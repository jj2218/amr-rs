//! Runtime CPU dispatch: the only module allowed to use `unsafe`.
//!
//! Each kernel listed in `kernels!` has one portable body (`#[inline(always)]`, in its own
//! module) and gets two copies: the portable build and an AVX2 build of the *same* source
//! (`#[target_feature(enable = "avx2")]`), picked per call from `is_x86_feature_detected!`
//! (std detects once and caches the answer in an atomic). Same source, so the same wrapping
//! integer results: the AVX2 copy only lets the compiler use 256-bit lanes. An `avx2 other`
//! entry runs another instance instead (a different lane width or a `const` flag): equal
//! results then rest on that body's contract, e.g. `first_min` is lane-width independent, and
//! the diff tests run on both paths.
//!
//! `--cfg portable` compiles the dispatch out, so the diff tests and the conformance
//! vectors run on both paths (`just test` and `just lint` cover both). Detection cost per call:
//! one cached atomic load and a branch, against kernels of 90 ns and up. The SSE2 kernels
//! (x86_64 baseline: both paths, except `vq_subvec_sse2`, the portable body of a dispatched
//! kernel) each have an array form for other targets; `--cfg amr_nb_arrays` selects those on
//! x86 too (`just test` runs the vectors with it).
//!
//! The frame decoder has no AVX2 copy (amr-wb's has, +4% there): here it measured 1.001×, as
//! its hot filters are SSE2 intrinsics already.
#![allow(unsafe_code)] // the approved scope: calls into #[target_feature] copies, below

use crate::basic_op::Flag;
use crate::encoder::codebook::{self as cb, RrTm};
use crate::encoder::gain_q::{
    mr475_search_portable, mr795_search_portable, qua_gain_exact_portable,
};
use crate::encoder::lpc::vq_pair_portable;
use crate::encoder::pitch_cl::convolve_portable;

/// Whether the AVX2 copies may run on this CPU.
#[cfg(all(target_arch = "x86_64", not(portable)))]
#[inline(always)]
fn avx2() -> bool {
    std::arch::is_x86_feature_detected!("avx2")
}

/// The first of the given paths.
#[cfg(all(target_arch = "x86_64", not(portable)))]
macro_rules! first {
    ($p:path $(, $rest:path)*) => {
        $p
    };
}

/// `pub fn name(args) = portable_body;` → the dispatching `name`. `= body, avx2 other;` gives
/// the AVX2 copy another instance of the same source (a `const` flag where the best form
/// depends on the lane width).
macro_rules! kernels {
    ($($(#[$m:meta])* $vis:vis fn $name:ident($($arg:ident: $ty:ty),* $(,)?) $(-> $ret:ty)?
        = $body:path $(, avx2 $avx:path)?;)*) => {$(
        $(#[$m])*
        #[allow(clippy::too_many_arguments)] // the C signatures
        $vis fn $name($($arg: $ty),*) $(-> $ret)? {
            #[cfg(all(target_arch = "x86_64", not(portable)))]
            if avx2() {
                #[target_feature(enable = "avx2")]
                fn avx2_copy($($arg: $ty),*) $(-> $ret)? {
                    first!($($avx,)? $body)($($arg),*)
                }
                // SAFETY: `avx2_copy` only requires AVX2, which `avx2()` just confirmed.
                return unsafe { avx2_copy($($arg),*) };
            }
            $body($($arg),*)
        }
    )*};
}

/// 8 × 8 transpose of i16 rows (`cb::transpose8_portable`'s three interleave rounds) as SSE2
/// `punpck` intrinsics: SSE2 is x86_64's baseline, so both paths run this one source. Written
/// over arrays, SLP rebuilt the network as 64 element loads (`vpinsrw`, ~2 uops each).
#[cfg(all(target_arch = "x86_64", not(amr_nb_arrays)))]
#[inline(always)]
pub(crate) fn transpose8(r: [[i16; 8]; 8]) -> [[i16; 8]; 8] {
    use core::arch::x86_64::{
        __m128i, _mm_unpackhi_epi16, _mm_unpackhi_epi32, _mm_unpackhi_epi64, _mm_unpacklo_epi16,
        _mm_unpacklo_epi32, _mm_unpacklo_epi64,
    };
    // The intrinsics need no runtime check; `[i16; 8]` and `__m128i` are 16 plain bytes each.
    // SAFETY: SSE2 is in the x86_64 baseline; the transmutes are between valid bit patterns.
    unsafe {
        let v: [__m128i; 8] = core::mem::transmute(r);
        let (lo, hi) = (_mm_unpacklo_epi16, _mm_unpackhi_epi16);
        let s = [
            lo(v[0], v[1]),
            hi(v[0], v[1]),
            lo(v[2], v[3]),
            hi(v[2], v[3]),
            lo(v[4], v[5]),
            hi(v[4], v[5]),
            lo(v[6], v[7]),
            hi(v[6], v[7]),
        ];
        let (lo2, hi2) = (_mm_unpacklo_epi32, _mm_unpackhi_epi32);
        let u = [
            lo2(s[0], s[2]),
            hi2(s[0], s[2]),
            lo2(s[1], s[3]),
            hi2(s[1], s[3]),
            lo2(s[4], s[6]),
            hi2(s[4], s[6]),
            lo2(s[5], s[7]),
            hi2(s[5], s[7]),
        ];
        let (lo4, hi4) = (_mm_unpacklo_epi64, _mm_unpackhi_epi64);
        let t = [
            lo4(u[0], u[4]),
            hi4(u[0], u[4]),
            lo4(u[1], u[5]),
            hi4(u[1], u[5]),
            lo4(u[2], u[6]),
            hi4(u[2], u[6]),
            lo4(u[3], u[7]),
            hi4(u[3], u[7]),
        ];
        core::mem::transmute(t)
    }
}

/// Other targets: the array form.
#[cfg(any(not(target_arch = "x86_64"), amr_nb_arrays))]
pub(crate) use cb::transpose8_portable as transpose8;

/// `norm_corr`'s lags on the fast paths (`pitch_cl::norm_lags_portable` is the contract and the
/// oracle) in SSE2: the 40-sample v stays in 5 registers from lag to lag; per lag the range
/// test (`pcmpgtw`), both dots as `pmaddwd` (one shared reduction), then the update
/// v'[j] = sat(((temp h[j]) >> h_fac) as i16 + v[j - 1]) from `pmulhw`/`pmullw` halves and
/// `paddsw` over v shifted one lane, v'[0] = temp >> scaling. (Saturating or not makes no
/// difference here: |prod| < 2^15 and |v| <= NORM_V_MAX, so a sum past ±32767 fails the next
/// lag's range test either way; no update follows the last lag.) Over arrays the update and
/// the dots ran as two out-of-line passes through memory per lag.
#[cfg(all(target_arch = "x86_64", not(amr_nb_arrays)))]
#[inline(always)]
pub(crate) fn norm_lags(
    p: &crate::encoder::pitch_cl::NormLags,
    v: &[i16; 40],
    exc_rev: &[i16],
    out: &mut [[i32; 2]],
) -> bool {
    use crate::encoder::pitch_cl::NORM_V_MAX;
    use core::arch::x86_64::{
        __m128i, _mm_add_epi32, _mm_adds_epi16, _mm_cmpgt_epi16, _mm_cvtsi32_si128,
        _mm_cvtsi128_si32, _mm_insert_epi16, _mm_loadu_si128, _mm_madd_epi16, _mm_movemask_epi8,
        _mm_mulhi_epi16, _mm_mullo_epi16, _mm_or_si128, _mm_set1_epi16, _mm_shuffle_epi32,
        _mm_sll_epi16, _mm_slli_si128, _mm_srl_epi16, _mm_srli_si128, _mm_unpackhi_epi32,
        _mm_unpacklo_epi32,
    };
    // Loads read 8 i16 at offsets 8q (q < 5) of the 40-sample arrays (h with h[0] = 0 in a
    // local copy).
    // SAFETY: SSE2 is in the x86_64 baseline; the accesses are in bounds (above).
    unsafe {
        let ld = |a: &[i16; 40], q: usize| _mm_loadu_si128(a.as_ptr().add(8 * q) as *const __m128i);
        let mut h0 = *p.h;
        h0[0] = 0; // lane 0 of the update is temp >> scaling, not a product
        let x: [__m128i; 5] = core::array::from_fn(|q| ld(p.xn, q));
        let h: [__m128i; 5] = core::array::from_fn(|q| ld(&h0, q));
        let mut v: [__m128i; 5] = core::array::from_fn(|q| ld(v, q));
        let (hi_max, lo_min) = (
            _mm_set1_epi16(NORM_V_MAX as i16),
            _mm_set1_epi16(-NORM_V_MAX as i16),
        );
        let lsh = _mm_cvtsi32_si128(16 - p.h_fac as i32);
        let rsh = _mm_cvtsi32_si128(p.h_fac as i32);
        let n = out.len();
        for (i, o) in out.iter_mut().enumerate() {
            let mut range = _mm_cvtsi32_si128(0);
            let (mut sx, mut sv) = (_mm_cvtsi32_si128(0), _mm_cvtsi32_si128(0));
            for q in 0..5 {
                range = _mm_or_si128(range, _mm_cmpgt_epi16(v[q], hi_max));
                range = _mm_or_si128(range, _mm_cmpgt_epi16(lo_min, v[q]));
                sx = _mm_add_epi32(sx, _mm_madd_epi16(x[q], v[q]));
                sv = _mm_add_epi32(sv, _mm_madd_epi16(v[q], v[q]));
            }
            if _mm_movemask_epi8(range) != 0 {
                return false;
            }
            // (sx0+sx2, sv0+sv2, sx1+sx3, sv1+sv3) -> lanes 0, 1: the two sums
            let t = _mm_add_epi32(_mm_unpacklo_epi32(sx, sv), _mm_unpackhi_epi32(sx, sv));
            let t = _mm_add_epi32(t, _mm_shuffle_epi32::<0xee>(t));
            let (s, s2) = (
                _mm_cvtsi128_si32(t),
                _mm_cvtsi128_si32(_mm_shuffle_epi32::<0x55>(t)),
            );
            if p.xx as i128 * s2 as i128 >= 1 << 60 {
                return false;
            }
            *o = [s, s2];
            if i + 1 < n {
                let temp = exc_rev[exc_rev.len() - 1 - i];
                if (temp as i32).abs() * p.hm >= 1 << (p.h_fac + 15) {
                    return false;
                }
                // ((temp h) >> h_fac) as i16 = (hi << (16 - h_fac)) | (lo >>> h_fac)
                let t = _mm_set1_epi16(temp);
                let prev = v;
                for q in 0..5 {
                    let (hi, lo) = (_mm_mulhi_epi16(t, h[q]), _mm_mullo_epi16(t, h[q]));
                    let prod = _mm_or_si128(_mm_sll_epi16(hi, lsh), _mm_srl_epi16(lo, rsh));
                    let mut sh = _mm_slli_si128::<2>(prev[q]);
                    if q > 0 {
                        sh = _mm_or_si128(sh, _mm_srli_si128::<14>(prev[q - 1]));
                    }
                    v[q] = _mm_adds_epi16(prod, sh);
                }
                v[0] = _mm_insert_epi16::<0>(v[0], (temp >> p.scaling) as i32);
            }
        }
        true
    }
}

/// Other targets: the array form.
#[cfg(any(not(target_arch = "x86_64"), amr_nb_arrays))]
pub(crate) use crate::encoder::pitch_cl::norm_lags_portable as norm_lags;

/// Four correlations at once (`dots4_portable` is the contract and the oracle):
/// sum_n x[n] w[off + k + n] for k = 0..4, wrapping, each x chunk loaded once for the four,
/// `pmaddwd` per window, the four sums reduced together (one transpose-add). For `comp_corr`
/// (auto-vectorised, LLVM kept x sign-extended to i32 across its lag loop and re-packed it
/// every lag: `vpackssdw` + `vextracti128` per 8 samples) and `autocorr`.
#[cfg(all(target_arch = "x86_64", not(amr_nb_arrays)))]
#[inline(always)]
pub(crate) fn dots4<const L: usize>(x: &[i16; L], w: &[i16], off: usize) -> [i32; 4] {
    use core::arch::x86_64::{
        __m128i, _mm_add_epi32, _mm_loadu_si128, _mm_madd_epi16, _mm_setzero_si128,
        _mm_storeu_si128, _mm_unpackhi_epi32, _mm_unpackhi_epi64, _mm_unpacklo_epi32,
        _mm_unpacklo_epi64,
    };
    const { assert!(L.is_multiple_of(8)) };
    let w = &w[off..off + L + 3];
    // Loads read 8 i16 at x[c..] and w[k + c..] with c + 8 <= L, k <= 3 (w: L + 3 long).
    // SAFETY: SSE2 is in the x86_64 baseline; the accesses are in bounds (above).
    unsafe {
        let ld = |p: *const i16| _mm_loadu_si128(p as *const __m128i);
        let mut acc = [_mm_setzero_si128(); 4];
        for c in (0..L).step_by(8) {
            let xv = ld(x.as_ptr().add(c));
            for (k, a) in acc.iter_mut().enumerate() {
                *a = _mm_add_epi32(*a, _mm_madd_epi16(xv, ld(w.as_ptr().add(k + c))));
            }
        }
        let (l32, h32) = (_mm_unpacklo_epi32, _mm_unpackhi_epi32);
        let u0 = _mm_add_epi32(l32(acc[0], acc[1]), h32(acc[0], acc[1]));
        let u1 = _mm_add_epi32(l32(acc[2], acc[3]), h32(acc[2], acc[3]));
        let r = _mm_add_epi32(_mm_unpacklo_epi64(u0, u1), _mm_unpackhi_epi64(u0, u1));
        let mut out = [0i32; 4];
        _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, r);
        out
    }
}

/// `dots4`' array form (its non-x86 body and oracle).
#[cfg_attr(all(target_arch = "x86_64", not(amr_nb_arrays)), allow(dead_code))]
pub(crate) fn dots4_portable<const L: usize>(x: &[i16; L], w: &[i16], off: usize) -> [i32; 4] {
    core::array::from_fn(|k| {
        x.iter()
            .zip(&w[off + k..off + k + L])
            .fold(0i32, |s, (&a, &b)| s.wrapping_add(a as i32 * b as i32))
    })
}

/// Other targets: the array form.
#[cfg(any(not(target_arch = "x86_64"), amr_nb_arrays))]
pub(crate) use dots4_portable as dots4;

/// `vq_subvec_portable` (the contract and the oracle) in SSE2 i16 lanes, 8 codewords per
/// step: t = (w (r - c)) >> 15 as i16 is (`pmulhw` << 1) | (`pmullw` >>> 15), and each
/// codeword's sum of t^2 over its components is `pmaddwd` of the components interleaved
/// pairwise (a pair wraps only at (-32768)^2, as the wrapping i32 sum does). SSE2 has no
/// `pmulld`: in i32 lanes the search was 3.7x its AVX2 copy. The portable build's body.
#[cfg(all(target_arch = "x86_64", not(amr_nb_arrays)))]
#[inline(always)]
pub(crate) fn vq_subvec_sse2<const D: usize>(
    lsf_r1: &mut [i16],
    cols: &[&[i16]; D],
    wf1: &[i16],
) -> i16 {
    use core::arch::x86_64::{
        __m128i, _mm_add_epi32, _mm_and_si128, _mm_andnot_si128, _mm_cmpgt_epi32, _mm_loadu_si128,
        _mm_madd_epi16, _mm_mulhi_epi16, _mm_mullo_epi16, _mm_or_si128, _mm_set1_epi16,
        _mm_set1_epi32, _mm_setr_epi32, _mm_setzero_si128, _mm_slli_epi16, _mm_srli_epi16,
        _mm_storeu_si128, _mm_sub_epi16, _mm_unpackhi_epi16, _mm_unpacklo_epi16,
    };
    const { assert!(D == 3 || D == 4) };
    let n = cols[0].len();
    assert!(n.is_multiple_of(8) && cols.iter().all(|c| c.len() == n));
    let (mut best, mut bidx) = ([0i32; 8], [0i32; 8]);
    // Loads read 8 i16 at cols[k][b..b + 8] (b + 8 <= n, checked above); stores write the
    // two `[i32; 8]` arrays' halves.
    // SAFETY: SSE2 is in the x86_64 baseline; the accesses are in bounds (above).
    unsafe {
        let r: [__m128i; D] = core::array::from_fn(|k| _mm_set1_epi16(lsf_r1[k]));
        let w: [__m128i; D] = core::array::from_fn(|k| _mm_set1_epi16(wf1[k]));
        let (mut blo, mut bhi) = (_mm_set1_epi32(i32::MAX), _mm_set1_epi32(i32::MAX));
        let (mut ilo, mut ihi) = (_mm_setzero_si128(), _mm_setzero_si128());
        let (mut idx_lo, mut idx_hi) = (_mm_setr_epi32(0, 1, 2, 3), _mm_setr_epi32(4, 5, 6, 7));
        let eight = _mm_set1_epi32(8);
        let sel = |m: __m128i, a: __m128i, b: __m128i| {
            _mm_or_si128(_mm_and_si128(m, a), _mm_andnot_si128(m, b))
        };
        for b in (0..n).step_by(8) {
            let t: [__m128i; D] = core::array::from_fn(|k| {
                let c = _mm_loadu_si128(cols[k].as_ptr().add(b) as *const __m128i);
                let d = _mm_sub_epi16(r[k], c);
                let (hi, lo) = (_mm_mulhi_epi16(w[k], d), _mm_mullo_epi16(w[k], d));
                _mm_or_si128(_mm_slli_epi16::<1>(hi), _mm_srli_epi16::<15>(lo))
            });
            let t3 = if D == 4 {
                t[D - 1]
            } else {
                _mm_setzero_si128()
            };
            let (p0, p1) = (
                (
                    _mm_unpacklo_epi16(t[0], t[1]),
                    _mm_unpackhi_epi16(t[0], t[1]),
                ),
                (_mm_unpacklo_epi16(t[2], t3), _mm_unpackhi_epi16(t[2], t3)),
            );
            let dlo = _mm_add_epi32(_mm_madd_epi16(p0.0, p0.0), _mm_madd_epi16(p1.0, p1.0));
            let dhi = _mm_add_epi32(_mm_madd_epi16(p0.1, p0.1), _mm_madd_epi16(p1.1, p1.1));
            // strict <: a lane keeps its first minimum
            let (mlo, mhi) = (_mm_cmpgt_epi32(blo, dlo), _mm_cmpgt_epi32(bhi, dhi));
            (blo, bhi) = (sel(mlo, dlo, blo), sel(mhi, dhi, bhi));
            (ilo, ihi) = (sel(mlo, idx_lo, ilo), sel(mhi, idx_hi, ihi));
            (idx_lo, idx_hi) = (_mm_add_epi32(idx_lo, eight), _mm_add_epi32(idx_hi, eight));
        }
        _mm_storeu_si128(best.as_mut_ptr() as *mut __m128i, blo);
        _mm_storeu_si128(best.as_mut_ptr().add(4) as *mut __m128i, bhi);
        _mm_storeu_si128(bidx.as_mut_ptr() as *mut __m128i, ilo);
        _mm_storeu_si128(bidx.as_mut_ptr().add(4) as *mut __m128i, ihi);
    }
    let index = crate::encoder::layout::first_min(&best, &bidx);
    for k in 0..D {
        lsf_r1[k] = cols[k][index as usize];
    }
    index as i16
}

/// Other targets: the array form.
#[cfg(any(not(target_arch = "x86_64"), amr_nb_arrays))]
#[inline(always)]
pub(crate) fn vq_subvec_sse2<const D: usize>(
    lsf_r1: &mut [i16],
    cols: &[&[i16]; D],
    wf1: &[i16],
) -> i16 {
    crate::encoder::lpc::vq_subvec_portable::<D, 16>(lsf_r1, cols, wf1)
}

/// `residu`'s block of 8 outputs (`filt::residu_block_portable` is the contract and the
/// oracle) over the window w = x[n..n + 18]: wrapping sums 0x800 + sum_k a[k] w[10 - k + i],
/// the taps in pairs (k, k + 1) as `pmaddwd` over windows 10 - k and 9 - k interleaved, the
/// last tap against zeros.
#[cfg(all(target_arch = "x86_64", not(amr_nb_arrays)))]
#[inline(always)]
pub(crate) fn residu_block(a: &[i16; 11], w: &[i16; 18]) -> [i32; 8] {
    // The windows and coefficient pairs below are laid out for order 10.
    const { assert!(crate::post::M == 10) };
    use core::arch::x86_64::{
        __m128i, _mm_add_epi32, _mm_loadu_si128, _mm_madd_epi16, _mm_set1_epi32, _mm_setzero_si128,
        _mm_shuffle_epi32, _mm_storeu_si128, _mm_unpackhi_epi16, _mm_unpacklo_epi16,
    };
    let pair = |lo: i16, hi: i16| (lo as u16 as u32 | (hi as u16 as u32) << 16) as i32;
    // Loads read 8 i16 at w[j..j + 8] (j <= 10) and a[0..8].
    // SAFETY: SSE2 is in the x86_64 baseline; the accesses are in bounds (above).
    unsafe {
        let ld = |p: *const i16| _mm_loadu_si128(p as *const __m128i);
        let a0 = ld(a.as_ptr());
        let c = [
            _mm_shuffle_epi32::<0x00>(a0),
            _mm_shuffle_epi32::<0x55>(a0),
            _mm_shuffle_epi32::<0xaa>(a0),
            _mm_shuffle_epi32::<0xff>(a0),
            _mm_set1_epi32(pair(a[8], a[9])),
            _mm_set1_epi32(pair(a[10], 0)),
        ];
        let (mut lo, mut hi) = (_mm_set1_epi32(0x800), _mm_set1_epi32(0x800));
        for (p, &c) in c.iter().enumerate() {
            let k = 2 * p;
            let u = ld(w.as_ptr().add(10 - k));
            let v = if k < 10 {
                ld(w.as_ptr().add(9 - k))
            } else {
                _mm_setzero_si128()
            };
            lo = _mm_add_epi32(lo, _mm_madd_epi16(_mm_unpacklo_epi16(u, v), c));
            hi = _mm_add_epi32(hi, _mm_madd_epi16(_mm_unpackhi_epi16(u, v), c));
        }
        let mut out = [0i32; 8];
        _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, lo);
        _mm_storeu_si128(out.as_mut_ptr().add(4) as *mut __m128i, hi);
        out
    }
}

/// Other targets: the array form.
#[cfg(any(not(target_arch = "x86_64"), amr_nb_arrays))]
pub(crate) use crate::post::filt::residu_block_portable as residu_block;

/// `pred_lt_3or6`'s block of 8 outputs (`pred_lt::block_portable` is the contract and the
/// oracle): wrapping sums 0x4000 + sum_t src[i + t] k[t], the taps in pairs as `pmaddwd`
/// (windows t and t + 1 interleaved against the broadcast pair (k[t], k[t + 1])). Over arrays
/// the pair form went scalar and the per-tap form is `pmullw`/`pmulhw` + unpacks per tap.
#[cfg(all(target_arch = "x86_64", not(amr_nb_arrays)))]
#[inline(always)]
pub(crate) fn pred_block(src: &[i16; 27], k: &[i16; 20]) -> [i32; 8] {
    // The loads below are laid out for these sizes.
    const { assert!(crate::excitation::pred_lt::BLOCK_IN == 27) };
    use core::arch::x86_64::{
        __m128i, _mm_add_epi32, _mm_loadu_si128, _mm_madd_epi16, _mm_set1_epi32, _mm_shuffle_epi32,
        _mm_storeu_si128, _mm_unpackhi_epi16, _mm_unpacklo_epi16,
    };
    // Loads read 8 i16 at src[t..t + 8] (t <= 19 < 27 - 7) and k[0..8], k[8..16], k[12..20].
    // SAFETY: SSE2 is in the x86_64 baseline; the accesses are in bounds (above).
    unsafe {
        let ld = |p: *const i16| _mm_loadu_si128(p as *const __m128i);
        let (k0, k1, k2) = (
            ld(k.as_ptr()),
            ld(k.as_ptr().add(8)),
            ld(k.as_ptr().add(12)),
        );
        let c = [
            _mm_shuffle_epi32::<0x00>(k0),
            _mm_shuffle_epi32::<0x55>(k0),
            _mm_shuffle_epi32::<0xaa>(k0),
            _mm_shuffle_epi32::<0xff>(k0),
            _mm_shuffle_epi32::<0x00>(k1),
            _mm_shuffle_epi32::<0x55>(k1),
            _mm_shuffle_epi32::<0xaa>(k1),
            _mm_shuffle_epi32::<0xff>(k1),
            _mm_shuffle_epi32::<0xaa>(k2),
            _mm_shuffle_epi32::<0xff>(k2),
        ];
        let (mut lo, mut hi) = (_mm_set1_epi32(0x4000), _mm_set1_epi32(0x4000));
        for (p, &c) in c.iter().enumerate() {
            let (a, b) = (ld(src.as_ptr().add(2 * p)), ld(src.as_ptr().add(2 * p + 1)));
            lo = _mm_add_epi32(lo, _mm_madd_epi16(_mm_unpacklo_epi16(a, b), c));
            hi = _mm_add_epi32(hi, _mm_madd_epi16(_mm_unpackhi_epi16(a, b), c));
        }
        let mut out = [0i32; 8];
        _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, lo);
        _mm_storeu_si128(out.as_mut_ptr().add(4) as *mut __m128i, hi);
        out
    }
}

/// Other targets: the array form.
#[cfg(any(not(target_arch = "x86_64"), amr_nb_arrays))]
pub(crate) use crate::excitation::pred_lt::block_portable as pred_block;

/// `syn_fast`'s blocks of 8 outputs (`filt::syn_blocks_portable` is the contract and the
/// oracle) in SSE2. Per block: the x term and the terms on the 10 outputs before it as 6
/// `pmaddwd` pairs per half; then per output the y[n-1], y[n-2] terms as scalar `imul`s (the
/// chain), and its terms in later outputs (k >= 3) as one `pmaddwd` per half off the chain.
/// Output m reads its lane right after the update of output m - 3, the last one it needs:
/// read from the latest accumulators, the (zero) updates of m - 2 and m - 1 put a `pmaddwd` on
/// the chain (1.2× slower). The outputs and the history stay in registers from block
/// to block (stored as elements and reloaded as vectors they missed store forwarding). Over
/// arrays, LLVM spilled the lanes (2× slower) and scalarised pair-interleaved columns.
#[cfg(all(target_arch = "x86_64", not(amr_nb_arrays)))]
#[inline(always)]
pub(crate) fn syn_blocks(
    t: &crate::post::filt::SynTab,
    old: &[i16; 16],
    x: &[[i16; 8]],
    y: &mut [[i16; 8]],
) -> u32 {
    use crate::post::filt::syn_out_of_range;
    use core::arch::x86_64::{
        __m128i, _mm_and_si128, _mm_cvtsi128_si32, _mm_loadu_si128, _mm_madd_epi16, _mm_or_si128,
        _mm_packs_epi32, _mm_set1_epi16, _mm_set1_epi32, _mm_shuffle_epi32, _mm_slli_epi32,
        _mm_storeu_si128, _mm_sub_epi32, _mm_unpackhi_epi16, _mm_unpacklo_epi16,
        _mm_unpacklo_epi32, _mm_unpacklo_epi64,
    };
    // Every load and store touches 16 bytes inside an array of at least that size (`[i32; 8]`
    // halves, `[i16; 8]`, `[i16; 16]` halves).
    // SAFETY: SSE2 is in the x86_64 baseline; the accesses are in bounds (above).
    unsafe {
        let ld = |p: *const i32| _mm_loadu_si128(p as *const __m128i);
        let ld16 = |p: *const i16| _mm_loadu_si128(p as *const __m128i);
        let mul = |c: &[i32; 8], h: usize, w: __m128i| _mm_madd_epi16(ld(c.as_ptr().add(h)), w);
        let (one, a0x) = (_mm_set1_epi16(1), _mm_set1_epi32(t.a0x));
        let low = _mm_set1_epi32(0xffff);
        // (u, v) as one i32 lane pair from two broadcasts, for `pmaddwd`.
        let pair =
            |u: __m128i, v: __m128i| _mm_or_si128(_mm_and_si128(u, low), _mm_slli_epi32::<16>(v));
        // w[q]: (old[2q], old[2q + 1]) broadcast
        let (o0, o1) = (ld16(old.as_ptr()), ld16(old.as_ptr().add(8)));
        let mut w = [
            _mm_shuffle_epi32::<0x00>(o0),
            _mm_shuffle_epi32::<0x55>(o0),
            _mm_shuffle_epi32::<0xaa>(o0),
            _mm_shuffle_epi32::<0xff>(o0),
            _mm_shuffle_epi32::<0x00>(o1),
        ];
        let (mut p2, mut p1) = (old[8] as i32, old[9] as i32);
        let mut bad = 0u32;
        for (xb, yo) in x.iter().zip(y) {
            let xv = ld16(xb.as_ptr());
            let mut lo = _mm_madd_epi16(_mm_unpacklo_epi16(xv, one), a0x);
            let mut hi = _mm_madd_epi16(_mm_unpackhi_epi16(xv, one), a0x);
            for (c, &w) in t.colp[..4].iter().zip(&w[..4]) {
                lo = _mm_sub_epi32(lo, mul(c, 0, w));
                hi = _mm_sub_epi32(hi, mul(c, 4, w));
            }
            // Output 0 doesn't read pair 4 (the previous block's last two outputs: k = 1, 2
            // there, on the chain), so it goes in after output 0's lane is read.
            let lo_all = _mm_sub_epi32(lo, mul(&t.colp[4], 0, w[4]));
            let hi = _mm_sub_epi32(hi, mul(&t.colp[4], 4, w[4]));
            // Output m from its accumulator (in lane 0); returns it broadcast.
            let mut out = |lane0: __m128i| -> __m128i {
                let s = _mm_cvtsi128_si32(lane0);
                let s = s
                    .wrapping_sub(t.a2.wrapping_mul(p2))
                    .wrapping_sub(t.a1.wrapping_mul(p1));
                bad |= syn_out_of_range(s);
                let v = s >> 12;
                (p2, p1) = (p1, v);
                _mm_set1_epi32(v)
            };
            let u = |m: usize| &t.upd[m];
            let b0 = out(lo);
            let lo3 = _mm_sub_epi32(lo_all, mul(u(0), 0, b0));
            let h0 = _mm_sub_epi32(hi, mul(u(0), 4, b0));
            let b1 = out(_mm_shuffle_epi32::<1>(lo_all));
            let h1 = _mm_sub_epi32(h0, mul(u(1), 4, b1));
            let b2 = out(_mm_shuffle_epi32::<2>(lo_all));
            let h2 = _mm_sub_epi32(h1, mul(u(2), 4, b2));
            let b3 = out(_mm_shuffle_epi32::<3>(lo3));
            let h3 = _mm_sub_epi32(h2, mul(u(3), 4, b3));
            let b4 = out(h1);
            let h4 = _mm_sub_epi32(h3, mul(u(4), 4, b4));
            let b5 = out(_mm_shuffle_epi32::<1>(h2));
            let b6 = out(_mm_shuffle_epi32::<2>(h3));
            let b7 = out(_mm_shuffle_epi32::<3>(h4));
            let (l, h) = (_mm_unpacklo_epi32, _mm_unpacklo_epi64);
            let ys = _mm_packs_epi32(h(l(b0, b1), l(b2, b3)), h(l(b4, b5), l(b6, b7)));
            _mm_storeu_si128(yo.as_mut_ptr() as *mut __m128i, ys);
            // next history pairs: (old[8], old[9]), (y0, y1), .., (y6, y7), each from its own
            // outputs (from the packed `ys`, all of them waited for y7)
            w = [w[4], pair(b0, b1), pair(b2, b3), pair(b4, b5), pair(b6, b7)];
        }
        bad
    }
}

/// Other targets: the array form.
#[cfg(any(not(target_arch = "x86_64"), amr_nb_arrays))]
pub(crate) use crate::post::filt::syn_blocks_portable as syn_blocks;

kernels! {
    /// `cor_h` into the 5-track modes' track-major layout (see `cor_h_tm_portable`): the AVX2
    /// copy keeps the multiply form of the sign scaling, the portable build the ±32767 form.
    pub(crate) fn cor_h_tm(h: &[i16], sign: &[i16; cb::L_CODE], rr: &mut RrTm, ov: &mut Flag)
        = cb::cor_h_tm_portable::<5, 8, 8, true>, avx2 cb::cor_h_tm_portable::<5, 8, 8, false>;
    /// `cor_h_tm` for two blocks only (MR475/MR515's track pairs, see `cor_h_blocks_portable`).
    pub(crate) fn cor_h_blocks(h: &[i16], sign: &[i16; cb::L_CODE], rr: &mut RrTm,
        blocks: &[[usize; 2]; 2], ov: &mut Flag) = cb::cor_h_blocks_portable;
    /// MR122's pulse search (see `search_10i40_portable`; `Pair::best`'s i32 form with AVX2).
    pub(crate) fn search_10i40(dn: &[i16; cb::L_CODE], rr: &RrTm, ipos: &mut [i16; 10],
        pos_max: &[i16; 5], codvec: &mut [usize; 10])
        = cb::search_10i40_portable::<false>, avx2 cb::search_10i40_portable::<true>;
    /// MR102's `cor_h` + pulse search (see `cor_h_search_8i40_portable`): rows of 10 lanes
    /// portable, 16 with AVX2.
    pub(crate) fn cor_h_search_8i40(h: &[i16], sign: &[i16; cb::L_CODE], dn: &[i16; cb::L_CODE],
        ipos: &mut [i16; 8], pos_max: &[i16; 4], codvec: &mut [usize; 8], ov: &mut Flag)
        = cb::cor_h_search_8i40_portable::<10, true, false>,
        avx2 cb::cor_h_search_8i40_portable::<16, false, true>;
    /// The step-5 searches (see `search_2i40_portable`, `search_3i40_portable`,
    /// `search_4i40_portable`).
    pub(crate) fn search_2i40(dn: &[i16; cb::L_CODE], rr: &RrTm, pairs: &[[usize; 2]])
        -> [usize; 2] = cb::search_2i40_portable;
    pub(crate) fn search_3i40(dn: &[i16; cb::L_CODE], dn2: &[i16; cb::L_CODE], rr: &RrTm)
        -> [usize; 3] = cb::search_3i40_portable;
    pub(crate) fn search_4i40(dn: &[i16; cb::L_CODE], dn2: &[i16; cb::L_CODE], rr: &RrTm)
        -> [usize; 4] = cb::search_4i40_portable;
    /// `cor_h_x`'s correlations (see `corr_hx_portable`).
    pub(crate) fn corr_hx(h: &[i16], x: &[i16]) -> [i32; cb::L_CODE] = cb::corr_hx_portable;
    /// The filtered codevector (see `filtered_code_portable`).
    pub(crate) fn filtered_code(h: &[i16], codvec: &[usize], sign: &[i16], y: &mut [i16],
        ov: &mut Flag) = cb::filtered_code_portable;
    /// `Convolve` (see `convolve_portable`).
    pub(crate) fn convolve(x: &[i16], h: &[i16], y: &mut [i16]) = convolve_portable;
    /// `Qua_gain`'s search without saturation (see `qua_gain_exact_portable`).
    pub(crate) fn qua_gain_exact(cols: &[&[i16]; 2], gcode0: i16, coeff: &[i16; 5],
        coeff_lo: &[i16; 5], gp_limit: i16) -> Option<usize>
        = qua_gain_exact_portable::<8>;
    /// `MR475_gain_quant`'s search (see `mr475_search_portable`).
    pub(crate) fn mr475_search(coeff: &[i16; 10], coeff_lo: &[i16; 10], gcode0: &[i16; 2],
        gp_limit: i16) -> Option<usize>
        = mr475_search_portable::<16>, avx2 mr475_search_portable::<8>;
    /// `MR795_gain_code_quant3`'s search (see `mr795_search_portable`).
    pub(crate) fn mr795_search(coeff: &[i16; 5], coeff_lo: &[i16; 5], l_tmp0: &[i32; 3],
        gcode0: i16, g_pitch: &[i16; 3]) -> Option<usize>
        = mr795_search_portable::<8>, avx2 mr795_search_portable::<16>;
    /// MR122's `Vq_subvec` / `Vq_subvec_s` in lanes (see `vq_pair_portable`).
    pub(crate) fn vq_pair(lsf_r1: &[i16], lsf_r2: &[i16], cols: &[&[i16]; 4], wf1: &[i16],
        wf2: &[i16]) -> Option<usize>
        = vq_pair_portable::<false, 16>, avx2 vq_pair_portable::<false, 32>;
    pub(crate) fn vq_pair_s(lsf_r1: &[i16], lsf_r2: &[i16], cols: &[&[i16]; 4], wf1: &[i16],
        wf2: &[i16]) -> Option<usize>
        = vq_pair_portable::<true, 16>, avx2 vq_pair_portable::<true, 32>;
    /// `Vq_subvec3`, `Vq_subvec4` (see `vq_subvec_portable`).
    pub(crate) fn vq_subvec3(lsf_r1: &mut [i16], cols: &[&[i16]; 3], wf1: &[i16]) -> i16
        = vq_subvec_sse2::<3>, avx2 crate::encoder::lpc::vq_subvec_portable::<3, 32>;
    pub(crate) fn vq_subvec4(lsf_r1: &mut [i16], cols: &[&[i16]; 4], wf1: &[i16]) -> i16
        = vq_subvec_sse2::<4>, avx2 crate::encoder::lpc::vq_subvec_portable::<4, 32>;
}
