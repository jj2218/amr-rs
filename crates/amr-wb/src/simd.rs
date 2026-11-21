//! Runtime CPU dispatch: the only module allowed to use `unsafe`.
//!
//! The frame decoder has one source and two copies: the portable build and an AVX2 build of
//! the *same* source (`#[target_feature(enable = "avx2")]`), picked per frame from
//! `is_x86_feature_detected!` (std caches the answer in an atomic). Same source, so the same
//! wrapping integer results: the AVX2 copy only lets the compiler use 256-bit lanes. The
//! encoder's hot kernels (`kernels!`) are dispatched the same way per call: each has one
//! portable body (`#[inline(always)]`, `*_portable` in its module, with the helpers it runs)
//! and an AVX2 copy of it. A copy gets AVX2 only for what inlines into it. The encoder frame
//! has no AVX2 copy: it measured 0.988×, as its hot kernels don't inline into the frame.
//! `--cfg portable` compiles the dispatch out (`just test` runs both paths).
//!
//! Where LLVM's lanes stay poor, a kernel gets an AVX2 intrinsics body here (`dots4`, `norm_lags`,
//! `fir32`, `phase_conv`, `oversamp_interpolate`), dispatched per call the same way; its `*_portable` body is the contract and oracle.
#![allow(unsafe_code)] // the approved scope: #[target_feature] copies and AVX2 kernels, below

use crate::cnst::{L_FRAME, L_SUBFR, M, RxFrameType};
use crate::decoder::State;
use crate::encoder::kernels as k;
#[cfg(all(target_arch = "x86_64", not(portable)))]
use crate::synthesis::HIST_PAIRS;
use crate::synthesis::L_FIR;

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

/// `State::decode` on the best instruction set of this CPU.
pub(crate) fn decode(
    st: &mut State,
    mode: usize,
    prms: &[i16],
    frame_type: RxFrameType,
    out: &mut [i16; L_FRAME],
) {
    #[cfg(all(target_arch = "x86_64", not(portable)))]
    if avx2() {
        #[target_feature(enable = "avx2")]
        fn avx2_copy(
            st: &mut State,
            mode: usize,
            prms: &[i16],
            frame_type: RxFrameType,
            out: &mut [i16; L_FRAME],
        ) {
            st.decode(mode, prms, frame_type, out)
        }
        // SAFETY: `avx2_copy` only requires AVX2, which `avx2()` just confirmed.
        return unsafe { avx2_copy(st, mode, prms, frame_type, out) };
    }
    st.decode(mode, prms, frame_type, out)
}

kernels! {
    /// `ACELP_4t64_fx` (see `acelp4::acelp_4t64_fx_portable`).
    pub(crate) fn acelp_4t64_fx(dn: &mut [i16; L_SUBFR], cn: &[i16], hh: &[i16],
        code: &mut [i16; L_SUBFR], y: &mut [i16; L_SUBFR], nbbits: i16, ser_size: i16,
        index_out: &mut [i16; 8]) = k::acelp_4t64_fx_portable;
    /// `ACELP_2t64_fx` (see `acelp2::acelp_2t64_fx_portable`).
    pub(crate) fn acelp_2t64_fx(dn: &mut [i16; L_SUBFR], cn: &[i16], hh: &[i16],
        code: &mut [i16; L_SUBFR], y: &mut [i16; L_SUBFR], index: &mut i16)
        = k::acelp_2t64_fx_portable;
    /// `cor_h_x` (see `acelp2::cor_h_x_portable`).
    pub(crate) fn cor_h_x(h: &[i16], x: &[i16], dn: &mut [i16; L_SUBFR]) = k::cor_h_x_portable;
    /// `Convolve`, `Residu`, `Syn_filt` (see `filters`).
    pub(crate) fn convolve(x: &[i16], h: &[i16], y: &mut [i16], l: usize) = k::convolve_portable;
    pub(crate) fn residu(a: &[i16], x: &[i16], at: usize, y: &mut [i16], lg: usize)
        = k::residu_portable;
    pub(crate) fn syn_filt(a: &[i16], x: Option<&[i16]>, y: &mut [i16], mem: &mut [i16],
        update: bool) = k::syn_filt_portable;
    /// `Pitch_fr4` (see `pitch::pitch_fr4_portable`).
    pub(crate) fn pitch_fr4(exc: &[i16], at: usize, xn: &[i16], h: &[i16], t0_min: i16,
        t0_max: i16, pit_frac: &mut i16, i_subfr: usize, t0_fr2: i16, t0_fr1: i16,
        l_subfr: usize) -> i16 = k::pitch_fr4_portable;
    /// `Decim_12k8` (see `preprocess::decim_12k8_portable`).
    pub(crate) fn decim_12k8(sig16k: &[i16], sig12k8: &mut [i16],
        mem: &mut [i16; 2 * crate::encoder::kernels::NB_COEF_DOWN]) = k::decim_12k8_portable;
    /// `Pitch_med_ol` (see `pitch_ol::pitch_med_ol_portable`).
    pub(crate) fn pitch_med_ol(wsp: &[i16], at: usize, l_min: i16, l_max: i16, l_frame: usize,
        l_0: i16, gain: &mut i16, hp_wsp_mem: &mut [i16; 9], old_hp_wsp: &mut [i16],
        wght_flg: i16) -> i16 = k::pitch_med_ol_portable;
    /// `Autocorr`, `Az_isp`, `Qpisf_2s_46b`, `Qpisf_2s_36b` (see `lpc`).
    pub(crate) fn autocorr(x: &[i16], r_h: &mut [i16; M + 1], r_l: &mut [i16; M + 1])
        = k::autocorr_portable;
    pub(crate) fn az_isp(a: &[i16], isp: &mut [i16; M], old_isp: &[i16; M])
        = k::az_isp_portable;
    pub(crate) fn qpisf_2s_46b(isf1: &[i16; M], isf_q: &mut [i16; M], past_isfq: &mut [i16; M],
        indice: &mut [i16; 7], nb_surv: usize) = k::qpisf_2s_46b_portable;
    pub(crate) fn qpisf_2s_36b(isf1: &[i16; M], isf_q: &mut [i16; M], past_isfq: &mut [i16; M],
        indice: &mut [i16; 5], nb_surv: usize) = k::qpisf_2s_36b_portable;
}

/// `acelp4::dots4_portable` (the contract and the oracle): the 4 dots
/// `Σ_{j<N} h[j]·p[o + 4k + j]`, on AVX2 as `pmaddwd`s sharing each 16-term `h` load across the 4 positions, then
/// one transposed reduction (the scalar form reduced each dot alone; its 16-term group ran
/// scalar). Wrapping like the portable sums: `pmaddwd`'s pair sum wraps only at
/// (-32768)^2 · 2, as the wrapping i32 sum does.
#[inline(always)]
pub(crate) fn dots4<const N: usize>(h: &[i16; L_SUBFR], p: &[i16], o: usize) -> [i32; 4] {
    #[cfg(all(target_arch = "x86_64", not(portable)))]
    if avx2() {
        // SAFETY: `dots4_avx2` only requires AVX2, which `avx2()` just confirmed.
        return unsafe { dots4_avx2::<N>(h, p, o) };
    }
    k::dots4_portable::<N>(h, p, o)
}

#[cfg(all(target_arch = "x86_64", not(portable)))]
#[target_feature(enable = "avx2")]
#[inline]
fn dots4_avx2<const N: usize>(h: &[i16; L_SUBFR], p: &[i16], o: usize) -> [i32; 4] {
    use core::arch::x86_64::{
        __m128i, __m256i, _mm_add_epi32, _mm_storeu_si128, _mm256_add_epi32,
        _mm256_castsi256_si128, _mm256_extracti128_si256, _mm256_loadu_si256, _mm256_madd_epi16,
        _mm256_setzero_si256, _mm256_unpackhi_epi32, _mm256_unpackhi_epi64, _mm256_unpacklo_epi32,
        _mm256_unpacklo_epi64,
    };
    const { assert!(N.is_multiple_of(16) && N <= L_SUBFR) };
    let p = &p[o..o + 12 + N]; // checked once: the loads' slices below then check for free
    let ld = |q: &[i16; 16]| -> __m256i {
        // SAFETY: 16 i16 behind a reference; AVX2 is enabled in this function.
        unsafe { _mm256_loadu_si256(q.as_ptr() as *const __m256i) }
    };
    let mut acc = [_mm256_setzero_si256(); 4];
    for c in (0..N).step_by(16) {
        let hv = ld(h[c..c + 16].try_into().unwrap());
        for (k, a) in acc.iter_mut().enumerate() {
            let w = ld(p[4 * k + c..4 * k + c + 16].try_into().unwrap());
            *a = _mm256_add_epi32(*a, _mm256_madd_epi16(hv, w));
        }
    }
    let (l32, h32) = (_mm256_unpacklo_epi32, _mm256_unpackhi_epi32);
    let u0 = _mm256_add_epi32(l32(acc[0], acc[1]), h32(acc[0], acc[1]));
    let u1 = _mm256_add_epi32(l32(acc[2], acc[3]), h32(acc[2], acc[3]));
    let r = _mm256_add_epi32(_mm256_unpacklo_epi64(u0, u1), _mm256_unpackhi_epi64(u0, u1));
    let r = _mm_add_epi32(_mm256_castsi256_si128(r), _mm256_extracti128_si256::<1>(r));
    let mut out = [0i32; 4];
    // SAFETY: `out` holds 4 i32 (16 bytes).
    unsafe { _mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, r) };
    out
}

/// `synthesis::fir32_portable` (the contract and the oracle): a 32-tap FIR's fast path,
/// out[n] = clamp16((Σ_{k<32} c[k]·x[n + k] + 2^(SH-1)) >> SH) for n < NO (`x`: NO + 31 samples
/// or more). On AVX2 as tap pairs: x interleaved once into X2[i] = (x[i], x[i + 1]), then per
/// pair one `pmaddwd` of 8 outputs' X2 window with the broadcast (c[2m], c[2m + 1]) (the
/// transposed portable form multiplies in i32 lanes, one tap per pass). Exact where the
/// portable form is: |Σ| + 2^(SH-1) < 2^31 (wrapping adds then give the true sum).
#[inline(always)]
pub(crate) fn fir32<const NO: usize, const SH: i32>(
    x: &[i16],
    c: &[i16; L_FIR + 2],
    out: &mut [i16; NO],
) {
    const { assert!(NO.is_multiple_of(16) && NO <= 80 && SH >= 1 && SH <= 16) };
    #[cfg(all(target_arch = "x86_64", not(portable)))]
    if avx2() {
        // SAFETY: `fir32_avx2` only requires AVX2, which `avx2()` just confirmed.
        return unsafe { fir32_avx2::<NO, SH>(x, c, out) };
    }
    crate::synthesis::fir32_portable::<NO, SH>(x, c, out)
}

/// `synthesis::hist_terms_portable` (the contract and the oracle; only builds with the AVX2
/// copy use this layout): the terms of the N outputs before a block of 8. On AVX2 one
/// `pmaddwd` per history pair, the coefficient pairs of the 8 outputs against the broadcast
/// (h[2p], h[2p + 1]) (LLVM's form of the i32 one-tap-per-pass body: `vpmulld`, dec 1.065).
#[cfg(all(target_arch = "x86_64", not(portable)))]
#[inline(always)]
pub(crate) fn hist_terms<const N: usize>(
    acc: &mut [i32; 8],
    c: &[[i16; 16]; HIST_PAIRS],
    h: &[i16; N],
) {
    if avx2() {
        // SAFETY: `hist_terms_avx2` only requires AVX2, which `avx2()` just confirmed.
        return unsafe { hist_terms_avx2::<N>(acc, c, h) };
    }
    crate::synthesis::hist_terms_portable::<N>(acc, c, h)
}

#[cfg(all(target_arch = "x86_64", not(portable)))]
#[target_feature(enable = "avx2")]
#[inline]
fn hist_terms_avx2<const N: usize>(acc: &mut [i32; 8], c: &[[i16; 16]; HIST_PAIRS], h: &[i16; N]) {
    use core::arch::x86_64::{
        __m256i, _mm256_add_epi32, _mm256_loadu_si256, _mm256_madd_epi16, _mm256_set1_epi32,
        _mm256_setzero_si256, _mm256_slli_epi32, _mm256_storeu_si256, _mm256_sub_epi32,
    };
    const { assert!(N.is_multiple_of(2)) };
    let c = &c[..N / 2];
    let mut t = _mm256_setzero_si256();
    for (cp, hp) in c.iter().zip(h.as_chunks::<2>().0) {
        // SAFETY: 16 i16 behind a reference; AVX2 is enabled in this function.
        let cv = unsafe { _mm256_loadu_si256(cp.as_ptr() as *const __m256i) };
        let hv = _mm256_set1_epi32(tap_pair(hp[0], hp[1]));
        t = _mm256_add_epi32(t, _mm256_madd_epi16(cv, hv));
    }
    // SAFETY: `acc` holds 8 i32 (32 bytes), read and written whole.
    unsafe {
        let a = _mm256_loadu_si256(acc.as_ptr() as *const __m256i);
        let r = _mm256_sub_epi32(a, _mm256_slli_epi32::<1>(t));
        _mm256_storeu_si256(acc.as_mut_ptr() as *mut __m256i, r);
    }
}

/// Two taps as one `pmaddwd` operand: c0 in the low 16 bits, c1 in the high.
#[cfg(all(target_arch = "x86_64", not(portable)))]
#[inline(always)]
fn tap_pair(c0: i16, c1: i16) -> i32 {
    (c0 as u16 as u32 | (c1 as u16 as u32) << 16) as i32
}

#[cfg(all(target_arch = "x86_64", not(portable)))]
#[target_feature(enable = "avx2")]
#[inline]
fn fir32_avx2<const NO: usize, const SH: i32>(
    x: &[i16],
    c: &[i16; L_FIR + 2],
    out: &mut [i16; NO],
) {
    use core::arch::x86_64::{
        __m256i, _mm256_add_epi32, _mm256_loadu_si256, _mm256_madd_epi16, _mm256_packs_epi32,
        _mm256_permute2x128_si256, _mm256_permute4x64_epi64, _mm256_set1_epi32,
        _mm256_setzero_si256, _mm256_srai_epi32, _mm256_storeu_si256, _mm256_unpackhi_epi16,
        _mm256_unpacklo_epi16,
    };
    // X2 needs i < NO + 30: NB2 = NO/16 + 2 blocks of 16 pairs, reading xp[..16·NB2 + 1]
    const NX: usize = 128;
    let nb2 = NO / 16 + 2; // <= 7: 16·7 + 1 <= NX, x2 holds 112
    let x = &x[..NO + L_FIR + 1]; // checked: the samples the outputs read
    let mut xp = [0i16; NX];
    xp[..x.len()].copy_from_slice(x);
    let mut x2 = [0i32; 16 * 7];
    let mut acc = [_mm256_setzero_si256(); 10]; // NO / 8 used
    // Loads of 16 i16 at xp[16b] and xp[16b + 1] (b < nb2 <= 7: within NX), stores of 8 i32
    // at x2[16b] and x2[16b + 8] (within 112); reads of 8 i32 at x2[8q + 2m] (q < NO/8,
    // m < 16: at most NO - 8 + 30 + 8 <= 110 <= 112); stores of 16 i16 at out[16j], j < NO/16.
    // SAFETY: every access is in bounds (above); AVX2 is enabled in this function.
    unsafe {
        let ld = |p: *const i16| _mm256_loadu_si256(p as *const __m256i);
        for b in 0..nb2 {
            let (a, a1) = (ld(xp.as_ptr().add(16 * b)), ld(xp.as_ptr().add(16 * b + 1)));
            // per 128-bit half: pairs 0..3 | 8..11 (lo), 4..7 | 12..15 (hi)
            let (lo, hi) = (_mm256_unpacklo_epi16(a, a1), _mm256_unpackhi_epi16(a, a1));
            let q = x2.as_mut_ptr().add(16 * b) as *mut __m256i;
            _mm256_storeu_si256(q, _mm256_permute2x128_si256::<0x20>(lo, hi));
            _mm256_storeu_si256(q.add(1), _mm256_permute2x128_si256::<0x31>(lo, hi));
        }
        for m in 0..16 {
            let cp = _mm256_set1_epi32(tap_pair(c[2 * m], c[2 * m + 1]));
            for (q, a) in acc[..NO / 8].iter_mut().enumerate() {
                let w = _mm256_loadu_si256(x2.as_ptr().add(8 * q + 2 * m) as *const __m256i);
                *a = _mm256_add_epi32(*a, _mm256_madd_epi16(w, cp));
            }
        }
        let half = _mm256_set1_epi32(1 << (SH - 1));
        let r = |a: __m256i| _mm256_srai_epi32::<SH>(_mm256_add_epi32(a, half));
        for j in 0..NO / 16 {
            // packs saturates and interleaves per 128-bit half: (0..3, 8..11 | 4..7, 12..15)
            let p = _mm256_packs_epi32(r(acc[2 * j]), r(acc[2 * j + 1]));
            let p = _mm256_permute4x64_epi64::<0xd8>(p);
            _mm256_storeu_si256(out.as_mut_ptr().add(16 * j) as *mut __m256i, p);
        }
    }
}

/// `synthesis::oversamp_interpolate_portable` (the contract and the oracle): the four
/// interpolated phases of `oversamp_12k8_to_16k`. On AVX2 an unaligned load of 16 samples at
/// offset b is the 8 sample pairs (s[b + 2k], s[b + 2k + 1]): output m of phase p at tap pair
/// j reads the pair at 4m + p + 2j, so one load serves phase p (even pairs) and p + 2 (odd
/// pairs) at the same j, and one `pmaddwd` with the two phases' tap pairs alternating does both
/// (LLVM's lanes widened every tap to i32 and multiplied with `vpmulld`).
#[inline(always)]
pub(crate) fn oversamp_interpolate(
    signal: &[i16; crate::synthesis::OVERSAMP_IN],
    out: &mut [i16; crate::cnst::L_SUBFR16K],
) {
    #[cfg(all(target_arch = "x86_64", not(portable)))]
    if avx2() {
        // SAFETY: `oversamp_avx2` only requires AVX2, which `avx2()` just confirmed.
        return unsafe { oversamp_avx2(signal, out) };
    }
    crate::synthesis::oversamp_interpolate_portable(signal, out)
}

#[cfg(all(target_arch = "x86_64", not(portable)))]
#[target_feature(enable = "avx2")]
#[inline]
fn oversamp_avx2(
    signal: &[i16; crate::synthesis::OVERSAMP_IN],
    out: &mut [i16; crate::cnst::L_SUBFR16K],
) {
    use crate::tables::FIR_UP;
    use core::arch::x86_64::{
        __m256i, _mm256_add_epi32, _mm256_loadu_si256, _mm256_madd_epi16, _mm256_max_epi32,
        _mm256_min_epi32, _mm256_set1_epi32, _mm256_setr_epi32, _mm256_srai_epi32,
        _mm256_storeu_si256,
    };
    // The 24 taps of phase p as 12 pairs (c[2j], c[2j + 1]) in one i32 each.
    let pair = |p: usize, j: usize| {
        tap_pair(
            FIR_UP[24 * (p - 1) + 2 * j],
            FIR_UP[24 * (p - 1) + 2 * j + 1],
        )
    };
    // the furthest load: group 1 (p = 2), j = 11, block 3, 16 samples
    const { assert!(2 + 2 * 11 + 16 * 3 + 16 <= crate::synthesis::OVERSAMP_IN) };
    let mut res = [[0i32; 32]; 2];
    // Group g: phases p = g + 1 (even lanes) and p + 2 (odd lanes); lane k of block q is
    // output m = 4q + k/2. Loads of 16 i16 at g + 1 + 2j + 16q (j < 12, q < 4): at most
    // 2 + 22 + 48 + 15 = 87 < OVERSAMP_IN = 88; stores of 8 i32 at res[g][8q], q < 4.
    // SAFETY: every access is in bounds (above); AVX2 is enabled in this function.
    unsafe {
        for (g, res) in res.iter_mut().enumerate() {
            let p = g + 1;
            let mut acc = [_mm256_set1_epi32(0x2000); 4];
            for j in 0..12 {
                let (a, b) = (pair(p, j), pair(p + 2, j));
                let c = _mm256_setr_epi32(a, b, a, b, a, b, a, b);
                let base = signal.as_ptr().add(p + 2 * j);
                for (q, acc) in acc.iter_mut().enumerate() {
                    let w = _mm256_loadu_si256(base.add(16 * q) as *const __m256i);
                    *acc = _mm256_add_epi32(*acc, _mm256_madd_epi16(w, c));
                }
            }
            // shl_int32(s, 2) >> 16 = clamp(s, -2^29, 2^29 - 1) >> 14
            let (lo, hi) = (
                _mm256_set1_epi32(-(1 << 29)),
                _mm256_set1_epi32((1 << 29) - 1),
            );
            for (q, acc) in acc.iter().enumerate() {
                let v = _mm256_srai_epi32::<14>(_mm256_min_epi32(_mm256_max_epi32(*acc, lo), hi));
                _mm256_storeu_si256(res.as_mut_ptr().add(8 * q) as *mut __m256i, v);
            }
        }
    }
    for (g, res) in res.iter().enumerate() {
        for (k, &v) in res.iter().enumerate() {
            out[5 * (k / 2) + g + 1 + 2 * (k % 2)] = v as i16;
        }
    }
}

/// `excitation::phase_conv_portable` (the contract and the oracle): `phase_dispersion`'s
/// convolution. On AVX2 a row is 4 × (`pmulhrsw` of 16 taps with the broadcast sample, then
/// `paddsw`): `pmulhrsw` is `mult_int16_r` but for (-32768)², which the taps exclude (a
/// `const` assertion there); LLVM's lanes do each product in i32 (`pmaddwd`, shifts, packs).
#[inline(always)]
pub(crate) fn phase_conv(code: &mut [i16; L_SUBFR], imp: &[i16; L_SUBFR]) {
    #[cfg(all(target_arch = "x86_64", not(portable)))]
    if avx2() {
        // SAFETY: `phase_conv_avx2` only requires AVX2, which `avx2()` just confirmed.
        return unsafe { phase_conv_avx2(code, imp) };
    }
    crate::excitation::phase_conv_portable(code, imp)
}

#[cfg(all(target_arch = "x86_64", not(portable)))]
#[target_feature(enable = "avx2")]
#[inline]
fn phase_conv_avx2(code: &mut [i16; L_SUBFR], imp: &[i16; L_SUBFR]) {
    use core::arch::x86_64::{
        __m256i, _mm256_adds_epi16, _mm256_loadu_si256, _mm256_mulhrs_epi16, _mm256_set1_epi16,
        _mm256_storeu_si256,
    };
    let mut code2 = [0i16; 2 * L_SUBFR];
    // Loads of 16 i16 at imp[16k] and code[16k] (k < 4), at code2[i + 16k] (i < 64: within
    // 128) and code2[16k], code2[64 + 16k]; stores at the same code2 and code places.
    // SAFETY: every access is in bounds (above); AVX2 is enabled in this function.
    unsafe {
        let ld = |p: *const i16| _mm256_loadu_si256(p as *const __m256i);
        let t = [0, 1, 2, 3].map(|k| ld(imp.as_ptr().add(16 * k)));
        for (i, &ci) in code.iter().enumerate() {
            if ci != 0 {
                let c = _mm256_set1_epi16(ci);
                for (k, &tk) in t.iter().enumerate() {
                    let p = code2.as_mut_ptr().add(i + 16 * k);
                    let v = _mm256_adds_epi16(ld(p), _mm256_mulhrs_epi16(c, tk));
                    _mm256_storeu_si256(p as *mut __m256i, v);
                }
            }
        }
        for k in 0..4 {
            let v = _mm256_adds_epi16(
                ld(code2.as_ptr().add(16 * k)),
                ld(code2.as_ptr().add(64 + 16 * k)),
            );
            _mm256_storeu_si256(code.as_mut_ptr().add(16 * k) as *mut __m256i, v);
        }
    }
}

/// `pitch::norm_lags_portable` (the contract and the oracle): `Norm_Corr`'s fast lags. On
/// AVX2 the 64-sample excf stays in 4 registers from lag to lag: per lag the range test, both
/// dots as `pmaddwd` (one shared reduction), then the update e'[i] = sat(mult(t, h[i]) +
/// e[i - 1]) from `pmulhw`/`pmullw` halves ((t h) >> 15 = (hi << 1) | (lo >>> 15); its one
/// overflow, t = h = -32768, is the only way to get -32768 and saturates to 32767) and
/// `paddsw` over e shifted one lane. The scalar form ran the dots and the update as
/// separate passes through memory.
#[inline(always)]
pub(crate) fn norm_lags(
    p: &crate::encoder::kernels::Lags,
    samples: &[i16],
    excf: &mut [i16; L_SUBFR],
    out: &mut [[i32; 2]],
) -> usize {
    #[cfg(all(target_arch = "x86_64", not(portable)))]
    if avx2() {
        // SAFETY: `norm_lags_avx2` only requires AVX2, which `avx2()` just confirmed.
        return unsafe { norm_lags_avx2(p, samples, excf, out) };
    }
    k::norm_lags_portable(p, samples, excf, out)
}

#[cfg(all(target_arch = "x86_64", not(portable)))]
#[target_feature(enable = "avx2")]
#[inline]
fn norm_lags_avx2(
    p: &crate::encoder::kernels::Lags,
    samples: &[i16],
    excf: &mut [i16; L_SUBFR],
    out: &mut [[i32; 2]],
) -> usize {
    use crate::encoder::kernels::LAGS_V_MAX;
    use core::arch::x86_64::{
        __m256i, _mm_add_epi32, _mm_cvtsi128_si32, _mm_shuffle_epi32, _mm256_add_epi32,
        _mm256_adds_epi16, _mm256_alignr_epi8, _mm256_castsi256_si128, _mm256_cmpeq_epi16,
        _mm256_cmpgt_epi16, _mm256_extracti128_si256, _mm256_loadu_si256, _mm256_madd_epi16,
        _mm256_movemask_epi8, _mm256_mulhi_epi16, _mm256_mullo_epi16, _mm256_or_si256,
        _mm256_permute2x128_si256, _mm256_set1_epi16, _mm256_setzero_si256, _mm256_slli_epi16,
        _mm256_srli_epi16, _mm256_storeu_si256, _mm256_unpackhi_epi32, _mm256_unpacklo_epi32,
        _mm256_xor_si256,
    };
    let ld = |a: &[i16; 16]| -> __m256i {
        // SAFETY: 16 i16 behind a reference; AVX2 is enabled in this function.
        unsafe { _mm256_loadu_si256(a.as_ptr() as *const __m256i) }
    };
    let blocks = |a: &[i16; L_SUBFR]| -> [__m256i; 4] {
        core::array::from_fn(|q| ld(a[16 * q..16 * q + 16].try_into().unwrap()))
    };
    let (x, h) = (blocks(p.xn), blocks(p.h));
    let mut e = blocks(excf);
    let (hi_max, lo_min) = (
        _mm256_set1_epi16(LAGS_V_MAX),
        _mm256_set1_epi16(-LAGS_V_MAX),
    );
    let min16 = _mm256_set1_epi16(i16::MIN);
    let n = out.len();
    let mut done = n;
    for (i, o) in out.iter_mut().enumerate() {
        let mut range = _mm256_setzero_si256();
        let (mut sx, mut sv) = (_mm256_setzero_si256(), _mm256_setzero_si256());
        for q in 0..4 {
            range = _mm256_or_si256(range, _mm256_cmpgt_epi16(e[q], hi_max));
            range = _mm256_or_si256(range, _mm256_cmpgt_epi16(lo_min, e[q]));
            sx = _mm256_add_epi32(sx, _mm256_madd_epi16(x[q], e[q]));
            sv = _mm256_add_epi32(sv, _mm256_madd_epi16(e[q], e[q]));
        }
        if _mm256_movemask_epi8(range) != 0 {
            done = i;
            break;
        }
        // (sx0+sx2, sv0+sv2, sx1+sx3, sv1+sv3) per half, the halves added: lanes 0, 1
        let t = _mm256_add_epi32(_mm256_unpacklo_epi32(sx, sv), _mm256_unpackhi_epi32(sx, sv));
        let t = _mm_add_epi32(_mm256_castsi256_si128(t), _mm256_extracti128_si256::<1>(t));
        let t = _mm_add_epi32(t, _mm_shuffle_epi32::<0xee>(t));
        let (sxe, see) = (
            _mm_cvtsi128_si32(t),
            _mm_cvtsi128_si32(_mm_shuffle_epi32::<0x55>(t)),
        );
        if see as i64 > p.see_max {
            done = i;
            break;
        }
        *o = [sxe, see];
        if i + 1 < n {
            let tv = _mm256_set1_epi16(samples[samples.len() - 1 - i]);
            let prev = e;
            for q in 0..4 {
                let (hi, lo) = (_mm256_mulhi_epi16(tv, h[q]), _mm256_mullo_epi16(tv, h[q]));
                let m = _mm256_or_si256(_mm256_slli_epi16::<1>(hi), _mm256_srli_epi16::<15>(lo));
                let m = _mm256_xor_si256(m, _mm256_cmpeq_epi16(m, min16)); // -32768 → 32767
                // e shifted one lane up: [prev[q-1][15], prev[q][0..15]] (zero into lane 0)
                let below = if q > 0 {
                    prev[q - 1]
                } else {
                    _mm256_setzero_si256()
                };
                let carry = _mm256_permute2x128_si256::<0x21>(below, prev[q]);
                let sh = _mm256_alignr_epi8::<14>(prev[q], carry);
                e[q] = _mm256_adds_epi16(m, sh);
            }
        }
    }
    for q in 0..4 {
        let d: &mut [i16; 16] = (&mut excf[16 * q..16 * q + 16]).try_into().unwrap();
        // SAFETY: 16 i16 behind a reference; AVX2 is enabled in this function.
        unsafe { _mm256_storeu_si256(d.as_mut_ptr() as *mut __m256i, e[q]) };
    }
    done
}
