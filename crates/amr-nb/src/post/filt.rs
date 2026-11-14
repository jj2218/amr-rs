//! `Syn_filt`, `Residu`, `Weight_Ai` (opencore `common/src/{syn_filt,residu,weight_a}.cpp`),
//! with TS 26.073's saturating operators: opencore accumulates in wrapping `Word32` arithmetic,
//! which differs where a partial sum saturates, and its `Syn_filt` has no overflow output (the
//! reference sets the global flag, which `Decoder_amr` acts on).

use super::M;
use crate::basic_op::{l_mac, l_msu, l_mult, l_shl, max_abs, round, sum_abs};

/// Longest `lg` any caller passes (the tests use up to a full frame).
const MAX_LG: usize = 160;

/// `Syn_filt`: LP synthesis filter 1/A(z), `a` in Q12 (`M + 1` coefficients). Returns whether
/// any of the reference's operators saturated (its global `Overflow`).
///
/// Computes `y[..lg]` from `x`; past outputs come from `mem` (`M` samples, oldest first). If
/// `update`, `mem` receives `y[lg - M..lg]`. (opencore's unrolled loop drops an odd last sample.)
/// Precondition (as in C): `lg >= M`. C callers that filter in place (`x == y`) or whose `mem`
/// aliases `x`/`y` get identical results by passing copies of `x` and `mem` (C reads each input
/// sample before writing the output at that index, and reads `mem` before writing anything).
pub(crate) fn syn_filt(
    a: &[i16],
    x: &[i16],
    y: &mut [i16],
    lg: usize,
    mem: &mut [i16],
    update: bool,
) -> bool {
    debug_assert!((M..=MAX_LG).contains(&lg));
    let a = &a[..=M];
    let len = lg;
    // The C's L_mult/L_msu chain saturates nowhere if |a0| max|x| + (sum |a_j|) max|y| < 2^30
    // (y over mem and the outputs). Checked after the plain pass with the outputs it produced:
    // up to the first saturation they are right, so if even their maximum keeps the bound,
    // there was none; and if no output left the 16-bit range (`syn_fast`), `round` and `L_shl`
    // didn't saturate either.
    // Every sample is 16-bit, so the coefficients alone often settle it (sum |a| < 2^15,
    // 11 values) without scanning x, y and mem.
    let (a0, sa) = ((a[0] as i64).abs(), sum_abs(&a[1..]));
    let ok = syn_fast(a, &x[..len], &mut y[..len], &mem[..M])
        && ((a0 + sa) << 15 < 1 << 30 || {
            let ym = max_abs(&y[..len]).max(max_abs(&mem[..M]));
            a0 * max_abs(&x[..len]) + sa * ym < 1 << 30
        });
    let ov = if ok {
        false
    } else {
        syn_exact(a, &x[..len], &mut y[..len], &mem[..M])
    };
    if update {
        mem[..M].copy_from_slice(&y[lg - M..lg]);
    }
    ov
}

/// `syn_filt` without saturation on the critical path: y[n] = s >> 12 unclamped, the
/// out-of-range test (`syn_out_of_range`) only collected in a flag. Returns false if any
/// sample would have saturated (then `y` is garbage and `syn_exact` redoes it). Blocks of 8
/// outputs (`simd::syn_blocks`): the x term and the terms on the 10 outputs before the block
/// with k >= 3 are one vector pass, then the outputs run serially with only the y[n-1], y[n-2] terms (k = 1, 2)
/// on the chain; each output's terms in the block's later outputs (k >= 3) are one vector
/// update off the chain. 16 scalar `imul`s per block instead of 44.
fn syn_fast(a: &[i16], x: &[i16], y: &mut [i16], mem: &[i16]) -> bool {
    const B: usize = 8;
    let t = syn_tab(a.try_into().unwrap());
    // old[..M]: the M outputs before the first block, padded for whole-vector loads.
    let mut old = [0i16; 2 * B];
    old[..M].copy_from_slice(&mem[..M]);
    // Whole blocks in place; a partial last block through zero-padded x (a false alarm there
    // only costs the exact pass).
    let (xc, xr) = x.as_chunks::<B>();
    let (yc, yr) = y.as_chunks_mut::<B>();
    let mut bad = crate::simd::syn_blocks(&t, &old, xc, yc);
    if !xr.is_empty() {
        // the M outputs before the last block: the end of mem ++ the whole blocks' outputs
        let mut h = [0i16; M + MAX_LG];
        let n = xc.len() * B;
        h[..M].copy_from_slice(&mem[..M]);
        h[M..M + n].copy_from_slice(yc.as_flattened());
        old[..M].copy_from_slice(&h[n..n + M]);
        let mut xt = [[0i16; B]];
        xt[0][..xr.len()].copy_from_slice(xr);
        let mut yt = [[0i16; B]];
        bad |= crate::simd::syn_blocks(&t, &old, &xt, &mut yt);
        yr.copy_from_slice(&yt[0][..yr.len()]);
    }
    bad == 0
}

/// `syn_fast`'s coefficient tables for `a` (see `SynTab`).
fn syn_tab(a: &[i16; M + 1]) -> SynTab {
    // A pair of i16 in one lane, as `pmaddwd` reads it.
    let pair = |lo: i16, hi: i16| (lo as u16 as u32 | (hi as u16 as u32) << 16) as i32;
    // col(j, i): coefficient of old[j] = y[n0 - M + j] in output n0 + i, for k = i + M - j in
    // 3..=M (k = 1, 2 run on the chain).
    let col = |j: usize, i: usize| {
        let k = i + M - j;
        if (3..=M).contains(&k) { a[k] } else { 0 }
    };
    SynTab {
        colp: core::array::from_fn(|q| {
            core::array::from_fn(|i| pair(col(2 * q, i), col(2 * q + 1, i)))
        }),
        upd: core::array::from_fn(|m| {
            core::array::from_fn(|i| if i >= m + 3 { pair(a[i - m], 0) } else { 0 })
        }),
        a0x: pair(a[0], 0x800),
        a1: a[1] as i32,
        a2: a[2] as i32,
    }
}

/// One `syn_fast` call's coefficients as `simd::syn_blocks` takes them, each lane a pair of i16 (low,
/// high) for `pmaddwd`: `colp[q][i]` the coefficients of old[2q], old[2q + 1] in output i
/// (0 where k is outside 3..=M), `upd[m][i]` output m's in output i (k = i - m >= 3), `a0x`
/// (a0, 0x800) against (x, 1).
pub(crate) struct SynTab {
    pub(crate) colp: [[i32; 8]; 5],
    pub(crate) upd: [[i32; 8]; 5],
    pub(crate) a0x: i32,
    pub(crate) a1: i32,
    pub(crate) a2: i32,
}

/// Flag-free only for s in [-2^27 + 2048, 2^27 - 2]: with the C's s' = s - 0x800,
/// `L_shl(2 s', 3)` saturates (and flags) below -2^27 + 2048 although the output is the same
/// -32768; the top end is conservative (exact pass). 1 if `s` is outside.
#[inline(always)]
pub(crate) fn syn_out_of_range(s: i32) -> u32 {
    (s.wrapping_add(134_217_728 - 2048) as u32 >= 0x0fff_ffff - 2048) as u32
}

/// `simd::syn_blocks`' array form (its non-x86 body, and its oracle): blocks of 8 outputs
/// from `old[..M]` (the M outputs before the first block) and `x`, wrapping i32 sums as
/// `pmaddwd` makes them. Outputs as i16 (exact unless the flag is set); returns the OR of
/// `syn_out_of_range`.
#[cfg_attr(all(target_arch = "x86_64", not(amr_nb_arrays)), allow(dead_code))]
pub(crate) fn syn_blocks_portable(
    t: &SynTab,
    old: &[i16; 16],
    x: &[[i16; 8]],
    y: &mut [[i16; 8]],
) -> u32 {
    let (lo, hi) = (|w: i32| w as i16 as i32, |w: i32| (w >> 16) as i16 as i32);
    let mut old: [i16; M] = old[..M].try_into().unwrap();
    let mut bad = 0u32;
    for (xb, yo) in x.iter().zip(y) {
        let mut acc: [i32; 8] =
            core::array::from_fn(|i| (xb[i] as i32 * lo(t.a0x)).wrapping_add(hi(t.a0x)));
        for (q, cq) in t.colp.iter().enumerate() {
            let (o0, o1) = (old[2 * q] as i32, old[2 * q + 1] as i32);
            for (s, &c) in acc.iter_mut().zip(cq) {
                *s = s.wrapping_sub((lo(c) * o0).wrapping_add(hi(c) * o1));
            }
        }
        let (mut p2, mut p1) = (old[M - 2] as i32, old[M - 1] as i32);
        for m in 0..8 {
            let s = acc[m]
                .wrapping_sub(t.a2.wrapping_mul(p2))
                .wrapping_sub(t.a1.wrapping_mul(p1));
            bad |= syn_out_of_range(s);
            let v = (s >> 12).clamp(i16::MIN as i32, i16::MAX as i32) as i16;
            yo[m] = v;
            if m < 5 {
                for (s, &c) in acc.iter_mut().zip(&t.upd[m]) {
                    *s = s.wrapping_sub(lo(c) * v as i32);
                }
            }
            (p2, p1) = (p1, v as i32);
        }
        old.copy_within(8.., 0);
        old[M - 8..].copy_from_slice(yo);
    }
    bad
}

/// `syn_filt`'s exact loop: the C's chain `L_mult`, `L_msu` (j = 1..=M), `L_shl(s, 3)`,
/// `round`, saturating. Returns whether anything saturated.
#[cold]
fn syn_exact(a: &[i16], x: &[i16], y: &mut [i16], mem: &[i16]) -> bool {
    let mut ov = false;
    let mut buf = [0i16; M + MAX_LG];
    buf[..M].copy_from_slice(&mem[..M]);
    for (n, &xn) in x.iter().enumerate() {
        let mut s = l_mult(xn, a[0], &mut ov);
        for j in 1..=M {
            s = l_msu(s, a[j], buf[M + n - j], &mut ov);
        }
        let t = l_shl(s, 3, &mut ov);
        ov |= t != s << 3; // L_shl saturated (our l_shl doesn't flag)
        buf[M + n] = round(t, &mut ov);
    }
    y.copy_from_slice(&buf[M..M + x.len()]);
    ov
}

/// `Residu`: LP analysis filter A(z), `a` in Q12 (`M + 1` coefficients).
///
/// `x` starts `M` samples before C's `input_ptr` (so `x[M + n]` is C `input_ptr[n]`) and must
/// hold `M + lg` samples. All `lg` outputs (opencore's unrolled loop skips the first `lg % 4`).
pub(crate) fn residu(a: &[i16], x: &[i16], y: &mut [i16], lg: usize) {
    let a = &a[..=M];
    let first = 0;
    // TS 26.073 saturates the L_mac chain, `L_shl(s, 3)` and `round` (opencore wraps and
    // truncates: wrong on full-scale input). If no partial sum can reach 2^30 (2^31 after
    // L_mult's doubling), that is exactly clamp((sum + 0x800) >> 12) on the exact sum.
    if sum_abs(a) * max_abs(&x[first..M + lg]) >= 1 << 30 {
        let ov = &mut false;
        for n in first..lg {
            let mut s = l_mult(x[M + n], a[0], ov);
            for j in 1..=M {
                s = l_mac(s, a[j], x[M + n - j], ov);
            }
            y[n] = round(l_shl(s, 3, ov), ov);
        }
        return;
    }
    // Exact i32 sums (|sum| < 2^30): the C's tap order doesn't matter. Blocks of outputs: 8
    // wide (`simd::residu_block`), then 4, then single outputs.
    let mut n = first;
    let a11: &[i16; M + 1] = a.try_into().unwrap();
    while n + 8 <= lg {
        let w: &[i16; M + 8] = x[n..n + M + 8].try_into().unwrap();
        store_q12(&crate::simd::residu_block(a11, w), &mut y[n..n + 8]);
        n += 8;
    }
    if n + 4 <= lg {
        block::<4>(a, x, y, n);
        n += 4;
    }
    while n < lg {
        block::<1>(a, x, y, n);
        n += 1;
    }
}

/// `residu`'s tails: blocks of B = 4 or 1 outputs n..n + B (taps outer, so the lanes
/// vectorise; 8 wide is `simd::residu_block`).
fn block<const B: usize>(a: &[i16], x: &[i16], y: &mut [i16], n: usize) {
    store_q12(&block_acc::<B>(a, x, n), &mut y[n..n + B]);
}

/// `residu`'s outputs from its sums: clamp(s >> 12).
fn store_q12(acc: &[i32], y: &mut [i16]) {
    for (o, &s) in y.iter_mut().zip(acc) {
        *o = (s >> 12).clamp(i16::MIN as i32, i16::MAX as i32) as i16;
    }
}

/// `block`'s sums 0x800 + sum_k a[k] x[M + n + i - k] (wrapping).
fn block_acc<const B: usize>(a: &[i16], x: &[i16], n: usize) -> [i32; B] {
    let mut acc = [0x800i32; B];
    for (k, &c) in a.iter().enumerate() {
        let xs: &[i16; B] = x[M + n - k..][..B].try_into().unwrap();
        for (s, &v) in acc.iter_mut().zip(xs) {
            *s = s.wrapping_add(c as i32 * v as i32);
        }
    }
    acc
}

/// `simd::residu_block`'s array form over its window `w` = x[n..n + M + 8].
#[cfg_attr(all(target_arch = "x86_64", not(amr_nb_arrays)), allow(dead_code))]
pub(crate) fn residu_block_portable(a: &[i16; M + 1], w: &[i16; M + 8]) -> [i32; 8] {
    block_acc::<8>(a, w, 0)
}

/// `Weight_Ai`: spectral expansion `a_exp[i] = round(L_mult(a[i], fac[i - 1]))`, `a_exp[0] =
/// a[0]`, saturating as TS 26.073 (opencore truncates: -32768 * -32768 wraps).
pub(crate) fn weight_ai(a: &[i16], fac: &[i16], a_exp: &mut [i16]) {
    let ov = &mut false;
    a_exp[0] = a[0];
    for i in 1..=M {
        a_exp[i] = round(l_mult(a[i], fac[i - 1], ov), ov);
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_c {
    use super::*;
    use crate::post::test_gen::{coefs, lpc, signal};
    use crate::testutil::{N, Rng};
    use opencore_sys as _;

    // Against the TS 26.073 reference (opencore wraps and skips odd/remainder samples).
    unsafe extern "C" {
        #[link_name = "ref_Syn_filt"]
        fn ref_syn_filt(a: *mut i16, x: *mut i16, y: *mut i16, lg: i16, mem: *mut i16, upd: i16);
        #[link_name = "ref_Residu"]
        fn ref_residu(a: *mut i16, x: *mut i16, y: *mut i16, lg: i16);
        #[link_name = "ref_Weight_Ai"]
        fn ref_weight_ai(a: *mut i16, fac: *const i16, a_exp: *mut i16);
    }

    #[allow(non_snake_case)]
    unsafe fn Weight_Ai(a: *mut i16, fac: *const i16, a_exp: *mut i16) {
        let _l = ref26073_sys::lock();
        unsafe { ref_weight_ai(a, fac, a_exp) }
    }

    /// The reference's `Syn_filt` (holds the reference lock; returns its overflow flag).
    #[allow(non_snake_case)]
    unsafe fn Syn_filt(
        a: *mut i16,
        x: *mut i16,
        y: *mut i16,
        lg: i16,
        mem: *mut i16,
        upd: i16,
    ) -> bool {
        let _l = ref26073_sys::lock();
        ref26073_sys::clear_overflow();
        unsafe { ref_syn_filt(a, x, y, lg, mem, upd) };
        ref26073_sys::overflow()
    }

    #[allow(non_snake_case)]
    unsafe fn Residu(a: *mut i16, x: *mut i16, y: *mut i16, lg: i16) {
        let _l = ref26073_sys::lock();
        unsafe { ref_residu(a, x, y, lg) }
    }

    fn check_syn_filt(
        a: [i16; M + 1],
        x: [i16; 160],
        y0: [i16; 160],
        mem0: [i16; M],
        lg: usize,
        update: bool,
    ) {
        let (mut y, mut mem) = (y0, mem0);
        let ov = syn_filt(&a, &x, &mut y, lg, &mut mem, update);
        let (mut ac, mut xc, mut yc, mut memc) = (a, x, y0, mem0);
        let ovc = unsafe {
            Syn_filt(
                ac.as_mut_ptr(),
                xc.as_mut_ptr(),
                yc.as_mut_ptr(),
                lg as i16,
                memc.as_mut_ptr(),
                update as i16,
            )
        };
        // The overflow flag too: the decoder re-synthesises on it.
        assert_eq!(
            (y, mem, ov),
            (yc, memc, ovc),
            "Syn_filt({a:?}, {x:?}, {mem0:?}, lg={lg})"
        );
    }

    /// `simd::syn_blocks` (SSE2 on x86_64) is its array form `syn_blocks_portable`: same
    /// range flag always, same outputs when it is clear; any coefficients (the `-32768` pair
    /// products too), full-scale history and input, 1 to 20 blocks.
    #[test]
    fn syn_blocks_is_the_array_form() {
        let mut g = Rng::new(11);
        let mut clear = 0;
        for case in 0..N {
            let mut a = coefs(&mut g);
            if case % 8 == 0 {
                a = [-32768; M + 1];
            }
            let t = syn_tab(&a);
            let mut old = [0i16; 16];
            signal(&mut g, &mut old[..M]);
            let nb = 1 + (g.u64() % 20) as usize;
            let mut x = [[0i16; 8]; 20];
            signal(&mut g, x.as_flattened_mut());
            let (mut y, mut yp) = ([[0i16; 8]; 20], [[0i16; 8]; 20]);
            let bad = crate::simd::syn_blocks(&t, &old, &x[..nb], &mut y[..nb]);
            let badp = syn_blocks_portable(&t, &old, &x[..nb], &mut yp[..nb]);
            assert_eq!(bad, badp, "{a:?} {old:?}");
            if bad == 0 {
                clear += 1;
                assert_eq!(y, yp, "{a:?} {old:?}");
            }
        }
        assert!(clear > N / 4, "{clear} of {N} cases without the flag");
    }

    /// `syn_blocks_portable` (the non-x86 body, run on x86 by this test only) is `syn_exact`
    /// where neither flags anything: encoder-like coefficients, 1 to 20 blocks.
    #[test]
    fn syn_blocks_portable_is_the_exact_loop() {
        let mut g = Rng::new(12);
        let mut clear = 0;
        for _ in 0..N {
            let a = lpc(&mut g);
            let mut old = [0i16; 16];
            for v in &mut old[..M] {
                *v = (g.u64() as i16) >> (2 + g.u64() % 8);
            }
            let nb = 1 + (g.u64() % 20) as usize;
            let mut x = [[0i16; 8]; 20];
            for v in x.as_flattened_mut() {
                *v = (g.u64() as i16) >> (4 + g.u64() % 8);
            }
            let mut y = [[0i16; 8]; 20];
            let bad = syn_blocks_portable(&syn_tab(&a), &old, &x[..nb], &mut y[..nb]);
            let mut ye = [0i16; 8 * 20];
            let ov = syn_exact(&a, x[..nb].as_flattened(), &mut ye[..8 * nb], &old[..M]);
            if bad == 0 && !ov {
                clear += 1;
                assert_eq!(y[..nb].as_flattened(), &ye[..8 * nb], "{a:?} {old:?}");
            }
        }
        assert!(clear > N / 2, "{clear} of {N} cases without a flag");
    }

    /// `simd::residu_block` is its array form on any coefficients and window, full scale too
    /// (on x86_64; elsewhere the two are one function and `residu_vs_c` is the oracle).
    #[test]
    fn residu_block_is_the_array_form() {
        let mut g = Rng::new(13);
        for case in 0..N {
            let (mut a, mut w) = ([0i16; M + 1], [0i16; M + 8]);
            if case % 8 == 0 {
                (a, w) = ([-32768; M + 1], [-32768; M + 8]);
            } else {
                signal(&mut g, &mut a);
                signal(&mut g, &mut w);
            }
            let r = crate::simd::residu_block(&a, &w);
            assert_eq!(r, residu_block_portable(&a, &w), "{a:?} {w:?}");
        }
    }

    #[test]
    fn post_syn_filt_vs_c() {
        let mut g = Rng::new(1);
        // The reference's Syn_filt holds lg + M samples in a 80-sample buffer: lg <= 70.
        let lgs = [10usize, 11, 12, 22, 40, 41, 69, 70];
        for case in 0..N / 4 {
            let lg = lgs[case % lgs.len()];
            let a = coefs(&mut g);
            let mut x = [0i16; 160];
            signal(&mut g, &mut x);
            let mut y0 = [0i16; 160];
            signal(&mut g, &mut y0);
            let mut mem0 = [0i16; M];
            signal(&mut g, &mut mem0);
            check_syn_filt(a, x, y0, mem0, lg, g.u64() & 1 == 0);
        }
    }

    /// A chain that saturates mid-way while every output stays 16-bit, so `syn_fast` passes and
    /// only the bound can send it to the exact loop: a0 = 32767, a1 = -32768, a2 = 32767, a
    /// history near `level`, and x chosen per sample to hold the output there (the second
    /// partial sum then passes 2^31 for the larger levels).
    #[test]
    fn post_syn_filt_saturating_chain_vs_c() {
        let mut g = Rng::new(5);
        let mut a = [0i16; M + 1];
        (a[0], a[1], a[2]) = (32767, -32768, 32767);
        for case in 0..2000 {
            let lg = [16usize, 24, 40, 64][case % 4]; // whole blocks: syn_fast pads with x = 0
            let level = 10_000 + (g.u64() % 22_000) as i32;
            let mem0 = [level as i16; M];
            let mut y = [0i32; M + 160];
            y[..M].fill(level);
            let mut x = [0i16; 160];
            for n in 0..lg {
                let s: i32 = (1..=M).map(|j| a[j] as i32 * y[M + n - j]).sum();
                let t = level + (g.u64() % 64) as i32 - 32;
                x[n] = ((t * 4096 + s) / a[0] as i32) as i16;
                y[M + n] = (x[n] as i32 * a[0] as i32 - s + 2048) >> 12;
            }
            let mut yf = [0i16; 160];
            assert!(syn_fast(&a, &x[..lg], &mut yf[..lg], &mem0), "case {case}");
            check_syn_filt(a, x, [0; 160], mem0, lg, g.u64() & 1 == 0);
        }
    }

    /// Stable filters (the codec's kind) with outputs up to full scale: the unsaturated fast
    /// path decides most cases, the exact fallback the rest (both counted, so a change of
    /// the generator can't silently drop either).
    #[test]
    fn post_syn_filt_stable_vs_c() {
        let mut g = Rng::new(3);
        let lgs = [10usize, 11, 16, 22, 40, 41, 69, 70]; // lg <= 70, as above
        let (mut fast, mut exact) = (0, 0);
        for case in 0..N / 4 {
            let lg = lgs[case % lgs.len()];
            let a = lpc(&mut g);
            let shift = (g.u64() % 8) as u32;
            let x: [i16; 160] = core::array::from_fn(|_| (g.u64() as i16) >> shift);
            let mut mem0 = [0i16; M];
            signal(&mut g, &mut mem0);
            let mut y = [0i16; 160];
            if syn_fast(&a, &x[..lg], &mut y[..lg], &mem0) {
                fast += 1;
            } else {
                exact += 1;
            }
            check_syn_filt(a, x, [0; 160], mem0, lg, g.u64() & 1 == 0);
        }
        assert!(
            fast > N / 16 && exact > N / 64,
            "fast {fast}, exact {exact}"
        );
    }

    /// The saturation boundary: an accumulator of exactly `0x07ff_ffff` (opencore's off-by-one
    /// gave `MIN_16` there; the reference 32767). With a = [4096, 0, -1, 0..], x = mem = 2047 the
    /// output is a fixed point at 2047, and `x[p] = 32767` makes the accumulator hit the
    /// boundary exactly at sample `p` (in both the `mem` and the `y` history loops).
    #[test]
    fn post_syn_filt_boundary_vs_c() {
        let mut a = [0i16; M + 1];
        a[0] = 4096;
        a[2] = -1;
        for p in [0usize, 1, 8, 9, 10, 11, 38, 39] {
            let mut x = [2047i16; 160];
            x[p] = 32767;
            check_syn_filt(a, x, [0; 160], [2047; M], 40, true);
        }
    }

    #[test]
    fn post_residu_vs_c() {
        let mut g = Rng::new(2);
        let lgs = [0usize, 4, 5, 7, 40, 41, 43, 160];
        for case in 0..N / 4 {
            let lg = lgs[case % lgs.len()];
            let a = coefs(&mut g);
            let mut x = [0i16; M + 160];
            signal(&mut g, &mut x);
            let mut y0 = [0i16; 160];
            signal(&mut g, &mut y0);
            let mut y = y0;
            residu(&a, &x, &mut y, lg);
            let (mut ac, mut yc) = (a, y0);
            unsafe {
                Residu(
                    ac.as_mut_ptr(),
                    x[M..].as_mut_ptr(),
                    yc.as_mut_ptr(),
                    lg as i16,
                )
            };
            assert_eq!(y, yc, "case {case}: Residu({a:?}, lg={lg})");
        }
    }

    #[test]
    fn post_weight_ai_vs_c() {
        let mut g = Rng::new(3);
        for case in 0..N {
            let mut a = [0i16; M + 1];
            let mut fac = [0i16; M];
            a.iter_mut().for_each(|v| *v = g.i16());
            fac.iter_mut().for_each(|v| *v = g.i16());
            let mut r = [0i16; M + 1];
            weight_ai(&a, &fac, &mut r);
            let mut c = [0i16; M + 1];
            unsafe { Weight_Ai(a.as_mut_ptr(), fac.as_ptr(), c.as_mut_ptr()) };
            assert_eq!(r, c, "case {case}: Weight_Ai({a:?}, {fac:?})");
        }
    }
}
