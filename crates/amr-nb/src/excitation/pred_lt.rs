//! Long-term (adaptive codebook) prediction with fractional lag (`common/pred_lt.cpp`).

use crate::basic_op::{Flag, MAX_32, l_mac, max_abs, round};

const UP_SAMP_MAX: i16 = 6;
const L_INTER10: usize = super::L_INTERPOL - 1;
const FIR_SIZE: usize = UP_SAMP_MAX as usize * L_INTER10 + 1;
/// Outputs per block, and the input window a block reads.
const B: usize = 8;
pub(crate) const BLOCK_IN: usize = B + 2 * L_INTER10 - 1;

/// `inter_6_pred_lt`: 1/6-resolution interpolation filter.
const INTER_6_PRED_LT: [i16; FIR_SIZE] = [
    29443, 28346, 25207, 20449, 14701, 8693, 3143, -1352, -4402, -5865, -5850, -4673, -2783, -672,
    1211, 2536, 3130, 2991, 2259, 1170, 0, -1001, -1652, -1868, -1666, -1147, -464, 218, 756, 1060,
    1099, 904, 550, 135, -245, -514, -634, -602, -451, -231, 0, 191, 308, 340, 296, 198, 78, -36,
    -120, -163, -165, -132, -79, -19, 34, 73, 91, 89, 70, 38, 0,
];

/// The 20 taps per phase `frac` (0..=6, the indices the C's loop can form): `k[9 - i]` =
/// c1[6i] (c1 = table + frac), `k[10 + i]` = c2[6i] (c2 = table + 6 - frac).
const PHASE_TAPS: [[i16; 2 * L_INTER10]; UP_SAMP_MAX as usize + 1] = {
    let mut t = [[0i16; 2 * L_INTER10]; UP_SAMP_MAX as usize + 1];
    let mut f = 0;
    while f <= UP_SAMP_MAX as usize {
        let mut i = 0;
        while i < L_INTER10 {
            t[f][L_INTER10 - 1 - i] = INTER_6_PRED_LT[f + UP_SAMP_MAX as usize * i];
            t[f][L_INTER10 + i] =
                INTER_6_PRED_LT[UP_SAMP_MAX as usize - f + UP_SAMP_MAX as usize * i];
            i += 1;
        }
        f += 1;
    }
    t
};

/// Sum of |tap| per phase: a window with 2 max|x| SUM_ABS + 0x8000 <= MAX_32 can't saturate.
const SUM_ABS: [i64; UP_SAMP_MAX as usize + 1] = {
    let mut t = [0i64; UP_SAMP_MAX as usize + 1];
    let mut f = 0;
    while f < t.len() {
        let mut i = 0;
        while i < 2 * L_INTER10 {
            let c = PHASE_TAPS[f][i] as i64;
            t[f] += if c < 0 { -c } else { c };
            i += 1;
        }
        f += 1;
    }
    t
};

/// Whether the 20-tap sums over `w` can't saturate (then the plain sums are exact).
#[inline(always)]
fn safe(w: &[i16], sum_abs: i64) -> bool {
    2 * max_abs(w) * sum_abs + 0x8000 <= MAX_32 as i64
}

/// One output exactly as TS 26.073: `L_mac` chain in the C's order (x0[-i] c1, x1[i] c2), `round`.
#[cold]
fn exact(w: &[i16], k: &[i16; 2 * L_INTER10]) -> i16 {
    let ov = &mut false;
    let mut s = 0;
    for i in 0..L_INTER10 {
        s = l_mac(s, w[L_INTER10 - 1 - i], k[L_INTER10 - 1 - i], ov);
        s = l_mac(s, w[L_INTER10 + i], k[L_INTER10 + i], ov);
    }
    round(s, ov)
}

/// `Pred_lt_3or6`. C's `exc` pointer is `exc_buf[exc_off..]`; the function reads back to
/// `exc_off - T0 - L_INTERPOL`, so the buffer needs that much history (the decoder's
/// `old_exc` with `exc_off = PIT_MAX + L_INTERPOL`). Overwrites `l_subfr` (rounded down to
/// even) samples starting at `exc_off`. TS 26.073 saturates the `L_mac` chain and `round`
/// (opencore wraps): plain 32-bit sums where a window can't saturate, the exact chain elsewhere.
pub(crate) fn pred_lt_3or6(
    exc_buf: &mut [i16],
    exc_off: usize,
    t0: i16,
    mut frac: i16,
    l_subfr: i16,
    flag3: i16,
    _ov: &mut Flag,
) {
    // Index (into exc_buf) of C's pX0.
    let mut x0 = exc_off as isize - t0 as isize;
    frac = frac.wrapping_neg();
    if flag3 != 0 {
        frac <<= 1;
    }
    if frac < 0 {
        frac += UP_SAMP_MAX;
        x0 -= 1;
    }
    // Sample n is a 20-tap FIR over exc_buf[b-9..=b+10] with b = x0 + n: the left half uses
    // the c1 phase taps (reversed), the right half the c2 phase taps. The C accumulates with
    // wrapping (not saturating) adds, which are associative, so one dot product in any order
    // is bit-exact and vectorises. Each sample reads only earlier samples (T0 >= PIT_MIN > 10).
    let k = &PHASE_TAPS[frac as usize];
    let sum_abs = SUM_ABS[frac as usize];
    let base = x0 as usize - (L_INTER10 - 1);
    let len = l_subfr as usize;
    let sample = |w: &[i16]| {
        w[..2 * L_INTER10]
            .iter()
            .zip(k)
            .fold(0x0000_4000i32, |s, (&x, &c)| {
                s.wrapping_add(x as i32 * c as i32)
            })
    };
    // Blocks of B outputs (`simd::pred_block`), read in place: the inputs end before the outputs.
    // Output n reads up to exc_buf[x0 + n + 10] and x0 >= exc_off - t0 - 1, so a block only
    // reads finished samples if t0 >= B + 10; otherwise (out-of-range lag) go per sample.
    let mut n = 0;
    // With t0 >= len + 10 no output is read back, so one scan of the whole input decides
    // every block at once (the codec's usual case); shorter lags check per block.
    let all_safe = t0 as usize >= len + L_INTER10
        && safe(&exc_buf[base..base + len + 2 * L_INTER10 - 1], sum_abs);
    if t0 as usize >= B + L_INTER10 {
        while n + B <= len {
            // Read in place (a stack copy made each block's 20 unaligned window loads miss
            // store forwarding): the inputs end before exc_off + n, the outputs start there.
            let (past, out) = exc_buf.split_at_mut(exc_off + n);
            let src: &[i16; BLOCK_IN] = past[base + n..][..BLOCK_IN].try_into().unwrap();
            if all_safe || safe(src, sum_abs) {
                let acc = crate::simd::pred_block(src, k);
                for (o, a) in out[..B].iter_mut().zip(acc) {
                    *o = (a >> 15) as i16;
                }
            } else {
                for (j, o) in out[..B].iter_mut().enumerate() {
                    *o = exact(&src[j..], k);
                }
            }
            n += B;
        }
    }
    // One at a time, each written before the next is read (TS 26.073; opencore unrolls by two,
    // reading both samples first, and drops an odd last one: different for lags < 12).
    let one = |w: &[i16]| {
        if safe(&w[..2 * L_INTER10], sum_abs) {
            (sample(w) >> 15) as i16
        } else {
            exact(w, k)
        }
    };
    for n in n..len {
        exc_buf[exc_off + n] = one(&exc_buf[base + n..]);
    }
}

/// One block of 8 outputs (wrapping sums, so any tap order): `simd::pred_block`'s array form
/// (its non-x86 body and its oracle).
#[cfg_attr(all(target_arch = "x86_64", not(amr_nb_arrays)), allow(dead_code))]
pub(crate) fn block_portable(src: &[i16; BLOCK_IN], k: &[i16; 2 * L_INTER10]) -> [i32; B] {
    let mut acc = [0x0000_4000i32; B];
    for (t, &c) in k.iter().enumerate() {
        let x: &[i16; B] = src[t..][..B].try_into().unwrap();
        for (a, &v) in acc.iter_mut().zip(x) {
            *a = a.wrapping_add(v as i32 * c as i32);
        }
    }
    acc
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::super::{L_INTERPOL, L_SUBFR, PIT_MAX};
    use super::*;
    use crate::testutil::{N, Rng};
    use opencore_sys as _;

    unsafe extern "C" {
        #[link_name = "ref_Pred_lt_3or6"]
        fn c_pred_lt_3or6(exc: *mut i16, t0: i16, frac: i16, l_subfr: i16, flag3: i16);
    }

    const OFF: usize = PIT_MAX + L_INTERPOL;
    const LEN: usize = L_SUBFR + PIT_MAX + L_INTERPOL + 1; // + 1 for L_subfr = 41

    /// `simd::pred_block` is `block_portable` on any window and phase taps, full scale too
    /// (all -32768: `pmaddwd`'s pair sum 2^31 wraps like the wrapping sums).
    #[test]
    fn pred_block_is_the_array_form() {
        let mut g = Rng::new(3);
        for case in 0..N {
            let src: [i16; BLOCK_IN] = if case % 8 == 0 {
                [-32768; BLOCK_IN] // with the -32768 taps: each pmaddwd pair wraps to -2^31
            } else {
                core::array::from_fn(|_| g.i16())
            };
            let k = if case % 4 == 0 {
                [-32768; 2 * L_INTER10]
            } else {
                PHASE_TAPS[case % PHASE_TAPS.len()]
            };
            assert_eq!(crate::simd::pred_block(&src, &k), block_portable(&src, &k));
        }
    }

    #[test]
    fn pred_lt_3or6_vs_ref() {
        let mut g = Rng::new(20);
        for n in 0..N / 4 {
            // Mostly realistic excitation amplitudes, sometimes full-range (wrapping sums).
            let amp = if n % 4 == 0 { 0xffff } else { 0x1fff };
            let buf: [i16; LEN] = core::array::from_fn(|_| {
                if amp == 0xffff {
                    g.i16()
                } else {
                    (g.u64() & amp) as i16 - (amp as i16 >> 1)
                }
            });
            let flag3 = (g.u64() % 2) as i16;
            // Decoder ranges: T0 in PIT_MIN_MR122..=PIT_MAX; frac -1..=1 (1/3 res) or
            // -2..=3 (1/6 res, Dec_lag6).
            // Plus some out-of-range short lags (10..18: still in bounds here) for the
            // per-sample path, where outputs feed later outputs of the same block.
            let t0 = if n % 8 == 1 {
                10 + (g.u64() % 8) as i16
            } else {
                18 + (g.u64() % 126) as i16
            };
            let frac = if flag3 != 0 {
                (g.u64() % 3) as i16 - 1
            } else {
                (g.u64() % 6) as i16 - 2
            };
            let l_subfr = if n % 16 == 0 { 41 } else { 40 };
            let mut r = buf;
            let mut c = buf;
            let mut ro = false;
            pred_lt_3or6(&mut r, OFF, t0, frac, l_subfr, flag3, &mut ro);
            {
                let _l = ref26073_sys::lock();
                unsafe { c_pred_lt_3or6(c.as_mut_ptr().add(OFF), t0, frac, l_subfr, flag3) };
            }
            assert_eq!(
                r, c,
                "Pred_lt_3or6(T0={t0}, frac={frac}, L={l_subfr}, flag3={flag3})"
            );
        }
    }
}
