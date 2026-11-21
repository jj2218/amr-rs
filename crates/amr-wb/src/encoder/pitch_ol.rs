//! Open-loop pitch (vo-amrwbenc `p_med_ol.c`, `hp_wsp.c`), with the reference's saturation.
//! The lag search runs on the weighted speech decimated by 2 (`OPL_DECIM`).

use super::ops::*;
use crate::enc_tables::{CORRWEIGHT, HP_WSP_A as A, HP_WSP_B as B};
use crate::excitation::median5;
pub(crate) use crate::simd::pitch_med_ol;

/// `Init_Hp_wsp` is a zeroed `[i16; 9]`. `scale_mem_Hp_wsp`: rescales the high-pass memory
/// by `2^exp` (the outputs in double precision, the inputs rounded).
pub(crate) fn scale_mem_hp_wsp(mem: &mut [i16; 9], exp: i16) {
    for i in (0..6).step_by(2) {
        (mem[i], mem[i + 1]) = l_extract(l_shl(l_comp(mem[i], mem[i + 1]), exp));
    }
    for v in &mut mem[6..] {
        *v = round(l_shl(l_deposit_h(*v), exp));
    }
}

/// `Hp_wsp`: 3rd-order high-pass of `wsp` into `hp_wsp`.
pub(crate) fn hp_wsp(wsp: &[i16], hp_wsp: &mut [i16], mem: &mut [i16; 9]) {
    let [
        mut y3_hi,
        mut y3_lo,
        mut y2_hi,
        mut y2_lo,
        mut y1_hi,
        mut y1_lo,
        mut x0,
        mut x1,
        mut x2,
    ] = *mem;
    for (out, &x) in hp_wsp.iter_mut().zip(wsp) {
        let x3 = x2;
        x2 = x1;
        x1 = x0;
        x0 = x;
        let mut l = 16384i32;
        l = l_mac(l, y1_lo, A[1]);
        l = l_mac(l, y2_lo, A[2]);
        l = l_mac(l, y3_lo, A[3]);
        l = l_shr(l, 15);
        l = l_mac(l, y1_hi, A[1]);
        l = l_mac(l, y2_hi, A[2]);
        l = l_mac(l, y3_hi, A[3]);
        l = l_mac(l, x0, B[0]);
        l = l_mac(l, x1, B[1]);
        l = l_mac(l, x2, B[2]);
        l = l_mac(l, x3, B[3]);
        l = l_shl(l, 2);
        (y3_hi, y3_lo) = (y2_hi, y2_lo);
        (y2_hi, y2_lo) = (y1_hi, y1_lo);
        (y1_hi, y1_lo) = l_extract(l);
        *out = round(l_shl(l, 1));
    }
    *mem = [y3_hi, y3_lo, y2_hi, y2_lo, y1_hi, y1_lo, x0, x1, x2];
}

/// `Pitch_med_ol`: the open-loop lag in (`l_min`, `l_max`] maximising the weighted
/// correlation of `wsp[at..at + l_frame]` (history before `at`), and its normalised
/// correlation on the high-passed signal (`gain`). `old_hp_wsp` holds `l_max` past high-passed
/// samples, then room for `l_frame`.
#[inline(always)]
pub(crate) fn pitch_med_ol_portable(
    wsp: &[i16],
    at: usize,
    l_min: i16,
    l_max: i16,
    l_frame: usize,
    l_0: i16,
    gain: &mut i16,
    hp_wsp_mem: &mut [i16; 9],
    old_hp_wsp: &mut [i16],
    wght_flg: i16,
) -> i16 {
    let mut ww = 198usize;
    let mut we = (98 + l_max - l_0) as usize;
    let mut max = i32::MIN;
    let mut tm = 0i16;
    let x = &wsp[at..at + l_frame];
    // every partial sum is at most sqrt(Σx²·Σw²) in magnitude over the window w of all lags
    // (Cauchy–Schwarz): below 2^30 the chains are plain sums
    let e = |v: &[i16]| v.iter().map(|&s| s as i64 * s as i64).sum::<i64>();
    let w = &wsp[at - l_max as usize..at + l_frame - 1 - l_min as usize];
    const B: i128 = (1 << 30) - 1;
    let plain = e(x) as i128 * e(w) as i128 <= B * B;
    for i in (l_min + 1..=l_max).rev() {
        let y = &wsp[at - i as usize..at - i as usize + l_frame];
        let mut r0 = if plain { dot(x, y) << 1 } else { 0 };
        if !plain {
            for j in 0..l_frame {
                r0 = l_mac(r0, x[j], y[j]);
            }
        }
        let (hi, lo) = l_extract(r0);
        r0 = mpy_32_16(hi, lo, CORRWEIGHT[ww]);
        ww -= 1;
        if l_0 > 0 && wght_flg > 0 {
            let (hi, lo) = l_extract(r0);
            r0 = mpy_32_16(hi, lo, CORRWEIGHT[we]);
            we = we.wrapping_sub(1);
        }
        if r0 >= max {
            max = r0;
            tm = i;
        }
    }
    let lm = l_max as usize;
    {
        let (_, hp) = old_hp_wsp.split_at_mut(lm);
        hp_wsp(x, &mut hp[..l_frame], hp_wsp_mem);
    }
    let hp = &old_hp_wsp[..lm + l_frame];
    let tm_u = tm as usize;
    let (mut r0, mut r1, mut r2) = (0i32, 1i32, 1i32);
    for j in lm..lm + l_frame {
        r0 = l_mac(r0, hp[j], hp[j - tm_u]);
        r1 = l_mac(r1, hp[j - tm_u], hp[j - tm_u]);
        r2 = l_mac(r2, hp[j], hp[j]);
    }
    let exp_r0 = norm_l(r0);
    let r0 = l_shl(r0, exp_r0);
    let mut exp_r1 = norm_l(r1);
    let r1 = l_shl(r1, exp_r1);
    let exp_r2 = norm_l(r2);
    let r2 = l_shl(r2, exp_r2);
    let mut r1 = l_mult(round(r1), round(r2));
    let i = norm_l(r1);
    r1 = l_shl(r1, i);
    exp_r1 = sub(62, add(add(exp_r1, exp_r2), i));
    isqrt_n(&mut r1, &mut exp_r1);
    let r0 = l_mult(round(r0), round(r1));
    let exp_r0 = add(sub(31, exp_r0), exp_r1);
    *gain = round(l_shl(r0, exp_r0));
    old_hp_wsp.copy_within(l_frame..l_frame + lm, 0);
    tm
}

/// `Σ x·y` over blocks of 16 (wrapping: exact when the caller's bound holds).
#[inline(always)]
fn dot(x: &[i16], y: &[i16]) -> i32 {
    let (xb, xr) = x.as_chunks::<16>();
    let (yb, yr) = y.as_chunks::<16>();
    let mut acc = [0i32; 16];
    for (a, b) in xb.iter().zip(yb) {
        for k in 0..16 {
            acc[k] = acc[k].wrapping_add(a[k] as i32 * b[k] as i32);
        }
    }
    let s = xr
        .iter()
        .zip(yr)
        .fold(0i32, |s, (&a, &b)| s.wrapping_add(a as i32 * b as i32));
    acc.iter().fold(s, |s, &a| s.wrapping_add(a))
}

/// `Med_olag`: pushes `prev_ol_lag` into the 5 past lags and returns their median.
pub(crate) fn med_olag(prev_ol_lag: i16, old_ol_lag: &mut [i16; 5]) -> i16 {
    old_ol_lag.copy_within(0..4, 1);
    old_ol_lag[0] = prev_ol_lag;
    median5(old_ol_lag, 2)
}

#[cfg(all(test, oracle))]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::encoder::filters::vs_ref::signal;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_Pitch_med_ol(
            wsp: *mut i16,
            l_min: i16,
            l_max: i16,
            l_frame: i16,
            l_0: i16,
            gain: *mut i16,
            hp_wsp_mem: *mut i16,
            old_hp_wsp: *mut i16,
            wght_flg: i16,
        ) -> i16;
        fn refwb_Med_olag(prev: i16, old: *mut i16) -> i16;
        fn refwb_scale_mem_Hp_wsp(mem: *mut i16, exp: i16);
    }

    /// Sequences of half and whole frames (64, 128) over signal classes (the 12-bit scaled
    /// weighted speech and louder), random old lags and weighting flags, memory rescaling.
    #[test]
    fn pitch_med_ol_matches_reference() {
        let mut r = Rng::new(41);
        for seq in 0..60u64 {
            let (mut mem, mut memc) = ([0i16; 9], [0i16; 9]);
            let (mut old, mut oldc) = ([0i16; 243], [0i16; 243]);
            let (mut lags, mut lagsc) = ([40i16; 5], [40i16; 5]);
            for call in 0..100 {
                let l_frame = if r.u64().is_multiple_of(2) { 64 } else { 128 };
                let mut w = signal(&mut r, seq + call / 20, 115 + l_frame);
                if seq % 2 == 0 {
                    for v in &mut w {
                        *v >>= 3;
                    }
                }
                // old_T0_med: a median of past lags in (17, 115]
                let l_0 = 18 + (r.u64() % 98) as i16;
                let flg = (r.u64() % 2) as i16;
                let (mut g, mut gc) = (0i16, 0i16);
                let t = pitch_med_ol(
                    &w, 115, 17, 115, l_frame, l_0, &mut g, &mut mem, &mut old, flg,
                );
                let exp = (r.u64() % 7) as i16 - 3;
                scale_mem_hp_wsp(&mut mem, exp);
                let m = med_olag(t, &mut lags);
                with_ref(|| {
                    // SAFETY: wsp has 115 samples of history; the memories their sizes.
                    unsafe {
                        let tc = refwb_Pitch_med_ol(
                            w.as_mut_ptr().add(115),
                            17,
                            115,
                            l_frame as i16,
                            l_0,
                            &mut gc,
                            memc.as_mut_ptr(),
                            oldc.as_mut_ptr(),
                            flg,
                        );
                        assert_eq!((t, g), (tc, gc), "Pitch_med_ol {seq} {call}");
                        refwb_scale_mem_Hp_wsp(memc.as_mut_ptr(), exp);
                        assert_eq!((mem, old), (memc, oldc), "Hp_wsp memories {seq} {call}");
                        assert_eq!(m, refwb_Med_olag(tc, lagsc.as_mut_ptr()), "Med_olag");
                    }
                });
            }
        }
    }
}
