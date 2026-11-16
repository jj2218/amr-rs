//! Closed-loop pitch (vo-amrwbenc `pitch_f4.c`, `g_pitch.c`) and pitch-gain clipping
//! (`gpclip.c`, with CR 26.173-0032's 6.60/8.85 thresholds, which vo predates), with the
//! reference's saturation.

use super::filters::convolve;
use super::ops::*;
use crate::basic_op::dot_product12;
use crate::cnst::{L_SUBFR, M, NBBITS_7K, NBBITS_9K, PIT_MIN};
use crate::enc_tables::INTER4_1;

const UP_SAMP: i16 = 4;
const L_INTERPOL1: usize = 4;

/// `Pitch_fr4`: the closed-loop lag in [`t0_min`, `t0_max`] (integer search on the normalised
/// correlation, then fractions: 1/4 below `t0_fr2`, 1/2 below `t0_fr1`); `exc` is the
/// excitation with its past before `at`. Returns the lag; `pit_frac` its fraction.
pub(crate) fn pitch_fr4(
    exc: &[i16],
    at: usize,
    xn: &[i16],
    h: &[i16],
    t0_min: i16,
    t0_max: i16,
    pit_frac: &mut i16,
    i_subfr: usize,
    t0_fr2: i16,
    t0_fr1: i16,
    l_subfr: usize,
) -> i16 {
    let t_min = t0_min - L_INTERPOL1 as i16;
    let t_max = t0_max + L_INTERPOL1 as i16;
    // corr_v[t - t_min] for t in [t_min, t_max]
    let mut corr_v = [0i16; 40];
    norm_corr(exc, at, xn, h, l_subfr, t_min, t_max, &mut corr_v);
    let c = |t: i16| (t - t_min) as usize;
    let mut max = corr_v[c(t0_min)];
    let mut t0 = t0_min;
    for i in t0_min + 1..=t0_max {
        if corr_v[c(i)] >= max {
            max = corr_v[c(i)];
            t0 = i;
        }
    }
    if i_subfr == 0 && t0 >= t0_fr1 {
        *pit_frac = 0;
        return t0;
    }
    let (mut step, mut fraction) = (1, -3);
    if (i_subfr == 0 && t0 >= t0_fr2) || t0_fr2 == PIT_MIN {
        step = 2;
        fraction = -2;
    }
    if t0 == t0_min {
        fraction = 0;
    }
    let mut max = interpol_4(&corr_v, c(t0), fraction);
    let mut i = fraction + step;
    while i <= 3 {
        let temp = interpol_4(&corr_v, c(t0), i);
        if temp > max {
            max = temp;
            fraction = i;
        }
        i += step;
    }
    if fraction < 0 {
        fraction += UP_SAMP;
        t0 -= 1;
    }
    *pit_frac = fraction;
    t0
}

/// `Norm_Corr`: the normalised correlation of `xn` with the filtered past excitation at each
/// lag in [`t_min`, `t_max`] (into `corr_norm[t - t_min]`).
fn norm_corr(
    exc: &[i16],
    at: usize,
    xn: &[i16],
    h: &[i16],
    l_subfr: usize,
    t_min: i16,
    t_max: i16,
    corr_norm: &mut [i16],
) {
    let mut k = at - t_min as usize;
    let mut excf = [0i16; L_SUBFR];
    convolve(&exc[k..], h, &mut excf, l_subfr);
    let mut l = 1i32;
    for &x in &xn[..l_subfr] {
        l = l_mac(l, x, x);
    }
    let exp = add(sub(30, norm_l(l)), 2);
    let scale = negate(shr(exp, 1));
    for t in t_min..=t_max {
        let mut l = 1i32;
        for i in 0..l_subfr {
            l = l_mac(l, xn[i], excf[i]);
        }
        let exp = norm_l(l);
        let corr = extract_h(l_shl(l, exp));
        let exp_corr = sub(30, exp);
        let mut l = 1i32;
        for &e in &excf[..l_subfr] {
            l = l_mac(l, e, e);
        }
        let exp = norm_l(l);
        let mut l = l_shl(l, exp);
        let mut exp_norm = sub(30, exp);
        isqrt_n(&mut l, &mut exp_norm);
        let norm = extract_h(l);
        let l = l_mult(corr, norm);
        corr_norm[(t - t_min) as usize] = round(l_shl(l, add(add(exp_corr, exp_norm), scale)));
        if t != t_max {
            k -= 1;
            for i in (1..l_subfr).rev() {
                excf[i] = add(mult(exc[k], h[i]), excf[i - 1]);
            }
            excf[0] = mult(exc[k], h[0]);
        }
    }
}

/// `Interpol_4`: `x` interpolated at `pos + frac/4` (frac in -4..=3).
fn interpol_4(x: &[i16], pos: usize, frac: i16) -> i16 {
    let (mut pos, mut frac) = (pos, frac);
    if frac < 0 {
        frac += UP_SAMP;
        pos -= 1;
    }
    let x = &x[pos + 1 - L_INTERPOL1..];
    // vo's table is phase-major: row 3 - frac holds that phase's 8 taps
    let k = (UP_SAMP - 1 - frac) as usize;
    let mut s = 0i32;
    for (&v, &c) in x[..2 * L_INTERPOL1].iter().zip(&INTER4_1[8 * k..8 * k + 8]) {
        s = l_mac(s, v, c);
    }
    round(l_shl(s, 1))
}

/// `G_pitch`: the pitch gain (Q14, in [0, 1.2]) of `y1` for target `xn`; `g_coeff` the
/// correlations (yy, exp_yy, xy, exp_xy) for the gain quantiser.
pub(crate) fn g_pitch(xn: &[i16], y1: &[i16], g_coeff: &mut [i16; 4]) -> i16 {
    let (mut exp_yy, mut exp_xy) = (0, 0);
    let yy = extract_h(dot_product12(y1, y1, &mut exp_yy));
    let xy = extract_h(dot_product12(xn, y1, &mut exp_xy));
    *g_coeff = [yy, exp_yy, xy, exp_xy];
    if xy < 0 {
        return 0;
    }
    let gain = div_s(xy >> 1, yy);
    let i = sub(exp_xy, exp_yy);
    shl(gain, i).min(19661)
}

const DIST_ISF_MAX_IO: i16 = 384;
const DIST_ISF_MAX: i16 = 307;
const DIST_ISF_THRES: i16 = 154;
const GAIN_PIT_THRES: i16 = 14746;
const GAIN_PIT_MIN: i16 = 9830;

/// 6.60 and 8.85 kbit/s: CR 26.173-0032's clipping thresholds.
fn low_rate(ser_size: i16) -> bool {
    ser_size == NBBITS_7K || ser_size == NBBITS_9K
}

/// `Init_gp_clip`.
pub(crate) const GP_CLIP_INIT: [i16; 2] = [DIST_ISF_MAX, GAIN_PIT_MIN];

/// `Gp_clip`: whether to clip the pitch gain (resonant LPC and high past pitch gain).
pub(crate) fn gp_clip(ser_size: i16, mem: &[i16; 2]) -> bool {
    if low_rate(ser_size) {
        let thres = add(
            14746,
            mult(1638, extract_l(l_mult(mem[0], 16384 / DIST_ISF_MAX_IO))),
        );
        mem[1] > thres
    } else {
        mem[0] < DIST_ISF_THRES && mem[1] > GAIN_PIT_THRES
    }
}

/// `Gp_clip_test_isf`: tracks the minimum ISF distance.
pub(crate) fn gp_clip_test_isf(ser_size: i16, isf: &[i16], mem: &mut [i16; 2]) {
    let mut dist_min = sub(isf[1], isf[0]);
    for i in 2..M - 1 {
        dist_min = dist_min.min(sub(isf[i], isf[i - 1]));
    }
    let dist = extract_h(l_mac(l_mult(26214, mem[0]), 6554, dist_min));
    let max = if low_rate(ser_size) {
        DIST_ISF_MAX_IO
    } else {
        DIST_ISF_MAX
    };
    mem[0] = dist.min(max);
}

/// `Gp_clip_test_gain_pit`: tracks the quantised pitch gain.
pub(crate) fn gp_clip_test_gain_pit(ser_size: i16, gain_pit: i16, mem: &mut [i16; 2]) {
    let l = if low_rate(ser_size) {
        l_mac(l_mult(32113, mem[1]), 655, gain_pit)
    } else {
        l_mac(l_mult(29491, mem[1]), 3277, gain_pit)
    };
    mem[1] = extract_h(l).max(GAIN_PIT_MIN);
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::cnst::{NBBITS_12K, NBBITS_24K, PIT_FR1_8B, PIT_FR1_9B, PIT_FR2, PIT_MAX};
    use crate::encoder::filters::vs_ref::signal;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_Pitch_fr4(
            exc: *mut i16,
            xn: *mut i16,
            h: *mut i16,
            t0_min: i16,
            t0_max: i16,
            pit_frac: *mut i16,
            i_subfr: i16,
            t0_fr2: i16,
            t0_fr1: i16,
            l_subfr: i16,
        ) -> i16;
        fn refwb_G_pitch(xn: *mut i16, y1: *mut i16, g_coeff: *mut i16, l: i16) -> i16;
        fn refwb_Gp_clip(ser_size: i16, mem: *mut i16) -> i16;
        fn refwb_Gp_clip_test_isf(ser_size: i16, isf: *mut i16, mem: *mut i16);
        fn refwb_Gp_clip_test_gain_pit(ser_size: i16, gain: i16, mem: *mut i16);
    }

    /// Lag ranges over [PIT_MIN, PIT_MAX] (the C's 16-wide windows), both resolutions and
    /// subframe kinds; excitation and targets from the signal classes (12-bit and louder),
    /// impulse responses in Q15. Then the gain on the same vectors.
    #[test]
    fn pitch_fr4_and_g_pitch_match_reference() {
        let mut r = Rng::new(71);
        for n in 0..20_000u64 {
            let at = PIT_MAX + 1 + 4;
            let mut exc = signal(&mut r, n, at + L_SUBFR + 1);
            let shift = (r.u64() % 4) as u32;
            for v in &mut exc {
                *v >>= shift;
            }
            let mut xn = signal(&mut r, n / 3, L_SUBFR);
            let mut h = signal(&mut r, n / 5, L_SUBFR);
            h[0] = h[0].saturating_abs().max(1000);
            let t0_min = PIT_MIN + (r.u64() % (PIT_MAX as u64 - 15 - PIT_MIN as u64 + 1)) as i16;
            let t0_max = t0_min + 15;
            let i_subfr = [0usize, 64, 128, 192][(r.u64() % 4) as usize];
            let (fr2, fr1) = if r.u64().is_multiple_of(2) {
                (PIT_MIN, PIT_FR1_8B)
            } else {
                (PIT_FR2, PIT_FR1_9B)
            };
            let mut frac = 0;
            let t = pitch_fr4(
                &exc, at, &xn, &h, t0_min, t0_max, &mut frac, i_subfr, fr2, fr1, L_SUBFR,
            );
            let mut y = [0i16; L_SUBFR];
            convolve(&exc[at..], &h, &mut y, L_SUBFR);
            let mut g = [0i16; 4];
            let gain = g_pitch(&xn, &y, &mut g);
            with_ref(|| {
                // SAFETY: exc holds PIT_MAX + 5 past samples (t_max + 1 back) and a subframe.
                unsafe {
                    let mut fc = 0;
                    let tc = refwb_Pitch_fr4(
                        exc.as_mut_ptr().add(at),
                        xn.as_mut_ptr(),
                        h.as_mut_ptr(),
                        t0_min,
                        t0_max,
                        &mut fc,
                        i_subfr as i16,
                        fr2,
                        fr1,
                        L_SUBFR as i16,
                    );
                    assert_eq!((t, frac), (tc, fc), "Pitch_fr4 {n}");
                    let mut gc = [0i16; 4];
                    let gainc = refwb_G_pitch(xn.as_mut_ptr(), y.as_mut_ptr(), gc.as_mut_ptr(), 64);
                    assert_eq!((gain, g), (gainc, gc), "G_pitch {n}");
                }
            });
        }
    }

    #[test]
    fn gp_clip_matches_reference() {
        let mut r = Rng::new(72);
        let sizes = [NBBITS_7K, NBBITS_9K, NBBITS_12K, NBBITS_24K, 35];
        let (mut m, mut mc) = (GP_CLIP_INIT, GP_CLIP_INIT);
        for n in 0..100_000 {
            let ser = sizes[(r.u64() % 5) as usize];
            let mut isf: [i16; M] =
                std::array::from_fn(|i| 200 * i as i16 + (r.u64() % 600) as i16);
            if n % 3 == 0 {
                isf[5] = isf[4] + (r.u64() % 100) as i16;
            }
            if (n / 50) % 2 == 0 {
                // evenly spread: the distance memory climbs to its rate's ceiling
                isf = std::array::from_fn(|i| 1000 * i as i16);
            }
            let g = (r.u64() % 20000) as i16;
            gp_clip_test_isf(ser, &isf, &mut m);
            let clip = gp_clip(ser, &m);
            gp_clip_test_gain_pit(ser, g, &mut m);
            with_ref(|| {
                // SAFETY: 16 ISFs and the 2-word memory.
                unsafe {
                    refwb_Gp_clip_test_isf(ser, isf.as_mut_ptr(), mc.as_mut_ptr());
                    assert_eq!(
                        clip,
                        refwb_Gp_clip(ser, mc.as_mut_ptr()) != 0,
                        "Gp_clip {n}"
                    );
                    refwb_Gp_clip_test_gain_pit(ser, g, mc.as_mut_ptr());
                }
            });
            assert_eq!(m, mc, "gp_clip memory {n}");
        }
    }
}
