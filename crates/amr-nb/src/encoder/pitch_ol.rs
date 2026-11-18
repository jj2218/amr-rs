//! Open-loop pitch analysis: `comp_corr`, `Pitch_ol` (with its `Lag_max`), the MR102 weighted
//! variant `Pitch_ol_wgh`, `hp_max` (VAD complex-signal measure) and the `ol_ltp` dispatcher.
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `enc/src/{calc_cor,pitch_ol,p_ol_wgh,
//! hp_max,ol_ltp}.cpp`.
//!
//! Signals are passed as slices that start `PIT_MAX` samples before the C's `signal[0]`.

use super::lpc::l_abs;
use super::vad::{
    Vad, vad_complex_detection_update, vad_tone_detection, vad_tone_detection_update,
};
use crate::basic_op::{
    Flag, MAX_16, MAX_32, MIN_16, MIN_32, div_s, l_add, l_mac, l_mac_dot, l_mac_energy, l_msu,
    l_shl, l_shr, l_sub, mpy_32, negate, norm_l, round, shl, shr, sub,
};
use crate::cnst::{CodecMode, L_FRAME, PIT_MAX};
use crate::gains::gmed_n;
use crate::math::inv_sqrt;
use crate::tables::CORRWEIGHT;

const THRESHOLD: i32 = 27853;
const PIT_MIN: usize = 20;
const PIT_MIN_MR122: usize = 18;
const L_FRAME_BY2: usize = 80;
/// Buffer of `PIT_MAX` history samples plus up to one frame.
const SIG_LEN: usize = PIT_MAX + L_FRAME;

/// Correlations indexed by lag (the C's `corr[-lag]`).
pub(crate) type Corr = [i32; PIT_MAX + 1];

#[inline(always)]
pub(super) fn split(x: i32) -> (i16, i16) {
    let hi = (x >> 16) as i16;
    (hi, ((x >> 1) - ((hi as i32) << 15)) as i16)
}

/// `comp_corr`: `corr[lag] = 2 * sum(s[n] * s[n - lag])` for `lag_min..=lag_max` (in groups of
/// four lags, so a few lags below `lag_min` may be written too). `s[PIT_MAX]` is `scal_sig[0]`.
/// SSE2 on both paths (`simd::dots4`; no AVX2 copy: the auto-vectorised one lost to it).
pub(crate) fn comp_corr(
    s: &[i16],
    l_frame: usize,
    lag_max: usize,
    lag_min: usize,
    corr: &mut Corr,
) {
    match l_frame {
        80 => comp_corr_lags::<80>(s, lag_max, lag_min, corr),
        _ => comp_corr_lags::<160>(s, lag_max, lag_min, corr),
    }
}

/// `comp_corr`'s lags in groups of four from lag_max down, each group one `simd::dots4`.
fn comp_corr_lags<const L: usize>(s: &[i16], lag_max: usize, lag_min: usize, corr: &mut Corr) {
    let groups = ((lag_max - lag_min) >> 2) + 1;
    let x: &[i16; L] = s[PIT_MAX..PIT_MAX + L].try_into().unwrap();
    for g in 0..groups {
        // lags lag_max - 4g - k at s[PIT_MAX - lag_max + 4g + k..]
        let d = crate::simd::dots4::<L>(x, s, PIT_MAX - lag_max + 4 * g);
        for (k, v) in d.into_iter().enumerate() {
            corr[lag_max - 4 * g - k] = v.wrapping_shl(1);
        }
    }
}

/// `comp_corr`'s lags as plain arrays: `comp_corr_lags`' oracle.
#[cfg(test)]
fn comp_corr_n_portable<const L: usize>(
    s: &[i16],
    lag_max: usize,
    lag_min: usize,
    corr: &mut Corr,
) {
    let x: &[i16; L] = s[PIT_MAX..PIT_MAX + L].try_into().unwrap();
    let groups = ((lag_max - lag_min) >> 2) + 1;
    for lag in (lag_max + 1 - 4 * groups..=lag_max).rev() {
        let y: &[i16; L] = s[PIT_MAX - lag..PIT_MAX - lag + L].try_into().unwrap();
        let t = x
            .iter()
            .zip(y)
            .fold(0i32, |t, (&a, &b)| t.wrapping_add(a as i32 * b as i32));
        corr[lag] = t.wrapping_shl(1);
    }
}

/// `comp_corr` over the scaled signal `s[..n]`, saturating like TS 26.073: every partial sum is
/// bounded by its energy, so with 2 sum s^2 <= MAX_32 the (vectorised) plain sums are exact;
/// a loud signal can stay above that after the `>> 3` scaling (opencore wraps there): then the
/// C's `L_mac` chain per lag.
fn comp_corr_exact(
    s: &[i16],
    n: usize,
    l_frame: usize,
    lag_max: usize,
    lag_min: usize,
    corr: &mut Corr,
) {
    let e: i64 = s[..n].iter().map(|&v| 2 * v as i64 * v as i64).sum();
    if e <= MAX_32 as i64 {
        return comp_corr(s, l_frame, lag_max, lag_min, corr);
    }
    let ov = &mut false;
    for lag in lag_min..=lag_max {
        corr[lag] = (0..l_frame).fold(0, |t, j| l_mac(t, s[PIT_MAX + j], s[PIT_MAX + j - lag], ov));
    }
}

/// `Lag_max` of `pitch_ol.cpp`: best lag in `lag_min..=lag_max` and its normalised correlation.
/// Also returns the lag's correlation and energy (`rmax`, `r0`, used by VAD2).
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn lag_max(
    vad: &mut Vad,
    corr: &Corr,
    s: &[i16],
    scal_fac: i16,
    scal_flag: bool,
    l_frame: usize,
    lag_max: usize,
    lag_min: usize,
    cor_max: &mut i16,
    dtx: bool,
    ov: &mut Flag,
) -> (usize, i32, i32) {
    let mut max = MIN_32;
    let mut p_max = lag_max;
    for i in (lag_min..=lag_max).rev() {
        if corr[i] >= max {
            max = corr[i];
            p_max = i;
        }
    }
    // the C's L_mac chain (saturating; opencore wraps)
    let mut t0 = l_mac_energy(&s[PIT_MAX - p_max..PIT_MAX - p_max + l_frame], ov);
    let r0 = t0;
    if dtx && let Vad::One(v) = vad {
        vad_tone_detection(v, max, t0, ov);
    }
    t0 = inv_sqrt(t0, ov);
    if scal_flag {
        t0 = if t0 > 0x3fff_ffff { MAX_32 } else { t0 << 1 };
    }
    let (max_h, max_l) = split(max);
    let (ener_h, ener_l) = split(t0);
    t0 = mpy_32(max_h, max_l, ener_h, ener_l, ov);
    *cor_max = if scal_flag {
        t0 = l_shr(t0, scal_fac, ov);
        if t0 > 0xffff {
            MAX_16
        } else if t0 < -0x10000 {
            MIN_16
        } else {
            (t0 >> 1) as i16
        }
    } else {
        t0 as i16
    };
    (p_max, max, r0)
}

/// `Pitch_ol`: open-loop pitch lag over `l_frame` samples; `sig[PIT_MAX]` is `signal[0]`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(crate) fn pitch_ol(
    vad: &mut Vad,
    mode: CodecMode,
    sig: &[i16],
    pit_min: usize,
    l_frame: usize,
    idx: i16,
    dtx: bool,
    ov: &mut Flag,
) -> i16 {
    let pit_max = PIT_MAX;
    if dtx && let Vad::One(v) = vad {
        vad_tone_detection_update(v, matches!(mode, CodecMode::MR475 | CodecMode::MR515), ov);
    }
    let n = pit_max + l_frame;
    let mut t0: i32 = 0;
    for &x in &sig[..n] {
        t0 = t0.wrapping_add((x as i32 * x as i32) << 1);
        if t0 < 0 {
            t0 = MAX_32;
            break;
        }
    }
    let mut scaled = [0i16; SIG_LEN];
    let scal_fac: i16 = if t0 == MAX_32 {
        for i in 0..n {
            scaled[i] = sig[i] >> 3;
        }
        3
    } else if t0 < 1_048_576 {
        for i in 0..n {
            scaled[i] = ((sig[i] as i32) << 3) as i16;
        }
        -3
    } else {
        scaled[..n].copy_from_slice(&sig[..n]);
        0
    };
    let mut corr: Corr = [0; PIT_MAX + 1];
    comp_corr_exact(&scaled, n, l_frame, pit_max, pit_min, &mut corr);
    let scal_flag = mode == CodecMode::MR122;

    let (mut max1, mut max2, mut max3) = (0i16, 0i16, 0i16);
    let j = 4 * pit_min;
    let (mut p_max1, mut rmax1, mut r01) = lag_max(
        vad, &corr, &scaled, scal_fac, scal_flag, l_frame, pit_max, j, &mut max1, dtx, ov,
    );
    let i = j - 1;
    let j = pit_min << 1;
    let (p_max2, rmax2, r02) = lag_max(
        vad, &corr, &scaled, scal_fac, scal_flag, l_frame, i, j, &mut max2, dtx, ov,
    );
    let i = j - 1;
    let (p_max3, rmax3, r03) = lag_max(
        vad, &corr, &scaled, scal_fac, scal_flag, l_frame, i, pit_min, &mut max3, dtx, ov,
    );
    if dtx
        && idx == 1
        && let Vad::One(v) = vad
    {
        let mut corr_hp_max = 0;
        hp_max(
            &corr,
            &scaled,
            l_frame,
            pit_max,
            pit_min,
            &mut corr_hp_max,
            ov,
        );
        vad_complex_detection_update(v, corr_hp_max);
    }
    if ((((max1 as i32) * THRESHOLD) >> 15) as i16) < max2 {
        max1 = max2;
        p_max1 = p_max2;
        (rmax1, r01) = (rmax2, r02);
    }
    if ((((max1 as i32) * THRESHOLD) >> 15) as i16) < max3 {
        p_max1 = p_max3;
        (rmax1, r01) = (rmax3, r03);
    }
    if dtx && let Vad::Two(v) = vad {
        v.l_rmax = l_add(v.l_rmax, rmax1, ov);
        v.l_r0 = l_add(v.l_r0, r01, ov);
    }
    p_max1 as i16
}

/// `hp_max`: maximum of the high-passed correlation, normalised by the high-passed energy.
pub(super) fn hp_max(
    corr: &Corr,
    s: &[i16],
    l_frame: usize,
    lag_max: usize,
    lag_min: usize,
    cor_hp_max: &mut i16,
    ov: &mut Flag,
) {
    let mut max = MIN_32;
    for i in (lag_min + 1..lag_max).rev() {
        let t0 = l_shl(corr[i], 1, ov);
        let l_temp = l_sub(t0, corr[i + 1], ov);
        let t0 = l_abs(l_sub(l_temp, corr[i - 1], ov));
        if t0 >= max {
            max = t0;
        }
    }
    let mut t0: i32 = 0;
    for n in 0..l_frame {
        t0 = l_mac(t0, s[PIT_MAX + n], s[PIT_MAX + n], ov);
    }
    let mut t1: i32 = 0;
    for n in 0..l_frame {
        t1 = l_mac(t1, s[PIT_MAX + n], s[PIT_MAX + n - 1], ov);
    }
    let l_temp = l_shl(t0, 1, ov);
    let t1 = l_shl(t1, 1, ov);
    let t0 = l_abs(l_sub(l_temp, t1, ov));
    let shift1 = norm_l(max) - 1;
    let max16 = (l_shl(max, shift1, ov) >> 16) as i16;
    let shift2 = norm_l(t0);
    let t016 = (l_shl(t0, shift2, ov) >> 16) as i16;
    let cor_max = if t016 != 0 { div_s(max16, t016) } else { 0 };
    let shift = shift1 - shift2;
    *cor_hp_max = if shift >= 0 {
        shr(cor_max, shift, ov)
    } else {
        shl(cor_max, negate(shift), ov)
    };
}

/// `pitchOLWghtState`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PitchOlWghtState {
    pub old_t0_med: i16,
    pub ada_w: i16,
    pub wght_flg: i16,
}

impl Default for PitchOlWghtState {
    fn default() -> Self {
        Self {
            old_t0_med: 40,
            ada_w: 0,
            wght_flg: 0,
        }
    }
}

/// `Lag_max` of `p_ol_wgh.cpp`: weighted lag search (lag weighting, optional neighbourhood
/// weighting around `old_lag`) plus the open-loop gain flag.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn lag_max_wgh(
    vad: &mut Vad,
    corr: &Corr,
    s: &[i16],
    l_frame: usize,
    lag_max: usize,
    lag_min: usize,
    old_lag: i16,
    cor_max: &mut i16,
    wght_flg: i16,
    gain_flg: &mut i16,
    dtx: bool,
    ov: &mut Flag,
) -> usize {
    // The C weights lag i (from lag_max down) by CORRWEIGHT[250 - (lag_max - i)], then (if
    // wght_flg) by CORRWEIGHT[123 + i - old_lag], each an `L_Extract` + `Mpy_32_16`, and keeps
    // the first maximum (>=, so the lowest lag of a tie). With weights >= 0 and lo in
    // 0..=0x7fff, `Mpy_32_16` can't saturate (<= 2 * 32767^2 + 2 * 32767): it is
    // 2 hi w + 2 ((lo w) >> 15) in i32, so every lag is one lane, then a max and a first
    // position (as `lag_max`).
    let wmul = |x: i32, w: i16| {
        let (hi, lo) = (x >> 16, (x >> 1) & 0x7fff);
        2 * hi * w as i32 + 2 * ((lo * w as i32) >> 15)
    };
    let mut t = [0i32; PIT_MAX + 1];
    let w1 = &CORRWEIGHT[250 + lag_min - lag_max..=250];
    for ((t, &c), &w) in t[lag_min..=lag_max]
        .iter_mut()
        .zip(&corr[lag_min..=lag_max])
        .zip(w1)
    {
        *t = wmul(c, w);
    }
    if wght_flg > 0 {
        let w2 = &CORRWEIGHT[(123 + lag_min as isize - old_lag as isize) as usize..];
        for (t, &w) in t[lag_min..=lag_max].iter_mut().zip(w2) {
            *t = wmul(*t, w);
        }
    }
    let tl = &t[lag_min..=lag_max];
    let max = tl.iter().copied().fold(MIN_32, i32::max);
    let p_max = lag_min + tl.iter().position(|&v| v == max).unwrap_or(0);
    // the C's two L_mac chains, as the exact helpers
    let (x, xl) = (
        &s[PIT_MAX..PIT_MAX + l_frame],
        &s[PIT_MAX - p_max..PIT_MAX - p_max + l_frame],
    );
    let t0 = l_mac_dot(0, x, xl, ov);
    let t1 = l_mac_energy(xl, ov);
    if dtx {
        match vad {
            Vad::One(v) => {
                vad_tone_detection_update(v, false, ov);
                vad_tone_detection(v, t0, t1, ov);
            }
            Vad::Two(v) => {
                v.l_rmax = l_add(v.l_rmax, t0, ov);
                v.l_r0 = l_add(v.l_r0, t1, ov);
            }
        }
    }
    let temp = round(t1, ov);
    let t1 = l_msu(t0, temp, 13107, ov);
    *gain_flg = round(t1, ov);
    *cor_max = 0;
    p_max
}

/// `Pitch_ol_wgh`: MR102 open-loop pitch with lag weighting; `sig[PIT_MAX]` is `signal[0]`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(crate) fn pitch_ol_wgh(
    st: &mut PitchOlWghtState,
    vad: &mut Vad,
    sig: &[i16],
    pit_min: usize,
    l_frame: usize,
    old_lags: &mut [i16; 5],
    ol_gain_flg: &mut [i16; 2],
    idx: usize,
    dtx: bool,
    ov: &mut Flag,
) -> i16 {
    let pit_max = PIT_MAX;
    let n = pit_max + l_frame;
    let mut t0: i32 = 0;
    for &x in &sig[..n] {
        t0 = l_mac(t0, x, x, ov);
    }
    let mut scaled = [0i16; SIG_LEN];
    if l_sub(t0, MAX_32, ov) == 0 {
        for i in 0..n {
            scaled[i] = shr(sig[i], 3, ov);
        }
    } else if l_sub(t0, 1_048_576, ov) < 0 {
        for i in 0..n {
            scaled[i] = shl(sig[i], 3, ov);
        }
    } else {
        scaled[..n].copy_from_slice(&sig[..n]);
    }
    let mut corr: Corr = [0; PIT_MAX + 1];
    comp_corr_exact(&scaled, n, l_frame, pit_max, pit_min, &mut corr);
    let mut max1 = 0;
    let p_max1 = lag_max_wgh(
        vad,
        &corr,
        &scaled,
        l_frame,
        pit_max,
        pit_min,
        st.old_t0_med,
        &mut max1,
        st.wght_flg,
        &mut ol_gain_flg[idx],
        dtx,
        ov,
    ) as i16;
    if ol_gain_flg[idx] > 0 {
        old_lags.copy_within(0..4, 1);
        old_lags[0] = p_max1;
        st.old_t0_med = gmed_n(old_lags, 5);
        st.ada_w = 32767;
    } else {
        st.old_t0_med = p_max1;
        st.ada_w = ((st.ada_w as i32 * 29491) >> 15) as i16;
    }
    st.wght_flg = (sub(st.ada_w, 9830, ov) >= 0) as i16;
    if dtx
        && idx == 1
        && let Vad::One(v) = vad
    {
        let mut corr_hp_max = 0;
        hp_max(
            &corr,
            &scaled,
            l_frame,
            pit_max,
            pit_min,
            &mut corr_hp_max,
            ov,
        );
        vad_complex_detection_update(v, corr_hp_max);
    }
    p_max1
}

/// `ol_ltp`: open-loop pitch lag for half a frame (`idx` = 0/1), or the whole frame for
/// MR475/MR515. `wsp[PIT_MAX]` is the C's `wsp[0]` at the analysed position.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(crate) fn ol_ltp(
    st: &mut PitchOlWghtState,
    vad: &mut Vad,
    mode: CodecMode,
    wsp: &[i16],
    old_lags: &mut [i16; 5],
    ol_gain_flg: &mut [i16; 2],
    idx: usize,
    dtx: bool,
    ov: &mut Flag,
) -> i16 {
    if mode != CodecMode::MR102 {
        *ol_gain_flg = [0, 0];
    }
    match mode {
        CodecMode::MR475 | CodecMode::MR515 => {
            pitch_ol(vad, mode, wsp, PIT_MIN, L_FRAME, idx as i16, dtx, ov)
        }
        m if m <= CodecMode::MR795 => {
            pitch_ol(vad, mode, wsp, PIT_MIN, L_FRAME_BY2, idx as i16, dtx, ov)
        }
        CodecMode::MR102 => pitch_ol_wgh(
            st,
            vad,
            wsp,
            PIT_MIN,
            L_FRAME_BY2,
            old_lags,
            ol_gain_flg,
            idx,
            dtx,
            ov,
        ),
        _ => pitch_ol(
            vad,
            mode,
            wsp,
            PIT_MIN_MR122,
            L_FRAME_BY2,
            idx as i16,
            dtx,
            ov,
        ),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{N, Rng};

    /// `simd::dots4` is `dots4_portable` (on x86_64): any x and window, all -32768 too
    /// (`pmaddwd`'s pair sum 2^31 wraps like the wrapping sums), both lengths in use.
    #[test]
    fn dots4_is_the_array_form() {
        let mut g = Rng::new(29);
        for case in 0..N {
            let full = case % 8 == 0;
            let w: [i16; 256] = core::array::from_fn(|_| if full { -32768 } else { g.i16() });
            let x: [i16; 240] = core::array::from_fn(|_| if full { -32768 } else { g.i16() });
            let off = (g.u64() % 13) as usize;
            let x80: &[i16; 80] = x[..80].try_into().unwrap();
            assert_eq!(
                crate::simd::dots4(x80, &w, off),
                crate::simd::dots4_portable(x80, &w, off)
            );
            assert_eq!(
                crate::simd::dots4(&x, &w, off),
                crate::simd::dots4_portable(&x, &w, off)
            );
        }
    }

    /// `comp_corr_lags` (on `simd::dots4`) is `comp_corr_n_portable`: full-scale signals (the
    /// i32 sums wrap), both frame lengths, any lag range the codec uses.
    #[test]
    fn comp_corr_lags_is_the_array_form() {
        let mut g = Rng::new(23);
        for case in 0..N / 8 {
            let sh = (g.u64() % 16) as u32;
            let s: [i16; SIG_LEN] = if case % 16 == 0 {
                [-32768; SIG_LEN] // pmaddwd's (-32768)^2 pairs wrap to -2^31
            } else {
                core::array::from_fn(|_| g.i16() >> sh)
            };
            let lag_min = 18 + (g.u64() % 42) as usize; // PIT_MIN_MR122 = 18
            let lag_max = lag_min + (g.u64() % (PIT_MAX - lag_min + 1) as u64) as usize;
            let (mut a, mut b) = ([7i32; PIT_MAX + 1], [7i32; PIT_MAX + 1]);
            if case % 2 == 0 {
                comp_corr_lags::<80>(&s, lag_max, lag_min, &mut a);
                comp_corr_n_portable::<80>(&s, lag_max, lag_min, &mut b);
            } else {
                comp_corr_lags::<160>(&s, lag_max, lag_min, &mut a);
                comp_corr_n_portable::<160>(&s, lag_max, lag_min, &mut b);
            }
            assert_eq!(a, b, "lags {lag_min}..={lag_max}");
        }
    }
}
