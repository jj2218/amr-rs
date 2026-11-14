//! Closed-loop pitch search and its helpers: `Convolve`, `Interpol_3or6`, `Pitch_fr` (with
//! `Norm_Corr`, `searchFrac`, `getRange`), `Enc_lag3`, `Enc_lag6`, `G_pitch`, `q_gain_pitch`,
//! the tonal-stability detector (`ton_stab`) and the `cl_ltp` driver.
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `enc/src/{convolve,inter_36,pitch_fr,
//! enc_lag3,enc_lag6,g_pitch,q_gain_p,ton_stab,cl_ltp}.cpp`.
#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1
#![allow(clippy::too_many_arguments)] // signatures mirror the C functions

use crate::basic_op::{
    Flag, MAX_16, MAX_32, abs_s, add, div_s, l_add, l_mac, l_mult, l_shl, max_abs, mpy_32, norm_l,
    round, shr, sub, sum_abs,
};
use crate::cnst::{L_SUBFR, M, Mode, PIT_MAX};
use crate::excitation::pred_lt_3or6;
use crate::math::inv_sqrt;
use crate::tables::{INTER_6, QUA_GAIN_PITCH};

const L_INTER_SRCH: i16 = 4;
const UP_SAMP_MAX: i16 = 6;
const L_FRAME_BY2: i16 = 80;
const N_FRAME: usize = 7;
/// Pitch gain clipping (0.95).
pub(crate) const GP_CLIP: i16 = 15565;
const NB_QUA_PITCH: usize = 16;

pub(crate) use crate::simd::convolve;

/// `Convolve` with `L = L_SUBFR`: `y[n] = sat(sum_{i<=n} x[i] h[n-i] >> 12)`, saturating like
/// TS 26.073 (`L_mac` chain, `L_shl(s, 3)`; opencore wraps). Portable body of the dispatched
/// `convolve` (`crate::simd`).
#[inline(always)]
pub(crate) fn convolve_portable(x: &[i16], h: &[i16], y: &mut [i16]) {
    const L: usize = L_SUBFR;
    // Partial sums stay below 2^30 (no L_mac saturation) if max|x| * sum|h| does; then the
    // result is exactly clamp(sum >> 12). Otherwise the C loop.
    if max_abs(&x[..L]) * sum_abs(&h[..L]) >= 1 << 30 {
        let ov = &mut false;
        for n in 0..L {
            let s = (0..=n).fold(0, |s, i| l_mac(s, x[i], h[n - i], ov));
            y[n] = (l_shl(s, 3, ov) >> 16) as i16;
        }
        return;
    }
    // Exact i32 sums may be reordered: outputs in blocks of 8 lanes, each adding x[i] times a
    // window of h zero-padded in front (hp[L + k] = h[k]), so every block is one vector loop.
    let mut hp = [0i16; 2 * L];
    hp[L..].copy_from_slice(&h[..L]);
    for (b, yb) in y[..L].as_chunks_mut::<8>().0.iter_mut().enumerate() {
        let n0 = 8 * b;
        let mut acc = [0i32; 8];
        for (i, &xi) in x[..n0 + 8].iter().enumerate() {
            let w: &[i16; 8] = hp[L + n0 - i..L + n0 - i + 8].try_into().unwrap();
            for (a, &hv) in acc.iter_mut().zip(w) {
                *a = a.wrapping_add(xi as i32 * hv as i32);
            }
        }
        for (yv, a) in yb.iter_mut().zip(acc) {
            *yv = (a >> 12).clamp(i16::MIN as i32, i16::MAX as i32) as i16;
        }
    }
}

/// `Interpol_3or6`: interpolates the normalised correlation `corr` at `corr[lag] + frac/3or6`.
pub(super) fn interpol_3or6(corr: &[i16], lag: usize, mut frac: i16, flag3: bool) -> i16 {
    if flag3 {
        frac <<= 1;
    }
    let mut x = lag;
    if frac < 0 {
        frac += UP_SAMP_MAX;
        x -= 1;
    }
    let c1 = frac as usize;
    let c2 = (UP_SAMP_MAX - frac) as usize;
    // TS 26.073's saturating L_mac chain and round (opencore wraps: saturated correlations
    // times the taps can pass 2^30).
    let ov = &mut false;
    let mut s: i32 = 0;
    for i in 0..L_INTER_SRCH as usize {
        s = l_mac(s, corr[x - i], INTER_6[c1 + 6 * i], ov);
        s = l_mac(s, corr[x + 1 + i], INTER_6[c2 + 6 * i], ov);
    }
    round(s, ov)
}

/// `Norm_Corr`: normalised correlation between `xn` and the filtered past excitation for lags
/// `t_min..=t_max`, written to `corr_norm[lag - t_min]`. `exc[exc_off]` is the C's `exc[0]`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn norm_corr(
    exc: &[i16],
    exc_off: usize,
    xn: &[i16],
    h: &[i16],
    t_min: usize,
    t_max: usize,
    corr_norm: &mut [i16],
    ov: &mut Flag,
) {
    const L: usize = L_SUBFR;
    let mut k = exc_off - t_min;
    let mut excf = [0i16; L];
    convolve(&exc[k..], h, &mut excf);
    let mut scaled = [0i16; L];
    let mut s: i64 = 0; // exact: the C's saturating 2*sum only matters as "> 2^26"
    for j in 0..L {
        scaled[j] = excf[j] >> 2;
        s += excf[j] as i64 * excf[j] as i64;
    }
    let (mut s_excf, h_fac, scaling) = if s <= (67_108_864 >> 1) {
        (excf, 12, 0)
    } else {
        (scaled, 14, 2)
    };
    let xn: &[i16; L] = xn[..L].try_into().unwrap();
    let h: &[i16; L] = h[..L].try_into().unwrap();
    // The fast i32 dots are exact while no partial sum reaches 2^30 (C's L_mac doubles and
    // saturates at 2^31). <e,e> < 2^30 is guaranteed by 40 max|e|^2 < 2^30; then
    // sum |x e| <= sqrt(<x,x> <e,e>) (Cauchy-Schwarz) bounds every partial of <x,e>.
    let xx: i64 = xn.iter().map(|&v| v as i64 * v as i64).sum();
    let hm = max_abs(h) as i32;
    // Every lag on the fast paths (the codec's usual case): `simd::norm_lags` keeps the
    // filtered excitation in registers from lag to lag. Else the loop below, lag by lag.
    let n = t_max + 1 - t_min;
    let mut sums = [[0i32; 2]; 40];
    let lags = NormLags {
        xn,
        h,
        hm,
        h_fac,
        scaling,
        xx,
    };
    if n <= 40 && crate::simd::norm_lags(&lags, &s_excf, &exc[k + 1 - n..k], &mut sums[..n]) {
        for (c, &[s, s2]) in corr_norm[..n].iter_mut().zip(&sums[..n]) {
            *c = norm_of(s << 1, s2 << 1, ov);
        }
        return;
    }
    for i in t_min..=t_max {
        let em = max_abs(&s_excf);
        let fast = if em <= NORM_V_MAX as i64 {
            let (s, s2) = dots(xn, &s_excf);
            ((xx as i128 * s2 as i128) < 1 << 60).then_some((s << 1, s2 << 1))
        } else {
            None
        };
        let (s, s2) = if let Some(f) = fast {
            f
        } else {
            let mut o = false;
            let s2 = s_excf.iter().fold(0, |a, &e| l_mac(a, e, e, &mut o));
            let s = xn
                .iter()
                .zip(&s_excf)
                .fold(0, |a, (&x, &e)| l_mac(a, x, e, &mut o));
            (s, s2)
        };
        corr_norm[i - t_min] = norm_of(s, s2, ov);
        if i != t_max {
            k -= 1;
            next_excf(&mut s_excf, h, hm, exc[k], h_fac, scaling);
        }
    }
}

/// One lag's normalised correlation from <x, v> and <v, v> (both doubled, as the C's L_mac).
#[inline(always)]
fn norm_of(s: i32, s2: i32, ov: &mut Flag) -> i16 {
    let s2 = inv_sqrt(s2, ov);
    let norm_h = (s2 >> 16) as i16;
    let norm_l = ((s2 >> 1) - ((norm_h as i32) << 15)) as i16;
    let corr_h = (s >> 16) as i16;
    let corr_l = ((s >> 1) - ((corr_h as i32) << 15)) as i16;
    // extract_h(L_shl(s, 16)): saturates (opencore truncates).
    mpy_32(corr_h, corr_l, norm_h, norm_l, ov).clamp(i16::MIN as i32, i16::MAX as i32) as i16
}

/// `norm_corr`'s per-call inputs for `simd::norm_lags`: the target, the impulse response, its
/// max |h|, the update's shifts and <x, x>.
pub(crate) struct NormLags<'a> {
    pub(crate) xn: &'a [i16; L_SUBFR],
    pub(crate) h: &'a [i16; L_SUBFR],
    pub(crate) hm: i32,
    pub(crate) h_fac: u32,
    pub(crate) scaling: u32,
    pub(crate) xx: i64,
}

/// Largest |v| with 40 |v|^2 < 2^30: `norm_corr`'s i32 dots are exact up to it.
pub(crate) const NORM_V_MAX: i32 = 5181;
const _: () = {
    let (v, l) = (NORM_V_MAX as i64, L_SUBFR as i64);
    assert!(l * v * v < 1 << 30 && l * (v + 1) * (v + 1) >= 1 << 30);
};

/// `simd::norm_lags`' array form (its non-x86 body and oracle): for lag i = 0.. out.len(),
/// out[i] = (<x, v_i>, <v_i, v_i>) as wrapping i32 dots, v_0 = `v`, v_{i+1} = `next_excf`'s
/// fast update of v_i by exc_rev[len - 1 - i] (the samples run backwards). False, `out`
/// partly written, as soon as a lag leaves the fast paths (max |v| > NORM_V_MAX, <x,x><v,v>
/// >= 2^60, |temp| max|h| >= 2^(h_fac + 15)).
#[cfg_attr(all(target_arch = "x86_64", not(amr_nb_arrays)), allow(dead_code))]
pub(crate) fn norm_lags_portable(
    p: &NormLags,
    v: &[i16; L_SUBFR],
    exc_rev: &[i16],
    out: &mut [[i32; 2]],
) -> bool {
    let mut v = *v;
    let n = out.len();
    for (i, o) in out.iter_mut().enumerate() {
        if max_abs(&v) > NORM_V_MAX as i64 {
            return false;
        }
        let (s, s2) = dots(p.xn, &v);
        if p.xx as i128 * s2 as i128 >= 1 << 60 {
            return false;
        }
        *o = [s, s2];
        if i + 1 < n {
            let temp = exc_rev[exc_rev.len() - 1 - i];
            if (temp as i32).abs() * p.hm >= 1 << (p.h_fac + 15) {
                return false;
            }
            next_excf(&mut v, p.h, p.hm, temp, p.h_fac, p.scaling);
        }
    }
    true
}

/// `norm_corr`'s two correlations (x·v, v·v), as two separate reductions: one loop with two
/// accumulators stays 80 scalar `imul`s, each reduction alone is `pmaddwd`.
#[inline(never)]
fn dots(xn: &[i16; L_SUBFR], v: &[i16; L_SUBFR]) -> (i32, i32) {
    let wsum = |a: i32, b: i32| a.wrapping_add(b);
    let s = xn
        .iter()
        .zip(v)
        .map(|(&x, &e)| x as i32 * e as i32)
        .fold(0, wsum);
    let s2 = v.iter().map(|&e| e as i32 * e as i32).fold(0, wsum);
    (s, s2)
}

/// `norm_corr`'s update to the next lag (not inlined either: same speed, and the lag loop stays
/// small). v[j] = add(extract_h(L_shl(L_mult(temp, h[j]), 15 - h_fac)), v[j - 1]),
/// v[0] = temp >> scaling: saturating like TS 26.073 (opencore wraps). The product clamp is
/// the L_shl saturation: extract_h is then i16::MAX/MIN.
#[inline(never)]
fn next_excf(
    v: &mut [i16; L_SUBFR],
    h: &[i16; L_SUBFR],
    hm: i32,
    temp: i16,
    h_fac: u32,
    scaling: u32,
) {
    let lim = 1i32 << (h_fac + 15);
    if (temp as i32).abs() * hm < lim {
        // No product reaches the L_shl clamp (hm = max |h|): no i32 clamp per lane.
        *v = core::array::from_fn(|j| {
            if j == 0 {
                temp >> scaling
            } else {
                (((temp as i32 * h[j] as i32) >> h_fac) as i16).saturating_add(v[j - 1])
            }
        });
        return;
    }
    *v = core::array::from_fn(|j| {
        if j == 0 {
            temp >> scaling
        } else {
            let p = (temp as i32 * h[j] as i32).clamp(-lim, lim - 1);
            ((p >> h_fac) as i16).saturating_add(v[j - 1])
        }
    });
}

/// `searchFrac`: fractional refinement around `lag` (in `corr`, indexed from `t_min`).
pub(super) fn search_frac(
    lag: &mut i16,
    frac: &mut i16,
    last_frac: i16,
    corr: &[i16],
    t_min: i16,
    flag3: bool,
) {
    let at = (*lag - t_min) as usize;
    let mut max = interpol_3or6(corr, at, *frac, flag3);
    for i in *frac + 1..=last_frac {
        let corr_int = interpol_3or6(corr, at, i, flag3);
        if corr_int > max {
            max = corr_int;
            *frac = i;
        }
    }
    if !flag3 {
        if *frac == -3 {
            *frac = 3;
            *lag -= 1;
        }
    } else {
        if *frac == -2 {
            *frac = 1;
            *lag -= 1;
        }
        if *frac == 2 {
            *frac = -1;
            *lag += 1;
        }
    }
}

/// `getRange`: `(t0_min, t0_max)` search range around `t0`.
pub(super) fn get_range(
    t0: i16,
    delta_low: i16,
    delta_range: i16,
    pitmin: i16,
    pitmax: i16,
) -> (i16, i16) {
    let mut t0_min = t0.wrapping_sub(delta_low).max(pitmin);
    let mut temp = t0_min + delta_range;
    if temp > pitmax {
        temp = pitmax;
        t0_min = pitmax - delta_range;
    }
    (t0_min, temp)
}

/// `Pitch_frState`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct PitchFrState {
    pub t0_prev_subframe: i16,
}

/// Per-mode `Pitch_fr` parameters (`mode_dep_parm`).
struct ModeParm {
    max_frac_lag: i16,
    flag3: bool,
    first_frac: i16,
    last_frac: i16,
    delta_int_low: i16,
    delta_int_range: i16,
    delta_frc_low: i16,
    delta_frc_range: i16,
    pit_min: i16,
}

const fn parm(
    max_frac_lag: i16,
    flag3: bool,
    ff: i16,
    lf: i16,
    dil: i16,
    dir: i16,
    dfl: i16,
    dfr: i16,
    pit_min: i16,
) -> ModeParm {
    ModeParm {
        max_frac_lag,
        flag3,
        first_frac: ff,
        last_frac: lf,
        delta_int_low: dil,
        delta_int_range: dir,
        delta_frc_low: dfl,
        delta_frc_range: dfr,
        pit_min,
    }
}

const MODE_DEP_PARM: [ModeParm; 8] = [
    parm(84, true, -2, 2, 5, 10, 5, 9, 20),
    parm(84, true, -2, 2, 5, 10, 5, 9, 20),
    parm(84, true, -2, 2, 3, 6, 5, 9, 20),
    parm(84, true, -2, 2, 3, 6, 5, 9, 20),
    parm(84, true, -2, 2, 3, 6, 5, 9, 20),
    parm(84, true, -2, 2, 3, 6, 10, 19, 20),
    parm(84, true, -2, 2, 3, 6, 5, 9, 20),
    parm(94, false, -3, 3, 3, 6, 5, 9, 18),
];

/// `Pitch_fr`: closed-loop fractional pitch search for one subframe. Returns
/// `(lag, frac, resu3, index)`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(crate) fn pitch_fr(
    st: &mut PitchFrState,
    mode: Mode,
    t_op: &[i16; 2],
    exc: &[i16],
    exc_off: usize,
    xn: &[i16],
    h: &[i16],
    i_subfr: i16,
    ov: &mut Flag,
) -> (i16, i16, bool, i16) {
    let p = &MODE_DEP_PARM[mode as usize];
    let mut frac = p.first_frac;
    let mut last_frac = p.last_frac;
    let mut delta_search = true;
    let pit_max = PIT_MAX as i16;
    let is_lowrate = matches!(mode, Mode::MR475 | Mode::MR515);
    let (t0_min, t0_max) =
        if (i_subfr == 0 || i_subfr == L_FRAME_BY2) && (!is_lowrate || i_subfr != L_FRAME_BY2) {
            delta_search = false;
            let fo = if i_subfr == 0 { 0 } else { 1 };
            get_range(
                t_op[fo],
                p.delta_int_low,
                p.delta_int_range,
                p.pit_min,
                pit_max,
            )
        } else {
            get_range(
                st.t0_prev_subframe,
                p.delta_frc_low,
                p.delta_frc_range,
                p.pit_min,
                pit_max,
            )
        };
    let t_min = t0_min - L_INTER_SRCH;
    let t_max = t0_max + L_INTER_SRCH;
    let mut corr = [0i16; 40];
    norm_corr(
        exc,
        exc_off,
        xn,
        h,
        t_min as usize,
        t_max as usize,
        &mut corr,
        ov,
    );
    let c = |i: i16| corr[(i - t_min) as usize];
    let mut max = c(t0_min);
    let mut lag = t0_min;
    for i in t0_min + 1..=t0_max {
        if c(i) >= max {
            max = c(i);
            lag = i;
        }
    }
    let delta_modes = matches!(mode, Mode::MR475 | Mode::MR515 | Mode::MR59 | Mode::MR67);
    if !delta_search && lag > p.max_frac_lag {
        frac = 0;
    } else if delta_search && delta_modes {
        let mut tmp_lag = st.t0_prev_subframe;
        if tmp_lag - t0_min > 5 {
            tmp_lag = t0_min + 5;
        }
        if t0_max - tmp_lag > 4 {
            tmp_lag = t0_max - 4;
        }
        if lag == tmp_lag || lag == tmp_lag - 1 {
            search_frac(&mut lag, &mut frac, last_frac, &corr, t_min, p.flag3);
        } else if lag == tmp_lag - 2 {
            frac = 0;
            search_frac(&mut lag, &mut frac, last_frac, &corr, t_min, p.flag3);
        } else if lag == tmp_lag + 1 {
            last_frac = 0;
            search_frac(&mut lag, &mut frac, last_frac, &corr, t_min, p.flag3);
        } else {
            frac = 0;
        }
    } else {
        search_frac(&mut lag, &mut frac, last_frac, &corr, t_min, p.flag3);
    }
    let index = if p.flag3 {
        enc_lag3(
            lag,
            frac,
            st.t0_prev_subframe,
            t0_min,
            t0_max,
            delta_search,
            delta_modes,
        )
    } else {
        enc_lag6(lag, frac, t0_min, delta_search)
    };
    st.t0_prev_subframe = lag;
    (lag, frac, p.flag3, index)
}

/// `Enc_lag3`: encodes a 1/3-resolution pitch lag.
pub(crate) fn enc_lag3(
    t0: i16,
    t0_frac: i16,
    t0_prev: i16,
    t0_min: i16,
    t0_max: i16,
    delta_flag: bool,
    flag4: bool,
) -> i16 {
    if !delta_flag {
        if t0 - 85 <= 0 {
            (t0 << 1) + t0 - 58 + t0_frac
        } else {
            t0 + 112
        }
    } else if !flag4 {
        let i = t0 - t0_min;
        i + (i << 1) + 2 + t0_frac
    } else {
        let mut tmp_lag = t0_prev;
        if tmp_lag - t0_min - 5 > 0 {
            tmp_lag = t0_min + 5;
        }
        if t0_max - tmp_lag - 4 > 0 {
            tmp_lag = t0_max - 4;
        }
        let uplag = t0 + (t0 << 1) + t0_frac;
        let i = tmp_lag - 2;
        let tmp_ind = i + (i << 1);
        if tmp_ind - uplag >= 0 {
            t0 - tmp_lag + 5
        } else {
            let i = tmp_lag + 1;
            let i = i + (i << 1);
            if i > uplag {
                uplag - tmp_ind + 3
            } else {
                t0 - tmp_lag + 11
            }
        }
    }
}

/// `Enc_lag6`: encodes a 1/6-resolution pitch lag (MR122).
pub(crate) fn enc_lag6(t0: i16, t0_frac: i16, t0_min: i16, delta_flag: bool) -> i16 {
    if !delta_flag {
        if t0 <= 94 {
            (t0 << 3) - (t0 << 1) - 105 + t0_frac
        } else {
            t0 + 368
        }
    } else {
        let temp = t0 - t0_min;
        (temp << 3) - (temp << 1) + 3 + t0_frac
    }
}

/// `G_pitch`: adaptive-codebook gain; `g_coeff` receives `<y1,y1>` and `-2<xn,y1>` style terms
/// (mantissa, exponent pairs). Follows TS 26.073 `g_pitch.c`, global `Overflow` included: each
/// sum restarts with `ov` cleared and is redone on `y1 / 4` if it saturated. (opencore's
/// optimised version tests for i32 wrap instead and denormalises the scaled `<xn,y1>` by 4
/// instead of 2: wrong gains on near-full-scale input, `tests/conformance.rs`.)
pub(crate) fn g_pitch(
    mode: Mode,
    xn: &[i16],
    y1: &[i16],
    g_coeff: &mut [i16],
    l_subfr: usize,
    ov: &mut Flag,
) -> i16 {
    /// `s = 1; for i { s = L_mac(s, a[i] >> sa, b[i] >> sb) }`, saturating, setting `ov`.
    /// Fast path (unshifted calls): if 1 + 2 max|a| sum|b| <= MAX_32, no product or partial
    /// sum saturates and the chain is 1 + 2 <a, b> (a plain, vectorising i32 dot).
    fn mac(a: &[i16], sa: u32, b: &[i16], sb: u32, ov: &mut Flag) -> i32 {
        if sa == 0 && sb == 0 && 2 * max_abs(a) * sum_abs(b) < MAX_32 as i64 {
            let d = a
                .iter()
                .zip(b)
                .fold(0i32, |s, (&x, &y)| s.wrapping_add(x as i32 * y as i32));
            return 1 + 2 * d;
        }
        a.iter()
            .zip(b)
            .fold(1, |s, (&x, &y)| l_add(s, l_mult(x >> sa, y >> sb, ov), ov))
    }
    let (xn, y1) = (&xn[..l_subfr], &y1[..l_subfr]);
    *ov = false;
    let mut s = mac(y1, 0, y1, 0, ov);
    let (yy, mut exp_yy);
    if !*ov {
        exp_yy = norm_l(s);
        yy = round(l_shl(s, exp_yy, ov), ov);
    } else {
        s = mac(y1, 2, y1, 2, ov);
        exp_yy = norm_l(s);
        yy = round(l_shl(s, exp_yy, ov), ov);
        exp_yy -= 4;
    }
    *ov = false;
    s = mac(xn, 0, y1, 0, ov);
    let (mut xy, exp_xy);
    if !*ov {
        let e = norm_l(s);
        xy = round(l_shl(s, e, ov), ov);
        exp_xy = e;
    } else {
        s = mac(xn, 0, y1, 2, ov);
        let e = norm_l(s);
        xy = round(l_shl(s, e, ov), ov);
        exp_xy = e - 2;
    }
    g_coeff[0] = yy;
    g_coeff[1] = 15 - exp_yy;
    g_coeff[2] = xy;
    g_coeff[3] = 15 - exp_xy;
    if xy < 4 {
        return 0;
    }
    xy >>= 1;
    let mut gain = div_s(xy, yy);
    gain = shr(gain, exp_xy - exp_yy, ov);
    gain = gain.min(19661);
    if mode == Mode::MR122 {
        gain &= !3;
    }
    gain
}

/// `q_gain_pitch`: scalar pitch-gain quantisation. For MR795 also fills the 3 candidate
/// gains/indices in `cand` (`(gain_cand, gain_cind)`).
pub(crate) fn q_gain_pitch(
    mode: Mode,
    gp_limit: i16,
    gain: &mut i16,
    cand: Option<(&mut [i16; 3], &mut [i16; 3])>,
    ov: &mut Flag,
) -> i16 {
    let q = &QUA_GAIN_PITCH;
    let mut err_min = abs_s(sub(*gain, q[0], ov));
    let mut index = 0usize;
    for i in 1..NB_QUA_PITCH {
        if q[i] <= gp_limit {
            let err = abs_s(sub(*gain, q[i], ov));
            if err < err_min {
                err_min = err;
                index = i;
            }
        }
    }
    if mode == Mode::MR795 {
        let ii = if index == 0 {
            0
        } else if index == NB_QUA_PITCH - 1 || q[index + 1] > gp_limit {
            index - 2
        } else {
            index - 1
        };
        let (gain_cand, gain_cind) = cand.expect("MR795 needs candidate outputs");
        for i in 0..3 {
            gain_cind[i] = (ii + i) as i16;
            gain_cand[i] = q[ii + i];
        }
        *gain = q[index];
    } else if mode == Mode::MR122 {
        *gain = q[index] & !3;
    } else {
        *gain = q[index];
    }
    index as i16
}

/// `tonStabState`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct TonStabState {
    pub gp: [i16; N_FRAME],
    pub count: i16,
}

/// `check_lsp`: detects resonances in the LSPs (tonal signals) for pitch-gain clipping.
pub(crate) fn check_lsp(st: &mut TonStabState, lsp: &[i16; M]) -> bool {
    let mut dist_min1 = MAX_16;
    for i in 3..M - 2 {
        dist_min1 = dist_min1.min(lsp[i].wrapping_sub(lsp[i + 1]));
    }
    let mut dist_min2 = MAX_16;
    for i in 1..3 {
        dist_min2 = dist_min2.min(lsp[i].wrapping_sub(lsp[i + 1]));
    }
    let dist_th = if lsp[1] > 32000 {
        600
    } else if lsp[1] > 30500 {
        800
    } else {
        1100
    };
    if dist_min1 < 1500 || dist_min2 < dist_th {
        st.count = st.count.wrapping_add(1);
    } else {
        st.count = 0;
    }
    if st.count >= 12 {
        st.count = 12;
        true
    } else {
        false
    }
}

/// `check_gp_clipping`.
pub(crate) fn check_gp_clipping(st: &TonStabState, g_pitch: i16, ov: &mut Flag) -> bool {
    let mut sum = shr(g_pitch, 3, ov);
    for &g in &st.gp {
        sum = add(sum, g, ov);
    }
    sum > GP_CLIP
}

/// `update_gp_clipping`.
pub(crate) fn update_gp_clipping(st: &mut TonStabState, g_pitch: i16) {
    st.gp.copy_within(1.., 0);
    st.gp[N_FRAME - 1] = g_pitch >> 3;
}

/// Outputs of `cl_ltp` for one subframe.
pub(crate) struct ClLtpOut {
    pub t0: i16,
    pub t0_frac: i16,
    pub gain_pit: i16,
    pub gp_limit: i16,
}

/// `cl_ltp`: closed-loop pitch search, adaptive excitation (written to `exc[exc_off..]`), its
/// filtered version `yl`, pitch gain, target update `xn2`/`res2`. Appends the lag index (and for
/// MR122 the pitch-gain index) to `ana` at `*ana_pos`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(crate) fn cl_ltp(
    pitch_st: &mut PitchFrState,
    ton_st: &TonStabState,
    mode: Mode,
    frame_offset: i16,
    t_op: &[i16; 2],
    h1: &[i16],
    exc: &mut [i16],
    exc_off: usize,
    res2: &mut [i16; L_SUBFR],
    xn: &[i16; L_SUBFR],
    lsp_flag: bool,
    xn2: &mut [i16; L_SUBFR],
    yl: &mut [i16; L_SUBFR],
    g_coeff: &mut [i16; 6],
    ana: &mut [i16],
    ana_pos: &mut usize,
    ov: &mut Flag,
) -> ClLtpOut {
    let (t0, t0_frac, resu3, index) =
        pitch_fr(pitch_st, mode, t_op, exc, exc_off, xn, h1, frame_offset, ov);
    ana[*ana_pos] = index;
    *ana_pos += 1;
    pred_lt_3or6(exc, exc_off, t0, t0_frac, L_SUBFR as i16, resu3 as i16, ov);
    convolve(&exc[exc_off..], h1, yl);
    let mut gain_pit = g_pitch(mode, xn, yl, g_coeff, L_SUBFR, ov);
    let mut gpc_flag = false;
    let mut gp_limit = MAX_16;
    if lsp_flag && gain_pit > GP_CLIP {
        gpc_flag = check_gp_clipping(ton_st, gain_pit, ov);
    }
    if matches!(mode, Mode::MR475 | Mode::MR515) {
        gain_pit = gain_pit.min(13926);
        if gpc_flag {
            gp_limit = GP_CLIP;
        }
    } else {
        if gpc_flag {
            gp_limit = GP_CLIP;
            gain_pit = GP_CLIP;
        }
        if mode == Mode::MR122 {
            ana[*ana_pos] = q_gain_pitch(Mode::MR122, gp_limit, &mut gain_pit, None, ov);
            *ana_pos += 1;
        }
    }
    // sub(v, extract_h(L_shl(L_mult(u, gain_pit), 1))), saturating like TS 26.073: the clamp
    // to +-2^29 is the L_shl saturation (opencore truncates and wraps).
    let scaled =
        |u: i16| ((u as i32 * gain_pit as i32).clamp(-(1 << 29), (1 << 29) - 1) >> 14) as i16;
    for i in 0..L_SUBFR {
        xn2[i] = xn[i].saturating_sub(scaled(yl[i]));
        res2[i] = res2[i].saturating_sub(scaled(exc[exc_off + i]));
    }
    ClLtpOut {
        t0,
        t0_frac,
        gain_pit,
        gp_limit,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::{N, Rng};

    /// `simd::norm_lags` is `norm_lags_portable`: same verdict, same sums up to it, on 1 to 40
    /// lags with excitations, targets and impulse responses from small to full scale (the
    /// range, energy and update limits all trip).
    #[test]
    fn norm_lags_is_the_array_form() {
        let mut g = Rng::new(19);
        let mut fast = 0;
        for case in 0..N {
            let sc = |g: &mut Rng, lo: u64| (g.u64() % 16).max(lo) as u32;
            let (sx, sv, sh, se) = (sc(&mut g, 0), sc(&mut g, 2), sc(&mut g, 0), sc(&mut g, 0));
            let xn: [i16; L_SUBFR] = core::array::from_fn(|_| g.i16() >> sx);
            let h: [i16; L_SUBFR] = core::array::from_fn(|_| g.i16() >> sh);
            let v: [i16; L_SUBFR] = core::array::from_fn(|_| g.i16() >> sv);
            let exc: [i16; 40] = core::array::from_fn(|_| g.i16() >> se);
            let (h_fac, scaling) = if case % 2 == 0 { (12, 0) } else { (14, 2) };
            let xx = xn.iter().map(|&v| v as i64 * v as i64).sum();
            let p = NormLags {
                xn: &xn,
                h: &h,
                hm: max_abs(&h) as i32,
                h_fac,
                scaling,
                xx,
            };
            let n = 1 + (g.u64() % 40) as usize;
            let (mut a, mut b) = ([[0i32; 2]; 40], [[0i32; 2]; 40]);
            let ra = crate::simd::norm_lags(&p, &v, &exc[..n - 1], &mut a[..n]);
            let rb = norm_lags_portable(&p, &v, &exc[..n - 1], &mut b[..n]);
            assert_eq!(ra, rb, "case {case}");
            if ra {
                fast += 1;
                assert_eq!(a, b, "case {case}");
            }
        }
        assert!(fast > N / 8 && fast < N - N / 8, "{fast} of {N} fast");
    }

    /// `dots` and `next_excf` vs the C's loops as written, on full-scale inputs (i32 sums and
    /// i16 adds wrap) and small ones. `pitch_fr_matches_c` checks the whole search vs the C.
    #[test]
    fn norm_corr_steps_are_the_c_loops() {
        let mut g = Rng::new(7);
        for case in 0..N / 10 {
            let shift = if case % 3 == 0 {
                0
            } else {
                (g.u64() % 15) as u32
            };
            let mut v = [0i16; L_SUBFR];
            let mut xn = [0i16; L_SUBFR];
            let mut h = [0i16; L_SUBFR];
            for ((a, b), c) in v.iter_mut().zip(&mut xn).zip(&mut h) {
                (*a, *b, *c) = (g.i16() >> shift, g.i16() >> shift, g.i16() >> shift);
            }
            let (mut s, mut s2) = (0i32, 0i32);
            for j in 0..L_SUBFR {
                s = s.wrapping_add(xn[j] as i32 * v[j] as i32);
                s2 = s2.wrapping_add(v[j] as i32 * v[j] as i32);
            }
            assert_eq!(dots(&xn, &v), (s, s2));
            let temp = g.i16();
            let (h_fac, scaling) = if case & 1 == 0 { (12, 0) } else { (14, 2) };
            // TS 26.073: add(extract_h(L_shl(L_mult(temp, h[j]), 15 - h_fac)), v[j - 1]).
            let o = &mut false;
            let mut want = v;
            for j in (1..L_SUBFR).rev() {
                let s = (l_shl(l_mult(temp, h[j], o), 15 - h_fac as i16, o) >> 16) as i16;
                want[j] = add(s, want[j - 1], o);
            }
            want[0] = temp >> scaling;
            next_excf(&mut v, &h, max_abs(&h) as i32, temp, h_fac, scaling);
            assert_eq!(v, want);
        }
    }
}
