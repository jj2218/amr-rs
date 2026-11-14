//! Gain decoding, prediction, error concealment and background-noise control.
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0): `common/src/{gc_pred,gmed_n}.cpp` and
//! `dec/src/{d_gain_p,d_gain_c,dec_gain,ec_gains,c_g_aver,ex_ctrl,bgnscd}.cpp`.
#![allow(dead_code)] // integrated into the decoder later
#![allow(clippy::needless_range_loop)] // keep the C loop/index structure

use crate::basic_op::{
    Flag, MAX_16, MAX_32, MIN_16, MIN_32, add, div_s, l_add, l_add_energy, l_mac, l_msu, l_mult,
    l_shl, l_shr, l_sub, mpy_32_16, mult, norm_l, norm_s, round, shl, shr, shr_r, sub,
};
use crate::math::{log2, log2_norm, pow2};
use crate::tables::{
    QUA_GAIN_CODE, QUA_GAIN_PITCH, TABLE_GAIN_HIGHRATES, TABLE_GAIN_LOWRATES, TABLE_GAIN_MR475,
};

/// C `enum Mode` (passed as an `int`). Local until a shared mode type exists.
pub(crate) use crate::cnst::Mode::{self, *};

use crate::cnst::{L_FRAME, L_SUBFR, M};

// ---------------------------------------------------------------------------------------------
// gmed_n.cpp
// ---------------------------------------------------------------------------------------------

const NMAX: usize = 9;

/// `gmed_n`: median of `ind[..n]` (n <= 9, odd in practice).
pub(crate) fn gmed_n(ind: &[i16], n: i16) -> i16 {
    let n = n as usize;
    let mut tmp = [0i16; NMAX];
    let mut tmp2 = [0i16; NMAX];
    tmp2[..n].copy_from_slice(&ind[..n]);
    // `ix` is declared outside the loops in C and keeps its value between iterations.
    let mut ix: i16 = 0;
    for i in 0..n {
        let mut max: i16 = -32767;
        for j in 0..n {
            if tmp2[j] >= max {
                max = tmp2[j];
                ix = j as i16;
            }
        }
        tmp2[ix as usize] = -32768;
        tmp[i] = ix;
    }
    let median_index = tmp[n >> 1];
    ind[median_index as usize]
}

// ---------------------------------------------------------------------------------------------
// gc_pred.cpp
// ---------------------------------------------------------------------------------------------

const NPRED: usize = 4;
const MEAN_ENER_MR122: i32 = 783741;
const MIN_ENERGY: i16 = -14336;
const MIN_ENERGY_MR122: i16 = -2381;
const PRED: [i16; NPRED] = [5571, 4751, 2785, 1556];
const PRED_MR122: [i16; NPRED] = [44, 37, 22, 12];

/// `gc_predState`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct GcPredState {
    pub past_qua_en: [i16; 4],
    pub past_qua_en_mr122: [i16; 4],
}

impl Default for GcPredState {
    fn default() -> Self {
        Self {
            past_qua_en: [MIN_ENERGY; 4],
            past_qua_en_mr122: [MIN_ENERGY_MR122; 4],
        }
    }
}

impl GcPredState {
    /// `gc_pred_reset`
    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }
}

/// `gc_pred`: returns `(exp_gcode0, frac_gcode0, exp_en, frac_en)`.
/// `exp_en`/`frac_en` are only computed for MR795 (C leaves them untouched otherwise; 0 here).
pub(crate) fn gc_pred(
    st: &GcPredState,
    mode: Mode,
    code: &[i16],
    ov: &mut Flag,
) -> (i16, i16, i16, i16) {
    let (mut exp_en, mut frac_en) = (0i16, 0i16);

    // TS 26.073's `L_mac` chain over code² (terms >= 0): min(2 sum, MAX_32). (opencore sums
    // c² >> 3 and wraps on `<< 4`: wrong for large codebook energies.)
    let sum: i64 = code[..L_SUBFR]
        .iter()
        .map(|&c| 2 * c as i64 * c as i64)
        .sum();
    let mut ener_code = sum.min(MAX_32 as i64) as i32;

    let (exp_gcode0, frac_gcode0) = if mode == MR122 {
        ener_code = (round(ener_code, ov) as i32 * 26214) << 1;
        let (exp, frac) = log2(ener_code, ov);
        let mut l_temp1 = ((exp as i32) - 30).wrapping_shl(16);
        ener_code = l_temp1.wrapping_add((frac as i32) << 1);
        let mut ener = MEAN_ENER_MR122;
        for i in 0..NPRED {
            l_temp1 = (st.past_qua_en_mr122[i] as i32 * PRED_MR122[i] as i32) << 1;
            ener = l_add(ener, l_temp1, ov);
        }
        l_temp1 = l_sub(ener, ener_code, ov);
        let exp_gcode0 = (l_temp1 >> 17) as i16;
        let l_temp2 = (exp_gcode0 as i32) << 15;
        l_temp1 >>= 2;
        (exp_gcode0, l_temp1.wrapping_sub(l_temp2) as i16)
    } else {
        let exp_code = norm_l(ener_code);
        ener_code = l_shl(ener_code, exp_code, ov);
        let (exp, frac) = log2_norm(ener_code, exp_code);
        let mut l_temp2 = (exp as i32 * -24660) << 1;
        let mut l_tmp = (frac as i32 * -24660) >> 15;
        if l_tmp & 0x0001_0000 != 0 {
            l_tmp |= 0xffff_0000u32 as i32;
        }
        l_tmp <<= 1;
        l_tmp = l_add(l_tmp, l_temp2, ov);
        if mode == MR102 {
            l_temp2 = 16678 << 7;
            l_tmp = l_add(l_tmp, l_temp2, ov);
        } else if mode == MR795 {
            frac_en = (ener_code >> 16) as i16;
            exp_en = -11 - exp_code;
            l_temp2 = 17062 << 7;
            l_tmp = l_add(l_tmp, l_temp2, ov);
        } else if mode == MR74 {
            l_temp2 = 32588 << 6;
            l_tmp = l_add(l_tmp, l_temp2, ov);
        } else if mode == MR67 {
            l_temp2 = 32268 << 6;
            l_tmp = l_add(l_tmp, l_temp2, ov);
        } else {
            l_temp2 = 16678 << 7;
            l_tmp = l_add(l_tmp, l_temp2, ov);
        }
        if l_tmp > 0x001f_ffff {
            *ov = true;
            l_tmp = MAX_32;
        } else if l_tmp < -2097152 {
            *ov = true;
            l_tmp = MIN_32;
        } else {
            l_tmp <<= 10;
        }
        for i in 0..4 {
            l_temp2 = (PRED[i] as i32 * st.past_qua_en[i] as i32) << 1;
            l_tmp = l_add(l_tmp, l_temp2, ov);
        }
        let gcode0 = (l_tmp >> 16) as i16;
        if mode == MR74 {
            l_tmp = (gcode0 as i32 * 5439) << 1;
        } else {
            l_tmp = (gcode0 as i32 * 5443) << 1;
        }
        if l_tmp < 0 {
            l_tmp = !((!l_tmp) >> 8);
        } else {
            l_tmp >>= 8;
        }
        let exp_gcode0 = (l_tmp >> 16) as i16;
        let l_temp1 = if l_tmp < 0 {
            !((!l_tmp) >> 1)
        } else {
            l_tmp >> 1
        };
        l_temp2 = (exp_gcode0 as i32) << 15;
        (exp_gcode0, l_sub(l_temp1, l_temp2, ov) as i16)
    };
    (exp_gcode0, frac_gcode0, exp_en, frac_en)
}

/// `gc_pred_update`
pub(crate) fn gc_pred_update(st: &mut GcPredState, qua_ener_mr122: i16, qua_ener: i16) {
    st.past_qua_en[3] = st.past_qua_en[2];
    st.past_qua_en_mr122[3] = st.past_qua_en_mr122[2];
    st.past_qua_en[2] = st.past_qua_en[1];
    st.past_qua_en_mr122[2] = st.past_qua_en_mr122[1];
    st.past_qua_en[1] = st.past_qua_en[0];
    st.past_qua_en_mr122[1] = st.past_qua_en_mr122[0];
    st.past_qua_en_mr122[0] = qua_ener_mr122;
    st.past_qua_en[0] = qua_ener;
}

/// `(av >> 2) | 0xc000` truncated back to Word16, as in C for negative sums.
fn quarter(av: i16) -> i16 {
    if av < 0 {
        ((av >> 2) as i32 | 0xc000) as i16
    } else {
        av >> 2
    }
}

/// `gc_pred_average_limited`: returns `(ener_avg_MR122, ener_avg)`.
pub(crate) fn gc_pred_average_limited(st: &GcPredState, ov: &mut Flag) -> (i16, i16) {
    let mut av_pred_en: i16 = 0;
    for i in 0..NPRED {
        av_pred_en = add(av_pred_en, st.past_qua_en_mr122[i], ov);
    }
    av_pred_en = quarter(av_pred_en);
    if av_pred_en < MIN_ENERGY_MR122 {
        av_pred_en = MIN_ENERGY_MR122;
    }
    let ener_avg_mr122 = av_pred_en;

    av_pred_en = 0;
    for i in 0..NPRED {
        av_pred_en = add(av_pred_en, st.past_qua_en[i], ov);
    }
    av_pred_en = quarter(av_pred_en);
    if av_pred_en < MIN_ENERGY {
        av_pred_en = MIN_ENERGY;
    }
    (ener_avg_mr122, av_pred_en)
}

// ---------------------------------------------------------------------------------------------
// d_gain_p.cpp, d_gain_c.cpp
// ---------------------------------------------------------------------------------------------

/// `d_gain_pitch` (table: `qua_gain_pitch`; `index` < 16).
pub(crate) fn d_gain_pitch(mode: Mode, index: i16) -> i16 {
    let mut gain = QUA_GAIN_PITCH[index as usize];
    if mode == MR122 {
        gain &= 0xFFFCu16 as i16;
    }
    gain
}

/// `d_gain_code` (table: `qua_gain_code`): returns `gain_code`.
pub(crate) fn d_gain_code(
    pred_state: &mut GcPredState,
    mode: Mode,
    index: i16,
    code: &[i16],
    ov: &mut Flag,
) -> i16 {
    let (exp, frac, _exp_inn_en, _frac_inn_en) = gc_pred(pred_state, mode, code, ov);
    let index = index & 31;
    let tbl_tmp = (index + (index << 1)) as usize;
    let p = &QUA_GAIN_CODE[tbl_tmp..tbl_tmp + 3];
    let temp = sub(mode as i16, MR122 as i16, ov);
    let gain_code = if temp == 0 {
        let mut gcode0 = pow2(exp, frac, ov) as i16;
        gcode0 = shl(gcode0, 4, ov);
        shl(mult(gcode0, p[0], ov), 1, ov)
    } else {
        let gcode0 = pow2(14, frac, ov) as i16;
        let mut l_tmp = l_mult(p[0], gcode0, ov);
        l_tmp = l_shr(l_tmp, sub(9, exp, ov), ov);
        (l_tmp >> 16) as i16
    };
    gc_pred_update(pred_state, p[1], p[2]);
    gain_code
}

// ---------------------------------------------------------------------------------------------
// dec_gain.cpp
// ---------------------------------------------------------------------------------------------

const MR475_VQ_SIZE: i16 = 256;

/// `Dec_gain`: returns `(gain_pit, gain_cod)`. Modes MR475/515/59/67/74/102 only (C passes
/// NULL `exp_en` to `gc_pred`, so MR795 would crash there). Tables come from `tables.rs`.
pub(crate) fn dec_gain(
    pred_state: &mut GcPredState,
    mode: Mode,
    index: i16,
    code: &[i16],
    even_subfr: i16,
    ov: &mut Flag,
) -> (i16, i16) {
    let gain_pit;
    let g_code;
    let qua_ener_mr122;
    let qua_ener;
    let mut index = shl(index, 2, ov);
    if mode == MR102 || mode == MR74 || mode == MR67 {
        let p = &TABLE_GAIN_HIGHRATES[index as usize..];
        gain_pit = p[0];
        g_code = p[1];
        qua_ener_mr122 = p[2];
        qua_ener = p[3];
    } else if mode == MR475 {
        index = index.wrapping_add((1 ^ even_subfr) << 1);
        if index > MR475_VQ_SIZE * 4 - 2 {
            index = MR475_VQ_SIZE * 4 - 2;
        }
        let p = &TABLE_GAIN_MR475[index as usize..];
        gain_pit = p[0];
        g_code = p[1];
        let (mut exp, frac) = log2(g_code as i32, ov);
        exp = exp.wrapping_sub(12);
        let temp1 = shr_r(frac, 5, ov);
        let temp2 = shl(exp, 10, ov);
        qua_ener_mr122 = add(temp1, temp2, ov);
        let mut l_tmp = mpy_32_16(exp, frac, 24660, ov);
        l_tmp = l_shl(l_tmp, 13, ov);
        qua_ener = round(l_tmp, ov);
    } else {
        let p = &TABLE_GAIN_LOWRATES[index as usize..];
        gain_pit = p[0];
        g_code = p[1];
        qua_ener_mr122 = p[2];
        qua_ener = p[3];
    }
    let (exp, frac, _, _) = gc_pred(pred_state, mode, code, ov);
    let gcode0 = pow2(14, frac, ov) as i16;
    let mut l_tmp = l_mult(g_code, gcode0, ov);
    let temp1 = 10i16.wrapping_sub(exp);
    l_tmp = l_shr(l_tmp, temp1, ov);
    let gain_cod = (l_tmp >> 16) as i16;
    gc_pred_update(pred_state, qua_ener_mr122, qua_ener);
    (gain_pit, gain_cod)
}

// ---------------------------------------------------------------------------------------------
// ec_gains.cpp
// ---------------------------------------------------------------------------------------------

/// `ec_gain_pitchState`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct EcGainPitchState {
    pub pbuf: [i16; 5],
    pub past_gain_pit: i16,
    pub prev_gp: i16,
}

impl Default for EcGainPitchState {
    fn default() -> Self {
        Self {
            pbuf: [1640; 5],
            past_gain_pit: 0,
            prev_gp: 16384,
        }
    }
}

impl EcGainPitchState {
    /// `ec_gain_pitch_reset`
    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }
}

/// `ec_gain_codeState`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct EcGainCodeState {
    pub gbuf: [i16; 5],
    pub past_gain_code: i16,
    pub prev_gc: i16,
}

impl Default for EcGainCodeState {
    fn default() -> Self {
        Self {
            gbuf: [1; 5],
            past_gain_code: 0,
            prev_gc: 1,
        }
    }
}

impl EcGainCodeState {
    /// `ec_gain_code_reset`
    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }
}

/// `ec_gain_code`: returns `gain_code`. `state` in 0..=6.
pub(crate) fn ec_gain_code(
    st: &EcGainCodeState,
    pred_state: &mut GcPredState,
    state: i16,
    ov: &mut Flag,
) -> i16 {
    const CDOWN: [i16; 7] = [32767, 32112, 32112, 32112, 32112, 32112, 22937];
    let mut tmp = gmed_n(&st.gbuf, 5);
    if sub(tmp, st.past_gain_code, ov) > 0 {
        tmp = st.past_gain_code;
    }
    tmp = mult(tmp, CDOWN[state as usize], ov);
    let (qua_ener_mr122, qua_ener) = gc_pred_average_limited(pred_state, ov);
    gc_pred_update(pred_state, qua_ener_mr122, qua_ener);
    tmp
}

/// `ec_gain_code_update` (`gain_code` is in/out).
pub(crate) fn ec_gain_code_update(
    st: &mut EcGainCodeState,
    bfi: i16,
    prev_bf: i16,
    gain_code: &mut i16,
    ov: &mut Flag,
) {
    if bfi == 0 {
        if prev_bf != 0 && sub(*gain_code, st.prev_gc, ov) > 0 {
            *gain_code = st.prev_gc;
        }
        st.prev_gc = *gain_code;
    }
    st.past_gain_code = *gain_code;
    for i in 1..5 {
        st.gbuf[i - 1] = st.gbuf[i];
    }
    st.gbuf[4] = *gain_code;
}

/// `ec_gain_pitch`: returns `gain_pitch`. `state` in 0..=6.
pub(crate) fn ec_gain_pitch(st: &EcGainPitchState, state: i16, ov: &mut Flag) -> i16 {
    const PDOWN: [i16; 7] = [32767, 32112, 32112, 26214, 9830, 6553, 6553];
    let mut tmp = gmed_n(&st.pbuf, 5);
    if sub(tmp, st.past_gain_pit, ov) > 0 {
        tmp = st.past_gain_pit;
    }
    mult(tmp, PDOWN[state as usize], ov)
}

/// `ec_gain_pitch_update` (`gain_pitch` is in/out).
pub(crate) fn ec_gain_pitch_update(
    st: &mut EcGainPitchState,
    bfi: i16,
    prev_bf: i16,
    gain_pitch: &mut i16,
    ov: &mut Flag,
) {
    if bfi == 0 {
        if prev_bf != 0 && sub(*gain_pitch, st.prev_gp, ov) > 0 {
            *gain_pitch = st.prev_gp;
        }
        st.prev_gp = *gain_pitch;
    }
    st.past_gain_pit = *gain_pitch;
    if sub(st.past_gain_pit, 16384, ov) > 0 {
        st.past_gain_pit = 16384;
    }
    for i in 1..5 {
        st.pbuf[i - 1] = st.pbuf[i];
    }
    st.pbuf[4] = st.past_gain_pit;
}

// ---------------------------------------------------------------------------------------------
// c_g_aver.cpp
// ---------------------------------------------------------------------------------------------

const L_CBGAINHIST: usize = 7;

/// `Cb_gain_averageState`
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct CbGainAverageState {
    pub cb_gain_history: [i16; L_CBGAINHIST],
    pub hang_var: i16,
    pub hang_count: i16,
}

impl CbGainAverageState {
    /// `Cb_gain_average_reset`
    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }
}

/// The C's `diff` loop in `Cb_gain_average`: per coefficient `sub`, `abs_s`, `norm_s`, `shl`,
/// `div_s`, `shr`/`shl`, summed with `add`. The 10 terms are independent, so they run as
/// branchless lanes (the C's form: `lsp_dev_sum_ref` in the tests). Exact:
/// - t = abs_s(sub(aver, lsp)) = min(|aver − lsp|, MAX_16) (`sub` sets ov past i16).
/// - aver ≤ 0 (`div_s` of a non-positive denominator) or t = 0: the term is 0.
/// - Else, with e1 = ⌊log2 t⌋ and e2 = ⌊log2 aver⌋ (f32 exponents, exact below 2^24):
///   shl(t, norm_s(t) − 1) = n, the top 14 bits of t (a right shift by 1 at e1 = 14), and
///   shl(aver, norm_s(aver)) = d = aver·2^(14 − e2), so n < 2^14 ≤ d and `div_s` is
///   q = ⌊2^15 n / d⌋ ∈ [2^13, 2^15). One f64 divide is exact: the quotient's distance to the
///   next integer (≥ 1/d) is far above f64's rounding error.
/// - The shift 2 + shift1 − shift2 = 1 + e2 − e1 ∈ [−13, 15], so the term is ⌊q·2^(e1 − e2 − 1)⌋
///   (an exact f32 product), clamped at MAX_16 (`shl`'s saturation; q > 0).
/// - The terms are ≥ 0: the `add` chain is min(Σ, MAX_16), setting ov iff Σ > MAX_16.
fn lsp_dev_sum(aver: &[i16; M], lsp: &[i16; M], ov: &mut Flag) -> i16 {
    // 12 lanes (whole SSE2 vectors); the pad lanes have aver = lsp = 0, so they add 0
    const L: usize = 12;
    // Float → int without `as` (its saturation fixups are scalar per lane): floor by the 2^52
    // magic add (round to nearest) and a compare, the integer from the low bits of x + 2^52.
    const MAGIC: f64 = 4_503_599_627_370_496.0; // 2^52; exact for 0 <= x < 2^31
    let (mut av, mut ls) = ([0i32; L], [0i32; L]);
    for i in 0..M {
        (av[i], ls[i]) = (aver[i] as i32, lsp[i] as i32);
    }
    let mut r = [0i32; L];
    let mut sub_ov = 0i32;
    for i in 0..L {
        let d = av[i] - ls[i];
        sub_ov |= ((d + 0x8000) as u32 > 0xffff) as i32;
        let t = d.abs().min(MAX_16 as i32);
        let b1 = (t as f32).to_bits() as i32;
        let b2 = (av[i].max(1) as f32).to_bits() as i32;
        let n = ((b1 >> 10) & 0x1fff) | 0x2000;
        let den = ((b2 >> 9) & 0x3fff) | 0x4000;
        // 2^(e1 − e2 − 1) from the exponent fields, (e1 + 127) − (e2 + 127) − 1 + 127;
        // clamped for the masked lanes (live ones are in [−15, 13])
        let k = ((b1 >> 23) - (b2 >> 23)).clamp(-14, 14);
        let p = f32::from_bits(((k + 126) << 23) as u32) as f64;
        let x = (n << 15) as f64 / den as f64;
        let m = (x + MAGIC) - MAGIC;
        let q = m - if m > x { 1.0 } else { 0.0 };
        let y = q * p;
        let m = (y + MAGIC) - MAGIC;
        let v = m - if m > y { 1.0 } else { 0.0 };
        let v = ((v + MAGIC).to_bits() as i32).min(MAX_16 as i32);
        r[i] = if av[i] > 0 && t > 0 { v } else { 0 };
    }
    let sum: i32 = r.iter().sum();
    *ov |= (sub_ov != 0) | (sum > MAX_16 as i32);
    sum.min(MAX_16 as i32) as i16
}

/// `Cb_gain_average`: returns the mixed codebook gain.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cb_gain_average(
    st: &mut CbGainAverageState,
    mode: Mode,
    gain_code: i16,
    lsp: &[i16; M],
    lsp_aver: &[i16; M],
    bfi: i16,
    prev_bf: i16,
    pdfi: i16,
    prev_pdf: i16,
    in_background_noise: i16,
    voiced_hangover: i16,
    ov: &mut Flag,
) -> i16 {
    let mut cb_gain_mix = gain_code;

    for i in 0..L_CBGAINHIST - 1 {
        st.cb_gain_history[i] = st.cb_gain_history[i + 1];
    }
    st.cb_gain_history[L_CBGAINHIST - 1] = gain_code;

    let diff = lsp_dev_sum(lsp_aver, lsp, ov);

    if diff > 5325 {
        st.hang_var = st.hang_var.wrapping_add(1);
    } else {
        st.hang_var = 0;
    }
    if st.hang_var > 10 {
        st.hang_count = 0;
    }

    if mode <= MR67 || mode == MR102 {
        let low_rate = mode == MR475 || mode == MR515 || mode == MR59;
        let tmp_diff = if ((pdfi != 0 && prev_pdf != 0) || bfi != 0 || prev_bf != 0)
            && voiced_hangover > 1
            && in_background_noise != 0
            && low_rate
        {
            diff.wrapping_sub(4506)
        } else {
            diff.wrapping_sub(3277)
        };
        let tmp1 = if tmp_diff > 0 { tmp_diff } else { 0 };
        let mut bg_mix = if 2048 < tmp1 { 8192 } else { shl(tmp1, 2, ov) };
        if st.hang_count < 40 || diff > 5325 {
            bg_mix = 8192;
        }

        let mut l_sum = l_mult(6554, st.cb_gain_history[2], ov);
        for i in 3..L_CBGAINHIST {
            l_sum = l_mac(l_sum, 6554, st.cb_gain_history[i], ov);
        }
        let mut cb_gain_mean = round(l_sum, ov);

        if (bfi != 0 || prev_bf != 0) && in_background_noise != 0 && low_rate {
            l_sum = l_mult(4681, st.cb_gain_history[0], ov);
            for i in 1..L_CBGAINHIST {
                l_sum = l_mac(l_sum, 4681, st.cb_gain_history[i], ov);
            }
            cb_gain_mean = round(l_sum, ov);
        }

        l_sum = l_mult(bg_mix, cb_gain_mix, ov);
        l_sum = l_mac(l_sum, 8192, cb_gain_mean, ov);
        l_sum = l_msu(l_sum, bg_mix, cb_gain_mean, ov);
        cb_gain_mix = round(l_shl(l_sum, 2, ov), ov);
    }

    st.hang_count = st.hang_count.wrapping_add(1);
    cb_gain_mix
}

// ---------------------------------------------------------------------------------------------
// ex_ctrl.cpp
// ---------------------------------------------------------------------------------------------

/// `Ex_ctrl`: scales `excitation[..40]` in place; always returns 0.
pub(crate) fn ex_ctrl(
    excitation: &mut [i16],
    exc_energy: i16,
    ex_energy_hist: &[i16],
    voiced_hangover: i16,
    prev_bfi: i16,
    careful_flag: i16,
    ov: &mut Flag,
) -> i16 {
    let mut exc_energy = exc_energy;
    let mut avg_energy = gmed_n(ex_energy_hist, 9);
    // shr(add(..), 1): the add saturates (TS 26.073; opencore sums in int).
    let mut prev_energy = shr(add(ex_energy_hist[7], ex_energy_hist[8], ov), 1, ov);
    if ex_energy_hist[8] < prev_energy {
        prev_energy = ex_energy_hist[8];
    }
    if exc_energy < avg_energy && exc_energy > 5 {
        let mut test_energy = shl(prev_energy, 2, ov);
        if voiced_hangover < 7 || prev_bfi != 0 {
            test_energy = sub(test_energy, prev_energy, ov);
        }
        if avg_energy > test_energy {
            avg_energy = test_energy;
        }
        let exp = norm_s(exc_energy);
        exc_energy = shl(exc_energy, exp, ov);
        exc_energy = div_s(16383, exc_energy);
        let mut t0 = l_mult(avg_energy, exc_energy, ov);
        t0 = l_shr(t0, sub(20, exp, ov), ov);
        if t0 > 32767 {
            t0 = 32767;
        }
        let mut scale_factor = t0 as i16;
        if careful_flag != 0 && scale_factor > 3072 {
            scale_factor = 3072;
        }
        for e in excitation[..L_SUBFR].iter_mut() {
            t0 = l_mult(scale_factor, *e, ov);
            t0 = l_shr(t0, 11, ov);
            *e = t0 as i16;
        }
    }
    0
}

// ---------------------------------------------------------------------------------------------
// bgnscd.cpp
// ---------------------------------------------------------------------------------------------

const L_ENERGYHIST: usize = 60;
const FRAMEENERGYLIMIT: i16 = 17578;
const LOWERNOISELIMIT: i16 = 20;
const UPPERNOISELIMIT: i16 = 1953;

/// `Bgn_scdState`
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct BgnScdState {
    pub frame_energy_hist: [i16; L_ENERGYHIST],
    pub bg_hangover: i16,
}

impl Default for BgnScdState {
    fn default() -> Self {
        Self {
            frame_energy_hist: [0; L_ENERGYHIST],
            bg_hangover: 0,
        }
    }
}

impl BgnScdState {
    /// `Bgn_scd_reset`
    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }
}

/// `Bgn_scd`: background-noise detection on one frame (`speech[..160]`, `ltp_gain_hist[..9]`).
/// Returns `inbgNoise` (0/1); `voiced_hangover` is in/out.
pub(crate) fn bgn_scd(
    st: &mut BgnScdState,
    ltp_gain_hist: &[i16],
    speech: &[i16],
    voiced_hangover: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let s = l_add_energy(&speech[..L_FRAME], ov);
    let curr_energy = if s > 0x1fff_ffff {
        MAX_16
    } else {
        (s >> 14) as i16
    };

    let mut frame_energy_min: i16 = 32767;
    for i in (0..L_ENERGYHIST).rev() {
        if st.frame_energy_hist[i] < frame_energy_min {
            frame_energy_min = st.frame_energy_hist[i];
        }
    }
    let l_temp = (frame_energy_min as i32) << 4;
    let noise_floor = if l_temp != (l_temp as i16) as i32 {
        if l_temp > 0 { MAX_16 } else { MIN_16 }
    } else {
        l_temp as i16
    };

    let mut max_energy = st.frame_energy_hist[0];
    for i in (1..=L_ENERGYHIST - 5).rev() {
        if max_energy < st.frame_energy_hist[i] {
            max_energy = st.frame_energy_hist[i];
        }
    }
    let mut max_energy_last_part = st.frame_energy_hist[2 * L_ENERGYHIST / 3];
    for i in 2 * L_ENERGYHIST / 3 + 1..L_ENERGYHIST {
        if max_energy_last_part < st.frame_energy_hist[i] {
            max_energy_last_part = st.frame_energy_hist[i];
        }
    }

    if max_energy > LOWERNOISELIMIT
        && curr_energy < FRAMEENERGYLIMIT
        && curr_energy > LOWERNOISELIMIT
        && (curr_energy < noise_floor || max_energy_last_part < UPPERNOISELIMIT)
    {
        if st.bg_hangover as i32 + 1 > 30 {
            st.bg_hangover = 30;
        } else {
            st.bg_hangover += 1;
        }
    } else {
        st.bg_hangover = 0;
    }

    let inbg_noise: i16 = if st.bg_hangover > 1 { 1 } else { 0 };

    for i in 0..L_ENERGYHIST - 1 {
        st.frame_energy_hist[i] = st.frame_energy_hist[i + 1];
    }
    st.frame_energy_hist[L_ENERGYHIST - 1] = curr_energy;

    let ltp_limit: i16 = if st.bg_hangover > 15 {
        16383
    } else if st.bg_hangover > 8 {
        15565
    } else {
        13926
    };

    let mut prev_voiced = gmed_n(&ltp_gain_hist[4..], 5) > ltp_limit;
    if st.bg_hangover > 20 {
        prev_voiced = gmed_n(ltp_gain_hist, 9) > ltp_limit;
    }

    if prev_voiced {
        *voiced_hangover = 0;
    } else {
        let temp = voiced_hangover.wrapping_add(1);
        *voiced_hangover = if temp > 10 { 10 } else { temp };
    }
    inbg_noise
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_c;
