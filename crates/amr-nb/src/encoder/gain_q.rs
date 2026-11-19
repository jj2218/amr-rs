//! Gain quantisation: `gainQuant` and everything below it (`calc_filt_energies`,
//! `calc_unfilt_energies`, `calc_target_energy`, `G_code`, `q_gain_code`, `Qua_gain`,
//! `MR795_gain_quant`, `gain_adapt`, `MR475_gain_quant`, `MR475_update_unq_pred`).
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `enc/src/{gain_q,calc_en,g_code,
//! q_gain_c,qua_gain,qgain795,g_adapt,qgain475}.cpp`.
#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1

use super::layout::{Soa, cols, first_min, soa};
use super::pitch_cl::q_gain_pitch;
use crate::basic_op::{
    Flag, MAX_32, abs_s, add, div_s, l_add, l_extract, l_mac, l_mac_dot, l_mac_energy, l_mult,
    l_shl, l_shr, l_sub, mac_32, mac_32_16, mpy_32_16, mult, negate, norm_l, round, shl, shr,
    shr_r, sub,
};
use crate::cnst::{CodecMode, L_SUBFR};
use crate::gains::{GcPredState, gc_pred, gc_pred_update, gmed_n};
use crate::math::{log2, pow2, sqrt_l_exp};
pub(crate) use crate::simd::{mr475_search, mr795_search, qua_gain_exact};
use crate::tables::{QUA_GAIN_CODE, TABLE_GAIN_HIGHRATES, TABLE_GAIN_LOWRATES, TABLE_GAIN_MR475};

const NB_QUA_CODE: usize = 32;
const LTPG_MEM_SIZE: usize = 5;
const MIN_QUA_ENER: i16 = -5443;
const MIN_QUA_ENER_MR122: i16 = -32768;
const MAX_QUA_ENER: i16 = 3037;
const MAX_QUA_ENER_MR122: i16 = 18284;

/// `GainAdaptState` (`repr(C)`: the C diff passes it by pointer).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
#[repr(C)]
pub(crate) struct GainAdaptState {
    pub onset: i16,
    pub prev_alpha: i16,
    pub prev_gc: i16,
    pub ltpg_mem: [i16; LTPG_MEM_SIZE],
}
// The C struct: 3 + LTPG_MEM_SIZE Word16s.
const _: () = assert!(size_of::<GainAdaptState>() == 2 * (3 + 5) && LTPG_MEM_SIZE == 5);

/// `gainQuantState`. `gain_idx_ptr` (an index into `ana`) is kept by the caller.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct GainQuantState {
    pub sf0_exp_gcode0: i16,
    pub sf0_frac_gcode0: i16,
    pub sf0_exp_target_en: i16,
    pub sf0_frac_target_en: i16,
    pub sf0_exp_coeff: [i16; 5],
    pub sf0_frac_coeff: [i16; 5],
    pub gc_pred_st: GcPredState,
    pub gc_pred_unq_st: GcPredState,
    pub adapt_st: GainAdaptState,
    /// Position in `ana` of the MR475 joint gain index (written in the odd subframe).
    pub gain_idx_pos: usize,
}

/// `calc_unfilt_energies`: energies of the LP residual, excitation, their correlation and the
/// LTP residual; returns the LTP prediction gain `ltpg`.
pub(super) fn calc_unfilt_energies(
    res: &[i16],
    exc: &[i16],
    code: &[i16],
    gain_pit: i16,
    frac_en: &mut [i16; 4],
    exp_en: &mut [i16; 4],
    ov: &mut Flag,
) -> i16 {
    let (res, exc, code) = (&res[..L_SUBFR], &exc[..L_SUBFR], &code[..L_SUBFR]);
    // The C's saturating L_mac chains (opencore wraps: wrong for full-scale res/exc).
    let s1 = l_mac_energy(res, ov);
    let s2 = l_mac_energy(exc, ov);
    let s3 = l_mac_dot(0, exc, code, ov);
    // The LTP residual t = sub(res, round(L_shl(L_mult(exc, gain_pit), 1))) per sample as
    // selects (vectorises; L_shl never flags, the others do).
    let mut sat = false;
    let t: [i16; L_SUBFR] = core::array::from_fn(|i| {
        let p = exc[i] as i32 * gain_pit as i32;
        let l = if p == 0x4000_0000 { MAX_32 } else { p << 1 };
        let x = if (l << 1) >> 1 != l {
            (l >> 31) ^ MAX_32
        } else {
            l << 1
        };
        let r = if x > 0x7fff_7fff {
            32767
        } else {
            x.wrapping_add(0x8000) >> 16
        };
        let d = res[i] as i32 - r;
        sat |= (p == 0x4000_0000) | (x > 0x7fff_7fff) | !(-32768..=32767).contains(&d);
        d.clamp(-32768, 32767) as i16
    });
    *ov |= sat;
    let s4 = l_mac_energy(&t, ov);
    if s1 < 400 {
        frac_en[0] = 0;
        exp_en[0] = -15;
    } else {
        let exp = norm_l(s1);
        frac_en[0] = (l_shl(s1, exp, ov) >> 16) as i16;
        exp_en[0] = 15 - exp;
    }
    let exp = norm_l(s2);
    frac_en[1] = (l_shl(s2, exp, ov) >> 16) as i16;
    exp_en[1] = 15 - exp;
    let exp = norm_l(s3);
    frac_en[2] = (l_shl(s3, exp, ov) >> 16) as i16;
    exp_en[2] = 2 - exp;
    let exp = norm_l(s4);
    let ltp_res_en = (l_shl(s4, exp, ov) >> 16) as i16;
    let mut exp = 15 - exp;
    frac_en[3] = ltp_res_en;
    exp_en[3] = exp;
    if ltp_res_en > 0 && frac_en[0] != 0 {
        let pred_gain = div_s(shr(frac_en[0], 1, ov), ltp_res_en);
        exp = sub(exp, exp_en[0], ov);
        let mut l_temp = (pred_gain as i32) << 16;
        l_temp = l_shr(l_temp, exp.wrapping_add(3), ov);
        let (ltpg_exp, ltpg_frac) = log2(l_temp, ov);
        let l_temp = crate::basic_op::l_comp(ltpg_exp.wrapping_sub(27), ltpg_frac, ov);
        round(l_shl(l_temp, 13, ov), ov)
    } else {
        0
    }
}

/// `calc_filt_energies`: the coefficients of the gain-quantiser error function; for MR475 and
/// MR795 also the optimum unquantised codebook gain `(frac, exp)`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn calc_filt_energies(
    mode: CodecMode,
    xn: &[i16],
    xn2: &[i16],
    y1: &[i16],
    y2: &[i16],
    g_coeff: &[i16],
    frac_coeff: &mut [i16; 5],
    exp_coeff: &mut [i16; 5],
    cod_gain_frac: &mut i16,
    cod_gain_exp: &mut i16,
    ov: &mut Flag,
) {
    frac_coeff[0] = g_coeff[0];
    exp_coeff[0] = g_coeff[1];
    frac_coeff[1] = negate(g_coeff[2]);
    exp_coeff[1] = add(g_coeff[3], 1, ov);
    let special = matches!(mode, CodecMode::MR795 | CodecMode::MR475);
    let init = if special { 0 } else { 1 };
    let (xn, xn2, y1) = (&xn[..L_SUBFR], &xn2[..L_SUBFR], &y1[..L_SUBFR]);
    let scaled_y2: [i16; L_SUBFR] = core::array::from_fn(|i| y2[i] >> 3);
    // The C's three L_mac chains. |tmp| <= 4096: s1's terms are >= 0 and sum below 2^31, so
    // it never saturates (a plain sum). s2, s3 run as the chains: a plain-sum fast path under
    // a 2 max|x| sum|tmp| < MAX_32 bound measured no faster in any mode.
    let s1 = scaled_y2.iter().fold(init, |s: i32, &t| {
        s.wrapping_add((t as i32 * t as i32) << 1)
    });
    let mut chain = |v: &[i16]| {
        v.iter()
            .zip(&scaled_y2)
            .fold(init, |s, (&x, &t)| l_mac(s, x, t, ov))
    };
    let (s2, s3) = (chain(xn), chain(y1));
    let exp = norm_l(s1);
    frac_coeff[2] = (l_shl(s1, exp, ov) >> 16) as i16;
    exp_coeff[2] = -3 - exp;
    let exp = norm_l(s2);
    frac_coeff[3] = negate((l_shl(s2, exp, ov) >> 16) as i16);
    exp_coeff[3] = 7 - exp;
    let exp = norm_l(s3);
    frac_coeff[4] = (l_shl(s3, exp, ov) >> 16) as i16;
    exp_coeff[4] = 7 - exp;
    if special {
        // Saturating chain like TS 26.073 (opencore wraps; xn2 can be full scale).
        let s1 = xn2
            .iter()
            .zip(&scaled_y2)
            .fold(0, |s, (&x, &t)| l_mac(s, x, t, ov));
        let exp = norm_l(s1);
        let frac = (l_shl(s1, exp, ov) >> 16) as i16;
        let exp = 6 - exp;
        if frac <= 0 {
            *cod_gain_frac = 0;
            *cod_gain_exp = 0;
        } else {
            *cod_gain_frac = div_s(shr(frac, 1, ov), frac_coeff[2]);
            *cod_gain_exp = exp.wrapping_sub(exp_coeff[2]).wrapping_sub(14);
        }
    }
}

/// `calc_target_energy`: energy of the target `xn` as `(exp, frac)`. TS 26.073's `L_mac` chain
/// over x² (terms >= 0, so it is min(2 sum, MAX_32)). opencore drops the doubling, which only
/// shifts the exponent (MR475 uses exponent differences) but saturates at 2^31 instead of
/// 2^30 and wraps beyond: wrong on full-scale input.
pub(super) fn calc_target_energy(xn: &[i16], ov: &mut Flag) -> (i16, i16) {
    let sum: i64 = xn[..L_SUBFR].iter().map(|&x| 2 * x as i64 * x as i64).sum();
    if sum > MAX_32 as i64 {
        *ov = true;
    }
    let s = sum.min(MAX_32 as i64) as i32;
    let exp = norm_l(s);
    (16 - exp, (l_shl(s, exp, ov) >> 16) as i16)
}

/// `G_code`: innovative codebook gain (MR122). TS 26.073 `g_code.c` with its saturating basic
/// operators (opencore's rewrite wraps: wrong gains on full-scale input). Once per subframe.
pub(super) fn g_code(xn2: &[i16], y2: &[i16]) -> i16 {
    let ov = &mut false;
    let scal: [i16; L_SUBFR] = core::array::from_fn(|i| y2[i] >> 1);
    let s = l_mac_dot(1, &xn2[..L_SUBFR], &scal, ov); // the C's L_mac chain, exact
    let exp_xy = norm_l(s);
    let xy = (l_shl(s, exp_xy, ov) >> 16) as i16;
    if xy <= 0 {
        return 0;
    }
    let s = l_mac_energy(&scal, ov);
    let exp_yy = norm_l(s);
    let yy = (l_shl(s, exp_yy, ov) >> 16) as i16;
    let gain = div_s(xy >> 1, yy);
    let i = exp_xy + 5 - exp_yy;
    shl(shr(gain, i, ov), 1, ov)
}

/// `q_gain_code`: scalar quantisation of the codebook gain (MR122, MR795 fallback).
pub(super) fn q_gain_code(
    mode: CodecMode,
    exp_gcode0: i16,
    frac_gcode0: i16,
    gain: &mut i16,
    qua_ener_mr122: &mut i16,
    qua_ener: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let mr122 = mode == CodecMode::MR122;
    let g_q0 = if mr122 { *gain >> 1 } else { *gain };
    let mut gcode0 = pow2(exp_gcode0, frac_gcode0, ov) as i16;
    gcode0 = shl(gcode0, if mr122 { 4 } else { 5 }, ov);
    // Saturating like TS 26.073 `q_gain_c.c` (opencore wraps).
    let o = &mut false;
    let mut err_of = |q: i16| -> i16 { abs_s(sub(g_q0, mult(gcode0, q, o), o)) };
    let mut err_min = err_of(QUA_GAIN_CODE[0]);
    let mut index = 0usize;
    for i in 1..NB_QUA_CODE {
        let err = err_of(QUA_GAIN_CODE[3 * i]);
        if err < err_min {
            err_min = err;
            index = i;
        }
    }
    let p = &QUA_GAIN_CODE[3 * index..];
    let temp = mult(gcode0, p[0], ov);
    *gain = if mr122 { shl(temp, 1, ov) } else { temp };
    *qua_ener_mr122 = p[1];
    *qua_ener = p[2];
    index as i16
}

/// Normalises 5 (or 10) error-function coefficients to a common exponent `e_max` as
/// double-precision `(hi, lo)` pairs.
pub(super) fn normalise_coeffs(
    frac: &[i16],
    exp_max: &[i16],
    e_max: i16,
    coeff: &mut [i16],
    coeff_lo: &mut [i16],
    ov: &mut Flag,
) {
    for i in 0..frac.len() {
        let j = e_max.wrapping_sub(exp_max[i]);
        let l_tmp = l_shr((frac[i] as i32) << 16, j, ov);
        (coeff[i], coeff_lo[i]) = l_extract(l_tmp, ov);
    }
}

/// `Qua_gain`: joint pitch/codebook gain VQ (MR59, MR67, MR74, MR102 and MR515).
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn qua_gain(
    mode: CodecMode,
    exp_gcode0: i16,
    frac_gcode0: i16,
    frac_coeff: &[i16; 5],
    exp_coeff: &[i16; 5],
    gp_limit: i16,
    gain_pit: &mut i16,
    gain_cod: &mut i16,
    qua_ener_mr122: &mut i16,
    qua_ener: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let (table, cols): (&[i16], [&[i16]; 2]) =
        if matches!(mode, CodecMode::MR102 | CodecMode::MR74 | CodecMode::MR67) {
            (&TABLE_GAIN_HIGHRATES, cols(&HIGHRATES_SOA))
        } else {
            (&TABLE_GAIN_LOWRATES, cols(&LOWRATES_SOA))
        };
    let gcode0 = pow2(14, frac_gcode0, ov) as i16;
    let exp_code = exp_gcode0.wrapping_sub(11);
    let mut exp_max = [0i16; 5];
    exp_max[0] = exp_coeff[0].wrapping_sub(13);
    exp_max[1] = exp_coeff[1].wrapping_sub(14);
    let temp = shl(exp_code, 1, ov).wrapping_add(15);
    exp_max[2] = add(exp_coeff[2], temp, ov);
    exp_max[3] = add(exp_coeff[3], exp_code, ov);
    exp_max[4] = add(exp_coeff[4], exp_code.wrapping_add(1), ov);
    let e_max = exp_max.iter().copied().max().unwrap().wrapping_add(1);
    let mut coeff = [0i16; 5];
    let mut coeff_lo = [0i16; 5];
    normalise_coeffs(frac_coeff, &exp_max, e_max, &mut coeff, &mut coeff_lo, ov);
    let index = match qua_gain_exact(&cols, gcode0, &coeff, &coeff_lo, gp_limit) {
        Some(i) => i,
        None => qua_gain_search(table, gcode0, &coeff, &coeff_lo, gp_limit, ov),
    };
    let p = &table[shl(index as i16, 2, ov) as usize..];
    *gain_pit = p[0];
    let g_code = p[1];
    *qua_ener_mr122 = p[2];
    *qua_ener = p[3];
    let l_tmp = l_mult(g_code, gcode0, ov);
    let l_tmp = l_shr(l_tmp, 10i16.wrapping_sub(exp_gcode0), ov);
    *gain_cod = (l_tmp >> 16) as i16;
    index as i16
}

/// `TABLE_GAIN_HIGHRATES` / `_LOWRATES` columns g_pitch, g_code (component-major: rows in lanes).
static HIGHRATES_SOA: Soa<2, 128> = soa(&TABLE_GAIN_HIGHRATES, 4);
static LOWRATES_SOA: Soa<2, 64> = soa(&TABLE_GAIN_LOWRATES, 4);

/// The `Qua_gain` search without saturation, `W` table rows per step in i32 lanes (wrapping,
/// branchless). `None` if, in a row the C evaluates (`g_pitch <= gp_limit`), any `mult`,
/// `Mpy_32_16` or `L_add` would saturate (the C then also sets the overflow flag): the caller
/// runs the faithful search. Otherwise no C operation saturates, so the result is the C's.
#[inline(always)]
pub(crate) fn qua_gain_exact_portable<const W: usize>(
    cols: &[&[i16]; 2],
    gcode0: i16,
    coeff: &[i16; 5],
    coeff_lo: &[i16; 5],
    gp_limit: i16,
) -> Option<usize> {
    // First tier: from the columns' maxima, every |n| (g_pitch, g_code = g_c g0 >> 15 and the
    // products) and hence every term and partial sum is bounded; within MAX_32 (and n <= 0x7fff,
    // so no `mult` saturates and no `Mpy_32_16` product is 2^30) the lanes skip their flags.
    let amax = |c: &[i16]| c.iter().fold(0i16, |m, &v| m.max(v.saturating_abs())) as i64;
    let (gp, gc) = (amax(cols[0]), (amax(cols[1]) * (gcode0 as i64).abs()) >> 15);
    let n = [(gp * gp) >> 15, gp, (gc * gc) >> 15, gc, (gc * gp) >> 15];
    let bound: i64 = (0..5)
        .map(|k| {
            2 * (coeff[k] as i64).abs() * n[k] + 2 * (((coeff_lo[k] as i64).abs() * n[k]) >> 15)
        })
        .sum();
    if bound <= MAX_32 as i64 && n.iter().all(|&v| v <= 0x7fff) {
        qua_gain_lanes::<W, false>(cols, gcode0, coeff, coeff_lo, gp_limit)
    } else {
        qua_gain_lanes::<W, true>(cols, gcode0, coeff, coeff_lo, gp_limit)
    }
}

/// `qua_gain_exact_portable`'s lanes; `CHECK`: flag saturating operations (else proven absent).
#[inline(always)]
fn qua_gain_lanes<const W: usize, const CHECK: bool>(
    cols: &[&[i16]; 2],
    gcode0: i16,
    coeff: &[i16; 5],
    coeff_lo: &[i16; 5],
    gp_limit: i16,
) -> Option<usize> {
    // Mpy_32_16(hi, lo, n) = 2*hi*n + 2*((lo*n) >> 15); flags hi*n == 2^30 and overflow.
    #[inline(always)]
    fn mpy<const CHECK: bool>(hi: i16, lo: i16, n: i32, bad: &mut bool) -> i32 {
        let p = hi as i32 * n;
        let lp = p << 1;
        let r2 = ((lo as i32 * n) >> 15) << 1;
        let t = lp.wrapping_add(r2);
        if CHECK {
            *bad |= (p == 0x4000_0000) | (((lp ^ t) & (r2 ^ t)) < 0);
        }
        t
    }
    #[inline(always)]
    fn add<const CHECK: bool>(a: i32, b: i32, bad: &mut bool) -> i32 {
        let t = a.wrapping_add(b);
        if CHECK {
            *bad |= ((a ^ t) & (b ^ t)) < 0;
        }
        t
    }
    let g0 = gcode0 as i32;
    let n = cols[0].len();
    assert!(n.is_multiple_of(W) && cols[1].len() == n);
    let mut dist = [MAX_32; W];
    let mut idx = [0i32; W];
    let mut bad_any = false;
    for b in (0..n).step_by(W) {
        let gp: &[i16; W] = cols[0][b..b + W].try_into().unwrap();
        let gc: &[i16; W] = cols[1][b..b + W].try_into().unwrap();
        for l in 0..W {
            let g_pitch = gp[l] as i32;
            let g_code = (gc[l] as i32 * g0) >> 15;
            let g2_pitch = (g_pitch * g_pitch) >> 15;
            let g2_code = (g_code * g_code) >> 15;
            let g_pit_cod = (g_code * g_pitch) >> 15;
            let mut bad = CHECK && (g_code | g2_pitch | g2_code | g_pit_cod) > 0x7fff;
            let mut d = mpy::<CHECK>(coeff[0], coeff_lo[0], g2_pitch, &mut bad);
            d = add::<CHECK>(
                d,
                mpy::<CHECK>(coeff[1], coeff_lo[1], g_pitch, &mut bad),
                &mut bad,
            );
            d = add::<CHECK>(
                d,
                mpy::<CHECK>(coeff[2], coeff_lo[2], g2_code, &mut bad),
                &mut bad,
            );
            d = add::<CHECK>(
                d,
                mpy::<CHECK>(coeff[3], coeff_lo[3], g_code, &mut bad),
                &mut bad,
            );
            d = add::<CHECK>(
                d,
                mpy::<CHECK>(coeff[4], coeff_lo[4], g_pit_cod, &mut bad),
                &mut bad,
            );
            let valid = g_pitch <= gp_limit as i32;
            bad_any |= valid & bad;
            let better = valid & (d < dist[l]);
            dist[l] = if better { d } else { dist[l] };
            idx[l] = if better { (b + l) as i32 } else { idx[l] };
        }
    }
    if bad_any {
        return None;
    }
    // Lane l saw rows l, l + W, ... in order, keeping the first minimum.
    Some(first_min(&dist, &idx) as usize)
}

/// The faithful `Qua_gain` search (saturating operators).
fn qua_gain_search(
    table: &[i16],
    gcode0: i16,
    coeff: &[i16; 5],
    coeff_lo: &[i16; 5],
    gp_limit: i16,
    ov: &mut Flag,
) -> usize {
    let mut dist_min = MAX_32;
    let mut index = 0usize;
    for (i, row) in table.as_chunks::<4>().0.iter().enumerate() {
        let g_pitch = row[0];
        if g_pitch <= gp_limit {
            let g_code = mult(row[1], gcode0, ov);
            let g2_pitch = mult(g_pitch, g_pitch, ov);
            let g2_code = mult(g_code, g_code, ov);
            let g_pit_cod = mult(g_code, g_pitch, ov);
            let mut l_tmp = mpy_32_16(coeff[0], coeff_lo[0], g2_pitch, ov);
            let t = mpy_32_16(coeff[1], coeff_lo[1], g_pitch, ov);
            l_tmp = l_add(l_tmp, t, ov);
            let t = mpy_32_16(coeff[2], coeff_lo[2], g2_code, ov);
            l_tmp = l_add(l_tmp, t, ov);
            let t = mpy_32_16(coeff[3], coeff_lo[3], g_code, ov);
            l_tmp = l_add(l_tmp, t, ov);
            let t = mpy_32_16(coeff[4], coeff_lo[4], g_pit_cod, ov);
            l_tmp = l_add(l_tmp, t, ov);
            if l_tmp < dist_min {
                dist_min = l_tmp;
                index = i;
            }
        }
    }
    index
}

/// `gain_adapt`: adaptive balance factor `alpha` for the MR795 codebook-gain search.
pub(super) fn gain_adapt(st: &mut GainAdaptState, ltpg: i16, gain_cod: i16, ov: &mut Flag) -> i16 {
    let mut adapt = if ltpg <= 2721 {
        0
    } else if ltpg <= 5443 {
        1
    } else {
        2
    };
    let tmp = shr_r(gain_cod, 1, ov);
    if tmp > st.prev_gc && gain_cod > 200 {
        st.onset = 8;
    } else if st.onset != 0 {
        st.onset -= 1;
    }
    if st.onset != 0 && adapt < 2 {
        adapt += 1;
    }
    st.ltpg_mem[0] = ltpg;
    let filt = gmed_n(&st.ltpg_mem, 5);
    let mut result = if adapt == 0 {
        if filt > 5443 {
            0
        } else if filt < 0 {
            16384
        } else {
            let f = shl(filt, 2, ov);
            16384 - mult(24660, f, ov)
        }
    } else {
        0
    };
    if st.prev_alpha == 0 {
        result = shr(result, 1, ov);
    }
    st.prev_alpha = result;
    st.prev_gc = gain_cod;
    st.ltpg_mem.copy_within(0..LTPG_MEM_SIZE - 1, 1);
    result
}

/// `MR795_gain_code_quant3`: joint search over 3 pitch-gain candidates × 32 codebook gains.
/// Returns `(gain_pit, gain_pit_ind, gain_cod, gain_cod_ind)`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn mr795_gain_code_quant3(
    exp_gcode0: i16,
    gcode0: i16,
    g_pitch_cand: &[i16; 3],
    g_pitch_cind: &[i16; 3],
    frac_coeff: &[i16; 5],
    exp_coeff: &[i16; 5],
    qua_ener_mr122: &mut i16,
    qua_ener: &mut i16,
    ov: &mut Flag,
) -> (i16, i16, i16, i16) {
    let exp_code = exp_gcode0.wrapping_sub(10);
    let mut exp_max = [0i16; 5];
    exp_max[0] = exp_coeff[0].wrapping_sub(13);
    exp_max[1] = exp_coeff[1].wrapping_sub(14);
    exp_max[2] = exp_coeff[2]
        .wrapping_add(shl(exp_code, 1, ov))
        .wrapping_add(15);
    exp_max[3] = exp_coeff[3].wrapping_add(exp_code);
    exp_max[4] = exp_coeff[4].wrapping_add(exp_code.wrapping_add(1));
    let e_max = add(exp_max.iter().copied().max().unwrap(), 1, ov);
    let mut coeff = [0i16; 5];
    let mut coeff_lo = [0i16; 5];
    normalise_coeffs(frac_coeff, &exp_max, e_max, &mut coeff, &mut coeff_lo, ov);
    let mut l_tmp0 = [0i32; 3];
    for (j, &g_pitch) in g_pitch_cand.iter().enumerate() {
        let g2_pitch = mult(g_pitch, g_pitch, ov);
        let l = mpy_32_16(coeff[0], coeff_lo[0], g2_pitch, ov);
        l_tmp0[j] = mac_32_16(l, coeff[1], coeff_lo[1], g_pitch, ov);
    }
    let k = match mr795_search(&coeff, &coeff_lo, &l_tmp0, gcode0, g_pitch_cand) {
        Some(k) => k,
        None => mr795_search_c(&coeff, &coeff_lo, &l_tmp0, gcode0, g_pitch_cand, ov),
    };
    let (pit_ind, cod_ind) = (k / NB_QUA_CODE, k % NB_QUA_CODE);
    let p = &QUA_GAIN_CODE[3 * cod_ind..];
    *qua_ener_mr122 = p[1];
    *qua_ener = p[2];
    let l_tmp = l_mult(p[0], gcode0, ov);
    let l_tmp = l_shr(l_tmp, 9i16.wrapping_sub(exp_gcode0), ov);
    (
        (g_pitch_cand[pit_ind]),
        g_pitch_cind[pit_ind],
        (l_tmp >> 16) as i16,
        cod_ind as i16,
    )
}

/// `QUA_GAIN_CODE`'s gain column (component-major: rows in lanes).
static QUA_CODE_SOA: Soa<1, NB_QUA_CODE> = soa(&QUA_GAIN_CODE, 3);

/// `MR795_gain_code_quant3`'s search after the per-candidate `l_tmp0`: the 3 x 32 distances
/// (candidate `j` outer, code `i` inner, as the C) in `W`-row i32 lanes, wrapping; returns
/// `32 * j + i` of the first minimum. `None` if any `L_mult` or `L_mac` in it would saturate:
/// the caller runs the C's loop. The `mult`s of a `lo` half never saturate (`lo` in
/// [0, 32767]), nor `mult(code, gcode0)` for `gcode0 >= 0`, but that one is checked anyway.
#[inline(always)]
pub(crate) fn mr795_search_portable<const W: usize>(
    coeff: &[i16; 5],
    coeff_lo: &[i16; 5],
    l_tmp0: &[i32; 3],
    gcode0: i16,
    g_pitch: &[i16; 3],
) -> Option<usize> {
    /// `L_mac(acc, a, b)`, flagging `a * b == 2^30` and overflow.
    #[inline(always)]
    fn mac(acc: i32, a: i32, b: i32, bad: &mut bool) -> i32 {
        let p = a * b;
        let p2 = p << 1;
        let t = acc.wrapping_add(p2);
        *bad |= (p == 0x4000_0000) | (((acc ^ t) & (p2 ^ t)) < 0);
        t
    }
    /// `L_mac(acc, x, 1)` (`x` a `mult` of a `lo` half: never 2^30).
    #[inline(always)]
    fn mac1(acc: i32, x: i32, bad: &mut bool) -> i32 {
        let p2 = x << 1;
        let t = acc.wrapping_add(p2);
        *bad |= ((acc ^ t) & (p2 ^ t)) < 0;
        t
    }
    /// `L_Extract(L_mult(a, b))` (`a * b == 2^30` flagged).
    #[inline(always)]
    fn extract(a: i32, b: i32, bad: &mut bool) -> (i32, i32) {
        let p = a * b;
        *bad |= p == 0x4000_0000;
        let l = p << 1;
        let hi = l >> 16;
        (hi, (l >> 1) - (hi << 15))
    }
    debug_assert!(coeff_lo.iter().all(|&v| v >= 0));
    let c = coeff.map(|v| v as i32);
    let cl = coeff_lo.map(|v| v as i32);
    let g0 = gcode0 as i32;
    let mut best = [MAX_32; W];
    let mut bidx = [0i32; W];
    let mut bad_any = false;
    for j in 0..3 {
        let (gp, acc0) = (g_pitch[j] as i32, l_tmp0[j]);
        for b in (0..NB_QUA_CODE).step_by(W) {
            let gc: &[i16; W] = QUA_CODE_SOA[0][b..b + W].try_into().unwrap();
            for l in 0..W {
                let g = (gc[l] as i32 * g0) >> 15;
                let mut bad = g > 0x7fff;
                let (h2, l2) = extract(g, g, &mut bad);
                let (hp, lp) = extract(g, gp, &mut bad);
                let mut t = mac(acc0, c[2], h2, &mut bad);
                t = mac1(t, (c[2] * l2) >> 15, &mut bad);
                t = mac1(t, (cl[2] * h2) >> 15, &mut bad);
                t = mac(t, c[3], g, &mut bad);
                t = mac1(t, (cl[3] * g) >> 15, &mut bad);
                t = mac(t, c[4], hp, &mut bad);
                t = mac1(t, (c[4] * lp) >> 15, &mut bad);
                t = mac1(t, (cl[4] * hp) >> 15, &mut bad);
                bad_any |= bad;
                let lt = t < best[l];
                best[l] = if lt { t } else { best[l] };
                bidx[l] = if lt {
                    (NB_QUA_CODE * j + b + l) as i32
                } else {
                    bidx[l]
                };
            }
        }
    }
    if bad_any {
        return None;
    }
    Some(first_min(&best, &bidx) as usize)
}

/// The C's `MR795_gain_code_quant3` search loop (saturating operators); returns `32 * j + i`.
pub(crate) fn mr795_search_c(
    coeff: &[i16; 5],
    coeff_lo: &[i16; 5],
    l_tmp0: &[i32; 3],
    gcode0: i16,
    g_pitch_cand: &[i16; 3],
    ov: &mut Flag,
) -> usize {
    let mut dist_min = MAX_32;
    let mut k = 0usize;
    for j in 0..3 {
        let g_pitch = g_pitch_cand[j];
        for i in 0..NB_QUA_CODE {
            let g_code = mult(QUA_GAIN_CODE[3 * i], gcode0, ov);
            let (g2_code_h, g2_code_l) = l_extract(l_mult(g_code, g_code, ov), ov);
            let (g_pit_cod_h, g_pit_cod_l) = l_extract(l_mult(g_code, g_pitch, ov), ov);
            let mut l_tmp = mac_32(l_tmp0[j], coeff[2], coeff_lo[2], g2_code_h, g2_code_l, ov);
            l_tmp = mac_32_16(l_tmp, coeff[3], coeff_lo[3], g_code, ov);
            l_tmp = mac_32(l_tmp, coeff[4], coeff_lo[4], g_pit_cod_h, g_pit_cod_l, ov);
            if l_tmp < dist_min {
                dist_min = l_tmp;
                k = NB_QUA_CODE * j + i;
            }
        }
    }
    k
}

/// `MR795_gain_code_quant_mod`: codebook-gain requantisation with the energy-matching criterion.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn mr795_gain_code_quant_mod(
    gain_pit: i16,
    exp_gcode0: i16,
    gcode0: i16,
    frac_en: &[i16; 4],
    exp_en: &[i16; 4],
    alpha: i16,
    gain_cod_unq: i16,
    gain_cod: &mut i16,
    qua_ener_mr122: &mut i16,
    qua_ener: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let mut coeff = [0i16; 5];
    let mut coeff_lo = [0i16; 5];
    let mut exp_coeff = [0i16; 5];
    let gain_code = shl(*gain_cod, 10i16.wrapping_sub(exp_gcode0), ov);
    let g2_pitch = mult(gain_pit, gain_pit, ov);
    let one_alpha = add(32767i16.wrapping_sub(alpha), 1, ov);
    let mut l_t1 = l_mult(alpha, frac_en[1], ov);
    l_t1 = l_shl(l_t1, 1, ov);
    let tmp = (l_t1 >> 16) as i16;
    l_t1 = l_mult(tmp, g2_pitch, ov);
    exp_coeff[1] = exp_en[1].wrapping_sub(15);
    let t = l_mult(alpha, frac_en[2], ov);
    let tmp = (l_shl(t, 1, ov) >> 16) as i16;
    coeff[2] = mult(tmp, gain_pit, ov);
    let exp = exp_gcode0.wrapping_sub(10);
    exp_coeff[2] = add(exp_en[2], exp, ov);
    let t = l_mult(alpha, frac_en[3], ov);
    coeff[3] = (l_shl(t, 1, ov) >> 16) as i16;
    let exp = shl(exp_gcode0, 1, ov).wrapping_sub(7);
    exp_coeff[3] = add(exp_en[3], exp, ov);
    coeff[4] = mult(one_alpha, frac_en[3], ov);
    exp_coeff[4] = add(exp_coeff[3], 1, ov);
    let l_tmp = l_mult(alpha, frac_en[0], ov);
    let (mut l_t0, exp) = sqrt_l_exp(l_tmp, ov);
    let exp = exp.wrapping_add(47);
    exp_coeff[0] = exp_en[0].wrapping_sub(exp);
    let mut e_max = exp_coeff[0].wrapping_add(31);
    for i in 1..=4 {
        if exp_coeff[i] > e_max {
            e_max = exp_coeff[i];
        }
    }
    l_t1 = l_shr(l_t1, e_max.wrapping_sub(exp_coeff[1]), ov);
    for i in 2..=4 {
        let l_tmp = l_shr(
            (coeff[i] as i32) << 16,
            e_max.wrapping_sub(exp_coeff[i]),
            ov,
        );
        (coeff[i], coeff_lo[i]) = l_extract(l_tmp, ov);
    }
    let exp = e_max.wrapping_sub(31);
    let tmp = exp.wrapping_sub(exp_coeff[0]);
    let s = shr(tmp, 1, ov);
    l_t0 = l_shr(l_t0, s, ov);
    if (tmp & 1) != 0 {
        (coeff[0], coeff_lo[0]) = l_extract(l_t0, ov);
        l_t0 = mpy_32_16(coeff[0], coeff_lo[0], 23170, ov);
    }
    let mut dist_min = MAX_32;
    let mut index = 0usize;
    for i in 0..NB_QUA_CODE {
        let g_code = mult(QUA_GAIN_CODE[3 * i], gcode0, ov);
        if g_code >= gain_code {
            break;
        }
        let (g2_code_h, g2_code_l) = l_extract(l_mult(g_code, g_code, ov), ov);
        let tmp = sub(g_code, gain_cod_unq, ov);
        let (d2_code_h, d2_code_l) = l_extract(l_mult(tmp, tmp, ov), ov);
        let mut l_tmp = mac_32_16(l_t1, coeff[2], coeff_lo[2], g_code, ov);
        l_tmp = mac_32(l_tmp, coeff[3], coeff_lo[3], g2_code_h, g2_code_l, ov);
        let (sq, exp) = sqrt_l_exp(l_tmp, ov);
        let s = shr(exp, 1, ov);
        l_tmp = l_shr(sq, s, ov);
        let d = l_sub(l_tmp, l_t0, ov);
        let tmp = round(d, ov);
        l_tmp = l_mult(tmp, tmp, ov);
        l_tmp = mac_32(l_tmp, coeff[4], coeff_lo[4], d2_code_h, d2_code_l, ov);
        if l_tmp < dist_min {
            dist_min = l_tmp;
            index = i;
        }
    }
    let p = &QUA_GAIN_CODE[3 * index..];
    *qua_ener_mr122 = p[1];
    *qua_ener = p[2];
    let l_tmp = l_mult(p[0], gcode0, ov);
    let l_tmp = l_shr(l_tmp, 9i16.wrapping_sub(exp_gcode0), ov);
    *gain_cod = (l_tmp >> 16) as i16;
    index as i16
}

/// `MR795_gain_quant`: pitch gain scalar + codebook gain (3 candidates), then the adaptive
/// requantisation of the codebook gain. Appends 2 indices to `ana`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn mr795_gain_quant(
    adapt_st: &mut GainAdaptState,
    res: &[i16],
    exc: &[i16],
    code: &[i16],
    frac_coeff: &[i16; 5],
    exp_coeff: &[i16; 5],
    exp_code_en: i16,
    frac_code_en: i16,
    exp_gcode0: i16,
    frac_gcode0: i16,
    cod_gain_frac: i16,
    cod_gain_exp: i16,
    gp_limit: i16,
    gain_pit: &mut i16,
    gain_cod: &mut i16,
    qua_ener_mr122: &mut i16,
    qua_ener: &mut i16,
    ana: &mut [i16],
    ana_pos: &mut usize,
    ov: &mut Flag,
) {
    let mut frac_en = [0i16; 4];
    let mut exp_en = [0i16; 4];
    let mut cand = [0i16; 3];
    let mut cind = [0i16; 3];
    q_gain_pitch(
        CodecMode::MR795,
        gp_limit,
        gain_pit,
        Some((&mut cand, &mut cind)),
        ov,
    );
    let gcode0 = pow2(14, frac_gcode0, ov) as i16;
    let (gp, gain_pit_index, gc, mut gain_cod_index) = mr795_gain_code_quant3(
        exp_gcode0,
        gcode0,
        &cand,
        &cind,
        frac_coeff,
        exp_coeff,
        qua_ener_mr122,
        qua_ener,
        ov,
    );
    *gain_pit = gp;
    *gain_cod = gc;
    let ltpg = calc_unfilt_energies(res, exc, code, *gain_pit, &mut frac_en, &mut exp_en, ov);
    let alpha = gain_adapt(adapt_st, ltpg, *gain_cod, ov);
    if frac_en[0] != 0 && alpha > 0 {
        frac_en[3] = frac_code_en;
        exp_en[3] = exp_code_en;
        let exp = sub(cod_gain_exp, exp_gcode0, ov).wrapping_add(10);
        let gain_cod_unq = shl(cod_gain_frac, exp, ov);
        gain_cod_index = mr795_gain_code_quant_mod(
            *gain_pit,
            exp_gcode0,
            gcode0,
            &frac_en,
            &exp_en,
            alpha,
            gain_cod_unq,
            gain_cod,
            qua_ener_mr122,
            qua_ener,
            ov,
        );
    }
    ana[*ana_pos] = gain_pit_index;
    ana[*ana_pos + 1] = gain_cod_index;
    *ana_pos += 2;
}

/// `MR475_quant_store_results`: decodes the chosen MR475 gains and updates the predictor.
pub(super) fn mr475_quant_store_results(
    pred_st: &mut GcPredState,
    p: &[i16],
    gcode0: i16,
    exp_gcode0: i16,
    gain_pit: &mut i16,
    gain_cod: &mut i16,
    ov: &mut Flag,
) {
    *gain_pit = p[0];
    let g_code = p[1];
    let l_tmp = (g_code as i32 * gcode0 as i32).wrapping_shl(1);
    let l_tmp = l_shr(l_tmp, 10i16.wrapping_sub(exp_gcode0), ov);
    *gain_cod = (l_tmp >> 16) as i16;
    let (exp, frac) = log2(g_code as i32, ov);
    let exp = exp.wrapping_sub(12);
    let tmp = shr_r(frac, 5, ov);
    let qua_ener_mr122 = tmp.wrapping_add(((exp as i32) << 10) as i16);
    let l_tmp = mpy_32_16(exp, frac, 24660, ov).wrapping_shl(13);
    let qua_ener = (l_tmp.wrapping_add(0x8000) >> 16) as i16;
    gc_pred_update(pred_st, qua_ener_mr122, qua_ener);
}

/// `MR475_update_unq_pred`: updates the unquantised predictor with the optimum gain.
pub(super) fn mr475_update_unq_pred(
    pred_st: &mut GcPredState,
    exp_gcode0: i16,
    frac_gcode0: i16,
    mut cod_gain_exp: i16,
    mut cod_gain_frac: i16,
    ov: &mut Flag,
) {
    let (qua_ener, qua_ener_mr122);
    if cod_gain_frac <= 0 {
        qua_ener = MIN_QUA_ENER;
        qua_ener_mr122 = MIN_QUA_ENER_MR122;
    } else {
        let frac_gcode0 = pow2(14, frac_gcode0, ov) as i16;
        if cod_gain_frac >= frac_gcode0 {
            cod_gain_frac >>= 1;
            cod_gain_exp = cod_gain_exp.wrapping_add(1);
        }
        let frac = div_s(cod_gain_frac, frac_gcode0);
        let tmp = cod_gain_exp.wrapping_sub(exp_gcode0).wrapping_sub(1);
        let (exp, frac) = log2(frac as i32, ov);
        let exp = exp.wrapping_add(tmp);
        let q = shr_r(frac, 5, ov).wrapping_add(((exp as i32) << 10) as i16);
        if q > MAX_QUA_ENER_MR122 {
            qua_ener = MAX_QUA_ENER;
            qua_ener_mr122 = MAX_QUA_ENER_MR122;
        } else {
            let l_tmp = mpy_32_16(exp, frac, 24660, ov);
            let l_tmp = l_shl(l_tmp, 13, ov);
            qua_ener = round(l_tmp, ov);
            qua_ener_mr122 = q;
        }
    }
    gc_pred_update(pred_st, qua_ener_mr122, qua_ener);
}

/// `TABLE_GAIN_MR475` component-major: g_pitch, g_code of subframe 0, then of subframe 1.
static MR475_SOA: Soa<4, 256> = soa(&TABLE_GAIN_MR475, 4);

/// Every gain in `TABLE_GAIN_MR475` is >= 0 (`mr475_search_portable` relies on it).
const _: () = {
    let mut i = 0;
    while i < TABLE_GAIN_MR475.len() {
        assert!(TABLE_GAIN_MR475[i] >= 0);
        i += 1;
    }
};

/// `MR475_gain_quant`'s search over the 256 rows, `MR475_W` rows in lanes (8 with AVX2, 16 on
/// SSE2, as measured). Wrapping i32 sums, as the C's `+`; a strict `<` minimum per lane, then
/// the first minimum over the lanes.
/// No `Mpy_32_16` here can saturate: every multiplier `n` is in [0, 32767] (table gains >= 0,
/// `gcode0 = Pow2(14, .) >= 0`), so `hi * n != 2^30`, and `lo >= 0` (bits 1..15 of a word), so
/// `2*hi*n + 2*((lo*n) >> 15)` stays below 2^31 when positive and can't wrap when negative.
#[inline(always)]
pub(crate) fn mr475_search_portable<const MR475_W: usize>(
    coeff: &[i16; 10],
    coeff_lo: &[i16; 10],
    gcode0: &[i16; 2],
    gp_limit: i16,
) -> Option<usize> {
    // First tier: every term is >= 0-weighted (gains, gcode0, the lo parts >= 0), so each
    // partial sum is bounded by sum |c| n_max over the table's column maxima: if that fits,
    // no add overflows for any row and the lanes skip their per-add flags.
    let nmax = |sf: usize| {
        let (gp, gc) = (
            MR475_MAX[2 * sf],
            (MR475_MAX[2 * sf + 1] * gcode0[sf] as i64) >> 15,
        );
        [(gp * gp) >> 15, gp, (gc * gc) >> 15, gc, (gc * gp) >> 15]
    };
    let bound: i64 = (0..2)
        .flat_map(|sf| {
            let n = nmax(sf);
            (0..5).map(move |k| {
                2 * (coeff[5 * sf + k] as i64).abs() * n[k]
                    + 2 * (((coeff_lo[5 * sf + k] as i64).abs() * n[k]) >> 15)
            })
        })
        .sum();
    if bound <= MAX_32 as i64 {
        mr475_lanes::<MR475_W, false>(coeff, coeff_lo, gcode0, gp_limit)
    } else {
        mr475_lanes::<MR475_W, true>(coeff, coeff_lo, gcode0, gp_limit)
    }
}

/// Column maxima (|.|) of `TABLE_GAIN_MR475` (g_pitch, g_code of sf0, then sf1).
const MR475_MAX: [i64; 4] = {
    let mut m = [0i64; 4];
    let mut i = 0;
    while i < TABLE_GAIN_MR475.len() {
        let v = TABLE_GAIN_MR475[i] as i64;
        let v = if v < 0 { -v } else { v };
        if v > m[i % 4] {
            m[i % 4] = v;
        }
        i += 1;
    }
    m
};

/// `mr475_search_portable`'s lanes; `CHECK`: flag each add's overflow (else proven absent).
#[inline(always)]
fn mr475_lanes<const MR475_W: usize, const CHECK: bool>(
    coeff: &[i16; 10],
    coeff_lo: &[i16; 10],
    gcode0: &[i16; 2],
    gp_limit: i16,
) -> Option<usize> {
    // The C's Mpy_32_16/Mac_32_16 chain, term by term in its order (sf0's five, then sf1's),
    // each as L_mac(s, hi, n) then L_mac(s, mult(lo, n), 1). Wrapping, with a per-lane flag for
    // any add that overflowed (where the C saturated): `None` if a row the C compares had one.
    #[inline(always)]
    fn add<const CHECK: bool>(a: i32, b: i32, bad: &mut bool) -> i32 {
        let t = a.wrapping_add(b);
        if CHECK {
            *bad |= ((a ^ t) & (b ^ t)) < 0;
        }
        t
    }
    debug_assert!(gcode0.iter().chain(coeff_lo).all(|&v| v >= 0));
    let (gp_limit, g0) = (gp_limit as i32, gcode0.map(|g| g as i32));
    let mut best = [MAX_32; MR475_W];
    let mut bidx = [0i32; MR475_W];
    let mut bad_any = false;
    for b in (0..256).step_by(MR475_W) {
        let mut dist = [0i32; MR475_W];
        let mut valid = [true; MR475_W];
        let mut bad = [false; MR475_W];
        for sf in 0..2 {
            let gp: &[i16; MR475_W] = MR475_SOA[2 * sf][b..b + MR475_W].try_into().unwrap();
            let gc: &[i16; MR475_W] = MR475_SOA[2 * sf + 1][b..b + MR475_W].try_into().unwrap();
            let (c, cl) = (&coeff[5 * sf..5 * sf + 5], &coeff_lo[5 * sf..5 * sf + 5]);
            for l in 0..MR475_W {
                // gains, gcode0 >= 0: none of these mults saturates
                let g_pitch = gp[l] as i32;
                let g_code = (gc[l] as i32 * g0[sf]) >> 15;
                let n = [
                    (g_pitch * g_pitch) >> 15,
                    g_pitch,
                    (g_code * g_code) >> 15,
                    g_code,
                    (g_code * g_pitch) >> 15,
                ];
                for k in 0..5 {
                    // |hi * n| < 2^30: L_mult doesn't saturate
                    dist[l] = add::<CHECK>(dist[l], (c[k] as i32 * n[k]) << 1, &mut bad[l]);
                    dist[l] =
                        add::<CHECK>(dist[l], ((cl[k] as i32 * n[k]) >> 15) << 1, &mut bad[l]);
                }
                // The C tests sf0 as a Word16 `g_pitch - gp_limit <= 0`: that only differs from
                // this for a negative limit, where sf1's test fails for every row anyway (gains
                // >= 0). One comparison for both keeps the lanes fast (a per-sf select: +27%).
                valid[l] &= g_pitch <= gp_limit;
            }
        }
        for l in 0..MR475_W {
            bad_any |= valid[l] & bad[l];
            let lt = valid[l] & (dist[l] < best[l]);
            best[l] = if lt { dist[l] } else { best[l] };
            bidx[l] = if lt { (b + l) as i32 } else { bidx[l] };
        }
    }
    if bad_any {
        return None;
    }
    Some(first_min(&best, &bidx) as usize)
}

/// The faithful `MR475_gain_quant` search (saturating operators), when the lanes decline.
#[cold]
fn mr475_search_exact(
    coeff: &[i16; 10],
    coeff_lo: &[i16; 10],
    gcode0: &[i16; 2],
    gp_limit: i16,
    ov: &mut Flag,
) -> usize {
    let (mut dist_min, mut index) = (MAX_32, 0usize);
    for (i, row) in TABLE_GAIN_MR475.as_chunks::<4>().0.iter().enumerate() {
        let terms = |sf: usize, ov: &mut Flag| -> [i16; 5] {
            let g_pitch = row[2 * sf];
            let g_code = mult(row[2 * sf + 1], gcode0[sf], ov);
            [
                mult(g_pitch, g_pitch, ov),
                g_pitch,
                mult(g_code, g_code, ov),
                g_code,
                mult(g_code, g_pitch, ov),
            ]
        };
        let n0 = terms(0, ov);
        let mut l_tmp = mpy_32_16(coeff[0], coeff_lo[0], n0[0], ov);
        for k in 1..5 {
            l_tmp = mac_32_16(l_tmp, coeff[k], coeff_lo[k], n0[k], ov);
        }
        if row[0] <= gp_limit && row[2] <= gp_limit {
            let n1 = terms(1, ov);
            for k in 0..5 {
                l_tmp = mac_32_16(l_tmp, coeff[5 + k], coeff_lo[5 + k], n1[k], ov);
            }
            if l_tmp < dist_min {
                dist_min = l_tmp;
                index = i;
            }
        }
    }
    index
}

/// `MR475_gain_quant`: joint VQ of the gains of two subframes (MR475).
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn mr475_gain_quant(
    pred_st: &mut GcPredState,
    sf0_exp_gcode0: i16,
    sf0_frac_gcode0: i16,
    sf0_exp_coeff: &[i16; 5],
    sf0_frac_coeff: &[i16; 5],
    sf0_exp_target_en: i16,
    mut sf0_frac_target_en: i16,
    sf1_code_nosharp: &[i16],
    sf1_exp_gcode0: i16,
    sf1_frac_gcode0: i16,
    sf1_exp_coeff: &[i16; 5],
    sf1_frac_coeff: &[i16; 5],
    sf1_exp_target_en: i16,
    mut sf1_frac_target_en: i16,
    gp_limit: i16,
    sf0_gain_pit: &mut i16,
    sf0_gain_cod: &mut i16,
    sf1_gain_pit: &mut i16,
    sf1_gain_cod: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let sf0_gcode0 = pow2(14, sf0_frac_gcode0, ov) as i16;
    let sf1_gcode0 = pow2(14, sf1_frac_gcode0, ov) as i16;
    let mut exp_max = [0i16; 10];
    let e = sf0_exp_gcode0.wrapping_sub(11);
    exp_max[0] = sf0_exp_coeff[0].wrapping_sub(13);
    exp_max[1] = sf0_exp_coeff[1].wrapping_sub(14);
    exp_max[2] = sf0_exp_coeff[2].wrapping_add(15i16.wrapping_add(e.wrapping_shl(1)));
    exp_max[3] = sf0_exp_coeff[3].wrapping_add(e);
    exp_max[4] = sf0_exp_coeff[4].wrapping_add(1i16.wrapping_add(e));
    let e = sf1_exp_gcode0.wrapping_sub(11);
    exp_max[5] = sf1_exp_coeff[0].wrapping_sub(13);
    exp_max[6] = sf1_exp_coeff[1].wrapping_sub(14);
    exp_max[7] = sf1_exp_coeff[2].wrapping_add(15i16.wrapping_add(e.wrapping_shl(1)));
    exp_max[8] = sf1_exp_coeff[3].wrapping_add(e);
    exp_max[9] = sf1_exp_coeff[4].wrapping_add(1i16.wrapping_add(e));
    let e = sf0_exp_target_en.wrapping_sub(sf1_exp_target_en);
    // Word16 `>>=` on the promoted int (shift counts stay well below 32 in practice).
    if e > 0 {
        sf1_frac_target_en = ((sf1_frac_target_en as i32) >> (e as u32).min(31)) as i16;
    } else {
        sf0_frac_target_en = ((sf0_frac_target_en as i32) >> ((-(e as i32)) as u32).min(31)) as i16;
    }
    let mut e = 0i16;
    let tmp = shr_r(sf1_frac_target_en, 1, ov);
    if tmp > sf0_frac_target_en {
        e = 1;
    } else {
        let tmp = ((sf0_frac_target_en as i32 + 3) >> 2) as i16;
        if tmp > sf1_frac_target_en {
            e = -1;
        }
    }
    for v in &mut exp_max[..5] {
        *v = v.wrapping_add(e);
    }
    let mut e = exp_max[0];
    for i in (1..10).rev() {
        if exp_max[i] > e {
            e = exp_max[i];
        }
    }
    e = e.wrapping_add(1);
    let mut coeff = [0i16; 10];
    let mut coeff_lo = [0i16; 10];
    for i in 0..10 {
        let f = if i < 5 {
            sf0_frac_coeff[i]
        } else {
            sf1_frac_coeff[i - 5]
        };
        let l_tmp = l_shr((f as i32) << 16, e.wrapping_sub(exp_max[i]), ov);
        coeff[i] = (l_tmp >> 16) as i16;
        coeff_lo[i] = ((l_tmp >> 1) - ((l_tmp >> 16) << 15)) as i16;
    }
    let gcode0 = [sf0_gcode0, sf1_gcode0];
    let index = mr475_search(&coeff, &coeff_lo, &gcode0, gp_limit)
        .unwrap_or_else(|| mr475_search_exact(&coeff, &coeff_lo, &gcode0, gp_limit, ov));
    let tmp = index << 2;
    mr475_quant_store_results(
        pred_st,
        &TABLE_GAIN_MR475[tmp..],
        sf0_gcode0,
        sf0_exp_gcode0,
        sf0_gain_pit,
        sf0_gain_cod,
        ov,
    );
    let (sf1_exp_gcode0, sf1_frac_gcode0, _, _) =
        gc_pred(pred_st, CodecMode::MR475, sf1_code_nosharp, ov);
    let sf1_gcode0 = pow2(14, sf1_frac_gcode0, ov) as i16;
    mr475_quant_store_results(
        pred_st,
        &TABLE_GAIN_MR475[tmp + 2..],
        sf1_gcode0,
        sf1_exp_gcode0,
        sf1_gain_pit,
        sf1_gain_cod,
        ov,
    );
    index as i16
}

/// Per-subframe gain outputs of `gainQuant`.
pub(crate) struct Gains {
    pub gain_pit: i16,
    pub gain_cod: i16,
}

/// `gainQuant`: quantises the pitch and codebook gains of one subframe (for MR475 the pair of
/// subframes, at the odd one; the even one only stores its analysis). Appends indices to `ana`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(crate) fn gain_quant(
    st: &mut GainQuantState,
    mode: CodecMode,
    res: &[i16],
    exc: &[i16],
    code: &[i16],
    xn: &[i16],
    xn2: &[i16],
    y1: &[i16],
    y2: &[i16],
    g_coeff: &[i16],
    even_subframe: bool,
    gp_limit: i16,
    sf0_gain_pit: &mut i16,
    sf0_gain_cod: &mut i16,
    gains: &mut Gains,
    ana: &mut [i16],
    ana_pos: &mut usize,
    ov: &mut Flag,
) {
    let mut frac_coeff = [0i16; 5];
    let mut exp_coeff = [0i16; 5];
    let (mut cod_gain_frac, mut cod_gain_exp) = (0i16, 0i16);
    if mode == CodecMode::MR475 {
        if even_subframe {
            st.gain_idx_pos = *ana_pos;
            *ana_pos += 1;
            st.gc_pred_unq_st = st.gc_pred_st;
            let (e, f, _, _) = gc_pred(&st.gc_pred_unq_st, mode, code, ov);
            st.sf0_exp_gcode0 = e;
            st.sf0_frac_gcode0 = f;
            calc_filt_energies(
                mode,
                xn,
                xn2,
                y1,
                y2,
                g_coeff,
                &mut st.sf0_frac_coeff,
                &mut st.sf0_exp_coeff,
                &mut cod_gain_frac,
                &mut cod_gain_exp,
                ov,
            );
            gains.gain_cod = shl(cod_gain_frac, cod_gain_exp.wrapping_add(1), ov);
            (st.sf0_exp_target_en, st.sf0_frac_target_en) = calc_target_energy(xn, ov);
            mr475_update_unq_pred(
                &mut st.gc_pred_unq_st,
                st.sf0_exp_gcode0,
                st.sf0_frac_gcode0,
                cod_gain_exp,
                cod_gain_frac,
                ov,
            );
        } else {
            let (exp_gcode0, frac_gcode0, _, _) = gc_pred(&st.gc_pred_unq_st, mode, code, ov);
            calc_filt_energies(
                mode,
                xn,
                xn2,
                y1,
                y2,
                g_coeff,
                &mut frac_coeff,
                &mut exp_coeff,
                &mut cod_gain_frac,
                &mut cod_gain_exp,
                ov,
            );
            let (exp_en, frac_en) = calc_target_energy(xn, ov);
            let st0 = (
                st.sf0_exp_gcode0,
                st.sf0_frac_gcode0,
                st.sf0_exp_coeff,
                st.sf0_frac_coeff,
                st.sf0_exp_target_en,
                st.sf0_frac_target_en,
            );
            ana[st.gain_idx_pos] = mr475_gain_quant(
                &mut st.gc_pred_st,
                st0.0,
                st0.1,
                &st0.2,
                &st0.3,
                st0.4,
                st0.5,
                code,
                exp_gcode0,
                frac_gcode0,
                &exp_coeff,
                &frac_coeff,
                exp_en,
                frac_en,
                gp_limit,
                sf0_gain_pit,
                sf0_gain_cod,
                &mut gains.gain_pit,
                &mut gains.gain_cod,
                ov,
            );
        }
    } else {
        let (exp_gcode0, frac_gcode0, exp_en, frac_en) = gc_pred(&st.gc_pred_st, mode, code, ov);
        let (mut qua_ener_mr122, mut qua_ener) = (0i16, 0i16);
        if mode == CodecMode::MR122 {
            gains.gain_cod = g_code(xn2, y2);
            ana[*ana_pos] = q_gain_code(
                mode,
                exp_gcode0,
                frac_gcode0,
                &mut gains.gain_cod,
                &mut qua_ener_mr122,
                &mut qua_ener,
                ov,
            );
            *ana_pos += 1;
        } else {
            calc_filt_energies(
                mode,
                xn,
                xn2,
                y1,
                y2,
                g_coeff,
                &mut frac_coeff,
                &mut exp_coeff,
                &mut cod_gain_frac,
                &mut cod_gain_exp,
                ov,
            );
            if mode == CodecMode::MR795 {
                mr795_gain_quant(
                    &mut st.adapt_st,
                    res,
                    exc,
                    code,
                    &frac_coeff,
                    &exp_coeff,
                    exp_en,
                    frac_en,
                    exp_gcode0,
                    frac_gcode0,
                    cod_gain_frac,
                    cod_gain_exp,
                    gp_limit,
                    &mut gains.gain_pit,
                    &mut gains.gain_cod,
                    &mut qua_ener_mr122,
                    &mut qua_ener,
                    ana,
                    ana_pos,
                    ov,
                );
            } else {
                ana[*ana_pos] = qua_gain(
                    mode,
                    exp_gcode0,
                    frac_gcode0,
                    &frac_coeff,
                    &exp_coeff,
                    gp_limit,
                    &mut gains.gain_pit,
                    &mut gains.gain_cod,
                    &mut qua_ener_mr122,
                    &mut qua_ener,
                    ov,
                );
                *ana_pos += 1;
            }
        }
        gc_pred_update(&mut st.gc_pred_st, qua_ener_mr122, qua_ener);
    }
}

#[cfg(test)]
pub(crate) mod vs_ref;
