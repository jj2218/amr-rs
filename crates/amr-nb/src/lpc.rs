//! LPC / LSF / LSP helpers of the decoder: LSF<->LSP conversion, LSF reordering, LSP -> A(z),
//! LPC interpolation, LSF interpolation, LSF dequantisation (`D_plsf_3`, `D_plsf_5`) and LSP
//! averaging.
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `common/src/{lsp_lsf,reorder,lsp_az,
//! int_lpc}.cpp` and `dec/src/{int_lsf,d_plsf,d_plsf_3,d_plsf_5,lsp_avg}.cpp`.
#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1

use crate::basic_op::{
    Flag, add, l_add, l_extract, l_mac, l_msu, l_mult, l_shl, l_shr, l_shr_r, l_sub, mpy_32_16,
    mult, neg_shr, negate, round, shl, sub,
};
use crate::tables::{
    DICO1_LSF_3, DICO1_LSF_5, DICO2_LSF_3, DICO2_LSF_5, DICO3_LSF_3, DICO3_LSF_5, DICO4_LSF_5,
    DICO5_LSF_5, MEAN_LSF_3, MEAN_LSF_5, MR515_3_LSF, MR795_1_LSF, PAST_RQ_INIT, PRED_FAC_3, SLOPE,
    TABLE,
};

pub(crate) use crate::cnst::{M, MP1, Mode};
/// `LSF_GAP` in cnst.h.
const LSF_GAP: i16 = 205;
/// `LSP_PRED_FAC_MR122` in cnst.h.
const LSP_PRED_FAC_MR122: i16 = 21299;
/// `EXPCONST` in lsp_avg.h.
const EXPCONST: i16 = 5243;

const DICO1_SIZE: i16 = 256;
const DICO2_SIZE: i16 = 512;
const DICO3_SIZE: i16 = 512;
const MR515_3_SIZE: i16 = 128;
const MR795_1_SIZE: i16 = 512;

// ---------------------------------------------------------------------------------------------
// lsp_lsf.cpp
// ---------------------------------------------------------------------------------------------

/// `Lsf_lsp`: LSF (0..=16383) to LSP (cosine domain). Precondition: `lsf[i] >> 8` in 0..=63.
/// TS 26.073's operators (the final `add` saturates; opencore wraps).
pub(crate) fn lsf_lsp(lsf: &[i16], lsp: &mut [i16], m: usize, ov: &mut Flag) {
    for i in 0..m {
        let ind = (lsf[i] >> 8) as usize;
        let offset = lsf[i] & 0x00ff;
        let l_tmp = l_mult(sub(TABLE[ind + 1], TABLE[ind], ov), offset, ov);
        lsp[i] = add(TABLE[ind], l_shr(l_tmp, 9, ov) as i16, ov);
    }
}

/// `Lsp_lsf`: LSP (cosine domain) to LSF, with TS 26.073's saturating `sub`, `L_shl`, `round`
/// and `add` (opencore wraps: different for LSPs at the -1 edge).
pub(crate) fn lsp_lsf(lsp: &[i16], lsf: &mut [i16], m: usize, ov: &mut Flag) {
    let mut ind: usize = 63;
    for i in (0..m).rev() {
        while sub(TABLE[ind], lsp[i], ov) < 0 {
            ind -= 1;
        }
        let l_tmp = l_mult(sub(lsp[i], TABLE[ind], ov), SLOPE[ind], ov);
        let v = round(l_shl(l_tmp, 3, ov), ov);
        lsf[i] = add(v, shl((ind) as i16, 8, ov), ov);
    }
}

// ---------------------------------------------------------------------------------------------
// reorder.cpp
// ---------------------------------------------------------------------------------------------

/// `Reorder_lsf`: enforce a minimum distance between consecutive LSFs (saturating `add`).
pub(crate) fn reorder_lsf(lsf: &mut [i16], min_dist: i16, n: usize, ov: &mut Flag) {
    let mut lsf_min = min_dist;
    for v in lsf.iter_mut().take(n) {
        if sub(*v, lsf_min, ov) < 0 {
            *v = lsf_min;
        }
        lsf_min = add(*v, min_dist, ov);
    }
}

// ---------------------------------------------------------------------------------------------
// lsp_az.cpp
// ---------------------------------------------------------------------------------------------

/// `Get_lsp_pol`: the polynomial F1(z) or F2(z) from every second LSP, TS 26.073's operators
/// (`Mpy_32_16`, saturating `L_shl`, `L_add`, `L_sub`, `L_msu`; opencore wraps).
fn get_lsp_pol(lsp: &[i16], f: &mut [i32; 6], ov: &mut Flag) {
    if !get_lsp_pol_fast(lsp, f) {
        get_lsp_pol_exact(lsp, f, ov);
    }
}

/// `get_lsp_pol` in one word per coefficient: the `L_Extract` pair times `x` is
/// `Mpy_32_16` = 2 (((f >> 1) x) >> 15) exactly (hi x 2^15 is a multiple of 2^15). Exact while
/// every intermediate stays below 2^30 in magnitude (then no operator saturates, `hi` is never
/// -32768), tracked as the OR of one's-complement magnitudes; returns false otherwise.
fn get_lsp_pol_fast(lsp: &[i16], f: &mut [i32; 6]) -> bool {
    let mut g = [0i64; 6];
    g[0] = 1 << 24;
    g[1] = -((lsp[0] as i64) << 10);
    let mut acc = 0i64;
    let (mut l, mut fi) = (2usize, 2usize);
    for i in 2..=5usize {
        g[fi] = g[fi - 2];
        for _ in 1..i {
            let t0 = (((g[fi - 1] >> 1) * lsp[l] as i64) >> 15) << 2;
            let a = g[fi] + g[fi - 2];
            let v = a - t0;
            acc |= (t0 ^ (t0 >> 63)) | (a ^ (a >> 63)) | (v ^ (v >> 63));
            g[fi] = v;
            fi -= 1;
        }
        g[fi] -= (lsp[l] as i64) << 10;
        acc |= g[fi] ^ (g[fi] >> 63);
        fi += i;
        l += 2;
    }
    *f = g.map(|v| v as i32);
    acc < 1 << 30
}

/// `Get_lsp_pol` as TS 26.073's operator chain: `get_lsp_pol`'s fallback.
#[cold]
fn get_lsp_pol_exact(lsp: &[i16], f: &mut [i32; 6], ov: &mut Flag) {
    let mut l = 0usize;
    f[0] = l_mult(4096, 2048, ov);
    f[1] = l_msu(0, lsp[l], 512, ov);
    l += 2;
    let mut fi = 2usize;
    for i in 2..=5usize {
        f[fi] = f[fi - 2];
        for _ in 1..i {
            let (hi, lo) = l_extract(f[fi - 1], ov);
            let t0 = l_shl(mpy_32_16(hi, lo, lsp[l], ov), 1, ov);
            f[fi] = l_add(f[fi], f[fi - 2], ov);
            f[fi] = l_sub(f[fi], t0, ov);
            fi -= 1;
        }
        f[fi] = l_msu(f[fi], lsp[l], 512, ov);
        fi += i;
        l += 2;
    }
}

/// `Lsp_Az`: LSPs (`M` values) to LPC coefficients `a[0..=M]` (Q12).
pub(crate) fn lsp_az(lsp: &[i16], a: &mut [i16], ov: &mut Flag) {
    let mut f1 = [0i32; 6];
    let mut f2 = [0i32; 6];
    get_lsp_pol(&lsp[0..], &mut f1, ov);
    get_lsp_pol(&lsp[1..], &mut f2, ov);
    for i in (1..=5usize).rev() {
        f1[i] = l_add(f1[i], f1[i - 1], ov);
        f2[i] = l_sub(f2[i], f2[i - 1], ov);
    }
    a[0] = 4096;
    let mut j = 10usize;
    for i in 1..=5usize {
        a[i] = l_shr_r(l_add(f1[i], f2[i], ov), 13, ov) as i16;
        a[j] = l_shr_r(l_sub(f1[i], f2[i], ov), 13, ov) as i16;
        j -= 1;
    }
}

// ---------------------------------------------------------------------------------------------
// int_lpc.cpp
// ---------------------------------------------------------------------------------------------

/// `(x >> 1) + (y >> 1)` per element, truncated to Word16 as in C.
fn half_sum(a: &[i16], b: &[i16], out: &mut [i16; M]) {
    for i in 0..M {
        out[i] = (a[i] >> 1).wrapping_add(b[i] >> 1);
    }
}

/// `Int_lpc_1and3`: interpolated A(z) for subframes 1 and 3, quantised for 2 and 4.
/// `az` holds 4 * MP1 coefficients.
pub(crate) fn int_lpc_1and3(
    lsp_old: &[i16],
    lsp_mid: &[i16],
    lsp_new: &[i16],
    az: &mut [i16],
    ov: &mut Flag,
) {
    let mut lsp = [0i16; M];
    half_sum(lsp_old, lsp_mid, &mut lsp);
    lsp_az(&lsp, &mut az[0..], ov);
    lsp_az(lsp_mid, &mut az[MP1..], ov);
    half_sum(lsp_mid, lsp_new, &mut lsp);
    lsp_az(&lsp, &mut az[2 * MP1..], ov);
    lsp_az(lsp_new, &mut az[3 * MP1..], ov);
}

/// `Int_lpc_1and3_2`: only the interpolated A(z) for subframes 1 and 3 (at `az[0]`, `az[2*MP1]`).
pub(crate) fn int_lpc_1and3_2(
    lsp_old: &[i16],
    lsp_mid: &[i16],
    lsp_new: &[i16],
    az: &mut [i16],
    ov: &mut Flag,
) {
    let mut lsp = [0i16; M];
    half_sum(lsp_old, lsp_mid, &mut lsp);
    lsp_az(&lsp, &mut az[0..], ov);
    half_sum(lsp_mid, lsp_new, &mut lsp);
    lsp_az(&lsp, &mut az[2 * MP1..], ov);
}

/// The three interpolated LSP vectors of `Int_lpc_1to3{,_2}` (weights 3/4+1/4, 1/2+1/2, 1/4+3/4).
fn int_1to3_step(lsp_old: &[i16], lsp_new: &[i16], step: usize, lsp: &mut [i16; M]) {
    for i in 0..M {
        lsp[i] = match step {
            0 => {
                let temp = lsp_old[i].wrapping_sub(lsp_old[i] >> 2);
                temp.wrapping_add(lsp_new[i] >> 2)
            }
            1 => (lsp_new[i] >> 1).wrapping_add(lsp_old[i] >> 1),
            _ => {
                let temp = lsp_new[i].wrapping_sub(lsp_new[i] >> 2);
                temp.wrapping_add(lsp_old[i] >> 2)
            }
        };
    }
}

/// `Int_lpc_1to3`: interpolated A(z) for subframes 1..3, quantised for subframe 4.
pub(crate) fn int_lpc_1to3(lsp_old: &[i16], lsp_new: &[i16], az: &mut [i16], ov: &mut Flag) {
    let mut lsp = [0i16; M];
    for step in 0..3 {
        int_1to3_step(lsp_old, lsp_new, step, &mut lsp);
        lsp_az(&lsp, &mut az[step * MP1..], ov);
    }
    lsp_az(lsp_new, &mut az[3 * MP1..], ov);
}

/// `Int_lpc_1to3_2`: interpolated A(z) for subframes 1..3 only.
pub(crate) fn int_lpc_1to3_2(lsp_old: &[i16], lsp_new: &[i16], az: &mut [i16], ov: &mut Flag) {
    let mut lsp = [0i16; M];
    for step in 0..3 {
        int_1to3_step(lsp_old, lsp_new, step, &mut lsp);
        lsp_az(&lsp, &mut az[step * MP1..], ov);
    }
}

// ---------------------------------------------------------------------------------------------
// int_lsf.cpp
// ---------------------------------------------------------------------------------------------

/// `Int_lsf`: LSF interpolation for subframe start `i_subfr` (0, 40, 80, 120; others: no-op).
pub(crate) fn int_lsf(
    lsf_old: &[i16],
    lsf_new: &[i16],
    i_subfr: i16,
    lsf_out: &mut [i16],
    ov: &mut Flag,
) {
    match i_subfr {
        0 => {
            for i in (0..M).rev() {
                let temp1 = neg_shr(lsf_old[i] as i32, 2) as i16;
                let temp2 = neg_shr(lsf_new[i] as i32, 2) as i16;
                lsf_out[i] = add(lsf_old[i].wrapping_sub(temp1), temp2, ov);
            }
        }
        40 => {
            for i in (0..M).rev() {
                let temp1 = neg_shr(lsf_old[i] as i32, 1) as i16;
                let temp2 = neg_shr(lsf_new[i] as i32, 1) as i16;
                lsf_out[i] = temp1.wrapping_add(temp2);
            }
        }
        80 => {
            for i in (0..M).rev() {
                let temp1 = neg_shr(lsf_old[i] as i32, 2) as i16;
                let temp2 = neg_shr(lsf_new[i] as i32, 2) as i16;
                lsf_out[i] = add(temp1, lsf_new[i].wrapping_sub(temp2), ov);
            }
        }
        120 => {
            for i in (0..M).rev() {
                lsf_out[i] = lsf_new[i];
            }
        }
        _ => {}
    }
}

// ---------------------------------------------------------------------------------------------
// d_plsf.cpp, d_plsf_3.cpp, d_plsf_5.cpp
// ---------------------------------------------------------------------------------------------

/// `D_plsfState`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct DPlsfState {
    /// Past quantised prediction residual.
    pub past_r_q: [i16; M],
    /// Past dequantised LSFs.
    pub past_lsf_q: [i16; M],
}

impl DPlsfState {
    /// `D_plsf_reset` (with `mean_lsf_5`, as the decoder passes it).
    pub(crate) fn reset(&mut self) {
        self.past_r_q = [0; M];
        self.past_lsf_q = MEAN_LSF_5;
    }
}

/// `Init_D_plsf_3`: load `past_r_q` from `past_rq_init[index * M..]` (index 0..=7).
pub(crate) fn init_d_plsf_3(st: &mut DPlsfState, index: i16) {
    let o = index as usize * M;
    st.past_r_q.copy_from_slice(&PAST_RQ_INIT[o..o + M]);
}

/// `D_plsf_3`: decode the 3 split-VQ LSF indices (all modes except MR122) into `lsp1_q`.
pub(crate) fn d_plsf_3(
    st: &mut DPlsfState,
    mode: Mode,
    bfi: i16,
    indice: &[i16],
    lsp1_q: &mut [i16],
    ov: &mut Flag,
) {
    const ALPHA: i16 = 29491;
    const ONE_ALPHA: i16 = 3277;
    let mut lsf1_r = [0i16; M];
    let mut lsf1_q = [0i16; M];

    if bfi != 0 {
        for i in 0..M {
            let temp = mult(st.past_lsf_q[i], ALPHA, ov);
            let index = mult(MEAN_LSF_3[i], ONE_ALPHA, ov);
            lsf1_q[i] = add(index, temp, ov);
        }
        if mode != Mode::MRDTX {
            for i in 0..M {
                let mut temp = mult(st.past_r_q[i], PRED_FAC_3[i], ov);
                temp = add(MEAN_LSF_3[i], temp, ov);
                st.past_r_q[i] = sub(lsf1_q[i], temp, ov);
            }
        } else {
            for i in 0..M {
                let temp = add(MEAN_LSF_3[i], st.past_r_q[i], ov);
                st.past_r_q[i] = sub(lsf1_q[i], temp, ov);
            }
        }
    } else {
        let index_limit_2: i16 = (DICO2_SIZE - 1) * 3;
        let p_cb2: &[i16] = &DICO2_LSF_3;
        let (p_cb1, p_cb3, index_limit_1, index_limit_3): (&[i16], &[i16], i16, i16) =
            if mode == Mode::MR475 || mode == Mode::MR515 {
                (
                    &DICO1_LSF_3,
                    &MR515_3_LSF,
                    (DICO1_SIZE - 1) * 3,
                    (MR515_3_SIZE - 1) * 4,
                )
            } else if mode == Mode::MR795 {
                (
                    &MR795_1_LSF,
                    &DICO3_LSF_3,
                    (MR795_1_SIZE - 1) * 3,
                    (DICO3_SIZE - 1) * 4,
                )
            } else {
                (
                    &DICO1_LSF_3,
                    &DICO3_LSF_3,
                    (DICO1_SIZE - 1) * 3,
                    (DICO3_SIZE - 1) * 4,
                )
            };

        let mut index = indice[0];
        let mut temp = index.wrapping_add(index.wrapping_shl(1));
        if temp > index_limit_1 {
            temp = index_limit_1;
        }
        let p = temp as usize;
        lsf1_r[0..3].copy_from_slice(&p_cb1[p..p + 3]);

        index = indice[1];
        if mode == Mode::MR475 || mode == Mode::MR515 {
            index = index.wrapping_shl(1);
        }
        temp = index.wrapping_add(index.wrapping_shl(1));
        if temp > index_limit_2 {
            temp = index_limit_2;
        }
        let p = temp as usize;
        lsf1_r[3..6].copy_from_slice(&p_cb2[p..p + 3]);

        index = indice[2];
        temp = index.wrapping_shl(2);
        if temp > index_limit_3 {
            temp = index_limit_3;
        }
        let p = temp as usize;
        lsf1_r[6..10].copy_from_slice(&p_cb3[p..p + 4]);

        if mode != Mode::MRDTX {
            for i in 0..M {
                let mut temp = mult(st.past_r_q[i], PRED_FAC_3[i], ov);
                temp = add(MEAN_LSF_3[i], temp, ov);
                lsf1_q[i] = add(lsf1_r[i], temp, ov);
                st.past_r_q[i] = lsf1_r[i];
            }
        } else {
            for i in 0..M {
                let temp = add(MEAN_LSF_3[i], st.past_r_q[i], ov);
                lsf1_q[i] = add(lsf1_r[i], temp, ov);
                st.past_r_q[i] = lsf1_r[i];
            }
        }
    }

    reorder_lsf(&mut lsf1_q, LSF_GAP, M, ov);
    st.past_lsf_q = lsf1_q;
    lsf_lsp(&lsf1_q, lsp1_q, M, ov);
}

/// `D_plsf_5`: decode the 5 split-VQ LSF indices of MR122 into two LSP vectors.
pub(crate) fn d_plsf_5(
    st: &mut DPlsfState,
    bfi: i16,
    indice: &[i16],
    lsp1_q: &mut [i16],
    lsp2_q: &mut [i16],
    ov: &mut Flag,
) {
    const ALPHA: i32 = 31128;
    const ONE_ALPHA: i32 = 1639;
    let mut lsf1_r = [0i16; M];
    let mut lsf2_r = [0i16; M];
    let mut lsf1_q = [0i16; M];
    let mut lsf2_q = [0i16; M];

    if bfi != 0 {
        for i in 0..M {
            let temp = ((st.past_lsf_q[i] as i32 * ALPHA) >> 15) as i16;
            let sign = ((MEAN_LSF_5[i] as i32 * ONE_ALPHA) >> 15) as i16;
            lsf1_q[i] = add(sign, temp, ov);
            lsf2_q[i] = lsf1_q[i];
            let temp = ((st.past_r_q[i] as i32 * LSP_PRED_FAC_MR122 as i32) >> 15) as i16;
            let temp = add(MEAN_LSF_5[i], temp, ov);
            st.past_r_q[i] = sub(lsf2_q[i], temp, ov);
        }
    } else {
        let p = shl(indice[0], 2, ov) as usize;
        let d = &DICO1_LSF_5[p..p + 4];
        (lsf1_r[0], lsf1_r[1], lsf2_r[0], lsf2_r[1]) = (d[0], d[1], d[2], d[3]);

        let p = shl(indice[1], 2, ov) as usize;
        let d = &DICO2_LSF_5[p..p + 4];
        (lsf1_r[2], lsf1_r[3], lsf2_r[2], lsf2_r[3]) = (d[0], d[1], d[2], d[3]);

        let sign = indice[2] & 1;
        let i = neg_shr(indice[2] as i32, 1) as i16;
        let p = shl(i, 2, ov) as usize;
        let d = &DICO3_LSF_5[p..p + 4];
        if sign == 0 {
            (lsf1_r[4], lsf1_r[5], lsf2_r[4], lsf2_r[5]) = (d[0], d[1], d[2], d[3]);
        } else {
            lsf1_r[4] = negate(d[0]);
            lsf1_r[5] = negate(d[1]);
            lsf2_r[4] = negate(d[2]);
            lsf2_r[5] = negate(d[3]);
        }

        let p = shl(indice[3], 2, ov) as usize;
        let d = &DICO4_LSF_5[p..p + 4];
        (lsf1_r[6], lsf1_r[7], lsf2_r[6], lsf2_r[7]) = (d[0], d[1], d[2], d[3]);

        let p = shl(indice[4], 2, ov) as usize;
        let d = &DICO5_LSF_5[p..p + 4];
        (lsf1_r[8], lsf1_r[9], lsf2_r[8], lsf2_r[9]) = (d[0], d[1], d[2], d[3]);

        for i in 0..M {
            let mut temp = mult(st.past_r_q[i], LSP_PRED_FAC_MR122, ov);
            temp = add(MEAN_LSF_5[i], temp, ov);
            lsf1_q[i] = add(lsf1_r[i], temp, ov);
            lsf2_q[i] = add(lsf2_r[i], temp, ov);
            st.past_r_q[i] = lsf2_r[i];
        }
    }

    reorder_lsf(&mut lsf1_q, LSF_GAP, M, ov);
    reorder_lsf(&mut lsf2_q, LSF_GAP, M, ov);
    st.past_lsf_q = lsf2_q;
    lsf_lsp(&lsf1_q, lsp1_q, M, ov);
    lsf_lsp(&lsf2_q, lsp2_q, M, ov);
}

// ---------------------------------------------------------------------------------------------
// lsp_avg.cpp
// ---------------------------------------------------------------------------------------------

/// `lsp_avgState`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct LspAvgState {
    /// Averaged LSPs (`lsp_meanSave`).
    pub lsp_mean_save: [i16; M],
}

impl LspAvgState {
    /// `lsp_avg_reset` (with `mean_lsf_5`, as the decoder passes it).
    pub(crate) fn reset(&mut self) {
        self.lsp_mean_save = MEAN_LSF_5;
    }
}

/// `lsp_avg`: first-order recursive average of the LSPs.
pub(crate) fn lsp_avg(st: &mut LspAvgState, lsp: &[i16], ov: &mut Flag) {
    for i in 0..M {
        let mut l_tmp = (st.lsp_mean_save[i] as i32) << 16;
        l_tmp = l_msu(l_tmp, EXPCONST, st.lsp_mean_save[i], ov);
        l_tmp = l_mac(l_tmp, EXPCONST, lsp[i], ov);
        st.lsp_mean_save[i] = round(l_tmp, ov);
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::testutil::{N, Rng};
    use opencore_sys as _;
    use std::ffi::c_void;

    #[repr(C)]
    struct CommonAmrTbls {
        ptrs: [*const c_void; 34],
    }

    #[repr(C)]
    #[derive(Default)]
    struct CDPlsfState {
        past_r_q: [i16; M],
        past_lsf_q: [i16; M],
    }

    #[repr(C)]
    #[derive(Default)]
    struct CLspAvgState {
        lsp_mean_save: [i16; M],
    }

    // The TS 26.073 reference, in the shape of the opencore functions these tests were written
    // against (`o`/`f`: its overflow flag; table and mean arguments unused: the reference has
    // its own). opencore only for the table set the tests pass around. The reference's flag is
    // global and sticky, so the asserts compare values and state only.
    unsafe extern "C" {
        fn get_const_tbls(t: *mut std::ffi::c_void); // opaque: gains declares it too
        fn ref_Lsf_lsp(lsf: *const i16, lsp: *mut i16, m: i16);
        fn ref_Lsp_lsf(lsp: *const i16, lsf: *mut i16, m: i16);
        fn ref_Reorder_lsf(lsf: *mut i16, min_dist: i16, n: i16);
        fn ref_Lsp_Az(lsp: *const i16, a: *mut i16);
        fn ref_Int_lpc_1and3(o: *const i16, m: *const i16, n: *const i16, az: *mut i16);
        fn ref_Int_lpc_1and3_2(o: *const i16, m: *const i16, n: *const i16, az: *mut i16);
        fn ref_Int_lpc_1to3(o: *const i16, n: *const i16, az: *mut i16);
        fn ref_Int_lpc_1to3_2(o: *const i16, n: *const i16, az: *mut i16);
        fn ref_Int_lsf(o: *const i16, n: *const i16, i_subfr: i16, out: *mut i16);
        fn ref_D_plsf_reset(st: *mut CDPlsfState) -> i32;
        fn ref_Init_D_plsf_3(st: *mut CDPlsfState, index: i16);
        fn ref_D_plsf_3(st: *mut CDPlsfState, mode: i32, bfi: i16, ind: *const i16, q: *mut i16);
        fn ref_D_plsf_5(
            st: *mut CDPlsfState,
            bfi: i16,
            indice: *const i16,
            lsp1_q: *mut i16,
            lsp2_q: *mut i16,
        ) -> i32;
        fn ref_lsp_avg_reset(st: *mut CLspAvgState) -> i32;
        fn ref_lsp_avg(st: *mut CLspAvgState, lsp: *const i16);
    }
    use crate::testutil::with_ref;
    const NO: *mut i32 = std::ptr::null_mut();
    #[allow(non_snake_case)]
    unsafe fn Lsf_lsp(lsf: *const i16, lsp: *mut i16, m: i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_Lsf_lsp(lsf, lsp, m) })
    }
    #[allow(non_snake_case)]
    unsafe fn Lsp_lsf(lsp: *const i16, lsf: *mut i16, m: i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_Lsp_lsf(lsp, lsf, m) })
    }
    #[allow(non_snake_case)]
    unsafe fn Reorder_lsf(lsf: *mut i16, min_dist: i16, n: i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_Reorder_lsf(lsf, min_dist, n) })
    }
    #[allow(non_snake_case)]
    unsafe fn Lsp_Az(lsp: *const i16, a: *mut i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_Lsp_Az(lsp, a) })
    }
    #[allow(non_snake_case)]
    unsafe fn Int_lpc_1and3(
        o: *const i16,
        m: *const i16,
        n: *const i16,
        az: *mut i16,
        f: *mut i32,
    ) {
        with_ref(f, || unsafe { ref_Int_lpc_1and3(o, m, n, az) })
    }
    #[allow(non_snake_case)]
    unsafe fn Int_lpc_1and3_2(
        o: *const i16,
        m: *const i16,
        n: *const i16,
        az: *mut i16,
        f: *mut i32,
    ) {
        with_ref(f, || unsafe { ref_Int_lpc_1and3_2(o, m, n, az) })
    }
    #[allow(non_snake_case)]
    unsafe fn Int_lpc_1to3(o: *const i16, n: *const i16, az: *mut i16, f: *mut i32) {
        with_ref(f, || unsafe { ref_Int_lpc_1to3(o, n, az) })
    }
    #[allow(non_snake_case)]
    unsafe fn Int_lpc_1to3_2(o: *const i16, n: *const i16, az: *mut i16, f: *mut i32) {
        with_ref(f, || unsafe { ref_Int_lpc_1to3_2(o, n, az) })
    }
    #[allow(non_snake_case)]
    unsafe fn Int_lsf(o: *const i16, n: *const i16, i_subfr: i16, out: *mut i16, f: *mut i32) {
        with_ref(f, || unsafe { ref_Int_lsf(o, n, i_subfr, out) })
    }
    #[allow(non_snake_case)]
    unsafe fn D_plsf_reset(st: *mut CDPlsfState, _mean_lsf_5: *const i16) -> i16 {
        with_ref(NO, || unsafe { ref_D_plsf_reset(st) }) as i16
    }
    #[allow(non_snake_case)]
    unsafe fn Init_D_plsf_3(st: *mut CDPlsfState, index: i16, _past_rq_init: *const i16) {
        with_ref(NO, || unsafe { ref_Init_D_plsf_3(st, index) })
    }
    #[allow(non_snake_case)]
    unsafe fn D_plsf_3(
        st: *mut CDPlsfState,
        mode: i32,
        bfi: i16,
        indice: *const i16,
        _t: *mut CommonAmrTbls,
        lsp1_q: *mut i16,
        o: *mut i32,
    ) {
        with_ref(o, || unsafe { ref_D_plsf_3(st, mode, bfi, indice, lsp1_q) })
    }
    #[allow(non_snake_case)]
    unsafe fn D_plsf_5(
        st: *mut CDPlsfState,
        bfi: i16,
        indice: *const i16,
        _t: *mut CommonAmrTbls,
        lsp1_q: *mut i16,
        lsp2_q: *mut i16,
        o: *mut i32,
    ) {
        with_ref(o, || unsafe {
            ref_D_plsf_5(st, bfi, indice, lsp1_q, lsp2_q)
        });
    }
    unsafe fn lsp_avg_reset(st: *mut CLspAvgState, _mean_lsf_5: *const i16) -> i16 {
        with_ref(NO, || unsafe { ref_lsp_avg_reset(st) }) as i16
    }
    unsafe fn c_lsp_avg(st: *mut CLspAvgState, lsp: *const i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_lsp_avg(st, lsp) })
    }

    const MEAN_LSF_5_IDX: usize = 12;
    const PAST_RQ_INIT_IDX: usize = 15;

    fn tbls() -> CommonAmrTbls {
        let mut t = CommonAmrTbls {
            ptrs: [std::ptr::null(); 34],
        };
        unsafe { get_const_tbls((&raw mut t).cast()) };
        t
    }

    const MODES: [Mode; 9] = [
        Mode::MR475,
        Mode::MR515,
        Mode::MR59,
        Mode::MR67,
        Mode::MR74,
        Mode::MR795,
        Mode::MR102,
        Mode::MR122,
        Mode::MRDTX,
    ];

    /// Random LSF vector in the `Lsf_lsp` domain (0..=16383); sorted (realistic) or not.
    fn rand_lsf(g: &mut Rng) -> [i16; M] {
        let mut lsf = [0i16; M];
        for v in &mut lsf {
            *v = match g.u64() % 8 {
                0 => [0, 16383, 255, 256][(g.u64() % 4) as usize],
                _ => (g.u64() % 16384) as i16,
            };
        }
        if !g.u64().is_multiple_of(4) {
            lsf.sort();
        }
        lsf
    }

    /// Realistic LSP vector (via `Lsf_lsp` of a sorted LSF vector).
    fn rand_lsp(g: &mut Rng) -> [i16; M] {
        let mut lsf = rand_lsf(g);
        lsf.sort();
        let mut lsp = [0i16; M];
        lsf_lsp(&lsf, &mut lsp, M, &mut false);
        lsp
    }

    #[test]
    fn lpc_lsf_lsp_vs_ref() {
        let mut g = Rng::new(1);
        for _ in 0..N {
            let lsf = rand_lsf(&mut g);
            let (mut r, mut c) = ([0i16; M], [0i16; M]);
            let (mut ro, mut co) = (false, 0);
            lsf_lsp(&lsf, &mut r, M, &mut ro);
            unsafe { Lsf_lsp(lsf.as_ptr(), c.as_mut_ptr(), M as i16, &mut co) };
            assert_eq!((r, ro), (c, ro), "Lsf_lsp({lsf:?})");

            // Lsp_lsf accepts any Word16 input (table[0] == 32767 stops the search).
            let mut lsp = [0i16; M];
            if g.u64().is_multiple_of(2) {
                lsp = r;
            } else {
                lsp.iter_mut().for_each(|v| *v = g.i16());
            }
            let (mut r, mut c) = ([0i16; M], [0i16; M]);
            lsp_lsf(&lsp, &mut r, M, &mut ro);
            unsafe { Lsp_lsf(lsp.as_ptr(), c.as_mut_ptr(), M as i16, &mut co) };
            assert_eq!((r, ro), (c, ro), "Lsp_lsf({lsp:?})");
        }
    }

    #[test]
    fn lpc_reorder_lsf_vs_ref() {
        let mut g = Rng::new(2);
        for _ in 0..N {
            let mut r = [0i16; M];
            if g.u64().is_multiple_of(2) {
                r.iter_mut().for_each(|v| *v = g.i16());
            } else {
                r = rand_lsf(&mut g);
            }
            let mut c = r;
            let min_dist = match g.u64() % 4 {
                0 => g.i16(),
                _ => (g.u64() % 1024) as i16,
            };
            let n = (g.u64() % (M as u64 + 1)) as usize;
            let (mut ro, mut co) = (false, 0);
            reorder_lsf(&mut r, min_dist, n, &mut ro);
            unsafe { Reorder_lsf(c.as_mut_ptr(), min_dist, n as i16, &mut co) };
            assert_eq!((r, ro), (c, ro), "Reorder_lsf(min_dist={min_dist}, n={n})");
        }
    }

    /// Valid LSPs, and every 4th case arbitrary ones (large and wrapping polynomial coefficients).
    #[test]
    fn lpc_lsp_az_vs_ref() {
        let mut g = Rng::new(3);
        for k in 0..N {
            let lsp = if k % 4 == 0 {
                core::array::from_fn(|_| g.i16())
            } else {
                rand_lsp(&mut g)
            };
            let (mut r, mut c) = ([0x5555i16; MP1 + 1], [0x5555i16; MP1 + 1]);
            let (mut ro, mut co) = (false, 0);
            lsp_az(&lsp, &mut r, &mut ro);
            unsafe { Lsp_Az(lsp.as_ptr(), c.as_mut_ptr(), &mut co) };
            assert_eq!((r, ro), (c, ro), "Lsp_Az({lsp:?})");
        }
    }

    #[test]
    fn lpc_int_lpc_vs_ref() {
        let mut g = Rng::new(4);
        for _ in 0..N / 4 {
            let (old, mid, new) = (rand_lsp(&mut g), rand_lsp(&mut g), rand_lsp(&mut g));
            let fresh = || ([0x5555i16; 4 * MP1], [0x5555i16; 4 * MP1], false, 0i32);

            let (mut r, mut c, mut ro, mut co) = fresh();
            int_lpc_1and3(&old, &mid, &new, &mut r, &mut ro);
            unsafe {
                Int_lpc_1and3(
                    old.as_ptr(),
                    mid.as_ptr(),
                    new.as_ptr(),
                    c.as_mut_ptr(),
                    &mut co,
                )
            };
            assert_eq!((r, ro), (c, ro), "Int_lpc_1and3");

            let (mut r, mut c, mut ro, mut co) = fresh();
            int_lpc_1and3_2(&old, &mid, &new, &mut r, &mut ro);
            unsafe {
                Int_lpc_1and3_2(
                    old.as_ptr(),
                    mid.as_ptr(),
                    new.as_ptr(),
                    c.as_mut_ptr(),
                    &mut co,
                )
            };
            assert_eq!((r, ro), (c, ro), "Int_lpc_1and3_2");

            let (mut r, mut c, mut ro, mut co) = fresh();
            int_lpc_1to3(&old, &new, &mut r, &mut ro);
            unsafe { Int_lpc_1to3(old.as_ptr(), new.as_ptr(), c.as_mut_ptr(), &mut co) };
            assert_eq!((r, ro), (c, ro), "Int_lpc_1to3");

            let (mut r, mut c, mut ro, mut co) = fresh();
            int_lpc_1to3_2(&old, &new, &mut r, &mut ro);
            unsafe { Int_lpc_1to3_2(old.as_ptr(), new.as_ptr(), c.as_mut_ptr(), &mut co) };
            assert_eq!((r, ro), (c, ro), "Int_lpc_1to3_2");
        }
    }

    #[test]
    fn lpc_int_lsf_vs_ref() {
        let mut g = Rng::new(5);
        for _ in 0..N {
            let (mut old, mut new) = ([0i16; M], [0i16; M]);
            if g.u64().is_multiple_of(2) {
                old.iter_mut().for_each(|v| *v = g.i16());
                new.iter_mut().for_each(|v| *v = g.i16());
            } else {
                (old, new) = (rand_lsf(&mut g), rand_lsf(&mut g));
            }
            let i_subfr = [0, 40, 80, 120, 20, -40, 160][(g.u64() % 7) as usize];
            let (mut r, mut c) = ([0x5555i16; M], [0x5555i16; M]);
            let (mut ro, mut co) = (false, 0);
            int_lsf(&old, &new, i_subfr, &mut r, &mut ro);
            unsafe { Int_lsf(old.as_ptr(), new.as_ptr(), i_subfr, c.as_mut_ptr(), &mut co) };
            assert_eq!((r, ro), (c, ro), "Int_lsf(i_subfr={i_subfr})");
        }
    }

    fn same_state(r: &DPlsfState, c: &CDPlsfState) -> bool {
        r.past_r_q == c.past_r_q && r.past_lsf_q == c.past_lsf_q
    }

    #[test]
    fn lpc_d_plsf_reset_init_vs_ref() {
        let t = tbls();
        let mut g = Rng::new(6);
        for _ in 0..10_000 {
            let mut r = DPlsfState::default();
            r.past_r_q.iter_mut().for_each(|v| *v = g.i16());
            r.past_lsf_q.iter_mut().for_each(|v| *v = g.i16());
            let mut c = CDPlsfState {
                past_r_q: r.past_r_q,
                past_lsf_q: r.past_lsf_q,
            };
            if g.u64().is_multiple_of(2) {
                r.reset();
                let ret = unsafe { D_plsf_reset(&mut c, t.ptrs[MEAN_LSF_5_IDX].cast()) };
                assert_eq!(ret, 0);
            } else {
                let index = (g.u64() % 8) as i16;
                init_d_plsf_3(&mut r, index);
                unsafe { Init_D_plsf_3(&mut c, index, t.ptrs[PAST_RQ_INIT_IDX].cast()) };
            }
            assert!(same_state(&r, &c), "D_plsf_reset/Init_D_plsf_3");
        }
    }

    #[test]
    fn lpc_d_plsf_3_vs_ref() {
        let mut t = tbls();
        let mut g = Rng::new(7);
        let mut seen = [[0usize; 2]; 9];
        for seq in 0..200 {
            let mut r = DPlsfState::default();
            let mut c = CDPlsfState::default();
            r.reset();
            unsafe { D_plsf_reset(&mut c, t.ptrs[MEAN_LSF_5_IDX].cast()) };
            if seq % 2 == 1 {
                let index = (g.u64() % 8) as i16;
                init_d_plsf_3(&mut r, index);
                unsafe { Init_D_plsf_3(&mut c, index, t.ptrs[PAST_RQ_INIT_IDX].cast()) };
            }
            // A sticky mode per sequence for some runs, random switching for others.
            let sticky = MODES[(g.u64() % 9) as usize];
            for call in 0..1000 {
                let mode = if seq % 3 == 0 {
                    MODES[(g.u64() % 9) as usize]
                } else {
                    sticky
                };
                // Bad-frame bursts plus isolated bad frames.
                let bfi = i16::from(g.u64().is_multiple_of(5) || (call / 50) % 7 == 3);
                // The indices the mode's bitstream carries (8/8/7, 9/9/9 or 8/9/9 bits): the
                // reference reads its tables unchecked (opencore's clamp has no spec meaning).
                let bits: [u32; 3] = match mode {
                    Mode::MR475 | Mode::MR515 => [8, 8, 7],
                    Mode::MR795 => [9, 9, 9],
                    _ => [8, 9, 9],
                };
                let indice: [i16; 3] = std::array::from_fn(|k| (g.u64() % (1 << bits[k])) as i16);
                let (mut rl, mut cl) = ([0i16; M], [0i16; M]);
                let (mut ro, mut co) = (false, 0);
                d_plsf_3(&mut r, mode, bfi, &indice, &mut rl, &mut ro);
                unsafe {
                    D_plsf_3(
                        &mut c,
                        mode as i32,
                        bfi,
                        indice.as_ptr(),
                        &mut t,
                        cl.as_mut_ptr(),
                        &mut co,
                    )
                };
                assert_eq!(
                    (rl, ro),
                    (cl, ro),
                    "D_plsf_3 seq {seq} call {call} mode {mode:?} bfi {bfi} {indice:?}"
                );
                assert!(same_state(&r, &c), "D_plsf_3 state seq {seq} call {call}");
                seen[mode as usize][bfi as usize] += 1;
            }
        }
        assert!(
            seen.iter().flatten().all(|&n| n > 0),
            "uncovered mode/bfi: {seen:?}"
        );
    }

    #[test]
    fn lpc_d_plsf_5_vs_ref() {
        let mut t = tbls();
        let mut g = Rng::new(8);
        for seq in 0..200 {
            let mut r = DPlsfState::default();
            let mut c = CDPlsfState::default();
            r.reset();
            unsafe { D_plsf_reset(&mut c, t.ptrs[MEAN_LSF_5_IDX].cast()) };
            for call in 0..1000 {
                let bfi = i16::from(g.u64().is_multiple_of(5) || (call / 50) % 7 == 3);
                // Field ranges: dico1 128, dico2 256, dico3 256 (+ sign bit), dico4 256, dico5 64.
                let indice: [i16; 5] = [
                    (g.u64() % 128) as i16,
                    (g.u64() % 256) as i16,
                    (g.u64() % 512) as i16,
                    (g.u64() % 256) as i16,
                    (g.u64() % 64) as i16,
                ];
                let (mut r1, mut r2, mut c1, mut c2) = ([0i16; M], [0i16; M], [0i16; M], [0i16; M]);
                let (mut ro, mut co) = (false, 0);
                d_plsf_5(&mut r, bfi, &indice, &mut r1, &mut r2, &mut ro);
                unsafe {
                    D_plsf_5(
                        &mut c,
                        bfi,
                        indice.as_ptr(),
                        &mut t,
                        c1.as_mut_ptr(),
                        c2.as_mut_ptr(),
                        &mut co,
                    )
                };
                assert_eq!(
                    (r1, r2, ro),
                    (c1, c2, ro),
                    "D_plsf_5 seq {seq} call {call} bfi {bfi} {indice:?}"
                );
                assert!(same_state(&r, &c), "D_plsf_5 state seq {seq} call {call}");
            }
        }
    }

    #[test]
    fn lpc_lsp_avg_vs_ref() {
        let t = tbls();
        let mut g = Rng::new(9);
        for seq in 0..200 {
            let mut r = LspAvgState::default();
            let mut c = CLspAvgState::default();
            if seq % 4 == 0 {
                // Garbage start state, then reset.
                r.lsp_mean_save.iter_mut().for_each(|v| *v = g.i16());
                c.lsp_mean_save = r.lsp_mean_save;
            }
            r.reset();
            let ret = unsafe { lsp_avg_reset(&mut c, t.ptrs[MEAN_LSF_5_IDX].cast()) };
            assert_eq!(ret, 0);
            assert_eq!(r.lsp_mean_save, c.lsp_mean_save, "lsp_avg_reset");
            for call in 0..1000 {
                let lsp = if g.u64().is_multiple_of(2) {
                    rand_lsp(&mut g)
                } else {
                    std::array::from_fn(|_| g.i16())
                };
                let (mut ro, mut co) = (false, 0);
                lsp_avg(&mut r, &lsp, &mut ro);
                unsafe { c_lsp_avg(&mut c, lsp.as_ptr(), &mut co) };
                assert_eq!(ro, co != 0, "lsp_avg ov seq {seq} call {call}");
                assert_eq!(
                    r.lsp_mean_save, c.lsp_mean_save,
                    "lsp_avg seq {seq} call {call}"
                );
            }
        }
    }
}
