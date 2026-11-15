//! LPC parameters: ISF dequantisation (`qpisf_2s.cpp`, `qisf_ns.cpp`), ISF → ISP
//! (`isp_isf.cpp`), ISP → A(z) (`isp_az.cpp`), interpolation (`interpolate_isp.cpp`), the 6.60
//! kbit/s high-band ISF extrapolation (`isf_extrapolation.cpp`) and weighting
//! (`weight_amrwb_lpc.cpp`).

use crate::basic_op::*;
use crate::cnst::{ISF_GAP, L_MEANBUF, M, M16K};
use crate::tables::*;

const MU: i16 = 10923;
const ALPHA: i16 = 29491;
const ONE_ALPHA: i16 = (32768 - ALPHA as i32) as i16;

/// `Reorder_isf`: keeps neighbouring ISFs at least `min_dist` apart (all but the last).
pub(crate) fn reorder_isf(isf: &mut [i16], min_dist: i16) {
    let mut isf_min = min_dist;
    let n = isf.len();
    for v in &mut isf[..n - 1] {
        if *v < isf_min {
            *v = isf_min;
        }
        isf_min = add_int16(*v, min_dist);
    }
}

/// The bad-frame branch shared by both `Dpisf_2s_*`: ISFs from the past and the mean.
fn conceal_isf(isf_q: &mut [i16; M], past_isfq: &mut [i16; M], isfold: &[i16; M], isf_buf: &[i16]) {
    let mut ref_isf = [0i16; M];
    for (i, r) in ref_isf.iter_mut().enumerate() {
        let mut l = mul_16by16_to_int32(MEAN_ISF[i], 8192);
        for j in 0..L_MEANBUF {
            l = mac_16by16_to_int32(l, isf_buf[j * M + i], 8192);
        }
        *r = amr_wb_round(l);
    }
    for i in 0..M {
        isf_q[i] = add_int16(
            mult_int16(ALPHA, isfold[i]),
            mult_int16(ONE_ALPHA, ref_isf[i]),
        );
    }
    for i in 0..M {
        let tmp = add_int16(ref_isf[i], mult_int16(past_isfq[i], MU));
        past_isfq[i] = sub_int16(isf_q[i], tmp) >> 1;
    }
}

/// Shifts the ISF history (`isf_buf`, `L_MEANBUF` frames) and stores `isf_q` first.
fn push_isf_buf(isf_buf: &mut [i16], isf_q: &[i16; M]) {
    isf_buf.copy_within(0..(L_MEANBUF - 1) * M, M);
    isf_buf[..M].copy_from_slice(isf_q);
}

/// `Dpisf_2s_46b`: 46-bit ISF dequantisation (modes above 6.60 kbit/s), decoder side.
pub(crate) fn dpisf_2s_46b(
    indice: &[i16],
    isf_q: &mut [i16; M],
    past_isfq: &mut [i16; M],
    isfold: &[i16; M],
    isf_buf: &mut [i16],
    bfi: bool,
) {
    if !bfi {
        let ix = |k: usize| indice[k] as usize;
        for i in 0..9 {
            isf_q[i] = DICO1_ISF[ix(0) * 9 + i];
        }
        for i in 0..7 {
            isf_q[i + 9] = DICO2_ISF[ix(1) * 7 + i];
        }
        // plain `+=` on int16 in the C (wrapping)
        for i in 0..3 {
            isf_q[i] = isf_q[i].wrapping_add(DICO21_ISF[ix(2) * 3 + i]);
            isf_q[i + 3] = isf_q[i + 3].wrapping_add(DICO22_ISF[ix(3) * 3 + i]);
            isf_q[i + 6] = isf_q[i + 6].wrapping_add(DICO23_ISF[ix(4) * 3 + i]);
            isf_q[i + 9] = isf_q[i + 9].wrapping_add(DICO24_ISF[ix(5) * 3 + i]);
            isf_q[i + 12] = isf_q[i + 12].wrapping_add(DICO25_ISF[ix(6) * 4 + i]);
        }
        isf_q[15] = isf_q[15].wrapping_add(DICO25_ISF[ix(6) * 4 + 3]);
        for i in 0..M {
            let tmp = isf_q[i];
            isf_q[i] = isf_q[i].wrapping_add(MEAN_ISF[i]);
            isf_q[i] = isf_q[i].wrapping_add(((MU as i32 * past_isfq[i] as i32) >> 15) as i16);
            past_isfq[i] = tmp;
        }
        push_isf_buf(isf_buf, isf_q);
    } else {
        conceal_isf(isf_q, past_isfq, isfold, isf_buf);
    }
    reorder_isf(isf_q, ISF_GAP);
}

/// `Dpisf_2s_36b`: 36-bit ISF dequantisation (6.60 kbit/s), decoder side.
pub(crate) fn dpisf_2s_36b(
    indice: &[i16],
    isf_q: &mut [i16; M],
    past_isfq: &mut [i16; M],
    isfold: &[i16; M],
    isf_buf: &mut [i16],
    bfi: bool,
) {
    if !bfi {
        let ix = |k: usize| indice[k] as usize;
        for i in 0..9 {
            isf_q[i] = DICO1_ISF[ix(0) * 9 + i];
        }
        for i in 0..7 {
            isf_q[i + 9] = add_int16(DICO2_ISF[ix(1) * 7 + i], DICO23_ISF_36B[ix(4) * 7 + i]);
        }
        for i in 0..5 {
            isf_q[i] = add_int16(isf_q[i], DICO21_ISF_36B[ix(2) * 5 + i]);
        }
        for i in 0..4 {
            isf_q[i + 5] = add_int16(isf_q[i + 5], DICO22_ISF_36B[ix(3) * 4 + i]);
        }
        for i in 0..M {
            let tmp = isf_q[i];
            isf_q[i] = add_int16(tmp, MEAN_ISF[i]);
            isf_q[i] = add_int16(isf_q[i], mult_int16(MU, past_isfq[i]));
            past_isfq[i] = tmp;
        }
        push_isf_buf(isf_buf, isf_q);
    } else {
        conceal_isf(isf_q, past_isfq, isfold, isf_buf);
    }
    reorder_isf(isf_q, ISF_GAP);
}

/// `Disf_ns`: comfort-noise (SID) ISF dequantisation.
pub(crate) fn disf_ns(indice: &[i16; 5], isf_q: &mut [i16; M]) {
    let ix = |k: usize| indice[k] as usize;
    isf_q[0] = DICO1_ISF_NOISE[ix(0) * 2];
    isf_q[1] = DICO1_ISF_NOISE[ix(0) * 2 + 1];
    for i in 0..3 {
        isf_q[i + 2] = DICO2_ISF_NOISE[ix(1) * 3 + i];
        isf_q[i + 5] = DICO3_ISF_NOISE[ix(2) * 3 + i];
    }
    for i in 0..4 {
        isf_q[i + 8] = DICO4_ISF_NOISE[ix(3) * 4 + i];
        isf_q[i + 12] = DICO5_ISF_NOISE[ix(4) * 4 + i];
    }
    for i in 0..M {
        isf_q[i] = add_int16(isf_q[i], MEAN_ISF_NOISE[i]);
    }
    reorder_isf(isf_q, ISF_GAP);
}

/// `Isf_isp`: ISF (normalised frequencies) to ISP (cosine domain), `m` = `isf.len()`. In place
/// when `isf` and `isp` are one buffer (the C calls it so): pass the same values.
pub(crate) fn isf_isp(isf: &[i16], isp: &mut [i16]) {
    let m = isp.len();
    isp[..m - 1].copy_from_slice(&isf[..m - 1]);
    isp[m - 1] = shl_int16(isf[m - 1], 1);
    isf_to_cos(isp);
}

/// The second half of `Isf_isp`, in place (the cosine table lookup).
///
/// The C reads `table[ind + 1]` with `ind = isf >> 7`: past the 129-entry table at 16384,
/// which comfort noise reaches (`CN_dithering` clamps ISF 14 to 16384), there times offset 0.
/// Here `ind` ≤ 127 and `offset` ≤ 128: the same value at 16384 (`table[127]` + the full
/// step), unchanged below. Values outside [0, 16384] (undefined in the C, never seen) are
/// clamped in release builds and fail the debug assertion (tests, fuzzing).
fn isf_to_cos(isp: &mut [i16]) {
    for v in isp.iter_mut() {
        debug_assert!(
            (0..=16384).contains(v),
            "ISF/ISP {v} outside the table's domain"
        );
        let x = (*v).clamp(0, 16384);
        let ind = ((x >> 7) as usize).min(127);
        let offset = x - ((ind as i16) << 7);
        let l = mul_16by16_to_int32(TABLE[ind + 1] - TABLE[ind], offset);
        *v = add_int16(TABLE[ind], (l >> 8) as i16);
    }
}

/// `Isf_isp(HfIsf, HfIsf, M16k)`: the in-place call.
pub(crate) fn isf_isp_inplace(isf: &mut [i16]) {
    let m = isf.len();
    isf[m - 1] = shl_int16(isf[m - 1], 1);
    isf_to_cos(isf);
}

/// `Get_isp_pol`: polynomial coefficients from every second ISP of `isp`, starting at `isp[0]`.
/// `SCALE`: 9 (12.8 kHz, `f[0]` = 0x00800000) or 7 (16 kHz, 0x00200000). Saturating as
/// TS 26.173 (`L_sub`, `L_add`, `L_msu`); opencore wraps.
fn get_isp_pol<const SCALE: u32>(isp: &[i16], f: &mut [i32], n: usize) {
    let half = 1i16 << (SCALE - 1); // L_msu(f, isp, half) = f - isp << SCALE
    f[0] = 1 << (14 + SCALE);
    f[1] = mul_16by16_to_int32(isp[0], -half);
    let mut fi = 2usize;
    let mut ii = 2usize;
    for i in 2..=n {
        f[fi] = f[fi - 2];
        for _ in 1..i {
            let t0 = shl_int32(mpy_32_16(f[fi - 1], isp[ii]), 1);
            f[fi] = sub_int32(f[fi], t0);
            f[fi] = add_int32(f[fi], f[fi - 2]);
            fi -= 1;
        }
        f[fi] = msu_16by16_from_int32(f[fi], isp[ii], half);
        fi += i;
        ii += 2;
    }
}

/// `Isp_Az`: ISP to LP coefficients `a[0..=m]` (Q12), `m` = `isp.len()` (16 or 20).
/// `adaptive_scaling`: scale `a` down when the coefficients would overflow (12.8 kHz, DTX).
pub(crate) fn isp_az(isp: &[i16], a: &mut [i16], adaptive_scaling: bool) {
    let m = isp.len();
    let nc = m >> 1;
    let mut f1 = [0i32; M16K / 2 + 1];
    let mut f2 = [0i32; M16K / 2];
    if nc > 8 {
        get_isp_pol::<7>(isp, &mut f1, nc);
        for v in &mut f1[..=nc] {
            *v = shl_int32(*v, 2);
        }
        get_isp_pol::<7>(&isp[1..], &mut f2, nc - 1);
        for v in &mut f2[..nc] {
            *v = shl_int32(*v, 2);
        }
    } else {
        get_isp_pol::<9>(isp, &mut f1, nc);
        get_isp_pol::<9>(&isp[1..], &mut f2, nc - 1);
    }
    // F2(z)·(1 - z^-2), then F1·(1 + isp[m-1]) and F2·(1 - isp[m-1]), saturating as TS
    // 26.173 (opencore wraps: the comfort-noise filters differed on extreme ISFs)
    for i in (2..nc).rev() {
        f2[i] = sub_int32(f2[i], f2[i - 2]);
    }
    let last = isp[m - 1];
    for i in 0..nc {
        f1[i] = add_int32(f1[i], mpy_32_16(f1[i], last));
        f2[i] = sub_int32(f2[i], mpy_32_16(f2[i], last));
    }
    a[0] = 4096;
    let mut tmax = 1i32;
    let rnd = |t: i32, q: i32| {
        (t.wrapping_shr(q as u32)).wrapping_add(t.wrapping_shr((q - 1) as u32) & 1) as i16
    };
    let mut j = m - 1;
    for i in 1..nc {
        let t0 = add_int32(f1[i], f2[i]);
        let t1 = t0.wrapping_sub((t0 < 0) as i32);
        tmax |= t1 ^ (t1 >> 31);
        a[i] = rnd(t0, 12);
        let t0 = sub_int32(f1[i], f2[i]);
        let t1 = t0.wrapping_sub((t0 < 0) as i32);
        tmax |= t1 ^ (t1 >> 31);
        a[j] = rnd(t0, 12);
        j -= 1;
    }
    let mut q = if adaptive_scaling {
        4 - normalize_amr_wb(tmax)
    } else {
        0
    };
    let q_sug;
    if q > 0 {
        q_sug = 12 + q as i32;
        let mut j = m - 1;
        for i in 1..nc {
            a[i] = rnd(add_int32(f1[i], f2[i]), q_sug);
            a[j] = rnd(sub_int32(f1[i], f2[i]), q_sug);
            j -= 1;
        }
        a[0] >>= q;
    } else {
        q_sug = 12;
        q = 0;
    }
    let t0 = add_int32(f1[nc], mpy_32_16(f1[nc], last));
    a[nc] = rnd(t0, q_sug);
    a[m] = shr_rnd(isp[m - 1], 3 + q);
}

/// `interpolate_isp`: A(z) of subframes 0..2 from ISPs interpolated with `frac`, and of
/// subframe 3 from `isp_new`. `az` holds 4 × (M + 1).
pub(crate) fn interpolate_isp(
    isp_old: &[i16; M],
    isp_new: &[i16; M],
    frac: &[i16],
    az: &mut [i16],
) {
    let mut isp = [0i16; M];
    for k in 0..3 {
        let fac_new = frac[k];
        let fac_old = add_int16(sub_int16(32767, fac_new), 1);
        for i in 0..M {
            let l = mul_16by16_to_int32(isp_old[i], fac_old);
            isp[i] = amr_wb_round(mac_16by16_to_int32(l, isp_new[i], fac_new));
        }
        isp_az(&isp, &mut az[k * (M + 1)..], false);
    }
    isp_az(isp_new, &mut az[3 * (M + 1)..], false);
}

/// `weight_amrwb_lpc`: `ap[i] = a[i]·gamma^i`, `m` = `a.len() - 1`, rounded as TS 26.173's
/// `round(L_mult(·))` (opencore's `(x·y + 0x4000) >> 15` wraps for -32768²).
pub(crate) fn weight_amrwb_lpc(a: &[i16], ap: &mut [i16], gamma: i16) {
    let m = a.len() - 1;
    ap[0] = a[0];
    let mut fac = gamma;
    for i in 1..m {
        ap[i] = amr_wb_round(mul_16by16_to_int32(a[i], fac));
        fac = amr_wb_round(mul_16by16_to_int32(fac, gamma));
    }
    ap[m] = amr_wb_round(mul_16by16_to_int32(a[m], fac));
}

const INV_LENGTH: i16 = 2731;

/// `isf_extrapolation`: extends the 16 ISFs of `hf_isf` to 20 (6.60 kbit/s high band) and
/// converts them to ISPs in place.
pub(crate) fn isf_extrapolation(hf_isf: &mut [i16; M16K]) {
    let mut isf_diff = [0i16; M - 2];
    let mut isf_corr = [0i32; 3];
    hf_isf[M16K - 1] = hf_isf[M - 1];
    for i in 1..M - 1 {
        isf_diff[i - 1] = sub_int16(hf_isf[i], hf_isf[i - 1]);
    }
    let mut l = 0i32;
    for i in 3..M - 1 {
        l = mac_16by16_to_int32(l, isf_diff[i - 1], INV_LENGTH);
    }
    let mut mean = amr_wb_round(l);
    let mut tmp = 0i16;
    for &d in &isf_diff {
        if d > tmp {
            tmp = d;
        }
    }
    let exp = norm_s(tmp);
    for d in &mut isf_diff {
        *d = shl_int16(*d, exp);
    }
    mean = shl_int16(mean, exp);
    for (k, lag) in [2usize, 3, 4].into_iter().enumerate() {
        for i in 7..M - 2 {
            let tmp2 = sub_int16(isf_diff[i], mean);
            let tmp3 = sub_int16(isf_diff[i - lag], mean);
            let (hi, lo) = int32_to_dpf(mul_16by16_to_int32(tmp2, tmp3));
            isf_corr[k] = add_int32(isf_corr[k], mpy_dpf_32(hi, lo, hi, lo));
        }
    }
    let mut max_corr = if isf_corr[0] > isf_corr[1] { 0 } else { 1 };
    if isf_corr[2] > isf_corr[max_corr] {
        max_corr = 2;
    }
    max_corr += 1;
    for i in M - 1..M16K - 1 {
        let t = sub_int16(hf_isf[i - 1 - max_corr], hf_isf[i - 2 - max_corr]);
        hf_isf[i] = add_int16(hf_isf[i - 1], t);
    }
    let mut tmp = add_int16(hf_isf[4], hf_isf[3]);
    tmp = sub_int16(hf_isf[2], tmp);
    tmp = mult_int16(tmp, 5461);
    tmp = tmp.wrapping_add(20390);
    if tmp > 19456 {
        tmp = 19456;
    }
    tmp = sub_int16(tmp, hf_isf[M - 2]);
    let mut tmp2 = sub_int16(hf_isf[M16K - 2], hf_isf[M - 2]);
    let exp2 = norm_s(tmp2);
    let exp = norm_s(tmp) - 1;
    // `int16 <<= n` in C: promoted, shifted (x86 masks the count to 5 bits), truncated.
    tmp = (tmp as i32).wrapping_shl(exp as u32) as i16;
    tmp2 = (tmp2 as i32).wrapping_shl(exp2 as u32) as i16;
    let coeff = div_16by16(tmp, tmp2);
    let exp = exp2 - exp;
    let mut diff = [0i16; M16K - M];
    for i in M - 1..M16K - 1 {
        let t = mult_int16(sub_int16(hf_isf[i], hf_isf[i - 1]), coeff);
        diff[i - (M - 1)] = shl_int16(t, exp);
    }
    for i in M..M16K - 1 {
        let (a, b) = (diff[i - (M - 1)], diff[i - M]);
        let t = a.wrapping_add(b).wrapping_sub(1280);
        if t < 0 {
            if a > b {
                diff[i - M] = 1280i16.wrapping_sub(a);
            } else {
                diff[i - (M - 1)] = 1280i16.wrapping_sub(b);
            }
        }
    }
    for i in M - 1..M16K - 1 {
        hf_isf[i] = add_int16(hf_isf[i - 1], diff[i - (M - 1)]);
    }
    for v in &mut hf_isf[..M16K - 1] {
        *v = mult_int16(*v, 26214);
    }
    isf_isp_inplace(hf_isf);
}
