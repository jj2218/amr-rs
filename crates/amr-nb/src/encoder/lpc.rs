//! Encoder LP analysis and LSP quantisation: `Autocorr`, `Lag_window`, `Levinson`, `lpc`,
//! `Az_lsp`, `Lsf_wt`, `Q_plsf_3`, `Q_plsf_5` and the `lsp` frame driver.
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `enc/src/{autocorr,lag_wind,levinson,
//! lpc,div_32}.cpp` and `common/src/{az_lsp,lsfwt,q_plsf,q_plsf_3,q_plsf_5,lsp}.cpp`.
#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1
#![allow(clippy::too_many_arguments)] // signatures mirror the C functions

pub(crate) use crate::simd::{vq_pair, vq_pair_s, vq_subvec3, vq_subvec4};

use super::layout::{Soa, cols, first_min, soa};
use crate::basic_op::{
    Flag, MAX_32, MIN_16, abs_s, add, div_s, l_extract, l_mac, l_msu, l_mult, l_shl, l_shr, mpy_32,
    mpy_32_16, negate, norm_l, norm_s, round, shl, shr, sub,
};
use crate::cnst::{CodecMode, M, MP1};
use crate::lpc::{
    int_lpc_1and3, int_lpc_1and3_2, int_lpc_1to3, int_lpc_1to3_2, lsf_lsp, lsp_lsf, reorder_lsf,
};
use crate::tables::{
    DICO1_LSF_3, DICO1_LSF_5, DICO2_LSF_3, DICO2_LSF_5, DICO3_LSF_3, DICO3_LSF_5, DICO4_LSF_5,
    DICO5_LSF_5, GRID, LAG_H, LAG_L, LSP_INIT_DATA, MEAN_LSF_3, MEAN_LSF_5, MR515_3_LSF,
    MR795_1_LSF, PAST_RQ_INIT, PRED_FAC_3, WINDOW_160_80, WINDOW_200_40, WINDOW_232_8,
};

/// LP analysis window length.
pub(crate) const L_WINDOW: usize = 240;
const LSF_GAP: i16 = 205;
const LSP_PRED_FAC_MR122: i32 = 21299;
const GRID_POINTS: usize = 60;
const NC: usize = M / 2;
const PAST_RQ_INIT_SIZE: usize = 8;

/// `(hi, lo)` double-precision split of a 32-bit value (the C's inline `L_Extract`).
#[inline(always)]
pub(super) fn split(x: i32) -> (i16, i16) {
    let hi = (x >> 16) as i16;
    (hi, ((x >> 1) - ((hi as i32) << 15)) as i16)
}

/// `L_abs` (opencore's branch-free variant: `L_abs(MIN_32) == MAX_32`).
#[inline(always)]
pub(crate) fn l_abs(x: i32) -> i32 {
    let y = x.wrapping_sub((x < 0) as i32);
    y ^ (y >> 31)
}

/// `Div_32`: `L_num / (denom_hi, denom_lo)`, both positive, `L_num < denom`.
pub(crate) fn div_32(l_num: i32, denom_hi: i16, denom_lo: i16, ov: &mut Flag) -> i32 {
    let approx = div_s(0x3fff, denom_hi);
    let mut result = mpy_32_16(denom_hi, denom_lo, approx, ov);
    result = 0x7fff_ffff - result;
    let (hi, lo) = split(result);
    result = mpy_32_16(hi, lo, approx, ov);
    let (hi, lo) = split(result);
    let (n_hi, n_lo) = split(l_num);
    result = mpy_32(n_hi, n_lo, hi, lo, ov);
    l_shl(result, 2, ov)
}

/// `Autocorr`: windowed autocorrelation `r[0..=m]` as `(r_h, r_l)`; returns the norm.
pub(crate) fn autocorr(
    x: &[i16],
    m: usize,
    r_h: &mut [i16],
    r_l: &mut [i16],
    wind: &[i16; L_WINDOW],
) -> i16 {
    let mut y: [i16; L_WINDOW] =
        std::array::from_fn(|i| ((x[i] as i32 * wind[i] as i32).wrapping_add(0x4000) >> 15) as i16);
    // TS 26.073 sums 2*t^2 with `L_mac` (saturating) and divides y by 4 while the sum is
    // MAX_32. The terms are >= 0 (-32768^2 saturates alone), so that is: the exact sum >= 2^31
    // (even sums never equal MAX_32), one vectorisable sum per round. (opencore re-sums with
    // wrapping after a shift and tests for <= 0: wrong when a clipped window still overflows.)
    let energy = |y: &[i16; L_WINDOW]| -> u64 {
        y.iter()
            .map(|&t| (t as i32 * t as i32) as u32 as u64)
            .sum::<u64>()
            << 1
    };
    let mut exact = energy(&y);
    let mut overfl_shft: i16 = 0;
    while exact >= 1 << 31 {
        overfl_shft += 4;
        y.iter_mut().for_each(|v| *v >>= 2);
        exact = energy(&y);
    }
    let mut sum = (exact as i32).wrapping_add(1);
    let norm = norm_l(sum);
    sum <<= norm;
    (r_h[0], r_l[0]) = split(sum);
    // y zero-padded: every lag sums over the full window length (no remainder loops).
    assert!(m <= 16, "autocorr: LP order {m} > the padding");
    let mut yp = [0i16; L_WINDOW + 16];
    yp[..L_WINDOW].copy_from_slice(&y);
    // lags 1..=m four at a time (`simd::dots4`)
    for g in (1..=m).step_by(4) {
        let d = crate::simd::dots4::<L_WINDOW>(&y, &yp, g);
        for (i, s) in (g..=m).zip(d) {
            (r_h[i], r_l[i]) = split(s.wrapping_shl((norm + 1) as u32));
        }
    }
    norm - overfl_shft
}

/// `Lag_window`: applies the lag window to `r[1..=m]`.
pub(crate) fn lag_window(m: usize, r_h: &mut [i16], r_l: &mut [i16], ov: &mut Flag) {
    for i in 1..=m {
        let x = mpy_32(r_h[i], r_l[i], LAG_H[i - 1], LAG_L[i - 1], ov);
        (r_h[i], r_l[i]) = split(x);
    }
}

/// `LevinsonState`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct LevinsonState {
    pub old_a: [i16; MP1],
}

impl Default for LevinsonState {
    fn default() -> Self {
        let mut old_a = [0; MP1];
        old_a[0] = 4096;
        Self { old_a }
    }
}

/// `Levinson`: LPC coefficients `a[0..=M]` (Q12) and the first 4 reflection coefficients.
pub(crate) fn levinson(
    st: &mut LevinsonState,
    rh: &[i16],
    rl: &[i16],
    a: &mut [i16],
    rc: &mut [i16; 4],
    ov: &mut Flag,
) {
    let mut ah = [0i16; M + 1];
    let mut al = [0i16; M + 1];
    let mut anh = [0i16; M + 1];
    let mut anl = [0i16; M + 1];

    let t1 = ((rh[1] as i32) << 16).wrapping_add((rl[1] as i32) << 1);
    let t2 = l_abs(t1);
    let mut t0 = div_32(t2, rh[0], rl[0], ov);
    if t1 > 0 {
        t0 = -t0; // L_negate: div_32's result is never MIN_32
    }
    let (mut kh, mut kl) = split(t0);
    rc[0] = round(t0, ov);
    t0 >>= 4;
    (ah[1], al[1]) = split(t0);
    t0 = l_abs(mpy_32(kh, kl, kh, kl, ov));
    t0 = 0x7fff_ffff - t0;
    let (hi, lo) = split(t0);
    t0 = mpy_32(rh[0], rl[0], hi, lo, ov);
    let mut alp_exp = norm_l(t0);
    t0 <<= alp_exp;
    let (mut alp_h, mut alp_l) = split(t0);

    for i in 2..=M {
        let mut t0: i32 = 0;
        for j in 1..i {
            t0 = t0.wrapping_add((rh[j] as i32 * al[i - j] as i32) >> 15);
            t0 = t0.wrapping_add((rl[j] as i32 * ah[i - j] as i32) >> 15);
            t0 = t0.wrapping_add(rh[j] as i32 * ah[i - j] as i32);
        }
        t0 = t0.wrapping_shl(5);
        let t1 = ((rh[i] as i32) << 16).wrapping_add((rl[i] as i32) << 1);
        t0 = t0.wrapping_add(t1);
        let t1 = l_abs(t0);
        let mut t2 = div_32(t1, alp_h, alp_l, ov);
        if t0 > 0 {
            t2 = l_negate(t2);
        }
        t2 = l_shl(t2, alp_exp, ov);
        (kh, kl) = split(t2);
        if i < 5 {
            rc[i - 1] = (t2.wrapping_add(0x8000) >> 16) as i16;
        }
        if abs_s(kh) > 32750 {
            a[..MP1].copy_from_slice(&st.old_a);
            *rc = [0; 4];
            return;
        }
        for j in 1..i {
            let mut t0 = (kh as i32 * al[i - j] as i32) >> 15;
            t0 = t0.wrapping_add((kl as i32 * ah[i - j] as i32) >> 15);
            t0 = t0.wrapping_add(kh as i32 * ah[i - j] as i32);
            t0 = t0.wrapping_add(((ah[j] as i32) << 15).wrapping_add(al[j] as i32));
            anh[j] = (t0 >> 15) as i16;
            anl[j] = t0.wrapping_sub((anh[j] as i32) << 15) as i16;
        }
        anh[i] = (t2 >> 20) as i16;
        anl[i] = ((t2 >> 5) - ((anh[i] as i32) << 15)) as i16;

        let mut t0 = l_abs(mpy_32(kh, kl, kh, kl, ov));
        t0 = 0x7fff_ffff - t0;
        let (hi, lo) = split(t0);
        t0 = (alp_h as i32 * lo as i32) >> 15;
        t0 = t0.wrapping_add((alp_l as i32 * hi as i32) >> 15);
        t0 = t0.wrapping_add(alp_h as i32 * hi as i32);
        t0 = t0.wrapping_shl(1);
        let j = norm_l(t0);
        t0 = t0.wrapping_shl(j as u32);
        (alp_h, alp_l) = split(t0);
        alp_exp += j;
        ah[1..=i].copy_from_slice(&anh[1..=i]);
        al[1..=i].copy_from_slice(&anl[1..=i]);
    }
    a[0] = 4096;
    for i in 1..=M {
        let t0 = ((ah[i] as i32) << 15).wrapping_add(al[i] as i32);
        a[i] = (t0.wrapping_add(0x2000) >> 14) as i16;
        st.old_a[i] = a[i];
    }
}

/// `L_negate`.
#[inline(always)]
pub(super) fn l_negate(x: i32) -> i32 {
    if x == i32::MIN { MAX_32 } else { -x }
}

/// `lpc`: LP analysis of one frame into `a` (4 × `MP1`; only the subframes the mode analyses
/// are written: 1 and 3 for MR122, 3 otherwise). `x` is the frame's analysis window start
/// (`p_window`) and `x_12k2` the MR122 one (`p_window_12k2`), both `L_WINDOW` long.
pub(crate) fn lpc(
    st: &mut LevinsonState,
    mode: CodecMode,
    x: &[i16],
    x_12k2: &[i16],
    a: &mut [i16; 4 * MP1],
    ov: &mut Flag,
) {
    let mut rc = [0i16; 4];
    let mut r_low = [0i16; MP1];
    let mut r_high = [0i16; MP1];
    if mode == CodecMode::MR122 {
        autocorr(x_12k2, M, &mut r_high, &mut r_low, &WINDOW_160_80);
        lag_window(M, &mut r_high, &mut r_low, ov);
        levinson(st, &r_high, &r_low, &mut a[MP1..2 * MP1], &mut rc, ov);
        autocorr(x_12k2, M, &mut r_high, &mut r_low, &WINDOW_232_8);
        lag_window(M, &mut r_high, &mut r_low, ov);
        levinson(st, &r_high, &r_low, &mut a[3 * MP1..], &mut rc, ov);
    } else {
        autocorr(x, M, &mut r_high, &mut r_low, &WINDOW_200_40);
        lag_window(M, &mut r_high, &mut r_low, ov);
        levinson(st, &r_high, &r_low, &mut a[3 * MP1..], &mut rc, ov);
    }
}

/// `Chebps`: the Chebyshev polynomial `f` (order `n`) at `x`, as TS 26.073's operator chain.
/// Fast form: the C keeps b1 as an `L_Extract` pair (hi, lo) and forms `Mpy_32_16` =
/// 2 (hi x + (lo x >> 15)) = (b x >> 16) << 1 (hi x 2^15 is a multiple of 2^15; b is even), and
/// `L_mac(t, b2_h, 0x8000)` + `L_msu(t, b2_l, 1)` subtract b2: one multiply per step. Exact while
/// no partial sum reaches 2^31 (checked per step, in i64); the final `extract_h(L_shl(t, 6))`
/// is then `clamp(t >> 10)`. Otherwise the operator chain (`chebps_exact`).
pub(super) fn chebps(x: i16, f: &[i16; NC + 1], n: usize) -> i16 {
    // Cheap bound first: |m| <= 2|b| + 4, |c| <= 2^29, and the C subtracts b2 as hi << 16
    // (|hi << 16| < |b2| + 2^16) then lo << 1: if every b stays below 2^29 - 2^17, no partial
    // sum reaches 2^31 at any step. Tracked as the OR of one's-complement |b| (an
    // upper bound of their maximum), off the serial chain; c - b2 is formed off it too.
    let xx = x as i64;
    let mut b2: i64 = 1 << 24;
    let mut b: i64 = ((x as i64) << 10) + ((f[1] as i64) << 14);
    let mut acc = b ^ (b >> 63);
    for &fi in &f[2..n] {
        let c = ((fi as i64) << 14) - b2;
        b2 = b;
        b = (((b * xx) >> 16) << 2) + c;
        acc |= b ^ (b >> 63);
    }
    let t = (((b * xx) >> 16) << 1) + ((f[n] as i64) << 13) - b2;
    if acc < (1 << 29) - (1 << 17) {
        return (t >> 10).clamp(i16::MIN as i64, i16::MAX as i64) as i16;
    }
    chebps_checked(x, f, n)
}

/// `chebps` with the bound checked per step (exact while no partial sum reaches 2^31).
#[cold]
fn chebps_checked(x: i16, f: &[i16; NC + 1], n: usize) -> i16 {
    // The C subtracts b2 as hi << 16 first: |hi << 16| < |b2| + 2^16.
    const LIM: i64 = (1 << 31) - (1 << 16);
    let xx = x as i64;
    let mut b2: i64 = 1 << 24;
    let mut b: i64 = ((x as i64) << 10) + ((f[1] as i64) << 14);
    let mut ok = true;
    for &fi in &f[2..n] {
        let m = ((b * xx) >> 16) << 2; // L_shl(Mpy_32_16, 1)
        let c = (fi as i64) << 14;
        ok &= m.abs() + b2.abs() + c.abs() < LIM;
        (b2, b) = (b, m - b2 + c);
    }
    let m = ((b * xx) >> 16) << 1;
    let c = (f[n] as i64) << 13;
    ok &= m.abs() + b2.abs() + c.abs() < LIM;
    if !ok {
        return chebps_exact(x, f, n);
    }
    ((m - b2 + c) >> 10).clamp(i16::MIN as i64, i16::MAX as i64) as i16
}

/// `Chebps` as TS 26.073's operator chain (saturating; opencore wraps): `chebps`' fallback.
#[cold]
pub(super) fn chebps_exact(x: i16, f: &[i16; NC + 1], n: usize) -> i16 {
    let ov = &mut false;
    let (mut b2_h, mut b2_l) = (256i16, 0i16);
    let mut t0 = l_mac(l_mult(x, 512, ov), f[1], 8192, ov);
    let (mut b1_h, mut b1_l) = l_extract(t0, ov);
    for &fi in &f[2..n] {
        t0 = l_shl(mpy_32_16(b1_h, b1_l, x, ov), 1, ov);
        t0 = l_mac(t0, b2_h, MIN_16, ov);
        t0 = l_msu(t0, b2_l, 1, ov);
        t0 = l_mac(t0, fi, 8192, ov);
        let (b0_h, b0_l) = l_extract(t0, ov);
        (b2_h, b2_l, b1_h, b1_l) = (b1_h, b1_l, b0_h, b0_l);
    }
    t0 = mpy_32_16(b1_h, b1_l, x, ov);
    t0 = l_mac(t0, b2_h, MIN_16, ov);
    t0 = l_msu(t0, b2_l, 1, ov);
    t0 = l_mac(t0, f[n], 4096, ov);
    (l_shl(t0, 6, ov) >> 16) as i16
}

/// `Az_lsp`: LPC coefficients `a[0..=M]` to LSPs; keeps `old_lsp` if fewer than M roots are
/// found. TS 26.073's operators (saturating; opencore wraps in the polynomials and the root
/// interpolation).
pub(crate) fn az_lsp(a: &[i16], lsp: &mut [i16], old_lsp: &[i16]) {
    let ov = &mut false;
    let mut f1 = [0i16; NC + 1];
    let mut f2 = [0i16; NC + 1];
    f1[0] = 1024;
    f2[0] = 1024;
    for i in 0..NC {
        let x = (l_mac(l_mult(a[i + 1], 8192, ov), a[M - i], 8192, ov) >> 16) as i16;
        f1[i + 1] = sub(x, f1[i], ov);
        let y = (l_msu(l_mult(a[i + 1], 8192, ov), a[M - i], 8192, ov) >> 16) as i16;
        f2[i + 1] = add(y, f2[i], ov);
    }
    let mut nf = 0usize;
    let mut ip = false;
    let mut coef = &f1;
    let mut xlow = GRID[0];
    let mut ylow = chebps(xlow, coef, NC);
    let mut j = 0usize;
    while nf < M && j < GRID_POINTS {
        j += 1;
        let mut xhigh = xlow;
        let mut yhigh = ylow;
        xlow = GRID[j];
        ylow = chebps(xlow, coef, NC);
        if l_mult(ylow, yhigh, ov) <= 0 {
            for _ in 0..4 {
                let xmid = add(shr(xlow, 1, ov), shr(xhigh, 1, ov), ov);
                let ymid = chebps(xmid, coef, NC);
                if l_mult(ylow, ymid, ov) <= 0 {
                    yhigh = ymid;
                    xhigh = xmid;
                } else {
                    ylow = ymid;
                    xlow = xmid;
                }
            }
            let x = sub(xhigh, xlow, ov);
            let mut y = sub(yhigh, ylow, ov);
            let xint = if y == 0 {
                xlow
            } else {
                let sign = y;
                y = abs_s(y);
                let exp = norm_s(y);
                y = shl(y, exp, ov);
                y = div_s(16383, y);
                y = l_shr(l_mult(x, y, ov), sub(20, exp, ov), ov) as i16;
                if sign < 0 {
                    y = negate(y);
                }
                sub(xlow, l_shr(l_mult(ylow, y, ov), 11, ov) as i16, ov)
            };
            lsp[nf] = xint;
            xlow = xint;
            nf += 1;
            ip = !ip;
            coef = if ip { &f2 } else { &f1 };
            ylow = chebps(xlow, coef, NC);
        }
    }
    if nf < M {
        lsp[..M].copy_from_slice(&old_lsp[..M]);
    }
}

/// `Lsf_wt`: LSF weighting factors for the quantiser.
pub(crate) fn lsf_wt(lsf: &[i16; M], wf: &mut [i16; M]) {
    wf[0] = lsf[1];
    for i in 1..M - 1 {
        wf[i] = lsf[i + 1].wrapping_sub(lsf[i - 1]);
    }
    wf[M - 1] = 16384i16.wrapping_sub(lsf[M - 2]);
    for w in wf.iter_mut() {
        let wgt_fct = *w;
        let temp = wgt_fct.wrapping_sub(1843);
        let wgt_fct = if temp > 0 {
            1843i16.wrapping_sub(((temp as i32 * 6242) >> 15) as i16)
        } else {
            3427i16.wrapping_sub(((wgt_fct as i32 * 28160) >> 15) as i16)
        };
        *w = ((wgt_fct as i32) << 3) as i16;
    }
}

pub(crate) static DICO1_SOA: Soa<3, 256> = soa(&DICO1_LSF_3, 3);
pub(crate) static DICO2_SOA: Soa<3, 512> = soa(&DICO2_LSF_3, 3);
/// `Vq_subvec3`'s `use_half`: every other codeword of `DICO2_LSF_3`.
pub(crate) static DICO2_HALF_SOA: Soa<3, 256> = soa(&DICO2_LSF_3, 6);
pub(crate) static MR795_1_SOA: Soa<3, 512> = soa(&MR795_1_LSF, 3);
pub(crate) static DICO3_SOA: Soa<4, 512> = soa(&DICO3_LSF_3, 4);
pub(crate) static MR515_3_SOA: Soa<4, 128> = soa(&MR515_3_LSF, 4);

/// `Vq_subvec3` / `Vq_subvec4` (`D` = 3, 4; `use_half` is a [`Soa`] of its own): the weighted
/// distance of `VQ_W` consecutive codewords in lanes (32: the AVX2 copy; 16: other targets).
/// Wrapping i32 sums, as the C; a strict `<` minimum per lane, then [`first_min`]: the C's
/// first minimum. Overwrites `lsf_r1[..D]` with the codeword. (x86 without AVX2 runs
/// `simd::vq_subvec_sse2`, i16 lanes: SSE2 has no `pmulld`.)
#[inline(always)]
#[cfg_attr(
    all(target_arch = "x86_64", portable, not(amr_nb_arrays)),
    allow(dead_code) // that build runs `simd::vq_subvec_sse2` only
)]
pub(crate) fn vq_subvec_portable<const D: usize, const VQ_W: usize>(
    lsf_r1: &mut [i16],
    cols: &[&[i16]; D],
    wf1: &[i16],
) -> i16 {
    let n = cols[0].len();
    assert!(n.is_multiple_of(VQ_W) && cols.iter().all(|c| c.len() == n));
    let r: [i16; D] = std::array::from_fn(|k| lsf_r1[k]);
    let w: [i32; D] = std::array::from_fn(|k| wf1[k] as i32);
    let mut best = [MAX_32; VQ_W];
    let mut bidx = [0i32; VQ_W];
    for b in (0..n).step_by(VQ_W) {
        let mut dist = [0i32; VQ_W];
        for k in 0..D {
            let c: &[i16; VQ_W] = cols[k][b..b + VQ_W].try_into().unwrap();
            for l in 0..VQ_W {
                let t = ((w[k] * r[k].wrapping_sub(c[l]) as i32) >> 15) as i16 as i32;
                dist[l] = dist[l].wrapping_add(t * t);
            }
        }
        for l in 0..VQ_W {
            let lt = dist[l] < best[l];
            best[l] = if lt { dist[l] } else { best[l] };
            bidx[l] = if lt { (b + l) as i32 } else { bidx[l] };
        }
    }
    let index = first_min(&best, &bidx);
    for k in 0..D {
        lsf_r1[k] = cols[k][index as usize];
    }
    index as i16
}

/// `Q_plsfState`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct QPlsfState {
    pub past_rq: [i16; M],
}

/// `Q_plsf_3`: quantises one LSP vector (all modes but MR122; `MRDTX` searches the predictor
/// initialisation and returns its index in `pred_init_i`).
pub(crate) fn q_plsf_3(
    st: &mut QPlsfState,
    mode: CodecMode,
    lsp1: &[i16],
    lsp1_q: &mut [i16],
    indice: &mut [i16],
    pred_init_i: &mut i16,
    ov: &mut Flag,
) {
    let mut lsf1 = [0i16; M];
    let mut wf1 = [0i16; M];
    let mut lsf_p = [0i16; M];
    let mut lsf_r1 = [0i16; M];
    let mut lsf1_q = [0i16; M];
    lsp_lsf(lsp1, &mut lsf1, M, ov);
    lsf_wt(&lsf1, &mut wf1);
    if mode != CodecMode::MRDTX {
        for i in 0..M {
            let temp = ((st.past_rq[i] as i32 * PRED_FAC_3[i] as i32) >> 15) as i16;
            lsf_p[i] = MEAN_LSF_3[i].wrapping_add(temp);
            lsf_r1[i] = lsf1[i].wrapping_sub(lsf_p[i]);
        }
    } else {
        *pred_init_i = 0;
        let mut l_min = 0x7fff_ffffi32;
        let mut temp_p = [0i16; M];
        let mut temp_r1 = [0i16; M];
        for j in 0..PAST_RQ_INIT_SIZE {
            let mut err: i32 = 0;
            for i in 0..M {
                temp_p[i] = MEAN_LSF_3[i].wrapping_add(PAST_RQ_INIT[j * M + i]);
                temp_r1[i] = lsf1[i].wrapping_sub(temp_p[i]);
                let l_temp = temp_r1[i] as i32 * temp_r1[i] as i32;
                err = err.wrapping_add(l_temp.wrapping_shl(1));
            }
            if err < l_min {
                l_min = err;
                lsf_r1 = temp_r1;
                lsf_p = temp_p;
                st.past_rq.copy_from_slice(&PAST_RQ_INIT[j * M..j * M + M]);
                *pred_init_i = j as i16;
            }
        }
    }
    let (d1, d2, d3) = match mode {
        CodecMode::MR475 | CodecMode::MR515 => {
            (cols(&DICO1_SOA), cols(&DICO2_HALF_SOA), cols(&MR515_3_SOA))
        }
        CodecMode::MR795 => (cols(&MR795_1_SOA), cols(&DICO2_SOA), cols(&DICO3_SOA)),
        _ => (cols(&DICO1_SOA), cols(&DICO2_SOA), cols(&DICO3_SOA)),
    };
    indice[0] = vq_subvec3(&mut lsf_r1[0..], &d1, &wf1[0..]);
    indice[1] = vq_subvec3(&mut lsf_r1[3..], &d2, &wf1[3..]);
    indice[2] = vq_subvec4(&mut lsf_r1[6..], &d3, &wf1[6..]);
    for i in 0..M {
        lsf1_q[i] = lsf_r1[i].wrapping_add(lsf_p[i]);
        st.past_rq[i] = lsf_r1[i];
    }
    reorder_lsf(&mut lsf1_q, LSF_GAP, M, ov);
    lsf_lsp(&lsf1_q, lsp1_q, M, ov);
}

/// `Vq_subvec`: joint 2+2-dimensional VQ of two LSF residual pairs (with early exits).
pub(super) fn vq_subvec(
    lsf_r1: &mut [i16],
    lsf_r2: &mut [i16],
    dico: &[i16],
    wf1: &[i16],
    wf2: &[i16],
    dico_size: usize,
) -> i16 {
    let mut dist_min = MAX_32;
    let mut index = 0usize;
    let aux1 = lsf_r1[0] as i32 * wf1[0] as i32;
    let aux2 = lsf_r1[1] as i32 * wf1[1] as i32;
    let aux3 = lsf_r2[0] as i32 * wf2[0] as i32;
    let aux4 = lsf_r2[1] as i32 * wf2[1] as i32;
    let term = |aux: i32, w: i16, d: i16| -> i32 {
        let t = (aux.wrapping_sub(w as i32 * d as i32) >> 15) as i16;
        t as i32 * t as i32
    };
    for i in 0..dico_size {
        let d = &dico[4 * i..4 * i + 4];
        let mut dist = term(aux1, wf1[0], d[0]);
        if dist >= dist_min {
            continue;
        }
        dist = dist.wrapping_add(term(aux2, wf1[1], d[1]));
        if dist >= dist_min {
            continue;
        }
        dist = dist.wrapping_add(term(aux3, wf2[0], d[2]));
        if dist >= dist_min {
            continue;
        }
        dist = dist.wrapping_add(term(aux4, wf2[1], d[3]));
        if dist < dist_min {
            dist_min = dist;
            index = i;
        }
    }
    let d = &dico[4 * index..];
    (lsf_r1[0], lsf_r1[1], lsf_r2[0], lsf_r2[1]) = (d[0], d[1], d[2], d[3]);
    index as i16
}

/// `Vq_subvec_s`: as `vq_subvec`, but also tries the negated codeword (sign in index bit 0).
pub(super) fn vq_subvec_s(
    lsf_r1: &mut [i16],
    lsf_r2: &mut [i16],
    dico: &[i16],
    wf1: &[i16],
    wf2: &[i16],
    dico_size: usize,
) -> i16 {
    let mut dist_min = MAX_32;
    let mut index = 0usize;
    let mut sign = false;
    let terms = |r: i16, w: i16, d: i16| -> (i32, i32) {
        let t1 = ((w as i32 * r.wrapping_sub(d) as i32) >> 15) as i16;
        let t2 = ((w as i32 * r.wrapping_add(d) as i32) >> 15) as i16;
        (t1 as i32 * t1 as i32, t2 as i32 * t2 as i32)
    };
    for i in 0..dico_size {
        let d = &dico[4 * i..4 * i + 4];
        let (a1, a2) = terms(lsf_r1[0], wf1[0], d[0]);
        let (b1, b2) = terms(lsf_r1[1], wf1[1], d[1]);
        let mut dist1 = a1.wrapping_add(b1);
        let mut dist2 = a2.wrapping_add(b2);
        if dist1 >= dist_min && dist2 >= dist_min {
            continue;
        }
        let (c1, c2) = terms(lsf_r2[0], wf2[0], d[2]);
        let (e1, e2) = terms(lsf_r2[1], wf2[1], d[3]);
        dist1 = dist1.wrapping_add(c1).wrapping_add(e1);
        dist2 = dist2.wrapping_add(c2).wrapping_add(e2);
        if dist1 < dist_min {
            dist_min = dist1;
            index = i;
            sign = false;
        }
        if dist2 < dist_min {
            dist_min = dist2;
            index = i;
            sign = true;
        }
    }
    let d = &dico[4 * index..];
    if sign {
        (lsf_r1[0], lsf_r1[1], lsf_r2[0], lsf_r2[1]) = (
            d[0].wrapping_neg(),
            d[1].wrapping_neg(),
            d[2].wrapping_neg(),
            d[3].wrapping_neg(),
        );
        (index * 2 + 1) as i16
    } else {
        (lsf_r1[0], lsf_r1[1], lsf_r2[0], lsf_r2[1]) = (d[0], d[1], d[2], d[3]);
        (index * 2) as i16
    }
}

static DICO1_5_SOA: Soa<4, 128> = soa(&DICO1_LSF_5, 4);
static DICO2_5_SOA: Soa<4, 256> = soa(&DICO2_LSF_5, 4);
static DICO3_5_SOA: Soa<4, 256> = soa(&DICO3_LSF_5, 4);
static DICO4_5_SOA: Soa<4, 256> = soa(&DICO4_LSF_5, 4);
static DICO5_5_SOA: Soa<4, 64> = soa(&DICO5_LSF_5, 4);

/// `vq_subvec` (`SIGNED` false) or `vq_subvec_s` (true) with `VQ_W` codewords in lanes; the
/// codeword is not written. The C's early exits skip a row once a partial distance reaches the
/// running minimum; while no sum wraps, the partial sums only grow, so that is "the full distance
/// is not below it" and the C's pick is the first minimum of the full distances (over 2i + sign
/// for `vq_subvec_s`): a strict `<` minimum per lane, then the first over the lanes. Every
/// term below 2^29 means no sum of four wraps; `None` if one is not (the caller runs the C's loop).
#[inline(always)]
pub(crate) fn vq_pair_portable<const SIGNED: bool, const VQ_W: usize>(
    lsf_r1: &[i16],
    lsf_r2: &[i16],
    cols: &[&[i16]; 4],
    wf1: &[i16],
    wf2: &[i16],
) -> Option<usize> {
    let n = cols[0].len();
    assert!(n.is_multiple_of(VQ_W) && cols.iter().all(|c| c.len() == n));
    let r = [lsf_r1[0], lsf_r1[1], lsf_r2[0], lsf_r2[1]];
    let w = [wf1[0], wf1[1], wf2[0], wf2[1]].map(|v| v as i32);
    let aux: [i32; 4] = std::array::from_fn(|k| r[k] as i32 * w[k]);
    let mut best = [MAX_32; VQ_W];
    let mut bidx = [0i32; VQ_W];
    let mut all = [0i32; VQ_W];
    for b in (0..n).step_by(VQ_W) {
        let (mut d1, mut d2) = ([0i32; VQ_W], [0i32; VQ_W]);
        for k in 0..4 {
            let c: &[i16; VQ_W] = cols[k][b..b + VQ_W].try_into().unwrap();
            for l in 0..VQ_W {
                let (t1, t2) = if SIGNED {
                    (
                        ((w[k] * r[k].wrapping_sub(c[l]) as i32) >> 15) as i16 as i32,
                        ((w[k] * r[k].wrapping_add(c[l]) as i32) >> 15) as i16 as i32,
                    )
                } else {
                    (
                        ((aux[k].wrapping_sub(w[k] * c[l] as i32)) >> 15) as i16 as i32,
                        0,
                    )
                };
                d1[l] = d1[l].wrapping_add(t1 * t1);
                all[l] |= t1 * t1;
                if SIGNED {
                    d2[l] = d2[l].wrapping_add(t2 * t2);
                    all[l] |= t2 * t2;
                }
            }
        }
        for l in 0..VQ_W {
            let step = if SIGNED { 2 } else { 1 };
            let lt = d1[l] < best[l];
            best[l] = if lt { d1[l] } else { best[l] };
            bidx[l] = if lt { (step * (b + l)) as i32 } else { bidx[l] };
            if SIGNED {
                let lt = d2[l] < best[l];
                best[l] = if lt { d2[l] } else { best[l] };
                bidx[l] = if lt {
                    (2 * (b + l) + 1) as i32
                } else {
                    bidx[l]
                };
            }
        }
    }
    if all.iter().fold(0, |a, &v| a | v) >= 1 << 29 {
        return None;
    }
    let index = first_min(&best, &bidx);
    Some(index as usize)
}

/// A `Q_plsf_5` split: the C's codebook, its columns, signed (`Vq_subvec_s`).
pub(super) type Split5 = (&'static [i16], [&'static [i16]; 4], bool);

/// `Q_plsf_5`'s five splits, in order.
pub(super) fn splits5() -> [Split5; 5] {
    [
        (&DICO1_LSF_5, cols(&DICO1_5_SOA), false),
        (&DICO2_LSF_5, cols(&DICO2_5_SOA), false),
        (&DICO3_LSF_5, cols(&DICO3_5_SOA), true),
        (&DICO4_LSF_5, cols(&DICO4_5_SOA), false),
        (&DICO5_LSF_5, cols(&DICO5_5_SOA), false),
    ]
}

/// One `Q_plsf_5` split: the lanes, or the C's loop if they can't decide (see
/// `vq_pair_portable`). Overwrites the residual pairs with the (signed) codeword.
pub(super) fn vq5(
    r1: &mut [i16],
    r2: &mut [i16],
    dico: &[i16],
    cols: [&[i16]; 4],
    wf1: &[i16],
    wf2: &[i16],
    signed: bool,
) -> i16 {
    let n = cols[0].len();
    let lanes = if signed {
        vq_pair_s(r1, r2, &cols, wf1, wf2)
    } else {
        vq_pair(r1, r2, &cols, wf1, wf2)
    };
    match lanes {
        None if signed => vq_subvec_s(r1, r2, dico, wf1, wf2, n),
        None => vq_subvec(r1, r2, dico, wf1, wf2, n),
        Some(i) => {
            let (row, neg) = if signed {
                (i >> 1, i & 1 == 1)
            } else {
                (i, false)
            };
            let c = cols.map(|c| if neg { c[row].wrapping_neg() } else { c[row] });
            (r1[0], r1[1], r2[0], r2[1]) = (c[0], c[1], c[2], c[3]);
            i as i16
        }
    }
}

/// `Q_plsf_5`: quantises the two LSP vectors of an MR122 frame (5 indices).
pub(crate) fn q_plsf_5(
    st: &mut QPlsfState,
    lsp1: &[i16],
    lsp2: &[i16],
    lsp1_q: &mut [i16],
    lsp2_q: &mut [i16],
    indice: &mut [i16],
    ov: &mut Flag,
) {
    let mut lsf1 = [0i16; M];
    let mut lsf2 = [0i16; M];
    let mut wf1 = [0i16; M];
    let mut wf2 = [0i16; M];
    let mut lsf_p = [0i16; M];
    let mut lsf_r1 = [0i16; M];
    let mut lsf_r2 = [0i16; M];
    let mut lsf1_q = [0i16; M];
    let mut lsf2_q = [0i16; M];
    lsp_lsf(lsp1, &mut lsf1, M, ov);
    lsp_lsf(lsp2, &mut lsf2, M, ov);
    lsf_wt(&lsf1, &mut wf1);
    lsf_wt(&lsf2, &mut wf2);
    for i in 0..M {
        lsf_p[i] =
            MEAN_LSF_5[i].wrapping_add(((st.past_rq[i] as i32 * LSP_PRED_FAC_MR122) >> 15) as i16);
        lsf_r1[i] = lsf1[i].wrapping_sub(lsf_p[i]);
        lsf_r2[i] = lsf2[i].wrapping_sub(lsf_p[i]);
    }
    let (r1, r2) = (&mut lsf_r1, &mut lsf_r2);
    for (k, (dico, c, signed)) in splits5().into_iter().enumerate() {
        let j = 2 * k;
        indice[k] = vq5(
            &mut r1[j..],
            &mut r2[j..],
            dico,
            c,
            &wf1[j..],
            &wf2[j..],
            signed,
        );
    }
    for i in 0..M {
        lsf1_q[i] = lsf_r1[i].wrapping_add(lsf_p[i]);
        lsf2_q[i] = lsf_r2[i].wrapping_add(lsf_p[i]);
        st.past_rq[i] = lsf_r2[i];
    }
    reorder_lsf(&mut lsf1_q, LSF_GAP, M, ov);
    reorder_lsf(&mut lsf2_q, LSF_GAP, M, ov);
    lsf_lsp(&lsf1_q, lsp1_q, M, ov);
    lsf_lsp(&lsf2_q, lsp2_q, M, ov);
}

/// `lspState`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct LspState {
    pub lsp_old: [i16; M],
    pub lsp_old_q: [i16; M],
    pub q_st: QPlsfState,
}

impl Default for LspState {
    fn default() -> Self {
        Self {
            lsp_old: LSP_INIT_DATA,
            lsp_old_q: LSP_INIT_DATA,
            q_st: QPlsfState::default(),
        }
    }
}

/// `lsp`: A(z) → LSPs, interpolation of the unquantised filters into `az`, quantisation (unless
/// `used_mode` is `MRDTX`) into `az_q`; writes the LSF indices to `ana[..3]` (or `[..5]`) and
/// returns how many it wrote.
pub(crate) fn lsp(
    st: &mut LspState,
    req_mode: CodecMode,
    used_mode: CodecMode,
    az: &mut [i16; 4 * MP1],
    az_q: &mut [i16; 4 * MP1],
    lsp_new: &mut [i16; M],
    ana: &mut [i16],
    ov: &mut Flag,
) -> usize {
    let mut lsp_new_q = [0i16; M];
    let mut written = 0;
    if req_mode == CodecMode::MR122 {
        let mut lsp_mid = [0i16; M];
        let mut lsp_mid_q = [0i16; M];
        az_lsp(&az[MP1..2 * MP1], &mut lsp_mid, &st.lsp_old);
        az_lsp(&az[3 * MP1..], lsp_new, &lsp_mid);
        int_lpc_1and3_2(&st.lsp_old, &lsp_mid, lsp_new, az, ov);
        if used_mode != CodecMode::MRDTX {
            q_plsf_5(
                &mut st.q_st,
                &lsp_mid,
                lsp_new,
                &mut lsp_mid_q,
                &mut lsp_new_q,
                ana,
                ov,
            );
            int_lpc_1and3(&st.lsp_old_q, &lsp_mid_q, &lsp_new_q, az_q, ov);
            written = 5;
        }
    } else {
        az_lsp(&az[3 * MP1..], lsp_new, &st.lsp_old);
        int_lpc_1to3_2(&st.lsp_old, lsp_new, az, ov);
        if used_mode != CodecMode::MRDTX {
            let mut pred_init_i = 0;
            q_plsf_3(
                &mut st.q_st,
                req_mode,
                lsp_new,
                &mut lsp_new_q,
                ana,
                &mut pred_init_i,
                ov,
            );
            int_lpc_1to3(&st.lsp_old_q, &lsp_new_q, az_q, ov);
            written = 3;
        }
    }
    st.lsp_old = *lsp_new;
    if used_mode != CodecMode::MRDTX {
        st.lsp_old_q = lsp_new_q;
    }
    written
}

#[cfg(all(test, oracle))]
pub(crate) mod vs_ref;
