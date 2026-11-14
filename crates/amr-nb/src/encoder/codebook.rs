//! Algebraic (fixed) codebook search: `cbsearch` and the per-mode searches `code_2i40_9bits`,
//! `code_2i40_11bits`, `code_3i40_14bits`, `code_4i40_17bits`, `code_8i40_31bits`,
//! `code_10i40_35bits`, with their helpers `cor_h_x`, `cor_h_x2`, `cor_h`, `set_sign`,
//! `set_sign12k2` and `search_10and8i40` (as `search_10i40`/`search_8i40_portable`).
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `enc/src/{cbsearch,c2_9pf,c2_11pf,c3_14pf,
//! c4_17pf,c8_31pf,c1035pf,s10_8pf,set_sign,cor_h,cor_h_x,cor_h_x2}.cpp`.
//!
//! The C's `h` points 40 samples into a zeroed buffer (`hvec`), so `h[-k]` reads 0; the ports
//! read zero-padded copies of `h` where they index it that way.
#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1

use super::lpc::l_abs;
use crate::basic_op::{
    Flag, MAX_16, MAX_32, MIN_16, add, l_add, l_mac, l_mac_energy_from, l_msu, l_mult, l_shl,
    max_abs, mult, negate, norm_l, round, shl, sum_abs,
};
use crate::cnst::{L_SUBFR, Mode};
use crate::math::inv_sqrt;
use crate::tables::{GRAY, STARTPOS, TRACKTABLE};

pub(crate) const L_CODE: usize = 40;
const NB_TRACK: usize = 5;
const STEP: usize = 5;
const NB_TRACK_MR102: usize = 4;
const STEP_MR102: usize = 4;

pub(crate) use crate::simd::{
    cor_h_blocks, cor_h_search_8i40, cor_h_tm, filtered_code, search_2i40, search_2i40_sat,
    search_3i40, search_4i40, search_10i40,
};

/// `cor_h_x`: backward-filtered target `dn = <x, h shifted>`, normalised (`sf` headroom).
pub(super) fn cor_h_x(h: &[i16], x: &[i16], dn: &mut [i16; L_CODE], sf: i16, ov: &mut Flag) {
    cor_h_x2(h, x, dn, sf, NB_TRACK, STEP, ov);
}

/// `cor_h_x2`'s y32[i] = the C's saturating chain `L_mac(s, x[j], h[j - i])` over j >= i. With
/// max|x| * sum|h| < 2^30 no partial sum saturates: then twice the (vectorised) plain sums.
fn corr_hx_exact(h: &[i16], x: &[i16], ov: &mut Flag) -> [i32; L_CODE] {
    if max_abs(&x[..L_CODE]) * sum_abs(&h[..L_CODE]) < 1 << 30 {
        crate::simd::corr_hx(h, x).map(|v| v << 1)
    } else {
        core::array::from_fn(|i| (i..L_CODE).fold(0, |s, j| l_mac(s, x[j], h[j - i], ov)))
    }
}

/// `cor_h_x`'s correlations y[i] = sum over j >= i of x[j] * h[j - i] (wrapping), in blocks of
/// 8 outputs (lanes) with x zero-padded, so every tap is one vector step (as in `convolve`).
/// The C's per-output triangle stayed scalar (359 `imul`s).
#[inline(always)]
pub(crate) fn corr_hx_portable(h: &[i16], x: &[i16]) -> [i32; L_CODE] {
    let h: &[i16; L_CODE] = h[..L_CODE].try_into().unwrap();
    let mut xp = [0i16; L_CODE + 8];
    xp[..L_CODE].copy_from_slice(&x[..L_CODE]);
    let mut y = [0i32; L_CODE];
    for (b, yb) in y.as_chunks_mut::<8>().0.iter_mut().enumerate() {
        let i0 = 8 * b;
        let mut acc = [0i32; 8];
        for (m, &hm) in h[..L_CODE - i0].iter().enumerate() {
            let w: &[i16; 8] = xp[i0 + m..i0 + m + 8].try_into().unwrap();
            for (a, &xv) in acc.iter_mut().zip(w) {
                *a = a.wrapping_add(hm as i32 * xv as i32);
            }
        }
        *yb = acc;
    }
    y
}

/// `cor_h_x2`: as `cor_h_x` with a configurable track layout (MR102). TS 26.073's saturating
/// `L_abs`/`L_add`/`round` (opencore wraps: wrong on full-scale input).
pub(super) fn cor_h_x2(
    h: &[i16],
    x: &[i16],
    dn: &mut [i16; L_CODE],
    sf: i16,
    nb_track: usize,
    step: usize,
    ov: &mut Flag,
) {
    let y32 = corr_hx_exact(h, x, ov);
    let mut tot: i32 = 5;
    for k in 0..nb_track {
        let max = (k..L_CODE)
            .step_by(step)
            .map(|i| l_abs(y32[i]))
            .fold(0, i32::max);
        tot = l_add(tot, max >> 1, ov);
    }
    let j = norm_l(tot) - sf;
    for i in 0..L_CODE {
        dn[i] = round(l_shl(y32[i], j, ov), ov);
    }
}

/// `cor_h`'s normalised impulse response `h2` (shared by `cor_h` and `cor_h_tm`), as TS 26.073
/// `cor_h.c` with its saturating operators. (opencore's rewrite takes the `h >> 1` branch only
/// on wrap, not at `extract_h(s) == 32767`, uses 32440 for `mult(32767, 32440)` and truncates
/// `h * k` instead of saturating: wrong on full-scale input.)
#[inline(always)]
fn cor_h_h2(h: &[i16], ov: &mut Flag) -> [i16; L_CODE] {
    let h = &h[..L_CODE];
    // 2 + sum 2 h^2, saturating: the terms are >= 0 (-32768^2 saturates alone).
    let sum: i64 = 2 + h.iter().map(|&v| 2 * v as i64 * v as i64).sum::<i64>();
    let s = sum.min(MAX_32 as i64) as i32;
    if (s >> 16) == 32767 {
        return core::array::from_fn(|i| h[i] >> 1);
    }
    let k = (l_shl(inv_sqrt(s >> 1, ov), 7, ov) >> 16) as i16;
    let k = mult(k, 32440, ov);
    core::array::from_fn(|i| round(l_shl(l_mult(h[i], k, ov), 9, ov), ov))
}

/// Whether `cor_h`'s sums can saturate for this `h2`: every partial sum of the matrix (and its
/// `round`) is bounded by e = sum h2^2 (Cauchy-Schwarz), so below 2e + 0x8000 <= MAX_32 the
/// wrapping kernels are exact.
#[inline(always)]
fn rr_may_saturate(h2: &[i16; L_CODE]) -> bool {
    let e: i64 = h2.iter().map(|&v| v as i64 * v as i64).sum();
    2 * e + 0x8000 > MAX_32 as i64
}

/// `cor_h`'s matrix exactly as TS 26.073 builds it (saturating `L_mac` chains). The rare
/// fallback of the fast kernels (`rr_may_saturate`).
#[cold]
fn cor_h_exact(h2: &[i16; L_CODE], sign: &[i16; L_CODE], ov: &mut Flag) -> [[i16; L_CODE]; L_CODE] {
    let mut m = [[0i16; L_CODE]; L_CODE];
    let mut s: i32 = 0;
    for k in 0..L_CODE {
        s = l_mac(s, h2[k], h2[k], ov);
        let i = L_CODE - 1 - k;
        m[i][i] = round(s, ov);
    }
    for dec in 1..L_CODE {
        let mut s: i32 = 0;
        for k in 0..L_CODE - dec {
            let (j, i) = (L_CODE - 1 - k, L_CODE - 1 - dec - k);
            s = l_mac(s, h2[k], h2[k + dec], ov);
            let v = mult(round(s, ov), mult(sign[i], sign[j], ov), ov);
            m[j][i] = v;
            m[i][j] = v;
        }
    }
    m
}

/// [`cor_h_exact`] into a track-major `Tm`: the blocks `keep` selects (by track pair) and every
/// diagonal.
#[cold]
fn cor_h_exact_tm<const T: usize, const N: usize, const W: usize>(
    h2: &[i16; L_CODE],
    sign: &[i16; L_CODE],
    rr: &mut Tm<T, N, W>,
    keep: impl Fn(usize, usize) -> bool,
    ov: &mut Flag,
) {
    let m = cor_h_exact(h2, sign, ov);
    for i in 0..L_CODE {
        for j in 0..L_CODE {
            if keep(i % T, j % T) || i == j {
                rr.b[i % T][j % T][i / T][j / T] = m[i][j];
            }
        }
        rr.diag[i % T][i / T] = m[i][i];
    }
}

/// `rr` in track-major order: `T` tracks of `N` positions (`i = t + T * n`), rows padded to
/// `W >= N` lanes with zeros: `b[ta][tb][ia][ib] = rr[ta + T * ia][tb + T * ib]`,
/// `diag[t][n] = rr[t + T * n][t + T * n]`. The searches read every candidate set as one
/// contiguous row.
#[derive(Clone)]
pub(crate) struct Tm<const T: usize, const N: usize, const W: usize> {
    pub(crate) b: [[[[i16; W]; N]; T]; T],
    pub(crate) diag: [[i16; W]; T],
}

/// The 5-track modes, MR475–MR795 and MR122 (5 tracks of 8).
pub(crate) type RrTm = Tm<NB_TRACK, 8, 8>;
/// MR102 (4 tracks of 10) in rows of `W`: 16 with AVX2 (whole 16-lane rows), 10 portable
/// (SSE2 holds 16 lanes' sums in no register file: `cor_h_tm` 2090 vs 1000 ns).
pub(crate) type RrTm4<const W: usize> = Tm<NB_TRACK_MR102, 10, W>;

impl<const T: usize, const N: usize, const W: usize> Tm<T, N, W> {
    pub(crate) const ZERO: Self = Self {
        b: [[[[0; W]; N]; T]; T],
        diag: [[0; W]; T],
    };
    /// `rr[i][j]`.
    #[inline(always)]
    pub(crate) fn at(&self, i: usize, j: usize) -> i16 {
        self.b[i % T][j % T][i / T][j / T]
    }
    /// `rr[q][t + T * n]` for n = 0..N: position `q` against every position of track `t`.
    #[inline(always)]
    fn row(&self, q: usize, t: usize) -> &[i16; W] {
        &self.b[q % T][t][q / T]
    }
}

/// `cor_h` writing a track-major `Tm` (T tracks of N, T * N = 40). Every cell, either
/// triangle, is the wrapping sum rr[m][j] = sum over k <= 39 - m of h2[k] * h2[k - δ],
/// δ = j - m (h2 is 0 outside 0..40; for j > m substitute k -> k + δ in the C's sum). Wrapping
/// sums don't depend on the order, so step k once, advance the sums, and write row m = 39 - k
/// whole (no mirror pass).
///
/// At step k only δ in k-39..=k still matter (later rows need δ >= k - 39, products vanish
/// for δ > k): N per residue r of δ mod T, and residues never mix. `acc[r]` holds those N
/// (lane p: δ = r + T(q + p) - 40, q = ceil((k + 1 - r) / T)), shifted down a lane when q
/// advances (k ≡ r mod T; padding lanes stay 0). Then the product's source is the same N lanes
/// of h2 every step, and the lanes are row m's cells of track tb = (r - k - 1) mod T
/// (j = tb + T ib). Steps run unrolled by T (`tm_step`), so every index is a constant and the
/// T sums stay in registers (as loops over r or with moving offsets, index math and
/// store-forwarding stalls cost 1.5–3×, RESULTS.md).
/// Portable body of the dispatched `cor_h_tm` and of MR102's `cor_h_search_8i40` (`NARROW`:
/// see the sign scaling).
#[inline(always)]
pub(crate) fn cor_h_tm_portable<
    const T: usize,
    const N: usize,
    const W: usize,
    const NARROW: bool,
>(
    h: &[i16],
    sign: &[i16; L_CODE],
    rr: &mut Tm<T, N, W>,
    ov: &mut Flag,
) {
    const { assert!(T * N == L_CODE && N <= W && T <= 5) };
    let h2 = cor_h_h2(h, ov);
    if rr_may_saturate(&h2) {
        return cor_h_exact_tm(&h2, sign, rr, |_, _| true, ov);
    }
    // hd[t][p] = h2[39 - (t + Tp)]; sg[t][p] = sign[t + Tp].
    let mut hd = [[0i16; W]; T];
    let mut sg = [[0i16; W]; T];
    for x in 0..L_CODE {
        hd[x % T][x / T] = h2[L_CODE - 1 - x];
        sg[x % T][x / T] = sign[x];
    }
    let mut acc = [[0i32; W]; T];
    macro_rules! steps {
        ($pm:literal, $kb:ident: $($k:literal)*) => {$(
            tm_step::<T, N, W, $k, 0, $pm>($kb, &h2, &hd, &sg, sign, &mut acc, rr);
            tm_step::<T, N, W, $k, 1, $pm>($kb, &h2, &hd, &sg, sign, &mut acc, rr);
            tm_step::<T, N, W, $k, 2, $pm>($kb, &h2, &hd, &sg, sign, &mut acc, rr);
            tm_step::<T, N, W, $k, 3, $pm>($kb, &h2, &hd, &sg, sign, &mut acc, rr);
            tm_step::<T, N, W, $k, 4, $pm>($kb, &h2, &hd, &sg, sign, &mut acc, rr);
        )*};
    }
    // set_sign's and set_sign12k2's signs are all ±32767: the sign product then needs no
    // multiply, which pays with narrow lanes (SSE2); with AVX2 the multiplies are cheaper.
    if NARROW && sign.iter().all(|&s| s == 32767 || s == -32767) {
        // sg as masks: -1 for a negative sign.
        for t in sg.iter_mut().flatten() {
            *t >>= 15;
        }
        for kb in 0..N {
            steps!(true, kb: 0 1 2 3 4);
        }
    } else {
        for kb in 0..N {
            steps!(false, kb: 0 1 2 3 4);
        }
    }
    let mut d: i32 = 0;
    for k in 0..L_CODE {
        d = d.wrapping_add(h2[k] as i32 * h2[k] as i32);
        let (m, d0) = (L_CODE - 1 - k, (d.wrapping_add(0x4000) >> 15) as i16);
        rr.b[m % T][m % T][m / T][m / T] = d0;
        rr.diag[m % T][m / T] = d0;
    }
}

/// `cor_h_tm` for the given blocks `rr.b[ta][tb]` only (and every diagonal): the 2-pulse
/// searches read one block per track pair (rows on the second pulse's track, lanes on the
/// first's). In `cor_h_tm`'s scheme a block takes all its cells from one residue,
/// r = (tb − ta) mod 5 (`tm_step`: tb − ta ≡ R), and writes its row m = ta + 5n at step
/// k = 39 − m: only that residue's sums advance, and only the block's 8 rows are scaled.
/// The other blocks are left as they are.
#[inline(always)]
pub(crate) fn cor_h_blocks_portable(
    h: &[i16],
    sign: &[i16; L_CODE],
    rr: &mut RrTm,
    blocks: &[[usize; 2]; 2],
    ov: &mut Flag,
) {
    const T: usize = NB_TRACK;
    let h2 = cor_h_h2(h, ov);
    if rr_may_saturate(&h2) {
        return cor_h_exact_tm(&h2, sign, rr, |ta, tb| blocks.contains(&[ta, tb]), ov);
    }
    // hd[t][p] = h2[39 - (t + Tp)] (as `cor_h_tm`)
    let mut hd = [[0i16; 8]; T];
    for x in 0..L_CODE {
        hd[x % T][x / T] = h2[L_CODE - 1 - x];
    }
    // blocks with equal residues share one accumulation (half the MR475/MR515 subframes)
    let [r0, r1] = blocks.map(|[ta, tb]| (tb + T - ta) % T);
    let (ta0, ta1) = (blocks[0][0], blocks[1][0]);
    if r0 == r1 {
        cor_h_residue(r0, &h2, &hd, sign, &[ta0, ta1], rr);
    } else {
        cor_h_residue(r0, &h2, &hd, sign, &[ta0], rr);
        cor_h_residue(r1, &h2, &hd, sign, &[ta1], rr);
    }
    let mut d: i32 = 0;
    for k in 0..L_CODE {
        d = d.wrapping_add(h2[k] as i32 * h2[k] as i32);
        let (m, d0) = (L_CODE - 1 - k, (d.wrapping_add(0x4000) >> 15) as i16);
        rr.b[m % T][m % T][m / T][m / T] = d0;
        rr.diag[m % T][m / T] = d0;
    }
}

/// Residue `r` (as a constant: the step's row of `hd` and the lane shift are fixed) of
/// `cor_h_blocks_portable`, for the blocks (ta, ta + r) of the given `tas`.
#[inline(always)]
fn cor_h_residue<const N: usize>(
    r: usize,
    h2: &[i16; L_CODE],
    hd: &[[i16; 8]; NB_TRACK],
    sign: &[i16; L_CODE],
    tas: &[usize; N],
    rr: &mut RrTm,
) {
    match r {
        0 => cor_h_block::<0, N>(h2, hd, sign, tas, rr),
        1 => cor_h_block::<1, N>(h2, hd, sign, tas, rr),
        2 => cor_h_block::<2, N>(h2, hd, sign, tas, rr),
        3 => cor_h_block::<3, N>(h2, hd, sign, tas, rr),
        _ => cor_h_block::<4, N>(h2, hd, sign, tas, rr),
    }
}

/// Residue R's sums over the 40 steps (unrolled by 5, as `tm_step`); for each block
/// (ta, ta + R), row m = ta + 5n written at step 39 − m.
#[inline(always)]
fn cor_h_block<const R: usize, const N: usize>(
    h2: &[i16; L_CODE],
    hd: &[[i16; 8]; NB_TRACK],
    sign: &[i16; L_CODE],
    tas: &[usize; N],
    rr: &mut RrTm,
) {
    const T: usize = NB_TRACK;
    let sgs: [[i32; 8]; N] =
        tas.map(|ta| core::array::from_fn(|l| sign[(ta + R) % T + T * l] as i32));
    let mut a = [0i32; 8];
    for kb in 0..8 {
        block_step::<R, 0, N>(kb, h2, hd, sign, &sgs, tas, &mut a, rr);
        block_step::<R, 1, N>(kb, h2, hd, sign, &sgs, tas, &mut a, rr);
        block_step::<R, 2, N>(kb, h2, hd, sign, &sgs, tas, &mut a, rr);
        block_step::<R, 3, N>(kb, h2, hd, sign, &sgs, tas, &mut a, rr);
        block_step::<R, 4, N>(kb, h2, hd, sign, &sgs, tas, &mut a, rr);
    }
}

/// Step k = 5 kb + KAPPA of `cor_h_block` (as `tm_step` for one residue).
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn block_step<const R: usize, const KAPPA: usize, const N: usize>(
    kb: usize,
    h2: &[i16; L_CODE],
    hd: &[[i16; 8]; NB_TRACK],
    sign: &[i16; L_CODE],
    sgs: &[[i32; 8]; N],
    tas: &[usize; N],
    a: &mut [i32; 8],
    rr: &mut RrTm,
) {
    const T: usize = NB_TRACK;
    if KAPPA == R {
        *a = core::array::from_fn(|l| if l + 1 < 8 { a[l + 1] } else { 0 });
    }
    let k = T * kb + KAPPA;
    let (hk, d) = (h2[k] as i32, &hd[(R + 2 * T - KAPPA - 1) % T]);
    for l in 0..8 {
        a[l] = a[l].wrapping_add(hk * d[l] as i32);
    }
    for (&ta, sg) in tas.iter().zip(sgs) {
        if ta == T - 1 - KAPPA {
            let sm = sign[L_CODE - 1 - k] as i32;
            let row = &mut rr.b[ta][(ta + R) % T][7 - kb];
            for l in 0..8 {
                let tmp1 = (a[l].wrapping_add(0x4000) >> 15) as i16;
                let tmp2 = ((sm * sg[l]) >> 15).min(32767) as i16; // mult: -32768² saturates
                row[l] = ((tmp1 as i32 * tmp2 as i32) >> 15) as i16;
            }
        }
    }
}

/// One residue R of `cor_h_tm`'s step k = T kb + KAPPA (constant indices, so the T residues'
/// sums stay in registers and interleave; R or KAPPA >= T: no step). `PM`: `sg` holds sign
/// masks (see `pm_scale`).
#[inline(always)]
fn tm_step<
    const T: usize,
    const N: usize,
    const W: usize,
    const KAPPA: usize,
    const R: usize,
    const PM: bool,
>(
    kb: usize,
    h2: &[i16; L_CODE],
    hd: &[[i16; W]; T],
    sg: &[[i16; W]; T],
    sign: &[i16; L_CODE],
    acc: &mut [[i32; W]; T],
    rr: &mut Tm<T, N, W>,
) {
    if KAPPA >= T || R >= T {
        return;
    }
    let k = T * kb + KAPPA;
    let a = &mut acc[R];
    if KAPPA == R {
        *a = core::array::from_fn(|l| if l + 1 < W { a[l + 1] } else { 0 });
    }
    let (ta, tb) = (T - 1 - KAPPA, (R + 2 * T - KAPPA - 1) % T);
    let (hk, sm) = (h2[k] as i32, sign[L_CODE - 1 - k]);
    for l in 0..W {
        a[l] = a[l].wrapping_add(hk * hd[tb][l] as i32);
    }
    let row = &mut rr.b[ta][tb][N - 1 - kb];
    for l in 0..W {
        let tmp1 = (a[l].wrapping_add(0x4000) >> 15) as i16;
        row[l] = if PM {
            pm_scale(tmp1, (sm >> 15) ^ sg[tb][l])
        } else {
            let tmp2 = ((sm as i32 * sg[tb][l] as i32) >> 15).min(32767) as i16; // mult: -32768² saturates
            ((tmp1 as i32 * tmp2 as i32) >> 15) as i16
        };
    }
}

/// `(t * tmp2) >> 15` for the sign product of two ±32767 signs, `tmp2 = (si * sj) >> 15`:
/// 32766 (same signs: t - ceil(t / 16384)) or -32767 (opposite: -t + floor(t / 32768)).
/// `opp` is -1 for opposite signs, else 0. Branch-free in i16 lanes, so it vectorises 8 wide:
/// ceil(t / 16384) = (t >> 14) + 1 + (-1 if t's low 14 bits are 0), all mod 2^16.
#[inline(always)]
pub(crate) fn pm_scale(t: i16, opp: i16) -> i16 {
    let zero = -(((t & 0x3fff) == 0) as i16);
    let same = t.wrapping_sub(t >> 14).wrapping_sub(1).wrapping_sub(zero);
    let opposite = (t >> 15).wrapping_sub(t);
    (same & !opp) | (opposite & opp)
}

/// `set_sign`: pulse signs from `dn` (made absolute); `dn2` marks the `8 - n` weakest positions
/// per track with -1.
pub(super) fn set_sign(
    dn: &mut [i16; L_CODE],
    sign: &mut [i16; L_CODE],
    dn2: &mut [i16; L_CODE],
    n: usize,
) {
    for i in (0..L_CODE).rev() {
        let mut val = dn[i];
        if val >= 0 {
            sign[i] = 32767;
        } else {
            sign[i] = -32767;
            val = negate(val);
            dn[i] = val;
        }
        dn2[i] = val;
    }
    let mut pos = 0;
    for i in 0..NB_TRACK {
        for _ in 0..8 - n {
            let mut min = 0x7fff;
            for j in (i..L_CODE).step_by(STEP) {
                if dn2[j] >= 0 && dn2[j] < min {
                    min = dn2[j];
                    pos = j;
                }
            }
            dn2[pos] = -1;
        }
    }
}

/// `set_sign12k2`'s `round(L_shl(L_mac(L_mult(k_cn, cn), k_dn, dn), 10))` per sample, for
/// k_cn, k_dn >= 0 (no L_mult saturates): t = k_cn cn + k_dn dn is exact in i32 and the chain
/// is round(sat(t << 11)), clamped where (t << 11) + 0x8000 leaves i32 (selects, so the 40
/// lanes vectorise; as `total_exc`). Unflagged: the encoder reads `ov` only in `g_pitch`,
/// which clears it first.
pub(super) fn sign_cor(
    k_cn: i16,
    cn: &[i16; L_CODE],
    k_dn: i16,
    dn: &[i16; L_CODE],
) -> [i16; L_CODE] {
    let (hi, lo) = ((MAX_32 - 0x8000) >> 11, i32::MIN >> 11);
    core::array::from_fn(|i| {
        let t = k_cn as i32 * cn[i] as i32 + k_dn as i32 * dn[i] as i32;
        let r = (t.wrapping_shl(11).wrapping_add(0x8000) >> 16) as i16;
        if t > hi {
            i16::MAX
        } else if t < lo {
            i16::MIN
        } else {
            r
        }
    })
}

/// `set_sign12k2`: pulse signs from a mix of `dn` and the LTP residual `cn`; best position per
/// track (`pos_max`) and the track order to search (`ipos`, 2 × `nb_track`).
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn set_sign12k2(
    dn: &mut [i16; L_CODE],
    cn: &[i16],
    sign: &mut [i16; L_CODE],
    pos_max: &mut [i16],
    nb_track: usize,
    ipos: &mut [i16],
    step: usize,
    ov: &mut Flag,
) {
    let mut en = [0i16; L_CODE];
    let cn: &[i16; L_CODE] = cn[..L_CODE].try_into().unwrap();
    // TS 26.073's operators throughout (opencore sums dn's energy wrapping, skips L_shl's
    // saturation in k_dn and keeps its L_mac quirk). 40 samples, MR122 only.
    let s = l_mac_energy_from(256, cn, ov);
    let k_cn = (l_shl(inv_sqrt(s, ov), 5, ov) >> 16) as i16;
    let t = l_mac_energy_from(256, dn, ov);
    let k_dn = (l_shl(inv_sqrt(t, ov), 5, ov) >> 16) as i16;
    let fast = k_cn >= 0 && k_dn >= 0;
    let cors = sign_cor(k_cn, cn, k_dn, dn);
    for i in 0..L_CODE {
        let val = dn[i];
        let cor = if fast {
            cors[i]
        } else {
            let s = l_mac(l_mult(k_cn, cn[i], ov), k_dn, val, ov);
            round(l_shl(s, 10, ov), ov)
        };
        let neg = cor < 0;
        sign[i] = if neg { -32767 } else { 32767 };
        en[i] = if neg { negate(cor) } else { cor };
        dn[i] = if neg { negate(val) } else { val };
    }
    let mut max_of_all = -1;
    let mut pos = 0usize;
    for i in 0..nb_track {
        let mut max = -1;
        for j in (i..L_CODE).step_by(step) {
            if en[j] > max {
                max = en[j];
                pos = j;
            }
        }
        pos_max[i] = pos as i16;
        if max > max_of_all {
            max_of_all = max;
            ipos[0] = i as i16;
        }
    }
    let mut pos = ipos[0] as usize;
    ipos[nb_track] = pos as i16;
    for i in 1..nb_track {
        pos += 1;
        if pos >= nb_track {
            pos = 0;
        }
        ipos[i] = pos as i16;
        ipos[i + nb_track] = pos as i16;
    }
}

/// Pitch sharpening of `v[T0..]` (the codebook searches' `if (T0 < L_CODE)` loops).
pub(super) fn sharpen(v: &mut [i16], t0: i16, sharp: i16, ov: &mut Flag) {
    for i in t0.max(0) as usize..L_CODE {
        let temp = mult(v[i - t0 as usize], sharp, ov);
        v[i] = add(v[i], temp, ov);
    }
}

/// One pulse-position loop of the depth-first searches, for 8 first pulses at once: lane `l`
/// runs the C's loop over its track `t`'s 8 positions (`d` row n: `dn[t + 5n]`, `diag` row
/// n: `rr[t + 5n][t + 5n]`, see `Rows`) for its own correlation `ps0[l]` and energy
/// `alp0[l]`; `x[n][l]` sums the lane's cross terms `rr[q][t + 5n]` over its already placed
/// pulses `q`; `sd`/`sc` are the shifts of the diagonal/cross energy terms. TS 26.073's
/// operators: `add` (saturating `ps`), `mult` (`sq`), `round`, and the compare
/// `L_msu(L_mult(alp, sq1), sq, alp_16) > 0`, i.e. `alp * sq1 > sq * alp_16` (a saturated
/// difference keeps its sign): all exact here. The energy is a saturating `L_mac` chain, but
/// with the callers' `alp0` and shifts (`STEP5_SHIFTS`, proven below) no partial sum leaves
/// the i32 range: the wrapping total is the C's. Returns each lane's `(sq, ps, alp, ix)`.
#[inline(always)]
pub(crate) fn lanes_next(
    ps0: &[i16; 8],
    alp0: &[i32; 8],
    d: &impl Rows,
    diag: &impl Rows,
    x: &[[i32; 8]; 8],
    sd: u32,
    sc: u32,
) -> Lanes {
    debug_assert!(
        STEP5_SHIFTS[1..]
            .iter()
            .any(|&(a, b, _)| (a, b) == (sd, sc)),
        "lanes_next: shifts outside STEP5_SHIFTS' proof"
    );
    let (mut sq, mut ps, mut alp, mut ix) = ([-1i16; 8], [0i16; 8], [1i16; 8], [0i16; 8]);
    for n in 0..8 {
        let (d, diag) = (d.row(n), diag.row(n));
        // Whole 8-lane rows per operation, selects as masks: the form SLP vectorises.
        let ps1: [i16; 8] = core::array::from_fn(|l| ps0[l].saturating_add(d[l]));
        let sq1: [i16; 8] =
            core::array::from_fn(|l| ((ps1[l] as i32 * ps1[l] as i32) >> 15).min(32767) as i16);
        let a16: [i16; 8] = core::array::from_fn(|l| {
            let alp1 = alp0[l]
                .wrapping_add((diag[l] as i32) << sd)
                .wrapping_add(x[n][l] << sc);
            (alp1.wrapping_add(0x8000) >> 16) as i16 // `round`: never saturates here (below)
        });
        let m: [i16; 8] = core::array::from_fn(|l| {
            -((alp[l] as i32 * sq1[l] as i32 > sq[l] as i32 * a16[l] as i32) as i16)
        });
        let sel = |a: i16, b: i16, m: i16| (a & m) | (b & !m);
        sq = core::array::from_fn(|l| sel(sq1[l], sq[l], m[l]));
        ps = core::array::from_fn(|l| sel(ps1[l], ps[l], m[l]));
        alp = core::array::from_fn(|l| sel(a16[l], alp[l], m[l]));
        ix = core::array::from_fn(|l| sel(n as i16, ix[l], m[l]));
    }
    Lanes { sq, ps, alp, ix }
}

/// The step-5 searches' shifts (the C's weights) per pulse p = 1..=3: (sd, sc) for
/// `lanes_next`'s diagonal and cross terms, and `.2`, the shift that makes the next pulse's
/// alp0 from this pulse's `alp` (index 0: the first pulse's `rr` diagonal << 14). The
/// `const` assertion below proves the chains with these shifts; both searches read them here.
pub(crate) const STEP5_SHIFTS: [(u32, u32, u32); 4] =
    [(0, 0, 14), (14, 15, 14), (12, 13, 16), (12, 13, 0)];

/// The step-5 searches' energy chains never leave the i32 range, for any 16-bit `rr`: they
/// never saturate, wrapping sums are the C's, and `round` never saturates (<= 0x7fff_7fff).
/// Pulse p's chain: alp0 (from the previous `alp`, or `rr` for p = 1), the diagonal term,
/// and p cross terms (one per earlier pulse), each at its 16-bit extremes; `alp` is `round`
/// of the chain. Every partial sum is a sub-sum, so the extremes bound it too. (Before this
/// was proven, a per-row guard, its |term| sums and an exact scalar fallback ran on every
/// block: RESULTS loop #20.) Three pulses after the first (`step5_block` asserts P <= 4).
const _: () = {
    let (lo, hi) = (-32768i64, 32767i64);
    let mut alp = (lo, hi); // pulse 1 starts from the rr diagonal
    let mut p = 1;
    while p <= 3 {
        let (sd, sc, _) = STEP5_SHIFTS[p];
        let sh0 = STEP5_SHIFTS[p - 1].2;
        let c = p as i64;
        let lo_sum = (alp.0 << sh0) + (lo << sd) + c * (lo << sc);
        let hi_sum = (alp.1 << sh0) + (hi << sd) + c * (hi << sc);
        assert!(lo_sum >= i32::MIN as i64 && hi_sum <= 0x7fff_7fff);
        alp = ((lo_sum + 0x8000) >> 16, (hi_sum + 0x8000) >> 16);
        p += 1;
    }
};

/// `lanes_next`'s per-position inputs: one value for all lanes (`[i16; 8]`, by position `n`)
/// or one per half block (`[&[i16; 8]; 2]`, `step5_halves`).
pub(crate) trait Rows {
    fn row(&self, n: usize) -> [i16; 8];
}

impl Rows for [i16; 8] {
    #[inline(always)]
    fn row(&self, n: usize) -> [i16; 8] {
        [self[n]; 8]
    }
}

/// Two half blocks' tracks: lanes `0..4` from the first, `4..8` from the second.
impl Rows for [&[i16; 8]; 2] {
    #[inline(always)]
    fn row(&self, n: usize) -> [i16; 8] {
        let mut r = [self[1][n]; 8];
        r[..4].fill(self[0][n]);
        r
    }
}

/// Per-lane results of `lanes_next` (`ix`: position index on the track, 0..8).
#[derive(Clone, Copy, Default)]
pub(crate) struct Lanes {
    pub(crate) sq: [i16; 8],
    pub(crate) ps: [i16; 8],
    pub(crate) alp: [i16; 8],
    pub(crate) ix: [i16; 8],
}

/// 8 × 8 transpose as the three interleave rounds (16-, 32-, 64-bit units) of `punpck*`
/// (`simd::transpose8` runs them as SSE2 intrinsics on x86_64).
// Used on other targets; on x86_64 only `transpose8_is_the_transpose` reads it.
#[cfg_attr(all(target_arch = "x86_64", not(test)), allow(dead_code))]
#[inline(always)]
pub(crate) fn transpose8_portable(r: [[i16; 8]; 8]) -> [[i16; 8]; 8] {
    let lo = |a: [i16; 8], b: [i16; 8]| [a[0], b[0], a[1], b[1], a[2], b[2], a[3], b[3]];
    let hi = |a: [i16; 8], b: [i16; 8]| [a[4], b[4], a[5], b[5], a[6], b[6], a[7], b[7]];
    let lo2 = |a: [i16; 8], b: [i16; 8]| [a[0], a[1], b[0], b[1], a[2], a[3], b[2], b[3]];
    let hi2 = |a: [i16; 8], b: [i16; 8]| [a[4], a[5], b[4], b[5], a[6], a[7], b[6], b[7]];
    let lo4 = |a: [i16; 8], b: [i16; 8]| [a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3]];
    let hi4 = |a: [i16; 8], b: [i16; 8]| [a[4], a[5], a[6], a[7], b[4], b[5], b[6], b[7]];
    let s = [
        lo(r[0], r[1]),
        hi(r[0], r[1]),
        lo(r[2], r[3]),
        hi(r[2], r[3]),
        lo(r[4], r[5]),
        hi(r[4], r[5]),
        lo(r[6], r[7]),
        hi(r[6], r[7]),
    ];
    let u = [
        lo2(s[0], s[2]),
        hi2(s[0], s[2]),
        lo2(s[1], s[3]),
        hi2(s[1], s[3]),
        lo2(s[4], s[6]),
        hi2(s[4], s[6]),
        lo2(s[5], s[7]),
        hi2(s[5], s[7]),
    ];
    [
        lo4(u[0], u[4]),
        hi4(u[0], u[4]),
        lo4(u[1], u[5]),
        hi4(u[1], u[5]),
        lo4(u[2], u[6]),
        hi4(u[2], u[6]),
        lo4(u[3], u[7]),
        hi4(u[3], u[7]),
    ]
}

/// One block of a step-5 search: the first pulse on track `ipos[0]` in 8 lanes (`i0 = ipos[0]
/// + 5l`), then lane by lane the C's depth-first loops for pulses 1..P on tracks `ipos[1..]`
/// (`lanes_next`), with the C's shifts: pulse 1 from `rr[i0][i0] << 14` with 14/15, pulse 2
/// from `alp << 14` with 12/13, pulse 3 from `alp << 16` with 12/13. Returns the last pulse's
/// lanes and every pulse's position index per lane (`ix[0][l] = l`). Reads `rr[i0][t + 5n]` as
/// `rr[t + 5n][i0]` (one row per n, lanes = i0): `cor_h_tm` writes `rr` symmetric, as the C's.
#[inline(always)]
pub(crate) fn step5_block<const P: usize>(
    ipos: &[usize; P],
    dt: &[[i16; 8]; NB_TRACK],
    rr: &RrTm,
) -> (Lanes, [[i16; 8]; P]) {
    const { assert!(P >= 2, "a block places at least one pulse after the first") };
    const {
        assert!(
            P <= 4,
            "STEP5_SHIFTS' proof covers three pulses after the first"
        )
    };
    let t0 = ipos[0];
    let mut ix = [[0i16; 8]; P];
    ix[0] = core::array::from_fn(|l| l as i16);
    let mut o = Lanes {
        sq: [0; 8],
        ps: dt[t0],
        alp: [0; 8],
        ix: [0; 8],
    };
    let mut alp0: [i32; 8] = core::array::from_fn(|l| (rr.diag[t0][l] as i32) << STEP5_SHIFTS[0].2);
    for p in 1..P {
        let t = ipos[p];
        // x[n][l] = rr[i0][t + 5n] + the earlier pulses' rows, transposed.
        let mut x: [[i32; 8]; 8] =
            core::array::from_fn(|n| core::array::from_fn(|l| rr.b[t][t0][n][l] as i32));
        for q in 1..p {
            let rows = core::array::from_fn(|l| rr.b[ipos[q]][t][ix[q][l] as usize]);
            let c = crate::simd::transpose8(rows);
            for n in 0..8 {
                for l in 0..8 {
                    x[n][l] += c[n][l] as i32;
                }
            }
        }
        let (sd, sc, sh) = STEP5_SHIFTS[p];
        o = lanes_next(&o.ps, &alp0, &dt[t], &rr.diag[t], &x, sd, sc);
        ix[p] = o.ix;
        alp0 = core::array::from_fn(|l| (o.alp[l] as i32) << sh);
    }
    (o, ix)
}

/// `dn` track-major: `[t][n] = dn[t + 5n]` (`step5_block`'s rows).
#[inline(always)]
pub(crate) fn dn_tm(dn: &[i16; L_CODE]) -> [[i16; 8]; NB_TRACK] {
    core::array::from_fn(|t| core::array::from_fn(|n| dn[t + NB_TRACK * n]))
}

/// The C's outer compare of the step-5 searches: whether `(sq, alp)` beats `(psk, alpk)`.
#[inline(always)]
fn beats(psk: i16, alpk: i16, sq: i16, alp: i16) -> bool {
    // L_msu(L_mult(alpk, sq), psk, alp) > 0: a saturated difference keeps its sign. (L_mult
    // saturates only for -32768^2; sq, psk >= -1.)
    alpk as i32 * sq as i32 > psk as i32 * alp as i32
}

/// `search_2i40` (MR475/MR515) and MR59's search (`SAT`: saturating `ps`): 2 pulses, one block
/// per track pair `[t0, t1]`, the 8 first pulses in lanes (`step5_block`), then the C's chain
/// over the lanes' winners. Portable body of the dispatched `search_2i40`/`search_2i40_sat`.
#[inline(always)]
pub(crate) fn search_2i40_portable<const SAT: bool>(
    dn: &[i16; L_CODE],
    rr: &RrTm,
    pairs: &[[usize; 2]],
    _ov: &mut Flag,
) -> [usize; 2] {
    let dt = dn_tm(dn);
    let mut codvec = [0usize, 1];
    let (mut psk, mut alpk) = (-1i16, 1i16);
    // All blocks first (independent, latency-bound lane chains overlap), then the C's chain.
    let mut outs = [(Lanes::default(), [[0i16; 8]; 2]); 8];
    assert!(
        pairs.len() <= outs.len(),
        "search_2i40: at most 8 track pairs (MR59)"
    );
    for (out, ipos) in outs.iter_mut().zip(pairs) {
        *out = step5_block(ipos, &dt, rr);
    }
    for (ipos, (o, ix)) in pairs.iter().zip(&outs) {
        for l in 0..8 {
            if beats(psk, alpk, o.sq[l], o.alp[l]) {
                (psk, alpk) = (o.sq[l], o.alp[l]);
                codvec = [ipos[0] + STEP * l, ipos[1] + STEP * ix[1][l] as usize];
            }
        }
    }
    codvec
}

/// `search_3i40` (MR67): 3 pulses, first pulses with `dn2 < 0` skipped (as `search_2i40`).
#[inline(always)]
pub(crate) fn search_3i40_portable(
    dn: &[i16; L_CODE],
    dn2: &[i16; L_CODE],
    rr: &RrTm,
    ov: &mut Flag,
) -> [usize; 3] {
    let dt = dn_tm(dn);
    let mut codvec = [0usize, 1, 2];
    let (mut psk, mut alpk) = (-1i16, 1i16);
    // The 12 blocks are independent until the C's compare chain over their lanes: all made
    // first (their latency-bound lane chains overlap), then the chain in the C's order.
    let mut blocks = [([0usize; 3], Lanes::default(), [[0i16; 8]; 3]); 12];
    let mut nb = 0;
    for track1 in (1..4).step_by(2) {
        for track2 in (2..5).step_by(2) {
            let mut ipos = [0usize, track1, track2];
            for _ in 0..3 {
                let (o, ix) = step5_block(&ipos, &dt, rr);
                blocks[nb] = (ipos, o, ix);
                nb += 1;
                ipos.rotate_right(1);
            }
        }
    }
    for (ipos, o, ix) in &blocks {
        let (ipos, o, ix) = (*ipos, o, ix);
        {
            {
                for l in 0..8 {
                    if dn2[ipos[0] + STEP * l] < 0 {
                        continue;
                    }
                    let (sq, alp) = (o.sq[l], o.alp[l]);
                    let s = l_mult(alpk, sq, ov);
                    if l_msu(s, psk, alp, ov) > 0 {
                        (psk, alpk) = (sq, alp);
                        codvec = core::array::from_fn(|p| ipos[p] + STEP * ix[p][l] as usize);
                    }
                }
            }
        }
    }
    codvec
}

/// Up to 4 first pulses of one step-5 block: its track order and the first pulses' indices
/// on track `ipos[0]` (`v[..n]`, increasing).
#[derive(Clone, Copy)]
pub(crate) struct Half<const P: usize> {
    pub(crate) ipos: [usize; P],
    pub(crate) v: [usize; 4],
    pub(crate) n: usize,
}

/// `step5_block` for two half blocks at once: lanes `0..4` are `h[0]`'s first pulses, `4..8`
/// `h[1]`'s, each half with its own track order (so `d`, `diag` and every pulse's cross row
/// are per lane: rows gathered per lane and transposed, including the first pulse's). Lanes
/// past a half's `n` hold arbitrary valid positions; the caller ignores them.
#[inline(always)]
pub(crate) fn step5_halves<const P: usize>(
    h: &[Half<P>; 2],
    dt: &[[i16; 8]; NB_TRACK],
    rr: &RrTm,
) -> (Lanes, [[i16; 8]; P]) {
    const {
        assert!(
            P <= 4,
            "STEP5_SHIFTS' proof covers three pulses after the first"
        )
    };
    let mut ix = [[0i16; 8]; P];
    ix[0] = core::array::from_fn(|l| (h[l / 4].v[l % 4] & 7) as i16);
    let (d0, g0) = (
        [&dt[h[0].ipos[0]], &dt[h[1].ipos[0]]],
        [&rr.diag[h[0].ipos[0]], &rr.diag[h[1].ipos[0]]],
    );
    let mut o = Lanes {
        sq: [0; 8],
        ps: core::array::from_fn(|l| d0[l / 4][ix[0][l] as usize & 7]),
        alp: [0; 8],
        ix: [0; 8],
    };
    let mut alp0: [i32; 8] =
        core::array::from_fn(|l| (g0[l / 4][ix[0][l] as usize & 7] as i32) << STEP5_SHIFTS[0].2);
    for p in 1..P {
        let t = [h[0].ipos[p], h[1].ipos[p]];
        let mut x = [[0i32; 8]; 8];
        for q in 0..p {
            let b = [&rr.b[h[0].ipos[q]][t[0]], &rr.b[h[1].ipos[q]][t[1]]];
            let rows = core::array::from_fn(|l| b[l / 4][ix[q][l] as usize & 7]);
            let c = crate::simd::transpose8(rows);
            for n in 0..8 {
                for l in 0..8 {
                    x[n][l] += c[n][l] as i32;
                }
            }
        }
        let d = [&dt[t[0]], &dt[t[1]]];
        let diag = [&rr.diag[t[0]], &rr.diag[t[1]]];
        let (sd, sc, sh) = STEP5_SHIFTS[p];
        o = lanes_next(&o.ps, &alp0, &d, &diag, &x, sd, sc);
        ix[p] = o.ix;
        alp0 = core::array::from_fn(|l| (o.alp[l] as i32) << sh);
    }
    (o, ix)
}

/// `search_4i40` (MR74, MR795): 4 pulses, first pulses with `dn2 < 0` skipped. `set_sign`
/// leaves 4 to 8 of a track's 8 positions (4 unless a track has 32767s), so the blocks' valid
/// first pulses are packed: each block splits into halves of up to 4, two halves per
/// `step5_halves` (normally 4 passes instead of 8). Halves keep the C's order, and so does
/// the outer compare chain over their lanes.
#[inline(always)]
pub(crate) fn search_4i40_portable(
    dn: &[i16; L_CODE],
    dn2: &[i16; L_CODE],
    rr: &RrTm,
    _ov: &mut Flag,
) -> [usize; 4] {
    let dt = dn_tm(dn);
    let none = Half {
        ipos: [0usize, 1, 2, 3],
        v: [0; 4],
        n: 0,
    };
    let mut halves = [none; 16];
    let mut nh = 0;
    for track in 3..5 {
        let mut ipos = [0usize, 1, 2, track];
        for _ in 0..4 {
            let mut cur = Half { ipos, ..none };
            for l in 0..8 {
                if dn2[ipos[0] + STEP * l] >= 0 {
                    cur.v[cur.n] = l;
                    cur.n += 1;
                    if cur.n == 4 {
                        halves[nh] = cur;
                        nh += 1;
                        cur.n = 0;
                    }
                }
            }
            if cur.n > 0 {
                halves[nh] = cur;
                nh += 1;
            }
            ipos.rotate_right(1);
        }
    }
    let mut codvec = [0usize, 1, 2, 3];
    let (mut psk, mut alpk) = (-1i16, 1i16);
    // All passes first (independent, latency-bound lane chains overlap), then the C's chain.
    let pairs = halves[..nh.next_multiple_of(2)].as_chunks::<2>().0;
    let mut outs = [(Lanes::default(), [[0i16; 8]; 4]); 8];
    for (out, pair) in outs.iter_mut().zip(pairs) {
        *out = step5_halves(pair, &dt, rr);
    }
    for (pair, (o, ix)) in pairs.iter().zip(&outs) {
        for (k, h) in pair.iter().enumerate() {
            for l in 4 * k..4 * k + h.n {
                if beats(psk, alpk, o.sq[l], o.alp[l]) {
                    (psk, alpk) = (o.sq[l], o.alp[l]);
                    codvec = core::array::from_fn(|p| h.ipos[p] + STEP * ix[p][l] as usize);
                }
            }
        }
    }
    codvec
}

/// `y[i] = round(sum_k h[i - codvec[k]] * sign[k])` (L_mac accumulation), the filtered codevector.
/// Blocks of 8 outputs in lanes over a zero-padded `h`, each `L_mac` exact (`L_mult` then a
/// saturating `L_add`, as selects). The C's first term is `L_mult` in some modes, `L_mac(0, ..)`
/// in others: the same function.
#[inline(always)]
pub(crate) fn filtered_code_portable(
    h: &[i16],
    codvec: &[usize],
    sign: &[i16],
    y: &mut [i16],
    ov: &mut Flag,
) {
    let mut hp = [0i16; 2 * L_CODE];
    hp[L_CODE..].copy_from_slice(&h[..L_CODE]);
    let mut sat = [0i32; 8];
    for (b, yb) in y[..L_CODE].as_chunks_mut::<8>().0.iter_mut().enumerate() {
        let mut s = [0i32; 8];
        for (&c, &sg) in codvec.iter().zip(sign) {
            let w: &[i16; 8] = hp[L_CODE + 8 * b - c..][..8].try_into().unwrap();
            for l in 0..8 {
                let p = w[l] as i32 * sg as i32;
                let big = p == 0x4000_0000; // L_mult saturates
                let q = if big { MAX_32 } else { p << 1 };
                let sum = s[l].wrapping_add(q);
                let o = ((s[l] ^ sum) & (q ^ sum)) < 0;
                sat[l] |= -((o | big) as i32);
                s[l] = if o { (s[l] >> 31) ^ MAX_32 } else { sum };
            }
        }
        for l in 0..8 {
            let r = s[l].wrapping_add(0x8000);
            let o = (!s[l] & r) < 0;
            sat[l] |= -(o as i32);
            yb[l] = (if o { MAX_32 } else { r } >> 16) as i16;
        }
    }
    *ov |= sat.iter().any(|&f| f != 0);
}

/// `code_2i40_9bits` (MR475, MR515): 2 pulses; returns the position index, `sign` the signs.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn code_2i40_9bits(
    sub_nr: usize,
    x: &[i16],
    h: &mut [i16],
    t0: i16,
    pitch_sharp: i16,
    code: &mut [i16],
    y: &mut [i16],
    sign: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let l_temp = (pitch_sharp as i32) << 1;
    let sharp = if l_temp != l_temp as i16 as i32 {
        *ov = true;
        if pitch_sharp > 0 { MAX_16 } else { MIN_16 }
    } else {
        l_temp as i16
    };
    if (t0 as usize) < L_CODE {
        sharpen(h, t0, sharp, ov);
    }
    let mut dn = [0i16; L_CODE];
    let mut dn2 = [0i16; L_CODE];
    let mut dn_sign = [0i16; L_CODE];
    let mut rr = RrTm::ZERO;
    cor_h_x(h, x, &mut dn, 1, ov);
    set_sign(&mut dn, &mut dn_sign, &mut dn2, 8);
    let pairs = [0, 1].map(|t| {
        let i = (sub_nr << 1) + (t << 3);
        [STARTPOS[i] as usize, STARTPOS[i + 1] as usize]
    });
    cor_h_blocks(h, &dn_sign, &mut rr, &pairs.map(|[t0, t1]| [t1, t0]), ov);
    let codvec = search_2i40(&dn, &rr, &pairs, ov);
    // build_code
    let pt = &TRACKTABLE[sub_nr * 5..];
    code[..L_CODE].fill(0);
    let mut indx: i16 = 0;
    let mut rsign: i16 = 0;
    let mut sgn = [0i16; 2];
    for k in 0..2 {
        let i = codvec[k];
        let mut index = ((i as i32 * 6554) >> 15) as i16;
        let track = i as i16 - 5 * index;
        let first = pt[track as usize];
        let track = if k == 0 {
            if first != 0 {
                index += 64;
            }
            0
        } else {
            index <<= 3;
            1
        };
        if dn_sign[i] > 0 {
            code[i] = 8191;
            sgn[k] = 32767;
            rsign += 1 << track;
        } else {
            code[i] = -8192;
            sgn[k] = -32768;
        }
        indx = indx.wrapping_add(index);
    }
    *sign = rsign;
    filtered_code(h, &codvec, &sgn, y, ov);
    if (t0 as usize) < L_CODE {
        sharpen(code, t0, sharp, ov);
    }
    indx
}

/// `code_2i40_11bits` (MR59): 2 pulses.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn code_2i40_11bits(
    x: &[i16],
    h: &mut [i16],
    t0: i16,
    pitch_sharp: i16,
    code: &mut [i16],
    y: &mut [i16],
    sign: &mut i16,
    ov: &mut Flag,
) -> i16 {
    const START_POS1: [usize; 2] = [1, 3];
    const START_POS2: [usize; 4] = [0, 1, 2, 4];
    let sharp = ((pitch_sharp as i32) << 1) as i16;
    if (t0 as usize) < L_CODE {
        sharpen(h, t0, sharp, ov);
    }
    let mut dn = [0i16; L_CODE];
    let mut dn2 = [0i16; L_CODE];
    let mut dn_sign = [0i16; L_CODE];
    let mut rr = RrTm::ZERO;
    cor_h_x(h, x, &mut dn, 1, ov);
    set_sign(&mut dn, &mut dn_sign, &mut dn2, 8);
    cor_h_tm(h, &dn_sign, &mut rr, ov);
    let pairs: [[usize; 2]; 8] = core::array::from_fn(|k| [START_POS1[k / 4], START_POS2[k % 4]]);
    let codvec = search_2i40_sat(&dn, &rr, &pairs, ov);
    code[..L_CODE].fill(0);
    let mut indx: i16 = 0;
    let mut rsign: i16 = 0;
    let mut sgn = [0i16; 2];
    for k in 0..2 {
        let i = codvec[k];
        let mut index = ((i as i32 * 6554) >> 15) as i16;
        let track_raw = i as i16 - (((index << 3) + (index << 1)) >> 1);
        let mut track = track_raw;
        match track_raw {
            0 => {
                track = 1;
                index <<= 6;
            }
            1 => {
                if k == 0 {
                    track = 0;
                    index <<= 1;
                } else {
                    track = 1;
                    index = (index << 6) + 16;
                }
            }
            2 => {
                track = 1;
                index = (index << 6) + 32;
            }
            3 => {
                track = 0;
                index = (index << 1) + 1;
            }
            4 => {
                track = 1;
                index = (index << 6) + 48;
            }
            _ => {}
        }
        if dn_sign[i] > 0 {
            code[i] = 8191;
            sgn[k] = 32767;
            let t = shl(1, track, ov);
            rsign = add(rsign, t, ov);
        } else {
            code[i] = -8192;
            sgn[k] = -32768;
        }
        indx = add(indx, index, ov);
    }
    *sign = rsign;
    filtered_code(h, &codvec, &sgn, y, ov);
    if (t0 as usize) < L_CODE {
        sharpen(code, t0, sharp, ov);
    }
    indx
}

/// `code_3i40_14bits` (MR67): 3 pulses.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn code_3i40_14bits(
    x: &[i16],
    h: &mut [i16],
    t0: i16,
    pitch_sharp: i16,
    code: &mut [i16],
    y: &mut [i16],
    sign: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let sharp = ((pitch_sharp as i32) << 1) as i16;
    if (t0 as usize) < L_CODE {
        sharpen(h, t0, sharp, ov);
    }
    let mut dn = [0i16; L_CODE];
    let mut dn2 = [0i16; L_CODE];
    let mut dn_sign = [0i16; L_CODE];
    let mut rr = RrTm::ZERO;
    cor_h_x(h, x, &mut dn, 1, ov);
    set_sign(&mut dn, &mut dn_sign, &mut dn2, 6);
    cor_h_tm(h, &dn_sign, &mut rr, ov);
    let codvec = search_3i40(&dn, &dn2, &rr, ov);
    code[..L_CODE].fill(0);
    let mut indx: i16 = 0;
    let mut rsign: i16 = 0;
    let mut sgn = [0i16; 3];
    for k in 0..3 {
        let i = codvec[k];
        let mut index = ((i as i32 * 6554) >> 15) as i16;
        let mut track = i as i16 - index * 5;
        match track {
            1 => index <<= 4,
            2 => index <<= 8,
            3 => {
                track = 1;
                index = (index << 4) + 8;
            }
            4 => {
                track = 2;
                index = (index << 8) + 128;
            }
            _ => {}
        }
        if dn_sign[i] > 0 {
            code[i] = 8191;
            sgn[k] = 32767;
            rsign = rsign.wrapping_add(1 << track);
        } else {
            code[i] = -8192;
            sgn[k] = -32768;
        }
        indx = indx.wrapping_add(index);
    }
    *sign = rsign;
    filtered_code(h, &codvec, &sgn, y, ov);
    if (t0 as usize) < L_CODE {
        sharpen(code, t0, sharp, ov);
    }
    indx
}

/// `code_4i40_17bits` (MR74, MR795): 4 pulses.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn code_4i40_17bits(
    x: &[i16],
    h: &mut [i16],
    t0: i16,
    pitch_sharp: i16,
    code: &mut [i16],
    y: &mut [i16],
    sign: &mut i16,
    ov: &mut Flag,
) -> i16 {
    let sharp = ((pitch_sharp as i32) << 1) as i16;
    if (t0 as usize) < L_CODE {
        sharpen(h, t0, sharp, ov);
    }
    let mut dn = [0i16; L_CODE];
    let mut dn2 = [0i16; L_CODE];
    let mut dn_sign = [0i16; L_CODE];
    let mut rr = RrTm::ZERO;
    cor_h_x(h, x, &mut dn, 1, ov);
    set_sign(&mut dn, &mut dn_sign, &mut dn2, 4);
    cor_h_tm(h, &dn_sign, &mut rr, ov);
    let codvec = search_4i40(&dn, &dn2, &rr, ov);
    code[..L_CODE].fill(0);
    let mut indx: i16 = 0;
    let mut rsign: i16 = 0;
    let mut sgn = [0i16; 4];
    for k in 0..4 {
        let i = codvec[k];
        let index0 = ((i as i32 * 6554) >> 15) as i16;
        let mut track = i as i16 - index0 * 5;
        let mut index = GRAY[index0 as usize];
        match track {
            1 => index <<= 3,
            2 => index <<= 6,
            3 => index <<= 10,
            4 => {
                track = 3;
                index = (index << 10) + 512;
            }
            _ => {}
        }
        if dn_sign[i] > 0 {
            code[i] = 8191;
            sgn[k] = 32767;
            rsign = rsign.wrapping_add(1 << track);
        } else {
            code[i] = -8192;
            sgn[k] = -32768;
        }
        indx = indx.wrapping_add(index);
    }
    *sign = rsign;
    filtered_code(h, &codvec, &sgn, y, ov);
    if t0 < L_CODE as i16 {
        sharpen(code, t0, sharp, ov);
    }
    indx
}

/// One pulse pair's N × N candidates (`[j][i]`, rows of `W` lanes; lanes `N..W` are padding):
/// correlation `ps`, its square `sq` and the energy `al` as the C compares it (the first pair's
/// 32-bit `alp2`, later pairs' `alp_16`). Built only for `chain` (the fallback, and tests).
pub(crate) struct Cands<const N: usize, const W: usize> {
    pub(crate) ps: [[i16; W]; N],
    pub(crate) sq: [[i16; W]; N],
    pub(crate) al: [[i32; W]; N],
}

/// The C's compare chain over `c` in j-major order: (j, i, sq, alp, ps) of the winner.
pub(crate) fn chain<const N: usize, const W: usize>(
    c: &Cands<N, W>,
) -> (usize, usize, i16, i16, i16) {
    let (mut sq, mut alp, mut ps, mut best) = (-1i16, 1i16, 0i16, (0, 0));
    for j in 0..N {
        for i in 0..N {
            let (sq2, alp2) = (c.sq[j][i], c.al[j][i]);
            if sq2 as i32 * alp as i32 > (sq as i32).wrapping_mul(alp2) {
                (sq, alp, ps, best) = (sq2, alp2 as i16, c.ps[j][i], (j, i));
            }
        }
    }
    (best.0, best.1, sq, alp, ps)
}

/// One pulse pair of `search_10and8i40`, first pulse at `ta + T * j` (rows), second at
/// `tb + T * i` (lanes): candidate (j, i) has `ps = d1[j] + tps[i]` and, if `wide` (the first
/// pair), `alp2 = (((base[j] + rr[j][i]) >> 2) + tal[i]) >> 1`, else
/// `alp_16 = (base[j] + (rr[j][i] << sh) + tal[i]) >> 16` (wrapping, as the C).
pub(crate) struct Pair<'a, const N: usize, const W: usize> {
    pub(crate) d1: &'a [i16; W],
    pub(crate) tps: [i16; W],
    pub(crate) base: [i32; W],
    pub(crate) tal: [i32; W],
    pub(crate) rr: &'a [[i16; W]; N],
    pub(crate) sh: u32,
    pub(crate) wide: bool,
}

impl<const N: usize, const W: usize> Pair<'_, N, W> {
    /// Row j's `ps`, `sq` and energy, whole rows (plain loops, the branch outside: the form
    /// that vectorises and inlines).
    #[inline(always)]
    fn row(&self, j: usize) -> ([i16; W], [i16; W], [i32; W]) {
        let (ps1, rj, b) = (self.d1[j], &self.rr[j], self.base[j]);
        let (mut ps, mut sq, mut al) = ([0i16; W], [0i16; W], [0i32; W]);
        for i in 0..W {
            ps[i] = ps1.wrapping_add(self.tps[i]);
            sq[i] = ((ps[i] as i32 * ps[i] as i32) >> 15) as i16;
        }
        if self.wide {
            for i in 0..W {
                al[i] = (((b + rj[i] as i32) >> 2) + self.tal[i]) >> 1;
            }
        } else {
            for i in 0..W {
                let alp2 = b.wrapping_add((rj[i] as i32) << self.sh);
                al[i] = (alp2.wrapping_add(self.tal[i]) >> 16) as i16 as i32;
            }
        }
        (ps, sq, al)
    }

    /// All candidates, for `chain`.
    pub(crate) fn cands(&self) -> Cands<N, W> {
        let mut c = Cands {
            ps: [[0; W]; N],
            sq: [[0; W]; N],
            al: [[0; W]; N],
        };
        for j in 0..N {
            (c.ps[j], c.sq[j], c.al[j]) = self.row(j);
        }
        c
    }

    /// `chain(&self.cands())`, faster: while every energy is in 1..=32767 and every `sq` in
    /// 0..=32767, `sq2 * alp > sq * alp2` is the ratio order `sq2 / alp2 > sq / alp` (no wrap,
    /// no truncation), so the first maximum is found per lane as the rows are made
    /// (independent chains over j, nothing stored), then across lanes, ties to the lowest
    /// (j, i). Any lane outside that range: `chain`. Padding lanes are checked too: read only
    /// up to N, LLVM drops them and splits the rows 8 + 2 (MR102: 1220 vs 920 ns). They can
    /// force a needless (still exact) fallback; on the vectors no pick falls back at all.
    /// They are never chosen. `L32`: the AVX2 form (`best_i32`), else `best_i16`.
    #[inline(always)]
    pub(crate) fn best<const L32: bool>(&self) -> (usize, usize, i16, i16, i16) {
        if L32 {
            self.best_i32()
        } else {
            self.best_i16()
        }
    }

    /// `best`'s portable form: the per-lane chain in i16 lanes (`bad` per element), the
    /// cross-lane pick a compare loop. SSE2 has no `pmulld`/`blendv`/`pcmpgtq`, and its 16
    /// xmm registers hold 8 i32 lanes in two: `best_i32` there was 1.07× this chain and the
    /// f64 pick 1.05–1.18× the loop.
    #[inline(always)]
    fn best_i16(&self) -> (usize, usize, i16, i16, i16) {
        let mut bad = [0i16; W];
        let mut mark = |sq: &[i16; W], al: &[i32; W]| {
            for i in 0..W {
                bad[i] |= -(((sq[i] < 0) | (al[i] < 1) | (al[i] > 32767)) as i16);
            }
        };
        let (_, mut bsq, al) = self.row(0);
        mark(&bsq, &al);
        let (mut bal, mut bj): ([i16; W], _) = (core::array::from_fn(|i| al[i] as i16), [0i16; W]);
        for j in 1..N {
            let (_, sq, al) = self.row(j);
            mark(&sq, &al);
            for i in 0..W {
                let al = al[i] as i16;
                let m = -((sq[i] as i32 * bal[i] as i32 > bsq[i] as i32 * al as i32) as i16);
                bsq[i] = (sq[i] & m) | (bsq[i] & !m);
                bal[i] = (al & m) | (bal[i] & !m);
                bj[i] = (j as i16 & m) | (bj[i] & !m);
            }
        }
        if bad.iter().fold(0, |a, &b| a | b) != 0 {
            return chain(&self.cands());
        }
        let mut b = 0;
        for i in 1..N {
            let (l, r) = (bsq[i] as i32 * bal[b] as i32, bsq[b] as i32 * bal[i] as i32);
            if l > r || (l == r && bj[i] < bj[b]) {
                b = i;
            }
        }
        let j = bj[b] as usize;
        (j, b, bsq[b], bal[b], self.d1[j].wrapping_add(self.tps[b]))
    }

    /// Row j in i32 lanes for `best_i32`: `sq` unclamped ((ps * ps) >> 15: 32768 for
    /// ps = -32768, where the C's i16 has -32768) and the energy.
    #[inline(always)]
    fn row32(&self, j: usize) -> ([i32; W], [i32; W]) {
        let (ps1, rj, b) = (self.d1[j], &self.rr[j], self.base[j]);
        let (mut sq, mut al) = ([0i32; W], [0i32; W]);
        for i in 0..W {
            let ps = ps1.wrapping_add(self.tps[i]) as i32;
            sq[i] = (ps * ps) >> 15;
        }
        if self.wide {
            for i in 0..W {
                al[i] = (((b + rj[i] as i32) >> 2) + self.tal[i]) >> 1;
            }
        } else {
            for i in 0..W {
                al[i] = b
                    .wrapping_add((rj[i] as i32) << self.sh)
                    .wrapping_add(self.tal[i])
                    >> 16;
            }
        }
        (sq, al)
    }

    /// `best`'s AVX2 form (in 128-bit lanes the i32 vectors spill and SSE2 has no `pcmpgtq`:
    /// slower there). The chain is i32 lanes (no i32 -> i16 mask packing) holding the values
    /// as made; each product masks its operands to 15 bits, which LLVM keeps (it can't prove
    /// them redundant), so the products are `pmaddwd`, not `vpmulld` (stored masked, the
    /// masks are dropped as redundant and the 10-cycle multiply bounds the chain). The range
    /// check is one OR per lane: `sq | al | (al - 1)` has a bit above 14 set exactly when `sq`
    /// is outside 0..=32767 or `al` outside 1..=32767. Across lanes, no compare chain: sq / al
    /// as an f64 orders exactly (distinct ratios of 15-bit values differ by a relative
    /// 1 / (sq * al') > 2^-30, f64 rounds within 2^-53, equal ones round equal), positive
    /// doubles order as their bits, so the low 16 bits (2^-36 relative) carry the tie rule
    /// (lowest j, then lane) and one i64 max finds the C's winner.
    #[inline(always)]
    fn best_i32(&self) -> (usize, usize, i16, i16, i16) {
        const M: i32 = 0x7fff;
        let (mut bsq, mut bal) = self.row32(0);
        let mut bj = [0i32; W];
        let mut lim: [i32; W] = core::array::from_fn(|i| bsq[i] | bal[i] | bal[i].wrapping_sub(1));
        for j in 1..N {
            let (sq, al) = self.row32(j);
            for i in 0..W {
                lim[i] |= sq[i] | al[i] | al[i].wrapping_sub(1);
                let m = -(((sq[i] & M) * (bal[i] & M) > (bsq[i] & M) * (al[i] & M)) as i32);
                bsq[i] = (sq[i] & m) | (bsq[i] & !m);
                bal[i] = (al[i] & m) | (bal[i] & !m);
                bj[i] = (j as i32 & m) | (bj[i] & !m);
            }
        }
        if lim.iter().fold(0, |a, &b| a | b) & !M != 0 {
            return chain(&self.cands());
        }
        let mut key = [0i64; W];
        for i in 0..W {
            let r = (bsq[i] as f64 / bal[i] as f64).to_bits() as i64;
            let k = (r & !0xffff) | (0xffff - (bj[i] as i64 * W as i64 + i as i64));
            key[i] = if i < N { k } else { 0 };
        }
        let mut w = W.next_power_of_two();
        while w > 1 {
            w /= 2;
            for i in 0..w {
                let o = if i + w < W { key[i + w] } else { 0 };
                key[i] = key[i].max(o);
            }
        }
        let t = 0xffff - (key[0] & 0xffff) as usize;
        let (j, b) = (t / W, t % W);
        (
            j,
            b,
            bsq[b] as i16,
            bal[b] as i16,
            self.d1[j].wrapping_add(self.tps[b]),
        )
    }
}

/// MR122's `search_10and8i40` (10 pulses, 5 tracks of 8). Portable body of the dispatched
/// `search_10i40` (see `search_pairs`).
#[inline(always)]
pub(crate) fn search_10i40_portable<const L32: bool>(
    dn: &[i16; L_CODE],
    rr: &RrTm,
    ipos: &mut [i16; 10],
    pos_max: &[i16; NB_TRACK],
    codvec: &mut [usize; 10],
) {
    search_pairs::<5, 8, 8, 10, L32>(dn, rr, ipos, pos_max, codvec);
}

/// MR102's `cor_h` and `search_10and8i40` (8 pulses, 4 tracks of 10), one kernel so the
/// layout's row width can differ per path. Portable body of the dispatched `cor_h_search_8i40`
/// (W = 10 with the ±32767 sign form; AVX2: 16, the multiply form).
#[allow(clippy::too_many_arguments)] // cor_h's and the search's inputs
#[inline(always)]
pub(crate) fn cor_h_search_8i40_portable<const W: usize, const NARROW: bool, const L32: bool>(
    h: &[i16],
    sign: &[i16; L_CODE],
    dn: &[i16; L_CODE],
    ipos: &mut [i16; 8],
    pos_max: &[i16; NB_TRACK_MR102],
    codvec: &mut [usize; 8],
    ov: &mut Flag,
) {
    let mut rr = RrTm4::<W>::ZERO;
    cor_h_tm_portable::<NB_TRACK_MR102, 10, W, NARROW>(h, sign, &mut rr, ov);
    search_8i40_portable::<W, L32>(dn, &rr, ipos, pos_max, codvec);
}

/// MR102's `search_10and8i40` (see `search_pairs`).
#[inline(always)]
pub(crate) fn search_8i40_portable<const W: usize, const L32: bool>(
    dn: &[i16; L_CODE],
    rr: &RrTm4<W>,
    ipos: &mut [i16; 8],
    pos_max: &[i16; NB_TRACK_MR102],
    codvec: &mut [usize; 8],
) {
    search_pairs::<4, 10, W, 8, L32>(dn, rr, ipos, pos_max, codvec);
}

/// The Q15 weights of `search_10and8i40`'s pair `p` (1..=4, pulses 2p and 2p + 1): `rrv`'s
/// diagonal and cross terms, the first pulse's diagonal and cross terms, the second pulse's
/// `rrv` and cross term; all powers of two (`pair_safe` shifts by them).
const PAIR_W: [[i64; 6]; 5] = [
    [0; 6],
    [4096, 8192, 2048, 4096, 16384, 4096],
    [4096, 8192, 1024, 2048, 8192, 2048],
    [2048, 4096, 512, 1024, 8192, 1024],
    [2048, 4096, 256, 512, 4096, 512],
];

/// Whether pair `p`'s L_mac chains stay below `MAX_32 - 0x8000` for every candidate (so
/// neither they nor their `round` saturate): bounded by absolute values, from `alp0_abs` (the
/// energy carried in), the earlier pulses `prior`, and the pair's tracks `ta`, `tb`. Reads the
/// track-major rows directly (`rr.b[q % T][t][q / T]` is position q against track t).
/// First tier: the same bound from per-block maxima (`bm[t][u]` = max |rr| over block (t, u),
/// `dm[t]` over track t's diagonal, made once per search): each term bounds its per-position
/// counterpart, so a pass here is a pass of the precise check below.
fn pair_safe<const T: usize, const N: usize, const W: usize>(
    rr: &Tm<T, N, W>,
    (bm, dm): &([[i64; T]; T], [i64; T]),
    p: usize,
    alp0_abs: i64,
    prior: &[usize],
    ta: usize,
    tb: usize,
) -> bool {
    let [rd, rc, ad, ac, ar, a2c] = PAIR_W[p];
    let lim = MAX_32 as i64 - 0x8000;
    let (sa, sb) = prior
        .iter()
        .fold((0, 0), |(a, b), &q| (a + bm[q % T][ta], b + bm[q % T][tb]));
    let srrv = 2 * (rd * dm[tb] + rc * sb);
    let a2 = alp0_abs + 2 * (ad * dm[ta] + ac * sa) + 2 * ar * ((srrv + 0x8000) >> 16);
    if srrv <= lim && a2 + 2 * a2c * bm[ta][tb] <= lim {
        return true;
    }
    // per position n of a track: sum over the earlier pulses of |rr[q][t + T n]|
    let cross_sum = |t: usize| -> [i32; N] {
        let mut s = [0i32; N];
        for &q in prior {
            let row = &rr.b[q % T][t][q / T];
            for n in 0..N {
                s[n] += (row[n] as i32).abs();
            }
        }
        s
    };
    let (sa, sb) = (cross_sum(ta), cross_sum(tb));
    // The weights are powers of two >= 256: the per-position maxima in i32 lanes as shifts by
    // log2(weight / 256) (< 2^25: 16 x 2^15 + 32 x 8 x 2^15, doubled), exact, scaled back once.
    let w = |v: i64| v.trailing_zeros() - 8;
    let (rd, rc, ad, ac) = (w(rd), w(rc), w(ad), w(ac));
    let (mut srrv, mut a1) = (0i32, 0i32);
    for n in 0..N {
        let db = (rr.diag[tb][n] as i32).abs();
        let da = (rr.diag[ta][n] as i32).abs();
        srrv = srrv.max(((db << rd) + (sb[n] << rc)) << 1);
        a1 = a1.max(((da << ad) + (sa[n] << ac)) << 1);
    }
    // max |c| over the block as one flat reduction (padding lanes are 0; per row and column
    // LLVM gathered columns)
    let (lo, hi) = rr.b[ta][tb]
        .as_flattened()
        .iter()
        .fold((0i16, 0i16), |(l, h), &c| (l.min(c), h.max(c)));
    let cross = (hi as i64).max(-(lo as i64));
    let (srrv, a1) = ((srrv as i64) << 8, (a1 as i64) << 8);
    let rrv = (srrv + 0x8000) >> 16;
    let a2 = alp0_abs + a1 + 2 * ar * rrv + 2 * a2c * cross;
    srrv <= lim && a2 <= lim
}

/// `search_10and8i40` exactly (scalar, TS 26.073's saturating operators): the fallback of
/// `search_pairs` where a chain may saturate. Same loops as the C (opencore `s10_8pf.cpp`).
#[cold]
fn search_pairs_exact<const T: usize, const N: usize, const W: usize, const P: usize>(
    dn: &[i16; L_CODE],
    rr: &Tm<T, N, W>,
    ipos: &mut [i16; P],
    pos_max: &[i16; T],
    codvec: &mut [usize; P],
) {
    let ov = &mut false;
    let r = |a: usize, b: usize| rr.at(a, b);
    let i0 = pos_max[ipos[0] as usize] as usize;
    let (mut psk, mut alpk) = (-1i16, 1i16);
    for (i, c) in codvec.iter_mut().enumerate() {
        *c = i;
    }
    for _ in 1..T {
        let mut pos = [0usize; P];
        pos[0] = i0;
        pos[1] = pos_max[ipos[1] as usize] as usize;
        let mut ps0 = add(dn[pos[0]], dn[pos[1]], ov);
        let mut alp0 = l_mult(r(pos[0], pos[0]), 2048, ov);
        alp0 = l_mac(alp0, r(pos[1], pos[1]), 2048, ov);
        alp0 = l_mac(alp0, r(pos[0], pos[1]), 4096, ov);
        let (mut sq, mut alp) = (-1i16, 1i16);
        for p in 1..T {
            let [rd, rc, ad, ac, ar, a2c] = PAIR_W[p].map(|w| w as i16);
            let (ta, tb) = (ipos[2 * p] as usize, ipos[2 * p + 1] as usize);
            let mut rrv = [0i16; L_CODE];
            for ib in (tb..L_CODE).step_by(T) {
                let mut s = l_mult(r(ib, ib), rd, ov);
                for &q in &pos[..2 * p] {
                    s = l_mac(s, r(q, ib), rc, ov);
                }
                rrv[ib] = round(s, ov);
            }
            let (mut ps, mut best) = (0i16, (ta, tb));
            (sq, alp) = (-1, 1);
            for ia in (ta..L_CODE).step_by(T) {
                let ps1 = add(ps0, dn[ia], ov);
                let mut alp1 = l_mac(alp0, r(ia, ia), ad, ov);
                for &q in &pos[..2 * p] {
                    alp1 = l_mac(alp1, r(q, ia), ac, ov);
                }
                for ib in (tb..L_CODE).step_by(T) {
                    let ps2 = add(ps1, dn[ib], ov);
                    let mut alp2 = l_mac(alp1, rrv[ib], ar, ov);
                    alp2 = l_mac(alp2, r(ia, ib), a2c, ov);
                    let sq2 = mult(ps2, ps2, ov);
                    let alp_16 = round(alp2, ov);
                    if l_msu(l_mult(alp, sq2, ov), sq, alp_16, ov) > 0 {
                        (sq, ps, alp, best) = (sq2, ps2, alp_16, (ia, ib));
                    }
                }
            }
            (pos[2 * p], pos[2 * p + 1]) = best;
            ps0 = ps;
            alp0 = l_mult(alp, 16384, ov);
        }
        if l_msu(l_mult(alpk, sq, ov), psk, alp, ov) > 0 {
            (psk, alpk) = (sq, alp);
            *codvec = pos;
        }
        let first = ipos[1];
        ipos.copy_within(2..P, 1);
        ipos[P - 1] = first;
    }
}

/// `search_10and8i40` (P = 2T pulses on T tracks of N) on the track-major `rr`: per pulse
/// pair, the N × N candidates are made in rows of W lanes and `Pair::best` finds the C's
/// winner. The C's shifts per pair: pair 1 32-bit `alp2`; pair 2 `s >> 2`, 11/12, `<< 14`;
/// pair 3 `s >> 3`, 10/11, `<< 14`; pair 4 (MR122 only) `s >> 3`, 9/10, `<< 13`.
#[inline(always)]
fn search_pairs<const T: usize, const N: usize, const W: usize, const P: usize, const L32: bool>(
    dn: &[i16; L_CODE],
    rr: &Tm<T, N, W>,
    ipos: &mut [i16; P],
    pos_max: &[i16; T],
    codvec: &mut [usize; P],
) {
    const { assert!(P == 2 * T && T * N == L_CODE && N <= W && T <= 5) };
    // TS 26.073 saturates `add`, the L_mac chains and `round`; this kernel's plain sums are its
    // results only where nothing can saturate (`pair_safe`), else the exact scalar search.
    let ipos_in = *ipos;
    let dn_sum: i32 = (0..T)
        .map(|t| {
            (0..N)
                .map(|n| (dn[t + T * n] as i32).abs())
                .max()
                .unwrap_or(0)
        })
        .sum();
    if 2 * dn_sum > 32767 {
        *ipos = ipos_in;
        return search_pairs_exact::<T, N, W, P>(dn, rr, ipos, pos_max, codvec);
    }
    let mut dnt = [[0i16; W]; T];
    for (x, &d) in dn.iter().enumerate() {
        dnt[x % T][x / T] = d;
    }
    let i0 = pos_max[ipos[0] as usize] as usize;
    // The C's T - 1 iterations (one per rotation of `ipos`) are independent until the final
    // comparison, and each is a serial chain of pairs (a pair's rows depend on the winners
    // before it): made level by level, the iterations are independent chains in flight
    // (MR122 AVX2 0.75 -> 0.66 of the old kernel with `best_i32`; flat with `best_i16`).
    // Rows `..T - 1` are used (`[_; T - 1]` isn't a valid generic length).
    let mut ips = [[0i16; P]; T];
    for ip in ips.iter_mut().take(T - 1) {
        *ip = *ipos;
        let pos = ipos[1];
        ipos.copy_within(2..P, 1);
        ipos[P - 1] = pos;
    }
    let mut index = [[0usize; P]; T];
    let (mut alp, mut ps, mut sq) = ([0i16; T], [0i16; T], [0i16; T]);
    let r = |a: usize, b: usize| rr.at(a, b) as i32;
    let maxabs = |v: &[i16]| {
        let (l, h) = v
            .iter()
            .fold((0i16, 0i16), |(l, h), &c| (l.min(c), h.max(c)));
        (h as i64).max(-(l as i64))
    };
    let bounds = (
        core::array::from_fn(|t| core::array::from_fn(|u| maxabs(rr.b[t][u].as_flattened()))),
        core::array::from_fn(|t| maxabs(&rr.diag[t])),
    );
    for ip in ips.iter().take(T - 1) {
        let i1 = pos_max[ip[1] as usize] as usize;
        let alp0_abs = 2
            * (2048 * (r(i0, i0) as i64).abs()
                + 2048 * (r(i1, i1) as i64).abs()
                + 4096 * (r(i0, i1) as i64).abs());
        if !pair_safe(
            rr,
            &bounds,
            1,
            alp0_abs,
            &[i0, i1],
            ip[2] as usize,
            ip[3] as usize,
        ) {
            *ipos = ipos_in;
            return search_pairs_exact::<T, N, W, P>(dn, rr, ipos, pos_max, codvec);
        }
    }
    for (it, ipos) in ips.iter().enumerate().take(T - 1) {
        let i1 = pos_max[ipos[1] as usize] as usize;
        let index = &mut index[it];
        (index[0], index[1]) = (i0, i1);

        // Pulses 2 and 3.
        let (ta, tb) = (ipos[2] as usize, ipos[3] as usize);
        let ps0 = (dn[i0] as i32 + dn[i1] as i32) as i16;
        let mut alp0 = r(i0, i0) << 12;
        alp0 = alp0.wrapping_add(r(i1, i1) << 12);
        alp0 = alp0.wrapping_add(r(i0, i1) << 13);
        alp0 = alp0.wrapping_add(0x8000);
        let s0 = alp0 >> 12;
        let (mut tps, mut tal, mut alp1) = ([0i16; W], [0i32; W], [0i32; W]);
        let (q0b, q1b) = (rr.row(i0, tb), rr.row(i1, tb));
        let (q0a, q1a) = (rr.row(i0, ta), rr.row(i1, ta));
        for n in 0..W {
            let s = (rr.diag[tb][n] as i32 >> 1) + q0b[n] as i32 + q1b[n] as i32;
            tps[n] = ps0.wrapping_add(dnt[tb][n]);
            tal[n] = ((s + 2) >> 2) as i16 as i32;
            alp1[n] = ((s0 + rr.diag[ta][n] as i32) >> 1) + q0a[n] as i32 + q1a[n] as i32;
        }
        let pair = Pair {
            d1: &dnt[ta],
            tps,
            base: alp1,
            tal,
            rr: &rr.b[ta][tb],
            sh: 0,
            wide: true,
        };
        let (j, i);
        (j, i, _, alp[it], ps[it]) = pair.best::<L32>();
        (index[2], index[3]) = (ta + T * j, tb + T * i);
    }

    // Later pulses, in pairs: (first pulse, its C shifts, earlier pulses).
    for p in 2..T {
        for (it, ip) in ips.iter().enumerate().take(T - 1) {
            let alp0_abs = 2 * 16384 * (alp[it] as i64).abs();
            let (ta, tb) = (ip[2 * p] as usize, ip[2 * p + 1] as usize);
            if !pair_safe(rr, &bounds, p, alp0_abs, &index[it][..2 * p], ta, tb) {
                *ipos = ipos_in;
                return search_pairs_exact::<T, N, W, P>(dn, rr, ipos, pos_max, codvec);
            }
        }
        let (sh_s, sh_d) = if p == 2 { (2, 11) } else { (3, 13 - p as u32) };
        let sh_t = if p == 4 { 13 } else { 14 };
        for (it, ipos) in ips.iter().enumerate().take(T - 1) {
            let (ta, tb) = (ipos[2 * p] as usize, ipos[2 * p + 1] as usize);
            let prev = &index[it][..2 * p];
            let alp0 = ((alp[it] as i32) << 15).wrapping_add(0x8000);
            let (mut s, mut a1) = ([0i32; W], [0i32; W]);
            for n in 0..W {
                s[n] = rr.diag[tb][n] as i32 >> 1;
                a1[n] = alp0.wrapping_add((rr.diag[ta][n] as i32) << sh_d);
            }
            for &q in prev {
                let (qb, qa) = (rr.row(q, tb), rr.row(q, ta));
                for n in 0..W {
                    s[n] += qb[n] as i32;
                    a1[n] = a1[n].wrapping_add((qa[n] as i32) << (sh_d + 1));
                }
            }
            let (mut tps, mut tal) = ([0i16; W], [0i32; W]);
            for n in 0..W {
                tps[n] = ps[it].wrapping_add(dnt[tb][n]);
                tal[n] = (((s[n] + (1 << (sh_s - 1))) >> sh_s) as i16 as i32) << sh_t;
            }
            let pair = Pair {
                d1: &dnt[ta],
                tps,
                base: a1,
                tal,
                rr: &rr.b[ta][tb],
                sh: sh_d + 1,
                wide: false,
            };
            let (j, i);
            (j, i, sq[it], alp[it], ps[it]) = pair.best::<L32>();
            index[it][2 * p] = ta + T * j;
            index[it][2 * p + 1] = tb + T * i;
        }
    }
    let (mut psk, mut alpk) = (-1i16, 1i16);
    for (i, c) in codvec.iter_mut().enumerate() {
        *c = i;
    }
    for it in 0..T - 1 {
        if alpk as i32 * sq[it] as i32 > psk as i32 * alp[it] as i32 {
            psk = sq[it];
            alpk = alp[it];
            *codvec = index[it];
        }
    }
}

/// `compress10` of `c8_31pf.cpp`.
pub(super) fn compress10(a: i16, b: i16, c: i16) -> i16 {
    let ia = a >> 1;
    let ib = (b >> 1) * 5;
    let ic = (c >> 1) * 25;
    let indx = (ib + ic + ia) << 3;
    indx + ((a & 1) + ((b & 1) << 1) + ((c & 1) << 2))
}

/// `code_8i40_31bits` (MR102): 8 pulses; writes 7 indices.
pub(super) fn code_8i40_31bits(
    x: &[i16],
    cn: &[i16],
    h: &[i16],
    cod: &mut [i16],
    y: &mut [i16],
    indx: &mut [i16],
    ov: &mut Flag,
) {
    const NB_PULSE: usize = 8;
    let mut ipos = [0i16; NB_PULSE];
    let mut pos_max = [0i16; NB_TRACK_MR102];
    let mut codvec = [0usize; NB_PULSE];
    let mut dn = [0i16; L_CODE];
    let mut sign = [0i16; L_CODE];
    cor_h_x2(h, x, &mut dn, 2, NB_TRACK_MR102, STEP_MR102, ov);
    set_sign12k2(
        &mut dn,
        cn,
        &mut sign,
        &mut pos_max,
        NB_TRACK_MR102,
        &mut ipos,
        STEP_MR102,
        ov,
    );
    cor_h_search_8i40(h, &sign, &dn, &mut ipos, &pos_max, &mut codvec, ov);
    // build_code
    let mut sign_indx = [-1i16; NB_TRACK_MR102];
    let mut pos_indx = [-1i16; NB_PULSE];
    let mut sgn = [0i16; NB_PULSE];
    cod[..L_CODE].fill(0);
    for k in 0..NB_PULSE {
        let i = codvec[k];
        let pos_index = (i >> 2) as i16;
        let track = i & 3;
        let sign_index = if sign[i] > 0 {
            cod[i] = (cod[i] as i32 + 8191) as i16;
            sgn[k] = 32767;
            0
        } else {
            cod[i] = (cod[i] as i32 - 8191) as i16;
            sgn[k] = -32768;
            1
        };
        if pos_indx[track] < 0 {
            pos_indx[track] = pos_index;
            sign_indx[track] = sign_index;
        } else if ((sign_index ^ sign_indx[track]) & 1) == 0 {
            if pos_indx[track] <= pos_index {
                pos_indx[track + NB_TRACK_MR102] = pos_index;
            } else {
                pos_indx[track + NB_TRACK_MR102] = pos_indx[track];
                pos_indx[track] = pos_index;
                sign_indx[track] = sign_index;
            }
        } else if pos_indx[track] <= pos_index {
            pos_indx[track + NB_TRACK_MR102] = pos_indx[track];
            pos_indx[track] = pos_index;
            sign_indx[track] = sign_index;
        } else {
            pos_indx[track + NB_TRACK_MR102] = pos_index;
        }
    }
    filtered_code(h, &codvec, &sgn, y, ov);
    // compress_code
    indx[..NB_TRACK_MR102].copy_from_slice(&sign_indx);
    indx[NB_TRACK_MR102] = compress10(pos_indx[0], pos_indx[4], pos_indx[1]);
    indx[NB_TRACK_MR102 + 1] = compress10(pos_indx[2], pos_indx[6], pos_indx[5]);
    let mut ia = pos_indx[3] >> 1;
    if (pos_indx[7] >> 1) & 1 == 1 {
        ia = 4 - ia;
    }
    let mut ib = (pos_indx[7] >> 1) * 5;
    ib += ia;
    ib <<= 5;
    ib += 12;
    let ic = (((ib as i32 * 1311) >> 15) as i16) << 2;
    indx[NB_TRACK_MR102 + 2] = (pos_indx[3] & 1) + ((pos_indx[7] & 1) << 1) + ic;
}

/// `code_10i40_35bits` (MR122): 10 pulses; writes 10 indices.
pub(super) fn code_10i40_35bits(
    x: &[i16],
    cn: &[i16],
    h: &[i16],
    cod: &mut [i16],
    y: &mut [i16],
    indx: &mut [i16],
    ov: &mut Flag,
) {
    const NB_PULSE: usize = 10;
    let mut ipos = [0i16; NB_PULSE];
    let mut pos_max = [0i16; NB_TRACK];
    let mut codvec = [0usize; NB_PULSE];
    let mut dn = [0i16; L_CODE];
    let mut sign = [0i16; L_CODE];
    let mut rr = RrTm::ZERO;
    cor_h_x(h, x, &mut dn, 2, ov);
    set_sign12k2(
        &mut dn,
        cn,
        &mut sign,
        &mut pos_max,
        NB_TRACK,
        &mut ipos,
        STEP,
        ov,
    );
    cor_h_tm(h, &sign, &mut rr, ov);
    search_10i40(&dn, &rr, &mut ipos, &pos_max, &mut codvec);
    // build_code
    cod[..L_CODE].fill(0);
    indx[..NB_TRACK].fill(-1);
    let mut sgn = [0i16; NB_PULSE];
    for k in 0..NB_PULSE {
        let i = codvec[k];
        let mut index = ((i as i32 * 6554) >> 15) as i16;
        let track = i - (index + (index << 2)) as usize;
        if sign[i] > 0 {
            cod[i] = cod[i].wrapping_add(4096);
            sgn[k] = 8192;
        } else {
            cod[i] = cod[i].wrapping_sub(4096);
            sgn[k] = -8192;
            index += 8;
        }
        let temp = indx[track];
        if temp < 0 {
            indx[track] = index;
        } else if ((index ^ temp) & 8) == 0 {
            if temp <= index {
                indx[track + 5] = index;
            } else {
                indx[track + 5] = temp;
                indx[track] = index;
            }
        } else if (temp & 7) <= (index & 7) {
            indx[track + 5] = temp;
            indx[track] = index;
        } else {
            indx[track + 5] = index;
        }
    }
    // y: the C's L_mult/L_mac chain over the ten +-8192-signed copies of h, then round
    // (saturating; opencore's shortcut, sum(+-64 h) >> 8, wraps for ten aligned pulses).
    filtered_code(h, &codvec, &sgn, y, ov);
    // q_p
    for (n, v) in indx[..10].iter_mut().enumerate() {
        let tmp = *v;
        *v = if n < 5 {
            (tmp & 0x8) | GRAY[(tmp & 0x7) as usize]
        } else {
            GRAY[(tmp & 0x7) as usize]
        };
    }
}

/// `cbsearch`: fixed-codebook search for one subframe. `h` is modified (pitch sharpening).
/// Appends the codebook indices to `ana` at `*ana_pos`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(crate) fn cbsearch(
    x: &[i16],
    h: &mut [i16],
    t0: i16,
    pitch_sharp: i16,
    gain_pit: i16,
    res2: &[i16],
    code: &mut [i16; L_SUBFR],
    y: &mut [i16; L_SUBFR],
    ana: &mut [i16],
    ana_pos: &mut usize,
    mode: Mode,
    sub_nr: usize,
    ov: &mut Flag,
) {
    let mut index = 0i16;
    let p = *ana_pos;
    match mode {
        Mode::MR475 | Mode::MR515 => {
            ana[p] = code_2i40_9bits(sub_nr, x, h, t0, pitch_sharp, code, y, &mut index, ov);
            ana[p + 1] = index;
            *ana_pos += 2;
        }
        Mode::MR59 => {
            ana[p] = code_2i40_11bits(x, h, t0, pitch_sharp, code, y, &mut index, ov);
            ana[p + 1] = index;
            *ana_pos += 2;
        }
        Mode::MR67 => {
            ana[p] = code_3i40_14bits(x, h, t0, pitch_sharp, code, y, &mut index, ov);
            ana[p + 1] = index;
            *ana_pos += 2;
        }
        Mode::MR74 | Mode::MR795 => {
            ana[p] = code_4i40_17bits(x, h, t0, pitch_sharp, code, y, &mut index, ov);
            ana[p + 1] = index;
            *ana_pos += 2;
        }
        Mode::MR102 => {
            let sharp = shl(pitch_sharp, 1, ov);
            sharpen_from(h, t0, sharp, ov);
            code_8i40_31bits(x, res2, h, code, y, &mut ana[p..], ov);
            *ana_pos += 7;
            sharpen_from(code, t0, sharp, ov);
        }
        _ => {
            let sharp = shl(gain_pit, 1, ov);
            for i in t0.max(0) as usize..L_SUBFR {
                let temp = ((h[i - t0 as usize] as i32 * sharp as i32) >> 15) as i16;
                h[i] = add(h[i], temp, ov);
            }
            code_10i40_35bits(x, res2, h, code, y, &mut ana[p..], ov);
            *ana_pos += 10;
            sharpen_from(code, t0, sharp, ov);
        }
    }
}

/// The unconditional `for (i = T0; i < L_SUBFR; i++)` sharpening loops of `cbsearch`.
pub(super) fn sharpen_from(v: &mut [i16], t0: i16, sharp: i16, ov: &mut Flag) {
    for i in t0.max(0) as usize..L_SUBFR {
        let temp = mult(v[i - t0 as usize], sharp, ov);
        v[i] = add(v[i], temp, ov);
    }
}
