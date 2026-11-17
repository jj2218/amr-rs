//! LP analysis and ISF quantisation (vo-amrwbenc `autocorr.c`, `lag_wind.c`, `levinson.c`,
//! `az_isp.c`, `isp_isf.c`, `qpisf_2s.c`, `qisf_ns.c`), with the reference's saturation.
//! ISF → ISP, A(z) and interpolation are the decoder's (`crate::lpc`).

use super::ops::*;
use crate::cnst::M;
use crate::enc_tables::*;
pub(crate) use crate::simd::{autocorr, az_isp, qpisf_2s_36b, qpisf_2s_46b};

/// LPC analysis window (`L_WINDOW`).
pub(crate) const L_WINDOW: usize = 384;
const NC: usize = M / 2;
const GRID_POINTS: usize = 100;

/// `Autocorr` (order 16): windowed autocorrelations `r` in double precision (r_h, r_l).
#[inline(always)]
pub(crate) fn autocorr_portable(x: &[i16], r_h: &mut [i16; M + 1], r_l: &mut [i16; M + 1]) {
    let mut y = [0i16; L_WINDOW];
    for i in 0..L_WINDOW {
        y[i] = mult_r(x[i], VO_WINDOW[i]);
    }
    // both energy chains have terms >= 0: each is min(its plain sum, MAX_32)
    let sat = |s: i64| s.min(i32::MAX as i64) as i32;
    let sum = sat(y
        .iter()
        .map(|&v| l_shr(l_mult(v, v), 8) as i64)
        .sum::<i64>()
        + (16 << 16));
    let shift = sub(4, norm_l(sum) >> 1).max(0);
    for v in &mut y {
        *v = shr_r(*v, shift);
    }
    let see: i64 = y.iter().map(|&v| v as i64 * v as i64).sum();
    let sum = sat(1 + 2 * see);
    let norm = norm_l(sum);
    (r_h[0], r_l[0]) = l_extract(l_shl(sum, norm));
    // every partial lag sum is at most Σy² in magnitude (Cauchy–Schwarz): if 2Σy² < 2^31 the
    // chains are plain sums
    let plain = 2 * see < 1 << 31;
    for i in 1..=M {
        let s = if plain {
            let mut s = 0i32;
            for j in 0..L_WINDOW - i {
                s = s.wrapping_add(y[j] as i32 * y[j + i] as i32);
            }
            s << 1
        } else {
            let mut s = 0i32;
            for j in 0..L_WINDOW - i {
                s = l_mac(s, y[j], y[j + i]);
            }
            s
        };
        (r_h[i], r_l[i]) = l_extract(l_shl(s, norm));
    }
}

/// `Lag_window`.
pub(crate) fn lag_window(r_h: &mut [i16; M + 1], r_l: &mut [i16; M + 1]) {
    for i in 1..=M {
        (r_h[i], r_l[i]) = l_extract(mpy_32(r_h[i], r_l[i], VOLAG_H[i - 1], VOLAG_L[i - 1]));
    }
}

/// `Levinson`: A(z) (Q12) and reflection coefficients from the autocorrelations; `mem` keeps
/// the last stable A(z) (16) and two reflection coefficients.
pub(crate) fn levinson(
    rh: &[i16; M + 1],
    rl: &[i16; M + 1],
    a: &mut [i16],
    rc: &mut [i16; M],
    mem: &mut [i16; M + 2],
) {
    let (mut ah, mut al) = ([0i16; M + 1], [0i16; M + 1]);
    let (mut anh, mut anl) = ([0i16; M + 1], [0i16; M + 1]);
    let t1 = l_comp(rh[1], rl[1]);
    let mut t0 = div_32(l_abs(t1), rh[0], rl[0]);
    if t1 > 0 {
        t0 = l_negate(t0);
    }
    let (mut kh, mut kl) = l_extract(t0);
    rc[0] = kh;
    (ah[1], al[1]) = l_extract(l_shr(t0, 4));
    let t0 = l_sub(i32::MAX, l_abs(mpy_32(kh, kl, kh, kl)));
    let (hi, lo) = l_extract(t0);
    let t0 = mpy_32(rh[0], rl[0], hi, lo);
    let mut alp_exp = norm_l(t0);
    let (mut alp_h, mut alp_l) = l_extract(l_shl(t0, alp_exp));
    for i in 2..=M {
        let mut t0 = 0i32;
        for j in 1..i {
            t0 = l_add(t0, mpy_32(rh[j], rl[j], ah[i - j], al[i - j]));
        }
        t0 = l_add(l_shl(t0, 4), l_comp(rh[i], rl[i]));
        let mut t2 = div_32(l_abs(t0), alp_h, alp_l);
        if t0 > 0 {
            t2 = l_negate(t2);
        }
        t2 = l_shl(t2, alp_exp);
        (kh, kl) = l_extract(t2);
        rc[i - 1] = kh;
        if abs_s(kh) > 32750 {
            a[0] = 4096;
            a[1..=M].copy_from_slice(&mem[..M]);
            rc[0] = mem[M];
            rc[1] = mem[M + 1];
            return;
        }
        for j in 1..i {
            let t0 = l_add(mpy_32(kh, kl, ah[i - j], al[i - j]), l_comp(ah[j], al[j]));
            (anh[j], anl[j]) = l_extract(t0);
        }
        (anh[i], anl[i]) = l_extract(l_shr(t2, 4));
        let t0 = l_sub(i32::MAX, l_abs(mpy_32(kh, kl, kh, kl)));
        let (hi, lo) = l_extract(t0);
        let t0 = mpy_32(alp_h, alp_l, hi, lo);
        let j = norm_l(t0);
        (alp_h, alp_l) = l_extract(l_shl(t0, j));
        alp_exp = add(alp_exp, j);
        ah[1..=i].copy_from_slice(&anh[1..=i]);
        al[1..=i].copy_from_slice(&anl[1..=i]);
    }
    a[0] = 4096;
    for i in 1..=M {
        a[i] = round(l_shl(l_comp(ah[i], al[i]), 1));
        mem[i - 1] = a[i];
    }
    mem[M] = rc[0];
    mem[M + 1] = rc[1];
}

/// `Az_isp`: the ISPs of A(z) (roots of F1, F2 by Chebyshev evaluation on a grid); `old_isp`
/// if fewer than 15 roots are found.
#[inline(always)]
pub(crate) fn az_isp_portable(a: &[i16], isp: &mut [i16; M], old_isp: &[i16; M]) {
    let (mut f1, mut f2) = ([0i16; NC + 1], [0i16; NC]);
    for i in 0..NC {
        let t0 = l_mult(a[i], 16384);
        f1[i] = round(l_mac(t0, a[M - i], 16384));
        f2[i] = round(l_msu(t0, a[M - i], 16384));
    }
    f1[NC] = a[NC];
    for i in 2..NC {
        f2[i] = add(f2[i], f2[i - 2]);
    }
    let mut nf = 0;
    let mut ip = 0;
    let cheb = |x: i16, ip: i32| {
        if ip == 0 {
            chebps2(x, &f1, NC)
        } else {
            chebps2(x, &f2, NC - 1)
        }
    };
    // the series at the grid points don't depend on the scan: both tables first, in lanes
    let grid = [cheb_grid(&f1, NC), cheb_grid(&f2, NC - 1)];
    let mut xlow = VOGRID[0];
    let mut ylow = grid[0][0];
    let mut j = 0;
    while nf < M - 1 && j < GRID_POINTS {
        j += 1;
        let mut xhigh = xlow;
        let mut yhigh = ylow;
        xlow = VOGRID[j];
        ylow = grid[ip as usize][j];
        if l_mult(ylow, yhigh) <= 0 {
            for _ in 0..2 {
                let xmid = add(xlow >> 1, xhigh >> 1);
                let ymid = cheb(xmid, ip);
                if l_mult(ylow, ymid) <= 0 {
                    yhigh = ymid;
                    xhigh = xmid;
                } else {
                    ylow = ymid;
                    xlow = xmid;
                }
            }
            let x = sub(xhigh, xlow);
            let y = sub(yhigh, ylow);
            let xint = if y == 0 {
                xlow
            } else {
                let sign = y;
                let y = abs_s(y);
                let exp = norm_s(y);
                let y = div_s(16383, shl(y, exp));
                let t0 = l_shr(l_mult(x, y), sub(20, exp));
                let mut y = extract_l(t0);
                if sign < 0 {
                    y = negate(y);
                }
                let t0 = l_shr(l_mult(ylow, y), 11);
                sub(xlow, extract_l(t0))
            };
            isp[nf] = xint;
            xlow = xint;
            nf += 1;
            ip ^= 1;
            ylow = cheb(xlow, ip);
        }
    }
    if nf < M - 1 {
        isp.copy_from_slice(old_isp);
    } else {
        isp[M - 1] = shl(a[M], 3);
    }
}

/// `Chebps2`: the Chebyshev series of order `n` at `x`, Q14 (saturated, never -32768).
#[inline(always)]
fn chebps2(x: i16, f: &[i16], n: usize) -> i16 {
    chebps2_lanes(&[x], f, n)[0]
}

/// `chebps2` at `L` points: `chebps2_plain`, else (some |u| over its bound) the C's chain.
#[inline(always)]
fn chebps2_lanes<const L: usize>(x: &[i16; L], f: &[i16], n: usize) -> [i16; L] {
    chebps2_plain(x, f, n).unwrap_or_else(|| chebps2_chain(x, f, n))
}

/// Largest |u| for which `chebps2_plain`'s recursion is the C's chain: a level's partial sums
/// are within 6·U + 2^28 + 2^17 (`2·Mpy_32_16` <= 2·U, `2^15·b2_h` <= U + 2^15, `4096·f_i`
/// <= 2^27, doubled by `L_shl`, then `2·b2_l` < 2^17), the last level's within 4·U + 2^17 +
/// 2^27.
const CHEB_U_MAX: i64 = 313_000_000;
const _: () = assert!(6 * CHEB_U_MAX + (1 << 28) + (1 << 17) < 1 << 31);

/// `chebps2` without saturation: each (hi, lo) `b` of the C, split by `L_Extract` from an even
/// `t`, is one word u = t >> 1 = hi·2^15 + lo, and `Mpy_32_16(hi, lo, x)` is 2·((u·x) >> 15).
/// A level `Mpy_32_16(b1, x)`, `L_mac(b2_h, -16384)`, `L_mac(f_i, 2048)`, `L_shl(1)`,
/// `L_msu(b2_l, 1)` is then u0 = 2·((u1·x) >> 15) - u2 + 4096·f_i (the Chebyshev
/// recursion), the last t = 2·((u1·x) >> 15) - 2·u2 + 4096·f_n, and `extract_h(L_shl(t, 6))`
/// = clamp16(t >> 10) (saturation is a monotone clamp). Exact while every |u| <=
/// `CHEB_U_MAX` (checked over all lanes and levels; `None` otherwise).
#[inline(always)]
fn chebps2_plain<const L: usize>(x: &[i16; L], f: &[i16], n: usize) -> Option<[i16; L]> {
    let m = |u: i64, x: i16| 2 * ((u * x as i64) >> 15);
    let u2_0 = 4096 * f[0] as i64;
    let mut u2 = [u2_0; L];
    let mut u1: [i64; L] = std::array::from_fn(|i| m(u2_0, x[i]) + 4096 * f[1] as i64);
    let mut big = [0i64; L];
    for &fi in &f[2..n] {
        for i in 0..L {
            let u0 = m(u1[i], x[i]) - u2[i] + 4096 * fi as i64;
            big[i] = big[i].max(u0.abs());
            (u2[i], u1[i]) = (u1[i], u0);
        }
    }
    if big.iter().any(|&b| b > CHEB_U_MAX) {
        return None;
    }
    Some(std::array::from_fn(|i| {
        let t = m(u1[i], x[i]) - 2 * u2[i] + 2048 * 2 * f[n] as i64;
        let cheb = (t >> 10).clamp(-32768, 32767) as i16;
        if cheb == -32768 { -32767 } else { cheb }
    }))
}

/// `chebps2` at `L` points as the C's chain, the same operations lane by lane (coefficients
/// outer, lanes inner, so the lanes vectorise).
#[inline(always)]
fn chebps2_chain<const L: usize>(x: &[i16; L], f: &[i16], n: usize) -> [i16; L] {
    let (h, l) = l_extract(l_mult(f[0], 4096));
    let (mut b2_h, mut b2_l) = ([h; L], [l; L]);
    let (mut b1_h, mut b1_l) = ([0i16; L], [0i16; L]);
    for i in 0..L {
        let t0 = l_mac(l_shl(mpy_32_16(b2_h[i], b2_l[i], x[i]), 1), f[1], 4096);
        (b1_h[i], b1_l[i]) = l_extract(t0);
    }
    for &fi in &f[2..n] {
        for i in 0..L {
            let mut t0 = mpy_32_16(b1_h[i], b1_l[i], x[i]);
            t0 = l_mac(t0, b2_h[i], -16384);
            t0 = l_mac(t0, fi, 2048);
            t0 = l_shl(t0, 1);
            t0 = l_msu(t0, b2_l[i], 1);
            let (b0_h, b0_l) = l_extract(t0);
            (b2_h[i], b2_l[i]) = (b1_h[i], b1_l[i]);
            (b1_h[i], b1_l[i]) = (b0_h, b0_l);
        }
    }
    std::array::from_fn(|i| {
        let mut t0 = mpy_32_16(b1_h[i], b1_l[i], x[i]);
        t0 = l_mac(t0, b2_h[i], -32768);
        t0 = l_msu(t0, b2_l[i], 1);
        t0 = l_mac(t0, f[n], 2048);
        let cheb = extract_h(l_shl(t0, 6));
        if cheb == -32768 { -32767 } else { cheb }
    })
}

/// Both polynomials' series at every grid point (`az_isp`'s scan reads them in order).
#[inline(always)]
fn cheb_grid(f: &[i16], n: usize) -> [i16; GRID_POINTS + 1] {
    const P: usize = (GRID_POINTS + 1).next_multiple_of(16);
    let x: [i16; P] = std::array::from_fn(|j| VOGRID[j.min(GRID_POINTS)]);
    let mut y = [0i16; GRID_POINTS + 1];
    for (b, xb) in x.as_chunks::<16>().0.iter().enumerate() {
        let yb = chebps2_lanes(xb, f, n);
        let k = (16 * b + 16).min(GRID_POINTS + 1) - 16 * b;
        y[16 * b..16 * b + k].copy_from_slice(&yb[..k]);
    }
    y
}

/// `Isp_isf` (order `m` = `isp.len()`): ISPs to normalised frequencies.
pub(crate) fn isp_isf(isp: &[i16], isf: &mut [i16]) {
    let m = isp.len();
    let mut ind = 127usize;
    for i in (0..m).rev() {
        if i >= m - 2 {
            ind = 127;
        }
        while TABLE[ind] < isp[i] {
            ind -= 1;
        }
        let l = l_mult(sub(isp[i], TABLE[ind]), SLOPE[ind]);
        isf[i] = round(l_shl(l, 4));
        isf[i] = add(isf[i], shl(ind as i16, 7));
    }
    isf[m - 1] = shr(isf[m - 1], 1);
}

/// 1/3 in Q15 (`MU`): the ISF prediction factor.
const MU: i16 = 10923;
const N_SURV_MAX: usize = 4;

/// A VQ codebook of `N` `D`-vectors (`flat`, the generated table) and its component-major
/// copy (`cols[j][i]`: entry i's component j), so the distances run in lanes over the entries;
/// `max` the largest |component|.
pub(crate) struct Book<const D: usize, const N: usize> {
    flat: &'static [i16],
    cols: [[i16; N]; D],
    max: i32,
}

impl<const D: usize, const N: usize> Book<D, N> {
    const fn new(flat: &'static [i16]) -> Self {
        assert!(flat.len() == D * N);
        let mut cols = [[0i16; N]; D];
        let mut max = 0;
        let mut i = 0;
        while i < N {
            let mut j = 0;
            while j < D {
                let v = flat[i * D + j];
                cols[j][i] = v;
                let a = (v as i32).abs();
                if a > max {
                    max = a;
                }
                j += 1;
            }
            i += 1;
        }
        Self { flat, cols, max }
    }

    /// Every entry's distance to `x`, the reference's chain `Σ L_mult(t, t)` with
    /// `t = sub(x, c)`, if no step can saturate: with `b = max|x| + max|c|`, `D·b² < 2^30`
    /// keeps every `t` exact and the sum of `2t²` below 2^31. `None`: the chain.
    #[inline(always)]
    fn dists(&self, x: &[i16]) -> Option<[i32; N]> {
        let x: &[i16; D] = x[..D].try_into().unwrap();
        let b = x.iter().map(|&v| (v as i32).abs()).max().unwrap_or(0) as i64 + self.max as i64;
        if D as i64 * b * b >= 1 << 30 {
            return None;
        }
        let mut d = [0i32; N];
        for (j, col) in self.cols.iter().enumerate() {
            let xj = x[j] as i32;
            for (v, &c) in d.iter_mut().zip(col) {
                let t = xj - c as i32;
                *v += t * t;
            }
        }
        Some(d.map(|v| v << 1))
    }

    /// The reference's distance chain for entry `i`.
    fn chain(&self, x: &[i16], i: usize) -> i32 {
        let mut dist = 0i32;
        for j in 0..D {
            let t = sub(x[j], self.flat[i * D + j]);
            dist = l_mac(dist, t, t);
        }
        dist
    }

    /// `Sub_VQ`: the nearest entry to `x` (written to `x`); (index, distance). The first
    /// minimum of the distances (the C's strict `<` from `MAX_32`).
    #[inline(always)]
    pub(crate) fn sub_vq(&self, x: &mut [i16]) -> (i16, i32) {
        let (index, dist_min) = match self.dists(x) {
            Some(d) => {
                let m = d.iter().copied().min().unwrap_or(i32::MAX);
                (d.iter().position(|&v| v == m).unwrap_or(0), m)
            }
            None => {
                let mut dist_min = i32::MAX;
                let mut index = 0;
                for i in 0..N {
                    let dist = self.chain(x, i);
                    if dist < dist_min {
                        dist_min = dist;
                        index = i;
                    }
                }
                (index, dist_min)
            }
        };
        x[..D].copy_from_slice(&self.flat[index * D..index * D + D]);
        (index as i16, dist_min)
    }

    /// `VQ_stage1`: the `surv` nearest entries, best first (ties: the lower index).
    #[inline(always)]
    fn vq_stage1(&self, x: &[i16], index: &mut [i16], surv: usize) {
        let mut dist_min = [i32::MAX; N_SURV_MAX];
        for (i, v) in index[..surv].iter_mut().enumerate() {
            *v = i as i16;
        }
        let fast = self.dists(x);
        for i in 0..N {
            let dist = match &fast {
                Some(d) => d[i],
                None => self.chain(x, i),
            };
            // dist_min is ascending: not below its last, no insertion
            if dist >= dist_min[surv - 1] {
                continue;
            }
            for k in 0..surv {
                if dist < dist_min[k] {
                    for l in (k + 1..surv).rev() {
                        dist_min[l] = dist_min[l - 1];
                        index[l] = index[l - 1];
                    }
                    dist_min[k] = dist;
                    index[k] = i as i16;
                    break;
                }
            }
        }
    }
}

static B1: Book<9, 256> = Book::new(&DICO1_ISF);
static B2: Book<7, 256> = Book::new(&DICO2_ISF);
static B21: Book<3, 64> = Book::new(&DICO21_ISF);
static B22: Book<3, 128> = Book::new(&DICO22_ISF);
static B23: Book<3, 128> = Book::new(&DICO23_ISF);
static B24: Book<3, 32> = Book::new(&DICO24_ISF);
static B25: Book<4, 32> = Book::new(&DICO25_ISF);
static B21_36: Book<5, 128> = Book::new(&DICO21_ISF_36B);
static B22_36: Book<4, 128> = Book::new(&DICO22_ISF_36B);
static B23_36: Book<7, 64> = Book::new(&DICO23_ISF_36B);
static BN1: Book<2, 64> = Book::new(&DICO1_ISF_NOISE);
static BN2: Book<3, 64> = Book::new(&DICO2_ISF_NOISE);
static BN3: Book<3, 64> = Book::new(&DICO3_ISF_NOISE);
static BN4: Book<4, 32> = Book::new(&DICO4_ISF_NOISE);
static BN5: Book<4, 32> = Book::new(&DICO5_ISF_NOISE);

/// The prediction residual `isf1 - mean - MU·past` of both quantisers.
fn residual(isf1: &[i16; M], past_isfq: &[i16; M]) -> [i16; M] {
    std::array::from_fn(|i| sub(sub(isf1[i], MEAN_ISF[i]), mult(MU, past_isfq[i])))
}

/// One split of a two-stage VQ: for each first-stage survivor of `isf` (dimension `dim1`,
/// codebook `dico1`), the second-stage splits `subs` (codebook, dimension); the best total
/// error's indices go to `ind1` and `inds` (unchanged if every total is `MAX_32`).
#[inline(always)]
fn two_stage<const D: usize>(
    isf: &[i16],
    book1: &Book<D, 256>,
    stage2: impl Fn(&mut [i16; 9], &mut [i16; 3]) -> i32,
    nb_surv: usize,
    ind1: &mut i16,
    inds: &mut [i16],
) {
    let mut surv1 = [0i16; N_SURV_MAX];
    book1.vq_stage1(isf, &mut surv1, nb_surv);
    let mut distance = i32::MAX;
    for &s in &surv1[..nb_surv] {
        let mut res2: [i16; 9] = [0; 9];
        for i in 0..D {
            res2[i] = sub(isf[i], book1.flat[i + s as usize * D]);
        }
        let mut tmp_ind = [0i16; 3];
        let temp = stage2(&mut res2, &mut tmp_ind);
        if temp < distance {
            distance = temp;
            *ind1 = s;
            let n = inds.len();
            inds.copy_from_slice(&tmp_ind[..n]);
        }
    }
}

/// `Qpisf_2s_46b`: 46-bit predictive two-stage ISF quantisation; `isf_q` the quantised ISFs.
#[inline(always)]
pub(crate) fn qpisf_2s_46b_portable(
    isf1: &[i16; M],
    isf_q: &mut [i16; M],
    past_isfq: &mut [i16; M],
    indice: &mut [i16; 7],
    nb_surv: usize,
) {
    let isf = residual(isf1, past_isfq);
    let (i0, rest) = indice.split_at_mut(1);
    let (i1, rest) = rest.split_at_mut(1);
    two_stage(
        &isf,
        &B1,
        |x, ind| {
            let (a, e0) = B21.sub_vq(&mut x[0..]);
            let (b, e1) = B22.sub_vq(&mut x[3..]);
            let (c, e2) = B23.sub_vq(&mut x[6..]);
            *ind = [a, b, c];
            l_add(l_add(e0, e1), e2)
        },
        nb_surv,
        &mut i0[0],
        &mut rest[..3],
    );
    two_stage(
        &isf[9..],
        &B2,
        |x, ind| {
            let (a, e0) = B24.sub_vq(&mut x[0..]);
            let (b, e1) = B25.sub_vq(&mut x[3..]);
            ind[..2].copy_from_slice(&[a, b]);
            l_add(e0, e1)
        },
        nb_surv,
        &mut i1[0],
        &mut rest[3..5],
    );
    let old = *isf_q;
    crate::lpc::dpisf_2s_46b(&indice[..], isf_q, past_isfq, &old, &mut [0; 3 * M], false);
}

/// `Qpisf_2s_36b`: 36-bit predictive two-stage ISF quantisation (6.60 kbit/s).
#[inline(always)]
pub(crate) fn qpisf_2s_36b_portable(
    isf1: &[i16; M],
    isf_q: &mut [i16; M],
    past_isfq: &mut [i16; M],
    indice: &mut [i16; 5],
    nb_surv: usize,
) {
    let isf = residual(isf1, past_isfq);
    let (i0, rest) = indice.split_at_mut(1);
    let (i1, rest) = rest.split_at_mut(1);
    two_stage(
        &isf,
        &B1,
        |x, ind| {
            let (a, e0) = B21_36.sub_vq(&mut x[0..]);
            let (b, e1) = B22_36.sub_vq(&mut x[5..]);
            ind[..2].copy_from_slice(&[a, b]);
            l_add(e0, e1)
        },
        nb_surv,
        &mut i0[0],
        &mut rest[..2],
    );
    two_stage(
        &isf[9..],
        &B2,
        |x, ind| {
            let (a, e0) = B23_36.sub_vq(&mut x[0..]);
            ind[0] = a;
            e0
        },
        nb_surv,
        &mut i1[0],
        &mut rest[2..3],
    );
    let old = *isf_q;
    crate::lpc::dpisf_2s_36b(&indice[..], isf_q, past_isfq, &old, &mut [0; 3 * M], false);
}

/// `Qisf_ns`: comfort-noise ISF quantisation (SID frames).
pub(crate) fn qisf_ns(isf1: &[i16; M], isf_q: &mut [i16; M], indice: &mut [i16; 5]) {
    for i in 0..M {
        isf_q[i] = sub(isf1[i], MEAN_ISF_NOISE[i]);
    }
    indice[0] = BN1.sub_vq(&mut isf_q[0..]).0;
    indice[1] = BN2.sub_vq(&mut isf_q[2..]).0;
    indice[2] = BN3.sub_vq(&mut isf_q[5..]).0;
    indice[3] = BN4.sub_vq(&mut isf_q[8..]).0;
    indice[4] = BN5.sub_vq(&mut isf_q[12..]).0;
    crate::lpc::disf_ns(indice, isf_q);
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::encoder::filters::vs_ref::signal;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_Autocorr(x: *mut i16, m: i16, r_h: *mut i16, r_l: *mut i16);
        fn refwb_Lag_window(r_h: *mut i16, r_l: *mut i16);
        fn refwb_Levinson(rh: *mut i16, rl: *mut i16, a: *mut i16, rc: *mut i16, mem: *mut i16);
        fn refwb_Az_isp(a: *mut i16, isp: *mut i16, old_isp: *mut i16);
        fn refwb_Isp_isf(isp: *mut i16, isf: *mut i16, m: i16);
        fn refwb_Qpisf_2s_46b(
            isf1: *mut i16,
            isf_q: *mut i16,
            past: *mut i16,
            ind: *mut i16,
            nb: i16,
        );
        fn refwb_Qpisf_2s_36b(
            isf1: *mut i16,
            isf_q: *mut i16,
            past: *mut i16,
            ind: *mut i16,
            nb: i16,
        );
        fn refwb_Qisf_ns(isf1: *mut i16, isf_q: *mut i16, ind: *mut i16);
    }

    /// The ISF quantisers on arbitrary vectors (any i16, a past prediction too): the
    /// distance chains saturate (the lanes' bound fails) and every tie rule shows.
    #[test]
    fn isf_quantisers_on_arbitrary_isfs() {
        let mut r = Rng::new(32);
        for n in 0..20_000u64 {
            let scale = (n % 16) as u32;
            let isf: [i16; M] = std::array::from_fn(|_| r.i16() >> scale);
            let past: [i16; M] = std::array::from_fn(|_| r.i16() >> (n % 7));
            let nb = 1 + (n % 4) as usize;
            let (mut p46, mut p36) = (past, past);
            let (mut p46c, mut p36c) = (past, past);
            let (mut q, mut qc) = ([0i16; M], [0i16; M]);
            let (mut i7, mut i7c, mut i5, mut i5c) = ([0i16; 7], [0i16; 7], [0i16; 5], [0i16; 5]);
            let mut isfc = isf;
            qpisf_2s_46b(&isf, &mut q, &mut p46, &mut i7, nb);
            with_ref(|| {
                // SAFETY: buffers of the sizes the reference reads and writes.
                unsafe {
                    refwb_Qpisf_2s_46b(
                        isfc.as_mut_ptr(),
                        qc.as_mut_ptr(),
                        p46c.as_mut_ptr(),
                        i7c.as_mut_ptr(),
                        nb as i16,
                    );
                    assert_eq!((q, p46, i7), (qc, p46c, i7c), "Qpisf_2s_46b {n}");
                    qpisf_2s_36b(&isf, &mut q, &mut p36, &mut i5, nb);
                    refwb_Qpisf_2s_36b(
                        isfc.as_mut_ptr(),
                        qc.as_mut_ptr(),
                        p36c.as_mut_ptr(),
                        i5c.as_mut_ptr(),
                        nb as i16,
                    );
                    assert_eq!((q, p36, i5), (qc, p36c, i5c), "Qpisf_2s_36b {n}");
                    qisf_ns(&isf, &mut q, &mut i5);
                    refwb_Qisf_ns(isfc.as_mut_ptr(), qc.as_mut_ptr(), i5c.as_mut_ptr());
                    assert_eq!((q, i5), (qc, i5c), "Qisf_ns {n}");
                }
            });
        }
    }

    /// The analysis chain on signal classes (white, square, quiet, loud), sequences of 300
    /// frames (Levinson's fallback memory): Autocorr → Lag_window → Levinson → Az_isp →
    /// Isp_isf → both quantisers (their predictor memories) and the noise quantiser.
    #[test]
    fn lp_analysis_matches_reference() {
        let mut r = Rng::new(31);
        for seq in 0..40u64 {
            let (mut mem, mut memc) = ([0i16; M + 2], [0i16; M + 2]);
            let mut old_isp = ISP_INIT;
            let (mut p46, mut p46c, mut p36, mut p36c) =
                ([0i16; M], [0i16; M], [0i16; M], [0i16; M]);
            for call in 0..300 {
                let mut x = signal(&mut r, seq * 7 + call / 50, L_WINDOW);
                if call % 13 == 0 {
                    // a pure tone: near-singular autocorrelation (Levinson's unstable branch)
                    let f = (r.u64() % 100) as f64 / 400.0;
                    for (n, v) in x.iter_mut().enumerate() {
                        *v = (20000.0 * (std::f64::consts::TAU * f * n as f64).sin()) as i16;
                    }
                }
                if call % 11 == 3 {
                    // a full-scale square: the scaled window's energy chain still saturates
                    let p = 1 + (r.u64() % 40) as usize;
                    for (n, v) in x.iter_mut().enumerate() {
                        *v = if (n / p).is_multiple_of(2) {
                            i16::MAX
                        } else {
                            i16::MIN
                        };
                    }
                }
                let (mut rh, mut rl) = ([0i16; M + 1], [0i16; M + 1]);
                let (mut rhc, mut rlc) = ([0i16; M + 1], [0i16; M + 1]);
                let (mut a, mut ac) = ([0i16; M + 1], [0i16; M + 1]);
                let (mut rc, mut rcc) = ([0i16; M], [0i16; M]);
                let (mut isp, mut ispc) = ([0i16; M], [0i16; M]);
                let (mut isf, mut isfc) = ([0i16; M], [0i16; M]);
                autocorr(&x, &mut rh, &mut rl);
                lag_window(&mut rh, &mut rl);
                levinson(&rh, &rl, &mut a, &mut rc, &mut mem);
                az_isp(&a, &mut isp, &old_isp);
                isp_isf(&isp, &mut isf);
                let nb = 1 + (r.u64() % 4) as usize;
                let (mut q, mut qc) = ([0i16; M], [0i16; M]);
                let (mut i7, mut i7c) = ([0i16; 7], [0i16; 7]);
                let (mut i5, mut i5c) = ([0i16; 5], [0i16; 5]);
                with_ref(|| {
                    // SAFETY: buffers of the sizes the reference reads and writes.
                    unsafe {
                        refwb_Autocorr(
                            x.as_mut_ptr(),
                            M as i16,
                            rhc.as_mut_ptr(),
                            rlc.as_mut_ptr(),
                        );
                        refwb_Lag_window(rhc.as_mut_ptr(), rlc.as_mut_ptr());
                        assert_eq!((rh, rl), (rhc, rlc), "Autocorr+Lag_window {seq} {call}");
                        refwb_Levinson(
                            rhc.as_mut_ptr(),
                            rlc.as_mut_ptr(),
                            ac.as_mut_ptr(),
                            rcc.as_mut_ptr(),
                            memc.as_mut_ptr(),
                        );
                        assert_eq!((a, rc, mem), (ac, rcc, memc), "Levinson {seq} {call}");
                        refwb_Az_isp(ac.as_mut_ptr(), ispc.as_mut_ptr(), old_isp.as_mut_ptr());
                        assert_eq!(isp, ispc, "Az_isp {seq} {call}");
                        refwb_Isp_isf(ispc.as_mut_ptr(), isfc.as_mut_ptr(), M as i16);
                        assert_eq!(isf, isfc, "Isp_isf {seq} {call}");
                        qpisf_2s_46b(&isf, &mut q, &mut p46, &mut i7, nb);
                        refwb_Qpisf_2s_46b(
                            isfc.as_mut_ptr(),
                            qc.as_mut_ptr(),
                            p46c.as_mut_ptr(),
                            i7c.as_mut_ptr(),
                            nb as i16,
                        );
                        assert_eq!((q, p46, i7), (qc, p46c, i7c), "Qpisf_2s_46b {seq} {call}");
                        qpisf_2s_36b(&isf, &mut q, &mut p36, &mut i5, nb);
                        refwb_Qpisf_2s_36b(
                            isfc.as_mut_ptr(),
                            qc.as_mut_ptr(),
                            p36c.as_mut_ptr(),
                            i5c.as_mut_ptr(),
                            nb as i16,
                        );
                        assert_eq!((q, p36, i5), (qc, p36c, i5c), "Qpisf_2s_36b {seq} {call}");
                        qisf_ns(&isf, &mut q, &mut i5);
                        refwb_Qisf_ns(isfc.as_mut_ptr(), qc.as_mut_ptr(), i5c.as_mut_ptr());
                        assert_eq!((q, i5), (qc, i5c), "Qisf_ns {seq} {call}");
                    }
                });
                old_isp = isp;
            }
        }
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref_direct {
    use super::*;
    use crate::encoder::filters::vs_ref::lpc;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_Levinson(rh: *mut i16, rl: *mut i16, a: *mut i16, rc: *mut i16, mem: *mut i16);
        fn refwb_Az_isp(a: *mut i16, isp: *mut i16, old_isp: *mut i16);
    }

    /// Arbitrary autocorrelations (often not positive definite: Levinson's unstable-filter
    /// fallback) and arbitrary A(z) (fewer than 15 ISPs found: Az_isp's fallback).
    #[test]
    fn levinson_and_az_isp_fallbacks_match_reference() {
        let mut r = Rng::new(32);
        let (mut mem, mut memc) = ([0i16; M + 2], [0i16; M + 2]);
        for n in 0..20_000 {
            let mut rh = [0i16; M + 1];
            let mut rl = [0i16; M + 1];
            rh[0] = 0x4000 | (r.u64() as i16 & 0x3fff);
            for i in 0..=M {
                if i > 0 {
                    rh[i] = ((r.i16() as i32 * rh[0] as i32) >> (15 + r.u64() % 3)) as i16;
                }
                rl[i] = (r.u64() & 0x7fff) as i16;
            }
            let (mut a, mut ac) = ([0i16; M + 1], [0i16; M + 1]);
            let (mut rc, mut rcc) = ([0i16; M], [0i16; M]);
            levinson(&rh, &rl, &mut a, &mut rc, &mut mem);
            let mut q: [i16; M + 1] = lpc(&mut r, M).try_into().unwrap();
            if n % 4 == 0 {
                // full-scale A(z) (some all one sign): the Chebyshev recursion past the plain
                // form's bound (`chebps2_plain`), the series saturating
                let kind = n / 4 % 3;
                for v in &mut q[1..] {
                    let big = 32767 - (r.u64() % 64) as i16;
                    *v = match kind {
                        0 => big,
                        1 => -big,
                        _ => r.i16(),
                    };
                }
            }
            let mut old: [i16; M] = std::array::from_fn(|_| r.i16());
            let (mut isp, mut ispc) = ([0i16; M], [0i16; M]);
            az_isp(&q, &mut isp, &old);
            with_ref(|| {
                // SAFETY: buffers of the sizes the reference reads and writes.
                unsafe {
                    refwb_Levinson(
                        rh.as_mut_ptr(),
                        rl.as_mut_ptr(),
                        ac.as_mut_ptr(),
                        rcc.as_mut_ptr(),
                        memc.as_mut_ptr(),
                    );
                    refwb_Az_isp(q.as_mut_ptr(), ispc.as_mut_ptr(), old.as_mut_ptr());
                }
            });
            assert_eq!((a, rc, mem), (ac, rcc, memc), "Levinson");
            assert_eq!(isp, ispc, "Az_isp");
        }
    }
}
