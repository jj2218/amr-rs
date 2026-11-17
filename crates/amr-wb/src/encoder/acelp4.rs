//! The 4-track algebraic codebook search of 8.85..23.85 kbit/s (vo-amrwbenc `c4t64fx.c`, with
//! CR 26.173-0036's impulse-response scaling, which vo predates) and the pulse index coding
//! (`q_pulse.c`), with the reference's saturation.

use super::ops::*;
use crate::basic_op::dot_product12;
use crate::cnst::L_SUBFR;
use crate::enc_tables::TIPOS;
pub(crate) use crate::simd::acelp_4t64_fx;

const NB_PULSE_MAX: usize = 24;
const NB_TRACK: usize = 4;
const STEP: usize = 4;
const NB_POS: usize = 16;
const MSIZE: usize = 256;
const NB_MAX: i16 = 8;
const NPMAXPT: usize = NB_PULSE_MAX.div_ceil(NB_TRACK);

/// `ACELP_4t64_fx`: `nbbits` (20, 36, 44, 52, 64, 72 or 88) pulses' worth of search on four
/// interleaved tracks; `code` (Q9), `y` (filtered, Q9) and the track indices `index`.
#[inline(always)]
pub(crate) fn acelp_4t64_fx_portable(
    dn: &mut [i16; L_SUBFR],
    cn: &[i16],
    hh: &[i16],
    code: &mut [i16; L_SUBFR],
    y: &mut [i16; L_SUBFR],
    nbbits: i16,
    ser_size: i16,
    index_out: &mut [i16; 8],
) {
    let mut nbpos = [0i16; 10];
    let (nbiter, alp0, nb_pulse): (usize, i16, usize) = match nbbits {
        20 => {
            nbpos[..2].copy_from_slice(&[4, 8]);
            (4, 8192, 4)
        }
        36 => {
            nbpos[..3].copy_from_slice(&[4, 8, 8]);
            (4, 4096, 8)
        }
        44 => {
            nbpos[..4].copy_from_slice(&[4, 6, 8, 8]);
            (4, 4096, 10)
        }
        52 => {
            nbpos[..4].copy_from_slice(&[4, 6, 8, 8]);
            (4, 4096, 12)
        }
        64 => {
            nbpos[..6].copy_from_slice(&[4, 4, 6, 6, 8, 8]);
            (3, 3277, 16)
        }
        72 => {
            nbpos[..7].copy_from_slice(&[2, 3, 4, 5, 6, 7, 8]);
            (3, 3072, 18)
        }
        88 => {
            nbpos.copy_from_slice(&[2, 2, 3, 4, 5, 6, 7, 8, 8, 8]);
            (if ser_size > 462 { 1 } else { 2 }, 2048, 24)
        }
        _ => (0, 0, 0),
    };
    let mut codvec = [0i16; NB_PULSE_MAX];
    for (i, c) in codvec[..nb_pulse].iter_mut().enumerate() {
        *c = i as i16;
    }
    let mut exp = 0;
    let mut s = dot_product12(cn, cn, &mut exp);
    isqrt_n(&mut s, &mut exp);
    let k_cn = round(l_shl(s, add(exp, 5)));
    let mut s = dot_product12(&dn[..], &dn[..], &mut exp);
    isqrt_n(&mut s, &mut exp);
    let k_dn = mult_r(alp0, round(l_shl(s, add(exp, 5 + 3))));
    let mut dn2 = [0i16; L_SUBFR];
    // k_cn = round(L_shl(s, ·)) of a positive isqrt_n, k_dn = mult_r(alp0 > 0, that form):
    // both >= 0, so no `L_mult` saturates, t = k_cn·cn + k_dn·dn is exact in i32, and the chain
    // extract_h(L_shl(L_mac(L_mult(k_cn, cn), k_dn, dn), 8)) = extract_h(sat(sat(2t)·256)) is
    // clamp16(t >> 7) (saturation is a monotone clamp).
    debug_assert!(k_cn >= 0 && k_dn >= 0, "acelp_4t64: negative sign weights");
    for i in 0..L_SUBFR {
        let t = k_cn as i32 * cn[i] as i32 + k_dn as i32 * dn[i] as i32;
        dn2[i] = (t >> 7).clamp(-32768, 32767) as i16;
    }
    let (mut sign, mut vec) = ([0i16; L_SUBFR], [0i16; L_SUBFR]);
    for i in 0..L_SUBFR {
        let neg = dn2[i] < 0;
        sign[i] = if neg { -32768 } else { 32767 };
        vec[i] = if neg { 32767 } else { -32768 };
        dn[i] = if neg { dn[i].saturating_neg() } else { dn[i] };
        dn2[i] = if neg { dn2[i].saturating_neg() } else { dn2[i] };
    }
    // The NB_MAX best positions per track get dn2 = k - NB_MAX (k-th best), the best in
    // pos_max. The C's NB_MAX passes of first-maximum picks (every dn2 is now >= 0, picked
    // ones negative) rank by (dn2 descending, position ascending): one rank count per position.
    let mut pos_max = [0usize; NB_TRACK];
    for t in 0..NB_TRACK {
        let d: [i16; NB_POS] = std::array::from_fn(|n| dn2[t + STEP * n]);
        let mut rank = [0i16; NB_POS];
        for n in 0..NB_POS {
            for m in 0..NB_POS {
                rank[n] += (d[m] > d[n] || (d[m] == d[n] && m < n)) as i16;
            }
        }
        for n in 0..NB_POS {
            if rank[n] < NB_MAX {
                dn2[t + STEP * n] = rank[n] - NB_MAX;
            }
            if rank[n] == 0 {
                pos_max[t] = t + STEP * n;
            }
        }
    }
    // h and -h after L_SUBFR zeros (the C's h_buf), scaled by h_shift
    let mut l = 0i32;
    for &v in &hh[..L_SUBFR] {
        l = l_mac(l, v, v);
    }
    let val = extract_h(l);
    let mut h_shift = 0;
    if nb_pulse >= 12 && val > 1024 {
        h_shift = 1;
    }
    if val > 0x6000 {
        h_shift = 2;
    }
    let mut hb = [0i16; 2 * L_SUBFR];
    let mut hb_inv = [0i16; 2 * L_SUBFR];
    for i in 0..L_SUBFR {
        hb[L_SUBFR + i] = shr(hh[i], h_shift);
        hb_inv[L_SUBFR + i] = negate(hb[L_SUBFR + i]);
    }
    let h = &hb[L_SUBFR..];
    let h_abs: i64 = h.iter().map(|&v| (v as i32).abs() as i64).sum();
    // h - ix (or h_inv - ix): the pulse at ix filtered
    let pulse = |ix: usize, sign: &[i16]| -> &[i16] {
        if sign[ix] < 0 {
            &hb_inv[L_SUBFR - ix..]
        } else {
            &hb[L_SUBFR - ix..]
        }
    };
    let mut rrixix = [[0i16; NB_POS]; NB_TRACK];
    let mut rrixiy = [[0i16; MSIZE]; NB_TRACK];
    // every chain's partial sums are at most 0x8000 + 2Σh² in magnitude (Cauchy–Schwarz on
    // the diagonals): below 2^31 none saturates and the same loops run on plain adds
    let hh: i64 = h.iter().map(|&v| v as i64 * v as i64).sum();
    if 0x8000 + 2 * hh <= i32::MAX as i64 {
        cor_h::<true>(h, &mut rrixix, &mut rrixiy);
    } else {
        cor_h::<false>(h, &mut rrixix, &mut rrixiy);
    }
    // rrixiy[k][16i + j] = mult(rrixiy, sign or vec of y): 32767 for equal signs, -32768
    // otherwise; in i32 lanes with mult's clamp, rr_max per lane.
    let mut rr_max = [0i32; NB_TRACK];
    for k in 0..NB_TRACK {
        let ty = (k + 1) % NB_TRACK;
        // per row sign of x: negate where the signs differ (mult by -32768), else mult by 32767
        let sy: [bool; NB_POS] = std::array::from_fn(|j| sign[ty + STEP * j] < 0);
        let flip = [sy, sy.map(|b| !b)];
        let mut m = [0i16; NB_POS];
        for (i, row) in rrixiy[k].as_chunks_mut::<NB_POS>().0.iter_mut().enumerate() {
            sign_row(row, &flip[(sign[k + STEP * i] < 0) as usize], &mut m);
        }
        rr_max[k] = m.into_iter().max().unwrap_or(0) as i32;
    }
    let (mut psk, mut alpk) = (-1i16, 1i16);
    let mut ind = [0i16; NPMAXPT * NB_TRACK];
    let mut ipos = [0usize; NB_PULSE_MAX];
    for k in 0..nbiter {
        for i in 0..nb_pulse {
            ipos[i] = TIPOS[k * 4 + i] as usize;
        }
        let (start, mut ps, mut alp);
        if nbbits == 20 {
            start = 0;
            ps = 0;
            alp = 0;
            vec = [0; L_SUBFR];
        } else if nbbits == 36 || nbbits == 44 {
            start = 2;
            let (ix, iy) = (pos_max[ipos[0]], pos_max[ipos[1]]);
            ind[0] = ix as i16;
            ind[1] = iy as i16;
            ps = add(dn[ix], dn[iy]);
            let (i, j) = (ix >> 2, iy >> 2);
            let mut s = l_mult(rrixix[ipos[0]][i], 4096);
            s = l_mac(s, rrixix[ipos[1]][j], 4096);
            s = l_mac(s, rrixiy[ipos[0]][(i << 4) + j], 8192);
            alp = round(s);
            let (p0, p1) = (pulse(ix, &sign), pulse(iy, &sign));
            for i in 0..L_SUBFR {
                vec[i] = add(p0[i], p1[i]);
            }
            if nbbits == 44 {
                ipos[8] = 0;
                ipos[9] = 1;
            }
        } else {
            start = 4;
            let p: [usize; 4] = std::array::from_fn(|n| pos_max[ipos[n]]);
            for n in 0..4 {
                ind[n] = p[n] as i16;
            }
            ps = add(add(add(dn[p[0]], dn[p[1]]), dn[p[2]]), dn[p[3]]);
            let q: [&[i16]; 4] = std::array::from_fn(|n| pulse(p[n], &sign));
            for i in 0..L_SUBFR {
                vec[i] = add(add(add(q[0][i], q[1][i]), q[2][i]), q[3][i]);
            }
            let mut l = 0i32;
            for &v in &vec {
                l = l_mac(l, v, v);
            }
            alp = round(l_shr(l, 3));
            if nbbits == 72 {
                ipos[16] = 0;
                ipos[17] = 1;
            }
        }
        let mut st = 0;
        let mut j = start;
        while j < nb_pulse {
            let (mut cor_x, mut cor_y) = ([0i16; NB_POS], [0i16; NB_POS]);
            let fast = HVec::new(&vec, h_abs);
            cor_h_vec(h, &fast, &vec, ipos[j], &sign, &rrixix, &mut cor_x);
            cor_h_vec(h, &fast, &vec, ipos[j + 1], &sign, &rrixix, &mut cor_y);
            let (ix, iy) = search_ixiy(
                nbpos[st],
                ipos[j],
                ipos[j + 1],
                &mut ps,
                &mut alp,
                &dn[..],
                &dn2,
                &cor_x,
                &cor_y,
                &rrixiy,
                &rr_max,
            );
            ind[j] = ix as i16;
            ind[j + 1] = iy as i16;
            let (p0, p1) = (pulse(ix, &sign), pulse(iy, &sign));
            for i in 0..L_SUBFR {
                vec[i] = add(vec[i], add(p0[i], p1[i]));
            }
            j += 2;
            st += 1;
        }
        let ps = mult(ps, ps);
        let s = l_msu(l_mult(alpk, ps), psk, alp);
        if s > 0 {
            psk = ps;
            alpk = alp;
            codvec[..nb_pulse].copy_from_slice(&ind[..nb_pulse]);
            *y = vec;
        }
    }
    let mut ind = [-1i16; NPMAXPT * NB_TRACK];
    code.fill(0);
    for v in y.iter_mut() {
        *v = shr_r(*v, 3);
    }
    let val = shr(512, h_shift);
    for k in 0..nb_pulse {
        let i = codvec[k] as usize;
        let mut index = (i >> 2) as i16;
        let track = i & 3;
        if sign[i] > 0 {
            code[i] = add(code[i], val);
            codvec[k] += 2 * L_SUBFR as i16;
        } else {
            code[i] = sub(code[i], val);
            index += NB_POS as i16;
        }
        let mut n = track * NPMAXPT;
        while ind[n] >= 0 {
            n += 1;
        }
        ind[n] = index;
    }
    let mut k = 0;
    let lo = |v: i32| v as i16;
    match nbbits {
        20 => {
            for t in 0..NB_TRACK {
                index_out[t] = lo(quant_1p_n1(ind[k], 4));
                k += NPMAXPT;
            }
        }
        36 => {
            for t in 0..NB_TRACK {
                index_out[t] = lo(quant_2p_2n1(ind[k], ind[k + 1], 4));
                k += NPMAXPT;
            }
        }
        44 => {
            for t in 0..NB_TRACK {
                index_out[t] = lo(if t < 2 {
                    quant_3p_3n1(ind[k], ind[k + 1], ind[k + 2], 4)
                } else {
                    quant_2p_2n1(ind[k], ind[k + 1], 4)
                });
                k += NPMAXPT;
            }
        }
        52 => {
            for t in 0..NB_TRACK {
                index_out[t] = lo(quant_3p_3n1(ind[k], ind[k + 1], ind[k + 2], 4));
                k += NPMAXPT;
            }
        }
        64 => {
            for t in 0..NB_TRACK {
                let l = quant_4p_4n(&ind[k..k + 4], 4);
                index_out[t] = lo((l >> 14) & 3);
                index_out[t + NB_TRACK] = lo(l & 0x3fff);
                k += NPMAXPT;
            }
        }
        72 => {
            for t in 0..NB_TRACK {
                if t < 2 {
                    let l = quant_5p_5n(&ind[k..k + 5], 4);
                    index_out[t] = lo((l >> 10) & 0x3ff);
                    index_out[t + NB_TRACK] = lo(l & 0x3ff);
                } else {
                    let l = quant_4p_4n(&ind[k..k + 4], 4);
                    index_out[t] = lo((l >> 14) & 3);
                    index_out[t + NB_TRACK] = lo(l & 0x3fff);
                }
                k += NPMAXPT;
            }
        }
        88 => {
            for t in 0..NB_TRACK {
                let l = quant_6p_6n_2(&ind[k..k + 6], 4);
                index_out[t] = lo((l >> 11) & 0x7ff);
                index_out[t + NB_TRACK] = lo(l & 0x7ff);
                k += NPMAXPT;
            }
        }
        _ => {}
    }
}

/// The C's `rrixix` and `rrixiy` chains; `PLAIN`: their adds without saturation (exact when
/// the caller's bound holds).
#[inline(always)]
fn cor_h<const PLAIN: bool>(
    h: &[i16],
    rrixix: &mut [[i16; NB_POS]; NB_TRACK],
    rrixiy: &mut [[i16; MSIZE]; NB_TRACK],
) {
    let mac = |c: i32, a: i16, b: i16| {
        if PLAIN {
            c.wrapping_add((a as i32 * b as i32).wrapping_mul(2))
        } else {
            l_mac(c, a, b)
        }
    };
    let mut cor = 0x8000i32;
    for i in 0..NB_POS {
        for t in (0..NB_TRACK).rev() {
            let v = h[4 * i + 3 - t];
            cor = mac(cor, v, v);
            rrixix[t][NB_POS - 1 - i] = extract_h(cor);
        }
    }
    for k in 0..NB_POS {
        let mut p = (MSIZE - 1 - k * NB_POS) as isize;
        let mut cor = 0x8000i32;
        let (mut a, mut b) = (0usize, 1 + STEP * k);
        for _ in k + 1..NB_POS {
            for t in [2usize, 1, 0] {
                cor = mac(cor, h[a], h[b]);
                a += 1;
                b += 1;
                rrixiy[t][p as usize] = extract_h(cor);
            }
            cor = mac(cor, h[a], h[b]);
            a += 1;
            b += 1;
            rrixiy[3][(p - NB_POS as isize) as usize] = extract_h(cor);
            p -= NB_POS as isize + 1;
        }
        for t in [2usize, 1, 0] {
            cor = mac(cor, h[a], h[b]);
            a += 1;
            b += 1;
            rrixiy[t][p as usize] = extract_h(cor);
        }
    }
    let mut pos = MSIZE - 1;
    for k in 0..NB_POS {
        let (mut p3, mut p) = (pos as isize, pos as isize - 1);
        let mut cor = 0x8000i32;
        let (mut a, mut b) = (0usize, 3 + STEP * k);
        for _ in k + 1..NB_POS {
            cor = mac(cor, h[a], h[b]);
            a += 1;
            b += 1;
            rrixiy[3][p3 as usize] = extract_h(cor);
            for t in [2usize, 1, 0] {
                cor = mac(cor, h[a], h[b]);
                a += 1;
                b += 1;
                rrixiy[t][p as usize] = extract_h(cor);
            }
            p3 -= NB_POS as isize + 1;
            p -= NB_POS as isize + 1;
        }
        cor = mac(cor, h[a], h[b]);
        rrixiy[3][p3 as usize] = extract_h(cor);
        pos -= 1;
    }
}

/// One row of `acelp_4t64`'s sign pass: `mult(c, -32768) = sat(-c)` where `neg`, else
/// `mult(c, 32767) = c - (c > 0) + (c == -32768)`; the running |row| maxima in `m` (results are
/// within ±32767). In i16 lanes: widened to i32, LLVM vectorised across the unrolled rows
/// with strided gathers instead.
#[inline(always)]
fn sign_row(row: &mut [i16; NB_POS], neg: &[bool; NB_POS], m: &mut [i16; NB_POS]) {
    for j in 0..NB_POS {
        let c = row[j];
        let pos = c - (c > 0) as i16 + (c == i16::MIN) as i16;
        let v = if neg[j] { c.saturating_neg() } else { pos };
        row[j] = v;
        m[j] = m[j].max(v.abs());
    }
}

/// `vec` zero-padded to two subframes, if `cor_h_vec`'s chains can't saturate on it:
/// every partial sum of `Σ 2·h[j-pos]·vec[j]` is at most `2·Σ|h|·max|vec|` in magnitude, so
/// below 2^31 the saturating chain is the plain sum (in any order).
struct HVec(Option<[i16; 2 * L_SUBFR]>);

impl HVec {
    #[inline(always)]
    fn new(vec: &[i16; L_SUBFR], h_abs: i64) -> Self {
        let max = vec.iter().map(|&v| (v as i32).abs()).max().unwrap_or(0) as i64;
        if 2 * h_abs * max >= 1 << 31 {
            return Self(None);
        }
        let mut p = [0i16; 2 * L_SUBFR];
        p[..L_SUBFR].copy_from_slice(vec);
        Self(Some(p))
    }
}

/// `cor_h_vec`: the correlation of `h` with `vec` at the track's positions, signed, plus
/// the track's `rrixix` (the plain fixed-length sums when `fast` allows, else the C's chain).
#[inline(always)]
fn cor_h_vec(
    h: &[i16],
    fast: &HVec,
    vec: &[i16; L_SUBFR],
    track: usize,
    sign: &[i16; L_SUBFR],
    rrixix: &[[i16; NB_POS]; NB_TRACK],
    cor: &mut [i16; NB_POS],
) {
    let h: &[i16; L_SUBFR] = h[..L_SUBFR].try_into().unwrap();
    let mut l = [0i32; NB_POS];
    match &fast.0 {
        // Position `pos` has `64 - pos` terms before the padding, so positions 4g..4g + 4 need
        // at most `64 - 16g`: each group's dots run only that far (fixed lengths, exact).
        Some(p) => {
            cor_h_vec_group::<64>(h, p, track, 0, &mut l);
            cor_h_vec_group::<48>(h, p, track, 1, &mut l);
            cor_h_vec_group::<32>(h, p, track, 2, &mut l);
            cor_h_vec_group::<16>(h, p, track, 3, &mut l);
        }
        None => {
            for (i, l) in l.iter_mut().enumerate() {
                let pos = track + STEP * i;
                for j in pos..L_SUBFR {
                    *l = l_mac(*l, h[j - pos], vec[j]);
                }
            }
        }
    }
    for i in 0..NB_POS {
        let pos = track + STEP * i;
        let corr = round(l_shl(l[i], 1));
        cor[i] = add(mult(corr, sign[pos]), rrixix[track][i]);
    }
}

/// `cor_h_vec`'s plain dots `2·Σ_{j<N} h[j]·p[pos + j]` for the 4 positions `i = 4g..4g + 4`
/// (`pos = track + 4i`); `N >= 64 - 16g` covers every term before `p`'s zero padding.
#[inline(always)]
fn cor_h_vec_group<const N: usize>(
    h: &[i16; L_SUBFR],
    p: &[i16; 2 * L_SUBFR],
    track: usize,
    g: usize,
    l: &mut [i32; NB_POS],
) {
    const { assert!(N <= L_SUBFR) };
    debug_assert!(
        N + 16 * g >= L_SUBFR,
        "cor_h_vec_group: the dots stop before the padding"
    );
    let h: &[i16; N] = h[..N].try_into().unwrap();
    for i in 4 * g..4 * g + 4 {
        let pos = track + STEP * i;
        let v: &[i16; N] = p[pos..pos + N].try_into().unwrap();
        let mut s = 0i32;
        for j in 0..N {
            s = s.wrapping_add(h[j] as i32 * v[j] as i32);
        }
        l[i] = s << 1;
    }
}

/// `search_ixiy`: the best pair of positions on tracks `track_x`, `track_y` (pulse 1 among
/// the `nb_pos_ix` preselected positions); updates `ps` and `alp`.
///
/// The C's chain in its order, its saturation removed where exact: the comparison
/// `L_msu(L_mult(alpk, sq), sqk, alp_16) > 0` is `alpk·sq > sqk·alp_16` (both products at
/// most 2^30: i32; `sq >= 0` and `sqk >= -1` keep each `L_mult` exact, saturation keeps the
/// sign), and `ps2`, `alp2` are plain sums when the bounds below hold (else the chain).
#[inline(always)]
#[allow(clippy::too_many_arguments)] // the C's
fn search_ixiy(
    nb_pos_ix: i16,
    track_x: usize,
    track_y: usize,
    ps: &mut i16,
    alp: &mut i16,
    dn: &[i16],
    dn2: &[i16; L_SUBFR],
    cor_x: &[i16; NB_POS],
    cor_y: &[i16; NB_POS],
    rrixiy: &[[i16; MSIZE]; NB_TRACK],
    rr_max: &[i32; NB_TRACK],
) -> (usize, usize) {
    let p2 = &rrixiy[track_x];
    let thres_ix = sub(nb_pos_ix, NB_MAX);
    let alp0 = l_add(l_deposit_h(*alp), 0x8000);
    let amax = |v: &[i16; NB_POS]| v.iter().map(|&c| (c as i32).abs()).max().unwrap_or(0);
    let dx: [i16; NB_POS] = std::array::from_fn(|n| dn[track_x + STEP * n]);
    let dy: [i16; NB_POS] = std::array::from_fn(|n| dn[track_y + STEP * n]);
    // every partial ps sum within |ps| + max|dx| + max|dy| (no i16 saturation, ps2 != -32768),
    // every partial alp sum within |alp0| + 8192·(max|cor_x| + max|cor_y|) + 16384·max|p2|
    let ps_ok = (*ps as i32).abs() + amax(&dx) + amax(&dy) <= 32767;
    let alp_ok = (alp0 as i64).abs()
        + 8192 * (amax(cor_x) + amax(cor_y)) as i64
        + 16384 * (rr_max[track_x] as i64)
        < 1 << 31;
    if !(ps_ok && alp_ok) {
        return search_ixiy_chain(
            nb_pos_ix, track_x, track_y, ps, alp, dn, dn2, cor_x, cor_y, rrixiy,
        );
    }
    if let Some((best, alpk)) = search_ixiy_lanes(
        thres_ix, track_x, *ps, alp0, dn2, &dx, &dy, cor_x, cor_y, p2,
    ) {
        let (ix, iy) = best.map_or((track_x, track_y), |(nx, ny)| {
            (track_x + STEP * nx, track_y + STEP * ny)
        });
        *ps = add(*ps, add(dn[ix], dn[iy]));
        *alp = alpk as i16;
        return (ix, iy);
    }
    let (mut sqk, mut alpk) = (-1i32, 1i32);
    let (mut ix, mut iy) = (track_x, track_y);
    for nx in 0..NB_POS {
        if dn2[track_x + STEP * nx] >= thres_ix {
            continue;
        }
        let ps1 = *ps as i32 + dx[nx] as i32;
        let alp1 = alp0 + ((cor_x[nx] as i32) << 13);
        let row: &[i16; NB_POS] = p2[NB_POS * nx..NB_POS * (nx + 1)].try_into().unwrap();
        let mut pos = None;
        for ny in 0..NB_POS {
            let ps2 = ps1 + dy[ny] as i32;
            let alp2 = alp1 + ((cor_y[ny] as i32) << 13) + ((row[ny] as i32) << 14);
            let alp_16 = alp2 >> 16;
            let sq = (ps2 * ps2) >> 15;
            if alpk * sq > sqk * alp_16 {
                sqk = sq;
                alpk = alp_16;
                pos = Some(ny);
            }
        }
        if let Some(p) = pos {
            ix = track_x + STEP * nx;
            iy = track_y + STEP * p;
        }
    }
    *ps = add(*ps, add(dn[ix], dn[iy]));
    *alp = alpk as i16;
    (ix, iy)
}

/// `search_ixiy`'s plain chain (the caller's bounds hold) in lanes: with every `alp_16 > 0`
/// the chain keeps the first maximum of `sq/alp_16` in (x, y) order (from `-1/1`, any
/// candidate wins first), so each row of 16 y lanes is computed at once, each lane keeps its
/// first best over the rows, and the lanes reduce in order (ties: the lower x, then the lower
/// y). `None` if some `alp_16 <= 0` (then the scalar chain). Returns the winner's (x, y)
/// indices on their tracks (`None`: no candidate) and `alpk`.
#[inline(always)]
#[allow(clippy::too_many_arguments)]
fn search_ixiy_lanes(
    thres_ix: i16,
    track_x: usize,
    ps: i16,
    alp0: i32,
    dn2: &[i16; L_SUBFR],
    dx: &[i16; NB_POS],
    dy: &[i16; NB_POS],
    cor_x: &[i16; NB_POS],
    cor_y: &[i16; NB_POS],
    p2: &[i16; MSIZE],
) -> Option<(Option<(usize, usize)>, i32)> {
    let (mut bsq, mut bal, mut bx) = ([-1i32; NB_POS], [1i32; NB_POS], [NB_POS; NB_POS]);
    for nx in 0..NB_POS {
        if dn2[track_x + STEP * nx] >= thres_ix {
            continue;
        }
        let ps1 = ps as i32 + dx[nx] as i32;
        let alp1 = alp0 + ((cor_x[nx] as i32) << 13);
        let row: &[i16; NB_POS] = p2[NB_POS * nx..NB_POS * (nx + 1)].try_into().unwrap();
        let (mut sq, mut al) = ([0i32; NB_POS], [0i32; NB_POS]);
        for n in 0..NB_POS {
            let ps2 = ps1 + dy[n] as i32;
            al[n] = (alp1 + ((cor_y[n] as i32) << 13) + ((row[n] as i32) << 14)) >> 16;
            sq[n] = (ps2 * ps2) >> 15;
        }
        if al.iter().any(|&a| a <= 0) {
            return None;
        }
        for n in 0..NB_POS {
            // 15-bit, al > 0: the cross products fit an i32
            let better = sq[n] * bal[n] > bsq[n] * al[n];
            bsq[n] = if better { sq[n] } else { bsq[n] };
            bal[n] = if better { al[n] } else { bal[n] };
            bx[n] = if better { nx } else { bx[n] };
        }
    }
    let (mut sqk, mut alpk, mut kx, mut ky) = (-1i32, 1i32, NB_POS, NB_POS);
    for n in 0..NB_POS {
        if bx[n] == NB_POS {
            continue;
        }
        let (l, r) = (bsq[n] * alpk, sqk * bal[n]);
        if l > r || l == r && bx[n] < kx {
            (sqk, alpk, kx, ky) = (bsq[n], bal[n], bx[n], n);
        }
    }
    if kx == NB_POS {
        return Some((None, 1));
    }
    Some((Some((kx, ky)), alpk))
}

/// `search_ixiy` as the C's chain.
#[allow(clippy::too_many_arguments)] // the C's
fn search_ixiy_chain(
    nb_pos_ix: i16,
    track_x: usize,
    track_y: usize,
    ps: &mut i16,
    alp: &mut i16,
    dn: &[i16],
    dn2: &[i16; L_SUBFR],
    cor_x: &[i16; NB_POS],
    cor_y: &[i16; NB_POS],
    rrixiy: &[[i16; MSIZE]; NB_TRACK],
) -> (usize, usize) {
    let p2 = &rrixiy[track_x];
    let thres_ix = sub(nb_pos_ix, NB_MAX);
    let alp0 = l_add(l_deposit_h(*alp), 0x8000);
    let (mut sqk, mut alpk) = (-1i16, 1i16);
    let (mut ix, mut iy) = (track_x, track_y);
    let mut q = 0;
    for (nx, x) in (track_x..L_SUBFR).step_by(STEP).enumerate() {
        let ps1 = add(*ps, dn[x]);
        let alp1 = l_mac(alp0, cor_x[nx], 4096);
        if dn2[x] < thres_ix {
            let mut pos = None;
            for (ny, y) in (track_y..L_SUBFR).step_by(STEP).enumerate() {
                let ps2 = add(ps1, dn[y]);
                let alp2 = l_mac(l_mac(alp1, cor_y[ny], 4096), p2[q], 8192);
                q += 1;
                let alp_16 = extract_h(alp2);
                let sq = mult(ps2, ps2);
                let s = l_msu(l_mult(alpk, sq), sqk, alp_16);
                if s > 0 {
                    sqk = sq;
                    alpk = alp_16;
                    pos = Some(y);
                }
            }
            if let Some(p) = pos {
                ix = x;
                iy = p;
            }
        } else {
            q += NB_POS;
        }
    }
    *ps = add(*ps, add(dn[ix], dn[iy]));
    *alp = alpk;
    (ix, iy)
}

/// `quant_1p_N1`: one pulse, N+1 bits.
pub(crate) fn quant_1p_n1(pos: i16, n: i16) -> i32 {
    let mask = (1i16 << n) - 1;
    let mut index = (pos & mask) as i32;
    if pos & NB_POS as i16 != 0 {
        index += 1 << n;
    }
    index
}

/// `quant_2p_2N1`: two pulses, 2N+1 bits.
pub(crate) fn quant_2p_2n1(pos1: i16, pos2: i16, n: i16) -> i32 {
    let mask = (1i16 << n) - 1;
    let nbp = NB_POS as i16;
    let (a, b) = ((pos1 & mask) as i32, (pos2 & mask) as i32);
    if ((pos2 ^ pos1) & nbp) == 0 {
        let mut index = if pos1 <= pos2 {
            (a << n) + b
        } else {
            (b << n) + a
        };
        if pos1 & nbp != 0 {
            index += 1 << (2 * n);
        }
        index
    } else if a <= b {
        let mut index = (b << n) + a;
        if pos2 & nbp != 0 {
            index += 1 << (2 * n);
        }
        index
    } else {
        let mut index = (a << n) + b;
        if pos1 & nbp != 0 {
            index += 1 << (2 * n);
        }
        index
    }
}

/// `quant_3p_3N1`: three pulses, 3N+1 bits.
pub(crate) fn quant_3p_3n1(pos1: i16, pos2: i16, pos3: i16, n: i16) -> i32 {
    let nb_pos = 1i16 << (n - 1);
    let (p, q, r) = if ((pos1 ^ pos2) & nb_pos) == 0 {
        (pos1, pos2, pos3)
    } else if ((pos1 ^ pos3) & nb_pos) == 0 {
        (pos1, pos3, pos2)
    } else {
        (pos2, pos3, pos1)
    };
    quant_2p_2n1(p, q, n - 1) + (((p & nb_pos) as i32) << n) + (quant_1p_n1(r, n) << (2 * n))
}

/// `quant_4p_4N1`: four pulses, 4N+1 bits.
pub(crate) fn quant_4p_4n1(pos1: i16, pos2: i16, pos3: i16, pos4: i16, n: i16) -> i32 {
    let nb_pos = 1i16 << (n - 1);
    let (p, q, r, s) = if ((pos1 ^ pos2) & nb_pos) == 0 {
        (pos1, pos2, pos3, pos4)
    } else if ((pos1 ^ pos3) & nb_pos) == 0 {
        (pos1, pos3, pos2, pos4)
    } else {
        (pos2, pos3, pos1, pos4)
    };
    quant_2p_2n1(p, q, n - 1) + (((p & nb_pos) as i32) << n) + (quant_2p_2n1(r, s, n) << (2 * n))
}

/// Splits `pos` by the bit `nb_pos`: (A without it, B with it).
fn split(pos: &[i16], nb_pos: i16) -> ([i16; 6], usize, [i16; 6], usize) {
    let (mut a, mut b, mut i, mut j) = ([0; 6], [0; 6], 0, 0);
    for &p in pos {
        if p & nb_pos == 0 {
            a[i] = p;
            i += 1;
        } else {
            b[j] = p;
            j += 1;
        }
    }
    (a, i, b, j)
}

/// `quant_4p_4N`: four pulses, 4N bits.
pub(crate) fn quant_4p_4n(pos: &[i16], n: i16) -> i32 {
    let n_1 = n - 1;
    let (a, i, b, _) = split(&pos[..4], 1 << n_1);
    let index = match i {
        0 => (1 << (4 * n - 3)) + quant_4p_4n1(b[0], b[1], b[2], b[3], n_1),
        1 => (quant_1p_n1(a[0], n_1) << (3 * n_1 + 1)) + quant_3p_3n1(b[0], b[1], b[2], n_1),
        2 => (quant_2p_2n1(a[0], a[1], n_1) << (2 * n_1 + 1)) + quant_2p_2n1(b[0], b[1], n_1),
        3 => (quant_3p_3n1(a[0], a[1], a[2], n_1) << n) + quant_1p_n1(b[0], n_1),
        _ => quant_4p_4n1(a[0], a[1], a[2], a[3], n_1),
    };
    index + (((i & 3) as i32) << (4 * n - 2))
}

/// `quant_5p_5N`: five pulses, 5N bits.
pub(crate) fn quant_5p_5n(pos: &[i16], n: i16) -> i32 {
    let n_1 = n - 1;
    let (a, i, b, _) = split(&pos[..5], 1 << n_1);
    let top = 1 << (5 * n - 1);
    let q3 = |p: &[i16; 6]| quant_3p_3n1(p[0], p[1], p[2], n_1) << (2 * n + 1);
    match i {
        0 => top + q3(&b) + quant_2p_2n1(b[3], b[4], n),
        1 => top + q3(&b) + quant_2p_2n1(b[3], a[0], n),
        2 => top + q3(&b) + quant_2p_2n1(a[0], a[1], n),
        3 => q3(&a) + quant_2p_2n1(b[0], b[1], n),
        4 => q3(&a) + quant_2p_2n1(a[3], b[0], n),
        _ => q3(&a) + quant_2p_2n1(a[3], a[4], n),
    }
}

/// `quant_6p_6N_2`: six pulses, 6N-2 bits.
pub(crate) fn quant_6p_6n_2(pos: &[i16], n: i16) -> i32 {
    let n_1 = n - 1;
    let (a, i, b, _) = split(&pos[..6], 1 << n_1);
    let top = 1 << (6 * n - 5);
    let (index, i) = match i {
        0 => (
            top + (quant_5p_5n(&b, n_1) << n) + quant_1p_n1(b[5], n_1),
            0,
        ),
        1 => (
            top + (quant_5p_5n(&b, n_1) << n) + quant_1p_n1(a[0], n_1),
            1,
        ),
        2 => (
            top + (quant_4p_4n(&b, n_1) << (2 * n_1 + 1)) + quant_2p_2n1(a[0], a[1], n_1),
            2,
        ),
        3 => (
            (quant_3p_3n1(a[0], a[1], a[2], n_1) << (3 * n_1 + 1))
                + quant_3p_3n1(b[0], b[1], b[2], n_1),
            3,
        ),
        4 => (
            (quant_4p_4n(&a, n_1) << (2 * n_1 + 1)) + quant_2p_2n1(b[0], b[1], n_1),
            2,
        ),
        5 => ((quant_5p_5n(&a, n_1) << n) + quant_1p_n1(b[0], n_1), 1),
        _ => ((quant_5p_5n(&a, n_1) << n) + quant_1p_n1(a[5], n_1), 0),
    };
    index + ((i & 3) << (6 * n - 4))
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::encoder::acelp2::vs_ref::search_input;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_ACELP_4t64_fx(
            dn: *mut i16,
            cn: *mut i16,
            h: *mut i16,
            code: *mut i16,
            y: *mut i16,
            nbbits: i16,
            ser_size: i16,
            index: *mut i16,
        );
        fn refwb_quant_4p_4N(pos: *mut i16, n: i16) -> i32;
        fn refwb_quant_5p_5N(pos: *mut i16, n: i16) -> i32;
        fn refwb_quant_6p_6N_2(pos: *mut i16, n: i16) -> i32;
        fn refwb_quant_3p_3N1(p1: i16, p2: i16, p3: i16, n: i16) -> i32;
    }

    /// Every bit allocation (both 88-bit iteration counts), impulse responses from quiet to
    /// loud (all three `h_shift`s, CR 0036's 2 included), wild ones, zero residuals.
    #[test]
    fn acelp_4t64_matches_reference() {
        let mut r = Rng::new(91);
        let cases = [
            (20, 177),
            (36, 253),
            (44, 285),
            (52, 317),
            (64, 365),
            (72, 397),
            (88, 461),
            (88, 477),
        ];
        for n in 0..16_000u64 {
            let (mut h, x, mut cn) = search_input(&mut r, n);
            if n % 3 == 1 {
                for v in &mut h {
                    *v = v.saturating_mul(1 + (r.u64() % 6) as i16);
                }
            }
            if n % 61 == 7 {
                // full-scale alternating h: the rrixiy chains saturate negative (-32768 entries)
                for (i, v) in h.iter_mut().enumerate() {
                    *v = if i % 2 == 0 {
                        32767
                    } else {
                        -32768 + (r.u64() % 64) as i16
                    };
                }
            }
            let (nb, ser) = cases[(n % 8) as usize];
            let mut dn = [0i16; L_SUBFR];
            crate::encoder::acelp2::cor_h_x(&h, &x, &mut dn);
            let mut dnc = dn;
            // y: the caller's buffer, kept if no candidate ever improves
            let mut y: [i16; L_SUBFR] = std::array::from_fn(|_| r.i16());
            let (mut code, mut ix) = ([0i16; L_SUBFR], [0i16; 8]);
            let (mut codec, mut yc, mut ixc) = ([0i16; L_SUBFR], y, [0i16; 8]);
            acelp_4t64_fx(&mut dn, &cn, &h, &mut code, &mut y, nb, ser, &mut ix);
            with_ref(|| {
                // SAFETY: 64-sample vectors, 8 indices.
                unsafe {
                    refwb_ACELP_4t64_fx(
                        dnc.as_mut_ptr(),
                        cn.as_mut_ptr(),
                        h.as_mut_ptr(),
                        codec.as_mut_ptr(),
                        yc.as_mut_ptr(),
                        nb,
                        ser,
                        ixc.as_mut_ptr(),
                    )
                }
            });
            assert_eq!(
                (dn, code, y, ix),
                (dnc, codec, yc, ixc),
                "ACELP_4t64_fx {nb} {n}"
            );
        }
    }

    /// The pulse index coders on random position sets (any signs, repeats).
    #[test]
    fn quant_pulses_match_reference() {
        let mut r = Rng::new(92);
        for _ in 0..100_000 {
            let mut p: [i16; 6] = std::array::from_fn(|_| (r.u64() % 32) as i16);
            let k = 2 + (r.u64() % 3) as i16;
            let q: [i16; 6] = std::array::from_fn(|i| p[i] & ((2 << k) - 1) | (p[i] & 16));
            let mut qq = q;
            with_ref(|| {
                // SAFETY: up to 6 positions.
                unsafe {
                    assert_eq!(quant_4p_4n(&p, 4), refwb_quant_4p_4N(p.as_mut_ptr(), 4));
                    assert_eq!(quant_5p_5n(&p, 4), refwb_quant_5p_5N(p.as_mut_ptr(), 4));
                    assert_eq!(quant_6p_6n_2(&p, 4), refwb_quant_6p_6N_2(p.as_mut_ptr(), 4));
                    assert_eq!(
                        quant_3p_3n1(q[0], q[1], q[2], k),
                        refwb_quant_3p_3N1(qq[0], qq[1], qq[2], k)
                    );
                    let _ = &mut qq;
                }
            });
        }
    }
}
