//! The target–impulse-response correlation (vo-amrwbenc `cor_h_x.c`) and the 2-track
//! algebraic codebook search of 6.60 kbit/s (`c2t64fx.c`), with the reference's saturation.

use super::ops::*;
use crate::basic_op::dot_product12;
use crate::cnst::L_SUBFR;
pub(crate) use crate::simd::{acelp_2t64_fx, cor_h_x};

/// `cor_h_x`: `dn[i] = Σ x[j]·h[j-i]`, scaled together so the codebook search stays in 12 bits.
#[inline(always)]
pub(crate) fn cor_h_x_portable(h: &[i16], x: &[i16], dn: &mut [i16; L_SUBFR]) {
    let y32 = corr_hx(
        h[..L_SUBFR].try_into().unwrap(),
        x[..L_SUBFR].try_into().unwrap(),
    );
    let mut l_tot = 1i32;
    for k in 0..4 {
        let mut l_max = 0i32;
        for i in (k..L_SUBFR).step_by(4) {
            l_max = l_max.max(l_abs(y32[i]));
        }
        let l_max = l_shr(l_max, 2);
        l_tot = l_add(l_tot, l_max);
        l_tot = l_add(l_tot, l_shr(l_max, 1));
    }
    let j = sub(norm_l(l_tot), 4);
    for i in 0..L_SUBFR {
        dn[i] = round(l_shl(y32[i], j));
    }
}

/// `cor_h_x`'s chains `1 + Σ_{j>=i} 2·x[j]·h[j-i]`: plain fixed-length dots over x
/// zero-padded to two subframes when exact (every partial sum is at most
/// `sqrt(Σx²·Σh²)` in magnitude, Cauchy–Schwarz: below 2^30 no step saturates), else the chain.
#[inline(always)]
fn corr_hx(h: &[i16; L_SUBFR], x: &[i16; L_SUBFR]) -> [i32; L_SUBFR] {
    let e = |v: &[i16; L_SUBFR]| v.iter().map(|&s| s as i64 * s as i64).sum::<i64>();
    const B: i128 = (1 << 30) - 1;
    if e(x) as i128 * e(h) as i128 <= B * B {
        let mut xp = [0i16; 2 * L_SUBFR];
        xp[..L_SUBFR].copy_from_slice(x);
        return std::array::from_fn(|i| {
            let w: &[i16; L_SUBFR] = xp[i..i + L_SUBFR].try_into().unwrap();
            let mut s = 0i32;
            for j in 0..L_SUBFR {
                s = s.wrapping_add(h[j] as i32 * w[j] as i32);
            }
            1 + (s << 1)
        });
    }
    std::array::from_fn(|i| {
        let mut l = 1i32;
        for j in i..L_SUBFR {
            l = l_mac(l, x[j], h[j - i]);
        }
        l
    })
}

const NB_POS2: usize = 32;
const MSIZE2: usize = 1024;

/// `ACELP_2t64_fx`: two pulses on two interleaved tracks of 32 positions (12-bit index);
/// `code` the codevector (Q9), `y` it filtered (Q9). `dn` gets the signs folded in.
#[inline(always)]
pub(crate) fn acelp_2t64_fx_portable(
    dn: &mut [i16; L_SUBFR],
    cn: &[i16],
    hh: &[i16],
    code: &mut [i16; L_SUBFR],
    y: &mut [i16; L_SUBFR],
    index: &mut i16,
) {
    let alp = 8192;
    let mut exp = 0;
    let mut s = dot_product12(cn, cn, &mut exp);
    isqrt_n(&mut s, &mut exp);
    let k_cn = round(l_shl(s, add(exp, 5)));
    let mut s = dot_product12(&dn[..], &dn[..], &mut exp);
    isqrt_n(&mut s, &mut exp);
    let k_dn = mult_r(alp, round(l_shl(s, add(exp, 5 + 3))));
    let mut dn2 = [0i16; L_SUBFR];
    for i in 0..L_SUBFR {
        let s = l_mac(l_mult(k_cn, cn[i]), k_dn, dn[i]);
        dn2[i] = extract_h(l_shl(s, 8));
    }
    let (mut sign, mut vec) = ([0i16; L_SUBFR], [0i16; L_SUBFR]);
    for i in 0..L_SUBFR {
        if dn2[i] >= 0 {
            sign[i] = 32767;
            vec[i] = -32768;
        } else {
            sign[i] = -32768;
            vec[i] = 32767;
            dn[i] = negate(dn[i]);
        }
    }
    // h_buf: [zeros, h, (unused), -h] as the C's h_buf[4·L_SUBFR]; code vectors read
    // h - ix (zeros before h) or h_inv - ix
    let mut h = [0i16; 2 * L_SUBFR];
    let mut h_inv = [0i16; 2 * L_SUBFR];
    for i in 0..L_SUBFR {
        h[L_SUBFR + i] = hh[i];
        h_inv[L_SUBFR + i] = negate(hh[i]);
    }
    let hp = &h[L_SUBFR..];
    let mut rrixix = [[0i16; NB_POS2]; 2];
    let mut rrixiy = [0i16; MSIZE2];
    // every chain's partial sums are at most 0x10000 + 2Σh² in magnitude (Cauchy–Schwarz on
    // the diagonals): below 2^31 none saturates and the same loops run on plain adds
    let hh: i64 = hp.iter().map(|&v| v as i64 * v as i64).sum();
    if 0x1_0000 + 2 * hh <= i32::MAX as i64 {
        cor_h2::<true>(hp, &mut rrixix, &mut rrixiy);
    } else {
        cor_h2::<false>(hp, &mut rrixix, &mut rrixiy);
    }
    // the odd positions' sign and its inverse (vec) contiguous: rows of 32 lanes
    let sy: [i16; NB_POS2] = std::array::from_fn(|n| sign[2 * n + 1]);
    let vy: [i16; NB_POS2] = std::array::from_fn(|n| vec[2 * n + 1]);
    for (n, row) in rrixiy.as_chunks_mut::<NB_POS2>().0.iter_mut().enumerate() {
        let psign = if sign[2 * n] < 0 { &vy } else { &sy };
        for (r, &s) in row.iter_mut().zip(psign) {
            *r = mult(*r, s);
        }
    }
    let (ix, iy) =
        search_2t(dn, &rrixix, &rrixiy).unwrap_or_else(|| search_2t_chain(dn, &rrixix, &rrixiy));
    code.fill(0);
    let (mut i0, mut i1) = ((ix >> 1) as i16, (iy >> 1) as i16);
    let p0: &[i16] = if sign[ix] > 0 {
        code[ix] = 512;
        &h[L_SUBFR - ix..]
    } else {
        code[ix] = -512;
        i0 += NB_POS2 as i16;
        &h_inv[L_SUBFR - ix..]
    };
    let p1: &[i16] = if sign[iy] > 0 {
        code[iy] = 512;
        &h[L_SUBFR - iy..]
    } else {
        code[iy] = -512;
        i1 += NB_POS2 as i16;
        &h_inv[L_SUBFR - iy..]
    };
    *index = add(shl(i0, 6), i1);
    for i in 0..L_SUBFR {
        y[i] = shr_r(add(p0[i], p1[i]), 3);
    }
}

/// The C's `rrixix` and `rrixiy` chains of the 2-track search; `PLAIN`: their adds without
/// saturation (exact when the caller's bound holds).
#[inline(always)]
fn cor_h2<const PLAIN: bool>(
    hp: &[i16],
    rrixix: &mut [[i16; NB_POS2]; 2],
    rrixiy: &mut [i16; MSIZE2],
) {
    let mac = |c: i32, a: i16, b: i16| {
        if PLAIN {
            c.wrapping_add((a as i32 * b as i32).wrapping_mul(2))
        } else {
            l_mac(c, a, b)
        }
    };
    let mut cor = 0x0001_0000i32;
    for i in 0..NB_POS2 {
        cor = mac(cor, hp[2 * i], hp[2 * i]);
        rrixix[1][NB_POS2 - 1 - i] = extract_h(cor);
        cor = mac(cor, hp[2 * i + 1], hp[2 * i + 1]);
        rrixix[0][NB_POS2 - 1 - i] = extract_h(cor);
    }
    for r in rrixix.iter_mut() {
        for v in r.iter_mut() {
            *v >>= 1;
        }
    }
    for k in 0..NB_POS2 {
        let (mut p1, mut p0) = (
            (MSIZE2 - 1 - k * NB_POS2) as isize,
            (MSIZE2 - 2 - k) as isize,
        );
        let mut cor = 0x8000i32;
        let (mut a, mut b) = (0usize, 1 + 2 * k);
        for _ in k + 1..NB_POS2 {
            cor = mac(cor, hp[a], hp[b]);
            a += 1;
            b += 1;
            rrixiy[p1 as usize] = extract_h(cor);
            cor = mac(cor, hp[a], hp[b]);
            a += 1;
            b += 1;
            rrixiy[p0 as usize] = extract_h(cor);
            p1 -= NB_POS2 as isize + 1;
            p0 -= NB_POS2 as isize + 1;
        }
        cor = mac(cor, hp[a], hp[b]);
        rrixiy[p1 as usize] = extract_h(cor);
    }
}

/// The 2-track pair search as the C's chain: (ix, iy).
#[inline(always)]
fn search_2t_chain(
    dn: &[i16; L_SUBFR],
    rrixix: &[[i16; NB_POS2]; 2],
    rrixiy: &[i16; MSIZE2],
) -> (usize, usize) {
    let (mut psk, mut alpk) = (-1i16, 1i16);
    let (mut ix, mut iy) = (0usize, 1usize);
    let mut p2 = 0;
    for (n0, i0) in (0..L_SUBFR).step_by(2).enumerate() {
        let ps1 = dn[i0];
        let alp1 = rrixix[0][n0];
        let mut pos = None;
        for (n1, i1) in (1..L_SUBFR).step_by(2).enumerate() {
            let ps2 = add(ps1, dn[i1]);
            let alp2 = add(alp1, add(rrixix[1][n1], rrixiy[p2]));
            p2 += 1;
            let sq = mult(ps2, ps2);
            let s = l_msu(l_mult(alpk, sq), psk, alp2);
            if s > 0 {
                psk = sq;
                alpk = alp2;
                pos = Some(i1);
            }
        }
        if let Some(p) = pos {
            ix = i0;
            iy = p;
        }
    }
    (ix, iy)
}

/// The 2-track pair search in lanes: as `search_ixiy` (the 4-track search), with every
/// `alp2 > 0` the C's chain keeps the first maximum of `sq/alp2` in (i0, i1) order; rows of 32
/// i1 lanes, each lane's first best over the rows, the lanes reduced in order. `None` if some
/// `alp2 <= 0` (then the chain decides).
#[inline(always)]
fn search_2t(
    dn: &[i16; L_SUBFR],
    rrixix: &[[i16; NB_POS2]; 2],
    rrixiy: &[i16; MSIZE2],
) -> Option<(usize, usize)> {
    let dn1: [i16; NB_POS2] = std::array::from_fn(|n| dn[2 * n + 1]);
    let (mut bsq, mut bal, mut bx) = ([-1i32; NB_POS2], [1i32; NB_POS2], [0usize; NB_POS2]);
    for n0 in 0..NB_POS2 {
        let ps1 = dn[2 * n0];
        let alp1 = rrixix[0][n0];
        let row: &[i16; NB_POS2] = rrixiy[NB_POS2 * n0..NB_POS2 * (n0 + 1)].try_into().unwrap();
        let (mut sq, mut al) = ([0i32; NB_POS2], [0i32; NB_POS2]);
        let mut min_al = i32::MAX;
        for n in 0..NB_POS2 {
            let ps2 = ps1.saturating_add(dn1[n]);
            al[n] = alp1.saturating_add(rrixix[1][n].saturating_add(row[n])) as i32;
            sq[n] = mult(ps2, ps2) as i32;
            min_al = min_al.min(al[n]);
        }
        if min_al <= 0 {
            return None;
        }
        for n in 0..NB_POS2 {
            // 15-bit, non-negative: the cross products fit an i32
            if sq[n] * bal[n] > bsq[n] * al[n] {
                bsq[n] = sq[n];
                bal[n] = al[n];
                bx[n] = n0;
            }
        }
    }
    // every lane took row 0's candidate at least (sq >= 0 > -alp): reduce in order
    let (mut sqk, mut alpk, mut kx, mut ky) = (bsq[0], bal[0], bx[0], 0);
    for n in 1..NB_POS2 {
        let (l, r) = (bsq[n] * alpk, sqk * bal[n]);
        if l > r || l == r && bx[n] < kx {
            (sqk, alpk, kx, ky) = (bsq[n], bal[n], bx[n], n);
        }
    }
    Some((2 * kx, 2 * ky + 1))
}

#[cfg(all(test, oracle))]
#[allow(unsafe_code)]
pub(crate) mod vs_ref {
    use super::*;
    use crate::encoder::filters::vs_ref::signal;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_cor_h_x(h: *mut i16, x: *mut i16, dn: *mut i16);
        fn refwb_ACELP_2t64_fx(
            dn: *mut i16,
            cn: *mut i16,
            h: *mut i16,
            code: *mut i16,
            y: *mut i16,
            index: *mut i16,
        );
    }

    /// A codebook-search input: impulse response (Q12, decaying, or wild), target, residual.
    pub(crate) fn search_input(
        r: &mut Rng,
        n: u64,
    ) -> ([i16; L_SUBFR], [i16; L_SUBFR], [i16; L_SUBFR]) {
        let mut h = [0i16; L_SUBFR];
        let wild = n.is_multiple_of(5);
        for (i, v) in h.iter_mut().enumerate() {
            *v = if wild {
                r.i16()
            } else {
                r.i16() >> (3 + i.min(48) / 4)
            };
        }
        if !wild {
            h[0] = 4096;
        }
        let x: [i16; L_SUBFR] = signal(r, n, L_SUBFR).try_into().unwrap();
        let mut cn: [i16; L_SUBFR] = signal(r, n / 3, L_SUBFR).try_into().unwrap();
        if n.is_multiple_of(7) {
            cn = [0; L_SUBFR];
        }
        (h, x, cn)
    }

    #[test]
    fn cor_h_x_and_acelp_2t64_match_reference() {
        let mut r = Rng::new(81);
        for n in 0..20_000u64 {
            let (mut h, mut x, mut cn) = search_input(&mut r, n);
            let (mut dn, mut dnc) = ([0i16; L_SUBFR], [0i16; L_SUBFR]);
            cor_h_x(&h, &x, &mut dn);
            let (mut code, mut y, mut ix) = ([0i16; L_SUBFR], [0i16; L_SUBFR], 0i16);
            let (mut codec, mut yc, mut ixc) = ([0i16; L_SUBFR], [0i16; L_SUBFR], 0i16);
            with_ref(|| {
                // SAFETY: 64-sample vectors, one index.
                unsafe {
                    refwb_cor_h_x(h.as_mut_ptr(), x.as_mut_ptr(), dnc.as_mut_ptr());
                    assert_eq!(dn, dnc, "cor_h_x {n}");
                    acelp_2t64_fx(&mut dn, &cn, &h, &mut code, &mut y, &mut ix);
                    refwb_ACELP_2t64_fx(
                        dnc.as_mut_ptr(),
                        cn.as_mut_ptr(),
                        h.as_mut_ptr(),
                        codec.as_mut_ptr(),
                        yc.as_mut_ptr(),
                        &mut ixc,
                    );
                }
            });
            assert_eq!(
                (dn, code, y, ix),
                (dnc, codec, yc, ixc),
                "ACELP_2t64_fx {n}"
            );
        }
    }
}
