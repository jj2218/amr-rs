//! The target–impulse-response correlation (vo-amrwbenc `cor_h_x.c`) and the 2-track
//! algebraic codebook search of 6.60 kbit/s (`c2t64fx.c`), with the reference's saturation.

use super::ops::*;
use crate::basic_op::dot_product12;
use crate::cnst::L_SUBFR;

/// `cor_h_x`: `dn[i] = Σ x[j]·h[j-i]`, scaled together so the codebook search stays in 12 bits.
pub(crate) fn cor_h_x(h: &[i16], x: &[i16], dn: &mut [i16; L_SUBFR]) {
    let mut y32 = [0i32; L_SUBFR];
    let mut l_tot = 1i32;
    for k in 0..4 {
        let mut l_max = 0i32;
        for i in (k..L_SUBFR).step_by(4) {
            let mut l = 1i32;
            for j in i..L_SUBFR {
                l = l_mac(l, x[j], h[j - i]);
            }
            y32[i] = l;
            l_max = l_max.max(l_abs(l));
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

const NB_POS2: usize = 32;
const MSIZE2: usize = 1024;

/// `ACELP_2t64_fx`: two pulses on two interleaved tracks of 32 positions (12-bit index);
/// `code` the codevector (Q9), `y` it filtered (Q9). `dn` gets the signs folded in.
pub(crate) fn acelp_2t64_fx(
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
    let mut cor = 0x0001_0000i32;
    for i in 0..NB_POS2 {
        cor = l_mac(cor, hp[2 * i], hp[2 * i]);
        rrixix[1][NB_POS2 - 1 - i] = extract_h(cor);
        cor = l_mac(cor, hp[2 * i + 1], hp[2 * i + 1]);
        rrixix[0][NB_POS2 - 1 - i] = extract_h(cor);
    }
    for r in rrixix.iter_mut() {
        for v in r.iter_mut() {
            *v >>= 1;
        }
    }
    let mut rrixiy = [0i16; MSIZE2];
    for k in 0..NB_POS2 {
        let (mut p1, mut p0) = (
            (MSIZE2 - 1 - k * NB_POS2) as isize,
            (MSIZE2 - 2 - k) as isize,
        );
        let mut cor = 0x8000i32;
        let (mut a, mut b) = (0usize, 1 + 2 * k);
        for _ in k + 1..NB_POS2 {
            cor = l_mac(cor, hp[a], hp[b]);
            a += 1;
            b += 1;
            rrixiy[p1 as usize] = extract_h(cor);
            cor = l_mac(cor, hp[a], hp[b]);
            a += 1;
            b += 1;
            rrixiy[p0 as usize] = extract_h(cor);
            p1 -= NB_POS2 as isize + 1;
            p0 -= NB_POS2 as isize + 1;
        }
        cor = l_mac(cor, hp[a], hp[b]);
        rrixiy[p1 as usize] = extract_h(cor);
    }
    let mut p0 = 0;
    for i in (0..L_SUBFR).step_by(2) {
        let psign = if sign[i] < 0 { &vec } else { &sign };
        for j in (1..L_SUBFR).step_by(2) {
            rrixiy[p0] = mult(rrixiy[p0], psign[j]);
            p0 += 1;
        }
    }
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

#[cfg(test)]
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
