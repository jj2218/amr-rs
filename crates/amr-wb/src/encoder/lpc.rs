//! LP analysis and ISF quantisation (vo-amrwbenc `autocorr.c`, `lag_wind.c`, `levinson.c`,
//! `az_isp.c`, `isp_isf.c`, `qpisf_2s.c`, `qisf_ns.c`), with the reference's saturation.
//! ISF → ISP, A(z) and interpolation are the decoder's (`crate::lpc`).

use super::ops::*;
use crate::cnst::M;
use crate::enc_tables::*;

/// LPC analysis window (`L_WINDOW`).
pub(crate) const L_WINDOW: usize = 384;
const NC: usize = M / 2;
const GRID_POINTS: usize = 100;

/// `Autocorr` (order 16): windowed autocorrelations `r` in double precision (r_h, r_l).
pub(crate) fn autocorr(x: &[i16], r_h: &mut [i16; M + 1], r_l: &mut [i16; M + 1]) {
    let mut y = [0i16; L_WINDOW];
    for i in 0..L_WINDOW {
        y[i] = mult_r(x[i], VO_WINDOW[i]);
    }
    let mut sum = l_deposit_h(16);
    for &v in &y {
        sum = l_add(sum, l_shr(l_mult(v, v), 8));
    }
    let shift = sub(4, norm_l(sum) >> 1).max(0);
    for v in &mut y {
        *v = shr_r(*v, shift);
    }
    let mut sum = 1i32;
    for &v in &y {
        sum = l_mac(sum, v, v);
    }
    let norm = norm_l(sum);
    (r_h[0], r_l[0]) = l_extract(l_shl(sum, norm));
    for i in 1..=M {
        let mut s = 0i32;
        for j in 0..L_WINDOW - i {
            s = l_mac(s, y[j], y[j + i]);
        }
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
pub(crate) fn az_isp(a: &[i16], isp: &mut [i16; M], old_isp: &[i16; M]) {
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
    let mut xlow = VOGRID[0];
    let mut ylow = cheb(xlow, ip);
    let mut j = 0;
    while nf < M - 1 && j < GRID_POINTS {
        j += 1;
        let mut xhigh = xlow;
        let mut yhigh = ylow;
        xlow = VOGRID[j];
        ylow = cheb(xlow, ip);
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
fn chebps2(x: i16, f: &[i16], n: usize) -> i16 {
    let (mut b2_h, mut b2_l) = l_extract(l_mult(f[0], 4096));
    let t0 = l_mac(l_shl(mpy_32_16(b2_h, b2_l, x), 1), f[1], 4096);
    let (mut b1_h, mut b1_l) = l_extract(t0);
    for &fi in &f[2..n] {
        let mut t0 = mpy_32_16(b1_h, b1_l, x);
        t0 = l_mac(t0, b2_h, -16384);
        t0 = l_mac(t0, fi, 2048);
        t0 = l_shl(t0, 1);
        t0 = l_msu(t0, b2_l, 1);
        let (b0_h, b0_l) = l_extract(t0);
        (b2_h, b2_l) = (b1_h, b1_l);
        (b1_h, b1_l) = (b0_h, b0_l);
    }
    let mut t0 = mpy_32_16(b1_h, b1_l, x);
    t0 = l_mac(t0, b2_h, -32768);
    t0 = l_msu(t0, b2_l, 1);
    t0 = l_mac(t0, f[n], 2048);
    let cheb = extract_h(l_shl(t0, 6));
    if cheb == -32768 { -32767 } else { cheb }
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

/// `Sub_VQ`: the nearest of `dico`'s `dim`-vectors to `x` (written to `x`); (index, distance).
pub(crate) fn sub_vq(x: &mut [i16], dico: &[i16], dim: usize, dico_size: usize) -> (i16, i32) {
    let mut dist_min = i32::MAX;
    let mut index = 0;
    for i in 0..dico_size {
        let mut dist = 0i32;
        for j in 0..dim {
            let t = sub(x[j], dico[i * dim + j]);
            dist = l_mac(dist, t, t);
        }
        if dist < dist_min {
            dist_min = dist;
            index = i;
        }
    }
    x[..dim].copy_from_slice(&dico[index * dim..index * dim + dim]);
    (index as i16, dist_min)
}

/// `VQ_stage1`: the `surv` nearest entries of `dico`, best first.
fn vq_stage1(
    x: &[i16],
    dico: &[i16],
    dim: usize,
    dico_size: usize,
    index: &mut [i16],
    surv: usize,
) {
    let mut dist_min = [i32::MAX; N_SURV_MAX];
    for (i, v) in index[..surv].iter_mut().enumerate() {
        *v = i as i16;
    }
    for i in 0..dico_size {
        let mut dist = 0i32;
        for j in 0..dim {
            let t = sub(x[j], dico[i * dim + j]);
            dist = l_mac(dist, t, t);
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

/// The prediction residual `isf1 - mean - MU·past` of both quantisers.
fn residual(isf1: &[i16; M], past_isfq: &[i16; M]) -> [i16; M] {
    std::array::from_fn(|i| sub(sub(isf1[i], MEAN_ISF[i]), mult(MU, past_isfq[i])))
}

/// One split of a two-stage VQ: for each first-stage survivor of `isf` (dimension `dim1`,
/// codebook `dico1`), the second-stage splits `subs` (codebook, dimension); the best total
/// error's indices go to `ind1` and `inds` (unchanged if every total is `MAX_32`).
fn two_stage(
    isf: &[i16],
    dico1: &[i16],
    dim1: usize,
    size1: usize,
    subs: &[(&[i16], usize)],
    nb_surv: usize,
    ind1: &mut i16,
    inds: &mut [i16],
) {
    let mut surv1 = [0i16; N_SURV_MAX];
    vq_stage1(isf, dico1, dim1, size1, &mut surv1, nb_surv);
    let mut distance = i32::MAX;
    for &s in &surv1[..nb_surv] {
        let mut stage2: [i16; 9] = [0; 9];
        for i in 0..dim1 {
            stage2[i] = sub(isf[i], dico1[i + s as usize * dim1]);
        }
        let mut temp = 0i32;
        let mut tmp_ind = [0i16; 3];
        let mut at = 0;
        for (k, &(dico, dim)) in subs.iter().enumerate() {
            let (ix, err) = sub_vq(&mut stage2[at..], dico, dim, dico.len() / dim);
            tmp_ind[k] = ix;
            temp = if k == 0 { err } else { l_add(temp, err) };
            at += dim;
        }
        if temp < distance {
            distance = temp;
            *ind1 = s;
            inds.copy_from_slice(&tmp_ind[..subs.len()]);
        }
    }
}

/// `Qpisf_2s_46b`: 46-bit predictive two-stage ISF quantisation; `isf_q` the quantised ISFs.
pub(crate) fn qpisf_2s_46b(
    isf1: &[i16; M],
    isf_q: &mut [i16; M],
    past_isfq: &mut [i16; M],
    indice: &mut [i16; 7],
    nb_surv: usize,
) {
    let isf = residual(isf1, past_isfq);
    let (i0, rest) = indice.split_at_mut(1);
    let (i1, rest) = rest.split_at_mut(1);
    let subs1: [(&[i16], usize); 3] = [(&DICO21_ISF, 3), (&DICO22_ISF, 3), (&DICO23_ISF, 3)];
    two_stage(
        &isf,
        &DICO1_ISF,
        9,
        256,
        &subs1,
        nb_surv,
        &mut i0[0],
        &mut rest[..3],
    );
    let subs2: [(&[i16], usize); 2] = [(&DICO24_ISF, 3), (&DICO25_ISF, 4)];
    two_stage(
        &isf[9..],
        &DICO2_ISF,
        7,
        256,
        &subs2,
        nb_surv,
        &mut i1[0],
        &mut rest[3..5],
    );
    let old = *isf_q;
    crate::lpc::dpisf_2s_46b(&indice[..], isf_q, past_isfq, &old, &mut [0; 3 * M], false);
}

/// `Qpisf_2s_36b`: 36-bit predictive two-stage ISF quantisation (6.60 kbit/s).
pub(crate) fn qpisf_2s_36b(
    isf1: &[i16; M],
    isf_q: &mut [i16; M],
    past_isfq: &mut [i16; M],
    indice: &mut [i16; 5],
    nb_surv: usize,
) {
    let isf = residual(isf1, past_isfq);
    let (i0, rest) = indice.split_at_mut(1);
    let (i1, rest) = rest.split_at_mut(1);
    let subs1: [(&[i16], usize); 2] = [(&DICO21_ISF_36B, 5), (&DICO22_ISF_36B, 4)];
    two_stage(
        &isf,
        &DICO1_ISF,
        9,
        256,
        &subs1,
        nb_surv,
        &mut i0[0],
        &mut rest[..2],
    );
    let subs2: [(&[i16], usize); 1] = [(&DICO23_ISF_36B, 7)];
    two_stage(
        &isf[9..],
        &DICO2_ISF,
        7,
        256,
        &subs2,
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
    let books: [(&[i16], usize, usize); 5] = [
        (&DICO1_ISF_NOISE, 2, 0),
        (&DICO2_ISF_NOISE, 3, 2),
        (&DICO3_ISF_NOISE, 3, 5),
        (&DICO4_ISF_NOISE, 4, 8),
        (&DICO5_ISF_NOISE, 4, 12),
    ];
    for (k, &(dico, dim, at)) in books.iter().enumerate() {
        indice[k] = sub_vq(&mut isf_q[at..], dico, dim, dico.len() / dim).0;
    }
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
        for _ in 0..20_000 {
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
