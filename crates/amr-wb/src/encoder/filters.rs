//! The encoder's small filters (vo-amrwbenc `scale.c`, `preemph.c`, `deemph.c`, `residu.c`,
//! `syn_filt.c`, `convolve.c`, `updt_tar.c`, `lp_dec2.c`, `hp400.c`), with the reference's
//! saturation where vo's plain C arithmetic wraps.

use super::ops::*;
use crate::cnst::{L_SUBFR16K, M, M16K};
use crate::enc_tables::LP_DEC2_H_FIR;

/// `Scale_sig`: `x = round(x << exp)`, saturating.
pub(crate) fn scale_sig(x: &mut [i16], exp: i16) {
    for v in x {
        *v = round(l_shl(l_deposit_h(*v), exp));
    }
}

/// `Preemph`: `x[i] -= mu·x[i-1]` (`mem` = x[-1]).
pub(crate) fn preemph(x: &mut [i16], mu: i16, mem: &mut i16) {
    let last = x[x.len() - 1];
    for i in (1..x.len()).rev() {
        x[i] = round(l_msu(l_deposit_h(x[i]), x[i - 1], mu));
    }
    x[0] = round(l_msu(l_deposit_h(x[0]), *mem, mu));
    *mem = last;
}

/// `Preemph2`: as [`preemph`], the output doubled.
pub(crate) fn preemph2(x: &mut [i16], mu: i16, mem: &mut i16) {
    let last = x[x.len() - 1];
    for i in (1..x.len()).rev() {
        x[i] = round(l_shl(l_msu(l_deposit_h(x[i]), x[i - 1], mu), 1));
    }
    x[0] = round(l_shl(l_msu(l_deposit_h(x[0]), *mem, mu), 1));
    *mem = last;
}

/// `Deemph2`: `y[i] = x[i]/2 + mu·y[i-1]` (`mem` = y[-1]).
pub(crate) fn deemph2(x: &mut [i16], mu: i16, mem: &mut i16) {
    let mut prev = *mem;
    for v in x.iter_mut() {
        *v = round(l_mac(l_mult(*v, 16384), prev, mu));
        prev = *v;
    }
    *mem = prev;
}

/// `Residu` (order 16): `y[i] = Σ a[j]·x[i-j]`; `x` holds 16 samples of history before
/// position `at`, the outputs are `y[..lg]` for `x[at..at + lg]`.
pub(crate) fn residu(a: &[i16], x: &[i16], at: usize, y: &mut [i16], lg: usize) {
    for i in 0..lg {
        let n = at + i;
        let mut s = l_mult(x[n], a[0]);
        for j in 1..=M {
            s = l_mac(s, a[j], x[n - j]);
        }
        y[i] = round(l_shl(s, 4));
    }
}

/// `Syn_filt` (order `a.len() - 1`, at most 20, `lg` at most 80): `y = x / A(z)`; `mem` the
/// last outputs, updated if `update`. `x` and `y` may be one buffer: pass `x = None`.
pub(crate) fn syn_filt(a: &[i16], x: Option<&[i16]>, y: &mut [i16], mem: &mut [i16], update: bool) {
    let m = a.len() - 1;
    let lg = y.len();
    let mut yy = [0i16; L_SUBFR16K + M16K];
    yy[..m].copy_from_slice(&mem[..m]);
    let s = sub(norm_s(a[0]), 2);
    let a0 = a[0] >> 1;
    for i in 0..lg {
        let xi = x.map_or(y[i], |x| x[i]);
        let mut l = l_mult(xi, a0);
        for j in 1..=m {
            l = l_msu(l, a[j], yy[m + i - j]);
        }
        yy[m + i] = round(l_shl(l, add(3, s)));
        y[i] = yy[m + i];
    }
    if update {
        mem[..m].copy_from_slice(&yy[lg..lg + m]);
    }
}

/// `Convolve`: `y[n] = Σ_{i<=n} x[i]·h[n-i]`, rounded.
pub(crate) fn convolve(x: &[i16], h: &[i16], y: &mut [i16], l: usize) {
    for n in 0..l {
        let mut s = 0i32;
        for i in 0..=n {
            s = l_mac(s, x[i], h[n - i]);
        }
        y[n] = round(s);
    }
}

/// `Updt_tar`: `x2 = x - gain·y` (gain Q14).
pub(crate) fn updt_tar(x: &[i16], x2: &mut [i16], y: &[i16], gain: i16, l: usize) {
    for i in 0..l {
        let t = l_msu(l_mult(x[i], 16384), y[i], gain);
        x2[i] = extract_h(l_shl(t, 1));
    }
}

/// `LP_Decim2`: 5-tap low-pass and decimation by 2, in place (`x[..l/2]` out).
pub(crate) fn lp_decim2(x: &mut [i16], l: usize, mem: &mut [i16; 3]) {
    let mut buf = [0i16; crate::cnst::L_FRAME + 3];
    buf[..3].copy_from_slice(mem);
    buf[3..3 + l].copy_from_slice(&x[..l]);
    mem.copy_from_slice(&x[l - 3..l]);
    for (j, i) in (0..l).step_by(2).enumerate() {
        let mut s = 0i32;
        for k in 0..5 {
            s = l_mac(s, buf[i + k], LP_DEC2_H_FIR[k]);
        }
        x[j] = round(s);
    }
}

/// `HP400_12k8`: 2nd-order 400 Hz high-pass (output /16), in place; the reference's
/// saturating chain.
pub(crate) fn hp400_12k8(signal: &mut [i16], mem: &mut [i16; 6]) {
    use crate::enc_tables::{HP400_A as A, HP400_B as B};
    let [mut y2_hi, mut y2_lo, mut y1_hi, mut y1_lo, mut x0, mut x1] = *mem;
    for s in signal.iter_mut() {
        let x2 = x1;
        x1 = x0;
        x0 = *s;
        let mut l = 16384i32;
        l = l_mac(l, y1_lo, A[1]);
        l = l_mac(l, y2_lo, A[2]);
        l = l_shr(l, 15);
        l = l_mac(l, y1_hi, A[1]);
        l = l_mac(l, y2_hi, A[2]);
        l = l_mac(l, x0, B[0]);
        l = l_mac(l, x1, B[1]);
        l = l_mac(l, x2, B[2]);
        l = l_shl(l, 1);
        y2_hi = y1_hi;
        y2_lo = y1_lo;
        (y1_hi, y1_lo) = l_extract(l);
        *s = round(l);
    }
    *mem = [y2_hi, y2_lo, y1_hi, y1_lo, x0, x1];
}

#[cfg(test)]
#[allow(unsafe_code)]
pub(crate) mod vs_ref {
    use super::*;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_Scale_sig(x: *mut i16, lg: i16, exp: i16);
        fn refwb_Preemph(x: *mut i16, mu: i16, lg: i16, mem: *mut i16);
        fn refwb_Preemph2(x: *mut i16, mu: i16, lg: i16, mem: *mut i16);
        fn refwb_Deemph2(x: *mut i16, mu: i16, lg: i16, mem: *mut i16);
        fn refwb_Residu(a: *mut i16, m: i16, x: *mut i16, y: *mut i16, lg: i16);
        fn refwb_Syn_filt(
            a: *mut i16,
            m: i16,
            x: *mut i16,
            y: *mut i16,
            lg: i16,
            mem: *mut i16,
            update: i16,
        );
        fn refwb_Convolve(x: *mut i16, h: *mut i16, y: *mut i16, l: i16);
        fn refwb_Updt_tar(x: *mut i16, x2: *mut i16, y: *mut i16, gain: i16, l: i16);
        fn refwb_LP_Decim2(x: *mut i16, l: i16, mem: *mut i16);
        fn refwb_HP400_12k8(x: *mut i16, lg: i16, mem: *mut i16);
    }

    /// A signal class: random, loud random, full-scale square, quiet.
    pub(crate) fn signal(r: &mut Rng, class: u64, n: usize) -> Vec<i16> {
        let mut t = 0u64;
        (0..n)
            .map(|_| {
                t += 1;
                match class % 4 {
                    0 => r.i16(),
                    1 => [i16::MAX, i16::MIN][((t / (1 + class % 7)) % 2) as usize],
                    2 => r.i16() >> (r.u64() % 14),
                    _ => (r.u64() % 64) as i16 - 32,
                }
            })
            .collect()
    }

    /// Stable-ish LPC filters in Q12 (a[0] = 4096) and wild ones (any coefficients).
    pub(crate) fn lpc(r: &mut Rng, m: usize) -> Vec<i16> {
        let mut a = vec![4096i16; m + 1];
        let wild = r.u64().is_multiple_of(4);
        for (k, c) in a.iter_mut().enumerate().skip(1) {
            *c = if wild {
                r.i16()
            } else {
                r.i16() >> (2 + k / 3)
            };
        }
        a
    }

    #[test]
    fn small_filters_match_reference() {
        let mut r = Rng::new(21);
        for seq in 0..400u64 {
            let (mut m1, mut c1, mut m2, mut c2, mut m3, mut c3) = (0i16, 0i16, 0, 0, 0, 0);
            let (mut md, mut cd) = ([0i16; 3], [0i16; 3]);
            let (mut mh, mut ch) = ([0i16; 6], [0i16; 6]);
            for call in 0..20 {
                let x = signal(&mut r, seq + call, 256);
                let mu = r.i16();
                let exp = (r.u64() % 25) as i16 - 12;
                let lg = [64usize, 80, 32, 256][(r.u64() % 4) as usize];
                with_ref(|| {
                    // SAFETY: every buffer holds lg (<= 256) samples, memories their sizes.
                    unsafe {
                        let (mut a, mut b) = (x.clone(), x.clone());
                        scale_sig(&mut a[..lg], exp);
                        refwb_Scale_sig(b.as_mut_ptr(), lg as i16, exp);
                        assert_eq!(a, b, "Scale_sig {seq} {call}");
                        let (mut a, mut b) = (x.clone(), x.clone());
                        preemph(&mut a[..lg], mu, &mut m1);
                        refwb_Preemph(b.as_mut_ptr(), mu, lg as i16, &mut c1);
                        assert_eq!((a, m1), (b, c1), "Preemph");
                        let (mut a, mut b) = (x.clone(), x.clone());
                        preemph2(&mut a[..lg], mu, &mut m2);
                        refwb_Preemph2(b.as_mut_ptr(), mu, lg as i16, &mut c2);
                        assert_eq!((a, m2), (b, c2), "Preemph2");
                        let (mut a, mut b) = (x.clone(), x.clone());
                        deemph2(&mut a[..lg], mu, &mut m3);
                        refwb_Deemph2(b.as_mut_ptr(), mu, lg as i16, &mut c3);
                        assert_eq!((a, m3), (b, c3), "Deemph2");
                        let (mut a, mut b) = (x.clone(), x.clone());
                        let l = lg & !1;
                        lp_decim2(&mut a, l, &mut md);
                        refwb_LP_Decim2(b.as_mut_ptr(), l as i16, cd.as_mut_ptr());
                        assert_eq!((a, md), (b, cd), "LP_Decim2");
                        let (mut a, mut b) = (x.clone(), x.clone());
                        hp400_12k8(&mut a[..lg], &mut mh);
                        refwb_HP400_12k8(b.as_mut_ptr(), lg as i16, ch.as_mut_ptr());
                        assert_eq!((a, mh), (b, ch), "HP400");
                        // Residu: 16 history samples, then lg
                        let mut q = lpc(&mut r, M);
                        let mut y = [0i16; 256];
                        let mut yc = [0i16; 256];
                        let n = lg.min(240);
                        let mut xc = x.clone();
                        residu(&q, &x, 16, &mut y, n);
                        refwb_Residu(
                            q.as_mut_ptr(),
                            16,
                            xc.as_mut_ptr().add(16),
                            yc.as_mut_ptr(),
                            n as i16,
                        );
                        assert_eq!(y, yc, "Residu");
                        // Syn_filt, orders 16 and 20, both updates, in place and not
                        let m = if r.u64().is_multiple_of(2) { M } else { M16K };
                        let n = lg.min(80);
                        let mut q = lpc(&mut r, m);
                        let mut mem: Vec<i16> = signal(&mut r, seq, m);
                        let mut memc = mem.clone();
                        let upd = r.u64().is_multiple_of(2);
                        let (mut y, mut yc) = ([0i16; 80], [0i16; 80]);
                        if r.u64().is_multiple_of(2) {
                            y[..n].copy_from_slice(&x[..n]);
                            yc[..n].copy_from_slice(&x[..n]);
                            syn_filt(&q, None, &mut y[..n], &mut mem, upd);
                            let p = yc.as_mut_ptr();
                            refwb_Syn_filt(
                                q.as_mut_ptr(),
                                m as i16,
                                p,
                                p,
                                n as i16,
                                memc.as_mut_ptr(),
                                upd as i16,
                            );
                        } else {
                            syn_filt(&q, Some(&x[..n]), &mut y[..n], &mut mem, upd);
                            refwb_Syn_filt(
                                q.as_mut_ptr(),
                                m as i16,
                                xc.as_mut_ptr(),
                                yc.as_mut_ptr(),
                                n as i16,
                                memc.as_mut_ptr(),
                                upd as i16,
                            );
                        }
                        assert_eq!((y, mem), (yc, memc), "Syn_filt");
                        // Convolve, Updt_tar on 64
                        let mut h = signal(&mut r, seq + 1, 64);
                        let (mut y, mut yc) = ([0i16; 64], [0i16; 64]);
                        convolve(&x[..64], &h, &mut y, 64);
                        refwb_Convolve(xc.as_mut_ptr(), h.as_mut_ptr(), yc.as_mut_ptr(), 64);
                        assert_eq!(y, yc, "Convolve");
                        let g = r.i16();
                        updt_tar(&x[..64], &mut y, &h, g, 64);
                        refwb_Updt_tar(xc.as_mut_ptr(), yc.as_mut_ptr(), h.as_mut_ptr(), g, 64);
                        assert_eq!(y, yc, "Updt_tar");
                    }
                });
            }
        }
    }
}
