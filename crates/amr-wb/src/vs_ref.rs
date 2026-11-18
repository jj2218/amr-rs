//! Unit differential tests against the TS 26.173 reference C (`ref26173-sys`): the operators,
//! the kernels with exact fast paths (each with an input class that drives its guard false,
//! so a guard-always-true mutant fails), and the stateful filters over call sequences
//! (outputs and memories compared after every call).
#![allow(unsafe_code)]

use crate::basic_op::*;
use crate::excitation::{agc2_amr_wb, pred_lt4};
use crate::lpc::{isf_isp, isp_az};
use crate::synthesis::*;
use crate::testutil::{Rng, with_ref};

unsafe extern "C" {
    fn refwb_shl(a: i16, b: i16) -> i16;
    fn refwb_L_shl(a: i32, b: i16) -> i32;
    fn refwb_L_shr(a: i32, b: i16) -> i32;
    fn refwb_L_shr_r(a: i32, b: i16) -> i32;
    fn refwb_round(a: i32) -> i16;
    fn refwb_norm_s(a: i16) -> i16;
    fn refwb_norm_l(a: i32) -> i16;
    fn refwb_Pow2(exponant: i16, fraction: i16) -> i32;
    fn refwb_Isf_isp(isf: *mut i16, isp: *mut i16, m: i16);
    fn refwb_Pred_lt4(exc: *mut i16, t0: i16, frac: i16, l_subfr: i16);
    fn refwb_Syn_filt_32(
        a: *mut i16,
        m: i16,
        exc: *mut i16,
        q: i16,
        hi: *mut i16,
        lo: *mut i16,
        lg: i16,
    );
    fn refwb_Syn_filt(
        a: *mut i16,
        m: i16,
        x: *mut i16,
        y: *mut i16,
        lg: i16,
        mem: *mut i16,
        update: i16,
    );
    fn refwb_Deemph_32(hi: *mut i16, lo: *mut i16, y: *mut i16, mu: i16, l: i16, mem: *mut i16);
    fn refwb_HP50_12k8(signal: *mut i16, lg: i16, mem: *mut i16);
    fn refwb_HP400_12k8(signal: *mut i16, lg: i16, mem: *mut i16);
    fn refwb_Oversamp_16k(sig12k8: *mut i16, lg: i16, sig16k: *mut i16, mem: *mut i16);
    fn refwb_Filt_6k_7k(signal: *mut i16, lg: i16, mem: *mut i16);
    fn refwb_Filt_7k(signal: *mut i16, lg: i16, mem: *mut i16);
    fn refwb_Dot_product12(x: *mut i16, y: *mut i16, lg: i16, exp: *mut i16) -> i32;
    fn refwb_Isp_Az(isp: *mut i16, a: *mut i16, m: i16, adaptive_scaling: i16);
    fn refwb_agc2(sig_in: *mut i16, sig_out: *mut i16, l_trm: i16);
}

const N: usize = 50_000;

#[test]
fn operators_match_reference() {
    let mut r = Rng::new(1);
    with_ref(|| {
        for _ in 0..N {
            let (a, x, n) = (r.i16(), r.i32(), r.shift());
            // SAFETY: plain value functions.
            unsafe {
                assert_eq!(shl_int16(a, n), refwb_shl(a, n), "shl({a}, {n})");
                assert_eq!(shl_int32(x, n), refwb_L_shl(x, n), "L_shl({x}, {n})");
                assert_eq!(shr_int32(x, n), refwb_L_shr(x, n), "L_shr({x}, {n})");
                assert_eq!(l_shr_r(x, n), refwb_L_shr_r(x, n), "L_shr_r({x}, {n})");
                assert_eq!(amr_wb_round(x), refwb_round(x), "round({x})");
                assert_eq!(norm_s(a), refwb_norm_s(a), "norm_s({a})");
                assert_eq!(normalize_amr_wb(x), refwb_norm_l(x), "norm_l({x})");
                let (e, f) = ((r.u64() % 64) as i16 - 16, a & 0x7fff);
                assert_eq!(power_of_2(e, f), refwb_Pow2(e, f), "Pow2({e}, {f})");
            }
        }
    });
}

/// Energies (the codec's only use) and dot products, quiet to saturating (`loud`, -32768
/// included), 64 and 80 samples.
#[test]
fn dot_product12_matches_reference() {
    let mut r = Rng::new(6);
    for k in 0..N / 10 {
        let n = if k % 2 == 0 { 64 } else { 80 };
        let shift = (r.u64() % 16) as u32;
        let loud = k % 5 == 0;
        let mut x: Vec<i16> = (0..n).map(|_| sample(&mut r, loud) >> shift).collect();
        let mut y: Vec<i16> = (0..n).map(|_| sample(&mut r, loud)).collect();
        with_ref(|| {
            for same in [true, false] {
                let (mut e, mut ce) = (0i16, 0i16);
                let got = if same {
                    dot_product12(&x, &x, &mut e)
                } else {
                    dot_product12(&x, &y, &mut e)
                };
                let yp = if same { x.as_mut_ptr() } else { y.as_mut_ptr() };
                // SAFETY: both buffers hold n samples.
                let want = unsafe { refwb_Dot_product12(x.as_mut_ptr(), yp, n as i16, &mut ce) };
                assert_eq!((got, e), (want, ce), "Dot_product12 case {k} same {same}");
            }
        });
    }
}

/// ISPs (cosines, Q15) of order 16 and 20, both scalings: realistic (decreasing), clustered at
/// ±1 (the polynomial chains saturate) and arbitrary.
#[test]
fn isp_az_matches_reference() {
    let mut r = Rng::new(8);
    for k in 0..N / 10 {
        let m = if k % 2 == 0 { 16 } else { 20 };
        let mut isp: Vec<i16> = match k % 3 {
            0 => {
                let mut v: Vec<i16> = (0..m).map(|_| r.i16()).collect();
                v.sort_unstable_by(|a, b| b.cmp(a));
                v
            }
            1 => (0..m)
                .map(|_| match r.u64() % 4 {
                    0 => i16::MIN,
                    1 => i16::MAX,
                    2 => i16::MIN + (r.u64() % 2048) as i16,
                    _ => i16::MAX - (r.u64() % 2048) as i16,
                })
                .collect(),
            _ => (0..m).map(|_| r.i16()).collect(),
        };
        let adaptive = r.u64().is_multiple_of(2);
        let mut a = [0i16; 21];
        let mut c = [0i16; 21];
        isp_az(&isp, &mut a[..=m], adaptive);
        // SAFETY: isp holds m values, c m + 1.
        with_ref(|| unsafe {
            refwb_Isp_Az(isp.as_mut_ptr(), c.as_mut_ptr(), m as i16, adaptive as i16)
        });
        assert_eq!(a, c, "Isp_Az case {k} m {m} adaptive {adaptive}");
    }
}

/// The 6.60/8.85 kbit/s gain control: quiet to full-scale (saturating energies and outputs)
/// signals, silent ones.
#[test]
fn agc2_matches_reference() {
    let mut r = Rng::new(10);
    for k in 0..N / 10 {
        let (si, so) = ((r.u64() % 16) as u32, (r.u64() % 16) as u32);
        let (li, lo) = (k % 3 == 0, k % 5 == 0);
        let mut x: Vec<i16> = (0..64).map(|_| sample(&mut r, li) >> si).collect();
        let mut y: Vec<i16> = (0..64).map(|_| sample(&mut r, lo) >> so).collect();
        if k % 7 == 0 {
            x.fill(0);
        }
        let mut c = y.clone();
        agc2_amr_wb(&x, &mut y);
        // SAFETY: both hold 64 samples.
        with_ref(|| unsafe { refwb_agc2(x.as_mut_ptr(), c.as_mut_ptr(), 64) });
        assert_eq!(y, c, "agc2 case {k}");
    }
}

/// A signal sample: mostly realistic, sometimes edge values; `loud` pins most at full scale.
fn sample(r: &mut Rng, loud: bool) -> i16 {
    if loud && !r.u64().is_multiple_of(4) {
        if r.u64() & 1 == 0 { i16::MAX } else { i16::MIN }
    } else {
        r.i16()
    }
}

#[test]
fn pred_lt4_matches_reference() {
    let mut r = Rng::new(2);
    const PAST: usize = 231 + 17;
    for k in 0..N / 10 {
        let loud = k % 2 == 0; // full-scale past excitation: row 1's guard fails
        let mut buf = [0i16; PAST + 65];
        for v in &mut buf[..PAST] {
            *v = sample(&mut r, loud);
        }
        let t0 = 34 + (r.u64() % 198) as i16;
        let frac = (r.u64() % 4) as i16;
        if k % 4 == 1 {
            // a full-scale window with the row's signs at a random output: the chain saturates
            let (mut f, mut st) = (-frac, PAST as isize - t0 as isize);
            if f < 0 {
                f += 4;
                st -= 1;
            }
            let start = (st - 15) as usize;
            let row = &crate::tables::INTER4_2[(3 - f as usize) * 32..][..32];
            let n = (r.u64() % 16) as usize;
            for (j, &c) in row.iter().enumerate() {
                if start + n + j < PAST {
                    buf[start + n + j] = if c < 0 { -32768 } else { 32767 };
                }
            }
        }
        let mut c = buf;
        pred_lt4(&mut buf, PAST, t0, frac, 65);
        // SAFETY: the reference reads exc[-t0-16..] and writes exc[0..65], all inside `c`.
        with_ref(|| unsafe { refwb_Pred_lt4(c.as_mut_ptr().add(PAST), t0, frac, 65) });
        assert_eq!(buf, c, "Pred_lt4 t0 {t0} frac {frac} loud {loud}");
    }
}

/// LP coefficients: a[0] = 4096 >> q (adaptive scaling), then either realistic or large ones
/// (Σ|a| past the fast paths' bounds).
fn lpc(r: &mut Rng, m: usize, large: bool) -> Vec<i16> {
    let mut a = vec![4096 >> (r.u64() % 4); m + 1];
    for v in &mut a[1..] {
        *v = if large {
            r.i16()
        } else {
            (r.u64() % 8192) as i16 - 4096
        };
    }
    a
}

#[test]
fn syn_filt_32_matches_reference() {
    let mut r = Rng::new(3);
    const M: usize = 16;
    for k in 0..N / 10 {
        let large = k % 2 == 0;
        let mut a = lpc(&mut r, M, large);
        let q = (r.u64() % 9) as i16;
        let mut exc: Vec<i16> = (0..64).map(|_| sample(&mut r, large)).collect();
        let mut hi = [0i16; M + 64];
        let mut lo = [0i16; M + 64];
        for j in 0..M {
            hi[j] = r.i16();
            lo[j] = if k % 7 == 0 {
                r.i16()
            } else {
                (r.u64() % 4096) as i16
            };
        }
        let (mut ch, mut cl) = (hi, lo);
        syn_filt_32(&a, M, &exc, q, &mut hi, &mut lo, 64);
        // SAFETY: the reference reads 16 past samples before the output pointers.
        with_ref(|| unsafe {
            refwb_Syn_filt_32(
                a.as_mut_ptr(),
                M as i16,
                exc.as_mut_ptr(),
                q,
                ch.as_mut_ptr().add(M),
                cl.as_mut_ptr().add(M),
                64,
            )
        });
        assert_eq!((hi, lo), (ch, cl), "Syn_filt_32 case {k} large {large}");
    }
    // the low halves' chain saturating alone (quiet high halves): the fast path's 2·Lmax·A bound.
    // A stable filter, (1 + 0.17/z)^16 (Σ|a[1..]| = 11.3 in Q12), so the outputs stay small.
    let mut p = [0f64; M + 1];
    p[0] = 1.0;
    for _ in 0..M {
        for k in (1..=M).rev() {
            p[k] += 0.17 * p[k - 1];
        }
    }
    let mut a: [i16; M + 1] = std::array::from_fn(|k| (p[k] * 4096.0).round() as i16);
    let mut exc = [0i16; 64];
    let (mut hi, mut lo) = ([0i16; M + 64], [0i16; M + 64]);
    lo[..M].fill(-32768);
    let (mut ch, mut cl) = (hi, lo);
    syn_filt_32(&a, M, &exc, 0, &mut hi, &mut lo, 64);
    // SAFETY: as above.
    with_ref(|| unsafe {
        refwb_Syn_filt_32(
            a.as_mut_ptr(),
            M as i16,
            exc.as_mut_ptr(),
            0,
            ch.as_mut_ptr().add(M),
            cl.as_mut_ptr().add(M),
            64,
        )
    });
    assert_eq!((hi, lo), (ch, cl), "Syn_filt_32, saturating low halves");
}

#[test]
fn wb_syn_filt_matches_reference() {
    let mut r = Rng::new(4);
    for k in 0..N / 10 {
        let m = if k % 3 == 0 { 20 } else { 16 };
        let large = k % 2 == 0;
        let mut a = lpc(&mut r, m, large);
        let mut x: Vec<i16> = (0..80).map(|_| sample(&mut r, large)).collect();
        let mut mem: Vec<i16> = (0..m).map(|_| r.i16()).collect();
        let (mut cx, mut cy, mut cmem) = (x.clone(), vec![0i16; 80], mem.clone());
        wb_syn_filt(&a, &mut x, &mut mem);
        // SAFETY: buffers of the lengths passed.
        with_ref(|| unsafe {
            refwb_Syn_filt(
                a.as_mut_ptr(),
                m as i16,
                cx.as_mut_ptr(),
                cy.as_mut_ptr(),
                80,
                cmem.as_mut_ptr(),
                1,
            )
        });
        assert_eq!(
            (x, mem),
            (cy, cmem),
            "Syn_filt m {m} case {k} large {large}"
        );
    }
}

/// Stateful filters: 40 calls per sequence on loud and realistic subframes, outputs and
/// memories compared after each call.
#[test]
fn filters_match_reference_over_sequences() {
    let mut r = Rng::new(5);
    for seq in 0..200 {
        let loud = seq % 2 == 0;
        // full-scale square waves: a high-pass's undershoot plus the next edge saturates HP50
        // (random loud samples never reach it)
        let (square, hold, mut t) = (seq % 4 == 1, 32 + r.u64() % 224, 0u64);
        let (mut m50, mut c50) = ([0i16; 6], [0i16; 6]);
        let (mut m400, mut c400) = ([0i16; 6], [0i16; 6]);
        let (mut mos, mut cos) = ([0i16; 24], [0i16; 24]);
        let (mut mbp, mut cbp) = ([0i16; 30], [0i16; 30]);
        let (mut mlp, mut clp) = ([0i16; 30], [0i16; 30]);
        let (mut mde, mut cde) = (0i16, 0i16);
        for call in 0..40 {
            let mut next = |r: &mut Rng| {
                t += 1;
                match square {
                    true if (t / hold) % 2 == 0 => i16::MAX,
                    true => i16::MIN,
                    false => sample(r, loud),
                }
            };
            let s64: [i16; 64] = std::array::from_fn(|_| next(&mut r));
            let s80: [i16; 80] = std::array::from_fn(|_| next(&mut r));
            let lo: [i16; 64] = std::array::from_fn(|_| (r.u64() % 4096) as i16);
            let ctx = format!("sequence {seq} call {call}");
            with_ref(|| {
                // SAFETY: every buffer has the length passed; memories as the reference sizes them.
                unsafe {
                    let (mut a, mut b) = (s64, s64);
                    highpass_50hz_at_12k8(&mut a, &mut m50);
                    refwb_HP50_12k8(b.as_mut_ptr(), 64, c50.as_mut_ptr());
                    assert_eq!((a, m50), (b, c50), "HP50 {ctx}");
                    let (mut a, mut b) = (s64, s64);
                    highpass_400hz_at_12k8(&mut a, &mut m400);
                    refwb_HP400_12k8(b.as_mut_ptr(), 64, c400.as_mut_ptr());
                    assert_eq!((a, m400), (b, c400), "HP400 {ctx}");
                    let (mut a, mut b) = ([0i16; 80], [0i16; 80]);
                    let mut s = s64;
                    oversamp_12k8_to_16k(&s64, &mut a, &mut mos);
                    refwb_Oversamp_16k(s.as_mut_ptr(), 64, b.as_mut_ptr(), cos.as_mut_ptr());
                    assert_eq!((a, mos), (b, cos), "Oversamp_16k {ctx}");
                    let (mut a, mut b) = (s80, s80);
                    band_pass_6k_7k(&mut a, &mut mbp);
                    refwb_Filt_6k_7k(b.as_mut_ptr(), 80, cbp.as_mut_ptr());
                    assert_eq!((a, mbp), (b, cbp), "Filt_6k_7k {ctx}");
                    let (mut a, mut b) = (s80, s80);
                    low_pass_filt_7k(&mut a, &mut mlp);
                    refwb_Filt_7k(b.as_mut_ptr(), 80, clp.as_mut_ptr());
                    assert_eq!((a, mlp), (b, clp), "Filt_7k {ctx}");
                    let (mut a, mut b) = ([0i16; 64], [0i16; 64]);
                    let (mut h, mut l) = (s64, lo);
                    deemphasis_32(&s64, &lo, &mut a, 22282, &mut mde);
                    refwb_Deemph_32(
                        h.as_mut_ptr(),
                        l.as_mut_ptr(),
                        b.as_mut_ptr(),
                        22282,
                        64,
                        &mut cde,
                    );
                    assert_eq!((a, mde), (b, cde), "Deemph_32 {ctx}");
                }
            });
        }
    }
}

/// Every ISF in [0, 16384) (the last one, doubled, in [0, 8192)). At 16384 the C reads
/// `table[129]`, past its table (undefined, so not called here): there the value is
/// `table[128]` (`lpc::isf_to_cos`), checked directly.
#[test]
fn isf_isp_matches_reference_on_its_whole_domain() {
    let mut isp = [0i16; 16];
    let mut top = [16384; 16];
    top[15] = 8192; // doubled
    isf_isp(&top, &mut isp);
    assert_eq!(isp, [crate::tables::TABLE[128]; 16]);
    for base in (0..16384i32).step_by(15) {
        let mut isf: [i16; 16] = std::array::from_fn(|i| (base + i as i32).min(16383) as i16);
        isf[15] = (base / 2) as i16;
        let mut isp = [0i16; 16];
        isf_isp(&isf, &mut isp);
        let mut c = [0i16; 16];
        // SAFETY: 16 ISFs in, 16 ISPs out (the table read at 16384 is the C's own, see above).
        with_ref(|| unsafe { refwb_Isf_isp(isf.as_mut_ptr(), c.as_mut_ptr(), 16) });
        assert_eq!(isp, c, "Isf_isp from {base}");
    }
}
