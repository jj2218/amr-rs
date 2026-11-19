//! Differential tests of `lpc` against the TS 26.073 reference (`ref26073-sys`;
//! helpers in `encoder::vs_ref`).
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::cnst::MP1;
use crate::cnst::{CodecMode, M};
use crate::encoder::layout::cols;
use crate::encoder::lpc::{
    DICO1_SOA, DICO2_HALF_SOA, DICO2_SOA, DICO3_SOA, MR515_3_SOA, MR795_1_SOA, vq_subvec3,
    vq_subvec4,
};
use crate::encoder::lpc::{
    L_WINDOW, LevinsonState, QPlsfState, autocorr, az_lsp, lag_window, levinson, q_plsf_3, q_plsf_5,
};
use crate::encoder::vs_ref::CFlag;
use crate::encoder::vs_ref::refc::*;
use crate::encoder::vs_ref::*;
use crate::tables::{DICO1_LSF_3, DICO2_LSF_3, DICO3_LSF_3, MR515_3_LSF, MR795_1_LSF};
use crate::tables::{WINDOW_160_80, WINDOW_200_40, WINDOW_232_8};
use crate::testutil::with_ref;
use crate::testutil::{N, Rng};

/// Autocorr alone around its overflow test: energies just below, at and above 2^31 (the C
/// stops adding at the first negative sum), -32768 samples, loud and quiet signals.
#[test]
fn autocorr_matches_c() {
    let wins: [&[i16; L_WINDOW]; 3] = [&WINDOW_200_40, &WINDOW_160_80, &WINDOW_232_8];
    let mut g = Rng::new(29);
    let mut n_ov = 0;
    for n in 0..N / 20 {
        let w = wins[n % 3];
        let mut x = [0i16; L_WINDOW];
        match n % 4 {
            0 => signal(&mut g, &mut x),
            1 => x
                .iter_mut()
                .for_each(|v| *v = [-32768, 32767][(g.u64() & 1) as usize]),
            _ => {
                // Windowed energy 2 * sum t^2 within +-2^18 of 2^31: amplitude by
                // bisection, then single samples nudged towards the target.
                let energy = |x: &[i16; L_WINDOW]| -> i64 {
                    (0..L_WINDOW)
                        .map(|i| {
                            let t = ((x[i] as i32 * w[i] as i32 + 0x4000) >> 15) as i64;
                            2 * t * t
                        })
                        .sum()
                };
                let target = (1i64 << 31) + (g.u64() % (1 << 19)) as i64 - (1 << 18);
                let signs: [i16; L_WINDOW] = std::array::from_fn(|_| (g.u64() & 1) as i16 * 2 - 1);
                let (mut lo, mut hi) = (1i16, 32767i16);
                while hi - lo > 1 {
                    let a = lo + (hi - lo) / 2;
                    x = std::array::from_fn(|i| a * signs[i]);
                    if energy(&x) < target { lo = a } else { hi = a }
                }
                x = std::array::from_fn(|i| lo * signs[i]);
                for _ in 0..500 {
                    let i = (g.u64() % L_WINDOW as u64) as usize;
                    let mut y = x;
                    y[i] = y[i].saturating_add((g.u64() % 201) as i16 - 100);
                    if (energy(&y) - target).abs() < (energy(&x) - target).abs() {
                        x = y;
                    }
                }
            }
        }
        let (mut rh, mut rl, mut ch, mut cl) = ([0i16; MP1], [0i16; MP1], [0i16; MP1], [0i16; MP1]);
        let rn = autocorr(&x, M, &mut rh, &mut rl, w);
        let mut o = 0;
        let cn = unsafe {
            c_autocorr(
                x.as_mut_ptr(),
                M as i16,
                ch.as_mut_ptr(),
                cl.as_mut_ptr(),
                w.as_ptr(),
                &mut o,
            )
        };
        n_ov += (rn < 0) as usize;
        assert_eq!((rn, rh, rl), (cn, ch, cl), "Autocorr {x:?}");
    }
    assert!(n_ov > N / 80, "overflow class too rare: {n_ov}");
}

unsafe extern "C" {
    fn ref_Autocorr(x: *mut i16, m: i16, rh: *mut i16, rl: *mut i16, w: *const i16) -> i16;
    fn ref_Lag_window(m: i16, rh: *mut i16, rl: *mut i16);
    fn ref_Levinson(st: *mut i16, rh: *mut i16, rl: *mut i16, a: *mut i16, rc: *mut i16) -> i16;
    fn ref_Az_lsp(a: *mut i16, lsp: *mut i16, old: *mut i16);
    fn ref_Q_plsf_3(
        st: *mut [i16; M],
        mode: i32,
        lsp: *mut i16,
        lsp_q: *mut i16,
        ind: *mut i16,
        pi: *mut i16,
    );
    fn ref_Q_plsf_5(
        st: *mut [i16; M],
        l1: *mut i16,
        l2: *mut i16,
        q1: *mut i16,
        q2: *mut i16,
        ind: *mut i16,
    );
    #[link_name = "Chebps_Wrapper"]
    fn c_chebps(x: i16, f: *mut i16, n: i16, o: *mut CFlag) -> i16;
    #[link_name = "Test_Vq_subvec3"]
    fn c_vq_subvec3(
        r: *mut i16,
        d: *const i16,
        wf: *mut i16,
        n: i16,
        half: i32,
        o: *mut CFlag,
    ) -> i16;
    #[link_name = "Test_Vq_subvec4"]
    fn c_vq_subvec4(r: *mut i16, d: *const i16, wf: *mut i16, n: i16, o: *mut CFlag) -> i16;
}

unsafe fn c_autocorr(
    x: *mut i16,
    m: i16,
    rh: *mut i16,
    rl: *mut i16,
    w: *const i16,
    o: *mut CFlag,
) -> i16 {
    unsafe { call(o, || ref_Autocorr(x, m, rh, rl, w)) }
}

unsafe fn c_lag_window(m: i16, rh: *mut i16, rl: *mut i16, o: *mut CFlag) {
    unsafe { call(o, || ref_Lag_window(m, rh, rl)) }
}

unsafe fn c_levinson<S>(
    st: *mut S,
    rh: *mut i16,
    rl: *mut i16,
    a: *mut i16,
    rc: *mut i16,
    o: *mut CFlag,
) -> i16 {
    unsafe { call(o, || ref_Levinson(st.cast(), rh, rl, a, rc)) }
}

unsafe fn c_az_lsp(a: *mut i16, lsp: *mut i16, old: *mut i16, o: *mut CFlag) {
    with_ref(o, || unsafe { ref_Az_lsp(a, lsp, old) })
}

unsafe fn c_q_plsf_3(
    st: *mut [i16; M],
    mode: i32,
    lsp: *mut i16,
    lsp_q: *mut i16,
    ind: *mut i16,
    pi: *mut i16,
    o: *mut CFlag,
) {
    with_ref(o, || unsafe { ref_Q_plsf_3(st, mode, lsp, lsp_q, ind, pi) })
}

unsafe fn c_q_plsf_5(
    st: *mut [i16; M],
    l1: *mut i16,
    l2: *mut i16,
    q1: *mut i16,
    q2: *mut i16,
    ind: *mut i16,
    o: *mut CFlag,
) {
    with_ref(o, || unsafe { ref_Q_plsf_5(st, l1, l2, q1, q2, ind) })
}

/// Random descending LSP vectors (cosine domain, as Az_lsp produces them).
pub(crate) fn lsp(g: &mut Rng) -> [i16; M] {
    let mut v: [i16; M] = std::array::from_fn(|_| ((g.u64() % 65535) as i32 - 32767) as i16);
    v.sort_unstable_by(|a, b| b.cmp(a));
    v
}

/// A C codebook, its columns, its size, `use_half`.
type VqCase<const D: usize> = (&'static [i16], [&'static [i16]; D], i16, bool);

/// Autocorr → Lag_window → Levinson (stateful over a signal sequence) on all three windows.
#[test]
fn lp_analysis_matches_c() {
    let wins: [&[i16; L_WINDOW]; 3] = [&WINDOW_200_40, &WINDOW_160_80, &WINDOW_232_8];
    let mut g = Rng::new(3);
    let (mut rs, mut cs) = (LevinsonState::default(), LevinsonState::default().old_a);
    for _ in 0..N / 10 {
        let mut x = [0i16; L_WINDOW];
        signal(&mut g, &mut x);
        let w = wins[(g.u64() % 3) as usize];
        let (mut rh, mut rl, mut ch, mut cl) = ([0i16; MP1], [0i16; MP1], [0i16; MP1], [0i16; MP1]);
        let (mut ov, mut o) = (false, 0);
        let rn = autocorr(&x, M, &mut rh, &mut rl, w);
        let cn = unsafe {
            c_autocorr(
                x.as_mut_ptr(),
                M as i16,
                ch.as_mut_ptr(),
                cl.as_mut_ptr(),
                w.as_ptr(),
                &mut o,
            )
        };
        assert_eq!((rn, rh, rl), (cn, ch, cl), "Autocorr {x:?}");
        lag_window(M, &mut rh, &mut rl, &mut ov);
        unsafe { c_lag_window(M as i16, ch.as_mut_ptr(), cl.as_mut_ptr(), &mut o) };
        assert_eq!((rh, rl), (ch, cl), "Lag_window");
        let (mut ra, mut ca, mut rrc, mut crc) = ([0i16; MP1], [0i16; MP1], [0i16; 4], [0i16; 4]);
        levinson(&mut rs, &rh, &rl, &mut ra, &mut rrc, &mut ov);
        unsafe {
            c_levinson(
                &mut cs,
                ch.as_mut_ptr(),
                cl.as_mut_ptr(),
                ca.as_mut_ptr(),
                crc.as_mut_ptr(),
                &mut o,
            )
        };
        assert_eq!((ra, rrc, rs.old_a), (ca, crc, cs), "Levinson {rh:?} {rl:?}");
    }
}

/// `Chebps` (b1 as one word instead of the C's hi/lo split) vs the C at every x: Az_lsp's
/// own f1/f2 (realistic), and full-scale random coefficients (wrapping sums, saturated
/// outputs). Az_lsp's diff alone can't see it: a low-bit change rarely moves a root.
#[test]
fn chebps_matches_c() {
    use crate::basic_op::{l_extract, l_mac, l_msu, l_mult, mpy_32_16};
    // Whether TS 26.073's chain saturates anywhere (L_shl doesn't flag: checked by hand),
    // or ends on opencore's off-by-one bound (t0 = 0x01ff_ffff, which it saturates to -32768).
    // There opencore differs from the reference; that class is checked through `Az_lsp`
    // against the reference (`az_lsp_matches_c`).
    fn spec_only(x: i16, f: &[i16; 6], n: usize) -> bool {
        let o = &mut false;
        let shl_sat = |t: i32, k: u32| (t << k) >> k != t;
        let (mut b2_h, mut b2_l) = (256i16, 0i16);
        let mut t0 = l_mac(l_mult(x, 512, o), f[1], 8192, o);
        let (mut b1_h, mut b1_l) = l_extract(t0, o);
        let mut sat = false;
        for &fi in &f[2..n] {
            let m = mpy_32_16(b1_h, b1_l, x, o);
            sat |= shl_sat(m, 1);
            t0 = l_mac(m << 1, b2_h, i16::MIN, o);
            t0 = l_msu(t0, b2_l, 1, o);
            t0 = l_mac(t0, fi, 8192, o);
            let (h, l) = l_extract(t0, o);
            (b2_h, b2_l, b1_h, b1_l) = (b1_h, b1_l, h, l);
        }
        t0 = mpy_32_16(b1_h, b1_l, x, o);
        t0 = l_mac(t0, b2_h, i16::MIN, o);
        t0 = l_msu(t0, b2_l, 1, o);
        t0 = l_mac(t0, f[n], 4096, o);
        *o || sat || shl_sat(t0, 6) || t0 == 0x01ff_ffff
    }
    let mut g = Rng::new(31);
    let (mut compared, mut skipped) = (0, 0);
    for n in 0..64 {
        let f: [i16; 6] = if n % 2 == 0 {
            let a = crate::post::test_gen::coefs(&mut g);
            let mut f = [1024i16; 6];
            for i in 0..5 {
                f[i + 1] = (((a[i + 1] as i32 + a[M - i] as i32) >> 2) as i16).wrapping_sub(f[i]);
            }
            f
        } else {
            std::array::from_fn(|_| g.i16())
        };
        for x in i16::MIN..=i16::MAX {
            // The fast forms against the operator chain, saturating classes included.
            assert_eq!(
                super::super::lpc::chebps(x, &f, 5),
                super::super::lpc::chebps_exact(x, &f, 5),
                "x {x} f {f:?}"
            );
            if spec_only(x, &f, 5) {
                skipped += 1;
                continue;
            }
            let r = super::super::lpc::chebps(x, &f, 5);
            let (mut cf, mut o) = (f, 0);
            let c = unsafe { c_chebps(x, cf.as_mut_ptr(), 5, &mut o) };
            assert_eq!(r, c, "x {x} f {f:?}");
            compared += 1;
        }
    }
    assert!(
        compared > 1_000_000 && skipped > 1000,
        "compared: {compared}, reference-only: {skipped}"
    );
}

#[test]
fn az_lsp_matches_c() {
    let mut g = Rng::new(4);
    for _ in 0..N {
        let mut a = crate::post::test_gen::coefs(&mut g);
        let mut old = [0i16; M];
        old.iter_mut().for_each(|v| *v = g.i16());
        let (mut r, mut c) = ([0i16; M], [0i16; M]);
        az_lsp(&a, &mut r, &old);
        let mut o = 0;
        unsafe { c_az_lsp(a.as_mut_ptr(), c.as_mut_ptr(), old.as_mut_ptr(), &mut o) };
        assert_eq!(r, c, "{a:?}");
    }
}

#[test]
fn lsf_quantisers_match_c() {
    let mut g = Rng::new(5);
    let (mut rs, mut cs) = (QPlsfState::default(), [0i16; M]);
    for _ in 0..N / 20 {
        let mut l1 = lsp(&mut g);
        let mut l2 = lsp(&mut g);
        let (mut ov, mut o) = (false, 0);
        if g.u64().is_multiple_of(3) {
            let (mut rq1, mut rq2, mut cq1, mut cq2) = ([0i16; M], [0i16; M], [0i16; M], [0i16; M]);
            let (mut ri, mut ci) = ([0i16; 5], [0i16; 5]);
            q_plsf_5(&mut rs, &l1, &l2, &mut rq1, &mut rq2, &mut ri, &mut ov);
            unsafe {
                c_q_plsf_5(
                    &mut cs,
                    l1.as_mut_ptr(),
                    l2.as_mut_ptr(),
                    cq1.as_mut_ptr(),
                    cq2.as_mut_ptr(),
                    ci.as_mut_ptr(),
                    &mut o,
                )
            };
            assert_eq!(
                (rq1, rq2, ri, rs.past_rq),
                (cq1, cq2, ci, cs),
                "Q_plsf_5 {l1:?} {l2:?}"
            );
        } else {
            let mode = CodecMode::ALL[(g.u64() % 9) as usize];
            if mode == CodecMode::MR122 {
                continue;
            }
            let (mut rq, mut cq, mut ri, mut ci) = ([0i16; M], [0i16; M], [0i16; 3], [0i16; 3]);
            let (mut rp, mut cp) = (0i16, 0i16);
            q_plsf_3(&mut rs, mode, &l1, &mut rq, &mut ri, &mut rp, &mut ov);
            unsafe {
                c_q_plsf_3(
                    &mut cs,
                    mode as i32,
                    l1.as_mut_ptr(),
                    cq.as_mut_ptr(),
                    ci.as_mut_ptr(),
                    &mut cp,
                    &mut o,
                )
            };
            assert_eq!(
                (rq, ri, rp, rs.past_rq),
                (cq, ci, cp, cs),
                "Q_plsf_3 {mode:?} {l1:?}"
            );
        }
    }
}

/// MR122's split VQs in lanes (`vq_pair`, `vq_pair_s`; `vq5` with its fallback) vs the C's
/// loops with early exits (`vq_subvec`, `vq_subvec_s`, checked against the C via
/// `Q_plsf_5`): realistic residuals, full-scale ones (wrapping sums: the lanes must decline),
/// zero/tiny weights (ties everywhere: first index).
#[test]
fn vq_pair_matches_loops() {
    use crate::encoder::lpc::{splits5, vq_pair, vq_pair_s, vq_subvec, vq_subvec_s, vq5};
    let mut g = Rng::new(30);
    let (mut declined, mut picked) = (0, 0);
    for n in 0..N / 4 {
        let class = n % 4;
        let v = |g: &mut Rng| -> [i16; 2] {
            std::array::from_fn(|_| match class {
                1 => g.u64() as i16,
                _ => (g.u64() % 6000) as i16 - 3000,
            })
        };
        let w = |g: &mut Rng| -> [i16; 2] {
            std::array::from_fn(|_| match class {
                1 => g.u64() as i16,
                2 => (g.u64() % 3) as i16,
                3 => [0, 32767, -32768][g.u64() as usize % 3],
                _ => 3000 + (g.u64() % 25000) as i16,
            })
        };
        let (r1, r2, w1, w2) = (v(&mut g), v(&mut g), w(&mut g), w(&mut g));
        for (dico, c, signed) in splits5() {
            let size = c[0].len();
            let (mut a1, mut a2, mut b1, mut b2) = (r1, r2, r1, r2);
            let want = if signed {
                vq_subvec_s(&mut a1, &mut a2, dico, &w1, &w2, size)
            } else {
                vq_subvec(&mut a1, &mut a2, dico, &w1, &w2, size)
            };
            let lanes = if signed {
                vq_pair_s(&r1, &r2, &c, &w1, &w2)
            } else {
                vq_pair(&r1, &r2, &c, &w1, &w2)
            };
            match lanes {
                Some(i) => {
                    picked += 1;
                    assert_eq!(
                        i as i16, want,
                        "{signed} {size} {r1:?} {r2:?} {w1:?} {w2:?}"
                    );
                }
                None => declined += 1,
            }
            let got = vq5(&mut b1, &mut b2, dico, c, &w1, &w2, signed);
            assert_eq!((got, b1, b2), (want, a1, a2), "{signed} {size}");
        }
    }
    assert!(
        declined > N / 50 && picked > N / 2,
        "classes: {declined} {picked}"
    );
}

/// The split-VQ searches in lanes vs the C on every codebook: full-scale residuals and
/// weights (wrapping distance sums), zero and tiny weights (ties everywhere: first index).
#[test]
fn vq_subvec_matches_c() {
    let mut g = Rng::new(26);
    let (mut ties, mut wraps) = (0, 0);
    for n in 0..N / 4 {
        let class = n % 5;
        let mut r: [i16; 4] = std::array::from_fn(|k| match class {
            0 => g.u64() as i16,
            // Opposite the first codeword (r - d wraps near -32768): its sum wraps.
            4 => DICO1_LSF_3[k % 3].wrapping_add(-32768 + (g.u64() % 64) as i16),
            _ => (g.u64() % 4000) as i16 - 2000,
        });
        let mut wf: [i16; 4] = std::array::from_fn(|_| match class {
            0 => g.u64() as i16,
            1 => (g.u64() % 3) as i16,
            2 => [0, -32768, 32767][g.u64() as usize % 3],
            4 => [-32768, 32767][g.u64() as usize % 2],
            _ => 3000 + (g.u64() % 20000) as i16,
        });
        // Class coverage: all-zero weights (every distance 0: index 0), and a wrapping
        // distance sum (exact sum >= 2^31) on the first codeword.
        ties += (wf[..3] == [0, 0, 0]) as usize;
        let exact: i64 = (0..3)
            .map(|k| {
                let t =
                    ((wf[k] as i32 * r[k].wrapping_sub(DICO1_LSF_3[k]) as i32) >> 15) as i16 as i64;
                t * t
            })
            .sum();
        wraps += (exact >= 1 << 31) as usize;
        let three: [VqCase<3>; 4] = [
            (&DICO1_LSF_3, cols(&DICO1_SOA), 256, false),
            (&DICO2_LSF_3, cols(&DICO2_SOA), 512, false),
            (&DICO2_LSF_3, cols(&DICO2_HALF_SOA), 256, true),
            (&MR795_1_LSF, cols(&MR795_1_SOA), 512, false),
        ];
        for (d, c, size, half) in three {
            let (mut rr, mut o) = (r, 0);
            let i = vq_subvec3(&mut rr, &c, &wf);
            let mut rs = r; // the portable build's SSE2 body, on every path
            assert_eq!(
                crate::simd::vq_subvec_sse2(&mut rs, &c, &wf),
                i,
                "sse2 {size}"
            );
            assert_eq!(rs, rr);
            let ci = unsafe {
                c_vq_subvec3(
                    r.as_mut_ptr(),
                    d.as_ptr(),
                    wf.as_mut_ptr(),
                    size,
                    half as i32,
                    &mut o,
                )
            };
            assert_eq!((i, &rr[..3]), (ci, &r[..3]), "{size} {half} {wf:?}");
        }
        let four: [VqCase<4>; 2] = [
            (&DICO3_LSF_3, cols(&DICO3_SOA), 512, false),
            (&MR515_3_LSF, cols(&MR515_3_SOA), 128, false),
        ];
        for (d, c, size, _) in four {
            let (mut rr, mut o) = (r, 0);
            let i = vq_subvec4(&mut rr, &c, &wf);
            let mut rs = r;
            assert_eq!(
                crate::simd::vq_subvec_sse2(&mut rs, &c, &wf),
                i,
                "sse2 {size}"
            );
            assert_eq!(rs, rr);
            let ci =
                unsafe { c_vq_subvec4(r.as_mut_ptr(), d.as_ptr(), wf.as_mut_ptr(), size, &mut o) };
            assert_eq!((i, rr), (ci, r), "{size} {wf:?}");
        }
    }
    assert!(
        ties > N / 1000 && wraps > N / 1000,
        "classes: {ties} {wraps}"
    );
}
