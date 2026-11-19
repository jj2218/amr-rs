//! Differential tests of `gain_q` against the TS 26.073 reference (`ref26073-sys`;
//! helpers in `encoder::vs_ref`).
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::basic_op::Flag;
use crate::cnst::{CodecMode, L_SUBFR};
use crate::encoder::gain_q::qua_gain;
use crate::encoder::gain_q::{
    GainAdaptState, calc_filt_energies, mr475_gain_quant, mr795_gain_quant,
};
use crate::encoder::gain_q::{g_code, q_gain_code};
use crate::encoder::vs_ref::CFlag;
use crate::encoder::vs_ref::CommonAmrTbls;
use crate::encoder::vs_ref::refc::*;
use crate::encoder::vs_ref::tables;
use crate::encoder::vs_ref::*;
use crate::gains::GcPredState;
use crate::testutil::{N, Rng};
use ref26073_sys::lock;

#[test]
fn g_code_matches_ref() {
    let _l = lock();
    let mut g = Rng::new(1);
    for _ in 0..N {
        let mut xn2 = [0i16; L_SUBFR];
        let mut y2 = [0i16; L_SUBFR];
        signal(&mut g, &mut xn2);
        signal(&mut g, &mut y2);
        if y2.iter().all(|&v| v >> 1 == 0) {
            continue; // <y2, y2> = 0: the reference's div_s aborts
        }
        let c = unsafe { ref_G_code(xn2.as_ptr(), y2.as_ptr()) };
        assert_eq!(g_code(&xn2, &y2), c, "xn2 {xn2:?} y2 {y2:?}");
    }
}

unsafe extern "C" {
    fn ref_G_code(xn2: *const i16, y2: *const i16) -> i16;
    fn ref_q_gain_code(
        mode: i32,
        exp: i16,
        frac: i16,
        gain: *mut i16,
        qe122: *mut i16,
        qe: *mut i16,
    ) -> i16;
    fn ref_Qua_gain(
        mode: i32,
        exp_gcode0: i16,
        frac_gcode0: i16,
        frac: *mut i16,
        exp: *mut i16,
        gp_limit: i16,
        gp: *mut i16,
        gc: *mut i16,
        q1: *mut i16,
        q2: *mut i16,
    ) -> i16;
    fn ref_MR475_gain_quant(
        pred: *mut [i16; 8],
        e0: i16,
        f0: i16,
        ec0: *mut i16,
        fc0: *mut i16,
        te0: i16,
        tf0: i16,
        code: *mut i16,
        e1: i16,
        f1: i16,
        ec1: *mut i16,
        fc1: *mut i16,
        te1: i16,
        tf1: i16,
        gp_limit: i16,
        p0: *mut i16,
        c0: *mut i16,
        p1: *mut i16,
        c1: *mut i16,
    ) -> i16;
    fn ref_calc_filt_energies(
        mode: i32,
        xn: *mut i16,
        xn2: *mut i16,
        y1: *mut i16,
        y2: *mut i16,
        g_coeff: *mut i16,
        frac: *mut i16,
        exp: *mut i16,
        cgf: *mut i16,
        cge: *mut i16,
    );
    fn ref_MR795_gain_quant(
        st: *mut GainAdaptState,
        res: *mut i16,
        exc: *mut i16,
        code: *mut i16,
        frac_coeff: *mut i16,
        exp_coeff: *mut i16,
        exp_code_en: i16,
        frac_code_en: i16,
        exp_gcode0: i16,
        frac_gcode0: i16,
        l_subfr: i16,
        cod_gain_frac: i16,
        cod_gain_exp: i16,
        gp_limit: i16,
        gain_pit: *mut i16,
        gain_cod: *mut i16,
        qua_ener_mr122: *mut i16,
        qua_ener: *mut i16,
        anap: *mut *mut i16,
    );
}

unsafe fn c_qua_gain(
    mode: i32,
    exp_gcode0: i16,
    frac_gcode0: i16,
    frac: *mut i16,
    exp: *mut i16,
    gp_limit: i16,
    gp: *mut i16,
    gc: *mut i16,
    q1: *mut i16,
    q2: *mut i16,
    _t: *mut crate::encoder::vs_ref::CommonAmrTbls,
    o: *mut CFlag,
) -> i16 {
    unsafe {
        call(o, || {
            ref_Qua_gain(
                mode,
                exp_gcode0,
                frac_gcode0,
                frac,
                exp,
                gp_limit,
                gp,
                gc,
                q1,
                q2,
            )
        })
    }
}

unsafe fn c_mr475_gain_quant(
    pred: *mut [i16; 8],
    e0: i16,
    f0: i16,
    ec0: *mut i16,
    fc0: *mut i16,
    te0: i16,
    tf0: i16,
    code: *mut i16,
    e1: i16,
    f1: i16,
    ec1: *mut i16,
    fc1: *mut i16,
    te1: i16,
    tf1: i16,
    gp_limit: i16,
    p0: *mut i16,
    c0: *mut i16,
    p1: *mut i16,
    c1: *mut i16,
    o: *mut CFlag,
) -> i16 {
    unsafe {
        call(o, || {
            ref_MR475_gain_quant(
                pred, e0, f0, ec0, fc0, te0, tf0, code, e1, f1, ec1, fc1, te1, tf1, gp_limit, p0,
                c0, p1, c1,
            )
        })
    }
}

unsafe fn c_calc_filt_energies(
    mode: i32,
    xn: *mut i16,
    xn2: *mut i16,
    y1: *mut i16,
    y2: *mut i16,
    g_coeff: *mut i16,
    frac: *mut i16,
    exp: *mut i16,
    cgf: *mut i16,
    cge: *mut i16,
    o: *mut CFlag,
) {
    unsafe {
        call(o, || {
            ref_calc_filt_energies(mode, xn, xn2, y1, y2, g_coeff, frac, exp, cgf, cge)
        })
    }
}

/// `Qua_gain` inputs: realistic coefficient exponents, sometimes extreme (saturating) ones.
pub(crate) struct GainIn {
    pub(crate) mode: CodecMode,
    pub(crate) exp_gcode0: i16,
    pub(crate) frac_gcode0: i16,
    pub(crate) frac: [i16; 5],
    pub(crate) exp: [i16; 5],
    pub(crate) gp_limit: i16,
}

pub(crate) fn gain_in(g: &mut Rng) -> GainIn {
    let class = g.u64() % 4;
    let wild = class == 0;
    let mut frac = [0i16; 5];
    let mut exp = [0i16; 5];
    for (f, e) in frac.iter_mut().zip(&mut exp) {
        *f = g.u64() as i16;
        *e = if wild {
            (g.u64() % 40) as i16 - 5
        } else {
            10 + (g.u64() % 12) as i16
        };
    }
    let exp_gcode0 = (g.u64() % if wild { 30 } else { 16 }) as i16;
    if class == 1 {
        // Saturating: equal normalised exponents (each coeff ~ +-16384), large same-sign
        // fractions, so the five Mpy_32_16 terms overflow the 32-bit sum.
        let (t, ec) = (5 + (g.u64() % 8) as i16, exp_gcode0 - 11);
        exp = [t + 13, t + 14, t - 2 * ec - 15, t - ec, t - ec - 1];
        let neg = g.u64() & 1 == 1;
        for f in &mut frac {
            let m = 28000 + (g.u64() % 4768) as i16;
            *f = if neg { -m } else { m };
        }
    }
    GainIn {
        mode: mode_of(g),
        exp_gcode0,
        frac_gcode0: (g.u64() % 32768) as i16,
        frac,
        exp,
        gp_limit: [32767, 16383, (g.u64() % 32768) as i16][g.u64() as usize % 3],
    }
}

pub(crate) fn qua_gain_rust(x: &GainIn) -> [i16; 5] {
    let (mut gp, mut gc, mut q1, mut q2, mut ov) = (0, 0, 0, 0, false);
    let i = qua_gain(
        x.mode,
        x.exp_gcode0,
        x.frac_gcode0,
        &x.frac,
        &x.exp,
        x.gp_limit,
        &mut gp,
        &mut gc,
        &mut q1,
        &mut q2,
        &mut ov,
    );
    let _ = ov;
    [i, gp, gc, q1, q2]
}

pub(crate) fn qua_gain_c(x: &GainIn, t: &mut CommonAmrTbls) -> [i16; 5] {
    let (mut gp, mut gc, mut q1, mut q2, mut o) = (0, 0, 0, 0, 0);
    let (mut frac, mut exp) = (x.frac, x.exp);
    let i = unsafe {
        c_qua_gain(
            x.mode as i32,
            x.exp_gcode0,
            x.frac_gcode0,
            frac.as_mut_ptr(),
            exp.as_mut_ptr(),
            x.gp_limit,
            &mut gp,
            &mut gc,
            &mut q1,
            &mut q2,
            t,
            &mut o,
        )
    };
    let _ = o;
    [i, gp, gc, q1, q2]
}

#[test]
fn q_gain_code_matches_ref() {
    let _l = lock();
    let mut g = Rng::new(2);
    for _ in 0..N {
        let mode = if g.u64() & 1 == 0 {
            CodecMode::MR122
        } else {
            CodecMode::MR795
        };
        let exp = (g.u64() % 20) as i16;
        let frac = g.u64() as i16 & 0x7fff;
        let gain0 = g.i16();
        let (mut gc, mut a, mut b) = (gain0, 0, 0);
        let ci = unsafe { ref_q_gain_code(mode as i32, exp, frac, &mut gc, &mut a, &mut b) };
        let (mut gr, mut ra, mut rb, mut ov): (i16, i16, i16, Flag) = (gain0, 0, 0, false);
        let ri = q_gain_code(mode, exp, frac, &mut gr, &mut ra, &mut rb, &mut ov);
        assert_eq!(
            (ri, gr, ra, rb),
            (ci, gc, a, b),
            "{mode:?} exp {exp} frac {frac} gain {gain0}"
        );
    }
}

#[test]
fn qua_gain_matches_c() {
    let mut g = Rng::new(24);
    let mut t = tables();
    for _ in 0..N {
        let x = gain_in(&mut g);
        assert_eq!(
            qua_gain_rust(&x),
            qua_gain_c(&x, &mut t),
            "{:?} {:?}",
            x.frac,
            x.exp
        );
    }
}

/// `MR475_gain_quant` (its search runs in lanes) vs the C: both subframes' `Qua_gain`-style
/// inputs (incl. the extreme-exponent class: wrapping sums), random predictor state and
/// target energies, MR475-like innovations. Ties: equal distances in the wild class.
#[test]
fn mr475_gain_quant_matches_c() {
    let mut g = Rng::new(27);
    let (mut ties, mut negative) = (0, 0);
    for n in 0..N {
        let (mut a, mut b) = (gain_in(&mut g), gain_in(&mut g));
        if n % 16 == 0 {
            // All distances 0: ties everywhere, the first row within gp_limit wins.
            (a.frac, b.frac) = ([0; 5], [0; 5]);
            ties += 1;
        }
        // MR475's innovation: 2 pulses of +-8191 (a full-band vector drives gc_pred's
        // exponent into the C's undefined L_shr by >= 32).
        let mut code = [0i16; 40];
        for _ in 0..2 {
            code[(g.u64() % 40) as usize] = if g.u64() & 1 == 1 { 8191 } else { -8191 };
        }
        let mut pred = [0i16; 8];
        for (i, p) in pred.iter_mut().enumerate() {
            *p = if i < 4 {
                -5443 + (g.u64() % 8480) as i16
            } else {
                g.u64() as i16
            };
        }
        let te = [(g.u64() % 30) as i16, (g.u64() % 30) as i16];
        let tf = [(g.u64() % 32768) as i16, (g.u64() % 32768) as i16];
        // Negative limits too (no row qualifies: index 0).
        let gp_limit = if n % 8 == 1 {
            g.u64() as i16
        } else {
            a.gp_limit
        };
        negative += (gp_limit < 0) as usize;
        let mut st = GcPredState {
            past_qua_en: pred[..4].try_into().unwrap(),
            past_qua_en_mr122: pred[4..].try_into().unwrap(),
        };
        let (mut r, mut ov) = ([0i16; 4], false);
        let ri = {
            let [p0, c0, p1, c1] = &mut r;
            mr475_gain_quant(
                &mut st,
                a.exp_gcode0,
                a.frac_gcode0,
                &a.exp,
                &a.frac,
                te[0],
                tf[0],
                &code,
                b.exp_gcode0,
                b.frac_gcode0,
                &b.exp,
                &b.frac,
                te[1],
                tf[1],
                gp_limit,
                p0,
                c0,
                p1,
                c1,
                &mut ov,
            )
        };
        let (mut c, mut o, mut cp) = ([0i16; 4], 0, pred);
        let (mut a2, mut b2, mut code2) = ((a.exp, a.frac), (b.exp, b.frac), code);
        let ci = unsafe {
            let [p0, c0, p1, c1] = &mut c;
            c_mr475_gain_quant(
                &mut cp,
                a.exp_gcode0,
                a.frac_gcode0,
                a2.0.as_mut_ptr(),
                a2.1.as_mut_ptr(),
                te[0],
                tf[0],
                code2.as_mut_ptr(),
                b.exp_gcode0,
                b.frac_gcode0,
                b2.0.as_mut_ptr(),
                b2.1.as_mut_ptr(),
                te[1],
                tf[1],
                gp_limit,
                p0,
                c0,
                p1,
                c1,
                &mut o,
            )
        };
        let rs: Vec<i16> = st
            .past_qua_en
            .iter()
            .chain(&st.past_qua_en_mr122)
            .copied()
            .collect();
        let _ = (ov, o);
        assert_eq!(
            (ri, r, &rs[..]),
            (ci, c, &cp[..]),
            "{:?} {:?}",
            a.exp,
            a.frac
        );
    }
    assert!(
        ties > N / 20 && negative > N / 40,
        "classes: {ties} {negative}"
    );
}

/// `MR795_gain_quant` (its `MR795_gain_code_quant3` search runs in lanes) vs the C:
/// `Qua_gain`-style inputs incl. the saturating class (the lanes decline, the C's loop
/// runs), a tie class (zero codebook-gain coefficients: per candidate all 32 distances
/// equal, the first wins), random adaptation state, residual, excitation and innovation.
#[test]
fn mr795_gain_quant_matches_ref() {
    let _l = ref26073_sys::lock();
    let mut g = Rng::new(29);
    let (mut ties, mut sat, mut unfilt) = (0, 0, 0);
    for n in 0..N {
        let mut x = gain_in(&mut g);
        if n % 16 == 0 {
            x.frac[2..].fill(0);
            ties += 1;
        } else if n % 4 == 1 {
            // Saturating: `gain_in`'s class 1 with MR795's exponent offsets (exp_code =
            // exp_gcode0 - 10): five same-sign terms of ~2^30 wrap the 32-bit sum.
            let (t, ec) = (5 + (g.u64() % 8) as i16, x.exp_gcode0 - 10);
            x.exp = [t + 13, t + 14, t - 2 * ec - 15, t - ec, t - ec - 1];
            let neg = g.u64() & 1 == 1;
            for f in &mut x.frac {
                let m = 28000 + (g.u64() % 4768) as i16;
                *f = if neg { -m } else { m };
            }
        }
        let mut v = [[0i16; 40]; 3];
        for s in v.iter_mut() {
            signal(&mut g, s);
        }
        if n % 8 == 2 {
            // calc_unfilt_energies' fallback: res at -32768 where exc is 0 gives an LTP
            // residual of -32768 (L_mult saturates in the C's L_mac chain), full-scale res
            // elsewhere pushes its energy past 2^31.
            for i in 0..40 {
                if g.u64().is_multiple_of(4) {
                    (v[0][i], v[1][i]) = (i16::MIN, 0);
                } else {
                    v[0][i] = g.i16();
                }
            }
            unfilt += 1;
        }
        let mut st = GainAdaptState {
            onset: (g.u64() % 3) as i16,
            prev_alpha: (g.u64() % 32768) as i16,
            prev_gc: (g.u64() % 32768) as i16,
            ltpg_mem: std::array::from_fn(|_| g.u64() as i16 >> 2),
        };
        let mut cst = st.clone();
        let gp_limit = [32767, 15565, 15565 + (g.u64() % 17203) as i16][n % 3];
        let gain_pit = (g.u64() % 32768) as i16;
        let (ece, fce) = ((g.u64() % 30) as i16 - 5, (g.u64() % 32768) as i16);
        let (cgf, cge) = ((g.u64() % 32768) as i16, (g.u64() % 30) as i16 - 5);
        let (mut r, mut ana, mut pos, mut ov) = ([gain_pit, 0, 0, 0], [0i16; 2], 0, false);
        {
            let [gp, gc, q1, q2] = &mut r;
            mr795_gain_quant(
                &mut st,
                &v[0],
                &v[1],
                &v[2],
                &x.frac,
                &x.exp,
                ece,
                fce,
                x.exp_gcode0,
                x.frac_gcode0,
                cgf,
                cge,
                gp_limit,
                gp,
                gc,
                q1,
                q2,
                &mut ana,
                &mut pos,
                &mut ov,
            );
        }
        let (mut c, mut cana, mut cv) = ([gain_pit, 0, 0, 0], [0i16; 2], v);
        let mut anap = cana.as_mut_ptr();
        ref26073_sys::clear_overflow();
        unsafe {
            let [gp, gc, q1, q2] = &mut c;
            let (mut f, mut e) = (x.frac, x.exp);
            let [a, b, d] = &mut cv;
            ref_MR795_gain_quant(
                &mut cst,
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                d.as_mut_ptr(),
                f.as_mut_ptr(),
                e.as_mut_ptr(),
                ece,
                fce,
                x.exp_gcode0,
                x.frac_gcode0,
                40,
                cgf,
                cge,
                gp_limit,
                gp,
                gc,
                q1,
                q2,
                &mut anap,
            );
        }
        sat += ref26073_sys::overflow() as usize;
        let _ = ov;
        assert_eq!((r, ana, &st), (c, cana, &cst), "{:?} {:?}", x.frac, x.exp);
    }
    assert!(
        ties > N / 20 && sat > N / 20 && unfilt > N / 10,
        "classes: {ties} {sat} {unfilt}"
    );
}

/// `calc_filt_energies` (its `L_mac` chains as plain sums under a no-overflow bound) vs the
/// C, all modes (MR475/MR795 also give the codebook gain): speech-like vectors, and
/// full-scale `xn`/`y1` against large `y2` (the bound fails, chains saturate: the flag).
#[test]
fn calc_filt_energies_matches_c() {
    let mut g = Rng::new(32);
    let mut sat = 0;
    for n in 0..N {
        let mut v = [[0i16; 40]; 4];
        for s in v.iter_mut() {
            signal(&mut g, s);
        }
        if n % 4 == 0 {
            let sign = if g.u64() & 1 == 0 { 1 } else { -1 };
            for i in 0..40 {
                let a = 32767 - (g.u64() % 64) as i16;
                (v[0][i], v[2][i]) = (sign * a, -sign * a);
                v[3][i] = if g.u64().is_multiple_of(8) {
                    g.i16()
                } else {
                    32767 - (g.u64() % 99) as i16
                };
            }
        }
        let mode = mode_of(&mut g);
        let gc: [i16; 4] = std::array::from_fn(|_| g.i16());
        let (mut rf, mut re, mut rcf, mut rce, mut ov) = ([0i16; 5], [0i16; 5], 0, 0, false);
        calc_filt_energies(
            mode, &v[0], &v[1], &v[2], &v[3], &gc, &mut rf, &mut re, &mut rcf, &mut rce, &mut ov,
        );
        let (mut cf, mut ce, mut ccf, mut cce, mut o, mut cv, mut cg) =
            ([0i16; 5], [0i16; 5], 0, 0, 0, v, gc);
        let [a, b, c, d] = &mut cv;
        unsafe {
            c_calc_filt_energies(
                mode as i32,
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                c.as_mut_ptr(),
                d.as_mut_ptr(),
                cg.as_mut_ptr(),
                cf.as_mut_ptr(),
                ce.as_mut_ptr(),
                &mut ccf,
                &mut cce,
                &mut o,
            )
        };
        sat += (o != 0) as usize;
        let _ = ov;
        assert_eq!((rf, re, rcf, rce), (cf, ce, ccf, cce), "{mode:?} {v:?}");
    }
    assert!(sat > N / 8, "saturating: {sat}");
}

/// `MR795_gain_code_quant3`'s lane search vs the C's loop (`mr795_search_c`, inside the C
/// diff above) on direct inputs: realistic coefficients (accepted), full-scale coefficients
/// and `l_tmp0` near +-2^31 (sums saturate: must decline), zero codebook-gain coefficients (hi
/// and lo)
/// with equal `l_tmp0` (all 96 distances tie: index 0). An accepted search must equal the
/// loop and the loop must not have saturated anywhere (else the lanes missed a flag).
#[test]
fn mr795_search_matches_loop() {
    use crate::encoder::gain_q::mr795_search_c;
    let mut g = Rng::new(33);
    let (mut declined, mut ties) = (0, 0);
    for n in 0..N {
        let class = n % 4;
        let mut coeff: [i16; 5] = std::array::from_fn(|_| match class {
            0 => g.i16(),
            _ => g.i16() >> (4 + g.u64() % 8),
        });
        let mut coeff_lo: [i16; 5] = std::array::from_fn(|_| (g.u64() % 32768) as i16);
        let mut l_tmp0: [i32; 3] = std::array::from_fn(|_| match class {
            0 => g.u64() as i32,
            _ => (g.u64() as i32) >> 8,
        });
        if class == 2 {
            coeff[2..].fill(0);
            coeff_lo[2..].fill(0);
            l_tmp0 = [l_tmp0[0]; 3];
            ties += 1;
        }
        let gcode0 = (g.u64() % 32768) as i16;
        let gp: [i16; 3] = std::array::from_fn(|_| (g.u64() % 32768) as i16);
        let mut ov = false;
        let want = mr795_search_c(&coeff, &coeff_lo, &l_tmp0, gcode0, &gp, &mut ov);
        match crate::simd::mr795_search(&coeff, &coeff_lo, &l_tmp0, gcode0, &gp) {
            Some(k) => {
                assert_eq!((k, ov), (want, false), "{coeff:?} {l_tmp0:?}");
                if class == 2 {
                    assert_eq!(k, 0);
                }
            }
            None => declined += 1,
        }
    }
    assert!(
        declined > N / 20 && declined < N / 2 && ties > N / 8,
        "classes: {declined} {ties}"
    );
}
