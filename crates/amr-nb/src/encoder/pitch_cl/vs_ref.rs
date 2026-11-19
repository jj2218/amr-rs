//! Differential tests of `pitch_cl` against the TS 26.073 reference (`ref26073-sys`;
//! helpers in `encoder::vs_ref`).
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::cnst::L_INTERPOL;
use crate::cnst::L_SUBFR;
use crate::cnst::{CodecMode, PIT_MAX};
use crate::encoder::pitch_cl::convolve;
use crate::encoder::pitch_cl::g_pitch;
use crate::encoder::pitch_cl::{PitchFrState, pitch_fr};
use crate::encoder::vs_ref::CFlag;
use crate::encoder::vs_ref::refc::*;
use crate::encoder::vs_ref::*;
use crate::testutil::{N, Rng};

#[test]
fn g_pitch_matches_c() {
    let mut g = Rng::new(7);
    for _ in 0..N {
        let (mut xn, mut y1) = ([0i16; L_SUBFR], [0i16; L_SUBFR]);
        signal(&mut g, &mut xn);
        signal(&mut g, &mut y1);
        let mode = mode_of(&mut g);
        let (mut rg, mut cg) = ([0i16; 6], [0i16; 6]);
        let (mut ov, mut o) = (false, 0);
        let r = g_pitch(mode, &xn, &y1, &mut rg, L_SUBFR, &mut ov);
        let c = unsafe {
            c_g_pitch(
                mode as i32,
                xn.as_mut_ptr(),
                y1.as_mut_ptr(),
                cg.as_mut_ptr(),
                L_SUBFR as i16,
                &mut o,
            )
        };
        let _ = (ov, o); // the reference's flag is global and sticky: values only
        assert_eq!((r, rg), (c, cg), "{mode:?} {xn:?} {y1:?}");
    }
}

unsafe extern "C" {
    fn ref_Convolve(x: *mut i16, h: *mut i16, y: *mut i16, l: i16);
    fn ref_Pitch_fr(
        st: *mut i16,
        mode: i32,
        t_op: *mut i16,
        exc: *mut i16,
        xn: *mut i16,
        h: *mut i16,
        l_subfr: i16,
        i_subfr: i16,
        frac: *mut i16,
        resu3: *mut i16,
        idx: *mut i16,
    ) -> i16;
    fn ref_G_pitch(mode: i32, xn: *mut i16, y1: *mut i16, g: *mut i16, l: i16) -> i16;
}

unsafe fn c_convolve(x: *mut i16, h: *mut i16, y: *mut i16, l: i16) {
    unsafe { call(std::ptr::null_mut(), || ref_Convolve(x, h, y, l)) }
}

unsafe fn c_pitch_fr(
    st: *mut i16,
    mode: i32,
    t_op: *mut i16,
    exc: *mut i16,
    xn: *mut i16,
    h: *mut i16,
    l_subfr: i16,
    i_subfr: i16,
    frac: *mut i16,
    resu3: *mut i16,
    idx: *mut i16,
    o: *mut CFlag,
) -> i16 {
    unsafe {
        call(o, || {
            ref_Pitch_fr(
                st, mode, t_op, exc, xn, h, l_subfr, i_subfr, frac, resu3, idx,
            )
        })
    }
}

unsafe fn c_g_pitch(
    mode: i32,
    xn: *mut i16,
    y1: *mut i16,
    g: *mut i16,
    l: i16,
    o: *mut CFlag,
) -> i16 {
    unsafe { call(o, || ref_G_pitch(mode, xn, y1, g, l)) }
}

/// `Pitch_fr` inputs: excitation history (`exc[EXC]` is the subframe start).
const EXC: usize = PIT_MAX + L_INTERPOL;

pub(crate) struct PitchIn {
    pub(crate) mode: CodecMode,
    pub(crate) t_op: [i16; 2],
    pub(crate) t0_prev: i16,
    pub(crate) i_subfr: i16,
    pub(crate) exc: [i16; EXC + 40],
    pub(crate) xn: [i16; 40],
    pub(crate) h: [i16; 40],
}

pub(crate) fn pitch_in(g: &mut Rng) -> PitchIn {
    let mut p = PitchIn {
        mode: mode_of(g),
        t_op: [20 + (g.u64() % 124) as i16, 20 + (g.u64() % 124) as i16],
        t0_prev: 20 + (g.u64() % 124) as i16,
        i_subfr: 40 * (g.u64() % 4) as i16,
        exc: [0; EXC + 40],
        xn: [0; 40],
        h: [0; 40],
    };
    signal(g, &mut p.exc);
    signal(g, &mut p.xn);
    impulse(g, &mut p.h);
    p.h[0] = 4096 + (g.u64() % 8192) as i16;
    p
}

/// (lag, frac, resu3, index, state, overflow)
pub(crate) fn pitch_fr_rust(p: &PitchIn) -> (i16, i16, bool, i16, i16, bool) {
    let mut st = PitchFrState {
        t0_prev_subframe: p.t0_prev,
    };
    let mut ov = false;
    let (lag, frac, resu3, idx) = pitch_fr(
        &mut st, p.mode, &p.t_op, &p.exc, EXC, &p.xn, &p.h, p.i_subfr, &mut ov,
    );
    let _ = ov;
    (lag, frac, resu3, idx, st.t0_prev_subframe, false)
}

pub(crate) fn pitch_fr_c(p: &PitchIn) -> (i16, i16, bool, i16, i16, bool) {
    let (mut st, mut t_op, mut exc, mut xn, mut h) = (p.t0_prev, p.t_op, p.exc, p.xn, p.h);
    let (mut frac, mut resu3, mut idx, mut o) = (0, 0, 0, 0);
    let lag = unsafe {
        c_pitch_fr(
            &mut st,
            p.mode as i32,
            t_op.as_mut_ptr(),
            exc.as_mut_ptr().add(EXC),
            xn.as_mut_ptr(),
            h.as_mut_ptr(),
            40,
            p.i_subfr,
            &mut frac,
            &mut resu3,
            &mut idx,
            &mut o,
        )
    };
    let _ = o; // the reference's flag is global and sticky: values only
    (lag, frac, resu3 != 0, idx, st, false)
}

#[test]
fn convolve_matches_c() {
    let mut g = Rng::new(22);
    for _ in 0..N {
        let (mut x, mut h) = ([0i16; 40], [0i16; 40]);
        signal(&mut g, &mut x);
        impulse(&mut g, &mut h);
        let mut y = [0i16; 40];
        convolve(&x, &h, &mut y);
        let mut cy = [0i16; 40];
        unsafe { c_convolve(x.as_mut_ptr(), h.as_mut_ptr(), cy.as_mut_ptr(), 40) };
        assert_eq!(y, cy, "{x:?} {h:?}");
    }
}

#[test]
fn pitch_fr_matches_c() {
    let mut g = Rng::new(25);
    for _ in 0..N / 10 {
        let mut p = pitch_in(&mut g);
        if g.u64().is_multiple_of(4) {
            // Full-scale target: <xn, excf> saturates while <excf, excf> stays small
            // (reaches norm_corr's Cauchy-Schwarz bound).
            p.xn = [[32767, -32768][g.u64() as usize & 1]; 40];
        }
        assert_eq!(
            pitch_fr_rust(&p),
            pitch_fr_c(&p),
            "{:?} {:?}",
            p.mode,
            p.t_op
        );
    }
}
