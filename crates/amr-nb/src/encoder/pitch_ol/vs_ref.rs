//! Differential tests of `pitch_ol` against the TS 26.073 reference (`ref26073-sys`;
//! helpers in `encoder::vs_ref`).
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::cnst::{CodecMode, L_FRAME, PIT_MAX};
use crate::encoder::pitch_ol::comp_corr;
use crate::encoder::pitch_ol::pitch_ol;
use crate::encoder::vad::{Vad, VadState};
use crate::encoder::vs_ref::whole::Gen;
use crate::encoder::vs_ref::*;
use crate::testutil::with_ref;
use crate::testutil::{N, Rng};
use std::ffi::c_void;

/// Pitch_ol with and without DTX (the VAD1 hooks: tone detection, complex detection).
#[test]
fn pitch_ol_matches_c() {
    let mut g = Rng::new(6);
    let mut sig = Gen::new(6);
    let (mut rv, mut cv) = (Vad::One(VadState::default()), VadState::default());
    let mut buf = [0i16; PIT_MAX + L_FRAME];
    for _ in 0..N / 40 {
        buf.copy_within(L_FRAME.., 0);
        if g.u64().is_multiple_of(4) {
            signal(&mut g, &mut buf[PIT_MAX..]);
        } else {
            buf[PIT_MAX..].copy_from_slice(&sig.frame());
        }
        let mode = mode_of(&mut g);
        let dtx = g.u64() & 1 == 1;
        let idx = (g.u64() & 1) as i16;
        let (pit_min, l_frame) = match mode {
            CodecMode::MR475 | CodecMode::MR515 => (20, L_FRAME),
            CodecMode::MR122 => (18, L_FRAME / 2),
            _ => (20, L_FRAME / 2),
        };
        let (mut ov, mut o) = (false, 0);
        let r = pitch_ol(&mut rv, mode, &buf, pit_min, l_frame, idx, dtx, &mut ov);
        let c = unsafe {
            c_pitch_ol(
                (&raw mut cv).cast(),
                mode as i32,
                buf.as_mut_ptr().add(PIT_MAX),
                pit_min as i16,
                PIT_MAX as i16,
                l_frame as i16,
                idx,
                dtx as i32,
                &mut o,
            )
        };
        let Vad::One(rvs) = &rv else { unreachable!() };
        assert_eq!((r, rvs), (c, &cv), "{mode:?} dtx {dtx}");
    }
}

unsafe extern "C" {
    fn ref_Pitch_ol(
        vad: *mut c_void,
        mode: i32,
        sig: *mut i16,
        pit_min: i16,
        pit_max: i16,
        l_frame: i16,
        idx: i16,
        dtx: i32,
    ) -> i16;
    #[link_name = "comp_corr"]
    fn c_comp_corr(sig: *mut i16, l_frame: i16, lag_max: i16, lag_min: i16, corr: *mut i32);
}

#[allow(clippy::too_many_arguments)]
unsafe fn c_pitch_ol(
    vad: *mut c_void,
    mode: i32,
    sig: *mut i16,
    pit_min: i16,
    pit_max: i16,
    l_frame: i16,
    idx: i16,
    dtx: i32,
    o: *mut CFlag,
) -> i16 {
    with_ref(o, || unsafe {
        ref_Pitch_ol(vad, mode, sig, pit_min, pit_max, l_frame, idx, dtx)
    })
}

/// Pitch search signal: PIT_MAX history + up to one frame, scaled as `pitch_ol` does.
pub(crate) fn pitch_signal(g: &mut Rng) -> [i16; PIT_MAX + L_FRAME] {
    let mut sig = [0i16; PIT_MAX + L_FRAME];
    signal(g, &mut sig);
    sig
}

/// C `comp_corr` output (`out[m]` = lag `lag_max - m`) in Rust's `corr[lag]` layout.
pub(crate) fn c_comp_corr_at(
    sig: &mut [i16],
    l: usize,
    max: usize,
    min: usize,
) -> [i32; PIT_MAX + 1] {
    let mut out = [0i32; PIT_MAX + 4];
    let p = sig[PIT_MAX..].as_mut_ptr();
    unsafe {
        c_comp_corr(
            p,
            l as i16,
            max as i16,
            min as i16,
            out.as_mut_ptr().add(max),
        )
    };
    let mut corr = [0i32; PIT_MAX + 1];
    for m in 0..((max - min) & !3) + 4 {
        corr[max - m] = out[m];
    }
    corr
}

#[test]
fn comp_corr_matches_c() {
    let mut g = Rng::new(23);
    for _ in 0..N / 10 {
        let mut sig = pitch_signal(&mut g);
        let l = [L_FRAME, L_FRAME / 2][g.u64() as usize & 1];
        let min = [18, 20][g.u64() as usize & 1];
        let mut corr = [0i32; PIT_MAX + 1];
        comp_corr(&sig, l, PIT_MAX, min, &mut corr);
        assert_eq!(
            corr,
            c_comp_corr_at(&mut sig, l, PIT_MAX, min),
            "l {l} min {min}"
        );
    }
}
