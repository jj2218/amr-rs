//! Differential tests of `vad` against the TS 26.073 reference (`ref26073-sys`;
//! helpers in `encoder::vs_ref`).
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::encoder::vad::{VadState, vad1};
use crate::encoder::vs_ref::whole::Gen;
use crate::encoder::vs_ref::*;
use crate::testutil::Rng;
use crate::testutil::with_ref;

#[test]
fn vad1_matches_c() {
    for seed in 0..20 {
        let mut sig = Gen::new(seed);
        let mut g = Rng::new(seed);
        let (mut r, mut c) = (VadState::default(), CVad1::default());
        unsafe { vad1_reset(&mut c) };
        let mut buf = [0i16; 200];
        for f in 0..1000 {
            buf.copy_within(160.., 0);
            buf[40..].copy_from_slice(&sig.frame());
            // The pitch/tone/complex inputs come from the pitch analysis; vary them too.
            let (tone, corr) = (g.i16(), g.i16());
            r.tone = tone;
            c.tone = tone;
            r.best_corr_hp = corr;
            c.best_corr_hp = corr;
            let (mut ov, mut o) = (false, 0);
            let rv = vad1(&mut r, &buf, &mut ov);
            let cv = unsafe { c_vad1(&mut c, buf.as_mut_ptr().add(40), &mut o) };
            assert!(
                rv == cv && same1(&r, &c),
                "seed {seed} frame {f}: {r:?} vs {c:?}"
            );
        }
    }
}

unsafe extern "C" {
    fn ref_vad1_reset(st: *mut CVad1) -> i32;
    fn ref_vad1(st: *mut CVad1, in_buf: *mut i16) -> i16;
}

#[repr(C)]
#[derive(Default, Debug, PartialEq)]
struct CVad1 {
    pub(crate) bckr_est: [i16; 9],
    pub(crate) ave_level: [i16; 9],
    pub(crate) old_level: [i16; 9],
    pub(crate) sub_level: [i16; 9],
    pub(crate) a_data5: [[i16; 2]; 3],
    pub(crate) a_data3: [i16; 5],
    pub(crate) burst_count: i16,
    pub(crate) hang_count: i16,
    pub(crate) stat_count: i16,
    pub(crate) vadreg: i16,
    pub(crate) pitch: i16,
    pub(crate) tone: i16,
    pub(crate) complex_high: i16,
    pub(crate) complex_low: i16,
    pub(crate) oldlag_count: i16,
    pub(crate) oldlag: i16,
    pub(crate) complex_hang_count: i16,
    pub(crate) complex_hang_timer: i16,
    pub(crate) best_corr_hp: i16,
    pub(crate) speech_vad_decision: i16,
    pub(crate) complex_warning: i16,
    pub(crate) sp_burst_count: i16,
    pub(crate) corr_hp_fast: i16,
}

unsafe fn vad1_reset(st: *mut CVad1) -> i16 {
    with_ref(std::ptr::null_mut(), || unsafe { ref_vad1_reset(st) }) as i16
}

unsafe fn c_vad1(st: *mut CVad1, in_buf: *mut i16, o: *mut CFlag) -> i16 {
    with_ref(o, || unsafe { ref_vad1(st, in_buf) })
}

fn same1(r: &VadState, c: &CVad1) -> bool {
    r.bckr_est == c.bckr_est
        && r.ave_level == c.ave_level
        && r.old_level == c.old_level
        && r.sub_level == c.sub_level
        && r.a_data5 == c.a_data5
        && r.a_data3 == c.a_data3
        && [
            r.burst_count,
            r.hang_count,
            r.stat_count,
            r.vadreg,
            r.pitch,
            r.tone,
        ] == [
            c.burst_count,
            c.hang_count,
            c.stat_count,
            c.vadreg,
            c.pitch,
            c.tone,
        ]
        && [
            r.complex_high,
            r.complex_low,
            r.oldlag_count,
            r.oldlag,
            r.complex_hang_count,
        ] == [
            c.complex_high,
            c.complex_low,
            c.oldlag_count,
            c.oldlag,
            c.complex_hang_count,
        ]
        && [r.complex_hang_timer, r.best_corr_hp, r.speech_vad_decision]
            == [c.complex_hang_timer, c.best_corr_hp, c.speech_vad_decision]
        && [r.complex_warning, r.sp_burst_count, r.corr_hp_fast]
            == [c.complex_warning, c.sp_burst_count, c.corr_hp_fast]
}
