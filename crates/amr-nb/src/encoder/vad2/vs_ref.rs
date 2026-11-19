//! Differential tests of `vad2` against the TS 26.073 reference (`ref26073-sys`;
//! helpers in `encoder::vs_ref`).
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::encoder::vad2::{Vad2State, ltp_flag_update, r_fft, vad2};
use crate::encoder::vs_ref::whole::Gen;
use crate::encoder::vs_ref::*;
use crate::testutil::with_ref;
use crate::testutil::{N, Rng};

#[test]
fn vad2_matches_c() {
    for seed in 0..20 {
        let mut sig = Gen::new(seed + 50);
        let mut g = Rng::new(seed);
        let (mut r, mut c) = (Vad2State::default(), CVad2::default());
        for f in 0..2000 {
            let frame = sig.frame();
            for half in frame.as_chunks::<80>().0 {
                let ltp = g.u64() & 1 == 1;
                r.ltp_flag = ltp;
                c.ltp_flag = ltp as i32;
                let mut h = *half;
                let (mut ov, mut o) = (false, 0);
                let rv = vad2(&h, &mut r, &mut ov);
                let cv = unsafe { c_vad2(h.as_mut_ptr(), &mut c, &mut o) };
                assert!(
                    rv as i16 == cv && same2(&r, &c),
                    "seed {seed} frame {f}: {r:?} vs {c:?}"
                );
            }
            let (rmax, r0) = (g.i32(), g.i32());
            (r.l_rmax, r.l_r0, c.l_rmax, c.l_r0) = (rmax, r0, rmax, r0);
            let mode = mode_of(&mut g);
            let (mut ov, mut o) = (false, 0);
            ltp_flag_update(&mut r, mode, &mut ov);
            unsafe { LTP_flag_update(&mut c, mode as i16, &mut o) };
            assert!(same2(&r, &c), "seed {seed} frame {f}: LTP_flag_update");
        }
    }
}

unsafe extern "C" {
    fn ref_vad2(farray: *mut i16, st: *mut CVad2) -> i16;
    fn ref_r_fft(farray: *mut i16);
    fn ref_LTP_flag_update(st: *mut CVad2, mode: i16);
}

#[repr(C)]
#[derive(Default, Debug, PartialEq)]
struct CVad2 {
    pub(crate) pre_emp_mem: i16,
    pub(crate) update_cnt: i16,
    pub(crate) hyster_cnt: i16,
    pub(crate) last_update_cnt: i16,
    pub(crate) ch_enrg_long_db: [i16; 16],
    pub(crate) l_frame_cnt: i32,
    pub(crate) l_ch_enrg: [i32; 16],
    pub(crate) l_ch_noise: [i32; 16],
    pub(crate) last_normb_shift: i16,
    pub(crate) tsnr: i16,
    pub(crate) hangover: i16,
    pub(crate) burstcount: i16,
    pub(crate) fupdate_flag: i16,
    pub(crate) neg_snr_var: i16,
    pub(crate) neg_snr_bias: i16,
    pub(crate) shift_state: i16,
    pub(crate) l_r0: i32,
    pub(crate) l_rmax: i32,
    pub(crate) ltp_flag: i32,
}

unsafe fn c_vad2(farray: *mut i16, st: *mut CVad2, o: *mut CFlag) -> i16 {
    with_ref(o, || unsafe { ref_vad2(farray, st) })
}

unsafe fn c_r_fft(farray: *mut i16, o: *mut CFlag) {
    with_ref(o, || unsafe { ref_r_fft(farray) })
}

#[allow(non_snake_case)]
unsafe fn LTP_flag_update(st: *mut CVad2, mode: i16, o: *mut CFlag) {
    with_ref(o, || unsafe { ref_LTP_flag_update(st, mode) })
}

fn same2(r: &Vad2State, c: &CVad2) -> bool {
    [r.pre_emp_mem, r.update_cnt, r.hyster_cnt, r.last_update_cnt]
        == [c.pre_emp_mem, c.update_cnt, c.hyster_cnt, c.last_update_cnt]
        && r.ch_enrg_long_db == c.ch_enrg_long_db
        && r.l_frame_cnt == c.l_frame_cnt
        && r.l_ch_enrg == c.l_ch_enrg
        && r.l_ch_noise == c.l_ch_noise
        && [
            r.last_normb_shift,
            r.tsnr,
            r.hangover,
            r.burstcount,
            r.fupdate_flag,
        ] == [
            c.last_normb_shift,
            c.tsnr,
            c.hangover,
            c.burstcount,
            c.fupdate_flag,
        ]
        && [r.neg_snr_var, r.neg_snr_bias, r.shift_state]
            == [c.neg_snr_var, c.neg_snr_bias, c.shift_state]
        && (r.l_r0, r.l_rmax, r.ltp_flag as i32) == (c.l_r0, c.l_rmax, c.ltp_flag)
}

#[test]
fn r_fft_matches_c() {
    let mut g = Rng::new(9);
    for _ in 0..N / 10 {
        let mut x = [0i16; 128];
        signal(&mut g, &mut x);
        let mut c = x;
        let mut ov = false;
        r_fft(&mut x, &mut ov);
        let mut o = 0;
        unsafe { c_r_fft(c.as_mut_ptr(), &mut o) };
        assert_eq!(x, c);
    }
}
