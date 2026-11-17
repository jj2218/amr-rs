//! Differential tests of the encoder modules against the C oracle (opencore-amr 0.1.6).
#![allow(unsafe_code)]

use crate::basic_op::Flag;
use crate::cnst::{L_FRAME, L_SUBFR, M, Mode, PIT_MAX};
use crate::testutil::{N, Rng};
use opencore_sys as _;

type CFlag = i32;

/// A signal of random amplitude (`>> 0..16`), sprinkled with edge values.
fn signal(g: &mut Rng, buf: &mut [i16]) {
    let scale = (g.u64() % 16) as u32;
    for v in buf {
        *v = if g.u64().is_multiple_of(32) {
            g.i16()
        } else {
            (g.u64() as i16) >> scale
        };
    }
}

pub(super) fn mode_of(g: &mut Rng) -> Mode {
    Mode::ALL[(g.u64() % 8) as usize]
}

/// Against the TS 26.073 reference (opencore's rewrites of both wrap on full-scale input, and
/// its `G_code` shifts by >= 32, which is UB).
mod gains {
    use super::*;
    use crate::encoder::gain_q::{g_code, q_gain_code};
    use ref26073_sys::lock;

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
    }

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

    #[test]
    fn q_gain_code_matches_ref() {
        let _l = lock();
        let mut g = Rng::new(2);
        for _ in 0..N {
            let mode = if g.u64() & 1 == 0 {
                Mode::MR122
            } else {
                Mode::MR795
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
                "{mode} exp {exp} frac {frac} gain {gain0}"
            );
        }
    }
}

const _: usize = M;

/// The TS 26.073 reference (`ref26073-sys`) behind the opencore-shaped signatures these tests
/// were written against: the tables argument is unused, `o` receives the reference's global
/// overflow flag (cleared before the call). Each call holds the reference lock.
#[allow(clippy::too_many_arguments)]
mod refc {
    use super::CFlag;
    use ref26073_sys::{clear_overflow, lock, overflow};

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
        fn ref_cor_h_x(h: *mut i16, x: *mut i16, dn: *mut i16, sf: i16);
        fn ref_cor_h_x2(h: *mut i16, x: *mut i16, dn: *mut i16, sf: i16, nb: i16, step: i16);
        fn ref_cor_h(h: *mut i16, sign: *mut i16, rr: *mut [i16; 40]);
        fn ref_set_sign12k2(
            dn: *mut i16,
            cn: *mut i16,
            sign: *mut i16,
            pos_max: *mut i16,
            nb_track: i16,
            ipos: *mut i16,
            step: i16,
        );
        fn ref_Autocorr(x: *mut i16, m: i16, rh: *mut i16, rl: *mut i16, w: *const i16) -> i16;
        fn ref_Lag_window(m: i16, rh: *mut i16, rl: *mut i16);
        fn ref_Levinson(st: *mut i16, rh: *mut i16, rl: *mut i16, a: *mut i16, rc: *mut i16)
        -> i16;
        fn ref_Pre_Process(st: *mut i16, signal: *mut i16, lg: i16) -> i32;
        fn ref_search_10and8i40(
            nb_pulse: i16,
            step: i16,
            nb_tracks: i16,
            dn: *mut i16,
            rr: *mut [i16; 40],
            ipos: *mut i16,
            pos_max: *mut i16,
            codvec: *mut i16,
        );
        fn ref_G_pitch(mode: i32, xn: *mut i16, y1: *mut i16, g: *mut i16, l: i16) -> i16;
        fn ref_cbsearch(
            x: *mut i16,
            h: *mut i16,
            t0: i16,
            sharp: i16,
            gain_pit: i16,
            res2: *mut i16,
            code: *mut i16,
            y: *mut i16,
            anap: *mut *mut i16,
            mode: i32,
            sub_nr: i16,
        ) -> i32;
    }

    /// Runs `f` on the reference with a cleared flag, storing the flag in `o` (if not null).
    unsafe fn call<T>(o: *mut CFlag, f: impl FnOnce() -> T) -> T {
        let _l = lock();
        clear_overflow();
        let r = f();
        if !o.is_null() {
            unsafe { *o = overflow() as CFlag };
        }
        r
    }

    pub(super) unsafe fn c_convolve(x: *mut i16, h: *mut i16, y: *mut i16, l: i16) {
        unsafe { call(std::ptr::null_mut(), || ref_Convolve(x, h, y, l)) }
    }
    pub(super) unsafe fn c_pitch_fr(
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
    pub(super) unsafe fn c_qua_gain(
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
        _t: *mut super::modules::CommonAmrTbls,
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
    pub(super) unsafe fn c_mr475_gain_quant(
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
                    pred, e0, f0, ec0, fc0, te0, tf0, code, e1, f1, ec1, fc1, te1, tf1, gp_limit,
                    p0, c0, p1, c1,
                )
            })
        }
    }
    pub(super) unsafe fn c_calc_filt_energies(
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
    pub(super) unsafe fn c_cor_h_x(h: *mut i16, x: *mut i16, dn: *mut i16, sf: i16, o: *mut CFlag) {
        unsafe { call(o, || ref_cor_h_x(h, x, dn, sf)) }
    }
    pub(super) unsafe fn c_cor_h_x2(
        h: *mut i16,
        x: *mut i16,
        dn: *mut i16,
        sf: i16,
        nb: i16,
        step: i16,
        o: *mut CFlag,
    ) {
        unsafe { call(o, || ref_cor_h_x2(h, x, dn, sf, nb, step)) }
    }
    pub(super) unsafe fn ref_search_10and8i40_w(
        nb_pulse: i16,
        step: i16,
        nb_tracks: i16,
        dn: *mut i16,
        rr: *mut [i16; 40],
        ipos: *mut i16,
        pos_max: *mut i16,
        codvec: *mut i16,
        o: *mut CFlag,
    ) {
        unsafe {
            call(o, || {
                ref_search_10and8i40(nb_pulse, step, nb_tracks, dn, rr, ipos, pos_max, codvec)
            })
        }
    }
    pub(super) unsafe fn c_autocorr(
        x: *mut i16,
        m: i16,
        rh: *mut i16,
        rl: *mut i16,
        w: *const i16,
        o: *mut CFlag,
    ) -> i16 {
        unsafe { call(o, || ref_Autocorr(x, m, rh, rl, w)) }
    }
    pub(super) unsafe fn c_lag_window(m: i16, rh: *mut i16, rl: *mut i16, o: *mut CFlag) {
        unsafe { call(o, || ref_Lag_window(m, rh, rl)) }
    }
    pub(super) unsafe fn c_levinson<S>(
        st: *mut S,
        rh: *mut i16,
        rl: *mut i16,
        a: *mut i16,
        rc: *mut i16,
        o: *mut CFlag,
    ) -> i16 {
        unsafe { call(o, || ref_Levinson(st.cast(), rh, rl, a, rc)) }
    }
    /// `Pre_Process` (its state: 6 Word16s, as `Pre_ProcessState`).
    pub(super) unsafe fn c_pre_process<S>(st: *mut S, signal: *mut i16, lg: i16) {
        unsafe {
            call(std::ptr::null_mut(), || {
                ref_Pre_Process(st.cast(), signal, lg)
            })
        };
    }
    pub(super) unsafe fn c_set_sign12k2(
        dn: *mut i16,
        cn: *mut i16,
        sign: *mut i16,
        pos_max: *mut i16,
        nb_track: i16,
        ipos: *mut i16,
        step: i16,
        o: *mut CFlag,
    ) {
        unsafe {
            call(o, || {
                ref_set_sign12k2(dn, cn, sign, pos_max, nb_track, ipos, step)
            })
        }
    }
    pub(super) unsafe fn c_cor_h(h: *mut i16, sign: *mut i16, rr: *mut [i16; 40], o: *mut CFlag) {
        unsafe { call(o, || ref_cor_h(h, sign, rr)) }
    }
    pub(super) unsafe fn c_g_pitch(
        mode: i32,
        xn: *mut i16,
        y1: *mut i16,
        g: *mut i16,
        l: i16,
        o: *mut CFlag,
    ) -> i16 {
        unsafe { call(o, || ref_G_pitch(mode, xn, y1, g, l)) }
    }
    pub(super) unsafe fn c_cbsearch(
        x: *mut i16,
        h: *mut i16,
        t0: i16,
        sharp: i16,
        gain_pit: i16,
        res2: *mut i16,
        code: *mut i16,
        y: *mut i16,
        anap: *mut *mut i16,
        mode: i32,
        sub_nr: i16,
        _t: *mut super::modules::CommonAmrTbls,
        o: *mut CFlag,
    ) {
        unsafe {
            call(o, || {
                ref_cbsearch(x, h, t0, sharp, gain_pit, res2, code, y, anap, mode, sub_nr);
            })
        }
    }
}

/// Whole-encoder differential: Rust `Encoder` vs the TS 26.073 reference encoder
/// (`Speech_Encode_Frame` + `sid_sync`, as its `coder.c` drives them) on synthetic signals that
/// reach the saturation and DTX paths the 3GPP vectors don't. (opencore wraps on many of them.)
mod whole {
    use super::*;
    use crate::cnst::L_FRAME;
    use crate::decoder::MAX_SERIAL_SIZE;
    use crate::encoder::{Dtx, Encoder};
    use std::ffi::c_void;
    use std::ptr::null_mut;

    unsafe extern "C" {
        fn ref_Speech_Encode_Frame_init(st: *mut *mut c_void, dtx: i32, id: *const u8) -> i32;
        fn ref_Speech_Encode_Frame_exit(st: *mut *mut c_void);
        fn ref_Speech_Encode_Frame(
            st: *mut c_void,
            mode: i32,
            speech: *mut i16,
            serial: *mut i16,
            used_mode: *mut i32,
        ) -> i32;
        fn ref_sid_sync_init(st: *mut *mut c_void) -> i32;
        fn ref_sid_sync_exit(st: *mut *mut c_void);
        fn ref_sid_sync(st: *mut c_void, mode: i32, tx_type: *mut i32);
    }

    /// Signal generator: per segment a random mix of tone, noise, silence and clipping.
    pub(crate) struct Gen {
        g: Rng,
        phase: f64,
        freq: f64,
        amp: f64,
        noise: i32,
        clip: bool,
        left: usize,
    }

    impl Gen {
        pub(crate) fn new(seed: u64) -> Self {
            Self {
                g: Rng::new(seed),
                phase: 0.0,
                freq: 0.0,
                amp: 0.0,
                noise: 0,
                clip: false,
                left: 0,
            }
        }
        pub(crate) fn frame(&mut self) -> [i16; L_FRAME] {
            if self.left == 0 {
                let g = &mut self.g;
                self.left = 1 + (g.u64() % 40) as usize;
                self.freq = (g.u64() % 3800) as f64 / 8000.0;
                self.amp = [0.0, 30.0, 500.0, 8000.0, 40000.0][(g.u64() % 5) as usize];
                self.noise = [0, 4, 200, 3000, 32767][(g.u64() % 5) as usize];
                self.clip = g.u64().is_multiple_of(8);
            }
            self.left -= 1;
            std::array::from_fn(|_| {
                self.phase += self.freq * std::f64::consts::TAU;
                let n = if self.noise > 0 {
                    (self.g.u64() % (2 * self.noise as u64 + 1)) as i32 - self.noise
                } else {
                    0
                };
                let v = (self.amp * self.phase.sin()) as i32 + n;
                let v = if self.clip { v.signum() * 32767 } else { v };
                v.clamp(-32768, 32767) as i16
            })
        }
    }

    fn run(seed: u64, dtx: bool, frames: usize) {
        let mut sig = Gen::new(seed);
        let mut g = Rng::new(seed ^ 0xabcdef);
        let _l = ref26073_sys::lock();
        let mut rust = Encoder::new(if dtx { Dtx::Vad1 } else { Dtx::Off });
        let (mut enc, mut sid) = (null_mut(), null_mut());
        unsafe {
            ref_Speech_Encode_Frame_init(&mut enc, dtx as i32, c"encoder".as_ptr().cast());
            ref_sid_sync_init(&mut sid);
        }
        let mut mode = mode_of(&mut g);
        for f in 0..frames {
            if g.u64().is_multiple_of(16) {
                mode = mode_of(&mut g);
            }
            let pcm = sig.frame();
            assert!(
                pcm.iter().any(|&s| s != 0x0008),
                "a homing frame: coder.c would reset"
            );
            let mut speech = pcm; // C overwrites its input
            let mut ets = [0i16; MAX_SERIAL_SIZE + 2];
            let (mut used, mut tx) = (0, 0);
            unsafe {
                ref_Speech_Encode_Frame(
                    enc,
                    mode as i32,
                    speech.as_mut_ptr(),
                    ets[1..].as_mut_ptr(),
                    &mut used,
                );
                ref_sid_sync(sid, used, &mut tx);
            }
            ets[0] = tx as i16;
            let r = rust.encode(mode, &pcm);
            let tx = match r.frame_type {
                crate::RxFrameType::SpeechGood => 0,
                crate::RxFrameType::SidFirst => 1,
                crate::RxFrameType::SidUpdate => 2,
                _ => 3,
            };
            assert_eq!(
                tx, ets[0],
                "seed {seed} dtx {dtx} frame {f} {mode}: tx type"
            );
            assert_eq!(
                r.bits[..],
                ets[1..=MAX_SERIAL_SIZE],
                "seed {seed} dtx {dtx} frame {f} {mode}: bits"
            );
        }
        unsafe {
            ref_Speech_Encode_Frame_exit(&mut enc);
            ref_sid_sync_exit(&mut sid);
        }
    }

    #[test]
    fn encoder_matches_ref_without_dtx() {
        for seed in 0..40 {
            run(seed, false, 500);
        }
    }

    #[test]
    fn encoder_matches_ref_with_dtx() {
        for seed in 100..140 {
            run(seed, true, 500);
        }
    }
}

/// VAD option 1 and 2, the FFT and `LTP_flag_update`: long call sequences, full state compared
/// after every call.
mod vads {
    use super::whole::Gen;
    use super::*;
    use crate::encoder::vad::{VadState, vad1};
    use crate::encoder::vad2::{Vad2State, ltp_flag_update, r_fft, vad2};

    #[repr(C)]
    #[derive(Default, Debug, PartialEq)]
    struct CVad1 {
        bckr_est: [i16; 9],
        ave_level: [i16; 9],
        old_level: [i16; 9],
        sub_level: [i16; 9],
        a_data5: [[i16; 2]; 3],
        a_data3: [i16; 5],
        burst_count: i16,
        hang_count: i16,
        stat_count: i16,
        vadreg: i16,
        pitch: i16,
        tone: i16,
        complex_high: i16,
        complex_low: i16,
        oldlag_count: i16,
        oldlag: i16,
        complex_hang_count: i16,
        complex_hang_timer: i16,
        best_corr_hp: i16,
        speech_vad_decision: i16,
        complex_warning: i16,
        sp_burst_count: i16,
        corr_hp_fast: i16,
    }

    #[repr(C)]
    #[derive(Default, Debug, PartialEq)]
    struct CVad2 {
        pre_emp_mem: i16,
        update_cnt: i16,
        hyster_cnt: i16,
        last_update_cnt: i16,
        ch_enrg_long_db: [i16; 16],
        l_frame_cnt: i32,
        l_ch_enrg: [i32; 16],
        l_ch_noise: [i32; 16],
        last_normb_shift: i16,
        tsnr: i16,
        hangover: i16,
        burstcount: i16,
        fupdate_flag: i16,
        neg_snr_var: i16,
        neg_snr_bias: i16,
        shift_state: i16,
        l_r0: i32,
        l_rmax: i32,
        ltp_flag: i32,
    }

    // The TS 26.073 reference (same state layouts), in the shape of the opencore functions
    // these tests were written against (`o`: its overflow flag, sticky).
    unsafe extern "C" {
        fn ref_vad1_reset(st: *mut CVad1) -> i32;
        fn ref_vad1(st: *mut CVad1, in_buf: *mut i16) -> i16;
        fn ref_vad2(farray: *mut i16, st: *mut CVad2) -> i16;
        fn ref_r_fft(farray: *mut i16);
        fn ref_LTP_flag_update(st: *mut CVad2, mode: i16);
    }
    use crate::testutil::with_ref;
    unsafe fn vad1_reset(st: *mut CVad1) -> i16 {
        with_ref(std::ptr::null_mut(), || unsafe { ref_vad1_reset(st) }) as i16
    }
    unsafe fn c_vad1(st: *mut CVad1, in_buf: *mut i16, o: *mut CFlag) -> i16 {
        with_ref(o, || unsafe { ref_vad1(st, in_buf) })
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
}

/// Stateless and short-state modules on random (often out-of-range) inputs.
mod modules {
    use super::refc::{c_autocorr, c_cbsearch, c_g_pitch, c_lag_window, c_levinson};

    // The TS 26.073 reference, in the shape of the opencore functions these tests were written
    // against (`o`: its overflow flag, sticky).
    unsafe extern "C" {
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
    }
    use crate::testutil::with_ref;
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

    unsafe extern "C" {
        #[link_name = "cbsearch"]
        fn oc_cbsearch(
            x: *mut i16,
            h: *mut i16,
            t0: i16,
            sharp: i16,
            gain_pit: i16,
            res2: *mut i16,
            code: *mut i16,
            y: *mut i16,
            anap: *mut *mut i16,
            mode: i32,
            sub_nr: i16,
            t: *mut CommonAmrTbls,
            o: *mut CFlag,
        );
    }
    use super::whole::Gen;
    use super::*;
    use crate::cnst::{L_FRAME, MP1, PIT_MAX};
    use crate::encoder::codebook::cbsearch;
    use crate::encoder::lpc::{
        L_WINDOW, LevinsonState, QPlsfState, autocorr, az_lsp, lag_window, levinson, q_plsf_3,
        q_plsf_5,
    };
    use crate::encoder::pitch_cl::g_pitch;
    use crate::encoder::pitch_ol::pitch_ol;
    use crate::encoder::vad::{Vad, VadState};
    use crate::tables::{WINDOW_160_80, WINDOW_200_40, WINDOW_232_8};
    use std::ffi::c_void;

    #[repr(C)]
    pub(super) struct CommonAmrTbls {
        ptrs: [*const c_void; 34],
    }

    unsafe extern "C" {
        fn get_const_tbls(t: *mut c_void);
        #[link_name = "Chebps_Wrapper"]
        fn c_chebps(x: i16, f: *mut i16, n: i16, o: *mut CFlag) -> i16;
    }

    pub(super) fn tables() -> CommonAmrTbls {
        let mut t = CommonAmrTbls {
            ptrs: [std::ptr::null(); 34],
        };
        unsafe { get_const_tbls((&raw mut t).cast()) };
        t
    }

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
                    let signs: [i16; L_WINDOW] =
                        std::array::from_fn(|_| (g.u64() & 1) as i16 * 2 - 1);
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
            let (mut rh, mut rl, mut ch, mut cl) =
                ([0i16; MP1], [0i16; MP1], [0i16; MP1], [0i16; MP1]);
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
            let (mut rh, mut rl, mut ch, mut cl) =
                ([0i16; MP1], [0i16; MP1], [0i16; MP1], [0i16; MP1]);
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
            let (mut ra, mut ca, mut rrc, mut crc) =
                ([0i16; MP1], [0i16; MP1], [0i16; 4], [0i16; 4]);
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
                    f[i + 1] =
                        (((a[i + 1] as i32 + a[M - i] as i32) >> 2) as i16).wrapping_sub(f[i]);
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

    /// Random descending LSP vectors (cosine domain, as Az_lsp produces them).
    fn lsp(g: &mut Rng) -> [i16; M] {
        let mut v: [i16; M] = std::array::from_fn(|_| ((g.u64() % 65535) as i32 - 32767) as i16);
        v.sort_unstable_by(|a, b| b.cmp(a));
        v
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
                let (mut rq1, mut rq2, mut cq1, mut cq2) =
                    ([0i16; M], [0i16; M], [0i16; M], [0i16; M]);
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
                let mode = Mode::ALL[(g.u64() % 9) as usize];
                if mode == Mode::MR122 {
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
                    "Q_plsf_3 {mode} {l1:?}"
                );
            }
        }
    }

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
                Mode::MR475 | Mode::MR515 => (20, L_FRAME),
                Mode::MR122 => (18, L_FRAME / 2),
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
            assert_eq!((r, rvs), (c, &cv), "{mode} dtx {dtx}");
        }
    }

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
            assert_eq!((r, rg), (c, cg), "{mode} {xn:?} {y1:?}");
        }
    }

    /// `cbsearch` inputs.
    pub(super) struct CbIn {
        pub(super) mode: Mode,
        x: [i16; L_SUBFR],
        res2: [i16; L_SUBFR],
        hvec: [i16; 2 * L_SUBFR],
        t0: i16,
        sharp: i16,
        gain_pit: i16,
        sub_nr: usize,
    }

    /// (index words written, indices, code, filtered code, sharpened h)
    type CbOut = (
        usize,
        [i16; 10],
        [i16; L_SUBFR],
        [i16; L_SUBFR],
        [i16; 2 * L_SUBFR],
    );

    pub(super) fn cb_in(g: &mut Rng, mode: Mode) -> CbIn {
        let (mut x, mut res2) = ([0i16; L_SUBFR], [0i16; L_SUBFR]);
        signal(g, &mut x);
        signal(g, &mut res2);
        let mut hvec = [0i16; 2 * L_SUBFR];
        // An impulse response: decaying, first tap large.
        let scale = (g.u64() % 4) as u32;
        for (i, v) in hvec[L_SUBFR..].iter_mut().enumerate() {
            *v = (((g.u64() as i16) >> (scale + (i as u32 / 6).min(12))) as i32)
                .clamp(-32768, 32767) as i16;
        }
        hvec[L_SUBFR] = 4096 + (g.u64() % 8192) as i16;
        CbIn {
            mode,
            x,
            res2,
            hvec,
            t0: 18 + (g.u64() % 126) as i16,
            sharp: (g.u64() % 13018) as i16,
            gain_pit: (g.u64() % 19662) as i16,
            sub_nr: (g.u64() % 4) as usize,
        }
    }

    pub(super) fn cb_rust(c: &CbIn) -> CbOut {
        let mut rh = c.hvec;
        let (mut rc, mut ry, mut ra) = ([0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; 10]);
        let (mut ov, mut pos) = (false, 0);
        cbsearch(
            &c.x,
            &mut rh[L_SUBFR..],
            c.t0,
            c.sharp,
            c.gain_pit,
            &c.res2,
            &mut rc,
            &mut ry,
            &mut ra,
            &mut pos,
            c.mode,
            c.sub_nr,
            &mut ov,
        );
        (pos, ra, rc, ry, rh)
    }

    pub(super) fn cb_c(c: &CbIn, t: &mut CommonAmrTbls) -> CbOut {
        let (mut x, mut res2, mut ch) = (c.x, c.res2, c.hvec);
        let (mut cc, mut cy, mut ca) = ([0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; 10]);
        let mut anap = ca.as_mut_ptr();
        let mut o = 0;
        unsafe {
            c_cbsearch(
                x.as_mut_ptr(),
                ch.as_mut_ptr().add(L_SUBFR),
                c.t0,
                c.sharp,
                c.gain_pit,
                res2.as_mut_ptr(),
                cc.as_mut_ptr(),
                cy.as_mut_ptr(),
                &mut anap,
                c.mode as i32,
                c.sub_nr as i16,
                t,
                &mut o,
            )
        };
        let n = unsafe { anap.offset_from(ca.as_ptr()) } as usize;
        (n, ca, cc, cy, ch)
    }

    /// `cb_c` on opencore (for `kernel_timing`).
    pub(super) fn cb_oc(c: &CbIn, t: &mut CommonAmrTbls) -> CbOut {
        let (mut x, mut res2, mut ch) = (c.x, c.res2, c.hvec);
        let (mut cc, mut cy, mut ca) = ([0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; 10]);
        let mut anap = ca.as_mut_ptr();
        let mut o = 0;
        unsafe {
            oc_cbsearch(
                x.as_mut_ptr(),
                ch.as_mut_ptr().add(L_SUBFR),
                c.t0,
                c.sharp,
                c.gain_pit,
                res2.as_mut_ptr(),
                cc.as_mut_ptr(),
                cy.as_mut_ptr(),
                &mut anap,
                c.mode as i32,
                c.sub_nr as i16,
                t,
                &mut o,
            )
        };
        let n = unsafe { anap.offset_from(ca.as_ptr()) } as usize;
        (n, ca, cc, cy, ch)
    }

    /// cbsearch for every mode: indices, code vector, filtered code vector, sharpened h.
    #[test]
    fn cbsearch_matches_c() {
        let mut g = Rng::new(8);
        let mut t = tables();
        for _ in 0..N / 20 {
            let mode = mode_of(&mut g);
            let c = cb_in(&mut g, mode);
            assert_eq!(
                cb_rust(&c),
                cb_c(&c, &mut t),
                "{mode} t0 {} sharp {}",
                c.t0,
                c.sharp
            );
        }
    }
}

/// Hot kernels, tested directly (they are rewritten for speed; the module tests above only
/// reach them through realistic callers).
mod kernels {
    use super::modules::cb_oc;
    use super::modules::{CbIn, CommonAmrTbls, cb_in, cb_rust, tables};
    use super::refc::*;
    use super::*;

    // opencore, for `kernel_timing` (`BENCH=kern`): the speed baseline stays the C we beat.
    unsafe extern "C" {
        #[link_name = "Convolve"]
        fn oc_convolve(x: *mut i16, h: *mut i16, y: *mut i16, l: i16);
        #[link_name = "cor_h"]
        fn oc_cor_h(h: *mut i16, sign: *mut i16, rr: *mut [i16; 40], o: *mut CFlag);
        #[link_name = "MR795_gain_quant"]
        fn oc_mr795_gain_quant(
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
            t: *mut CommonAmrTbls,
            o: *mut CFlag,
        );
    }
    use crate::cnst::L_INTERPOL;
    use crate::encoder::codebook::{
        Pair, RrTm, RrTm4, Tm, chain, cor_h_search_8i40, cor_h_tm, cor_h_tm_portable, cor_h_x,
        filtered_code, search_8i40_portable, search_10i40, search_10i40_portable, set_sign,
        set_sign12k2,
    };
    use crate::encoder::gain_q::{
        GainAdaptState, calc_filt_energies, mr475_gain_quant, mr795_gain_quant, qua_gain,
    };
    use crate::encoder::layout::cols;
    use crate::encoder::lpc::{
        DICO1_SOA, DICO2_HALF_SOA, DICO2_SOA, DICO3_SOA, MR515_3_SOA, MR795_1_SOA, vq_subvec3,
        vq_subvec4,
    };
    use crate::encoder::pitch_cl::{PitchFrState, convolve, pitch_fr};
    use crate::encoder::pitch_ol::comp_corr;
    use crate::gains::GcPredState;
    use crate::tables::{DICO1_LSF_3, DICO2_LSF_3, DICO3_LSF_3, MR515_3_LSF, MR795_1_LSF};

    unsafe extern "C" {
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
        #[link_name = "search_10and8i40"]
        fn c_search_10and8i40(
            nb_pulse: i16,
            step: i16,
            nb_tracks: i16,
            dn: *mut i16,
            rr: *mut [i16; 40],
            ipos: *mut i16,
            pos_max: *mut i16,
            codvec: *mut i16,
            o: *mut CFlag,
        );
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
        #[link_name = "Syn_filt"]
        fn c_syn_filt(a: *mut i16, x: *mut i16, y: *mut i16, lg: i16, mem: *mut i16, update: i16);
        #[link_name = "comp_corr"]
        fn c_comp_corr(sig: *mut i16, l_frame: i16, lag_max: i16, lag_min: i16, corr: *mut i32);
    }

    /// An impulse response: random, decaying at a random rate, sometimes full scale.
    fn impulse(g: &mut Rng, h: &mut [i16]) {
        let scale = (g.u64() % 6) as u32;
        let full = g.u64().is_multiple_of(8);
        for (i, v) in h.iter_mut().enumerate() {
            let shift = if full {
                0
            } else {
                scale + (i as u32 / 6).min(12)
            };
            *v = (g.u64() as i16) >> shift;
        }
    }

    /// `rr` in a track-major layout (padding lanes zero).
    fn tm<const T: usize, const N: usize, const W: usize>(rr: &[[i16; 40]; 40]) -> Tm<T, N, W> {
        let mut t = Tm::ZERO;
        for i in 0..40 {
            for j in 0..40 {
                t.b[i % T][j % T][i / T][j / T] = rr[i][j];
            }
            t.diag[i % T][i / T] = rr[i][i];
        }
        t
    }

    /// MR102's `cor_h` (both instances of `cor_h_search_8i40`: 10 lanes and the ±32767 form,
    /// 16 lanes and the multiply form).
    #[test]
    fn cor_h_tm4_matches_c() {
        let mut g = Rng::new(21);
        for _ in 0..N {
            let mut h = [0i16; 40];
            impulse(&mut g, &mut h);
            let mut sign = [0i16; 40];
            signs(&mut g, &mut sign);
            let (mut crr, mut o) = ([[0i16; 40]; 40], 0);
            unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), crr.as_mut_ptr(), &mut o) };
            let (mut r10, mut r16, mut ov10, mut ov16) = (RrTm4::ZERO, RrTm4::ZERO, false, false);
            cor_h_tm_portable::<4, 10, 10, true>(&h, &sign, &mut r10, &mut ov10);
            cor_h_tm_portable::<4, 10, 16, false>(&h, &sign, &mut r16, &mut ov16);
            let (c10, c16) = (tm::<4, 10, 10>(&crr), tm::<4, 10, 16>(&crr));
            assert_eq!(
                (r10.b, r10.diag, ov10 && o != 0),
                (c10.b, c10.diag, ov10 && o != 0),
                "{h:?} {sign:?}"
            );
            assert_eq!(
                (r16.b, r16.diag, ov16 && o != 0),
                (c16.b, c16.diag, ov16 && o != 0),
                "{h:?} {sign:?}"
            );
        }
    }

    /// Pulse signs: ±32767 as `set_sign12k2` makes them, sometimes arbitrary.
    fn signs(g: &mut Rng, sign: &mut [i16; 40]) {
        let any = g.u64().is_multiple_of(4);
        for s in sign {
            *s = match (any, g.u64() & 1) {
                (true, _) => g.i16(),
                (false, 0) => 32767,
                _ => -32767,
            };
        }
    }

    #[test]
    fn cor_h_tm_matches_c() {
        let mut g = Rng::new(26);
        for _ in 0..N / 4 {
            let mut h = [0i16; 40];
            impulse(&mut g, &mut h);
            let mut sign = [0i16; 40];
            signs(&mut g, &mut sign);
            let (mut rr, mut ov) = (RrTm::ZERO, false);
            cor_h_tm(&h, &sign, &mut rr, &mut ov);
            let (mut crr, mut o) = ([[0i16; 40]; 40], 0);
            unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), crr.as_mut_ptr(), &mut o) };
            let c = tm(&crr);
            let _ = (ov, o);
            assert_eq!((rr.b, rr.diag), (c.b, c.diag), "{h:?} {sign:?}");
        }
    }

    /// `cor_h_blocks` (dispatched) vs the C's `cor_h`: two random blocks (any tracks, equal ones
    /// too; half the cases with equal residues, one shared accumulation), with every diagonal;
    /// the other blocks stay untouched (a 0x5555 fill). Signs ±32767 (as `set_sign`'s) and of
    /// any value.
    #[test]
    fn cor_h_blocks_matches_c() {
        let mut g = Rng::new(27);
        for _ in 0..N / 4 {
            let mut h = [0i16; 40];
            impulse(&mut g, &mut h);
            let mut sign = [0i16; 40];
            signs(&mut g, &mut sign);
            let mut blocks: [[usize; 2]; 2] =
                core::array::from_fn(|_| [(g.u64() % 5) as usize, (g.u64() % 5) as usize]);
            if g.u64() & 1 == 0 {
                // equal residues (one shared accumulation), as in half the codec's calls
                let r = (blocks[0][1] + 5 - blocks[0][0]) % 5;
                blocks[1][1] = (blocks[1][0] + r) % 5;
            }
            let mut rr = RrTm::ZERO;
            rr.b = [[[[0x5555; 8]; 8]; 5]; 5];
            let mut ov = false;
            crate::simd::cor_h_blocks(&h, &sign, &mut rr, &blocks, &mut ov);
            let (mut crr, mut o) = ([[0i16; 40]; 40], 0);
            unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), crr.as_mut_ptr(), &mut o) };
            let c = tm(&crr);
            let _ = (ov, o);
            assert_eq!(rr.diag, c.diag, "{h:?} {sign:?}");
            for ta in 0..5 {
                for tb in 0..5 {
                    let want = if blocks.contains(&[ta, tb]) || ta == tb {
                        c.b[ta][tb]
                    } else {
                        [[0x5555; 8]; 8]
                    };
                    if blocks.contains(&[ta, tb]) || ta != tb {
                        assert_eq!(rr.b[ta][tb], want, "block {ta} {tb}: {h:?} {sign:?}");
                    } else {
                        // a diagonal block not asked for: only its diagonal is written
                        for n in 0..8 {
                            assert_eq!(rr.b[ta][ta][n][n], c.b[ta][ta][n][n]);
                        }
                    }
                }
            }
        }
    }

    /// `pm_scale` is cor_h's sign scaling for ±32767 signs, for every `tmp1`.
    #[test]
    fn pm_scale_is_exact() {
        use crate::encoder::codebook::pm_scale;
        for t in i16::MIN..=i16::MAX {
            for (si, sj) in [
                (32767, 32767),
                (-32767, -32767),
                (32767, -32767),
                (-32767, 32767),
            ] {
                let tmp2 = ((si * sj) >> 15) as i16;
                let want = ((t as i32 * tmp2 as i32) >> 15) as i16;
                assert_eq!(pm_scale(t, ((si ^ sj) >> 15) as i16), want, "{t} {si} {sj}");
            }
        }
    }

    /// `search_10and8i40`'s inputs: `dn`, square `rr`, `ipos`, `pos_max`.
    type SearchIn<const T: usize, const P: usize> =
        ([i16; 40], [[i16; 40]; 40], [i16; P], [i16; T]);

    /// Search inputs for `T` tracks (MR122: 5, MR102: 4): `rr` from `cor_h` with the encoder's
    /// kind of `h` and signs (the timing class: `Pair::best` falls back on <1% of picks, 0 on the
    /// vectors) or, if `odd`, of any `h` and signs (non-positive energies), sometimes a random
    /// symmetric matrix (energies of any sign: the exact-chain fallback), or a slowly decaying `h`
    /// with one sign for all pulses and small `dn` (large positive `rr` everywhere: the pair
    /// chains saturate while `dn` keeps the fast path, so `pair_safe` must decline), and `dn` from
    /// a few values whose sums tie and wrap to -32768; tracks and maxima as `set_sign12k2` leaves them.
    fn search_in<const T: usize, const P: usize>(g: &mut Rng, odd: bool) -> SearchIn<T, P> {
        let mut rr = [[0i16; 40]; 40];
        let flat = odd && g.u64().is_multiple_of(4);
        if flat {
            let mut h = [0i16; 40];
            let (top, step) = (8192 + (g.u64() % 24576) as i16, (g.u64() % 400) as i16);
            for (i, v) in h.iter_mut().enumerate() {
                *v = top - step * i as i16 / 2;
            }
            let mut sign = [if g.u64() & 1 == 0 { 32767 } else { -32767 }; 40];
            let mut o = 0;
            unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), rr.as_mut_ptr(), &mut o) };
        } else if odd && g.u64().is_multiple_of(4) {
            let shift = (g.u64() % 16) as u32;
            for i in 0..40 {
                for j in 0..=i {
                    rr[i][j] = g.i16() >> shift;
                    rr[j][i] = rr[i][j];
                }
            }
        } else {
            let mut h = [0i16; 40];
            let mut sign = [0i16; 40];
            if odd {
                impulse(g, &mut h);
                signs(g, &mut sign);
            } else {
                // As the encoder's: decaying, first tap large, ±32767 signs.
                let scale = 2 + (g.u64() % 4) as u32;
                for (i, v) in h.iter_mut().enumerate() {
                    *v = (g.u64() as i16) >> (scale + (i as u32 / 6).min(12));
                }
                h[0] = 4096 + (g.u64() % 8192) as i16;
                for s in &mut sign {
                    *s = if g.u64() & 1 == 0 { 32767 } else { -32767 };
                }
            }
            let mut o = 0;
            unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), rr.as_mut_ptr(), &mut o) };
        }
        let mut dn = [0i16; 40];
        if flat {
            for d in &mut dn {
                *d = g.i16() >> (6 + g.u64() % 4);
            }
        } else if !odd {
            // the timing class: `dn` as `cor_h_x` normalises it (2 sum of the tracks' max |dn|
            // within 32767: the fast path; full-scale values took the exact search every time)
            for d in &mut dn {
                *d = g.i16() >> (4 + g.u64() % 4);
            }
        } else if odd && g.u64().is_multiple_of(4) {
            for d in &mut dn {
                *d = [0, 1, 8192, 16384, -16384, -8192][(g.u64() % 6) as usize];
            }
        } else {
            signal(g, &mut dn);
        }
        let (t, n) = (T as u64, 40 / T as u64);
        let a = (g.u64() % t) as i16;
        let ipos = core::array::from_fn(|k| (a + k as i16) % t as i16);
        let pos_max = core::array::from_fn(|p| p as i16 + (t * (g.u64() % n)) as i16);
        (dn, rr, ipos, pos_max)
    }

    /// `search_10i40` (MR122) and `search_8i40_portable` (MR102, both widths) vs the C's
    /// `search_10and8i40`, each with both `Pair::best` forms (`L32`).
    fn search_matches_c<const T: usize, const P: usize>(
        seed: u64,
        rust: impl Fn(&[i16; 40], &[[i16; 40]; 40], &mut [i16; P], &[i16; T], &mut [usize; P]),
    ) {
        let mut g = Rng::new(seed);
        for _ in 0..N / 4 {
            let (mut dn, mut rr, mut ipos, mut pos_max) = search_in::<T, P>(&mut g, true);
            let (mut rpos, mut cv) = (ipos, [0usize; P]);
            rust(&dn, &rr, &mut rpos, &pos_max, &mut cv);
            let (mut ccv, mut o) = ([0i16; P], 0);
            unsafe {
                ref_search_10and8i40_w(
                    P as i16,
                    T as i16,
                    T as i16,
                    dn.as_mut_ptr(),
                    rr.as_mut_ptr(),
                    ipos.as_mut_ptr(),
                    pos_max.as_mut_ptr(),
                    ccv.as_mut_ptr(),
                    &mut o,
                )
            };
            assert_eq!(
                (cv.map(|v| v as i16), rpos),
                (ccv, ipos),
                "{dn:?} {pos_max:?}"
            );
        }
    }

    #[test]
    fn search_10i40_matches_c() {
        search_matches_c::<5, 10>(27, |dn, rr, ip, pm, cv| {
            search_10i40(dn, &tm(rr), ip, pm, cv)
        });
        search_matches_c::<5, 10>(27, |dn, rr, ip, pm, cv| {
            search_10i40_portable::<false>(dn, &tm(rr), ip, pm, cv)
        });
        search_matches_c::<5, 10>(27, |dn, rr, ip, pm, cv| {
            search_10i40_portable::<true>(dn, &tm(rr), ip, pm, cv)
        });
    }

    #[test]
    fn search_8i40_matches_c() {
        search_matches_c::<4, 8>(29, |dn, rr, ip, pm, cv| {
            search_8i40_portable::<10, false>(dn, &tm::<4, 10, 10>(rr), ip, pm, cv)
        });
        search_matches_c::<4, 8>(29, |dn, rr, ip, pm, cv| {
            search_8i40_portable::<16, true>(dn, &tm::<4, 10, 16>(rr), ip, pm, cv)
        });
    }

    /// The step-5 searches' `step5_block` (8 first pulses in lanes, later pulses' rows
    /// transposed, plain wrapping energy sums: `STEP5_SHIFTS`' proof) is the C's nested loops
    /// (`for i = start; i < L_CODE; i += STEP` over the square `rr`) with TS 26.073's saturating
    /// operators, for every lane: random symmetric `rr` of any magnitude, a ±full-scale class
    /// (the chains at their proven limits, 4 pulses), `dn` near the i16 limits, and a tie-heavy
    /// class (few distinct values).
    #[test]
    fn step5_block_is_the_c_loops() {
        use crate::basic_op::{add, l_mac, l_msu, l_mult, mult, round};
        use crate::encoder::codebook::{dn_tm, step5_block};
        let mut g = Rng::new(29);
        for case in 0..N {
            let ties = case % 2 == 0;
            // the 16-bit extremes (np = 4, not `ties`): the chains at their range limits
            let full = case % 6 == 5;
            let shift = (g.u64() % 16) as u32;
            let mut rr = [[0i16; 40]; 40];
            for i in 0..40 {
                for j in 0..=i {
                    rr[i][j] = if ties {
                        [0, 1, 4096, -4096][(g.u64() % 4) as usize]
                    } else if full {
                        [32767, 32766, -32768, -32767][(g.u64() % 4) as usize]
                    } else {
                        g.i16() >> shift
                    };
                    rr[j][i] = rr[i][j];
                }
            }
            let mut dn = [0i16; 40];
            for d in &mut dn {
                *d = if ties {
                    [0, 100, -100, 16384][(g.u64() % 4) as usize]
                } else if g.u64().is_multiple_of(3) {
                    [32767, -32768, 16384][(g.u64() % 3) as usize]
                } else {
                    g.i16()
                };
            }
            let ipos: [usize; 4] = core::array::from_fn(|_| (g.u64() % 5) as usize);
            // The C's loop for one pulse: weights 1/4, 1/2 (pulse 1) or 1/16, 1/8, the cross
            // terms from the latest pulse back.
            let c_loop = |ps0: i16, alp0: i32, t: usize, prev: &[usize], p: usize, o: &mut bool| {
                let (wd, wc) = if p == 1 { (8192, 16384) } else { (2048, 4096) };
                let (mut sq, mut alp, mut ps, mut ix) = (-1i16, 1i16, 0i16, t);
                for i in (t..40).step_by(5) {
                    let ps1 = add(ps0, dn[i], o);
                    let mut alp1 = l_mac(alp0, rr[i][i], wd, o);
                    for &q in prev.iter().rev() {
                        alp1 = l_mac(alp1, rr[q][i], wc, o);
                    }
                    let sq1 = mult(ps1, ps1, o);
                    let alp_16 = round(alp1, o);
                    if l_msu(l_mult(alp, sq1, o), sq, alp_16, o) > 0 {
                        (sq, ps, alp, ix) = (sq1, ps1, alp_16, i);
                    }
                }
                (sq, ps, alp, ix)
            };
            let dt = dn_tm(&dn);
            let rt: RrTm = tm(&rr);
            let np = 2 + case % 3;
            let mut o = false;
            let (got, gix) = match np {
                2 => {
                    let (l, ix) = step5_block(&[ipos[0], ipos[1]], &dt, &rt);
                    (l, ix.to_vec())
                }
                3 => {
                    let (l, ix) = step5_block(&[ipos[0], ipos[1], ipos[2]], &dt, &rt);
                    (l, ix.to_vec())
                }
                _ => {
                    let (l, ix) = step5_block(&ipos, &dt, &rt);
                    (l, ix.to_vec())
                }
            };
            for l in 0..8 {
                let i0 = ipos[0] + 5 * l;
                let mut pulses = vec![i0];
                let (mut ps, mut alp0) = (dn[i0], (rr[i0][i0] as i32) << 14);
                let mut last = (0, 0, 0, 0);
                for p in 1..np {
                    last = c_loop(ps, alp0, ipos[p], &pulses, p, &mut o);
                    pulses.push(last.3);
                    ps = last.1;
                    alp0 = (last.2 as i32) << if p == 1 { 14 } else { 16 };
                }
                let lane: Vec<usize> = (0..np).map(|p| ipos[p] + 5 * gix[p][l] as usize).collect();
                assert_eq!(
                    (got.sq[l], got.ps[l], got.alp[l], lane),
                    (last.0, last.1, last.2, pulses),
                    "case {case} lane {l} {ipos:?}"
                );
            }
        }
    }

    /// `simd::transpose8` (SSE2 intrinsics on x86_64) is the transpose, and so is the array form.
    #[test]
    fn transpose8_is_the_transpose() {
        use crate::encoder::codebook::transpose8_portable;
        let mut g = Rng::new(32);
        for _ in 0..1000 {
            let r: [[i16; 8]; 8] = core::array::from_fn(|_| core::array::from_fn(|_| g.i16()));
            let want: [[i16; 8]; 8] = core::array::from_fn(|n| core::array::from_fn(|l| r[l][n]));
            assert_eq!(crate::simd::transpose8(r), want);
            assert_eq!(transpose8_portable(r), want);
        }
    }

    /// `search_4i40` (MR74/MR795; packed half blocks, `step5_halves`) is the C's four nested
    /// loops with the `dn2 >= 0` skip and the outer compare chain, both instances (dispatched
    /// and portable). Classes: random `rr`/`dn`, tie-heavy (few distinct values), and tracks
    /// with 5..=8 cells of |dn| = 32767, which `set_sign` (its stale `pos`) leaves with 5..=8
    /// valid first pulses: blocks split into a full and a partial half, pairs then mix blocks.
    #[test]
    fn search_4i40_is_the_c_loops() {
        use crate::encoder::codebook::search_4i40_portable;
        let mut g = Rng::new(31);
        let (mut wide, mut ties_n) = (0, 0);
        for case in 0..N {
            let ties = case % 3 == 0;
            let shift = (g.u64() % 16) as u32;
            let mut rr = [[0i16; 40]; 40];
            for i in 0..40 {
                for j in 0..=i {
                    rr[i][j] = if ties {
                        [0, 1, 4096, -4096][(g.u64() % 4) as usize]
                    } else {
                        g.i16() >> shift
                    };
                    rr[j][i] = rr[i][j];
                }
            }
            let mut dn = [0i16; 40];
            for d in &mut dn {
                *d = if ties {
                    [0, 100, -100, 16384][(g.u64() % 4) as usize]
                } else {
                    g.i16() >> (g.u64() % 4)
                };
            }
            if case % 3 == 1 {
                // A track with k > 4 cells of |dn| = 32767 keeps k first pulses (only 8 - k
                // removals find a minimum below 32767).
                for t in 0..5 {
                    let k = 5 + (g.u64() % 4) as usize;
                    for n in 0..8 {
                        if n < k {
                            dn[t + 5 * n] = [32767, -32767, -32768][(g.u64() % 3) as usize];
                        } else if dn[t + 5 * n].unsigned_abs() >= 32767 {
                            dn[t + 5 * n] = 0;
                        }
                    }
                }
            }
            let (mut sign, mut dn2) = ([0i16; 40], [0i16; 40]);
            set_sign(&mut dn, &mut sign, &mut dn2, 4);
            let valid: Vec<usize> = (0..5)
                .map(|t| (0..8).filter(|n| dn2[t + 5 * n] >= 0).count())
                .collect();
            wide += valid.iter().any(|&v| v > 4) as usize;
            ties_n += ties as usize;
            // The C's loops with TS 26.073's saturating operators.
            use crate::basic_op::{add, l_mac, l_msu, l_mult, mult, round};
            let o = &mut false;
            let mut want = [0usize, 1, 2, 3];
            let (mut psk, mut alpk) = (-1i16, 1i16);
            for track in 3..5 {
                let mut ipos = [0usize, 1, 2, track];
                for _ in 0..4 {
                    for i0 in (ipos[0]..40).step_by(5) {
                        if dn2[i0] < 0 {
                            continue;
                        }
                        let mut pulses = vec![i0];
                        let (mut ps0, mut alp0) = (dn[i0], l_mult(rr[i0][i0], 8192, o));
                        let (mut sq, mut alp) = (-1i16, 1i16);
                        for p in 1..4 {
                            let (wd, wc) = if p == 1 { (8192, 16384) } else { (2048, 4096) };
                            let mut ps = 0i16;
                            let mut ix = ipos[p];
                            (sq, alp) = (-1, 1);
                            for i in (ipos[p]..40).step_by(5) {
                                let ps1 = add(ps0, dn[i], o);
                                let mut alp1 = l_mac(alp0, rr[i][i], wd, o);
                                for &q in pulses.iter().rev() {
                                    alp1 = l_mac(alp1, rr[q][i], wc, o);
                                }
                                let sq1 = mult(ps1, ps1, o);
                                let alp_16 = round(alp1, o);
                                if l_msu(l_mult(alp, sq1, o), sq, alp_16, o) > 0 {
                                    (sq, ps, alp, ix) = (sq1, ps1, alp_16, i);
                                }
                            }
                            pulses.push(ix);
                            ps0 = ps;
                            alp0 = if p == 1 {
                                l_mult(alp, 8192, o)
                            } else {
                                (alp as i32) << 16
                            };
                        }
                        if l_msu(l_mult(alpk, sq, o), psk, alp, o) > 0 {
                            (psk, alpk) = (sq, alp);
                            want = [pulses[0], pulses[1], pulses[2], pulses[3]];
                        }
                    }
                    ipos.rotate_right(1);
                }
            }
            let rt: RrTm = tm(&rr);
            let got = crate::simd::search_4i40(&dn, &dn2, &rt);
            assert_eq!(got, want, "case {case} valid {valid:?}");
            let got = search_4i40_portable(&dn, &dn2, &rt);
            assert_eq!(got, want, "case {case} portable, valid {valid:?}");
        }
        assert!(
            wide >= N / 4 && ties_n >= N / 4,
            "classes: wide {wide}, ties {ties_n}"
        );
    }

    /// `lanes_next` is the C's loop (TS 26.073 operators) for every lane, on its callers'
    /// domain: per pulse level the start energy as they make it (pulse 1: rr << 14; pulse 2:
    /// alp << 14, alp any i16; pulse 3: alp << 16, |alp| <= 18432) with 1, 2 or 3 cross terms
    /// of any 16-bit value, `d` and `diag` any 16-bit rows (two half blocks' rows, as
    /// `step5_halves` passes them). The `const` assertion after `lanes_next` proves no chain
    /// leaves the i32 range there; the C side here sums in i64 and asserts that too.
    #[test]
    fn lanes_next_is_the_c_loop() {
        use crate::basic_op::{add, mult, round};
        use crate::encoder::codebook::lanes_next;
        let mut g = Rng::new(30);
        for case in 0..N {
            let level = case % 3;
            let full = case % 5 == 0; // the 16-bit extremes only
            let v16 = |g: &mut Rng| {
                if full {
                    [i16::MIN, i16::MAX][(g.u64() % 2) as usize]
                } else {
                    g.i16() >> (g.u64() % 16)
                }
            };
            let ps0: [i16; 8] = core::array::from_fn(|_| g.i16());
            let alp0: [i32; 8] = core::array::from_fn(|_| match level {
                0 | 1 => (v16(&mut g) as i32) << 14,
                _ => ((v16(&mut g) as i32).clamp(-18432, 18431)) << 16,
            });
            let (sd, sc) = if level == 0 { (14, 15) } else { (12, 13) };
            let dh: [[i16; 8]; 2] = core::array::from_fn(|_| core::array::from_fn(|_| g.i16()));
            let gh: [[i16; 8]; 2] = core::array::from_fn(|_| core::array::from_fn(|_| v16(&mut g)));
            let (d, diag) = ([&dh[0], &dh[1]], [&gh[0], &gh[1]]);
            let x: [[i32; 8]; 8] = core::array::from_fn(|_| {
                core::array::from_fn(|_| (0..=level).map(|_| v16(&mut g) as i32).sum())
            });
            let got = lanes_next(&ps0, &alp0, &d, &diag, &x, sd, sc);
            let mut o = false;
            for l in 0..8 {
                let (mut sq, mut alp, mut ps, mut ix) = (-1i16, 1i16, 0i16, 0i16);
                for n in 0..8 {
                    let ps1 = add(ps0[l], dh[l / 4][n], &mut o);
                    let alp1 =
                        alp0[l] as i64 + ((gh[l / 4][n] as i64) << sd) + ((x[n][l] as i64) << sc);
                    assert!(i32::try_from(alp1).is_ok(), "level {level}: {alp1}");
                    let sq1 = mult(ps1, ps1, &mut o);
                    let alp_16 = round(alp1 as i32, &mut o);
                    if alp as i32 * sq1 as i32 > sq as i32 * alp_16 as i32 {
                        (sq, ps, alp, ix) = (sq1, ps1, alp_16, n as i16);
                    }
                }
                assert_eq!(
                    (got.sq[l], got.ps[l], got.alp[l], got.ix[l]),
                    (sq, ps, alp, ix),
                    "case {case} lane {l}"
                );
            }
        }
    }

    /// `Pair::best` (per-lane argmax as the rows are made) is the compare chain over the pair's
    /// candidates, both forms (`L32`): random rows (energies of any sign and size: the
    /// fallback), a tie-heavy class (few distinct values), both energy forms, random padding
    /// lanes; and an edge class for the range check: valid wide rows with one cell at an
    /// energy of 0, -1, 1, 32767 or 32768, or ps = -32768 (sq = -32768 in the C), also where
    /// only they could differ (that cell first, all others sq = 0 or it also sq = 0).
    fn best_is_the_chain<const R: usize, const W: usize>(seed: u64) {
        let mut g = Rng::new(seed);
        let (mut fast, mut edges, mut edge_fast) = (0, 0, 0);
        for case in 0..N {
            if case % 4 == 1 {
                let ok = best_edge_case::<R, W>(&mut g, edges % 4);
                (edges, edge_fast) = (edges + 1, edge_fast + ok as usize);
                continue;
            }
            let odd = g.u64().is_multiple_of(4);
            let wide = g.u64() & 1 == 0;
            let pad_ps = g.i16();
            let (mut d1, mut tps, mut base, mut tal) = ([0i16; W], [0i16; W], [0i32; W], [0i32; W]);
            let mut rr = [[0i16; W]; R];
            for i in 0..W {
                if i >= R {
                    // Padding as `search_pairs` makes it: rr 0, tal 0, tps = ps + a zero dn
                    // lane (d1 is per row: unused here).
                    (tps[i], base[i], tal[i]) = (pad_ps, base[0], 0);
                } else if odd {
                    (d1[i], tps[i]) = (g.i16(), g.i16());
                    // The wide pair's adds don't wrap in the C either: `alp1` is a few i16
                    // terms (|alp1| < 2^18), `tal` an i16. The other pairs wrap: any i32.
                    (base[i], tal[i]) = if wide {
                        (g.i32() >> (13 + g.u64() % 19), g.i16() as i32)
                    } else {
                        (g.i32() >> (g.u64() % 32), g.i32() >> (g.u64() % 32))
                    };
                    for row in &mut rr {
                        row[i] = g.i16();
                    }
                } else {
                    d1[i] = few(&mut g, &[0, 3000, 6000]) as i16;
                    tps[i] = few(&mut g, &[0, 3000, 9000]) as i16;
                    // Energies in 1..=32767 either way (wide: (base + rr) / 8 + tal / 2).
                    (base[i], tal[i]) = if wide {
                        (
                            few(&mut g, &[8, 4000, 40000]),
                            few(&mut g, &[2, 1000, 8000]),
                        )
                    } else {
                        (few(&mut g, &[1 << 16, 1 << 28]), few(&mut g, &[0, 1 << 20]))
                    };
                    for row in &mut rr {
                        row[i] = few(&mut g, &[0, 100, 2000]) as i16;
                    }
                }
            }
            let p = Pair::<R, W> {
                d1: &d1,
                tps,
                base,
                tal,
                rr: &rr,
                sh: 10 + (g.u64() % 3) as u32,
                wide,
            };
            let c = p.cands();
            let ok = (0..R)
                .all(|j| (0..R).all(|i| c.sq[j][i] >= 0 && (1..=32767).contains(&c.al[j][i])));
            fast += ok as usize;
            let want = chain(&c);
            assert_eq!(p.best::<false>(), want, "{:?} {:?}", c.sq, c.al);
            assert_eq!(p.best::<true>(), want, "{:?} {:?}", c.sq, c.al);
        }
        assert!(fast > 3 * N / 8, "{fast}");
        // Both sides of the range check: valid edge cells (kind 2, energy 1 or 32767) take the
        // lanes, the others fall back.
        assert!(
            edge_fast > edges / 20 && edge_fast < edges / 4,
            "{edge_fast} of {edges}"
        );
    }

    /// One of `v`, uniformly.
    fn few(g: &mut Rng, v: &[i32]) -> i32 {
        v[(g.u64() % v.len() as u64) as usize]
    }

    /// One case of `best_is_the_chain`'s edge class; whether every candidate is in the fast
    /// range. Wide energies `(((base + rr) >> 2) + tal) >> 1` with base 800, so tal 0 or 65534
    /// and rr pick the energy exactly. Kinds: 0 ps = -32768 anywhere, 1 on all of row 0 (one tps)
    /// with every other sq 0 (lanes reading it as sq 0 take it first), 2 an edge energy
    /// anywhere, 3 energy 0 and sq 0 at (0, 0) (nothing beats 0 / 0 in a lane; the C skips it).
    fn best_edge_case<const R: usize, const W: usize>(g: &mut Rng, kind: usize) -> bool {
        let (mut d1, mut tps, mut tal) = ([0i16; W], [0i16; W], [0i32; W]);
        let mut rr = [[0i16; W]; R];
        let ps: &[i32] = if kind == 1 {
            &[0, 100]
        } else {
            &[0, 3000, 9000]
        };
        let pad_ps = if kind == 1 { 0 } else { few(g, ps) as i16 };
        for i in 0..W {
            if i >= R {
                tps[i] = pad_ps;
                continue;
            }
            d1[i] = few(g, ps) as i16;
            tps[i] = if kind == 1 && i > 0 {
                tps[0]
            } else {
                few(g, ps) as i16
            };
            tal[i] = few(g, &[0, 65534]);
            // Energies 1, 100, 1000 (tal 0) or 32767, 32757 (tal 65534): rr - 800 = 8 al - tal.
            let v: &[i32] = if tal[i] == 0 {
                &[8, 800, 8000]
            } else {
                &[0, -80]
            };
            for row in &mut rr {
                row[i] = (few(g, v) - 800) as i16;
            }
        }
        let (j0, i0) = match kind {
            1 | 3 => (0, 0),
            _ => ((g.u64() % R as u64) as usize, (g.u64() % R as u64) as usize),
        };
        match kind {
            // ps = d1[j0] + tps[i0] = -32768.
            0 | 1 => d1[j0] = (-32768 - tps[i0] as i32) as i16,
            // al 0, -1 or 1 (tal 0); 32768 or 32767 (tal 65534).
            2 => {
                let v: &[i32] = if tal[i0] == 0 { &[0, -8, 8] } else { &[8, 0] };
                rr[j0][i0] = (few(g, v) - 800) as i16;
            }
            _ => {
                (tal[i0], d1[j0]) = (0, -tps[i0]);
                rr[j0][i0] = -800;
            }
        }
        let p = Pair::<R, W> {
            d1: &d1,
            tps,
            base: [800; W],
            tal,
            rr: &rr,
            sh: 0,
            wide: true,
        };
        let c = p.cands();
        let ok =
            (0..R).all(|j| (0..R).all(|i| c.sq[j][i] >= 0 && (1..=32767).contains(&c.al[j][i])));
        let want = chain(&c);
        assert_eq!(p.best::<false>(), want, "{:?} {:?}", c.sq, c.al);
        assert_eq!(p.best::<true>(), want, "{:?} {:?}", c.sq, c.al);
        ok
    }

    /// MR122's 8 × 8 candidates, MR102's 10 × 10 in 16 lanes (AVX2) and 10 (portable).
    #[test]
    fn best_is_the_chain_8_10() {
        best_is_the_chain::<8, 8>(28);
        best_is_the_chain::<10, 16>(30);
        best_is_the_chain::<10, 10>(32);
    }

    /// `cor_h_x` (5 tracks) and `cor_h_x2` (MR102's 4) vs the C, both `sf`; full-scale classes
    /// (`signal` edge values, `impulse`'s full-scale `h`) make the sums wrap and `l_shl`
    /// saturate.
    #[test]
    fn cor_h_x_matches_c() {
        let mut g = Rng::new(31);
        for case in 0..N / 2 {
            let (mut x, mut h) = ([0i16; 40], [0i16; 40]);
            signal(&mut g, &mut x);
            impulse(&mut g, &mut h);
            if case % 7 == 0 {
                x.iter_mut()
                    .for_each(|v| *v = if g.u64() & 1 == 0 { 32767 } else { -32768 });
            }
            let sf = 1 + (case & 1) as i16;
            let (mut dn, mut ov) = ([0i16; 40], false);
            let (mut cdn, mut co) = ([0i16; 40], 0 as CFlag);
            if case % 3 == 0 {
                crate::encoder::codebook::cor_h_x2(&h, &x, &mut dn, sf, 4, 4, &mut ov);
                unsafe {
                    c_cor_h_x2(
                        h.as_mut_ptr(),
                        x.as_mut_ptr(),
                        cdn.as_mut_ptr(),
                        sf,
                        4,
                        4,
                        &mut co,
                    )
                };
            } else {
                cor_h_x(&h, &x, &mut dn, sf, &mut ov);
                unsafe {
                    c_cor_h_x(
                        h.as_mut_ptr(),
                        x.as_mut_ptr(),
                        cdn.as_mut_ptr(),
                        sf,
                        &mut co,
                    )
                };
            }
            assert_eq!((dn, ov), (cdn, co != 0), "case {case}: {x:?} {h:?}");
        }
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

    /// Pitch search signal: PIT_MAX history + up to one frame, scaled as `pitch_ol` does.
    fn pitch_signal(g: &mut Rng) -> [i16; PIT_MAX + L_FRAME] {
        let mut sig = [0i16; PIT_MAX + L_FRAME];
        signal(g, &mut sig);
        sig
    }

    /// C `comp_corr` output (`out[m]` = lag `lag_max - m`) in Rust's `corr[lag]` layout.
    fn c_comp_corr_at(sig: &mut [i16], l: usize, max: usize, min: usize) -> [i32; PIT_MAX + 1] {
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

    /// `Pitch_fr` inputs: excitation history (`exc[EXC]` is the subframe start).
    const EXC: usize = PIT_MAX + L_INTERPOL;
    struct PitchIn {
        mode: Mode,
        t_op: [i16; 2],
        t0_prev: i16,
        i_subfr: i16,
        exc: [i16; EXC + 40],
        xn: [i16; 40],
        h: [i16; 40],
    }

    fn pitch_in(g: &mut Rng) -> PitchIn {
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
    fn pitch_fr_rust(p: &PitchIn) -> (i16, i16, bool, i16, i16, bool) {
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

    fn pitch_fr_c(p: &PitchIn) -> (i16, i16, bool, i16, i16, bool) {
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
    fn pitch_fr_matches_c() {
        let mut g = Rng::new(25);
        for _ in 0..N / 10 {
            let mut p = pitch_in(&mut g);
            if g.u64().is_multiple_of(4) {
                // Full-scale target: <xn, excf> saturates while <excf, excf> stays small
                // (reaches norm_corr's Cauchy-Schwarz bound).
                p.xn = [[32767, -32768][g.u64() as usize & 1]; 40];
            }
            assert_eq!(pitch_fr_rust(&p), pitch_fr_c(&p), "{} {:?}", p.mode, p.t_op);
        }
    }

    /// `Qua_gain` inputs: realistic coefficient exponents, sometimes extreme (saturating) ones.
    struct GainIn {
        mode: Mode,
        exp_gcode0: i16,
        frac_gcode0: i16,
        frac: [i16; 5],
        exp: [i16; 5],
        gp_limit: i16,
    }

    fn gain_in(g: &mut Rng) -> GainIn {
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

    fn qua_gain_rust(x: &GainIn) -> [i16; 5] {
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

    fn qua_gain_c(x: &GainIn, t: &mut CommonAmrTbls) -> [i16; 5] {
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
                mode, &v[0], &v[1], &v[2], &v[3], &gc, &mut rf, &mut re, &mut rcf, &mut rce,
                &mut ov,
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
            assert_eq!((rf, re, rcf, rce), (cf, ce, ccf, cce), "{mode} {v:?}");
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

    /// `sign_cor` (`set_sign12k2`'s lanes) vs the operator chain for any gains >= 0: the clamps
    /// at both ends, which `set_sign12k2`'s normalised gains don't reach.
    #[test]
    fn sign_cor_is_the_operator_chain() {
        use crate::basic_op::{l_mac, l_mult, l_shl, round};
        const L_CODE: usize = 40;
        let mut g = Rng::new(41);
        let mut clamped = 0;
        for _ in 0..N / 10 {
            let (kc, kd) = (g.i16().max(0), g.i16().max(0));
            let cn: [i16; L_CODE] = core::array::from_fn(|_| g.i16());
            let dn: [i16; L_CODE] = core::array::from_fn(|_| g.i16());
            let o = &mut false;
            let want: [i16; L_CODE] = core::array::from_fn(|i| {
                round(l_shl(l_mac(l_mult(kc, cn[i], o), kd, dn[i], o), 10, o), o)
            });
            assert_eq!(super::super::codebook::sign_cor(kc, &cn, kd, &dn), want);
            clamped += want
                .iter()
                .filter(|&&v| v == i16::MAX || v == i16::MIN)
                .count();
        }
        assert!(clamped > 1000, "clamped: {clamped}");
    }

    /// `set_sign12k2` (its per-sample pass as selects, the `cn` energy as an exact sum) vs the
    /// C, both track layouts: speech-like `dn`/`cn`; full-scale ones (-32768, the energy sum
    /// and `L_mac` saturate, the flag); constant ones (equal `en`: ties in the track maxima).
    #[test]
    fn set_sign12k2_matches_c() {
        let mut g = Rng::new(30);
        let (mut sat, mut ties, mut lmac) = (0, 0, 0);
        for n in 0..N {
            let (mut dn, mut cn) = ([0i16; 40], [0i16; 40]);
            if n % 8 == 3 {
                // L_mac saturating (only reachable this way: |k_cn * cn << 1| < 2^20, so it
                // needs k_dn * dn ~ 2^30, and the energy is even): sum dn^2 = 2^31 wraps the
                // energy to 256, k_dn = Inv_sqrt(256) >> 11 = 32767, at dn = 32767 that is
                // 2 * k_dn * dn = 2^31 - 131070; cn = 3..10 alone keeps k_cn near 32767, so
                // k_cn * cn << 1 > 131070 tips the sum over.
                let mut rest = 1i64 << 31;
                for v in dn.iter_mut() {
                    let r = rest.isqrt().min(32767);
                    *v = r as i16;
                    rest -= r * r;
                }
                assert_eq!((rest, dn[0], dn[1]), (0, 32767, 32767));
                cn[(g.u64() % 2) as usize] = 3 + (g.u64() % 8) as i16;
            } else {
                match n % 4 {
                    0 => {
                        const V: [i16; 4] = [32767, -32768, 30000, -30000];
                        for v in dn.iter_mut().chain(&mut cn) {
                            *v = if g.u64() & 1 == 0 {
                                V[(g.u64() % 4) as usize]
                            } else {
                                g.i16()
                            };
                        }
                    }
                    1 => {
                        (dn, cn) = (
                            [(g.u64() % 64) as i16 - 32; 40],
                            [(g.u64() % 64) as i16; 40],
                        );
                        ties += 1;
                    }
                    _ => {
                        signal(&mut g, &mut dn);
                        signal(&mut g, &mut cn);
                    }
                }
            }
            let (nb, step) = if n & 4 == 0 { (5, 5) } else { (4, 4) };
            let (mut rd, mut rs, mut rp, mut ri, mut ov) =
                (dn, [0i16; 40], [0i16; 5], [0i16; 10], false);
            set_sign12k2(&mut rd, &cn, &mut rs, &mut rp, nb, &mut ri, step, &mut ov);
            let (mut cd, mut cc, mut cs, mut cp, mut ci, mut o) =
                (dn, cn, [0i16; 40], [0i16; 5], [0i16; 10], 0);
            unsafe {
                c_set_sign12k2(
                    cd.as_mut_ptr(),
                    cc.as_mut_ptr(),
                    cs.as_mut_ptr(),
                    cp.as_mut_ptr(),
                    nb as i16,
                    ci.as_mut_ptr(),
                    step as i16,
                    &mut o,
                )
            };
            sat += (o != 0) as usize;
            // In that class only L_mac can saturate (energies, round and negate don't).
            lmac += (n % 8 == 3 && o != 0) as usize;
            assert_eq!(
                (rd, rs, rp, ri, ov),
                (cd, cs, cp, ci, o != 0),
                "{dn:?} {cn:?}"
            );
        }
        assert!(
            sat > N / 20 && ties > N / 8 && lmac > N / 10,
            "classes: {sat} {ties} {lmac}"
        );
    }

    /// `filtered_code` in lanes vs the C's chain step by step (`L_mult`, `L_mac`, `round`: the
    /// ported basic ops, each checked against the C): full-scale `h` with -32768 taps, signs
    /// +-32767/-32768 or small sets, pulses stacked on few positions (saturation, and `L_mac`'s
    /// "accumulator equals the product" quirk); the overflow flag too.
    #[test]
    fn filtered_code_matches_chain() {
        use crate::basic_op::{l_mac, l_mult, round};
        let mut g = Rng::new(28);
        let (mut n_ov, mut n_quirk) = (0, 0);
        const V: [i16; 6] = [32767, -32768, 16384, -16384, 21846, -21846];
        for n in 0..N {
            let class = n % 3;
            let mut h = [0i16; 40];
            for v in &mut h {
                *v = match class {
                    0 => g.u64() as i16 >> (g.u64() % 8),
                    _ => V[g.u64() as usize % 6],
                };
            }
            let np = 1 + (g.u64() % 10) as usize;
            let (mut cv, mut sg) = ([0usize; 10], [0i16; 10]);
            for k in 0..np {
                cv[k] = (g.u64() % if class == 2 { 3 } else { 40 }) as usize;
                sg[k] = match class {
                    0 => [32767, -32768][g.u64() as usize % 2],
                    _ => V[g.u64() as usize % 6],
                };
            }
            let (mut y, mut ov) = ([0i16; 40], false);
            filtered_code(&h, &cv[..np], &sg[..np], &mut y, &mut ov);
            let (mut yr, mut ovr) = ([0i16; 40], false);
            for i in 0..40 {
                let at = |c: usize| if i >= c { h[i - c] } else { 0 };
                let mut s = l_mult(at(cv[0]), sg[0], &mut ovr);
                for k in 1..np {
                    n_quirk += (s == at(cv[k]) as i32 * sg[k] as i32 && s != 0) as usize;
                    s = l_mac(s, at(cv[k]), sg[k], &mut ovr);
                }
                yr[i] = round(s, &mut ovr);
            }
            n_ov += ovr as usize;
            assert_eq!((y, ov), (yr, ovr), "{h:?} {cv:?} {sg:?}");
        }
        // n_quirk: L_var3 == product, where opencore skipped the saturation
        assert!(
            n_ov > N / 10 && n_quirk > 100,
            "classes too rare: {n_ov} {n_quirk}"
        );
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

    /// A C codebook, its columns, its size, `use_half`.
    type VqCase<const D: usize> = (&'static [i16], [&'static [i16]; D], i16, bool);

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
                    let t = ((wf[k] as i32 * r[k].wrapping_sub(DICO1_LSF_3[k]) as i32) >> 15) as i16
                        as i64;
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
                let ci = unsafe {
                    c_vq_subvec4(r.as_mut_ptr(), d.as_ptr(), wf.as_mut_ptr(), size, &mut o)
                };
                assert_eq!((i, rr), (ci, r), "{size} {wf:?}");
            }
        }
        assert!(
            ties > N / 1000 && wraps > N / 1000,
            "classes: {ties} {wraps}"
        );
    }

    /// ns per call of `f` over 64 rotating inputs (best of 20 runs of 20k calls). Only the
    /// kernels whose name contains `$KERN`, if set.
    fn time(name: &str, mut f: impl FnMut(usize)) {
        if std::env::var("KERN").is_ok_and(|k| !name.contains(&k)) {
            return;
        }
        let best = (0..20)
            .map(|_| {
                let t = std::time::Instant::now();
                for i in 0..20_000 {
                    f(i & 63);
                }
                t.elapsed().as_secs_f64() / 20_000.0 * 1e9
            })
            .fold(f64::INFINITY, f64::min);
        eprintln!("{name:<24} {best:8.1} ns");
    }

    /// Kernel timings, Rust vs C, ns per call (`BENCH=kern just bench`). Whole-codec timing is
    /// too coarse for one kernel, and callgrind's instruction and mispredict counts mislead.
    #[test]
    #[ignore]
    fn kernel_timing() {
        use std::hint::black_box as bb;
        let mut g = Rng::new(99);
        let mut xs = [[0i16; 40]; 64];
        let mut hs = [[0i16; 40]; 64];
        for (x, h) in xs.iter_mut().zip(&mut hs) {
            // the encoder's kind of input (an excitation, a decaying impulse response): the
            // fast path; full-scale `signal`/`impulse` took the exact loop every call
            for v in x.iter_mut() {
                *v = g.i16() >> (3 + g.u64() % 4);
            }
            let scale = 2 + (g.u64() % 4) as u32;
            for (i, v) in h.iter_mut().enumerate() {
                *v = g.i16() >> (scale + (i as u32 / 6).min(12));
            }
            h[0] = 4096 + (g.u64() % 8192) as i16;
        }
        let mut y = [0i16; 40];
        time("convolve", |i| convolve(bb(&xs[i]), bb(&hs[i]), bb(&mut y)));
        time("convolve C", |i| unsafe {
            oc_convolve(
                bb(xs[i].as_ptr() as *mut i16),
                bb(hs[i].as_ptr() as *mut i16),
                bb(y.as_mut_ptr()),
                40,
            )
        });
        let mut syn = Vec::new();
        while syn.len() < 64 {
            // the encoder's kind of input: LPC residual-like x and past outputs, no saturation
            // (the exact fallback is rare in the codec)
            let a = crate::post::test_gen::lpc(&mut g);
            let x: [i16; 40] = std::array::from_fn(|_| (g.u64() as i16) >> (6 + g.u64() % 4));
            let mut mem = [0i16; M];
            let mut y = [0i16; 40];
            crate::post::filt::syn_filt(&a, &x, &mut y, 40, &mut mem, true);
            // and within `syn_filt`'s bound (a0 max|x| + sum|a_j| max|y| < 2^30): 1.6% of the
            // old inputs took the exact loop, the decoder's calls 39 in 714k
            let am = |v: &[i16]| v.iter().map(|&c| (c as i64).abs()).max().unwrap_or(0);
            let sa: i64 = a[1..].iter().map(|&c| (c as i64).abs()).sum();
            if (a[0] as i64).abs() * am(&x) + sa * am(&y) < 1 << 30 {
                syn.push((a, x, mem));
            }
        }
        let mut ys = [0i16; 40];
        time("syn_filt 40", |i| {
            let (a, x, mut mem) = syn[i];
            crate::post::filt::syn_filt(bb(&a), bb(&x), bb(&mut ys), 40, bb(&mut mem), true);
        });
        time("syn_filt 40 C", |i| unsafe {
            let (mut a, mut x, mut mem) = syn[i];
            c_syn_filt(
                bb(a.as_mut_ptr()),
                bb(x.as_mut_ptr()),
                bb(ys.as_mut_ptr()),
                40,
                bb(mem.as_mut_ptr()),
                1,
            )
        });
        let mut sigs = [[0i16; PIT_MAX + L_FRAME]; 8];
        for sig in &mut sigs {
            *sig = pitch_signal(&mut g);
        }
        {
            unsafe extern "C" {
                #[link_name = "Pred_lt_3or6"]
                fn c_pred_lt(exc: *mut i16, t0: i16, frac: i16, l: i16, flag3: i16, o: *mut i32);
            }
            // the encoder's call: l_subfr 41, lags over the whole range, both resolutions,
            // frac -2..=2 (all but the 1/6 resolution's frac 3)
            const OFF: usize = PIT_MAX + 11;
            let args: [(i16, i16, i16); 64] = std::array::from_fn(|i| {
                let t0 = 20 + (g.u64() % 124) as i16;
                (t0, (g.u64() % 5) as i16 - 2, (i & 1) as i16)
            });
            let mut buf = [0i16; OFF + 41];
            // an excitation's range (full-scale `sigs` failed the window guard: the exact loop)
            let excs = sigs.map(|s| s.map(|v| v >> 4));
            time("pred_lt_3or6", |i| {
                let (t0, frac, flag3) = args[i];
                buf[..PIT_MAX + 11].copy_from_slice(&excs[i & 7][..PIT_MAX + 11]);
                crate::excitation::pred_lt_3or6(
                    bb(&mut buf),
                    OFF,
                    t0,
                    frac,
                    41,
                    flag3,
                    &mut Flag::default(),
                )
            });
            time("pred_lt_3or6 C", |i| unsafe {
                let (t0, frac, flag3) = args[i];
                buf[..PIT_MAX + 11].copy_from_slice(&excs[i & 7][..PIT_MAX + 11]);
                let mut o = 0;
                c_pred_lt(bb(buf.as_mut_ptr().add(OFF)), t0, frac, 41, flag3, &mut o)
            });
        }
        let mut corr = [0i32; PIT_MAX + 1];
        time("comp_corr 160", |i| {
            comp_corr(bb(&sigs[i & 7]), L_FRAME, PIT_MAX, 20, bb(&mut corr))
        });
        time("comp_corr 160 C", |i| {
            bb(c_comp_corr_at(bb(&mut sigs[i & 7]), L_FRAME, PIT_MAX, 20));
        });
        let pins: Vec<PitchIn> = (0..64).map(|_| pitch_in(&mut g)).collect();
        time("pitch_fr", |i| {
            bb(pitch_fr_rust(bb(&pins[i])));
        });
        time("pitch_fr C", |i| {
            bb(pitch_fr_c(bb(&pins[i])));
        });
        let mut lsf = [[0i16; 10]; 64];
        let mut wf = [[0i16; 10]; 64];
        for (l, w) in lsf.iter_mut().zip(&mut wf) {
            for (a, b) in l.iter_mut().zip(w) {
                (*a, *b) = (
                    (g.u64() % 2000) as i16 - 1000,
                    3000 + (g.u64() % 20000) as i16,
                );
            }
        }
        time("vq_subvec3 512", |i| {
            let mut r = lsf[i];
            bb(vq_subvec3(&mut r[3..], bb(&cols(&DICO2_SOA)), &wf[i][3..]));
        });
        time("vq_subvec4 512", |i| {
            let mut r = lsf[i];
            bb(vq_subvec4(&mut r[6..], bb(&cols(&DICO3_SOA)), &wf[i][6..]));
        });
        let (mut dn, mut ov) = ([0i16; 40], false);
        time("cor_h_x", |i| {
            cor_h_x(bb(&hs[i]), bb(&xs[i]), bb(&mut dn), 1, &mut ov)
        });
        let (mut sg, mut dn2) = ([0i16; 40], [0i16; 40]);
        time("set_sign", |i| {
            let mut d = xs[i];
            set_sign(bb(&mut d), &mut sg, &mut dn2, 4);
            bb(&dn2);
        });
        // LP analysis on speech-like windows: autocorr, and az_lsp on the LPC it leads to.
        let (wins, mut azs) = ([[0i16; 240]; 64], [[0i16; 11]; 64]);
        let mut wins = wins;
        let mut lev = crate::encoder::lpc::LevinsonState::default();
        for (x, a) in wins.iter_mut().zip(&mut azs) {
            signal(&mut g, x);
            let (mut rh, mut rl, mut rc, mut ov) = ([0i16; 11], [0i16; 11], [0i16; 4], false);
            crate::encoder::lpc::autocorr(x, 10, &mut rh, &mut rl, &crate::tables::WINDOW_200_40);
            crate::encoder::lpc::lag_window(10, &mut rh, &mut rl, &mut ov);
            crate::encoder::lpc::levinson(&mut lev, &rh, &rl, a, &mut rc, &mut ov);
        }
        let (mut rh, mut rl) = ([0i16; 11], [0i16; 11]);
        time("autocorr", |i| {
            bb(crate::encoder::lpc::autocorr(
                bb(&wins[i]),
                10,
                &mut rh,
                &mut rl,
                &crate::tables::WINDOW_200_40,
            ));
        });
        let mut lsp = [0i16; 10];
        time("az_lsp", |i| {
            crate::encoder::lpc::az_lsp(bb(&azs[i]), bb(&mut lsp), &crate::tables::LSP_INIT_DATA)
        });
        let gins: Vec<(GainIn, GainIn)> = (0..64)
            .map(|_| (gain_in(&mut g), gain_in(&mut g)))
            .collect();
        let mut code = [0i16; 40];
        (code[7], code[30]) = (8191, -8191);
        time("mr475_gain_quant", |i| {
            let (a, b) = &gins[i];
            let mut st = crate::gains::GcPredState::default();
            let (mut r, mut ov) = ([0i16; 4], false);
            let [p0, c0, p1, c1] = &mut r;
            bb(crate::encoder::gain_q::mr475_gain_quant(
                &mut st,
                a.exp_gcode0,
                a.frac_gcode0,
                &a.exp,
                &a.frac,
                10,
                20000,
                bb(&code),
                b.exp_gcode0,
                b.frac_gcode0,
                &b.exp,
                &b.frac,
                10,
                18000,
                a.gp_limit,
                p0,
                c0,
                p1,
                c1,
                &mut ov,
            ));
        });
        let cv = [3usize, 9, 21, 37];
        let sgn = [32767i16, -32768, 32767, -32768];
        time("filtered_code 4", |i| {
            filtered_code(bb(&hs[i]), bb(&cv), &sgn, bb(&mut y), &mut ov)
        });
        let mut t = tables();
        for mode in [Mode::MR59, Mode::MR67, Mode::MR74, Mode::MR102, Mode::MR122] {
            let cins: Vec<CbIn> = (0..64).map(|_| cb_in(&mut g, mode)).collect();
            time(&format!("cbsearch {mode}"), |i| {
                bb(cb_rust(bb(&cins[i])));
            });
            time(&format!("cbsearch {mode} C"), |i| {
                bb(cb_oc(bb(&cins[i]), &mut t));
            });
        }
        let gins: Vec<GainIn> = (0..64)
            .map(|_| GainIn {
                mode: Mode::MR102,
                ..gain_in(&mut g)
            })
            .collect();
        time("qua_gain 102", |i| {
            bb(qua_gain_rust(bb(&gins[i])));
        });
        time("qua_gain 102 C", |i| {
            bb(qua_gain_c(bb(&gins[i]), &mut t));
        });
        // MR795_gain_quant on encoder-like inputs (realistic exponents; the lanes accept).
        let gins: Vec<(GainIn, [[i16; 40]; 3])> = (0..64)
            .map(|_| {
                let x = GainIn {
                    exp: std::array::from_fn(|_| 10 + (g.u64() % 12) as i16),
                    exp_gcode0: (g.u64() % 16) as i16,
                    ..gain_in(&mut g)
                };
                let mut v = [[0i16; 40]; 3];
                v.iter_mut().for_each(|s| signal(&mut g, s));
                (x, v)
            })
            .collect();
        time("mr795_gain_quant", |i| {
            let (x, v) = &gins[i];
            let (mut st, mut r, mut ana, mut pos) =
                (GainAdaptState::default(), [12000i16, 0, 0, 0], [0i16; 2], 0);
            let [gp, gc, q1, q2] = &mut r;
            mr795_gain_quant(
                &mut st,
                bb(&v[0]),
                &v[1],
                &v[2],
                bb(&x.frac),
                &x.exp,
                12,
                20000,
                x.exp_gcode0,
                x.frac_gcode0,
                20000,
                12,
                32767,
                gp,
                gc,
                q1,
                q2,
                &mut ana,
                &mut pos,
                &mut ov,
            );
            bb(r);
        });
        time("mr795_gain_quant C", |i| {
            let (x, v) = &gins[i];
            let (mut st, mut r, mut ana, mut o, mut v) = (
                GainAdaptState::default(),
                [12000i16, 0, 0, 0],
                [0i16; 2],
                0,
                *v,
            );
            let mut anap = ana.as_mut_ptr();
            let (mut f, mut e) = (x.frac, x.exp);
            let [gp, gc, q1, q2] = &mut r;
            let [a, b, d] = &mut v;
            unsafe {
                oc_mr795_gain_quant(
                    &mut st,
                    a.as_mut_ptr(),
                    b.as_mut_ptr(),
                    d.as_mut_ptr(),
                    bb(f.as_mut_ptr()),
                    e.as_mut_ptr(),
                    12,
                    20000,
                    x.exp_gcode0,
                    x.frac_gcode0,
                    40,
                    20000,
                    12,
                    32767,
                    gp,
                    gc,
                    q1,
                    q2,
                    &mut anap,
                    &mut t,
                    &mut o,
                )
            };
            bb(r);
        });
        let mut signs = [[0i16; 40]; 64];
        for s in signs.iter_mut().flatten() {
            *s = if g.u64() & 1 == 0 { 32767 } else { -32767 };
        }
        let (mut rr, mut ov, mut o) = ([[0i16; 40]; 40], false, 0);
        let mut rt = RrTm::ZERO;
        time("cor_h_tm", |i| {
            cor_h_tm(bb(&hs[i]), bb(&signs[i]), bb(&mut rt), &mut ov)
        });
        // MR475/MR515's calls: blocks of different residues (sub_nr 1), of equal ones (sub_nr 0)
        time("cor_h_blocks", |i| {
            let b = [[[3, 0], [4, 2]], [[2, 0], [3, 1]]][i & 1];
            crate::simd::cor_h_blocks(bb(&hs[i]), bb(&signs[i]), bb(&mut rt), &b, &mut ov)
        });
        fn c_search<const T: usize, const P: usize>(
            (dn, rr, ip, pm): &SearchIn<T, P>,
            o: &mut CFlag,
        ) -> [i16; P] {
            let (mut dn, mut rr, mut ip, mut pm, mut cv) = (*dn, *rr, *ip, *pm, [0i16; P]);
            unsafe {
                c_search_10and8i40(
                    P as i16,
                    T as i16,
                    T as i16,
                    dn.as_mut_ptr(),
                    rr.as_mut_ptr(),
                    ip.as_mut_ptr(),
                    pm.as_mut_ptr(),
                    cv.as_mut_ptr(),
                    o,
                )
            };
            cv
        }
        let sins: Vec<_> = (0..64)
            .map(|_| {
                let s = search_in::<5, 10>(&mut g, false);
                let t = tm(&s.1);
                (s, t)
            })
            .collect();
        time("search_10i40", |i| {
            let ((dn, _, ipos, pm), rt) = bb(&sins[i]);
            let (mut ip, mut cv) = (*ipos, [0usize; 10]);
            search_10i40(dn, rt, &mut ip, pm, &mut cv);
            bb(cv);
        });
        time("search_10i40 C", |i| {
            bb(c_search(&bb(&sins[i]).0, &mut o));
        });
        // The step-5 searches (MR67, MR74/MR795) on the same kind of rr, dn2 as set_sign marks
        // it (no C entry point of their own: `cbsearch MR67/MR74` time them with cor_h).
        let s5: Vec<_> = sins
            .iter()
            .map(|((dn, _, _, _), rt)| {
                let (mut d, mut sign, mut dn2) = (*dn, [0i16; 40], [0i16; 40]);
                crate::encoder::codebook::set_sign(&mut d, &mut sign, &mut dn2, 6);
                (d, dn2, rt)
            })
            .collect();
        time("search_3i40", |i| {
            let (dn, dn2, rt) = bb(&s5[i]);
            bb(crate::simd::search_3i40(dn, dn2, rt));
        });
        time("search_4i40", |i| {
            let (dn, dn2, rt) = bb(&s5[i]);
            bb(crate::simd::search_4i40(dn, dn2, rt));
        });
        // MR102's cor_h + search (one kernel) on the search inputs' h and signs.
        let sins: Vec<_> = (0..64)
            .map(|_| {
                let (h, sign) = (hs[g.u64() as usize % 64], signs[g.u64() as usize % 64]);
                let (dn, _, ipos, pm) = search_in::<4, 8>(&mut g, false);
                (h, sign, dn, ipos, pm)
            })
            .collect();
        time("cor_h_search_8i40", |i| {
            let (h, sign, dn, ipos, pm) = bb(&sins[i]);
            let (mut ip, mut cv) = (*ipos, [0usize; 8]);
            cor_h_search_8i40(h, sign, dn, &mut ip, pm, &mut cv, &mut ov);
            bb(cv);
        });
        time("cor_h_search_8i40 C", |i| {
            let (mut h, mut sign, dn, ipos, pm) = *bb(&sins[i]);
            let mut rr = [[0i16; 40]; 40];
            unsafe { oc_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), rr.as_mut_ptr(), &mut o) };
            bb(c_search(&(dn, rr, ipos, pm), &mut o));
        });
        time("cor_h C", |i| unsafe {
            oc_cor_h(
                bb(hs[i].as_ptr() as *mut i16),
                bb(signs[i].as_ptr() as *mut i16),
                bb(rr.as_mut_ptr()),
                &mut o,
            )
        });
    }
}

/// `Pre_Process` vs the C, from reset and from random states (any `lo`, as the one-word form
/// holds for every integer lo) on full-scale and speech-like input. Classes counted: outputs
/// whose accumulator truncates y1_hi (|l| >= 2^27; the form's sext28).
mod pre_process {
    use super::super::{PreProcessState, pre_process};
    use super::refc::c_pre_process;
    use super::*;

    #[repr(C)]
    #[derive(Debug, PartialEq, Eq)]
    struct CState([i16; 6]);

    unsafe extern "C" {}

    fn c_of(s: &PreProcessState) -> CState {
        CState([s.y2_hi, s.y2_lo, s.y1_hi, s.y1_lo, s.x0, s.x1])
    }

    #[test]
    fn pre_process_vs_c() {
        let mut g = Rng::new(77);
        let mut trunc = 0;
        for seq in 0..N / 50 {
            let mut st = if seq.is_multiple_of(2) {
                PreProcessState::default()
            } else {
                PreProcessState {
                    y2_hi: g.i16(),
                    y2_lo: g.i16(),
                    y1_hi: g.i16(),
                    y1_lo: g.i16(),
                    x0: g.i16(),
                    x1: g.i16(),
                }
            };
            let mut cst = c_of(&st);
            for call in 0..8 {
                // a length-1 first call on a random state: one input pair is kept as is
                let lg = if call == 0 && seq % 4 == 1 {
                    1
                } else {
                    [L_FRAME, L_FRAME, 1, 2, 3, 40][call % 6]
                };
                let mut sig = [0i16; L_FRAME];
                if seq.is_multiple_of(3) {
                    // full-scale square-ish waves drive the filter state to its range
                    let a = if g.u64().is_multiple_of(2) {
                        i16::MAX
                    } else {
                        i16::MIN
                    };
                    let p = 1 + (g.u64() % 40) as usize;
                    for (i, v) in sig.iter_mut().enumerate() {
                        *v = if (i / p).is_multiple_of(2) { a } else { !a };
                    }
                } else {
                    signal(&mut g, &mut sig);
                }
                let mut csig = sig;
                // count truncating accumulators with the C's own two-word recursion
                let mut s2 = st.clone();
                for i in 0..lg {
                    let l = s2.y1_hi as i32 * 7807
                        + ((s2.y1_lo as i32 * 7807) >> 15)
                        + s2.y2_hi as i32 * -3733
                        + ((s2.y2_lo as i32 * -3733) >> 15)
                        + s2.x1 as i32 * 1899
                        + s2.x0 as i32 * -3798
                        + sig[i] as i32 * 1899;
                    trunc += (l >> 12 != (l >> 12) as i16 as i32) as u32;
                    (s2.y2_hi, s2.y2_lo) = (s2.y1_hi, s2.y1_lo);
                    s2.y1_hi = (l >> 12) as i16;
                    s2.y1_lo = l.wrapping_shl(3).wrapping_sub((s2.y1_hi as i32) << 15) as i16;
                    (s2.x1, s2.x0) = (s2.x0, sig[i]);
                }
                pre_process(&mut st, &mut sig[..lg]);
                unsafe { c_pre_process(&mut cst, csig.as_mut_ptr(), lg as i16) };
                assert_eq!((sig, &c_of(&st)), (csig, &cst), "seq {seq} call {call}");
            }
        }
        assert!(trunc > 100, "truncating class: {trunc}");
    }
}
