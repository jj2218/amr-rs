//! The encoder's differential tests against the TS 26.073 reference (`ref26073-sys`): the whole
//! encoder and `Pre_Process`, and the helpers each module's own `vs_ref` shares (signals, impulse
//! responses, opencore's tables, the reference in opencore's function shapes: `refc`). opencore
//! where a test says so, and as `kern`'s speed baseline.
#![allow(unsafe_code)]

use crate::cnst::{CodecMode, L_FRAME, M};
use crate::encoder::gain_q::GainAdaptState;
use crate::testutil::{N, Rng};
use opencore_sys as _;
use std::ffi::c_void;

pub(crate) type CFlag = i32;

/// A signal of random amplitude (`>> 0..16`), sprinkled with edge values.
pub(crate) fn signal(g: &mut Rng, buf: &mut [i16]) {
    let scale = (g.u64() % 16) as u32;
    for v in buf {
        *v = if g.u64().is_multiple_of(32) {
            g.i16()
        } else {
            (g.u64() as i16) >> scale
        };
    }
}

pub(crate) fn mode_of(g: &mut Rng) -> CodecMode {
    CodecMode::ALL[(g.u64() % 8) as usize]
}

const _: usize = M;

/// The TS 26.073 reference (`ref26073-sys`) behind the opencore-shaped signatures these tests
/// were written against: the tables argument is unused, `o` receives the reference's global
/// overflow flag (cleared before the call). Each call holds the reference lock.
#[allow(clippy::too_many_arguments)]
pub(crate) mod refc {
    use super::CFlag;
    use ref26073_sys::{clear_overflow, lock, overflow};

    unsafe extern "C" {
        pub(crate) fn ref_cor_h(h: *mut i16, sign: *mut i16, rr: *mut [i16; 40]);
        pub(crate) fn ref_Pre_Process(st: *mut i16, signal: *mut i16, lg: i16) -> i32;
        pub(crate) fn ref_cbsearch(
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
    pub(crate) unsafe fn call<T>(o: *mut CFlag, f: impl FnOnce() -> T) -> T {
        let _l = lock();
        clear_overflow();
        let r = f();
        if !o.is_null() {
            unsafe { *o = overflow() as CFlag };
        }
        r
    }

    /// `Pre_Process` (its state: 6 Word16s, as `Pre_ProcessState`).
    pub(crate) unsafe fn c_pre_process<S>(st: *mut S, signal: *mut i16, lg: i16) {
        unsafe {
            call(std::ptr::null_mut(), || {
                ref_Pre_Process(st.cast(), signal, lg)
            })
        };
    }
    pub(crate) unsafe fn c_cor_h(h: *mut i16, sign: *mut i16, rr: *mut [i16; 40], o: *mut CFlag) {
        unsafe { call(o, || ref_cor_h(h, sign, rr)) }
    }
    pub(crate) unsafe fn c_cbsearch(
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
        _t: *mut super::CommonAmrTbls,
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
pub(crate) mod whole {
    use super::*;
    use crate::cnst::L_FRAME;
    use crate::decoder::MAX_SERIAL_SIZE;
    use crate::encoder::{Dtx, Encoder};
    use std::ffi::c_void;
    use std::ptr::null_mut;

    unsafe extern "C" {
        pub(crate) fn ref_Speech_Encode_Frame_init(
            st: *mut *mut c_void,
            dtx: i32,
            id: *const u8,
        ) -> i32;
        pub(crate) fn ref_Speech_Encode_Frame_exit(st: *mut *mut c_void);
        pub(crate) fn ref_Speech_Encode_Frame(
            st: *mut c_void,
            mode: i32,
            speech: *mut i16,
            serial: *mut i16,
            used_mode: *mut i32,
        ) -> i32;
        pub(crate) fn ref_sid_sync_init(st: *mut *mut c_void) -> i32;
        pub(crate) fn ref_sid_sync_exit(st: *mut *mut c_void);
        pub(crate) fn ref_sid_sync(st: *mut c_void, mode: i32, tx_type: *mut i32);
    }

    /// Signal generator: per segment a random mix of tone, noise, silence and clipping.
    pub(crate) struct Gen {
        pub(crate) g: Rng,
        pub(crate) phase: f64,
        pub(crate) freq: f64,
        pub(crate) amp: f64,
        pub(crate) noise: i32,
        pub(crate) clip: bool,
        pub(crate) left: usize,
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

    pub(crate) fn run(seed: u64, dtx: bool, frames: usize) {
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
            let r = rust.encode(crate::Mode::ALL[mode as usize], &pcm);
            let tx = match r.frame_type {
                crate::RxFrameType::SpeechGood => 0,
                crate::RxFrameType::SidFirst => 1,
                crate::RxFrameType::SidUpdate => 2,
                _ => 3,
            };
            assert_eq!(
                tx, ets[0],
                "seed {seed} dtx {dtx} frame {f} {mode:?}: tx type"
            );
            assert_eq!(
                r.bits[..],
                ets[1..=MAX_SERIAL_SIZE],
                "seed {seed} dtx {dtx} frame {f} {mode:?}: bits"
            );
        }
        unsafe {
            ref_Speech_Encode_Frame_exit(&mut enc);
            ref_sid_sync_exit(&mut sid);
        }
    }

    #[test]
    pub(crate) fn encoder_matches_ref_without_dtx() {
        for seed in 0..40 {
            run(seed, false, 500);
        }
    }

    #[test]
    pub(crate) fn encoder_matches_ref_with_dtx() {
        for seed in 100..140 {
            run(seed, true, 500);
        }
    }
}

unsafe extern "C" {
    #[link_name = "cbsearch"]
    pub(crate) fn oc_cbsearch(
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

#[repr(C)]
pub(crate) struct CommonAmrTbls {
    pub(crate) ptrs: [*const c_void; 34],
}

unsafe extern "C" {
    pub(crate) fn get_const_tbls(t: *mut c_void);
}

pub(crate) fn tables() -> CommonAmrTbls {
    let mut t = CommonAmrTbls {
        ptrs: [std::ptr::null(); 34],
    };
    unsafe { get_const_tbls((&raw mut t).cast()) };
    t
}

// opencore, for `kernel_timing` (`BENCH=kern`): the speed baseline stays the C we beat.
unsafe extern "C" {
    #[link_name = "Convolve"]
    pub(crate) fn oc_convolve(x: *mut i16, h: *mut i16, y: *mut i16, l: i16);
    #[link_name = "cor_h"]
    pub(crate) fn oc_cor_h(h: *mut i16, sign: *mut i16, rr: *mut [i16; 40], o: *mut CFlag);
    #[link_name = "MR795_gain_quant"]
    pub(crate) fn oc_mr795_gain_quant(
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

unsafe extern "C" {
    #[link_name = "search_10and8i40"]
    pub(crate) fn c_search_10and8i40(
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
    #[link_name = "Syn_filt"]
    pub(crate) fn c_syn_filt(
        a: *mut i16,
        x: *mut i16,
        y: *mut i16,
        lg: i16,
        mem: *mut i16,
        update: i16,
    );
}

/// An impulse response: random, decaying at a random rate, sometimes full scale.
pub(crate) fn impulse(g: &mut Rng, h: &mut [i16]) {
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

/// Pulse signs: ±32767 as `set_sign12k2` makes them, sometimes arbitrary.
pub(crate) fn signs(g: &mut Rng, sign: &mut [i16; 40]) {
    let any = g.u64().is_multiple_of(4);
    for s in sign {
        *s = match (any, g.u64() & 1) {
            (true, _) => g.i16(),
            (false, 0) => 32767,
            _ => -32767,
        };
    }
}

/// `Pre_Process` vs the C, from reset and from random states (any `lo`, as the one-word form
/// holds for every integer lo) on full-scale and speech-like input. Classes counted: outputs
/// whose accumulator truncates y1_hi (|l| >= 2^27; the form's sext28).
pub(crate) mod pre_process {
    use super::super::{PreProcessState, pre_process};
    use super::refc::c_pre_process;
    use super::*;

    #[repr(C)]
    #[derive(Debug, PartialEq, Eq)]
    pub(crate) struct CState([i16; 6]);

    fn c_of(s: &PreProcessState) -> CState {
        CState([s.y2_hi, s.y2_lo, s.y1_hi, s.y1_lo, s.x0, s.x1])
    }

    #[test]
    pub(crate) fn pre_process_vs_ref() {
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
