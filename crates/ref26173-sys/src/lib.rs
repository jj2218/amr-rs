//! The 3GPP TS 26.173 AMR-WB reference C code, the unit-test oracle: every symbol it defines is
//! prefixed `refwb_` (see build.rs). Tests declare the internal functions they call themselves,
//! with the reference signatures (`Word16 = i16`, `Word32 = i32`, `Flag = i32`); the decoder's
//! entry points are declared here.
//!
//! The reference keeps its overflow flag in one global (`refwb_Overflow`), shared by every
//! thread: hold [`lock`] while calling into it (tests run in parallel).
//!
//! It includes CR 26.173-0027..29 (`dtx_vad_hist`: a NO_DATA frame during speech is DTX only
//! after 8 frames without VAD), which the Rust port deliberately leaves out (a live patent).
//! Differential tests must not send NO_DATA frames that the rule would treat differently.
#![allow(unsafe_code, non_upper_case_globals)]

use std::ffi::c_void;
use std::sync::{Mutex, MutexGuard};

unsafe extern "C" {
    /// The reference's global overflow flag (`Flag Overflow` in basicop2.c).
    pub static mut refwb_Overflow: i32;

    /// `dec_main.c`: allocates and resets a decoder state.
    pub fn refwb_Init_decoder(state: *mut *mut c_void);
    pub fn refwb_Close_decoder(state: *mut c_void);
    pub fn refwb_Reset_decoder(state: *mut c_void, reset_all: i16);
    /// One frame: `prms` is the serial bit array (±127 per bit), `frame_type` `RX_*`.
    pub fn refwb_decoder(
        mode: i16,
        prms: *mut i16,
        synth16k: *mut i16,
        frame_length: *mut i16,
        state: *mut c_void,
        frame_type: i16,
    );
    /// `homing.c`
    pub fn refwb_decoder_homing_frame_test(prms: *mut i16, mode: i16) -> i16;
    pub fn refwb_decoder_homing_frame_test_first(prms: *mut i16, mode: i16) -> i16;
}

static LOCK: Mutex<()> = Mutex::new(());

/// Serialises use of the reference code (its global flag); a panicking holder doesn't poison it.
pub fn lock() -> MutexGuard<'static, ()> {
    LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

/// Reads the global overflow flag (hold [`lock`]).
pub fn overflow() -> bool {
    // SAFETY: a plain int the C code reads and writes; callers hold `lock`.
    unsafe { refwb_Overflow != 0 }
}

/// Clears the global overflow flag (hold [`lock`]).
pub fn clear_overflow() {
    // SAFETY: as `overflow`.
    unsafe { refwb_Overflow = 0 }
}

/// The reference decoder's frame loop (`decoder.c` `main`) as it runs on `.awb` (MIME) and
/// ITU input: SID, NO_DATA and lost frames take the last mode (`Read_serial`'s `prev_mode`),
/// homing per TS 26.173, the 2 LSBs cleared (14-bit output).
pub struct Decoder {
    st: *mut c_void,
    mode_old: i16,
    reset_flag: i16,
    reset_flag_old: i16,
}

/// `RX_NO_DATA` and `RX_SPEECH_LOST` (`dtx.h`): frames without a mode of their own.
const RX_SPEECH_LOST: i16 = 2;
const RX_NO_DATA: i16 = 7;
/// `RX_SID_FIRST`..=`RX_SID_BAD`.
const RX_SID_FIRST: i16 = 4;
const RX_SID_BAD: i16 = 6;

impl Decoder {
    pub fn new() -> Self {
        let mut st = std::ptr::null_mut();
        // SAFETY: Init_decoder allocates and writes the state pointer.
        unsafe { refwb_Init_decoder(&mut st) };
        Self {
            st,
            mode_old: 0,
            reset_flag: 0,
            reset_flag_old: 1,
        }
    }

    /// Decodes one frame; `prms` holds the frame's bits (35 for SID frames, none read for
    /// NO_DATA), at least 477 words (the reference reads up to the largest mode).
    pub fn decode(
        &mut self,
        mut mode: i16,
        frame_type: i16,
        prms: &mut [i16],
        out: &mut [i16; 320],
    ) {
        assert!(prms.len() >= 477);
        let _g = lock();
        let p = prms.as_mut_ptr();
        // SAFETY: `prms` holds every bit the reference reads, `out` 320 samples, `st` is live.
        unsafe {
            if frame_type == RX_NO_DATA || frame_type == RX_SPEECH_LOST {
                mode = self.mode_old;
                self.reset_flag = 0;
            } else {
                if (RX_SID_FIRST..=RX_SID_BAD).contains(&frame_type) {
                    mode = self.mode_old; // Read_serial: `*mode = st->prev_mode`
                }
                self.mode_old = mode;
                if self.reset_flag_old == 1 {
                    self.reset_flag = refwb_decoder_homing_frame_test_first(p, mode);
                }
            }
            if self.reset_flag != 0 && self.reset_flag_old != 0 {
                out.fill(0x0008); // EHF_MASK
            } else {
                let mut len = 0;
                refwb_decoder(mode, p, out.as_mut_ptr(), &mut len, self.st, frame_type);
                for s in out.iter_mut() {
                    *s &= !3;
                }
            }
            if self.reset_flag_old == 0 {
                self.reset_flag = refwb_decoder_homing_frame_test(p, mode);
            }
            if self.reset_flag != 0 {
                refwb_Reset_decoder(self.st, 1);
            }
        }
        self.reset_flag_old = self.reset_flag;
    }
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for Decoder {
    fn drop(&mut self) {
        // SAFETY: `st` came from Init_decoder and is freed once.
        unsafe { refwb_Close_decoder(self.st) };
    }
}
