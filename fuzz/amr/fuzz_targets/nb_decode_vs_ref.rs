//! Decoder differential: the Rust decoder against the TS 26.073 reference decoder
//! (`Speech_Decode_Frame`, as its `decoder.c` drives it, from `ref26073-sys`), frame by frame.
//! Input per frame: one byte (RX frame type `% 8` of the eight, mode `>> 3 % 8`) and 31 bytes
//! of serial bits (244 used). Every frame type is reachable, bad and NO_DATA frames included;
//! any PCM difference is a crash.
//!
//! The reference's homing is in `decoder.c`, outside `Speech_Decode_Frame`: comparing stops at
//! the first frame that would home a freshly reset Rust decoder (output 160 × 0x0008).
#![no_main]

use amr_nb::{Decoder, MAX_SERIAL_SIZE, Mode, RxFrameType};
use libfuzzer_sys::fuzz_target;
use std::ffi::c_void;
use std::ptr::null_mut;

const MODES: [Mode; 8] = [
    Mode::MR475,
    Mode::MR515,
    Mode::MR59,
    Mode::MR67,
    Mode::MR74,
    Mode::MR795,
    Mode::MR102,
    Mode::MR122,
];

unsafe extern "C" {
    fn ref_Speech_Decode_Frame_init(st: *mut *mut c_void, id: *const u8) -> i32;
    fn ref_Speech_Decode_Frame_exit(st: *mut *mut c_void);
    fn ref_Speech_Decode_Frame(
        st: *mut c_void,
        mode: i32,
        serial: *mut i16,
        frame_type: i32,
        synth: *mut i16,
    ) -> i32;
}

const TYPES: [RxFrameType; 8] = [
    RxFrameType::SpeechGood,
    RxFrameType::SpeechDegraded,
    RxFrameType::Onset,
    RxFrameType::SpeechBad,
    RxFrameType::SidFirst,
    RxFrameType::SidUpdate,
    RxFrameType::SidBad,
    RxFrameType::NoData,
];

/// The reference decoder's state (freed on drop, also when an assertion fails).
struct Ref(*mut c_void);

impl Drop for Ref {
    fn drop(&mut self) {
        unsafe { ref_Speech_Decode_Frame_exit(&mut self.0) }
    }
}

fuzz_target!(|data: &[u8]| {
    let _l = ref26073_sys::lock();
    let mut c = Ref(null_mut());
    unsafe { ref_Speech_Decode_Frame_init(&mut c.0, c"Decoder".as_ptr().cast()) };
    let mut rust = Decoder::new();
    // decoder.c passes its previous mode with NO_DATA; the reference never initialises it
    // (sp_dec.h), the Rust decoder starts from MR475 (`Decoder::new`): the harness follows.
    let mut prev_mode = Mode::MR475;
    for (i, chunk) in data.as_chunks::<32>().0.iter().enumerate() {
        let ft = TYPES[usize::from(chunk[0] % 8)];
        let mode = MODES[usize::from(chunk[0] >> 3 & 7)];
        let mut bits: [i16; MAX_SERIAL_SIZE] =
            std::array::from_fn(|k| i16::from(chunk[1 + k / 8] >> (k % 8) & 1));
        let mut probe = [0i16; 160];
        Decoder::new().decode(ft, mode, &bits, &mut probe);
        if probe == [0x0008; 160] {
            break; // a homing frame: decoder.c would reset the reference
        }
        let mut r = [0i16; 160];
        rust.decode(ft, mode, &bits, &mut r);
        let m = if ft == RxFrameType::NoData {
            prev_mode
        } else {
            mode
        };
        prev_mode = m;
        let mut cpcm = [0i16; 160];
        unsafe {
            ref_Speech_Decode_Frame(
                c.0,
                m as i32,
                bits.as_mut_ptr(),
                ft as i32,
                cpcm.as_mut_ptr(),
            )
        };
        assert_eq!(r, cpcm, "frame {i}: {ft:?} {mode:?}");
    }
});
