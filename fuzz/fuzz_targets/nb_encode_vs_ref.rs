//! Encoder differential: the Rust encoder against the TS 26.073 reference encoder
//! (`Speech_Encode_Frame` + `sid_sync`, as its `coder.c` drives them, from `ref26073-sys`),
//! frame by frame, on arbitrary PCM. Input: one byte of flags (bit 0: DTX with VAD1), then per
//! frame one mode byte (`% 8`) and 160 little-endian samples. The TX frame type and the
//! serial bits must be equal; any difference is a crash.
//!
//! The reference's homing is in `coder.c`, outside the functions called here: comparing stops
//! after the first encoder homing frame (160 × 0x0008), where the Rust encoder resets.
#![no_main]

use amr_nb::{Dtx, Encoder, MAX_SERIAL_SIZE, Mode, RxFrameType};
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

/// The reference encoder's state (freed on drop, also when an assertion fails).
struct Ref(*mut c_void, *mut c_void);

impl Drop for Ref {
    fn drop(&mut self) {
        unsafe {
            ref_Speech_Encode_Frame_exit(&mut self.0);
            ref_sid_sync_exit(&mut self.1);
        }
    }
}

fuzz_target!(|data: &[u8]| {
    let Some((&flags, frames)) = data.split_first() else {
        return;
    };
    let dtx = flags & 1 != 0;
    let _l = ref26073_sys::lock();
    let mut c = Ref(null_mut(), null_mut());
    unsafe {
        ref_Speech_Encode_Frame_init(&mut c.0, i32::from(dtx), c"encoder".as_ptr().cast());
        ref_sid_sync_init(&mut c.1);
    }
    let mut rust = Encoder::new(if dtx { Dtx::Vad1 } else { Dtx::Off });
    for (i, chunk) in frames.as_chunks::<{ 1 + 320 }>().0.iter().enumerate() {
        let mode = MODES[usize::from(chunk[0] % 8)];
        let pcm: [i16; 160] =
            std::array::from_fn(|k| i16::from_le_bytes([chunk[1 + 2 * k], chunk[2 + 2 * k]]));
        let r = rust.encode(mode, &pcm);
        let mut speech = pcm; // C overwrites its input
        let mut serial = [0i16; MAX_SERIAL_SIZE];
        let (mut used, mut tx) = (0, 0);
        unsafe {
            ref_Speech_Encode_Frame(
                c.0,
                mode as i32,
                speech.as_mut_ptr(),
                serial.as_mut_ptr(),
                &mut used,
            );
            ref_sid_sync(c.1, used, &mut tx);
        }
        let rtx = match r.frame_type {
            RxFrameType::SpeechGood => 0,
            RxFrameType::SidFirst => 1,
            RxFrameType::SidUpdate => 2,
            _ => 3,
        };
        assert_eq!((rtx, r.bits), (tx, serial), "frame {i}: {mode:?}");
        if pcm == [0x0008; 160] {
            break; // the Rust encoder homes here; coder.c would reset the reference
        }
    }
});
