//! AMR-WB decoder differential: the Rust `amr_wb::Decoder` against the TS 26.173 reference
//! decoder (`ref26173_sys::Decoder`, its `decoder.c` frame loop with homing), frame by frame.
//! Input per frame: one byte (RX frame type `% 8`, mode `>> 3 % 9`) and 60 bytes of serial
//! bits (477 used; bits past the frame's size are zero, as a storage reader leaves them).
//! Any PCM difference is a crash.
//!
//! The one deliberate difference (CR 26.173-0027..29, NO_DATA during speech; see the amr-wb
//! crate docs) is kept out by construction, as in `crates/amr-wb/tests/common/stress.rs`:
//! a NO_DATA frame outside DTX (no SID frame since the last good speech frame) becomes
//! SPEECH_LOST.
#![no_main]

use amr_wb::{Decoder, L_FRAME, MAX_SERIAL_SIZE, Mode, RxFrameType, SID_BITS};
use libfuzzer_sys::fuzz_target;

const TYPES: [RxFrameType; 8] = [
    RxFrameType::SpeechGood,
    RxFrameType::SpeechProbablyDegraded,
    RxFrameType::SpeechLost,
    RxFrameType::SpeechBad,
    RxFrameType::SidFirst,
    RxFrameType::SidUpdate,
    RxFrameType::SidBad,
    RxFrameType::NoData,
];

fuzz_target!(|data: &[u8]| {
    use RxFrameType::*;
    let mut rust = Decoder::new();
    let mut c = ref26173_sys::Decoder::new();
    let mut in_dtx = false;
    for (i, chunk) in data.as_chunks::<61>().0.iter().enumerate() {
        let mut ft = TYPES[usize::from(chunk[0] % 8)];
        let mode = Mode::from_index(usize::from(chunk[0] >> 3) % 9).unwrap();
        if ft == NoData && !in_dtx {
            ft = SpeechLost;
        }
        in_dtx = match ft {
            SidFirst | SidUpdate | SidBad => true,
            NoData | SpeechBad | SpeechLost => in_dtx,
            SpeechGood | SpeechProbablyDegraded => false,
        };
        let nbits = match ft {
            SidFirst | SidUpdate | SidBad => SID_BITS,
            NoData | SpeechLost => 0,
            _ => mode.bits(),
        };
        let bits: [i16; MAX_SERIAL_SIZE] = std::array::from_fn(|k| {
            if k < nbits {
                i16::from(chunk[1 + k / 8] >> (k % 8) & 1)
            } else {
                0
            }
        });
        let mut r = [0i16; L_FRAME];
        rust.decode(ft, mode, &bits, &mut r);
        let mut prms = bits.map(|b| if b != 0 { 127 } else { -127 });
        let mut cpcm = [0i16; L_FRAME];
        c.decode(mode as i16, ft as i16, &mut prms, &mut cpcm);
        assert_eq!(r, cpcm, "frame {i}: {ft:?} {mode:?}");
    }
});
