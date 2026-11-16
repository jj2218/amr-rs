//! AMR-WB encoder differential: the Rust `amr_wb::Encoder` against the TS 26.173 reference
//! encoder (`ref26173_sys::Encoder`, its `coder.c` frame loop with homing and `Write_serial`'s
//! frame types), frame by frame, on arbitrary PCM. Input: one byte of flags (bit 0: DTX; bits
//! 1..=3: each frame is encoded 1 + that many times in a row; bits 4..=7: the samples are
//! shifted right by that much), then per frame one mode byte (`% 9`) and 320 little-endian
//! samples. Repeats and quiet frames let inputs of a few frames reach the end of the DTX
//! hangover (8 frames without voice) and the SID frames. The TX frame type and the serial bits
//! must be equal; any difference is a crash.
#![no_main]

use amr_wb::{Encoder, L_FRAME16K, Mode, RxFrameType, SID_BITS};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Some((&flags, frames)) = data.split_first() else {
        return;
    };
    let dtx = flags & 1 != 0;
    let reps = 1 + usize::from((flags >> 1) & 7);
    let shift = flags >> 4;
    let mut rust = Encoder::new(dtx);
    let mut c = ref26173_sys::Encoder::new(dtx);
    let mut want = [0i16; 480];
    for (i, chunk) in frames
        .as_chunks::<{ 1 + 2 * L_FRAME16K }>()
        .0
        .iter()
        .enumerate()
    {
        let m = chunk[0] % 9;
        let pcm: [i16; L_FRAME16K] = std::array::from_fn(|k| {
            i16::from_le_bytes([chunk[1 + 2 * k], chunk[2 + 2 * k]]) >> shift
        });
        for r in 0..reps {
            let w = c.encode(i16::from(m), &pcm, &mut want);
            let f = rust.encode(Mode::from_index(usize::from(m)).unwrap(), &pcm);
            let (tx, n) = match f.frame_type {
                RxFrameType::SpeechGood => (0, f.mode.bits()),
                RxFrameType::SidFirst => (1, SID_BITS),
                RxFrameType::SidUpdate => (2, SID_BITS),
                _ => (3, SID_BITS),
            };
            let got: Vec<i16> = [0x6b21, tx, i16::from(m)]
                .into_iter()
                .chain(f.bits[..n].iter().map(|&b| if b != 0 { 127 } else { -127 }))
                .collect();
            assert_eq!(got[..], want[..w], "frame {i} (rep {r}): mode {m}");
        }
    }
});
