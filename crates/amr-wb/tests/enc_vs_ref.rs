//! Whole-encoder differential test against the TS 26.173 reference encoder on stress input:
//! real speech (the vector inputs) scaled up to clipping, white noise at every level, held
//! tones and full-scale square waves (resonant LPC, peaky impulse responses), silence and
//! near-silence (the DTX paths), encoder homing frames, per-frame mode switching, DTX on and
//! off. Every output word is compared.
mod common;

use amr_wb::{Encoder, Mode};
use common::stress::Rng;
use common::{cod_words, encoder_cases};

/// An input class with its parameters, held for a run of frames.
#[derive(Clone, Copy, Debug)]
enum Class {
    /// Real speech times a gain (saturating).
    Speech(i32),
    /// White noise of amplitude 2^k.
    Noise(u32),
    /// A held tone (Hz, amplitude), optionally sweeping.
    Tone(f64, f64, f64),
    /// A full-scale square wave of half-period n samples.
    Square(usize),
    /// Near-silence of amplitude 2^k (0: digital silence).
    Quiet(u32),
    /// The encoder homing frame.
    Homing,
    /// Speech with bursts of random samples.
    Bursts,
}

fn pick(r: &mut Rng) -> Class {
    match r.below(7) {
        0 => Class::Speech([1, 2, 4, 16][r.below(4) as usize]),
        1 => Class::Noise(r.below(16) as u32),
        2 => Class::Tone(
            50.0 + r.below(7900) as f64,
            [32767.0, 20000.0, 3000.0, 100.0][r.below(4) as usize],
            if r.chance(30) {
                r.below(2000) as f64
            } else {
                0.0
            },
        ),
        3 => Class::Square(1 + r.below(200) as usize),
        4 => Class::Quiet(r.below(4) as u32),
        5 => Class::Homing,
        _ => Class::Bursts,
    }
}

/// One 320-sample frame of `class` at time `t` (samples; advanced).
fn frame(r: &mut Rng, class: Class, speech: &[i16], t: &mut usize) -> [i16; 320] {
    let mut f = [0i16; 320];
    let start = *t;
    *t += 320;
    for (i, v) in f.iter_mut().enumerate() {
        let n = start + i;
        *v = match class {
            Class::Speech(g) => (speech[n % speech.len()] as i32 * g).clamp(-32768, 32767) as i16,
            Class::Noise(k) => ((r.next() as i16) as i32 >> (15 - k)) as i16,
            Class::Tone(fr, a, sweep) => {
                let s = n as f64 / 16000.0;
                (a * (std::f64::consts::TAU * (fr + sweep * s) * s).sin()) as i16
            }
            Class::Square(half) => {
                if (n / half).is_multiple_of(2) {
                    32767
                } else {
                    -32768
                }
            }
            Class::Quiet(k) => (r.below(1 << k) as i16) - (1 << k) / 2,
            Class::Homing => 8,
            Class::Bursts => {
                if r.chance(5) {
                    r.next() as i16
                } else {
                    speech[n % speech.len()]
                }
            }
        };
    }
    f
}

/// Runs `frames` frames of a stream (classes held for runs of frames, modes switching)
/// through both encoders; returns the first differing frame.
fn stream(seed: u64, frames: usize, dtx: bool, speech: &[i16]) -> Option<String> {
    let mut r = Rng::seeded(seed);
    let mut rust = Encoder::new(dtx);
    let mut c = ref26173_sys::Encoder::new(dtx);
    let mut t = 0;
    let mut class = pick(&mut r);
    let mut mode = r.below(9) as u8;
    let mut want = [0i16; 480];
    for n in 0..frames {
        if r.chance(8) {
            class = pick(&mut r);
        }
        if r.chance(10) {
            mode = r.below(9) as u8;
        }
        let f = frame(&mut r, class, speech, &mut t);
        let w = c.encode(mode as i16, &f, &mut want);
        let got = cod_words(&rust.encode(Mode::from_index(mode as usize).unwrap(), &f));
        if got[..] != want[..w] {
            return Some(format!(
                "seed {seed} dtx {dtx}: frame {n} ({class:?}, mode {mode}) differs"
            ));
        }
    }
    None
}

#[test]
fn stress_streams_match_reference() {
    let speech: Vec<i16> = encoder_cases()
        .into_iter()
        .filter(|c| c.name.contains("T00") || c.name.contains("T05"))
        .flat_map(|c| c.input)
        .collect();
    assert!(!speech.is_empty());
    const N: u64 = 64;
    let mut failing = Vec::new();
    for seed in 0..N {
        if let Some(f) = stream(seed, 300, seed % 2 == 0, &speech) {
            failing.push(f);
        }
    }
    println!(
        "encoder stress streams: {}/{N} bit-exact vs the reference",
        N as usize - failing.len()
    );
    assert!(failing.is_empty(), "{failing:#?}");
}

/// Speech with pauses, 600 frames: 30 frames of T00's speech (repeated), then 30 of silence
/// (the DTX periods), with the mode changing every 60 frames (each mode).
fn speech_with_pauses() -> Vec<(Mode, [i16; 320])> {
    let speech: Vec<i16> = encoder_cases()
        .into_iter()
        .find(|c| c.name.contains("T00"))
        .unwrap()
        .input;
    let chunks = speech.as_chunks::<320>().0;
    (0..600)
        .map(|n| {
            let pcm = if (n / 30) % 2 == 0 {
                chunks[n % chunks.len()]
            } else {
                [0; 320]
            };
            (Mode::ALL[(n / 60) % 9], pcm)
        })
        .collect()
}

/// The encoder's output through `.awb` (RFC 4867 storage) and back: speech, SID_FIRST,
/// SID_UPDATE and NO_DATA frames (DTX on, speech with pauses) keep their type and bits, but
/// SID_FIRST's comfort-noise bits, which are written as zeros.
#[test]
fn encoder_output_round_trips_through_awb() {
    use amr_wb::rfc4867::{Frame, read_storage, write_storage};
    use amr_wb::{RxFrameType, SID_BITS};
    let mut enc = Encoder::new(true);
    let mut frames = Vec::new();
    for (mode, pcm) in speech_with_pauses() {
        let f = enc.encode(mode, &pcm);
        // what is sent: nothing of a NO_DATA frame
        frames.push(if f.frame_type == RxFrameType::NoData {
            Frame::no_data()
        } else {
            f
        });
    }
    let kinds: std::collections::HashSet<_> = frames.iter().map(|f| f.frame_type).collect();
    assert_eq!(kinds.len(), 4, "every frame type: {kinds:?}");
    let back: Vec<Frame> = read_storage(&write_storage(&frames).unwrap())
        .unwrap()
        .collect::<Result<_, _>>()
        .unwrap();
    for f in &mut frames {
        if f.frame_type == RxFrameType::SidFirst {
            f.bits[..SID_BITS].fill(0);
        }
    }
    assert_eq!(back, frames);
}

/// Our `.awb` writer vs the TS 26.173 reference's (`Write_serial`, format 2 = MIME, as its
/// `coder -mime` writes): the reference encoder's frames through its writer, ours through
/// `write_storage`, byte-equal (DTX on, speech with pauses, every mode; SID_FIRST frames
/// whose comfort-noise analysis is not zero, which both send as zeros).
#[test]
#[allow(unsafe_code)]
fn awb_matches_the_reference_writer() {
    use amr_wb::RxFrameType;
    use amr_wb::rfc4867::{STORAGE_MAGIC, write_storage};
    use std::ffi::{c_long, c_void};
    unsafe extern "C" {
        fn tmpfile() -> *mut c_void;
        fn ftell(f: *mut c_void) -> c_long;
        fn rewind(f: *mut c_void);
        fn fread(p: *mut c_void, size: usize, n: usize, f: *mut c_void) -> usize;
        fn fclose(f: *mut c_void) -> i32;
        fn refwb_Init_write_serial(st: *mut *mut c_void) -> i16;
        fn refwb_Close_write_serial(st: *mut c_void) -> i16;
        fn refwb_Write_serial(
            fp: *mut c_void,
            prms: *mut i16,
            coding_mode: i16,
            mode: i16,
            st: *mut c_void,
            format: i16,
        );
    }
    const MRDTX: i16 = 9;
    let input = speech_with_pauses();
    let mut enc = Encoder::new(true);
    let frames: Vec<_> = input.iter().map(|(m, pcm)| enc.encode(*m, pcm)).collect();
    let first_cn = frames
        .iter()
        .filter(|f| f.frame_type == RxFrameType::SidFirst && f.bits.contains(&1))
        .count();
    assert!(first_cn > 0, "no SID_FIRST frame with comfort-noise bits");

    let mut c = ref26173_sys::Encoder::new(true);
    let mut words = [0i16; 480];
    // SAFETY: `Write_serial` reads at most 477 bits of `prms`; the state and file are live
    // until closed. It touches no global (no lock needed; `encode` takes it).
    let theirs = unsafe {
        let mut tx = std::ptr::null_mut();
        refwb_Init_write_serial(&mut tx);
        let fp = tmpfile();
        assert!(!fp.is_null());
        for (mode, pcm) in &input {
            let n = c.encode(*mode as i16, pcm, &mut words);
            // `coder`'s mode out: MRDTX for SID and NO_DATA frames (TX type 1..=3).
            let coding_mode = if words[1] == 0 { *mode as i16 } else { MRDTX };
            let mut prms = [0i16; 477];
            prms[..n - 3].copy_from_slice(&words[3..n]);
            refwb_Write_serial(fp, prms.as_mut_ptr(), coding_mode, *mode as i16, tx, 2);
        }
        let len = ftell(fp) as usize;
        rewind(fp);
        let mut out = vec![0u8; len];
        assert_eq!(fread(out.as_mut_ptr().cast(), 1, len, fp), len);
        fclose(fp);
        refwb_Close_write_serial(tx);
        out
    };
    let ours = write_storage(&frames).unwrap();
    assert!(
        ours[STORAGE_MAGIC.len()..] == theirs[..],
        "the .awb bytes differ"
    );
    println!(
        "awb vs reference Write_serial: {} frames byte-equal ({first_cn} SID_FIRST with CN bits)",
        frames.len()
    )
}
