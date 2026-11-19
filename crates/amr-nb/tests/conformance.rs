//! Decoder conformance against the official 3GPP TS 26.074 test vectors (`just vectors`).
#![allow(unsafe_code)] // FFI to the C oracle

mod common;

use common::{
    COD_FRAME, Case, EncCase, L_FRAME, MAX_SERIAL_SIZE, run_all, run_all_enc, run_all_sync,
};
use std::ffi::c_void;
use std::ptr::null_mut;

unsafe extern "C" {
    fn GSMInitDecode(state: *mut *mut c_void, id: *mut i8) -> i16;
    fn GSMDecodeFrameExit(state: *mut *mut c_void);
    fn Speech_Decode_Frame_reset(state: *mut c_void) -> i16;
    fn AMRDecode(
        state: *mut c_void,
        frame_type: i32,
        bits: *mut u8,
        pcm: *mut i16,
        format: i32,
    ) -> i16;
    fn AMREncodeInit(enc: *mut *mut c_void, sid: *mut *mut c_void, dtx: i32) -> i16;
    fn AMREncodeExit(enc: *mut *mut c_void, sid: *mut *mut c_void);
    fn AMREncodeReset(enc: *mut c_void, sid: *mut c_void) -> i16;
    fn AMREncode(
        enc: *mut c_void,
        sid: *mut c_void,
        mode: i32,
        input: *mut i16,
        out: *mut u8,
        frame_type: *mut i32,
        format: i16,
    ) -> i16;
    static prmnofsf: [i16; 8];
    static bitno_MR475: [i16; 17];
    static bitno_MR515: [i16; 19];
    static bitno_MR59: [i16; 19];
    static bitno_MR67: [i16; 19];
    static bitno_MR74: [i16; 19];
    static bitno_MR795: [i16; 23];
    static bitno_MR102: [i16; 39];
    static bitno_MR122: [i16; 57];
}
use opencore_sys as _;

const ETS: i32 = 0; // decoder input format
const AMR_TX_ETS: i16 = 2; // encoder output format

/// Decoder homing frame per mode: the encoder's output for an encoder homing frame
/// (160 × 0x0008) from the reset state (TS 26.073 §homing), plus the number of bits that
/// make up its first subframe.
fn dhf() -> Vec<([i16; MAX_SERIAL_SIZE], usize)> {
    let bitno: [&[i16]; 8] = unsafe {
        [
            &bitno_MR475,
            &bitno_MR515,
            &bitno_MR59,
            &bitno_MR67,
            &bitno_MR74,
            &bitno_MR795,
            &bitno_MR102,
            &bitno_MR122,
        ]
    };
    (0..8)
        .map(|mode| {
            let (mut enc, mut sid) = (null_mut(), null_mut());
            let mut ehf = [0x0008i16; L_FRAME];
            let mut out = [0i16; 1 + MAX_SERIAL_SIZE + 1];
            let mut ft = 0;
            unsafe {
                AMREncodeInit(&mut enc, &mut sid, 0);
                AMREncode(
                    enc,
                    sid,
                    mode,
                    ehf.as_mut_ptr(),
                    out.as_mut_ptr().cast(),
                    &mut ft,
                    AMR_TX_ETS,
                );
                AMREncodeExit(&mut enc, &mut sid);
            }
            let first: i16 = bitno[mode as usize][..unsafe { prmnofsf[mode as usize] } as usize]
                .iter()
                .sum();
            (out[1..=MAX_SERIAL_SIZE].try_into().unwrap(), first as usize)
        })
        .collect()
}

/// C oracle with the reference decoder's homing logic (TS 26.073 `decoder.c`), which opencore
/// lacks. Proves the harness (file parsing, TX→RX mapping, homing) before the Rust decoder uses it.
#[test]
fn c_oracle_with_homing_is_bit_exact() {
    let dhf = dhf();
    run_all("c-oracle", &common::OPENCORE_DEFECTS, |case: &Case| {
        c_decode(case, &dhf)
    });
}

fn c_decode(case: &Case, dhf: &[([i16; MAX_SERIAL_SIZE], usize)]) -> Vec<i16> {
    {
        let mut st = null_mut();
        unsafe { GSMInitDecode(&mut st, null_mut()) };
        let mut out = Vec::with_capacity(case.frames.len() * L_FRAME);
        let (mut reset_old, mut prev_mode) = (true, 0i16);
        for f in &case.frames {
            let mode = if f.rx_type == common::rx::NO_DATA {
                prev_mode
            } else {
                f.mode
            };
            prev_mode = mode;
            let (pattern, first) = &dhf[mode as usize];
            let mut reset = reset_old && f.bits[..*first] == pattern[..*first];
            let mut pcm = [0i16; L_FRAME];
            if reset && reset_old {
                pcm = [0x0008; L_FRAME];
            } else {
                let mut ets = [0i16; MAX_SERIAL_SIZE + 2];
                ets[0] = f.rx_type;
                ets[1..=MAX_SERIAL_SIZE].copy_from_slice(&f.bits);
                ets[MAX_SERIAL_SIZE + 1] = f.mode;
                unsafe { AMRDecode(st, 0, ets.as_mut_ptr().cast(), pcm.as_mut_ptr(), ETS) };
            }
            if !reset_old {
                reset = f.bits == *pattern;
            }
            if reset {
                unsafe { Speech_Decode_Frame_reset(st) };
            }
            reset_old = reset;
            out.extend_from_slice(&pcm);
        }
        unsafe { GSMDecodeFrameExit(&mut st) };
        out
    }
}

fn rx_type(t: i16) -> amr_nb::RxFrameType {
    use amr_nb::RxFrameType::*;
    [
        SpeechGood,
        SpeechDegraded,
        Onset,
        SpeechBad,
        SidFirst,
        SidUpdate,
        SidBad,
        NoData,
    ][t as usize]
}

fn mode(m: i16) -> amr_nb::Mode {
    use amr_nb::Mode::*;
    // NO_DATA frames carry -1; the decoder ignores the mode for them.
    [MR475, MR515, MR59, MR67, MR74, MR795, MR102, MR122]
        .get(m as usize)
        .copied()
        .unwrap_or(MR475)
}

/// The Rust decoder (homing built in) against the official vectors.
#[test]
fn rust_decoder_is_bit_exact() {
    run_all("rust", &[], rust_decode);
}

fn rust_decode(case: &Case) -> Vec<i16> {
    let mut dec = amr_nb::Decoder::new();
    let mut out = Vec::with_capacity(case.frames.len() * L_FRAME);
    let mut pcm = [0i16; L_FRAME];
    for f in &case.frames {
        dec.decode(rx_type(f.rx_type), mode(f.mode), &f.bits, &mut pcm);
        out.extend_from_slice(&pcm);
    }
    out
}

/// The vectors packed per RFC 4867 and read back before decoding: each path must stay
/// bit-exact (the vectors carry speech, SID_FIRST, SID_UPDATE and NO_DATA frames).
#[test]
fn rust_decoder_via_rfc4867_is_bit_exact() {
    use amr_nb::rfc4867::{self, Packing};
    type Path = fn(&[rfc4867::Frame]) -> Vec<rfc4867::Frame>;
    let paths: [(&str, Path); 3] = [
        (".amr storage", |f| {
            rfc4867::read_storage(&rfc4867::write_storage(f).unwrap())
                .unwrap()
                .collect::<Result<_, _>>()
                .unwrap()
        }),
        ("RTP octet-aligned", |f| rtp(f, Packing::OctetAligned)),
        ("RTP bandwidth-efficient", |f| {
            rtp(f, Packing::BandwidthEfficient)
        }),
    ];
    /// Packets of 3 frames (the last may be shorter), no codec mode request.
    fn rtp(frames: &[rfc4867::Frame], p: Packing) -> Vec<rfc4867::Frame> {
        frames
            .chunks(3)
            .flat_map(|c| {
                let pkt = rfc4867::write_rtp(None, c, p).unwrap();
                rfc4867::read_rtp(&pkt, p).unwrap().1.collect::<Vec<_>>()
            })
            .collect()
    }
    for (label, path) in paths {
        run_all(&format!("rust via {label}"), &[], |case: &Case| {
            let serial: Vec<rfc4867::Frame> = case
                .frames
                .iter()
                .map(|f| rfc4867::Frame {
                    frame_type: rx_type(f.rx_type),
                    mode: mode(f.mode),
                    bits: f.bits,
                })
                .collect();
            let mut dec = amr_nb::Decoder::new();
            let mut out = Vec::with_capacity(case.frames.len() * L_FRAME);
            let mut pcm = [0i16; L_FRAME];
            for f in path(&serial) {
                dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
                out.extend_from_slice(&pcm);
            }
            out
        });
    }
}

/// Encoder homing frame (TS 26.073): 160 samples of 0x0008.
fn is_ehf(pcm: &[i16]) -> bool {
    pcm.iter().all(|&s| s == 0x0008)
}

/// C encoder oracle with the reference encoder's homing (opencore lacks it). Proves the encoder
/// harness (file parsing, mode schedule, `.COD` layout) before the Rust encoder uses it.
/// opencore builds VAD option 1 only, so the DTX2 (VAD2) cases are left out.
#[test]
fn c_encoder_oracle_with_homing_is_bit_exact() {
    run_all_enc(
        "c-encoder-oracle (no VAD2)",
        &common::OPENCORE_DEFECTS,
        |c| c.vad != 2,
        |case: &EncCase| c_encode(case.vad != 0, &case.modes, &case.input),
    );
}

/// The C encoder with homing over `input` (linear, whole frames), one requested mode per frame.
fn c_encode(dtx: bool, modes: &[i16], input: &[i16]) -> Vec<[i16; COD_FRAME]> {
    {
        {
            let (mut enc, mut sid) = (null_mut(), null_mut());
            unsafe { AMREncodeInit(&mut enc, &mut sid, dtx as i32) };
            let mut out = Vec::new();
            for (pcm, &mode) in input.as_chunks::<L_FRAME>().0.iter().zip(modes) {
                let mut speech = *pcm; // C overwrites its input
                let mut ets = [0i16; MAX_SERIAL_SIZE + 2];
                let mut ft = 0;
                unsafe {
                    AMREncode(
                        enc,
                        sid,
                        mode as i32,
                        speech.as_mut_ptr(),
                        ets.as_mut_ptr().cast(),
                        &mut ft,
                        AMR_TX_ETS,
                    );
                }
                let mut cod = [0i16; COD_FRAME];
                cod[..MAX_SERIAL_SIZE + 2].copy_from_slice(&ets);
                out.push(cod);
                if is_ehf(pcm) {
                    unsafe { AMREncodeReset(enc, sid) };
                }
            }
            unsafe { AMREncodeExit(&mut enc, &mut sid) };
            out
        }
    }
}

/// `.COD` words of one encoder output frame.
fn cod_frame(f: &amr_nb::rfc4867::Frame) -> [i16; COD_FRAME] {
    use amr_nb::RxFrameType::*;
    let mut cod = [0i16; COD_FRAME];
    cod[0] = match f.frame_type {
        SpeechGood => 0,
        SidFirst => 1,
        SidUpdate => 2,
        NoData => 3,
        t => panic!("unexpected encoder frame type {t:?}"),
    };
    cod[1..=MAX_SERIAL_SIZE].copy_from_slice(&f.bits);
    cod[1 + MAX_SERIAL_SIZE] = if f.frame_type == NoData {
        -1
    } else {
        f.mode as i16
    };
    cod
}

/// The Rust encoder (homing built in) against the official vectors.
#[test]
fn rust_encoder_is_bit_exact() {
    run_all_enc(
        "rust-encoder",
        &[],
        |_| true,
        |case: &EncCase| rust_encode(case.vad, &case.modes, &case.input),
    );
}

/// The Rust encoder over `input` (linear, whole frames), one requested mode per frame.
fn rust_encode(vad: u8, modes: &[i16], input: &[i16]) -> Vec<[i16; COD_FRAME]> {
    let mut enc = amr_nb::Encoder::new(
        [amr_nb::Dtx::Off, amr_nb::Dtx::Vad1, amr_nb::Dtx::Vad2][vad as usize],
    );
    input
        .as_chunks::<L_FRAME>()
        .0
        .iter()
        .zip(modes)
        .map(|(pcm, &m)| cod_frame(&enc.encode(mode(m), pcm)))
        .collect()
}

/// Bit and frame synchronisation (TS 26.074 §8), linear and companded: Rust and the C oracle.
#[test]
fn rust_encoder_sync_is_bit_exact() {
    run_all_sync("rust-encoder", &common::SYNC_ALL, |m, input| {
        rust_encode(0, &vec![m; input.len() / L_FRAME], input)
    });
}

/// The C oracle on the sync sequences: its frame-sync counts are pinned (opencore wraps on
/// them: RESULTS.md defect 10); bit sync as amr-nb's.
#[test]
fn c_encoder_oracle_sync() {
    run_all_sync(
        "c-encoder-oracle",
        &common::OPENCORE_SYNC_EXACT,
        |m, input| c_encode(false, &vec![m; input.len() / L_FRAME], input),
    );
}

#[path = "../examples/amr.rs"]
#[allow(dead_code)] // its `main`: the test calls `run`, as the command line does
mod amr;

/// `examples/amr.rs` as a user runs it: every linear constant-mode encoder input, plain and
/// DTX with VAD option 1 (the CLI's `dtx`), through `amr encode <in> <out.amr> <rate> [dtx]`
/// and `amr decode`, equals the decoder vector of that `.COD`.
#[test]
fn example_cli_round_trips_to_the_decoder_vectors() {
    const RATES: [&str; 8] = ["4.75", "5.15", "5.9", "6.7", "7.4", "7.95", "10.2", "12.2"];
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("amr-cli");
    std::fs::create_dir_all(&dir).unwrap();
    let outs: std::collections::HashMap<String, Vec<i16>> = common::cases()
        .into_iter()
        .filter(|c| c.law == common::Law::Linear)
        .map(|c| (c.name, c.expected))
        .collect();
    let cases: Vec<EncCase> = common::encoder_cases_where(|n| outs.contains_key(n))
        .into_iter()
        .filter(|c| c.vad < 2 && c.modes.iter().all(|&m| m == c.modes[0]))
        .collect();
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let failing: Vec<String> = std::thread::scope(|s| {
        let jobs: Vec<_> = (0..threads)
            .map(|k| {
                let (cases, dir, outs) = (&cases, &dir, &outs);
                s.spawn(move || {
                    let mut failing = Vec::new();
                    for (i, case) in cases.iter().enumerate().skip(k).step_by(threads) {
                        let path = |ext: &str| dir.join(format!("{i}.{ext}")).display().to_string();
                        let (inp, file, out) = (path("raw"), path("amr"), path("out"));
                        let bytes: Vec<u8> =
                            case.input.iter().flat_map(|v| v.to_le_bytes()).collect();
                        std::fs::write(&inp, bytes).unwrap();
                        let mut args = vec!["encode", &inp, &file, RATES[case.modes[0] as usize]];
                        if case.vad == 1 {
                            args.push("dtx");
                        }
                        amr::run(&args).unwrap();
                        amr::run(&["decode", &file, &out]).unwrap();
                        let got: Vec<i16> = std::fs::read(&out)
                            .unwrap()
                            .as_chunks::<2>()
                            .0
                            .iter()
                            .map(|&b| i16::from_le_bytes(b))
                            .collect();
                        if got != outs[&case.name] {
                            failing.push(case.name.clone());
                        }
                    }
                    failing
                })
            })
            .collect();
        jobs.into_iter().flat_map(|j| j.join().unwrap()).collect()
    });
    println!(
        "examples/amr.rs: {}/{} encoder inputs (plain and DTX) round-trip to their decoder vectors bit-exact",
        cases.len() - failing.len(),
        cases.len()
    );
    assert!(failing.is_empty(), "differ: {failing:?}");
}
