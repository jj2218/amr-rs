//! Codec speed: Rust vs opencore (C, -O3) on the 3GPP NODTX vectors, per mode: the decoder on
//! the `.COD` streams, the encoder on the `.INP` speech. `BENCH=dec` or `BENCH=enc` runs one
//! side only. Run via `just bench` (pins one core). Prints µs/frame (best of BENCH_REPS,
//! default 15) and the speedup. `just profile` runs it once under callgrind.
#![allow(unsafe_code)] // FFI to the C baseline

#[path = "../tests/common/mod.rs"]
mod common;

use std::ffi::c_void;
use std::hint::black_box;
use std::ptr::null_mut;
use std::time::Instant;

unsafe extern "C" {
    fn GSMInitDecode(state: *mut *mut c_void, id: *mut i8) -> i16;
    fn GSMDecodeFrameExit(state: *mut *mut c_void);
    fn AMRDecode(
        state: *mut c_void,
        frame_type: i32,
        bits: *mut u8,
        pcm: *mut i16,
        format: i32,
    ) -> i16;
    fn AMREncodeInit(enc: *mut *mut c_void, sid: *mut *mut c_void, dtx: i32) -> i16;
    fn AMREncodeExit(enc: *mut *mut c_void, sid: *mut *mut c_void);
    fn AMREncode(
        enc: *mut c_void,
        sid: *mut c_void,
        mode: i32,
        input: *mut i16,
        out: *mut u8,
        frame_type: *mut i32,
        format: i16,
    ) -> i16;
}
use opencore_sys as _;

const MODES: [&str; 8] = ["475", "515", "59", "67", "74", "795", "102", "122"];

fn reps() -> usize {
    std::env::var("BENCH_REPS")
        .ok()
        .and_then(|r| r.parse().ok())
        .unwrap_or(15)
}

/// `BENCH_MODE=<mode>` (e.g. 102) runs only that mode (`just bench`, `just profile`).
fn wanted(name: &str) -> bool {
    std::env::var("BENCH_MODE").map_or(true, |m| m == name)
}

/// Exits on an unknown `BENCH_MODE` (it would run nothing) or `BENCH` (it would run both sides).
fn check_mode() {
    if let Ok(m) = std::env::var("BENCH_MODE")
        && !MODES.contains(&m.as_str())
    {
        eprintln!(
            "bench: unknown BENCH_MODE '{m}' (modes: {})",
            MODES.join(" ")
        );
        std::process::exit(2);
    }
    if let Ok(b) = std::env::var("BENCH")
        && !["", "dec", "enc"].contains(&b.as_str())
    {
        eprintln!("bench: unknown BENCH '{b}' (dec, enc; unset: both)");
        std::process::exit(2);
    }
}

fn best<F: FnMut()>(mut f: F) -> f64 {
    (0..reps())
        .map(|_| {
            let t = Instant::now();
            f();
            t.elapsed().as_secs_f64()
        })
        .fold(f64::INFINITY, f64::min)
}

/// `AMREncode` output format: ETS (one Word16 per bit), the cheapest packing.
const AMR_TX_ETS: i16 = 2;

const MODE: [amr_nb::Mode; 8] = [
    amr_nb::Mode::MR475,
    amr_nb::Mode::MR515,
    amr_nb::Mode::MR59,
    amr_nb::Mode::MR67,
    amr_nb::Mode::MR74,
    amr_nb::Mode::MR795,
    amr_nb::Mode::MR102,
    amr_nb::Mode::MR122,
];

fn main() {
    check_mode();
    let which = std::env::var("BENCH").unwrap_or_default();
    if which != "enc" {
        decoder();
    }
    if which != "dec" {
        encoder();
    }
}

/// Prints one table row.
fn row(name: &str, rust: f64, c: f64, frames: usize) {
    let (ur, uc) = (rust / frames as f64 * 1e6, c / frames as f64 * 1e6);
    println!("{name:<6} {ur:>10.3} {uc:>10.3} {:>7.2}x", uc / ur);
}

/// The "all" row (every mode timed: not with `BENCH_MODE`, nor in a stages build).
fn total(tr: f64, tc: f64, tn: usize) {
    if tn > 0 && std::env::var("BENCH_MODE").is_err() {
        row("all", tr, tc, tn);
        println!("({tn} frames/rep, best of {})", reps());
    }
}

/// `--cfg codec_stages` (`BENCH=stages`): µs per frame per stage (mean of all reps, timer
/// overhead included), then the frame's total.
#[cfg(codec_stages)]
fn stages_row(name: &str, n: usize, frames: usize) {
    let ns = amr_nb::stages::take();
    let per = (frames * reps()) as f64 * 1e3;
    let cols: String = ns[..n]
        .iter()
        .map(|&v| format!(" {:>13.2}", v as f64 / per))
        .collect();
    let total = ns.iter().sum::<u64>() as f64 / per;
    println!("stages {name:<5}{cols}  = {total:.2}");
}

/// The stage table's header.
#[cfg(codec_stages)]
fn stages_head(side: &str, names: &[&str]) {
    let cols: String = names.iter().map(|n| format!(" {n:>13}")).collect();
    println!("stages {side}  {cols}");
}

fn decoder() {
    let cases = common::cases();
    #[cfg(codec_stages)]
    stages_head("dec", &amr_nb::stages::NAMES_DEC);
    println!("decoder {:>9} {:>10} {:>8}", "rust µs", "c µs", "speedup");
    let (mut tr, mut tc, mut tn) = (0.0, 0.0, 0usize);
    for (m, name) in MODES.iter().enumerate().filter(|(_, n)| wanted(n)) {
        let set: Vec<_> = cases
            .iter()
            .filter(|c| c.name.starts_with("nodtx/") && c.name.ends_with(&format!("_{name}")))
            .collect();
        let frames: usize = set.iter().map(|c| c.frames.len()).sum();
        let mode = MODE[m];
        #[cfg(codec_stages)]
        amr_nb::stages::take();
        let rust = best(|| {
            for c in &set {
                let mut d = amr_nb::Decoder::new();
                let mut pcm = [0i16; common::L_FRAME];
                for f in &c.frames {
                    d.decode(amr_nb::RxFrameType::SpeechGood, mode, &f.bits, &mut pcm);
                    black_box(&pcm);
                }
            }
        });
        #[cfg(codec_stages)]
        {
            stages_row(name, amr_nb::stages::NAMES_DEC.len(), frames);
            continue;
        }
        // Pre-built ETS frames so the C timing excludes Rust-side formatting.
        let ets: Vec<Vec<[i16; common::MAX_SERIAL_SIZE + 2]>> = set
            .iter()
            .map(|c| {
                c.frames
                    .iter()
                    .map(|f| {
                        let mut e = [0i16; common::MAX_SERIAL_SIZE + 2];
                        e[1..=common::MAX_SERIAL_SIZE].copy_from_slice(&f.bits);
                        e[common::MAX_SERIAL_SIZE + 1] = m as i16;
                        e
                    })
                    .collect()
            })
            .collect();
        let mut ets = ets;
        let c = best(|| {
            for frames in ets.iter_mut() {
                let mut st = null_mut();
                let mut pcm = [0i16; common::L_FRAME];
                unsafe { GSMInitDecode(&mut st, null_mut()) };
                for e in frames.iter_mut() {
                    unsafe { AMRDecode(st, 0, e.as_mut_ptr().cast(), pcm.as_mut_ptr(), 0) };
                    black_box(&pcm);
                }
                unsafe { GSMDecodeFrameExit(&mut st) };
            }
        });
        row(name, rust, c, frames);
        (tr, tc, tn) = (tr + rust, tc + c, tn + frames);
    }
    total(tr, tc, tn);
}

/// Encoder (DTX off) on the NODTX speech inputs, one case per input file and mode.
fn encoder() {
    // NODTX speech `T<nn>_<mode>` only (not T21's mode switching): reads only those files.
    let cases = common::encoder_cases_where(|n| n.starts_with("nodtx/T") && n.contains('_'));
    #[cfg(codec_stages)]
    stages_head("enc", &amr_nb::stages::NAMES_ENC);
    println!("encoder {:>9} {:>10} {:>8}", "rust µs", "c µs", "speedup");
    let (mut tr, mut tc, mut tn) = (0.0, 0.0, 0usize);
    for (m, name) in MODES.iter().enumerate().filter(|(_, n)| wanted(n)) {
        let inputs: Vec<&[[i16; common::L_FRAME]]> = cases
            .iter()
            .filter(|c| c.name.starts_with("nodtx/T") && c.name.ends_with(&format!("_{name}")))
            .map(|c| c.input.as_chunks::<{ common::L_FRAME }>().0)
            .collect();
        let frames: usize = inputs.iter().map(|i| i.len()).sum();
        #[cfg(codec_stages)]
        amr_nb::stages::take();
        let rust = best(|| {
            for input in &inputs {
                let mut e = amr_nb::Encoder::new(amr_nb::Dtx::Off);
                for pcm in input.iter() {
                    black_box(e.encode(MODE[m], pcm));
                }
            }
        });
        #[cfg(codec_stages)]
        {
            stages_row(name, amr_nb::stages::NAMES_ENC.len(), frames);
            continue;
        }
        let c = best(|| {
            for input in &inputs {
                let (mut enc, mut sid) = (null_mut(), null_mut());
                let mut out = [0i16; common::MAX_SERIAL_SIZE + 2];
                let mut ft = 0;
                unsafe { AMREncodeInit(&mut enc, &mut sid, 0) };
                for pcm in input.iter() {
                    let mut speech = *pcm; // C overwrites its input
                    unsafe {
                        AMREncode(
                            enc,
                            sid,
                            m as i32,
                            speech.as_mut_ptr(),
                            out.as_mut_ptr().cast(),
                            &mut ft,
                            AMR_TX_ETS,
                        )
                    };
                    black_box(&out);
                }
                unsafe { AMREncodeExit(&mut enc, &mut sid) };
            }
        });
        row(name, rust, c, frames);
        (tr, tc, tn) = (tr + rust, tc + c, tn + frames);
    }
    total(tr, tc, tn);
}
