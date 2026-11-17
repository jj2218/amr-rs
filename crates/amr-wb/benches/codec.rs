//! Codec speed on the 3GPP TS 26.174 single-mode vectors (`T_<mode>/T00..T20`), per mode.
//! Decoder: Rust vs opencore (C, -O3), both the same frame loop (homing test, decode, 14-bit
//! output) on pre-parsed serial bits. `BENCH=enc`: encoder, Rust vs vo-amrwbenc (C, -O3) on
//! the vectors' input PCM. Run via `CRATE=amr-wb just bench` (pins one core). Prints µs/frame
//! (best of BENCH_REPS, default 15) and the speedup; `BENCH_MODE=<mode>` (e.g. 1265) runs one
//! mode.
#![allow(unsafe_code)] // FFI to the C baseline

#[path = "../tests/common/mod.rs"]
mod common;

use amr_wb::{Decoder, Encoder, L_FRAME16K, Mode};
use std::hint::black_box;
use std::time::Instant;

const MODES: [&str; 9] = [
    "660", "885", "1265", "1425", "1585", "1825", "1985", "2305", "2385",
];

fn reps() -> usize {
    std::env::var("BENCH_REPS")
        .ok()
        .and_then(|r| r.parse().ok())
        .unwrap_or(15)
}

/// `BENCH_MODE=<mode>` runs only that mode.
fn wanted(name: &str) -> bool {
    std::env::var("BENCH_MODE").map_or(true, |m| m == name)
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

fn row(name: &str, rust: f64, c: f64, frames: usize) {
    let (ur, uc) = (rust / frames as f64 * 1e6, c / frames as f64 * 1e6);
    println!("{name:<6} {ur:>10.3} {uc:>10.3} {:>7.2}x", uc / ur);
}

fn main() {
    if let Ok(m) = std::env::var("BENCH_MODE")
        && !MODES.contains(&m.as_str())
    {
        eprintln!(
            "bench: unknown BENCH_MODE '{m}' (modes: {})",
            MODES.join(" ")
        );
        std::process::exit(2);
    }
    if std::env::var("BENCH").as_deref() == Ok("enc") {
        encoder();
        return;
    }
    let cases = common::cases();
    #[cfg(codec_stages)]
    println!(
        "stages dec  {}",
        amr_wb::stages::NAMES
            .iter()
            .map(|n| format!(" {n:>11}"))
            .collect::<String>()
    );
    println!("decoder {:>9} {:>10} {:>8}", "rust µs", "c µs", "speedup");
    let (mut tr, mut tc, mut tn) = (0.0, 0.0, 0usize);
    for name in MODES.iter().filter(|n| wanted(n)) {
        let dir = format!("T_{name}/");
        let mut set: Vec<_> = cases.iter().filter(|c| c.name.starts_with(&dir)).collect();
        set.retain(|c| !c.name.contains("Dtx"));
        let frames: usize = set.iter().map(|c| c.frames.len()).sum();
        #[cfg(codec_stages)]
        amr_wb::stages::take();
        let rust = best(|| {
            for c in &set {
                let mut d = Decoder::new();
                let mut pcm = [0i16; L_FRAME16K];
                for f in &c.frames {
                    d.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
                    black_box(&pcm);
                }
            }
        });
        // `--cfg codec_stages` (`BENCH=stages`): µs per frame per stage (mean of all reps)
        #[cfg(codec_stages)]
        {
            let ns = amr_wb::stages::take();
            let per = (frames * reps()) as f64 * 1e3;
            let cols: String = ns[..amr_wb::stages::NAMES.len()]
                .iter()
                .map(|&v| format!(" {:>11.2}", v as f64 / per))
                .collect();
            println!("stages {name:<4}{cols}");
            continue;
        }
        // the C reads prms without writing them: one copy outside the timing
        let mut prms: Vec<Vec<_>> = set
            .iter()
            .map(|c| c.frames.iter().map(|f| f.prms).collect())
            .collect();
        let c = best(|| {
            for (case, p) in set.iter().zip(prms.iter_mut()) {
                // SAFETY: the shim allocates its state; each prms holds 477 words.
                unsafe {
                    let st = opencore_sys::oc_wb_init();
                    let mut pcm = [0i16; L_FRAME16K];
                    for (f, bits) in case.frames.iter().zip(p.iter_mut()) {
                        opencore_sys::oc_wb_decode(
                            st,
                            f.mode as i16,
                            f.frame_type as i16,
                            bits.as_mut_ptr(),
                            pcm.as_mut_ptr(),
                        );
                        black_box(&pcm);
                    }
                    opencore_sys::oc_wb_exit(st);
                }
            }
        });
        row(name, rust, c, frames);
        (tr, tc, tn) = (tr + rust, tc + c, tn + frames);
    }
    if std::env::var("BENCH_MODE").is_err() {
        row("all", tr, tc, tn);
    }
    // `BENCH=stages`: the encoder's stages too
    #[cfg(codec_stages)]
    encoder();
}

/// `BENCH=enc`: the encoders on the single-mode vectors' input PCM (T00..T20, no DTX), Rust vs
/// vo-amrwbenc, both through their frame loops (homing test, 14-bit input, serial output).
fn encoder() {
    let cases = common::enc_cases();
    #[cfg(codec_stages)]
    println!(
        "stages enc  {}",
        amr_wb::stages::NAMES_ENC
            .iter()
            .map(|n| format!(" {n:>9}"))
            .collect::<String>()
    );
    println!("encoder {:>9} {:>10} {:>8}", "rust µs", "c µs", "speedup");
    let (mut tr, mut tc, mut tn) = (0.0, 0.0, 0usize);
    for (m, name) in MODES.iter().enumerate().filter(|(_, n)| wanted(n)) {
        let dir = format!("T_{name}/");
        let set: Vec<_> = cases
            .iter()
            .filter(|c| c.name.starts_with(&dir) && !c.dtx && !c.sync)
            .collect();
        let frames: usize = set.iter().map(|c| c.input.len() / L_FRAME16K).sum();
        let mode = Mode::from_index(m).unwrap();
        #[cfg(codec_stages)]
        amr_wb::stages::take();
        let rust = best(|| {
            for c in &set {
                let mut e = Encoder::new(false);
                for f in c.input.as_chunks::<L_FRAME16K>().0 {
                    black_box(e.encode(mode, f));
                }
            }
        });
        // `--cfg codec_stages`: µs per frame per encoder stage (mean of all reps), then the
        // per-frame total
        #[cfg(codec_stages)]
        {
            let ns = amr_wb::stages::take();
            let per = (frames * reps()) as f64 * 1e3;
            let cols: String = ns[..amr_wb::stages::NAMES_ENC.len()]
                .iter()
                .map(|&v| format!(" {:>9.2}", v as f64 / per))
                .collect();
            let total = ns.iter().sum::<u64>() as f64 / per;
            println!("stages {name:<4}{cols}  = {total:.2}");
            continue;
        }
        let c = best(|| {
            for c in &set {
                // SAFETY: the shim allocates its state, copies the 320 samples, writes at
                // most 480 words.
                unsafe {
                    let st = opencore_sys::vo_wbenc_init();
                    let mut out = [0i16; 480];
                    for f in c.input.as_chunks::<L_FRAME16K>().0 {
                        black_box(opencore_sys::vo_wbenc_encode(
                            st,
                            m as i32,
                            0,
                            f.as_ptr(),
                            out.as_mut_ptr(),
                        ));
                        black_box(&out);
                    }
                    opencore_sys::vo_wbenc_exit(st);
                }
            }
        });
        row(name, rust, c, frames);
        (tr, tc, tn) = (tr + rust, tc + c, tn + frames);
    }
    if std::env::var("BENCH_MODE").is_err() {
        row("all", tr, tc, tn);
    }
}
