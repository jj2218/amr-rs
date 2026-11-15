//! Decoder speed: Rust vs opencore (C, -O3) on the 3GPP TS 26.174 single-mode vectors
//! (`T_<mode>/T00..T20`), per mode. Both run the same frame loop (homing test, decode, 14-bit
//! output) on pre-parsed serial bits. Run via `CRATE=amr-wb just bench` (pins one core).
//! Prints µs/frame (best of BENCH_REPS, default 15) and the speedup; `BENCH_MODE=<mode>`
//! (e.g. 1265) runs one mode. `BENCH=enc` has nothing to run (no encoder yet).
#![allow(unsafe_code)] // FFI to the C baseline

#[path = "../tests/common/mod.rs"]
mod common;

use amr_wb::{Decoder, L_FRAME16K};
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
        println!("amr-wb: no encoder to benchmark");
        return;
    }
    let cases = common::cases();
    println!("decoder {:>9} {:>10} {:>8}", "rust µs", "c µs", "speedup");
    let (mut tr, mut tc, mut tn) = (0.0, 0.0, 0usize);
    for name in MODES.iter().filter(|n| wanted(n)) {
        let dir = format!("T_{name}/");
        let mut set: Vec<_> = cases.iter().filter(|c| c.name.starts_with(&dir)).collect();
        set.retain(|c| !c.name.contains("Dtx"));
        let frames: usize = set.iter().map(|c| c.frames.len()).sum();
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
}
