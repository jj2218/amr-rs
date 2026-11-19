//! Conformance against the 3GPP TS 26.174 encoder vectors (every `.inp` with its `.cod`):
//! proving the harness with the TS 26.173 reference and vo-amrwbenc.
mod common;

use amr_wb::L_FRAME;
use common::{EncCase, encoder_cases};

/// Encodes a case frame by frame with `enc(mode, dtx, frame, out) -> words`.
fn encode(
    case: &EncCase,
    mut enc: impl FnMut(u8, &[i16; 320], &mut [i16; 480]) -> usize,
) -> Vec<i16> {
    let mut out = Vec::with_capacity(case.expected.len());
    let mut buf = [0i16; 480];
    for (n, f) in case.input.as_chunks::<320>().0.iter().enumerate() {
        let w = enc(case.mode(n), f, &mut buf);
        out.extend_from_slice(&buf[..w]);
    }
    out
}

fn encode_ref(case: &EncCase) -> Vec<i16> {
    let mut e = ref26173_sys::Encoder::new(case.dtx);
    encode(case, |m, f, o| e.encode(m as i16, f, o))
}

#[allow(unsafe_code)]
fn encode_vo(case: &EncCase) -> Vec<i16> {
    // SAFETY: the shim allocates its state, copies the 320 samples, writes at most 480 words.
    unsafe {
        let st = opencore_sys::vo_wbenc_init();
        let out = encode(case, |m, f, o| {
            opencore_sys::vo_wbenc_encode(st, m as i32, case.dtx as i32, f.as_ptr(), o.as_mut_ptr())
                as usize
        });
        opencore_sys::vo_wbenc_exit(st);
        out
    }
}

/// The frame (counting the `.cod` frame headers) of the first differing word, if any.
fn first_diff(got: &[i16], want: &[i16]) -> Option<usize> {
    let i = got
        .iter()
        .zip(want)
        .position(|(a, b)| a != b)
        .or((got.len() != want.len()).then(|| got.len().min(want.len())))?;
    Some(
        frames(&want[..(i + 1).min(want.len())])
            .len()
            .saturating_sub(1),
    )
}

/// Splits `.cod` words into frames (each starts with 0x6b21; a ±127 bit word never equals it).
fn frames(w: &[i16]) -> Vec<&[i16]> {
    let starts: Vec<usize> = (0..w.len()).filter(|&i| w[i] == 0x6b21).collect();
    starts
        .iter()
        .zip(starts.iter().skip(1).chain([&w.len()]))
        .map(|(&a, &b)| &w[a..b])
        .collect()
}

/// Runs `encode` on every case; prints one summary line and the failing names. Frame-sync
/// cases (TS 26.174, as TS 26.074 §8.2) compare the 2nd frame after the last decoder homing
/// frame, the output for an all-0x0008 input frame.
fn run(label: &str, encode: impl Fn(&EncCase) -> Vec<i16> + Sync, threads: usize) -> Vec<String> {
    let cases = encoder_cases();
    let dhf: Vec<Vec<i16>> = (0..9)
        .map(|m| {
            encode(&EncCase {
                name: String::new(),
                input: vec![8; 320],
                modes: vec![m],
                dtx: false,
                sync: false,
                expected: Vec::new(),
            })
        })
        .collect();
    let check = |c: &EncCase| {
        let mut got = encode(c);
        if c.sync {
            let f = frames(&got);
            let at = f.iter().rposition(|&f| f == dhf[c.modes[0] as usize]);
            got = at
                .and_then(|i| f.get(i + 2))
                .map_or(Vec::new(), |f| f.to_vec());
        }
        first_diff(&got, &c.expected).map(|f| format!("{} (frame {f})", c.name))
    };
    // cases are independent: in `threads` interleaved shares (failures listed in case order)
    let mut failing: Vec<(usize, String)> = std::thread::scope(|sc| {
        let workers: Vec<_> = (0..threads)
            .map(|t| {
                let (cases, check) = (&cases, &check);
                sc.spawn(move || {
                    (t..cases.len())
                        .step_by(threads)
                        .filter_map(|i| check(&cases[i]).map(|f| (i, f)))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| w.join().unwrap())
            .collect()
    });
    failing.sort();
    let failing: Vec<String> = failing.into_iter().map(|(_, f)| f).collect();
    println!(
        "{label}: {}/{} TS 26.174 encoder vectors bit-exact",
        cases.len() - failing.len(),
        cases.len()
    );
    for f in &failing {
        println!("  {label} failing: {f}");
    }
    failing
}

/// Bit synchronisation (as TS 26.074 §8.1): `bitsync.inp` holds bit-shifted homing-frame
/// triplets; per mode, a fresh encoder on each triplet emits a decoder homing frame for
/// exactly the unshifted one (the first). Returns the failing modes.
fn bit_sync(label: &str, encode: impl Fn(&EncCase) -> Vec<i16>) -> Vec<String> {
    let bits = common::read_bitsync();
    let mut failing = Vec::new();
    for m in 0..9u8 {
        let case = |input: Vec<i16>| EncCase {
            name: String::new(),
            input,
            modes: vec![m],
            dtx: false,
            sync: false,
            expected: Vec::new(),
        };
        let dhf = encode(&case(vec![8; 320]));
        let hits: Vec<usize> = bits
            .chunks(3 * 320)
            .enumerate()
            .filter(|(_, t)| frames(&encode(&case(t.to_vec()))).contains(&&dhf[..]))
            .map(|(i, _)| i)
            .collect();
        if hits != [0] {
            failing.push(format!("mode {m}: homing frame after triplets {hits:?}"));
        }
    }
    println!(
        "{label} bit sync: {}/9 modes ({} triplets each)",
        9 - failing.len(),
        bits.len() / 960
    );
    failing
}

/// The harness itself: the reference must reproduce every vector (one thread: the reference
/// keeps a global overflow flag).
#[test]
fn reference_encoder_on_vectors() {
    let mut failing = run("ref26173 encoder (C reference)", encode_ref, 1);
    failing.extend(bit_sync("ref26173 encoder (C reference)", encode_ref));
    assert!(
        failing.is_empty(),
        "the reference must reproduce every vector: {failing:?}"
    );
}

/// vo-amrwbenc's count is reported (it predates CR 26.173-0032 and -0036), not asserted.
#[test]
fn vo_amrwbenc_on_vectors() {
    let _ = run("vo-amrwbenc (C baseline)", encode_vo, 1);
    let _ = bit_sync("vo-amrwbenc (C baseline)", encode_vo);
}

fn encode_rust(case: &EncCase) -> Vec<i16> {
    let mut e = amr_wb::Encoder::new(if case.dtx {
        amr_wb::Dtx::On
    } else {
        amr_wb::Dtx::Off
    });
    encode(case, |m, f, o| {
        let w = common::cod_words(&e.encode(amr_wb::Mode::from_index(m as usize).unwrap(), f));
        o[..w.len()].copy_from_slice(&w);
        w.len()
    })
}

/// The Rust encoder on every TS 26.174 vector, and bit synchronisation.
#[test]
fn rust_encoder_on_vectors() {
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get().min(8));
    let mut failing = run("amr-wb encoder (Rust)", encode_rust, threads);
    failing.extend(bit_sync("amr-wb encoder (Rust)", encode_rust));
    assert!(failing.is_empty(), "{} failing: {failing:?}", failing.len());
}

#[path = "../examples/awb.rs"]
#[allow(dead_code)] // its `main`: the test calls `run`, as the command line does
mod awb;

/// `examples/awb.rs` as a user runs it: every constant-mode encoder input, plain and DTX,
/// through `awb encode <in> <out.awb> <rate> [dtx]` and `awb decode`, equals the decoder
/// vector of that `.cod` (220 cases: T22 switches modes, which the CLI can't), plus the frame
/// the CLI pads a final partial one into (the reference drops it).
#[test]
fn example_cli_round_trips_to_the_decoder_vectors() {
    const RATES: [&str; 9] = [
        "6.60", "8.85", "12.65", "14.25", "15.85", "18.25", "19.85", "23.05", "23.85",
    ];
    let dir = std::path::Path::new(env!("CARGO_TARGET_TMPDIR")).join("awb-cli");
    std::fs::create_dir_all(&dir).unwrap();
    let t = common::vectors_dir().join("T");
    let cases: Vec<_> = encoder_cases()
        .into_iter()
        .filter(|c| !c.sync && c.modes.len() == 1 && t.join(&c.name).with_extension("out").exists())
        .collect();
    let threads = std::thread::available_parallelism().map_or(4, |n| n.get());
    let failing: Vec<String> = std::thread::scope(|s| {
        let jobs: Vec<_> = (0..threads)
            .map(|k| {
                let (cases, dir, t) = (&cases, &dir, &t);
                s.spawn(move || {
                    let mut failing = Vec::new();
                    for (i, case) in cases.iter().enumerate().skip(k).step_by(threads) {
                        let path = |ext: &str| dir.join(format!("{i}.{ext}")).display().to_string();
                        let (inp, file, out) = (path("raw"), path("awb"), path("out"));
                        let bytes: Vec<u8> =
                            case.input.iter().flat_map(|v| v.to_le_bytes()).collect();
                        std::fs::write(&inp, bytes).unwrap();
                        let mut args = vec!["encode", &inp, &file, RATES[case.modes[0] as usize]];
                        if case.dtx {
                            args.push("dtx");
                        }
                        awb::run(&args).unwrap();
                        awb::run(&["decode", &file, &out]).unwrap();
                        let got: Vec<i16> = std::fs::read(&out)
                            .unwrap()
                            .as_chunks::<2>()
                            .0
                            .iter()
                            .map(|&b| i16::from_le_bytes(b))
                            .collect();
                        let want = common::read_i16(&t.join(&case.name).with_extension("out"));
                        // the reference drops a final partial frame, the CLI pads it
                        let padded = got.len() == case.input.len().div_ceil(L_FRAME) * L_FRAME;
                        if !padded || got[..want.len()] != want[..] {
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
        "examples/awb.rs: {}/{} encoder inputs (plain and DTX) round-trip to their decoder vectors bit-exact",
        cases.len() - failing.len(),
        cases.len()
    );
    assert!(failing.is_empty(), "differ: {failing:?}");
}
