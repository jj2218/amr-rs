//! Conformance against the 3GPP TS 26.174 encoder vectors (every `.inp` with its `.cod`):
//! proving the harness with the TS 26.173 reference and vo-amrwbenc.
mod common;

use common::{EncCase, enc_cases};

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
    let cases = enc_cases();
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
    let mut e = amr_wb::Encoder::new(case.dtx);
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
