//! Conformance against the 3GPP TS 26.174 decoder vectors (every `.cod` with its `.out`):
//! the Rust decoder, and (proving the harness) the TS 26.173 reference and opencore.
mod common;

use amr_wb::{Decoder, L_FRAME16K, RxFrameType};
use common::{Case, cases};

fn decode_rust(case: &Case) -> Vec<i16> {
    let mut dec = Decoder::new();
    let mut out = Vec::with_capacity(case.frames.len() * L_FRAME16K);
    let mut pcm = [0i16; L_FRAME16K];
    for f in &case.frames {
        dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
        out.extend_from_slice(&pcm);
    }
    out
}

#[allow(unsafe_code)]
fn decode_opencore(case: &Case) -> Vec<i16> {
    // SAFETY: the shim allocates its state; prms holds every bit it reads (477 words).
    unsafe {
        let st = opencore_sys::oc_wb_init();
        let mut out = Vec::new();
        let mut pcm = [0i16; L_FRAME16K];
        for f in &case.frames {
            let mut prms = f.prms;
            opencore_sys::oc_wb_decode(
                st,
                f.mode as i16,
                f.frame_type as i16,
                prms.as_mut_ptr(),
                pcm.as_mut_ptr(),
            );
            out.extend_from_slice(&pcm);
        }
        opencore_sys::oc_wb_exit(st);
        out
    }
}

fn decode_ref(case: &Case) -> Vec<i16> {
    let mut dec = ref26173_sys::Decoder::new();
    let mut out = Vec::new();
    let mut pcm = [0i16; L_FRAME16K];
    for f in &case.frames {
        let mut prms = f.prms;
        dec.decode(f.mode as i16, f.frame_type as i16, &mut prms, &mut pcm);
        out.extend_from_slice(&pcm);
    }
    out
}

/// First differing sample as (frame, sample), if any.
fn first_diff(got: &[i16], want: &[i16]) -> Option<(usize, usize)> {
    if got.len() != want.len() {
        return Some((got.len().min(want.len()) / L_FRAME16K, 0));
    }
    got.iter()
        .zip(want)
        .position(|(a, b)| a != b)
        .map(|i| (i / L_FRAME16K, i % L_FRAME16K))
}

/// Runs `decode` on every case; prints one summary line and the failing names.
fn run(label: &str, decode: impl Fn(&Case) -> Vec<i16>) -> Vec<String> {
    let cases = cases();
    let mut failing = Vec::new();
    for c in &cases {
        if let Some((f, s)) = first_diff(&decode(c), &c.expected) {
            failing.push(format!("{} (frame {f}, sample {s})", c.name));
        }
    }
    println!(
        "{label}: {}/{} TS 26.174 decoder vectors bit-exact",
        cases.len() - failing.len(),
        cases.len()
    );
    for f in &failing {
        println!("  {label} failing: {f}");
    }
    failing
}

#[test]
fn decoder_vectors() {
    let failing = run("amr-wb decoder", decode_rust);
    assert!(failing.is_empty(), "{} vectors differ", failing.len());
}

/// The harness itself: the reference must reproduce every vector. opencore's count is
/// reported (it has known deviations), not asserted.
#[test]
fn c_decoders_on_vectors() {
    let failing = run("ref26173 (C reference)", decode_ref);
    assert!(
        failing.is_empty(),
        "the reference must reproduce every vector"
    );
    let _ = run("opencore (C baseline)", decode_opencore);
}

/// No vector sends NO_DATA right after a speech frame, the case CR 26.173-0027 changed (the
/// one deliberate difference from the reference); were one added, the oracle would diverge.
#[test]
fn no_vector_has_no_data_mid_speech() {
    for c in cases() {
        let mut prev_speech = false;
        for (n, f) in c.frames.iter().enumerate() {
            let speech = matches!(
                f.frame_type,
                RxFrameType::SpeechGood
                    | RxFrameType::SpeechProbablyDegraded
                    | RxFrameType::SpeechBad
                    | RxFrameType::SpeechLost
            );
            if f.frame_type == RxFrameType::NoData {
                assert!(
                    !prev_speech,
                    "{}: NO_DATA after speech at frame {n}",
                    c.name
                );
            } else {
                prev_speech = speech;
            }
        }
    }
}
