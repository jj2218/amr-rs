//! Conformance against the 3GPP TS 26.174 decoder vectors (every `.cod` with its `.out`):
//! the Rust decoder, and (proving the harness) the TS 26.173 reference and opencore.
mod common;

use amr_wb::{Decoder, L_FRAME, RxFrameType, rfc4867};
use common::{Case, cases};

fn decode_rust(case: &Case) -> Vec<i16> {
    let mut dec = Decoder::new();
    let mut out = Vec::with_capacity(case.frames.len() * L_FRAME);
    let mut pcm = [0i16; L_FRAME];
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
        let mut pcm = [0i16; L_FRAME];
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
    let mut pcm = [0i16; L_FRAME];
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
        return Some((got.len().min(want.len()) / L_FRAME, 0));
    }
    got.iter()
        .zip(want)
        .position(|(a, b)| a != b)
        .map(|i| (i / L_FRAME, i % L_FRAME))
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

/// Every vector through an RFC 4867 `.awb` file: written, read back, decoded. Speech
/// frames marked "probably degraded" have no packed form (none in the set: asserted).
#[test]
fn awb_files_decode_like_the_vectors() {
    let decode_awb = |case: &Case| {
        let frames: Vec<rfc4867::Frame> = case
            .frames
            .iter()
            .map(|f| rfc4867::Frame {
                frame_type: f.frame_type,
                mode: f.mode,
                bits: f.bits,
            })
            .collect();
        let file = rfc4867::write_storage(&frames).expect("every vector frame is packable");
        let mut dec = Decoder::new();
        let mut out = Vec::with_capacity(frames.len() * L_FRAME);
        let mut pcm = [0i16; L_FRAME];
        let mut n = 0;
        for f in rfc4867::read_storage(&file).unwrap() {
            let f = f.unwrap();
            dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
            out.extend_from_slice(&pcm);
            n += 1;
        }
        assert_eq!(n, frames.len(), "{}", case.name);
        out
    };
    let failing = run("amr-wb via .awb (RFC 4867)", decode_awb);
    assert!(failing.is_empty(), "{} vectors differ", failing.len());
}

/// A SID frame's mode is ignored (the previous mode is used, as the reference reads `.awb`
/// files): the DTX vectors decode unchanged with every SID frame's mode rotated.
#[test]
fn sid_mode_is_ignored() {
    let mut sids = 0;
    for mut case in cases().into_iter().filter(|c| c.name.contains("Dtx")) {
        for f in &mut case.frames {
            if matches!(
                f.frame_type,
                RxFrameType::SidFirst | RxFrameType::SidUpdate | RxFrameType::SidBad
            ) {
                f.mode = amr_wb::Mode::from_index((f.mode as usize + 4) % 9).unwrap();
                sids += 1;
            }
        }
        assert_eq!(
            first_diff(&decode_rust(&case), &case.expected),
            None,
            "{}",
            case.name
        );
    }
    assert!(sids > 100, "only {sids} SID frames");
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
