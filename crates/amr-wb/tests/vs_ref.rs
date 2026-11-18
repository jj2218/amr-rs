//! Whole-decoder differential tests against the TS 26.173 reference C on stress streams the
//! vectors don't have: corrupted and random bits (loud, clipping signals), bad, lost and
//! degraded frames, SID and NO_DATA sequences, mode switching on every frame
//! (`common::stress`). One summary line per stream class.
//!
//! The one deliberate difference (CR 26.173-0027..29, see the crate docs) is kept out by
//! construction; `no_data_mid_speech_is_the_one_difference` pins it.
//! `WB_STRESS_STREAMS=<n>` runs more streams per class (default 40).
mod common;

use amr_wb::{Decoder, L_FRAME, RxFrameType};
use common::stress::{Mix, Rng, pool, stream};
use common::{Frame, cases};

/// Decodes `frames` with both decoders; the first differing frame, if any.
pub fn first_divergence(frames: &[Frame]) -> Option<usize> {
    let mut rust = Decoder::new();
    let mut c = ref26173_sys::Decoder::new();
    let (mut a, mut b) = ([0i16; L_FRAME], [0i16; L_FRAME]);
    for (n, f) in frames.iter().enumerate() {
        rust.decode(f.frame_type, f.mode, &f.bits, &mut a);
        let mut prms = f.prms;
        c.decode(f.mode as i16, f.frame_type as i16, &mut prms, &mut b);
        if a != b {
            return Some(n);
        }
    }
    None
}

#[test]
fn stress_streams_match_reference() {
    let n: u64 = std::env::var("WB_STRESS_STREAMS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(40);
    let pool = pool();
    let mut failing = 0;
    for k in 0..8u64 {
        let mix = Mix {
            bad_frames: k & 1 != 0,
            dtx: k & 2 != 0,
            random_bits: k & 4 != 0,
        };
        let mut bad = Vec::new();
        for seed in 1..=n {
            let frames = stream(&mut Rng::seeded(seed * 8 + k), &pool, 120, mix);
            if let Some(at) = first_divergence(&frames) {
                let ctx: Vec<String> = frames[at.saturating_sub(2)..=at]
                    .iter()
                    .map(|f| format!("{:?}/{:?}", f.frame_type, f.mode))
                    .collect();
                bad.push(format!(
                    "seed {}: frame {at} ({})",
                    seed * 8 + k,
                    ctx.join(", ")
                ));
            }
        }
        println!(
            "amr-wb vs ref26173, {mix:?}: {}/{n} streams bit-exact",
            n as usize - bad.len()
        );
        for b in bad.iter().take(3) {
            println!("  differs: {b}");
        }
        failing += bad.len();
    }
    assert_eq!(failing, 0, "streams differ from the reference");
}

/// The documented difference, pinned. NO_DATA frames right after speech change only the DTX
/// hangover bookkeeping, which shows at the next SID_FIRST: the reference (CR 26.173-0029)
/// counted them as speech and adds a hangover there (comfort noise from the averaged
/// history); this crate follows the pre-CR rule, which counted them as DTX. The control
/// stream (no NO_DATA before the SID) must match. If the pinned one stops differing, the
/// CR's rule has crept in.
#[test]
fn no_data_mid_speech_is_the_one_difference() {
    let cases = cases();
    let case = cases
        .iter()
        .find(|c| c.name == "T_1265/T00_1265.cod")
        .unwrap();
    let mut frames: Vec<Frame> = case.frames[..60].to_vec();
    frames[36].frame_type = RxFrameType::SidFirst;
    for f in &mut frames[37..] {
        f.frame_type = RxFrameType::NoData;
    }
    assert_eq!(
        first_divergence(&frames),
        None,
        "control: speech, SID_FIRST, NO_DATA"
    );
    for f in &mut frames[30..36] {
        f.frame_type = RxFrameType::NoData;
    }
    let d = first_divergence(&frames);
    assert!(
        matches!(d, Some(n) if n >= 36),
        "NO_DATA mid-speech must differ from the CR 0029 reference at or after the SID (got {d:?})"
    );
}

/// Sets `n` bits of `f` at `at` to `v` (MSB first), as the decoder reads parameters.
fn set_parm(f: &mut Frame, at: usize, n: usize, v: u32) {
    for k in 0..n {
        let b = (v >> (n - 1 - k) & 1) as i16;
        f.bits[at + k] = b;
        f.prms[at + k] = if b != 0 { 127 } else { -127 };
    }
}

/// Isf_isp reads `table[ind + 1]`; DTX can make ISF 14 exactly 16384 (CN_dithering's clamp),
/// `ind` = 128 on a 129-entry table. Speech frames with the largest ISF 14 the 46-bit
/// quantiser gives (16216) fill the DTX history; the SID_FIRST after them averages it, and the
/// dither (on from the earlier SID_UPDATE) pushes it past 16384 within a few frames.
#[test]
fn dtx_isf_at_the_table_end() {
    let cases = cases();
    let case = cases
        .iter()
        .find(|c| c.name == "T_1265/T00_1265.cod")
        .unwrap();
    let speech = &case.frames[..60];
    let mut frames: Vec<Frame> = speech[..10].to_vec();
    let mut sid = speech[10].clone();
    sid.frame_type = RxFrameType::SidFirst;
    frames.push(sid.clone());
    sid.frame_type = RxFrameType::SidUpdate;
    set_parm(&mut sid, 34, 1, 1); // CN dithering on
    frames.push(sid.clone());
    let mut nd = sid.clone();
    nd.frame_type = RxFrameType::NoData;
    frames.extend(std::iter::repeat_n(nd.clone(), 5));
    for f in &speech[12..52] {
        let mut f = f.clone();
        set_parm(&mut f, 1 + 8, 8, 208); // dico2_isf row with the largest ISF 14
        set_parm(&mut f, 1 + 8 + 8 + 6 + 7 + 7 + 5, 5, 23); // dico25_isf, likewise
        frames.push(f);
    }
    sid.frame_type = RxFrameType::SidFirst;
    frames.push(sid);
    frames.extend(std::iter::repeat_n(nd, 40));
    assert_eq!(first_divergence(&frames), None);
}
