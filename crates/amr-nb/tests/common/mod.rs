//! Shared conformance-test plumbing: locating and parsing the 3GPP TS 26.074 test vectors.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

pub mod g711;
pub use g711::Law;

/// Samples per 20 ms frame.
pub const L_FRAME: usize = 160;
/// Serial bits per frame in the ETS format.
pub const MAX_SERIAL_SIZE: usize = 244;
/// Words per frame in a 26.074 `.COD` file: `[tx_type, 244 bits, mode, 0, 0, 0, 0]`.
const COD_WORDS: usize = 250;

/// RX frame types (3GPP TS 26.073 `enum RXFrameType`).
pub mod rx {
    pub const SPEECH_GOOD: i16 = 0;
    pub const SPEECH_DEGRADED: i16 = 1;
    pub const ONSET: i16 = 2;
    pub const SPEECH_BAD: i16 = 3;
    pub const SID_FIRST: i16 = 4;
    pub const SID_UPDATE: i16 = 5;
    pub const SID_BAD: i16 = 6;
    pub const NO_DATA: i16 = 7;
}

/// One decoder input frame.
#[derive(Clone)]
pub struct Frame {
    pub rx_type: i16,
    pub bits: [i16; MAX_SERIAL_SIZE],
    /// Codec mode (0..=7); -1 in the file for NO_DATA frames.
    pub mode: i16,
}

/// One `.COD` → `.OUT` conformance case.
pub struct Case {
    pub name: String,
    /// Sample encoding of `expected`: the decoder's linear output is compressed before comparing.
    pub law: Law,
    pub frames: Vec<Frame>,
    pub expected: Vec<i16>,
}

/// `.cache/3gpp` of this checkout, or of the main checkout when run from a git worktree.
pub fn vectors_dir() -> PathBuf {
    let ws = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let local = ws.join(".cache/3gpp");
    if local.join(".unpacked").exists() {
        return local;
    }
    let common = std::process::Command::new("git")
        .args(["rev-parse", "--path-format=absolute", "--git-common-dir"])
        .current_dir(&ws)
        .output()
        .ok()
        .filter(|o| o.status.success())
        .map(|o| PathBuf::from(String::from_utf8_lossy(&o.stdout).trim()));
    let main = common.and_then(|c| c.parent().map(|p| p.join(".cache/3gpp")));
    match main {
        Some(d) if d.join(".unpacked").exists() => d,
        _ => panic!("3GPP test vectors not found: run `just vectors`"),
    }
}

fn read_i16(path: &Path) -> Vec<i16> {
    let bytes = std::fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    bytes
        .as_chunks::<2>()
        .0
        .iter()
        .map(|c| i16::from_le_bytes(*c))
        .collect()
}

/// Maps the encoder's TX frame type (as stored in `.COD`) to the decoder's RX type,
/// as the reference decoder does for encoder-generated bitstreams.
fn tx_to_rx(tx: i16) -> i16 {
    match tx {
        0 => rx::SPEECH_GOOD,
        1 => rx::SID_FIRST,
        2 => rx::SID_UPDATE,
        3 => rx::NO_DATA,
        4 => rx::SPEECH_DEGRADED,
        5 => rx::SPEECH_BAD,
        6 => rx::SID_BAD,
        7 => rx::ONSET,
        _ => panic!("invalid TX frame type {tx}"),
    }
}

pub fn read_cod(path: &Path) -> Vec<Frame> {
    let w = read_i16(path);
    assert_eq!(
        w.len() % COD_WORDS,
        0,
        "{}: not a whole number of frames",
        path.display()
    );
    w.as_chunks::<COD_WORDS>()
        .0
        .iter()
        .map(|f| Frame {
            rx_type: tx_to_rx(f[0]),
            bits: f[1..=MAX_SERIAL_SIZE].try_into().unwrap(),
            mode: f[1 + MAX_SERIAL_SIZE],
        })
        .collect()
}

/// Cases in the SHA-pinned vector set, linear + A-law + µ-law (each law has a third). Asserted,
/// so a name filter can't silently drop a set again (it once skipped all of DTX2).
pub const DECODER_CASES: usize = 3 * 191;
pub const ENCODER_CASES: usize = 3 * 199;

/// All decoder cases (linear, A-law, µ-law), sorted by name.
/// Sets: NODTX (T00..T22 per mode + T21 mode switching), DTX (VAD1) and DTX2 (VAD2).
pub fn cases() -> Vec<Case> {
    let root = vectors_dir();
    let mut out = Vec::new();
    for set in ["nodtx", "dtx", "dtx2"] {
        let mut stack = vec![root.join(set)];
        while let Some(dir) = stack.pop() {
            for e in std::fs::read_dir(&dir).unwrap().flatten() {
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                    continue;
                }
                let stem = p.file_stem().unwrap().to_string_lossy().into_owned();
                if p.extension().is_some_and(|x| x == "COD")
                    && let Some(v) = parse_stem(&stem)
                {
                    let out_path = p.with_extension("OUT");
                    if out_path.exists() {
                        out.push(Case {
                            name: format!("{set}/{stem}"),
                            law: v.law,
                            frames: read_cod(&p),
                            expected: read_i16(&out_path),
                        });
                    }
                }
            }
        }
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    assert_eq!(out.len(), DECODER_CASES, "decoder conformance cases found");
    out
}

/// The companded cases opencore (the C oracle) gets wrong: its arithmetic wraps where TS 26.073
/// saturates (10.2/12.2 kbit/s on near-full-scale input; the 3GPP reference reproduces them, so
/// does amr-nb). Pinned for the C oracle's runs by their first mismatch: a changed mismatch
/// fails the run. RESULTS.md, defects found in the C.
pub const OPENCORE_DEFECTS: [&str; 10] = [
    "nodtx/T03A102: first mismatch frame 37 sample 82 (got 44, expected 47), rx_type 0 mode 6",
    "dtx/DTX2U122: first mismatch frame 461 word 90 (got 1, expected 0), tx_type 0 vs 0",
    "dtx/DTX3U122: first mismatch frame 1255 word 91 (got 1, expected 0), tx_type 0 vs 0",
    "nodtx/T01U122: first mismatch frame 20 word 91 (got 1, expected 0), tx_type 0 vs 0",
    "nodtx/T03A102: first mismatch frame 37 word 155 (got 0, expected 1), tx_type 0 vs 0",
    "nodtx/T03A122: first mismatch frame 121 word 193 (got 1, expected 0), tx_type 0 vs 0",
    "nodtx/T03U122: first mismatch frame 39 word 207 (got 0, expected 1), tx_type 0 vs 0",
    "nodtx/T05A122: first mismatch frame 179 word 140 (got 1, expected 0), tx_type 0 vs 0",
    "nodtx/T15U122: first mismatch frame 49 word 141 (got 1, expected 0), tx_type 0 vs 0",
    "nodtx/T20U102: first mismatch frame 563 word 163 (got 0, expected 1), tx_type 0 vs 0",
];

/// Splits case results into (bit-exact flags, pinned mismatches, failure messages): a pinned
/// case must still fail with exactly its pinned message. `kind` is `" sample "` for decoder
/// messages, `" word "` for encoder ones.
fn triage<'a>(
    kind: &str,
    pins: &[&str],
    names: impl Iterator<Item = &'a str>,
    results: &[Result<(), String>],
) -> (Vec<bool>, usize, Vec<String>) {
    let (mut known, mut failures) = (0, Vec::new());
    let exact = names
        .zip(results)
        .map(|(name, r)| {
            let pinned = pins
                .iter()
                .find(|d| d.starts_with(&format!("{name}: ")) && d.contains(kind));
            match (r, pinned) {
                (Ok(()), None) => return true,
                (Err(e), Some(d)) if e == d => known += 1,
                (Ok(()), Some(d)) => failures.push(format!("now bit-exact, unpin: {d}")),
                (Err(e), _) => failures.push(e.clone()),
            }
            false
        })
        .collect();
    (exact, known, failures)
}

/// A parsed `.COD` vector name.
struct VectorName {
    /// Encoder input stem (`T08A` for `T08A795`).
    input: String,
    /// Linear input stem, which names the `.MOD` file (`T21` for `T21A`).
    linear: String,
    law: Law,
    /// Fixed mode index; `None` for the T21 mode-switching case.
    mode: Option<usize>,
}

/// Vector names: `T<nn>[_<mode>]`, `DTX<n>[_<mode>]` (VAD1) or `DT<nn>_<mode>` (VAD2);
/// A-law/µ-law ones put `A`/`U` in place of the `_` (e.g. `T08A795`, `DT22U475`, `T21A`).
fn parse_stem(stem: &str) -> Option<VectorName> {
    let prefix = ["DTX", "DT", "T"]
        .into_iter()
        .find(|p| stem.starts_with(p))?;
    let digits = stem[prefix.len()..]
        .find(|c: char| !c.is_ascii_digit())
        .unwrap_or(stem.len() - prefix.len());
    let (linear, rest) = stem.split_at(prefix.len() + digits);
    if digits == 0 {
        return None;
    }
    let (law, mode) = match rest.chars().next() {
        None => (Law::Linear, ""),
        Some('_') => (Law::Linear, &rest[1..]),
        Some('A') => (Law::A, &rest[1..]),
        Some('U') => (Law::Mu, &rest[1..]),
        _ => return None,
    };
    let input = match law {
        Law::Linear => linear.to_string(),
        _ => stem[..linear.len() + 1].to_string(),
    };
    let mode = match mode {
        "" => None,
        m => Some(MODE_NAMES.iter().position(|n| *n == m)?),
    };
    Some(VectorName {
        input,
        linear: linear.to_string(),
        law,
        mode,
    })
}

/// The one summary line per run, e.g. `rust: 572/573 cases bit-exact (252063 frames; linear
/// 191/191, A-law 190/191, µ-law 191/191; 1 pinned opencore defect)`.
fn report(
    label: &str,
    laws: impl Iterator<Item = Law>,
    exact: &[bool],
    known: usize,
    frames: usize,
) {
    let mut n = [(0, 0); 3];
    for (law, &ok) in laws.zip(exact) {
        let e = &mut n[law as usize];
        e.0 += ok as usize;
        e.1 += 1;
    }
    let by_law = ["linear", "A-law", "µ-law"]
        .iter()
        .zip(n)
        .filter(|(_, (_, t))| *t > 0)
        .map(|(l, (p, t))| format!("{l} {p}/{t}"))
        .collect::<Vec<_>>()
        .join(", ");
    let known = match known {
        0 => String::new(),
        1 => "; 1 pinned opencore defect".into(),
        k => format!("; {k} pinned opencore defects"),
    };
    eprintln!(
        "{label}: {}/{} cases bit-exact ({frames} frames; {by_law}{known})",
        exact.iter().filter(|&&e| e).count(),
        exact.len()
    );
}

/// Compares decoder output (linear; compressed here for A-law/µ-law cases) with the expected
/// `.OUT`; returns a one-line diagnosis on mismatch.
pub fn compare(case: &Case, got: &[i16]) -> Result<(), String> {
    let got: Vec<i16> = got.iter().map(|&s| case.law.compress(s)).collect();
    if got.len() != case.expected.len() {
        return Err(format!(
            "{}: {} samples, expected {}",
            case.name,
            got.len(),
            case.expected.len()
        ));
    }
    match got.iter().zip(&case.expected).position(|(a, b)| a != b) {
        None => Ok(()),
        Some(i) => {
            let f = i / L_FRAME;
            Err(format!(
                "{}: first mismatch frame {f} sample {} (got {}, expected {}), rx_type {} mode {}",
                case.name,
                i % L_FRAME,
                got[i],
                case.expected[i],
                case.frames[f].rx_type,
                case.frames[f].mode
            ))
        }
    }
}

/// Runs `decode` over every case; panics with a summary listing each failing case (`pins`:
/// mismatches expected, as `OPENCORE_DEFECTS`).
pub fn run_all(label: &str, pins: &[&str], mut decode: impl FnMut(&Case) -> Vec<i16>) {
    let cases = cases();
    assert!(!cases.is_empty(), "no conformance cases found");
    let results: Vec<_> = cases.iter().map(|c| compare(c, &decode(c))).collect();
    let (exact, known, failures) =
        triage(" sample ", pins, cases.iter().map(|c| &*c.name), &results);
    report(
        label,
        cases.iter().map(|c| c.law),
        &exact,
        known,
        cases.iter().map(|c| c.frames.len()).sum(),
    );
    assert!(
        failures.is_empty(),
        "{label}: {} failing cases:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Words per encoder output frame in a `.COD` file: `[tx_type, 244 bits, mode, 0, 0, 0, 0]`.
pub const COD_FRAME: usize = COD_WORDS;

/// One `.INP` → `.COD` encoder conformance case.
pub struct EncCase {
    pub name: String,
    /// DTX: 0 off, 1 with VAD option 1 (`DTX*`), 2 with VAD option 2 (`DT<nn>_<mode>`).
    pub vad: u8,
    pub law: Law,
    /// Input PCM, linear (A-law/µ-law files expanded), a whole number of frames.
    pub input: Vec<i16>,
    /// Requested mode (0..=7) per frame.
    pub modes: Vec<i16>,
    /// Expected `.COD` frames.
    pub expected: Vec<[i16; COD_FRAME]>,
}

const MODE_NAMES: [&str; 8] = ["475", "515", "59", "67", "74", "795", "102", "122"];

/// Encoder cases, each in linear, A-law and µ-law: NODTX `T<nn>_<mode>` (with `T<nn>.INP`),
/// T21 mode switching (with `T21.MOD`), DTX `DTX<n>_<mode>` (VAD1, with `DTX<n>.INP`) and DTX2
/// `DT<nn>_<mode>` (VAD2, with `DT<nn>.INP`), sorted by name.
pub fn encoder_cases() -> Vec<EncCase> {
    encoder_cases_where(|_| true)
}

/// The encoder cases whose name (`<set>/<stem>`) `keep` selects; only those files are read.
pub fn encoder_cases_where(keep: impl Fn(&str) -> bool) -> Vec<EncCase> {
    let root = vectors_dir();
    // Keyed by (set, stem): a stem may recur in another set or subdirectory.
    let mut inputs = std::collections::HashMap::<(&str, String), PathBuf>::new();
    let mut cods = Vec::new();
    for set in ["nodtx", "dtx", "dtx2"] {
        let mut stack = vec![root.join(set)];
        while let Some(dir) = stack.pop() {
            for e in std::fs::read_dir(&dir).unwrap().flatten() {
                let p = e.path();
                if p.is_dir() {
                    stack.push(p);
                    continue;
                }
                let stem = p.file_stem().unwrap().to_string_lossy().into_owned();
                match p.extension().and_then(|x| x.to_str()) {
                    Some("INP") => {
                        if let Some(prev) = inputs.insert((set, stem), p.clone()) {
                            assert!(
                                read_i16(&prev) == read_i16(&p),
                                "{} and {} differ",
                                prev.display(),
                                p.display()
                            );
                        }
                    }
                    Some("COD") => {
                        if let Some(v) = parse_stem(&stem) {
                            cods.push((set, stem, v, p.clone()))
                        }
                    }
                    _ => {}
                }
            }
        }
    }
    assert_eq!(cods.len(), ENCODER_CASES, "encoder conformance cases found");
    let mut out = Vec::new();
    for (set, stem, v, cod) in cods {
        if !keep(&format!("{set}/{stem}")) {
            continue;
        }
        let inp = &inputs[&(set, v.input)];
        let input: Vec<i16> = read_i16(inp).into_iter().map(|s| v.law.expand(s)).collect();
        let frames = input.len() / L_FRAME;
        let modes = match v.mode {
            Some(m) => vec![m as i16; frames],
            None => {
                let path = inp.with_file_name(format!("{}.MOD", v.linear));
                let text = std::fs::read_to_string(&path)
                    .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
                text.split_whitespace()
                    .map(|t| {
                        MODE_NAMES
                            .iter()
                            .position(|n| t == format!("MR{n}"))
                            .expect("mode name") as i16
                    })
                    .collect()
            }
        };
        let w = read_i16(&cod);
        out.push(EncCase {
            name: format!("{set}/{stem}"),
            vad: match set {
                "dtx" => 1,
                "dtx2" => 2,
                _ => 0,
            },
            law: v.law,
            input,
            modes,
            expected: w.as_chunks::<COD_FRAME>().0.to_vec(),
        });
    }
    out.sort_by(|a, b| a.name.cmp(&b.name));
    out
}

/// Runs `encode` over every encoder case `keep` selects (it returns the `.COD` frames); panics
/// with a summary listing each failing case and its first mismatching frame.
pub fn run_all_enc(
    label: &str,
    pins: &[&str],
    keep: impl Fn(&EncCase) -> bool,
    mut encode: impl FnMut(&EncCase) -> Vec<[i16; COD_FRAME]>,
) {
    let cases: Vec<EncCase> = encoder_cases().into_iter().filter(keep).collect();
    assert!(!cases.is_empty(), "no encoder conformance cases found");
    let results: Vec<_> = cases
        .iter()
        .map(|c| compare_cod(&c.name, &encode(c), &c.expected))
        .collect();
    let (exact, known, failures) = triage(" word ", pins, cases.iter().map(|c| &*c.name), &results);
    report(
        label,
        cases.iter().map(|c| c.law),
        &exact,
        known,
        cases.iter().map(|c| c.expected.len()).sum(),
    );
    assert!(
        failures.is_empty(),
        "{label}: {} failing cases:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// Compares encoder output with the expected `.COD` frames; a one-line diagnosis on mismatch.
fn compare_cod(
    name: &str,
    got: &[[i16; COD_FRAME]],
    expected: &[[i16; COD_FRAME]],
) -> Result<(), String> {
    if got.len() != expected.len() {
        return Err(format!(
            "{name}: {} frames, expected {}",
            got.len(),
            expected.len()
        ));
    }
    let Some(f) = got.iter().zip(expected).position(|(a, b)| a != b) else {
        return Ok(());
    };
    let w = got[f]
        .iter()
        .zip(&expected[f])
        .position(|(a, b)| a != b)
        .unwrap();
    Err(format!(
        "{name}: first mismatch frame {f} word {w} (got {}, expected {}), tx_type {} vs {}",
        got[f][w], expected[f][w], got[f][0], expected[f][0]
    ))
}

/// Frame-synchronisation cases (TS 26.074 §8.2) in the pinned set: 8 modes × 160 shifts × 3 laws.
pub const SYNC_CASES: usize = 3 * 8 * 160;

/// Frame-sync cases per law (linear, A, µ) and mode (4.75 … 12.2), 160 each: all bit-exact.
pub const SYNC_ALL: [[usize; 8]; 3] = [[160; 8]; 3];

/// opencore's bit-exact frame-sync cases per law and mode (`SYNC_ALL` for amr-nb and the 3GPP
/// reference): its arithmetic wraps on the full-scale sync sequences. Pinned: a change fails.
pub const OPENCORE_SYNC_EXACT: [[usize; 8]; 3] = [
    [129, 145, 91, 90, 85, 69, 57, 72],
    [129, 145, 91, 90, 85, 69, 57, 72],
    [134, 144, 98, 91, 91, 76, 64, 68],
];

/// One frame-synchronisation case: `SEQS*.INP` (3 homing frames, then two identical frames)
/// delayed by `shift` zero samples; `expected` is the second encoder output frame after the
/// last decoder homing frame (`S<shift><law><mode>.COD`).
pub struct SyncCase {
    pub name: String,
    pub law: Law,
    pub mode: i16,
    /// Linear input, padded with zeros to whole frames.
    pub input: Vec<i16>,
    pub expected: [i16; COD_FRAME],
}

fn sync_dir() -> PathBuf {
    vectors_dir().join("nodtx")
}

const LAWS: [(Law, &str); 3] = [(Law::Linear, "_"), (Law::A, "A"), (Law::Mu, "U")];

pub fn sync_cases() -> Vec<SyncCase> {
    let dir = sync_dir();
    let mut out = Vec::new();
    for (law, tag) in LAWS {
        for (mode, name) in MODE_NAMES.iter().enumerate() {
            // 12.2 has its own sequence: no lookahead, so its 3rd homing frame ends early.
            let seq = match (mode, tag) {
                (7, "_") => "SEQS_122".to_string(),
                (7, t) => format!("SEQS{t}122"),
                (_, "_") => "SEQS".to_string(),
                (_, t) => format!("SEQS{t}"),
            };
            let seq: Vec<i16> = read_i16(&dir.join(format!("SYNC_INP/{seq}.INP")))
                .into_iter()
                .map(|s| law.expand(s))
                .collect();
            for shift in 0..L_FRAME {
                let stem = format!("S{shift:03}{tag}{name}");
                let w = read_i16(&dir.join(format!("SYNC_{name}/{stem}.COD")));
                let mut input = vec![0; shift];
                input.extend_from_slice(&seq);
                input.resize(input.len().div_ceil(L_FRAME) * L_FRAME, 0);
                out.push(SyncCase {
                    name: format!("sync/{stem}"),
                    law,
                    mode: mode as i16,
                    input,
                    expected: w.try_into().expect("one .COD frame"),
                });
            }
        }
    }
    assert_eq!(out.len(), SYNC_CASES, "sync cases found");
    out
}

/// The frame the sync procedure records: the 2nd output frame after the last decoder homing
/// frame (`dhf`), if any.
pub fn sync_frame(out: &[[i16; COD_FRAME]], dhf: &[i16; COD_FRAME]) -> Option<[i16; COD_FRAME]> {
    let last = out.iter().rposition(|f| f == dhf)?;
    out.get(last + 2).copied()
}

/// Runs `encode(mode, linear input)` (fresh encoder, DTX off) over every frame-synchronisation
/// case (bit-exact counts per law and mode must be `expect`), and the bit-synchronisation check
/// (§8.1): of the 13 bit-shifted homing triplets in `BITS*.INP`, exactly the unshifted one makes
/// the encoder emit a decoder homing frame.
pub fn run_all_sync(
    label: &str,
    expect: &[[usize; 8]; 3],
    mut encode: impl FnMut(i16, &[i16]) -> Vec<[i16; COD_FRAME]>,
) {
    let dhf: Vec<_> = (0..8).map(|m| encode(m, &[0x0008; L_FRAME])[0]).collect();
    let cases = sync_cases();
    let mut failures = Vec::new();
    let mut counts = [[0; 8]; 3];
    for c in &cases {
        let got = sync_frame(&encode(c.mode, &c.input), &dhf[c.mode as usize]);
        if got.is_none() {
            failures.push(format!("{}: no 2nd frame after a homing frame", c.name));
        }
        counts[c.law as usize][c.mode as usize] += (got == Some(c.expected)) as usize;
    }
    let exact: usize = counts.iter().flatten().sum();
    let pinned = if exact < cases.len() {
        " (as pinned)"
    } else {
        ""
    };
    eprintln!(
        "{label} frame sync: {exact}/{} cases bit-exact{pinned} ({} frames)",
        cases.len(),
        cases.len(),
    );
    if counts != *expect {
        failures.push(format!(
            "frame sync bit-exact counts {counts:?}, expected {expect:?}"
        ));
    }

    let mut found = Vec::new(); // (law, the unshifted triplet alone found)
    for (law, tag) in LAWS {
        let file = if tag == "_" {
            "BITS".into()
        } else {
            format!("BITS{tag}")
        };
        let bits: Vec<i16> = read_i16(&sync_dir().join(format!("SYNC_INP/{file}.INP")))
            .into_iter()
            .map(|s| law.expand(s))
            .collect();
        for mode in 0..8 {
            let hits: Vec<usize> = bits
                .chunks(3 * L_FRAME)
                .enumerate()
                .filter(|(_, t)| encode(mode, t).contains(&dhf[mode as usize]))
                .map(|(i, _)| i)
                .collect();
            if hits != [0] {
                failures.push(format!("bit sync {file} mode {mode}: homing frame after triplets {hits:?}, expected [0]"));
            }
            found.push((law, hits == [0]));
        }
    }
    let ok: Vec<bool> = found.iter().map(|f| f.1).collect();
    report(
        &format!("{label} bit sync"),
        found.iter().map(|f| f.0),
        &ok,
        0,
        3 * 8 * 13 * 3,
    );
    assert!(
        failures.is_empty(),
        "{label}: {} failing sync cases:\n{}",
        failures.len(),
        failures.join("\n")
    );
}
