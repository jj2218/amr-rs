//! Shared conformance-test plumbing: locating and parsing the 3GPP TS 26.174 test vectors.
#![allow(dead_code)]

pub mod stress;

use amr_wb::{MAX_SERIAL_SIZE, Mode, RxFrameType};
use std::path::{Path, PathBuf};

/// Serial bits per frame type: modes 0..=8, then SID (`nb_of_bits`).
const NB_BITS: [usize; 10] = [132, 177, 253, 285, 317, 365, 397, 461, 477, 35];
/// `TXRXFLAG` values of the default (TS 26.173) bitstream format.
const TX_FRAME_TYPE: i16 = 0x6b21;
const RX_FRAME_TYPE: i16 = 0x6b20;

/// One decoder input frame.
#[derive(Clone)]
pub struct Frame {
    pub frame_type: RxFrameType,
    pub mode: Mode,
    /// The frame's bits, 0/1 (the rest 0).
    pub bits: [i16; MAX_SERIAL_SIZE],
    /// The same bits as the C decoders read them (BIT_0 = -127, BIT_1 = 127).
    pub prms: [i16; MAX_SERIAL_SIZE],
}

/// One `.cod` → `.out` conformance case.
pub struct Case {
    pub name: String,
    pub frames: Vec<Frame>,
    pub expected: Vec<i16>,
}

/// `.cache/3gpp/ts26174` of this checkout, or of the main checkout when run from a worktree;
/// `TS26174_DIR` overrides it (an older vector version, unpacked alike, to attribute a CR).
pub fn vectors_dir() -> PathBuf {
    if let Some(d) = std::env::var_os("TS26174_DIR") {
        return PathBuf::from(d);
    }
    let ws = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let local = ws.join(".cache/3gpp/ts26174");
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
    let main = common.and_then(|c| c.parent().map(|p| p.join(".cache/3gpp/ts26174")));
    match main {
        Some(d) if d.join(".unpacked").exists() => d,
        _ => panic!("3GPP TS 26.174 test vectors not found: run `just vectors`"),
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

fn rx_type(n: i16) -> RxFrameType {
    use RxFrameType::*;
    [
        SpeechGood,
        SpeechProbablyDegraded,
        SpeechLost,
        SpeechBad,
        SidFirst,
        SidUpdate,
        SidBad,
        NoData,
    ][n as usize]
}

/// Parses a `.cod` file in the default format: per frame `TXRXFLAG`, frame type, mode, then
/// one word per bit (±127). TX frame types map to RX types as the reference decoder does.
pub fn read_cod(path: &Path) -> Vec<Frame> {
    let w = read_i16(path);
    let mut frames = Vec::new();
    let mut i = 0;
    while i + 3 <= w.len() {
        let (flag, ft, mode) = (w[i], w[i + 1], w[i + 2]);
        let frame_type = match flag {
            TX_FRAME_TYPE => rx_type([0, 4, 5, 7][ft as usize]),
            RX_FRAME_TYPE => rx_type(ft),
            _ => panic!("{}: bad TXRXFLAG {flag:#x} at word {i}", path.display()),
        };
        use RxFrameType::*;
        let n = match frame_type {
            SidFirst | SidUpdate | SidBad | NoData => NB_BITS[9],
            _ => NB_BITS[mode as usize],
        };
        let raw = &w[i + 3..i + 3 + n];
        let mut f = Frame {
            frame_type,
            mode: Mode::from_index(mode as usize).expect("mode"),
            bits: [0; MAX_SERIAL_SIZE],
            prms: [-127; MAX_SERIAL_SIZE],
        };
        for (k, &b) in raw.iter().enumerate() {
            f.bits[k] = (b == 127) as i16;
            f.prms[k] = b;
        }
        frames.push(f);
        i += 3 + n;
    }
    assert_eq!(i, w.len(), "{}: trailing words", path.display());
    frames
}

/// Decoder cases in the SHA-pinned TS 26.174 set (j00): T00..T20 per mode, T22 mode switching,
/// Dtx1..Dtx5. Asserted, so a filter can't silently drop a set.
pub const DECODER_CASES: usize = 221;

/// All decoder cases (`.cod` with a matching `.out`), sorted by name.
pub fn cases() -> Vec<Case> {
    let root = vectors_dir().join("T");
    let mut cods = Vec::new();
    let mut stack = vec![root.clone()];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                stack.push(p);
            } else if p.extension().is_some_and(|x| x.eq_ignore_ascii_case("cod"))
                && p.with_extension("out").exists()
            {
                cods.push(p);
            }
        }
    }
    cods.sort();
    let cases: Vec<Case> = cods
        .iter()
        .map(|p| Case {
            name: p.strip_prefix(&root).unwrap().display().to_string(),
            frames: read_cod(p),
            expected: read_i16(&p.with_extension("out")),
        })
        .collect();
    assert_eq!(cases.len(), DECODER_CASES, "decoder vector count");
    cases
}

/// Mode directory suffixes in mode order (`T_660` .. `T_2385`).
pub const MODE_NAMES: [&str; 9] = [
    "660", "885", "1265", "1425", "1585", "1825", "1985", "2305", "2385",
];

/// One `.inp` → `.cod` encoder conformance case.
#[derive(Clone)]
pub struct EncCase {
    pub name: String,
    pub input: Vec<i16>,
    /// The mode per frame (the last one repeats, as the reference's mode file reader does).
    pub modes: Vec<u8>,
    pub dtx: bool,
    /// Frame sync: `expected` is the 2nd output frame after the last decoder homing frame.
    pub sync: bool,
    /// The `.cod` words: per frame 0x6b21, TX frame type, mode, one ±127 word per bit.
    pub expected: Vec<i16>,
}

impl EncCase {
    pub fn mode(&self, frame: usize) -> u8 {
        self.modes[frame.min(self.modes.len() - 1)]
    }
}

/// Encoder cases in the SHA-pinned TS 26.174 set (j00): the 221 decoder `.cod` files (T00..T23
/// per mode but T22, T22 mode switching, Dtx1..Dtx5) and 320 frame-sync streams per mode.
pub const ENCODER_CASES: usize = 221 + 320 * 9;

/// The `.inp` (and mode list) for a `.cod` path relative to `T/`.
fn enc_case(root: &Path, rel: &Path) -> EncCase {
    let stem = rel.file_stem().unwrap().to_str().unwrap();
    let dir = rel.parent().unwrap();
    let dtx = dir.starts_with("DTX_test_vectors");
    let (input, modes) = if stem == "T22" {
        let m = std::fs::read_to_string(root.join("T_inp/T22.mod")).unwrap();
        (
            root.join("T_inp/T22.inp"),
            m.split_whitespace().map(|v| v.parse().unwrap()).collect(),
        )
    } else {
        let (seq, mode) = stem.split_once('_').unwrap();
        let mode = MODE_NAMES.iter().position(|&n| n == mode).unwrap() as u8;
        let inp = dir
            .parent()
            .unwrap()
            .join("T_inp")
            .join(format!("{seq}.inp"));
        (root.join(inp), vec![mode])
    };
    EncCase {
        name: rel.display().to_string(),
        input: read_i16(&input),
        modes,
        dtx,
        sync: false,
        expected: read_i16(&root.join(rel)),
    }
}

/// All encoder cases, sorted by name: `T/**/*.cod` with an input, then `S/Sync_<mode>`
/// (`syncNNN`: `seqsync.inp`, homing frames then speech, delayed by NNN zero samples).
pub fn encoder_cases() -> Vec<EncCase> {
    let root = vectors_dir().join("T");
    let mut cods = Vec::new();
    let mut stack = vec![root.clone()];
    while let Some(d) = stack.pop() {
        for e in std::fs::read_dir(&d).unwrap() {
            let p = e.unwrap().path();
            if p.is_dir() {
                stack.push(p);
            } else if p.extension().is_some_and(|x| x.eq_ignore_ascii_case("cod")) {
                cods.push(p.strip_prefix(&root).unwrap().to_path_buf());
            }
        }
    }
    cods.sort();
    let mut cases: Vec<EncCase> = cods.iter().map(|r| enc_case(&root, r)).collect();
    let s = vectors_dir().join("S");
    let seq = read_i16(&s.join("Sync_inp/seqsync.inp"));
    for (m, name) in MODE_NAMES.iter().enumerate() {
        for k in 0..320 {
            let rel = format!("Sync_{name}/sync{k:03}_{name}.cod");
            let mut input = vec![0; k];
            input.extend_from_slice(&seq);
            input.resize(input.len().div_ceil(320) * 320, 0);
            cases.push(EncCase {
                name: rel.clone(),
                input,
                modes: vec![m as u8],
                dtx: false,
                sync: true,
                expected: read_i16(&s.join(&rel)),
            });
        }
    }
    assert_eq!(cases.len(), ENCODER_CASES, "encoder vector count");
    cases
}

/// `S/Sync_inp/bitsync.inp`: bit-shifted homing-frame triplets (3 frames each).
pub fn read_bitsync() -> Vec<i16> {
    let w = read_i16(&vectors_dir().join("S/Sync_inp/bitsync.inp"));
    assert_eq!(w.len() % 960, 0, "bitsync.inp: whole triplets");
    w
}

/// An encoder frame as a `.cod` vector stores it (the reference's default format): 0x6b21, the
/// TX frame type (`TX_*`), the requested mode, then one ±127 word per bit (a SID's bits, also
/// for SID_FIRST and NO_DATA frames, as the reference writes them).
#[allow(dead_code)] // not every test binary encodes
pub fn cod_words(f: &amr_wb::rfc4867::Frame) -> Vec<i16> {
    let (tx, n) = match f.frame_type {
        RxFrameType::SpeechGood => (0, f.mode.bits()),
        RxFrameType::SidFirst => (1, amr_wb::SID_BITS),
        RxFrameType::SidUpdate => (2, amr_wb::SID_BITS),
        RxFrameType::NoData => (3, amr_wb::SID_BITS),
        t => panic!("not an encoder frame type: {t:?}"),
    };
    let mut w = vec![0x6b21, tx, f.mode as i16];
    w.extend(f.bits[..n].iter().map(|&b| if b != 0 { 127 } else { -127 }));
    w
}
