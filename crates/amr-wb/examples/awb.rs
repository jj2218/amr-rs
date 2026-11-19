//! Encode raw PCM to an `.awb` file, or decode one back.
//!
//! ```text
//! cargo run --release --example awb -- encode speech.raw speech.awb [mode] [dtx]
//! cargo run --release --example awb -- decode speech.awb speech.raw
//! ```
//!
//! Raw PCM is 16 kHz, 16-bit little-endian mono, e.g. from
//! `ffmpeg -i in.wav -f s16le -ar 16000 -ac 1 speech.raw`. `mode` is the bit rate in kbit/s
//! (6.60, 8.85, 12.65, 14.25, 15.85, 18.25, 19.85, 23.05, 23.85; default 12.65); `dtx` turns
//! on discontinuous transmission. A final partial frame is padded with silence. Decoding stops
//! at a truncated frame or a reserved frame type (10–13), as RFC 4867 has a file end there.

use amr_wb::{Decoder, Dtx, Encoder, L_FRAME, Mode, rfc4867};
use std::error::Error;

fn mode(rate: &str) -> Option<Mode> {
    Some(match rate {
        "6.60" | "6.6" => Mode::MD66,
        "8.85" => Mode::MD885,
        "12.65" => Mode::MD1265,
        "14.25" => Mode::MD1425,
        "15.85" => Mode::MD1585,
        "18.25" => Mode::MD1825,
        "19.85" => Mode::MD1985,
        "23.05" => Mode::MD2305,
        "23.85" => Mode::MD2385,
        _ => return None,
    })
}

fn encode(input: &str, output: &str, mode: Mode, dtx: Dtx) -> Result<usize, Box<dyn Error>> {
    let bytes = std::fs::read(input)?;
    let samples: Vec<i16> = bytes
        .as_chunks::<2>()
        .0
        .iter()
        .map(|&b| i16::from_le_bytes(b))
        .collect();
    let mut enc = Encoder::new(dtx);
    let frames: Vec<rfc4867::Frame> = samples
        .chunks(L_FRAME)
        .map(|chunk| {
            let mut pcm = [0i16; L_FRAME];
            pcm[..chunk.len()].copy_from_slice(chunk);
            enc.encode(mode, &pcm)
        })
        .collect();
    std::fs::write(output, rfc4867::write_storage(&frames)?)?;
    Ok(frames.len())
}

fn decode(input: &str, output: &str) -> Result<usize, Box<dyn Error>> {
    let file = std::fs::read(input)?;
    let mut dec = Decoder::new();
    let mut out = Vec::new();
    let mut n = 0;
    for frame in rfc4867::read_storage(&file)? {
        let f = frame?; // truncated or a reserved frame type: the stream ends
        let mut pcm = [0i16; L_FRAME];
        dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
        out.extend(pcm.iter().flat_map(|s| s.to_le_bytes()));
        n += 1;
    }
    std::fs::write(output, out)?;
    Ok(n)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let frames = run(&args.iter().map(String::as_str).collect::<Vec<_>>())?;
    eprintln!("{frames} frames ({:.2} s)", frames as f64 * 0.02);
    Ok(())
}

/// The command line without the program name; the number of frames written. Public for the
/// conformance test, which runs every encoder vector through it (`tests/enc_conformance.rs`).
pub fn run(args: &[&str]) -> Result<usize, Box<dyn Error>> {
    let usage = "usage: awb encode <in.raw> <out.awb> [mode] [dtx] | awb decode <in.awb> <out.raw>";
    Ok(match *args {
        ["encode", input, output, ref rest @ ..] => {
            let (mut rate, mut dtx) = (Mode::MD1265, Dtx::Off);
            for &opt in rest {
                match (opt, mode(opt)) {
                    ("dtx", _) => dtx = Dtx::On,
                    (_, Some(m)) => rate = m,
                    _ => return Err(format!("unknown option {opt}: {usage}").into()),
                }
            }
            encode(input, output, rate, dtx)?
        }
        ["decode", input, output] => decode(input, output)?,
        _ => return Err(usage.into()),
    })
}
