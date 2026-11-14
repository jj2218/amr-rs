//! Encode raw PCM to an `.amr` file, or decode one back.
//!
//! ```text
//! cargo run --release --example amr -- encode speech.raw speech.amr [mode] [dtx]
//! cargo run --release --example amr -- decode speech.amr speech.raw
//! ```
//!
//! Raw PCM is 8 kHz, 16-bit little-endian mono, e.g. from
//! `ffmpeg -i in.wav -f s16le -ar 8000 -ac 1 speech.raw`. `mode` is the bit rate in kbit/s
//! (4.75, 5.15, 5.9, 6.7, 7.4, 7.95, 10.2, 12.2; default 12.2); `dtx` turns on discontinuous
//! transmission with VAD option 1. A final partial frame is padded with silence. Decoding
//! conceals frames of other codecs (types 9–11) like lost ones.

use amr_nb::{Decoder, Dtx, Encoder, MAX_SERIAL_SIZE, Mode, RxFrameType, rfc4867};
use std::error::Error;

fn mode(rate: &str) -> Option<Mode> {
    Some(match rate {
        "4.75" => Mode::MR475,
        "5.15" => Mode::MR515,
        "5.9" => Mode::MR59,
        "6.7" => Mode::MR67,
        "7.4" => Mode::MR74,
        "7.95" => Mode::MR795,
        "10.2" => Mode::MR102,
        "12.2" => Mode::MR122,
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
        .chunks(160)
        .map(|chunk| {
            let mut pcm = [0i16; 160];
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
        let mut pcm = [0i16; 160];
        match frame {
            Ok(f) => dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm),
            // A frame of another codec (types 9..=11): conceal it, keep the timing.
            Err(rfc4867::Error::UnsupportedFrameType(_)) => {
                let bits = [0; MAX_SERIAL_SIZE];
                dec.decode(RxFrameType::NoData, Mode::MR475, &bits, &mut pcm);
            }
            Err(e) => return Err(e.into()), // truncated: the stream ends
        }
        out.extend(pcm.iter().flat_map(|s| s.to_le_bytes()));
        n += 1;
    }
    std::fs::write(output, out)?;
    Ok(n)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let usage = "usage: amr encode <in.raw> <out.amr> [mode] [dtx] | amr decode <in.amr> <out.raw>";
    let frames = match args.iter().map(String::as_str).collect::<Vec<_>>()[..] {
        ["encode", input, output, ref rest @ ..] => {
            let (mut rate, mut dtx) = (Mode::MR122, Dtx::Off);
            for &opt in rest {
                match (opt, mode(opt)) {
                    ("dtx", _) => dtx = Dtx::Vad1,
                    (_, Some(m)) => rate = m,
                    _ => return Err(format!("unknown option {opt}: {usage}").into()),
                }
            }
            encode(input, output, rate, dtx)?
        }
        ["decode", input, output] => decode(input, output)?,
        _ => return Err(usage.into()),
    };
    eprintln!("{frames} frames ({:.2} s)", frames as f64 * 0.02);
    Ok(())
}
