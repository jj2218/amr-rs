//! Pure-Rust, bit-exact AMR-NB (3GPP TS 26.073) speech codec: encoder, decoder, DTX with
//! VAD option 1 or 2, and RFC 4867 packing (`.amr` files and RTP payloads).
//!
//! Bit-exact on the official 3GPP TS 26.074 test vectors (linear, A-law, µ-law, bit and
//! frame synchronisation) and differentially tested against the TS 26.073 reference C. Ported
//! from opencore-amr 0.1.6 (Apache-2.0); see NOTICE.
//!
//! Audio is 8 kHz, 16-bit mono PCM in frames of 160 samples (20 ms). The [`Encoder`] turns a
//! frame into an [`rfc4867::Frame`] (frame type, mode and serial bits); the [`Decoder`] turns
//! one back into PCM; [`rfc4867`] converts frames to and from bytes.
//!
//! ```
//! use amr_nb::{Decoder, Dtx, Encoder, Mode, rfc4867};
//!
//! // Encode one second of a 440 Hz tone at 12.2 kbit/s and store it as an `.amr` file.
//! let pcm: Vec<i16> = (0..8000)
//!     .map(|n| (8000.0 * (n as f64 * 440.0 * std::f64::consts::TAU / 8000.0).sin()) as i16)
//!     .collect();
//! let mut enc = Encoder::new(Dtx::Off);
//! let frames: Vec<rfc4867::Frame> = pcm
//!     .chunks_exact(160)
//!     .map(|chunk| enc.encode(Mode::MR122, chunk.try_into().unwrap()))
//!     .collect();
//! let file = rfc4867::write_storage(&frames)?;
//! assert_eq!(file.len(), 6 + 50 * 32); // magic, then 32 octets per 12.2 kbit/s frame
//!
//! // Read it back and decode.
//! let mut dec = Decoder::new();
//! let mut out = Vec::new();
//! for frame in rfc4867::read_storage(&file)? {
//!     let f = frame?;
//!     let mut samples = [0i16; 160];
//!     dec.decode(f.frame_type, f.mode, &f.bits, &mut samples);
//!     out.extend_from_slice(&samples);
//! }
//! assert_eq!(out.len(), pcm.len());
//! # Ok::<(), rfc4867::Error>(())
//! ```

#![deny(missing_docs)]

mod basic_op;
mod cnst;
mod math;

// Codec modules, ported from opencore.
mod decoder;
mod dtx;
mod encoder;
mod excitation;
mod gains;
mod lpc;
mod post;
mod simd;
// Stage timers: public only in `--cfg codec_stages` builds (`BENCH=stages just bench`).
#[cfg(codec_stages)]
#[doc(hidden)]
pub mod stages;
#[cfg(not(codec_stages))]
#[allow(dead_code)] // the names and `take` are read only by the stages bench
mod stages;
mod tables;

pub mod rfc4867;

pub use cnst::Mode;
pub use decoder::{Decoder, MAX_SERIAL_SIZE, RxFrameType};
pub use encoder::{Dtx, Encoder};

#[cfg(test)]
mod testutil;
