//! Pure-Rust, bit-exact AMR-WB (3GPP TS 26.173, ITU-T G.722.2) speech decoder and encoder,
//! with RFC 4867 packing (`.awb` files and RTP payloads).
//!
//! Bit-exact on the official 3GPP TS 26.174 decoder and encoder test vectors and
//! differentially tested against the TS 26.173 reference C. The decoder is ported from
//! opencore-amr 0.1.6, the encoder from vo-amrwbenc 0.1.3 (both Apache-2.0; see NOTICE), each
//! corrected to the specification where it differs.
//!
//! Audio is 16 kHz, 16-bit mono PCM in frames of 320 samples (20 ms), 14-bit as the reference
//! codec uses it. The [`Decoder`] turns one frame of serial bits (TS 26.173 parameter order,
//! one bit per `i16`) into PCM, the [`Encoder`] PCM into such bits; [`rfc4867`] reads such
//! frames from bytes and writes them.
//!
//! **One deliberate difference from TS 26.173:** NO_DATA frames are handled as before CR
//! 26.173-0027..29 (2007), as opencore does: a NO_DATA frame during speech starts the
//! comfort-noise hangover at once, where the current specification first waits for 8 frames
//! without voice activity. The two agree on every conformance vector; they differ only when
//! NO_DATA frames arrive in the middle of speech (frame loss signalled as NO_DATA).
//!
//! ```
//! use amr_wb::{Decoder, Dtx, Encoder, L_FRAME, Mode, rfc4867};
//!
//! // Encode one second of a 440 Hz tone at 12.65 kbit/s and store it as an `.awb` file.
//! let pcm: Vec<i16> = (0..16000)
//!     .map(|n| (8000.0 * (n as f64 * 440.0 * std::f64::consts::TAU / 16000.0).sin()) as i16)
//!     .collect();
//! let mut enc = Encoder::new(Dtx::Off);
//! let frames: Vec<rfc4867::Frame> = pcm.as_chunks::<L_FRAME>().0.iter()
//!     .map(|frame| enc.encode(Mode::MD1265, frame))
//!     .collect();
//! let file = rfc4867::write_storage(&frames)?;
//!
//! // Read it back and decode.
//! let mut dec = Decoder::new();
//! let mut out = Vec::new();
//! for frame in rfc4867::read_storage(&file)? {
//!     let f = frame?;
//!     let mut samples = [0i16; L_FRAME];
//!     dec.decode(f.frame_type, f.mode, &f.bits, &mut samples);
//!     out.extend_from_slice(&samples);
//! }
//! assert_eq!(out.len(), pcm.len());
//! # Ok::<(), rfc4867::Error>(())
//! ```
#![deny(missing_docs)]
mod basic_op;
mod bits;
mod cnst;
mod decoder;
mod dtx;
mod enc_tables;
mod encoder;
mod excitation;
mod gains;
mod lpc;
pub mod rfc4867;
mod simd;
// Stage timers: public only in `--cfg codec_stages` builds (`BENCH=stages just bench`).
#[cfg(codec_stages)]
#[doc(hidden)]
pub mod stages;
#[cfg(not(codec_stages))]
#[allow(dead_code)] // the names and `take` are read only by the stages bench
mod stages;
mod synthesis;
mod tables;
#[cfg(test)]
mod testutil;
#[cfg(all(test, oracle))]
mod vs_ref;

pub use cnst::{L_FRAME, MAX_SERIAL_SIZE, Mode, RxFrameType, SID_BITS};
pub use decoder::Decoder;
pub use encoder::{Dtx, Encoder};

/// The README's example compiles and runs.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
struct ReadmeDoctest;
