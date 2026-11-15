//! Pure-Rust, bit-exact AMR-WB (3GPP TS 26.173, ITU-T G.722.2) speech decoder, with RFC 4867
//! packing (`.awb` files and RTP payloads).
//!
//! Bit-exact on the official 3GPP TS 26.174 decoder test vectors and differentially tested
//! against the TS 26.173 reference C. Ported from opencore-amr 0.1.6 (Apache-2.0); see NOTICE.
//!
//! Audio is 16 kHz, 16-bit mono PCM in frames of 320 samples (20 ms), 14-bit as the reference
//! decoder outputs it. The [`Decoder`] turns one frame of serial bits (TS 26.173 parameter
//! order, one bit per `i16`) into PCM; [`rfc4867`] reads such frames from bytes (and writes
//! them).
//!
//! **One deliberate difference from TS 26.173:** NO_DATA frames are handled as before CR
//! 26.173-0027..29 (2007), as opencore does: a NO_DATA frame during speech starts the
//! comfort-noise hangover at once, where the current specification first waits for 8 frames
//! without voice activity. The two agree on every conformance vector; they differ only when
//! NO_DATA frames arrive in the middle of speech (frame loss signalled as NO_DATA).
#![deny(missing_docs)]
mod basic_op;
mod bits;
mod cnst;
mod decoder;
mod dtx;
mod excitation;
mod gain;
mod lpc;
pub mod rfc4867;
mod synthesis;
mod tables;
#[cfg(test)]
mod testutil;
#[cfg(test)]
mod vs_c;

pub use cnst::{L_FRAME16K, MAX_SERIAL_SIZE, Mode, RxFrameType, SID_BITS};
pub use decoder::Decoder;
