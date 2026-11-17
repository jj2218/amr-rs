//! The AMR-WB encoder (3GPP TS 26.190), ported from vo-amrwbenc 0.1.3 (Apache-2.0) with the
//! TS 26.173 reference's semantics (vo's `vo_*` macros wrap where the reference saturates),
//! including CR 26.173-0032 (pitch-gain clipping at 6.60/8.85) and -0036 (`c4t64fx`
//! saturation), which vo-amrwbenc predates.
// Faithful ports: index loops and the C argument lists.
#![allow(clippy::needless_range_loop, clippy::too_many_arguments)]

mod acelp2;
mod acelp4;
mod bits;
mod coder;
mod dtx;
mod filters;
#[cfg(test)]
mod kern;
mod lpc;
mod ops;
mod pitch;
mod pitch_ol;
mod preprocess;
mod q_gain;
mod vad;

/// The portable bodies of the dispatched kernels (`crate::simd`).
pub(crate) mod kernels {
    pub(crate) use super::acelp2::{acelp_2t64_fx_portable, cor_h_x_portable};
    pub(crate) use super::acelp4::acelp_4t64_fx_portable;
    pub(crate) use super::filters::{convolve_portable, residu_portable, syn_filt_portable};
    pub(crate) use super::lpc::{
        autocorr_portable, az_isp_portable, qpisf_2s_36b_portable, qpisf_2s_46b_portable,
    };
    pub(crate) use super::pitch::pitch_fr4_portable;
    pub(crate) use super::pitch_ol::pitch_med_ol_portable;
    pub(crate) use super::preprocess::{NB_COEF_DOWN, decim_12k8_portable};
}

use crate::cnst::{EHF_MASK, L_FRAME16K, MAX_SERIAL_SIZE, Mode, RxFrameType};
use crate::rfc4867::Frame;
use bits::BIT_1;
use coder::CoderState;
use dtx::MRDTX;

/// Transmitted frame type (TS 26.173 `TX_*`, `dtx.h`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum TxFrameType {
    /// A speech frame.
    Speech = 0,
    /// The first SID frame of a DTX period (no comfort-noise parameters used).
    SidFirst,
    /// A SID update: comfort-noise parameters ([`crate::SID_BITS`] bits).
    SidUpdate,
    /// Nothing to send (DTX pause).
    NoData,
}

/// What a decoder receives for a frame sent without errors.
impl From<TxFrameType> for RxFrameType {
    fn from(t: TxFrameType) -> Self {
        match t {
            TxFrameType::Speech => RxFrameType::SpeechGood,
            TxFrameType::SidFirst => RxFrameType::SidFirst,
            TxFrameType::SidUpdate => RxFrameType::SidUpdate,
            TxFrameType::NoData => RxFrameType::NoData,
        }
    }
}

/// AMR-WB encoder: 16 kHz, 16-bit PCM frames of 320 samples to serial bits, as the TS 26.173
/// reference encoder's frame loop (`coder.c`): the input is taken at 14 bits (the 2 LSBs
/// cleared), an encoder homing frame resets the encoder after it is coded, and with DTX the
/// SID frames are scheduled as its `Write_serial` does. Output frames use the representation
/// of [`crate::rfc4867::Frame`], ready for packing or decoding.
///
/// ```
/// use amr_wb::{Decoder, Encoder, L_FRAME16K, Mode, RxFrameType};
///
/// let mut enc = Encoder::new(true); // with DTX
/// let mut dec = Decoder::new();
/// let tone: [i16; L_FRAME16K] = std::array::from_fn(|n| ((n as f64 * 0.3).sin() * 8000.0) as i16);
/// let f = enc.encode(Mode::MD1265, &tone);
/// assert_eq!(f.frame_type, RxFrameType::SpeechGood);
/// let mut out = [0i16; L_FRAME16K];
/// dec.decode(f.frame_type, f.mode, &f.bits, &mut out);
/// // silence: after the VAD's hangover, comfort noise (SID) and NO_DATA frames
/// let types: Vec<_> = (0..30).map(|_| enc.encode(Mode::MD1265, &[0; L_FRAME16K]).frame_type).collect();
/// assert!(types.contains(&RxFrameType::SidFirst) && types[29] == RxFrameType::NoData);
/// ```
#[derive(Clone)]
pub struct Encoder {
    st: Box<CoderState>,
    dtx: bool,
    // `TX_State` (bits.c), without its `sid_handover_debt`: never raised in the reference
    sid_update_counter: i16,
    prev_ft: TxFrameType,
}

impl Encoder {
    /// An encoder in its initial (homed) state; `dtx` enables voice activity detection and
    /// discontinuous transmission (comfort noise).
    pub fn new(dtx: bool) -> Self {
        Self {
            st: Box::new(CoderState::new()),
            dtx,
            sid_update_counter: 3,
            prev_ft: TxFrameType::Speech,
        }
    }

    /// Back to the initial state (the DTX setting is kept).
    pub fn reset(&mut self) {
        self.st.reset(true);
        self.sid_update_counter = 3;
        self.prev_ft = TxFrameType::Speech;
    }

    /// Encodes one frame at `mode`. The frame's serial bits, one per element (0 or 1, zero past
    /// the frame's size) in TS 26.173 parameter order, as [`crate::Decoder::decode`] reads
    /// them: `mode.bits()` for a speech frame, [`crate::SID_BITS`] for a SID frame. With DTX
    /// the frame can be `SidFirst`, `SidUpdate` or `NoData` (`mode` is then the mode
    /// indication); as in the reference's output, `SidFirst` and `NoData` frames also carry
    /// the comfort-noise analysis in their first [`crate::SID_BITS`] bits, which are not
    /// transmitted (a SID_FIRST frame's parameters are unused: [`crate::rfc4867`] writes them as
    /// zeros; a NO_DATA frame is not sent).
    pub fn encode(&mut self, mode: Mode, speech: &[i16; L_FRAME16K]) -> Frame {
        crate::stages::start();
        let homing = speech.iter().all(|&s| s == EHF_MASK);
        let signal = speech.map(|s| s & !3);
        let mut prms = [0i16; MAX_SERIAL_SIZE];
        let (used, n) = self.st.coder(mode as i16, &signal, &mut prms, self.dtx);
        if homing {
            self.st.reset(true);
        }
        let mut bits = [0i16; MAX_SERIAL_SIZE];
        for (b, &p) in bits.iter_mut().zip(&prms[..n]) {
            *b = (p == BIT_1) as i16;
        }
        let frame = Frame {
            frame_type: self.frame_type(used == MRDTX).into(),
            mode,
            bits,
        };
        crate::stages::mark(crate::stages::E_TAIL);
        frame
    }

    /// `Write_serial`'s TX frame type for a SID (`MRDTX`) or speech frame.
    fn frame_type(&mut self, sid: bool) -> TxFrameType {
        use TxFrameType::*;
        let ft = if sid {
            self.sid_update_counter -= 1;
            if self.prev_ft == Speech {
                self.sid_update_counter = 3;
                SidFirst
            } else if self.sid_update_counter == 0 {
                self.sid_update_counter = 8;
                SidUpdate
            } else {
                NoData
            }
        } else {
            self.sid_update_counter = 8;
            Speech
        };
        self.prev_ft = ft;
        ft
    }
}
