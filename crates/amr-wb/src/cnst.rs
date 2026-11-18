//! Codec constants (`pvamrwbdecoder_cnst.h`, `dtx.h`), the codec mode and the RX frame type.

/// Samples per 12.8 kHz frame.
pub(crate) const L_FRAME12K8: usize = 256;
/// Samples per 16 kHz subframe.
pub(crate) const L_SUBFR16K: usize = 80;
/// Samples per 12.8 kHz subframe.
pub(crate) const L_SUBFR: usize = 64;
pub(crate) const NB_SUBFR: usize = 4;
/// LPC order (12.8 kHz).
pub(crate) const M: usize = 16;
/// LPC order of the 16 kHz high band.
pub(crate) const M16K: usize = 20;
pub(crate) const PIT_SHARP: i16 = 27853;
pub(crate) const PIT_MIN: i16 = 34;
pub(crate) const PIT_FR2: i16 = 128;
pub(crate) const PIT_FR1_9B: i16 = 160;
pub(crate) const PIT_FR1_8B: i16 = 92;
pub(crate) const PIT_MAX: usize = 231;
pub(crate) const L_INTERPOL: usize = 16 + 1;
pub(crate) const PREEMPH_FAC: i16 = 22282;
pub(crate) const Q_MAX: i16 = 8;
pub(crate) const RANDOM_INITSEED: i16 = 21845;
pub(crate) const L_MEANBUF: usize = 3;
/// Encoder homing frame pattern, the decoder's output for a decoder homing frame.
pub(crate) const EHF_MASK: i16 = 0x0008;
pub(crate) const ISF_GAP: i16 = 128;

/// Bits per frame per mode (`AMR_WB_COMPRESSED`, `nb_of_bits`): modes 0..=8, then SID.
pub(crate) const NBBITS_7K: i16 = 132;
pub(crate) const NBBITS_9K: i16 = 177;
pub(crate) const NBBITS_12K: i16 = 253;
pub(crate) const NBBITS_14K: i16 = 285;
pub(crate) const NBBITS_16K: i16 = 317;
pub(crate) const NBBITS_18K: i16 = 365;
pub(crate) const NBBITS_20K: i16 = 397;
pub(crate) const NBBITS_23K: i16 = 461;
pub(crate) const NBBITS_24K: i16 = 477;

/// Samples per 20 ms output frame (16 kHz).
pub const L_FRAME: usize = 320;
/// The most serial bits a frame carries (23.85 kbit/s).
pub const MAX_SERIAL_SIZE: usize = 477;

/// AMR-WB codec mode (bit rate), TS 26.173 numbering (0..=8; 9 is SID).
#[allow(clippy::upper_case_acronyms)] // the standard's names
#[repr(i16)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Mode {
    /// 6.60 kbit/s
    MD66 = 0,
    /// 8.85 kbit/s
    MD885,
    /// 12.65 kbit/s
    MD1265,
    /// 14.25 kbit/s
    MD1425,
    /// 15.85 kbit/s
    MD1585,
    /// 18.25 kbit/s
    MD1825,
    /// 19.85 kbit/s
    MD1985,
    /// 23.05 kbit/s
    MD2305,
    /// 23.85 kbit/s
    MD2385,
}

impl Mode {
    /// All modes, slowest first.
    pub const ALL: [Mode; 9] = [
        Mode::MD66,
        Mode::MD885,
        Mode::MD1265,
        Mode::MD1425,
        Mode::MD1585,
        Mode::MD1825,
        Mode::MD1985,
        Mode::MD2305,
        Mode::MD2385,
    ];

    /// The mode with TS 26.173 number `n` (0..=8).
    pub fn from_index(n: usize) -> Option<Mode> {
        Self::ALL.get(n).copied()
    }

    /// Serial bits in a speech frame of this mode.
    pub fn bits(self) -> usize {
        crate::tables::AMR_WB_COMPRESSED[self as usize] as usize
    }
}

/// Serial bits in a SID frame.
pub const SID_BITS: usize = 35;

/// Received frame type (TS 26.173 `RX_*`, `dtx.h`), as given to the decoder.
#[repr(i16)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RxFrameType {
    /// Speech frame received correctly.
    SpeechGood = 0,
    /// Speech frame with possible bit errors (decoded, its gains concealed).
    SpeechProbablyDegraded,
    /// Speech frame lost (concealed; the mode is the previous one).
    SpeechLost,
    /// Speech frame with errors (concealed).
    SpeechBad,
    /// First SID frame of a DTX period.
    SidFirst,
    /// SID update with comfort-noise parameters.
    SidUpdate,
    /// Corrupted SID frame.
    SidBad,
    /// Nothing received (DTX pause; the mode is the previous one). See the crate docs for the
    /// one deliberate difference from TS 26.173 here.
    NoData,
}
