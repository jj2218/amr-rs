//! Shared codec constants (`cnst.h`) and the codec mode (`mode.h`).

/// LPC order.
pub(crate) const M: usize = 10;
/// LPC order + 1.
pub(crate) const MP1: usize = M + 1;
/// Samples per 20 ms frame (8 kHz).
pub const L_FRAME: usize = 160;
/// Samples per subframe.
pub(crate) const L_SUBFR: usize = 40;
/// Maximum pitch lag.
pub(crate) const PIT_MAX: usize = 143;
/// Length of the pitch interpolation filter + 1.
pub(crate) const L_INTERPOL: usize = 11;

/// AMR-NB codec mode (bit rate), TS 26.073 numbering (0..=7).
#[allow(clippy::upper_case_acronyms)] // the standard's names
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Mode {
    /// 4.75 kbit/s
    MR475 = 0,
    /// 5.15 kbit/s
    MR515,
    /// 5.90 kbit/s
    MR59,
    /// 6.70 kbit/s
    MR67,
    /// 7.40 kbit/s
    MR74,
    /// 7.95 kbit/s
    MR795,
    /// 10.2 kbit/s
    MR102,
    /// 12.2 kbit/s
    MR122,
}

impl Mode {
    /// All modes, slowest first.
    pub const ALL: [Mode; 8] = [
        Mode::MR475,
        Mode::MR515,
        Mode::MR59,
        Mode::MR67,
        Mode::MR74,
        Mode::MR795,
        Mode::MR102,
        Mode::MR122,
    ];

    /// The mode with TS 26.073 number `n` (0..=7).
    pub fn from_index(n: usize) -> Option<Mode> {
        Self::ALL.get(n).copied()
    }

    /// Serial bits in a speech frame of this mode.
    pub fn bits(self) -> usize {
        crate::tables::NUMOFBITS[self as usize] as usize
    }
}

/// Serial bits in a SID frame: its comfort-noise parameters (TS 26.101 table 3; the STI and
/// the mode indication that follow them in a packed frame are not serial bits).
pub const SID_BITS: usize = 35;

/// The codec's internal mode, TS 26.073 `enum Mode`: the speech modes and `MRDTX`, the used
/// mode of comfort-noise (SID) frames. Ordered as in C (the codec compares modes, e.g.
/// `mode <= MR74`).
#[allow(clippy::upper_case_acronyms)] // the standard's names
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) enum CodecMode {
    MR475 = 0,
    MR515,
    MR59,
    MR67,
    MR74,
    MR795,
    MR102,
    MR122,
    MRDTX,
}

impl CodecMode {
    /// All modes in discriminant order: `ALL[i]` is mode `i` (speech modes are `ALL[..8]`).
    pub(crate) const ALL: [CodecMode; 9] = [
        CodecMode::MR475,
        CodecMode::MR515,
        CodecMode::MR59,
        CodecMode::MR67,
        CodecMode::MR74,
        CodecMode::MR795,
        CodecMode::MR102,
        CodecMode::MR122,
        CodecMode::MRDTX,
    ];
}

impl From<Mode> for CodecMode {
    fn from(m: Mode) -> Self {
        CodecMode::ALL[m as usize]
    }
}
