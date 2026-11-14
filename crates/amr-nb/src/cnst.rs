//! Shared codec constants (`cnst.h`) and the codec mode (`mode.h`).

/// LPC order.
pub(crate) const M: usize = 10;
/// LPC order + 1.
pub(crate) const MP1: usize = M + 1;
/// Samples per frame.
pub(crate) const L_FRAME: usize = 160;
/// Samples per subframe.
pub(crate) const L_SUBFR: usize = 40;
/// Maximum pitch lag.
pub(crate) const PIT_MAX: usize = 143;
/// Length of the pitch interpolation filter + 1.
pub(crate) const L_INTERPOL: usize = 11;

/// AMR-NB codec mode (bit rate), TS 26.073 `enum Mode`. Ordered as in C (the codec compares
/// modes, e.g. `mode <= MR74`).
#[allow(clippy::upper_case_acronyms)] // the standard's names
#[repr(i32)]
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
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
    /// Comfort-noise (SID) parameters.
    MRDTX,
}

impl Mode {
    /// All modes in discriminant order: `ALL[i]` is mode `i` (speech modes are `ALL[..8]`).
    pub(crate) const ALL: [Mode; 9] = [
        Mode::MR475,
        Mode::MR515,
        Mode::MR59,
        Mode::MR67,
        Mode::MR74,
        Mode::MR795,
        Mode::MR102,
        Mode::MR122,
        Mode::MRDTX,
    ];
}

impl std::fmt::Display for Mode {
    /// The C enumerator name, e.g. `MR122`.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Debug::fmt(self, f)
    }
}
