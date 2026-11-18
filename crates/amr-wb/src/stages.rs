//! Decoder and encoder stage timers for time attribution (callgrind counts instructions, `perf`
//! is not available): `--cfg codec_stages` (`CRATE=amr-wb BENCH=stages just bench`) times the
//! stretch since the previous mark and adds it to the stage being marked; otherwise every call is a no-op.
//! Decoder stages are `D_*` (`NAMES_DEC`), encoder stages `E_*` (`NAMES_ENC`), in frame order.

/// Stage slots.
pub const N: usize = 11;
/// Decoder stage indices (`mark` after each stage).
pub const D_PARAMS: usize = 0;
pub const D_EXC: usize = 1;
pub const D_LF: usize = 2;
pub const D_OVERSAMP: usize = 3;
pub const D_HF_GAIN: usize = 4;
pub const D_HF_SYN: usize = 5;
pub const D_HF_FIR: usize = 6;
pub const D_TAIL: usize = 7;
/// The decoder's stage names.
pub const NAMES_DEC: [&str; 8] = [
    "isf+dtx", "exc", "lf_synth", "oversamp", "hf_gain", "hf_syn", "hf_fir", "tail",
];
/// Encoder stage indices (`mark` after each stage).
pub const E_PRE: usize = 0;
pub const E_VAD: usize = 1;
pub const E_LPC: usize = 2;
pub const E_OL: usize = 3;
pub const E_ISFQ: usize = 4;
pub const E_SUB_PRE: usize = 5;
pub const E_CL: usize = 6;
pub const E_CB: usize = 7;
pub const E_GAIN: usize = 8;
pub const E_SUB_POST: usize = 9;
pub const E_TAIL: usize = 10;
/// The encoder's stage names.
pub const NAMES_ENC: [&str; 11] = [
    "pre",
    "vad+dtx",
    "lpc",
    "wsp+ol",
    "isf_q+res",
    "target+h",
    "cl_pitch",
    "codebook",
    "gains",
    "update+hf",
    "tail",
];

#[cfg(codec_stages)]
mod imp {
    use std::cell::Cell;
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
    use std::time::Instant;

    static NS: [AtomicU64; super::N] = [const { AtomicU64::new(0) }; super::N];
    thread_local! { static LAST: Cell<Option<Instant>> = const { Cell::new(None) }; }

    pub(crate) fn start() {
        LAST.set(Some(Instant::now()));
    }

    pub(crate) fn mark(stage: usize) {
        let now = Instant::now();
        if let Some(t) = LAST.get() {
            NS[stage].fetch_add((now - t).as_nanos() as u64, Relaxed);
        }
        LAST.set(Some(now));
    }

    /// Nanoseconds per stage since the last `take`.
    pub fn take() -> [u64; super::N] {
        core::array::from_fn(|i| NS[i].swap(0, Relaxed))
    }
}

#[cfg(not(codec_stages))]
mod imp {
    #[inline(always)]
    pub(crate) fn start() {}
    #[inline(always)]
    pub(crate) fn mark(_: usize) {}
}

#[cfg(codec_stages)]
pub use imp::take;
pub(crate) use imp::{mark, start};
