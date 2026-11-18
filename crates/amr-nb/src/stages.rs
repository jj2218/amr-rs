//! Decoder and encoder stage timers for time attribution (callgrind counts instructions, `perf`
//! is not available): `--cfg codec_stages` (`BENCH=stages just bench`) times the stretch since
//! the previous mark and adds it to the stage being marked; otherwise every call is a no-op.
//! Decoder stages are `D_*` (`NAMES_DEC`), encoder stages `E_*` (`NAMES_ENC`), in frame order.

/// Stage slots.
pub const N: usize = 9;
/// Decoder stage indices (`mark` after each stage).
pub const D_LSP: usize = 0;
pub const D_PARAMS: usize = 1;
pub const D_EXC: usize = 2;
pub const D_SYNTH: usize = 3;
pub const D_TAIL: usize = 4;
pub const D_POSTFILTER: usize = 5;
pub const D_POSTPROC: usize = 6;
/// The decoder's stage names.
pub const NAMES_DEC: [&str; 7] = [
    "bits+lsp",
    "pitch+cb+gain",
    "ph_disp",
    "synth+exc",
    "bgn+tail",
    "postfilter",
    "post_proc",
];
/// Encoder stage indices (`mark` after each stage).
pub const E_PRE: usize = 0;
pub const E_LPC: usize = 1;
pub const E_OL: usize = 2;
pub const E_SUB_PRE: usize = 3;
pub const E_CL: usize = 4;
pub const E_CB: usize = 5;
pub const E_GAIN: usize = 6;
pub const E_SUB_POST: usize = 7;
pub const E_TAIL: usize = 8;
/// The encoder's stage names.
pub const NAMES_ENC: [&str; 9] = [
    "pre+vad", "lpc+lsp", "ol_pitch", "sub_pre", "cl_pitch", "codebook", "gains", "sub_post",
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
