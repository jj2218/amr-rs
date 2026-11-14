//! Encoder and decoder stage timers for time attribution (callgrind counts instructions, `perf` is not
//! available): `--cfg amr_nb_stages` (`BENCH=stages just bench`) times the stretch since the
//! previous mark and adds it to the stage being marked; otherwise every call is a no-op.

/// Stage indices, in frame order (`mark` after each stage).
pub const PRE: usize = 0;
pub const LPC: usize = 1;
pub const OL: usize = 2;
pub const SUB_PRE: usize = 3;
pub const CL: usize = 4;
pub const CB: usize = 5;
pub const GAIN: usize = 6;
pub const SUB_POST: usize = 7;
pub const TAIL: usize = 8;
/// The decoder's stages (the same slots, `NAMES_DEC`).
pub const D_LSP: usize = 0;
pub const D_PARAMS: usize = 1;
pub const D_EXC: usize = 2;
pub const D_SYNTH: usize = 3;
pub const D_TAIL: usize = 4;
pub const D_POSTFILTER: usize = 5;
pub const D_POSTPROC: usize = 6;
/// The decoder's stage names.
pub const NAMES_DEC: [&str; 9] = [
    "bits+lsp",
    "pitch+cb+gain",
    "ph_disp",
    "synth+exc",
    "bgn+tail",
    "postfilter",
    "post_proc",
    "-",
    "-",
];
/// Their names, as `BENCH=stages` prints them.
pub const NAMES: [&str; 9] = [
    "pre+vad", "lpc+lsp", "ol_pitch", "sub_pre", "cl_pitch", "codebook", "gains", "sub_post",
    "tail",
];

#[cfg(amr_nb_stages)]
mod imp {
    use std::cell::Cell;
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
    use std::time::Instant;

    static NS: [AtomicU64; 9] = [const { AtomicU64::new(0) }; 9];
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
    pub fn take() -> [u64; 9] {
        core::array::from_fn(|i| NS[i].swap(0, Relaxed))
    }
}

#[cfg(not(amr_nb_stages))]
mod imp {
    #[inline(always)]
    pub(crate) fn start() {}
    #[inline(always)]
    pub(crate) fn mark(_: usize) {}
}

#[cfg(amr_nb_stages)]
pub use imp::take;
pub(crate) use imp::{mark, start};
