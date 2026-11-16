//! Decoder stage timers for time attribution (callgrind counts instructions, `perf` is not
//! available): `--cfg codec_stages` (`BENCH=stages CRATE=amr-wb just bench`) times the stretch
//! since the previous mark and adds it to the stage being marked; otherwise every call is a no-op.

/// Stage indices (`mark` after each stage).
pub const PARAMS: usize = 0;
pub const EXC: usize = 1;
pub const LF: usize = 2;
pub const OVERSAMP: usize = 3;
pub const HF_GAIN: usize = 4;
pub const HF_SYN: usize = 5;
pub const HF_FIR: usize = 6;
pub const TAIL: usize = 7;
/// Their names, as `BENCH=stages` prints them.
pub const NAMES: [&str; 9] = [
    "isf+dtx", "exc", "lf_synth", "oversamp", "hf_gain", "hf_syn", "hf_fir", "tail", "-",
];

#[cfg(codec_stages)]
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
