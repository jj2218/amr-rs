//! The 3GPP TS 26.073 AMR-NB reference C code (VAD option 1), the unit-test oracle: every
//! symbol it defines is prefixed `ref_` (see build.rs). Tests declare the functions they call
//! themselves, with the reference signatures (`Word16 = i16`, `Word32 = i32`, `Flag = i32`).
//!
//! The reference keeps its overflow flag in one global (`ref_Overflow`), shared by every
//! thread: hold [`lock`] while calling into it (tests run in parallel).
#![allow(unsafe_code, non_upper_case_globals)]

use std::sync::{Mutex, MutexGuard};

unsafe extern "C" {
    /// The reference's global overflow flag (`Flag Overflow` in basicop2.c).
    pub static mut ref_Overflow: i32;
}

static LOCK: Mutex<()> = Mutex::new(());

/// Serialises use of the reference code (its global flag); a panicking holder doesn't poison it.
pub fn lock() -> MutexGuard<'static, ()> {
    LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

/// Reads the global overflow flag (hold [`lock`]).
pub fn overflow() -> bool {
    // SAFETY: a plain int the C code reads and writes; callers hold `lock`.
    unsafe { ref_Overflow != 0 }
}

/// Clears the global overflow flag (hold [`lock`]).
pub fn clear_overflow() {
    // SAFETY: as `overflow`.
    unsafe { ref_Overflow = 0 }
}
