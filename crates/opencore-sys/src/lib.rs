//! FFI to opencore-amr (AMR-NB), used only as a test and benchmark oracle.
//!
//! Public interface (`interf_dec.h`/`interf_enc.h`) plus internal functions that the Rust
//! port is differentially tested against. Signatures mirror the C headers exactly.
#![allow(unsafe_code, missing_docs, non_snake_case)]

use std::ffi::c_void;

unsafe extern "C" {
    pub fn Decoder_Interface_init() -> *mut c_void;
    pub fn Decoder_Interface_exit(state: *mut c_void);
    pub fn Decoder_Interface_Decode(state: *mut c_void, input: *const u8, out: *mut i16, bfi: i32);

    pub fn Encoder_Interface_init(dtx: i32) -> *mut c_void;
    pub fn Encoder_Interface_exit(state: *mut c_void);
    /// NB: despite the `const` in C, this overwrites `speech` (pre-processing in place).
    pub fn Encoder_Interface_Encode(
        state: *mut c_void,
        mode: i32,
        speech: *const i16,
        out: *mut u8,
        force_speech: i32,
    ) -> i32;
}
