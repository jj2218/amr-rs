//! FFI to opencore-amr, used only as a test and benchmark oracle.
//!
//! Public interfaces (`interf_dec.h`/`interf_enc.h`, `dec_if.h`) plus internal functions that
//! the Rust ports are differentially tested against. Signatures mirror the C headers exactly.
//! Features: `nb` (AMR-NB codec), `wb` (AMR-WB decoder).
#![allow(unsafe_code, missing_docs, non_snake_case)]

use std::ffi::c_void;

#[cfg(feature = "nb")]
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

#[cfg(feature = "wb")]
unsafe extern "C" {
    /// RFC 4867 storage-format frames (ToC byte first), as `amrwb/wrapper.cpp`.
    pub fn D_IF_init() -> *mut c_void;
    pub fn D_IF_exit(state: *mut c_void);
    pub fn D_IF_decode(state: *mut c_void, input: *const u8, out: *mut i16, bfi: i32);

    /// `shim_wb.c`: the same frame loop on the serial parameter array (one `i16` per bit,
    /// ±127), as the TS 26.174 vectors store it. `frame_type` is `RX_*` (`dtx.h`).
    pub fn oc_wb_init() -> *mut c_void;
    pub fn oc_wb_exit(state: *mut c_void);
    pub fn oc_wb_decode(
        state: *mut c_void,
        mode: i16,
        frame_type: i16,
        prms: *mut i16,
        out: *mut i16,
    );
}
