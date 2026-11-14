//! Excitation decoding: algebraic codebooks, pitch lag decoding, long-term prediction
//! (`Pred_lt_3or6`) and anti-sparseness phase dispersion (`ph_disp`).
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `dec/{d2_9pf,d2_11pf,d3_14pf,d4_17pf,
//! d8_31pf,d1035pf,dec_lag3,dec_lag6,ph_disp}.cpp` and `common/pred_lt.cpp`.
#![allow(dead_code)] // wired into the decoder later
#![allow(unused_imports)] // re-exports used once the decoder lands

mod codebook;
mod lag;
mod ph_disp;
pub(crate) mod pred_lt;

pub(crate) use codebook::{
    dec_8i40_31bits, dec_10i40_35bits, decode_2i40_9bits, decode_2i40_11bits, decode_3i40_14bits,
    decode_4i40_17bits,
};
pub(crate) use lag::{dec_lag3, dec_lag6};
pub(crate) use ph_disp::{
    PhDispState, ph_disp, ph_disp_lock, ph_disp_release, ph_disp_reset, total_exc,
};
pub(crate) use pred_lt::pred_lt_3or6;

pub(crate) use crate::cnst::{L_INTERPOL, L_SUBFR, Mode, PIT_MAX};
