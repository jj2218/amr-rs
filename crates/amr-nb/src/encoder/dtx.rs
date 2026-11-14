//! Encoder-side DTX: `dtx_buffer`, `tx_dtx_handler`, `dtx_enc` (SID parameters) and `sid_sync`
//! (TX frame type sequencing).
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `enc/src/{dtx_enc,sid_sync}.cpp`.

use super::lpc::{QPlsfState, q_plsf_3};
use crate::basic_op::{Flag, MAX_16, MAX_32, MIN_16, add, l_add, neg_shr, sub};
use crate::cnst::{L_FRAME, M, Mode};
use crate::gains::GcPredState;
use crate::lpc::{lsf_lsp, lsp_lsf, reorder_lsf};
use crate::math::log2;
use crate::tables::LSP_INIT_DATA;

const DTX_HIST_SIZE: usize = 8;
const DTX_HANG_CONST: i16 = 7;
const DTX_ELAPSED_FRAMES_THRESH: i16 = 24 + 7 - 1;
const LSF_GAP: i16 = 205;

/// `dtx_encState`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct DtxEncState {
    pub lsp_hist: [i16; M * DTX_HIST_SIZE],
    pub log_en_hist: [i16; DTX_HIST_SIZE],
    pub hist_ptr: i16,
    pub log_en_index: i16,
    pub init_lsf_vq_index: i16,
    pub lsp_index: [i16; 3],
    pub dtx_hangover_count: i16,
    pub dec_ana_elapsed_count: i16,
}

impl Default for DtxEncState {
    /// `dtx_enc_reset`.
    fn default() -> Self {
        let mut lsp_hist = [0; M * DTX_HIST_SIZE];
        for h in lsp_hist.as_chunks_mut::<M>().0 {
            h.copy_from_slice(&LSP_INIT_DATA);
        }
        Self {
            lsp_hist,
            log_en_hist: [0; DTX_HIST_SIZE],
            hist_ptr: 0,
            log_en_index: 0,
            init_lsf_vq_index: 0,
            lsp_index: [0; 3],
            dtx_hangover_count: DTX_HANG_CONST,
            dec_ana_elapsed_count: 32767,
        }
    }
}

/// `dtx_enc`: computes (when due) and writes the 5 SID parameters to `ana`.
pub(crate) fn dtx_enc(
    st: &mut DtxEncState,
    compute_sid_flag: bool,
    q_st: &mut QPlsfState,
    pred_state: &mut GcPredState,
    ana: &mut [i16],
    ov: &mut Flag,
) {
    // TS 26.073 (current): only when asked. (opencore, from an older version, also recomputes
    // while `log_en_index == 0`, i.e. after digital silence; the TS 26.074 vectors pass both.)
    if compute_sid_flag {
        let mut log_en: i16 = 0;
        let mut l_lsp = [0i32; M];
        for i in (0..DTX_HIST_SIZE).rev() {
            let temp = neg_shr(st.log_en_hist[i] as i32, 2) as i16;
            log_en = add(log_en, temp, ov);
            for j in (0..M).rev() {
                l_lsp[j] = l_add(l_lsp[j], st.lsp_hist[i * M + j] as i32, ov);
            }
        }
        log_en = neg_shr(log_en as i32, 1) as i16;
        let mut lsp = [0i16; M];
        for j in (0..M).rev() {
            lsp[j] = neg_shr(l_lsp[j], 3) as i16;
        }
        st.log_en_index = log_en.wrapping_add(2560);
        st.log_en_index = st.log_en_index.wrapping_add(128);
        st.log_en_index = neg_shr(st.log_en_index as i32, 8) as i16;
        st.log_en_index = st.log_en_index.clamp(0, 63);
        log_en = ((st.log_en_index as i32) << 8) as i16;
        log_en = sub(log_en, 11560, ov);
        log_en = log_en.clamp(-14436, 0);
        pred_state.past_qua_en = [log_en; 4];
        let log_en = ((5443 * log_en as i32) >> 15) as i16;
        pred_state.past_qua_en_mr122 = [log_en; 4];
        let mut lsf = [0i16; M];
        lsp_lsf(&lsp, &mut lsf, M, ov);
        reorder_lsf(&mut lsf, LSF_GAP, M, ov);
        lsf_lsp(&lsf, &mut lsp, M, ov);
        let mut lsp_q = [0i16; M];
        q_plsf_3(
            q_st,
            Mode::MRDTX,
            &lsp,
            &mut lsp_q,
            &mut st.lsp_index,
            &mut st.init_lsf_vq_index,
            ov,
        );
    }
    ana[0] = st.init_lsf_vq_index;
    ana[1..4].copy_from_slice(&st.lsp_index);
    ana[4] = st.log_en_index;
}

/// `dtx_buffer`: stores the frame's LSPs and log energy in the DTX history.
pub(crate) fn dtx_buffer(st: &mut DtxEncState, lsp_new: &[i16; M], speech: &[i16], ov: &mut Flag) {
    st.hist_ptr += 1;
    if st.hist_ptr as usize == DTX_HIST_SIZE {
        st.hist_ptr = 0;
    }
    let p = st.hist_ptr as usize * M;
    st.lsp_hist[p..p + M].copy_from_slice(lsp_new);
    let mut l_frame_en: i32 = 0;
    for &x in &speech[..L_FRAME] {
        l_frame_en = l_frame_en.wrapping_add((x as i32 * x as i32) << 1);
        if l_frame_en < 0 {
            l_frame_en = MAX_32;
            break;
        }
    }
    let (log_en_e, log_en_m) = log2(l_frame_en, ov);
    let l_temp = (log_en_e as i32) << 10;
    let mut log_en = if l_temp != l_temp as i16 as i32 {
        *ov = true;
        if log_en_e > 0 { MAX_16 } else { MIN_16 }
    } else {
        l_temp as i16
    };
    log_en = log_en.wrapping_add(log_en_m >> 5);
    log_en = log_en.wrapping_sub(8521);
    st.log_en_hist[st.hist_ptr as usize] = log_en >> 1;
}

/// `tx_dtx_handler`: DTX hangover; switches `used_mode` to `MRDTX`. Returns whether a new SID
/// may be computed.
pub(crate) fn tx_dtx_handler(
    st: &mut DtxEncState,
    vad_flag: bool,
    used_mode: &mut Mode,
    ov: &mut Flag,
) -> bool {
    st.dec_ana_elapsed_count = add(st.dec_ana_elapsed_count, 1, ov);
    if vad_flag {
        st.dtx_hangover_count = DTX_HANG_CONST;
        false
    } else if st.dtx_hangover_count == 0 {
        st.dec_ana_elapsed_count = 0;
        *used_mode = Mode::MRDTX;
        true
    } else {
        st.dtx_hangover_count -= 1;
        if add(st.dec_ana_elapsed_count, st.dtx_hangover_count, ov) < DTX_ELAPSED_FRAMES_THRESH {
            *used_mode = Mode::MRDTX;
        }
        false
    }
}

/// TX frame type (`enum TXFrameType`, the subset `sid_sync` produces).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum TxFrameType {
    SpeechGood,
    SidFirst,
    SidUpdate,
    NoData,
}

/// `sid_syncState`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct SidSyncState {
    pub sid_update_rate: i16,
    pub sid_update_counter: i16,
    pub sid_handover_debt: i16,
    pub prev_ft: TxFrameType,
}

impl Default for SidSyncState {
    fn default() -> Self {
        Self {
            sid_update_rate: 8,
            sid_update_counter: 3,
            sid_handover_debt: 0,
            prev_ft: TxFrameType::SpeechGood,
        }
    }
}

impl SidSyncState {
    /// `sid_sync_reset` (keeps the update rate).
    pub(crate) fn reset(&mut self) {
        self.sid_update_counter = 3;
        self.sid_handover_debt = 0;
        self.prev_ft = TxFrameType::SpeechGood;
    }
}

/// `sid_sync`: TX frame type of a frame encoded with `mode` (the used mode).
pub(crate) fn sid_sync(st: &mut SidSyncState, mode: Mode) -> TxFrameType {
    let ft = if mode == Mode::MRDTX {
        st.sid_update_counter -= 1;
        if st.prev_ft == TxFrameType::SpeechGood {
            st.sid_update_counter = 3;
            TxFrameType::SidFirst
        } else if st.sid_handover_debt > 0 && st.sid_update_counter > 2 {
            st.sid_handover_debt -= 1;
            TxFrameType::SidUpdate
        } else if st.sid_update_counter == 0 {
            st.sid_update_counter = st.sid_update_rate;
            TxFrameType::SidUpdate
        } else {
            TxFrameType::NoData
        }
    } else {
        st.sid_update_counter = st.sid_update_rate;
        TxFrameType::SpeechGood
    };
    st.prev_ft = ft;
    ft
}
