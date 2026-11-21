//! The AMR-NB encoder: `cod_amr` (`cod_amr.cpp`) with its subframe helpers (`pre_big`,
//! `subframePreProc`, `subframePostProc`), the input high-pass filter (`Pre_Process`), the frame
//! wrapper (`sp_enc.cpp`), `Prm2bits`, `sid_sync` and the TS 26.073 encoder homing, which
//! opencore lacks.
//! Bit-exact port of opencore-amr 0.1.6 (Apache-2.0); conformance-tested against the official
//! 3GPP TS 26.074 test vectors (`just conform`).

#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1

pub(crate) mod codebook;
mod dtx;
pub(crate) mod gain_q;
#[cfg(all(test, oracle))]
mod kern;
pub(crate) mod layout;
pub(crate) mod lpc;
pub(crate) mod pitch_cl;
pub(crate) mod pitch_ol;
mod vad;
mod vad2;

#[cfg(all(test, oracle))]
pub(crate) mod vs_ref;

use crate::basic_op::{Flag, add, extract_h, l_mult, l_shl, round, sub};
use crate::cnst::{CodecMode, L_FRAME, L_INTERPOL, L_SUBFR, M, MP1, Mode, PIT_MAX};
use crate::decoder::MAX_SERIAL_SIZE;
use crate::dtx::RxFrameType;
use crate::excitation::pred_lt_3or6;
use crate::post::filt::{residu, syn_filt, weight_ai};
use crate::rfc4867::Frame;
use crate::tables::{
    BITNO_MR59, BITNO_MR67, BITNO_MR74, BITNO_MR102, BITNO_MR122, BITNO_MR475, BITNO_MR515,
    BITNO_MR795, BITNO_MRDTX, PRMNO,
};
use codebook::cbsearch;
use dtx::{DtxEncState, SidSyncState, TxFrameType, dtx_buffer, dtx_enc, sid_sync, tx_dtx_handler};
use gain_q::{GainQuantState, Gains, gain_quant};
use lpc::{L_WINDOW, LevinsonState, LspState, lpc, lsp};
use pitch_cl::{PitchFrState, TonStabState, check_lsp, cl_ltp, convolve, update_gp_clipping};
use pitch_ol::{PitchOlWghtState, ol_ltp};
use vad::{LOOKAHEAD, Vad, VadState, vad_pitch_detection, vad1};
use vad2::{Vad2State, ltp_flag_update, vad2};

const L_TOTAL: usize = 320;
const L_NEXT: usize = 40;
const L_FRAME_BY2: usize = 80;
const MAX_PRM_SIZE: usize = 57;
const SHARPMAX: i16 = 13017;
const SHARPMIN: i16 = 0;
/// Offsets into `old_speech` of the C's `speech`, `new_speech`, `p_window`, `p_window_12k2`.
const SPEECH: usize = L_TOTAL - L_FRAME - L_NEXT;
const NEW_SPEECH: usize = L_TOTAL - L_FRAME;
const P_WINDOW: usize = L_TOTAL - L_WINDOW;
const P_WINDOW_12K2: usize = P_WINDOW - L_NEXT;
/// Offset into `old_exc` of the C's `exc`.
const EXC: usize = PIT_MAX + L_INTERPOL;
/// Encoder homing frame sample value (TS 26.073 `EHF_MASK`).
const EHF: i16 = 0x0008;

const GAMMA1: [i16; M] = [
    30802, 28954, 27217, 25584, 24049, 22606, 21250, 19975, 18777, 17650,
];
const GAMMA1_12K2: [i16; M] = [
    29491, 26542, 23888, 21499, 19349, 17414, 15672, 14105, 12694, 11425,
];
const GAMMA2: [i16; M] = [19661, 11797, 7078, 4247, 2548, 1529, 917, 550, 330, 198];

/// `Pre_ProcessState`: the input high-pass filter (cut-off 80 Hz) and /2 scaling.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct PreProcessState {
    y2_hi: i16,
    y2_lo: i16,
    y1_hi: i16,
    y1_lo: i16,
    x0: i16,
    x1: i16,
}

/// `Pre_Process` (in place).
///
/// The C keeps y1, y2 as `L_Extract(L)`, (hi, lo) with u = (hi << 15) + lo = L >> 1, where
/// L = sat(16 l). hi*a + ((lo*a) >> 15) is exactly (u*a) >> 15 (hi*a*2^15 is a multiple of
/// 2^15), so the serial chain has one multiply per tap instead of the split and two (|l| <
/// 2^30: no overflow).
fn pre_process(st: &mut PreProcessState, signal: &mut [i16]) {
    let mut u1 = ((st.y1_hi as i32) << 15) + st.y1_lo as i32;
    let mut u2 = ((st.y2_hi as i32) << 15) + st.y2_lo as i32;
    let mut x_n_2 = st.x1;
    let mut x_n_1 = st.x0;
    for s in signal.iter_mut() {
        let x = x_n_2 as i32 * 1899 + x_n_1 as i32 * -3798 + *s as i32 * 1899;
        (x_n_2, x_n_1) = (x_n_1, *s);
        let l = ((u1 as i64 * 7807) >> 15) as i32 + ((u2 as i64 * -3733) >> 15) as i32 + x;
        // The C's sum is 2l (|l| < 2^30: no saturation before) and `L_shl(.., 3)` saturates:
        // output round(16l), state 16l >> 1 (opencore wraps: wrong on clipped input).
        let lt = if l >= 1 << 27 {
            i32::MAX
        } else if l < -(1 << 27) {
            i32::MIN
        } else {
            l << 4
        };
        *s = round(lt, &mut false);
        (u2, u1) = (u1, lt >> 1);
    }
    ((st.y1_hi, st.y1_lo), (st.y2_hi, st.y2_lo)) = crate::post::post_pro::one_word_pairs(
        signal.len(),
        u1,
        u2,
        (st.y1_hi, st.y1_lo),
        (st.y2_hi, st.y2_lo),
    );
    st.x1 = x_n_2;
    st.x0 = x_n_1;
}

/// `cod_amrState` (sub-states inline; the C's buffer pointers are the offset constants above).
#[derive(Clone)]
struct CodAmr {
    old_speech: [i16; L_TOTAL],
    old_wsp: [i16; L_FRAME + PIT_MAX],
    old_lags: [i16; 5],
    ol_gain_flg: [i16; 2],
    old_exc: [i16; L_FRAME + PIT_MAX + L_INTERPOL],
    /// `hvec[L_SUBFR..]` (the C's `h1`); `hvec[..L_SUBFR]` is always zero.
    h1: [i16; L_SUBFR],
    lpc_st: LevinsonState,
    lsp_st: LspState,
    cl_ltp_st: PitchFrState,
    gain_quant_st: GainQuantState,
    pitch_ol_wght_st: PitchOlWghtState,
    ton_stab_st: TonStabState,
    vad_st: Vad,
    dtx: bool,
    dtx_mode: Dtx,
    dtx_enc_st: DtxEncState,
    mem_syn: [i16; M],
    mem_w0: [i16; M],
    mem_w: [i16; M],
    /// `mem_err[..M]` is the error filter memory, `mem_err[M..]` the C's `error`.
    mem_err: [i16; M + L_SUBFR],
    sharp: i16,
    overflow: Flag,
}

impl CodAmr {
    /// `cod_amr_init` + `cod_amr_reset`.
    fn new(dtx_mode: Dtx) -> Self {
        Self {
            old_speech: [0; L_TOTAL],
            old_wsp: [0; L_FRAME + PIT_MAX],
            old_lags: [40; 5],
            ol_gain_flg: [0; 2],
            old_exc: [0; L_FRAME + PIT_MAX + L_INTERPOL],
            h1: [0; L_SUBFR],
            lpc_st: LevinsonState::default(),
            lsp_st: LspState::default(),
            cl_ltp_st: PitchFrState::default(),
            gain_quant_st: GainQuantState::default(),
            pitch_ol_wght_st: PitchOlWghtState::default(),
            ton_stab_st: TonStabState::default(),
            vad_st: match dtx_mode {
                Dtx::Vad2 => Vad::Two(Vad2State::default()),
                _ => Vad::One(VadState::default()),
            },
            dtx: dtx_mode != Dtx::Off,
            dtx_mode,
            dtx_enc_st: DtxEncState::default(),
            mem_syn: [0; M],
            mem_w0: [0; M],
            mem_w: [0; M],
            mem_err: [0; M + L_SUBFR],
            sharp: SHARPMIN,
            overflow: false,
        }
    }

    /// `cod_amr_reset`. C clears only the history parts of some buffers (and not `h1` or
    /// `ol_gain_flg`); everything it leaves is rewritten before it is read, so a full reset is
    /// equivalent.
    fn reset(&mut self) {
        *self = Self::new(self.dtx_mode);
    }
}

/// `pre_big`: weighted speech `wsp` for half a frame (two subframes).
fn pre_big(
    mode: CodecMode,
    a_t: &[i16; 4 * MP1],
    frame_offset: usize,
    old_speech: &[i16; L_TOTAL],
    mem_w: &mut [i16; M],
    old_wsp: &mut [i16],
) {
    let g1 = if mode <= CodecMode::MR795 {
        &GAMMA1
    } else {
        &GAMMA1_12K2
    };
    let mut a_off = if frame_offset > 0 { 2 * MP1 } else { 0 };
    let mut off = frame_offset;
    for _ in 0..2 {
        let mut ap1 = [0i16; MP1];
        let mut ap2 = [0i16; MP1];
        weight_ai(&a_t[a_off..], g1, &mut ap1);
        weight_ai(&a_t[a_off..], &GAMMA2, &mut ap2);
        let wsp = &mut old_wsp[PIT_MAX + off..PIT_MAX + off + L_SUBFR];
        residu(&ap1, &old_speech[SPEECH + off - M..], wsp, L_SUBFR);
        let x: [i16; L_SUBFR] = (&*wsp).try_into().unwrap();
        syn_filt(&ap2, &x, wsp, L_SUBFR, mem_w, true);
        a_off += MP1;
        off += L_SUBFR;
    }
}

/// `subframePreProc`: impulse response `h1`, target `xn`, LP residual `res2` (also copied to
/// `exc`) for one subframe. `speech` starts `M` samples before the subframe.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
fn subframe_pre_proc(
    mode: CodecMode,
    a: &[i16],
    aq: &[i16],
    speech: &[i16],
    mem_err: &mut [i16; M + L_SUBFR],
    mem_w0: &[i16; M],
    exc: &mut [i16],
    h1: &mut [i16; L_SUBFR],
    xn: &mut [i16; L_SUBFR],
    res2: &mut [i16; L_SUBFR],
) {
    let g1 = if matches!(mode, CodecMode::MR122 | CodecMode::MR102) {
        &GAMMA1_12K2
    } else {
        &GAMMA1
    };
    let mut ap1 = [0i16; MP1];
    let mut ap2 = [0i16; MP1];
    weight_ai(a, g1, &mut ap1);
    weight_ai(a, &GAMMA2, &mut ap2);
    let mut ai_zero = [0i16; L_SUBFR];
    ai_zero[..MP1].copy_from_slice(&ap1);
    let mut zero = [0i16; M];
    syn_filt(aq, &ai_zero, h1, L_SUBFR, &mut zero, false);
    let h: [i16; L_SUBFR] = *h1;
    syn_filt(&ap2, &h, h1, L_SUBFR, &mut zero, false);
    residu(aq, speech, res2, L_SUBFR);
    exc[..L_SUBFR].copy_from_slice(res2);
    let (mem, error) = mem_err.split_at_mut(M);
    let mut mem_copy: [i16; M] = (&*mem).try_into().unwrap();
    syn_filt(aq, &exc[..L_SUBFR], error, L_SUBFR, &mut mem_copy, false);
    residu(&ap1, &mem_err[..], xn, L_SUBFR);
    let x: [i16; L_SUBFR] = *xn;
    let mut w0 = *mem_w0;
    syn_filt(&ap2, &x, xn, L_SUBFR, &mut w0, false);
}

/// `subframePostProc`: final excitation, synthesis, filter memory updates; updates `sharp`.
#[allow(clippy::too_many_arguments)] // mirrors the C signature
fn subframe_post_proc(
    speech: &[i16],
    mode: CodecMode,
    i_subfr: usize,
    gain_pit: i16,
    gain_code: i16,
    aq: &[i16],
    synth: &mut [i16; L_FRAME],
    xn: &[i16; L_SUBFR],
    code: &[i16; L_SUBFR],
    y1: &[i16; L_SUBFR],
    y2: &[i16; L_SUBFR],
    mem_syn: &mut [i16; M],
    mem_err: &mut [i16; M + L_SUBFR],
    mem_w0: &mut [i16; M],
    exc: &mut [i16],
    sharp: &mut i16,
) {
    let (temp_shift, k_shift, pitch_fac) = if mode != CodecMode::MR122 {
        (1, 2, gain_pit)
    } else {
        (2, 4, gain_pit >> 1)
    };
    *sharp = gain_pit.min(SHARPMAX);
    // Saturating like TS 26.073 `spstproc.c` (opencore wraps: wrong on full-scale input).
    let ov = &mut false;
    let e: &mut [i16; L_SUBFR] = (&mut exc[i_subfr..i_subfr + L_SUBFR]).try_into().unwrap();
    crate::excitation::total_exc::<true>(e, code, pitch_fac, gain_code, temp_shift);
    syn_filt(
        aq,
        &exc[i_subfr..i_subfr + L_SUBFR],
        &mut synth[i_subfr..i_subfr + L_SUBFR],
        L_SUBFR,
        mem_syn,
        true,
    );
    for (j, i) in (L_SUBFR - M..L_SUBFR).enumerate() {
        mem_err[j] = sub(speech[i_subfr + i], synth[i_subfr + i], ov);
        let temp = extract_h(l_shl(l_mult(y1[i], gain_pit, ov), 1, ov));
        let k = extract_h(l_shl(l_mult(y2[i], gain_code, ov), k_shift, ov));
        mem_w0[j] = sub(xn[i], add(temp, k, ov), ov);
    }
}

impl CodAmr {
    /// `cod_amr`: encodes one pre-processed frame into the parameters `ana`; returns the used
    /// mode (`MRDTX` when DTX decides on a SID/no-data frame).
    fn encode(
        &mut self,
        mode: CodecMode,
        new_speech: &[i16; L_FRAME],
        ana: &mut [i16; MAX_PRM_SIZE],
    ) -> CodecMode {
        let ov = &mut self.overflow;
        let mut a_t = [0i16; 4 * MP1];
        let mut aq_t = [0i16; 4 * MP1];
        let mut lsp_new = [0i16; M];
        let mut synth = [0i16; L_FRAME];
        let mut ana_pos = 0usize;
        self.old_speech[NEW_SPEECH..].copy_from_slice(new_speech);
        let mut used_mode = mode;
        let compute_sid_flag = if self.dtx {
            let vad_flag = match &mut self.vad_st {
                Vad::One(v) => vad1(v, &self.old_speech[NEW_SPEECH - LOOKAHEAD..], ov) != 0,
                Vad::Two(v) => {
                    let first = vad2(&self.old_speech[NEW_SPEECH..], v, ov);
                    vad2(&self.old_speech[NEW_SPEECH + L_FRAME_BY2..], v, ov) || first
                }
            };
            tx_dtx_handler(&mut self.dtx_enc_st, vad_flag, &mut used_mode, ov)
        } else {
            false
        };
        crate::stages::mark(crate::stages::E_PRE);
        lpc(
            &mut self.lpc_st,
            mode,
            &self.old_speech[P_WINDOW..],
            &self.old_speech[P_WINDOW_12K2..],
            &mut a_t,
            ov,
        );
        ana_pos += lsp(
            &mut self.lsp_st,
            mode,
            used_mode,
            &mut a_t,
            &mut aq_t,
            &mut lsp_new,
            &mut ana[..],
            ov,
        );
        dtx_buffer(
            &mut self.dtx_enc_st,
            &lsp_new,
            &self.old_speech[NEW_SPEECH..],
            ov,
        );
        let mut lsp_flag = false;
        if used_mode == CodecMode::MRDTX {
            dtx_enc(
                &mut self.dtx_enc_st,
                compute_sid_flag,
                &mut self.lsp_st.q_st,
                &mut self.gain_quant_st.gc_pred_st,
                &mut ana[ana_pos..],
                ov,
            );
            self.old_exc[..EXC].fill(0);
            self.mem_w0 = [0; M];
            self.mem_err[..M].fill(0);
            // (C also clears `zero` and `hvec[..L_SUBFR]`, which are always zero here.)
            self.lsp_st = LspState::default();
            self.lsp_st.lsp_old = lsp_new;
            self.lsp_st.lsp_old_q = lsp_new;
            self.cl_ltp_st = PitchFrState::default();
            self.sharp = SHARPMIN;
        } else {
            lsp_flag = check_lsp(&mut self.ton_stab_st, &self.lsp_st.lsp_old);
        }
        if self.dtx
            && let Vad::Two(v) = &mut self.vad_st
        {
            v.l_rmax = 0;
            v.l_r0 = 0;
        }

        crate::stages::mark(crate::stages::E_LPC);
        let mut t_op = [0i16; L_FRAME / L_FRAME_BY2];
        for sub_nr in 0..2 {
            let i_subfr = sub_nr * L_FRAME_BY2;
            pre_big(
                mode,
                &a_t,
                i_subfr,
                &self.old_speech,
                &mut self.mem_w,
                &mut self.old_wsp,
            );
            if mode != CodecMode::MR475 && mode != CodecMode::MR515 {
                t_op[sub_nr] = ol_ltp(
                    &mut self.pitch_ol_wght_st,
                    &mut self.vad_st,
                    mode,
                    &self.old_wsp[i_subfr..],
                    &mut self.old_lags,
                    &mut self.ol_gain_flg,
                    sub_nr,
                    self.dtx,
                    ov,
                );
            }
        }
        if mode == CodecMode::MR475 || mode == CodecMode::MR515 {
            t_op[0] = ol_ltp(
                &mut self.pitch_ol_wght_st,
                &mut self.vad_st,
                mode,
                &self.old_wsp,
                &mut self.old_lags,
                &mut self.ol_gain_flg,
                1,
                self.dtx,
                ov,
            );
            t_op[1] = t_op[0];
        }
        if self.dtx {
            match &mut self.vad_st {
                Vad::One(v) => vad_pitch_detection(v, &t_op, ov),
                Vad::Two(v) => ltp_flag_update(v, mode, ov),
            }
        }

        crate::stages::mark(crate::stages::E_OL);
        if used_mode != CodecMode::MRDTX {
            self.subframes(
                used_mode, &a_t, &aq_t, &t_op, lsp_flag, &mut synth, ana, ana_pos,
            );
            self.old_exc.copy_within(L_FRAME.., 0);
        }
        self.old_wsp.copy_within(L_FRAME.., 0);
        self.old_speech.copy_within(L_FRAME.., 0);
        used_mode
    }

    /// The subframe loop of `cod_amr` (speech frames only).
    #[allow(clippy::too_many_arguments)] // the C function's locals
    fn subframes(
        &mut self,
        mode: CodecMode,
        a_t: &[i16; 4 * MP1],
        aq_t: &[i16; 4 * MP1],
        t_op: &[i16; 2],
        lsp_flag: bool,
        synth: &mut [i16; L_FRAME],
        ana: &mut [i16; MAX_PRM_SIZE],
        mut ana_pos: usize,
    ) {
        let ov = &mut self.overflow;
        let mut xn = [0i16; L_SUBFR];
        let mut xn2 = [0i16; L_SUBFR];
        let mut code = [0i16; L_SUBFR];
        let mut y1 = [0i16; L_SUBFR];
        let mut y2 = [0i16; L_SUBFR];
        let mut g_coeff = [0i16; 6];
        let mut res = [0i16; L_SUBFR];
        let mut res2;
        // MR475: state of the even subframe, for redoing it once the joint gains are known.
        let mut xn_sf0 = [0i16; L_SUBFR];
        let mut y2_sf0 = [0i16; L_SUBFR];
        let mut code_sf0 = [0i16; L_SUBFR];
        let mut h1_sf0 = [0i16; L_SUBFR];
        let mut mem_syn_save = [0i16; M];
        let mut mem_w0_save = [0i16; M];
        let mut mem_err_save = [0i16; M];
        let mut sharp_save = 0i16;
        let (mut t0_sf0, mut t0_frac_sf0, mut i_subfr_sf0) = (0i16, 0i16, 0usize);
        let (mut gain_pit_sf0, mut gain_code_sf0) = (0i16, 0i16);

        let mut even = false;
        for sub_nr in 0..4 {
            let i_subfr = sub_nr * L_SUBFR;
            let a = &a_t[sub_nr * MP1..];
            let mut aq_off = sub_nr * MP1;
            even = !even;
            if even && mode == CodecMode::MR475 {
                mem_syn_save = self.mem_syn;
                mem_w0_save = self.mem_w0;
                mem_err_save.copy_from_slice(&self.mem_err[..M]);
                sharp_save = self.sharp;
            }
            let speech = &self.old_speech[SPEECH + i_subfr - M..];
            let mem_w0 = if mode != CodecMode::MR475 {
                &self.mem_w0
            } else {
                &mem_w0_save
            };
            subframe_pre_proc(
                mode,
                a,
                &aq_t[aq_off..],
                speech,
                &mut self.mem_err,
                mem_w0,
                &mut self.old_exc[EXC + i_subfr..],
                &mut self.h1,
                &mut xn,
                &mut res,
            );
            crate::stages::mark(crate::stages::E_SUB_PRE);
            if mode == CodecMode::MR475 && even {
                h1_sf0 = self.h1;
            }
            res2 = res;
            let cl = cl_ltp(
                &mut self.cl_ltp_st,
                &self.ton_stab_st,
                mode,
                i_subfr as i16,
                t_op,
                &self.h1,
                &mut self.old_exc,
                EXC + i_subfr,
                &mut res2,
                &xn,
                lsp_flag,
                &mut xn2,
                &mut y1,
                &mut g_coeff,
                &mut ana[..],
                &mut ana_pos,
                ov,
            );
            crate::stages::mark(crate::stages::E_CL);
            if sub_nr == 0 && self.ol_gain_flg[0] > 0 {
                self.old_lags[1] = cl.t0;
            }
            if sub_nr == 3 && self.ol_gain_flg[1] > 0 {
                self.old_lags[0] = cl.t0;
            }
            cbsearch(
                &xn2,
                &mut self.h1,
                cl.t0,
                self.sharp,
                cl.gain_pit,
                &res2,
                &mut code,
                &mut y2,
                &mut ana[..],
                &mut ana_pos,
                mode,
                sub_nr,
                ov,
            );
            crate::stages::mark(crate::stages::E_CB);
            let mut g = Gains {
                gain_pit: cl.gain_pit,
                gain_cod: 0,
            };
            gain_quant(
                &mut self.gain_quant_st,
                mode,
                &res,
                &self.old_exc[EXC + i_subfr..],
                &code,
                &xn,
                &xn2,
                &y1,
                &y2,
                &g_coeff,
                even,
                cl.gp_limit,
                &mut gain_pit_sf0,
                &mut gain_code_sf0,
                &mut g,
                &mut ana[..],
                &mut ana_pos,
                ov,
            );
            crate::stages::mark(crate::stages::E_GAIN);
            update_gp_clipping(&mut self.ton_stab_st, g.gain_pit);
            let speech = &self.old_speech[SPEECH..];
            if mode != CodecMode::MR475 {
                subframe_post_proc(
                    speech,
                    mode,
                    i_subfr,
                    g.gain_pit,
                    g.gain_cod,
                    &aq_t[aq_off..],
                    synth,
                    &xn,
                    &code,
                    &y1,
                    &y2,
                    &mut self.mem_syn,
                    &mut self.mem_err,
                    &mut self.mem_w0,
                    &mut self.old_exc[EXC..],
                    &mut self.sharp,
                );
                crate::stages::mark(crate::stages::E_SUB_POST);
            } else if even {
                i_subfr_sf0 = i_subfr;
                xn_sf0 = xn;
                y2_sf0 = y2;
                code_sf0 = code;
                t0_sf0 = cl.t0;
                t0_frac_sf0 = cl.t0_frac;
                subframe_post_proc(
                    speech,
                    mode,
                    i_subfr,
                    g.gain_pit,
                    g.gain_cod,
                    &aq_t[aq_off..],
                    synth,
                    &xn,
                    &code,
                    &y1,
                    &y2,
                    &mut mem_syn_save,
                    &mut self.mem_err,
                    &mut mem_w0_save,
                    &mut self.old_exc[EXC..],
                    &mut self.sharp,
                );
                crate::stages::mark(crate::stages::E_SUB_POST);
                self.sharp = sharp_save;
            } else {
                self.mem_err[..M].copy_from_slice(&mem_err_save);
                pred_lt_3or6(
                    &mut self.old_exc,
                    EXC + i_subfr_sf0,
                    t0_sf0,
                    t0_frac_sf0,
                    L_SUBFR as i16,
                    1,
                    ov,
                );
                convolve(&self.old_exc[EXC + i_subfr_sf0..], &h1_sf0, &mut y1);
                aq_off -= MP1;
                subframe_post_proc(
                    speech,
                    mode,
                    i_subfr_sf0,
                    gain_pit_sf0,
                    gain_code_sf0,
                    &aq_t[aq_off..],
                    synth,
                    &xn_sf0,
                    &code_sf0,
                    &y1,
                    &y2_sf0,
                    &mut self.mem_syn,
                    &mut self.mem_err,
                    &mut self.mem_w0,
                    &mut self.old_exc[EXC..],
                    &mut sharp_save,
                );
                crate::stages::mark(crate::stages::E_SUB_POST);
                aq_off += MP1;
                let speech_sf = &self.old_speech[SPEECH + i_subfr - M..];
                subframe_pre_proc(
                    mode,
                    a,
                    &aq_t[aq_off..],
                    speech_sf,
                    &mut self.mem_err,
                    &self.mem_w0,
                    &mut self.old_exc[EXC + i_subfr..],
                    &mut self.h1,
                    &mut xn,
                    &mut res,
                );
                crate::stages::mark(crate::stages::E_SUB_PRE);
                pred_lt_3or6(
                    &mut self.old_exc,
                    EXC + i_subfr,
                    cl.t0,
                    cl.t0_frac,
                    L_SUBFR as i16,
                    1,
                    ov,
                );
                convolve(&self.old_exc[EXC + i_subfr..], &self.h1, &mut y1);
                let speech = &self.old_speech[SPEECH..];
                subframe_post_proc(
                    speech,
                    mode,
                    i_subfr,
                    g.gain_pit,
                    g.gain_cod,
                    &aq_t[aq_off..],
                    synth,
                    &xn,
                    &code,
                    &y1,
                    &y2,
                    &mut self.mem_syn,
                    &mut self.mem_err,
                    &mut self.mem_w0,
                    &mut self.old_exc[EXC..],
                    &mut self.sharp,
                );
                crate::stages::mark(crate::stages::E_SUB_POST);
            }
        }
    }
}

/// `Prm2bits`: parameters to serial bits (MSB first per parameter).
fn prm2bits(mode: CodecMode, prm: &[i16; MAX_PRM_SIZE], bits: &mut [i16; MAX_SERIAL_SIZE]) {
    let bitno: &[i16] = match mode {
        CodecMode::MR475 => &BITNO_MR475,
        CodecMode::MR515 => &BITNO_MR515,
        CodecMode::MR59 => &BITNO_MR59,
        CodecMode::MR67 => &BITNO_MR67,
        CodecMode::MR74 => &BITNO_MR74,
        CodecMode::MR795 => &BITNO_MR795,
        CodecMode::MR102 => &BITNO_MR102,
        CodecMode::MR122 => &BITNO_MR122,
        CodecMode::MRDTX => &BITNO_MRDTX,
    };
    let mut pos = 0;
    for (&value, &n) in prm.iter().zip(&bitno[..PRMNO[mode as usize] as usize]) {
        let n = n as usize;
        for k in 0..n {
            bits[pos + k] = (value >> (n - 1 - k)) & 1;
        }
        pos += n;
    }
}

/// Discontinuous transmission (TS 26.093) and the voice activity detector (TS 26.094) it uses.
/// The 3GPP reference selects the VAD at compile time; here it is per encoder.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Dtx {
    /// Every frame is a speech frame.
    #[default]
    Off,
    /// DTX with VAD option 1.
    Vad1,
    /// DTX with VAD option 2.
    Vad2,
}

/// AMR-NB encoder (TS 26.073), bit-exact with the 3GPP reference, including encoder homing.
///
/// Input is 160 samples of 16-bit PCM at 8 kHz per frame (the codec uses the 13 most significant
/// bits). Output frames use the same representation as [`crate::rfc4867::Frame`], ready for
/// packing or decoding.
///
/// With DTX, silence becomes comfort-noise (SID) frames and `NoData` frames, which a
/// transmitter doesn't send:
///
/// ```
/// use amr_nb::{Dtx, Encoder, Mode, RxFrameType};
///
/// let mut enc = Encoder::new(Dtx::Vad1);
/// let silence = [0i16; 160];
/// let types: Vec<RxFrameType> = (0..50)
///     .map(|_| enc.encode(Mode::MR122, &silence).frame_type)
///     .collect();
/// assert_eq!(types[0], RxFrameType::SpeechGood); // the VAD's hangover
/// assert!(types.contains(&RxFrameType::SidFirst));
/// assert_eq!(types[49], RxFrameType::NoData);
/// ```
#[derive(Clone)]
pub struct Encoder {
    pre: PreProcessState,
    cod: Box<CodAmr>,
    sid: SidSyncState,
}

impl Encoder {
    /// A new encoder; `dtx` selects discontinuous transmission (comfort noise) and its voice
    /// activity detector.
    pub fn new(dtx: Dtx) -> Self {
        Self {
            pre: PreProcessState::default(),
            cod: Box::new(CodAmr::new(dtx)),
            sid: SidSyncState::default(),
        }
    }

    /// Back to the initial state (the DTX setting is kept).
    pub fn reset(&mut self) {
        self.pre = PreProcessState::default();
        self.cod.reset();
        self.sid.reset();
    }

    /// Encodes one frame at `mode`. The frame's serial bits, one per element (0 or 1, zero past
    /// the frame's size) in TS 26.073 parameter order, as [`crate::Decoder::decode`] reads
    /// them: `mode.bits()` for a speech frame, [`crate::SID_BITS`] for a SID frame. With DTX
    /// the frame can be `SidFirst`, `SidUpdate` or `NoData` (`mode` is then the mode
    /// indication); as in the reference's output, `SidFirst` and `NoData` frames also carry
    /// the comfort-noise analysis in their first [`crate::SID_BITS`] bits, which are not
    /// transmitted (a SID_FIRST frame's parameters are unused: [`crate::rfc4867`] writes them as
    /// zeros; a NO_DATA frame is not sent).
    pub fn encode(&mut self, mode: Mode, pcm: &[i16; L_FRAME]) -> Frame {
        crate::stages::start();
        let homing = pcm.iter().all(|&s| s == EHF);
        let mut speech = pcm.map(|s| s & !7);
        pre_process(&mut self.pre, &mut speech);
        let mut prm = [0i16; MAX_PRM_SIZE];
        let used_mode = self.cod.encode(mode.into(), &speech, &mut prm);
        let mut bits = [0i16; MAX_SERIAL_SIZE];
        prm2bits(used_mode, &prm, &mut bits);
        let frame_type = match sid_sync(&mut self.sid, used_mode) {
            TxFrameType::SpeechGood => RxFrameType::SpeechGood,
            TxFrameType::SidFirst => RxFrameType::SidFirst,
            TxFrameType::SidUpdate => RxFrameType::SidUpdate,
            TxFrameType::NoData => RxFrameType::NoData,
        };
        if homing {
            self.reset();
        }
        crate::stages::mark(crate::stages::E_TAIL);
        Frame {
            frame_type,
            mode,
            bits,
        }
    }
}

/// The DTX setting; the rest is the codec's internal memory.
impl std::fmt::Debug for Encoder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Encoder")
            .field("dtx", &self.cod.dtx_mode)
            .finish_non_exhaustive()
    }
}

impl Default for Encoder {
    /// An encoder without DTX.
    fn default() -> Self {
        Self::new(Dtx::Off)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testutil::Rng;

    /// TS 26.073 encoder homing: after an encoder homing frame the encoder continues exactly like
    /// a freshly created one (the 3GPP vectors can't show this: their homing frames come first).
    #[test]
    fn homing_frame_resets_the_encoder() {
        let mut g = Rng::new(1);
        for dtx in [Dtx::Off, Dtx::Vad1, Dtx::Vad2] {
            for mode in Mode::ALL {
                let mut used = Encoder::new(dtx);
                for _ in 0..20 {
                    let pcm: [i16; L_FRAME] = std::array::from_fn(|_| g.i16() >> 2);
                    used.encode(mode, &pcm);
                }
                used.encode(mode, &[EHF; L_FRAME]);
                let mut fresh = Encoder::new(dtx);
                for _ in 0..20 {
                    let pcm: [i16; L_FRAME] = std::array::from_fn(|_| g.i16() >> 2);
                    assert_eq!(
                        used.encode(mode, &pcm),
                        fresh.encode(mode, &pcm),
                        "{dtx:?} {mode:?}"
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod debug {
    use super::{Dtx, Encoder};
    use crate::Decoder;

    #[test]
    fn debug_shows_the_dtx_setting_and_hides_the_state() {
        assert_eq!(
            format!("{:?}", Encoder::new(Dtx::Vad2)),
            "Encoder { dtx: Vad2, .. }"
        );
        assert_eq!(format!("{:?}", Decoder::new()), "Decoder { .. }");
    }
}
