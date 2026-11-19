//! The AMR-NB decoder: `Decoder_amr` (`dec_amr.cpp`), the frame wrapper (`sp_dec.cpp`) and the
//! TS 26.073 decoder homing, which opencore lacks.
//! Bit-exact port of opencore-amr 0.1.6 (Apache-2.0); conformance-tested against the official
//! 3GPP TS 26.074 test vectors (`just conform`).

use crate::basic_op::{Flag, MAX_16, MIN_16, add, l_mac_energy, l_mult, l_shr, mult, round, shl};
use crate::cnst::{CodecMode, L_FRAME, L_INTERPOL, L_SUBFR, M, MP1, Mode, PIT_MAX};
use crate::dtx::{DtxDecDeps, DtxDecState, DtxState, build_cn_param, dtx_dec};
use crate::excitation::{
    PhDispState, dec_8i40_31bits, dec_10i40_35bits, dec_lag3, dec_lag6, decode_2i40_9bits,
    decode_2i40_11bits, decode_3i40_14bits, decode_4i40_17bits, ph_disp, ph_disp_lock,
    ph_disp_release, pred_lt_3or6, total_exc,
};
use crate::gains::{
    BgnScdState, CbGainAverageState, EcGainCodeState, EcGainPitchState, GcPredState, bgn_scd,
    cb_gain_average, d_gain_code, d_gain_pitch, dec_gain, ec_gain_code, ec_gain_code_update,
    ec_gain_pitch, ec_gain_pitch_update, ex_ctrl,
};
use crate::lpc::{
    DPlsfState, LspAvgState, d_plsf_3, d_plsf_5, int_lpc_1and3, int_lpc_1to3, int_lsf, lsf_lsp,
    lsp_avg,
};
use crate::post::filt::syn_filt;
use crate::post::post_pro::{PostProcessState, post_process};
use crate::post::pstfilt::{PostFilterState, post_filter};
use crate::tables;

pub use crate::dtx::RxFrameType;

/// Serial bits per frame (the largest mode, MR122, uses all 244).
pub const MAX_SERIAL_SIZE: usize = 244;
const MAX_PRM_SIZE: usize = 57;
const AZ_SIZE: usize = 4 * MP1;
const EXC_ENERGY_HIST_LEN: usize = 9;
const LTP_GAIN_HISTORY_LEN: usize = 9;
const SHARPMAX: i16 = 13017;
const SHARPMIN: i16 = 0;
const PIT_MIN: i16 = 20;
const PIT_MIN_MR122: i16 = 18;
const L_FRAME_BY2: usize = 80;
/// `old_exc` length and the offset of the current excitation (`exc`) inside it.
const OLD_EXC_LEN: usize = L_SUBFR + PIT_MAX + L_INTERPOL;
const EXC: usize = PIT_MAX + L_INTERPOL;
/// Encoder homing frame sample value (TS 26.073 `EHF_MASK`).
const EHF: i16 = 0x0008;

/// `Decoder_amrState`.
#[derive(Clone)]
struct DecoderAmr {
    old_exc: [i16; OLD_EXC_LEN],
    lsp_old: [i16; M],
    mem_syn: [i16; M],
    sharp: i16,
    old_t0: i16,
    prev_bf: i16,
    prev_pdf: i16,
    state: i16,
    exc_energy_hist: [i16; EXC_ENERGY_HIST_LEN],
    t0_lag_buff: i16,
    in_background_noise: i16,
    voiced_hangover: i16,
    ltp_gain_history: [i16; LTP_GAIN_HISTORY_LEN],
    background_state: BgnScdState,
    nodata_seed: i16,
    cb_gain_aver_state: CbGainAverageState,
    lsp_avg_st: LspAvgState,
    lsf_state: DPlsfState,
    ec_gain_p_st: EcGainPitchState,
    ec_gain_c_st: EcGainCodeState,
    pred_state: GcPredState,
    ph_disp_st: PhDispState,
    dtx_decoder_state: DtxDecState,
    overflow: Flag,
}

const LSP_INIT: [i16; M] = [
    30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000,
];

/// `(x << 1)` saturated to 16 bits, as the C code spells it out for `pit_sharp`.
fn sat_shl1(x: i16) -> i16 {
    let t = (x as i32) << 1;
    if t != t as i16 as i32 {
        if x > 0 { MAX_16 } else { MIN_16 }
    } else {
        t as i16
    }
}

impl DecoderAmr {
    /// `Decoder_amr_init` (everything it resets is reset again by `reset(MR475)`).
    fn new() -> Self {
        let mut s = Self {
            old_exc: [0; OLD_EXC_LEN],
            lsp_old: [0; M],
            mem_syn: [0; M],
            sharp: 0,
            old_t0: 0,
            prev_bf: 0,
            prev_pdf: 0,
            state: 0,
            exc_energy_hist: [0; EXC_ENERGY_HIST_LEN],
            t0_lag_buff: 40,
            in_background_noise: 0,
            voiced_hangover: 0,
            ltp_gain_history: [0; LTP_GAIN_HISTORY_LEN],
            background_state: BgnScdState::default(),
            nodata_seed: 0,
            cb_gain_aver_state: CbGainAverageState::default(),
            lsp_avg_st: LspAvgState::default(),
            lsf_state: DPlsfState::default(),
            ec_gain_p_st: EcGainPitchState::default(),
            ec_gain_c_st: EcGainCodeState::default(),
            pred_state: GcPredState::default(),
            ph_disp_st: PhDispState::default(),
            dtx_decoder_state: DtxDecState::new(),
            overflow: false,
        };
        s.reset(CodecMode::MR475);
        s
    }

    /// `Decoder_amr_reset`. With `MRDTX`, speech history needed by comfort noise is kept.
    fn reset(&mut self, mode: CodecMode) {
        let dtx = mode == CodecMode::MRDTX;
        // C zeroes only the history part (`PIT_MAX + L_INTERPOL`), not the current subframe.
        self.old_exc[..EXC].fill(0);
        if !dtx {
            self.mem_syn = [0; M];
        }
        self.sharp = SHARPMIN;
        self.old_t0 = 40;
        self.overflow = false;
        if !dtx {
            self.lsp_old = LSP_INIT;
        }
        self.prev_bf = 0;
        self.prev_pdf = 0;
        self.state = 0;
        self.t0_lag_buff = 40;
        self.in_background_noise = 0;
        self.voiced_hangover = 0;
        if !dtx {
            self.exc_energy_hist = [0; EXC_ENERGY_HIST_LEN];
        }
        self.ltp_gain_history = [0; LTP_GAIN_HISTORY_LEN];
        self.cb_gain_aver_state.reset();
        if !dtx {
            self.lsp_avg_st.reset();
        }
        self.lsf_state.reset();
        self.ec_gain_p_st.reset();
        self.ec_gain_c_st.reset();
        if !dtx {
            self.pred_state.reset();
        }
        self.background_state.reset();
        self.nodata_seed = 21845;
        self.ph_disp_st.reset();
        if !dtx {
            self.dtx_decoder_state.reset();
        }
    }

    /// `Decoder_amr`: decodes one frame of parameters into `synth`; `a_t` receives the four
    /// subframe LP filters for the post-filter.
    fn decode(
        &mut self,
        mode: CodecMode,
        parm: &[i16],
        frame_type: RxFrameType,
        synth: &mut [i16; L_FRAME],
        a_t: &mut [i16; AZ_SIZE],
    ) {
        use CodecMode::*;
        let mut ov = self.overflow;
        let new_dtx_state = self.dtx_decoder_state.rx_dtx_handler(frame_type, &mut ov);

        if new_dtx_state != DtxState::Speech {
            self.overflow = ov;
            self.reset(MRDTX);
            let mut ov = self.overflow;
            dtx_dec(
                &mut self.dtx_decoder_state,
                DtxDecDeps {
                    mem_syn: &mut self.mem_syn,
                    lsf_state: &mut self.lsf_state,
                    pred_state: &mut self.pred_state,
                    aver_state: &mut self.cb_gain_aver_state,
                },
                new_dtx_state,
                mode,
                parm,
                synth,
                a_t,
                &mut ov,
            );
            lsf_lsp(&self.lsf_state.past_lsf_q, &mut self.lsp_old, M, &mut ov);
            lsp_avg(&mut self.lsp_avg_st, &self.lsf_state.past_lsf_q, &mut ov);
            self.dtx_decoder_state.dtx_global_state = new_dtx_state;
            self.overflow = ov;
            return;
        }

        // Bad / lost frame handling.
        let mut parm_buf = [0i16; MAX_PRM_SIZE + 1];
        parm_buf[..parm.len().min(MAX_PRM_SIZE + 1)]
            .copy_from_slice(&parm[..parm.len().min(MAX_PRM_SIZE + 1)]);
        let (mut bfi, mut pdfi) = (0i16, 0i16);
        match frame_type {
            RxFrameType::SpeechBad | RxFrameType::NoData | RxFrameType::Onset => {
                bfi = 1;
                if matches!(frame_type, RxFrameType::NoData | RxFrameType::Onset) {
                    let m = mode as usize;
                    build_cn_param(
                        &mut self.nodata_seed,
                        tables::PRMNO[m] as usize,
                        bitno(mode),
                        &mut parm_buf,
                        &tables::WINDOW_200_40,
                        &mut ov,
                    );
                }
            }
            RxFrameType::SpeechDegraded => pdfi = 1,
            _ => {}
        }
        if bfi != 0 {
            self.state += 1;
        } else if self.state == 6 {
            self.state = 5;
        } else {
            self.state = 0;
        }
        if self.state > 6 {
            self.state = 6;
        }
        match self.dtx_decoder_state.dtx_global_state {
            DtxState::Dtx => {
                self.state = 5;
                self.prev_bf = 0;
            }
            DtxState::DtxMute => {
                self.state = 5;
                self.prev_bf = 1;
            }
            DtxState::Speech => {}
        }

        // LSP decoding and interpolation.
        let prev_lsf = self.lsf_state.past_lsf_q;
        let mut lsp_new = [0i16; M];
        let mut p = 0usize; // index into parm_buf (the C `parm` pointer)
        if mode != MR122 {
            d_plsf_3(
                &mut self.lsf_state,
                mode,
                bfi,
                &parm_buf,
                &mut lsp_new,
                &mut ov,
            );
            p += 3;
            int_lpc_1to3(&self.lsp_old, &lsp_new, a_t, &mut ov);
        } else {
            let mut lsp_mid = [0i16; M];
            d_plsf_5(
                &mut self.lsf_state,
                bfi,
                &parm_buf,
                &mut lsp_mid,
                &mut lsp_new,
                &mut ov,
            );
            p += 5;
            int_lpc_1and3(&self.lsp_old, &lsp_mid, &lsp_new, a_t, &mut ov);
        }
        self.lsp_old = lsp_new;
        crate::stages::mark(crate::stages::D_LSP);

        let mut code = [0i16; L_SUBFR];
        let mut excp = [0i16; L_SUBFR];
        let mut exc_enhanced = [0i16; L_SUBFR];
        let mut lsf_i = [0i16; M];
        let mut t0 = 0i16; // C: uninitialised; Dec_lag6 only reads it for relative lags
        let mut t0_frac;
        let mut index_mr475 = 0i16;
        let mut gain_pit = 0i16;
        let mut gain_code;
        let mut even_subfr = 0i16;
        let mut subfr_nr = -1i16;

        for i_subfr in (0..L_FRAME).step_by(L_SUBFR) {
            subfr_nr += 1;
            even_subfr = 1 - even_subfr;
            let mut pit_flag = i_subfr as i16;
            if i_subfr == L_FRAME_BY2 && mode != MR475 && mode != MR515 {
                pit_flag = 0;
            }

            // Pitch lag and adaptive codebook.
            let mut index = parm_buf[p];
            p += 1;
            if mode != MR122 {
                let flag4 = matches!(mode, MR475 | MR515 | MR59 | MR67) as i16;
                let (delta_frc_low, delta_frc_range) =
                    if mode == MR795 { (10, 19) } else { (5, 9) };
                let mut t0_min = (self.old_t0 - delta_frc_low).max(PIT_MIN);
                let mut t0_max = t0_min + delta_frc_range;
                if t0_max > PIT_MAX as i16 {
                    t0_max = PIT_MAX as i16;
                    t0_min = t0_max - delta_frc_range;
                }
                (t0, t0_frac) =
                    dec_lag3(index, t0_min, t0_max, pit_flag, self.old_t0, flag4, &mut ov);
                self.t0_lag_buff = t0;
                if bfi != 0 {
                    if self.old_t0 < PIT_MAX as i16 {
                        self.old_t0 += 1;
                    }
                    t0 = self.old_t0;
                    t0_frac = 0;
                    if self.in_background_noise != 0
                        && self.voiced_hangover > 4
                        && matches!(mode, MR475 | MR515 | MR59)
                    {
                        t0 = self.t0_lag_buff;
                    }
                }
                pred_lt_3or6(
                    &mut self.old_exc,
                    EXC,
                    t0,
                    t0_frac,
                    L_SUBFR as i16,
                    1,
                    &mut ov,
                );
            } else {
                t0_frac = dec_lag6(
                    index,
                    PIT_MIN_MR122,
                    PIT_MAX as i16,
                    pit_flag,
                    &mut t0,
                    &mut ov,
                );
                if !(bfi == 0 && (pit_flag == 0 || index < 61)) {
                    self.t0_lag_buff = t0;
                    t0 = self.old_t0;
                    t0_frac = 0;
                }
                pred_lt_3or6(
                    &mut self.old_exc,
                    EXC,
                    t0,
                    t0_frac,
                    L_SUBFR as i16,
                    0,
                    &mut ov,
                );
            }

            // Innovative (algebraic) codebook.
            let mut pit_sharp;
            match mode {
                MR475 | MR515 => {
                    index = parm_buf[p];
                    let sign = parm_buf[p + 1];
                    p += 2;
                    decode_2i40_9bits(subfr_nr, sign, index, &mut code, &mut ov);
                    pit_sharp = sat_shl1(self.sharp);
                }
                MR59 => {
                    index = parm_buf[p];
                    let sign = parm_buf[p + 1];
                    p += 2;
                    decode_2i40_11bits(sign, index, &mut code);
                    pit_sharp = sat_shl1(self.sharp);
                }
                MR67 => {
                    index = parm_buf[p];
                    let sign = parm_buf[p + 1];
                    p += 2;
                    decode_3i40_14bits(sign, index, &mut code);
                    pit_sharp = sat_shl1(self.sharp);
                }
                MR74 | MR795 => {
                    index = parm_buf[p];
                    let sign = parm_buf[p + 1];
                    p += 2;
                    decode_4i40_17bits(sign, index, &mut code);
                    pit_sharp = sat_shl1(self.sharp);
                }
                MR102 => {
                    dec_8i40_31bits(&parm_buf[p..p + 7], &mut code, &mut ov);
                    p += 7;
                    pit_sharp = sat_shl1(self.sharp);
                }
                _ => {
                    // MR122: pitch gain is coded per subframe before the codebook.
                    index = parm_buf[p];
                    p += 1;
                    gain_pit = if bfi != 0 {
                        ec_gain_pitch(&self.ec_gain_p_st, self.state, &mut ov)
                    } else {
                        d_gain_pitch(mode, index)
                    };
                    ec_gain_pitch_update(
                        &mut self.ec_gain_p_st,
                        bfi,
                        self.prev_bf,
                        &mut gain_pit,
                        &mut ov,
                    );
                    dec_10i40_35bits(&parm_buf[p..p + 10], &mut code);
                    p += 10;
                    pit_sharp = sat_shl1(gain_pit);
                }
            }
            // Pitch sharpening of the innovation.
            for i in (t0.max(0) as usize)..L_SUBFR {
                let temp = mult(code[i - t0 as usize], pit_sharp, &mut ov);
                code[i] = add(code[i], temp, &mut ov);
            }

            // Gains.
            if mode == MR475 {
                if even_subfr != 0 {
                    index_mr475 = parm_buf[p];
                    p += 1;
                }
                if bfi == 0 {
                    (gain_pit, gain_code) = dec_gain(
                        &mut self.pred_state,
                        mode,
                        index_mr475,
                        &code,
                        even_subfr,
                        &mut ov,
                    );
                } else {
                    gain_pit = ec_gain_pitch(&self.ec_gain_p_st, self.state, &mut ov);
                    gain_code = ec_gain_code(
                        &self.ec_gain_c_st,
                        &mut self.pred_state,
                        self.state,
                        &mut ov,
                    );
                }
                ec_gain_pitch_update(
                    &mut self.ec_gain_p_st,
                    bfi,
                    self.prev_bf,
                    &mut gain_pit,
                    &mut ov,
                );
                ec_gain_code_update(
                    &mut self.ec_gain_c_st,
                    bfi,
                    self.prev_bf,
                    &mut gain_code,
                    &mut ov,
                );
                pit_sharp = gain_pit.min(SHARPMAX);
            } else if mode <= MR74 || mode == MR102 {
                index = parm_buf[p];
                p += 1;
                if bfi == 0 {
                    (gain_pit, gain_code) = dec_gain(
                        &mut self.pred_state,
                        mode,
                        index,
                        &code,
                        even_subfr,
                        &mut ov,
                    );
                } else {
                    gain_pit = ec_gain_pitch(&self.ec_gain_p_st, self.state, &mut ov);
                    gain_code = ec_gain_code(
                        &self.ec_gain_c_st,
                        &mut self.pred_state,
                        self.state,
                        &mut ov,
                    );
                }
                ec_gain_pitch_update(
                    &mut self.ec_gain_p_st,
                    bfi,
                    self.prev_bf,
                    &mut gain_pit,
                    &mut ov,
                );
                ec_gain_code_update(
                    &mut self.ec_gain_c_st,
                    bfi,
                    self.prev_bf,
                    &mut gain_code,
                    &mut ov,
                );
                pit_sharp = gain_pit.min(SHARPMAX);
                if mode == MR102 && self.old_t0 > (L_SUBFR + 5) as i16 {
                    pit_sharp >>= 2;
                }
            } else {
                index = parm_buf[p];
                p += 1;
                if mode == MR795 {
                    gain_pit = if bfi != 0 {
                        ec_gain_pitch(&self.ec_gain_p_st, self.state, &mut ov)
                    } else {
                        d_gain_pitch(mode, index)
                    };
                    ec_gain_pitch_update(
                        &mut self.ec_gain_p_st,
                        bfi,
                        self.prev_bf,
                        &mut gain_pit,
                        &mut ov,
                    );
                    index = parm_buf[p];
                    p += 1;
                }
                gain_code = if bfi == 0 {
                    d_gain_code(&mut self.pred_state, mode, index, &code, &mut ov)
                } else {
                    ec_gain_code(
                        &self.ec_gain_c_st,
                        &mut self.pred_state,
                        self.state,
                        &mut ov,
                    )
                };
                ec_gain_code_update(
                    &mut self.ec_gain_c_st,
                    bfi,
                    self.prev_bf,
                    &mut gain_code,
                    &mut ov,
                );
                pit_sharp = if mode == MR795 {
                    gain_pit.min(SHARPMAX)
                } else {
                    gain_pit
                };
            }

            if mode != MR475 || even_subfr == 0 {
                self.sharp = gain_pit.min(SHARPMAX);
            }
            pit_sharp = shl(pit_sharp, 1, &mut ov);
            if pit_sharp > 16384 {
                for (e, &x) in excp.iter_mut().zip(&self.old_exc[EXC..]) {
                    let temp = mult(x, pit_sharp, &mut ov);
                    let mut t = l_mult(temp, gain_pit, &mut ov);
                    if mode == MR122 {
                        t >>= 1;
                    }
                    *e = round(t, &mut ov);
                }
            }

            if bfi == 0 {
                self.ltp_gain_history.copy_within(1.., 0);
                self.ltp_gain_history[8] = gain_pit;
            }
            if (self.prev_bf != 0 || bfi != 0)
                && self.in_background_noise != 0
                && matches!(mode, MR475 | MR515 | MR59)
            {
                if gain_pit > 12288 {
                    gain_pit = ((gain_pit - 12288) >> 1) + 12288;
                }
                gain_pit = gain_pit.min(14745);
            }

            // Codebook gain smoothing in stationary background noise.
            int_lsf(
                &prev_lsf,
                &self.lsf_state.past_lsf_q,
                i_subfr as i16,
                &mut lsf_i,
                &mut ov,
            );
            let mut gain_code_mix = cb_gain_average(
                &mut self.cb_gain_aver_state,
                mode,
                gain_code,
                &lsf_i,
                &self.lsp_avg_st.lsp_mean_save,
                bfi,
                self.prev_bf,
                pdfi,
                self.prev_pdf,
                self.in_background_noise,
                self.voiced_hangover,
                &mut ov,
            );
            if mode > MR67 && mode != MR102 {
                gain_code_mix = gain_code;
            }

            // Total excitation.
            let (pitch_fac, tmp_shift) = if mode <= MR102 {
                (gain_pit, 1)
            } else {
                (gain_pit >> 1, 2)
            };
            exc_enhanced[..L_SUBFR].copy_from_slice(&self.old_exc[EXC..EXC + L_SUBFR]);
            total_exc::<true>(
                (&mut self.old_exc[EXC..EXC + L_SUBFR]).try_into().unwrap(),
                &code,
                pitch_fac,
                gain_code,
                tmp_shift,
            );

            // Phase dispersion (anti-sparseness).
            ph_disp_release(&mut self.ph_disp_st);
            if matches!(mode, MR475 | MR515 | MR59)
                && self.voiced_hangover > 3
                && self.in_background_noise != 0
                && bfi != 0
            {
                ph_disp_lock(&mut self.ph_disp_st);
            }
            crate::stages::mark(crate::stages::D_PARAMS);
            ph_disp(
                &mut self.ph_disp_st,
                mode,
                &mut exc_enhanced,
                gain_code_mix,
                gain_pit,
                &mut code,
                pitch_fac,
                tmp_shift,
                &mut ov,
            );
            crate::stages::mark(crate::stages::D_EXC);

            // Excitation energy (for error concealment).
            let l_temp = l_mac_energy(&exc_enhanced, &mut ov) >> 1;
            let (sq, exp) = crate::math::sqrt_l_exp(l_temp, &mut ov);
            let l_temp = l_shr(sq, (exp >> 1) + 15, &mut ov);
            let exc_energy = (l_temp >> 2) as i16;

            if matches!(mode, MR475 | MR515 | MR59)
                && self.voiced_hangover > 5
                && self.in_background_noise != 0
                && self.state < 4
                && ((pdfi != 0 && self.prev_pdf != 0) || bfi != 0 || self.prev_bf != 0)
            {
                let careful_flag = (pdfi != 0 && bfi == 0) as i16;
                ex_ctrl(
                    &mut exc_enhanced,
                    exc_energy,
                    &self.exc_energy_hist,
                    self.voiced_hangover,
                    self.prev_bf,
                    careful_flag,
                    &mut ov,
                );
            }
            if !(self.in_background_noise != 0 && (bfi != 0 || self.prev_bf != 0) && self.state < 4)
            {
                self.exc_energy_hist.copy_within(1.., 0);
                self.exc_energy_hist[8] = exc_energy;
            }

            // Synthesis.
            let az = &a_t[subfr_nr as usize * MP1..(subfr_nr as usize + 1) * MP1];
            let out = &mut synth[i_subfr..i_subfr + L_SUBFR];
            if pit_sharp > 16384 {
                for (e, &x) in excp.iter_mut().zip(&exc_enhanced) {
                    *e = add(*e, x, &mut ov);
                }
                crate::post::agc::agc2(&exc_enhanced, &mut excp, L_SUBFR, &mut ov);
                ov = syn_filt(az, &excp, out, L_SUBFR, &mut self.mem_syn, false);
            } else {
                ov = syn_filt(az, &exc_enhanced, out, L_SUBFR, &mut self.mem_syn, false);
            }
            // TS 26.073: on overflow (the reference's global flag, which `Syn_filt` sets) scale
            // the excitation down and re-synthesise. (opencore's Syn_filt reports nothing: this
            // branch was unreachable there.)
            if ov {
                for x in self.old_exc.iter_mut() {
                    *x >>= 2;
                }
                for x in exc_enhanced.iter_mut() {
                    *x >>= 2;
                }
                syn_filt(az, &exc_enhanced, out, L_SUBFR, &mut self.mem_syn, true);
            } else {
                self.mem_syn.copy_from_slice(&out[L_SUBFR - M..]);
            }
            self.old_exc.copy_within(L_SUBFR.., 0);
            self.old_t0 = t0;
            crate::stages::mark(crate::stages::D_SYNTH);
        }

        self.in_background_noise = bgn_scd(
            &mut self.background_state,
            &self.ltp_gain_history,
            synth,
            &mut self.voiced_hangover,
            &mut ov,
        );
        self.dtx_decoder_state
            .activity_update(&self.lsf_state.past_lsf_q, synth, &mut ov);
        self.prev_bf = bfi;
        self.prev_pdf = pdfi;
        lsp_avg(&mut self.lsp_avg_st, &self.lsf_state.past_lsf_q, &mut ov);
        self.dtx_decoder_state.dtx_global_state = new_dtx_state;
        self.overflow = ov;
    }
}

/// Bits per parameter for `mode` (`bitno_ptr[mode]`).
fn bitno(mode: CodecMode) -> &'static [i16] {
    match mode {
        CodecMode::MR475 => &tables::BITNO_MR475,
        CodecMode::MR515 => &tables::BITNO_MR515,
        CodecMode::MR59 => &tables::BITNO_MR59,
        CodecMode::MR67 => &tables::BITNO_MR67,
        CodecMode::MR74 => &tables::BITNO_MR74,
        CodecMode::MR795 => &tables::BITNO_MR795,
        CodecMode::MR102 => &tables::BITNO_MR102,
        CodecMode::MR122 => &tables::BITNO_MR122,
        CodecMode::MRDTX => &tables::BITNO_MRDTX,
    }
}

/// `Bits2prm`: serial bits (one per `i16`) to codec parameters for `mode`.
fn bits2prm(
    mode: CodecMode,
    bits: &[i16; MAX_SERIAL_SIZE],
    prm: &mut [i16; MAX_PRM_SIZE + 1],
) -> usize {
    let mut b = 0;
    let widths = bitno(mode);
    for (p, &w) in prm.iter_mut().zip(widths) {
        let mut v: i16 = 0;
        for &bit in &bits[b..b + w as usize] {
            v = (v << 1) | bit;
        }
        *p = v;
        b += w as usize;
    }
    widths.len()
}

/// Speech-mode parameters of the decoder homing frame per mode (TS 26.073 §5.3): the encoder's
/// output for an encoder homing frame from the reset state. Verified against opencore's encoder
/// by `homing::dhf_params_match_opencore_encoder`.
#[rustfmt::skip]
const DHF_PARAMS: [&[i16]; 8] = include!("dhf_params.in");

/// A bit-exact AMR-NB decoder (3GPP TS 26.073), including decoder homing.
///
/// Input is one frame of serial bits (one bit per `i16`, as in the 3GPP test vectors), its
/// receive frame type and codec mode; output is 160 samples of 13-bit PCM (left-aligned in `i16`).
///
/// A lost frame is decoded as [`RxFrameType::NoData`] (or [`RxFrameType::SpeechBad`] when its
/// bits arrived damaged): the decoder conceals it from its state.
///
/// ```
/// use amr_nb::{Decoder, MAX_SERIAL_SIZE, Mode, RxFrameType};
///
/// let mut dec = Decoder::new();
/// let mut pcm = [0i16; 160];
/// let bits = [0i16; MAX_SERIAL_SIZE]; // one bit per element, as received
/// dec.decode(RxFrameType::SpeechGood, Mode::MR122, &bits, &mut pcm);
/// // The next frame was lost in transit: conceal it.
/// dec.decode(RxFrameType::NoData, Mode::MR122, &[0; MAX_SERIAL_SIZE], &mut pcm);
/// ```
#[derive(Clone)]
pub struct Decoder {
    amr: DecoderAmr,
    post_state: PostFilterState,
    post_hp_state: PostProcessState,
    prev_mode: Mode,
    homed: bool,
}

/// `Decoder { .. }`: the state is the codec's internal memory.
impl std::fmt::Debug for Decoder {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Decoder").finish_non_exhaustive()
    }
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Decoder {
    /// A decoder in the reset (homed) state.
    pub fn new() -> Self {
        let mut d = Self {
            amr: DecoderAmr::new(),
            post_state: PostFilterState::default(),
            post_hp_state: PostProcessState::default(),
            prev_mode: Mode::MR475,
            homed: true,
        };
        d.reset();
        d
    }

    /// `Speech_Decode_Frame_reset`: back to the initial state.
    pub fn reset(&mut self) {
        self.amr.reset(CodecMode::MR475);
        self.post_state.reset();
        self.post_hp_state.reset();
        self.prev_mode = Mode::MR475;
        self.homed = true;
    }

    /// Decodes one frame. `bits` holds the frame's serial bits, one per element (0 or 1), in
    /// TS 26.073 parameter order: `mode.bits()` of them for speech frames, [`SID_BITS`] for
    /// SID frames, none read for [`RxFrameType::NoData`]. `mode` is ignored for NO_DATA
    /// frames: the last frame's is used (for a SID frame, `mode` is its mode indication), as
    /// the TS 26.073 reference does.
    ///
    /// [`SID_BITS`]: crate::SID_BITS
    pub fn decode(
        &mut self,
        frame_type: RxFrameType,
        mode: Mode,
        bits: &[i16; MAX_SERIAL_SIZE],
        out: &mut [i16; L_FRAME],
    ) {
        crate::stages::start();
        let mode = if frame_type == RxFrameType::NoData {
            self.prev_mode
        } else {
            mode
        };
        let mut prm = [0i16; MAX_PRM_SIZE + 1];
        let dhf = DHF_PARAMS[mode as usize];
        let first_sf = tables::PRMNOFSF[mode as usize] as usize;
        let is_dhf = |prm: &[i16], n: usize| prm[..n] == dhf[..n];

        bits2prm(mode.into(), bits, &mut prm);
        // Homed: a homing frame (checked up to the end of subframe 1) produces the EHF output.
        let mut reset = self.homed && is_dhf(&prm, first_sf);
        if reset {
            out.fill(EHF);
        } else {
            // Reuse the homing check's parameters when decode_frame would parse the same.
            self.decode_frame(frame_type, mode.into(), bits, &prm, out);
        }
        if !self.homed {
            reset = is_dhf(&prm, dhf.len());
        }
        if reset {
            self.reset();
        } else {
            self.prev_mode = mode;
        }
        self.homed = reset;
    }

    /// `GSMFrameDecode` (opencore `sp_dec.cpp`).
    fn decode_frame(
        &mut self,
        frame_type: RxFrameType,
        mode: CodecMode,
        bits: &[i16; MAX_SERIAL_SIZE],
        parsed: &[i16; MAX_PRM_SIZE + 1],
        synth: &mut [i16; L_FRAME],
    ) {
        let mut parm = *parsed;
        // `parsed` is bits2prm(mode, bits): valid unless a SID frame needs MRDTX's layout.
        if matches!(frame_type, RxFrameType::SidBad | RxFrameType::SidUpdate) {
            bits2prm(CodecMode::MRDTX, bits, &mut parm);
        }
        let mut az = [0i16; AZ_SIZE];
        self.amr.decode(mode, &parm, frame_type, synth, &mut az);
        crate::stages::mark(crate::stages::D_TAIL);
        let mut ov = self.amr.overflow;
        post_filter(&mut self.post_state, mode, synth, &az, &mut ov);
        crate::stages::mark(crate::stages::D_POSTFILTER);
        post_process(&mut self.post_hp_state, synth, L_FRAME, &mut ov);
        self.amr.overflow = ov;
        for s in synth.iter_mut() {
            *s &= !7; // 13-bit output (0xfff8)
        }
        crate::stages::mark(crate::stages::D_POSTPROC);
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod homing {
    use super::*;
    use opencore_sys as _;
    use std::ffi::c_void;
    use std::ptr::null_mut;

    unsafe extern "C" {
        fn AMREncodeInit(enc: *mut *mut c_void, sid: *mut *mut c_void, dtx: i32) -> i16;
        fn AMREncodeExit(enc: *mut *mut c_void, sid: *mut *mut c_void);
        fn AMREncode(
            enc: *mut c_void,
            sid: *mut c_void,
            mode: i32,
            input: *mut i16,
            out: *mut u8,
            frame_type: *mut i32,
            format: i16,
        ) -> i16;
    }

    /// DHF parameters from opencore's encoder (ETS output format = 2).
    fn dhf_from_encoder(mode: CodecMode) -> Vec<i16> {
        let (mut enc, mut sid) = (null_mut(), null_mut());
        let mut ehf = [EHF; L_FRAME];
        let mut out = [0i16; MAX_SERIAL_SIZE + 2];
        let mut ft = 0;
        unsafe {
            AMREncodeInit(&mut enc, &mut sid, 0);
            AMREncode(
                enc,
                sid,
                mode as i32,
                ehf.as_mut_ptr(),
                out.as_mut_ptr().cast(),
                &mut ft,
                2,
            );
            AMREncodeExit(&mut enc, &mut sid);
        }
        let bits: [i16; MAX_SERIAL_SIZE] = out[1..=MAX_SERIAL_SIZE].try_into().unwrap();
        let mut prm = [0i16; MAX_PRM_SIZE + 1];
        let n = bits2prm(mode, &bits, &mut prm);
        prm[..n].to_vec()
    }

    #[test]
    fn dhf_params_match_opencore_encoder() {
        for mode in &CodecMode::ALL[..8] {
            assert_eq!(
                DHF_PARAMS[*mode as usize],
                dhf_from_encoder(*mode),
                "{mode:?}"
            );
        }
    }

    /// Regenerates `dhf_params.in`: `cargo test -p amr-nb print_dhf_params -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn print_dhf_params() {
        println!(
            "// @generated by decoder::homing::print_dhf_params (opencore encoder output for the EHF)"
        );
        println!("[");
        for mode in &CodecMode::ALL[..8] {
            println!("    &{:?},", dhf_from_encoder(*mode));
        }
        println!("]");
    }
}

/// The decoder against the TS 26.073 reference decoder (`Speech_Decode_Frame`, as its
/// `decoder.c` drives it) on random frames: every RX frame type (bad and degraded speech, SID
/// first/update/bad, NO_DATA: the concealment paths), random modes and bits, mode switches.
/// Random bits are never a decoder homing frame, so the reference's homing logic is idle.
#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::testutil::Rng;
    use std::ffi::c_void;
    use std::ptr::null_mut;

    unsafe extern "C" {
        fn ref_Speech_Decode_Frame_init(st: *mut *mut c_void, id: *const u8) -> i32;
        fn ref_Speech_Decode_Frame_exit(st: *mut *mut c_void);
        fn ref_Speech_Decode_Frame(
            st: *mut c_void,
            mode: i32,
            serial: *mut i16,
            frame_type: i32,
            synth: *mut i16,
        ) -> i32;
    }

    fn run(seed: u64, frames: usize) {
        use RxFrameType::*;
        let _l = ref26073_sys::lock();
        let mut g = Rng::new(seed);
        let mut dec = Decoder::new();
        let mut st = null_mut();
        unsafe { ref_Speech_Decode_Frame_init(&mut st, c"Decoder".as_ptr().cast()) };
        let types = [
            SpeechGood,
            SpeechGood,
            SpeechGood,
            SpeechGood,
            SpeechDegraded,
            Onset,
            SpeechBad,
            SidFirst,
            SidUpdate,
            SidBad,
            NoData,
        ];
        let (mut mode, mut prev_mode) = (Mode::MR122, Mode::MR122);
        let mut ft = SpeechGood;
        for f in 0..frames {
            if g.u64().is_multiple_of(8) {
                mode = Mode::ALL[(g.u64() % 8) as usize];
            }
            if g.u64().is_multiple_of(4) {
                ft = types[(g.u64() % types.len() as u64) as usize];
            }
            let mut bits: [i16; MAX_SERIAL_SIZE] = core::array::from_fn(|_| (g.u64() & 1) as i16);
            let mut r = [0i16; L_FRAME];
            dec.decode(ft, mode, &bits, &mut r);
            let m = if ft == NoData { prev_mode } else { mode };
            prev_mode = m;
            let mut c = [0i16; L_FRAME];
            unsafe {
                ref_Speech_Decode_Frame(st, m as i32, bits.as_mut_ptr(), ft as i32, c.as_mut_ptr())
            };
            assert_eq!(r, c, "seed {seed} frame {f}: {ft:?} {mode:?}");
        }
        unsafe { ref_Speech_Decode_Frame_exit(&mut st) };
    }

    /// DTX muting with a low comfort-noise energy: a SID_UPDATE (random bits: log_en index
    /// <= 10 gives log_en <= 0), then NO_DATA past `DTX_MAX_EMPTY_THRESH` (mute, log_en -=
    /// 256: the reference goes negative, opencore and this port once clamped at 0). Found by
    /// `nb_decode_vs_ref` after 30 min.
    #[test]
    fn decoder_matches_ref_through_dtx_mute() {
        let _l = ref26073_sys::lock();
        let mut g = Rng::new(61);
        for seed in 0..20 {
            let mut dec = Decoder::new();
            let mut st = null_mut();
            unsafe { ref_Speech_Decode_Frame_init(&mut st, c"Decoder".as_ptr().cast()) };
            let mode = Mode::ALL[(g.u64() % 8) as usize];
            let mut sid: [i16; MAX_SERIAL_SIZE] = core::array::from_fn(|_| (g.u64() & 1) as i16);
            for f in 0..130 {
                let ft = if f == 0 {
                    RxFrameType::SidUpdate
                } else {
                    RxFrameType::NoData
                };
                let mut r = [0i16; L_FRAME];
                dec.decode(ft, mode, &sid, &mut r);
                let mut c = [0i16; L_FRAME];
                unsafe {
                    ref_Speech_Decode_Frame(
                        st,
                        mode as i32,
                        sid.as_mut_ptr(),
                        ft as i32,
                        c.as_mut_ptr(),
                    )
                };
                assert_eq!(r, c, "seed {seed} frame {f}: {ft:?} {mode:?}");
            }
            unsafe { ref_Speech_Decode_Frame_exit(&mut st) };
        }
    }

    #[test]
    fn decoder_matches_ref_on_random_frames() {
        for seed in 0..20 {
            run(seed, 2000);
        }
    }
}
