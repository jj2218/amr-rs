//! Decoder DTX / comfort noise: `dtx_dec.cpp` (state machine + activity update; `dtx_dec`
//! itself lands with the decoder core), `b_cn_cod.cpp`, `a_refl.cpp`.
//! Bit-exact port of opencore-amr 0.1.6 (Apache-2.0).

use crate::basic_op::{
    Flag, MAX_32, abs_s, add, div_s, l_add, l_add_energy, l_msu, l_mult, l_shl, l_shr_r, l_sub,
    norm_l, round, shl,
};
use crate::math::log2;

use crate::cnst::{L_FRAME, L_SUBFR, M};
pub(crate) const DTX_HIST_SIZE: usize = 8;
const DTX_HANG_CONST: i16 = 7;
const DTX_MAX_EMPTY_THRESH: i16 = 50;
const DTX_ELAPSED_FRAMES_THRESH: i16 = 24 + 7 - 1;
const PN_INITIAL_SEED: i32 = 0x7081_6958;

/// Receive frame type (TS 26.073 `enum RXFrameType`): how the channel delivered the frame.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[repr(i32)]
pub enum RxFrameType {
    /// Speech frame received correctly.
    SpeechGood = 0,
    /// Speech frame with possible bit errors.
    SpeechDegraded,
    /// Onset: first speech frame after DTX (no speech parameters).
    Onset,
    /// Speech frame with errors (concealed).
    SpeechBad,
    /// First SID frame of a DTX period.
    SidFirst,
    /// SID update with comfort-noise parameters.
    SidUpdate,
    /// Corrupted SID frame.
    SidBad,
    /// Nothing received (DTX pause or frame loss).
    NoData,
}

/// `enum DTXStateType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
#[repr(i32)]
pub(crate) enum DtxState {
    #[default]
    Speech = 0,
    Dtx,
    DtxMute,
}

/// `dtx_decState`.
#[derive(Clone, Debug)]
pub(crate) struct DtxDecState {
    pub since_last_sid: i16,
    pub true_sid_period_inv: i16,
    pub log_en: i16,
    pub old_log_en: i16,
    pub l_pn_seed_rx: i32,
    pub lsp: [i16; M],
    pub lsp_old: [i16; M],
    pub lsf_hist: [i16; M * DTX_HIST_SIZE],
    pub lsf_hist_ptr: i16,
    pub lsf_hist_mean: [i16; M * DTX_HIST_SIZE],
    pub log_pg_mean: i16,
    pub log_en_hist: [i16; DTX_HIST_SIZE],
    pub log_en_hist_ptr: i16,
    pub log_en_adjust: i16,
    pub dtx_hangover_count: i16,
    pub dec_ana_elapsed_count: i16,
    pub sid_frame: i16,
    pub valid_data: i16,
    pub dtx_hangover_added: i16,
    pub dtx_global_state: DtxState,
    pub data_updated: i16,
}

const LSP_INIT: [i16; M] = [
    30000, 26000, 21000, 15000, 8000, 0, -8000, -15000, -21000, -26000,
];

impl DtxDecState {
    /// A freshly reset state.
    pub fn new() -> Self {
        let mut s = Self {
            since_last_sid: 0,
            true_sid_period_inv: 0,
            log_en: 0,
            old_log_en: 0,
            l_pn_seed_rx: 0,
            lsp: [0; M],
            lsp_old: [0; M],
            lsf_hist: [0; M * DTX_HIST_SIZE],
            lsf_hist_ptr: 0,
            lsf_hist_mean: [0; M * DTX_HIST_SIZE],
            log_pg_mean: 0,
            log_en_hist: [0; DTX_HIST_SIZE],
            log_en_hist_ptr: 0,
            log_en_adjust: 0,
            dtx_hangover_count: 0,
            dec_ana_elapsed_count: 0,
            sid_frame: 0,
            valid_data: 0,
            dtx_hangover_added: 0,
            dtx_global_state: DtxState::Speech,
            data_updated: 0,
        };
        s.reset();
        s
    }

    /// `dtx_dec_reset`
    pub fn reset(&mut self) {
        self.since_last_sid = 0;
        self.true_sid_period_inv = 1 << 13;
        self.log_en = 3500;
        self.old_log_en = 3500;
        self.l_pn_seed_rx = PN_INITIAL_SEED;
        self.lsp = LSP_INIT;
        self.lsp_old = LSP_INIT;
        self.lsf_hist_ptr = 0;
        self.log_pg_mean = 0;
        self.log_en_hist_ptr = 0;
        let first = [
            1384, 2077, 3420, 5108, 6742, 8122, 9863, 11092, 12714, 13701,
        ];
        for h in self.lsf_hist.as_chunks_mut::<M>().0 {
            h.copy_from_slice(&first);
        }
        self.lsf_hist_mean = [0; M * DTX_HIST_SIZE];
        self.log_en_hist = [self.log_en; DTX_HIST_SIZE];
        self.log_en_adjust = 0;
        self.dtx_hangover_count = DTX_HANG_CONST;
        self.dec_ana_elapsed_count = 32767;
        self.sid_frame = 0;
        self.valid_data = 0;
        self.dtx_hangover_added = 0;
        self.dtx_global_state = DtxState::Dtx;
        self.data_updated = 0;
    }

    /// `dtx_dec_activity_update`
    pub fn activity_update(&mut self, lsf: &[i16; M], frame: &[i16; L_FRAME], ov: &mut Flag) {
        self.lsf_hist_ptr += M as i16;
        if self.lsf_hist_ptr == 80 {
            self.lsf_hist_ptr = 0;
        }
        let p = self.lsf_hist_ptr as usize;
        self.lsf_hist[p..p + M].copy_from_slice(lsf);

        let l_frame_en = l_add_energy(frame, ov);
        let (mut log_en_e, mut log_en_m) = log2(l_frame_en, ov);
        let t = (log_en_e as i32) << 10;
        log_en_e = if t != t as i16 as i32 {
            *ov = true;
            if log_en_e > 0 { i16::MAX } else { i16::MIN }
        } else {
            t as i16
        };
        log_en_m >>= 5; // arithmetic shift == the C `~((~x) >> 5)` for negatives
        let log_en = log_en_e.wrapping_add(log_en_m).wrapping_sub(7497 + 1024);

        self.log_en_hist_ptr += 1;
        if self.log_en_hist_ptr == DTX_HIST_SIZE as i16 {
            self.log_en_hist_ptr = 0;
        }
        self.log_en_hist[self.log_en_hist_ptr as usize] = log_en;
    }

    /// `rx_dtx_handler`
    pub fn rx_dtx_handler(&mut self, frame_type: RxFrameType, ov: &mut Flag) -> DtxState {
        use RxFrameType as R;
        let in_dtx = matches!(self.dtx_global_state, DtxState::Dtx | DtxState::DtxMute);
        let new_state = if matches!(frame_type, R::SidFirst | R::SidUpdate | R::SidBad)
            || in_dtx && matches!(frame_type, R::NoData | R::SpeechBad | R::Onset)
        {
            let mut s = DtxState::Dtx;
            if self.dtx_global_state == DtxState::DtxMute
                && matches!(frame_type, R::SidBad | R::SidFirst | R::Onset | R::NoData)
            {
                s = DtxState::DtxMute;
            }
            self.since_last_sid += 1;
            if frame_type != R::SidUpdate && self.since_last_sid > DTX_MAX_EMPTY_THRESH {
                s = DtxState::DtxMute;
            }
            s
        } else {
            self.since_last_sid = 0;
            DtxState::Speech
        };

        if self.data_updated == 0 && frame_type == R::SidUpdate {
            self.dec_ana_elapsed_count = 0;
        }
        self.dec_ana_elapsed_count = add(self.dec_ana_elapsed_count, 1, ov);
        self.dtx_hangover_added = 0;

        let enc_state = if matches!(
            frame_type,
            R::SidFirst | R::SidUpdate | R::SidBad | R::Onset | R::NoData
        ) {
            if frame_type == R::NoData && new_state == DtxState::Speech {
                DtxState::Speech
            } else {
                DtxState::Dtx
            }
        } else {
            DtxState::Speech
        };

        if enc_state == DtxState::Speech {
            self.dtx_hangover_count = DTX_HANG_CONST;
        } else if self.dec_ana_elapsed_count > DTX_ELAPSED_FRAMES_THRESH {
            self.dtx_hangover_added = 1;
            self.dec_ana_elapsed_count = 0;
            self.dtx_hangover_count = 0;
        } else if self.dtx_hangover_count == 0 {
            self.dec_ana_elapsed_count = 0;
        } else {
            self.dtx_hangover_count -= 1;
        }

        if new_state != DtxState::Speech {
            self.sid_frame = 0;
            self.valid_data = 0;
            match frame_type {
                R::SidFirst => self.sid_frame = 1,
                R::SidUpdate => {
                    self.sid_frame = 1;
                    self.valid_data = 1;
                }
                R::SidBad => {
                    self.sid_frame = 1;
                    self.dtx_hangover_added = 0;
                }
                _ => {}
            }
        }
        new_state
    }
}

const LSF_HIST_MEAN_SCALE: [i16; M] = [20000, 20000, 20000, 20000, 20000, 18000, 16384, 8192, 0, 0];
const DTX_LOG_EN_ADJUST: [i16; 9] = [-1023, -878, -732, -586, -440, -294, -148, 0, 0];
const LSF_GAP: i16 = 205;

/// Decoder state owned by other groups that `dtx_dec` reads and updates.
pub(crate) struct DtxDecDeps<'a> {
    pub mem_syn: &'a mut [i16; M],
    pub lsf_state: &'a mut crate::lpc::DPlsfState,
    pub pred_state: &'a mut crate::gains::GcPredState,
    pub aver_state: &'a mut crate::gains::CbGainAverageState,
}

/// `dtx_dec`: comfort-noise synthesis for one frame (SID / no-data / mute).
/// Covered end to end by the 3GPP DTX/DTX2 conformance vectors.
#[allow(clippy::too_many_arguments)]
pub(crate) fn dtx_dec(
    st: &mut DtxDecState,
    deps: DtxDecDeps<'_>,
    new_state: DtxState,
    mode: crate::cnst::CodecMode,
    parm: &[i16],
    synth: &mut [i16; L_FRAME],
    a_t: &mut [i16; 4 * (M + 1)],
    ov: &mut Flag,
) {
    use crate::basic_op::{MAX_16, MIN_16, l_mac, l_shr, mult, shr, sub};
    use crate::cnst::CodecMode;
    use crate::lpc::{d_plsf_3, init_d_plsf_3, lsf_lsp, lsp_az, lsp_lsf, reorder_lsf};
    let DtxDecDeps {
        mem_syn,
        lsf_state,
        pred_state,
        aver_state,
    } = deps;
    let m = mode as usize;

    if st.dtx_hangover_added != 0 && st.sid_frame != 0 {
        st.log_en_adjust = DTX_LOG_EN_ADJUST[m];
        let mut ptr = st.lsf_hist_ptr as usize + M;
        if ptr == 80 {
            ptr = 0;
        }
        let cur = st.lsf_hist_ptr as usize;
        st.lsf_hist.copy_within(cur..cur + M, ptr);
        let mut ptr = st.log_en_hist_ptr as usize + 1;
        if ptr == DTX_HIST_SIZE {
            ptr = 0;
        }
        st.log_en_hist[ptr] = st.log_en_hist[st.log_en_hist_ptr as usize];

        st.log_en = 0;
        let mut l_lsf = [0i32; M];
        for i in (0..DTX_HIST_SIZE).rev() {
            st.log_en = add(st.log_en, st.log_en_hist[i] >> 3, ov);
            for j in (0..M).rev() {
                l_lsf[j] = l_add(l_lsf[j], st.lsf_hist[i * M + j] as i32, ov);
            }
        }
        let mut lsf = [0i16; M];
        for j in (0..M).rev() {
            lsf[j] = (l_lsf[j] >> 3) as i16;
        }
        lsf_lsp(&lsf, &mut st.lsp, M, ov);
        st.log_en = sub(st.log_en, st.log_en_adjust, ov);

        st.lsf_hist_mean = st.lsf_hist;
        for i in (0..M).rev() {
            let mut l_mean: i32 = 0;
            for j in (0..8).rev() {
                l_mean = l_add(l_mean, st.lsf_hist_mean[i + j * M] as i32, ov);
            }
            let lsf_mean = (l_mean >> 3) as i16;
            for j in (0..8).rev() {
                let k = i + j * M;
                let mut v = sub(st.lsf_hist_mean[k], lsf_mean, ov);
                v = mult(v, LSF_HIST_MEAN_SCALE[i], ov);
                let negative = v < 0;
                v = abs_s(v);
                if v > 655 {
                    v = 655 + ((v - 655) >> 2);
                }
                if v > 1310 {
                    v = 1310;
                }
                st.lsf_hist_mean[k] = if negative { -v } else { v };
            }
        }
    }

    if st.sid_frame != 0 {
        st.lsp_old = st.lsp;
        st.old_log_en = st.log_en;
        if st.valid_data != 0 {
            let tmp_int_length = st.since_last_sid.min(32);
            st.since_last_sid = 0;
            let t = (tmp_int_length as i32) << 10;
            let temp = if t != t as i16 as i32 {
                *ov = true;
                if tmp_int_length > 0 { MAX_16 } else { MIN_16 }
            } else {
                t as i16
            };
            st.true_sid_period_inv = if tmp_int_length >= 2 {
                div_s(1 << 10, temp)
            } else {
                1 << 14
            };

            init_d_plsf_3(lsf_state, parm[0]);
            d_plsf_3(lsf_state, CodecMode::MRDTX, 0, &parm[1..], &mut st.lsp, ov);
            lsf_state.past_r_q = [0; M];

            let log_en_index = parm[4];
            st.log_en = if !(-64..=63).contains(&log_en_index) {
                if log_en_index > 0 { MAX_16 } else { MIN_16 }
            } else {
                log_en_index << (11 - 2)
            };
            st.log_en = sub(st.log_en, 2560 * 2, ov);
            if log_en_index == 0 {
                st.log_en = MIN_16;
            }
            if st.data_updated == 0 || st.dtx_global_state == DtxState::Speech {
                st.lsp_old = st.lsp;
                st.old_log_en = st.log_en;
            }
        }
        let ma_pred_init = ((st.log_en >> 1).wrapping_sub(9000)).clamp(-14436, 0);
        pred_state.past_qua_en = [ma_pred_init; 4];
        let mr122 = ((ma_pred_init as i32 * 5443) >> 15) as i16;
        pred_state.past_qua_en_mr122 = [mr122; 4];
    }

    let adj = DTX_LOG_EN_ADJUST[m];
    let scaled = if adj > 1023 {
        MAX_16
    } else if adj < -1024 {
        MIN_16
    } else {
        ((((adj as i32) << 5) * 3277) >> 15) as i16
    };
    let temp = scaled >> 5;
    st.log_en_adjust = add(((st.log_en_adjust as i32 * 29491) >> 15) as i16, temp, ov);

    let mut int_fac = shl(st.since_last_sid.wrapping_add(1), 10, ov);
    int_fac = mult(int_fac, st.true_sid_period_inv, ov);
    int_fac = if int_fac > 1024 {
        16384
    } else if int_fac < -2048 {
        MIN_16
    } else {
        int_fac << 4
    };
    let mut l_log_en_int = l_mult(int_fac, st.log_en, ov);
    let mut lsp_int = [0i16; M];
    for i in (0..M).rev() {
        lsp_int[i] = mult(int_fac, st.lsp[i], ov);
    }
    int_fac = 16384i16.wrapping_sub(int_fac);
    l_log_en_int = l_mac(l_log_en_int, int_fac, st.old_log_en, ov);
    for i in (0..M).rev() {
        lsp_int[i] = add(lsp_int[i], mult(int_fac, st.lsp_old[i], ov), ov);
        let t = (lsp_int[i] as i32) << 1;
        lsp_int[i] = if t != t as i16 as i32 {
            *ov = true;
            if lsp_int[i] > 0 { MAX_16 } else { MIN_16 }
        } else {
            t as i16
        };
    }

    let mut lsf_variab_factor = st.log_pg_mean.wrapping_sub(2457);
    lsf_variab_factor = 4096i16.wrapping_sub(mult(lsf_variab_factor, 9830, ov));
    lsf_variab_factor = if lsf_variab_factor > 4095 {
        MAX_16
    } else if lsf_variab_factor < 0 {
        0
    } else {
        lsf_variab_factor << 3
    };
    let lsf_variab_index = pseudonoise(&mut st.l_pn_seed_rx, 3) as usize;

    let mut lsf_int = [0i16; M];
    lsp_lsf(&lsp_int, &mut lsf_int, M, ov);
    let mut lsf_int_variab = lsf_int;
    for i in (0..M).rev() {
        let d = mult(
            lsf_variab_factor,
            st.lsf_hist_mean[i + lsf_variab_index * M],
            ov,
        );
        lsf_int_variab[i] = add(lsf_int_variab[i], d, ov);
    }
    reorder_lsf(&mut lsf_int, LSF_GAP, M, ov);
    reorder_lsf(&mut lsf_int_variab, LSF_GAP, M, ov);
    lsf_state.past_lsf_q = lsf_int;
    let mut lsp_int_variab = [0i16; M];
    lsf_lsp(&lsf_int, &mut lsp_int, M, ov);
    lsf_lsp(&lsf_int_variab, &mut lsp_int_variab, M, ov);
    let (mut acoeff, mut acoeff_variab) = ([0i16; M + 1], [0i16; M + 1]);
    lsp_az(&lsp_int, &mut acoeff, ov);
    lsp_az(&lsp_int_variab, &mut acoeff_variab, ov);
    for sub_az in a_t.as_chunks_mut::<{ M + 1 }>().0 {
        sub_az.copy_from_slice(&acoeff);
    }

    let mut refl = [0i16; M];
    a_refl(acoeff[1..].try_into().unwrap(), &mut refl, ov);
    let mut pred_err = MAX_16;
    for &r in &refl {
        let t = (r as i32 * r as i32) >> 15;
        let temp = if t <= 0x7fff {
            MAX_16 - t as i16
        } else {
            *ov = true;
            0
        };
        pred_err = mult(pred_err, temp, ov);
    }
    let (log_pg_e, log_pg_m) = log2(pred_err as i32, ov);
    let mut log_pg = shl(log_pg_e.wrapping_sub(15), 12, ov);
    log_pg = shr(
        sub(0, add(log_pg, shr(log_pg_m, 15 - 12, ov), ov), ov),
        1,
        ov,
    );
    st.log_pg_mean = add(mult(29491, st.log_pg_mean, ov), mult(3277, log_pg, ov), ov);

    l_log_en_int = l_shr(l_log_en_int, 10, ov);
    l_log_en_int = l_add(l_log_en_int, 4 * 65536, ov);
    l_log_en_int = l_sub(l_log_en_int, l_shl(log_pg as i32, 4, ov), ov);
    l_log_en_int = l_add(l_log_en_int, l_shl(st.log_en_adjust as i32, 5, ov), ov);
    let log_en_int_e = (l_log_en_int >> 16) as i16;
    let log_en_int_m = l_shr(l_sub(l_log_en_int, (log_en_int_e as i32) << 16, ov), 1, ov) as i16;
    let level = crate::math::pow2(log_en_int_e, log_en_int_m, ov) as i16;

    let mut ex = [0i16; L_SUBFR];
    for i in 0..4 {
        build_cn_code(&mut st.l_pn_seed_rx, &mut ex, ov);
        for e in ex.iter_mut().rev() {
            *e = mult(level, *e, ov);
        }
        let out = &mut synth[i * L_SUBFR..(i + 1) * L_SUBFR];
        crate::post::filt::syn_filt(&acoeff_variab, &ex, out, L_SUBFR, mem_syn, true);
    }

    aver_state.hang_var = 20;
    aver_state.hang_count = 0;

    if new_state == DtxState::DtxMute {
        let tmp_int_length = match st.since_last_sid {
            n if n > 32 => 32,
            n if n <= 0 => 8,
            n => n,
        };
        let t = (tmp_int_length as i32) << 10;
        let temp = if t != t as i16 as i32 {
            *ov = true;
            if tmp_int_length > 0 { MAX_16 } else { MIN_16 }
        } else {
            t as i16
        };
        st.true_sid_period_inv = div_s(1 << 10, temp);
        st.since_last_sid = 0;
        st.lsp_old = st.lsp;
        st.old_log_en = st.log_en;
        // subtract 1/8 in Q11 (-6/8 dB): TS 26.073 saturates and lets it go negative
        // (opencore clamps at 0: a decoder that mutes twice with a low CN energy then differs)
        st.log_en = sub(st.log_en, 256, ov);
    }

    if st.sid_frame != 0 && (st.valid_data != 0 || st.dtx_hangover_added != 0) {
        st.since_last_sid = 0;
        st.data_updated = 1;
    }
}

/// `pseudonoise`: `no_bits` bits from the 31-bit LFSR in `shift_reg`.
pub(crate) fn pseudonoise(shift_reg: &mut i32, no_bits: i16) -> i16 {
    let mut noise_bits: i16 = 0;
    for _ in 0..no_bits {
        let sn = (*shift_reg & 1) ^ ((*shift_reg >> 28) & 1);
        noise_bits = (noise_bits << 1) | (*shift_reg & 1) as i16;
        *shift_reg >>= 1;
        if sn & 1 != 0 {
            *shift_reg |= 0x4000_0000;
        }
    }
    noise_bits
}

/// `build_CN_code`: random 10-pulse innovation for comfort noise.
pub(crate) fn build_cn_code(seed: &mut i32, cod: &mut [i16; L_SUBFR], ov: &mut Flag) {
    cod.fill(0);
    for k in 0..10i16 {
        let i = pseudonoise(seed, 2);
        let temp = l_mult(i, 10, ov) as i16;
        let i = add(temp >> 1, k, ov) as usize;
        let j = pseudonoise(seed, 1);
        cod[i] = if j > 0 { 4096 } else { -4096 };
    }
}

/// `build_CN_param`: random parameters for a lost/no-data frame.
pub(crate) fn build_cn_param(
    seed: &mut i16,
    n_param: usize,
    param_size_table: &[i16],
    parm: &mut [i16],
    window_200_40: &[i16],
    ov: &mut Flag,
) {
    let t = l_mult(*seed, 31821, ov) >> 1;
    *seed = l_add(t, 13849, ov) as i16;
    let start = (*seed & 0x7f) as usize;
    for i in 0..n_param {
        let mask = !(0xFFFFi32 << param_size_table[i]) as i16;
        parm[i] = window_200_40[start + i] & mask;
    }
}

/// `A_Refl`: LPC coefficients `a` to reflection coefficients `refl` (zeroed if unstable).
pub(crate) fn a_refl(a: &[i16; M], refl: &mut [i16; M], ov: &mut Flag) {
    a_refl_impl(a, refl, ov);
}

/// Returns `true` if the inner loop bailed out early. In C that path then copies elements of the
/// uninitialised local `bState` into `aState` (undefined behaviour; confirmed with valgrind);
/// here `b_state` is zero-initialised, so the result is deterministic.
fn a_refl_impl(a: &[i16; M], refl: &mut [i16; M], ov: &mut Flag) -> bool {
    let mut undefined_in_c = false;
    let mut a_state = *a;
    let mut b_state = [0i16; M];
    for i in (0..M).rev() {
        if abs_s(a_state[i]) >= 4096 {
            refl.fill(0);
            break;
        }
        refl[i] = shl(a_state[i], 3, ov);
        let t = l_mult(refl[i], refl[i], ov);
        let mut acc = l_sub(MAX_32, t, ov);
        let norm_shift = norm_l(acc);
        let scale = 15 - norm_shift;
        acc = l_shl(acc, norm_shift, ov);
        let norm_prod = round(acc, ov);
        let mult = div_s(16384, norm_prod);
        for j in 0..i {
            let mut acc = (a_state[j] as i32) << 16;
            acc = l_msu(acc, refl[i], a_state[i - j - 1], ov);
            let temp = round(acc, ov);
            let mut t = l_mult(mult, temp, ov);
            t = l_shr_r(t, scale, ov);
            let abs = {
                let y = t.wrapping_sub((t < 0) as i32);
                y ^ (y >> 31)
            };
            if abs > 32767 {
                refl.fill(0);
                undefined_in_c = true;
                break; // C quirk: leaves only the inner loop; the outer loop continues
            }
            b_state[j] = t as i16;
        }
        a_state[..i].copy_from_slice(&b_state[..i]);
    }
    undefined_in_c
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::testutil::{N, Rng};
    use opencore_sys as _;

    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq)]
    struct CDtxDecState {
        since_last_sid: i16,
        true_sid_period_inv: i16,
        log_en: i16,
        old_log_en: i16,
        l_pn_seed_rx: i32,
        lsp: [i16; M],
        lsp_old: [i16; M],
        lsf_hist: [i16; M * DTX_HIST_SIZE],
        lsf_hist_ptr: i16,
        lsf_hist_mean: [i16; M * DTX_HIST_SIZE],
        log_pg_mean: i16,
        log_en_hist: [i16; DTX_HIST_SIZE],
        log_en_hist_ptr: i16,
        log_en_adjust: i16,
        dtx_hangover_count: i16,
        dec_ana_elapsed_count: i16,
        sid_frame: i16,
        valid_data: i16,
        dtx_hangover_added: i16,
        dtx_global_state: i32,
        data_updated: i16,
    }

    impl From<&DtxDecState> for CDtxDecState {
        fn from(s: &DtxDecState) -> Self {
            Self {
                since_last_sid: s.since_last_sid,
                true_sid_period_inv: s.true_sid_period_inv,
                log_en: s.log_en,
                old_log_en: s.old_log_en,
                l_pn_seed_rx: s.l_pn_seed_rx,
                lsp: s.lsp,
                lsp_old: s.lsp_old,
                lsf_hist: s.lsf_hist,
                lsf_hist_ptr: s.lsf_hist_ptr,
                lsf_hist_mean: s.lsf_hist_mean,
                log_pg_mean: s.log_pg_mean,
                log_en_hist: s.log_en_hist,
                log_en_hist_ptr: s.log_en_hist_ptr,
                log_en_adjust: s.log_en_adjust,
                dtx_hangover_count: s.dtx_hangover_count,
                dec_ana_elapsed_count: s.dec_ana_elapsed_count,
                sid_frame: s.sid_frame,
                valid_data: s.valid_data,
                dtx_hangover_added: s.dtx_hangover_added,
                dtx_global_state: s.dtx_global_state as i32,
                data_updated: s.data_updated,
            }
        }
    }

    // The TS 26.073 reference, in the shape of the opencore functions these tests were written
    // against (`o`: its overflow flag, sticky; table arguments unused).
    unsafe extern "C" {
        fn ref_dtx_dec_reset(st: *mut CDtxDecState) -> i32;
        fn ref_dtx_dec_activity_update(st: *mut CDtxDecState, lsf: *mut i16, frame: *mut i16);
        fn ref_rx_dtx_handler(st: *mut CDtxDecState, frame_type: i32) -> i32;
        fn ref_pseudonoise(reg: *mut i32, n: i16) -> i16;
        fn ref_build_CN_code(seed: *mut i32, cod: *mut i16);
        fn ref_build_CN_param(seed: *mut i16, n: i16, sizes: *const i16, parm: *mut i16);
        fn ref_A_Refl(a: *mut i16, refl: *mut i16);
    }
    use crate::testutil::with_ref;
    const NO: *mut i32 = std::ptr::null_mut();
    unsafe fn dtx_dec_reset(st: *mut CDtxDecState) -> i16 {
        with_ref(NO, || unsafe { ref_dtx_dec_reset(st) }) as i16
    }
    unsafe fn dtx_dec_activity_update(
        st: *mut CDtxDecState,
        lsf: *mut i16,
        frame: *mut i16,
        o: *mut i32,
    ) {
        with_ref(o, || unsafe { ref_dtx_dec_activity_update(st, lsf, frame) })
    }
    unsafe fn rx_dtx_handler(st: *mut CDtxDecState, frame_type: i32, o: *mut i32) -> i32 {
        with_ref(o, || unsafe { ref_rx_dtx_handler(st, frame_type) })
    }
    unsafe fn c_pseudonoise(reg: *mut i32, n: i16) -> i16 {
        with_ref(NO, || unsafe { ref_pseudonoise(reg, n) })
    }
    #[allow(non_snake_case)]
    unsafe fn build_CN_code(seed: *mut i32, cod: *mut i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_build_CN_code(seed, cod) })
    }
    #[allow(non_snake_case)]
    unsafe fn build_CN_param(
        seed: *mut i16,
        n: i16,
        sizes: *const i16,
        parm: *mut i16,
        _win: *const i16,
        o: *mut i32,
    ) {
        with_ref(o, || unsafe { ref_build_CN_param(seed, n, sizes, parm) })
    }
    #[allow(non_snake_case)]
    unsafe fn A_Refl(a: *mut i16, refl: *mut i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_A_Refl(a, refl) })
    }

    const TYPES: [RxFrameType; 8] = [
        RxFrameType::SpeechGood,
        RxFrameType::SpeechDegraded,
        RxFrameType::Onset,
        RxFrameType::SpeechBad,
        RxFrameType::SidFirst,
        RxFrameType::SidUpdate,
        RxFrameType::SidBad,
        RxFrameType::NoData,
    ];

    #[test]
    fn dtx_state_machine_vs_ref() {
        let mut g = Rng::new(7);
        for run in 0..200 {
            let mut r = DtxDecState::new();
            let mut c: CDtxDecState = (&r).into();
            unsafe { dtx_dec_reset(&mut c) };
            assert_eq!(CDtxDecState::from(&r), c, "reset");
            for step in 0..1000 {
                // Bursty frame-type sequences so DTX/mute/hangover paths are all reached.
                let ft = if g.u64().is_multiple_of(4) {
                    TYPES[(g.u64() % 8) as usize]
                } else {
                    TYPES[(run + step / 37) % 8]
                };
                if g.u64().is_multiple_of(3) {
                    r.data_updated = (g.u64() % 2) as i16;
                    c.data_updated = r.data_updated;
                }
                let (mut ro, mut co) = (false, 0);
                let rs = r.rx_dtx_handler(ft, &mut ro);
                let cs = unsafe { rx_dtx_handler(&mut c, ft as i32, &mut co) };
                r.dtx_global_state = rs;
                c.dtx_global_state = cs;
                assert_eq!(
                    (rs as i32, ro),
                    (cs, co != 0),
                    "rx_dtx_handler run {run} step {step}"
                );

                let mut lsf = [0i16; M];
                lsf.iter_mut().for_each(|x| *x = g.i16() & 0x3fff);
                let mut frame = [0i16; L_FRAME];
                let amp = 1i16 << (g.u64() % 16);
                frame.iter_mut().for_each(|x| *x = g.i16() % amp.max(1));
                let (mut ro, mut co) = (false, 0);
                r.activity_update(&lsf, &frame, &mut ro);
                unsafe {
                    dtx_dec_activity_update(&mut c, lsf.as_mut_ptr(), frame.as_mut_ptr(), &mut co)
                };
                assert_eq!(ro, co != 0, "activity_update overflow");
                assert_eq!(CDtxDecState::from(&r), c, "state run {run} step {step}");
            }
        }
    }

    #[test]
    fn comfort_noise_builders_vs_ref() {
        let mut g = Rng::new(11);
        let sizes: &[i16] = &crate::tables::BITNO_MR122;
        for _ in 0..N / 10 {
            let mut rs = g.i32();
            let mut cs = rs;
            let n = (g.u64() % 12) as i16;
            assert_eq!(pseudonoise(&mut rs, n), unsafe {
                c_pseudonoise(&mut cs, n)
            });
            assert_eq!(rs, cs);

            let (mut rc, mut cc) = ([0i16; L_SUBFR], [0i16; L_SUBFR]);
            let (mut ro, mut co) = (false, 0);
            build_cn_code(&mut rs, &mut rc, &mut ro);
            unsafe { build_CN_code(&mut cs, cc.as_mut_ptr(), &mut co) };
            assert_eq!((rc, rs, ro), (cc, cs, co != 0));

            let mut rseed = g.i16();
            let mut cseed = rseed;
            let n = (g.u64() % (sizes.len() as u64 + 1)) as usize;
            let (mut rp, mut cp) = ([0i16; 57], [0i16; 57]);
            let (mut ro, mut co) = (false, 0);
            build_cn_param(
                &mut rseed,
                n,
                sizes,
                &mut rp,
                &crate::tables::WINDOW_200_40,
                &mut ro,
            );
            unsafe {
                build_CN_param(
                    &mut cseed,
                    n as i16,
                    sizes.as_ptr(),
                    cp.as_mut_ptr(),
                    crate::tables::WINDOW_200_40.as_ptr(),
                    &mut co,
                )
            };
            assert_eq!((rp, rseed, ro), (cp, cseed, co != 0));
        }
    }

    /// C's `A_Refl` reads its uninitialised `bState` stack array when the inner loop bails out
    /// early (valgrind: "conditional jump depends on uninitialised value"), so C is undefined
    /// for such inputs. Compare on the inputs the decoder produces (LPC from valid LSPs via
    /// `Lsp_Az`, as in `dtx_dec`) and on the deterministic `|a[i]| >= 4096` rejection path.
    #[test]
    fn a_refl_vs_ref() {
        let mut g = Rng::new(13);
        let mut undefined = 0;
        for case in 0..N {
            let mut a = [0i16; M];
            if case % 8 == 0 {
                a.iter_mut().for_each(|x| *x = g.i16() % 4096);
                a[M - 1] = if g.u64() & 1 == 0 {
                    4096
                } else {
                    -4096 - (g.i16() & 0xfff)
                };
            } else {
                // Random LSPs, descending in the cosine domain with a minimum spacing.
                let lsp = loop {
                    let mut l = [0i16; M];
                    l.iter_mut()
                        .for_each(|x| *x = ((g.u64() % 64000) as i32 - 32000) as i16);
                    l.sort_unstable_by(|a, b| b.cmp(a));
                    if l.windows(2).all(|w| w[0] as i32 - w[1] as i32 >= 300) {
                        break l;
                    }
                };
                let mut az = [0i16; M + 1];
                let mut o = false;
                crate::lpc::lsp_az(&lsp, &mut az, &mut o);
                a.copy_from_slice(&az[1..]);
            }
            let (mut rr, mut cr) = ([0i16; M], [0i16; M]);
            let (mut ro, mut co) = (false, 0);
            let mut ac = a;
            if a_refl_impl(&a, &mut rr, &mut ro) {
                undefined += 1; // C reads uninitialised memory here: nothing to compare against
                continue;
            }
            unsafe { A_Refl(ac.as_mut_ptr(), cr.as_mut_ptr(), &mut co) };
            assert_eq!((rr, ro), (cr, co != 0), "A_Refl({a:?})");
        }
        eprintln!("A_Refl: {undefined}/{N} inputs hit C's undefined path (skipped)");
        assert!(undefined < N / 10, "too few comparable cases");
    }
}
