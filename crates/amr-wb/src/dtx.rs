//! Discontinuous transmission, decoder side (`dtx_decoder_amr_wb.cpp`): the RX DTX state
//! machine and comfort-noise generation from SID parameters.
//!
//! TS 26.173 as of v7.0.0 plus CR 26.173-0030: CR 0026 (March 2007: late SID_UPDATE frames
//! don't mute, a division-by-zero guard) and CR 0030 (2015: the energy history starts at
//! 3500 >> 3), both newer than opencore 0.1.6. NO_DATA handling stays opencore's, before
//! CR 26.173-0027..29 (September 2007): every NO_DATA frame counts as DTX for the hangover
//! logic. The CRs' `dtx_vad_hist` rule (a NO_DATA frame during speech is DTX only after 8
//! frames without voice activity) is deliberately absent: it is covered by a live patent.

use crate::basic_op::*;
use crate::bits::Bits;
use crate::cnst::{ISF_GAP, L_FRAME, M, RANDOM_INITSEED, RxFrameType};
use crate::excitation::noise_gen_amrwb;
use crate::lpc::disf_ns;

pub(crate) const SPEECH: i16 = 0;
pub(crate) const DTX: i16 = 1;
pub(crate) const DTX_MUTE: i16 = 2;

const DTX_MAX_EMPTY_THRESH: i16 = 50;
const DTX_HIST_SIZE: usize = 8;
const DTX_ELAPSED_FRAMES_THRESH: i16 = 24 + 7 - 1;
const DTX_HANG_CONST: i16 = 7;
const ISF_DITH_GAP: i16 = 448;
const ISF_FACTOR_LOW: i16 = 256;
const ISF_FACTOR_STEP: i16 = 2;
const GAIN_FACTOR: i16 = 75;

/// `dtx_decState`.
#[derive(Clone)]
pub(crate) struct DtxDecState {
    pub since_last_sid: i16,
    pub true_sid_period_inv: i16,
    pub log_en: i16,
    pub old_log_en: i16,
    pub isf: [i16; M],
    pub isf_old: [i16; M],
    pub cng_seed: i16,
    pub isf_hist: [i16; M * DTX_HIST_SIZE],
    pub log_en_hist: [i16; DTX_HIST_SIZE],
    pub hist_ptr: i16,
    pub dtx_hangover_count: i16,
    pub dec_ana_elapsed_count: i16,
    pub sid_frame: bool,
    pub valid_data: bool,
    pub dtx_hangover_added: bool,
    pub dtx_global_state: i16,
    pub data_updated: bool,
    pub dither_seed: i16,
    pub cn_dith: bool,
}

impl DtxDecState {
    pub(crate) fn new(isf_init: &[i16; M]) -> Self {
        let mut st = Self {
            since_last_sid: 0,
            true_sid_period_inv: 0,
            log_en: 0,
            old_log_en: 0,
            isf: [0; M],
            isf_old: [0; M],
            cng_seed: 0,
            isf_hist: [0; M * DTX_HIST_SIZE],
            log_en_hist: [0; DTX_HIST_SIZE],
            hist_ptr: 0,
            dtx_hangover_count: 0,
            dec_ana_elapsed_count: 0,
            sid_frame: false,
            valid_data: false,
            dtx_hangover_added: false,
            dtx_global_state: SPEECH,
            data_updated: false,
            dither_seed: 0,
            cn_dith: false,
        };
        st.reset(isf_init);
        st
    }

    /// `dtx_dec_amr_wb_reset`.
    pub(crate) fn reset(&mut self, isf_init: &[i16; M]) {
        self.since_last_sid = 0;
        self.true_sid_period_inv = 1 << 13;
        self.log_en = 3500;
        self.old_log_en = 3500;
        self.cng_seed = RANDOM_INITSEED;
        self.hist_ptr = 0;
        self.isf = *isf_init;
        self.isf_old = *isf_init;
        for h in self.isf_hist.as_chunks_mut::<M>().0 {
            h.copy_from_slice(isf_init);
        }
        // CR 0030: the history holds log energies scaled down by 2^3 (opencore: 3500)
        self.log_en_hist = [3500 >> 3; DTX_HIST_SIZE];
        self.dtx_hangover_count = DTX_HANG_CONST;
        self.dec_ana_elapsed_count = 32767;
        self.sid_frame = false;
        self.valid_data = false;
        self.dtx_hangover_added = false;
        self.dtx_global_state = SPEECH;
        self.data_updated = false;
        self.dither_seed = RANDOM_INITSEED;
        self.cn_dith = false;
    }

    /// `dtx_dec_amr_wb`: comfort-noise excitation `exc2` (one frame) and ISFs.
    pub(crate) fn dec(
        &mut self,
        exc2: &mut [i16; L_FRAME],
        new_state: i16,
        isf: &mut [i16; M],
        prms: &mut Bits,
    ) {
        if self.dtx_hangover_added && self.sid_frame {
            // average the history: the frame before the hangover is the reference
            let hp = self.hist_ptr as usize;
            let mut ptr = hp + 1;
            if ptr == DTX_HIST_SIZE {
                ptr = 0;
            }
            self.isf_hist.copy_within(hp * M..hp * M + M, ptr * M);
            self.log_en_hist[ptr] = self.log_en_hist[hp];
            self.log_en = 0;
            let mut l_isf = [0i32; M];
            for i in 0..DTX_HIST_SIZE {
                self.log_en = add_int16(self.log_en, self.log_en_hist[i]);
                for j in 0..M {
                    l_isf[j] = add_int32(l_isf[j], self.isf_hist[i * M + j] as i32);
                }
            }
            self.log_en = (self.log_en >> 1).wrapping_add(1024);
            if self.log_en < 0 {
                self.log_en = 0;
            }
            for j in 0..M {
                self.isf[j] = (l_isf[j] >> 3) as i16;
            }
        }
        if self.sid_frame {
            self.isf_old = self.isf;
            self.old_log_en = self.log_en;
            if self.valid_data {
                let len = self.since_last_sid.min(32);
                self.true_sid_period_inv = if len >= 2 {
                    div_16by16(1 << 10, shl_int16(len, 10))
                } else {
                    1 << 14
                };
                let mut ind = [0i16; 5];
                for (k, n) in [6, 6, 6, 5, 5].into_iter().enumerate() {
                    ind[k] = prms.parm(n);
                }
                disf_ns(&ind, &mut self.isf);
                let log_en_index = prms.parm(6);
                self.cn_dith = prms.bit() != 0;
                self.log_en = shl_int16(log_en_index, 15 - 6);
                self.log_en = mult_int16(self.log_en, 12483);
                if !self.data_updated || self.dtx_global_state == SPEECH {
                    self.isf_old = self.isf;
                    self.old_log_en = self.log_en;
                }
            }
        }
        if self.sid_frame && self.valid_data {
            self.since_last_sid = 0;
        }
        // interpolate ISFs and energy between the last two SIDs
        let mut int_fac = shl_int16(self.since_last_sid, 10);
        int_fac = mult_int16(int_fac, self.true_sid_period_inv);
        if int_fac > 1024 {
            int_fac = 1024;
        }
        int_fac = shl_int16(int_fac, 4);
        let mut l_log_en_int = mul_16by16_to_int32(int_fac, self.log_en);
        for i in 0..M {
            isf[i] = mult_int16(int_fac, self.isf[i]);
        }
        let int_fac = 16384 - int_fac;
        l_log_en_int = mac_16by16_to_int32(l_log_en_int, int_fac, self.old_log_en);
        for i in 0..M {
            isf[i] = add_int16(isf[i], mult_int16(int_fac, self.isf_old[i]));
            isf[i] = shl_int16(isf[i], 1);
        }
        if self.cn_dith {
            cn_dithering(isf, &mut l_log_en_int, &mut self.dither_seed);
        }
        l_log_en_int >>= 9;
        let mut e = extract_h(l_log_en_int);
        let m = (sub_int32(l_log_en_int, l_deposit_h(e)) >> 1) as i16;
        e = e.wrapping_add(15);
        let mut level32 = power_of_2(e, m);
        let exp0 = normalize_amr_wb(level32);
        level32 = level32.wrapping_shl(exp0 as u32);
        let exp0 = 15 - exp0;
        let level = (level32 >> 16) as i16;

        // white noise scaled to the level
        for v in exc2.iter_mut() {
            *v = noise_gen_amrwb(&mut self.cng_seed) >> 4;
        }
        let mut exp = 0i16;
        let mut ener32 = dot_product12(exc2, exc2, &mut exp);
        one_ov_sqrt_norm(&mut ener32, &mut exp);
        let gain = mult_int16(level, extract_h(ener32));
        let exp = exp.wrapping_add(exp0).wrapping_add(4);
        for v in exc2.iter_mut() {
            *v = shl_int16(mult_int16(*v, gain), exp);
        }

        if new_state == DTX_MUTE {
            let mut len = self.since_last_sid.min(32);
            if len <= 0 {
                len = 8; // CR 0026: opencore divides by zero (its div_16by16 then gives 0)
            }
            self.true_sid_period_inv = div_16by16(1 << 10, shl_int16(len, 10));
            self.since_last_sid = 0;
            self.old_log_en = self.log_en;
            // -3/8 dB per muted frame, below 0 too: TS 26.173 has no floor here (opencore
            // clamped at 0, so long mutes stayed louder)
            self.log_en = sub_int16(self.log_en, 64);
        }
        if self.sid_frame && (self.valid_data || self.dtx_hangover_added) {
            self.since_last_sid = 0;
            self.data_updated = true;
        }
    }

    /// `dtx_dec_amr_wb_activity_update`: ISF and energy history from a speech frame.
    pub(crate) fn activity_update(&mut self, isf: &[i16; M], exc: &[i16]) {
        self.hist_ptr += 1;
        if self.hist_ptr as usize == DTX_HIST_SIZE {
            self.hist_ptr = 0;
        }
        let hp = self.hist_ptr as usize;
        self.isf_hist[hp * M..hp * M + M].copy_from_slice(isf);
        let mut l = 0i32;
        for &v in &exc[..L_FRAME] {
            l = mac_16by16_to_int32(l, v, v);
        }
        l >>= 1;
        let (mut e, mut m) = (0i16, 0i16);
        amrwb_log_2(l, &mut e, &mut m);
        let log_en = shl_int16(e, 7).wrapping_add(m >> 8).wrapping_sub(1024);
        self.log_en_hist[hp] = log_en;
    }

    /// `rx_amr_wb_dtx_handler`: the new DTX state for a received frame type.
    pub(crate) fn rx_handler(&mut self, frame_type: RxFrameType) -> i16 {
        use RxFrameType::*;
        let in_dtx = self.dtx_global_state == DTX || self.dtx_global_state == DTX_MUTE;
        let new_state = if matches!(frame_type, SidFirst | SidUpdate | SidBad)
            || (in_dtx && matches!(frame_type, NoData | SpeechBad | SpeechLost))
        {
            let mut s = DTX;
            if self.dtx_global_state == DTX_MUTE
                && matches!(frame_type, SidBad | SidFirst | SpeechLost | NoData)
            {
                s = DTX_MUTE;
            }
            self.since_last_sid = add_int16(self.since_last_sid, 1);
            // CR 0026: a late SID_UPDATE (since_last_sid is reset only after it) doesn't mute
            if frame_type != SidUpdate && self.since_last_sid > DTX_MAX_EMPTY_THRESH {
                s = DTX_MUTE;
            }
            s
        } else {
            self.since_last_sid = 0;
            SPEECH
        };
        if !self.data_updated && frame_type == SidUpdate {
            self.dec_ana_elapsed_count = 0;
        }
        self.dec_ana_elapsed_count = add_int16(self.dec_ana_elapsed_count, 1);
        self.dtx_hangover_added = false;
        // opencore (pre CR 26.173-0027): every NO_DATA frame counts as DTX here
        let enc_dtx = matches!(frame_type, SidFirst | SidUpdate | SidBad | NoData);
        if !enc_dtx {
            self.dtx_hangover_count = DTX_HANG_CONST;
        } else if self.dec_ana_elapsed_count > DTX_ELAPSED_FRAMES_THRESH {
            self.dtx_hangover_added = true;
            self.dec_ana_elapsed_count = 0;
            self.dtx_hangover_count = 0;
        } else if self.dtx_hangover_count == 0 {
            self.dec_ana_elapsed_count = 0;
        } else {
            self.dtx_hangover_count -= 1;
        }
        if new_state != SPEECH {
            self.sid_frame = false;
            self.valid_data = false;
            match frame_type {
                SidFirst => self.sid_frame = true,
                SidUpdate => {
                    self.sid_frame = true;
                    self.valid_data = true;
                }
                SidBad => {
                    self.sid_frame = true;
                    self.dtx_hangover_added = false;
                }
                _ => {}
            }
        }
        new_state
    }
}

/// `CN_dithering`: dithers the comfort-noise ISFs and energy.
fn cn_dithering(isf: &mut [i16; M], l_log_en_int: &mut i32, dither_seed: &mut i16) {
    let rand2 = |seed: &mut i16| {
        let a = noise_gen_amrwb(seed) >> 1;
        let b = noise_gen_amrwb(seed) >> 1;
        a.wrapping_add(b)
    };
    let r = rand2(dither_seed);
    *l_log_en_int = add_int32(*l_log_en_int, mul_16by16_to_int32(r, GAIN_FACTOR));
    if *l_log_en_int < 0 {
        *l_log_en_int = 0;
    }
    let mut dither_fac = ISF_FACTOR_LOW;
    let r = rand2(dither_seed);
    let temp = add_int16(isf[0], mult_int16_r(r, dither_fac));
    isf[0] = if temp < ISF_GAP { ISF_GAP } else { temp };
    for i in 1..M - 1 {
        dither_fac = add_int16(dither_fac, ISF_FACTOR_STEP);
        let r = rand2(dither_seed);
        let temp = add_int16(isf[i], mult_int16_r(r, dither_fac));
        let temp1 = sub_int16(temp, isf[i - 1]);
        isf[i] = if temp1 < ISF_DITH_GAP {
            isf[i - 1].wrapping_add(ISF_DITH_GAP)
        } else {
            temp
        };
    }
    if isf[M - 2] > 16384 {
        isf[M - 2] = 16384;
    }
}
