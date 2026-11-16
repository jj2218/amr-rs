//! Voice activity detection (vo-amrwbenc `wb_vad.c`, constants from `wb_vad_c.h`), with the
//! reference's saturation: a 12-band filter bank, SNR against a background estimate, adaptive
//! threshold, hangover.

use super::ops::*;

const FRAME_LEN: usize = 256;
const COMPLEN: usize = 12;
const UNIRSHFT: i16 = 7;
const SCALE: f64 = 128.0;
const MAX_16: f64 = 32767.0;
/// `(Word16)(x)`: C's truncating conversion of a constant expression.
const fn w16(x: f64) -> i16 {
    x as i16
}
const TONE_THR: i16 = w16(0.65 * MAX_16);
const SP_EST_COUNT: i16 = 80;
const SP_ACTIVITY_COUNT: i16 = 25;
const ALPHA_SP_UP: i16 = w16((1.0 - 0.85) * MAX_16);
const ALPHA_SP_DOWN: i16 = w16((1.0 - 0.85) * MAX_16);
const NOM_LEVEL: i16 = 2050;
const SPEECH_LEVEL_INIT: i16 = NOM_LEVEL;
const MIN_SPEECH_LEVEL1: i16 = w16(NOM_LEVEL as f64 * 0.063);
const MIN_SPEECH_LEVEL2: i16 = w16(NOM_LEVEL as f64 * 0.2);
const MIN_SPEECH_SNR: i16 = 4096;
const ALPHA_UP1: i16 = w16((1.0 - 0.95) * MAX_16);
const ALPHA_DOWN1: i16 = w16((1.0 - 0.936) * MAX_16);
const ALPHA_UP2: i16 = w16((1.0 - 0.985) * MAX_16);
const ALPHA_DOWN2: i16 = w16((1.0 - 0.943) * MAX_16);
const ALPHA3: i16 = w16((1.0 - 0.95) * MAX_16);
const ALPHA4: i16 = w16((1.0 - 0.9) * MAX_16);
const ALPHA5: i16 = w16((1.0 - 0.5) * MAX_16);
const THR_MIN: i16 = w16(1.6 * SCALE);
const THR_HIGH: i16 = w16(6.0 * SCALE);
const THR_LOW: i16 = w16(1.7 * SCALE);
const NO_P1: i16 = 31744;
const NO_P2: i16 = 19786;
/// `MAX_16*(float)(a)/(float)(b)`: the C evaluates these slopes in single precision.
const fn slope(a: i16, b: i16) -> i16 {
    (32767.0f32 * a as f32 / b as f32) as i16
}
const NO_SLOPE: i16 = slope(THR_LOW - THR_HIGH, NO_P2 - NO_P1);
const SP_CH_MIN: i16 = w16(-0.75 * SCALE);
const SP_CH_MAX: i16 = w16(0.75 * SCALE);
const SP_P1: i16 = 22527;
const SP_P2: i16 = 17832;
const SP_SLOPE: i16 = slope(SP_CH_MAX - SP_CH_MIN, SP_P2 - SP_P1);
const HANG_HIGH: i16 = 12;
const HANG_LOW: i16 = 2;
const HANG_P1: i16 = THR_LOW;
const HANG_P2: i16 = w16(4.0 * SCALE);
const HANG_SLOPE: i16 = slope(HANG_LOW - HANG_HIGH, HANG_P2 - HANG_P1);
const BURST_HIGH: i16 = 8;
const BURST_P1: i16 = THR_HIGH;
const BURST_P2: i16 = THR_LOW;
const BURST_SLOPE: i16 = slope(3 - BURST_HIGH, BURST_P2 - BURST_P1);
const STAT_COUNT: i16 = 20;
const STAT_THR_LEVEL: i16 = 184;
const STAT_THR: i16 = 1000;
const NOISE_MIN: i16 = 40;
const NOISE_MAX: i16 = 20000;
const NOISE_INIT: i16 = 150;
const VAD_POW_LOW: i32 = 30000;
const POW_TONE_THR: i32 = 686080;
const COEFF3: i16 = 13363;
const COEFF5_1: i16 = 21955;
const COEFF5_2: i16 = 6390;

/// `VadVars`.
#[derive(Clone, Debug, PartialEq)]
#[repr(C)]
pub(crate) struct VadState {
    bckr_est: [i16; COMPLEN],
    ave_level: [i16; COMPLEN],
    old_level: [i16; COMPLEN],
    sub_level: [i16; COMPLEN],
    a_data5: [[i16; 2]; 5],
    a_data3: [i16; 6],
    burst_count: i16,
    hang_count: i16,
    stat_count: i16,
    vadreg: i16,
    tone_flag: i16,
    sp_est_cnt: i16,
    sp_max: i16,
    sp_max_cnt: i16,
    speech_level: i16,
    prev_pow_sum: i32,
}

impl VadState {
    /// `wb_vad_reset` (`stat_count` keeps its value there; it starts at 0 here).
    pub(crate) fn new() -> Self {
        Self {
            bckr_est: [NOISE_INIT; COMPLEN],
            ave_level: [NOISE_INIT; COMPLEN],
            old_level: [NOISE_INIT; COMPLEN],
            sub_level: [0; COMPLEN],
            a_data5: [[0; 2]; 5],
            a_data3: [0; 6],
            burst_count: 0,
            hang_count: 0,
            stat_count: 0,
            vadreg: 0,
            tone_flag: 0,
            sp_est_cnt: 0,
            sp_max: 0,
            sp_max_cnt: 0,
            speech_level: SPEECH_LEVEL_INIT,
            prev_pow_sum: 0,
        }
    }

    /// `wb_vad_reset`: every field but `stat_count`.
    pub(crate) fn reset(&mut self) {
        let stat_count = self.stat_count;
        *self = Self::new();
        self.stat_count = stat_count;
    }

    /// `wb_vad_tone_detection`: shifts in the flag `p_gain > 0.65`.
    pub(crate) fn tone_detection(&mut self, p_gain: i16) {
        self.tone_flag >>= 1;
        if p_gain > TONE_THR {
            self.tone_flag |= 0x4000;
        }
    }

    /// `wb_vad`: the VAD decision for one frame (1 = speech).
    pub(crate) fn vad(&mut self, in_buf: &[i16; FRAME_LEN]) -> i16 {
        let mut l = 0i32;
        for &v in in_buf {
            l = l_mac(l, v, v);
        }
        let pow_sum = l_add(l, self.prev_pow_sum);
        self.prev_pow_sum = l;
        if pow_sum < POW_TONE_THR {
            self.tone_flag &= 0x1fff;
        }
        let mut level = [0i16; COMPLEN];
        self.filter_bank(in_buf, &mut level);
        let flag = self.vad_decision(&level, pow_sum);
        let mut l = 0i32;
        for &v in &level[1..] {
            l = l_add(l, v as i32);
        }
        self.estimate_speech(extract_h(l_shl(l, 12)));
        flag
    }

    fn filter_bank(&mut self, input: &[i16; FRAME_LEN], level: &mut [i16; COMPLEN]) {
        let mut b = [0i16; FRAME_LEN];
        for (d, &s) in b.iter_mut().zip(input) {
            *d = s >> 1;
        }
        let d5 = &mut self.a_data5;
        let d3 = &mut self.a_data3;
        for i in 0..FRAME_LEN / 2 {
            filter5(&mut b, 2 * i, 2 * i + 1, &mut d5[0]);
        }
        for i in 0..FRAME_LEN / 4 {
            filter5(&mut b, 4 * i, 4 * i + 2, &mut d5[1]);
            filter5(&mut b, 4 * i + 1, 4 * i + 3, &mut d5[2]);
        }
        for i in 0..FRAME_LEN / 8 {
            filter5(&mut b, 8 * i, 8 * i + 4, &mut d5[3]);
            filter5(&mut b, 8 * i + 2, 8 * i + 6, &mut d5[4]);
            filter3(&mut b, 8 * i + 3, 8 * i + 7, &mut d3[0]);
        }
        for i in 0..FRAME_LEN / 16 {
            filter3(&mut b, 16 * i, 16 * i + 8, &mut d3[1]);
            filter3(&mut b, 16 * i + 4, 16 * i + 12, &mut d3[2]);
            filter3(&mut b, 16 * i + 6, 16 * i + 14, &mut d3[3]);
        }
        for i in 0..FRAME_LEN / 32 {
            filter3(&mut b, 32 * i, 32 * i + 16, &mut d3[4]);
            filter3(&mut b, 32 * i + 8, 32 * i + 24, &mut d3[5]);
        }
        // (band, count1, count2, ind_m, ind_a, scale)
        const BANDS: [(usize, usize, usize, usize, usize, i16); 12] = [
            (11, FRAME_LEN / 4 - 48, FRAME_LEN / 4, 4, 1, 14),
            (10, FRAME_LEN / 8 - 24, FRAME_LEN / 8, 8, 7, 15),
            (9, FRAME_LEN / 8 - 24, FRAME_LEN / 8, 8, 3, 15),
            (8, FRAME_LEN / 8 - 24, FRAME_LEN / 8, 8, 2, 15),
            (7, FRAME_LEN / 16 - 12, FRAME_LEN / 16, 16, 14, 16),
            (6, FRAME_LEN / 16 - 12, FRAME_LEN / 16, 16, 6, 16),
            (5, FRAME_LEN / 16 - 12, FRAME_LEN / 16, 16, 4, 16),
            (4, FRAME_LEN / 16 - 12, FRAME_LEN / 16, 16, 12, 16),
            (3, FRAME_LEN / 32 - 6, FRAME_LEN / 32, 32, 8, 17),
            (2, FRAME_LEN / 32 - 6, FRAME_LEN / 32, 32, 24, 17),
            (1, FRAME_LEN / 32 - 6, FRAME_LEN / 32, 32, 16, 17),
            (0, FRAME_LEN / 32 - 6, FRAME_LEN / 32, 32, 0, 17),
        ];
        for (k, c1, c2, m, a, scale) in BANDS {
            level[k] = level_calculation(&b, &mut self.sub_level[k], c1, c2, m, a, scale);
        }
    }

    fn update_cntrl(&mut self, level: &[i16; COMPLEN]) {
        if (self.tone_flag & 0x7c00) == 0x7c00 || (self.vadreg & 0x7f80) == 0 {
            self.stat_count = STAT_COUNT;
        } else {
            let mut stat_rat = 0i16;
            for i in 0..COMPLEN {
                let (num, denom) = if level[i] > self.ave_level[i] {
                    (level[i], self.ave_level[i])
                } else {
                    (self.ave_level[i], level[i])
                };
                let num = num.max(STAT_THR_LEVEL);
                let denom = denom.max(STAT_THR_LEVEL);
                let exp = norm_s(denom);
                let temp = div_s(num >> 1, shl(denom, exp));
                stat_rat = add(stat_rat, shr(temp, sub(8, exp)));
            }
            if stat_rat > STAT_THR {
                self.stat_count = STAT_COUNT;
            } else if (self.vadreg & 0x4000) != 0 && self.stat_count != 0 {
                self.stat_count = sub(self.stat_count, 1);
            }
        }
        let alpha = if self.stat_count == STAT_COUNT {
            32767
        } else if (self.vadreg & 0x4000) == 0 {
            ALPHA5
        } else {
            ALPHA4
        };
        for i in 0..COMPLEN {
            self.ave_level[i] = add(
                self.ave_level[i],
                mult_r(alpha, sub(level[i], self.ave_level[i])),
            );
        }
    }

    fn hangover_addition(&mut self, low_power: bool, hang_len: i16, burst_len: i16) -> i16 {
        if low_power {
            self.burst_count = 0;
            self.hang_count = 0;
            return 0;
        }
        if (self.vadreg & 0x4000) != 0 {
            self.burst_count = add(self.burst_count, 1);
            if self.burst_count >= burst_len {
                self.hang_count = hang_len;
            }
            return 1;
        }
        self.burst_count = 0;
        if self.hang_count > 0 {
            self.hang_count = sub(self.hang_count, 1);
            return 1;
        }
        0
    }

    fn noise_estimate_update(&mut self, level: &[i16; COMPLEN]) {
        self.update_cntrl(level);
        let (alpha_up, alpha_down, bckr_add) = if (0x7800 & self.vadreg) == 0 {
            (ALPHA_UP1, ALPHA_DOWN1, 2)
        } else if self.stat_count == 0 {
            (ALPHA_UP2, ALPHA_DOWN2, 2)
        } else {
            (0, ALPHA3, 0)
        };
        for i in 0..COMPLEN {
            let temp = sub(self.old_level[i], self.bckr_est[i]);
            let b = &mut self.bckr_est[i];
            if temp < 0 {
                *b = add(-2, add(*b, mult_r(alpha_down, temp))).max(NOISE_MIN);
            } else {
                *b = add(bckr_add, add(*b, mult_r(alpha_up, temp))).min(NOISE_MAX);
            }
        }
        self.old_level = *level;
    }

    fn vad_decision(&mut self, level: &[i16; COMPLEN], pow_sum: i32) -> i16 {
        let mut snr_sum = 0i32;
        for i in 0..COMPLEN {
            let exp = norm_s(self.bckr_est[i]);
            let temp = shl(self.bckr_est[i], exp);
            let temp = div_s(level[i] >> 1, temp);
            let temp = shl(temp, sub(exp, UNIRSHFT - 1));
            snr_sum = l_mac(snr_sum, temp, temp);
        }
        let mut l = 0i32;
        for &b in &self.bckr_est[1..] {
            l = l_add(l, b as i32);
        }
        let noise_level = extract_h(l_shl(l, 12));
        let temp = shl(mult(noise_level, MIN_SPEECH_SNR), 3);
        if self.speech_level < temp {
            self.speech_level = temp;
        }
        let ilog2_noise_level = ilog2(noise_level);
        let ilog2_speech_level = ilog2(sub(self.speech_level, temp));
        let temp = add(mult(NO_SLOPE, sub(ilog2_noise_level, NO_P1)), THR_HIGH);
        let temp2 = add(SP_CH_MIN, mult(SP_SLOPE, sub(ilog2_speech_level, SP_P1)))
            .clamp(SP_CH_MIN, SP_CH_MAX);
        let vad_thr = add(temp, temp2).max(THR_MIN);
        self.vadreg >>= 1;
        if snr_sum > l_mult(vad_thr, 512 * COMPLEN as i16) {
            self.vadreg |= 0x4000;
        }
        let low_power = pow_sum < VAD_POW_LOW;
        self.noise_estimate_update(level);
        let hang_len = add(mult(HANG_SLOPE, sub(vad_thr, HANG_P1)), HANG_HIGH).max(HANG_LOW);
        let burst_len = add(mult(BURST_SLOPE, sub(vad_thr, BURST_P1)), BURST_HIGH);
        self.hangover_addition(low_power, hang_len, burst_len)
    }

    fn estimate_speech(&mut self, in_level: i16) {
        if sub(self.sp_est_cnt, self.sp_max_cnt) > SP_EST_COUNT - SP_ACTIVITY_COUNT {
            self.sp_est_cnt = 0;
            self.sp_max = 0;
            self.sp_max_cnt = 0;
        }
        self.sp_est_cnt = add(self.sp_est_cnt, 1);
        if ((self.vadreg & 0x4000) != 0 || in_level > self.speech_level)
            && in_level > MIN_SPEECH_LEVEL1
        {
            if in_level > self.sp_max {
                self.sp_max = in_level;
            }
            self.sp_max_cnt = add(self.sp_max_cnt, 1);
            if self.sp_max_cnt >= SP_ACTIVITY_COUNT {
                let tmp = self.sp_max >> 1;
                let alpha = if tmp > self.speech_level {
                    ALPHA_SP_UP
                } else {
                    ALPHA_SP_DOWN
                };
                if tmp > MIN_SPEECH_LEVEL2 {
                    self.speech_level = add(
                        self.speech_level,
                        mult_r(alpha, sub(tmp, self.speech_level)),
                    );
                }
                self.sp_max = 0;
                self.sp_max_cnt = 0;
                self.sp_est_cnt = 0;
            }
        }
    }
}

/// `ilog2`: -1024·log2(mant·2^-15)-ish (TS 26.173's approximation).
pub(crate) fn ilog2(mant: i16) -> i16 {
    let mant = mant.max(1);
    let ex = norm_s(mant);
    let mut mant = shl(mant, ex);
    for _ in 0..3 {
        mant = mult(mant, mant);
    }
    let l = l_mult(mant, mant);
    let ex2 = norm_l(l);
    let mant = extract_h(l_shl(l, ex2));
    let res = add(shl(add(ex, 16), 10), shl(ex2, 6));
    sub(add(res, 127), mant >> 8)
}

/// `filter5`: 5th-order half-band low/high-pass pair on `b[i0]`, `b[i1]`.
fn filter5(b: &mut [i16], i0: usize, i1: usize, data: &mut [i16; 2]) {
    let temp0 = sub(b[i0], mult(COEFF5_1, data[0]));
    let temp1 = add(data[0], mult(COEFF5_1, temp0));
    data[0] = temp0;
    let temp0 = sub(b[i1], mult(COEFF5_2, data[1]));
    let temp2 = add(data[1], mult(COEFF5_2, temp0));
    data[1] = temp0;
    b[i0] = extract_h(l_shl(l_add(temp1 as i32, temp2 as i32), 15));
    b[i1] = extract_h(l_shl(l_sub(temp1 as i32, temp2 as i32), 15));
}

/// `filter3`: 3rd-order half-band pair.
fn filter3(b: &mut [i16], i0: usize, i1: usize, data: &mut i16) {
    let temp1 = sub(b[i1], mult(COEFF3, *data));
    let temp2 = add(*data, mult(COEFF3, temp1));
    *data = temp1;
    b[i1] = extract_h(l_shl(l_sub(b[i0] as i32, temp2 as i32), 15));
    b[i0] = extract_h(l_shl(l_add(b[i0] as i32, temp2 as i32), 15));
}

/// `level_calculation`: the band's level (sum of |x|), with the previous frame's tail.
fn level_calculation(
    data: &[i16],
    sub_level: &mut i16,
    count1: usize,
    count2: usize,
    ind_m: usize,
    ind_a: usize,
    scale: i16,
) -> i16 {
    let mut l1 = 0i32;
    for i in count1..count2 {
        l1 = l_mac(l1, 1, abs_s(data[ind_m * i + ind_a]));
    }
    let mut l2 = l_add(l1, l_shl(*sub_level as i32, sub(16, scale)));
    *sub_level = extract_h(l_shl(l1, scale));
    for i in 0..count1 {
        l2 = l_mac(l2, 1, abs_s(data[ind_m * i + ind_a]));
    }
    extract_h(l_shl(l2, scale))
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::encoder::filters::vs_ref::signal;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_wb_vad_init(st: *mut *mut VadState) -> i16;
        fn refwb_wb_vad_exit(st: *mut *mut VadState);
        fn refwb_wb_vad(st: *mut VadState, input: *mut i16) -> i16;
        fn refwb_wb_vad_tone_detection(st: *mut VadState, p_gain: i16);
        fn refwb_ilog2(mant: i16) -> i16;
    }

    #[test]
    fn ilog2_matches_reference() {
        with_ref(|| {
            for m in i16::MIN..=i16::MAX {
                // SAFETY: a value function.
                assert_eq!(ilog2(m), unsafe { refwb_ilog2(m) }, "ilog2 {m}");
            }
        });
    }

    /// Sequences of frames switching between silence, noise and loud classes (the VAD's
    /// hangover, stationarity and speech-level paths), random pitch gains for the tone flags;
    /// the whole state after every frame (`VadVars` mirrored as `#[repr(C)]`).
    #[test]
    fn wb_vad_matches_reference() {
        let mut r = Rng::new(51);
        for seq in 0..30u64 {
            let mut st = VadState::new();
            let mut c: *mut VadState = std::ptr::null_mut();
            // SAFETY: init allocates the state, exit frees it.
            with_ref(|| unsafe { refwb_wb_vad_init(&mut c) });
            for call in 0..400u64 {
                let phase = (call / (5 + seq % 40)) % 5;
                let mut x: [i16; FRAME_LEN] = match phase {
                    0 => [0; FRAME_LEN],
                    1 => signal(&mut r, 3, FRAME_LEN).try_into().unwrap(),
                    2 => signal(&mut r, 2, FRAME_LEN).try_into().unwrap(),
                    _ => signal(&mut r, seq + call, FRAME_LEN).try_into().unwrap(),
                };
                if phase == 1 && seq % 3 == 0 {
                    for (n, v) in x.iter_mut().enumerate() {
                        *v = (8000.0 * (n as f64 * 0.3).sin()) as i16;
                    }
                }
                let g = if r.u64().is_multiple_of(3) {
                    30000
                } else {
                    r.i16()
                };
                let f = st.vad(&x);
                st.tone_detection(g);
                with_ref(|| {
                    // SAFETY: a live state, a 256-sample frame.
                    unsafe {
                        assert_eq!(f, refwb_wb_vad(c, x.as_mut_ptr()), "wb_vad {seq} {call}");
                        refwb_wb_vad_tone_detection(c, g);
                        assert_eq!(st, *c, "VadVars {seq} {call}");
                    }
                });
            }
            // SAFETY: frees the state from init.
            with_ref(|| unsafe { refwb_wb_vad_exit(&mut c) });
        }
    }
}
