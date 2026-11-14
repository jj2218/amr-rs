//! Voice activity detector, option 1 (`vad1`, TS 26.094): filter-bank sub-band levels, noise
//! estimation, tone/pitch/complex-signal detection and the hangover logic.
//!
//! Bit-exact port of opencore-amr 0.1.6 (Apache-2.0) `enc/src/vad1.cpp`.
#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1

use crate::basic_op::{
    Flag, abs_s, add, div_s, l_add, l_mac, l_msu, l_shl, l_sub, mult, mult_r, norm_s, round, shl,
    shr, sub,
};

const FRAME_LEN: usize = 160;
pub(crate) const COMPLEN: usize = 9;
const INV_COMPLEN: i16 = 3641;
/// Samples of look-ahead before the frame that the power measure includes.
pub(crate) const LOOKAHEAD: usize = 40;
const UNIRSHFT: i16 = 6;
const TONE_THR: i16 = 21298;
const ALPHA_UP1: i16 = 1638;
const ALPHA_DOWN1: i16 = 2097;
const ALPHA_UP2: i16 = 491;
const ALPHA_DOWN2: i16 = 1867;
const ALPHA3: i16 = 1638;
const ALPHA4: i16 = 3276;
const ALPHA5: i16 = 16383;
const VAD_THR_HIGH: i16 = 1260;
const VAD_THR_LOW: i16 = 720;
const VAD_P1: i16 = 0;
const VAD_SLOPE: i16 = -2808;
const STAT_COUNT: i16 = 20;
const CAD_MIN_STAT_COUNT: i16 = 5;
const STAT_THR_LEVEL: i16 = 184;
const STAT_THR: i16 = 1000;
const NOISE_MIN: i16 = 40;
const NOISE_MAX: i16 = 16000;
const NOISE_INIT: i16 = 150;
const HANG_NOISE_THR: i16 = 100;
const BURST_LEN_HIGH_NOISE: i16 = 4;
const HANG_LEN_HIGH_NOISE: i16 = 7;
const BURST_LEN_LOW_NOISE: i16 = 5;
const HANG_LEN_LOW_NOISE: i16 = 4;
const VAD_POW_LOW: i32 = 15000;
const POW_PITCH_THR: i32 = 343_040;
const POW_COMPLEX_THR: i32 = 15000;
const COEFF3: i16 = 13363;
const COEFF5_1: i16 = 21955;
const COEFF5_2: i16 = 6390;
const LTHRESH: i16 = 4;
const NTHRESH: i16 = 4;
const CVAD_THRESH_ADAPT_HIGH: i16 = 19660;
const CVAD_THRESH_ADAPT_LOW: i16 = 16383;
const CVAD_THRESH_IN_NOISE: i16 = 21298;
const CVAD_THRESH_HANG: i16 = 22936;
const CVAD_HANG_LIMIT: i16 = 100;
const CVAD_HANG_LENGTH: i16 = 250;
const CVAD_LOWPOW_RESET: i16 = 13106;
const CVAD_MIN_CORR: i16 = 13106;
const CVAD_ADAPT_SLOW: i16 = 655;
const CVAD_ADAPT_FAST: i16 = 2621;
const CVAD_ADAPT_REALLY_FAST: i16 = 6553;

/// The VAD option of a DTX encoder (the reference code chooses at compile time).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Vad {
    One(VadState),
    Two(super::vad2::Vad2State),
}

/// `vadState1` (field order as in C).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct VadState {
    pub bckr_est: [i16; COMPLEN],
    pub ave_level: [i16; COMPLEN],
    pub old_level: [i16; COMPLEN],
    pub sub_level: [i16; COMPLEN],
    pub a_data5: [[i16; 2]; 3],
    pub a_data3: [i16; 5],
    pub burst_count: i16,
    pub hang_count: i16,
    pub stat_count: i16,
    pub vadreg: i16,
    pub pitch: i16,
    pub tone: i16,
    pub complex_high: i16,
    pub complex_low: i16,
    pub oldlag_count: i16,
    pub oldlag: i16,
    pub complex_hang_count: i16,
    pub complex_hang_timer: i16,
    pub best_corr_hp: i16,
    pub speech_vad_decision: i16,
    pub complex_warning: i16,
    pub sp_burst_count: i16,
    pub corr_hp_fast: i16,
}

impl Default for VadState {
    /// `vad1_reset`.
    fn default() -> Self {
        Self {
            bckr_est: [NOISE_INIT; COMPLEN],
            ave_level: [NOISE_INIT; COMPLEN],
            old_level: [NOISE_INIT; COMPLEN],
            sub_level: [0; COMPLEN],
            a_data5: [[0; 2]; 3],
            a_data3: [0; 5],
            burst_count: 0,
            hang_count: 0,
            stat_count: 0,
            vadreg: 0,
            pitch: 0,
            tone: 0,
            complex_high: 0,
            complex_low: 0,
            oldlag_count: 0,
            oldlag: 0,
            complex_hang_count: 0,
            complex_hang_timer: 0,
            best_corr_hp: CVAD_LOWPOW_RESET,
            speech_vad_decision: 0,
            complex_warning: 0,
            sp_burst_count: 0,
            corr_hp_fast: CVAD_LOWPOW_RESET,
        }
    }
}

#[inline(always)]
pub(super) fn mul15(a: i16, b: i16) -> i16 {
    ((a as i32 * b as i32) >> 15) as i16
}

pub(super) fn first_filter_stage(
    inp: &[i16],
    out: &mut [i16; FRAME_LEN],
    data: &mut [i16; 2],
    ov: &mut Flag,
) {
    let mut data0 = data[0];
    let mut data1 = data[1];
    for i in 0..FRAME_LEN / 4 {
        let mut temp0 = mul15(COEFF5_1, data0);
        let mut temp1 = inp[4 * i] >> 2;
        temp0 = sub(temp1, temp0, ov);
        temp1 = mul15(COEFF5_1, temp0);
        temp1 = add(data0, temp1, ov);
        let mut temp3 = mul15(COEFF5_2, data1);
        let mut temp2 = inp[4 * i + 1] >> 2;
        temp3 = sub(temp2, temp3, ov);
        temp2 = mul15(COEFF5_2, temp3);
        temp2 = add(data1, temp2, ov);
        out[4 * i] = add(temp1, temp2, ov);
        out[4 * i + 1] = sub(temp1, temp2, ov);
        temp1 = mul15(COEFF5_1, temp0);
        temp2 = inp[4 * i + 2] >> 2;
        data0 = sub(temp2, temp1, ov);
        temp1 = mul15(COEFF5_1, data0);
        temp1 = add(temp0, temp1, ov);
        data1 = mul15(COEFF5_2, temp3);
        temp2 = inp[4 * i + 3] >> 2;
        data1 = sub(temp2, data1, ov);
        temp2 = mul15(COEFF5_2, data1);
        temp2 = add(temp3, temp2, ov);
        out[4 * i + 2] = add(temp1, temp2, ov);
        out[4 * i + 3] = sub(temp1, temp2, ov);
    }
    data[0] = data0;
    data[1] = data1;
}

pub(super) fn filter5(buf: &mut [i16], i0: usize, i1: usize, data: &mut [i16; 2], ov: &mut Flag) {
    let mut temp0 = mult(COEFF5_1, data[0], ov);
    temp0 = sub(buf[i0], temp0, ov);
    let mut temp1 = mult(COEFF5_1, temp0, ov);
    temp1 = add(data[0], temp1, ov);
    data[0] = temp0;
    temp0 = mult(COEFF5_2, data[1], ov);
    temp0 = sub(buf[i1], temp0, ov);
    let mut temp2 = mult(COEFF5_2, temp0, ov);
    temp2 = add(data[1], temp2, ov);
    data[1] = temp0;
    temp0 = add(temp1, temp2, ov);
    buf[i0] = shr(temp0, 1, ov);
    temp0 = sub(temp1, temp2, ov);
    buf[i1] = shr(temp0, 1, ov);
}

pub(super) fn filter3(buf: &mut [i16], i0: usize, i1: usize, data: &mut i16, ov: &mut Flag) {
    let mut temp1 = mult(COEFF3, *data, ov);
    temp1 = sub(buf[i1], temp1, ov);
    let mut temp2 = mult(COEFF3, temp1, ov);
    temp2 = add(*data, temp2, ov);
    *data = temp1;
    temp1 = sub(buf[i0], temp2, ov);
    buf[i1] = shr(temp1, 1, ov);
    temp1 = add(buf[i0], temp2, ov);
    buf[i0] = shr(temp1, 1, ov);
}

#[allow(clippy::too_many_arguments)] // mirrors the C signature
pub(super) fn level_calculation(
    data: &[i16],
    sub_level: &mut i16,
    count1: usize,
    count2: usize,
    ind_m: usize,
    ind_a: usize,
    scale: i16,
    ov: &mut Flag,
) -> i16 {
    let mut l_temp1: i32 = 0;
    for i in count1..count2 {
        l_temp1 = l_mac(l_temp1, 1, abs_s(data[ind_m * i + ind_a]), ov);
    }
    let shifted = l_shl(*sub_level as i32, sub(16, scale, ov), ov);
    let mut l_temp2 = l_add(l_temp1, shifted, ov);
    *sub_level = (l_shl(l_temp1, scale, ov) >> 16) as i16;
    for i in 0..count1 {
        l_temp2 = l_mac(l_temp2, 1, abs_s(data[ind_m * i + ind_a]), ov);
    }
    (l_shl(l_temp2, scale, ov) >> 16) as i16
}

pub(super) fn filter_bank(
    st: &mut VadState,
    inp: &[i16],
    level: &mut [i16; COMPLEN],
    ov: &mut Flag,
) {
    let mut tmp = [0i16; FRAME_LEN];
    first_filter_stage(inp, &mut tmp, &mut st.a_data5[0], ov);
    for i in 0..FRAME_LEN / 4 {
        filter5(&mut tmp, 4 * i, 4 * i + 2, &mut st.a_data5[1], ov);
        filter5(&mut tmp, 4 * i + 1, 4 * i + 3, &mut st.a_data5[2], ov);
    }
    for i in 0..FRAME_LEN / 8 {
        filter3(&mut tmp, 8 * i, 8 * i + 4, &mut st.a_data3[0], ov);
        filter3(&mut tmp, 8 * i + 2, 8 * i + 6, &mut st.a_data3[1], ov);
        filter3(&mut tmp, 8 * i + 3, 8 * i + 7, &mut st.a_data3[4], ov);
    }
    for i in 0..FRAME_LEN / 16 {
        filter3(&mut tmp, 16 * i, 16 * i + 8, &mut st.a_data3[2], ov);
        filter3(&mut tmp, 16 * i + 4, 16 * i + 12, &mut st.a_data3[3], ov);
    }
    let sl = &mut st.sub_level;
    level[8] = level_calculation(
        &tmp,
        &mut sl[8],
        FRAME_LEN / 4 - 8,
        FRAME_LEN / 4,
        4,
        1,
        15,
        ov,
    );
    level[7] = level_calculation(
        &tmp,
        &mut sl[7],
        FRAME_LEN / 8 - 4,
        FRAME_LEN / 8,
        8,
        7,
        16,
        ov,
    );
    level[6] = level_calculation(
        &tmp,
        &mut sl[6],
        FRAME_LEN / 8 - 4,
        FRAME_LEN / 8,
        8,
        3,
        16,
        ov,
    );
    level[5] = level_calculation(
        &tmp,
        &mut sl[5],
        FRAME_LEN / 8 - 4,
        FRAME_LEN / 8,
        8,
        2,
        16,
        ov,
    );
    level[4] = level_calculation(
        &tmp,
        &mut sl[4],
        FRAME_LEN / 8 - 4,
        FRAME_LEN / 8,
        8,
        6,
        16,
        ov,
    );
    level[3] = level_calculation(
        &tmp,
        &mut sl[3],
        FRAME_LEN / 16 - 2,
        FRAME_LEN / 16,
        16,
        4,
        16,
        ov,
    );
    level[2] = level_calculation(
        &tmp,
        &mut sl[2],
        FRAME_LEN / 16 - 2,
        FRAME_LEN / 16,
        16,
        12,
        16,
        ov,
    );
    level[1] = level_calculation(
        &tmp,
        &mut sl[1],
        FRAME_LEN / 16 - 2,
        FRAME_LEN / 16,
        16,
        8,
        16,
        ov,
    );
    level[0] = level_calculation(
        &tmp,
        &mut sl[0],
        FRAME_LEN / 16 - 2,
        FRAME_LEN / 16,
        16,
        0,
        16,
        ov,
    );
}

pub(super) fn update_cntrl(st: &mut VadState, level: &[i16; COMPLEN], ov: &mut Flag) {
    if st.complex_warning != 0 && st.stat_count < CAD_MIN_STAT_COUNT {
        st.stat_count = CAD_MIN_STAT_COUNT;
    }
    if (st.pitch & 0x6000) == 0x6000 || (st.tone & 0x7c00) == 0x7c00 || (st.vadreg & 0x7f80) == 0 {
        st.stat_count = STAT_COUNT;
    } else {
        let mut stat_rat: i16 = 0;
        for i in 0..COMPLEN {
            let (mut num, mut denom) = if level[i] > st.ave_level[i] {
                (level[i], st.ave_level[i])
            } else {
                (st.ave_level[i], level[i])
            };
            num = num.max(STAT_THR_LEVEL);
            denom = denom.max(STAT_THR_LEVEL);
            let exp = norm_s(denom);
            denom = shl(denom, exp, ov);
            let mut temp = shr(num, 1, ov);
            temp = div_s(temp, denom);
            let s = sub(8, exp, ov);
            let t = shr(temp, s, ov);
            stat_rat = add(stat_rat, t, ov);
        }
        if stat_rat > STAT_THR {
            st.stat_count = STAT_COUNT;
        } else if (st.vadreg & 0x4000) != 0 && st.stat_count != 0 {
            st.stat_count = sub(st.stat_count, 1, ov);
        }
    }
    let alpha = if st.stat_count == STAT_COUNT {
        32767
    } else if (st.vadreg & 0x4000) == 0 {
        ALPHA5
    } else {
        ALPHA4
    };
    for i in 0..COMPLEN {
        let mut temp = sub(level[i], st.ave_level[i], ov);
        temp = mult_r(alpha, temp, ov);
        st.ave_level[i] = add(st.ave_level[i], temp, ov);
    }
}

pub(super) fn hangover_addition(
    st: &mut VadState,
    noise_level: i16,
    low_power: bool,
    ov: &mut Flag,
) -> i16 {
    let (burst_len, hang_len) = if noise_level > HANG_NOISE_THR {
        (BURST_LEN_HIGH_NOISE, HANG_LEN_HIGH_NOISE)
    } else {
        (BURST_LEN_LOW_NOISE, HANG_LEN_LOW_NOISE)
    };
    if low_power {
        st.burst_count = 0;
        st.hang_count = 0;
        st.complex_hang_count = 0;
        st.complex_hang_timer = 0;
        return 0;
    }
    if st.complex_hang_timer > CVAD_HANG_LIMIT && st.complex_hang_count < CVAD_HANG_LENGTH {
        st.complex_hang_count = CVAD_HANG_LENGTH;
    }
    if st.complex_hang_count != 0 {
        st.burst_count = BURST_LEN_HIGH_NOISE;
        st.complex_hang_count = sub(st.complex_hang_count, 1, ov);
        return 1;
    } else if (st.vadreg & 0x3ff0) == 0 && st.corr_hp_fast > CVAD_THRESH_IN_NOISE {
        return 1;
    }
    if (st.vadreg & 0x4000) != 0 {
        st.burst_count = add(st.burst_count, 1, ov);
        if st.burst_count >= burst_len {
            st.hang_count = hang_len;
        }
        return 1;
    }
    st.burst_count = 0;
    if st.hang_count > 0 {
        st.hang_count = sub(st.hang_count, 1, ov);
        return 1;
    }
    0
}

pub(super) fn noise_estimate_update(st: &mut VadState, level: &[i16; COMPLEN], ov: &mut Flag) {
    update_cntrl(st, level, ov);
    let mut bckr_add = 2;
    let (alpha_up, alpha_down) =
        if (0x7800 & st.vadreg) == 0 && (st.pitch & 0x7800) == 0 && st.complex_hang_count == 0 {
            (ALPHA_UP1, ALPHA_DOWN1)
        } else if st.stat_count == 0 && st.complex_hang_count == 0 {
            (ALPHA_UP2, ALPHA_DOWN2)
        } else {
            bckr_add = 0;
            (0, ALPHA3)
        };
    for i in 0..COMPLEN {
        let mut temp = sub(st.old_level[i], st.bckr_est[i], ov);
        if temp < 0 {
            temp = mult_r(alpha_down, temp, ov);
            temp = add(st.bckr_est[i], temp, ov);
            st.bckr_est[i] = add(-2, temp, ov).max(NOISE_MIN);
        } else {
            temp = mult_r(alpha_up, temp, ov);
            temp = add(st.bckr_est[i], temp, ov);
            st.bckr_est[i] = add(bckr_add, temp, ov).min(NOISE_MAX);
        }
    }
    st.old_level = *level;
}

pub(super) fn complex_estimate_adapt(st: &mut VadState, low_power: bool, ov: &mut Flag) {
    let alpha = if st.best_corr_hp < st.corr_hp_fast {
        if st.corr_hp_fast < CVAD_THRESH_ADAPT_HIGH {
            CVAD_ADAPT_FAST
        } else {
            CVAD_ADAPT_REALLY_FAST
        }
    } else if st.corr_hp_fast < CVAD_THRESH_ADAPT_HIGH {
        CVAD_ADAPT_FAST
    } else {
        CVAD_ADAPT_SLOW
    };
    let mut l_tmp = (st.corr_hp_fast as i32) << 16;
    l_tmp = l_msu(l_tmp, alpha, st.corr_hp_fast, ov);
    l_tmp = l_mac(l_tmp, alpha, st.best_corr_hp, ov);
    st.corr_hp_fast = round(l_tmp, ov).max(CVAD_MIN_CORR);
    if low_power {
        st.corr_hp_fast = CVAD_MIN_CORR;
    }
}

pub(super) fn complex_vad(st: &mut VadState, low_power: bool, ov: &mut Flag) -> i16 {
    st.complex_high = shr(st.complex_high, 1, ov);
    st.complex_low = shr(st.complex_low, 1, ov);
    if !low_power {
        if st.corr_hp_fast > CVAD_THRESH_ADAPT_HIGH {
            st.complex_high |= 0x4000;
        }
        if st.corr_hp_fast > CVAD_THRESH_ADAPT_LOW {
            st.complex_low |= 0x4000;
        }
    }
    if st.corr_hp_fast > CVAD_THRESH_HANG {
        st.complex_hang_timer = add(st.complex_hang_timer, 1, ov);
    } else {
        st.complex_hang_timer = 0;
    }
    ((st.complex_high & 0x7f80) == 0x7f80 || (st.complex_low & 0x7fff) == 0x7fff) as i16
}

pub(super) fn vad_decision(
    st: &mut VadState,
    level: &[i16; COMPLEN],
    pow_sum: i32,
    ov: &mut Flag,
) -> i16 {
    let mut l_temp: i32 = 0;
    for i in 0..COMPLEN {
        let exp = norm_s(st.bckr_est[i]);
        let mut temp = shl(st.bckr_est[i], exp, ov);
        temp = div_s(shr(level[i], 1, ov), temp);
        let s = sub(exp, UNIRSHFT - 1, ov);
        temp = shl(temp, s, ov);
        l_temp = l_mac(l_temp, temp, temp, ov);
    }
    let mut snr_sum = (l_shl(l_temp, 6, ov) >> 16) as i16;
    snr_sum = mult(snr_sum, INV_COMPLEN, ov);
    let mut l_temp: i32 = 0;
    for i in 0..COMPLEN {
        l_temp = l_add(l_temp, st.bckr_est[i] as i32, ov);
    }
    let noise_level = (l_shl(l_temp, 13, ov) >> 16) as i16;
    let mut temp1 = sub(noise_level, VAD_P1, ov);
    temp1 = mult(VAD_SLOPE, temp1, ov);
    let vad_thr = add(temp1, VAD_THR_HIGH, ov).max(VAD_THR_LOW);
    st.vadreg = shr(st.vadreg, 1, ov);
    if snr_sum > vad_thr {
        st.vadreg |= 0x4000;
    }
    let low_power = l_sub(pow_sum, VAD_POW_LOW, ov) < 0;
    complex_estimate_adapt(st, low_power, ov);
    st.complex_warning = complex_vad(st, low_power, ov);
    noise_estimate_update(st, level, ov);
    st.speech_vad_decision = hangover_addition(st, noise_level, low_power, ov);
    st.speech_vad_decision
}

/// `vad_complex_detection_update`.
pub(crate) fn vad_complex_detection_update(st: &mut VadState, best_corr_hp: i16) {
    st.best_corr_hp = best_corr_hp;
}

/// `vad_tone_detection`.
pub(crate) fn vad_tone_detection(st: &mut VadState, t0: i32, t1: i32, ov: &mut Flag) {
    let temp = round(t1, ov);
    if temp > 0 && l_msu(t0, temp, TONE_THR, ov) > 0 {
        st.tone |= 0x4000;
    }
}

/// `vad_tone_detection_update`.
pub(crate) fn vad_tone_detection_update(st: &mut VadState, one_lag_per_frame: bool, ov: &mut Flag) {
    st.tone = shr(st.tone, 1, ov);
    if one_lag_per_frame {
        st.tone = shr(st.tone, 1, ov);
        st.tone |= 0x2000;
    }
}

/// `vad_pitch_detection`.
pub(crate) fn vad_pitch_detection(st: &mut VadState, t_op: &[i16; 2], ov: &mut Flag) {
    let mut lagcount: i16 = 0;
    for &t in t_op {
        let temp = abs_s(sub(st.oldlag, t, ov));
        if temp < LTHRESH {
            lagcount += 1;
        }
        st.oldlag = t;
    }
    st.pitch = shr(st.pitch, 1, ov);
    if add(st.oldlag_count, lagcount, ov) >= NTHRESH {
        st.pitch |= 0x4000;
    }
    st.oldlag_count = lagcount;
}

/// `vad1`: VAD decision for one frame. `buf` holds `LOOKAHEAD` samples before the frame, then
/// the frame (`FRAME_LEN` samples).
pub(crate) fn vad1(st: &mut VadState, buf: &[i16], ov: &mut Flag) -> i16 {
    let mut level = [0i16; COMPLEN];
    let mut pow_sum: i32 = 0;
    for &x in &buf[..FRAME_LEN] {
        pow_sum = l_mac(pow_sum, x, x, ov);
    }
    if l_sub(pow_sum, POW_PITCH_THR, ov) < 0 {
        st.pitch &= 0x3fff;
    }
    if l_sub(pow_sum, POW_COMPLEX_THR, ov) < 0 {
        st.complex_low &= 0x3fff;
    }
    filter_bank(st, &buf[LOOKAHEAD..], &mut level, ov);
    vad_decision(st, &level, pow_sum, ov)
}
