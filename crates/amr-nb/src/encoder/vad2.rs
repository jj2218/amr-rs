//! Voice activity detector, option 2 (`vad2`, TS 26.094 §B): 16-channel FFT energy estimator
//! with noise tracking, plus `LTP_flag_update` and the real FFT (`r_fft`) it uses.
//!
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `common/src/{vad2,r_fft}.cpp` and
//! `enc/src/lflg_upd.cpp` (opencore ships these but builds VAD1 only).
#![allow(clippy::needless_range_loop)] // index loops mirror the C code 1:1

use crate::basic_op::{
    Flag, abs_s, add, div_s, l_add, l_extract, l_mac, l_msu, l_mult, l_negate, l_shl, l_shr,
    l_shr_r, l_sub, mpy_32_16, mult, mult_r, norm_s, round, shl, shr, shr_r, sub,
};
use crate::cnst::CodecMode;
use crate::math::{log2, pow2};
use crate::tables::{
    BURSTCOUNT_TABLE, CH_TBL_SH, HANGOVER_TABLE, II_TABLE, PHS_TBL, VM_TBL, VM_THRESHOLD_TABLE,
};

const FRM_LEN: usize = 80;
const DELAY: usize = 24;
const FFT_LEN: usize = 128;
const NUM_CHAN: usize = 16;
const UPDATE_THLD: i16 = 35;
const HYSTER_CNT_THLD: i16 = 6;
const UPDATE_CNT_THLD: i16 = 50;
const MIN_NOISE_ENRG_0: i32 = 32;
const INE_NOISE_0: i32 = 8192;
const FRACTIONAL_BITS_0: i16 = 9;
const HIGH_ALPHA: i16 = 29491;
const LOW_ALPHA: i16 = 22938;
const ALPHA_RANGE: i16 = HIGH_ALPHA - LOW_ALPHA;
const DEV_THLD: i16 = 7168;
const PRE_EMP_FAC: i16 = -26214;
const CEE_SM_FAC: i16 = 18022;
const ONE_MINUS_CEE_SM_FAC: i16 = 14746;
const CNE_SM_FAC: i16 = 3277;
const ONE_MINUS_CNE_SM_FAC: i16 = 29491;
const FFT_HEADROOM: i16 = 2;
/// Per shift state: `{NOISE_FLOOR_CHAN, MIN_CHAN_ENRG, INE_NOISE, FRACTIONAL_BITS}`.
const NOISE_FLOOR_CHAN: [i32; 2] = [512, 16];
const MIN_CHAN_ENRG: [i32; 2] = [32, 1];
const INE_NOISE: [i32; 2] = [8192, 256];
const FBITS: [i16; 2] = [9, 4];
/// `STATE_1_TO_0_SHIFT_R`, `STATE_0_TO_1_SHIFT_R`.
const STATE_CHANGE_SHIFT_R: [i16; 2] = [4 - 9, 9 - 4];
const ENRG_NORM_SHIFT: [i16; 2] = [9 - 1 + 2, 4 - 1 + 2];

const CH_TBL: [[usize; 2]; NUM_CHAN] = [
    [2, 3],
    [4, 5],
    [6, 7],
    [8, 9],
    [10, 11],
    [12, 13],
    [14, 16],
    [17, 19],
    [20, 22],
    [23, 26],
    [27, 30],
    [31, 35],
    [36, 41],
    [42, 48],
    [49, 55],
    [56, 63],
];

/// `vadState2` (field order as in C).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct Vad2State {
    pub pre_emp_mem: i16,
    pub update_cnt: i16,
    pub hyster_cnt: i16,
    pub last_update_cnt: i16,
    pub ch_enrg_long_db: [i16; NUM_CHAN],
    pub l_frame_cnt: i32,
    pub l_ch_enrg: [i32; NUM_CHAN],
    pub l_ch_noise: [i32; NUM_CHAN],
    pub last_normb_shift: i16,
    pub tsnr: i16,
    pub hangover: i16,
    pub burstcount: i16,
    pub fupdate_flag: i16,
    pub neg_snr_var: i16,
    pub neg_snr_bias: i16,
    pub shift_state: i16,
    pub l_r0: i32,
    pub l_rmax: i32,
    pub ltp_flag: bool,
}

/// `c_fft`: in-place complex FFT of 64 points (interleaved re/im), scaled by 1/2 per stage.
pub(super) fn c_fft(x: &mut [i16; FFT_LEN], ov: &mut Flag) {
    const SIZE: i16 = FFT_LEN as i16;
    let mut j: i16 = 0;
    let mut i: i16 = 0;
    while i < SIZE - 2 {
        if j > i {
            x.swap(i as usize, j as usize);
            x.swap(i as usize + 1, j as usize + 1);
        }
        let mut k = SIZE / 2;
        while j >= k {
            j = sub(j, k, ov);
            k = shr(k, 1, ov);
        }
        j = add(j, k, ov);
        i += 2;
    }
    for (stage, &ii) in II_TABLE.iter().enumerate() {
        let jj = shl(2, stage as i16, ov);
        let kk = shl(jj, 1, ov);
        let ii2 = shl(ii, 1, ov);
        let mut ji: i16 = 0;
        let mut j = 0;
        while j < jj {
            let mut k = j;
            while k < SIZE {
                let kj = add(k, jj, ov) as usize;
                let (k_, ji_) = (k as usize, ji as usize);
                let mut ftmp_real = l_mult(x[kj], PHS_TBL[ji_], ov);
                ftmp_real = l_msu(ftmp_real, x[kj + 1], PHS_TBL[ji_ + 1], ov);
                let mut ftmp_imag = l_mult(x[kj + 1], PHS_TBL[ji_], ov);
                ftmp_imag = l_mac(ftmp_imag, x[kj], PHS_TBL[ji_ + 1], ov);
                let tmp1 = round(ftmp_real, ov);
                let tmp2 = round(ftmp_imag, ov);
                let t = sub(x[k_], tmp1, ov);
                x[kj] = shr(t, 1, ov);
                let t = sub(x[k_ + 1], tmp2, ov);
                x[kj + 1] = shr(t, 1, ov);
                let t = add(x[k_], tmp1, ov);
                x[k_] = shr(t, 1, ov);
                let t = add(x[k_ + 1], tmp2, ov);
                x[k_ + 1] = shr(t, 1, ov);
                k += kk;
            }
            ji = add(ji, ii2, ov);
            j += 2;
        }
    }
}

/// `r_fft`: in-place real FFT of 128 points (packed: `x[0]` DC, `x[1]` Nyquist).
pub(super) fn r_fft(x: &mut [i16; FFT_LEN], ov: &mut Flag) {
    c_fft(x, ov);
    let (f1r, f2r) = (x[0], x[1]);
    x[0] = add(f1r, f2r, ov);
    x[1] = sub(f1r, f2r, ov);
    let mut i = 2;
    while i <= FFT_LEN / 2 {
        let j = FFT_LEN - i;
        let ftmp1_real = add(x[i], x[j], ov);
        let ftmp1_imag = sub(x[i + 1], x[j + 1], ov);
        let ftmp2_real = add(x[i + 1], x[j + 1], ov);
        let ftmp2_imag = sub(x[j], x[i], ov);
        let l1r = (ftmp1_real as i32) << 16;
        let l1i = (ftmp1_imag as i32) << 16;
        let mut t = l_mac(l1r, ftmp2_real, PHS_TBL[i], ov);
        t = l_msu(t, ftmp2_imag, PHS_TBL[i + 1], ov);
        x[i] = round(l_shr(t, 1, ov), ov);
        let mut t = l_mac(l1i, ftmp2_imag, PHS_TBL[i], ov);
        t = l_mac(t, ftmp2_real, PHS_TBL[i + 1], ov);
        x[i + 1] = round(l_shr(t, 1, ov), ov);
        let mut t = l_mac(l1r, ftmp2_real, PHS_TBL[j], ov);
        t = l_mac(t, ftmp2_imag, PHS_TBL[j + 1], ov);
        x[j] = round(l_shr(t, 1, ov), ov);
        let mut t = l_negate(l1i);
        t = l_msu(t, ftmp2_imag, PHS_TBL[j], ov);
        t = l_mac(t, ftmp2_real, PHS_TBL[j + 1], ov);
        x[j + 1] = round(l_shr(t, 1, ov), ov);
        i += 2;
    }
}

/// `fn10Log10`: 10 log10 of a value with `fbits` fractional bits (Q8 result).
fn fn10_log10(l_input: i32, fbits: i16, ov: &mut Flag) -> i16 {
    let (integer, fraction) = log2(l_input, ov);
    let integer = sub(integer, fbits, ov);
    let l = mpy_32_16(integer, fraction, 24660, ov);
    l_shr_r(l, 6, ov) as i16
}

/// `block_norm`: normalises `inp` to `headroom` bits below full scale; returns the shift.
fn block_norm(inp: &[i16], out: &mut [i16], headroom: i16, ov: &mut Flag) -> i16 {
    let max = inp.iter().map(|&v| abs_s(v)).max().unwrap();
    if max != 0 {
        let scnt = sub(norm_s(max), headroom, ov);
        for (o, &v) in out.iter_mut().zip(inp) {
            *o = shl(v, scnt, ov);
        }
        scnt
    } else {
        out.fill(0);
        sub(16, headroom, ov)
    }
}

/// `vad2`: VAD option 2 decision for 80 samples (half a frame).
pub(crate) fn vad2(farray: &[i16], st: &mut Vad2State, ov: &mut Flag) -> bool {
    let mut input_buffer = [0i16; FRM_LEN];
    let mut data_buffer = [0i16; FFT_LEN];
    let mut ch_snr = [0i16; NUM_CHAN];
    let mut ch_enrg_db = [0i16; NUM_CHAN];

    st.l_frame_cnt = l_add(st.l_frame_cnt, 1, ov);
    let normb_shift = block_norm(&farray[..FRM_LEN], &mut input_buffer, FFT_HEADROOM, ov);
    let s = sub(st.last_normb_shift, normb_shift, ov);
    st.pre_emp_mem = shr_r(st.pre_emp_mem, s, ov);
    st.last_normb_shift = normb_shift;
    let t = mult(PRE_EMP_FAC, st.pre_emp_mem, ov);
    data_buffer[DELAY] = add(input_buffer[0], t, ov);
    for j in 1..FRM_LEN {
        let t = mult(PRE_EMP_FAC, input_buffer[j - 1], ov);
        data_buffer[DELAY + j] = add(input_buffer[j], t, ov);
    }
    st.pre_emp_mem = input_buffer[FRM_LEN - 1];
    r_fft(&mut data_buffer, ov);

    let mut state_change = false;
    if st.shift_state == 0 {
        if normb_shift <= -FFT_HEADROOM + 2 {
            state_change = true;
            st.shift_state = 1;
        }
    } else if normb_shift >= -FFT_HEADROOM + 5 {
        state_change = true;
        st.shift_state = 0;
    }
    let ss = st.shift_state as usize;
    if state_change {
        for e in st.l_ch_enrg.iter_mut() {
            *e = l_shr(*e, STATE_CHANGE_SHIFT_R[ss], ov);
        }
    }
    let (alpha, one_m_alpha) = if st.l_frame_cnt == 1 {
        (32767, 0)
    } else {
        (CEE_SM_FAC, ONE_MINUS_CEE_SM_FAC)
    };
    for i in 0..NUM_CHAN {
        let mut lenrg: i32 = 0;
        for j in CH_TBL[i][0]..=CH_TBL[i][1] {
            lenrg = l_mac(lenrg, data_buffer[2 * j], data_buffer[2 * j], ov);
            lenrg = l_mac(lenrg, data_buffer[2 * j + 1], data_buffer[2 * j + 1], ov);
        }
        let tmp = shl(normb_shift, 1, ov);
        let tmp = sub(tmp, ENRG_NORM_SHIFT[ss], ov);
        lenrg = l_shr_r(lenrg, tmp, ov);
        let tmp = mult(alpha, CH_TBL_SH[i], ov);
        let (hi, lo) = l_extract(lenrg, ov);
        let ltmp = mpy_32_16(hi, lo, tmp, ov);
        let (hi, lo) = l_extract(st.l_ch_enrg[i], ov);
        let ltmp1 = mpy_32_16(hi, lo, one_m_alpha, ov);
        st.l_ch_enrg[i] = l_add(ltmp, ltmp1, ov);
        if st.l_ch_enrg[i] < MIN_CHAN_ENRG[ss] {
            st.l_ch_enrg[i] = MIN_CHAN_ENRG[ss];
        }
    }
    let mut ltce: i32 = 0;
    for i in 0..NUM_CHAN {
        ltce = l_add(ltce, st.l_ch_enrg[i], ov);
    }
    let mut lpeak: i32 = 0;
    for i in 2..NUM_CHAN {
        if l_sub(st.l_ch_enrg[i], lpeak, ov) > 0 {
            lpeak = st.l_ch_enrg[i];
        }
    }
    let (hi, lo) = l_extract(ltce, ov);
    let ltmp = mpy_32_16(hi, lo, 20480, ov);
    let p2a_flag = l_sub(lpeak, ltmp, ov) > 0;

    if st.l_frame_cnt <= 4 {
        for i in 0..NUM_CHAN {
            st.l_ch_noise[i] = if p2a_flag || st.l_ch_enrg[i] < INE_NOISE[ss] {
                INE_NOISE_0
            } else if st.shift_state == 1 {
                l_shr(st.l_ch_enrg[i], STATE_CHANGE_SHIFT_R[0], ov)
            } else {
                st.l_ch_enrg[i]
            };
        }
    }

    let mut vm_sum: i16 = 0;
    for i in 0..NUM_CHAN {
        ch_enrg_db[i] = fn10_log10(st.l_ch_enrg[i], FBITS[ss], ov);
        let ch_noise_db = fn10_log10(st.l_ch_noise[i], FRACTIONAL_BITS_0, ov);
        ch_snr[i] = sub(ch_enrg_db[i], ch_noise_db, ov);
        let tmp = mult(21845, ch_snr[i], ov);
        let ch_snrq = shr_r(tmp, 6, ov);
        let j = ch_snrq.clamp(0, 89) as usize;
        vm_sum = add(vm_sum, VM_TBL[j], ov);
    }

    let xt;
    if st.l_frame_cnt <= 4 || st.fupdate_flag != 0 {
        let tce_db: i16 = 14320;
        st.neg_snr_var = 0;
        st.neg_snr_bias = 0;
        let mut ltne: i32 = 0;
        for i in 0..NUM_CHAN {
            ltne = l_add(ltne, st.l_ch_noise[i], ov);
        }
        let tne_db = fn10_log10(ltne, FRACTIONAL_BITS_0, ov);
        xt = sub(tce_db, tne_db, ov);
        st.tsnr = xt;
    } else {
        let mut ltmp1: i32 = 0;
        for i in 0..NUM_CHAN {
            let mut ltmp2 = l_mult(ch_snr[i], 10885, ov);
            ltmp2 = l_shr(ltmp2, 8, ov);
            let (hi, lo) = l_extract(ltmp2, ov);
            let hi = add(hi, 3, ov);
            ltmp2 = pow2(hi, lo, ov);
            ltmp1 = l_add(ltmp1, ltmp2, ov);
        }
        xt = fn10_log10(ltmp1, 4 + 3, ov);
        if xt > st.tsnr {
            let a = l_mult(29491, st.tsnr, ov);
            let b = l_mult(3277, xt, ov);
            st.tsnr = round(l_add(a, b, ov), ov);
        } else {
            let tmp = mult(20480, st.tsnr, ov);
            if sub(xt, tmp, ov) > 0 {
                let a = l_mult(32702, st.tsnr, ov);
                let b = l_mult(66, xt, ov);
                st.tsnr = round(l_add(a, b, ov), ov);
            }
        }
    }
    let tmp = mult(st.tsnr, 10923, ov);
    let tsnrq = shr(tmp, 8, ov).clamp(0, 19) as usize;

    if xt < 0 {
        let mut ltmp1 = l_mult(xt, xt, ov);
        ltmp1 = l_shl(ltmp1, 7, ov);
        let tmp = round(ltmp1, ov);
        let a = l_mult(32440, st.neg_snr_var, ov);
        let b = l_mult(328, tmp, ov);
        st.neg_snr_var = round(l_add(a, b, ov), ov).min(1024);
        let mut tmp = sub(st.neg_snr_var, 166, ov);
        tmp = shl(tmp, 4, ov);
        tmp = mult_r(tmp, 24576, ov);
        st.neg_snr_bias = if tmp < 0 { 0 } else { shr(tmp, 8, ov) };
    }

    let tmp = add(VM_THRESHOLD_TABLE[tsnrq], st.neg_snr_bias, ov);
    let ivad = if vm_sum > tmp {
        st.burstcount = add(st.burstcount, 1, ov);
        if st.burstcount > BURSTCOUNT_TABLE[tsnrq] {
            st.hangover = HANGOVER_TABLE[tsnrq];
        }
        true
    } else {
        st.burstcount = 0;
        st.hangover = sub(st.hangover, 1, ov);
        if st.hangover <= 0 {
            st.hangover = 0;
            false
        } else {
            true
        }
    };

    let mut ch_enrg_dev: i16 = 0;
    if st.l_frame_cnt == 1 {
        st.ch_enrg_long_db = ch_enrg_db;
    } else {
        for i in 0..NUM_CHAN {
            let tmp = abs_s(sub(st.ch_enrg_long_db[i], ch_enrg_db[i], ov));
            ch_enrg_dev = add(ch_enrg_dev, tmp, ov);
        }
    }

    let tmp = sub(st.tsnr, xt, ov);
    let (alpha, one_m_alpha) = if tmp <= 0 || st.tsnr <= 0 {
        (HIGH_ALPHA, (32768 - HIGH_ALPHA as i32) as i16)
    } else if tmp > st.tsnr {
        (LOW_ALPHA, (32768 - LOW_ALPHA as i32) as i16)
    } else {
        let t = div_s(tmp, st.tsnr);
        let t = mult(ALPHA_RANGE, t, ov);
        let alpha = sub(HIGH_ALPHA, t, ov);
        (alpha, sub(32767, alpha, ov))
    };
    for i in 0..NUM_CHAN {
        let a = l_mult(one_m_alpha, ch_enrg_db[i], ov);
        let b = l_mult(alpha, st.ch_enrg_long_db[i], ov);
        st.ch_enrg_long_db[i] = round(l_add(a, b, ov), ov);
    }

    let mut update_flag = false;
    st.fupdate_flag = 0;
    if vm_sum <= UPDATE_THLD {
        if st.burstcount == 0 {
            update_flag = true;
            st.update_cnt = 0;
        }
    } else if l_sub(ltce, NOISE_FLOOR_CHAN[ss], ov) > 0
        && ch_enrg_dev < DEV_THLD
        && !p2a_flag
        && !st.ltp_flag
    {
        st.update_cnt = add(st.update_cnt, 1, ov);
        if st.update_cnt >= UPDATE_CNT_THLD {
            update_flag = true;
            st.fupdate_flag = 1;
        }
    }
    if st.update_cnt == st.last_update_cnt {
        st.hyster_cnt = add(st.hyster_cnt, 1, ov);
    } else {
        st.hyster_cnt = 0;
    }
    st.last_update_cnt = st.update_cnt;
    if st.hyster_cnt > HYSTER_CNT_THLD {
        st.update_cnt = 0;
    }
    if update_flag {
        let tmp = if st.shift_state == 1 {
            STATE_CHANGE_SHIFT_R[0]
        } else {
            0
        };
        for i in 0..NUM_CHAN {
            let l1 = l_shr(st.l_ch_enrg[i], tmp, ov);
            let (hi, lo) = l_extract(l1, ov);
            let ltmp = mpy_32_16(hi, lo, CNE_SM_FAC, ov);
            let (hi, lo) = l_extract(st.l_ch_noise[i], ov);
            let ltmp1 = mpy_32_16(hi, lo, ONE_MINUS_CNE_SM_FAC, ov);
            st.l_ch_noise[i] = l_add(ltmp, ltmp1, ov);
            if st.l_ch_noise[i] <= MIN_NOISE_ENRG_0 {
                st.l_ch_noise[i] = MIN_NOISE_ENRG_0;
            }
        }
    }
    ivad
}

/// `LTP_flag_update`: whether the open-loop pitch gain of the frame exceeds a mode threshold.
pub(crate) fn ltp_flag_update(st: &mut Vad2State, mode: CodecMode, ov: &mut Flag) {
    let thresh = match mode {
        CodecMode::MR475 | CodecMode::MR515 => 18022,
        CodecMode::MR102 => 19660,
        _ => 21299,
    };
    let (hi, lo) = l_extract(st.l_r0, ov);
    let ltmp = mpy_32_16(hi, lo, thresh, ov);
    st.ltp_flag = st.l_rmax > ltmp;
}

#[cfg(all(test, oracle))]
pub(crate) mod vs_ref;
