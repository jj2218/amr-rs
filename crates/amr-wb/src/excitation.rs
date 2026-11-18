//! Excitation: adaptive codebook (`pred_lt4.cpp`, `lagconceal.cpp`), algebraic codebook
//! (`dec_acelp_*_in_64.cpp`, `dec_alg_codebook.cpp`), and the per-subframe helpers
//! (`preemph_amrwb_dec`, `pit_shrp`, `voice_factor`, `phase_dispersion`, `noise_gen_amrwb`,
//! `scale_signal`, `agc2_amr_wb`, `median5`).

use crate::basic_op::*;
use crate::cnst::L_SUBFR;
use crate::tables::{INTER4_2, PH_IMP_LOW, PH_IMP_MID};

const UP_SAMP: i16 = 4;
const L_INTERPOL2: usize = 16;

/// `Pred_lt4`: adaptive-codebook excitation `buf[pos..pos + l_subfr]` from the past
/// excitation at lag `t0` + `frac`/4, interpolated with a 32-tap FIR. In place: later samples
/// read earlier outputs (`t0` >= 34). Saturating as TS 26.173 (`L_mac` chain, `L_shl`,
/// `round`); opencore sums plain products and truncates `>> 14` to 16 bits, wrapping on loud
/// excitation (conformance vectors expose it).
///
/// Computed as clamp16((Σ c·x + 2^13) >> 14) (plain sums, `simd::fir32`'s form), exact for
/// every input by `CHAIN_IS_CLAMPED_SUM` (below): the chain's L_shl + round equal it.
#[inline(always)] // into simd::decode's AVX2 copy
pub(crate) fn pred_lt4(buf: &mut [i16], pos: usize, t0: i16, frac: i16, l_subfr: usize) {
    let mut start = pos as isize - t0 as isize;
    let mut frac = -frac;
    if frac < 0 {
        frac += UP_SAMP;
        start -= 1;
    }
    start -= L_INTERPOL2 as isize - 1;
    let row = (UP_SAMP - 1 - frac) as usize;
    let fir: &[i16; 2 * L_INTERPOL2] = INTER4_2[row * 2 * L_INTERPOL2..][..2 * L_INTERPOL2]
        .try_into()
        .unwrap();
    let start = start as usize;
    // Blocks of outputs on `fir32`'s tap pairs: a block of B outputs at n0 reads
    // buf[start + n0..start + n0 + B + 31], all before pos + n0 when pos - start >= B + 31
    // (pos - start = t0 + 15 or 16 >= 49): no output reads one of its own block.
    let mut n0 = 0;
    let d = pos - start;
    loop {
        let (src, dst) = buf.split_at_mut(pos + n0);
        let x = &src[start + n0..];
        if d >= 63 && n0 + 32 <= l_subfr {
            crate::simd::fir32::<32, 14>(x, fir, (&mut dst[..32]).try_into().unwrap());
            n0 += 32;
        } else if d >= 47 && n0 + 16 <= l_subfr {
            crate::simd::fir32::<16, 14>(x, fir, (&mut dst[..16]).try_into().unwrap());
            n0 += 16;
        } else {
            break;
        }
    }
    for n in n0..l_subfr {
        let x: &[i16; 2 * L_INTERPOL2] = buf[start + n..][..2 * L_INTERPOL2].try_into().unwrap();
        buf[pos + n] = (dot(x, fir).wrapping_add(1 << 13) >> 14).clamp(-32768, 32767) as i16;
    }
}

/// `Pred_lt4`'s chain (`L_mac` per tap from 0, `L_shl(·, 1)`, `round`) is the clamped plain
/// sum for every input when each `INTER4_2` row has `Σ|c| < 32768 + 16384`. Plain: S = Σ 2·c·x,
/// |S| <= 65536·Σ|c| < 2^32 - 2^14, so S/2 + 2^13 fits an i32. Without saturation the chain is S
/// and round(L_shl(S, 1)) = clamp16((S/2 + 2^13) >> 14). A chain first saturates (to ±2^31)
/// at tap k only if 65536·Σ_{j<=k}|c_j| >= 2^31, i.e. that prefix >= 32768; the taps after k
/// then sum below 16384 in |c|, so they move the chain, and S beyond the prefix, by under
/// 65536·16384 = 2^30: both stay at or beyond ±2^30 with the prefix's sign, where both forms
/// give 32767 (S >= 2^30) or -32768 (S <= -2^30).
const CHAIN_IS_CLAMPED_SUM: () = {
    let mut r = 0;
    while r < 4 {
        let (mut s, mut k) = (0i32, 0);
        while k < 2 * L_INTERPOL2 {
            s += (INTER4_2[r * 2 * L_INTERPOL2 + k] as i32).abs();
            k += 1;
        }
        assert!(
            s < 32768 + 16384,
            "an INTER4_2 row's chain may saturate unclamped"
        );
        r += 1;
    }
};
const _: () = CHAIN_IS_CLAMPED_SUM;

/// Plain `Σ x·y` (wrapping; callers bound it below 2^31). Fixed length: one vector loop.
#[inline(always)]
pub(crate) fn dot<const N: usize>(x: &[i16; N], y: &[i16; N]) -> i32 {
    x.iter()
        .zip(y)
        .fold(0i32, |s, (&a, &b)| s.wrapping_add(a as i32 * b as i32))
}

const NB_POS: i16 = 16;

fn dec_1p_n1(index: i32, n: i16, offset: i16, pos: &mut [i16]) {
    let mask = (1i32 << n) - 1;
    let mut pos1 = ((index & mask) + offset as i32) as i16;
    if (index >> n) & 1 == 1 {
        pos1 = pos1.wrapping_add(NB_POS);
    }
    pos[0] = pos1;
}

fn dec_2p_2n1(index: i32, n: i16, offset: i16, pos: &mut [i16]) {
    let mask = sub_int16(shl_int16(1, n), 1) as i32;
    let mut pos1 = add_int32(shr_int32(index, n) & mask, offset as i32) as i16;
    let tmp = shl_int16(n, 1);
    let i = (index >> tmp) & 1;
    let mut pos2 = add_int16((index & mask) as i16, offset);
    if pos2 < pos1 {
        if i == 1 {
            pos1 = pos1.wrapping_add(NB_POS);
        } else {
            pos2 = pos2.wrapping_add(NB_POS);
        }
    } else if i == 1 {
        pos1 = pos1.wrapping_add(NB_POS);
        pos2 = pos2.wrapping_add(NB_POS);
    }
    pos[0] = pos1;
    pos[1] = pos2;
}

fn dec_3p_3n1(index: i32, n: i16, offset: i16, pos: &mut [i16]) {
    let mask = (1i32 << (2 * n - 1)) - 1;
    let idx = index & mask;
    let mut j = offset;
    if (index >> ((n << 1) - 1)) & 1 != 0 {
        j = j.wrapping_add(1 << (n - 1));
    }
    dec_2p_2n1(idx, n - 1, j, pos);
    let mask = (1i32 << (n + 1)) - 1;
    let idx = (index >> (n << 1)) & mask;
    dec_1p_n1(idx, n, offset, &mut pos[2..]);
}

fn dec_4p_4n1(index: i32, n: i16, offset: i16, pos: &mut [i16]) {
    let tmp = (n << 1) - 1;
    let mask = (1i32 << tmp) - 1;
    let idx = index & mask;
    let mut j = offset;
    if (index >> tmp) & 1 != 0 {
        j = j.wrapping_add(1 << (n - 1));
    }
    dec_2p_2n1(idx, n - 1, j, pos);
    let mask = (1i32 << ((n << 1) + 1)) - 1;
    let idx = (index >> (n << 1)) & mask;
    dec_2p_2n1(idx, n, offset, &mut pos[2..]);
}

fn dec_4p_4n(index: i32, n: i16, offset: i16, pos: &mut [i16]) {
    let n_1 = n - 1;
    let j = offset.wrapping_add(1 << n_1);
    match (index >> ((n << 2) - 2)) & 3 {
        0 => {
            if (index >> ((n_1 << 2) + 1)) & 1 != 0 {
                dec_4p_4n1(index, n_1, j, pos);
            } else {
                dec_4p_4n1(index, n_1, offset, pos);
            }
        }
        1 => {
            dec_1p_n1(index >> (3 * n_1 + 1), n_1, offset, pos);
            dec_3p_3n1(index, n_1, j, &mut pos[1..]);
        }
        2 => {
            dec_2p_2n1(index >> ((n_1 << 1) + 1), n_1, offset, pos);
            dec_2p_2n1(index, n_1, j, &mut pos[2..]);
        }
        _ => {
            dec_3p_3n1(index >> (n_1 + 1), n_1, offset, pos);
            dec_1p_n1(index, n_1, j, &mut pos[3..]);
        }
    }
}

fn dec_5p_5n(index: i32, n: i16, offset: i16, pos: &mut [i16]) {
    let n_1 = n - 1;
    let j = add_int16(offset, shl_int16(1, n_1));
    let idx = index >> ((n << 1) + 1);
    if (index >> (5 * n - 1)) & 1 != 0 {
        dec_3p_3n1(idx, n_1, j, pos);
    } else {
        dec_3p_3n1(idx, n_1, offset, pos);
    }
    dec_2p_2n1(index, n, offset, &mut pos[3..]);
}

fn dec_6p_6n_2(index: i32, n: i16, offset: i16, pos: &mut [i16]) {
    let n_1 = n - 1;
    let j = offset.wrapping_add(1 << n_1);
    let (mut offset_a, mut offset_b) = (j, j);
    if (index >> (6 * n - 5)) & 1 == 0 {
        offset_a = offset;
    } else {
        offset_b = offset;
    }
    match (index >> (6 * n - 4)) & 3 {
        0 => {
            dec_5p_5n(index >> n, n_1, offset_a, pos);
            dec_1p_n1(index, n_1, offset_a, &mut pos[5..]);
        }
        1 => {
            dec_5p_5n(index >> n, n_1, offset_a, pos);
            dec_1p_n1(index, n_1, offset_b, &mut pos[5..]);
        }
        2 => {
            dec_4p_4n(index >> (2 * n_1 + 1), n_1, offset_a, pos);
            dec_2p_2n1(index, n_1, offset_b, &mut pos[4..]);
        }
        _ => {
            dec_3p_3n1(index >> (3 * n_1 + 1), n_1, offset, pos);
            dec_3p_3n1(index, n_1, j, &mut pos[3..]);
        }
    }
}

/// `dec_acelp_2p_in_64`: 12-bit algebraic codebook, 2 pulses (6.60 kbit/s).
pub(crate) fn dec_acelp_2p_in_64(index: i16, code: &mut [i16; L_SUBFR]) {
    code.fill(0);
    let i = ((index >> 5) & 0x003E) as usize;
    code[i] = if (index >> 6) & 32 == 0 { 512 } else { -512 };
    let i = (((index & 0x001F) << 1) + 1) as usize;
    code[i] = if index & 32 == 0 { 512 } else { -512 };
}

fn add_pulses(pos: &[i16], track: usize, code: &mut [i16; L_SUBFR]) {
    for &p in pos {
        let i = (((p & (NB_POS - 1)) << 2) as usize) + track;
        code[i] = if p & NB_POS == 0 {
            code[i].wrapping_add(512)
        } else {
            code[i].wrapping_sub(512)
        };
    }
}

/// `dec_acelp_4p_in_64`: algebraic codebook, 4 tracks, `nbbits` 20..88 (8.85..23.85 kbit/s).
#[allow(clippy::needless_range_loop)] // k is the track, as in the C
pub(crate) fn dec_acelp_4p_in_64(index: &[i16; 8], nbbits: i16, code: &mut [i16; L_SUBFR]) {
    code.fill(0);
    let mut pos = [0i16; 6];
    let joined = |k: usize, s: u32| ((index[k] as i32) << s) + index[k + 4] as i32;
    match nbbits {
        20 => {
            for k in 0..4 {
                dec_1p_n1(index[k] as i32, 4, 0, &mut pos);
                add_pulses(&pos[..1], k, code);
            }
        }
        36 => {
            for k in 0..4 {
                dec_2p_2n1(index[k] as i32, 4, 0, &mut pos);
                add_pulses(&pos[..2], k, code);
            }
        }
        44 => {
            for k in 0..2 {
                dec_3p_3n1(index[k] as i32, 4, 0, &mut pos);
                add_pulses(&pos[..3], k, code);
            }
            for k in 2..4 {
                dec_2p_2n1(index[k] as i32, 4, 0, &mut pos);
                add_pulses(&pos[..2], k, code);
            }
        }
        52 => {
            for k in 0..4 {
                dec_3p_3n1(index[k] as i32, 4, 0, &mut pos);
                add_pulses(&pos[..3], k, code);
            }
        }
        64 => {
            for k in 0..4 {
                dec_4p_4n(joined(k, 14), 4, 0, &mut pos);
                add_pulses(&pos[..4], k, code);
            }
        }
        72 => {
            for k in 0..2 {
                dec_5p_5n(joined(k, 10), 4, 0, &mut pos);
                add_pulses(&pos[..5], k, code);
            }
            for k in 2..4 {
                dec_4p_4n(joined(k, 14), 4, 0, &mut pos);
                add_pulses(&pos[..4], k, code);
            }
        }
        88 => {
            for k in 0..4 {
                dec_6p_6n_2(joined(k, 11), 4, 0, &mut pos);
                add_pulses(&pos[..6], k, code);
            }
        }
        _ => {}
    }
}

/// `preemph_amrwb_dec`: `x[i] -= mu·x[i-1]`, backwards, without memory.
pub(crate) fn preemph_amrwb_dec(x: &mut [i16], mu: i16) {
    for i in (1..x.len()).rev() {
        x[i] = amr_wb_round(msu_16by16_from_int32((x[i] as i32) << 16, x[i - 1], mu));
    }
}

/// `Pit_shrp`: `x[i] += sharp·x[i - pit_lag]`.
pub(crate) fn pit_shrp(x: &mut [i16], pit_lag: i16, sharp: i16) {
    for i in pit_lag.max(0) as usize..x.len() {
        let l = mac_16by16_to_int32((x[i] as i32) << 16, x[i - pit_lag as usize], sharp);
        x[i] = amr_wb_round(l);
    }
}

/// `noise_gen_amrwb`: 16-bit linear congruential generator.
#[inline(always)]
pub(crate) fn noise_gen_amrwb(seed: &mut i16) -> i16 {
    *seed = fxp_mac_16by16(*seed, 31821, 13849) as i16;
    *seed
}

/// `voice_factor`: voicing in [-1, 1] (Q15) from the pitch and code energies.
pub(crate) fn voice_factor(
    exc: &[i16],
    q_exc: i16,
    gain_pit: i16,
    code: &[i16],
    gain_code: i16,
) -> i16 {
    let (mut exp1, mut exp2) = (0i16, 0i16);
    let mut ener1 = extract_h(dot_product12(exc, exc, &mut exp1));
    exp1 = sub_int16(exp1, q_exc.wrapping_shl(1));
    let l = mul_16by16_to_int32(gain_pit, gain_pit);
    let exp = normalize_amr_wb(l);
    let tmp = (l.wrapping_shl(exp as u32) >> 16) as i16;
    ener1 = mult_int16(ener1, tmp);
    exp1 = exp1.wrapping_sub(exp + 10);
    let mut ener2 = extract_h(dot_product12(code, code, &mut exp2));
    let exp = norm_s(gain_code);
    let tmp = shl_int16(gain_code, exp);
    let tmp = mult_int16(tmp, tmp);
    ener2 = mult_int16(ener2, tmp);
    exp2 = exp2.wrapping_sub(exp << 1);
    let i = exp1.wrapping_sub(exp2) as i32;
    // TS 26.173's `shr`: counts past 15 give 0 (opencore's `>>=` is masked by x86 to 5 bits,
    // so a zero pitch excitation, 40 apart in exponent, kept energy)
    let shr = |v: i16, n: i32| v >> n.min(15);
    if i >= 0 {
        ener1 >>= 1;
        ener2 = shr(ener2, i + 1);
    } else {
        ener1 = shr(ener1, 1 - i);
        ener2 >>= 1;
    }
    let tmp = ener1.wrapping_sub(ener2);
    let ener1 = ener1.wrapping_add(ener2.wrapping_add(1));
    if tmp >= 0 {
        div_16by16(tmp, ener1)
    } else {
        negate_int16(div_16by16(negate_int16(tmp), ener1))
    }
}

const PITCH_0_9: i16 = 14746;
const PITCH_0_6: i16 = 9830;

/// `phase_dispersion`: smooths sparse codevectors at low pitch gain; `disp_mem` holds the
/// previous state, code gain and six pitch gains.
pub(crate) fn phase_dispersion(
    gain_code: i16,
    gain_pit: i16,
    code: &mut [i16; L_SUBFR],
    mode: i16,
    disp_mem: &mut [i16; 8],
) {
    let mut state: i16 = if gain_pit < PITCH_0_6 {
        0
    } else if gain_pit < PITCH_0_9 {
        1
    } else {
        2
    };
    disp_mem.copy_within(2..7, 3);
    disp_mem[2] = gain_pit;
    let prev_gain_code = disp_mem[1];
    if sub_int16(gain_code, prev_gain_code) > shl_int16(prev_gain_code, 1) {
        if state < 2 {
            state += 1;
        }
    } else {
        let j = disp_mem[2..8].iter().filter(|&&g| g < PITCH_0_6).count();
        if j > 2 {
            state = 0;
        }
        if state as i32 > disp_mem[0] as i32 + 1 {
            state -= 1;
        }
    }
    disp_mem[1] = gain_code;
    disp_mem[0] = state;
    let state = state.wrapping_add(mode);
    let imp = match state {
        0 => &PH_IMP_LOW,
        1 => &PH_IMP_MID,
        _ => return,
    };
    crate::simd::phase_conv(code, imp);
}

// `mult_int16_r` saturates only for (-32768)²; with no -32768 tap it is `pmulhrsw` exactly
// (`simd::phase_conv`'s AVX2 body).
const _: () = {
    let mut j = 0;
    while j < L_SUBFR {
        assert!(PH_IMP_LOW[j] != i16::MIN && PH_IMP_MID[j] != i16::MIN);
        j += 1;
    }
};

/// `phase_dispersion`'s circular convolution of `code` with `imp` (taps without -32768), in
/// the C's order: one saturating row per non-zero sample into `code2`, then the wrapped half
/// added on, saturating. The contract and the oracle of `simd::phase_conv`.
#[inline(always)]
pub(crate) fn phase_conv_portable(code: &mut [i16; L_SUBFR], imp: &[i16; L_SUBFR]) {
    let mut code2 = [0i16; 2 * L_SUBFR];
    for i in 0..L_SUBFR {
        if code[i] != 0 {
            for j in 0..L_SUBFR {
                code2[i + j] = add_int16(code2[i + j], mult_int16_r(code[i], imp[j]));
            }
        }
    }
    for i in 0..L_SUBFR {
        code[i] = add_int16(code2[i], code2[i + L_SUBFR]);
    }
}

const L_LTPHIST: usize = 5;
const ONE_PER_3: i16 = 10923;
const ONE_PER_LTPHIST: i16 = 6554;

/// `Init_Lagconc`.
pub(crate) fn init_lagconc(lag_hist: &mut [i16; L_LTPHIST]) {
    lag_hist.fill(64);
}

/// The C's `insertion_sort` (ascending).
fn insertion_sort(a: &mut [i16; L_LTPHIST]) {
    for n in 0..L_LTPHIST {
        let x = a[n];
        let mut i = n as isize - 1;
        while i >= 0 && x < a[i as usize] {
            a[i as usize + 1] = a[i as usize];
            i -= 1;
        }
        a[(i + 1) as usize] = x;
    }
}

/// `lagconceal`: pitch lag for a bad frame from the lag and gain histories.
pub(crate) fn lagconceal(
    gain_hist: &[i16],
    lag_hist: &[i16; L_LTPHIST],
    t0: &mut i16,
    old_t0: i16,
    seed: &mut i16,
    unusable_frame: bool,
) {
    let last_gain = gain_hist[4];
    let sec_last_gain = gain_hist[3];
    let last_lag = lag_hist[0];
    let min_lag = *lag_hist.iter().min().unwrap();
    let max_lag = *lag_hist.iter().max().unwrap();
    let min_gain = *gain_hist[..L_LTPHIST].iter().min().unwrap();
    let lag_dif = sub_int16(max_lag, min_lag);
    let random_lag = |seed: &mut i16| {
        let mut h = *lag_hist;
        insertion_sort(&mut h);
        let dif = sub_int16(h[4], h[2]).min(40);
        let d = noise_gen_amrwb(seed);
        let d2 = mult_int16(dif >> 1, d);
        let tmp = add_int16(add_int16(h[2], h[3]), h[4]);
        add_int16(mult_int16(tmp, ONE_PER_3), d2)
    };
    let clamp = |t: i16| {
        if t > max_lag {
            max_lag
        } else if t < min_lag {
            min_lag
        } else {
            t
        }
    };
    if unusable_frame {
        *t0 = if min_gain > 8192 && lag_dif < 10 {
            old_t0
        } else if last_gain > 8192 && sec_last_gain > 8192 {
            lag_hist[0]
        } else {
            random_lag(seed)
        };
        *t0 = clamp(*t0);
    } else {
        let mut mean_lag = 0i16;
        for &l in lag_hist {
            mean_lag = add_int16(mean_lag, l);
        }
        let mean_lag = mult_int16(mean_lag, ONE_PER_LTPHIST);
        let t = *t0 as i32;
        let tmp = t - max_lag as i32;
        let tmp2 = t - last_lag as i32;
        let (min_l, max_l) = (min_lag as i32, max_lag as i32);
        let keep = (lag_dif < 10 && t > min_l - 5 && tmp < 5)
            || (last_gain > 8192 && sec_last_gain > 8192 && tmp2 + 10 > 0 && tmp2 < 10)
            || (min_gain < 6554 && last_gain == min_gain && t > min_l && t < max_l)
            || (lag_dif < 70 && t > min_l && t < max_l)
            || (t > mean_lag as i32 && t < max_l);
        if !keep {
            *t0 = if (min_gain > 8192 && lag_dif < 10) || (last_gain > 8192 && sec_last_gain > 8192)
            {
                lag_hist[0]
            } else {
                random_lag(seed)
            };
            *t0 = clamp(*t0);
        }
    }
}

/// `scale_signal`: `x = round(L_shl(x << 16, exp))`, TS 26.173's `Scale_sig`. For `exp < 0`
/// that is `(x + 2^(e-1)) >> e` in 32 bits; opencore saturates the addition in 16 bits (so
/// 32764..=32767 >> 3 gave 4095, not 4096) and masks the count to 4 bits.
pub(crate) fn scale_signal(x: &mut [i16], exp: i16) {
    if exp > 0 {
        for v in x.iter_mut() {
            *v = amr_wb_round(shl_int32((*v as i32) << 16, exp));
        }
    } else if exp < 0 {
        let e = (-(exp as i32)).min(31) as u32;
        for v in x.iter_mut() {
            *v = amr_wb_round(((*v as i32) << 16) >> e);
        }
    }
}

/// `L_shl(l, K)`, saturating, as compares and selects (vectorises; `shl_int32` widens to i64).
#[inline(always)]
fn shl_sat<const K: u32>(l: i32) -> i32 {
    if l > i32::MAX >> K {
        i32::MAX
    } else if l < i32::MIN >> K {
        i32::MIN
    } else {
        l << K
    }
}

/// The decoder's excitation mix per sample, as TS 26.173 (every op saturating):
/// round(L_shl(L_mac(L_shl(L_mult(c, g), 5), e, gp), 1)). Branch-free over fixed arrays, so
/// LLVM runs it in lanes (`simd::decode`'s AVX2 copy too).
#[inline(always)]
pub(crate) fn mix_exc(c: &[i16; L_SUBFR], g: i16, e: &[i16; L_SUBFR], gp: i16) -> [i16; L_SUBFR] {
    let mut out = [0i16; L_SUBFR];
    for ((o, &c), &e) in out.iter_mut().zip(c).zip(e) {
        let l = shl_sat::<5>(mul_16by16_to_int32(c, g));
        *o = amr_wb_round(shl_sat::<1>(l.saturating_add(mul_16by16_to_int32(e, gp))));
    }
    out
}

/// The enhancer's code tilt: round(L_msu(L_msu(c[i] << 16, c[i + 1], t), c[i - 1], t)), the
/// missing neighbours at the ends as 0 (an `L_msu` of 0 is the identity).
#[inline(always)]
pub(crate) fn tilt_code(c: &[i16; L_SUBFR], t: i16) -> [i16; L_SUBFR] {
    let mut p = [0i16; L_SUBFR + 2];
    p[1..=L_SUBFR].copy_from_slice(c);
    let mut out = [0i16; L_SUBFR];
    for (o, w) in out.iter_mut().zip(p.windows(3)) {
        let l = msu_16by16_from_int32((w[1] as i32) << 16, w[2], t);
        *o = amr_wb_round(msu_16by16_from_int32(l, w[0], t));
    }
    out
}

#[cfg(test)]
#[test]
fn mix_exc_and_tilt_code_are_the_chains() {
    let mut r = crate::testutil::Rng::new(11);
    let v = |r: &mut crate::testutil::Rng| match r.u64() % 4 {
        0 => [i16::MIN, i16::MAX, -1, 0][(r.u64() % 4) as usize],
        _ => r.i16() >> (r.u64() % 12),
    };
    for _ in 0..20_000 {
        let c: [i16; L_SUBFR] = std::array::from_fn(|_| v(&mut r));
        let e: [i16; L_SUBFR] = std::array::from_fn(|_| v(&mut r));
        let (g, gp, t) = (v(&mut r), v(&mut r), v(&mut r));
        let want: [i16; L_SUBFR] = std::array::from_fn(|i| {
            let mut l = shl_int32(mul_16by16_to_int32(c[i], g), 5);
            l = mac_16by16_to_int32(l, e[i], gp);
            amr_wb_round(shl_int32(l, 1))
        });
        assert_eq!(mix_exc(&c, g, &e, gp), want, "mix g {g} gp {gp}");
        let want: [i16; L_SUBFR] = std::array::from_fn(|i| {
            let mut l = (c[i] as i32) << 16;
            if i + 1 < L_SUBFR {
                l = msu_16by16_from_int32(l, c[i + 1], t);
            }
            if i > 0 {
                l = msu_16by16_from_int32(l, c[i - 1], t);
            }
            amr_wb_round(l)
        });
        assert_eq!(tilt_code(&c, t), want, "tilt t {t}");
    }
}

/// `Scale_sig(x, -3)` of one value: `(x + 4) >> 3` without 16-bit saturation.
#[inline(always)]
pub(crate) fn shr3_round(x: i16) -> i16 {
    ((x as i32 + 4) >> 3) as i16
}

/// `agc2_amr_wb`: scales `sig_out` to the energy of `sig_in`.
pub(crate) fn agc2_amr_wb(sig_in: &[i16], sig_out: &mut [i16]) {
    // The C's energy chains (2·t² from t = v >> 2, then L_mac of the others): terms >= 0
    // and |t| <= 8192, so each is min(Σ 2·t², MAX_32), a plain sum.
    let energy = |x: &[i16]| {
        let e: u64 = x
            .iter()
            .map(|&v| ((v >> 2) as i32 * (v >> 2) as i32) as u64)
            .sum();
        (2 * e).min(MAX_32 as u64) as i32
    };
    let s = energy(sig_out);
    if s == 0 {
        return;
    }
    let mut exp = normalize_amr_wb(s) - 1;
    // `L_shl` as TS 26.173: exp = -1 (energy >= 2^30) halves s; opencore's `s << exp` became
    // a shift by 31 on x86
    let gain_out = amr_wb_round(shl_int32(s, exp));
    let s = energy(sig_in);
    let g0 = if s == 0 {
        0
    } else {
        let i = normalize_amr_wb(s);
        let gain_in = amr_wb_round(shl_int32(s, i));
        exp -= i;
        let mut s = div_16by16(gain_out, gain_in) as i32;
        s = shl_int32(s, 7);
        s = shr_int32(s, exp);
        s = one_ov_sqrt(s);
        amr_wb_round(shl_int32(s, 9))
    };
    for v in sig_out.iter_mut() {
        // extract_h(L_shl(v·g0, 3)), the saturation as a clamp (each sample on its own)
        let p = (*v as i64 * g0 as i64) << 3;
        *v = (p.clamp(i32::MIN as i64, i32::MAX as i64) >> 16) as i16;
    }
}

/// `median5`: the median of `x[c-2..=c+2]`.
pub(crate) fn median5(x: &[i16], c: usize) -> i16 {
    let (mut x1, mut x2, mut x3, mut x4, mut x5) = (x[c - 2], x[c - 1], x[c], x[c + 1], x[c + 2]);
    if x2 < x1 {
        std::mem::swap(&mut x1, &mut x2);
    }
    if x3 < x1 {
        std::mem::swap(&mut x1, &mut x3);
    }
    if x4 < x1 {
        std::mem::swap(&mut x1, &mut x4);
    }
    if x5 < x1 {
        x5 = x1;
    }
    if x3 < x2 {
        std::mem::swap(&mut x2, &mut x3);
    }
    if x4 < x2 {
        std::mem::swap(&mut x2, &mut x4);
    }
    if x5 < x2 {
        x5 = x2;
    }
    if x4 < x3 {
        x3 = x4;
    }
    if x5 < x3 {
        x3 = x5;
    }
    x3
}
