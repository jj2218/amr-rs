//! Synthesis of one subframe (`synthesis_amr_wb.cpp`): the 12.8 kHz LP synthesis in double
//! precision (`Syn_filt_32`), de-emphasis (`deemphasis_32`), 50 Hz high-pass, oversampling to
//! 16 kHz (`oversamp_12k8_to_16k`), and the 6.4–7 kHz band from shaped noise (`wb_syn_filt`,
//! `highpass_400Hz_at_12k8`, `band_pass_6k_7k`, `low_pass_filt_7k`).
//!
//! The C's sums of products (`fxp_mac_16by16`) wrap, so their order is free; each filter here
//! computes one output at a time instead of the C's unrolled blocks.

use crate::basic_op::*;
use crate::cnst::{L_SUBFR, L_SUBFR16K, M, M16K, NBBITS_7K, NBBITS_24K, PREEMPH_FAC};
use crate::dtx::SPEECH;
use crate::excitation::{dot, noise_gen_amrwb};
use crate::lpc::{isf_extrapolation, isp_az, weight_amrwb_lpc};
use crate::tables::{FIR_6K_7K, FIR_7K, FIR_UP, HP_GAIN};

const L_FILT: usize = 12;
const L_FIR: usize = 30;

/// Synthesis filter memories (part of `Decoder_State`).
#[derive(Clone)]
pub(crate) struct SynthState {
    pub mem_syn_hi: [i16; M],
    pub mem_syn_lo: [i16; M],
    pub mem_deemph: i16,
    pub mem_sig_out: [i16; 6],
    pub mem_oversamp: [i16; 2 * L_FILT],
    pub mem_syn_hf: [i16; M16K],
    pub mem_hf: [i16; L_FIR],
    pub mem_hf3: [i16; L_FIR],
    pub mem_hp400: [i16; 6],
    pub seed2: i16,
}

impl Default for SynthState {
    fn default() -> Self {
        Self {
            mem_syn_hi: [0; M],
            mem_syn_lo: [0; M],
            mem_deemph: 0,
            mem_sig_out: [0; 6],
            mem_oversamp: [0; 2 * L_FILT],
            mem_syn_hf: [0; M16K],
            mem_hf: [0; L_FIR],
            mem_hf3: [0; L_FIR],
            mem_hp400: [0; 6],
            seed2: 0,
        }
    }
}

/// `Syn_filt_32`: 1/A(z) on `exc` (`>> q_new`) into (hi, lo) halves; `hi`/`lo` hold `m` past
/// samples, then the `lg` new ones. As TS 26.173: the input and output scaling follow `a[0]`
/// (which `Isp_Az`'s adaptive scaling lowers for comfort noise), and the sums saturate.
/// opencore assumed `a[0]` = 4096 and summed plainly (louder comfort noise after rescaling).
///
/// Fast path (exact): with `|lo| <= Lmax` (4095 for every output, `lo` is the low 12 bits)
/// and `|hi|, |exc| <= 32768`, no partial sum of the chains reaches
/// `2·Lmax·A/4096 + 1 + 65536·(|a0| + A)` (A = Σ|a[1..=16]|); below 2^31 the chains are
/// plain dot products.
#[allow(clippy::needless_range_loop)] // i is the sample: it indexes exc and the outputs
pub(crate) fn syn_filt_32(
    a: &[i16],
    m: usize,
    exc: &[i16],
    q_new: i16,
    hi: &mut [i16],
    lo: &mut [i16],
    lg: usize,
) {
    let s = norm_s(a[0]) - 2;
    let a0 = shl_int16(a[0], -(4 + q_new));
    if m == M {
        let ar: [i16; M] = std::array::from_fn(|k| a[M - k]);
        let sum_a: i64 = ar.iter().map(|&c| (c as i64).abs()).sum();
        let lmax = lo[..M]
            .iter()
            .map(|&v| (v as i64).abs())
            .max()
            .unwrap()
            .max(4095);
        let bound = 2 * lmax * sum_a / 4096 + 1 + 65536 * ((a0 as i64).abs() + sum_a);
        if bound < 1 << 31 {
            for i in 0..lg {
                let n = M + i;
                let lw: &[i16; M] = lo[i..n].try_into().unwrap();
                let hw: &[i16; M] = hi[i..n].try_into().unwrap();
                let mut l = dot(lw, &ar).wrapping_shl(1).wrapping_neg() >> 12;
                l += 2 * exc[i] as i32 * a0 as i32;
                l -= dot(hw, &ar).wrapping_shl(1);
                l = shl_int32(l, 3 + s);
                hi[n] = extract_h(l);
                lo[n] = ((l >> 4) - ((hi[n] as i32) << 12)) as i16;
            }
            return;
        }
    }
    for i in 0..lg {
        let n = m + i;
        let mut l = 0i32;
        for j in 1..=m {
            l = msu_16by16_from_int32(l, lo[n - j], a[j]);
        }
        l >>= 12;
        l = mac_16by16_to_int32(l, exc[i], a0);
        for j in 1..=m {
            l = msu_16by16_from_int32(l, hi[n - j], a[j]);
        }
        l = shl_int32(l, 3 + s);
        hi[n] = extract_h(l);
        lo[n] = msu_16by16_from_int32(l >> 4, hi[n], 2048) as i16;
    }
}

/// `wb_syn_filt` (in place, `update` = 1): 1/A(z) of order `m` = `a.len() - 1` on `xy`. As
/// TS 26.173's `Syn_filt`: scaled by `a[0]` (lowered by `Isp_Az`'s adaptive scaling for comfort
/// noise) and saturating; opencore assumed `a[0]` = 4096 and summed plainly.
///
/// Fast path (exact): with 16-bit signals no partial sum reaches `65536·(|a0| + Σ|a[1..]|)`;
/// below 2^31 the chain is a plain dot product.
pub(crate) fn wb_syn_filt(a: &[i16], xy: &mut [i16], mem: &mut [i16]) {
    let m = a.len() - 1;
    let lg = xy.len();
    let s = norm_s(a[0]) - 2;
    let a0 = a[0] >> 1;
    let mut yy = [0i16; M16K + L_SUBFR16K];
    yy[..m].copy_from_slice(mem);
    let sum_a: i64 = a[1..].iter().map(|&c| (c as i64).abs()).sum();
    if 65536 * ((a0 as i64).abs() + sum_a) < 1 << 31 {
        match m {
            M => wb_syn_fast::<M>(a, a0, s, xy, &mut yy),
            M16K => wb_syn_fast::<M16K>(a, a0, s, xy, &mut yy),
            _ => unreachable!("order {m}"),
        }
        mem.copy_from_slice(&xy[lg - m..]);
        return;
    }
    for n in 0..lg {
        let mut l = mul_16by16_to_int32(xy[n], a0);
        for j in 1..=m {
            l = msu_16by16_from_int32(l, a[j], yy[m + n - j]);
        }
        let y = amr_wb_round(shl_int32(l, 3 + s));
        yy[m + n] = y;
        xy[n] = y;
    }
    mem.copy_from_slice(&xy[lg - m..]);
}

/// `wb_syn_filt`'s fast path at order `N` (`yy` holds the memory, then the outputs).
fn wb_syn_fast<const N: usize>(a: &[i16], a0: i16, s: i16, xy: &mut [i16], yy: &mut [i16]) {
    let ar: [i16; N] = std::array::from_fn(|k| a[N - k]);
    for n in 0..xy.len() {
        let w: &[i16; N] = yy[n..n + N].try_into().unwrap();
        let l = 2 * xy[n] as i32 * a0 as i32 - dot(w, &ar).wrapping_shl(1);
        let y = amr_wb_round(shl_int32(l, 3 + s));
        yy[N + n] = y;
        xy[n] = y;
    }
}

/// `deemphasis_32`: `y[i] = x[i] + mu·y[i-1]` on the (hi, lo) signal. Both additions
/// saturate as in TS 26.173 (`L_mac`); opencore adds plainly and wraps on loud signals.
pub(crate) fn deemphasis_32(x_hi: &[i16], x_lo: &[i16], y: &mut [i16], mu: i16, mem: &mut i16) {
    let mut prev = *mem;
    for i in 0..y.len() {
        let mut l = mac_16by16_to_int32((x_hi[i] as i32) << 16, x_lo[i], 8);
        l = shl_int32(l, 3);
        l = mac_16by16_to_int32(l, prev, mu >> 1);
        l = shl_int32(l, 1);
        y[i] = amr_wb_round(l);
        prev = y[i];
    }
    *mem = prev;
}

/// `highpass_50Hz_at_12k8`: 2nd-order IIR, double-precision feedback, in place. The output
/// sum saturates as in TS 26.173 (an `L_mac` chain, then `L_shl`); opencore wraps there (plain
/// products, `<< 2`), which conformance vectors expose.
pub(crate) fn highpass_50hz_at_12k8(signal: &mut [i16], mem: &mut [i16; 6]) {
    const A1: i16 = 16211;
    const A2: i16 = -8021;
    const B: [i16; 3] = [4053, -8106, 4053];
    let [mut y2_hi, mut y2_lo, mut y1_hi, mut y1_lo, mut x0, mut x1] = *mem;
    for s in signal.iter_mut() {
        let x2 = x1;
        x1 = x0;
        x0 = *s;
        // low halves first (cannot overflow), then the high halves and the input taps
        let mut l = fxp_mac_16by16(y2_lo, A2, fxp_mac_16by16(y1_lo, A1, 8192)) >> 14;
        for (a, b) in [(y1_hi, A1), (y2_hi, A2), (x0, B[0]), (x1, B[1]), (x2, B[2])] {
            l = mac_16by16_to_int32(l, a, b);
        }
        let l = shl_int32(l, 2);
        y2_hi = y1_hi;
        y2_lo = y1_lo;
        (y1_hi, y1_lo) = int32_to_dpf(l);
        *s = amr_wb_shl1_round(l);
    }
    *mem = [y2_hi, y2_lo, y1_hi, y1_lo, x0, x1];
}

/// `highpass_400Hz_at_12k8`, in place.
pub(crate) fn highpass_400hz_at_12k8(signal: &mut [i16], mem: &mut [i16; 6]) {
    let [mut y2_hi, mut y2_lo, mut y1_hi, mut y1_lo, mut x0, mut x1] = *mem;
    for s in signal.iter_mut() {
        let mut l1 = fxp_mac_16by16(y1_lo, 29280, 8192);
        let mut l2 = fxp_mul_16by16(y1_hi, 29280);
        l1 = fxp_mac_16by16(y2_lo, -14160, l1);
        l2 = fxp_mac_16by16(y2_hi, -14160, l2);
        let x2 = x1;
        x1 = x0;
        x0 = *s;
        l2 = fxp_mac_16by16(x2, 915, l2);
        l2 = fxp_mac_16by16(x1, -1830, l2);
        l2 = fxp_mac_16by16(x0, 915, l2);
        l1 = (l1 >> 13).wrapping_add(l2.wrapping_shl(2));
        y2_hi = y1_hi;
        y2_lo = y1_lo;
        *s = (l1.wrapping_add(0x8000) >> 16) as i16;
        y1_hi = (l1 >> 16) as i16;
        y1_lo = (l1.wrapping_sub((y1_hi as i32) << 16) >> 1) as i16;
    }
    *mem = [y2_hi, y2_lo, y1_hi, y1_lo, x0, x1];
}

/// `oversamp_12k8_to_16k` (`AmrWbUp_samp`, `AmrWbInterpol`): 64 samples at 12.8 kHz to 80 at
/// 16 kHz; every fifth output copies an input, the others are 24-tap interpolations.
pub(crate) fn oversamp_12k8_to_16k(
    sig12k8: &[i16; L_SUBFR],
    sig16k: &mut [i16],
    mem: &mut [i16; 2 * L_FILT],
) {
    let mut signal = [0i16; 2 * L_FILT + L_SUBFR];
    signal[..2 * L_FILT].copy_from_slice(mem);
    signal[2 * L_FILT..].copy_from_slice(sig12k8);
    for (j, out) in sig16k[..L_SUBFR16K].iter_mut().enumerate() {
        let i = (j * 6554) >> 13;
        let phase = j % 5; // the C's `frac` counts down 5, 4, 3, 2, 1
        *out = if phase == 0 {
            signal[L_FILT + i]
        } else {
            let fir = &FIR_UP[(phase - 1) * 24..phase * 24];
            let x = &signal[L_FILT + i - 11..L_FILT + i + 13];
            let mut sum = 0x2000i32;
            for (&a, &c) in x.iter().zip(fir) {
                sum = fxp_mac_16by16(a, c, sum);
            }
            (shl_int32(sum, 2) >> 16) as i16
        };
    }
    mem.copy_from_slice(&signal[L_SUBFR..L_SUBFR + 2 * L_FILT]);
}

/// A 31-tap FIR over `x` (30 past samples, then 80 new) into `out`: TS 26.173's `L_mac`
/// chain and `round`. Fast path (exact): while `2·max|x|·Σ|c| < 2^31` no partial sum
/// saturates and the chain is the plain dot product. opencore summed plainly (wrapping on loud
/// input) and, at 7 kHz, pre-added the end taps' samples in 16 bits.
/// `x` ends in one zero of slack: the fast path runs 32 taps (the last coefficient 0), one
/// whole vector width.
fn fir31(x: &[i16; L_FIR + L_SUBFR16K + 1], fir: &[i16; L_FIR + 1], out: &mut [i16; L_SUBFR16K]) {
    let abs_sum: i64 = fir.iter().map(|&c| (c as i64).abs()).sum();
    let max = x.iter().map(|&v| (v as i64).abs()).max().unwrap();
    let fir32: [i16; L_FIR + 2] = std::array::from_fn(|k| if k <= L_FIR { fir[k] } else { 0 });
    let fast = 2 * max * abs_sum < 1 << 31;
    for (n, o) in out.iter_mut().enumerate() {
        let l = if fast {
            let w: &[i16; L_FIR + 2] = x[n..n + L_FIR + 2].try_into().unwrap();
            dot(w, &fir32).wrapping_shl(1)
        } else {
            let w = &x[n..n + L_FIR + 1];
            w.iter()
                .zip(fir)
                .fold(0, |l, (&a, &c)| mac_16by16_to_int32(l, a, c))
        };
        *o = amr_wb_round(l);
    }
}

/// `band_pass_6k_7k`: 31-tap FIR on `signal >> 2`, in place.
pub(crate) fn band_pass_6k_7k(signal: &mut [i16; L_SUBFR16K], mem: &mut [i16; L_FIR]) {
    // the table holds taps 0..=29; the filter is symmetric (tap 30 = tap 0)
    let fir: [i16; L_FIR + 1] = std::array::from_fn(|k| FIR_6K_7K[k.min(L_FIR - k)]);
    let mut x = [0i16; L_FIR + L_SUBFR16K + 1];
    x[..L_FIR].copy_from_slice(mem);
    for (d, &s) in x[L_FIR..].iter_mut().zip(signal.iter()) {
        *d = s >> 2;
    }
    fir31(&x, &fir, signal);
    mem.copy_from_slice(&x[L_SUBFR16K..L_SUBFR16K + L_FIR]);
}

/// `low_pass_filt_7k`: 31-tap FIR, in place (23.85 kbit/s).
pub(crate) fn low_pass_filt_7k(signal: &mut [i16; L_SUBFR16K], mem: &mut [i16; L_FIR]) {
    let mut x = [0i16; L_FIR + L_SUBFR16K + 1];
    x[..L_FIR].copy_from_slice(mem);
    x[L_FIR..L_FIR + L_SUBFR16K].copy_from_slice(signal);
    fir31(&x, &FIR_7K, signal);
    mem.copy_from_slice(&x[L_SUBFR16K..L_SUBFR16K + L_FIR]);
}

/// `synthesis_amr_wb`: one subframe from the excitation `exc` (12.8 kHz, modified in place as
/// the C does) to 80 samples at 16 kHz in `synth16k`.
#[allow(clippy::too_many_arguments)] // the C's signature
pub(crate) fn synthesis_amr_wb(
    aq: &[i16],
    exc: &mut [i16],
    q_new: i16,
    synth16k: &mut [i16],
    prms: i16,
    hf_isf: &mut [i16; M16K],
    nb_bits: i16,
    new_dtx_state: i16,
    st: &mut SynthState,
    vad_hist: i16,
    bfi: bool,
) {
    let mut synth_hi = [0i16; M + L_SUBFR];
    let mut synth_lo = [0i16; M + L_SUBFR];
    let mut synth = [0i16; L_SUBFR];
    let mut hf = [0i16; L_SUBFR16K];
    synth_hi[..M].copy_from_slice(&st.mem_syn_hi);
    synth_lo[..M].copy_from_slice(&st.mem_syn_lo);
    syn_filt_32(aq, M, exc, q_new, &mut synth_hi, &mut synth_lo, L_SUBFR);
    st.mem_syn_hi.copy_from_slice(&synth_hi[L_SUBFR..]);
    st.mem_syn_lo.copy_from_slice(&synth_lo[L_SUBFR..]);
    deemphasis_32(
        &synth_hi[M..],
        &synth_lo[M..],
        &mut synth,
        PREEMPH_FAC,
        &mut st.mem_deemph,
    );
    highpass_50hz_at_12k8(&mut synth, &mut st.mem_sig_out);
    oversamp_12k8_to_16k(&synth, synth16k, &mut st.mem_oversamp);

    // high band: white noise scaled to the excitation energy
    for v in hf.iter_mut() {
        *v = noise_gen_amrwb(&mut st.seed2) >> 3;
    }
    for v in exc[..L_SUBFR].iter_mut() {
        *v = crate::excitation::shr3_round(*v);
    }
    let q_new = q_new.wrapping_sub(3);
    let (mut exp_ener, mut exp) = (0i16, 0i16);
    let ener = extract_h(dot_product12(
        &exc[..L_SUBFR],
        &exc[..L_SUBFR],
        &mut exp_ener,
    ));
    exp_ener = exp_ener.wrapping_sub(q_new.wrapping_shl(1));
    let mut tmp = extract_h(dot_product12(&hf, &hf, &mut exp));
    if tmp > ener {
        tmp >>= 1;
        exp += 1;
    }
    let mut l = l_deposit_h(div_16by16(tmp, ener));
    exp = exp.wrapping_sub(exp_ener);
    one_ov_sqrt_norm(&mut l, &mut exp);
    let l = shl_int32(l, exp.wrapping_add(1));
    let tmp = (l >> 16) as i16;
    for v in hf.iter_mut() {
        *v = (fxp_mul_16by16(*v, tmp) >> 15) as i16;
    }

    // tilt of the 12.8 kHz synthesis above 400 Hz sets the high-band gain
    highpass_400hz_at_12k8(&mut synth, &mut st.mem_hp400);
    let mut l = mac_16by16_to_int32(1, synth[0], synth[0]);
    let mut l2 = 1i32;
    for i in 1..L_SUBFR {
        l = mac_16by16_to_int32(l, synth[i], synth[i]);
        l2 = mac_16by16_to_int32(l2, synth[i], synth[i - 1]);
    }
    let exp = normalize_amr_wb(l);
    let ener = (l.wrapping_shl(exp as u32) >> 16) as i16;
    let tmp = (l2.wrapping_shl(exp as u32) >> 16) as i16;
    let fac = if tmp > 0 { div_16by16(tmp, ener) } else { 0 };
    let gain1 = 32767 - fac;
    let gain2 = shl_int16(mult_int16(gain1, 20480), 1);
    let mut tmp = if vad_hist > 0 {
        gain2.wrapping_sub(1)
    } else {
        gain1.wrapping_sub(1)
    };
    if tmp != 0 {
        tmp = tmp.wrapping_add(1);
    }
    if tmp < 3277 {
        tmp = 3277;
    }
    if nb_bits >= NBBITS_24K && !bfi {
        let g = HP_GAIN[prms as usize];
        for v in hf.iter_mut() {
            *v = mult_int16(*v, g).wrapping_shl(1);
        }
    } else {
        for v in hf.iter_mut() {
            *v = mult_int16(*v, tmp);
        }
    }

    let mut ap = [0i16; M16K + 1];
    if nb_bits <= NBBITS_7K && new_dtx_state == SPEECH {
        let mut hf_a = [0i16; M16K + 1];
        isf_extrapolation(hf_isf);
        isp_az(hf_isf, &mut hf_a, false);
        weight_amrwb_lpc(&hf_a, &mut ap, 29491);
        wb_syn_filt(&ap, &mut hf, &mut st.mem_syn_hf);
    } else {
        weight_amrwb_lpc(&aq[..M + 1], &mut ap[..M + 1], 19661);
        wb_syn_filt(&ap[..M + 1], &mut hf, &mut st.mem_syn_hf[M16K - M..]);
    }
    band_pass_6k_7k(&mut hf, &mut st.mem_hf);
    if nb_bits >= NBBITS_24K {
        low_pass_filt_7k(&mut hf, &mut st.mem_hf3);
    }
    for (s, &h) in synth16k[..L_SUBFR16K].iter_mut().zip(&hf) {
        *s = add_int16(*s, h);
    }
}
