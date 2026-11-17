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
use crate::excitation::noise_gen_amrwb;
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
/// Fast path (exact): with `|lo| <= Lmax` (4095 for every output, `lo` is the low 12 bits),
/// `|hi| <= H` and `|exc| <= E`, no partial sum of the chains reaches
/// `2·Lmax·A` (low halves) or `2·Lmax·A/4096 + 1 + 2·E·|a0| + 2·H·A` (A = Σ|a[1..=16]|); below
/// 2^31 the chains are plain
/// sums. H is checked after the fast pass, over the history and the outputs it made: the first
/// output that differed from the chains would have read only exact ones, all within H.
pub(crate) fn syn_filt_32(
    a: &[i16],
    m: usize,
    exc: &[i16],
    q_new: i16,
    hi: &mut [i16],
    lo: &mut [i16],
    lg: usize,
) {
    syn_filt_32_then(a, m, exc, q_new, hi, lo, lg, |_, _, _| {});
}

/// `syn_filt_32`, calling `post(n0, hi, lo)` on each block of B outputs of the fast path (for
/// the filters that follow it, so their chains run beside this one). Returns false if the
/// exact chains made the outputs instead (then `post` saw void values).
#[allow(clippy::needless_range_loop)] // i is the sample: it indexes exc and the outputs
#[allow(clippy::too_many_arguments)] // syn_filt_32's, and the block hook
#[inline(always)]
fn syn_filt_32_then(
    a: &[i16],
    m: usize,
    exc: &[i16],
    q_new: i16,
    hi: &mut [i16],
    lo: &mut [i16],
    lg: usize,
    post: impl FnMut(usize, &[i32; B], &[i32; B]),
) -> bool {
    let s = norm_s(a[0]) - 2;
    let a0 = shl_int16(a[0], -(4 + q_new));
    if m == M && lg.is_multiple_of(B) {
        let sum_a: i64 = a[1..=M].iter().map(|&c| (c as i64).abs()).sum();
        let lmax = max_abs(&lo[..M]).max(4095);
        let fixed = 2 * lmax * sum_a / 4096 + 1 + 2 * max_abs(&exc[..lg]) * (a0 as i64).abs();
        // (the history's H alone already too large: no fast pass; the check after it decides)
        if 2 * lmax * sum_a < 1 << 31 && fixed + 2 * max_abs(&hi[..M]) * sum_a < 1 << 31 {
            syn32_fast(a, a0, s, exc, hi, lo, lg, post);
            if fixed + 2 * max_abs(&hi[..M + lg]) * sum_a < 1 << 31 {
                return true;
            }
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
    false
}

/// max |v| (0 for none).
#[inline(always)]
fn max_abs(v: &[i16]) -> i64 {
    v.iter().map(|&x| (x as i32).abs()).max().unwrap_or(0) as i64
}

/// Outputs per block of the fast synthesis filters.
const B: usize = 8;

/// `c[j][i]`: 2·(the coefficient of y[n0 - N + j] in output n0 + i) for the N outputs before
/// a block, 0 where that sample is not among output i's N predecessors outside the block.
#[inline(always)]
fn hist_coefs<const N: usize>(a: &[i16]) -> [[i32; B]; N] {
    std::array::from_fn(|j| {
        std::array::from_fn(|i| if j >= i { 2 * a[N + i - j] as i32 } else { 0 })
    })
}

/// `acc[i] -= Σ_j c[j][i]·h[j]`: the terms of the N outputs before a block, one pass of B
/// lanes per history sample (wrapping: the callers' bounds keep every partial sum in range).
#[inline(always)]
fn hist_terms<const N: usize>(acc: &mut [i32; B], c: &[[i32; B]; N], h: &[i16; N]) {
    for (cj, &hj) in c.iter().zip(h) {
        for (s, &cji) in acc.iter_mut().zip(cj) {
            *s = s.wrapping_sub(cji.wrapping_mul(hj as i32));
        }
    }
}

/// `syn_filt_32`'s fast path (`lg` a multiple of B): blocks of B outputs, the terms of the M
/// outputs before the block for both halves as `hist_terms`, then per output only the terms
/// on earlier outputs of its block, the y[n-1] ones last (the serial chain). The same plain
/// sums as the C's chains reordered (exact under the caller's bound); the low-half sum is
/// complete before its `>> 12`. `post` gets each block's outputs.
#[allow(clippy::too_many_arguments)] // syn_filt_32's, and the block hook
#[inline(always)]
fn syn32_fast(
    a: &[i16],
    a0: i16,
    s: i16,
    exc: &[i16],
    hi: &mut [i16],
    lo: &mut [i16],
    lg: usize,
    mut post: impl FnMut(usize, &[i32; B], &[i32; B]),
) {
    let c = hist_coefs::<M>(a);
    let a2: [i32; B] = std::array::from_fn(|k| 2 * a[k] as i32);
    for n0 in (0..lg).step_by(B) {
        let mut acc_lo = [0i32; B];
        let mut acc_hi: [i32; B] = std::array::from_fn(|i| 2 * exc[n0 + i] as i32 * a0 as i32);
        hist_terms(&mut acc_lo, &c, lo[n0..n0 + M].try_into().unwrap());
        hist_terms(&mut acc_hi, &c, hi[n0..n0 + M].try_into().unwrap());
        let (mut yh, mut yl) = ([0i32; B], [0i32; B]);
        for m in 0..B {
            let (mut sl, mut sh) = (acc_lo[m], acc_hi[m]);
            for k in (1..=m).rev() {
                sl = sl.wrapping_sub(a2[k].wrapping_mul(yl[m - k]));
                sh = sh.wrapping_sub(a2[k].wrapping_mul(yh[m - k]));
            }
            let l = shl_int32((sl >> 12).wrapping_add(sh), 3 + s);
            yh[m] = l >> 16;
            yl[m] = (l >> 4) - (yh[m] << 12);
        }
        for m in 0..B {
            hi[M + n0 + m] = yh[m] as i16;
            lo[M + n0 + m] = yl[m] as i16;
        }
        post(n0, &yh, &yl);
    }
}

/// `wb_syn_filt` (in place, `update` = 1): 1/A(z) of order `m` = `a.len() - 1` on `xy`. As
/// TS 26.173's `Syn_filt`: scaled by `a[0]` (lowered by `Isp_Az`'s adaptive scaling for comfort
/// noise) and saturating; opencore assumed `a[0]` = 4096 and summed plainly.
///
/// Fast path (exact): with `|x| <= X` and `|y| <= Y` no partial sum reaches
/// `2·X·|a0| + 2·Y·Σ|a[1..]|`; below 2^31 the chain is a plain sum. Y is checked after the fast
/// pass, over the memory and the outputs it made (as in `syn_filt_32`).
#[inline(always)] // into simd::decode's AVX2 copy
pub(crate) fn wb_syn_filt(a: &[i16], xy: &mut [i16], mem: &mut [i16]) {
    let m = a.len() - 1;
    let lg = xy.len();
    let s = norm_s(a[0]) - 2;
    let a0 = a[0] >> 1;
    let mut yy = [0i16; M16K + L_SUBFR16K];
    yy[..m].copy_from_slice(mem);
    let sum_a: i64 = a[1..].iter().map(|&c| (c as i64).abs()).sum();
    let fixed = 2 * max_abs(xy) * (a0 as i64).abs();
    if lg.is_multiple_of(B) && lg <= L_SUBFR16K && fixed + 2 * max_abs(mem) * sum_a < 1 << 31 {
        let mut x = [0i16; L_SUBFR16K];
        x[..lg].copy_from_slice(xy);
        match m {
            M => wb_syn_fast::<M>(a, a0, s, xy, &mut yy),
            M16K => wb_syn_fast::<M16K>(a, a0, s, xy, &mut yy),
            _ => unreachable!("order {m}"),
        }
        // as in `syn_filt_32`: the bound over the outputs made
        if fixed + 2 * max_abs(xy) * sum_a < 1 << 31 {
            mem.copy_from_slice(&xy[lg - m..]);
            return;
        }
        xy.copy_from_slice(&x[..lg]);
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

/// `wb_syn_filt`'s fast path at order `N` (`yy` holds the memory, then the outputs; `xy` a
/// multiple of B long): blocks of B outputs as in `syn32_fast`.
pub(crate) fn wb_syn_fast<const N: usize>(
    a: &[i16],
    a0: i16,
    s: i16,
    xy: &mut [i16],
    yy: &mut [i16],
) {
    let c = hist_coefs::<N>(a);
    let a2: [i32; B] = std::array::from_fn(|k| 2 * a[k] as i32);
    for (n0, xb) in (0..).step_by(B).zip(xy.as_chunks_mut::<B>().0) {
        let mut acc: [i32; B] = std::array::from_fn(|i| 2 * xb[i] as i32 * a0 as i32);
        hist_terms(&mut acc, &c, yy[n0..n0 + N].try_into().unwrap());
        let mut y = [0i32; B];
        for m in 0..B {
            let mut l = acc[m];
            for k in (1..=m).rev() {
                l = l.wrapping_sub(a2[k].wrapping_mul(y[m - k]));
            }
            y[m] = amr_wb_round(shl_int32(l, 3 + s)) as i32;
        }
        for m in 0..B {
            xb[m] = y[m] as i16;
            yy[N + n0 + m] = y[m] as i16;
        }
    }
}

/// `deemphasis_32`: `y[i] = x[i] + mu·y[i-1]` on the (hi, lo) signal. Both additions
/// saturate as in TS 26.173 (`L_mac`); opencore adds plainly and wraps on loud signals.
pub(crate) fn deemphasis_32(x_hi: &[i16], x_lo: &[i16], y: &mut [i16], mu: i16, mem: &mut i16) {
    for i in 0..y.len() {
        y[i] = deemph_step(x_hi[i], x_lo[i], mu, mem);
    }
}

/// One `deemphasis_32` sample; `mem` is the previous output.
#[inline(always)]
fn deemph_step(x_hi: i16, x_lo: i16, mu: i16, mem: &mut i16) -> i16 {
    let mut l = mac_16by16_to_int32((x_hi as i32) << 16, x_lo, 8);
    l = shl_int32(l, 3);
    // round(L_shl(L_mac(l, mem, mu >> 1), 1)) as one clamp: saturation is a monotone clamp,
    // so clamp(2·clamp(v)) = clamp(2v), and the rounded high half of clamp(w) + 0x8000 is that
    // of w + 0x8000 clamped to 16 bits (mu >> 1 != -32768: L_mult is 2·a·b).
    let v = l as i64 + 2 * (*mem as i64 * (mu >> 1) as i64);
    *mem = ((2 * v + 0x8000) >> 16).clamp(-32768, 32767) as i16;
    *mem
}

/// `highpass_50Hz_at_12k8`: 2nd-order IIR, double-precision feedback, in place. The output
/// sum saturates as in TS 26.173 (an `L_mac` chain, then `L_shl`); opencore wraps there (plain
/// products, `<< 2`), which conformance vectors expose.
pub(crate) fn highpass_50hz_at_12k8(signal: &mut [i16], mem: &mut [i16; 6]) {
    let mut m = *mem;
    for s in signal.iter_mut() {
        *s = hp50_step(&mut m, *s);
    }
    *mem = m;
}

/// One `highpass_50Hz_at_12k8` sample; `mem` = [y2_hi, y2_lo, y1_hi, y1_lo, x0, x1].
#[inline(always)]
fn hp50_step(mem: &mut [i16; 6], x: i16) -> i16 {
    const A1: i16 = 16211;
    const A2: i16 = -8021;
    const B: [i16; 3] = [4053, -8106, 4053];
    let [y2_hi, y2_lo, y1_hi, y1_lo, x0, x1] = *mem;
    let (x2, x1, x0) = (x1, x0, x);
    // The C's chain: the low halves (|.| <= lo_max), then L_mac of y1_hi, y2_hi, x0, x1, x2.
    // Its first four partial sums are bounded by p3_max < 2^31 (y_lo in [0, 2^15), 16-bit rest):
    // they cannot saturate. Then L_shl(sat(sat(p3 + t4) + t5), 2) = clamp(4·(p3 + t4 + t5)):
    // clamp(4·clamp(u)) = clamp(4u), and if p3 + t4 saturates, |t5| < 3/4·2^31 keeps
    // 4·(±2^31 + t5) past the same end. So the step is one plain sum (any order: the y1 terms,
    // the feedback chain, last) and one clamp.
    const _: () = {
        let lo_max = ((1 << 15) * (A1 as i64 - A2 as i64) + 8192) >> 14;
        let p3_max = lo_max + 65536 * (A1 as i64 - A2 as i64 + B[0] as i64);
        assert!(p3_max < 1 << 31 && 65536 * (B[2] as i64).abs() < 3 << 29);
    };
    let off = 2 * (y2_hi as i64 * A2 as i64 + x0 as i64 * B[0] as i64)
        + 2 * (x1 as i64 * B[1] as i64 + x2 as i64 * B[2] as i64);
    let lo = fxp_mac_16by16(y2_lo, A2, fxp_mac_16by16(y1_lo, A1, 8192)) >> 14;
    let s = lo as i64 + off + 2 * (y1_hi as i64 * A1 as i64);
    let l = (4 * s).clamp(i32::MIN as i64, i32::MAX as i64) as i32;
    let (hi, lo) = int32_to_dpf(l);
    *mem = [y1_hi, y1_lo, hi, lo, x0, x1];
    amr_wb_shl1_round(l)
}

/// `highpass_400Hz_at_12k8`, in place.
#[cfg(test)]
pub(crate) fn highpass_400hz_at_12k8(signal: &mut [i16], mem: &mut [i16; 6]) {
    let mut m = *mem;
    for s in signal.iter_mut() {
        *s = hp400_step(&mut m, *s);
    }
    *mem = m;
}

/// One `highpass_400Hz_at_12k8` sample; `mem` as `hp50_step`'s.
#[inline(always)]
fn hp400_step(mem: &mut [i16; 6], x: i16) -> i16 {
    let [y2_hi, y2_lo, y1_hi, y1_lo, x0, x1] = *mem;
    let mut l1 = fxp_mac_16by16(y1_lo, 29280, 8192);
    let mut l2 = fxp_mul_16by16(y1_hi, 29280);
    l1 = fxp_mac_16by16(y2_lo, -14160, l1);
    l2 = fxp_mac_16by16(y2_hi, -14160, l2);
    let (x2, x1, x0) = (x1, x0, x);
    l2 = fxp_mac_16by16(x2, 915, l2);
    l2 = fxp_mac_16by16(x1, -1830, l2);
    l2 = fxp_mac_16by16(x0, 915, l2);
    l1 = (l1 >> 13).wrapping_add(l2.wrapping_shl(2));
    let hi = (l1 >> 16) as i16;
    let lo = (l1.wrapping_sub((hi as i32) << 16) >> 1) as i16;
    *mem = [y1_hi, y1_lo, hi, lo, x0, x1];
    (l1.wrapping_add(0x8000) >> 16) as i16
}

/// The 12.8 kHz synthesis: `syn_filt_32` on `exc` (its memories in `st`), then `deemphasis_32`, `highpass_50Hz_at_12k8` (into `synth`) and
/// `highpass_400Hz_at_12k8` (into `synth400`), one sample at a time: on `syn_filt_32`'s fast
/// path per block of its outputs, so the four serial chains run side by side.
#[inline(always)]
fn lf_synthesis(
    aq: &[i16],
    exc: &[i16],
    q_new: i16,
    synth: &mut [i16; L_SUBFR],
    synth400: &mut [i16; L_SUBFR],
    st: &mut SynthState,
) {
    let mut hi = [0i16; M + L_SUBFR];
    let mut lo = [0i16; M + L_SUBFR];
    hi[..M].copy_from_slice(&st.mem_syn_hi);
    lo[..M].copy_from_slice(&st.mem_syn_lo);
    // (deemphasis, 50 Hz, 400 Hz) memories
    type Mems = (i16, [i16; 6], [i16; 6]);
    let step = |(md, m50, m400): &mut Mems, x_hi: i16, x_lo: i16| {
        let s50 = hp50_step(m50, deemph_step(x_hi, x_lo, PREEMPH_FAC, md));
        (s50, hp400_step(m400, s50))
    };
    let mems: Mems = (st.mem_deemph, st.mem_sig_out, st.mem_hp400);
    let mut mm = mems;
    let fused = syn_filt_32_then(
        aq,
        M,
        exc,
        q_new,
        &mut hi,
        &mut lo,
        L_SUBFR,
        |n0, yh, yl| {
            for m in 0..B {
                (synth[n0 + m], synth400[n0 + m]) = step(&mut mm, yh[m] as i16, yl[m] as i16);
            }
        },
    );
    if !fused {
        mm = mems;
        for i in 0..L_SUBFR {
            (synth[i], synth400[i]) = step(&mut mm, hi[M + i], lo[M + i]);
        }
    }
    (st.mem_deemph, st.mem_sig_out, st.mem_hp400) = mm;
    st.mem_syn_hi.copy_from_slice(&hi[L_SUBFR..]);
    st.mem_syn_lo.copy_from_slice(&lo[L_SUBFR..]);
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
    // Output j = 5m + p (m < 16) reads signal[1 + 4m + d + t] (t < 24, d = p - 1) for p > 0
    // and copies signal[L_FILT + 4m] for p = 0. Polyphase: with the input split by index mod 4,
    // each tap of a phase is one pass over 16 contiguous lanes (m). Plain sums, as the C.
    const N: usize = L_SUBFR16K / 5;
    let ph: [[i16; 22]; 4] = std::array::from_fn(|r| std::array::from_fn(|q| signal[4 * q + r]));
    for p in 1..5 {
        let fir = &FIR_UP[(p - 1) * 24..p * 24];
        let mut acc = [0x2000i32; N];
        for (t, &c) in fir.iter().enumerate() {
            let r = p + t; // 1 + d + t
            let w: &[i16; N] = ph[r % 4][r / 4..r / 4 + N].try_into().unwrap();
            for (s, &v) in acc.iter_mut().zip(w) {
                *s = fxp_mac_16by16(v, c, *s);
            }
        }
        for (m, &s) in acc.iter().enumerate() {
            sig16k[5 * m + p] = (shl_int32(s, 2) >> 16) as i16;
        }
    }
    for m in 0..N {
        sig16k[5 * m] = signal[L_FILT + 4 * m];
    }
    mem.copy_from_slice(&signal[L_SUBFR..L_SUBFR + 2 * L_FILT]);
}

/// A 31-tap FIR over `x` (30 past samples, then 80 new) into `out`: TS 26.173's `L_mac`
/// chain and `round`. Fast path (exact): while `2·max|x|·Σ|c| < 2^31` no partial sum
/// saturates and the chain is the plain dot product. opencore summed plainly (wrapping on loud
/// input) and, at 7 kHz, pre-added the end taps' samples in 16 bits.
/// `x` ends in one zero of slack: the fast path runs 32 taps (the last coefficient 0), one
/// whole vector width.
#[inline(always)] // into simd::decode's AVX2 copy
fn fir31(x: &[i16; L_FIR + L_SUBFR16K + 1], fir: &[i16; L_FIR + 1], out: &mut [i16; L_SUBFR16K]) {
    let abs_sum: i64 = fir.iter().map(|&c| (c as i64).abs()).sum();
    let max = x.iter().map(|&v| (v as i64).abs()).max().unwrap();
    let fir32: [i16; L_FIR + 2] = std::array::from_fn(|k| if k <= L_FIR { fir[k] } else { 0 });
    if 2 * max * abs_sum < 1 << 31 {
        // transposed: the outputs in lanes, one pass per tap
        let mut acc = [0i32; L_SUBFR16K];
        for (k, &c) in fir32.iter().enumerate() {
            let w: &[i16; L_SUBFR16K] = x[k..k + L_SUBFR16K].try_into().unwrap();
            for (s, &v) in acc.iter_mut().zip(w) {
                *s = s.wrapping_add(v as i32 * c as i32);
            }
        }
        for (o, &l) in out.iter_mut().zip(&acc) {
            *o = amr_wb_round(l.wrapping_shl(1));
        }
        return;
    }
    for (n, o) in out.iter_mut().enumerate() {
        let w = &x[n..n + L_FIR + 1];
        let l = w
            .iter()
            .zip(fir)
            .fold(0, |l, (&a, &c)| mac_16by16_to_int32(l, a, c));
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

/// `noise_gen_amrwb(seed) >> 3` for each of `hf` in turn, as 8 interleaved generators: the LCG
/// s' = 31821·s + 13849 (mod 2^16) jumps 8 steps as s' = A8·s + C8, so after the first 8
/// samples each lane steps on its own (no serial multiply chain).
fn noise_block(hf: &mut [i16; L_SUBFR16K], seed: &mut i16) {
    const J: (i16, i16) = {
        let (mut a, mut c, mut k) = (1i16, 0i16, 0);
        while k < 8 {
            (a, c) = (
                a.wrapping_mul(31821),
                c.wrapping_mul(31821).wrapping_add(13849),
            );
            k += 1;
        }
        (a, c)
    };
    let mut lane = [0i16; 8];
    for v in lane.iter_mut() {
        *v = noise_gen_amrwb(seed);
    }
    for (b, blk) in hf.as_chunks_mut::<8>().0.iter_mut().enumerate() {
        if b > 0 {
            for s in lane.iter_mut() {
                *s = s.wrapping_mul(J.0).wrapping_add(J.1);
            }
        }
        for (o, &s) in blk.iter_mut().zip(&lane) {
            *o = s >> 3;
        }
    }
    *seed = lane[7];
}

#[cfg(test)]
#[test]
fn noise_block_is_the_generator() {
    for seed0 in [0i16, 1, -1, 21845, i16::MIN, i16::MAX, 12345] {
        let (mut s, mut t) = (seed0, seed0);
        let mut hf = [0i16; L_SUBFR16K];
        noise_block(&mut hf, &mut s);
        let want: [i16; L_SUBFR16K] = std::array::from_fn(|_| noise_gen_amrwb(&mut t) >> 3);
        assert_eq!((hf, s), (want, t), "seed {seed0}");
    }
}

/// The high-band gain's energy and neighbour correlation of `synth`: the C's chains
/// `L_mac(1, s[0], s[0])`, then per i >= 1 `L_mac(l, s[i], s[i])` and `L_mac(l2, s[i], s[i-1])`
/// (from 1). The energy's terms are >= 0: it is min(1 + Σ 2·s², MAX_32). Below MAX_32 no
/// partial sum of the correlation saturates either (|2·s[i]·s[i-1]| <= s[i]² + s[i-1]²), so it
/// is a plain sum; past it, the chain.
fn tilt_energies(synth: &[i16; L_SUBFR]) -> (i32, i32) {
    let e: u64 = synth.iter().map(|&v| (v as i32 * v as i32) as u64).sum();
    if 2 * e < MAX_32 as u64 {
        let c = synth
            .windows(2)
            .map(|w| w[0] as i32 * w[1] as i32)
            .sum::<i32>();
        ((1 + 2 * e) as i32, 1 + 2 * c)
    } else {
        let l2 = (1..L_SUBFR).fold(1, |l2, i| mac_16by16_to_int32(l2, synth[i], synth[i - 1]));
        (MAX_32, l2)
    }
}

#[cfg(test)]
#[test]
fn tilt_energies_match_the_chains() {
    let mut r = crate::testutil::Rng::new(7);
    for k in 0..20_000 {
        // quiet to full scale: the shift; k % 3 == 0 pins samples to the extremes
        let shift = (r.u64() % 12) as u32;
        let s: [i16; L_SUBFR] = std::array::from_fn(|_| match (k % 3, r.u64() % 4) {
            (0, 0) => i16::MIN,
            (0, 1) => i16::MAX,
            _ => r.i16() >> shift,
        });
        let mut l = mac_16by16_to_int32(1, s[0], s[0]);
        let mut l2 = 1;
        for i in 1..L_SUBFR {
            l = mac_16by16_to_int32(l, s[i], s[i]);
            l2 = mac_16by16_to_int32(l2, s[i], s[i - 1]);
        }
        assert_eq!(tilt_energies(&s), (l, l2), "case {k}");
    }
}

/// `synthesis_amr_wb`: one subframe from the excitation `exc` (12.8 kHz, modified in place as
/// the C does) to 80 samples at 16 kHz in `synth16k`.
#[allow(clippy::too_many_arguments)] // the C's signature
#[inline(always)] // into simd::decode's AVX2 copy
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
    let mut synth = [0i16; L_SUBFR];
    let mut hf = [0i16; L_SUBFR16K];
    let mut synth400 = [0i16; L_SUBFR];
    lf_synthesis(aq, exc, q_new, &mut synth, &mut synth400, st);
    crate::stages::mark(crate::stages::LF);
    oversamp_12k8_to_16k(&synth, synth16k, &mut st.mem_oversamp);
    crate::stages::mark(crate::stages::OVERSAMP);

    // high band: white noise scaled to the excitation energy
    noise_block(&mut hf, &mut st.seed2);
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
    let (l, l2) = tilt_energies(&synth400);
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

    crate::stages::mark(crate::stages::HF_GAIN);
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
    crate::stages::mark(crate::stages::HF_SYN);
    band_pass_6k_7k(&mut hf, &mut st.mem_hf);
    if nb_bits >= NBBITS_24K {
        low_pass_filt_7k(&mut hf, &mut st.mem_hf3);
    }
    for (s, &h) in synth16k[..L_SUBFR16K].iter_mut().zip(&hf) {
        *s = add_int16(*s, h);
    }
    crate::stages::mark(crate::stages::HF_FIR);
}
