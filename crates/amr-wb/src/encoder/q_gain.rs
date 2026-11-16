//! Joint pitch / code gain quantisation (vo-amrwbenc `q_gain2.c`), with the reference's
//! saturation: MA-predicted code gain, a table search minimising the weighted error.

use super::ops::*;
use crate::basic_op::dot_product12;
use crate::enc_tables::{Q_GAIN2_PRED as PRED, T_QUA_GAIN6B, T_QUA_GAIN7B};

const MEAN_ENER: i16 = 30;
const RANGE: usize = 64;
const NB_QUA_GAIN7B: usize = 128;

/// `Init_Q_gain2`: the past quantised energies.
pub(crate) const Q_GAIN2_INIT: [i16; 4] = [-14336; 4];

/// `Q_gain2`: the gain codebook index (`nbits` 6 or 7); `gain_pit` (Q14, in: the pitch gain
/// found, out: quantised), `gain_cod` the code gain (Q16). `g_coeff` from `G_pitch`; `q_xn`
/// the scaling of `xn`, `y1`.
pub(crate) fn q_gain2(
    xn: &[i16],
    y1: &[i16],
    q_xn: i16,
    y2: &[i16],
    code: &[i16],
    g_coeff: &[i16; 4],
    nbits: i16,
    gain_pit: &mut i16,
    gain_cod: &mut i32,
    gp_clip: bool,
    past_qua_en: &mut [i16; 4],
) -> i16 {
    let (t_qua_gain, min_ind, size): (&[i16], usize, usize) = if nbits == 6 {
        (&T_QUA_GAIN6B, 0, if gp_clip { RANGE - 16 } else { RANGE })
    } else {
        let j = NB_QUA_GAIN7B - RANGE - if gp_clip { 27 } else { 0 };
        let mut min_ind = 0;
        for i in 0..j {
            // the C: p = t_qua_gain7b + RANGE (a quarter in), p += 2
            if *gain_pit > T_QUA_GAIN7B[RANGE + 2 * i] {
                min_ind += 1;
            }
        }
        (&T_QUA_GAIN7B, min_ind, RANGE)
    };
    let mut coeff = [0i16; 5];
    let mut coeff_lo = [0i16; 5];
    let mut exp_coeff = [0i16; 5];
    coeff[0] = g_coeff[0];
    exp_coeff[0] = g_coeff[1];
    coeff[1] = negate(g_coeff[2]);
    exp_coeff[1] = add(g_coeff[3], 1);
    let mut exp = 0;
    coeff[2] = extract_h(dot_product12(y2, y2, &mut exp));
    exp_coeff[2] = add(sub(exp, 18), shl(q_xn, 1));
    coeff[3] = extract_h(l_negate(dot_product12(xn, y2, &mut exp)));
    exp_coeff[3] = add(sub(exp, 9 - 1), q_xn);
    coeff[4] = extract_h(dot_product12(y1, y2, &mut exp));
    exp_coeff[4] = add(sub(exp, 9 - 1), q_xn);
    let mut exp_code = 0;
    let l = dot_product12(code, code, &mut exp_code);
    let exp_code = sub(exp_code, 18 + 6 + 31);
    let (e, frac) = log2(l);
    let e = add(e, exp_code);
    let mut l = mpy_32_16(e, frac, -24660);
    l = l_mac(l, MEAN_ENER, 8192);
    l = l_shl(l, 10);
    for k in 0..4 {
        l = l_mac(l, PRED[k], past_qua_en[k]);
    }
    let gcode0 = extract_h(l);
    let l = l_shr(l_mult(gcode0, 5443), 8);
    let (exp_gcode0, frac) = l_extract(l);
    let gcode0 = extract_l(pow2(14, frac));
    let exp_gcode0 = sub(exp_gcode0, 14);
    let exp_code = add(exp_gcode0, 4);
    let exp_max = [
        sub(exp_coeff[0], 13),
        sub(exp_coeff[1], 14),
        add(exp_coeff[2], add(15, shl(exp_code, 1))),
        add(exp_coeff[3], exp_code),
        add(exp_coeff[4], add(1, exp_code)),
    ];
    let mut e_max = exp_max[0];
    for &e in &exp_max[1..] {
        if e > e_max {
            e_max = e;
        }
    }
    for i in 0..5 {
        let j = add(sub(e_max, exp_max[i]), 2);
        let l = l_shr(l_deposit_h(coeff[i]), j);
        (coeff[i], coeff_lo[i]) = l_extract(l);
        coeff_lo[i] = shr(coeff_lo[i], 3);
    }
    let mut dist_min = i32::MAX;
    let mut index = 0;
    for i in 0..size {
        let g_pitch = t_qua_gain[2 * (min_ind + i)];
        let g_code = mult_r(t_qua_gain[2 * (min_ind + i) + 1], gcode0);
        let g2_pitch = mult_r(g_pitch, g_pitch);
        let g_pit_cod = mult_r(g_code, g_pitch);
        let (g2_code, g2_code_lo) = l_extract(l_mult(g_code, g_code));
        let mut l = l_shr(l_mult(coeff[2], g2_code_lo), 3);
        l = l_mac(l, coeff_lo[0], g2_pitch);
        l = l_mac(l, coeff_lo[1], g_pitch);
        l = l_mac(l, coeff_lo[2], g2_code);
        l = l_mac(l, coeff_lo[3], g_code);
        l = l_mac(l, coeff_lo[4], g_pit_cod);
        l = l_shr(l, 12);
        l = l_mac(l, coeff[0], g2_pitch);
        l = l_mac(l, coeff[1], g_pitch);
        l = l_mac(l, coeff[2], g2_code);
        l = l_mac(l, coeff[3], g_code);
        l = l_mac(l, coeff[4], g_pit_cod);
        if l < dist_min {
            dist_min = l;
            index = i;
        }
    }
    let index = index + min_ind;
    *gain_pit = t_qua_gain[2 * index];
    let g_code = t_qua_gain[2 * index + 1];
    *gain_cod = l_shl(l_mult(g_code, gcode0), add(exp_gcode0, 4));
    let (e, frac) = log2(g_code as i32);
    let l = mpy_32_16(sub(e, 11), frac, 24660);
    let qua_ener = extract_l(l_shr(l, 3));
    past_qua_en.copy_within(0..3, 1);
    past_qua_en[0] = qua_ener;
    index as i16
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::cnst::L_SUBFR;
    use crate::encoder::filters::vs_ref::signal;
    use crate::encoder::pitch::g_pitch;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_Q_gain2(
            xn: *mut i16,
            y1: *mut i16,
            q_xn: i16,
            y2: *mut i16,
            code: *mut i16,
            g_coeff: *mut i16,
            l_subfr: i16,
            nbits: i16,
            gain_pit: *mut i16,
            gain_cod: *mut i32,
            gp_clip: i16,
            mem: *mut i16,
        ) -> i16;
    }

    /// Sequences (the predictor memory) over both table sizes, clipping on and off, targets
    /// and filtered vectors from the signal classes, sparse codevectors (Q9 pulses).
    #[test]
    fn q_gain2_matches_reference() {
        let mut r = Rng::new(101);
        for seq in 0..100u64 {
            let (mut mem, mut memc) = (Q_GAIN2_INIT, Q_GAIN2_INIT);
            for call in 0..200u64 {
                let mut xn = signal(&mut r, seq + call, L_SUBFR);
                let mut y1 = signal(&mut r, seq + call / 2, L_SUBFR);
                let mut y2 = signal(&mut r, call, L_SUBFR);
                let mut code = vec![0i16; L_SUBFR];
                for _ in 0..1 + r.u64() % 24 {
                    code[(r.u64() % 64) as usize] +=
                        if r.u64().is_multiple_of(2) { 512 } else { -512 };
                }
                let mut g = [0i16; 4];
                let gp = g_pitch(&xn, &y1, &mut g);
                let q_xn = (r.u64() % 12) as i16;
                let nbits = 6 + (r.u64() % 2) as i16;
                let clip = r.u64().is_multiple_of(3);
                let (mut gain_pit, mut gain_cod) = (gp, 0i32);
                let (mut gpc, mut gcc) = (gp, 0i32);
                let ix = q_gain2(
                    &xn,
                    &y1,
                    q_xn,
                    &y2,
                    &code,
                    &g,
                    nbits,
                    &mut gain_pit,
                    &mut gain_cod,
                    clip,
                    &mut mem,
                );
                let ixc = with_ref(|| {
                    // SAFETY: 64-sample vectors, the 4 correlations and memory words.
                    unsafe {
                        refwb_Q_gain2(
                            xn.as_mut_ptr(),
                            y1.as_mut_ptr(),
                            q_xn,
                            y2.as_mut_ptr(),
                            code.as_mut_ptr(),
                            g.as_mut_ptr(),
                            64,
                            nbits,
                            &mut gpc,
                            &mut gcc,
                            clip as i16,
                            memc.as_mut_ptr(),
                        )
                    }
                });
                assert_eq!(
                    (ix, gain_pit, gain_cod, mem),
                    (ixc, gpc, gcc, memc),
                    "Q_gain2 {seq} {call} nbits {nbits} clip {clip} q_xn {q_xn} gp {gp} g {g:?}"
                );
            }
        }
    }
}
