//! Gain decoding (`dec_gain2_amr_wb.cpp`): pitch gain and MA-predicted code gain, with the
//! frame-erasure concealment.

use crate::basic_op::*;
use crate::excitation::median5;
use crate::tables::*;

const MEAN_ENER: i16 = 30;

/// Gain decoder memory (`dec_gain[23]` in the C): past quantised energies `[0..4]`, past pitch
/// gain `[4]`, past code gain `[5]`, previous code gain `[6]`, pitch-gain buffer `[7..12]`,
/// code-gain buffer `[12..17]`, the pitch gains `lagconceal` reads `[17..22]`.
pub(crate) type GainMem = [i16; 23];
pub(crate) const PBUF2: usize = 17;

/// `dec_gain2_amr_wb_init`.
pub(crate) fn dec_gain2_amr_wb_init(mem: &mut GainMem) {
    mem[..4].fill(-14336);
    mem[4..22].fill(0);
    mem[22] = 21845;
}

/// `dec_gain2_amr_wb`: (gain_pit Q14, gain_cod Q16) for one subframe.
#[allow(clippy::too_many_arguments)] // the C's signature
pub(crate) fn dec_gain2_amr_wb(
    index: i16,
    nbits: i16,
    code: &[i16],
    bfi: bool,
    prev_bfi: bool,
    state: i16,
    unusable_frame: bool,
    vad_hist: i16,
    mem: &mut GainMem,
) -> (i16, i32) {
    const PAST_GAIN_PIT: usize = 4;
    const PAST_GAIN_CODE: usize = 5;
    const PREV_GC: usize = 6;
    const PBUF: usize = 7;
    const GBUF: usize = 12;
    let mut exp = 0i16;
    let mut l = dot_product12(code, code, &mut exp);
    exp = exp.wrapping_sub(24);
    one_ov_sqrt_norm(&mut l, &mut exp);
    let gcode_inov = extract_h(shl_int32(l, exp.wrapping_sub(3)));
    let st = state as usize;

    if bfi {
        mem[PAST_GAIN_PIT] = median5(mem, PBUF + 2).min(15565);
        let gain_pit = if unusable_frame {
            mult_int16(PDOWN_UNUSABLE[st], mem[PAST_GAIN_PIT])
        } else {
            mult_int16(PDOWN_USABLE[st], mem[PAST_GAIN_PIT])
        };
        let tmp = median5(mem, GBUF + 2);
        mem[PAST_GAIN_CODE] = if vad_hist > 2 {
            tmp
        } else if unusable_frame {
            mult_int16(CDOWN_UNUSABLE[st], tmp)
        } else {
            mult_int16(CDOWN_USABLE[st], tmp)
        };
        // the mean of the past energies, -3 dB, as TS 26.173 (opencore: sum / 8, and its
        // pipelined copy keeps q[3] instead of shifting it)
        let q = &mut mem[..4];
        let sum = q[0] as i32 + q[1] as i32 + q[2] as i32 + q[3] as i32;
        let qua_ener = ((sum >> 2) as i16).wrapping_sub(3072).max(-14336);
        q.copy_within(0..3, 1);
        q[0] = qua_ener;
        mem.copy_within(GBUF + 1..GBUF + 5, GBUF);
        mem.copy_within(PBUF + 1..PBUF + 5, PBUF);
        mem[GBUF + 4] = mem[PAST_GAIN_CODE];
        mem[PBUF + 4] = mem[PAST_GAIN_PIT];
        return (
            gain_pit,
            mul_16by16_to_int32(mem[PAST_GAIN_CODE], gcode_inov),
        );
    }

    let mut l = shl_int32(l_deposit_h(MEAN_ENER), 8);
    for k in 0..4 {
        l = mac_16by16_to_int32(l, PRED[k], mem[k]);
    }
    let gcode0 = extract_h(l);
    let l = (gcode0 as i32 * 5443) >> 7;
    let (exp_gcode0, frac) = int32_to_dpf(l);
    let gcode0 = power_of_2(14, frac) as i16;
    let exp_gcode0 = exp_gcode0.wrapping_sub(14);
    let p = if nbits == 6 {
        &T_QUA_GAIN6B[..]
    } else {
        &T_QUA_GAIN7B[..]
    };
    let gain_pit = p[(index as usize) << 1];
    let g_code = p[((index as usize) << 1) + 1];
    let mut gain_cod = shl_int32(
        mul_16by16_to_int32(g_code, gcode0),
        exp_gcode0.wrapping_add(4),
    );
    if prev_bfi {
        let l = mul_16by16_to_int32(mem[PREV_GC], 5120);
        if gain_cod > l && gain_cod > 6_553_600 {
            gain_cod = l;
        }
    }
    mem[PAST_GAIN_CODE] = amr_wb_round(shl_int32(gain_cod, 3));
    mem[PAST_GAIN_PIT] = gain_pit;
    mem[PREV_GC] = mem[PAST_GAIN_CODE];
    // shift the histories (opencore's pipelined copy leaves [b1, b1, b2, b3], not the
    // [b1, b2, b3, b4] of TS 26.173: wrong concealment gains after a bad frame)
    for buf in [GBUF, PBUF, PBUF2] {
        mem.copy_within(buf + 1..buf + 5, buf);
    }
    mem[GBUF + 4] = mem[PAST_GAIN_CODE];
    mem[PBUF + 4] = mem[PAST_GAIN_PIT];
    mem[PBUF2 + 4] = mem[PAST_GAIN_PIT];
    let (hi, lo) = int32_to_dpf(gain_cod);
    let gain_cod = shl_int32(mul_32by16(hi, lo, gcode_inov), 3);
    mem.copy_within(0..3, 1);
    let (mut e, mut f) = (0i16, 0i16);
    amrwb_log_2(g_code as i32, &mut e, &mut f);
    let e = e.wrapping_sub(11);
    mem[0] = (mul_32by16(e, f, 24660) >> 3) as i16;
    (gain_pit, gain_cod)
}
