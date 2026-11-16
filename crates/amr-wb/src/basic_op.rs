//! Basic operators and math, named as opencore's AMR-WB decoder (`pvamrwbdecoder_basic_op*.h`,
//! `pvamrwb_math_op.cpp`, `normalize_amr_wb.cpp`) with TS 26.173's semantics: no overflow
//! flag (the decoder reads none). Where opencore's versions differ from the ETSI operators
//! they stand for (`round`, the shifts, `norm_s`/`norm_l`), these follow the specification;
//! each such spot says what opencore did. C's plain `int` arithmetic is `wrapping_*`.

use crate::tables::{LOG2_NORM_TABLE, TABLE_ISQRT, TABLE_POW2};

pub(crate) const MAX_16: i16 = i16::MAX;
pub(crate) const MAX_32: i32 = i32::MAX;

#[inline(always)]
pub(crate) fn negate_int16(a: i16) -> i16 {
    if a == i16::MIN { MAX_16 } else { -a }
}

// The shifts are TS 26.173's `shl`, `L_shl`, `L_shr` for every count: opencore masks the
// counts (4 bits for 16-bit shifts and negative 32-bit ones, 5 bits for `shr_int32`), so a
// large right shift of a small energy kept bits the spec shifts out (bad-frame high band).

/// Saturating 16-bit shift left; negative counts shift right (by at most 15).
#[inline(always)]
pub(crate) fn shl_int16(a: i16, b: i16) -> i16 {
    if b < 0 {
        a >> (-(b as i32)).min(15)
    } else {
        let r = (a as i32) << b.min(16);
        if r != r as i16 as i32 {
            (a >> 15) ^ MAX_16
        } else {
            r as i16
        }
    }
}

/// Saturating 32-bit shift left; negative counts shift right (by at most 31).
#[inline(always)]
pub(crate) fn shl_int32(a: i32, b: i16) -> i32 {
    if b > 0 {
        let r = (a as i64) << b.min(32);
        if r != r as i32 as i64 {
            (a >> 31) ^ MAX_32
        } else {
            r as i32
        }
    } else {
        a >> (-(b as i32)).min(31)
    }
}

/// Arithmetic 32-bit shift right (by at most 31); negative counts shift left, saturating.
#[inline(always)]
pub(crate) fn shr_int32(a: i32, b: i16) -> i32 {
    if b >= 0 {
        a >> b.min(31)
    } else {
        shl_int32(a, (-(b as i32)).min(32) as i16)
    }
}

#[inline(always)]
pub(crate) fn add_int16(a: i16, b: i16) -> i16 {
    a.saturating_add(b)
}

#[inline(always)]
pub(crate) fn sub_int16(a: i16, b: i16) -> i16 {
    a.saturating_sub(b)
}

#[inline(always)]
pub(crate) fn mult_int16(a: i16, b: i16) -> i16 {
    let p = (a as i32 * b as i32) >> 15;
    if (p >> 15) != (p >> 31) {
        ((p >> 31) ^ MAX_16 as i32) as i16
    } else {
        p as i16
    }
}

#[inline(always)]
pub(crate) fn mult_int16_r(a: i16, b: i16) -> i16 {
    let p = (a as i32 * b as i32 + 0x4000) >> 15;
    if (p >> 15) != (p >> 31) {
        ((p >> 31) ^ MAX_16 as i32) as i16
    } else {
        p as i16
    }
}

#[inline(always)]
pub(crate) fn add_int32(a: i32, b: i32) -> i32 {
    a.saturating_add(b)
}

#[inline(always)]
pub(crate) fn sub_int32(a: i32, b: i32) -> i32 {
    a.saturating_sub(b)
}

/// `L_mult`: `2·a·b`, saturating only for `a = b = -32768`.
#[inline(always)]
pub(crate) fn mul_16by16_to_int32(a: i16, b: i16) -> i32 {
    let p = a as i32 * b as i32;
    if p != 0x4000_0000 { p << 1 } else { MAX_32 }
}

#[inline(always)]
pub(crate) fn mac_16by16_to_int32(acc: i32, a: i16, b: i16) -> i32 {
    acc.saturating_add(mul_16by16_to_int32(a, b))
}

#[inline(always)]
pub(crate) fn msu_16by16_from_int32(acc: i32, a: i16, b: i16) -> i32 {
    acc.saturating_sub(mul_16by16_to_int32(a, b))
}

/// Rounds to the high half, saturating: TS 26.173's `round`. opencore only special-cases
/// `MAX_32` and wraps `[0x7fff8000, 0x7ffffffe]` to -32768.
#[inline(always)]
pub(crate) fn amr_wb_round(a: i32) -> i16 {
    (a.saturating_add(0x8000) >> 16) as i16
}

/// `round(L_shl(a, 1))`, saturating both (opencore's version wraps near `MAX_32 / 2`).
#[inline(always)]
pub(crate) fn amr_wb_shl1_round(a: i32) -> i16 {
    amr_wb_round(shl_int32(a, 1))
}

/// `(hi, lo)` double precision times 16 bits (wrapping, as the C).
#[inline(always)]
pub(crate) fn mul_32by16(hi: i16, lo: i16, n: i16) -> i32 {
    ((hi as i32 * n as i32).wrapping_add((lo as i32 * n as i32) >> 15)).wrapping_shl(1)
}

/// Plain `acc + a·b` (wrapping).
#[inline(always)]
pub(crate) fn fxp_mac_16by16(a: i16, b: i16, acc: i32) -> i32 {
    acc.wrapping_add(a as i32 * b as i32)
}

#[inline(always)]
pub(crate) fn fxp_mul_16by16(a: i16, b: i16) -> i32 {
    a as i32 * b as i32
}

#[inline(always)]
pub(crate) fn extract_h(a: i32) -> i16 {
    (a >> 16) as i16
}

#[inline(always)]
pub(crate) fn l_deposit_h(a: i16) -> i32 {
    (a as i32) << 16
}

/// `norm_l`: the left shift that normalises `x` (0 for 0; negative `x` as `!x`), TS 26.173's.
/// opencore's `normalize_amr_wb` agrees for positive `x` only (27 for 0), and its `norm_s` is
/// that minus 16: 11 for a zero code gain, which skewed the voicing after lost frames.
#[inline(always)]
pub(crate) fn normalize_amr_wb(x: i32) -> i16 {
    if x == 0 {
        0
    } else {
        ((if x < 0 { !x } else { x }).leading_zeros() - 1) as i16
    }
}

/// `norm_s`: as [`normalize_amr_wb`], for 16 bits.
#[inline(always)]
pub(crate) fn norm_s(x: i16) -> i16 {
    if x == 0 {
        0
    } else {
        ((if x < 0 { !x } else { x }).leading_zeros() - 1) as i16
    }
}

/// Arithmetic shift right with rounding (count masked to 4 bits).
pub(crate) fn shr_rnd(a: i16, b: i16) -> i16 {
    let mut r = a >> (b & 0xf);
    if b != 0 && (a as i32 & (1i32 << (b - 1))) != 0 {
        r = r.wrapping_add(1);
    }
    r
}

/// `a / b` in Q15 for `0 <= a <= b`; 0 outside that range.
pub(crate) fn div_16by16(a: i16, b: i16) -> i16 {
    if a > b || a < 0 {
        return 0;
    }
    if a == 0 {
        return 0;
    }
    if a == b {
        return MAX_16;
    }
    let mut out: i16 = 0;
    let mut num = a as i32;
    let den = b as i32;
    let (den2, den4) = (den << 1, den << 2);
    for _ in 0..5 {
        out <<= 3;
        num <<= 3;
        if num >= den4 {
            num -= den4;
            out |= 4;
        }
        if num >= den2 {
            num -= den2;
            out |= 2;
        }
        if num >= den {
            num -= den;
            out |= 1;
        }
    }
    out
}

/// `1/sqrt(frac · 2^exp)` as (frac, exp), in place.
pub(crate) fn one_ov_sqrt_norm(frac: &mut i32, exp: &mut i16) {
    if *frac <= 0 {
        *exp = 0;
        *frac = 0x7fff_ffff;
        return;
    }
    if *exp & 1 == 1 {
        *frac >>= 1;
    }
    *exp = negate_int16(((*exp as i32 - 1) >> 1) as i16);
    *frac >>= 9;
    let i = extract_h(*frac) - 16;
    *frac >>= 1;
    let a = (*frac as i16) & 0x7fff;
    let i = i as usize;
    *frac = l_deposit_h(TABLE_ISQRT[i]);
    let tmp = TABLE_ISQRT[i] - TABLE_ISQRT[i + 1];
    *frac = msu_16by16_from_int32(*frac, tmp, a);
}

pub(crate) fn one_ov_sqrt(x: i32) -> i32 {
    let e = normalize_amr_wb(x);
    let mut x = x.wrapping_shl(e as u32);
    let mut e = 31 - e;
    one_ov_sqrt_norm(&mut x, &mut e);
    shl_int32(x, e)
}

/// `L_shr_r`: arithmetic shift right with rounding; counts past 31 give 0, negative counts
/// shift left saturating.
#[inline(always)]
pub(crate) fn l_shr_r(x: i32, n: i16) -> i32 {
    if n > 31 {
        return 0;
    }
    let r = shr_int32(x, n);
    if n > 0 && (x >> (n - 1)) & 1 != 0 {
        r + 1
    } else {
        r
    }
}

/// `2^(exponant + fraction/2^15)` (`fraction` in [0, 32767]): TS 26.173's `Pow2`. opencore
/// computed it at half scale and shifted with plain C shifts, x86-masked for exponents past
/// 29 (loud comfort noise came out wrong).
pub(crate) fn power_of_2(exponant: i16, fraction: i16) -> i32 {
    let l = (fraction as i32) << 6;
    let i = (l >> 16) as usize;
    let a = ((l >> 1) as i16) & 0x7fff;
    let tmp = TABLE_POW2[i] - TABLE_POW2[i + 1];
    let l = msu_16by16_from_int32(l_deposit_h(TABLE_POW2[i]), tmp, a);
    l_shr_r(l, sub_int16(30, exponant))
}

/// Normalised `Σ x·y` (+1, saturating), as (mantissa, exponent). `lg` a multiple of 8.
pub(crate) fn dot_product12(x: &[i16], y: &[i16], exp: &mut i16) -> i32 {
    debug_assert!(x.len() == y.len() && x.len().is_multiple_of(8));
    let sum = if std::ptr::eq(x, y) {
        // an energy: the terms are >= 0 (L_mult saturates only up), so the chain is
        // min(1 + Σ 2·a², MAX_32), a plain sum
        let e: u64 = x.iter().map(|&a| (a as i32 * a as i32) as u64).sum();
        (1 + 2 * e).min(MAX_32 as u64) as i32
    } else {
        x.iter()
            .zip(y)
            .fold(1i32, |sum, (&a, &b)| mac_16by16_to_int32(sum, a, b))
    };
    let sft = normalize_amr_wb(sum);
    *exp = 30 - sft;
    sum.wrapping_shl(sft as u32)
}

pub(crate) fn lg2_normalized(x: i32, exp: i16, exponent: &mut i16, fraction: &mut i16) {
    if x <= 0 {
        *exponent = 0;
        *fraction = 0;
        return;
    }
    *exponent = 30 - exp;
    let x = x >> 9;
    let i = (extract_h(x) - 32) as usize;
    let a = ((x >> 1) as i16) & 0x7fff;
    let y = l_deposit_h(LOG2_NORM_TABLE[i]);
    let tmp = LOG2_NORM_TABLE[i] - LOG2_NORM_TABLE[i + 1];
    *fraction = extract_h(msu_16by16_from_int32(y, tmp, a));
}

pub(crate) fn amrwb_log_2(x: i32, exponent: &mut i16, fraction: &mut i16) {
    let e = normalize_amr_wb(x);
    lg2_normalized(shl_int32(x, e), e, exponent, fraction);
}

/// Splits into (hi, lo) double precision.
#[inline(always)]
pub(crate) fn int32_to_dpf(x: i32) -> (i16, i16) {
    let hi = (x >> 16) as i16;
    let lo = (x.wrapping_sub((hi as i32) << 16) >> 1) as i16;
    (hi, lo)
}

/// `Mpy_32_16`: a 32-bit value (split as `L_Extract` does) times 16 bits, saturating.
/// opencore's `fxp_mul32_by_16b(x, n) << 1` is the same product without the saturation, and
/// it adds the products plainly where TS 26.173 uses `L_add`.
#[inline(always)]
pub(crate) fn mpy_32_16(x: i32, n: i16) -> i32 {
    let (hi, lo) = int32_to_dpf(x);
    mac_16by16_to_int32(mul_16by16_to_int32(hi, n), mult_int16(lo, n), 1)
}

pub(crate) fn mpy_dpf_32(hi1: i16, lo1: i16, hi2: i16, lo2: i16) -> i32 {
    let l = mul_16by16_to_int32(hi1, hi2);
    let l = mac_16by16_to_int32(l, mult_int16(hi1, lo2), 1);
    mac_16by16_to_int32(l, mult_int16(lo1, hi2), 1)
}
