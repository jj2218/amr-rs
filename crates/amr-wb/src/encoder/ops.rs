//! The TS 26.173 basic operators under their ETSI names (`basicop2.c`, `oper_32b.c`,
//! `math_op.c`), as the encoder's C reads them: aliases of `crate::basic_op` (spec semantics,
//! no overflow flag) plus the few the decoder never needed.

#[allow(unused_imports)] // the full set, as the C has it
pub(crate) use crate::basic_op::{
    add_int16 as add, add_int32 as l_add, amr_wb_round as round, div_16by16 as div_s, extract_h,
    int32_to_dpf as l_extract, l_deposit_h, mac_16by16_to_int32 as l_mac, mpy_dpf_32 as mpy_32,
    msu_16by16_from_int32 as l_msu, mul_16by16_to_int32 as l_mult, mult_int16 as mult,
    mult_int16_r as mult_r, negate_int16 as negate, norm_s, normalize_amr_wb as norm_l,
    one_ov_sqrt_norm as isqrt_n, power_of_2 as pow2, shl_int16 as shl, shl_int32 as l_shl,
    shr_int32 as l_shr, sub_int16 as sub, sub_int32 as l_sub,
};

/// `shr`: arithmetic shift right (saturating left shift for negative counts).
#[inline(always)]
pub(crate) fn shr(a: i16, b: i16) -> i16 {
    if b < 0 { shl(a, -b) } else { a >> b.min(15) }
}

/// `shr_r`: shift right with rounding.
pub(crate) fn shr_r(a: i16, b: i16) -> i16 {
    if b > 15 {
        return 0;
    }
    let r = shr(a, b);
    if b > 0 && (a >> (b - 1)) & 1 != 0 {
        r + 1
    } else {
        r
    }
}

#[inline(always)]
pub(crate) fn abs_s(a: i16) -> i16 {
    if a == i16::MIN { i16::MAX } else { a.abs() }
}

#[inline(always)]
pub(crate) fn l_abs(a: i32) -> i32 {
    if a == i32::MIN { i32::MAX } else { a.abs() }
}

#[inline(always)]
pub(crate) fn l_negate(a: i32) -> i32 {
    if a == i32::MIN { i32::MAX } else { -a }
}

#[inline(always)]
pub(crate) fn extract_l(a: i32) -> i16 {
    a as i16
}

/// `Mpy_32_16(hi, lo, n)`: double precision times 16 bits.
#[inline(always)]
pub(crate) fn mpy_32_16(hi: i16, lo: i16, n: i16) -> i32 {
    l_mac(l_mult(hi, n), mult(lo, n), 1)
}

/// `L_Comp(hi, lo)`: `hi·2^16 + lo·2`.
#[inline(always)]
pub(crate) fn l_comp(hi: i16, lo: i16) -> i32 {
    l_mac((hi as i32) << 16, lo, 1)
}

/// `Div_32(L_num, denom)`: `L_num / denom` in Q31 (0 <= L_num < denom, denom normalised).
pub(crate) fn div_32(l_num: i32, denom_hi: i16, denom_lo: i16) -> i32 {
    let approx = div_s(0x3fff, denom_hi);
    let l = l_sub(i32::MAX, mpy_32_16(denom_hi, denom_lo, approx));
    let (hi, lo) = l_extract(l);
    let (hi, lo) = l_extract(mpy_32_16(hi, lo, approx));
    let (n_hi, n_lo) = l_extract(l_num);
    l_shl(mpy_32(n_hi, n_lo, hi, lo), 2)
}

/// `Log2(L_x)`: (exponent, fraction).
pub(crate) fn log2(l_x: i32) -> (i16, i16) {
    let (mut e, mut f) = (0, 0);
    crate::basic_op::amrwb_log_2(l_x, &mut e, &mut f);
    (e, f)
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_shr(a: i16, b: i16) -> i16;
        fn refwb_Div_32(n: i32, hi: i16, lo: i16) -> i32;
        fn refwb_Log2(x: i32, e: *mut i16, f: *mut i16);
        fn refwb_div_s(a: i16, b: i16) -> i16;
        fn refwb_Mpy_32_16(hi: i16, lo: i16, n: i16) -> i32;
    }

    #[test]
    fn operators_match_reference() {
        let mut r = Rng::new(3);
        with_ref(|| {
            for _ in 0..200_000 {
                let (a, b, x) = (r.i16(), r.i16(), r.i32());
                let n = r.shift();
                // SAFETY: plain value functions.
                unsafe {
                    assert_eq!(shr(a, n), refwb_shr(a, n), "shr {a} {n}");
                    let (hi, lo) = l_extract(x);
                    assert_eq!(mpy_32_16(hi, lo, a), refwb_Mpy_32_16(hi, lo, a));
                    if x > 0 {
                        let (mut e, mut f) = (0, 0);
                        refwb_Log2(x, &mut e, &mut f);
                        assert_eq!(log2(x), (e, f), "Log2 {x}");
                        // Div_32: numerator below a normalised denominator
                        let d = x | 0x4000_0000;
                        let (dh, dl) = l_extract(d);
                        let num = (r.u64() % d as u64) as i32;
                        assert_eq!(div_32(num, dh, dl), refwb_Div_32(num, dh, dl), "Div_32");
                    }
                    let (p, q) = (abs_s(a), abs_s(b).max(1));
                    let (p, q) = (p.min(q), p.max(q));
                    assert_eq!(div_s(p, q), refwb_div_s(p, q), "div_s {p} {q}");
                }
            }
        });
    }
}
