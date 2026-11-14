//! Fixed-point math helpers: `Inv_sqrt`, `Log2`, `Log2_norm`, `Pow2`, `sqrt_l_exp`.
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `common/src/{inv_sqrt,log2,log2_norm,pow2,sqrt_l}.cpp`.
#![allow(dead_code)] // used progressively as decoder modules land

use crate::basic_op::{Flag, l_msu, l_mult, l_shl, l_shr_r, norm_l};
use crate::tables::{INV_SQRT_TBL, LOG2_TBL, POW2_TBL, SQRT_L_TBL};

/// `Inv_sqrt`: 1/sqrt(L_x), L_x > 0, result in Q30-ish as per the C code.
pub(crate) fn inv_sqrt(l_x: i32, _ov: &mut Flag) -> i32 {
    if l_x <= 0 {
        return 0x3fff_ffff;
    }
    let mut exp = norm_l(l_x);
    let mut x = l_x << exp;
    exp = 30 - exp;
    if exp & 1 == 0 {
        x >>= 1;
    }
    exp >>= 1;
    exp += 1;
    x >>= 9;
    let i = (x >> 16) as i16 - 16;
    let a = ((x >> 1) as i16) & 0x7fff;
    let i = i as usize;
    let mut y = (INV_SQRT_TBL[i] as i32) << 16;
    let tmp = INV_SQRT_TBL[i].wrapping_sub(INV_SQRT_TBL[i + 1]);
    y = y.wrapping_sub((tmp as i32 * a as i32) << 1);
    y >> exp
}

/// `Log2_norm`: returns `(exponent, fraction)` for an already-normalised `l_x`.
pub(crate) fn log2_norm(l_x: i32, exp: i16) -> (i16, i16) {
    if l_x <= 0 {
        return (0, 0);
    }
    let exponent = 30 - exp;
    let x = l_x >> 10;
    let i = ((x >> 15) as i16 - 32) as usize;
    let a = x & 0x7fff;
    let mut y = (LOG2_TBL[i] as i32) << 16;
    let tmp = LOG2_TBL[i].wrapping_sub(LOG2_TBL[i + 1]);
    y = y.wrapping_sub((tmp as i32 * a) << 1);
    (exponent, (y >> 16) as i16)
}

/// `Log2`: returns `(exponent, fraction)`.
pub(crate) fn log2(l_x: i32, _ov: &mut Flag) -> (i16, i16) {
    let exp = norm_l(l_x);
    log2_norm(l_x.wrapping_shl(exp as u32), exp)
}

/// `Pow2`: 2^(exponent.fraction).
pub(crate) fn pow2(exponent: i16, fraction: i16, ov: &mut Flag) -> i32 {
    let mut x = l_mult(fraction, 32, ov);
    let i = (((x >> 16) as i16) & 31) as usize;
    let a = ((x >> 1) & 0x7fff) as i16;
    x = (POW2_TBL[i] as i32) << 16;
    let tmp = POW2_TBL[i].wrapping_sub(POW2_TBL[i + 1]);
    x = l_msu(x, tmp, a, ov);
    l_shr_r(x, 30i16.wrapping_sub(exponent), ov)
}

/// `sqrt_l_exp`: returns `(sqrt, exp)`.
pub(crate) fn sqrt_l_exp(l_x: i32, ov: &mut Flag) -> (i32, i16) {
    if l_x <= 0 {
        return (0, 0);
    }
    let e = norm_l(l_x) & !1;
    let mut x = l_shl(l_x, e, ov);
    x >>= 10;
    let mut i = ((x >> 15) as i16) & 63;
    let a = (x as i16) & 0x7fff;
    if i > 15 {
        i -= 16;
    }
    let i = i as usize;
    let mut y = (SQRT_L_TBL[i] as i32) << 16;
    let tmp = SQRT_L_TBL[i].wrapping_sub(SQRT_L_TBL[i + 1]);
    y = l_msu(y, tmp, a, ov);
    (y, e)
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_c {
    use super::*;
    use crate::testutil::with_ref;
    use crate::testutil::{N, Rng};

    // The TS 26.073 reference, in the shape of the opencore functions these tests were written
    // against (`o`: its overflow flag).
    unsafe extern "C" {
        fn ref_Inv_sqrt(x: i32) -> i32;
        fn ref_Log2(x: i32, e: *mut i16, f: *mut i16);
        fn ref_Log2_norm(x: i32, exp: i16, e: *mut i16, f: *mut i16);
        fn ref_Pow2(e: i16, f: i16) -> i32;
        fn ref_sqrt_l_exp(x: i32, e: *mut i16) -> i32;
    }
    #[allow(non_snake_case)]
    unsafe fn Inv_sqrt(x: i32, o: *mut i32) -> i32 {
        with_ref(o, || unsafe { ref_Inv_sqrt(x) })
    }
    #[allow(non_snake_case)]
    unsafe fn Log2(x: i32, e: *mut i16, f: *mut i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_Log2(x, e, f) })
    }
    #[allow(non_snake_case)]
    unsafe fn Log2_norm(x: i32, exp: i16, e: *mut i16, f: *mut i16) {
        with_ref(std::ptr::null_mut(), || unsafe {
            ref_Log2_norm(x, exp, e, f)
        })
    }
    #[allow(non_snake_case)]
    unsafe fn Pow2(e: i16, f: i16, o: *mut i32) -> i32 {
        with_ref(o, || unsafe { ref_Pow2(e, f) })
    }
    unsafe fn c_sqrt_l_exp(x: i32, e: *mut i16, o: *mut i32) -> i32 {
        with_ref(o, || unsafe { ref_sqrt_l_exp(x, e) })
    }

    #[test]
    fn math_vs_c() {
        let mut g = Rng::new(42);
        for _ in 0..N {
            // Positive inputs (the functions' domain) plus the documented <= 0 branches.
            let x = if g.u64().is_multiple_of(16) {
                g.i32().min(0)
            } else {
                (g.i32() & 0x7fff_ffff).max(1)
            };
            let (mut o1, mut c1) = (false, 0);
            assert_eq!(
                inv_sqrt(x, &mut o1),
                unsafe { Inv_sqrt(x, &mut c1) },
                "Inv_sqrt({x})"
            );

            let (mut e, mut f, mut c) = (0, 0, 0);
            unsafe { Log2(x, &mut e, &mut f, &mut c) };
            assert_eq!(log2(x, &mut o1), (e, f), "Log2({x})");

            let n = norm_l(x);
            let xn = x.wrapping_shl(n as u32);
            unsafe { Log2_norm(xn, n, &mut e, &mut f) };
            assert_eq!(log2_norm(xn, n), (e, f), "Log2_norm({xn}, {n})");

            let (mut q, mut qc) = (0, 0);
            let mut qo = false;
            let r = sqrt_l_exp(x, &mut qo);
            let c = unsafe { c_sqrt_l_exp(x, &mut q, &mut qc) };
            assert_eq!((r, qo), ((c, q), qc != 0), "sqrt_l_exp({x})");

            // Pow2 domain: exponent 0..=30, fraction 0..=32767.
            let (pe, pf) = ((g.u64() % 31) as i16, (g.i16() & 0x7fff));
            let (mut po, mut pc) = (false, 0);
            let r = pow2(pe, pf, &mut po);
            let c = unsafe { Pow2(pe, pf, &mut pc) };
            assert_eq!((r, po), (c, pc != 0), "Pow2({pe}, {pf})");
        }
    }
}
