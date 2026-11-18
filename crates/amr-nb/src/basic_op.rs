//! Fixed-point basic operators, bit-exact with opencore-amr's C equivalents
//! (`basic_op_c_equivalent.h`, `basic_op.h` and the `common/src/*.cpp` operator files).
//!
//! Semantics deliberately follow opencore, quirks included (e.g. `shl` by a large negative
//! amount yields 0 even for negative input; `l_mac` checks overflow with `> 0`). Shifts use
//! wrapping counts, matching the x86 behaviour of the C build. `ov` is opencore's
//! `Flag *pOverflow`: operators only ever set it, never clear it.
//! Derived from opencore-amr 0.1.6 (Apache-2.0).

/// opencore's `Flag`.
pub type Flag = bool;

pub const MAX_32: i32 = i32::MAX;
pub const MIN_32: i32 = i32::MIN;
pub const MAX_16: i16 = i16::MAX;
pub const MIN_16: i16 = i16::MIN;

/// `add_16`
#[inline(always)]
pub fn add(var1: i16, var2: i16, ov: &mut Flag) -> i16 {
    let sum = var1 as i32 + var2 as i32;
    if sum > 0x7fff {
        *ov = true;
        MAX_16
    } else if sum < -32768 {
        *ov = true;
        MIN_16
    } else {
        sum as i16
    }
}

/// `sub`
#[inline(always)]
pub fn sub(var1: i16, var2: i16, ov: &mut Flag) -> i16 {
    let diff = var1 as i32 - var2 as i32;
    if (diff.wrapping_add(32768) as u32) > 0xffff {
        *ov = true;
        if diff > 0x7fff { MAX_16 } else { MIN_16 }
    } else {
        diff as i16
    }
}

/// `shl`, as TS 26.073 `basicop2.c` (a negative count is `shr`, so it keeps the sign; a count
/// above 15 saturates any non-zero input). Never sets `ov`.
#[inline(always)]
pub fn shl(var1: i16, var2: i16, ov: &mut Flag) -> i16 {
    if var2 < 0 {
        shr(var1, -(var2.max(-16)), ov)
    } else if var2 > 15 {
        if var1 == 0 { 0 } else { (var1 >> 15) ^ MAX_16 }
    } else {
        let out = (var1 as i32).wrapping_shl(var2 as u32) as i16;
        if (out as i32).wrapping_shr(var2 as u32) != var1 as i32 {
            (var1 >> 15) ^ MAX_16
        } else {
            out
        }
    }
}

/// `shr`
#[inline(always)]
pub fn shr(var1: i16, var2: i16, ov: &mut Flag) -> i16 {
    if var2 == 0 {
        var1
    } else if var2 > 0 {
        var1 >> var2.min(15)
    } else {
        // TS 26.073: shl(var1, -var2) with the count clamped to 16; above 15 any non-zero
        // input saturates (and sets the flag).
        let s = (-(var2 as i32)).min(16) as u32;
        let result = (var1 as i32) << s;
        if result != result as i16 as i32 {
            *ov = true;
            if var1 > 0 { MAX_16 } else { MIN_16 }
        } else {
            result as i16
        }
    }
}

/// `shr_r`
#[inline(always)]
pub fn shr_r(var1: i16, var2: i16, ov: &mut Flag) -> i16 {
    if var2 > 15 {
        return 0;
    }
    let mut out = shr(var1, var2, ov);
    if var2 > 0 && (var1 as i32 & (1i32 << (var2 - 1))) != 0 {
        out = out.wrapping_add(1);
    }
    out
}

/// `mult`
#[inline(always)]
pub fn mult(var1: i16, var2: i16, ov: &mut Flag) -> i16 {
    let product = (var1 as i32 * var2 as i32) >> 15;
    if product > 0x7fff {
        *ov = true;
        MAX_16
    } else {
        product as i16
    }
}

/// `mult_r`
#[inline(always)]
pub fn mult_r(var1: i16, var2: i16, ov: &mut Flag) -> i16 {
    let mut p = var1 as i32 * var2 as i32;
    p += 0x4000;
    p >>= 15;
    p |= -(p & 0x10000);
    if p > 0x7fff {
        *ov = true;
        MAX_16
    } else if p < -32768 {
        *ov = true;
        MIN_16
    } else {
        p as i16
    }
}

/// `L_mult`
#[inline(always)]
pub fn l_mult(var1: i16, var2: i16, ov: &mut Flag) -> i32 {
    let p = var1 as i32 * var2 as i32;
    if p != 0x4000_0000 {
        p << 1
    } else {
        *ov = true;
        MAX_32
    }
}

/// `L_add`
#[inline(always)]
pub fn l_add(l_var1: i32, l_var2: i32, ov: &mut Flag) -> i32 {
    // Branchless: the overflow test is data-dependent and mispredicts as a short-circuit `&&`.
    let (_, o) = l_var1.overflowing_add(l_var2);
    *ov |= o;
    l_var1.saturating_add(l_var2)
}

/// `L_sub`
#[inline(always)]
pub fn l_sub(l_var1: i32, l_var2: i32, ov: &mut Flag) -> i32 {
    let diff = l_var1.wrapping_sub(l_var2);
    if (l_var1 ^ l_var2) < 0 && ((diff ^ l_var1) & MIN_32) != 0 {
        *ov = true;
        if l_var1 < 0 { MIN_32 } else { MAX_32 }
    } else {
        diff
    }
}

/// `L_mac` = `L_add(l_var3, L_mult(var1, var2))`, as TS 26.073 (opencore skips the
/// saturation when `l_var3` equals the product, and returns `MAX_32` for -32768²).
#[inline(always)]
pub fn l_mac(l_var3: i32, var1: i16, var2: i16, ov: &mut Flag) -> i32 {
    let p = l_mult(var1, var2, ov);
    l_add(l_var3, p, ov)
}

/// Energy: `xs.fold(0, |s, x| l_mac(s, x, x, ov))`, bit-exact. Fast path: the terms are
/// non-negative, so if no `x` is -32768 (the `0x4000_0000` product case) and the exact total
/// fits in 32 bits, no partial sum saturated and the plain (vectorisable) sum is the answer.
/// Otherwise the faithful fold runs.
#[inline]
pub fn l_mac_energy(xs: &[i16], ov: &mut Flag) -> i32 {
    l_mac_energy_from(0, xs, ov)
}

/// [`l_mac_energy`] starting from `init >= 0` (`xs.fold(init, ..)`).
#[inline]
pub fn l_mac_energy_from(init: i32, xs: &[i16], ov: &mut Flag) -> i32 {
    debug_assert!(init >= 0);
    let (sum, any_min) = xs.iter().fold((init as i64, false), |(s, m), &x| {
        (s + 2 * (x as i32 * x as i32) as i64, m | (x == MIN_16))
    });
    if !any_min && sum <= MAX_32 as i64 {
        sum as i32
    } else {
        xs.iter().fold(init, |s, &x| l_mac(s, x, x, ov))
    }
}

/// max |v| (0 for an empty slice), for the fast-path bounds: folded as i16 min/max
/// (`pminsw`/`pmaxsw`; i64 or i32 abs/max of widened lanes stayed scalar or half-width).
#[inline(always)]
pub fn max_abs(v: &[i16]) -> i64 {
    let (lo, hi) = v
        .iter()
        .fold((0i16, 0i16), |(l, h), &c| (l.min(c), h.max(c)));
    (hi as i64).max(-(lo as i64))
}

/// sum |v|, for the fast-path bounds: an i32 sum (exact up to 65,536 elements), widened once.
#[inline(always)]
pub fn sum_abs(v: &[i16]) -> i64 {
    debug_assert!(v.len() <= 1 << 16);
    v.iter().map(|&c| (c as i32).abs()).sum::<i32>() as i64
}

/// `a.zip(b).fold(init, |s, (x, y)| l_mac(s, x, y, ov))`, bit-exact. Fast path: if
/// `|init| + 2 sum |x y| <= MAX_32`, no partial sum (nor product) saturates and the plain
/// (vectorisable) sum is the answer; otherwise the faithful fold.
#[inline]
pub fn l_mac_dot(init: i32, a: &[i16], b: &[i16], ov: &mut Flag) -> i32 {
    let (sum, abs) = a.iter().zip(b).fold((0i64, 0i64), |(s, m), (&x, &y)| {
        let p = x as i64 * y as i64;
        (s + p, m + p.abs())
    });
    if (init as i64).abs() + 2 * abs <= MAX_32 as i64 {
        (init as i64 + 2 * sum) as i32
    } else {
        a.iter().zip(b).fold(init, |s, (&x, &y)| l_mac(s, x, y, ov))
    }
}

/// Energy as `Bgn_scd` and `dtx_dec_activity_update` compute it: an `L_add` chain over
/// x² << 1 per sample, `MAX_32` (with no overflow flag) for −32768's 0x4000_0000 product. The
/// terms are ≥ 0, so a prefix sum first passes `MAX_32` exactly when the total does and then
/// stays saturated: the chain is min(Σ, MAX_32), with ov iff Σ > MAX_32.
/// Σ = Σ 2x² − #(x = −32768) (2^31 vs `MAX_32`): one vectorising i64 sum.
#[inline]
pub fn l_add_energy(xs: &[i16], ov: &mut Flag) -> i32 {
    let (sq, mins) = xs.iter().fold((0i64, 0i64), |(s, m), &x| {
        (s + (x as i32 * x as i32) as i64, m + (x == MIN_16) as i64)
    });
    let sum = 2 * sq - mins;
    *ov |= sum > MAX_32 as i64;
    sum.min(MAX_32 as i64) as i32
}

/// `L_msu`
#[inline(always)]
pub fn l_msu(l_var3: i32, var1: i16, var2: i16, ov: &mut Flag) -> i32 {
    let p = l_mult(var1, var2, ov);
    l_sub(l_var3, p, ov)
}

/// `L_shl`, as TS 26.073 (a negative count is `L_shr`; above 31 saturates). Never sets `ov`.
#[inline(always)]
pub fn l_shl(l_var1: i32, var2: i16, _ov: &mut Flag) -> i32 {
    if var2 > 31 {
        if l_var1 == 0 {
            0
        } else {
            (l_var1 >> 31) ^ MAX_32
        }
    } else if var2 > 0 {
        let out = l_var1.wrapping_shl(var2 as u32);
        if out.wrapping_shr(var2 as u32) != l_var1 {
            (l_var1 >> 31) ^ MAX_32
        } else {
            out
        }
    } else {
        l_var1 >> (-(var2 as i32)).min(31)
    }
}

/// `L_shr`, as TS 26.073 (31 or more gives the sign; a negative count is `L_shl`). Never sets
/// `ov`.
#[inline(always)]
pub fn l_shr(l_var1: i32, var2: i16, ov: &mut Flag) -> i32 {
    if var2 >= 0 {
        l_var1 >> var2.min(31)
    } else {
        l_shl(l_var1, -(var2.max(-32)), ov)
    }
}

/// `L_shr_r`
#[inline(always)]
pub fn l_shr_r(l_var1: i32, var2: i16, ov: &mut Flag) -> i32 {
    if var2 > 31 {
        return 0;
    }
    let mut r = l_shr(l_var1, var2, ov);
    if var2 > 0 && (l_var1 & (1i32 << (var2 - 1))) != 0 {
        r = r.wrapping_add(1);
    }
    r
}

/// `negate`
/// The C idiom `x < 0 ? ~((~x) >> n) : x >> n`: an arithmetic shift that rounds negative values
/// towards zero (for an `i16` input the result fits an `i16`).
#[inline(always)]
pub fn neg_shr(x: i32, n: u32) -> i32 {
    if x < 0 { !((!x) >> n) } else { x >> n }
}

#[inline(always)]
pub fn negate(var1: i16) -> i16 {
    if var1 == MIN_16 { MAX_16 } else { -var1 }
}

/// `L_negate`
#[inline(always)]
pub fn l_negate(l_var1: i32) -> i32 {
    if l_var1 == MIN_32 { MAX_32 } else { -l_var1 }
}

/// `abs_s` (opencore's branch-free variant: `abs_s(MIN_16) == MAX_16`).
#[inline(always)]
pub fn abs_s(var1: i16) -> i16 {
    let y = var1.wrapping_sub((var1 < 0) as i16);
    y ^ (y >> 15)
}

/// `norm_s`. Negative inputs normalise `~var1` like TS 26.073 `basicop2.c` (opencore's `-var1`
/// gives one less for -2^k, e.g. 14 instead of 15 for -1).
#[inline(always)]
pub fn norm_s(var1: i16) -> i16 {
    if var1 == 0 {
        return 0;
    }
    ((var1 ^ (var1 >> 15)).leading_zeros() as i16) - 1
}

/// `norm_l`. Negative inputs normalise `~l_var1` like TS 26.073 `basicop2.c` (opencore's
/// `-l_var1` gives one less for -2^k, e.g. 30 instead of 31 for -1).
#[inline(always)]
pub fn norm_l(l_var1: i32) -> i16 {
    if l_var1 == 0 {
        return 0;
    }
    ((l_var1 ^ (l_var1 >> 31)).leading_zeros() as i16) - 1
}

/// `pv_round` (`round` in the ETSI code).
#[inline(always)]
pub fn round(l_var1: i32, ov: &mut Flag) -> i16 {
    (l_add(l_var1, 0x8000, ov) >> 16) as i16
}

/// `div_s`
pub fn div_s(var1: i16, var2: i16) -> i16 {
    if var1 > var2 || var1 < 0 {
        return 0;
    }
    if var1 == 0 {
        return 0;
    }
    if var1 == var2 {
        return MAX_16;
    }
    // The C's 15-step restoring division computes floor(var1 * 2^15 / var2) for
    // 0 < var1 < var2; one hardware divide gives the same quotient.
    (((var1 as i32) << 15) / var2 as i32) as i16
}

/// `extract_h`
#[inline(always)]
pub fn extract_h(l_var1: i32) -> i16 {
    (l_var1 >> 16) as i16
}

/// `L_Extract`: returns `(hi, lo)`.
#[inline(always)]
pub fn l_extract(l_var: i32, _ov: &mut Flag) -> (i16, i16) {
    let temp = l_var >> 16;
    (temp as i16, ((l_var >> 1) - (temp << 15)) as i16)
}

/// `L_Comp`
#[inline(always)]
pub fn l_comp(hi: i16, lo: i16, ov: &mut Flag) -> i32 {
    l_mac((hi as i32) << 16, lo, 1, ov)
}

/// `Mpy_32`. Never sets `ov`.
#[inline(always)]
pub fn mpy_32(hi1: i16, lo1: i16, hi2: i16, lo2: i16, _ov: &mut Flag) -> i32 {
    let mut lp = hi1 as i32 * hi2 as i32;
    lp = if lp != 0x4000_0000 { lp << 1 } else { MAX_32 };
    let sat_add = |a: i32, p32: i32| -> i32 {
        let s = a.wrapping_add(p32 << 1);
        if (a ^ p32) > 0 && (s ^ a) < 0 {
            if a < 0 { MIN_32 } else { MAX_32 }
        } else {
            s
        }
    };
    lp = sat_add(lp, (hi1 as i32 * lo2 as i32) >> 15);
    sat_add(lp, (lo1 as i32 * hi2 as i32) >> 15)
}

/// `Mpy_32_16`
#[inline(always)]
pub fn mpy_32_16(hi: i16, lo: i16, var2: i16, ov: &mut Flag) -> i32 {
    let mut lp = hi as i32 * var2 as i32;
    if lp != 0x4000_0000 {
        lp <<= 1;
    } else {
        *ov = true;
        lp = MAX_32;
    }
    let r = (lo as i32 * var2 as i32) >> 15;
    let s = lp.wrapping_add(r << 1);
    if (lp ^ r) > 0 && (s ^ lp) < 0 {
        *ov = true;
        if lp < 0 { MIN_32 } else { MAX_32 }
    } else {
        s
    }
}

/// `Mac_32`
#[inline(always)]
pub fn mac_32(l_var3: i32, hi1: i16, lo1: i16, hi2: i16, lo2: i16, ov: &mut Flag) -> i32 {
    let mut acc = l_mac(l_var3, hi1, hi2, ov);
    let p = mult(hi1, lo2, ov);
    acc = l_mac(acc, p, 1, ov);
    let p = mult(lo1, hi2, ov);
    l_mac(acc, p, 1, ov)
}

/// `Mac_32_16`
#[inline(always)]
pub fn mac_32_16(l_var3: i32, hi: i16, lo: i16, var2: i16, ov: &mut Flag) -> i32 {
    let acc = l_mac(l_var3, hi, var2, ov);
    let p = mult(lo, var2, ov);
    l_mac(acc, p, 1, ov)
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    //! Differential tests against the TS 26.073 reference (`ref26073-sys`); opencore only where
    //! the reference has no counterpart (`div_s` out of contract: the reference aborts).
    use super::*;
    use crate::testutil::{N, Rng};
    use opencore_sys as _;
    use ref26073_sys::{clear_overflow, lock, overflow};

    unsafe extern "C" {
        fn ref_add(a: i16, b: i16) -> i16;
        fn ref_sub(a: i16, b: i16) -> i16;
        fn ref_shl(a: i16, b: i16) -> i16;
        fn ref_shr(a: i16, b: i16) -> i16;
        fn ref_shr_r(a: i16, b: i16) -> i16;
        fn ref_mult(a: i16, b: i16) -> i16;
        fn ref_mult_r(a: i16, b: i16) -> i16;
        fn ref_L_mult(a: i16, b: i16) -> i32;
        fn ref_L_add(a: i32, b: i32) -> i32;
        fn ref_L_sub(a: i32, b: i32) -> i32;
        fn ref_L_mac(a: i32, b: i16, c: i16) -> i32;
        fn ref_L_msu(a: i32, b: i16, c: i16) -> i32;
        fn ref_L_shl(a: i32, b: i16) -> i32;
        fn ref_L_shr(a: i32, b: i16) -> i32;
        fn ref_L_shr_r(a: i32, b: i16) -> i32;
        fn ref_L_Comp(hi: i16, lo: i16) -> i32;
        fn ref_L_Extract(a: i32, hi: *mut i16, lo: *mut i16);
        fn ref_norm_s(a: i16) -> i16;
        fn ref_norm_l(a: i32) -> i16;
        fn ref_round(a: i32) -> i16;
        fn ref_div_s(a: i16, b: i16) -> i16;
        fn ref_negate(a: i16) -> i16;
        fn ref_abs_s(a: i16) -> i16;
        fn ref_L_negate(a: i32) -> i32;
        fn ref_Mpy_32(a: i16, b: i16, c: i16, d: i16) -> i32;
        fn ref_Mpy_32_16(a: i16, b: i16, c: i16) -> i32;
        fn ref_Mac_32(x: i32, a: i16, b: i16, c: i16, d: i16) -> i32;
        fn ref_Mac_32_16(x: i32, a: i16, b: i16, c: i16) -> i32;
        // opencore, for `div_s` outside its contract
        #[link_name = "div_s"]
        fn oc_div_s(a: i16, b: i16) -> i16;
    }

    /// Calls the reference with a cleared global flag; returns (result, flag).
    fn rf<T>(f: impl FnOnce() -> T) -> (T, bool) {
        clear_overflow();
        let r = f();
        (r, overflow())
    }

    /// Runs `N` random cases; `f` returns (rust, c) results incl. overflow flags.
    fn diff<T: PartialEq + std::fmt::Debug>(
        name: &str,
        mut f: impl FnMut(&mut Rng) -> (T, T, String),
    ) {
        let _g = lock();
        let mut rng = Rng::new(name.len() as u64 * 7919);
        for _ in 0..N {
            let (r, c, args) = f(&mut rng);
            assert_eq!(r, c, "{name}({args}): rust vs c");
        }
    }

    /// `$flag`: whether the Rust op tracks the overflow flag (`shl`, `L_shl`, `L_shr` don't:
    /// no caller reads it after them), compared with the reference's global.
    macro_rules! op2 {
        ($test:ident, $rs:expr, $c:ident, $a:ident, $b:ident, $flag:expr) => {
            #[test]
            fn $test() {
                diff(stringify!($test), |g| {
                    let (a, b) = (g.$a(), g.$b());
                    let mut ro = false;
                    let r = $rs(a, b, &mut ro);
                    let (c, co) = rf(|| unsafe { $c(a, b) });
                    ((r, ro && $flag), (c, co && $flag), format!("{a}, {b}"))
                });
            }
        };
    }
    op2!(add_vs_ref, add, ref_add, i16, i16, true);
    op2!(sub_vs_ref, sub, ref_sub, i16, i16, true);
    op2!(shl_vs_ref, shl, ref_shl, i16, shift, false);
    op2!(shr_vs_ref, shr, ref_shr, i16, shift, true);
    op2!(shr_r_vs_ref, shr_r, ref_shr_r, i16, shift, true);
    op2!(mult_vs_ref, mult, ref_mult, i16, i16, true);
    op2!(mult_r_vs_ref, mult_r, ref_mult_r, i16, i16, true);
    op2!(l_mult_vs_ref, l_mult, ref_L_mult, i16, i16, true);
    op2!(l_add_vs_ref, l_add, ref_L_add, i32, i32, true);
    op2!(l_sub_vs_ref, l_sub, ref_L_sub, i32, i32, true);
    op2!(l_shl_vs_ref, l_shl, ref_L_shl, i32, shift, false);
    op2!(l_shr_vs_ref, l_shr, ref_L_shr, i32, shift, false);
    op2!(l_shr_r_vs_ref, l_shr_r, ref_L_shr_r, i32, shift, false);
    op2!(l_comp_vs_ref, l_comp, ref_L_Comp, i16, i16, false);

    #[test]
    fn mac_msu_vs_ref() {
        let mut n = 0u32;
        diff("mac", |g| {
            let (b, c) = (g.i16(), g.i16());
            // Every 4th case: L_var3 == product, where opencore's `> 0` test skipped the
            // saturation (3 * product overflows for |product| > 2^31 / 3).
            n += 1;
            let a = if n.is_multiple_of(4) {
                b as i32 * c as i32
            } else {
                g.i32()
            };
            let (mut r1, mut r2) = (false, false);
            let r = (l_mac(a, b, c, &mut r1), r1, l_msu(a, b, c, &mut r2), r2);
            let (m, c1) = rf(|| unsafe { ref_L_mac(a, b, c) });
            let (s, c2) = rf(|| unsafe { ref_L_msu(a, b, c) });
            (r, (m, c1, s, c2), format!("{a}, {b}, {c}"))
        });
    }

    /// `l_add_energy` vs the C's loop (`Bgn_scd`, `dtx_dec_activity_update`) with its `L_add`.
    /// Classes counted: −32768 samples (the `MAX_32` term), the chain saturating, the total
    /// exactly `MAX_32` (no overflow), unsaturated totals above 0x1fff_ffff (`Bgn_scd`'s
    /// `MAX_16` branch).
    #[test]
    fn l_add_energy_vs_ref() {
        let _g = lock();
        let mut g = Rng::new(62);
        let mut cls = [0u32; 4];
        for case in 0..N / 4 {
            let amp = [0x7fff, 0x3fff, 0x0fff, 0x1fff][case % 4];
            let mut x: [i16; 160] = core::array::from_fn(|_| (g.i16() as i32 & amp) as i16);
            if case % 8 == 0 {
                x[(g.u64() % 160) as usize] = MIN_16;
            }
            if case % 16 == 1 {
                // one −32768, the rest 0: the total is MAX_32 and nothing overflows; or the
                // −32768 first, then one more term
                x = [0; 160];
                x[(g.u64() % 160) as usize] = MIN_16;
                if g.u64() & 1 == 0 {
                    x[(g.u64() % 160) as usize] = g.i16();
                }
            }
            let mut c = 0i32;
            clear_overflow();
            for i in (0..160).rev() {
                let p = x[i] as i32 * x[i] as i32;
                let t = if p != 0x4000_0000 { p << 1 } else { MAX_32 };
                c = unsafe { ref_L_add(c, t) };
            }
            let co = overflow() as i32;
            cls[0] += x.contains(&MIN_16) as u32;
            cls[1] += (co != 0) as u32;
            cls[2] += (c == MAX_32 && co == 0) as u32;
            cls[3] += (co == 0 && c > 0x1fff_ffff && c < MAX_32) as u32;
            let mut ro = case % 3 == 0;
            let r = l_add_energy(&x, &mut ro);
            assert_eq!((r, ro), (c, co != 0 || case % 3 == 0), "case {case}: {x:?}");
        }
        assert!(cls.iter().all(|&n| n > 500), "classes {cls:?}");
    }

    #[test]
    fn l_mac_energy_vs_fold() {
        let mut g = Rng::new(42);
        for n in 0..N {
            let len = (g.u64() % 41) as usize;
            let amp = [0x7fff, 0x3fff, 0x0fff][n % 3];
            let mut xs: Vec<i16> = (0..len).map(|_| (g.i16() as i32 & amp) as i16).collect();
            match n % 8 {
                // opencore's l_var3 == product quirk: [a, a, 2a] (an ordinary sum in TS 26.073).
                0 => xs.extend([13378, -13378, 26756]),
                1 if len > 0 => xs[len / 2] = MIN_16,
                _ => {}
            }
            let (mut fo, mut eo) = (n % 5 == 0, n % 5 == 0);
            let init = [0, 256, 1][n % 3];
            let f = xs.iter().fold(init, |s, &x| l_mac(s, x, x, &mut fo));
            let e = l_mac_energy_from(init, &xs, &mut eo);
            assert_eq!((e, eo), (f, fo), "{init} {xs:?}");
        }
    }

    #[test]
    fn unary_vs_ref() {
        diff("unary", |g| {
            let (a, l) = (g.i16(), g.i32());
            let mut ro = false;
            let r = (
                negate(a),
                abs_s(a),
                norm_s(a),
                norm_l(l),
                l_negate(l),
                round(l, &mut ro),
                ro,
            );
            let (rc, co) = rf(|| unsafe { ref_round(l) });
            let c = unsafe {
                (
                    ref_negate(a),
                    ref_abs_s(a),
                    ref_norm_s(a),
                    ref_norm_l(l),
                    ref_L_negate(l),
                    rc,
                    co,
                )
            };
            (r, c, format!("{a}, {l}"))
        });
    }

    /// In contract (0 <= a <= b, b > 0) against the reference; outside it (the reference
    /// aborts) against opencore, whose behaviour the port keeps.
    #[test]
    fn div_s_vs_ref() {
        diff("div_s", |g| {
            let (mut a, mut b) = (g.i16(), g.i16());
            let inside = g.u64() % 4 != 0;
            if inside {
                (a, b) = (a.saturating_abs(), b.saturating_abs().max(1));
                if a > b {
                    std::mem::swap(&mut a, &mut b);
                }
            }
            let inside = inside || (0 <= a && a <= b && b > 0);
            let c = unsafe {
                if inside {
                    ref_div_s(a, b)
                } else {
                    oc_div_s(a, b)
                }
            };
            (div_s(a, b), c, format!("{a}, {b}"))
        });
    }

    #[test]
    fn extract_vs_ref() {
        diff("extract", |g| {
            let l = g.i32();
            let mut o = false;
            let (mut hi, mut lo) = (0, 0);
            unsafe { ref_L_Extract(l, &mut hi, &mut lo) };
            (l_extract(l, &mut o), (hi, lo), format!("{l}"))
        });
    }

    #[test]
    fn mpy_mac_32_vs_ref() {
        diff("mpy32", |g| {
            let (x, a, b, c, d) = (
                g.i32(),
                g.i16(),
                g.i16() & 0x7fff,
                g.i16(),
                g.i16() & 0x7fff,
            );
            let (mut o1, mut o2, mut o3, mut o4) = (false, false, false, false);
            let r = (
                mpy_32(a, b, c, d, &mut o1),
                mpy_32_16(a, b, c, &mut o2),
                o2,
                mac_32(x, a, b, c, d, &mut o3),
                o3,
                mac_32_16(x, a, b, c, &mut o4),
                o4,
            );
            let m = unsafe { ref_Mpy_32(a, b, c, d) };
            let (m16, q2) = rf(|| unsafe { ref_Mpy_32_16(a, b, c) });
            let (mc, q3) = rf(|| unsafe { ref_Mac_32(x, a, b, c, d) });
            let (mc16, q4) = rf(|| unsafe { ref_Mac_32_16(x, a, b, c) });
            let cc = (m, m16, q2, mc, q3, mc16, q4);
            let _ = o1;
            (r, cc, format!("{x}, {a}, {b}, {c}, {d}"))
        });
    }
}
