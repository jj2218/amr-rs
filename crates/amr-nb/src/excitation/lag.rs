//! Pitch lag decoding (`dec_lag3.cpp`, `dec_lag6.cpp`).

use crate::basic_op::{Flag, mult, sub};

/// `Dec_lag3`: 1/3-resolution lag decoding (all modes but MR122). Returns `(T0, T0_frac)`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn dec_lag3(
    index: i16,
    t0_min: i16,
    t0_max: i16,
    i_subfr: i16,
    t0_prev: i16,
    flag4: i16,
    ov: &mut Flag,
) -> (i16, i16) {
    let t0;
    let t0_frac;
    if i_subfr == 0 {
        if index < 197 {
            let mut tmp_lag = index.wrapping_add(2);
            tmp_lag = mult(tmp_lag, 10923, ov);
            let mut i = tmp_lag.wrapping_add(19);
            t0 = i;
            i <<= 1;
            i = i.wrapping_add(t0);
            tmp_lag = index.wrapping_sub(i);
            t0_frac = tmp_lag.wrapping_add(58);
        } else {
            t0 = index.wrapping_sub(112);
            t0_frac = 0;
        }
    } else if flag4 == 0 {
        let mut i = index.wrapping_add(2);
        i = ((i as i32 * 10923) >> 15) as i16;
        i = i.wrapping_sub(1);
        t0 = i.wrapping_add(t0_min);
        i = i.wrapping_add(i << 1);
        let tmp_lag = index.wrapping_sub(2);
        t0_frac = tmp_lag.wrapping_sub(i);
    } else {
        let mut tmp_lag = t0_prev;
        let mut i = sub(tmp_lag, t0_min, ov);
        if i > 5 {
            tmp_lag = t0_min.wrapping_add(5);
        }
        i = t0_max.wrapping_sub(tmp_lag);
        if i > 4 {
            tmp_lag = t0_max.wrapping_sub(4);
        }
        if index < 4 {
            i = tmp_lag.wrapping_sub(5);
            t0 = i.wrapping_add(index);
            t0_frac = 0;
        } else if index < 12 {
            i = index.wrapping_sub(5);
            i = ((i as i32 * 10923) >> 15) as i16;
            i = i.wrapping_sub(1);
            t0 = i.wrapping_add(tmp_lag);
            i = i.wrapping_add(i << 1);
            tmp_lag = index.wrapping_sub(9);
            t0_frac = tmp_lag.wrapping_sub(i);
        } else {
            i = index.wrapping_sub(12);
            i = i.wrapping_add(tmp_lag);
            t0 = i.wrapping_add(1);
            t0_frac = 0;
        }
    }
    (t0, t0_frac)
}

/// `Dec_lag6`: 1/6-resolution lag decoding (MR122). `t0` is in/out as in C (the previous
/// subframe's integer lag is read when `i_subfr != 0`); returns `T0_frac`.
pub(crate) fn dec_lag6(
    index: i16,
    pit_min: i16,
    pit_max: i16,
    i_subfr: i16,
    t0: &mut i16,
    _ov: &mut Flag,
) -> i16 {
    if i_subfr == 0 {
        if index < 463 {
            let mut i = index.wrapping_add(5);
            i = ((i as i32 * 5462) >> 15) as i16;
            i = i.wrapping_add(17);
            *t0 = i;
            i <<= 1;
            i = i.wrapping_add(*t0);
            i <<= 1;
            i = index.wrapping_sub(i);
            i.wrapping_add(105)
        } else {
            *t0 = index.wrapping_sub(368);
            0
        }
    } else {
        let mut t0_min = t0.wrapping_sub(5);
        if t0_min < pit_min {
            t0_min = pit_min;
        }
        let mut t0_max = t0_min.wrapping_add(9);
        if t0_max > pit_max {
            t0_max = pit_max;
            t0_min = t0_max.wrapping_sub(9);
        }
        let mut i = index.wrapping_add(5);
        i = ((i as i32 * 5462) >> 15) as i16;
        i = i.wrapping_sub(1);
        *t0 = i.wrapping_add(t0_min);
        i = i.wrapping_add(i << 1);
        i <<= 1;
        let k = index.wrapping_sub(3);
        k.wrapping_sub(i)
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_c {
    use super::*;
    use crate::testutil::{N, Rng};
    use opencore_sys as _;

    // The TS 26.073 reference, in the shape of the opencore functions these tests were
    // written against (`ov`/`o`: its overflow flag; table arguments unused).
    unsafe extern "C" {
        fn ref_Dec_lag3(
            index: i16,
            t0_min: i16,
            t0_max: i16,
            i_subfr: i16,
            t0_prev: i16,
            t0: *mut i16,
            t0_frac: *mut i16,
            flag4: i16,
        );
        fn ref_Dec_lag6(
            index: i16,
            pit_min: i16,
            pit_max: i16,
            i_subfr: i16,
            t0: *mut i16,
            t0_frac: *mut i16,
        );
    }
    #[allow(clippy::too_many_arguments)]
    unsafe fn c_dec_lag3(
        index: i16,
        t0_min: i16,
        t0_max: i16,
        i_subfr: i16,
        t0_prev: i16,
        t0: *mut i16,
        t0_frac: *mut i16,
        flag4: i16,
        o: *mut i32,
    ) {
        crate::testutil::with_ref(o, || unsafe {
            ref_Dec_lag3(index, t0_min, t0_max, i_subfr, t0_prev, t0, t0_frac, flag4)
        })
    }
    unsafe fn c_dec_lag6(
        index: i16,
        pit_min: i16,
        pit_max: i16,
        i_subfr: i16,
        t0: *mut i16,
        t0_frac: *mut i16,
        o: *mut i32,
    ) {
        crate::testutil::with_ref(o, || unsafe {
            ref_Dec_lag6(index, pit_min, pit_max, i_subfr, t0, t0_frac)
        })
    }

    /// The decoder's search window around the previous lag (as in `Dec_lag3`'s callers).
    fn window(t0_prev: i16, pit_min: i16, pit_max: i16) -> (i16, i16) {
        let mut t0_min = (t0_prev - 5).max(pit_min);
        let mut t0_max = t0_min + 9;
        if t0_max > pit_max {
            t0_max = pit_max;
            t0_min = t0_max - 9;
        }
        (t0_min, t0_max)
    }

    #[test]
    fn dec_lag3_vs_c() {
        let mut g = Rng::new(10);
        for _ in 0..N {
            // i_subfr == 0: 8-bit absolute index; others: 5-bit (flag4 = 0) or 4-bit
            // (flag4 = 1, MR475/MR515/MR59/MR67) relative index around the previous lag.
            let i_subfr = [0i16, 40, 80, 120][(g.u64() % 4) as usize];
            let flag4 = (g.u64() % 2) as i16;
            let index = match (i_subfr, flag4) {
                (0, _) => (g.u64() % 256) as i16,
                (_, 0) => (g.u64() % 32) as i16,
                _ => (g.u64() % 16) as i16,
            };
            let t0_prev = 20 + (g.u64() % 124) as i16;
            let (t0_min, t0_max) = window(t0_prev, 20, 143);
            // T0_prev is a free input in Dec_lag3: also exercise values outside the window.
            let t0_prev = if g.u64().is_multiple_of(4) {
                20 + (g.u64() % 124) as i16
            } else {
                t0_prev
            };
            let mut ro = false;
            let r = dec_lag3(index, t0_min, t0_max, i_subfr, t0_prev, flag4, &mut ro);
            let (mut t0, mut fr, mut co) = (0, 0, 0);
            unsafe {
                c_dec_lag3(
                    index, t0_min, t0_max, i_subfr, t0_prev, &mut t0, &mut fr, flag4, &mut co,
                )
            };
            assert_eq!(
                (r, ro),
                ((t0, fr), co != 0),
                "Dec_lag3({index},{t0_min},{t0_max},{i_subfr},{t0_prev},{flag4})"
            );
        }
    }

    #[test]
    fn dec_lag6_vs_c() {
        let mut g = Rng::new(11);
        for _ in 0..N {
            // MR122: 9-bit absolute index for i_subfr == 0, 6-bit relative otherwise.
            let i_subfr = [0i16, 40, 80, 120][(g.u64() % 4) as usize];
            let index = if i_subfr == 0 {
                (g.u64() % 512) as i16
            } else {
                (g.u64() % 64) as i16
            };
            let t0_in = 18 + (g.u64() % 126) as i16;
            let (mut rt0, mut ro) = (t0_in, false);
            let rf = dec_lag6(index, 18, 143, i_subfr, &mut rt0, &mut ro);
            let (mut t0, mut fr, mut co) = (t0_in, 0, 0);
            unsafe { c_dec_lag6(index, 18, 143, i_subfr, &mut t0, &mut fr, &mut co) };
            assert_eq!(
                (rt0, rf, ro),
                (t0, fr, co != 0),
                "Dec_lag6({index},{i_subfr},{t0_in})"
            );
        }
    }
}
