//! `agc`, `agc2` and their energy helpers (opencore `dec/src/agc.cpp`), with TS 26.073's
//! saturating basic operators throughout (opencore's plain C arithmetic in `agc` wraps and
//! truncates).

use crate::basic_op::{
    Flag, MAX_16, MAX_32, MIN_16, MIN_32, add, div_s, l_mac_energy, l_mult, l_shl, l_shr, mult,
    norm_l, round, sub,
};
use crate::cnst::L_FRAME;
use crate::math::inv_sqrt;

/// `agcState`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct AgcState {
    pub(crate) past_gain: i16,
}

impl Default for AgcState {
    fn default() -> Self {
        Self { past_gain: 4096 }
    }
}

impl AgcState {
    /// `agc_reset`.
    pub(crate) fn reset(&mut self) {
        self.past_gain = 4096;
    }
}

/// `energy_old` (static in C): energy of `in_[..l_trm] >> 2`.
fn energy_old(in_: &[i16], l_trm: usize, ov: &mut Flag) -> i32 {
    let mut t = [0i16; L_FRAME];
    for (t, &v) in t.iter_mut().zip(&in_[..l_trm]) {
        *t = v >> 2;
    }
    l_mac_energy(&t[..l_trm], ov)
}

/// `energy_new` (static in C): energy of `in_[..l_trm]` scaled by 1/16, falling back to
/// `energy_old` (with the overflow flag restored first) if the sum saturates.
fn energy_new(in_: &[i16], l_trm: usize, ov: &mut Flag) -> i32 {
    let ov_save = *ov;
    let s = l_mac_energy(&in_[..l_trm], ov);
    if s != MAX_32 {
        s >> 4
    } else {
        *ov = ov_save;
        energy_old(in_, l_trm, ov)
    }
}

/// `agc`: scales `sig_out[..l_trm]` to the energy of `sig_in[..l_trm]`, smoothing the gain with
/// `agc_fac` across calls.
pub(crate) fn agc(
    st: &mut AgcState,
    sig_in: &[i16],
    sig_out: &mut [i16],
    agc_fac: i16,
    l_trm: usize,
    ov: &mut Flag,
) {
    let mut s = energy_new(sig_out, l_trm, ov);
    if s == 0 {
        st.past_gain = 0;
        return;
    }
    let mut exp = norm_l(s) - 1;
    let mut l_temp = l_shl(s, exp, ov);
    let gain_out = round(l_temp, ov);

    s = energy_new(sig_in, l_trm, ov);
    let g0: i16 = if s == 0 {
        0
    } else {
        let i = norm_l(s);
        l_temp = l_shl(s, i, ov);
        let gain_in = round(l_temp, ov);
        exp -= i;
        s = div_s(gain_out, gain_in) as i32;
        s = l_shl(s, 7, ov);
        s = l_shr(s, exp, ov);
        s = inv_sqrt(s, ov);
        let i = round(l_shl(s, 9, ov), ov);
        mult(i, sub(32767, agc_fac, ov), ov)
    };

    let mut gain = st.past_gain;
    if agc_fac >= 0 && g0 >= 0 {
        // Every new gain is in -32767..=32767 (gain fac >> 15 >= -32767 for fac <= 32767): `mult`
        // can't saturate and `add` only upwards, a min (flagged). The serial recursion alone,
        // then the scaling as lanes: `l_mult` can't saturate (gain > -32768), `L_shl(.., 3)`
        // clamps 16 v gain, i.e. v gain at +-2^27, then `>> 16`.
        let mut gains = [0i16; L_FRAME];
        let mut sat = false;
        for g in &mut gains[..l_trm] {
            let s = ((gain as i32 * agc_fac as i32) >> 15) + g0 as i32;
            sat |= s > MAX_16 as i32;
            gain = s.min(MAX_16 as i32) as i16;
            *g = gain;
        }
        *ov |= sat;
        for (v, &g) in sig_out[..l_trm].iter_mut().zip(&gains) {
            let p = *v as i32 * g as i32;
            *v = (p.clamp(-(1 << 27), (1 << 27) - 1) >> 12) as i16;
        }
    } else {
        for v in &mut sig_out[..l_trm] {
            gain = mult(gain, agc_fac, ov);
            gain = add(gain, g0, ov);
            *v = (l_shl(l_mult(*v, gain, ov), 3, ov) >> 16) as i16;
        }
    }
    st.past_gain = gain;
}

/// `agc2`: like `agc` but without smoothing or state (gain applied directly, saturating).
pub(crate) fn agc2(sig_in: &[i16], sig_out: &mut [i16], l_trm: usize, ov: &mut Flag) {
    let mut s = energy_new(sig_out, l_trm, ov);
    if s == 0 {
        return;
    }
    let mut exp = norm_l(s) - 1;
    let mut l_temp = l_shl(s, exp, ov);
    let gain_out = round(l_temp, ov);

    s = energy_new(sig_in, l_trm, ov);
    let g0 = if s == 0 {
        0
    } else {
        let i = norm_l(s);
        l_temp = l_shl(s, i, ov);
        let gain_in = round(l_temp, ov);
        exp -= i;
        let temp = div_s(gain_out, gain_in);
        s = temp as i32;
        s = if s > 0x00FF_FFFF {
            MAX_32
        } else if s < -16_777_216 {
            MIN_32
        } else {
            s << 7
        };
        s = l_shr(s, exp, ov);
        s = inv_sqrt(s, ov);
        l_temp = if s > 0x003F_FFFF {
            MAX_32
        } else if s < -4_194_304 {
            MIN_32
        } else {
            s << 9
        };
        round(l_temp, ov)
    };

    for i in (0..l_trm).rev() {
        let l_temp = l_mult(sig_out[i], g0, ov);
        sig_out[i] = if l_temp > 0x0FFF_FFFF {
            MAX_16
        } else if l_temp < -268_435_456 {
            MIN_16
        } else {
            (l_temp >> 13) as i16
        };
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::post::test_gen::signal;
    use crate::testutil::Rng;
    use opencore_sys as _;

    /// `#[repr(C)]` mirror of the C state struct (one `Word16`).
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct CState(i16);

    unsafe extern "C" {
        // agc, agc2, agc_reset against the TS 26.073 reference (values and state: its flag is
        // global); the energy helpers (static there) against opencore.
        #[link_name = "ref_agc"]
        fn c_agc(st: *mut CState, i: *mut i16, o: *mut i16, fac: i16, l: i16) -> i32;
        #[link_name = "ref_agc2"]
        fn c_agc2(i: *mut i16, o: *mut i16, l: i16);
        #[link_name = "ref_agc_reset"]
        fn agc_reset(st: *mut CState) -> i16;
        fn energy_old_Wrapper(i: *mut i16, l: i16, ov: *mut i32) -> i32;
        fn energy_new_Wrapper(i: *mut i16, l: i16, ov: *mut i32) -> i32;
    }

    /// Signal pairs: random, all-zero (the `s == 0` branches) and full-scale (energy saturation).
    fn sig(g: &mut Rng, buf: &mut [i16]) {
        match g.u64() % 10 {
            0 => buf.fill(0),
            1 => buf
                .iter_mut()
                .for_each(|v| *v = [i16::MIN, i16::MAX][(g.u64() & 1) as usize]),
            _ => signal(g, buf),
        }
    }

    const LENS: [usize; 5] = [40, 40, 160, 1, 13];

    #[test]
    fn post_energy_vs_ref() {
        let mut g = Rng::new(20);
        for case in 0..50_000 {
            let l = LENS[case % LENS.len()];
            let mut x = [0i16; 160];
            sig(&mut g, &mut x);
            for new in [false, true] {
                let mut ov = g.u64() & 1 == 0;
                let mut cov = ov as i32;
                let (r, c) = if new {
                    (energy_new(&x, l, &mut ov), unsafe {
                        energy_new_Wrapper(x.as_mut_ptr(), l as i16, &mut cov)
                    })
                } else {
                    (energy_old(&x, l, &mut ov), unsafe {
                        energy_old_Wrapper(x.as_mut_ptr(), l as i16, &mut cov)
                    })
                };
                assert_eq!((r, ov), (c, cov != 0), "case {case} new={new}");
            }
        }
    }

    #[test]
    fn post_agc_vs_ref() {
        let mut g = Rng::new(21);
        for seq in 0..20 {
            let mut st = AgcState { past_gain: g.i16() };
            let mut cst = CState(st.past_gain);
            if seq % 2 == 0 {
                st.reset();
                assert_eq!(unsafe { agc_reset(&mut cst) }, 0);
                assert_eq!(st.past_gain, cst.0);
            }
            for call in 0..1000 {
                let l = LENS[call % LENS.len()];
                let fac = if seq % 3 == 0 { g.i16() } else { 29491 };
                let (mut i, mut o) = ([0i16; 160], [0i16; 160]);
                sig(&mut g, &mut i);
                sig(&mut g, &mut o);
                let (mut ci, mut co) = (i, o);
                let mut ov = false;
                agc(&mut st, &i, &mut o, fac, l, &mut ov);
                {
                    let _l = ref26073_sys::lock();
                    unsafe { c_agc(&mut cst, ci.as_mut_ptr(), co.as_mut_ptr(), fac, l as i16) };
                }
                assert_eq!((o, st.past_gain), (co, cst.0), "seq {seq} call {call}");
                assert_eq!(i, ci);
            }
        }
    }

    #[test]
    fn post_agc2_vs_ref() {
        let mut g = Rng::new(22);
        for case in 0..50_000 {
            let l = LENS[case % LENS.len()];
            let (mut i, mut o) = ([0i16; 160], [0i16; 160]);
            sig(&mut g, &mut i);
            sig(&mut g, &mut o);
            let (mut ci, mut co) = (i, o);
            let mut ov = false;
            agc2(&i, &mut o, l, &mut ov);
            {
                let _l = ref26073_sys::lock();
                unsafe { c_agc2(ci.as_mut_ptr(), co.as_mut_ptr(), l as i16) };
            }
            assert_eq!(o, co, "case {case}");
        }
    }
}
