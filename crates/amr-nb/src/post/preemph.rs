//! `preemphasis` (opencore `dec/src/preemph.cpp`).

use super::L_FRAME;
use crate::basic_op::{Flag, MAX_16, MIN_16};

/// `preemphasisState`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct PreemphasisState {
    pub(crate) mem_pre: i16,
}

impl PreemphasisState {
    /// `preemphasis_reset`.
    pub(crate) fn reset(&mut self) {
        self.mem_pre = 0;
    }
}

/// `preemphasis`: in-place filter `signal[i] -= g * signal[i - 1]` over `l` samples, with
/// `signal[-1]` taken from (and the last input sample saved to) `st.mem_pre`; `l` in 1..=L_FRAME.
pub(crate) fn preemphasis(
    st: &mut PreemphasisState,
    signal: &mut [i16],
    g: i16,
    l: usize,
    ov: &mut Flag,
) {
    // The C runs backwards, so every output reads the unmodified x[i - 1]: independent
    // lanes, `mult`'s and `sub`'s saturation as selects, their overflow in one flag.
    let x = &mut signal[..l];
    debug_assert!((1..=L_FRAME).contains(&l));
    let mut prev = [0i16; L_FRAME];
    prev[0] = st.mem_pre;
    prev[1..l].copy_from_slice(&x[..l - 1]);
    st.mem_pre = x[l - 1];
    let mut sat = 0u32;
    for (o, &p) in x.iter_mut().zip(&prev[..l]) {
        let m0 = (g as i32 * p as i32) >> 15; // mult: > MAX_16 only for (-32768)^2
        let m = m0.min(0x7fff);
        let d = *o as i32 - m;
        let r = d.clamp(MIN_16 as i32, MAX_16 as i32);
        sat |= (m0 != m) as u32 | (r != d) as u32;
        *o = r as i16;
    }
    *ov |= sat != 0;
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_c {
    use super::*;
    use crate::basic_op::{mult, sub};
    use crate::post::test_gen::signal;
    use crate::testutil::Rng;
    use opencore_sys as _;

    /// `#[repr(C)]` mirror of the C state struct (one `Word16`).
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct CState(i16);

    // The TS 26.073 reference, in the shape of the opencore functions these tests were
    // written against (`ov`/`o`: its overflow flag; table arguments unused).
    unsafe extern "C" {
        fn ref_preemphasis(st: *mut CState, sig: *mut i16, g: i16, l: i16) -> i32;
        fn ref_preemphasis_reset(st: *mut CState) -> i32;
    }
    unsafe fn c_preemphasis(st: *mut CState, sig: *mut i16, g: i16, l: i16, ov: *mut i32) {
        crate::testutil::with_ref(ov, || unsafe { ref_preemphasis(st, sig, g, l) });
    }
    unsafe fn preemphasis_reset(st: *mut CState) -> i16 {
        crate::testutil::with_ref(std::ptr::null_mut(), || unsafe {
            ref_preemphasis_reset(st)
        }) as i16
    }

    #[test]
    fn post_preemphasis_vs_c() {
        let mut g = Rng::new(10);
        let (mut mult_sat, mut sub_sat) = (0u32, 0u32);
        for seq in 0..20 {
            let mut st = PreemphasisState { mem_pre: g.i16() };
            let mut cst = CState(st.mem_pre);
            if seq % 2 == 0 {
                st.reset();
                assert_eq!(unsafe { preemphasis_reset(&mut cst) }, 0);
                assert_eq!(st.mem_pre, cst.0);
            }
            let gain = g.i16();
            for call in 0..1000 {
                let l = [40usize, 40, 160, 1, 2, 7][call % 6];
                let mut sig = [0i16; 160];
                signal(&mut g, &mut sig);
                let mut gg = if call % 4 == 0 { g.i16() } else { gain };
                if call % 5 == 0 {
                    // mult's only saturation, (-32768)^2, then sub's (x - 32767 < -32768)
                    gg = i16::MIN;
                    for v in &mut sig[..l] {
                        if g.u64().is_multiple_of(3) {
                            *v = i16::MIN;
                        }
                    }
                }
                // observed saturations of each operator (the scalar reference ops)
                for i in 0..l {
                    let (mut om, mut os) = (false, false);
                    let prev = if i == 0 { st.mem_pre } else { sig[i - 1] };
                    let p = mult(gg, prev, &mut om);
                    sub(sig[i], p, &mut os);
                    mult_sat += om as u32;
                    sub_sat += os as u32;
                }
                let mut csig = sig;
                let mut ov = g.u64().is_multiple_of(8);
                let mut cov = ov as i32;
                preemphasis(&mut st, &mut sig, gg, l, &mut ov);
                unsafe { c_preemphasis(&mut cst, csig.as_mut_ptr(), gg, l as i16, &mut cov) };
                let _ = (ov, cov); // the reference's flag is global: values and state only
                assert_eq!((sig, st.mem_pre), (csig, cst.0), "seq {seq} call {call}");
            }
        }
        assert!(
            mult_sat > 100 && sub_sat > 100,
            "classes: {mult_sat} {sub_sat}"
        );
    }
}
