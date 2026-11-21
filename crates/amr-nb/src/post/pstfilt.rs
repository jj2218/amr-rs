//! `Post_Filter`: adaptive formant post-filter (opencore `dec/src/pstfilt.cpp`).

use super::agc::{AgcState, agc};
use super::filt::{residu, syn_filt, weight_ai};
use super::preemph::{PreemphasisState, preemphasis};
use super::{AGC_FAC, CodecMode, L_FRAME, L_SUBFR, M, MP1, MU};
use crate::basic_op::{Flag, MAX_32, div_s, l_mac, l_mult};

/// `L_H`: length of the truncated impulse response.
const L_H: usize = 22;

const GAMMA3_MR122: [i16; M] = [22938, 16057, 11240, 7868, 5508, 3856, 2699, 1889, 1322, 925];
const GAMMA3: [i16; M] = [18022, 9912, 5451, 2998, 1649, 907, 499, 274, 151, 83];
const GAMMA4_MR122: [i16; M] = [
    24576, 18432, 13824, 10368, 7776, 5832, 4374, 3281, 2461, 1846,
];
const GAMMA4: [i16; M] = [22938, 16057, 11240, 7868, 5508, 3856, 2699, 1889, 1322, 925];

/// `Post_FilterState`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct PostFilterState {
    pub(crate) res2: [i16; L_SUBFR],
    pub(crate) mem_syn_pst: [i16; M],
    pub(crate) preemph_state: PreemphasisState,
    pub(crate) agc_state: AgcState,
    pub(crate) synth_buf: [i16; M + L_FRAME],
}

impl Default for PostFilterState {
    fn default() -> Self {
        Self {
            res2: [0; L_SUBFR],
            mem_syn_pst: [0; M],
            preemph_state: PreemphasisState::default(),
            agc_state: AgcState::default(),
            synth_buf: [0; M + L_FRAME],
        }
    }
}

impl PostFilterState {
    /// `Post_Filter_reset`.
    pub(crate) fn reset(&mut self) {
        self.mem_syn_pst = [0; M];
        self.res2 = [0; L_SUBFR];
        self.synth_buf = [0; M + L_FRAME];
        self.agc_state.reset();
        self.preemph_state.reset();
    }
}

/// C's `L_mult`-like accumulation in the tilt computation: on the single overflowing product
/// it sets the flag and *stops* (the saturated term is not added).
fn autocorr_h(h: &[i16; L_H], lag: usize, ov: &mut Flag) -> i32 {
    // Fast path: if 2 sum|h[i] h[i + lag]| <= MAX_32, no partial sum of the L_add chain
    // saturates and no product is 0x4000_0000 (the C's break), so the plain sum is exact.
    let (sum, abs) = h[..L_H - lag]
        .iter()
        .zip(&h[lag..])
        .fold((0i32, 0i64), |(s, a), (&x, &y)| {
            let p = x as i32 * y as i32;
            (s.wrapping_add(p << 1), a + (p.unsigned_abs() as i64))
        });
    if 2 * abs <= MAX_32 as i64 {
        return sum;
    }
    autocorr_h_chain(h, lag, ov)
}

/// `autocorr_h`'s chain as TS 26.073: `L_mult(h[0], h[lag])`, then `L_mac` upwards, saturating.
/// (opencore adds from the top and stops at a -32768^2 product: different once it saturates.)
#[cold]
fn autocorr_h_chain(h: &[i16; L_H], lag: usize, ov: &mut Flag) -> i32 {
    let mut l_tmp = l_mult(h[0], h[lag], ov);
    for i in 1..L_H - lag {
        l_tmp = l_mac(l_tmp, h[i], h[i + lag], ov);
    }
    l_tmp
}

/// `Post_Filter`: filters one frame `syn[..L_FRAME]` in place, using the four subframe LP
/// filters `az_4[..4 * MP1]`.
pub(crate) fn post_filter(
    st: &mut PostFilterState,
    mode: CodecMode,
    syn: &mut [i16],
    az_4: &[i16],
    ov: &mut Flag,
) {
    let mut ap3 = [0i16; MP1];
    let mut ap4 = [0i16; MP1];
    // syn_work = &synth_buf[M]
    st.synth_buf[M..M + L_FRAME].copy_from_slice(&syn[..L_FRAME]);

    for (sf, i_subfr) in (0..L_FRAME).step_by(L_SUBFR).enumerate() {
        let az = &az_4[sf * MP1..(sf + 1) * MP1];
        if mode == CodecMode::MR122 || mode == CodecMode::MR102 {
            weight_ai(az, &GAMMA3_MR122, &mut ap3);
            weight_ai(az, &GAMMA4_MR122, &mut ap4);
        } else {
            weight_ai(az, &GAMMA3, &mut ap3);
            weight_ai(az, &GAMMA4, &mut ap4);
        }

        // Residu(Ap3, &syn_work[i_subfr], res2, L_SUBFR): input starts M samples earlier.
        residu(&ap3, &st.synth_buf[i_subfr..], &mut st.res2, L_SUBFR);

        // Impulse response of Ap3/Ap4: C filters `h` in place with `mem = &h[M + 1]` (zeros,
        // update = 0); equivalent to filtering a copy (see `syn_filt`).
        // Filtered over 24 samples (x zero past 22, the filter causal: the first 22 outputs
        // are the C's): whole blocks of 8, no tail block (its history rebuild cost more).
        let mut x = [0i16; L_H.next_multiple_of(8)];
        x[..MP1].copy_from_slice(&ap3);
        let mut mem = [0i16; M]; // the C's &h[M + 1..]: zeros
        let mut hp = [0i16; L_H.next_multiple_of(8)];
        syn_filt(&ap4, &x, &mut hp, x.len(), &mut mem, false);
        let h: [i16; L_H] = hp[..L_H].try_into().unwrap();

        let temp1 = (autocorr_h(&h, 0, ov) >> 16) as i16;
        let mut temp2 = (autocorr_h(&h, 1, ov) >> 16) as i16;
        if temp2 <= 0 {
            temp2 = 0;
        } else {
            let mut l_tmp = (temp2 as i32 * MU as i32) >> 15;
            if l_tmp & 0x0001_0000 != 0 {
                l_tmp |= 0xffff_0000u32 as i32;
            }
            temp2 = l_tmp as i16;
            temp2 = div_s(temp2, temp1);
        }

        preemphasis(&mut st.preemph_state, &mut st.res2, temp2, L_SUBFR, ov);
        syn_filt(
            &ap4,
            &st.res2,
            &mut syn[i_subfr..i_subfr + L_SUBFR],
            L_SUBFR,
            &mut st.mem_syn_pst,
            true,
        );
        agc(
            &mut st.agc_state,
            &st.synth_buf[M + i_subfr..M + i_subfr + L_SUBFR],
            &mut syn[i_subfr..i_subfr + L_SUBFR],
            AGC_FAC,
            L_SUBFR,
            ov,
        );
    }
    // memmove(&syn_work[-M], &syn_work[L_FRAME - M], M)
    st.synth_buf.copy_within(L_FRAME..L_FRAME + M, 0);
}

#[cfg(all(test, oracle))]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::post::test_gen::{coefs, signal};
    use crate::testutil::{N, Rng};
    use opencore_sys as _;

    /// `#[repr(C)]` mirror of `Post_FilterState` (all `Word16`, no padding).
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct CState {
        res2: [i16; L_SUBFR],
        mem_syn_pst: [i16; M],
        mem_pre: i16,
        past_gain: i16,
        synth_buf: [i16; M + L_FRAME],
    }

    impl CState {
        fn of(s: &PostFilterState) -> Self {
            Self {
                res2: s.res2,
                mem_syn_pst: s.mem_syn_pst,
                mem_pre: s.preemph_state.mem_pre,
                past_gain: s.agc_state.past_gain,
                synth_buf: s.synth_buf,
            }
        }
    }

    /// The TS 26.073 reference's `Post_FilterState`: its preemphasis and agc states (one
    /// `Word16` each) behind pointers.
    #[repr(C)]
    struct RState {
        res2: [i16; L_SUBFR],
        mem_syn_pst: [i16; M],
        preemph: *mut i16,
        agc: *mut i16,
        synth_buf: [i16; M + L_FRAME],
    }

    /// `CState` as a reference state (`subs` holds mem_pre and past_gain).
    fn to_ref(c: &CState, subs: &mut [i16; 2]) -> RState {
        *subs = [c.mem_pre, c.past_gain];
        let p = subs.as_mut_ptr();
        RState {
            res2: c.res2,
            mem_syn_pst: c.mem_syn_pst,
            preemph: p,
            agc: p.wrapping_add(1),
            synth_buf: c.synth_buf,
        }
    }

    fn from_ref(r: &RState, subs: &[i16; 2]) -> CState {
        CState {
            res2: r.res2,
            mem_syn_pst: r.mem_syn_pst,
            mem_pre: subs[0],
            past_gain: subs[1],
            synth_buf: r.synth_buf,
        }
    }

    unsafe extern "C" {
        fn ref_Post_Filter(st: *mut RState, mode: i32, syn: *mut i16, az: *mut i16) -> i32;
        fn ref_Post_Filter_reset(st: *mut RState) -> i32;
    }

    #[allow(non_snake_case)]
    fn Post_Filter_reset(c: &mut CState) -> i32 {
        let _l = ref26073_sys::lock();
        let mut subs = [0i16; 2];
        let mut r = to_ref(c, &mut subs);
        let ret = unsafe { ref_Post_Filter_reset(&mut r) };
        *c = from_ref(&r, &subs);
        ret
    }

    const MODES: [CodecMode; 9] = [
        CodecMode::MR475,
        CodecMode::MR515,
        CodecMode::MR59,
        CodecMode::MR67,
        CodecMode::MR74,
        CodecMode::MR795,
        CodecMode::MR102,
        CodecMode::MR122,
        CodecMode::MRDTX,
    ];

    /// `autocorr_h`'s bounded plain sum vs the C's loop (`autocorr_h_chain`), both lags.
    /// Classes, counted: the fast path, a declined bound where the chain still doesn't
    /// saturate, a saturating chain, and the break on a (-32768)^2 product.
    #[test]
    fn post_autocorr_h_vs_chain() {
        let mut g = Rng::new(31);
        let (mut fast, mut declined, mut sat, mut brk) = (0, 0, 0, 0);
        for n in 0..N {
            let shift = [0u32, 1, 2, 4, 8][n % 5];
            let mut h: [i16; L_H] = core::array::from_fn(|_| g.i16() >> shift);
            if n % 7 == 0 {
                // same sign everywhere: the lag-1 sum grows monotonically
                h.iter_mut()
                    .for_each(|v| *v = v.saturating_abs() >> (n % 3));
            }
            if n % 11 == 0 {
                let i = (g.u64() % (L_H as u64 - 1)) as usize;
                (h[i], h[i + 1]) = (i16::MIN, i16::MIN);
            }
            for lag in [0, 1] {
                let (mut o, mut oc) = (false, false);
                let r = autocorr_h(&h, lag, &mut o);
                let c = autocorr_h_chain(&h, lag, &mut oc);
                assert_eq!((r, o), (c, oc), "lag {lag} {h:?}");
                let abs: i64 = (0..L_H - lag)
                    .map(|i| (h[i] as i32 * h[i + lag] as i32).unsigned_abs() as i64)
                    .sum();
                let has_min = (0..L_H - lag).any(|i| h[i] == i16::MIN && h[i + lag] == i16::MIN);
                match (2 * abs <= MAX_32 as i64, oc) {
                    (true, _) => fast += 1,
                    (false, false) => declined += 1,
                    (false, true) if has_min => brk += 1,
                    (false, true) => sat += 1,
                }
            }
        }
        assert!(
            fast > N / 4 && declined > N / 100 && sat > N / 100 && brk > N / 100,
            "classes: fast {fast}, declined {declined}, sat {sat}, break {brk}"
        );
    }

    #[test]
    fn post_filter_vs_ref() {
        let mut g = Rng::new(30);
        for seq in 0..12 {
            // Start from reset, or from a random (dirty) state.
            let mut st = PostFilterState::default();
            st.res2.iter_mut().for_each(|v| *v = g.i16());
            st.mem_syn_pst.iter_mut().for_each(|v| *v = g.i16());
            st.synth_buf.iter_mut().for_each(|v| *v = g.i16());
            st.preemph_state.mem_pre = g.i16();
            st.agc_state.past_gain = g.i16();
            let mut cst = CState::of(&st);
            if seq % 2 == 0 {
                st.reset();
                assert_eq!(Post_Filter_reset(&mut cst), 0);
                assert_eq!(CState::of(&st), cst);
            }
            for frame in 0..1000 {
                // Fixed mode per sequence for most, switching every frame for the rest.
                let mode = if seq < 9 {
                    MODES[seq % 9]
                } else {
                    MODES[(g.u64() % 9) as usize]
                };
                let mut syn = [0i16; L_FRAME];
                match g.u64() % 12 {
                    0 => {}
                    _ => signal(&mut g, &mut syn),
                }
                let mut az = [0i16; 4 * MP1];
                for sf in 0..4 {
                    az[sf * MP1..(sf + 1) * MP1].copy_from_slice(&coefs(&mut g));
                }
                let mut csyn = syn;
                let mut ov = false;
                post_filter(&mut st, mode, &mut syn, &az, &mut ov);
                {
                    let _l = ref26073_sys::lock();
                    let mut subs = [0i16; 2];
                    let mut r = to_ref(&cst, &mut subs);
                    unsafe {
                        ref_Post_Filter(&mut r, mode as i32, csyn.as_mut_ptr(), az.as_mut_ptr())
                    };
                    cst = from_ref(&r, &subs);
                }
                assert_eq!(syn, csyn, "seq {seq} frame {frame} mode {mode:?}: output");
                assert_eq!(CState::of(&st), cst, "seq {seq} frame {frame}: state");
            }
        }
    }
}
