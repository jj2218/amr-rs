//! `Post_Process`: output high-pass filter and up-scaling (opencore `dec/src/post_pro.cpp`).

use crate::basic_op::{Flag, l_shl, round};

const B: [i16; 3] = [7699, -15398, 7699];
const A: [i16; 3] = [8192, 15836, -7667];

/// `Post_ProcessState`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct PostProcessState {
    pub(crate) y2_hi: i16,
    pub(crate) y2_lo: i16,
    pub(crate) y1_hi: i16,
    pub(crate) y1_lo: i16,
    pub(crate) x0: i16,
    pub(crate) x1: i16,
}

impl PostProcessState {
    /// `Post_Process_reset`.
    pub(crate) fn reset(&mut self) {
        *self = Self::default();
    }
}

/// The (hi, lo) state pairs y1, y2 after `n` samples of a one-word IIR (`post_process`,
/// `pre_process`): the pairs the filter computed decode from u as (u >> 15, u & 0x7fff); pairs
/// it didn't compute (n < 2) are the input pairs, kept as they are (as in the C, even with a
/// lo outside [0, 0x7fff]).
pub(crate) fn one_word_pairs(
    n: usize,
    u1: i32,
    u2: i32,
    y1: (i16, i16),
    y2: (i16, i16),
) -> ((i16, i16), (i16, i16)) {
    let pair = |u: i32| ((u >> 15) as i16, (u & 0x7fff) as i16);
    match n {
        0 => (y1, y2),
        1 => (pair(u1), y1),
        _ => (pair(u1), pair(u2)),
    }
}

/// `Post_Process`: 2nd-order high-pass IIR (cut-off 60 Hz) with x2 gain, in place on
/// `signal[..lg]`.
pub(crate) fn post_process(
    st: &mut PostProcessState,
    signal: &mut [i16],
    lg: usize,
    ov: &mut Flag,
) {
    if !post_process_fast(st, &mut signal[..lg], ov) {
        post_process_exact(st, signal, lg, ov);
    }
}

/// Longest `lg` the fast path takes (the decoder passes `L_FRAME`).
const MAX_LG: usize = 160;

/// `post_process` while every y = (x and feedback terms) stays in [-2^28, 2^28): then the C's
/// `l_shl(y, 3)` can't saturate, so its state u = l_tmp >> 1 is exactly 4y and the serial
/// chain is y -> imul -> sar -> add (the C's chain adds the saturation test and the shifts).
/// The outputs `round(l_shl(8y, 1))` don't feed back: a separate pass of independent lanes.
/// Returns false, with `st`, `signal` and `ov` untouched, if some y leaves the range.
fn post_process_fast(st: &mut PostProcessState, sig: &mut [i16], ov: &mut Flag) -> bool {
    const LIM: i32 = 1 << 28;
    let n = sig.len();
    if n > MAX_LG {
        return false;
    }
    let (c_a1, c_a2) = (A[1] as i64, A[2] as i64);
    let (c_b0, c_b1, c_b2) = (B[0] as i32, B[1] as i32, B[2] as i32);
    let mut u1 = ((st.y1_hi as i32) << 15) + st.y1_lo as i32;
    let mut u2 = ((st.y2_hi as i32) << 15) + st.y2_lo as i32;
    let (mut x1, mut x0) = (st.x1, st.x0);
    let mut ys = [0i32; MAX_LG];
    let mut bad = false;
    for (&v, yo) in sig.iter().zip(&mut ys) {
        let x2 = x1;
        x1 = x0;
        x0 = v;
        let x = (x0 as i32 * c_b0)
            .wrapping_add(x1 as i32 * c_b1)
            .wrapping_add(x2 as i32 * c_b2);
        // The u2 term and x are known a sample early: add the u1 term last.
        let rest = (((u2 as i64 * c_a2) >> 15) as i32).wrapping_add(x);
        let y = rest.wrapping_add(((u1 as i64 * c_a1) >> 15) as i32);
        bad |= !(-LIM..LIM).contains(&y);
        *yo = y;
        u2 = u1;
        u1 = y.wrapping_shl(2);
    }
    if bad {
        return false;
    }
    // round(l_shl(8y, 1)): 16y saturates outside [-2^27, 2^27); round's l_add overflows (and
    // sets ov) for t > MAX_32 - 0x8000, a saturated t or 16y just below 2^31.
    let mut o = false;
    for (s, &y) in sig.iter_mut().zip(&ys) {
        let t = if y >= 1 << 27 {
            i32::MAX
        } else if y < -(1 << 27) {
            i32::MIN
        } else {
            y << 4
        };
        let r = t > i32::MAX - 0x8000;
        o |= r;
        *s = if r {
            i16::MAX
        } else {
            ((t + 0x8000) >> 16) as i16
        };
    }
    *ov |= o;
    (st.x1, st.x0) = (x1, x0);
    ((st.y1_hi, st.y1_lo), (st.y2_hi, st.y2_lo)) =
        one_word_pairs(n, u1, u2, (st.y1_hi, st.y1_lo), (st.y2_hi, st.y2_lo));
    true
}

/// `post_process`'s exact loop (any y: `l_shl`'s saturation on the chain).
fn post_process_exact(st: &mut PostProcessState, signal: &mut [i16], lg: usize, ov: &mut Flag) {
    // The C keeps y1, y2 as (hi, lo) pairs with u = t >> 1 = (hi << 15) + lo, lo in
    // [0, 0x7fff]; hi*a + ((lo*a) >> 15) is exactly (u*a) >> 15 (hi*a is a multiple of 2^15
    // in u*a), and as an i64 product it equals the C's wrapping i32 sum mod 2^32. One
    // multiply per tap on the serial chain instead of the split and two.
    let (c_a1, c_a2) = (A[1] as i64, A[2] as i64);
    let (c_b0, c_b1, c_b2) = (B[0] as i32, B[1] as i32, B[2] as i32);
    let mut u1 = ((st.y1_hi as i32) << 15) + st.y1_lo as i32;
    let mut u2 = ((st.y2_hi as i32) << 15) + st.y2_lo as i32;
    let (mut x1, mut x0) = (st.x1, st.x0);
    for s in &mut signal[..lg] {
        let x2 = x1;
        x1 = x0;
        x0 = *s;
        let x = (x0 as i32 * c_b0)
            .wrapping_add(x1 as i32 * c_b1)
            .wrapping_add(x2 as i32 * c_b2);
        let y = (((u1 as i64 * c_a1) >> 15) as i32).wrapping_add(((u2 as i64 * c_a2) >> 15) as i32);
        let l_tmp = l_shl(y.wrapping_add(x), 3, ov);
        *s = round(l_shl(l_tmp, 1, ov), ov);
        u2 = u1;
        u1 = l_tmp >> 1;
    }
    (st.x1, st.x0) = (x1, x0);
    ((st.y1_hi, st.y1_lo), (st.y2_hi, st.y2_lo)) =
        one_word_pairs(lg, u1, u2, (st.y1_hi, st.y1_lo), (st.y2_hi, st.y2_lo));
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_c {
    use super::*;
    use crate::post::test_gen::signal;
    use crate::testutil::Rng;
    use opencore_sys as _;

    /// `#[repr(C)]` mirror of `Post_ProcessState`.
    #[repr(C)]
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct CState([i16; 6]);

    impl CState {
        fn of(s: &PostProcessState) -> Self {
            Self([s.y2_hi, s.y2_lo, s.y1_hi, s.y1_lo, s.x0, s.x1])
        }
    }

    // The TS 26.073 reference, in the shape of the opencore functions these tests were
    // written against (`ov`/`o`: its overflow flag; table arguments unused).
    unsafe extern "C" {
        fn ref_Post_Process(st: *mut CState, sig: *mut i16, lg: i16) -> i32;
        fn ref_Post_Process_reset(st: *mut CState) -> i32;
    }
    #[allow(non_snake_case)]
    unsafe fn Post_Process(st: *mut CState, sig: *mut i16, lg: i16, ov: *mut i32) {
        crate::testutil::with_ref(ov, || unsafe { ref_Post_Process(st, sig, lg) });
    }
    #[allow(non_snake_case)]
    unsafe fn Post_Process_reset(st: *mut CState) -> i16 {
        crate::testutil::with_ref(std::ptr::null_mut(), || unsafe {
            ref_Post_Process_reset(st)
        }) as i16
    }

    /// Random and reset states (any lo), random signals, plus full-scale square waves of
    /// random amplitude. Classes counted: the fast path, its fallback (some y outside
    /// [-2^28, 2^28)), and fast calls with a full-scale output (+ and -).
    #[test]
    fn post_process_vs_c() {
        let mut g = Rng::new(40);
        let (mut fast, mut exact, mut sat_hi, mut sat_lo) = (0, 0, 0, 0);
        for seq in 0..20 {
            let mut st = PostProcessState {
                y2_hi: g.i16(),
                y2_lo: g.i16(),
                y1_hi: g.i16(),
                y1_lo: g.i16(),
                x0: g.i16(),
                x1: g.i16(),
            };
            let mut cst = CState::of(&st);
            if seq % 2 == 0 {
                st.reset();
                assert_eq!(unsafe { Post_Process_reset(&mut cst) }, 0);
                assert_eq!(CState::of(&st), cst);
            }
            for call in 0..1000 {
                // a length-1 first call on a random state: one input pair is kept as is
                let lg = if call == 0 && seq % 4 == 1 {
                    1
                } else {
                    [160usize, 160, 40, 1, 7][call % 5]
                };
                let mut sig = [0i16; 160];
                if call % 3 == 2 {
                    let (amp, half) = ((g.u64() % 32768) as i16, 1 + (g.u64() % 80) as usize);
                    for (i, v) in sig.iter_mut().enumerate() {
                        *v = if (i / half) % 2 == 0 { amp } else { -amp };
                    }
                } else {
                    signal(&mut g, &mut sig);
                }
                let mut csig = sig;
                let (mut st_f, mut sig_f, mut ov_f) = (st, sig, false);
                if post_process_fast(&mut st_f, &mut sig_f[..lg], &mut ov_f) {
                    fast += 1;
                    sat_hi += sig_f[..lg].contains(&i16::MAX) as u32;
                    sat_lo += sig_f[..lg].contains(&i16::MIN) as u32;
                } else {
                    exact += 1;
                }
                let mut ov = g.u64().is_multiple_of(8);
                let mut cov = ov as i32;
                post_process(&mut st, &mut sig, lg, &mut ov);
                unsafe { Post_Process(&mut cst, csig.as_mut_ptr(), lg as i16, &mut cov) };
                let _ = (ov, cov); // the reference's flag is global: values and state only
                assert_eq!((sig, CState::of(&st)), (csig, cst), "seq {seq} call {call}");
            }
        }
        assert!(
            fast > 10_000 && exact > 2000 && sat_hi > 2000 && sat_lo > 2000,
            "classes: fast {fast}, exact {exact}, sat_hi {sat_hi}, sat_lo {sat_lo}"
        );
    }

    /// The fast path's edges, one sample from a state that puts y exactly there (x = 0,
    /// u2 = 0, u1 solved for y = (u1 * a1) >> 15): round's overflow below 2^27 (16y + 0x8000
    /// past MAX_32 without 16y saturating), 16y's saturation at +-2^27, the fallback at +-2^28.
    #[test]
    fn post_process_edges_vs_c() {
        let mut n = 0;
        for edge in [(1 << 27) - 0x800, 1 << 27, -(1 << 27), 1 << 28, -(1 << 28)] {
            for d in -40i64..40 {
                let y = edge as i64 + d;
                // u0 = floor(y 2^15 / a1): (u0 * a1) >> 15 is y or y - 1 (a1 > 0); u0 + 1 reaches
                // y and beyond, so the three cover the edge
                let u0 = (y << 15).div_euclid(A[1] as i64);
                for u1 in u0 - 1..=u0 + 1 {
                    let st = PostProcessState {
                        y1_hi: (u1 >> 15) as i16,
                        y1_lo: (u1 & 0x7fff) as i16,
                        ..Default::default()
                    };
                    let cst = CState::of(&st);
                    for lg in [1usize, 3] {
                        let (mut sig, mut csig) = ([0i16; 3], [0i16; 3]);
                        let (mut ov, mut cov) = (false, 0);
                        let (mut s, mut c) = (st, cst);
                        post_process(&mut s, &mut sig, lg, &mut ov);
                        unsafe { Post_Process(&mut c, csig.as_mut_ptr(), lg as i16, &mut cov) };
                        let _ = (ov, cov);
                        assert_eq!((sig, CState::of(&s)), (csig, c), "u1 {u1}");
                        n += 1;
                    }
                }
            }
        }
        assert_eq!(n, 5 * 80 * 3 * 2);
    }
}
