//! Adaptive phase dispersion / anti-sparseness processing (`dec/ph_disp.cpp`).

use super::{L_SUBFR, Mode};
use crate::basic_op::{Flag, MAX_32, MIN_32, add, l_add, l_mult, l_shl, round};
use crate::tables::{PH_IMP_LOW, PH_IMP_LOW_MR795, PH_IMP_MID, PH_IMP_MID_MR795};
use Mode::{MR74, MR102, MR122, MR795};

const PHDGAINMEMSIZE: usize = 5;
const PHDTHR1LTP: i16 = 9830;
const PHDTHR2LTP: i16 = 14746;
const ONFACTPLUS1: i32 = 16384;
const ONLENGTH: i16 = 2;

/// `ph_dispState`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(crate) struct PhDispState {
    pub(crate) gain_mem: [i16; PHDGAINMEMSIZE],
    pub(crate) prev_state: i16,
    pub(crate) prev_cb_gain: i16,
    pub(crate) lock_full: i16,
    pub(crate) onset: i16,
}

impl PhDispState {
    pub(crate) fn reset(&mut self) {
        ph_disp_reset(self);
    }
}

/// `ph_disp_reset` (always returns 0: the C NULL check has no Rust equivalent).
pub(crate) fn ph_disp_reset(state: &mut PhDispState) -> i16 {
    *state = PhDispState::default();
    0
}

/// `ph_disp_lock`
pub(crate) fn ph_disp_lock(state: &mut PhDispState) {
    state.lock_full = 1;
}

/// `ph_disp_release`
pub(crate) fn ph_disp_release(state: &mut PhDispState) {
    state.lock_full = 0;
}

/// `ph_disp`. `mode` is the C `enum Mode` value (MR475 = 0 .. MR122 = 7). The
/// `common_amr_tbls` argument is replaced by the `PH_IMP_*` tables.
#[allow(clippy::too_many_arguments)]
pub(crate) fn ph_disp(
    state: &mut PhDispState,
    mode: Mode,
    x: &mut [i16],
    cb_gain: i16,
    ltp_gain: i16,
    inno: &mut [i16],
    pitch_fac: i16,
    tmp_shift: i16,
    ov: &mut Flag,
) {
    state.gain_mem[4] = state.gain_mem[3];
    state.gain_mem[3] = state.gain_mem[2];
    state.gain_mem[2] = state.gain_mem[1];
    state.gain_mem[1] = state.gain_mem[0];
    state.gain_mem[0] = ltp_gain;

    let mut imp_nr: i16 = if ltp_gain < PHDTHR2LTP {
        if ltp_gain > PHDTHR1LTP { 1 } else { 0 }
    } else {
        2
    };

    let mut l_temp = (state.prev_cb_gain as i32 * ONFACTPLUS1) << 1;
    if l_temp > 0x1fff_ffff {
        *ov = true;
        l_temp = MAX_32;
    } else if l_temp < -536_870_912 {
        *ov = true;
        l_temp = MIN_32;
    } else {
        l_temp <<= 2;
    }
    let tmp1 = round(l_temp, ov);
    if cb_gain > tmp1 {
        state.onset = ONLENGTH;
    } else if state.onset > 0 {
        state.onset -= 1;
    }

    if state.onset == 0 {
        let i1 = state.gain_mem.iter().filter(|&&g| g < PHDTHR1LTP).count();
        if i1 > 2 {
            imp_nr = 0;
        }
    }
    if imp_nr as i32 > state.prev_state as i32 + 1 && state.onset == 0 {
        imp_nr -= 1;
    }
    if imp_nr < 2 && state.onset > 0 {
        imp_nr += 1;
    }
    if cb_gain < 10 {
        imp_nr = 2;
    }
    if state.lock_full == 1 {
        imp_nr = 0;
    }
    state.prev_state = imp_nr;
    state.prev_cb_gain = cb_gain;

    if mode != MR122 && mode != MR102 && mode != MR74 && imp_nr < 2 {
        let mut inno_sav = [0i16; L_SUBFR];
        let mut ps_poss = [0usize; L_SUBFR];
        let mut nze = 0usize;
        for i in 0..L_SUBFR {
            if inno[i] != 0 {
                ps_poss[nze] = i;
                nze += 1;
            }
            inno_sav[i] = inno[i];
            inno[i] = 0;
        }
        let imp2: &[i16; 2 * L_SUBFR] = if mode == MR795 {
            if imp_nr == 0 {
                &IMP2_LOW_MR795
            } else {
                &IMP2_MID_MR795
            }
        } else if imp_nr == 0 {
            &IMP2_LOW
        } else {
            &IMP2_MID
        };
        // Circular convolution: inno[j] += imp[(j - ppos) mod L], i.e. a window of imp ++ imp
        // starting at L - ppos: one fixed-length saturating pass per pulse (vectorised; the
        // C's head/tail split gave variable-length loops that stayed scalar).
        let inno: &mut [i16; L_SUBFR] = (&mut inno[..L_SUBFR]).try_into().unwrap();
        let mut o = false;
        for &ppos in &ps_poss[..nze] {
            let w: &[i16; L_SUBFR] = imp2[L_SUBFR - ppos..][..L_SUBFR].try_into().unwrap();
            o |= add_scaled(inno, w, inno_sav[ppos]);
        }
        *ov |= o;
    }

    let x: &mut [i16; L_SUBFR] = (&mut x[..L_SUBFR]).try_into().unwrap();
    total_exc::<false>(
        x,
        inno[..L_SUBFR].try_into().unwrap(),
        pitch_fac,
        cb_gain,
        tmp_shift,
    );
}

/// The total excitation, per sample (independent lanes, branchless so they vectorise):
/// `x = round(l_shl(L_mac(l_mult(x, pf), c, g), shift))` with TS 26.073's operators (both
/// callers, `ph_disp` and the decoder; `MAC` is kept for the call sites and no longer differs).
/// The decoder's `shift` is 1 or 2; any shift below 16 is exact (it goes through `l_shl`, whose
/// branch on it is loop-invariant; tested -40..=3). No overflow output: the decoder resets the
/// flag before `Syn_filt`, its only reader, as the reference does.
#[inline]
pub(crate) fn total_exc<const MAC: bool>(
    x: &mut [i16; L_SUBFR],
    c: &[i16; L_SUBFR],
    pf: i16,
    g: i16,
    shift: i16,
) {
    debug_assert!(shift < 16);
    if pf != i16::MIN && g != i16::MIN && (0..=2).contains(&shift) {
        // Fast form (every caller's case): no L_mult saturates, so |x pf|, |c g| < 2^30 and
        // t = x pf + c g is exact in i32. The chain is then round(sat(t << k)), k = 1 + shift:
        // clamped where (t << k) + 0x8000 leaves i32.
        let k = 1 + shift as u32;
        let (hi, lo) = ((MAX_32 - 0x8000) >> k, i32::MIN >> k);
        for (xi, &ci) in x.iter_mut().zip(c) {
            let t = *xi as i32 * pf as i32 + ci as i32 * g as i32;
            let r = (t.wrapping_shl(k).wrapping_add(0x8000) >> 16) as i16;
            *xi = if t > hi {
                i16::MAX
            } else if t < lo {
                i16::MIN
            } else {
                r
            };
        }
        return;
    }
    for (xi, &ci) in x.iter_mut().zip(c) {
        let p = *xi as i32 * pf as i32;
        let pm = p == 0x4000_0000;
        let a = if pm { MAX_32 } else { p.wrapping_shl(1) };
        // L_mac = L_add(a, L_mult(c, g)), saturating, as TS 26.073 (both call sites; opencore's
        // L_mac skips the saturation when a equals the product, and one site uses a wrapping
        // `L_add(a, c * g << 1)`).
        let r = ci as i32 * g as i32;
        let rm = r == 0x4000_0000;
        let b = if rm { MAX_32 } else { r << 1 };
        let sum = a.wrapping_add(b);
        let f = ((a ^ sum) & (b ^ sum)) < 0;
        let s = if f { (a >> 31) ^ MAX_32 } else { sum };
        let t = l_shl(s, shift, &mut false);
        let ro = t > MAX_32 - 0x8000;
        *xi = if ro {
            i16::MAX
        } else {
            ((t + 0x8000) >> 16) as i16
        };
    }
}

/// `imp ++ imp`, so a circular shift of `imp` is one contiguous window.
const fn doubled(imp: &[i16; L_SUBFR]) -> [i16; 2 * L_SUBFR] {
    let mut d = [0i16; 2 * L_SUBFR];
    let mut i = 0;
    while i < 2 * L_SUBFR {
        d[i] = imp[i % L_SUBFR];
        i += 1;
    }
    d
}
const IMP2_LOW: [i16; 2 * L_SUBFR] = doubled(&PH_IMP_LOW);
const IMP2_MID: [i16; 2 * L_SUBFR] = doubled(&PH_IMP_MID);
const IMP2_LOW_MR795: [i16; 2 * L_SUBFR] = doubled(&PH_IMP_LOW_MR795);
const IMP2_MID_MR795: [i16; 2 * L_SUBFR] = doubled(&PH_IMP_MID_MR795);

/// `dst[j] = sat(dst[j] + (c * imp[j]) >> 15)`; true if any add saturated.
#[inline(never)]
fn add_scaled(dst: &mut [i16; L_SUBFR], imp: &[i16; L_SUBFR], c: i16) -> bool {
    let mut o = 0u16;
    for (d, &h) in dst.iter_mut().zip(imp) {
        let t = ((c as i32 * h as i32) >> 15) as i16;
        let (s, f) = d.overflowing_add(t);
        o |= f as u16;
        *d = if f {
            if t < 0 { i16::MIN } else { i16::MAX }
        } else {
            s
        };
    }
    o != 0
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_c {
    use super::*;
    use crate::testutil::{N, Rng};
    use opencore_sys as _;

    /// Mirror of C `ph_dispState`.
    #[repr(C)]
    #[derive(Default)]
    struct CState {
        gain_mem: [i16; PHDGAINMEMSIZE],
        prev_state: i16,
        prev_cb_gain: i16,
        lock_full: i16,
        onset: i16,
    }

    /// Opaque storage for C `CommonAmrTbls` (34 table pointers), filled by `get_const_tbls`.
    #[repr(C)]
    struct CTbls([*const i16; 34]);

    // The TS 26.073 reference, in the shape of the opencore functions these tests were written
    // against (`o`: its overflow flag; tables unused). opencore only for its table set.
    unsafe extern "C" {
        fn ref_ph_disp_reset(s: *mut CState) -> i32;
        fn ref_ph_disp_lock(s: *mut CState);
        fn ref_ph_disp_release(s: *mut CState);
        fn ref_ph_disp(
            s: *mut CState,
            mode: i32,
            x: *mut i16,
            cb_gain: i16,
            ltp_gain: i16,
            inno: *mut i16,
            pitch_fac: i16,
            tmp_shift: i16,
        );
        fn get_const_tbls(t: *mut std::ffi::c_void); // opaque: also declared in lpc/gains tests
    }
    use crate::testutil::with_ref;
    unsafe fn c_ph_disp_reset(s: *mut CState) -> i16 {
        with_ref(std::ptr::null_mut(), || unsafe { ref_ph_disp_reset(s) }) as i16
    }
    unsafe fn c_ph_disp_lock(s: *mut CState) {
        with_ref(std::ptr::null_mut(), || unsafe { ref_ph_disp_lock(s) })
    }
    unsafe fn c_ph_disp_release(s: *mut CState) {
        with_ref(std::ptr::null_mut(), || unsafe { ref_ph_disp_release(s) })
    }
    #[allow(clippy::too_many_arguments)]
    unsafe fn c_ph_disp(
        s: *mut CState,
        mode: i32,
        x: *mut i16,
        cb_gain: i16,
        ltp_gain: i16,
        inno: *mut i16,
        pitch_fac: i16,
        tmp_shift: i16,
        _tbls: *mut CTbls,
        o: *mut i32,
    ) {
        with_ref(o, || unsafe {
            ref_ph_disp(s, mode, x, cb_gain, ltp_gain, inno, pitch_fac, tmp_shift)
        })
    }

    fn same(r: &PhDispState, c: &CState) -> bool {
        r.gain_mem == c.gain_mem
            && r.prev_state == c.prev_state
            && r.prev_cb_gain == c.prev_cb_gain
            && r.lock_full == c.lock_full
            && r.onset == c.onset
    }

    /// Innovation vector: sparse pulses (as the codebooks produce) or dense random.
    fn inno(g: &mut Rng) -> [i16; L_SUBFR] {
        let mut v = [0i16; L_SUBFR];
        match g.u64() % 4 {
            0 => v.iter_mut().for_each(|s| *s = g.i16()),
            1 => {} // all zero
            _ => {
                for _ in 0..1 + g.u64() % 10 {
                    let amp = [8191i16, -8192, 4096, -4096, 16382, -16382][(g.u64() % 6) as usize];
                    v[(g.u64() % 40) as usize] = amp;
                }
            }
        }
        v
    }

    /// Gain value: realistic ranges around the thresholds plus full-range edge values.
    fn gain(g: &mut Rng, centre: i16) -> i16 {
        match g.u64() % 4 {
            0 => g.i16(),
            1 => (g.u64() % 20) as i16, // below the cbGain < 10 threshold
            _ => (centre as i32 + (g.u64() % 16384) as i32 - 8192).clamp(0, 32767) as i16,
        }
    }

    #[test]
    fn ph_disp_vs_c() {
        let mut g = Rng::new(30);
        let mut tbls = CTbls([core::ptr::null(); 34]);
        unsafe { get_const_tbls((&raw mut tbls).cast()) };
        let mut r = PhDispState::default();
        let mut c = CState::default();
        assert_eq!(ph_disp_reset(&mut r), unsafe { c_ph_disp_reset(&mut c) });
        for n in 0..N {
            match g.u64() % 64 {
                0 => {
                    r.reset();
                    unsafe { c_ph_disp_reset(&mut c) };
                }
                1..=4 => {
                    ph_disp_lock(&mut r);
                    unsafe { c_ph_disp_lock(&mut c) };
                }
                5..=12 => {
                    ph_disp_release(&mut r);
                    unsafe { c_ph_disp_release(&mut c) };
                }
                _ => {}
            }
            assert!(same(&r, &c), "state after lock/release/reset, call {n}");

            // Modes MR475..MR122 (0..=7); the slow-varying mode exercises state history.
            let mode = Mode::ALL[(g.u64() % 8) as usize];
            let cb_gain = gain(&mut g, 2000);
            let ltp_gain = gain(&mut g, 12000);
            let pitch_fac = gain(&mut g, 8000);
            let tmp_shift = [1i16, 2, 0][(g.u64() % 3) as usize];
            let x0: [i16; L_SUBFR] = core::array::from_fn(|_| {
                if g.u64().is_multiple_of(2) {
                    g.i16()
                } else {
                    (g.u64() % 8192) as i16 - 4096
                }
            });
            let inno0 = inno(&mut g);

            let (mut rx, mut ri, mut ro) = (x0, inno0, false);
            ph_disp(
                &mut r, mode, &mut rx, cb_gain, ltp_gain, &mut ri, pitch_fac, tmp_shift, &mut ro,
            );
            let (mut cx, mut ci, mut co) = (x0, inno0, 0);
            unsafe {
                c_ph_disp(
                    &mut c,
                    mode as i32,
                    cx.as_mut_ptr(),
                    cb_gain,
                    ltp_gain,
                    ci.as_mut_ptr(),
                    pitch_fac,
                    tmp_shift,
                    &mut tbls,
                    &mut co,
                )
            };
            let ctx = format!(
                "call {n}: mode={mode} cb={cb_gain} ltp={ltp_gain} pf={pitch_fac} sh={tmp_shift}"
            );
            assert_eq!(rx, cx, "x, {ctx}");
            assert_eq!(ri, ci, "inno, {ctx}");
            let _ = co; // the reference's flag is global: values and state only
            assert!(same(&r, &c), "state, {ctx}: {r:?}");
        }
    }

    // The model's operators: TS 26.073's (the reference's basic operators, `o` its flag).
    unsafe extern "C" {
        fn ref_L_mult(a: i16, b: i16) -> i32;
        fn ref_L_mac(a: i32, b: i16, c: i16) -> i32;
        fn ref_L_shl(a: i32, b: i16) -> i32;
        fn ref_round(a: i32) -> i16;
    }
    #[allow(non_snake_case)]
    unsafe fn shim_L_mult(a: i16, b: i16, o: *mut i32) -> i32 {
        crate::testutil::with_ref(o, || unsafe { ref_L_mult(a, b) })
    }
    #[allow(non_snake_case)]
    unsafe fn shim_L_mac(a: i32, b: i16, c: i16, o: *mut i32) -> i32 {
        crate::testutil::with_ref(o, || unsafe { ref_L_mac(a, b, c) })
    }
    #[allow(non_snake_case)]
    unsafe fn shim_L_shl(a: i32, b: i16, o: *mut i32) -> i32 {
        crate::testutil::with_ref(o, || unsafe { ref_L_shl(a, b) })
    }
    unsafe fn pv_round(a: i32, o: *mut i32) -> i16 {
        crate::testutil::with_ref(o, || unsafe { ref_round(a) })
    }

    /// `total_exc` (both forms) vs the C operators composed per sample, as `ph_disp` and
    /// `Decoder_amr` run them. Edge-heavy inputs; classes counted: a 0x4000_0000 product in
    /// `L_mult` and in `L_mac`, `L_mac`'s unsaturated equal-operands overflow (accumulator ==
    /// product, the sum wraps), the `L_add`/`L_mac` saturation, `L_shl`'s, round's overflow.
    #[test]
    fn total_exc_vs_c() {
        let mut g = Rng::new(77);
        let pick = |g: &mut Rng| -> i16 {
            match g.u64() % 6 {
                0 => [-32768, 32767, 0, 1, -1, 16384][(g.u64() % 6) as usize],
                1 => (g.u64() % 512) as i16 - 256,
                _ => g.i16(),
            }
        };
        let mut cls = [0u32; 6];
        for case in 0..N {
            let mac = case % 2 == 0;
            let shift = [1i16, 2, 1, 2, 0, -1, 3, -40][(g.u64() % 8) as usize];
            let (pf, gn) = (pick(&mut g), pick(&mut g));
            let mut x: [i16; L_SUBFR] = core::array::from_fn(|_| pick(&mut g));
            let mut c: [i16; L_SUBFR] = core::array::from_fn(|_| pick(&mut g));
            // the equal-operands quirk: a = 2 x pf == r = c g, and the sum a + 2r = 6 x pf
            // wraps (x = c near 32767 needs pf > 2^31 / (6 * 32767) ~ 10923)
            if mac && case % 8 == 0 && (10_923..16384).contains(&pf) {
                let gq = 2 * pf;
                for i in 0..L_SUBFR / 2 {
                    (x[i], c[i]) = (32767 - i as i16, 32767 - i as i16);
                }
                let mut xr = x;
                let (r, co) = run_c(&mut xr, &c, pf, gq, shift, mac, &mut cls);
                let mut xs = x;
                total_exc::<true>(&mut xs, &c, pf, gq, shift);
                let _ = co; // the reference's flag is global: values only
                assert_eq!(xs, r, "case {case} (opencore's equal-operands quirk)");
            }
            let mut xr = x;
            let (r, co) = run_c(&mut xr, &c, pf, gn, shift, mac, &mut cls);
            if mac {
                total_exc::<true>(&mut x, &c, pf, gn, shift)
            } else {
                total_exc::<false>(&mut x, &c, pf, gn, shift)
            };
            let _ = co;
            assert_eq!(x, r, "case {case}: mac {mac} pf {pf} g {gn} shift {shift}");
        }
        assert!(cls.iter().all(|&n| n > 1000), "classes {cls:?}");
    }

    /// The C operators per sample; counts the classes into `cls`.
    fn run_c(
        x: &mut [i16; L_SUBFR],
        c: &[i16; L_SUBFR],
        pf: i16,
        g: i16,
        shift: i16,
        mac: bool,
        cls: &mut [u32; 6],
    ) -> ([i16; L_SUBFR], bool) {
        let mut o = 0i32;
        for (xi, &ci) in x.iter_mut().zip(c) {
            let p = *xi as i32 * pf as i32;
            let r = ci as i32 * g as i32;
            cls[0] += (p == 0x4000_0000) as u32;
            unsafe {
                let a = shim_L_mult(*xi, pf, &mut o);
                // TS 26.073 accumulates with L_mac at both call sites (`mac` no longer matters)
                if mac {
                    cls[1] += (r == 0x4000_0000) as u32;
                    cls[2] += (a == r && a.checked_add(r << 1).is_none()) as u32;
                }
                let s = shim_L_mac(a, ci, g, &mut o);
                cls[3] += (a.checked_add(r.wrapping_shl(1)).is_none()) as u32;
                let t = shim_L_shl(s, shift, &mut o);
                cls[4] += (shift > 0 && s.checked_mul(1 << shift).is_none_or(|v| v != t)) as u32;
                cls[5] += (t > i32::MAX - 0x8000) as u32;
                *xi = pv_round(t, &mut o);
            }
        }
        (*x, o != 0)
    }
}
