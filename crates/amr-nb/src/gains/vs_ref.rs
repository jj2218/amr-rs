//! Differential tests of the gains group against the TS 26.073 reference (`ref_*`); opencore
//! only for its table set and two state helpers (see the `extern` block).
use super::*;
use crate::testutil::{N, Rng};
use opencore_sys as _;
use ref26073_sys::lock;

unsafe extern "C" {
    fn ref_gc_pred_reset(st: *mut CGcPred) -> i16;
    fn ref_gc_pred(
        st: *mut CGcPred,
        mode: i32,
        code: *mut i16,
        exp_gcode0: *mut i16,
        frac_gcode0: *mut i16,
        exp_en: *mut i16,
        frac_en: *mut i16,
    );
    fn ref_gc_pred_update(st: *mut CGcPred, qua_ener_mr122: i16, qua_ener: i16);
    fn ref_gc_pred_average_limited(st: *mut CGcPred, a: *mut i16, b: *mut i16);
    fn ref_d_gain_code(st: *mut CGcPred, mode: i32, index: i16, code: *mut i16, g: *mut i16);
    fn ref_Dec_gain(
        st: *mut CGcPred,
        mode: i32,
        index: i16,
        code: *mut i16,
        even_subfr: i16,
        gain_pit: *mut i16,
        gain_cod: *mut i16,
    );
}

#[repr(C)]
#[derive(Default)]
struct CGcPred {
    past_qua_en: [i16; 4],
    past_qua_en_mr122: [i16; 4],
}

#[repr(C)]
#[derive(Default)]
struct CEcPitch {
    pbuf: [i16; 5],
    past_gain_pit: i16,
    prev_gp: i16,
}

#[repr(C)]
#[derive(Default)]
struct CEcCode {
    gbuf: [i16; 5],
    past_gain_code: i16,
    prev_gc: i16,
}

#[repr(C)]
#[derive(Default)]
struct CCbAver {
    cb_gain_history: [i16; 7],
    hang_var: i16,
    hang_count: i16,
}

#[repr(C)]
struct CBgnScd {
    frame_energy_hist: [i16; 60],
    bg_hangover: i16,
}

/// `CommonAmrTbls`: 34 table pointers.
type CTbls = [*const i16; 34];
const TBL_QUA_GAIN_PITCH: usize = 18;

// The TS 26.073 reference, in the shape of the opencore functions these tests were written
// against (`ov`/`o`: its overflow flag, sticky; tables unused). opencore only for its table
// set and the gc_pred/Dec_gain state helpers above.
unsafe extern "C" {
    #[link_name = "get_const_tbls"]
    fn c_get_const_tbls(t: *mut std::ffi::c_void); // opaque: lpc.rs declares it too
    #[link_name = "div_s"]
    fn oc_div_s(a: i16, b: i16) -> i16;
    #[link_name = "Cb_gain_average"]
    fn oc_cb_gain_average(
        st: *mut CCbAver,
        mode: i32,
        gain_code: i16,
        lsp: *mut i16,
        lsp_aver: *mut i16,
        bfi: i16,
        prev_bf: i16,
        pdfi: i16,
        prev_pdf: i16,
        in_bg_noise: i16,
        voiced_hangover: i16,
        ov: *mut i32,
    ) -> i16;
    fn ref_gmed_n(ind: *mut i16, n: i16) -> i16;
    fn ref_d_gain_pitch(mode: i32, index: i16) -> i16;
    fn ref_ec_gain_code_reset(st: *mut CEcCode) -> i32;
    fn ref_ec_gain_code(st: *mut CEcCode, pred: *mut CGcPred, state: i16, gain: *mut i16);
    fn ref_ec_gain_code_update(st: *mut CEcCode, bfi: i16, prev_bf: i16, g: *mut i16);
    fn ref_ec_gain_pitch_reset(st: *mut CEcPitch) -> i32;
    fn ref_ec_gain_pitch(st: *mut CEcPitch, state: i16, g: *mut i16);
    fn ref_ec_gain_pitch_update(st: *mut CEcPitch, bfi: i16, prev_bf: i16, g: *mut i16);
    fn ref_Cb_gain_average_reset(st: *mut CCbAver) -> i16;
    fn ref_Cb_gain_average(
        st: *mut CCbAver,
        mode: i32,
        gain_code: i16,
        lsp: *mut i16,
        lsp_aver: *mut i16,
        bfi: i16,
        prev_bf: i16,
        pdfi: i16,
        prev_pdf: i16,
        in_bg_noise: i16,
        voiced_hangover: i16,
    ) -> i16;
    fn ref_sub(a: i16, b: i16) -> i16;
    fn ref_add(a: i16, b: i16) -> i16;
    fn ref_shr(a: i16, b: i16) -> i16;
    fn ref_shl(a: i16, b: i16) -> i16;
    fn ref_norm_s(a: i16) -> i16;
    fn ref_div_s(a: i16, b: i16) -> i16;
    fn ref_abs_s(a: i16) -> i16;
    fn ref_Ex_ctrl(
        exc: *mut i16,
        exc_energy: i16,
        hist: *mut i16,
        voiced_hangover: i16,
        prev_bfi: i16,
        careful: i16,
    ) -> i16;
    fn ref_Bgn_scd_reset(st: *mut CBgnScd) -> i16;
    fn ref_Bgn_scd(
        st: *mut CBgnScd,
        ltp_gain_hist: *mut i16,
        speech: *mut i16,
        voiced_hangover: *mut i16,
    ) -> i16;
}

use crate::testutil::with_ref;
const NO: *mut i32 = std::ptr::null_mut();
unsafe fn c_gmed_n(ind: *mut i16, n: i16) -> i16 {
    with_ref(NO, || unsafe { ref_gmed_n(ind, n) })
}
unsafe fn c_d_gain_pitch(mode: i32, index: i16, _tbl: *const i16) -> i16 {
    with_ref(NO, || unsafe { ref_d_gain_pitch(mode, index) })
}
unsafe fn c_ec_gain_code_reset(st: *mut CEcCode) -> i16 {
    with_ref(NO, || unsafe { ref_ec_gain_code_reset(st) }) as i16
}
unsafe fn c_ec_gain_code(
    st: *mut CEcCode,
    pred: *mut CGcPred,
    state: i16,
    g: *mut i16,
    o: *mut i32,
) {
    with_ref(o, || unsafe { ref_ec_gain_code(st, pred, state, g) })
}
unsafe fn c_ec_gain_code_update(
    st: *mut CEcCode,
    bfi: i16,
    prev_bf: i16,
    g: *mut i16,
    o: *mut i32,
) {
    with_ref(o, || unsafe {
        ref_ec_gain_code_update(st, bfi, prev_bf, g)
    })
}
unsafe fn c_ec_gain_pitch_reset(st: *mut CEcPitch) -> i16 {
    with_ref(NO, || unsafe { ref_ec_gain_pitch_reset(st) }) as i16
}
unsafe fn c_ec_gain_pitch(st: *mut CEcPitch, state: i16, g: *mut i16, o: *mut i32) {
    with_ref(o, || unsafe { ref_ec_gain_pitch(st, state, g) })
}
unsafe fn c_ec_gain_pitch_update(
    st: *mut CEcPitch,
    bfi: i16,
    prev_bf: i16,
    g: *mut i16,
    o: *mut i32,
) {
    with_ref(o, || unsafe {
        ref_ec_gain_pitch_update(st, bfi, prev_bf, g)
    })
}
unsafe fn c_cb_gain_average_reset(st: *mut CCbAver) -> i16 {
    with_ref(NO, || unsafe { ref_Cb_gain_average_reset(st) })
}
#[allow(clippy::too_many_arguments)]
unsafe fn c_cb_gain_average(
    st: *mut CCbAver,
    mode: i32,
    gain_code: i16,
    lsp: *mut i16,
    lsp_aver: *mut i16,
    bfi: i16,
    prev_bf: i16,
    pdfi: i16,
    prev_pdf: i16,
    in_bg_noise: i16,
    voiced_hangover: i16,
    o: *mut i32,
) -> i16 {
    // A non-positive average LSP makes its div_s leave the contract (the reference aborts):
    // those calls on opencore, whose result the port keeps.
    let aver = unsafe { std::slice::from_raw_parts(lsp_aver, 10) };
    if aver.iter().any(|&a| a <= 0) {
        return unsafe {
            oc_cb_gain_average(
                st,
                mode,
                gain_code,
                lsp,
                lsp_aver,
                bfi,
                prev_bf,
                pdfi,
                prev_pdf,
                in_bg_noise,
                voiced_hangover,
                o,
            )
        };
    }
    with_ref(o, || unsafe {
        ref_Cb_gain_average(
            st,
            mode,
            gain_code,
            lsp,
            lsp_aver,
            bfi,
            prev_bf,
            pdfi,
            prev_pdf,
            in_bg_noise,
            voiced_hangover,
        )
    })
}
unsafe fn c_sub(a: i16, b: i16, o: *mut i32) -> i16 {
    with_ref(o, || unsafe { ref_sub(a, b) })
}
unsafe fn add_16(a: i16, b: i16, o: *mut i32) -> i16 {
    with_ref(o, || unsafe { ref_add(a, b) })
}
unsafe fn c_shr(a: i16, b: i16, o: *mut i32) -> i16 {
    with_ref(o, || unsafe { ref_shr(a, b) })
}
unsafe fn c_norm_s(a: i16) -> i16 {
    with_ref(NO, || unsafe { ref_norm_s(a) })
}
/// `div_s` in its contract (0 <= a <= b, b > 0) on the reference; outside it (the reference
/// aborts, the spec defines nothing) on opencore, whose result the port keeps.
unsafe fn c_div_s(a: i16, b: i16) -> i16 {
    if 0 <= a && a <= b && b > 0 {
        with_ref(NO, || unsafe { ref_div_s(a, b) })
    } else {
        unsafe { oc_div_s(a, b) }
    }
}
unsafe fn shim_shl(a: i16, b: i16, o: *mut i32) -> i16 {
    with_ref(o, || unsafe { ref_shl(a, b) })
}
unsafe fn shim_abs_s(a: i16) -> i16 {
    with_ref(NO, || unsafe { ref_abs_s(a) })
}
unsafe fn c_ex_ctrl(
    exc: *mut i16,
    exc_energy: i16,
    hist: *mut i16,
    voiced_hangover: i16,
    prev_bfi: i16,
    careful: i16,
    o: *mut i32,
) -> i16 {
    with_ref(o, || unsafe {
        ref_Ex_ctrl(exc, exc_energy, hist, voiced_hangover, prev_bfi, careful)
    })
}
unsafe fn c_bgn_scd_reset(st: *mut CBgnScd) -> i16 {
    with_ref(NO, || unsafe { ref_Bgn_scd_reset(st) })
}
unsafe fn c_bgn_scd(
    st: *mut CBgnScd,
    ltp_gain_hist: *mut i16,
    speech: *mut i16,
    voiced_hangover: *mut i16,
    o: *mut i32,
) -> i16 {
    with_ref(o, || unsafe {
        ref_Bgn_scd(st, ltp_gain_hist, speech, voiced_hangover)
    })
}

fn c_tbls() -> CTbls {
    let mut t: CTbls = [std::ptr::null(); 34];
    unsafe { c_get_const_tbls((&raw mut t).cast()) };
    t
}

fn chk_pred(r: &GcPredState, c: &CGcPred, ctx: &str) {
    assert_eq!(r.past_qua_en, c.past_qua_en, "past_qua_en: {ctx}");
    assert_eq!(
        r.past_qua_en_mr122, c.past_qua_en_mr122,
        "past_qua_en_MR122: {ctx}"
    );
}

/// Uniform value in `lo..=hi`.
fn range(g: &mut Rng, lo: i32, hi: i32) -> i16 {
    (lo + (g.u64() % (hi - lo + 1) as u64) as i32) as i16
}

/// An innovative-codebook-like vector: sparse pulses, small noise, or full-range values.
fn code_vec(g: &mut Rng) -> [i16; 40] {
    let mut c = [0i16; 40];
    match g.u64() % 5 {
        0 => {
            for _ in 0..(g.u64() % 11) {
                let amp = [4096, 8191, 2048, 8192, 16383][(g.u64() % 5) as usize];
                c[(g.u64() % 40) as usize] = if g.u64() & 1 == 0 { amp } else { -amp };
            }
        }
        1 => c.iter_mut().for_each(|x| *x = range(g, -64, 64)),
        2 => c.iter_mut().for_each(|x| *x = range(g, -2000, 2000)),
        3 => {} // all zero
        _ => c.iter_mut().for_each(|x| *x = g.i16()),
    }
    c
}

/// A plausible quantised-energy pair (sometimes full range).
fn qua_pair(g: &mut Rng) -> (i16, i16) {
    if g.u64().is_multiple_of(8) {
        (g.i16(), g.i16())
    } else {
        (range(g, -5000, 3000), range(g, -16000, 10000))
    }
}

#[test]
fn gmed_n_vs_ref() {
    let mut g = Rng::new(1);
    for _ in 0..N {
        let n = (1 + g.u64() % 9) as i16;
        let mut v = [0i16; 9];
        let small = g.u64() & 1 == 0;
        for x in v.iter_mut() {
            *x = if small { range(&mut g, -3, 3) } else { g.i16() };
            if g.u64().is_multiple_of(10) {
                *x = i16::MIN;
            }
        }
        let r = gmed_n(&v, n);
        let c = unsafe { c_gmed_n(v.as_mut_ptr(), n) };
        assert_eq!(r, c, "gmed_n({v:?}, {n})");
    }
}

/// Against the reference (opencore's `ener_code` wraps for large codebook energies).
#[test]
fn gc_pred_vs_ref() {
    let _g = lock();
    let mut g = Rng::new(2);
    let mut r = GcPredState::default();
    let mut c = CGcPred::default();
    unsafe { ref_gc_pred_reset(&mut c) };
    chk_pred(&r, &c, "reset");
    for k in 0..N {
        if k % 3000 == 0 {
            r.reset();
            unsafe { ref_gc_pred_reset(&mut c) };
        }
        let mode = CodecMode::ALL[(g.u64() % 8) as usize];
        let mut code = code_vec(&mut g);
        let mut ro = false;
        let out = gc_pred(&r, mode, &code, &mut ro);
        let (mut e, mut f, mut ee, mut fe) = (0, 0, 0, 0);
        unsafe {
            ref_gc_pred(
                &mut c,
                mode as i32,
                code.as_mut_ptr(),
                &mut e,
                &mut f,
                &mut ee,
                &mut fe,
            )
        };
        let ctx = format!("call {k} mode {mode:?} code {code:?} state {r:?}");
        assert_eq!((out.0, out.1), (e, f), "gc_pred: {ctx}");
        if mode == MR795 {
            assert_eq!((out.2, out.3), (ee, fe), "gc_pred en: {ctx}");
        }

        let (ra, rb) = gc_pred_average_limited(&r, &mut ro);
        let (mut ca, mut cb) = (0, 0);
        unsafe { ref_gc_pred_average_limited(&mut c, &mut ca, &mut cb) };
        assert_eq!((ra, rb), (ca, cb), "average_limited: {ctx}");

        let (q1, q2) = qua_pair(&mut g);
        gc_pred_update(&mut r, q1, q2);
        unsafe { ref_gc_pred_update(&mut c, q1, q2) };
        chk_pred(&r, &c, &ctx);
    }
}

#[test]
fn d_gain_pitch_vs_ref() {
    let t = c_tbls();
    for mode in CodecMode::ALL {
        for index in 0..16 {
            let c = unsafe { c_d_gain_pitch(mode as i32, index, t[TBL_QUA_GAIN_PITCH]) };
            assert_eq!(
                d_gain_pitch(mode, index),
                c,
                "d_gain_pitch({mode:?}, {index})"
            );
        }
    }
}

#[test]
fn d_gain_code_vs_ref() {
    let _g = lock();
    let mut g = Rng::new(3);
    let mut r = GcPredState::default();
    let mut c = CGcPred::default();
    unsafe { ref_gc_pred_reset(&mut c) };
    for k in 0..N {
        if k % 2500 == 0 {
            r.reset();
            unsafe { ref_gc_pred_reset(&mut c) };
        }
        // Decoder uses MR795/MR122; cover all speech modes (index is masked with 31 in C).
        let mode = if g.u64() & 1 == 0 {
            [MR795, MR122][(g.u64() % 2) as usize]
        } else {
            CodecMode::ALL[(g.u64() % 8) as usize]
        };
        // The reference indexes its table with the raw index: in range only.
        let index = range(&mut g, 0, 31);
        let mut code = code_vec(&mut g);
        let mut ro = false;
        let rg = d_gain_code(&mut r, mode, index, &code, &mut ro);
        let mut cg = 0;
        unsafe { ref_d_gain_code(&mut c, mode as i32, index, code.as_mut_ptr(), &mut cg) };
        let ctx = format!("call {k} mode {mode:?} index {index}");
        assert_eq!(rg, cg, "d_gain_code: {ctx}");
        chk_pred(&r, &c, &ctx);
        // Occasionally inject arbitrary predictor history (as ec_gain_code would).
        if g.u64().is_multiple_of(20) {
            let (q1, q2) = qua_pair(&mut g);
            gc_pred_update(&mut r, q1, q2);
            unsafe { ref_gc_pred_update(&mut c, q1, q2) };
        }
    }
}

#[test]
fn dec_gain_vs_ref() {
    let _g = lock();
    let mut g = Rng::new(4);
    let mut r = GcPredState::default();
    let mut c = CGcPred::default();
    unsafe { ref_gc_pred_reset(&mut c) };
    // Dec_gain passes NULL exp_en/frac_en to gc_pred: MR795 (and MR122) are not valid here.
    let modes = [MR475, MR515, MR59, MR67, MR74, MR102];
    for k in 0..N {
        if k % 2000 == 0 {
            r.reset();
            unsafe { ref_gc_pred_reset(&mut c) };
        }
        let mode = modes[(g.u64() % 6) as usize];
        let max_index = match mode {
            MR475 => 255,
            MR515 | MR59 => 63,
            _ => 127,
        };
        let index = range(&mut g, 0, max_index);
        let even = (g.u64() & 1) as i16;
        let mut code = code_vec(&mut g);
        let mut ro = false;
        let (rp, rc) = dec_gain(&mut r, mode, index, &code, even, &mut ro);
        let (mut cp, mut cc) = (0, 0);
        unsafe {
            ref_Dec_gain(
                &mut c,
                mode as i32,
                index,
                code.as_mut_ptr(),
                even,
                &mut cp,
                &mut cc,
            )
        };
        let ctx = format!("call {k} mode {mode:?} index {index} even {even}");
        assert_eq!((rp, rc), (cp, cc), "Dec_gain: {ctx}");
        chk_pred(&r, &c, &ctx);
    }
}

#[test]
fn ec_gains_vs_ref() {
    let mut g = Rng::new(5);
    let (mut rp, mut rc, mut rpred) = (
        EcGainPitchState::default(),
        EcGainCodeState::default(),
        GcPredState::default(),
    );
    let (mut cp, mut cc, mut cpred) = (CEcPitch::default(), CEcCode::default(), CGcPred::default());
    let mut prev_bf = 0i16;
    for k in 0..N {
        if k % 4000 == 0 {
            rp.reset();
            rc.reset();
            rpred.reset();
            unsafe {
                c_ec_gain_pitch_reset(&mut cp);
                c_ec_gain_code_reset(&mut cc);
                ref_gc_pred_reset(&mut cpred);
            }
        }
        // Bad frames come in bursts.
        let bfi = if prev_bf != 0 {
            (!g.u64().is_multiple_of(3)) as i16
        } else {
            g.u64().is_multiple_of(6) as i16
        };
        let state = range(&mut g, 0, 6);
        let ctx = format!("call {k} bfi {bfi} prev_bf {prev_bf} state {state}");
        let (mut ro, mut co) = (false, 0);
        let (mut gp, mut gc);
        if bfi != 0 {
            gp = ec_gain_pitch(&rp, state, &mut ro);
            let mut c_gp = 0;
            unsafe { c_ec_gain_pitch(&mut cp, state, &mut c_gp, &mut co) };
            assert_eq!(gp, c_gp, "ec_gain_pitch: {ctx}");
            gc = ec_gain_code(&rc, &mut rpred, state, &mut ro);
            let mut c_gc = 0;
            unsafe { c_ec_gain_code(&mut cc, &mut cpred, state, &mut c_gc, &mut co) };
            assert_eq!(gc, c_gc, "ec_gain_code: {ctx}");
        } else {
            gp = if g.u64().is_multiple_of(10) {
                g.i16()
            } else {
                range(&mut g, 0, 20000)
            };
            gc = if g.u64().is_multiple_of(10) {
                g.i16()
            } else {
                range(&mut g, 0, 12000)
            };
            let (q1, q2) = qua_pair(&mut g);
            gc_pred_update(&mut rpred, q1, q2);
            unsafe { ref_gc_pred_update(&mut cpred, q1, q2) };
        }
        let (mut c_gp, mut c_gc) = (gp, gc);
        ec_gain_pitch_update(&mut rp, bfi, prev_bf, &mut gp, &mut ro);
        ec_gain_code_update(&mut rc, bfi, prev_bf, &mut gc, &mut ro);
        unsafe {
            c_ec_gain_pitch_update(&mut cp, bfi, prev_bf, &mut c_gp, &mut co);
            c_ec_gain_code_update(&mut cc, bfi, prev_bf, &mut c_gc, &mut co);
        }
        let _ = (ro, co); // the reference's flag is global: values and state only
        assert_eq!((gp, gc), (c_gp, c_gc), "updates: {ctx}");
        assert_eq!(
            (rp.pbuf, rp.past_gain_pit, rp.prev_gp),
            (cp.pbuf, cp.past_gain_pit, cp.prev_gp),
            "pitch state: {ctx}"
        );
        assert_eq!(
            (rc.gbuf, rc.past_gain_code, rc.prev_gc),
            (cc.gbuf, cc.past_gain_code, cc.prev_gc),
            "code state: {ctx}"
        );
        chk_pred(&rpred, &cpred, &ctx);
        prev_bf = bfi;
    }
}

/// `Cb_gain_average`'s `diff` loop as the C runs it: its operators, coefficient by coefficient.
/// Counts the classes into `cls`.
fn lsp_dev_sum_c(aver: &[i16; 10], lsp: &[i16; 10], cls: &mut [u32; 7]) -> (i16, bool) {
    let mut o = 0i32;
    let mut diff = 0i16;
    for i in 0..10 {
        unsafe {
            let mut so = 0;
            let mut t = shim_abs_s(c_sub(aver[i], lsp[i], &mut so));
            let s1 = c_norm_s(t) - 1;
            cls[0] += (so != 0) as u32; // sub saturates
            cls[1] += (aver[i] <= 0) as u32; // div_s of a non-positive denominator
            cls[2] += (t == 0) as u32;
            cls[3] += (t >= 16384) as u32; // n = t >> 1
            o |= so;
            t = shim_shl(t, s1, &mut o);
            let s2 = c_norm_s(aver[i]);
            let d = shim_shl(aver[i], s2, &mut o);
            let q = c_div_s(t, d);
            let shift = 2 + s1 - s2;
            let r = if shift >= 0 {
                c_shr(q, shift, &mut o)
            } else {
                let r = shim_shl(q, -shift, &mut o);
                let live = aver[i] > 0 && t > 0;
                cls[4] += (live && r == MAX_16) as u32; // shl saturates
                cls[5] += (live && r < MAX_16) as u32;
                r
            };
            let mut ao = 0;
            diff = add_16(diff, r, &mut ao);
            cls[6] += (ao != 0) as u32; // the sum saturates
            o |= ao;
        }
    }
    (diff, o != 0)
}

/// `lsp_dev_sum` (lanes) vs the C's loop. Each coefficient: aver mostly positive (all
/// magnitudes, powers of two ± 1), lsp at distances that cross every `norm_s` boundary, some
/// full-range (sub saturates); some calls with one live lane (a saturated `shl` term alone).
/// Every class is counted.
#[test]
fn lsp_dev_sum_vs_ref() {
    let mut g = Rng::new(61);
    let mut cls = [0u32; 7];
    let mut lone_sat = 0u32;
    let pow = |g: &mut Rng| -> i32 { (1i32 << (g.u64() % 16)) + (g.u64() % 3) as i32 - 1 };
    for case in 0..N {
        let big = case % 4 == 0; // large deviations: the sum saturates
        let mut aver = [0i16; 10];
        let mut lsp = [0i16; 10];
        for i in 0..10 {
            let a = match g.u64() % 8 {
                0 => g.i16() as i32,
                1 => pow(&mut g),
                2 => range(&mut g, 1, 300) as i32,
                _ => range(&mut g, 1, 32767) as i32,
            }
            .clamp(-32768, 32767);
            let dev = match g.u64() % 8 {
                0 => g.i16() as i32,
                1 => 0,
                2 if big => range(&mut g, -32768, 32767) as i32,
                _ => pow(&mut g) * if g.u64() & 1 == 0 { 1 } else { -1 },
            };
            aver[i] = a as i16;
            lsp[i] = if g.u64().is_multiple_of(16) {
                g.i16()
            } else {
                (a - dev).clamp(-32768, 32767) as i16
            };
        }
        // one live lane: a saturated term alone decides the sum
        let lone = case % 8 == 1;
        if lone {
            let keep = (g.u64() % 10) as usize;
            for i in (0..10).filter(|&i| i != keep) {
                lsp[i] = aver[i];
            }
        }
        let n7 = cls[4];
        let (cd, co) = lsp_dev_sum_c(&aver, &lsp, &mut cls);
        lone_sat += (lone && cls[4] > n7) as u32;
        let mut ro = false;
        let rd = lsp_dev_sum(&aver, &lsp, &mut ro);
        let _ = (ro, co); // the reference's flag is global: values only
        assert_eq!(rd, cd, "case {case}: aver {aver:?} lsp {lsp:?}");
    }
    assert!(
        cls.iter().all(|&n| n > 1000) && lone_sat > 1000,
        "classes {cls:?} {lone_sat}"
    );
}

#[test]
fn cb_gain_average_vs_ref() {
    let mut g = Rng::new(6);
    let mut r = CbGainAverageState::default();
    let mut c = CCbAver::default();
    unsafe { c_cb_gain_average_reset(&mut c) };
    let mut lsp_aver = [0i16; 10];
    let mut noise = 0i32;
    let mut mode = MR475;
    // Coverage: [bgMix < 8192 mixing used, hangVar > 10, bad-frame mean in background noise]
    let mut cov = [0usize; 3];
    for k in 0..N {
        if k % 5000 == 0 {
            r.reset();
            unsafe { c_cb_gain_average_reset(&mut c) };
        }
        // Regimes of ~60 frames: fixed mode, stationary/non-stationary spectrum.
        if k % 60 == 0 {
            mode = CodecMode::ALL[(g.u64() % 9) as usize];
            noise = [0, 30, 300, 3000, 20000][(g.u64() % 5) as usize];
            for (i, x) in lsp_aver.iter_mut().enumerate() {
                *x = if g.u64().is_multiple_of(16) {
                    g.i16()
                } else {
                    (30000 - 6000 * i as i32 + range(&mut g, -2000, 2000) as i32) as i16
                };
            }
        }
        let mut lsp = [0i16; 10];
        for i in 0..10 {
            let d = if noise == 0 {
                0
            } else {
                (g.u64() % (2 * noise as u64 + 1)) as i32 - noise
            };
            lsp[i] = (lsp_aver[i] as i32 + d).clamp(-32768, 32767) as i16;
        }
        let gain_code = if g.u64().is_multiple_of(10) {
            g.i16()
        } else {
            range(&mut g, 0, 8000)
        };
        let flags: Vec<i16> = (0..5).map(|_| (g.u64() & 1) as i16).collect();
        let (bfi, prev_bf, pdfi, prev_pdf, inbg) =
            (flags[0], flags[1], flags[2], flags[3], flags[4]);
        let vh = range(&mut g, 0, 10);
        let mixing = r.hang_count >= 39 && (mode <= MR67 || mode == MR102);
        let mut ro = false;
        let rr = cb_gain_average(
            &mut r, mode, gain_code, &lsp, &lsp_aver, bfi, prev_bf, pdfi, prev_pdf, inbg, vh,
            &mut ro,
        );
        let mut co = 0;
        let cr = unsafe {
            c_cb_gain_average(
                &mut c,
                mode as i32,
                gain_code,
                lsp.as_mut_ptr(),
                lsp_aver.as_mut_ptr(),
                bfi,
                prev_bf,
                pdfi,
                prev_pdf,
                inbg,
                vh,
                &mut co,
            )
        };
        let ctx =
            format!("call {k} mode {mode:?} lsp {lsp:?} aver {lsp_aver:?} flags {flags:?} vh {vh}");
        let _ = (ro, co);
        assert_eq!(rr, cr, "Cb_gain_average: {ctx}");
        assert_eq!(
            (r.cb_gain_history, r.hang_var, r.hang_count),
            (c.cb_gain_history, c.hang_var, c.hang_count),
            "state: {ctx}"
        );
        cov[0] += (mixing && r.hang_var == 0 && rr != gain_code) as usize;
        cov[1] += (r.hang_var > 10) as usize;
        cov[2] += ((bfi | prev_bf) != 0 && inbg != 0 && mode <= MR59) as usize;
    }
    assert!(cov.iter().all(|&n| n > 100), "coverage {cov:?}");
}

#[test]
fn ex_ctrl_vs_ref() {
    let mut g = Rng::new(7);
    let mut scaled = 0;
    for k in 0..N {
        let wide = g.u64().is_multiple_of(8);
        let mut hist = [0i16; 9];
        for x in hist.iter_mut() {
            *x = if wide {
                g.i16()
            } else {
                range(&mut g, 0, 4000)
            };
        }
        let exc_energy = if wide {
            g.i16()
        } else {
            range(&mut g, 0, 4000)
        };
        let mut exc = [0i16; 40];
        exc.iter_mut().for_each(|x| *x = g.i16());
        let (vh, prev_bfi, careful) = (
            range(&mut g, 0, 10),
            (g.u64() & 1) as i16,
            (g.u64() & 1) as i16,
        );
        let mut rexc = exc;
        let before = exc;
        let mut ro = false;
        let rr = ex_ctrl(&mut rexc, exc_energy, &hist, vh, prev_bfi, careful, &mut ro);
        let mut co = 0;
        let cr = unsafe {
            c_ex_ctrl(
                exc.as_mut_ptr(),
                exc_energy,
                hist.as_mut_ptr(),
                vh,
                prev_bfi,
                careful,
                &mut co,
            )
        };
        let ctx = format!(
            "case {k} energy {exc_energy} hist {hist:?} vh {vh} bfi {prev_bfi} careful {careful}"
        );
        let _ = (ro, co);
        assert_eq!((rr, rexc), (cr, exc), "Ex_ctrl: {ctx}");
        scaled += (rexc != before) as usize;
    }
    assert!(scaled > N / 10, "scaling branch hit {scaled} times");
}

#[test]
fn bgn_scd_vs_ref() {
    let mut g = Rng::new(8);
    let mut r = BgnScdState::default();
    let mut c = CBgnScd {
        frame_energy_hist: [0; 60],
        bg_hangover: 0,
    };
    unsafe { c_bgn_scd_reset(&mut c) };
    let mut ltp = [0i16; 9];
    let (mut rvh, mut cvh) = (0i16, 0i16);
    let (mut amp, mut ltp_hi) = (0i32, 0i32);
    let mut speech = [0i16; 160];
    // Coverage: [in background noise, bgHangover > 20, voiced (hangover reset), hangover = 10]
    let mut cov = [0usize; 4];
    for k in 0..N / 4 {
        if k % 5000 == 0 {
            r.reset();
            unsafe { c_bgn_scd_reset(&mut c) };
        }
        // Regimes of ~40 frames: silence, noise of various levels, loud speech, clipping.
        if k % 40 == 0 {
            amp = [0, 5, 30, 60, 150, 400, 1000, 5000, 32767][(g.u64() % 9) as usize];
            ltp_hi = [8000, 14000, 16000, 20000][(g.u64() % 4) as usize];
        }
        let a = (amp + amp * range(&mut g, -3, 3) as i32 / 10).clamp(0, 32767);
        for x in speech.iter_mut() {
            *x = if a == 0 { 0 } else { range(&mut g, -a, a) };
        }
        if amp == 32767 && g.u64() & 1 == 0 {
            speech.iter_mut().for_each(|x| *x = i16::MIN);
        }
        ltp.copy_within(1.., 0);
        ltp[8] = if g.u64().is_multiple_of(20) {
            g.i16()
        } else {
            range(&mut g, 0, ltp_hi)
        };
        let mut ro = false;
        let rr = bgn_scd(&mut r, &ltp, &speech, &mut rvh, &mut ro);
        let mut co = 0;
        let cr = unsafe {
            c_bgn_scd(
                &mut c,
                ltp.as_mut_ptr(),
                speech.as_mut_ptr(),
                &mut cvh,
                &mut co,
            )
        };
        let ctx = format!("call {k} amp {a} ltp {ltp:?}");
        let _ = (ro, co);
        assert_eq!((rr, rvh), (cr, cvh), "Bgn_scd: {ctx}");
        assert_eq!(r.frame_energy_hist, c.frame_energy_hist, "hist: {ctx}");
        assert_eq!(r.bg_hangover, c.bg_hangover, "bgHangover: {ctx}");
        cov[0] += (rr != 0) as usize;
        cov[1] += (r.bg_hangover > 20) as usize;
        cov[2] += (rvh == 0) as usize;
        cov[3] += (rvh == 10) as usize;
    }
    assert!(cov.iter().all(|&n| n > 100), "coverage {cov:?}");
}
