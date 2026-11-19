//! Kernel timings, Rust vs C, ns per call (`BENCH=kern just bench`, `KERN=<substring>`).
//! Whole-codec timing is too coarse for one kernel, and callgrind's instruction and mispredict
//! counts mislead.
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::basic_op::Flag;
use crate::cnst::{CodecMode, L_FRAME, M, PIT_MAX};
use crate::encoder::codebook::{
    RrTm, cor_h_search_8i40, cor_h_tm, cor_h_x, filtered_code, search_10i40, set_sign,
};
use crate::encoder::gain_q::{GainAdaptState, mr795_gain_quant};
use crate::encoder::layout::cols;
use crate::encoder::lpc::{DICO2_SOA, DICO3_SOA, vq_subvec3, vq_subvec4};
use crate::encoder::pitch_cl::convolve;
use crate::encoder::pitch_ol::comp_corr;
use crate::encoder::vs_ref::tables;
use crate::encoder::vs_ref::*;
use crate::encoder::{codebook::vs_ref::*, gain_q::vs_ref::*};
use crate::encoder::{pitch_cl::vs_ref::*, pitch_ol::vs_ref::*};
use crate::testutil::Rng;

/// ns per call of `f` over 64 rotating inputs (best of 20 runs of 20k calls). Only the
/// kernels whose name contains `$KERN`, if set.
fn time(name: &str, mut f: impl FnMut(usize)) {
    if std::env::var("KERN").is_ok_and(|k| !name.contains(&k)) {
        return;
    }
    let best = (0..20)
        .map(|_| {
            let t = std::time::Instant::now();
            for i in 0..20_000 {
                f(i & 63);
            }
            t.elapsed().as_secs_f64() / 20_000.0 * 1e9
        })
        .fold(f64::INFINITY, f64::min);
    eprintln!("{name:<24} {best:8.1} ns");
}

/// Kernel timings, Rust vs C, ns per call (`BENCH=kern just bench`). Whole-codec timing is
/// too coarse for one kernel, and callgrind's instruction and mispredict counts mislead.
#[test]
#[ignore]
fn kernel_timing() {
    use std::hint::black_box as bb;
    let mut g = Rng::new(99);
    let mut xs = [[0i16; 40]; 64];
    let mut hs = [[0i16; 40]; 64];
    for (x, h) in xs.iter_mut().zip(&mut hs) {
        // the encoder's kind of input (an excitation, a decaying impulse response): the
        // fast path; full-scale `signal`/`impulse` took the exact loop every call
        for v in x.iter_mut() {
            *v = g.i16() >> (3 + g.u64() % 4);
        }
        let scale = 2 + (g.u64() % 4) as u32;
        for (i, v) in h.iter_mut().enumerate() {
            *v = g.i16() >> (scale + (i as u32 / 6).min(12));
        }
        h[0] = 4096 + (g.u64() % 8192) as i16;
    }
    let mut y = [0i16; 40];
    time("convolve", |i| convolve(bb(&xs[i]), bb(&hs[i]), bb(&mut y)));
    time("convolve C", |i| unsafe {
        oc_convolve(
            bb(xs[i].as_ptr() as *mut i16),
            bb(hs[i].as_ptr() as *mut i16),
            bb(y.as_mut_ptr()),
            40,
        )
    });
    let mut syn = Vec::new();
    while syn.len() < 64 {
        // the encoder's kind of input: LPC residual-like x and past outputs, no saturation
        // (the exact fallback is rare in the codec)
        let a = crate::post::test_gen::lpc(&mut g);
        let x: [i16; 40] = std::array::from_fn(|_| (g.u64() as i16) >> (6 + g.u64() % 4));
        let mut mem = [0i16; M];
        let mut y = [0i16; 40];
        crate::post::filt::syn_filt(&a, &x, &mut y, 40, &mut mem, true);
        // and within `syn_filt`'s bound (a0 max|x| + sum|a_j| max|y| < 2^30): 1.6% of the
        // old inputs took the exact loop, the decoder's calls 39 in 714k
        let am = |v: &[i16]| v.iter().map(|&c| (c as i64).abs()).max().unwrap_or(0);
        let sa: i64 = a[1..].iter().map(|&c| (c as i64).abs()).sum();
        if (a[0] as i64).abs() * am(&x) + sa * am(&y) < 1 << 30 {
            syn.push((a, x, mem));
        }
    }
    let mut ys = [0i16; 40];
    time("syn_filt 40", |i| {
        let (a, x, mut mem) = syn[i];
        crate::post::filt::syn_filt(bb(&a), bb(&x), bb(&mut ys), 40, bb(&mut mem), true);
    });
    time("syn_filt 40 C", |i| unsafe {
        let (mut a, mut x, mut mem) = syn[i];
        c_syn_filt(
            bb(a.as_mut_ptr()),
            bb(x.as_mut_ptr()),
            bb(ys.as_mut_ptr()),
            40,
            bb(mem.as_mut_ptr()),
            1,
        )
    });
    let mut sigs = [[0i16; PIT_MAX + L_FRAME]; 8];
    for sig in &mut sigs {
        *sig = pitch_signal(&mut g);
    }
    {
        unsafe extern "C" {
            #[link_name = "Pred_lt_3or6"]
            fn c_pred_lt(exc: *mut i16, t0: i16, frac: i16, l: i16, flag3: i16, o: *mut i32);
        }
        // the encoder's call: l_subfr 41, lags over the whole range, both resolutions,
        // frac -2..=2 (all but the 1/6 resolution's frac 3)
        const OFF: usize = PIT_MAX + 11;
        let args: [(i16, i16, i16); 64] = std::array::from_fn(|i| {
            let t0 = 20 + (g.u64() % 124) as i16;
            (t0, (g.u64() % 5) as i16 - 2, (i & 1) as i16)
        });
        let mut buf = [0i16; OFF + 41];
        // an excitation's range (full-scale `sigs` failed the window guard: the exact loop)
        let excs = sigs.map(|s| s.map(|v| v >> 4));
        time("pred_lt_3or6", |i| {
            let (t0, frac, flag3) = args[i];
            buf[..PIT_MAX + 11].copy_from_slice(&excs[i & 7][..PIT_MAX + 11]);
            crate::excitation::pred_lt_3or6(
                bb(&mut buf),
                OFF,
                t0,
                frac,
                41,
                flag3,
                &mut Flag::default(),
            )
        });
        time("pred_lt_3or6 C", |i| unsafe {
            let (t0, frac, flag3) = args[i];
            buf[..PIT_MAX + 11].copy_from_slice(&excs[i & 7][..PIT_MAX + 11]);
            let mut o = 0;
            c_pred_lt(bb(buf.as_mut_ptr().add(OFF)), t0, frac, 41, flag3, &mut o)
        });
    }
    let mut corr = [0i32; PIT_MAX + 1];
    time("comp_corr 160", |i| {
        comp_corr(bb(&sigs[i & 7]), L_FRAME, PIT_MAX, 20, bb(&mut corr))
    });
    time("comp_corr 160 C", |i| {
        bb(c_comp_corr_at(bb(&mut sigs[i & 7]), L_FRAME, PIT_MAX, 20));
    });
    let pins: Vec<PitchIn> = (0..64).map(|_| pitch_in(&mut g)).collect();
    time("pitch_fr", |i| {
        bb(pitch_fr_rust(bb(&pins[i])));
    });
    time("pitch_fr C", |i| {
        bb(pitch_fr_c(bb(&pins[i])));
    });
    let mut lsf = [[0i16; 10]; 64];
    let mut wf = [[0i16; 10]; 64];
    for (l, w) in lsf.iter_mut().zip(&mut wf) {
        for (a, b) in l.iter_mut().zip(w) {
            (*a, *b) = (
                (g.u64() % 2000) as i16 - 1000,
                3000 + (g.u64() % 20000) as i16,
            );
        }
    }
    time("vq_subvec3 512", |i| {
        let mut r = lsf[i];
        bb(vq_subvec3(&mut r[3..], bb(&cols(&DICO2_SOA)), &wf[i][3..]));
    });
    time("vq_subvec4 512", |i| {
        let mut r = lsf[i];
        bb(vq_subvec4(&mut r[6..], bb(&cols(&DICO3_SOA)), &wf[i][6..]));
    });
    let (mut dn, mut ov) = ([0i16; 40], false);
    time("cor_h_x", |i| {
        cor_h_x(bb(&hs[i]), bb(&xs[i]), bb(&mut dn), 1, &mut ov)
    });
    let (mut sg, mut dn2) = ([0i16; 40], [0i16; 40]);
    time("set_sign", |i| {
        let mut d = xs[i];
        set_sign(bb(&mut d), &mut sg, &mut dn2, 4);
        bb(&dn2);
    });
    // LP analysis on speech-like windows: autocorr, and az_lsp on the LPC it leads to.
    let (wins, mut azs) = ([[0i16; 240]; 64], [[0i16; 11]; 64]);
    let mut wins = wins;
    let mut lev = crate::encoder::lpc::LevinsonState::default();
    for (x, a) in wins.iter_mut().zip(&mut azs) {
        signal(&mut g, x);
        let (mut rh, mut rl, mut rc, mut ov) = ([0i16; 11], [0i16; 11], [0i16; 4], false);
        crate::encoder::lpc::autocorr(x, 10, &mut rh, &mut rl, &crate::tables::WINDOW_200_40);
        crate::encoder::lpc::lag_window(10, &mut rh, &mut rl, &mut ov);
        crate::encoder::lpc::levinson(&mut lev, &rh, &rl, a, &mut rc, &mut ov);
    }
    let (mut rh, mut rl) = ([0i16; 11], [0i16; 11]);
    time("autocorr", |i| {
        bb(crate::encoder::lpc::autocorr(
            bb(&wins[i]),
            10,
            &mut rh,
            &mut rl,
            &crate::tables::WINDOW_200_40,
        ));
    });
    let mut lsp = [0i16; 10];
    time("az_lsp", |i| {
        crate::encoder::lpc::az_lsp(bb(&azs[i]), bb(&mut lsp), &crate::tables::LSP_INIT_DATA)
    });
    let gins: Vec<(GainIn, GainIn)> = (0..64)
        .map(|_| (gain_in(&mut g), gain_in(&mut g)))
        .collect();
    let mut code = [0i16; 40];
    (code[7], code[30]) = (8191, -8191);
    time("mr475_gain_quant", |i| {
        let (a, b) = &gins[i];
        let mut st = crate::gains::GcPredState::default();
        let (mut r, mut ov) = ([0i16; 4], false);
        let [p0, c0, p1, c1] = &mut r;
        bb(crate::encoder::gain_q::mr475_gain_quant(
            &mut st,
            a.exp_gcode0,
            a.frac_gcode0,
            &a.exp,
            &a.frac,
            10,
            20000,
            bb(&code),
            b.exp_gcode0,
            b.frac_gcode0,
            &b.exp,
            &b.frac,
            10,
            18000,
            a.gp_limit,
            p0,
            c0,
            p1,
            c1,
            &mut ov,
        ));
    });
    let cv = [3usize, 9, 21, 37];
    let sgn = [32767i16, -32768, 32767, -32768];
    time("filtered_code 4", |i| {
        filtered_code(bb(&hs[i]), bb(&cv), &sgn, bb(&mut y), &mut ov)
    });
    let mut t = tables();
    for mode in [
        CodecMode::MR59,
        CodecMode::MR67,
        CodecMode::MR74,
        CodecMode::MR102,
        CodecMode::MR122,
    ] {
        let cins: Vec<CbIn> = (0..64).map(|_| cb_in(&mut g, mode)).collect();
        time(&format!("cbsearch {mode:?}"), |i| {
            bb(cb_rust(bb(&cins[i])));
        });
        time(&format!("cbsearch {mode:?} C"), |i| {
            bb(cb_oc(bb(&cins[i]), &mut t));
        });
    }
    let gins: Vec<GainIn> = (0..64)
        .map(|_| GainIn {
            mode: CodecMode::MR102,
            ..gain_in(&mut g)
        })
        .collect();
    time("qua_gain 102", |i| {
        bb(qua_gain_rust(bb(&gins[i])));
    });
    time("qua_gain 102 C", |i| {
        bb(qua_gain_c(bb(&gins[i]), &mut t));
    });
    // MR795_gain_quant on encoder-like inputs (realistic exponents; the lanes accept).
    let gins: Vec<(GainIn, [[i16; 40]; 3])> = (0..64)
        .map(|_| {
            let x = GainIn {
                exp: std::array::from_fn(|_| 10 + (g.u64() % 12) as i16),
                exp_gcode0: (g.u64() % 16) as i16,
                ..gain_in(&mut g)
            };
            let mut v = [[0i16; 40]; 3];
            v.iter_mut().for_each(|s| signal(&mut g, s));
            (x, v)
        })
        .collect();
    time("mr795_gain_quant", |i| {
        let (x, v) = &gins[i];
        let (mut st, mut r, mut ana, mut pos) =
            (GainAdaptState::default(), [12000i16, 0, 0, 0], [0i16; 2], 0);
        let [gp, gc, q1, q2] = &mut r;
        mr795_gain_quant(
            &mut st,
            bb(&v[0]),
            &v[1],
            &v[2],
            bb(&x.frac),
            &x.exp,
            12,
            20000,
            x.exp_gcode0,
            x.frac_gcode0,
            20000,
            12,
            32767,
            gp,
            gc,
            q1,
            q2,
            &mut ana,
            &mut pos,
            &mut ov,
        );
        bb(r);
    });
    time("mr795_gain_quant C", |i| {
        let (x, v) = &gins[i];
        let (mut st, mut r, mut ana, mut o, mut v) = (
            GainAdaptState::default(),
            [12000i16, 0, 0, 0],
            [0i16; 2],
            0,
            *v,
        );
        let mut anap = ana.as_mut_ptr();
        let (mut f, mut e) = (x.frac, x.exp);
        let [gp, gc, q1, q2] = &mut r;
        let [a, b, d] = &mut v;
        unsafe {
            oc_mr795_gain_quant(
                &mut st,
                a.as_mut_ptr(),
                b.as_mut_ptr(),
                d.as_mut_ptr(),
                bb(f.as_mut_ptr()),
                e.as_mut_ptr(),
                12,
                20000,
                x.exp_gcode0,
                x.frac_gcode0,
                40,
                20000,
                12,
                32767,
                gp,
                gc,
                q1,
                q2,
                &mut anap,
                &mut t,
                &mut o,
            )
        };
        bb(r);
    });
    let mut signs = [[0i16; 40]; 64];
    for s in signs.iter_mut().flatten() {
        *s = if g.u64() & 1 == 0 { 32767 } else { -32767 };
    }
    let (mut rr, mut ov, mut o) = ([[0i16; 40]; 40], false, 0);
    let mut rt = RrTm::ZERO;
    time("cor_h_tm", |i| {
        cor_h_tm(bb(&hs[i]), bb(&signs[i]), bb(&mut rt), &mut ov)
    });
    // MR475/MR515's calls: blocks of different residues (sub_nr 1), of equal ones (sub_nr 0)
    time("cor_h_blocks", |i| {
        let b = [[[3, 0], [4, 2]], [[2, 0], [3, 1]]][i & 1];
        crate::simd::cor_h_blocks(bb(&hs[i]), bb(&signs[i]), bb(&mut rt), &b, &mut ov)
    });
    fn c_search<const T: usize, const P: usize>(
        (dn, rr, ip, pm): &SearchIn<T, P>,
        o: &mut CFlag,
    ) -> [i16; P] {
        let (mut dn, mut rr, mut ip, mut pm, mut cv) = (*dn, *rr, *ip, *pm, [0i16; P]);
        unsafe {
            c_search_10and8i40(
                P as i16,
                T as i16,
                T as i16,
                dn.as_mut_ptr(),
                rr.as_mut_ptr(),
                ip.as_mut_ptr(),
                pm.as_mut_ptr(),
                cv.as_mut_ptr(),
                o,
            )
        };
        cv
    }
    let sins: Vec<_> = (0..64)
        .map(|_| {
            let s = search_in::<5, 10>(&mut g, false);
            let t = tm(&s.1);
            (s, t)
        })
        .collect();
    time("search_10i40", |i| {
        let ((dn, _, ipos, pm), rt) = bb(&sins[i]);
        let (mut ip, mut cv) = (*ipos, [0usize; 10]);
        search_10i40(dn, rt, &mut ip, pm, &mut cv);
        bb(cv);
    });
    time("search_10i40 C", |i| {
        bb(c_search(&bb(&sins[i]).0, &mut o));
    });
    // The step-5 searches (MR67, MR74/MR795) on the same kind of rr, dn2 as set_sign marks
    // it (no C entry point of their own: `cbsearch MR67/MR74` time them with cor_h).
    let s5: Vec<_> = sins
        .iter()
        .map(|((dn, _, _, _), rt)| {
            let (mut d, mut sign, mut dn2) = (*dn, [0i16; 40], [0i16; 40]);
            crate::encoder::codebook::set_sign(&mut d, &mut sign, &mut dn2, 6);
            (d, dn2, rt)
        })
        .collect();
    time("search_3i40", |i| {
        let (dn, dn2, rt) = bb(&s5[i]);
        bb(crate::simd::search_3i40(dn, dn2, rt));
    });
    time("search_4i40", |i| {
        let (dn, dn2, rt) = bb(&s5[i]);
        bb(crate::simd::search_4i40(dn, dn2, rt));
    });
    // MR102's cor_h + search (one kernel) on the search inputs' h and signs.
    let sins: Vec<_> = (0..64)
        .map(|_| {
            let (h, sign) = (hs[g.u64() as usize % 64], signs[g.u64() as usize % 64]);
            let (dn, _, ipos, pm) = search_in::<4, 8>(&mut g, false);
            (h, sign, dn, ipos, pm)
        })
        .collect();
    time("cor_h_search_8i40", |i| {
        let (h, sign, dn, ipos, pm) = bb(&sins[i]);
        let (mut ip, mut cv) = (*ipos, [0usize; 8]);
        cor_h_search_8i40(h, sign, dn, &mut ip, pm, &mut cv, &mut ov);
        bb(cv);
    });
    time("cor_h_search_8i40 C", |i| {
        let (mut h, mut sign, dn, ipos, pm) = *bb(&sins[i]);
        let mut rr = [[0i16; 40]; 40];
        unsafe { oc_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), rr.as_mut_ptr(), &mut o) };
        bb(c_search(&(dn, rr, ipos, pm), &mut o));
    });
    time("cor_h C", |i| unsafe {
        oc_cor_h(
            bb(hs[i].as_ptr() as *mut i16),
            bb(signs[i].as_ptr() as *mut i16),
            bb(rr.as_mut_ptr()),
            &mut o,
        )
    });
}
