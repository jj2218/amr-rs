//! Kernel timings (`CRATE=amr-wb BENCH=kern just bench`, `KERN=<substring>`; with `AB=<ref>` the
//! min of 3 runs per side and B/A): ns per call of the encoder's dispatched kernels, on the
//! realistic classes of the diff tests' inputs only (decaying impulse responses, stable Q12
//! filters, moderate signals), so the fast paths are what is timed (amr-nb's harness once timed
//! exact fallbacks on full-scale inputs). Ranks kernels below the resolution of a whole-codec
//! A/B; stage timers (`BENCH=stages`) say which kernels matter.

use super::acelp2::vs_ref::search_input;
use super::filters::vs_ref::{lpc, signal};
use super::lpc::L_WINDOW;
use crate::cnst::{L_SUBFR, M};
use crate::simd::{
    acelp_2t64_fx, acelp_4t64_fx, autocorr, az_isp, convolve, cor_h_x, residu, syn_filt,
};
use crate::testutil::Rng;
use std::hint::black_box as bb;

/// ns per call of `f` over 64 rotating inputs (best of 20 runs of 5k calls), if its name
/// contains `$KERN` (or that is unset).
fn time(name: &str, mut f: impl FnMut(usize)) {
    if std::env::var("KERN").is_ok_and(|k| !name.contains(&k)) {
        return;
    }
    let best = (0..20)
        .map(|_| {
            let t = std::time::Instant::now();
            for i in 0..5_000 {
                f(i & 63);
            }
            t.elapsed().as_secs_f64() / 5_000.0 * 1e9
        })
        .fold(f64::INFINITY, f64::min);
    eprintln!("{name:<24} {best:8.1} ns");
}

#[test]
#[ignore]
fn kernel_timing() {
    let mut r = Rng::new(7);
    // codebook inputs: n not a multiple of 5 (wild h) or 7 (zero cn)
    let cases: Vec<_> = (0..)
        .filter(|n: &u64| !n.is_multiple_of(5) && !n.is_multiple_of(7))
        .take(64)
        .map(|n| search_input(&mut r, n))
        .collect();
    let dns: Vec<[i16; L_SUBFR]> = cases
        .iter()
        .map(|(h, x, _)| {
            let mut dn = [0i16; L_SUBFR];
            cor_h_x(h, x, &mut dn);
            dn
        })
        .collect();
    // stable filters (`lpc` draws a wild one 1 in 4) and moderate signals
    let filters: Vec<Vec<i16>> = std::iter::repeat_with(|| lpc(&mut r, M))
        .filter(|a| a[1..].iter().map(|&c| (c as i32).abs()).sum::<i32>() < 4 * 4096)
        .take(64)
        .collect();
    let sigs: Vec<Vec<i16>> = (0..64).map(|_| signal(&mut r, 2, L_WINDOW)).collect();

    let (mut code, mut y, mut dn, mut ix) =
        ([0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; 8]);
    for (nb, ser, label) in [
        (36, 253, "acelp_4t64 8.85"),
        (64, 365, "acelp_4t64 15.85"),
        (88, 461, "acelp_4t64 23.05"),
    ] {
        time(label, |i| {
            dn = dns[i];
            acelp_4t64_fx(
                bb(&mut dn),
                &cases[i].2,
                &cases[i].0,
                &mut code,
                &mut y,
                nb,
                ser,
                &mut ix,
            );
        });
    }
    let mut index = 0i16;
    time("acelp_2t64 6.60", |i| {
        dn = dns[i];
        acelp_2t64_fx(
            bb(&mut dn),
            &cases[i].2,
            &cases[i].0,
            &mut code,
            &mut y,
            &mut index,
        );
    });
    time("cor_h_x", |i| {
        cor_h_x(bb(&cases[i].0), &cases[i].1, &mut dn)
    });
    time("convolve", |i| {
        convolve(bb(&cases[i].1), &cases[i].0, &mut y, L_SUBFR)
    });
    time("residu", |i| {
        residu(bb(&filters[i]), &sigs[i], M, &mut y, L_SUBFR)
    });
    // the synthesis filter's input is an excitation, not a signal: 1/8 scale
    let quiet: Vec<[i16; L_SUBFR]> = sigs
        .iter()
        .map(|s| std::array::from_fn(|j| s[j] >> 3))
        .collect();
    let mut mem = [0i16; M];
    time("syn_filt", |i| {
        mem.fill(0);
        syn_filt(bb(&filters[i]), Some(&quiet[i]), &mut y, &mut mem, false);
    });
    let (mut r_h, mut r_l) = ([0i16; M + 1], [0i16; M + 1]);
    time("autocorr", |i| autocorr(bb(&sigs[i]), &mut r_h, &mut r_l));
    let (mut isp, old) = ([0i16; M], [0i16; M]);
    time("az_isp", |i| az_isp(bb(&filters[i]), &mut isp, &old));
}
