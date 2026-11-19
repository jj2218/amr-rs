//! Differential tests of `codebook` against the TS 26.073 reference (`ref26073-sys`;
//! helpers in `encoder::vs_ref`).
#![allow(unsafe_code, clippy::too_many_arguments)]

use crate::cnst::CodecMode;
use crate::cnst::L_SUBFR;
use crate::encoder::codebook::Tm;
use crate::encoder::codebook::cbsearch;
use crate::encoder::codebook::{Pair, chain};
use crate::encoder::codebook::{
    RrTm, RrTm4, cor_h_tm, cor_h_tm_portable, cor_h_x, filtered_code, search_8i40_portable,
    search_10i40, search_10i40_portable, set_sign, set_sign12k2,
};
use crate::encoder::vs_ref::CFlag;
use crate::encoder::vs_ref::refc::*;
use crate::encoder::vs_ref::tables;
use crate::encoder::vs_ref::*;
use crate::testutil::{N, Rng};

/// cbsearch for every mode: indices, code vector, filtered code vector, sharpened h.
#[test]
fn cbsearch_matches_c() {
    let mut g = Rng::new(8);
    let mut t = tables();
    for _ in 0..N / 20 {
        let mode = mode_of(&mut g);
        let c = cb_in(&mut g, mode);
        assert_eq!(
            cb_rust(&c),
            cb_c(&c, &mut t),
            "{mode:?} t0 {} sharp {}",
            c.t0,
            c.sharp
        );
    }
}

unsafe extern "C" {
    fn ref_cor_h_x(h: *mut i16, x: *mut i16, dn: *mut i16, sf: i16);
    fn ref_cor_h_x2(h: *mut i16, x: *mut i16, dn: *mut i16, sf: i16, nb: i16, step: i16);
    fn ref_set_sign12k2(
        dn: *mut i16,
        cn: *mut i16,
        sign: *mut i16,
        pos_max: *mut i16,
        nb_track: i16,
        ipos: *mut i16,
        step: i16,
    );
    fn ref_search_10and8i40(
        nb_pulse: i16,
        step: i16,
        nb_tracks: i16,
        dn: *mut i16,
        rr: *mut [i16; 40],
        ipos: *mut i16,
        pos_max: *mut i16,
        codvec: *mut i16,
    );
}

unsafe fn c_cor_h_x(h: *mut i16, x: *mut i16, dn: *mut i16, sf: i16, o: *mut CFlag) {
    unsafe { call(o, || ref_cor_h_x(h, x, dn, sf)) }
}

unsafe fn c_cor_h_x2(
    h: *mut i16,
    x: *mut i16,
    dn: *mut i16,
    sf: i16,
    nb: i16,
    step: i16,
    o: *mut CFlag,
) {
    unsafe { call(o, || ref_cor_h_x2(h, x, dn, sf, nb, step)) }
}

unsafe fn ref_search_10and8i40_w(
    nb_pulse: i16,
    step: i16,
    nb_tracks: i16,
    dn: *mut i16,
    rr: *mut [i16; 40],
    ipos: *mut i16,
    pos_max: *mut i16,
    codvec: *mut i16,
    o: *mut CFlag,
) {
    unsafe {
        call(o, || {
            ref_search_10and8i40(nb_pulse, step, nb_tracks, dn, rr, ipos, pos_max, codvec)
        })
    }
}

unsafe fn c_set_sign12k2(
    dn: *mut i16,
    cn: *mut i16,
    sign: *mut i16,
    pos_max: *mut i16,
    nb_track: i16,
    ipos: *mut i16,
    step: i16,
    o: *mut CFlag,
) {
    unsafe {
        call(o, || {
            ref_set_sign12k2(dn, cn, sign, pos_max, nb_track, ipos, step)
        })
    }
}

pub(crate) fn cb_rust(c: &CbIn) -> CbOut {
    let mut rh = c.hvec;
    let (mut rc, mut ry, mut ra) = ([0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; 10]);
    let (mut ov, mut pos) = (false, 0);
    cbsearch(
        &c.x,
        &mut rh[L_SUBFR..],
        c.t0,
        c.sharp,
        c.gain_pit,
        &c.res2,
        &mut rc,
        &mut ry,
        &mut ra,
        &mut pos,
        c.mode,
        c.sub_nr,
        &mut ov,
    );
    (pos, ra, rc, ry, rh)
}

/// `rr` in a track-major layout (padding lanes zero).
pub(crate) fn tm<const T: usize, const N: usize, const W: usize>(
    rr: &[[i16; 40]; 40],
) -> Tm<T, N, W> {
    let mut t = Tm::ZERO;
    for i in 0..40 {
        for j in 0..40 {
            t.b[i % T][j % T][i / T][j / T] = rr[i][j];
        }
        t.diag[i % T][i / T] = rr[i][i];
    }
    t
}

/// `search_10i40` (MR122) and `search_8i40_portable` (MR102, both widths) vs the C's
/// `search_10and8i40`, each with both `Pair::best` forms (`L32`).
fn search_matches_c<const T: usize, const P: usize>(
    seed: u64,
    rust: impl Fn(&[i16; 40], &[[i16; 40]; 40], &mut [i16; P], &[i16; T], &mut [usize; P]),
) {
    let mut g = Rng::new(seed);
    for _ in 0..N / 4 {
        let (mut dn, mut rr, mut ipos, mut pos_max) = search_in::<T, P>(&mut g, true);
        let (mut rpos, mut cv) = (ipos, [0usize; P]);
        rust(&dn, &rr, &mut rpos, &pos_max, &mut cv);
        let (mut ccv, mut o) = ([0i16; P], 0);
        unsafe {
            ref_search_10and8i40_w(
                P as i16,
                T as i16,
                T as i16,
                dn.as_mut_ptr(),
                rr.as_mut_ptr(),
                ipos.as_mut_ptr(),
                pos_max.as_mut_ptr(),
                ccv.as_mut_ptr(),
                &mut o,
            )
        };
        assert_eq!(
            (cv.map(|v| v as i16), rpos),
            (ccv, ipos),
            "{dn:?} {pos_max:?}"
        );
    }
}

/// MR102's `cor_h` (both instances of `cor_h_search_8i40`: 10 lanes and the ±32767 form,
/// 16 lanes and the multiply form).
#[test]
fn cor_h_tm4_matches_c() {
    let mut g = Rng::new(21);
    for _ in 0..N {
        let mut h = [0i16; 40];
        impulse(&mut g, &mut h);
        let mut sign = [0i16; 40];
        signs(&mut g, &mut sign);
        let (mut crr, mut o) = ([[0i16; 40]; 40], 0);
        unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), crr.as_mut_ptr(), &mut o) };
        let (mut r10, mut r16, mut ov10, mut ov16) = (RrTm4::ZERO, RrTm4::ZERO, false, false);
        cor_h_tm_portable::<4, 10, 10, true>(&h, &sign, &mut r10, &mut ov10);
        cor_h_tm_portable::<4, 10, 16, false>(&h, &sign, &mut r16, &mut ov16);
        let (c10, c16) = (tm::<4, 10, 10>(&crr), tm::<4, 10, 16>(&crr));
        assert_eq!(
            (r10.b, r10.diag, ov10 && o != 0),
            (c10.b, c10.diag, ov10 && o != 0),
            "{h:?} {sign:?}"
        );
        assert_eq!(
            (r16.b, r16.diag, ov16 && o != 0),
            (c16.b, c16.diag, ov16 && o != 0),
            "{h:?} {sign:?}"
        );
    }
}

#[test]
fn cor_h_tm_matches_c() {
    let mut g = Rng::new(26);
    for _ in 0..N / 4 {
        let mut h = [0i16; 40];
        impulse(&mut g, &mut h);
        let mut sign = [0i16; 40];
        signs(&mut g, &mut sign);
        let (mut rr, mut ov) = (RrTm::ZERO, false);
        cor_h_tm(&h, &sign, &mut rr, &mut ov);
        let (mut crr, mut o) = ([[0i16; 40]; 40], 0);
        unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), crr.as_mut_ptr(), &mut o) };
        let c = tm(&crr);
        let _ = (ov, o);
        assert_eq!((rr.b, rr.diag), (c.b, c.diag), "{h:?} {sign:?}");
    }
}

/// `cor_h_blocks` (dispatched) vs the C's `cor_h`: two random blocks (any tracks, equal ones
/// too; half the cases with equal residues, one shared accumulation), with every diagonal;
/// the other blocks stay untouched (a 0x5555 fill). Signs ±32767 (as `set_sign`'s) and of
/// any value.
#[test]
fn cor_h_blocks_matches_c() {
    let mut g = Rng::new(27);
    for _ in 0..N / 4 {
        let mut h = [0i16; 40];
        impulse(&mut g, &mut h);
        let mut sign = [0i16; 40];
        signs(&mut g, &mut sign);
        let mut blocks: [[usize; 2]; 2] =
            core::array::from_fn(|_| [(g.u64() % 5) as usize, (g.u64() % 5) as usize]);
        if g.u64() & 1 == 0 {
            // equal residues (one shared accumulation), as in half the codec's calls
            let r = (blocks[0][1] + 5 - blocks[0][0]) % 5;
            blocks[1][1] = (blocks[1][0] + r) % 5;
        }
        let mut rr = RrTm::ZERO;
        rr.b = [[[[0x5555; 8]; 8]; 5]; 5];
        let mut ov = false;
        crate::simd::cor_h_blocks(&h, &sign, &mut rr, &blocks, &mut ov);
        let (mut crr, mut o) = ([[0i16; 40]; 40], 0);
        unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), crr.as_mut_ptr(), &mut o) };
        let c = tm(&crr);
        let _ = (ov, o);
        assert_eq!(rr.diag, c.diag, "{h:?} {sign:?}");
        for ta in 0..5 {
            for tb in 0..5 {
                let want = if blocks.contains(&[ta, tb]) || ta == tb {
                    c.b[ta][tb]
                } else {
                    [[0x5555; 8]; 8]
                };
                if blocks.contains(&[ta, tb]) || ta != tb {
                    assert_eq!(rr.b[ta][tb], want, "block {ta} {tb}: {h:?} {sign:?}");
                } else {
                    // a diagonal block not asked for: only its diagonal is written
                    for n in 0..8 {
                        assert_eq!(rr.b[ta][ta][n][n], c.b[ta][ta][n][n]);
                    }
                }
            }
        }
    }
}

/// `pm_scale` is cor_h's sign scaling for ±32767 signs, for every `tmp1`.
#[test]
fn pm_scale_is_exact() {
    use crate::encoder::codebook::pm_scale;
    for t in i16::MIN..=i16::MAX {
        for (si, sj) in [
            (32767, 32767),
            (-32767, -32767),
            (32767, -32767),
            (-32767, 32767),
        ] {
            let tmp2 = ((si * sj) >> 15) as i16;
            let want = ((t as i32 * tmp2 as i32) >> 15) as i16;
            assert_eq!(pm_scale(t, ((si ^ sj) >> 15) as i16), want, "{t} {si} {sj}");
        }
    }
}

#[test]
fn search_10i40_matches_c() {
    search_matches_c::<5, 10>(27, |dn, rr, ip, pm, cv| {
        search_10i40(dn, &tm(rr), ip, pm, cv)
    });
    search_matches_c::<5, 10>(27, |dn, rr, ip, pm, cv| {
        search_10i40_portable::<false>(dn, &tm(rr), ip, pm, cv)
    });
    search_matches_c::<5, 10>(27, |dn, rr, ip, pm, cv| {
        search_10i40_portable::<true>(dn, &tm(rr), ip, pm, cv)
    });
}

#[test]
fn search_8i40_matches_c() {
    search_matches_c::<4, 8>(29, |dn, rr, ip, pm, cv| {
        search_8i40_portable::<10, false>(dn, &tm::<4, 10, 10>(rr), ip, pm, cv)
    });
    search_matches_c::<4, 8>(29, |dn, rr, ip, pm, cv| {
        search_8i40_portable::<16, true>(dn, &tm::<4, 10, 16>(rr), ip, pm, cv)
    });
}

/// The step-5 searches' `step5_block` (8 first pulses in lanes, later pulses' rows
/// transposed, plain wrapping energy sums: `STEP5_SHIFTS`' proof) is the C's nested loops
/// (`for i = start; i < L_CODE; i += STEP` over the square `rr`) with TS 26.073's saturating
/// operators, for every lane: random symmetric `rr` of any magnitude, a ±full-scale class
/// (the chains at their proven limits, 4 pulses), `dn` near the i16 limits, and a tie-heavy
/// class (few distinct values).
#[test]
fn step5_block_is_the_c_loops() {
    use crate::basic_op::{add, l_mac, l_msu, l_mult, mult, round};
    use crate::encoder::codebook::{dn_tm, step5_block};
    let mut g = Rng::new(29);
    for case in 0..N {
        let ties = case % 2 == 0;
        // the 16-bit extremes (np = 4, not `ties`): the chains at their range limits
        let full = case % 6 == 5;
        let shift = (g.u64() % 16) as u32;
        let mut rr = [[0i16; 40]; 40];
        for i in 0..40 {
            for j in 0..=i {
                rr[i][j] = if ties {
                    [0, 1, 4096, -4096][(g.u64() % 4) as usize]
                } else if full {
                    [32767, 32766, -32768, -32767][(g.u64() % 4) as usize]
                } else {
                    g.i16() >> shift
                };
                rr[j][i] = rr[i][j];
            }
        }
        let mut dn = [0i16; 40];
        for d in &mut dn {
            *d = if ties {
                [0, 100, -100, 16384][(g.u64() % 4) as usize]
            } else if g.u64().is_multiple_of(3) {
                [32767, -32768, 16384][(g.u64() % 3) as usize]
            } else {
                g.i16()
            };
        }
        let ipos: [usize; 4] = core::array::from_fn(|_| (g.u64() % 5) as usize);
        // The C's loop for one pulse: weights 1/4, 1/2 (pulse 1) or 1/16, 1/8, the cross
        // terms from the latest pulse back.
        let c_loop = |ps0: i16, alp0: i32, t: usize, prev: &[usize], p: usize, o: &mut bool| {
            let (wd, wc) = if p == 1 { (8192, 16384) } else { (2048, 4096) };
            let (mut sq, mut alp, mut ps, mut ix) = (-1i16, 1i16, 0i16, t);
            for i in (t..40).step_by(5) {
                let ps1 = add(ps0, dn[i], o);
                let mut alp1 = l_mac(alp0, rr[i][i], wd, o);
                for &q in prev.iter().rev() {
                    alp1 = l_mac(alp1, rr[q][i], wc, o);
                }
                let sq1 = mult(ps1, ps1, o);
                let alp_16 = round(alp1, o);
                if l_msu(l_mult(alp, sq1, o), sq, alp_16, o) > 0 {
                    (sq, ps, alp, ix) = (sq1, ps1, alp_16, i);
                }
            }
            (sq, ps, alp, ix)
        };
        let dt = dn_tm(&dn);
        let rt: RrTm = tm(&rr);
        let np = 2 + case % 3;
        let mut o = false;
        let (got, gix) = match np {
            2 => {
                let (l, ix) = step5_block(&[ipos[0], ipos[1]], &dt, &rt);
                (l, ix.to_vec())
            }
            3 => {
                let (l, ix) = step5_block(&[ipos[0], ipos[1], ipos[2]], &dt, &rt);
                (l, ix.to_vec())
            }
            _ => {
                let (l, ix) = step5_block(&ipos, &dt, &rt);
                (l, ix.to_vec())
            }
        };
        for l in 0..8 {
            let i0 = ipos[0] + 5 * l;
            let mut pulses = vec![i0];
            let (mut ps, mut alp0) = (dn[i0], (rr[i0][i0] as i32) << 14);
            let mut last = (0, 0, 0, 0);
            for p in 1..np {
                last = c_loop(ps, alp0, ipos[p], &pulses, p, &mut o);
                pulses.push(last.3);
                ps = last.1;
                alp0 = (last.2 as i32) << if p == 1 { 14 } else { 16 };
            }
            let lane: Vec<usize> = (0..np).map(|p| ipos[p] + 5 * gix[p][l] as usize).collect();
            assert_eq!(
                (got.sq[l], got.ps[l], got.alp[l], lane),
                (last.0, last.1, last.2, pulses),
                "case {case} lane {l} {ipos:?}"
            );
        }
    }
}

/// `simd::transpose8` (SSE2 intrinsics on x86_64) is the transpose, and so is the array form.
#[test]
fn transpose8_is_the_transpose() {
    use crate::encoder::codebook::transpose8_portable;
    let mut g = Rng::new(32);
    for _ in 0..1000 {
        let r: [[i16; 8]; 8] = core::array::from_fn(|_| core::array::from_fn(|_| g.i16()));
        let want: [[i16; 8]; 8] = core::array::from_fn(|n| core::array::from_fn(|l| r[l][n]));
        assert_eq!(crate::simd::transpose8(r), want);
        assert_eq!(transpose8_portable(r), want);
    }
}

/// `search_4i40` (MR74/MR795; packed half blocks, `step5_halves`) is the C's four nested
/// loops with the `dn2 >= 0` skip and the outer compare chain, both instances (dispatched
/// and portable). Classes: random `rr`/`dn`, tie-heavy (few distinct values), and tracks
/// with 5..=8 cells of |dn| = 32767, which `set_sign` (its stale `pos`) leaves with 5..=8
/// valid first pulses: blocks split into a full and a partial half, pairs then mix blocks.
#[test]
fn search_4i40_is_the_c_loops() {
    use crate::encoder::codebook::search_4i40_portable;
    let mut g = Rng::new(31);
    let (mut wide, mut ties_n) = (0, 0);
    for case in 0..N {
        let ties = case % 3 == 0;
        let shift = (g.u64() % 16) as u32;
        let mut rr = [[0i16; 40]; 40];
        for i in 0..40 {
            for j in 0..=i {
                rr[i][j] = if ties {
                    [0, 1, 4096, -4096][(g.u64() % 4) as usize]
                } else {
                    g.i16() >> shift
                };
                rr[j][i] = rr[i][j];
            }
        }
        let mut dn = [0i16; 40];
        for d in &mut dn {
            *d = if ties {
                [0, 100, -100, 16384][(g.u64() % 4) as usize]
            } else {
                g.i16() >> (g.u64() % 4)
            };
        }
        if case % 3 == 1 {
            // A track with k > 4 cells of |dn| = 32767 keeps k first pulses (only 8 - k
            // removals find a minimum below 32767).
            for t in 0..5 {
                let k = 5 + (g.u64() % 4) as usize;
                for n in 0..8 {
                    if n < k {
                        dn[t + 5 * n] = [32767, -32767, -32768][(g.u64() % 3) as usize];
                    } else if dn[t + 5 * n].unsigned_abs() >= 32767 {
                        dn[t + 5 * n] = 0;
                    }
                }
            }
        }
        let (mut sign, mut dn2) = ([0i16; 40], [0i16; 40]);
        set_sign(&mut dn, &mut sign, &mut dn2, 4);
        let valid: Vec<usize> = (0..5)
            .map(|t| (0..8).filter(|n| dn2[t + 5 * n] >= 0).count())
            .collect();
        wide += valid.iter().any(|&v| v > 4) as usize;
        ties_n += ties as usize;
        // The C's loops with TS 26.073's saturating operators.
        use crate::basic_op::{add, l_mac, l_msu, l_mult, mult, round};
        let o = &mut false;
        let mut want = [0usize, 1, 2, 3];
        let (mut psk, mut alpk) = (-1i16, 1i16);
        for track in 3..5 {
            let mut ipos = [0usize, 1, 2, track];
            for _ in 0..4 {
                for i0 in (ipos[0]..40).step_by(5) {
                    if dn2[i0] < 0 {
                        continue;
                    }
                    let mut pulses = vec![i0];
                    let (mut ps0, mut alp0) = (dn[i0], l_mult(rr[i0][i0], 8192, o));
                    let (mut sq, mut alp) = (-1i16, 1i16);
                    for p in 1..4 {
                        let (wd, wc) = if p == 1 { (8192, 16384) } else { (2048, 4096) };
                        let mut ps = 0i16;
                        let mut ix = ipos[p];
                        (sq, alp) = (-1, 1);
                        for i in (ipos[p]..40).step_by(5) {
                            let ps1 = add(ps0, dn[i], o);
                            let mut alp1 = l_mac(alp0, rr[i][i], wd, o);
                            for &q in pulses.iter().rev() {
                                alp1 = l_mac(alp1, rr[q][i], wc, o);
                            }
                            let sq1 = mult(ps1, ps1, o);
                            let alp_16 = round(alp1, o);
                            if l_msu(l_mult(alp, sq1, o), sq, alp_16, o) > 0 {
                                (sq, ps, alp, ix) = (sq1, ps1, alp_16, i);
                            }
                        }
                        pulses.push(ix);
                        ps0 = ps;
                        alp0 = if p == 1 {
                            l_mult(alp, 8192, o)
                        } else {
                            (alp as i32) << 16
                        };
                    }
                    if l_msu(l_mult(alpk, sq, o), psk, alp, o) > 0 {
                        (psk, alpk) = (sq, alp);
                        want = [pulses[0], pulses[1], pulses[2], pulses[3]];
                    }
                }
                ipos.rotate_right(1);
            }
        }
        let rt: RrTm = tm(&rr);
        let got = crate::simd::search_4i40(&dn, &dn2, &rt);
        assert_eq!(got, want, "case {case} valid {valid:?}");
        let got = search_4i40_portable(&dn, &dn2, &rt);
        assert_eq!(got, want, "case {case} portable, valid {valid:?}");
    }
    assert!(
        wide >= N / 4 && ties_n >= N / 4,
        "classes: wide {wide}, ties {ties_n}"
    );
}

/// `lanes_next` is the C's loop (TS 26.073 operators) for every lane, on its callers'
/// domain: per pulse level the start energy as they make it (pulse 1: rr << 14; pulse 2:
/// alp << 14, alp any i16; pulse 3: alp << 16, |alp| <= 18432) with 1, 2 or 3 cross terms
/// of any 16-bit value, `d` and `diag` any 16-bit rows (two half blocks' rows, as
/// `step5_halves` passes them). The `const` assertion after `lanes_next` proves no chain
/// leaves the i32 range there; the C side here sums in i64 and asserts that too.
#[test]
fn lanes_next_is_the_c_loop() {
    use crate::basic_op::{add, mult, round};
    use crate::encoder::codebook::lanes_next;
    let mut g = Rng::new(30);
    for case in 0..N {
        let level = case % 3;
        let full = case % 5 == 0; // the 16-bit extremes only
        let v16 = |g: &mut Rng| {
            if full {
                [i16::MIN, i16::MAX][(g.u64() % 2) as usize]
            } else {
                g.i16() >> (g.u64() % 16)
            }
        };
        let ps0: [i16; 8] = core::array::from_fn(|_| g.i16());
        let alp0: [i32; 8] = core::array::from_fn(|_| match level {
            0 | 1 => (v16(&mut g) as i32) << 14,
            _ => ((v16(&mut g) as i32).clamp(-18432, 18431)) << 16,
        });
        let (sd, sc) = if level == 0 { (14, 15) } else { (12, 13) };
        let dh: [[i16; 8]; 2] = core::array::from_fn(|_| core::array::from_fn(|_| g.i16()));
        let gh: [[i16; 8]; 2] = core::array::from_fn(|_| core::array::from_fn(|_| v16(&mut g)));
        let (d, diag) = ([&dh[0], &dh[1]], [&gh[0], &gh[1]]);
        let x: [[i32; 8]; 8] = core::array::from_fn(|_| {
            core::array::from_fn(|_| (0..=level).map(|_| v16(&mut g) as i32).sum())
        });
        let got = lanes_next(&ps0, &alp0, &d, &diag, &x, sd, sc);
        let mut o = false;
        for l in 0..8 {
            let (mut sq, mut alp, mut ps, mut ix) = (-1i16, 1i16, 0i16, 0i16);
            for n in 0..8 {
                let ps1 = add(ps0[l], dh[l / 4][n], &mut o);
                let alp1 =
                    alp0[l] as i64 + ((gh[l / 4][n] as i64) << sd) + ((x[n][l] as i64) << sc);
                assert!(i32::try_from(alp1).is_ok(), "level {level}: {alp1}");
                let sq1 = mult(ps1, ps1, &mut o);
                let alp_16 = round(alp1 as i32, &mut o);
                if alp as i32 * sq1 as i32 > sq as i32 * alp_16 as i32 {
                    (sq, ps, alp, ix) = (sq1, ps1, alp_16, n as i16);
                }
            }
            assert_eq!(
                (got.sq[l], got.ps[l], got.alp[l], got.ix[l]),
                (sq, ps, alp, ix),
                "case {case} lane {l}"
            );
        }
    }
}

/// MR122's 8 × 8 candidates, MR102's 10 × 10 in 16 lanes (AVX2) and 10 (portable).
#[test]
fn best_is_the_chain_8_10() {
    best_is_the_chain::<8, 8>(28);
    best_is_the_chain::<10, 16>(30);
    best_is_the_chain::<10, 10>(32);
}

/// `cor_h_x` (5 tracks) and `cor_h_x2` (MR102's 4) vs the C, both `sf`; full-scale classes
/// (`signal` edge values, `impulse`'s full-scale `h`) make the sums wrap and `l_shl`
/// saturate.
#[test]
fn cor_h_x_matches_c() {
    let mut g = Rng::new(31);
    for case in 0..N / 2 {
        let (mut x, mut h) = ([0i16; 40], [0i16; 40]);
        signal(&mut g, &mut x);
        impulse(&mut g, &mut h);
        if case % 7 == 0 {
            x.iter_mut()
                .for_each(|v| *v = if g.u64() & 1 == 0 { 32767 } else { -32768 });
        }
        let sf = 1 + (case & 1) as i16;
        let (mut dn, mut ov) = ([0i16; 40], false);
        let (mut cdn, mut co) = ([0i16; 40], 0 as CFlag);
        if case % 3 == 0 {
            crate::encoder::codebook::cor_h_x2(&h, &x, &mut dn, sf, 4, 4, &mut ov);
            unsafe {
                c_cor_h_x2(
                    h.as_mut_ptr(),
                    x.as_mut_ptr(),
                    cdn.as_mut_ptr(),
                    sf,
                    4,
                    4,
                    &mut co,
                )
            };
        } else {
            cor_h_x(&h, &x, &mut dn, sf, &mut ov);
            unsafe {
                c_cor_h_x(
                    h.as_mut_ptr(),
                    x.as_mut_ptr(),
                    cdn.as_mut_ptr(),
                    sf,
                    &mut co,
                )
            };
        }
        assert_eq!((dn, ov), (cdn, co != 0), "case {case}: {x:?} {h:?}");
    }
}

/// `sign_cor` (`set_sign12k2`'s lanes) vs the operator chain for any gains >= 0: the clamps
/// at both ends, which `set_sign12k2`'s normalised gains don't reach.
#[test]
fn sign_cor_is_the_operator_chain() {
    use crate::basic_op::{l_mac, l_mult, l_shl, round};
    const L_CODE: usize = 40;
    let mut g = Rng::new(41);
    let mut clamped = 0;
    for _ in 0..N / 10 {
        let (kc, kd) = (g.i16().max(0), g.i16().max(0));
        let cn: [i16; L_CODE] = core::array::from_fn(|_| g.i16());
        let dn: [i16; L_CODE] = core::array::from_fn(|_| g.i16());
        let o = &mut false;
        let want: [i16; L_CODE] = core::array::from_fn(|i| {
            round(l_shl(l_mac(l_mult(kc, cn[i], o), kd, dn[i], o), 10, o), o)
        });
        assert_eq!(super::super::codebook::sign_cor(kc, &cn, kd, &dn), want);
        clamped += want
            .iter()
            .filter(|&&v| v == i16::MAX || v == i16::MIN)
            .count();
    }
    assert!(clamped > 1000, "clamped: {clamped}");
}

/// `set_sign12k2` (its per-sample pass as selects, the `cn` energy as an exact sum) vs the
/// C, both track layouts: speech-like `dn`/`cn`; full-scale ones (-32768, the energy sum
/// and `L_mac` saturate, the flag); constant ones (equal `en`: ties in the track maxima).
#[test]
fn set_sign12k2_matches_c() {
    let mut g = Rng::new(30);
    let (mut sat, mut ties, mut lmac) = (0, 0, 0);
    for n in 0..N {
        let (mut dn, mut cn) = ([0i16; 40], [0i16; 40]);
        if n % 8 == 3 {
            // L_mac saturating (only reachable this way: |k_cn * cn << 1| < 2^20, so it
            // needs k_dn * dn ~ 2^30, and the energy is even): sum dn^2 = 2^31 wraps the
            // energy to 256, k_dn = Inv_sqrt(256) >> 11 = 32767, at dn = 32767 that is
            // 2 * k_dn * dn = 2^31 - 131070; cn = 3..10 alone keeps k_cn near 32767, so
            // k_cn * cn << 1 > 131070 tips the sum over.
            let mut rest = 1i64 << 31;
            for v in dn.iter_mut() {
                let r = rest.isqrt().min(32767);
                *v = r as i16;
                rest -= r * r;
            }
            assert_eq!((rest, dn[0], dn[1]), (0, 32767, 32767));
            cn[(g.u64() % 2) as usize] = 3 + (g.u64() % 8) as i16;
        } else {
            match n % 4 {
                0 => {
                    const V: [i16; 4] = [32767, -32768, 30000, -30000];
                    for v in dn.iter_mut().chain(&mut cn) {
                        *v = if g.u64() & 1 == 0 {
                            V[(g.u64() % 4) as usize]
                        } else {
                            g.i16()
                        };
                    }
                }
                1 => {
                    (dn, cn) = (
                        [(g.u64() % 64) as i16 - 32; 40],
                        [(g.u64() % 64) as i16; 40],
                    );
                    ties += 1;
                }
                _ => {
                    signal(&mut g, &mut dn);
                    signal(&mut g, &mut cn);
                }
            }
        }
        let (nb, step) = if n & 4 == 0 { (5, 5) } else { (4, 4) };
        let (mut rd, mut rs, mut rp, mut ri, mut ov) =
            (dn, [0i16; 40], [0i16; 5], [0i16; 10], false);
        set_sign12k2(&mut rd, &cn, &mut rs, &mut rp, nb, &mut ri, step, &mut ov);
        let (mut cd, mut cc, mut cs, mut cp, mut ci, mut o) =
            (dn, cn, [0i16; 40], [0i16; 5], [0i16; 10], 0);
        unsafe {
            c_set_sign12k2(
                cd.as_mut_ptr(),
                cc.as_mut_ptr(),
                cs.as_mut_ptr(),
                cp.as_mut_ptr(),
                nb as i16,
                ci.as_mut_ptr(),
                step as i16,
                &mut o,
            )
        };
        sat += (o != 0) as usize;
        // In that class only L_mac can saturate (energies, round and negate don't).
        lmac += (n % 8 == 3 && o != 0) as usize;
        assert_eq!(
            (rd, rs, rp, ri, ov),
            (cd, cs, cp, ci, o != 0),
            "{dn:?} {cn:?}"
        );
    }
    assert!(
        sat > N / 20 && ties > N / 8 && lmac > N / 10,
        "classes: {sat} {ties} {lmac}"
    );
}

/// `filtered_code` in lanes vs the C's chain step by step (`L_mult`, `L_mac`, `round`: the
/// ported basic ops, each checked against the C): full-scale `h` with -32768 taps, signs
/// +-32767/-32768 or small sets, pulses stacked on few positions (saturation, and `L_mac`'s
/// "accumulator equals the product" quirk); the overflow flag too.
#[test]
fn filtered_code_matches_chain() {
    use crate::basic_op::{l_mac, l_mult, round};
    let mut g = Rng::new(28);
    let (mut n_ov, mut n_quirk) = (0, 0);
    const V: [i16; 6] = [32767, -32768, 16384, -16384, 21846, -21846];
    for n in 0..N {
        let class = n % 3;
        let mut h = [0i16; 40];
        for v in &mut h {
            *v = match class {
                0 => g.u64() as i16 >> (g.u64() % 8),
                _ => V[g.u64() as usize % 6],
            };
        }
        let np = 1 + (g.u64() % 10) as usize;
        let (mut cv, mut sg) = ([0usize; 10], [0i16; 10]);
        for k in 0..np {
            cv[k] = (g.u64() % if class == 2 { 3 } else { 40 }) as usize;
            sg[k] = match class {
                0 => [32767, -32768][g.u64() as usize % 2],
                _ => V[g.u64() as usize % 6],
            };
        }
        let (mut y, mut ov) = ([0i16; 40], false);
        filtered_code(&h, &cv[..np], &sg[..np], &mut y, &mut ov);
        let (mut yr, mut ovr) = ([0i16; 40], false);
        for i in 0..40 {
            let at = |c: usize| if i >= c { h[i - c] } else { 0 };
            let mut s = l_mult(at(cv[0]), sg[0], &mut ovr);
            for k in 1..np {
                n_quirk += (s == at(cv[k]) as i32 * sg[k] as i32 && s != 0) as usize;
                s = l_mac(s, at(cv[k]), sg[k], &mut ovr);
            }
            yr[i] = round(s, &mut ovr);
        }
        n_ov += ovr as usize;
        assert_eq!((y, ov), (yr, ovr), "{h:?} {cv:?} {sg:?}");
    }
    // n_quirk: L_var3 == product, where opencore skipped the saturation
    assert!(
        n_ov > N / 10 && n_quirk > 100,
        "classes too rare: {n_ov} {n_quirk}"
    );
}

/// `cbsearch` inputs.
pub(crate) struct CbIn {
    pub(crate) mode: CodecMode,
    pub(crate) x: [i16; L_SUBFR],
    pub(crate) res2: [i16; L_SUBFR],
    pub(crate) hvec: [i16; 2 * L_SUBFR],
    pub(crate) t0: i16,
    pub(crate) sharp: i16,
    pub(crate) gain_pit: i16,
    pub(crate) sub_nr: usize,
}

/// (index words written, indices, code, filtered code, sharpened h)
type CbOut = (
    usize,
    [i16; 10],
    [i16; L_SUBFR],
    [i16; L_SUBFR],
    [i16; 2 * L_SUBFR],
);

pub(crate) fn cb_in(g: &mut Rng, mode: CodecMode) -> CbIn {
    let (mut x, mut res2) = ([0i16; L_SUBFR], [0i16; L_SUBFR]);
    signal(g, &mut x);
    signal(g, &mut res2);
    let mut hvec = [0i16; 2 * L_SUBFR];
    // An impulse response: decaying, first tap large.
    let scale = (g.u64() % 4) as u32;
    for (i, v) in hvec[L_SUBFR..].iter_mut().enumerate() {
        *v = (((g.u64() as i16) >> (scale + (i as u32 / 6).min(12))) as i32).clamp(-32768, 32767)
            as i16;
    }
    hvec[L_SUBFR] = 4096 + (g.u64() % 8192) as i16;
    CbIn {
        mode,
        x,
        res2,
        hvec,
        t0: 18 + (g.u64() % 126) as i16,
        sharp: (g.u64() % 13018) as i16,
        gain_pit: (g.u64() % 19662) as i16,
        sub_nr: (g.u64() % 4) as usize,
    }
}

fn cb_c(c: &CbIn, t: &mut CommonAmrTbls) -> CbOut {
    let (mut x, mut res2, mut ch) = (c.x, c.res2, c.hvec);
    let (mut cc, mut cy, mut ca) = ([0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; 10]);
    let mut anap = ca.as_mut_ptr();
    let mut o = 0;
    unsafe {
        c_cbsearch(
            x.as_mut_ptr(),
            ch.as_mut_ptr().add(L_SUBFR),
            c.t0,
            c.sharp,
            c.gain_pit,
            res2.as_mut_ptr(),
            cc.as_mut_ptr(),
            cy.as_mut_ptr(),
            &mut anap,
            c.mode as i32,
            c.sub_nr as i16,
            t,
            &mut o,
        )
    };
    let n = unsafe { anap.offset_from(ca.as_ptr()) } as usize;
    (n, ca, cc, cy, ch)
}

/// `cb_c` on opencore (for `kernel_timing`).
pub(crate) fn cb_oc(c: &CbIn, t: &mut CommonAmrTbls) -> CbOut {
    let (mut x, mut res2, mut ch) = (c.x, c.res2, c.hvec);
    let (mut cc, mut cy, mut ca) = ([0i16; L_SUBFR], [0i16; L_SUBFR], [0i16; 10]);
    let mut anap = ca.as_mut_ptr();
    let mut o = 0;
    unsafe {
        oc_cbsearch(
            x.as_mut_ptr(),
            ch.as_mut_ptr().add(L_SUBFR),
            c.t0,
            c.sharp,
            c.gain_pit,
            res2.as_mut_ptr(),
            cc.as_mut_ptr(),
            cy.as_mut_ptr(),
            &mut anap,
            c.mode as i32,
            c.sub_nr as i16,
            t,
            &mut o,
        )
    };
    let n = unsafe { anap.offset_from(ca.as_ptr()) } as usize;
    (n, ca, cc, cy, ch)
}

/// `search_10and8i40`'s inputs: `dn`, square `rr`, `ipos`, `pos_max`.
pub(crate) type SearchIn<const T: usize, const P: usize> =
    ([i16; 40], [[i16; 40]; 40], [i16; P], [i16; T]);

/// Search inputs for `T` tracks (MR122: 5, MR102: 4): `rr` from `cor_h` with the encoder's
/// kind of `h` and signs (the timing class: `Pair::best` falls back on <1% of picks, 0 on the
/// vectors) or, if `odd`, of any `h` and signs (non-positive energies), sometimes a random
/// symmetric matrix (energies of any sign: the exact-chain fallback), or a slowly decaying `h`
/// with one sign for all pulses and small `dn` (large positive `rr` everywhere: the pair
/// chains saturate while `dn` keeps the fast path, so `pair_safe` must decline), and `dn` from
/// a few values whose sums tie and wrap to -32768; tracks and maxima as `set_sign12k2` leaves them.
pub(crate) fn search_in<const T: usize, const P: usize>(g: &mut Rng, odd: bool) -> SearchIn<T, P> {
    let mut rr = [[0i16; 40]; 40];
    let flat = odd && g.u64().is_multiple_of(4);
    if flat {
        let mut h = [0i16; 40];
        let (top, step) = (8192 + (g.u64() % 24576) as i16, (g.u64() % 400) as i16);
        for (i, v) in h.iter_mut().enumerate() {
            *v = top - step * i as i16 / 2;
        }
        let mut sign = [if g.u64() & 1 == 0 { 32767 } else { -32767 }; 40];
        let mut o = 0;
        unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), rr.as_mut_ptr(), &mut o) };
    } else if odd && g.u64().is_multiple_of(4) {
        let shift = (g.u64() % 16) as u32;
        for i in 0..40 {
            for j in 0..=i {
                rr[i][j] = g.i16() >> shift;
                rr[j][i] = rr[i][j];
            }
        }
    } else {
        let mut h = [0i16; 40];
        let mut sign = [0i16; 40];
        if odd {
            impulse(g, &mut h);
            signs(g, &mut sign);
        } else {
            // As the encoder's: decaying, first tap large, ±32767 signs.
            let scale = 2 + (g.u64() % 4) as u32;
            for (i, v) in h.iter_mut().enumerate() {
                *v = (g.u64() as i16) >> (scale + (i as u32 / 6).min(12));
            }
            h[0] = 4096 + (g.u64() % 8192) as i16;
            for s in &mut sign {
                *s = if g.u64() & 1 == 0 { 32767 } else { -32767 };
            }
        }
        let mut o = 0;
        unsafe { c_cor_h(h.as_mut_ptr(), sign.as_mut_ptr(), rr.as_mut_ptr(), &mut o) };
    }
    let mut dn = [0i16; 40];
    if flat {
        for d in &mut dn {
            *d = g.i16() >> (6 + g.u64() % 4);
        }
    } else if !odd {
        // the timing class: `dn` as `cor_h_x` normalises it (2 sum of the tracks' max |dn|
        // within 32767: the fast path; full-scale values took the exact search every time)
        for d in &mut dn {
            *d = g.i16() >> (4 + g.u64() % 4);
        }
    } else if odd && g.u64().is_multiple_of(4) {
        for d in &mut dn {
            *d = [0, 1, 8192, 16384, -16384, -8192][(g.u64() % 6) as usize];
        }
    } else {
        signal(g, &mut dn);
    }
    let (t, n) = (T as u64, 40 / T as u64);
    let a = (g.u64() % t) as i16;
    let ipos = core::array::from_fn(|k| (a + k as i16) % t as i16);
    let pos_max = core::array::from_fn(|p| p as i16 + (t * (g.u64() % n)) as i16);
    (dn, rr, ipos, pos_max)
}

/// `Pair::best` (per-lane argmax as the rows are made) is the compare chain over the pair's
/// candidates, both forms (`L32`): random rows (energies of any sign and size: the
/// fallback), a tie-heavy class (few distinct values), both energy forms, random padding
/// lanes; and an edge class for the range check: valid wide rows with one cell at an
/// energy of 0, -1, 1, 32767 or 32768, or ps = -32768 (sq = -32768 in the C), also where
/// only they could differ (that cell first, all others sq = 0 or it also sq = 0).
fn best_is_the_chain<const R: usize, const W: usize>(seed: u64) {
    let mut g = Rng::new(seed);
    let (mut fast, mut edges, mut edge_fast) = (0, 0, 0);
    for case in 0..N {
        if case % 4 == 1 {
            let ok = best_edge_case::<R, W>(&mut g, edges % 4);
            (edges, edge_fast) = (edges + 1, edge_fast + ok as usize);
            continue;
        }
        let odd = g.u64().is_multiple_of(4);
        let wide = g.u64() & 1 == 0;
        let pad_ps = g.i16();
        let (mut d1, mut tps, mut base, mut tal) = ([0i16; W], [0i16; W], [0i32; W], [0i32; W]);
        let mut rr = [[0i16; W]; R];
        for i in 0..W {
            if i >= R {
                // Padding as `search_pairs` makes it: rr 0, tal 0, tps = ps + a zero dn
                // lane (d1 is per row: unused here).
                (tps[i], base[i], tal[i]) = (pad_ps, base[0], 0);
            } else if odd {
                (d1[i], tps[i]) = (g.i16(), g.i16());
                // The wide pair's adds don't wrap in the C either: `alp1` is a few i16
                // terms (|alp1| < 2^18), `tal` an i16. The other pairs wrap: any i32.
                (base[i], tal[i]) = if wide {
                    (g.i32() >> (13 + g.u64() % 19), g.i16() as i32)
                } else {
                    (g.i32() >> (g.u64() % 32), g.i32() >> (g.u64() % 32))
                };
                for row in &mut rr {
                    row[i] = g.i16();
                }
            } else {
                d1[i] = few(&mut g, &[0, 3000, 6000]) as i16;
                tps[i] = few(&mut g, &[0, 3000, 9000]) as i16;
                // Energies in 1..=32767 either way (wide: (base + rr) / 8 + tal / 2).
                (base[i], tal[i]) = if wide {
                    (
                        few(&mut g, &[8, 4000, 40000]),
                        few(&mut g, &[2, 1000, 8000]),
                    )
                } else {
                    (few(&mut g, &[1 << 16, 1 << 28]), few(&mut g, &[0, 1 << 20]))
                };
                for row in &mut rr {
                    row[i] = few(&mut g, &[0, 100, 2000]) as i16;
                }
            }
        }
        let p = Pair::<R, W> {
            d1: &d1,
            tps,
            base,
            tal,
            rr: &rr,
            sh: 10 + (g.u64() % 3) as u32,
            wide,
        };
        let c = p.cands();
        let ok =
            (0..R).all(|j| (0..R).all(|i| c.sq[j][i] >= 0 && (1..=32767).contains(&c.al[j][i])));
        fast += ok as usize;
        let want = chain(&c);
        assert_eq!(p.best::<false>(), want, "{:?} {:?}", c.sq, c.al);
        assert_eq!(p.best::<true>(), want, "{:?} {:?}", c.sq, c.al);
    }
    assert!(fast > 3 * N / 8, "{fast}");
    // Both sides of the range check: valid edge cells (kind 2, energy 1 or 32767) take the
    // lanes, the others fall back.
    assert!(
        edge_fast > edges / 20 && edge_fast < edges / 4,
        "{edge_fast} of {edges}"
    );
}

/// One of `v`, uniformly.
fn few(g: &mut Rng, v: &[i32]) -> i32 {
    v[(g.u64() % v.len() as u64) as usize]
}

/// One case of `best_is_the_chain`'s edge class; whether every candidate is in the fast
/// range. Wide energies `(((base + rr) >> 2) + tal) >> 1` with base 800, so tal 0 or 65534
/// and rr pick the energy exactly. Kinds: 0 ps = -32768 anywhere, 1 on all of row 0 (one tps)
/// with every other sq 0 (lanes reading it as sq 0 take it first), 2 an edge energy
/// anywhere, 3 energy 0 and sq 0 at (0, 0) (nothing beats 0 / 0 in a lane; the C skips it).
fn best_edge_case<const R: usize, const W: usize>(g: &mut Rng, kind: usize) -> bool {
    let (mut d1, mut tps, mut tal) = ([0i16; W], [0i16; W], [0i32; W]);
    let mut rr = [[0i16; W]; R];
    let ps: &[i32] = if kind == 1 {
        &[0, 100]
    } else {
        &[0, 3000, 9000]
    };
    let pad_ps = if kind == 1 { 0 } else { few(g, ps) as i16 };
    for i in 0..W {
        if i >= R {
            tps[i] = pad_ps;
            continue;
        }
        d1[i] = few(g, ps) as i16;
        tps[i] = if kind == 1 && i > 0 {
            tps[0]
        } else {
            few(g, ps) as i16
        };
        tal[i] = few(g, &[0, 65534]);
        // Energies 1, 100, 1000 (tal 0) or 32767, 32757 (tal 65534): rr - 800 = 8 al - tal.
        let v: &[i32] = if tal[i] == 0 {
            &[8, 800, 8000]
        } else {
            &[0, -80]
        };
        for row in &mut rr {
            row[i] = (few(g, v) - 800) as i16;
        }
    }
    let (j0, i0) = match kind {
        1 | 3 => (0, 0),
        _ => ((g.u64() % R as u64) as usize, (g.u64() % R as u64) as usize),
    };
    match kind {
        // ps = d1[j0] + tps[i0] = -32768.
        0 | 1 => d1[j0] = (-32768 - tps[i0] as i32) as i16,
        // al 0, -1 or 1 (tal 0); 32768 or 32767 (tal 65534).
        2 => {
            let v: &[i32] = if tal[i0] == 0 { &[0, -8, 8] } else { &[8, 0] };
            rr[j0][i0] = (few(g, v) - 800) as i16;
        }
        _ => {
            (tal[i0], d1[j0]) = (0, -tps[i0]);
            rr[j0][i0] = -800;
        }
    }
    let p = Pair::<R, W> {
        d1: &d1,
        tps,
        base: [800; W],
        tal,
        rr: &rr,
        sh: 0,
        wide: true,
    };
    let c = p.cands();
    let ok = (0..R).all(|j| (0..R).all(|i| c.sq[j][i] >= 0 && (1..=32767).contains(&c.al[j][i])));
    let want = chain(&c);
    assert_eq!(p.best::<false>(), want, "{:?} {:?}", c.sq, c.al);
    assert_eq!(p.best::<true>(), want, "{:?} {:?}", c.sq, c.al);
    ok
}
