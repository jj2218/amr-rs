//! Algebraic codebook decoders (`d2_9pf`, `d2_11pf`, `d3_14pf`, `d4_17pf`, `d8_31pf`, `d1035pf`).
//! Table-pointer arguments (`startPos_ptr`, `dgray_ptr`) are replaced by the tables themselves.

use super::L_SUBFR;
use crate::basic_op::{Flag, add, l_mult, l_shr, mult, negate, shl};
use crate::tables::{DGRAY, STARTPOS};

/// Writes the pulse amplitudes of `pos` into the zeroed `cod`, as all the 2/3/4-pulse
/// decoders do (`cod[pos[j]] = i * 16383 - 8192`, i.e. +8191 or -8192).
fn place_pulses(pos: &[i16], mut sign: i16, cod: &mut [i16]) {
    for p in pos {
        let i = sign & 1;
        cod[*p as usize] = i.wrapping_mul(16383).wrapping_sub(8192);
        sign >>= 1;
    }
}

/// `decode_2i40_9bits` (MR475/MR515).
pub(crate) fn decode_2i40_9bits(
    sub_nr: i16,
    sign: i16,
    mut index: i16,
    cod: &mut [i16],
    ov: &mut Flag,
) {
    let mut pos = [0i16; 2];
    let mut j = index & 64;
    j >>= 3;
    let mut i = index & 7;
    let mut k = shl(sub_nr, 1, ov);
    k = k.wrapping_add(j);
    pos[0] = i.wrapping_mul(5).wrapping_add(STARTPOS[k as usize]);
    k += 1;
    index >>= 3;
    i = index & 7;
    pos[1] = i.wrapping_mul(5).wrapping_add(STARTPOS[k as usize]);
    cod[..L_SUBFR].fill(0);
    place_pulses(&pos, sign, cod);
}

/// `decode_2i40_11bits` (MR59).
pub(crate) fn decode_2i40_11bits(sign: i16, mut index: i16, cod: &mut [i16]) {
    let mut pos = [0i16; 2];
    let mut j = index & 1;
    index >>= 1;
    let mut i = index & 7;
    pos[0] = i * 5 + j * 2 + 1;
    index >>= 3;
    j = index & 3;
    index >>= 2;
    i = index & 7;
    pos[1] = if j == 3 { i * 5 + 4 } else { i * 5 + j };
    cod[..L_SUBFR].fill(0);
    place_pulses(&pos, sign, cod);
}

/// `decode_3i40_14bits` (MR67).
pub(crate) fn decode_3i40_14bits(sign: i16, mut index: i16, cod: &mut [i16]) {
    let mut pos = [0i16; 3];
    let mut i = index & 7;
    pos[0] = i * 5;
    index >>= 3;
    let mut j = index & 1;
    index >>= 1;
    i = index & 7;
    pos[1] = i * 5 + j * 2 + 1;
    index >>= 3;
    j = index & 1;
    index >>= 1;
    i = index & 7;
    pos[2] = i * 5 + j * 2 + 2;
    cod[..L_SUBFR].fill(0);
    place_pulses(&pos, sign, cod);
}

/// `decode_4i40_17bits` (MR74/MR795).
pub(crate) fn decode_4i40_17bits(sign: i16, mut index: i16, cod: &mut [i16]) {
    let mut pos = [0i16; 4];
    let mut i = DGRAY[(index & 7) as usize];
    pos[0] = i * 5;
    index >>= 3;
    i = DGRAY[(index & 7) as usize];
    pos[1] = i * 5 + 1;
    index >>= 3;
    i = DGRAY[(index & 7) as usize];
    pos[2] = i * 5 + 2;
    index >>= 3;
    let j = index & 1;
    index >>= 1;
    i = DGRAY[(index & 7) as usize];
    pos[3] = i * 5 + 3 + j;
    cod[..L_SUBFR].fill(0);
    place_pulses(&pos, sign, cod);
}

/// `decompress10` (static in d8_31pf.cpp).
fn decompress10(
    mut msbs: i16,
    lsbs: i16,
    index1: usize,
    index2: usize,
    index3: usize,
    pos_indx: &mut [i16],
    ov: &mut Flag,
) {
    if msbs > 124 {
        msbs = 124;
    }
    let mut ia = mult(msbs, 1311, ov);
    let mut t = l_mult(ia, 25, ov);
    ia = (msbs as i32).wrapping_sub(t >> 1) as i16;
    let mut ib = mult(ia, 6554, ov);
    t = l_mult(ib, 5, ov);
    ib = ia.wrapping_sub((t >> 1) as i16);
    ib = shl(ib, 1, ov);
    let mut ic = lsbs.wrapping_sub((lsbs >> 2) << 2);
    pos_indx[index1] = ib.wrapping_add(ic & 1);
    ib = mult(ia, 6554, ov);
    ib = shl(ib, 1, ov);
    pos_indx[index2] = ib.wrapping_add(ic >> 1);
    ib = lsbs >> 2;
    ic = mult(msbs, 1311, ov);
    ic = shl(ic, 1, ov);
    pos_indx[index3] = add(ib, ic, ov);
}

/// `NB_TRACK_MR102` (cnst.h).
const NB_TRACK_MR102: usize = 4;

/// `decompress_code` (static in d8_31pf.cpp).
fn decompress_code(indx: &[i16], sign_indx: &mut [i16], pos_indx: &mut [i16], ov: &mut Flag) {
    sign_indx[..NB_TRACK_MR102].copy_from_slice(&indx[..NB_TRACK_MR102]);
    let mut msbs = indx[NB_TRACK_MR102] >> 3;
    let mut lsbs = indx[NB_TRACK_MR102] & 7;
    decompress10(msbs, lsbs, 0, 4, 1, pos_indx, ov);
    msbs = indx[NB_TRACK_MR102 + 1] >> 3;
    lsbs = indx[NB_TRACK_MR102 + 1] & 7;
    decompress10(msbs, lsbs, 2, 6, 5, pos_indx, ov);
    msbs = indx[NB_TRACK_MR102 + 2] >> 2;
    lsbs = indx[NB_TRACK_MR102 + 2] & 3;
    let mut t = l_mult(msbs, 25, ov);
    let mut ia = l_shr(t, 1, ov) as i16;
    ia = ia.wrapping_add(12);
    let msbs0_24 = ia >> 5;
    ia = mult(msbs0_24, 6554, ov);
    ia &= 1;
    let mut ib = mult(msbs0_24, 6554, ov);
    t = l_mult(ib, 5, ov);
    ib = msbs0_24.wrapping_sub((t >> 1) as i16);
    if ia == 1 {
        ib = 4i16.wrapping_sub(ib);
    }
    ib = shl(ib, 1, ov);
    ia = lsbs & 1;
    pos_indx[3] = add(ib, ia, ov);
    ia = mult(msbs0_24, 6554, ov);
    ia = shl(ia, 1, ov);
    pos_indx[7] = ia.wrapping_add(lsbs >> 1);
}

/// `dec_8i40_31bits` (MR102): `index` holds the 7 codebook parameters.
pub(crate) fn dec_8i40_31bits(index: &[i16], cod: &mut [i16], ov: &mut Flag) {
    const POS_CODE: i16 = 8191;
    const NEG_CODE: i16 = 8191;
    let mut linear_signs = [0i16; NB_TRACK_MR102];
    let mut linear_codewords = [0i16; 8];
    cod[..L_SUBFR].fill(0);
    decompress_code(index, &mut linear_signs, &mut linear_codewords, ov);
    for j in 0..NB_TRACK_MR102 {
        let pos1 = (linear_codewords[j] << 2).wrapping_add(j as i16);
        let mut sign = if linear_signs[j] == 0 {
            POS_CODE
        } else {
            -NEG_CODE
        };
        if pos1 < L_SUBFR as i16 {
            cod[pos1 as usize] = sign;
        }
        let pos2 = (linear_codewords[j + 4] << 2).wrapping_add(j as i16);
        if pos2 < pos1 {
            sign = negate(sign);
        }
        if pos2 < L_SUBFR as i16 {
            cod[pos2 as usize] = cod[pos2 as usize].wrapping_add(sign);
        }
    }
}

/// `dec_10i40_35bits` (MR122): `index` holds the 10 codebook parameters.
pub(crate) fn dec_10i40_35bits(index: &[i16], cod: &mut [i16]) {
    const NB_TRACK: usize = 5;
    cod[..L_SUBFR].fill(0);
    for j in 0..NB_TRACK {
        let tmp = index[j];
        let mut i = DGRAY[(tmp & 7) as usize];
        i = i.wrapping_mul(5);
        let pos1 = i.wrapping_add(j as i16);
        i = (tmp >> 3) & 1;
        let mut sign: i16 = if i == 0 { 4096 } else { -4096 };
        cod[pos1 as usize] = sign;
        i = DGRAY[(index[j + 5] & 7) as usize];
        i = i.wrapping_mul(5);
        let pos2 = i.wrapping_add(j as i16);
        if pos2 < pos1 {
            sign = negate(sign);
        }
        cod[pos2 as usize] = cod[pos2 as usize].wrapping_add(sign);
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::testutil::{N, Rng};
    use opencore_sys as _;

    // The TS 26.073 reference, in the shape of the opencore functions these tests were
    // written against (`ov`/`o`: its overflow flag; table arguments unused).
    unsafe extern "C" {
        fn ref_decode_2i40_9bits(sub: i16, sign: i16, index: i16, cod: *mut i16);
        fn ref_decode_2i40_11bits(sign: i16, index: i16, cod: *mut i16);
        fn ref_decode_3i40_14bits(sign: i16, index: i16, cod: *mut i16);
        fn ref_decode_4i40_17bits(sign: i16, index: i16, cod: *mut i16);
        fn ref_dec_8i40_31bits(index: *mut i16, cod: *mut i16);
        fn ref_dec_10i40_35bits(index: *mut i16, cod: *mut i16);
    }
    use crate::testutil::with_ref;
    unsafe fn c_decode_2i40_9bits(
        sub: i16,
        sign: i16,
        index: i16,
        _sp: *const i16,
        cod: *mut i16,
        o: *mut i32,
    ) {
        with_ref(o, || unsafe {
            ref_decode_2i40_9bits(sub, sign, index, cod)
        })
    }
    unsafe fn c_decode_2i40_11bits(sign: i16, index: i16, cod: *mut i16) {
        with_ref(std::ptr::null_mut(), || unsafe {
            ref_decode_2i40_11bits(sign, index, cod)
        })
    }
    unsafe fn c_decode_3i40_14bits(sign: i16, index: i16, cod: *mut i16) {
        with_ref(std::ptr::null_mut(), || unsafe {
            ref_decode_3i40_14bits(sign, index, cod)
        })
    }
    unsafe fn c_decode_4i40_17bits(sign: i16, index: i16, _dgray: *const i16, cod: *mut i16) {
        with_ref(std::ptr::null_mut(), || unsafe {
            ref_decode_4i40_17bits(sign, index, cod)
        })
    }
    unsafe fn c_dec_8i40_31bits(index: *mut i16, cod: *mut i16, o: *mut i32) {
        with_ref(o, || unsafe { ref_dec_8i40_31bits(index, cod) })
    }
    unsafe fn c_dec_10i40_35bits(index: *mut i16, cod: *mut i16, _dgray: *const i16) {
        with_ref(std::ptr::null_mut(), || unsafe {
            ref_dec_10i40_35bits(index, cod)
        })
    }

    /// Random value with `n` bits (the parameter's width in the bitstream).
    fn bits(g: &mut Rng, n: u32) -> i16 {
        (g.u64() & ((1u64 << n) - 1)) as i16
    }

    /// Output buffer pre-filled with junk, so a missed zeroing would show.
    fn junk(g: &mut Rng) -> [i16; L_SUBFR] {
        core::array::from_fn(|_| g.i16())
    }

    #[test]
    fn small_codebooks_vs_ref() {
        let mut g = Rng::new(1);
        for _ in 0..N {
            // decode_2i40_9bits: subNr 0..=3, index 7 bits, sign 2 bits (MR475/MR515).
            let (sub, sign, idx) = ((g.u64() % 4) as i16, bits(&mut g, 2), bits(&mut g, 7));
            let mut r = junk(&mut g);
            let mut c = r;
            let (mut ro, mut co) = (false, 0);
            decode_2i40_9bits(sub, sign, idx, &mut r, &mut ro);
            unsafe {
                c_decode_2i40_9bits(sub, sign, idx, STARTPOS.as_ptr(), c.as_mut_ptr(), &mut co)
            };
            assert_eq!(
                (r, ro),
                (c, co != 0),
                "decode_2i40_9bits({sub},{sign},{idx})"
            );

            // decode_2i40_11bits: index 9 bits, sign 2 bits (MR59).
            let (sign, idx) = (bits(&mut g, 2), bits(&mut g, 9));
            let mut r = junk(&mut g);
            let mut c = r;
            decode_2i40_11bits(sign, idx, &mut r);
            unsafe { c_decode_2i40_11bits(sign, idx, c.as_mut_ptr()) };
            assert_eq!(r, c, "decode_2i40_11bits({sign},{idx})");

            // decode_3i40_14bits: index 11 bits, sign 3 bits (MR67).
            let (sign, idx) = (bits(&mut g, 3), bits(&mut g, 11));
            let mut r = junk(&mut g);
            let mut c = r;
            decode_3i40_14bits(sign, idx, &mut r);
            unsafe { c_decode_3i40_14bits(sign, idx, c.as_mut_ptr()) };
            assert_eq!(r, c, "decode_3i40_14bits({sign},{idx})");

            // decode_4i40_17bits: index 13 bits, sign 4 bits (MR74/MR795).
            let (sign, idx) = (bits(&mut g, 4), bits(&mut g, 13));
            let mut r = junk(&mut g);
            let mut c = r;
            decode_4i40_17bits(sign, idx, &mut r);
            unsafe { c_decode_4i40_17bits(sign, idx, DGRAY.as_ptr(), c.as_mut_ptr()) };
            assert_eq!(r, c, "decode_4i40_17bits({sign},{idx})");
        }
    }

    #[test]
    fn dec_8i40_31bits_vs_ref() {
        let mut g = Rng::new(2);
        for n in 0..N {
            // MR102: 4 sign bits, two 10-bit and one 7-bit position words. Every 8th case
            // forces MSBs >= 124 to hit the clamp in decompress10.
            let mut idx = [
                bits(&mut g, 1),
                bits(&mut g, 1),
                bits(&mut g, 1),
                bits(&mut g, 1),
                bits(&mut g, 10),
                bits(&mut g, 10),
                bits(&mut g, 7),
            ];
            if n % 8 == 0 {
                idx[4] |= 0x3e0;
                idx[5] |= 0x3e0;
            }
            let mut r = junk(&mut g);
            let mut c = r;
            let mut ci = idx;
            let (mut ro, mut co) = (false, 0);
            dec_8i40_31bits(&idx, &mut r, &mut ro);
            unsafe { c_dec_8i40_31bits(ci.as_mut_ptr(), c.as_mut_ptr(), &mut co) };
            assert_eq!((r, ro), (c, co != 0), "dec_8i40_31bits({idx:?})");
        }
    }

    #[test]
    fn dec_10i40_35bits_vs_ref() {
        let mut g = Rng::new(3);
        for _ in 0..N {
            // MR122: five 4-bit (sign + position) and five 3-bit position words.
            let idx: [i16; 10] = core::array::from_fn(|k| bits(&mut g, if k < 5 { 4 } else { 3 }));
            let mut r = junk(&mut g);
            let mut c = r;
            let mut ci = idx;
            dec_10i40_35bits(&idx, &mut r);
            unsafe { c_dec_10i40_35bits(ci.as_mut_ptr(), c.as_mut_ptr(), DGRAY.as_ptr()) };
            assert_eq!(r, c, "dec_10i40_35bits({idx:?})");
        }
    }
}
