//! Pre-processing (vo-amrwbenc `decim54.c`): 16 kHz → 12.8 kHz decimation.

use crate::enc_tables::FIR_DOWN1;

/// Filter half-length (`NB_COEF_DOWN`): `Decim_12k8` keeps 2·15 input samples of memory.
pub(crate) const NB_COEF_DOWN: usize = 15;
/// 4/5 in Q15 (`DOWN_FAC`).
const DOWN_FAC: i32 = 26215;
/// Input step per output sample, in quarters (5/4, `FAC5`).
const FAC5: usize = 5;

/// `Decim_12k8`: decimates `sig16k` (`lg` samples) to `sig12k8` (`lg·4/5` samples).
pub(crate) fn decim_12k8(sig16k: &[i16], sig12k8: &mut [i16], mem: &mut [i16; 2 * NB_COEF_DOWN]) {
    let lg = sig16k.len();
    let mut signal = [0i16; 320 + 2 * NB_COEF_DOWN];
    signal[..2 * NB_COEF_DOWN].copy_from_slice(mem);
    signal[2 * NB_COEF_DOWN..2 * NB_COEF_DOWN + lg].copy_from_slice(sig16k);
    let lg_down = ((lg as i32 * DOWN_FAC) >> 15) as usize;
    down_samp(&signal[..2 * NB_COEF_DOWN + lg], &mut sig12k8[..lg_down]);
    mem.copy_from_slice(&signal[lg..lg + 2 * NB_COEF_DOWN]);
}

/// `Down_samp` (on the whole buffer, memory first; the C's `sig` is `signal + 15`): output
/// `j` interpolates at 5j/4 with the 30-tap phase `frac` of
/// `fir_down1`. The reference's `Interpol` is an `L_mac` chain, `L_shl(·, 1)` and `round`; no
/// partial sum can saturate (asserted), so it is one plain sum, then clamp(4·sum) and round.
fn down_samp(sig: &[i16], sig_d: &mut [i16]) {
    const _: () = {
        let mut f = 0;
        while f < 4 {
            let (mut s, mut k) = (0i64, 0);
            while k < 2 * NB_COEF_DOWN {
                s += (FIR_DOWN1[30 * f + k] as i64).abs();
                k += 1;
            }
            assert!(2 * 32768 * s < 1 << 31);
            f += 1;
        }
    };
    let mut pos = 0;
    for out in sig_d {
        let (i, frac) = (pos >> 2, pos & 3);
        // the C's x = sig + i - NB_COEF_DOWN + 1 = signal + i + 1
        let x = &sig[i + 1..i + 1 + 2 * NB_COEF_DOWN];
        let fir = &FIR_DOWN1[30 * frac..30 * frac + 2 * NB_COEF_DOWN];
        let s: i32 = x.iter().zip(fir).map(|(&a, &b)| a as i32 * b as i32).sum();
        let l = (4 * s as i64).clamp(i32::MIN as i64, i32::MAX as i64) as i32;
        *out = (l.saturating_add(0x8000) >> 16) as i16;
        pos += FAC5;
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_Decim_12k8(sig16k: *mut i16, lg: i16, sig12k8: *mut i16, mem: *mut i16);
    }

    /// Frames (320) and the coder's lookahead calls (L_FILT16k = 15) over 1000-call sequences:
    /// random, loud (full scale) and square waves (the saturating round); outputs and memory.
    #[test]
    fn decim_12k8_matches_reference() {
        let mut r = Rng::new(11);
        for seq in 0..12 {
            let (mut m, mut c) = ([0i16; 30], [0i16; 30]);
            let mut t = 0u64;
            for call in 0..1000 {
                let lg = if call % 3 == 2 { 15 } else { 320 };
                let x: Vec<i16> = (0..lg)
                    .map(|_| {
                        t += 1;
                        match seq % 3 {
                            0 => r.i16(),
                            1 => [i16::MAX, i16::MIN][((t / 7) % 2) as usize],
                            _ => r.i16() >> (r.u64() % 12),
                        }
                    })
                    .collect();
                let (mut a, mut b) = ([0i16; 256], [0i16; 256]);
                let mut xc = x.clone();
                decim_12k8(&x, &mut a, &mut m);
                with_ref(|| {
                    // SAFETY: lg input samples, lg·4/5 outputs, the 30-word memory.
                    unsafe {
                        refwb_Decim_12k8(xc.as_mut_ptr(), lg as i16, b.as_mut_ptr(), c.as_mut_ptr())
                    }
                });
                assert_eq!((a, m), (b, c), "sequence {seq} call {call}");
            }
        }
    }
}
