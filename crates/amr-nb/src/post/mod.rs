//! Decoder post-processing group: LP synthesis/analysis filters (`Syn_filt`, `Residu`,
//! `Weight_Ai`), `preemphasis`, `agc`/`agc2`, the formant post-filter (`Post_Filter`) and the
//! output high-pass filter (`Post_Process`).
//! Bit-exact ports of opencore-amr 0.1.6 (Apache-2.0) `common/src/{syn_filt,residu,weight_a}.cpp`
//! and `dec/src/{preemph,agc,pstfilt,post_pro}.cpp`.

pub(crate) use crate::cnst::{CodecMode, L_FRAME, L_SUBFR, M, MP1};
pub(crate) mod agc;
pub(crate) mod filt;
pub(crate) mod post_pro;
pub(crate) mod preemph;
pub(crate) mod pstfilt;

// Post-filter constants (`pstfilt.cpp`, local to this group).
/// `MU`: post-filter tilt factor.
const MU: i16 = 26214;
/// `AGC_FAC`: AGC smoothing factor.
const AGC_FAC: i16 = 29491;

/// Shared input generators for this group's differential tests.
#[cfg(test)]
pub(crate) mod test_gen {
    use super::M;
    use crate::testutil::Rng;

    /// A stable LP filter in Q12 like the encoder's: step-up from random reflection
    /// coefficients (|k| < 0.9), redrawn until every coefficient fits 16 bits.
    pub(crate) fn lpc(g: &mut Rng) -> [i16; M + 1] {
        loop {
            let mut a = [0f64; M + 1];
            a[0] = 1.0;
            for i in 1..=M {
                let k = ((g.u64() % 1801) as f64 - 900.0) / 1000.0;
                let prev = a;
                for j in 1..i {
                    a[j] = prev[j] + k * prev[i - j];
                }
                a[i] = k;
            }
            if a.iter().all(|v| (v * 4096.0).abs() < 32767.0) {
                return a.map(|v| (v * 4096.0).round() as i16);
            }
        }
    }

    /// LPC-like coefficients: Q12 `a[0] = 4096` with taps of random magnitude, or (1 in 3)
    /// fully random 16-bit values (unstable filters, drives the saturation paths).
    pub(crate) fn coefs(g: &mut Rng) -> [i16; M + 1] {
        let mut a = [0i16; M + 1];
        if g.u64().is_multiple_of(3) {
            a.iter_mut().for_each(|v| *v = g.i16());
        } else {
            let scale = (g.u64() % 14) as u32;
            a[0] = 4096;
            for v in &mut a[1..] {
                *v = (g.u64() as i16) >> scale;
            }
        }
        a
    }

    /// Fills `buf` with a signal of random amplitude, sprinkled with edge values.
    pub(crate) fn signal(g: &mut Rng, buf: &mut [i16]) {
        let scale = (g.u64() % 16) as u32;
        for v in buf {
            *v = if g.u64().is_multiple_of(16) {
                g.i16()
            } else {
                (g.u64() as i16) >> scale
            };
        }
    }
}
