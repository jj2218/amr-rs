//! Encoder-side DTX (vo-amrwbenc `dtx.c`: `dtx_enc`, `dtx_buffer`, `tx_dtx_handler`), with the
//! reference's saturation: ISF and energy history, SID parameters, comfort-noise excitation.

use super::bits::Serial;
use super::lpc::qisf_ns;
use super::ops::*;
use crate::basic_op::dot_product12;
use crate::cnst::{L_FRAME12K8, M, RANDOM_INITSEED};
use crate::enc_tables::EN_ADJUST;
use crate::excitation::noise_gen_amrwb as random;

const DTX_HIST_SIZE: usize = 8;
const DTX_ELAPSED_FRAMES_THRESH: i16 = 24 + 7 - 1;
pub(crate) const DTX_HANG_CONST: i16 = 7;
const INV_MED_THRESH: i16 = 14564;
const GAIN_THR: i16 = 180;
/// `MRDTX`: the mode number of SID / no-data frames.
pub(crate) const MRDTX: i16 = 9;

/// `dtx_encState`.
#[derive(Clone, Debug, PartialEq)]
#[repr(C)]
pub(crate) struct DtxEncState {
    isf_hist: [i16; M * DTX_HIST_SIZE],
    log_en_hist: [i16; DTX_HIST_SIZE],
    hist_ptr: i16,
    log_en_index: i16,
    cng_seed: i16,
    pub(crate) dtx_hangover_count: i16,
    dec_ana_elapsed_count: i16,
    d: [i32; 28],
    sum_d: [i32; DTX_HIST_SIZE],
}

impl DtxEncState {
    pub(crate) fn new(isf_init: &[i16; M]) -> Self {
        let mut s = Self {
            isf_hist: [0; M * DTX_HIST_SIZE],
            log_en_hist: [0; DTX_HIST_SIZE],
            hist_ptr: 0,
            log_en_index: 0,
            cng_seed: 0,
            dtx_hangover_count: 0,
            dec_ana_elapsed_count: 0,
            d: [0; 28],
            sum_d: [0; DTX_HIST_SIZE],
        };
        s.reset(isf_init);
        s
    }

    /// `dtx_enc_reset` (it zeroes `sumD[0..7]` only; `sumD[7]` is rewritten before it is read).
    pub(crate) fn reset(&mut self, isf_init: &[i16; M]) {
        self.hist_ptr = 0;
        self.log_en_index = 0;
        for h in self.isf_hist.as_chunks_mut::<M>().0 {
            h.copy_from_slice(isf_init);
        }
        self.cng_seed = RANDOM_INITSEED;
        self.log_en_hist = [0; DTX_HIST_SIZE];
        self.dtx_hangover_count = DTX_HANG_CONST;
        self.dec_ana_elapsed_count = 32767;
        self.d = [0; 28];
        self.sum_d[..DTX_HIST_SIZE - 1].fill(0);
    }

    /// `dtx_enc`: the SID parameters (written to `prms`), the comfort-noise ISFs (`isf`) and
    /// excitation (`exc2`).
    pub(crate) fn enc(&mut self, isf: &mut [i16; M], exc2: &mut [i16], prms: &mut Serial) {
        let mut log_en = 0i16;
        for &h in &self.log_en_hist {
            log_en = add(log_en, h);
        }
        let mut isf_order = [0i16; 3];
        self.find_frame_indices(&mut isf_order);
        let mut l_isf = [0i32; M];
        aver_isf_history(&mut self.isf_hist, &isf_order, &mut l_isf);
        for j in 0..M {
            isf[j] = extract_l(l_shr(l_isf[j], 3));
        }
        let log_en = mult(add(shr(log_en, 2), 512), 21504);
        self.log_en_index = shr(log_en, 6).clamp(0, 63);
        let mut indice = [0i16; 5];
        let isf_in = *isf;
        qisf_ns(&isf_in, isf, &mut indice);
        for (k, n) in [6, 6, 6, 5, 5].into_iter().enumerate() {
            prms.parm(indice[k], n);
        }
        prms.parm(self.log_en_index, 6);
        let cn_dith = self.dithering_control();
        prms.parm(cn_dith, 1);
        let log_en = mult(shl(self.log_en_index, 15 - 6), 12483);
        let log_en_int_e = add(shr(log_en, 10), 16 - 1);
        let log_en_int_m = shl(log_en & 0x3ff, 5);
        let level32 = pow2(log_en_int_e, log_en_int_m);
        let exp0 = norm_l(level32);
        let level = extract_h(l_shl(level32, exp0));
        let exp0 = sub(15, exp0);
        for v in exc2[..L_FRAME12K8].iter_mut() {
            *v = shr(random(&mut self.cng_seed), 4);
        }
        let mut exp = 0;
        let mut ener32 = dot_product12(&exc2[..L_FRAME12K8], &exc2[..L_FRAME12K8], &mut exp);
        isqrt_n(&mut ener32, &mut exp);
        let gain = mult(level, extract_h(ener32));
        let exp = add(add(exp0, exp), 4);
        for v in exc2[..L_FRAME12K8].iter_mut() {
            *v = shl(mult(*v, gain), exp);
        }
    }

    /// `dtx_buffer`: stores the frame's ISFs and log energy of `enr` (residual, `L_FRAME12K8`).
    pub(crate) fn buffer(&mut self, isf_new: &[i16; M], enr: i32, codec_mode: i16) {
        self.hist_ptr = add(self.hist_ptr, 1);
        if self.hist_ptr == DTX_HIST_SIZE as i16 {
            self.hist_ptr = 0;
        }
        let p = self.hist_ptr as usize * M;
        self.isf_hist[p..p + M].copy_from_slice(isf_new);
        let (e, m) = log2(enr);
        let log_en = add(shl(e, 7), shr(m, 15 - 7));
        let log_en = sub(log_en, add(1024, EN_ADJUST[codec_mode as usize]));
        self.log_en_hist[self.hist_ptr as usize] = log_en;
    }

    /// `tx_dtx_handler`: the DTX hangover; sets `used_mode` to `MRDTX` once it has run out.
    pub(crate) fn tx_handler(&mut self, vad_flag: i16, used_mode: &mut i16) {
        self.dec_ana_elapsed_count = add(self.dec_ana_elapsed_count, 1);
        if vad_flag != 0 {
            self.dtx_hangover_count = DTX_HANG_CONST;
        } else if self.dtx_hangover_count == 0 {
            self.dec_ana_elapsed_count = 0;
            *used_mode = MRDTX;
        } else {
            self.dtx_hangover_count = sub(self.dtx_hangover_count, 1);
            if add(self.dec_ana_elapsed_count, self.dtx_hangover_count) < DTX_ELAPSED_FRAMES_THRESH
            {
                *used_mode = MRDTX;
            }
        }
    }

    fn find_frame_indices(&mut self, indices: &mut [i16; 3]) {
        let mut tmp = (DTX_HIST_SIZE - 1) as i16;
        let mut j = -1i16;
        for i in 0..DTX_HIST_SIZE - 1 {
            j = add(j, tmp);
            self.sum_d[i] = l_sub(self.sum_d[i], self.d[j as usize]);
            tmp -= 1;
        }
        self.sum_d.copy_within(0..DTX_HIST_SIZE - 1, 1);
        self.sum_d[0] = 0;
        let mut tmp = 0usize;
        let mut i = 27usize;
        while i >= 12 {
            tmp += 1;
            for j in (1..=tmp).rev() {
                self.d[i - j + 1] = self.d[i - j - tmp];
            }
            i -= tmp;
        }
        let hp = self.hist_ptr as usize;
        let mut ptr = hp as isize;
        for i in 1..DTX_HIST_SIZE {
            ptr -= 1;
            if ptr < 0 {
                ptr = (DTX_HIST_SIZE - 1) as isize;
            }
            let mut l = 0i32;
            for j in 0..M {
                let t = sub(
                    self.isf_hist[hp * M + j],
                    self.isf_hist[ptr as usize * M + j],
                );
                l = l_mac(l, t, t);
            }
            self.d[i - 1] = l;
            self.sum_d[0] = l_add(self.sum_d[0], l);
            self.sum_d[i] = l_add(self.sum_d[i], l);
        }
        let (mut summax, mut summin) = (self.sum_d[0], self.sum_d[0]);
        indices[0] = 0;
        indices[2] = 0;
        for i in 1..DTX_HIST_SIZE {
            if self.sum_d[i] > summax {
                indices[0] = i as i16;
                summax = self.sum_d[i];
            }
            if self.sum_d[i] < summin {
                indices[2] = i as i16;
                summin = self.sum_d[i];
            }
        }
        let mut summax2nd = -2147483647i32;
        indices[1] = -1;
        for i in 0..DTX_HIST_SIZE {
            if self.sum_d[i] > summax2nd && i as i16 != indices[0] {
                indices[1] = i as i16;
                summax2nd = self.sum_d[i];
            }
        }
        for ix in indices.iter_mut() {
            *ix = sub(self.hist_ptr, *ix);
            if *ix < 0 {
                *ix = add(*ix, DTX_HIST_SIZE as i16);
            }
        }
        let tmp = norm_l(summax);
        let summax = l_shl(summax, tmp);
        let summin = l_shl(summin, tmp);
        if l_mult(round(summax), INV_MED_THRESH) <= summin {
            indices[0] = -1;
        }
        let summax2nd = l_shl(summax2nd, tmp);
        if l_mult(round(summax2nd), INV_MED_THRESH) <= summin {
            indices[1] = -1;
        }
    }

    fn dithering_control(&self) -> i16 {
        let mut isf_diff = 0i32;
        for &s in &self.sum_d {
            isf_diff = l_add(isf_diff, s);
        }
        let mut cn_dith = (l_shr(isf_diff, 26) > 0) as i16;
        let mut mean = 0i16;
        for &h in &self.log_en_hist {
            mean = add(mean, h);
        }
        let mean = shr(mean, 3);
        let mut gain_diff = 0i16;
        for &h in &self.log_en_hist {
            gain_diff = add(gain_diff, abs_s(sub(h, mean)));
        }
        if gain_diff > GAIN_THR {
            cn_dith = 1;
        }
        cn_dith
    }
}

/// `aver_isf_history`: Σ of the 8 ISF vectors, the (up to two) outliers in `indices[0..2]`
/// replaced by the median one `indices[2]` (-1: none).
fn aver_isf_history(
    isf_old: &mut [i16; M * DTX_HIST_SIZE],
    indices: &[i16; 3],
    isf_aver: &mut [i32; M],
) {
    let mut isf_tmp = [0i16; 2 * M];
    for k in 0..2 {
        if indices[k] != -1 {
            let (a, b) = (indices[k] as usize * M, indices[2] as usize * M);
            for i in 0..M {
                isf_tmp[k * M + i] = isf_old[a + i];
                isf_old[a + i] = isf_old[b + i];
            }
        }
    }
    for j in 0..M {
        let mut l = 0i32;
        for i in 0..DTX_HIST_SIZE {
            l = l_add(l, isf_old[i * M + j] as i32);
        }
        isf_aver[j] = l;
    }
    for k in 0..2 {
        if indices[k] != -1 {
            let a = indices[k] as usize * M;
            isf_old[a..a + M].copy_from_slice(&isf_tmp[k * M..k * M + M]);
        }
    }
}

#[cfg(test)]
#[allow(unsafe_code)]
mod vs_ref {
    use super::*;
    use crate::testutil::{Rng, with_ref};

    unsafe extern "C" {
        fn refwb_dtx_enc_init(st: *mut *mut DtxEncState, isf_init: *const i16) -> i16;
        fn refwb_dtx_enc_exit(st: *mut *mut DtxEncState);
        fn refwb_dtx_enc(
            st: *mut DtxEncState,
            isf: *mut i16,
            exc2: *mut i16,
            prms: *mut *mut i16,
        ) -> i16;
        fn refwb_dtx_buffer(st: *mut DtxEncState, isf: *mut i16, enr: i32, mode: i16) -> i16;
        fn refwb_tx_dtx_handler(st: *mut DtxEncState, vad: i16, mode: *mut i16);
    }

    /// ISFs near each other (the median/outlier logic) or spread, energies from silence to
    /// loud, VAD patterns driving the hangover; SID encodings at random points. Whole state,
    /// SID bits, CN ISFs and excitation compared.
    #[test]
    fn dtx_enc_matches_reference() {
        let mut r = Rng::new(61);
        for seq in 0..40u64 {
            let isf_init = crate::enc_tables::ISF_INIT;
            let mut st = DtxEncState::new(&isf_init);
            let mut c: *mut DtxEncState = std::ptr::null_mut();
            // SAFETY: init allocates the state; sumD[7] is uninitialised in the C (never read
            // before written): aligned here for the state comparison.
            with_ref(|| unsafe {
                refwb_dtx_enc_init(&mut c, isf_init.as_ptr());
                (*c).sum_d[7] = st.sum_d[7];
            });
            let mut base: [i16; M] = std::array::from_fn(|i| 400 + 800 * i as i16);
            for call in 0..300u64 {
                let spread = if seq % 2 == 0 { 30 } else { 600 };
                if r.u64().is_multiple_of(10) {
                    base = std::array::from_fn(|i| 300 + 900 * i as i16 + (r.u64() % 200) as i16);
                }
                let mut isf: [i16; M] = std::array::from_fn(|i| {
                    base[i] + (r.u64() % spread) as i16 - (spread / 2) as i16
                });
                let enr = match r.u64() % 4 {
                    0 => 0,
                    1 => r.i32().max(1),
                    _ => (r.u64() % 1_000_000) as i32,
                };
                let mode = (r.u64() % 9) as i16;
                let vad = ((call / (3 + seq % 11)) % 3 == 0) as i16;
                let (mut m, mut mc) = (mode, mode);
                st.tx_handler(vad, &mut m);
                st.buffer(&isf, enr, mode);
                let sid = r.u64().is_multiple_of(4);
                let (mut bits, mut cbits) = ([0i16; 35], [0i16; 35]);
                let (mut exc, mut excc) = ([0i16; L_FRAME12K8], [0i16; L_FRAME12K8]);
                let mut isfc = isf;
                if sid {
                    let mut s = Serial::new(&mut bits);
                    st.enc(&mut isf, &mut exc, &mut s);
                }
                // SAFETY: a live state; 16 ISFs, 256 samples, 35 SID bits.
                with_ref(|| unsafe {
                    refwb_tx_dtx_handler(c, vad, &mut mc);
                    refwb_dtx_buffer(c, isfc.as_mut_ptr(), enr, mode);
                    if sid {
                        let mut p = cbits.as_mut_ptr();
                        refwb_dtx_enc(c, isfc.as_mut_ptr(), excc.as_mut_ptr(), &mut p);
                    }
                    assert_eq!(m, mc, "tx_dtx_handler {seq} {call}");
                    assert_eq!(
                        (isf, bits, exc),
                        (isfc, cbits, excc),
                        "dtx_enc {seq} {call}"
                    );
                    assert_eq!(st, *c, "dtx_encState {seq} {call}");
                });
            }
            // SAFETY: frees the state from init.
            with_ref(|| unsafe { refwb_dtx_enc_exit(&mut c) });
        }
    }
}
