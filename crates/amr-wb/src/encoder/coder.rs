//! The encoder's frame (vo-amrwbenc `voAMRWBEnc.c`: `Reset_encoder`, `coder`, `synthesis`),
//! with the reference's saturation: pre-processing, LP analysis, open-loop pitch, VAD/DTX,
//! then per subframe the closed-loop pitch, the algebraic codebook, the gains and (23.85
//! kbit/s) the high-band gain.

use super::acelp2::{acelp_2t64_fx, cor_h_x};
use super::acelp4::acelp_4t64_fx;
use super::bits::Serial;
use super::dtx::{DtxEncState, MRDTX};
use super::filters::*;
use super::lpc::*;
use super::ops::*;
use super::pitch::*;
use super::pitch_ol::*;
use super::preprocess::{NB_COEF_DOWN, decim_12k8};
use super::q_gain::{Q_GAIN2_INIT, q_gain2};
use super::vad::VadState;
use crate::basic_op::dot_product12;
use crate::cnst::*;
use crate::enc_tables::{HP_GAIN, INTERPOL_FRAC, ISF_INIT, ISP_INIT};
use crate::excitation::{noise_gen_amrwb as random, pit_shrp, pred_lt4, voice_factor};
use crate::lpc::{interpolate_isp as int_isp, isf_isp, isp_az, weight_amrwb_lpc as weight_a};
use crate::synthesis::{
    band_pass_6k_7k as filt_6k_7k, deemphasis_32, highpass_50hz_at_12k8, syn_filt_32, wb_syn_filt,
};

const L_NEXT: usize = 64;
const L_TOTAL: usize = 384;
const L_FILT: usize = 12;
const L_FILT16K: usize = NB_COEF_DOWN;
const OPL_DECIM: usize = 2;
const GAMMA1: i16 = 30147;
const TILT_FAC: i16 = 22282;
const GP_CLIP: i16 = 15565;
const MP1: usize = M + 1;
const OLD_WSP: usize = PIT_MAX / OPL_DECIM;
const OLD_EXC: usize = PIT_MAX + L_INTERPOL;

/// `nb_of_bits`: serial bits per mode (`MRDTX`: a SID frame).
pub(crate) const NB_OF_BITS: [i16; 10] = [
    NBBITS_7K, NBBITS_9K, NBBITS_12K, NBBITS_14K, NBBITS_16K, NBBITS_18K, NBBITS_20K, NBBITS_23K,
    NBBITS_24K, 35,
];

/// `Coder_State`.
#[derive(Clone)]
pub(crate) struct CoderState {
    mem_decim: [i16; 2 * L_FILT16K],
    mem_sig_in: [i16; 6],
    mem_preemph: i16,
    old_speech: [i16; L_TOTAL - L_FRAME],
    old_wsp: [i16; OLD_WSP],
    old_exc: [i16; OLD_EXC],
    mem_levinson: [i16; M + 2],
    ispold: [i16; M],
    ispold_q: [i16; M],
    past_isfq: [i16; M],
    mem_wsp: i16,
    mem_decim2: [i16; 3],
    mem_w0: i16,
    mem_syn: [i16; M],
    tilt_code: i16,
    old_wsp_max: i16,
    old_wsp_shift: i16,
    q_old: i16,
    q_max: [i16; 2],
    gp_clip: [i16; 2],
    qua_gain: [i16; 4],
    old_t0_med: i16,
    ol_gain: i16,
    ada_w: i16,
    ol_wght_flg: i16,
    old_ol_lag: [i16; 5],
    hp_wsp_mem: [i16; 9],
    old_hp_wsp: [i16; L_FRAME / OPL_DECIM + OLD_WSP],
    vad: VadState,
    dtx: DtxEncState,
    first_frame: bool,
    isfold: [i16; M],
    l_gc_thres: i32,
    mem_syn_hi: [i16; M],
    mem_syn_lo: [i16; M],
    mem_deemph: i16,
    mem_sig_out: [i16; 6],
    mem_hp400: [i16; 6],
    mem_syn_hf: [i16; M],
    mem_hf: [i16; 2 * L_FILT16K],
    mem_hf2: [i16; 2 * L_FILT16K],
    seed2: i16,
    vad_hist: i16,
    gain_alpha: i16,
}

impl CoderState {
    pub(crate) fn new() -> Self {
        let mut s = Self {
            mem_decim: [0; 2 * L_FILT16K],
            mem_sig_in: [0; 6],
            mem_preemph: 0,
            old_speech: [0; L_TOTAL - L_FRAME],
            old_wsp: [0; OLD_WSP],
            old_exc: [0; OLD_EXC],
            mem_levinson: [0; M + 2],
            ispold: ISP_INIT,
            ispold_q: ISP_INIT,
            past_isfq: [0; M],
            mem_wsp: 0,
            mem_decim2: [0; 3],
            mem_w0: 0,
            mem_syn: [0; M],
            tilt_code: 0,
            old_wsp_max: 0,
            old_wsp_shift: 0,
            q_old: 15,
            q_max: [15; 2],
            gp_clip: GP_CLIP_INIT,
            qua_gain: Q_GAIN2_INIT,
            old_t0_med: 40,
            ol_gain: 0,
            ada_w: 0,
            ol_wght_flg: 0,
            old_ol_lag: [40; 5],
            hp_wsp_mem: [0; 9],
            old_hp_wsp: [0; L_FRAME / OPL_DECIM + OLD_WSP],
            vad: VadState::new(),
            dtx: DtxEncState::new(&ISF_INIT),
            first_frame: true,
            isfold: ISF_INIT,
            l_gc_thres: 0,
            mem_syn_hi: [0; M],
            mem_syn_lo: [0; M],
            mem_deemph: 0,
            mem_sig_out: [0; 6],
            mem_hp400: [0; 6],
            mem_syn_hf: [0; M],
            mem_hf: [0; 2 * L_FILT16K],
            mem_hf2: [0; 2 * L_FILT16K],
            seed2: RANDOM_INITSEED,
            vad_hist: 0,
            gain_alpha: 32767,
        };
        s.reset(true);
        s
    }

    /// `Reset_encoder`: the excitation, synthesis and quantiser memories; with `reset_all`
    /// everything (as a fresh encoder).
    pub(crate) fn reset(&mut self, reset_all: bool) {
        self.old_exc = [0; OLD_EXC];
        self.mem_syn = [0; M];
        self.past_isfq = [0; M];
        self.mem_w0 = 0;
        self.tilt_code = 0;
        self.first_frame = true;
        self.gp_clip = GP_CLIP_INIT;
        self.l_gc_thres = 0;
        if reset_all {
            self.old_speech = [0; L_TOTAL - L_FRAME];
            self.old_wsp = [0; OLD_WSP];
            self.mem_decim2 = [0; 3];
            self.mem_decim = [0; 2 * L_FILT16K];
            self.mem_sig_in = [0; 6];
            self.mem_levinson = [0; M + 2];
            self.qua_gain = Q_GAIN2_INIT;
            self.hp_wsp_mem = [0; 9];
            self.ispold = ISP_INIT;
            self.ispold_q = ISP_INIT;
            self.mem_preemph = 0;
            self.mem_wsp = 0;
            self.q_old = 15;
            self.q_max = [15; 2];
            self.old_wsp_max = 0;
            self.old_wsp_shift = 0;
            self.old_t0_med = 40;
            self.ol_gain = 0;
            self.ada_w = 0;
            self.ol_wght_flg = 0;
            self.old_ol_lag = [40; 5];
            self.old_hp_wsp = [0; L_FRAME / OPL_DECIM + OLD_WSP];
            self.mem_syn_hf = [0; M];
            self.mem_syn_hi = [0; M];
            self.mem_syn_lo = [0; M];
            self.mem_sig_out = [0; 6];
            self.mem_hf = [0; 2 * L_FILT16K];
            self.mem_hp400 = [0; 6];
            self.isfold = ISF_INIT;
            self.mem_deemph = 0;
            self.seed2 = RANDOM_INITSEED;
            self.mem_hf2 = [0; 2 * L_FILT16K];
            self.gain_alpha = 32767;
            self.vad_hist = 0;
            self.vad.reset();
            self.dtx.reset(&ISF_INIT);
        }
    }

    /// The open-loop lag of half (or, at 6.60, all) the frame of `wsp` from `at`, and the
    /// lag-weighting state after it.
    fn open_loop(&mut self, wsp: &[i16], at: usize, l_frame: usize) -> i16 {
        let t = pitch_med_ol(
            wsp,
            at,
            (PIT_MIN as usize / OPL_DECIM) as i16,
            OLD_WSP as i16,
            l_frame,
            self.old_t0_med,
            &mut self.ol_gain,
            &mut self.hp_wsp_mem,
            &mut self.old_hp_wsp,
            self.ol_wght_flg,
        );
        if self.ol_gain > 19661 {
            self.old_t0_med = med_olag(t, &mut self.old_ol_lag);
            self.ada_w = 32767;
        } else {
            self.ada_w = mult(self.ada_w, 29491);
        }
        self.ol_wght_flg = (self.ada_w >= 26214) as i16;
        self.vad.tone_detection(self.ol_gain);
        t * OPL_DECIM as i16
    }

    /// `coder`: one 20 ms frame at `mode` (0..=8); with `allow_dtx` the DTX handler may turn
    /// it into `MRDTX` (a SID frame). Writes the serial bits (±127) to `prms`; returns the mode
    /// used and the number of bits.
    pub(crate) fn coder(
        &mut self,
        mode: i16,
        speech16k: &[i16; L_FRAME16K],
        prms: &mut [i16; MAX_SERIAL_SIZE],
        allow_dtx: bool,
    ) -> (i16, usize) {
        let mut mode = mode;
        let mut ser_size = NB_OF_BITS[mode as usize];
        let codec_mode = mode;
        let mut prms = Serial::new(prms);
        // old_speech: [past (L_TOTAL - L_FRAME), new frame]; speech, new_speech, the window
        let mut old_speech = [0i16; L_TOTAL];
        const NEW_SPEECH: usize = L_TOTAL - L_FRAME - L_FILT;
        const SPEECH: usize = L_TOTAL - L_FRAME - L_NEXT;
        let mut old_wsp = [0i16; L_FRAME + OLD_WSP];
        const WSP: usize = OLD_WSP;
        let mut old_exc = [0i16; (L_FRAME + 1) + OLD_EXC];
        const EXC: usize = OLD_EXC;
        old_speech[..L_TOTAL - L_FRAME].copy_from_slice(&self.old_speech);
        old_wsp[..OLD_WSP].copy_from_slice(&self.old_wsp);
        old_exc[..OLD_EXC].copy_from_slice(&self.old_exc);

        // 16 kHz → 12.8 kHz, then L_FILT lookahead samples from zeros (autocorrelation only)
        decim_12k8(
            speech16k,
            &mut old_speech[NEW_SPEECH..NEW_SPEECH + L_FRAME],
            &mut self.mem_decim,
        );
        let mut mem = self.mem_decim;
        decim_12k8(
            &[0; L_FILT16K],
            &mut old_speech[NEW_SPEECH + L_FRAME..],
            &mut mem,
        );
        highpass_50hz_at_12k8(
            &mut old_speech[NEW_SPEECH..NEW_SPEECH + L_FRAME],
            &mut self.mem_sig_in,
        );
        let mut mem = self.mem_sig_in;
        highpass_50hz_at_12k8(&mut old_speech[NEW_SPEECH + L_FRAME..], &mut mem);

        // preemphasis, scaled for precision (Q_new) by the frame's maximum
        let mu = PREEMPH_FAC >> 1;
        let ns = NEW_SPEECH;
        let pre = |x: i16, prev: i16| l_msu(l_mult(x, 16384), prev, mu);
        let mut l_max = l_abs(pre(old_speech[ns], self.mem_preemph));
        for i in 1..L_FRAME + L_FILT {
            let l = l_abs(pre(old_speech[ns + i], old_speech[ns + i - 1]));
            if l > l_max {
                l_max = l;
            }
        }
        let tmp = extract_h(l_max);
        let shift = if tmp == 0 {
            Q_MAX
        } else {
            sub(norm_s(tmp), 1).clamp(0, Q_MAX)
        };
        let q_new = shift.min(self.q_max[0]).min(self.q_max[1]);
        let mut exp = sub(q_new, self.q_old);
        self.q_old = q_new;
        self.q_max[1] = self.q_max[0];
        self.q_max[0] = shift;
        let tmp = old_speech[ns + L_FRAME - 1];
        for i in (1..L_FRAME + L_FILT).rev() {
            old_speech[ns + i] = round(l_shl(
                pre(old_speech[ns + i], old_speech[ns + i - 1]),
                q_new,
            ));
        }
        old_speech[ns] = round(l_shl(pre(old_speech[ns], self.mem_preemph), q_new));
        self.mem_preemph = tmp;
        scale_sig(&mut old_speech[..L_TOTAL - L_FRAME - L_FILT], exp);
        scale_sig(&mut old_exc[..OLD_EXC], exp);
        scale_sig(&mut self.mem_syn, exp);
        scale_sig(&mut self.mem_decim2, exp);
        scale_sig(std::slice::from_mut(&mut self.mem_wsp), exp);
        scale_sig(std::slice::from_mut(&mut self.mem_w0), exp);

        // VAD (on the new speech, slightly ahead), DTX
        let mut buf: [i16; L_FRAME] = old_speech[ns..ns + L_FRAME].try_into().unwrap();
        scale_sig(&mut buf, sub(1, q_new));
        let vad_flag = self.vad.vad(&buf);
        self.vad_hist = if vad_flag == 0 {
            add(self.vad_hist, 1)
        } else {
            0
        };
        if allow_dtx {
            self.dtx.tx_handler(vad_flag, &mut mode);
            ser_size = NB_OF_BITS[mode as usize];
        }
        if mode != MRDTX {
            prms.parm(vad_flag, 1);
        }

        // LP analysis (window centred on the 4th subframe), interpolated A(z) per subframe
        let (mut r_h, mut r_l) = ([0i16; M + 1], [0i16; M + 1]);
        let mut rc = [0i16; M];
        let mut a = [0i16; NB_SUBFR * MP1];
        let mut ispnew = [0i16; M];
        autocorr(&old_speech[..L_WINDOW], &mut r_h, &mut r_l);
        lag_window(&mut r_h, &mut r_l);
        levinson(&r_h, &r_l, &mut a, &mut rc, &mut self.mem_levinson);
        az_isp(&a, &mut ispnew, &self.ispold);
        int_isp(&self.ispold, &ispnew, &INTERPOL_FRAC, &mut a);
        self.ispold = ispnew;
        let mut isf = [0i16; M];
        isp_isf(&ispnew, &mut isf);
        gp_clip_test_isf(ser_size, &isf, &mut self.gp_clip);

        // weighted speech, decimated by 2 and scaled to 12 bits; open-loop pitch
        let mut ap = [0i16; MP1];
        for (k, p_a) in a.as_chunks::<MP1>().0.iter().enumerate() {
            weight_a(p_a, &mut ap, GAMMA1);
            let at = SPEECH + k * L_SUBFR;
            residu(
                &ap,
                &old_speech,
                at,
                &mut old_wsp[WSP + k * L_SUBFR..],
                L_SUBFR,
            );
        }
        deemph2(
            &mut old_wsp[WSP..WSP + L_FRAME],
            TILT_FAC,
            &mut self.mem_wsp,
        );
        let mut max = 0i16;
        for &w in &old_wsp[WSP..WSP + L_FRAME] {
            max = max.max(abs_s(w));
        }
        let tmp = max.max(self.old_wsp_max);
        self.old_wsp_max = max;
        let shift = sub(norm_s(tmp), 3).min(0);
        lp_decim2(&mut old_wsp[WSP..], L_FRAME, &mut self.mem_decim2);
        scale_sig(&mut old_wsp[WSP..WSP + L_FRAME / OPL_DECIM], shift);
        exp = add(exp, sub(shift, self.old_wsp_shift));
        self.old_wsp_shift = shift;
        scale_sig(&mut old_wsp[..OLD_WSP], exp);
        scale_sig(&mut self.old_hp_wsp[..OLD_WSP], exp);
        scale_mem_hp_wsp(&mut self.hp_wsp_mem, exp);
        let half = L_FRAME / 2 / OPL_DECIM;
        let t_op;
        let t_op2;
        if ser_size == NBBITS_7K {
            t_op = self.open_loop(&old_wsp, WSP, L_FRAME / OPL_DECIM);
            t_op2 = t_op;
        } else {
            t_op = self.open_loop(&old_wsp, WSP, half);
            t_op2 = self.open_loop(&old_wsp, WSP + half, half);
        }

        let mut exc2 = [0i16; L_FRAME];
        if mode == MRDTX {
            // comfort noise: buffer the ISFs and energy, send the SID parameters
            residu(
                &a[3 * MP1..],
                &old_speech,
                SPEECH,
                &mut old_exc[EXC..],
                L_FRAME,
            );
            for i in 0..L_FRAME {
                exc2[i] = shr(old_exc[EXC + i], q_new);
            }
            let mut l = 0i32;
            for &e in &exc2 {
                l = l_mac(l, e, e);
            }
            self.dtx.buffer(&isf, l_shr(l, 1), codec_mode);
            self.dtx.enc(&mut isf, &mut exc2, &mut prms);
            let mut ispnew_q = [0i16; M];
            isf_isp(&isf, &mut ispnew_q);
            let mut aq = [0i16; MP1];
            isp_az(&ispnew_q, &mut aq, false);
            for k in 0..NB_SUBFR {
                let at = k * L_SUBFR;
                self.synthesis(
                    &aq,
                    &mut exc2[at..at + L_SUBFR],
                    0,
                    &speech16k[at * 5 / 4..],
                );
            }
            self.isfold = isf;
            self.reset(false);
            self.old_speech.copy_from_slice(&old_speech[L_FRAME..]);
            self.old_wsp
                .copy_from_slice(&old_wsp[L_FRAME / OPL_DECIM..L_FRAME / OPL_DECIM + OLD_WSP]);
            return (mode, prms.len());
        }

        // ISF quantisation
        let isf_in = isf;
        if ser_size <= NBBITS_7K {
            let mut indice = [0i16; 5];
            qpisf_2s_36b(&isf_in, &mut isf, &mut self.past_isfq, &mut indice, 4);
            for (k, n) in [8, 8, 7, 7, 6].into_iter().enumerate() {
                prms.parm(indice[k], n);
            }
        } else {
            let mut indice = [0i16; 7];
            qpisf_2s_46b(&isf_in, &mut isf, &mut self.past_isfq, &mut indice, 4);
            for (k, n) in [8, 8, 6, 7, 7, 5, 5].into_iter().enumerate() {
                prms.parm(indice[k], n);
            }
        }
        // stability of the ISFs (for the noise enhancer)
        let mut l = 0i32;
        for i in 0..M - 1 {
            let t = sub(isf[i], self.isfold[i]);
            l = l_mac(l, t, t);
        }
        let tmp = mult(extract_h(l_shl(l, 8)), 26214);
        let stab_fac = shl(sub(20480, tmp), 1).max(0);
        self.isfold = isf;
        let mut ispnew_q = [0i16; M];
        isf_isp(&isf, &mut ispnew_q);
        if self.first_frame {
            self.first_frame = false;
            self.ispold_q = ispnew_q;
        }
        let mut aq = [0i16; NB_SUBFR * MP1];
        int_isp(&self.ispold_q, &ispnew_q, &INTERPOL_FRAC, &mut aq);
        self.ispold_q = ispnew_q;
        for (k, p_aq) in aq.as_chunks::<MP1>().0.iter().enumerate() {
            let at = k * L_SUBFR;
            residu(
                p_aq,
                &old_speech,
                SPEECH + at,
                &mut old_exc[EXC + at..],
                L_SUBFR,
            );
        }
        if vad_flag == 0 {
            for i in 0..L_FRAME {
                exc2[i] = shr(old_exc[EXC + i], q_new);
            }
            let mut l = 0i32;
            for &e in &exc2 {
                l = l_mac(l, e, e);
            }
            self.dtx.buffer(&isf, l_shr(l, 1), codec_mode);
        }

        let range = |t: i16| {
            let t0_min = sub(t, 8).max(PIT_MIN);
            let t0_max = add(t0_min, 15);
            if t0_max > PIT_MAX as i16 {
                (PIT_MAX as i16 - 15, PIT_MAX as i16)
            } else {
                (t0_min, t0_max)
            }
        };
        let (mut t0_min, mut t0_max) = range(t_op);
        for k in 0..NB_SUBFR {
            let i_subfr = k * L_SUBFR;
            let p_a = &a[k * MP1..(k + 1) * MP1];
            let p_aq = &aq[k * MP1..(k + 1) * MP1];
            let mut pit_flag = i_subfr;
            if i_subfr == 2 * L_SUBFR && ser_size > NBBITS_7K {
                pit_flag = 0;
                (t0_min, t0_max) = range(t_op2);
            }
            // target for the pitch search: the weighted error of the zero-input response
            let mut error = [0i16; M + L_SUBFR];
            for i in 0..M {
                error[i] = sub(old_speech[SPEECH + i + i_subfr - M], self.mem_syn[i]);
            }
            residu(
                p_aq,
                &old_speech,
                SPEECH + i_subfr,
                &mut old_exc[EXC + i_subfr..],
                L_SUBFR,
            );
            {
                let (mem, out) = error.split_at_mut(M);
                let mut m = [0i16; M];
                m.copy_from_slice(mem);
                syn_filt(
                    p_aq,
                    Some(&old_exc[EXC + i_subfr..EXC + i_subfr + L_SUBFR]),
                    out,
                    &mut m,
                    false,
                );
            }
            weight_a(p_a, &mut ap, GAMMA1);
            let mut xn = [0i16; L_SUBFR];
            residu(&ap, &error, M, &mut xn, L_SUBFR);
            deemph2(&mut xn, TILT_FAC, &mut self.mem_w0);

            // target in the residual domain for the innovation search
            let mut code = [0i16; M + L_SUBFR];
            code[M..M + L_SUBFR / 2].copy_from_slice(&xn[..L_SUBFR / 2]);
            let mut t = 0;
            preemph2(&mut code[M..M + L_SUBFR / 2], TILT_FAC, &mut t);
            weight_a(p_a, &mut ap, GAMMA1);
            {
                let (mem, out) = code.split_at_mut(M);
                let mut m = [0i16; M];
                m.copy_from_slice(mem);
                syn_filt(&ap, None, &mut out[..L_SUBFR / 2], &mut m, false);
            }
            let mut cn = [0i16; L_SUBFR];
            residu(p_aq, &code, M, &mut cn, L_SUBFR / 2);
            cn[L_SUBFR / 2..]
                .copy_from_slice(&old_exc[EXC + i_subfr + L_SUBFR / 2..EXC + i_subfr + L_SUBFR]);

            // impulse response of the weighted synthesis filter
            let mut error = [0i16; M + L_SUBFR];
            weight_a(p_a, &mut error[M..], GAMMA1);
            let mut h1 = [0i16; L_SUBFR];
            for i in 0..L_SUBFR {
                let mut l = l_mult(error[i + M], 16384);
                for j in 1..=M {
                    l = l_msu(l, p_aq[j], error[i + M - j]);
                }
                h1[i] = round(l_shl(l, 3));
                error[i + M] = h1[i];
            }
            let mut t = 0;
            deemph2(&mut h1, TILT_FAC, &mut t);
            let mut h2 = h1;
            scale_sig(&mut h2, -2);
            scale_sig(&mut xn, shift);
            scale_sig(&mut h1, add(1, shift));

            // closed-loop fractional pitch and its index
            let mut t0_frac = 0;
            let mut t0;
            if ser_size <= NBBITS_9K {
                t0 = pitch_fr4(
                    &old_exc,
                    EXC + i_subfr,
                    &xn,
                    &h1,
                    t0_min,
                    t0_max,
                    &mut t0_frac,
                    pit_flag,
                    PIT_MIN,
                    PIT_FR1_8B,
                    L_SUBFR,
                );
                if pit_flag == 0 {
                    let index = if t0 < PIT_FR1_8B {
                        sub(add(shl(t0, 1), shr(t0_frac, 1)), PIT_MIN * 2)
                    } else {
                        add(sub(t0, PIT_FR1_8B), (PIT_FR1_8B - PIT_MIN) * 2)
                    };
                    prms.parm(index, 8);
                    (t0_min, t0_max) = range(t0);
                } else {
                    let i = sub(t0, t0_min);
                    prms.parm(add(shl(i, 1), shr(t0_frac, 1)), 5);
                }
            } else {
                t0 = pitch_fr4(
                    &old_exc,
                    EXC + i_subfr,
                    &xn,
                    &h1,
                    t0_min,
                    t0_max,
                    &mut t0_frac,
                    pit_flag,
                    PIT_FR2,
                    PIT_FR1_9B,
                    L_SUBFR,
                );
                if pit_flag == 0 {
                    let index = if t0 < PIT_FR2 {
                        sub(add(shl(t0, 2), t0_frac), PIT_MIN * 4)
                    } else if t0 < PIT_FR1_9B {
                        add(
                            sub(add(shl(t0, 1), shr(t0_frac, 1)), PIT_FR2 * 2),
                            (PIT_FR2 - PIT_MIN) * 4,
                        )
                    } else {
                        add(
                            add(sub(t0, PIT_FR1_9B), (PIT_FR2 - PIT_MIN) * 4),
                            (PIT_FR1_9B - PIT_FR2) * 2,
                        )
                    };
                    prms.parm(index, 9);
                    (t0_min, t0_max) = range(t0);
                } else {
                    let i = sub(t0, t0_min);
                    prms.parm(add(shl(i, 2), t0_frac), 6);
                }
            }
            let clip_gain = gp_clip(ser_size, &self.gp_clip);

            // adaptive codebook vector (unfiltered, then low-passed) and its gains
            pred_lt4(&mut old_exc, EXC + i_subfr, t0, t0_frac, L_SUBFR + 1);
            let mut y1 = [0i16; L_SUBFR];
            let mut g_coeff = [0i16; 4];
            let mut dn = [0i16; L_SUBFR];
            let gain1 = if ser_size > NBBITS_9K {
                convolve(&old_exc[EXC + i_subfr..], &h1, &mut y1, L_SUBFR);
                let mut g = g_pitch(&xn, &y1, &mut g_coeff);
                if clip_gain && g > GP_CLIP {
                    g = GP_CLIP;
                }
                updt_tar(&xn, &mut dn, &y1, g, L_SUBFR);
                g
            } else {
                0
            };
            let mut code = [0i16; L_SUBFR];
            for i in 0..L_SUBFR {
                let e = EXC + i_subfr + i;
                let mut l = l_mult(5898, old_exc[e - 1]);
                l = l_mac(l, 20972, old_exc[e]);
                l = l_mac(l, 5898, old_exc[e + 1]);
                code[i] = round(l);
            }
            let mut y2 = [0i16; L_SUBFR];
            convolve(&code, &h1, &mut y2, L_SUBFR);
            let mut g_coeff2 = [0i16; 4];
            let mut gain2 = g_pitch(&xn, &y2, &mut g_coeff2);
            if clip_gain && gain2 > GP_CLIP {
                gain2 = GP_CLIP;
            }
            let mut xn2 = [0i16; L_SUBFR];
            updt_tar(&xn, &mut xn2, &y2, gain2, L_SUBFR);
            let mut select = false;
            if ser_size > NBBITS_9K {
                let mut l = 0i32;
                for &d in &dn {
                    l = l_mac(l, d, d);
                }
                for &x in &xn2 {
                    l = l_msu(l, x, x);
                }
                select = l <= 0;
                prms.parm(select as i16, 1);
            }
            let mut gain_pit;
            if !select {
                gain_pit = gain2;
                old_exc[EXC + i_subfr..EXC + i_subfr + L_SUBFR].copy_from_slice(&code);
                y1 = y2;
                g_coeff = g_coeff2;
            } else {
                gain_pit = gain1;
                xn2 = dn;
            }
            let cn_in = cn;
            updt_tar(
                &cn_in,
                &mut cn,
                &old_exc[EXC + i_subfr..],
                gain_pit,
                L_SUBFR,
            );
            scale_sig(&mut cn, shift);

            // the fixed-gain pitch contribution in the impulse response, the innovation
            let mut t = 0;
            preemph(&mut h2, self.tilt_code, &mut t);
            if t0_frac > 2 {
                t0 = add(t0, 1);
            }
            pit_shrp(&mut h2, t0, PIT_SHARP);
            cor_h_x(&h2, &xn2, &mut dn);
            let mut indice = [0i16; 8];
            // (y2 still holds the filtered pitch vector: a search that never improves on its
            // start leaves it there, as the C does)
            if ser_size <= NBBITS_7K {
                acelp_2t64_fx(&mut dn, &cn, &h2, &mut code, &mut y2, &mut indice[0]);
                prms.parm(indice[0], 12);
            } else {
                let (nbbits, widths): (i16, &[usize]) = match ser_size {
                    s if s <= NBBITS_9K => (20, &[5, 5, 5, 5]),
                    s if s <= NBBITS_12K => (36, &[9, 9, 9, 9]),
                    s if s <= NBBITS_14K => (44, &[13, 13, 9, 9]),
                    s if s <= NBBITS_16K => (52, &[13, 13, 13, 13]),
                    s if s <= NBBITS_18K => (64, &[2, 2, 2, 2, 14, 14, 14, 14]),
                    s if s <= NBBITS_20K => (72, &[10, 10, 2, 2, 10, 10, 14, 14]),
                    _ => (88, &[11; 8]),
                };
                acelp_4t64_fx(
                    &mut dn,
                    &cn,
                    &h2,
                    &mut code,
                    &mut y2,
                    nbbits,
                    ser_size,
                    &mut indice,
                );
                for (k, &w) in widths.iter().enumerate() {
                    prms.parm(indice[k], w);
                }
            }
            let mut t = 0;
            preemph(&mut code, self.tilt_code, &mut t);
            pit_shrp(&mut code, t0, PIT_SHARP);

            // gains
            let mut l_gain_code = 0i32;
            let nbits = if ser_size <= NBBITS_9K { 6 } else { 7 };
            let index = q_gain2(
                &xn,
                &y1,
                add(q_new, shift),
                &y2,
                &code,
                &g_coeff,
                nbits,
                &mut gain_pit,
                &mut l_gain_code,
                clip_gain,
                &mut self.qua_gain,
            );
            prms.parm(index, nbits as usize);
            gp_clip_test_gain_pit(ser_size, gain_pit, &mut self.gp_clip);
            let mut gain_code = round(l_shl(l_gain_code, q_new));

            // voicing → tilt of the next code; target memory; excitation and synthesis
            let mut exc2 = [0i16; L_SUBFR];
            exc2.copy_from_slice(&old_exc[EXC + i_subfr..EXC + i_subfr + L_SUBFR]);
            scale_sig(&mut exc2, shift);
            let voice_fac = voice_factor(&exc2, shift, gain_pit, &code, gain_code);
            self.tilt_code = add(shr(voice_fac, 2), 8192);
            let mut l = l_mult(gain_code, y2[L_SUBFR - 1]);
            l = l_shl(l, add(5, shift));
            l = l_negate(l);
            l = l_mac(l, xn[L_SUBFR - 1], 16384);
            l = l_msu(l, y1[L_SUBFR - 1], gain_pit);
            l = l_shl(l, sub(1, shift));
            self.mem_w0 = round(l);
            if ser_size >= NBBITS_24K {
                exc2.copy_from_slice(&old_exc[EXC + i_subfr..EXC + i_subfr + L_SUBFR]);
            }
            for i in 0..L_SUBFR {
                let e = &mut old_exc[EXC + i_subfr + i];
                let mut l = l_shl(l_mult(gain_code, code[i]), 5);
                l = l_mac(l, *e, gain_pit);
                *e = round(l_shl(l, 1));
            }
            let mut synth = [0i16; L_SUBFR];
            syn_filt(
                p_aq,
                Some(&old_exc[EXC + i_subfr..EXC + i_subfr + L_SUBFR]),
                &mut synth,
                &mut self.mem_syn,
                true,
            );

            if ser_size >= NBBITS_24K {
                // noise enhancer (code gain smoothed towards a threshold on stable, unvoiced
                // frames) and pitch enhancer, then the high band's gain
                let (gc, gc_lo) = l_extract(l_gain_code);
                let tmp = sub(16384, shr(voice_fac, 1));
                let fac = mult(stab_fac, tmp);
                let mut l = l_gain_code;
                if l < self.l_gc_thres {
                    l = l_add(l, mpy_32_16(gc, gc_lo, 6226));
                    if l > self.l_gc_thres {
                        l = self.l_gc_thres;
                    }
                } else {
                    l = mpy_32_16(gc, gc_lo, 27536);
                    if l < self.l_gc_thres {
                        l = self.l_gc_thres;
                    }
                }
                self.l_gc_thres = l;
                l_gain_code = mpy_32_16(gc, gc_lo, sub(32767, fac));
                let (gt, gt_lo) = l_extract(l);
                l_gain_code = l_add(l_gain_code, mpy_32_16(gt, gt_lo, fac));
                let tmp = add(shr(voice_fac, 3), 4096);
                let mut code2 = [0i16; L_SUBFR];
                code2[0] = round(l_msu(l_deposit_h(code[0]), code[1], tmp));
                for i in 1..L_SUBFR - 1 {
                    let l = l_msu(l_deposit_h(code[i]), code[i + 1], tmp);
                    code2[i] = round(l_msu(l, code[i - 1], tmp));
                }
                code2[L_SUBFR - 1] = round(l_msu(
                    l_deposit_h(code[L_SUBFR - 1]),
                    code[L_SUBFR - 2],
                    tmp,
                ));
                gain_code = round(l_shl(l_gain_code, q_new));
                for i in 0..L_SUBFR {
                    let mut l = l_shl(l_mult(code2[i], gain_code), 5);
                    l = l_mac(l, exc2[i], gain_pit);
                    exc2[i] = round(l_shl(l, 1));
                }
                let corr_gain =
                    self.synthesis(p_aq, &mut exc2, q_new, &speech16k[i_subfr * 5 / 4..]);
                prms.parm(corr_gain, 4);
            }
        }
        self.old_speech.copy_from_slice(&old_speech[L_FRAME..]);
        self.old_wsp
            .copy_from_slice(&old_wsp[L_FRAME / OPL_DECIM..L_FRAME / OPL_DECIM + OLD_WSP]);
        self.old_exc
            .copy_from_slice(&old_exc[L_FRAME..L_FRAME + OLD_EXC]);
        (mode, prms.len())
    }

    /// `synthesis`: the 12.8 kHz synthesis of `exc` (scaled by `q_new`; modified), its tilt
    /// and HF noise, and the quantised high-band gain correction against the input
    /// `synth16k` (6–7 kHz band). Returns the gain index.
    fn synthesis(&mut self, aq: &[i16], exc: &mut [i16], q_new: i16, synth16k: &[i16]) -> i16 {
        let mut synth_hi = [0i16; M + L_SUBFR];
        let mut synth_lo = [0i16; M + L_SUBFR];
        synth_hi[..M].copy_from_slice(&self.mem_syn_hi);
        synth_lo[..M].copy_from_slice(&self.mem_syn_lo);
        syn_filt_32(aq, M, exc, q_new, &mut synth_hi, &mut synth_lo, L_SUBFR);
        self.mem_syn_hi.copy_from_slice(&synth_hi[L_SUBFR..]);
        self.mem_syn_lo.copy_from_slice(&synth_lo[L_SUBFR..]);
        let mut synth = [0i16; L_SUBFR];
        deemphasis_32(
            &synth_hi[M..],
            &synth_lo[M..],
            &mut synth,
            PREEMPH_FAC,
            &mut self.mem_deemph,
        );
        highpass_50hz_at_12k8(&mut synth, &mut self.mem_sig_out);
        let mut hf_sp: [i16; L_SUBFR16K] = synth16k[..L_SUBFR16K].try_into().unwrap();
        // HF noise at the excitation's energy
        let mut hf = [0i16; L_SUBFR16K];
        for v in &mut hf {
            *v = shr(random(&mut self.seed2), 3);
        }
        scale_sig(exc, -3);
        let q_new = sub(q_new, 3);
        let mut exp_ener = 0;
        let ener = extract_h(dot_product12(exc, exc, &mut exp_ener));
        exp_ener = sub(exp_ener, add(q_new, q_new));
        let mut exp = 0;
        let mut tmp = extract_h(dot_product12(&hf, &hf, &mut exp));
        if tmp > ener {
            tmp >>= 1;
            exp = add(exp, 1);
        }
        let mut l = l_deposit_h(div_s(tmp, ener));
        exp = sub(exp, exp_ener);
        isqrt_n(&mut l, &mut exp);
        let tmp = extract_h(l_shl(l, add(exp, 1)));
        for v in &mut hf {
            *v = mult(*v, tmp);
        }
        // tilt of the synthesis → the estimated HF gain
        hp400_12k8(&mut synth, &mut self.mem_hp400);
        let mut l = 1i32;
        for &s in &synth {
            l = l_mac(l, s, s);
        }
        let exp = norm_l(l);
        let ener = extract_h(l_shl(l, exp));
        let mut l = 1i32;
        for i in 1..L_SUBFR {
            l = l_mac(l, synth[i], synth[i - 1]);
        }
        let tmp = extract_h(l_shl(l, exp));
        let fac = if tmp > 0 { div_s(tmp, ener) } else { 0 };
        let gain1 = sub(32767, fac);
        let gain2 = shl(mult(sub(32767, fac), 20480), 1);
        let (weight1, weight2) = if self.vad_hist > 0 {
            (0, 32767)
        } else {
            (32767, 0)
        };
        let mut tmp = add(mult(weight1, gain1), mult(weight2, gain2));
        if tmp != 0 {
            tmp = add(tmp, 1);
        }
        let hp_est_gain = tmp.max(3277);
        // the noise through the weighted LPC, both band-passed; the gain the input needs
        let mut ap = [0i16; MP1];
        weight_a(aq, &mut ap, 19661);
        wb_syn_filt(&ap, &mut hf, &mut self.mem_syn_hf);
        filt_6k_7k(&mut hf, &mut self.mem_hf);
        filt_6k_7k(&mut hf_sp, &mut self.mem_hf2);
        scale_sig(&mut hf_sp, -1);
        let mut exp_ener = 0;
        let ener = extract_h(dot_product12(&hf_sp, &hf_sp, &mut exp_ener));
        let mut exp = 0;
        let mut tmp = extract_h(dot_product12(&hf, &hf, &mut exp));
        if tmp > ener {
            tmp >>= 1;
            exp = add(exp, 1);
        }
        let mut l = l_deposit_h(div_s(tmp, ener));
        exp = sub(exp, exp_ener);
        isqrt_n(&mut l, &mut exp);
        let hp_calc_gain = extract_h(l_shl(l, exp));
        let hang = self.dtx.dtx_hangover_count;
        let l = l_shl(l_mult(hang, 4681), 15);
        self.gain_alpha = mult(self.gain_alpha, extract_h(l));
        if hang > 6 {
            self.gain_alpha = 32767;
        }
        let hp_est_gain = hp_est_gain >> 1;
        let hp_corr_gain = add(
            mult(hp_calc_gain, self.gain_alpha),
            mult(sub(32767, self.gain_alpha), hp_est_gain),
        );
        let mut dist_min = 32767;
        let mut ind = 0;
        for (i, &g) in HP_GAIN.iter().enumerate() {
            let d = sub(hp_corr_gain, g);
            let dist = mult(d, d);
            if dist_min > dist {
                dist_min = dist;
                ind = i as i16;
            }
        }
        ind
    }
}
