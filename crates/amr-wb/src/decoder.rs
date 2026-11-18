//! The decoder: state (`Decoder_State`), one frame (`pvDecoder_AmrWb`), homing
//! (`homing_amr_wb_dec.cpp`) and the frame loop of opencore's `amrwb/wrapper.cpp`.

use crate::basic_op::*;
use crate::bits::Bits;
use crate::cnst::*;
use crate::dtx::{DtxDecState, SPEECH};
use crate::excitation::*;
use crate::gains::{GainMem, PBUF2, dec_gain2_amr_wb, dec_gain2_amr_wb_init};
use crate::lpc::*;
use crate::synthesis::{SynthState, synthesis_amr_wb};
use crate::tables::*;

/// Past excitation kept between frames.
const OLD_EXC: usize = PIT_MAX + L_INTERPOL;

/// `Decoder_State`.
#[derive(Clone)]
pub(crate) struct State {
    old_exc: [i16; OLD_EXC],
    ispold: [i16; M],
    isfold: [i16; M],
    isf_buf: [i16; L_MEANBUF * M],
    past_isfq: [i16; M],
    tilt_code: i16,
    q_old: i16,
    qsubfr: [i16; 4],
    l_gc_thres: i32,
    syn: SynthState,
    seed: i16,
    old_t0: i16,
    lag_hist: [i16; 5],
    dec_gain: GainMem,
    seed3: i16,
    disp_mem: [i16; 8],
    prev_bfi: bool,
    state: i16,
    first_frame: bool,
    dtx: DtxDecState,
    vad_hist: i16,
}

impl State {
    fn new() -> Self {
        let mut st = Self {
            old_exc: [0; OLD_EXC],
            ispold: [0; M],
            isfold: [0; M],
            isf_buf: [0; L_MEANBUF * M],
            past_isfq: [0; M],
            tilt_code: 0,
            q_old: 0,
            qsubfr: [0; 4],
            l_gc_thres: 0,
            syn: SynthState::default(),
            seed: 0,
            old_t0: 0,
            lag_hist: [0; 5],
            dec_gain: [0; 23],
            seed3: 0,
            disp_mem: [0; 8],
            prev_bfi: false,
            state: 0,
            first_frame: true,
            dtx: DtxDecState::new(&ISF_INIT),
            vad_hist: 0,
        };
        st.reset(true);
        st
    }

    /// `pvDecoder_AmrWb_Reset`.
    fn reset(&mut self, reset_all: bool) {
        self.old_exc.fill(0);
        self.past_isfq.fill(0);
        self.old_t0 = 64;
        self.first_frame = true;
        self.l_gc_thres = 0;
        self.tilt_code = 0;
        self.disp_mem.fill(0);
        self.q_old = Q_MAX;
        self.qsubfr = [Q_MAX; 4];
        if reset_all {
            dec_gain2_amr_wb_init(&mut self.dec_gain);
            self.syn = SynthState {
                seed2: RANDOM_INITSEED,
                ..SynthState::default()
            };
            init_lagconc(&mut self.lag_hist);
            self.ispold = ISP_INIT;
            self.isfold = ISF_INIT;
            for b in self.isf_buf.as_chunks_mut::<M>().0 {
                b.copy_from_slice(&ISF_INIT);
            }
            self.seed = RANDOM_INITSEED;
            self.seed3 = RANDOM_INITSEED;
            self.state = 0;
            self.prev_bfi = false;
            self.dtx.reset(&ISF_INIT);
            self.vad_hist = 0;
        }
    }

    /// `pvDecoder_AmrWb`: one frame of `mode` (the previous mode for NO_DATA/lost frames).
    /// Inlined into `simd::decode`'s copies (one per instruction set).
    #[inline(always)]
    pub(crate) fn decode(
        &mut self,
        mode: usize,
        prms: &[i16],
        frame_type: RxFrameType,
        synth16k: &mut [i16; L_FRAME16K],
    ) {
        let mut prms = Bits::new(prms);
        let nb_bits = AMR_WB_COMPRESSED[mode];
        let mut exc2 = [0i16; L_FRAME];
        let mut isf = [0i16; M];
        let mut hf_isf = [0i16; M16K];
        let mut aq = [0i16; NB_SUBFR * (M + 1)];

        let new_dtx_state = self.dtx.rx_handler(frame_type);
        if new_dtx_state != SPEECH {
            self.dtx.dec(&mut exc2, new_dtx_state, &mut isf, &mut prms);
        }
        use RxFrameType::*;
        let (bfi, unusable_frame) = match frame_type {
            SpeechBad | SpeechProbablyDegraded => (true, false),
            NoData | SpeechLost => (true, true),
            _ => (false, false),
        };
        if bfi {
            self.state = (self.state + 1).min(6);
        } else {
            self.state >>= 1;
        }
        if self.dtx.dtx_global_state == crate::dtx::DTX {
            self.state = 5;
            self.prev_bfi = false;
        } else if self.dtx.dtx_global_state == crate::dtx::DTX_MUTE {
            self.state = 5;
            self.prev_bfi = true;
        }
        if new_dtx_state == SPEECH {
            let vad_flag = prms.bit();
            if !bfi {
                self.vad_hist = if vad_flag == 0 {
                    add_int16(self.vad_hist, 1)
                } else {
                    0
                };
            }
        }

        if new_dtx_state != SPEECH {
            let mut ispnew = [0i16; M];
            isf_isp(&isf, &mut ispnew);
            isp_az(&ispnew, &mut aq, true);
            let isf_tmp = self.isfold;
            for (j, i_subfr) in (0..L_FRAME).step_by(L_SUBFR).enumerate() {
                interp_hf_isf(&isf_tmp, &isf, j, &mut hf_isf);
                synthesis_amr_wb(
                    &aq[..M + 1],
                    &mut exc2[i_subfr..i_subfr + L_SUBFR],
                    0,
                    &mut synth16k[i_subfr * 5 / 4..],
                    1,
                    &mut hf_isf,
                    nb_bits,
                    new_dtx_state,
                    &mut self.syn,
                    self.vad_hist,
                    bfi,
                );
            }
            self.reset(false);
            self.isfold = isf;
            self.prev_bfi = bfi;
            self.dtx.dtx_global_state = new_dtx_state;
            return;
        }

        // past excitation, then this frame (+1 sample for the interpolation's look-ahead)
        let mut old_exc = [0i16; L_FRAME + 1 + OLD_EXC];
        old_exc[..OLD_EXC].copy_from_slice(&self.old_exc);
        let exc = OLD_EXC;

        // ISFs
        let mut ind = [0i16; 8];
        if nb_bits > NBBITS_7K {
            for (k, n) in [8, 8, 6, 7, 7, 5, 5].into_iter().enumerate() {
                ind[k] = prms.parm(n);
            }
            dpisf_2s_46b(
                &ind,
                &mut isf,
                &mut self.past_isfq,
                &self.isfold,
                &mut self.isf_buf,
                bfi,
            );
        } else {
            ind[0] = prms.parm(8);
            ind[1] = prms.parm(8);
            ind[2] = prms.parm(14);
            ind[3] = ind[2] & 0x007F;
            ind[2] >>= 7;
            ind[4] = prms.parm(6);
            dpisf_2s_36b(
                &ind,
                &mut isf,
                &mut self.past_isfq,
                &self.isfold,
                &mut self.isf_buf,
                bfi,
            );
        }
        let mut ispnew = [0i16; M];
        isf_isp(&isf, &mut ispnew);
        if self.first_frame {
            self.first_frame = false;
            self.ispold = ispnew;
        }
        interpolate_isp(&self.ispold, &ispnew, &INTERPOL_FRAC, &mut aq);
        self.ispold = ispnew;

        // stability factor from the ISF distance to the previous frame
        let mut l = 0i32;
        for i in 0..M - 1 {
            let t = sub_int16(isf[i], self.isfold[i]);
            l = mac_16by16_to_int32(l, t, t);
        }
        let t = extract_h(shl_int32(l, 8));
        let t = mult_int16(t, 26214);
        let stab_fac = shl_int16(20480i16.wrapping_sub(t), 1).max(0);
        let isf_tmp = self.isfold;
        self.isfold = isf;
        crate::stages::mark(crate::stages::D_PARAMS);

        let mut t0_min = 0i16;
        let mut q_new = 0i16;
        let mut code = [0i16; L_SUBFR];
        let mut excp = [0i16; L_SUBFR];
        let mut sub_exc2 = [0i16; L_SUBFR];
        for (sf, i_subfr) in (0..L_FRAME).step_by(L_SUBFR).enumerate() {
            // pitch lag
            let pit_flag = if i_subfr == 2 * L_SUBFR && nb_bits > NBBITS_7K {
                0
            } else {
                i_subfr
            };
            let (mut t0, mut t0_frac);
            if pit_flag == 0 {
                if nb_bits <= NBBITS_9K {
                    let index = prms.parm(8);
                    if index < (PIT_FR1_8B - PIT_MIN) * 2 {
                        t0 = PIT_MIN + (index >> 1);
                        t0_frac = sub_int16(index, shl_int16(sub_int16(t0, PIT_MIN), 1));
                        t0_frac = shl_int16(t0_frac, 1);
                    } else {
                        t0 = add_int16(index, PIT_FR1_8B - (PIT_FR1_8B - PIT_MIN) * 2);
                        t0_frac = 0;
                    }
                } else {
                    let mut index = prms.parm(9);
                    if index < (PIT_FR2 - PIT_MIN) * 4 {
                        t0 = PIT_MIN + (index >> 2);
                        t0_frac = sub_int16(index, shl_int16(sub_int16(t0, PIT_MIN), 2));
                    } else if index < ((PIT_FR2 - PIT_MIN) << 2) + ((PIT_FR1_9B - PIT_FR2) << 1) {
                        index -= (PIT_FR2 - PIT_MIN) << 2;
                        t0 = PIT_FR2 + (index >> 1);
                        t0_frac = sub_int16(index, shl_int16(sub_int16(t0, PIT_FR2), 1));
                        t0_frac = shl_int16(t0_frac, 1);
                    } else {
                        t0 = add_int16(
                            index,
                            PIT_FR1_9B - (PIT_FR2 - PIT_MIN) * 4 - (PIT_FR1_9B - PIT_FR2) * 2,
                        );
                        t0_frac = 0;
                    }
                }
                t0_min = (t0 - 8).max(PIT_MIN);
                if t0_min + 15 > PIT_MAX as i16 {
                    t0_min = PIT_MAX as i16 - 15;
                }
            } else if nb_bits <= NBBITS_9K {
                let index = prms.parm(5);
                t0 = t0_min + (index >> 1);
                t0_frac = sub_int16(index, shl_int16(t0 - t0_min, 1));
                t0_frac = shl_int16(t0_frac, 1);
            } else {
                let index = prms.parm(6);
                t0 = t0_min + (index >> 2);
                t0_frac = sub_int16(index, shl_int16(t0 - t0_min, 2));
            }
            if bfi {
                lagconceal(
                    &self.dec_gain[PBUF2..],
                    &self.lag_hist,
                    &mut t0,
                    self.old_t0,
                    &mut self.seed3,
                    unusable_frame,
                );
                t0_frac = 0;
            }
            pred_lt4(&mut old_exc, exc + i_subfr, t0, t0_frac, L_SUBFR + 1);

            // low-pass the adaptive codebook (always for 6.60 and 8.85 kbit/s)
            let select = if unusable_frame {
                1
            } else if nb_bits <= NBBITS_9K {
                0
            } else {
                prms.bit()
            };
            if select == 0 {
                // saturating as TS 26.173 (opencore sums plain products and wraps the << 1)
                let e = &old_exc[exc + i_subfr - 1..];
                for i in 0..L_SUBFR {
                    let l = mul_16by16_to_int32(5898, e[i]);
                    let l = mac_16by16_to_int32(l, 20972, e[i + 1]);
                    code[i] = amr_wb_round(mac_16by16_to_int32(l, 5898, e[i + 2]));
                }
                old_exc[exc + i_subfr..exc + i_subfr + L_SUBFR].copy_from_slice(&code);
            }

            // algebraic codebook
            if unusable_frame {
                for c in code.iter_mut() {
                    *c = noise_gen_amrwb(&mut self.seed) >> 3;
                }
            } else if nb_bits <= NBBITS_7K {
                dec_acelp_2p_in_64(prms.parm(12), &mut code);
            } else {
                let (widths, nbbits): (&[usize], i16) = if nb_bits <= NBBITS_9K {
                    (&[5, 5, 5, 5], 20)
                } else if nb_bits <= NBBITS_12K {
                    (&[9, 9, 9, 9], 36)
                } else if nb_bits <= NBBITS_14K {
                    (&[13, 13, 9, 9], 44)
                } else if nb_bits <= NBBITS_16K {
                    (&[13, 13, 13, 13], 52)
                } else if nb_bits <= NBBITS_18K {
                    (&[2, 2, 2, 2, 14, 14, 14, 14], 64)
                } else if nb_bits <= NBBITS_20K {
                    (&[10, 10, 2, 2, 10, 10, 14, 14], 72)
                } else {
                    (&[11, 11, 11, 11, 11, 11, 11, 11], 88)
                };
                for (k, &n) in widths.iter().enumerate() {
                    ind[k] = prms.parm(n);
                }
                dec_acelp_4p_in_64(&ind, nbbits, &mut code);
            }
            preemph_amrwb_dec(&mut code, self.tilt_code);
            let lag = if t0_frac > 2 { t0 + 1 } else { t0 };
            pit_shrp(&mut code, lag, PIT_SHARP);

            // gains
            let (index, nbits) = if nb_bits <= NBBITS_9K {
                (prms.parm(6), 6)
            } else {
                (prms.parm(7), 7)
            };
            let (gain_pit, mut l_gain_code) = dec_gain2_amr_wb(
                index,
                nbits,
                &code,
                bfi,
                self.prev_bfi,
                self.state,
                unusable_frame,
                self.vad_hist,
                &mut self.dec_gain,
            );

            // scale the excitation to the code gain's headroom
            let tmp = (*self.qsubfr.iter().min().unwrap()).min(Q_MAX);
            q_new = 0;
            let mut lt = l_gain_code;
            while lt < 0x0800_0000 && q_new < tmp {
                lt <<= 1;
                q_new += 1;
            }
            let gain_code = amr_wb_round(lt);
            scale_signal(
                &mut old_exc[i_subfr..i_subfr + OLD_EXC + L_SUBFR],
                q_new - self.q_old,
            );
            self.q_old = q_new;

            if !bfi {
                self.lag_hist.copy_within(0..4, 1);
                self.lag_hist[0] = t0;
                self.old_t0 = t0;
            }

            // adaptive part at Q_exc = -3 for the voicing and the 8.85/6.60 pitch enhancer
            for i in 0..L_SUBFR {
                let e = old_exc[exc + i_subfr + i];
                sub_exc2[i] = shr3_round(e);
            }
            let pit_sharp = if nb_bits <= NBBITS_9K {
                let ps = shl_int16(gain_pit, 1);
                if ps > 16384 {
                    for i in 0..L_SUBFR {
                        let t = mult_int16(sub_exc2[i], ps);
                        excp[i] = amr_wb_round(mul_16by16_to_int32(t, gain_pit) >> 1);
                    }
                }
                ps
            } else {
                0
            };
            let voice_fac = voice_factor(&sub_exc2, -3, gain_pit, &code, gain_code);
            self.tilt_code = (voice_fac >> 2) + 8192;

            // total excitation (kept for the next subframes' pitch)
            sub_exc2.copy_from_slice(&old_exc[exc + i_subfr..exc + i_subfr + L_SUBFR]);
            let mut max = 1i32;
            for i in 0..L_SUBFR {
                let mut l = shl_int32(mul_16by16_to_int32(code[i], gain_code), 5);
                l = mac_16by16_to_int32(l, old_exc[exc + i_subfr + i], gain_pit);
                let t = amr_wb_round(shl_int32(l, 1));
                old_exc[exc + i_subfr + i] = t;
                let t = t.wrapping_sub((t < 0) as i16);
                max |= (t ^ (t >> 15)) as i32;
            }
            let t = add_int16(norm_s(max as i16), q_new).wrapping_sub(1);
            self.qsubfr.copy_within(0..3, 1);
            self.qsubfr[0] = t;

            // phase dispersion and noise enhancer, for the synthesis only
            let disp_mode = if nb_bits <= NBBITS_7K {
                0
            } else if nb_bits <= NBBITS_9K {
                1
            } else {
                2
            };
            phase_dispersion(
                (l_gain_code >> 16) as i16,
                gain_pit,
                &mut code,
                disp_mode,
                &mut self.disp_mem,
            );
            let fac = mult_int16(stab_fac, 16384 - (voice_fac >> 1));
            // saturating as TS 26.173 (`Mpy_32_16`, `L_add`; opencore adds plainly)
            let mut lt = l_gain_code;
            if lt < self.l_gc_thres {
                lt = add_int32(lt, mpy_32_16(l_gain_code, 6226));
                if lt > self.l_gc_thres {
                    lt = self.l_gc_thres;
                }
            } else {
                lt = mpy_32_16(l_gain_code, 27536);
                if lt < self.l_gc_thres {
                    lt = self.l_gc_thres;
                }
            }
            self.l_gc_thres = lt;
            l_gain_code = mpy_32_16(l_gain_code, 32767 - fac);
            l_gain_code = add_int32(l_gain_code, mpy_32_16(lt, fac));

            // pitch-sharpened code (tilt from the voicing) plus adaptive part
            let tmp = (voice_fac >> 3) + 4096;
            let gain_code = amr_wb_round(shl_int32(l_gain_code, q_new));
            for i in 0..L_SUBFR {
                // one L_msu per neighbour, as TS 26.173 (opencore summed the neighbours into
                // a 16-bit argument first, which wraps on loud codevectors)
                let mut l = (code[i] as i32) << 16;
                if i + 1 < L_SUBFR {
                    l = msu_16by16_from_int32(l, code[i + 1], tmp);
                }
                if i > 0 {
                    l = msu_16by16_from_int32(l, code[i - 1], tmp);
                }
                l = mul_16by16_to_int32(amr_wb_round(l), gain_code);
                l = shl_int32(l, 5);
                l = mac_16by16_to_int32(l, sub_exc2[i], gain_pit);
                sub_exc2[i] = amr_wb_round(shl_int32(l, 1));
            }
            if nb_bits <= NBBITS_9K && pit_sharp > 16384 {
                for i in 0..L_SUBFR {
                    excp[i] = add_int16(excp[i], sub_exc2[i]);
                }
                agc2_amr_wb(&sub_exc2, &mut excp);
                sub_exc2 = excp;
            }

            if nb_bits <= NBBITS_7K {
                interp_hf_isf(&isf_tmp, &isf, sf, &mut hf_isf);
            } else {
                self.syn.mem_syn_hf[..M16K - M].fill(0);
            }
            let corr_gain = if nb_bits >= NBBITS_24K {
                prms.parm(4)
            } else {
                0
            };
            crate::stages::mark(crate::stages::D_EXC);
            synthesis_amr_wb(
                &aq[sf * (M + 1)..(sf + 1) * (M + 1)],
                &mut sub_exc2,
                q_new,
                &mut synth16k[i_subfr + (i_subfr >> 2)..],
                corr_gain,
                &mut hf_isf,
                nb_bits,
                new_dtx_state,
                &mut self.syn,
                self.vad_hist,
                bfi,
            );
        }

        self.old_exc
            .copy_from_slice(&old_exc[L_FRAME..L_FRAME + OLD_EXC]);
        scale_signal(&mut old_exc[exc..exc + L_FRAME], -q_new);
        self.dtx.activity_update(&isf, &old_exc[exc..exc + L_FRAME]);
        self.dtx.dtx_global_state = new_dtx_state;
        self.prev_bfi = bfi;
    }
}

/// The ISFs of subframe `j`: `isf_old` and `isf` mixed by `INTERPOL_FRAC[j]` (the 16 low ISFs
/// of the 6.60 kbit/s and DTX high band).
fn interp_hf_isf(isf_old: &[i16; M], isf: &[i16; M], j: usize, hf_isf: &mut [i16; M16K]) {
    for i in 0..M {
        let l = mul_16by16_to_int32(isf_old[i], sub_int16(32767, INTERPOL_FRAC[j]));
        hf_isf[i] = amr_wb_round(mac_16by16_to_int32(l, isf[i], INTERPOL_FRAC[j]));
    }
}

const MODE_24K: usize = 8;

/// `dhf_test`: whether the first `nparms` bits of `bits` equal the decoder homing frame.
fn dhf_test(bits: &[i16], mode: usize, nparms: i16) -> bool {
    const DHF: [&[i16]; 9] = [
        &DFH_M7K, &DFH_M9K, &DFH_M12K, &DFH_M14K, &DFH_M16K, &DFH_M18K, &DFH_M20K, &DFH_M23K,
        &DFH_M24K,
    ];
    let dhf = DHF[mode];
    let mut prms = Bits::new(bits);
    let mut param = [0i16; 32];
    let mut i = 0usize;
    let shift;
    if mode != MODE_24K {
        let mut j = 0i16;
        while nparms - 15 > j {
            param[i] = prms.parm(15);
            j += 15;
            i += 1;
        }
        let tmp = nparms - j;
        param[i] = prms.parm(tmp as usize);
        shift = 15 - tmp;
        param[i] = shl_int16(param[i], shift);
    } else {
        for (k, p) in param.iter_mut().enumerate().take(31) {
            *p = prms.parm(15);
            match k {
                10 => *p &= 0x61FF,
                17 => *p &= 0xE0FFu16 as i16, // the value has 15 bits: the same as the C's int mask
                24 => *p &= 0x7F0F,
                _ => {}
            }
        }
        param[31] = shl_int16(prms.parm(8), 7);
        i = 31;
        shift = 0;
    }
    let n = i;
    let mut j = 0i16;
    let mut i = 0;
    while i < n {
        j = param[i] ^ dhf[i];
        if j != 0 {
            break;
        }
        i += 1;
    }
    let mask = shl_int16(0x7fff >> shift, shift);
    j |= param[i] ^ (dhf[i] & mask);
    j == 0
}

/// AMR-WB decoder: one 20 ms frame of serial bits to 320 samples at 16 kHz.
///
/// The frame loop of the TS 26.173 reference decoder: NO_DATA and lost frames reuse the last
/// mode, decoder homing frames reset the decoder (and output the encoder homing pattern), and
/// the output is 14-bit (the 2 LSBs cleared).
#[derive(Clone)]
pub struct Decoder {
    st: Box<State>,
    mode_old: Mode,
    reset_flag: bool,
    reset_flag_old: bool,
}

impl Default for Decoder {
    fn default() -> Self {
        Self::new()
    }
}

impl Decoder {
    /// A decoder in its initial (homed) state.
    pub fn new() -> Self {
        Self {
            st: Box::new(State::new()),
            mode_old: Mode::MD66,
            reset_flag: false,
            reset_flag_old: true,
        }
    }

    /// Back to the initial state.
    pub fn reset(&mut self) {
        *self = Self::new();
    }

    /// Decodes one frame. `bits` holds the frame's serial bits, one per element (0 or 1), in
    /// TS 26.173 parameter order: `mode.bits()` of them for speech frames, [`SID_BITS`] for
    /// SID frames, none read for [`RxFrameType::NoData`] and [`RxFrameType::SpeechLost`].
    /// `mode` is ignored for SID, NO_DATA and lost frames: the last speech frame's is used (as
    /// the TS 26.173 reference does for `.awb` files; a SID's mode indication is not used).
    ///
    /// [`SID_BITS`]: crate::SID_BITS
    pub fn decode(
        &mut self,
        frame_type: RxFrameType,
        mode: Mode,
        bits: &[i16; MAX_SERIAL_SIZE],
        out: &mut [i16; L_FRAME16K],
    ) {
        crate::stages::start();
        use RxFrameType::*;
        let mode = if matches!(frame_type, NoData | SpeechLost) {
            self.reset_flag = false;
            self.mode_old
        } else {
            // SID frames: the previous mode, as the reference reads `.awb` and ITU files
            // (`Read_serial`: `*mode = st->prev_mode`); their mode indication is ignored.
            let mode = if matches!(frame_type, SidFirst | SidUpdate | SidBad) {
                self.mode_old
            } else {
                mode
            };
            self.mode_old = mode;
            if self.reset_flag_old {
                self.reset_flag = dhf_test(bits, mode as usize, PRMNOFSF[mode as usize]);
            }
            mode
        };
        if self.reset_flag && self.reset_flag_old {
            out.fill(EHF_MASK);
        } else {
            crate::simd::decode(&mut self.st, mode as usize, bits, frame_type, out);
        }
        for s in out.iter_mut() {
            *s &= !3;
        }
        if !self.reset_flag_old {
            self.reset_flag = dhf_test(bits, mode as usize, AMR_WB_COMPRESSED[mode as usize]);
        }
        if self.reset_flag {
            self.st.reset(true);
        }
        self.reset_flag_old = self.reset_flag;
        crate::stages::mark(crate::stages::D_TAIL);
    }
}
