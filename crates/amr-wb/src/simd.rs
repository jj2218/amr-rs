//! Runtime CPU dispatch: the only module allowed to use `unsafe`.
//!
//! The frame decoder has one source and two copies: the portable build and an AVX2 build of
//! the *same* source (`#[target_feature(enable = "avx2")]`), picked per frame from
//! `is_x86_feature_detected!` (std caches the answer in an atomic). Same source, so the same
//! wrapping integer results: the AVX2 copy only lets the compiler use 256-bit lanes. The
//! encoder's hot kernels (`kernels!`) are dispatched the same way per call: each has one
//! portable body (`#[inline(always)]`, `*_portable` in its module, with the helpers it runs)
//! and an AVX2 copy of it. A copy gets AVX2 only for what inlines into it.
//! `--cfg portable` compiles the dispatch out (`just test` runs both paths).
#![allow(unsafe_code)] // the approved scope: calls into #[target_feature] copies, below

use crate::cnst::{L_FRAME16K, L_SUBFR, M, RxFrameType};
use crate::decoder::State;
use crate::encoder::kernels as k;

/// `State::decode` on the best instruction set of this CPU.
pub(crate) fn decode(
    st: &mut State,
    mode: usize,
    prms: &[i16],
    frame_type: RxFrameType,
    out: &mut [i16; L_FRAME16K],
) {
    #[cfg(all(target_arch = "x86_64", not(portable)))]
    if std::arch::is_x86_feature_detected!("avx2") {
        #[target_feature(enable = "avx2")]
        fn avx2_copy(
            st: &mut State,
            mode: usize,
            prms: &[i16],
            frame_type: RxFrameType,
            out: &mut [i16; L_FRAME16K],
        ) {
            st.decode(mode, prms, frame_type, out)
        }
        // SAFETY: `avx2_copy` only requires AVX2, which was just detected.
        return unsafe { avx2_copy(st, mode, prms, frame_type, out) };
    }
    st.decode(mode, prms, frame_type, out)
}

/// `pub fn name(args) = portable_body;` → the dispatching `name` (amr-nb's macro).
macro_rules! kernels {
    ($($(#[$m:meta])* $vis:vis fn $name:ident($($arg:ident: $ty:ty),* $(,)?) $(-> $ret:ty)?
        = $body:path;)*) => {$(
        $(#[$m])*
        #[allow(clippy::too_many_arguments)] // the C signatures
        $vis fn $name($($arg: $ty),*) $(-> $ret)? {
            #[cfg(all(target_arch = "x86_64", not(portable)))]
            if std::arch::is_x86_feature_detected!("avx2") {
                #[target_feature(enable = "avx2")]
                fn avx2_copy($($arg: $ty),*) $(-> $ret)? {
                    $body($($arg),*)
                }
                // SAFETY: `avx2_copy` only requires AVX2, which was just detected.
                return unsafe { avx2_copy($($arg),*) };
            }
            $body($($arg),*)
        }
    )*};
}

kernels! {
    /// `ACELP_4t64_fx` (see `acelp4::acelp_4t64_fx_portable`).
    pub(crate) fn acelp_4t64_fx(dn: &mut [i16; L_SUBFR], cn: &[i16], hh: &[i16],
        code: &mut [i16; L_SUBFR], y: &mut [i16; L_SUBFR], nbbits: i16, ser_size: i16,
        index_out: &mut [i16; 8]) = k::acelp_4t64_fx_portable;
    /// `ACELP_2t64_fx` (see `acelp2::acelp_2t64_fx_portable`).
    pub(crate) fn acelp_2t64_fx(dn: &mut [i16; L_SUBFR], cn: &[i16], hh: &[i16],
        code: &mut [i16; L_SUBFR], y: &mut [i16; L_SUBFR], index: &mut i16)
        = k::acelp_2t64_fx_portable;
    /// `cor_h_x` (see `acelp2::cor_h_x_portable`).
    pub(crate) fn cor_h_x(h: &[i16], x: &[i16], dn: &mut [i16; L_SUBFR]) = k::cor_h_x_portable;
    /// `Convolve`, `Residu`, `Syn_filt` (see `filters`).
    pub(crate) fn convolve(x: &[i16], h: &[i16], y: &mut [i16], l: usize) = k::convolve_portable;
    pub(crate) fn residu(a: &[i16], x: &[i16], at: usize, y: &mut [i16], lg: usize)
        = k::residu_portable;
    pub(crate) fn syn_filt(a: &[i16], x: Option<&[i16]>, y: &mut [i16], mem: &mut [i16],
        update: bool) = k::syn_filt_portable;
    /// `Pitch_fr4` (see `pitch::pitch_fr4_portable`).
    pub(crate) fn pitch_fr4(exc: &[i16], at: usize, xn: &[i16], h: &[i16], t0_min: i16,
        t0_max: i16, pit_frac: &mut i16, i_subfr: usize, t0_fr2: i16, t0_fr1: i16,
        l_subfr: usize) -> i16 = k::pitch_fr4_portable;
    /// `Decim_12k8` (see `preprocess::decim_12k8_portable`).
    pub(crate) fn decim_12k8(sig16k: &[i16], sig12k8: &mut [i16],
        mem: &mut [i16; 2 * crate::encoder::kernels::NB_COEF_DOWN]) = k::decim_12k8_portable;
    /// `Pitch_med_ol` (see `pitch_ol::pitch_med_ol_portable`).
    pub(crate) fn pitch_med_ol(wsp: &[i16], at: usize, l_min: i16, l_max: i16, l_frame: usize,
        l_0: i16, gain: &mut i16, hp_wsp_mem: &mut [i16; 9], old_hp_wsp: &mut [i16],
        wght_flg: i16) -> i16 = k::pitch_med_ol_portable;
    /// `Autocorr`, `Az_isp`, `Qpisf_2s_46b`, `Qpisf_2s_36b` (see `lpc`).
    pub(crate) fn autocorr(x: &[i16], r_h: &mut [i16; M + 1], r_l: &mut [i16; M + 1])
        = k::autocorr_portable;
    pub(crate) fn az_isp(a: &[i16], isp: &mut [i16; M], old_isp: &[i16; M])
        = k::az_isp_portable;
    pub(crate) fn qpisf_2s_46b(isf1: &[i16; M], isf_q: &mut [i16; M], past_isfq: &mut [i16; M],
        indice: &mut [i16; 7], nb_surv: usize) = k::qpisf_2s_46b_portable;
    pub(crate) fn qpisf_2s_36b(isf1: &[i16; M], isf_q: &mut [i16; M], past_isfq: &mut [i16; M],
        indice: &mut [i16; 5], nb_surv: usize) = k::qpisf_2s_36b_portable;
}
