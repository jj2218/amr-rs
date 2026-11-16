//! Runtime CPU dispatch: the only module allowed to use `unsafe`.
//!
//! The frame decoder has one source and two copies: the portable build and an AVX2 build of
//! the *same* source (`#[target_feature(enable = "avx2")]`), picked per frame from
//! `is_x86_feature_detected!` (std caches the answer in an atomic). Same source, so the same
//! wrapping integer results: the AVX2 copy only lets the compiler use 256-bit lanes.
//! `--cfg portable` compiles the dispatch out (`just test` runs both paths).
#![allow(unsafe_code)] // the approved scope: calls into #[target_feature] copies, below

use crate::cnst::{L_FRAME16K, RxFrameType};
use crate::decoder::State;

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
