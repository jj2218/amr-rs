# amr-wb

Pure-Rust, bit-exact AMR-WB (3GPP TS 26.173, ITU-T G.722.2) speech decoder and encoder, with
RFC 4867 packing (`.awb` files, RTP payloads octet-aligned and bandwidth-efficient). No
`unsafe`, no dependencies. The decoder is ported from opencore-amr 0.1.6, the encoder from
vo-amrwbenc 0.1.3 (both Apache-2.0), each corrected to the specification wherever it differs.

```rust
use amr_wb::{Decoder, L_FRAME16K, rfc4867};

let mut dec = Decoder::new();
let mut pcm = [0i16; L_FRAME16K]; // 20 ms at 16 kHz
for frame in rfc4867::read_storage(&std::fs::read("speech.awb")?)? {
    let f = frame?; // frames before a damaged one are still delivered
    dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
}

let mut enc = amr_wb::Encoder::new(true); // with DTX (comfort noise in pauses)
let frames: Vec<_> = samples.as_chunks::<L_FRAME16K>().0.iter()
    .map(|pcm| enc.encode(amr_wb::Mode::MD1265, pcm))
    .collect();
std::fs::write("out.awb", rfc4867::write_storage(&frames)?)?;
```

## Results
Measured, with the recipe that reproduces each number, in
[`docs/RESULTS.md`](https://github.com/jj2218/amr-rs/blob/main/docs/RESULTS.md#amr-wb-amr-wb-decoder)
(decoder) and [its encoder section](https://github.com/jj2218/amr-rs/blob/main/docs/RESULTS.md#amr-wb-encoder-amr-wb-ts-26190).

- **Bit-exact on the whole 3GPP TS 26.174 decoder set** (v19.0.0): 221/221 cases, 119,660
  frames, every mode, mode switching and DTX. opencore-amr gets 189/221.
- **Encoder bit-exact on the whole TS 26.174 encoder set**: 3,101/3,101 cases (every mode,
  DTX, mode switching, frame and bit synchronisation), with the two encoder changes since
  vo-amrwbenc (CR 26.173-0032 and -0036). vo-amrwbenc gets 3,006/3,101.
- **Encoder 1.43× the speed of vo-amrwbenc** (1.34–1.54× per mode).
- **Differentially tested against the TS 26.173 reference C**: the whole decoder on 1,600
  stress streams (corrupted and random bits, bad, lost and degraded frames, SID and NO_DATA
  sequences, mode switching), the whole encoder on 64 stress streams (speech up to clipping,
  noise, tones, silence, homing, DTX, mode switching), and the filters and operators of both
  one by one, bit-exact.
- **`.awb` files**: every vector, packed per RFC 4867 and read back, decodes bit-exactly;
  the reader agrees with the TS 26.173 reference's on 20,000 random frames, and the writer's
  output for the encoder's frames (DTX, every mode) is byte-equal to the reference's.
- **≈2× the speed of opencore-amr's C decoder** (-O3; 1.85× at 6.60 kbit/s, 1.98–2.09× at the
  other modes), with runtime AVX2 dispatch and a portable fallback.
- **Fuzzed against the TS 26.173 reference**: the decoder on arbitrary frames (1 h), the
  encoder on arbitrary PCM (30 min, comfort noise included).
- **Ten defects of opencore's decoder fixed** (saturation, rounding, shift counts, the
  synthesis filters' scaling for comfort noise, concealment gains, comfort-noise muting): see
  the results page. One more, an out-of-bounds table read shared with the reference C, would
  be a panic in Rust and is handled exactly.

**One deliberate difference from TS 26.173:** NO_DATA frames are handled as before
CR 26.173-0027..29 (2007), as opencore does. A NO_DATA frame during speech starts the
comfort-noise hangover at once; the current specification first waits for 8 frames without
voice activity. Both agree on every conformance vector and differ only in the comfort noise
after NO_DATA frames that arrive mid-speech.

License: Apache-2.0 (see LICENSE and NOTICE).
