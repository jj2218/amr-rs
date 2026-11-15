# amr-wb

Pure-Rust, bit-exact AMR-WB (3GPP TS 26.173, ITU-T G.722.2) speech decoder, with RFC 4867
packing (`.awb` files, RTP payloads octet-aligned and bandwidth-efficient). No `unsafe`, no
dependencies. Ported from opencore-amr 0.1.6 (Apache-2.0) and corrected to the specification
wherever opencore differs from it.

```rust
use amr_wb::{Decoder, L_FRAME16K, rfc4867};

let mut dec = Decoder::new();
let mut pcm = [0i16; L_FRAME16K]; // 20 ms at 16 kHz
for frame in rfc4867::read_storage(&std::fs::read("speech.awb")?)? {
    let f = frame?; // frames before a damaged one are still delivered
    dec.decode(f.frame_type, f.mode, &f.bits, &mut pcm);
}
```

## Results
Measured, with the recipe that reproduces each number, in
[`docs/RESULTS.md`](https://github.com/jj2218/amr-rs/blob/main/docs/RESULTS.md#amr-wb-amr-wb-decoder).

- **Bit-exact on the whole 3GPP TS 26.174 decoder set** (v19.0.0): 221/221 cases, 119,660
  frames, every mode, mode switching and DTX. opencore-amr gets 189/221.
- **Differentially tested against the TS 26.173 reference C**: the whole decoder on 1,600
  stress streams (corrupted and random bits, bad, lost and degraded frames, SID and NO_DATA
  sequences, mode switching), and its filters and operators one by one, bit-exact.
- **`.awb` files**: every vector, packed per RFC 4867 and read back, decodes bit-exactly;
  the reader agrees with the TS 26.173 reference's on 20,000 random frames.
- **≈1.25× the speed of opencore-amr's C decoder** (-O3), scalar code.
- **Fuzzed against the TS 26.173 reference** (1 h, the decoder on arbitrary frames).
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
