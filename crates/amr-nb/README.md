# amr-nb

A pure-Rust AMR-NB speech codec (3GPP TS 26.073): encoder, decoder, discontinuous transmission
with VAD option 1 or 2, and RFC 4867 packing (`.amr` files, RTP payloads octet-aligned and
bandwidth-efficient). Ported from opencore-amr 0.1.6.

```rust
use amr_nb::{Decoder, Dtx, Encoder, Mode, rfc4867};

// Encode: 160 samples of 8 kHz, 16-bit PCM per 20 ms frame.
let mut enc = Encoder::new(Dtx::Off);
let frames: Vec<rfc4867::Frame> = pcm
    .chunks_exact(160)
    .map(|chunk| enc.encode(Mode::MR122, chunk.try_into().unwrap()))
    .collect();
std::fs::write("speech.amr", rfc4867::write_storage(&frames)?)?;

// Decode.
let mut dec = Decoder::new();
let mut out = [0i16; 160];
for frame in rfc4867::read_storage(&std::fs::read("speech.amr")?)? {
    let f = frame?; // examples/amr.rs conceals other codecs' frames instead
    dec.decode(f.frame_type, f.mode, &f.bits, &mut out);
}
```

A command-line example converts raw PCM to `.amr` and back:
`cargo run --release --example amr -- encode speech.raw speech.amr 12.2`.

## Results
Measured, with the recipe that reproduces each number, in
[`docs/RESULTS.md`](https://github.com/jj2218/amr-rs/blob/main/docs/RESULTS.md).

- **Bit-exact on the whole 3GPP TS 26.074 test-vector set**: linear PCM, A-law and µ-law,
  plus the bit and frame synchronisation sequences. Decoder 191/191 cases per law (573 cases,
  252,063 frames), also through `.amr` storage and both RTP packings; encoder 199/199 per law
  (no DTX, mode switching, DTX with VAD option 1 and 2); frame sync 3,840/3,840, bit sync 24/24.
- **Differentially tested against the TS 26.073 reference C**: its basic operators, math and
  encoder and decoder functions one by one, the whole encoder (80 synthetic signals × 500
  frames) and the whole decoder (20 × 2,000 random frames of every receive type), bit-exact.
  Fuzzed against it as well: the decoder target found one defect the port had inherited from
  opencore (fixed); after the fix, one-hour runs of the decoder and encoder targets (1.18 M and
  1.02 M runs) found no difference.
- **Faster than opencore-amr 0.1.6** compiled from source at `-O3`: decoder ≈3.06×
  (4.5 vs 13.6 µs per 20 ms frame), encoder ≈2.25× (24.2 vs 54.6 µs; per mode 2.08–2.39×).
  One pinned P-core of an i7-12700H (±5%), 2026-11-14; x86-64 baseline build with runtime AVX2
  dispatch. Without the AVX2 paths (`--cfg amr_nb_portable`): decoder 3.02×, encoder 1.78×.
- **Twelve defects found in the legacy C** (opencore-amr), each absent in this port, found by
  the vectors, differential tests, valgrind, inspection and fuzzing: no encoder/decoder homing;
  a `const` input buffer overwritten; an uninitialised read and two undefined-behaviour cases
  (signed overflow, an oversized shift); no PCM written for frame types 9–14; no input length
  check; the quality bit ignored; no VAD option 2; wrapping arithmetic where TS 26.073
  saturates; an unreachable overflow path; a comfort-noise energy clamp.
- **Safe Rust**: the only `unsafe` is `src/simd.rs` (runtime AVX2 dispatch, SSE2 kernels),
  plus FFI in the tests.

## Scope
Single channel. RTP without interleaving, CRCs or robust sorting (the SDP defaults); `.amr`
files without the multi-channel header. Frame types 9–14 (other codecs' SID, reserved) are
reported as errors.

## License
Apache-2.0 ([LICENSE](https://github.com/jj2218/amr-rs/blob/main/crates/amr-nb/LICENSE)).
Derived from opencore-amr 0.1.6 (Apache-2.0) and, through it, from the 3GPP TS 26.073 code;
see [NOTICE](https://github.com/jj2218/amr-rs/blob/main/crates/amr-nb/NOTICE).
