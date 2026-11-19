# Changelog

One section per release, newest first. The release commit in the repository carries the
section of each crate it releases.

## 0.4.2 (2026-11-19)
Added
- `Debug` for `Decoder` and `Encoder` (the encoder shows its `Dtx` setting; the state stays
  opaque), so a struct holding a codec can derive `Debug`.
- An example CLI, `examples/awb.rs`: `encode` raw 16 kHz PCM → `.awb` (any mode, optional
  DTX) and `decode` `.awb` → PCM. Every 3GPP encoder input through it equals its decoder vector.

Changed
- Faster decoder on x86-64 CPUs with AVX2: about 2.5× the speed of opencore-amr (was 2.3×;
  2.28–2.65× per mode). The output is unchanged (bit-exact on the 3GPP vectors and against
  the 3GPP reference).

## 0.4.1 (2026-11-18)
Changed
- Faster decoder on x86-64 CPUs with AVX2: about 2.3× the speed of opencore-amr (was 2.2×).
  The output is unchanged (bit-exact on the 3GPP vectors and against the 3GPP reference).
- The package description names DTX and RFC 4867 packing; the keywords include `amr`.
- README: a Scope section (single channel, the supported RTP options, Rust 1.88) and a License
  section with links that work on crates.io.

## 0.4.0 (2026-11-18)
One API shape for `amr-nb` and `amr-wb`. Breaking changes, with the migration:
- `Encoder::new` takes a `Dtx` (`Off`, `On`) instead of a `bool`. Migration:
  `Encoder::new(true)` → `Encoder::new(Dtx::On)`, `false` → `Dtx::Off`.
- `L_FRAME16K` is renamed `L_FRAME` (320 samples per frame). Migration: rename.

Added
- `Default` for `Encoder` (without DTX).

## 0.3.1 (2026-11-18)
Changed
- Faster on x86-64 CPUs with AVX2 (detected at run time): the decoder is about 2.2× the speed of
  opencore-amr (was 2.0×), the encoder about 1.6× vo-amrwbenc (was 1.43×; 1.51–1.61× per mode).
  The output is unchanged (bit-exact on the 3GPP vectors and against the 3GPP reference); the
  portable code's speed is unchanged.

## 0.3.0 (2026-11-17)
Changed
- The encoder is about 1.4× the speed of vo-amrwbenc (was 0.46×): 1.43× on average, 1.34–1.54×
  per mode. The output is unchanged (bit-exact on the whole 3GPP TS 26.174 encoder set and
  against the 3GPP reference). On x86-64 CPUs with AVX2 its hot kernels use it, detected at
  run time; other CPUs run the portable code.

## 0.2.0 (2026-11-16)
Added
- `Encoder`: AMR-WB (3GPP TS 26.173, ITU-T G.722.2) encoder in pure Rust, every mode, with
  voice activity detection and discontinuous transmission (comfort noise in pauses: SID and
  NO_DATA frames) and encoder homing. Ported from vo-amrwbenc 0.1.3 and corrected to the
  specification where it differs, including the two later encoder changes (CR 26.173-0032
  and -0036). Its frames go straight into `rfc4867::write_storage` / `write_rtp`.
- Bit-exact on the whole 3GPP TS 26.174 encoder set: 3,101/3,101 cases (every mode, DTX,
  mode switching, frame and bit synchronisation); vo-amrwbenc gets 3,006/3,101.
  Differentially tested against the 3GPP reference C (stress streams, every module) and
  fuzzed against it.
- Speed: the encoder is not yet tuned, about half vo-amrwbenc's (0.46× on average, 0.42–0.53×
  per mode).

Changed
- About 2× the speed of opencore-amr's decoder (was 1.25×): 1.85× at 6.60 kbit/s, 1.98–2.09×
  at the other modes. The output is unchanged (bit-exact on every 3GPP test vector and
  against the 3GPP reference). On x86-64 CPUs with AVX2 the decoder uses it, detected at
  run time; other CPUs run the portable code.
- `rfc4867`: a SID_FIRST frame's 35 comfort-noise bits are written as zeros, as the 3GPP
  reference's `.awb` writer does (decoders don't use them), so an `.awb` file from the
  encoder is byte-identical to the reference's.

## 0.1.1 (2026-11-15)
Fixed
- The decoder could panic (index out of bounds) on rare comfort-noise frames: a crafted
  stream could crash it. Now handled exactly; the output is unchanged (still bit-exact).
- SID (comfort-noise) frames are decoded with the previous frame's mode, as the 3GPP
  reference does when it reads `.awb` files; their mode indication is no longer used, so a
  corrupted one has no effect. The 3GPP test vectors decode unchanged.

Added
- `rfc4867`: RFC 4867 packing, `.awb` files and RTP payloads (octet-aligned and
  bandwidth-efficient), checked against the 3GPP reference's own file reader.

## 0.1.0 (2026-11-15)
Added
- First release: AMR-WB (3GPP TS 26.173, ITU-T G.722.2) decoder in pure Rust, ported from
  opencore-amr and corrected to the specification where opencore differs from it.
- Bit-exact on the whole 3GPP TS 26.174 decoder test-vector set (221 cases) and
  differentially tested against the 3GPP reference C; about 1.25× the speed of opencore-amr.
