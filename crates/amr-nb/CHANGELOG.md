# Changelog

One section per release, newest first. The release commit in the repository carries the
section of each crate it releases.

## 0.2.0 (2026-11-18)
One API shape for `amr-nb` and `amr-wb`. Breaking changes, with the migration:
- `Mode::MRDTX` is removed: it is not a bit rate. `Encoder::encode` no longer panics on it,
  and a frame can no longer carry it. Migration: pass a speech mode (DTX is set in
  `Encoder::new`); a SID frame's `mode` is its mode indication, as before.
- `Display` for `Mode` is removed (it printed the same as `Debug`). Migration: `{:?}`.

Added
- `Mode::ALL` (the 8 modes), `Mode::from_index`, `Mode::bits` (serial bits per speech frame).
- `L_FRAME` (160 samples per frame) and `SID_BITS` (35: a SID frame's comfort-noise bits, TS
  26.101 table 3).
- `Hash` for `Mode` and `RxFrameType`.

Fixed
- `rfc4867` writes a SID_FIRST frame's comfort-noise bits as zeros, as TS 26.101 4.2.3
  requires (it sent the encoder's analysis, as opencore-amr does). Decoders don't use them.

## 0.1.1 (2026-11-14)
Changed
- Documentation for users: the README's results summarise the measured conformance, speed
  and robustness, and the code comments no longer refer to development notes.
- `rust-version = "1.88"` (the oldest Rust that builds the library and the example).

## 0.1.0 (2026-11-14)
Added
- First release: AMR-NB (3GPP TS 26.073) encoder and decoder in pure Rust, with
  discontinuous transmission (VAD option 1 or 2) and RFC 4867 packing (`.amr` files, RTP
  payloads octet-aligned and bandwidth-efficient).
- Bit-exact on the whole 3GPP test-vector set, differentially tested against the 3GPP
  reference C, and faster than opencore-amr's C on encoding and decoding.
