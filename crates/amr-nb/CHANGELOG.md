# Changelog

One section per release, newest first. The release commit in the repository carries the
section of each crate it releases.

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
