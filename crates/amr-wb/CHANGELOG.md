# Changelog

One section per release, newest first. The release commit in the repository carries the
section of each crate it releases.

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
