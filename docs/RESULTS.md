# Results

amr-nb's measured results. Every number is reproducible with the recipe or command in its
heading; the oracles (3GPP vectors and reference C, opencore-amr) are downloaded and verified on
first use (see the README).

## Conformance (`just conform`)
The whole 3GPP TS 26.074 set: NODTX (all 8 modes, T00–T20, T22), T21 mode switching, DTX
(VAD1) and DTX2 (VAD2), each in linear PCM, A-law and µ-law (clause 9), plus the bit and frame
synchronisation sequences (clause 8). Bit-exact = every output sample (decoder) or every `.COD`
word (encoder: frame type, 244 bits, mode) identical. A-law/µ-law: the harness expands the
8-bit input with G.711 before encoding and compresses the decoder's output before comparing
(`tests/common/g711.rs`, a test-side helper).

| decoder | cases | frames | bit-exact (linear / A-law / µ-law) |
|---|---:|---:|---|
| amr-nb (Rust) | 573 | 252,063 | **191/191 / 191/191 / 191/191** |
| amr-nb via `.amr` storage (RFC 4867 §5) | 573 | 252,063 | **191/191 / 191/191 / 191/191** |
| amr-nb via RTP octet-aligned (3 frames/packet) | 573 | 252,063 | **191/191 / 191/191 / 191/191** |
| amr-nb via RTP bandwidth-efficient (3 frames/packet) | 573 | 252,063 | **191/191 / 191/191 / 191/191** |
| opencore 0.1.6 (C), raw, linear | 191 | 84,021 | 1/191 (no homing) |
| opencore 0.1.6 (C) + our homing shim | 573 | 252,063 | 191/191 / 190/191 / 191/191 (defect 10) |

| encoder | cases | frames | bit-exact (linear / A-law / µ-law) |
|---|---:|---:|---|
| amr-nb (Rust), all | 597 | 252,087 | **199/199 / 199/199 / 199/199** |
| … NODTX + T21 mode switching | 531 | 180,843 | 177/177 / 177/177 / 177/177 |
| … DTX, VAD option 1 | 33 | 47,268 | 11/11 / 11/11 / 11/11 |
| … DTX2, VAD option 2 | 33 | 23,976 | 11/11 / 11/11 / 11/11 |
| opencore 0.1.6 (C) + our homing shim | 564 | 228,111 | 188/188 / 185/188 / 182/188 (no VAD2 build; defect 10) |
| frame sync, amr-nb (§8.2, 8 modes × 160 shifts × 3 laws) | 3,840 | 3,840 | **3,840/3,840** |
| frame sync, opencore | 3,840 | 3,840 | 2,242/3,840 (defect 10) |
| bit sync, amr-nb and opencore (§8.1, 13 triplets × 8 modes × 3 laws) | 24 | 936 | **24/24** |

The 3GPP reference C (TS 26.073, `.cache/3gpp/ref26073`) reproduces every vector, the SYNC set
included (checked with its own `encoder`/`decoder` on the cases opencore misses and on all
3,840 sync cases). opencore's misses are its arithmetic (defect 10). The conformance run pins
them for the C oracle (exact first mismatch; per-mode sync counts), so a change either way fails
`conform`. Bit sync: of the 13 bit-shifted homing triplets, exactly the unshifted one yields a
decoder homing frame, in every mode and law.

The "via" rows pack every frame (speech, SID_FIRST, SID_UPDATE, NO_DATA) with
`amr_nb::rfc4867`, read it back and decode. The encoder vectors can't show encoder homing (raw
opencore also gets 188/188: their homing frames come first), so a unit test checks it:
after a homing frame the encoder continues exactly like a fresh one.

## Against the TS 26.073 reference (`cargo test -p amr-nb`)
The unit-test oracle is the 3GPP reference C (`crates/ref26073-sys`: built from the pinned
download, every symbol prefixed `ref_`). All green, values and state (the reference's overflow
flag is one sticky global, compared only where the port tracks the same thing):
- **Basic operators** (`basicop2.c`, `oper_32b.c`) and math (`Inv_sqrt`, `Log2`, `Pow2`,
  `sqrt_l_exp`) on edge-biased random input.
- **Encoder**: LP analysis (`Autocorr`, `Lag_window`, `Levinson`), `Az_lsp`, `Q_plsf_3`,
  `Q_plsf_5`, `Pre_Process`, `Pitch_ol`, VAD1, VAD2 (`vad2`, `r_fft`, `LTP_flag_update`),
  `Pitch_fr`, `G_pitch`, `cbsearch` (all modes), `search_10and8i40`, `cor_h`, `cor_h_x`,
  `set_sign12k2`, `Convolve`, `G_code`, `q_gain_code`, `Qua_gain`, `MR475_gain_quant`,
  `MR795_gain_quant`, `calc_filt_energies`, `gc_pred`.
- **Decoder**: `D_plsf_3`, `D_plsf_5`, LSP/LSF conversions, `Int_lpc`, `Int_lsf`, `lsp_avg`,
  `Dec_lag3`, `Dec_lag6`, every codebook decoder, `d_gain_code`, `Dec_gain`, `d_gain_pitch`,
  `ec_gain_*`, `Cb_gain_average`, `Ex_ctrl`, `Bgn_scd`, `ph_disp`, `Pred_lt_3or6`, `Syn_filt`,
  `Residu`, `Weight_Ai`, `preemphasis`, `agc`, `agc2`, `Post_Filter`, `Post_Process`, the DTX
  decoder (`rx_dtx_handler`, `dtx_dec_activity_update`, comfort-noise builders, `A_Refl`).
- **Whole encoder**: 80 synthetic signals × 500 frames (tones, noise, silence, clipping; mode
  switches; 40 with VAD1 DTX) vs `Speech_Encode_Frame` + `sid_sync`, bit-exact.
- **Whole decoder**: 20 × 2,000 random frames (every RX type incl. bad, degraded, SID bad and
  NO_DATA; random modes and bits) vs `Speech_Decode_Frame`, bit-exact.

opencore stays the oracle only where the reference has nothing to say: inputs outside a
contract (`div_s`, and `Cb_gain_average` with a non-positive average LSP, where the reference
aborts), functions that are `static` there (`Chebps` on its non-saturating inputs, `Vq_subvec*`,
the `agc` energies; each covered through its public caller against the reference), the RFC 4867
reader (compared via opencore's own unpacking) and the raw `comp_corr` kernel. And it is the
speed baseline (`just bench`, `BENCH=kern`).

## Packed formats, RFC 4867 (`cargo test --release -p amr-nb rfc4867`)
`.amr` storage and RTP payloads (octet-aligned, bandwidth-efficient; single channel, no
interleaving/CRC/robust sorting). Oracles, all green:
- **Writer vs opencore's encoder**: 4,000 frames (DTX on, all 8 modes, SID_FIRST/SID_UPDATE/
  NO_DATA) packed by us are byte-identical to `Encoder_Interface_Encode`'s output.
- **Reader vs opencore's `Decoder_Interface_Decode`**: identical PCM on that encoder stream and
  on 20,000 random frames (all speech modes, SID with random STI/mode indication, NO_DATA,
  random padding bits).
- **Round trip** (all 8 modes × good/bad, SID first/update/bad × 8 mode indications, NO_DATA;
  1–4 frames per packet), frame sizes vs opencore's tables, header layouts worked out by hand
  from the RFC, and every truncation/reserved-FT error path. `.amr` files are read as a stream
  (frames before a bad one are kept; types 9–11 are skipped), so memory stays O(1) per frame.

opencore has no bandwidth-efficient mode and no `.amr` header handling; those are checked by
round trip, hand-derived layouts and the conformance rows above.

## Decoder speed (`just bench`)
µs per 20 ms frame, best of 15, one pinned P-core of an i7-12700H (laptop: ±5%, and Rust and
C drift together, so compare A/B runs back to back). Rust: release profile (fat LTO), portable
x86-64 (SSE2 baseline, no `target-cpu`; the encoder picks AVX2 copies at runtime, `src/simd.rs`). C: opencore compiled from source at `-O3` (Debian's
distro build at `-O2` is slower still, ~15 µs).

| mode | Rust | C | speedup |
|---|---:|---:|---:|
| 4.75 | 4.6 | 14.0 | 3.03× |
| 5.15 | 4.5 | 13.7 | 3.06× |
| 5.90 | 4.5 | 13.8 | 3.05× |
| 6.70 | 4.5 | 13.9 | 3.08× |
| 7.40 | 4.3 | 13.2 | 3.09× |
| 7.95 | 4.5 | 14.0 | 3.10× |
| 10.2 | 4.3 | 12.8 | 3.02× |
| 12.2 | 4.5 | 13.5 | 3.02× |
| **all** | **4.5** | **13.6** | **≈3.06×** |

Measured 2026-11-14. Portable build (`--cfg amr_nb_portable`, SSE2 only): 3.02×.

## Encoder speed (`just bench`)
Same machine and builds as the decoder table; input: the NODTX speech vectors (`.INP`), DTX
off, 59,856 frames per rep (2026-11-12). The C drifts ±5% with temperature (54.6–63 µs over
the day), so compare the ratio, or A/B runs interleaved back to back.

| mode | Rust | C | speedup |
|---|---:|---:|---:|
| 4.75 | 22.2 | 50.6 | 2.28× |
| 5.15 | 18.4 | 39.7 | 2.15× |
| 5.90 | 21.2 | 48.3 | 2.28× |
| 6.70 | 27.5 | 63.0 | 2.29× |
| 7.40 | 25.2 | 58.0 | 2.30× |
| 7.95 | 26.6 | 63.7 | 2.39× |
| 10.2 | 26.1 | 58.8 | 2.25× |
| 12.2 | 26.5 | 55.0 | 2.08× |
| **all** | **24.2** | **54.6** | **≈2.25×** |

Measured 2026-11-14 (2.25–2.27× over that day's runs). Portable build: 1.78×.

## Fuzzing (`just fuzz [target] [seconds]`)
cargo-fuzz (libFuzzer) on nightly-2026-11-01, debug assertions and overflow checks on, ASan on
the Rust side. Corpus kept across runs in `.cache/fuzz/`. `rfc4867`: 2026-11-12, 600 s, corpus
continued. The codec targets diff against the TS 26.073 reference since 2026-11-14:
600 s each, new corpora, then 1 h each.

| target | checks | runs | cov | findings |
|---|---|---:|---:|---|
| `rfc4867` | `.amr`, RTP octet-aligned and bandwidth-efficient readers on arbitrary bytes: no panics, sane error order, accepted frames round-trip through both packings, canonical form never longer | 2,023,256 | 524 | none |
| `decode_vs_ref` | Rust decoder vs the TS 26.073 reference `Speech_Decode_Frame` (`ref26073-sys`), PCM equal frame by frame: every RX frame type (bad, degraded, onset, SID first/update/bad, NO_DATA), any mode and bits; stops at a homing frame (the reference homes in `decoder.c`) | 452,074 (600 s); 30 min: 1; after the fix 1,175,462 (1 h): 0 | 2,162 | defect 12 (DTX mute clamp, inherited from opencore), fixed; both saved crash inputs replay clean |
| `encode_vs_ref` | Rust encoder (DTX off / VAD1, per-frame mode from the input) vs the reference `Speech_Encode_Frame` + `sid_sync` on arbitrary PCM: TX type and serial bits equal | 550,312; 1,015,759 (1 h) | 5,518 | none |

## Robustness: defects found in the legacy C (each absent in the Rust port)
1. **No homing** (TS 26.073 encoder/decoder homing): raw opencore's decoder gets 1/191 on the
   official vectors (its encoder passes them only because their homing frames come first).
2. **`const` contract violation**: `Encoder_Interface_Encode(…, const short* speech, …)`
   overwrites the caller's buffer.
3. **Uninitialised memory read** in `A_Refl` (DTX comfort-noise path): branches on an
   uninitialised stack array for some filters (valgrind-confirmed; ~8% of random valid-LSP filters).
4. **Signed-overflow UB** in `L_add` as used by `gc_pred`: gcc removes the saturation check
   (not reachable by the decoder; the Rust port saturates as specified).
5. **Silent no-op on frame types 9–14**: `Decoder_Interface_Decode` writes no PCM and returns
   `void`, so the caller plays stale or uninitialised samples (test:
   `opencore_leaves_pcm_unwritten_for_unsupported_frame_types`). Rust: `Error::UnsupportedFrameType`.
6. **No length check** (by inspection of `interf_dec.h`/`wmf_to_ets`): `Decoder_Interface_Decode`
   takes no input length and reads the full frame size implied by the header byte, so a
   truncated frame over-reads the caller's buffer. Rust: `Error::Truncated`.
7. **Quality bit ignored**: frames marked damaged (Q = 0, RFC 4867 §4.3.2) decode identically
   to good ones instead of being concealed (test: `opencore_ignores_the_quality_bit`).
   Rust: maps them to `SpeechBad`/`SidBad`.
8. **Oversized shift (UB) in the MR122 codebook gain** (`G_code`): `gain >>= i - 1` with
   `i - 1 >= 32` when the target and the filtered codevector are nearly orthogonal and the
   codevector energy is large. x86 masks the count and returns the gain *unshifted* (e.g. 16,236
   where the denormalised gain is 0), so the result depends on the CPU and compiler. Found by
   the `G_code` differential test on random input (5 of 200,000 cases); not reached by the vectors,
   the 40,000-frame synthetic encoder differential or the fuzz run. Rust: gain 0.
9. **No VAD option 2**: opencore ships `vad2.cpp`/`r_fft.cpp` but neither builds nor wires them
   (they don't compile as shipped: missing `l_deposit_h.h`, operators that exist only `inline`),
   so it can't produce the DTX2 vectors. Rust: `Dtx::Vad2`, 11/11 DTX2 vectors bit-exact
   (oracle for the module tests: those two files, built with a compat shim in `opencore-sys`).

10. **Wrapping where TS 26.073 saturates.** opencore's speed rewrites replace the saturating
    basic operators with plain `int` arithmetic in many places, and its `norm_s`/`norm_l`, `shl`,
    `L_shl`/`L_shr` (counts beyond ±15/31) and `L_mac` (no saturation when `L_var3` equals the
    product) differ from `basicop2.c`. Visible on full-scale input: 10 of the A-law/µ-law
    vectors and 1,598 of the 3,840 frame-sync vectors fail (table above; the reference
    reproduces them all). Found in: `Residu`, `Syn_filt`, `Convolve`, `Pred_lt_3or6`,
    `Interpol_3or6`, `Norm_Corr`, `G_pitch` (also denormalises the scaled `<xn,y1>` by 4
    instead of 2), `G_code`, `q_gain_code`, `cl_ltp`'s target update, `subframePostProc`,
    `Pre_Process`, `Autocorr` (re-sums with wrapping after scaling), `calc_unfilt_energies`,
    `calc_filt_energies`, `calc_target_energy`, `gc_pred`, `cor_h` (scaling threshold, `h2`,
    `mult(32767, 32440)`), `cor_h_x`/`cor_h_x2`, `set_sign12k2`, the codebook searches' compare
    and energy chains, MR122's filtered codevector, `MR475_gain_quant`'s search, `agc` and the
    post-filter's energy chain.
11. **The decoder's overflow re-synthesis is unreachable**: opencore's `Syn_filt` reports no
    overflow, so `Decoder_amr` never takes the reference's "scale the excitation by 1/4 and
    re-synthesise" branch. Rust: `syn_filt` returns the reference's flag.
12. **DTX mute clamps the comfort-noise energy at 0**: entering DTX_MUTE (no SID update for
    50 frames), `dtx_dec` lowers `log_en` by 256 and opencore clamps it at 0 where TS 26.073
    saturates and lets it go negative; a low-energy SID followed by long silence then decodes
    to different comfort noise. The Rust port had inherited the clamp; `decode_vs_ref` found
    it after 30 min of fuzzing and it is fixed (test
    `decoder_matches_ref_through_dtx_mute`; the vectors never reach it).

Also differs from the current reference, not a bug: `dtx_enc` recomputes the SID parameters
while `log_en_index == 0` (an older reference version; the vectors pass either way, the Rust
port follows the current one). And opencore's unrolled loops drop the remainder (`Syn_filt` an
odd last sample, `Residu` the first `lg % 4`, `Pred_lt_3or6` an odd last one; the codec only
uses lengths 22 and 40).

The Rust codec's only `unsafe` is `src/simd.rs` (runtime AVX2 dispatch and SSE2 transposes),
plus FFI in the tests.
