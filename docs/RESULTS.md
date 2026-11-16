# Results

Measured results of amr-nb (the sections below) and amr-wb (its
[decoder](#amr-wb-amr-wb-decoder) and [encoder](#amr-wb-encoder-amr-wb-ts-26190) at the end).
Every number is reproducible with the recipe or command in its heading; the
oracles (3GPP vectors and reference C, opencore-amr) are downloaded and verified on first use
(see the README).

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

## amr-wb (AMR-WB decoder)
The AMR-WB decoder (3GPP TS 26.173, ITU-T G.722.2), ported from opencore-amr 0.1.6's decoder
and corrected to the specification wherever the two differ (below). No `unsafe`.

### Conformance (`just conform`)
The whole decoder set of 3GPP TS 26.174 v19.0.0 (j00, its latest version): every `.cod` with
an `.out`. T00–T20 in each of the 9 modes, T22 mode switching, and the DTX sequences (Dtx1–5 at
23.85 kbit/s, Dtx3 in every mode). Bit-exact = every output sample identical.

| decoder | cases | frames | bit-exact |
|---|---:|---:|---:|
| amr-wb (Rust) | 221 | 119,660 | **221/221** |
| 3GPP TS 26.173 reference C (`ref26173-sys`) | 221 | 119,660 | 221/221 |
| opencore 0.1.6 (C), its frame loop | 221 | 119,660 | 189/221 |

opencore fails the same 32 vectors (T02, T13, T20 in most modes, T05, T08, T14, Dtx2) on the
Rel-11 vectors (b00), so its misses are its own deviations (defects below), not later changes.
Since opencore, TS 26.173 changed the decoder in three places, none visible on the vectors (the
Rel-11 and current reference each reproduce both vector sets): CR 0026 (2007, late SID_UPDATE
and a division by zero in DTX_MUTE), CR 0030 (2015, the comfort-noise energy history's initial
value), both adopted, and CR 0027–0029 (2007), deliberately not adopted (next paragraph). The
other CRs (0031–0036) change the encoder only.

**The one deliberate difference from TS 26.173:** NO_DATA frames are handled as before CR
26.173-0027–0029, as opencore does: a NO_DATA frame during speech counts as DTX for the
hangover logic at once, where the CR first waits for 8 frames without voice activity. The CR's
rule is covered by a live patent (Nokia US 8090588, EP 2201565). It only changes the comfort
noise of the next SID_FIRST after NO_DATA frames that arrive mid-speech (frame loss signalled as
NO_DATA). No vector has NO_DATA mid-speech (asserted by the conformance test).

### Against the TS 26.173 reference (`cargo test -p amr-wb`)
The unit-test oracle is the 3GPP reference C (`crates/ref26173-sys`, built from the pinned
download, symbols prefixed `refwb_`).
- **Whole decoder** on 1,600 stress streams of 120 frames (`just conform` runs 200 per class,
  `cargo test` 40): real frames with flipped bits, random and all-ones bits (full-scale,
  clipping signals), bad, lost and degraded frames, SID_FIRST/SID_UPDATE/SID_BAD and NO_DATA
  sequences, mode switching on every frame. All bit-exact. NO_DATA is generated only in DTX
  (the excluded difference); `no_data_mid_speech_is_the_one_difference` pins it.
- **Kernels**: the basic operators and `Pow2`, `Pred_lt4`, `Syn_filt_32`, `Syn_filt` (both
  orders), and over 40-call sequences `HP50_12k8`, `HP400_12k8`, `Oversamp_16k`, `Filt_6k_7k`,
  `Filt_7k`, `Deemph_32` (outputs and memories). Each exact fast path (below) has an input class
  that drives its guard false; forcing each guard true fails these tests. `Isf_isp` on every
  ISF of its domain [0, 16384] (defect 11).
- **Extreme states**: a crafted stream (speech with the largest ISF 14 the quantiser gives,
  then comfort noise with dithering) drives ISF 14 to 16384, the table's end (defect 11):
  bit-exact (`dtx_isf_at_the_table_end`; it panicked before the fix).

### Packed formats, RFC 4867 (`cargo test --release -p amr-wb rfc4867`)
`.awb` storage and RTP payloads (octet-aligned, bandwidth-efficient; single channel, no
interleaving/CRC/robust sorting): `amr-wb::rfc4867`. Oracles, all green:
- **Every vector through `.awb`**: the 221 TS 26.174 decoder cases packed, read back and
  decoded: **221/221 bit-exact** (`awb_files_decode_like_the_vectors`).
- **Against the TS 26.173 reference's file reader** (`Read_serial`, its `decoder -mime`): it
  reads what we write (every mode × good/bad, SID first/update/bad × 9 mode indications,
  SPEECH_LOST, NO_DATA), and both readers agree on 20,000 random frames (every frame type,
  random Q, padding and mode indications).
- **Against the reference's file writer** (`Write_serial`, its `coder -mime`): 600 encoder
  frames (speech with pauses, DTX on, every mode) written by us and by the reference coder +
  writer are byte-equal (`awb_matches_the_reference_writer`). A SID_FIRST frame's 35
  comfort-noise bits are sent as zeros, as there (10 such frames carry non-zero analysis bits).
- **A SID frame's mode indication is not used**: the decoder takes the previous mode for SID
  frames, as the reference's `.awb` and ITU readers do (`prev_mode`), so a corrupted
  indication can't change the output. The vectors (whose own format carries a mode per SID
  frame) still decode 221/221, also with every SID frame's mode rotated (`sid_mode_is_ignored`;
  it fails if the decoder uses the SID's mode).
- **Round trip** in both RTP packings (1–4 frames per packet), sizes vs opencore's tables,
  header layouts worked out by hand from the RFC, every truncation/reserved-FT error path.
  Frame types 10–13 (reserved, size unknown) end a file and discard an RTP packet.

### Fuzzing (`just fuzz wb_<target> [seconds]`)
As amr-nb's (nightly-2026-11-01, ASan on the Rust side), 2026-11-15.

| target | checks | runs | cov | findings |
|---|---|---:|---:|---|
| `wb_decode_vs_ref` | Rust decoder vs the TS 26.173 reference's frame loop (`ref26173-sys`, homing included), PCM equal frame by frame: every RX frame type, mode and bits; NO_DATA only in DTX (the one difference, kept out; without that rule it fails in ~20 s) | 1,235,729 (1 h) | 2,120 | none |
| `wb_rfc4867` | `.awb`, RTP octet-aligned and bandwidth-efficient readers on arbitrary bytes: no panics, accepted frames round-trip through both packings, canonical form never longer | 589,309 (120 s) | 494 | none |

The fuzzer did not reach defect 11 (an extreme state: a history of maximal ISFs, then
dithered comfort noise); a bound per table lookup over the generated tables did.

### Decoder speed (`CRATE=amr-wb just bench`)
µs per 20 ms frame, best of 15, same machine and builds as amr-nb's tables, on the
single-mode vectors (93,258 frames per rep). Both run the same frame loop (homing test,
decode, 14-bit output) on pre-parsed serial bits.

| mode | Rust | C | speedup |
|---|---:|---:|---:|
| 6.60 | 19.3 | 35.5 | 1.85× |
| 8.85 | 16.4 | 33.6 | 2.04× |
| 12.65 | 15.6 | 30.9 | 1.98× |
| 14.25 | 15.6 | 31.2 | 2.00× |
| 15.85 | 15.6 | 31.8 | 2.03× |
| 18.25 | 16.1 | 32.7 | 2.03× |
| 19.85 | 16.1 | 33.1 | 2.05× |
| 23.05 | 16.4 | 34.2 | 2.08× |
| 23.85 | 18.0 | 37.7 | 2.09× |
| **all** | **16.6** | **33.4** | **≈2.0×** |

Measured 2026-11-15 (two runs: all 2.00× and 2.02×; 8.85 varies 1.90–2.04×). The frame decoder
runs in an AVX2 copy of the same source when the CPU has AVX2 (+4% over the portable build).
Exactness costs: the specification's saturating chains are serial where opencore's plain sums
vectorise (the straight port with the corrections ran at 0.76×). The gains come from exact
reformulations: bounds that prove a chain can't saturate (checked after a fast pass over the
outputs it made), saturations collapsed into one clamp (a clamp is monotone), energies in
closed form, the synthesis filters in blocks with only in-block terms serial, and the
follow-on filters run beside them per block.

### Robustness: defects found in opencore's AMR-WB decoder (each absent in the Rust port)
1. **Wrapping where TS 26.173 saturates** (as amr-nb's defect 10): `HP50_12k8`, `Pred_lt4`,
   the adaptive codebook's low-pass, `Deemph_32`, `Syn_filt_32`, `Syn_filt`, `Isp_Az` and its
   polynomials, the noise enhancer (`Mpy_32_16`, `L_add`), the code's tilt filter (it sums two
   neighbours into one 16-bit argument), `Filt_6k_7k` and `Filt_7k` (which also pre-adds the end
   taps' samples in 16 bits). On loud signals: 32 of the 221 vectors fail.
2. **`round` wraps** for `[0x7fff8000, 0x7ffffffe]` to −32768 (only `MAX_32` is special-cased);
   the same for its `round(x << 1)`.
3. **Shift counts masked**: `shl`, `L_shl`, `L_shr` take counts modulo 16 or 32 where the
   specification saturates them; in `voice_factor` and `agc2` plain C shifts by a negative
   count or one ≥ 32 (UB; x86 masks, `agc2` shifts by 31 instead of halving), and `Pow2`'s
   final shift likewise for loud comfort noise.
4. **`norm_s`/`norm_l` of 0 and negative values**: `norm_s(0)` is 11, not 0 (a zero code gain
   then skews the voicing after lost frames).
5. **`Scale_sig(x, -3)` saturates in 16 bits**: 32764..32767 become 4095, not 4096, in the
   excitation's energy and the high band's level.
6. **The synthesis filters assume `a[0]` = 4096**: `Isp_Az` lowers `a[0]` (adaptive scaling)
   for comfort noise, which the specification's `Syn_filt_32` and `Syn_filt` follow and
   opencore's ignore, so its comfort noise comes out louder after rescaling.
7. **Concealment gains**: the gain histories' pipelined shift leaves `[b1, b1, b2, b3, new]`
   instead of `[b1, b2, b3, b4, new]`, a bad frame's energy history keeps its oldest entry, and
   the mean of the past energies is `Σ/8` instead of `Σ/4`: wrong gains after every bad frame.
8. **DTX mute clamps the comfort-noise energy at 0** (as amr-nb's defect 12), where the
   specification lets it go negative.
9. **No CR 0026**: in DTX_MUTE right after a SID update opencore divides by zero
   (`div_16by16(1024, 0)` returns 0 there), and a late SID_UPDATE mutes.
10. **The quality bit is ignored**: `D_IF_decode` sets `quality = 1` (as its AMR-NB
   interface, amr-nb's defect 7), so frames marked damaged (Q = 0) decode like good ones instead
   of being concealed (pinned: `opencore_ignores_the_quality_bit`).
11. **`Isf_isp` reads past its cosine table** (also the TS 26.173 reference): at ISF 16384 it
   reads `table[129]` of 129 entries, times offset 0, so the value is right but the read is
   out of bounds; in Rust it panicked (index out of bounds). Comfort noise reaches it: after a
   SID_FIRST averages a history of maximal speech ISFs, the dithering pushes ISF 14 past
   16384 and clamps it there. Now indexed exactly (`ind ≤ 127`, `offset ≤ 128`, same value).

## amr-wb encoder (AMR-WB, TS 26.190)
The AMR-WB encoder, ported from vo-amrwbenc 0.1.3 (Apache-2.0) module by module, with the
saturation of the TS 26.173 reference (vo's `vo_*` macros and plain sums wrap where it
saturates), including CR 26.173-0032 (pitch-gain clipping at 6.60/8.85 kbit/s) and -0036
(`c4t64fx` impulse-response scaling), which vo-amrwbenc predates. No `unsafe`.

### Conformance (`cargo test --release -p amr-wb --test enc_conformance`)
The whole encoder set of 3GPP TS 26.174 v19.0.0: every `.inp` with its `.cod` (T00–T21 in each
mode, DTX sequences, mode switching: 221 cases), the frame-synchronisation set (2,880 cases,
each compared at the 2nd frame after its last decoder homing frame) and bit synchronisation
(`bitsync.inp`: per mode, exactly the unshifted homing triplet gives a homing frame), as
TS 26.074 §8 describes for AMR. Bit-exact = every output word (frame type, mode, bits) equal.

| encoder | cases | bit-exact | bit sync |
|---|---:|---:|---:|
| amr-wb (Rust) | 3,101 | **3,101/3,101** | **9/9** |
| 3GPP TS 26.173 reference C (`ref26173-sys`) | 3,101 | 3,101/3,101 | 9/9 |
| vo-amrwbenc 0.1.3 (C) | 3,101 | 3,006/3,101 | |

vo-amrwbenc fails 95, 90 of them also on the vectors from before the two CRs (its `vo_*`
macros and plain sums wrap where the specification saturates). The Rust encoder also passes the vectors in the overflow-checked debug build.

### Against the TS 26.173 reference (`cargo test -p amr-wb`)
- **Whole encoder** on 64 stress streams of 300 frames (`enc_vs_ref`): real speech scaled up to
  clipping, white noise at every level, held tones and sweeps up to full scale, full-scale
  square waves, silence and near-silence, encoder homing frames, speech with random bursts;
  mode switching, DTX on and off. All bit-exact. Classes are held for runs of frames: with a new
  class every frame, the CR 0036 branch (a peaky impulse response from a resonant LPC filter)
  was never reached (its mutant survived; held tones reach it).
- **Modules**, each against its `refwb_*` function: the operators (`Div_32`, `Log2`, `shr`,
  `div_s`, `Mpy_32_16`), `Decim_12k8`, `Scale_sig`, `Preemph`/`Preemph2`, `Deemph2`, `Residu`,
  `Syn_filt` (orders 16/20, in place, both updates), `Convolve`, `Updt_tar`, `LP_Decim2`,
  `HP400_12k8`, `Autocorr`→`Lag_window`→`Levinson`→`Az_isp`→`Isp_isf` on signal classes and on
  arbitrary autocorrelations and A(z) (Levinson's unstable-filter and Az_isp's too-few-roots
  fallbacks), `Qpisf_2s_46b`/`36b` and `Qisf_ns`, `Pitch_med_ol`/`Hp_wsp`/`Med_olag`, `wb_vad`
  (the whole `VadVars` after every frame) and `ilog2` on every input, the encoder DTX
  (`dtx_enc`, `dtx_buffer`, `tx_dtx_handler`, the whole state), `Pitch_fr4`, `G_pitch`,
  `Gp_clip` (CR 0032), `cor_h_x`, `ACELP_2t64_fx`, `ACELP_4t64_fx` (every bit allocation,
  CR 0036's scaling), the pulse coders, `Q_gain2`. Sequences for the stateful ones; a mutant of
  every special case fails.

### Fuzzing (`just fuzz wb_encode_vs_ref [seconds]`)
As the decoder's (nightly-2026-11-01, ASan on the Rust side), 2026-11-16.

| target | checks | runs | cov | findings |
|---|---|---:|---:|---|
| `wb_encode_vs_ref` | Rust encoder vs the TS 26.173 reference's frame loop (homing, `Write_serial` frame types) on arbitrary PCM, DTX on/off, any mode per frame (each frame repeated 1–8 times, scaled down 0–15 bits: DTX and SID frames reached; a DTX-only mutant fails in ~4 min) | 55,720 (30 min) | 2,895 | none |

### Encoder speed (`CRATE=amr-wb BENCH=enc just bench`)
µs per 20 ms frame, best of 15, on the single-mode vectors' input PCM (T00–T20, no DTX). Both
run their frame loop (homing test, 14-bit input, serial output); vo-amrwbenc at -O3, its C path
(no ARM assembly).

| mode | Rust | C | speedup |
|---|---:|---:|---:|
| 6.60 | 205.7 | 108.5 | 0.53× |
| 8.85 | 252.6 | 121.6 | 0.48× |
| 12.65 | 299.5 | 137.6 | 0.46× |
| 14.25 | 331.7 | 145.3 | 0.44× |
| 15.85 | 332.3 | 147.3 | 0.44× |
| 18.25 | 345.5 | 152.2 | 0.44× |
| 19.85 | 364.2 | 153.9 | 0.42× |
| 23.05 | 357.2 | 154.8 | 0.43× |
| 23.85 | 299.4 | 156.5 | 0.52× |
| **all** | **309.8** | **142.0** | **0.46×** |

Measured 2026-11-16 (a 3-rep run before: all 0.46×, 281 vs 130 µs).

The straight port (the reference's saturating chains, no fast paths yet) runs at about half
vo-amrwbenc's speed; the decoder's straight port was at 0.76× before its speed work.
