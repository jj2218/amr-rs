# Results

The scoreboard. Every number is reproducible with the `just` recipe in its heading.

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
3,840 sync cases). opencore's misses are its arithmetic (defect 10), and they were the Rust
port's too until 2026-11-12. The conformance run pins them for the C oracle (exact first
mismatch; per-mode sync counts), so a change either way fails `conform`. Bit sync: of the 13
bit-shifted homing triplets, exactly the unshifted one yields a decoder homing frame, in every
mode and law.

The "via" rows pack every frame (speech, SID_FIRST, SID_UPDATE, NO_DATA) with
`amr_nb::rfc4867`, read it back and decode. Until 2026-11-08 the decoder harness skipped the 11
DTX2 cases (a file-name filter) while this page claimed them; they pass. The encoder vectors
can't show encoder homing (raw opencore also gets 188/188: their homing frames come first), so
a unit test checks it: after a homing frame the encoder continues exactly like a fresh one.


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

2026-11-14, optimisation loop after 39 experiments (SSE2 `syn_fast`, `pred_lt`, `residu`, the
postfilter's impulse response in whole blocks; loop section below): 2.47× → 3.06× (portable 3.02×). Before: after winning back the exactness fixes' cost (the win-back log below the encoder's
log): 2.08× → 2.47×. Before the fixes (2026-11-11, opencore's wrapping semantics): ≈2.51×.

### Optimisation log
Each step kept `just conform` at 180/180 and its module's differential test green. Found with
`just profile` (callgrind; branch mispredicts via `--branch-sim=yes`), measured with `just bench`.
Rust µs/frame (all modes) after each step:

| step | µs | note |
|---|---:|---|
| faithful port | 12.27 | baseline (1.19×) |
| `pred_lt_3or6` as one 20-tap dot product | 11.95 | wrapping sums are associative, so tap order is free |
| `div_s` as one hardware divide | 11.02 | the 15-step restoring division = `(a << 15) / b` |
| branchless `l_add` / `l_mac` saturation | 9.14 | the `&&` overflow tests caused ~55% of simulated mispredicts |
| `l_mac_energy`: exact i64 sum, faithful fallback | 8.28 | agc + excitation energy; fallback keeps opencore's `> 0` quirk |
| `pred_lt` in 8-sample blocks from a snapshot | 7.85 | LLVM's own vectorisation needed alias checks that fail for short lags |
| `residu` in 8/4-output blocks | ≈7.7 | |
| `ph_disp`: contiguous convolution, local `ov` | ≈7.6 | |
| `syn_filt`: y[n-1] in a register | ≈7.3 | the recursion is latency-bound; store→load forwarding was on its critical path |
| parse parameters once per frame | 7.1 | within noise; removes duplicated work |
| `syn_filt` in 8-output blocks | flat | kept for the encoder |
| `ph_disp`: one 40-sample saturating pass per pulse | −1% | 2.03× → 2.04× |
| `pred_lt_3or6`: phase taps as a `const` table, in-place block reads | −2% | 2.05× → 2.09×; stack copy missed store forwarding |
| `preemphasis` and `autocorr_h` as lanes | −1.5% | 2.08× → 2.11×; the C runs backwards, so outputs are independent |
| `post_process` state as one word | −1.6% | 2.14× → 2.17×; hi·a + ((lo·a) >> 15) = (u·a) >> 15 |
| `post_process` fast path while \|y\| < 2^28 | −3% | 2.15× → 2.22×; `l_shl` can't saturate there |
| total excitation as branchless lanes (`total_exc`) | −7.5% | 2.19× → 2.38×; saturating operators as selects |
| `cb_gain_average` LSP deviations as lanes (`lsp_dev_sum`) | −2.5% | 2.38× → 2.44×; `norm_s`/`div_s` exact in floats |
| `bgn_scd` / DTX frame energies as one exact sum | −3% | 2.36× → 2.49× with the step above |

Tried, no gain, reverted: `syn_filt` with windowed pair dot products (latency-bound, same
time); branchless `add`/`sub`/`l_mult` (the loops still didn't vectorise); `post_process` state
in locals / reordered sums (LLVM already does this; the one-word state paid, see the log);
`agc`'s gain recursion run alone up to its fixed point (exact: once gain[n] = gain[n−1] it
stays), then one scaling pass: flat (9.94 vs 9.96 µs); `get_lsp_pol`'s `L_Extract` pair as one
word (`((f >> 1)·x) >> 15`, exact, as `chebps`): Rust −0.1–1% in 6 of 6 pairs, ratios split
(decoder 3 lower, encoder 3 higher), reverted.

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

2026-11-14, optimisation loop after 30 experiments (`pair_safe`, SSE2 `syn_fast`, `pred_lt`,
`residu`, `norm_corr`, `comp_corr`/`autocorr`, the step-5 guard proven dead, MR102's weighted
lag pick, gain searches' first tiers; loop section below): 1.83× → 2.25–2.27× after 39 (this
run's machine state 2.25×; portable 1.41× → 1.78×, the SSE2 VQs). Before: after the win-back (log below): 1.65× → 1.83× with TS 26.073's saturation
everywhere. Before the fixes (main after PR #31, opencore's wrapping semantics): ≈2.20×
(portable ≈1.69×; after the fixes 1.41×, win-back audit below). The log below is the perf phase before the fixes.

### Optimisation log
Each step kept `just conform` at 199/199 encoder vectors and the differential tests green
(new direct C diffs: `cor_h`, `Convolve`, `comp_corr`, `Qua_gain` with a saturating input class,
`Pitch_fr`, `cor_h_tm` vs `cor_h` (5 tracks, and MR102's 4 tracks at both row widths),
`search_10i40`/`search_8i40` vs `search_10and8i40` with fallback-forcing inputs; the per-lane
argmax (`Pair::best`) vs the compare chain on tie-heavy input; `pm_scale` exhaustively;
`Test_Vq_subvec3/4` on all six codebooks, `MR475_gain_quant`, `Autocorr` at the 2^31 energy
boundary; `filtered_code` vs the `L_mult`/`L_mac`/`round` chain and MR122's lane VQs vs the C's
early-exit loops, `MR795_gain_quant`, `set_sign12k2` (incl. its only `L_mac`-saturating input),
`Chebps` at every x, `calc_filt_energies`; each with saturating/full-scale and tie-heavy
classes and a mutant that must fail). Gains: speedup vs the C from interleaved `AB` runs (all
modes unless one is named; early steps in µs/frame). Kernel timings and the full rationale are
in each step's commit message (`git log -S '<kernel>' -- docs/RESULTS.md`).

| step | gain (AVX2; portable) | why it paid |
|---|---|---|
| faithful port | 1.21× | baseline |
| `convolve` in 8-output blocks | 46.2 → 45.0 µs | wrapping sums: order is free |
| `cor_h` one sweep over all diagonals | 45.1 → 44.7 µs | 8-lane chunks |
| `comp_corr` compile-time frame length | 45.0 → 43.9 µs | one `pmaddwd` reduction per lag |
| `qua_gain` exact i32 lanes + overflow flags | 45.0 → 43.4 µs | faithful fallback on any saturation |
| AVX2 dispatch (`simd.rs`): `cor_h`, `convolve`, `comp_corr`, `vq_subvec4` | 1.35× → 1.38× | `#[target_feature]` copies of one source |
| MR122 `cor_h_tm` + `search_10i40` per-lane argmax | 12.2 1.26× → 1.34× | track-major `rr`, candidates as 8-lane rows |
| `cor_h_tm` sign scaling by lane width | portable 12.2 1.14× → 1.22× | portable: multiply-free ±32767 form in i16 lanes |
| step-5 searches walk contiguous track-major rows | 1.40× → 1.42× | `next_pulse` over `cor_h_tm` |
| step-5 searches in lanes (`step5_block`) | 1.42× → 1.55×; 1.32× → 1.38× | 8 first pulses = 8 exact chains, 8×8 transposes |
| MR102 + MR122 generic pair search (`search_pairs`) | 1.54× → 1.61×; 1.39× → 1.42× | `Tm<T, N, W>`, 4-track `cor_h_tm` |
| `syn_filt` in 8-output blocks | 1.59× → 1.63×; 1.42× → 1.43× | vector `past_terms`, serial in-block terms, exact fallback |
| `norm_corr` helpers not inlined | 1.62× → 1.65×; 1.43× → 1.46× | inlined they were 80 scalar `imul`s per lag |
| `cor_h_x` in 8-output lanes (`corr_hx`) | 1.64× → 1.67× | the C's triangle was 359 scalar `imul`s |
| `vq_subvec3/4` over component-major codebooks | 1.69× → 1.74×; 1.46× → 1.48× | one load per component, no gathers |
| `MR475_gain_quant` in lanes | 4.75 1.55× → 1.89×; 1.39× → 1.51× | no `Mpy_32_16` there can saturate: no fallback |
| `filtered_code` 8-output blocks, `L_mac` as selects | all 1.72× → 1.79× (with MR475) | sign-bit overflow test vectorises |
| `autocorr` energy as one exact sum | 1.78× → 1.80×; 1.49× → 1.50× | "first negative sum" = "total ≥ 2^31" |
| MR122 `Vq_subvec(_s)` in lanes | 12.2 1.57× → 1.69×; 1.42× → 1.44× | early exits = first minimum of full distances |
| `Qua_gain` component-major, dispatched | 1.81× → 1.87× | AoS rows had blocked lanes |
| `MR795_gain_code_quant3` in lanes | 7.95 1.71× → 1.98×; 1.42× → 1.52× | saturating input declines to the C's loop |
| MR122 build code as shifted slices of `h`; `set_sign12k2` as selects | 12.2 ~1.65× → ~1.9×; 1.44× → 1.58× | `(h·±8192) >> 7` = ±64 h exactly |
| `chebps` b1 as one word | 1.96× → 1.98×; 1.55× → 1.56× | t0 even: `L_Extract` pair = `(t0·x) >> 16` |
| `calc_filt_energies` plain sums under a bound | 1.96× → 2.01×; 1.55× → 1.59× | 2 max\|x\| Σ\|tmp\| < MAX_32: no prefix saturates |
| `calc_unfilt_energies` (MR795) | 7.95 2.04× → 2.06× | residual as selects, s4 exact sum |
| `pred_lt_3or6` phase-tap table + `pre_process` one-word state | 2.01× → 2.03× | see the decoder log |
| MR475/MR515 only the `rr` blocks read (`cor_h_blocks`) | 5.15 1.85× → 1.95×, 4.75 1.98× → 2.07×; 1.52× → 1.70× | 2 of 25 blocks, one residue each |
| MR122/MR102 `Pair::best` AVX2 form (`best_i32`) | 12.2 1.94× → 2.21×, all 2.07× → 2.08× | masks at the multiply keep `pmaddwd`; f64 key argmax |
| MR74/MR795 `search_4i40` packed halves, SSE2 transposes | all 2.10–2.14× → 2.18–2.23×; 1.65× → 1.73× | 4 of 8 lanes valid: two blocks per pass |

Tried, no gain, reverted:
- `syn_fast` (per 8-output block, timed without the fallback: 23 ns): pairs of outputs in one
  i64 (SWAR: y·(a_k + 2^32 a_k+1) gives two outputs' terms in one `imul`, exact while
  Σ|a_1..a_8| ≤ 15.0 in Q12, 99.93% of the codec's calls): 44 → 24 `imul`s, 29 ns; the a1 term
  last on the chain (imul → sub → sar): 156 vs 150 ns per call; `old` read in place from
  `mem ++ y`, or passed in registers (two packed u64, no narrow-store/wide-load stall): flat;
  breaking the block-to-block dependency on `old` entirely (timing hack): flat, so it isn't
  latency between blocks; `past_terms` over a pair-interleaved `col` (`pmaddwd`'s shape):
  scalarised (72 `imul`s), 2× slower. llvm-mca (`llvm-mca-14`) put both loop bodies at ~60
  cycles; the measured ~100 includes `past_terms` (~130 instructions).
- `syn_filt2`: `subframe_pre_proc`'s two independent `syn_filt` pairs (same `aq`, same `ap2`)
  interleaved in one call (`syn_fast` generic over streams): Rust 40.0 vs 40.0 µs. `syn_fast`
  is bound by scalar `imul` throughput (44 per block of 8 on one port), not by chain latency.
- `pred_lt_3or6`'s block with tap pairs in `pmaddwd`'s shape (interleaved windows, `from_fn`
  rows, chunked taps): every form went fully scalar (389 vs 155 ns); its AVX2 copy uses
  `vpmulld` (162 vs 155 ns).
- `cor_h` reading the accumulators backwards, and as 16-bit-lane products (`pmulhw` form): slower.
- `convolve` as 40 full-length dot products: 2× slower than the faithful loop (458 ns).
- `comp_corr` with 4 lags per pass sharing x (the C's shape): 5920 ns, loses `pmaddwd`.
- `norm_corr`: snapshot update loop, fixed-length/split reductions: no change.
- `vq_subvec3/4` with 8 codewords in lanes over the AoS codebook: no gain (stride-3/4 gathers on
  SSE2; done since over component-major tables). With them: the product `(w*x) >> 15` from
  `pmulhw`/`pmullw` halves (LLVM kept `vpmulld`; 1.4× slower), components squared in pairs as one
  `pmaddwd`: no gain; 8 lanes: erratic (vq4 1.4× slower).
- `filtered_code` with all 40 accumulators in one lane array: 261 vs 180 ns (accumulators
  loaded/stored per pulse, `pmuludq` multiplies); blocks of 8 with `overflowing_add`: scalarised
  (`seto` per lane), 256 ns. The sign-bit overflow test vectorises: 87 ns.
- `syn_fast` blocks written in place (no zero-padded x/y copies, `memset`/`memcpy` gone; tail
  through a padded block): 193 → 213 ns (A/B). Its columns as overlapping slices of a
  zero-padded `a` (no 8×8 `col` per call): `past_terms` went scalar (64 `imul`s), 210 ns.
- `next_excf` as a plain offset loop: already vectorised (its 7 scalar `imul`s are the
  39-element remainder), `pitch_fr` flat.
- `az_lsp` evaluating both polynomials on the whole grid in lanes (122 lane evaluations instead
  of ~61 serial `chebps` calls): 1151 → 1184 ns.
- `MR475_gain_quant` lanes with `Mpy_32_16`'s saturation as selects: unreachable there (see
  the log), dropped; 16 lanes on AVX2 spilled (1020 vs 680 ns with 8).
- AVX2 copies (same source, `#[target_feature]`) with no gain, min of 3 interleaved `BENCH=kern`
  runs: `pitch_fr`, `vq_subvec3`, `set_sign`, `qua_gain`, the per-mode searches
  (`code_*`, `search_10and8i40`; `cbsearch` within ±2%): flat. `cor_h_x` +11%, `filtered_code` +36%.
- `next_pulse` branchless selection: within noise. Track-major `rr` copy + tree argmax
  (exact while all `alp > 0`, `sq >= 0`): search −15%, but the regrouping costs 374 ns; net noise
  (done since by writing the layout in `cor_h_tm`).
- `cor_h_tm` shapes (AVX2): lower triangle + mirror pass 1093 ns; all 79 offsets, row windows
  read at moving unaligned offsets 1105–1552 ns (store-forwarding stalls, index math per row);
  only the live 8 lanes per residue, shifted in memory 833 ns, residue loop outside 800 ns;
  unrolled by 5 with constant indices 520 ns (kept). Portable: the multiply-free sign form as
  an `i16` select 5300 ns (not vectorised; as an `i32` mask 1115 ns, kept for SSE2 only: 664 vs
  520 ns on AVX2), `(x*y)>>15` from `pmulhw`/`pmullw` halves: no change.
- `search_10i40`: candidates fused into the argmax (rows made on demand, no 64-entry buffer):
  3817 vs 2127 ns (the closure defeats vectorisation; done since with plain loops: `Pair::best`,
  kept); branch-free lane reduction: +2% (noise).
  On the vectors the exact-chain fallback is never taken (0 of 800,000 picks).
- `cor_h_tm` on SSE2 (16 xmm registers: the 5 residues' sums and `hd` spill): only the 15
  blocks ta <= tb in the steps + 10 8×8 transposes (unpack rounds) 922 vs 899 ns, AVX2 576 vs
  496 (the sums still advance every step, only the scaling is saved); residue loop outermost
  1088 ns; the product as a `pmaddwd`-shaped pair sum (`hd` interleaved with zeros) 1427 ns
  SSE2, 591 AVX2. The i16 ±32767 form on AVX2: 807 vs 496 ns (the multiplies stay there).
- `next_pulse` making the 8 candidates as vector rows (i32 lanes) before the chain: no faster
  than the scalar loop on contiguous rows (the chain extracts lanes, `pextrw`). Superseded by
  lanes across the first pulse (`step5_block`, kept).
- `lanes_next` written per element (one `for l` body with `if` selects, a `bool` overflow
  flag, an `i32` clamp for `add`): not vectorised at all (scalar `imul` + branches; MR59 +5%,
  the rest gained only from independent chains). Whole-row passes (`from_fn` per operation,
  i16 mask selects, per-lane flags, `saturating_add`) give `pmaddwd` lanes: MR59 −40%.
- MR102/MR122 searches: the square `rr` regrouped into the 4-track layout after `cor_h`:
  `cor_h` 559 → 1127 ns, 10.2 1.42× → 1.37× (writing it in `cor_h_tm`: 700). Candidates
  stored as rows, then `pick`: `search_8i40` no faster than the C until the energies were i16
  (`vpmulld` → `vpmaddwd`). W = 16 with only lanes `..N` read: LLVM drops the dead lanes and
  splits rows 8 + 2 (1220 ns; flagging all 16 lanes keeps whole vectors: 920). Cross-lane
  argmax as a branch-free halving tree: slower (1490 → 1660, 1220 → 1298 ns). W = 16 on the
  portable build: `cor_h_tm` 2090 ns (16 i32 lanes × 4 residues don't fit in 16 xmm), 10.2
  1.34× → 1.26×; W = 10 with AVX2: search 1070 vs 920 ns.
- Kernel timing inputs: the searches' random `h` and signs (not ±32767, full-scale `h`) make
  energies < 1 and took the exact-chain fallback on 6–21% of picks (0 on the vectors), hiding
  the pair searches' gain; the timing class now uses the encoder's kind of `h` and signs (<1%).
  `cbsearch`'s inputs (`cb_in`) still do it on ~12% of MR102 picks: its numbers understate.
- `pitch_ol`'s energy (stop at the first negative wrapped sum = "exact total ≥ 2^31") as one
  u64 sum: 62–70 vs 72 ns per call (the old loop is ~1 cycle per sample; the u64 widening eats
  the lanes), under 0.1% of a frame; `pitch_ol_wgh`'s with `l_mac_energy`: reverted with it.
- `Pair::best` (MR122 `search_10i40` AVX2 kernel, 1260 ns before): two chains per lane in
  i16 state (rows 0..4 / 4..8, merged): 1240 (LLVM packs i32 compare masks into i16 blends,
  ~120 instructions per row); i16 state compared as `pmulhw`/`pmullw` (high, unsigned low)
  halves: 1830 (SLP mixed lane groupings, 2-lane i64 masks); in i32 lanes two chains paid only
  while the products were `vpmulld` (890 vs 960; with `pmaddwd` one chain: 830); interleaving
  the iterations with the old form: flat on both paths. The i32 form on SSE2: 1.07× (no
  `pmulld`, two xmm per vector, spills), the f64 pick there 1.05× (MR122) to 1.18× (MR102,
  scalar `cmov` max: no `pcmpgtq`): portable keeps the i16 form.
- `cor_h_blocks` with the residue a runtime value (one generic step loop): 120 ns per block,
  `cor_h_tm` 0.64× (const residue + dispatch: 95 ns per block, 0.55×). Its portable build with
  `cor_h_tm`'s multiply-free ±32767 scaling: 5.15/4.75 +1% (3 of 3 each; 8 rows per block).
- `syn_filt`: the 8-lane history window in registers (`[i16; 8]` or `[i32; 8]`, shifted per
  sample): LLVM scalarises it (11 `imul`s per sample, register rotation), 158–185 vs 180 ns.
  The block form with its vector pass inlined: 52 scalar `imul`s per block, no faster (SLP has
  no seed: the lanes feed scalar chains); not inlined it vectorises. Blocks of 4: 1.07× the old
  time (more calls, the boundary more often); unpadded tail handling: 0.95 vs 0.90 (padded).
  An AVX2 copy: no gain over SSE2 (dropped).
- `norm_corr` helpers inlined again: `pitch_fr` 611 → 840–890 ns (scalar again).
- `search_4i40` packing: `transpose8` out of line (to stop SLP's element loads): slower than
  inlined (call, spills; the intrinsics fixed it). `lanes_next` with per-lane `d`/`diag` rows
  for every caller: MR59 +4% (the one-track callers keep broadcast rows, `Rows`).

Per-mode `cbsearch` kernel timings can mislead by ±5% between builds: the same search
machine code showed 7.40 +5% in the portable build's kernel run while that mode's interleaved
`AB_MODE=74` (now `BENCH_MODE`) line went −1.5%. Decide with the encoder A/B lines.

Codec facts behind the parked leads: the codebook searches (`cbsearch`, 3.4–6.3 µs per
subframe, ~40% of a frame) are serial compare chains; `set_sign12k2`'s `L_mac` saturates only
when Σ dn² = 2^31 exactly; `Az_lsp`'s diff can't see low-bit changes in `Chebps` (a direct
diff at every x can). General lessons from both optimisation logs (measurement, exact
fixed-point equivalences, vectorisation) are kept with the porting method, outside this page.

## Win-back log (README Next #1, after the exactness fixes)
Every experiment, kept or not, so none is retried. Ratios vs opencore from interleaved
`AB=<prev> just bench` runs (3 pairs; BENCH=enc / dec), start 2026-11-12: encoder 1.65×,
decoder 2.08×. Attribution (`just profile`, 8eafc9c → main): +42% encoder instructions,
spread over the new guards and saturating operators, not the exact fallbacks: on the
vectors those run on ≤ 0.3% of calls (`norm_corr` 0.3%, `syn_filt` 39 of 714k decoder calls,
the others ≈ 0).

| step | encoder | decoder | kept |
|---|---|---|---|
| guards' max\|v\| / sum\|v\| as i16 min/max and i32 sums (`max_abs`, `sum_abs`; i64 abs/max stayed scalar), `norm_corr`'s bound by Cauchy–Schwarz (<x,x><e,e> < 2^60) | 1.65× → 1.70× | 2.08× → 2.12× | yes |
| `g_pitch`'s chains as one i32 dot under 1 + 2 max\|a\| sum\|b\| <= MAX_32; `chebps` guarded by the OR of \|b\| < 2^29 − 2^17 (one OR per step, off the chain; per-step check as the fallback, its limit tightened by 2^16 for the C's hi-first subtraction) | 1.70× → 1.73× | – | yes |
| `subframe_post_proc`'s total excitation (`total_exc`) as branchless selects: round(sat(t << k)) clamps where t leaves [i32::MIN >> k, (i32::MAX − 0x8000) >> k] | 1.73× → 1.77× | – | yes |
| `get_lsp_pol` in one word per coefficient (`Mpy_32_16` of the `L_Extract` pair = 2(((f >> 1) x) >> 15)), exact while the OR of magnitudes < 2^30, else the operator chain | 1.77× → 1.78× | 2.14× → 2.22× | yes |
| one `total_exc` (`ph_disp.rs`) for decoder, `ph_disp` and encoder: guarded fast form (pf, g ≠ −32768, shift 0..=2: exact i32 t, two clamps as selects), dead overflow output dropped, `#[inline]` instead of `#[inline(never)]` (encoder's own copy deleted) | 1.78× → 1.77–1.78× (with the flag and `inline(never)`: −1%) | 2.22× → 2.38× | yes |
| `syn_filt` guard: coefficients-only first tier ((\|a0\| + sum\|a_j\|) 2^15 < 2^30: every sample is 16-bit), max scans only past it; new test class: a chain saturating mid-way with all outputs 16-bit (the guard's old bound loosened 8× had passed the tests) | 1.78× → 1.78–1.79× | 2.40× → 2.45× | yes |
| `agc`: the gain recursion alone (scalar; fac >= 0, g0 >= 0 keep gains in −32767..=32767: `mult` exact, `add` a min), then the scaling as lanes (`L_shl(.., 3)` = clamp of v·gain at ±2^27) | – | 2.43× → 2.49× | yes |
| `pair_safe` maxima in i32 lanes (weights / 256, exact), the 8×8 cross max as i16 min/max. First judged on the "all" line: flat (1.76–1.78×); per mode (`BENCH_MODE`): 12.2 1.53× → 1.63×, 10.2 flat. A one-mode change needs its mode's line | 12.2 1.53× → 1.63× | – | yes |
| `next_excf` plain 16-bit lanes when \|temp\|·max\|h\| can't reach the `L_shl` clamp, max\|h\| per call | `pitch_fr` kernel 1.07× slower (random full-scale inputs: mostly the clamped path) | – | no |
| the same with max\|h\| once per `norm_corr`, judged end to end | 1.78–1.80× → 1.80–1.81× | – | yes |
| `lanes_next` (step-5 searches) saturation bound in u32 lanes (xa clamped at 2^(31 − sc), \|alp0\| as the room) instead of u64 | all 1.80–1.81× → 1.81–1.82×; 7.40 flat, 5.90 1.83× → 1.85× (marginal) | – | yes |
| `set_sign12k2`'s per-sample chain as `sign_cor` (the `total_exc` clamp-select form, k = 11; its clamps tested directly, unreachable through the codec's normalised gains) | 12.2 1.63× → 1.66× | – | yes |
| `cor_h_x2`'s normalisation `round(L_shl(y, j))` as selects (`round_shl`) | flat by the ratio (1.82× → 1.81–1.82×; Rust −1%, the C as much) | – | no |
| `l_mac_energy` in closed form (min(init + 2 sum x², MAX_32), flagged past it: −32768 alone passes it; no fold) with squares summed in u32 pairs into u64 | flat (1.82–1.84× both) | flat (2.48–2.51× both) | no |
| whole-frame AVX2 dispatch (`Encoder::encode`'s body as a `kernels!` entry) | flat (1.82× → 1.81×): `CodAmr::encode`, `cbsearch`, `gain_quant` stay out-of-line baseline functions, so the copy reaches almost nothing | – | no |
| `pred_lt_3or6`: one guard scan over the whole input when t0 >= len + 10 (no output read back) instead of one per block (Σ\|taps\| is 53.5k–61.7k per phase: a coefficients-only tier can't decide) | Rust flat (ratio 1.83× → 1.81×: B's C ran faster) | 2.46× → 2.48× (marginal) | yes |
| (fix, found while reading `syn_fast`'s range check) `syn_filt` reported no overflow for s ∈ [−2^27, −2^27 + 2048), where the reference's `L_shl(2 s', 3)` saturates (same output, −32768) and the decoder re-synthesises; the fast range now starts at −2^27 + 2048, and the `syn_filt` tests compare the flag (the old bound fails two of them) | – | – | fix |
| `syn_fast`'s range flag once per block (min/max of the 8 outputs, conservative) instead of three scalar ops per output | – | flat (2.43–2.49× both) | no |
| `g_code` (MR122) via the exact helpers (`l_mac_dot`, `l_mac_energy`) instead of serial `l_mac` folds | 12.2 1.66–1.67× → 1.67–1.69× (marginal) | – | yes |
| `calc_filt_energies`' MR795/MR475 `<xn2, y2>` chain via `l_mac_dot` | 7.95 flat (1.92–1.96× both) | – | no |

**Audit of the pre-fix optimisations** (2026-11-13, perf/audit). Each guarded fast path vs its
plain exact form (the guard forced to the exact fallback, `#[cold]` dropped), every mode from
one interleaved run (`AB=HEAD BENCH_MODE=each just bench`: median of 3 per mode; B/A < 1: the
fast path pays). A/A noise: ±2% per mode.

| fast path (B = exact form only) | encoder B/A per mode (475 515 59 67 74 795 102 122) | decoder B/A | kept |
|---|---|---|---|
| `convolve` (8-output blocks under max\|x\| sum\|h\| < 2^30) | 0.62 0.67 0.72 0.75 0.74 0.74 0.75 0.74 | – | yes |
| `comp_corr` (plain sums under 2 sum s² <= MAX_32) | 0.56 0.51 0.56 0.63 0.60 0.62 0.62 0.62 | – | yes |
| `residu` (8/4-output blocks under sum\|a\| max\|x\| < 2^30) | 0.90 0.91 0.89 0.94 0.93 0.94 0.94 0.94 | 0.87–0.89 | yes |
| `syn_filt` (`syn_fast` blocks + two-tier guard) | 0.53 0.56 0.61 0.66 0.65 0.66 0.67 0.66 | 0.43–0.46 | yes |
| `pred_lt_3or6` (8-output blocks under the window guard) | 0.84 0.90 0.91 0.94 0.93 0.93 0.93 0.94 | 0.69–0.71 | yes |
| `cor_h` (`cor_h_tm`/`cor_h_blocks` under `rr_may_saturate`) | 0.61 0.56 0.68 0.73 0.71 0.73 0.82 0.73 | – | yes |
| `cor_h_x` (`corr_hx` lanes under max\|x\| sum\|h\| < 2^30) | 0.88 0.85 0.88 0.90 0.89 0.90 0.90 0.90 | – | yes |
| `qua_gain` (i32 lanes + overflow flags; MR515–MR102) | 1.00 0.91 0.92 0.88 0.87 1.00 0.88 1.00 | – | yes |
| `calc_filt_energies` (plain s2, s3 under 2 max\|x\| sum\|tmp\| < MAX_32) | 0.99 0.99 0.99 1.01 0.99 0.99 1.00 1.01 | – | **no: simplified to the chains** |
| step-5 searches (`step5_block`/`step5_halves` lanes, `step5_exact` fallback) | 0.94 0.93 0.79 0.56 0.58 0.59 1.01 1.00 | – | yes |
| `search_pairs` (MR102/MR122 under `pair_safe`) | 1.01 1.03 1.03 1.02 1.02 1.00 **0.65 0.53** | – | yes |
| `filtered_code` (8-output blocks, `L_mac` as selects) vs the scalar `l_mac` chain | 0.98 0.98 0.97 0.99 0.97 0.98 0.94 0.92 | – | yes |
| `MR795_gain_code_quant3` lanes | 1.00 1.02 1.00 1.02 1.01 **0.91** 1.01 1.00 | – | yes |
| `MR475_gain_quant` lanes | **0.73** 1.01 1.01 1.01 1.01 1.01 1.02 1.01 | – | yes |
| MR122 `Vq_subvec(_s)` lanes vs the C's early-exit loops | 1.00 1.01 1.00 1.01 0.99 1.00 1.01 **0.93** | – | yes |
| `post_process` (fast path while \|y\| < 2^28) | – | 0.96–0.99 | yes |

Bold: the only mode(s) running a mode-specific path (the others within ±3%: layout noise).
Only `calc_filt_energies`' guarded s2/s3 dots no longer paid (two 40-term dots per subframe,
≲ 1% of a frame, and the guard scans about as much data); s1 (terms >= 0) stays a plain sum
without a guard. The simplified code vs main: 0.98 1.00 0.98 1.00 0.99 0.99 1.00 1.00 (all
0.992): within noise, so possibly ≲ 0.5% lost, below what this bench resolves.

Portable build (`RUSTFLAGS="--cfg amr_nb_portable"`, one run after the audit): encoder
1.41× (4.75 1.41, 5.15 1.48, 5.90 1.44, 6.70 1.37, 7.40 1.43, 7.95 1.51, 10.2 1.37, 12.2
1.31; ≈1.69× before the saturation fixes), decoder 2.50× (2.46–2.53 per mode).

## Optimisation loop (perf/loop, open-ended)
Started 2026-11-13 from encoder ≈1.83×, decoder ≈2.47× (portable 1.41× / 2.50×). One row per
experiment, kept or not. Per-mode columns: B/A of the speedup vs opencore from
`AB=<previous kept state> BENCH=enc|dec BENCH_MODE=each just bench` (median of 3; > 1: B faster;
A/A noise ±2%),
modes 475 515 59 67 74 795 102 122.

| # | what | encoder B/A | decoder B/A | kept | why |
|---|---|---|---|---|---|
| 1 | `pair_safe` (MR102/MR122 pair-search guard, 12.6% of a 12.2 frame by instructions): the block's max\|c\| as one flat i16 min/max (row/column folds were gathered column-wise: 66 `movzwl`), weights as shifts (all powers of two; SSE2 has no `pmulld`); 552 → 308 instructions. Test: new `search_in` class (slowly decaying `h`, one sign, small `dn`): a guard-always-true mutant passed the old tests, fails now | 1.002 1.006 1.003 1.009 1.006 0.997 **1.067 1.038** (all 1.018) | – | yes | the guard read its block column by column |
| 2 | `pair_safe` first tier from per-block maxima (25 block and 5 diagonal max\|rr\| once per search; the bound's sums as scalars over the prior pulses' tracks), the per-position check only when it fails. Its `a2` half is untested by any input class (as before: it nears 2^31 only with every term at full scale), the `srrv` half is (mutant) | 1.008 1.004 1.039 1.000 1.002 1.001 **1.031 1.037** (all 1.016) | – | yes | 16 per-position checks per subframe → one table |
| 3 | `syn_fast` as SSE2 blocks (`simd::syn_blocks`): only the y[n−1], y[n−2] terms stay scalar `imul`s on the chain (16 per block of 8, was 44); each output's k ≥ 3 terms in the block's later outputs are one `pmaddwd` per half off the chain, read by output m right after output m − 3's update (reading the latest accumulators put a `pmaddwd` on the chain: 1.2× slower); the x term (x, 1)·(a0, 0x800) and the 10 earlier outputs as 6 `pmaddwd` pairs; outputs and history stay in registers across blocks (as element stores + vector reloads: no gain, store forwarding). Over arrays, LLVM spilled (2× slower). Kernel 150 → 120 ns. Oracle: `syn_blocks_is_the_array_form` + the `syn_filt` reference tests; two mutants fail | 1.061 1.036 1.014 1.024 1.042 1.047 1.025 1.030 (all 1.036) | 1.088 1.088 1.085 1.092 1.106 1.073 1.106 1.087 (all 1.092) | yes | multiplies off the chain; per-output state in registers |
| 4 | `syn_tab`'s rows as windows of a zero-padded coefficient array (instead of a range check per entry) | – | – | no | identical machine code (LLVM already built it so); kernel 1.00 |
| 5 | `syn_blocks`: the next block's history pairs from each output pair's own broadcasts (not from the packed outputs, which all waited for y7), and the pair (y6, y7) of the previous block (k = 1, 2 for output 0: zero coefficients) applied after output 0's lane is read. Kernel 120 → 94 ns; a history-rotation mutant fails | 1.018 1.025 1.022 1.016 1.017 1.007 1.015 1.004 (all 1.016) | 1.016 1.037 1.029 1.044 1.035 1.030 1.043 1.033 (all 1.033) | yes | a zero-coefficient vector term still carries its operand's latency |
| 6 | `comp_corr_n` `#[inline(always)]` (the 5.15 profile listed it outside the AVX2 copy) | – | – | no | kernel 1.01: `just profile` builds without inlining, so its symbols don't show what the dispatched copy contains |
| 7 | `post_process_fast`: the state's `<< 2` folded into the coefficients ((4y·a) >> 15 = (y·4a) >> 15 under the range check), first two samples peeled from the state words | – | 0.993 1.008 1.014 1.006 1.005 1.017 1.000 1.040 (all 1.008) | no | within noise (LLVM likely folds it already) |
| 8 | `pred_lt_3or6`'s block as SSE2 intrinsics (`simd::pred_block`): tap pairs as `pmaddwd` (windows t, t + 1 interleaved against the broadcast pair), 10 per half instead of 20 `pmullw`/`pmulhw` + unpack rounds; the array form stays as the non-x86 body and oracle (`pred_block_is_the_array_form` incl. the all −32768 wrap, a tap-pair mutant fails; review found the non-x86 alias naming a private module: `just lint` now type-checks aarch64). Kernel 0.92 (random input, mostly the exact path) | 1.012 1.024 1.022 1.008 1.004 1.005 1.017 1.014 (all 1.012) | 1.028 1.055 1.052 1.032 1.030 1.046 1.041 1.039 (all 1.036) | yes | the pair shape needs intrinsics (arrays went scalar) |
| 9 | `residu`'s 8-output block as SSE2 intrinsics (`simd::residu_block`): 11 taps as 6 `pmaddwd` pairs per half (the last against zeros); the array form stays as the non-x86 body and oracle (`residu_block_is_the_array_form`, full scale too; a swapped-pair mutant fails) | 1.030 1.023 1.025 1.020 1.006 1.020 1.012 1.015 (all 1.020) | 1.013 1.014 1.023 1.025 1.021 1.027 1.031 1.025 (all 1.020) | yes | same pair shape as `pred_lt`; every mode up |
| 10 | `convolve`'s fast path as SSE2 tap pairs (`pmaddwd`, the 20 x pairs broadcast once per call, `packs` for the clamp) | 1.008 1.010 1.017 1.004 1.013 1.003 1.010 1.012 (all 1.009) | – | no | within noise: the AVX2 copy's auto-vectorised 8-lane form was already close (kernel 0.86 on a fallback-heavy input) |
| 11 | `syn_filt2`: `subframe_pre_proc`'s two independent filter pairs (aq: h1 and error; ap2: h1 and xn) as one call each, blocks of both streams in lockstep sharing the table (the SSE2 block step as a per-stream state) | 1.005 1.006 1.017 0.993 1.006 1.019 1.008 1.012 (all 1.005) | – | no | flat: with 16 chain `imul`s the chain is latency-bound, yet a second stream didn't fill it (the core already overlaps the calls, or per-call setup dominates) |
| 12 | `norm_corr`'s lags in SSE2 (`simd::norm_lags`): the filtered excitation stays in 5 registers from lag to lag; range test, both dots (`pmaddwd`, one shared reduction), the update from `pmulhw`/`pmullw` halves + `paddsw` over the one-lane shift; any lag off the fast paths sends the call to the old loop. Oracle `norm_lags_is_the_array_form` (verdict and sums; a dropped lane-0 insert fails; `paddsw` → `paddw` is an equivalent mutant: the next range test rejects either) | 1.052 1.066 1.052 1.037 1.025 1.041 1.035 1.038 (all 1.042) | – | yes | two out-of-line passes through memory per lag → registers |
| 13 | `syn_blocks`' range test on the block's 8 sums at once (two `pcmpgtd` per half, gathered by the unpacks the outputs need anyway) instead of 4 scalar ops per output | 1.019 1.009 1.049 1.016 1.007 1.003 1.003 0.999 (all 1.014) | 0.997 1.008 1.005 0.989 0.989 1.020 0.987 1.005 (all 1.002) | no | flat: the scalar checks ran beside the chain for free (kernel 0.97) |
| 14 | `lag_max`'s argmax (C: from lag_max down, ties `>=`) as a max reduction + first-position scan (a tie mutant fails the unit test and the vectors) | 0.998 1.003 1.005 1.014 0.992 0.999 0.994 1.005 (all 1.003) | – | no | flat (~124 lags per frame); the C-shaped loop stays |
| 15 | `comp_corr` in SSE2 (`simd::comp_corr_lags`): 4 lags per pass sharing each x load, one transpose-add reduction per 4 lags. Auto-vectorised, LLVM had kept x sign-extended to i32 across the lag loop and re-packed it every lag (`vpackssdw` + `vextracti128` per 8 samples, 128-bit `vpmaddwd` only). Kernel 1328 → 645 ns; its AVX2 dispatch dropped (same SSE2 body on both paths); oracle `comp_corr_lags_is_the_array_form` incl. all −32768 and lag_min 18 (a lane-order mutant fails) | 1.029 1.034 1.014 1.017 1.003 1.003 1.018 1.015 (all 1.017) | – | yes | a hoisted sign-extension cost 3 ops per 8 samples per lag |
| 16 | pulse position → (track, index) as a const table instead of q % 5, q / 5 (`Tm::row`, `pair_safe`'s first tier; callgrind's heaviest lines in the pair search after the products) | 0.994 0.971 0.966 0.986 0.997 1.006 0.987 1.011 (all 0.988) | – | no | slower: division by a constant is a cheap multiply-shift off the critical path; the table load adds latency to the address. callgrind's instruction counts misled again |
| 17 | one SSE2 primitive `simd::dots4` (four shifted-window dots, one shared reduction) for `comp_corr` (its group loop back in plain Rust) and `autocorr`'s 10 lags (three groups over its padded window); 33 lines in, 49 out of `simd.rs`. Oracle `dots4_is_the_array_form` (all −32768 too; a window-offset mutant fails) | 1.008 0.990 1.009 0.996 1.026 1.013 0.998 1.020 (all 1.010) | – | yes | simpler (one kernel, two users) at equal or better speed |
| 18 | `SynFilt`: `syn_filt`'s coefficient tables made once per coefficient set and reused (decoder postfilter `ap4` ×2, encoder `subframe_pre_proc` `aq` ×2, `ap2` ×2) | 0.993 0.995 0.990 0.996 0.996 0.991 0.995 0.970 (all 0.989) | 0.979 0.977 0.983 0.976 0.968 0.988 0.965 0.977 (all 0.976) | no | slower: built inline the tables never leave registers; kept in a struct they are stored and reloaded per block |
| 19 | `l_add_energy` (`bgn_scd`'s frame energy; its line was callgrind's heaviest in the decoder, 2.6%) with squares summed in u32 pairs, widened once per pair | 1.002 1.019 1.028 1.016 1.016 0.992 0.996 0.998 (all 1.006) | 1.000 1.016 1.006 1.000 0.995 1.004 1.003 1.003 (all 1.004) | no | flat: 160 scalar i64 adds per frame are cheap in time |
| 20 | step-5 searches' saturation guard deleted: with the C's shifts and any 16-bit `rr`, every energy chain stays in [−2^31, 0x7fff7fff] (pulse 1 ≤ 2,147,418,112; pulse 2 |·| ≤ 1,207,959,552 so |alp| ≤ 18432; pulse 3 ≤ 2,147,389,440, ≥ −2^31 exactly), so neither `L_mac` nor `round` can saturate: the per-row bound in `lanes_next`, the \|term\| sums for it, the `round` clamp and the scalar fallback `step5_exact` were dead. The worst cases are a `const` assertion over the shift table both searches read (`STEP5_SHIFTS`; changing one shift fails to compile; P <= 4 asserted). `lanes_next`'s test now draws from its callers' domain; `step5_block`'s C-loop test gains a ±full-scale class (a first tier from block maxima, tried first, left a guard-always-true mutant alive: why the proof was found) | 1.012 0.999 1.035 1.050 1.090 1.070 1.005 0.998 (all 1.035) | – | yes | proving a guard unnecessary beat making it cheaper; −67 lines |
| – | (tooling) `search_10i40`/`cor_h_search_8i40` kernel timings: their input's full-scale `dn` sent every call to the exact search (2 Σ max\|dn\| > 32767), which the codec never takes (1 call in a whole profiled run); now `dn` as `cor_h_x` normalises it: 4962 → 1034 ns (C 2577), 5062 → 2053 ns (C 2681) | – | – | – | the kernels looked slower than the C |
| 21 | MR102's `lag_max_wgh`: the lags' weighting in lanes (`Mpy_32_16` of an `L_Extract` pair by a positive weight can't saturate: 2 hi w + 2 ((lo w) >> 15) in i32), then max + first position; its two serial `L_mac` chains as `l_mac_dot`/`l_mac_energy`. A tie mutant (`rposition`) fails the vectors and the whole-encoder reference tests | 0.999 0.999 0.996 0.990 1.003 1.009 **1.030** 1.004 (all 1.002) | – | yes | 10.2's weighted open-loop pick was 5.2% of its instructions, scalar |
| 22 | MR102's AVX2 copy with `cor_h_tm`'s multiply-free sign scaling (i16 lanes, `NARROW`; its 16-lane rows fit one `ymm` as i16, two as i32) | kernel 1.06 vs the C's 1.04 drift (realistic input) | – | no | the `vpmulld` form stays faster, as for 12.2's 8-lane rows |
| 23 | decoder `bits2prm`: bits packed into big-endian u64 words (when all are 0/1), each parameter one shifted window | – | 0.968 0.988 0.970 0.992 0.981 0.986 0.991 0.981 (all 0.978) | no | slower: the packing fold is a serial shift chain per bit |
| – | (tooling) `convolve`/`pred_lt_3or6` kernel timings on encoder-like inputs (excitation range, decaying impulse response): full-scale inputs took the exact loops every call. 424 → 98 ns (C 267), 265 → 82 ns (C 152); `syn_filt` (127 vs 172 ns, 1.35×) is now the weakest kernel against the C | – | – | – | kernel timings must take the codec's path |
| 24 | `syn_filt`'s bound from max \|x\|, \|y\| collected inside `syn_blocks` (`pminsw`/`pmaxsw` on the vectors it holds, `mem` from the first history load) instead of three `max_abs` scans after it. Ablation first: the bound check costs 28 of 127 ns (the coefficient tier fails for real LPC filters), the tables 20 ns | – | – (kernel 1.00) | no | the in-kernel min/max cost what the scans did; the ablation's 28 ns was not the scans |
| 25 | `syn_filt`'s tables built in SSE2 (each `colp`/`upd` row: windows of a zero-padded coefficient array, `punpcklwd`/`punpckhwd`) | – | – (kernel 1.07, slower) | no | stored, then reloaded per block; built inline, LLVM keeps part of it in registers (as exp 18) |
| – | (tooling) `syn_filt` kernel timing inputs filtered by its bound: 1.6% of them took the exact loop (the decoder: 39 in 714k calls), 127 → 107 ns (C 168). Exp 24's ablation (28 ns for the bound check) was mostly those fallbacks | – | – | – | – |
| 26 | `az_lsp`'s bisection two levels at a time (the first midpoint and both possible second ones evaluated concurrently: 2 `chebps` latencies per 2 steps instead of 2, 6 calls instead of 4) | – (kernel 0.97) | – | no | below the bar: the grid scan already overlaps under speculation, and the bisection less than its 4 chains suggested |
| – | (tooling) `BENCH=stages just bench`: the encoder's µs per frame per stage (`--cfg amr_nb_stages` timers, no-ops otherwise): codebook search 10.5–13.9 of 25–30 µs in 5.90–12.2 (6.70: 47%), callgrind had it near 20%; 4.75's gains 4 µs; lpc+lsp 3.3 (12.2: 5.7); closed-loop pitch 3.7–4.2 | – | – | – | time attribution without `perf` |
| 27 | step-5 searches: all blocks (or half-block passes) made first, then the C's compare chain over their lanes, in its order (the blocks are independent; the chain interleaved with them serialised their latency-bound lane chains). Kernels: `search_3i40` 0.90, `search_4i40` 0.96 | 1.020 1.010 **1.047** 1.005 1.016 1.013 0.992 1.000 (all 1.010) | – | yes | latency-bound blocks overlap once nothing between them depends on the chain |
| 28 | `lanes_next` over two blocks' streams in one row loop (`step5_blocks::<P, 2>`, `search_3i40` as 6 pairs): independent chains interleaved, 16 lanes per AVX2 register | – (kernel `search_3i40` 1.33, slower) | – | no | two streams' 8×8 i32 cross terms and lane state outgrow the registers |
| 29 | `MR475_gain_quant`'s lanes with a first tier: every term is non-negatively weighted, so all partial sums are bounded by Σ \|c\| n_max over the table's column maxima (a `const`); within MAX_32 the lanes skip their per-add overflow flags (4 ops on each of 20 adds per row). Kernel 0.77; an always-pass mutant fails the MR475 reference test | **1.028** 0.985 0.966 0.986 1.001 0.996 0.999 1.002 (all 0.996) | – | yes | 4.75's gain stage is 19% of its frame (`BENCH=stages`); the other modes don't run it (their ±3%: layout) |
| 30 | `Qua_gain`'s lanes (MR515–MR102) with the same first tier: the columns' maxima bound every n (≤ 0x7fff: no `mult` saturates, no `Mpy_32_16` product is 2^30) and every partial sum; within MAX_32 the per-operation flags go. Kernel 0.75; an always-pass mutant fails the `Qua_gain` reference test | 1.000 **1.010 1.015 1.023 1.039** 0.999 **1.012** 1.006 (all 1.012) | – | yes | every mode that runs it gains, the others are flat |
| 31 | `MR795_gain_code_quant3`'s lanes with the same first tier (|acc0| + the terms at the table's g maximum) | 1.010 0.996 1.000 1.003 0.993 **1.010** 0.997 0.995 (all 0.999) | – | no | kernel 0.85, but 7.95 +1.0%: within noise (4 calls per frame, ~0.4 µs each) |
| 32 | `az_lsp`'s grid values for both polynomials precomputed in i32 lanes (b x >> 16 as (b >> 16) x + ((b & 0xffff) x >> 16), exact under `chebps`' guard), with an AVX2 copy; the scan reads them. Stage split first (12.2: `az_lsp` ×2 2.3 µs, `q_plsf_5` 1.8, `lpc` 1.5 per frame) | 0.993 0.999 0.999 0.986 1.007 1.001 1.000 0.998 (all 0.997) | – | no | kernel 0.83 but flat end to end (12.2 too): the scalar scan already overlaps under speculation, and the lanes evaluate twice the points |
| – | (tooling) `BENCH=stages` covers the decoder too: per frame the postfilter is ~1.7 of ~5 µs (a third), parameters + pitch + code + gains ~1.05, synthesis ~0.9, `post_process` 0.55 | – | – | – | – |
| 33 | postfilter's 22-sample impulse response (`syn_filt` of ap3 through 1/ap4) filtered over 24 samples, x zero-padded (causal: the first 22 are the C's): whole blocks, no tail block. Stage split first: per 7.40 frame weight_ai + residu 0.6 µs, this filter 0.46, autocorr_h + preemphasis 0.34, the main syn_filt 0.5, agc 0.14 | – | 1.062 1.017 1.016 1.031 1.036 1.010 1.031 1.028 (all 1.027) | yes | the tail block's history rebuild and copies cost more than two extra outputs |
| 34 | `weight_ai` in lanes for factors >= 0 (every gamma table: neither L_mult nor round can saturate) | 1.018 1.000 1.007 0.998 1.009 1.007 0.996 0.998 (all 1.004) | 1.003 1.005 1.012 1.003 0.995 0.985 1.014 1.004 (all 1.002) | no | flat: 10 scalar products per call |
| 35 | `cor_h_tm`'s AVX2 copy with the multiply-free ±32767 sign scaling in i32 lanes (t − ceil(t/16384) or (t >> 15) − t, selected; the i16 form had lost to i16 repacking). Stage split before: 12.2's `q_plsf_5` VQs 1.34 µs per frame, likely table misses in a full frame (lane throughput ~0.3–0.6 µs) | – (kernel 1.21, slower) | – | no | two `vpmulld` beat the shifts and selects on AVX2 in i32 lanes too |
| 36 | `lanes_next`: when every row's energy is >= 1, a tournament over the 8 rows (cross-multiplied ratio compares are transitive then; the earlier row keeps ties) instead of the C's chain of 8 selections (depth 3 instead of 8) | – (kernels `search_3i40` 2.16, `search_4i40` 1.79: slower) | – | no | 8 rows × 4 lane vectors live at once spill; the chain keeps 4 |
| 37 | `lanes_next` as two interleaved chains (even rows, odd rows), merged by ratio (lower row on ties) when every energy is >= 1: half the chain latency, little extra state | – (kernels `search_3i40` 1.04, `search_4i40` 1.06: slower) | – | no | with exp 36 and 28: the step-5 searches are not limited by the row chain's latency |
| 38 | portable build: `Vq_subvec3/4` in SSE2 i16 lanes (`simd::vq_subvec_sse2`: t = (`pmulhw` << 1) \| (`pmullw` >>> 15), each codeword's Σ t² as `pmaddwd` over its components interleaved pairwise); the AVX2 copy keeps its i32 form. Stage split first: portable vs AVX2 per frame, lpc+lsp +1.9 µs (the VQs 3.7× slower without `pmulld`), codebook +4.5–7.4, gains +1.7. Kernels 0.68 / 0.53 (portable); oracle in `vq_subvec_matches_c` on every path; a dropped-low-bit mutant fails | portable: 1.009 1.009 0.998 **1.020 1.028 1.030 1.033** 0.993 (all 1.017) | – | yes | SSE2 has no `pmulld`: i16 lanes with `pmulhw`/`pmullw` halves |
| 39 | MR67's `search_3i40` with its valid first pulses packed: every lane its own track order (`step5_lanes`: `d`/`diag` and cross rows gathered per lane, transposed), 72 candidates in 9 passes instead of 12 blocks | – (kernel `search_3i40` 1.37, slower) | – | no | the per-lane gathers (2 extra transposes per pulse, row gathers by lane) cost more than 3 passes |
| 40 | `comp_corr`'s AVX2 copy with 256-bit `vpmaddwd` (a `#[target_feature]` body called from the dispatch: no `unsafe` call site) | 1.017 1.027 1.010 0.958 0.995 1.006 0.989 1.016 (all 0.995; 6.70 re-measured alone: flat) | – | no | kernel 0.53, but ~0.3 µs of a frame: within noise |
| 41 | decoder: each subframe post-filtered right after its synthesis (`post_filter_sf`, separate output buffer, the post-filter's overflow flag OR'd in at frame end), so its chains sit beside the next subframe's synthesis instead of after the whole frame. Bit-exact on every decoder vector path | – | 1.012 0.993 1.004 1.047 0.985 1.003 0.985 0.992 (all 1.004) | no | flat: adjacency didn't let the core overlap the two chains (as exp 11) |

## Fuzzing (`just fuzz [target] [seconds]`)
cargo-fuzz (libFuzzer) on nightly-2026-11-01, debug assertions and overflow checks on, ASan on
the Rust side. Corpus kept across runs in `.cache/fuzz/`. `rfc4867`: 2026-11-12, 600 s, corpus
continued. The codec targets diff against the TS 26.073 reference since 2026-11-14 (perf/loop;
they were against opencore, which wraps where the reference saturates): 600 s each, new corpora.
Their first run found a harness gap, not a codec bug: the reference leaves `decoder.c`'s
`prev_mode` (used for a NO_DATA frame) uninitialised; the Rust decoder starts from MR475.

| target | checks | runs | cov | findings |
|---|---|---:|---:|---|
| `rfc4867` | `.amr`, RTP octet-aligned and bandwidth-efficient readers on arbitrary bytes: no panics, sane error order, accepted frames round-trip through both packings, canonical form never longer | 2,023,256 | 524 | none |
| `decode_vs_ref` | Rust decoder vs the TS 26.073 reference `Speech_Decode_Frame` (`ref26073-sys`), PCM equal frame by frame: every RX frame type (bad, degraded, onset, SID first/update/bad, NO_DATA), any mode and bits; stops at a homing frame (the reference homes in `decoder.c`) | 452,074 (600 s); 30 min: 1; after the fix 1,175,462 (1 h): 0 | 2,162 | defect 12 (DTX mute clamp, inherited from opencore), fixed; both saved crash inputs replay clean |
| `encode_vs_ref` | Rust encoder (DTX off / VAD1, per-frame mode from the input) vs the reference `Speech_Encode_Frame` + `sid_sync` on arbitrary PCM: TX type and serial bits equal | 550,312; 1,015,759 (1 h) | 5,518 | none |

Earlier, against opencore (2026-11-12, 600 s): decoder 100,831 runs / cov 2,055, encoder
272,121 / 5,148, no findings; 2026-11-08 (pre-rewrite code): 1.84 M / 204 k / 1.22 M runs.

Known C defects (list below) kept out of the differential, each by construction:
- **1 no homing**: comparison stops at the first frame that could home the Rust decoder (the
  encoder target: after the first encoder homing frame).
- **2 `const` input overwritten**: the C encoder gets a copy of each frame.
- **3 `A_Refl` uninitialised read, 4 `L_add` UB**: oracle built with `opencore-sys/pin-ub`
  (`-ftrivial-auto-var-init=zero -fwrapv`). A 180 s run without it found no difference either.
- **5 frame types 9–14**: skipped (12–14 end the input: size unknown).
- **6 no length check**: C gets every frame zero-padded to 64 octets.
- **7 quality bit ignored**: compared frames have Q forced to 1. Sanity check: dropping this
  exclusion makes the target fail within seconds (`SidBad`, frame 0), so the oracle bites.

## Release 0.1 (2026-11-14)
`amr-nb` 0.1.0, exported from the development repository to `amr-rs` as one commit. In the
exported tree (78 files), before the commit:
- `just check` green: export guard, fmt, clippy (x86-64, portable, aarch64, fuzz targets), docs
  with `-D warnings`, every test suite against the oracles (131 unit tests and 8 conformance
  tests per build; doc tests 6/6), 248 s including `cargo package`.
- `cargo package -p amr-nb` clean: 44 files, 241 kB (`src/`, `examples/`, README, LICENSE,
  NOTICE; no tests or benches, which need the C oracles; path dev-dependencies stripped).
- The export guard (`scripts/export-guard.sh`) fails on each planted violation: development-only
  files and notes, a vector file, a `.cpp`, a C file that isn't one of the FFI shims, a
  reference to an unexported doc, a third-party copyright line outside NOTICE, a `just` recipe
  that doesn't exist. Its first real catch: this section, which named development-only files
  (the release was refused and the clone restored). The development repository's name is
  checked by the export script, which doesn't ship.

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
    post-filter's energy chain. The Rust port inherited them (it was differentially tested
    against opencore) until 2026-11-12; now the reference is the oracle.
11. **The decoder's overflow re-synthesis is unreachable**: opencore's `Syn_filt` reports no
    overflow, so `Decoder_amr` never takes the reference's "scale the excitation by 1/4 and
    re-synthesise" branch. Rust: `syn_filt` returns the reference's flag.
12. **DTX mute clamps the comfort-noise energy at 0**: entering DTX_MUTE (no SID update for
    50 frames), `dtx_dec` lowers `log_en` by 256 and opencore clamps it at 0 where TS 26.073
    saturates and lets it go negative; a low-energy SID followed by long silence then decodes
    to different comfort noise. The Rust port had inherited the clamp; `decode_vs_ref` found
    it after 30 min (2026-11-14, perf/loop) and it is fixed (test
    `decoder_matches_ref_through_dtx_mute`; the vectors never reach it).

Also differs from the current reference, not a bug: `dtx_enc` recomputes the SID parameters
while `log_en_index == 0` (an older reference version; the vectors pass either way, the Rust
port follows the current one). And opencore's unrolled loops drop the remainder (`Syn_filt` an
odd last sample, `Residu` the first `lg % 4`, `Pred_lt_3or6` an odd last one; the codec only
uses lengths 22 and 40).

The Rust codec's only `unsafe` is `src/simd.rs` (runtime AVX2 dispatch and SSE2 transposes),
plus FFI in the tests.
