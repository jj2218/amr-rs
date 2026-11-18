# amr-rs

Pure-Rust AMR speech codecs, bit-exact against the 3GPP specifications.

```sh
cargo add amr-nb   # AMR-NB, 8 kHz
cargo add amr-wb   # AMR-WB, 16 kHz
```

| crate | |
|---|---|
| [`amr-nb`](crates/amr-nb) | AMR-NB (3GPP TS 26.073): encoder, decoder, DTX, RFC 4867 packing (`.amr`, RTP). Bit-exact on the whole TS 26.074 vector set; faster than opencore-amr's C. [README: example, speed, scope](crates/amr-nb/README.md) |
| [`amr-wb`](crates/amr-wb) | AMR-WB (3GPP TS 26.173, ITU-T G.722.2): decoder, encoder, DTX, RFC 4867 packing (`.awb`, RTP). Bit-exact on the whole TS 26.174 vector set and against the TS 26.173 reference; faster than opencore-amr's decoder, and with AVX2 than vo-amrwbenc. [README: example, speed, scope](crates/amr-wb/README.md) |
| `opencore-sys`, `ref26073-sys`, `ref26173-sys` | test-only (not published): the C oracles, compiled from pinned downloads |

No dependencies. On x86-64 CPUs with AVX2 the hot kernels use it, detected
at run time; elsewhere the portable code runs, with the same output. Measured results, each
with the command that reproduces it: [`docs/RESULTS.md`](docs/RESULTS.md).

## Development
Needs a stable Rust toolchain (pinned in `rust-toolchain.toml`, installed by rustup), [`just`](https://github.com/casey/just),
a C compiler, `nm`/`objcopy` (binutils), `curl` and `unzip`.

```sh
just check     # the full gate: guard, fmt, clippy, docs, tests against the oracles
just conform   # bit-exactness vs the 3GPP vectors (TS 26.074, TS 26.174), one summary line each
just bench     # µs per frame, Rust vs opencore's C
just --list    # everything else (fuzz, profile, ...)
```

**Oracles.** The tests compare against the 3GPP test vectors (TS 26.074, TS 26.174), the 3GPP
reference C (TS 26.073, TS 26.173) and opencore-amr 0.1.6. `just check` downloads them on first use into the
gitignored `.cache/` and verifies each against a pinned SHA-256. The 3GPP material is "all
rights reserved": it is used as a test oracle only, never committed, never ported from. The
codecs are ported from opencore-amr (Apache-2.0).

**CI** can't fetch the 3GPP files (their server refuses GitHub's runners), so it checks what
builds without them; run `just check` locally before a pull request.

**Releases** are published from a private development repository, one commit per release.
Pull requests are welcome: an accepted change is applied there with the same paths, credited
with `Co-authored-by`, and ships with the next release; the pull request is then closed with a
link to it.

## License
Apache-2.0 ([LICENSE](LICENSE)); attribution in [`crates/amr-nb/NOTICE`](crates/amr-nb/NOTICE)
and [`crates/amr-wb/NOTICE`](crates/amr-wb/NOTICE).
