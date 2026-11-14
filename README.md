# amr-rs

Pure-Rust AMR speech codecs, bit-exact against the 3GPP specifications.

| crate | |
|---|---|
| [`amr-nb`](crates/amr-nb) | AMR-NB (3GPP TS 26.073): encoder, decoder, DTX, RFC 4867 packing. Bit-exact on the whole TS 26.074 vector set; decoder ≈3.06×, encoder ≈2.25× the speed of opencore-amr's C. [README](crates/amr-nb/README.md) |
| `opencore-sys`, `ref26073-sys` | test-only (not published): the C oracles, compiled from pinned downloads |

Measured results, each with the command that reproduces it: [`docs/RESULTS.md`](docs/RESULTS.md).

## Development
Needs a stable Rust toolchain (pinned in `rust-toolchain.toml`, installed by rustup), [`just`](https://github.com/casey/just),
a C compiler, `nm`/`objcopy` (binutils), `curl` and `unzip`.

```sh
just check     # the full gate: guard, fmt, clippy, docs, tests against the oracles
just conform   # bit-exactness vs the 3GPP TS 26.074 vectors, one summary line each
just bench     # µs per frame, Rust vs opencore's C
just --list    # everything else (fuzz, profile, ...)
```

**Oracles.** The tests compare against the 3GPP test vectors (TS 26.074), the 3GPP reference C
(TS 26.073) and opencore-amr 0.1.6. `just check` downloads them on first use into the
gitignored `.cache/` and verifies each against a pinned SHA-256. The 3GPP material is "all
rights reserved": it is used as a test oracle only, never committed, never ported from. The
codec is ported from opencore-amr (Apache-2.0).

**CI** can't fetch the 3GPP files (their server refuses GitHub's runners), so it checks what
builds without them; run `just check` locally before a pull request.

**Releases** are published from a private development repository, one commit per release.
Pull requests are welcome: an accepted change is applied there with the same paths, credited
with `Co-authored-by`, and ships with the next release; the pull request is then closed with a
link to it.

## License
Apache-2.0 ([LICENSE](LICENSE)); attribution in [`crates/amr-nb/NOTICE`](crates/amr-nb/NOTICE).
