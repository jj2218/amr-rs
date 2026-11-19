# Every recurring task is a recipe here, backed by a script when it has logic: `just --list`.
# The recipes below `rules` are the development repository's, unchanged.

default:
    @just --list

# The full gate: export guard, format, lint, tests against the oracles (downloaded into .cache/)
check:
    @bash scripts/check.sh

# Export guard: no third-party source, vectors or stray files in the repo
rules:
    @bash scripts/export-guard.sh

fmt-check:
    cargo fmt --all --check
    cargo fmt --manifest-path fuzz/Cargo.toml --check

lint: opencore vectors msrv
    cargo clippy --workspace --all-targets --locked -- -D warnings
    RUSTFLAGS="--cfg portable" CARGO_TARGET_DIR=target/portable cargo clippy -p amr-nb -p amr-wb --all-targets --locked -- -D warnings
    cargo clippy --manifest-path fuzz/Cargo.toml --locked -- -D warnings
    cargo clippy -q -p amr-nb -p amr-wb --lib --locked --target aarch64-unknown-linux-gnu -- -D warnings
    RUSTDOCFLAGS="-D warnings" cargo doc -q -p amr-nb -p amr-wb --no-deps --locked

# The library builds with its minimum supported Rust (`rust-version` in Cargo.toml; installs it)
msrv:
    #!/usr/bin/env bash
    set -euo pipefail
    v=$(sed -n 's/^rust-version = "\(.*\)"$/\1/p' Cargo.toml)
    rustup toolchain install "$v" --profile minimal --no-self-update >/dev/null 2>&1
    CARGO_TARGET_DIR=target/msrv cargo "+$v" check -q -p amr-nb -p amr-wb --lib --locked
    echo "msrv: amr-nb and amr-wb build with Rust $v"

# Needs the C oracle (opencore-sys); fetched automatically, cached in .cache/. Runs the codec
# tests three ways in parallel: runtime AVX2 dispatch (default), the portable path
# (src/simd.rs), the vectors with the SSE2 kernels' array forms (the non-x86 bodies; their unit
# oracles run in the other two). One line per run, a failing run's whole log (scripts/test.sh).
test: opencore vectors
    @bash scripts/test.sh

# The tests and conformance on aarch64: cross-built, each test binary under qemu-aarch64 (the portable code; correctness only, qemu timing means nothing; C-oracle-only tests skipped; ~15 min)
test-aarch64: opencore vectors
    @TARGET=aarch64 bash scripts/test.sh

# Download + SHA-verify the 3GPP test vectors (TS 26.074, TS 26.174) and reference C (TS 26.073, TS 26.173) into .cache/3gpp (never committed)
vectors:
    @bash scripts/vectors.sh

# Codec conformance. AMR-NB vs the whole 3GPP TS 26.074 set (linear, A-law, µ-law; bit + frame sync): decoder (direct and via RFC 4867), encoder; pinned vector deviations are counted, not failed. AMR-WB decoder and encoder vs the TS 26.174 set (plus the C reference, opencore and vo-amrwbenc on it), and vs the TS 26.173 reference on stress streams (decoder 1600, encoder 64)
conform: opencore vectors
    @cargo test -q --release -p amr-nb --test conformance -- --nocapture 2>&1 | grep -E "bit-exact|failing|mismatch|panicked|error" | sed "s/^\.*//" || true
    @cargo test -q --release -p amr-wb --test conformance -- --nocapture 2>&1 | grep -E "bit-exact|panicked|error" | sed "s/^\.*//" || true
    @cargo test -q --release -p amr-wb --test enc_conformance -- --nocapture 2>&1 | grep -E "bit-exact|bit sync|panicked|error" | sed "s/^\.*//" || true
    @WB_STRESS_STREAMS=200 cargo test -q --release -p amr-wb --test vs_ref -- --nocapture 2>&1 | grep -E "bit-exact|differs|panicked|error" | sed "s/^\.*//" || true
    @cargo test -q --release -p amr-wb --test enc_vs_ref -- --nocapture 2>&1 | grep -E "bit-exact|differs|panicked|error" | sed "s/^\.*//" || true

# Decoder + encoder speed, Rust vs opencore C (-O3), pinned to one core; µs/frame per mode (BENCH=dec|enc; BENCH=all: both codecs × sides, the "all" lines (before/after a sweep); BENCH=stages: encoder and decoder µs/frame per stage, time attribution; BENCH=kern: ns/call per kernel, KERN=<substring> filters; AB=<ref>: ref vs working tree (B=<ref>: vs that ref; 3 reps per run), BENCH_MODE=<mode>: one mode only (each: per-mode table); with kern: min of 3 + ratio)
bench: opencore vectors
    @bash scripts/bench.sh

# Disassemble one function (count, top instructions, listing in target/asm.s), e.g. 'cor_h_tm.*avx2_copy'
[positional-arguments]
asm symbol:
    @bash scripts/asm.sh "$1"

# Rank the Rust codec's functions by sampled cycles (perf on the shipped build, inline frames resolved; AVX2 dispatch and portable; CRATE=amr-wb, BENCH=dec|enc|all, BENCH_MODE=<mode>, TOP=<n>)
perf: opencore vectors
    @bash scripts/perf.sh

# Rank the Rust codec's functions by self cost (callgrind, no-inline copy; CRATE=amr-wb, BENCH=dec|enc, BENCH_MODE=<mode>: one mode)
profile: opencore vectors
    @bash scripts/profile.sh

# Fuzz (nightly + cargo-fuzz, installed on first use): per codec (nb_*, wb_*) the rfc4867 readers, the decoder and the encoder vs the 3GPP reference C
[positional-arguments]
fuzz target="all" seconds="60": vectors
    @bash scripts/fuzz.sh "$1" "$2"

# Fetch + SHA-verify opencore-amr 0.1.6 and vo-amrwbenc 0.1.3 into .cache/opencore (C oracles + bench baselines)
opencore:
    @bash scripts/opencore.sh

# Regenerate the codec crates' src/tables.rs from the opencore C tables (amr-wb's enc_tables.rs: vo-amrwbenc's)
gen-tables: opencore
    @python3 -I scripts/gen_tables.py nb .cache/opencore/opencore-amr-0.1.6/opencore/codecs_v2/audio/gsm_amr/amr_nb crates/amr-nb/src/tables.rs
    @python3 -I scripts/gen_tables.py wb .cache/opencore/opencore-amr-0.1.6/opencore/codecs_v2/audio/gsm_amr/amr_wb crates/amr-wb/src/tables.rs
    @python3 -I scripts/gen_tables.py wbenc .cache/opencore/vo-amrwbenc-0.1.3/amrwbenc crates/amr-wb/src/enc_tables.rs crates/amr-wb/src/tables.rs
