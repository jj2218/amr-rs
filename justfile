# Every recurring task is a recipe here, backed by a script when it has logic: `just --list`.

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
    RUSTFLAGS="--cfg amr_nb_portable" CARGO_TARGET_DIR=target/portable cargo clippy -p amr-nb --all-targets --locked -- -D warnings
    cargo clippy --manifest-path fuzz/Cargo.toml --locked -- -D warnings
    cargo clippy -q -p amr-nb --lib --locked --target aarch64-unknown-linux-gnu -- -D warnings
    RUSTDOCFLAGS="-D warnings" cargo doc -q -p amr-nb --no-deps --locked

# The library builds with its minimum supported Rust (`rust-version` in Cargo.toml; installs it)
msrv:
    #!/usr/bin/env bash
    set -euo pipefail
    v=$(sed -n 's/^rust-version = "\(.*\)"$/\1/p' Cargo.toml)
    rustup toolchain install "$v" --profile minimal --no-self-update >/dev/null 2>&1
    CARGO_TARGET_DIR=target/msrv cargo "+$v" check -q -p amr-nb --lib --locked
    echo "msrv: amr-nb builds with Rust $v"

# Codec tests vs the C oracles, twice: runtime AVX2 dispatch (default) and the portable path;
# then the vectors with the SSE2 kernels' array forms (the non-x86 bodies).
test: opencore vectors
    cargo test --workspace --locked
    RUSTFLAGS="--cfg amr_nb_portable" CARGO_TARGET_DIR=target/portable cargo test -q -p amr-nb --locked
    RUSTFLAGS="--cfg amr_nb_portable --cfg amr_nb_arrays" CARGO_TARGET_DIR=target/arrays cargo test -q --release -p amr-nb --locked --test conformance

# Download + SHA-verify the 3GPP TS 26.074 test vectors and TS 26.073 reference C into .cache/3gpp (never committed)
vectors:
    @bash scripts/vectors.sh

# Fetch + SHA-verify opencore-amr 0.1.6 into .cache/opencore (C oracle + speed baseline)
opencore:
    @bash scripts/opencore.sh

# Bit-exactness vs the whole 3GPP TS 26.074 set, decoder and encoder (one summary line each)
conform: opencore vectors
    @cargo test -q --release -p amr-nb --test conformance -- --nocapture 2>&1 | grep -E "bit-exact|failing|mismatch|panicked|error" | sed "s/^\.*//" || true

# Speed, Rust vs opencore C (-O3), µs/frame per mode, pinned core (BENCH_CPU=n, BENCH=dec|enc|stages|kern, AB=<ref>; see scripts/bench.sh)
bench: opencore vectors
    @bash scripts/bench.sh

# Rank the codec's functions by self cost (callgrind; BENCH=dec|enc, BENCH_MODE=<mode>)
profile: opencore vectors
    @bash scripts/profile.sh

# Disassemble one function of the release build, e.g. 'cor_h_tm.*avx2_copy'
[positional-arguments]
asm symbol:
    @bash scripts/asm.sh "$1"

# Fuzz (nightly + cargo-fuzz, installed on first use): rfc4867 readers, decoder and encoder vs the TS 26.073 reference
[positional-arguments]
fuzz target="all" seconds="60": vectors
    @bash scripts/fuzz.sh "$1" "$2"

# Regenerate crates/amr-nb/src/tables.rs from the opencore C tables
gen-tables: opencore
    @python3 -I scripts/gen_tables.py .cache/opencore/opencore-amr-0.1.6/opencore/codecs_v2/audio/gsm_amr/amr_nb crates/amr-nb/src/tables.rs
