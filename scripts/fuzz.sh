#!/usr/bin/env bash
# Coverage-guided fuzzing with cargo-fuzz (libFuzzer) on a pinned nightly; the only nightly
# use in the repo (`just check` and CI stay on stable). Targets: fuzz/fuzz_targets/*.rs.
# Usage: fuzz.sh [target|all] [seconds per target]. Corpus and crash inputs persist in the
# gitignored .cache/fuzz/{corpus,artifacts}/<target>/. Prints one line per target.
set -euo pipefail
cd "$(dirname "$0")/.."

NIGHTLY=nightly-2026-11-01
CARGO_FUZZ=0.13.2
want=${1:-all}
secs=${2:-60}

rustup toolchain list | grep -q "^$NIGHTLY" || rustup toolchain install "$NIGHTLY" --profile minimal
if [ "$(cargo fuzz --version 2>/dev/null)" != "cargo-fuzz $CARGO_FUZZ" ]; then
  cargo install cargo-fuzz --locked --version "$CARGO_FUZZ"
fi

if [ "$want" = all ]; then
  targets=$(basename -s .rs fuzz/fuzz_targets/*.rs)
else
  targets=$want
fi
cargo +"$NIGHTLY" fuzz build --fuzz-dir fuzz 2>&1 | grep -E "^(error|warning)" || true

status=0
for t in $targets; do
  corpus=.cache/fuzz/corpus/$t art=.cache/fuzz/artifacts/$t log=.cache/fuzz/$t.log
  mkdir -p "$corpus" "$art"
  if cargo +"$NIGHTLY" fuzz run --fuzz-dir fuzz "$t" "$corpus" -- \
      -max_total_time="$secs" -artifact_prefix="$art/" -print_final_stats=1 >"$log" 2>&1; then
    runs=$(grep -oP 'stat::number_of_executed_units:\s*\K\d+' "$log" || echo "?")
    cov=$(grep -oP '#\d+\s+DONE\s+cov: \K\d+' "$log" || grep -oP 'cov: \K\d+' "$log" | tail -1)
    echo "fuzz $t: ok, $runs runs in ${secs}s, cov $cov, corpus $(ls "$corpus" | wc -l)"
  else
    status=1
    echo "fuzz $t: FAILED (log $log)"
    grep -m1 -E "panicked at|ERROR: libFuzzer|SUMMARY" "$log" || true
    grep -A1 -m1 "panicked at" "$log" | tail -1 || true
    grep -oE "$art/[a-z]+-[0-9a-f]+" "$log" | tail -1 | sed 's/^/  reproduce: cargo +'"$NIGHTLY"' fuzz run --fuzz-dir fuzz '"$t"' /' || true
  fi
done
exit $status
