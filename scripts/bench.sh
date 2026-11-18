#!/usr/bin/env bash
# CRATE=<crate> (default amr-nb): which codec. Build its benchmark once and run it pinned to one core (less noise on laptops).
# BENCH=kern: ns per call of the hot encoder kernels, Rust vs C (an ignored lib test).
# AB=<git ref>: that ref (A) vs the working tree (B; B=<git ref>: that ref instead, e.g. two
# `git stash create` snapshots of consecutive steps), interleaved runs of the "all" line
# (BENCH=enc by default; BENCH_MODE=<mode> runs and compares only that mode, e.g. 122); the C drifts
# with temperature, so judge by the Rust/C ratio. BENCH_MODE=each: every mode in one run, then
# the median ratio per mode for A and B (B/A > 1: B faster).
# AB=<ref> BENCH=kern: per kernel, min ns over 3 interleaved runs of each side, and B/A.
set -euo pipefail
cd "$(dirname "$0")/.."
crate=${CRATE:-amr-nb}  # CRATE=amr-wb: the other codec
# BENCH=all: every codec × side, its "all" line only: a sweep's before/after (run on each tree).
if [ "${BENCH:-}" = all ]; then
    [ -z "${AB:-}" ] || { echo "bench: BENCH=all doesn't take AB= (run it on each tree)" >&2; exit 1; }
    for c in amr-nb amr-wb; do for b in dec enc; do
        printf '%-6s %s %s\n' "$c" "$b" "$(env -u BENCH_MODE CRATE=$c BENCH=$b bash "$0" | grep '^all ' || echo failed)"
    done; done
    exit
fi
# BENCH=stages: the encoder's and decoder's µs per frame per stage (`--cfg codec_stages` timers,
# mean of all reps, timer overhead included): time attribution, as `perf` is unavailable.
if [ "${BENCH:-}" = stages ]; then
    export RUSTFLAGS="${RUSTFLAGS:-} --cfg codec_stages" CARGO_TARGET_DIR="$PWD/target/stages"
    unset BENCH
    exe=$(cargo bench -p "$crate" --bench codec --no-run 2>&1 | sed -n 's/.*Executable .*(\(.*\)).*/\1/p' | tail -1)
    taskset -c "${BENCH_CPU:-2}" "$exe" | grep "^stages"
    exit
fi
# BENCH_MODE=all is the default (every mode, the "all" line).
[ "${BENCH_MODE:-}" = all ] && unset BENCH_MODE
# The portable build gets its own target dirs (as in `just test`), so switching paths doesn't
# rebuild the fat-LTO binaries each time.
ab=target/ab
case "${RUSTFLAGS:-}" in *portable*)
    export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$PWD/target/portable}"; ab=target/ab-portable ;; esac
build() { cargo bench -p "$crate" --bench codec --no-run 2>&1 | sed -n 's/.*Executable .*(\(.*\)).*/\1/p' | tail -1; }
kern() { taskset -c "${BENCH_CPU:-2}" cargo test -q --release -p "$crate" --lib kernel_timing \
    -- --ignored --nocapture 2>&1 | grep -E " ns$|error|panicked"; }
if [ -n "${AB:-}" ]; then
    # best of 3 reps per run: the A/B takes the median of 3 interleaved runs (15: ~50 min)
    export BENCH_REPS=${BENCH_REPS:-3}
    src=target/ab-src
    git worktree remove --force "$src" 2>/dev/null || rm -rf "$src"
    git worktree add -q --detach "$src" "$AB"
    trap 'git worktree remove --force "$src"' EXIT
    ln -s "$PWD/.cache" "$src/.cache"
    if [ -n "${B:-}" ]; then
        bsrc=target/ab-src-b
        git worktree remove --force "$bsrc" 2>/dev/null || rm -rf "$bsrc"
        git worktree add -q --detach "$bsrc" "$B"
        trap 'git worktree remove --force "$src"; git worktree remove --force "$bsrc"' EXIT
        ln -s "$PWD/.cache" "$bsrc/.cache"
    fi
    if [ "${BENCH:-}" = kern ]; then
        [ -z "${B:-}" ] || { echo "bench: B=<ref> is not supported with BENCH=kern (side B is the working tree)" >&2; exit 1; }
        for r in 1 2 3; do
            a=$(cd "$src" && CARGO_TARGET_DIR="$OLDPWD/$ab" kern) \
                || { echo "bench: $AB has no kernel_timing test (or it failed)" >&2; exit 1; }
            b=$(kern) || { echo "bench: kernel_timing failed in the working tree" >&2; exit 1; }
            sed 's/^/A /' <<<"$a"; sed 's/^/B /' <<<"$b"
        done | awk -v kf="${KERN:-}" '/error|panicked/ {print; bad=1; next}
            { side=$1; ns=$(NF-1); name=$0; sub(/^[AB] +/, "", name); sub(/ +[0-9.]+ ns$/, "", name)
              if (kf != "" && index(name, kf) == 0) next
              k=side SUBSEP name; if (!(k in m) || ns < m[k]) m[k]=ns
              if (!(name in seen)) { seen[name]=1; order[++n]=name } }
            END { if (bad) exit 1
              printf "%-22s %9s %9s  %5s\n", "kernel (min of 3)", "A ns", "B ns", "B/A"
              for (i=1; i<=n; i++) { a=m["A" SUBSEP order[i]]; b=m["B" SUBSEP order[i]]
                printf "%-22s %9.1f %9s  %5s\n", order[i], a, (b=="" ? "-" : sprintf("%.1f", b)),
                    (b=="" || a==0 ? "-" : sprintf("%.2f", b/a)) } }'
        exit
    fi
    abs() { case "$1" in /*) echo "$1" ;; *) echo "$2/$1" ;; esac; }
    if [ -n "${B:-}" ]; then
        b=$(abs "$(cd "$bsrc" && CARGO_TARGET_DIR="$OLDPWD/$ab-b" build)" "$PWD/$bsrc")
    else
        b=$(abs "$(build)" "$PWD")
    fi
    a=$(abs "$(cd "$src" && CARGO_TARGET_DIR="$OLDPWD/$ab" build)" "$PWD/$src")
    [ -f "$a" ] && [ -x "$a" ] && [ -f "$b" ] || { echo "bench: AB build failed (does ${B:+$B or }$AB have the codec bench?)" >&2; exit 1; }
    if [ "${BENCH_MODE:-}" = each ]; then
        for r in 1 2 3; do for v in A B; do
            exe=$a; [ $v = B ] && exe=$b
            BENCH=${BENCH:-enc} env -u BENCH_MODE taskset -c "${BENCH_CPU:-2}" "$exe" |
                awk -v v=$v '$4 ~ /x$/ {print v, $1, $2/$3}'
        done; done | awk '{ k=$1 SUBSEP $2; n[k]++; r[k, n[k]]=$3; if (!($2 in seen)) { seen[$2]=1; o[++m]=$2 } }
            function med(k,  i, j, t, a) { for (i=1; i<=n[k]; i++) a[i]=r[k, i]
                for (i=1; i<=n[k]; i++) for (j=i+1; j<=n[k]; j++) if (a[j]<a[i]) { t=a[i]; a[i]=a[j]; a[j]=t }
                return a[int((n[k]+1)/2)] }
            END { printf "%-5s %7s %7s  %6s\n", "mode", "A x", "B x", "B/A"
                for (i=1; i<=m; i++) { ra=med("A" SUBSEP o[i]); rb=med("B" SUBSEP o[i])
                    printf "%-5s %6.2fx %6.2fx  %6.3f\n", o[i], 1/ra, 1/rb, ra/rb } }'
        exit
    fi
    for r in 1 2 3; do for v in A B; do
        exe=$a; [ $v = B ] && exe=$b
        BENCH=${BENCH:-enc} taskset -c "${BENCH_CPU:-2}" "$exe" | { grep -E "^${BENCH_MODE:-all} " || { echo "bench: no '${BENCH_MODE:-all}' line" >&2; exit 1; }; } |
            awk -v v=$v '{printf "%s rust %s  c %s  %s\n", v, $2, $3, $4}'
    done; done
    exit
fi
if [ "${BENCH:-}" = kern ]; then
    kern
    exit
fi
exe=$(build)
[ -x "$exe" ] || { echo "bench: could not locate the benchmark binary" >&2; exit 1; }
taskset -c "${BENCH_CPU:-2}" "$exe"
