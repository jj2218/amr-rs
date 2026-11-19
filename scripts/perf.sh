#!/usr/bin/env bash
# Rank the Rust codec's functions by sampled CPU cycles: time, where `just profile` counts
# instructions. perf records the shipped bench binary (fat LTO) pinned to BENCH_CPU (a
# performance core on hybrid CPUs); each sample's inline chain (addr2line -i, debug info) is
# charged to its innermost codec function outside basic_op, so code inlined into the AVX2
# copies and Decoder::decode still ranks as its own function (operators count to their caller).
# CRATE= (amr-nb), BENCH=dec|enc (enc; all: both codecs × sides), BENCH_MODE=<mode>, TOP= (10),
# PERF_REPS= (3). Two builds: runtime AVX2 dispatch and portable (--cfg portable), as in bench.
# The header gives the Rust codec's share of the samples (the C baseline runs in the same binary).
# Detail: perf report -i target/perf.data (the last run); perf annotate -i target/perf.data <sym>.
set -euo pipefail
cd "$(dirname "$0")/.."
if [ "${BENCH:-}" = all ]; then
    for c in amr-nb amr-wb; do for b in dec enc; do CRATE=$c BENCH=$b bash "$0"; done; done
    exit
fi
crate=${CRATE:-amr-nb} side=${BENCH:-enc}
case $side in dec | enc) ;; *) echo "perf: BENCH=dec|enc|all" >&2; exit 1 ;; esac
# Debian's `perf` wrapper refuses a host kernel it has no package for; perf_<ver> works then.
perf=
for p in perf $(compgen -c perf_ | sort -uV); do
    "$p" stat -e task-clock true >/dev/null 2>&1 && { perf=$p; break; }
done
[ -n "$perf" ] || { echo "perf: no working perf (container: --cap-add=PERFMON, linux-perf; host: perf_event_paranoid <= 2)" >&2; exit 1; }
ev=cycles:u; [ -d /sys/bus/event_source/devices/cpu_core ] && ev=cpu_core/cycles/u  # hybrid: P-core PMU
for build in avx2 portable; do
    flags=${RUSTFLAGS:-} dir=$PWD/target
    [ $build = portable ] && flags="$flags --cfg portable" dir=$PWD/target/portable  # as just test/bench
    exe=$(RUSTFLAGS=$flags CARGO_TARGET_DIR=$dir cargo bench -p "$crate" --bench codec --no-run 2>&1 |
        sed -n 's/.*Executable .*(\(.*\)).*/\1/p' | tail -1)
    case $exe in /*) ;; *) exe=$PWD/$exe ;; esac
    [ -x "$exe" ] || { echo "perf: could not build the $build benchmark" >&2; exit 1; }
    BENCH=$side BENCH_REPS=${PERF_REPS:-3} taskset -c "${BENCH_CPU:-2}" \
        "$perf" record -q -e "$ev" -F 4999 -o target/perf.data "$exe" >/dev/null 2>&1
    "$perf" script -i target/perf.data -F ip,sym,symoff,dso 2>/dev/null |
        python3 -I -c '
import collections, re, subprocess, sys
exe, mod, head, top = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
base = {}
for line in subprocess.run(["nm", "-C", "--defined-only", exe], capture_output=True, text=True).stdout.splitlines():
    a, _, name = line.split(" ", 2)
    base.setdefault(name, int(a, 16))
hits, samples = collections.Counter(), 0
for line in sys.stdin:
    samples += 1
    m = re.match(r"\s*[0-9a-f]+ (.*)\+0x([0-9a-f]+) \((.*)\)$", line)
    if m and m.group(3) == exe and m.group(1) in base:
        hits[base[m.group(1)] + int(m.group(2), 16)] += 1
addrs = sorted(hits)
out = subprocess.run(["addr2line", "-a", "-i", "-f", "-C", "-e", exe], capture_output=True, text=True,
                     input="\n".join(f"{a:x}" for a in addrs)).stdout.splitlines()
cost, i = collections.Counter(), 0
while i < len(out):  # per address: "0x…", then (function, file:line) pairs, innermost first
    addr, i, fns = int(out[i], 16), i + 1, []
    while i < len(out) and not out[i].startswith("0x"):
        fns.append(out[i]); i += 2
    own = [f for f in fns if f.startswith(mod + "::") and "::basic_op::" not in f]
    if own:
        cost[re.sub(r"::\{closure#\d+\}", "", own[0])] += hits[addr]
total = sum(cost.values())
print(f"{head}: {100 * total / max(samples, 1):.0f}% of {samples} samples in the Rust codec, self cycles:")
for f, n in cost.most_common(top):
    print(f"{100 * n / max(total, 1):5.1f}%  {f}")
' "$exe" "${crate//-/_}" "$crate $side $build${BENCH_MODE:+ mode $BENCH_MODE}" "${TOP:-10}"
done
