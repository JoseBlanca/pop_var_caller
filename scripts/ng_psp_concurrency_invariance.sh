#!/bin/sh
# **The psp route's VCF does not depend on how many threads it was given**
# (doc/devel/ng/spec/run_streaming.md §12.2). The same stored cohort is called at 1, 2, 4 and 8
# threads and the four files are compared byte for byte.
#
# **A script and not a test, and the reason is the flag.** `--threads` builds rayon's *global*
# pool, which a process may build once: a unit test sweeping thread counts would build the pool
# on its first call and silently run every later count at the first one's width, reporting a
# sweep it did not do. One invocation a thread count is the only way to measure it.
#
#   scripts/ng_psp_concurrency_invariance.sh <reference.fa> <catalog.parquet> <out-dir> <psp>...
set -eu

if [ "$#" -lt 4 ]; then
  sed -n '2,14p' "$0"
  exit 2
fi

reference=$1; catalog=$2; out=$3
shift 3

root=$(cd "$(dirname "$0")/.." && pwd)
bin=""
for candidate in "$root/target-container/release/pop_var_caller_exp" \
                 "$root/target/release/pop_var_caller_exp"; do
  if [ -x "$candidate" ] && { [ -z "$bin" ] || [ "$candidate" -nt "$bin" ]; }; then
    bin=$candidate
  fi
done
if [ -z "$bin" ]; then
  echo "no release build of pop_var_caller_exp; build one first" >&2
  exit 1
fi

mkdir -p "$out"
stored=""
for psp in "$@"; do
  stored="$stored --psp $psp"
done

for threads in 1 2 4 8; do
  # **The same --output path every time, renamed after.** The VCF header echoes the command
  # line, so a different output name would make a different header for identical calls and the
  # comparison would fail on a filename.
  # shellcheck disable=SC2086
  "$bin" call-from-psps \
    --reference "$reference" --catalog "$catalog" $stored \
    --defaults --threads "$threads" --output "$out/run.vcf" > "$out/threads-$threads.log" 2>&1
  mv "$out/run.vcf" "$out/threads-$threads.vcf"
  mv "$out/run.parameters.toml" "$out/threads-$threads.parameters.toml"
  printf 'threads %s: ' "$threads"; grep -vc '^#' "$out/threads-$threads.vcf"
done

# **`##commandline` carries `--threads N`, so it is the one line that differs by construction.**
# Measured the hard way: without this the four files differ on that line alone and the comparison
# reports a failure that is the flag being recorded, not the calls moving.
for threads in 1 2 4 8; do
  grep -v '^##commandline=' "$out/threads-$threads.vcf" > "$out/threads-$threads.comparable"
done

status=0
for threads in 2 4 8; do
  if cmp -s "$out/threads-1.comparable" "$out/threads-$threads.comparable"; then
    echo "1 and $threads threads: IDENTICAL"
  else
    echo "1 and $threads threads: DIFFERENT"
    diff "$out/threads-1.comparable" "$out/threads-$threads.comparable" | head -20
    status=1
  fi
done
exit "$status"
