#!/bin/sh
# **The caller's speed and memory, measured the same way before and after a change to its
# floating-point maths** (doc/devel/implementation_plans/portable_float.md, steps A3 and B).
#
# Three parts, each writing under <output-dir>:
#
#   benches  — the five criterion benches, each run twice in a row so the second run's
#              "change" line shows run-to-run noise. Criterion keeps its own history under
#              the target directory; the logs here are what the report quotes.
#   runs     — call-from-alignments, generate-psps and call-from-psps on the given CRAMs,
#              each repeated <repeats> times, with wall seconds and peak resident memory
#              appended one line a run to runs.tsv.
#   fit      — estimate-parameters on the given psps, repeated <repeats> times, into fit.tsv.
#              Kept apart from `runs` because one fit costs as much as many calling runs.
#
# Usage:
#
#   scripts/portable_float_baseline.sh benches <output-dir>
#   scripts/portable_float_baseline.sh runs <output-dir> <repeats> <reference.fa> \
#       <catalog.parquet> <regions.bed> <alignment.cram>...
#   scripts/portable_float_baseline.sh fit <output-dir> <repeats> <reference.fa> \
#       <catalog.parquet> <sample.psp>...
#
# Run it natively on macOS, and inside the dev container for Linux:
#
#   DEV_EXTRA_MOUNT=<dir holding the crams and catalog> \
#     ./scripts/dev.sh scripts/portable_float_baseline.sh runs ...
#
# **Native parallelism.** The two calling commands take `--threads 0`, every core the machine
# offers; estimate-parameters has no thread flag and uses every core. generate-psps has no thread
# flag and walks its samples one after another, so it runs as one invocation over the whole
# cohort, as the identity oracle runs it. The benches keep the parallelism each one sets for
# itself.
#
# **Peak memory** comes from the operating system's account of the finished process:
# `/usr/bin/time -l` on macOS (bytes), and `getrusage` of the child through python3 on Linux,
# whose container has no `/usr/bin/time` (kilobytes). Both are written out in megabytes.
set -eu

if [ "$#" -lt 2 ]; then
  sed -n '2,${/^#/!q;p;}' "$0"
  exit 2
fi

mode=$1; out=$2
shift 2
root=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$out"

# Wall seconds and peak megabytes of one command, printed as "seconds<TAB>megabytes". The
# command's own output goes to the log file named first.
measure() {
  log=$1; shift
  if [ "$(uname)" = Darwin ]; then
    if ! /usr/bin/time -l "$@" > "$log" 2> "$log.time"; then
      echo "$1 failed; see $log and $log.time" >&2
      exit 1
    fi
    bytes=$(awk '/maximum resident set size/ {print $1}' "$log.time")
    seconds=$(awk '/ real / {print $1}' "$log.time")
    awk -v s="$seconds" -v b="$bytes" 'BEGIN { printf "%s\t%.1f\n", s, b / 1048576 }'
  else
    python3 - "$log" "$@" <<'PY'
import resource, subprocess, sys, time
log = sys.argv[1]
start = time.monotonic()
with open(log, "w") as handle:
    status = subprocess.run(sys.argv[2:], stdout=handle, stderr=subprocess.STDOUT).returncode
seconds = time.monotonic() - start
if status != 0:
    sys.exit(f"{sys.argv[2]} exited {status}; see {log}")
kilobytes = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
print(f"{seconds:.2f}\t{kilobytes / 1024:.1f}")
PY
  fi
}

record_build() {
  # What was measured, so two runs weeks apart can be told apart: when, where, which compiler,
  # which commit, and — for a binary — its checksum and time of building. `git` may be unusable
  # inside the container, whose mount does not include a worktree's git directory.
  {
    echo "date: $(date '+%Y-%m-%d %H:%M:%S %z')"
    echo "machine: $(uname -sm)"
    echo "rustc: $(rustc --version 2>/dev/null || echo unknown)"
    echo "commit: $(git -C "$root" rev-parse HEAD 2>/dev/null || echo unknown)"
    if [ -n "$2" ]; then
      if command -v sha256sum > /dev/null 2>&1; then
        sum=$(sha256sum "$2" | cut -d' ' -f1)
      else
        sum=$(shasum -a 256 "$2" | cut -d' ' -f1)
      fi
      echo "binary: $2"
      echo "binary sha256: $sum"
      echo "binary modified: $(date -r "$2" '+%Y-%m-%d %H:%M:%S %z')"
    fi
  } > "$1"
}

# The newest release build of pop_var_caller that runs on this machine. The container builds to
# target-container and a native build to target (CLAUDE.md), and a Linux binary does not run on
# macOS nor the reverse, so both are tried. Never run `cargo bench` into these directories between
# a build and a measurement: see `benches` below.
find_binary() {
  found=""
  for candidate in "$root/target-container/release/pop_var_caller" \
                   "$root/target/release/pop_var_caller"; do
    if [ -x "$candidate" ] && "$candidate" --version > /dev/null 2>&1 &&
       { [ -z "$found" ] || [ "$candidate" -nt "$found" ]; }; then
      found=$candidate
    fi
  done
  if [ -z "$found" ]; then
    echo "no release build of pop_var_caller that runs on this machine; build one first" >&2
    exit 1
  fi
  echo "$found"
}

case "$mode" in
  benches)
    # **The benches build into their own directory.** `cargo bench` also rebuilds the package's
    # binary, with the bench profile and whatever features the bench asks for, into the same
    # `release` directory the release build uses — and `find_binary` then picks that build for the
    # next `runs` or `fit`. Measured on 2026-09-14: a fit and two output comparisons ran on it
    # before anyone noticed.
    CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$root/target}-bench"
    export CARGO_TARGET_DIR
    # In the container rustc is killed for memory above two parallel jobs (CLAUDE.md); a native
    # build may override this with CARGO_BUILD_JOBS.
    CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-2}"
    export CARGO_BUILD_JOBS
    record_build "$out/build.txt" ""
    for bench in ng_site_quality_perf ng_ssr_delimiter_perf ng_generic_pileup_perf ng_psp_perf; do
      for pass in 1 2; do
        echo "=== $bench, pass $pass"
        cargo bench --bench "$bench" > "$out/$bench.$pass.log" 2>&1
      done
    done
    for pass in 1 2; do
      echo "=== ng_joint_fit_perf, pass $pass"
      cargo bench --features bench-fixtures --bench ng_joint_fit_perf \
        > "$out/ng_joint_fit_perf.$pass.log" 2>&1
    done
    ;;
  runs)
    if [ "$#" -lt 5 ]; then
      echo "runs needs <repeats> <reference.fa> <catalog.parquet> <regions.bed> <alignment.cram>..." >&2
      exit 2
    fi
    repeats=$1; reference=$2; catalog=$3; regions=$4
    shift 4
    bin=$(find_binary)
    alignments=""
    for cram in "$@"; do alignments="$alignments --alignment $cram"; done

    record_build "$out/build.txt" "$bin"
    echo "# binary: $bin (see build.txt)" > "$out/runs.tsv"
    echo "# regions: $regions ($(wc -l < "$regions" | tr -d ' ') lines); alignments: $#" >> "$out/runs.tsv"
    printf 'repeat\tcommand\tseconds\tpeak_megabytes\n' >> "$out/runs.tsv"
    repeat=1
    while [ "$repeat" -le "$repeats" ]; do
      work=$out/work
      rm -rf "$work"
      mkdir -p "$work/psps"
      echo "=== repeat $repeat of $repeats"
      # shellcheck disable=SC2086
      result=$(measure "$work/alignments.log" "$bin" call-from-alignments \
        --reference "$reference" --catalog "$catalog" $alignments --regions "$regions" \
        --defaults --threads 0 --output "$work/from_alignments.vcf")
      printf '%s\tcall-from-alignments\t%s\n' "$repeat" "$result" >> "$out/runs.tsv"
      # shellcheck disable=SC2086
      result=$(measure "$work/generate.log" "$bin" generate-psps \
        --reference "$reference" --catalog "$catalog" $alignments --regions "$regions" \
        --output-dir "$work/psps" --force)
      printf '%s\tgenerate-psps\t%s\n' "$repeat" "$result" >> "$out/runs.tsv"
      stored=""
      for psp in "$work"/psps/*.psp; do stored="$stored --psp $psp"; done
      # shellcheck disable=SC2086
      result=$(measure "$work/psps.log" "$bin" call-from-psps \
        --reference "$reference" --catalog "$catalog" $stored \
        --defaults --threads 0 --output "$work/from_psps.vcf")
      printf '%s\tcall-from-psps\t%s\n' "$repeat" "$result" >> "$out/runs.tsv"
      repeat=$((repeat + 1))
    done
    cat "$out/runs.tsv"
    ;;
  fit)
    if [ "$#" -lt 4 ]; then
      echo "fit needs <repeats> <reference.fa> <catalog.parquet> <sample.psp>..." >&2
      exit 2
    fi
    repeats=$1; reference=$2; catalog=$3
    shift 3
    bin=$(find_binary)
    stored=""
    for psp in "$@"; do stored="$stored --psp $psp"; done

    record_build "$out/build.txt" "$bin"
    echo "# binary: $bin (see build.txt)" > "$out/fit.tsv"
    echo "# psps: $#" >> "$out/fit.tsv"
    printf 'repeat\tcommand\tseconds\tpeak_megabytes\n' >> "$out/fit.tsv"
    repeat=1
    while [ "$repeat" -le "$repeats" ]; do
      echo "=== repeat $repeat of $repeats"
      # shellcheck disable=SC2086
      result=$(measure "$out/fit.$repeat.log" "$bin" estimate-parameters \
        --reference "$reference" --catalog "$catalog" $stored \
        --output "$out/fit.$repeat.parameters.toml" --force)
      printf '%s\testimate-parameters\t%s\n' "$repeat" "$result" >> "$out/fit.tsv"
      repeat=$((repeat + 1))
    done
    cat "$out/fit.tsv"
    ;;
  *)
    sed -n '2,${/^#/!q;p;}' "$0"
    exit 2
    ;;
esac
