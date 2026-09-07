#!/usr/bin/env bash
# **The hidden-duplication filter on real reads** — plan steps D1 and D2 of
# doc/devel/ng/impl_plan/hidden_paralog_filter.md.
#
# One cohort of alignment files, called three ways, and the numbers each run reports:
#
#   off      --paralog-fdr 0            the run before this work, byte for byte
#   drop     --paralog-fdr <target>     flagged records leave the file
#   tag      --paralog-fdr <target> --paralog-filter-tag
#                                       flagged records stay, on the hiddenParalog filter
#
# **Why all three every time.** The off run is what the other two are measured against — how many
# records went, and whether anything but the filter's own lines moved. The tag run is where a
# flagged record can still be read: a dropped record leaves no trace of its ratio in the file, so
# anything counted per flagged record has to be taken from the tagged run.
#
# Usage:
#
#   scripts/ng_paralog_filter_runs.sh <reference.fa> <catalog.parquet> <regions.bed> \
#       <output-dir> <alignment.cram>...
#
# Environment:
#
#   NG_PARALOG_FDR   the target false-discovery rate for the two on runs (default 0.01,
#                    spec §3.6's own value)
#   NG_THREADS       threads for the calling pass (default 4)
#
# It expects a release build of `pop_var_caller_exp` in `target-container/release` or
# `target/release` and takes whichever is newer, because a machine with no container runtime
# builds to the second (see CLAUDE.md).
#
# **bash, not sh, and python3 for the measuring.** The container's `/bin/sh` is dash, and it
# carries no GNU `time`, so wall and peak resident memory are read with `getrusage` over the
# child — which reports the peak the kernel recorded rather than a sample of it.
set -euo pipefail

if (( $# < 5 )); then
  # Every comment line from the shebang to the first line that is not one, so the usage cannot
  # fall off the end of a fixed range as this header grows.
  sed -n '2,${/^#/!q;p;}' "$0"
  exit 2
fi

reference=$1; catalog=$2; regions=$3; out=$4
shift 4

fdr=${NG_PARALOG_FDR:-0.01}
threads=${NG_THREADS:-4}

root=$(cd "$(dirname "$0")/.." && pwd)
bin=""
for candidate in "$root/target-container/release/pop_var_caller_exp" \
                 "$root/target/release/pop_var_caller_exp"; do
  if [[ -x "$candidate" ]] && { [[ -z "$bin" ]] || [[ "$candidate" -nt "$bin" ]]; }; then
    bin=$candidate
  fi
done
if [[ -z "$bin" ]]; then
  echo "no release build of pop_var_caller_exp; build one first" >&2
  exit 1
fi

rm -rf "$out"
mkdir -p "$out"

alignments=()
samples=0
for cram in "$@"; do
  alignments+=(--alignment "$cram")
  samples=$((samples + 1))
done

echo "cohort: $samples sample(s); target false-discovery rate $fdr; $threads thread(s)"
# **Which build produced the numbers.** `$bin` is chosen by modification time between two target
# directories, so the output of this script cannot be attributed to a build unless it says which.
echo "binary: $bin ($(date -r "$bin" '+%Y-%m-%d %H:%M:%S'))"

# **Wall and peak resident, from the kernel rather than from a poll.** `ru_maxrss` under
# `RUSAGE_CHILDREN` is the largest resident size any waited-for child reached, so it is the peak
# and not a sample that may have missed it.
#
# **Each call is its own python process and waits for exactly one child**, which is what makes the
# mark after it that child's own peak rather than the largest of several. The assertion below is
# how that stays true: a future edit that ran two things under one `measure` would report the
# larger of the two for both, silently.
measure() {
  local where=$1; shift
  python3 - "$where" "$@" <<'PYTHON'
import resource, subprocess, sys, time

where, command = sys.argv[1], sys.argv[2:]
before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
assert before == 0, (
    f"this process has already waited for a child (mark {before}), so ru_maxrss after the next "
    "one is the larger of the two and not this run's peak"
)
started = time.monotonic()
with open(where, "wb") as log:
    finished = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
wall = time.monotonic() - started
peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
# **`ru_maxrss` is kilobytes on Linux and bytes on macOS**, and this script is meant to be run in
# the dev container, which is Linux. A host run would otherwise print a number a thousand times
# too large without saying so.
if sys.platform == "darwin":
    peak //= 1024
elif not sys.platform.startswith("linux"):
    sys.exit(f"ru_maxrss units are not known for {sys.platform}; run this in the dev container")
print(f"wall {wall:.1f}s, peak resident {peak / 1024:.1f} MB")
sys.exit(finished.returncode)
PYTHON
}

run_one() {
  local side=$1; shift
  # **Every run writes to the same basename.** The header records the parameters file beside the
  # VCF by name (`##parametersFile=`), so a run written as `off.vcf` and one written as
  # `drop.vcf` would differ in the header whatever the filter did — which would make the off run
  # uncomparable with the standing pre-filter baseline. Each writes `run.vcf` and is moved after.
  printf '%s: ' "$side"
  measure "$out/$side.log" \
    "$bin" call-from-alignments \
    --reference "$reference" --catalog "$catalog" "${alignments[@]}" --regions "$regions" \
    --defaults --threads "$threads" "$@" --output "$out/run.vcf" || {
      echo "the $side run failed:" >&2; tail -40 "$out/$side.log" >&2; exit 1; }
  mv "$out/run.vcf" "$out/$side.vcf"
  mv "$out/run.parameters.toml" "$out/$side.parameters.toml"
}

echo
echo "=== the three runs ==="
run_one off --paralog-fdr 0
run_one drop --paralog-fdr "$fdr"
run_one tag --paralog-fdr "$fdr" --paralog-filter-tag

# **The off run is the one compared against the standing pre-filter baseline**, so it is the only
# one whose hash means anything; the other two carry the filter's own lines by design.
grep -v '^##commandline=' "$out/off.vcf" > "$out/off.comparable"

echo
echo "=== records written ==="
declare -A written
for side in off drop tag; do
  # `grep -vc` exits 1 on a count of zero, which inside `$( )` would not trip `set -e`.
  written[$side]=$(grep -vc '^#' "$out/$side.vcf" || true)
  printf '%-5s %s\n' "$side" "${written[$side]}"
done
# **A run that wrote nothing satisfies every check below**, so say so rather than printing three
# clean zeros. It is not necessarily a fault — an interval with no variants in it does this — but
# nothing measured on such a run means anything.
if [[ ${written[off]} -eq 0 ]]; then
  echo "WARNING: the off run wrote no records at all — nothing below is a measurement of anything"
fi
tagged=$(awk -F'\t' '!/^#/ { n = split($7, ids, ";"); for (i = 1; i <= n; i++) if (ids[i] == "hiddenParalog") { print; next } }' \
  "$out/tag.vcf" | wc -l | tr -d ' ')
printf 'of the tag run, carrying hiddenParalog: %s\n' "$tagged"
if [[ $tagged -eq 0 && ${written[off]} -gt 0 ]]; then
  echo "NOTE: nothing was flagged, so every relation below is the empty case"
fi
printf 'the off run, sha256 of everything but ##commandline: '
shasum -a 256 "$out/off.comparable" | cut -d' ' -f1

echo
echo "=== what the filter told the operator ==="
for side in drop tag; do
  echo "--- $side ---"
  # **An empty report is a failure, not a quiet section.** `sed` prints nothing and exits 0 when
  # the marker is absent, so a reworded first line would silently take π, the cut, the
  # convergence flag and the drop count out of every report written from this output.
  report=$(sed -n '/^hidden-duplication filter:/,$p' "$out/$side.log")
  if [[ -z $report ]]; then
    echo "the $side run printed no filter report — has the report's first line been reworded?" >&2
    exit 1
  fi
  printf '%s\n' "$report"
done

echo
echo "=== the header's calibration line ==="
for side in off drop tag; do
  printf '%-5s ' "$side"
  line=$(grep '^##paralogFilter=' "$out/$side.vcf" || true)
  # The off run must not carry it and the other two must; either way round is a broken run.
  if [[ $side == off ]]; then
    [[ -z $line ]] || { echo "the off run carries a calibration line: $line" >&2; exit 1; }
    echo "(none — the filter did not run, which is right)"
  else
    [[ -n $line ]] || { echo "the $side run carries no calibration line" >&2; exit 1; }
    printf '%s\n' "$line"
  fi
done

echo
echo "=== the three file relations, on real reads ==="
python3 "$root/scripts/ng_paralog_filter_relations.py" \
  "$out/off.vcf" "$out/drop.vcf" "$out/tag.vcf"
