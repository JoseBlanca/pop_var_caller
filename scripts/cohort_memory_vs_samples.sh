#!/usr/bin/env bash
#
# cohort_memory_vs_samples.sh — how a cohort run's peak memory grows with the
# number of samples in it.
#
# The caller opens one `.psp` per sample and holds all of them for the whole
# run, so any cost paid per open file is multiplied by the cohort size. This
# script runs the same calling job at several cohort sizes and records peak
# resident memory, which separates the two things that make up the peak:
#
#   the slope     memory paid per sample — it is what decides whether a
#                 thousand-sample cohort fits on a machine;
#   the intercept memory the run needs whatever the cohort size.
#
# Both come out of a straight-line fit over the rows this writes.
#
# Usage:
#   scripts/cohort_memory_vs_samples.sh \
#       --psp-dir <dir of .psp files> --reference <reference.fa> \
#       --out-dir <where results go> [--sizes 1,5,10,25,50] [--threads 4]
#
# Portability: every path is an argument; the binary is taken as the newer of
# `target/release` and `target-container/release` (a machine with a container
# runtime builds into the second, one without into the first); inputs are
# checked before the first long run; results go to a file; and the machine and
# its core count are written into the results header, because a memory figure
# without the machine on it is not a measurement.

set -euo pipefail

PSP_DIR=""; REFERENCE=""; OUT_DIR=""; SIZES="1,5,10,25,50"; THREADS=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --psp-dir)   PSP_DIR="$2";   shift 2 ;;
        --reference) REFERENCE="$2"; shift 2 ;;
        --out-dir)   OUT_DIR="$2";   shift 2 ;;
        --sizes)     SIZES="$2";     shift 2 ;;
        --threads)   THREADS="$2";   shift 2 ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

for required in PSP_DIR REFERENCE OUT_DIR; do
    [[ -n "${!required}" ]] || { echo "missing --${required//_/-}" >&2; exit 2; }
done
[[ -d "$PSP_DIR" ]]       || { echo "no such directory: $PSP_DIR" >&2; exit 1; }
[[ -f "$REFERENCE" ]]     || { echo "no such file: $REFERENCE" >&2; exit 1; }
[[ -f "$REFERENCE.fai" ]] || { echo "reference index missing: $REFERENCE.fai" >&2; exit 1; }

# `mapfile` would be the obvious way to read this and does not exist in bash
# 3.2, which is what macOS ships.
PSPS=()
while IFS= read -r line; do PSPS+=("$line"); done \
    < <(find "$PSP_DIR" -maxdepth 1 -name '*.psp' | sort)
(( ${#PSPS[@]} )) || { echo "no .psp files in $PSP_DIR" >&2; exit 1; }

newest_binary() {
    local name="$1" newest="" candidate
    for candidate in "target/release/$name" "target-container/release/$name" \
                     "target/release/examples/$name" "target-container/release/examples/$name"; do
        [[ -x "$candidate" ]] || continue
        [[ -z "$newest" || "$candidate" -nt "$newest" ]] && newest="$candidate"
    done
    printf '%s' "$newest"
}

CALLER="$(newest_binary pop_var_caller)"
[[ -n "$CALLER" ]] || {
    echo "pop_var_caller not built. Build it first:" >&2
    echo "  ./scripts/dev.sh cargo build --release   # with a container runtime" >&2
    echo "  cargo build --release                    # without one" >&2
    exit 1
}

HERE="$(cd "$(dirname "$0")" && pwd)"
PEAK_RSS="$HERE/peak_rss.sh"
[[ -x "$PEAK_RSS" ]] || { echo "missing $PEAK_RSS" >&2; exit 1; }

mkdir -p "$OUT_DIR"
RESULTS="$OUT_DIR/memory_vs_samples.tsv"
{
    echo "# cohort peak memory against sample count"
    echo "# host	$(hostname)"
    echo "# uname	$(uname -srm)"
    echo "# cores	$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo unknown) logical"
    echo "# threads	$THREADS"
    echo "# psp_dir	$(cd "$PSP_DIR" && pwd)"
    echo "# available	${#PSPS[@]} .psp files"
    echo "# reference	$REFERENCE"
    echo "# caller	$CALLER"
    printf 'n_samples\tpeak_rss_mb\twall_s\tvcf_records\n'
} > "$RESULTS"

IFS=',' read -ra SIZE_LIST <<< "$SIZES"
for n in "${SIZE_LIST[@]}"; do
    if (( n > ${#PSPS[@]} )); then
        echo ">> skipping N=$n: only ${#PSPS[@]} files available" >&2
        continue
    fi
    subset=("${PSPS[@]:0:n}")
    vcf="$OUT_DIR/calls_n${n}.vcf"
    timing="$OUT_DIR/time_n${n}.log"

    echo ">> N=$n ..." >&2
    "$PEAK_RSS" "$timing" "$CALLER" var-calling \
        --reference "$REFERENCE" --output "$vcf" --threads "$THREADS" \
        "${subset[@]}" >"$OUT_DIR/call_n${n}.log" 2>&1
    read -r rss wall < "$timing"

    records=$(grep -cv '^#' "$vcf" 2>/dev/null || echo 0)
    printf "%s\t%s\t%s\t%s\n" "$n" "$rss" "$wall" "$records" >> "$RESULTS"
    echo ">> N=$n: peak ${rss} MB, wall ${wall} s, $records records" >&2
    rm -f "$vcf"
done

echo >&2
cat "$RESULTS"
