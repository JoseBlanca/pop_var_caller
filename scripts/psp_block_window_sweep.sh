#!/usr/bin/env bash
#
# psp_block_window_sweep.sh — how much of a cohort run's peak memory follows
# the per-sample store, and how much is a floor the store cannot touch.
#
# The `.psp` block is the cohort reader's decode unit: a reader inflates a whole
# block before it can hand out a record, so the block's span in reference
# coordinates sets how much every open sample forces the reader to hold. This
# script rewrites one cohort's `.psp` files at several block spans, runs the
# caller on each, and records peak resident memory and wall time.
#
# What the result says: the part of the peak that moves with the span is the
# part a store redesign can address. The part that does not move is the floor —
# the per-sample columns assembled for a cohort chunk, the merger's projections,
# the genotype fit — and no encoding change reaches it.
#
# Usage:
#   scripts/psp_block_window_sweep.sh \
#       --psp-dir  <dir of .psp files>     \
#       --reference <reference.fa>          \
#       --out-dir  <where results go>       \
#       [--n-samples 50] [--threads 4]      \
#       [--windows 5000,20000,80000]        \
#       [--keep-rechunked]
#
# Portability, because this runs on more than one machine:
#   - every path is an argument; nothing is hard-coded;
#   - the binaries are looked for in BOTH `target/release` and
#     `target-container/release`, taking whichever is newer — a machine with a
#     container runtime builds into the second, a machine without builds into
#     the first, and a script that checks only one silently runs a stale build;
#   - inputs are checked before the first long run, not part-way through;
#   - results are written to a file, so a run that takes hours survives a closed
#     terminal;
#   - the machine, its cores and the resolved paths go in the results header. A
#     memory figure without the machine on it is not a measurement.

set -euo pipefail

PSP_DIR=""
REFERENCE=""
OUT_DIR=""
N_SAMPLES=50
THREADS=4
WINDOWS="5000,20000,80000"
KEEP_RECHUNKED=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --psp-dir)        PSP_DIR="$2"; shift 2 ;;
        --reference)      REFERENCE="$2"; shift 2 ;;
        --out-dir)        OUT_DIR="$2"; shift 2 ;;
        --n-samples)      N_SAMPLES="$2"; shift 2 ;;
        --threads)        THREADS="$2"; shift 2 ;;
        --windows)        WINDOWS="$2"; shift 2 ;;
        --keep-rechunked) KEEP_RECHUNKED=1; shift ;;
        *) echo "unknown argument: $1" >&2; exit 2 ;;
    esac
done

for required in PSP_DIR REFERENCE OUT_DIR; do
    if [[ -z "${!required}" ]]; then
        echo "missing --${required//_/-}" >&2
        exit 2
    fi
done

# ---------------------------------------------------------------------------
# Check every input before the first long run.
# ---------------------------------------------------------------------------

[[ -d "$PSP_DIR" ]]       || { echo "no such directory: $PSP_DIR" >&2; exit 1; }
[[ -f "$REFERENCE" ]]     || { echo "no such file: $REFERENCE" >&2; exit 1; }
[[ -f "$REFERENCE.fai" ]] || { echo "reference index missing: $REFERENCE.fai" >&2; exit 1; }

# `mapfile` would be the obvious way to read this list, and it does not exist in
# bash 3.2 — which is what macOS ships and what this was first run on.
PSPS=()
while IFS= read -r line; do
    PSPS+=("$line")
done < <(find "$PSP_DIR" -maxdepth 1 -name '*.psp' | sort)

if (( ${#PSPS[@]} < N_SAMPLES )); then
    echo "asked for $N_SAMPLES samples, found ${#PSPS[@]} .psp files in $PSP_DIR" >&2
    exit 1
fi
SUBSET=("${PSPS[@]:0:N_SAMPLES}")

# Take the newer of the two build trees for each binary. Neither is
# authoritative: the container build writes to `target-container`, a direct
# cargo build to `target`, and a machine can have both.
newest_binary() {
    local name="$1" newest="" candidate
    for candidate in "target/release/$name" "target-container/release/$name" \
                     "target/release/examples/$name" "target-container/release/examples/$name"; do
        [[ -x "$candidate" ]] || continue
        if [[ -z "$newest" || "$candidate" -nt "$newest" ]]; then
            newest="$candidate"
        fi
    done
    printf '%s' "$newest"
}

CALLER="$(newest_binary pop_var_caller)"
RECHUNK="$(newest_binary psp_rechunk)"

if [[ -z "$CALLER" ]]; then
    echo "pop_var_caller not built. Build it first:" >&2
    echo "  ./scripts/dev.sh cargo build --release          # with a container runtime" >&2
    echo "  cargo build --release                           # without one" >&2
    exit 1
fi
if [[ -z "$RECHUNK" ]]; then
    echo "psp_rechunk not built. Build it first:" >&2
    echo "  cargo build --release --example psp_rechunk" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

HERE="$(cd "$(dirname "$0")" && pwd)"
PEAK_RSS="$HERE/peak_rss.sh"
[[ -x "$PEAK_RSS" ]] || { echo "missing $PEAK_RSS" >&2; exit 1; }

# ---------------------------------------------------------------------------
# The header, then one row per block span.
# ---------------------------------------------------------------------------

RESULTS="$OUT_DIR/block_window_sweep.tsv"
{
    echo "# psp block-window sweep"
    echo "# host	$(hostname)"
    echo "# uname	$(uname -srm)"
    echo "# cores	$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo unknown) logical"
    echo "# threads	$THREADS"
    echo "# samples	$N_SAMPLES"
    echo "# psp_dir	$(cd "$PSP_DIR" && pwd)"
    echo "# reference	$REFERENCE"
    echo "# caller	$CALLER"
    echo "# rechunk	$RECHUNK"
    printf 'window_bp\tpsp_mb\tpeak_rss_mb\twall_s\tvcf_records\n'
} > "$RESULTS"

echo ">> sweeping block windows: $WINDOWS" >&2
echo ">> results: $RESULTS" >&2

IFS=',' read -ra WINDOW_LIST <<< "$WINDOWS"
for window in "${WINDOW_LIST[@]}"; do
    stage="$OUT_DIR/psp_${window}"
    if [[ ! -d "$stage" ]]; then
        echo ">> rechunking $N_SAMPLES samples at ${window} bp ..." >&2
        mkdir -p "$stage"
        "$RECHUNK" "$window" "$stage" "${SUBSET[@]}" >"$OUT_DIR/rechunk_${window}.log" 2>&1
    else
        echo ">> reusing existing $stage" >&2
    fi

    psp_mb=$(du -sk "$stage" | awk '{printf "%.1f", $1 / 1024}')
    vcf="$OUT_DIR/calls_${window}.vcf"
    timing="$OUT_DIR/time_${window}.log"

    echo ">> calling at ${window} bp ..." >&2
    "$PEAK_RSS" "$timing" "$CALLER" var-calling \
        --reference "$REFERENCE" \
        --output "$vcf" \
        --threads "$THREADS" \
        "$stage"/*.psp \
        >"$OUT_DIR/call_${window}.log" 2>&1
    read -r rss wall < "$timing"

    records=$(grep -cv '^#' "$vcf" || true)

    printf '%s\t%s\t%s\t%s\t%s\n' \
        "$window" "$psp_mb" "$rss" "$wall" "$records" >> "$RESULTS"
    echo ">> ${window} bp: peak ${rss} MB, wall ${wall} s, psp ${psp_mb} MB, $records records" >&2

    if (( ! KEEP_RECHUNKED )); then
        rm -rf "$stage"
    fi
done

echo >&2
echo ">> done. $RESULTS" >&2
cat "$RESULTS"
