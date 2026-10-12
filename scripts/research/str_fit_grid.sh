#!/usr/bin/env bash
# Run the repeat-tract fitting bench over a grid of cells, one method per process.
#
# The bench (examples/ng_str_fit_bench.rs) draws one motif length's repeat classes from a known model
# and scores each fitting method on it. One run of it is one cell — a motif length, a number of
# samples, a read depth, a truth and a seed — and one method. This script runs a grid of them and
# leaves one rows file per run in OUT, which scripts/research/str_fit_summarise.py turns into tables
# (doc/devel/ng/research/str_fit_research_plan.md §2, §3.1).
#
# Usage (every setting is an environment variable; the defaults are the plan's three required
# corners, for the counting fit, on both default truths' PCR variant):
#
#   OUT=tmp/str_fit_grid METHODS="A B" SEEDS="1 2 3" PERIODS="1 2" \
#     CELLS="1:3 2000:3 1:300" ./scripts/research/str_fit_grid.sh
#
#   CELLS     samples:depth pairs, or "full" for the plan's grid (1, 10, 60, 500, 2000 samples ×
#             3, 10, 30, 300 reads)
#   METHODS   method names the bench knows (A, B, B_within, …)
#   TRUTH     pcr, pcr_free, or a TOML model file
#   LOCI      tomato2 (as tomato2 has them) or a number a class
#   LOCI_AT_MOST, FOLDS, HELDOUT_LOCI
#             passed to the bench as --loci-at-most, --folds, --heldout-loci when set
#   EXTRA     any further bench arguments
#
# **Build first**; the script refuses to run without a binary and takes the newer of
# target/release and target-container/release (CLAUDE.md, "Container vs. host"):
#   cargo build --release --features bench-fixtures --example ng_str_fit_bench
#
# A run whose rows file already exists is skipped, so an interrupted grid resumes where it stopped.
# Each run's peak memory is read by /usr/bin/time (-l on macOS, -v on Linux) and appended to its
# rows file as `method = <name>, parameter = peak_rss_mib`: the whole process, the drawn cohort
# included, so a method's own share is that minus the draw's (run the same cell with METHODS=none
# to see the draw alone).

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT="${OUT:-$REPO_ROOT/tmp/str_fit_grid}"
METHODS="${METHODS:-B}"
SEEDS="${SEEDS:-1 2 3}"
PERIODS="${PERIODS:-1 2}"
CELLS="${CELLS:-1:3 2000:3 1:300}"
TRUTH="${TRUTH:-pcr}"
LOCI="${LOCI:-tomato2}"

if [[ "$CELLS" == "full" ]]; then
    CELLS=""
    for samples in 1 10 60 500 2000; do
        for depth in 3 10 30 300; do
            CELLS+="$samples:$depth "
        done
    done
fi

binary=""
for candidate in "$REPO_ROOT/target-container/release/examples/ng_str_fit_bench" \
                 "$REPO_ROOT/target/release/examples/ng_str_fit_bench"; do
    if [[ -x "$candidate" && ( -z "$binary" || "$candidate" -nt "$binary" ) ]]; then
        binary="$candidate"
    fi
done
if [[ -z "$binary" ]]; then
    echo "not built: no ng_str_fit_bench in target/ or target-container/; build it with" >&2
    echo "  cargo build --release --features bench-fixtures --example ng_str_fit_bench" >&2
    exit 1
fi
echo "using $binary" >&2

if [[ "$(uname)" == "Darwin" ]]; then
    time_flag="-l"
else
    time_flag="-v"
fi

# The peak resident memory /usr/bin/time reports, in MiB.
peak_mib() {
    local report="$1"
    if [[ "$time_flag" == "-l" ]]; then
        # macOS: "<bytes>  maximum resident set size"
        awk '/maximum resident set size/ { printf "%.1f", $1 / 1048576 }' "$report"
    else
        # Linux: "Maximum resident set size (kbytes): <kib>"
        awk -F': ' '/Maximum resident set size/ { printf "%.1f", $2 / 1024 }' "$report"
    fi
}

mkdir -p "$OUT"
for cell in $CELLS; do
    samples="${cell%%:*}"
    depth="${cell##*:}"
    for period in $PERIODS; do
        for seed in $SEEDS; do
            for method in $METHODS; do
                name="p${period}_n${samples}_d${depth}_$(basename "$TRUTH" .toml)_s${seed}_${method}"
                rows="$OUT/$name.tsv"
                if [[ -s "$rows" ]]; then
                    echo "skip $name: done" >&2
                    continue
                fi
                arguments=(--period "$period" --samples "$samples" --depth "$depth"
                           --seed "$seed" --truth "$TRUTH" --loci "$LOCI"
                           --out "$rows" --truth-out "$OUT/${name%_*}.truth.tsv")
                if [[ "$method" == "none" ]]; then
                    arguments+=(--methods "" --folds 0)
                else
                    arguments+=(--methods "$method")
                fi
                [[ -n "${LOCI_AT_MOST:-}" ]] && arguments+=(--loci-at-most "$LOCI_AT_MOST")
                [[ -n "${FOLDS:-}" ]] && arguments+=(--folds "$FOLDS")
                [[ -n "${HELDOUT_LOCI:-}" ]] && arguments+=(--heldout-loci "$HELDOUT_LOCI")
                # shellcheck disable=SC2206
                [[ -n "${EXTRA:-}" ]] && arguments+=($EXTRA)
                echo "run $name" >&2
                started=$(date +%s)
                if /usr/bin/time "$time_flag" -o "$OUT/$name.time" \
                        "$binary" "${arguments[@]}" 2> "$OUT/$name.log"; then
                    peak="$(peak_mib "$OUT/$name.time")"
                    cell_columns="$(sed -n 2p "$rows" | cut -f1-6)"
                    printf '%s\t%s\tall\tpeak_rss_mib\tNA\t%s\tNA\n' \
                        "$cell_columns" "$method" "${peak:-NA}" >> "$rows"
                    echo "  done in $(( $(date +%s) - started )) s, peak ${peak:-NA} MiB" >&2
                else
                    echo "  FAILED after $(( $(date +%s) - started )) s; see $OUT/$name.log" >&2
                fi
            done
        done
    done
done
