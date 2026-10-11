#!/usr/bin/env bash
# **The stuttering benchmark's first pass**: 16 runs, PCR-free and PCR libraries, 1 and 50 samples,
# 5, 15, 30 and 60 reads a position (doc/devel/ng/spec/str_stutter_benchmark.md §4.1; owner,
# 2026-10-09).
#
#   scripts/ng_str_stutter_first_pass.sh <grid-dir> <bin-dir>
#
# Run it inside the dev container. <bin-dir> holds a frozen copy of pvc,
# ng_str_stutter_benchmark and ng_typed_region_dump, so the tree can be rebuilt while the grid runs.
#
# Each cohort — one design and sample count, at the simulator's default seed — runs at 30× first,
# and its fit is reused at 5×, 15× and 60× (`FIT_FROM`; the samples are the same at every depth).
# A run that fails is logged and the grid moves on; rerunning the script resumes every run where
# it stopped. Runs are `<grid-dir>/<design>_<samples>s_<depth>x`; progress in `<grid-dir>/grid.log`.
#
# **A simulator built since 2026-10-11 lays 500 STR bundles by default**, and the segments then sit
# elsewhere along the genome: the cohorts are not the 2026-10-10 first pass's
# (doc/devel/reports/implementations/str_stutter_benchmark_first_pass_2026-10-10.md). Those come
# back only with `bundles=0`, which this script does not pass. A fit from a run made before bundles
# is refused by `FIT_FROM`, because the settings line it compares now names `bundles`.
set -uo pipefail

grid=${1:?usage: $0 <grid-dir> <bin-dir>}
bins=${2:?usage: $0 <grid-dir> <bin-dir>}
root=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$grid"
grid=$(cd "$grid" && pwd)
bins=$(cd "$bins" && pwd)
log=$grid/grid.log

# Four HaplotypeCaller processes of 2 GB fit the 16 GB container beside everything else.
export GATK_JOBS=${GATK_JOBS:-4} GATK_HEAP=${GATK_HEAP:-2g} BIN_DIR=$bins

run() {
    local dir=$1 fit_from=$2
    shift 2
    echo "$(date -u +%FT%TZ) start $(basename "$dir") ${fit_from:+(fit from $(basename "$fit_from"))}" >> "$log"
    if FIT_FROM=$fit_from "$root/scripts/ng_str_stutter_benchmark.sh" "$dir" "$@" > "$dir.log" 2>&1; then
        echo "$(date -u +%FT%TZ) done  $(basename "$dir")" >> "$log"
    else
        echo "$(date -u +%FT%TZ) FAILED $(basename "$dir") — see $dir.log" >> "$log"
        return 1
    fi
}

for cohort in "pcr 50" "pcr 1" "pcr_free 50" "pcr_free 1"; do
    read -r design samples <<< "$cohort"
    first=$grid/${design}_${samples}s_30x
    run "$first" "" design="$design" samples="$samples" depth=30 || continue
    for depth in 5 15 60; do
        run "$grid/${design}_${samples}s_${depth}x" "$first" \
            design="$design" samples="$samples" depth="$depth"
    done
done
echo "$(date -u +%FT%TZ) grid finished" >> "$log"
