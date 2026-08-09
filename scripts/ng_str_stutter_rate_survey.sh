#!/usr/bin/env bash
# Fit the STR slippage rate for many libraries, N at a time — run this ON RICK.
#
# `ng_str_stutter_rate` fits, per (library, motif period, reference repeat count), how often a read
# shows a length other than its allele's. This drives it over a list of alignments and keeps one
# result file per library, so the run resumes rather than restarting.
#
# ## Run it in three stages, not one
#
# **The per-library cost is not known**, and that is the whole reason for staging. The fits are the
# expensive half, they scale with how many distinct locus shapes a stratum holds, and nothing has
# yet timed one library over a whole chromosome. Launching 2,475 of anything whose unit cost you
# have not measured is how the last survey spent two days.
#
#   LIMIT=1   ./scripts/ng_str_stutter_rate_survey.sh ...   # one library — read the elapsed time
#   LIMIT=20  ./scripts/ng_str_stutter_rate_survey.sh ...   # the twenty, one per project
#   LIMIT=0   ./scripts/ng_str_stutter_rate_survey.sh ...   # everything on the list
#
# Each stage reuses the previous stage's results — a library already fitted is not fitted again —
# so the stages cost nothing extra beyond the first.
#
# Usage:
#   ./ng_str_stutter_rate_survey.sh OUT.tsv REF LIBRARIES.txt
#
# where LIBRARIES.txt is one alignment path per line. On rick, and note **no `./scripts/dev.sh`**:
# that box has no container runtime (CLAUDE.md, "Container vs. host").
#
#   cargo build --release --example ng_str_stutter_rate
#   ./target/release/examples/ng_str_stutter_rate --self-check     # once per build, ~7 min
#   LIMIT=1 JOBS=20 ./scripts/ng_str_stutter_rate_survey.sh ~/tmp/stutter_rate.tsv \
#       /home/joxi/refs/S_lycopersicum_chromosomes.4.00.fa ~/tmp/libraries.txt
#
# Knobs:
#   JOBS=20        libraries fitted at once. Each is a separate single-threaded process, so this is
#                  a core count. Each also holds its own windowed reference and its own stratum
#                  tables, so watch memory on the first batch before trusting a large value.
#   LIMIT=1        how many libraries from the list to attempt this run. 0 means all of them.
#   CONTIGS=SL4.0ch01   the walk. One chromosome gives ~644k STR loci, which is ample per stratum.
#   MIN_LOCI=500   fewest loci a stratum needs before it is fitted at all. Below it the row is not
#                  emitted, because a rate fitted on a handful of loci is noise wearing a number.
#   MAX_REPEATS=30 strata above this are skipped. The cost of a fit grows with the stratum's
#                  distinct shapes, and the longest tracts are both the thinnest and the slowest.

set -euo pipefail

# ---------------------------------------------------------------------------
# What counts as a finished library
# ---------------------------------------------------------------------------
#
# **A result is finished when it holds a fitted row, not when the file exists.** This is the lesson
# the copy-floor survey paid for: 1,300 walks wrote a well-formed table with nothing in it, and
# every later pass skipped them because a file was there. Used by the worker (skip), by the summary
# (count) and by the merge (fold in) — checking only one of the three leaves the other two able to
# report a library that was never measured.
#
# awk rather than grep, because `grep -q -v` with several `-e` patterns does not behave the same in
# every grep, and a portability bug here reads as "no results" rather than as an error.
has_rows() {
    [[ -s "$1" ]] || return 1
    awk -F'\t' '
        /^#/ { next }
        $1 == "read_group" { next }
        { found = 1; exit }
        END { exit !found }
    ' "$1"
}

# ---------------------------------------------------------------------------
# The per-library worker, re-entered through `--one`
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "--one" ]]; then
    cram=$2
    key=$(basename "$cram" | tr -c 'A-Za-z0-9._-' '_')
    result="$RATE_WORK/$key.tsv"
    status="$RATE_WORK/$key.status"
    has_rows "$result" && exit 0
    rm -f "$result"
    if [[ ! -f "$cram" ]]; then echo "missing" > "$status"; exit 0; fi

    read -r -a walk_args <<< "$RATE_WALK_ARGS"
    started=$(date +%s)
    if ! "$RATE_BIN" "${walk_args[@]}" --min-loci "$RATE_MIN_LOCI" \
            --max-repeats "$RATE_MAX_REPEATS" "$RATE_REF" "$cram" \
            > "$result.partial" 2> "$RATE_WORK/$key.err"; then
        rm -f "$result.partial"
        echo "failed: $(tail -n 1 "$RATE_WORK/$key.err" 2>/dev/null | tr '\t' ' ')" > "$status"
        exit 0
    fi
    elapsed=$(( $(date +%s) - started ))
    # **The elapsed time travels with the result**, because the first question after the first batch
    # is whether the rest is affordable, and a time in a terminal scrollback is a time nobody has.
    if ! has_rows "$result.partial"; then
        rm -f "$result.partial"
        echo "fitted nothing in ${elapsed}s — no stratum reached $RATE_MIN_LOCI loci" > "$status"
        exit 0
    fi
    printf '#elapsed_seconds\t%s\n' "$elapsed" >> "$result.partial"
    mv "$result.partial" "$result"
    rm -f "$status"
    echo "  done in ${elapsed}s: $(basename "$cram")" >&2
    exit 0
fi

OUT=${1:?"usage: $0 OUT.tsv REF LIBRARIES.txt"}
REF=${2:?"usage: $0 OUT.tsv REF LIBRARIES.txt"}
LIST=${3:?"usage: $0 OUT.tsv REF LIBRARIES.txt"}

JOBS=${JOBS:-20}
LIMIT=${LIMIT:-1}
CONTIGS=${CONTIGS:-SL4.0ch01}
MIN_LOCI=${MIN_LOCI:-500}
MAX_REPEATS=${MAX_REPEATS:-30}
WORK=${WORK:-$OUT.work}

[[ -f "$REF" ]] || { echo "reference not found: $REF" >&2; exit 1; }
[[ -f "$LIST" ]] || { echo "library list not found: $LIST" >&2; exit 1; }

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SOURCE="$REPO_ROOT/examples/ng_str_stutter_rate.rs"

# Two target trees, because the container build points CARGO_TARGET_DIR at the second and a host
# build does not. Take the newer and say which — a survey running a months-old binary is how the
# last one produced 1,300 unusable results.
BIN=""
for candidate in "$REPO_ROOT/target/release/examples/ng_str_stutter_rate" \
                 "$REPO_ROOT/target-container/release/examples/ng_str_stutter_rate"; do
    [[ -x "$candidate" ]] || continue
    if [[ -z "$BIN" || "$candidate" -nt "$BIN" ]]; then BIN=$candidate; fi
done
[[ -n "$BIN" ]] || {
    echo "not built: no ng_str_stutter_rate in target/ or target-container/" >&2
    echo "  cargo build --release --example ng_str_stutter_rate" >&2
    exit 1
}
if [[ -f "$SOURCE" && "$SOURCE" -nt "$BIN" ]]; then
    echo "$BIN is older than $SOURCE — rebuild before running:" >&2
    echo "  cargo build --release --example ng_str_stutter_rate" >&2
    exit 1
fi
echo "using $BIN" >&2

mkdir -p "$WORK"

# **The work directory belongs to one set of walk settings.** Results from two different walks merge
# into a table that compares nothing, and narrowing the walk is the obvious thing to reach for when
# a run is taking too long — so make the mistake impossible rather than document it.
# A BED to walk instead of whole contigs — overrides CONTIGS. The rate needs loci per stratum, so
# a slice is only for smoke-testing the driver; a real run wants the chromosome.
REGIONS=${REGIONS:-}
if [[ -n "$REGIONS" ]]; then
    [[ -f "$REGIONS" ]] || { echo "regions BED not found: $REGIONS" >&2; exit 1; }
    WALK_ARGS="--regions $REGIONS"
else
    WALK_ARGS="--contigs $CONTIGS"
fi

PARAMS="walk=${REGIONS:-$CONTIGS} min_loci=$MIN_LOCI max_repeats=$MAX_REPEATS"
STAMP="$WORK/.rate-params"
if [[ -f "$STAMP" ]] && [[ "$(cat "$STAMP")" != "$PARAMS" ]]; then
    echo "$WORK was built with different settings:" >&2
    echo "  it holds:  $(cat "$STAMP")" >&2
    echo "  you asked: $PARAMS" >&2
    echo "Use a different OUT path, or delete $WORK." >&2
    exit 1
fi
echo "$PARAMS" > "$STAMP"

# A read loop rather than `mapfile`, which is bash 4 and absent on macOS's bash 3.2 — this script
# is meant to be testable off the machine it runs on.
all_libraries=()
while IFS= read -r line; do
    [[ -n "${line// /}" ]] && all_libraries+=("$line")
done < "$LIST"
if ((LIMIT > 0)) && ((LIMIT < ${#all_libraries[@]})); then
    libraries=("${all_libraries[@]:0:LIMIT}")
else
    libraries=("${all_libraries[@]}")
fi

echo "fitting ${#libraries[@]} of ${#all_libraries[@]} library/libraries over ${REGIONS:-$CONTIGS}, $JOBS at a time" >&2
echo "  work directory: $WORK (delete it to force a re-run)" >&2

export RATE_WORK="$WORK" RATE_REF="$REF" RATE_BIN="$BIN" RATE_WALK_ARGS="$WALK_ARGS"
export RATE_MIN_LOCI="$MIN_LOCI" RATE_MAX_REPEATS="$MAX_REPEATS"
export -f has_rows
SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")

printf '%s\0' "${libraries[@]}" \
    | xargs -0 -n 1 -P "$JOBS" "$SELF" --one \
    || echo "  (some libraries failed; their reasons are in $WORK/*.status)" >&2

# ---------------------------------------------------------------------------
# What came back
# ---------------------------------------------------------------------------
list_results() {
    find "$WORK" -maxdepth 1 -name '*.tsv' -size +0 | sort | while read -r f; do
        if has_rows "$f"; then echo "$f"; fi
    done
}
results=$(list_results)
count=$(echo "$results" | grep -c . || true)
((count > 0)) || { echo "no library could be fitted; nothing written" >&2; exit 1; }

# **The timing summary is the point of the first stage.** It is printed before the merge so that a
# run launched at LIMIT=1 answers its question without anyone reading a table.
echo >&2
echo "$results" | while read -r f; do awk -F'\t' '$1=="#elapsed_seconds"{print $2}' "$f"; done \
    | sort -n | awk -v jobs="$JOBS" -v total="${#all_libraries[@]}" '
        {v[NR]=$1; s+=$1}
        END {
          if (NR == 0) exit
          printf "elapsed per library: min %ds, median %ds, max %ds, mean %.0fs (n=%d)\n",
                 v[1], v[int((NR+1)/2)], v[NR], s/NR, NR > "/dev/stderr"
          printf "at %d at a time, all %d libraries would take about %.1f hours\n",
                 jobs, total, (s/NR) * total / jobs / 3600 > "/dev/stderr"
        }'
echo >&2

{
    printf '#survey\tcontigs=%s\tmin_loci=%s\tmax_repeats=%s\tattempted=%s\tfitted=%s\n' \
        "${REGIONS:-$CONTIGS}" "$MIN_LOCI" "$MAX_REPEATS" "${#libraries[@]}" "$count"
    for s in "$WORK"/*.status; do
        [[ -f "$s" ]] || continue
        printf '#missing\t%s\t%s\n' "$(basename "$s" .status)" "$(cat "$s")"
    done
    printf '#rg_columns\tlibrary_key\trg_id\tsample\tlibrary\tfile\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '$1=="#rg"{n=split($6,p,"/"); $2=p[n]"::"$3; print}' "$f"
    done
    printf 'library_key\tperiod\trepeats\tloci\treads\tentries\tslip_rate\tgain_share\tstep_decay\tstart_spread\tguard_share\theterozygosity\tidentified\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '
            $1 == "#rg" { n = split($6, p, "/"); stable[$2] = p[n] "::" $3; next }
            /^#/ { next }
            $1 == "read_group" { next }
            { $1 = stable[$1]; print }
        ' "$f"
    done
} > "$OUT"

rows=$(grep -vc '^#' "$OUT" || true)
echo "wrote $OUT — $count library/libraries, $((rows - 1)) fitted strata" >&2
echo "re-run with a larger LIMIT to add more; finished libraries are not refitted." >&2
