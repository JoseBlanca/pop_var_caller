#!/usr/bin/env bash
# Survey where tracts start to stutter, across many libraries — run this ON RICK.
#
# ng routes a locus to the STR path when it is **likely to stutter**, not merely when it contains a
# repeat (doc/devel/ng/spec/parameter_prepass_ssr.md §5.1), and the per-period copy floors are where
# that line is drawn. Whether a tract stutters is a property of the **library** — PCR amplification
# stutters more than a PCR-free preparation — and every number ng has today comes from a single
# library per species. This runs the measurement across the archive so the floors rest on the axis
# that actually drives them.
#
# Usage:
#   ./ng_str_library_survey.sh OUT.tsv REF CRAM [CRAM ...]
#
# Example on rick:
#   ./scripts/ng_str_library_survey.sh ~/tmp/stutter_by_library.tsv \
#       /home/joxi/refs/S_lycopersicum_chromosomes.4.00.fa \
#       /media/tomato25_bams/crams/*/*.cram
#
# ## It is built for an archive that is being written to while it runs
#
# The walk takes hours, so a file that was complete when the run started can be appended to before
# the walk reaches it, and one that was mid-write at the start can be finished long before the run
# ends. **A gate at the front cannot express either.** So:
#
#   * **one walk per file**, not one per batch — a file's failure costs that file and nothing else,
#     and the region-typing pass it repeats is a fraction of the read fetching it does not;
#   * **each file is checked immediately before its walk and again immediately after**, and a file
#     whose size or mtime moved in between has its result thrown away. This is the check that
#     matters: a half-written CRAM does not reliably *error*, it can simply run out of records part
#     way through the region, and that library would come back with thinner strata and a copy floor
#     read off them with nothing saying why;
#   * **files that were not ready are retried in later passes**, so a library still being written at
#     the start is picked up once it settles rather than lost;
#   * **results are kept in a work directory as they are produced**, so re-running the same command
#     resumes instead of starting over — which is also how you add libraries to a finished survey.
#
# Knobs, all env-overridable:
#   CONTIGS=SL4.0ch01   the walk. ~90 Mb gives ~200k loci, which settles these curves; the cost
#                       scales with the number of libraries, so a wider walk buys precision you do
#                       not need. Set to "" to walk the whole genome.
#   MIN_COPIES=2        type from this many copies at every period. **Not optional for this
#                       question**: at ng's defaults region typing emits nothing below
#                       [6,4,4,3,3,3], so every curve would start exactly where the floor is meant
#                       to be decided and the measurement would be censored at the wrong place.
#   JOBS=8              files walked at once. **The knob that decides whether this finishes.** One
#                       file over 90 Mb takes ~10 minutes, so 2,475 of them sequentially is 17 days
#                       and at 32-way is half a day. Each walk is independent — that is what the
#                       per-file design buys — so set this to about the core count, backing off if
#                       the archive's disks are the bottleneck rather than the CPU.
#   REGIONS=            a BED to walk instead of whole contigs. **The other way to make a large
#                       survey finish**: these curves need thousands of loci per stratum, not the
#                       ~200k that a whole chromosome gives, so a 10 Mb slice is ~9x faster and
#                       still ample. Overrides CONTIGS.
#   PASSES=3            how many times to come back for files that were not ready.
#   WAIT=900            seconds between passes. Long enough that a file being written has plausibly
#                       finished; the survey is hours long, so this costs nothing.
#   WORK=OUT.work       where per-file results accumulate. Delete it to force a full re-run.

set -euo pipefail

# ---------------------------------------------------------------------------
# The per-file worker, re-entered through `--one`
# ---------------------------------------------------------------------------
#
# Walking one file is its own invocation so the parent can drive many at once with `xargs -P`. The
# worker reports back through `$WORK/<key>.status` rather than stdout, because several of them are
# writing at the same time and interleaved lines would be unreadable and unparseable.
if [[ "${1:-}" == "--one" ]]; then
    cram=$2
    key=$(basename "$cram" | tr -c 'A-Za-z0-9._-' '_')
    result="$SURVEY_WORK/$key.tsv"
    status="$SURVEY_WORK/$key.status"
    stamp_of() { stat -c '%s %Y' "$1" 2>/dev/null || stat -f '%z %m' "$1" 2>/dev/null || echo gone; }
    index_of() {
        local candidate
        for candidate in "$cram.crai" "$cram.bai" "${cram%.*}.crai" "${cram%.*}.bai"; do
            [[ -f "$candidate" ]] && { echo "$candidate"; return 0; }
        done
        return 1
    }
    [[ -s "$result" ]] && exit 0
    if [[ ! -f "$cram" ]]; then echo "missing" > "$status"; exit 0; fi
    if ! index=$(index_of); then echo "no index yet" > "$status"; exit 0; fi
    if [[ "$cram" -nt "$index" ]]; then
        echo "index older than the file — still growing, or not re-indexed" > "$status"
        exit 0
    fi
    if [[ "$SURVEY_QUICKCHECK" == "1" ]] && ! samtools quickcheck "$cram" 2>/dev/null; then
        echo "no EOF block — still being written" > "$status"
        exit 0
    fi
    before=$(stamp_of "$cram")
    read -r -a survey_args <<< "$SURVEY_ARGS"
    if ! "$SURVEY_BIN" "${survey_args[@]}" "$SURVEY_REF" "$cram" \
            > "$result.partial" 2> "$SURVEY_WORK/$key.err"; then
        rm -f "$result.partial"
        echo "walk failed: $(tail -n 1 "$SURVEY_WORK/$key.err" 2>/dev/null | tr '\t' ' ')" > "$status"
        exit 0
    fi
    # **The check that catches the silent case.** A file appended to during its own walk yields a
    # result that looks like a thin library rather than an error, so the stamp — not the exit code —
    # decides whether the result is trustworthy.
    if [[ "$before" != "$(stamp_of "$cram")" ]]; then
        rm -f "$result.partial"
        echo "changed during its own walk — result discarded" > "$status"
        exit 0
    fi
    mv "$result.partial" "$result"
    rm -f "$status"
    exit 0
fi

OUT=${1:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
REF=${2:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
shift 2
(($# > 0)) || { echo "no CRAMs given" >&2; exit 1; }

# `${CONTIGS-...}` without the colon, so an explicitly empty CONTIGS means "the whole genome"
# rather than falling back to the default — which is what the usage above promises.
CONTIGS=${CONTIGS-SL4.0ch01}
REGIONS=${REGIONS:-}
MIN_COPIES=${MIN_COPIES:-2}
PASSES=${PASSES:-3}
WAIT=${WAIT:-900}
WORK=${WORK:-$OUT.work}
JOBS=${JOBS:-8}

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BIN="$REPO_ROOT/target/release/examples/ng_str_stutter_by_library"

[[ -x "$BIN" ]] || {
    echo "not built: $BIN" >&2
    echo "  cargo build --release --example ng_str_stutter_by_library" >&2
    exit 1
}
[[ -f "$REF" ]] || { echo "reference not found: $REF" >&2; exit 1; }

mkdir -p "$WORK"

# ---------------------------------------------------------------------------
# The work directory belongs to one set of walk parameters
# ---------------------------------------------------------------------------
#
# **Resuming with a different walk would silently mix incomparable results.** The per-file results
# are keyed by filename alone, so a run over `SL4.0ch01` and a later one over a 10 Mb BED would land
# side by side in the same directory and merge into one table whose libraries were measured over
# different loci — a survey that looks whole and compares nothing. Since narrowing the walk is the
# obvious thing to reach for when a run is taking too long, this is a mistake worth making
# impossible rather than documenting.
PARAMS="contigs=${CONTIGS:-} regions=${REGIONS:-} min_copies=$MIN_COPIES"
STAMP="$WORK/.survey-params"
existing=$(find "$WORK" -maxdepth 1 -name '*.tsv' -size +0 2>/dev/null | wc -l | tr -d ' ')
if [[ -f "$STAMP" ]]; then
    recorded=$(cat "$STAMP")
    if [[ "$recorded" != "$PARAMS" ]]; then
        echo "this work directory was built with different walk parameters:" >&2
        echo "  it holds: $recorded" >&2
        echo "  you asked: $PARAMS" >&2
        echo "Results from two different walks are not comparable. Use a different OUT path, or" >&2
        echo "delete $WORK to start over." >&2
        exit 1
    fi
elif ((existing > 0)); then
    # Written by a version of this script that did not stamp, so the parameters are unknown.
    echo "$WORK holds $existing result(s) but no record of how they were walked." >&2
    echo "If they came from the same walk you are asking for now, re-run with" >&2
    echo "  SURVEY_ADOPT_EXISTING=1" >&2
    echo "Otherwise delete $WORK. Mixing two walks makes a table that compares nothing." >&2
    [[ "${SURVEY_ADOPT_EXISTING:-0}" == "1" ]] || exit 1
    echo "$PARAMS" > "$STAMP"
else
    echo "$PARAMS" > "$STAMP"
fi

walk_args=()
if [[ -n "$REGIONS" ]]; then
    [[ -f "$REGIONS" ]] || { echo "regions BED not found: $REGIONS" >&2; exit 1; }
    walk_args=(--regions "$REGIONS")
elif [[ -n "$CONTIGS" ]]; then
    walk_args=(--contigs "$CONTIGS")
fi
walk_args+=(--min-copies "$MIN_COPIES")

command -v samtools >/dev/null && HAVE_QUICKCHECK=1 || HAVE_QUICKCHECK=0
((HAVE_QUICKCHECK)) || echo "warning: samtools not on PATH — no EOF check, only indexes and stamps" >&2

key_of() { basename "$1" | tr -c 'A-Za-z0-9._-' '_'; }

export SURVEY_WORK="$WORK" SURVEY_REF="$REF" SURVEY_BIN="$BIN"
export SURVEY_QUICKCHECK="$HAVE_QUICKCHECK"
export SURVEY_ARGS="${walk_args[*]}"
SELF=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")

declare -a pending=("$@")
total=${#pending[@]}
echo "surveying $total file(s) over ${REGIONS:-${CONTIGS:-the whole genome}}, from $MIN_COPIES copies, $JOBS at a time" >&2
echo "  work directory: $WORK (delete it to force a full re-run)" >&2

declare -a deferred_files=() deferred_reasons=()

for ((pass = 1; pass <= PASSES; pass++)); do
    ((${#pending[@]} > 0)) || break
    ((pass == 1)) || {
        echo "pass $pass: waiting ${WAIT}s for ${#pending[@]} unfinished file(s) to settle" >&2
        sleep "$WAIT"
    }
    # Clear last pass's verdicts so a file that has since become readable is not still carrying one.
    for cram in "${pending[@]}"; do rm -f "$WORK/$(key_of "$cram").status"; done

    echo "pass $pass: ${#pending[@]} file(s) to walk" >&2
    # **`xargs -P` rather than a hand-rolled job pool**, because it is the portable way to keep N
    # slots full: a `wait`-on-every-N loop stalls the whole group on its slowest member, and these
    # files differ several-fold in size.
    #
    # **The file is appended as the last argument rather than substituted with `-I`.** `-I` implies
    # one argument per run, so passing `-n 1` beside it made GNU xargs warn that the two are
    # mutually exclusive — harmless, since `-I` won and one file per worker is what was wanted, but
    # `-I` also has a history of interacting badly with `-P` across implementations. Appending needs
    # neither flag and the worker's signature (`--one <cram>`) already puts the file last.
    printf '%s\0' "${pending[@]}" \
        | xargs -0 -n 1 -P "$JOBS" "$SELF" --one \
        || echo "  (some workers reported failures; their reasons are in the work directory)" >&2

    declare -a still_pending=()
    deferred_files=()
    deferred_reasons=()
    done_count=0
    for cram in "${pending[@]}"; do
        key=$(key_of "$cram")
        if [[ -s "$WORK/$key.tsv" ]]; then
            done_count=$((done_count + 1))
            continue
        fi
        still_pending+=("$cram")
        deferred_files+=("$cram")
        deferred_reasons+=("$(cat "$WORK/$key.status" 2>/dev/null || echo 'no result and no reason recorded')")
    done
    echo "  pass $pass done: $done_count walked, ${#still_pending[@]} outstanding" >&2
    pending=("${still_pending[@]:-}")
    # `${arr[@]:-}` yields one empty element on an empty array under `set -u`; drop it.
    ((${#pending[@]} == 1)) && [[ -z "${pending[0]}" ]] && pending=()
done

surveyed=$(find "$WORK" -maxdepth 1 -name '*.tsv' -size +0 | wc -l | tr -d ' ')
((surveyed > 0)) || { echo "no file could be surveyed; nothing written" >&2; exit 1; }

# Merge. **The numeric read_group is minted per invocation**, so rows are re-keyed onto
# `(file, rg_id)` — the stable identity, since the SAM specification makes `@RG ID` unique within
# its file — before the per-file results are joined.
echo "merging $surveyed result(s)" >&2
{
    # **The output says what it does not contain.** A survey missing a third of the archive because
    # those CRAMs were mid-write looks exactly like a survey of an archive that size, so the
    # exclusions travel with the numbers rather than in a terminal scrollback.
    printf '#survey\tcontigs=%s\tmin_copies=%s\tgiven=%s\tsurveyed=%s\tmissing=%s\tpasses=%s\n' \
        "${CONTIGS:-ALL}" "$MIN_COPIES" "$total" "$surveyed" "${#deferred_files[@]}" "$PASSES"
    for ((i = 0; i < ${#deferred_files[@]}; i++)); do
        printf '#missing\t%s\t%s\n' "${deferred_files[$i]}" "${deferred_reasons[$i]}"
    done
    # Three passes so the file reads top to bottom — every `#rg`, then every `#floor`, then the
    # rows — rather than interleaving one library's three blocks with the next library's.
    # The results are one small file per library, so re-reading them twice more costs nothing.
    results=$(find "$WORK" -maxdepth 1 -name '*.tsv' -size +0 | sort)

    printf '#rg_columns\tlibrary_key\trg_id\tsample\tlibrary\tlibrary_origin\texperiment\texperiment_origin\tplatform\tloci\treads\tmean_tract_bases_per_read\tfile\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '
            $1 == "#rg" { n = split($13, path, "/"); $2 = path[n] "::" $3; print }
        ' "$f"
    done

    printf '#floor_columns\tlibrary_key\tperiod\timplied_floor\tcriterion\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '
            $1 == "#rg" { n = split($13, path, "/"); stable[$2] = path[n] "::" $3; next }
            $1 == "#floor" { $2 = stable[$2]; print }
        ' "$f"
    done

    printf 'library_key\tperiod\trepeats\tloci\treads\toff_ref_reads\toff_ref_share\tnot_whole_reads\tguard_share\tend_bucket_reads\n'
    echo "$results" | while read -r f; do
        awk -F'\t' -v OFS='\t' '
            $1 == "#rg" { n = split($13, path, "/"); stable[$2] = path[n] "::" $3; next }
            /^#/ { next }
            $1 == "library_key" || $1 == "read_group" { next }
            { $1 = stable[$1]; print }
        ' "$f"
    done
} > "$OUT"

libraries=$(grep -c '^#rg	' "$OUT" || true)
# Exclude the column-header line, which is not a comment and is not a stratum either.
rows=$(grep -v '^#' "$OUT" | grep -vc '^library_key	' || true)
echo "wrote $OUT — $libraries library/libraries, $rows stratum rows" >&2
if ((${#deferred_files[@]} > 0)); then
    echo "${#deferred_files[@]} file(s) never became readable and are listed as #missing in the output:" >&2
    for ((i = 0; i < ${#deferred_files[@]}; i++)); do
        echo "  $(basename "${deferred_files[$i]}") — ${deferred_reasons[$i]}" >&2
    done
    echo "Re-run the same command later to pick them up; finished files are not re-walked." >&2
fi
echo >&2
echo "The blocks in it:" >&2
echo "  #survey/#missing  what was walked and what was not, with the reason" >&2
echo "  #rg               one line per read group: which library, and how much data it gave" >&2
echo "  #floor            the copy floor each period's own data implies, per library — the answer" >&2
echo "  rows              the per-(library, period, repeat count) curves behind the floors" >&2
echo >&2
echo "Before comparing libraries, join to rick_sample_manifest.sh on the \`@RG\` id: read length" >&2
echo "is a confound this tool cannot see, and mixing 100 bp with 150 bp libraries makes stutter" >&2
echo "appear to start at different tract lengths for purely geometric reasons." >&2
