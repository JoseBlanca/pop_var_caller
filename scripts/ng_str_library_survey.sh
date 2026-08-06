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
#   PASSES=3            how many times to come back for files that were not ready.
#   WAIT=900            seconds between passes. Long enough that a file being written has plausibly
#                       finished; the survey is hours long, so this costs nothing.
#   WORK=OUT.work       where per-file results accumulate. Delete it to force a full re-run.

set -euo pipefail

OUT=${1:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
REF=${2:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
shift 2
(($# > 0)) || { echo "no CRAMs given" >&2; exit 1; }

# `${CONTIGS-...}` without the colon, so an explicitly empty CONTIGS means "the whole genome"
# rather than falling back to the default — which is what the usage above promises.
CONTIGS=${CONTIGS-SL4.0ch01}
MIN_COPIES=${MIN_COPIES:-2}
PASSES=${PASSES:-3}
WAIT=${WAIT:-900}
WORK=${WORK:-$OUT.work}

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BIN="$REPO_ROOT/target/release/examples/ng_str_stutter_by_library"

[[ -x "$BIN" ]] || {
    echo "not built: $BIN" >&2
    echo "  cargo build --release --example ng_str_stutter_by_library" >&2
    exit 1
}
[[ -f "$REF" ]] || { echo "reference not found: $REF" >&2; exit 1; }

mkdir -p "$WORK"
contig_args=()
[[ -n "$CONTIGS" ]] && contig_args=(--contigs "$CONTIGS")

command -v samtools >/dev/null && HAVE_QUICKCHECK=1 || HAVE_QUICKCHECK=0
((HAVE_QUICKCHECK)) || echo "warning: samtools not on PATH — no EOF check, only indexes and stamps" >&2

# `stat` differs between GNU and BSD and this script is written on one and run on the other.
stamp_of() {
    stat -c '%s %Y' "$1" 2>/dev/null || stat -f '%z %m' "$1" 2>/dev/null || echo "gone"
}

index_of() {
    local cram=$1 candidate
    for candidate in "$cram.crai" "$cram.bai" "${cram%.*}.crai" "${cram%.*}.bai"; do
        [[ -f "$candidate" ]] && { echo "$candidate"; return 0; }
    done
    return 1
}

# Is this file ready to be walked *right now*? Two states an in-progress file can be in, and each
# needs its own test: still being written (no EOF block, which `samtools quickcheck` detects), or
# written but not yet indexed — or indexed and then appended to, which `quickcheck` passes.
ready_reason() {
    local cram=$1 index
    [[ -f "$cram" ]] || { echo "missing"; return 1; }
    index=$(index_of "$cram") || { echo "no index yet"; return 1; }
    [[ "$cram" -nt "$index" ]] && { echo "index older than the file — still growing, or not re-indexed"; return 1; }
    if ((HAVE_QUICKCHECK)) && ! samtools quickcheck "$cram" 2>/dev/null; then
        echo "no EOF block — still being written"
        return 1
    fi
    return 0
}

# A stable key per file, so a resumed run recognises what it already has and the merge can attribute
# rows without depending on the numeric read-group id, which is minted per invocation.
key_of() { basename "$1" | tr -c 'A-Za-z0-9._-' '_'; }

declare -a pending=("$@")
total=${#pending[@]}
echo "surveying $total file(s) over ${CONTIGS:-the whole genome}, from $MIN_COPIES copies" >&2
echo "  work directory: $WORK (delete it to force a full re-run)" >&2

declare -a deferred_files=() deferred_reasons=()

for ((pass = 1; pass <= PASSES; pass++)); do
    ((${#pending[@]} > 0)) || break
    ((pass == 1)) || {
        echo "pass $pass: waiting ${WAIT}s for ${#pending[@]} unfinished file(s) to settle" >&2
        sleep "$WAIT"
    }
    declare -a still_pending=()
    deferred_files=()
    deferred_reasons=()
    index=0
    for cram in "${pending[@]}"; do
        index=$((index + 1))
        key=$(key_of "$cram")
        result="$WORK/$key.tsv"
        # Resume: anything already walked cleanly is left alone, which is what makes re-running the
        # same command add libraries rather than redo them.
        [[ -s "$result" ]] && continue

        if ! reason=$(ready_reason "$cram"); then
            still_pending+=("$cram")
            deferred_files+=("$cram")
            deferred_reasons+=("$reason")
            echo "  [$index/${#pending[@]}] deferring $(basename "$cram") — $reason" >&2
            continue
        fi

        before=$(stamp_of "$cram")
        echo "  [$index/${#pending[@]}] walking $(basename "$cram")" >&2
        if ! "$BIN" "${contig_args[@]}" --min-copies "$MIN_COPIES" "$REF" "$cram" \
                > "$result.partial" 2> "$WORK/$key.err"; then
            rm -f "$result.partial"
            still_pending+=("$cram")
            deferred_files+=("$cram")
            deferred_reasons+=("walk failed: $(tail -n 1 "$WORK/$key.err" 2>/dev/null | tr '\t' ' ')")
            echo "    failed — will retry in a later pass" >&2
            continue
        fi
        # **The check that catches the silent case.** A file appended to during its own walk yields
        # a result that looks like a thin library rather than an error, so the stamp decides whether
        # the result is trustworthy — not whether the command succeeded.
        after=$(stamp_of "$cram")
        if [[ "$before" != "$after" ]]; then
            rm -f "$result.partial"
            still_pending+=("$cram")
            deferred_files+=("$cram")
            deferred_reasons+=("changed during its own walk — result discarded")
            echo "    changed while being walked; discarding and retrying later" >&2
            continue
        fi
        mv "$result.partial" "$result"
    done
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
