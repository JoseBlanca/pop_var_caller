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
#   ./ng_str_library_survey.sh OUT.tsv REF /media/tomato25_bams/crams/*/*.cram
#
# Example on rick:
#   ./scripts/ng_str_library_survey.sh ~/tmp/stutter_by_library.tsv \
#       /home/joxi/refs/S_lycopersicum_chromosomes.4.00.fa \
#       /media/tomato25_bams/crams/*/*.cram
#
# Knobs, all env-overridable:
#   CONTIGS=SL4.0ch01   the walk. ~90 Mb gives ~200k loci, which settles these curves; the cost
#                       scales with the number of libraries, so a wider walk buys precision you do
#                       not need. Set to "" to walk the whole genome.
#   MIN_COPIES=2        type from this many copies at every period. **Not optional for this
#                       question**: at ng's defaults region typing emits nothing below
#                       [6,4,4,3,3,3], so every curve would start exactly where the floor is meant
#                       to be decided and the measurement would be censored at the wrong place.
#   BATCH=20            CRAMs per invocation. Each one holds its files open for the whole walk, so
#                       this bounds open handles and memory rather than changing any answer — the
#                       loci walked are identical across batches given the same reference, contigs
#                       and MIN_COPIES, so the pieces join exactly.

set -euo pipefail

OUT=${1:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
REF=${2:?"usage: $0 OUT.tsv REF CRAM [CRAM ...]"}
shift 2
(($# > 0)) || { echo "no CRAMs given" >&2; exit 1; }

# `${CONTIGS-...}` without the colon, so an explicitly empty CONTIGS means "the whole genome"
# rather than falling back to the default — which is what the usage above promises.
CONTIGS=${CONTIGS-SL4.0ch01}
MIN_COPIES=${MIN_COPIES:-2}
BATCH=${BATCH:-20}

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BIN="$REPO_ROOT/target/release/examples/ng_str_stutter_by_library"

# Fail on the setup rather than half way through a six-hour walk.
[[ -x "$BIN" ]] || {
    echo "not built: $BIN" >&2
    echo "  cargo build --release --example ng_str_stutter_by_library" >&2
    exit 1
}
[[ -f "$REF" ]] || { echo "reference not found: $REF" >&2; exit 1; }
for cram in "$@"; do
    [[ -f "$cram" ]] || { echo "not a file: $cram" >&2; exit 1; }
done

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

contig_args=()
[[ -n "$CONTIGS" ]] && contig_args=(--contigs "$CONTIGS")

total=$#
echo "surveying $total library file(s) over ${CONTIGS:-the whole genome}, from $MIN_COPIES copies" >&2

batch=0
while (($# > 0)); do
    batch=$((batch + 1))
    chunk=("${@:1:BATCH}")
    if (($# > BATCH)); then shift "$BATCH"; else shift $#; fi
    echo "  batch $batch: ${#chunk[@]} file(s)" >&2
    "$BIN" "${contig_args[@]}" --min-copies "$MIN_COPIES" "$REF" "${chunk[@]}" \
        > "$WORK/batch.$batch.tsv"
done

# Merge. **The numeric read_group is minted per run**, so a plain concatenation would collide two
# batches' group 0. `(file, rg_id)` is the stable identity — the SAM specification makes `@RG ID`
# unique within its file — so each batch's rows are re-keyed onto it before they are joined.
echo "  merging $batch batch(es)" >&2
{
    printf '#survey\tcontigs=%s\tmin_copies=%s\tfiles=%s\tbatches=%s\n' \
        "${CONTIGS:-ALL}" "$MIN_COPIES" "$total" "$batch"
    printf '#rg_columns\tlibrary_key\trg_id\tsample\tlibrary\tlibrary_origin\texperiment\texperiment_origin\tplatform\tloci\treads\tmean_tract_bases_per_read\tfile\n'
    printf '#floor_columns\tlibrary_key\tperiod\timplied_floor\tcriterion\n'
    printf 'library_key\tperiod\trepeats\tloci\treads\toff_ref_reads\toff_ref_share\tnot_whole_reads\tguard_share\tend_bucket_reads\n'

    for f in "$WORK"/batch.*.tsv; do
        awk -F'\t' -v OFS='\t' '
            # Build numeric id -> stable key from this batch own #rg table.
            $1 == "#rg" {
                # $2 numeric id, $3 rg_id, $13 file
                n = split($13, path, "/")
                key = path[n] "::" $3
                stable[$2] = key
                $1 = "#rg"; $2 = key
                print
                next
            }
            $1 == "#floor" { $2 = stable[$2]; print; next }
            # Drop each batch own column headers; the merged ones are written once above.
            /^#/ { next }
            $1 == "library_key" || $1 == "read_group" { next }
            { $1 = stable[$1]; print }
        ' "$f"
    done
} > "$OUT"

rows=$(grep -vc '^#' "$OUT" || true)
libraries=$(grep -c '^#rg\b' "$OUT" || true)
echo "wrote $OUT — $libraries library/libraries, $rows stratum rows" >&2
echo >&2
echo "The three blocks in it:" >&2
echo "  #rg    one line per read group: which library, and how much data it contributed" >&2
echo "  #floor the copy floor each period's own data implies, per library — the answer" >&2
echo "  rows   the per-(library, period, repeat count) curves the floors were read off" >&2
echo >&2
echo "Before comparing libraries, join to rick_sample_manifest.sh on the \`@RG\` id: read length" >&2
echo "is a confound this tool cannot see, and mixing 100 bp with 150 bp libraries makes stutter" >&2
echo "appear to start at different tract lengths for purely geometric reasons." >&2
