#!/usr/bin/env bash
# Why does ng miss the truth variants it misses — is it genotyping them wrong,
# or is it never looking?
#
#   benchmarks/giab/src/ng_missed_sites_probe.sh [COVERAGE] [SAMPLE ...]
#
# ng does not build loci inside tandem repeats yet, so a truth variant inside a
# tract is out of its reach by construction and scores as a false negative in
# the accuracy dashboard exactly like a variant it looked at and got wrong.
# Those two are worth very different amounts, and the accuracy table cannot
# tell them apart. This does.
#
# For each (sample, class) it takes the truth variants ng emitted nothing for,
# writes them as a one-base-per-site BED, and runs ng again over just those
# bases. ng's own run report then says how many loci it built there — so
#
#   loci_built ≈ 0   the sites are ground ng does not call; the recall gap is
#                    the unbuilt repeat-tract path, not the genotyper.
#   loci_built ≈ N   ng looked at the sites and did not call them; that IS the
#                    genotyper, and worth chasing.
#
# The production caller's missed sites go through the same probe, and note what
# that row means: the locus count is always ng's, so for the production caller
# it says how much of ITS residual miss list also lies on ground ng cannot
# reach. That is the comparison worth having — whether the two callers are
# failing on the same ground or on different ground.
#
# Output: results/per_sample/ng_missed_sites.tsv, one row per
# (coverage, caller, sample, class), read by freebayes_comparison_dashboard.py.
#
# Args:
#   COVERAGE   bam/ coverage subdir (default: 300x — the depth where a miss is
#              least likely to be sampling and most likely to be the caller)
#   SAMPLE...  subset of {HG002,HG003,HG004} (default: all three)
#
# Env overrides:
#   NG_BIN, NG_CATALOG   as in run_ng_per_sample.sh
#   CALLERS              space-separated result subdirs to probe
#                        (default: "ng high-recall")

set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"

COVERAGE="${1:-300x}"
shift || true
SAMPLES=("$@")
if (( ${#SAMPLES[@]} == 0 )); then
    SAMPLES=(HG002 HG003 HG004)
fi
CALLERS="${CALLERS:-ng high-recall}"

REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BAM_DIR="$BENCH_DIR/per_sample/bam/$COVERAGE"
BED_DIR="$BENCH_DIR/per_sample/bed"
TRUTH_DIR="$BENCH_DIR/per_sample/vcf"
RESULTS_DIR="$BENCH_DIR/results/per_sample"
NG_CATALOG="${NG_CATALOG:-${REFERENCE}.repeats.parquet}"
TSV_OUT="${TSV_OUT:-$RESULTS_DIR/ng_missed_sites.tsv}"

command -v bcftools >/dev/null || { echo "bcftools not on PATH" >&2; exit 1; }

discover_bin() {
    if [[ -z "${NG_BIN:-}" ]]; then
        local candidate
        for candidate in \
            "$PROJECT_ROOT/target-container/release/pop_var_caller_exp" \
            "$PROJECT_ROOT/target/release/pop_var_caller_exp"; do
            if [[ -x "$candidate" ]] && "$candidate" --version >/dev/null 2>&1; then
                NG_BIN="$candidate"; break
            fi
        done
    fi
    [[ -n "${NG_BIN:-}" && -x "${NG_BIN}" ]] || {
        echo "no pop_var_caller_exp binary found; set NG_BIN=<path>" >&2; exit 1; }
}

align_file() {
    if [[ "$COVERAGE" != "300x" ]]; then
        echo "$BAM_DIR/${1}.${COVERAGE}.seed42.bam"; return
    fi
    case "$1" in
        HG002) echo "$BAM_DIR/HG002_reads_selected_100_rg.cram" ;;
        *)     echo "$BAM_DIR/${1}_bench_azar_merged_100.sorted.m5.bam" ;;
    esac
}

discover_bin
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# Same normalisation as the accuracy dashboard, so "missed" means the same
# thing here as it does in the table this explains: BED-restricted,
# left-aligned, biallelic-split, class-filtered, matched on POS+REF+ALT.
normalize() {
    local src="$1" cls="$2" bed="$3" out="$4"; shift 4
    bcftools view "$@" -T "$bed" "$src" -Ou \
      | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
      | bcftools view -v "$cls" -Ou \
      | bcftools sort -Oz -o "$out" 2>/dev/null
    bcftools index -f -t "$out"
}

printf 'coverage\tcaller\tsample\tclass\tmissed_sites\tloci_built\tcalled_bp\ttract_bp\n' > "$TSV_OUT"

echo "binary    : $NG_BIN"
echo "dataset   : per_sample / $COVERAGE"
echo "callers   : $CALLERS"
echo "samples   : ${SAMPLES[*]}"
echo "tsv       : $TSV_OUT"
echo
printf '%-14s %-8s %-7s %13s %11s\n' caller sample class missed_sites loci_built
printf '%s\n' "----------------------------------------------------------------"

for caller in $CALLERS; do
    for sample in "${SAMPLES[@]}"; do
        bed="$BED_DIR/${sample}_bench_azar_merged_100.bed"
        truth="$TRUTH_DIR/${sample}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz"
        query="$RESULTS_DIR/$COVERAGE/$caller/${sample}.vcf"
        aln="$(align_file "$sample")"
        [[ -f "$query" && -f "$aln" ]] || { echo "skip $caller/$sample (missing input)" >&2; continue; }

        for cls in snps indels; do
            normalize "$truth" "$cls" "$bed" "$TMP/t.vcf.gz" -f PASS
            normalize "$query" "$cls" "$bed" "$TMP/q.vcf.gz"
            rm -rf "$TMP/isec"
            bcftools isec -p "$TMP/isec" "$TMP/t.vcf.gz" "$TMP/q.vcf.gz" >/dev/null

            # 0000 = truth-only = the sites this caller emitted nothing for.
            # One BED line per distinct position; a position carrying two ALTs
            # is one place to look, not two.
            bcftools query -f '%CHROM\t%POS\n' "$TMP/isec/0000.vcf" \
              | sort -u -k1,1 -k2,2n \
              | awk '{ printf "%s\t%d\t%d\n", $1, $2-1, $2 }' > "$TMP/missed.bed"
            missed=$(wc -l < "$TMP/missed.bed" | tr -d ' ')
            if (( missed == 0 )); then
                printf '%s\t%s\t%s\t%s\t0\t0\t0\t0\n' \
                    "$COVERAGE" "$caller" "$sample" "$cls" >> "$TSV_OUT"
                printf '%-14s %-8s %-7s %13d %11s\n' "$caller" "$sample" "$cls" 0 "-"
                continue
            fi

            # Ask ng to call over exactly those bases and read its run report.
            # The walk covers a typed region whole, so the base counts below
            # are of the expanded regions, not of the requested bases; the
            # locus count is the number that answers the question.
            report="$TMP/report.txt"
            "$NG_BIN" call-from-alignments \
                --reference "$REFERENCE" \
                --catalog "$NG_CATALOG" \
                --alignment "$aln" \
                --regions "$TMP/missed.bed" \
                --output "$TMP/probe.vcf" \
                --defaults > "$report" 2>&1 || true

            loci=$(awk '/^loci called:/ { print $3; exit }' "$report")
            called_bp=$(awk '/^  called:/ { print $2; exit }' "$report")
            tract_bp=$(awk '/repeat tracts this caller has not built yet:/ { print $(NF-2); exit }' "$report")
            loci="${loci:-0}"; called_bp="${called_bp:-0}"; tract_bp="${tract_bp:-0}"

            printf '%s\t%s\t%s\t%s\t%d\t%s\t%s\t%s\n' \
                "$COVERAGE" "$caller" "$sample" "$cls" \
                "$missed" "$loci" "$called_bp" "$tract_bp" >> "$TSV_OUT"
            printf '%-14s %-8s %-7s %13d %11s\n' "$caller" "$sample" "$cls" "$missed" "$loci"
        done
    done
done

echo
echo "tsv: $TSV_OUT"
