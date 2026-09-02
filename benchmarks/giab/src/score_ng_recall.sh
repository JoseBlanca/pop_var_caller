#!/usr/bin/env bash
# Precision and recall for one directory of ng calls, per sample and class.
#
#   benchmarks/giab/src/score_ng_recall.sh <coverage> <results subdir>
#
# The method is `accuracy_dashboard.py`'s, in shell: restrict both truth and
# query to the sample's own confident BED, left-align and split multi-allelics,
# filter to one class, and intersect on POS+REF+ALT. TP and FN are counted on
# the truth side, FP on the query side. The truth set is FILTER PASS; the query
# keeps its own FILTER, as the dashboard does.
#
# Two directories under `results/per_sample/<coverage>/` can be scored the same
# way and compared — which is what a routing change needs.
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"
COVERAGE="${1:?coverage, e.g. 30x}"
SUBDIR="${2:?results subdirectory, e.g. ng}"

REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BED_DIR="$BENCH_DIR/per_sample/bed"
TRUTH_DIR="$BENCH_DIR/per_sample/vcf"
QUERY_DIR="$BENCH_DIR/results/per_sample/$COVERAGE/$SUBDIR"

printf 'sample\tclass\ttp\tfp\tfn\trecall\tprecision\n'
for sample in HG002 HG003 HG004; do
    bed="$BED_DIR/${sample}_bench_azar_merged_100.bed"
    truth="$TRUTH_DIR/${sample}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz"
    query="$QUERY_DIR/${sample}.vcf"
    [[ -f "$query" ]] || { echo "missing $query" >&2; exit 1; }
    for cls in snps indels; do
        work="$(mktemp -d)"
        bcftools view -f PASS -T "$bed" "$truth" -Ou \
            | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
            | bcftools view -v "$cls" -Oz -o "$work/truth.vcf.gz"
        bcftools index -t "$work/truth.vcf.gz"
        bcftools view -T "$bed" "$query" -Ou \
            | bcftools norm -f "$REFERENCE" -m -any -Ou 2>/dev/null \
            | bcftools view -v "$cls" -Oz -o "$work/query.vcf.gz"
        bcftools index -t "$work/query.vcf.gz"
        bcftools isec -p "$work/isec" "$work/truth.vcf.gz" "$work/query.vcf.gz" >/dev/null
        fn=$(grep -vc '^#' "$work/isec/0000.vcf" || true)
        fp=$(grep -vc '^#' "$work/isec/0001.vcf" || true)
        tp=$(grep -vc '^#' "$work/isec/0002.vcf" || true)
        recall=$(awk -v t="$tp" -v f="$fn" 'BEGIN{ if (t+f) printf "%.4f", t/(t+f); else print "." }')
        precision=$(awk -v t="$tp" -v f="$fp" 'BEGIN{ if (t+f) printf "%.4f", t/(t+f); else print "." }')
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$sample" "$cls" "$tp" "$fp" "$fn" "$recall" "$precision"
        rm -rf "$work"
    done
done
