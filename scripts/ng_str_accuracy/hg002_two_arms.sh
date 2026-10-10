#!/usr/bin/env bash
# One binary run two ways on GIAB HG002's tandem-repeat benchmark, at 30x and 50x, scored for
# genotype accuracy and per-tract verdict flips — how the edge-deletion rule of 2026-10-10 was
# measured (doc/devel/reports/implementations/str_accuracy_2026-10-10.md §3).
#
#   DEV_EXTRA_MOUNT=<main checkout>/benchmarks scripts/dev.sh \
#       scripts/ng_str_accuracy/hg002_two_arms.sh <out-dir> <pvc-binary> <ENV_SWITCH>
#
# Arm `base` runs with ENV_SWITCH unset, arm `experiment` with ENV_SWITCH=1 — an experiment build
# puts its change behind an environment variable so one binary is both arms. Writes into
# <out-dir>: each arm's VCF, `genotype.tsv` (right genotypes by period class and depth, the
# number to compare), `sweep.tsv` (allele-exact precision and recall by QUAL), per-arm verdicts,
# and `flips_<depth>.txt` (tracts whose verdict changed, and which way). The tract ground,
# `tier.bed`, is typed with this checkout's `ng_typed_region_dump`, built first if missing.
#
# Env: BENCHMARKS (default: the main checkout's benchmarks/, which holds the data; mount it).
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
out=${1:?usage: $0 <out-dir> <binary> <ENV_SWITCH>}; bin=${2:?binary}; switch=${3:?environment switch}
data=${BENCHMARKS:-$(cd "$root/.." && pwd)/pop_var_caller/benchmarks}
ref=$data/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
truth=$data/ssr_hg002/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz
mkdir -p "$out"
if [[ ! -s $out/tier.bed ]]; then
    dump=$root/target-container/release/examples/ng_typed_region_dump
    [[ -x $dump ]] || (cd "$root" && cargo build --release --example ng_typed_region_dump)
    sort -k1,1 -k2,2n "$data/ssr_hg002/regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed" | cut -f1-3 > "$out/tier_sorted.bed"
    NG_REFERENCE_CHECK=skip "$dump" "$ref" "$out/tier_sorted.bed" calling 2>/dev/null \
        | awk -F'\t' '!/^#/ && $4=="ssr_locus" { print $1"\t"($2-1)"\t"$3"\t"$6 }' > "$out/tier.bed"
fi
for cov in 30x 50x; do
    bam=$data/ssr_hg002/bam/$cov/HG002_TR_v1.0.1_Tier_${cov}.bam
    for arm in base experiment; do
        vcf=$out/${arm}_${cov}.vcf
        if [[ ! -s $vcf ]]; then
            if [[ $arm == experiment ]]; then export "$switch"=1; else unset "$switch"; fi
            NG_REFERENCE_CHECK=skip "$bin" call-from-alignments --reference "$ref" --catalog "$ref.repeats.parquet" \
                --alignment "$bam" --regions "$out/tier_sorted.bed" --output "$vcf" --defaults > "$out/${arm}_${cov}.log" 2>&1
            unset "$switch"
        fi
        python3 "$root/benchmarks/lib/tract_qual_experiment.py" --reference "$ref" --truth "$truth" --query "$vcf" \
            --confident-bed "$out/tier_sorted.bed" --tract-bed "$out/tier.bed" \
            --arm $arm --ground tandem_repeat_tier --depth $cov --sample HG002 --genotype-sample HG002 \
            --calibration-out "$out/calibration.tsv" --sweep-out "$out/sweep.tsv" --genotype-out "$out/genotype.tsv" \
            --verdicts-out "$out/verdicts_${arm}_${cov}.tsv" > "$out/score_${arm}_${cov}.log" 2>&1
        echo "$cov $arm scored"
    done
    python3 "$root/benchmarks/lib/tract_verdict_flips.py" --baseline "$out/verdicts_base_${cov}.tsv" \
        --arm "$out/verdicts_experiment_${cov}.tsv" > "$out/flips_${cov}.txt" 2>&1 || true
done
echo "done: $out"
