#!/usr/bin/env bash
# Two binaries on GIAB HG002's tandem-repeat benchmark, at 30x and 50x, scored for
# genotype accuracy and per-tract verdict flips — how the edge-deletion rule of 2026-10-10 was
# measured (doc/devel/reports/implementations/str_accuracy_2026-10-10.md §3).
#
#   DEV_EXTRA_MOUNT=<main checkout>/benchmarks scripts/dev.sh \
#       scripts/ng_str_accuracy/hg002_two_binaries.sh <out-dir> <base-binary> <experiment-binary>
#
# Arm `base` is the first binary, arm `experiment` the second (2026-10-11; `hg002_two_arms.sh` is
# the same with one binary and an environment switch). Writes into
# <out-dir>: each arm's VCF, `genotype.tsv` (right genotypes by period class and depth, the
# number to compare), `sweep.tsv` (allele-exact precision and recall by QUAL), per-arm verdicts,
# and `flips_<depth>.txt` (tracts whose verdict changed, and which way). The tract ground,
# `tier.bed`, is typed with this checkout's `ng_typed_region_dump`, built first if missing.
#
# Env: BENCHMARKS (default: the main checkout's benchmarks/, which holds the data; mount it).
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
out=${1:?usage: $0 <out-dir> <base-binary> <experiment-binary>}; base_bin=${2:?base binary}; experiment_bin=${3:?experiment binary}
# TIER_KIND: which typed regions are the scored ground — ssr_locus (default) or ssr_bundle. A bundle
# has no single period (the dump writes `.`); PERIOD_FOR_BUNDLES sets the one it is scored under
# (default 2, so bundles count as period 2 and above). The ground is typed with this checkout's
# dump: a build that retypes bundles as ordinary sequence must be scored on a tier.bed made from
# main's typing, copied into <out-dir> first.
tier_kind=${TIER_KIND:-ssr_locus}; bundle_period=${PERIOD_FOR_BUNDLES:-2}
data=${BENCHMARKS:-$(cd "$root/.." && pwd)/pop_var_caller/benchmarks}
ref=$data/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
truth=$data/ssr_hg002/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz
mkdir -p "$out"
if [[ ! -s $out/tier.bed ]]; then
    dump=$root/target-container/release/examples/ng_typed_region_dump
    [[ -x $dump ]] || (cd "$root" && cargo build --release --example ng_typed_region_dump)
    sort -k1,1 -k2,2n "$data/ssr_hg002/regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed" | cut -f1-3 > "$out/tier_sorted.bed"
    NG_REFERENCE_CHECK=skip "$dump" "$ref" "$out/tier_sorted.bed" calling 2>/dev/null \
        | awk -F'\t' -v k="$tier_kind" -v p="$bundle_period" '!/^#/ && $4==k { print $1"\t"($2-1)"\t"$3"\t"($6=="." ? p : $6) }' > "$out/tier.bed"
fi
for cov in 30x 50x; do
    bam=$data/ssr_hg002/bam/$cov/HG002_TR_v1.0.1_Tier_${cov}.bam
    for arm in base experiment; do
        vcf=$out/${arm}_${cov}.vcf
        if [[ ! -s $vcf ]]; then
            if [[ $arm == experiment ]]; then bin=$experiment_bin; else bin=$base_bin; fi
            NG_REFERENCE_CHECK=skip "$bin" call-from-alignments --reference "$ref" --catalog "$ref.repeats.parquet" \
                --alignment "$bam" --regions "$out/tier_sorted.bed" --output "$vcf" --defaults > "$out/${arm}_${cov}.log" 2>&1
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
