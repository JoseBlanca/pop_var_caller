#!/usr/bin/env bash
# Copy one run of the stuttering benchmark's first pass into this worktree, so a build can be
# called and scored on it without the first pass's 9 callers' outputs or its fit.
#
#   scripts/ng_str_accuracy/copy_run.sh <run>        e.g. pcr_1s_15x, pcr_free_50s_60x
#
# Copies the reads, the truth, the reference and its catalog, the true-stutter parameters, and
# GATK's and HipSTR's calls (so they are scored beside ng), plus main's calls as of 4b528afe
# (`ng_main_true_stutter.vcf`, the "before" of doc/devel/reports/implementations/
# str_accuracy_2026-10-10.md). Nothing ng wrote is copied: psps are made fresh by the build
# being measured (`call_and_score.sh`).
#
# Env: FIRST_PASS (default: the tract-ceiling worktree's tmp/stutter_bench/first_pass, where the
# 16 runs of 2026-10-10 live), STR_RUNS_DIR (default: <repo>/tmp/str_accuracy/runs).
# Sizes: a one-sample run is 25–190 MB; a 50-sample run 0.9 GB at 5x to 8.9 GB at 60x.
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
src=${FIRST_PASS:-$(cd "$root/.." && pwd)/pop_var_caller-tract-ceiling/tmp/stutter_bench/first_pass}/${1:?usage: $0 <run>}
dst=${STR_RUNS_DIR:-$root/tmp/str_accuracy/runs}/$1
[[ -d $src ]] || { echo "no run at $src" >&2; exit 1; }
mkdir -p "$dst"
for f in reference.fa reference.fa.fai reference.dict reference.fa.repeats.parquet segments.bed confident.bed \
         segment_alleles.tsv segment_genotypes.tsv truth.vcf libraries.tsv stutter_rates.tsv simulate.log \
         ng_defaults.parameters.toml true_stutter.parameters.toml gatk.vcf.gz gatk.vcf.gz.tbi \
         hipstr.vcf.gz hipstr.vcf.gz.tbi hipstr_fixed_stutter.vcf.gz hipstr_fixed_stutter.vcf.gz.tbi \
         ng_main_true_stutter.vcf; do
    [[ -e $src/$f ]] && cp "$src/$f" "$dst/"
done
[[ -d $dst/bams ]] || cp -R "$src/bams" "$dst/"
echo "copied $1 to $dst"
