#!/usr/bin/env bash
# run_prior_comparison.sh — re-call the existing HG002 .psp files under BOTH genotype
# priors, so the two arms differ by nothing but `PVC_SSR_MARGINALIZED_PRIOR`.
#
# Answers the first half of Q5 in doc/devel/ng/spec/calling_priors.md: does the
# marginalized genotype prior genotype repeat tracts better than the plug-in one?
# Result: doc/devel/reports/ssr_prior_hg002_single_sample_2026-08-18.md
#
# The pileups are NOT rebuilt — run_ours_coverages.sh made them, and re-using them is
# what keeps the two arms differing only in the prior. The native 300x rung is skipped
# by default: at that depth the read likelihood swamps any prior.
#
# Run inside the dev container (the binary is a Linux build):
#   ./scripts/dev.sh bash benchmarks/ssr_hg002/src/run_prior_comparison.sh
# Env: BIN, THREADS, COVERAGES.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPO="$(cd "$ROOT/../.." && pwd)"
BIN="${BIN:-$REPO/target-container/release/pop_var_caller}"
# The Tier-restricted catalog, not the genome-wide one: ssr-call's burn-in has to
# sample loci the BAM actually covers (README, "PREP GOTCHAS").
CAT="$ROOT/catalog/HG002_Tier_restricted.cat"
THREADS="${THREADS:-8}"
COVERAGES="${COVERAGES:-50 30 20 15}"
OUT="$ROOT/results/ours_prior_comparison"
mkdir -p "$OUT/plugin" "$OUT/marginalized" "$OUT/logs"

for cov in $COVERAGES; do
  psp="$ROOT/results/ours/psp/HG002_${cov}x.ssr.psp"
  if [[ ! -f "$psp" ]]; then
    echo "!! missing $psp — run run_ours_coverages.sh first" >&2
    continue
  fi

  echo "[plug-in]      ${cov}x"
  "$BIN" ssr-call --catalog "$CAT" --output "$OUT/plugin/HG002_${cov}x.ssr.vcf" \
    --threads "$THREADS" "$psp" > "$OUT/logs/plugin_${cov}x.log" 2>&1

  echo "[marginalized] ${cov}x"
  PVC_SSR_MARGINALIZED_PRIOR=1 "$BIN" ssr-call --catalog "$CAT" \
    --output "$OUT/marginalized/HG002_${cov}x.ssr.vcf" \
    --threads "$THREADS" "$psp" > "$OUT/logs/marginalized_${cov}x.log" 2>&1

  plugin_n=$(grep -vc '^#' "$OUT/plugin/HG002_${cov}x.ssr.vcf")
  marg_n=$(grep -vc '^#' "$OUT/marginalized/HG002_${cov}x.ssr.vcf")
  echo "   ${cov}x records: plug-in=$plugin_n marginalized=$marg_n"
done
echo "done — score with: cd $ROOT && uv run src/prior_genotype_accuracy.py"
