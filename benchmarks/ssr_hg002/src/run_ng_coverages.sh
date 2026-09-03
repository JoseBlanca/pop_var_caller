#!/usr/bin/env bash
# run_ng_coverages.sh — run ng over the HG002 tandem-repeat benchmark at every
# coverage rung, so its repeat-tract calls can be scored against GIAB's
# assembly-based truth.
#
#   bash benchmarks/ssr_hg002/src/run_ng_coverages.sh [COVERAGE ...]
#
# The sibling of `run_ours_coverages.sh`, `run_hipstr_coverages.sh` and
# `run_freebayes_coverages.sh`, and it differs from all three in one way worth
# knowing: ng goes from alignments to VCF in ONE process, so there is no pileup
# stage and nothing on disk between the BAM and the calls.
#
# WHY THIS BENCHMARK AND NOT `giab/per_sample`
#
# The trio's per-sample benchmark holds 100 random confident intervals — about
# 4,200 bases of repeat tract a sample, on which ng writes about 50 tract
# records. That is enough to say whether tracts are called at all and far too
# few to measure how often a record at QUAL 200 is wrong. This benchmark is
# 50,000 tandem-repeat intervals over 6.1 Mb with 36,497 truth records, and ng
# writes about 6,400 tract records on it. It is one sample, so it says nothing
# about a cohort.
#
# The truth is assembly-based, so it does not share a short-read caller's
# stutter-versus-allele failure — which is what makes it usable as truth here.
# See `benchmarks/ssr_hg002/README.txt`.
#
# TWO FILES A COVERAGE, as the other runners write:
#   HG002_<cov>.raw.vcf   every record ng emitted, no QUAL floor
#   HG002_<cov>.vcf       QUAL >= MIN_QUAL, the comparable set
# A QUAL sweep needs the first; a headline number reads the second.
#
# Env: NG_BIN, NG_CATALOG, SSR_HG002_ROOT, GIAB_ROOT, MIN_QUAL, OUT_DIR, DRY_RUN.
set -euo pipefail

SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${SSR_HG002_ROOT:-$(cd "$SRC_DIR/.." && pwd)}"
BENCHMARKS="$(cd "$ROOT/.." && pwd)"
REPO="$(cd "$BENCHMARKS/.." && pwd)"

GIAB_ROOT="${GIAB_ROOT:-$BENCHMARKS/giab}"
REFERENCE="$GIAB_ROOT/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
NG_CATALOG="${NG_CATALOG:-${REFERENCE}.repeats.parquet}"
OUT_DIR="${OUT_DIR:-$ROOT/results/ng}"
MIN_QUAL="${MIN_QUAL:-30}"
REGIONS_SOURCE="$ROOT/regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed"

COVERAGES=("$@")
if (( ${#COVERAGES[@]} == 0 )); then
    COVERAGES=(30x 50x)
fi

discover_bin() {
    if [[ -z "${NG_BIN:-}" ]]; then
        local candidate
        for candidate in \
            "$REPO/target-container/release/pop_var_caller_exp" \
            "$REPO/target/release/pop_var_caller_exp"; do
            if [[ -x "$candidate" ]] && "$candidate" --version >/dev/null 2>&1; then
                NG_BIN="$candidate"
                break
            fi
        done
    fi
    [[ -n "${NG_BIN:-}" && -x "${NG_BIN}" ]] || {
        echo "no pop_var_caller_exp binary; build it or set NG_BIN" >&2
        exit 1
    }
}

discover_bin
mkdir -p "$OUT_DIR"

# ng walks its regions in order, and the shipped Tier BED is not sorted.
REGIONS="$OUT_DIR/tier_sorted.bed"
if [[ ! -s "$REGIONS" ]]; then
    sort -k1,1 -k2,2n "$REGIONS_SOURCE" | cut -f1-3 > "$REGIONS"
fi

echo "binary    : $NG_BIN"
echo "reference : $REFERENCE"
echo "catalog   : $NG_CATALOG"
echo "regions   : $REGIONS ($(wc -l < "$REGIONS") intervals)"
echo "model     : --defaults (nothing fitted)"
echo "out dir   : $OUT_DIR"
echo

# header lines pass through; data rows kept iff QUAL is numeric and >= min.
FILTER_AWK='/^#/ { print; next } $6 != "." && ($6 + 0) >= min { print }'

for cov in "${COVERAGES[@]}"; do
    bam="$ROOT/bam/${cov}/HG002_TR_v1.0.1_Tier_${cov}.bam"
    [[ -f "$bam" ]] || { echo "!! missing $bam — skipping $cov" >&2; continue; }
    raw="$OUT_DIR/HG002_${cov}.raw.vcf"
    gated="$OUT_DIR/HG002_${cov}.vcf"
    log="$OUT_DIR/HG002_${cov}.log"

    echo "=================================================================="
    echo "coverage  : $cov"
    echo "alignment : $bam"
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "DRY-RUN"
        continue
    fi
    t0=$(date +%s)
    NG_REFERENCE_CHECK="${NG_REFERENCE_CHECK:-skip}" "$NG_BIN" call-from-alignments \
        --reference "$REFERENCE" \
        --catalog "$NG_CATALOG" \
        --alignment "$bam" \
        --regions "$REGIONS" \
        --output "$raw" \
        --defaults 2>&1 | tee "$log"
    t1=$(date +%s)
    awk -v min="$MIN_QUAL" "$FILTER_AWK" "$raw" > "$gated"
    echo
    echo "elapsed                     : $((t1 - t0)) s"
    echo "records (raw)               : $(grep -vc '^#' "$raw" || true)"
    echo "records (QUAL >= $MIN_QUAL) : $(grep -vc '^#' "$gated" || true)"
    echo "tract records (raw)         : $(grep -c 'STR;' "$raw" || true)"
    echo
done

echo "done. VCFs under: $OUT_DIR"
