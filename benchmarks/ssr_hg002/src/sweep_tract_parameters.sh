#!/usr/bin/env bash
# One ng run over the HG002 tandem-repeat benchmark at a stated set of
# repeat-tract parameters, scored for genotype accuracy.
#
#   sweep_tract_parameters.sh <label> [options]
#
#     --share S            one slip share for every stratum (default: the run's own,
#                          i.e. HipSTR's shipped 0.10 through the defaults file)
#     --base B --slope S   a slip share that rises with tract length instead
#     --shorter S          the share of slips that shorten the tract (default 0.50)
#     --fall-off F         how fast two-repeat slips fall off (default 0.05)
#     --rows FILE          slippage rows written elsewhere, e.g. a fit's output
#     --outlier W          the repeat-tract outlier weight (default: the run's own 0.01)
#     --concentration C    the fallback length-spectrum concentration — how much
#                          prior belief, in chromosomes, is spread over a tract's
#                          candidate lengths (default: the run's own 1.0). **This
#                          is the only thing on the tract path that moves the
#                          het/hom balance**: the genotype row is a marginalised
#                          Dirichlet-multinomial, so with K candidates the share
#                          of prior mass on heterozygous genotypes is
#                          (K-1)C / (K(C+1)) — 42% at K=6, C=1, rising toward
#                          (K-1)/K as C grows and falling to nothing as C -> 0.
#     --coverage C         30x (default) or 50x
#     --out DIR            where to put the run (default: tmp/tract_sweep/<label>)
#
# This is the instrument behind
# `doc/devel/reports/ng_tract_genotype_improvement_2026-09-02.md` §2. Its result:
# no flat stutter setting beats the shipped one, a fitted per-stratum set is
# worth about a fifth of a point, and the outlier weight — the bound on how far
# one read may pull a genotype, inherited at 0.01 and never measured — is worth
# about four tenths.
#
# **Three of these settings are the same dial.** The slip share, the outlier
# weight and the length-spectrum concentration all trade a spurious heterozygote
# against a collapsed one, and each is at or beside its own peak on this
# benchmark. A fourth point on any of them is not worth a run; what is left is
# not a mis-set balance between the two error classes.
#
# **A parameters file written before a shipped default moved cannot be replayed
# unedited**: the run refuses a value that disagrees with the built-in while
# claiming `defaulted`. Pass `--outlier` with the value that file carries, or
# with the current default, and say which in the label.
#
# THE CONTROL, and run it first after any change here: with no options at all
# the run must be BYTE-IDENTICAL to the `--defaults` run, because it is handed
# that run's own parameters file unedited. It is.
#
# A second, weaker control is worth knowing about: `--share 0.10 --shorter 0.50
# --fall-off 0.05` writes the shipped numbers out as rows, and that run differs
# from `--defaults` by 2 genotypes in 3,648. A supplied row rebuilds the
# part-repeat shares as a twentieth of the whole-repeat mass where the shipped
# model states them as 0.01 each, so a sweep's rows are comparable with each
# other but carry a 0.05-point offset against a `--defaults` baseline.
#
# Env: NG_BIN, SSR_HG002_ROOT, GIAB_ROOT, TRACT_GROUND_DIR.
set -euo pipefail

# **The tooling comes from this script's own tree and the data from
# `SSR_HG002_ROOT`, and they are deliberately two roots.** The benchmark data is
# gitignored and lives only in the primary checkout, so a run driven from a
# worktree reads the data over there and must still use the binary, the scorer
# and the ground built over here — resolving both from one root gives the
# worktree's changes no effect, silently.
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_REPO="$(cd "$SRC_DIR/../../.." && pwd)"
ROOT="${SSR_HG002_ROOT:-$(cd "$SRC_DIR/.." && pwd)}"
BENCHMARKS="$(cd "$ROOT/.." && pwd)"
REPO="$SCRIPT_REPO"
GIAB_ROOT="${GIAB_ROOT:-$BENCHMARKS/giab}"
REFERENCE="$GIAB_ROOT/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
NG_BIN="${NG_BIN:-$SCRIPT_REPO/target-container/release/pop_var_caller_exp}"
SCORER="$SCRIPT_REPO/benchmarks/lib/tract_qual_experiment.py"
PYTHON=$(command -v python3 > /dev/null && echo python3 || echo "uv run --no-project python")

LABEL="${1:?a label for this setting, e.g. outlier0.10}"; shift
SHARE=""; BASE=""; SLOPE=""; SHORTER=0.50; FALLOFF=0.05; ROWS=""; OUTLIER=""
CONCENTRATION=""; COVERAGE=30x; OUT=""
while (( $# )); do
    case "$1" in
        --share) SHARE="$2"; shift 2 ;;
        --base) BASE="$2"; shift 2 ;;
        --slope) SLOPE="$2"; shift 2 ;;
        --shorter) SHORTER="$2"; shift 2 ;;
        --fall-off) FALLOFF="$2"; shift 2 ;;
        --rows) ROWS="$2"; shift 2 ;;
        --outlier) OUTLIER="$2"; shift 2 ;;
        --concentration) CONCENTRATION="$2"; shift 2 ;;
        --coverage) COVERAGE="$2"; shift 2 ;;
        --out) OUT="$2"; shift 2 ;;
        *) echo "unknown option: $1" >&2; exit 2 ;;
    esac
done
OUT="${OUT:-$REPO/tmp/tract_sweep/$LABEL}"
mkdir -p "$OUT"

# The tract ground and the confident regions, as `run_tract_qual_experiment.sh`
# builds them.
#
# **Overridable, because the benchmark data and the checkout can be different
# trees.** `SSR_HG002_ROOT` points at the data, which on this project lives only
# in the primary checkout; a run driven from a worktree wants the ground the
# worktree built, not the one beside the data.
GROUND_DIR="${TRACT_GROUND_DIR:-$REPO/tmp/tract_qual/ground}"
CONFIDENT="$GROUND_DIR/tier_sorted.bed"
TRACTS="$GROUND_DIR/tier.bed"
for needed in "$CONFIDENT" "$TRACTS"; do
    [[ -s "$needed" ]] || {
        echo "missing $needed — run benchmarks/lib/run_tract_qual_experiment.sh tandem_repeat_tier first" >&2
        exit 1
    }
done

DEFAULTS_PARAMETERS="$ROOT/results/ng/HG002_${COVERAGE}.raw.parameters.toml"
[[ -f "$DEFAULTS_PARAMETERS" ]] || {
    echo "missing $DEFAULTS_PARAMETERS — run src/run_ng_coverages.sh $COVERAGE first" >&2
    exit 1
}

# **The parameters file is the defaults run's own, edited.** Reusing what the
# caller wrote means the six groups of numbers this sweep is not about cannot
# drift from what the baseline scored under.
if [[ -n "$OUTLIER" ]]; then
    # Its own commentary requires a value somebody typed to say `supplied`, so
    # that a run cannot report a number we chose as one it inherited.
    sed "s|^repeat_tract_outlier_weight = .*|repeat_tract_outlier_weight = { value = $OUTLIER, warrant = \"supplied\" }|" \
        "$DEFAULTS_PARAMETERS" > "$OUT/parameters.toml"
    grep -q "value = $OUTLIER" "$OUT/parameters.toml" || {
        echo "the outlier weight was not set" >&2; exit 1; }
else
    cp "$DEFAULTS_PARAMETERS" "$OUT/parameters.toml"
fi

if [[ -n "$CONCENTRATION" ]]; then
    # Same rule as the outlier weight: a value somebody typed says `supplied`,
    # so a run cannot report a number we chose as one it inherited.
    sed -i.bak "s|^fallback_length_spectrum_concentration = .*|fallback_length_spectrum_concentration = { value = $CONCENTRATION, warrant = \"supplied\" }|" \
        "$OUT/parameters.toml"
    rm -f "$OUT/parameters.toml.bak"
    grep -q "value = $CONCENTRATION, warrant = \"supplied\"" "$OUT/parameters.toml" || {
        echo "the length-spectrum concentration was not set" >&2; exit 1; }
fi

if [[ -n "$SHARE" || -n "$BASE" || -n "$ROWS" ]]; then
    if [[ -z "$ROWS" ]]; then
        ROWS="$OUT/rows.toml"
        if [[ -n "$SHARE" ]]; then
            $PYTHON "$SRC_DIR/tract_slippage_rows.py" --share "$SHARE" \
                --shorter "$SHORTER" --fall-off "$FALLOFF" --out "$ROWS"
        else
            $PYTHON "$SRC_DIR/tract_slippage_rows.py" --base "$BASE" --slope "$SLOPE" \
                --shorter "$SHORTER" --fall-off "$FALLOFF" --out "$ROWS"
        fi
    fi
    # The empty array goes and the rows are appended at the END of the file — an
    # array-of-tables closes the table it sits in.
    grep -v '^slippage_by_stratum_and_group = \[\]$' "$OUT/parameters.toml" > "$OUT/parameters.tmp"
    mv "$OUT/parameters.tmp" "$OUT/parameters.toml"
    cat "$ROWS" >> "$OUT/parameters.toml"
fi

NG_REFERENCE_CHECK="${NG_REFERENCE_CHECK:-skip}" "$NG_BIN" call-from-alignments \
    --reference "$REFERENCE" \
    --catalog "${REFERENCE}.repeats.parquet" \
    --alignment "$ROOT/bam/${COVERAGE}/HG002_TR_v1.0.1_Tier_${COVERAGE}.bam" \
    --regions "$CONFIDENT" \
    --output "$OUT/calls.vcf" \
    --parameters "$OUT/parameters.toml" > "$OUT/run.log" 2>&1

$PYTHON "$SCORER" \
    --reference "$REFERENCE" \
    --truth "$ROOT/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz" \
    --query "$OUT/calls.vcf" \
    --confident-bed "$CONFIDENT" --tract-bed "$TRACTS" \
    --arm "$LABEL" --ground tandem_repeat_tier --depth "$COVERAGE" --sample HG002 \
    --genotype-sample HG002 \
    --calibration-out "$(dirname "$OUT")/calibration.tsv" \
    --sweep-out "$(dirname "$OUT")/sweep.tsv" \
    --genotype-out "$(dirname "$OUT")/genotype.tsv"

echo "[$LABEL @ $COVERAGE] tract records: $(grep -c 'STR;' "$OUT/calls.vcf" || true)"
