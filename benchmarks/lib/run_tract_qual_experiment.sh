#!/usr/bin/env bash
# The tract QUAL experiment, end to end — `doc/devel/ng/spec/calling_loop_ssr.md` §3.3.
#
#   benchmarks/lib/run_tract_qual_experiment.sh <ground> [out dir]
#
# `<ground>` is one of `per_sample`, `tandem_repeat_tier`, `simulator` or `all`.
# The output directory defaults to `tmp/tract_qual` in the repository, and ends
# holding two tables written by `tract_qual_experiment.py`:
#
#   calibration.tsv  records binned by QUAL against the share truly variant
#   sweep.tsv        precision and recall as a QUAL threshold sweeps
#
# Both are appended to, so a second ground adds rows rather than replacing them;
# delete the directory to start over.
#
# WHAT THE THREE GROUNDS ARE FOR, AND WHY THERE ARE THREE
#
#   per_sample          The GIAB trio's 100 random confident intervals, at 30x
#                       and 50x — the ground every standing ng number was
#                       measured on, so a figure here is read against C4's.
#                       It holds about 4,200 bases of repeat tract a sample and
#                       ng writes about 50 tract records on it, which is too few
#                       to see a one-in-a-thousand error rate. It is here for
#                       continuity, not for the calibration.
#
#   tandem_repeat_tier  GIAB's HG002 tandem-repeat benchmark, 50,000 Tier
#                       intervals with assembly-based truth, at the same two
#                       depths. About 6,400 tract records — this is where the
#                       calibration can actually be read. One sample.
#
#   simulator           `examples/ng_tract_simulator`: tracts whose genotypes we
#                       chose, sequenced under a slippage we set. The only
#                       ground where the truth is exact, where the slippage can
#                       be moved away from what the caller assumes, and where
#                       the fitted-against-defaulted split has two sides — no
#                       command fits a parameters file yet, so every cell on
#                       either GIAB ground is `Defaulted`.
#
# ARMS
#
#   ng                  the inherited site-quality fold as it stands (spec
#                       §3.3's arm A), `--defaults`
#   ng_fitted           the same caller handed the slippage the reads were drawn
#                       under — simulator only, since only there is it known
#   production          the existing caller on the same ground (arm C): the
#                       `high-recall` preset on `per_sample`, and `ssr-call` on
#                       the tandem-repeat ground
#
# Env: NG_BIN, NG_EXAMPLE_DIR, GIAB_ROOT, SSR_HG002_ROOT, DEPTHS, SIM_TRACTS,
#      SIM_SAMPLES, SIM_SLIPS.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCHMARKS="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$BENCHMARKS/.." && pwd)"

GROUND="${1:?ground: per_sample | tandem_repeat_tier | simulator | all}"
OUT_DIR="${2:-$REPO/tmp/tract_qual}"

GIAB_ROOT="${GIAB_ROOT:-$BENCHMARKS/giab}"
SSR_HG002_ROOT="${SSR_HG002_ROOT:-$BENCHMARKS/ssr_hg002}"
REFERENCE="$GIAB_ROOT/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
NG_BIN="${NG_BIN:-$REPO/target-container/release/pop_var_caller_exp}"
NG_EXAMPLE_DIR="${NG_EXAMPLE_DIR:-$REPO/target-container/release/examples}"
SCORER="$HERE/tract_qual_experiment.py"

# **The scorer needs a Python and the two environments have different ones.**
# The development container ships `python3` and no `uv`; the macOS host is the
# other way round, and this repository's rule is that host Python goes through
# `uv run`. Nothing here needs a package, so either interpreter serves.
if [[ -n "${PYTHON_RUNNER:-}" ]]; then
    read -r -a PYTHON <<< "$PYTHON_RUNNER"
elif command -v python3 > /dev/null; then
    PYTHON=(python3)
else
    PYTHON=(uv run --no-project python)
fi
DEPTHS="${DEPTHS:-30x 50x}"

CALIBRATION="$OUT_DIR/calibration.tsv"
SWEEP="$OUT_DIR/sweep.tsv"
GROUND_DIR="$OUT_DIR/ground"
mkdir -p "$GROUND_DIR"

# The typed regions a run routes on, as a `chrom start end period` BED of the
# repeat tracts alone — the ground both arms are scored on.
#
# `calling` is the floors `call-from-alignments` actually routes on. The dump's
# other setting, `catalog`, is what the file was stored at, and the two differ
# by a lot: on HG002's confident intervals the calling floors call 335 stretches
# a repeat tract where the storage floors call many more.
build_tract_ground() {
    local regions="$1" out="$2"
    if [[ -s "$out" ]]; then
        echo "[ground] reusing $out ($(wc -l < "$out") tracts)"
        return
    fi
    echo "[ground] typing $regions"
    NG_REFERENCE_CHECK=skip "$NG_EXAMPLE_DIR/ng_typed_region_dump" \
        "$REFERENCE" "$regions" calling 2>/dev/null \
        | awk -F'\t' '!/^#/ && $4=="ssr_locus" { print $1"\t"($2-1)"\t"$3"\t"$6 }' \
        > "$out"
    echo "[ground] $out holds $(wc -l < "$out") tracts"
}

score() {
    local arm="$1" ground="$2" depth="$3" sample="$4"
    local truth="$5" query="$6" confident="$7" tracts="$8"
    [[ -f "$query" ]] || { echo "!! missing $query — skipping" >&2; return; }
    "${PYTHON[@]}" "$SCORER" \
        --reference "$9" --truth "$truth" --query "$query" \
        --confident-bed "$confident" --tract-bed "$tracts" \
        --arm "$arm" --ground "$ground" --depth "$depth" --sample "$sample" \
        --calibration-out "$CALIBRATION" --sweep-out "$SWEEP"
}

# ---------------------------------------------------------------------------

run_per_sample() {
    local bed_dir="$GIAB_ROOT/per_sample/bed"
    local truth_dir="$GIAB_ROOT/per_sample/vcf"
    for sample in HG002 HG003 HG004; do
        local confident="$bed_dir/${sample}_bench_azar_merged_100.bed"
        local truth="$truth_dir/${sample}_GRCh38_1_22_v4.2.1_benchmark.selected_100.vcf.gz"
        local tracts="$GROUND_DIR/per_sample_${sample}.bed"
        build_tract_ground "$confident" "$tracts"
        for depth in $DEPTHS; do
            local results="$GIAB_ROOT/results/per_sample/$depth"
            # ng's ungated file, so the sweep starts below the runner's own
            # QUAL 30 floor rather than at it.
            score ng per_sample "$depth" "$sample" "$truth" \
                "$results/ng/${sample}.raw.vcf" "$confident" "$tracts" "$REFERENCE"
            score production per_sample "$depth" "$sample" "$truth" \
                "$results/high-recall/${sample}.vcf" "$confident" "$tracts" "$REFERENCE"
        done
    done
}

run_tandem_repeat_tier() {
    local confident="$GROUND_DIR/tier_sorted.bed"
    if [[ ! -s "$confident" ]]; then
        sort -k1,1 -k2,2n \
            "$SSR_HG002_ROOT/regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed" \
            | cut -f1-3 > "$confident"
    fi
    local truth="$SSR_HG002_ROOT/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz"
    local tracts="$GROUND_DIR/tier.bed"
    build_tract_ground "$confident" "$tracts"
    for depth in $DEPTHS; do
        score ng tandem_repeat_tier "$depth" HG002 "$truth" \
            "$SSR_HG002_ROOT/results/ng/HG002_${depth}.raw.vcf" \
            "$confident" "$tracts" "$REFERENCE"
        score production tandem_repeat_tier "$depth" HG002 "$truth" \
            "$SSR_HG002_ROOT/results/ours/vcf/HG002_${depth}.ssr.vcf" \
            "$confident" "$tracts" "$REFERENCE"
    done
}

# The simulator arm, at a ladder of true slippage levels.
#
# `SIM_SLIPS` names the share of reads that report a length other than their
# allele's. **0.10 is the number the caller assumes when nothing is fitted**, so
# that rung is the case where the shipped model is exactly right and every other
# rung is a case where it is not — which is the risk §3.3 states.
run_simulator() {
    local tracts_n="${SIM_TRACTS:-4000}"
    local samples="${SIM_SAMPLES:-3}"
    local slips="${SIM_SLIPS:-0.02 0.10 0.25}"
    for slip in $slips; do
        for depth in $DEPTHS; do
            local reads="${depth%x}"
            local dir="$OUT_DIR/sim_slip${slip}_${depth}"
            echo "[simulator] slip=$slip depth=$depth -> $dir"
            mkdir -p "$dir"
            "$NG_EXAMPLE_DIR/ng_tract_simulator" "$dir" \
                "tracts=$tracts_n" "samples=$samples" "depth=$reads" \
                "slip_share=$slip" > "$dir/simulator.log" 2>&1
            samtools faidx "$dir/reference.fa"
            local alignments=()
            for bam in "$dir"/sim*.bam; do
                samtools index "$bam"
                alignments+=(--alignment "$bam")
            done
            "$NG_BIN" repeat-catalog --reference "$dir/reference.fa" \
                --output "$dir/reference.fa.repeats.parquet" > "$dir/catalog.log" 2>&1

            # **The fixture's own oracle, and the run stops if it fails.** The
            # region typing has to call exactly the tracts the simulator laid
            # down: one more is a repeat the flanks grew by accident, one fewer
            # is a tract the routing merged into a cluster, and either way the
            # ground the scorer restricts to is not the ground the reads were
            # drawn on.
            NG_REFERENCE_CHECK=skip "$NG_EXAMPLE_DIR/ng_typed_region_dump" \
                "$dir/reference.fa" "$dir/confident.bed" calling 2>/dev/null \
                | awk -F'\t' '!/^#/ && $4=="ssr_locus" { print $1"\t"($2-1)"\t"$3"\t"$6 }' \
                > "$dir/typed_tracts.bed"
            if ! diff -q <(sort -k2,2n "$dir/typed_tracts.bed") \
                        <(sort -k2,2n "$dir/tracts.bed") > /dev/null; then
                echo "!! the typed tracts are not the injected tracts in $dir" >&2
                diff <(sort -k2,2n "$dir/typed_tracts.bed") \
                     <(sort -k2,2n "$dir/tracts.bed") | head -10 >&2
                exit 1
            fi

            NG_REFERENCE_CHECK=skip "$NG_BIN" call-from-alignments \
                --reference "$dir/reference.fa" \
                --catalog "$dir/reference.fa.repeats.parquet" \
                "${alignments[@]}" --regions "$dir/confident.bed" \
                --output "$dir/ng.vcf" --defaults > "$dir/ng.log" 2>&1

            # The same reads scored under the model that made them. The empty
            # slippage table is deleted and the true rows appended at the end,
            # because an array-of-tables closes the table it sits in.
            grep -v '^slippage_by_stratum_and_group = \[\]$' \
                "$dir/ng.parameters.toml" > "$dir/fitted.parameters.toml"
            cat "$dir/slippage_rows.toml" >> "$dir/fitted.parameters.toml"
            NG_REFERENCE_CHECK=skip "$NG_BIN" call-from-alignments \
                --reference "$dir/reference.fa" \
                --catalog "$dir/reference.fa.repeats.parquet" \
                "${alignments[@]}" --regions "$dir/confident.bed" \
                --output "$dir/ng_fitted.vcf" \
                --parameters "$dir/fitted.parameters.toml" \
                > "$dir/ng_fitted.log" 2>&1

            for arm in ng ng_fitted; do
                score "$arm" "simulator_slip${slip}" "$depth" pooled \
                    "$dir/truth.vcf" "$dir/${arm}.vcf" \
                    "$dir/confident.bed" "$dir/tracts.bed" "$dir/reference.fa"
            done
        done
    done
}

case "$GROUND" in
    per_sample) run_per_sample ;;
    tandem_repeat_tier) run_tandem_repeat_tier ;;
    simulator) run_simulator ;;
    all) run_per_sample; run_tandem_repeat_tier; run_simulator ;;
    *) echo "unknown ground: $GROUND" >&2; exit 2 ;;
esac

echo
echo "calibration : $CALIBRATION"
echo "sweep       : $SWEEP"
