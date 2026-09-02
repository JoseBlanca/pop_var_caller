#!/usr/bin/env bash
# Run ng — the experimental caller — over the GIAB `per_sample` dataset, so its
# calls can be scored beside the production caller's and freebayes'.
#
#   benchmarks/giab/src/run_ng_per_sample.sh [COVERAGE] [SAMPLE ...]
#
# Mirrors run_ours_per_sample.sh and run_freebayes_per_sample.sh: the per_sample
# dataset holds three independent single-sample callsets (HG002, HG003, HG004),
# each with its OWN random 100-region BED and its OWN GIAB truth VCF, so every
# sample is called on its own, restricted to its confident regions.
#
# ng goes from alignments to VCF in ONE process — there is no .psp stage — so a
# run is a single `call-from-alignments` invocation per sample. Three things
# about that command shape this script:
#
#   * It needs a tandem-repeat catalog built from the same reference. The
#     default is `<reference>.repeats.parquet`, which is where
#     `pop_var_caller_exp repeat-catalog` writes it; this script builds it once
#     if it is not there, which on GRCh38 is the slow part of a first run.
#   * It needs either a fitted parameters file or `--defaults`. No command
#     writes a fitted file yet, so every run here is `--defaults`: no
#     base-quality calibration, no contamination, no inbreeding. The genotypes
#     are what the reads alone say under those assumptions.
#   * It writes a `<output-stem>.parameters.toml` beside the VCF recording the
#     numbers it scored with.
#
# Fairness with the other two callers. ng applies no QUAL floor of its own,
# while the production caller defaults to `--min-qual 30` and freebayes gets a
# post-call QUAL >= 30 gate in its runner. So this script writes BOTH, exactly
# as the freebayes runner does:
#   {sample}.raw.vcf  — every record ng emitted (no QUAL filter)
#   {sample}.vcf      — QUAL >= MIN_QUAL (the headline, comparable set)
#
# Args:
#   COVERAGE   bam/ coverage subdir to use (default: 300x)
#   SAMPLE...  subset of {HG002,HG003,HG004} to run (default: all three)
#
# Env overrides:
#   NG_BIN       binary to invoke (default: auto-detect host/container build)
#   NG_CATALOG   repeat catalog path (default: <reference>.repeats.parquet)
#   MIN_QUAL     QUAL floor for the gated output (default 30)
#   DRY_RUN=1    print the commands instead of running them

set -euo pipefail

# benchmarks/giab/src -> benchmarks/giab -> benchmarks -> repo root
SRC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BENCH_DIR="$(cd "$SRC_DIR/.." && pwd)"
PROJECT_ROOT="$(cd "$BENCH_DIR/../.." && pwd)"

COVERAGE="${1:-300x}"
shift || true
SAMPLES=("$@")
if (( ${#SAMPLES[@]} == 0 )); then
    SAMPLES=(HG002 HG003 HG004)
fi

MIN_QUAL="${MIN_QUAL:-30}"
REFERENCE="$BENCH_DIR/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
BAM_DIR="$BENCH_DIR/per_sample/bam/$COVERAGE"
BED_DIR="$BENCH_DIR/per_sample/bed"
OUT_DIR="$BENCH_DIR/results/per_sample/$COVERAGE/ng"
NG_CATALOG="${NG_CATALOG:-${REFERENCE}.repeats.parquet}"

# --- sample -> alignment filename / BED filename ---------------------------
# Same two naming schemes as the other two runners (kept in sync):
#   300x  — HG002 is a CRAM; HG003/HG004 are the `.m5.bam` reheadered copies
#           (ng does not need @SQ M5, but those are the files we keep around).
#   else  — downsampled tiers (5x..50x): `{sample}.{cov}.seed42.bam`.
align_file() {
    local dir="$BAM_DIR"
    if [[ "$COVERAGE" != "300x" ]]; then
        local bam="$dir/${1}.${COVERAGE}.seed42.bam"
        if [[ -f "$bam" ]]; then echo "$bam"; return 0; fi
        echo "missing alignment: $bam" >&2; return 1
    fi
    case "$1" in
        HG002) echo "$dir/HG002_reads_selected_100_rg.cram" ;;
        HG003|HG004)
            local m5="$dir/${1}_bench_azar_merged_100.sorted.m5.bam"
            local orig="$dir/${1}_bench_azar_merged_100.sorted.bam"
            if [[ -f "$m5" ]]; then echo "$m5"
            elif [[ -f "$orig" ]]; then echo "$orig"
            else
                echo "missing 300x BAM for $1: $m5 (or $orig)" >&2
                return 1
            fi ;;
        *) echo "unknown sample: $1" >&2; return 1 ;;
    esac
}
bed_file() {
    case "$1" in
        HG002) echo "HG002_bench_azar_merged_100.bed" ;;
        HG003) echo "HG003_bench_azar_merged_100.bed" ;;
        HG004) echo "HG004_bench_azar_merged_100.bed" ;;
        *) echo "unknown sample: $1" >&2; return 1 ;;
    esac
}

# --- binary discovery (container build first, then host) -------------------
discover_bin() {
    if [[ -z "${NG_BIN:-}" ]]; then
        local candidate
        for candidate in \
            "$PROJECT_ROOT/target-container/release/pop_var_caller_exp" \
            "$PROJECT_ROOT/target/release/pop_var_caller_exp"; do
            if [[ -x "$candidate" ]] && "$candidate" --version >/dev/null 2>&1; then
                NG_BIN="$candidate"
                break
            fi
        done
    fi
    if [[ -z "${NG_BIN:-}" || ! -x "${NG_BIN}" ]]; then
        echo "no pop_var_caller_exp binary found." >&2
        echo "build with: ./scripts/dev.sh cargo build --release --bin pop_var_caller_exp" >&2
        echo "or set NG_BIN=<path>" >&2
        exit 1
    fi
}

preflight() {
    local f
    for f in "$@"; do
        [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
    done
}

run() {
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        printf 'DRY-RUN:'; printf ' %q' "$@"; printf '\n'
        return 0
    fi
    "$@"
}

record_count() {
    local vcf="$1"
    [[ -f "$vcf" ]] || { echo "?"; return; }
    grep -vc '^#' "$vcf" || true
}

# header lines pass through; data rows kept iff QUAL (col 6) is numeric and
# >= min. `$6 + 0` coerces; non-numeric ('.') falls to 0 and is dropped.
FILTER_AWK='/^#/ { print; next } $6 != "." && ($6 + 0) >= min { print }'

discover_bin
preflight "$REFERENCE" "${REFERENCE}.fai"
mkdir -p "$OUT_DIR"

# Build the catalog once if it is missing. On GRCh38 this is minutes, and every
# later run of every coverage tier reads the file instead of rescanning.
if [[ ! -f "$NG_CATALOG" && "${DRY_RUN:-0}" != "1" ]]; then
    echo "no repeat catalog at $NG_CATALOG — building it (once)"
    "$NG_BIN" repeat-catalog --reference "$REFERENCE" --output "$NG_CATALOG"
    echo
fi

echo "binary    : $NG_BIN"
echo "dataset   : per_sample / $COVERAGE"
echo "reference : $REFERENCE"
echo "catalog   : $NG_CATALOG"
echo "model     : --defaults (nothing fitted; no calibration, contamination or inbreeding)"
echo "min QUAL  : $MIN_QUAL (gated output; raw kept alongside)"
echo "samples   : ${SAMPLES[*]}"
echo "out dir   : $OUT_DIR"
echo

for sample in "${SAMPLES[@]}"; do
    aln="$(align_file "$sample")"
    bed="$BED_DIR/$(bed_file "$sample")"
    raw_vcf="$OUT_DIR/${sample}.raw.vcf"
    vcf="$OUT_DIR/${sample}.vcf"
    log="$OUT_DIR/${sample}.log"

    preflight "$aln" "$bed"

    echo "=================================================================="
    echo "sample    : $sample"
    echo "alignment : $aln"
    echo "regions   : $bed ($(wc -l < "$bed") intervals)"
    echo "raw vcf   : $raw_vcf"
    echo "vcf       : $vcf (QUAL >= $MIN_QUAL)"
    echo

    t0=$(date +%s)
    # ng prints its run report — what it called, and what ground it declined to
    # speak for — on stdout, so both streams go to the log and to the terminal.
    echo "[ng] $sample -> $raw_vcf"
    run "$NG_BIN" call-from-alignments \
        --reference "$REFERENCE" \
        --catalog "$NG_CATALOG" \
        --alignment "$aln" \
        --regions "$bed" \
        --output "$raw_vcf" \
        --defaults 2>&1 | tee "$log"
    t1=$(date +%s)

    [[ "${DRY_RUN:-0}" == "1" ]] && { echo; continue; }

    # Apply the QUAL gate to produce the comparable headline set.
    awk -v min="$MIN_QUAL" "$FILTER_AWK" "$raw_vcf" > "$vcf"

    echo
    echo "elapsed               : $((t1 - t0)) s"
    echo "records (raw)         : $(record_count "$raw_vcf")"
    echo "records (QUAL >= $MIN_QUAL) : $(record_count "$vcf")"
    echo
done

echo "done. VCFs under: $OUT_DIR"
