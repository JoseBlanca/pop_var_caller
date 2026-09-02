#!/usr/bin/env bash
# ng runner — shared across benchmarks.
#
#   benchmarks/lib/run_ng.sh <bench.config.sh> [single|cohort]
#
# single : one alignment file -> single-sample VCF, restricted to the
#          benchmark BED.
# cohort : every alignment file -> one joint VCF in one invocation. There is
#          no per-sample intermediate — ng holds every sample's file open and
#          walks them all at one frontier — so unlike run_ours.sh there is no
#          .psp stage to skip on a re-run.
#
# ng is the experimental caller (`pop_var_caller_exp call-from-alignments`).
# Three things about it shape this script:
#
#   * It needs a tandem-repeat catalog built from the same reference, and
#     refuses one built from another. The default path is
#     `<reference>.repeats.parquet`; where the reference sits on a read-only
#     mount, set NG_CATALOG to somewhere writable. The catalog is built once
#     if missing (about 100 s on GRCh38, at one thread).
#   * It needs either a fitted parameters file or `--defaults`. No command
#     writes a fitted file yet, so the default here is `--defaults`: no
#     base-quality calibration, no contamination, no inbreeding.
#   * It does NOT call inside tandem repeats. Every repeat tract in the BED is
#     counted as ground it cannot speak for and named in the run report, which
#     this script tees to the log — read it before reading the recall.
#
# Unlike run_ours.sh, ng takes a BED, so it is restricted to the benchmark
# regions exactly as GATK and freebayes are.
#
# A post-call QUAL >= MIN_QUAL filter is applied inline with awk, matching what
# run_freebayes.sh does, so a benchmark's fairness setting applies to ng too.
#
# Env overrides (see common.sh for the rest):
#   NG_BIN        binary (default: auto-detect container then host build)
#   NG_CATALOG    repeat catalog (default: <reference>.repeats.parquet)
#   NG_PARAMETERS fitted parameters file; without it the run is --defaults
#   REFERENCE     FASTA (.fai sibling required)
#   PLOIDY        --ploidy (default 2)
#   MIN_QUAL      QUAL floor applied pre-write (default 30)
#   EXTRA_ARGS    appended verbatim to the ng command line
#   DRY_RUN=1     print commands instead of running them

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

CONFIG="${1:-}"
MODE="${2:-single}"
bench_load_config "$CONFIG"

EXTRA_ARGS=${EXTRA_ARGS:-}
OUT_DIR="$OUT_ROOT/ng"
NG_CATALOG="${NG_CATALOG:-${REFERENCE}.repeats.parquet}"

# --- binary discovery ------------------------------------------------------
# Mirrors bench_discover_ours_bin: container build first (canonical per
# CLAUDE.md), then host build, verifying each actually runs here — a Linux ELF
# under target-container/ is +x but unusable on a macOS host.
discover_ng_bin() {
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

ensure_catalog() {
    [[ -f "$NG_CATALOG" ]] && return 0
    if [[ "${DRY_RUN:-0}" == "1" ]]; then
        echo "(dry run: would build the repeat catalog at $NG_CATALOG)"
        return 0
    fi
    echo "no repeat catalog at $NG_CATALOG — building it (once)"
    "$NG_BIN" repeat-catalog --reference "$REFERENCE" --output "$NG_CATALOG"
    echo
}

# The model flag pair: a fitted file if one was named, the compiled-in defaults
# otherwise. They are mutually exclusive at the CLI.
model_args() {
    if [[ -n "${NG_PARAMETERS:-}" ]]; then
        printf '%s\n%s\n' --parameters "$NG_PARAMETERS"
    else
        printf '%s\n' --defaults
    fi
}
model_label() {
    if [[ -n "${NG_PARAMETERS:-}" ]]; then echo "$NG_PARAMETERS"
    else echo "--defaults (nothing fitted)"; fi
}

# header lines pass through; data rows kept iff QUAL (col 6) is numeric and
# >= min. `$6 + 0` coerces; non-numeric falls to 0 and is dropped.
FILTER_AWK='/^#/ { print; next } $6 != "." && ($6 + 0) >= min { print }'

# ng writes its VCF to a path it is given rather than to stdout, so the QUAL
# gate is a second pass over that file. The ungated file is kept beside the
# gated one under a `.raw.vcf` name, the way run_freebayes_per_sample.sh does,
# so the low-QUAL tail can still be inspected.
gate_output() {
    local raw="$1" out="$2"
    awk -v min="$MIN_QUAL" "$FILTER_AWK" "$raw" > "$out"
}

run_single() {
    local cram base raw_vcf out_vcf log
    cram="$(bench_single_cram)"
    base="$(bench_sample_base "$cram")"
    bench_preflight "$cram" "$BED" "$REFERENCE" "${REFERENCE}.fai"

    raw_vcf="$OUT_DIR/single_${base}.raw.vcf"
    out_vcf="$OUT_DIR/single_${base}.vcf"
    log="$OUT_DIR/single_${base}.log"
    mkdir -p "$OUT_DIR"
    ensure_catalog

    echo "binary    : $NG_BIN"
    echo "mode      : single ($BENCH_NAME) — alignments -> VCF in one process"
    echo "sample    : $base"
    echo "input     : $cram"
    echo "reference : $REFERENCE"
    echo "catalog   : $NG_CATALOG"
    echo "model     : $(model_label)"
    echo "regions   : $BED ($(wc -l < "$BED") intervals)"
    echo "ploidy    : $PLOIDY"
    echo "min QUAL  : $MIN_QUAL (gated output; raw kept as $(basename "$raw_vcf"))"
    echo "output    : $out_vcf"
    echo

    local t0 t1
    t0=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "" -- "$NG_BIN" call-from-alignments \
        --reference "$REFERENCE" \
        --catalog "$NG_CATALOG" \
        --alignment "$cram" \
        --regions "$BED" \
        --ploidy "$PLOIDY" \
        --output "$raw_vcf" \
        $(model_args) \
        $EXTRA_ARGS 2>&1 | tee "$log"
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0

    gate_output "$raw_vcf" "$out_vcf"
    echo
    echo "elapsed: $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf") (QUAL >= $MIN_QUAL) of $(bench_record_count "$raw_vcf") written"
}

run_cohort() {
    bench_list_crams
    bench_preflight "$BED" "$REFERENCE" "${REFERENCE}.fai"

    local raw_vcf="$OUT_DIR/cohort.raw.vcf"
    local out_vcf="$OUT_DIR/cohort.vcf"
    local log="$OUT_DIR/cohort.log"
    mkdir -p "$OUT_DIR"
    ensure_catalog

    # One --alignment flag per sample.
    local aln_args=() c
    for c in "${BENCH_CRAMS[@]}"; do
        aln_args+=(--alignment "$c")
    done

    echo "binary    : $NG_BIN"
    echo "mode      : cohort ($BENCH_NAME) — every sample in one process"
    echo "samples   : ${#BENCH_CRAMS[@]}"
    echo "reference : $REFERENCE"
    echo "catalog   : $NG_CATALOG"
    echo "model     : $(model_label)"
    echo "regions   : $BED ($(wc -l < "$BED") intervals)"
    echo "ploidy    : $PLOIDY"
    echo "min QUAL  : $MIN_QUAL (gated output; raw kept as $(basename "$raw_vcf"))"
    echo "output    : $out_vcf"
    echo

    local t0 t1
    t0=$(bench_now)
    # shellcheck disable=SC2086
    bench_run "" -- "$NG_BIN" call-from-alignments \
        --reference "$REFERENCE" \
        --catalog "$NG_CATALOG" \
        "${aln_args[@]}" \
        --regions "$BED" \
        --ploidy "$PLOIDY" \
        --output "$raw_vcf" \
        $(model_args) \
        $EXTRA_ARGS 2>&1 | tee "$log"
    t1=$(bench_now)
    [[ "${DRY_RUN:-0}" == "1" ]] && return 0

    gate_output "$raw_vcf" "$out_vcf"
    echo
    echo "elapsed: $((t1 - t0)) s"
    echo "records: $(bench_record_count "$out_vcf") (QUAL >= $MIN_QUAL) of $(bench_record_count "$raw_vcf") written"
    echo "samples in vcf: $(bench_sample_count "$out_vcf")"
}

discover_ng_bin

case "$MODE" in
    single) run_single ;;
    cohort) run_cohort ;;
    *) echo "unknown mode: $MODE (expected single|cohort)" >&2; exit 2 ;;
esac
