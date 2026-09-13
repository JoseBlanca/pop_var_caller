#!/usr/bin/env bash
# Run the tomato1 performance benchmark suite that feeds perf_dashboard.py
# (freebayes and GATK). The production caller's scripts (`perf_ours_*`) were
# deleted with that caller in promotion Milestone D; the dashboard still reads
# the results they wrote.
#
# Split by where each tool lives:
#   - HOST  (uv run --script): freebayes.
#   - CONTAINER (scripts/dev.sh, reference genome bind-mounted): GATK,
#     which isn't on the host. The genome lives outside the project tree,
#     so it's mounted via DEV_EXTRA_MOUNT.
#
# Sections and the scripts that feed them:
#   §1 one intermediate, 4 threads:      perf_gatk_gvcf_4t
#   §2 scaling CRAM→VCF:                 perf_freebayes
#   §3 scaling intermediate→VCF:         perf_gatk_joint
#
# NOTE: GATK direct multi-sample HaplotypeCaller (CRAM→VCF) is deliberately
# NOT measured — its wall is super-linear in N (re-assembly over pooled deep
# coverage), so it's impractical past a few samples. GATK's scalable route
# is the GVCF flow (§4, perf_gatk_joint).
#
# §4 assumes the per-sample intermediates already exist:
#   - GATK: results/gatk/cohort/gvcf/*.g.vcf.gz
# build them first with the cohort drivers if missing (see the warnings
# this script prints).
#
# Usage:
#   benchmarks/tomato1/scripts/run_perf_all.sh [all|host|container]
#
#   all        run everything (default)
#   host       only the host scripts (freebayes)
#   container  only the GATK scripts (inside scripts/dev.sh)
#
# Env:
#   GENOMES            dir bind-mounted into the container for GATK's
#                      reference (default: $HOME/genomes)
#   SIZES, THREADS, REFERENCE, GATK_BIN, JAVA_HEAP, ... forwarded to the
#                      per-caller scripts (see each script's header).

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
TEST_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
GENOMES="${GENOMES:-$HOME/genomes}"
MODE="${1:-all}"

case "$MODE" in
    all|host|container) ;;
    *) echo "usage: $0 [all|host|container]" >&2; exit 2 ;;
esac

HOST_SCRIPTS=(
    perf_freebayes.py            # §2
)
CONTAINER_SCRIPTS=(
    perf_gatk_gvcf_4t.py         # §1
    perf_gatk_joint.py           # §3
)

declare -a RESULTS=()

record() { RESULTS+=("$1"); }

run_host() {
    local script="$1" sizes_override="${2:-}"
    echo
    echo "=== [host] $script ==="
    local rc=0
    if [[ -n "$sizes_override" ]]; then
        SIZES="$sizes_override" uv run --script "$SCRIPT_DIR/$script" || rc=$?
    else
        uv run --script "$SCRIPT_DIR/$script" || rc=$?
    fi
    if [[ $rc -eq 0 ]]; then record "ok   host/$script"
    else record "FAIL host/$script (exit $rc)"; fi
}

run_container() {
    local script="$1"
    echo
    echo "=== [container] $script ==="
    if DEV_EXTRA_MOUNT="$GENOMES" "$PROJECT_ROOT/scripts/dev.sh" \
            python3 "benchmarks/tomato1/scripts/$script"; then
        record "ok   container/$script"
    else
        record "FAIL container/$script (exit $?)"
    fi
}

# --- Pre-flight warning for §4 GATK (GVCFs must be pre-built) -----------
gvcf_dir="$TEST_DIR/results/gatk/cohort/gvcf"
if [[ "$MODE" != host ]] && ! compgen -G "$gvcf_dir/*.g.vcf.gz" >/dev/null; then
    echo "WARNING: no *.g.vcf.gz in $gvcf_dir — §4 perf_gatk_joint will have no inputs." >&2
    echo "         build them first: DEV_EXTRA_MOUNT=\$HOME/genomes ./scripts/dev.sh \\" >&2
    echo "             benchmarks/tomato1/scripts/build_cohort_intermediates.sh gvcf" >&2
fi

# --- Run -----------------------------------------------------------------
if [[ "$MODE" == all || "$MODE" == host ]]; then
    for s in "${HOST_SCRIPTS[@]}"; do run_host "$s"; done
fi
if [[ "$MODE" == all || "$MODE" == container ]]; then
    for s in "${CONTAINER_SCRIPTS[@]}"; do run_container "$s"; done
fi

# --- Summary -------------------------------------------------------------
echo
echo "=== summary ==="
printf '%s\n' "${RESULTS[@]}"
echo
echo "view the dashboard with:"
echo "  uvx marimo edit --sandbox benchmarks/tomato1/scripts/perf_dashboard.py"

# Non-zero overall exit if anything failed.
if printf '%s\n' "${RESULTS[@]}" | grep -q '^FAIL'; then
    exit 1
fi
