#!/usr/bin/env bash
# Rerun pvc's benchmarks and regenerate the figures quoted about it.
#
#   benchmarks/suite/run_suite.sh                     everything not yet done, pvc at this commit
#   benchmarks/suite/run_suite.sh --stop              finish the running step, then stop
#   benchmarks/suite/run_suite.sh --status            steps done and the next one, per benchmark
#   benchmarks/suite/run_suite.sh --only giab_panel --depth 30x --caller pvc
#   benchmarks/suite/run_suite.sh --rerun-competitors --only hg002_tr --caller hipstr
#   benchmarks/suite/run_suite.sh --slice --only giab_panel        the small test configuration
#   benchmarks/suite/run_suite.sh --summary-only      rebuild results.tsv and summary.md
#
# Other options: --threads N (default 6, the laptop's fast cores), --rescore, --redo-step NAME,
# --allow-dirty (run with uncommitted changes; the commit is recorded as <sha>-dirty),
# --data-root DIR (the checkout holding benchmarks/ and its data; from a worktree, the main
# checkout, found automatically).
#
# This part runs on the host. It records the machine, samples the host's load every 30 s for as
# long as the suite runs, and starts the suite inside the dev container, whose virtual machine gets
# as many cores as --threads and BENCH_MEM of memory (default 32g). Results go to
# benchmark_results/ in this checkout; see benchmarks/suite/README.md.
#
# Stopping: --stop from another terminal ends the run after the step that is running; Ctrl-C ends
# it at once and that step is redone next time. Either way, rerunning the same command resumes.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/../.." && pwd)"
RESULTS="$PROJECT_DIR/benchmark_results"

threads=6
allow_dirty=0
data_root=""
pass=()
slice=0
stop=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stop) stop=1 ;;
        --allow-dirty) allow_dirty=1 ;;
        --threads) threads="$2"; pass+=(--threads "$2"); shift ;;
        --data-root) data_root="$2"; shift ;;
        --slice) slice=1; pass+=(--slice) ;;
        *) pass+=("$1") ;;
    esac
    shift
done
[[ $slice == 1 ]] && RESULTS="$RESULTS/slice"
mkdir -p "$RESULTS"

if [[ $stop == 1 ]]; then
    touch "$RESULTS/STOP"
    echo "asked the running suite to stop after its current step ($RESULTS/STOP)"
    exit 0
fi
rm -f "$RESULTS/STOP"

# The data live in the checkout that has them; a worktree has only the tracked files.
if [[ -z "$data_root" ]]; then
    if [[ -d "$PROJECT_DIR/benchmarks/giab/per_sample/bam" ]]; then
        data_root="$PROJECT_DIR"
    else
        data_root="$(git -C "$PROJECT_DIR" worktree list --porcelain | awk '/^worktree /{print $2; exit}')"
    fi
fi
[[ -d "$data_root/benchmarks" ]] || { echo "no benchmarks/ under --data-root $data_root" >&2; exit 1; }

commit="$(git -C "$PROJECT_DIR" rev-parse --short=8 HEAD)"
if [[ -n "$(git -C "$PROJECT_DIR" status --porcelain --untracked-files=no)" ]]; then
    if [[ $allow_dirty == 0 ]]; then
        echo "the checkout has uncommitted changes; commit them, or pass --allow-dirty" >&2
        exit 1
    fi
    commit="$commit-dirty"
fi

mem="${BENCH_MEM:-32g}"
model="$(sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -m)"
cores="$(sysctl -n hw.ncpu 2>/dev/null || nproc)"
fast="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo '?')"
slow="$(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null || echo '?')"
ram_gb=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1073741824 ))
machine="$model, $cores cores ($fast fast, $slow slower), ${ram_gb} GB; Linux VM with $threads cores, $mem"

# vcfeval, pinned; not in the dev image.
"$PROJECT_DIR/benchmarks/suite/fetch_tools.sh"

# The host's load, every 30 s, for as long as the suite runs. The virtual machine cannot see it.
(
    while true; do
        load="$(sysctl -n vm.loadavg 2>/dev/null | awk '{print $2}')"
        [[ -z "$load" ]] && load="$(awk '{print $1}' /proc/loadavg)"
        echo "$(date +%s) $load" >> "$RESULTS/host_load.tsv"
        sleep 30
    done
) &
sampler=$!
trap 'kill $sampler 2>/dev/null || true' EXIT

extra=()
[[ "$data_root" != "$PROJECT_DIR" ]] && extra=(DEV_EXTRA_MOUNT="$data_root")

status=0
env DEV_CPUS="$threads" DEV_MEM="$mem" ${extra[@]+"${extra[@]}"} \
    "$PROJECT_DIR/scripts/dev.sh" python3 "$PROJECT_DIR/benchmarks/suite/suite.py" \
    --data-root "$data_root" --pvc-commit "$commit" --machine "$machine" \
    ${pass[@]+"${pass[@]}"} || status=$?
exit $status
