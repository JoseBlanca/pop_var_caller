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
        --stop-at)
            # The next time the clock reads HH:MM, local time (the virtual machine runs on UTC).
            at="$(date -j -f "%H:%M:%S" "$2:00" +%s 2>/dev/null)" || { echo "--stop-at takes HH:MM" >&2; exit 1; }
            [[ "$at" -le "$(date +%s)" ]] && at=$(( at + 86400 ))
            pass+=(--stop-at-epoch "$at"); echo "no step will start after $(date -r "$at")"; shift ;;
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

# The host's load every 5 s, for as long as the suite runs; the virtual machine cannot see it. Two
# figures a line: the one-minute load average, and the cores' worth of CPU used by every process
# except the suite's own virtual machine ("other CPU"). The load average counts the suite too, so
# it is the other CPU that says whether the machine was quiet. The suite's virtual machine is any
# `com.apple.Virtualization.VirtualMachine` process that did not exist before the suite started;
# a container another session starts during the run would be counted as the suite's.
vm_name='com.apple.Virtualization.VirtualMachine'
before_vms=" $( (pgrep -f "$vm_name" || true) | tr "\n" " ") "
(
    while true; do
        load="$(sysctl -n vm.loadavg 2>/dev/null | awk '{print $2}')"
        [[ -z "$load" ]] && load="$(awk '{print $1}' /proc/loadavg)"
        other="$(ps -A -o pid=,%cpu=,comm= | awk -v before="$before_vms" -v vm="$vm_name" -v me=$$ '
            { if (index($3, vm) && index(before, " " $1 " ") == 0) next; t += $2 }
            END { printf "%.2f", t / 100 }')"
        echo "$(date +%s) $load $other" >> "$RESULTS/host_load.tsv"
        sleep 5
    done
) &
sampler=$!
trap 'kill $sampler 2>/dev/null || true' EXIT

extra=()
[[ "$data_root" != "$PROJECT_DIR" ]] && extra=(DEV_EXTRA_MOUNT="$data_root")

status=0
# caffeinate keeps the Mac from sleeping while the suite runs.
caffeinate -is env DEV_CPUS="$threads" DEV_MEM="$mem" ${extra[@]+"${extra[@]}"} \
    "$PROJECT_DIR/scripts/dev.sh" python3 "$PROJECT_DIR/benchmarks/suite/suite.py" \
    --data-root "$data_root" --pvc-commit "$commit" --machine "$machine" \
    ${pass[@]+"${pass[@]}"} || status=$?
exit $status
