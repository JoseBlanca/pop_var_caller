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

threads=""
profile=""
allow_dirty=0
data_root=""
pass=()
slice=0
stop=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stop) stop=1 ;;
        --allow-dirty) allow_dirty=1 ;;
        --threads) threads="$2"; shift ;;
        --profile) profile="$2"; shift ;;
        --data-root) data_root="$2"; shift ;;
        --slice) slice=1; pass+=(--slice) ;;
        --stop-at)
            # The next time the clock reads HH:MM, local time (the virtual machine runs on UTC).
            at="$(date -j -f "%H:%M:%S" "$2:00" +%s 2>/dev/null || date -d "$2" +%s 2>/dev/null)" \
                || { echo "--stop-at takes HH:MM" >&2; exit 1; }
            [[ "$at" -le "$(date +%s)" ]] && at=$(( at + 86400 ))
            pass+=(--stop-at-epoch "$at"); echo "no step will start after $(date -d "@$at" 2>/dev/null || date -r "$at")"; shift ;;
        *) pass+=("$1") ;;
    esac
    shift
done
# A machine profile sets what differs between machines: cores and memory for the suite, where the
# tomato CRAMs and their reference are, one more directory to mount. benchmarks/suite/profiles/
# <name>.sh, chosen by --profile or else by this machine's short hostname; none on the laptop.
[[ -z "$profile" ]] && profile="$(hostname -s 2>/dev/null || hostname)"
profile_file="$PROJECT_DIR/benchmarks/suite/profiles/$profile.sh"
if [[ -f "$profile_file" ]]; then
    # shellcheck source=/dev/null
    source "$profile_file"
    echo "machine profile: $profile_file"
fi
threads="${threads:-${BENCH_THREADS:-6}}"
pass+=(--threads "$threads")
[[ -n "${BENCH_TOMATO_CRAMS:-}" ]] && pass+=(--tomato-crams "$BENCH_TOMATO_CRAMS")
[[ -n "${BENCH_TOMATO_REFERENCE:-}" ]] && pass+=(--tomato-reference "$BENCH_TOMATO_REFERENCE")
[[ -n "${BENCH_TOMATO_SIZES:-}" ]] && pass+=(--tomato-sizes "$BENCH_TOMATO_SIZES")
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

# pvc's results are labelled by the last commit that changed pvc's own code (the crate, its
# vendored libraries and the examples the simulated benchmark runs), not by the branch's head, so
# a change to the suite's scripts alone does not make pvc look new and rerun.
pvc_paths=(src Cargo.toml Cargo.lock vendor examples build.rs)
commit="$(git -C "$PROJECT_DIR" log -1 --format=%h --abbrev=8 -- "${pvc_paths[@]}")"
if [[ -n "$(git -C "$PROJECT_DIR" status --porcelain --untracked-files=no -- "${pvc_paths[@]}")" ]]; then
    if [[ $allow_dirty == 0 ]]; then
        echo "the checkout has uncommitted changes; commit them, or pass --allow-dirty" >&2
        exit 1
    fi
    commit="$commit-dirty"
fi

mem="${BENCH_MEM:-32g}"
if [[ "$(uname -s)" == Darwin ]]; then
    model="$(sysctl -n machdep.cpu.brand_string)"
    cores="$(sysctl -n hw.ncpu)"
    fast="$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || echo '?')"
    slow="$(sysctl -n hw.perflevel1.physicalcpu 2>/dev/null || echo '?')"
    ram_gb=$(( $(sysctl -n hw.memsize) / 1073741824 ))
    machine="$model, $cores cores ($fast fast, $slow slower), ${ram_gb} GB; Linux VM with $threads cores, $mem"
else
    model="$(awk -F': ' '/^model name/{print $2; exit}' /proc/cpuinfo)"
    cores="$(nproc --all)"
    ram_gb=$(( $(awk '/^MemTotal/{print $2}' /proc/meminfo) / 1048576 ))
    machine="$(hostname -s), $model, $cores logical cores, ${ram_gb} GB; the suite limited to $threads cores and $mem"
fi

# vcfeval, pinned; not in the dev image.
"$PROJECT_DIR/benchmarks/suite/fetch_tools.sh"

# The host's load every 5 s, for as long as the suite runs; the virtual machine cannot see it. Two
# figures a line: the one-minute load average, and the cores' worth of CPU used by every process
# except the suite's own virtual machine ("other CPU"). The load average counts the suite too, so
# it is the other CPU that says whether the machine was quiet. The suite's virtual machine is any
# `com.apple.Virtualization.VirtualMachine` process that did not exist before the suite started;
# a container another session starts during the run would be counted as the suite's.
#
# **On Linux** the containers' processes are ordinary host processes, so "the suite's own" are
# instead every process of this user inside a podman container (a `libpod` cgroup): the dev
# container and the DeepVariant jobs. Other users' containers stay counted as other load.
# `ps`'s %cpu is a lifetime average; the busy figure therefore comes from /proc's tick counters,
# sampled 5 s apart, on Linux.
vm_name='com.apple.Virtualization.VirtualMachine'
before_vms=" $( (pgrep -f "$vm_name" || true) | tr "\n" " ") "
linux_other() {  # cores used by others' processes over the last interval, from /proc ticks
    local me; me="$(id -u)"
    for d in /proc/[0-9]*; do
        local pid="${d#/proc/}"
        [[ -r "$d/stat" ]] || continue
        if [[ "$(stat -c %u "$d" 2>/dev/null)" == "$me" ]] && grep -q libpod "$d/cgroup" 2>/dev/null; then
            continue
        fi
        # utime + stime, fields 14 and 15 after the command name.
        awk -v p="$pid" '{sub(/^.*\) /, ""); print p, $12 + $13}' "$d/stat" 2>/dev/null
    done
}
(
    hz="$(getconf CLK_TCK 2>/dev/null || echo 100)"
    prev=""
    while true; do
        load="$(sysctl -n vm.loadavg 2>/dev/null | awk '{print $2}')"
        [[ -z "$load" ]] && load="$(awk '{print $1}' /proc/loadavg)"
        if [[ "$(uname -s)" == Darwin ]]; then
            other="$(ps -A -o pid=,%cpu=,comm= | awk -v before="$before_vms" -v vm="$vm_name" '
                { if (index($3, vm) && index(before, " " $1 " ") == 0) next; t += $2 }
                END { printf "%.2f", t / 100 }')"
        else
            now="$(linux_other)"
            other="$(awk -v hz="$hz" 'NR==FNR {old[$1]=$2; next} ($1 in old) {t += $2 - old[$1]}
                END {printf "%.2f", t / hz / 5}' <(printf '%s\n' "$prev") <(printf '%s\n' "$now"))"
            prev="$now"
        fi
        echo "$(date +%s) $load $other" >> "$RESULTS/host_load.tsv"
        sleep 5
    done
) &
sampler=$!
trap 'kill $sampler 2>/dev/null || true' EXIT

# DeepVariant and GLnexus run in their own official images, started on this host by the job runner
# (host_jobs.sh), since the dev container cannot start containers; where the machine cannot run
# them it says so and the suite skips those steps.
BENCH_MEM="$mem" "$PROJECT_DIR/benchmarks/suite/host_jobs.sh" "$RESULTS" "$PROJECT_DIR" "$data_root" \
    "${DEV_GENOMES_DIR:-$HOME/genomes}" "${BENCH_EXTRA_MOUNT:-}" &
host_jobs=$!
trap 'kill $sampler $host_jobs 2>/dev/null || true' EXIT

extra=()
if [[ "$data_root" != "$PROJECT_DIR" ]]; then
    extra=(DEV_EXTRA_MOUNT="$data_root")
    [[ -n "${BENCH_EXTRA_MOUNT:-}" ]] && { echo "dev.sh mounts one extra directory: the data root and the profile's BENCH_EXTRA_MOUNT cannot both be mounted" >&2; exit 1; }
elif [[ -n "${BENCH_EXTRA_MOUNT:-}" ]]; then
    extra=(DEV_EXTRA_MOUNT="$BENCH_EXTRA_MOUNT")
fi
# The memory and core limits are applied to podman too (DEV_LIMIT_PODMAN), not only to Apple's VM.
extra+=(DEV_LIMIT_PODMAN=1)
keep_awake=()
command -v caffeinate > /dev/null && keep_awake=(caffeinate -is)

status=0
# caffeinate keeps the Mac from sleeping while the suite runs (macOS only).
${keep_awake[@]+"${keep_awake[@]}"} env DEV_CPUS="$threads" DEV_MEM="$mem" ${extra[@]+"${extra[@]}"} \
    "$PROJECT_DIR/scripts/dev.sh" python3 "$PROJECT_DIR/benchmarks/suite/suite.py" \
    --data-root "$data_root" --pvc-commit "$commit" --machine "$machine" \
    ${pass[@]+"${pass[@]}"} || status=$?
exit $status
