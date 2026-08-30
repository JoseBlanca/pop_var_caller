#!/usr/bin/env bash
#
# peak_rss.sh — run a command and print its peak resident memory and wall time.
#
#   scripts/peak_rss.sh <result-file> <command> [args...]
#
# Writes one line to <result-file>:   <peak_rss_mb> <wall_seconds>
# The command's own stdout and stderr pass through untouched, so the caller can
# redirect them as it likes — which is why the measurement goes to a file rather
# than to stdout, where it would be interleaved with the command's own output.
#
# Three ways of getting the number, because no single one is available
# everywhere this project runs:
#
#   GNU time (-v)        Linux with the `time` package installed. Reports
#                        kilobytes.
#   BSD time (-l)        macOS. Reports bytes. Reading one as the other is a
#                        thousand-fold error in the only number this produces,
#                        so the two are never mixed.
#   /proc/PID/status     the fallback, and the one that actually gets used
#                        inside this project's dev container, which ships no
#                        time(1) at all. `VmHWM` is the kernel's own high-water
#                        mark for resident memory: it only ever rises, so
#                        sampling it and keeping the largest reading gives the
#                        true peak even at a slow sampling rate. The last
#                        sample must land before the process exits, hence the
#                        20 ms interval.

set -uo pipefail

if (( $# < 2 )); then
    echo "usage: peak_rss.sh <result-file> <command> [args...]" >&2
    exit 2
fi
RESULT_FILE="$1"; shift

start=$(date +%s)

if /usr/bin/time -v true 2>/dev/null; then
    log=$(mktemp)
    /usr/bin/time -v -o "$log" "$@"
    status=$?
    rss_mb=$(awk '/Maximum resident set size/ {printf "%.1f", $NF / 1024; exit}' "$log")
    rm -f "$log"

elif /usr/bin/time -l true 2>/dev/null; then
    log=$(mktemp)
    { /usr/bin/time -l "$@" ; } 2>"$log"
    status=$?
    rss_mb=$(awk '/maximum resident set size/ {printf "%.1f", $1 / 1048576; exit}' "$log")
    rm -f "$log"

elif [[ -r /proc/self/status ]]; then
    "$@" &
    pid=$!
    peak_kb=0
    while kill -0 "$pid" 2>/dev/null; do
        # VmHWM is the peak, not the current size, so a missed sample costs
        # nothing as long as one lands while the process is alive.
        hwm=$(awk '/^VmHWM:/ {print $2; exit}' "/proc/$pid/status" 2>/dev/null || true)
        if [[ -n "${hwm:-}" ]] && (( hwm > peak_kb )); then
            peak_kb=$hwm
        fi
        # `sleep` here is the shell's, not a subshell of the measured process.
        sleep 0.02
    done
    wait "$pid"
    status=$?
    rss_mb=$(awk -v kb="$peak_kb" 'BEGIN {printf "%.1f", kb / 1024}')

else
    echo "peak_rss.sh: no way to measure peak memory on this machine" >&2
    exit 1
fi

end=$(date +%s)
echo "${rss_mb:-0} $((end - start))" > "$RESULT_FILE"
exit "$status"
