#!/usr/bin/env bash
# Rebuild the project's development container image.
#
# `dev.sh` and `claude.sh` only build the image when it doesn't already exist,
# so changes to the Containerfile (or to Cargo.toml/Cargo.lock, which the
# image pre-fetches) won't be picked up automatically. Run this script
# whenever you've edited any of those and want the next `dev.sh` /
# `claude.sh` invocation to use the updated image.
#
# Usage:
#   ./scripts/rebuild-image.sh              # rebuild, reusing layer cache
#   ./scripts/rebuild-image.sh --no-cache   # full rebuild from scratch
#   ./scripts/rebuild-image.sh --pull       # also refresh the FROM base image
#
# Any extra arguments are forwarded to the container runtime's `build`.
#
# # Watching a rebuild, and why that needed fixing (2026-08-15)
#
# **The build is slow enough that you have to be able to see where it is.** A
# rebuild after a base-image change ran for over an hour with nothing on screen,
# and the reason it was invisible was the caller: `rebuild-image.sh | tail -40`
# holds every line until the pipeline ends, which is exactly when the output
# stops being useful. **Do not pipe this script into `tail`, `head` or `grep`.**
#
# So it now does two things itself. It asks for `--progress plain`, which
# streams one line a step instead of the `auto` default's TTY-oriented display
# that collapses to nothing when stdout is not a terminal. And it tees
# everything to a log under `tmp/`, so a build started in one shell can be
# followed from another:
#
#     tail -f tmp/rebuild-image.log
#
# # If it is slow, look at the builder before you look at the Containerfile
#
# Apple's `container` runs builds inside a **separate long-lived builder VM**,
# created once with its own defaults — measured 2026-08-15 as **2 CPUs and
# 2 GB**, on an 18-core, 64 GB host. `dev.sh`'s `--cpus` / `--memory` do not
# reach it: those are `container run` flags and this is a `container build`.
# Three layers here compile Rust from source (`flamegraph`, `samply`,
# `cargo-show-asm`) and a fourth fetches the project's whole crate graph, so two
# cores is the whole story of a multi-hour rebuild. To resize it:
#
#     container builder stop && container builder delete
#     container builder start --cpus 8 --memory 16g
#
# That is a one-time change to the machine, not to this repository, which is why
# it is written here rather than done here.
set -euo pipefail

IMAGE="${IMAGE:-pop-var-caller-dev}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
LOG="${REBUILD_LOG:-$PROJECT_DIR/tmp/rebuild-image.log}"

# Pick the same runtime dev.sh does: podman (Linux dev box) or Apple's
# `container` CLI (macOS). Both accept `build -t … -f … <context>`, and both
# take `--progress plain`.
if command -v podman >/dev/null 2>&1; then
    RUNTIME=podman
elif command -v container >/dev/null 2>&1; then
    RUNTIME=container
else
    echo "Error: neither 'podman' nor 'container' (Apple) found on PATH." >&2
    exit 1
fi

mkdir -p "$(dirname "$LOG")"
echo "rebuilding $IMAGE with $RUNTIME; following along: tail -f $LOG" >&2
if [ "$RUNTIME" = container ]; then
    echo "builder: $(container builder status 2>/dev/null | tail -1)" >&2
fi
started=$(date +%s)

# **Every line is stamped with how long the build has been going.** `--progress
# plain` says which step is running but not for how long, and "which layer is
# the hour" is the only question a slow rebuild raises.
#
# The stamping is a bash read loop rather than `awk`, and rather than `ts`: `ts`
# is not installed here, and macOS ships the one-true-awk, whose `systime()`
# does not exist — reaching for it killed a build at four seconds with
# `awk: calling undefined function systime`. `SECONDS` costs no process a line,
# where `$(date)` in the loop would fork one, and `printf` in bash flushes as it
# goes so `tail -f` is live rather than a block at the end.
#
# `exec` is gone on purpose: the stamping and the tee have to outlive the build
# to flush the log, and `PIPESTATUS` is what carries the build's exit code past
# them.
set +e
"$RUNTIME" build \
    --progress plain \
    -t "$IMAGE" \
    -f "$PROJECT_DIR/Containerfile" \
    "$@" \
    "$PROJECT_DIR" 2>&1 \
    | { SECONDS=0; while IFS= read -r line; do printf '[%5ds] %s\n' "$SECONDS" "$line"; done; } \
    | tee "$LOG"
status=${PIPESTATUS[0]}
set -e
echo "build finished in $(( $(date +%s) - started ))s (exit $status); log in $LOG" >&2
if [ "$status" -ne 0 ]; then
    echo "the slowest steps: grep for the largest jumps in the [   Ns] column" >&2
fi
exit "$status"
