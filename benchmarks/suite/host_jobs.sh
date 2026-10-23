#!/usr/bin/env bash
# The host job runner: runs DeepVariant and GLnexus in their own official images for the suite.
#
#   benchmarks/suite/host_jobs.sh <results-dir> <project-dir> <data-root> <genomes-dir> [<extra-dir>]
#
# run_suite.sh starts it in the background and stops it on exit. The suite itself runs inside the
# dev container, which cannot start other containers; for a DeepVariant or GLnexus step it writes a
# request into <results-dir>/host_jobs/ and waits, and this runner executes it on the host with the
# host's own podman or Docker, then writes the result beside the request.
#
# It runs only the two images below, only with the project, the data root and the genomes directory
# mounted (the data and genomes read-only), and only with a working directory inside the project.
# Inside the image the command runs under the suite's own measuring (lib/measure_cli.py), so its
# time and memory are measured as every other caller's are.
#
# **Where it can run them:** an x86-64 machine with podman or Docker. Anywhere else (this Mac:
# arm64, and Apple's `container` is not used for x86 images) it writes `capabilities.json` saying
# why not, answers every request with that reason, and the suite skips those steps.
set -uo pipefail

results=$1
project=$2
data_root=$3
genomes=$4
extra_dir=${5:-}
memory="${BENCH_MEM:-32g}"
jobs="$results/host_jobs"
mkdir -p "$jobs"

image_of() {  # the only images this runner will run (bash 3.2 has no associative arrays)
    case "$1" in
        deepvariant) echo "docker.io/google/deepvariant:1.10.0" ;;
        glnexus) echo "ghcr.io/dnanexus-rnd/glnexus:v1.4.1" ;;
        *) echo "" ;;
    esac
}

arch="$(uname -m)"
runtime=""
if [[ "$arch" == x86_64 ]]; then
    if command -v podman > /dev/null; then runtime=podman
    elif command -v docker > /dev/null; then runtime=docker
    fi
fi
if [[ -n "$runtime" ]]; then
    available=true
    reason="x86-64 with $runtime"
elif [[ "$arch" != x86_64 ]]; then
    available=false
    reason="this machine is $arch; DeepVariant's images are built for x86-64 only"
else
    available=false
    reason="no podman or Docker on this x86-64 machine"
fi
printf '{"available": %s, "reason": "%s", "arch": "%s", "runtime": "%s", "images": {"deepvariant": "%s", "glnexus": "%s"}}\n' \
    "$available" "$reason" "$arch" "$runtime" "$(image_of deepvariant)" "$(image_of glnexus)" \
    > "$jobs/capabilities.json"

running=""
cleanup() { [[ -n "$running" && -n "$runtime" ]] && "$runtime" rm -f "$running" > /dev/null 2>&1; exit 0; }
trap cleanup INT TERM EXIT

field() {  # field <request.json> <key>
    python3 -c 'import json,sys; v=json.load(open(sys.argv[1]))[sys.argv[2]]; print(v)' "$1" "$2"
}

inside() {  # is the path inside the project?
    case "$(cd "$1" 2>/dev/null && pwd)/" in "$project"/*) return 0 ;; *) return 1 ;; esac
}

while true; do
    for request in "$jobs"/*.request; do
        [[ -e "$request" ]] || continue
        id="$(basename "$request" .request)"
        answer="$jobs/$id.answer"
        if [[ "$available" != true ]]; then
            echo "{\"ok\": false, \"unavailable\": true, \"error\": \"$reason\"}" > "$answer.tmp"
            mv "$answer.tmp" "$answer"; rm -f "$request"; continue
        fi
        image_key="$(field "$request" image)"
        image="$(image_of "$image_key")"
        workdir="$(field "$request" workdir)"
        result="$(field "$request" result)"
        threads="$(field "$request" threads)"
        label="$(field "$request" label)"
        command="$(field "$request" command)"
        if [[ -z "$image" ]] || ! inside "$workdir"; then
            echo "{\"ok\": false, \"error\": \"refused: image '$image_key' or working directory outside the project\"}" > "$answer.tmp"
            mv "$answer.tmp" "$answer"; rm -f "$request"; continue
        fi
        user=()
        [[ "$runtime" == docker ]] && user=(--user "$(id -u):$(id -g)")
        mounts=(-v "$project:$project" -v "$data_root:$data_root:ro")
        [[ -d "$genomes" && "$genomes" != "$data_root" ]] && mounts+=(-v "$genomes:$genomes:ro")
        [[ -n "$extra_dir" && -d "$extra_dir" ]] && mounts+=(-v "$extra_dir:$extra_dir:ro")
        running="pvcsuite-$id"
        # Measured by the suite's own code where the image has Python; otherwise run as it is, and
        # the step records no measurement (GLnexus's image may lack Python).
        mkdir -p "$workdir/logs"
        "$runtime" run --rm --name "$running" --cpus "$threads" --memory "$memory" ${user[@]+"${user[@]}"} \
            "${mounts[@]}" -w "$workdir" \
            -e MCLI="$project/benchmarks/suite/lib/measure_cli.py" -e RESULT="$result" \
            -e LOGDIR="$workdir/logs" -e LABEL="$label" -e CMD="$command" \
            "$image" sh -c 'if command -v python3 > /dev/null 2>&1; then
                    exec python3 "$MCLI" "$RESULT" "$LOGDIR" "$LABEL" -- "$CMD"
                else
                    sh -c "$CMD" > "$LOGDIR/$LABEL.stdout" 2> "$LOGDIR/$LABEL.stderr"
                fi' > "$jobs/$id.runner.log" 2>&1
        status=$?
        running=""
        echo "{\"ok\": $([[ $status == 0 ]] && echo true || echo false), \"exit\": $status, \"runner_log\": \"$jobs/$id.runner.log\"}" > "$answer.tmp"
        mv "$answer.tmp" "$answer"
        rm -f "$request"
    done
    sleep 2
done
