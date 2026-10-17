#!/usr/bin/env bash
# Call copied benchmark runs with one build, under a label, and score every label so far.
#
#   scripts/ng_str_accuracy/call_and_score.sh <label> <pvc-binary> <run>...
#
# Run it inside the dev container. For each run: makes the psps (once per PSP_DIR — default
# `psps`; give each build that changes how reads are measured its own, e.g. PSP_DIR=psps_<label>),
# calls `call-from-psps` with the simulator's TRUE stutter (`true_stutter.parameters.toml`, one
# stutter group per library — not a fit), with `--explain-loci` over the segments, writes
# `<label>.vcf`, then scores GATK, main-as-of-4b528afe and every `x_*` label in the run with
# `score_callers.py`. **Label your builds `x_<something>`**: only those are scored.
#
# Copy a build's binary somewhere stable first (e.g. tmp/str_accuracy/bin/<label>) so a later
# `cargo build` cannot change what an unfinished run is measuring.
#
# Env — set inside the container, e.g. `scripts/dev.sh env PSP_DIR=psps_x <this script> ...`, since
# dev.sh does not pass the host environment through: PSP_DIR, STR_RUNS_DIR (default
# <repo>/tmp/str_accuracy/runs), THREADS (default 8), and
# anything the binary itself reads (an experiment's switch, say).
set -euo pipefail
root=$(cd "$(dirname "$0")/../.." && pwd)
runs=${STR_RUNS_DIR:-$root/tmp/str_accuracy/runs}
label=${1:?usage: $0 <label> <binary> <run>...}; ng=${2:?binary}; shift 2
psps=${PSP_DIR:-psps}
for run in "$@"; do
    out=$runs/$run
    if [[ ! -f $out/$psps.done ]]; then
        args=()
        for bam in "$out"/bams/*.bam; do args+=(--alignment "$bam"); done
        "$ng" generate-psps --reference "$out/reference.fa" --catalog "$out/reference.fa.repeats.parquet" \
            "${args[@]}" --regions "$out/confident.bed" --output-dir "$out/$psps" --force > "$out/$psps.log" 2>&1
        touch "$out/$psps.done"
    fi
    start=$SECONDS
    "$ng" call-from-psps --reference "$out/reference.fa" --catalog "$out/reference.fa.repeats.parquet" \
        --psp "$out/$psps" --parameters "$out/true_stutter.parameters.toml" --threads "${THREADS:-8}" \
        --explain-loci "$out/segments.bed" "$out/$label.explain.tsv" \
        --output "$out/$label.vcf" --log "$out/$label.pvc.log" --force > "$out/$label.log" 2>&1
    echo "$run $label called in $((SECONDS - start)) s"
    ls "$out"/*.vcf | xargs -n1 basename | sed 's/\.vcf$//' | grep -E '^(x_|ng_main_true_stutter$)' > "$out/callers.list"
    SKIP_VCFEVAL=1 python3 "$root/scripts/ng_str_accuracy/score_callers.py" "$out"
done
