#!/usr/bin/env bash
# The neighbour-evidence pass on the tomato cohort: how much missing data it fills at 100, 200 and
# 300 accessions, and how often a fill is wrong (tomato_fill.py says how).
#
# usage: tomato_fill_test.sh <cohort.tsv> <reference.fa> <work dir> [sizes, default 100,200,300]
#
# cohort.tsv: one line per alignment file, `sample<TAB>path`, in the order accessions join the
# cohort; the first N samples are the cohort of size N, so each cohort holds the smaller ones. A
# sample with two files has two lines. The suite's `inputs/tomato/catalogue.tsv` (sample, project,
# path) cut to columns 1 and 3 is one; it lists samples by name, so shuffle it with a fixed seed
# first for cohorts that mix projects, and with them depths:
#   cut -f1,3 catalogue.tsv | shuf --random-source=<(yes) > cohort.tsv
#
# Psps made already: give each sample one line, `sample<TAB>path.psp`, and CATALOG, the repeat
# catalog they were made with. Nothing is generated then, REGIONS is unused, and calling covers the
# ground the psps cover — the whole genome, if that is what they hold.
#
# Environment: REGIONS (a BED; default all of SL4.0ch01), THREADS (default 16), PVC (the binary;
# default the newer of target/release/pvc and target-container/release/pvc),
# PYTHON (default `uv run python`).
#
# Per sample, generate-psps once, shared by every size. Per size: estimate-parameters, then
# call-from-psps twice, with `--refine-genotypes off` and `on`. How often a fill is wrong comes from
# the fill check the refined run writes beside its VCF (`refined.fill_check.tsv`): it hides
# genotypes the reads settle and scores what the pass would write there. Until 2026-11-03 this script
# hid 2% of the confidently called SNP genotypes in the pass-off VCF and refined that with
# `pvc refine-genotypes`; that command is removed, and a VCF cannot be fed back into calling.
# Every step is skipped when its output is there, so a run that stops resumes.
set -euo pipefail

cohort=$1
reference=$2
work=$3
sizes=${4:-100,200,300}
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../../../../.." && pwd)
threads=${THREADS:-16}
read -r -a python <<< "${PYTHON:-uv run python}"

if [[ -z ${PVC:-} ]]; then
    PVC=""
    for candidate in "$repo/target/release/pvc" "$repo/target-container/release/pvc"; do
        if [[ -x $candidate && ( -z $PVC || $candidate -nt $PVC ) ]]; then
            PVC=$candidate
        fi
    done
fi
[[ -n $PVC ]] || { echo "no pvc binary: cargo build --release" >&2; exit 1; }
echo "pvc: $PVC ($("$PVC" --version))"

mkdir -p "$work/psps" "$work/logs"
if [[ -z ${REGIONS:-} ]]; then
    length=$(awk '$1 == "SL4.0ch01" { print $2 }' "$reference.fai")
    printf 'SL4.0ch01\t0\t%s\n' "$length" > "$work/regions.bed"
    REGIONS=$work/regions.bed
fi
catalog=${CATALOG:-$work/catalog.parquet}
[[ -e $catalog ]] || "$PVC" repeat-catalog --reference "$reference" --output "$catalog" \
    --log "$work/logs/catalog.log" --force

largest=$(tr ',' '\n' <<< "$sizes" | sort -n | tail -1)
mapfile -t samples < <(cut -f1 "$cohort" | awk '!seen[$0]++' | head -n "$largest")
(( ${#samples[@]} == largest )) || {
    echo "cohort.tsv names ${#samples[@]} samples, fewer than $largest" >&2; exit 1; }

# One psp per sample, in a directory of its own: generate-psps names it by the files' @RG SM, which
# need not be the name cohort.tsv gives. generate-psps is one thread, so THREADS samples run at once.
psp_of() {
    local sample=$1 name path
    while IFS=$'\t' read -r name path; do
        if [[ $name == "$sample" && $path == *.psp ]]; then
            echo "$path"
            return
        fi
    done < "$cohort"
    ls "$work/psps/$sample"/*.psp
}
make_psp() {
    local sample=$1 name path
    local alignments=()
    while IFS=$'\t' read -r name path; do
        [[ $name == "$sample" ]] && alignments+=(--alignment "$path")
    done < "$cohort"
    "$PVC" generate-psps --reference "$reference" --catalog "$catalog" "${alignments[@]}" \
        --regions "$REGIONS" --output-dir "$work/psps/$sample" --force \
        > "$work/logs/psp.$sample.log" 2>&1
    touch "$work/psps/$sample/done"
}
for sample in "${samples[@]}"; do
    [[ -e $work/psps/$sample/done ]] && continue
    awk -F"\t" -v s="$sample" '$1 == s && $2 ~ /\.psp$/ { found = 1 } END { exit !found }' "$cohort" && continue
    while (( $(jobs -rp | wc -l) >= threads )); do
        wait -n
    done
    make_psp "$sample" &
done
while (( $(jobs -rp | wc -l) > 0 )); do
    wait -n
done

for size in ${sizes//,/ }; do
    run=$work/n$size
    mkdir -p "$run"
    psps=()
    for sample in "${samples[@]:0:size}"; do
        psps+=(--psp "$(psp_of "$sample")")
    done
    [[ -e $run/params.toml ]] || "$PVC" estimate-parameters --reference "$reference" \
        --catalog "$catalog" "${psps[@]}" --threads "$threads" --output "$run/params.toml" \
        --log "$run/params.log" --force
    [[ -e $run/calls.vcf.gz ]] || "$PVC" call-from-psps --reference "$reference" \
        --catalog "$catalog" "${psps[@]}" --parameters "$run/params.toml" \
        --refine-genotypes off --threads "$threads" --output "$run/calls.vcf.gz" \
        --log "$run/calls.log" --force
    [[ -e $run/refined.vcf.gz ]] || "$PVC" call-from-psps --reference "$reference" \
        --catalog "$catalog" "${psps[@]}" --parameters "$run/params.toml" \
        --refine-genotypes on --threads "$threads" --output "$run/refined.vcf.gz" \
        --log "$run/refined.log" --force
    {
        echo "== $size accessions"
        grep -h "neighbour-evidence pass" "$run/refined.log" || true
        "${python[@]}" "$here/tomato_fill.py" missing "$run/calls.vcf.gz" "$run/refined.vcf.gz"
        cat "$run/refined.fill_check.tsv"
    } | tee "$run/summary.txt"
done
