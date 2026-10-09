#!/usr/bin/env bash
# **One configuration of the stuttering benchmark, end to end**
# (doc/devel/ng/spec/str_stutter_benchmark.md §9): simulate a cohort, map it with minimap2, call it
# with ng both ways, GATK and HipSTR, and leave every output where the scorer finds it.
#
#   scripts/ng_str_stutter_benchmark.sh <out-dir> [simulator key=value ...]
#
# Run it inside the dev container, which has minimap2, samtools, GATK (/opt/gatk/gatk) and HipSTR:
#
#   ./scripts/dev.sh scripts/ng_str_stutter_benchmark.sh tmp/stutter_bench/pcr_50_30 \
#       design=pcr samples=50 depth=30
#
# The key=value settings go to `ng_str_stutter_benchmark` unchanged (its module documentation
# lists them). Every stage writes `<stage>.done` when it finishes, and a rerun skips the stages
# already done, so an interrupted run resumes; delete a marker to redo its stage and the ones
# after it. Wall-clock seconds per stage go to `timings.tsv`. **The callers' times are not a
# performance comparison**: GATK's per-sample step runs several samples at once (GATK_JOBS) to
# keep a 50-sample run to hours, which is not how the other callers are run.
#
# What it writes into <out-dir>, beside the simulator's own files:
#   bams/<sample>.bam          minimap2 -ax sr, sorted and indexed; one read group a sample
#   typing_check.tsv           the regions the calling typing makes of the reference
#   ng_defaults.vcf            call-from-alignments --defaults
#   ng_fitted.vcf              call-from-psps with that fit
#   ng_fit.parameters.toml     the fit (estimate-parameters over the psps)
#   ng_*.explain.tsv           --explain-loci over every STR segment
#   gatk.vcf.gz                HaplotypeCaller -ERC GVCF, CombineGVCFs, GenotypeGVCFs
#   hipstr.vcf.gz              HipSTR over the STR segments, all samples jointly, its defaults
#   hipstr_fixed_stutter.vcf.gz  the same with one fixed stutter model and no read floor
#   score_*.tsv                the scorer's tables (scripts/ng_str_stutter_score.py)
#
# Settings (environment):
#   THREADS     threads for minimap2, samtools and ng (default: all)
#   GATK_JOBS   HaplotypeCaller processes at once (default: THREADS / 2, at least 1)
#   GATK_HEAP   Java heap per GATK process (default 4g)
#   KEEP_FASTQ  1 keeps each sample's FASTQ after mapping (default: deleted, ~230 MB a sample at 60×)
set -euo pipefail

out=${1:?usage: $0 <out-dir> [simulator key=value ...]}
shift
root=$(cd "$(dirname "$0")/.." && pwd)
threads=${THREADS:-$(nproc)}
gatk_jobs=${GATK_JOBS:-$(( threads / 2 > 0 ? threads / 2 : 1 ))}
gatk_heap=${GATK_HEAP:-4g}
gatk=${GATK_BIN:-/opt/gatk/gatk}
hipstr=${HIPSTR_BIN:-HipSTR}

# The newer of the two release trees: the container builds into target-container/, a machine
# without a container runtime into target/ (CLAUDE.md, "Container vs. host").
newest() {
    local best="" candidate
    for candidate in "$root/target-container/release/$1" "$root/target/release/$1"; do
        if [[ -x "$candidate" && ( -z "$best" || "$candidate" -nt "$best" ) ]]; then
            best=$candidate
        fi
    done
    [[ -n "$best" ]] || { echo "no built $1; build it with cargo build --release" >&2; exit 1; }
    echo "$best"
}
ng=$(newest pop_var_caller)
simulator=$(newest examples/ng_str_stutter_benchmark)
typing_dump=$(newest examples/ng_typed_region_dump)
for tool in minimap2 samtools bcftools "$gatk" "$hipstr"; do
    command -v "$tool" > /dev/null || { echo "$tool not found; run this in the dev container" >&2; exit 1; }
done

mkdir -p "$out"
out=$(cd "$out" && pwd)
reference=$out/reference.fa
catalog=$out/reference.fa.repeats.parquet
[[ -f "$out/timings.tsv" ]] || printf 'stage\tseconds\n' > "$out/timings.tsv"

# stage <name> <command...>: run once, record its time, mark it done.
stage() {
    local name=$1
    shift
    if [[ -f "$out/$name.done" ]]; then
        echo "[skip] $name"
        return 0
    fi
    echo "[run]  $name"
    local started=$SECONDS
    "$@"
    printf '%s\t%d\n' "$name" $(( SECONDS - started )) >> "$out/timings.tsv"
    touch "$out/$name.done"
}

simulate() {
    "$simulator" "$out" "$@" > "$out/simulate.log"
    samtools faidx "$reference"
    "$ng" repeat-catalog --reference "$reference" --output "$catalog" > "$out/catalog.log" 2>&1
}

# **The typed STR segments must be the planted ones**, exactly: a segment the typing misses is
# scored for ng on the generic path alone, and one it invents has no truth. Stops the run if not.
check_typing() {
    NG_REFERENCE_CHECK=skip "$typing_dump" "$reference" "$out/confident.bed" calling \
        > "$out/typing_check.tsv"
    local typed planted
    typed=$(awk -F'\t' '$4 == "ssr_locus" { print $1 "\t" ($2 - 1) "\t" $3 "\t" $6 }' \
        "$out/typing_check.tsv" | sort)
    planted=$(sort "$out/segments.bed")
    if [[ "$typed" != "$planted" ]]; then
        echo "!! the typed STR segments are not the planted ones:" >&2
        diff <(echo "$typed") <(echo "$planted") | head -10 >&2
        exit 1
    fi
    echo "       typing: $(echo "$planted" | wc -l) STR segments, all as planted"
}

samples() {
    tail -n +2 "$out/libraries.tsv" | cut -f1
}

map_reads() {
    mkdir -p "$out/bams"
    local sample
    for sample in $(samples); do
        minimap2 -ax sr -t "$threads" -R "@RG\tID:$sample\tSM:$sample\tLB:$sample\tPL:ILLUMINA" \
            "$reference" "$out/${sample}_R1.fastq.gz" "$out/${sample}_R2.fastq.gz" 2>> "$out/minimap2.log" \
            | samtools sort -@ "$threads" -o "$out/bams/$sample.bam" - 2>> "$out/samtools.log"
        samtools index "$out/bams/$sample.bam"
        [[ "${KEEP_FASTQ:-0}" == 1 ]] || rm -f "$out/${sample}_R1.fastq.gz" "$out/${sample}_R2.fastq.gz"
    done
}

alignments() {
    local sample
    for sample in $(samples); do
        printf -- '--alignment\n%s\n' "$out/bams/$sample.bam"
    done
}

ng_defaults() {
    mapfile -t args < <(alignments)
    "$ng" call-from-alignments --reference "$reference" --catalog "$catalog" "${args[@]}" \
        --regions "$out/confident.bed" --defaults --threads "$threads" \
        --explain-loci "$out/segments.bed" "$out/ng_defaults.explain.tsv" \
        --output "$out/ng_defaults.vcf" > "$out/ng_defaults.log" 2>&1
}

ng_psps() {
    mapfile -t args < <(alignments)
    "$ng" generate-psps --reference "$reference" --catalog "$catalog" "${args[@]}" \
        --regions "$out/confident.bed" --output-dir "$out/psps" --force \
        > "$out/ng_generate_psps.log" 2>&1
}

ng_fit() {
    "$ng" estimate-parameters --reference "$reference" --catalog "$catalog" --psp "$out/psps" \
        --output "$out/ng_fit.parameters.toml" --force > "$out/ng_estimate_parameters.log" 2>&1
}

# The fit has its own name: ng writes `<output>.parameters.toml` beside every VCF, and refuses
# to overwrite the file it was handed.
ng_fitted() {
    "$ng" call-from-psps --reference "$reference" --catalog "$catalog" --psp "$out/psps" \
        --parameters "$out/ng_fit.parameters.toml" --threads "$threads" \
        --explain-loci "$out/segments.bed" "$out/ng_fitted.explain.tsv" \
        --output "$out/ng_fitted.vcf" > "$out/ng_call_from_psps.log" 2>&1
}

gatk_one() {
    local sample=$1
    "$gatk" --java-options "-Xmx$gatk_heap" HaplotypeCaller -R "$reference" \
        -I "$out/bams/$sample.bam" -O "$out/gatk/$sample.g.vcf.gz" -ERC GVCF \
        > "$out/gatk/$sample.log" 2>&1
}

run_gatk() {
    mkdir -p "$out/gatk"
    [[ -f "$out/reference.dict" ]] || \
        "$gatk" CreateSequenceDictionary -R "$reference" > "$out/gatk/dict.log" 2>&1
    export -f gatk_one
    export gatk gatk_heap reference out
    samples | xargs -P "$gatk_jobs" -I{} bash -c 'gatk_one "$1"' _ {}
    local variants=()
    local sample
    for sample in $(samples); do
        variants+=(--variant "$out/gatk/$sample.g.vcf.gz")
    done
    "$gatk" --java-options "-Xmx$gatk_heap" CombineGVCFs -R "$reference" "${variants[@]}" \
        -O "$out/gatk/combined.g.vcf.gz" > "$out/gatk/combine.log" 2>&1
    "$gatk" --java-options "-Xmx$gatk_heap" GenotypeGVCFs -R "$reference" \
        -V "$out/gatk/combined.g.vcf.gz" -O "$out/gatk.vcf.gz" > "$out/gatk/genotype.log" 2>&1
}

# HipSTR's region file: chrom, 1-based start, end, period, copies, name — the STR segments.
#
# **Twice, because its defaults call nothing at low total depth**: HipSTR learns each locus's
# stutter from the cohort's reads and skips a locus with fewer than 100 reads over all samples
# (`--min-reads`), which is every locus of one sample below about 100×. The second run is its own
# answer for that case: a fixed stutter model for every locus (`--def-stutter-model`: in-frame
# slips 5% up and 5% down, geometric step 0.95 as its VCF reports — its help says 0.9) and no
# read floor.
hipstr_regions() {
    awk -F'\t' -v OFS='\t' '{ print $1, $2 + 1, $3, $4, ($3 - $2) / $4, sprintf("seg%05d", NR - 1) }' \
        "$out/segments.bed" > "$out/hipstr_regions.bed"
}

run_hipstr() {
    local name=$1
    shift
    [[ -s "$out/hipstr_regions.bed" ]] || hipstr_regions
    local bams names
    bams=$(samples | sed "s#^#$out/bams/#; s#\$#.bam#" | paste -sd, -)
    names=$(samples | paste -sd, -)
    "$hipstr" --bams "$bams" --bam-samps "$names" --bam-libs "$names" --fasta "$reference" \
        --regions "$out/hipstr_regions.bed" --str-vcf "$out/$name.vcf.gz" \
        --log "$out/$name.log" "$@" > /dev/null 2>&1
}

stage simulate simulate "$@"
stage check_typing check_typing
stage map_reads map_reads
stage ng_defaults ng_defaults
stage ng_psps ng_psps
stage ng_fit ng_fit
stage ng_fitted ng_fitted
stage gatk run_gatk
stage hipstr run_hipstr hipstr
stage hipstr_fixed_stutter run_hipstr hipstr_fixed_stutter --def-stutter-model --min-reads 1
stage score python3 "$root/scripts/ng_str_stutter_score.py" "$out"
echo "done: $out"
