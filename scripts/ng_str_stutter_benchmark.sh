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
#   typing_check.tsv           the regions the calling typing makes of the reference (checked
#                              against segments.bed and bundles.bed before anything is mapped)
#   bundle_runs_check.tsv      the caller's reference runs of every typed STR bundle (checked
#                              against the simulator's stutter sites, bundle_runs.tsv)
#   ng_defaults.vcf            call-from-alignments --defaults
#   ng_fitted.vcf              call-from-psps with that fit
#   ng_true_stutter.vcf        call-from-psps with the defaults and the simulator's true stutter
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
#   FIT_FROM    a run directory of the same cohort at another depth whose fit is reused (see ng_fit)
#   REDO        stage names to run again although done, space-separated (e.g. REDO=score)
#   MAIN_BIN    a build of main: adds main's defaults, fitted and true-stutter runs (ng_main_*)
#   BIN_DIR     a directory holding pvc, ng_str_stutter_benchmark, ng_typed_region_dump and
#               ng_str_bundle_runs, used instead of the newest build — so a grid can run on a
#               frozen copy while the tree is rebuilt
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
if [[ -n "${BIN_DIR:-}" ]]; then
    ng=$BIN_DIR/pvc
    simulator=$BIN_DIR/ng_str_stutter_benchmark
    typing_dump=$BIN_DIR/ng_typed_region_dump
    bundle_runs=$BIN_DIR/ng_str_bundle_runs
else
    ng=$(newest pvc)
    simulator=$(newest examples/ng_str_stutter_benchmark)
    typing_dump=$(newest examples/ng_typed_region_dump)
    bundle_runs=$(newest examples/ng_str_bundle_runs)
fi
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
    if [[ -f "$out/$name.done" && " ${REDO:-} " != *" $name "* ]]; then
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
# scored for ng on the generic path alone, and one it invents has no truth. **And every planted STR
# bundle must be typed as one bundle over the laid hull, a separated one with the drawn members, and
# the caller's reference runs there must be the simulator's stutter sites**
# (scripts/ng_str_stutter_check_bundles.py says exactly what matches). Stops the run if not.
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
    if [[ -f "$out/bundles.bed" ]]; then
        if ! NG_REFERENCE_CHECK=skip "$bundle_runs" "$reference" "$out/confident.bed" \
                > "$out/bundle_runs_check.tsv" 2> "$out/bundle_runs_check.log"; then
            echo "ng_str_bundle_runs failed; its log, $out/bundle_runs_check.log:" >&2
            tail -n 20 "$out/bundle_runs_check.log" >&2
            exit 1
        fi
        python3 "$root/scripts/ng_str_stutter_check_bundles.py" "$out"
    fi
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
        --output "$out/ng_defaults.vcf" --log "$out/ng_defaults.pvc.log" \
        > "$out/ng_defaults.log" 2>&1
    # **The numbers the defaults run scored with are in its log**, between two marker lines;
    # ng_str_stutter_true_parameters.py starts from them.
    sed -n '/^# ---- parameters file begins ----$/,/^# ---- parameters file ends ----$/p' \
        "$out/ng_defaults.pvc.log" | sed '1d;$d' > "$out/ng_defaults.parameters.toml"
}

ng_psps() {
    mapfile -t args < <(alignments)
    "$ng" generate-psps --reference "$reference" --catalog "$catalog" "${args[@]}" \
        --regions "$out/confident.bed" --output-dir "$out/psps" --force \
        > "$out/ng_generate_psps.log" 2>&1
}

# **One fit a cohort** (owner, 2026-10-09): the fit is the slow step, and a cohort — one design,
# sample count and seed — is the same samples at every depth, because the simulator draws the
# genotypes and libraries before the depth is used. With FIT_FROM=<run dir of the same cohort>,
# that run's fit is copied after checking the two runs' settings differ in depth alone.
ng_fit() {
    if [[ -n "${FIT_FROM:-}" ]]; then
        local mine theirs
        mine=$(grep '^settings' "$out/simulate.log" | sed 's/depth: [0-9.]*, //; s/threads: [0-9]*, //')
        theirs=$(grep '^settings' "$FIT_FROM/simulate.log" | sed 's/depth: [0-9.]*, //; s/threads: [0-9]*, //')
        if [[ "$mine" != "$theirs" ]]; then
            echo "!! FIT_FROM=$FIT_FROM is another cohort:" >&2
            printf '   %s\n   %s\n' "$theirs" "$mine" >&2
            exit 1
        fi
        cp "$FIT_FROM/ng_fit.parameters.toml" "$out/ng_fit.parameters.toml"
        echo "fit copied from $FIT_FROM" > "$out/ng_estimate_parameters.log"
        return
    fi
    "$ng" estimate-parameters --reference "$reference" --catalog "$catalog" --psp "$out/psps" \
        --output "$out/ng_fit.parameters.toml" --force > "$out/ng_estimate_parameters.log" 2>&1
}

# **Each calling run names its own log**: the default, `pvc_calls.log`, is one name a directory,
# and a second run into `$out` would be refused on the first one's.
ng_fitted() {
    "$ng" call-from-psps --reference "$reference" --catalog "$catalog" --psp "$out/psps" \
        --parameters "$out/ng_fit.parameters.toml" --threads "$threads" \
        --explain-loci "$out/segments.bed" "$out/ng_fitted.explain.tsv" \
        --output "$out/ng_fitted.vcf" --log "$out/ng_fitted.pvc.log" \
        > "$out/ng_call_from_psps.log" 2>&1
}

# **ng with the true stutter** (owner, 2026-10-09): the defaults run's parameters file with its
# stutter tables replaced by the simulator's (scripts/ng_str_stutter_true_parameters.py) — no fit.
# Against the fitted run it says whether the fit can be skipped; against --defaults, what a right
# stutter model buys.
ng_true_stutter() {
    python3 "$root/scripts/ng_str_stutter_true_parameters.py" "$out" > "$out/ng_true_stutter_parameters.log"
    "$ng" call-from-psps --reference "$reference" --catalog "$catalog" --psp "$out/psps" \
        --parameters "$out/true_stutter.parameters.toml" --threads "$threads" \
        --explain-loci "$out/segments.bed" "$out/ng_true_stutter.explain.tsv" \
        --output "$out/ng_true_stutter.vcf" --log "$out/ng_true_stutter.pvc.log" \
        > "$out/ng_true_stutter.log" 2>&1
}

# **The same calls made by main** (owner, 2026-10-09), when MAIN_BIN names a build of main: the
# comparison that says what calling each STR segment both ways changed. Main fits with the same
# code and reads this branch's parameters files unchanged, so its fitted and true-stutter runs use
# this run's two files. It cannot read this branch's catalog (a catalog records the build that
# wrote it) or its psps (format 2.0), so it writes its own. Main lacks the compound-motif fix
# (db7388ec): it types the AAAAAC, AAAAAG and AAAAAT hexanucleotides as ordinary sequence, which
# the scorer's `motif` column lets the report leave out.
main_psps() {
    mapfile -t args < <(alignments)
    "$MAIN_BIN" repeat-catalog --reference "$reference" --output "$out/main_catalog.parquet" --force \
        > "$out/main_catalog.log" 2>&1
    "$MAIN_BIN" generate-psps --reference "$reference" --catalog "$out/main_catalog.parquet" \
        "${args[@]}" --regions "$out/confident.bed" --output-dir "$out/psps_main" --force \
        > "$out/main_generate_psps.log" 2>&1
}

ng_main_defaults() {
    mapfile -t args < <(alignments)
    "$MAIN_BIN" call-from-alignments --reference "$reference" --catalog "$out/main_catalog.parquet" \
        "${args[@]}" --regions "$out/confident.bed" --defaults --threads "$threads" \
        --output "$out/ng_main_defaults.vcf" > "$out/ng_main_defaults.log" 2>&1
}

# ng_main_called_with <name> <parameters file>
ng_main_called_with() {
    "$MAIN_BIN" call-from-psps --reference "$reference" --catalog "$out/main_catalog.parquet" \
        --psp "$out/psps_main" --parameters "$2" --threads "$threads" \
        --output "$out/$1.vcf" > "$out/$1.log" 2>&1
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
# **Not the STR bundles** (bundles.bed): a HipSTR region carries one period and HipSTR fits one
# stutter model to it, so a bundle given as its hull would be called with every member's stutter
# priced at one member's period, and its calls would measure that mismatch rather than HipSTR. It is
# scored on segments only, as before; GATK and ng are scored on the bundles.
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
stage ng_true_stutter ng_true_stutter
if [[ -n "${MAIN_BIN:-}" ]]; then
    stage main_psps main_psps
    stage ng_main_defaults ng_main_defaults
    stage ng_main_true_stutter ng_main_called_with ng_main_true_stutter "$out/true_stutter.parameters.toml"
fi
stage ng_fit ng_fit
stage ng_fitted ng_fitted
if [[ -n "${MAIN_BIN:-}" ]]; then
    stage ng_main_fitted ng_main_called_with ng_main_fitted "$out/ng_fit.parameters.toml"
fi
stage gatk run_gatk
stage hipstr run_hipstr hipstr
stage hipstr_fixed_stutter run_hipstr hipstr_fixed_stutter --def-stutter-model --min-reads 1
stage score python3 "$root/scripts/ng_str_stutter_score.py" "$out"
echo "done: $out"
