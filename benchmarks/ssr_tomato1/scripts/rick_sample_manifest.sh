#!/usr/bin/env bash
# Emit a per-sample manifest for a directory of CRAMs/BAMs — run this ON RICK, before slicing.
#
# We need to choose which samples to bring across, and two properties decide it:
#
#   * read length — the complete-vs-partial frontier IS a read-length effect, so mixing a 100 bp
#     project with a 150 bp one makes stutter appear to start at different tract lengths for purely
#     geometric reasons. Either restrict to one length or stratify by it; you cannot ignore it.
#   * library / experiment — the grouping the whole question turns on. One project can mix
#     preparations, so the unit is the experiment (the @RG LB/ID), not the project directory.
#
# Read length is sampled from the first records rather than measured over the file, so this stays
# fast on full-size CRAMs. Sampling is fine for a modal length; it would not be fine for a mean.
#
# Usage:
#   ./rick_sample_manifest.sh IN_DIR [REF] > manifest.tsv
#
#   REF is needed only if the CRAMs cannot find their reference via the header URI.
#
# Example:
#   for d in /media/tomato25_bams/crams/*/; do
#       ./rick_sample_manifest.sh "$d" /home/joxi/refs/S_lycopersicum_chromosomes.4.00.fa
#   done > ~/tmp/tomato_manifest.tsv

set -euo pipefail

IN_DIR=${1:?"usage: $0 IN_DIR [REF]"}
REF=${2:-}
SAMPLE_READS=${SAMPLE_READS:-2000}

[[ -d "$IN_DIR" ]] || { echo "not a directory: $IN_DIR" >&2; exit 1; }

ref_args=()
[[ -n "$REF" ]] && ref_args=(-T "$REF")

shopt -s nullglob
files=("$IN_DIR"/*.cram "$IN_DIR"/*.bam)
if (( ${#files[@]} == 0 )); then
    echo "no *.cram or *.bam in $IN_DIR" >&2
    exit 1
fi

# One header line, emitted only when stdout is not being appended to an existing manifest.
if [[ -z "${MANIFEST_NO_HEADER:-}" ]]; then
    printf 'dir\tfile\trun\tsample\tlibrary\tplatform\tread_len_mode\tread_len_spread\tn_sampled\n'
fi

for f in "${files[@]}"; do
    base=$(basename "$f")

    # First @RG only: these are single-read-group files, and a second would mean the sample was
    # merged across libraries — which the `library` column would then misreport, so flag it.
    rg=$(samtools view -H "${ref_args[@]}" "$f" 2>/dev/null | grep -m1 '^@RG' || true)
    n_rg=$(samtools view -H "${ref_args[@]}" "$f" 2>/dev/null | grep -c '^@RG' || true)

    field() { printf '%s' "$rg" | tr '\t' '\n' | grep -m1 "^$1:" | cut -d: -f2- || true; }
    run=$(field ID); sm=$(field SM); lb=$(field LB); pl=$(field PL)
    [[ "$n_rg" -gt 1 ]] && lb="${lb}+MERGED(${n_rg}RG)"

    # Modal read length over a sample of records, plus how many distinct lengths appear: a wide
    # spread means trimmed reads, which changes spanning behaviour as surely as a shorter read does.
    lens=$(samtools view "${ref_args[@]}" "$f" 2>/dev/null | head -n "$SAMPLE_READS" \
           | awk '{print length($10)}' | sort -n | uniq -c | sort -rn || true)
    mode=$(printf '%s' "$lens" | head -1 | awk '{print $2}')
    distinct=$(printf '%s' "$lens" | grep -c . || true)
    n=$(printf '%s' "$lens" | awk '{s+=$1} END {print s+0}')

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$(basename "$IN_DIR")" "$base" "${run:-?}" "${sm:-?}" "${lb:-?}" "${pl:-?}" \
        "${mode:-?}" "${distinct:-?}" "${n:-0}"
done
