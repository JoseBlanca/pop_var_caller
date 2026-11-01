#!/bin/bash
# **Relations between runs on transformed reads**: the companion of run_invariance_oracle.sh. That
# script reruns the same reads in ways that must not change a byte; this one changes the reads in
# ways whose effect on the calls is known, and checks the calls moved only as they should.
#
#   doubled  every read fed twice, each copy under a new name. No QUAL may fall and no genotype
#            with GQ >= 20 may change. **Only at depth.** A candidate allele needs 2 reads, so at
#            a few reads a position every single misread becomes a candidate once doubled; on
#            the 3x tomato extract that added 109,698 records and set aside 15,635 confident
#            calls (2026-11-01), while GIAB at 30x held. DOUBLED=measure reports it without
#            failing, DOUBLED=skip leaves it out; the default is check.
#   empty    one more sample, with a read group and no reads. Every other sample's calls, and every
#            site field, must stay the same.
#   strand   every read's strand flag flipped, and its mate's. Genotypes and SPPEN, the strand and
#            read-position penalty, must stay the same. Run on copies of both sides without the
#            pairs whose fragment is shorter than a read — see below.
#   alone    each sample called by itself, against the same sample in the cohort. A measurement:
#            the differences are counted and characterised, never a failure.
#
# **What flipping the strand means in a BAM.** The record stores the read's bases and qualities
# already turned to the reference's forward strand, so "the same fragment, read off the other
# strand" is the record with its 0x10 bit flipped and nothing else: its bases, qualities, CIGAR and
# position are the reference-strand view either way. The mate's 0x20 bit is flipped with it, so
# the pair stays consistent. One thing does change: a pair whose fragment is shorter than a read
# has its bases past the adaptor dropped by a test that depends on which read is forward. On the
# tomato extract those pairs alone changed 238 genotypes, so both sides of the strand relation drop
# them first (`samtools view -e`), and how many were dropped is printed beside the result.
#
# Every run is call-from-alignments with --defaults, so no fit can differ between the two sides,
# and with --min-site-quality 0 --hidden-dup-tag, so a call that got worse shows up as a worse call
# rather than as an absence. The comparisons are in read_transform_compare.py.
#
# Every relation runs and reports; the script exits 1 at the end if any pass/fail one broke.
#
# Usage:
#
#   scripts/read_transform_oracle.sh <reference.fa> <catalog.parquet> <regions.bed> \
#       <output-dir> <alignment.bam>...
#
# Needs samtools. PVC_BIN picks the binary; without it, the newest of
# target{,-container}/{soak,release}/pvc. THREADS (default 6) sets --threads.
set -eu

if [ "$#" -lt 5 ]; then
  sed -n '2,${/^#/!q;p;}' "$0"
  exit 2
fi

reference=$1; catalog=$2; regions=$3; out=$4
shift 4
alignments=("$@")
threads=${THREADS:-6}
doubled_mode=${DOUBLED:-check}
case $doubled_mode in check|measure|skip) ;; *) echo "DOUBLED is check, measure or skip" >&2; exit 2;; esac

root=$(cd "$(dirname "$0")/.." && pwd)
bin=${PVC_BIN:-}
if [ -z "$bin" ]; then
  for candidate in "$root"/target{,-container}/{soak,release}/pvc; do
    if [ -x "$candidate" ] && { [ -z "$bin" ] || [ "$candidate" -nt "$bin" ]; }; then
      bin=$candidate
    fi
  done
fi
if [ -z "$bin" ]; then
  echo "no build of pvc; build one first (cargo build --profile soak --bin pvc)" >&2
  exit 1
fi
echo "binary: $bin"
mkdir -p "$out/bam"
compare=(python3 "$root/scripts/read_transform_compare.py")

failures=0
# check <label> <relation> <base.vcf> <other.vcf> [sample]
check() {
  local label=$1; shift
  echo "--- $label"
  if "${compare[@]}" "$@"; then
    echo "HELD       $label"
  else
    echo "BROKEN     $label"
    failures=$((failures + 1))
  fi
}

call() {
  local name=$1; shift
  "$bin" call-from-alignments --reference "$reference" --catalog "$catalog" \
    $(for f in "$@"; do printf -- '--alignment\n%s\n' "$f"; done) \
    --defaults --regions "$regions" --min-site-quality 0 --hidden-dup-tag \
    --threads "$threads" --output "$out/$name.vcf" --log "$out/$name.log" --force \
    > "$out/$name.stdout" 2>&1
}

# The sample a BAM's first read group names.
sample_of() {
  samtools view -H "$1" | awk -F'\t' '/^@RG/ { for (i = 2; i <= NF; i++)
    if ($i ~ /^SM:/) { print substr($i, 4); exit } }'
}

# transform <awk-program> <input.bam> <output.bam>: the program sees every record line.
transform() {
  samtools view -h "$2" |
    awk -F'\t' -v OFS='\t' "/^@/ { print; next } $1" |
    samtools view -b -o "$3" -
  samtools index "$3"
}

# full_length <input.bam> <output.bam>: the pairs whose fragment is at least a read long.
full_length() {
  samtools view -b -e '(tlen >= qlen) || (tlen <= -qlen)' -o "$2" "$1"
  samtools index "$2"
}

doubled=(); full_length_base=(); swapped=()
for bam in "${alignments[@]}"; do
  name=$(basename "$bam" .bam)
  if [ "$doubled_mode" != skip ]; then
    transform '{ print; $1 = $1 "_copy"; print }' "$bam" "$out/bam/$name.doubled.bam"
    doubled+=("$out/bam/$name.doubled.bam")
  fi
  full_length "$bam" "$out/bam/$name.full_length.bam"
  # mawk has no bit operations: flip 0x10 and 0x20 by arithmetic.
  transform '{ f = $2
    f += (int(f / 16) % 2) ? -16 : 16
    if (int(f / 1) % 2) f += (int(f / 32) % 2) ? -32 : 32
    $2 = f; print }' "$out/bam/$name.full_length.bam" "$out/bam/$name.swapped.bam"
  full_length_base+=("$out/bam/$name.full_length.bam"); swapped+=("$out/bam/$name.swapped.bam")
done

# The empty sample: the first BAM's header with its read groups renamed, and no records.
samtools view -H "${alignments[0]}" |
  sed -E '/^@RG/ { s/\tID:[^\t]*/\tID:EMPTY/; s/\tSM:[^\t]*/\tSM:EMPTY/; s/\tLB:[^\t]*/\tLB:EMPTY/; }' |
  samtools view -b -o "$out/bam/empty.bam" -
samtools index "$out/bam/empty.bam"

echo "=== calling ==="
call base "${alignments[@]}"
[ "$doubled_mode" = skip ] || call doubled "${doubled[@]}"
call full_length "${full_length_base[@]}"
call swapped "${swapped[@]}"
call with_empty "${alignments[@]}" "$out/bam/empty.bam"
for bam in "${alignments[@]}"; do
  call "alone_$(sample_of "$bam")" "$bam"
done

echo "=== relations ==="
case $doubled_mode in
  check) check "every read twice" doubled "$out/base.vcf" "$out/doubled.vcf" ;;
  measure)
    echo "--- every read twice (a measurement: DOUBLED=measure)"
    "${compare[@]}" doubled "$out/base.vcf" "$out/doubled.vcf" || true ;;
esac
check "strands swapped, short fragments dropped on both sides" strand "$out/full_length.vcf" \
  "$out/swapped.vcf"
short=0
for bam in "${alignments[@]}"; do
  short=$((short + $(samtools view "$bam" |
    awk '{ t = $9 < 0 ? -$9 : $9; if (t < length($10)) n++ } END { print n + 0 }')))
done
echo "    ($short records dropped from both sides: a fragment shorter than the read, whose" \
  "adaptor test reads the strand, or no fragment length at all)"
check "an empty sample added" empty "$out/base.vcf" "$out/with_empty.vcf" EMPTY
for bam in "${alignments[@]}"; do
  sample=$(sample_of "$bam")
  echo "--- $sample alone against $sample in the cohort (a measurement)"
  "${compare[@]}" alone "$out/base.vcf" "$out/alone_$sample.vcf" "$sample"
done

echo
if [ "$failures" -eq 0 ]; then
  echo "every relation held"
else
  echo "$failures relation(s) broken"
  exit 1
fi
