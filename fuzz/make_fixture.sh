#!/bin/bash
# **Builds the fuzz targets' fixture**: a run small enough to repeat thousands of times a second,
# made of real reads so the fuzzer starts from inputs that reach deep into the caller.
#
# Two 10 kb stretches of GRCh38 (chosen for their repeats and variants in the GIAB benchmark
# regions) become a two-contig reference, ctgA and ctgB. The 5x GIAB reads of HG002 and HG003 whose
# pair lies wholly inside a stretch are moved onto it. Then pvc builds everything else a run reads:
# the repeat catalog, the psps, and a parameters file fitted from them.
#
# Usage, inside the dev container:
#
#   fuzz/make_fixture.sh <GRCh38.fa> <giab_trio-dir> <pvc-binary>
#
# Writes fuzz/fixture/ and the targets' seed corpora under fuzz/corpus/.
set -eu

reference=$1; giab=$2; pvc=$3
here=$(cd "$(dirname "$0")" && pwd)
fixture=$here/fixture
rm -rf "$fixture"; mkdir -p "$fixture"

windows=("chr1 110600000 110610000 ctgA" "chr15 78836000 78846000 ctgB")

: > "$fixture/reference.fa"
for window in "${windows[@]}"; do
  set -- $window
  samtools faidx "$reference" "$1:$2-$3" | sed "1s/.*/>$4/" >> "$fixture/reference.fa"
done
samtools faidx "$fixture/reference.fa"

for sample in HG002 HG003; do
  bam=$giab/bam/5x/$sample/$sample.5x.seed42.bam
  {
    samtools view -H "$bam" | grep -v '^@SQ\|^@PG'
    awk -F'\t' '{ printf "@SQ\tSN:%s\tLN:%s\n", $1, $2 }' "$fixture/reference.fa.fai"
    for window in "${windows[@]}"; do
      set -- $window
      # A pair is kept when both reads start inside the window and the read ends inside it
      # (150-base reads, so a start more than 160 bases from the end is safe); positions move
      # by the window's start.
      samtools view "$bam" "$1:$2-$3" | awk -F'\t' -v OFS='\t' -v s="$2" -v e="$3" -v c="$4" '
        $7 == "=" && $4 >= s && $4 <= e - 160 && $8 >= s && $8 <= e - 160 {
          $3 = c; $4 = $4 - s + 1; $8 = $8 - s + 1; print }'
    done
  } | samtools sort -o "$fixture/$sample.bam" -
  samtools index "$fixture/$sample.bam"
done

printf 'ctgA\t0\t10001\nctgB\t0\t10001\n' > "$fixture/regions.bed"
"$pvc" repeat-catalog --reference "$fixture/reference.fa" --output "$fixture/catalog.parquet" --force
"$pvc" generate-psps --reference "$fixture/reference.fa" --catalog "$fixture/catalog.parquet" \
  --alignment "$fixture/HG002.bam" --alignment "$fixture/HG003.bam" \
  --output-dir "$fixture/psps" --force
"$pvc" estimate-parameters --reference "$fixture/reference.fa" --catalog "$fixture/catalog.parquet" \
  --psp "$fixture"/psps/HG002*.psp --psp "$fixture"/psps/HG003*.psp \
  --output "$fixture/parameters.toml" --force
rm -f "$fixture"/*.log "$fixture"/psps/*.log

# Seeds. The BAM target takes the decompressed stream (it compresses each input itself, so the
# fuzzer mutates records rather than deflate blocks); a 30-record cut is a second, small seed.
corpus=$here/corpus
rm -rf "$corpus"; mkdir -p "$corpus"/{bam,cram,fasta,bed,parameters,psp}
bgzip -dc "$fixture/HG002.bam" > "$corpus/bam/HG002"
{ samtools view -H "$fixture/HG002.bam"; samtools view "$fixture/HG002.bam" | head -30; } |
  samtools view -b -o - - | bgzip -dc > "$corpus/bam/HG002_30_records"
samtools view -C -T "$fixture/reference.fa" -o "$corpus/cram/HG002" "$fixture/HG002.bam"
# The FASTA target's input is the FASTA and its index joined by one line holding only `@@FAI`.
{ cat "$fixture/reference.fa"; echo '@@FAI'; cat "$fixture/reference.fa.fai"; } > "$corpus/fasta/fixture"
cp "$fixture/regions.bed" "$corpus/bed/fixture"
printf 'ctgA\t100\t200\tname\nctgB\t5000\t9000\n# comment\n' > "$corpus/bed/four_columns"
cp "$fixture/parameters.toml" "$corpus/parameters/fitted"
cp "$fixture"/psps/HG002*.psp "$corpus/psp/HG002"
echo "fixture in $fixture, seeds in $corpus"
