#!/bin/sh
# pvc on GIAB HG002 with and without --call-satellites (plan call_satellites.md, step 6).
# 300x over the 1,000 small-variant regions; 30x and 10x over the Tandem-Repeat Tier regions.
set -eu
W=/Users/jose/devel/pop_var_caller-call-satellites
P=$W/target/release/pvc
B=/Users/jose/devel/pop_var_caller/benchmarks
REF=$B/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
OUT=$W/tmp/giab/pvc
mkdir -p "$OUT"

run() { # name alignment regions extra...
  name=$1; alignment=$2; regions=$3; shift 3
  /usr/bin/time -l "$P" call-from-alignments --reference "$REF" --alignment "$alignment" \
    --regions "$regions" --defaults --threads 6 --output "$OUT/$name.vcf" \
    --log "$OUT/$name.log" --force "$@" > "$OUT/$name.out" 2>&1
  echo "$name: $(grep -vc '^#' "$OUT/$name.vcf") records, $(grep 'maximum resident' "$OUT/$name.out")"
}

for setting in default call_satellites; do
  extra=""
  [ "$setting" = call_satellites ] && extra="--call-satellites"
  run "300x_$setting" "$B/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram" \
    "$B/human_genome_bottle/HG002_bench_azar_merged_1000.bed" $extra
  for depth in 30x 10x; do
    run "tr${depth}_$setting" "$B/ssr_hg002/bam/$depth/HG002_TR_v1.0.1_Tier_$depth.bam" \
      "$W/tmp/giab/tr_regions.bed" $extra
  done
done
