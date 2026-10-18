#!/bin/sh
# GATK 4.6 HaplotypeCaller over the Tandem-Repeat set's satellite ground, padded 150 bp, at 30x and
# 10x (plan call_satellites.md, step 6). Run inside the dev container.
set -eu
W=/Users/jose/devel/pop_var_caller-call-satellites/tmp/giab
B=/Users/jose/devel/pop_var_caller/benchmarks
REF=/Users/jose/genomes/h_sapiens/gca_grch38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
mkdir -p $W/gatk
for depth in 30x 10x; do
  /opt/gatk/gatk --java-options "-Xmx6g" HaplotypeCaller -R $REF \
    -I $B/ssr_hg002/bam/$depth/HG002_TR_v1.0.1_Tier_$depth.bam \
    -L $W/sat_tr.sorted.bed --interval-padding 150 -stand-call-conf 0 \
    -O $W/gatk/tr${depth}.vcf.gz > $W/gatk/tr${depth}.log 2>&1
  echo "$depth: $(zcat $W/gatk/tr${depth}.vcf.gz | grep -vc '^#') records"
done
