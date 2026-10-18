#!/bin/sh
# Score pvc (default, --call-satellites) and GATK against GIAB HG002 truth over satellite ground,
# with the project's exact matcher (benchmarks/lib/compare_to_truth.sh's method) and with RTG
# vcfeval in genotype and allele (--squash-ploidy) modes. Plan call_satellites.md, step 6.
# Run inside the dev container (Java, bcftools). Writes $W/scores.tsv.
set -eu
W=/Users/jose/devel/pop_var_caller-call-satellites/tmp/giab
B=/Users/jose/devel/pop_var_caller/benchmarks
REF=$B/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
RTG="java -Xmx6g -jar /Users/jose/devel/pop_var_caller-call-satellites/tmp/tools/rtg-tools-3.13/RTG.jar"
S=$W/score
rm -rf $S && mkdir -p $S/truth
TSV=$W/scores.tsv
printf 'dataset\tcaller\tmeasure\tclass\ttp\tfp\tfn\n' > $TSV

# A VCF bgzipped and indexed, under $S, named $2.
indexed() {
  case "$1" in
    *.gz) bcftools view -Oz -o "$S/$2.vcf.gz" "$1" ;;
    *) bcftools view -Oz -o "$S/$2.vcf.gz" "$1" ;;
  esac
  bcftools index -t -f "$S/$2.vcf.gz"
  echo "$S/$2.vcf.gz"
}

# Records of a VCF split into SNPs and the rest, as "snps others".
classes() {
  bcftools view -H "$1" 2>/dev/null | awk -F'\t' '{
    n = split($5, alts, ","); snp = (length($4) == 1);
    for (i = 1; i <= n; i++) if (length(alts[i]) != 1) snp = 0;
    if (snp) s++; else o++ } END { print s + 0, o + 0 }'
}

# The exact matcher's normalisation: BED-restrict, left-align, split, keep one class.
normalize() { # src class out bed [extra view args]
  src=$1; cls=$2; out=$3; bed=$4; shift 4
  bcftools view "$@" -T "$bed" "$src" -Ou \
    | bcftools norm -f "$REF" -m -any -Ou 2>/dev/null \
    | bcftools view -v "$cls" -Oz -o "$out"
  bcftools index -t -f "$out"
}

score() { # dataset truth bed caller calls
  dataset=$1; truth=$2; bed=$3; caller=$4; calls=$5
  tag=${dataset}_$caller

  set -- $(classes "$(bcftools view -T "$bed" -Oz -o $S/$tag.on_sat.vcf.gz "$calls" && echo $S/$tag.on_sat.vcf.gz)")
  printf '%s\t%s\tcalls\tsnp\t%s\t.\t.\n%s\t%s\tcalls\tother\t%s\t.\t.\n' \
    "$dataset" "$caller" "$1" "$dataset" "$caller" "$2" >> $TSV

  for cls in snps indels; do
    normalize "$truth" $cls $S/$tag.truth.$cls.vcf.gz "$bed" -f PASS,.
    normalize "$calls" $cls $S/$tag.calls.$cls.vcf.gz "$bed"
    bcftools isec -p $S/$tag.isec.$cls $S/$tag.truth.$cls.vcf.gz $S/$tag.calls.$cls.vcf.gz \
      > /dev/null 2>&1
    fn=$(grep -vc '^#' $S/$tag.isec.$cls/0000.vcf || true)
    fp=$(grep -vc '^#' $S/$tag.isec.$cls/0001.vcf || true)
    tp=$(grep -vc '^#' $S/$tag.isec.$cls/0003.vcf || true)
    printf '%s\t%s\texact\t%s\t%s\t%s\t%s\n' "$dataset" "$caller" $cls $tp $fp $fn >> $TSV
  done

  for mode in genotype allele; do
    extra=""
    [ $mode = allele ] && extra="--squash-ploidy"
    out=$S/$tag.vcfeval.$mode
    $RTG vcfeval -b "$truth" -c "$calls" -t $W/grch38.sdf --bed-regions "$bed" \
      --sample HG002 -o $out $extra > $out.log 2>&1
    set -- $(classes $out/tp-baseline.vcf.gz) $(classes $out/fp.vcf.gz) $(classes $out/fn.vcf.gz)
    printf '%s\t%s\tvcfeval_%s\tsnp\t%s\t%s\t%s\n' "$dataset" "$caller" $mode $1 $3 $5 >> $TSV
    printf '%s\t%s\tvcfeval_%s\tother\t%s\t%s\t%s\n' "$dataset" "$caller" $mode $2 $4 $6 >> $TSV
  done
}

T300=$(indexed $B/human_genome_bottle/results/HG002_bottle_reference_1000rg.vcf.gz truth/giab_1000)
TTR=$(indexed $B/ssr_hg002/truth/HG002_GRCh38_TandemRepeats_v1.0.1_50000.vcf.gz truth/tr)

score 300x "$T300" $W/sat_1000.sorted.bed gatk \
  "$(indexed $B/human_genome_bottle/results/gatk/single_HG002_reads_selected_1000_rg.vcf gatk_300x)"
for setting in default call_satellites; do
  score 300x "$T300" $W/sat_1000.sorted.bed pvc_$setting "$(indexed $W/pvc/300x_$setting.vcf pvc_300x_$setting)"
done
for depth in 30x 10x; do
  score tr$depth "$TTR" $W/sat_tr.sorted.bed gatk "$(indexed $W/gatk/tr$depth.vcf.gz gatk_tr$depth)"
  for setting in default call_satellites; do
    score tr$depth "$TTR" $W/sat_tr.sorted.bed pvc_$setting \
      "$(indexed $W/pvc/tr${depth}_$setting.vcf pvc_tr${depth}_$setting)"
  done
done
cat $TSV
