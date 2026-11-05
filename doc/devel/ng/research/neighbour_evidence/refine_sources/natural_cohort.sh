#!/bin/bash
# usage: [DEPTHS=30:10,5:10,2:10,1:10] [OUT=dir] natural_cohort.sh <name> <inbreeding>   — 40-sample msprime cohort, 10 each at 30/20/10/5x,
# 4 Mb over 2 contigs: simulate, map, catalogue, psps, fitted parameters. Runs inside the container.
set -euo pipefail
P=$(cd "$(dirname "$0")/../../../../../.." && pwd)  # the project root
name=$1; F=$2
D=${OUT:-$P/tmp/refine_sources}/$name; mkdir -p $D/sim $D/map $D/psps
uv=$P/benchmark_results/inputs/uv
if [ ! -s $D/sim/truth.vcf.gz ]; then
  UV_CACHE_DIR=$uv/cache UV_PYTHON_INSTALL_DIR=$uv/python uv run --no-project --python 3.12 --with msprime==1.3.3 \
    python $P/scripts/msprime_haplotypes.py $D/sim/haplotypes.tsv --samples 40 --contigs 2 --contig-bp 3000000 --inbreeding $F
  $P/target-container/release/examples/ng_str_stutter_benchmark $D/sim genome_bp=4000000 contigs=2 samples=40 \
    depths=${DEPTHS:-30:10,20:10,10:10,5:10} haplotypes=$D/sim/haplotypes.tsv inbreeding=$F threads=6 > $D/sim/simulate.log
  samtools faidx $D/sim/reference.fa
  bgzip -f $D/sim/truth.vcf && tabix -p vcf $D/sim/truth.vcf.gz
fi
for i in $(seq -f %03g 0 39); do
  s=sim$i
  [ -s $D/map/$s.bam.bai ] && continue
  minimap2 -ax sr -t 6 -R "@RG\tID:$s\tSM:$s\tLB:$s\tPL:ILLUMINA" $D/sim/reference.fa $D/sim/${s}_R1.fastq.gz $D/sim/${s}_R2.fastq.gz 2>/dev/null \
    | samtools sort -@ 2 -o $D/map/$s.bam - && samtools index $D/map/$s.bam && rm -f $D/sim/${s}_R?.fastq.gz
done
pvc=$P/target-container/release/pvc
[ -s $D/repeats.parquet ] || $pvc repeat-catalog --reference $D/sim/reference.fa --output $D/repeats.parquet --log $D/catalog.log --force
ls $D/map/*.bam | xargs -P 6 -I{} sh -c "f={}; s=\$(basename \$f .bam); [ -s $D/psps/\$s.psp ] || $pvc generate-psps --reference $D/sim/reference.fa --catalog $D/repeats.parquet --alignment \$f --regions $D/sim/confident.bed --output-dir $D/psps --force > /dev/null 2>&1"
psps=(); for f in $D/psps/sim*.psp; do psps+=(--psp $f); done
[ -s $D/fitted.parameters.toml ] || $pvc estimate-parameters --reference $D/sim/reference.fa --catalog $D/repeats.parquet "${psps[@]}" --output $D/fitted.parameters.toml --log $D/fit.log --threads 6 --force > /dev/null
echo "$name ready"
