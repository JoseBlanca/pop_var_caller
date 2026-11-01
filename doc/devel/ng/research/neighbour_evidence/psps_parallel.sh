#!/usr/bin/env bash
# psps for every CRAM in <cram dir> into <psp dir>, <jobs> at a time. Runs in the dev container.
set -euo pipefail
WT=/Users/jose/devel/pop_var_caller-neighbour-hmm
T1=/Users/jose/devel/pop_var_caller/benchmarks/tomato1
REF=/Users/jose/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa
CAT=$WT/tmp/tomato/S_lycopersicum_chromosomes.4.00.repeats.parquet
mkdir -p "$2"
ls "$1"/*.bench.cram | xargs -P "$3" -I{} sh -c \
  "$WT/target-container/release/pvc generate-psps --quiet --reference $REF --catalog $CAT --regions $T1/regions.bed --alignment {} --output-dir $2/\$(basename {} .bench.cram)"
