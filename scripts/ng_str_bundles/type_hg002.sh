#!/usr/bin/env bash
# Type HG002 tandem-repeat benchmark regions with this checkout's ng_typed_region_dump (built first:
# `cargo build --release --example ng_typed_region_dump`), writing <checkout>/tmp/bundles/typed.tsv.
#   DEV_EXTRA_MOUNT=<main checkout>/benchmarks scripts/dev.sh bash scripts/ng_str_bundles/type_hg002.sh
set -euo pipefail
w=$(cd "$(dirname "$0")/../.." && pwd)
data=/Users/jose/devel/pop_var_caller/benchmarks
ref=$data/giab/ref_genome_GRCh38/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
mkdir -p $w/tmp/bundles
sort -k1,1 -k2,2n $data/ssr_hg002/regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed | cut -f1-3 > $w/tmp/bundles/tier_sorted.bed
NG_REFERENCE_CHECK=skip $w/target-container/release/examples/ng_typed_region_dump $ref $w/tmp/bundles/tier_sorted.bed calling 2>/dev/null > $w/tmp/bundles/typed.tsv
