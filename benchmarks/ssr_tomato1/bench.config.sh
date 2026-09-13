#!/usr/bin/env bash
# Benchmark config: ssr_tomato1 — SSR (microsatellite) genotyping on the same
# ~63-sample S. lycopersicum cohort as the SNP `tomato1` bench. Parallel to the
# SNP benches. The production STR caller that this config first drove, and its
# runner `run_ours_ssr.sh`, were deleted in promotion Milestone D; the config
# stays for `run_hipstr.sh`, and its results stay under results/ours/. Every
# value is env-overridable.
#
# NOTE: the cohort CRAMs are pre-sliced to an SSR-targeted region set.
#
# BED here is NOT a caller flag — it is the slice the CRAMs were cut to
# (ssr_regions.bed: ~15k catalog SSRs ±500 bp, ch00 excluded; built by
# scripts/pick_ssr_regions.py, sliced on rick via lib/slice_crams_on_rick.sh).
# run_hipstr.sh uses it to restrict HipSTR's loci to the covered set so both
# callers genotype the same loci. Regenerate the crams if you change it.
BED="${BED:-ssr_regions.bed}"

BENCH_NAME="ssr_tomato1"

# Reference: the SAME SL4.0 the cohort CRAMs were aligned to. The catalog's
# whole-reference MD5 is validated against every .ssr.psp, so the catalog and
# the pileups must use this exact FASTA. A sibling .fai is required.
REFERENCE="${REFERENCE:-$HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa}"

# SSR-targeted cohort CRAMs (gitignored local data), sliced on rick to
# ssr_regions.bed and copied here. Self-contained — no longer borrows the SNP
# bench's CRAMs (those covered only ~774 SSRs; this set covers ~15k).
CRAM_DIR="${CRAM_DIR:-$BENCH_DIR/crams}"
CRAM_GLOB="${CRAM_GLOB:-*.bench.cram}"
CRAM_SUFFIX="${CRAM_SUFFIX:-.bench.cram}"
# First cohort CRAM by name; overridable. (single mode is unused for SSR — the
# cohort pre-pass needs multiple samples — but bench_single_cram still expects it.)
SINGLE_CRAM="${SINGLE_CRAM:-}"

PLOIDY="${PLOIDY:-2}"
THREADS="${THREADS:-4}"

# CATALOG defaults to $OUT_ROOT/ours/$BENCH_NAME.ssr.catalog inside run_hipstr.sh
# (OUT_ROOT is resolved after this file is sourced); set CATALOG to override.
