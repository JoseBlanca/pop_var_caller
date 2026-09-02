# Shared benchmark runners

Caller runners + helpers shared across benchmarks. Each benchmark
(`benchmarks/<name>/`) contributes only a `bench.config.sh` describing
its paths and knobs; the scripts here hold all the common scaffolding
(binary discovery, preflight, timing, record counting).

## Layout

```
benchmarks/
  lib/
    common.sh             # sourced by the runners; config + helpers
    run_ours.sh           # pop_var_caller   (single | cohort)
    run_gatk.sh           # GATK             (single | cohort)
    run_freebayes.sh      # freebayes        (single | cohort)
    run_ng.sh             # ng, the experimental caller (single | cohort)
    prepare_reference.sh  # build .fai (+ .dict for GATK) for a reference
    compare_to_truth.sh   # precision/recall/F1 vs a truth VCF (accuracy)
  <name>/
    bench.config.sh       # per-benchmark paths/knobs (see comments inside)
    crams/                # input CRAMs (+ .crai)
    results/              # caller outputs land here (results/<caller>/...)
```

## Usage

```sh
# 1. one-time: make sure the reference has its index siblings
benchmarks/lib/prepare_reference.sh benchmarks/<name>/bench.config.sh

# 2. run a caller — mode is `single` (one sample) or `cohort` (all)
benchmarks/lib/run_ours.sh      benchmarks/<name>/bench.config.sh single
benchmarks/lib/run_gatk.sh      benchmarks/<name>/bench.config.sh cohort
benchmarks/lib/run_freebayes.sh benchmarks/<name>/bench.config.sh single
benchmarks/lib/run_ng.sh        benchmarks/<name>/bench.config.sh cohort

# 3. (accuracy benchmarks) score caller VCFs against a truth set.
#    With no VCF args it scores the four single-sample outputs.
benchmarks/lib/compare_to_truth.sh benchmarks/<name>/bench.config.sh
```

`DRY_RUN=1` makes the runners print the exact command they would run
(shell-quoted) instead of executing — handy for checking what a config
resolves to without needing every tool installed. Every config value is
also overridable from the environment (`REFERENCE=… THREADS=8 …`).

## Benchmarks

- **tomato1** — 63-accession *S. lycopersicum* cohort, CRAMs pre-sliced to
  the 80 intervals of `regions.bed` (8 Mb of SL4.0) at about three reads a
  position. Multi-sample and no truth set: it is where cohort behaviour and
  cross-caller agreement are checked, not accuracy. The perf experiments
  (`tomato1/scripts/perf_*.py`) build their PSP/GVCF inputs via the cohort
  runners.
- **human_genome_bottle** — single-sample GIAB HG002, CRAM restricted to
  the 1000-region (~5 Mb) benchmark BED on GRCh38. Accuracy benchmark:
  compare each caller against the GIAB truth VCF with `compare_to_truth.sh`.

## Notes

- GATK lives at `/opt/gatk/gatk` inside the dev container; override with
  `GATK_BIN=…`. `pop_var_caller` is auto-detected from
  `target-container/release` then `target/release` (override with
  `POP_VAR_CALLER_BIN=…`).
- `pop_var_caller` has no BED/region flag, so `run_ours.sh` processes the
  whole CRAM — fine here because the benchmark CRAMs are already
  pre-sliced to the region set. GATK, freebayes and ng restrict via the BED.
- **ng** is `pop_var_caller_exp call-from-alignments`, auto-detected the same
  way (override with `NG_BIN=`). It needs a tandem-repeat catalog built from
  the same reference: `run_ng.sh` builds one at `NG_CATALOG` if it is missing,
  which takes about 100 s on GRCh38. Where the reference sits on a read-only
  mount, the benchmark config points `NG_CATALOG` somewhere writable.
- **Two limits to hold against ng's numbers.** It does not call inside tandem
  repeats — every tract in the BED is counted as ground it cannot speak for and
  named in the run report, which `run_ng.sh` tees to the log — and nothing is
  fitted, because no command writes a fitted parameters file yet, so every run
  is `--defaults`. Both cost recall, not precision.
- `compare_to_truth.sh` scores **allele** concordance (POS+REF+ALT), not
  genotype concordance; for GT-level eval use rtg vcfeval or hap.py.
