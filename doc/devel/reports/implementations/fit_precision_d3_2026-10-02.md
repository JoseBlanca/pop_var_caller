# Fit precision, step D3 — the tool that fits chosen strata of a real cohort on every sample and on the grown subset

**Date:** 2026-10-02. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step D3.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §4.5 item 3. **Branch:** `fit-precision`.

## 1. What was built

**`examples/ng_ssr_subset_against_every_sample.rs`**, for the owner to run on kimura. It opens a cohort exactly as
`estimate-parameters` does, takes its repeat-tract evidence, and fits chosen strata twice — on every sample
(`fit_stratum`) and as a run now does (`fit_strata_on_sample_subsets`) — printing, a stratum at a time:

- the samples the subset reached and both fits' times;
- the slippage level, the shorter share, the fall-off and the concentration of each fit, with their standard
  errors;
- how far apart each number is, in the standard deviation the difference has when the subset's samples are part of
  the whole, `√(subset error² − whole error²)`: about one is what drawing a subset does.

```text
cargo run --release --example ng_ssr_subset_against_every_sample -- \
    --reference ref.fa --psp cohort/ --output unused.toml --inbreeding 0.9
```

- **The flags are `estimate-parameters`'s.** `--output` is required by them and not written to.
- **The homozygote excess** is `--inbreeding` for every sample, or 0. The SNP/indel fit that measures it per sample
  is not run (hours on kimura); both fits of a stratum use the same value.
- **Which strata**: `STRATA=1:8,2:10`, or five spread across the cohort's sizes — the largest and those an eighth,
  a quarter, a half and three quarters of the way down — so both the costly strata and the ones whose subset
  grows are shown.
- **The subset fit runs first and prints at once**; a line before the fit on every sample gives the likelihood
  table it will hold. On kimura that fit is the long part: spec §1.1 records seven strata taking 11 h 20 min under
  the old five-round limit.
- When the subset grew to every sample, both fits read the same samples and "apart" is given in the whole fit's
  errors instead. The last line totals both fits' time.
- `--inbreeding` is checked as `estimate-parameters` checks it.
- `FIRST_SUBSET`, `MIN_SAMPLES_A_GROUP` and `TARGET` override 256, 8 and 0.02, to try the tool on a small cohort.

To make that possible without copying a hundred lines of checks, two pieces of the run are now functions of their
own, each called by `estimate-parameters` as before:

- `cli::estimate_parameters::open_the_censuses` — the psps opened as a cohort, every census judged against the
  run's settings, the reference and catalog checked against the psps, the censuses assembled — with
  `contig_id_of`;
- `run::the_tract_strata_of_a_cohort` (and `refuse_another_selection`, the census-digest check `fit_a_cohort`
  makes first) — every stratum's evidence.

No fitted number moves.

## 2. What was run

On the four-accession tomato cohort's psps (`tmp/fit_precision/d3_smoke4.log`), with `FIRST_SUBSET=2
MIN_SAMPLES_A_GROUP=1 TARGET=0.5 STRATA=1:10,1:11 --inbreeding 0.9`: both strata stopped at a subset of 2 samples,
and their levels were −0.76 and +1.69 of the difference's spread from the fit on all four. With the defaults on
four samples the subset holds every sample, and the two fits print the same numbers to the last digit
(`d3_smoke.log`). Four samples show the plumbing works and nothing about what a subset of 256 in 2,169 loses — that
is the kimura run, checkpoint D's.

Review ([fit_precision_d3_2026-10-02.md](../reviews/fit_precision_d3_2026-10-02.md)): the two functions split out
change no behaviour, checked line by line and by `run::census_fit` (9 tests) and `cli::estimate_parameters` (21).

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite on the reviewed code: 5,025 passed, 3 failed (the pre-existing `examples/ng_generic_loci_dump.rs` × 2
  and `examples/ng_ssr_loci_dump.rs` × 1), 9 ignored. The review's fixes touched only the example.
- The oracle cohort (`scripts/promote_ng_oracle.sh`): every checksum matches the baseline.
