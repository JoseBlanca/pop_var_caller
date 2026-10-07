# Fit precision, step D2 — a large cohort's repeat-tract strata read from a growing subset of samples

**Date:** 2026-10-02. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step D2.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §4.4, as amended with the owner at this step.
**Branch:** `fit-precision`. Review: [fit_precision_d2_2026-10-02.md](../reviews/fit_precision_d2_2026-10-02.md),
fixes: [fixes_applied_fit_precision_d2_2026-10-02.md](../reviews/fixes_applied_fit_precision_d2_2026-10-02.md).

## 1. What was built

Each repeat-tract stratum is fitted on the reads of its samples, and the cost grows with the samples: at kimura's
2,169 samples the largest stratum's likelihood table alone is 5.9 GiB (spec §4.4). This step fits each stratum of
a large cohort on a subset of its samples, grown only until the stratum's slippage level is measured well enough.

- **`fit_strata_on_sample_subsets(strata, homozygote_excess, sample_order, config)`**, called by
  `run::census_fit` with step D1's order of the cohort's samples. A cohort of at most `SampleSubsets::first`
  samples (256) goes to `fit_strata` unchanged, so its output is the same bytes as before.
- **Above 256, each stratum** (`fit_on_growing_subsets`):
  1. is refused or not on its **whole** evidence, before any subset is drawn, by the refusal floor as before;
  2. takes the first 256 samples of the order, and every slippage group with reads in the stratum is topped up
     to 8 of its readers (`MIN_SAMPLES_A_GROUP`), or all of its own (`grow_subset`). A larger subset holds the
     smaller one's samples, doubles it, and takes every sample once it would hold more than three quarters of
     them (`next_subset_size`);
  3. a subset with fewer tracts with reads than the refusal floor is not fitted; it grows;
  4. the first subset fitted is fitted from every starting point; each later one by one walk from the last
     answer, which decides only whether to grow;
  5. the subset stops growing once every live slippage group that still has readers outside it has a level
     with a standard error below 0.02 of itself (`LEVEL_RELATIVE_ERROR_TARGET`, `level_is_measured_to`), or
     once it holds every sample. **The subset the answer is taken from is fitted from every starting point and
     from the last answer as well, the best winning.**
- **The answer's evidence counts are the subset's** (tracts with reads, reads crossing), so the curves across
  strata weigh the stratum by the reads it was fitted on. The substitution counts are still summed over every
  sample, and the guard's counts are the whole stratum's.
- **`StratumFit::samples_fitted_on`**: the subset's size; `None` when no subset was drawn. The run's log gains
  one line: the fewest, the median and the most samples a stratum was read from, and how many grew to every
  sample. Writing it in the parameters file is added to plan step E2.
- `SsrFitConfig::subsets` (`SampleSubsets`) holds the three settings, so tests can draw small cohorts;
  `fit_strata` ignores it.
- A stratum that grows to every sample is fitted on its own evidence, not on a copy of it.
- Under `--str-param-estimates-at-once N > 1` the strata go to the pool, each on one thread: a stratum's subset
  fits follow one another. The default, one stratum at a time, is unaffected.

## 2. Decided with the owner

The review measured the rule as first built and found six things the spec had left open or had wrong; the owner
accepted all six (2026-10-02), and spec §4.4 is amended with them:

1. the answer's subset is fitted from every starting point and from the last answer — on one drawn stratum of
   thirteen classes, one walk from the last answer had ended 1.3 to 41 log-likelihood units below three fresh
   starts on the same samples, in all four comparisons;
2. every slippage group is topped up to 8 readers, not only one with none;
3. the target is judged only on groups that still have readers outside the subset;
4. every sample is taken once doubling would pass three quarters of them;
5. a subset thinner than the refusal floor grows;
6. the parameters file is to record how many samples each stratum was fitted on (plan step E2).

Two consequences are recorded in the spec: a well-read stratum's curve weight, and the one-thread-a-stratum
schedule above.

## 3. What was measured

**The subset against every sample** (`the_subset_against_every_sample`, ignored; log
`tmp/fit_precision/d2_subset_measure.log`). Strata drawn at a level of 0.05, fitted on every sample and on the
grown subset. "Apart" is the subset's level minus the whole fit's, in the standard deviation that difference has
when the subset's samples are part of the whole: √(subset error² − whole error²).

| stratum drawn | draws | subset reached | apart | time, subset against whole |
|---|---|---|---|---|
| 3 classes, 300 tracts × 1,024 samples × 3 reads | 3 | 256 each | −0.34, −0.20, −1.09 | 0.20 to 0.37 |
| 3 classes, 1,000 × 1,024 × 3 | 3 | 256 each | +0.98, −1.30, −0.43 | 0.16 to 0.25 |
| 3 classes, 300 × 1,024 × 30 | 3 | 256 each | 0.00, −1.94, −0.55 | 0.19 to 0.32 |
| 13 classes, 200 × 600 × 3 | 3 | 256 each | −3.71, −0.99, −2.58 | 0.28 to 0.54 |

- **Every stratum stopped at its first subset of 256**, its level's error under 2% of the level, in 16 to 54% of
  the time the whole fit took.
- **At three classes the subset's level differs from the whole fit's by what drawing a quarter of the samples
  explains**: within 1.94 of that difference's standard deviation in all nine draws.
- **At thirteen classes the subset's level was below the whole fit's in all three draws**, by 1.0 to 3.7 of that
  standard deviation. Both were above the truth of 0.05, and the subset was the closer: 0.0543 to 0.0546 against
  0.0548 to 0.0564. The 256-point average's error a tract grows with the cohort
  ([its measurement](fit_precision_quadrature_average_2026-10-01.md)), which may be why the fit on more samples
  sits further from the truth; that is not measured here.
- These runs were made with the rule as first built. On these strata the final rule takes the same path — one
  slippage group, its readers in the first 256 well over 8, 256 below three quarters of the cohort, and no
  growth — so the numbers are the final code's.
- **Cost of the final refit, from the review's measurement**: a stratum that grows to every sample of 2,169 takes
  about twice the walking it took before D2; one that stops at 512, under half.

**The four-accession oracle cohort** (`scripts/promote_ng_oracle.sh`): four samples, so no subset is drawn —
every checksum matches the baseline.

## 4. Tests

| test | what it shows |
|---|---|
| `a_cohort_no_larger_than_the_first_subset_is_fitted_whole` | with a first subset of 20, twenty samples give `fit_strata`'s outcome and no subset; twenty-one take a subset |
| `a_subset_grows_until_the_level_is_measured_or_every_sample_is_in` | 64 samples, first subset 8: an unreachable target takes all 64 and four walks, scoring at least what three starts on every sample score; a target the first fit meets stops at 8 from three walks, with fewer reads crossing |
| `the_subset_doubles_until_three_quarters_of_the_cohort` | 256, 512, 1,024, then 2,169; at 64, 16, 32, then 64 |
| `a_subset_tops_up_a_thin_group_and_keeps_its_samples_as_it_grows` | three readers of a group in the first eight: five more added; doubling keeps every sample |
| `the_floor_is_judged_on_the_whole_stratum_and_a_thin_subset_grows` | forty tracts each read by one sample: fitted on the 8 samples that read 8 tracts, not refused; six such tracts refused below the floor |
| `the_subset_is_the_same_samples_whatever_order_they_arrive_in` | indices shuffled and the evidence relabelled: the same names, three readers of a missing group added |
| `the_two_schedules_give_the_same_bits_on_subsets` | one stratum at a time and two at once, on subsets: the same bits |
| `the_level_is_judged_only_in_groups_that_can_still_grow` | two groups at 1% and 5%: the 2% target met only when the second has no readers left outside; a growing group without an error has not met it |
| `one_read_makes_a_sample_a_reader_of_its_group` | one read in a bucket makes a reader; empty buckets do not |
| `the_subsets_summary_counts_the_sizes` | the log line's fewest, median, most and every-sample count |
| `the_subset_against_every_sample` (ignored) | §3 |

Seven mutations of the final code, each run against the module's tests and restored (`tmp/fit_precision/d2_mut/`):

| mutation | result |
|---|---|
| the floor not judged on the whole stratum | killed (`the_floor_is_judged_…`) |
| a subset thinner than the floor fitted anyway | survived at first — the test's top-up already reached 8 tracts; with the test topping up to one reader, killed |
| `<` for `<=` at the first subset's size | killed (`a_cohort_no_larger_than_the_first_subset_…`) |
| a group topped up only when it holds none | killed (two tests) |
| each subset drawn afresh, not added to | survived — and changes nothing found: with every group topped up to 8, a group's kept readers are always the earliest of its readers in the order, so the subsets nest either way; the cumulative flags are a guard |
| the answer's subset not refitted from every start | killed (`a_subset_grows_…`) |
| the target asked of every live group | killed (`the_level_is_judged_only_in_groups_that_can_still_grow`) |

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Module (`cargo test --release --lib parameter_estimation::joint::ssr_fit`): 71 passed, 4 ignored.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 5,025 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 9 ignored; the cross-platform
  checksums pass unchanged.
- The oracle cohort (`scripts/promote_ng_oracle.sh`): every checksum matches the baseline.
