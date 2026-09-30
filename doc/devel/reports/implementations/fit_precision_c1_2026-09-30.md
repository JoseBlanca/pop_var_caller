# Fit precision, step C1 — each repeat-tract stratum's fit carries the standard errors of its numbers

**Date:** 2026-09-30. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step C1.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §4.2, and §3.2's rule for an error wider than its
range. **Branch:** `fit-precision`. **Review:** [fit_precision_c1_2026-09-30.md](../reviews/fit_precision_c1_2026-09-30.md);
**fixes:** [fixes_applied_fit_precision_c1_2026-09-30.md](../reviews/fixes_applied_fit_precision_c1_2026-09-30.md).
This report describes the step as committed, after the review's fixes.

## 1. What was built

The repeat-tract fit estimates, for each stratum (every tract of one motif length and reference repeat
count), three slippage numbers a slippage group — how often a read slips, which way, how fast two-unit slips
fall off — and a length spectrum and a concentration describing the stratum's alleles. **Each stratum's fit
now also says how precisely its own tracts determine each of those numbers, or why they do not.** No number
moves.

- **`StratumFit::standard_errors: Option<StratumErrors>`** — per slippage group the three numbers' errors
  (`SlippageErrors`), per allele class its share's, and the concentration's, each a `StratumError`:
  - `Estimated(error)`, on the number's own scale;
  - `NotIdentified` — once the stratum's other numbers are accounted for, the log-likelihood is not curved
    downwards in its direction (flat, or a saddle);
  - `NotPlaced` — its error on the scale the climb moves it on is above `NOT_PLACED` = 3 units of logit or
    log, so a one-error interval spans odds or a size about 400 times apart; or, for a number in `[0, 1]` (a
    slippage number or a share), its error on its own scale is wider than that whole range — spec §3.2's rule;
  - `NoShare` — a length class with no share at all, held at zero.

  The field is `None` only on fixtures built by hand. **The errors are the stratum's own fit's**: after
  `fit_strata` draws the curves, `slippage` beside them may hold a blend, which has no error (spec §5.2).
- **How** (`standard_errors_at`, spec §4.2):
  - the second derivatives of the stratum's **total** log-likelihood — the mean a tract that the climb
    maximises, times every tract it is a mean over — by central differences (`curvature_of`), at the answer
    the winning walk returned;
  - over every number the stratum fits, on the climb's scales: logit for the slippage numbers, log for the
    concentration, and for the spectrum the log of each class's share over the largest class's (the shares
    sum to one, so one class is not a number of its own); one step, 10⁻², on every coordinate
    (`CURVATURE_STEP`), those scales being relative already;
  - the negative of that curvature inverted by the SNP/indel fit's `invert_identified`, which drops a number
    keeping less than 10⁻⁸ of its own curvature once the others are accounted for;
  - carried to each number's own scale by its derivative: `p(1 − p)` for a slippage number, the
    concentration itself, and the shares' slopes in the log-ratios for each share.
- **Why `NotPlaced`** (found in review). A golden-section climb moves a number a bounded span a round, so a
  number whose likelihood keeps rising towards an end of its range stops where five rounds took it. There the
  likelihood is flat, and `p(1 − p)` turns a logit error of about 1,000 into a tiny one: a shorter share of
  0.9999997 came back ± 0.0003 where moving it to the true 0.83 costs 0.646 log-likelihood units, and a
  concentration run off to 8 × 10⁶ came back with a curvature that was rounding noise. Numbers the tracts
  determine sat at 0.02 to 0.5 on the climb's scale, stuck ones at 1,000 to 2,400 (review measurements). This
  is the repeat-tract counterpart of spec §3.2's "wider than its range", which is also applied as written: on
  the oracle cohort, before it was, the largest class's share of one stratum came back with an error of 2.85,
  carried from other classes' log-ratios that were not placed.
- **Both thread schedules compute the same bits**, now checked by the schedule tests.
- **The run's log** prints, once every stratum has its answer: over the strata fitted on their own tracts, the
  median and largest error of the level and of the concentration as a percentage of the number, of the two
  shares, the fall-off and the class shares as they are, and the numbers without an error counted by why.
- The period-pooled length spectrum (the tract prior's middle rung) is fitted without errors; nothing reads
  them there.

## 2. Choices the spec left open

1. **The step**, 10⁻² on every coordinate of the climb's scales. At a third of it the errors of a drawn stratum
   agree to 2 in 100 (tested).
2. **The spectrum's log-ratios are taken against the largest class**; in exact arithmetic the choice changes no
   error.
3. **`NOT_PLACED` = 3**, a threshold between the two populations the review measured. For the owner at
   checkpoint C.
4. **The drop order.** When two numbers cannot be told apart, `invert_identified` drops the later in the
   layout — slippage numbers come first, then the classes, then the concentration — and the kept one's error
   is taken with the other held.
5. **`invert_identified` is shared** with the SNP/indel fit, its visibility widened to the joint-fit module.

## 3. What was measured

### 3.1 Cost

On drawn strata at the production span (13 allele classes, 16 numbers, 513 evaluations), one stratum at a time
(`tmp/fit_precision/c1/probe_cost.log`, before the review's fixes, which leave the number of evaluations unchanged):

| stratum | the climb, 3 starting points | the errors | share |
|---|---|---|---|
| 200 tracts × 20 samples × 3 reads | 14.72 s | 1.69 s | 11.5% |
| 1,000 × 20 × 3 | 62.82 s | 6.92 s | 11.0% |
| 500 × 63 × 3 | 96.83 s | 10.55 s | 10.9% |

The correctness reviewer measured 10.8 to 12.3% on others, thin ones included, with 19 likelihood-table and
501 quadrature rebuilds for the 513 evaluations. **The oracle cohort's repeat-tract fit time cannot show it**:
four runs of the unchanged stratum-fit code on this branch took 1 min 24 s to 1 min 42 s.

### 3.2 Whether the errors mean what they say — for C2

C2's coverage test is where that is measured. Seen so far:

- At 3 allele classes the level's error halves as the tracts quadruple, 0.00164 at 400 tracts to 0.00082 at
  1,600, and the fitted levels lie 2.3 and 1.0 of their errors from the truth. Over 20 redrawn strata of that
  regime, the fitted numbers' spread was 0.83 to 1.22 times the reported errors (correctness reviewer); an
  independent rebuild of the curvature agreed with the code's errors to 0.3%.
- **At 13 classes, the fitted level lay 3.0 to 4.5 of its errors above the truth on all three strata of §3.1**
  (0.0570 ± 0.0023, 0.0548 ± 0.0011, 0.0525 ± 0.0008 against 0.05). Either the errors are too small there or
  the climb, which stops at five rounds, stops short; C2 and C3 will tell.

### 3.3 What moves, and the four-accession oracle cohort

**No number.** `scripts/promote_ng_oracle.sh` (`tmp/fit_precision/oracle_c1g/`): every checksum matches its
baseline, the parameters file included. Its log's new line, over the 15 strata fitted on their own tracts:
the level ± a median 25% of itself (largest 100%), the shorter share ± 0.106 (largest 0.206), the fall-off ±
0.070 (largest 0.155), the concentration ± 31% of itself (largest 86%), the class shares ± 0.011 (largest
0.672); of their 255 numbers, 76 are not placed and 54 not identified. Four samples at about 3 reads place
half a stratum's numbers at most.

## 4. Tests

| test | what it shows |
|---|---|
| `central_differences_give_a_quadratics_second_derivatives` | on a quadratic of four variables, every second derivative to 10⁻⁸, in `1 + 2p²` = 33 evaluations |
| `the_errors_are_carried_to_each_numbers_own_scale` | a given inverse carried to the slippage numbers, the concentration and each share (the largest class's included) exactly; a dropped log-ratio leaves its share `NotIdentified` and the others with theirs |
| `a_stratums_errors_shrink_as_its_tracts_grow` | every number has an error; the level's halves (ratio 0.50) from 400 to 1,600 tracts — the mean's curvature would leave it unchanged |
| `a_stratums_errors_do_not_depend_on_the_difference_step` | at the default step and a third of it every error agrees to 2 in 100 |
| `a_number_the_climb_left_at_the_end_of_its_reach_is_not_placed` | three thin drawn strata: a fall-off stopped at 1.3 × 10⁻⁷, a shorter share at 0.9999997 and a concentration at 8.0 × 10⁶ are each `NotPlaced`, the stratum's level keeping its error |
| `tracts_without_reads_change_no_error` | 300 tracts no sample read, added, change no error beyond 10⁻⁶ of itself |
| `the_curvatures_centre_is_the_fitted_answer` | the layout's centre set back is the fitted answer to 10⁻¹² |
| `the_summary_gives_each_kinds_errors_and_counts_the_missing_by_why` | the log's line: percentages, medians, largest, the missing by why |
| `the_two_ways_of_spending_the_pool_give_the_same_bits`, `any_number_of_strata_at_once_gives_the_same_bits` | now compare every error, and require each fitted stratum to carry them |

- Repeat-tract fit module (`cargo test --release --lib parameter_estimation::joint::ssr_fit`): 42 passed.
- Four mutations the review found surviving now fail a test (`tmp/fit_precision/c1/mut/`): no errors in the
  several-at-once schedule; the tract count over tracts with reads; the spectrum not renormalised; the
  summary's level not a percentage.
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 4,992 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 5 ignored.

## 5. Deviations and what is left

1. **`StratumError` with its reasons, not the plan's `Option<f64>`**, so a number without an error says why;
   `StratumError::value` gives the `Option`. `NotPlaced` applies spec §3.2's range rule here, on the climb's
   scale.
2. **Deferred from review**: a test that the off-diagonal curvature matters (all of it set to zero survived,
   moving class shares up to 18%); the identification floor's calibration for finite differences; moving the
   new code to its own file and `invert_identified` to a shared module; starting a stratum's errors in the
   several-at-once schedule as soon as its own walks finish.
3. The errors reach the parameters file at E2; C2 measures their coverage.
