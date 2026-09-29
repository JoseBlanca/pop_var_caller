# Fit precision, step B1 — the SNP/indel fit stops when every parameter is within a tenth of its error of the maximum

**Date:** 2026-09-29. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step B1.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §2 and §3.3 as amended at this step, §3.4.
**Branch:** `fit-precision`. **Review:** [fit_precision_b1_2026-09-29.md](../reviews/fit_precision_b1_2026-09-29.md);
**fixes:** [fixes_applied_fit_precision_b1_2026-09-29.md](../reviews/fixes_applied_fit_precision_b1_2026-09-29.md).
**Why the rule is not the one first specified:** [fit_precision_b1_stopped_2026-09-29.md](fit_precision_b1_stopped_2026-09-29.md).

## 1. What was built

The SNP/indel fit used to stop when no parameter moved by more than a thousandth of itself over a cycle. On
kimura's 2,169 samples that never held, and every start ran to the 200-pass limit. **It now stops when every
parameter is within a tenth of its own standard error of the likelihood's maximum.**

- **When it judges.** Once a cycle gains less than 10⁻⁴ log-likelihood units a position (the old log-likelihood
  rule, kept as the trigger), every later cycle's first pass also sums the information and each parameter's score,
  at the cycle's starting parameters (`InformationSums::slope`, new).
- **How far each parameter is from the maximum** (`newton_step`, `fit/standard_errors.rs`). The two sums describe
  the log-likelihood near the parameters as a quadratic, `g·d − ½ dᵀ I d`, and its maximum within the intervals
  the fit keeps the parameters in is the distance. Two solves are used, the same matrix the errors come from:
  - the whole matrix up to 20 samples and 188 parameters;
  - above that, the arrow of blocks: each sample's own parameters given the cohort's, then the cohort's from
    what the samples leave of its slope.

  **The step stays inside the bounds** (`bounded_step`), by the active-set rule for a box:
  - a parameter whose step would carry it past an end stops there, and its distance is the way there;
  - the rest are solved again with its move taken out of their scores;
  - a stopped parameter whose slope at the solution points back into its interval is freed.
- **When it stops** (`fit/settled.rs`). A parameter is settled when its distance is below `settled_fraction`,
  0.1, of its standard error from the same pass. One without an error or a distance is settled by definition.
  **The fit has converged at the end of a cycle whose first pass found every parameter settled.**
- **What the run's log says.** Each cycle's progress line gives how many parameters are not yet within a tenth of
  an error of the maximum and names the furthest, with its distance in errors. Before judging starts, it gives
  the threshold at which it will. The start line compares its final pass with the average *plain* pass, the
  in-loop information passes left out.

`JointFitConfig::stillness` is now `settled_fraction`. `largest_relative_move` and its scale floors are gone.

## 2. How the rule came to be Newton's

Built as spec §2 first worded it, the rule projected each parameter's remaining distance from its last two moves.
On the four-accession oracle cohort it stopped every start at 36 to 57 passes, 74 log-likelihood units below a fit
run to 1,000 passes; the invariant share was 28 errors from where that fit put it. The stop report has the
measurements. The owner approved the Newton distance (2026-09-29). On 15 drawn cohorts, before the review's bound
fix, every fit the Newton rule called converged lay within 0.088 of an error of a 600-pass fit.

Review then found that a fit whose maximum lies on a bound never settled. Two cases:
- a clean error rate the golden section leaves 4.4 × 10⁻¹⁰ above its floor;
- a density shape the fit walks towards its floor.

The bounded step fixes both (§4).

## 3. What moves

- **The oracle cohort** (four tomato accessions at about 3 reads, `tmp/fit_precision/oracle_b1f/`): **nothing.**
  - Every start still runs to the 200-pass limit (198, 198, 200 passes). The fit's two checksums are unchanged:
    `fitted.parameters.toml` `393e7cd3…`, `from_fit.comparable` `1f16dec8…`.
  - Its log now says why. Start 3 at pass 170: 5 of 20 parameters not yet within a tenth of an error, the furthest
    the duplicated share at 0.38 errors. Start 2 at pass 176: 7, the invariant share at 2.85.
  - A fit run to 1,000 passes gains a further 1.5 units and is still gaining (0.2 units from pass 700 to 999). This
    cohort's density is barely determined: the first shape's error at the maximum is ± 36.6.
  - Two window-coverage checksums differ from `scripts/promote_ng_oracle.baseline`. That comes from the merge of
    main (its psp format 1.1 counts the reads the depth cap discarded into window depth). It is re-recorded in its
    own commit.
- **The cross-platform checksums** (`cli::cross_platform_digests`), re-recorded with the explanation in their doc.
  - Two samples over 600 bases give almost no parameter a standard error, and those without one are settled by
    definition, so every start now converges after 18 passes where it ran to the limit.
  - The density's first shape stops at 1.546, where 200 passes took it to 0.288.
  - Two lines of the parameters file move: the ordinary-site prior's reference concentration 39.089 → 20.047, and
    its alternative total 6.60 × 10⁻³ → 2.08 × 10⁻².
  - In the calls three records' quality moves: QUAL 541.9 → 549.5 and 542.7 → 550.3 at the two SNPs, 7.4 → 14.3 at
    the repeat tract, and one sample's GQ 73 → 71 at each SNP. No genotype moves.
  - The digests were `9d421df0…`/`7a570c44…`; they are now `22dda760…`/`2c0a8946…`.

## 4. What was measured

### 4.1 The 200-cohort coverage test

`the_errors_mean_what_they_say`, `tmp/fit_precision/b1/coverage_b1c.log`, 2,139 s, against A8's run:

| regime | fits converged (A8) | passes a start, median (A8) | within one error, by kind | flat point, errors (A8) |
|---|---|---|---|---|
| 2 samples × 3 reads | 8 of 200 (146) | 198 (161) | 0.618–1.000 | −0.31 to +0.28 (−0.35 to +0.32) |
| 4 × 3 | 182 (198) | 93 (84) | 0.660–0.810 | −0.064 to +0.027 (−0.108 to +0.029) |
| 4 × 30 | 200 (200) | 42 (45) | 0.655–0.690 | −0.035 to 0 (−0.041 to 0) |
| 20 × 3 | 200 (200) | 18 (24) | 0.650–0.735 | −0.009 to 0 (−0.002 to 0) |
| 20 × 30 | 200 (200) | 15 (18) | 0.665–0.722 | −0.008 to 0 (−0.002 to 0) |

- **The errors still mean what they say** in every regime: the share within one error per kind differs from A8's
  by at most 0.02 (the invariant share at 2 samples), and by at most 0.015 elsewhere.
- **At 20 samples the fits take about a quarter fewer passes.** They stop within 0.009 errors of the likelihood's
  flat point, against 0.002, inside the tenth the rule allows.
- **At 2 and 4 samples × 3 reads fewer fits converge.** The rule no longer calls a fit converged while it crawls
  across a flat likelihood: those fits run to the limit and say so. On the reviewed commit the tests reviewer re-ran
  three such cohorts to 2,000 passes: one converged at 219, one at 876, and one had not converged at 1,999.
- The test's convergence check is re-recorded from this run: at 4 samples × 3 reads, at least 175 of 200.

### 4.2 Fits whose maximum lies on a bound

`a_fit_whose_maximum_lies_on_a_bound_converges`, the tests reviewer's six cohorts:

| cohort | what sits on a bound | passes |
|---|---|---|
| 20 samples × 3 reads, 4 × 3, 4 × 30 | the homozygote excesses, at 0 | 21, 57, 27 |
| 20 × 3, 4 × 8 | the clean error rates, drawn at 10⁻⁷ below their floor of 10⁻⁶ | 15, 45 |
| 20 × 3 | the density's second shape, drawn at 150 above its bound of 50 | 36 |

All six converge. On the reviewed commit, which held a parameter only exactly on an end, the last three ran to the
200-pass limit.

### 4.3 The stop against a long fit

`a_settled_fit_lands_within_a_tenth_of_an_error_of_a_long_one`: at 20 samples × 3 reads the fit settles in 18
passes, its furthest parameter 0.0043 errors from a 600-pass fit. At 4 samples × 8 reads with the duplicated class
it settles in 105 passes, 0.063 errors away.

### 4.4 Cost

The tests reviewer measured a pass that also sums the information at 1.48 to 1.80 times a plain pass, at 4 to 64
samples (one thread, median of 9 interleaved runs). About 30% of a small cohort's passes carry it once judging
starts, so a whole fit takes 11 to 13% longer. `newton_step` itself takes 1.6 ms at 2,000 one-library samples.

### 4.5 Determinism

The summed slope, the Newton step and a whole judged fit give the same bits at 1, 4 and 8 threads (review probe:
3, 25 and 4 samples with the duplicated class, over 13 chunks). The committed thread-width test now also runs a
fit whose stop is decided by an in-loop information pass. The solve uses only + − × ÷ and sqrt.

## 5. Tests and validation

| test | what it shows |
|---|---|
| `settled::tests` (6) | the settled count at and below the fraction, without an error or a distance, at a fraction of zero, a NaN distance; the furthest parameter |
| `on_the_whole_matrix_the_newton_step_is_the_dense_solve`, `on_the_arrow_…` | free steps equal a dense solve to 10⁻⁹ |
| `a_parameter_carried_past_an_end_stops_there` | a share, a shape and an excess, both matrices, both ends: stops at the end, the rest re-solved with its move taken out |
| `a_parameter_on_an_end_stepping_inwards_is_free` | an end with an inward step leaves the dense solve |
| `a_rate_the_golden_section_leaves_beside_an_end_stops_there` | the M-step's resting point beside either end gives a distance below 10⁻⁹ |
| `at_one_sample_the_newton_step_leaves_the_held_excess_out` | the excess held fixed is not solved for, on both matrices |
| `check_a_pass_against_its_scores` (extended) | the summed slope equals the scores summed by hand to 6.4 × 10⁻¹⁶ |
| `a_fit_whose_maximum_lies_on_a_bound_converges` | §4.2 |
| `a_settled_fit_lands_within_a_tenth_of_an_error_of_a_long_one` | §4.3 |
| `a_census_of_many_chunks_fits_to_the_same_bits_at_any_pool_width` (extended) | a judged stop is the same bits at one thread and four |

- The fit module (`cargo test --release --lib parameter_estimation::joint::fit`): 94 passed, 1 ignored.
- `cargo fmt --check` and `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- The full suite: in the commit message; the three pre-existing failures remain.

## 6. Deviations and what is left

1. **The rule is not the first-specified one** (§2): the owner approved the change and the spec is amended.
2. **The bounded step** is review's fix. It makes "held at a bound" hold for a parameter approaching one, and the
   spec amendment says so.
3. **Judging starts at the old log-likelihood rule** and never stops once started; the review found no fixture on
   which starting earlier changes a verdict, only the cost.
4. **Left for later** (fixes report §5):
   - a parameter without information or not identified inside the solve's tests, to B2;
   - one shared type for the errors' inversion and the solve, for the owner (about half a day);
   - renaming `log_likelihood_stillness`;
   - a loop-level test that stopping needs every parameter settled, with B3's per-cycle record.
5. **For B2:** at 2 and 4 samples a full Newton step from part-way through a fit points the right way but is
   several times too long (review: it lands 39 to 331 log-likelihood units below the 600-pass fit). "The current
   value plus the step" as a start's projected endpoint needs a guard there.
