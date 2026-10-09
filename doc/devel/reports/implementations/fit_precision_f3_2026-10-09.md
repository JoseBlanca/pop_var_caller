# Fit precision, step F3 — the shapes' step capped, and a start returns its best point

**Date:** 2026-10-09 · **Branch:** `fit-precision-followup` · **Plan:**
[fit_precision.md](../../implementation_plans/fit_precision.md), step F3 and the two decisions recorded beside it.

## Why

The owner's traced kimura run (2,169 samples; report `tmp/fit_trace_kimura_2026-10-08.md` in the main checkout)
found the SNP/indel fit going round a 12-pass loop. Over 11 passes the log-likelihood climbs 1.58 million units
while the allele-frequency shape b grows about twofold a pass to 34.4; then the plain pass steadying a SQUAREM jump
leaves b at its lower bound, 0.02, and the 1.58 million is lost and kept. No pass limit finishes it, and the fit
returned start 2's last point, 128,600 units below the best point any start reached.

The shapes' update (`step_beta_shapes`) is one Newton step with the complete-data curvature, which in b shrinks
about as a/b² (2.35 × 10⁻⁴ against a/b² = 2.30 × 10⁻⁴ a position at a = 0.27, b = 34.4), so a modest slope asks
for a step larger than b. That this step made the loss is inferred from the code: the trace records where each
pass left the parameters, not which update in the pass moved them.

## What was built

- **No Beta shape moves more than a factor of two in a pass** (`MAX_SHAPE_FACTOR_A_PASS`), for the frequency
  density's shapes and the duplicated class's carrier shapes. The cap is a box the Newton step is kept inside, as
  the shapes' bounds already were: where one shape's step crosses the cap, it is held there and the other takes the
  step the same curvature gives with it held (`the_shapes_newton_step`, now taking a box a shape). A shape that is
  not a positive number moves nothing.
- **A start that stops at the pass limit returns the best point it measured**, a cycle's start or a jump that held,
  when its last point is more than `JUMP_SLACK` (one unit) below it, with one more final pass there; the start's
  log line says so and by how much. A start that converged or agreed returns where it did.
- **Dropped, owner, 2026-10-09: going back to the best point mid-fit.** Built and measured first, then removed;
  the plan records why.

## Measured

**The kimura loop's step, capped** (`no_shape_moves_more_than_a_factor_of_two_in_a_pass`): at a = 0.27, b = 34.4,
a slope pulling b a little down asks, uncapped, for a step past the floor (b = 0.02); capped, b halves to 17.2 and a
takes its step with b held there.

**The four-accession tomato cohort** (`scripts/promote_ng_oracle.sh`, 20 regions, 7,451 records): every checksum
matches the baseline, fitted file and calls included. Its three starts run to the pass limit and none returns a
best point more than a unit above its last.

**The 21-sample cross-platform fixture moves**, and only through the best-point return: switching the return off
reproduces its recorded checksums exactly, switching the cap off does not. All three starts stop at the limit below
their best, by 15.3, 1,004.3 and 535.8 units; the winner, start 1, now returns −5,099.1 against −5,114.4. 24 of
the file's 549 lines move, every library's multiplier among them (read group 0 from 4.85 to 7.66). Its calls: the
same 22 records and genotypes; `QUAL` lower by 0.7 to 78.4 at every site; 9 sample GQs move, by −4 to +8.
Re-recorded in the Linux container; the same on macOS (arm64).

**What remains: another update that is not uphill.** On the cohort of
`a_stretch_some_samples_carry_twice_is_not_read_as_heterozygosity` (30 drawn samples, 4,000 positions, 3 reads),
with the fallback off, every start's cycles start lower than the one before in 19 of 39 cycles; after pass 20 the
largest such drop per start is 1.7 to 23.2 units against median gains of 0.07 to 9.2. Holding the frequency
density's shapes fixed leaves this almost unchanged (largest drops 3.5 to 23.0). None of the six starts (two fits,
three starts each) converges in 120 passes. Before this step the same cohort's starts returned 1.5 to 16.4 units
below their best (measured in the step's review).

**The fallback, before it was dropped**, on that cohort: judged against the best point, it replayed one cycle every
three passes to the pass limit (37 times in one start); judged against the cycle before, 19 times a start, with
five of six starts returning the same log-likelihood to the printed digits as without it and one 0.4 units higher.

## Tests

- `no_shape_moves_more_than_a_factor_of_two_in_a_pass`: the kimura case above; a step inside the factor is the
  Newton step unchanged; from b = 10 a step asking for 45 stops at 20; from (1, 1) a step asking for (4, 1) stops
  at exactly 2 in the first shape; a shape that is not a number moves nothing.
- `the_shapes_step_follows_the_slope_by_the_complete_data_curvature`, unchanged in substance, now pins the step
  inside the bounds alone (`newton_step_within_the_bounds`).
- `every_start_at_the_pass_limit_returns_the_best_point_it_reached`: on the 30-sample cohort, without and with the
  duplicated class, every start at the limit returns within `JUMP_SLACK` of the best log-likelihood any of its
  passes entered with. Switching the return off fails it (start 2 returns −12,292.7 against a best of −12,267.3).
- `a_start_returns_the_best_point_it_reached` (`fit/information.rs`): the same property on three small drawn
  cohorts — a guard only; at that size no start returns below its best.
- `the_log_says_when_a_start_returned_its_best_point`.

## Deviations

- The decision said the shapes' step is "taken on the log scale and capped at a factor of 2". Built: the Newton
  step on the natural scale, kept inside a box of a factor of two either way — a cap on the log scale, the step
  itself unchanged.
- The review fanned out to two reviewers covering all nine categories between them.
