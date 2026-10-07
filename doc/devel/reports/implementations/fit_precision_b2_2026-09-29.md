# Fit precision, step B2 — a later start stops as soon as it is heading where an earlier one converged

**Date:** 2026-09-29. **Plan:** [fit_precision.md](../../implementation_plans/fit_precision.md) step B2.
**Spec:** [fit_precision.md](../../ng/spec/fit_precision.md) §3.4, with §2 and §3.3 as amended at step B1.
**Branch:** `fit-precision`. **Review:** [fit_precision_b2_2026-09-29.md](../reviews/fit_precision_b2_2026-09-29.md);
**fixes:** [fixes_applied_fit_precision_b2_2026-09-29.md](../reviews/fixes_applied_fit_precision_b2_2026-09-29.md).
This report describes the step as committed, after the review's fixes.

## 1. What was built

The SNP/indel fit runs three starting points one after another and keeps the one with the best
log-likelihood. On kimura all three reached the same answer, each after 200 passes. **A later start now
stops as soon as it is heading where an earlier start converged.**

- **What it is judged against** (`fit_jointly`, `settled::is_the_new_yardstick`). The best-scoring start
  that has converged so far: its log-likelihood, values and standard errors from its final pass
  (`EarlierAnswer::of`). A start stopped at the pass limit is never the yardstick — its errors need not
  be the maximum's (spec §3.4's second trap) — so until one converges, later starts run to their own end.
- **When it stops** (`maximise`, `settled::agrees`). On every pass that sums the information (step B1),
  the start's *projected endpoint* — each parameter's value plus its Newton distance to the maximum — is
  compared with the earlier answer. **When the endpoint lies within half of the earlier answer's standard
  error of that answer, on every parameter the earlier answer gives an error** (`agreement_fraction`,
  default `AGREEMENT_FRACTION` = 0.5), the start stops at once, without the rest of its cycle.
  - A start that the same pass finds **settled** is not stopped by agreement: it finishes its cycle,
    converges and competes.
  - A parameter the earlier answer has no error for is not compared; if it has none at all, nothing
    agrees with it. A parameter without a Newton distance is taken where it is.
- **What it returns.** The parameters that pass judged, with that pass's statistics and information. It
  runs no final pass: the judging pass already scored those parameters, and a start that agreed never
  wins, so its per-position posteriors are never read.
- **Who wins** (`settled::is_the_new_best`). The best log-likelihood among the starts that did not
  agree. The winner's statistics are asserted to come from a final pass.
- **What the log says.** "stopped, heading where start 1 converged (its projected endpoint within 0.5 of
  that start's standard error on every parameter it gives one) after 20 pass(es)", and "no last pass of
  its own, as it does not compete".
- `StartOutcome::converged` became `ended: StartEnd` — `Converged`, `AtTheLimit` or `Agreed { with }`.
  `JointFit::converged` is the winner's `ended == Converged`. The earlier answer travels in `WhichStart`
  beside the start's number and point.

## 2. Choices the spec left open

1. **A start that agreed does not compete to be the answer.** Spec §3.4 says it "contributes nothing
   new". Its stop point is part-way — 0.04 to 1.09 errors from the answer it agreed with in §4 — so
   letting it win would return a less finished fit. **What it costs is in §4.2, for the owner.**
2. **It stops mid-cycle**, right after the judging pass, and runs no final pass.
3. **A start found settled on the same pass converges rather than agreeing.** Built the other way first,
   review showed it took a start that had just met the settled test out of the competition; on the
   cross-platform fixture that alone moved the checksum.
4. **The yardstick is the best-scoring converged start so far**, not necessarily the first: if start 1
   runs to the limit and start 2 converges, start 3 is judged against start 2. Spec §3.4's decided
   sentence says "the best earlier start"; its trap speaks of "the first start" only.

## 3. What moves

- **The cross-platform checksums: nothing.** On the fixture (two samples over 600 bases) every start
  converges after 18 passes, as under B1; no start is stopped by agreement, since its judging pass at
  pass 16 already finds it settled. Both checksums pass unchanged.
- **The oracle cohort** (four tomato accessions at about 3 reads, `tmp/fit_precision/oracle_b2f/`):
  **nothing.** No start converges within 200 passes (198, 198, 200, as under B1), so none is judged against
  another; `scripts/promote_ng_oracle.sh` reports every checksum matching its baseline.
- **Drawn cohorts:** the answer moves where a later start that would have won alone now agrees (§4.2).

## 4. What was measured

A probe (`tmp/fit_precision/b2/probe2.rs`, not committed; output `probe2_keep.log`) drew 54 cohorts, nine
regimes × six seeds, and ran each twice: every start on its own (step B1's behaviour), and as `fit_jointly`
now runs them. Passes count each start's final pass.

### 4.1 Passes saved, and whether an agreeing start was really heading there

| regime | cohorts with an agreeing start | agreed at pass | passes, all six cohorts: B1 → B2 | agreeing start's own end, errors from the answer |
|---|---|---|---|---|
| 20 samples × 3 reads | 6 | 10–13 | 399 → 276 | 0.008–0.053 |
| 20 × 8, duplicated class | 6 | 10–13 | 491 → 312 | 0.003–0.052 |
| 8 × 3, duplicated class | 5 | 22–56 | 2,369 → 1,545 | 0.001–0.087 |
| 4 × 3 | 5 | 10–20 | 1,720 → 1,129 | 0.001–0.034 |
| 4 × 8, duplicated class | 3 | 10–65 | 2,714 → 2,256 | 0.024–0.087; one 0.33 |
| 2 × 3 | 0 | — | 3,590 → 3,590 | — |
| 2 × 20, duplicated class | 4 | 10–13 | 1,728 → 1,505 | 0.002–0.109 |
| 1 × 3 | 5 | 10–28 | 1,446 → 748 | 0.011–0.132 |
| 1 × 30 | 5 | 10 | 303 → 231 | 0.016–0.035 |

- **78 starts agreed, in 39 of the 54 cohorts.** Over those 39 the passes fell from 6,737 to 3,569, 47%
  fewer; over all 54, from 14,760 to 11,592. Where no start converges — every 2 × 3 cohort — nothing
  changes.
- **No start agreed while heading elsewhere.** Run on its own, every agreeing start that converges ends
  within 0.132 errors of the answer it agreed with, on the parameters compared. The one further, 0.33, is
  a start that alone runs to the pass limit (below). B1's review found that at 2 and 4 samples a Newton
  step from part-way through a fit is several times too long; no agreement here came from it. The
  correctness reviewer's 54 cohorts on other seeds (1 to 8 samples) found the same: 70 agreeing starts,
  their own ends at most 0.194 errors from the answer.
- **At one and two samples few parameters are compared**: 2 to 5 of 8 at 1 × 3, 5 of 8 at 1 × 30, 5 or 6
  of 14 at 2 × 20 with the duplicated class. The rest have no error in the earlier answer.

### 4.2 What not competing costs

In 33 of the 39 cohorts the answer is the one B1 returned or moves by at most 0.041 errors and 0.003
log-likelihood units. The cohorts where B2 returns a lower answer than B1 by more than 0.01 units:

| cohort | B1 − B2, log-likelihood | apart, errors | why |
|---|---|---|---|
| 4 × 8, duplicated, seed 4 | 0.165 | 0.33 | start 3 alone runs to the limit at a higher point than start 1's converged answer, which B1's settled test therefore stopped short of |
| 2 × 20, duplicated, seed 4 | 0.150 | 0.08 on the compared parameters | start 3's own end differs on parameters the rule does not compare — the correctness reviewer found the carrier Beta's second shape at 1.152 against the answer's 1.703 |
| 1 × 3, seed 4 | 0.124 | 0.13 | as above: 2 to 5 of 8 parameters compared |
| 1 × 3, seed 5 | 0.080 | 0.03 | as above |
| 2 × 20, duplicated, seed 5 | 0.012 | 0.09 | as above |
| 4 × 8, duplicated, seed 0 | 0.010 | 0.04 | |

**All are at one to four samples**, where the answer is least determined; at 8 and 20 samples the largest
loss is 0.0011 units. The alternative — letting an agreeing start compete — would not recover them, since
an agreeing start stops part-way (its stop point scores below the answer in 77 of 78 cases). Recovering
them means comparing the parameters without an error as well, or not stopping a start at the smallest
cohorts; both are design changes for the owner. A guard proposed in review — refuse to agree when the
judging pass already scores above the answer — would have applied to 1 of the 78 agreeing starts.

### 4.3 Determinism

The agreement decision reads the same information pass and Newton distance as step B1's settled test;
`agrees` uses only + − × ÷ and comparisons. The correctness reviewer found every start's ending, passes,
log-likelihood and parameters the same bits at 1, 4 and 8 threads on a 20-sample, 16-chunk cohort where
starts 2 and 3 agree. The committed `a_fitted_cohort_writes_the_same_bytes_at_any_pool_width` runs 1, 4 and
7 threads.

## 5. Tests and validation

| test | what it shows |
|---|---|
| `settled::tests::a_start_agrees_when_its_endpoint_is_within_half_an_error_on_every_parameter` | 0.4 and 0.25 errors agree, 0.6 does not; a parameter without an error in the earlier answer is not compared |
| `settled::tests::the_endpoint_and_not_the_value_is_compared` | a value 5 errors off with an endpoint on the answer agrees; the reverse does not |
| `settled::tests::nothing_agrees_with_an_answer_that_gives_no_error_or_at_a_nan_endpoint`, `a_distance_list_of_another_length_is_refused`, `at_a_fraction_of_zero_no_start_agrees` | the edges of `agrees` |
| `settled::tests::the_yardstick_is_the_best_converged_start`, `the_answer_is_the_best_start_that_did_not_agree` | the loop's two decisions |
| `information::tests::a_later_start_heading_where_the_first_converged_stops_there` | 4 samples × 3 reads: start 2 alone converges in 58 passes; judged against start 1's answer it agrees after 20, its value then still 1.09 errors from the answer, so only its endpoint could agree; against the answer moved 3 errors on the invariant share, it converges in its own 58 |
| `information::tests::a_start_both_settled_and_agreeing_converges` | settled and agreeing on one pass converges; never settled, it agrees at once; at an agreement fraction of zero it runs pass for pass as alone |
| `information::tests::a_start_at_the_pass_limit_is_no_yardstick` | through `fit_jointly`, on 4 samples where every start runs to the limit: the best of them is returned, pass for pass |
| `information::tests::later_starts_that_agree_leave_the_first_converged_answer` | through `fit_jointly`, on 20 samples where start 2 alone scores higher: start 1's answer is returned |
| `standard_errors::tests::a_parameter_the_newton_step_cannot_solve_has_no_distance` | carried from B1's review: a mismapped rate with no information and a share the others mimic have no distance, and every other distance is the dense solve without them, on both matrices |

- Fit module (`cargo test --release --lib parameter_estimation::joint::fit`): 106 passed, 1 ignored.
- **Each of five mutations fails a test** (`tmp/fit_precision/b2/mut/`): a start at the limit made the
  yardstick (`a_start_at_the_pass_limit_is_no_yardstick`, and the predicate's test); the value compared for
  the endpoint (`a_later_start_heading_where_the_first_converged_stops_there`); agreement before the settled
  test (`a_start_both_settled_and_agreeing_converges`); agreed starts competing, and the first converged answer
  kept in place of the best (the predicates' tests).
- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 4,979 passed, 3 failed (the
  pre-existing `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 5 ignored.

## 6. For the owner, and what is left

1. **Whether the losses of §4.2 are acceptable** — up to 0.165 log-likelihood units at one to four
   samples, none above 0.0011 at 8 and 20. The recommendation is to accept them for checkpoint B: kimura
   (2,169 samples) is where the passes matter, and there every parameter has an error.
2. **Left for step B3:** the per-start record in `JointFit`; `StartEnd` becomes public there. Deferred
   from review: a start-number type in place of a bare `usize`, and one value for the judging pass's
   verdict in `maximise` (fixes report).
