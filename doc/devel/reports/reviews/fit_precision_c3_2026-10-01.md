# Code Review: fit_precision_c3
**Date:** 2026-10-01
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step C3 — the repeat-tract climb stops when a Newton step puts every number within a tenth of its
error of the peak
**Status:** Request-changes (applied: [fixes_applied_fit_precision_c3_2026-10-01.md](fixes_applied_fit_precision_c3_2026-10-01.md))

---

## 1. Scope

- **What was reviewed:** the diff `91995ed5..f7aa3174` (a review object of the working tree): in `ssr_fit.rs`
  the walk's stopping rule (`walk_until_settled`, `StratumClimb::judge`, `newton_distances_on_the_natural_scale`,
  `furthest_in_errors`, `HeldPoint`, the reuse of the winning walk's errors), `StratumFit::ending` and `walks` and
  their readers, `DEFAULT_MAX_ROUNDS`, `SsrFitConfig::settled_fraction`; the spec §4.3 amendment; the oracle
  baseline's two fit lines; the implementation report and the stopped report.
- **Categories:** three reviewers, each in its own worktree detached at `f7aa3174` — *correctness and cost*;
  *design and claims*; *tests and errors*. Findings:
  `tmp/review_2026-10-01_fit_precision_c3/{correctness_cost,design_claims,reliability_errors}.md`, evidence
  beside them.

## 2. Verdict

Request-changes. The rule, the Newton distance and the reuse of the winning walk's errors are correct and give
the same bits on both schedules and at any pool width. But half the judgement — the shares' and the
concentration's distances — could be wrong with every test passing (M1). The round limit of 30 was justified by
a range a drawn walk contradicts (M2). The spec and the report describe behaviour the code does not have (M3,
M4).

## 3. Execution status

- Module tests on the reviewed tree: 51 passed, 3 ignored (all three reviewers).
- Mutations: 22 planned and 8 run (correctness), 13 run (design), 28 run (tests and errors), each restored and
  the tree checked clean. Against the author's suite, the tests-and-errors reviewer's 28 left 14 survivors, 5 of
  which cannot change behaviour on any input a fixture makes.
- Experiments: the returned score against a fresh score, and the reused errors against fresh ones, bit for bit
  on 24 walks and both schedules; Newton distances against known moves of the level, a share and the
  concentration; two thirteen-class draws refitted with every walk's record printed; the oracle calls before and
  after compared field by field.
- Every number in the reports, the spec amendment, the baseline comment and the doc comments was checked by at
  least one reviewer; the wrong ones are findings below.

## 4. Open questions and assumptions

1. Whether a judgement where no number has an error should count as settled (Mi1) — the spec says yes.

## 5. Top 3 priorities

1. **M1** — test the shares' and the concentration's distances, and the kept-coordinate indexing.
2. **M2, M3** — re-measure the round limit, and make the spec say what the code does.
3. **M4** — correct the oracle comparison the checkpoint is decided on.

## 6. Findings

### Major

**M1: the Newton distance's carriage to the shares and the concentration, and its indexing, are untested.**
**Categories:** correctness and cost; tests and errors; design and claims. **Confidence:** High.
Dropping the shares' cross terms (the largest class's distance then reads 0, against 0.217 measured), leaving the
concentration on the log scale (about 2× off), leaving the shares on the log-ratio scale, using the inverse's
diagonal only, leaving the shares out of `furthest_in_errors`, and indexing the slope by column instead of by
kept coordinate all passed the 51 tests. `a_point_one_error_from_the_peak_has_not_settled` moves the level only.

**M2: the round limit's reason is contradicted by a drawn walk.** **Categories:** design and claims; correctness
and cost. **Confidence:** High. The doc, the spec and the report said walks that settled took 3 to 29 rounds at
thirteen classes and that none settled after 29. In the thirteen-class draw 2, a walk settled at round 36. And the
29 was lengthened by the put-off itself: the same walk judged at every trigger settled at 17.

**M3: the §4.3 amendment disagrees with the code.** **Categories:** design and claims; correctness and cost.
**Confidence:** High on the mismatch. A losing round inside the put-off was not judged, and could end the walk as
having lost a round at a point never judged; the spec says every losing round is judged. The winning walk's
errors are reused whenever its last round was judged at the returned point; the spec and the doc said "when it
settled". The reviewer also asked whether the put-off and the stop on a losing round were the owner's choices.

**M4: report §2.2's oracle comparison is miscounted, and two mechanisms are wrong.** **Categories:** design and
claims. **Confidence:** High. Five strata, not four, differ on some number by 0.37 to 233 errors at 0.000 to 0.054
units short. Stratum 2:9's differing share is 0.031, not near zero. In 1:15 no walk reached 40 rounds. 1:8's
"equals the reference" holds by construction, since its winner is the reference climb.

**M5: nothing tests that the walk's trigger is in total units at the call that matters.** **Categories:** tests and
errors. **Confidence:** High. Passing one tract instead of the stratum's count in `climb_from` passed every test.

**M6: the score `HeldPoint` holds is never checked against the point it holds.** **Categories:** tests and errors.
**Confidence:** Medium. Returning the right bracket's score with the left bracket's point passed the suite, and
changes where a losing round lands.

**M7: the ignored comparison asserts nothing.** **Categories:** tests and errors. **Confidence:** High. The step's
acceptance measurement prints and passes whatever the rule does.

### Minor

- **Mi1:** a judgement where no number has an error reports `Settled`, and the record cannot tell it from a real
  one (*correctness and cost; tests and errors*). On the oracle 130 of the winners' 255 numbers carry no error.
- **Mi2:** a walk whose start scores `NaN` never stops as having lost a round; it runs every round.
- **Mi3:** the baseline comment and report §2.4 say AF, GQ and QUAL moved "by at most 0.5"; sample GQs moved by up
  to 2, AF by up to 0.036.
- **Mi4:** a round at thirteen classes is 307 evaluations, not 288, so a judgement is about 1.7 rounds, not 1.8
  (*design and claims; correctness and cost; tests and errors*).
- **Mi5:** report §4 item 3, "53 judgements where judging at every trigger took 104, with the same answers": the
  put-off cost 22 more rounds, and draw 0's answer moved (0.040 errors against 0.083).
- **Mi6:** report §2.1, "in the fifth draw no walk settled": two of its three walks settled (*design and claims;
  tests and errors*).
- **Mi7:** `DEFAULT_MAX_ROUNDS`'s doc blames the old limit of 5 for a 17.6-unit shortfall that the new rule at 40
  rounds still has (17.4).
- **Mi8:** stale text the step should have updated: the plan's C3 entry and checkpoint C, spec §4.1, §4.2's
  "578 at p = 17", §4.5 item 2's twenty rounds, §6's time, and `DEFAULT_REFUSAL_FLOOR`'s "83% settle", measured
  under the old rule.
- **Mi9:** `StratumFit::ending` and `walks` can disagree, and hand-built fixtures carry a settled fit with no walks.
- **Mi10:** `the_best_walk` and `WalkedStratum` hold the same three things in two orders, as tuples.
- **Mi11:** `Walked` repeats `Climb` field for field.
- **Mi12:** `StratumDistances` mirrors `StratumErrors` by position, paired by `zip`.
- **Mi13:** "climb" and "walk" name the same one-start unit in different types.
- **Mi14:** `PeriodLengthSpectrum::converged` stays a `bool` named for the old rule.
- **Mi15:** `SsrFitConfig::settled_fraction`'s default (0.1) is not stated where a caller reads it.
- **Mi16:** an empty `starting_points` panics without naming the field (code that predates the step).
- **Mi17:** `climb_endings_summary` has no test.
- **Mi18:** no test sets `settled_fraction`, so a judge reading the constant instead passes.
- **Mi19:** the two-schedule tests compare no ending or walk record, and run at 2 rounds, where the errors are
  rarely reused (*tests and errors; correctness and cost*).
- **Mi20:** no fixture has two slippage groups, so a judge reading every group's distance from group 0 passes.
- **Mi21:** after a round that loses with a better point part-way through, the next projection uses that partial
  gain; neither the spec nor a test says which gain it should use.

### Nits

`settled()` should be `is_settled()`; `WalkRecord`'s `ending` and `rounds` undocumented; "one a starting point";
`LostARound`'s first sentence carries three conditions; `remaining_gain_target` cites "plan step C3's report";
the spec's "not placed" undefined; `one_round` re-scores the point its last golden section scored (one evaluation
in 307); the comparison test's knobs (a reference of 0 rounds panics, a class count other than 3 or 13 runs
nothing); a test doc's "the default five"; `HeldPoint::offer` on exact ties.

## 7. Out of scope observations

- The 256-point average the likelihood stands on — measured separately
  ([fit_precision_quadrature_average_2026-10-01.md](../implementations/fit_precision_quadrature_average_2026-10-01.md)).
- Strata where another start's walk climbs higher than the settled winner (two on the oracle cohort): the
  starting points, not the stopping rule.

## 8. Missing tests to add now

`newton_distances_are_carried_to_each_numbers_scale_like_the_errors` (M1); the furthest number over every kind and
a `NaN` distance (M1); `the_held_point_scores_what_it_holds` (M6); assertions in the comparison (M7); the settled
fraction read from the config (Mi18); the climb summary (Mi17); the winning walk's errors against a fresh
curvature in both schedules at the default limit (Mi19); a two-group stratum (Mi20); a `NaN` start (Mi2).

## 9. What's good

- The judgement takes its slope from the curvature's own evaluations, so the Newton step costs nothing beyond
  the curvature C1 already computes.
- The `Climbing` trait lets the stopping rule be tested on scripted gains, without a likelihood.
- Reusing the winning walk's errors gives the same bits as a fresh curvature, on both schedules — checked by
  experiment and by a mutation that perturbs them by 10⁻⁹.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib parameter_estimation::joint::ssr_fit`;
`NG_FIT_PRECISION_CLIMB_DRAWS=2 NG_FIT_PRECISION_CLIMB_CLASSES=3 scripts/dev.sh cargo test --release --lib
the_climbs_stop_loses_nothing_against_a_longer_climb -- --ignored --nocapture`.
