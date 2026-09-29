# Code Review: fit_precision_b2
**Date:** 2026-09-29
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step B2 — a later start of the SNP/indel fit stops as soon as it is heading where an earlier one converged
**Status:** Request-changes (applied: [fixes_applied_fit_precision_b2_2026-09-29.md](fixes_applied_fit_precision_b2_2026-09-29.md))

---

## 1. Scope

- **What was reviewed:** the diff `7e36285f..d9f2adab` (a review object, not on the branch; `ca11a2d8` between
  them only re-records the oracle baseline): `fit.rs` (`fit_jointly`'s start loop, `WhichStart`, `maximise`,
  `StartOutcome`, `JointFitConfig::agreement_fraction`), `fit/settled.rs` (`AGREEMENT_FRACTION`, `StartEnd`,
  `EarlierAnswer`, `agrees`), the new tests in `fit/information.rs` and `fit/standard_errors.rs`,
  `cli/cross_platform_digests.rs` (the fit checksum re-recorded), and the implementation report.
- **Out of scope:** the repeat-tract fit; B1's rule itself, except where B2 reads it.
- **Categories:** three reviewers, each in its own worktree —
  *correctness and cost* (reliability of behaviour, float_portability, unsafe_concurrency, extras);
  *design and claims* (naming, idiomatic, defaults, module_structure, smells, refactor_safety, diff against
  intent, every number in the prose); *tests and errors* (reliability: tests, mutation testing; errors).
  Findings: `tmp/review_2026-09-29_fit_precision_b2/{correctness_cost,design_claims,reliability_errors}.md`,
  evidence beside them.

## 2. Verdict

Request-changes: neither of spec §3.4's two traps has a test that fails when it is broken (B1), and the order
of the agreement and settled checks takes a start that has just converged out of the competition (M1).

## 3. Execution status

- The reviewers ran the fit module and the cross-platform tests (101 passed, 1 ignored), 31 mutations between
  them (10, 7 and 14), and probes on 54 new drawn cohorts (1 to 8 samples, 3 to 30 reads), one-sample cohorts,
  B1's build on the cross-platform fixture with a per-pass trace, and determinism at 1, 4 and 8 threads.
- The oracle, the full suite, fmt and clippy were checked against the author's logs, not re-run.
- Needs verification: none open.

## 4. Open questions and assumptions

1. Whether an agreeing start should compete (the report's choice 1) — for the owner at checkpoint B, with
   M3's measurement.

## 5. Top 3 priorities

1. **B1** — tests that fail when a start at the pass limit becomes the yardstick, and when the value is
   compared instead of the projected endpoint.
2. **M1** — check the settled test before agreement, so a start that has converged competes.
3. **M2/M3** — correct the report's account of what not competing costs, with the cases measured.

## 6. Findings

### Blocker

**B1: `fit.rs:1883`, `fit.rs:2210` — Neither of spec §3.4's two traps has a test that fails when it is broken.**
**Categories:** all three (convergent). **Confidence:** High.
Letting a start stopped at the pass limit become the yardstick passes every test; on cohorts where no start
converges, later starts then "agree" with it after 19 to 70 passes and the returned answer moves (2 samples × 3
reads, seed 0: −3582.5268 against −3582.4637; 4 × 3: −6423.406 against −6423.398; another 4 × 3 cohort
−5987.3107 against −5987.2868), and the log says "heading where start 1 converged" of a start that never did.
Comparing the value instead of the projected endpoint also passes: the integration test's 20-sample cohort agrees
at pass 10 either way; at 4 × 3 agreement comes 12 to 24 passes later. Fix: a `fit_jointly`-level test on a
cohort where every start runs to the limit; an integration cohort whose stop point's value is still over half an
error from the answer.

### Major

**M1: `fit.rs:2219` — Agreement was checked before the settled test, so a start that had just converged was
recorded as agreed and lost the win; that alone moved the cross-platform checksum.**
**Categories:** design and claims; tests and errors; correctness and cost (case analysis). **Confidence:** High.
On the fixture, starts 2 and 3 agree on pass 16, the same pass that settles them under B1; with the settled test
first the fixture writes B1's bytes exactly. The report presented the order as a tie-break. Fix: agree only
when some parameter is not settled.

**M2: report §3, §6 — The explanation of why the agreeing starts stopped above the answer is wrong.**
**Categories:** design and claims; tests and errors; correctness and cost. **Confidence:** High.
At pass 16 all three starts lie within 0.015 units of one another (−44.8858, −44.8737, −44.8710); each then loses
about 0.07 over its last cycle, the climb not being monotone here (`JUMP_SLACK` = 1.0). The report compared
starts 2 and 3 at pass 16 with start 1 after its final pass, and put choice 1 to the owner on that comparison.

**M3: report §4.1 — What not letting an agreeing start compete costs is larger than reported, through a
mechanism the report does not describe, and the 4 × 8 case is presented as within the rule's allowance.**
**Categories:** correctness and cost; design and claims. **Confidence:** High.
Only parameters the earlier answer gives an error are compared; at 2 samples × 20 reads with the duplicated
class 6 of 14 are. On one such cohort start 3 agrees at 0.076 errors on the compared parameters while its own
end has the carrier Beta's second shape at 1.152 against 1.703, and scores 0.150 units above what B2 returns.
The 4 × 8 case (0.165 units, 0.33 errors) is B2 returning a lower answer: the half error bounds the endpoint at
the stop, not where the start would have ended. Fix: report both, with the losses measured, for the owner;
say "on every parameter the earlier answer gives an error" in the doc and the log.

**M4: `fit.rs:2304`, `:2346-2360` — An agreed start's statistics carry no per-position posteriors, and only
`competes` keeps such a start from winning.**
**Categories:** design and claims; correctness and cost; tests and errors. **Confidence:** Medium.
Flipping choice 1 would hand empty posteriors to the contamination fit. The `statistics` doc still said "the
final pass's … the per-position posteriors". Fix: correct the doc; assert the winner's statistics come from a
final pass.

**M5: tests — Nothing exercises `fit_jointly` choosing the yardstick and the winner.**
**Categories:** tests and errors. **Confidence:** High.
The only test running `fit_jointly` with agreement is the cross-platform checksum. "The first converged answer
instead of the best" and "agreed starts compete" need two pure predicates to test cheaply. Fix: the predicates,
their unit tests, and two `fit_jointly`-level tests.

### Minor

**Mi1: `settled.rs:112-131` — `agrees` asserted one of three length pairs**; a short `distances` would make a
start agree more easily. *(all three)*

**Mi2: `agrees` is true when the earlier answer gives no parameter an error**, at any fraction, so "at zero no
start agrees" was false in that case. *(all three)*

**Mi3: report §2 choice 4 — the yardstick being the best converged start departs from the second trap's literal
"the first start"**, while the report said it followed it. *(design and claims)*

**Mi4: report §4.1 — the passes column left out the cohorts in which no start converged** (4 × 3: 1,443 → 1,013
over all five, not 848 → 418). *(design and claims; correctness and cost)*

**Mi5: "in the fifth significant figure" is wrong** (the 7th and 6th). *(design and claims)*

**Mi6: the earlier answer built in two places, its values in three.** *(design and claims)*

**Mi7: `ended`, `agreed_with` and `agreed_at` kept consistent by hand.** *(design and claims)*

**Mi8: a start number counted from one travels as a bare `usize` beside a zero-based `index`.** *(design and
claims)*

**Mi9: `settled.rs`'s module doc did not describe the agreement rule.** *(design and claims)*

**Mi10: nothing measured at one sample**, the low end of the range. *(design and claims)*

**Mi11: `best.expect`'s message did not name the invariant** (the first start never agrees). *(tests and
errors)*

**Mi12: report §4.2 — "1, 4 and 8 threads"**: the committed test runs 1, 4 and 7. *(design and claims)*

### Nits

`StartEnd` could be `StartEnding`; `Agreed { with }` could be `with_start`; `let slot = 1` where
`cohort::P_INVARIANT` names it; `pub` fields on a `pub(super)` struct; the log closure's `to_owned` and
`Option<String>`; no test gives `agrees` a NaN; the oracle baseline's re-record belongs to its own commit.

## 7. Out of scope observations

- `fit.rs` `JUMP_SLACK` = 1.0: on the two-sample fixture every start loses 0.07 to 0.15 log-likelihood units within
  one cycle (pass 17 → 18: −44.8626 → −45.0164), so a converged answer can score below a point the same start
  passed through. A question about the accelerator, not B2's.

## 8. Missing tests to add now

`a_start_at_the_pass_limit_is_no_yardstick`, `later_starts_that_agree_leave_the_first_converged_answer`,
`a_start_both_settled_and_agreeing_converges`, an endpoint-only agreement in
`a_later_start_heading_where_the_first_converged_stops_there`, `the_yardstick_is_the_best_converged_start`,
`the_answer_is_the_best_start_that_did_not_agree`, `nothing_agrees_with_an_answer_that_gives_no_error`. Code for
most is in `reliability_errors.md`.

## 9. What's good

- No false agreement in 70 + 78 agreeing starts over 108 drawn cohorts: an over-long Newton step at small cohorts
  overshoots, so a start disagrees rather than agrees.
- The agreement decision is the same bits at 1, 4 and 8 threads with agreeing starts.
- The new `a_parameter_the_newton_step_cannot_solve_has_no_distance` kills a stale-index mutant in
  `whole_matrix_solve` that no older test caught.
- An agreed start reusing its judging pass saves a final pass at no cost to any output.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`;
`scripts/dev.sh cargo test --release --lib cli::cross_platform_digests`;
the probe, `tmp/fit_precision/b2/probe2.rs` pasted into the tests module, `NG_PROBE_SEEDS=6`.
