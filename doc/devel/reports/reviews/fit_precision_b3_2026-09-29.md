# Code Review: fit_precision_b3
**Date:** 2026-09-29
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step B3 — each start of the SNP/indel fit records how it ended, and why when it did not converge
**Status:** Approve-with-changes (applied: [fixes_applied_fit_precision_b3_2026-09-29.md](fixes_applied_fit_precision_b3_2026-09-29.md))

---

## 1. Scope

- **What was reviewed:** the diff `0b9459a0..fc0a89a0` (a review object, not on the branch): `fit.rs`
  (`StartRecord`, `JointFit::starts`, `how_far_from_settled`, the start's log line, `PassSummary::not_settled`),
  `fit/settled.rs` (`StartEnding`, `FurthestFromSettled`), the tests in `fit/information.rs`, and the
  implementation report.
- **Categories:** three reviewers, each in its own worktree — *correctness and cost*; *design and claims*;
  *tests and errors*. Findings: `tmp/review_2026-09-29_fit_precision_b3/{correctness_cost,design_claims,reliability_errors}.md`,
  evidence beside them.

## 2. Verdict

Approve-with-changes: the record is right in every case built, no output byte moves, and the reason is computed
at the values each start returned; but the log names parameters in a way the owner cannot read on kimura (M1),
and the tests check neither the record's content nor where the per-pass verdict lands (M2, M3).

## 3. Execution status

- The reviewers ran the fit module and the cross-platform tests (110 passed, 1 ignored), 37 mutations between
  them (4, 7 and 26), one-sample and duplicated-class cohorts, determinism at 1, 4 and 8 threads on six fixtures,
  and the reason's cost on a drawn cohort of kimura's size (2,169 samples, 2,653 read groups).
- Needs verification: none open.

## 4. Open questions and assumptions

None.

## 5. Top 3 priorities

1. **M1** — name the furthest parameter as the reader knows it: samples by name, words not trace keys.
2. **M2** — a loop-level fixture whose judged passes find exactly one unsettled, and assertions on which pass
   carries each verdict.
3. **M3** — check the recorded reason against an independent recomputation at the returned values.

## 6. Findings

### Major

**M1: `fit.rs:1907-1918`, `:2431` — The start's log line names parameters by the trace's keys and positions.**
**Categories:** design and claims; correctness and cost; tests and errors (convergent). **Confidence:** High.
`density_a`, `duplicated_share`, `clean_error_3`, `hom_excess_1734`: at kimura's 2,169 samples nothing maps
`hom_excess_1922` (a drawn cohort of that size named it) back to a sample, while the standard-error line in the
same log already says "allele-frequency shape a" and `fit_jointly` holds the sample names. "0.39 errors" beside a
sequencing-error rate reads as a count of errors. Fix: one plain-name function; "standard errors".

**M2: `information.rs:4118` — The loop-level test cannot fail for "stop with one unsettled" or for the verdict
written on the wrong pass.** **Categories:** all three. **Confidence:** High.
Its fixture judges twice (3, then 0), so stopping at ≤ 1 or ≤ 2 gives the same result; stopping at ≤ 1 survives
the whole suite while moving fits (21 passes for 24; 115 for 121). Writing each verdict on the first trace row
also survives. Fix: 20 samples at seed `+3` (3, then 1, then 0), with assertions on the passes carrying verdicts.

**M3: `information.rs:4054` — `each_start_is_recorded_as_it_ended` checks the reason exists and agrees with
itself, not that it is right.** **Categories:** all three. **Confidence:** High.
Surviving mutations, each changing the record: fraction 0 ("17 of 17"), total + 5, the first parameter's name,
the distance doubled, the reason from the last judged cycle (0.23497 for 0.23444 errors). Fix: recompute from a
fresh pass at the returned parameters at a non-default fraction.

### Minor

**Mi1: `fit.rs:2431`, `:2334` — `names[at]` without asserting the name list and the error list agree in length.**
*(tests and errors)*

**Mi2: a start that did not converge and has no parameter to judge prints nothing in place of the reason.**
*(design and claims; tests and errors)*

**Mi3: `fit.rs:1911` — a start that reaches the limit before any cycle judged it, already settled, would print
"0 of N parameter(s) not yet within 0.1 … at 0.00 errors" after "stopped at the pass limit".** *(correctness)*

**Mi4: the verdict is attached with `trace.last_mut()`**, right only while `step` adds one trace row a call.
*(design and claims)*

**Mi5: `FurthestFromSettled::parameter` is a `String` packing kind and index**, which step E2 would have to
parse. *(design and claims)*

**Mi6: an agreed start repeats the inversion and Newton solve its judging pass just did; a winner at the limit
inverts a third time.** *(design and claims)*

**Mi7: the report's "1.6 ms" is the inversion alone**: the whole reason costs 5.45 ms a start at kimura's size,
against 524 ms for one pass over 2,000 positions there. *(correctness and cost)*

**Mi8: the oracle's start line runs to about 450 characters.** *(design and claims)*

### Nits

`named()` builds every name to use one; public types in the private `settled` module, re-exported (precedent
exists; rustdoc is quiet); a NaN distance would print "NaN errors"; the report's header links point at files
that did not exist yet.

## 7. Out of scope observations

- One pass's information sums differ in their last digits depending on whether it keeps per-position posteriors
  — deterministic, and older than B3; an agreed start's reason and an at-limit start's come from sums rounded
  differently. *(correctness and cost)*

## 8. Missing tests to add now

`a_fit_stops_only_on_a_pass_that_found_every_parameter_settled` on the seed-`+3` fixture;
`a_start_at_the_limit_reports_the_settled_test_at_the_values_it_returned`; a unit test of the log clause. Code
in `reliability_errors.md`.

## 9. What's good

- The reason is computed at the values each start returned: bit-identical to a fresh pass there, and different
  from the last judged cycle's (correctness probe).
- `fit.starts` and the trace verdicts are the same at 1, 4 and 8 threads on six fixtures.
- The whole reason costs 5.45 ms a start at kimura's size.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib -- parameter_estimation::joint::fit cli::cross_platform_digests`.
