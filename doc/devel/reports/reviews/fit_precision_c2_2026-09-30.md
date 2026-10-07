# Code Review: fit_precision_c2
**Date:** 2026-09-30
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step C2 — a measurement test of whether a stratum's errors mean what they say
**Status:** Request-changes (applied: [fixes_applied_fit_precision_c2_2026-09-30.md](fixes_applied_fit_precision_c2_2026-09-30.md))

---

## 1. Scope

- **What was reviewed:** the diff `33072de2..c5010ff2` (a review object): the ignored test
  `a_stratums_errors_mean_what_they_say` and its helpers in `ssr_fit.rs`, and the implementation report, above all
  its §3 diagnosis of the thirteen-class results.
- **Categories:** three reviewers, each in its own worktree — *correctness of the diagnosis*; *design and claims*;
  *tests and errors*. Findings: `tmp/review_2026-09-30_fit_precision_c2/{correctness_cost,design_claims,reliability_errors}.md`,
  evidence beside them.

## 2. Verdict

Request-changes: the test fails on unchanged code at the draw count its bounds claim to hold at (B1) and passes
when it measures nothing (M1); the report's §3 is right about the cause but overstates what more points achieve
and leaves the owner no recommendation (M3, M4).

## 3. Execution status

- A harness refitting the test's own draws reproduced all twelve three-class lines to the last digit; 12 code
  mutants and 13 exact scalings of every error; the test at 50 to 100 draws; the log-likelihood at the truth
  with the integral at 256 to 65,536 points; five-class draws; three-class draws refitted at 4,096 points.

## 4. Open questions and assumptions

1. The integral (M3): for the owner.

## 5. Top 3 priorities

1. **B1** — hold the bounds only at the default count.
2. **M1, M2** — the test must fail when nothing is measured, and check that every number has an error.
3. **M3, M4** — correct §3 and give the owner a recommendation.

## 6. Findings

### Blocker

**B1: `hold_to_…` bounds at ≥ 50 draws — the test fails on unchanged code at 50 draws.**
**Categories:** tests and errors; correctness. **Confidence:** High.
"3 classes, 300 tracts x 20 samples x 30 reads, shorter share: the estimates sit +0.358 errors from the truth"
against a bound of 0.35; 19 of the counts from 50 to 100 fail. The shorter share sits +0.22 errors high in every
run, and the mean's own noise is 0.15 at 50 draws.

### Major

**M1: the test passes when every stratum is refused.** *(tests and errors; design and claims)* Both `continue`
paths skip a draw uncounted, and only the kinds present are checked: `fit_stratum` forced to `None` passed in
0.15 s.

**M2: the bounds cannot see errors 10 to 20% off, and never check numbers without an error.** *(tests and errors;
correctness)* Every error scaled by 0.88 to 1.23 passes at 3 reads; doubling the off-diagonal curvature survived
while turning 6 to 12 concentration errors and 26 to 27 share errors into "not identified".

**M3: report §3 — right cause, overstated remedy.** *(correctness; design and claims)* The draw is exactly the
fit's model, so the 256-point integral is the only approximation between them, and at the truth it is 3.7 units a
tract low. But 4,096 points is still 0.45 a tract low; the 1,024-point climb did not converge; the defect starts by
five classes; it reaches the three-class shares' errors (×1.29 to ×1.56 at 4,096 points). §3 did not address the
slippage level, the number C1 left for C2.

**M4: §3 gives the owner a decision without a recommendation, and its framing is half right.** *(design and
claims)* The point count is a setting, not the model; the real reason it is the owner's is that raising it moves
every stratum and reverses a recorded decision — and C1's errors, C3's check, D's target and E2 all stand on it.

### Minor

- The shares' bounds sit under their measured figures and so accept the known under-coverage.
- A kind's line printed after its assertion hides the rest of a failing run.
- `kind.contains("length")` picks the bounds; magic numbers (0.4, 20, the bounds) inline; an `allow` without its
  reason.
- A non-Unicode `NG_FIT_PRECISION_STRATUM_DRAWS` read as unset.

### Nits

`group` holds errors; `StratumCoverage` duplicates `Tally` in `fit/information.rs` with a different spread; the
three-class regimes share seeds; 13 classes at 30 reads not run.

## 7. Out of scope observations

- The 256-point integral (M3) — the fit's model, for the owner.

## 8. Missing tests to add now

Skipped draws counted and required to be none; every number with an error at three classes; the spread bounded;
a five-class regime and a check of the integral at the truth against 65,536 points (both for after the owner's
decision).

## 9. What's good

- The coverage arithmetic — distances, shares within one and two, spread, the reference class — is right, and
  the report's §2 table matches the logs exactly.
- Seeding makes a 40-minute measurement reproducible to the last digit.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib a_stratums_errors_mean_what_they_say -- --ignored --nocapture`.
