# Code Review: fit_precision_c1
**Date:** 2026-09-30
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step C1 — each repeat-tract stratum's fit carries the standard errors of its numbers
**Status:** Request-changes (applied: [fixes_applied_fit_precision_c1_2026-09-30.md](fixes_applied_fit_precision_c1_2026-09-30.md))

---

## 1. Scope

- **What was reviewed:** the diff `efc325e2..57618353` (a review object, not on the branch): `ssr_fit.rs` (the
  types `StratumErrors` and `SlippageErrors`, `StratumFit::standard_errors`, `the_best_walk`, `curvature_of`,
  `CurvatureLayout`, `standard_errors_at`, `errors_on_the_natural_scale`, `standard_errors_summary`, both thread
  schedules, the tests), `invert_identified`'s visibility, the test fixtures elsewhere, and the implementation
  report.
- **Categories:** three reviewers, each in its own worktree — *correctness and cost*; *design and claims*;
  *tests and errors*. Findings: `tmp/review_2026-09-30_fit_precision_c1/{correctness_cost,design_claims,reliability_errors}.md`,
  evidence beside them.

## 2. Verdict

Request-changes. The arithmetic is right: rebuilt independently (plain coordinates, own differences, own
inverse), all seven errors of a drawn stratum matched to 0.3%, and over 20 redrawn strata the spread of the
fitted numbers was 0.83 to 1.22 times the reported error. But a number the climb runs off to the end of its reach
gets an error hundreds of times too small (M1), a runaway concentration gets rounding noise (M2), and the tests do
not guard the parallel schedule's errors (B1).

## 3. Execution status

- 47 mutations between the reviewers (10, 6 and 31); independent rebuilds of the curvature; 20 redrawn strata;
  thin and runaway strata; costs at 13 classes; both schedules at 1, 4 and 7 strata at once and pool widths 1, 4
  and 8.
- The full suite, fmt and clippy were checked against the author's logs, not re-run.

## 4. Open questions and assumptions

1. What an error of a number the tracts do not place should be (M1, M2, M3): recommended by all three reviewers,
   no error, by a reason, as spec §3.2's "wider than its range" does for the SNP/indel fit.

## 5. Top 3 priorities

1. **M1/M2** — no error for a number whose error on the climb's own scale says the tracts do not place it.
2. **B1** — compare the errors in the schedule tests.
3. **M3/M4** — say why an error is absent, and correct the report's cost claim.

## 6. Findings

### Blocker

**B1: `ssr_fit.rs` `every_fitted_number` — the schedule tests never read the errors, and the report says they do.**
**Categories:** tests and errors; correctness and cost. **Confidence:** High.
Returning no errors in the several-strata-at-once branch, computing them at another step there, from the wrong
stratum's tracts, or at a losing walk's answer, all pass every test. On the code as it stands the errors are the
same bits in every schedule (both reviewers' probes).

### Major

**M1: `errors_on_the_natural_scale` — a slippage number the climb left at the end of its reach gets an error
hundreds of times too small.** **Categories:** correctness and cost; tests and errors. **Confidence:** High.
A golden section moves a number at most a fixed span a round, so a number whose likelihood keeps rising towards 0
or 1 stops where five rounds took it (the fall-off at exactly 1.3178 × 10⁻⁷ on 3 of 5 thin drawn strata). There
the likelihood is flat, and `p(1 − p)` turns a logit error of about 1,000 into one of 10⁻⁴: a shorter share of
0.9999997 ± 0.0003, where moving it to the true 0.83 costs 0.646 log-likelihood units. At one sample and one read
a tract the level came back 0.059 ± 5.20. Fix: no error where the error on the climb's scale exceeds about 3
(numbers the tracts determine sat at 0.02 to 0.5, stuck ones at 1,000 to 2,400).

**M2: a concentration run off to the climb's reach gets a curvature that is rounding noise.**
**Categories:** correctness and cost; design and claims. **Confidence:** High.
Drawn at 50, fitted at 8.0 × 10⁶: an error of 3.7 × 10⁹ at the default step, none at a third of it, 3.4 × 10⁹ at
three times it. The oracle's "largest 247.649" times the concentration has the same signature. Fix: the same rule
on the log scale.

**M3: four reasons for a missing error are one `None`, and the log calls all of them "not curved downwards".**
**Categories:** design and claims. **Confidence:** High.
A number dropped as not identified, a class held at zero share, a variance not positive, and `p(1 − p)` = 0.

**M4: report §3.3 — "+21% on the oracle" is run-to-run noise, and its mechanism does not hold.**
**Categories:** design and claims; correctness and cost. **Confidence:** High.
The stratum fit's code was unchanged before C1, and four earlier oracle runs of it took 1 min 24 s to 1 min 42 s
(B3's own first run 1 min 42 s, the same as C1's). The climb rebuilds the quadrature on 3,780 of its 4,608
evaluations too, so the rebuilds do not favour the errors on thin strata: on the thinnest stratum drawn the errors
cost 12.3% of the climb.

**M5: after `fit_strata` the errors sit beside blended numbers they are not the errors of.**
**Categories:** design and claims. **Confidence:** Medium.
The smoothing overwrites `StratumFit::slippage` with a blend; the errors stay the own fit's. Spec §5.2 gives a
blended number no error.

**M6: nothing checks that the curvature is taken at the fitted answer, that tracts without reads leave it alone,
or what the summary prints.** **Categories:** tests and errors; correctness and cost. **Confidence:** High.
Mis-centring the spectrum moved the errors up to 45%; removing the spectrum's renormalisation 10 to 34%; the
tract count taken over tracts with reads only; a wrong scale or count in the summary — all survive.

### Minor

**Mi1: `CURVATURE_STEP`'s doc — "a ten-thousandth" holds only near zero on the climb's scale**, and scaling by
`max(1, |x|)` has no reason on logit and log scales. *(design and claims)*

**Mi2: the 10⁻⁸ identification floor was calibrated on exact scores; here the curvature is a finite
difference.** *(design and claims, medium)*

**Mi3: the drop order** — when two numbers cannot be told apart, the later is dropped (slippage first, then the
classes, then the concentration), and the kept one's error is taken with the other held; not recorded in the
report. *(design and claims)*

**Mi4: the log line** mixes relative and absolute errors in one format, never summarises the class shares, and
prints 247.649. *(design and claims)*

**Mi5: `invert_identified` shared by widening a private submodule of the SNP/indel fit**; the ~300 new lines could
sit in their own file. *(design and claims)*

**Mi6: the field doc names a path that does not exist** (a fit read back from a parameters file). *(design and
claims)*

**Mi7: two `.expect` calls without a PANIC-FREE note; an empty `starting_points` panics.** *(tests and errors)*

**Mi8: the several-at-once schedule starts the errors after every walk**, lengthening its tail by about 11% of the
longest climb. *(correctness and cost)*

### Nits

`Result` as an either-type; `which: usize` with a `_ =>` fall-through; `fit_pooled` no longer pools; spec §4.2's
p = 17 and 578 evaluations should be 16 and 513 (the code is right).

## 7. Out of scope observations

- The fitted level sat 3.0 to 4.5 errors above the truth on three drawn 13-class strata (author's measurement) —
  for step C2.

## 8. Missing tests to add now

`every_fitted_number` with the errors; a number left at the climb's reach; tracts without reads; the centre's
round trip; the summary; the off-diagonal curvature (all-off-diagonal-zero survived, moving class shares up to
18%).

## 9. What's good

- The curvature of the total, not the mean, and the delta method, confirmed by an independent rebuild to 0.3%.
- The errors' spread over redrawn strata matches them (0.83 to 1.22) at 3 classes.
- Grouping the evaluations by slippage moves: 19 table rebuilds for 513 evaluations.
- The errors cost 10.8 to 12.3% of the climb at 13 classes, as spec §4.2 estimated.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib parameter_estimation::joint::ssr_fit`.
