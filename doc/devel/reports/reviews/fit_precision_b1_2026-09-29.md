# Code Review: fit_precision_b1
**Date:** 2026-09-29
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step B1 — the SNP/indel fit stops when every parameter is within a tenth of its error of the maximum, by a Newton step
**Status:** Request-changes (applied: [fixes_applied_fit_precision_b1_2026-09-29.md](fixes_applied_fit_precision_b1_2026-09-29.md))

---

## 1. Scope

- **What was reviewed:** the diff `c37e14e9..f72349b5` (a review object, not on the branch): `fit.rs`, the new
  `fit/settled.rs`, `fit/standard_errors.rs` (`newton_step`, its solves, `in_coordinate_order`),
  `fit/information.rs` (`InformationSums::slope`, the new tests), the spec amendment (§2, §3.3, §3.4), the plan's
  note under B1, and the stop report `fit_precision_b1_stopped_2026-09-29.md`.
- **Out of scope:** the repeat-tract fit; the merge of main (`c37e14e9`), checked separately (both cross-platform
  checksums unchanged on it alone).
- **Categories:** three reviewers, each in its own worktree —
  *correctness and cost* (reliability of behaviour, float_portability, unsafe_concurrency, extras/cost);
  *design and claims* (naming, idiomatic, defaults, module_structure, smells, refactor_safety, diff against
  intent, every number in the prose); *tests and errors* (reliability: tests, mutation testing; errors).
  Findings files: `tmp/review_2026-09-29_fit_precision_b1/{correctness_cost,design_claims,reliability_errors}.md`,
  evidence beside them.

## 2. Verdict

Request-changes: the rule never settles a fit whose maximum lies on a parameter's bound (M1, M2), and nothing
tested that case (M3).

## 3. Execution status

- The reviewers ran the fit module's tests (90 passed, 1 ignored), 44 mutations between them (34, 7 and 3), probes
  on drawn cohorts (Newton steps on real scores, bound cohorts, determinism at 1, 4 and 8 threads, pass cost), and
  re-derived the stop report's tables from the logs.
- `cli::cross_platform_digests` was expected to fail on this commit (fitted numbers move) and was not reported.
- Needs verification: none open.

## 4. Open questions and assumptions

1. With a parameter held at a bound, each other parameter is judged by its error from the whole inversion
   against a distance solved with the held one fixed (correctness, design). The spec does not say which is meant;
   the pairing can only make settling easier. Carried as Mi4.

## 5. Top 3 priorities

1. **M1/M2** — make the Newton step respect the parameters' bounds properly (a box-constrained step), so fits whose
   maximum lies on a bound converge.
2. **M3** — a fit-level test with maxima on bounds, and unit tests at both ends and for cohort-level parameters.
3. **M4** — the progress line must say how far the fit is from stopping and which parameter holds it.

## 6. Findings

### Major

**M1: `standard_errors.rs:662` — A clean error rate the golden section leaves beside its floor is never held.**
**Categories:** tests and errors. **Confidence:** High.
`maximise_error_rate`'s golden section stops within about 10⁻⁹ of an end and returns the bracket's middle; a rate
whose maximum is its floor came back at 1.000437 × 10⁻⁶, 4.4 × 10⁻¹⁰ above the floor of 10⁻⁶, twice. The hold
tolerance for that interval was 2 × 10⁻¹⁰, so the rate's step past the floor was never zeroed. Cohorts drawn with a
clean rate of 10⁻⁷: 20 samples × 3 reads, 198 passes, 61 parameters unsettled; 4 × 8, 198 passes, 9 unsettled.
Fix: a step that would cross an end stops there (M2's fix); or make the golden section return the end.

**M2: `standard_errors.rs:600-626` — A parameter the fit walks towards a bound is never held; its step overshoots.**
**Categories:** tests and errors; correctness (cross-category). **Confidence:** High.
Only a parameter already on an end was held. On the A7 cohort (second shape drawn at 150), `density_a` drifts from
0.105 to 0.059 towards its floor of 0.02 with a step of −1.48; the fit is unsettled at 198 passes (8 parameters)
and at 999 (6). Holding the crossing steps with a distance of `end − value` alone broke another cohort; the fix is
the active-set rule for a box, choosing holds by the slope's sign.

**M3: fit-level tests — Holding at a bound has no test that can fail.**
**Categories:** tests and errors; design and claims. **Confidence:** High.
Setting the hold tolerance to zero leaves all 90 fit-module tests green while three excess-at-zero cohorts go from
converging (21, 57, 27 passes) to the limit. Deleting the upper-end clause survives; so do dropping the cohort
parameters' holds from either solve and swapping `a` and `b` in `values_in_layout` (design reviewer's probe kills
all three).

**M4: `fit.rs:2158-2175` — The progress line does not say how far the fit is from stopping.**
**Categories:** tests and errors; design and claims; correctness and cost. **Confidence:** High.
It printed "judged this cycle"; the count not settled was computed and dropped, and no parameter was named. A fit
that runs to the limit — the kimura case B1 exists for — gives no reason.

### Minor

**Mi1: `fit.rs:2232` — the start line's "last pass against the average pass" is diluted** by the information passes
now inside the loop: 1.30 printed where the direct ratio is 1.48 (4 samples), 1.56 against 1.70 (20).
*(correctness and cost)*

**Mi2: spec §2 amendment, `standard_errors.rs:575-577` — the distance cap cannot change a verdict**, since a
reported error is never wider than its interval; the spec's example for it (6.3 errors on the duplicated share)
is 0.031 against an interval of 0.05, not cut; and "moved by 0.11" lacks its unit (0.11 of the 1,000-pass fit's
error). *(all three)*

**Mi3: the Newton tests never cover a parameter without information, one not identified, or one sample.**
*(tests and errors)*

**Mi4: `arrow_step`/`whole_matrix_step` repeat the parameter selection of `of_the_blocks`/`of_the_whole_matrix`**
(four copies of one entry rule); with nothing held they match exactly. *(design and claims)*

**Mi5: leftovers of the replaced rule** — `log_likelihood_stillness` now only starts the judging (rename);
`maximisation`'s floors and a stale comment about "reporting convergence"; spec §8 still lists `MAX_CONTRACTION`
and `ERROR_REFRESH_CYCLES`. *(design and claims)*

**Mi6: the stop report states the pre-decision state** ("not committed") and, per the design reviewer, attributes
the parameters left behind to the invariant share where the fixed non-reference share was (medium confidence).
*(design and claims)*

**Mi7: `settled_fraction`'s doc is imprecise at zero; a negative or NaN fraction silently never stops.**
*(design and claims)*

**Mi8: `step(&x, false)` passes a bare boolean where `PassKeeps` exists.** *(design and claims)*

**Mi9: `values_in_layout` lists the cohort slots by position, not through `cohort::`.** *(design and claims)*

**Mi10: the judging gate, and "stop only when every parameter is settled", have no test that fails** when changed
(starting judged, the gate resetting each cycle, stopping with one parameter unsettled: 15 passes against 18).
*(design and claims; tests and errors)*

**Mi11: `settled.rs:42` — `assert_eq!` in production code without the invariant that keeps it from firing.**
*(tests and errors)*

**Mi12: the thread-width test runs 3 passes, so it never reaches an information pass inside the loop.**
*(correctness and cost)*

### Nits

`not_settled` returns a count (`count_not_settled`); `held_none` means "none newly held this round";
`AT_A_BOUND` lacks its noun; `judging` could be `is_judging`; the settled-by-definition justification in
`settled.rs` is loose for an error wider than its range.

## 7. Out of scope observations

- The oracle's two window-coverage checksums differ from `scripts/promote_ng_oracle.baseline` since the merge of
  main: main's psp format 1.1 counts the reads the depth cap discarded into window depth, and did not re-record
  the baseline. Re-record in its own commit.
- B2: at 2 and 4 samples the full Newton step points the right way but is several times too long (it lands 39 to
  331 log-likelihood units below the 600-pass fit, from 0.05 to 11.9 before it). "Current value plus step" as a
  projected endpoint needs care there. *(correctness and cost)*

## 8. Missing tests to add now

`a_parameter_pushed_past_either_bound_is_held_there`, `a_clean_rate_the_golden_section_leaves_at_its_floor_is_held`,
`a_fit_whose_maximum_lies_on_a_bound_converges`, `at_one_sample_the_newton_step_leaves_the_held_excess_out`, a
cohort-parameter hold on both matrices (design reviewer's probe), `a_parameter_the_newton_step_cannot_solve_has_no_distance`,
`a_one_sample_fit_converges`, and the thread-width test run to convergence. Code for most is in the findings files.

## 9. What's good

- `arrow_step` is the exact solve of the arrow system, verified by derivation and against a dense solve.
- The summed slope, the Newton step and a whole judged fit give the same bits at 1, 4 and 8 threads (3, 25 and 4
  samples with the duplicated class, over 13 chunks).
- On real scores at 20 and 40 samples a Newton step from pass 6 cuts the log-likelihood gap to a 600-pass fit from
  0.045 to 0.004 and from 0.066 to 0.004 units: the distance it measures is the real one.
- `newton_step` costs 1.6 ms at 2,000 one-library samples, negligible beside a pass.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`;
`scripts/dev.sh cargo test --release --lib the_errors_mean_what_they_say -- --ignored`;
the oracle, `tmp/fit_precision/b1/run_oracle_b1f.sh`.
