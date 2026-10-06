# Fix Application Report: fit_precision_b1_2026-09-29.md

**Date:** 2026-09-29
**Source review:** `doc/devel/reports/reviews/fit_precision_b1_2026-09-29.md`
**Source state reviewed against:** `f72349b5` (a review object of the working tree, parent `c37e14e9`)
**Execution mode:** non-interactive, within the owner's approval of the Newton rule
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 4
- Minors: 12
- Nits: 5 (grouped)

### Outcome totals
- Applied: 4 (M3, M4, Mi7, Mi11)
- Applied with adaptation: 5 (M1, M2, Mi1, Mi2, Mi12)
- Already fixed: 1 (Mi9, by M1/M2's rewrite)
- Deferred: 5 (Mi3 in part, Mi4, Mi8, Mi10, the Nits)
- Disputed: 1 (Mi6 in part)
- Applied in part: Mi5 (the stale comment and spec §8; the rename deferred), Mi3 (one sample added)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --release --lib parameter_estimation::joint::fit` → 0, 94 passed, 1 ignored
- `cargo test --all-targets --all-features --no-fail-fast` → the numbers are in the commit message and the
  implementation report; the three pre-existing failures (`examples/ng_generic_loci_dump.rs` × 2,
  `examples/ng_ssr_loci_dump.rs` × 1) remain
- `cargo doc --no-deps`, `cargo audit` → not run (no public API or dependency change)
- Performance check → not applicable (no bench covers the stop rule); the information pass's cost was measured
  in review (1.48 to 1.80 times a plain pass)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | a rate the golden section leaves beside its floor is never held | Apply | Applied with adaptation | `fit/standard_errors.rs` | Pass |
| M2 | Major | a parameter walked towards a bound is never held | Apply | Applied with adaptation | `fit/standard_errors.rs` | Pass |
| M3 | Major | holding at a bound has no test that can fail | Apply | Applied | `fit/standard_errors.rs`, `fit/information.rs` | Pass |
| M4 | Major | the progress line does not say how far the fit is from stopping | Apply | Applied | `fit.rs`, `fit/settled.rs` | Pass |
| Mi1 | Minor | the last-pass cost ratio is diluted | Apply | Applied with adaptation | `fit.rs` | Pass |
| Mi2 | Minor | the distance cap changes no verdict; the spec's example for it is wrong | Apply | Applied with adaptation | `fit/standard_errors.rs`, spec | Pass |
| Mi3 | Minor | Newton tests miss no-information, not-identified, one sample | Apply in part | Applied in part / Deferred | `fit/standard_errors.rs` | Pass |
| Mi4 | Minor | the solves repeat the errors' parameter selection | Defer | Deferred | None | N/A |
| Mi5 | Minor | leftovers of the replaced rule | Apply in part | Applied in part / Deferred | `fit.rs`, spec | Pass |
| Mi6 | Minor | the stop report's state and one attribution | Apply / Dispute | Applied / Disputed | stop report | N/A |
| Mi7 | Minor | `settled_fraction`'s doc at zero, negative, NaN | Apply | Applied | `fit.rs` | Pass |
| Mi8 | Minor | a bare boolean where `PassKeeps` exists | Defer | Deferred | None | N/A |
| Mi9 | Minor | `values_in_layout` by position | Already fixed | Already fixed | — | Pass |
| Mi10 | Minor | the judging gate and "every parameter settled" untested at the loop | Defer | Deferred | None | N/A |
| Mi11 | Minor | an assertion without its invariant | Apply | Applied | `fit/settled.rs` | Pass |
| Mi12 | Minor | the thread-width test never reaches an in-loop information pass | Apply | Applied with adaptation | `fit.rs` | Pass |
| Nits | Nit | names | Defer | Deferred | None | N/A |

## 3. Questions asked and answers

None. M1 and M2 needed a design choice the owner's approval already covered: the owner approved "a parameter at a
bound whose step points outward is held there, the rest solved again"; the fix makes that rule hold for a parameter
*approaching* its bound as well, which is the rule's intent, and the spec amendment says so.

## 4. Per-finding log

### M1 and M2 — parameters at or approaching a bound
- **Final status:** Applied with adaptation (one fix for both).
- **Implementation:** `newton_step` now takes the Newton step **within the parameters' bounds** (`bounded_step`):
  the maximum of the quadratic model `g·d − ½ dᵀ I d` in the box, by the active-set rule — a free parameter whose
  step would carry it past an end is fixed there, its distance the way there, and the rest solved again with its
  move taken out of their scores (`information_times`); a fixed one whose slope at the solution points back in is
  freed. `AT_A_BOUND`, `pushes_past_its_bound` and the interval cap are gone. The solves work on one flat layout
  (`Layout`), and take the slope as an argument.
- **Adaptation:** M1's suggested fix (the golden section returning the end) was not needed: a rate
  4.4 × 10⁻¹⁰ above its floor with its step pointing out is fixed at the floor with a distance of 4.4 × 10⁻¹⁰.
  The M-steps' arithmetic is untouched.
- **Verification:** `a_fit_whose_maximum_lies_on_a_bound_converges` (the reviewer's six cohorts) — every one
  converges: excesses at 0, 21, 57 and 27 passes; clean rates below their floor, 15 and 45; the second shape above
  its bound, 36. The reviewer measured the last three at 198 passes, not converged, on `f72349b5`.
- **Tests:** `a_parameter_carried_past_an_end_stops_there` (a cohort share, a density shape and a sample's excess,
  both matrices, both ends reached, a quarter of the step inside and exactly on the end),
  `a_parameter_on_an_end_stepping_inwards_is_free`, `a_rate_the_golden_section_leaves_beside_an_end_stops_there`,
  and the fit-level test above; the two superseded tests removed.

### M3 — holding at a bound untested
- **Final status:** Applied, by M1/M2's tests (the fit-level test fails when bounds are ignored: the reviewer's
  M8 turned three of its cohorts into runs to the limit).

### M4 — the progress line
- **Final status:** Applied. Once judged, the line says how many of the parameters are not yet within the settled
  fraction of an error of the maximum, and the furthest by name with its distance in errors
  (`settled::furthest`, tested); before, the threshold at which judging starts. Plan step B3's per-start record
  stays B3's.

### Mi1 — the cost ratio
- **Final status:** Applied with adaptation. The alternation counts the passes that also sum the information and
  their time; the start line compares the final pass with the average *plain* pass. Wording changed to "the
  average plain pass before it"; its test updated.

### Mi2 — the cap
- **Final status:** Applied with adaptation: the cap is removed (the bounded step never leaves an interval), and the
  spec text for it replaced by the bounded step's; "0.11" now carries its unit where it survives, and the "34 times"
  names the invariant share.

### Mi3 — Newton tests' coverage
- **Final status:** Applied in part: `at_one_sample_the_newton_step_leaves_the_held_excess_out` (both matrices).
  **Deferred:** a no-information row and a not-identified parameter in the solve (the reviewer's test 5), to B2,
  which reads the step as an endpoint.

### Mi4 — the repeated selection
- **Final status:** Deferred. With nothing fixed the selections match exactly; with a parameter fixed at an end the
  others are judged by their marginal error against a distance solved with it fixed, which can only make settling
  easier. A shared "one inverse, `errors()` and `solve()`" type is about half a day and must keep the errors' bits;
  recorded in the implementation report for the owner.

### Mi5 — leftovers
- **Final status:** Applied in part: `maximisation`'s doc and its stale comment rewritten; spec §8 notes that
  `MAX_CONTRACTION` and `ERROR_REFRESH_CYCLES` are no longer used. **Deferred:** renaming
  `log_likelihood_stillness` (a public config field; its doc now says it only starts the judging).

### Mi6 — the stop report
- **Final status:** Applied: a note at its top says it is superseded and points at what was built.
  **Disputed:** the attribution. The sweep printed `coordinate {judged}` after incrementing `judged`, so
  coordinates count from 1 among parameters with an error: coordinate 2 is the invariant share. The first probe on
  the 4 × 8 cohort shows it directly — the parameter 0.2341 errors apart was 0.87595 against 0.88517, invariant-share
  values (`tmp/fit_precision/b1/probe1.log`, `trajectory.py`).

### Mi7 — `settled_fraction`'s doc
- **Final status:** Applied: at zero, below it or not a number, the fit never stops before `max_passes` unless no
  parameter has both an error and a distance.

### Mi8, Mi10, Nits
- **Final status:** Deferred: small and local; Mi10 needs the per-cycle verdict in the trace, which B3 adds.

### Mi9 — `values_in_layout`
- **Final status:** Already fixed: `Layout::values_of` writes through the `cohort::` constants.

### Mi11 — the assertion
- **Final status:** Applied: a comment names why the two lengths agree.

### Mi12 — the thread-width test
- **Final status:** Applied with adaptation: the test also runs a fit with a settled fraction of 10⁹, so its first
  judged cycle — after an in-loop information pass — stops it, and requires the same bits at one thread and four.
  A 3-sample cohort does not converge within 200 passes at 0.1, so the wide fraction is what makes the stop
  reachable.

## 5. Deferred findings to carry forward
- Mi3 (no-information and not-identified parameters in the solve) — to B2.
- Mi4 (one inverse type for errors and solves) — for the owner.
- Mi5 (rename `log_likelihood_stillness`), Mi8, Nits — small.
- Mi10 — with B3's per-cycle record.

## 6. Disputed findings to return to reviewer
- Mi6's attribution — above.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no bench covers the stop rule; the information pass's cost was measured in review.

## 10. Commands run
- `scripts/dev.sh cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`
- `scripts/dev.sh cargo test --all-targets --all-features --no-fail-fast`
- `scripts/dev.sh cargo test --release --lib the_errors_mean_what_they_say -- --ignored`
- the oracle, `tmp/fit_precision/b1/run_oracle_b1f.sh`

## 11. Command results
In the implementation report, [fit_precision_b1_2026-09-29.md](../implementations/fit_precision_b1_2026-09-29.md) §5.

## 12. Notes
- B2 inherits a finding from the correctness review: at 2 and 4 samples the full Newton step is several times too
  long (right direction); "current value plus step" as a projected endpoint needs care there.
