# Code Review: fit_precision_f3
**Date:** 2026-10-09
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step F3 as first built — the Beta shapes' step capped; a mid-fit fallback to the best point; a start returning its best point
**Status:** Request-changes (applied, and one change dropped by the owner — see [fixes_applied_fit_precision_f3_2026-10-09.md](fixes_applied_fit_precision_f3_2026-10-09.md))

---

## 1. Scope

- **Reviewed:** commit `8fce5fe4` (temporary, folded into the step's commit) against `b94dc466`, branch
  `fit-precision-followup`.
- **In-scope files:** `src/parameter_estimation/joint/fit.rs` (`maximise`, `StartOutcome`, the start's log line,
  `step_beta_shapes`, the new tests), `src/parameter_estimation/joint/fit/information.rs` (the new property test);
  the plan section for consistency.
- **Categories:** two reviewers in their own worktrees — reliability, errors, refactor_safety and float_portability
  with mutation testing (11 mutations); naming, idiomatic, smells, defaults and extras with every stated number
  checked. Per-category files: `tmp/review_2026-10-09_f3_fit_never_keeps_a_loss/`.

## 2. Verdict

Request-changes. The cap does what the decision asked; the fallback, as built, can never let a start converge on an
existing test cohort; and no test reached either new branch.

## 3. Execution status

- `cargo test --release --lib joint::fit`: 118 passed. The four new or changed tests pass.
- Mutations: 11 run; 1 caught (removing the cap), 8 survived having changed the fit's output, 2 changed nothing.
- Findings labelled "Needs verification": 0.

## 4. Open questions and assumptions

1. Whether the fallback should stay — answered by the owner (2026-10-09): dropped (B1, M1).

## 5. Top 3 priorities

1. B1 — no test runs the fallback or the best-point return, though an existing 2-second fixture reaches both.
2. M1 — the fallback replays one cycle every three passes to the pass limit on that fixture.
3. M2 — going back to a jumped-to point retraces the step that lost.

## 6. Findings

### Blocker

**B1: fit.rs `maximise`; fit/information.rs — no test fails when the fallback or the best-point return is broken.**
`a_stretch_some_samples_carry_twice_is_not_read_as_heterozygosity` (30 samples, 4,000 positions, 3 reads, 120
passes) reaches both; the parent commit returns each of its six starts 1.5 to 16.4 units below its best. Returning
the best parameters with the last point's statistics survived every test. **Categories:** reliability, extras.

### Major

**M1: a plain cycle from the best point that loses is replayed exactly, every three passes.** On that cohort one
start fell back 37 times, 108 of its 118 passes on one repeated cycle; with the fallback off all six starts return
the same log-likelihoods. With the cap removed, five of six starts are unchanged: something other than the shapes'
update makes the alternation settle 1.4 to 3.7 units below a point it visited. **Categories:** reliability.

**M2: the go-back target can be a jumped-to point**, and on kimura that is where the losing step began, so going
back there and taking a plain cycle retraces it: the 12-pass loop would become a 3-pass stall. **Categories:** smells,
claims.

**M3: the cap clamps each shape of the joint step on its own** — the failure the A7 review found at a bound. On the
test's fixture the cap stops b at 17.2 and leaves a where it fits b = 0.02. *Fix:* project at the cap as at a bound.
**Categories:** smells.

### Minor

- The `maximise` doc still says the plain alternation never lowers the log-likelihood; the module doc, the
  `JUMP_SLACK` doc, the `trace_the_returned_fit` and `fit_trace` docs omit the new behaviour. **Categories:** smells,
  claims.
- A start can say "converged" while returning a best point it never judged settled. **Categories:** reliability.
- The last pass's time covers two passes when the best is returned. **Categories:** reliability, smells.
- The cap test asserts through the constant, so a factor of 3, or a cap on the second shape only, passes.
  **Categories:** reliability.
- `f64::clamp` panics on a shape that is not a number. **Categories:** errors.
- The best point is an anonymous tuple; `returned_its_best` reads as a yes/no; `plain_cycle` is not a question.
  **Categories:** idiomatic, naming.
- The cap test's doc comment sat on another test. **Categories:** smells, reliability.
- The decision's "step taken on the log scale" and "where the previous cycle started" are not what was built, and
  neither departure was recorded. **Categories:** extras, claims.

### Nits

The log figures carry no unit; `maximise` is about 295 lines; `MAX_SHAPE_FACTOR_A_PASS` gives no source; the bounds
clamp after the factor clamp cannot bind.

## 7. Out of scope observations

- Some update other than the Beta shapes' is not uphill on the 30-sample cohort (M1). Not traced; the next step.

## 8. Missing tests to add now

- `every_start_at_the_pass_limit_returns_the_best_point_it_reached` on the 30-sample cohort, every start, both
  fits (B1).
- The cap's factor pinned by value on the first shape (Minor).

## 9. What's good

- The new `StartOutcome` fields are destructured exhaustively where the outcome is consumed.
- The kimura figures quoted in the plan and the docs all check against the owner's report.

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --release --lib joint::fit`
