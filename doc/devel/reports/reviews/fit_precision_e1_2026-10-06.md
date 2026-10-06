# Code Review: fit_precision_e1
**Date:** 2026-10-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step E1 — `Estimate` gains `standard_error: Option<f64>`, and the fits' errors are carried onto it
**Status:** Approve-with-changes

---

## 1. Scope

- **Reviewed:** the diff of review object `141e5a7c` against `a1fea02a` (branch `fit-precision`), plan step E1 of
  [fit_precision.md](../../implementation_plans/fit_precision.md); spec [fit_precision.md](../../ng/spec/fit_precision.md) §5.1.
- **In scope:** [src/parameter_estimation/mod.rs](../../../../src/parameter_estimation/mod.rs),
  [src/parameter_estimation/joint/fit.rs](../../../../src/parameter_estimation/joint/fit.rs),
  [src/run/census_fit.rs](../../../../src/run/census_fit.rs); nine further files that only add `standard_error: None`
  to `Estimate` literals.
- **Out of scope:** the parameters file (plan step E2), the repeat-tract slippage errors (they live on `StratumFit`).
- **Dispatched**, two reviewers, each in its own worktree at `141e5a7c`:
  - reliability, float portability, and whether the diff matches its intent — the step moves numbers that reach a later
    output;
  - errors, defaults, refactor safety, naming, idiomatic Rust, smells — a public type gains a field.

## 2. Verdict

Approve-with-changes.

## 3. Execution status

- `cargo fmt --check`: clean. `cargo clippy --all-targets --all-features -- -D warnings`: clean.
- `cargo test --release --lib -- parameter_estimation::joint::fit:: run::census_fit parameter_estimation::tests`:
  122 passed, 1 ignored.
- Full suite (`cargo test --all-targets --all-features --no-fail-fast`): 5,025 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2 and `examples/ng_ssr_loci_dump.rs` × 1), 9 ignored.
- Mutations: 13 run across both reviewers; 10 killed; 3 changed no behaviour (M2, M3); none survived with a behaviour change.
- Needs verification: none.

## 4. Open questions and assumptions

1. Whether a counted rate with no mismatch should carry no error or an upper bound such as 3/n (M1). Resolved in the
   fixes as *no error*, the reading of spec §5.1's "None where nothing determined the value"; recorded for the owner
   at checkpoint E.

## 5. Top 3 priorities

1. M1 — a substitution rate counted with no mismatch gets an error of exactly zero.
2. M2 — nothing tests the error each library's rate carries out of `census_fit`.
3. M3 — that error has no reader yet; E2 needs a route for it into the base-quality calibration.

## 6. Findings

### Major

**M1: [src/run/census_fit.rs:392](../../../../src/run/census_fit.rs#L392) — a rate counted with no mismatch gets an
error of zero.** Confidence High. **Categories:** reliability, smells (convergent). √(p(1 − p)/n) is zero at p = 0
and p = 1, and `ErrorRate` admits both. The fixture's own stratum (460 bases, none mismatching) produced `Some(0.0)`;
E2 would write *known exactly* beside a rate from a count that found one outcome only. **Fix:** `None` when no base
mismatched or every one did; a test for each.

**M2: [src/run/census_fit.rs:346](../../../../src/run/census_fit.rs#L346) — no test sees the error carried onto each
library's `Estimate<ErrorRate>`.** Confidence High. Dropping it and giving every library the first library's error both
passed every test: in the one fixture reaching this code all three read groups' errors are `None`. **Fix:** move the
conversion into the fit, where the two-library test can check it against a fit whose libraries' errors differ.

**M3: [src/run/census_fit.rs:346](../../../../src/run/census_fit.rs#L346) — the library error has no reader.**
Confidence High. It goes to `RunParameters::assemble`, whose `ReadGroupCalibration::from_fitted_rate` keeps only the
value and warrant, and `parameters_file_of` passes `ReadsBehindEachCalibration::nothing_was_fitted`. A mutation to
`Some(-1.0)` passed 1,198 tests across `run::`, `calling::parameters_file` and `cli::`. **Fix (E2):** carry the
multiplier's error, the rate's error over the mean minted error, beside the multiplier.

### Minor

**Mi1: [src/parameter_estimation/joint/fit.rs:292](../../../../src/parameter_estimation/joint/fit.rs#L292) — two
public maps that must hold the same read groups.** `clean_rate_standard_error` beside `noise`, nothing enforcing it, and
the reader's `.get(..).flatten()` turns a missing group into *no error*. **Fix:** private field, read only through a
method that panics on a mismatch.

**Mi2: [src/parameter_estimation/mod.rs:136](../../../../src/parameter_estimation/mod.rs#L136) — the doc's reasons for
`None` are incomplete**: *not told apart* and *wider than its range* are missing, both reachable; the comment at the
homozygote excess says only *at one sample*; "derived from others get None" is contradicted by the inbreeding
coefficient, which keeps the excess's error.

**Mi3: [src/parameter_estimation/mod.rs:140](../../../../src/parameter_estimation/mod.rs#L140) — a bare
`Option<f64>` admits a negative, zero or NaN error.** **Fix:** a checked type.

### Nits

- [census_fit.rs:1030](../../../../src/run/census_fit.rs#L1030): "about three fifths of the rate" is 0.577, the rate
  over √3.
- [census_fit.rs:560](../../../../src/run/census_fit.rs#L560): the map is rebound mutable to set one error.
- [fit.rs:2075](../../../../src/parameter_estimation/joint/fit.rs#L2075): another bare `0` for *ordinary position*.

## 7. Out of scope observations

- [census_fit.rs:506](../../../../src/run/census_fit.rs#L506): the doc of `parameters_file_of` says this route has no
  per-library minted totals, while `parameters_from_the_fit` builds them. For E2, which needs them.
- [to_run_parameters.rs:355](../../../../src/calling/parameters_file/to_run_parameters.rs#L355) and `:503`: `None` is
  right for a version-1 file; E2 must read the new key there or a run that reads a file and writes it back loses it.
- The two-library test pairs `fit.rates.keys()` (name order) with sample indices, which holds for its names `s0`…`s5`
  and not for `s10`.

## 8. Missing tests to add now

- `a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own`: extend with a zero-mismatch and
  an all-mismatch stratum, both without an error (M1).
- `a_samples_two_libraries_come_back_at_their_own_rates`: assert each library's carried error is its own (M2).

## 9. What's good

- The library error is found through the sample and section it belongs to
  ([fit.rs:2071](../../../../src/parameter_estimation/joint/fit.rs#L2071)); the two-library test killed the wrong-class,
  wrong-library and wrong-sample mutations.
- The inbreeding coefficient takes the excess's error with a one-line reason
  ([census_fit.rs:476](../../../../src/run/census_fit.rs#L476)).

## 10. Commands to re-verify

- `./scripts/dev.sh cargo test --release --lib -- parameter_estimation::joint::fit:: run::census_fit parameter_estimation::tests`
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings`

Per-category files: `tmp/review_2026-10-06_fit_precision_e1/`.
