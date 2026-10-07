# Fix Application Report: fit_precision_e1_2026-10-06.md

**Date:** 2026-10-06
**Source review:** `doc/devel/reports/reviews/fit_precision_e1_2026-10-06.md`
**Source state reviewed against:** review object `141e5a7c` (parent `a1fea02a`, branch `fit-precision`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3
- Minors: 3
- Nits: 3

### Outcome totals
- Applied: 2 (M1, Mi2)
- Applied with adaptation: 2 (M2, Mi1)
- Deferred: 2 (M3 to E2, Mi3)
- Nits: one applied (the √3 wording), two not

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --release --lib -- parameter_estimation::joint::fit:: run::census_fit parameter_estimation::tests` → 0, 122 passed, 1 ignored
- `cargo test --all-targets --all-features --no-fail-fast` → the three pre-existing failures only (the step's implementation report gives the counts)
- `cargo doc`, `cargo audit` → not run; no dependency or public doc link changed beyond the fields' own docs
- Performance check → not applicable: no hot path changed

### Unresolved high-priority findings
- M3 — deferred to plan step E2, whose subject it is.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed |
|---|---|---|---|---|---|
| M1 | Major | zero error on a one-outcome count | Apply | Applied | `src/run/census_fit.rs`, `src/parameter_estimation/mod.rs` |
| M2 | Major | library error untested | Apply | Applied with adaptation | `src/parameter_estimation/joint/fit.rs`, `src/run/census_fit.rs` |
| M3 | Major | library error has no reader | Defer | Deferred | None |
| Mi1 | Minor | two maps that must agree | Apply | Applied with adaptation (with M2) | `src/parameter_estimation/joint/fit.rs` |
| Mi2 | Minor | incomplete `None` reasons | Apply | Applied | `src/parameter_estimation/mod.rs`, `fit.rs` |
| Mi3 | Minor | bare `Option<f64>` | Defer | Deferred | None |

## 3. Questions asked and answers

None asked. M1's choice between *no error* and an upper bound is recorded for the owner at checkpoint E.

## 4. Per-finding log

### M1 — a rate counted with one outcome only gets an error of zero
- **Final status:** Applied.
- **Implementation:** the binomial error is computed only when some but not every compared base mismatched; otherwise
  `None`. The `Estimate` doc names the case. Spec §5.1 already says `None` where nothing determined the value; the spec
  is not edited without the owner.
- **Tests:** `a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own` now asserts the
  fixture's own zero rate (460 bases, none mismatching) and an added stratum of 5 bases all mismatching carry no error.
  Mutation `<` → `<=` on the upper condition: killed by that test.

### M2 and Mi1 — the library error, tested and kept with its rate
- **Final status:** Applied with adaptation. The review suggested a helper in `census_fit`; the conversion moved into the
  fit instead, as `JointFit::sequencing_error_rates`, which returns each read group's ordinary-position rate as an
  `Estimate<ErrorRate>` with its error. The error map is private and the method panics if a read group has rates and no
  error entry, which also settles Mi1. `census_fit` calls it.
- **Tests:** `a_samples_two_libraries_come_back_at_their_own_rates` asserts each library's returned error equals its own
  slot's and its value equals the fit's rate. Mutation dropping the error in the method: killed by that test.

### M3 — the library error has no reader
- **Final status:** Deferred to plan step E2, which writes the multiplier's error. Noted there: the multiplier is the
  rate over the mean minted error, so its error is the rate's over the same; and `parameters_file_of`'s doc on minted
  totals needs reconciling with `parameters_from_the_fit`.

### Mi2 — the reasons for `None`
- **Final status:** Applied: the `Estimate` field doc lists every reason, including not told apart, wider than its range
  and the one-outcome count, and says a value that is another under a second name keeps that value's error; the field doc
  in `fit.rs` and the excess's comment no longer imply one reason.

### Mi3 — a checked type for the error
- **Final status:** Deferred. Spec §5.1 fixes the field as `Option<f64>`, and every value written into it comes from a
  square root of a positive number or the fit's `StandardError::Estimated`, which holds only placed, finite errors. A
  checked type would be a spec change for the owner.

### Nits
- The √3 wording: applied. The mutable rebinding and the bare `0`: not changed, being one test line and a call whose
  argument the surrounding code already spells the same way.

## 5–8.
Deferred: M3 (E2), Mi3 (owner). Disputed, failed validation, blocked: none.

## 9. Performance check
Skipped — no `Apply` touched a hot path.

## 10–11. Commands run and results
- `./scripts/dev.sh cargo clippy --all-targets --all-features -- -D warnings` → 0
- `./scripts/dev.sh cargo test --release --lib -- parameter_estimation::joint::fit:: run::census_fit parameter_estimation::tests` → 0, 122 passed
- two mutations planted and restored (`tmp/fit_precision/e1/mut/`), each killed

## 12. Notes
The two-library test pairs name order with sample index; it holds for its six names and is left as the review recorded it.
