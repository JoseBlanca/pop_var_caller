# Fix Application Report: fit_precision_a3_2026-09-27.md

**Date:** 2026-09-27
**Source review:** `doc/devel/reports/reviews/fit_precision_a3_2026-09-27.md`
**Source state reviewed against:** `9943a969` (review-only commit), branch `fit-precision`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 6
- Minors: 6
- Nits: 2

### Outcome totals
- Applied: 11
- Applied with adaptation: 2
- Deferred: 2
- Disputed: 0
- Failed validation: 0
- Awaiting user answer: 0 (two items carried to checkpoint A)

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --all-targets --all-features --no-fail-fast` → counts in the commit message
  (`tmp/fit_precision/suite_a3b.log`); the three pre-existing failures only; both
  `cli::cross_platform_digests` tests pass unchanged
- `cargo test --release --lib parameter_estimation::joint::fit` → 0, 51 passed
- Performance check → not applicable (nothing calls the errors outside tests)

### Unresolved high-priority findings
- M1 — the block approximation understates the samples' errors at small cohorts; A4 measures, the
  owner decides at checkpoint A.

## 2. Findings table

| ID | Severity | Title | Decision | Final status | Files changed | Follow-up |
|---|---|---|---|---|---|---|
| B1 | Blocker | no errors, or millions, at 1–2 samples | Apply | Applied with adaptation | `standard_errors.rs`, `information.rs` | Checkpoint A |
| M1 | Major | the gap is chance only at large cohorts | Apply | Applied (prose); decision deferred | test doc, report | Checkpoint A, A4 |
| M2 | Major | write-back untested | Apply | Applied | `standard_errors.rs` | No |
| M3 | Major | a failed sample block drops the sample | Apply | Applied (superseded by B1's rule) | `standard_errors.rs` | No |
| M4 | Major | threshold pinned from one side | Apply | Applied | `standard_errors.rs` | No |
| M5 | Major | `c == 0` untested | Apply | Applied | `standard_errors.rs` | No |
| M6 | Major | comparison doc overstates | Apply | Applied | `information.rs` | No |
| Mi1 | Minor | non-finite diagonal untested | Apply | Applied | `standard_errors.rs` | No |
| Mi2 | Minor | misplaced doc | Apply | Applied (function rewritten) | `standard_errors.rs` | No |
| Mi3 | Minor | tolerances keyed on strings | Apply | Applied | `information.rs` | No |
| Mi4 | Minor | `SampleSolve` naming | Apply | Applied | `standard_errors.rs` | No |
| Mi5 | Minor | spec shows the rejected formula | Defer | Deferred | None | Owner |
| Mi6 | Minor | wrong numbers in prose | Apply | Applied | test doc, report, module doc | No |

## 3. Questions asked and answers

None asked. Carried to checkpoint A: M1 (with A4's measurement), Mi5.

## 4. Per-finding log

### B1 — the smallest cohorts
- **Final status:** Applied with adaptation. `invert_identified` builds the Cholesky factor in order
  and drops a parameter whose pivot keeps less than `IDENTIFIED_SHARE` = 10⁻⁸ of its reference
  curvature, then refactors without it — in the cohort's reduced block (reference: the cohort's
  diagonal before the samples' share is taken) and in each sample's own block. `StandardError` says
  why an error is absent: `NoInformation`, `HeldFixed`, `NotIdentified`, or `WiderThanItsRange`.
- **Adaptation:** the reviewer's fix alone left, at two samples, errors of 3.6 on the invariant share
  and 720, 1,148 and 1,591 on three Beta shapes — over parameters the full matrix does not invert at
  all. An error wider than the whole interval the fit keeps its parameter in is therefore reported as
  `WiderThanItsRange` (keeping the width for reporting), and the parameter stays in the inversion so
  its uncertainty still widens the others'. Not in the spec; recorded for checkpoint A. A warning
  when a parameter is dropped is left to A5, which prints the counts by reason.
- **Tests:** `the_smallest_cohorts_keep_their_error_rates` (one and two samples, 30,000 positions: every
  error rate keeps its error, no share's error reaches 1); `a_parameter_the_others_mimic_is_dropped_alone`;
  `an_error_wider_than_its_range_is_not_an_error`.

### M1 — the reading of the gap
- **Final status:** Applied for the prose; the decision is deferred. The comparison's doc and the report
  say the gap is chance at 20 samples and 8 reads, and not at 4 samples or 3 reads, where the blocks'
  errors are too small, with the review's measurements. A4's coverage is to be measured there.

### M2 — write-back
- **Final status:** Applied. `a_parameter_without_information_mid_block_leaves_the_others_in_place`
  and `a_diagonal_that_is_not_finite_is_no_information` remove parameters from the middle of the
  cohort's and a sample's block. A planted slip — each cohort error written to its position among the
  kept parameters — fails both (`tmp/fit_precision/mutant_a3.log`); the file was restored and checked.

### M3 — a failed sample block
- **Final status:** Applied, by B1's rule: a sample's parameter it cannot tell apart is dropped alone
  and the sample stays in the cohort's sum. `a_sample_parameter_it_cannot_tell_apart_is_dropped_alone`.

### M4 — the threshold
- **Final status:** Applied. `the_threshold_sits_between_poor_information_and_rounding` keeps
  `1 − R²` = 10⁻⁶ and drops 10⁻¹⁰; `a_parameter_in_small_units_is_still_identified` checks it is
  relative to the parameter's own curvature.

### M5 — no cohort information
- **Final status:** Applied. `without_cohort_information_the_samples_keep_their_own_errors`.

### M6 — the comparison's doc
- **Final status:** Applied. It now says it pins each kind's measured range, not the algebra, which the
  unit tests pin.

### Mi1–Mi4, Mi6
- Applied as listed in the table. Mi3: a `Kind` enum; Mi4: `SampleBlockInverse` with
  `own_inverse_times_cross`; Mi6: every figure corrected — +0.000 to +0.001; −0.166 at 4 samples, a
  sixth; the threshold described by what it keeps (`1 − R²`), not by a multiple of rounding; the
  "slightly smaller" rule is gone with the rule itself.

### Mi5, Nits
- Deferred to the owner: the spec's §3.2 text, and 63 samples (the plan names 4 and 20).
- `error_from_variance` at a variance of exactly 0 returns `NotIdentified`; unreachable, left.

## 5. Deferred findings to carry forward
- M1's decision (checkpoint A, with A4).
- Mi5 (owner).

## 6.–8.
None disputed, failed or blocked.

## 9. Performance check
- Not triggered.

## 10. Commands run
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit -- --nocapture`
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit::standard_errors` (with
  the planted slip, then restored)
- `scripts/dev.sh bash -c 'cargo fmt --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all-targets --all-features --no-fail-fast'`

## 11. Command results
- fit module → 0, 51 passed; `standard_errors` → 0, 13 passed; with the planted slip → 2 failed
- fmt, clippy, full suite → see §1 and the commit message

## 12. Notes
- Two temporary probes in `information.rs` (printing the two-sample and four-sample errors) were
  applied from a backup copy and restored from it; the restore was checked against the tree.
