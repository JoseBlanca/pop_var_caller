# Fix Application Report: fit_precision_c1_2026-09-30.md

**Date:** 2026-09-30
**Source review:** `doc/devel/reports/reviews/fit_precision_c1_2026-09-30.md`
**Source state reviewed against:** `57618353` (a review object of the working tree, parent `efc325e2`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 6
- Minors: 8
- Nits: grouped

### Outcome totals
- Applied: 7 (B1, M1, M2, M3, M4, M6, Mi6)
- Applied with adaptation: 3 (M5, Mi1, Mi4)
- Applied in part: 1 (Mi7)
- Deferred: 5 (Mi2, Mi3 recorded in the report, Mi5, Mi8, the off-diagonal test)
- Nits: deferred

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --release --lib parameter_estimation::joint::ssr_fit` → 0, 42 passed
- `cargo test --all-targets --all-features --no-fail-fast` → 101, 4,992 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1), 5 ignored
- `cargo doc --no-deps`, `cargo audit` → not run (no dependency change)
- Performance check → not applicable (no bench covers the stratum errors; their cost is measured in the
  implementation report)

### Unresolved high-priority findings
- None. `NOT_PLACED`'s value goes to the owner at checkpoint C.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| B1 | Blocker | schedule tests never read the errors | Apply | Applied | `ssr_fit.rs` tests | Pass; mutation killed |
| M1 | Major | a number stuck at the climb's reach gets a tiny error | Apply | Applied | `ssr_fit.rs` | Pass |
| M2 | Major | a runaway concentration gets rounding noise | Apply | Applied | `ssr_fit.rs` | Pass |
| M3 | Major | four reasons folded into one `None` | Apply | Applied | `ssr_fit.rs` | Pass |
| M4 | Major | "+21%" is noise | Apply | Applied | impl report | N/A |
| M5 | Major | errors beside blended numbers | Apply | Applied with adaptation | `ssr_fit.rs` doc | N/A |
| M6 | Major | centring, empty tracts, summary untested | Apply | Applied | `ssr_fit.rs` tests | Pass; mutations killed |
| Mi1 | Minor | the step's doc and `max(1, |x|)` | Apply | Applied with adaptation | `ssr_fit.rs` | Pass |
| Mi2 | Minor | the identification floor for finite differences | Defer | Deferred | — | N/A |
| Mi3 | Minor | the drop order | Apply | Recorded (report §2) | impl report | N/A |
| Mi4 | Minor | the log line | Apply | Applied with adaptation | `ssr_fit.rs` | Pass |
| Mi5 | Minor | module structure | Defer | Deferred | — | N/A |
| Mi6 | Minor | the field doc's non-existent path | Apply | Applied | `ssr_fit.rs` | N/A |
| Mi7 | Minor | `.expect` without a note; empty starting points | Apply in part | Applied in part | `ssr_fit.rs` | Pass |
| Mi8 | Minor | errors start after every walk | Defer | Deferred | — | N/A |
| Nits | Nit | | Defer | Deferred | — | N/A |

## 3. Questions asked and answers

None during the run.

## 4. Per-finding log

### B1 — the schedules' errors
- `every_fitted_number` lays out every error too (`every_error`, a negative code for each reason) and requires a
  fitted stratum to carry them. Returning no errors in the several-at-once schedule now fails three tests.

### M1, M2, M3 — the reasons, and the rule for a number not placed
- `StratumError { Estimated, NotIdentified, NotPlaced, NoShare }` replaces `Option<f64>`; `value()` gives the
  option. A coordinate whose error on the climb's scale exceeds `NOT_PLACED` = 3 is `NotPlaced`, the threshold
  the correctness reviewer recommended from the two populations it measured (0.02 to 0.5 against 1,000 to
  2,400); applied to the slippage numbers, the concentration, and each non-largest class's log-ratio.
- **Also found while applying it**: the oracle's summary then showed a length-class share with an error of 2.85 —
  the largest class's, carried from other classes' log-ratios that were not placed. A number in `[0, 1]` whose
  error exceeds that range is now `NotPlaced` too, spec §3.2's rule as written; tested in
  `the_errors_are_carried_to_each_numbers_own_scale`.
- Test: `a_number_the_climb_left_at_the_end_of_its_reach_is_not_placed`, on the reviewer's three strata — a
  fall-off at 1.3 × 10⁻⁷, a shorter share at 0.9999997, a concentration at 8.0 × 10⁶, each `NotPlaced`.

### M4 — the oracle's time
- The claim is removed; the report says four runs of the unchanged code took 1 min 24 s to 1 min 42 s, and gives
  the drawn strata's and the reviewer's measurements.

### M5 — blended numbers
- **Adaptation:** the field's doc says the errors are the own fit's, not the blend's, and that the own fit's level
  is the one `LevelProvenance::slipped_reads` is counted from. Keeping the own fit's numbers beside them waits for
  E2, which writes the error into the origin blocks that already describe the own fit (spec §5.2).

### M6 — tests
- `tracts_without_reads_change_no_error`, `the_curvatures_centre_is_the_fitted_answer`,
  `the_summary_gives_each_kinds_errors_and_counts_the_missing_by_why`. The tract count over tracts with reads,
  the spectrum not renormalised, and the summary's level not a percentage each now fail one
  (`tmp/fit_precision/c1/mut/`).

### Mi1 — the step
- **Adaptation:** one step, 10⁻², on every coordinate, the scales being relative already; the doc no longer
  claims a ten-thousandth.

### Mi4 — the log line
- Percentages for the level and the concentration, the class shares summarised, the missing counted by why.

### Mi6, Mi7
- The doc names the fixtures; the `expect` on the spectrum's largest class carries a PANIC-FREE note. The empty
  `starting_points` panic is older than C1 and left.

### Deferred
- Mi2: whether the 10⁻⁸ floor suits a finite-difference curvature — C2's coverage will show a floor that keeps
  noise; Mi5: moving the code; Mi8: starting a stratum's errors as its walks finish; the off-diagonal test.

## 5. Deferred findings to carry forward
- Mi2, Mi5, Mi8, the off-diagonal test, the Nits.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — no bench covers the stratum errors.

## 10. Commands run
- `scripts/dev.sh cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::ssr_fit`
- `scripts/dev.sh cargo test --all-targets --all-features --no-fail-fast`
- four mutations (`tmp/fit_precision/c1/mut/run.sh`); the oracle (`tmp/fit_precision/c1/run_oracle_c1f.sh`)

## 11. Command results
In the implementation report, [fit_precision_c1_2026-09-30.md](../implementations/fit_precision_c1_2026-09-30.md) §4.

## 12. Notes
None.
