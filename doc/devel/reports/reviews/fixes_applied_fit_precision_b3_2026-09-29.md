# Fix Application Report: fit_precision_b3_2026-09-29.md

**Date:** 2026-09-29
**Source review:** `doc/devel/reports/reviews/fit_precision_b3_2026-09-29.md`
**Source state reviewed against:** `fc0a89a0` (a review object of the working tree, parent `0b9459a0`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 3
- Minors: 8
- Nits: grouped

### Outcome totals
- Applied: 8 (M1, M2, M3, Mi1, Mi2, Mi3, Mi4, Mi7)
- Applied with adaptation: 1 (Mi8)
- Deferred: 1 (Mi6)
- Disputed: 1 (Mi5)
- Nits: applied in part

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `cargo test --release --lib -- parameter_estimation::joint::fit cli::cross_platform_digests` → 0, 113 passed,
  1 ignored; both checksums unchanged
- `cargo test --all-targets --all-features --no-fail-fast` → 101, 4,984 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1), 5 ignored
- `cargo doc --no-deps`, `cargo audit` → not run (no dependency change; the new public types are documented)
- Performance check → not applicable (logging and a record; the reason's cost was measured in review, 5.45 ms a
  start at kimura's size)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Initial decision | Final status | Files changed | Validation |
|---|---|---|---|---|---|---|
| M1 | Major | trace keys and positions in the log | Apply | Applied | `fit.rs`, `fit/standard_errors.rs`, `fit/settled.rs` | Pass |
| M2 | Major | loop test blind to "stop at one" and to the verdict's pass | Apply | Applied | `fit/information.rs` | Pass; mutations below |
| M3 | Major | the record's reason not checked | Apply | Applied | `fit/information.rs` | Pass; mutations below |
| Mi1 | Minor | `names[at]` unasserted | Apply | Applied | `fit.rs` | Pass |
| Mi2 | Minor | no reason printed when none can be computed | Apply | Applied | `fit/settled.rs`, `fit.rs` | Pass |
| Mi3 | Minor | a never-judged, settled start at the limit | Apply | Applied | `fit/settled.rs` | Pass |
| Mi4 | Minor | verdict attached by `last_mut` | Apply | Applied | `fit.rs` | Pass |
| Mi5 | Minor | `parameter` as a `String` for E2 | Dispute | Disputed | — | N/A |
| Mi6 | Minor | repeated inversion | Defer | Deferred | — | N/A |
| Mi7 | Minor | the cost figure | Apply | Applied | impl report | N/A |
| Mi8 | Minor | the line's length and "errors" | Apply | Applied with adaptation | `fit/settled.rs`, `fit.rs` | Pass |
| Nits | Nit | | Apply in part | Applied in part / Deferred | | |

## 3. Questions asked and answers

None.

## 4. Per-finding log

### M1 — names
- **Implementation:** `Parameters::plain_names` — the cohort's in `COHORT_PARAMETER_NAMES` (now `pub(super)`),
  the standard-error line's words; "homozygote excess of <sample>"; "error rate at ordinary positions of
  <sample>" or "of <sample>'s library k" where the sample has several. `WhichStart` carries the sample names;
  the start's reason and the per-cycle progress line both use them, and say "standard errors" and "where the
  likelihood peaks". Test: `the_log_names_each_parameter_by_its_sample`.

### M2 — the loop test
- **Implementation:** the fixture moved to 20 samples at seed `+3`, whose judged passes find 3, 1 and 0
  unsettled at passes 10, 13 and 16; the test requires a judged pass with exactly one, no verdict on the first
  pass, verdicts at least three passes apart, and the converged cycle finished after its judging pass.
- **Verification:** stopping at ≤ 1 unsettled and writing every verdict on pass 1 each fail it
  (`tmp/fit_precision/b3/mut/`).

### M3 — the record's reason
- **Implementation:** `a_start_at_the_limit_reports_the_settled_test_at_the_values_it_returned` — at a settled
  fraction of 0.15, recomputes errors, Newton distances, the count, the furthest and its distance from a fresh
  pass (posteriors kept, as the final pass keeps them) at the returned parameters; they match, the distance
  bit for bit on this cohort (0.23443559761649738).
- **Verification:** the first parameter's name, the distance doubled, a fraction of 0 and a total 5 too high
  each fail it. The reason read from the last judged cycle would move the distance to 0.23497 (the tests
  reviewer's measurement), past the test's 10⁻⁹.

### Mi1, Mi4
- `how_far_from_settled` asserts one name a parameter; the verdict is written on the trace row whose pass is
  the judging pass's number.

### Mi2, Mi3, Mi8 — the clause
- `settled::describe_short_of_settled` writes the clause, tested for all four cases: converged (nothing), some
  unsettled, every parameter within the fraction (a start stopped before any cycle judged it), and nothing to
  judge. **Adaptation (Mi8):** the line keeps the timing clause after it; the wording is shorter per item and
  says "standard errors".

### Mi5 — a `String` for step E2 (disputed)
- Spec §5.2 has the file record each start's outcome — converged, at the limit or agreed, with its passes — not
  which parameter held it; E2 does not parse this field. It is now a plain name for a reader.

### Mi6 — the repeated inversion (deferred)
- 5.45 ms a start at kimura's size against 524 ms a pass (review); at most three starts a fit.

### Mi7
- The report now gives 5.45 ms for the whole computation, and says whose measurement it is.

### Nits
- A NaN distance would print as "NaN standard errors": deferred, as no fixture produces one. `named()` building
  every name, and the types' home in `settled`: deferred. The report's header links now resolve.

## 5. Deferred findings to carry forward
- Mi6; the re-exported types' home (Nit).

## 6. Disputed findings to return to reviewer
- Mi5 — above.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — logging and a record; the reason's cost measured in review.

## 10. Commands run
- `scripts/dev.sh cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --release --lib -- parameter_estimation::joint::fit cli::cross_platform_digests`
- `scripts/dev.sh cargo test --all-targets --all-features --no-fail-fast`
- six mutations (`tmp/fit_precision/b3/mut/run.sh`); the oracle (`tmp/fit_precision/b3/run_oracle_b3f.sh`)

## 11. Command results
In the implementation report, [fit_precision_b3_2026-09-29.md](../implementations/fit_precision_b3_2026-09-29.md) §3.

## 12. Notes
None.
