# Fix Application Report: fit_precision_c2_2026-09-30.md

**Date:** 2026-09-30
**Source review:** `doc/devel/reports/reviews/fit_precision_c2_2026-09-30.md`
**Source state reviewed against:** `c5010ff2` (a review object of the working tree, parent `33072de2`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 1
- Majors: 4
- Minors: 4
- Nits: grouped

### Outcome totals
- Applied: 4 (B1, M1, M3, M4)
- Applied in part: 2 (M2, the Minors)
- Deferred: the five-class regime and the integral check (after the owner's decision), the Nits

### Validation summary
- `cargo fmt --check` → 0, clean
- `cargo clippy --all-targets --all-features -- -D warnings` → 0, clean
- `a_stratums_errors_mean_what_they_say` at 100 draws → 0, passed, every figure the same as before the fixes
- `cargo test --all-targets --all-features --no-fail-fast` → 101, 4,992 passed, 3 failed (the pre-existing
  `examples/ng_generic_loci_dump.rs` × 2, `examples/ng_ssr_loci_dump.rs` × 1), 6 ignored
- Performance check → not applicable (a test-only step)

### Unresolved high-priority findings
- M3's subject — the 256-point integral — is for the owner (report §4).

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| B1 | Blocker | fails at 50 draws | Applied | `ssr_fit.rs` test |
| M1 | Major | passes when nothing is measured | Applied | `ssr_fit.rs` test |
| M2 | Major | cannot see errors 10–20% off; missing errors unchecked | Applied in part | `ssr_fit.rs` test |
| M3 | Major | §3 overstates the remedy | Applied | impl report |
| M4 | Major | no recommendation | Applied | impl report |
| Minors | Minor | | Applied in part | `ssr_fit.rs` test |

## 3. Questions asked and answers

None; the owner's decision on the integral is asked at the end of the run.

## 4. Per-finding log

### B1
- The bounds apply only at the default count of draws (`STRATUM_DRAWS`, 100), where they were measured, in
  `hold_to_the_three_class_bounds`; the mean-distance bound is 0.4 (the shorter share sits +0.22 high in every run).

### M1
- Draws that give nothing to tally are counted and required to be none; each kind must have one estimate a draw
  (two for the other length shares).

### M2
- **Applied:** every three-class estimate must have an error (the doubled-off-diagonal mutant turned dozens into
  "not identified"); the non-share kinds' spread is held to 0.85 to 1.2. **Not applied:** tighter bounds on the
  shares, whose under-coverage is the integral's (M3). The doc says what the bounds cannot see.

### M3, M4
- Report §3 rewritten: what the fit averages over, in words; the reviewer's table of the log-likelihood at the
  truth and at the fit against the number of points; that no count tried is known to suffice; that the defect
  starts by five classes and reaches the three-class shares; the level. §4 gives the recommendation: a separate
  investigation of the average before step D, continuing with C3 now.

### Minors
- **Applied:** every line printed before any bound is checked; a non-Unicode count refused. **Deferred:** a kind
  enum, named bounds, the `allow`'s reason.

## 5. Deferred findings to carry forward
- A five-class regime and a test of the integral at the truth — after the owner's decision on it.
- The Minors and Nits not applied.

## 6. Disputed findings to return to reviewer
None.

## 7. Failed-validation findings
None.

## 8. Blocked-by-context-mismatch findings
None.

## 9. Performance check
Skipped — a test-only step.

## 10. Commands run
- `scripts/dev.sh cargo fmt`, `cargo clippy --all-targets --all-features -- -D warnings`
- `scripts/dev.sh cargo test --release --lib a_stratums_errors_mean_what_they_say -- --ignored --nocapture`
- `scripts/dev.sh cargo test --all-targets --all-features --no-fail-fast`

## 11. Command results
In the implementation report, [fit_precision_c2_2026-09-30.md](../implementations/fit_precision_c2_2026-09-30.md) §5.

## 12. Notes
None.
