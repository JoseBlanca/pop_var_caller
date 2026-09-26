# Fix Application Report: estimation_memory_issue1_2026-09-26.md

**Date:** 2026-09-26
**Source review:** `doc/devel/reports/reviews/estimation_memory_issue1_2026-09-26.md`
**Source state reviewed against:** `a75855a9`, branch `evidence-release`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

### Review totals
- Blockers: 0
- Majors: 0
- Minors: 5
- Nits: 6

### Outcome totals
- Applied: 5 Minors, 4 Nits
- Deferred: 1 Nit (one declaration of the three counts)
- Won't fix: 1 Nit (`drop(evidence)` kept for the intent)

### Validation summary
- `cargo fmt --check` → clean
- `cargo clippy --all-targets --all-features -- -D warnings` → clean
- `cargo test --all-targets --all-features --no-fail-fast` → see §10
- Performance check → not applicable (no hot-path change)

### Unresolved high-priority findings
- None.

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| Mi1 | Minor | handover of counts untested | Applied | `src/run/census_fit.rs` |
| Mi2 | Minor | nothing-compared and nonzero rate untested end to end | Applied | `src/run/census_fit.rs` |
| Mi3 | Minor | loop variable named `stratum` | Applied | `src/run/census_fit.rs` |
| Mi4 | Minor | field doc says "read off the evidence" | Applied | `src/run/census_fit.rs` |
| Mi5 | Minor | report names the wrong status block | Applied | the implementation report |
| Nit | Nit | type doc sentence and unsourced figure | Applied | `ssr_fit.rs` |
| Nit | Nit | `stratum` field undocumented | Applied | `ssr_fit.rs` |
| Nit | Nit | tautological assertion | Applied (deleted) | `ssr_fit.rs` |
| Nit | Nit | closure argument name | Applied | `census_fit.rs` |
| Nit | Nit | counts declared on both types | Deferred | — |
| Nit | Nit | `drop(evidence)` | Won't fix: shows the intent at no cost | — |

## 3. Questions asked and answers

None.

## 4. Per-finding log

- **Mi1.** `a_cohort_of_censuses_is_fitted_both_halves` now asserts the counts' strata equal the
  outcomes' strata in order, that at least one has bases compared, and that the second fit's counts
  equal the first's.
- **Mi2.** New test `a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own`.
  The fixture builder `a_fitted_cohorts_parameters` became `a_fitted_cohorts_parameters_with(extra)`,
  which appends `extra` to the fit's counts before assembly; the old name calls it with none. The
  test counts the added stratum's rates against the fixture's own stratum's, so it does not depend
  on how many read groups the fixture has. Adapted from the review's sketch: `key.stratum.period`
  is a field, not a method.
- **Mi3, Mi4.** Renamed, reworded.
- **Mi5.** Report corrected in place, with a note saying so.
- **Nits.** As in the table.

Validation after all fixes: `cargo test --lib -- census_fit substitution_counts cross_platform` →
12 passed, 0 failed (the new test is the twelfth).

## 5. Deferred findings to carry forward
- Nit: `StratumEvidence` could hold a `StratumSubstitutionCounts` rather than repeat its fields —
  reaches every place the evidence is built.

## 6–9.
None.

## 10. Commands run
- `./scripts/dev.sh cargo test --lib -- census_fit substitution_counts cross_platform` → 0, 12 passed
- full validation: recorded in the commit's report line (below)

## 12. Notes
None.
