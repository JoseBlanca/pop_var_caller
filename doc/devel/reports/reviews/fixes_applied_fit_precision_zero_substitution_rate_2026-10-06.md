# Fix Application Report: fit_precision_zero_substitution_rate_2026-10-06.md

**Date:** 2026-10-06
**Source review:** `doc/devel/reports/reviews/fit_precision_zero_substitution_rate_2026-10-06.md`
**Source state reviewed against:** review object `d696f374` (parent `479e5ddb`)
**Execution mode:** non-interactive
**Overall status:** Partial — M1 awaits the owner

---

| ID | Final status | What was done |
|---|---|---|
| M1 | Ask | Put to the owner: refuse a zero or one in a file, or convert a zero with a count on reading. Not changed here; the behaviour for a version-1 file is the one it had before this change. |
| Mi1 | Applied | The doc says what a zero does — a read with one disagreeing base is explained by no tract length and stops counting for any genotype — and its examples say "no mismatch in 500 bases", "no mismatch in 34". |
| Mi2 | Applied in code; spec Ask | `StratumEvidence::substitution_rate`'s doc points to the rule. Spec `parameter_prepass_ssr.md` §4.2 not edited without the owner. |
| Mi3 | Applied | The full branch is `== bases_compared`, with a `debug_assert!` that no more are mismatching than compared; an impossible count gives a ratio above one, which no rate type accepts. |
| Mi4 | Applied | The baseline note says the 52 rows gain the error they lacked. |
| Mi5 | Applied | `a_count_with_one_outcome_only_takes_half_a_count_of_the_other` holds the rule; the run-level test renamed `every_stratum_with_bases_compared_gets_a_rate_and_its_binomial_error`. |
| Nit | Applied | `bases_compared`, `bases_compared_f64`. |

Validation: `cargo fmt --check` and clippy `-D warnings` clean; full suite in the implementation report. No number moved
since the review object: the fixes are docs, names, tests and the unreachable `>` case.
