# Fix Application Report: fit_precision_d3_2026-10-02.md

**Date:** 2026-10-02
**Source review:** `doc/devel/reports/reviews/fit_precision_d3_2026-10-02.md`
**Source state reviewed against:** `836c696d` (a review object of the working tree, parent `349c0532`)
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

- Review: 0 Blockers, 2 Majors, 3 Minors, Nits grouped. Applied: all.
- Validation: `cargo fmt --check` → 0; `cargo clippy --all-targets --all-features -- -D warnings` → 0; the tool run
  on the four-sample cohort with a two-sample subset (`tmp/fit_precision/d3_smoke4.log`). The library is unchanged
  since the reviewed code's full suite (5,025 passed, the 3 pre-existing failures) and oracle (every checksum
  matches).

## 2. Findings table

| ID | Severity | Title | Final status | Files changed |
|---|---|---|---|---|
| M1 | Major | "apart" when the subset is every sample | Applied | the example |
| M2 | Major | default strata; silent long fits | Applied | the example |
| Mi1 | Minor | `--inbreeding` unchecked | Applied | the example |
| Mi2 | Minor | the reading of "apart" | Applied | the example |
| Mi3 | Minor | no time total | Applied | the example |
| Nits | Nit | | Applied | the example |

## 3. Questions asked and answers

None.

## 4. Per-finding log

- **M1:** when the subset grew to every sample, the difference is printed in the whole fit's errors and labelled so.
- **M2:** the default strata are five spread across the cohort's sizes; the subset fit runs first and prints at
  once; a line before the fit on every sample gives its likelihood table; the doc states what kimura's fits cost.
- **Mi1:** `--inbreeding` checked with `InbreedingF::try_new`.
- **Mi2:** the doc says about one in twenty numbers is expected beyond two by chance, that the measure takes the
  errors at their word (so `--inbreeding` on an inbred cohort), what "-" means, and which lines answer spec §4.5
  item 3.
- **Mi3:** a last line totals both fits' time.
- **Nits:** a refused subset fit prints its reason without dumping the other fit; the fit on every sample runs only
  after the subset fit succeeds; the doc sentence rewritten; "every sample, no subset drawn" when none was.

## 5–8. Deferred, disputed, failed, blocked
None.

## 9. Performance check
Not applicable: an example.
