# Fix Application Report: estimation_memory_issue3_2026-09-26.md

**Date:** 2026-09-26
**Source review:** `doc/devel/reports/reviews/estimation_memory_issue3_2026-09-26.md`
**Source state reviewed against:** `dbdba774`, branch `contamination-two-pass`
**Execution mode:** non-interactive
**Overall status:** Completed

---

## 1. Executive summary

- Review totals: 0 Blockers, 1 Major, 5 Minors, 5 Nits.
- Applied: the Major, 4 Minors, 3 Nits. Deferred: Mi5 and the two performance savings. Won't fix:
  2 Nits (carried-over code).
- Unresolved high-priority findings: none.

## 2. Findings table

| ID | Final status | What was done |
|---|---|---|
| M1 | Applied | `two_passes_agree_on_deep_reads_and_a_second_alternative_allele`: 30 samples at 150 reads a position, with every third sample's alternative reads at every fourth position recorded as `G` (a new fixture variant, `structured_panel_with_a_second_allele`, which the existing panel now calls with no `G`). It asserts some position's common allele is `G`, some `G` reads are left out, and some depth code is clamped by the cap, then compares. `a_dosage_whose_weights_all_underflow_is_twice_the_panel_frequency` pins the fallback |
| Mi1 | Applied | Report corrected: per-library lists and dosages both counted, at 52,525 and at the current 31,758 markers |
| Mi2 | Applied | Report states the full-suite figures |
| Mi3 | Applied | `the_same_markers` destructures `Marker` |
| Mi4 | Applied | `two_passes_agree_on_a_panel_below_the_coverage_floor`: 1 and 5 samples, no markers either way |
| Mi5 | Deferred | A scratch type owning the pair is a reasonable refactor, but each pair has one call site per pass, both covered by the differential tests; not needed for this step |
| Nits | 3 applied (renamed `fill_the_reads_of_each_library` with parameter `markers`; the no-op `allow` removed; the 64 MB deviation recorded in the report), 2 won't fix (`+=` and `.expect` carried over verbatim; the test helper's `..Default::default()` is the file's pattern) | |
| Perf | Deferred | Plan §6 fixes the shape (the per-library pass unchanged); folding the dosages into it is a follow-up, recorded in the report |

## 3. Questions asked and answers

None.

## 4. Verification of the new test against the defects it exists for

Each defect put in by hand, the file restored byte for byte afterwards (`cmp`):

- the common-allele filter in `pool_one_sample` made always true → the deep test fails ("how many
  markers"), the other three differential tests pass;
- the depth cap's clamp dropped in `pool_one_sample` → the deep test fails the same way.

## 10. Commands run

- `./scripts/dev.sh cargo test --lib -- contamination cross_platform` → 40 passed, 0 failed.
- `./scripts/dev.sh cargo test --lib -- two_passes` under each defect above → 1 failed, 3 passed.
