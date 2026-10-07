# Code Review: fit_precision_errors_pinned_cross_platform
**Date:** 2026-10-06
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** a second cross-platform digest over the fixture with sequencing errors, so the fit's standard errors are pinned
**Status:** Request-changes

---

## 1. Scope

Review object `755d29b6` against `2ffca34f`: `src/cli/test_fixtures.rs` (the fixture variant) and
`src/cli/cross_platform_digests.rs` (the new test and the thread-count test over both fixtures). One reviewer, in its
own worktree: reliability, float portability, naming, smells, and the prose of the doc comments. Test-only change.

## 2. Verdict

Request-changes: one Major.

## 3. Execution status

All three digest tests pass at `755d29b6` (Linux container, release). Mutations: 9 considered, 8 built, 1 simulated —
4 killed, 3 survived with a proven internal change, 2 changed no behaviour.

## 6. Findings

### Major

**M1: [cross_platform_digests.rs:166](../../../../src/cli/cross_platform_digests.rs#L166) — the errors are pinned only
at 20 samples or fewer.** Confidence High. Above `FULL_MATRIX_SAMPLES` (20) the errors come from each sample's own block
of the information, a path a two-sample cohort never reaches: re-rounding the block sums changed 55,865 of them and moved
neither checksum. That path is every large cohort's. **Fix:** a cohort of 21 samples or more with its own checksums, or a
test-only way to force the blocks.

### Minor

- **Mi1** `:170` — the cohort's score for the mismapped share does not reach the file (computed with std `exp`: 1,668
  calls differed; doubled: no byte moved); the doc implies the whole inversion is pinned.
- **Mi2** `test_fixtures.rs:262` — "never at a variant site" is false: an error lands on the other sample's site, and
  two inside the tract.
- **Mi3** `test_fixtures.rs:367` — "no reference position collects the same wrong base" is false: two positions do.
- **Mi4** `:170` — "near 7" is off for one library (5.59).
- **Mi5** `:348` — the pool-width comment points at a join-order defect three chunks cannot show.
- **Mi6** `:315` — "the platforms this project runs on": no x86_64 measurement.

### Nits

The hash seeded from the name's first byte only; the substitution guard compares reference coordinates that the tract's
deletion shifts (correct here by placement); unnamed hash constants; a `Vec` a read for three bases; a boolean
parameter.

## 8. Missing tests

`a_fit_with_standard_errors_through_the_blocks_writes_the_same_bytes_on_every_platform` (M1).

## 9. What's good

Each new checksum is guarded by an assertion that the file carries what it exists to pin — every multiplier with an
error — so a fixture that stopped producing errors fails loudly instead of pinning nothing.
