# Code Review: fit_precision_a8
**Date:** 2026-09-28 (finished 2026-09-29)
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A8 of `fit_precision.md` — the whole information matrix at small cohorts
**Status:** Request-changes (no Blocker; one design limit for the owner)

---

## 1. Scope

- **What was reviewed:** the diff `11da7a62..c6ec817d` (design, reliability) and `11da7a62..e3641146`
  (correctness and cost; `e3641146` is `c6ec817d` with the first two reviews' fixes). Both are
  review-only commits on no branch.
- **In-scope files:** [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
  (`FULL_MATRIX_SAMPLES`, `FullInformation`, `InformationSums::for_a_cohort_of`, the coverage test),
  [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)
  (`StandardErrors::of`, `of_the_whole_matrix`, `of_the_blocks`), one line of
  [fit.rs](../../../../src/parameter_estimation/joint/fit.rs), the draft implementation report.
- **Categories dispatched**, three reviewers, each in its own worktree:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan step A8; every figure (1 mutation, a probe of the identification order);
  - *reliability and errors* — tests and mutation (24 mutations, three review tests);
  - *correctness, float portability, cost* — an independent derivation of the layout, bit identity at
    1, 4 and 8 threads, no fitted number moved against the parent, the final pass's time and memory
    at 20 samples of 1, 4 and 16 libraries (1 mutation).
  A first dispatch of the first two was cut off by a restart of the session before either wrote
  findings; both were dispatched again, reading the first attempt's saved logs.

## 2. Verdict

**Request-changes, for tests and docs; the code is right.** The correctness reviewer derived the
upper-triangle layout, `row_start`, `entry` and the rows of each sample for samples of one to four
libraries, and found that the whole-matrix errors are what spec §3.2's outer-product estimator means.
**No fitted number moves**: the checksum tests pass unchanged, and three drawn fits taking the
whole-matrix path give the same digests of everything but the errors on this step and its parent.
**The same bits at 1, 4 and 8 threads**, whole matrix included, on 20 samples of one to four libraries
over 128 chunks, with and without the duplicated class. Only `+` and `×` in the new code, plus the
existing inversion.

## 3. Execution status

- Each reviewer: `cargo test --release --lib parameter_estimation::joint::fit` — 76 passed at
  `c6ec817d`, 79 at `e3641146`, 1 ignored (plus each reviewer's probes).
- Mutations: reliability 24 (21 caught by committed tests, 2 only by its probes, 1 no behaviour
  change); design 1 (survived all committed tests, caught by its probe); correctness 1 (caught).
- The coverage test at 3 cohorts a regime (reliability); its 200-cohort run is the author's.
- Not run: the full suite (the author's).

## 4. Open questions and assumptions

- **Whether to cap the whole matrix by parameters as well as samples** (Major, below) is the owner's,
  since the plan sets the limit on samples.

## 5. Findings

### Major

**M1: [information.rs](../../../../src/parameter_estimation/joint/fit/information.rs) `FULL_MATRIX_SAMPLES`
— the limit is on samples, not on the matrix's size.** **Categories:** design, correctness and cost
(convergent). Samples of many libraries make the matrix grow as the square of the parameters, and the
final pass holds 129 copies: 20 samples of 4 libraries, 188 parameters, 18.3 MB; of 16, 668
parameters, 231 MB (the pass's peak measured 261 MB higher, the blocks included), 1.33 to 1.64 times
the blocks' pass; of 50, 2,028 parameters, 2.1 GB. **Recommendation (both):** also cap the parameters;
the correctness reviewer proposes 188.

**M2: [standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)
`of_the_whole_matrix` — nothing tests the order the parameters are inverted in.** **Categories:**
design (Minor), reliability (Major). Taking the cohort's parameters first survives every committed
test: the only drop test makes one cohort parameter mimic another, which drops the later slot under
either order. The order decides which is dropped when a cohort parameter mimics a sample's.

**M3: `of_the_whole_matrix` — nothing checks that each sample slot is judged against its own
interval.** **Category:** reliability. Judging every slot by the clean rate's interval survives, since
every fixture's errors are near 10⁻³; a mismapped rate's error of 0.3 would then be reported wider
than its range.

### Minor

- **m1** (design) — the module doc of `standard_errors.rs` contradicts the code: two samples'
  parameters "are never paired"; the order is "the cohort's first"; at two samples the matrix "does
  not invert at all", where the whole-matrix path now drops those parameters as not identified. The
  drops follow the blocks' rule only on an arrow.
- **m2** (design) — spec §3.2 still reads "a block approximation, not the full matrix". The spec is
  the owner's.
- **m3** (design, reliability) — `InformationSums::absorb` and `FullInformation::absorb` check only in
  debug builds that two chunks agree on the whole matrix's presence and size; a release build would
  silently report errors from part of the positions. `Statistics::absorb` uses `unreachable!`.
- **m4** (design) — the inverted parameters are a bare `(Option<usize>, usize, usize)` tuple.
- **m5** (design) — the coverage test's "the whole matrix beats the blocks by 0.08" compares shares
  over different sets of estimates; its bounds were widened (0.78 to 0.80) before the 200-cohort run.
- **m6** (design) — `whole_matrix_by_hand` was inserted between `entries_with_scale` and its doc.
- **m7** (reliability) — the coverage test's new assertion checks each reported error's size, not the
  reason for a missing one.
- **m8** (reliability) — `a_sample_without_reads_carries_no_information` builds blocks only; nothing
  runs a sample without reads through the whole matrix.
- **m9** (reliability, cross-category) — two test docs still call the blocks "what the fit reports".
- **m10** (correctness) — the draft report's §2 and §5 are stale.

### Nit

- **n1** (design) — "full" and "whole" name the same thing; `FULL_MATRIX_SAMPLES` does not say it is a
  limit; `sample_parameters` returns matrix rows; a scratch row inside a summed type forces a
  hand-written `Debug`; the `blocks` tally's doc overstates what it is.

## 6. Figures checked

All correct: 68 parameters and 2,346 products, 724 for the blocks, 228 and 26,106, 26 parameters in the
dense test, 1.23–1.67 and 0.89–1.02 from the A4 report. Corrected: 129 matrices at the peak, not 128;
the A4 scatter figures apply to the error rates and the mismapped share, not to every estimate; the
gap at 20 samples is 10.6%, the mismapped rates at 3 reads (1.067 against 0.965).
