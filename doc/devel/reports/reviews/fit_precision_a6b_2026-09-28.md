# Code Review: fit_precision_a6b
**Date:** 2026-09-28
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A6, second part, of `fit_precision.md` — each library's rates get their own standard errors
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the diff `373a2b7e..09f58fd6` (`09f58fd6`, "WIP: fit precision A6 second
  part, for review", a review-only commit on no branch).
- **In-scope files:** [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
  (the row layout, `PositionScores`, `fill_read_slopes`, the branch scorers, `InformationSums`),
  [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)
  (rows of each sample's length, `bounds_of_own`, `named`, `described`),
  [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) (wiring, the whole-fit test), the
  A6b implementation report.
- **Categories dispatched**, three reviewers, each in its own worktree:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan A6 and spec §3.2; the prose (8 mutations);
  - *reliability and errors* — tests and mutation (22 mutations, probes with three-library and
    no-library samples);
  - *correctness, float portability, cost* — an independent derivation, bit identity against the
    parent, coverage over drawn two-library cohorts, the final pass's cost (16 mutations). Its first
    dispatch could not run shell commands (a tool error) and reviewed nothing; it was dispatched
    again (`tmp/review_2026-09-28_fit_precision_a6b/correctness_cost_attempt1_not_run.md`).

## 2. Verdict

**Request-changes**, for one defect in the file itself and for tests that could not fail. The maths
is right: the correctness reviewer derived each library's score and the inversion over blocks of
different sizes independently. **One-library errors are the same bits as before** on 10 of 10 fits
(4 samples and 1, with and without the duplicated class, one at 25 reads; 1 and 4 threads), and a
two-library cohort gives the same bytes at 1, 4 and 8 threads. **The errors of two-library samples mean
what they say**, measured over 240 drawn cohorts:

| cohort | clean rates within 1 / 2 errors | mismapped within 1 / 2 | spread in errors (clean / mismapped) |
|---|---|---|---|
| 4 samples × 2 libraries, 3 and 6 reads | 65.4% / 93.3% | 66.7% / 93.3% | 1.05 / 1.06 |
| 8 samples × 2 libraries, 3 and 6 reads | 67.5% / 95.6% | 67.9% / 93.4% | 0.98 / 1.04 |
| 4 samples × 2 libraries, 15 and 25 reads | 65.8% / 94.2% | 69.2% / 96.0% | 1.03 / 0.96 |
| control: 4 samples × 1 library, 9 reads | 65.0% / 95.8% | 60.4% / 92.9% | 0.99 / 1.14 |

(20,000 positions a cohort, 60 cohorts a regime, all converged, no error absent; nominal 68.3% and
95.4%.)

## 3. Execution status

- Each reviewer: `cargo test --release --lib parameter_estimation::joint::fit` — 65 passed, 1 ignored.
- Mutations: design 8 run (2 survived with a behaviour change); reliability 22 run (8 survived with a
  behaviour change, 0 no change); correctness 16 run (2 survived with a behaviour change).
- Cost of the final (information-collecting) pass, relative to a plain pass in the same process, 4
  rounds × 5 repetitions, 60,000 positions, one thread: one library a sample (16 samples) from
  1.53–1.56 to 1.69–1.71 times; two libraries a sample (8 × 2) from 1.37–1.45 to 1.71–1.79 times (the
  parent computed no slopes for them).
- Not re-run by reviewers: the oracle, the full suite.

## 4. Open questions and assumptions

None needing the owner.

## 5. Top 3 priorities

1. **M1** — the double-encoded file.
2. **M2, M3, M4** — tests that cannot fail on a sibling's error, a third library, the slots' intervals.
3. **M5, M6** — a sample holding no library; the row-skip rule for several libraries.

## 6. Findings

### Major

**M1: [information.rs](../../../../src/parameter_estimation/joint/fit/information.rs) — the whole
file is double-encoded.** **Categories:** all three (convergent). Every non-ASCII character (—, §,
×, ⁻, ψ, Σ) is stored as UTF-8 of its UTF-8 bytes read as Latin-1: 117 lines, including three test
messages; about 230 of the diff's lines only re-encode old text. Line 1644, added after the damage, is
correct. Cause: a perl edit whose replacement held a wide character, which made perl encode the whole
buffer on output. **Fix:** decode once, re-encoding only the characters above U+00FF.

**M2: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs)
`a_samples_two_libraries_come_back_at_their_own_rates` — cannot tell a library's error from its
sibling's.** **Categories:** design, reliability, correctness. Reading every library's error from the
first library's slot still passes (largest distance 2.22 errors, under the bound of 3); the largest
real distance is 1.91, so any error at least 0.64 times the reported one passes; the mismapped rates
are not checked at all. **Fix:** assert something only each library's own error satisfies.

**M3: [information.rs](../../../../src/parameter_estimation/joint/fit/information.rs) — no pass runs a
sample of three libraries.** **Category:** reliability. Dropping the third library from the
segregating branch (a slope 9.6% off its finite difference) passes the suite. **Fix:** a
finite-difference test over a sample of three libraries.

**M4: [standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)
`bounds_of_own`, [information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
`sample::class_of` — which interval each slot is judged against is untested.** **Categories:**
design (Minor there), reliability, correctness (Minor there). Judging the excess as a clean rate, swapping
a further library's classes, or capping the slot at the call site all pass; an error between 0.2 and
0.45 flips between an error and "wider than its range". `class_of` writes the one-library row size
as a bare `3`. **Fix:** a layout test, and a check through the inversion.

**M5: a sample holding no library is never exercised, and is reachable.** **Category:** reliability.
The census keeps every sample, so a sample with repeat-tract sections only gets an empty library list.
On the reviewed code it behaves as the report says, but a mutation that panics and one that counts a
read group that does not exist in the log both survive. **Fix:** a test.

**M6: `InformationSums::add_position`'s row-skip rule is untested for several libraries.**
**Category:** reliability. Skipping a row when any score is zero instead of all passes the suite and
widens a mismapped-rate error from 0.0160 to 0.0229 with 7 of 600 positions holding an empty library.
**Fix:** a by-hand comparison of a pass's sums on a cohort of multi-library samples.

### Minor

- **Mi1: `TOLERANCE`'s doc** says the largest disagreement any test reaches is 1.7 × 10⁻⁷; it is now
  6.76 × 10⁻⁷ (a homozygote excess in the two-library fixture).
- **Mi2: rows of different lengths are zipped with no check** in `add_position`, `absorb`, `named`
  and `described`.
- **Mi3: `JointFit::standard_errors`' doc** still says "each sample's two error rates" before
  stating the new layout.
- **Mi4: the ranged test uses `cohort::DUPLICATED_SHARE`** as the count of cohort coordinates, where
  `cohort_coordinates` exists.
- **Mi5: the report's §2 cross-reference** to "§5, the oracle" should be §6; and it states no cost
  for one library a sample, where the final pass is about 10% longer (§3).

### Nits

The report's 5.7 × 10⁻⁵ at a step of 10⁻⁵ is in `tmp/fit_precision/a6b_tests.log`, not the
`a6b_fd_*.log` files the report names.

## 7. Out of scope observations

None.

## 8. Missing tests to add now

- `a_sample_of_three_libraries_has_every_slope_right` (M3).
- `a_sample_holding_no_library_carries_no_information` (M5).
- `each_slot_of_a_samples_row_is_its_own_parameter_with_its_own_bounds` and a width check through the
  inversion (M4).
- A mixed-library cohort in `a_pass_sums_the_scores_multiplied_pairwise` (M6).
- An own-error assertion in the two-library whole-fit test (M2).

## 9. What's good

- The row layout that keeps a one-library sample's three slots made the one-library errors
  bit-identical without special cases.
- The arrow inversion needed only its block size read from the sums: it was already written over the
  parameters a block keeps.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`
- Per-category files: `tmp/review_2026-09-28_fit_precision_a6b/{design_claims,reliability_errors,correctness_cost}.md`;
  probes under `…/correctness_evidence/` and `…/evidence_*`.
