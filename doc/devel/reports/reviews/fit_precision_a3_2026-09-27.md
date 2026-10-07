# Code Review: fit_precision_a3
**Date:** 2026-09-27
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A3 of `fit_precision.md` — standard errors from the information blocks
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the diff `0988d9ac..9943a969` (`9943a969`, a review-only commit on no branch).
- **In-scope files:** [src/parameter_estimation/joint/fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)
  (new), [src/parameter_estimation/joint/fit.rs](../../../../src/parameter_estimation/joint/fit.rs)
  (`fits_homozygote_excess`), the test-only comparison in
  [src/parameter_estimation/joint/fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs),
  the A3 implementation report.
- **Categories dispatched**, three reviewers, each in its own worktree:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan A3 and spec §3.2, §3.6 item 2; the prose (15 mutations);
  - *reliability* — tests and mutation (23 mutations);
  - *correctness and statistics* — the arrow algebra, the rules for an absent error, the
    interpretation of the comparison, float portability; probes at 1 and 2 samples and at 3,000 to
    300,000 positions (3 mutations).

## 2. Verdict

**Request-changes.** The arrow algebra is right (derived independently; every figure reproduced), and
portable. But the rule for a parameter the data cannot determine broke at the smallest cohorts: with
the duplicated class on, one or two samples gave no error at all — the error rates' included, against
spec §3.2 — or, where a rounding pivot passed the threshold, errors in the millions. And the author's
reading of the block-versus-full gap as chance holds only at large cohorts.

## 3. Execution status

- Each reviewer: `cargo test --release --lib parameter_estimation::joint::fit` — 43 passed at `9943a969`.
- Mutations: design 15 run (6 survived: 4 with a behaviour change, 1 proven equivalent, 1 argued
  unreachable); reliability 23 run (12 killed, 7 survived with a behaviour change, 3 no behaviour
  change, 1 run only against the reviewer's added tests); correctness 3 run, 3 killed.
- Findings labelled *Needs verification*: 0.

## 4. Open questions and assumptions

1. **The block approximation at small cohorts and low depth** (M1): the samples' block errors are too
   small there. Which matrix matches the spread of the estimates is A4's measurement; the choice after
   it is spec question 3 at checkpoint A.
2. **The spec's text** (§3.2) shows the sketch the plan compared and did not ship, and describes block
   sizes that were not built (A2's open question). Owner.

## 5. Top 3 priorities

1. **B1** — drop only the parameter the data cannot determine, and raise the threshold.
2. **M1** — correct the reading of the block-versus-full gap; make A4 measure where it is real.
3. **M2** — test where each error is written back when a parameter in the middle of a block is left out.

## 6. Findings

### Blocker

#### B1: src/parameter_estimation/joint/fit/standard_errors.rs — at one or two samples, no errors at all, or errors in the millions
- **Categories:** correctness
- **Confidence:** High
- **Problem:** with the duplicated class on, at one sample the four frequency-density parameters cannot
  all be told apart (three genotypes, four numbers), and at two the duplicated share and carrier shapes
  cannot. The cohort's reduced block is then exactly singular, one pivot holds rounding, and at 3,000
  and 30,000 positions its inversion fails: every error, the clean error rates' included, came back
  `None`. At 300,000 positions and two samples the rounding pivot kept 3.2 × 10⁻¹² of its diagonal and
  passed the 10⁻¹² threshold: a share got an error of 13.5 on [0, 1], the duplicated share 158, a
  carrier shape 2.0 × 10⁶.
- **Suggested fix:** when a pivot fails, mark only that parameter and refactor without it, inside a
  sample's own block too; a threshold of about 10⁻⁸, judged against the cohort's diagonal before the
  samples' share is taken; tests on drawn cohorts of one and two samples.

### Major

#### M1: the comparison's doc and the A3 report — the gap between blocks and full matrix is chance only at large cohorts
- **Categories:** correctness; design (the 4-sample reading unsupported)
- **Confidence:** High
- **Problem:** two samples at one position share its unknown frequency and class; their expected score
  product is minus the expected posterior covariance of their complete-data scores, not zero. Raising
  the positions separates chance from a real term: at 20 samples and 8 reads the mismapped rates' median
  gap fell from −0.145 to −0.014 (3,000 to 30,000 positions); at 4 samples it went −0.166, −0.135,
  −0.109 (up to 100,000); at 4 samples and 3 reads the clean rates' went −0.331 to −0.186. At 4 samples
  every pair's clean-rate scores correlated at about −0.012. So at small cohorts the blocks' errors are
  too small.
- **Suggested fix:** say so; A4 must measure coverage at two to four samples, at 3 reads, with 30,000
  positions or more.

#### M2: src/parameter_estimation/joint/fit/standard_errors.rs — nothing tests where an error is written back
- **Categories:** reliability
- **Confidence:** High
- **Problem:** every fixture removed parameters only from the end of a block, where a parameter's
  position among the kept ones equals its slot. Writing back by position (M13, M14) survived all tests.
  A read group with no mismapped positions leaves slot 1 of 3 empty — the realistic case.
- **Suggested fix:** tests removing a parameter in the middle of the cohort's and a sample's block.

#### M3: standard_errors.rs — a sample whose own block fails to invert: untested, and its effect not "slight"
- **Categories:** design, reliability (convergent)
- **Confidence:** High
- **Problem:** the whole sample was dropped, and left out of what the samples explain of the cohort's
  information; the cohort's errors are then too small by an amount that grows with the sample's tie to
  the cohort (3.2% on one test matrix), and the sample's independent parameters are lost. No test.

#### M4: standard_errors.rs — the pivot threshold pinned from one side only
- **Categories:** design, reliability (convergent)
- **Confidence:** High
- **Problem:** a threshold of 10⁻² (a pair correlated at 0.999 no longer inverts) and one compared with
  another entry's diagonal survived.

#### M5: standard_errors.rs — no cohort information (`c == 0`) untested
- **Categories:** reliability, design (convergent)
- **Confidence:** High
- **Problem:** a mutant returning no error for every sample survived.

#### M6: information.rs — the comparison test does not catch what its doc says
- **Categories:** design, reliability, correctness (convergent)
- **Confidence:** High
- **Problem:** "a regression that moved a block or dropped a term would move these by far more": with
  the cohort's uncertainty carried into each sample's error dropped, or `B_s` transposed, it passes
  (the 4-sample excess moved from −0.022…+0.089 to −0.302…−0.175, inside ±0.35). Only the unit tests
  catch these.

### Minor

#### Mi1: standard_errors.rs — a non-finite diagonal untested
- **Categories:** reliability · **Confidence:** High

#### Mi2: standard_errors.rs — `inverse_of_positive_definite`'s doc attached to `SINGULAR_PIVOT`
- **Categories:** design, correctness, reliability (convergent) · **Confidence:** High

#### Mi3: information.rs — tolerances keyed on display strings
- **Categories:** design · **Confidence:** High
- **Problem:** renaming a kind silently widens its bound to ±0.35.

#### Mi4: standard_errors.rs — `SampleSolve` / `solved` name an action, not a value
- **Categories:** design · **Confidence:** High

#### Mi5: spec §3.2 still shows the rejected formula
- **Categories:** design · **Confidence:** High

#### Mi6: wrong or unsupported numbers in the prose
- **Categories:** design, reliability, correctness (convergent) · **Confidence:** High
- "+0.001 for each kind" (clean rates +0.000); "about a seventh below at both sizes" (−0.166 at 4
  samples, a sixth); "about ten thousand times a double's rounding" (10⁻¹² is 4,504 times machine
  epsilon); "slightly smaller" (unmeasured); the 4-sample explanation.

### Nits

- `error_from_variance` accepting a variance of exactly 0 has no reachable input.
- Spec §3.6 item 2 asks for 4, 20 and 63 samples; 63 is not run.

## 7. Out of scope observations

- The Beta-shape update of the maximisation (A1 finding) means the density-shape scores are not
  centred on zero at the fitted shapes; A4's coverage will show it if it matters.

## 8. Missing tests to add now

- Drawn cohorts of one and two samples, the duplicated class on: every error rate keeps its error, no
  error is made of rounding (B1).
- A parameter without information in the middle of the cohort's and a sample's block (M2).
- A sample parameter mimicked by another in its own block, dropped alone (M3).
- The threshold from both sides, and a parameter in small units (M4).
- No cohort information (M5).
- A non-finite diagonal (Mi1).

## 9. What's good

- The arrow inverse is exact and cheap — no matrix larger than 8 × 8 at any cohort size — and the unit
  tests check it against a dense inversion of the whole arrow.
- The one-sample rule is one function shared by the maximisation, the parameters' provenance and the
  errors.
- The comparison reports the spec's own sketch beside the shipped formula, per kind, so the choice is
  measured rather than asserted.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit -- --nocapture`
- `scripts/dev.sh cargo fmt --check`; `cargo clippy --all-targets --all-features -- -D warnings`;
  `cargo test --all-targets --all-features --no-fail-fast`

Per-category files: `tmp/review_2026-09-27_fit_precision_a3/` (`design.md`, `reliability.md`,
`correctness.md`).
