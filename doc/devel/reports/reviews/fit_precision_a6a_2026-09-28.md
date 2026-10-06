# Code Review: fit_precision_a6a
**Date:** 2026-09-28
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step A6, first part, of `fit_precision.md` — each library's reads under its own error rates
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the diff `8df657ca..e529839d` (`e529839d`, "WIP: fit precision A6 first part,
  for review", a review-only commit on no branch).
- **In-scope files:** [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) (the per-position
  reader, `one_position`'s likelihood and tallies, `fit_jointly`, the generator),
  [fit/information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
  (`fill_read_slopes`, tests), [fit/standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs)
  (the interim reason), the A6a implementation report, `PROJECT_STATUS.md`.
- **Categories dispatched**, three reviewers, each in its own worktree:
  - *design and claims* — naming, idiomatic, module structure, refactor safety, smells, defaults;
    conformance to plan A6 and checkpoint A's decision 1; the prose (8 mutations);
  - *reliability and errors* — tests and mutation (33 mutations, with probe tests drawing cohorts of
    several libraries a sample);
  - *correctness, float portability, concurrency, cost* — an independent derivation of the
    per-library likelihood and tallies; bit identity against the parent; the cost of a pass (9
    mutations).

## 2. Verdict

**Request-changes.** The arithmetic is right: the correctness reviewer derived the per-library
likelihood, the tallies and each library's Poisson centre independently and they match; one-library
fits are byte-identical to the parent's on nine drawn cohorts across the three reviewers (1 and 4
threads, duplicated class on, a thirteen-chunk census, depths stored as ranges), and two
multi-library cohorts give the same bytes at 1, 3 and 9 threads. What is wrong is the evidence: the
step's central behaviour — each library's tallies get its own reads — is pinned in one of four
branches, and the interim treatment of the errors understates every other error in a cohort that
holds a sample of several libraries.

## 3. Execution status

- Each reviewer: `cargo test --release --lib parameter_estimation::joint::fit` — 60 passed, 1 ignored,
  at `e529839d`.
- Mutations: design 8 run (5 killed, 2 survived with a behaviour change, 1 no change); reliability 33
  run (14 killed, 15 survived with a behaviour change, 4 no change); correctness 9 run (3 killed, 5
  survived with a behaviour change, 1 no change).
- Cost (correctness reviewer; 8 samples, 40,000 positions, 10 passes, 4 threads, median of 5, three
  interleaved rounds): one library a sample, the review commit 2.3 to 3.6% slower than the parent;
  two libraries at 1.5 reads a position against one at 3 reads, 8.8% slower for the same reads.
- Not re-run by reviewers: the oracle's checksums, the per-pass times.
- Findings labelled *Needs verification*: 0.

## 4. Open questions and assumptions

None needing the owner.

## 5. Top 3 priorities

1. **B1** — a conservation test for every branch's tallies.
2. **M1** — the interim errors of a cohort holding a sample of several libraries.
3. **M2, M3, M4** — tests that can fail on the per-library likelihood, depth handling and counts.

## 6. Findings

### Blocker

**B1: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) `one_position` — three of the four
branches' per-library crediting is untested.** **Categories:** reliability, correctness (convergent).
Crediting a sample's reads to its first library in the fixed-alternative, segregating or duplicated
branch passes all 60 tests; each moves the tallies (on a probe, one library's segregating candidate
count went from 27.5 to 3.6 reads and its sibling's from 34.3 to 58.2; in a whole fit, one clean rate
doubled, 0.00082 to 0.00164, against 0.001 drawn). The whole-fit test cannot see it: 92% of its
positions are invariant and the duplicated class is off. **Fix:** a conservation test on
`one_position` — over one position each library's credited candidate and other-base reads, summed
over classes, branches and genotypes, equal that library's own disagreeing reads (up to the 10⁻¹²
skips) — with the duplicated class on.

### Major

**M1: [standard_errors.rs](../../../../src/parameter_estimation/joint/fit/standard_errors.rs),
[information.rs](../../../../src/parameter_estimation/joint/fit/information.rs) `fill_read_slopes` —
in a cohort holding a sample of several libraries every other error is computed as if its rates were
known.** **Categories:** design, correctness, reliability (convergent). Zero slopes drop those rates
out of the inversion. Measured by zeroing every sample's rate scores on one-library cohorts of 4,000
positions: the mismapped share's error falls to 0.26 of its value at 4 samples and 3 reads, 0.70 at 8
samples and 3 reads, 0.92 at 8 and 8 reads, 0.90 at 20 and 3; the homozygote excesses' by at most 1%.
The report says the excess "keeps its error", which is true only to that 1%. The plan's split clause
says a multi-library sample's *errors* are reported absent. **Fix:** report absent every error of a
cohort that holds a sample of several libraries, for the life of this interim, or measure and state the
understatement.

**M2: [information.rs](../../../../src/parameter_estimation/joint/fit/information.rs)
`a_samples_libraries_are_scored_each_under_its_own_rates` — does not rule out scoring every library
under one library's rates.** **Categories:** design, reliability. Scoring each library under its
sample's *last* library's rates passes all 60 tests; swapping two libraries' rates passes the whole-fit
test and is caught only by the `> 10` bound on one slope (−3.88 under the swap). **Fix:** assert what
correct scoring implies exactly: an empty second library leaves the likelihood the first library's alone
whatever the second's rates; and swapping two libraries' reads together with their rates leaves the
likelihood unchanged.

**M3: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) the per-position reader — no test
puts a sample of several libraries where a depth is a range, or whose libraries differ in depth; nine
mutants survive.** **Categories:** reliability, correctness. Which library's coverage centres the
weights, which slot the weights go in, the per-library coverage function itself, and the rate each
branch's expected reference reads use. No fixture has one library walked where its sibling is not.
**Fix:** a test on a drawn two-library cohort at different depths under a cap above 124: each library's
coverage near its drawn depth, and each ranged library's weights equal to `fill_depth_weights` at its
own coverage. **Context:** with the shipped cap of 124 no stored depth is a range, so this code moves
nothing today (`src/run/gatherer.rs`).

**M4: [fit.rs](../../../../src/parameter_estimation/joint/fit.rs) `one_position` — a sample's count of
positions with reads is untested for several libraries.** **Categories:** all three. Counting only
the first library passed the suite and changed a count from 400 to 330 (and 2,848 to 2,335 on another
probe); counting the deepest library instead of the sum lowered positions with two reads from 203 to
157. The count is written to the parameters file and is the per-sample estimates' `observations`.
**Fix:** assert it against the drawn reads.

### Minor

- **Mi1: the report's "up to 482 of kimura's 2,651 read groups"** counts the read groups after each
  sample's first; every library of a multi-library sample moves, so the figure lies between 483 and
  964. (The plan carries the same figure; not the implementer's to edit.)
- **Mi2: the report says per-library coverage moves numbers "above the census cap of 124"**; a stored
  depth is a range only when the run's cap is above 124, and the shipped cap is 124, so it moves
  nothing today.
- **Mi3: the empty-library branch in `one_position` is reachable and unreported.** A sample holding
  no ordinary-position section (a tracts-only census) panicked before at `group_index[s][0]`; now the
  fit completes and reports a homozygote excess from no reads, marked fitted. Untested.
- **Mi4: removing the `with_several_libraries_not_scored` call in `fit_jointly` passes the suite.**
- **Mi5: the interim rule is decided in four places**, and `named` re-derives it instead of reading
  the marked rows, so the log and the trace can drift.
- **Mi6: stale doc** — `sample::CLEAN_ERROR_RATE` still says "Its first read group's error rate".
- **Mi7: the log's new reason reads badly** — "(6 a sample of several libraries, not scored yet)" —
  and "scored" means nothing to the reader.
- **Mi8: the report states no cost**; measured 2.3 to 3.6% slower for one library a sample (§3).
- **Mi9: `Estimate::observations` of a library's noise is the cohort's position count**, now that a
  library's rates are fitted from its own reads alone (pre-existing semantics, sharper now).
- **Mi10: `Scratch::new` takes three counts by position**, and callers take the library count from
  two sources.

### Nits

"Library" and "read group" name one index in turn; a few `sample` bindings now hold a library; the
`section == 0` special case on `sample_invariant` is not needed for bit identity (measured).

## 7. Out of scope observations

- A read group no sample holds gets zero tallies and is reported fitted, with the run's position count
  as its observations (pre-existing).
- The position count is read from the first sample, so a first sample with no ordinary section gives a
  fit of no positions (pre-existing, not probed).
- A single sample with two libraries fits (clean rates 0.00164 and 0.00772 against 0.002 and 0.008
  drawn); no test covers the shape.

## 8. Missing tests to add now

- `each_librarys_tallies_hold_its_own_reads_in_every_branch` (B1).
- `an_empty_library_leaves_its_samples_likelihood_to_its_sibling` and
  `swapping_two_libraries_reads_and_rates_leaves_the_likelihood` (M2).
- `each_library_is_weighted_around_its_own_depth` (M3).
- The two-library whole-fit test asserting positions with reads (M4) and the interim errors (Mi4).

## 9. What's good

- The first-library-assigns, later-libraries-add shape keeps one-library arithmetic byte-identical,
  verified on nine cohorts at two thread counts.
- The generator draws a library's read count before the sample's genotype, so every existing fixture
  is the same cohort.

## 10. Commands to re-verify

- `scripts/dev.sh cargo test --release --lib parameter_estimation::joint::fit`
- Per-category files: `tmp/review_2026-09-28_fit_precision_a6a/{design_claims,reliability_errors,correctness_cost}.md`;
  the correctness reviewer's probes in `…/correctness_evidence/`.
