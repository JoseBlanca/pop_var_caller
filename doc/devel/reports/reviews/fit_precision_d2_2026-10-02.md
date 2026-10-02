# Code Review: fit_precision_d2
**Date:** 2026-10-02
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step D2 — a large cohort's repeat-tract strata fitted on a growing subset of samples
**Status:** Request-changes (applied: [fixes_applied_fit_precision_d2_2026-10-02.md](fixes_applied_fit_precision_d2_2026-10-02.md))

---

## 1. Scope

- **What was reviewed:** the diff `cbb27a9c..3e5b9dff` (a review object of the working tree): in `ssr_fit.rs` the
  subset settings, `StratumFit::samples_fitted_on`, `climb_from_parameters`, `StratumEvidence::of_samples` and
  `readers_by_group`, `subset_members`, `level_measured_to`, `parameters_of_fit`, `fit_from`,
  `fit_on_growing_subsets`, `fit_strata_on_sample_subsets`, `fit_strata_with`, `subsets_summary`, seven tests; the
  caller in `run/census_fit.rs`; fixtures.
- **Categories:** two reviewers, each in its own worktree — *correctness and tests*; *design and claims*.
  Findings: `tmp/review_2026-10-02_fit_precision_d2/{correctness_tests,design_claims}.md`.

## 2. Verdict

Request-changes. The loop terminates, the indices line up, the two schedules give the same bits and a cohort of
256 or fewer is untouched; but the refusal floor on the new path was untested, the subsets did not nest, a
stratum could be answered from fewer tracts than the floor, and one walk from the last answer settled up to 41
log-likelihood units below three fresh starts at thirteen classes.

## 3. Execution status

- Module tests on the reviewed tree: 68 passed (both reviewers).
- Mutations: 17 (correctness and tests; 11 survived, 9 of them changing output on some input) and 6 (design and
  claims; 3 survived, each changing behaviour on some input).
- A probe refitting one drawn stratum's subsets of 16 to 256 samples both from three starts and from the last
  answer, at three and thirteen classes.

## 4. Open questions and assumptions

1. Six design choices for the owner (M1, M3, M6, Mi7–Mi9) — answered 2026-10-02: all accepted.

## 5. Top 3 priorities

1. **B1, M3** — test the floor on the subset path, and keep a thin subset from being fitted.
2. **M1** — fit the answer's subset from every starting point.
3. **M2** — make the subsets nest.

## 6. Findings

### Blocker

**B1: the refusal floor on the subset path is untested.** *(correctness and tests)* Deleting it, or judging it on
the first subset (the trap spec §4.4 names), passed every test.

### Major

- **M1:** one walk from the last answer can settle below three fresh starts — at thirteen classes, 1.3 to 41
  units below in all four comparisons, once 0.79 of an error off in the level. *(design and claims)*
- **M2:** the subsets do not nest: a group's added readers left the subset when it doubled (8 readers to 1).
  *(correctness and tests)*
- **M3:** a stratum can be answered from fewer tracts than the floor (4 against 8). *(correctness and tests)*
- **M4:** the ≤ 256 test compares `fit_strata` with itself, and the 256/257 boundary is untested. *(both)*
- **M5:** the warm start, the doubling and the stop at every sample are untested: starting cold, quadrupling and
  `>` for `>=` survive. *(correctness and tests)*
- **M6:** the evidence counts now mean the subset's, but their docs and the parameters file say the stratum's;
  `samples_fitted_on` reaches only the log. *(design and claims)*

### Minor

- **Mi1:** the new section split `refused_before_any_walk` from its doc comment (clippy: empty line after doc
  comment), and a negated float comparison. *(both)*
- **Mi2:** "a subset none of whose samples read the stratum is passed over" describes a path that cannot run; the
  test named for it passes by the top-up. *(both)*
- **Mi3:** the subset is a full copy beside the whole evidence, about 1.4 GB on kimura's largest stratum when it
  grows to every sample. *(design and claims)*
- **Mi4:** `sample_order: &[usize]` is unchecked; `SsrFitConfig::subsets` is silently ignored by `fit_strata`.
  *(both)*
- **Mi5:** when a sample counts as a group's reader is untested at low depth. *(correctness and tests)*
- **Mi6:** `subsets_summary` has no test. *(correctness and tests)*
- **Mi7:** a group is topped up only when the first samples hold none of its reads. *(design and claims)*
- **Mi8:** every live group's level must reach the target, even one whose readers are all in. *(design and claims)*
- **Mi9:** at 2,169 samples the 2,048 step is 94% of the cohort. *(design and claims)*
- **Mi10:** under several strata at once, each stratum's walks run in sequence on one thread. *(design and claims)*
- **Mi11:** a subset-fitted stratum's weight in the curves changes, and spec §6 does not say so. *(design and
  claims)*

### Nits

`SampleSubsets::first`, `level_measured_to`, `fit_from` and `cohort: usize` named loosely; `first.max(1)` turns 0
into 1; the subset rule not logged; plan labels in doc comments; the summary's upper median; no progress line in
the parallel arm; a zero level with zero error untested; the guard counts `of_samples` keeps are never read; after
growth `walks` holds only the last subset's walks.

## 7. Out of scope observations

- The curve blend multiplies one group's level by reads crossing summed over every group — code that predates the
  step, moot while every read group is pooled.

## 8. Missing tests to add now

The floor on the subset path; nesting; a thin subset growing; the 256/257 boundary; the growth schedule; the final
refit; groups judged only while growing; one read making a reader; the summary.

## 9. What's good

- At 256 samples or fewer the new entry returns `fit_strata` itself, so the old path cannot drift.
- The sample indices of `sample_order`, the evidence and the homozygote excess line up by construction, both
  coming from `cohort.sample_names()` in one place.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib parameter_estimation::joint::ssr_fit`.
