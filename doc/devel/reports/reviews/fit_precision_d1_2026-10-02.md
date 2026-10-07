# Code Review: fit_precision_d1
**Date:** 2026-10-02
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step D1 — the fixed order a repeat-tract stratum takes its samples in
**Status:** Approve-with-changes (applied: [fixes_applied_fit_precision_d1_2026-10-02.md](fixes_applied_fit_precision_d1_2026-10-02.md))

---

## 1. Scope

- **What was reviewed:** the diff `9854283e..c15c624f` (a review object of the working tree): the new module
  `src/parameter_estimation/joint/sample_order.rs` (`SAMPLE_ORDER_SEED`, `sample_order`, five tests) and its
  declaration in `joint/mod.rs`.
- **Categories:** one reviewer in its own worktree, covering correctness, tests, naming, defaults and idiom — the
  step is one pure function with no callers. Findings: `tmp/review_2026-10-02_fit_precision_d1/review.md`.

## 2. Verdict

Approve-with-changes. The order depends on the names and the seed only; the one real exposure is how a name's
bytes reach the hash (Mi1).

## 3. Execution status

- Module tests on the reviewed tree: 5 passed (Linux aarch64, dev container). x86_64 and the macOS host not run.
- Seven mutations, each against the five tests and a temporary 300-name test: five killed; dropping the name
  tie-break survives and changes nothing a test can build (it needs a 64-bit hash collision); an unstable sort
  survives the five tests and matters only for duplicate names, which a run refuses.

## 4. Open questions and assumptions

None.

## 5. Top 3 priorities

1. **Mi1** — feed the name's bytes to xxh3 directly.
2. **Mi3** — what D2 must pass into the fit.

## 6. Findings

### Minor

- **Mi1:** the name reaches xxh3 through `str`'s `Hash`, which appends a `0xFF` by an unstable default that the
  standard library says is "not yet decided"; the other candidate it documents, a length prefix, would make the
  bytes platform-dependent. Only the pinned test guards it. Fix: `xxh3_64_with_seed(name.as_bytes(), seed)`,
  pinned order re-recorded.
- **Mi2:** the module doc says the hash is over the name's bytes; it was over the bytes and a `0xFF`.
- **Mi3:** sample names reach neither `fit_strata` nor `StratumEvidence`. D2 should compute the order once in
  `census_fit.rs` from `cohort.sample_names()` — the same indices as `SampleTractReads.sample` and
  `homozygote_excess` — and pass the order, not the names; it will also need each sample's rank.

### Nits

The pinned test's doc reads as if platforms may differ; "the oracle cohort" is an internal label; the doc
describes duplicate names a run refuses (`RunError::SampleAppearsTwice`); the duplicate-name test is too small to
tell a stable sort from an unstable one.

## 7. Out of scope observations

- The census's `hash_position` and `repeat_catalog::strata::hash_locus` feed names through `str`'s `Hash` too
  (Mi1's exposure); changing them would move every census and catalog selection.

## 8. Missing tests to add now

None beyond the fixes.

## 9. What's good

- The pinned order is a real guard: feeding the bytes without the `0xFF` fails it.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib parameter_estimation::joint::sample_order`.
