# Code Review: fit_precision_d3
**Date:** 2026-10-02
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** plan step D3 — the tool comparing chosen strata fitted on every sample and on the grown subset
**Status:** Approve-with-changes (applied: [fixes_applied_fit_precision_d3_2026-10-02.md](fixes_applied_fit_precision_d3_2026-10-02.md))

---

## 1. Scope

- **What was reviewed:** the diff `349c0532..836c696d` (a review object of the working tree):
  `cli::estimate_parameters::open_the_censuses` and `contig_id_of`, split out of `fit_and_assemble`;
  `run::refuse_another_selection` and `the_tract_strata_of_a_cohort`, split out of `fit_a_cohort`;
  `examples/ng_ssr_subset_against_every_sample.rs`.
- **Categories:** one reviewer in its own worktree — behaviour of the split, the tool's correctness, cost and
  prose. Findings: `tmp/review_2026-10-02_fit_precision_d3/review.md`.

## 2. Verdict

Approve-with-changes. The split changes no behaviour: the same checks in the same order, the selection checked
before the SNP/indel fit, the same tract count. The tool computes the right numbers on the right samples; two
things about how it reports them needed changing.

## 3. Execution status

`run::census_fit` 9 passed, `cli::estimate_parameters` 21 passed; clippy clean on the example; the tool run on the
four-sample tomato cohort with a first subset of 2 (670 s).

## 4. Open questions and assumptions

None.

## 5. Top 3 priorities

1. **M1** — "apart" when the subset grew to every sample.
2. **M2** — the default strata and the silent fit on every sample.

## 6. Findings

### Major

- **M1:** when the subset grows to every sample, the subset's error can exceed the whole fit's by a hair, and
  "apart" divides by nearly zero.
- **M2:** the default picks the five most expensive strata, and nothing is printed while each runs; on kimura that
  could be over a day silent. The checkpoint asks for strata "chosen to include the largest".

### Minor

- **Mi1:** `--inbreeding` is not checked as `estimate-parameters` checks it.
- **Mi2:** the doc's reading of "apart" overstates one large value, does not say what "-" means, and does not map
  to spec §4.5 item 3.
- **Mi3:** no total of the time saved, which spec §4.5 item 3 asks for.

### Nits

When one fit fails the other's result is dumped whole; a chosen stratum below the floor still spends a fit on
every sample first; "One line a stratum and number" is hard to parse; the header says "subset of N samples" when no
subset was drawn.

## 7. Out of scope observations

None.

## 8. Missing tests to add now

None: an example, run by hand.

## 9. What's good

- The tool opens the cohort through the same function `estimate-parameters` uses, so it cannot drift from the run.

## 10. Commands to re-verify

`scripts/dev.sh cargo test --release --lib run::census_fit`; `... cli::estimate_parameters`.
