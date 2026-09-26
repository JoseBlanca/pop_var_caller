# Issue 1: the repeat-tract evidence is released when the fit returns

*Implementation report, 2026-09-26. Plan:
[estimation_memory.md](../../implementation_plans/estimation_memory.md) §4. Branch
`evidence-release`.*

## Plan

`estimate-parameters` kept every sample's repeat-tract evidence (`Vec<StratumEvidence>`) alive
after the per-stratum fit, inside `CohortFit`, until the parameters file had been assembled. The
only thing read from it after the fit is each stratum's substitution rate, which is two counts a
stratum. This step keeps those two counts and drops the evidence as soon as `fit_strata` returns.

It does not lower the command's peak: the fit reads the evidence throughout, so the evidence and
the fit's likelihood tables are resident together whatever happens afterwards (plan §4). What it
frees is the evidence during the assembly of the parameters file — 0.7 GiB at 100 tomato samples
on kimura, and growing with every sample.

## Assumptions

- **Where the new type lives.** The plan names it (`StratumSubstitutionCounts`) but not its home. It
  is in `ssr_fit.rs`, beside `StratumEvidence`, because it is the part of that type that outlives
  it, and `StratumEvidence::substitution_counts()` builds it.
- **One division, not two.** `StratumEvidence::substitution_rate()` now delegates to the new type's
  `substitution_rate()`, so the rate the parameters file records and the rate a harness prints
  from the evidence (`examples/ng_joint_records_walk.rs`) are the same expression.

## Changes made

- [ssr_fit.rs](../../../../src/parameter_estimation/joint/ssr_fit.rs): new
  `StratumSubstitutionCounts { stratum, bases_compared, mismatching_bases }` with
  `substitution_rate()`; new `StratumEvidence::substitution_counts()`.
- [census_fit.rs](../../../../src/run/census_fit.rs): `CohortFit::tract_evidence` is replaced by
  `CohortFit::substitution_counts: Vec<StratumSubstitutionCounts>`, built in `fit_a_cohort` right
  after `fit_strata`, after which the evidence is dropped explicitly. `parameters_from_the_fit`
  reads the counts where it read the evidence; its loop is otherwise unchanged.
- No other reader of `tract_evidence` existed: `grep -rn tract_evidence src examples tests benches`
  found only `census_fit.rs` (the other hits are an unrelated function in
  `summarise_condition.rs`).

## Tests added

After review (see [the review](../reviews/estimation_memory_issue1_2026-09-26.md)), two more:

- `census_fit::tests::a_cohort_of_censuses_is_fitted_both_halves` now also asserts one count a
  stratum in the fit's order, at least one stratum with bases compared, and the same counts when
  the cohort is fitted twice.
- `census_fit::writing_the_parameters_file::a_stratum_with_nothing_compared_gets_no_rate_and_one_with_mismatches_gets_its_own`
  — adds two strata's counts to a real fit: one with nothing compared gets no rate (not a fitted
  zero), one with 3 disagreeing in 4,000 gets that rate for every read group, over 4,000
  observations.

From the first pass:

- `ssr_fit::tests::substitution_counts_give_the_rate_the_evidence_gives` — the counts carry the
  stratum and both counts, give the rate the evidence gives (3 in 4,000), and give `None` where no
  base was compared, as the evidence does.
- The existing tests of `parameters_from_the_fit` and `parameters_file_of`, and the cross-platform
  checksum test (`cli::cross_platform_digests`), cover the loop that reads the counts. **The
  checksum test passes unchanged** (`FITTED_PARAMETERS_MD5` and `CALLS_MD5` not touched).

## Validation

In the dev container (Apple `container`, arm64 Linux):

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo test --all-targets --all-features --no-fail-fast` | 4,865 passed, 3 failed |

The library target alone: 4,749 passed, 0 failed, 4 ignored; that includes the checksum test.

**The three failures are in two probe examples this step does not touch**, and neither file
references `ssr_fit`, `census_fit` or the new type:

- `examples/ng_generic_loci_dump.rs`: `a_deletion_across_a_region_boundary_keeps_the_support_a_single_region_walk_gives`
  and `only_the_rows_that_departed_from_the_reference_carry_a_chain_id`;
- `examples/ng_ssr_loci_dump.rs`: `a_cap_above_the_depth_is_invisible_and_a_cap_below_it_bites`,
  which panics building its fixture's catalog (`FlankExceedsBundleThreshold { flank_bp: 30,
  bundle_threshold: 15 }`), before any fit.

Both examples' failures are already recorded in `PROJECT_STATUS.md`, in the block "Step 5/6 — STR
observations through a run" (corrected after review: this report first placed them in the census
block).

## Tradeoffs and follow-ups

- The peak is unchanged by design; issue 2 bounds the tables and a smaller evidence layout is the
  next step after this plan (plan §8).
