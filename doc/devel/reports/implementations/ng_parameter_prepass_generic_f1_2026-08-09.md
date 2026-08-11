# ng step 4, the SNP/indel path — F1: the two ways in

**Date:** 2026-08-09. **Plan:** `doc/devel/ng/impl_plan/parameter_prepass_generic.md`, F1.
**Design:** `doc/devel/ng/arch/parameter_prepass_generic.md` §1.1.

## What was built

`GenericEstimationConfig`, `estimate_generic_parameters(loci, config)` and
`GenericAccumulators::estimate(config)`, in a new file `generic/estimate.rs`.

**The first is literally the second over an accumulator fed by the stream.** Not a parallel
implementation with a test asserting agreement — the same call. The plan's "so the two cannot
diverge" is therefore structural, and the test that asserts them equal is confirming the
wiring rather than holding the property up.

**`GenericEstimationConfig::accumulators()`** mints a shard's accumulator from the config, so
every shard of a sample is handed the same `edges` object — which `merge` requires and proves
by `Arc::ptr_eq` rather than by comparing lengths.

**`GenericAccumulators::covered_positions()`**, the warrant a supplied `F` carries. **From the
windowed table and not the read-group one**: a site covered by two libraries enters the
read-group table twice, and its positions would be counted twice with it.

## The three debts Milestone E recorded are discharged

- **`fit_coupled` has a production caller.** It was called only by its own test.
- **`take_supplied_inbreeding` has one**, and its warrant is asserted against a number a
  reader can compute: the fixture's loci are one base each, so covered positions equal the
  site count.
- **`fallback_error_rates` is exercised non-empty**, which took a specific fixture shape.
  The ladder is *fitted here → borrowed → supplied → defaulted*, so a supplied rate is
  consulted only when a group is too thin to fit **and** no sibling qualifies to lend. Three
  read groups covering disjoint thirds of the sites, each one site below the floor, is that
  shape — while their union is comfortably above the floor the genotype frequencies are
  checked against. Two groups supplied and one not, so the last two rungs are both reached in
  one sample.

## What it is proven against

- **`the_two_entry_points_answer_identically`** — asserted as whole `GenericSampleParameters`
  values rather than field by field, so a field added later is covered without revisiting it.
  **What it cannot say**: the two share `estimate`, so it proves the stream-walking half feeds
  the accumulator the same loci, not that the reduction is right. That is F2's question.
- **`a_failed_walk_propagates_rather_than_fitting_the_prefix`** — the stream yields the whole
  sample and *then* fails, so the prefix is above every floor and would have fitted. That is
  what makes it a test of the propagation rather than of the floors.
- **`a_supplied_inbreeding_coefficient_is_taken_without_a_fit`** — value, provenance, warrant,
  and that no runs model is reported.
- **`a_supplied_rate_is_used_where_no_group_can_fit_or_lend`** — the `Supplied` and
  `Defaulted` rungs, from the config field, in one sample.

## Recorded deviations

- **A new error variant, `LocusGeneration { sample, cause }`.** Arch §1.1 mandates both
  `Result<_, ParameterEstimationError>` and that a `LocusGenerationError` propagates; a
  variant is the only way to satisfy both. **Its cause is a rendered `String` and that is a
  real loss**: `LocusGenerationError` is `Debug + Error` and neither `Clone` nor `PartialEq`,
  which this enum is and its tests rely on, so carrying it would mean dropping both from every
  variant. What is lost is `source()` chaining — a caller can read what failed but cannot match
  on which stage. The field is deliberately **not** named `source`, which would claim a chain
  it does not provide (and which `thiserror` would try to treat as one).
- **A new file rather than a sixth job for `generic/mod.rs`**, which Milestone C's review
  already flagged as doing four. The arch module table has no entry-point row either way, so
  this is the same shape as `depth_bins.rs` was: an owner call, recorded rather than taken.
- **`library_shares` widened from private to `pub(super)`.** The runs model needs each
  library's share paired with its *settled* rate — including a borrowed or supplied one, since
  the alternative is a library scored against a rate nobody chose — and the shares are
  computed in `coupled_fit.rs`. The pairing is by read group and never by position, which is
  the same discipline `noise_from` keeps and for the same reason: a rule with two libraries'
  rates swapped between their groups is still a probability over the cell space, so none of
  the scoring rule's identities can see it.

## What F1 does not establish

- **Nothing here checks what the reduction computes.** Every number still comes from the fits
  proven in Milestones D and E; F1 only says the loci reach them. F2 is where a directly
  filled accumulator is refitted and the answer checked against known parameters.
- **`read_admission` is recorded and not read.** Its consumer is the `SampleSummary` assembly,
  which belongs to the cohort gather and is out of this plan's scope (spec §7). It is on the
  config now rather than later because the config is where a caller states it, and a field
  added afterwards would be one every existing caller had to be revisited to fill.
- **The runs model is run with `RunsModelStarts::default()` and the config cannot change it.**
  Arch §1.1 does not list a field for it and nothing in either cohort needs one.
