# ng cohort merge — A3: the reach walk, closing loci only

*Implementation report, 2026-08-17. Step A3 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §4.1, §4.2, §9, §11 and
[arch](../../ng/arch/cohort_merge.md) §3.*

## 1. Plan

`LocusCloser`: one walk over the samples' observations merged by position, closing each
locus as its reach stops growing. Every closed locus comes out — judging them is A4.

## 2. Assumptions and recorded deviations

- **`ClosedLocus::members` is a `Vec<SampleMembers>`, not the `&'a [SampleObservationRef]`
  the arch sketches.** An `Iterator` cannot yield an item borrowing from the iterator
  itself, so a slice into closer-owned scratch is not expressible. The shape that landed
  keeps the arch's intent — *the members are borrowed from what the walk was given, never
  copied out of it* — and improves on it in one way: since each sample's observations
  inside a locus are **contiguous**, one entry per covering sample carries a borrowed
  slice, rather than one entry per observation. At a locus every sample covers, that is
  k entries instead of one per observation, and it is already the grouping B3 needs.
- **`LocusCloser::over` takes `&[&'a [SampleLocusObservations]]` and copies the k slice
  references in.** Borrowing the outer array as well would tie every caller to keeping
  that array alive beside the data, which the cache's scoped `with_observations` (arch
  §4) does not do. One k-pointer allocation per walk.
- **`ClosedLocus` carries `non_reference_reads`.** The plan's A3 line does not name it,
  but the arch's walk keeps three running values — start, reach and the non-reference
  total — and summing it here is what makes A4 a verdict rather than a second pass. A4
  reads it; nothing else changes.
- **The verdict is not here.** A3's contract is explicit that no locus is judged, so
  `ClosedLocus` has no `verdict` field; A4 adds it, which is what makes A4 bisectable on
  its own.

**One precondition, and the review changed how it is held:** each sample's slice must be
in coordinate order. The first draft documented it and checked nothing, on the grounds
that a check costs a comparison per observation on the hot path. The review showed what
violating it does — a locus closed over ground it does not cover, silently — and the check
now runs, at release level, inside the walk. One comparison against per-observation work
that already includes a base-by-base sequence comparison.

## 3. Changes made

**[`src/ng/run/cohort_merge/close.rs`](../../../../src/ng/run/cohort_merge/close.rs)**
(new) — `LocusCloser`, `ClosedLocus`, `SampleMembers`.

The walk keeps three running values and takes heads while the next one begins inside the
reach:

- the argmin over per-sample heads, ties to the lowest sample index — the tie-break is
  what makes the members' order a function of the data (spec §9);
- the reach extends *as the walk scans*, which is what makes one pass enough however
  late a deletion widens the locus;
- a change of contig closes the locus whatever the positions say, since no observation
  crosses one.

Members are recovered without a second scan: the walk records where each cursor stood
when the locus opened, and what it consumed by the close is that sample's run.

## 4. Tests added

Twenty-three, in `close.rs` — fourteen written with the step, nine added by the review.
The three the plan names, first:

- `two_adjacent_snps_are_two_loci` — spec §4.1's own example; neither observation covers
  the other's base.
- `a_deletion_and_a_snp_inside_it_are_one_locus` — the other example, with the two in
  *different* samples, which is the case the merge exists for.
- `a_late_widening_pulls_in_a_following_observation` — 10–12, then a deletion at 11
  reaching 20, then a SNP at 15. The third is inside the locus only because the second
  widened it: when the walk passed position 10 the reach was 12, and 15 is well beyond
  it. A fixed reach closes 10–12 and starts a new locus at 15.

Then the properties the spec asks for: `the_loci_do_not_depend_on_which_sample_carried_what`
(the same five observations dealt out three ways, including one where the widener and
what it pulls in sit in different samples); `closed_loci_are_disjoint_and_ascending`
(and every observation a member exactly once); `a_locus_never_crosses_a_contig`;
`closing_is_uncapped`; `an_insertion_does_not_widen_a_locus`;
`the_non_reference_total_is_summed_across_the_locus` (with nine reference reads present,
so summing `num_obs` instead answers 18 rather than 9);
`an_all_reference_locus_still_closes`; `nothing_to_walk_yields_nothing` (k = 0 and empty
samples); and the two member-shape tests.

## 4a. What the review changed

The step went through `rust-code-review` (five categories, two agents mutation-testing in
worktrees) and `apply-code-review-fixes` before it was committed —
[review](../reviews/ng_cohort_merge_a3_2026-08-17.md),
[fixes](../reviews/fixes_applied_ng_cohort_merge_a3_2026-08-17.md). Five things below are
the review's, and the first two are defects the first draft shipped:

- **`span()` called `region.len()`**, which A2 had documented one step earlier as
  overflowing at the coordinate ceiling — a debug panic and **0 in release**. `span()` is
  what A4's width verdict reads, so the widest locus expressible would have been judged
  narrow enough to build. It now subtracts before adding.
- **The coordinate-order precondition was documented three times and checked nowhere**,
  and violating it does not fail loudly: it closes a locus over ground it does not cover.
  There is now a release-level check inside the walk, and a test.
- **Two mutants survived the original suite** — members emitted in consumption order, and
  `saturating_add` replaced by `+`. Both now have tests that kill them.
- **The walk's layout follows the sibling merge**: keys beside the cursors rather than read
  through them, scratch reused instead of a `Vec` cloned per locus, and the member vector
  sized exactly during the walk.
- **Two doc comments named mechanisms the code does not use**, including one crediting the
  tie-break with a determinism that actually comes from the member collection.

## 5. Validation

In the container: `cargo fmt --check` clean; `cargo clippy --lib --all-features --
-D warnings` clean; `cargo test --lib ng::run::cohort_merge` **25 passed, 0 failed** (16
before the review's fixes). Full suite figures are in the commit message.

**The walk was checked against production's, not merely argued to match it.** The review's
refactor_safety agent ported production's `CohortSpanFold` + `derive_is_kept` grouping into
a temporary test and ran it against `LocusCloser` over **5,000 random cohorts**: loci and
membership identical in all 5,000, with **5,025 members joining at exactly the running
reach**, so the boundary the two rules share is genuinely exercised rather than dodged.

## 6. Tradeoffs and follow-ups

- **The argmin is still O(k) per observation** — the keys are contiguous now, but the scan
  is a scan. A heap would make it O(log k); the sequencing agreed with the review is to
  take that only after a measurement, and the first thing that runs this at scale is the
  region-width sweep (spec §14 Q1).
- **One allocation per closed locus** for the members, which is what `Iterator` costs: an
  item cannot borrow a buffer its iterator owns. If a caller wants it back, the shape is
  internal iteration, and A4 and B1 are where a real caller decides.
- **Partition-invariance and a k = 3,000 fixture belong at E4**, over the whole builder
  arrangement, rather than as unit tests here.
