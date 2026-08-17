# Fixes applied — ng cohort merge, B2

*2026-08-17, against [the review](ng_cohort_merge_b2_2026-08-17.md). All in
`src/ng/run/cohort_merge/build.rs`; nothing else in the tree was touched.*

## Behaviour changed

- **A sample's records are checked to be disjoint and ascending**, once per sample, in
  `LocusReferenceBases::over` — before any read is consulted. The composition's own backstop
  is reached only by a read named at both overlapping records, so the same producer defect
  was loud or silent depending on which reads the data carried; the silent half lost the
  sample's whole evidence at that locus. Carries the "must become a `RunError` on the psp
  path" note its two siblings carry.
- **Reads removed as evidence are counted.** `alleles_of_sample` answers how many it
  removed; `AlleleTable::reads_removed_as_evidence()` sums them over the locus. Nothing
  consumes it yet — where the loss surfaces is B3's or C1's — but it is no longer
  unrecoverable.
- **`REFERENCE_ALLELE`** names the reference's position in the table, and `over` asserts the
  reference landed there.

## Behaviour kept, claims corrected

- **The two branches of the derivation differ on one shape**, and the doc said they did not.
  A fragment whose mates overlap and disagree is two sequences of one record under one chain
  id: at a sole record both are emitted, across several the fragment is removed. Kept — with
  one record there is nothing to compose across, so each mate's sequence is already complete
  evidence, while composing across records needs the fragment as a unit — and now stated,
  with `a_read_showing_two_things_at_a_sole_record_still_contributes_both` pinning it beside
  the removal test.
- **"Known to have covered the whole locus" is not what the code decides**, and both the
  module header and `AlleleTable` said it was. The criterion is presence at every record
  *that sample* minted; ground outside those records is written from the reference because
  the sample minted nothing there. The two things the derivation does not claim are now
  stated at it, the fragment whose mates flank an unsequenced insert included.
- **`index_of` answers identity and only identity.** Its doc offered B3 a route —
  recompose and look the bytes up — that cannot recover the per-allele moments, which are
  summed per observation.
- **The chain-id assertion cannot fire on the STR path**, and the comment now gives the
  reason rather than calling it a shape that "somehow spanned two records": two STR records
  of one sample never chain, segments being the reference's own partition. The
  `num_obs == 0` escape is documented as guarding a state no producer emits.

## Tests

86 in the module, from 80 at review time. The two that answer the Blockers were each run
against the mutation they exist to kill:

- `an_allele_composed_across_records_closes_on_the_reference_it_consumed_not_its_length` —
  an insertion and a deletion inside multi-record compositions. Closing from `composed.len()`
  instead of the reference consumed now **fails** it (`ATGTG` for `ATGTGA`, `ACTTA` for
  `ACTT`); before, that mutation left all 80 tests green.
- `two_samples_that_each_hold_several_records_are_derived_independently` — **and its first
  version did not kill its mutant**, because the two samples used different read ids, where
  a carried-over buffer composes the same allele again and the table is unchanged. The id
  space is per file, so two samples sharing id 7 is the ordinary case; with that, dropping
  the buffer's reset now **fails** the test.
- `a_read_showing_two_things_at_a_sole_record_still_contributes_both`, the two
  overlapping-records refusals (with and without a shared read), the out-of-coordinate-order
  refusal, and `the_reads_removed_as_evidence_are_counted`.
- The depth-cap test now asserts the **sameness** of the table with and without a capped
  read, since the count it sets is never read here; before, it was its sibling plus one
  inert line.

## Shape, not behaviour

`ReadSighting` replaces the `(ChainId, u32, u32)` triple, its field order being the sort
order; `slice::chunk_by` replaces the hand-rolled run-grouping; `AlleleIndex` became
`AlleleLookup`; placements are built once per record rather than once per (read, record);
the working buffer is cleared on entry.

## Not applied

- **`Arc<[u8]>`** for the alleles, which would stop each distinct allele's bytes being held
  twice: it changes the type arch §4 fixes for `CohortObservation::alleles`, and B3 would
  hand the saving back when it converted. Recorded at the field.
- **Composing once per distinct read pattern** rather than per read — measured by the review
  at 39% of a cross-record locus's cost, at 1,000 samples × 300 reads. It also yields the
  read multiplicity per allele, which B3 needs, so it belongs there.
- **Hoisting the working buffer above the locus** — plan C1 owns it.
- **The psp obligations in arch §5.** Three assertions now instruct a later step to convert
  them and §5 names none; a design-doc edit, so it is offered to the owner rather than made.

## Validation

In the container:

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — 86 passed, 0 failed.
- `cargo test --lib` — 3,709 passed, 0 failed, 11 ignored (3,703 at review time).
</content>
