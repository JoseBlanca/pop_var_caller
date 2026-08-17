# ng cohort merge — B3: `CohortObservation` and per-sample support

*Implementation report, 2026-08-17. Step B3 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[arch](../../ng/arch/cohort_merge.md) §4 and [spec](../../ng/spec/cohort_merge.md) §4.2, with
the owner's ruling of 2026-08-17 on how a record's quality sums are divided.*

## 1. Plan

Per sample, support against the allele table B2 built: reads and their quality sums summed
where two of that sample's own observations reached the same allele, never merged across
alleles, and a sample with no coverage keeping no support at all — which stays a different
fact from reference-only support.

## 2. The ruling this implements, and the one number that is exact

**Read counts are exact.** Every read is named (B0), so it lands on exactly one allele of the
locus and is counted there.

**The five quality sums are not, wherever a locus spans several of a sample's records.** The
mint stores them already summed over the reads behind one observed sequence — `q_sum`,
`num_fwd`, `mapq_sum`, `mapq_sum_sq`, `placed_left` — so when those reads take different
paths across the locus the sums have to be divided, and nothing recorded says how much
belongs to which read. The owner ruled for a proportional split (2026-08-17), after asking
what freebayes does. What it does is not available to us: it holds one object per read all
the way to the likelihood (`freebayes/src/Sample.h`, a vector of per-read observations per
allele), so the question never arises there — that is the storage cost B0 declined when it
chose to carry one identifier per read rather than five numbers.

So, per read composed across records:

- **the quality takes the weakest sighting's mean read** — the allele is wrong if any of its
  pieces is, so a read's evidence for it is limited by the piece it saw least well;
- **strand, placement and mapping quality take the read's pooled share**, `(Σ sums)/(Σ read
  counts)` over the sequences it was seen in. Mapping quality and strand are properties of
  the read, identical at each sighting, so the pooling estimates one number rather than
  mixing two;
- **the divided counts are rounded once per allele**, not once per read.

At a sample with one record — every STR locus and the ordinary generic one — nothing is
divided and every sum is the mint's own.

## 3. A divergence from production, deliberate and worth knowing

Production faces this in `project_compound_scalars`
([`per_group_merger.rs`](../../../../src/var_calling/per_group_merger.rs)) and takes the
**`min`** over the constituents' mean `q_sum`. Its own plan says the opposite:

> `S_s(C) ≈ |chain_id ∩| × min_over_constituents(q_sum_at_record_i / num_obs_at_record_i)`
> … with `min_over_constituents` because … **the cross-record compound's effective quality
> cannot exceed any single constituent's**
> ([`cohort_per_group_merger.md`](../../implementation_plans/cohort_per_group_merger.md),
> step 3)

**`q_sum` is a sum of `ln P(error)`** ([`pileup_record.rs:47`](../../../../src/pileup_record.rs)),
so it is negative and the number **nearer zero is the worse read**. `min` therefore picks the
constituent the read saw *best*, which is the opposite of the sentence justifying it: it
makes an allele spanning several records look better evidenced than any single piece of it.
ng takes the maximum — the weakest sighting — and pins it with
`the_weakest_sighting_sets_a_composed_reads_quality`. Production is frozen and nothing here
changes it; the divergence is recorded at the code and raised for the owner.

## 4. Two departures from the architecture's declarations

Both are recorded at the types, and neither changes what the step delivers.

- **`per_sample` is sparse.** Arch §4 says `Vec<SampleSupport>` "indexed by the run's sample
  order", and its next sentence says a sample with no coverage "has no support at all, which
  is a different fact from reference-only support and stays one". Those two readings pull
  apart at k = 3,000 samples where three cover a locus. The code keeps only the covering
  samples, in ascending order, each naming its own `sample` index — the shape the walk hands
  over ([`SampleMembers`](../../../../src/ng/run/cohort_merge/close.rs)) — so the
  distinction is structural rather than resting on a zeroed row.
- **`per_allele` holds a dedicated `AlleleSupport`, not `SequenceObservation`.** Three of
  that type's fields cannot be right per allele: `bases` (the table already holds them),
  `read_witness` (always `Complete` here) and `read_group` — an allele's support aggregates
  over read groups, so no single value is true, and writing one would be a fabricated field
  rather than an approximation. **The consequence is that the read-group axis does not
  survive into a cohort observation**; it survives where it is used, in the census the
  parameter pre-pass fits. Flagged for the owner rather than assumed.

## 5. Changes made

All in [`build.rs`](../../../../src/ng/run/cohort_merge/build.rs):

- **`CohortObservation`** — `over(&ClosedLocus)`: the locus's ground, its alleles, and the
  covering samples' support.
- **`SampleSupport`** — the sample it belongs to, its row against the allele table, and three
  counts: the reads that said nothing (carried through from the records), the reads removed
  as evidence, and **the reads whose sums were divided rather than measured**, which is what
  tells a later step how approximate the row is.
- **`AlleleSupport`** — reads and the five sums, per allele.
- **`AlleleTable::assemble`** — the table and the support in **one pass**: interning an
  allele answers which allele it is, so a sample's support accumulates as its alleles are
  derived. The alternative — deriving again and looking the bytes up — would compose every
  read's allele a second time, which is the expensive half of this module.
  `AlleleTable::over` is now that walk with the support discarded.
- **`ShownBy`** — what backed an allele the derivation emitted, so the caller need not
  compose it again to attribute support: every read behind one sequence of a sole record, or
  one read with its sighting at each of the sample's records.
- **`share_of_one_read`**, **`AlleleSupportTally`**, **`round_to_u32`/`round_to_u64`** — the
  division, the accumulation in `f64` so that rounding happens once, and the saturating
  conversion back.

## 6. What the tests pin

11 tests added (86 → 97 in the module) — 8 with the step and 3 from
[the review](../reviews/ng_cohort_merge_b3_2026-08-17.md), which found that **the case this
whole division exists for had no test**: one observation's reads splitting onto two alleles.

- **one record is exact** — four reads' numbers land whole on the SNP's allele and two on the
  reference's;
- **summed within an allele, never across** — the same bases from two read groups are one
  allele and their sums add, while a different allele keeps its own;
- **the division, with its numbers spelled out** — the composed read's quality is the weaker
  of −2.0 and −0.5 per read, and its four pooled shares are 2/3, 170/3, 9,700/3 and 1/3;
- **the weakest sighting governs** — −1.0 against production's −6.0, five nats apart;
- **rounding happens once per allele** — three reads holding half a left-placed read each
  make 2, where rounding per read would claim all three started left;
- **a sample with no coverage has no entry**, walked through `LocusCloser` rather than
  hand-built;
- **every row is as wide as the final table**, including a sample derived before a later one
  introduced an allele;
- **the reads that said nothing are carried through** and summed across a sample's records.

**Two of my own expectations were wrong and the code was right both times**, which is why
they are worth naming: I had the sign of the quality comparison backwards in one fixture's
prose (−0.5 is the *worse* read, not −2.0), and claimed a left-placed count of 1 where the
fixture gives 1.5 → 2. Both were caught by running the tests, and the second turned out to be
the sharper demonstration of rounding once per allele.

## 7. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean. (`--all-targets` is red on this
  branch and was before this work.)
- `cargo test --lib ng::run::cohort_merge` — 97 passed, 0 failed.
- `cargo test --lib` — **3,720 passed, 0 failed, 11 ignored** (3,709 before this step).

## 8. What B3 does not do

- **Nothing consumes the counts yet.** Where the removed reads and the divided-sum count
  surface is C1's (what a region reports) and the emission step's (spec §13).
- **No builder, no cache, no threads** — C1 onwards.
- **The cross-record derivation still composes once per read.** The B2 review measured a
  prototype composing once per distinct read *pattern* at 39% cheaper on a locus of 1,000
  samples × 300 reads; it belongs here or at C1, and it needs the read multiplicity this step
  now has.
</content>
