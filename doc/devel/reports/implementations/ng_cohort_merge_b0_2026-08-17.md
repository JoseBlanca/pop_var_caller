# ng cohort merge — B0: every read is named, the reference-matching ones included

*Implementation report, 2026-08-17. Step B0 of [the plan](../../ng/impl_plan/cohort_merge.md),
added by the owner's ruling on spec §14 question 2. The code is in the generic locus generator,
one module upstream of the cohort merge, because that is where a read's identity is decided.*

## 1. The ruling this implements

> **Either we know the read covered the whole locus, and its allele is elongated with what it
> showed; or we know it did not cover it, and it is removed as evidence. Not being able to decide
> which is an error that must never happen.** — the owner, 2026-08-17

**The third case was the ordinary one.** A cohort locus can span several of a sample's records —
one sample's SNP at 12, another's deletion covering 10–14 — and a read that showed the SNP had to
be placed at every other position of that locus before its allele could be written. The generic
mint recorded a chain id only for a read that *disagreed* with the reference, on the argument that
"a chain id marks which haplotype a read came from, and the reference is the default — a default
needs no tag". So a read that covered position 14 and agreed with the reference there, and a read
that never reached 14, were the same absence.

**What the id is now for is different from what it was for.** It used to answer *which haplotype*;
it now also answers *was this read here*. That is why the old rule could not be kept.

## 2. Changes made

- [`open_record.rs`](../../../../src/ng/locus_generation/pileup/open_record.rs) — the general
  fold records `state.chain_id` for every read it folds, not only for the ones that departed.
- [`fast_column.rs`](../../../../src/ng/locus_generation/pileup/fast_column.rs) — the same for
  the ordinary-column fast lane, which had its own copy of the rule.
- **`read_agreed_with_reference` is deleted**, with the test that pinned its one subtlety
  (a witness with a hole never counted as agreeing). It existed only to gate the chain id, and
  with the gate gone it had no caller. Its reasoning is in this commit's parent if it is ever
  wanted again.
- The field doc on `KeyedObservation::chain_ids` now carries the ruling, what the id answers, and
  that the cost — an identifier per read per position, genome-wide — was accepted rather than
  avoided.

**In the parity harness against production** ([`parity.rs`](../../../../src/ng/locus_generation/pileup/parity.rs)),
which is what says ng's walk still computes what production's does:

- chain ids on a **reference-matching** observation are cleared on **both** sides before the
  loci are compared, since production still withholds the ids of its REF bucket. The same rule
  applied to both sides is this harness's own standing requirement; what stays compared is the
  ids of the observations that departed, where the two walkers still agree.
- `locus_chain_ids` — the locus-level comparison — is restricted the same way.
- **A new assertion on ng's own side**: an observation with reads carries at least one id. That is
  the property B0 buys, and production has no counterpart to compare it against.

## 3. What this surfaced, and it is not small

**A read that straddles a walk boundary has two identities.** The tiling test
`two_adjacent_regions_concatenate_into_the_single_region_walk` splits one region at position 50
and asserts the two halves emit what one walk does. It now fails on ids alone: the read the
whole-region walk calls 1 throughout is **1 before the cut and 3 after it**, because the split
walk meets it twice and allocates an id each time. Every other byte is identical.

**It is safe, and only because of where walks are cut.** The merge compares ids within one
sample's stream inside one cohort locus; a segment is never cut and no *observation* crosses one
([`run_streaming.md`](../../ng/spec/run_streaming.md) §4.3, the owner's rule of 2026-08-09), so no
chain of overlapping observations does either and the two identities of a straddling read always
land in different loci. **A run that ever cut inside a segment would break the merge's
read-linking silently** — a dependency the cohort merge now has on a rule made for another reason,
recorded here and at the comparison. **A second dependency of the same kind** rides with it: §6.2
makes segment independence conditional on every sample's segmentation being the same one, so two
samples segmented differently would put one sample's cut inside another's locus.

The two tiling tests therefore compare the evidence **up to the numbering of ids, per locus**:
identical bytes everywhere else, and inside each locus one consistent renaming, so a walk that
merged two reads into one identity or split one into two still fails. The helper carrying that
comparison is in [`pileup/mod.rs`](../../../../src/ng/locus_generation/pileup/mod.rs) with the
measured example in its doc.

## 4. What B0 does *not* do

- **The depth cap is left as it is.** It acts per position, so a read counted at one position and
  capped at the next is absent there for a reason that is not coverage. I raised this as needing
  identities; it does not: a capped read cannot be elongated either, because what it showed there
  was never recorded, so it takes the same branch as a read that never covered — removed as
  evidence. The outcome is decided; what is lost is a little depth at capped positions in
  multi-record loci, which is a high-depth effect. **The correction is mine and it is recorded
  because I told the owner otherwise.**
- **The STR path needs no ids and has none.** An STR locus is one record, so `ReadWitness` already
  says whether a read spanned it, and the rule is decidable there today.
- **Nothing about the psp encoding**, which does not exist yet. The per-file cost of carrying an
  id per read per position lands there when it is written, and the owner accepted it.

## 5. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::locus_generation` — 362 passed, 0 failed, 1 ignored.
- `cargo test --lib` — **3,679 passed, 0 failed, 11 ignored** (3,680 before; the deleted
  predicate took its one test with it).
- The differential against production's walker is part of that run and green:
  `every_divergence_from_production_is_one_of_the_six_named_classes`,
  `ng_agrees_with_production_where_production_fabricated_nothing`.

## 6. Tradeoffs and follow-ups

- **The memory and file cost is real, and one number for its size already existed.** Production's
  own walker records, beside the `allele_index == 0` rule this step departs from, that the ids it
  drops are about **96.6% of all chain ids on real cohorts**
  ([`pileup/walker/open_record.rs:155`](../../../../src/pileup/walker/open_record.rs)). ng's old
  rule was per read rather than per bucket, so the two do not withhold exactly the same set — but
  to within that, keeping them all is **about thirty times the ids** either caller used to carry. It is a `u64` per read per position, and **the price is not only the eight
  bytes**: an observation that used to hold an empty `Vec` now holds a populated one, so it is an
  allocation per observation per position as well. In the direct path both are transient; in the
  per-sample files the ids are the dominant term. Sizing it belongs with the measurement spec §8
  already owes for what one observation costs.
- **Two design documents now disagree with the code**: spec §14 Q2 still records the question as
  open with production's leaning, and `locus_generation`'s spec still describes the chain id as a
  tag the reference does not need. Both need the ruling written in; that is a design edit and not
  this run's to make silently.
