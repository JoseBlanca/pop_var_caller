# ng cohort merge — A2: the two derivations, on the observation types

*Implementation report, 2026-08-17. Step A2 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[arch](../../ng/arch/cohort_merge.md) §2 and [spec](../../ng/spec/cohort_merge.md)
§4.1, §4.3, §11.*

## 1. Plan

The cohort merge's walk needs two numbers per observation — how far it reaches along
the reference, and how many of its reads showed something other than the reference —
and neither is a fact about merging. Put both on the observation types, put the
comparison they share in one place, and move the second caller onto it.

## 2. Assumptions and recorded deviations

**One decision the design documents do not make, and it changes what the keep rule
counts.** `non_reference_reads` sums over the **complete** observations only. A
`ReadWitness::Partial` observation's `bases` cover only the stretch its read witnessed,
so comparing them against the locus's whole `reference_bases` reports every partial as
non-reference — including a read that agreed with the reference over every base it
actually saw. Arch §2 says the sum is "`matches_reference` over the members"; spec §4.1
says the comparison "is exactly what the census writer does today", and the census
writer makes it over `complete_observations()`. The second is the more specific
instruction and the only one that yields a defensible number, so it is what landed.

**What it costs, stated plainly:** a variant witnessed only by partial reads never
reaches the keep threshold, so nothing is emitted over that locus at all. Scoring a
partial needs a censored likelihood that does not exist yet
([`locus_generation.md`](../../ng/spec/locus_generation.md) §3), so a partial cannot be
told from a short allele today either way. **Raised at Checkpoint A** — it is the
owner's to confirm, and it is cheap to change.

**`reach()` is `max(end, start)`, not production's `start + max(span, 1) − 1`.** The two
agree on every input. Production's form needs the span, and `GenomeRegion::len()`
computes `end + 1` before subtracting — which panics in a debug build on a region ending
at the top of the coordinate space, despite that method's doc promising to saturate.
This was **found by writing the test**, not predicted: the first draft used production's
form and `the_reach_is_productions_reach` failed with *attempt to add with overflow* at
[`types.rs:94`](../../../../src/ng/types.rs). The defect in `GenomeRegion::len` is left
alone — it is a shared type, unreachable from real coordinates, and out of this step's
blast radius — but a test now pins that `reach` itself answers there, with a comment
saying why production's expression cannot be checked against at that input.

## 3. Changes made

**[`src/ng/locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs)** —
three methods:

- `SequenceObservation::matches_reference(&self, reference_bases: &[u8]) -> bool`, the
  one definition of "non-reference" in the codebase. It is a byte comparison, which is
  all it has ever been; what is new is that there is one of it. Its doc states the
  precondition the census's use depends on: the slice handed in must cover the same
  stretch the observation's `bases` do.
- `SampleLocusObservations::reach(&self) -> Position` — the last reference base the
  locus covers, named rather than read off `region.end` because the merge groups by it
  and someone will want to compare ng's rule with production's in one place.
- `SampleLocusObservations::non_reference_reads(&self) -> u32` — `num_obs` summed,
  saturating, over the complete observations whose bases differ.

**[`src/ng/parameter_estimation/joint/census.rs`](../../../../src/ng/parameter_estimation/joint/census.rs)**
— `add_generic`'s inline `*observation.bases == *locus.reference_bases` becomes
`observation.matches_reference(&locus.reference_bases)`. Behaviour identical; the
comment says which subset makes the whole locus's reference bases the right stretch.

## 4. Tests added

Six, all in `locus_generation`'s `mod tests`:

- `the_reach_is_productions_reach` — a SNP and a deletion, each checked against
  production's arithmetic recomputed in the test rather than against a copied answer.
- `a_locus_at_the_coordinate_ceiling_reaches_its_own_end` — the saturating edge spec §11
  names, and the one input where production's expression cannot be evaluated (see §2).
- `an_inverted_region_reaches_its_own_start` — `GenomeRegion` has public fields and no
  constructor enforcing `start <= end`; a reach behind its own start would make the A3
  walk close every locus at once.
- `matches_reference_compares_the_bases_it_is_given` — reference, SNP and deletion
  against one stretch, plus the same deletion against its own.
- `non_reference_reads_sums_only_the_reads_that_differ` — 40 reference reads beside 3
  and 2 non-reference ones. The 40 are what make it discriminating: summing every
  observation gives 45, inverting the predicate gives 40, and only the right rule gives
  5.
- `a_partial_that_agreed_with_the_reference_is_not_counted_against_it` — the §2 decision,
  as a fixture: a read that stopped after two bases having matched both.
- `a_locus_with_no_observations_has_no_non_reference_reads` — no coverage is not an
  error.

**The census's own tests are the other half of this step's evidence.** They were not
touched: `cargo test --lib ng::parameter_estimation::joint` is 103 passed, 0 failed.

## 5. Validation

In the container: `cargo fmt --check` clean; `cargo clippy --lib --all-features --
-D warnings` clean; `cargo test --lib ng::locus_generation::tests` 26 passed, 0 failed;
`cargo test --lib ng::parameter_estimation::joint` 103 passed, 0 failed, 501.66s. Full
suite figures are in the commit message.

## 6. Tradeoffs and follow-ups

- **The partial decision above** — Checkpoint A.
- **`GenomeRegion::len()` overflows in debug at the coordinate ceiling** while its doc
  says it saturates. Not touched here; worth its own small fix.
- Nothing calls `reach` or `non_reference_reads` yet. A3 walks on the first, A4 judges
  on the second.
