# ng step 6 at a repeat tract — implementation plan


> **⚠ 2026-08-26 — the repeat tract's genotype prior no longer indexes its seed by the
> cohort's modal repeat count**, so every sentence below that justifies carrying the mode
> *for the prior* has lost its reason. The seed is now the fitted **length spectrum** the
> joint repeat fit produces per stratum, indexed by whole-repeat offset from the
> **reference** tract length, which every locus already knows
> ([`../spec/population_diversity.md`](../spec/population_diversity.md) §4.2; built by step
> E2e of [`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md)). Handing the mode
> where the reference length belongs is now a measured error: it moves 0.595 of the prior's
> mass off the reference length onto 0.091 on `seed_ssr`'s own fixture. **Whether selection
> should still carry the mode is undecided and is this document's to settle** — it may have
> uses of its own; what it no longer has is this one.
*Draft, 2026-08-24. Turns the settled design —
[`../spec/candidate_alleles_ssr.md`](../spec/candidate_alleles_ssr.md) and
[`../arch/candidate_alleles_ssr.md`](../arch/candidate_alleles_ssr.md) — into build order. **Not a
place for new design.***

*Follows [`candidate_alleles.md`](candidate_alleles.md), which builds the shared module this one
extends. Its first milestone is a change to the merge, because the repeat-tract path cannot be
built without a field `CohortObservation` does not carry.*

---

## Scope

**In:**

- `src/ng/run/cohort_merge/close.rs` and `build.rs` — `ClosedLocus` carries the locus kind, and
  `CohortObservation` gains `pub kind: LocusKind` (Milestone A). Without it there is no motif, so
  no period, so no repeat count and no ladder.
- `src/ng/calling/allele_candidates/ssr.rs` — `RepeatLadder`, `SsrSelectionConfig` with
  `DEFAULT_MAX_OFF_GRID_SHARE`, `SsrLocusSelection`, `select_ssr`.
- `src/ng/calling/allele_candidates/mod.rs` — `SelectionScratch` gains the ladder's buffers.
- `examples/ng_candidate_selection_probe.rs` — a repeat-tract arm, so the ladder is exercised on
  real reads through the merge (Milestone E).

**Out (later plans):**

- **The ordinary path and everything shared** — [`candidate_alleles.md`](candidate_alleles.md),
  which must have merged first.
- **How a rung's prior weight divides between two spellings** —
  [`calling_prior.md`](calling_prior.md)'s follow-on, against
  [`../spec/calling_priors.md`](../spec/calling_priors.md) §5.2. Spec §5 hands it the fact it was
  waiting for and a leaning; this plan produces `repeat_counts` and stops there.
- **The merge's variability filter at a tract too long for a read to span** — already booked in
  the merge's own code ([`build.rs:1025-1034`](../../../../src/ng/run/cohort_merge/build.rs)) and
  left there. Spec §8.
- **Wiring the STR generator into a production run's `GeneratorSet`** — nothing in `src/` builds
  one with the `Ssr` slot filled today; only tests and examples do. Milestone E fills it in the
  probe, which is all this plan needs. A real run's wiring belongs with
  [`calling_loop.md`](calling_loop.md)'s builder wiring.

## Principles (how the order was chosen)

- **The blocker first, and alone.** Milestone A is a change to a merged, tested module and touches
  no calling code; it lands and is reviewed on its own so that a defect in it cannot be confused
  with a defect in the selector.
- **Types first, then implementation**, within each milestone (project rule).
- **The algorithmic heart before the plumbing.** The ladder and nomination are built and proven on
  hand-built loci before anything assembles an `SsrLocusSelection`.
- **Reuse over rewrite.** The `±1` rescue and its `occupied` test are ported unchanged; the shared
  support rule is called, not rewritten (arch §5's table names every row).
- **Verify against ground truth, and here there is a real oracle at both ends.** Production's
  selector is the differential: with its three replaced rules switched in, this one must reproduce
  production's candidate set on the tomato panel; with them switched out, it must move by the
  amounts spec §4.1 and §5 measured on HG002 against the assembly truth set. **Unlike the ordinary
  path, this half has a failing state at both ends** (spec §10).
- **Isolate the silent steps.** Three steps produce a quietly wrong answer rather than a crash —
  **B1** the rung keying, **C2** the `±1` rescue, **D2** periodicity's units — and each lands as
  its own commit with its oracle green before and after.
- **Container builds.** All `cargo` via `./scripts/dev.sh`; a native host build at completion.

## Preconditions

- **[`candidate_alleles.md`](candidate_alleles.md) merged**, all four milestones: the shared
  config and constants, `SelectionVerdict`, `UnmatchedSupport`, `AlleleRemap`, `LocusSelection`,
  `SelectionScratch`, and the fold and ranking `select_ssr` reuses.
- **The merge, complete** ([`cohort_merge.md`](cohort_merge.md)), including the closer's exhaustive
  branch on the locus kind ([`close.rs:110-125`](../../../../src/ng/run/cohort_merge/close.rs)),
  which is what makes Milestone A a carry-through rather than a new derivation.
- **The STR locus generator**: `SsrGenerator` and `SsrDetail` with its `Motif`
  ([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs);
  `Motif::period` at [`types.rs:1138`](../../../../src/ng/types.rs)).
- **The measurement rig**: `benchmarks/ssr_hg002/` with its Tier catalog, its 5×–300× read ladder
  and its assembly-derived truth set; production's `ssr-pileup` and
  `examples/ssr_slip_dump`, which produced every number spec §4.1 and §5 quote.
- **Production's selector as the differential's other end**:
  [`candidate_set.rs`](../../../../src/ssr/cohort/candidate_set.rs) and
  [`rung_ladder.rs`](../../../../src/ssr/cohort/rung_ladder.rs), read-only — `src/ssr/` is frozen.

## Branch and merge

> **2026-09-02:** superseded for execution — Milestones B–E now run on the STR loop plan's
> branch and worktree ([`calling_loop_ssr.md`](calling_loop_ssr.md) §Branch), in parallel
> with the observations plan that delivers Milestone A.

- **Branch** `ng-candidate-alleles-ssr`, from `main` after `ng-candidate-alleles` has merged. No
  worktree.
- **Conflict surface:** `calling/allele_candidates/mod.rs` (the scratch gains fields) and the two
  merge files of Milestone A. Nothing else in flight edits them.

---

## The steps

### Milestone A — the merge carries the locus kind

> **2026-09-02 — ownership note.** This milestone is now delivered by the observations plan
> ([`run_ssr_observations.md`](run_ssr_observations.md) Milestone A, against
> [`../spec/run_ssr_observations.md`](../spec/run_ssr_observations.md) §4), which runs in
> parallel with the STR loop plan and merges this change to `main` first. Do not build it
> twice: when that checkpoint has merged, flip A1 here and start at Milestone B.

**A1. `ClosedLocus` and `CohortObservation` carry `LocusKind`.**  ☐
The closer already has the kind in scope where it builds a `ClosedLocus`
([`close.rs:713-721`](../../../../src/ng/run/cohort_merge/close.rs)) and drops it; carry it
through, and have `CohortObservation::over` copy it onto the assembled locus. **The clone is one
per *built* locus, not per closed one**, and `LocusKind::Ssr` boxes two flanks — say so in the doc
comment, per the repo's clone-audit rule, and note that the generic variant clones nothing.
Nothing else changes: no verdict moves, no test's expectation changes except the constructors that
now name a kind.
*Source:* arch §1; spec §2. *Depends:* none.

> **Checkpoint A:** the merge's own tests are unchanged in outcome and the field is populated on
> both kinds. **Reviewed and merged before any calling code is written**, so that a defect here
> cannot be confused with a defect in the selector. Pause for review.

### Milestone B — the ladder

**B1. `RepeatLadder`, keyed by repeat count. — own commit, do not bundle.**  ☐
Built once per locus from the merge's allele table and `SsrDetail::motif`: each rung holds the
merge-table indices of the sequences at that repeat count, ascending, and the ladder holds the
cohort's modal repeat count — most-supported rung, ties to the shorter, and **not the reference's
count**. The key is `bases.len() / motif.period()`, floored, and both offsets are measured in
units, never bases.
**Why this one is isolated:** the genotype prior indexes its seed by the same integer
([`../spec/calling_priors.md`](../spec/calling_priors.md) §5), so a keying that disagrees puts a
candidate's prior mass on the wrong rung and every genotype at the locus shifts, silently.
*Oracle:* two sequences of one length land on one rung and stay distinct inside it; a sequence
whose length is not a whole number of units lands on the floored rung and is counted once; the
mode is the rung with the most reads and the shorter of a tie.
*Depends:* A1. *Source:* arch §2.1; spec §3.

**B2. `SsrSelectionConfig` and the per-sample length histogram.**  ☐
The config — the shared `CandidateSelectionConfig`, `max_off_grid_share` with
`DEFAULT_MAX_OFF_GRID_SHARE = 0.10` marked inherited-and-never-measured, and `ploidy` taken from
the caller rather than a constant — and the per-sample fold from the merge's rows into a length
histogram over spanning reads. `SelectionScratch` gains both buffers.
*Depends:* B1. *Source:* arch §2.2, §3.1; spec §7, Q1.

> **Checkpoint B:** the ladder and the histograms are proven on hand-built loci, with nothing
> nominated yet. Pause for review.

### Milestone C — nomination

**C1. The per-sample bar over rungs, and top-`ploidy`.**  ☐
A repeat count is nominated when the sample's reads at it clear the **shared** support rule against
that sample's spanning reads — `MinAltReads::reached_by` with the rung's read total as the
numerator, not a second predicate. The top `ploidy` by support are promoted, ties to the shorter
length. Production's `is_clear_peak` is not called and not ported.
*Depends:* B2. *Source:* arch §3.1; spec §4.

**C2. The `±1` rescue and the union. — own commit, do not bundle.**  ☐
When a sample resolved fewer than `ploidy` counts, promote each promoted count's `±1` neighbours
**where some sample's reads reached that length** — production's `occupied` test
([`candidate_set.rs:221`](../../../../src/ssr/cohort/candidate_set.rs)), ported unchanged. The
cohort's promoted set is the union across samples.
**Why this one is isolated:** dropping the `occupied` test invents a length, and firing the rescue
on a sample that did resolve its ploidy widens every locus — both are extra candidates, not
crashes.
*Oracle:* a sample resolving one length nominates its occupied neighbours and not its unoccupied
ones; a sample resolving two nominates no neighbour at all; **and the test spec §13 names as the
one production cannot pass — a sample with 150 reads at count 10 and 150 at 11 nominates both.**
*Depends:* C1. *Source:* arch §3.1; spec §4, §4.1.

> **Checkpoint C:** nomination is complete, and the adjacent-length heterozygote test is green.
> Pause for review.

### Milestone D — admission, periodicity, and the outputs

**D1. Sequence admission within a promoted rung.**  ☐
Every sequence on a promoted rung faces the shared support rule, asked of the sequence: no
privileged representative, no recurrence term. Then the shared cap and truncation, the shared
leftover, and the reference tract admitted first and exempt from both.
*Depends:* C2. *Source:* arch §3.2; spec §5.

**D2. Periodicity, per sample. — own commit, do not bundle.**  ☐
A sample is non-periodic when more than `max_off_grid_share` of its spanning reads sit at a length
that is not a whole number of motif units from the ladder's mode; the locus is `NotPeriodic` only
when **no** sample is periodic. A `NotPeriodic` locus returns the reference tract alone, that
verdict, an empty leftover and a one-entry `repeat_counts`.
**Why this one is isolated:** production measures this offset in bases while keying its ladder in
units ([`candidate_set.rs:114-145`](../../../../src/ssr/cohort/candidate_set.rs) against
[`rung_ladder.rs:301`](../../../../src/ssr/cohort/rung_ladder.rs)); repeating that mismatch
refuses periodic loci and admits non-periodic ones, and nothing fails.
*Oracle:* a locus whose reads are all a whole number of units from the mode is periodic whatever
their lengths in bases; one sample being periodic saves a locus every other sample fails.
*Depends:* D1. *Source:* arch §3.3; spec §7, §3.

**D3. What the prior takes: `repeat_counts` and `modal_repeat_count`.**  ☐
`SsrLocusSelection` returns both, parallel to the candidates with the reference at index 0.
`fill_ssr_seed` takes exactly this slice and must not recompute it — one producer for an integer
two modules have to agree on.
*Depends:* D2. *Source:* arch §2.3; spec §3.

> **Checkpoint D:** `select_ssr` is complete and proven on hand-built loci, including the
> `NotPeriodic` and reference-only outcomes. Pause for review.

### Milestone E — the two ends of the differential, on real data

**E1. Production's rules as a test-only arm.**  ☐
In `#[cfg(test)]`, a re-implementation of production's three replaced rules — the clear-peak test,
the cohort-summed depth gate, the same-length sibling bar — driving the same fold. **Not a field of
`SsrSelectionConfig`:** the shipping binary carries one rule, and a configuration nobody should set
does not belong in it (arch, *Test & bench shape*). With the arm switched in, `select_ssr` must
reproduce production's candidate set on a committed tomato fixture.
*Depends:* D3. *Source:* arch *Test & bench shape*; spec §10.

**E2. The probe's repeat-tract arm, and the HG002 numbers.**  ☐
Give `examples/ng_candidate_selection_probe.rs` a `GeneratorSet` with the `Ssr` slot filled, as
`examples/ng_ssr_loci_dump.rs` already does, so repeat tracts reach the merge and `select_ssr`. Then
re-run spec §4.1's and §5's scoring against `benchmarks/ssr_hg002/truth/` **through the shipped
module** rather than through the offline scorers that produced those tables. The numbers to
reproduce, at 300×: both alleles of an adjacent-length heterozygote offered at about 97%, both
spellings of one length at about 86%, and about 1.26 candidate sequences per tract. A difference is
a defect and must be traced.
*Depends:* E1. *Source:* spec §4.1, §5, §13.

> **Checkpoint E:** both ends of the differential green — production reproduced with its rules
> switched in, and the measured improvement reproduced with them switched out. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | the merge's existing tests unchanged in outcome; the kind populated on both variants |
| B | hand-built loci: one rung for two same-length sequences, kept distinct; the floored key; the mode as most-supported-then-shorter |
| C | the occupied-neighbour rescue in both directions; **the adjacent-length heterozygote nominating both** — the case production cannot |
| D | the sequence bar with no privileged representative; `NotPeriodic` per sample with one periodic sample saving the locus; `repeat_counts` parallel to the candidates |
| E | **external, both ends**: production's candidate set reproduced with its rules switched in on a tomato fixture; spec §4.1's and §5's HG002 recall reproduced with them switched out |

**Property tests from Milestone D on**, inherited from the shared plan and re-run on this path:
the surviving list is a subset of the merge's table, contains the reference, never exceeds the cap;
adding a sample that shows only the reference tract changes neither the candidate set nor any other
sample's nomination; and every rule gives the same answer for one sample alone that it gives for
that sample inside a cohort of sixty-three.

## Out of scope (next plans)

- **Dividing a rung's prior weight between two spellings** —
  [`../spec/calling_priors.md`](../spec/calling_priors.md) §5.2 and its plan.
- **Discovery at a repeat tract, and the cap's refuse-don't-evict rule** —
  [`calling_bakeoffs.md`](calling_bakeoffs.md), against
  [`../spec/candidate_alleles.md`](../spec/candidate_alleles.md) §6.3.
- **Wiring the STR generator into a production run** — [`calling_loop.md`](calling_loop.md)'s
  builder wiring; Milestone E fills it in the probe only.
- **The three measurements that would move the soft numbers** — spec Q1 (the periodicity threshold,
  never measured by anyone), Q2 (the cap at a tract in a large cohort), Q3 (the whole path is one
  human sample). All are parameters; no code waits on them.
