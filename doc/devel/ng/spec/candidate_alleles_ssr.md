# ng step 6 at a repeat tract — choosing the tract sequences a locus is called over


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
*Design spec, 2026-08-24. **No code yet — this settles the design.** Second of two documents on
candidate selection. It inherits everything in
[`candidate_alleles.md`](candidate_alleles.md) — the support rule, the cap and its ranking, the
verdict, the reference at index 0, where it runs — and replaces only what a repeat tract needs
replaced. **Read that one first**; this document does not restate it. The types and signatures
are [`../arch/candidate_alleles_ssr.md`](../arch/candidate_alleles_ssr.md).*

*Reads on: [`locus_generation_ssr.md`](locus_generation_ssr.md) — the observations this narrows;
[`cohort_merge.md`](cohort_merge.md) §4.2. Read by: [`calling_priors.md`](calling_priors.md) §5,
which lays the prior's mass over the same ladder this builds, and
[`read_likelihoods.md`](read_likelihoods.md) §4.5, whose junk term is spread over the lengths this
candidate set can reach.*

*Production's equivalent is
[`ssr/cohort/candidate_set.rs`](../../../../src/ssr/cohort/candidate_set.rs) with
[`rung_ladder.rs`](../../../../src/ssr/cohort/rung_ladder.rs). Everything said about them is a
record of what they do, not a proposal to change them — `src/ssr/` is frozen production.*

---

## 1. What this is

**At an ordinary locus the merge has already unified every sample's observations into one flat
table, and selection is a bar and a cap over it. At a repeat tract that is not enough, because the
alleles are not unordered.** A tract of 11 repeats is adjacent to one of 10 and one of 12 and far
from one of 4, and both the genotype prior and the stutter model are written on that ordering.
Selection here therefore has structure: sequences are grouped by how many repeats they carry, each
sample nominates which groups to promote, and the sequences inside a promoted group each face the
support rule.

### 1.1 Goals, beyond the ones inherited

1. **Produce a candidate set the genotype prior can index.** The prior's mass falls off
   geometrically with distance from the cohort's modal repeat count
   ([`calling_priors.md`](calling_priors.md) §5.1), so the grouping by repeat count is the prior's
   own index and not an implementation choice.
2. **Offer both spellings of one repeat length where the reads show both.** An interrupted repeat
   — two tract sequences of the same length differing by an interior base — is a real allele class
   and the read likelihood already separates the two by about 28 Phred per distinguishing base
   ([`read_likelihoods.md`](read_likelihoods.md) §4.6).

### 1.2 Non-goals

Everything in [`candidate_alleles.md`](candidate_alleles.md) §1.2, plus: **this does not model
stutter.** Which tract lengths a polymerase can slip to is the read likelihood's
([`read_likelihoods.md`](read_likelihoods.md) §4.2); selection only decides which observed
sequences are worth scoring.

### 1.3 Vocabulary

- **tract sequence** — the bases a read showed across the repeat tract, as the locus generator
  minted them ([`SequenceObservation`](../../../../src/ng/locus_generation/mod.rs)). An ng repeat
  allele **is** a sequence, never a repeat count: two alleles of the same length can differ inside
  ([`locus_generation_ssr.md`](locus_generation_ssr.md) §3).
- **rung** — all the observed tract sequences carrying the same number of repeats, grouped. The
  set of occupied rungs at a locus is its **ladder**. Production's word, and it is used here
  because [`calling_priors.md`](calling_priors.md) §5.1 already writes the prior on it.
- **spanning reads** — a sample's reads that crossed the whole tract and so name a length. A read
  that ran out inside it does not, and is held apart (§8).

---

## 2. What this needs that does not exist yet

**The merge closes repeat-tract loci but the assembled observation cannot key a ladder.** Closing
handles them deliberately — a tract is exempt from the span cap, and the exemption is exhaustive
on the locus kind so a new kind is a compile error rather than a silent default
([`close.rs:110-125`](../../../../src/ng/run/cohort_merge/close.rs)). But the assembled
`CohortObservation` carries only the region, the allele table and the per-sample rows
([`build.rs:922-929`](../../../../src/ng/run/cohort_merge/build.rs)): **it drops the locus kind,
and with it the motif.** Without the motif there is no period, and without a period there are no
repeat counts and no ladder.

**That is one field on the merge's own type and it is the merge's to add**, not this step's. It is
stated here because it is the scheduling fact: the ordinary path of
[`candidate_alleles.md`](candidate_alleles.md) is buildable today and this document is not.

**A second thing is already in place and worth knowing**, because an earlier reading of this
module said otherwise: the merge's per-sample rows are keyed on `(allele, read group)`, not pooled
([`SupportedAllele`, `build.rs:1089`](../../../../src/ng/run/cohort_merge/build.rs)). The stutter
model is fitted per read group, so the STR emission needs them apart, and it has them. **The
support rule pools them back** — it counts reads, and a read is a read whichever lane it came from
— which is deliberate and is the one place pooling is correct
([`SampleSupport::pooled_support_for`, `:1070`](../../../../src/ng/run/cohort_merge/build.rs)).

---

## 3. The ladder — ported, and why

**A rung is keyed by repeat count: the tract sequence's length in bases divided by the motif's,
rounded down.** Production keys it exactly so
([`rung_ladder.rs:291-320`](../../../../src/ssr/cohort/rung_ladder.rs)) and ng keeps it, because
the prior's geometric seed is indexed by the same quantity and the two must land on the same
integer ([`calling_priors.md`](calling_priors.md) §5.1). **Inside a rung, sequences stay separate
by exact bytes** — the grouping is for nomination and for the prior, never a merge of evidence.

**One inconsistency in production is not ported.** Its periodicity gate measures a read's tract in
**bases** while the ladder keys in **units**
([`candidate_set.rs:114-145`](../../../../src/ssr/cohort/candidate_set.rs) against
[`rung_ladder.rs:301`](../../../../src/ssr/cohort/rung_ladder.rs)), so a 13-base read at a
dinucleotide is judged out of frame by the gate and then pooled onto rung 6 with the 12-base
reads, where it counts toward that rung's totals. ng measures both in units. Recorded as a fact
about the code, not as a criticism of a caller that works.

---

## 4. Which rungs a sample nominates

**Decision: a sample nominates every repeat count whose spanning reads clear the inherited support
rule, and the two best-supported of those are promoted — with the `±1` neighbours added when it
resolved fewer than its ploidy.**

```text
a sample s nominates repeat count L  ⟺
        spanning_reads_s(L) ≥ max( 2,  ceil( 0.05 × spanning_reads_s ) )

promote the top `ploidy` by support; if fewer than `ploidy` cleared it,
also promote each nominated L's occupied ±1 neighbours
```

The `±1` rescue is production's and is kept for its reason: a sample whose two copies the reads
could not resolve is more likely to carry a neighbour of what it did resolve than anything else,
and "occupied" keeps it honest — a neighbour is added only where some sample's reads actually
reached that length
([`candidate_set.rs:221`](../../../../src/ssr/cohort/candidate_set.rs)). **Nothing here invents a
length.**

### 4.1 What is replaced, and what it costs — measured

**Production nominates a repeat count only if its read count exceeds *both* neighbouring counts by
more than three reads** (`is_clear_peak`,
[`rung_ladder.rs:274-288`](../../../../src/ssr/cohort/rung_ladder.rs), `prominence = 3` at
[`candidate_set.rs:83`](../../../../src/ssr/cohort/candidate_set.rs)). Two consequences follow
from the predicate itself, before any data:

- **A sample heterozygous for two adjacent repeat counts has no clear peak.** With the two lengths
  at similar read counts, neither exceeds the other by three. The rescue does not save it, because
  it iterates over the peaks that *were* found and there are none
  ([`candidate_set.rs:239-258`](../../../../src/ssr/cohort/candidate_set.rs)) — so such a sample
  nominates nothing at all and the locus falls back to whatever other samples and the reference
  supply.
- **It is an absolute read count with no depth term**, so it binds only at the shallow end: a
  sample needs four reads at one length with both neighbours empty to nominate anything, and at
  300 reads a clear peak is never in doubt.

**And a one-repeat difference is the commonest heterozygote at a tract**, so this is not a corner.

*Measured 2026-08-24 on GIAB HG002 (`benchmarks/ssr_hg002/`), real reads through production's own
Stage-1 pileup over the 13,272-tract Tier catalog, scored against the assembly-derived truth
genotypes. **Of HG002's heterozygous tracts, 51% carry their two copies one repeat apart.** The
table gives the fraction of tracts where **both** true repeat counts are offered as candidates.*

| HG002 coverage | tracts scored | median spanning reads | production's rule | the rule of §4 |
|---|---|---|---|---|
| 5× | 179 | 4 | 33% | 100% |
| 10× | 3,858 | 8 | 35% | 98% |
| 20× | 11,365 | 16 | 45% | 98% |
| 30× | 12,593 | 23 | 50% | 98% |
| 50× | 12,846 | 39 | 60% | 98% |
| 300× | 12,905 | 229 | 78% | 97% |

*The 5× row rests on 179 tracts — the rest were refused by the depth gate of §6 — and is
indicative only. Homozygous tracts score 99.9–100% under both rules, which is the control: a
broken truth join would fail those too, and the rule fails only where the mechanism above predicts
it will.*

**The replacement is cheaper as well as better, which is unusual enough to state.** At 300×
production offers 1.65 candidate repeat counts per tract of which 0.52 are neither true allele nor
reference; the rule of §4 offers 1.28 of which 0.15 are. At 30× the two are level on cost — 1.28
candidates each — and 48 points apart on recall. **This is not a trade.**

**Why tomato could not have found this, and it matters because tomato is where this project's
repeat-tract work has been done.** The panel is an extreme selfer — apparent `F_IS` about 0.82
([`calling_priors.md`](calling_priors.md) §5.1) — so heterozygotes of any kind are rare and
adjacent-length ones rarer still. **The one measurement that looks like it contradicts §4.1 was
made there**: production's own drop attribution, on the tomato panel in July 2026, put "the
variant was never nominated as a candidate" at 3 of 108 confident real misses
([`ssr_emission_drop_attribution_2026-07-08.md`](../../reports/ssr_emission_drop_attribution_2026-07-08.md)).
**That 3% is a fact about a selfing panel, not about the mechanism.** HG002 is outbred and about
72 tandem repeats in 100 are heterozygous there
([`calling_priors.md`](calling_priors.md) §5.3).

---

## 5. Two spellings of one repeat length

**Decision: a second sequence on an occupied rung faces the same support rule as everything else —
no representative is privileged and no recurrence bar applies.**

```text
a tract sequence q survives  ⟺  ∃ sample s :
        reads_s(q) ≥ max( 2, ceil( 0.05 × spanning_reads_s ) )
```

That is the inherited rule of [`candidate_alleles.md`](candidate_alleles.md) §3 asked of the
*sequence* rather than of the rung. The rungs decide which lengths are in play (§4); within a
promoted length, every spelling stands on its own reads.

**Production does something different and it cannot work below three samples.** It promotes the
rung's best-supported sequence unconditionally and makes any sibling clear all three of: 8 reads,
**3 distinct samples**, and 10% of that rung's reads
([`candidate_set.rs:169-191`](../../../../src/ssr/cohort/candidate_set.rs), constants at
[`:86-88`](../../../../src/ssr/cohort/candidate_set.rs)). The three-sample term is an absolute
constant with no cohort-size clamp, so **at one and two samples no second spelling can ever be
promoted** — the mechanism is simply absent. Production clamps the *other* recurrence constant for
small cohorts and not this one
([`driver.rs:368`](../../../../src/ssr/cohort/driver.rs)), so the omission is visible in its own
code.

*Measured 2026-08-24, same HG002 runs, scored at **sequence** level: each haplotype's tract
sequence reconstructed from the reference tract plus the phased truth VCF's edits inside it, then
looked for among the candidates each rule offers. "Shown by some read" is the ceiling the evidence
allows — a rule cannot offer a sequence no read carried.*

At 300×, of 11,978 tracts scored:

| class | tracts | production | the rule of §5 | shown by some read |
|---|---|---|---|---|
| homozygous | 11,283 | 99.8% | 99.8% | 99.9% |
| heterozygous, different repeat length | 399 | 72.4% | 86.0% | 90.0% |
| **heterozygous, same repeat length, different spelling** | **296** | **35.1%** | **86.1%** | **93.6%** |

**The same-length class is 43% of HG002's heterozygous tracts — 296 of 695 — and production offers
both spellings at about a third of them.** The third it does reach are the ones where one copy
happens to equal the reference tract, which is seeded free; the other copy is covered by the rung's
representative, and the representative is the only sequence per rung a single-sample run can get.
Production's figure is flat across the whole coverage ladder — 33% at 10×, 38% at 30×, 35% at 300×
— because the limit is structural rather than evidential.

**The rule of §5 costs fewer candidates, not more:** 1.26 sequences per tract at 300× against
production's 1.57.

**Why the share is 5% here too.** Sweeping it at 300×: 2% reaches 88.5% on the same-length class
at 1.46 candidates per tract, 5% reaches 86.1% at 1.26, 10% reaches 85.8% at 1.22, 20% reaches
84.5% at 1.20. **The knee is between 2% and 5%, and 5% is chosen so that one number governs both
paths** ([`candidate_alleles.md`](candidate_alleles.md) §3.3 settles the same value on the
ordinary path). It costs 2.4 points of same-length recall against 2% at 300× and nothing at all at
30× and below, where the floor binds.

**This answers a question the genotype prior was holding open.** [`calling_priors.md`](calling_priors.md)
§5.2 asks how a rung's prior weight divides between two alleles that sit on it, and could not be
answered while it was unknown whether both are ever selected. They are, at about 86% of the tracts
where the reads allow it. **Leaning, and it is that document's decision, not this one's: divide the
rung's weight evenly**, because the support rule has already established both are real and
weighting the division by reads would count the same evidence the likelihood is already using.

---

## 6. The gates that are not ported

**Two of production's three locus gates read the cohort, and one is measured badly broken at a
single sample.** [`candidate_alleles.md`](candidate_alleles.md) §3.2 makes reading the cohort
disqualifying; this section is what that costs in ported code.

**The depth gate is dropped entirely.** `min_cohort_depth = 10` is a **sum over the cohort**
([`cohort_depth`, `candidate_set.rs:94-100`](../../../../src/ssr/cohort/candidate_set.rs),
constant at [`:84`](../../../../src/ssr/cohort/candidate_set.rs)), so the same tract is refused
alone and admitted in company. *Measured on HG002 as a single sample, 13,272 tracts:*

| coverage | median spanning reads | tracts refused as `LowDepth` |
|---|---|---|
| 5× | 4 | 13,091 (98.6%) |
| 10× | 8 | 9,359 (70.5%) |
| 20× | 16 | 1,649 (12.4%) |
| 30× | 23 | 369 (2.8%) |
| 50× | 39 | 90 (0.7%) |
| 300× | 229 | 22 (0.2%) |

At 63 tomato accessions at 3 reads a position the same gate passes almost everything, because the
sum is about 190. **Depth is asked once, upstream, per sample, by the merge's keep rule**
([`cohort_merge.md`](cohort_merge.md) §4.3), and asking it again here with a cohort denominator is
the drift the merge itself retired in August. There is no depth verdict
([`candidate_alleles.md`](candidate_alleles.md) §6.2).

**The three sibling constants are dropped**, replaced by §5's rule: `min_same_length_reads = 8`,
`min_same_length_samples = 3` and `min_same_length_fraction = 0.10`. Two of the three read the
cohort and the middle one is dead below three samples.

**The cap of 24 and its refusal are dropped**, replaced by a cap of 32 with
truncation ([`candidate_alleles.md`](candidate_alleles.md) §4). Production's is checked against a
set that grows monotonically with cohort size and refuses the locus outright
([`candidate_set.rs:272-276`](../../../../src/ssr/cohort/candidate_set.rs)), so at a large panel it
refuses loci for being polymorphic rather than for being bad.

---

## 7. The periodicity verdict

> **⚠ 2026-09-02 — the grid is anchored on the reference tract's length, not on the ladder's
> mode.** Owner's decision, taken at step D2 of
> [`../impl_plan/candidate_alleles_ssr.md`](../impl_plan/candidate_alleles_ssr.md). Read
> literally, "a whole number of motif units from the ladder's mode" anchors the grid at **zero**
> — the mode is a repeat count, so its length in bases is a whole number of units and cancels out
> of the subtraction — and a zero-anchored grid refuses a real class of tract. The catalog trims
> every tract back to whole motif copies **at both ends**, but a length-changing interruption
> inside puts the two ends out of phase, so the tract's own reference length is then not a
> multiple of the period; such a tract is admitted whenever the break is late enough to clear the
> purity floor of 0.8. Measured through the catalog's own `minimal_trim` and `recompute_purity`:
> 49 bases of an `AT` repeat with one extra base 40 bases in trims to 49 bases at **purity
> 0.816**, and 49 is odd — so every read at that tract's reference length is off a zero-anchored
> grid and the locus is refused. Production avoids this by anchoring on the commonest observed
> length in bases; the reference tract's length is preferred to that because it is a property of
> the locus rather than of the reads, so it cannot move with depth, and because it is the
> quantity the genotype prior was re-indexed onto on 2026-08-27.

**`NotPeriodic` survives, and it is the one verdict this path adds.** A stretch the catalog called
a repeat tract whose reads do not actually vary in whole motif units is not a tract this caller's
model describes: the stutter distribution is written on whole-repeat and part-repeat regimes
([`read_likelihoods.md`](read_likelihoods.md) §4.2) and the prior's ladder is written on repeat
counts. Genotyping it against that model would produce a confident answer from a model that does
not apply.

**The measure is production's, with its denominator corrected to be per sample.** Production asks
what fraction of the **cohort's** reads sit at a length that is not a whole number of motif units
from the modal length, and rejects above `max_out_of_frame_frac = 0.10`
([`candidate_set.rs:114-145`](../../../../src/ssr/cohort/candidate_set.rs), constant at
[`:85`](../../../../src/ssr/cohort/candidate_set.rs)). **The 10% is inherited and unmeasured**;
what changes is that the fraction is taken of each sample's own spanning reads and the locus is
`NotPeriodic` only when **no** sample is periodic by it — the same "one sample suffices" shape as
everything else here, for the same reason.

**A `NotPeriodic` locus still yields a candidate table**: the reference tract alone, verdict
`NotPeriodic`. What the run does with it is emission's.

---

## 8. What this path does *not* owe, and what it decides for others

**No leftover pool.** [`candidate_alleles.md`](candidate_alleles.md) §5 exists because the
SNP/indel likelihood needs the error mass of reads matching no candidate. **This path needs none:**
a read no allele explains is already carried by the junk term, spread uniformly over the tract
lengths the stutter model can reach from the candidate set
([`read_likelihoods.md`](read_likelihoods.md) §4.5, and `reachable_length_count` in
[`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §4.1). A read at an unselected length
lands there by construction.

**But the candidate set is that term's denominator, which is a coupling worth naming.** Widen the
set and the junk floor per length falls; narrow it and the floor rises. This is why a candidate may
not be added between two passes of the frequency loop
([`calling_em_loop.md`](calling_em_loop.md) §4's table) and why a discovery round that would
overflow the cap is refused rather than allowed to evict
([`candidate_alleles.md`](candidate_alleles.md) §6.3).

**Partial reads are not read here.** A read that ran out inside a tract gives a lower bound on the
length, not a length, and the merge holds it on its own axis
([`PartialObservation`, `build.rs:1130`](../../../../src/ng/run/cohort_merge/build.rs)). It does
not nominate a rung and does not count toward any denominator, which is the same treatment
`reads_compared_with_reference` already gives it upstream. **One consequence belongs to the merge,
not here, and is already booked there:** a sample carrying an allele too long for a read to span
shows no complete observation, so the merge's variability filter reads it as *nothing varied*, and
"one line of the rule has to change for that to mean what it says, and only on this path"
([`build.rs:1025-1034`](../../../../src/ng/run/cohort_merge/build.rs)).

---

## 9. One sample and a thousand, three reads and three hundred

**At one sample everything in §4 and §5 still works**, which is the point of removing the three
cohort-denominated constants: the ladder is that sample's histogram, its two best-supported
lengths are promoted, and every spelling on them faces the same rule. §5's table is a
single-sample measurement.

**At three reads a position the ladder is nearly empty and the reference carries the locus.** A
sample with 3 spanning reads clears the rule at one length at most, so it nominates at most one
rung; the candidate set is the reference plus that. This is the same collapse
[`candidate_alleles.md`](candidate_alleles.md) §7 describes and it has the same answer — the union
across samples is what recovers a length one accession could not resolve. On the tomato panel the
share is inert and the rule is a count of 2.

**At three hundred reads the share does the work.** §5's sweep is that measurement: at 300× moving
the share from 2% to 5% cuts candidates per tract from 1.46 to 1.26.

**At a thousand samples the ladder widens with the cohort**, as the ordinary path's table does
([`candidate_alleles.md`](candidate_alleles.md) §4.2), and the cap of six is what bounds it. **A
repeat tract is where that cap is most likely to bind**, because a tract genuinely carries more
alleles than a SNP does — this is the one place the extrapolation of that document's Q2 is most
likely to come true, and §12's Q2 here records it as the same open question seen from this side.

---

## 10. Reuse map

| what | existing code | how ng reuses it |
|---|---|---|
| the rung ladder — keying by repeat count, sequences kept apart inside a rung | [`build_rungs`, `rung_ladder.rs:291`](../../../../src/ssr/cohort/rung_ladder.rs) | **ported**, with the frame inconsistency of §3 corrected |
| per-sample nomination and top-`ploidy` selection | [`assemble_candidates`, `candidate_set.rs:194`](../../../../src/ssr/cohort/candidate_set.rs) | **shape ported, predicate replaced** (§4.1) |
| the `±1` neighbour rescue and its `occupied` test | [`candidate_set.rs:239-258`](../../../../src/ssr/cohort/candidate_set.rs), [`:221`](../../../../src/ssr/cohort/candidate_set.rs) | **ported unchanged** |
| the union across samples, deduplicated by exact bytes | [`candidate_set.rs:223-227`](../../../../src/ssr/cohort/candidate_set.rs) | ported; the linear `contains` scan is replaced by a set, since the cap is now checked after a bar rather than after an unbounded union |
| the periodicity measure | [`is_periodic`, `candidate_set.rs:114`](../../../../src/ssr/cohort/candidate_set.rs) | shape ported, denominator moved per sample (§7); the 10% is inherited and unmeasured |
| the reference tract seeded as allele 0 before any gate | [`candidate_set.rs:200-202`](../../../../src/ssr/cohort/candidate_set.rs) | ported; it is the inherited invariant ([`candidate_alleles.md`](candidate_alleles.md) §6.1) |
| the depth gate, the three sibling constants, the cap of 24 | [`candidate_set.rs:80-90`](../../../../src/ssr/cohort/candidate_set.rs) | **not ported** (§6) |

**The oracle is a differential, not a parity test, and that is deliberate.** Three of production's
rules are being replaced on purpose, so byte-identical output is impossible by construction and a
parity test with an escape clause has no failing state. Instead: build the selector so the peak
test, the depth gate and the sibling constants are switchable; with them switched **in**, require
it to reproduce production's candidate set on the tomato panel exactly; then switch them **out**
and report what moved, against the numbers in §4.1 and §5. That has a failing state at both ends.

---

## 11. Deferred, with a recommended home

- **How a rung's prior weight divides between two spellings** — to
  [`calling_priors.md`](calling_priors.md) §5.2, whose question it is. §5 gives it the fact it was
  waiting for and a leaning.
- **The merge's variability filter at a tract too long for a read to span** — already booked at
  [`build.rs:1025-1034`](../../../../src/ng/run/cohort_merge/build.rs) and left there (§8).
- **Carrying the locus kind and motif on `CohortObservation`** — to
  [`cohort_merge.md`](cohort_merge.md), whose type it is (§2). This document cannot be built until
  it lands.
- **A purity adjustment on a candidate's stutter level** — not this step's at all; recorded at
  [`read_likelihoods.md`](read_likelihoods.md) §4.6 as deferred there, and repeated here only so
  that nobody reads §3's rung keying as the place it should live.

---

## 12. Resolved decisions and open questions

**Resolved.**

1. **The ladder is ported**, keyed by repeat count, because the genotype prior is written on the
   same index (§3).
2. **The peak test is not ported.** Measured: production offers both alleles of an
   adjacent-length heterozygote at 33–78% depending on coverage, against 97–100% for the support
   rule, with fewer candidates (§4.1).
3. **Both spellings of one length face the same support rule**, with no privileged representative
   and no recurrence bar. Measured: 35% against 86% at 300× on the class that is 43% of HG002's
   heterozygous tracts (§5).
4. **The share is 5%, one number for both paths** (§5).
5. **The depth gate, the three sibling constants and the cap of 24 are dropped** (§6).
6. **`NotPeriodic` survives, per sample** (§7).

**Open.**

- **Q1 — the periodicity threshold of 10% is inherited from production and has never been
  measured**, here or there. It decides how many catalog tracts are refused as not really
  periodic, and §4.1's measurements were all taken with the gate switched off, so nothing here
  constrains it. **Leaning: ship 10% and mark it soft.** *What would settle it:* on HG002, sweep
  it and count tracts refused against true heterozygous tracts lost, using the same rig as §5.
- **Q2 — how many alleles may a repeat tract be called over? SETTLED at 32** (owner's decision,
  2026-08-24), against the sibling document's six. A tract carries more real alleles than a SNP
  does and the ladder widens with the cohort, so this is where
  [`candidate_alleles.md`](candidate_alleles.md) §4.2's extrapolation is most likely to come true;
  six was inherited from a SNP/indel setting with nothing behind it for tracts.

  **32 costs 528 diploid genotypes a sample a locus** — a locus over `A` alleles has `A(A+1)/2`,
  so six is 21 and 32 is 528, and the calling loop scores every sample against every genotype.
  That is the price, and it is the reason the number is 32 rather than HipSTR's effective ceiling:
  64 alleles would be 2,080.

  **What decided it is that HipSTR has no such limit at all.** Its `gen_candidate_seqs`
  ([`HaplotypeGenerator.cpp:175-235`](../../../../HipSTR/src/SeqAlignment/HaplotypeGenerator.cpp))
  keeps every tract sequence passing an admission test — one sample showing it with ≥2 reads *and*
  ≥20% of that sample's reads, or a cohort-level 5% test — and **never ranks or truncates**. The
  1,000 of `MAX_TOTAL_HAPLOTYPES` is a later safety check on the product across flank and repeat
  blocks, and when it trips HipSTR abandons the locus. So HipSTR controls the count purely by
  strict admission; 32 with truncation plus §4.1's missing-genotype rule keeps more loci callable
  than either that or production's refusal at 24.

  **Still not measured, and the measurement is still owed:** the only STR benchmark with a truth
  set is one sample. *What would confirm 32:* the tomato panel's tracts routed through the merge
  once §2's field exists, histogramming candidates per tract at 1, 4, 16 and 63 accessions. If
  that histogram's tail reaches 32 at 63 accessions, the number is too low for a large panel.

  **What the other STR callers allow, since it is the one piece of evidence available before that
  run** (read 2026-08-24 from the vendored trees; the owner asked for it against this question):

  | caller | budget at one tract | above it |
  |---|---|---|
  | **ng, as specified** | **6 alleles counting the reference — 5 alternatives** | truncate to the best six (§4.1) |
  | production's repeat-tract path | 24 candidates ([`candidate_set.rs:272-276`](../../../../src/ssr/cohort/candidate_set.rs)) | **refuse the locus**, every sample no-called |
  | **HipSTR** | **1,000 haplotypes** (`MAX_TOTAL_HAPLOTYPES`, [`genotyper_bam_processor.h:110`](../../../../HipSTR/src/genotyper_bam_processor.h)) | **abandon the locus** — "Aborting genotyping of the locus as too many candidate haplotypes were found" ([`seq_stutter_genotyper.cpp:609`](../../../../HipSTR/src/seq_stutter_genotyper.cpp)) |

  **HipSTR's thousand is not a thousand tract lengths, and the conversion is the useful part.** Its
  haplotype is a combination across blocks — left flank × the repeat block × right flank — and the
  cap is on the product (`ncombs_` is `nopts_` multiplied over the blocks,
  [`Haplotype.cpp:124-140`](../../../../HipSTR/src/SeqAlignment/Haplotype.cpp)). Each flank
  contributes at most `MAX_FLANK_HAPLOTYPES = 4` and contributes 1 when the flank is invariant, so
  the repeat block's own budget runs from **1,000 tract sequences** where both flanks are fixed
  down to **62** where both are maximally variable. **So ng's five alternatives sit somewhere
  between 12 and 200 times tighter than HipSTR's**, and roughly four times tighter than
  production's 24.

  **That comparison is weaker than it looks in one direction and stronger in another.** Weaker:
  HipSTR's cap and production's are both set where they are *never meant to bind* — both refuse the
  locus when they do, which is the behaviour §4.1 rejects, so their numbers are safety valves and
  not statements about how many alleles a tract carries. Stronger: ng's six is a **working part**
  at a tract in a way it is not at a SNP, and it is the only one of the three that stays useful
  when it binds. **The number still has to be chosen against a measurement rather than against
  these three**, and the measurement above is what would do it. *The cap's value is a field of
  `SsrSelectionConfig::shared`, so a different default for this path costs no code — only a second
  constant and a reason.*
- **Q3 — everything in §4.1 and §5 is one human sample.** The cohort axis is untested on this path
  entirely, because `benchmarks/ssr_hg002/` is single-sample by construction and the tomato panel
  has no repeat-tract truth set. **Leaning: ship as specified** — every rule here is per-sample by
  design, so the cohort enters only through the union, which cannot lose an allele. *What would
  settle it:* a multi-sample repeat-tract truth set, which this project does not have.

---

## 13. How we know it works

- **Unit, the ladder:** two sequences of the same length land on one rung and stay distinct inside
  it; a sequence whose length is not a whole number of motif units lands on the floored rung and
  is counted there once, not twice.
- **Unit, nomination:** a sample with 150 reads at length 10 and 150 at length 11 nominates both —
  the case production's peak test cannot nominate at all, and the single most important test in
  this document.
- **Unit, the rescue:** a sample resolving one length nominates its occupied `±1` neighbours and
  not its unoccupied ones.
- **Unit, spellings:** a rung carrying a pure sequence at 200 reads and an interrupted one at 20,
  in a run of one sample, yields both — production yields only the first, and the test asserts the
  difference rather than the value.
- **Unit, no cohort term:** adding a sample that shows only the reference tract changes neither
  the candidate set nor any other sample's nomination. This is the property that fails first if a
  cohort denominator creeps back in.
- **Unit, one sample:** every rule above produces the same answer with a cohort of one that it
  produces for that sample inside a cohort of sixty-three.
- **Unit, `NotPeriodic`:** a locus whose reads are all a non-unit offset from the mode yields the
  reference tract alone with that verdict, and one sample being periodic is enough to avoid it.
- **Differential, the oracle:** §10's two runs — reproduce production exactly with the old rules
  switched in, then report the move.
- **Regression, the definition of done:** HG002 at 30× through the whole caller with real
  candidates, scoring at least the 98% adjacent-length recall of §4.1 and the 86% same-length
  recall of §5. The harness for both already exists and produced every number in this document on
  2026-08-24: production's `ssr-pileup` into `examples/ssr_slip_dump`, scored against
  `benchmarks/ssr_hg002/truth/`.
