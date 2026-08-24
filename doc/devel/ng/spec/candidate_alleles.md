# ng step 6 — choosing the alleles a locus is called over

*Design spec, 2026-08-24. **No code yet — this settles the design.** First of two documents on
candidate selection. This one carries the rules both paths share **and** the ordinary SNP/indel
path; the repeat-tract path is
[`candidate_alleles_ssr.md`](candidate_alleles_ssr.md), which inherits everything here and
replaces only what a repeat tract needs replaced. The types and signatures are
[`../arch/candidate_alleles.md`](../arch/candidate_alleles.md).*

*Reads on: [`cohort_merge.md`](cohort_merge.md) §4.2 — the allele table this narrows. Read by:
[`read_likelihoods.md`](read_likelihoods.md) §3.3 (the pooled leftover this creates),
[`calling_priors.md`](calling_priors.md) §4 (which divides its concentration by how many
alternatives survive here), and [`calling_em_loop.md`](calling_em_loop.md) §4 (which named this
gap and whose discovery rounds re-enter through this document's admission rule).*

*Production's nearest equivalents are
[`per_group_merger.rs`](../../../../src/var_calling/per_group_merger.rs) (SNP/indel) and
[`ssr/cohort/candidate_set.rs`](../../../../src/ssr/cohort/candidate_set.rs) (repeat tracts).
Everything said about those files is a record of what they do, not a proposal to change them —
`src/ssr/` and `src/var_calling/` are frozen production.*

---

## 1. What this is

**When the caller reaches a locus it has to decide the short list of sequences it will consider
as possibly present.** The cohort merge has already collected every sequence any sample's reads
showed over that stretch of genome and unified them into one table
([`cohort_merge.md`](cohort_merge.md) §4.2); it narrows nothing. Everything downstream is defined
over the narrowed list: the read likelihood scores each observation against each candidate, the
genotype prior lays its mass over the genotypes those candidates make, and the VCF's `ALT` column
is what survived. **This document says how the narrowing is done.**

Until it exists the caller cannot run end to end on real data: the calling loop is built and
proven, but only over candidate sets handed to it from a test fixture
([`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md), which records this as its one
blocker, twice).

### 1.1 Goals

1. **Narrow the merge's table to the sequences worth the arithmetic**, on evidence that is a
   property of the data and not of the run's composition.
2. **Degrade across the committed range** — one sample to several thousand, a few reads a
   position to several hundred (`CLAUDE.md`; the same goal every sibling spec carries). Nothing
   in the admission rule may branch on either axis.
3. **Produce the pooled leftover the SNP/indel read likelihood needs**
   ([`read_likelihoods.md`](read_likelihoods.md) §3.3's `q_sum_other`), which has no other
   producer.
4. **Be a pure function of one locus and the run's frozen parameters**, so where it runs is a
   cost decision rather than a correctness one.

### 1.2 Non-goals, and what this does not do

- **It does not invent a sequence.** Only what some sample's reads actually reached can become a
  candidate. A reader coming from GATK or HipSTR will expect otherwise: GATK's candidates are
  events read off assembled haplotypes
  ([`EventMap.java:240`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/haplotype/EventMap.java)),
  and HipSTR admits alleles named in a reference panel with no read support at all
  ([`HaplotypeGenerator.cpp:190-204`](../../../../HipSTR/src/SeqAlignment/HaplotypeGenerator.cpp)).
  ng does neither, on either path.
- **It does not discover alleles hidden under stutter.** That is the calling loop's, it runs
  between whole runs of the frequency loop, and it ships off
  ([`calling_em_loop.md`](calling_em_loop.md) §4). What it *does* owe this document is stated in
  §6.3: a discovered sequence enters through this document's rule, not through one of its own.
- **It does not prune after genotyping.** Dropping alleles no sample's best genotype used is the
  loop's, and it renumbers every allele id, so it needs a remapping this step never mints
  ([`calling/mod.rs:170-186`](../../../../src/ng/calling/mod.rs)).
- **It does not decide what the VCF says**, including whether a locus with no surviving
  alternative is emitted at all. That belongs to emission (steps 10–13).
- **It does not re-ask a question the merge already asked.** The merge decides whether a locus is
  worth building at all ([`cohort_merge.md`](cohort_merge.md) §4.3); this step decides which of
  its alleles are worth calling, and never revisits the first question — see §6.2.

### 1.3 Vocabulary

Three terms, all used below before anything is built on them.

- **compared reads** — a sample's reads whose whole sequence over the locus was read off and
  compared with the reference. Reads that stopped inside the locus, and reads that covered it and
  showed nothing, are in neither the numerator nor the denominator of anything here
  ([`SampleLocusObservations::reads_compared_with_reference`](../../../../src/ng/locus_generation/mod.rs)).
  At a cohort locus this is the sum of a sample's rows in the merge's table, because the merge
  admits only complete observations onto alleles.
- **the leftover** — for one sample, the reads whose sequence was in the merge's table and did
  not survive selection, with their summed log error probability. The read likelihood consumes
  the mass; §5 explains why the count has to travel with it.
- **the verdict** — what selection says happened at this locus, beyond the list itself: whether
  the list is everything that cleared the bar or was cut to fit the cap.

---

## 2. Where it sits, and why the placement is free

**Input:** one assembled `CohortObservation` — the locus's span, its unified allele table with
the reference at index 0, and per sample one row per `(allele, read group)` carrying a read count
and five quality sums
([`build.rs:922`](../../../../src/ng/run/cohort_merge/build.rs),
[`:965`](../../../../src/ng/run/cohort_merge/build.rs),
[`:1089`](../../../../src/ng/run/cohort_merge/build.rs)).

**Output:** `CandidateAlleles` — the surviving sequences with the reference at index 0, built and
merged ([`calling/mod.rs:86`](../../../../src/ng/calling/mod.rs)) — plus the verdict of §6 and,
per sample, the leftover of §5.

**Decision: it runs inside the merge's builder, on the locus the builder has just assembled.**
The calling loop already runs there ([`../arch/calling_em_loop.md`](../arch/calling_em_loop.md)
§7), so selection and the loop are the same function on the same thread, and no narrowed locus
has to travel anywhere before it is called.

**It commutes, and here is why rather than an assertion that it does.** Selection reads one
`CohortObservation` and the run's frozen parameters. It reads no other locus, no accumulated
state, no clock and no random source, so its result is identical wherever it is called. The
placement is therefore a cost decision, and there are two costs pulling the same way: the buffer
between the builders and the ordered output holds narrowed loci rather than whole tables, and the
leftover of §5 is a sum over rows the builder is holding at that instant — run selection later
and every dropped allele's per-sample rows have to travel to it intact, to be added up and thrown
away.

---

## 3. The admission rule

**An alternative survives if some single sample's reads lent it enough — where "enough" is a read
count at low depth and a share of that sample's own reads at high depth.**

```text
an alternative allele a survives  ⟺  ∃ sample s :
        reads_s(a)  ≥  max( floor,  ceil( share × compared_reads_s ) )

floor = 2 reads          share = 5 in 100 of that sample's compared reads
```

**This is the merge's own keep rule, asked one level down.** The merge asks each sample whether
its *non-reference reads pooled* reach `max(floor, ceil(share × compared reads))` and builds the
locus if any one does ([`cohort_merge.md`](cohort_merge.md) §4.3); the type that answers it is
[`MinAltReads`](../../../../src/ng/run/cohort_merge/mod.rs), `required_of`
([`:451`](../../../../src/ng/run/cohort_merge/mod.rs)) and `reached_by`
([`:466`](../../../../src/ng/run/cohort_merge/mod.rs)). Selection asks the identical question of
each *allele* separately. **The type is reused rather than copied**, so the two rules cannot drift
apart and a sweep of one is a sweep of both.

**The two halves bind at opposite ends of the depth range, which is the whole reason the rule has
two halves.** At 3 compared reads the bar is 2 reads; at 300 it is 15. §3.3 measures what the
share half buys at depth.

### 3.1 Why a bar at all, when production has none

**Production's ordinary path admits every distinct byte sequence any sample showed**, with no
count, fraction, quality or strand test — including buckets with zero observations
([`per_group_merger.rs:1105-1128`](../../../../src/var_calling/per_group_merger.rs)). Its
`--min-alt-obs-per-sample`, default 2
([`variant_caller.rs:380`](../../../../src/var_calling/variant_caller.rs),
[`var_calling/mod.rs:72`](../../../../src/var_calling/mod.rs)), is a *record* survival test: if one
alternative has two reads, a second alternative with a single read rides into the genotyping beside
it. bcftools has no SNP allele bar either — any base seen at base quality 1 or more becomes an
`ALT` ([`bam2bcf.c:996-1002`](../../../../bcftools/bam2bcf.c)) — and GATK has none at generation
time, relying on ranked caps afterwards.

**Two costs make a bar worth having here, and the first is the larger.**

- **The genotype prior divides by the alternative count.** `α_alt(a) = θ / (number of
  alternatives)` ([`calling_priors.md`](calling_priors.md) §4, production's own
  [`genetics.rs:214`](../../../../src/genetics.rs)). Admit five alternatives of which four are
  sequencing error and the real one starts with a fifth of the concentration it should have.
  Production has this defect; ng would inherit it.
- **Genotypes grow faster than alleles.** A diploid locus has `C(A+1, 2)` genotypes — 3 at two
  alleles, 21 at six, 55 at ten — and the loop scores every sample against every genotype
  ([`calling_em_loop.md`](calling_em_loop.md) §8).

**freebayes and HipSTR both have a bar of this exact shape**, which is where the shape comes from:
freebayes wants at least 2 reads **and** at least 5% of one sample's observations, taking the first
sample that passes
([`AlleleParser.cpp:3913-3937`](../../../../freebayes/src/AlleleParser.cpp), defaults `-C 2`,
`-F 0.05`); HipSTR wants at least 2 reads **and** at least 20% of one sample's spanning reads
([`HaplotypeGenerator.cpp:181`](../../../../HipSTR/src/SeqAlignment/HaplotypeGenerator.cpp)).

### 3.2 One sample suffices, and no term reads the cohort

**This is the rule the rest of the design hangs on: no term of the admission bar may read the
cohort.** The bar is asked of each sample separately against that sample's own reads, and one
sample reaching it admits the allele for everyone. The cohort enters selection in exactly two
places — the union of what the samples nominated, and the tie-break the cap uses when it binds
(§4) — and neither can remove an allele a sample's own reads earned.

**Why, in one sentence: otherwise a sample's candidate list depends on who else is in the run.**
That is the same ground on which the repeat-tract junk term was rewritten to stop reading the
cohort ([`read_likelihoods.md`](read_likelihoods.md) §4.5, the owner's decision of 2026-08-19),
and the same ground on which the merge retired its cohort-summed keep rule
([`cohort_merge.md`](cohort_merge.md) §4.3 and §7.3, on a measurement: the old rule discarded 997
loci in 1,000 at one tomato accession and only 439 at 63).

**Two alternatives were live and both lost.**

| | why it lost |
|---|---|
| a bar on the **cohort's total** reads for the allele (production's ranking key, freebayes' `--min-alternate-total`) | it grows with the cohort while a rare allele's evidence does not. At 63 tomato accessions a 2% cohort share asks for about 14 reads, which one or two carriers at 3 reads a position cannot supply; at 1,000 samples it asks for about 200 |
| a bar requiring **two or more samples** | it cannot admit a private allele, and finding those is what a large panel is for. Recurrence is still used — as the cap's tie-break (§4), where it can only reorder, never exclude |

### 3.3 The numbers, and what they cost

**Both are inherited from the merge and both were, until now, unmeasured against a truth set** —
[`cohort_merge.md`](cohort_merge.md) §4.3 books that debt in its own last paragraph. This section
pays it for the allele-level rule.

*Measured 2026-08-24 on the GIAB trio (HG002/3/4) over the 100 benchmark intervals of
`benchmarks/giab/per_sample/bed/HG002_bench_azar_merged_100.bed`, 572 kb, scored against the
v4.2.1 truth VCFs. "Lost to the bar" means some sample's reads did show the allele and the bar
rejected it; of the 920 true alternative alleles carried by at least one sample inside those
regions, 6 fall where the merge built no locus and 1 no sample showed, and those are the merge's
and the reads', not this step's.*

At **30×** — 5,777 alternatives in the merge's tables:

| bar | alternatives kept | true alleles lost to the bar |
|---|---|---|
| 2 reads, no share | 3,091 | 1 |
| 2 reads or 2% | 3,091 | 1 |
| **2 reads or 5%** | **2,977** | **1** |
| 2 reads or 10% | 1,601 | 2 |
| 3 reads or 2% | 1,539 | **5** |

At **300×** — 15,474 alternatives:

| bar | alternatives kept | true alleles lost to the bar |
|---|---|---|
| 2 reads, no share | 10,793 | 2 |
| 2 reads or 2% | 5,596 | 2 |
| **2 reads or 5%** | **2,308** | **2** |
| 2 reads or 10% | 1,273 | 4 |

**Two things follow.** The share is free up to 5 in 100 at both depths while removing three
quarters of the table at 300×, which is why the share here is 5% where the merge's is 2% — and it
changes nothing below about 40 compared reads a sample, so a tomato-depth run sees the identical
rule it would have seen at 2%. And **the floor is the expensive knob**: raising it from 2 to 3
loses five true alleles where raising the share to 10% loses two for the same reduction in table
size. **The floor stays at 2 and should be defended there.**

**The share is what controls the allele count at depth, not the cap.** On the same 300× trio, a
count-only bar leaves 471 of 7,478 loci carrying more than three alternatives — 6.3% — while a
2-in-100 share brings that to 16 loci, 0.2%, and 5 in 100 to a single locus. The cap of §4 is
therefore guarding against something the share has already almost eliminated at that depth.

**At tomato depth the share is inert and the bar is a count**, which is worth stating plainly
because it is the thin end of the committed range. On 63 tomato accessions at about 3 reads a
position over 400 kb (53,935 built loci, same date), a share of 0% and a share of 2% differ by 4
loci in 53,935. §7 says what that costs.

**A third bar was considered and is not built.** freebayes carries `--min-alternate-qsum`,
defaulted to 0 ([`Parameters.cpp:484`](../../../../freebayes/src/Parameters.cpp)), and the merge
already stores `q_sum` per `(sample, allele, read group)`, so a bar on summed read quality would
cost nothing to compute. It is left out because it is a second unmeasured number and the read
count already carries most of the signal; **it is the first thing to add if the run of §12's Q1
shows the count bar admitting error alleles at high depth.** Recorded here so nobody invents a
producer for it.

### 3.4 The numerator is supplied, not fixed

**The rule above is a predicate, and its numerator is an argument.** Selection supplies "reads of
this sample that showed this allele". A discovery round supplies a different, narrower numerator —
reads of this sample that the converged model is currently explaining as slippage and which imply
this sequence — against the same denominator, that sample's reads at the locus, and the same two
constants.

**This retires a second bar that already exists in the architecture.**
[`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §6.2 types a
`DiscoveryBar { min_reads: 2, min_spanning_read_share: 0.15 }`, HipSTR's high-depth human
setting, marked inherited and soft. Its own spec then argues that a single pair of numbers cannot
be right at 3 reads a position and at 300 ([`calling_em_loop.md`](calling_em_loop.md) §4.1) —
which is the argument that produced the two-part rule above. **Two different pairs of numbers
admitting alleles into one table at one locus is the thing to avoid**, so `DiscoveryBar` should
go and discovery should call this rule. That is an edit to a sibling document and is not made
here.

---

## 4. The cap, and what happens above it

**A locus is called over at most six alleles including the reference. Above that the list is cut
to the best six; the locus is never refused.**

**Six is inherited.** It is production's `DEFAULT_MAX_ALLELES_PER_RECORD`
([`per_group_merger.rs:57`](../../../../src/var_calling/per_group_merger.rs)) and GATK's
`--max-alternate-alleles` default
([`GenotypeCalculationArgumentCollection.java:29`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/genotyper/GenotypeCalculationArgumentCollection.java)),
and ng's own documents already do their arithmetic at that width — "21 at six alleles"
([`calling_em_loop.md`](calling_em_loop.md) §1.3). **It has never been measured and is soft.**
§4.2 says how often it binds, which is the fact that decides whether its value matters.

### 4.1 Truncate, never refuse — and the ranking

**Truncation, because refusal throws away the good alleles with the bad.** A locus carrying two
obvious variants and six noise sequences is refused whole under production's repeat-tract rule
(over 24 candidates, every sample no-called,
[`candidate_set.rs:272-276`](../../../../src/ssr/cohort/candidate_set.rs)) and under HipSTR's
(over 1,000 haplotypes, the locus is abandoned). GATK, freebayes, bcftools' indel path and
production's own ordinary path all truncate instead. **Truncation is only defensible because §5's
leftover exists** — the cut alleles' reads keep their error mass in the arithmetic instead of
being silently reassigned to the reference, which is what GATK does
(`AlleleLikelihoods.marginalize` takes a maximum over the collapsed alleles and drops the rest;
`AlleleSubsettingUtils.subsetAlleles` deletes the genotype entries and subtracts the maximum).

**The ranking: the largest share of one sample's compared reads the allele took, maximised over
samples.** Ties break on how many samples cleared the bar, then on the cohort's read total, then
on the bases themselves so the result is deterministic at any worker count.

**This is not production's ranking, and the difference is what it does at scale.** Production
ranks by the cohort's raw read total
([`enforce_max_alleles`, `per_group_merger.rs:1434`](../../../../src/var_calling/per_group_merger.rs)),
so at a thousand samples the alleles pushed into the leftover are the private ones: a real allele
carried by one sample at 30× scores 15 reads, while a systematic mismapping artefact at 1% in 800
samples scores 240. A within-sample share inverts that — 0.5 against 0.01. **Two other callers
normalise per sample before pooling for the same reason**: bcftools divides each sample's quality
sums by that sample's total before adding across the cohort
([`bam2bcf.c:968-975`](../../../../bcftools/bam2bcf.c)), and HipSTR pools `count/sample_reads`
rather than counts
([`HaplotypeGenerator.cpp:183`](../../../../HipSTR/src/SeqAlignment/HaplotypeGenerator.cpp)).

**The ranking degrades correctly across the depth range without a branch**, which is the reason
for the tie-break order. At 300 reads a sample the shares separate cleanly and the first key
decides. At 3 reads every admitted allele has a share near 0.67, the first key ties, and how many
samples showed it decides — which is the only signal there is at that depth.

**No second, harder ceiling is needed.** Production carries one at 64 alleles, where the group is
refused outright, because its genotype enumeration uses a `u64` bitmask
([`per_group_merger.rs:936-940`](../../../../src/var_calling/per_group_merger.rs)). ng has no
bitmask, and with the bar of §3 applied before the cap no locus reaches the loop above six. What
remains is `admit`'s own refusal at 65,536 alleles, which exists so a discovered allele cannot
wrap onto the reference's index
([`calling/mod.rs:170-186`](../../../../src/ng/calling/mod.rs)) and is a backstop, not a policy.

### 4.2 How often it binds — measured

*Same runs as §3.3, with the bar at 2 reads or 2%.*

| cohort | built loci | loci above six alleles | loci above four |
|---|---|---|---|
| 63 tomato accessions, ~3 reads a position, 400 kb | 53,935 | 23 (1 in 2,300) | 65 (0.12%) |
| GIAB trio, 30×, 572 kb | 4,177 | 0 | 0 |
| GIAB trio, 300×, 572 kb | 7,478 | 0 | 4 (0.05%) |

**So at these cohort sizes the cap is a safety valve, not a working part**, and its exact value
carries little. **But what it guards against grows with the cohort**, which is why it is here.
Holding the allele table fixed at the 63-accession one and varying only how many samples the bar
is asked of:

| samples asked | alternatives admitted per locus | loci with none | loci above six | most at one locus |
|---|---|---|---|---|
| 1 | 0.02 | 98.1% | 0 | 2 |
| 4 | 0.11 | 88.9% | 0 | 4 |
| 16 | 0.22 | 78.9% | 3 | 14 |
| 63 | 0.77 | 27.4% | 23 | 14 |

Part of that growth is real segregating variation, not error — [`cohort_merge.md`](cohort_merge.md)
§7.3 makes the same point about the same panel — but the prior of §3.1 dilutes either way, so the
cap has to be there. **Whether it becomes a working part at a thousand samples is an
extrapolation from this table, not a measurement** (§12, Q2).

**And the ranking's advantage is, honestly, unmeasurable at these sizes.** Ranking by
within-sample share and production's ranking by cohort read total keep different alleles at 17 of
53,935 tomato loci and at none of the trio's. The argument for the share ranking is a
thousand-sample argument and no thousand-sample cohort exists here (§12, Q2).

---

## 5. The leftover, and why it carries a count

**Every allele the bar or the cap removes keeps its reads' error mass in the arithmetic.** The
SNP/indel genotype likelihood carries a term for reads matching no candidate —
`q_sum_other`, the pooled error mass, which is the same under every genotype, cancels in
genotyping, and is kept because the data likelihood feeds emission and `QUAL`
([`read_likelihoods.md`](read_likelihoods.md) §3.3, and the field
`GenericSampleEvidence::unmatched_q_sum` in
[`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §2.1). **Nothing upstream produces
it, because nothing upstream drops anything. Selection is what creates the pool, so selection owes
it.**

**It costs nothing to produce.** The merge already stores `q_sum` per `(sample, allele, read
group)` ([`build.rs:1245`](../../../../src/ng/run/cohort_merge/build.rs)), so a sample's leftover
is the sum over the alleles selection dropped. There is no new upstream producer and no second
pass over anything.

**Decision: the leftover carries a per-sample read count beside the mass, and that is what makes
truncation defensible.** The mass alone is identical under every genotype, so a sample whose true
allele was cut is scored confidently against a set that does not contain it, with nothing
per-sample saying so. A count of that sample's reads in the pool is what lets a later step
no-call it. It costs one integer per sample per locus. **Truncation (§4.1) and this count are one
decision, not two: if the count is ever dropped, refusing the locus becomes the correct policy
again.**

**Production produces the same thing and shows the trap.** Its cap pools the dropped alleles'
per-sample support into `other_scalars`, whose `q_sum` enters every genotype's likelihood
([`pool_dropped_other_scalars`, `per_group_merger.rs:1555`](../../../../src/var_calling/per_group_merger.rs);
consumed at [`:1979-1986`](../../../../src/var_calling/per_group_merger.rs)). **And when
contamination is configured the loop recomputes likelihoods from the per-allele table and never
reads it** ([`posterior_engine.rs:1480`](../../../../src/var_calling/posterior_engine.rs) and the
substitution around [`:2429`](../../../../src/var_calling/posterior_engine.rs)), so the pool
silently vanishes on that path. ng's own likelihood has the same contamination mixture
([`read_likelihoods.md`](read_likelihoods.md) §3.6) and so has the same omission available.
Flagged as a trap, not as a change to production.

**Measured size** (2026-08-24, bar at 2 reads or 2%): the leftover takes 0.36% of the reads on the
63-accession tomato panel and 1.05% of the trio's at 300×.

### 5.1 What is *not* in the leftover

- **Partial reads.** A read that stopped inside the locus does not say what the sample carries, it
  says the sample carries *at least* this, and it is held on its own axis with its own bases and
  witnessed stretch ([`SampleSupport::partials`](../../../../src/ng/run/cohort_merge/build.rs),
  `PartialObservation` at [`:1130`](../../../../src/ng/run/cohort_merge/build.rs)). Selection does
  not read them, does not count them toward any bar, and does not fold them into the pool: they
  are not reads matching no candidate, they are reads matching a *set* of candidates, which
  [`read_likelihoods.md`](read_likelihoods.md) §5.3 scores. **Selection can empty that set** — a
  partial whose only compatible allele was dropped — and that is recorded here so §5.3's author
  knows the state is reachable; what to do about it is theirs.
- **Reads that produced no observation, and reads removed as evidence.** Both are counted by the
  merge and neither carries a quality sum
  ([`build.rs:1036-1052`](../../../../src/ng/run/cohort_merge/build.rs)), so neither can join a
  pool of error mass. They were never in the allele table and selection did not drop them.
- **The reference's own reads.** The reference allele is never dropped (§6.1).

---

## 6. The verdict, the reference, and what the loop may do afterwards

### 6.1 The reference is allele 0 and is never a candidate for removal

It is placed before any sample's evidence is read and is exempt from both the bar and the cap. The
type makes it impossible to build a table without it and impossible to move it
([`calling/mod.rs:86-118`](../../../../src/ng/calling/mod.rs)); production, freebayes and GATK all
force it to index 0 for the same reason — every downstream branch tests against it.

### 6.2 The verdict — two cases, and a deliberate absence

```text
Selected                    the list is everything that cleared the bar
Truncated { dropped: u16 }  the cap bound; `dropped` alternatives were cut
```

**There is no depth verdict, and its absence is a decision.** The architecture sketched
`Admission { Ok, LowDepth, NotPeriodic, TooManyAlleles }`
([`../arch/ng_step_interfaces.md`](../arch/ng_step_interfaces.md) §3) and it was never built —
today's `CandidateAlleles` has no verdict field at all. Two of those four do not survive what this
document settles. `TooManyAlleles` named a refusal and §4.1 chose truncation. **`LowDepth` would
re-ask the merge's keep rule with a different denominator**, and production's version of it is a
sum over the cohort — measured refusing 98.6% of repeat tracts at one sample sequenced to 5× and
0.2% at 300× ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §6). Depth is asked once,
upstream, per sample. `NotPeriodic` survives and belongs only to repeat tracts
([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §7).

**A locus can select down to the reference alone, and that is a first-class outcome.** The merge
builds a locus when some sample's non-reference reads *pooled* reach its rule; two reads split one
and one across two alternatives clear that and clear neither allele bar. Measured, this is more
than one built locus in four and the fraction is the same on both benchmarks: 27.4% on tomato,
27.3% on the trio at 30×, 28.0% at 300×. **It is `Selected` with an empty alternative list, it is
counted, and what the run does with it is emission's.** It is not an error and must not become
one.

### 6.3 What the calling loop may do afterwards

The loop owns three things this document does not: when a discovery round runs, where its
numerator comes from, and the prune with its id remapping
([`calling_em_loop.md`](calling_em_loop.md) §4). Two rules bind them to this document.

- **A discovered sequence enters through §3's rule**, with the numerator of §3.4. The
  architecture already asks for this in words — "discovered sequences enter `CandidateAlleles`
  through the same admission as selected ones"
  ([`../arch/calling_em_loop.md`](../arch/calling_em_loop.md) §6.2) — and until now there was no
  admission to enter through.
- **When the table is already at the cap, a discovery round's proposal is refused and counted, not
  admitted with a re-ranking.** Evicting an allele mid-run changes the repeat-tract likelihood's
  junk denominator, which is defined over the candidate set
  ([`read_likelihoods.md`](read_likelihoods.md) §4.5), and so changes every genotype's
  likelihood — the exact invalidation the "only between whole runs of the loop" rule exists to
  prevent ([`calling_em_loop.md`](calling_em_loop.md) §4's table). There is also no ranking
  argument for preferring an allele found under stutter to one that cleared the bar on a sample's
  own reads.

---

## 7. One sample and a thousand, three reads and three hundred

**No term of the rule branches on either axis** (§3.2), so this section is about what the one rule
gives at each end rather than about four rules.

**At three reads a position the bar is a count and nothing can repair that inside selection.**
There is no fraction of three reads between one third and two thirds, so every conjunctive bar in
every caller here is inert at that depth: freebayes' 5%, HipSTR's 20%, ours at 5%. The count is 2,
and a heterozygous carrier at 3 compared reads shows two or more alternative reads half the time —
so **half of that sample's heterozygous sites are not offered as candidates from that sample**.

**What recovers them is the union across samples, not the bar**, and that is the positive argument
for choosing the list once across the cohort. With 63 accessions and an allele carried by six of
them, the chance that no carrier clears the bar is under 2 in 100.

**At one sample there is no union, and the loss above is what 3× single-sample sequencing costs.**
It is a property of the data, not of this design, and the two repairs available are both worse: a
floor of 1 mints an allele from every sequencing error, and a floor that changes at `n = 1` is a
constant with no measurement behind it and makes a sample's calls depend on the run's composition.
**What protects a low-depth call is not the bar.** At three reads a position the read likelihoods
are nearly flat across genotypes ([`calling_em_loop.md`](calling_em_loop.md) §7), so an allele
admitted on two reads gets a posterior the prior can overrule and contributes almost nothing to
the cohort's expected copies — and at one sample the prior's cohort term is exactly zero, which is
why one sample at three reads is the weakest corner of the committed range. The false-allele rate
there is §12's Q1.

**At three hundred reads the share does the work** and the count is never the binding half: §3.3
measured a count-only bar leaving 6.3% of loci above three alternatives at 300× against 0.2% once
a 2-in-100 share is added.

**At a thousand samples the rule asks the same question it asks at one**, but more samples are
more independent chances for one of them to reach the bar, so what is *admitted* grows even where
the locus does not (§4.2's table). That is what the cap and its ranking are for, and it is the one
part of this design whose behaviour at that size is extrapolated rather than measured.

---

## 8. Cross-cutting concerns

**Cost.** One pass over the locus's rows, which the builder is already holding: for each sample,
one walk of its `supported` rows accumulating per allele. `alleles × samples` work, no allocation
beyond the surviving table itself, and it runs before the loop whose per-pass cost is
`samples × genotypes` — so selection is cheap in the place it matters most, by making the second
number smaller.

**Memory.** Selection *reduces* what the builder passes on: the narrowed table plus one leftover
pair per sample, against the full table plus every dropped allele's rows. That is the second
reason for the placement of §2.

**Errors.** A locus that selects to the reference alone is not an error (§6.2). A locus above the
cap is not an error (§4.1). The failure modes that *are* errors are caller bugs — a support row
naming an allele outside the table, a non-finite quality sum, a compared-read count of zero on a
sample that has rows — and they are assertions, structural ones held in release, which is the
convention [`../arch/read_likelihoods.md`](../arch/read_likelihoods.md) §1.1 already sets.

**Determinism.** The output must be byte-identical at any worker count, which the merge already
requires of itself ([`cohort_merge.md`](cohort_merge.md) §9). Selection is a pure function of
one locus (§2), so the only exposure is the cap's ranking, and its final tie-break on the allele's
bases closes it: no two distinct alleles can compare equal. **The ranking must not rest on a
stable sort of equal keys** — production's does, which is deterministic but arbitrary
([`per_group_merger.rs:1434`](../../../../src/var_calling/per_group_merger.rs)) — because the row
order it would inherit is the order samples were walked in.

---

## 9. Reuse map

| what | existing code | how ng reuses it |
|---|---|---|
| the two-part support rule and its constants | [`MinAltReads`, `cohort_merge/mod.rs:424`](../../../../src/ng/run/cohort_merge/mod.rs) | **the type itself, unchanged**, with a different numerator. The share moves from 2% to 5% for this rule only; the merge keeps its own (§3.3) |
| the allele table and its reference invariant | [`CandidateAlleles`, `calling/mod.rs:86`](../../../../src/ng/calling/mod.rs) | built and merged; selection fills it through `admit` |
| the evidence | [`CohortObservation`, `cohort_merge/build.rs:922`](../../../../src/ng/run/cohort_merge/build.rs) | read, not changed |
| the cap's value | [`DEFAULT_MAX_ALLELES_PER_RECORD`, `per_group_merger.rs:57`](../../../../src/var_calling/per_group_merger.rs) | the number is inherited and declared inherited; the ranking is **not** ported (§4.1) |
| pooling the dropped alleles' error mass | [`pool_dropped_other_scalars`, `per_group_merger.rs:1555`](../../../../src/var_calling/per_group_merger.rs) | the shape is the same; ng adds the read count (§5) and does not repeat the contamination omission |

**There is no parity oracle for this path, and inventing one would be worse than saying so.**
Production has no per-allele bar at all (§3.1), so there is nothing to be byte-identical with. The
repeat-tract path gets a differential instead
([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) §10).

---

## 10. Deferred, with a recommended home

- **A bar on summed read quality** — to §12's Q1. Recorded in §3.3 with what would trigger it, so
  nobody invents a producer for a number this step does not currently need.
- **What a locus with no surviving alternative becomes in the output** — to the emission step's
  spec, which owns what a calling run writes. This document fixes only that the state exists, is
  legal, and is counted (§6.2).
- **Retiring `DiscoveryBar`** — to [`../arch/calling_em_loop.md`](../arch/calling_em_loop.md)
  §6.2, whose type it is. §3.4 gives the reason; the edit is not made here.
- **The rule for a partial read whose compatible set selection emptied** — to
  [`read_likelihoods.md`](read_likelihoods.md) §5.3, which owns how a partial is scored. §5.1
  records that selection can create the state.
- **Whether the emitted depth counts the leftover's reads** — to emission. Production's does not,
  a documented convention that makes its `DP` not the position's depth
  ([`posterior_engine.rs:877`](../../../../src/var_calling/posterior_engine.rs)); ng should decide
  rather than inherit by accident.

---

## 11. Resolved decisions and open questions

**Resolved.**

1. **The bar** — `max(2 reads, 5% of that sample's compared reads)`, asked per sample, one sample
   suffices (§3). It beat a cohort-total bar (drifts with cohort size, measured) and a
   two-sample bar (cannot admit a private allele).
2. **The cap** — six alleles including the reference, truncate rather than refuse, ranked by the
   largest within-sample share (§4). Refusal loses the good alleles with the bad; production's
   cohort-total ranking truncates the private alleles first at scale.
3. **Cohort or sample** — one list per locus across every sample; the per-sample reasoning lives
   entirely in the bar (§3.2).
4. **The leftover** — per sample, summed from the merge's own `q_sum`, carrying a read count
   beside the mass (§5). The count and truncation stand or fall together.
5. **Where it runs** — the merge's builder; it commutes because selection is a pure function of
   one locus, so the choice is made on memory (§2).
6. **The verdict** — `Selected` / `Truncated { dropped }`, no depth verdict (§6.2).

**Open.**

- **Q1 — the false-allele rate at one sample and three reads a position.** The number that decides
  whether a floor of 2 is defensible at the weakest corner of the committed range, and it is
  computable before it is measured: the parameters pre-pass fits a per-read-group error rate, and
  the chance that two of a sample's three reads carry the same wrong base at the same position
  follows from it. **Leaning: the floor stays at 2.** *What would settle it:* the arithmetic from
  the fitted rate first, then HG002 downsampled to 3× as a single sample, counting alternatives
  admitted outside the truth set's confident regions. Confirm before lowering the floor.
- **Q2 — does the cap become a working part at a thousand samples, and does the ranking then
  matter?** §4.2's growth table says the admitted count rises with the cohort and extrapolates to
  the cap binding, but 63 accessions is as far as the measurement goes, and at that size the two
  rankings disagree at 17 loci in 53,935. **Leaning: keep both as specified** — the cap costs
  almost nothing when it does not bind, and the ranking's tie-break order is what makes it degrade
  correctly with depth regardless. *What would settle it:* a cohort of several hundred samples, or
  a resampling study that holds the allele table fixed and grows the sample set past 63.
- **Q3 — the share of 5% is measured on one human trio over 572 kb.** It is free there at both 30×
  and 300× and it is inert below about 40 compared reads a sample, so no low-coverage run is
  exposed to it. **Leaning: ship 5%.** *What would settle it:* the same scoring on a second
  high-depth cohort, or on the trio over a wider region set.

---

## 12. How we know it works

- **Unit, the bar:** at 3 compared reads the rule asks 2; at 300 it asks 15; a sample showing 2 of
  3 passes and a sample showing 2 of 300 fails. One sample passing admits the allele for a cohort
  in which no other sample shows it at all.
- **Unit, the cap:** a locus with eight alternatives clearing the bar yields five plus the
  reference, `Truncated { dropped: 3 }`, and the five kept are the five with the largest
  within-sample share. A locus with six or fewer yields `Selected`.
- **Unit, determinism:** two alternatives with identical shares, identical sample counts and
  identical read totals are ordered by their bases, and reversing the input row order does not
  change the output.
- **Unit, the leftover:** the dropped alleles' `q_sum` and read counts sum into the per-sample
  pool exactly; a sample with no dropped alleles has a pool of zero and no branch is needed to get
  it.
- **Unit, the reference:** a locus where every alternative fails the bar yields a table of length
  one, `Selected`, and is counted.
- **Property:** the surviving list is a subset of the merge's table, contains the reference, and
  its length never exceeds the cap — for any input, including tables built from random support.
- **Property, scale-freedom:** adding a sample that shows only reference reads changes neither the
  surviving list nor any other sample's leftover. This is §3.2's principle as a test, and it is
  the one that fails first if a cohort term creeps into the bar.
- **Regression, the definition of done:** the GIAB trio and the tomato panel run end to end
  through the calling loop with real candidates rather than fixtures — which is the blocker
  [`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md) records, and clearing it is what
  this document is for.
- **The measurement harness exists**: `examples/ng_candidate_selection_probe.rs` walks a real
  cohort, applies the bar and the cap, and with `NG_SELECT_DUMP` writes a per-allele table for
  truth scoring. Every number in §3.3, §4.2 and §5 came from it on 2026-08-24.
