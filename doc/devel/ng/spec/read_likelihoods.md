# ng — the read likelihood models

*Design spec draft, 2026-08-19. **No code yet — this settles the design.** Second of three documents
on variant calling; the first is the **genotype prior**
([`calling_priors.md`](calling_priors.md)) and the third is the **EM loop** that ties them together,
which does not exist yet. This document defines what the prior gets multiplied by: the caller adds
the two in log space and normalises, and the result is the posterior ng emits.*

*Reads on: [`cohort_merge.md`](cohort_merge.md) — what a cohort observation is, which is this
document's input; [`parameter_prepass.md`](parameter_prepass.md) and its two path siblings
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) and
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) — every parameter this model consumes and
none of which it fits; [`alignment.md`](alignment.md) §5.1 and §5.2 — the sequence comparison this
model composes, and the stutter distribution the two documents share;
[`locus_generation_ssr.md`](locus_generation_ssr.md) §3 — what an STR allele is and what a read that
saw only part of a tract claims. Production's equivalent code is
[`per_group_merger.rs`](../../../../src/var_calling/per_group_merger.rs) and
[`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs) (SNP/indel) and
[`ssr/cohort/read_model/`](../../../../src/ssr/cohort/read_model/) (STR). Everything said about
those files is a record of what they do, not a proposal to change them — `src/ssr/` and
`src/var_calling/` are frozen production.*

---

## 1. What this is

**Given that this sample's two copies of this locus are such-and-such, how probable are the reads we
actually saw here?** That is the whole question. The answer is one number per candidate genotype per
sample, in log space, and it is written `Lg` — the **genotype likelihood**. It is built out of a
smaller number, `Lr` — the **read likelihood**, the probability that one copy of one allele produced
one observed sequence.

The caller adds `Lg` to the log-prior of [`calling_priors.md`](calling_priors.md) and normalises
across genotypes. The prior says what is likely before any read is looked at; this says what the
reads change.

*(Notation: `Lr` and `Lg`, never `Q`. `Q` means Phred throughout this project and a second meaning
would collide with it. Where production's code calls the read likelihood `Qᵣ`, this document renames
it on quotation.)*

There are two of these models, because there are two kinds of locus and reads go wrong at them in
different ways:

- **At an ordinary site**, a read can show a base the individual does not carry. One rate covers it.
- **At a repeat tract**, a read can also show a whole repeat more or fewer than the individual
  carries, because the copying enzyme slipped. That is **stutter**, it happens thousands of times
  more often than a base misread, and a model with one substitution rate has nowhere to put it.

**This model fits nothing: it is a pure function of the numbers it is handed and the reads it is
shown.** Every number in it — the per-read-group error rate; the four numbers the parameters fit measures per
stratum, which are how often a read slips, which way, how far, and the substitution rate inside a
tract; and the contamination fractions — is fitted before calling starts
([`parameter_prepass.md`](parameter_prepass.md) §1). §7 lists each one and where it comes from, and
what this document decides is the *shape* the numbers are used in.

**"Fits nothing" is not the same as "everything stays frozen", and an earlier version of this
sentence blurred the two.** Whether the caller *re-estimates* a parameter between its own iterations
and hands this model a new value is the EM loop document's decision, not this one's, and it is a real
decision rather than a formality: **production's STR caller re-fits two of the seven stutter
parameters per locus**, in up to three rounds, shrunk toward the pre-pass's per-stratum value
([`em.rs:69`](../../../../src/ssr/cohort/em.rs)). §6.1 sorts every parameter into what must stay
frozen, what is genuinely a candidate for re-estimation, and what the caller re-estimates without
this model ever seeing it — which is what the EM document will need.

### 1.1 Goals

1. **One formula for both paths, differing in one term.** Both are the same formula (§2), and
   what changes is the **emission** — the answer to *how does one copy of one allele produce one
   observed sequence*, which is a base-quality term at an ordinary site and a slippage term at a
   repeat tract.
2. **Exactly computable from what the merge hands over.** A cohort observation is a distinct observed
   sequence with a count, the ids of the reads that showed it, and a few *summed* per-read
   quantities (§1.4). The ids let each read's haplotype be built; what is summed away is each read's
   quality. **A likelihood that treats an observation's reads as interchangeable when their
   qualities differ is wrong by an amount nobody measured.** §2.3 makes that a contract and §12
   makes it a test.
3. **A read carries a base quality, and the pre-pass fits an error rate. They are different
   things, and the model uses both.** A base quality is the sequencer's claim about one read, and it
   is the only number that says this read is better than that one. The fitted rate covers everything
   else that makes a read show an allele the sample does not carry — a mismapped read, a chimera,
   contamination — and it was measured on this data rather than claimed by the machine. §3.2 keeps
   both: the base qualities say which reads differ, the fitted rate says how large the errors are.
4. **Degrade across the committed range** — one sample to several thousand, a few reads a position
   to several hundred (`CLAUDE.md`, *What this caller has to work on*). §6 answers both ends of both
   axes, and the answer at one sample is *unchanged*, not *degraded*: the read likelihood is a
   per-sample quantity and no part of it should need a cohort. Production's STR model breaks that in
   one place — it spreads its outlier weight over the distinct sequences the whole cohort showed —
   and §4.5 is about it.
5. **Ploidy-generic and multi-allelic from the first line of code**, for the same reason the prior
   is: the biallelic-diploid shortcut is what production had to unpick, and the general form costs a
   loop bound.
6. **A pure function.** Same evidence and same parameters give the same numbers at any thread count,
   in any order.

### 1.2 Non-goals, and what this document does not do

- **It does not fit anything.** See above. If a number in this document has no named source in the
  pre-pass, that is a defect in this document.
- **It does not choose the candidate alleles.** That is candidate generation
  ([`ng_proposal.md`](ng_proposal.md) step 6). This model scores the alleles it is handed.
- **It does not run the EM, and it does not compute a posterior or a genotype.** Those are the third
  document's. This one defines a function; that one decides when it is called and with what.
- **It does not decide whether a site is emitted or what its QUAL is.** Emission is a separate
  question from genotyping ([`ng_proposal.md`](ng_proposal.md) step 11). What this document owes
  emission is a **data likelihood** that is comparable between loci — which is why §3.3 keeps a term
  that cancels in genotyping and §4.5's open question matters more to emission than to calls.
- **It does not align reads.** The comparison of two sequences under a per-base error rate belongs
  to the alignment module ([`alignment.md`](alignment.md) §5.1) and this model composes it. §7 draws
  the line, because [`alignment.md`](alignment.md) §5.2 already draws it from the other side and the
  two documents must not disagree.
- **It does not use phasing, linkage, or the reads at any neighbouring locus.** Each locus is scored
  on its own reads. `chain_ids` — the id of every read folded into an observation — is what lets the
  merge build each read's haplotype across the several records one locus can span, and what a later
  step would chain loci with. **Nothing in this model reads it**, and §1.4 says what that does and
  does not cost.

### 1.3 Vocabulary

Six terms do real work below. Three are this project's and three are borrowed.

- **observation** — one distinct sequence the reads showed at a locus, with the number of reads that
  showed it. **Not one read**, and not one allele: the identity is
  `(bases, read witness, read group)`, so the same sequence seen by a read that spanned the whole
  locus and by one that saw part of it is two observations
  ([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs)).
- **complete** and **partial** observation — a read that showed the whole locus pins what the
  sample carries there; a read that entered the locus and ran off its own end proves only that the
  locus is *at least* what it showed. The statistician's word for the second is **censored**: the
  value is known to exceed a threshold, not known. §5 is about them.
- **stutter** — the copying enzyme adding or dropping whole repeats before sequencing, so a
  read reports a tract one repeat longer or shorter than the DNA it came from. At tomato dinucleotide
  tracts a read is nearly five times as likely to have lost a repeat as to have gained one — 2,438
  reads against 501 ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §3).
- **the anatomy of a repeat locus**, fixed here and used nowhere else in this document with any other
  words. The **motif** is the short sequence that repeats — `CA` — and its length in bases is the
  **period**. One copy of the motif inside the run is a **repeat**, and how many there are is the
  **repeat count**. The run of them is the **tract**, and the unique sequence either side of it is
  the **flank**. *(These are the repo's own words: `Motif`, `period()`, `left_flank`/`right_flank`,
  and `tract` throughout. "Repeat unit" is a third name for the motif and this document does not use
  it; nor does it say "track", which is a different word.)*
- **a whole-repeat change** and **a part-repeat change** — a read's length differs from the allele's
  by a whole number of repeats, or it does not. The first is slippage; the second is an ordinary
  small insertion or deletion, is rarer, and gets its own parameters. **These are what HipSTR calls
  *in frame* and *out of frame*, and this document does not use those words** — *frame* is borrowed
  from coding sequence, it says nothing about repeats to a reader who has not met it, and it was read
  here as meaning *inside the tract* against *in the flanks*, which is a different distinction
  entirely. The parameter table in §4.2 carries HipSTR's field names for whoever ports from it.
- **stratum** — one *(motif period, repeat count)* cell. The pre-pass fits the four slippage numbers
  per read group per stratum, because how much a tract slips depends on how many repeats it has more
  than on anything else: 9 reads in 10,000 below four repeats against 2 in 100 at six or more
  ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4, §5).
- **read group** — one `@RG`, one lane, one chemistry. It is the grain at which noise is a
  property of the data ([`parameter_prepass.md`](parameter_prepass.md) §1.1), and §2.3 explains why
  the merge must not collapse it.

### 1.4 What the evidence actually is, and why that decides the shape

The cohort merge hands over, per sample per locus, a list of
observations, and each carries
([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs)):

| field | what it holds |
|---|---|
| `bases` | the sequence those reads showed, projected onto the locus span |
| `read_witness` | whole locus, or which stretch of it |
| `read_group` | which lane |
| `num_obs` | how many reads showed it |
| `chain_ids` | **the id of every read folded here**, one per read or per read pair, the reference-matching ones named too — **generic path only; empty on the STR path**, which does not phase ([`locus_generation/mod.rs:358`](../../../../src/ng/locus_generation/mod.rs)) |
| `q_sum` | the **sum over those reads of the logarithm of each read's error probability** |
| `num_fwd`, `mapq_sum`, `mapq_sum_sq`, `placed_left` | strand, mapping-quality and read-position moments, for the site filters |

**On the generic path every read is named, so a read's haplotype can still be built.** `chain_ids`
carries the id of every read folded into an observation, the ones matching the reference included, so
one cohort locus spanning several of a sample's records can be walked record by record and each
read's own stretch of sequence stitched back together. **On the STR path the list is empty**, and
costs nothing there: a tract is one record, so there is nothing to stitch. That is what the merge's open question about compound alleles turns
on — whether a sample's two nearby changes sat on the same molecule or on opposite ones
([`cohort_merge.md`](cohort_merge.md) §14, question 2).

**What the merge does not keep is each read's quality.** Base and mapping quality arrive summed over
the reads folded together, and the list of names is carried beside them rather than being part of
what separates one observation from another — so nothing says which named read carried which
quality.

**So the constraint on this model is about quality, not about reads.** It may use `num_obs` and
`q_sum` freely. It may not use a **non-linear** function of a per-read quality, because `q_sum`
recovers only the **geometric mean** of the reads' error probabilities — `exp(q_sum / num_obs)` — and
the geometric mean is not the arithmetic mean. On an observation half at Phred 30 and half at Phred
20, the geometric mean is 3.2 in a thousand and the arithmetic mean is 5.5 in a thousand: a factor of
1.7.

**Which mean a term wants depends on the term, and the difference is not always a defect.**
Production's contamination path substitutes the geometric mean for every read's own error
([`posterior_engine.rs:1531`](../../../../src/var_calling/posterior_engine.rs)), and **that is the
right substitution for the term that dominates**, which is worth stating so that nobody "fixes" it. A
read the genotype cannot explain contributes `log ε` plus a constant, so adding the logs and taking
the geometric mean are the same arithmetic — identical to the last bit with contamination switched
off, and 0.14 nats apart under this configuration: 20 reads none of which the genotype explains, a 3%
contamination fraction, and the contaminant carrying that allele at 1 in 1,000. **The gap grows with
how often the contaminant carries it** — 1.13 nats at 1 in 100, 1.89 at 1 in 2. The term that wants the **arithmetic** mean is the one for reads the genotype *does* explain,
`log(1 − ε)`, and it is worth 0.05 nats over those same 20 reads — a fifth of a Phred. **ng does not
carry the substitution because §3.3 drops the term it would have appeared in, not because it is
wrong.**

§2.3 turns this into a contract and §3.3 shows the SNP/indel formula that satisfies it exactly — it
is not a coincidence that freebayes' expression has the shape it has.

---

## 2. One formula, two emissions

### 2.1 The shared formula

**A read arrives one of three ways: it copied one of the individual's own allele copies faithfully
enough, it copied one of them and something went wrong, or it did not come from this individual's
copies of this locus at all.** The third possibility is what stops one strange read from ruling out
a genotype, and every caller has some version of it. Writing that down for a genotype `g` whose copy
counts are `k_a` over a ploidy `P`:

```text
log Lg(g)  =  Σ_o  n_o · log[ (1 − λ) · Σ_a (k_a / P) · Lr(o | a)  +  λ · U(o) ]
```

- `o` ranges over the sample's observations at this locus; `n_o` is how many reads showed `o`.
- `k_a / P` is the chance a read was copied from a copy of allele `a` — the reads are assumed to be
  drawn from the individual's copies in proportion to how many there are.
- `Lr(o | a)` is the **emission**: the probability that one copy of allele `a` produces observation
  `o`. This is the only part that differs between the two paths.
- `λ` is the chance a read did **not** come from this individual's copies of this locus, and `U(o)`
  is what such a read shows.

**The third term is literally two likelihoods added, weighted by how often each happens** — one for
the read having come from this sample, one for it having come from somewhere else — which is
production's shape and freebayes' shape both. Written out at a contamination fraction `c`, it is
`(1 − c) × P(read | this sample's genotype) + c × P(read | the contaminating population)`
([`posterior_engine.rs:1555`](../../../../src/var_calling/posterior_engine.rs)). §2.2 says what
plays the part of `λ` and `U` on each path.

**Two quite different things put a read in that third term, and only one of them cancels.** The
distinction is worth pinning down before either path is read, because the same algebra behaves
differently in the two cases:

- **A read no candidate allele could have produced at all** — junk, a read from a paralogous tract,
  a chimera. Its emission is zero under *every* allele, so the bracket collapses to `λ·U(o)` for
  every genotype alike, the term is the same number in every row, and it **drops out** when the
  caller normalises. That is what stops one strange read from ruling a genotype out and turning the
  sample's whole row into nothing.
- **A read the candidates *can* produce, arriving from a contaminating individual.** Here the
  bracket still contains the genotype through `P(read | this sample's genotype)`, so the term
  **does not cancel** — and it must not. What it does is soften the read's pull: a handful of reads
  carrying an allele this sample does not have is no longer forced to be either an error or
  evidence of heterozygosity, because there is a third explanation with a measured weight.

### 2.2 What each path plugs in

| | SNP/indel path | STR path |
|---|---|---|
| what an observation is | the projected allele sequence over the locus span | the tract sequence the read showed |
| `Lr(o \| a)` | the read is right with probability `1 − ε`, and wrong with probability `ε`, shared out over what a wrong read could show (§3.5) | how likely a slip of this size is, times how well the letters match once the candidate is stretched to the read's length (§4.2, §4.3) |
| `ε` from | the read's own quality, rescaled to the read group's fitted rate (§3.2) | the fitted substitution rate for this read group and stratum ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.1) |
| `λ` from | this read group's contamination fraction, or zero (§3.6) | a fixed outlier weight (§4.5) |
| `U(o)` | the contaminating population's frequency for this allele's class | uniform over what could have been seen (§4.5) |
| what the third term carries | contamination | reads no allele explains, **and** contamination — three components, §4.5.1 |

**The SNP/indel path uses the third term for contamination; the STR path uses it for reads no allele
explains, and for contamination when that is switched on.** The SNP/indel path needs no term for
unexplained reads: one is charged as an error, which is a large negative number rather than a zero,
so no genotype's likelihood collapses and nothing has to be floored. Contamination at a repeat tract
is modelled by no existing caller — not production, HipSTR or GangSTR — and §4.5.1 specifies it,
builds it, and turns it **on by default**, for the reason §3.6 gives.

### 2.3 The aggregation contract

**The formula in §2.1 has a logarithm outside a sum, and that is where an implementation quietly goes
wrong.** `log` of a mixture is not a sum of `log`s, so an observation's `n_o` reads may be pooled into
one term **only if every one of them would have got the same number**. The reads in one observation
already share their bases, their witness and their read group, so the only thing that can differ
between them is their quality — which is why this contract is about quality and about nothing else
(§1.4). Two rules follow, and both are requirements on other components, not on this one:

1. **The merge must keep `read_group` as part of an observation's identity.** Two reads showing the
   same sequence from two lanes have different error rates and must not be pooled. The evidence type
   already keys on it; [`cohort_merge.md`](cohort_merge.md) §4.2 describes a sample's moments being
   summed "where two of its own observations projected onto the same allele", and this document
   requires that the summing stop at the read-group boundary. On a single-library sample — which is
   most of them — this costs nothing. **Built 2026-08-23**
   ([`calling_prerequisites.md`](../impl_plan/calling_prerequisites.md) B1): the merge's rows are
   one per `(allele, read group)`, ascending. The sentence §4.2 of the merge's own spec still
   carries is the one that was corrected here, and it is corrected there too.
2. **Anything else that varies read by read must enter as a term that is exactly additive in
   `q_sum`, or not at all.** §3.3's formula is built to satisfy this and §12's ninth test pins it:
   compute the likelihood from a list of made-up reads and from the aggregate of those same reads
   and require the two to agree **bit for bit**. This one *is* bitwise: §3.3's formula sums the same
   terms in the same order either way, with no round trip through probability space.

**What this contract costs, stated plainly.** A read's mapping quality cannot enter as its own
`λ` — the natural and correct treatment, since a mismapped read is exactly a read that came from
somewhere else — because `log[(1 − λ_r)x + λ_r U]` is not additive in `λ_r` and the merge holds only
the sum and sum-of-squares of mapping quality. So mapping error is folded into the read's error
probability instead, at the moment the observation is minted, which is what production does
([`open_record.rs:793`](../../../../src/pileup/walker/open_record.rs)). §10 records the per-read
mixture as deferred and says what it would cost to enable.

### 2.4 Why this is one document and not two

The two emissions have no mathematics in common. The formula around them does, and so do five decisions that would
drift if they were written twice: how observations are aggregated (§2.3), the copy-weighted mixture
that makes the model ploidy-generic, the treatment of reads that saw only part of the locus (§5), the
no-allocation contract (§8), and determinism. Production writes a copy-weighted mixture wrapped in one logarithm per observation twice over —
[`posterior_engine.rs:1475`](../../../../src/var_calling/posterior_engine.rs) and
[`ssr/cohort/likelihood.rs:75`](../../../../src/ssr/cohort/likelihood.rs) — without either half
knowing about the other. **They are not the same expression**, and §2.2 says how they differ: only
the STR one carries a term for reads no allele explains, only the SNP/indel one carries
contamination. §3 and §4 are each self-contained enough that a
coder implementing one need not read the other.

---

## 3. The SNP/indel emission

### 3.1 What a read's error probability is, and what every other caller does with it

**Two different numbers both answer to the name "error rate", and the owner's question is the right
one: they are not the same quantity.**

- **A base quality** is the sequencing instrument's own claim about one base of one read. It is per
  read, it is the only thing that distinguishes one read from another, and it describes the
  sequencing chemistry and nothing else.
- **The pre-pass's fitted rate `ε`** is the rate at which an admitted read shows an allele the
  individual does not carry, *whatever put it there*
  ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2): the sequencer misreading a
  base, but also a polymerase misincorporating during amplification, two fragments joining into a
  chimera, a base damaged before sequencing, DNA from another individual, a read from a paralogous
  locus mismapped here. It is per read group, it is measured from the data rather than claimed, and
  it cannot tell one read from another.

**The first is precise and incomplete; the second is complete and undifferentiated.** So the question
is not which to believe. It is how to use both.

**What the four reference callers do.** All four vendored callers use the read's own base quality in
the likelihood, and all four refuse to take it as reported. What differs is the device.

| | uses per-read base quality | how it refuses to trust it as reported |
|---|---|---|
| **freebayes** | yes — per observation | combines it with mapping quality into one number, `1 − (1 − ε_base)(1 − ε_map)` ([`DataLikelihood.cpp:77`](../../../../freebayes/src/DataLikelihood.cpp)); then multiplies the whole error term by a **read-dependence factor**, default 0.9, so the error evidence from many reads is discounted about a tenth ([`Parameters.cpp:475`](../../../../freebayes/src/Parameters.cpp)) |
| **GATK** | yes — every cell of the pair-HMM is `1 − ε(Q)` on a match and `ε(Q)/3` on a mismatch ([`LoglessPairHMM.java:89`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/LoglessPairHMM.java)) | **refits the quality from data** (below); squashes qualities below Phred 18 to the minimum ([`PairHMM.java:31`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMM.java)); and caps how much any one read may discriminate between two alleles at Phred 45, on the grounds that the read might have come from elsewhere in the genome ([`LikelihoodEngineArgumentCollection.java:109`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/LikelihoodEngineArgumentCollection.java), applied at [`AlleleLikelihoods.java:441`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/genotyper/AlleleLikelihoods.java)) |
| **bcftools / samtools** | yes | recomputes each base's quality from a local realignment (BAQ — on by default in its *partial* form,
over problematic regions only; `--full-BAQ` applies it everywhere), **caps it at Phred 60** and at Phred 30 to 50 under the per-platform profiles, lowers a base toward its neighbours' quality, and downweights the *n*-th read supporting the same base and strand geometrically — the MAQ dependency coefficient ([`mpileup.c:1382`](../../../../bcftools/mpileup.c), [`errmod.c:77`](../../../../htslib/errmod.c)) |
| **HipSTR** | yes, in its read-to-haplotype alignment, capped at Phred 41 ([`base_quality.h`](../../../../HipSTR/src/base_quality.h)) | — |
| **GangSTR** | **no** — its enclosing-read likelihood is a geometric stutter term with no quality in it at all ([`enclosing_class.cpp:63`](../../../../GangSTR/src/enclosing_class.cpp)) | — |
| **ours, today** | yes — per read, the worse of the lowest base quality across the locus window and the mapping quality ([`open_record.rs:793`](../../../../src/pileup/walker/open_record.rs)) | nothing; the reported quality is used as reported |

**GATK's device is the interesting one, because it is the fitted-rate idea done properly.** Base
quality recalibration counts, over the whole genome, how often a base actually mismatched, split by
**read group × reported quality × sequence context × position in the read**
([`StandardCovariateList.java`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/recalibration/covariates/StandardCovariateList.java)),
and writes the empirical rate back onto each base. So GATK ends up with a per-read number whose
*scale* is measured and whose *shape* is the instrument's. **What it costs is a catalogue of known
variants**: without one, every real heterozygous site counts as a mismatch and the recalibration
learns that the sequencer is far worse than it is. GATK makes that catalogue a required argument,
not an optional one
([`BaseRecalibrator.java:108`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/walkers/bqsr/BaseRecalibrator.java)).
**A tomato panel does not have one**, and neither does any of the species this caller is meant to
serve, which is why ng cannot simply adopt it.

### 3.2 The decision: keep the read's own number, set its scale from the fit

**ng charges each read the error probability the read carries, rescaled by one number per read group
so that the average over that read group's admitted reads equals the rate the pre-pass measured.**

```text
ε_read(after)  =  ε_read(as minted) × scale(read group)

                       fitted error rate for this read group
scale  =  ────────────────────────────────────────────────────────────────
          geometric mean, over that group's admitted reads, of ε as minted
```

> **Written per read, applied per observation — the caller never holds a per-read `ε` at all.**
> The merge keeps, for each allele in each read group at each locus, how many reads support it and
> the **sum of their log error probabilities**; the reads themselves are gone from there on (§1.4).
> From those two numbers exactly one average is recoverable — `exp(q_sum / num_obs)`, the geometric
> mean — and that is what §3.3 charges. So the line above states an intent, not a mechanism: in the
> code the scale is **one addition of `log scale` per observation**, and nothing is multiplied read
> by read.
>
> **The two are the same arithmetic**, which is why the line above is still true of the result:
> `exp(Σ ln(s·ε) / n) = s · exp(Σ ln ε / n)`. Scaling every read and scaling their geometric mean
> are one operation, so the aggregation costs nothing.
>
> **And this is the sharper form of the correction recorded below (owner, 2026-08-24).** An earlier
> version of this section asked for the **arithmetic** mean of the per-read `ε`. That was not a
> worse choice than the geometric mean — **it was not a choice at all**, because nowhere in the
> model does a per-read `ε` survive to be averaged that way. Reading the sentence above as a
> per-read mechanism is what made it look reachable.

**Why this and not either half alone.** The base qualities carry the only information that
distinguishes one read from another, and at three reads a position that distinction is the whole
call: a locus with one alternative read at Phred 40 and one at Phred 13 are not the same evidence,
and a single fitted rate says they are. The fitted rate carries the only information about the
error sources the instrument cannot see — mismapping, chimeras, index hopping, damage — and it is
the number that was actually measured on this data rather than asserted by the machine. The scale
keeps the shape from the first and the size from the second.

**This is GATK's benefit without GATK's requirement, and the reason it is available to us is the
pre-pass's central decision.** GATK needs a variant catalogue because it identifies mismatches by
masking known sites; the pre-pass instead sums over the unknown genotype rather than choosing one
([`parameter_prepass.md`](parameter_prepass.md) §3), so it measures the same rate with no catalogue
at all. The recalibration is coarser than GATK's — one scale per read group against GATK's four-way
covariate table — and that is the honest trade: the covariates GATK conditions on need per-base
counts this pre-pass does not keep, and [`parameter_prepass_generic.md`](parameter_prepass_generic.md)
§2 works through why a quality axis in its histograms is not affordable.

**It cannot double-count, and the reason is that the scale is a ratio.** Production's minted per-read
number already contains mapping error, because it takes the worse of base and mapping quality. If
mapping error is already in the per-read numbers, their mean is already nearer the fitted rate, so
the scale comes out nearer one and adds nothing further. Whatever the per-read number already
captures, the scale supplies only the remainder.

**What it asks of the pre-pass — two numbers per read group.** A running sum of the per-read
**log** error probability, and the count of reads it was summed over. Two scalars per read group,
no new traversal, and **already carried**: an observation's `q_sum` is that sum over its own reads
and `num_obs` is that count, so the accumulator adds up numbers the walk has already produced
rather than minting anything.

> **The average is the geometric mean, and an earlier version of this section asked for the
> arithmetic one (corrected 2026-08-24, owner).** It said "a running sum of the per-read error
> probability", which is a sum of `ε` and not of `ln ε`, and nothing carries that — the walk sums
> the logarithms and throws the individual reads away, and `Σ ε` cannot be recovered from `Σ ln ε`.
> Supplying it would have meant a second accumulation at fold time and a new field on every
> observation.
>
> **Taking the geometric mean instead is not a concession, it is the self-consistent choice**, and
> the reason is what the scale is applied *to*. The model charges an observation
> `exp(q_sum / num_obs)` — a geometric mean — and production does the same, clamped
> ([`posterior_engine.rs:1536`](../../../../src/var_calling/posterior_engine.rs); production has no
> recalibration at all, so there is nothing there to copy but the quantity). A scale built from an
> arithmetic mean and applied to a geometric one would not make the calibrated property hold in the
> model's own terms, so paying for the arithmetic sum would buy an inexactness rather than remove
> one.
>
> **The two are not the same number, and they are twenty-five to forty-four times apart on real
> reads** (measured 2026-08-24, `examples/ng_minted_error_means.rs`; report
> [`ng_prereq_closeout_two_averages_2026-08-24.md`](../../reports/implementations/ng_prereq_closeout_two_averages_2026-08-24.md)).
>
> | | tomato, 63 accessions, 2.5× to 28.6× | HG002, 100 benchmark regions, 300× |
> |---|---|---|
> | read-positions | 5,485,730,235 | 172,616,054 |
> | geometric mean of `ε` | 5.982 × 10⁻⁴ (Phred 32.2) | 2.905 × 10⁻⁴ (Phred 35.4) |
> | arithmetic mean of `ε` | 1.505 × 10⁻² (Phred 18.2) | 1.282 × 10⁻² (Phred 18.9) |
> | ratio | **25.2** | **44.1** |
>
> Per read group on tomato the ratio runs 22.7 to 37.0, median 24.4 over 63 accessions — **no read
> group anywhere near one**. Building the scale from the arithmetic mean would therefore have
> divided every charged error by 25 to 44: **14 to 16 Phred**, every read treated as that much
> cleaner than the pre-pass measured it to be.
>
> **And the arithmetic mean is not measuring the chemistry.** A read's minted error is the *worse*
> of its base and mapping qualities, and the mate-overlap rule silences the losing mate of an
> overlapping pair by giving it base quality Phred 0 — an error probability of exactly one, on a
> read that still counts. `ln 1 = 0`, so such a read adds nothing to the log sum and a whole unit to
> the probability sum. Measured: **9 read-positions in 1,000 on HG002 carry `ε = 1`, and they are
> 73% of Σ ε**; on tomato, 7 in 1,000 and 47%. So the arithmetic mean is mostly a measurement of how
> often mates overlap. That is a second, independent reason for the geometric mean, and it does not
> depend on the self-consistency argument above.

Two requirements on them, and the second is easy to get wrong:

- **The per-read quantity summed must be computed by the same function the locus generator mints
  with**, or the scale calibrates against a different definition of "how wrong is this read" than the
  one it is applied to. §12's tenth test calls that function from both sides on the same read.
- **The sum must run over exactly the reads the surviving error-rate estimate was fitted from.** The
  error rate is one of the three quantities the pre-pass fits **twice** — once after each sample's
  walk from the read-group histogram, and once at the gather from the census sites — with which
  estimate is kept still open ([`parameter_prepass.md`](parameter_prepass.md) §1.3, §4.1). *An
  earlier version of this section named only the histogram, which was wrong.* **So the accumulator is
  per route and both routes need one**, exactly as both routes are built today: the histogram route
  sums over the sites its histogram counts, the census route
  ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.1) over the census sites its
  fit reads. Whichever route §4.1's comparison keeps, its accumulator goes with it and the other is
  deleted with the other fit. **A scale whose numerator and denominator come from different site sets
  is not a calibration**, because the two sets do not have the same quality profile — the census
  sites are chosen and the histogram's are not.

  **The per-position depth cap does divide them, and by 3 parts in 100 on the deepest real sample
  we have.** The histogram route thins every position to at most
  [`MAX_BINNED_DEPTH`](../../../../src/ng/parameter_estimation/generic/depth_bins.rs) = 124 reads
  before fitting; the accumulator thins nothing. *An earlier version of this paragraph said the two
  were undivided, and gave the per-site argument as the reason. The per-site argument is right and
  it is not the whole of it.*

  **Per site, the cap is harmless**: the draw is hypergeometric on counts and never looks at a
  read's quality, so the mean log error over the kept reads has the same expectation as over all of
  them. **Across sites, it re-weights.** A 500-read position casts 500 votes in the denominator and
  124 in the population the numerator was fitted from, so the two averages weight deep positions
  differently — and deep positions are not a random sample of the genome, they are where reads pile
  up from elsewhere and where mapping quality collapses.

  **The size, measured rather than argued** (`examples/ng_minted_error_means.rs`, 2026-08-24). On
  HG002's 100 benchmark regions at 300×, where the cap fires at essentially every position — the fit
  sees 70,288,390 of 172,616,054 read-positions, 41 in 100 — the denominator's geometric mean is
  2.9055 × 10⁻⁴ against 2.9862 × 10⁻⁴ with each position thinned to the cap first. **The
  re-weighting moves it by 2.7%, which is 0.12 Phred.** On the 63-accession tomato cohort at 10× to
  28.6× it moves it by nothing at all: on the deepest accession 228,468,065 of 228,492,796
  read-positions are under the cap, and the mean moves by a factor of 1.0000.

  **So the divergence is real, bounded at the top of the depth range, and unmeasured beyond 300×.**
  Two ways to close it, and the choice is the owner's: thin the accumulator at the same cap, which
  makes the two counts identical and costs a multiply per site; or leave it and accept 3 parts in
  100 at 300×, on the argument that the population the *scale* is applied to at calling time is
  every read and not a thinned subsample. **Nothing decides this until the scale has a consumer**,
  which is [`calling_read_likelihoods.md`](../impl_plan/calling_read_likelihoods.md) A2.

  **The census route cannot supply either number as it stands**, and that is a fact about its
  records rather than about this requirement: its per-position unit is a depth code and a sparse
  list of non-reference allele counts
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md)), with no quality in
  it at all. So its accumulator waits on §4.1's comparison — if that route wins, the records gain a
  quality field; if the histogram route wins, nothing is owed.

**This does not reopen the pre-pass's decision about base qualities, and it is worth saying why.**
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2 decides that base qualities are
not modelled *inside the fit*, and its argument is exact: letting `ε` vary by base quality destroys
the `(depth, alternative count)` sufficient statistic the whole estimator rests on. Nothing here
varies `ε` inside the fit. The scale is one scalar computed after the fit, from an accumulator that
sits beside the fitted object rather than inside it, and it is used only downstream. Both statements
are true at once.

*(One sentence elsewhere needs its scope narrowed rather than its meaning changed.
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.1 rules out the third of
`verifyBamID2`'s sources of upward bias in the contamination estimate on the grounds that "ng fits
one error rate and does not read base qualities at all". That stays true of **the fit**, which is
what biases the estimate, and stops being true of **the caller**. The edit is not made here.)*

**At the ends of the range.** The fitted rate is pooled over hundreds of millions of base
observations, so it is equally well determined at 3 reads a position and at 300, and it exists at
one sample as readily as at a thousand — it is a property of a read group, and one sample has read
groups. **Where the pre-pass emits no rate** — too little data, or a sample whose noise is off the
end of its ladder, which two of five real alignments were
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1) — the scale is 1 and the
qualities are used as reported. **That must be visible in the run's output**, because a run
calibrated against a measurement and a run trusting the instrument are otherwise indistinguishable.

### 3.3 The formula

**A read either shows something the genotype can produce, in which case it is charged only for which
copy it came from, or it does not, in which case it is charged for being wrong.** For a genotype `g`
with copy counts `k_a` over ploidy `P`, at a locus whose allele table has `A` alleles:

```text
log Lg(g)  =   Σ         n_o · log( k_{a(o)} / P )                     ← the reads it explains
            o : k_{a(o)} > 0

            +  Σ         [ q_sum_o  +  n_o · ( log scale_{r(o)} − log m(a(o), g) ) ]
            o : k_{a(o)} = 0                                            ← the reads it calls errors

            +  q_sum_other                                              ← the pooled leftover
```

- `a(o)` is the allele the observation's bases match, `r(o)` its read group.
- `q_sum_o` is the observation's summed log error, straight off the evidence.
- `scale` is §3.2's per-read-group calibration; in log space it is one addition per observation.
- `m(a, g)` is how many things a wrong read could have shown, which §3.5 settles.
- `q_sum_other` is the pool of reads matching no allele in the table. It is the same for every
  genotype so it cancels in genotyping, and it is **kept** because the data likelihood also feeds
  emission and QUAL, where an absolute value is compared between loci.

**A sample with no reads gives zero for every genotype and needs no branch** — an empty sum is
zero — so such a sample is decided by the prior alone, which is the right answer and not a special
case.

**Two approximations live in this formula and both have a size.**

- **A read that the genotype explains is not charged for being right.** The exact term is
  `log(k_a/P) + log(1 − ε)` and the second half is dropped, because recovering it needs the
  *arithmetic* mean of `ε` over the observation's reads and `q_sum` gives the geometric mean (§1.4).
  The omission is at most `n·ε` nats for a genotype explaining `n` reads at rate `ε`, and what
  matters is the *difference* between two genotypes, which is the difference in how many reads each
  explains. At 20 reads with 2 of them alternative, at an error rate of 1 in a thousand, that
  difference is 0.002 nats — a hundredth of a Phred. At 300 reads and a poor library at 1 in 20, comparing a heterozygote
  explaining all 300 against a reference homozygote explaining 285, it is 0.75 nats, about 3 Phred.
  **So it is negligible at good chemistry and small but real at bad chemistry and high depth**, and
  it always favours the genotype that explains more reads.
- **Every read supporting one allele is charged the same, whatever its own quality**, once it is on
  the error side. That is exactly what `q_sum` encodes and it is not an approximation at all: the
  sum of logs is the log of the product, which is what a labelled-read likelihood wants.

**How this compares with freebayes, term by term.** freebayes' default branch is the closest thing
in any other caller to the formula above, and the two agree on more than they differ on. Setting
them side by side is the quickest way to see what is genuinely new here and what is inherited
([`DataLikelihood.cpp:135`](../../../../freebayes/src/DataLikelihood.cpp),
[`:138`](../../../../freebayes/src/DataLikelihood.cpp),
[`:162`](../../../../freebayes/src/DataLikelihood.cpp)):

| | freebayes, default branch | this document |
|---|---|---|
| a read the genotype explains | `log(k_a / P)`, no `(1 − ε)` factor | **the same**, and §3.3's first approximation is the same omission |
| a read the genotype cannot explain | `log[ 1 − (1 − ε_base)(1 − ε_map) ]` | `log[ max(ε_base, ε_map) ] + log scale − log m` |
| combining base and mapping error | a **sum** — the read is wrong if either is | the **worse of the two**, which is production's mint |
| which base's quality, on a multi-base allele | the **average** over the allele's bases, and the **minimum** across alleles merged into one haplotype ([`Allele.cpp:158`](../../../../freebayes/src/Allele.cpp), [`AlleleParser.cpp:3154`](../../../../freebayes/src/AlleleParser.cpp)) | the **minimum** across the whole locus window ([`open_record.rs:944`](../../../../src/pileup/walker/open_record.rs)) |
| where the error mass goes | nowhere — the full `ε` is charged | divided by 3 for a plain substitution (§3.5) |
| calibration against a measured rate | none | the per-read-group scale (§3.2) |
| pooled-error discount | the read-dependence factor, 0.9 | none — deferred (§10) |
| multinomial coefficient | none in this branch; present in `--legacy-gls` | none (§3.4) |

**Two of those rows are differences of size rather than of kind, and both are small.** Summing the
base and mapping error rather than taking the worse of the two differs by at most a factor of 2 —
0.69 nats, 3 Phred — and only when the two are within a few times of each other; where one dominates,
which is the ordinary case, they agree to within a percent. Taking the minimum quality across the
locus window rather than the average within an allele is the more conservative of the two, charging a
larger error probability and so making an unexplained read cheaper to dismiss, which favours the
reference; how much depends on how variable quality is across a locus's bases and has not been
measured here.

**The rows that are differences of kind are the scale, the divisor and the discount** — the first is
new, the second follows the pre-pass's own noise model rather than freebayes, and the third is
something freebayes has that we do not.

### 3.4 The multinomial coefficient production inherited — dropped

**Production's SNP/indel likelihood carries a term this formula does not**
([`per_group_merger.rs:1948`](../../../../src/var_calling/per_group_merger.rs)):
`log(N!) − Σ_a log(n_a!)` summed over the alleles the genotype carries, where `N` is how many reads
those alleles account for. It comes from freebayes — but from freebayes' `--legacy-gls` branch
([`DataLikelihood.cpp:155`](../../../../freebayes/src/DataLikelihood.cpp)), and **freebayes' default
branch has no such term**, returning a plain sum of per-observation logarithms
([`:162`](../../../../freebayes/src/DataLikelihood.cpp); the default is set at
[`Parameters.cpp:459`](../../../../freebayes/src/Parameters.cpp)).

**It is not a normalisation that either model wants.** Over labelled reads — which is what the
evidence is, since each read carries its own error probability — there is no combinatorial
coefficient at all. Over unlabelled counts there is one, but it runs over *every* allele including
the ones called errors, is therefore the same for every genotype, and cancels. Production's runs over
a subset that changes with the genotype, so it is neither, and it does not cancel.

**What it does, with its size.** For a diploid heterozygote at `n` reads split with a fraction `f` on
one allele, the coefficient adds about `n · H(f)` nats, where `H` is the binary entropy; for a
homozygote it adds nothing. Worked at 20 reads with 2 of them alternative at Phred 30, comparing a
reference homozygote against a heterozygote and including the genotype prior at a diversity of 1 in
a thousand. **Both rows charge a wrong read the full `ε`, which is production's setting and not this
document's, so that the coefficient is the only thing that differs between them:**

| | reference homozygote favoured by |
|---|---|
| without the coefficient | 6.96 nats — **30 Phred** |
| with the coefficient (production) | 1.71 nats — **7 Phred** |

**Two alternative reads in twenty at Phred 30 is the canonical false-positive shape**, and the
coefficient turns a clear reference call into a coin flip. At a balanced heterozygote the coefficient
instead adds 12.1 nats at 20 reads and 205 nats at 300 — 53 and 890 Phred of extra confidence in a
call that was already going to be made.

**It contributes to the recorded depth inflation of QUAL without being the mechanism behind it**, and
the difference matters because dropping the coefficient will not cure that defect. The recorded
diagnosis is the *error* term: false-positive sites carry a persistent alternative fraction of about
20 in 100, so their alternative-read count rises linearly with coverage and the heterozygote's
advantage is dominated by `N_alt · (−ln ε)`, linear in depth — the false-positive QUAL median running
1 → 3 → 150 from 5× to full
([`../reports/qual_fp_depth_inflation_2026-06-10.md`](../../reports/qual_fp_depth_inflation_2026-06-10.md)).
The coefficient's `n·H(f)` sits on top of that and is the smaller part at a 20-in-100 fraction.
**The shipped cure was an allele-balance penalty at VCF encode** — which is where the parenthesis
below says allele balance belongs.

*(One property is genuinely lost by dropping it, and it should be named rather than discovered later.
Because `n·H(f)` is largest at `f = ½`, the coefficient makes an unbalanced heterozygote score worse
than a balanced one — a 60:40 split at 300 reads loses 6.0 nats, 26 Phred, against 50:50. That is a
real allele-balance test, and it is arriving through a term with no probabilistic reading. Allele
balance belongs in the site filters, where production already computes it
([`vcf/qual_refine.rs`](../../../../src/vcf/qual_refine.rs)), and where it can be
calibrated. It does not belong in the likelihood by accident.)*

**Decision: ng does not carry the coefficient.** This changes genotypes against production, so it is
not free: §12's thirteenth test runs both forms over the same merged records and reports which
genotypes move, so the change is attributed rather than asserted.

### 3.5 Where a wrong read's probability goes

A read that the genotype cannot produce is wrong. **But wrong how?** If the individual carries `A` at
this base and the read shows `C`, the chance of that particular misread is not the chance of any
misread — there were three bases it could have gone to.

The three callers answer differently. GATK divides by three
([`LoglessPairHMM.java:90`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/utils/pairhmm/LoglessPairHMM.java));
freebayes and production divide by nothing; production's own STR substitution term divides by three
([`pair_hmm.rs`](../../../../src/ssr/cohort/pair_hmm.rs)), so the two halves of production disagree
with each other. The pre-pass's noise model is explicit that the denominator is three — *"three bases
to go wrong into, one to come back to"*
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2).

**The size is `log 3` per wrong read — 1.10 nats, 4.8 Phred** — and it does not cancel, because how
many reads a genotype calls wrong varies by genotype. Dividing by three makes a wrong read *less*
probable, so calling a read wrong costs more and the divisor **favours the heterozygote**: at §3.4's
2-in-20 example it moves the comparison 2.2 nats, about 10 Phred, that way — the same order as the
coefficient and in the opposite direction.

**Together the two changes take that configuration from production's 7 Phred for the reference
homozygote to 21 Phred**, which is where the two specified choices land it: 4.76 nats, with the
dropped coefficient worth +5.25 nats towards the reference and the divisor worth −2.20 back.

**Decision: `m(a, g) = 3` when the observation differs from every allele the genotype carries by
substitution at exactly one position, and `m = 1` otherwise.** The first case is a plain substitution
and three is the physical fact. The second case includes insertions, deletions and multi-position
differences, where there is no finite set of things a wrong read could have shown and any divisor
would be invented; leaving the mass unspread is the conservative choice, and conservative here means
favouring the reference, which is the direction a caller should err in when its model runs out.

**`m` is a property of the allele pair, not of the read**, so it is computed once per
`(allele, genotype)` and costs nothing per read. The classification uses the projected sequences the
merge already unified, which is why unification by exact byte match is sound only because indels were
left-aligned upstream ([`cohort_merge.md`](cohort_merge.md) §4.2).

### 3.6 Contamination

**Some of a sample's reads come from someone else's DNA**, and the pre-pass estimates how many —
**per read group**, from the census sites
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5). **Per read group and not per
sample** (owner, 2026-08-19): contamination is DNA from another individual mixed into a *library*,
and two libraries of one sample can carry different amounts of it, so it takes the grain
[`parameter_prepass.md`](parameter_prepass.md) §1.1 gives every property of the chemistry. *(That is
finer than per sample, not a departure from it: the earlier ruling that contamination is a property
of the sample rather than of the locus (§3.6 below) is about a different axis and both hold.)*

**The estimator produces that grain, and did not until 2026-08-20.** ng fits a contamination fraction
inside the joint route
([`joint/contamination.rs`](../../../../src/ng/parameter_estimation/joint/contamination.rs), wired
into that fit, with an `estimate-contamination` subcommand that hands the same estimator to somebody
comparing methods). It is still off by default. **It gave one number per sample and now gives one per
read group**, as do [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.1 and §3.4;
only [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5's route defers an estimator at
all, and that is the other of the two routes §4.1 there is comparing.

**So this section no longer records an approximation.** It used to say: *a per-sample estimate is
applied to every read group of that sample* — exact where a second seedling was in the tube before
the libraries were split off it, wrong where a neighbouring library on the run hopped its index into
one of them. Both causes are real, only the second can make two libraries of one plant differ, and
only the read-group grain can express it. The fraction this document consumes is now fitted at the
grain it consumes.

**An objection this section used to make was wrong, and is kept here so nobody makes it again.** It
argued that splitting a sample's reads by read group *"would hand every fit less data"*, citing the
measured cost of partitioning a panel — about **+0.015** on every sample's fraction, enough to put 41
to 47 of 50 clean samples over a 1% threshold. **That measurement is about splitting the *panel* to
estimate allele *frequencies* from a twelfth of it, and nothing in this change splits the panel.**
Ancestry belongs to the individual, so each sample's coordinates and its fitted frequency are still
computed from all of its reads and from every sample in the cohort; the only thing that takes the
finer grain is the read count one fraction is fitted from
([`../reports/contamination_grain_decomposition_2026-08-20.md`](../reports/contamination_grain_decomposition_2026-08-20.md)
decomposes it parameter by parameter, with the code behind each row).

**And the split costs less than the depth it appears to give away**
([`../reports/contamination_read_group_grain_2026-08-20.md`](../reports/contamination_read_group_grain_2026-08-20.md)).
A plant with two libraries, one carrying 6% stray reads and one clean, returns **0.0628 and 0.0008**
per library against **0.0307 for both** at the old grain. And a library holding three reads a position
returns 0.026 when it is the plant's only library, 0.046 when it is half of a six-read plant and 0.057
when it is a quarter of a twelve-read plant, against a planted 0.060 — the same reads in the library
being measured each time. What limits the estimate is how well the plant's genotype and the panel's
frequencies are known, and neither of those changes grain. **A plant sequenced from one library
returns the identical number at either grain**, which is every sample of every benchmark cohort here.

**What this document must not lose.** `c` below can now be four different things and a consumer cannot
act on them alike: fitted from that read group's own reads; fitted from every read of the plant and
copied onto it, which is what a plant with one library gets and what the sample grain gives; or a
number near zero because **nothing could be measured**, which is not the same claim as *measured
clean*. The parameters carry, beside each fraction, how many markers that read group had a read at,
how many reads it had there, and which of the first two it was. **A library with too little evidence
returns a fraction near zero rather than a refusal** — the likelihood barely moves with `c` and the
search keeps zero, which is the right default for a value this term multiplies — and those counts are
what tell it from a library that was measured and found clean.

That is §2.1's third term with a real `U`:

```text
log Lg(g)  =  Σ_o  n_o · log[ (1 − c) · own(o | g)  +  c · q(o) ]
```

- `c` is this read group's contamination fraction, zero when the parameters fit emits none.
- `own(o | g)` is what §3.3's formula computes per read, as a probability rather than a logarithm:
  `k_a/P` for an explained observation, `ε̄ / m` for a wrong one.
- `q(o)` is **the contaminating population's frequency of the allele this observation shows, at this
  locus** — read off the caller's own frequency estimate for the samples sequenced in the same batch
  as this one, and recomputed every iteration as that estimate moves (owner, 2026-08-24).

> **`q(o)` used to be three numbers — the frequency of an allele *class*, reference against
> substitution against insertion-or-deletion — averaged over the census sites and emitted by a
> side-pass in the parameter pre-pass. That is deleted (owner, 2026-08-24), and the reason is worth
> keeping, because the old shape looks reasonable.**
>
> **The class split is not a modelling requirement. It is what a caller does when it has no
> per-locus frequency to use.** Production splits by class because it has no population frequency
> for an arbitrary alternative allele and has to fall back on a class average
> ([`var_calling/contamination_estimation.rs`](../../../../src/var_calling/contamination_estimation.rs),
> `q_b_per_batch`, a three-entry simplex per sequencing batch). **This caller has the frequency**:
> it is the same per-locus allele frequency the genotype prior reads and the loop re-estimates, so
> the honest `q(o)` is that number and the classes disappear with the ignorance that motivated them.
>
> **Three consequences, all of them simplifications.** The pre-pass owes nothing here — no
> side-pass, nothing new in the census, and the census's fifth allele code (which lumps indels with
> `N` and spanning deletions and so cannot cleanly answer "how often is a contaminant read an
> indel") never has to. There is no new estimator to design, which matters because ng's
> contamination fit produces no per-read posterior that a read came from the contaminant, so
> production's estimator could not have been ported anyway. And the answer is per locus rather than
> per run, which is strictly more information at the place it is used.
>
> **One consequence that is not a simplification, and it contradicts a sentence below.** §6.1's
> first tier says the contaminated likelihood is computed once per sample per locus and read
> unchanged by every iteration. **That stops being true**: the fraction `c` stays frozen, and
> `q(o)` moves with the loop. This does not reopen the ruling that *contamination* never enters the
> loop — that ruling is about `c`, which is a property of the library and of nothing else. It is the
> second half of the mixture that is a property of the locus, and a per-locus quantity is what a
> per-locus loop is for.
>
> **Which samples' frequency.** The batch, not the whole cohort: a contaminant is most likely a
> neighbour on the same run, and the run already declares who was sequenced beside whom
> ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §1.6, whose default is one
> batch holding everything — so a run that declares no batches gets the cohort frequency and loses
> nothing it had).
>
> **Owner:** the read likelihood builds this; nothing is owed by the parameter pre-pass.

**At `c = 0` this is §3.3 exactly, and that is worth checking rather than asserting.** For an
observation the genotype explains, `own` is `k_a/P` and the term is `n_o·log(k_a/P)` — §3.3's first
line. For one it cannot, `own` is `ε̄/m`, and because `ε̄` is `exp(q_sum/n_o)` by construction,
`n_o·log ε̄` **is** `q_sum_o` — §3.3's second line. **The identity is exact algebra and the arithmetic
is not**: §8 evaluates the mixture in probability space and takes one logarithm, so the round trip
through `exp(q_sum/n_o)` and back costs a few units in the last place. **So there is no second code
path and no discontinuity to guard, and the agreement is to a stated tolerance rather than bitwise**
— a few ulp, which §12's eleventh test fixes. *(Production does keep a `c = 0` branch, and this is
why it can afford to: its two forms differ algebraically, so it has a real discontinuity to jump
rather than a rounding one.)* *(Production has one, and falls back to its precomputed
non-contaminated value whenever `c` is zero rather than relying on the two agreeing
([`posterior_engine.rs:1509`](../../../../src/var_calling/posterior_engine.rs)) — because its `own`
carries an extra `(1 − ε)` factor and divides the error mass by the allele count less one, floored at
one, so its two forms
genuinely differ. ng's do not.)*

**What turning it on does cost, with its size.** Once `c` is above zero the logarithm sits outside a
sum, so the reads pooled into one observation can no longer each carry their own error probability
and §1.4's geometric mean is substituted for them. Measured on 20 reads none of which the genotype
explains, half at Phred 30 and half at Phred 20 — a deliberately wide spread — at a 3% contamination
fraction, with the contaminant carrying that allele at 1 in 1,000: **0.14 nats, six-tenths of a
Phred**. It is zero at `c = 0`, it grows with `c`, and **it grows with the contaminant's own frequency
for the allele — 1.13 nats at 1 in 100 and 1.89 at 1 in 2.** So it is small where a contaminant allele
is rare and not where it is common.

**Decision: on by default, wherever the parameters fit emits a contamination fraction above its own
floor.** The estimator that produces it is built and currently switched off (§3.6 above), so what
this decides is what the caller does with the number once that estimator is turned on; where no
fraction is emitted, `c` is zero and the formula is §3.3's at every locus. **A read group whose
fraction rests on almost no evidence also arrives near zero** rather than as a refusal (§3.6), so the
gate is on the value and the evidence counts beside it, never on the value alone. *An earlier version of this section said off by default, "exactly as production has it",
which was a habit rather than an argument.* Three reasons to switch it on. The fraction is one of the parameters fit's named
outputs, and an output the caller declines to consume is one that should not have been produced. The failure it prevents — a contaminated sample called
heterozygous for its contaminant's allele — is a genotype error we know how to avoid. And it is free
where it is not needed: a read group the pre-pass puts at or below its floor gets `c = 0` and the
formula above *is* §3.3, so a clean cohort is untouched by the default.

**The one direction to watch is an overestimated fraction**, which suppresses real heterozygotes by
attributing their alternative reads to the contaminant. `verifyBamID2`'s recorded failure runs the
other way — a true 10% contamination returned as 2.9% when the frequencies came from the wrong
population ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.2) — so the known
failure mode is the safe one here. **The run's output must still carry the fraction used, per
sample**, because a genotype computed at `c = 0.03` and one at `c = 0` are otherwise
indistinguishable.

**The repeat-tract path defaults on too**, for the reason that decides both: contamination is a
property of the sample and not of the marker, so a library carrying another individual's DNA carries
it at repeat tracts as well. What differs there is only how well the *second* half of the mixture is
known — the locus's own allele frequency here, a stand-in built from the prior's seed there — and
§4.5.1 bounds what that costs.

**At one sample there is no contamination estimate** — it is a comparison between samples
([`parameter_prepass.md`](parameter_prepass.md) §1.1) — so `c` is absent and §3.3 is what runs. That
is *emit it as absent*, not *silently fit a zero*.

**The fraction belongs here and not to the EM loop; the frequency beside it does belong to the
loop.** `c` *is* estimated from the cohort, but it is estimated **before** calling and frozen: one
constant per library, unchanged however the caller iterates. **`q(o)` is the other half of the
mixture and it is not a constant** — it is the contaminating population's frequency for the allele
this observation shows, at this locus, which is the loop's own estimate and moves with it
(corrected 2026-08-24; the paragraph this replaces said both halves were frozen and that the
contaminated likelihood was computed once per sample per locus).

**So the two halves sit in different tiers, and that is the whole of the change.** The fraction is
§6.1's first tier, frozen before the loop starts, exactly as production fills its mixture table once
and never touches it again ([`posterior_engine.rs:2361`](../../../../src/var_calling/posterior_engine.rs) — the
mixture branch's own comment; `let did_mixture` is on 2367).
The frequency is recomputed per iteration alongside the frequency the genotype prior already reads,
which costs a lookup rather than a fit. **Production's arrangement is not available to us and does
not need to be**: it freezes both because its second half is a three-entry class average with
nothing per-locus in it, and ours is the locus's own number.

**Decided, and it closes a door rather than leaving one open (owner, 2026-08-19): contamination
never enters the EM loop.** How contaminated a sample is is a property of that sample — of its
library and of how it was handled — and not of any locus, so there is nothing about it for a
per-locus loop to re-estimate.

**The one half that could in principle have varied is the contaminant's allele frequency**, since a
contaminating read shows whatever the contaminating population carries *at this locus*, and the
cohort's own per-locus frequency is exactly what the EM rewrites each iteration. **That route is
closed by the same ruling.** Tying a sample-level quantity to a number the loop keeps rewriting would
make the contamination term move locus by locus and pass by pass, for a refinement nobody has
measured, and it would cost the model the property that makes it cheap — that its numbers survive
every pass of the caller's loop. Production uses one frequency per allele *class* — reference,
substitution alternative, insertion-or-deletion alternative — averaged over the census sites, and
freebayes likewise uses per read-group constants
([`Contamination.h`](../../../../freebayes/src/Contamination.h)). ng follows both.

### 3.7 Positions that are not the kind of position the model assumes — and why the repair is not here

**Everything so far assumes that every position of the genome goes wrong at the same rate** — one
rate per read group, applied wherever that read group's reads land (§3.2). **Some positions do not,
and not because of the sequencer.** A stretch the sample carries twice and the reference once
collects two copies' reads at a single place, so a fixed share of them disagree with the reference at
every depth; a read from a similar sequence elsewhere in the genome lands here; the local sequence is
hard to read or hard to map.

**How much that misses has been measured, on HG002's confident regions.** At the 550,976 positions
the benchmark records no variant of any kind, **818 carry three or more alternative reads where one
rate predicts 29**. Those 818 have to be explained somehow, and the only genotype in the model that
can explain them is the heterozygote — so the surplus came back as heterozygosity, at **1.41 times
the benchmark's count** ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1).

**Decision: this model does not carry a class for them, and the reason is that it cannot see the one
thing that would separate them from real heterozygotes** (owner, 2026-08-19). A read likelihood is
handed one locus. From the reads at that locus alone, a position where half the reads disagree
because the sample carries two copies and a position where half the reads disagree because the sample
is heterozygous are **the same distribution** — the parameters fit says so in as many words, which is
why it refuses to widen its own noisy class toward a half
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1). Adding a latent class here
would ask the likelihood to tell apart two things its evidence cannot distinguish, and it would pay
for the attempt at every ordinary locus.

**What can tell them apart is depth, and depth is not this model's to see.** A collapsed duplication
carries about twice the sample's own median coverage for that GC content. Depth is roughly
independent of allele frequency, so selecting on it does not select on the thing being measured —
which is the same argument [`calling_priors.md`](calling_priors.md) §4.1 makes for choosing census
positions, and which the project has already proved the hard way from the other side: **a class that
identifies these positions by how heterozygous they look took the tomato panel's median accession
from 0.867 heterozygous positions per kilobase to 0.064 — 93% of them — to clear an artefact worth
about 11%** ([`../reports/duplicated_class_on_real_reads_2026-08-14.md`](../reports/duplicated_class_on_real_reads_2026-08-14.md)).

**Home: the site filter that runs after genotyping** — step 11a of
[`ng_proposal.md`](ng_proposal.md), where production already hard-drops on a hidden-paralog
likelihood ratio and where coverage is available. That step sees the whole genome's depth and this
model sees one locus, so the signal that works lives there and not here.

**What this model consumes instead is the single rate the parameters fit already emits**, the
share-weighted average of its clean and noisy rates — the chance a read disagrees with the reference
at a position drawn at random. That is the number §3.2 calibrates against, so nothing above changes.
**What it costs is stated rather than hidden:** at a position of the noisy kind the model is
over-confident, and the 818 in 550,976 above are the size of the class on one high-coverage human
sample. Whether that reaches the *calls* — as against the parameter fit, where it was measured — is
what the site filter's own validation has to answer.

*(Two different classes have been proposed in this project and it is worth not confusing them. The
**clean/noisy pair** is fitted and kept: three numbers per read group, and the parameters fit uses
them internally to stop the surplus above landing on heterozygosity
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.1). The **duplicated class**
is a separate, more aggressive one that was measured on real reads and **ships off**, for exactly the
reason this section gives ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2).
What this section decides is narrower than either: that the *read likelihood* carries neither.)*

### 3.8 What this model gives up at three reads, and what it must not break at three hundred

**At 3 reads a position** the likelihood is nearly flat and the prior decides. Under the model as
specified, one alternative read of three at Phred 30 favours a heterozygote over a reference
homozygote by 5.93 nats — 26 Phred — while the prior at a diversity of 1 in a thousand favours the
homozygote by 6.91 nats, so the call is a reference homozygote by **0.98 nats, about 4 Phred**.
**Move the per-read error rate to 1 in 20 and the same three reads favour the heterozygote by only
2.01 nats, and the call is reference by 4.89 nats — 21 Phred.**

Read that pair the right way round, because the direction surprises people: **a worse error rate
makes the reference call more confident**, since the alternative read is cheaper to explain away. So
at this depth the model's answer turns on the prior and on the *scale* of the error rate, and barely
at all on its per-read resolution — which is the argument for §3.2's calibration being the half that
matters here, and the per-read half being what matters at depth. It is also why a mis-set error rate
at 3 reads a position moves calls silently in whichever direction it is wrong.

**At 300 reads a position** — the GIAB trio's depth is about 313 reads a sample, where a true
heterozygote shows about 156 alternative reads ([`cohort_merge.md`](cohort_merge.md) §4.3) — the
likelihood is overwhelming and the prior is irrelevant. What breaks at this end is not the arithmetic
but the assumption that reads are independent given the genotype: 300 reads multiply any
misspecification 300 times. The two devices the reference callers use against exactly this are
freebayes' read-dependence factor, which discounts the pooled error evidence by about a tenth
([`Parameters.cpp:475`](../../../../freebayes/src/Parameters.cpp)), and samtools' dependency
coefficient, which downweights the *n*-th read supporting the same base and strand geometrically
([`errmod.c:77`](../../../../htslib/errmod.c)). **ng has neither, and neither does production.**
Recorded as deferred (§10) rather than adopted, because a flat
discount is a knob with no reading and nothing here has measured what it should be set to.

---

## 4. The STR emission — Model A

### 4.1 Why this model, and what the comparison actually measured

**Three read likelihood models were implemented behind one interface and scored against synthetic
data with known genotypes.** The winner is the one this document specifies, and the numbers rather
than the verdict are what is worth carrying
([report](../../reports/implementations/ssr_stutter_scoring_model_bakeoff_2026-06-24.md), 16 samples,
3 loci, depth 60, parameters taken from the truth so the models are compared and not the estimators):

| | genotype concordance | calibration error, messy data | time to score |
|---|---|---|---|
| **Model A** — whole-repeat and part-repeat slips priced apart, HipSTR's shape | 1.0000 on all three data sets | **0.0001** | 34 ms |
| **Model B** — stutter summed outside a substitution alignment, the incumbent | 1.0000 on all three data sets | 0.0166 | 8,522 ms |
| **Model C** — one alignment with a cheap whole-repeat gap and a stiff single-base gap | 0.0208 on part-repeat data | 0.9769 | 682 ms |

**Read the table carefully, because the headline is not concordance.** A and B tie at perfect
concordance, so the comparison did not separate them on accuracy at that depth. What separates them
is calibration — how well the model's stated confidence matches how often it is right, where A is
**166 times closer** on the messiest data — and cost, where A is **250 times cheaper** on the same
data. The timings are from a debug build, so the ratio is trustworthy and the absolute numbers are
not. The cost matters because in production's STR pileup the read likelihood is about three quarters
of self-time.

**Model C's failure is a warning this document inherits.** Folding stutter into the alignment — one
score, two gap prices — let part-repeat reads pull the genotype: 2 loci in 100 called correctly on
part-repeat data. That is the reason §4.2 and §4.3 are two separate factors and not one alignment,
and [`alignment.md`](alignment.md) §4.2 records the same result from the aligner's side.

**The rule at a repeat tract is that the mapper's alignment is not trusted, and the read is aligned
again.** What is taken from it is only where the read landed and how much of it was
clipped, which decides whether the read reaches this locus at all. The tract's length is then read off
an alignment this pipeline computed. **So what this model scores is always a sequence somebody
measured, and never a repeat count read off the reference**; the shortcut that would have counted
repeats instead of measuring them was considered and parked.

**Some reads are allowed to skip the alignment**
([`ssr.rs:858`](../../../../src/ng/locus_generation/ssr.rs)). If a read's CIGAR carries only matches —
no indel and no clip — its aligned span brackets the whole window, and its bases over that window are
**byte-identical to the reference**, then it copied the reference across the tract. Its repeat count
is the reference's and its tract sits where the reference's tract sits, so aligning it would spend a
full dynamic-programming pass rediscovering a span already known. **This is an exact shortcut and not
an approximation:** it returns the span the alignment would have returned, so the observation reaching
§4.2 is the same sequence with the same witness either way.

**The exception is narrow, and two cases that look like it are not.** A read carrying so much as one
substitution inside the tract is aligned. So is any read at a locus where the reference's own repeat
run carries on past the tract region typing drew, because there the measurement answers *at least this
long* rather than a length. And the base-quality gate runs on the skipped reads like every other:
skipping the alignment skips *where the tract is*, not *whether the read is good enough to count*.

**It is worth having because most reads at most tracts are that read.** On tomato at three reads a
position it takes 4 reads in every 9 out of the aligner — 57,507 of 132,069 over 24 spans of
chromosome 1. The saving shrinks as a panel diverges from the reference, and the correctness argument
does not depend on its size.

### 4.2 The stutter distribution

**How likely is it that a copy of an allele `L` bases long produced a read showing `L + Δ` bases?**
That is the whole of this section, and it is a genetics model — a statement about how a polymerase
slips — rather than an alignment algorithm. **It belongs to this document.**
[`alignment.md`](alignment.md) §5.2 sets out the same distribution, because a candidate repeat-aware
aligner there would consume it too, and says outright that it is *"not an alignment algorithm"*. §7
below fixes the ownership and asks for that section to be repointed here; the edit is not made here.

**The distribution has two regimes and that split is its defining structure.**

- **A whole-repeat change** — `Δ` is a whole number of repeats. This is slippage: the common event, its size
  measured in repeats.
- **A part-repeat change** — `Δ` is not a whole number of repeats. This is an ordinary small insertion or
  deletion, or an interruption; it is rarer and its size is measured in base pairs.

Each regime splits again by direction, because losing repeats is more common than gaining them — at
tomato dinucleotides by a factor of 4.9, 2,438 reads against 501
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §3).

**Seven numbers say what a read shows.** Five of them split the reads copied from one allele into the
kinds of length they can carry — the allele's own, a whole repeat longer or shorter, part of a repeat
longer or shorter — and those five sum to one. The other two say, within each kind of change, how
much of it is a single step rather than several. HipSTR's field names are in brackets, for anyone
reading the two side by side:

| name | what it is |
|---|---|
| `same_length_share` [`log_equal_`] | the share of reads showing the allele's own length — neither longer nor shorter |
| `whole_repeat_longer_share`, `whole_repeat_shorter_share` [`in_up_`, `in_down_`] | the share of reads a whole repeat longer, and the share a whole repeat shorter — at any size |
| `whole_repeat_one_step_share` [`in_geom_`] | **of the reads that slipped by whole repeats**, the share that moved by exactly one |
| `part_repeat_longer_share`, `part_repeat_shorter_share` [`out_up_`, `out_down_`] | the same two shares, for reads longer or shorter by part of a repeat |
| `part_repeat_one_step_share` [`out_geom_`] | **of the reads that changed by part of a repeat**, the share that moved by exactly one base |

The five shares are one distribution over what a read can show, so the same-length share is
whatever the four change shares leave:
```text
same_length_share = 1 − whole_repeat_longer_share − whole_repeat_shorter_share
                     − part_repeat_longer_share  − part_repeat_shorter_share
```

**The distribution.** Write `Δ` for the read's length change in base pairs and `p` for the repeat
motif's length:

```text
a whole-repeat change  (Δ divisible by p), n = Δ/p repeats:
    n = 0   →  same_length_share
    n > 0   →  whole_repeat_longer_share  · whole_repeat_one_step_share · (1 − whole_repeat_one_step_share)^(n − 1)
    n < 0   →  whole_repeat_shorter_share · whole_repeat_one_step_share · (1 − whole_repeat_one_step_share)^(|n| − 1)

a part-repeat change   (Δ not divisible by p), e = Δ − Δ/p  (truncated division):
    e > 0   →  part_repeat_longer_share  · part_repeat_one_step_share · (1 − part_repeat_one_step_share)^(e − 1)
    e < 0   →  part_repeat_shorter_share · part_repeat_one_step_share · (1 − part_repeat_one_step_share)^(|e| − 1)
```

**Read a one-step share as: of the reads that moved at all, this fraction moved by exactly one step —
one repeat on the whole-repeat branch, one base on the part-repeat branch.** A
larger value concentrates the mass on one-repeat slips; HipSTR ships 0.95, and that is exactly how to
read it — **nineteen slips in twenty are one repeat**.

**It is a share of the slips, not a multiplier applied at each step, and the two are complements.** A
parameter expressed as the chance of *carrying on* to the next step is `1 − this`; production's
`StutterShape` stores that complement and calls it `decay`
([`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs) converts, `geom = 1 − decay`).
**Getting the two the wrong way round inverts the size distribution — large slips become the common
ones — and nothing crashes**, which is why this name says *share* rather than anything that could be
read as a rate of decline. It is also why §12's fourth test exists.

**Two ways to check which one a number is, without reading anyone's source.** The mean slip size is
`1 / this share`, so 0.95 means an average slip of 1.05 repeats and 0.05 would mean an average of
twenty. And HipSTR's own model calls the quantity the *no further step* probability, deriving the
per-step multiplier from it as the complement — the field is called `geom` after the distribution
rather than after the number, which is how this gets read backwards.

**Why part-repeat sizes are re-indexed.** The part-repeat geometric is indexed by `Δ − Δ/p` rather
than by `Δ` so that its support has no gaps: at period 3 the part-repeat values 1, 2, 4, 5, 7 map to
1, 2, 3, 4, 5. Without it the geometric would be evaluated at indices that skip the multiples,
distorting the distribution. It is **not** about double counting — the two regimes are disjoint by
construction, since a change is part-repeat precisely when it is not a multiple of the period.

**Which regime applies is decided by arithmetic and never by what was inserted.** Inserting a `C`
into a poly-A tract is classified as one repeat of slippage exactly as inserting an `A` is. The
composition is caught, but through §4.3's substitution term rather than through the part-repeat
branch. **At period 1 this catches roughly three of every four single-base insertions**, and
[`alignment.md`](alignment.md) §4.2 works through what it costs.

**A slip can land in more than one place.** In a pure repeat, adding a repeat anywhere gives the same
sequence and there is nothing to choose. In an interrupted repeat the placements give genuinely
different sequences, so the model enumerates them and averages with equal weight. **ng enumerates for
whole-repeat slips only and resizes a part-repeat change at the tract's end in a single placement**,
which is production's split and is stated here because otherwise a coder has to guess; it is a
simplification of the same class as the fixed part-repeat share below. **A slip that
cannot be reached from the allele at all — contracting away more repeats than exist — scores zero.**
That leaves the allele's distribution summing to less than one, which under-scores short alleles
relative to long ones. **The size depends entirely on that share and is not always small.** The
shortest tract the copy floors admit is four repeats, at hexamers
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5.1.1), so the unreachable tail there is a
contraction of four repeats or more. At a slippage level of 2 in 100 and HipSTR's shipped one-step share
of 0.95, that tail is about **2 parts in a million**; at a one-step share of 0.5 — a stress value, not
a fitted one — it is about **2 parts in a thousand**, a thousandfold difference. *(An earlier version
called 0.5 "production's fallback decay". It is not: production's only coded 0.5 fallback is
`DEFAULT_G0_FALLBACK_P` ([`param_estimation.rs:167`](../../../../src/ssr/cohort/param_estimation.rs)),
the **genotype prior's** pseudocount decay, a different model. Production codes no fallback for the
stutter one-step share.)* **So it must be computed
and reported per candidate rather than assumed negligible**, and it is one more reason the step
chance is a fitted parameter and not a default. It is the same shape of defect that made Model C
collapse: a model that quietly loses mass on some candidates and not others is comparing them on
different scales.

**Slips beyond a cutoff score zero and the read falls to §4.5's outlier term.** Production applies
one constant, 10, to the *repeat* count on the whole-repeat branch and to the re-indexed *base-pair* count on the part-repeat branch
([`param_estimation.rs:21`](../../../../src/ssr/cohort/param_estimation.rs)) — one number on two
scales. **ng carries two cutoffs, named for what they count**, because 10 repeats at a hexamer is 60 base
pairs and 10 base pairs is not the same claim.

**Their values, because §1.2 says a number with no named source is a defect in this document.**
Neither is fitted and neither has a source in the parameters fit, so both are **inherited from
production at 10 and declared inherited rather than measured** — `max_whole_repeat_slip = 10`
repeats, `max_part_repeat_slip = 10` base pairs. Production's own comment calls its 10 "a provisional
choice", so these are two named constants awaiting a measurement, not findings. **What should set
them is the mass they discard**, which the builder already computes per candidate for §12's fifth
test.

**Where the parameters come from.** All of them from the pre-pass, per read group per stratum
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.1): the level — how often a read slips at
all — supplies the four direction shares once split by the fitted direction ratio, and the fitted
fall-off supplies the two one-step shares. **Two of the seven are placeholders in production and must be
recorded as such rather than mistaken for fitted values**: the part-repeat share is a fixed 5% of
the whole-repeat share ([`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs)), and the
whole-repeat and part-repeat one-step shares are tied to one number where HipSTR keeps them independent.
The follow-up that would replace the first — binning part-repeat reads separately, as HipSTR does —
is recorded in the model comparison's report (§4.1) and in [`alignment.md`](alignment.md) §5.2,
**not** in the parameters fit's own specification;
§10 gives the second a home.

**Licence.** HipSTR is GPL-2 and this project is not, so **the vendored tree is not a source to
implement from** — the rule [`locus_generation_ssr.md`](locus_generation_ssr.md) §3 already carries.
**The source is the manuscript HipSTR itself asks to be cited by** — Willems, Zielinski, Yuan,
Gordon, Gymrek and Erlich, *Genome-wide profiling of heritable and de novo STR variations*, **Nature
Methods 2017** (`nmeth.4267`), named in that project's own README. **Whether its text carries enough
to reimplement the distribution is not something this repository can establish**, so the operative
source for a coder is §4.2 above, which sets the distribution out in full precisely so that nobody
has cause to open that tree.

**Do not treat production's port as a cleared second source until someone establishes that it is
one.** [`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs) implements this distribution
and §9's reuse map ports from it — but its own header says it *"mirrors
`HipSTR/src/stutter_model.cpp::log_stutter_pmf`"*, and a test is named for mirroring the HipSTR
formula. That wording may describe only what the code is equivalent to, and it may describe where the
code came from; nobody now remembers which. **Provenance is what a licence question turns on, so it
has to be settled rather than assumed** — and if it was written from the source, the port is a
rewrite from the paper rather than a copy forward. Either way the header should end up saying which.

*(For contrast, the alignment module's per-base-quality emissions and flank gap costs came from
Dindel — Albers et al. 2011 — and there is no ambiguity there, because **Dindel's source is not
vendored in this repository at all**, so the publication is the only place they could have come
from.)*

### 4.3 The substitution term

**Once the candidate has been stretched to the read's length, the only thing left to explain is
which letters differ.** That comparison is the alignment module's, not this one's
([`alignment.md`](alignment.md) §5.1): two equal-length sequences scored under one flat error rate,
`(1 − ε)` per matching base and `ε/3` per mismatching one, with gaps confined to a couple of bases at
either end.

**Each of the two factors answers one question and is unable to answer the other's, and that is
deliberate.** The stutter factor decides **how long** the read is: it has already stretched or
trimmed the candidate to the read's own length before the comparison starts. The substitution factor
then decides **which letters**, over two sequences that are the same length by construction — so
there is nothing left for it to insert or delete.

**The reason to keep them apart is that otherwise there would be two explanations for one
observation.** Suppose the substitution factor could also delete a few bases in the middle of a
tract. Then a read one repeat shorter than the candidate could be explained twice over: as a
polymerase slip, or as a deletion during sequencing. Both would fit, and **nothing in the data would
say which** — so raising the slippage rate and lowering the sequencing-indel rate would describe the
reads exactly as well as the reverse. The parameters fit measures those as two separate numbers
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.1), and it can only do that because each
one explains something the other cannot.

*(The alignment module enforces the same separation on its own side, for callers that hand it two
sequences of different lengths: a gap may only sit within a couple of bases of either end, which is
slop for a fuzzy extraction boundary rather than an allowance for real indels
([`alignment.md`](alignment.md) §5.1). **This composition never reaches that branch**, since the
stutter factor equalises the lengths first — which is exactly why that document warns the branch will
ship untested if nobody drives it deliberately.)*

**The error rate is the pre-pass's STR substitution rate for this read group and stratum, and it is
deliberately not the SNP path's.** The two are separate parameters that are never tied
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §1.1): each is the error parameter of its own
model, absorbing whatever that model cannot otherwise explain, and forcing one number to carry both
would make each model wrong in a way neither could report. Where a stratum barely slips the two
models describe nearly the same thing and the two rates must agree to within a quarter-Phred, which
is that document's §4.5 and is a test rather than an aspiration.

**This path uses the fitted rate and not the reads' own qualities, and the reason is not that the
qualities are missing — it is that they are the wrong quantity.** The evidence does carry a summed
per-read error on this path too, filled by the STR generator and recorded there as unconsumed
([`locus_generation_ssr.md`](locus_generation_ssr.md) §3). But **what it carries is a per-*read*
number — the chance that this read is wrong somewhere — and the substitution term needs a per-*base*
rate**, applied once for each of the tract's twenty or forty bases. Those are different quantities
and one does not stand in for the other: a read's error probability is roughly the sum of its bases',
so using it per base would overcharge by the tract's length.

**On the SNP/indel path the per-read number is exactly what is wanted**, because an observation there
is one allele and the question is whether the read shows the wrong one. Inside a tract the question
is per base. **That is the whole of the difference between §3.2's answer and this one**, and it is
why the two paths disagree without either being inconsistent.

Two further reasons point the same way. The fitted rate is measured **inside tracts**, where reported
base quality is systematically worse and therefore least trustworthy
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1). And it is fitted per stratum, so it
already conditions on the one covariate that matters most here. §10 keeps the door open for a
per-base treatment if the evidence's shape ever changes.

### 4.4 Which stratum a candidate belongs to

**A read's chance of slipping is a property of the tract it was copied from, and that is the
candidate allele, not the reference.** So the stutter parameters are looked up by the **candidate's**
motif period and repeat count, and they therefore differ between the candidates at one locus. A
candidate of 6 repeats and one of 12 at the same tract are drawn from different strata and slip at
measurably different rates — slippage rises about 1.3-fold per repeat count over the measured range,
and at tomato dinucleotides reaches 15.0% at 12 to 15 repeats
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3, §4).

Production already does this, by a linear function of the candidate's repeat count rather than a table
lookup ([`em.rs:385`](../../../../src/ssr/cohort/em.rs)). ng uses the table, because the table is
what the pre-pass emits.

**The consequence a coder must not miss: the stutter parameters cannot be hoisted out of the
candidate loop.** They are per `(read group, candidate)`, not per locus. What *can* be hoisted is the
lookup, which is a small table indexed by period and repeat count.

**A stratum whose numbers were borrowed rather than fitted is used exactly as a fitted one**, with no
branch and no down-weighting — a borrowed value is the best estimate available and treating it as
second-class would mean inventing a discount. **But the provenance travels**: the locus's output
carries the weakest provenance of any parameter that entered it, so a genotype resting on a direction
split borrowed from two repeat counts away is distinguishable in the run's output from one resting on
a fit. At the bottom of the repeat-count range most direction splits *will* be borrowed and that is
the pre-pass's design working: reaching enough slipped reads to fit one there takes about
310,000 loci for a direction split and about 880,000 for a fall-off, against tomato's 1.73 million in
total
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.5).

### 4.5 The third term on this path: reads nothing explains, and reads from someone else

**Some reads at a repeat tract are explained by no allele at all** — a read from a paralogous tract
elsewhere, a chimera, a somatic length in a long tract, a read whose delimiter anchored wrongly.
Without somewhere to put them, one such read drives every candidate genotype's likelihood to zero.
§2.1's third term is where they go: with weight `λ` a read is drawn from a junk distribution instead
of from the individual's alleles.

**Production sets `λ` to 0.01 and spreads it uniformly over `D`, the number of distinct sequences the
whole cohort showed at that locus** ([`em.rs:393`](../../../../src/ssr/cohort/em.rs)). **ng inherits
the 0.01 and declares it inherited**: it is not fitted, the parameters fit has no source for it, and
§1.2 requires that be said rather than left blank. **`D` is the defect.** A single sample showing two distinct sequences gets a junk floor of 0.005; a 63-accession
panel showing twenty gets 0.0005, ten times lower. So a sample's own genotype likelihood changes when
an unrelated sample is added to the run, which is not a property a per-sample likelihood may have.

**Be precise about how much it moves, because the honest answer is "less than it looks".** For a read
that no candidate explains, the junk term is the same under every genotype and cancels when the
caller normalises. It does not cancel for a read that the candidates explain *weakly* — one whose
emission is near the floor — where the floor sets how much that read's weak preference counts. And it
does not cancel at all in the **data likelihood**, which emission and QUAL read as an absolute number
compared between loci, so a locus called in a large cohort and the same locus called alone are scored
on different scales. **The size at genotyping is unmeasured** and open question 4 (§11) is how to
measure it.

**Decision: the junk distribution is a property of the locus, not of the cohort** (owner,
2026-08-19). It is **uniform over the tract lengths the stutter model's support reaches from any
candidate at this locus** — a count this model computes from the candidate set alone, without asking
what any sample showed. That keeps the read likelihood a per-sample function at every cohort size,
which goal 4 requires and which nothing else in this document breaks. *An earlier version left this
as a leaning while §6 and §9 already described it as decided; the inconsistency was the reviewer's
finding, and this is the side it is resolved on.* **Open question 4 (§11) survives as the size
measurement only** — how much production's cohort-wide denominator was actually moving.

#### 4.5.1 Contamination at a repeat tract — built, on by default, and measured

**A contaminating read at a tract is not junk: it shows a length that is a real allele in some
individual, and often a common one.** Today it has to be explained as slippage where the stutter
model can reach it — inflating the apparent slip rate — or fall to the outlier floor where it cannot.
The first is how a contaminated sample gets called heterozygous for its contaminant's allele, and it
is the same failure the SNP path's mixture exists to prevent (§3.6).

**No caller models it — not ng, not production, not HipSTR, not GangSTR** (checked: `contamin`
appears nowhere in GangSTR's source). **It is nonetheless not a new term**, only the third term with
its two ways of not coming from this sample kept apart:

```text
log Lg(g)  =  Σ_o  n_o · log[ (1 − λ − c) · Σ_a (k_a/P) · Lr(o | a)
                            +      λ      · uniform over the reachable lengths
                            +      c      · seed(o) ]

λ    the outlier weight — reads no allele explains (§4.5), unchanged
c    this read group's contamination fraction, from the pre-pass
seed the prior's geometric decay away from the cohort's modal repeat count at this
     locus, normalised to a distribution over lengths
```

**Three components rather than two, and folding them together would be wrong.** A junk read can show
anything, so its distribution is flat and its term cancels between genotypes (§4.5); a contaminating
read shows a *plausible* length, so its distribution is peaked and its term does not cancel. Adding
`c` to `λ` and keeping one flat distribution would give the contaminant's reads the junk treatment,
which is what the model already does and what this section exists to stop. `1 − λ − c` is floored
positive; at the values in play — `λ` of 1 in 100 and `c` of 1 to 3 in 100 — it is nowhere near the
floor.

**Why the fraction may be carried over from ordinary sites.** How contaminated a sample is is a
property of its DNA and not of the locus, so the per-read-group estimate the pre-pass fits from the
census sites applies at a tract though no tract is in that census
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2, §5). **One caveat, and it wants measuring rather than
assuming.** The fraction is a share of *reads*, and at a repeat tract a read has to earn its place
before it counts: it must carry both ends of the tract with a few bases of unique sequence on each
side. A read that runs out inside the tract says only *at least this long* (§5), and one that falls
entirely inside says nothing and is thrown away. **How often a read manages that depends on how long
the tract is** — so if the contaminating individual's allele is longer than this sample's, its reads
qualify less often and the contaminant is under-represented among the reads the model scores; if it
is shorter, over-represented. **The direction is clear in each case; what nobody has measured is
which case is commoner across a genome, or by how much.**

**Where the second half of the mixture comes from, and why it is harder here.** The mixture needs two
things: how *often* a read comes from the contaminant, which is the fraction above, and **what such a
read shows**. At an ordinary site the second one is three numbers — how often the contaminating
population carries the reference, a substitution, an insertion or deletion — and the parameters fit
measures them. **At a tract a contaminating read shows a tract length, so what is needed is a list of
how common each length is at this locus**, and no two tracts have the same list.

**Two requirements, and the obvious source meets only one.** The list has to be **specific to this
locus**, since the lengths in play differ at every tract; and it has to be **fixed before the
caller's loop starts**, because contamination is frozen and must not move from one pass to the next
(§3.6). The cohort's own fitted allele frequencies at the locus are the natural answer to the first
— they are precisely how common each length is in this panel — and they fail the second, because
they are what the caller rewrites on every pass.

**So the list comes from the genotype prior's starting shape**
([`calling_priors.md`](calling_priors.md) §5.1): mass falling away geometrically from the commonest
length the cohort actually showed at that tract, at a rate the parameters fit measured over a group
of loci. It is specific to the locus, because the length it is centred on is that locus's own; and it
is computed once from the observations rather than from the caller's estimates, so it does not move
while the caller iterates. **Both requirements met, and neither by accident.**

**What is weaker here than at ordinary sites, and how much weaker.** At a SNP both halves of the
mixture are measured: the fraction, and the frequency of the allele the observation shows. At a tract
only the fraction is, and the seed stands in for the second half — so the term assumes the contaminant was
drawn from the same population as the cohort. **That assumption is bounded by `c` and it is not the
one `verifyBamID2` was built to remove.** An earlier version of this section cited that tool's
recorded failure — a true 10% contamination returned as 2.9% when the allele frequencies came from
the wrong population
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.2) — as though it applied
here. It does not: that failure is in *estimating* the fraction, which happens at ordinary sites by
the proper method and is already done by the time this term runs. What a wrong seed can spoil is only
the *application*, and its whole weight is `c`.

**So work through what a wrong seed actually does.** If the contaminant's lengths are unlike the
cohort's mode, its reads simply are not covered by the seed and fall to the outlier floor instead —
which is where they go today, so nothing is lost. If the seed puts mass on the mode and the sample's
own modal reads are partly attributed to the contaminant, at most `c` of their weight moves, and it
moves equally for every genotype that carries the mode. **The realistic failure is therefore small
and bounded, which is why this is on rather than off.**

**It costs one number inside a logarithm that is already being taken**, and it costs this path
nothing in exactness — the SNP mixture forces per-read qualities to collapse into one number per
observation (§1.4), and this path already uses a flat fitted rate per stratum, so it gives up nothing
it was using.

**Decision: on by default wherever the pre-pass emits a fraction above its floor, the same rule as
§3.6, and kept on unless a measurement shows it harming calls** (owner, 2026-08-19). **The reason is
that contamination is a property of the sample and not of the marker.** If a library carries another
individual's DNA, that DNA is at the repeat tracts as much as at the SNPs, and a caller that corrects
for it at one and not the other is treating one number as two. The weaker half of the mixture here is
the seed, not the fraction, and the paragraph above bounds what a wrong seed can cost.

**Built behind the model seam all the same**, so the measurement can be run against a build with it
off. §12's item 17 is that measurement and open question 7 (§11) says what would switch it off.

### 4.6 Two alleles of the same length

**ng stores an STR allele as the sequence the reads showed, not as a repeat count**, because two
alleles of the same length can differ by an interior substitution — an interrupted repeat — and a
count cannot tell them apart ([`locus_generation_ssr.md`](locus_generation_ssr.md) §3).

**The emission handles this correctly and needs no decision.** §4.2's stutter factor depends only on
the length change, so both same-length alleles get the same one; §4.3's substitution factor compares
the read's letters against each candidate's own letters, so a read carrying the interruption scores
higher against the interrupted allele and lower against the pure one, by `log(3(1−ε)/ε)` per
distinguishing base — at an error rate of 1 in 200, 6.4 nats or 28 Phred. **So the read likelihood
separates them and the genotype prior does not**: the prior's geometric seed is indexed by repeat
count, so both alleles land on the same rung and the sibling document leaves how to divide that rung
open as its own third question ([`calling_priors.md`](calling_priors.md) §5.2). Nothing here needs to
wait for that answer.

**One thing does need care.** The placement enumeration of §4.2 is what makes an interrupted
candidate's slip land in several distinct sequences, and the interruption also changes how much the
tract slips at all — production carries a purity factor on the level for exactly that
([`em.rs:387`](../../../../src/ssr/cohort/em.rs)). The pre-pass fits its strata over loci of mixed
purity and **emits no purity of its own** — `SsrSegment::purity_fraction()` is region typing's
measure, which that document uses to say when its rate comparison is valid, and §4.4's list of what
the walk emits does not contain it
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.5), so **the level a candidate gets is the
stratum's, unadjusted for that candidate's own purity**, and the adjustment production applies has no
fitted source here. Recorded as deferred (§10) rather than invented.

### 4.7 What this model gives up at three reads, and at three hundred

**At 3 reads a locus** — tomato's depth — a tract of six repeats or more slips at 2.0%
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5 — **that rate is HG002's, pooled over
periods; tomato's own has not been quoted here**), so **about one locus in seventeen shows a single
slipped read and the other sixteen show none.** The stutter parameters therefore barely move
a genotype at this depth: what decides it is the substitution term, the candidate set, and the prior.
That is worth knowing before anyone attributes a tomato result to the stutter model.

**At 300 reads a locus** the same tract shows about 6 slipped reads, and the stutter parameters are
what separates a genuine heterozygote whose two alleles differ by two repeats from a homozygote read
with two repeats of stutter. **This is where the model earns its keep, and it is also where the two
placeholders of §4.2 bite hardest**: at 300 reads a locus carries about 0.3 part-repeat reads at
production's fixed 5% share of a 2% level, so the part-repeat parameters remain thinly exercised
even at depth, while the whole-repeat ones are measured 6 reads at a time per locus and thousands of loci
at a time in the fit.

**Below the copy floors the model describes the data badly and says so.** Tracts under four repeats
are 19.9% of loci, produce 1.7% of all slippage, and **58.5% of even that is not a whole-repeat
change** — an ordinary small indel, which this model has nowhere to put
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5). Those tracts are excluded from the STR
path upstream by the copy floors, and this section exists so that nobody lowers a floor without
knowing that the model on the other side of it is mostly describing the wrong thing.

### 4.8 A property this model should have and has never been tested for

**Where a tract does not slip, sending it down the STR path should cost time and not accuracy**
(owner, 2026-08-19). The two paths would then differ only in what they cost to run, and the choice of
which one a locus goes down would stop being a choice about how well it gets called. **Nothing has
tested that**, and the property is worth writing down because the things that would break it are
knowable in advance.

**One of the two halves is already required elsewhere.** Where a stratum barely slips, this path's
noise model *is* the SNP/indel path's, so their two substitution rates have to meet to within a
quarter-Phred or one of them is wrong — that is the parameters fit's own requirement, not an
aspiration ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.5). If the rates agree and the
genotypes still do not, the difference is in the models rather than in the numbers.

**Four things could put the two paths apart at zero slippage, and the first is structural.**

- **One path charges error per base and the other per read.** The STR path scores the whole tract
  base by base, so a read differing at two bases of a 20-base tract pays `(ε/3)²`. The SNP/indel path
  charges a read one error probability however many bases it differs at, because its per-read number
  is the worst base quality across the locus (§3.1). At an error rate of 1 in 200 those are
  **1,970 times apart at two mismatches — 33 Phred — and 61 Phred at three**. **Neither is obviously
  the right one:** per base is the more faithful model, per read is the approximation that survives
  the merge's aggregation (§2.3). At one mismatch the two agree to within a Phred, which is why this
  has never shown up.
- **The outlier term has no counterpart on the SNP/indel path** (§2.2). It reweights every read by
  about 1 in 100, which mostly cancels between genotypes and does not exactly.
- **The candidate sets are built differently** — a rung ladder against the observed alleles — so the
  two paths need not even be choosing between the same genotypes.
- **The priors differ deliberately** and are not this document's to reconcile: the STR prior is
  centred on the cohort's modal length and the SNP/indel prior privileges the reference, for reasons
  [`calling_priors.md`](calling_priors.md) §5 argues at length. **So the property has to be stated
  about accuracy and not about agreement** — the two paths should call *correctly* at the same rate,
  not identically.

**Open question 8 (§11) carries it, and §12's item 18 is the measurement.**

---

## 5. Reads that saw only part of the locus

### 5.1 What a partial observation claims, and why it has been thrown away until now

A read that entered the locus and ran off its own end shows a prefix or a suffix of what the sample
carries. It does not say what the sample carries; it says the sample carries **at least** this.
Production discards these outright at repeat tracts, counting them and moving on
([`locus_tally.rs:91`](../../../../src/ssr/pileup/locus_tally.rs)). **HipSTR does not** — its own
documentation says a read that does not fully extend across the repeat still gives a lower bound and
that HipSTR uses it. So the censored term below has a precedent in a working caller, and what is new
here is carrying it on the generic path as well. ng's locus generators keep them,
and [`locus_generation_ssr.md`](locus_generation_ssr.md) §3 calls that *"the one place this path is
new rather than a port"* — and then hands the modelling here: *"whether and how a likelihood uses a
lower bound — the censored term — is the caller's to design"*.

**This document designs it and turns it on.** Two reasons. The evidence is not small: on chromosome 1
of one tomato accession the generator recorded 7,085 such observations (`SRR7279503`; what share of
all observations that is has not been measured, and whoever builds this should measure it first).
And it is the only route to alleles longer than a read, which is otherwise a whole class of variation
this caller cannot see.

**The danger is precise and must be designed against: a partial scored as though it were complete
mis-scores as a *short* allele**, because its bases are a prefix of the truth. That is why
`complete_observations()` exists as a deliberate guard on the evidence type
([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs)) and why the sections below
are careful that a lower bound produces a *less discriminating* likelihood, never a differently
biased one.

### 5.2 Scoring a read that ran out inside a tract

**This section says how to put a number on *at least this long*.** §4.2's stutter distribution
answers a question a partial read cannot ask: *how probable is it that this allele produced a tract
of exactly this length?* A read that ran out gives no exact length. What it gives is a floor — the
tract is at least as long as the stretch the read got through — so the question becomes **how
probable is it that this allele produced a tract of at least that length**, and the answer is the
stutter distribution added up over every length at or above the floor.

**That sum is closed-form and finite**, because §4.2's two cutoffs bound the distribution's support,
so nothing is approximated by taking it. The only approximation in this section is a separate one and
the last part of this section says what it costs.

A partial observation shows `ℓ` base pairs of tract and stops. Under §4.2, for a candidate `a` of
length `L`:

```text
Lr(partial of length ℓ | a)  =  Σ            P_stutter(Δ | a) · subst( the ℓ bases | first ℓ bases of a stretched by Δ )
                                Δ : L + Δ ≥ ℓ
```

**In words: sum the stutter distribution over every length change that would still leave the tract at
least as long as what the read saw, weighting each by how well the letters matched.**

**The factorised form, which is what to implement.** For a pure repeat every reachable stretching
agrees on the first `ℓ` bases, because stretching appends or trims at the tract's end, so the
substitution factor comes out of the sum:

```text
Lr(partial of length ℓ | a)  ≈  subst( the ℓ bases | first ℓ bases of a's own tiling )
                                 ×  P( tract length ≥ ℓ | a )
```

**and the tail probability is a closed form**, because a geometric's tail is a geometric. Writing
`n₀` for the smallest whole-repeat count satisfying `L + n₀·p ≥ ℓ` and `e₀` for the smallest
part-repeat index satisfying the same:

```text
P( length ≥ ℓ | a )  =  Σ over the whole-repeat terms from n₀ upward, capped at the repeat cutoff
                      + Σ over the part-repeat terms from e₀ upward, capped at the base-pair cutoff
```

Both sums are finite, because §4.2's cutoffs truncate the support, so **the tail is an exact finite
sum and not an approximation at all.** Only the factorisation is approximate.

**The size of the factorisation, and when to pay for the exact form.** It is exact for a pure tract.
For an interrupted candidate the placements give different first-`ℓ` bases, and the error is bounded
by the spread of the substitution factor across placements — at most `log(3(1−ε)/ε)` per
distinguishing base inside the witnessed prefix, the same 6.4 nats §4.6 quotes at an error rate of 1
in 200. **So: use the factorised form on pure candidates and the exact sum on interrupted ones.** The
exact sum costs one substitution evaluation per reachable length change, which is at most the repeat
cutoff plus the base-pair cutoff, and interrupted candidates are the minority.

**Two properties this form has and a coder should check they still hold after any change.** A partial
is always *less* discriminating than a complete observation of the same bases, because a tail
probability varies less between candidates than a point probability does. And a partial that could
have come from every candidate contributes the same to every genotype and drops out, which is what
should happen to a read that saw too little to say anything.

### 5.3 Scoring a read that saw only part of an ordinary locus

A partial observation here covers a run of positions inside the locus span and nothing outside it.
Its bases cannot be compared against a whole-span allele — that comparison would report a read
agreeing with the reference over everything it saw as non-reference, which is the trap the evidence
type's own documentation names.

**The rule: an allele is compatible with a partial observation when the allele's projection,
restricted to the positions the read witnessed, equals the read's bases.** Then, for genotype `g`:

```text
term(o | g)  =  Σ  ( k_a / P )   over the alleles a in g that are compatible with o
```

and if no allele in `g` is compatible, the observation is charged as an error exactly as in §3.3,
with `m = 1`, because a multi-position difference has no finite set of wrong outcomes to divide by.

**This is exactly aggregable**, which is why it is the rule chosen: the witnessed run is already part
of an observation's identity ([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs)
— `read_witness` carries the offset and the number of positions covered), so every read pooled into
one observation witnessed the same stretch and gets the same compatibility verdict.

**What it gives up.** A partial that is compatible with two of the genotype's alleles contributes
`(k_a + k_b)/P`, which is 1 for a diploid heterozygote — no information, correctly. A read that
witnessed one base of a ten-base deletion says almost nothing and is scored as saying almost nothing.
**The model gains nothing at ordinary sites, where nearly every read spans the single position, and
gains at exactly the loci where reads run out: long deletions and wide overlapping-variant groups.**

### 5.4 Where partial reads exist at all, and what that means for the merge

**Two separate things get confused here**, so they are separated first. **Scoring partial reads at a
locus the merge already built** is what §5.2 and §5.3 specify. **Letting partial reads rescue a locus
the merge would otherwise discard** is a change to the merge's own rule, and §5.4.2 concludes it
should not be made on the generic path.

**Correction, 2026-08-21: an earlier version of this section said the first of those "needs nothing
from any other component: the evidence is there". It is not there.** The merge's collation drops
every observation whose witness is not `Complete` before projecting it onto an allele
([`cohort_merge/build.rs:1351`](../../../../src/ng/run/cohort_merge/build.rs)), and its projection
panics rather than pad one ([`:323`](../../../../src/ng/run/cohort_merge/build.rs)). So a cohort
observation as built carries no partial reads at all, and §5.2 and §5.3 have nothing to read.
**This is a second requirement on the merge, alongside the read-group one of §2.3, and §7's table
carries it:** a partial observation must survive collation, keyed by the stretch it witnessed, and
projected over that stretch rather than the whole locus span. It changes no locus's existence — that
is §5.4.2's separate question, still answered *no* — only what evidence a built locus hands on.
**Where it bites hardest is the repeat path**, where §5.4.1's bottom row shows over half the reads at
a 60-base tract are partial and an allele longer than a read can only ever be witnessed partially.

#### 5.4.1 A read is partial only where the locus is wide on the reference

**The witness is measured in the locus's reference positions**
([`witness.rs`](../../../../src/ng/locus_generation/witness.rs)), and that single fact decides which
variant classes have partial reads at all. For reads of length `L` over a locus covering `W`
reference bases, with read starts uniform, the share of overlapping reads that are partial is
`(2W − 2)/(L + W − 1)`. At 150-base reads:

| locus reference width | share of overlapping reads that are partial |
|---:|---:|
| 1 | **none** |
| 10 | 11 in 100 |
| 30 | 32 in 100 |
| 50 — the widest generic locus ng builds | 49 in 100 |
| 100 — the widest repeat tract the segmentation admits | 80 in 100 |

**So on the SNP/indel path partial reads barely exist, and where they exist they are mostly
reference reads.** Three consequences, and each follows from the table rather than from judgement:

- **A single-base substitution has none.** Its locus is one reference base wide.
- **An insertion has none either, and that is not obvious.** An insertion's reference span is its
  anchor base, however long the inserted sequence
  ([`cohort_merge.md`](cohort_merge.md) §3.1) — width stays 1, so every read covering the anchor is
  complete. *(A read that ran out inside a long inserted sequence is a different problem — a
  truncated allele sequence, not a partial witness — and it belongs to read preparation, not here.)*
- **A deletion has them, and they carry the reference allele preferentially.** A read carrying the
  deletion crosses every deleted reference position without spending a single read base, so it
  reaches the far side of a wide locus far more cheaply than a read carrying the reference does.
  **The alternative evidence at a deletion is therefore biased towards complete and the partial reads
  towards reference.**

**The repeat path is where the table's bottom row lives, and it is the opposite case.** A tract's
locus is as wide as the tract — up to 100 bases — so over half of the overlapping reads are partial
at a 60-base tract, and **an allele longer than a read can only ever be witnessed partially.** That
is not a corner: it is the whole class of expanded repeats.

#### 5.4.2 So the merge's rule stays as it is on the generic path

The merge decides whether to build a locus at all by asking whether some single sample showed enough
reads disagreeing with the reference, and it excludes partial reads from both halves of that question
— because comparing a partial's shorter bases against the whole locus's reference bases would call a
faithful partial non-reference
([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs), `non_reference_reads`).

**An earlier version of this section recommended changing that, and §5.4.1 says why it was wrong.**
The classes where the rule could lose a real variant for want of partial evidence are the classes
that have no partial reads (substitutions, insertions) or whose partial reads carry the reference
(deletions). **The rule is a cost bet and it is a good one** — on one 100 kb tomato interval,
assembling the surviving loci cost 170 ms of a 425 ms single-threaded merge
([`cohort_merge.md`](cohort_merge.md) §4.3) — and nothing in §5.4.1 argues for paying more of it.

**The filter applies to repeat tracts too** (owner, 2026-08-19): a catalogue tract no sample varied
at is discarded like any other locus, and being in the catalogue buys it nothing. **One line of the
rule has to change for that to mean what it says, and only on this path.**

A read that runs off its own end inside a tract stays `Partial` even when its witnessed run covers
every position of the locus — the saturation is a placement convention and not a claim to have seen
the whole tract, which
[`witness.rs`](../../../../src/ng/locus_generation/witness.rs) states outright: *"a reach that
saturates says the read ran out"*. **So a sample carrying an allele too long for a read to span shows
no complete observation at all**, and the filter, which counts complete observations only, reads that
as *nothing varied here*. **Not varying and carrying an allele too long to measure would be the same
answer**, and the second is precisely the class §5's censored term exists to reach.

**The repair is the restricted comparison §5.3 already needs:** count a partial as non-reference when
its bases differ from the reference **over the stretch it witnessed**, and count it in the
denominator on the same terms, so the rule's two halves keep counting the same reads. The witness
carries the offset and the length, so nothing new has to be computed. It changes no verdict at a
tract any sample spans, because a spanning sample reaches the bar on its complete reads alone; it
changes the verdict only where every read carrying the allele ran out, which is the one case at
issue.

---

## 6. What varies per locus, what does not, and what happens at one sample and at a thousand

**As specified here, nothing in this model varies with the cohort, and that is a property to defend
rather than an accident.** The read likelihood is a per-sample question — *given this sample's
genotype, how probable are this sample's reads* — and every parameter it reads is fitted per read
group, per sample, or per stratum. Adding a sample to a run must not change any other sample's
likelihood. **§6.1's second tier is the one thing that would break it**, and that is why it is set
out rather than left to be discovered.

| what | grain | changes with the locus? | changes with the cohort? |
|---|---|---|---|
| per-read-group error rate and its calibration scale | read group | no | no |
| slippage level, direction split, fall-off, STR substitution rate | read group × stratum | **yes** — through the candidate's stratum (§4.4) | no — **as this document specifies them**, and the first three are a live candidate for per-locus re-estimation, which would change that answer (§6.1) |
| contamination fraction and the contaminating population's frequencies | read group | no | fitted from the cohort, then **frozen**; absent at one sample |
| the allele table and the candidate set | locus | yes | yes — the merge unifies across samples |
| the outlier term's spread | — | should be per locus | **is per cohort today, and §4.5 is why that is wrong** |

### 6.1 What must stay frozen, what the EM loop may re-estimate, and what it re-estimates elsewhere

**Three tiers, and only the middle one is open.**

**Tier one — frozen for the whole run, and this document requires it.** The per-read-group error rate
and its calibration scale, the STR substitution rate, the contamination fraction and the contaminating
population's frequencies. **Contamination is frozen by a ruling and for its own reason** — how
contaminated a sample is is a property of that sample and not of any locus, so a per-locus loop has
nothing about it to re-estimate (§3.6, owner, 2026-08-19). **For the rest the reason is not tidiness,
it is selection.** The caller only ever sees
the loci that survived the merge's variability filter — positions selected precisely for carrying
non-reference reads. An error rate re-fitted from those loci would be an error rate measured on the
sites least representative of the genome, and it would come back enormous. That is the same class of
defect the pre-pass exists to remove, arriving from the other end
([`parameter_prepass.md`](parameter_prepass.md) §2). **A caller that re-estimates its own noise rate
from the sites it kept is measuring its own selection.**

**Tier two — frozen per run as this document specifies them, and a genuine candidate for per-locus
re-estimation.** The stutter shape — the direction split and the fall-off — and the slippage level.
**The case for adapting them is real, and both reference implementations do it — HipSTR far more
thoroughly than this document had it.** Its expectation-maximization re-estimates **all six stutter
parameters per locus**, jointly with the genotypes: both direction shares and the one-step share, in
each of the whole-repeat and part-repeat branches, and therefore the same-length share as well since
that is what the four direction shares leave. It starts from a fixed guess and iterates until nothing
moves ([`em_stutter_genotyper.cpp`](../../../../HipSTR/src/em_stutter_genotyper.cpp), the model
rebuilt at each pass; *an earlier version of this paragraph said three parameters, which was wrong*).
**Production adapts less: the shape and the level, in up to three rounds**, each shrunk toward the
parameters fit's per-stratum value by 50 pseudo-counts for the shape and 20 slipped reads for the
level, with the rounds settable to zero to keep the seed
([`em.rs:69`](../../../../src/ssr/cohort/em.rs)). **ng as specified here adapts none of them.** So the
three implementations sit at three points on one axis, and this document picks the frozen end
knowingly rather than by omission. The case for it is that a tract can behave unlike
its stratum for reasons the stratum cannot see — an interruption, a nearby indel, somatic instability
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.6).

**Decision: ship frozen, build the seam, and settle it on a measurement rather than on this
paragraph** (owner, 2026-08-19). *An earlier version argued the case from three properties. Only one
of them was a reason; the other two were hypotheses dressed as objections, and they are recast below
as things to measure.*

**The one that is a constraint rather than a hypothesis.** The stutter parameters arrive through the
per-call context and nothing in §4.2 asks where they came from, so the EM loop can turn per-locus
adaptation on without touching this model. **That has to stay true** — it is what lets the two be
compared in one build — and it is the whole of what this document requires.

**What was wrongly offered as an objection: that a per-locus refit makes the read likelihood
cohort-dependent.** It does, and **the model is already cohort-dependent in the same way**: the
per-stratum numbers are fitted from the whole cohort's walk. A refit at the locus differs in *how
much evidence stands behind each number*, not in kind, and depending on the evidence is what an
estimator is for. **The genuine concern is narrower**: parameters re-estimated from the very reads
being genotyped can fit that locus's noise rather than its chemistry, and how badly is a number
nobody has.

**The axis along which it would fail is reads at the locus, not cohort size.** At 300 reads on one
sample there is plenty at a tract to fit from; at 3 reads across 63 accessions there is not, and the
shrinkage exists to refuse the attempt. So a single deep sample may well recover a locus's own
parameters where a wide shallow panel cannot — which is the opposite of how the cohort-dependence
framing made it sound.

**The hypothesis to test, stated so it can be wrong** (owner): most loci will land close to their
stratum's fitted values, and a minority will sit a long way off — and it is that minority the
adaptation exists for. **The measurement is §12's item 19**, on the two benchmarks the STR path
already has: HG002's tandem-repeat benchmark, where the truth is an assembly rather than another
caller, and tomato's recurrence-based standard
([`silver_standard.py`](../../../../benchmarks/ssr_tomato1/scripts/silver_standard.py)), where the
truth is weaker but the depth is 3 reads a position and the panel is 63 accessions — between them
they cover both ends of the axis that matters.

*Leaning, and it is a leaning about order rather than about the answer: measure adaptation against
§4.2's two placeholders first* — a part-repeat share fixed at 5% and two one-step shares tied to one
value are known misspecifications of this same emission, and un-nailing them is already scheduled
work, so a per-locus refit measured before them would be fitting around a defect.

**Tier three — re-estimated every iteration, and this model never sees it.** The per-locus allele
frequencies. They enter the **prior** as expected allele copies summed over the other samples
([`calling_priors.md`](calling_priors.md) §6), and no term of §2.1 reads them. That separation is
worth keeping deliberately: it is what lets the read likelihood be computed once per
`(sample, observation, candidate)` and reused across every iteration of the caller's loop, which is
the difference between a cheap EM and an expensive one.

**One sample.** Everything above except contamination is available and unchanged. Contamination is a
comparison between samples and does not exist at one, so `c` is absent and §3.3's formula runs — the
non-contaminated path, not a contaminated path with a fitted zero. **The single-sample case is
therefore the *simple* case for this model, not the weak one.** Where it is weak is upstream and
downstream: the candidate set is thinner and the prior carries more of the call
([`calling_priors.md`](calling_priors.md) §6).

**A thousand samples.** Nothing here grows with the cohort. What grows is the number of alleles the
merge unifies at a locus, which is capped, and the number of genotypes, which grows as
`C(A + P − 1, P)` — 21 at six alleles and a diploid, 126 at a tetraploid. The per-sample cost is
`observations × genotypes` multiply-adds and nothing else, so the model is linear in cohort size with
no shared state, which is what lets the EM document parallelise it however it likes.

**Three reads a position.** §3.8 and §4.7 give each path's answer. In one sentence: the likelihood is
nearly flat, the prior decides, and what matters about this model there is the *scale* of its error
rates rather than their per-read resolution.

**Three hundred reads a position.** The likelihood is overwhelming and the prior is irrelevant. What
matters is that the model's misspecifications are multiplied 300 times: the positions §3.7 hands to
the site filter, and §4.2's two placeholder parameters, are the named ones, and §3.8 records that ng has no general
overdispersion device where freebayes and samtools both have one.

---

## 7. Who owns what

**This boundary is already recorded in two other documents from their own side, and the three must
say the same thing.**

| the thing | owner | this model's relation to it |
|---|---|---|
| measuring a read's repeat tract — the ruler | [`alignment.md`](alignment.md) §4.2 | consumes the measurement; never re-measures |
| comparing two equal-length sequences under one flat error rate | [`alignment.md`](alignment.md) §5.1 | **composes** it as §4.3's factor |
| the stutter distribution — how likely a length change is | **this document**, §4.2 | [`alignment.md`](alignment.md) §5.2 sets it out because a candidate aligner there would consume it, and states it is not an alignment algorithm. **That section should be repointed here**, the way [`cohort_merge.md`](cohort_merge.md) repointed [`run_streaming.md`](run_streaming.md) §10, **and its *in frame* / *out of frame* wording moved to §1.3's**; neither edit is made here |
| every parameter | [`parameter_prepass.md`](parameter_prepass.md) and its two path siblings | reads them frozen; **fits nothing** |
| the evidence — observations, counts, summed moments | [`cohort_merge.md`](cohort_merge.md) | consumes it, and placed **two** requirements on it: **keep read group in the identity** (§2.3) — **built 2026-08-23**, rows are one per `(allele, read group)` — and **let a partial observation survive collation**, keyed and projected over the stretch it witnessed, because the built merge discards it and §5 has nothing to score without it (§5.4, corrected 2026-08-21), which is still owed. Neither changes which loci get built — an earlier draft asked for that third change and §5.4.2 withdrew it |
| which alleles are candidates | candidate generation, [`ng_proposal.md`](ng_proposal.md) step 6 | scores what it is handed |
| the genotype prior | [`calling_priors.md`](calling_priors.md) | added to this in log space by the third document |
| when this is called, with which candidates, and how the results are combined | the EM loop document | defines the function only |
| **which parameters are re-estimated between the caller's iterations** | the EM loop document | this document sorts them into three tiers (§6.1) and specifies the frozen behaviour as the default; it does not decide the middle tier |

**One interface correction, and it matters before any code is written.** The step sketch in
[`ng_step_interfaces.md`](../arch/ng_step_interfaces.md) §7 proposes
`read_log_lik(read, allele, params)` — one read at a time. **There are no reads at this point in the
pipeline** (§1.4). The interface is one distinct observation with a count, and the outer call returns
a whole genotype row, because §2.1's junk term puts a logarithm around a sum over alleles and a
per-allele function cannot express it.

---

## 8. Cross-cutting concerns

**Numerics.** Everything in log space at the genotype level; the inner mixture of §2.1 is evaluated
in probability space and logged once, which is where the single `log` per observation per genotype
goes. Every probability that could reach zero before a logarithm is floored, so an impossible
observation yields a finite, very negative number and one read cannot turn a sample's whole row into
`NaN`. Production floors the per-allele error rate into `[1e-12, 0.5]`
([`contamination_estimation.rs:1449`](../../../../src/var_calling/contamination_estimation.rs)) and
the stutter geometrics into `(0.01, 0.99)`
([`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs)); ng carries both, as named
constants with their reasons rather than as bare numbers.

**Determinism.** The model is a pure function of (frozen parameters, observations, candidate
alleles). No RNG, no clock, no thread-dependent iteration, no accumulation whose result depends on
order — the sum over observations must run in a fixed order, which the merge already guarantees by
sorting observations by bytes. This is the same requirement byte-identical output at any worker count
already imposes ([`run_streaming.md`](run_streaming.md) §12).

**Cost and memory.** Per sample per locus, `observations × genotypes` inner terms on the SNP/indel
path. On the STR path the emission is evaluated once per `(observation, candidate)` and reused across
every genotype containing that candidate — **the caching is not an optimisation, it is what makes the
cost `observations × candidates` instead of `observations × genotypes`**, a factor of 10 at six
candidates and a diploid. **Nothing may allocate inside the per-sample loop.** The caller hands in
scratch sized by candidate count and observation count and the model fills it; production lifted
exactly these buffers out of its own iteration after a profile put the allocator's self-time at about
16% of cycles ([`posterior_engine.rs:1874`](../../../../src/var_calling/posterior_engine.rs)). The
STR path additionally needs a placement buffer and a substitution-alignment buffer, both per worker.

**Errors.** The model has no failure mode of its own that is not a caller bug: a genotype table that
disagrees with the allele count, a negative count, a parameter outside its declared range, a
candidate whose stratum is not in the table. These are assertions, and the structural ones must hold
in release — production asserts the analogous invariant because a scratch array too short for the allele count
would otherwise be indexed out of bounds silently
([`per_group_merger.rs:1963`](../../../../src/var_calling/per_group_merger.rs)).

**Provenance.** The model does not branch on it (§4.4). It **propagates** it: the weakest provenance
of any parameter that entered a locus travels onto that locus's output, so a call resting on borrowed
or defaulted parameters is distinguishable from one resting on fitted ones without re-running
anything.

---

## 9. Reuse map

| what | production code | how ng reuses it |
|---|---|---|
| the SNP/indel closed form | [`per_group_merger.rs:1948`](../../../../src/var_calling/per_group_merger.rs) | **shape ported, two terms changed** — the multinomial coefficient dropped (§3.4), the error mass divided by three where the alleles differ by one substitution (§3.5) |
| the per-read error probability | [`open_record.rs:793`](../../../../src/pileup/walker/open_record.rs), [`:944`](../../../../src/pileup/walker/open_record.rs) | ported — the worse of the lowest base quality across the window and the mapping quality — and then rescaled per read group (§3.2) |
| the read-group calibration scale | — | **new**; production has nothing to port because it never fitted an error rate |
| the contamination mixture | [`posterior_engine.rs:1475`](../../../../src/var_calling/posterior_engine.rs) | ported as §3.6, with its geometric-mean substitution recorded as a known cost rather than inherited silently |
| the STR stutter distribution | [`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs) | **ported**, with the two placeholder parameters named as placeholders, the slip cutoff split into two constants — one counting repeats, one counting bases — and *in frame* / *out of frame* renamed (§1.3, §4.2) |
| the substitution comparison | [`pair_hmm.rs`](../../../../src/ssr/cohort/pair_hmm.rs) | belongs to the alignment module; composed, not ported here (§4.3, §7) |
| the placement enumeration for interrupted candidates | [`stutter.rs`](../../../../src/ssr/cohort/stutter.rs) | ported |
| the copy-weighted mixture with an outlier floor | [`likelihood.rs:75`](../../../../src/ssr/cohort/likelihood.rs) | **ported as §2.1's shared formula, for both paths**, with the outlier spread made per-locus rather than per-cohort (§4.5) |
| the censored term | — | **new** (§5); production discards partial observations on both paths |
| the swappable model interface | [`read_model/mod.rs`](../../../../src/ssr/cohort/read_model/mod.rs) | shape ported — the seam is what let production swap models with zero test re-baselining, and §4.5.1 needs it |
| Model B, the comparator | [`classic.rs`](../../../../src/ssr/cohort/read_model/classic.rs) | **test-only**, exactly as production keeps it: an independent implementation is worth more as an oracle than as an alternative |

---

## 10. Deferred, with a recommended home

- **A per-read mismapping mixture.** The right treatment of a mismapped read is §2.1's third term with
  a per-read weight, not a larger error probability. It cannot be expressed against the current
  evidence (§2.3), and enabling it means either splitting an observation's identity by mapping-quality
  band — which multiplies the table — or carrying per-read mapping quality to the caller. GATK's cap
  on per-read discrimination is the cheap approximation of the same idea and could be adopted first.
  **Home:** here, once someone measures how much mapping quality varies within an observation; the
  evidence already carries `mapq_sum_sq`, so the spread is free to report today.
- **A general overdispersion device** — freebayes' read-dependence factor or samtools' dependency
  coefficient (§3.8). **Home:** here, but only once the site filter of §3.7 exists, because
  the two repair the same symptom and one of them names its mechanism.
- **A per-base error rate on the STR path.** The evidence carries a per-read error probability, which is the wrong unit for a per-base substitution term (§4.3); carrying per-base qualities to the caller would make the question askable. The moments as filled are unconsumed
  ([`locus_generation_ssr.md`](locus_generation_ssr.md) §3), and HipSTR does use them. **Home:** here,
  behind the model seam, if a measurement shows the fitted per-stratum rate losing calls the qualities
  would have made.
- **A part-repeat stutter estimator**, replacing production's fixed 5% share, and **untying the
  whole-repeat and part-repeat one-step shares**, which production ties and HipSTR does not (§4.2).
  **Home: unowned, and that is the finding.** The follow-up is named in the model comparison's report
  and in [`alignment.md`](alignment.md) §5.2, but
  [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) has never claimed it — its part-repeat reads
  go to a guard bucket it calls a diagnostic rather than a parameter. **Somebody must claim it before
  it can be scheduled**, and the parameters fit is the right owner. This document's obligation is to keep
  the two parameters separate in the model so the estimator has somewhere to land.
- **A purity adjustment on the slippage level** for an interrupted candidate (§4.6). **Home:** the
  interrupted-repeat work, which has to say what an interruption is for genotyping before it can say
  what it does to slippage.
- **Alleles longer than a read**, by counting reads that fall entirely inside a tract — GangSTR's
  approach. **Home:** already claimed, and **deferred rather than rejected** on our data
  ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §7); §5's censored term is the cheaper half
  of the same capability and is being built.

---

## 11. Open questions

**Q1 — what the calibration scale does on a sample whose noise is off the pre-pass's ladder.** Two of
five real alignments ask for a noisier class than the model covers and are refused, at 0.42% and
0.49% of sites, and what they are asking for fits a duplication the reference does not carry
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1). Such a sample gets a one-rate
answer, and §3.2 would then scale that library's charged errors by a ratio computed against a rate that is
itself absorbing an artefact. *Leaning: scale anyway and report it,* because the alternative — trusting
the instrument on the sample we have most reason to distrust — is worse, and because the refusal is
already reported. **Settled by:** fitting the scale on those two tomato accessions and on three that
were not refused, and comparing the four numbers; if the refused samples' scales sit far from the
others, the scale needs a floor and a ceiling rather than a free ratio.

**Q2 — closed, 2026-08-19 (owner).** It asked whether the read likelihood should carry a second class
of position for the ones that disagree with the reference far more than one rate predicts. **It should
not**, and §3.7 gives the reason: from one locus's reads a duplicated position and a heterozygous
position are the same distribution, so a latent class here would be asked to separate two things its
evidence cannot separate, and every ordinary locus would pay for the attempt. The signal that does
separate them is depth, which the site filter has and this model does not.

**What moved rather than closed is the work.** These positions are still real — 818 in 550,976 on
HG002 — and still have to be dealt with. **Home:** step 11a's artifact filter
([`ng_proposal.md`](ng_proposal.md)), whose own validation has to show what they cost the calls; the
1.41 figure is a fact about a parameter fit, and a fit that comes back wrong and a genotype that
comes back wrong are different claims.

**Q3 — is dividing the error mass by three right for a substitution and by nothing for everything
else?** §3.5 sizes the choice at 1.10 nats, 4.8 Phred, per wrongly-explained read, and the two halves
of production already disagree about it. *Leaning: as specified* — three where the physics gives
three, one where nothing does. **Settled by:** GIAB HG002, single sample, at 5× and at full
depth, with the divisor forced to one everywhere as the comparator — the same instrument
[`calling_priors.md`](calling_priors.md) §2.2 uses, so the two documents' effects are read off one
bench; the indel half of the truth set is what would show the second
branch mattering.

**Q4 — how much does the outlier floor's cohort-wide denominator actually move a genotype?** §4.5
shows it is a tenfold difference in the floor between one sample and a 63-accession panel, and argues
it cancels for reads no candidate explains and not for reads the candidates explain weakly.
*Leaning: fix it regardless of size*, because a per-sample likelihood that depends on other samples
violates goal 4 and the fix is cheap. **Settled by:** re-calling the HG002 tandem-repeat bundle and
one tomato interval with the cohort-wide denominator and with a per-locus one, and counting genotypes
that move; if the answer is zero, the fix is still made and the measurement is what lets it be made
without a re-baselining argument.

**Q5 — closed on the generic path, 2026-08-19, and reopened as a different question on the repeat
path.** It asked whether partial reads should count towards the merge's decision to build a locus.
**On the SNP/indel path: no**, and §5.4.1 gives the mechanism rather than a preference — a
substitution's locus is one reference base wide and has no partial reads, an insertion's reference
span is its anchor base and has none either, and a deletion's partial reads carry the reference
preferentially because a deletion-carrying read crosses the deleted positions without spending read
bases. There is no class left for the change to rescue, and the rule it would loosen is a cost bet
worth 170 ms of a 425 ms merge on one 100 kb tomato interval.

**On the repeat path the filter applies too** (owner, 2026-08-19) — a catalogue tract no sample
varied at is discarded like any other locus. **What is left is not a question but a consequence, and
§5.4.2 carries it:** because a read that ran out stays `Partial` however much of the tract it
covered, a sample whose allele is too long for a read to span shows no complete observation, and a
filter counting complete observations only cannot tell that sample from one that did not vary. The
fix is to test a partial against the reference over the stretch it witnessed, in both halves of the
rule. **Home:** [`cohort_merge.md`](cohort_merge.md) §4.3, as a one-rule amendment on the repeat
path. **Worth measuring afterwards rather than before:** how many catalogue tracts change verdict,
which the STR benchmark bundle answers, since alleles that long are common there and rare at
tomato's tract lengths.

**Q6 — closed, 2026-08-19.** It asked whether the STR substitution term should use the reads' own
qualities the way §3.2 does. It should not, and the reason is a unit mismatch rather than a
preference: the evidence carries a per-**read** error probability and the term needs a per-**base**
rate (§4.3). The two are not interchangeable, and the fitted per-stratum rate is the per-base
quantity, measured inside tracts where the reported qualities are least trustworthy. **What remains
open is a different question with its own home:** whether a per-base treatment would be worth the
evidence shape it needs. §10 records it as deferred rather than open, because nothing can settle it
until someone decides to carry per-base qualities to the caller, and nothing else in this document
wants them.

**Q7 — does modelling contamination at a repeat tract help, and is the cohort's own modal length a
good enough stand-in for what a contaminant shows?** §4.5.1 builds it and turns it on. Two things
are unknown and only the first is cheap. *Does it move genotypes at all* — at the 1 to 3 in 100
contamination that gets flagged in practice, against an outlier weight already set at 1 in 100, the
term may be swamped by the floor it sits beside. *And is the second distribution right* — the seed is
centred on this cohort's modal repeat count, so a contaminant from elsewhere is mis-attributed toward
the mode, the direction in which this term hurts rather than helps. **Leaning: keep it on, and switch it off only
on evidence of harm** — contamination is a property of the sample rather than of the marker, so
correcting for it at ordinary sites and not at repeat tracts would treat one number as two.
**Settled by:** §12's item 17, in that order — the simulated panel first, because it is runnable
today and gives exact truth, and the planted GIAB mixture second, because it needs reads we do not
hold.

**Q8 — where a tract does not slip, does the STR path call as well as the SNP/indel path would?**
§4.8 states the property and names the four things that could break it, the largest being that one
path charges a read's error per base and the other per read — 33 Phred apart at two mismatches over a
20-base tract, and 61 at three. *Leaning: it should hold, and if it does not the per-base/per-read
split is where to look first*, because it is structural, it is the only one of the four whose size is
already known, and it is the one that grows with tract length. **Settled by:** §12's item 18.
**Consequence if it fails:** the copy floors stop being a tuning choice about which loci are worth
the compute and become a correctness boundary, which is a much stronger claim than anything
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5.1 makes for them.

---

## 12. How we know it works

**Unit tests, each pinning a property rather than a value.** The first three are properties
production already tests and ng's port should carry across
([`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs)).

1. **A read matching its allele is dominated by the same-length term.** A read identical to a candidate scores at
   least `same_length_share × (1 − ε)^length` and within 5% of it.
2. **Direction and size are ordered as fitted.** With a contraction-biased direction split, a read
   one repeat short outscores one a repeat long; and one repeat outscores two.
3. **A whole repeat beats a single stray base.** A read a whole repeat longer scores above one a
   single base longer — **whenever the part-repeat share times its one-step share is below the
   whole-repeat product**. *Not simply "at any part-repeat share below the whole-repeat one", which
   suffices only while the two one-step shares are tied:* a whole-repeat share of 0.02 at a one-step
   share of 0.1 gives 0.002, against a part-repeat share of 0.019 at 0.95 giving 0.018 — nine times
   higher, both inside §8's clamps. Production ties them and §10 schedules untying them, so the test
   must survive that.
4. **The stutter distribution sums to one.** Over the full untruncated support, for every period from
   1 to 6 and for direction splits from symmetric to 5:1, the seven parameters produce a
   distribution summing to 1 to floating-point tolerance. **No production test pins this**, and it is
   the test that catches a mis-set same-length share, a one-step share read as its complement, and an
   part-repeat re-indexing off by one — all three of which are silent.
5. **Truncation removes a stated mass, and the builder reports it.** With the cutoffs applied, sum
   the distribution over its whole reachable support and subtract from one; the builder's *reported*
   loss must equal that to floating-point tolerance, for every candidate length, every period, and
   one-step shares across the clamped range. **The test pins that the loss is computed and surfaced,
   not that it is small** — *an earlier version compared it against "a named bound" and named none,
   which is unrunnable.* Its size is a property of the parameters, running from 2 in a million to 2 in
   a thousand across that range (§4.2).
6. **The junk term cancels for a read nothing explains.** For an observation whose emission is zero
   under every candidate, the difference between any two genotypes' log-likelihoods is bit-for-bit
   what it is with that observation removed.
7. **Ploidy generality.** At ploidy 2 the copy mixture reproduces a biallelic calculation done by
   hand. At ploidy 4, with every observation matching one allele, a genotype carrying two copies each
   of two alleles scores **between** the two homozygous quadruples; where the observations are split
   between those alleles it scores **above both**, which is the whole reason a mixed genotype is
   callable. *An earlier version asserted the between-ness at every split, which is false:* one read
   on each allele, emission 0.9 on a match and 0.001 otherwise, gives the mixed genotype −1.59 nats
   against −7.01 for either homozygote. Pin both cases; the second is the one a wrong copy weighting
   breaks.
8. **Order independence.** Permuting the observations, and permuting the candidates, changes no
   genotype's log-likelihood by a single bit.
9. **The aggregation contract holds exactly.** Build a list of individual reads with different error
   probabilities; compute the SNP/indel likelihood by looping over those reads; compute it again from
   the counts and summed logs the merge would have produced; require the two to agree **bit for
   bit**. This is the test that would have caught §1.4's geometric-mean substitution, and it is the
   reason §3.3's formula is shaped the way it is.
10. **The calibration scale reproduces the fitted rate.** Given a read group's minted per-read error
    probabilities and the pre-pass's fitted rate, the scaled probabilities' mean equals the fitted
    rate to floating-point tolerance. A second assertion checks the definitional requirement of §3.2:
    the quantity the pre-pass averaged and the quantity the locus generator mints are computed by
    the *same function*, checked by calling it from both sides on the same read.
11. **The contamination mixture is the plain formula at zero.** With every contamination fraction
    set to zero, the mixture of §3.6 returns §3.3's log-likelihood **to within a few units in the
    last place** — the tolerance stated as a constant in the test, not chosen per case — for every
    genotype, at several allele counts and quality spreads. *Not bitwise: §3.6's identity is exact
    algebra, and §8's evaluation order puts an `exp`/`log` round trip between the two forms.* This is what lets contamination default
    on: it pins that a clean sample is untouched, and it fails the moment anyone reintroduces
    production's extra `(1 − ε)` factor or its allele-count divisor into `own`.
12. **The censored tail is a complement.** For every candidate and every period, the probability that
    the tract is at least `ℓ` plus the probability that it is shorter must sum to **the truncated
    distribution's own total mass** — one minus the loss test 5 reports — to floating-point tolerance.
    **Not to 1**: §4.2 makes slips beyond the cutoffs score zero, and that loss runs from 2 in a
    million to 2 in a thousand depending on the one-step share, orders above any tolerance. *An
    earlier version said 1, which a correct implementation must fail.* And where the constraint admits exactly one length change, the censored
    likelihood equals the complete likelihood at that change, bit for bit — which is the test of the
    tail arithmetic rather than of the tolerance.
13. **A partial never out-discriminates a complete observation.** For the same bases and a stated
    parameter set, the log-likelihood ratio between any two candidates is no larger in magnitude for
    the partial than for the complete observation. Stated on a fixed parameter set rather than as a
    universal claim, because it is a property of decreasing-tailed distributions and this document
    has not proved it for every parameter combination.

**The change measurements, which are not tests but must be run before adoption.**

13a. **The two averages of a read's error — DONE, 2026-08-24.** Per read group on a real cohort,
    the geometric and the arithmetic mean of the minted per-read error over the sites the
    error-rate histogram counts, and the ratio between them. Not a unit test, because the property
    is about a walk over real reads and no fixture can stand in for one:
    `examples/ng_minted_error_means.rs`. **Result: 25.2 on the 63-accession tomato cohort and 44.1
    on HG002 at 300×** — §3.2 carries the table and what it means. The same run is what checks
    §3.2's site-set requirement on real data: the reads the accumulator counts and the reads the
    walk emitted at those loci agreed exactly, 172,616,054 on HG002 and on all 63 tomato
    accessions.

14. **The dropped multinomial coefficient.** Compute both forms over the same merged records and
    report which genotypes move and by how much, on GIAB HG002 at 5× and at full depth. §3.4 argues
    from arithmetic; this is what turns the argument into a number. **A change to §3.3 that does not
    move anything on this bench has not been tested by anything that matters.**
15. **The end-to-end regression, SNP/indel path.** GIAB single-sample genotype accuracy at true
    variants, and the count of true homozygous-variant sites called heterozygous — the same two
    numbers [`calling_priors.md`](calling_priors.md) §12 uses, so the prior's effect and the
    likelihood's are read off one bench and can be attributed apart.
16. **The end-to-end regression, STR path.** The GIAB HG002 tandem-repeat bundle, scored on genotype
    accuracy given detection (`benchmarks/ssr_hg002/src/prior_genotype_accuracy.py`). **It is a
    regression guard, not a discriminating test**: production's STR caller reaches 64% of
    truth-variant loci at 50× and 28% at 15×, and the loci it reaches are the well-covered ones where
    the read likelihood is sharp and every model agrees. It will catch a change that breaks STR
    genotyping and it will not say whether this model is the right one at low depth. That gap is the
    state of the evidence, not an oversight.

17. **Does contamination at a repeat tract earn its place?** Two rungs, and the first is runnable
    today. **Simulated:** the generator that chose Model A already builds a cohort with known
    genotypes ([`sim.rs`](../../../../src/ssr/cohort/sim.rs)); add a contaminating individual at a
    known fraction and score genotype concordance and calibration with §4.5.1 on and off, sweeping
    the fraction over 0, 1, 3 and 10 in 100 and the depth over the committed range. It is the same
    harness that produced §4.1's table, so the comparison is on proven ground. **Planted on real
    reads:** mix reads from an unrelated GIAB sample into HG002 at those same fractions and score
    against the HG002 tandem-repeat benchmark. **That second rung needs a sample we do not hold** —
    `benchmarks/ssr_hg002/` is HG002 alone — so it is an acquisition to schedule rather than a run to
    queue, and the contaminant must be **unrelated**: a trio member shares half its alleles and would
    measure the easiest case as though it were the general one. **What would switch it off:** the term is on by
    default, so this measurement is looking for harm rather than for benefit — a loss of concordance
    at a fraction the pre-pass would flag, or any movement at all at a fraction of zero, where
    §12's item 11 says there must be none.

18. **Does the STR path lose anything at a tract that does not slip?** Take the loci whose fitted
    slippage level is at or near zero — short tracts, and any stratum the parameters fit puts below a
    named threshold — and call them **both ways**: down the STR path, and down the SNP/indel path by
    treating the tract as an ordinary stretch of genome. **Score each against truth rather than
    against each other**, because §4.8 makes the property about accuracy and the two priors differ on
    purpose. Two data sets: the simulator, where truth is exact and the slippage level can be set to
    zero outright; and HG002's tandem-repeat benchmark, where the short tracts are real. **What it
    would show:** if the STR path is the weaker one, report the loci that moved and how many bases
    their reads differ at — §4.8 predicts the gap is concentrated where a read differs at two bases
    or more.

19. **Do the stutter parameters want to be fitted per locus?** Fit them per locus and compare against
    the frozen per-stratum values, on both STR benchmarks and at several depths.
    **First, the distribution of the disagreement** — for every locus with enough reads to fit,
    how far its own fitted values sit from its stratum's. The hypothesis is that most land close and
    a minority sit a long way off; **report the shape of that distribution, not its mean**, because
    the whole case for adaptation rests on the tail.
    **Then, whether it changes calls, and where.** Score both ways against truth: HG002's
    tandem-repeat benchmark, where truth is an assembly and depth can be subsampled down a ladder;
    and tomato's recurrence-based standard
    ([`silver_standard.py`](../../../../benchmarks/ssr_tomato1/scripts/silver_standard.py)), at 3
    reads a position over 63 accessions. **The depth ladder is the experiment**, not a robustness
    check: adaptation should help where a locus has reads to spare and hurt where it does not, so
    report accuracy against reads-at-the-locus and find where the two curves cross.
    **And sweep the shrinkage**, since it is the knob that trades the two off — production's 50
    pseudo-counts for the shape and 20 slipped reads for the level are inherited values nobody here
    has measured. **What would settle it:** adaptation is adopted if it beats the frozen values at
    some depth and does not lose to them at any, with the shrinkage that achieves both named.
