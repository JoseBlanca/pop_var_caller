# ng step 13 — how sure the caller is: the genotype quality and the site quality

*Design spec, 2026-08-25. **No code yet — this settles the design.** One document, one module:
`src/ng/calling/quality/`.*

***What this is written against.** `SampleGenotypeCall`, `LocusInference`, `CandidateAlleles`,
`Phred` and `SpectrumSeed` are in `main`. `CallingScratch`, `FrozenParameters`, `LocusEvidence`
and the missing-genotype variant of `SampleGenotypeCall` are **step A1 of the calling loop**,
in flight on branch `ng-calling-loop` and not yet merged — this document names them without line
numbers for that reason.*

*Reads on: [`calling_em_loop.md`](calling_em_loop.md) §2 (the loop whose last pass produces two of
these three numbers), [`calling_priors.md`](calling_priors.md) §4 (the fitted concentration this
reuses as a prior on how many variant copies a cohort holds), and
[`cohort_merge.md`](cohort_merge.md) §4.2 (the per-allele read counts the artifact tests read).
Read by: whoever builds emission and site filtering — steps 11a and 11b of
[`ng_proposal.md`](ng_proposal.md) — which have no document yet and which own the threshold this
document deliberately does not set.*

*Production's equivalents are
[`posterior_engine.rs:3371`](../../../../src/var_calling/posterior_engine.rs) and
[`:3466`](../../../../src/var_calling/posterior_engine.rs) (the two model numbers) and
[`vcf/qual_refine.rs`](../../../../src/vcf/qual_refine.rs) (the artifact correction). Everything
said about those files is a record of what they do, not a proposal to change them — `src/vcf/`,
`src/ssr/` and `src/var_calling/` are frozen production.*

---

## 1. What this is

**A caller that says "this sample is A/T here" has told you half of what you need.** The other
half is how much to believe it, and a VCF carries that in two places: `GQ`, one per sample, saying
how sure the caller is of *that sample's genotype*; and `QUAL`, one per site, saying how sure it is
that *anything at all varies here*. They answer different questions and they can disagree — a site
can be certainly polymorphic while one thin sample's own genotype is a coin flip, and a site where
every sample is confidently called hom-ref is certain in every sample and not a variant at all.

**This document says how ng computes both, and where.** The arithmetic is ported from production,
where it was built, measured and twice repaired against the GIAB human benchmark in June 2026. What
is new here is one input ng has and production does not (§5.4), and a placement that makes
impossible a defect production shipped and had to fix (§3.5).

### 1.1 Goals

1. **A per-sample genotype quality** that says how much of the sample's own posterior sits on the
   genotype that was called.
2. **A per-site quality** that says how unlikely it is that *no* sample in the cohort carries a
   non-reference allele, given every sample's reads — a marginal over the cohort, not a product of
   per-sample hom-ref probabilities.
3. **A site quality that does not grow more confident at a false site as depth grows.** This is
   the property production had to add and it is the reason there is a third number rather than
   two: an artifact's reads recur at a steady fraction of the depth, so a model that only asks
   *are there more variant reads than error would give?* gets **more** sure of a false site the
   deeper you sequence it.
4. **Degrade across the committed range** — one sample to several thousand, a few reads a position
   to several hundred (`CLAUDE.md`; the goal every sibling spec carries). §7 and §9 answer both
   ends, and §13's open questions say where the evidence stops.
5. **One quality per site at every moment**, so the number a filter tests and the number a file
   carries cannot be two different numbers. §3.5 records the production defect this exists to
   make impossible.

### 1.2 Non-goals, and what this does not do

- **It does not decide whether a site is emitted.** That is step 11b's, and the threshold is its
  document's to set. What this produces is the number such a threshold would be applied to. Where
  the threshold sits relative to this document's stage is the one ordering constraint this
  document does impose (§3.5).
- **It does not filter, drop or annotate for artifacts beyond the two tests in §6.** ng's
  production sibling also hard-drops on a hidden-paralog likelihood ratio and on DUST; those are
  step 11a's, they are not quality numbers, and nothing here replaces them.
- **It does not recalibrate.** The quality *is* the model's posterior, corrected for two named
  artifact shapes. There is no held-out calibration set, no BQSR-style upstream recalibration, and
  no post-hoc mapping from raw score to observed error rate. That matches freebayes, HipSTR and
  GangSTR; GATK is the field's outlier here ([`ng_proposal.md`](ng_proposal.md) step 13).
- **It does not report an interval or an expansion probability.** GangSTR's `REPCI`, `STDERR` and
  `QEXP` have no counterpart here. Deferred with a home in §12.
- **It does not cover repeat tracts.** §8 says why, and where that work goes.
- **It does not change a genotype.** The two corrections in §6 move `QUAL` only. A sample's called
  genotype, its `GQ`, and the cohort's fitted allele frequencies are what the loop produced and
  are not revisited. Production draws the same line
  ([`qual_refine.rs:35-45`](../../../../src/vcf/qual_refine.rs)), and it is worth stating because
  the obvious next thought — *if the site is an artifact, shouldn't the genotypes change too?* —
  is a different design with a different failure mode.

### 1.3 Vocabulary

Four terms; the last two are this document's.

- **Phred** — a probability written as `−10·log₁₀(p)`. A quality of 30 means a one-in-a-thousand
  chance of being wrong. ng's newtype is [`Phred`](../../../../src/ng/types.rs) (`types.rs:429`),
  and it **refuses** an infinite value rather than coercing one, which §5.3 has to answer for.
- **cohort allele count** — how many copies of a non-reference allele the whole cohort carries,
  summed over samples: `0` at an invariant site, `2N` if every one of `N` diploid samples is
  homozygous for it. The site quality is a statement about this one number being zero.
- **the primary alternative allele** — at a locus with more than one alternative, the one the most
  reads across the cohort reached. The two artifact tests are defined against it and treat the
  locus as if it had two alleles. Exact when it does; an approximation when it does not (§6.4).
- **an artifact** — a systematic reason reads carry a sequence the sample does not: reads mapped
  from a paralogous copy elsewhere in the genome, a recurrent context error of the chemistry, a
  misaligned indel. What makes it artifact rather than noise is that it recurs at a *fraction of
  the depth* rather than at a fixed small rate, which is what §6's tests are built to see.

---

## 2. The three numbers, in one place

| | what it says | scale | computed from | where (§3) |
|---|---|---|---|---|
| **genotype quality** | how much of this sample's posterior sits on the genotype called | Phred, capped at 99 | the sample's posterior row over candidate genotypes | in the loop's final pass |
| **site quality, before correction** | how unlikely it is that no sample carries a non-reference allele | Phred, capped at 9999 | the whole `samples × genotypes` likelihood table, plus the run's fitted spectrum | in the worker, once the loop has stopped |
| **the two artifact penalties** | how much of that confidence the shape of the variant reads does not support | Phred, subtracted | nine pooled counts (§6.3) | in the ordered output stage |

The site quality a file carries is the second minus the third, floored at zero.

---

## 3. Where each number is computed, and why

**This section is the design decision the rest of the document hangs on**, and it is not the
obvious one. The obvious placement — compute every quality downstream, after calling, from the
called variants — is what production does for the artifact correction, and it is where its one
shipped defect came from. It is also, for two of the three numbers, not available: the inputs are
gone by then.

### 3.1 The genotype quality is taken during the loop's final pass

**Not after it.** `CallingScratch::posterior_row` (step A1) is **one** genotype-length buffer, and
each sample in turn is scored into it. When the last sample has been scored, the earlier samples'
posteriors no longer exist anywhere. A quality computed "after the loop" would
have to keep a `samples × genotypes` posterior table — a second buffer the size of the largest one
already allocated, held for the whole locus, to produce one number per sample.

**Decision:** this module owns the arithmetic and the cap as a function of one posterior row;
[`calling_loop.md`](../impl_plan/calling_loop.md)'s step C3, which already walks every sample once
more to take the argmax, calls it as it goes. *Rejected:* widening the scratch to keep the
posteriors, for the memory above; and letting the loop inline the formula, which puts the cap in
two places the first time anything else needs a genotype quality.

**And it is taken once, not every pass — which is a cost point, not a correctness one.** A quality
computed on pass 3 of 50 is a perfectly good quality; it is simply thrown away when pass 4 rewrites
the posterior row. Computing it every pass would multiply one logarithm per sample by the pass
count and keep only the last result. There is no argument to make beyond that.

### 3.2 The site quality's baseline is computed in the worker

Its first step collapses the whole `samples × genotypes` log-likelihood table into a small
per-sample kernel (§5.2). That table is `CallingScratch.lg_table`, per-worker scratch, overwritten
at the next locus — so a downstream stage would need it carried on every called locus in flight.
At 3,000 samples and six candidate alleles that is 21 genotypes a sample, about half a megabyte per
locus, multiplied by the run's look-ahead ([`run_streaming.md`](run_streaming.md) §3.5). Against a
project whose memory thesis is trading resident bytes for sample-count scaling, that is the wrong
direction.

**The second reason is cost, and it is the stronger one.** The convolution in §5.2 is quadratic in
sample count (§9). The run's skeleton yields segments **in genome order from a pool of workers**
([`run_streaming.md`](run_streaming.md) §3.5), so anything downstream of the yield is serial. A
quadratic-in-cohort-size computation on the serial path would put the whole run behind one thread
at exactly the cohort sizes where it costs most.

**Decision:** the baseline is computed in the worker, once the loop has stopped, from the table the
loop already built. It needs nothing the worker does not already hold —
`FrozenParameters` carries both the ploidy and the fitted spectrum (its `ploidy` and `prior_seed`
fields, step A1).

### 3.3 So does the nine-number artifact summary

The two artifact tests need per-allele read counts and the called genotypes. The read counts live
in the merge's `SampleSupport`, borrowed through `LocusEvidence`, and are gone once the locus is
released.

**But their input reduces to nine numbers, and that is what makes the stage in §3.4 affordable.**
Reading production's refinement ([`qual_refine.rs:83-115`](../../../../src/vcf/qual_refine.rs)),
everything the two tests consume is pooled across samples except one accumulated scalar:

- from the reference allele: how many reads, how many on the forward strand, how many placed in the
  left half of their read;
- the same three for the primary alternative allele;
- the total reads at the locus across every allele;
- how many alternative-allele reads the **called genotypes** lead you to expect — a sum over samples
  of `(copies of the primary alternative in that sample's call ÷ ploidy) × that sample's depth`;
- which allele the primary alternative is.

**Nine numbers per locus, independent of cohort size.** Everything cohort-shaped is summed away in
the worker, where the evidence and the genotypes are both in hand.

### 3.4 The two penalties, and the final quality, are the first stage on the output stream

Everything left is arithmetic on those nine numbers plus the baseline: two binomial tail
probabilities and a subtraction, a few dozen floating-point operations per locus. Putting it on the
serial ordered stream costs nothing measurable and buys two things.

**It is the part most likely to move.** The two tests are where the field disagrees, where the
tuning lives, and where §13's open question sits. A stage that can be run over the same called
stream at two settings, without re-running the loop, is the shape this repository keeps needing —
the same reason [`calling_priors.md`](calling_priors.md) put two priors behind one seam.

**And it is where the emission threshold has to sit** (§3.5). Keeping the correction and the gate
adjacent, in that order, is what makes the ordering checkable rather than conventional.

**Decision:** the module exposes the penalties and the corrected quality as functions of the
summary; the *stream wiring* — which stage runs where, and what the emission threshold is — belongs
to step 11's document, not this one.

### 3.5 One quality field, mutated in place — and the defect that makes this a rule

Production keeps the corrected value nowhere: it recomputes it at VCF-encode time from the engine's
baseline ([`record_encode.rs:260`](../../../../src/vcf/record_encode.rs)). Between the correction
shipping (`68d7a181`, 2026-06-11) and the repair (`71e2338a`, 2026-06-27), the `--min-qual`
emission gate compared the **engine baseline** while the **corrected** number went into the QUAL
column — so sites were emitted `PASS` carrying a written `QUAL` of 0. On the GIAB HG002 benchmark
that was 40 false positives at 30× and 64 at 50×; routing both through one function took them to 14
and 14, for at most one true positive lost.

**Decision: the called locus carries exactly one quality at every moment.** The worker writes the
baseline into it; the stage of §3.4 overwrites it with the corrected value and records the two
penalties beside it. There is never a second quality field for anything to read by mistake, and
**nothing between the worker and the stage may read the field** — an invariant worth an assertion
rather than a comment.

**The baseline is not lost by this.** It is `corrected + allele_balance + strand_and_position`,
except where the subtraction floored at zero — and a site whose penalties exceed its baseline is
one no threshold would keep, so the arithmetic that is unrecoverable is exactly the arithmetic
nobody needs. *Rejected:* keeping both numbers on the record, which is the shape that failed above;
and recomputing the correction at write time, which is the same shape again.

---

## 4. The genotype quality

**What it is:** the caller scored every candidate genotype for this sample and normalised those
scores into probabilities. One of them won. The genotype quality is how much probability is *not*
on the winner, written as a Phred.

```
GQ = min( cap, −10·log₁₀( 1 − p_best ) )
```

Three details, all inherited:

- **`p_best` is nudged below one before the logarithm.** A sample whose reads make every genotype
  but one impossible produces `p_best = 1.0` exactly, and `log₁₀(0)` is `−∞`. Production clamps to
  one unit in the last place ([`posterior_engine.rs:3391`](../../../../src/var_calling/posterior_engine.rs));
  ng does the same, and the reason it is a clamp rather than an error is §5.3's.
- **The cap is 99**, matching GATK and bcftools so that downstream tools see the range they expect
  ([`DEFAULT_MAX_GQ_PHRED`](../../../../src/var_calling/posterior_engine.rs)). **Inherited, never
  measured against anything** — it is a convention, and marking it soft is the point of saying so.
- **The argmax is over the posterior, and ties go to the lower genotype index.** The table's order
  is fixed ([`genotype_table.rs`](../../../../src/ng/calling/genotype_table.rs)), so this
  is deterministic; a fold that keeps the first strict maximum gives it.

**A sample with no genotype has no genotype quality.** `SampleGenotypeCall` gained a missing variant
for the sample the candidate step declared uncallable ([`calling_em_loop.md`](calling_em_loop.md)
§5.0); such a sample was set aside before the first pass and has no posterior row at all. Its `GQ`
is absent, not zero — the two mean different things in a VCF and conflating them is how a missing
call becomes an uncertain one.

---

## 5. The site quality, before correction

### 5.1 What it is a statement about

**Not "is this sample variant" summed up.** The question is: *given every sample's reads, how
unlikely is it that the cohort carries no copy of any non-reference allele here?*

```
QUAL = −10·log₁₀ P( cohort allele count = 0 | every sample's reads )
```

This matters, and production's own comment records why the earlier formula was replaced
([`posterior_engine.rs:3446`](../../../../src/var_calling/posterior_engine.rs)). The obvious
alternative — multiply each sample's probability of being hom-ref — is not a normalised posterior,
and every hom-ref sample you add multiplies in another factor below one. Under it, `QUAL` grows with
cohort size at a site nobody carries. The marginal does not: adding a hom-ref sample to a
sparse-variant prior adds essentially no evidence either way, so the number stays bounded by what
the few non-hom-ref samples actually justify. This is also what GATK means by `P(AC = 0 | data)`.

**At a locus with several alternative alleles there is one quality, over their union.** "Any
non-reference allele" is the collapse; a triallelic site does not get three qualities. Inherited
from production, and it is what a VCF's `QUAL` column means — *is there a variant here* — as against
`AC`/`AF`, which are per-allele and stay per-allele.

**The samples it runs over are the samples the locus was called on.** A sample set aside as
uncallable (§4) contributed no likelihood row and does not enter the count axis. That is the same
denominator [`calling_em_loop.md`](calling_em_loop.md) §9 uses for the cohort's expected copies, and
using a different one here would make two numbers on the same record disagree about who was in the
cohort.

### 5.2 How it is computed

Four steps, in the order the code runs them, ported from
[`posterior_engine.rs:3466`](../../../../src/var_calling/posterior_engine.rs):

1. **Collapse.** Each sample's row over candidate genotypes becomes a row over *how many
   non-reference copies that genotype has* — `ploidy + 1` entries, three for a diploid, by
   log-sum-exp over the genotypes that share a count. This is the only step that reads the full
   table, and it is why the baseline is computed in the worker (§3.2).
2. **Convolve.** Fold the samples one at a time into a running distribution over the cohort allele
   count: after `s` samples the array holds the log-probability of each total from 0 to `s·ploidy`.
3. **Apply the prior** on the cohort allele count (§5.4).
4. **Normalise** over every possible total and read off the entry at zero.

### 5.3 Three traps, all of them load-bearing

**Run the fold in the linear domain, not the log domain.** The mathematically identical log-domain
version spends its whole inner loop in `exp`/`ln`, and production measured it at **88% of that
path's own time at 200 samples** before rewriting it as a plain multiply-add convolution with a
per-sample rescaling ([`posterior_engine.rs:3588`](../../../../src/var_calling/posterior_engine.rs)).
Port the rewritten version. The scheme: divide each sample's kernel by its own maximum, fold, then
renormalise by the result's maximum, carrying both divisions in a running log scale — which keeps
every live value in `(0, 1]` and leaves one `exp` per kernel entry and one `ln` per sample.

**Track the zero entry separately, in logs.** The entry the whole calculation exists to read is the
one that underflows first: `log P(count = 0) = Σ_s (that sample's log-likelihood of zero copies)`,
a trivial running sum, and at a confident cohort it goes far below anything the rescaled linear
buffer can hold. Read it back from the linear array and a strongly-supported variant reports
`QUAL = +∞`. Keep the exact log-domain value and override the array's entry with it.

**And ng's `Phred` refuses infinity where production clamps it.** `Phred::from_log_prob` returns
`DomainError::PhredInfinite` for a log-probability of `−∞`, deliberately, and its own doc comment
says the consumer's answer is *"cap at its own ceiling and carry on"*
([`types.rs:454-476`](../../../../src/ng/types.rs)). **This module is that consumer and names the
ceiling**: `MAX_SITE_QUALITY = 9999`, production's `QUAL_MAX`
([`record_encode.rs:44`](../../../../src/vcf/record_encode.rs)), inherited and soft. A `NaN`
is not capped — it is a bug in the arithmetic above and must surface as one, which `Phred::try_new`
already does.

### 5.4 The prior on the cohort allele count — and this is where ng differs

The count needs a prior, and the shape is a Beta-Binomial: a Dirichlet on the split between
reference and non-reference frequency, integrated over that frequency. It takes two concentrations.

**Production's are two constants copied from GATK** — `α_ref = 10`, `α_alt = 0.01` for a SNP and
`0.00125` for an indel, each carrying "Revisit against the cohort calibration set" in its doc
comment ([`posterior_engine.rs:112`, `:119`, `:128`](../../../../src/var_calling/posterior_engine.rs)).
Nobody revisited them.

**ng already holds the same two numbers, fitted.** `SpectrumSeed { alpha_ref, alpha_alt_total }`
([`genotype_prior/mod.rs:494`](../../../../src/ng/calling/genotype_prior/mod.rs)) is the run's
frequency spectrum, reached through `FrozenParameters::prior_seed`, and it lands at
`(1, θ)` on a neutral panel — `α_ref = 1` being the neutral `1/p` density written as a Dirichlet
([`calling_priors.md`](calling_priors.md) §2.3, §4). **Decision: the site quality's prior is the
run's seed.** Two constants that were guessed become one number that was measured, on the same
axis, from the same cohort. *Rejected:* keeping production's pair for comparability — a spec that
imports a value marked *revisit this* into a system that has already fitted it is preserving an
accident.

**It is not a free swap, and the size is worth carrying.** The prior sets how much read evidence a
site must produce before its quality climbs off zero. That is the Phred of the prior odds against
the site being polymorphic, and with `n = ploidy × samples` chromosomes it is a two-term closed
form anyone can re-derive:

```
P(count = 0) = B(α_alt, α_ref + n) / B(α_alt, α_ref)      B = the Beta function
prior toll   = −10·log₁₀( (1 − P(count = 0)) / P(count = 0) )   in Phred
```

| samples (diploid) | production `(10, 0.01)` | ng at `θ` = 1 in 1,000 | ng at `θ` = 1 in 100 |
|---:|---:|---:|---:|
| 1 | 27 | 28 | 18 |
| 3 | 23 | 26 | 16 |
| 63 | 16 | 23 | 13 |
| 1,000 | 13 | 21 | 11 |
| 3,000 | 12 | 20 | 10 |

**Both benchmarks in this repository sit at or near the middle column**, which is the case ng will
be run on first. Human diversity is about one variant per kilobase; tomato1's **fitted** `θ` is
6 in 10,000 over its 52 chromosomes ([`calling_priors.md`](calling_priors.md) §4.1), which is about
2 Phred stricter than the middle column at every cohort size in the table. The right-hand column is
the diverse end that same section names — `θ` of 1 in 100 — carried here because §1.1's fourth goal
asks for both ends of the range, not for the corner we happen to test on.

Read the table as: at one sample the two priors agree within a Phred, and they part company as the
cohort grows — at 3,000 samples ng's fitted seed asks for about **8 Phred more** evidence than
production's constants before a site is called confident. That is the property worth having, and
the diverse column shows why: on a panel ten times as polymorphic a site being variant is genuinely
more likely before any read is looked at, by about 10 Phred, and two fixed constants cannot know
that. It is also why §14's oracle runs ng's arithmetic under production's pair *before* it runs it
under the seed — the port and the improvement are separable, and only one of them is a port.

---

## 6. The artifact correction

### 6.1 Why there is one at all

**A caller's genotype likelihood treats each variant-supporting read as an independent rare error.**
Twice the reads at the same variant fraction is twice the evidence, so the site quality grows about
linearly with the number of variant reads. At a real site that is right. At an artifact site —
reads mapped in from a paralogous copy, a recurrent context error — the wrong reads recur at a
*steady fraction of the depth*, so their number grows with coverage too, and **the caller gets more
confident about a false site the deeper you sequence it.**

Measured on GIAB HG002 by production in June 2026: the median quality of its false-positive SNPs
went 1 → 3 → 150 from 5× to 301×, while freebayes held its own false positives near zero at every
depth. A caller whose mistakes look more certain with more data is mis-calibrated, whatever its
mistakes cost.

**The fix is to judge the *shape* of the variant evidence rather than its amount**, because shape
is the thing that gets *clearer* with depth at an artifact. Two tests, each a Phred subtracted from
the baseline, summed as independent evidence that the site is an artifact.

### 6.2 The two tests

**Allele balance — does the variant-read fraction match what the called genotypes imply?** A single
heterozygote should show about half its reads carrying the variant; a cohort where two samples in
sixty carry one copy each should show about that fraction overall. The expected fraction is read
from the **called genotypes**, not from the fitted allele frequency — the frequency adapts to the
artifact and would excuse it, the genotypes do not. The penalty is the two-sided binomial
improbability of the observed split against that expectation.

Two guards, both inherited and both load-bearing:

- **Only a deficit is penalised.** These artifacts present *fewer* variant reads than a real call at
  that frequency would. An excess is a different phenomenon and this test says nothing about it.
- **Skipped when the genotypes expect nearly all variant reads** (expected fraction ≥ 0.9). A
  homozygous-variant sample's handful of reference reads is sequencing error; a binomial against a
  probability near one reads that as a deficit and charges for it.

**Strand and read position — are the variant reads a fair sample of the site's reads?** A real
variant's reads come off both strands and sit at varied places within the read; artifact reads
often pile on one strand or at one end. The penalty is the improbability of the variant reads'
forward-strand fraction, and of their placed-left fraction, against the *reference* reads' own
fractions at the same site — the larger of the two. Using the reference reads as the expectation
rather than a fixed one half is what makes it robust at a site whose coverage is one-sided for
innocent reasons.

**And this test is ramped in, because at two or three variant reads it has no power.** Three reads
land on one strand by chance often enough that the test called genuine low-coverage heterozygotes
biased and charged them a flat 10–17 Phred — harmless where the baseline is in the hundreds,
lethal at 5×. Production scaled the penalty by a ramp on variant-read count: nothing at or below
3, rising linearly to full at 7 or more
([`qual_refine.rs:180`](../../../../src/vcf/qual_refine.rs), commit `959c5c41`, 2026-06-27). On
GIAB that restored recall (5×: 0.640 → 0.706; 10×: 0.885 → 0.913) while holding the medium-depth
false-positive floor at 14 rather than the 28 and 43 that simply removing the test gives.

**The ramp endpoints are `(3, 7)` and they are soft.** They were read off the GIAB alternative-read
distributions — killed real heterozygotes carried 2–3 variant reads, medium-depth artifacts 5 or
more — at one sample, on human data, in June 2026. §13's Q1 is whether they hold on a cohort. In ng
they are named `pub const`s with that provenance in the doc comment, not production's
`PVC_BIAS_RAMP` environment variable: this repository's convention is typed configuration, and an
environment variable that silently changes a quality is the shape §3.5 exists to prevent.

**The allele-balance test needs no such ramp**, and that is a measurement rather than a symmetry:
production's per-record decomposition found it at zero for true heterozygotes at every depth while
growing with depth for false positives. It is the strand test that had the power problem.

### 6.3 What the tests are given

The nine numbers of §3.3, computed in the worker while the evidence and the genotypes are both in
hand. Two rules about who is counted:

- **A sample with no genotype is in neither the counts nor the expectation.** It has no called
  genotype to derive an expected variant-read count from, so counting its reads in the observed
  total while it contributes nothing to the expectation would manufacture an apparent excess — and
  since only a deficit is penalised, it would quietly *weaken* the test rather than break it
  loudly. Excluding it also keeps this cohort the same cohort as §5.1's.
- **The primary alternative allele is chosen on pooled read count across the whole cohort**, before
  any per-sample work, so every sample's counts are gathered against one allele.

### 6.4 What the correction gives up

**It is a two-allele approximation.** Exact at a biallelic site, which is most of them; at a site
with three or more alternatives it tests only the primary one and the reference. A site whose
artifact lives in a minor alternative is not seen. Inherited from production, and recorded here
because it is the kind of limit that reads as a bug when someone meets it in the field.

**It does not help indels at depth.** Production's measurement: high-depth indel false positives
stay confident under both tests, because they present at a variant fraction near one half —
alignment artifacts that look exactly like real heterozygotes to a test built on read fractions.
The mechanism and its fix are upstream of quality, in how indel-supporting reads are assigned —
the same root cause the genotype-concordance investigation found, where true homozygous-variant
indels carry a persistent 15–20% reference-read fraction at every depth because stutter and shifted
reads collapse into the reference allele
([`gt_concordance_vs_giab_2026-07-04.md`](../../reports/gt_concordance_vs_giab_2026-07-04.md)).
Nothing this document does will move them.

---

## 7. One sample, and three thousand

**At one sample** the cohort allele count runs 0, 1, 2 and the convolution is three terms; the site
quality reduces to a closed form, which is production's own single-sample test
([`posterior_engine.rs:4845`](../../../../src/var_calling/posterior_engine.rs)) and is ng's
hand-checkable oracle (§14). Nothing branches on cohort size: one sample is the general case with
`N = 1`.

**At three thousand** the axis is 6,000 long and the fold is quadratic (§9). Nothing about the
answer degrades — the rescaling scheme in §5.3 is what keeps a 6,000-long distribution numerically
alive — but the cost is real and unmeasured above 200 samples (§13's Q3).

**At three reads a position** both artifact tests are near-silent by construction: the strand test
is ramped to nothing below four variant reads, and a binomial on two reads out of five cannot
distinguish much from anything. That is the intended behaviour — a test with no power should
charge nothing — and it means a thin cohort's quality is essentially the model's own, corrected
only where the reads can support a correction.

**At three hundred reads a position** both tests are at full strength and the baseline is large, so
the correction is the part that decides whether a site survives a threshold. This is the depth the
whole correction was built for.

---

## 8. Repeat tracts are not in this document

**Nothing in ng can score a tract yet.** The repeat-tract read-likelihood row and the repeat-tract
candidate path are both unwritten, and [`calling_loop.md`](../impl_plan/calling_loop.md) names them
as its own blocker. A quality for a tract has nothing to be a quality *of*.

**And the mechanism is genuinely different, not merely untested here.** Production's repeat-tract
caller does not correct a posterior for artifact shape; it runs a per-locus **emission model** with
three swappable arms — a heuristic, a likelihood-ratio test on whether a second allele earns its
extra parameter, and a port of freebayes' formula — in `src/ssr/cohort/`, alongside a posterior
genotype quality. Strand and read-position bias are not what goes wrong at a tract; slippage is,
and slippage is already inside the read likelihood.

**Home:** a sibling document, `calling_quality_ssr.md`, written the way
[`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) is written against
[`candidate_alleles.md`](candidate_alleles.md) — inheriting the genotype quality of §4 unchanged
(it is a property of a posterior, not of a locus kind) and replacing §5 and §6. **Not before** the
two blockers above close.

---

## 9. Cost, memory, determinism

**Cost.** The genotype quality is one logarithm per sample. The artifact tests are a few dozen
operations per locus. **The site quality's fold is the only term that grows with cohort size, and
it grows quadratically**: about `ploidy·(ploidy+1)/2 · N²` multiply-adds, which is roughly 27
million at 3,000 diploid samples against 12 thousand at tomato's 63. Production measured this path
only at 200 samples. §13's Q3 carries it, and §3.2's placement — in the worker, in parallel, off
the ordered output path — is what keeps the answer to that question from being a run-wide serial
bottleneck whichever way it lands.

**Memory.** Per worker, four buffers, all of them living in `CallingScratch` beside the ones already
there rather than being allocated per locus: the collapsed per-sample kernels,
`samples × (ploidy + 1)` doubles; the two linear count-axis buffers the fold alternates between;
and the log-domain result. At 3,000 diploid samples that is about **216 KB** in total — 72 KB of
kernels and three arrays of roughly 6,000 doubles. Per locus carried downstream: the nine numbers
of §3.3 and the two penalties, all scalars. **Nothing cohort-shaped crosses the worker boundary**,
which is the property §3.2 and §3.3 were chosen for.

**Determinism.** Same inputs, same bits, at any worker count. The fold runs over samples in the
run's fixed sample order — the same requirement, for the same reason, that
[`calling_em_loop.md`](calling_em_loop.md) §8 puts on the M-step sum, and a reordered
floating-point sum here is quietly different output at another worker count rather than a crash.
The transcendentals use `std`'s exact `exp`/`ln` rather than the approximate math backend, as
production's do, so the quality does not depend on which backend a build selected.

---

## 10. The types

```rust
// src/ng/calling/quality/mod.rs

/// The pooled read counts the two artifact tests read, and the one number the
/// called genotypes contribute. Nine scalars, whatever the cohort size — every
/// per-sample quantity is summed away where the evidence lives (§3.3).
pub struct ArtifactTestCounts {
    /// Which allele the tests treat as *the* alternative: the non-reference
    /// allele the most reads across the cohort reached.
    pub primary_alternative: AlleleId,
    pub reference_reads: f64,
    pub reference_forward_reads: f64,
    pub reference_placed_left_reads: f64,
    pub alternative_reads: f64,
    pub alternative_forward_reads: f64,
    pub alternative_placed_left_reads: f64,
    /// Every allele's reads, summed over the samples the locus was called on.
    pub total_reads: f64,
    /// How many alternative-allele reads the called genotypes lead you to
    /// expect: `Σ_s (copies of the primary alternative in sample s's call
    /// ÷ ploidy) × sample s's depth`.
    pub genotype_expected_alternative_reads: f64,
}

/// How much of the site quality each artifact test took away. Recorded beside
/// the corrected quality so the uncorrected one stays recoverable as their sum
/// without a second quality field for anything to read by mistake (§3.5).
pub struct ArtifactPenalties {
    pub allele_balance: Phred,
    pub strand_and_read_position: Phred,
}

/// Phred cap on a per-sample genotype quality. GATK's and bcftools'
/// convention, inherited, never measured — see §4.
pub const MAX_GENOTYPE_QUALITY: f32 = 99.0;

/// Phred cap on a site quality, and the answer to `Phred`'s refusal of an
/// infinite value (§5.3). Production's `QUAL_MAX`, inherited.
pub const MAX_SITE_QUALITY: f32 = 9999.0;

/// The variant-read count at or below which the strand and read-position test
/// is charged nothing, and the count at which it is charged in full; linear
/// between. Read off the GIAB HG002 alternative-read distributions at one
/// sample in June 2026 — soft, and §13's Q1 is whether they hold on a cohort.
pub const BIAS_RAMP_NO_POWER_BELOW: f64 = 3.0;
pub const BIAS_RAMP_FULL_POWER_AT: f64 = 7.0;
```

`LocusInference` gains two fields: the site quality (`Phred`, written twice — §3.5) and the
`ArtifactTestCounts`. The penalties join it when the stage of §3.4 runs.

---

## 11. Reuse map

| what | production code | how ng reuses it |
|---|---|---|
| genotype quality from a posterior row, with the below-one clamp and the cap | [`posterior_engine.rs:3371`](../../../../src/var_calling/posterior_engine.rs) | **ported unchanged**; called per sample from the loop's final pass rather than from a loop over a retained posterior table (§3.1) |
| the site quality: collapse, fold, prior, normalise | [`posterior_engine.rs:3466`](../../../../src/var_calling/posterior_engine.rs) | **arithmetic ported; the prior's two numbers replaced by the run's fitted spectrum** (§5.4) |
| the linear-domain fold with per-sample rescaling and the exact zero term | [`posterior_engine.rs:3588`](../../../../src/var_calling/posterior_engine.rs) | ported with both numerical devices — they are the two traps of §5.3 |
| `ln Γ` | [`posterior_engine.rs:3406`](../../../../src/var_calling/posterior_engine.rs), [`qual_refine.rs:219`](../../../../src/vcf/qual_refine.rs) | production carries two copies of the same Lanczos approximation; ng writes one, in this module |
| the two artifact tests and their guards | [`qual_refine.rs:52`](../../../../src/vcf/qual_refine.rs) | ported, over the nine-number summary rather than over a VCF record (§3.3) |
| the power ramp on the strand test | [`qual_refine.rs:180`](../../../../src/vcf/qual_refine.rs) | ported; endpoints become typed constants rather than an environment variable (§6.2) |
| the binomial two-sided tail | [`qual_refine.rs:281`](../../../../src/vcf/qual_refine.rs) | **the closed-form incomplete-beta half only** — §13's Q2 |
| one function for the written value and the gated value | [`record_encode.rs:260`](../../../../src/vcf/record_encode.rs) | **the property is kept, the mechanism is not**: ng stores one quality instead of recomputing a correction at write time (§3.5) |

**Parity oracle, and it is one arm of a two-arm test.** Run ng's site quality with production's
`(α_ref = 10, α_alt = 0.01)` on the same likelihood table and it must reproduce production's number.
Then run it with the fitted seed and report what moved — the same *differential rather than parity*
shape [`calling_em_loop.md`](calling_em_loop.md) §10 uses for the repeat-tract loop, and for the
same reason: two arms that differ by a recorded decision must not be compared as though one were a
bug in the other.

---

## 12. Deferred, with a recommended home

- **The emission threshold, and what a site is dropped for.** Steps 11a and 11b of
  [`ng_proposal.md`](ng_proposal.md), which have no document. This one supplies the number and
  pins one thing about the consumer: **the threshold reads the field the stage of §3.4 wrote**,
  never a value recomputed from the baseline. **Home:** step 11's spec, when it is written.
- **Uncertainty beyond a point estimate** — a confidence interval on a repeat count, an expansion
  probability, the `REPCI`/`STDERR`/`QEXP` family GangSTR emits. Nothing here produces them and no
  consumer has asked. **Home:** the repeat-tract sibling of §8, where they would mean something.
- **Bias annotations as VCF fields in their own right.** freebayes and GATK publish strand and
  position statistics for downstream tools to filter on, where ng would publish only the two
  penalties. Whether the raw counts of §3.3 should reach the file too is the VCF writer's decision,
  not this module's. **Home:** the VCF output document, unwritten.
- **A genotype quality that knows the site is an artifact.** §1.2 excludes it: the correction moves
  the site quality only. Whether an artifact site's per-sample qualities should also fall is a real
  question with a different failure mode — it changes calls, not just their confidence. **Home:**
  its own investigation, if a measurement ever asks for it.

---

## 13. Open questions

**Q1 — do the artifact tests hold on a cohort, and at three reads a position?** OPEN, and it is the
one that would change a shipped number. Everything in §6 was measured on **one human sample** at
5× to 301×, in June 2026, by the production caller. Two things are untested: the allele-balance
test's expectation comes from the *called genotypes across the cohort*, so on a 63-sample panel at
three reads a position it is reading a quantity nobody has checked it on; and the ramp endpoints
`(3, 7)` are counts of one sample's reads, where a cohort's pooled variant-read count crosses 7
for reasons that have nothing to do with one sample's power. *Leaning:* port unchanged and measure
before moving anything — the endpoints are the kind of value that invites tuning before there is
evidence. **Settled by:** running the tomato panel (63 accessions, about three reads a position)
and the GIAB trio through the stage of §3.4 with the tests on and off, and comparing what each
test charges against what the truth set says — the same decomposition production's per-record dump
produced, over a cohort instead of a sample.

**Q2 — one binomial-tail implementation or two?** OPEN, low stakes, decide before coding.
Production has both: an exact discrete sum, and a closed-form incomplete-beta tail it switches to
above 2,000 reads, keeping the sum below that **only to hold its own output byte-identical** at the
depths it had already validated ([`qual_refine.rs:260-269`](../../../../src/vcf/qual_refine.rs)).
The two agree to floating point, not bitwise. ng has no byte-identity obligation to production —
its oracle is a differential (§11) — and at cohort scale the pooled read count is above 2,000 in
any case, so the discrete sum would be dead code on the runs that matter. *Leaning:* the
incomplete-beta path only, one implementation, and state the tolerance the differential is checked
to. **Settled by:** the coder, at the keyboard, against the differential's measured agreement.

**Q3 — what does the quadratic fold cost at the top of the cohort range?** OPEN. §9's arithmetic
says 27 million multiply-adds per locus at 3,000 diploid samples; production profiled the path only
at 200, where it is four orders of magnitude smaller, and no run in this repository has exercised
it near the committed ceiling. *Leaning:* it is affordable, because §3.2 puts it in the worker
where it parallelises with the loop it follows, and because the fold is a vectorisable multiply-add
rather than transcendentals. **Settled by:** timing one locus's site quality at 63, 200, 1,000 and
3,000 samples — a benchmark over a synthesised likelihood table, needing no data. **If it is not
affordable**, the honest lever is the count axis: a cohort where the fitted spectrum puts
essentially no mass above a few dozen copies does not need 6,000 entries. That is a truncation with
a stated error bound, not an approximation to be reached for casually, and it should not be built
before the measurement.

**Q4 — should the site quality's prior be the run's spectrum, or the locus's own class?**
RESOLVED as the run's spectrum (§5.4), and recorded because the alternative is tempting. ng fits
`α_alt` once for the run; production splits its constant by allele class, charging an indel eight
times less prior probability than a SNP. **Rejected:** ng's seed is measured and production's class
split is a GATK-inherited constant, so adopting the split would mean re-introducing a guessed ratio
to modulate a fitted number. If the spectrum should differ by class, the place to find that out is
the parameter fit that produced it, not here.

---

## 14. How we know it works

1. **The single-sample closed form.** At one diploid sample the site quality has a three-term
   closed form; assert the computed value against it to floating point. Production holds this test
   ([`posterior_engine.rs:4845`](../../../../src/var_calling/posterior_engine.rs)) and it is the
   one place the whole four-step calculation is checkable by hand.
2. **The differential against production**, both arms (§11): with production's two constants the
   numbers must agree; with the fitted seed, the movement is reported and must match §5.4's table
   in sign and rough size at the cohort sizes tested. A silent agreement in the second arm is a
   failure — it would mean the seed is not reaching the prior.
3. **Adding hom-ref samples does not inflate the site quality.** Take a locus, append twenty
   samples whose reads are all reference, and require the quality not to rise. This is the property
   the marginal formula exists for (§5.1) and the one the rejected product-of-hom-ref formula
   fails; production holds it as
   `site_qual_does_not_inflate_when_no_sample_has_strong_alt_evidence`
   ([`posterior_engine.rs:4885`](../../../../src/var_calling/posterior_engine.rs)).
4. **A confident cohort is large and finite, never infinite.** Drive the evidence hard enough that
   the zero entry underflows the linear buffer and require a capped quality rather than an error or
   an infinity — the §5.3 trap, and the one that will bite whoever skips the exact log-domain zero
   term.
5. **Two different invariances, and they are not the same test.** *Worker count* must change
   nothing **bitwise**: the fold runs over samples in the run's fixed order (§9), so the same
   cohort at 1 and at 16 workers gives identical bits, and a test that allows a tolerance here
   would pass a fold that had picked up the thread order. *Permuting the cohort itself* legitimately
   moves the last bits, because the fold's summation order genuinely changes; assert agreement to a
   tolerance instead — production's proptest uses `1e-6` on the site quality
   ([`posterior_engine.rs:6387`](../../../../src/var_calling/posterior_engine.rs)). Writing the
   second test where the first belongs is how a run-order dependency ships unnoticed.
6. **The false-positive depth curve, which is the whole reason §6 exists.** On GIAB HG002 across
   the depth ladder, the median quality of false-positive SNPs must **not** rise with depth. This
   is the failing state the correction was built against — production's went 1 → 3 → 150 from 5× to
   301× without it — and it is the only test here that would catch the correction being wired in
   but inert.
7. **Recall is not paid for it at low depth.** On the same ladder, at 5× and 10×, the true
   positives surviving a fixed threshold must not fall against the uncorrected arm. Production's
   ramp exists because this test failed once (141 true SNPs at 5×, 80 at 10×); ng inherits the fix
   and must inherit the test that proves it.
8. **One quality, and it is the corrected one.** Assert that nothing reads the quality field
   between the worker and the stage, and that the value a threshold sees is the value the record
   carries — the §3.5 invariant, and the shape of production's shipped defect.
