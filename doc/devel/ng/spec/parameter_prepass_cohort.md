# ng — the parameter pre-pass: the cohort gather

*Design spec, 2026-07-31. **No code yet — this settles the design.** Companion arch and plan docs do
not exist yet. One of five documents covering ng step 4; the shared framing is in
[`parameter_prepass.md`](parameter_prepass.md), which walks each sample once and produces that
sample's summary. **This document is about what happens afterwards, reading those summaries and
nothing else.** It fills `CohortEstimator`
([`../arch/ng_step_interfaces.md:351`](../arch/ng_step_interfaces.md)), whose signature is unchanged
by any of this. `src/ssr/` and `src/pileup/` are frozen production: everything said about them here
is a record, not a change.*

---

## 1. Scope — what this document is, and when it runs

Some quantities a caller needs cannot be got out of one sample, however deeply it is sequenced.
**How variable is this population? Is this library contaminated? Which of these samples are
relatives? Do these two libraries share a chemistry?** Each is a *comparison between* samples, so no
amount of data about any one of them contains the answer.

This document is about those, and it is the pass in which they are computed.

**It is not, however, "the cross-sample half" of a clean split.** The per-sample rates — error rate,
heterozygosity, distance from the reference, slippage — might also be fitted here, if they turn out
to be better fitted from the censuses than from the per-sample histograms
([`parameter_prepass.md`](parameter_prepass.md) §4.1 sets that comparison up and does not resolve it).
Nothing below depends on how it resolves. What is settled is that **the quantities in §3–§7 can only
be computed here**, because there is nowhere else the samples meet.

**It runs inside cohort variant calling, not at the end of the walk.** That walk finishes by producing each
sample's summary; this step reads every sample's, once, and touches no read. The two passes are one
statistical design and two separate runs, and they are described separately because of length,
not architecture.

**How a summary gets from one to the other is neither document's business.** It may be held in
memory, serialized per sample, or folded into whatever the pipeline already persists — that is the
pipeline's decision and is settled elsewhere. What both passes require is a property rather than a
mechanism: **this step must reach every sample's summary without walking the reads again**, and the
summary must carry enough identity for §2's compatibility check to be possible.

**Goals.**

1. Estimate the parameters that need the whole panel: **the cohort's diversity** (§3), **the
   frequency spectrum** (§4), **contamination** (§5), **relatedness** (§6) and **which read groups
   share a chemistry** (§7).
2. Do it **without a burn-in** — no second traversal of the reads, ever. Everything here reads the per-sample
   summaries. This is the constraint that shaped those accumulators, and it is the one thing
   in this document that must not be given up quietly.
3. Emit them, with the same provenance and uncertainty marking the walk uses ([`parameter_prepass.md`](parameter_prepass.md) §6), into the
   `ModelParams` the caller consumes as frozen inputs.

**Non-goals.**

- **Genotyping.** As upstream, nothing here calls a variant.
- **Designing the caller's priors.** This step estimates quantities. What the cohort caller builds
  out of them — how a diversity number becomes a genotype prior, what form that prior takes — is the
  caller's design and not settled here. The distinction is easy to lose and §3 is where it usually
  goes.
- **Re-deriving a per-sample parameter that has already been fitted.** The per-base error rate,
  heterozygosity, distance from the reference, slippage and inbreeding all belong to one sample or
  one read group. This step consumes them and does not revisit them. **It may, however, be the pass
  in which some of them are *first* computed**: whether the rates are fitted from the per-sample
  histograms or from the censuses is an open comparison
  ([`parameter_prepass.md`](parameter_prepass.md) §4.1), and a census is cross-sample data that arrives
  here anyway. Nothing in this document depends on which way that resolves — it consumes the
  parameters wherever they were produced. Two things do belong here regardless: a read group too thin
  to fit on its own (§7), and every quantity that is a comparison between samples.

**It does not:**

- read a single read, or open an alignment file;
- require every sample to be present at once beyond what a small gather needs — the summaries are
  small ([`parameter_prepass.md`](parameter_prepass.md) §6) and the arithmetic over them is not the expensive part of anything;
- decide the order in which samples were processed, or depend on it (§8).

---

## 2. What it receives

Everything below reads four things out of each sample's summary. They come from the walk and are
restated here only so this document can be read on its own.

| from the per-sample walk | what it is | used by |
|---|---|---|
| per-read-group parameters | the error rate and the stutter parameters — the chemistry, each with its provenance and the count of observations behind it | §7 (grouping), §5 (contamination) |
| per-sample parameters | heterozygosity and the homozygous-non-reference rate, likewise with provenance and evidence count. **These are per sample, not per read group**, because they are counted over whole sites and a read-group histogram has split them ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §1) | §3 (diversity), §5 |
| the sample's inbreeding coefficient | `F_autozygosity`, from the windowed accumulator | §3 (diversity) |
| **the census sites** | reads supporting **each allele** (A/C/G/T, plus a bucket for anything else) at a fixed set of scattered positions, **the same positions in every sample**, kept **per read group** so that chemistry stays separable | §4, §5, §6 — everything cross-sample. All of these want the *sample's* counts, which is the read groups at a position summed — exact, since they are raw counts at one place |
| **the STR census sites** | this sample's reads at each whole-repeat offset from the reference tract length, plus a non-whole-repeat bucket, at a fixed set of STR loci | §3 — the STR diversity |

**The histograms those first two rows were fitted from do not arrive here at all.** `F_autozygosity`, the
heterozygosity and the homozygous-non-reference rate all come out of one sample's windowed histogram
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4), which is reduced to exactly those
numbers at the end of that sample's walk and then dropped
([`parameter_prepass.md`](parameter_prepass.md) §1.3). **For those three, this step receives the
estimates and not the evidence behind them**, which is why nothing here can refit them and nothing
here tries. The censuses are the exception and are the reason this document exists: they *are*
evidence, they arrive raw, and §3–§6 fit from them directly.

**The census sites are what make this document possible.** A summary is a per-sample marginal, and a
marginal cannot say *which samples carried an allele together*. That correspondence is the whole
content of a frequency spectrum, of a contamination estimate, and of a relatedness matrix.
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) keeps it, at a few megabytes per sample, precisely so this step exists without a burn-in.

**And it keeps allele identity, not just a count of non-reference reads** — without which two samples
carrying *different* alternative alleles at one site would look like they shared one. **Deciding each
site's allele set is this step's job**, because only this step sees every sample: it unions what the
samples observed, works out which allele is major (the reference need not be), and drops any site
where the "other" bucket is large enough to mean the position is not a clean substitution. STR loci
do not appear here at all — they are in a **second census of their own**, holding tract lengths
rather than base counts, because a hypervariable multi-allelic locus does not belong in a
substitution-rate statistic ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)
§1.1). Both sets arrive in every sample's summary and are read by different sections here: the
generic one by §3–§6, the STR one by §3's second diversity.

**Trap: the positions must be identical across samples, and that is an assumption to check, not to
trust.** They are chosen by a rule that is a pure function of position, so they are identical by
construction — but only for runs that agree on **all four** inputs to that rule ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)):

| must match | why a mismatch is silent |
|---|---|
| the selection **seed** | a different seed selects a disjoint set; both summaries look well-formed |
| the **reference** digest | coordinates mean different things, so "the same position" is not |
| the **analysed region set** digest | this is the one that will actually happen. A sample run genome-wide and one run under a `--regions` BED select from different domains, and their overlap is neither the intended sample nor a random one |
| the **region-typing** parameters | they decide which loci are STR loci, so they set the boundary the two censuses partition on. Different copy floors put a locus in the generic set for one sample and the STR set for another |

The region set is the likeliest of the four to differ by accident, because a BED feels like a
runtime convenience rather than part of the data's identity. It is not: it defines what population
of sites the estimate describes. **This step must refuse to gather across summaries that disagree on
any of the four** (§8) — a spectrum built from mismatched sites is not approximately wrong, it is
meaningless, and averaging it produces a number that looks reasonable.

---

## 3. The cohort's diversity

**How variable is this population at an ordinary site?** The number that answers it is the **expected
heterozygosity**: how often two copies of such a site, drawn at random from the panel, differ. It is a
property of the population and it is measurable, so this step measures it.

**"At an ordinary site" is load-bearing, not a hedge.** Repeat tracts mutate orders of magnitude
faster than bases do, so the population's diversity at STR loci is a **different number**, not a
correction to this one. A consumer that applies the number below to a repeat tract will badly
understate how many alleles to expect there.

**So this step estimates two diversities, from two censuses.** The generic set, which excludes
repeat tracts, gives the number below. The STR set gives the STR one
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §1.1). They are computed the
same way and must be emitted as **separate parameters that a consumer cannot confuse** — today's STR
path uses a SNP-scale constant for its own, never measured, which is the failure this separation
exists to prevent recurring.

**The STR one is a gene diversity, not a heterozygosity.** A repeat tract is multi-allelic by length,
so "how often do two copies differ" is the right question but the two-allele arithmetic below is not.
Use the chance that two copies drawn at random from the panel carry different **lengths**, which
reduces to the same thing where a locus happens to be biallelic. The inbreeding correction applies
identically: a selfing crop's STR loci look less diverse than the population is, for the same reason
and by the same factor.

**What is done with it afterwards is not this step's concern.** A cohort caller needs some measure of
how variable a population is before it can judge whether a rare allele is surprising, and how it
turns that measure into a genotype prior is its own design. What is worth recording is only that
today the number is not measured at all: the STR path takes a fixed `SFS_THETA = 0.01`
([`src/ssr/cohort/freebayes_emit.rs:42`](../../../../src/ssr/cohort/freebayes_emit.rs)), freebayes'
default, commented *"Fixed, not a per-run knob"*, while the interfaces list a population diversity
among the panel parameters this step is supposed to estimate
([`../arch/ng_step_interfaces.md:347-349`](../arch/ng_step_interfaces.md)). Nothing estimates it.
This section supplies the number; the caller decides what it is worth.

**Each sample contributes its own observed heterozygosity, and inbreeding is not an optional
correction on it — it is load-bearing.** Inbreeding suppresses observed heterozygosity, so a selfing
crop's samples look far less diverse than the population is. The relation is the textbook one:

```text
observed heterozygosity  =  expected heterozygosity × (1 − F)

  Hobs = Hexp (1 − F)          so this step computes
  Hexp = mean over samples of   Hobs(sample) / (1 − F(sample))
```

**and it is the ratio estimator ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6.1) for `F_hom_excess` read backwards.** That estimator is `F = 1 − Hobs/Hexp`
— the same equation solved for the other unknown. Which raises a trap worth stating plainly:

> **Do not take an inbreeding coefficient from the ratio estimator and then compute expected
> heterozygosity from it.** That
> is circular: the ratio estimator *needs* an expected heterozygosity to produce `F_hom_excess`, so
> feeding its answer back in returns whatever you assumed. The runs-of-homozygosity estimator has no such problem,
> because it reads `F_autozygosity` off the **genomic distribution** of heterozygosity — long homozygous stretches
> against ordinary ones — and never needs a population expectation. **This is a constraint the gather
> places on the open choice of inbreeding estimator** ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §11): the runs estimator is what makes
> the cohort's diversity estimable at all.

**Build it from observed heterozygosity, not from the non-reference rate, and that is deliberate.**
Diversity means polymorphism *within the panel*, not distance from the reference accession. A site
where every sample is homozygous for the alternative allele is not a variant at all — it is a place
where the reference carries the odd allele — and it contributes nothing to observed heterozygosity,
so it contributes nothing here. Estimating diversity from "how often we see a non-reference allele"
would count every quirk of the reference accession as cohort polymorphism, and on a crop reference
that is a large number.

**Soft: how good this number is on a selfing crop is unmeasured.** Inbreeding, drift in a small
effective population and selfing all pull a population away from textbook expectations, and
`Hobs = Hexp (1 − F)` is a textbook relation. **Settled by:** comparing it against the allele counts
from §4's spectrum on the tomato cohort — the two are independent routes to the same quantity, so
they should agree. Until then the number is an improvement on an unmeasured constant, which is a low
bar and the honest claim.

---

## 4. The frequency spectrum

**Across the sites that vary, how many copies of the allele does the panel carry?** Most variants are
rare and few are common, and the shape of that fall-off is a richer description of a population's
variation than any single number.

**Why the per-sample rates cannot give it.** Those rates already *are* a frequency
spectrum — for a sample of two chromosomes: every site carries none, one or two alternative copies.
Recovering the panel's spectrum from that means inverting a downsampling, and downsampling is
many-to-one, so it cannot be inverted. What survives is the one number §3 estimates.

**The census sites can, which is why the walk keeps them.** The method built for our coverage estimates
the spectrum from genotype likelihoods by expectation-maximization **without calling genotypes**
(ANGSD's `realSFS`, *Bioinformatics* 2015, and its stochastic-EM successor, *Genetics* 2022) — the
same "sum over the genotype" principle as [`parameter_prepass.md`](parameter_prepass.md) §3, one level up — and it is unbiased at low
coverage, where a spectrum built from genotype calls is not. What it consumes is a per-site
likelihood over the **sample allele frequency across all individuals**: exactly what the census sites
preserve and no summary can.

**Precision.** A target of roughly ten thousand sites variable across the panel puts a couple of
thousand in the singleton class under a neutral shape and tens in the high-frequency tail — ample for
a prior, and far more than ample for §3's single number. [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5 sets the selection threshold
that delivers it.

**A third route exists, and is recorded only so it is not mistaken for a cheap win.** The *spatial*
pattern of heterozygosity carries more than its average: its density along a diploid genome tracks
the local time to the most recent common ancestor of the two copies, which is what PSMC (Li & Durbin
2011) infers a population-size history from, and a population-size history determines an expected
spectrum (Polanski & Kimmel 2003). The windowed accumulator of [`parameter_prepass_generic.md`](parameter_prepass_generic.md) is already that input. **But the two
routes are known to disagree:** Beichman, Phung & Lohmueller (*G3* 2017) ran both on the same human
populations and found the spectrum-fitted model predicted the observed spectrum **9 log-likelihood
units** from optimal while whole-genome-derived models came in at **152 or worse**, attributing the
gap to unmodelled demography rather than to statistical power. Our case is harder than theirs on
three counts — about 3 reads per plant against the 10–20× that route wants, selfing, and a
domestication bottleneck, which is precisely the unmodelled demography blamed there. Since the census
sites reach the spectrum directly, this route is a curiosity rather than a fallback.

---

## 5. Contamination

**The fitted error rate is not only sequencing error.** A read showing an allele the individual does
not carry may be a misread base, DNA from another individual in the library, or a read from a
paralogous locus mismapped here. **At a site the sample is homozygous for, nothing in one sample
separates them**, so the fitted rate is their sum, and
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2 says so.

> **The criterion below is unchanged, and the estimator now has a home** —
> [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4, the per-site route's generic
> path, which fits the population frequencies this criterion needs over the same loci in the same pass
> (owner, 2026-08-12). It adds a third signature to the two below: a contaminant allele sits in
> **few** reads where a real heterozygote's two alleles are balanced — which needs depth to mean
> anything, since at three reads a site one alternative read is a third. **Not on the STR loci**, where
> stutter lands on the population's own alleles by construction (§4.1 there).

**The cohort separates them, and the criterion is simple.** Contaminant reads carry *real segregating
alleles* and land preferentially at sites polymorphic in the panel. Sequencing errors do not — they
fall uniformly, and the allele they produce has no reason to match anything. So the question is
whether a sample's low-level alternative reads are concentrated on the sites §4 finds variable, and
whether they carry the alleles the panel carries there.

**That is a question about the census sites, asked one sample at a time.** The same object §4 reads
for the spectrum, read again per sample against the spectrum it produced.

**Deferred, with a home: the estimator itself.** The evidence is settled here — it is the census
sites plus §4's spectrum — and the arithmetic on top is not. **Home:** its own spec, or this one when
someone builds it; nothing upstream changes when it lands, which is the property that makes
deferring it safe.

**Why it matters even before it is built.** An uncorrected contamination rate is absorbed into the
per-read-group error rate, where it inflates the estimate and — because contamination concentrates on
polymorphic sites while error does not — makes a single flat rate understate the alternative-read
rate exactly where variants are. The cost is not uniform, which is what makes it worth separating.

---

## 6. Relatedness

**Which of these samples are relatives?** Nothing in the spec's earlier drafts provided this, and the
census sites give it for free: pairwise comparison at positions every sample has in common is the
standard input to any kinship or identity-by-descent estimate.

**It is listed here because the data now exists, not because this step needs it.** No parameter in
this document depends on a relatedness matrix. What does depend on it is downstream: a cohort caller
that assumes samples are independent draws is wrong in a cohort of siblings, and an inbreeding
coefficient interpreted as ancestry is misleading when the "population" is a family.

**Deferred, with a home: the estimator, and the decision about whether the caller should use it.**
**Home:** whichever spec takes on cohort structure. Recorded here so that the next person to want it
finds the data already gathered rather than proposing a new pass.

---

## 7. Sample groups, and read groups too thin to fit

**Two libraries prepared the same way share a chemistry, and pooling them buys precision for free.**
The walk fits every **chemistry** parameter per read group because that is the finest grain available
and the safest default, and it explicitly leaves the fold to whoever knows better
([`parameter_prepass.md`](parameter_prepass.md) §5). This is that place. The sample's own quantities —
heterozygosity, distance from the reference — are never per read group in the first place, so nothing
here folds them.

**What can be decided here that a per-sample pass cannot.** Whether two read groups agree is a
comparison, so it needs both.

**How close counts as agreeing is open, and the obvious answer is no longer available.** An earlier
version pooled read groups whose fitted values overlapped within their uncertainty intervals; the
walk no longer emits one ([`parameter_prepass.md`](parameter_prepass.md) §6), and for a good reason —
these are priors. What is left is a comparison of the values themselves, and the yardstick that fits
the rest of the design is the one §3 of that document uses to set the scan spacing: **two read groups
agree if the gap between them is smaller than a caller could feel.** That makes the grouping rule and
the scan spacing the same soft number, measured once.

*The obvious trap:* that yardstick ignores how much data stood behind each fit, so a read group with
80,000 reads and one with 8 could be pooled on a coincidence. The evidence count arrives with every
parameter (§2) and should gate the comparison — but how is not settled, and §10 keeps it open.

**The panel-pooled fallback lives here for the same reason — but it is the rung *below* a borrow
the sample can make for itself.** A library below `MIN_SITES_TO_FIT` in a sample whose other
libraries were fitted takes their mean, within the sample and needing no other sample
([`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §5.4). Only when
**every** library of a sample is thin is there nothing to borrow from, and only then does a step
holding every read group have something to supply. An earlier version of this paragraph put the
whole read-group borrow here.

**Soft, and unvalidated on our data.** Neither benchmark cohort exercises this: tomato declares one
read group per sample and HG002 one for the whole file ([`parameter_prepass.md`](parameter_prepass.md) §5). The grouping rule is therefore
chosen from first principles and cannot be tested on anything we hold — the same gap [`parameter_prepass.md`](parameter_prepass.md) §10
records for the read-group grain, and it is tested the same way, on the simulator.

---

## 8. Cross-cutting concerns

**Memory and time.** This step reads per-sample summaries, never reads. The parameters are
kilobytes each; the two censuses are one to two megabytes per read group, so roughly a hundred
megabytes across a 50-sample cohort at one read group per sample
([`parameter_prepass.md`](parameter_prepass.md) §6). The arithmetic is a single pass over them. Nothing here is on any critical path, which
is why the walk is allowed to pay for the census sites and this step is not allowed to ask for a
second traversal.

**Errors.** Three failures are distinct and must not be conflated:

- **a summary from a different run** — a different selection seed, reference digest, region-typing parameters, **analysed
  region set** or format version. **Refuse to gather**, loudly, naming which of them differs.
  Silently mixing site sets produces a spectrum that is meaningless rather than noisy (§2), and the
  region set is the one that will differ in practice.
- **a sample missing a parameter** — it was too thin to fit. Not an error: it borrows (§7) and is
  marked.
- **too few samples for a panel quantity.** A two-sample "cohort" has no usable spectrum. Emit the
  quantity as absent rather than as a number backed by almost nothing, because a consumer will use a
  number.

**Determinism.** The gather must not depend on the order the summaries are read. Every quantity here
is a sum or a fit over an unordered set, so this is achievable by construction — sort by sample
identifier before any floating-point reduction, and do not let a parallel gather change the summation
order.

**Concurrency.** Not needed. This is a small single-threaded pass; making it parallel would trade
determinism for a saving nobody will measure.

---

## 9. Deferred, with a recommended home

- **The contamination estimator** (§5). Evidence settled here, arithmetic deferred. **Home:** its own
  spec or a later revision of this one.
- **The relatedness estimator** (§6), and whether the caller should consume it. **Home:** whichever
  spec takes on cohort structure.
- **The grouping rule for sample groups** (§7) — how close two read groups must be to be pooled, and
  how the evidence behind each fit enters that comparison,
  and whether grouping is hierarchical. **Home:** here, once a multi-library cohort exists to test it
  on.
- **Anything requiring per-site data beyond the census sites** — for example a spectrum stratified by
  genomic context, or linkage between sites. The census sites are scattered *precisely so that sites
  are near-independent* ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)), which makes them the wrong object for any question about
  linkage. **Home:** the cohort caller, which sees every site.

---

## 10. Open questions

1. **Is the diversity from §3 consistent with the spectrum from §4?** — OPEN, and it is the best
   check either has, because they are independent routes to the same population. *Leaning:* they will
   disagree on the tomato cohort, because §3's relation is neutral-equilibrium and a domesticated
   selfing crop is neither. **Settled by:** computing both on tomato and comparing. **This is the
   first thing to run once both exist.**
2. **How many samples does a usable spectrum need?** — OPEN. The singleton class is the best-resolved
   and the most informative, and it is also the one that shrinks fastest as the panel shrinks.
   *Leaning:* set a floor and emit the spectrum as absent below it, rather than emitting a noisy one.
   **Settled by:** subsampling the tomato cohort and watching where the spectrum stops being stable.
3. **Should contamination be estimated per read group or per sample?** — OPEN. Contamination and
   index hopping are properties of a library and a run, which argues for the read group; a sample
   swap is a property of the sample. *Leaning:* per read group, because it matches where the error
   rate that absorbs it is fitted ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §3). **Settled by:** whether a multi-library sample ever
   shows one library contaminated and another clean, which needs a cohort we do not have.

---

## 11. How we know it works

1. **The two diversity routes agree on synthetic data.** Simulate a cohort at a known diversity;
   §3's `Hexp` and §4's spectrum must both recover it. This is the primary test and it needs no real
   data. Run it at more than one inbreeding level, since §3's correction is where the two can come
   apart.
2. **The spectrum recovers a known shape.** Simulate under a known demography, and the spectrum
   fitted from the census sites must match the spectrum computed from the true genotypes, within its
   own error. This is the test that the sampling in [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) is unbiased — if the census sites
   were selected on the data, this is where it shows.
3. **Mismatched summaries are refused.** Gather two summaries differing in the selection seed, in the
   reference digest, in the analysed region set, and in the region-typing parameters — four separate
   cases, and the third is the one that will happen in practice. Each must fail and name what
   differs, never average (§8). The fourth is the subtlest: it moves loci between the two censuses
   sets rather than changing either set's size, so both summaries look entirely well-formed.
4. **The gather is order-independent** — same summaries, same parameters, whatever order they are
   read in.
5. **A cohort too small for a panel quantity reports it as absent**, not as a number.
