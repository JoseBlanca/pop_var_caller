# ng — the genotype prior

*Design spec draft, 2026-08-18, amended 2026-08-19 (§4.1: the prior's starting point is read off the
pre-pass's fitted spectrum rather than fixed at the neutral shape). First of three
documents on variant calling, written in pieces at the owner's direction; the siblings are the
**read likelihood models** (one for the SNP/indel path, one with stutter for the STR path) and
the **EM loop** that ties them together. Neither exists yet, and this document is written so it
can be read before them.*

> **⛔ §5 is superseded, and the code no longer matches it.**
> [`population_diversity.md`](population_diversity.md) §4 replaces the repeat tract's prior seed:
> the shape is no longer constructed as a geometric decay from the cohort's modal repeat count, and
> the total is no longer scaled to reproduce a measured repeat gene diversity. It is the **length
> spectrum and concentration the joint repeat fit already produces per stratum**, mapped onto the
> locus's candidate lengths by their offset from the *reference* tract length. Step E2e of
> [`../impl_plan/calling_loop.md`](../impl_plan/calling_loop.md) built that and deleted the
> construction, `SeedDecayPerRepeat`, and `SsrSeedOutcome` with its `DiversityUnreachable` refusal.
> **§5 and its open question Q2 are kept as the record of what was replaced and why** — Q2 asked
> what to do at a locus the geometry could not hold, which is a question the fitted pair cannot
> raise. §4, the ordinary-site half, is **not** superseded: `population_diversity.md` §3 supplies
> its inputs and changes none of its rules.
> *(The status line above said "no code yet" until 2026-08-26; §4 and §5 were both built.)*

*Reads on: [`cohort_merge.md`](cohort_merge.md) — what a cohort observation is, the input to
calling; [`parameter_prepass.md`](parameter_prepass.md) and
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) — the frozen parameters this prior
consumes and does not fit. Production's equivalent code is
[`posterior_engine.rs`](../../../../src/var_calling/posterior_engine.rs) (SNP/indel) and
[`ssr/cohort/em.rs`](../../../../src/ssr/cohort/em.rs) (STR); the shared mathematics already sits
in [`genetics.rs`](../../../../src/genetics.rs). Everything said about those files is a record of
what they do, not a proposal to change them — `src/ssr/` and `src/var_calling/` are frozen
production.*

---

## 1. What this is

**Before any read is examined, some genotypes are more likely than others, and the prior is where
that belief is written down.** At a locus with a set of alleles, for one sample, it produces one
log-probability per possible genotype. The caller multiplies that by what the reads say — the read
likelihood, the sibling document — and normalises; the result is the posterior this caller emits.

The belief has exactly two sources, and both are measured before calling starts:

- **How variable the population is.** In a nearly-invariant population almost every sample is
  homozygous for the reference allele, and a single read carrying something else is more likely to
  be an error than a variant. In a diverse one it is not.
- **How homozygous this individual runs.** A selfing tomato accession is homozygous nearly
  everywhere; an outbred human is not. This is the inbreeding coefficient `F`, and it is a property
  of the sample, not of the locus.

Both come from the parameter pre-pass ([`parameter_prepass.md`](parameter_prepass.md) §1) as frozen
inputs. **The prior fits nothing.** The one thing it does learn — how common each allele is *at this
locus* — it learns from the other samples in the cohort during the EM loop, and §6 says exactly how.

**One term is used throughout and is worth pinning down here.** The prior over a locus's allele
frequencies is a Dirichlet, and a Dirichlet is described by one positive number per allele — `α_ref`,
`α_alt`, one for each. Those numbers are its **concentration**, and the way to read them is as
**chromosomes the prior behaves as though it had already seen**. Their ratio is the frequency it
expects, `α_a / Σα`; their sum is how much conviction that is, since `Var(p) = p(1 − p)/(Σα + 1)`, so
a larger sum is a tighter belief about the same expectation. A concentration of `(1, 0.005)` says
*the alternative allele is expected at about one in two hundred, held with one chromosome's worth of
conviction*. **The name is unhelpful and standard** — it sounds like a description of how
concentrated the distribution is, which is true only of the sum — but it is what the literature and
production's code call it ([`genetics.rs:214`](../../../../src/genetics.rs)), so this document keeps
it. Reading it as chromosome counts is also what makes §6's arithmetic obvious: the cohort's observed
allele copies are added straight onto these numbers, because they are the same unit.

### 1.1 Goals

1. **One prior for both paths, two starting points.** The mathematics of turning a population's
   variability into a genotype probability is the same for a SNP and for a repeat tract. What
   differs is the frequency spectrum each path faces. At an ordinary site most alternative alleles
   are rare, so the one chromosome the reference happens to record is almost always the common one,
   and the starting concentration can put its mass there. Repeat tracts mutate orders of magnitude
   faster: their alternative alleles are not rare, often no length is clearly the major one, and
   the reference's length is one draw among several common ones. So the STR path starts its mass
   somewhere else entirely (§5).
2. **Marginalize, never plug in.** Do not estimate the allele frequency and then use the estimate
   as if it were the truth. §2 gives the mechanism and the measured size of the difference; §3 is
   the shape that follows.
3. **Degrade across the committed range** — one sample to several thousand, a few reads a position
   to several hundred (`CLAUDE.md`, *What this caller has to work on*). Both the starting point and
   the per-locus term cover both ends with no special case: §4.1's projection returns the neutral
   shape at one sample because there is nothing there to move it, and §6's cohort term is exactly
   zero at one sample and swamps the starting point at a thousand.
4. **Multi-allelic and ploidy-generic from the first line of code.** Not because polyploids are
   scheduled, but because the biallelic-diploid shortcut is what production had to unpick to get
   here, and the general form costs a loop bound (§3.3).
5. **Consume frozen parameters and stay a pure function.** Same inputs, same log-priors, at any
   thread count.

### 1.2 Non-goals, and what this document does not do

- **It does not compute read likelihoods and it does not run the EM.** Both are siblings. This
  document defines a function; the EM document decides when it is called and with what.
- **It does not decide whether a site is emitted, or what its QUAL is.** Emission is a separate
  question from genotyping, and production's STR work found the two behave independently
  ([`ng_proposal.md`](ng_proposal.md) step 11).
- **It does not estimate diversity, the frequency spectrum, or `F`.** Those are the pre-pass's
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3, §4;
  [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §6). It does *read* the fitted
  spectrum and rewrite it as a concentration (§4.1) — a change of representation, with nothing
  fitted here that the pre-pass has not already fitted.
- **It does not use relatedness.** The pre-pass estimates it
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6) and this prior ignores it —
  every sample is treated as an independent draw from the population. §10 records where a
  relatedness-aware prior would go.
- **It does not link loci.** No linkage disequilibrium, no phasing, no haplotype prior. Each locus
  is independent of every other.

---

## 2. Integrating over the frequency, not substituting an estimate for it

**The pre-pass supplies a concentration; what this document decides is what to do with it.**
Collapsing it to a single frequency and squaring that is one option, integrating over the whole
distribution is the other, and the two give different genotypes from identical inputs. The
difference between them is the largest measured calling defect in the project's history.
**Integrating is not by itself the repair, though.** It fixes nothing unless the distribution being
integrated over is the right one: §2.3 shows the wrong distribution returning the same wrong answer
as the plug-in route, at more expense.

### 2.1 The two ways to build a prior from a frequency

Hardy–Weinberg turns an allele frequency into a genotype probability. At a biallelic site with
alternative-allele frequency `p`, for a diploid:

```text
P(0/0) = (1 − p)²      P(0/1) = 2p(1 − p)      P(1/1) = p²
```

But `p` is not known — it is close to what the cohort is being asked about. There are two ways to
proceed:

- **Plug-in.** Estimate `p`, call the estimate `p̄`, and substitute it into the three formulas.
- **Marginalized.** Treat `p` as a quantity with a distribution rather than a value, and average
  the three formulas over that distribution, weighting each `p` by how plausible it is.

### 2.2 They differ by exactly the variance, and always in the same direction

The genotype probability is quadratic in `p`, and the average of a curve is not the curve at the
average. For the homozygous-alternative genotype:

```text
E[p²] = p̄² + Var(p)
```

Plug-in uses `p̄²`, so it **undercounts homozygotes by exactly the variance of `p`** — the same
term for hom-ref, by the same algebra — and since the three probabilities sum to one, that mass
goes somewhere:

```text
P(het) = 2·E[p(1 − p)] = 2p̄(1 − p̄) − 2·Var(p)
```

**Plug-in hands heterozygotes `2·Var(p)` — the sum of what it takes off the two homozygotes, and
twice what it takes off each.** This is
not a tuning accident. `Var(p)` is how badly the frequency is pinned down, so the bias is
negligible with a thousand samples and dominant with one sample at low depth — **it is largest
precisely in the corner this caller commits to supporting.**

**The measured consequence, on GIAB (human, three samples, each called on its own, at 5×).**
Production's plug-in prior gave SNP genotype accuracy at true variants of **83.6%**, with **214
sites where a sample carrying two copies of the variant was called heterozygous**. Replacing the
plug-in prior with the marginalized one took the same benchmark to **94.6%**, and those 214 sites
to **8** — with precision, recall and the emitted variant set byte-identical, because the prior
moves genotypes among emitted variants and does not change which sites are emitted. Report:
`doc/devel/reports/sfs_prior_giab_validation_2026-07-04.md`. **That measurement is a fact about
one corner** — one sample at a time, high-quality human data, 5× — and the range it covers is the
low-coverage single-sample end, which is the end it was aimed at.

### 2.3 The trap: marginalizing over the wrong distribution reproduces the bug

For a diploid biallelic site, averaging Hardy–Weinberg over a Dirichlet with concentrations
`(α_ref, α_alt)` gives a ratio that can be checked by hand:

```text
P(het) : P(hom-alt)  =  2·α_ref : (1 + α_alt)
```

Production's plug-in path regularised its frequency estimate with `α_ref = 10`, `α_alt = 0.01`
([`posterior_engine.rs:107`](../../../../src/var_calling/posterior_engine.rs),
[`:114`](../../../../src/var_calling/posterior_engine.rs)), which gave a prior ratio of 22:1 for
heterozygous **at the configuration that matters — one diploid sample's own two copies inside the
estimate**, putting `p̂` near `1/12` and the Hardy–Weinberg ratio `2(1 − p̂)/p̂` near 22. *(Read off
the pseudocounts alone, with no sample in them, `p̂` is about 1 in 10,000 and the ratio is near
2,000:1, which is the same failure an order of magnitude further along. The 22 is the number the
GIAB run actually met.)* **Marginalize over that same Dirichlet and the ratio is 20:1** — the same
wrong answer, computed more expensively.

The whole gain in §2.2 lives in one number: `α_ref = 1` rather than 10. That value is what the
neutral site frequency spectrum — density proportional to `1/p`, most polymorphic sites carrying a
rare allele — looks like written as a Dirichlet, and it puts the ratio at the defensible **2:1** as
`α_alt → 0`. Production fixed it there as a named constant with that reasoning
([`genetics.rs:179`](../../../../src/genetics.rs)).

**So the decision this document records is not "marginalize". It is "marginalize over the site
frequency spectrum", and the starting concentration is where that decision is actually made.**
§4 sets it for the SNP/indel path and §5 for the STR path, and they differ.

---

## 3. The shape: a Dirichlet-multinomial with an inbreeding branch

### 3.1 The random-mating term

**A genotype is a handful of allele copies drawn from a population whose composition we are unsure
of.** Averaging the multinomial draw over a Dirichlet prior on the composition has a closed form —
the Dirichlet-multinomial — so nothing is integrated numerically. For a genotype `g` whose copy
counts are `k_a` (how many copies of allele `a`, summing to the ploidy `m`):

```text
log P_random(g) = log C(m; k)  +  Σ_a [ lgamma(α_a + k_a) − lgamma(α_a) ]
```

`log C(m; k)` is the multinomial coefficient — how many orderings of the copies give this
genotype. A zero count contributes exactly zero to the sum, so the cost is one `lgamma` pair per
allele the genotype actually carries.

The genotype-independent term `lgamma(Σα + m) − lgamma(Σα)` is **omitted** from the expression
above: it is the same for every genotype, so it cancels when the caller normalises. The values are
therefore log-priors up to a shared additive constant, which is what a softmax consumes.
Production's primitive ([`genetics.rs:127`](../../../../src/genetics.rs)) takes flat arrays and
returns exactly this, and ng ports it unchanged (§9).

**The constant cancels in a row and it does not cancel in a mixture, and that distinction is the
one thing in this section that has already cost a defect.** §3.2 combines this term with a second
branch whose value, `α_a / Σα`, is a true probability. Adding a constant to one summand of a
`logsumexp` and not the other changes the answer, so **the random-mating branch has to be put back
on the probability scale before the two are mixed** — §3.2 says how, and what it costs not to.
Production does not, and its own defence is that its default inbreeding coefficient is zero, where
the second branch is not there to be mixed.

### 3.2 The inbreeding branch

`F` is the probability that an individual's two copies of a locus are **identical by descent** —
inherited from the same ancestral copy — rather than two independent draws from the population.
That makes the prior a two-branch mixture, not a correction term:

```text
with probability F      the copies are one copy counted twice, so the genotype is
                        homozygous for allele a with probability α_a / Σα
with probability 1 − F  the copies are independent draws, giving the term in §3.1
```

In log space, for a genotype that is homozygous for allele `a`:

```text
log P(g) = logsumexp( log(1 − F) + log P_random(g),  log F + log(α_a / Σα) )
```

and for any genotype that is not homozygous, only the first branch exists:

```text
log P(g) = log(1 − F) + log P_random(g)
```

**`log P_random(g)` here is the whole probability, including the `lgamma(Σα + m) − lgamma(Σα)`
term §3.1 drops.** The second branch is a true probability; mixing an unnormalised first branch
against it inflates the random-mating side by `Σα(Σα + 1)` at diploidy, and the inbreeding
coefficient then does a fraction of the work it should. **Measured, biallelic diploid at tomato's
fitted diversity of 6 in 10,000, reading the heterozygote to homozygous-alternative prior ratio at
`F = 0.8`:** the model says 0.222 at one sample and the uncorrected mixture gives 0.400 — it
travels 90% of the way from the outbred 2:1 to the right answer. **At cohort scale it barely moves
at all**, because the concentration handed in is the leave-one-out one (§6) so `Σα` grows with the
cohort and the inflation grows with its square: 50 samples reach 3.6% of the way and 1,000 samples
0.09%, so an inbreeding coefficient of 0.8 buys almost nothing.

**Production has this defect and it is live rather than latent.** Its engine mixes the two scales,
and its default coefficient of zero hides it — but the pipeline also hands the engine the
per-sample coefficients the diversity estimator *fitted*, as overrides, and those are not zero on
an inbred panel. **ng corrects it deliberately (owner, 2026-08-22), and it is the one place the
port departs from what it was ported from.**

This is otherwise production's mixture, from
[`posterior_engine.rs:3799`](../../../../src/var_calling/posterior_engine.rs) and
[`:3217`](../../../../src/var_calling/posterior_engine.rs), and the STR side's port of it
([`ssr/cohort/em.rs:290`](../../../../src/ssr/cohort/em.rs)).

**Why the mixture rather than the textbook Wright formulas.** The familiar
`P(het) = 2pq(1 − F)`, `P(hom-alt) = p² + Fpq` is biallelic and diploid. It exists in production
([`genetics.rs:66`](../../../../src/genetics.rs)) but **its calling engine never calls it** — the
only caller in the tree is the hidden-paralog filter's locus score
([`paralog/locus_score.rs:200`](../../../../src/paralog/locus_score.rs)), which is biallelic by
construction. The mixture above says the same thing at two alleles and keeps saying something at
three, four, or at a tetraploid. Since goal 4 is ploidy-generic code, the mixture is the form ng
builds, and the Wright formulas become a **test oracle** rather than a code path (§12).

### 3.3 What "homozygous" means above diploidy — a trap

The mixture's homozygous branch fires when **every** copy is the same allele. At ploidy 2 that is
the whole story. At ploidy 4 it is not: two of four copies can be identical by descent while the
other two are not, and the model above has no state for that — a tetraploid genotype `AABC` gets
the random-mating branch alone, as though `F` were zero for it.

**This is a known and accepted simplification, and it is not this document's to fix.** The
pre-pass already defers the population genetics above diploidy — several identity-by-descent
coefficients instead of one `F` — to a spec of its own
([`parameter_prepass.md`](parameter_prepass.md) §8), with the instruction that its definitions
degrade to the diploid ones at `P = 2`. **The coder should build the prior so the homozygous test
is one function** rather than an inlined comparison, so that spec has one place to change.

---

## 4. The starting point on the SNP/indel path

**At an ordinary site most alternative alleles are rare, which is why the reference allele is
usually the common one, and the starting concentration says so — but only just.**

The concentration has two numbers, and the alternative one is shared out across however many
alternative alleles the site carries:

```text
α_ref     = the reference allele's concentration
α_alt(a)  = (the alternative concentration) / (number of alternative alleles),
            floored at a tiny positive value
```

On a neutral panel those two numbers are `1` and `θ`, where `θ` is the cohort's **expected
heterozygosity at ordinary sites**, estimated by the pre-pass
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3) — not the non-reference rate,
which counts every quirk of the reference accession as cohort polymorphism. **They are not fixed
there**, and §4.1 says what sets them.

Three properties are worth stating because they are why this is the setting and not another:

- **`α_ref = 1` is where the fit lands on a neutral panel, not a tuning knob.** It is the `1/p`
  density written as a Dirichlet, and it puts the het:hom-alt ratio at 2:1 at every realistic
  diversity, because `α_alt` is always small (§2.3). Nothing here privileges the reference as such:
  it holds the large concentration only because it is one chromosome drawn from the population, and
  under a `1/p` spectrum a drawn chromosome carries the site's common allele nearly every time.
  Where that stops being true — §5 — the asymmetry goes with it.
- **No separate invariant-site mass is needed.** With `α_alt = θ` the hom-ref probability comes out
  at `1 − 3θ/2`, which is the genetically correct weight for a site being invariant. Production
  records this as the reason it dropped a separate monomorphic term
  ([`genetics.rs:179`](../../../../src/genetics.rs)).
- **Splitting `θ` across the alternative alleles keeps a site's total polymorphism independent of
  how many alleles it happens to carry** ([`genetics.rs:214`](../../../../src/genetics.rs)). A
  triallelic site is not twice as polymorphic as a biallelic one merely for having a third allele
  in the table.

**One sample is not the fallback case, and the spec used to imply it was.** The pre-pass's expected
heterozygosity is a mean over samples of each sample's observed heterozygosity divided by `(1 − F)`
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3), and a mean over one sample is
that genome's own heterozygosity corrected for its own inbreeding — a measurement of the
population's diversity, not a stand-in for one. **`θ` is fitted at every cohort size down to one.**

**The fallback is for having no estimate at all** — too few sites to fit, or a sample whose `F` is
unknown. Then the pre-pass falls back to a species-range value; production's is `1e-3`, roughly
human nucleotide diversity, marked as a weakly-informative fallback rather than a default for the
estimate ([`diversity.rs:78`](../../../../src/var_calling/diversity.rs)). **Soft: that number is a
human figure and a tomato panel is more diverse than a human one.** ng should carry the fallback as
a named, overridable parameter and record which of the two produced the value in the run's output,
because a run at the fallback and a run at a fitted `θ` are otherwise indistinguishable.

### 4.1 Where the two numbers come from: the fitted spectrum, projected

**The neutral `1/p` density and the neutral frequency spectrum are the same statement written
twice** — once at a locus, once across a panel. Under neutrality the expected number of sites at
which `k` of `2N` sampled chromosomes carry the alternative allele is `θ/k`; draw `2N` chromosomes
from a Dirichlet with `α_ref = 1`, `α_alt = θ` and `θ/k` is what comes out **in the limit of small
`θ`**. The gap is the panel's own chance that a site is polymorphic, about `θ · H(2N)` with `H` the
harmonic number: **3 in a thousand** at `θ = 6 in 10,000` over 52 chromosomes — the diversity fitted
on tomato1 — rising to about **8%** at `θ` of 1 in 100 across two thousand chromosomes. That is close
enough to treat the two as one
object and too far to call an identity, which is why §12's test builds its target as the family's
**exact expected spectrum, computed in closed form**, rather than by writing `θ/k` — and not by
sampling either, for the reason test 5 gives. So there is no choice to make between "the neutral shape" and
"the pre-pass's measured spectrum". They are one object at two sample sizes, and the concentration is
read off the spectrum rather than fixed.

**How it is read off.** The pre-pass estimates the panel's spectrum from the census sites,
regularized toward the neutral `θ/k` shape
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4.1). This document projects that
spectrum onto the two-parameter family above: compute what a candidate `(α_ref, α_alt)` predicts
for the sample allele counts at `2N` chromosomes, and take the pair whose prediction matches the
fitted spectrum. **That is a change of representation, not a second estimate**; nothing is fitted
here that the pre-pass has not already fitted.

**"Matches" has to name an objective, because different ones give different pairs.** Use the
**maximum-likelihood fit of the predicted class probabilities to the fitted spectrum's class
weights** — equivalently, minimise the Kullback–Leibler divergence from the fitted spectrum to the
predicted one — over **all** classes including the monomorphic ones, since those are what pin
`α_alt` against `θ` rather than leaving only its shape identified. The prediction is the family's
**exact expected spectrum** at `2N`, computed in closed form; nothing is simulated.

**The prediction must use §3.2's two-branch sampling, carrying the panel's `F`, and not independent
chromosomes.** A panel's `2N` chromosomes are not `2N` independent draws once its individuals are
inbred, and treating them as such biases `α_ref` **down**, with a fixed sign. On a simulated panel of
26 individuals whose population spectrum is exactly `Dirichlet(1, θ)`, an independent-chromosome
projection returns `α_ref = 0.914` at `F = 0.6` and `0.860` at `F = 0.9`, where the two-branch
version returns `1.000` at every `F` from 0 to 0.9. At tomato's fitted `F` of 0.8 to 0.9 that is a 12
to 14% error in the number the whole section is about. It costs nothing to avoid, and it keeps the
projection and the genotype prior inside one sampling model.

**The same requirement binds one step earlier, and the pre-pass carries it**
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4.1): the *neutral shape the estimate
is regularized toward* must be written in the panel's own sampling too, or the projection is matched
against a reference shape the panel could not have produced. **On tomato this is the dominant feature
of the spectrum rather than a correction to it** — in the cohort VCF over 26 accessions, 10,786 sites
carry the alternative allele on exactly two chromosomes against 5,142 on exactly one, doubletons
outnumbering singletons 2.1 to 1, which no independent-chromosome spectrum produces at any `θ`
because `θ/k` falls monotonically. Inbreeding is not a second-order effect on a selfer's spectrum; it
is its shape.

**The regularizer's own scale is `θ`, and that is what makes the single-sample case work.** The
neutral prior is `θ/k`, not a bare shape: its total mass is the `θ` of
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3, which is each sample's observed
heterozygosity from the windowed histogram
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §5) divided by `(1 − F)` and
averaged. **So the shape is theoretical and the scale is measured, and the two come from different
accumulators.**

**At one sample the projection returns `(1, θ)`, and no test of the cohort size makes it do so.** A census
site's whole content is which samples carried an allele together, and a panel of one has no such
correspondence: no site is variable across it, so the census contributes nothing and the spectrum
is its `θ/k` prior untouched — at that genome's own measured `θ`. **The single-sample case
therefore rests on the per-sample windowed histogram, not on the census sites**, which is worth
knowing when judging what a thin sample can be trusted for. At a thousand samples the census sites
outweigh the regularizer and the panel's own spectrum sets the pair.

**And in the middle the pre-pass hands over nothing at all, which this document has to answer for.**
Below a panel-size floor the spectrum is emitted **as absent** rather than as a thin estimate
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §4.1, §8; its Q3 puts the floor in the
tens), so every cohort from two samples up to that floor arrives here with no spectrum — while a
single sample arrives with one. **On absence the concentration is `(1, θ)`**: the same pair the
one-sample case produces, for the same reason, since a spectrum too thin to emit and a panel with
nothing to fit carry the same information about shape. That is a branch on *absence*, which the
pre-pass forces and which the code must have; it is not a branch on cohort size, and nothing
downstream may test `N`. **One formula covers both ends of the committed range**, and the middle is
covered by the absent case rather than by a third rule.

**What two parameters cannot hold.** A spectrum with a mode at intermediate frequency — an excess
of alleles at middling frequency, from balancing selection or from strong population structure —
projects onto the nearest monotone shape and the excess is lost. A bottlenecked panel's
flatter-than-neutral spectrum *is* representable, by `α_alt` rising above `θ`. **Whether a
domesticated selfer is where the family gets stretched is still untested.** A fit against an
independently-called VCF of 18 accessions returned `α_alt = 0.81 θ`, slightly below `θ` rather than
above — but an independent-chromosome projection on a *perfectly neutral* panel at `F = 0.9` returns
about `0.83 θ` by itself, so that number is consistent with being nothing but the inbreeding bias the
paragraph above requires the projection to remove. It says nothing about tomato either way until the
fit is redone with two-branch sampling. The richer alternative — carrying the spectrum as a
mixture over allele frequencies — was rejected because §6's cohort term stops being counts added
onto one Dirichlet and becomes a posterior over mixture weights, per sample per locus, and nothing
yet shows two parameters losing anything that moves a genotype (Q4, §11).

**Trap: the census sites are themselves loci this caller genotypes.** The spectrum is fitted from a
scattered set of sites that are later called, so each locus's starting concentration carries a
trace of that locus's own data — the double-counting §6 subtracts the sample's own copies to avoid.
At the
ten thousand census sites the pre-pass targets, any one site contributes about a ten-thousandth of
the spectrum, which is far below anything that moves a genotype. It is written down because it is
the same mechanism, not because it is a defect to fix.

**Trap, and this one has already cost a measurement: exclude a site on depth, never on how
heterozygous it looks.** Collapsed duplications — a stretch the sample carries twice and the
reference once — put about half their reads at odds with the reference in every sample, which is
what a heterozygote looks like, and they inflate the low-frequency classes the projection is most
sensitive to. The tempting repair is to identify those positions and drop them. **A repair that
identifies them by their apparent heterozygosity is selecting on the estimand**, which is the
ascertainment [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §3 forbids
when choosing the positions in the first place, and it fails in exactly the way that rule predicts:
ng's duplicated class takes the tomato panel's median accession from 0.867 heterozygous positions per
kilobase to **0.064** — 93% of them — where the artefact is worth about 11%, because across a cohort
a heterozygote and a duplication carrier give the same read counts and the only signal separating
them is some accession being homozygous for the non-reference allele, which an inbred panel has made
rare ([`../reports/duplicated_class_on_real_reads_2026-08-14.md`](../reports/duplicated_class_on_real_reads_2026-08-14.md)).
**Depth is a covariate rather than the estimand**, roughly independent of allele frequency, so
excluding a census position whose depth runs well above its sample's GC-corrected median removes
duplications without conditioning on what is being measured. It is a far weaker instrument — it
reaches only duplications that show a depth excess — and weak is the correct trade here, because the
strong instrument biases the answer and the weak one does not.

**Report which regime produced the pair.** The strength of the neutral regularizer — how many
sites' worth of pseudo-counts — is a named, overridable parameter, and the run's output must carry
it together with whether the concentration came out prior-dominated or data-dominated. This is the
same complaint the paragraph above makes about `θ`'s fallback, for the same reason: two runs that
used different information are otherwise indistinguishable in their output.

**And a panel-wide ratio is the wrong number to quote as reassurance.** At a regularizer worth ten
sites against 31,084 variable ones, the aggregate ratio is 3,100 to 1 — but in that same panel the
thinnest allele-count class held **2 sites** and was outweighed only **39 to 1**, about 13 to 1 once
scaled down to the ten thousand variable census sites the pre-pass targets. **The tail is where the
regularizer actually binds**, so the reported figure should be per class, not a single ratio.

### 4.2 SNP and indel share one θ, for now

Production ran different pseudocounts for SNP alternatives and for indel alternatives in its
plug-in path — `0.01` against `0.00125`, an 8:1 ratio inherited from GATK's per-class default
([`posterior_engine.rs:114`](../../../../src/var_calling/posterior_engine.rs),
[`:123`](../../../../src/var_calling/posterior_engine.rs)). The genetic argument behind it is real:
indels arise at a different rate from substitutions, so their spectra differ.

**ng uses a single `θ` for both classes in the first version, and the reason is that the pre-pass
measures one.** Its diversity comes from observed heterozygosity summed over a windowed histogram
([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §5), which does not separate
substitutions from short insertions and deletions. Splitting the prior by class before the estimate
is split by class would mean inventing the ratio, and the 8:1 above is inherited from another tool
and never measured here.

**Open (Q1, §11), with the measurement that settles it.**

---

## 5. The starting point on the STR path

**Everything §4 leans on is false at a repeat tract**, and each of the three failures is separate:

- **The reference length carries no presumption, because the spectrum is not the same shape.**
  §4's asymmetry is earned by rarity: at an ordinary site the alternatives are rare, so the one
  chromosome the reference records is almost always the common one. At a repeat tract the
  alternatives are not rare and often no length is clearly the major one, so the reference
  accession's length is one draw among several common ones. `α_ref = 1` against small alternatives
  asserts §4's spectrum at a locus that does not have it.
- **The alleles are not unordered.** A tract of 11 repeats is adjacent to one of 10 and 12, and far
  from one of 4. A Dirichlet that splits mass evenly across alternatives throws that away — and it
  is the only structure that makes a rare long allele believable at all.
- **The diversity is a different number.** Repeat tracts mutate orders of magnitude faster than
  bases do, so the cohort's STR diversity is a separate parameter, which the pre-pass estimates
  separately and requires be emitted so a consumer cannot confuse the two
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3). **Today's production STR path
  uses a fixed `SFS_THETA = 0.01`**, freebayes' default, commented *"Fixed, not a per-run knob"*
  ([`freebayes_emit.rs:42`](../../../../src/ssr/cohort/freebayes_emit.rs)) — a SNP-scale constant
  standing in for a quantity nobody measured. Not repeating that is the point of the separation.

**§4.1's projection does not reach here, and the reason is the same separation.** The spectrum the
pre-pass fits is built from the generic census, which excludes repeat tracts outright
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2), so there is no STR spectrum to
project and the concentration below is not read off one. What the STR path takes from the pre-pass
is a single number — the panel's STR gene diversity — and §5.1 says what it does with it.

### 5.1 The seed: where the mass sits, and how much of it there is

**Two questions, deliberately answered by two different parameters** — this is the one place ng
departs from production's shape rather than porting it:

```text
shape:  the mass falls off geometrically with distance from the cohort's modal repeat count
        weight(allele j) ∝ decay ^ |Δ_j|,   Δ_j = (repeat count of j) − (cohort modal count)

total:  the weights are scaled so that the prior's own implied gene diversity equals the
        cohort's STR gene diversity D from the pre-pass, which for a shape whose Simpson
        index is c means

              Σ α  =  D / (1 − c − D),        c = Σ_j (weight_j / Σ weights)²

        c is the shape's Simpson index: the chance that two copies drawn from the
        shape alone land on the same repeat count

        and not Σ α = D, which is a different quantity in different units
```

The **shape** is production's `G₀`, ported: a geometric decay in unit offset from the locus's modal
length, floored at a tiny positive value so a far allele — a long heterozygous copy the candidate
set nearly missed — keeps a non-zero prior and stays recoverable rather than falling into an
absorbing zero ([`allele_freq_prior.rs:25`](../../../../src/ssr/cohort/allele_freq_prior.rs)). The
decay is fitted per group of loci by the pre-pass.

The **total** is new here, and the reason is a defect production records against itself. `G₀` was
designed as a regulariser for a *plug-in* frequency estimate, where its total mass means smoothing
strength and nothing else. Its size is a by-product of the decay: at the fallback decay of `0.5`
([`param_estimation.rs:167`](../../../../src/ssr/cohort/param_estimation.rs)) the weights over a
mode and a few neighbours on each side sum to between 2 and 3. **Reused unchanged as a Dirichlet
concentration, that number stops being a smoothing knob and becomes a claim about how polymorphic
the locus is** — and how hard the prior will resist the reads. Whether 2.5 is the right claim for a
repeat tract is unmeasured; production carries the question as a deferred re-tune
("`G₀`-as-DM-concentration too-diffuse/tight"). The pre-pass's STR gene diversity is the measured
quantity that answers it, which is why ng sets the total from that and leaves `G₀` to do only the
job it was fitted for.

**Setting `Σα` to `D` itself would be a units error, and it was in this document until
2026-08-19.** Gene diversity is a probability — the chance two copies drawn at random carry
different lengths — while a concentration is a count of chromosomes (§1). What a Dirichlet with
total `A` and Simpson index `c` actually implies is `A(1 − c)/(A + 1)`, so `A = D` asserts
`D(1 − c)/(D + 1)`, which is always less than `D`. **Measured on 1,236 polymorphic tomato STR
loci** at the coded fallback decay of `0.5`: the median locus carries `D = 0.087` and the prior
would assert `0.030` — a paired median ratio of **0.40**, tenth percentile 0.22. The total that
reproduces the measurement is a median of **2.8 × D**, ninetieth percentile 8.5
([`../../../../benchmarks/ssr_tomato1/scripts/g0_total_vs_gene_diversity.py`](../../../../benchmarks/ssr_tomato1/scripts/g0_total_vs_gene_diversity.py)).
So the correction is one-directional, modest in the middle of the distribution and large in its
tail.

**A ceiling the total cannot reach, and this is the harder half.** `A(1 − c)/(A + 1)` rises to
`1 − c` and stops, so a locus whose measured `D` is at or above `1 − c` **has no concentration that
reproduces it at all** — the geometry itself cannot express that much diversity, and rescaling is
not a repair. At the fallback decay that is **119 of those 1,236 loci, about one in ten**; 242 at a
decay of `0.3` and 49 at `0.7`. **The knob that moves it is the decay**, which is fitted for a
different job entirely, so this is a defect in the *shape* rather than in the total and §5.1's
separation of the two questions does not resolve it. **Both figures are, if anything, optimistic:**
the measured `D` comes from called genotypes at about 3 reads a position on a panel whose apparent
`F_IS` is 0.82, and low-depth calling in a selfer books ambiguous sites as homozygous, which
understates allele diversity. Open (Q2, §11).

**One locus the rule does not reach: a tract with a single candidate length.** Its shape has a
Simpson index of exactly 1 and therefore a ceiling of 0, so `D ≥ 1 − c` would refuse every
monomorphic tract whatever the measurement, including a measurement of zero. There is nothing to
refuse — one length is one genotype, whose prior probability is 1 at any positive concentration —
so the builder seeds it at one chromosome and the rule starts at two candidate lengths.

**The one-in-ten figure is a tomato figure, and at the other end of the cohort range the refusal
is the rule rather than the exception.** A single diploid sample shows at most three lengths at a
tract, and the most those shapes can imply at the fallback decay is `0.444` at two lengths and
`0.625` at three — while the pre-pass, which fits this quantity at every cohort size down to one
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3), returns that genome's own
repeat diversity, about `0.72` on the GIAB HG002 benchmark where 72 tandem repeats in 100 are
heterozygous (§5.3). **So one outbred genome is refused at every tract**, and no decay rescues it:
the ceiling saturates at `1 − 5/27 = 0.815` at the fallback decay however many lengths a tract
carries. Whatever Q2 settles has to work at ten refusals in ten, not only at one in ten.

**And a second spelling of one length raises the ceiling**, which gives Q3 a size: an interrupted
repeat sitting on an occupied rung flattens the shape, so a two-length tract's ceiling goes from
`0.444` to `0.625` merely for the cohort having shown the interruption. Whether a locus is refused
therefore depends on how many spellings of one length the panel happened to carry, not only on how
much prior mass it collects (§5.2).

**Separating the two questions is still what makes the STR prior defensible:** the pre-pass
measures how variable repeat tracts are in this panel, and the geometry says where that variability
sits relative to the mode. Neither number is asked to do the other's job — but the geometry has to
be able to hold what the measurement reports, and at one locus in ten it cannot.

### 5.2 Two alleles of the same length

ng stores an STR allele as the **observed sequence, not a repeat count**, because two alleles of
the same length can differ by an interior substitution — an interrupted repeat — and a count cannot
tell them apart ([`locus_generation_ssr.md`](locus_generation_ssr.md) §3). The seed above is
indexed by repeat count, so two such alleles land on the same rung.

Production computes each candidate's offset from `sequence length / period`
([`allele_freq_prior.rs:25`](../../../../src/ssr/cohort/allele_freq_prior.rs)), so both alleles
receive the rung's full weight independently — a locus with an interruption variant collects more
total prior mass than one without, purely for having it. **Leaning: divide the rung's weight among
the alleles that sit on it**, so the total is a property of the locus and not of how many spellings
of one length the cohort happened to show. Open (Q3, §11), and it needs the interrupted-repeat work
to say how the division should be weighted.

### 5.3 What is and is not measured here

**The marginalized prior is adopted on the STR path by argument, not yet by measurement, and that
should be visible to whoever benchmarks it.** The GIAB result in §2.2 is a SNP result. Production's
STR comparison — tomato, 51 samples, against HipSTR — was **inconclusive rather than negative**:
the marginalized prior emitted about 30% fewer of the loci HipSTR called polymorphic (774 → 539) at
flat concordance (96.5% → 96.1%), but tomato is an extreme selfer (`F_IS ≈ 0.82`) and HipSTR is
blind to `F`, so "recall against HipSTR" penalises any `F`-aware prior and cannot separate a
conservative prior from a wrong one. Production therefore ships its marginalized STR prior behind
an environment toggle with the plug-in prior still the default
([`driver.rs:287`](../../../../src/ssr/cohort/driver.rs)). Report:
`doc/devel/reports/ssr_marginalized_prior_benchmark_2026-07-07.md`.

**The single-sample half of that question has since been answered, and the answer is a zero.**
`benchmarks/ssr_hg002/` holds the GIAB HG002 tandem-repeat benchmark v1.0.1 on GRCh38 — 36,497
phased genotypes from an **assembly** rather than from another caller, about **72% of them
heterozygous**, so the sample is outbred, §3.2's inbreeding branch is inert, and the comparison
isolates the frequency prior. Re-calling the same pileups under both priors at 50, 30, 20 and 15×
gives, at every locus both priors emitted, **the identical genotype — 0 differences out of the 732 to 1,679 loci
both emitted, depending on coverage.** The only effect is that the marginalized prior withholds about **1 locus
in 100**, and those are close to a coin flip: at 30× it drops 13 correct genotypes and 12 wrong
ones. Report:
[`../../reports/ssr_prior_hg002_single_sample_2026-08-18.md`](../../reports/ssr_prior_hg002_single_sample_2026-08-18.md).

**The reason it is a zero is §2.3's trap seen from the other side.** The SNP failure came from a
*reference-privileged* starting concentration, and the STR path never had one — `G₀` is centred on
the cohort's modal repeat count, which at one sample is that sample's own mode. There is no
reference allele being favoured, so there is nothing for marginalizing to undo. **What this buys
ng is that adopting the marginalized prior on the STR path costs nothing at the single-sample
end**, and it relocates the tomato result above: since one sample is neutral, that 30% is a
*cohort* effect — the leave-one-out term doing something — and has to be judged in a cohort rather
than blamed on the prior's form.

---

## 6. What varies per locus, what does not, and what happens at one sample

**The frequency spectrum is a property of the genome, not of a locus.** §4 and §5 set one starting
concentration for the whole run. What varies locus by locus is **what the other samples showed
there**, and it enters as counts added onto that starting point:

```text
for sample s at this locus:
    α'_s(a) = α_seed(a) + max(0, E[copies of a across the cohort] − E[copies of a in s])
```

The expectation is over the current genotype posteriors — no genotype is called to produce it,
which is what lets it be used at low coverage. The `max(0, …)` guards floating-point noise only:
the sample's own count is one non-negative addend of the total, so the true difference cannot be
negative ([`em.rs:278`](../../../../src/ssr/cohort/em.rs),
[`posterior_engine.rs:3200`](../../../../src/var_calling/posterior_engine.rs)).

**Subtracting the sample's own contribution is not a refinement, it is what makes the prior a
prior.** Without it, a sample's own reads would arrive twice — once through the likelihood and once
through the frequency they helped estimate. With one sample and no subtraction, a genuinely
homozygous-variant sample could only push the frequency estimate to a diluted value, and would then
be told that value made it heterozygous.

**It is not, however, the mechanism behind §2.2's 214 sites, and this document said so twice before
2026-08-19.** §2.3 has the counterfactual: production's plug-in ran `α_ref = 10`, and *with* the
subtraction in place its estimate is still about 1 in 10,000, which puts the heterozygous prior far
higher than 22:1 rather than repairing it — while at `α_ref = 1` *without* the subtraction the
failure disappears. Double counting is neither necessary nor sufficient for that measurement; the
starting concentration is. Subtracting the sample's own contribution is required because using a
sample's reads twice is wrong, which needs no measurement to justify.

**At one sample the cohort term is exactly zero**, and no branch is needed to make it so: the
cohort total and the sample's own count are the same number. What remains is the starting
concentration bent by that sample's `F` — a prior identical at every locus of the same allele
count and class, with every site-to-site distinction coming from the reads. That is the correct
answer, not a degraded one: **a single genome carries no information about how common an allele is
at a particular locus.** It does carry information about how variable the genome is *on average* —
that is what `θ` is, and §4.1 fits it from one sample as readily as from a thousand. The two are
different questions, and only the first needs other samples.

**At several thousand samples the cohort term swamps the starting concentration** and the prior
converges on the panel's own frequencies; there, the marginalized and plug-in priors agree, because
`Var(p)` has gone to zero (§2.2). Both ends of the committed range are one formula.

**Cost note for the EM document.** Because `α'_s` differs per sample, the prior cannot be computed
once per locus and shared, which is how production's homogeneous-`F` fast path works
([`posterior_engine.rs:3799`](../../../../src/var_calling/posterior_engine.rs)); its
leave-one-out path runs per sample and carries a deferred perf item for a large-cohort shared
approximation. That trade belongs to the EM document, not here.

---

## 7. The inbreeding coefficient: where it comes from, what it means

**Per sample, fitted by the pre-pass, frozen before calling, never iterated by the EM.** The prior
reads it and nothing writes it.

Four things about the number that a reader will otherwise get wrong:

- **It is genome-wide.** One value per sample, applied at every locus. A locus in a region of
  unusual ancestry gets the sample's average. Accepted, not modelled.
- **What the prior needs is the homozygote excess measured against the frequencies it actually
  uses, and those are pooled over the whole panel.** §6 builds every concentration from the cohort's
  expected allele copies, with no notion of subpopulation. A panel drawn from several subpopulations
  therefore over-predicts heterozygotes by exactly twice the variance of the allele frequency across
  those subpopulations — the Wahlund effect — and that distortion is algebraically identical to
  inbreeding:

  ```text
  P(AA) = p̄² + Var(p)      = p̄² + F_ST·p̄(1 − p̄)
  P(Aa) = 2p̄(1 − p̄) − 2Var(p) = 2p̄(1 − p̄)·(1 − F_ST)
  ```

  which is Wright's form with `F = F_ST`. So the coefficient this mixture wants is the **total**
  deficit against the pooled panel, and not autozygosity alone. **In Wright's hierarchy that is
  `F_IT`, the individual against the total, not `F_IS`**, which is the individual against its own
  subpopulation — a distinction with teeth here, because `F_IS` is roughly what runs of homozygosity
  already deliver, so calling the wanted quantity `F_IS` would say the substitution below costs
  nothing.
  **It should not be reported to a user as a pedigree statement about the accession**, and whatever
  writes it into the output owes that caveat. *(The `Var(p)` here is real variation across
  subpopulations. §2.2's `Var(p)` is uncertainty about a frequency we have not pinned down: the same
  algebra from a different source, and the two stay separately identified because the second shrinks
  as the cohort grows and the first does not.)*
- **What the pre-pass supplies is `F_autozygosity`, which under-corrects, and the reason is
  estimation rather than meaning.** The two quantities carry separate names
  ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §5): `F_autozygosity` is
  `F_ROH`, read off runs of homozygosity in the windowed accumulator, and `F_hom_excess` is
  `1 − Hobs/Hexp` against the fitted spectrum — **`F_IS` on an unstructured panel and `F_IT` the
  moment structure exists**, since its `Hexp` comes from the pooled frequencies. That sibling's table
  names it `F_IS` unconditionally, which is the mislabel this bullet exists to correct. The bullet
  above asks for the second and the caller is
  handed the first, because the second is `1 − Hobs/Hexp` and so passes any bias in observed
  heterozygosity through at full size — false heterozygotes from collapsed paralogs at five times
  tomato's rate of one per kilobase inflate it eight-fold while leaving the runs estimate unmoved —
  and because it needs a panel, so it does not exist at one sample. **What the substitution costs is
  the structure component of the deficit.** On a selfing panel most of the deficit really is
  autozygosity and the two run together, so the gap is small; on an outbred structured panel it would
  not be, and this prior would over-predict heterozygotes there. Open (Q6, §11).

**How much the substitution costs, and why it is nonetheless the right trade.** Wright's
hierarchical relation gives `1 − F_IT = (1 − F_IS)(1 − F_ST)`, where the coefficient measured
against pooled panel frequencies is `F_IT`, what runs of homozygosity read is `F_IS`, and the
structure is `F_ST`. The prior multiplies its heterozygote branch by `(1 − F)`, so substituting one
for the other over-predicts heterozygotes by

```text
(1 − F_autozygosity) / (1 − F_IT)  =  1 / (1 − F_ST)
```

**and the sample's own inbreeding cancels out of that ratio entirely.** At `F_ST = 0.20` — the
strongest divergence the pre-pass's own contamination harness simulates
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.2) — it is 25% too much
heterozygote prior, about 1 on the Phred scale, against the 60 Phred that two clean alternative
bases at Q30 supply. **The error the other choice risks is not bounded like that.** Where `F` is
near 1 the whole of `1 − F` is `Hobs/Hexp`, so the eight-fold inflation of observed heterozygosity
that collapsed paralogs produce moves the heterozygote branch eight-fold with it — and at one sample
there is no estimate to move. One horn costs about a Phred on a structured panel; the other can be
an order of magnitude wrong on data we hold, and undefined across the input range this caller
commits to. That asymmetry is the justification, and it is why the choice does not turn on which
quantity is conceptually the better match.

**`1/(1 − F_ST)` is a ceiling rather than an equality, and the `F_IS ≈ F_autozygosity` step leaks in
both directions.** Downward: 100 kb windows resolve runs of about 300 kb and longer, so autozygosity
old enough to have been broken into shorter tracts is invisible to the runs estimator while still
suppressing heterozygosity, which makes the true factor `1 / [(1 − F_ST)·(1 − F_old)]` — worse than
the ceiling, and unquantified here. Upward: runs of homozygosity also capture identity by descent
generated by a subpopulation's own recent drift, which is already part of `F_ST`, so `F_autozygosity`
can sit *above* `F_IS` and the real over-prediction falls between 1 and `1/(1 − F_ST)`. On a selfing
crop `F_old` is captured, because selfing regenerates long tracts every generation; on an old
bottlenecked outbred population it is the term with no number.
- **It must be estimable without a population expectation.** The pre-pass's ratio estimator
  `F = 1 − Hobs/Hexp` needs an expected heterozygosity, which is itself computed from `F`; the
  cohort gather states plainly that feeding one into the other is circular and that the
  runs-of-homozygosity estimator is what makes the pair estimable
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3). This prior inherits that
  constraint rather than restating it: **it consumes whichever estimator the pre-pass settles on,
  and its correctness depends on that choice being the non-circular one.**

**`F = 1` makes every heterozygote impossible, and the mixture survives that limit** — the
heterozygotes fall to the probability floor and the two homozygotes stand in the ratio
`α_ref : α_alt`, which production pins in a test
(`dirichlet_prior_full_inbreeding_concentrates_on_homozygotes`,
[`posterior_engine.rs:4341`](../../../../src/var_calling/posterior_engine.rs)). **What is capped is
the estimate, not the model:** production's inbreeding estimator clamps at `0.99`
([`inbreeding.rs:25`](../../../../src/paralog/inbreeding.rs)), so no sample whose `F` was **fitted**
reaches the caller carrying a prior that has ruled heterozygotes out.

**A sample whose `F` was supplied does.** Production's second door is `--inbreeding-coefficient`,
whose parser admits the closed `[0, 1]` and is pinned admitting `1.0`
([`parsers.rs:166`](../../../../src/pop_var_caller/cli/parsers.rs), test at
[`:392`](../../../../src/pop_var_caller/cli/parsers.rs)); the value goes to the engine as typed
([`pipeline.rs:343`](../../../../src/var_calling/pipeline.rs)). The cap is a line inside one
estimator, and the other route round it is a command-line flag. *(Corrected 2026-08-23 while
building `InbreedingF`'s half-open check — the earlier text claimed the guarantee for every
sample.)*

**ng has no such second door today, and the point of the newtype is to keep it that way.** Every
construction of `InbreedingF` from a raw number outside test code is the fitted path
([`runs.rs:634`](../../../../src/ng/parameter_estimation/generic/runs.rs)); a supplied coefficient
arrives already built, so the checked constructor is the only way in. **The exposure is
prospective:** whoever eventually gives ng a command line will reach for production's parser, whose
range is closed, and a flag that admitted `1.0` would hand the prior a value no fit can produce.
The half-open type is what makes that a compile-and-check problem rather than a silent one.

ng should carry that ceiling in a validated newtype — `InbreedingF` in `[0, 1)`, per the interface
conventions ([`ng_step_interfaces.md`](../arch/ng_step_interfaces.md)) — so it is a property of the
type rather than a line inside one estimator that a second estimator can forget. **The trap for the
coder:** production's cap lives in the estimator and its engine config accepts `1.0`
([`posterior_engine.rs:4348`](../../../../src/var_calling/posterior_engine.rs) sets it), so porting
the engine without the newtype ports a gap.

**What the newtype's range buys, stated exactly, because it is easy to over-read.** Excluding the
endpoint removes the mathematical limit and nothing more: it keeps `ln(1 − F)` finite. It is not a
numerically meaningful cap — the largest `f64` below one leaves `1 − F = 2⁻⁵³`, about **160 on the
Phred scale** against every heterozygote, where two clean alternative bases at Q30 supply **60**.
Production's `0.99` is 20 Phred, which evidence can overcome. So the type makes `F = 1`
unrepresentable and **every estimator still owes its own cap**; ng's fitted path clamps at `0.99`
for exactly that reason ([`calling_prerequisites.md`](../impl_plan/calling_prerequisites.md) A2).

---

## 8. Cross-cutting concerns

**Numerics.** Everything in log space. Every concentration is strictly positive — the alternative
concentrations are floored so `lgamma` stays finite when the estimated diversity is exactly zero, a
fully invariant cohort or an explicit `--diversity 0`
([`genetics.rs:187`](../../../../src/genetics.rs)). Probabilities that reach zero before a
logarithm are floored rather than allowed to become `−∞`
([`genetics.rs:18`](../../../../src/genetics.rs)), so a zero-probability genotype yields a finite,
very negative log-prior and a whole sample's row cannot become `NaN` on one impossible genotype.

**Determinism.** The prior is a pure function of (frozen parameters, expected allele copies,
allele table). No RNG, no clock, no thread-dependent iteration. The one place order matters is the
sum of expected copies across the cohort, which must be accumulated in a fixed sample order — that
is the EM document's contract, and it is the same requirement the merge already carries for
byte-identical output at any worker count ([`run_streaming.md`](run_streaming.md) §12, item 1).

**Cost and memory.** Per sample per locus: one `lgamma` per (allele, non-zero copy count) pair,
plus one `logsumexp` for each homozygous genotype. **Nothing may allocate inside the per-sample
loop** — the caller hands in scratch sized by allele count and genotype count, and the prior fills
it. Production lifted exactly these buffers out of its EM iteration for a measured reason: a
profile put the allocator's own self-time at about 16% of cycles before the lift
([`posterior_engine.rs:1874`](../../../../src/var_calling/posterior_engine.rs), and the profile
report it cites). Memory is `O(alleles + genotypes)` of scratch per worker, independent of cohort
size.

**Errors.** The prior has no failure mode of its own that is not a caller bug: a non-positive
concentration, a mis-shaped allele table, a genotype table that disagrees with the allele count.
These are assertions, not `Result`s, and the structural ones must hold in release — production's
primitive asserts them because a short coefficient array would otherwise let the iteration
**silently truncate** and corrupt every downstream genotype index without panicking
([`genetics.rs:127`](../../../../src/genetics.rs)).

---

## 9. Reuse map

**The mathematics of this document already exists in production and is already shared between its
two callers.** ng ports it; nothing here needs inventing.

| what | production code | how ng reuses it |
|---|---|---|
| the Dirichlet-multinomial log-priors | [`genetics.rs:127`](../../../../src/genetics.rs) | **ported as-is** — it already takes flat arrays specifically to avoid a back-reference into the engine |
| the SNP/indel starting concentration | `alpha_from_diversity`, [`genetics.rs:214`](../../../../src/genetics.rs) | **shape ported, source not** — production hard-codes `(1, θ)`; §4.1 reads the pair off the fitted spectrum instead, and `(1, θ)` is where that lands on a neutral panel |
| the projection of a spectrum onto `(α_ref, α_alt)` | — | **new** (§4.1); production has nothing to port, because it never fitted a spectrum |
| the alt floor | [`genetics.rs:187`](../../../../src/genetics.rs) | ported with its reasoning |
| the inbreeding mixture | [`posterior_engine.rs:3799`](../../../../src/var_calling/posterior_engine.rs), [`:3217`](../../../../src/var_calling/posterior_engine.rs) | ported as §3.2's two-branch form |
| leave-one-out concentration | [`em.rs:278`](../../../../src/ssr/cohort/em.rs) | ported; §6 is its contract |
| the STR geometric seed | [`allele_freq_prior.rs:25`](../../../../src/ssr/cohort/allele_freq_prior.rs) | **shape ported, total mass not** — §5.1 |
| the biallelic Wright formulas | [`genetics.rs:66`](../../../../src/genetics.rs) | **test oracle only**, not a code path (§3.2) |

**Parity oracle.** Production's own test module computes each Dirichlet-multinomial log-prior a
second way, from the rising factorial with no `lgamma` at all
([`genetics.rs`](../../../../src/genetics.rs) tests, `pochhammer_ln` / `dm_log_prior_oracle`). ng's
port should carry that oracle across: it is an independent implementation, not a golden value, so
it keeps checking after the constants move.

---

## 10. Deferred, with a recommended home

- **A structure-aware prior — per-sample allele frequencies instead of panel-pooled ones.** This is
  the proper repair for the Wahlund term §7 currently absorbs into `F`: give each sample the
  frequencies its own ancestry implies and there is no excess left for a coefficient to correct.
  **The machinery is already specified upstream, for another purpose.** Contamination estimation
  needs exactly this and has it — each individual's frequency at a marker as a **linear function of
  that individual's coordinates in a principal-component space**, fitted across the whole panel and
  maximised jointly with everything else
  ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.2, following
  `verifyBamID2`, Zhang et al., *Genome Research* 2020). One measured warning this prior would
  inherit: **borrow strength across the panel, never partition it** — per-subpopulation frequencies
  fitted from twelve samples each performed worse than pooled ones on that harness, adding about
  0.015 to every sample's contamination estimate. **That is a defect of partitioning, not of sample
  size**, so a larger cohort is not what makes partitioning safe; the smooth fit is.

  **The shape it would take here.** §6's cohort term stops being an unweighted leave-one-out sum and
  becomes a fitted frequency: at each locus, regress expected allele dosage on the panel's first `K`
  components across every sample — `K+1` coefficients, one small fit per locus per EM iteration —
  then predict this sample's own frequency with itself held out, which a linear fit gives in closed
  form rather than by refitting. The concentration is then `α_seed(a) + (effective cohort size) ·
  f_s(a)`. Two properties are easy to lose and must not be: `Σα'` still has to state how much
  evidence the cohort supplies, so the weights need normalising to an effective sample count rather
  than to whatever the regression returns; and `f_s` must be held in `[0, 1]`, since a linear
  predictor leaves the interval at loci with strong loadings.

  **Two gates, and neither is a round number chosen by taste.** *Enough samples:* the fit spends
  `K+1` parameters per locus, so the crossover against the pooled frequency sits near twenty samples
  per parameter — about 100 accessions at `verifyBamID2`'s default of four components. *Enough
  structure:* gate on **how many components are significant under the Tracy–Widom criterion**
  (Patterson, Price & Reich 2006), not on the proportion of variance they explain — that proportion
  moves with marker count, sample count and read depth, so two panels with identical structure give
  different numbers. Report the gate itself as the `F_ST` those components imply, because §7 states
  the gain as `1/(1 − F_ST)` and a threshold in the same units as the benefit needs no calibration:
  at `F_ST = 0.05` the refinement buys a fifth of a Phred, at 0.15 about 0.7.

  **What adopting it would fix beyond its own size.** Once each sample is judged against its own
  ancestry's frequencies, the Wahlund excess is out of the frequencies, so the coefficient §3.2's
  mixture needs becomes the within-subpopulation one — which is what `F_autozygosity` estimates. The
  approximation §7 accepts stops being an approximation, and Q6's structure half retires.

  **One leak to state before anyone builds it.** A sample's principal-component coordinates are
  estimated from data including its own genotypes, so holding it out of the per-locus regression does
  not fully hold it out of its own frequency. It is the double-counting mechanism of §2.2 again, small
  at a hundred samples and not to be discovered later.

  **Deferred, and the size is why.** The gain to genotyping is bounded by that same `1/(1 − F_ST)` —
  about 1 Phred at the strongest divergence the contamination harness simulates — while the measured
  payoff of individual-specific frequencies was on contamination, where a 3% contaminated sample came
  back at 0.5%. **Home:** here, as a by-product once contamination produces those frequencies and
  settles the component count (that document's open question 4). Not worth scheduling ahead of
  anything measured in whole Phreds.
- **A relatedness-aware prior.** The pre-pass estimates relatedness between sample pairs
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §6) and this prior treats every
  sample as an independent draw. Using it means the leave-one-out term is no longer a plain sum
  over the other samples. It is the sibling of the item above — one corrects *which* frequencies a
  sample is judged against, the other corrects *double-counting* between samples that share
  ancestry. **Home:** a spec of its own, after the three calling documents land.
- **Partial identity by descent above diploidy** (§3.3). **Home:** already claimed — the pre-pass's
  deferred *population genetics above diploidy* spec ([`parameter_prepass.md`](parameter_prepass.md)
  §8). This document's only obligation is to keep the homozygous test in one function.
- **A shared-prior fast path for large cohorts**, where the leave-one-out correction is small
  enough that one prior serves every sample. **Home:** the EM loop document, as a perf item with a
  measured threshold.
- **A prior carrying the spectrum as a mixture over allele frequencies**, rather than projected onto
  two concentrations (§4.1). It is the form that would hold a mode at intermediate frequency, and it
  costs §6 its closed form. **Home:** here, if Q4b's comparison shows the projection moving a
  genotype.

---

## 11. Open questions

**Q1 — one `θ` for SNPs and indels, or two?** *Leaning: one, for the first version.* The pre-pass
measures one number (§4.2), and the 8:1 SNP-to-indel ratio production carries is inherited from
GATK and never measured here. **Settled by:** whether the pre-pass's windowed histogram can be
split by class without a second traversal; if it can, fit both and compare indel false-positive
rate at matched recall on GIAB, where the truth set distinguishes the classes. **Confirm before
code** only in the sense that the concentration function must take the class as an argument even
if both classes pass the same value — otherwise splitting later touches every call site.

**Which function that is, settled 2026-08-23: the per-locus expansion, not the projection.** The
projection reads the *shape* of variation off the panel's allele counts, which the pre-pass fits
without separating the two classes; a class-specific *scale* belongs where the run's total is
shared out over a locus's alleles. When the estimate splits, the run holds one seed carrying both
totals and the per-locus expansion picks between them; carrying the argument at both ends would
apply the ratio twice. **What is still open under this question is a locus carrying one
alternative of each kind** — a substitution and a short indel at one position, which one class per
locus cannot express and which nothing in the caller can currently tell apart.

**Q2 — reopened 2026-08-19, and it is now a question about the shape rather than the total.** As
posed it asked whether the STR total concentration is the pre-pass's STR gene diversity. It is not,
in those words: gene diversity is a probability and a concentration is a count of chromosomes, and
`Σα = D` makes the prior assert about two-fifths of what was measured on tomato. §5.1 carries the
mapping that fixes it, `Σα = D/(1 − c − D)`, and that part is settled.

**What is open is the locus the geometry cannot hold at any total.** The prior's implied diversity
is bounded by `1 − c`, and on 1,236 polymorphic tomato loci **119 — about one in ten — measure at or
above that bound** at the coded fallback decay, 242 at `0.3` and 49 at `0.7`. **At one outbred
sample it is every locus** (§5.1), so the candidate answers below have to be judged at ten in ten
and not only at one in ten — which rules out any of them that is affordable only because it is
rare.

Three candidate answers, none measured: fit the decay against gene diversity as well as against
stutter, so the shape has to be able to hold what the panel shows; carry a floor of prior mass on the alleles the
geometry starves, which is the `G₀` floor doing a second job; or let the seed refuse such loci to
the reads entirely, which is honest and gives up the regulariser exactly where the locus is most
polymorphic. *Leaning: none — the choice needs the number of affected loci on a second panel, since
one in ten on a selfing crop at 3 reads a position may be a tomato figure rather than a caller
figure.* **Settled by:** re-running [`../../../../benchmarks/ssr_tomato1/scripts/g0_total_vs_gene_diversity.py`](../../../../benchmarks/ssr_tomato1/scripts/g0_total_vs_gene_diversity.py)
on the GIAB HG002 bundle, where depth is high and the genotypes are not a selfer's, and then the
surviving option against truth on Q5's benchmark.

**Q3 — how is a rung's weight divided between two alleles of the same length?** *Leaning: divide
it, so a locus's total prior mass does not grow with the number of spellings the cohort showed*
(§5.2). **Settled by:** the interrupted-repeat work, which has to say whether an interruption is a
separate allele for genotyping or a spelling of one; the branch that explored this
(`ssr-interruptions`) reached a two-phase spec and no code.

**Q4 — closed as posed, 2026-08-19, and replaced by two narrower questions.** It asked whether the
starting concentration should come from the pre-pass's measured spectrum *instead of* the neutral
`1/p` shape. It should not, because there is no *instead*: §4.1's correspondence makes the neutral
shape the value the measured spectrum takes on a neutral panel — to within `θ · H(2N)`, 3 in a
thousand at tomato's diversity — so reading the concentration off the spectrum subsumes the neutral
setting rather than replacing it. What survives:

- **Q4a — how strong is the neutral regularizer on the spectrum?** How many sites' worth of
  pseudo-counts hold the estimate at `θ/k` before the census sites move it. *Leaning: weak enough
  that the census sites dominate in aggregate* — 3,100 to 1 at a regularizer worth ten sites against
  31,084 variable ones — *but the aggregate is the wrong figure to lean on*: the thinnest allele-count
  class in that same panel was outweighed only 39 to 1, and the tail is where the regularizer binds
  (§4.1). **Settled by:** sweeping it on the tomato panel and watching where the fitted pair stops
  moving, with the ratio reported per class.
- **Q4b — does the two-parameter projection lose anything that moves a genotype?** §4.1 says what
  the family cannot hold — a spectrum with a mode at intermediate frequency. **Attempted 2026-08-19,
  held open, and the obstacle is not the one it looked like.** Fitted on tomato1, an independently
  called VCF of 18 accessions gives `α_ref = 1.01`, inside the band where the projection changes
  nothing; our own caller's VCF over the same accessions gives 1.72, and over all 26 gives 4.00 —
  which by §2.3's ratio would move the prior from 2:1 to 8:1, the direction of the historical bug,
  arrived at by fitting.

  **Neither number is a fact about tomato's frequency spectrum, and re-running the caller would not
  make one.** The inputs are two unfiltered VCFs from different sample sets — 26 accessions against
  18, the second missing eight samples — both predating the hidden-paralog filter by five weeks, and
  neither is the census-site fitting §4.1 specifies. But the deeper obstacle survives fixing all of
  that: **the spectrum is fitted by marginalizing over unknown genotypes from the observations, so
  filtering the caller's output changes nothing about it.** Removing the artefact means excluding the
  affected positions from the fit, and the tool built to identify them selects on apparent
  heterozygosity and so removes 93% of the panel's real heterozygosity to clear an 11% artefact
  (§4.1's second trap). **Leave the duplications in and the low-frequency classes are inflated; take
  them out with what exists and the fit is gutted.** Neither reading is the panel's spectrum.

  **Settled by:** the census-site spectrum once it exists, with duplication-suspect positions
  excluded **on depth** rather than on apparent heterozygosity (§4.1). Until then, no claim that the
  projection is neutral has been tested, and none should be made.

**Q5 — does the marginalized prior genotype STR loci better than the plug-in one?** **Half
resolved, 2026-08-18; the cohort half is open.**

**Resolved: at one sample it makes no difference.** On GIAB HG002 at 50, 30, 20 and 15×, every
locus both priors emitted got the identical genotype — 0 differences out of the 732 to 1,679
loci both emitted, depending on coverage — and the marginalized prior withheld about 1 locus in 100, split evenly between calls the
plug-in had right and calls it had wrong (§5.3;
[report](../../reports/ssr_prior_hg002_single_sample_2026-08-18.md)). **The reasoning that closes
it:** the marginalized prior's advantage on the SNP path came from replacing a reference-privileged
starting concentration (§2.3), and the STR seed is mode-centred, so there was never that pull to
undo. **Consequence for this document: nothing changes.** The marginalized prior is kept for both
paths — it is free here, it is the correct treatment of an unknown frequency, and the single-sample
end is where it could have cost something and does not.

**Open: does it help in a cohort?** At one sample the leave-one-out term of §6 is zero by
construction, so the run above could not test the mechanism the prior exists for. The only
multi-sample evidence is tomato — 51 samples, 30% fewer loci emitted than HipSTR called
polymorphic — and now that the single-sample end is known to be neutral, that 30% is a cohort
effect to be explained rather than a property of the prior's form. **Settled by:** a multi-sample
panel with truth that is not another caller. We do not have one; tomato's truth is HipSTR, and
HG002 is one sample. **This is the manager's item** — it is blocked on data, not on work, and until
that data exists the honest statement in any report is that the cohort behaviour of the STR prior
is unmeasured.

**A limitation of the instrument, worth carrying.** Production's STR caller detects 64% of
truth-variant HG002 loci at 50× and 28% at 15×, and essentially none at 5×. The loci it reaches are
the well-covered ones, where the read likelihood is sharp and any prior is a small term — so this
benchmark cannot reach the low-depth single-sample regime where the SNP effect was measured. A
caller with better low-depth recall would make the same question askable again.

**Q6 — how much does the prior lose by taking `F_autozygosity` where its mixture wants the total
homozygote excess against the panel's pooled frequencies?** **Half answered by §7's algebra:** the
structure component costs a factor `1/(1 − F_ST)` on the heterozygote prior, about 1 Phred at
`F_ST = 0.20`, and the sample's own inbreeding cancels out — so this horn is bounded and small
wherever `F_ST` is, while the alternative estimator is not bounded at all. **The decision is
settled; what stays open is one term.** `F_autozygosity` stands in for `F_IS`, and runs shorter than
about 300 kb are invisible to it, so old autozygosity suppresses heterozygosity without being
counted. **Settled by:** a panel with old inbreeding and little recent inbreeding — the case where
the runs estimator and the ratio should disagree most and neither tomato nor GIAB provides one.
**Manager's item, blocked on data rather than on work.** Until it exists the honest statement is
that the inbreeding term is bounded on structure and unquantified on old inbreeding.

---

## 12. How we know it works

**Unit tests, each pinning a property rather than a value.** The first four are properties
production already tests and ng's port should carry across:

1. **The 2:1 ratio holds across diversity.** At `F = 0`, biallelic diploid, the het:hom-alt prior
   ratio stays at 2:1 for every realistic `θ`. This is the §2.3 trap's tripwire: it fails the
   moment someone raises `α_ref`.
2. **The invariant mass tracks `θ`.** Hom-ref prior ≈ `1 − 3θ/2`; raising `θ` moves mass onto the
   variant genotypes.
3. **The full-inbreeding limit is well-behaved.** At `F = 1` every heterozygote sits at the
   probability floor and the two homozygotes stand in the ratio `α_ref : α_alt`. The estimator
   never delivers `F = 1` (§7), so this tests the mathematics at its edge, not a case the caller
   meets.
4. **The independent oracle.** Every log-prior matches a rising-factorial computation that uses no
   `lgamma` (§9).

Seven more are ng's own, and each pins a claim this document makes:

5. **A neutral spectrum projects to `(1, θ)`.** Build the target as the **exact expected spectrum**
   of `Dirichlet(1, θ)` at the panel size under test, in closed form. **Not** by writing the counts
   `θ/k` — that is the small-`θ` approximation (§4.1), which would put the test's own error at about
   `θ · H(2N)`, 3 in a thousand at tomato's diversity, and no honest tolerance could tell that from a
   wiring bug. **And not by drawing sites at random either:** Monte-Carlo noise falls as one over the
   square root of the site count, so no tolerance a fit could be held to would need a number of sites
   anybody can generate, and loosening it to fit destroys what the test is for.

   **What the projection must return is the pair to within the resolution the search was asked
   for, and closer when it is asked for more** — the second half is what separates a resolution
   from a bias, and floating-point equality is not available from any bounded search. Measured on
   the shipped one, over six panel sizes from one individual to 150, three diversities and two
   inbreeding coefficients: 0.25% on `α_ref` and 0.31% on `α_alt` at the shipped 1% resolution, and
   2 parts in 100,000 at a thousand-fold finer.
6. **The projection is invariant to inbreeding.** Build the exact expected spectrum of one population
   frequency density at `F = 0`, `0.6` and `0.9` under §3.2's two-branch sampling; the projection
   must return the same pair from all three. An independent-chromosome projection does not — `α_ref`
   falls to about 0.91 and 0.86 — so this test is what holds §4.1's two-branch requirement in place
   rather than leaving it as prose.
7. **One sample: the projection is the neutral shape.** With `n = 1` no census site is variable
   across the panel, so the fitted spectrum is its regularizer and the pair is `(1, θ)` — the same
   two numbers §4 sets, reached without any test of `n` — the only branch permitted is on the spectrum being absent.
8. **One sample: the cohort term is exactly zero.** With `n = 1` the leave-one-out concentration
   equals the starting concentration, bit for bit, at every locus — no tolerance, no branch.
9. **Monotone in cohort evidence.** Raising the expected copies of an allele across the cohort
   cannot lower its prior weight for a sample that did not contribute the rise.
10. **The STR seed implies the diversity the pre-pass measured — the prior's own `A(1 − c)/(A + 1)`
    recovers `D` to floating-point tolerance**, whatever the decay and however many alleles the locus
    carries (§5.1). **Not** that the concentration sums to `D`, which is the units error §5.1
    records; a test written that way passes on a prior asserting two-fifths of the measurement.
11. **A locus the geometry cannot hold is refused, not silently rescaled.** Where the measured `D` is
    at or above `1 − c` no total reproduces it, so the seed builder must say so — one locus in ten on
    tomato at the fallback decay, and **every locus on one outbred genome** (§5.1). What it does
    instead is Q2's to settle; what it must not do is return the closest total it can reach as though
    it had met the target. **The rule starts at two candidate lengths:** a single-length tract has a
    ceiling of exactly 0 and would be refused at any measurement, and there is nothing there to
    refuse.

**The end-to-end check, and the definition of done for the manager:** the GIAB single-sample 5×
regression of §2.2 — genotype accuracy at true variants and the count of true homozygous-variant
sites called heterozygous. Those two numbers are what the prior is for, and a change to this module
that does not move them has not been tested by anything that matters.

**The STR path's equivalent is the HG002 tandem-repeat bundle**, scored on genotype accuracy given
detection (`benchmarks/ssr_hg002/src/prior_genotype_accuracy.py`). It has been run against the two
priors and cannot separate them at one sample (§5.3, Q5), so **it is a regression guard rather than
a discriminating test**: it will catch a change that breaks STR genotyping, and it will not tell
anyone whether the prior is the right one in a cohort. That gap is the state of the evidence, not
an oversight to be papered over.
