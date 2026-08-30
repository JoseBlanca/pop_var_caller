# ng — the ordinary-site prior's two numbers, measured on the census

**Status:** spec, 2026-08-27. Follows the measurements in
[`../../reports/ng_ordinary_site_prior_moments_2026-08-27.md`](../../reports/ng_ordinary_site_prior_moments_2026-08-27.md),
which answers the six questions in
[`../research/ordinary_site_prior_moments.md`](../research/ordinary_site_prior_moments.md). **No
code yet.**

**What it changes:** where the SNP/indel genotype prior's two numbers come from.
**What it does not change:** the prior itself, the joint fit's model, the repeat-tract path, and
anything about how a genotype is called.

**Companions:** [`ordinary_site_seed.md`](ordinary_site_seed.md) (the seam this edits),
[`population_diversity.md`](population_diversity.md) §3 (where the numbers come from today),
[`calling_priors.md`](calling_priors.md) §4 (what the caller does with them),
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §3 (how census positions
are chosen).

---

## 1. What this is, in one paragraph

The prior for an ordinary site holds two numbers: a concentration for the reference allele and one
for the alternatives. Written the other way round they are **an expected allele frequency and a
total conviction**, and the pair is fixed by two properties of the population — its **mean
alternative-allele frequency** and its **heterozygosity**, how often two copies of a position drawn
at random from the population differ.

**Today the frequency is taken by a detour.** The joint pre-pass fits a curve describing the
population, the caller evaluates that curve into the allele-count classes a panel of `N` diploid
individuals has, and a search finds the two-parameter pair whose own predicted classes best match.

**This spec removes that detour, and it does so in two steps that are worth taking separately.**

---

## 2. Two changes, and the first does not depend on the second

**Step one: integrate the two numbers off the fitted curve in closed form.**

```text
E[2 f (1 − f)]  =  p_segregating · 2ab / ((a + b)(a + b + 1))     `expected_heterozygosity`, exists
E[f]            =  p_fixed_alt  +  p_segregating · a / (a + b)     one line, does not exist
```

A concentration pair is equivalent to those two moments — §6's identity inverts — so this map
reproduces the fitted curve's own first two moments **exactly, at every panel size**. It deletes the
projection and the search, needs no census pass and no inbreeding coefficient, and it is better than
what ships today in **34 of the 36 cells measured** (report §9.1). **What it discards is the curve's
third moment and above**, which is a modelling choice measured to reach no genotype (§7's closed
question).

**Step two, which is the rest of this document: average the two moments over the census positions
instead.** It costs an expectation-step change, a variance term (§3.1) and an inbreeding coefficient
(§3) — and it buys the one thing step one cannot give: **an error that goes to zero as the cohort
grows.** The curve's does not, because a curve converges on the best-fitting member of its family
and the census average has no family to converge to. Measured at 63 individuals on a population with
two frequency peaks: the census average returns 1.000 of the panel's own genotypes at every depth,
the integrated curve 1.05 to 1.07, and today's search 0.75 (report §9.1).

**Step one is not a stepping stone to step two — it is a smaller, independent repair**, and a run
that takes it is strictly better off than today whether or not step two ever happens. **The rest of
this spec is step two.**

---

## 3. The two estimators (step two)

**Both are averages over the census positions**, and both take the panel's alternative-copy count
`k` at a position among its `2N` chromosomes.

```text
mean allele frequency   f  =  mean over positions of   k / 2N

heterozygosity          π  =  mean over positions of   2 k (2N − k)
                                                      ─────────────  ÷  (1 − F/(2N − 1))
                                                       2N (2N − 1)
```

**`F` is the panel's inbreeding coefficient** — the probability that an individual's two copies of a
position are one ancestral copy counted twice. §4 is entirely about it.

**The `2N − 1` rather than `2N` in the denominator** is the finite-panel correction: it is what
makes the answer a property of the population rather than of the panel. Dropping it returns the
heterozygosity `1/(2N)` low, which is 50% at one individual (report §10).

### 3.1 `k` is never counted — it is an expectation, and the second formula needs its variance

**At three reads a position a heterozygote often shows only one of its two alleles and a
sequencing error often looks like a third.** So `k` is an expected count under the read model,
taken from the joint fit's converged per-position posteriors
(`JointFit::genotype_posterior`, three numbers a sample a position):

```text
E[k]     =  Σ over samples i of   P(het at i)  +  2 · P(both copies non-reference at i)

Var(k)  ≈=  Σ over samples i of ( [ P(het at i) + 4 · P(both at i) ]
                                − [ P(het at i) + 2 · P(both at i) ]² )
```

**The square is inside the sum, per sample**, and the parentheses are not decoration: squaring the
whole sum instead is a different and much larger number, and it is the first thing to check against
this section's own §9 test 2.

The first estimator is **linear** in `k`, so substituting `E[k]` is exact.

**The second is not, and this is the part an implementation will get wrong.** `k (2N − k)` is
quadratic, so

```text
E[k (2N − k)]  =  2N · E[k]  −  E[k]²  −  Var(k)
```

**Substituting `E[k]` into the formula and stopping there returns the heterozygosity high by
exactly `Var(k)`.** At one sample and three reads a position that is not a correction, it is most
of the answer: measured, **2.538 ± 0.165 times the truth without the term and 1.219 ± 0.152 with
it** (report §4.1). **At 63 samples the two agree to within 0.003 at three reads and to three
decimals from eight reads up**, so a test written on a cohort would not catch the term's absence —
**the test for this must run at one sample and three reads.**

**`Var(k)` above is the sum of the samples' own posterior variances, and it is not the whole
variance.** The exact quantity is `Σᵢ Varᵢ + ΣΣ_{i≠j} Covᵢⱼ`; the samples at a position are coupled
through the frequency they share, so the covariance is **positive** and the sum above is an
**under**-estimate — which makes the heterozygosity come back slightly **high**, in the same
direction as dropping the term altogether and much smaller.

**Its size has not been measured and this spec does not claim one.** An earlier draft credited it
with the shortfall the report sees at ten samples; that shortfall runs the other way and is ten
times larger, so it belongs to something else (report §4.2). What is known is the size of the whole
variance term at ten samples and three reads — between 1.6 and 2.2 parts in a hundred — which
bounds the residual above by that. **An implementation may not claim the term is exact.** Whether
computing the exact variance is worth it is §8's first open question.

### 3.2 Every kept position counts, unweighted

**Do not weight positions by one minus the fit's posterior that they are mismapped.** The fit
already scores each position under both the ordinary and the mismapped error rate and weights the
two by their posteriors, so a position whose reads look mismapped has already had its genotype
posteriors computed under an error rate that explains them without a heterozygote. Weighting a
second time removes real variation as well as artefact.

**Measured over eighteen cells** — one position in a hundred planted mismapped, at 6 and at 25
disagreeing reads in 100, three panel sizes and two depths (report §6, §6.1). The weighting is
worse in thirteen, better in four and unchanged in one, and every one of those differences is
smaller than its own error bar **except one, which is a loss**: at one sample and three reads with
the stronger plant it moves the estimate from 0.842 of the truth to **0.674**.

**And the plant itself barely registers**, which is the finding behind the recommendation: leaving
the planted positions in, unweighted, changes the estimate by at most 7.4 percentage points against
the same cohort with nothing planted, and by under 1 point at 63 samples.

---

## 4. The inbreeding coefficient, needed by step two and not by step one

**Without `F` the heterozygosity estimator is wrong by `1 − F/(2N − 1)`**, and the size at the two
ends of the committed cohort range is not the same order of thing:

| individuals | shortfall at `F = 0.8` |
|---:|---|
| 1 | **80%** |
| 2 | 27% |
| 3 | 16% |
| 10 | 4% |
| 63 | 0.6% |
| 1000 | 0.04% |

**Measured across four populations at nine panel sizes.** Of the 36 cells, 21 sit within one
standard error of the value the factor predicts, 33 within two and all 36 within three — which is
what a correct formula and 36 draws give (report §3.3).

### 4.1 Where it comes from: the runs estimator, not the fit's homozygote excess

**Two estimators of the coefficient already exist in this tree, and they read different things.**

- **`parameter_estimation::generic::runs`** fits a two-state model over genome windows and returns
  the **share of the analysable genome lying in runs of homozygosity** — stretches where an
  individual's two copies descend from one recent ancestor. It reads the *distribution* of
  heterozygosity along the genome. It returns an `InbreedingF`, works from **one genome**, and needs
  no population expectation.
- **`JointFit::hom_excess`** is how much less heterozygous an individual is than the fitted
  frequencies predict. It reads the *total*, against a population expectation the same fit produced.
  It returns a `HomozygoteExcess`, a deliberately separate type.

**Decision: the correction takes the runs estimator's coefficient.** Three reasons, and the first
is decisive:

1. **Taking the homozygote excess is circular, and this is not a new observation.**
   `parameter_prepass_generic.md` §6.3 states it in as many words: *"Do not take `F` from the ratio
   estimator and then compute the cohort's diversity from it… the ratio estimator needs a diversity
   to produce `F`, so feeding its `F` back in returns whatever was assumed."* This spec's §2
   estimator divides a diversity by `1 − F`; taking that `F` from a quantity fitted against the
   same fit's own diversity closes exactly that loop. **The joint fit is not the pure ratio
   estimator — at two samples and above it also sees how many samples carry the allele at each
   position, which is real information the ratio has not got, and it recovers 0.80 from a truth of
   0.8 there (report §3.5). But the direction of the dependence is the wrong way round**, and the
   runs estimator has no such dependence at all.
2. **It is what the derivation asks for.** §3's factor is the probability that an individual's two
   copies at a position are one ancestral copy — realized autozygosity, which is what a run of
   homozygosity *is*. The homozygote excess additionally absorbs population structure: a cohort that
   is really two subpopulations looks homozygote-excessive for reasons no individual's parents
   caused.
3. **It works at one sample**, which §4.2 shows nothing reading positions one at a time can do.

**And it is already the caller's choice.** `parameter_prepass_generic.md` §6.2 decides the runs
estimator is *"the one a caller reads"*, on three reasons of its own. **This spec adopts that
decision rather than making a second one.**

**What the homozygote excess is still for**: §6.3 of that document keeps it as a **labelled
diagnostic** — once the runs estimator has produced a coefficient and this spec's §2 has produced a
diversity, feed that diversity back through `F = 1 − Hobs/Hexp` and see whether the same coefficient
comes out. Disagreement in a predictable direction is population structure, and it costs one
division a sample.

**⚑ One integration question this spec cannot close.** The runs estimator lives in the per-sample
histogram route and walks genome *windows*; the estimators of §2 live in the joint route and walk
*census positions*. **A run on the joint route alone does not produce a runs coefficient today**, and
whether it should gain one, or whether the two routes' outputs are joined at the run's assembly, is
§8's fifth open question. Until that is settled an implementation may take `JointFit::hom_excess`
and **must say in the output that it did**, because the circularity above is then live.

**The panel's value is the plain mean over samples, unweighted**, whichever estimator supplied it.
With a per-individual `Fᵢ` the estimator's expectation is `π · (1 − F̄/(2N − 1))` with `F̄` the
unweighted mean, so a sample with more census positions covered must not count for more. **Nothing
has tested this on a panel of mixed coefficients** — every drawn panel in the report shares one — so
a weighted rule would have passed every arm there.

**⚑ A scope statement, for the fallback path only.** Where an implementation does take the
homozygote excess (see the integration question above), it absorbs whatever makes an individual less
heterozygous than the fitted frequencies predict — autozygosity, but also population structure and
residual mismapping. Dividing by `1 − F̄` then recovers the **whole** population's diversity rather
than a subpopulation's, which is the same quantity the repeat path's `F` is documented as being. The
runs estimator does not have this ambiguity, which is reason 2 above.

### 4.2 At one sample nothing that reads positions one at a time can estimate it

**The information is absent, not thin.** A single genome shows exactly two things across the
census — how often it is heterozygous, and how often both copies are non-reference — and both are
products of the population's spread of frequencies and `F`. Two populations can be built whose
single-genome census is drawn from **the identical distribution** while their diversities stand in
the ratio `1 − F`: at `F = 0.8`, one at 7.58 differences per kilobase and one at 1.52, both showing
a genome heterozygous at 0.00151515 of positions and homozygous non-reference at 0.00857576
(report §3.4, `examples/ng_prior_moment_one_sample_inbreeding.rs`). No estimator can prefer one
over the other.

**This is not a new defect.** Handed a genome from the selfing population, the fit that ships today
returns a homozygote excess of **exactly 0.000** and a diversity **one fifth** of the truth, at 3,
20 and 100 reads a position alike. `design_principles.md` §0 says *emit it as absent* is a
legitimate answer at one sample and *silently emit a fitted zero* is not; this is the second.

**Decision (owner, 2026-08-27). The coefficient is fitted by default at every cohort size, and a
user may override it — per sample, or with one value for the whole panel, including zero.** The
override is not a single-sample feature: a user who knows how their material was bred knows it
whatever the cohort size, and a fitted coefficient at three samples is worth overriding for the
same reason it is at one.

| what the user gave | what the run uses | what it reports |
|---|---|---|
| nothing | the fitted coefficient, per sample, averaged over the panel | the value, and that it was fitted |
| one value | that value for every sample | the value, and that it came from the user |
| one value per sample | those | the values, and that they came from the user |

**No refusal, and no cohort-size branch.** The earlier draft of this section refused a user-supplied
coefficient above one sample, on the grounds that two ways to set one number is worse than one. That
is wrong twice over: the run reports which it used, so there is no ambiguity to protect against; and
it would deny the override exactly where a user is most likely to be right about their own material.

**What the run owes at one sample depends on which estimator produced the coefficient.**

- **From the runs estimator** — the preferred source (§3.1) — a single genome is a perfectly good
  input, because runs of homozygosity are a within-genome signal. The run reports the coefficient
  and the window count it was fitted over, and warns if that count is near the estimator's own floor
  of 3,000 windows.
- **From the joint fit's homozygote excess**, a single genome returns **0.000 whatever the truth**
  (report §3.4), so the diversity is `π (1 − F)` reported as `π`. **A run in that state must say its
  coefficient came from one genome's totals and is therefore a floor**, and say what the diversity
  would be at a stated coefficient — a true diversity is `1/(1 − F)` times the reported one, so a
  user who knows their material is selfing can multiply.

**Where the second warning stops is a measurement and not a taste, and it has been made**: the
joint fit's coefficient goes from 0.000 at one sample to 0.833 at two, against a truth of 0.8, and
stays within 0.03 of the truth from three samples to sixty-three at both three and twenty reads a
position (report §3.5). **So the warning is for one sample and no other.**

---

## 5. Where the code lives

**The estimators live with the fit, not with the prior.** They are reductions over
`JointFit::genotype_posterior`, which is a census-shaped array the calling side has never seen and
must not learn about (`calling_priors.md` §0 keeps the prior free of a back-reference into its
producer).

- **`parameter_estimation::joint`** gains a function taking the converged posteriors, the sample
  count and the panel's `F`, and returning the two numbers. `JointFit` carries them as fields
  beside `expected_heterozygosity`.
- **`JointFit::expected_heterozygosity` stays**, and stays read off the density. It is the fit's own
  description of the population and other consumers exist for it; what stops is the *caller*
  reading its numbers off the curve. **A run must emit both and they must be comparable** — see
  §6.
- **`calling::run_parameters` is where the seam actually is**, and it is the adapter that goes.
  `FittedFrequencySpectrum::of` calls `FrequencyDensity::allele_count_classes` to build the class
  weights and `::view` hands them to `project_spectrum_seed` (`run_parameters.rs:158`, `:192`).
  That whole type exists to carry a borrowed projection between the two, and with the projection
  gone it has nothing to carry.
- **`calling::genotype_prior::seed_generic` loses the projection and the search — and it loses them
  at step one, before any of the rest of this spec is built.**
  `fit_spectrum_shape`, `fill_expected_spectrum` and `FittedSpectrum` have no consumer outside it
  and `run_parameters`; `FrequencyDensity::allele_count_classes` loses its only production caller.
  **They are deleted, not deprecated**, and what replaces `project_spectrum_seed` takes the two
  measured numbers and returns the pair by §6's identity.
- **What must not be deleted with them is the variable-census-site count.**
  `FittedFrequencySpectrum::of` computes it as one minus the two end classes, which is *how many
  positions came out variable across this panel* rather than how many segregate in the population,
  and its own comment records the two differing 6.6-fold at one individual on a tomato-like
  density. §6.2 needs the same quantity and it now comes from the posteriors directly — the count
  of positions that segregate, in the soft form §6.2 defines.

**The `genotype_posteriors` flag turns on.** It is off by default today because it weighs twelve
bytes a position a sample — 1.5 GB for fifty samples over two million positions. **Two ways to pay
it, and the choice is the implementer's**: keep the flag off and accumulate the two moments inside
the E-step's own per-position loop, where the posteriors are already in a scratch buffer and
nothing is stored; or turn the flag on and reduce afterwards. **The first is what the design wants**
— the moments are two running sums and the array exists only to be summed — and the second is a
legitimate first cut for getting the numbers to agree with the harnesses. Whichever ships, the
memory it costs is stated in the module's own documentation.

---

## 6. What the seed does with the two numbers

**Unchanged from [`ordinary_site_seed.md`](ordinary_site_seed.md) §3**, and this spec adds nothing
to it. A pair of expected frequency `f` and total `A` implies a heterozygosity of
`2 f (1 − f) · A / (A + 1)`, so the total that reproduces a measured `π` is

```text
A  =  π / (2 f (1 − f) − π)
```

and the seed is `(A (1 − f), A f)`. The three ways this can have no answer, and the requirement
that none of them be silent, are that document's §3.1 and stand as written.

### 6.1 The blend toward the neutral shape goes away

**`ordinary_site_seed.md` §4 interpolates the expected frequency between a neutral shape and the
panel's own, at a weight `w = N / (N + N₀)` that rises with the panel.** That blend exists because
the *search's* shape was untrustworthy at a small panel — three allele-count classes at one
individual carry almost no shape.

**The direct estimate has no such weakness, and the blend now costs rather than protects.** Handed
the population exactly, the search at one individual is exact to three decimals — and the blended
seed is **0.62× to 0.92× the truth** across four populations, purely from the pull toward the
neutral shape (report §9). What the blend was insuring against is not a risk the direct estimate
runs: at one individual the direct mean frequency is unbiased on all four populations (report §3.1),
so there is no small-panel shape to shrink away from.

**The blend also damped *noise* at a small panel, and there is none left to damp.** At one
individual and the shipped two-million-position census, the direct mean frequency's spread over
runs is **0.77% to 1.62% of the population's** across the four populations and both inbreeding
coefficients, and 19 runs in 20 land within 1.5% to 3.2% (report §7). It needs no inbreeding
coefficient either. **So a single-genome run's shape is not the noisy thing the blend was built
for.**

**Decision: no blend. The panel size does not appear in the seed at all.** That is the whole point
of the change — the prior is a statement about the population, and `N₀` was a constant that had to
be fitted for a defect this removes.

**What replaces the protection it gave** is §6.2's floor, which keys on the census rather than on
the panel.

### 6.2 What happens where the census is thin

**A thin census is a real risk and it is a different one from a small panel.** Both moments are
averages over positions, and their precision rests on how many positions **segregate** — a
position where the population carries one allele contributes zero to the heterozygosity and tells
you nothing about the frequency. A two-million-position census on a panel segregating at 1 in 200
carries about ten thousand segregating positions
(`parameter_prepass_census_sites.md` §6.1); one run over a small `--regions` BED may carry a few
hundred.

**Decision: the run reports the count and the spread; it does not silently fall back.**

- **A count of segregating positions travels with the two moments**, as the variable-site count
  already does for the spectrum. **It cannot be "positions whose expected alternative-copy count is
  above zero", and an earlier draft of this spec said that**: posteriors from reads are continuous,
  so that count is essentially every position and the run would report 100% segregating. It must be
  a **soft count** — `Σ over positions of P(the position segregates in this panel)`, which the fit
  can form from the same per-position quantities — or a count above a stated threshold, and the
  implementation must say which. The soft count is preferred: it needs no constant, and it degrades
  smoothly where the reads are thin instead of stepping.
- **The spread of each moment travels with it too.** It is computed across positions, and **it must
  be labelled as a floor rather than an interval**: census positions are scattered but linked, so a
  spread computed as though they were independent is too narrow by a factor
  `parameter_prepass_census_sites.md` §5 puts between 3 and 16.
- **Below a floor on the segregating count the run falls to the neutral shape at the measured
  diversity** — `ordinary_site_seed.md`'s middle rung, unchanged — and says so.

**Where the floor goes is not set here, and the reason is that nothing has measured it.** §7's
second open question names the experiment. **An implementation may not pick a number for it**; until
it is measured the run reports the count and takes no floor, which is the honest state and is
distinguishable in the output from a floor that never fires.

---

## 7. What the run reports

**A run that used different information must not look the same as one that did not** — the
requirement `calling_priors.md` §4 makes about production's silent fallback. Beyond what
`SeedRegime` already carries:

| what | why it must be there |
|---|---|
| the two measured moments | they are the run's estimate of the population and a user may want them |
| the fit's own `expected_heterozygosity`, read off the density | **two independent estimates of one quantity, and a large gap between them is the run telling you something** — the direct one is a census average, the curve's is a model fit, and they should agree |
| how many census positions segregated | a run returning a few hundred has thin moments and good error rates, and only this says which |
| the spread of each moment, labelled a floor | §6.2 |
| where the inbreeding coefficient came from — the runs estimator, the joint fit's homozygote excess, or the user | §3.1, and the middle one carries a circularity the output must not hide |
| for the runs estimator, how many windows it was fitted over | its own floor is 3,000, below which what it returns is its noise (`parameter_prepass_generic.md` §6.1) |

**The gap between the two heterozygosity estimates is worth naming rather than leaving to be
noticed.** They are computed from the same converged fit by two routes: one averages the
per-position posteriors, the other integrates the fitted curve.

**Measured at 63 samples across three populations and four depths, the curve's number sits between
1.1% below and 10.7% above the census's own**, and which population is the wide one is *not*
predicted by whether the curve can hold its shape — the widest cell, 10.7%, is the rare-allele
population, which it can hold exactly:

| | 3 reads | 8 | 20 | 100 |
|---|---:|---:|---:|---:|
| nearly all alleles rare (inside the family) | 1.040 | 1.063 | **1.107** | 1.069 |
| two peaks (outside it) | 1.039 | 1.042 | 1.041 | 1.054 |
| one peak (inside it) | **0.989** | 0.992 | 1.021 | 0.995 |

**So a threshold on this gap cannot be set at a tenth**: a converged, healthy fit already reaches
10.7% here. What the number is good for is being *printed* — two routes to one quantity, and a run
that shows them far apart is a run worth looking at. **Where the threshold goes is §7's fourth open
question**, and it needs a fit that genuinely failed to converge to calibrate against, which nothing
in this work produced.

---

## 8. Open questions

1. **Whether the exact `Var(k)` is worth computing, and how large the residual actually is.**
   §3.1's variance is a sum of per-sample variances and ignores the positive coupling between
   samples at a position. **Nothing has measured what that costs** — the report bounds it above by
   the whole variance term's size, 1.6 to 2.2 parts in a hundred at ten samples and three reads,
   and no lower. The exact quantity is available inside the expectation step, where the position's
   frequency is already being integrated over: `Σ_nodes weight(node) · Var(k | f = node)` plus the
   variance of `E[k | f]` across nodes. So this is a question about cost, not about reachability —
   and computing it once would settle the size at the same time.
2. **Where the thin-census floor goes.** Subsample a real cohort's census down and watch where the
   two moments stop being stable, which is the shape of experiment
   [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §10's third question already asks
   for. **It cannot be answered on drawn data alone**: what makes a real census thin is linkage,
   and drawn positions are independent.
3. **What the single-sample flag is called, and which of the two coefficients it sets.** §4.1 keeps
   the homozygote excess and the genotype prior's autozygosity coefficient as separate quantities,
   which the code already does. **A user at one sample plausibly wants to set both**, and whether
   one flag sets two things or there are two flags is an interface decision this spec does not take.
4. **Where the threshold on §6's two-route disagreement goes.** A converged, healthy fit already
   shows the curve's heterozygosity 10.7% above the census's own on one of the three populations
   measured, so a threshold at a tenth would fire on good runs. Calibrating it needs a fit that
   genuinely failed to converge, and nothing in this work produced one.
5. **How a joint-route run gets a runs coefficient.** §4.1 prefers it and the joint route does not
   produce one: the runs estimator walks genome windows in the per-sample histogram route, and this
   route walks census positions. Three shapes are possible — run the windows estimator alongside,
   fit runs from the census positions themselves, or join the two routes' outputs at the run's
   assembly — and **the second is worth checking before the others**, because the census is dense
   enough in principle: two million positions over 800 Mb is one position every 400 bases
   (`parameter_prepass_census_sites.md` §1), so a megabase run of homozygosity holds about 2,500
   census positions, and at a diversity of 7 per kilobase an equally long stretch outside one would
   hold roughly 17 heterozygotes against nearly none inside. **That is arithmetic, not a
   measurement**, and it says the signal is there rather than that a model reads it.
6. **Whether a panel of mixed inbreeding coefficients behaves as §4.1 says.** The unweighted mean is
   right by derivation, and every drawn panel behind this spec shares one coefficient across its
   individuals — so a weighted rule would have passed every measurement made.

**One thing deliberately left shut.** *A better curve* is not an alternative to any of this: the
prior can hold two moments, so the curve's shape was never reaching a genotype. That was measured
in July on the repeat path and came back a null
([`../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md`](../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md)),
and this spec is the explanation of that null from the other side.

---

## 9. How we will know it works

**Three tests that fail if the code is wrong, and each names the case that catches it.** The
recurring failure in this area is a test that cannot fail, so each entry says what it would miss.

1. **The estimators against known genotypes, at one individual and at a thousand.** Feed
   posteriors that are point masses on known genotypes; both moments must return the census's own
   exact values. **At one individual and `F = 0.8`**, because that is where the inbreeding factor is
   `1 − F` and a missing one is an 80% error; at a thousand, because that is where a `2N` in place
   of `2N − 1` is a 0.05% error and only the one-individual case would catch it.
2. **The variance term, at one sample and three reads a position.** Posteriors midway between
   genotypes — the shape low depth produces — must give a heterozygosity below what substituting
   the mean gives, by the variance. **A cohort test cannot catch a missing variance term**: at 63
   samples the two agree to three decimals (§3.1).
3. **The seed's implied heterozygosity is the measured one, at every panel size and shape.** The
   identity in §5 is exact, so this is a real check on the arithmetic rather than a restatement of
   it — and it is the test `seed_generic` already carries for the pinned total, extended to the
   new source of `f`.

**And the harnesses stay.** The three programs behind the report are the regression check for
anything that changes these numbers; a change that moves them and leaves the tests green is a
change whose effect nobody has looked at.
