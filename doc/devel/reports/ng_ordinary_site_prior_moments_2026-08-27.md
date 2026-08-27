# The genotype prior's two numbers, taken straight from the census — six questions answered

**Status:** measurement report, 2026-08-27. Answers the six questions in
[`../ng/research/ordinary_site_prior_moments.md`](../ng/research/ordinary_site_prior_moments.md).
Everything here is drawn cohorts; this checkout cannot rebuild the tomato census, which the plan
records as a standing limit.

**A short answer written for the plan's author**, carrying the conclusions without the
measurements, is
[`../ng/research/ordinary_site_prior_moments_answer.md`](../ng/research/ordinary_site_prior_moments_answer.md).

**The three programs behind it:**
[`examples/ng_prior_moment_estimators.rs`](../../../examples/ng_prior_moment_estimators.rs) — the
estimators on genotypes that are known, and what the caller's current path returns on the same
populations;
[`examples/ng_prior_moments_from_reads.rs`](../../../examples/ng_prior_moments_from_reads.rs) — the
same estimators on reads, through the joint fit;
[`examples/ng_prior_moment_one_sample_inbreeding.rs`](../../../examples/ng_prior_moment_one_sample_inbreeding.rs)
— two populations five-fold apart in diversity that one genome cannot tell apart.

---

## 1. The recommendation, first

**Two changes, and the small one is not a step toward the large one — it is worth making whether or
not the large one happens.**

1. **Delete the projection and the search, and integrate the two numbers off the fitted curve in
   closed form.** Two lines. It removes the whole measured defect, needs no census pass and no
   inbreeding coefficient, and is better than what ships today in 34 of the 36 cells measured. §9.1.
2. **Then decide whether to average over the census positions as well**, which is the plan's own
   proposal. It is the only route whose error goes to zero as the cohort grows — 1.000 against the
   panel's genotypes at 63 individuals on a population the curve's family cannot hold, where the
   curve settles at 1.05 and stays. **It costs three things the curve route does not need**: the
   fitting loop must accumulate two running sums; the heterozygosity needs a variance term without
   which it returns 2.5× the truth at one sample and three reads; and it needs the panel's
   inbreeding coefficient. **The decision turns on two measured facts** — above eight reads a
   position the choice moves no genotype at all, and at three reads it can move up to 3 calls in 100
   at segregating loci (§9.2). **So it is a decision about low-coverage calling and about how
   accurate the run's printed diversity has to be, not about calling in general.**

**The argument that started this work — that four numbers cannot be compressed to two without
biasing the prior — identified a real defect and named the wrong cause.** The prior's two numbers
*are* two exact integrals of the four-number curve; the bias belongs to how the caller currently
picks them, which is a divergence fit over the panel's allele-count classes rather than a moment
match. §9.1.

**The proposal**, written up in
[`ordinary_site_prior_moments.md`](../ng/research/ordinary_site_prior_moments.md) and called *the
plan* throughout this report, is to compute the prior's two numbers — the population's **mean
alternative-allele frequency** and its **heterozygosity**, how often two copies of a position drawn
at random from the population differ — from the census positions directly, rather than by fitting a
curve and searching for a two-parameter pair that reproduces it.

**The mean frequency comes back unbiased at every panel size and on every population**, as the plan
expected. **The heterozygosity does too — but only when the panel mates at random**, and that is
the finding this report leads with, because the caller's own benchmark cohort does not: tomato
accessions are self-pollinated, and the fit puts their inbreeding coefficient at 0.8 to 0.9.

Three things follow, in the order they matter:

- **The heterozygosity estimator needs a third input the plan does not name: the panel's
  inbreeding coefficient `F`.** Without it, a single self-pollinated tomato accession reports
  **one fifth** of its population's diversity, and a panel of two reports three quarters of it.
  The correction is a single factor, `1 − F/(2N − 1)`, derived rather than fitted, and from two
  samples up the joint fit already estimates the `F` it needs. §3.
- **The inbreeding coefficient is needed by the census route and not by the closed-form one**, and
  that is the largest practical difference between them. The curve's moments are properties of a
  population with no individuals in them, so nothing has to be undone; the census average works over
  realized genotypes and carries a factor of `1 − F/(2N − 1)`. **Where the coefficient is needed, it
  should not come from the joint fit.** The joint fit's homozygote excess is measured against a diversity the same fit
  produced, so dividing that diversity by `1 − F` closes a loop `parameter_prepass_generic.md` §6.3
  already warns about. **The runs-of-homozygosity estimator in `generic::runs` has no such
  dependence, reads the coefficient off where heterozygotes sit along the genome, and works from a
  single genome** — and §6.2 of that document already decided it is the one a caller reads. §3.1 of
  the accompanying spec. **The joint fit's own coefficient is 0.000 at one sample whatever the
  truth, and within 0.03 of it from two samples up** — so the two matter in different places. §3.4,
  §3.5.
- **The defect being removed is worth up to a quarter of the number.** Handed the population
  exactly, today's projection-and-search returns a mean frequency between **1.22× and 0.79× the
  truth at 200 individuals depending on the population's shape**, and the error grows steadily with
  the panel in all four. On fitted curves it is worst at **0.749×** — 63 individuals, a population
  with two frequency peaks — where integrating the same curve gives 1.05× and the census average
  gives 1.000×. §9, §9.1.

**What does not follow: an accuracy claim about called genotypes.** Seeding the caller from the
direct moments instead of from the search moves almost no genotype calls. §8 gives the number and
the control that says which of those cells could have seen movement at all.

---

## 2. What was measured, and what each number is scored against

**This is stated before the results because getting it wrong is how a measurement invents a
finding**, and it did so on the immediately preceding branch.

Two different yardsticks are used, never mixed:

| question | scored against | why |
|---|---|---|
| **Is the estimator biased?** | the moments of **the positions that cohort was itself drawn at** | those positions are a finite sample of the population, and their own moments sit a few percent away from it. Scoring against the population would put that scatter into every error column and call it estimator bias. |
| **How far from the population is the number a run prints?** | the population's own moments, in closed form | a run does not get to average over censuses, so this is what a user sees. |

**The populations.** Four, at diversities from 6 to 15 differences per 10,000 bases — which spans
this project's two benchmark cohorts, tomato at about 6 and a human panel at about 10, and reaches
half again above the higher of them. Two of them are shapes the joint
fit's own description of a population can hold exactly (one Beta over the positions that segregate,
plus a spike of positions carrying only the reference base and a spike carrying only a
non-reference one). **Two of them are not**, and that is deliberate: on a population drawn from the
fit's own family, an estimator that quietly assumes that family looks perfect.

How far outside is measured rather than asserted — the largest gap between a population's spread of
segregating frequencies and the closest single Beta's, over 400,000 draws from each:

| population | heterozygosity, per kb | mean frequency | gap from the closest single Beta |
|---|---:|---:|---:|
| nearly all alternative alleles rare, tomato-like — `Beta(0.2, 1)` | 0.61 | 0.00167 | 0.016 |
| where it varies, the reference base is the rare one — `Beta(3, 0.6)` | 0.87 | 0.00933 | 0.001 |
| **two peaks, off centre** — 0.7 of `Beta(30, 70)` and 0.3 of `Beta(90, 10)` | 1.38 | 0.00242 | **0.269** |
| **a lump at one intermediate frequency** — 0.55 of `Beta(0.25, 1.5)` and 0.45 of `Beta(42, 58)` | 1.51 | 0.00234 | **0.192** |

A gap of 0.27 means that at some frequency the two curves disagree about **27% of all segregating
positions**.

**The first two rows are 0.000 by construction, and the column reads 0.016 and 0.001 because this
statistic has a noise floor.** Moment-matching a Beta to a single Beta returns that same Beta
exactly, so those two populations have no gap at all; what is printed is the disagreement between
two finite samples of one distribution. Comparing two sets of 400,000 draws has a floor near 0.002
even when the distributions are identical, and 0.016 is eight times that — worth checking whether
drawing the two samples alternately from one generator correlates them, which this report has not
done. **The two rows that matter are unaffected**: 0.269 and 0.192 are more than a hundred times
that floor.

**Inbreeding.** `F` is the probability that an individual's two copies of a position are one
ancestral copy counted twice. **It is not swept everywhere, and where it is not is worth knowing:**
§3's estimator sweep runs at `F = 0` (random mating) and at `F = 0.8`, tomato's fitted range; §7's
one-genome arm runs both; §4 to §6 and §8 run at a single homozygote excess of 0.15, because every
arm there is a full joint fit and the ratio they report is insensitive to it; §9's comparison
against the current path runs at `F = 0` only, for a reason that section gives.

---

## 3. Question 1 — do the estimators recover the truth?

**Sweep:** nine panel sizes from 1 to 1,000 diploid individuals, four populations, two inbreeding
coefficients, 20,000 census positions and 24 independently drawn cohorts a cell. The positions and
every sample's genotype are drawn once a replicate and each panel-size arm reads the first `N`
samples of that same draw, so nothing moves across the arms for a reason unrelated to panel size.
No reads: these are the genotypes the cohorts were drawn with.

### 3.1 The mean allele frequency is unbiased everywhere

`mean over positions of k / 2N`, for `k` alternative copies among the panel's `2N` chromosomes.

**Across all eight population-by-inbreeding cells and all nine panel sizes, the largest departure
from the drawn positions' own mean frequency is 3.4%, at one individual, against a precision of
±1.5% for that cell.** From ten individuals up, the worst cell anywhere in the sweep is 0.45%; from
200 up it is 0.19%; at 1,000 individuals every cell is within 0.07%, against a precision of
±0.06%.

Inbreeding does not touch it, and cannot: `E[k / 2N]` is the frequency whatever the two copies
inside an individual are doing.

### 3.2 The heterozygosity is unbiased under random mating

`mean over positions of 2 k (2N − k) / (2N (2N − 1))` — Nei's average heterozygosity, where the
`2N − 1` rather than `2N` is the finite-panel correction that makes the answer a property of the
population rather than of the panel.

**At `F = 0`, across all four populations and all nine panel sizes, the largest departure is 4.5%,
at one individual, against a precision of ±5.9%.** From ten individuals up every cell is within
1.8%; at 1,000 individuals every cell is within 0.11%.

### 3.3 At `F = 0.8` it is not unbiased, and the size is not a detail

The same estimator, same populations, `F = 0.8`:

| individuals | tomato-like | reference base rare | two peaks | a lump |
|---:|---:|---:|---:|---:|
| 1 | **−81.7%** | **−81.1%** | **−77.7%** | **−80.7%** |
| 2 | −23.4% | −24.4% | −28.5% | −26.2% |
| 3 | −15.5% | −11.8% | −15.4% | −16.7% |
| 5 | −8.7% | −6.9% | −10.3% | −8.8% |
| 10 | −3.7% | −2.8% | −4.8% | −4.1% |
| 25 | −2.9% | −0.7% | −2.4% | −1.9% |
| 63 | −1.1% | −0.2% | −1.4% | −0.8% |
| 200 | +0.3% | −0.2% | −0.4% | −0.2% |
| 1000 | +0.0% | +0.2% | −0.1% | −0.0% |

Every cell above is measured to a precision of between ±0.05% and ±2.9%. **No entry from one individual to ten is
within reach of zero**; from 25 up two are, both on the reference-base-rare population — −0.67% ±
0.68% at 25 individuals and −0.23% ± 0.46% at 63 — which is what the factor predicts there, since
it puts the shortfall at 1.6% and 0.64%.

**Against what the factor `1 − F/(2N − 1)` predicts** — −80.0% at one individual, −26.7% at two,
−16.0% at three, −4.2% at ten, −0.64% at 63 and −0.04% at 1,000 — **21 of the 36 cells sit within
one standard error of the prediction, 33 within two, and all 36 within three.** For 36 draws from a
correct formula the expected counts are about 25, 34 and 36. The worst cell is 2.5 standard errors
out; there is no cell where the formula and the measurement part company.

**The mechanism, and it is arithmetic rather than a defect in anything.** The estimator asks how
often two chromosomes drawn at random from the panel differ. One pair in `2N − 1` is the two copies
inside a single individual, and those are the same ancestral copy with probability `F`, in which
case they never differ. So the estimator's expectation is the population's heterozygosity times

```text
1 − F / (2N − 1)
```

At one individual there is only that one pair, and the factor is `1 − F` exactly: **a single
self-pollinated genome counting its own heterozygous positions measures `π (1 − F)`, not `π`.** At
63 individuals the same factor is 0.994 and nobody would notice.

**Checked independently of the sweep**, at one known frequency of 0.30 with the panel redrawn
400,000 times and no census sampling in it at all — the estimator's mean divided by
`2 f (1 − f) · (1 − F/(2N−1))`:

| individuals | `F = 0` | `F = 0.8` |
|---:|---:|---:|
| 1 | 0.998 | 1.007 |
| 5 | 1.001 | 1.001 |
| 50 | 1.000 | 1.000 |

**Dividing by the factor restores the estimator.** At two individuals and above, every cell of
every population returns to within 5.0% of the drawn positions' own heterozygosity — the worst is
+5.00% ± 2.26% at three individuals — and within 0.26% at 1,000. At one individual the correction is a division by `1 − F = 0.2`, which multiplies
the estimate and its scatter alike by five: the bias comes back to between −8.5% and +11.3%, but
the run pins those to only ±8 to ±15%, so what the correction buys at one individual is a
recentred estimate with five times the spread.

### 3.4 At one sample nothing can supply the correction, and the caller already has this defect

**The plan writes the heterozygosity estimator with two inputs — the allele counts and the panel
size. It needs three.** From two samples up the third is available: the joint fit already estimates
a per-sample homozygote excess (`JointFit::hom_excess`), which is this coefficient.

**At one sample it is not available from the sample, and the reason is not that the estimate is
noisy — it is that the information is absent.** A single genome shows exactly two things across
the census: how often it is heterozygous, and how often both its copies are non-reference. Writing
`s` for the share of positions that segregate, `q` for the share fixed non-reference, and `m1`,
`m2` for the first two moments of the population's spread of frequencies over what segregates:

```text
heterozygous          (1 − F) · s · 2 (m1 − m2)
both copies non-ref   q + s · m2 + F · s · (m1 − m2)
```

Hold the spread fixed and both observables equal between an outbred population and a selfing one,
and the two equations solve exactly: `s_outbred = (1 − F) · s_selfing` and
`q_outbred = q_selfing + F · s_selfing · m1`. **The resulting pair of populations produce a single
genome's census from the identical distribution, position for position and read for read, while
their diversities stand in the ratio `1 − F`.**

Built at `F = 0.8` (`examples/ng_prior_moment_one_sample_inbreeding.rs`): a selfing population at
**7.58 differences per kilobase** and an outbred one at **1.52**, five-fold apart, both showing a
single genome heterozygous at 0.00151515 of positions and homozygous non-reference at 0.00857576 —
agreeing to one part in `10^18`, which is floating-point noise. Since the data distributions are
the same, no estimator can prefer one over the other.

**⚠ That holds for a census read one position at a time, and no further.** Both populations here
draw every position independently, so a selfed genome in this construction has its homozygosity
scattered evenly across the genome. **A real selfed genome does not**: it is a mosaic of long
stretches where both copies descend from one ancestor and carry almost no heterozygotes, and
stretches where they do not. **That mosaic is information one genome carries and this construction
deliberately removes**, and §3.5 is about reading it.

**What the shipped fit does there is the part that matters, and it is not neutral.** One genome,
200,000 census positions, `fit_jointly` asked for nothing special:

| depth | fitted diversity, selfing genome | over its truth | fitted diversity, outbred genome | over its truth | fitted homozygote excess, both |
|---:|---:|---:|---:|---:|---:|
| 3 | 0.001453 | **0.192×** | 0.001447 | 0.955× | 0.000 |
| 20 | 0.001346 | **0.178×** | 0.001495 | 0.987× | 0.000 |
| 100 | 0.001534 | **0.202×** | 0.001529 | 1.009× | 0.000 |

Two readings, and the second is the important one:

1. **The two rows at each depth agree to within their own sampling error**, as the argument says
   they must — the fit returns about 0.0015 whichever population the genome came from, at every
   depth. A genome is heterozygous at about 303 of these 200,000 positions, so each row carries
   roughly 6% of sampling scatter and the two together about 8%: the 0.4% and 0.3% gaps at three and
   a hundred reads are well inside that, and the 11% gap at twenty is 1.4 times it. **Depth does
   not help and cannot** — more reads pin the genome's genotypes better, and the genotypes are what
   carry no information about `F`.
2. **The fit resolves the ambiguity at the bottom of its own range**: it returns a homozygote
   excess of **0.000** for a genome drawn from a population with `F = 0.8`, at every depth, and
   therefore reports **a fifth** of that population's diversity, silently.

   **Whether that zero is a fitted interior answer or the search pinned at its lower endpoint is not
   settled here, and the distinction matters.** The maximisation is a golden-section search over
   `[0, 1]` (`fit.rs:2648`) and the type refuses anything below zero, whose own note says *"an
   unconstrained fit will go negative under a heterozygote excess"* — so an endpoint hit is a real
   possibility and would be a **diagnosable state a run could report**, which is a better outcome
   than the one this report assumes. Settling it needs the likelihood printed at a few coefficients,
   which nothing here did.

   Either way, `design_principles.md` §0's rule for the one-sample end stands: *emit it as absent*
   is a legitimate answer and *silently emit a fitted zero* is not, and what a user gets today is
   the second.

**So the inbreeding shortfall at one sample is not a cost of the proposal. It is a defect the
caller has today**, and the direct estimator's virtue here is that it makes the mechanism explicit:
the missing factor has a name, a value, and one place to be supplied from. **A single-genome run
must take `F` from the user, or report its diversity as a lower bound and say so.**

### 3.5 Two samples is enough for the cohort route, and one genome is not out of reach either

**Two questions the recommendation rests on, and the first has a sharp answer.** §3.4 shows one
genome cannot separate the coefficient from the diversity; it says nothing about how many samples it
takes before the fit actually finds it. Measured on one selfing population at `F = 0.8`, 50,000
census positions, the cohort drawn afresh at each size:

| individuals | fitted coefficient at 3 reads | at 20 reads | fitted diversity at 3 reads, over the truth | at 20 reads |
|---:|---:|---:|---:|---:|
| 1 | **0.000** | **0.000** | 0.246× | 0.198× |
| 2 | 0.833 | 0.825 | 0.911× | 1.052× |
| 3 | 0.805 | 0.829 | 1.038× | 1.011× |
| 5 | 0.809 | 0.797 | 1.014× | 0.989× |
| 10 | 0.814 | 0.813 | 1.090× | 1.126× |
| 25 | 0.804 | 0.813 | 1.089× | 1.091× |
| 63 | 0.803 | 0.802 | 1.157× | 1.186× |

**The cliff is between one sample and two, and there is no ramp.** At two samples the fit returns
0.833 and 0.825 against a truth of 0.8; from three up every cell is within 0.03 of it, at three
reads a position as much as at twenty. **Depth is not what makes the difference — a second genome
is.** The reason is visible in the drawer: the samples at a position share that position's
frequency, so how many of them carry the allele says what the frequency is, and the excess of
homozygotes over what that frequency predicts is the coefficient. One genome never learns the
frequency at any position, so it cannot form that excess.

**The diversity column is what a run reports today, uncorrected.** At one sample it is a fifth of
the truth; from two up it is within 19% and mostly within 10%, because the fitted curve is being
fitted against genotypes whose homozygote excess the fit has separately accounted for.

**And one genome is not the dead end §3.4 makes it look**, because the caller already carries a
second estimator that reads the coefficient off the genome's *structure* rather than off its total:
`parameter_estimation::generic::runs` fits a two-state model over genome windows and returns the
share lying in runs of homozygosity. **It needs no cohort, and it needs no population expectation**
— which is exactly the property `parameter_prepass_generic.md` §6.3 says makes a cohort's diversity
estimable at all, and §6.2 already decides it is *"the one a caller reads"*.

**What that costs and what it needs are documented and not small:** at least 3,000 windows, which a
tomato genome (8,004) and a human one (31,000) clear easily and a region-restricted run may not; and
a starting spread that spans the two states' separation, without which it returns a **silent zero**
— the same failure mode as the one this report found in the joint fit, from a different cause.

**Nothing in this report measured it**, and it changes which input the correction should take. §9 of
the accompanying spec now takes it as the preferred source and says why.

---

## 4. Question 2 — what does it cost to read from reads instead of genotypes?

**Nobody can count `k`.** At three reads a position a heterozygote often shows only one of its two
alleles and a sequencing error often looks like a third, so the allele count has to be an *expected*
count under the read model, taken from the joint fit's own per-position posteriors
(`JointFit::genotype_posterior`).

**Sweep:** three panel sizes (1, 10, 63 individuals), four depths (3, 8, 20, 100 reads a position),
10,000 census positions, three independently drawn cohorts a cell, homozygote excess 0.15 in every
drawn individual. Both arms — the genotypes the cohort was drawn with, and the posteriors the fit
returned from its reads — see **the same positions in the same cohort**, so the ratio between them
carries no census sampling at all. Each cell is quoted with one standard deviation over the
replicates divided by their number's square root; **a departure from 1.000 smaller than that is one
the run cannot see.**

**⚠ These are not §3's four populations, and the difference matters.** Every arm from §4 to §8 is a
full joint fit, so a two-million-position census is out of reach; what sets how well the fit
resolves a population is the number of positions that **segregate**, and at tomato's rate of 1 in
200 a ten-thousand-position census would carry only fifty. **The populations here segregate at 2 in
100 instead**, which puts about 200 in each drawn census — a fiftieth of a real run's ten thousand
rather than a two-hundredth. Their shapes are §2's, their diversities are not:

| population | shape | heterozygosity, per kb | mean frequency | used in |
|---|---|---:|---:|---|
| nearly all alternative alleles rare | one `Beta(0.2, 1)` — inside the fit's family | 3.03 | 0.00533 | §4, §6, §8 |
| two peaks, off centre | 0.7 of `Beta(30, 70)` and 0.3 of `Beta(90, 10)` — outside it | 6.89 | 0.01060 | §4, §5, §6, §8 |
| one peak at the same place | one `Beta(30, 70)` — inside it, and §5's control | 7.24 | 0.00622 | §5 |

**The diversity level is not what §4 to §8 measure.** What they measure is the *ratio* between
reading the moments off reads and reading them off the same positions' genotypes, and both arms see
the same positions. What the level does affect is how sharply the fit resolves the population, and
that is a limit on this whole half of the report — §12 says so.

### 4.1 The second formula needs a term the plan does not write, and at one sample it is most of the answer

The frequency estimator is linear in `k`, so substituting the posterior mean is exact. **The
heterozygosity estimator is not**: `k (2N − k)` is quadratic, and

```text
E[k (2N − k)]  =  2N · E[k]  −  E[k]²  −  Var(k)
```

Dropping `Var(k)` — substituting the posterior mean and evaluating the formula — returns the
heterozygosity **high** by exactly the variance the reads left behind. At one individual and three
reads a position that is not a correction, it is the bulk of the number:

| population | plain substitution | with the variance term |
|---|---:|---:|
| rare alleles, one individual, 3 reads | **2.538 ± 0.165×** | 1.219 ± 0.152× |
| two peaks, one individual, 3 reads | **1.858 ± 0.049×** | 0.916 ± 0.057× |
| rare alleles, one individual, 8 reads | 1.185 ± 0.051× | 1.024 ± 0.042× |
| rare alleles, one individual, 20 reads | 1.003 ± 0.002× | 0.977 ± 0.001× |
| rare alleles, 63 individuals, 3 reads | 0.959 ± 0.019× | 0.957 ± 0.019× |

**So an implementation that writes the plan's formula with an expected count substituted into it
returns two and a half times the truth at the corner this caller commits to.** The variance used
here is the sum of the samples' own posterior variances, which is exact only if the samples are
independent given the reads; they are not, and the residual is inside the numbers below.

### 4.2 What it costs, once that term is there

Heterozygosity from the posteriors divided by heterozygosity from the same positions' genotypes:

| individuals | 3 reads | 8 reads | 20 reads | 100 reads |
|---:|---:|---:|---:|---:|
| **rare alleles** | | | | |
| 1 | 1.219 ± 0.152 | 1.024 ± 0.042 | 0.977 ± 0.001 | 0.973 ± 0.001 |
| 10 | 0.903 ± 0.025 | 0.900 ± 0.020 | 0.950 ± 0.005 | 0.920 ± 0.048 |
| 63 | 0.957 ± 0.019 | 0.957 ± 0.011 | 0.984 ± 0.004 | 0.965 ± 0.006 |
| **two peaks** | | | | |
| 1 | 0.916 ± 0.057 | 0.929 ± 0.019 | 0.988 ± 0.003 | 0.983 ± 0.001 |
| 10 | 0.753 ± 0.036 | 0.862 ± 0.015 | 0.923 ± 0.016 | 0.917 ± 0.008 |
| 63 | 1.001 ± 0.002 | 1.000 ± 0.001 | 1.000 ± 0.000 | 0.999 ± 0.001 |

The mean frequency costs less everywhere: the worst cell in the whole sweep is 0.878 ± 0.017 (two
peaks, ten individuals, three reads) and every 63-individual cell is 0.985 or better.

**⚠ One deviation from what a real run would do, declared because it affects the `today` columns and
nothing else.** The shipped seam hands the run's own inbreeding coefficient to the search
(`run_parameters::project_seed`), and these arms hand it zero, on cohorts drawn with a homozygote
excess of 0.15. The reason is that the step before the search — evaluating the fitted curve into the
panel's allele-count classes — takes no coefficient at all
(`FrequencyDensity::allele_count_classes`), so passing a non-zero one to the search alone compares
it against a spectrum no such panel produces. **What it costs is a few percent on the `today`
columns**: `project_spectrum_seed`'s own documentation records the reference concentration moving
8.6% across `F = 0` to `0.6`, so at 0.15 it is smaller than that. §8's comparison between the two
seeds is unaffected in its conclusion — the control there shows that a **40%** difference in the
alternative concentration moves no calls, and this is a few percent.

Three readings:

- **At 63 individuals, reading from reads costs almost nothing, at three reads a position as much
  as at a hundred.** On the two-peaked population it is within 0.001 of 1.000 at every depth; on
  the rare-allele one it is between 1.6% and 4.3% low, and the four depths do not order themselves
  — 4.3% at three reads, 4.3% at eight, 1.6% at twenty, 3.5% at a hundred — so what is left there
  is not a depth effect.
- **At one individual the cost is worst at three reads and small by twenty** — 22% high, then 2%
  high, then 2% low — and the direction differs by population at three reads, which a systematic
  pull would not do. **"Small" is not "gone":** the 20- and 100-read cells are 0.977 ± 0.001 and
  0.973 ± 0.001, which are 23 and 27 standard errors from 1.000. At one sample there is no
  cross-sample coupling, so §4.1's variance term is **exact** there and cannot be the cause.
  **A 2 to 3% shortfall at one sample survives to any depth and this report does not explain it.**
- **Ten individuals is the worst panel size on both populations**, worse than 63 at every depth on
  both, and worse than one individual at every depth except three reads on the rare-allele
  population. **The shortfall to be explained is 9.7 points on the rare-allele population and 24.7
  on the two-peaked one, at three reads.**

  **It is not the variance term's residual, and an earlier draft of this report said it was.** Two
  reasons, either of which is enough. *Sign*: the exact quantity is
  `Var(k) = Σᵢ Varᵢ + ΣΣ_{i≠j} Covᵢⱼ`, the samples at a position are coupled through the frequency
  they share so the covariance is **positive**, and the estimator subtracts only the first sum —
  subtracting too little, which pushes the ratio **above** 1.000. These cells are below it. *Size*:
  the whole variance term at ten samples and three reads is the gap between the two posterior
  columns, 1.6 points on the rare-allele population and 2.2 on the two-peaked one. Setting it to
  anything at all moves those cells by at most 2 points against a shortfall of 10 and 25.

  **What does fit is §5's finding**, and ten samples at three reads is exactly the cell where §5
  measures it: the population outside the fit's family loses 16 points more than a matched
  population inside it. That accounts for the two-peaked column and not for the rare-allele one,
  **which is left unexplained here.**

---

## 5. Question 3 — does the fit's own population curve pull the answer toward itself?

**This is the question the plan names as the one that can kill the proposal.** The per-position
posteriors are computed under the fit's description of the population — one Beta over what
segregates plus two end masses — so a moment computed from them inherits a pull toward that
description. **On a population drawn from that family the pull is invisible**, which is why the
sweep carries a population outside it.

**The comparison needs a control, and the first draft of this measurement did not have one.** The
two-peaked population differs from the rare-allele one in two ways at once: it is outside the fit's
family, *and* its alternative alleles sit at frequencies of 0.3 and 0.9, where a handful of reads
settles a genotype easily, rather than piled up near zero where they do not. A gap between those
two could be either. The third population keeps the frequencies and drops the second peak — one
`Beta(30, 70)`, squarely inside the family, at a heterozygosity within 5% of the two-peaked one's.
**A pull toward the fit's own family has to show as a gap between those two.**

**Heterozygosity from the posteriors divided by heterozygosity from the same positions' own
genotypes**, the two-peaked population against its matched single-Beta control, and the gap between
them:

| individuals | depth | two peaks (outside) | one peak (inside) | gap |
|---:|---:|---:|---:|---:|
| 1 | 3 | 0.916 ± 0.057 | 0.944 ± 0.057 | −0.028 ± 0.081 |
| 1 | 8 | 0.929 ± 0.019 | 0.932 ± 0.028 | −0.003 ± 0.034 |
| 1 | 20 | 0.988 ± 0.003 | 0.983 ± 0.001 | +0.005 ± 0.003 |
| 1 | 100 | 0.983 ± 0.001 | 0.983 ± 0.000 | 0.000 ± 0.001 |
| 10 | 3 | 0.753 ± 0.036 | 0.912 ± 0.067 | **−0.159 ± 0.076** |
| 10 | 8 | 0.862 ± 0.015 | 0.945 ± 0.009 | **−0.083 ± 0.018** |
| 10 | 20 | 0.923 ± 0.016 | 0.920 ± 0.028 | +0.003 ± 0.032 |
| 10 | 100 | 0.917 ± 0.008 | 0.948 ± 0.013 | −0.031 ± 0.015 |
| 63 | 3 | 1.001 ± 0.002 | 0.998 ± 0.004 | +0.003 ± 0.004 |
| 63 | 8 | 1.000 ± 0.001 | 0.999 ± 0.001 | +0.001 ± 0.001 |
| 63 | 20 | 1.000 ± 0.000 | 1.000 ± 0.000 | 0.000 ± 0.000 |
| 63 | 100 | 0.999 ± 0.001 | 1.000 ± 0.000 | −0.001 ± 0.001 |

**The pull is confined to one corner, and it falls from 16 points at three reads to 8 at eight to
nothing at twenty — then returns at 3 points at a hundred, which nothing here explains.**

- **At ten individuals it is there and it is not small.** At three reads a position the population
  the fit's family cannot hold loses **16 percentage points more** of its heterozygosity than the
  matched population it can hold — 0.753 against 0.912 — and at eight reads, 8 points. Those are
  2.1 and 4.6 times their own uncertainty.
- **By twenty reads a position it is gone**: +0.003 ± 0.032. The −0.031 ± 0.015 at a hundred reads
  does not fit a monotone story with the zero at twenty, and with three replicates a cell it is
  the size of what these runs return by chance; **it is reported, not explained.**
- **At 63 individuals there is nothing to find at any depth.** Both populations come back within
  0.003 of their own genotypes' answer, at three reads as much as at a hundred.
- **At one individual there is nothing to find either**, but the run's precision there is ±0.08 at
  three reads, so what this shows is that the pull is smaller than eight parts in a hundred, not
  that it is absent.

**A perfect control does not exist, and this one is not it.** Being outside the fit's family *is* a
difference in where the alternative alleles sit — the second peak at frequency 0.9 is what puts the
population outside, so it cannot be removed while keeping the population outside. The two match on
heterozygosity to within 5% (6.89 against 7.24 differences per kilobase) and on where the main peak
sits, and they do not match on the mean alternative-allele frequency: 0.0106 against 0.00622,
because 30 positions in 100 of the two-peaked one's segregating set sit at frequency 0.9.

**Two reasons the gap is still better explained by the family than by that difference**, and both
are arguments rather than further measurements:

- **A position at frequency 0.9 is easier to read at three reads, not harder.** Most samples there
  carry two copies of the alternative and all three of their reads say so; the ambiguous case is a
  heterozygote, and the two populations have similar numbers of those — that is what matching on
  heterozygosity means. Read difficulty would push the two-peaked population's ratio **up**, and
  the measured gap is downward.
- **The fitted curve loses on the same population by the same order.** At ten individuals and three
  reads the curve's own heterozygosity comes back 0.740 of the genotypes' on the two-peaked
  population and 0.902 on the control (the `het: today /gt` column of §4's sweep, which the recorded output carries in full) — so the fit
  is genuinely describing the two-peaked population worse, which is the mechanism the pull is named
  for, seen from the other side.

**So the plan's kill condition is not met, and the reason is not that the pull vanishes
monotonically.** The condition reads: *"a pull toward the EM's Beta that does not shrink with
depth."* From three reads to twenty it shrinks by every measure — 16 points, then 8, then nothing —
and **at 63 individuals it does not exist at any depth**, which is the stronger of the two facts
because 63 is a cohort size at which a run has a panel at all. What does not fit a clean monotone
story is the −3.1 ± 1.5 point cell at a hundred reads, two standard errors from zero with a zero
beside it at twenty. **Three replicates a cell is not enough to say whether that is a real residual
or a draw**, and this report does not claim to know.

**What must be said in the honest form the plan asks for**, though: at ten samples and three reads
a position, the direct estimate is **unbiased given the read model rather than unbiased outright**.
The claim that survives without qualification is the one about genotypes (§3), and the claim about
posteriors carries the family the fit assumes with it.

---

## 6. Question 4 — what do mismapped positions do to the diversity?

A position where two stretches of genome the reference holds once both pile their reads up reads
part non-reference **in every sample**, which is heterozygosity's own signature. The joint fit
produces a per-position posterior that a position is that (`JointFit::noisy_posterior`).

**Planted at one position in a hundred**, with reads there disagreeing with the reference at 6 in
100 instead of 2 in 1,000 — the rate `examples/ng_joint_sample_count_sweep.rs` plants. The
heterozygosity is then computed twice, and **both are divided by the same panel's own genotypes**,
so these rows compare directly with §4's no-plant sweep:

| population, depth | individuals | over all positions | weighted by one minus the mismapped posterior | control |
|---|---:|---:|---:|---:|
| rare alleles, 3 reads | 1 | 1.207 ± 0.235 | 1.176 ± 0.264 | 1.011 |
| | 10 | 0.896 ± 0.021 | 0.880 ± 0.020 | 1.014 |
| | 63 | 0.957 ± 0.019 | 0.941 ± 0.028 | 1.015 |
| rare alleles, 20 reads | 1 | 0.976 ± 0.001 | 0.975 ± 0.003 | 0.990 |
| | 10 | 0.950 ± 0.005 | 0.942 ± 0.010 | 1.007 |
| | 63 | 0.984 ± 0.004 | 0.976 ± 0.015 | 1.009 |
| two peaks, 3 reads | 1 | 0.874 ± 0.107 | 0.815 ± 0.132 | 1.006 |
| | 10 | 0.754 ± 0.031 | 0.746 ± 0.035 | 1.000 |
| | 63 | 1.002 ± 0.002 | 1.001 ± 0.006 | 1.000 |
| two peaks, 20 reads | 1 | 0.990 ± 0.007 | 0.985 ± 0.008 | 0.989 |
| | 10 | 0.924 ± 0.016 | 0.929 ± 0.014 | 0.993 |
| | 63 | 1.000 ± 0.000 | 1.006 ± 0.004 | 0.994 |

**The control column is the genotypes over every position against the genotypes over the unplanted
ones**, and it sits between 0.989 and 1.015 throughout — so planting changed what the reads say and
not what is there, which is what makes the other two columns readable.

**Two things, and the first is a null:**

- **Planting barely moves the unweighted estimate.** Cell for cell against §4's no-plant sweep on
  the same populations at the same depths, the largest change anywhere is **4.2 percentage points**
  — the two-peaked population at one sample and three reads, 0.916 without plants against 0.874
  with them — and nine of the twelve cells move by 0.2 points or less. One position in a hundred,
  reading non-reference at 6 in 100 where a heterozygote reads about 500 in 1,000, is not enough to
  register.
- **The weighting is worse in nine of the twelve cells and better in three**, by margins smaller
  than the error bars in every one of them. It is not doing anything.

**⚠ An earlier version of this table divided by the population's own heterozygosity over unplanted
positions**, which is not a panel quantity — so the ratio kept the inbreeding factor of §3
uncancelled and the one-sample rows were 15% low before mismapping did anything. **Every number
above was 0.85 times its correct value at one sample.** The error was found by review, not by the
run.

### 6.1 A harder plant, because the first one was mild

At 6 disagreeing reads in 100 and three reads a position, most planted positions show no
disagreeing read at all, so §6's table is a weak test of anything. **Re-run on the two-peaked
population with the planted positions disagreeing at 25 reads in 100** — four times as strong, and
inside the ceiling of 0.45 the fit's own error-rate maximisation allows:

| depth | individuals | over all positions | weighted by one minus the mismapped posterior | control | §4's no-plant cell |
|---:|---:|---:|---:|---:|---:|
| 3 | 1 | 0.842 ± 0.068 | **0.674 ± 0.046** | 1.006 | 0.916 |
| 3 | 10 | 0.747 ± 0.034 | 0.743 ± 0.039 | 1.000 | 0.753 |
| 3 | 63 | 0.970 ± 0.008 | 0.970 ± 0.004 | 1.000 | 1.001 |
| 20 | 1 | 1.004 ± 0.012 | 0.983 ± 0.016 | 0.989 | 0.988 |
| 20 | 10 | 0.924 ± 0.014 | 0.930 ± 0.013 | 0.993 | 0.923 |
| 20 | 63 | 1.000 ± 0.000 | 1.006 ± 0.004 | 0.994 | 1.000 |

**Even at four times the strength the plant costs almost nothing, and the weighting costs a lot.**
The unweighted estimate differs from the same cell without plants by at most **7.4 percentage
points**, at one sample and three reads (0.842 against 0.916), and by 3.1 or less everywhere else.
**The weighting takes that worst cell 17 further points down**, 0.842 to 0.674 — and it is the only
change anywhere in the table larger than an error bar.

**So the answer to the plan's question is the opposite of the leaning it recorded.** The plan
expected the unweighted estimate to be badly high and the weighted one to rescue it. What happens
is that the fit's mismapped class is doing the work already — a position whose reads look mismapped
has its genotype posteriors computed under an error rate that explains them without a heterozygote
— and applying the same posterior a second time removes real variation along with artefact.
**Recommendation: compute both moments over every kept position, unweighted.**

Across the eighteen cells of §6 and §6.1 together, the weighting is worse in thirteen, better in
four and unchanged in one, and its single effect larger than an error bar is a loss.

**What this does not settle.** The planted positions here are ordinary positions with a raised
error rate, which is what the fit's mismapped class models
(`JointFit::noisy_posterior`). A duplication where two copies pile up and about half the reads
disagree in every sample is a different shape, and the fit has a separate class for it
(`JointFit::duplicated_posterior`); nothing here tests weighting by that one.

---

## 7. Question 5 — is the answer usable at one sample, not merely unbiased?

**Unbiased is not enough.** A prior seeded from a number with enormous spread is worse than one
seeded from a constant, and the constant a run falls back to when nothing could be fitted is
`ExpectedHeterozygosity::SPECIES_FALLBACK`, one difference per thousand bases.

**Measured at the shipped census size** — one genome at 2,000,000 census positions
(`parameter_prepass_census_sites.md` §5), the positions redrawn every run, so the spread carries
both the census's sampling and the genome's, over 40 runs a cell:

| population | F | heterozygosity: off the population | spread over runs | 19 runs in 20 land inside | the constant is off by | mean frequency: off the population | spread over runs | 19 runs in 20 land inside |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| rare alleles, tomato-like | 0 | +0.93% | 2.67% | ±5.2% | **+65%** | +0.27% | 1.45% | ±2.8% |
| | 0.8 | −0.71% | 8.30% | ±16.3% | **+65%** | −0.24% | 1.30% | ±2.5% |
| reference base rare | 0 | +0.03% | 2.09% | ±4.1% | +15% | −0.00% | 0.77% | ±1.5% |
| | 0.8 | +0.16% | 4.69% | ±9.2% | +15% | +0.15% | 0.91% | ±1.8% |
| two peaks | 0 | −0.23% | 1.83% | ±3.6% | −27% | −0.43% | 1.20% | ±2.4% |
| | 0.8 | +0.68% | 4.57% | ±9.0% | −27% | −0.15% | 1.58% | ±3.1% |
| a lump | 0 | −0.33% | 1.72% | ±3.4% | −34% | −0.28% | 1.10% | ±2.2% |
| | 0.8 | +0.14% | 4.06% | ±8.0% | −34% | −0.54% | 1.62% | ±3.2% |

**On seven of the eight rows one genome beats the constant by three-fold to twelve-fold**, and the
margin is widest where the constant is worst: on the tomato-like population the constant is 65%
high and 19 runs in 20 of the measurement land within 5.2% (random mating) or 16.3% (`F = 0.8`,
with `F` supplied).

**The eighth row is the one to watch.** On the population where the reference base is the rare one
at `F = 0.8`, the constant is 15% off and the measurement's 95% band is ±9.2% — a margin of 1.6,
which linkage alone could erase (see below). **Where the population happens to sit near one
difference per thousand bases, measuring it is worth much less**, and a run should be able to say
so from the spread it prints rather than leaving the user to know it.

**The mean frequency is the tighter of the two and needs no inbreeding coefficient**, so at one
genome it is usable without qualification: 19 runs in 20 land within 1.5% to 3.2% of the population
on all four shapes at both coefficients. **That is what lets the seed drop its blend toward a
neutral shape** — the blend exists to damp a small panel's noisy shape estimate, and at a
full-size census there is no such noise to damp.

**Two things this number is not.**

- **It assumes `F` is known.** The `F = 0.8` rows divide by `1 − F`, which multiplies the spread by
  five along with the estimate. Without `F` the answer is not wide, it is wrong by a factor of five
  (§3.4).
- **It assumes census positions are independent, and they are not.** The cited factor is on the
  **interval**, not on the variance: *"if a 100 kb stretch behaves as one draw, two million
  positions carry 8,000 draws and every such interval widens sixteen-fold; the true factor lies
  between about 3 and 16"* (`parameter_prepass_census_sites.md` §5), and `√(2,000,000 / 8,000)` is
  15.8. **At the top of that range the `F = 0.8` bands become ±261%, ±147%, ±144% and ±128%**, and
  at the bottom ±49%, ±28%, ±27% and ±24%. **At 16× the measurement is worse than the constant on
  all four populations at both coefficients; at 3× it still beats it on the tomato-like one.**
  Linkage does not touch the bias.

  **So how far one genome's measurement beats the constant is not settled by this report**, and it
  turns on a number nobody here measured: how far linkage reaches in the panel. What follows
  regardless is that **a run must print the spread it is claiming**, because the same estimate is
  worth a great deal on one population and nothing on another.

---

## 8. Question 6 — does any of it move a genotype?

**Seed the caller both ways and compare the calls.** At each cell of the read sweep, one sample's
genotype is called at 200,000 freshly drawn loci under the seed today's projection-and-search
returns and under the seed the direct moments return, with the caller's own genotype prior — the one that ships, which averages over the locus's unknown
allele frequencies rather than fixing them at an estimate (`MarginalizedDirichletPrior`) — and a
straightforward read likelihood. Counted over the loci that
segregate, since a locus where the population carries one allele is called the same way by
anything.

**The two seeds differ by up to 41% in the alternative concentration.** That is the second of the
prior's two numbers: the seed is a pair, a weight of belief on the reference allele and one shared
out over the alternatives, and it is the alternative half that carries how often a variant is
expected. The widest gap in the sweep is 1.32 × 10⁻² against 9.4 × 10⁻³, at 63 individuals on the
single-peak population.

**Calls that come out differently: 0.00% in 28 of the 36 cells.** The largest anywhere is
0.44 ± 0.44% of segregating loci; the next three are 0.17%, 0.14% and 0.13%. **The loci are redrawn
for every replicate and every panel size**, so those spreads are over independent draws — an
earlier version seeded them by depth alone and called the same 200,000 loci in every cell.

**A null is only worth reading if the comparison could have seen something, so every cell also runs
a control**: the same comparison with the direct seed's alternative concentration multiplied by
three.

| depth | control: calls moved by trebling the concentration | the two seeds, same cells |
|---:|---|---|
| 3 reads | 1.24% to 3.61% of segregating loci, non-zero in all nine | 0.00% in eight, 0.44% in one |
| 8 reads | 0.22% to 0.95%, non-zero in all nine | 0.00% in three, 0.02% to 0.17% in six |
| 20 reads | 0.00% in five cells, 0.01% or 0.02% in four | 0.00% in eight, 0.01% in one |
| 100 reads | 0.00% in every cell | 0.00% in every cell |

**So the null means one thing at three and eight reads and much less at twenty and a hundred.** At
three and eight reads the control moves calls in all eighteen cells, so the comparison there has
power — and the two routes' seeds, 40% apart at the widest, move at most 0.44 in 100 against a
control that moves up to 3.6. **Seven of the eight non-zero cells sit at those two depths**,
including the largest, so it is not that nothing moves; it is that what moves is between a tenth and
a twentieth of what the control moves.

**At a hundred reads nothing moves under any prior tested**, so that row says only that the reads
have taken the decision away entirely. **At twenty reads the control keeps a little power in four
cells of nine, and the two seeds moved a call in one of them** — 0.01%, in the same cell where the
control managed 0.02%. So where any power survived at twenty reads, the two routes did differ, by
about half of what a three-fold prior change achieves.

**Why the null is a step function in depth rather than a small probability.** At a fixed depth the
call is decided by a threshold on how many of the reads carry the alternative allele, and that count
is a whole number. One alternative read is about **750 times** likelier under a heterozygote than
under a homozygous reference at an error rate of 2 in 1,000 — a log-odds of 6.6 — and each reference
read is worth about 0.7 the other way, so the reads' verdict moves in steps of roughly **7** as the
alternative count changes by one. Shifting the prior's alternative concentration by 40% moves its
contribution by `ln 1.4` = **0.34**.

**So the prior moves the threshold or it does not, and every locus at that depth moves together.**
That is why 28 cells are exactly 0.00% and not scattered around some small rate: a 0.34 shift
almost never crosses a whole-number boundary, and when it does — as at 63 individuals and three
reads — a whole class of loci crosses at once and the cell reads 0.41%. Trebling the concentration
shifts the contribution by 1.1 instead, which crosses a boundary often enough to show at three and
eight reads and never at twenty, where the reads have already decided every locus.

**An earlier draft of this section argued the null from a probability — "0.34 in steps of 7, so
about one time in twenty" — which predicts 5 calls in 100 and would have contradicted the
measurement three paragraphs above it.**

**The consequence for the recommendation.** The case for the direct estimate is **not** that it
calls genotypes better. On these draws it calls them identically. The case is that it reports the
population correctly, costs no search, and stops the answer depending on the panel — and **a run
that reports a diversity 20% wrong is wrong in its output whether or not a genotype moved.**

---

## 9. What the detour costs today, on the same four populations

**Today the caller does not read the mean frequency off the census.** It takes the fitted
population curve, evaluates it into the `2N + 1` allele-count classes a panel of `N` diploid
individuals has, and searches for the two-parameter pair whose own predicted classes best match
([`population_diversity.md`](../ng/spec/population_diversity.md) §3.2). The number the caller uses
is that pair's mean frequency, blended toward a neutral shape by how much the panel has earned
([`ordinary_site_seed.md`](../ng/spec/ordinary_site_seed.md) §4).

**Handed the population exactly** — no census sampling, no fitting error, which is the best case
for the current path rather than a fair fight — the pair's mean frequency, over the population's
own:

| individuals | tomato-like | reference base rare | two peaks | a lump |
|---:|---:|---:|---:|---:|
| 1 | 0.999× | 0.999× | 1.000× | 1.000× |
| 10 | 1.097× | 0.918× | 0.924× | 1.004× |
| 63 | **1.177×** | **0.876×** | **0.817×** | 1.030× |
| 200 | **1.217×** | **0.861×** | **0.787×** | 1.043× |

Three things the plan did not have:

1. **The tomato-like column reproduces the plan's own figures** — 1.18× at 63 and 1.22× at 200 —
   through the shipped projection rather than a copy of it, which is what makes the other three
   columns worth reading.
2. **The error's sign is not fixed.** On a population where the reference base is the rare one at
   the positions that vary, and on the two-peaked one, the same machinery comes back **low** by the
   same order: 0.79× at 200 individuals. So "the recovered frequency is 1.2× the truth" is a fact
   about one shape, and what is general is that **it is wrong by up to a fifth in whichever
   direction the shape dictates, and more the larger the panel.**
3. **At one individual the search itself is exact** — three classes, two free after
   normalisation, against two parameters — but **the seed the caller uses is not**, because the
   blend toward the neutral shape takes a fifth of the shape there. The seed's mean frequency at one
   individual is 0.816×, 0.622×, 0.893× and 0.916× of the population's on the four shapes. That
   loss belongs to the blend, not to the projection, and the direct estimate does not incur it.

**The direct estimate on the same panels**, scored the same way against the population, sits
between **0.956× and 1.014×** at every panel size from ten up, on all four populations and at both
inbreeding coefficients (§3, right-hand columns).

**That comparison is not like for like, and the asymmetry favours the detour.** The detour's column
above has no sampling in it at all — it is handed the population exactly. The direct estimate's
figures carry the census's own scatter, whose standard error at 20,000 positions is ±1.1% to ±3.5%
depending on the population, so several of them are within a standard error or two of 1.000 and
none is more than 4.5% away. **What is not within reach of chance is the detour's column**: 1.217×
and 0.787× at 200 individuals, on a population handed to it exactly, with no sampling that could
have produced them.

### 9.1 Is the bias the compression, or the detour? It is the detour

**The argument that started this work was that the fit describes the population with four numbers,
the prior holds two, and no map from four to two can avoid biasing the prior.** The first half is
true and the second does not follow, and separating them changes what should be built.

**Two numbers is a real loss of information, and it is not a bias.** A concentration pair is
equivalent to the first two moments of the frequency distribution it stands for: a pair of expected
frequency `f` and total `A` has `E[f] = f` and `E[2f(1−f)] = 2f(1−f)·A/(A+1)`, and those two
equations invert. **So the two numbers the prior wants *are* two integrals of the curve**, and both
have a closed form over the fitted density:

```text
E[2 f (1 − f)]  =  p_segregating · 2ab / ((a + b)(a + b + 1))     already in the code
E[f]            =  p_fixed_alt  +  p_segregating · a / (a + b)     one line, not in the code
```

That map reproduces the density's own first two moments **exactly, at any panel size**. What it
discards is the third moment and above — which is a modelling choice, and one measured to reach no
genotype: swapping the repeat path's frequency prior for a richer spectrum in July changed nothing
at any quality threshold
([`ssr_marg_sfs_genotype_prior_2026-07-09.md`](ssr_marg_sfs_genotype_prior_2026-07-09.md)).

**What the caller does instead is not that map.** It evaluates the curve into the panel's `2N + 1`
allele-count classes and searches for the pair whose own predicted classes best match. That is a
divergence fit over a histogram, not a moment match, so its answer depends on how many classes there
are — which is why the error grows with the panel and why the panel appears in a quantity about the
population.

**Measured, all three routes on the same fitted curves**, each divided by the same panel's own
genotypes. Mean allele frequency, range over the four depths:

| population | individuals | today: project and search | the curve, integrated | the census average |
|---|---:|---:|---:|---:|
| nearly all alleles rare | 1 | 0.830–0.901 | 0.995–1.021 | 0.994–1.020 |
| | 10 | 0.947–1.015 | 0.974–0.997 | 0.969–0.986 |
| | 63 | 1.040–1.120 | 0.993–1.028 | 0.985–0.996 |
| two peaks (outside the family) | 1 | 0.866–0.880 | 0.986–1.023 | 0.985–1.022 |
| | 10 | 0.754–0.852 | 0.894–0.998 | 0.878–0.964 |
| | 63 | **0.749–0.766** | 1.050–1.068 | **1.000–1.000** |
| one peak (inside it) | 1 | 0.955–1.027 | 0.975–1.036 | 0.974–1.035 |
| | 10 | 0.842–0.909 | 0.943–1.017 | 0.917–0.955 |
| | 63 | 0.831–0.899 | 0.990–1.037 | 0.998–1.000 |

**Three readings, in the order they matter:**

1. **The search is the worst of the three in 34 of the 36 cells**, and its worst cell is 25% low —
   two-peaked population, 63 individuals, where both alternatives are within 7%. **So the bias is
   the detour and not the compression**, and the four-into-two argument, while it identified a real
   defect, named the wrong cause.
2. **Integrating the curve is within 7% everywhere and within 3% in two thirds of the cells**, with
   no census pass, no per-position posteriors, and — because the curve's moments are population
   quantities with no individuals in them — **no inbreeding correction at any cohort size.** It is
   two lines against the machinery §4 and §3 describe.
3. **Only the census average converges.** At 63 individuals it is 1.000 to three decimals on both
   populations the curve's family cannot hold exactly, where the curve settles 5 to 7% off and stays
   there at every depth. **The curve's error is a modelling error: more data does not reduce it,
   because the curve is converging on the best-fitting member of its family rather than on the
   population.** The census average has no family to converge to.

**The heterozygosity says the same thing less sharply**, because there the caller already takes the
curve's own integral rather than the search's (`ordinary_site_seed.md` §3). Comparing the curve's
integral against the census average over the same 36 cells: the two are within 0.2% of each other at
one sample, the curve is about 2 points closer at ten, and **at 63 individuals the census average is
closer in 9 of 12 cells** — exactly 1.000 on the two-peaked population where the curve is 4% high.
Same pattern, same reason.

### 9.2 Does the choice between the two surviving routes move a genotype?

**Not above eight reads a position, and at three reads it can — by as much as a three-fold change in
the prior does.** The two seeds compared here are the integrated curve against the census average,
at the same 36 cells and by the same method as §8:

| depth | calls moved, curve against census | for scale: trebling the concentration |
|---:|---|---|
| 3 reads | 0.00% in eight of nine cells, **2.98%** in the ninth | 1.24% to 3.61% |
| 8 reads | 0.00% to 0.08% | 0.22% to 0.95% |
| 20 reads | 0.00% in every cell | 0.00% to 0.02% |
| 100 reads | 0.00% in every cell | 0.00% |

**The size of the gap does not predict whether calls move — the depth does.** The two routes'
alternative concentrations never differ by more than 14% and usually by under 5%; the cell that
moved 3 calls in 100 differs by **6.7%**, and a cell differing by **14%** moved none. That is §8's
step-function mechanism seen from the other side: a prior shift either crosses a whole-read
threshold at a given depth or it does not.

**At 63 samples and three reads — tomato's corner — the two routes move between 0 and 0.45 calls in
100** at segregating loci. The 2.98% cell is at ten samples and three reads.

**So the extrapolation an earlier draft of this report made — that the two routes are too close to
reach a genotype — was not safe**, and the measurement replaces it.

**What this does to the recommendation is in §1.**

---

## 10. How these measurements were checked

**This project's failures are in measurements that were subtly wrong and invented a finding, not in
arithmetic.** Four checks were run against that, and each is a mutation of the source with the
result recorded — the source was restored immediately and the restore verified by comparison
against a saved copy.

| what was changed | what should have caught it | what it did |
|---|---|---|
| the estimator's denominator, `2N(2N − 1)` → `2N·2N` | the check at one known frequency, which shares no code with the sweep | ratios of 0.499, 0.900 and 0.990 at 1, 5 and 50 individuals instead of 1.000 |
| the inbreeding factor, `1 − F/(2N − 1)` → `1 − F/2N` | the same check | 0.336 instead of 1.007 at one individual, and 0.991 instead of 1.001 at five — the `F = 0` rows and the 50-individual row are the only ones it leaves alone |
| **every population's two Beta parameters swapped** | anything | **every heterozygosity was unchanged**, and every mean frequency moved — 0.00167 to 0.00433 on the tomato-like population |
| the variance term's sign in the read harness, `− Var(k)` → `+ Var(k)` | the plain-substitution column beside it | the corrected column returned 2.602 where the uncorrected one returned 1.666, which is impossible for a subtraction; the direct seed's alternative concentration went from 1.5 × 10⁻³ to 4.8 × 10⁻² and the share of genotype calls that move went from 0.00% to 5.04% |

**The third row is the one worth carrying forward.** `E[2 f (1 − f)]` under `Beta(a, b)` is
`2ab / ((a + b)(a + b + 1))`, which is **symmetric in `a` and `b`** — so a fixture set checked only
on its diversity cannot detect a transposed pair, however many populations it holds. What catches it
is the mean frequency, `a / (a + b)`, and only because none of the four populations is symmetric
about a frequency of one half. **A future population added to these programs must be checked for
that**, and the check is one line: does swapping its two shape parameters move its mean frequency?

**Four defects the review of this report found, and none of them was arithmetic:**

| what was wrong | how it read | what it actually was |
|---|---|---|
| **§6's mismapped tables divided a panel estimator by a population quantity** | one sample looked 15% worse under mismapping than it was | the inbreeding factor of §3 does not cancel between those two, so every one-sample cell was 0.85 times its correct value. Both columns now divide by the same panel's own genotypes, and §6 was rebuilt. |
| **§4.2 attributed the ten-sample shortfall to the variance term's residual** | a plausible mechanism, hedged as "an explanation and not a measurement" | the residual has the **opposite sign** and is **ten times too small** to account for it. Hedging the confidence did not make the direction less wrong. |
| **§8 argued its own null from a probability** | "0.34 in steps of 7, so about one in twenty" | that predicts 5 calls in 100 against a measurement of 0.44 in the worst cell. The mechanism is a step function, not a probability. |
| **§7 applied the linkage factor to the variance where its source applies it to the interval** | a four-fold understatement of how much linkage widens the one-sample band | `√(2,000,000 / 8,000)` is 15.8, and it multiplies the interval. |

**Two of those four are the same failure**: a mechanism asserted because it was the expected answer,
in a place where the numbers beside it said otherwise. Every counted figure in the draft the review
read was correct.

**Two defects these programs had, and how they surfaced:**

- One population was written with `0.9900` of positions invariant and `0.0250` fixed
  non-reference — 1.015 between them. Nothing crashed: the segregating share clamps at zero, so the
  program printed a population with a **negative** heterozygosity of −0.0033 and a drawer that never
  segregated. It is now refused by an assertion that names the population and both shares.
- The first read sweep ran one drawn cohort a cell and printed ratios with no spread beside them.
  Two of the readings it appeared to support — a depth trend in the read model's cost, and a gap
  between the two populations — did not survive replicates and error bars. **Every ratio in §4 to
  §6 is now three independently drawn cohorts with its own uncertainty printed**, and the
  conclusions that were drawn are only the ones larger than it.
- **Two generators were seeded so that things meant to be independent were not.** The two
  populations of §3.4 were keyed on their names' *lengths*, and "selfing" and "outbred" are both
  seven bytes, so they shared a stream — in a table whose whole reading is that two independent
  draws come back the same. And §8's genotype comparison was keyed on depth alone, so every
  replicate and panel size at one depth called the identical 200,000 loci. Both are fixed and both
  arms re-run; **neither changed a conclusion**, which is luck rather than diligence.

**What a reviewer should check first**, in order of how much rests on it:

1. **That §3.3's factor is derived correctly**, and that the check at one known frequency really
   shares no code with the sweep it is checking.
2. **That §5's argument survives its control being imperfect.** The control cannot match the
   population it controls in every respect — being outside the fit's family *is* a difference in
   where the alleles sit — and §5 says which respects it does and does not match. The two reasons
   it gives for preferring the family explanation are arguments, not measurements.
3. **That §8's control column licenses exactly the null beside it and no more** — it does not, at
   twenty and a hundred reads, and §8 says so.

**Four accidents the populations here share, each with the mutation that would exploit it**, found
by the review of this report and left as a warning to whoever extends these programs:

- **The inbreeding factor is checked against the model that generated it.** The drawer implements
  `F` as *one ancestral copy counted twice with probability `F`* and the correction is derived from
  exactly that. Replacing `F` by `√F` in **both** leaves every table byte-identical. So §3.3
  confirms the arithmetic, not the parameterisation.
- **Every panel is drawn with one `F` shared by all its individuals.** A rule that weighted the
  panel's coefficient by each sample's covered positions, instead of taking the plain mean, would
  pass every arm here and be wrong on a real panel mixing selfed and outcrossed accessions. (The
  plain mean is the correct rule: with per-individual `Fᵢ` the estimator's expectation is
  `π(1 − F̄/(2N − 1))` with `F̄` unweighted.)
- **A planted mismapped position carries real variation.** Deleting the guard that restricts the
  control column to unplanted positions moves every cell of §6 by less than its error bar, because
  the planted positions are a random one in a hundred drawn from the same population.
- **§6's mild plant cannot tell "the weighting does nothing" from "the weighting is not applied".**
  Replacing the per-position weight by a constant would make the weighted column identical to the
  unweighted one, which supports §6's conclusion more cleanly than the real measurement does.
  **Only §6.1's harder plant catches it.**

---

## 11. Reproducing every number here

All three programs are deterministic: the same command gives the same digits. **The one exception
is the `seconds` column**, which is wall clock — re-running the harder-plant configuration after an
unrelated edit reproduced every measured value to the last digit and moved only that column.

```text
./scripts/dev.sh cargo run --release --example ng_prior_moment_estimators
    §2's population table, §3's sweep, §7's one-genome table, §9's current-path table.
    20,000 positions and 24 replicates by default, and 2,000,000 positions with 40 replicates
    for the one-genome arm; 22 s.

./scripts/dev.sh cargo run --release --example ng_prior_moments_from_reads -- 10000 0.15 0.01 3 0,1 mismapped
    §4, §6 and §8, on the first two populations. 37 minutes.

./scripts/dev.sh cargo run --release --example ng_prior_moments_from_reads -- 10000 0.15 0.01 3 2 nomismapped
    §5's control population. 10 minutes.

./scripts/dev.sh cargo run --release --example ng_prior_moments_from_reads -- 10000 0.15 0.01 3 1 mismappedonly 0.25
    §6.1's harder plant. 6 minutes.

./scripts/dev.sh cargo run --release --example ng_prior_moments_from_reads -- 10000 0.15 0.01 3 0,1 nomismapped
./scripts/dev.sh cargo run --release --example ng_prior_moments_from_reads -- 10000 0.15 0.01 3 2 nomismapped
    §9.1's three-route comparison, which needs the column the two runs above predate. 35 minutes.

./scripts/dev.sh cargo run --release --example ng_prior_moment_one_sample_inbreeding
    §3.4 and §3.5. 200,000 positions for the two-population construction and 50,000 for the
    cohort-size sweep; 14 minutes.
```

**Which recorded output each table came from**, since three of these programs were extended while
the report was being written and an earlier draft quoted a superseded run:

```text
§2, §3, §9                estimators_run5.md
§7                        estimators_run5.md   (the mean-frequency columns exist only there)
§4, §5, §8                reads_v2_main.md and reads_v2_control.md
§6, §6.1                  reads_v2_main.md and reads_v2_hardplant.md
§9.1                      reads_v3_main.md and reads_v3_control.md
§3.4, §3.5                one_sample_run3.md
```

The arguments are, in order: census positions, the homozygote excess every drawn individual
carries, the share of positions planted mismapped, how many cohorts to draw a cell, which
populations to run, which halves to run (`mismapped`, `nomismapped` or `mismappedonly`), and how
often a read at a planted position disagrees with the reference. **Naming the populations is not
optional in the first command**: omitting it runs all three, and the recorded output holds two.

**And the one that is not a limit but a gap in the runs above:** `ng_prior_moment_one_sample_inbreeding.rs`
now carries a cohort-size sweep as well as the two-population construction, so its name is narrower
than its contents.

**Validation on this branch**, all through `./scripts/dev.sh`: `cargo test --lib` gives 4,874
passed / 0 failed / 14 ignored, matching the branch point; `cargo fmt --all -- --check` and
`cargo clippy --all-targets --all-features -- -D warnings` both exit 0; `cargo doc --no-deps --lib`
reports the same 27 unresolved intra-doc links it reported before this work, all pre-existing.

---

## 12. What this did not measure

- **A real panel.** Every number is a drawn cohort. Confirmation on the tomato census is
  [`parameter_prepass_cohort.md`](../ng/spec/parameter_prepass_cohort.md) §10's third question and
  stays open.
- **Linkage between census positions.** Positions are drawn independently here. Real census
  positions are scattered but correlated, which widens every spread quoted and leaves every bias
  alone ([`parameter_prepass_census_sites.md`](../ng/spec/parameter_prepass_census_sites.md) §1).
  **So every spread in this report is optimistic by a factor `parameter_prepass_census_sites.md` §5 puts between 3 and 16**, and
  every bias is not.
- **A census as thin, in segregating positions, as a real one is rich.** §4 to §8 give the fit
  about 200 segregating positions a cohort against a real run's ten thousand, because every arm
  there is a full joint fit. **The direction that biases is not obvious and was not measured**: a
  thinner census makes the fitted curve looser, which could make the pull §5 measures either
  larger (less data to resist the family) or smaller (a looser curve pulls less hard). §3, §7 and
  §9 do not have this limit — they run at 20,000 and 2,000,000 positions and involve no fit.
- **The runs-of-homozygosity estimator, which §3.5 now makes the preferred source of the inbreeding
  coefficient.** Nothing here ran it, and every genome drawn in this report has its homozygosity
  scattered evenly across positions rather than gathered into runs — so these programs could not
  have tested it even if they had called it. **Two things would have to be measured before the
  recommendation rests on it**: whether it recovers the coefficient at the depths this caller
  commits to, and whether a census of scattered positions is a substrate it can work on at all,
  since it was built to walk genome windows.
- **Whether the exact `Var(k)` closes §4's residual.** The variance used is the sum of the samples'
  own and ignores their coupling at a position. The exact quantity is reachable inside the fit's own
  integration over frequency; nothing here computed it.
- **The repeat-tract prior**, which reads a fitted length spectrum per stratum and takes no such
  detour.
- **More than one alternative allele at a position.** Every drawn position here is biallelic, so
  nothing here says how the alternative concentration should be shared out where a position carries
  two different alternative alleles.
- **Depths above 124 reads a position.** The harnesses subsample each position down to 124 reads
  before recording it, so the 100-read arm is the deepest this census can carry. §13 says why that
  is not merely a choice.

---

## 13. One defect found on the way, in code this work did not touch

**A census recorded with a depth cap of 255 makes the joint fit panic at any position holding more
than 124 reads.** Hit while building the read sweep, at 100 reads a position:

```text
thread panicked at src/ng/parameter_estimation/joint/fit.rs:903:
a stored depth code stands for 35 depths, more than the 32 this fit reserves room for
```

**The mechanism.** The census depth ladder is exact to 124 reads a position and widens above it;
its first widening bin stands for 35 different depths (`generic::depth_bins`). The fit reserves
room for 32 (`MAX_RECORDED_SPREAD`) and asserts rather than truncating. The census's own cap
narrows a bin's range to the cap before the fit sees it
(`DepthCap::denominator_for`), so **a cap at or below 124 makes every bin exact and the assertion
unreachable** — and `DepthCap::MAX` is 255, which does not.

**Why it has not bitten.** Nothing in `src/ng/run/` or `main.rs` constructs a `DepthCap` yet: the
census writer takes one as an argument and no run wires it. `examples/ng_joint_sample_count_sweep.rs`
passes `DepthCap::MAX` and draws at 8 reads a position, so it never reaches the ladder's widening
half.

**What it costs whoever wires the flag.** A user who sets a large cap on a deep run gets a panic
rather than a number, and the message names a code width rather than the cap they set. **Two
constants in two modules have to agree and nothing says so** — not the ladder's documentation, not
`MAX_RECORDED_SPREAD`'s, and not `DepthCap`'s, whose own comment explains only why the cap cannot
exceed 255. **Left as found**: this report's scope is the genotype prior's two numbers, and the fix
is a decision about which of the two constants gives way.
