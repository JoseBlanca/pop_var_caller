# An answer to the plan: the detour is the defect, not the compression

**Status:** findings, 2026-08-27. Written for whoever wrote
[`ordinary_site_prior_moments.md`](ordinary_site_prior_moments.md), which this answers. The
measurements and every number quoted here are in
[`../../reports/ng_ordinary_site_prior_moments_2026-08-27.md`](../../reports/ng_ordinary_site_prior_moments_2026-08-27.md);
the implementation this recommends is
[`../spec/ordinary_site_prior_moments.md`](../spec/ordinary_site_prior_moments.md).

**In one paragraph.** The plan is right that the caller's mean allele frequency is badly wrong and
right that the error grows with the panel — measured at **1.22× the truth at 200 individuals** on a
tomato-like population and **0.75× at 63** on one with two frequency peaks. It is wrong about the
cause. The defect is not that four numbers are being compressed into two; it is *how* the two are
chosen. **There is a map from the four to the two that is exact at every panel size, it is two
lines of arithmetic, and it makes most of the plan's work unnecessary.** What remains of the plan is
a smaller and still real advantage, and it now has to justify itself against that two-line fix
rather than against the broken thing it was compared with.

---

## 1. Why "four numbers cannot go into two without bias" is wrong

The plan's §1 says:

> **Two parameters cannot follow a distribution with a spike at frequency zero**, so the recovered
> frequency comes back 1.18× the truth at 63 individuals and 1.22× at 200.

**The measurement is right and the sentence joining the two halves is not.** Compressing four
numbers into two loses information. **Losing information is not the same as biasing what is kept**,
and here what is kept can be kept exactly.

### 1.1 The prior's two numbers *are* two integrals of the curve

The prior holds a pair of concentrations, `(α_ref, α_alt)`. Written the other way round they are an
expected allele frequency and a total conviction, and that pair determines exactly two properties of
the frequency distribution it stands for:

```text
E[f]           =  α_alt / (α_ref + α_alt)                        write this f
E[2 f (1 − f)] =  2 f (1 − f) · A / (A + 1)      where A = α_ref + α_alt
```

**Those two equations invert.** Given a wanted mean frequency `f` and a wanted heterozygosity `π`,
the total that produces them is `A = π / (2f(1−f) − π)` — which is
`seed_generic::total_for_diversity`, already in the tree. **So a concentration pair and a pair of
moments are the same object in two notations.**

And the fitted curve's own two moments have closed forms:

```text
E[2 f (1 − f)]  =  p_segregating · 2ab / ((a + b)(a + b + 1))
E[f]            =  p_fixed_alt  +  p_segregating · a / (a + b)
```

The first is `FrequencyDensity::expected_heterozygosity` and **has been in the code all along**. The
second does not exist and is one line.

**So the map from four numbers to two is: evaluate two integrals, then invert the identity above.**
It reproduces the fitted curve's first two moments *exactly*, for any curve in the family, at any
panel size. There is no approximation in it to be biased.

### 1.2 The spike at zero is not an obstacle to that map

The plan's mechanism is that a two-parameter Beta "has no spike to offer against a distribution that
has two, so it compromises". **That is true of a curve-fitting exercise and irrelevant to a moment
match.** A point mass at frequency zero contributes exactly `0` to `E[f]` and exactly `0` to
`E[2f(1−f)]`. It is not something the two numbers have to *represent*; it is something they have to
*integrate over*, and integrating over it is trivial. The same goes for the mass at frequency one:
it contributes its own weight to `E[f]` and zero to the heterozygosity, and both are one term.

**What is genuinely lost is the third moment and above** — the *shape* of the curve beyond its mean
and spread. That loss is a modelling decision, and it was measured on this project a month before
the plan was written: swapping the repeat path's frequency prior for a richer spectrum changed
nothing at any quality threshold
([`../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md`](../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md)).
**The plan cites that null and explains it correctly** — the prior can hold only two moments, so the
curve's shape was never reaching a genotype. It then does not apply the same reasoning one step
earlier, to the question of how the two moments should be obtained.

---

## 2. Where the bias actually comes from

**The caller does not take the map of §1.1.** It does this
([`../spec/population_diversity.md`](../spec/population_diversity.md) §3.2):

1. **Evaluate the curve into the `2N + 1` allele-count classes a panel of `N` diploid individuals
   has.** This step is exact — it is a change of representation and nothing is fitted.
2. **Search for the concentration pair whose own predicted classes best match those classes**, by
   minimising a divergence over the `2N + 1` numbers.

**Step 2 is a curve fit in class space, not a moment match, and that is the whole defect.**
Minimising a divergence over a histogram is not the same objective as reproducing two integrals, so
it returns a different pair — and **the histogram has `2N + 1` entries, so the objective itself
changes with the panel.**

**That is why the error grows with the panel, and it is the signature to reason from.** A moment
match cannot depend on how many individuals the run holds, because neither moment mentions them. A
class-space fit must. So *the panel appearing in the answer at all* was already proof that the route
was not a moment match — which is the observation the plan's own §1 makes ("the panel size should
not appear in the answer at all, and today it does") without following it to its cause.

**The spike does make the class-space fit worse**, exactly as the plan says: the family cannot put
mass at class 0, so the fit trades the ends against the middle, and the more classes there are the
more the middle outvotes the ends. **But that is a property of the fitting criterion, not of the
two-number representation.**

### 2.1 The three routes, measured side by side

Mean allele frequency, each route divided by the same panel's own genotypes, on cohorts drawn from
three populations and fitted with the shipped joint fit. Range over four read depths from 3 to 100:

| population | individuals | today: project and search | integrate the curve | average the census |
|---|---:|---:|---:|---:|
| nearly all alternative alleles rare | 1 | 0.830–0.901 | 0.995–1.021 | 0.994–1.020 |
| | 10 | 0.947–1.015 | 0.974–0.997 | 0.969–0.986 |
| | 63 | 1.040–1.120 | 0.993–1.028 | 0.985–0.996 |
| two frequency peaks, outside the family | 1 | 0.866–0.880 | 0.986–1.023 | 0.985–1.022 |
| | 10 | 0.754–0.852 | 0.894–0.998 | 0.878–0.964 |
| | 63 | **0.749–0.766** | 1.050–1.068 | **1.000–1.000** |
| one peak, inside the family | 1 | 0.955–1.027 | 0.975–1.036 | 0.974–1.035 |
| | 10 | 0.842–0.909 | 0.943–1.017 | 0.917–0.955 |
| | 63 | 0.831–0.899 | 0.990–1.037 | 0.998–1.000 |

**The search is the worst of the three in 34 of the 36 cells**, and its worst is a quarter low. The
closed-form integral of the *same four numbers* is within 7% everywhere. **If the compression were
the problem, that middle column could not exist.**

---

## 3. What to do: two changes, and the first does not depend on the second

### 3.1 Step one — integrate, do not search. Two lines.

- Add `E[f]` to `FrequencyDensity` beside `expected_heterozygosity`.
- Delete the projection into allele-count classes and the search that reads it
  (`fit_spectrum_shape`, `fill_expected_spectrum`, `FittedSpectrum`,
  `FrequencyDensity::allele_count_classes`, and the adapter `FittedFrequencySpectrum` in
  `calling::run_parameters` that exists only to carry a projection between them).
- Solve the pair from the two moments by the identity already in `total_for_diversity`.

**Three things fall out for free**, and each is worth naming because each is a thing the plan
proposed to solve some other way:

- **The panel size stops appearing in the answer.** Not because it was corrected for, but because
  neither integral mentions it.
- **The blend toward the neutral shape can go too.** That blend exists to damp the *search's*
  unreliable shape at a small panel, and with no search there is nothing to damp. It is not free
  today: at one individual it puts the seed's mean frequency at **0.62× to 0.92×** of the truth
  across four populations, on populations where the search itself was exact.
- **No inbreeding coefficient is needed anywhere.** The curve's moments are properties of a
  *population*; there are no particular individuals in them to be inbred. §4 below is why that
  matters.

**The search also costs 11.8 minutes at 3,200 individuals**, by the code's own measurement, and that
goes with it.

### 3.2 Step two — average over the census positions, if it earns it

This is the plan's proposal, and **it survives the above with one real advantage left**: it is the
only route whose error goes to zero as the cohort grows. A curve converges on the best-fitting
member of its family; a census average has no family to converge to. The table in §2.1 shows it
directly — at 63 individuals on the two-peaked population the census average is **1.000 at every
depth** while the integrated curve settles at 1.05–1.07 and stays there however much data arrives.

**What it costs is three pieces of machinery step one does not need:**

1. **The fitting loop must accumulate two running sums.** It already computes, at each census
   position, how likely each sample is to be heterozygous and how likely both its copies are to be
   non-reference; today those are used and discarded. (Or store them — `genotype_posteriors` exists
   — at 12 bytes a position a sample, which is 1.5 GB for fifty samples over two million positions.)
2. **The heterozygosity needs a term the plan's formula does not have.** `k(2N − k)` is quadratic in
   the allele count, so substituting the posterior *mean* count is not the mean of the formula:
   `E[k(2N − k)] = 2N·E[k] − E[k]² − Var(k)`. **Dropping `Var(k)` returns 2.538 ± 0.165 times the
   truth at one sample and three reads a position**, against 1.219 ± 0.152 with it. At 63 samples the
   two agree to three decimals from eight reads up — so a test written on a cohort cannot catch its
   absence.
3. **The heterozygosity must be divided by `1 − F/(2N − 1)`.** §4.

---

## 4. The plan's heterozygosity estimator is not unbiased, and this is the thing it most needs to hear

The plan states, as settled:

> At `N = 1` the diversity estimator collapses to `k(2 − k)` … **A single genome estimating
> population diversity by counting its own heterozygous sites is the estimator, not an
> approximation of it.**

**That is true only for a population that mates at random.** Nei's estimator asks how often two
chromosomes drawn from the panel differ. One pair in `2N − 1` is the two copies inside one
individual, and in a self-pollinated plant those are the same ancestral copy with probability `F`.
So its expectation is the population's heterozygosity times `1 − F/(2N − 1)`.

Measured across four populations and nine panel sizes at `F = 0.8` — tomato's fitted range — with 21
of 36 cells within one standard error of that factor's prediction, 33 within two and all 36 within
three:

| individuals | 1 | 2 | 3 | 10 | 63 | 1000 |
|---|---:|---:|---:|---:|---:|---:|
| shortfall | **−80%** | −27% | −16% | −4% | −0.6% | −0.04% |

**A single self-pollinated accession reports one fifth of its population's diversity.** The
correction is that one factor and it is exact. **Step one of §3 needs none of this**, because the
curve is a population and not a panel.

**Where the coefficient should come from, if step two is built:** not from the joint fit's
homozygote excess. That quantity is measured against a diversity the same fit produced, so dividing
that diversity by `1 − F` closes the loop
[`../spec/parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md) §6.3 already warns
about. The runs-of-homozygosity estimator in `parameter_estimation::generic::runs` reads the
coefficient off *where heterozygotes sit along the genome*, needs no population expectation, works
from one genome, and §6.2 of that document already decided it is the one a caller reads.

---

## 5. How the prior is applied in the calling loop, and how much precision is worth buying

**The short answer: far less than this work assumed.** The plan's sixth question asks whether any of
this moves a genotype, and its leaning is that calls "may move at three reads". They do not move
much at all.

**Measured**: one sample's genotype called at 200,000 drawn loci under two seeds, with the caller's
own genotype prior (`MarginalizedDirichletPrior`) and a straightforward read likelihood, counted over
the loci that segregate.

- **Two seeds differing by up to 41% in the alternative concentration moved 0.00% of calls in 28 of
  36 cells**, and 0.44% in the worst.
- **A control that trebles the concentration** moved 1.24% to 3.61% of calls at three reads a
  position, 0.22% to 0.95% at eight, and **nothing at all at twenty and a hundred**.

**Whether a call moves is a threshold effect, not a proportional one**, and that is the single most
useful thing to know about applying this prior. At a fixed depth the call is decided by a threshold
on how many reads carry the alternative allele, and that count is a whole number. One alternative
read is worth a log-odds of about 6.6 in favour of a heterozygote at an error rate of 2 in 1,000, so
the reads' verdict moves in steps of roughly 7, while a 40% shift in the prior moves its
contribution by `ln 1.4` = 0.34. **So a prior change either moves the threshold — and then every
borderline locus at that depth moves together — or it moves nothing.** That is why 28 of 36 cells
read exactly `0.00%` rather than scattering around some small rate.

**Measured directly, and it is why the last bullet of the plan's sixth question cannot be answered
with a rule of thumb.** Comparing the two routes of §3 against each other — the integrated curve
against the census average — over the same 36 cells:

| depth | calls moved, integrated curve against census average | for scale: trebling the concentration |
|---:|---|---|
| 3 reads | 0.00% in eight of nine cells, **2.98%** in the ninth | 1.24% to 3.61% |
| 8 reads | 0.00% to 0.08% | 0.22% to 0.95% |
| 20 reads | 0.00% in every cell | 0.00% to 0.02% |
| 100 reads | 0.00% in every cell | 0.00% |

**The two routes' concentrations never differ by more than 14% and usually by under 5%** — and the
cell that moved 3 calls in 100 differs by 6.7%, while a cell differing by 14% moved none. **So the
size of the gap does not predict whether calls move; the depth does.**

**Three consequences for the loop:**

- **Above twenty reads a position the prior is not participating in the decision at all.** Nothing
  tested here moves a call there, including a three-fold change. Effort spent on these two numbers
  is wasted on any high-coverage run.
- **At three reads a position the prior does decide borderline loci**, and the choice between the
  two routes of §3 can move up to **3 calls in 100 at segregating loci** — the same order as a
  three-fold change in the prior. **At 63 samples and three reads, which is tomato's corner, it is
  between 0 and 0.45 in 100.**
- **Being wrong by a quarter is a different matter from being wrong by a twentieth.** The search is
  wrong by up to a quarter (§2.1), and that is what §3.1 removes for two lines. **The remaining
  choice between the curve and the census is worth deciding on what the run's reported diversity
  should be, and on low-coverage calling, not on high-coverage calling.**

**What this does not cover, and should not be read as covering:** this is a genotype call, not the
caller. It has no candidate-allele selection, no quality gate and no cohort step, and every locus
drawn is biallelic. It isolates the one thing two seeds differ in, which is what the question asks,
and it says nothing about how a wrong prior interacts with the rest of the loop.

---

## 6. What the plan got right

Stated plainly, because most of it did:

- **Both estimators are unbiased at every panel size on true genotypes**, which is what the plan's
  first question asks. Largest departure across four populations and nine panel sizes: 3.4% for the
  mean frequency and 4.5% for the heterozygosity, both at one individual and both inside their own
  uncertainty. From ten individuals up, 0.45% and 1.8%.
- **The census is not ascertained and needs no correction for it.** Confirmed by construction.
- **The kill condition in the plan's §5 is not met.** The fit's own curve does pull the census
  averages toward itself at ten samples and low depth — a population the family cannot hold loses 16
  percentage points more of its heterozygosity than a matched population it can, at three reads —
  but the pull is 8 points at eight reads and gone at twenty, and **at 63 individuals it does not
  exist at any depth**.
- **The insistence that the sweep carry a population the Beta family cannot hold.** Without it the
  most useful column of §2.1 above would not exist: on the two shapes inside the family, all three
  routes look tolerable at 63 individuals and the difference between them is invisible.
- **Holding the drawn positions fixed across the panel-size arms.** Kept, and it matters.

---

## 7. What the plan missed, beyond §1

- **The inbreeding factor** (§4), which is 80% at the corner this caller commits to.
- **The variance term** (§3.2, item 2), without which an implementation of the plan's own formula
  returns two and a half times the truth at one sample and three reads.
- **The fourth question's answer is the opposite of its leaning.** The plan expects the unweighted
  heterozygosity to be badly high with mismapped positions present, and the weighted one to rescue
  it. Planting one position in a hundred at four times the disagreement rate of the sibling harness,
  the unweighted estimate is at most **7.4 percentage points** from the same cohort with nothing
  planted — and the weighting is worse in thirteen of eighteen cells, its one effect larger than an
  error bar being a **17-point loss** at one sample and three reads. The fit already models the
  class; weighting by its posterior a second time removes real variation along with artefact.
- **At one sample the coefficient cannot be estimated by anything that reads census positions one at
  a time.** Two populations can be built whose single-genome census is drawn from the identical
  distribution while their diversities stand in the ratio `1 − F`; built at `F = 0.8`, one at 7.58
  differences per kilobase and one at 1.52, both showing a genome heterozygous at 0.00151515 of
  positions. The shipped fit returns a homozygote excess of **0.000** for both, at 3, 20 and 100
  reads alike. **From two samples up it returns 0.83 against a truth of 0.8** — the cliff is between
  one sample and two, and depth does nothing about it.

---

## 8. What is still unmeasured

- **Any real panel.** Everything is drawn cohorts; this checkout cannot rebuild the tomato census,
  which the plan itself records as a standing limit.
- **Linkage between census positions.** Every position here is drawn independently. Real ones are
  linked, which widens every spread quoted by a factor
  [`../spec/parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md) §5 puts
  between 3 and 16 — **on the interval, not the variance** — and leaves every bias alone.
- **The runs-of-homozygosity estimator**, which §4 makes the preferred source of the inbreeding
  coefficient. Nothing here ran it, and every genome drawn in this work has its homozygosity
  scattered evenly across positions rather than gathered into runs, so these programs could not have
  tested it even if they had called it.
- **Anything but a single-sample, biallelic genotype call.** §5's comparison is the shipped genotype
  prior against a straightforward read likelihood, and nothing else in the loop. What a wrong prior
  does once candidate-allele selection, the quality gate and the cohort step are in front of it is
  not measured here.
