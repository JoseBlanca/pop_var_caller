# Two numbers from the census, with no curve in between — a research plan

**Status:** research plan, 2026-08-27. **No code and no measurement yet — this says what to
measure and what would settle it.** Written from a design conversation with the owner; the
proposal is the owner's.

**⚠ Answered, 2026-08-27, and §1's explanation of the defect is wrong.** The bias is real and the
figures in §1 are reproduced, but it does not come from compressing four numbers into two — that
map exists, is exact at every panel size, and is two lines. It comes from choosing the two numbers
by a divergence fit over the panel's allele-count classes. **Read
[`ordinary_site_prior_moments_answer.md`](ordinary_site_prior_moments_answer.md) before acting on
this document**; the measurements are in
[`../../reports/ng_ordinary_site_prior_moments_2026-08-27.md`](../../reports/ng_ordinary_site_prior_moments_2026-08-27.md)
and the implementation in [`../spec/ordinary_site_prior_moments.md`](../spec/ordinary_site_prior_moments.md).

**The objective, in one sentence.** Estimate the two numbers the ordinary-site genotype prior
needs **directly from the census positions**, without bias, at every cohort size from one sample
to thousands and every depth from about three reads a position to hundreds.

**Companions:** [`../spec/calling_priors.md`](../spec/calling_priors.md) §4 (what the caller does
with the two numbers), [`../spec/population_diversity.md`](../spec/population_diversity.md) §3
(where they come from today), [`../spec/ordinary_site_seed.md`](../spec/ordinary_site_seed.md)
(the seam as it now stands), [`../spec/parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md)
§3 (how census positions are chosen).

---

## 1. What the prior needs, and what it is getting

**The prior holds exactly two numbers**, and no design has ever changed that: a concentration for
the reference allele and one for the alternative. Written the other way round, they are **an
expected allele frequency and a total conviction**, and the pair is fixed by the population's
**first two moments** — the mean alternative-allele frequency `E[f]`, and the heterozygosity
`E[2f(1−f)]`, which is the same information as the variance.

**What it is getting today is those two numbers taken by a detour.** The joint pre-pass fits a
four-parameter description of the population — a spike of mass at frequency exactly zero, a spike
at exactly one, and a Beta over what segregates. The caller then projects that onto the panel's
`2N + 1` allele-count classes and *searches* for the two-parameter pair whose own predicted classes
best match.

**The detour has a bias that grows with the panel, and the mechanism is not subtle.** Two
parameters are fitted to `2N + 1` numbers. At one individual there are three classes, two free
after normalisation, against two parameters — the fit is exact. At sixty-three there are 127, and
a two-parameter Beta has no spike to offer against a distribution that has two, so it compromises.

**Measured, on drawn densities, before any of it is repaired:** the recovered expected frequency
comes back **1.18× the truth at 63 individuals and 1.22× at 200** on a tomato-like shape
(`examples/ng_spectrum_panel_floor.rs`). The seed's implied diversity used to carry the same error
— 9.9% low at 63, 18.6% on a human-like shape — and **that half is already fixed**: the conviction
is now solved from the measured heterozygosity rather than taken from the search
([`../spec/ordinary_site_seed.md`](../spec/ordinary_site_seed.md) §3). What remains biased is the
expected frequency.

**And there is a smell that says the same thing as the number.** The prior is a statement about the
*population*. The panel size should not appear in the answer at all, and today it does.

---

## 2. The proposal

**Estimate the two moments directly from the census positions, and never form the curve.**

```text
diversity        π  =  mean over positions of   2 k (2N − k) / (2N (2N − 1))
expected freq.   f  =  mean over positions of   k / 2N
```

for `k` alternative copies among the panel's `2N` chromosomes.

**Both are classical and both are unbiased at any panel size.** The first is Nei's average
heterozygosity, and the `2N(2N − 1)` denominator is the finite-panel correction — the thing that
makes the estimate a property of the population rather than of the panel. The second is unbiased
because `E[k / 2N | f] = f` at every `N`.

**One sample is not a special case, and this is the part worth checking by hand.** At `N = 1` the
denominator is `2N(2N − 1) = 2` and the diversity estimator collapses to `k(2 − k)`, which is 1 if
that genome is heterozygous at the position and 0 otherwise. Its expectation over positions is the
population's heterozygosity exactly. **A single genome estimating population diversity by counting
its own heterozygous sites is the estimator, not an approximation of it.**

**What this removes:** the Beta family and its inability to hold a spike; the projection onto the
panel; the search, which the code's own measurement puts at 11.8 minutes at 3,200 individuals; the
panel size in the answer; and the whole question of which curve to fit.

---

## 3. What is already settled — do not re-derive these

**The census is not ascertained, and the spec says so in as many words.** Positions are kept by
`hash(contig, p, seed) < threshold` over the analysed regions — *"Never select on the data.
'Positions that looked variable' is ascertainment: variability is a function of depth and error, so
selecting on it conditions on the quantity being measured"*
([`../spec/parameter_prepass_census_sites.md`](../spec/parameter_prepass_census_sites.md) §3). So
both estimators are unbiased for the analysed regions with no correction. **What they are not is
genome-wide** — the census domain is the reference intersected with any `--regions` BED — and that
is a statement of scope rather than a bias.

**Linkage costs precision, not correctness.** Census positions are scattered but not independent;
the same section records that treating them as independent makes *"precision figures that are
optimistic, not estimates that are wrong"*. That is exactly right for both moments: linkage inflates
the variance of the estimate and leaves its expectation alone. **So confidence intervals computed
under independence will be too narrow, and §4's question 5 must not use them.**

**The density does not go away.** It is not only an output — it is the joint fit's own prior over
allele frequency, which the expectation step needs to compute per-position posteriors at all. It
keeps being fitted; what stops is *reading the caller's two numbers off it*.

**A richer curve is not the answer, and that was measured rather than argued.** In July the STR
path's frequency prior was swapped for the Ewens `θ/k` spectrum and compared on 561 tomato loci:
indistinguishable at every `QUAL ≥ 10`, 84.1% / 0.19% against 84.0% / 0.18% at `Q ≥ 20`
([`../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md`](../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md)).
**This proposal explains that null from the other side**: the prior can hold only two moments, so
the curve's shape was never reaching a genotype.

---

## 4. The questions, in the order to answer them

### Q1 — Do the estimators recover the truth on drawn cohorts, across the committed range?

The positive control, and without it a clean-looking answer cannot be told from one with no
information in it.

**Sweep** cohort size `1, 2, 3, 5, 10, 25, 63, 200, 1000` crossed with depth `3, 8, 20, 100`
reads a sample, several population shapes including at least one the Beta family cannot hold.
**Report** the bias and the spread of both moments against the values the cohort was drawn with.
**Hold the drawn positions fixed across the cohort-size arms**, as
`examples/ng_joint_sample_count_sweep.rs` does, so nothing moves for a reason unrelated to panel
size.

**Settles it if:** the bias is within the draw-to-draw spread at every cell. **Leaning:** it will
be, on true genotypes — the estimators are unbiased by construction and this arm is checking the
implementation rather than the mathematics.

### Q2 — What does it cost to read from reads instead of genotypes?

**We cannot count `k`.** At three reads a position a heterozygote often looks homozygous and a
sequencing error often looks like a variant; the first deflates the diversity and the second
inflates it, and they do not cancel. So both moments must be computed over **expected** allele
counts under the read model, not counted ones.

**Two arms on the same draws:** the moments from the true genotypes (the oracle), and the moments
from the converged per-position posteriors. **The gap between them is what the read model
contributes**, and it is the number a user needs.

**Settles it if:** the gap is small and shrinks with depth. **Leaning:** it grows as depth falls
and is worst at the corner this caller commits to — one sample at three reads.

### Q3 — Does the EM's own Beta prior pull the answer back toward itself?

**This is the question that can kill the proposal, and it should be run early.** If the moments are
computed from posteriors that were themselves computed under a Beta prior over frequency, they
inherit shrinkage toward that Beta. At high depth the reads dominate and it vanishes; at three
reads it will not.

**The test has to use densities the Beta family cannot represent** — bimodal, or with a mass at
intermediate frequency — because on a density the family *can* hold, shrinkage toward it is
invisible. Compare the moments from posteriors against the moments from true genotypes on those
shapes, at each depth.

**Settles it if:** the pull is negligible at every depth. **If it is not**, the finding is real and
the options are worth naming now: iterate the moments to a fixed point; or widen the EM's frequency
prior deliberately; or accept the shrinkage and report it. **Leaning: it will be visible at three
reads and small at twenty**, and the honest outcome may be that the estimate is unbiased given the
read model rather than unbiased outright.

### Q4 — What do mismapped positions do to the diversity?

A position where two stretches of genome the reference holds once are both piling reads up reads
part-non-reference **in every sample**, which is heterozygosity's own signature. The joint fit
already produces a per-position posterior that a position is mismapped.

**Draw cohorts with a known share of such positions**, and compare the diversity computed with and
without weighting each position by one minus that posterior.

**Settles it if:** the weighted estimate recovers the truth and the unweighted one does not.
**Leaning:** unweighted will be badly high — the fit's own documentation records that ignoring the
class puts observed heterozygosity 50.6% above the truth on a fifty-sample selfing panel at three
reads a position.

### Q5 — Is the answer *usable* at one sample and three reads, not merely unbiased?

Unbiased is not enough: a prior seeded from a number with enormous spread is worse than one seeded
from a constant. **Report the spread of both moments over repeated draws at the hard corner**, and
say what fraction of runs would land far enough out to matter.

**Do not compute the interval under independence** — §3 says why it would be too narrow. Use the
spread across drawn cohorts.

**Settles it if:** the spread is small enough that the seed is better than the species-range
constant. **Leaning:** yes at any real census size, because one genome carries tens of thousands of
census positions; but this is the arm that decides whether a single-sample run should fall back.

### Q6 — Does any of it move a genotype?

The end-to-end question, and the one that says whether the rest was worth doing. **Seed the caller
both ways** — from the current search and from the direct moments — on the same drawn cohorts, and
compare the calls.

**Settles it if:** calls move measurably at the depths where the prior matters. **Leaning: they
will barely move at high depth and may move at three reads**, which is where the prior competes
with the reads. **A null here is a real result and must be reported as one** — it would say the
detour costs nothing in practice and the case for removing it is simplicity and eleven minutes of
run time rather than accuracy.

---

## 5. What would make the answer "no"

Stated in advance so the result is not argued backwards from:

- **Q3 shows a pull toward the EM's Beta that does not shrink with depth.** Then the estimate is
  not unbiased in any useful sense and the proposal needs the fixed-point iteration before it is
  worth anything.
- **Q5 shows the spread at one sample and three reads is wider than the distance between a good
  seed and the species-range constant.** Then the direct estimate is right on average and useless
  in the case this caller most needs it.
- **Q1 fails at large `N`.** That would mean the implementation, not the mathematics, and it is the
  cheapest failure to find — which is why it is first.

---

## 6. Out of scope

- **The repeat-tract prior**, which reads a fitted length spectrum per stratum and takes no such
  detour ([`../spec/population_diversity.md`](../spec/population_diversity.md) §4).
- **Changing the joint fit's model.** The density stays; only what reads it changes.
- **A richer prior family** (§3).
- **Implementation.** This plan produces numbers and a recommendation; the spec and the code follow
  it.

---

## 7. Deliverable

A report under [`../../reports/`](../../reports/) carrying the six answers, each with the
measurement behind it, and a recommendation: adopt, adopt with the fixed-point iteration, or leave
the seam as it is. **If the recommendation is to adopt, a spec for the implementation step follows
it** — the estimators, where they live, what happens where the census is thin, and how the run
reports which numbers it used.

**The standing limit, which applies to every arm above.** This checkout cannot rebuild the tomato
census, so all of it is drawn cohorts. Confirmation on a real panel is
[`../spec/parameter_prepass_cohort.md`](../spec/parameter_prepass_cohort.md) §10's third question
and stays open regardless of what this finds.
