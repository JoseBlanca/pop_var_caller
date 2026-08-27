# ng — the ordinary-site prior's seed: pin its diversity, shrink its shape

**Status:** design spec, 2026-08-26; **§1's diagnosis and §4 are superseded, 2026-08-27**, by
[`ordinary_site_prior_moments.md`](ordinary_site_prior_moments.md) and the measurements behind it
([`../../reports/ng_ordinary_site_prior_moments_2026-08-27.md`](../../reports/ng_ordinary_site_prior_moments_2026-08-27.md)).
**§3 stands and is shipped.**

> **⛔ Two of this document's four claims did not survive measurement, and the reader needs both
> before §1.**
>
> **§1's diagnosis is wrong.** It says the trouble is that a four-number curve with a spike at
> frequency zero cannot be compressed into two numbers without one. **It can, exactly, at every
> panel size**: the prior's two numbers are two *integrals* of the curve, and a spike at zero
> contributes nothing to either. What actually costs the accuracy is that the seed is chosen by a
> **curve fit in class space** — the density is evaluated into the panel's `2N + 1` allele-count
> classes and a pair is searched for that minimises a divergence against them. That objective has
> `2N + 1` terms, so its answer moves with the cohort; a moment match cannot. This document even
> notices the symptom — that the panel size appears in an answer about the population — and does
> not follow it to the cause.
>
> **§4's ramp is retired.** It exists to damp the small-panel noise of a search that is itself
> being deleted, and it was measured to cost 0.62× to 0.92× of the truth at one sample. With the
> two numbers integrated off the curve there is nothing left to blend.
>
> **§3 is untouched and correct.** Solving the total from the measured heterozygosity is half of
> the replacement, and the identity it defines is the one the new spec inverts. It is shipped.
>
> **What to build is [`ordinary_site_prior_moments.md`](ordinary_site_prior_moments.md) §2's step
> one**, which is this document's §3 plus one closed-form line for the expected frequency, minus
> the projection, the search, the ramp and the inbreeding coefficient at this seam.

It supersedes part of [`calling_priors.md`](calling_priors.md) §4.1, which specifies the seed as
*the pair whose predicted allele-count spectrum best matches the fitted one*, and part of
[`population_diversity.md`](population_diversity.md) §3.4, whose ordinary-site ladder switches
between two rungs. §3 below replaces the first; §4 replaced the switch with a ramp and is itself
now replaced.

**Companions:** [`population_diversity.md`](population_diversity.md) §3 (where the two fitted
numbers come from, and how they reach the caller — built as plan step E2f),
[`calling_priors.md`](calling_priors.md) §4 (what the caller does with the pair),
[`str_slippage_level_curve.md`](str_slippage_level_curve.md) §5.1 (the blend-by-precision device
this document borrows).

**It is the ordinary-site path only.** The repeat tract's seed is a length spectrum and a
concentration read off the joint repeat fit ([`population_diversity.md`](population_diversity.md)
§4, built as plan step E2e); nothing here touches it.

---

## 1. What this is, and the two problems it fixes

**The SNP/indel genotype prior starts from two numbers** — `α_ref` and `α_alt`, chromosomes'
worth of belief that a locus carries the reference allele and that it carries the alternative.
Their ratio is the allele frequency the prior expects; their sum is how much conviction that is.
Everything downstream reads exactly these two (`calling_priors.md` §1).

Since plan step E2f both of their inputs are fitted from the run's own data: the population's
allele-frequency density — two point masses and a Beta over what segregates — and the expected
heterozygosity read off it. The seed is then chosen by projecting that density into the panel's
`2N + 1` allele-count classes and searching for the pair whose own predicted classes best match.
**That last step is where both of this document's problems live, and they sit at opposite ends of
the cohort range.**

### 1.1 At the small end the shape is not merely noisy — it is wrong

Measured on drawn cohorts at a known truth (Beta(0.7, 2.5) over 10% of positions, 3,000
positions, 20 reads a sample, six draws at each size). **`a` below 1 is the rare-allele pile-up a
neutral population has**; the truth is 0.7.

| samples | fitted `a`, median (range) | error in the fitted diversity, median / worst |
|---:|---|---|
| 1 | **1.15** (1.11 – 1.26) | +13.6% / +25.8% |
| 2 | **1.28** (1.19 – 1.64) | −4.7% / −20.6% |
| 3 | **1.28** (1.03 – 1.60) | +3.8% / −8.5% |
| 5 | 0.89 (0.74 – 1.12) | +5.5% / +16.8% |
| 10 | 0.91 (0.70 – 1.13) | +3.2% / +7.3% |
| 25 | **0.74** (0.68 – 0.81) | −2.3% / +11.0% |
| 63 | **0.71** (0.62 – 0.92) | −3.6% / +14.3% |

**Two readings, and the second is the one that matters.** The diversity is usable at every size
but one — at a single sample it runs about 14% high. The *shape* is a different story: at one,
two and three samples it comes back **above 1 in every draw**, which is not a noisy version of the
truth but the opposite claim — *there is no excess of rare alleles here*. It is on the right side
from five samples and tight by twenty-five.

**⚠ These are drawn cohorts at 20 reads a sample.** Tomato's panel sits at about three. Nothing
here measures the small end at low depth, and §6 makes that part of the work.

### 1.2 At the large end the pair loses the diversity it was fitted from

**Even with the density known exactly**, compressing it to two numbers costs diversity, and costs
more the larger the panel. Measured at `F = 0` in
[`examples/ng_spectrum_panel_floor.rs`](../../../../examples/ng_spectrum_panel_floor.rs): the
heterozygosity the *fitted pair* implies, against the density's own.

| density | 1 individual | 63 | 200 |
|---|---:|---:|---:|
| tomato-like, Beta(0.20, 1.00), 4 in 1,000 segregating | −0.1% | **−9.9%** | −15.4% |
| human-like, Beta(0.35, 1.20) | −0.1% | **−18.6%** | −25.8% |
| flat over what segregates, Beta(1, 1) | +0.1% | −40.9% | −49.4% |
| middling frequencies, Beta(4, 4) | +0.0% | −53.9% | −61.0% |

**The projection is not what loses it.** The class weights carry the density's own heterozygosity
at every panel size — pinned by
`joint::fit::tests::the_classes_carry_the_densitys_heterozygosity_at_every_panel`. What loses it
is the search that follows: a two-parameter Dirichlet cannot hold mass piled at *invariant* **and**
a spread over what segregates, so it compromises, and the more classes it is fitted over the more
it trades the ends against the middle.

### 1.3 Why this is not answered by a better-fitting family

**Because that was measured, on this project's own data, and came back null.** In July the STR
path's frequency prior was swapped from the Dirichlet-multinomial to the **Ewens `θ/k` site
frequency spectrum** — the neutral spectrum itself — and compared on 561 silver-standard tomato
loci ([`../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md`](../../reports/ssr_marg_sfs_genotype_prior_2026-07-09.md)).
The two were indistinguishable at every `QUAL ≥ 10`: recall within 0.7 points, false positives
within 0.03%; at `Q ≥ 20`, 84.1% / 0.19% against 84.0% / 0.18%. Both beat the plug-in prior by 3
to 4 points, so the comparison could see differences — it could not see that one. The verdict
recorded then: **the SFS's edge lives in marginalising the frequency in the *emit* decision, not
in a better genotype frequency model.**

**And the cheap version is cheap for a reason that generalises**: the finite-`k` Ewens spectrum
**is** a symmetric Dirichlet(`θ/k`), so "use the real spectrum" was the same family at a different
base measure rather than a different family.

**What that probe did not cover**, stated so nobody reads it as more than it is: the STR path
only, tomato at about three reads a plant, 561 loci, with depth recorded as the wall — and two
shapes *inside* the Dirichlet family.

**So a richer family is deliberately out of scope here**, and §7 records what would reopen it.

---

## 2. Goals and non-goals

### Goals

1. **The seed's implied diversity is the diversity that was measured**, at every panel size and
   every shape — by construction rather than by luck.
2. **The shape moves smoothly from the neutral one to the panel's own**, as the panel earns it,
   with no size at which a run's answer jumps.
3. **How much of the shape came from the panel reaches the run's output.** Two runs that leaned
   differently on their data must be distinguishable in what they emit, which is the complaint
   `calling_priors.md` §4 makes about production's own fallback.
4. **An answer at one sample and at a thousand**, and at three reads a position as well as at
   three hundred (`doc/devel/specs/design_principles.md` §0).

### Non-goals

- **A richer prior family** (§1.3). Not a mixture of Betas, not a point mass at *invariant*, not
  the Ewens spectrum.
- **Changing the projection.** `FrequencyDensity::allele_count_classes` is exact and stays.
- **Changing the search.** `fill_expected_spectrum` and the two-parameter fit stay as they are;
  what changes is what is done with the pair they return.
- **The repeat-tract seed**, which is a different object on a different path.
- **The site quality**, where §1.3's probe says the spectrum's real edge lives. §7 keeps it.

---

## 3. Pin the diversity

**Decision: the pair is set from two numbers — the expected allele frequency and the measured
diversity — and the search supplies only the first.**

A Dirichlet(`α_ref`, `α_alt`) with total `A = α_ref + α_alt` and expected frequency
`f = α_alt / A` implies, for a diploid drawn from it,

```text
heterozygosity  =  2 f (1 − f) · A / (A + 1)
```

which is the Beta-binomial at one alternative copy in two draws. **So a pair is exactly `(f, A)`
in other clothes**, and fixing the heterozygosity to the measured `θ` determines `A` once `f` is
known:

```text
t = θ / (2 f (1 − f)),        A = t / (1 − t),        α_alt = A f,   α_ref = A (1 − f)
```

**Take `f` from the fit and `A` from the measurement, not the other way round.** The two
quantities are not equally trustworthy: `θ` is read straight off the fitted density with no panel
in it, and §1.1 shows it within a few per cent from two samples up, while `A` is the number the
point masses corrupt (§1.2). **The alternative — keep `A` from the fit and solve for `f`** — puts
the measured number into the quantity the fit gets right and the fitted number into the one it
gets wrong.

**This is a rescale, not a replacement**, and the size is worth knowing before anyone fears it: on
the tomato-like density at 63 individuals the fitted pair is `(0.1618, 3.180e-4)` and the pinned
one is about `(0.1829, 3.59e-4)` — the total rises by 13% and the expected frequency does not
move.

**The neutral rung is already pinned and needs no change.** At `(1, θ)` the implied heterozygosity
is `2θ / ((1 + θ)(2 + θ))`, which is `θ` to within 0.15% at `θ = 10⁻³`. **That is what makes §4
possible**: after this section both ends of the ramp imply the same diversity, so the ramp
interpolates *shape alone*.

### 3.1 Three ways this can have no answer, and none of them may be silent

**`2 f (1 − f) ≤ θ`.** The shape's own maximum implied diversity is below the measurement, so no
total reaches it. **Rule: fall to the neutral rung and mark the run.** Do not rescale toward the
ceiling, do not clamp. *This is the same failure the repeat-tract seed used to have — a shape
scaled to a measurement it could not reach — which fired at every tract at one outbred sample and
was deleted with the construction that caused it (`population_diversity.md` §4.2). It must not
return silently on the other path.*

**`θ > 1/2`.** `2 f (1 − f)` cannot exceed a half, so no pair implies it. `ExpectedHeterozygosity`
admits the whole of `[0, 1]`, so this is expressible. **Rule: refuse the run's assembly with a
message naming the fit**, because a heterozygosity above a half is not a thin estimate, it is a
fit that did not converge.

**`θ = 0`.** A cohort with no variation at all. `A` goes to zero and every entry with it.
**Rule: the seed is floored at `MIN_ALT_CONCENTRATION`**, the same floor every other seed builder
in this tree applies, and the run says the diversity was zero.

---

## 4. Shrink the shape

> **⛔ Superseded 2026-08-27** by
> [`ordinary_site_prior_moments.md`](ordinary_site_prior_moments.md) §6.1. The ramp damps the
> small-panel noise of the class-space search, and that search is being deleted — measured, the
> blend costs **0.62× to 0.92× of the truth at one sample**, and its own sweep put the best
> half-weight panel size at zero on every arm. **Kept as the record of what was tried.** The one
> paragraph below that survives is the last: the two rungs do not disagree only about shape, and
> the size of the disagreement in the *total* is what a reader should carry away.

**Decision: the expected frequency is interpolated between the neutral one and the panel's own,
by a weight that rises with how much the panel supports a shape. The diversity stays pinned
throughout (§3), so nothing but the shape moves.**

```text
f_neutral = θ / (1 + θ)                      the neutral rung's own expected frequency
f_fitted                                     from the search
ln f = (1 − w) · ln f_neutral + w · ln f_fitted
```

**In log space, because the two can be orders of magnitude apart** — on the tomato-like density at
63 individuals they are `6.1e-4` and `2.0e-3`, and a linear blend of numbers that small is
dominated by the larger one at every weight but zero.

Then `A` follows from §3's identity and the seed is `(A(1 − f), A f)`.

**`w = 0` is exactly the neutral rung and `w = 1` is exactly §3's pinned fit.** The two rungs
`population_diversity.md` §3.4 switches between are the two ends of this, which is what makes this
a ramp rather than a third rung: **the ladder's top two rungs stop being alternatives.** Its
bottom rung — no diversity fitted at all, so the species-range constant — is untouched and stays a
rung, because there is nothing to interpolate when there is no measurement.

### 4.1 Where the weight comes from

**It is set from measurement, not chosen**, and this is the part of the work that is a measurement
rather than an edit.

**Form:** a weight in `[0, 1]`, rising with the panel, and **one constant to fit** — the
half-weight point:

```text
w = N / (N + N₀)
```

`N` is the panel size in diploid individuals. **`N₀` is where the two ends are equally
trustworthy**, which is the criterion rather than a taste: below it the neutral shape is the
better guess and above it the panel's own is.

**What to measure, and the shape of the run is already established.**
[`examples/ng_joint_sample_count_sweep.rs`](../../../../examples/ng_joint_sample_count_sweep.rs)
refits *the same drawn positions* at 2, 3, 5, 10, 25 and 50 samples — holding the positions fixed
across the arms is what stops the answer moving for a reason that has nothing to do with panel
size, and the new sweep must do the same. For each cohort size, and **at several depths**:

1. draw a cohort from a known density;
2. fit it, project it, and take `f_fitted`;
3. compare `f_fitted` and `f_neutral` against the truth's own expected frequency;
4. `N₀` is where the two errors cross.

**Both axes, not one.** §1.1's table is at 20 reads a sample and tomato sits near three; a weight
fitted at high depth and applied at low would lean on a shape the reads never supported.
**If `N₀` turns out to depend on depth, the weight takes both** — and that is a finding, not a
complication.

**⚠ It has to be drawn data, and the spec says so rather than leaving it to be discovered.** This
checkout cannot rebuild tomato's census — the CRAMs are not in the repository — so the sweep runs
on cohorts drawn at known parameters across the committed range. That is what
`design_principles.md` §0 asks for in any case; what it is not is a confirmation on real data, and
§7 keeps that open.

### 4.2 What travels to the output

**The weight itself.** A run at `w = 0.1` and a run at `w = 0.9` used different information and
must not look the same in what they emit (goal 3). The existing `SeedRegime` has three variants
and no room for a continuum; **the fitted variant carries the weight**, or an implementer may
replace the top two with one variant that does — either is fine, and what is not fine is a run
whose output cannot say how much of its shape it borrowed.

**The two things already carried stay carried**: how far the fitted pair sits from the measurement
it was fitted to (`SpectrumMatch`), and whether the search ran out of range.

**And the refusals of §3.1 are regimes, not silences.** A run that fell to the neutral rung
because no total could reach its measured diversity must say so, distinguishably from one that
fell there because `w` was zero.

---

## 5. What this replaces, section by section

| document | what it says now | what this changes |
|---|---|---|
| `calling_priors.md` §4.1 | the seed is the pair whose predicted spectrum best matches the fitted one | the search supplies the *shape* only; the total comes from the diversity (§3) |
| `calling_priors.md` §4.1 | *"A spectrum makes the diversity moot: it carries its own scale"* | false after §3 — the diversity is read on **every** rung |
| `population_diversity.md` §3.4 | a ladder whose top two rungs are alternatives | the top two rungs are the ends of one ramp (§4); the bottom rung is unchanged |
| `population_diversity.md` §9 q3 | where to put a panel-size floor, *"confirm before code"* | **answered**: there is no floor, because there is no switch left to place one at. The measurement that closed it is §4.1's, not the divergence sweep §9 proposed — that statistic is smallest at the smallest panel and cannot locate a floor |
| `seed_generic.rs` `project_spectrum_seed` | *"a spectrum arrived — fit it … `diversity` is not read"* | the diversity is read on every path through this function |

---

## 6. How we know it works

1. **The seed's implied diversity is the measured one**, at every panel size from 1 to the
   projection's ceiling and at every shape in §1.2's grid — asserted, not sampled.
2. **The two ends are exactly the two rungs.** At `w = 0` the seed is `(1, θ)` to the same 0.15%
   §3 quotes; at `w = 1` it is §3's pinned fit. Both pinned by tests that fail if either end
   drifts.
3. **The middle is better than either end.** On drawn cohorts across the committed range, the
   blended seed's expected frequency is closer to the truth than the neutral shape below `N₀` and
   closer than the fitted shape above it — which is the claim `N₀` is fitted to, checked on draws
   the fit did not see.
4. **The weight is monotone in the panel and never leaves `[0, 1]`.**
5. **Each of §3.1's three failures is reachable by a test and reported distinguishably** — the
   unreachable diversity, the impossible one, and the zero.
6. **Cost is once per run.** Nothing here is per locus, per sample or per pass; the seed is
   projected once and frozen (`calling_priors.md` §2.3).
7. **Two runs that leaned differently on their panels emit different records** (goal 3).

**⚠ Built on 2026-08-26 (`ng-seed-shrinkage`), and three of these seven need correcting against
what was measured.** The design is unchanged; these are factual amendments to the checklist.

- **Item 2's 0.15% is the diversity's, not the pair's.** At `w = 0` the seed's implied
  heterozygosity is `θ` by construction, and the *pair* is `(1, θ)` moved up by about `3 θ` on each
  concentration — because the literal `(1, θ)` implies `2θ/((1+θ)(2+θ))`, short of `θ` by about
  `1.5 θ`, and making that up costs twice as much on the pair. Measured: **0.03% at `θ = 10⁻⁴`,
  0.18% at 6 × 10⁻⁴ and 3.07% at 10⁻²**
  (`seed_generic::projection_tests::a_neutral_panel_projects_to_one_and_theta`). So the two ends
  are the two rungs at every realistic diversity and visibly not at `θ = 10⁻²`.
- **Item 3 is not what the measurement supports, and it is not what `N₀` was fitted to.** Per drawn
  cohort the blended shape can never be worse than *both* ends — in log space it is a convex
  combination of the two errors — so "better than either end" is not a testable claim in that
  direction. What was measured over 42 held-out cells is: at least as good as both in some, between
  the two in the rest, worse than both in none. And the crossing this item assumes — the panel's
  own shape worse below `N₀` and better above — **does not exist**: the panel's own shape is
  exact at one individual and degrades monotonically with the panel. §4.1's criterion was replaced
  by *the value that puts the blended shape nearest the truth over the whole grid*. See
  [`../../reports/implementations/ng_seed_shrinkage_2026-08-26.md`](../../reports/implementations/ng_seed_shrinkage_2026-08-26.md) §5.
- **Item 6's second half is arithmetic, not a test.** The 399-prediction count is asserted; that the
  pin and the blend add no prediction is visible in the code and pinned by nothing, because there is
  nothing to instrument — both are closed-form arithmetic on the pair the search already returned.

---

## 7. Open questions, and one thing deliberately left shut

**Open.**

1. **Does `N₀` depend on depth?** §4.1 measures both axes. **Leaning: it does** — the shape is
   information about how many *chromosomes* carry the allele, and at three reads a sample a
   genotype is barely observed, so a panel earns its shape more slowly. **What would settle it:**
   the sweep of §4.1, read as a surface rather than a curve.
2. **Confirmation on a real cohort.** §4.1's sweep is drawn data by necessity. **What would settle
   it:** subsample the tomato panel and refit, which is
   [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §10's third question and needs the
   census this checkout cannot rebuild.
3. **Should the *inbreeding* used in the projection enter the weight?** The fitted pair moves 8.6%
   to 14.0% across `F = 0.6` to `0.9` (`seed_generic.rs`), and every measurement in this document
   is at `F = 0`. **Leaning: no** — inbreeding changes what the panel's classes look like, not how
   much the panel knows. **Confirm before relying on it.**

**Shut, and here is the key.** **A richer prior family is not the answer to §1.2**, on the
evidence of §1.3's null result and because the finite-`k` Ewens spectrum is the Dirichlet this
already uses. **What would reopen it** is the one place that probe said the spectrum's edge lives:
if a prior that kept the density's *point mass* — a locus is invariant with probability `p`, and
otherwise its frequency is Beta — moved the **site quality**, rather than the genotypes. That is a
bounded probe against [`calling_quality.md`](calling_quality.md) once the caller emits, and it
costs the calling loop a `logsumexp` per genotype and a mixing-weight update per sample per pass,
which is why it is not free and not now.
