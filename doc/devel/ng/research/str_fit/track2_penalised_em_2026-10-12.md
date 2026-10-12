# Track 2: penalised EM — every repeat class of a motif length fitted at once

**Status:** research report, 2026-10-12, branch `str-fit-d`. Avenue D of
[`../str_fit_research_plan.md`](../str_fit_research_plan.md). Measured on drawn data two ways: a
motif length composed with the per-class drawer `draw_stratum` along known curves with a deliberate
kink (§3), and phase 0's shared bench at its corners (§4). Nothing in the caller changes. Every time
quoted was measured while other fits compiled and ran on the same machine, so **every time is
provisional**. **§10 (2026-10-13) adds where the time goes and what was cut, standard errors and
their honesty on the bench, and what lies behind the overshoot at 1 sample × 3 reads.**

---

## 1. The question, and the short answer

**The question (plan §4, avenue D):** does borrowing strength across repeat count *during* the fit
beat fitting each class on its own and smoothing afterwards — at thin classes and at one sample —
without smoothing away a real kink such as HG002's 9-repeat homopolymers?

**The short answer:**

- **Yes at thin classes and at one sample, for every number; no loss elsewhere except at the top of
  the repeat range on the thinnest evidence.** On a drawn motif length (§3) at 1 sample × 3 reads,
  penalised counting roughly halves every error: the slip rate at thin classes (5 to 30 loci) 0.97 →
  0.63 relative error, the allele length variation 0.31 → 0.13 of chromosomes misplaced. From 10
  samples × 10 reads up the gain moves to the numbers a class measures worst — the shorter share,
  fall-off, between-locus share F and allele length variation, whose thin-class errors fall
  by 13 to 70%, most often by a third to a half — while the slip rate, which every read measures, gains little and the existing curve-and-blend
  matches it at middling classes. At 2,000 samples × 3 reads the two fits agree.
- **A well-measured kink survives.** The slip rate at 9 repeats, drawn at 0.53 of its curve, is kept
  as well as the counting fit keeps it from 447 slipped reads up (0.053 against 0.053 relative error
  at 10 × 10; 0.018 against 0.018 at 50 × 30), and it no longer drags the rest of the curve's
  smoothing down (§3.3). At 671 slipped reads in an inbred cohort it is pulled further than by the
  counting fit (0.047 against 0.032). Getting there took a heavier-tailed bend penalty than the
  blend's knee implies (§2.3).
- **On the shared bench's tomato-like truth** (§4), penalised counting predicts held-out loci almost
  as well as the truth itself at one sample: −0.04 to −0.05 log-likelihood units a locus, against the
  counting fit's −0.57 to −0.91 (at 3 reads, after the caller's curves) and −46 to −91 (at 300 reads).
  Much of that gap is the counting fit's zero-share defect and its F bound, not the curves; the runs
  that split them are listed in §8. At 2,000 samples × 3 reads (one-base runs, one seed) the gap is −1.35 against
  −428, again mostly the zero-share defect.
- **On the drawn grid a round costs what a round of the counting fit costs**, but it runs every
  class until the slowest settles:
  0.9 to 1.5 times the counting fit's time on the drawn grid at ±4, and 1.5 to 35 times on the
  bench at ±9, where its strength search and pairs alternations weigh more (§5; every time is
  provisional, the machine was shared).
- **F is reported as absent at one sample**, and **no allele length can be driven to a share of
  zero** (§2.5, §2.6).
- **Not settled:** whether choosing the smoothing from the E-step's expected counts smooths too
  little at one sample; a strength choice that sees all the numbers at once, without which weighing
  allele draws as clustered by locus misleads the stutter at 2,000 samples (§2.5, §3.4); standard
  errors.

---

## 2. What the fit does

### 2.1 Vocabulary

- **Class.** Every locus with one motif length and one reference repeat count. The fit here takes
  every class of one motif length at once.
- **The counting fit.** Today's prototype ([`composite_fit.rs`](../../../../../src/parameter_estimation/joint/ssr_fit/composite_fit.rs)),
  which fits one class at a time by expectation maximisation (EM): each round shares every read
  between "read exactly" and "slipped this far, this way", and every sample's alleles among the
  lengths, given the class's current numbers (the **E-step**); the class's new numbers then follow
  from those shared-out counts in closed form (the **M-step**).
- **Logit.** `ln(p / (1 − p))`, the scale on which a probability between 0 and 1 becomes a number
  that can be any size. A straight line on the logit scale rises like an exponential while the
  probability is small and bends to approach 1, so it never leaves the range a probability can take.
- **Bend.** The second difference of a curve at a repeat count: its value there, doubled and taken
  from the values one repeat either side. A straight line has no bends.
- **Penalty strength.** How heavily the fit charges a curve for its bends. Zero lets every class
  keep its own answer; a very large strength forces a straight line on the logit scale.
- **Marginal likelihood.** How probable the counts are under a given penalty strength, averaged over
  every curve that strength allows. It is highest at the smoothness the counts themselves support,
  which is what lets the fit choose the strength rather than be told it.

### 2.2 The model

**The E-step is the counting fit's, untouched.** Each round, every class's expected counts are
gathered at that class's current numbers, by the same code the counting fit runs
(`composite_fit::expected_counts`, split out of the counting fit with no number changing: §6).

**The M-step fits every class of the motif length together.** Each of the class's numbers is the
answer to a counting question, and the M-step answers each question for all the classes at once,
as a curve across repeat count that is charged for bending:

| number | the counts it is fitted to, per class | scale the curve is drawn on |
|---|---|---|
| slip rate | slipped reads, out of all reads | logit |
| shorter share | slipped reads that went shorter, out of slipped reads | logit |
| fall-off | "one repeat further" steps, out of all steps of the slipped reads (a slip of `k` repeats is `k − 1` of them and one stop) | logit |
| between-locus share F | pairs of a sample's two copies that drew one lineage from the locus, out of pairs not identical by descent | logit |
| allele length variation | allele draws at each length | each length's log-ratio to the reference length's share |

With no penalty each row reproduces the counting fit's closed-form answer exactly — checked by a
test (§6). With a penalty, a class whose counts say little is pulled toward the curve its
neighbours draw, by as much as its counts leave undecided; a class with many counts keeps its own
answer.

**One value per repeat count, not a spline with knots.** Each curve has one free value at every
repeat count from the lowest class with reads to the highest, and the penalty is on the bends
between consecutive values. This is the limiting case of the P-splines of Eilers and Marx (a
penalised spline whose knots sit at every repeat count), known as the Whittaker smoother. Repeat
counts are whole numbers and a motif length has a few dozen at most, so nothing is gained by fewer
knots, and nothing then has to choose where they go. A repeat count with no loci between two that
have them takes the curve's value there (the bend penalty draws a straight line through it on the
logit scale). Outside the fitted range a curve is held at its end value, never continued, as the
slippage curve already is ([`str_slippage_level_curve.md`](../../spec/str_slippage_level_curve.md) §6).

**The allele length variation: each length's own curve, one strength for all lengths.** For every
allele length other than the reference length, the logarithm of its share over the reference
length's share is a curve across repeat count; the shares at a class are then those ratios
normalised. Why this choice:

- **The reference length is present at every class**, so its share is never near zero and every
  ratio to it is well defined.
- **Nothing is assumed about how shares fall away from the reference length.** Each length is
  smoothed only along repeat count, never against its neighbouring lengths, so the two edge classes
  — which mean "this far or further" and so sit above the falling tail — are smoothed like any other
  length and are never forced onto a tail. A shape across lengths is track 3's question.
- **One penalty strength for all lengths** pools the evidence about how fast the allele spread
  changes with repeat count; eighteen separate strengths would each rest on a thin length's counts.
  The bends of each length are eased separately (below), so one length's kink does not ease the
  others.

The rejected alternative, a penalty on the difference between neighbouring classes' shares
themselves, is the same idea on the share scale; it would let a share become negative or fail to
sum to one between steps, which the log-ratio scale cannot.

### 2.3 How smooth each curve is, and why a kink survives

**The penalty strength is chosen from the data by the marginal likelihood**, in the
Laplace approximation that restricted maximum likelihood rests on (Wood 2011, *J. R. Stat. Soc. B*
73:3). For each curve and each candidate strength, the fit finds the best curve and scores the
strength by how well the curve fits, minus what its bends cost, minus a term that charges a weak
strength for the freedom it gives. The strength is searched over 10⁻⁴ to 10¹⁰ on a grid, from the
strongest to the weakest, then refined. Each of the five numbers gets its own strength.

Why this and not the two alternatives the brief named:

- **Cross-validation over whole loci** would refit the whole EM, every round's E-step included, once
  per held-out fold per candidate strength: about 5 folds × 30 strengths, so roughly 150 times the
  fit's cost. The marginal likelihood costs one small banded solve per candidate: choosing every
  strength, and every other part of the M-step, took 1 to 9% of the fit's time on the drawn grid
  (§5).
- **Generalised cross-validation** costs the same as the marginal likelihood but is known to choose
  too little smoothing now and then, by a wide margin (Wood 2011 §1); the marginal likelihood's
  optimum is flatter-bottomed and more stable.

**One reservation about the choice, unmeasured.** The marginal likelihood here is computed from the
*expected* counts of the E-step, as though they were observed. They carry more apparent information
than the reads do when the genotypes are uncertain, which is the case at one sample and three reads,
so the strength chosen there may be weaker than the reads would support. The shared harness's
held-out loci can measure whether that matters (§8).

**A strong, well-measured kink survives because the penalty on a bend stops growing as its square
beyond a knee.** A bend is measured in standard deviations of the curve's smoothness. Up to 2.5 of
them — the knee of the slippage curve's blend
([`str_slippage_level_curve.md`](../../spec/str_slippage_level_curve.md) §7.2) — it costs its square,
as in any penalised spline; beyond, it costs only `2 · ln(bend / 2.5)` more, so its weight in the
fit falls as `2 / bend²`. As a prior on the bends this is a normal distribution whose tails fall
only as `bend⁻²`: a bend of many deviations is rare but allowed, and a class whose counts insist on
one pays a bounded price for it. A thin class's noise cannot earn the same allowance, because the
penalty pulls a thin class's bend well inside the knee before its size is ever judged. The strength
is chosen under this same penalty, its tails included, and searched from the strongest down, so the
marginal likelihood compares "smooth with one kink" against "rough everywhere" on the counts.

**Why the tail costs 2 and not 6.25.** The first version charged `2.5² · ln(bend / 2.5)`, which
keeps the penalty's slope continuous at the knee. With that tail, at every drawn setting with a kink
except 1 sample × 3 reads, the marginal likelihood chose a slip-rate curve smoothed 60 to 100 times
less (strength 5 to 10 against 300 to 700 without the kink) rather than pay for one eased bend, and
the thin classes lost what they had borrowed: at 1 sample × 300 reads their slip rate's error went
from 0.023 without the kink to 0.049 with it. At a tail of 2 the bend at 9 repeats is eased in every
draw at 10 × 10, 50 × 30 and 1 × 300, the strength stays at 25 to 610 (all but one draw above 100),
and the thin classes are back to 0.024 (§3.3). Where no bend passes the knee the two tails give the
same answer to the third decimal (1 sample × 3 reads, and the two settings drawn without the kink at
10 × 10 and 1 × 300).

**Strengths are chosen afresh each round until the fit has nearly settled, then held**: until no
number moves by more than 10⁻⁵, or for at most 200 rounds. Both limits were needed. Choosing every
round left one draw in five at 50 samples × 30 reads alternating between two strengths and never
settling in 500 rounds. Holding after a fixed 30 rounds instead froze, at 2,000 samples (about 300
rounds), strengths chosen on counts still shaped by the flat starting shares. With the strengths
held, the remaining rounds are penalised EM at a fixed penalty, which climbs its objective every
round (Green 1990, *J. R. Stat. Soc. B* 52:443).

### 2.4 The between-locus share in an inbred cohort

At a homozygote excess near one, a sample's two copies are nearly always one draw and say almost
nothing about how much loci differ, so the counting fit reads F off pairs of samples at a locus
instead (`fit_class_by_counting_with_pairs`). The penalised fit does the same when asked
(`pair_alternations`): after penalised counting settles, each class's pairwise likelihood is
tabulated over F at its current stutter and allele length variation, the F curve is fitted to those
tables with the same penalty and strength choice, and penalised counting runs again with F held.
At most `pair_alternations` times (3 in the drawn grid, 10 on the bench, as method B uses), stopping
once no class's concentration moves by more than 1%.

### 2.5 Allele draws: never zero, and why they are not reweighted by locus

**One allele draw a class is spread evenly over the lengths** (1/19 of a draw each at ±9) before
the shares are chosen. Without a floor, a length the reads stop supporting gets a share near zero,
the next E-step gives it no draws, and it can never come back: EM has an absorbing boundary at zero.
The counting fit has this defect, and on the shared bench it loses 46 to 428 log-likelihood units a
held-out locus at 300 reads and at 2,000 samples (§4). A log-ratio curve alone does not prevent it, because the ridge
that keeps an unsupported log-ratio finite still lets it fall to a share of about 10⁻⁹. **Half a draw
a length (the Jeffreys prior) was tried first and was too much:** at one sample and 100 loci a
class, 19 half-draws moved enough mass into the tails to bias the reference length's share by −0.13
(one-base runs) and −0.20 (two-base runs) on the bench. At one draw a class the bias is −0.01 to
−0.03 (§4).

**An option, off by default: weigh each class's allele draws as clustered by locus**
(`weigh_allele_draws_by_locus`), dividing them by `1 + (m − 1) · F`, `m` being the allele draws a
locus holds (two per sample with reads) and `F` the class's between-locus share. The E-step counts
every sample's alleles as independent draws from the class, but samples at one locus share its
frequencies; under the model's own Dirichlet a locus's allele counts vary more than independent
draws by exactly that factor. Dividing every count of a class by one number leaves that class's own
shares unchanged. **It helps the allele length variation and, at many samples, wrecks the stutter**
(§3.4): at 50 × 30 it halves the thin classes' allele-length error; at 2,000 × 3 it leaves the
allele-length curve straight, slightly wrong at rich classes, and the stutter absorbs the misfit
(slip rate 12% wrong at rich classes against 1.4% without it). The strength is chosen from the
allele counts alone, which cannot see what the reads say about the spread through the stutter. So it
is off until the strength choice can see both.

### 2.6 F at one sample: absent

**When no locus of the motif length has two samples with reads, the between-locus share is
reported as absent** (`None`), with its curve and concentration. At one sample how much loci differ
cannot be told from that sample's own inbreeding, and the project reports such a number as absent
rather than as a fitted value. The genotype prior still needs a value, so the fit works at the one
read off the sample's own two copies and keeps it as `working_concentration`, which is not an
estimate. The bench's `D_working_F` reports that working value only so that the held-out score,
which needs every number, can be computed at one sample (§4).

### 2.7 The interface

```rust
pub fn fit_period_by_penalised_counting(
    classes: &[StratumEvidence],      // every class of one motif length, any order
    homozygote_excess: &[f64],        // by sample
    config: &SsrFitConfig,            // the allele window
    fit: &PenalisedFitConfig,         // start, rounds, tolerance, smoothing, knee, tail, floor, pairs
) -> PenalisedPeriodFit               // per class (None if it had no reads), the curves, cost
```

`PenalisedPeriodFit::curves` holds one `ParameterCurve` per number (and one per allele length),
each with its values at every repeat count of its range, its strength, and the repeat counts where a
bend was eased; `ParameterCurve::at(repeats)` holds flat outside the range. A class with no reads is
not fitted and gets `None`; the curves still give it a value. Each library is scored at its class's
stutter: library multipliers are not fitted here, as in the counting fit.

---


## 3. One drawn motif length, along known curves

### 3.1 What was drawn

**Every class of one motif length, 8 to 20 repeats, drawn along known curves with one deliberate
kink**, by `examples/str_fit_penalised_em.rs` through the per-class drawer `draw_stratum`. Every
number bends across repeat count, so none is a straight line on the logit scale (the shape the
penalty leaves alone), and the fall-off dips and rises again:

| repeats | 8 | 12 | 16 | 20 |
|---|---:|---:|---:|---:|
| slip rate | 0.0050 | 0.0243 | 0.0781 | 0.1641 |
| shorter share | 0.599 | 0.761 | 0.808 | 0.776 |
| fall-off | 0.231 | 0.184 | 0.199 | 0.287 |
| between-locus share F | 0.299 | 0.200 | 0.143 | 0.112 |
| ratio of each allele length's share to the next one in | 0.15 | 0.29 | 0.43 | 0.57 |

**The kink:** the slip rate at 9 repeats is drawn at 0.53 of its curve (0.0041 where the curve says
0.0077), as HG002's 9-repeat homopolymers sit 87% below the current curve. **Loci per class** fall as
real catalogs do: 2,000 at 8 repeats, 1,100 at 9, then 605, 333, 183, 101, 55, 30, 17, **none at
17**, 12, 8 and 5 at 20. Shorter alleles are 1.3 times as common as longer ones at the same
distance, and the two edge lengths carry everything at or past them, so they sit above the falling
tail. **Allele window ±4** (reads and alleles), not the shipped ±9, so that many draws fit in the
time; every method is scored at the same window.

**Settings**, every sample with the same homozygote excess:

| setting | excess | draws | slipped reads at the kink | between-locus share read from |
|---|---:|---:|---:|---|
| 1 sample × 3 reads | 0.3 | 10 | 13 | each sample's two copies |
| 10 samples × 10 reads | 0.3 | 5 | 447 | each sample's two copies |
| 50 samples × 30 reads | 0.3 | 5 | 6,712 | each sample's two copies |
| 1 sample × 300 reads | 0.3 | 5 | 1,342 | each sample's two copies |
| 50 samples × 3 reads, inbred | 0.95 | 3 | 671 | pairs of samples, 3 alternations |
| 2,000 samples × 3 reads | 0.3 | 1 | 26,849 | each sample's two copies |

The 1 × 3, 10 × 10 and 1 × 300 settings were also drawn without the kink.

**Three fits on the same draws:** the counting fit class by class (`counting`); the same, then the
slip rate drawn across repeat count by the existing curve and blended with each class's own answer
(`blend`, [`str_slippage_level_curve.md`](../../spec/str_slippage_level_curve.md) §4 and §7; the
other numbers have no such curve after the counting fit); and penalised counting (`penalised`), at
its defaults.

**What is scored**, per class, over the draws: the root-mean-square error of the slip rate relative
to its truth, of the shorter share, fall-off and between-locus share in their own units, and for the
allele length variation the mean **total-variation distance**, the share of chromosomes that would
have to change length to turn the fitted spread into the true one. Classes are grouped by how many
loci they hold: rich (8, 10 and 11 repeats: 333 to 2,000 loci), middling (12 to 14: 55 to 183), thin
(15, 16, 18 to 20: 5 to 30), the empty class at 17 (penalised only: the curve's value there), and the
kink at 9 (1,100 loci).

### 3.2 Results

Each cell is counting → **penalised**; for the slip rate, counting / blend → **penalised**. Lower
is better.

**Slip rate (relative error; counting / blend → penalised)**

| setting | rich | middling | thin | no loci (17) | the kink (9) |
|---|---|---|---|---|---|
| 1 × 3 | 0.389 / 0.340 → **0.210** | 0.743 / 0.600 → **0.389** | 0.967 / 0.983 → **0.625** | — / — → **0.522** | 0.861 / 0.764 → **0.809** |
| 10 × 10 | 0.056 / 0.057 → **0.053** | 0.080 / 0.054 → **0.055** | 0.098 / 0.075 → **0.071** | — / — → **0.075** | 0.053 / 0.064 → **0.053** |
| 50 × 30 | 0.008 / 0.007 → **0.008** | 0.017 / 0.014 → **0.016** | 0.026 / 0.024 → **0.017** | — / — → **0.008** | 0.018 / 0.018 → **0.018** |
| 1 × 300 | 0.021 / 0.019 → **0.020** | 0.036 / 0.020 → **0.027** | 0.054 / 0.037 → **0.024** | — / — → **0.017** | 0.032 / 0.040 → **0.035** |
| 50 × 3, inbred | 0.028 / 0.036 → **0.029** | 0.060 / 0.037 → **0.042** | 0.077 / 0.048 → **0.041** | — / — → **0.025** | 0.032 / 0.056 → **0.047** |
| 2,000 × 3 | 0.017 / 0.017 → **0.014** | 0.015 / 0.015 → **0.017** | 0.048 / 0.047 → **0.045** | — / — → **0.050** | 0.031 / 0.031 → **0.032** |

**Shorter share**

| setting | rich | middling | thin | no loci (17) | the kink (9) |
|---|---|---|---|---|---|
| 1 × 3 | 0.228 → **0.187** | 0.272 → **0.129** | 0.406 → **0.360** | — → **0.334** | 0.318 → **0.195** |
| 10 × 10 | 0.027 → **0.014** | 0.024 → **0.013** | 0.057 → **0.034** | — → **0.019** | 0.033 → **0.013** |
| 50 × 30 | 0.005 → **0.003** | 0.007 → **0.004** | 0.010 → **0.006** | — → **0.006** | 0.002 → **0.003** |
| 1 × 300 | 0.013 → **0.008** | 0.014 → **0.007** | 0.025 → **0.012** | — → **0.011** | 0.012 → **0.005** |
| 50 × 3, inbred | 0.022 → **0.012** | 0.029 → **0.007** | 0.029 → **0.021** | — → **0.011** | 0.016 → **0.009** |
| 2,000 × 3 | 0.008 → **0.007** | 0.014 → **0.002** | 0.009 → **0.005** | — → **0.007** | 0.002 → **0.003** |

**Fall-off**

| setting | rich | middling | thin | no loci (17) | the kink (9) |
|---|---|---|---|---|---|
| 1 × 3 | 0.174 → **0.103** | 0.161 → **0.108** | 0.311 → **0.171** | — → **0.146** | 0.390 → **0.100** |
| 10 × 10 | 0.012 → **0.011** | 0.029 → **0.012** | 0.037 → **0.023** | — → **0.011** | 0.014 → **0.006** |
| 50 × 30 | 0.004 → **0.003** | 0.005 → **0.003** | 0.007 → **0.004** | — → **0.004** | 0.006 → **0.003** |
| 1 × 300 | 0.008 → **0.005** | 0.014 → **0.005** | 0.019 → **0.013** | — → **0.006** | 0.012 → **0.005** |
| 50 × 3, inbred | 0.009 → **0.008** | 0.021 → **0.008** | 0.029 → **0.022** | — → **0.017** | 0.023 → **0.009** |
| 2,000 × 3 | 0.018 → **0.017** | 0.005 → **0.005** | 0.030 → **0.023** | — → **0.009** | 0.020 → **0.017** |

**Between-locus share F (— where penalised counting reports none: one sample)**

| setting | rich | middling | thin | no loci (17) | the kink (9) |
|---|---|---|---|---|---|
| 1 × 3 | 0.070 → **—** | 0.211 → **—** | 0.478 → **—** | — → **—** | 0.031 → **—** |
| 10 × 10 | 0.014 → **0.005** | 0.036 → **0.009** | 0.086 → **0.026** | — → **0.023** | 0.005 → **0.003** |
| 50 × 30 | 0.009 → **0.006** | 0.009 → **0.003** | 0.054 → **0.023** | — → **0.014** | 0.009 → **0.007** |
| 1 × 300 | 0.043 → **—** | 0.080 → **—** | 0.156 → **—** | — → **—** | 0.029 → **—** |
| 50 × 3, inbred | 0.009 → **0.007** | 0.017 → **0.010** | 0.058 → **0.038** | — → **0.025** | 0.006 → **0.006** |
| 2,000 × 3 | 0.011 → **0.011** | 0.016 → **0.015** | 0.035 → **0.026** | — → **0.012** | 0.011 → **0.011** |

**Allele length variation (share of chromosomes misplaced)**

| setting | rich | middling | thin | no loci (17) | the kink (9) |
|---|---|---|---|---|---|
| 1 × 3 | 0.029 → **0.016** | 0.078 → **0.031** | 0.314 → **0.130** | — → **0.110** | 0.018 → **0.011** |
| 10 × 10 | 0.015 → **0.012** | 0.040 → **0.023** | 0.139 → **0.086** | — → **0.081** | 0.013 → **0.012** |
| 50 × 30 | 0.013 → **0.012** | 0.042 → **0.040** | 0.105 → **0.090** | — → **0.084** | 0.008 → **0.008** |
| 1 × 300 | 0.024 → **0.014** | 0.076 → **0.034** | 0.274 → **0.116** | — → **0.095** | 0.020 → **0.008** |
| 50 × 3, inbred | 0.014 → **0.013** | 0.040 → **0.036** | 0.116 → **0.101** | — → **0.079** | 0.009 → **0.009** |
| 2,000 × 3 | 0.009 → **0.009** | 0.041 → **0.041** | 0.102 → **0.102** | — → **0.101** | 0.016 → **0.016** |

**Thin classes and one sample.** At 1 sample × 3 reads penalised counting roughly halves every
error: the slip rate at thin classes 0.97 → 0.63 relative error, the fall-off 0.31 → 0.17, the
allele length variation 0.31 → 0.13 of chromosomes misplaced. At the well-covered settings the gain
sits in the thin classes and in the numbers a class measures worst: at 50 × 30 the thin classes'
between-locus share 0.054 → 0.023 and fall-off 0.007 → 0.004, while the rich classes were already
measured as well as their loci allow and do not move. **The empty class at 17 repeats** gets a value
from the curve about as good as a thin class's (slip rate within 0.8 to 7.5% from 10 × 10 up).

**The slip rate is where penalised counting gains least**, because each class measures its own
slip rate from every read: even a 5-locus class at 10 × 10 has hundreds of reads behind it. The
blend is as good or better at middling classes (0.054 against 0.055 at 10 × 10; 0.020 against 0.027
at 1 × 300; 0.037 against 0.042 inbred), and penalised counting is better at thin ones (0.071 against
0.075; 0.024 against 0.037; 0.041 against 0.048).

**At 2,000 samples × 3 reads the two fits agree** within a draw's noise on every number (one draw):
each class holds millions of reads and its loci, not its reads, limit it. The allele-length curve's
strength came out at 1.9, almost no smoothing, so the thin classes' allele length variation is the
counting fit's (0.102 both).

### 3.3 The kink

**A well-measured kink survives, and does not cost the rest of the curve its smoothing.** The class
at 9 repeats keeps its own slip rate as well as the counting fit does wherever it is measured:
relative error 0.053 against 0.053 at 10 × 10 (447 slipped reads), 0.018 against 0.018 at 50 × 30
(6,712), 0.035 against 0.032 at 1 × 300 (1,342), 0.032 against 0.031 at 2,000 × 3 (26,849). At
10 × 10, 50 × 30, 1 × 300 and 2,000 × 3 the bend at 9 repeats was eased in every draw, and the one at
10 in all but one. The blend keeps the kink about as well (0.064, 0.018, 0.040, 0.031).

**The inbred setting pulls it most:** 0.047 against the counting fit's 0.032, on 671 slipped reads;
the blend's is 0.056.

**The kink no longer drags the whole curve's smoothing down.** The slip-rate strength is 290 to 680
without the kink and 25 to 610 with it (all but one draw above 100), and the thin classes' slip rate
is as good with the kink as without (1 × 300: 0.024 and 0.023; 10 × 10: 0.071 and 0.072). With the
first tail, 6.25 (§2.3), the kink had pulled the strength down to 5 to 10 and the thin classes' error
up to 0.049 at 1 × 300.

**At 1 sample × 3 reads the kink is not measured, by any method**: 13 slipped reads at 9 repeats,
and the counting fit's own answer there is 0.86 relative error. Penalised counting smooths it toward
its neighbours (0.81), which is the right behaviour for a class whose counts say almost nothing.

### 3.4 Where penalised counting is worse

- **The top of the range at 1 sample × 3 reads.** There the strengths come out at 10⁹ and above in
  most draws, a straight line on the logit scale, because the counts cannot show the curve bending.
  The drawn slip rate bends down at the top, so the straight line runs high: on average +36% at 18
  repeats, +42% at 19 and +49% at 20, where the counting fit is +60%, +17% and +6% on average but far
  noisier (root-mean-square 1.13, 0.90 and 0.61 against penalised 0.59, 0.69 and 0.81). A curve drawn
  through thin evidence leans on the classes that have evidence, and here those are at the bottom.
- **The kink class** in the inbred setting (above), and **the middling classes' slip rate** against
  the blend (§3.2).
- **With the clustering option on** (§2.5), at 2,000 × 3 (one draw): the slip rate at rich classes 0.123 relative
  error against 0.014 with it off (the counting fit's 0.017), the shorter share at thin classes
  0.103 against 0.005, while the allele length variation improves at middling and thin classes
  (0.041 → 0.017, 0.102 → 0.066) and worsens at rich ones (0.009 → 0.018). At 50 × 30 the same option
  halves the thin classes' allele-length error (0.090 → 0.050) and the middling classes' (0.040 →
  0.020), with every stutter number within 0.001 of the default's. It is off by default for the
  2,000-sample result.

## 4. On the shared bench, at the corners

**What was run:** the phase 0 bench (`examples/ng_str_fit_bench.rs`,
[`phase0_bench_2026-10-12.md`](phase0_bench_2026-10-12.md)) with penalised counting wired in as `D`
(`ThePenalisedFit`, the spread read off pairs of samples as method B reads it) and `D_working_F` (the
same fit, reporting the concentration it worked at). The tomato-like truth with PCR stutter (F 0.7,
homozygote excess 0.85 to 0.99 a sample), one-base runs (14 classes, 8 to 21 repeats) and two-base
runs (21 classes, 6 to 27 repeats), at most 100 loci a class, the shipped ±9 window, the bench's
five-fold held-out score; seeds 1 to 3 at one sample, seed 1 at 2,000 samples with 10 held-out loci a
class. B was rerun on the same seeded draws; the climb (A) is in phase 0's report. Recovery figures
are the bench summariser's: the median over classes of each class's root-mean-square error over the
seeds.

**Held-out prediction**: the log-likelihood of held-out loci under each fit minus under the truth,
a locus (zero is as good as the truth):

| cell | B, as fitted | B, after the caller's curves | D (working F at one sample) |
|---|---:|---:|---:|
| one-base runs, 1 sample × 3 reads | −8.0 | −0.57 | **−0.039** |
| two-base runs, 1 sample × 3 reads | −19.9 | −0.91 | **−0.050** |
| one-base runs, 1 sample × 300 reads | −46.0 | −46.0 | **−0.044** |
| two-base runs, 1 sample × 300 reads | −91.3 | −91.3 | **−0.050** |
| one-base runs, 2,000 samples × 3 reads | −428.3 | −428.3 | **−1.35** |
| two-base runs, 2,000 samples × 3 reads | −439.3 | −440.0 | not scored |

**Two things stand behind those differences, and these runs cannot split them.** At 300 reads and at
2,000 samples, most of B's loss is, by phase 0's account, the counting fit's zero-share defect,
which another track is fixing: a length driven to a share of exactly zero costs a held-out locus carrying it hundreds of
units. And at one sample B sits at its concentration bound (F ≈ 0 against a truth of 0.7) where D's
working F is nearer (root-mean-square 0.12 to 0.44), and the held-out score at one sample depends on
F. So these figures say penalised counting predicts held-out loci almost as well as the truth at
these corners; they do not yet say how much of that its curves give over B's.

**Recovery**, the numbers D reports:

| cell | number | B after curves | D |
|---|---|---:|---:|
| one-base, 1 × 3 | slip rate (relative) | 0.167 | **0.098** |
| | shorter share | 0.087 | **0.063** |
| | fall-off | 0.054 | **0.037** |
| | reference length's share | 0.068 | **0.044** |
| | share one repeat shorter | 0.036 | **0.015** |
| two-base, 1 × 3 | slip rate (relative) | 0.140 | **0.096** |
| | reference length's share | 0.070 | **0.044** |
| one-base, 1 × 300 | slip rate (relative) | 0.019 | **0.015** |
| | reference length's share | 0.065 | **0.033** |
| two-base, 1 × 300 | slip rate (relative) | 0.019 | **0.015** |
| | reference length's share | 0.063 | **0.040** |
| one-base, 2,000 × 3 | slip rate (relative) | 0.0067 | 0.0063 |
| | shorter share | 0.0023 | 0.0019 |
| | reference length's share | 0.029 | 0.030 |
| | between-locus share F | 0.030 | 0.029 |
| two-base, 2,000 × 3 | slip rate (relative) | 0.0034 | 0.0028 |
| | reference length's share | 0.052 | 0.052 |
| | between-locus share F | 0.054 | 0.055 |

At 2,000 samples the two agree on every number to within a few parts in a thousand. The two-base
held-out score for D was not run: its fit took 52 minutes under this load, so five held-out refits
would have taken about five hours.

**F is absent from D at one sample** (0 of 14 and 0 of 21 classes answered); B reports it at its
bound, median bias −0.70. **The reference length's share** carries D's one systematic lean, the mass
the floor of §2.5 spreads into the tails: median bias −0.011 to −0.032, against B's +0.001 to +0.018.

## 5. Cost

**Measured with other fits compiling and running on the same 18-core machine; every time is
provisional.** The load average ran from about 50 to 145 during these runs, so the same fit's wall
time varied up to four-fold between runs; the comparison that holds is within one run, where both
fits met the same load in turn.

**On the drawn grid at ±4, a round of penalised counting costs what a round of the counting fit
costs**: each round runs the counting fit's own E-step over every class, and the M-step (every
curve, every strength choice) took 1 to 2% of the fit's time in the final drawn grid, and 28% in the
inbred setting, where the pairs route's likelihood tables are built. **What differs is the number of
rounds.** The counting fit stops each class when that class settles; penalised counting runs every
class until the slowest settles:

| setting | counting fit, s | penalised, s | counting fit, class-rounds | penalised, rounds × classes |
|---|---:|---:|---:|---:|
| 1 sample × 3 reads | 16.1 | 14.7 | 2,079 | 169 × 12 |
| 10 samples × 10 reads | 6.5 | 7.6 | 513 | 41 × 12 |
| 50 samples × 30 reads | 13.4 | 13.5 | 572 | 39 × 12 |
| 1 sample × 300 reads | 2.4 | 3.2 | 519 | 45 × 12 |
| 50 samples × 3 reads, inbred | 131 | 158 | — | 342 × 12 |
| 2,000 samples × 3 reads | 536 | 821 | 1,991 | 243 × 12 |

(Means over the draws, final drawn grid, ±4 window; the counting fit's class-rounds are summed over
its classes; in the inbred setting its rounds are not comparable, being counted for the last
alternation only.)

**On the bench, at the shipped ±9 window and 14 or 21 classes, D was slower than this:** its fit on
all loci took 20 and 102 s at 1 sample × 3 reads (one- and two-base runs) against B's 2.7 and 2.9 s,
4.9 and 7.3 s at 300 reads against 1.3 and 0.6 s, and 992 and 3,120 s at 2,000 samples against 666
and 659 s. Where that time goes was not split; the drawn grid at ±4 put 98 to 99% of it in the
E-step, but at ±9 the allele-length curve has 18 log-ratios a repeat count and the strength search
refits it 30 to 55 times a round while strengths are being chosen, so the M-step's share is likely
larger there. The bench's tolerance (10⁻⁷) and pairs alternations (up to 10, each running to
tolerance) also add rounds that penalised counting runs for every class at once. Peak memory was
lower than B's: 983 MiB against 2,450 MiB for two-base runs at 2,000 samples.

**Memory** is the counting fit's: one class's expected counts at a time, and each curve's banded
matrix, a few hundred kilobytes at most for the allele lengths. The whole final drawn grid, the
2,000-sample draw included, peaked at 0.95 GiB resident.

**At the shipped window of ±9 a round's E-step costs about eight times what it costs at ±4**, for
both fits: 190 genotypes and 17 read buckets against 45 and 9. The overall ratio between the two
fits at ±9 is the bench's, above.

**Where the rounds could be saved**, untried: a class whose numbers have stopped moving could skip
its E-step and reuse its last counts until its curve moves it again.

## 6. Checks on the code

- **The counting fit is unchanged to the bit.** Its round was split out as `expected_counts` and its
  pairs as `SamplePairs`; six fits (allele spans 1, 3 and 4; homozygote excess 0.3, 0.4 and 0.95;
  with and without the pairs route) print the same bits before and after, and its four tests pass.
- **With the penalty near zero (and no floor), penalised counting is the counting fit class by
  class** (test `with_no_penalty_the_fit_is_the_counting_fit`: three classes, every number within
  10⁻⁴).
- **A kink drawn at 0.53 of its curve on 400,000 reads is kept within 5%**, while a 40-read point
  with no slipped read borrows to within 25% of its truth (test, five noisy draws).
- **The banded solve matches a dense one** (test).
- **Every run converged** in the final grid: the largest change in the last round was below 10⁻⁶.

## 7. What is uncertain, and what was not done

- **The drawn grid is one motif length with one set of curves**, at a ±4 window. The bench's
  tomato-like truth (§4) is the second, at ±9, at its corners only.
- **The strength is chosen from the E-step's expected counts** as though they were observed (§2.3);
  where genotypes are uncertain this may choose too little smoothing. Not measured.
- **The strength choice sees one number's counts at a time.** That is why the clustering option
  misleads it (§2.5): the reads' say on the allele spread, through the stutter, is not in the
  allele counts. A joint choice, or the clustering weight applied with a guard, was not tried.
- **No standard errors** (avenue F), so error honesty is not scored.
- **Library stutter multipliers are not fitted**: every library is scored at its class's stutter, as
  in the counting fit.
- **The pairs route** has no unit test of its own; it was run in the inbred drawn setting and on the
  bench.
- **The between-locus share's absence** is decided per motif length (no locus with two samples with
  reads), not per class.
- **One kink, at a rich class.** A kink at a middling class, and two kinks in one curve, were not
  drawn.
- **Seeds 2 and 3 at 2,000 samples on the bench** were not run (each fit there takes 10 to 15
  minutes and the held-out score an hour more under this machine's load).

## 8. What the shared harness should measure next

1. **The corners again after the counting fit's zero-share fix**, to split §4's held-out gap into
   what the fix gives B and what the curves give D; and with F held at the truth for both at one
   sample, to take F out of that comparison.
2. **The grid's middle cells** (10 and 60 samples × 3 to 30 reads), where the drawn grid puts the
   gain on the shorter share, fall-off, F and allele length variation rather than on the slip rate.
3. **A kink at increasing evidence**: one class's slip rate at 0.53 of its curve, its slipped reads
   swept from 100 to 10,000, reporting how much of the kink each method keeps (D, B after the
   caller's curves, the blend). The inbred setting here kept less of it than the counting fit did at
   671 slipped reads.
4. **Whether choosing the strength on expected counts smooths too little at 1 sample × 3 reads**:
   the held-out score at strengths fixed one and two natural-log units above the chosen ones.
5. **The clustering option with a guard**, at 60 and 500 samples, where it may help the allele
   length variation without the stutter paying.
6. **The real cohorts' held-out prediction** (tomato2, kimura), the only measure that works without a
   truth set (plan §2).
7. **Cost at kimura's scale** against the plan's budget of one hour on 64 GiB: per round D costs B's
   E-step, so the question is how many rounds 2,169 samples need (243 at 2,000 drawn samples here,
   against the counting fit's 128 to 308 a class).

## 9. How to reproduce

```sh
# the drawn motif length: [allele span] [draws] [setting] [excess] [pair alternations]
#                         [slip rate at 9 repeats as a share of its curve] [tail] [1: clustering on]
cargo build --release --features bench-fixtures --example str_fit_penalised_em
target/release/examples/str_fit_penalised_em 4 10 1x3 0.3 0 > rows.tsv
target/release/examples/str_fit_penalised_em 4 5 10x10 0.3 0 1.0 > no_kink.tsv
target/release/examples/str_fit_penalised_em 4 3 50x3 0.95 3 > inbred.tsv
target/release/examples/str_fit_penalised_em 4 5 50x30 0.3 0 0.53 2.0 1 > clustered.tsv

# the bench corners
cargo build --release --features bench-fixtures --example ng_str_fit_bench
OUT=tmp/bench_corners METHODS="B D D_working_F" SEEDS="1 2 3" PERIODS="1 2" CELLS="1:3 1:300" \
    LOCI_AT_MOST=100 ./scripts/research/str_fit_grid.sh
OUT=tmp/bench_corners METHODS="B D" SEEDS="1" PERIODS="1 2" CELLS="2000:3" \
    LOCI_AT_MOST=100 HELDOUT_LOCI=10 ./scripts/research/str_fit_grid.sh
uv run --no-project --with numpy python scripts/research/str_fit_summarise.py tmp/bench_corners
```

The drawn rows are seeded and reproduce exactly; the times do not.

| file | what |
|---|---|
| `src/parameter_estimation/joint/ssr_fit/penalised_fit.rs` | the fit, its curves, its tests |
| `src/parameter_estimation/joint/ssr_fit/composite_fit.rs` | `expected_counts` and `SamplePairs`, split out |
| `src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs` | `ThePenalisedFit`, methods D and D_working_F |
| `examples/str_fit_penalised_em.rs` | the drawn motif length of §3 |

---

## 10. Follow-up, 2026-10-13: where the time goes, standard errors, the overshoot

Three questions from the lead after the first report. Same branch; nothing merged from track 1.
**Every time here is provisional and the machine's load swung widely between runs** (one round of
E-steps over every class, run one class after another, took between 2 and 59 ms for the same cell on
different runs). So the comparisons below are made **within one run**, where both fits met the same
load in turn, or in units of that run's own single pass.

**What penalised counting depends on in `composite_fit.rs`**, for the integration branch:
`expected_counts` and `ClassCounts` (one round's E-step), `SamplePairs::of`, `::log_likelihood`,
`::log_likelihood_by_tract` and `::is_empty` (F off pairs of samples). This follow-up added
`log_likelihood_by_tract` and `CountedClass::rounds_in_all`; nothing the counting fit computes
changed (its four tests pass).

### 10.1 Where the time went, and what was cut

**How it was measured.** `examples/str_fit_penalised_cost.rs` draws one bench cell (tomato-like PCR
truth, at most 100 loci a class, the shipped ±9 window, seed 1) and runs, in one process: one round
of E-steps over every class (a *pass*); the counting fit as the bench runs it (B); and penalised
counting with its time split into E-step, strength search, pairs of samples and standard errors. A
*class-round* is one class's E-step in one round.

**Three causes, in order of size:**

1. **Lockstep rounds.** B stops each class when that class settles; penalised counting ran every
   class until the slowest settled. At 2,000 samples × 3 reads on two-base runs that was 1,051
   rounds × 21 classes, 22,071 class-rounds, against B's 8,501.
2. **E-steps one class at a time.** Each class's E-step spread only over its own loci. At one sample a
   class has 100 loci and one sample each, far too little to fill 18 cores; B fits whole classes side
   by side.
3. **The strength search**, re-run for every curve every round: 2 to 4 s a fit at the small cells,
   where the whole fit took 1 to 38 s; and a whole-range search at the start of every counting after
   each reading of F.

**What was changed:**

- **A curve whose strength moves by less than 0.1 natural-log units when chosen keeps it**
  (`strength_steady_within`), and only the fit's first round searches the whole range. The search's
  time fell 4.5- to 16-fold in one run (for example 2.31 → 0.14 s, two-base runs at 1 × 3), and the rounds fell
  where the strength had kept the numbers moving (two-base runs at 60 × 10: 534 → 288).
- **Classes' E-steps run side by side.**
- **A class whose numbers moved by less than the tolerance keeps its last expected counts**
  (`reuse_settled_counts`) until the curves move it again; every class's numbers are still chosen
  every round. Class-rounds counted fell 1.04- to 3.3-fold at the small cells (against the steady
  schedule without it, same run), and 1.6-fold (one-base)
  and 3.0-fold (two-base) at 2,000 × 3, where they are now below B's.

**Accuracy is unchanged:** on the drawn motif length of §3, the reuse rule leaves every number at
every setting identical to the third decimal; the steady strength moves the 1 × 3 setting within its
draws' noise (slip rate at middling classes 0.39 → 0.30 relative error, at the unmeasured kink class
0.81 → 0.92) and nothing else beyond 0.002.

**Before and after**, penalised counting's time as a multiple of B's in the same run. *Before* is the
version of §5 (strengths every round, classes one at a time, every class counted every round);
*after* is the current default.

| cell | B class-rounds | D class-rounds, before → after | D ÷ B time, before | D ÷ B time, after |
|---|---:|---:|---:|---:|
| one-base, 1 sample × 3 reads | 1,723 | 5,082 → 4,085 | 15 | 3.1 to 4.4 |
| two-base, 1 sample × 3 reads | 2,917 | 4,809 → 4,561 | 22 | 2.9 to 8.2 |
| one-base, 1 sample × 300 reads | 1,613 | 448 → 395 | 1.0 | 0.17 to 0.35 |
| two-base, 1 sample × 300 reads | 1,136 | 945 → 837 | 2.4 | 1.1 to 1.7 |
| one-base, 60 samples × 10 reads | 1,320 | 1,470 → 922 | 1.4 | 1.3 to 1.7 |
| two-base, 60 samples × 10 reads | 2,698 | 11,214 → 1,847 | 2.5 | 1.3 to 1.9 |
| one-base, 2,000 samples × 3 reads | 2,741 | 4,298 → 2,642 | 1.5 | 1.06 |
| two-base, 2,000 samples × 3 reads | 8,501 | 22,071 → 7,247 | 4.7 | 0.69 |

The *after* ranges are two runs under different loads; the lower of each pair is from a run with
the steady strengths and side-by-side classes but before the reuse rule. At 2,000 samples the *before* ratio is the
bench's (§4), the *after* ratio is in units of each run's single pass, because B and D ran in separate
processes there. The *before* class-rounds at 2,000 samples are from the steady-strength,
side-by-side run, the closest measured.

**What remains:**

- **At one sample, D still takes 3 to 8 times B's time**, though only seconds in all (0.3 to 8 s a
  fit). Its rounds are short (100 loci of one sample each), so the cost of handing each round's work
  to the threads is most of it: about 2 ms a class-round against B's 1.1 ms. B, which hands a whole
  class's fit to one thread, does not pay it every round. Not changed.
- **At 60 samples × 10 reads, reading F off pairs of samples is most of the time**: 29 to 36 s of
  35 to 48 s in the second run (its 57-point table of the pairwise likelihood for every class, twice).
  B's own search over the same likelihood (about 89 evaluations a class a reading) is in its 20 to 37
  s. A table only near the current value would cut it; not done.

### 10.2 Standard errors

**What they are** (`PenalisedClassErrors`, the doc comment states it in full). Write `A` for the
counting objective's curvature in the curves' unbounded values, `B` for the spread of its slope from
locus to locus, and `P` for the bend penalties at the strengths and bend weights fitted. A penalised
fit's Bayesian covariance would be `(A + P)⁻¹` if the objective were a true likelihood; it is not,
because it treats samples at one locus as independent, so `A` overstates what the reads say. The
covariance reported is `H⁻¹ (B + P) H⁻¹` with `H = A + P`: `B` in place of `A`, as the sandwich formula
does for an unpenalised fit (track 1), with the penalty's part kept, which is the allowance a
Bayesian interval makes for the curve's own bias. With `B = A` it is the usual Bayesian covariance.

- `A` is taken numerically, class by class: the change in the class's slope over a step of 0.001 in
  each of its numbers either way, the slope from the expected counts by Fisher's identity, exactly as
  track 1 computes the counting fit's.
- `B` sums each locus's slope times itself: the locus is the unit that repeats independently.
- F read off pairs of samples answers a second objective, so it gets the same formula over its own
  curve, its slope and curvature taken numerically locus by locus from the pairwise likelihood; the
  other numbers' errors are taken with F held, so they leave out F's uncertainty. The pseudo-draw
  floor is left out of `H` and `P`.
- At one sample F is absent, and so is its error; the other numbers are conditional on the working F.

**Honesty on the bench**: of the estimates with an error, the share whose truth lay within one and two
errors (honest: about 0.68 and 0.95). Tomato-like PCR truth, at most 100 loci a class, seeds 1 to 6
at the small cells, seed 1 at 2,000 samples (the 2,000-sample runs were made on the build before the
reuse rule, which on the drawn grid changes no number at the third decimal).

| cell | slip rate | shorter share | fall-off | F | shares ≥ 0.05 | shares 0.01–0.05 | shares < 0.01 |
|---|---|---|---|---|---|---|---|
| one-base, 1 × 3 | 0.69 / 0.98 | 0.56 / 0.90 | 0.58 / 0.93 | absent | 0.67 / 0.91 | 0.70 / 0.89 | 0.57 / 0.88 |
| two-base, 1 × 3 | 0.68 / 0.90 | 0.49 / 0.83 | 0.60 / 0.86 | absent | 0.66 / 0.90 | 0.71 / 0.93 | 0.56 / 0.88 |
| one-base, 1 × 300 | 0.67 / 0.92 | 0.57 / 0.85 | 0.51 / 0.90 | absent | 0.69 / 0.89 | 0.67 / 0.93 | 0.63 / 0.90 |
| two-base, 1 × 300 | 0.64 / 0.95 | 0.65 / 0.96 | 0.77 / 0.98 | absent | 0.60 / 0.92 | 0.64 / 0.88 | 0.58 / 0.88 |
| one-base, 60 × 10 | 0.62 / 0.92 | 0.63 / 0.99 | 0.77 / 0.96 | 0.63 / 0.95 | 0.67 / 0.94 | 0.53 / 0.68 | 0.27 / 0.34 |
| two-base, 60 × 10 | 0.72 / 0.95 | 0.85 / 0.98 | 0.63 / 0.94 | 0.66 / 0.94 | 0.64 / 0.89 | 0.41 / 0.52 | 0.20 / 0.25 |
| one-base, 2,000 × 3 | 0.57 / 0.86 | 0.43 / 0.86 | 0.64 / 0.86 | 0.64 / 1.00 | 0.57 / 0.90 | 0.55 / 0.67 | 0.23 / 0.31 |
| two-base, 2,000 × 3 | 0.81 / 0.90 | 0.81 / 1.00 | 0.76 / 1.00 | 0.57 / 0.90 | 0.66 / 0.90 | 0.36 / 0.52 | 0.22 / 0.26 |

Estimates a cell: 84 or 126 for each stutter number and F at the small cells, 14 or 21 at 2,000
samples; 294 to 1,296 shares a group.

**The stutter numbers, F and every share of 0.05 or more are honest within the noise of these
counts**, the shorter share the least so (0.43 to 0.57 within one error at four cells). **Shares
below 0.05 are not, from 60 samples up**: a share of 0.004 is reported as about 10⁻⁶ with an error of
about 10⁻⁵, 10 to 100 errors below the truth. At a between-locus share of 0.7, most loci are fixed at
one length, so a share of 0.004 at 100 loci is carried by a fraction of one locus. In most draws no
locus carries the length, the estimate falls nearly to zero, and an error read off the loci's own
spread cannot see a length that no locus carried. Track 1 reports the same limit for the counting
fit's sandwich (shares below 0.01 at 60 samples: 36 within one error, 50 within two); ours is
worse (20 to 27, and 25 to 34) because at 60 samples the floor of one allele draw a class is
negligible against about 12,000 draws, so nothing holds the estimate off zero. At one sample the
floor is large enough that small shares cover (0.56 to 0.63 / 0.88 to 0.90). **A share's error should
be flagged where the share times the loci is below a few loci**, as track 1 recommends; a model-based
error for such a share, or an interval on the log-ratio scale, is not tried.

**Cost of the errors**: `2k + 1` E-step passes a class (`k` numbers: 22 at ±9, with F) plus one pass
locus by locus, then dense matrices of up to about 480 rows. On the final build at the bench's small
cells the whole fit with errors took 0.1 to 4.0 s (the summariser's fit time, seeds 1 to 6); in the
cost runs, on an earlier build that took the classes one at a time, the errors took 2.8 to 21 s
against fits of 0.4 to 48 s.

### 10.3 The overshoot at 1 sample × 3 reads

**The question:** at 1 sample × 3 reads the drawn slip rate's curve ran high at the top of the
range (+49% at 20 repeats, §3.4). Is the strength chosen from expected counts too weak, too strong,
or is it the rule that holds a curve flat beyond the fitted range?

**It is not the hold-flat rule.** That rule applies only beyond the highest class with reads, and
the classes at 18, 19 and 20 repeats all have reads (12, 8 and 5 loci).

**It is not a wrong strength either: no strength removes it without adding an error elsewhere.** Every
curve was held at fixed strengths from 1 to 10⁶ on the same ten draws (`str_fit_penalised_em`, ninth
argument). Slip rate, root-mean-square relative error:

| strength | rich | middling | thin | at 18 | at 19 | at 20 |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 0.31 | 0.52 | 0.73 | 0.89 | 0.73 | 0.59 |
| 10 | 0.25 | 0.45 | 0.59 | 0.64 | 0.58 | 0.55 |
| 100 | 0.20 | 0.35 | 0.59 | 0.59 | 0.65 | 0.72 |
| 1,000 | 0.18 | 0.20 | 0.66 | 0.64 | 0.79 | 0.96 |
| 10⁴ | 0.19 | 0.19 | 0.72 | 0.69 | 0.88 | 1.08 |
| 10⁶ | 0.19 | 0.19 | 0.73 | 0.69 | 0.89 | 1.10 |
| chosen (mostly 10⁹) | 0.21 | 0.39 | 0.63 | 0.59 | 0.69 | 0.81 |

A strong curve suits the rich and middling classes (0.19) and fails the top (1.08 at 20 repeats); a
weak one does the reverse. The chosen strength lies between. **Two things combine at the top:**

- **the classes' own counts are biased high at one sample**: the counting fit alone, class by class,
  is +50 to +60% at 16 and 18 repeats on average; a weak curve follows them;
- **the drawn truth bends down on the logit scale at the top**, and a strong curve's straight line
  through the well-measured bottom classes cannot follow the bend it cannot see.

**On the bench's truth the same thing appears at the other end.** At 1 × 3 the top of the range is
nearly unbiased (within ±10% from 15 repeats up, both motif lengths), and the bottom class runs high: +86% at
8 repeats on one-base runs and +127% at 6 on two-base runs, where the slip rate is lowest (0.005 to
0.0075) and a class's few slipped reads cannot show where the curve bends. The errors cover it
within two errors in 6 draws of 6 on one-base runs, but in only 2 of 6 on two-base runs.

**The limit, stated:** with one strength for the whole curve, at a corner where the end classes'
counts are biased and too few to show a bend, penalised counting extrapolates the shape of the
classes that have evidence. The curve at the ends of the range at one sample × three reads is the
least trustworthy part of the fit, and its errors widen there (relative error 0.20 to 0.68
at the end classes against 0.08 to 0.16 mid-range on the bench) but do not always cover it. What
might help, untried: a strength that varies along the curve (an adaptive penalty), or a first-order
penalty near the ends so that a curve flattens rather than extends its slope where its evidence ends.

### 10.4 What the shared harness should measure next, revised

To §8's list:

1. **The integration branch's head-to-head** against the fixed counting fit, held-out and recovery,
   now that both have errors.
2. **Small shares' errors**: a model-based or log-ratio interval for a share carried by fewer than a
   few loci, scored as here.
3. **Cost at kimura's scale** with the reuse rule, against the plan's budget of one hour on 64 GiB:
   at 2,000 drawn samples D counted fewer class-rounds than B.
