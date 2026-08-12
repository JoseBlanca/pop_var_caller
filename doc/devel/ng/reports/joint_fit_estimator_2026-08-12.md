# The per-site route: what the selection does on real genomes, and which estimator to fit with

*Research report, 2026-08-12. Covers `spec/parameter_prepass_joint_fit.md`,
`parameter_prepass_joint_loci.md`, `parameter_prepass_joint_records.md` and their three arch
documents. **Two programs stand behind it**: `examples/ng_joint_loci_probe.rs`, which runs the
selection over a reference and its repeat catalog, and `examples/ng_joint_fit_harness.rs`, which
fits the estimator against a truth it drew itself. The first unit of code also landed —
`src/ng/parameter_estimation/joint/loci.rs`.*

---

## 1. What was asked and what came back

Three things, in the order they were done.

1. **An adversarial review of the six documents**, acting on what it found rather than listing it.
   Fifteen findings; the three that changed the design are §3, §4 and §5 below, and the rest are in
   §7.
2. **The locus selection built and measured on real genomes** — §2. It was the cheapest thing
   outstanding and it closed an open question three decisions were waiting on.
3. **The estimator settled and measured against a known truth** — §3 to §6.

---

## 2. The selection, measured on tomato and on GRCh38

`src/ng/parameter_estimation/joint/loci.rs`, driven by `examples/ng_joint_loci_probe.rs`.

### 2.1 The generic rule does what the spec says

Keep position `p` when `hash(contig, p, seed) < threshold`, over the analysed regions.

| | tomato SL4.00 | GRCh38 + hs38d1 |
|---|---:|---:|
| reference bases | 782,520,033 | 3,105,715,063 |
| positions kept for a target of 2,000,000 | 2,002,505 | 1,999,981 |
| gaps between neighbours, p10 / median / p90 | 41 / 271 / 900 bp | 155 / 1,019 / 3,385 bp |
| the geometric prediction at that density | 41 / 271 / 900 | 155 / 1,019 / 3,385 |
| per-contig counts, chi-square | 19.8 on 12 d.f. | 127.9 on 127 d.f. |
| a sweep over shards cut at arbitrary coordinates | identical | identical |

**The gap distribution matching a geometric to three figures is the sharper test of the hash's
uniformity than the contig counts are**, and it passes at both genomes and both densities.

### 2.2 Assembly gaps were not excluded, and on human that is 5% of the budget

Nothing in the rule looked at the reference's bases. A position inside an `N` run is kept by a rule
that cannot see it and covered by no read in any sample.

| reference | not `A`/`C`/`G`/`T` | kept positions landing there |
|---|---:|---:|
| tomato SL4.00 | 44,731 bases (0.01%) | **135** |
| GRCh38 + hs38d1 | 165,046,090 bases (5.31%) | **106,423** (5.32%) |

It is not merely wasted budget. The fit derives each sample's heterozygosity as a mean of genotype
posteriors over the kept loci, and a locus with no reads contributes its **prior** — the model's own
prediction — rather than evidence. **Fixed**: the domain is now the analysed regions intersected with
the unambiguous bases, the mask riding on the same forward pass over the FASTA that computes the
reference's digests, and the threshold's denominator is the masked length.

### 2.3 The STR set needs no sampling at all

At the STR path's calling floors `[8, 6, 6, 6, 5, 4]`, a 30 bp flank and a 100 bp satellite cap,
tomato holds **462,701 STR loci in 141 strata**. A cap above the largest stratum keeps every one of
them, which is under a quarter of the generic budget.

| cap | loci kept | strata capped, of 141 |
|---:|---:|---:|
| 100 | 8,699 | 68 |
| 1,000 | 41,271 | 21 |
| 20,000 | 157,752 | 3 |
| none | 462,701 | 0 |

One stratum — period 1 at 8 repeats — holds 217,812 loci, 47% of the total; **68 strata hold fewer
than a hundred each**, so a cap does nearly all its work on three strata and those three are where
the parameter is already best determined.

**Three consequences**, and the third is the one nobody had priced. Sampling is unnecessary, so the
per-stratum reweighting never has to fire; the thin strata still borrow from their neighbours, which
selection was never going to fix; and **holding 462,701 loci in every sample makes the STR records
the larger half of the cohort's memory bill rather than the smaller**, which is the first reason the
cap mechanism has ever had to exist.

*Still open: the same count on GRCh38, which the catalog file on hand cannot answer — it was written
in an older header format. It matters for the GIAB arm of the comparison and for nothing else.*

---

## 3. The estimator: a population frequency, not a panel count

**The specification wrote one model and named another**, and the two are not notations for the same
thing.

- **A count in the panel.** The `2N` chromosomes carry `c` copies of the allele, `c` drawn from a
  spectrum. Given `c` the samples' genotypes are not independent — they must add up to `c` — and the
  sum over configurations is a convolution. This is ANGSD's `realSFS` (Nielsen et al. 2012), and its
  virtue is that the conditional distribution of the genotypes given `c` contains no frequency at
  all.
- **A frequency in the population.** The population carries the allele at frequency `f`, drawn from a
  density; given `f` the samples are independent draws.

The document's formula was the second and its citation was the first.

**Decision: the frequency, and inbreeding forces it.** The cancellation that makes the count form
work needs each individual's genotype to be a pair of independent draws at `f` — that is what makes a
sample's weight `f^j (1−f)^(P−j)` times a combinatorial constant, so `f` divides out. Under an
inbreeding coefficient a diploid heterozygote's weight is `2f(1−f)(1−F)` and each homozygote's
carries an extra `F·f(1−f)`, which is not of that form. This route requires one inbreeding
coefficient per sample and one contamination fraction per sample, and contamination needs to know
*which* allele the population carries rather than only how many copies of it there are. Neither is
expressible in the count form. (It is why the estimators that fit inbreeding from genotype
likelihoods — ngsF, Vieira et al. 2013 — work per site rather than per count.)

**A claim made in review and withdrawn here.** An earlier statement of this finding said the
frequency form throws information away by treating the samples as independent. It does not: the
correlation between two samples at a locus is *induced* by the shared unknown frequency, and
integrating over it keeps all of it. The frequency form is an exact likelihood under its own
generative story, not an approximation to the count form.

**The real hazard is the other one, and it is why the density has four numbers.** A weight per allele
*count* is a distribution over an almost-observable quantity. A weight per *frequency* on the same
`2N + 1` grid is not — a frequency reaches the data only through the genotypes it draws — so
recovering it means undoing a binomial blur, and the maximiser of an unregularised mixing
distribution is discrete (Lindsay 1983). So what is fitted is

```text
π(f) = p_invariant·[f = 0] + p_fixed_alt·[f = 1] + (1 − p_invariant − p_fixed_alt)·Beta(f; a, b)
```

with the integral over `f` a fixed quadrature whose node count is accuracy and not freedom.

---

## 4. The harness, and what it measures

`examples/ng_joint_fit_harness.rs`. The method is `ng_multilib_key_harness.rs`'s, one level up: an
estimator maximises a sum over **the patterns a locus can show** — every sample's depth and
alternative-read count together — so replacing the observed pattern counts with each pattern's exact
probability under the truth gives the objective the estimator climbs with an infinite genome. Its
maximiser is the value the estimate converges to, and the gap from the truth is bias computed rather
than sampled.

Two arms, **one objective**: the exact arm enumerates the pattern space (feasible at two or three
samples) and the drawn arm supplies empirical weights instead (fifty samples). Three candidates:
the adopted four-number density, the same form with one free weight per frequency in a `2N + 1`
grid, and the count form with the convolution — the last fitted with inbreeding at zero, since it
cannot carry it.

**One diagnostic is load-bearing and was missing from the first version: the objective evaluated at
the truth.** A fitted score below it is an optimiser that stopped early, which is a different failure
from a biased estimator and is indistinguishable from one in a table of parameter errors.

---

## 5. What the fits return

### 5.1 Exactly, at two samples and three reads a site

| | truth | fitted | error |
|---|---:|---:|---:|
| clean error rate | 1.895 × 10⁻³ | 2.030 × 10⁻³ | **+7.1%** |
| noisy error rate | 5.29 × 10⁻² | 8.61 × 10⁻² | **+62.7%** |
| noisy-locus share | 0.0088 | 0.0040 | **−54.3%** |
| diversity `Hexp` | 1.536 × 10⁻³ | 1.487 × 10⁻³ | **−3.2%** |
| inbreeding, truth 0.6 | 0.600 | 0.620 | **+3.3%** |

**The last two rows are the estimator working and the middle two are not bias.** The objective at the
true values is −0.113529 and the fit reached −0.113531 — the truth is the higher point, so the
maximiser is where it should be. The middle two rows are a ridge so flat that a climb stops
2 × 10⁻⁶ nats per locus short with the noisy class's rate wrong by two thirds and its share wrong by
half, the two trading against each other. Over two million loci that gap is about 4 nats, so it is
information rather than nothing — it is just spread very thin at two samples.

### 5.2 At ten samples, from drawn data, all three candidates agree

40,000 loci, three reads a site.

| candidate | clean rate | noisy rate | noisy share | `Hexp` | inbreeding (truth 0.6) |
|---|---:|---:|---:|---:|---:|
| four-number density | −8.9% | −29.5% | +72.6% | **+2.4%** | +6.7% |
| free weight per frequency | −8.7% | −30.4% | +77.9% | **+2.4%** | +6.7% |
| count form (inbreeding 0) | −9.1% | −31.6% | +81.7% | **+2.7%** | — |

**The three descriptions of the frequency are interchangeable for every parameter a caller reads.**
`Hexp` differs between them by 0.3 percentage points and the error rates by under 2. What separates
them is the fitted description's own shape: the largest single component holds 25.6% of the
segregating mass under the free grid against 19.0% under the count form, so the free grid is the
spikier of the two — mildly, at ten samples, where it has 21 cells.

**So the case for the four-number density is not accuracy. It is that it is the only one of the three
that carries what this route has to emit**, plus a shape that cannot degenerate.

### 5.3 When the truth's density is one the fitted shape cannot reach

Ten samples, 40,000 loci, the segregating frequencies drawn from two components at 0.2 and 0.8 in
equal parts — the shape two diverged subpopulations leave behind, and one a single Beta cannot make.

| truth's density | clean rate | noisy rate | noisy share | `Hexp` |
|---|---:|---:|---:|---:|
| one bump (in the fitted family) | −8.9% | −29.6% | +72.9% | **+2.4%** |
| two bumps (outside it) | −10.9% | −33.9% | +96.4% | **+7.0%** |

**The misspecification costs about five percentage points on the diversity and about two on the error
rates.** That is the answer to *"is a single Beta enough shape for a landrace panel"*: enough for a
caller's prior, where an error of a few percent in a prior changes nothing a caller does; **not
enough to publish as a site-frequency spectrum**, and the spec says so. The replacement if that
changes is a mixture of two Betas — two more numbers and the same quadrature — rather than a weight
per allele count.

### 5.4 A second bump buys nothing, and the instability is somewhere else

§5.3 leaves an obvious next move: give the density a second bump so it can reach a two-subpopulation
truth. Three more fitted numbers — two shape numbers and the second bump's share. Ten samples, 40,000
loci, each shape fitted from three starting points and the best-scoring fit taken.

| truth's density | four numbers | seven numbers | likelihood gap |
|---|---:|---:|---:|
| one bump, inside the fitted family | `Hexp` **−0.8%** | **+0.4%** | 2.5 × 10⁻⁵ nats/locus |
| two bumps, outside it | `Hexp` **+4.9%** | **+5.8%** | 1 × 10⁻⁶ nats/locus |

**The extra numbers do not recover the two-bump truth better.** They cannot do worse in likelihood —
an extra parameter never can — and here they do not do better in any quantity a consumer reads.

**And the start-dependence they were expected to add is not what dominates.** Across the three starts,
`Hexp` spans 11.4% of its own value under four numbers and 12.8% under seven on the one-bump truth;
on the two-bump truth, 7.3% under four and **2.6%** under seven. So the second bump adds about a
percentage point on one truth and removes five on the other — noise, beside the spread that is
already there.

**What is already there is the two-class noise model, which is a mixture too.** The same three starts,
the same sizes, under both shapes:

| where the search began | clean error rate | `Hexp` |
|---|---:|---:|
| the two noise classes far apart (`ε_clean` 5 × 10⁻⁴, `ε_noisy` 2 × 10⁻¹) | **−1.3%** | −0.8% |
| a middling start (6 × 10⁻³ and 2 × 10⁻²) | −8.9% | +2.4% |
| **the two classes close together** (2 × 10⁻² and 6 × 10⁻²) | **−45.8%** | +10.6% |

**A start that puts the two noise classes near each other collapses them into one and reports
convergence**, costing 46% of the clean error rate — which is exactly the failure
[`parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md) §6.5 records for its own
two-state model, and it is now measured on this route rather than inherited.

**Two decisions follow.** The density keeps four numbers, and the spec's open question about its
shape is closed by measurement rather than by argument. And **the starting points are not an
implementation detail**: three spanning the separation between the classes, best score taken, is part
of the estimator.

### 5.5 How many samples each parameter needs

The adopted candidate, 40,000 drawn loci, three reads a site, one inbreeding coefficient shared by
the panel. Errors against the truth; the inbreeding column is the fitted value itself.

| samples | clean rate | noisy rate | noisy share | `Hexp` | `F`, truth 0 | `F`, truth 0.6 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | −35.0% | −62.6% | +474% | +47.2% | **0.097** | **0.530** |
| 2 | −22.0% | −56.2% | +438% | +2.8% | 0.000 | **0.437** |
| 5 | −21.7% | −59.0% | +383% | +14.2% | 0.000 | 0.607 |
| 10 | −8.9% | −29.6% | +73% | +2.4% | 0.039 | 0.640 |
| 25 | **−1.4%** | **−0.9%** | **+5.4%** | −12.2% | 0.000 | 0.613 |
| 50 | **−0.7%** | +3.2% | **−2.3%** | −2.5% | 0.000 | 0.601 |

**The noisy-locus class is the estimator's weak point and it needs about twenty-five samples.** Below
ten its share is wrong by a factor of four or five — not noisy, *wrong*, and in the same direction at
every panel size, which is what §5.1's flat ridge looks like when there is too little data to pin it.
Between ten and twenty-five it goes from wrong by three quarters to right within a twentieth. **The
error rates follow the same curve**, reaching one part in a hundred at twenty-five samples.

**Inbreeding needs about five samples, and at one it is not identified at all** — 0.097 where the
truth is zero and 0.530 where the truth is 0.6, pulled towards the middle from both sides because one
individual's heterozygote deficit and a density concentrated near zero are the same observation.
That is the measured version of what spec §6.1 claims, and it is why the route must emit
`F_hom_excess` as *not identified* at one sample rather than as a fitted number.

**The `Hexp` column says nothing about panel size and should not be read as if it did.** At 40,000
loci and a segregating share of 0.008, about 320 loci carry the whole of that estimate, so its scatter
here is the *locus budget* rather than the sample count. What sets `Hexp`'s precision is
[`parameter_prepass_joint_loci.md`](../spec/parameter_prepass_joint_loci.md) §4.3's downward budget
sweep, and this table is not it.

### 5.6 How the errors respond to the number of loci

Fifty samples, three reads a site, a fresh draw at each budget, three starts, best taken.

| positions | segregating | clean rate | noisy rate | noisy share | `Hexp` | seconds |
|---:|---:|---:|---:|---:|---:|---:|
| 5,000 | 39 | −1.1% | −2.0% | +15.9% | +14.0% | 5 |
| 20,000 | 149 | −1.4% | +2.8% | −13.3% | −3.5% | 17 |
| 80,000 | 627 | +0.2% | −0.4% | **−1.3%** | −7.0% | 66 |
| 320,000 | 2,614 | −0.4% | −1.3% | −1.0% | **+1.6%** | 273 |
| 1,280,000 | 10,208 | −0.0% | +0.2% | −0.9% | −0.9% | 1,103 |

**Three parameters, three different curves, and pooling them would have hidden all of it.**

- **The error rates are finished at five thousand positions.** They are fitted from *reads*, and five
  thousand positions at fifty samples and three reads is 750,000 read observations. Two hundred and
  fifty-six times the budget buys one further percentage point.
- **The noisy-locus share needs about eighty thousand.** It is a property of loci rather than of
  reads, so loci are what it wants: wrong by a sixth at five thousand, within a hundredth at eighty
  thousand.
- **The diversity is the slowest and it tracks the segregating count**, not the budget — several
  percent of scatter until a few thousand sites segregate, near 1% at ten thousand.
- **Inbreeding needs about twenty thousand.** Against a drawn truth of 0.600 the fit returns 0.563 at
  five thousand positions and 0.600, 0.607, 0.596, 0.603 at every budget above it. The same run at a
  truth of zero returns 0.052 and then 0.000 — so the low-budget failure is an inbreeding coefficient
  invented out of scatter, in both directions.

*The whole table repeats at an inbreeding of 0.6 with every column within a percentage point of the
outbred one, so none of these budgets is a function of how inbred the panel is.*

**So the estimates in this table are satisfied at 320,000 positions, a sixth of the two million the
census budget sets.** What holds the budget at two million is the parameter that is not in the table:
contamination wants about ten thousand segregating markers, and this panel yields 10,208 of them at
1.28 M positions. **Two million is a contamination budget.**

### 5.7 A bigger panel does not buy fewer loci

The same sweep at two hundred samples, four times the panel:

| positions | noisy share, 50 | at 200 | `Hexp`, 50 | at 200 | clean rate, 50 | at 200 |
|---:|---:|---:|---:|---:|---:|---:|
| 5,000 | +15.9% | **+16.0%** | +14.0% | −9.0% | −1.1% | +0.7% |
| 20,000 | −13.3% | **−12.4%** | −3.5% | +3.3% | −1.4% | −0.0% |
| 80,000 | −1.3% | +4.5% | −7.0% | −4.0% | +0.2% | +0.5% |
| 320,000 | −1.0% | −2.8% | +1.6% | +0.4% | −0.4% | −0.2% |

**The noisy-locus share does not feel the panel at all** — +16% at five thousand positions at both
sizes. It cannot: `w` is the share of *loci* that are noisy, so it is counted in loci. The larger
panel buys the read-driven parameters at small budgets instead — the clean rate from −1.1% to +0.7%
at five thousand positions, and the inbreeding coefficient from 0.052 to 0.004 against a truth of
zero.

**So a cohort of thousands does not relieve the locus budget**, and the measurement that might have
said otherwise says it does not either. Counting both — sites segregating in the *population*, which
is a property of the truth, and sites segregating in the *panel*, which is what contamination can
use:

| panel | share of the population's segregating sites that the panel sees |
|---:|---:|
| 50 | 77–82% |
| 200 | 82–85% |
| 1,000 | **89–92%** |

**Twenty times the panel buys about a seventh more usable markers.** Under this truth's density —
`Beta(0.3, 1.2)`, a neutral rare-allele pile-up — about a quarter of segregating sites sit below one
in a hundred and most of those are already visible to fifty samples. A population with a heavier rare
tail would answer differently, which is a reason to re-run this on the real cohort rather than carry
the number as a constant.

*A claim of §5.7's, refined by the same run.* At five thousand positions the noisy-locus share is
+16% at fifty, two hundred **and** a thousand samples — the panel buys nothing there. At twenty
thousand it is −13.3%, −12.4% and −4.0%, so a thousand-sample panel does help. **Samples are not
useless to it; the budget still has to be counted in loci.**

**The fit's cost is almost independent of the panel, which is why this sweep was affordable.** At a
locus, every sample with no alternative read contributes a factor depending only on its depth and the
inbreeding coefficient, so where the panel shares one coefficient the product over the samples is a
table of a dozen numbers raised to counts. Only the samples that showed an alternative read are
scored one at a time — about six in a thousand at three reads and an error rate of 0.002, plus
whoever really carries the variant. The fast path is checked against the slow one on every run and
agrees to twelve decimal places.

---

## 6. What is settled and what is not

**Settled by measurement.**

- The generic selection's arithmetic, scatter and order-independence, on two real references (§2.1).
- That assembly gaps had to be excluded, and what they cost when they are not (§2.2).
- That the STR set fits whole, and what that does to the memory bill (§2.3).
- That the estimator is consistent — the objective's maximum is at the truth (§5.1).
- That the choice between the three frequency descriptions does not move any parameter a caller reads
  (§5.2), so the choice rests on what each can carry rather than on accuracy.
- That the four-number density's shape misspecification costs ~5% of `Hexp` and ~2% of the error
  rates on a two-subpopulation panel (§5.3).

**Open, with the measurement named.**

- **How many samples each parameter needs** — the panel sweep at 1, 2, 5, 10, 25 and 50 samples. The
  noisy class's two numbers are the ones to watch; `Hexp` and inbreeding are already within a few
  percent at ten.
- **How many starting points, spanning what.** §5.1's flat ridge is the reason this matters: three
  starts spanning the separation between the clean and the noisy class are in the harness, and what
  is not yet known is whether more of them close the gap or whether the ridge simply needs more data.
- **Whether the alternative allele being summed over rather than chosen matters in practice.** The
  harness is biallelic today, so the test that would catch it — a cohort with no real variant
  anywhere, where a fit that picks the largest of three error counts inflates the rare classes — needs
  the four-allele emission first.
- **Everything about the STR path's per-locus length frequencies.** The design mirrors the generic
  path (a latent drawn from a fitted per-stratum prior, with one concentration number saying how
  monomorphic loci are), and no code has been written for it.
- **The two-class residual on real data**, which is the comparison's headline measurement and now has
  three arms rather than one (spec §8).

---

## 7. The rest of the review, and what was done about it

| finding | what was done |
|---|---|
| `Hobs` is derived from posteriors whose prior is the fitted density and `F`, so `1 − Hobs/Hexp` restates the fitted parameter rather than checking it; at three reads a site about a fifth of kept loci carry no usable data | spec §5.2 rewritten to say which half of the circularity opens; §3.2 now requires the two evidence counts to be emitted beside the rates |
| the decisive measurement was specified on HG002 alone, where the mechanism being tested does not exist | rewritten as three arms with what each can show — the GIAB trio, a drawn sweep in sample count, and the tomato cohort |
| the STR path reintroduced the per-locus parameter the generic path rejects | rewritten: the locus's length frequencies are a latent drawn from a fitted per-stratum Dirichlet, with the concentration saying how monomorphic loci are; the per-stratum model is its large-concentration limit |
| the third site class needs a per-sample coverage-by-window summary no document held, and its grain is (locus, sample) rather than the locus | records spec §4 now specifies the object and sizes it — 1.6 MB per sample on tomato, 80 MB across fifty; fit spec §2.2 states the grain |
| choosing each site's alternative allele from the data biases the rare end of the density | fit spec §3.1.1: the three non-reference bases are summed over with an equal prior, and the arch doc forbids an allele parameter in the likelihood's signature |
| neither the depth ladder nor the per-position depth cap travelled with the records | records spec §5: thirteen identity values, five of which say in what units the evidence was written down |
| the STR records were priced without a locus count | records spec §6 recomputed against 462,701 loci; the STR set is the larger half |
| `F_hom_excess`'s sign is load-bearing and was never constrained | spec §5.1 and the arch newtype constrain it to `[0, 1]`, with the reason |
| the route's behaviour at one sample was never stated | spec §6.1: it degenerates to the per-sample estimator; the class posterior, the inbreeding coefficient and contamination come back as not identified |
| `SelectionIdentity` could not derive `Eq`, and the `verifyBamID2` error table's units are unverified | both noted where they occur |
