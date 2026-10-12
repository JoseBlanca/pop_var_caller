# The counting fit: honest standard errors, and two defects fixed on the way

**Status:** track 1 of [`../str_fit_research_plan.md`](../str_fit_research_plan.md) (avenues B and
F), 2026-10-12, branch `str-fit-b`. Measured on classes drawn by the fit's own per-class drawing
helper (`ssr_fit::bench_fixtures::draw_stratum`) and on the shared fitting bench of phase 0
([`phase0_bench_2026-10-12.md`](phase0_bench_2026-10-12.md)). **Every time quoted is provisional**:
other agents' builds and fits ran on the same machine throughout, at load averages of 55 to 135 on
18 cores.

---

## 1. The answers

**1. Both error methods are honest wherever the fit is unbiased and the number is carried by more
than a handful of tracts.** On drawn classes with more than one sample and 100 or 200 draws a
setting, the truth lay within one reported error in 60 to 78 draws of 100 and within two in 86 to
98 of 100, for every number, by the bootstrap over tracts and by the sandwich formula alike (§4).
Honest errors give about 68 and 95. The settings run on 40 or 60 draws range from 50 to 80 and from
80 to 100, within their wider noise, apart from shares below 0.01 (§5.1) and one share biased by alleles past the window (§5.2). At one sample both
methods cover every number.

**2. The plain curvature is not honest for the allele-length shares and the concentration.** At 20
samples a tract and a between-locus share of two thirds it holds the truth within one error in 16
to 37 draws of 100: its share errors are about a quarter of the true spread, its concentration
errors about 40% of it. At 20 samples its stutter errors are nearly right; at 500 they are too small
as well (the fall-off within one error in 20 draws of 100).

**3. Where the errors are not honest, the cause is outside the error method** (§5):

- **a share that few tracts carry** — on the fitting bench, shares below 0.01 at about 70 tracts a
  class hold the truth within one error in only 31 to 36 estimates of 100;
- **bias in the fit**, which no error covers once the error is small: on the bench at 2,000
  samples the slip rate of long tracts came out 0.5 to 2.4% low, up to 3.3 errors from the truth.
  Alleles past the fit's outermost length are not the cause at the precision a drawn setting
  reaches; the cause is open (§5.2).

**4. Two defects of the counting fit, found by phase 0, are fixed** (§6, §7):

- **allele-length shares of exactly zero**, which a held-out locus at that length cannot survive:
  at 2,000 samples the held-out score went from 385 log-likelihood units a locus below the truth to
  0.56 above it;
- **at one sample the fit read the concentration from pairs of samples, of which there are none**,
  returned its search bound and fitted the stutter under it: a slip rate a third of the truth at
  three reads. It now reads the concentration from each sample's own two copies, and returns it as
  absent where those cannot place it either.

**5. Cost** (§8, provisional): the sandwich costs 0.01 to 0.3 times what the fit itself costs; a
bootstrap of 50 replicates 0.25 to 9 times.

**Recommendation (§9):** the sandwich as the counting fit's errors; the bootstrap kept as a check;
the shares' errors flagged where few tracts carry the length.

---

## 2. The two methods

**Why the errors need care.** The counting fit scores each sample at a tract as if it drew its
alleles from the class's allele-length variation independently of the other samples there. In
truth samples at one tract share that tract's allele frequencies. The estimates are still right on
average — each sample's own term is exact — but the objective's curvature counts every sample as
fresh evidence about the class, when twenty samples at one tract carry much less than twenty
samples at twenty tracts. **The unit that repeats independently is the tract**, and both methods are
built on it.

**The bootstrap** draws the class's tracts again with replacement, as many as it has, refits each
draw, and reports each number's standard deviation over the refits. Each refit starts at the full
fit's answer.

**The sandwich** is the formula a composite likelihood takes, after Godambe:

- **A**, the curvature: how fast the objective's slope changes as each number moves;
- **B**, the spread: the sum over tracts of each tract's slope times itself, which measures how much
  the slope varies from tract to tract;
- the errors are the square roots of the diagonal of `A⁻¹ B A⁻ᵀ`.

When the objective is the true likelihood, B equals A and the formula reduces to `A⁻¹`, the plain
inverse curvature — called **naive** here. The gap between the two is the price of treating samples
at a tract as independent.

**How the slopes are computed.** A tract's slope comes from the expected counts the fit's own round
already gathers: reads that slipped and did not, slips shorter and longer, repeats slipped, allele
draws at each length. This is Fisher's identity: the slope of the log of a sum over hidden
genotypes and slips equals the expected slope of the log-likelihood that knew them, which is a short
formula in the counts (`counting_fit_errors::counting_slope`). It is exact and costs one pass over
the reads; the test `the_slope_from_the_expected_counts_is_the_objectives_slope` checks it against
the numerical slope of the objective. The curvature is taken numerically: the change in the summed
slope over a small step of each number either way, two passes a number. A purely numerical slope
was the alternative; it would need two passes a number for the slopes as well, and its error would
enter B squared.

**The concentration from pairs of samples is a second objective** (the pairwise likelihood). Its
slope is taken numerically tract by tract, and the sandwich becomes the one for a set of estimating
equations: A is not symmetric, which `A⁻¹ B A⁻ᵀ` allows. The naive errors there take the counting
objective's curvature for the stutter and shares and the pairwise objective's for the
concentration.

**Where a number gets no error.** A stutter number fitted at an end of its range (a shorter share of
exactly one when no read slipped longer) is held there and reported as not placed: the formula
needs room either side. A concentration the fit returns as absent (§7) has no error, and neither
has the between-locus share.

---

## 3. The settings drawn

Each class is drawn with a slip rate of 0.08, a shorter share of 0.8, a fall-off of 0.25 and a
concentration of 0.5 (a between-locus share of two thirds, near tomato2's 0.67–0.71). **Each setting
changes one thing from the base**: 200 tracts, 20 samples, 3 reads a sample at every tract, a
homozygote excess of 0.4, three allele lengths (one repeat short, the reference, one long) at
shares 0.26, 0.48 and 0.26, the concentration read from pairs of samples, and the counting fit run
as on real cohorts (`examples/ng_str_counting_fit.rs`: tolerance 10⁻⁷, up to 2,000 rounds, up to 10
alternations).

The **±9 window** draws nineteen lengths, −9 to +9, falling off from the reference as 0.55 to the
power of the distance, plus 0.04 at −9 and 0.03 at −8; reads are recorded to ±8, as the census does,
so an allele at ±9 has no read bucket of its own. **Alleles past the window** are drawn to −12 (with
0.03, 0.02 and 0.02 added at −12, −11 and −10) and fitted at ±9, where the −9 length stands for "9
or more short".

| setting | tracts | samples | reads a sample | homozygote excess | allele lengths | draws | replicates |
|---|---|---|---|---|---|---|---|
| base | 200 | 20 | 3 | 0.4 | ±1 | 200 | 50 |
| tracts_20 | **20** | 20 | 3 | 0.4 | ±1 | 200 | 50 |
| tracts_2000 | **2,000** | 20 | 3 | 0.4 | ±1 | 100 | 30 |
| samples_500 | 200 | **500** | 3 | 0.4 | ±1 | 40 | 20 |
| reads_30 | 200 | 20 | **30** | 0.4 | ±1 | 200 | 50 |
| excess_095 | 200 | 20 | 3 | **0.95** | ±1 | 200 | 50 |
| window_9 | 200 | 20 | 3 | 0.4 | **±9** | 60 | 20 |
| window_9_alleles_beyond | 200 | 20 | 3 | 0.4 | **±9, drawn to −12** | 40 | 20 |
| one_sample_reads_3 | 2,000 | **1** | 3 | 0.4 | ±1 | 200 | 50 |
| one_sample_reads_30 | 2,000 | **1** | **30** | 0.4 | ±1 | 200 | 50 |
| one_sample_reads_30_inbred | 2,000 | **1** | **30** | **1.0** | ±1 | 200 | 50 |

Fewer draws and replicates at the slow settings are a matter of time on a shared machine. With 200
draws a coverage of 68 has a standard error of about 3 points; with 40 draws about 7.

The tracts_2000 row was measured before the two fixes of §6 and §7. At 2,000 tracts the share prior
moves a share by at most 1/2,000 and the concentration fix does not apply (20 samples), so the row
stands; every other row was measured after them.

---

## 4. Coverage on drawn classes

Each cell is **the percentage of draws whose truth lies within one / two reported errors** of the
estimate. Honest is about 68 / 95. "(n without)" counts draws where the method gave that number no
error; the percentages are over the rest.

### 4.1 Several samples a tract

The reference-length share stands for the shares; the other shares behave alike (appendix).

**Slip rate**

| setting | spread of estimates | mean sandwich error | bootstrap | sandwich | naive |
|---|---|---|---|---|---|
| base | 0.0042 | 0.0040 | 70 / 92 | 70 / 94 | 67 / 92 |
| tracts_20 | 0.0145 | 0.0124 | 62 / 91 | 64 / 92 | 60 / 94 |
| tracts_2000 | 0.0013 | 0.0013 | 61 / 92 | 67 / 95 | 66 / 93 |
| samples_500 | 0.0010 | 0.0010 | 60 / 92 | 62 / 95 | 52 / 88 |
| reads_30 | 0.0009 | 0.0009 | 70 / 97 | 70 / 98 | 69 / 98 |
| excess_095 | 0.0027 | 0.0030 | 74 / 96 | 73 / 98 | 74 / 98 |
| window_9 | 0.0034 | 0.0034 | 67 / 92 | 68 / 95 | 63 / 92 |
| window_9_alleles_beyond | 0.0028 | 0.0034 | 78 / 95 | 78 / 98 | 70 / 95 |

**Fall-off**

| setting | spread of estimates | mean sandwich error | bootstrap | sandwich | naive |
|---|---|---|---|---|---|
| base | 0.0531 | 0.0557 | 70 / 96 | 70 / 96 | 59 / 90 |
| tracts_20 | 0.1771 | 0.1840 | 60 / 86 | 68 / 95 (13 without) | 59 / 91 (13 without) |
| tracts_2000 | 0.0183 | 0.0176 | 62 / 95 | 63 / 95 | 51 / 88 |
| samples_500 | 0.0321 | 0.0363 | 68 / 92 | 72 / 95 | 20 / 52 |
| reads_30 | 0.0086 | 0.0092 | 74 / 96 | 74 / 96 | 75 / 96 |
| excess_095 | 0.0298 | 0.0305 | 70 / 96 | 70 / 97 | 70 / 97 |
| window_9 | 0.0212 | 0.0216 | 62 / 93 | 73 / 97 | 67 / 90 |
| window_9_alleles_beyond | 0.0204 | 0.0216 | 68 / 95 | 70 / 95 | 55 / 92 |

**Reference-length share**

| setting | spread of estimates | mean sandwich error | bootstrap | sandwich | naive |
|---|---|---|---|---|---|
| base | 0.0288 | 0.0291 | 64 / 96 | 66 / 96 | 19 / 38 |
| tracts_20 | 0.0848 | 0.0841 | 68 / 94 | 66 / 94 | 23 / 39 |
| tracts_2000 | 0.0089 | 0.0092 | 67 / 97 | 66 / 98 | 26 / 43 |
| samples_500 | 0.0318 | 0.0285 | 68 / 92 | 65 / 92 | 0 / 8 |
| reads_30 | 0.0284 | 0.0285 | 68 / 95 | 68 / 94 | 23 / 41 |
| excess_095 | 0.0289 | 0.0290 | 68 / 94 | 68 / 96 | 19 / 42 |
| window_9 | 0.0283 | 0.0253 | 58 / 93 | 57 / 95 | 5 / 25 |
| window_9_alleles_beyond | 0.0242 | 0.0250 | 68 / 98 | 70 / 95 | 8 / 32 |

**Between-locus share F**

| setting | spread of estimates | mean sandwich error | bootstrap | sandwich | naive |
|---|---|---|---|---|---|
| base | 0.0231 | 0.0247 | 74 / 96 | 74 / 95 | 34 / 63 |
| tracts_20 | 0.0792 | 0.0774 | 66 / 90 | 64 / 92 | 30 / 56 |
| tracts_2000 | 0.0075 | 0.0078 | 70 / 92 | 72 / 96 | 23 / 58 |
| samples_500 | 0.0205 | 0.0219 | 78 / 95 | 75 / 98 | 2 / 10 |
| reads_30 | 0.0228 | 0.0229 | 67 / 96 | 66 / 98 | 28 / 52 |
| excess_095 | 0.0217 | 0.0246 | 72 / 96 | 74 / 97 | 37 / 66 |
| window_9 | 0.0193 | 0.0200 | 65 / 92 | 73 / 97 | 23 / 57 |
| window_9_alleles_beyond | 0.0242 | 0.0201 | 62 / 90 | 62 / 90 | 28 / 45 |

**In every row the bootstrap and the sandwich sit near 68 / 95, within the noise of the row's
draws, and near each other**: their mean errors differ by at most 18% (the concentration at 20
tracts) and by under 10% in most rows. **The naive column sits at 5 to 37 within one error for the
share and the between-locus share at 20 samples, and at 0 to 2 at 500.** The naive gap is the size of
the sharing between samples: at a between-locus share of two thirds, the twenty samples of a tract
are worth about one and a half independent samples for the reference share (its naive error is 0.0079
against a spread of 0.0288 at the base setting, a variance thirteen times too small).

**At 20 tracts both methods still cover** (60–77 within one error, 86–95 within two), with the
bootstrap low on the fall-off (60 / 86). The fall-off is barely measured there: its spread over draws
is 0.18 around a truth of 0.25, and in 13 of 200 draws it was fitted at an end of its range, where
the sandwich gives it no error. In 6 of 200 draws no read slipped longer, the shorter share was
fitted at one, and neither method gives it an error.

### 4.2 One sample

| setting | number | spread of estimates | bootstrap | sandwich | naive |
|---|---|---|---|---|---|
| one_sample_reads_3 | level | 0.0118 | 64 / 92 | 64 / 92 | 64 / 92 |
| one_sample_reads_3 | shorter_share | 0.0503 | 68 / 92 | 68 / 96 | 68 / 96 |
| one_sample_reads_3 | fall_off | 0.0664 | 70 / 94 | 70 / 96 | 68 / 96 |
| one_sample_reads_3 | share_0 | 0.0120 | 64 / 94 | 65 / 96 | 65 / 96 |
| one_sample_reads_3 | between_locus_share | 0.0842 | 66 / 92 | 66 / 94 | 66 / 94 |
| one_sample_reads_30 | level | 0.0012 | 76 / 94 | 76 / 96 | 78 / 96 |
| one_sample_reads_30 | shorter_share | 0.0074 | 63 / 96 | 64 / 94 | 64 / 94 |
| one_sample_reads_30 | fall_off | 0.0137 | 63 / 94 | 62 / 94 | 64 / 94 |
| one_sample_reads_30 | share_0 | 0.0111 | 68 / 94 | 66 / 94 | 66 / 94 |
| one_sample_reads_30 | between_locus_share | 0.0207 | 64 / 94 | 64 / 93 | 64 / 93 |
| one_sample_reads_30_inbred | level | 0.0013 | 70 / 94 | 70 / 96 | 70 / 96 |
| one_sample_reads_30_inbred | shorter_share | 0.0070 | 64 / 96 | 66 / 95 | 65 / 94 |
| one_sample_reads_30_inbred | fall_off | 0.0128 | 64 / 94 | 64 / 94 | 64 / 93 |
| one_sample_reads_30_inbred | share_0 | 0.0115 | 63 / 93 | 64 / 94 | 64 / 94 |
| one_sample_reads_30_inbred | between_locus_share | absent | — (none of 200) | — (none of 200) | — (none of 200) |

These run the pairs entry point, `fit_class_by_counting_with_pairs`, which at one sample now reads
the concentration from the sample's own two copies (§7); before that change the same draws gave
the slip rates of §7's table.

**At one sample the composite likelihood is the exact likelihood** — there is no second sample at a
tract to share alleles with — so the naive errors match the sandwich's, as they do in every row
here. That is a check of the sandwich's code as much as a result. At an excess of one the sample's
two copies are always one draw and say nothing about how tracts differ: the concentration is absent
in all 200 draws, and so are its errors (§7).

---

## 5. Coverage on the fitting bench

The bench of phase 0 draws a whole motif length (period 1: 14 classes, 8 to 21 repeats) from the
tomato-like truth — slip rates from 0.005 to 0.27, homozygote excesses of 0.85–0.99 a sample, a
between-locus share of 0.7, alleles drawn to ±14 — with at most 100 loci a class. Each cell below
pools the estimates of every class and seed. Shares are grouped by their true value.

| cell | number | estimates | sandwich: within 1 / 2 errors |
|---|---|---|---|
| 60 samples, 3 reads, seeds 1–6 | slip rate, 8–15 / 16–21 repeats | 48 / 36 | 67 / 98 and 72 / 89 |
| | shorter share, 8–15 / 16–21 | 48 / 36 | 67 / 88 and 72 / 94 |
| | fall-off, 8–15 / 16–21 | 48 / 36 | 65 / 94 and 72 / 100 |
| | between-locus share | 84 | 63 / 93 |
| | shares, truth ≥ 0.05 | 294 | 64 / 87 |
| | shares, truth 0.01–0.05 | 402 | 50 / 67 |
| | shares, truth < 0.01 | 900 | 36 / 50 |
| 60 samples, 30 reads, seeds 1–6 | slip rate, 8–15 / 16–21 | 48 / 36 | 77 / 96 and 72 / 86 |
| | shorter share, 8–15 / 16–21 | 48 / 36 | 62 / 90 and 67 / 86 |
| | fall-off, 8–15 / 16–21 | 48 / 36 | 60 / 92 and 69 / 86 |
| | between-locus share | 84 | 73 / 94 |
| | shares, truth ≥ 0.05 | 294 | 66 / 90 |
| | shares, truth 0.01–0.05 | 402 | 51 / 68 |
| | shares, truth < 0.01 | 900 | 31 / 43 |
| 2,000 samples, 3 reads, seed 1 | slip rate, 8–15 / 16–21 | 8 / 6 | 50 / 88 and **17 / 33** |
| | shorter share, 8–15 / 16–21 | 8 / 6 | 38 / 88 and 67 / 67 |
| | fall-off, 8–15 / 16–21 | 8 / 6 | **25 / 62** and 67 / 83 |
| | between-locus share | 14 | 57 / 100 |
| | shares, truth ≥ 0.05 | 49 | 61 / 90 |
| | shares, truth 0.01–0.05 | 67 | 54 / 67 |
| | shares, truth < 0.01 | 150 | 33 / 47 |

### 5.1 Shares few tracts carry

**A share near zero at about 70 loci has a sampling distribution nothing like a bell.** At a
between-locus share of 0.7 most loci are fixed at one length, so a share of 0.005 is carried by a
fraction of one locus on average; the estimate is either near the prior's floor or a few times the
truth, and a symmetric ± error around either misses. Both methods share this, as both describe the
spread by one number.

On two seeds at 60 samples and three reads the bootstrap covers shares below 0.01 at 32 / 45 and
the sandwich at 37 / 52; on the same estimates both cover shares above 0.05 at 65 / 88, and every
stutter number and the between-locus share alike (on 12 to 28 estimates a number, the bootstrap
within two estimates of the sandwich).

On the drawn ±9 window at 200 tracts (appendix) the same holds: shares of 0.008 or less cover at
50–73 within one error and 60–83 within two, shares of 0.014 or more at 57–75 and 85–98. That puts the line near **3 tracts' worth of the length**
(share × tracts), measured on two settings only.

### 5.2 Bias at long tracts and 2,000 samples

**At 2,000 samples the stutter's errors are small enough that a small bias shows.** On that cell's
one seed the slip rate of the six classes of 16 to 21 repeats came out at 0.2636 to 0.2686 against a
truth of 0.27, 0.5 to 2.4% low and 0.9 to 3.3 errors from it; the fall-off at 8 to 10 repeats came
out at 0.080 to 0.135 against 0.10 to 0.14, 2.0 to 3.3 errors low. At 60 samples the same classes
cover (the rows above).

**Alleles past the fit's outermost length were the first suspect**, because the bench draws alleles
to ±14 and the fit stops at ±9, which then stands for "9 or more". A drawn setting tests it alone:
alleles drawn to −12, fitted at ±9 (window_9_alleles_beyond, 40 draws). There the slip rate comes
out at 0.0803 against 0.08, with a spread of 0.0028, and the stutter numbers cover within the noise of 40 draws (slip rate 78 / 98, shorter share 62 / 90, fall-off 70 / 95, sandwich).
**So the far alleles do not bias the stutter by as much as 0.003 at a slip rate of 0.08.** They do
move share mass: the −8 share comes out at 0.0107 against a truth of 0.0023, covered in 0 draws of
40, because an allele 10 to 12 repeats short sends its one-step slips to where an allele at −8
sends its exact reads.

**What causes the bias at 2,000 samples is not settled.** It is one seed of 14 classes, and the
candidates left are the bench's own model — slip rates near 0.27, per-sample homozygote excesses
drawn from 0.85–0.99, Poisson depths — none of which the drawn settings vary. A run of the
2,000-sample cell over several seeds, with the excess fixed and then with the slip rate lowered,
would split them; each run costs about an hour on this machine as it is loaded.

---

## 6. Shares of exactly zero (phase 0's first finding)

**The cause, locus by locus.** The bench's held-out score fits a method on four fifths of each
class's loci and scores the other fifth. A new option (`--locus-out`) writes each held-out locus's
score with the read length most of its reads report and the share the fit gave that length. Before
the fix, on period 1, seed 1:

| cell | held-out loci | at a length the fit gave a share of exactly 0 | their loss against the truth | the rest's loss |
|---|---|---|---|---|
| 1 sample, 300 reads | 971 | 64 | 52,263 units (817 a locus) | 21 units (0.02 a locus) |
| 2,000 samples, 3 reads | 276 | 10 | 94,608 units (9,461 a locus) | 11,701 units (44 a locus) |

At 2,000 samples the other 266 loci still lose 44 units each; a locus with a minority allele at a
zero-share length is not caught by looking at its main length alone, which is the likely remainder.
Before the fix the fit returned 74 shares of exactly zero over the 14 classes of that cell.

**Why a share reaches exactly zero and stays there.** The counting fit updates a share as the
length's expected allele draws over all draws. A length that no training locus carries gets a
vanishing share, the genotypes using it fall below the fit's posterior floor, it draws nothing, and
its share becomes zero — after which its prior is zero, so it can never draw again. A share of zero
is also what maximum likelihood says for a length no training locus carries, so the fix is a prior,
not a change to the rounds.

**The fix: a prior weight of one tract's worth of allele draws, spread evenly over the lengths**
(`composite_fit::SHARE_PRIOR_TRACTS`). Each round adds `draws ÷ tracts ÷ lengths` to every length's
draws. Counting it in tracts rather than in draws keeps its size against the evidence fixed whatever
the cohort: at 2,000 samples a tract holds about 4,000 draws and a prior of a few draws would be
lost. Its pull on a share is at most about `1 / tracts`:

| setting | reference-length share, before → after | spread of estimates |
|---|---|---|
| base, 200 tracts | 0.4800 → 0.4791 (truth 0.4762) | 0.029 |
| tracts_20 | 0.4773 → 0.4690 (truth 0.4762) | 0.085 |
| one sample, 30 reads, 2,000 tracts | 0.4774 → 0.4773 | 0.011 |

The largest shift, at 20 tracts, is a tenth of the spread of the estimates, and coverage did not
move (§4).

**The held-out loss after the fix** (same cells, same draws):

| cell | held-out score against the truth, before → after, a locus |
|---|---|
| 1 sample, 300 reads (575 loci scored both times) | −47.8 → −0.14; the 33 of them at a zero share −833 → −2.4 |
| 2,000 samples, 3 reads (276 loci) | **−385 → +0.56** |

At 2,000 samples every one of the 14 classes now sits within 2.7 units a locus of the truth
(before: 7 classes lost 81 to 1,836). The fit beating the truth by 0.56 a locus is within what one
seed's held-out loci allow; the Monte Carlo error of the figure is 0.0025.

---

## 7. The concentration at one sample (phase 0's second finding)

**What the fit did.** With the concentration read from pairs of samples, a class with no tract of
two samples has no pair: the pairwise likelihood is flat, its search returned its upper end (10⁴, a
between-locus share of 10⁻⁴), and the stutter and shares were then counted at that concentration.
That says a sample is heterozygous far more often than it is, so reads that differ are put down to
two alleles rather than to slips. On drawn classes, before the fix:

| one sample, 2,000 tracts | slip rate (truth 0.08) | shorter share (truth 0.8) | reference share (truth 0.476) |
|---|---|---|---|
| 3 reads | 0.024, spread 0.004 | 1.000 in every draw | 0.466 |
| 30 reads | 0.080 | 0.799 | 0.496, spread 0.013; 27 draws in 100 within one sandwich error |

**The rule now.** `fit_class_by_counting_with_pairs` reads the concentration from pairs only where
some tract has two samples; otherwise it keeps the estimate from each sample's own two copies,
which is what a single sample has. **A concentration is returned as absent (`None`) when the data
cannot place it** (`composite_fit::concentration_is_placed`):

- its estimator has nothing to read: no tract of two samples for the pairs; no sample with a
  homozygote excess below one for the own copies, since a sample whose copies are always one draw
  says nothing about how tracts differ; or
- the estimate stopped within 1% of an end of its search, 10⁻³ or 10⁴.

The fit still counts the stutter at some concentration (`concentration_counted_at`), which the
errors and the refits use; it is not reported. **Both error methods give the concentration and the
between-locus share no error when the fit gives no concentration.**

After the fix the one-sample rows of §4.2 recover every number, with honest errors, and at an excess
of one the concentration is absent in all 200 draws. On the fitting bench at one sample (excess
0.85–0.99) the concentration fitted on all loci is placed in 9 of 14 classes at 300 reads and in none
at 3 reads. The bench's held-out score needs a concentration, so a class without one is not scored: at 3 reads that
leaves 18 of 911 held-out loci scored. Whether the score should integrate over the concentration in
that case is a question for the bench.

---

## 8. Cost

**Measured one draw after another, three draws a setting, with every replicate of a bootstrap run in
parallel on the shared machine.** Seconds a class.

| setting | one fit | sandwich (× a fit) | bootstrap, 50 replicates (× a fit) |
|---|---|---|---|
| base | 11.2 | 0.30 (0.03×) | 20.0 (1.8×) |
| tracts_20 | 9.8 | 0.22 (0.02×) | 9.1 (0.9×) |
| reads_30 | 3.5 | 0.63 (0.18×) | 7.3 (2.1×) |
| excess_095 | 34.6 | 0.43 (0.01×) | 8.6 (0.25×) |
| one_sample_reads_3 | 7.3 | 0.24 (0.03×) | 22.2 (3.0×) |
| one_sample_reads_30 | 1.2 | 0.36 (0.30×) | 1.9 (1.6×) |
| window_9 (19 lengths) | 45.7 | 3.9 (0.09×) | 408 (8.9×) |

On the fitting bench's period 1 at 60 samples and 3 reads (all 14 classes, 19 lengths, one run
each): the fit alone 39–46 s, the fit with the sandwich 53–63 s, the fit with a 50-replicate
bootstrap 512–547 s.

**Why the bootstrap is cheaper than 50 fits.** A replicate starts at the full fit's answer, and with
the concentration from pairs a fit from scratch spends most of its time in the first counting from
flat shares and in the spread's searches, which a warm start mostly skips. At nineteen lengths a
replicate's counting takes more rounds and the advantage shrinks (8.9× for 50 replicates).

**Why the sandwich is cheap.** It is two passes over the reads for each number fitted, plus one:
8 numbers at three lengths, 24 at nineteen, against a fit that takes tens to hundreds of rounds.

---

## 9. Recommendation

1. **Report the sandwich as the counting fit's standard errors.** It is honest wherever the
   bootstrap is, its errors are within 18% of the bootstrap's in every drawn setting (under 10% in
   most), and it costs 0.01 to 0.3 times what the fit costs.
2. **Keep the bootstrap as a check**, run on a sample of classes when the fit changes, not on every
   class of every run: it costs one refit a replicate.
3. **Do not report the plain curvature** for the shares or the concentration.
4. **Flag the errors of shares few tracts carry.** Below about 0.01 at 70 to 200 tracts neither
   method's error means what it says (§5.1). A rule such as "no error where share × tracts is below
   3" would say so; the threshold rests on two settings and is not built.
5. **Keep both fixes** (§6, §7). Outside the cases they fix, they move the answers by less than
   the noise: coverage on every drawn setting is unchanged, and the reference share moves by a
   tenth of its spread at 20 tracts and less elsewhere.
6. **For the spec** (`str_counting_fit.md` §3, open question 1): the sandwich, the slope from the
   expected counts, the curvature from two passes a number. Open question 2 — whether the
   parameters file carries the shares' and the concentration's errors — is not answered here; both
   exist if it should.

---

## 10. What is not settled

- **Only one truth on drawn classes**: slip rate 0.08, concentration 0.5. The fitting bench adds the
  tomato-like curves, at 60 and 2,000 samples.
- **Few draws at the slow settings** (§3), and one seed at the bench's 2,000-sample cell.
- **The stutter's bias at long tracts and 2,000 samples** (§5.2): seen on one seed, cause open.
- **The bench's held-out score skips classes without a concentration** (§7).
- **Times are provisional** (§8).

---

## 11. How to rerun

```text
# coverage on drawn classes, one setting
cargo run --release --features bench-fixtures --example ng_str_counting_fit_errors -- \
  --draws 200 --replicates 50 --output tmp/coverage --setting base
# cost
cargo run --release --features bench-fixtures --example ng_str_counting_fit_errors -- \
  --cost --draws 2 --replicates 50 --output tmp/cost
# the bench: held-out with the per-locus file, and error honesty without held-out
target/release/examples/ng_str_fit_bench --period 1 --samples 2000 --depth 3 --seed 1 --truth pcr \
  --methods B --loci-at-most 100 --heldout-loci 20 --out rows.tsv --locus-out loci.tsv
target/release/examples/ng_str_fit_bench --period 1 --samples 60 --depth 3 --seed 1 --truth pcr \
  --methods B,B_bootstrap --loci-at-most 100 --folds 0 --out rows.tsv
```

Code: `src/parameter_estimation/joint/ssr_fit/counting_fit_errors.rs` (`sandwich_errors`,
`bootstrap_errors`) and `composite_fit.rs` (the fixes); the bench's method B reports the sandwich
errors, `B_bootstrap` the bootstrap's.

---

## Appendix: full coverage tables

### Settings with several samples a tract
#### tracts_20 — 20 tracts, 20 samples, 3 reads a sample, homozygote excess 0.4, window ±1, concentration from pairs of samples

Rounds of the last counting: median 85, most 1009. 200 draws, 50 bootstrap replicates each, 249 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0787 | 0.0145 | 62 / 91 | 64 / 92 | 60 / 94 | 0.0124 / 0.0124 / 0.0120 |
| shorter_share | 0.8000 | 0.7962 | 0.0978 | 67 / 91 (6 without) | 68 / 95 (6 without) | 66 / 96 (6 without) | 0.0942 / 0.0977 / 0.0908 |
| fall_off | 0.2500 | 0.2475 | 0.1771 | 60 / 86 | 68 / 95 (13 without) | 59 / 91 (13 without) | 0.1691 / 0.1840 / 0.1476 |
| share_-1 | 0.2619 | 0.2689 | 0.0779 | 62 / 92 | 66 / 92 | 24 / 40 | 0.0727 / 0.0741 / 0.0216 |
| share_0 | 0.4762 | 0.4690 | 0.0848 | 68 / 94 | 66 / 94 | 23 / 39 | 0.0840 / 0.0841 / 0.0240 |
| share_1 | 0.2619 | 0.2622 | 0.0789 | 66 / 88 | 67 / 90 | 28 / 49 | 0.0719 / 0.0725 / 0.0212 |
| concentration | 0.5000 | 0.5826 | 0.2052 | 77 / 94 | 74 / 94 | 30 / 60 | 0.2357 / 0.1990 / 0.0800 |
| between_locus_share | 0.6667 | 0.6419 | 0.0792 | 66 / 90 | 64 / 92 | 30 / 56 | 0.0791 / 0.0774 / 0.0310 |

#### base — 200 tracts, 20 samples, 3 reads a sample, homozygote excess 0.4, window ±1, concentration from pairs of samples

Rounds of the last counting: median 70, most 140. 200 draws, 50 bootstrap replicates each, 823 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0801 | 0.0042 | 70 / 92 | 70 / 94 | 67 / 92 | 0.0039 / 0.0040 / 0.0037 |
| shorter_share | 0.8000 | 0.7997 | 0.0285 | 68 / 93 | 70 / 94 | 66 / 94 | 0.0286 / 0.0284 / 0.0270 |
| fall_off | 0.2500 | 0.2487 | 0.0531 | 70 / 96 | 70 / 96 | 59 / 90 | 0.0561 / 0.0557 / 0.0441 |
| share_-1 | 0.2619 | 0.2617 | 0.0274 | 64 / 94 | 66 / 95 | 18 / 36 | 0.0254 / 0.0257 / 0.0071 |
| share_0 | 0.4762 | 0.4791 | 0.0288 | 64 / 96 | 66 / 96 | 19 / 38 | 0.0289 / 0.0291 / 0.0079 |
| share_1 | 0.2619 | 0.2592 | 0.0243 | 68 / 95 | 68 / 94 | 23 / 46 | 0.0251 / 0.0254 / 0.0070 |
| concentration | 0.5000 | 0.5004 | 0.0522 | 73 / 96 | 74 / 96 | 33 / 62 | 0.0558 / 0.0557 / 0.0217 |
| between_locus_share | 0.6667 | 0.6673 | 0.0231 | 74 / 96 | 74 / 95 | 34 / 63 | 0.0244 / 0.0247 / 0.0096 |

#### reads_30 — 200 tracts, 20 samples, 30 reads a sample, homozygote excess 0.4, window ±1, concentration from pairs of samples

Rounds of the last counting: median 2, most 3. 200 draws, 50 bootstrap replicates each, 174 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0800 | 0.0009 | 70 / 97 | 70 / 98 | 69 / 98 | 0.0009 / 0.0009 / 0.0010 |
| shorter_share | 0.8000 | 0.7999 | 0.0046 | 78 / 96 | 77 / 97 | 78 / 98 | 0.0049 / 0.0049 / 0.0049 |
| fall_off | 0.2500 | 0.2500 | 0.0086 | 74 / 96 | 74 / 96 | 75 / 96 | 0.0092 / 0.0092 / 0.0092 |
| share_-1 | 0.2619 | 0.2630 | 0.0244 | 70 / 96 | 70 / 95 | 21 / 40 | 0.0250 / 0.0249 / 0.0065 |
| share_0 | 0.4762 | 0.4744 | 0.0284 | 68 / 95 | 68 / 94 | 23 / 41 | 0.0284 / 0.0285 / 0.0074 |
| share_1 | 0.2619 | 0.2626 | 0.0229 | 72 / 96 | 74 / 96 | 25 / 43 | 0.0250 / 0.0248 / 0.0065 |
| concentration | 0.5000 | 0.5059 | 0.0515 | 68 / 96 | 66 / 94 | 28 / 52 | 0.0521 / 0.0521 / 0.0176 |
| between_locus_share | 0.6667 | 0.6648 | 0.0228 | 67 / 96 | 66 / 98 | 28 / 52 | 0.0227 / 0.0229 / 0.0077 |

#### excess_095 — 200 tracts, 20 samples, 3 reads a sample, homozygote excess 0.95, window ±1, concentration from pairs of samples

Rounds of the last counting: median 22, most 31. 200 draws, 50 bootstrap replicates each, 376 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0798 | 0.0027 | 74 / 96 | 73 / 98 | 74 / 98 | 0.0030 / 0.0030 / 0.0030 |
| shorter_share | 0.8000 | 0.7977 | 0.0192 | 72 / 95 | 70 / 95 | 72 / 94 | 0.0194 / 0.0193 / 0.0192 |
| fall_off | 0.2500 | 0.2526 | 0.0298 | 70 / 96 | 70 / 97 | 70 / 97 | 0.0305 / 0.0305 / 0.0306 |
| share_-1 | 0.2619 | 0.2623 | 0.0243 | 68 / 95 | 70 / 96 | 22 / 40 | 0.0253 / 0.0255 / 0.0070 |
| share_0 | 0.4762 | 0.4803 | 0.0289 | 68 / 94 | 68 / 96 | 19 / 42 | 0.0287 / 0.0290 / 0.0080 |
| share_1 | 0.2619 | 0.2574 | 0.0241 | 68 / 98 | 68 / 98 | 22 / 38 | 0.0252 / 0.0253 / 0.0070 |
| concentration | 0.5000 | 0.5059 | 0.0492 | 72 / 96 | 72 / 98 | 37 / 67 | 0.0566 / 0.0559 / 0.0237 |
| between_locus_share | 0.6667 | 0.6647 | 0.0217 | 72 / 96 | 74 / 97 | 37 / 66 | 0.0247 / 0.0246 / 0.0104 |

#### window_9 — 200 tracts, 20 samples, 3 reads a sample, homozygote excess 0.4, window ±9, concentration from pairs of samples

Rounds of the last counting: median 128, most 258. 60 draws, 20 bootstrap replicates each, 5344 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0802 | 0.0034 | 67 / 92 | 68 / 95 | 63 / 92 | 0.0033 / 0.0034 / 0.0031 |
| shorter_share | 0.8000 | 0.7969 | 0.0232 | 60 / 92 | 62 / 97 | 58 / 97 | 0.0224 / 0.0223 / 0.0213 |
| fall_off | 0.2500 | 0.2534 | 0.0212 | 62 / 93 | 73 / 97 | 67 / 90 | 0.0207 / 0.0216 / 0.0189 |
| share_-9 | 0.0386 | 0.0404 | 0.0122 | 72 / 95 | 75 / 95 | 83 / 100 | 0.0131 / 0.0141 / 0.0160 |
| share_-8 | 0.0303 | 0.0343 | 0.0134 | 62 / 87 | 67 / 87 | 78 / 97 | 0.0114 / 0.0124 / 0.0160 |
| share_-7 | 0.0041 | 0.0048 | 0.0030 | 68 / 80 | 73 / 80 | 23 / 47 | 0.0032 / 0.0033 / 0.0010 |
| share_-6 | 0.0075 | 0.0079 | 0.0041 | 72 / 83 | 73 / 83 | 22 / 50 | 0.0047 / 0.0045 / 0.0013 |
| share_-5 | 0.0137 | 0.0156 | 0.0079 | 58 / 85 | 62 / 87 | 22 / 38 | 0.0065 / 0.0067 / 0.0018 |
| share_-4 | 0.0249 | 0.0235 | 0.0089 | 57 / 88 | 58 / 88 | 12 / 38 | 0.0079 / 0.0082 / 0.0023 |
| share_-3 | 0.0453 | 0.0445 | 0.0130 | 60 / 92 | 62 / 92 | 22 / 35 | 0.0111 / 0.0115 / 0.0032 |
| share_-2 | 0.0823 | 0.0804 | 0.0161 | 63 / 90 | 63 / 93 | 18 / 40 | 0.0151 / 0.0154 / 0.0042 |
| share_-1 | 0.1497 | 0.1500 | 0.0210 | 65 / 90 | 68 / 95 | 17 / 40 | 0.0192 / 0.0202 / 0.0055 |
| share_0 | 0.2722 | 0.2708 | 0.0283 | 58 / 93 | 57 / 95 | 5 / 25 | 0.0253 / 0.0253 / 0.0068 |
| share_1 | 0.1497 | 0.1460 | 0.0232 | 65 / 85 | 65 / 85 | 18 / 43 | 0.0192 / 0.0200 / 0.0054 |
| share_2 | 0.0823 | 0.0823 | 0.0144 | 63 / 93 | 68 / 98 | 18 / 35 | 0.0151 / 0.0155 / 0.0042 |
| share_3 | 0.0453 | 0.0460 | 0.0115 | 60 / 90 | 67 / 95 | 23 / 38 | 0.0114 / 0.0117 / 0.0032 |
| share_4 | 0.0249 | 0.0253 | 0.0086 | 65 / 93 | 67 / 95 | 20 / 37 | 0.0083 / 0.0085 / 0.0024 |
| share_5 | 0.0137 | 0.0141 | 0.0057 | 65 / 93 | 70 / 93 | 25 / 40 | 0.0059 / 0.0063 / 0.0017 |
| share_6 | 0.0075 | 0.0062 | 0.0036 | 62 / 75 | 67 / 73 | 32 / 42 | 0.0036 / 0.0037 / 0.0011 |
| share_7 | 0.0041 | 0.0043 | 0.0038 | 55 / 63 | 53 / 63 | 20 / 40 | 0.0027 / 0.0028 / 0.0009 |
| share_8 | 0.0023 | 0.0019 | 0.0016 | 43 / 60 | 52 / 60 | 50 / 58 | 0.0011 / 0.0013 / 0.0011 |
| share_9 | 0.0013 | 0.0019 | 0.0019 | 53 / 68 | 50 / 70 | 48 / 73 | 0.0013 / 0.0012 / 0.0011 |
| concentration | 0.5000 | 0.5148 | 0.0447 | 65 / 95 | 72 / 98 | 23 / 55 | 0.0443 / 0.0460 / 0.0170 |
| between_locus_share | 0.6667 | 0.6607 | 0.0193 | 65 / 92 | 73 / 97 | 23 / 57 | 0.0193 / 0.0200 / 0.0074 |

#### samples_500 — 200 tracts, 500 samples, 3 reads a sample, homozygote excess 0.4, window ±1, concentration from pairs of samples

Rounds of the last counting: median 71, most 98. 40 draws, 20 bootstrap replicates each, 2075 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0800 | 0.0010 | 60 / 92 | 62 / 95 | 52 / 88 | 0.0009 / 0.0010 / 0.0007 |
| shorter_share | 0.8000 | 0.8008 | 0.0093 | 50 / 95 | 55 / 98 | 38 / 68 | 0.0083 / 0.0091 / 0.0054 |
| fall_off | 0.2500 | 0.2478 | 0.0321 | 68 / 92 | 72 / 95 | 20 / 52 | 0.0347 / 0.0363 / 0.0088 |
| share_-1 | 0.2619 | 0.2566 | 0.0224 | 62 / 98 | 62 / 98 | 12 / 20 | 0.0249 / 0.0252 / 0.0014 |
| share_0 | 0.4762 | 0.4813 | 0.0318 | 68 / 92 | 65 / 92 | 0 / 8 | 0.0281 / 0.0285 / 0.0016 |
| share_1 | 0.2619 | 0.2621 | 0.0275 | 60 / 92 | 75 / 88 | 2 / 5 | 0.0250 / 0.0250 / 0.0014 |
| concentration | 0.5000 | 0.5114 | 0.0478 | 75 / 95 | 78 / 98 | 2 / 10 | 0.0494 / 0.0501 / 0.0044 |
| between_locus_share | 0.6667 | 0.6623 | 0.0205 | 78 / 95 | 75 / 98 | 2 / 10 | 0.0214 / 0.0219 / 0.0019 |

#### window_9_alleles_beyond — 200 tracts, 20 samples, 3 reads a sample, homozygote excess 0.4, window ±9, concentration from pairs of samples

Rounds of the last counting: median 66, most 133. 40 draws, 20 bootstrap replicates each, 3716 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0803 | 0.0028 | 78 / 95 | 78 / 98 | 70 / 95 | 0.0034 / 0.0034 / 0.0031 |
| shorter_share | 0.8000 | 0.8001 | 0.0264 | 68 / 85 | 62 / 90 | 58 / 85 | 0.0217 / 0.0227 / 0.0213 |
| fall_off | 0.2500 | 0.2513 | 0.0204 | 68 / 95 | 70 / 95 | 55 / 92 | 0.0218 / 0.0216 / 0.0187 |
| share_-9 | 0.0679 | 0.0640 | 0.0146 | 78 / 85 | 80 / 95 | 62 / 80 | 0.0141 / 0.0148 / 0.0096 |
| share_-8 | 0.0023 | 0.0107 | 0.0021 | 0 / 8 | 0 / 5 | 90 / 100 | 0.0029 / 0.0027 / 0.0090 |
| share_-7 | 0.0041 | 0.0042 | 0.0033 | 60 / 72 | 52 / 70 | 22 / 38 | 0.0031 / 0.0029 / 0.0009 |
| share_-6 | 0.0075 | 0.0077 | 0.0050 | 68 / 80 | 68 / 85 | 20 / 38 | 0.0042 / 0.0043 / 0.0013 |
| share_-5 | 0.0137 | 0.0134 | 0.0059 | 68 / 90 | 70 / 92 | 30 / 58 | 0.0062 / 0.0062 / 0.0017 |
| share_-4 | 0.0248 | 0.0227 | 0.0099 | 52 / 80 | 55 / 82 | 10 / 22 | 0.0078 / 0.0081 / 0.0022 |
| share_-3 | 0.0452 | 0.0446 | 0.0128 | 52 / 95 | 50 / 95 | 8 / 30 | 0.0113 / 0.0117 / 0.0032 |
| share_-2 | 0.0821 | 0.0831 | 0.0175 | 70 / 92 | 72 / 95 | 25 / 42 | 0.0161 / 0.0157 / 0.0043 |
| share_-1 | 0.1493 | 0.1496 | 0.0197 | 75 / 90 | 75 / 92 | 32 / 48 | 0.0194 / 0.0203 / 0.0055 |
| share_0 | 0.2715 | 0.2665 | 0.0242 | 68 / 98 | 70 / 95 | 8 / 32 | 0.0249 / 0.0250 / 0.0068 |
| share_1 | 0.1493 | 0.1461 | 0.0207 | 68 / 90 | 62 / 92 | 18 / 50 | 0.0198 / 0.0200 / 0.0054 |
| share_2 | 0.0821 | 0.0869 | 0.0143 | 78 / 98 | 72 / 100 | 20 / 48 | 0.0165 / 0.0160 / 0.0043 |
| share_3 | 0.0452 | 0.0468 | 0.0116 | 72 / 98 | 75 / 98 | 25 / 42 | 0.0117 / 0.0120 / 0.0032 |
| share_4 | 0.0248 | 0.0247 | 0.0091 | 68 / 88 | 65 / 88 | 18 / 32 | 0.0084 / 0.0084 / 0.0023 |
| share_5 | 0.0137 | 0.0124 | 0.0061 | 55 / 85 | 62 / 88 | 8 / 22 | 0.0059 / 0.0058 / 0.0016 |
| share_6 | 0.0075 | 0.0070 | 0.0046 | 60 / 78 | 62 / 78 | 18 / 50 | 0.0039 / 0.0040 / 0.0012 |
| share_7 | 0.0041 | 0.0046 | 0.0036 | 50 / 70 | 55 / 70 | 32 / 42 | 0.0030 / 0.0031 / 0.0009 |
| share_8 | 0.0023 | 0.0021 | 0.0016 | 50 / 70 | 52 / 68 | 50 / 75 | 0.0013 / 0.0013 / 0.0013 |
| share_9 | 0.0025 | 0.0029 | 0.0026 | 55 / 70 | 55 / 70 | 50 / 65 | 0.0020 / 0.0021 / 0.0013 |
| concentration | 0.5000 | 0.5045 | 0.0561 | 62 / 92 | 65 / 88 | 28 / 45 | 0.0462 / 0.0457 / 0.0167 |
| between_locus_share | 0.6667 | 0.6656 | 0.0242 | 62 / 90 | 62 / 90 | 28 / 45 | 0.0200 / 0.0201 / 0.0074 |

(tracts_2000 measured before the fixes of §6 and §7; see §3.)
#### tracts_2000 — 2000 tracts, 20 samples, 3 reads a sample, homozygote excess 0.4, window ±1, concentration from pairs of samples

Rounds of the last counting: median 69, most 88. 100 draws, 30 bootstrap replicates each, 3177 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0799 | 0.0013 | 61 / 92 | 67 / 95 | 66 / 93 | 0.0012 / 0.0013 / 0.0012 |
| shorter_share | 0.8000 | 0.8003 | 0.0081 | 75 / 97 | 74 / 98 | 74 / 97 | 0.0092 / 0.0091 / 0.0086 |
| fall_off | 0.2500 | 0.2473 | 0.0183 | 62 / 95 | 63 / 95 | 51 / 88 | 0.0171 / 0.0176 / 0.0139 |
| share_-1 | 0.2619 | 0.2615 | 0.0078 | 67 / 96 | 69 / 98 | 19 / 42 | 0.0082 / 0.0082 / 0.0023 |
| share_0 | 0.4762 | 0.4769 | 0.0089 | 67 / 97 | 66 / 98 | 26 / 43 | 0.0091 / 0.0092 / 0.0025 |
| share_1 | 0.2619 | 0.2616 | 0.0090 | 64 / 91 | 64 / 93 | 16 / 34 | 0.0079 / 0.0081 / 0.0022 |
| concentration | 0.5000 | 0.5007 | 0.0169 | 70 / 92 | 71 / 97 | 23 / 58 | 0.0169 / 0.0176 / 0.0068 |
| between_locus_share | 0.6667 | 0.6664 | 0.0075 | 70 / 92 | 72 / 96 | 23 / 58 | 0.0075 / 0.0078 / 0.0030 |

### One sample
#### one_sample_reads_30_inbred — 2000 tracts, 1 samples, 30 reads a sample, homozygote excess 1, window ±1, concentration from pairs of samples

Rounds of the last counting: median 54, most 62. 200 draws, 50 bootstrap replicates each, 86 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0799 | 0.0013 | 70 / 94 | 70 / 96 | 70 / 96 | 0.0013 / 0.0013 / 0.0013 |
| shorter_share | 0.8000 | 0.8002 | 0.0070 | 64 / 96 | 66 / 95 | 65 / 94 | 0.0066 / 0.0067 / 0.0067 |
| fall_off | 0.2500 | 0.2513 | 0.0128 | 64 / 94 | 64 / 94 | 64 / 93 | 0.0122 / 0.0123 / 0.0122 |
| share_-1 | 0.2619 | 0.2624 | 0.0097 | 72 / 94 | 73 / 94 | 73 / 94 | 0.0099 / 0.0098 / 0.0098 |
| share_0 | 0.4762 | 0.4753 | 0.0115 | 63 / 93 | 64 / 94 | 64 / 94 | 0.0111 / 0.0112 / 0.0112 |
| share_1 | 0.2619 | 0.2623 | 0.0098 | 70 / 94 | 70 / 96 | 70 / 96 | 0.0099 / 0.0098 / 0.0098 |
| concentration | 0.5000 | NaN | NaN | — (none of 200) | — (none of 200) | — (none of 200) | 0.0000 / 0.0000 / 0.0000 |
| between_locus_share | 0.6667 | NaN | NaN | — (none of 200) | — (none of 200) | — (none of 200) | 0.0000 / 0.0000 / 0.0000 |

#### one_sample_reads_30 — 2000 tracts, 1 samples, 30 reads a sample, homozygote excess 0.4, window ±1, concentration from pairs of samples

Rounds of the last counting: median 62, most 68. 200 draws, 50 bootstrap replicates each, 92 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0798 | 0.0012 | 76 / 94 | 76 / 96 | 78 / 96 | 0.0013 / 0.0013 / 0.0013 |
| shorter_share | 0.8000 | 0.7993 | 0.0074 | 63 / 96 | 64 / 94 | 64 / 94 | 0.0070 / 0.0070 / 0.0070 |
| fall_off | 0.2500 | 0.2510 | 0.0137 | 63 / 94 | 62 / 94 | 64 / 94 | 0.0129 / 0.0131 / 0.0131 |
| share_-1 | 0.2619 | 0.2615 | 0.0092 | 72 / 94 | 71 / 94 | 70 / 94 | 0.0091 / 0.0092 / 0.0092 |
| share_0 | 0.4762 | 0.4773 | 0.0111 | 68 / 94 | 66 / 94 | 66 / 94 | 0.0104 / 0.0105 / 0.0105 |
| share_1 | 0.2619 | 0.2612 | 0.0100 | 60 / 94 | 61 / 94 | 62 / 94 | 0.0092 / 0.0092 / 0.0092 |
| concentration | 0.5000 | 0.5032 | 0.0468 | 64 / 94 | 64 / 94 | 64 / 94 | 0.0441 / 0.0443 / 0.0443 |
| between_locus_share | 0.6667 | 0.6659 | 0.0207 | 64 / 94 | 64 / 93 | 64 / 93 | 0.0194 / 0.0196 / 0.0196 |

#### one_sample_reads_30_own_copies — 2000 tracts, 1 samples, 30 reads a sample, homozygote excess 0.4, window ±1, concentration from each sample's two copies

Rounds of the last counting: median 62, most 68. 200 draws, 50 bootstrap replicates each, 92 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0798 | 0.0012 | 76 / 94 | 76 / 96 | 78 / 96 | 0.0013 / 0.0013 / 0.0013 |
| shorter_share | 0.8000 | 0.7993 | 0.0074 | 63 / 96 | 64 / 94 | 64 / 94 | 0.0070 / 0.0070 / 0.0070 |
| fall_off | 0.2500 | 0.2510 | 0.0137 | 63 / 94 | 62 / 94 | 64 / 94 | 0.0129 / 0.0131 / 0.0131 |
| share_-1 | 0.2619 | 0.2615 | 0.0092 | 72 / 94 | 71 / 94 | 70 / 94 | 0.0091 / 0.0092 / 0.0092 |
| share_0 | 0.4762 | 0.4773 | 0.0111 | 68 / 94 | 66 / 94 | 66 / 94 | 0.0104 / 0.0105 / 0.0105 |
| share_1 | 0.2619 | 0.2612 | 0.0100 | 60 / 94 | 61 / 94 | 62 / 94 | 0.0092 / 0.0092 / 0.0092 |
| concentration | 0.5000 | 0.5032 | 0.0468 | 64 / 94 | 64 / 94 | 64 / 94 | 0.0441 / 0.0443 / 0.0443 |
| between_locus_share | 0.6667 | 0.6659 | 0.0207 | 64 / 94 | 64 / 93 | 64 / 93 | 0.0194 / 0.0196 / 0.0196 |

#### one_sample_reads_3 — 2000 tracts, 1 samples, 3 reads a sample, homozygote excess 0.4, window ±1, concentration from pairs of samples

Rounds of the last counting: median 477, most 1036. 200 draws, 50 bootstrap replicates each, 511 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0814 | 0.0118 | 64 / 92 | 64 / 92 | 64 / 92 | 0.0108 / 0.0110 / 0.0109 |
| shorter_share | 0.8000 | 0.8023 | 0.0503 | 68 / 92 | 68 / 96 | 68 / 96 | 0.0514 / 0.0511 / 0.0505 |
| fall_off | 0.2500 | 0.2529 | 0.0664 | 70 / 94 | 70 / 96 | 68 / 96 | 0.0681 / 0.0681 / 0.0673 |
| share_-1 | 0.2619 | 0.2610 | 0.0105 | 64 / 92 | 66 / 94 | 66 / 94 | 0.0100 / 0.0102 / 0.0102 |
| share_0 | 0.4762 | 0.4775 | 0.0120 | 64 / 94 | 65 / 96 | 65 / 96 | 0.0114 / 0.0115 / 0.0115 |
| share_1 | 0.2619 | 0.2615 | 0.0099 | 67 / 96 | 64 / 96 | 64 / 96 | 0.0099 / 0.0100 / 0.0100 |
| concentration | 0.5000 | 0.4968 | 0.1910 | 68 / 94 | 66 / 94 | 66 / 95 | 0.1865 / 0.1821 / 0.1805 |
| between_locus_share | 0.6667 | 0.6786 | 0.0842 | 66 / 92 | 66 / 94 | 66 / 94 | 0.0811 / 0.0811 / 0.0804 |

#### one_sample_reads_3_own_copies — 2000 tracts, 1 samples, 3 reads a sample, homozygote excess 0.4, window ±1, concentration from each sample's two copies

Rounds of the last counting: median 477, most 1036. 200 draws, 50 bootstrap replicates each, 859 s. Cells: % of draws with the truth within 1 / 2 errors.

| number | truth | mean estimate | spread of estimates | bootstrap | sandwich | naive | mean error: bootstrap / sandwich / naive |
|---|---|---|---|---|---|---|---|
| level | 0.0800 | 0.0814 | 0.0118 | 64 / 92 | 64 / 92 | 64 / 92 | 0.0108 / 0.0110 / 0.0109 |
| shorter_share | 0.8000 | 0.8023 | 0.0503 | 68 / 92 | 68 / 96 | 68 / 96 | 0.0514 / 0.0511 / 0.0505 |
| fall_off | 0.2500 | 0.2529 | 0.0664 | 70 / 94 | 70 / 96 | 68 / 96 | 0.0681 / 0.0681 / 0.0673 |
| share_-1 | 0.2619 | 0.2610 | 0.0105 | 64 / 92 | 66 / 94 | 66 / 94 | 0.0100 / 0.0102 / 0.0102 |
| share_0 | 0.4762 | 0.4775 | 0.0120 | 64 / 94 | 65 / 96 | 65 / 96 | 0.0114 / 0.0115 / 0.0115 |
| share_1 | 0.2619 | 0.2615 | 0.0099 | 67 / 96 | 64 / 96 | 64 / 96 | 0.0099 / 0.0100 / 0.0100 |
| concentration | 0.5000 | 0.4968 | 0.1910 | 68 / 94 | 66 / 94 | 66 / 95 | 0.1865 / 0.1821 / 0.1805 |
| between_locus_share | 0.6667 | 0.6786 | 0.0842 | 66 / 92 | 66 / 94 | 66 / 94 | 0.0811 / 0.0811 / 0.0804 |
