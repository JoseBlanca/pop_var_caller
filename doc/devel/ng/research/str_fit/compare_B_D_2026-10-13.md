# The counting fit against penalised counting, head to head

**Status:** research report, 2026-10-13, branch `str-fit-compare` (tracks 1 and 2 merged). The
comparison the research plan's phase 2 asks for ([`../str_fit_research_plan.md`](../str_fit_research_plan.md)
§5) between avenue **B**, the counting fit class by class with the caller's curves drawn afterwards,
and avenue **D**, penalised counting, which fits every class of a motif length at once with each
number a smooth curve across repeat count. The climb **A**, which the caller ships, was run where one
fit finishes in about half an hour. **Every figure is on drawn data** from the fitting bench of
phase 0 ([`phase0_bench_2026-10-12.md`](phase0_bench_2026-10-12.md)); none is on a real cohort.

---

## 1. The answer

**Which method predicts held-out loci better depends on the number of samples, and the reason is one
detail of D, not penalised counting itself.** The score is the log-likelihood of the reads at loci the
fit never saw, minus what the true model scores on the same loci, in log-likelihood units a locus:
zero is as good as the truth. "D − B" is the difference between the two methods' scores, seed by seed
on the same draws, with its standard error over seeds; above zero, D predicts better. B is taken after
the caller's curves, its better form.

| samples | reads | one-base runs: D − B | two-base runs: D − B | better |
|---|---|---|---|---|
| 1 | 3 to 300 | +0.10 to +0.15 (± 0.01) | +0.14 to +0.16 (± 0.01) | **D**, four times closer to the truth |
| 10 | 3 | +0.045 ± 0.015 | +0.061 ± 0.010 | **D** |
| 10 | 10 to 300 | +0.015 to +0.025 (± 0.01) | +0.001 to +0.025 (± 0.02) | level, D slightly ahead |
| 60 | 3 to 300 | −0.12 to −0.16 (± 0.02) | −0.21 to −0.26 (± 0.03) | **B** |
| 500 | 3 and 10 | −0.50 to −0.56 (± 0.04) | −0.83 to −0.92 (± 0.05) | **B** |
| 2,000 | 3 | −0.62 (one seed) | −1.17 ± 0.10 (two seeds) | **B** |

At one sample both methods are scored with the class's own allele-length shares kept (§5); with the
fallback calling uses today, both lose most of what they fitted.

**Why D loses from 60 samples up: the floor it puts under the allele-length shares does not grow with
the cohort.** D adds one allele draw a class, spread over the 19 lengths, so that no length's share
reaches exactly zero. B adds one *locus's* worth of draws (track 1's fix). At one sample a locus holds
about two draws and the two floors are alike; at 60 samples a locus holds about 115 and D's floor is
a hundred times smaller. A length that no training locus carries then gets a share near 10⁻⁶ where B
gives it about 5 × 10⁻⁴, and a held-out locus fixed at that length costs several units. At 60 samples
and 10 reads (two-base runs, seed 1) D placed 140 of its 399 shares below a tenth of their truth; B
placed 7.

**With B's floor, D predicts better than B in every cell where that was run.** `D_tract_prior` is D
with B's prior on the shares added and nothing else changed (§3.3):

| cell | D − B | D with B's prior − B |
|---|---|---|
| one-base runs, 10 samples, 3 reads | +0.045 ± 0.015 | **+0.100 ± 0.014** |
| two-base runs, 10 samples, 3 reads | +0.061 ± 0.010 | **+0.149 ± 0.012** |
| one-base runs, 60 samples, 10 reads | −0.157 ± 0.018 | **+0.027 ± 0.004** |
| two-base runs, 60 samples, 10 reads | −0.258 ± 0.010 | **+0.040 ± 0.006** |
| one-base runs, 500 samples, 3 reads (seed 1) | −0.56 | **+0.010** |
| two-base runs, 500 samples, 3 reads (seed 1) | −0.83 | **+0.026** |
| one sample, 3 reads, both motif lengths (3 seeds) | | the same as D, within 0.008 a locus |

So the method the measurements favour is **penalised counting with the counting fit's prior on the
shares**. The prior's size, one locus's worth, was chosen in track 1 for B and was not tuned here.

**Recovering the true numbers, D does as well as B or better in nearly every cell** (§4): at one
sample and three reads its slip rate is about half as far from the truth (relative error 0.087
against 0.176, one-base runs); from 10 samples up the two are within a few per cent of each other on
one-base runs and D is ahead on two-base runs (0.0042 against 0.0087 at 2,000 samples). The
reference length's share is where D's floor shows again: no better than B's from 500 samples up.

**Both methods' standard errors are honest for the three stutter numbers and F** (§6): the truth
within one reported error in 43 to 86 estimates of 100 and within two in 64 to 100, around the 68 and
95 that honest errors give, with two exceptions: B's slip rate at 500 and 2,000 samples, biased by
B's prior as track 1 §12 found (within one error 0.38 to 0.43 times in 100 on one-base runs), and D's
shorter share at three reads (0.43 to 0.56). **Neither is honest for shares below 0.01 beyond one
sample** (within two errors 0.21 to 0.60 times in 100), D the worse from 60 samples up; at one sample
D's small shares are honest and B's are not.

**Cost and the budget** (§7). Alone on this Mac the two fits cost the same to within 40%: at 2,000
samples and three reads, 150 to 210 seconds a motif length at 100 loci a class. **Scaled to kimura
(2,169 samples, about 86,300 loci), neither fits in one hour: about 5 hours for B and 6 for D on this
machine's 18 cores**, if the cost grows with the sample-locus rows as it did here. The real counting
fit on kimura at ±6 took 15.5 minutes, which this per-row cost does not explain (§7.3), so the budget
answer is uncertain by up to that factor. **Memory fits in 64 GiB only if classes are fitted one at a
time**: the table both methods build to read F from pairs of samples takes 2.9 kB a sample-locus row,
about 28 GiB for one class of 5,000 loci at kimura's scale. D builds one class's table at a time
(peak 3.3 GiB at 2,000 samples and 400 loci a class); the bench's B builds every class's at once
(20.6 GiB on the same cell).

**The climb (A)**, run at one sample and at 10 samples × 3 reads, seed 1 (§8), predicts held-out
loci worse than D in every cell (−0.27 to −0.46 a locus against −0.03 to −0.14) and worse than B
whenever B keeps its own shares; it is better than B only under calling's fallback at one sample,
because it always reports an F. One fit took 450 to 3,200 seconds under load.

**Not run** (§10): the held-out score at 2,000 samples for 10, 30 and 300 reads and at 500 samples for
30 and 300 reads (recovery and errors were run at 500 × 30, 500 × 300 and 2,000 × 10), 2,000 × 30 and
2,000 × 300 entirely, and the climb beyond 10 samples.

---

## 2. What was compared

### 2.1 Vocabulary

- **Class.** Every locus with one motif length and one reference repeat count. One-base runs
  (period 1) have 14 classes, 8 to 21 repeats; two-base runs (period 2) 21 classes, 6 to 27 repeats.
- **Held-out score.** Each class's loci are split into five groups; each group is scored under a fit
  made on the other four. A locus's score is the log of the probability of all its samples' reads,
  its own allele frequencies averaged out. The truth's score on the same loci is subtracted.
- **Standard error over seeds.** The spread of a figure across independently drawn cohorts, divided
  by the square root of their number: "± 0.01" means the figure would move by about that if the
  experiment were repeated.
- **Paired difference.** D's score minus B's on the same drawn cohort and the same held-out loci,
  averaged over seeds. It removes the cohort-to-cohort noise the two share.
- **F, the between-locus share.** How far a locus's own allele frequencies depart from its class's
  typical shape: at 0 every locus has the class's shape, at 1 every locus is fixed at one length.
- **Calling's fallback.** When a class has no F of its own, calling seeds its genotype prior from a
  flat shape at a stated F and drops the class's own shares (track 1 §12.2).

### 2.2 The methods

| name | what it is |
|---|---|
| `B_raw` | the counting fit class by class, F read off pairs of samples, sandwich errors, shares' prior of one locus's worth (track 1) |
| `B_curves` | `B_raw` after the caller's curves across repeat count, which replace the three stutter numbers |
| `D` | penalised counting as track 2 left it: every class of a motif length at once, each number a penalised curve, F off pairs of samples, its own errors; F absent at one sample |
| `D_working_F` | D reporting at one sample the F it worked at, so that every class is scored at its own shares |
| `D_tract_prior` | D with the counting fit's shares prior added (`PenalisedFitConfig::allele_prior_tracts = 1`); a measuring variant, added on this branch |
| `A_raw`, `A_curves` | the climb as the caller runs it, before and after the curves |

### 2.3 The cells

All on the bench's default truth (`CohortModel::tomato_like`, PCR stutter, F 0.7, homozygote excess
0.85 to 0.99 a sample), at most 100 loci a class, the shipped ±9 allele window:

| samples | reads | seeds | held-out score |
|---|---|---|---|
| 1 | 3, 10, 30, 300 | 6 | every held-out locus, both fallbacks |
| 10, 60 | 3, 10, 30, 300 | 3 | every held-out locus |
| 500 | 3, 10 | 3 | 20 loci a class |
| 500 | 30, 300 | 3 | not run; recovery and errors only |
| 2,000 | 3 | one-base 2, two-base 3 | 20 loci a class, seed 1 (one-base) and seeds 1–2 (two-base); seed 3 recovery and errors only |
| 2,000 | 10 | 2 (seeds 2, 3) | not run; recovery and errors only |
| 2,000 | 30, 300 | — | not run |

Two more truths. **PCR-free** (`tomato_like(.., PcrFree)`, every slip rate a tenth): at 1 × 3, 1 × 300
and 60 × 10 with the held-out score, 3 seeds, both motif lengths, and at 2,000 × 3, seed 1, without
it. **Two-base runs with F = 0.4**: the default truth with its F replaced by 0.4 in every class
(`--between-locus-share 0.4`; tomato2's two-base classes give 0.29 to 0.50), period 2 only, at the
same cells.

---

## 3. Held-out prediction

### 3.1 By cohort size

The table of §1 comes from the appendix's *Held-out prediction, side by side*. In order of cohort:

- **At one sample** D sits 0.03 to 0.05 units a locus from the truth in every cell, truth and motif
  length; B after the curves 0.14 to 0.20. Both are scored with their own shares (§5).
- **At 10 samples** D is ahead by 0.05 to 0.06 a locus at three reads and level, within about one
  standard error, at 10 to 300 reads.
- **From 60 samples** B is ahead, and by more as the cohort grows: 0.12 to 0.26 a locus at 60
  samples, 0.5 to 0.9 at 500, 0.6 to 1.2 at 2,000. Two-base runs lose more than one-base runs at
  every size.

The other truths leave this unchanged where they were run. At 60 samples and 10 reads D − B is
−0.139 ± 0.016 (one-base, PCR-free), −0.256 ± 0.009 (two-base, PCR-free) and −0.219 ± 0.011
(two-base, F 0.4), against −0.157 and −0.258 on the default truth.

### 3.2 Where D's loss comes from

At 60 samples and 10 reads, two-base runs, seed 1, of the 399 shares fitted (21 classes × 19
lengths): D put 110 below 10⁻⁴ and 140 below a tenth of their truth; B none below 10⁻⁴ and 7 below a
tenth. At 10 samples the counts are 40 and 51 for D, 0 and 6 for B. Track 2 saw the same shares
(its §10.2: "a share of 0.004 is reported as about 10⁻⁶").

### 3.3 D with B's prior

`D_tract_prior` adds to each length's allele draws, every round, the class's draws divided by its
loci with reads and by the 19 lengths: the counting fit's rule (`SHARE_PRIOR_TRACTS`). D's own floor
stays. It was run at 1 × 3, 10 × 3 and 60 × 10 (3 seeds) and 500 × 3 (seed 1), both motif lengths,
with the figures in §1's second table. With it D overtakes B at 60 and 500 samples and doubles its
lead at 10; at one sample it changes nothing that matters.

**What this does not settle.** The 500-sample figure is one seed. The prior costs B about 0.8% of
slip rate at long, thin classes and 2,000 samples (track 1 §12.1); whether it costs D the same was
not measured, and nor were its effects on D's errors for small shares beyond the appendix's tables.

---

## 4. Recovery

The median over classes of each class's root-mean-square error over the seeds; the slip rate relative
to its truth, the rest absolute. Selected cells of the default truth; every cell is in the appendix.

| cell | number | B after curves | D |
|---|---|---|---|
| one-base, 1 × 3 | slip rate | 0.176 | **0.087** |
| | shorter share | 0.094 | **0.075** |
| | reference-length share | 0.060 | **0.042** |
| two-base, 1 × 3 | slip rate | 0.153 | **0.133** |
| | reference-length share | 0.087 | **0.049** |
| one-base, 1 × 300 | slip rate | 0.020 | 0.020 |
| two-base, 10 × 3 | slip rate | 0.078 | **0.045** |
| | F | 0.059 | **0.022** |
| one-base, 60 × 3 | slip rate | **0.019** | 0.022 |
| two-base, 60 × 10 | slip rate | 0.020 | **0.013** |
| two-base, 500 × 3 | slip rate | 0.013 | **0.010** |
| two-base, 500 × 10 | reference-length share | **0.070** | 0.076 |
| one-base, 2,000 × 3 | slip rate | 0.0083 | **0.0063** |
| two-base, 2,000 × 3 | slip rate | 0.0087 | **0.0042** |
| | reference-length share | **0.067** | 0.077 |

**D recovers the stutter as well or better in nearly every cell**, most at one sample and on two-base
runs, whose classes are thinner. F is recovered alike from 60 samples up; at 10 samples D's is better
on two-base runs. At one sample D reports no F; B reports one in the classes where it can place it
from a sample's own two copies (none to 9 of 14 or 21), and those are far off (0.22 to 0.48
root-mean-square, truth 0.7). The PCR-free truth gives the same order: at 2,000 × 3 the slip rate
0.017 (B) against 0.010 (D) on one-base runs.

## 5. One sample: the two fallbacks

At one sample no method can read F from pairs of samples. B places F in some classes from the
sample's own two copies; D never reports it. Calling's fallback, which drops a class's own shares when
it has no F, therefore decides most of the one-sample score:

| cell | B after curves: calling's fallback / own shares | D: calling's fallback / own shares | D at its working F | held-out loci whose class B left without F |
|---|---|---|---|---|
| one-base, 1 × 3 | −1.524 / −0.141 | −1.628 / **−0.037** | −0.038 | 91% |
| one-base, 1 × 300 | −0.942 / −0.157 | −1.653 / **−0.041** | −0.040 | 49% |
| two-base, 1 × 3 | −1.802 / −0.186 | −1.828 / **−0.047** | −0.048 | 96% |
| two-base, 1 × 300 | −1.252 / −0.196 | −1.885 / **−0.049** | −0.051 | 65% |

Six seeds each; standard errors 0.002 to 0.12, in the appendix's fallback table with 10 and 30 reads
and the other truths. **Under calling's fallback both methods are 0.9 to 1.9 units a locus from the
truth, and D is the worse because no class of it has an F. With the class's own shares kept, both
improve five- to fortyfold and D is four times closer to the truth than B.** D at its working F
scores as D with its shares kept: at one sample the F used hardly matters, as track 1 §13 found for
B. Track 1's case for a rung in calling's ladder that keeps a class's own shares when it has no F
applies to D as much as to B.

## 6. Error honesty

The share of estimates whose truth lies within one / two reported errors, pooled over classes and
seeds; honest errors give about 0.68 / 0.95. Ranges over the default truth's cells; every cell, with
the counts behind it, is in the appendix (*Error honesty, side by side*).

| number | B, 1 sample | D, 1 sample | B, 10 to 2,000 samples | D, 10 to 2,000 samples |
|---|---|---|---|---|
| slip rate | 0.59–0.71 / 0.90–0.97 | 0.60–0.73 / 0.87–0.98 | 0.38–0.86 / 0.68–0.98 | 0.48–0.76 / 0.81–0.98 |
| shorter share | 0.62–0.70 / 0.87–0.98 | 0.47–0.82 / 0.83–1.00 | 0.48–0.79 / 0.82–1.00 | 0.38–0.90 / 0.62–1.00 |
| fall-off | 0.54–0.76 / 0.86–0.98 | 0.36–0.77 / 0.85–0.98 | 0.40–0.76 / 0.82–1.00 | 0.45–0.83 / 0.69–1.00 |
| F | few classes | absent | 0.50–0.78 / 0.83–1.00 | 0.43–0.83 / 0.78–1.00 |
| shares ≥ 0.05 | 0.64–0.71 / 0.88–0.93 | 0.55–0.69 / 0.86–0.94 | 0.54–0.68 / 0.84–0.93 | 0.52–0.68 / 0.86–0.93 |
| shares 0.01–0.05 | 0.35–0.51 / 0.41–0.60 | 0.59–0.71 / 0.86–0.93 | 0.38–0.57 / 0.50–0.70 | 0.37–0.55 / 0.49–0.72 |
| shares < 0.01 | 0.18–0.26 / 0.26–0.35 | 0.56–0.63 / 0.87–0.93 | 0.25–0.38 / 0.35–0.52 | 0.18–0.45 / 0.21–0.60 |

**The stutter numbers and F are honest for both methods**, within the noise of 14 to 126 estimates a
cell, with three departures:

- **B's slip rate at many samples and three reads.** Within one error 0.43 / 0.68 at 2,000 samples and
  0.38 / 0.83 at 500 (one-base runs): the bias from B's shares prior that track 1 §12.1 traced, about
  0.8% of slip rate at long, thin classes, which the errors do not cover. D's slip rate there: 0.61 /
  0.86 and 0.48 / 0.88.
- **D's shorter share at three reads**: 0.43 to 0.56 within one error at 1, 10 and 60 samples and
  0.38 at 500, one-base runs; its errors there are too small by about a third.
- **D's fall-off at one sample and 10 reads**: 0.36 / 0.85.

**Small shares are not honest for either method beyond one sample**, as tracks 1 and 2 found: a share
of 0.005 at 100 loci and F 0.7 is carried by a fraction of one locus, and an error read off the loci's
spread cannot see a length no locus carried. D is the worse of the two from 60 samples up (shares
below 0.01 within two errors 0.21 to 0.36 times in 100 against B's 0.36 to 0.52), the floor again; at
one sample D's are honest and B's are not. Shares of 0.05 or more are honest for both everywhere.

## 7. Cost

### 7.1 Measured alone

Each fit alone on this Mac (Apple silicon, 18 cores, 64 GiB), seed 1, at most 100 loci a class unless
the cell says otherwise, without the held-out score: the fit on every locus, its standard errors
included. **Each started with the 1-minute load below 3** (the runner waited for it, and watched
memory every 10 s, ready to stop a fit if less than 8 GiB stayed free; none came near). An earlier
attempt was taken at a load of 18 to 20 and ended when the machine crashed; its timings are
discarded. "Draw alone" is the same cell run with no method: the drawn cohort's own memory.

| cell | rows with reads | B: fit, s / CPU, s / peak, MiB | D: fit, s / CPU, s / peak, MiB | draw alone, MiB |
|---|---|---|---|---|
| one-base, 1 × 3 | 911 | 0.1 / 1.2 / 10 | 0.3 / 2.5 / 16 | 7 |
| one-base, 1 × 300 | 971 | 0.3 / 4.4 / 10 | 0.1 / 0.7 / 15 | 7 |
| two-base, 1 × 3 | 828 | 0.1 / 1.4 / 10 | 0.3 / 2.0 / 22 | 7 |
| two-base, 1 × 300 | 878 | 0.05 / 0.8 / 10 | 0.2 / 0.8 / 22 | 7 |
| one-base, 60 × 10 | 58,259 | 2.6 / 42 / 237 | 2.5 / 40 / 45 | 16 |
| two-base, 60 × 10 | 52,679 | 2.5 / 41 / 115 | 2.8 / 39 / 52 | 15 |
| one-base, 500 × 3 | 461,435 | 40.5 / 659 / 1,167 | 39.3 / 665 / 348 | 72 |
| two-base, 500 × 3 | 417,234 | 37.5 / 629 / 589 | 57.3 / 944 / 252 | 66 |
| one-base, 2,000 × 3 | 1,845,379 | 154.7 / 2,636 / 5,287 | 147.4 / 2,502 / 984 | 268 |
| two-base, 2,000 × 3 | 1,668,546 | 147.7 / 2,465 / 3,634 | 210.4 / 3,506 / 955 | 243 |
| one-base, 2,000 × 3, 200 loci a class | 3,088,597 | 294.4 / 4,954 / 6,183 | 264.6 / 4,632 / 1,753 | 442 |
| one-base, 2,000 × 3, 400 loci a class | 5,115,915 | 511.5 / 8,776 / 20,596 | 475.5 / 8,390 / 3,293 | 728 |
| two-base, 2,000 × 3, every tomato2 locus (1,328) | 2,523,976 | 232.6 / 3,913 / 5,514 | 296.1 / 5,013 / 3,042 | 363 |

**At one sample both fits take under half a second; from 60 samples up they cost within 40% of each
other**, level on one-base runs and B ahead on two-base runs (148 against 210 s at 2,000 × 3). CPU time
is about 17 times the fit's wall time at the large cells: both use the whole machine.

**Memory is where they differ: B's peak is 3.6 to 6.2 GiB at 2,000 samples and 100 to 200 loci a
class and 20.6 GiB at 400; D's 1.0 to 3.3 GiB.** Both read F off pairs of samples, and for that each
builds, per class, a table of every sample's likelihood of every ordered pair of allele lengths: 19 ×
19 numbers of 8 bytes, 2.9 kB, a sample-locus row. D builds one class's table at a time; the bench's
B fits every class side by side (`TheCountingFit::fit_period` maps the classes in parallel) and so
holds them all at once. D's peaks match one class's table plus the drawn cohort: at 400 loci a class,
400 × 1,900 rows × 2.9 kB = 2.0 GiB plus 0.7 GiB for the draw, against 3.3 GiB measured.

### 7.2 The budget: kimura, 2,169 samples, about three reads

**The loci.** tomato2's catalog (`TOMATO2_LOCI`: 14,589 loci over 32 Mb, classes with at least 8
loci) scaled by 31.7, the ratio of tomato's whole-genome repeat-tract loci (462,701,
`spec/parameter_prepass_joint_loci.md`) to tomato2's, and capped at the pre-pass's 5,000 loci a class:
43,478 one-base, 32,867 two-base, 86,303 in all, within 0.5% of the 86,688 tracts the real kimura run
fitted. It leaves out classes in which tomato2 has fewer than 8 loci and the whole genome may have
more. At three reads 95% of samples have a read at a locus: 178 million sample-locus rows.

**Time.** The cost a row grew slowly with the loci a class: B 1.43, 1.60 and 1.72 CPU-milliseconds at
100, 200 and 400 loci a class (one-base runs), D 1.36, 1.50 and 1.64; at every tomato2 locus of the
two-base runs B 1.55 and D 1.99. At 1.72 ms for B and 1.99 ms for D, and the measured 17-fold
speed-up of 18 cores:

| motif length | rows | B | D |
|---|---|---|---|
| one-base | 89.6 million | 2.5 hours | 2.9 hours |
| two-base | 67.7 million | 1.9 hours | 2.2 hours |
| all | 178 million | **about 5.0 hours** | **about 5.8 hours** |

**By this measure neither fits within one hour on this machine**: each would need five to six times
less work a row, or a machine that much faster.

**Memory.** The drawn cohort costs 142 bytes a row here, 23.5 GiB for kimura's rows. One class of
5,000 loci adds a pair table of 5,000 × 2,061 rows × 2.9 kB = **27.7 GiB**. D, a class at a time,
would peak near 23.5 + 27.7 = 51 GiB, inside 64 GiB; less if each motif length's evidence were loaded
alone. B as the bench runs it would hold several such tables at once and would not fit; fitted a class
at a time it would peak like D. A table of the 190 unordered genotypes in place of the 361 ordered
pairs would halve the 27.7 GiB.

### 7.3 Against the real kimura run

The counting fit on kimura's own evidence (86,688 tracts, 2,169 samples, ±6 window) took 15.5 minutes
and 31.5 GiB ([`../str_fit_research_plan.md`](../str_fit_research_plan.md) §1.1). The ±9 window
doubles the genotypes (190 against 91), about a factor of two; the extrapolation above is still ten
times that run. Phase 0 found the bench slower per row than kimura (`phase0_bench_2026-10-12.md` §4)
without explaining it. They differ in machine, in the standard errors (computed here; 0.01 to 0.3
times the fit by track 1 §8) and possibly in the rounds the bench's truth needs; none was measured
here. **So the budget answer is "no, by five- to sixfold" on the bench's cost a row, and uncertain by
up to that factor until both fits are timed on kimura's own evidence at ±9.**

## 8. The climb (A)

Run at one sample (3 and 300 reads) and at 10 samples × 3 reads (one-base runs), seed 1, with the
held-out score. Each run took 1.9 to 2.4 hours at one sample and 3.7 hours at 10 samples under load,
the climb's six fits taking nearly all of it. Held-out score against the truth, a locus, on the same
draws as B and D:

| cell | A as fitted / after curves | B after curves: calling / own shares | D own shares | D at working F |
|---|---|---|---|---|
| one-base, 1 × 3 | −0.316 / −0.273 | −1.685 / −0.139 | −0.027 | −0.033 |
| one-base, 1 × 300 | −0.325 / −0.312 | −0.837 / −0.153 | −0.049 | −0.049 |
| two-base, 1 × 3 | −0.577 / −0.463 | −1.735 / −0.205 | −0.058 | −0.054 |
| two-base, 1 × 300 | −0.348 / −0.332 | −1.098 / −0.166 | −0.051 | −0.056 |
| one-base, 10 × 3 | −0.331 / −0.318 | −0.180 | −0.140 (D as fitted) | |

**The climb always reports an F** (biased low by 0.29 to 0.62 at one sample, phase 0 §3.1), so
calling's fallback never touches it; that is why it beats B under the fallback. Against D, or B with
its own shares, it is worse by 0.1 to 0.4 a locus. Its recovery at these cells is in the appendix
(one seed, so not paired with B's and D's six). One fit took 450 to 1,620 seconds at one sample and
3,180 at 10 samples, under load, against under a second (one sample) and under 4 seconds (10
samples) for B and D.

## 9. The merge, and what reproduced

`str-fit-b` (track 1) was merged first, then `str-fit-d` (track 2). Both had split the counting fit's
round out of `composite_fit.rs`. The merge keeps track 1's file and adds what penalised counting reads
from it: the round's counts as the public `ClassCounts`, `expected_counts` (the sum of track 1's
`expected_counts_by_tract` over a class), track 2's `SamplePairs` copied unchanged, and
`CountedClass::rounds_in_all`. Penalised counting's check against the counting fit now calls the
counting fit without its shares prior (`fit_class_by_counting_at_share_prior` at zero), so both
maximise the same objective as before.

**Every number checked reproduced exactly.** All 111 tests of `ssr_fit` pass. The bench's rows match
each track's own runs value for value, times aside: D at 1 × 3, 1 × 300 and 60 × 10, both motif
lengths, seeds 1–6 (36 runs, track 2's honesty runs); B at one sample, 3 and 300 reads, seeds 1–3, at
calling's fallback and with own shares kept (12 runs, track 1 §13). Track 2's drawn-grid example at
1 × 300 gives the same penalised rows; its counting-fit rows changed, as they now carry track 1's
fixes.

Added on this branch: `--both-fallbacks` (both held-out fallbacks in one run), cell names that say
when `--between-locus-share` replaced the truth's F, `D_tract_prior`, and the summariser's
`--compare` tables over several directories.

## 10. Not done, and uncertain

- **Held-out score not run** at 2,000 samples for 10, 30 and 300 reads, nor at 500 samples for 30 and
  300 reads: one held-out run at 2,000 samples took about three hours on the shared machine. A second
  one-base seed at 2,000 × 3 was started and stopped for time. Recovery and errors were run at
  500 × 30 and 500 × 300 (3 seeds) and 2,000 × 10 (seeds 2 and 3); 2,000 × 30 and 2,000 × 300 were
  not run at all.
- **The other truths at 2,000 samples** were run without the held-out score, seed 1 only.
- **The climb** was run at one sample and 10 samples × 3 reads, seed 1 only.
- **`D_tract_prior`** at 500 samples is one seed; its prior's size is borrowed from B, not tuned.
- **The budget** rests on the bench's cost a row, which disagrees with the real kimura run (§7.3).
- **Times in the appendix's *Cost* table were taken with up to six runs sharing the machine**, and
  other users' work; only §7.1's are quiet.
- **One family of drawn truths.** How either method does where real classes depart from smooth curves
  (a kink such as HG002's 9-repeat homopolymers) is track 2 §3's measurement, not repeated here.

## 11. How to rerun

```sh
cargo build --release --features bench-fixtures --example ng_str_fit_bench
# one cell, B, D and D with B's prior on the same draws, both fallbacks
target/release/examples/ng_str_fit_bench --period 2 --samples 60 --depth 10 --seed 1 --truth pcr \
    --loci-at-most 100 --methods B,D,D_tract_prior --both-fallbacks --out rows.tsv
# the F = 0.4 truth
target/release/examples/ng_str_fit_bench --period 2 --samples 60 --depth 10 --seed 1 --truth pcr \
    --between-locus-share 0.4 --loci-at-most 100 --methods B,D --both-fallbacks --out rows.tsv
# the tables
uv run --no-project --with numpy python scripts/research/str_fit_summarise.py DIR [DIR ...] \
    --compare B_raw,B_curves,D,D_working_F,D_tract_prior
```

The grid's job lists and runners and every rows file are in this worktree's `tmp/` (`tmp/grid`,
`tmp/grid_b`, `tmp/gridA`, `tmp/cost_quiet`), not committed; the draws are seeded and reproduce.

---

## Appendix: the full tables

Generated by `scripts/research/str_fit_summarise.py tmp/grid tmp/grid_b tmp/gridA --compare
A_raw,A_curves,B_raw,B_curves,D,D_working_F,D_tract_prior`; `—` where a method was not run in that
cell or gave no estimate. "A − B_curves" and the like are paired only on the seeds both ran.

### Held-out prediction, side by side

Log-likelihood a held-out locus under each fit minus under the truth (zero is as good as the truth; below zero is worse), at calling's fallback for a class without its own numbers. Mean over the seeds ± its standard error. The last columns are the paired difference, seed by seed, of D against B: above zero means D predicts better.

| cell | seeds | loci a seed | A_raw | A_curves | B_raw | B_curves | D | D_working_F | D_tract_prior | D − B_raw | A_raw − B_curves | A_curves − B_curves | D − B_curves | D_working_F − B_curves | D_tract_prior − B_curves |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| period 1, 1 samples, 10 reads, pcr | 6 | 971 | — | — | -1.087 ± 0.088 | -1.092 ± 0.088 | -1.650 ± 0.015 | -0.036 ± 0.008 | — | -0.564 ± 0.082 | — | — | -0.559 ± 0.082 | +1.056 ± 0.095 | — |
| period 1, 1 samples, 3 reads, pcr | 6 | 921 | -0.316 | -0.273 | -1.537 ± 0.043 | -1.524 ± 0.042 | -1.628 ± 0.021 | -0.038 ± 0.002 | -1.616 ± 0.043 | -0.091 ± 0.036 | +1.368 | +1.412 | -0.104 ± 0.037 | +1.486 ± 0.042 | -0.050 ± 0.045 |
| period 1, 1 samples, 3 reads, pcr_free | 3 | 921 | — | — | -1.510 ± 0.109 | -1.489 ± 0.113 | -1.657 ± 0.047 | -0.035 ± 0.001 | — | -0.148 ± 0.103 | — | — | -0.168 ± 0.108 | +1.454 ± 0.113 | — |
| period 1, 1 samples, 30 reads, pcr | 6 | 971 | — | — | -0.887 ± 0.083 | -0.868 ± 0.082 | -1.695 ± 0.020 | -0.030 ± 0.003 | — | -0.807 ± 0.073 | — | — | -0.826 ± 0.073 | +0.838 ± 0.084 | — |
| period 1, 1 samples, 300 reads, pcr | 6 | 971 | -0.325 | -0.312 | -0.954 ± 0.059 | -0.942 ± 0.058 | -1.653 ± 0.027 | -0.040 ± 0.002 | — | -0.699 ± 0.072 | +0.512 | +0.526 | -0.712 ± 0.072 | +0.902 ± 0.059 | — |
| period 1, 1 samples, 300 reads, pcr_free | 3 | 971 | — | — | -1.125 ± 0.145 | -0.975 ± 0.125 | -1.688 ± 0.036 | -0.036 ± 0.002 | — | -0.564 ± 0.177 | — | — | -0.713 ± 0.156 | +0.939 ± 0.125 | — |
| period 1, 10 samples, 10 reads, pcr | 3 | 971 | — | — | -0.219 ± 0.005 | -0.207 ± 0.006 | -0.191 ± 0.017 | — | — | +0.028 ± 0.014 | — | — | +0.016 ± 0.012 | — | — |
| period 1, 10 samples, 3 reads, pcr | 3 | 971 | -0.331 | -0.318 | -0.196 ± 0.020 | -0.186 ± 0.022 | -0.141 ± 0.007 | — | -0.087 ± 0.009 | +0.055 ± 0.013 | -0.151 | -0.138 | +0.045 ± 0.015 | — | +0.100 ± 0.014 |
| period 1, 10 samples, 30 reads, pcr | 3 | 971 | — | — | -0.221 ± 0.015 | -0.211 ± 0.014 | -0.186 ± 0.018 | — | — | +0.035 ± 0.009 | — | — | +0.025 ± 0.009 | — | — |
| period 1, 10 samples, 300 reads, pcr | 3 | 971 | — | — | -0.210 ± 0.003 | -0.206 ± 0.004 | -0.191 ± 0.013 | — | — | +0.019 ± 0.013 | — | — | +0.015 ± 0.014 | — | — |
| period 1, 2000 samples, 3 reads, pcr | 1 | 276 | — | — | -0.561 | -0.537 | -1.157 | — | — | -0.595 | — | — | -0.619 | — | — |
| period 1, 500 samples, 10 reads, pcr | 3 | 276 | — | — | -0.458 ± 0.036 | -0.461 ± 0.035 | -0.959 ± 0.045 | — | — | -0.501 ± 0.013 | — | — | -0.498 ± 0.019 | — | — |
| period 1, 500 samples, 3 reads, pcr | 3 | 276 | — | — | -0.441 ± 0.022 | -0.427 ± 0.016 | -0.989 ± 0.053 | — | -0.385 | -0.548 ± 0.035 | — | — | -0.562 ± 0.038 | — | +0.010 |
| period 1, 60 samples, 10 reads, pcr | 3 | 971 | — | — | -0.277 ± 0.025 | -0.264 ± 0.022 | -0.421 ± 0.011 | — | -0.237 ± 0.019 | -0.144 ± 0.020 | — | — | -0.157 ± 0.018 | — | +0.027 ± 0.004 |
| period 1, 60 samples, 10 reads, pcr_free | 3 | 971 | — | — | -0.279 ± 0.024 | -0.259 ± 0.024 | -0.399 ± 0.008 | — | — | -0.120 ± 0.016 | — | — | -0.139 ± 0.016 | — | — |
| period 1, 60 samples, 3 reads, pcr | 3 | 971 | — | — | -0.239 ± 0.004 | -0.227 ± 0.003 | -0.344 ± 0.010 | — | — | -0.105 ± 0.014 | — | — | -0.117 ± 0.013 | — | — |
| period 1, 60 samples, 30 reads, pcr | 3 | 971 | — | — | -0.250 ± 0.010 | -0.242 ± 0.010 | -0.391 ± 0.006 | — | — | -0.141 ± 0.016 | — | — | -0.149 ± 0.016 | — | — |
| period 1, 60 samples, 300 reads, pcr | 3 | 971 | — | — | -0.258 ± 0.038 | -0.256 ± 0.040 | -0.413 ± 0.047 | — | — | -0.154 ± 0.010 | — | — | -0.157 ± 0.009 | — | — |
| period 2, 1 samples, 10 reads, pcr | 6 | 878 | — | — | -1.614 ± 0.101 | -1.594 ± 0.101 | -1.903 ± 0.021 | -0.042 ± 0.005 | — | -0.289 ± 0.091 | — | — | -0.309 ± 0.090 | +1.552 ± 0.101 | — |
| period 2, 1 samples, 3 reads, pcr | 6 | 833 | -0.577 | -0.463 | -1.854 ± 0.025 | -1.802 ± 0.023 | -1.828 ± 0.020 | -0.048 ± 0.003 | -1.850 ± 0.039 | +0.026 ± 0.022 | +1.158 | +1.272 | -0.026 ± 0.021 | +1.754 ± 0.023 | -0.040 ± 0.031 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | 3 | 829 | — | — | -1.822 ± 0.059 | -1.778 ± 0.054 | -1.863 ± 0.020 | -0.049 ± 0.005 | — | -0.042 ± 0.076 | — | — | -0.085 ± 0.072 | +1.729 ± 0.058 | — |
| period 2, 1 samples, 3 reads, pcr_free | 3 | 830 | — | — | -1.980 ± 0.381 | -1.697 ± 0.110 | -1.880 ± 0.045 | -0.045 ± 0.009 | — | +0.100 ± 0.368 | — | — | -0.183 ± 0.085 | +1.652 ± 0.119 | — |
| period 2, 1 samples, 30 reads, pcr | 6 | 878 | — | — | -1.447 ± 0.114 | -1.412 ± 0.120 | -1.888 ± 0.026 | -0.046 ± 0.003 | — | -0.442 ± 0.121 | — | — | -0.476 ± 0.129 | +1.366 ± 0.119 | — |
| period 2, 1 samples, 300 reads, pcr | 6 | 878 | -0.348 | -0.332 | -1.269 ± 0.115 | -1.252 ± 0.116 | -1.885 ± 0.013 | -0.051 ± 0.003 | — | -0.615 ± 0.121 | +0.750 | +0.766 | -0.633 ± 0.122 | +1.201 ± 0.118 | — |
| period 2, 1 samples, 300 reads, pcr_F0.4 | 3 | 878 | — | — | -1.231 ± 0.043 | -1.219 ± 0.041 | -1.912 ± 0.028 | -0.053 ± 0.003 | — | -0.681 ± 0.017 | — | — | -0.693 ± 0.014 | +1.166 ± 0.039 | — |
| period 2, 1 samples, 300 reads, pcr_free | 3 | 878 | — | — | -1.605 ± 0.356 | -1.273 ± 0.149 | -1.856 ± 0.015 | -0.041 ± 0.002 | — | -0.251 ± 0.372 | — | — | -0.584 ± 0.164 | +1.232 ± 0.148 | — |
| period 2, 10 samples, 10 reads, pcr | 3 | 878 | — | — | -0.320 ± 0.018 | -0.299 ± 0.021 | -0.274 ± 0.010 | — | — | +0.046 ± 0.010 | — | — | +0.025 ± 0.013 | — | — |
| period 2, 10 samples, 3 reads, pcr | 3 | 878 | — | — | -0.288 ± 0.021 | -0.255 ± 0.020 | -0.193 ± 0.016 | — | -0.106 ± 0.009 | +0.095 ± 0.013 | — | — | +0.061 ± 0.010 | — | +0.149 ± 0.012 |
| period 2, 10 samples, 30 reads, pcr | 3 | 878 | — | — | -0.279 ± 0.022 | -0.265 ± 0.021 | -0.263 ± 0.022 | — | — | +0.015 ± 0.016 | — | — | +0.001 ± 0.019 | — | — |
| period 2, 10 samples, 300 reads, pcr | 3 | 878 | — | — | -0.284 ± 0.026 | -0.271 ± 0.026 | -0.268 ± 0.012 | — | — | +0.016 ± 0.019 | — | — | +0.003 ± 0.019 | — | — |
| period 2, 2000 samples, 3 reads, pcr | 2 | 388 | — | — | -0.744 ± 0.050 | -0.720 ± 0.029 | -1.889 ± 0.067 | — | — | -1.146 ± 0.117 | — | — | -1.169 ± 0.096 | — | — |
| period 2, 500 samples, 10 reads, pcr | 3 | 388 | — | — | -0.661 ± 0.009 | -0.650 ± 0.010 | -1.568 ± 0.054 | — | — | -0.907 ± 0.046 | — | — | -0.919 ± 0.048 | — | — |
| period 2, 500 samples, 3 reads, pcr | 3 | 388 | — | — | -0.575 ± 0.022 | -0.543 ± 0.018 | -1.374 ± 0.044 | — | -0.488 | -0.799 ± 0.047 | — | — | -0.831 ± 0.052 | — | +0.026 |
| period 2, 60 samples, 10 reads, pcr | 3 | 878 | — | — | -0.348 ± 0.030 | -0.332 ± 0.029 | -0.590 ± 0.027 | — | -0.292 ± 0.023 | -0.242 ± 0.011 | — | — | -0.258 ± 0.010 | — | +0.040 ± 0.006 |
| period 2, 60 samples, 10 reads, pcr_F0.4 | 3 | 878 | — | — | -0.391 ± 0.005 | -0.375 ± 0.007 | -0.594 ± 0.006 | — | — | -0.203 ± 0.010 | — | — | -0.219 ± 0.011 | — | — |
| period 2, 60 samples, 10 reads, pcr_free | 3 | 878 | — | — | -0.352 ± 0.031 | -0.326 ± 0.029 | -0.582 ± 0.031 | — | — | -0.231 ± 0.007 | — | — | -0.256 ± 0.009 | — | — |
| period 2, 60 samples, 3 reads, pcr | 3 | 878 | — | — | -0.321 ± 0.015 | -0.291 ± 0.018 | -0.502 ± 0.045 | — | — | -0.181 ± 0.032 | — | — | -0.211 ± 0.030 | — | — |
| period 2, 60 samples, 30 reads, pcr | 3 | 878 | — | — | -0.390 ± 0.043 | -0.377 ± 0.040 | -0.605 ± 0.021 | — | — | -0.215 ± 0.028 | — | — | -0.229 ± 0.026 | — | — |
| period 2, 60 samples, 300 reads, pcr | 3 | 878 | — | — | -0.382 ± 0.029 | -0.380 ± 0.029 | -0.641 ± 0.030 | — | — | -0.259 ± 0.017 | — | — | -0.261 ± 0.017 | — | — |

### Held-out prediction at one sample: the two fallbacks

Calling's fallback scores a class with no concentration of its own at a flat shape and the stated concentration; the other keeps the class's own shares at the stated concentration. D at one sample has no concentration in any class, so every class falls back; `D_working_F` is D scored at the concentration it worked at. The last column is the share of held-out loci whose class B left without a concentration. Mean over the seeds ± its standard error, a locus, against the truth; the paired column is seed by seed, above zero where D predicts better.

| cell | B_raw: calling / own | B_curves: calling / own | D: calling / own | D_working_F | D − B_curves, own shares | B classes on fallback |
|---|---|---|---|---|---|---|
| period 1, 1 samples, 10 reads, pcr | -1.087 ± 0.088 / -0.180 ± 0.013 | -1.092 ± 0.088 / -0.182 ± 0.015 | -1.650 ± 0.015 / -0.037 ± 0.007 | -0.036 ± 0.008 | +0.145 ± 0.014 | 0.60 |
| period 1, 1 samples, 3 reads, pcr | -1.537 ± 0.043 / -0.158 ± 0.011 | -1.524 ± 0.042 / -0.141 ± 0.012 | -1.628 ± 0.021 / -0.037 ± 0.004 | -0.038 ± 0.002 | +0.104 ± 0.011 | 0.91 |
| period 1, 1 samples, 3 reads, pcr_free | -1.510 ± 0.109 / -0.163 ± 0.002 | -1.489 ± 0.113 / -0.144 ± 0.006 | -1.657 ± 0.047 / -0.035 ± 0.002 | -0.035 ± 0.001 | +0.108 ± 0.005 | 0.89 |
| period 1, 1 samples, 30 reads, pcr | -0.887 ± 0.083 / -0.160 ± 0.012 | -0.868 ± 0.082 / -0.139 ± 0.011 | -1.695 ± 0.020 / -0.031 ± 0.003 | -0.030 ± 0.003 | +0.107 ± 0.011 | 0.43 |
| period 1, 1 samples, 300 reads, pcr | -0.954 ± 0.059 / -0.169 ± 0.006 | -0.942 ± 0.058 / -0.157 ± 0.006 | -1.653 ± 0.027 / -0.041 ± 0.002 | -0.040 ± 0.002 | +0.116 ± 0.007 | 0.49 |
| period 1, 1 samples, 300 reads, pcr_free | -1.125 ± 0.145 / -0.301 ± 0.148 | -0.975 ± 0.125 / -0.148 ± 0.001 | -1.688 ± 0.036 / -0.037 ± 0.002 | -0.036 ± 0.002 | +0.111 ± 0.001 | 0.50 |
| period 2, 1 samples, 10 reads, pcr | -1.614 ± 0.101 / -0.220 ± 0.016 | -1.594 ± 0.101 / -0.196 ± 0.017 | -1.903 ± 0.021 / -0.040 ± 0.005 | -0.042 ± 0.005 | +0.157 ± 0.013 | 0.81 |
| period 2, 1 samples, 3 reads, pcr | -1.854 ± 0.025 / -0.257 ± 0.013 | -1.802 ± 0.023 / -0.186 ± 0.010 | -1.828 ± 0.020 / -0.047 ± 0.003 | -0.048 ± 0.003 | +0.138 ± 0.007 | 0.96 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | -1.822 ± 0.059 / -0.232 ± 0.018 | -1.778 ± 0.054 / -0.169 ± 0.006 | -1.863 ± 0.020 / -0.046 ± 0.002 | -0.049 ± 0.005 | +0.123 ± 0.008 | 0.95 |
| period 2, 1 samples, 3 reads, pcr_free | -1.980 ± 0.381 / -0.492 ± 0.279 | -1.697 ± 0.110 / -0.209 ± 0.009 | -1.880 ± 0.045 / -0.042 ± 0.007 | -0.045 ± 0.009 | +0.168 ± 0.002 | 0.90 |
| period 2, 1 samples, 30 reads, pcr | -1.447 ± 0.114 / -0.229 ± 0.013 | -1.412 ± 0.120 / -0.192 ± 0.006 | -1.888 ± 0.026 / -0.046 ± 0.003 | -0.046 ± 0.003 | +0.147 ± 0.009 | 0.72 |
| period 2, 1 samples, 300 reads, pcr | -1.269 ± 0.115 / -0.213 ± 0.011 | -1.252 ± 0.116 / -0.196 ± 0.011 | -1.885 ± 0.013 / -0.049 ± 0.003 | -0.051 ± 0.003 | +0.147 ± 0.012 | 0.65 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | -1.231 ± 0.043 / -0.195 ± 0.005 | -1.219 ± 0.041 / -0.183 ± 0.003 | -1.912 ± 0.028 / -0.050 ± 0.003 | -0.053 ± 0.003 | +0.132 ± 0.003 | 0.63 |
| period 2, 1 samples, 300 reads, pcr_free | -1.605 ± 0.356 / -0.512 ± 0.233 | -1.273 ± 0.149 / -0.176 ± 0.022 | -1.856 ± 0.015 / -0.039 ± 0.002 | -0.041 ± 0.002 | +0.137 ± 0.020 | 0.66 |

### Recovery, side by side

Median over the classes of each class's root mean square error over the seeds: relative for the slip rate (`level`), absolute otherwise. A dash: the method gave no estimate (F at one sample is absent by design for B and D).

| cell | method | level | shorter_share | fall_off | between_locus_share | share_0 |
|---|---|---|---|---|---|---|
| period 1, 1 samples, 10 reads, pcr | B_raw | 0.0993 | 0.0527 | 0.0368 | 0.2239 | 0.0535 |
| period 1, 1 samples, 10 reads, pcr | B_curves | 0.0743 | 0.0295 | 0.0291 | 0.2239 | 0.0535 |
| period 1, 1 samples, 10 reads, pcr | D | 0.0617 | 0.0198 | 0.0224 | — | 0.0337 |
| period 1, 1 samples, 10 reads, pcr | D_working_F | 0.0617 | 0.0198 | 0.0224 | 0.1913 | 0.0337 |
| period 1, 1 samples, 3 reads, pcr | A_raw | 0.1139 | 0.0976 | 0.0530 | 0.6574 | 0.0700 |
| period 1, 1 samples, 3 reads, pcr | A_curves | 0.2125 | 0.1007 | 0.0234 | 0.6574 | 0.0700 |
| period 1, 1 samples, 3 reads, pcr | B_raw | 0.2175 | 0.1586 | 0.0962 | 0.2415 | 0.0602 |
| period 1, 1 samples, 3 reads, pcr | B_curves | 0.1757 | 0.0943 | 0.0602 | 0.2415 | 0.0602 |
| period 1, 1 samples, 3 reads, pcr | D | 0.0873 | 0.0752 | 0.0450 | — | 0.0420 |
| period 1, 1 samples, 3 reads, pcr | D_working_F | 0.0873 | 0.0752 | 0.0450 | 0.3459 | 0.0420 |
| period 1, 1 samples, 3 reads, pcr | D_tract_prior | 0.0961 | 0.0627 | 0.0383 | — | 0.0483 |
| period 1, 1 samples, 3 reads, pcr_free | B_raw | 0.6567 | 0.3692 | 0.1973 | 0.4731 | 0.0630 |
| period 1, 1 samples, 3 reads, pcr_free | B_curves | 0.6696 | 0.3186 | 0.1282 | 0.4731 | 0.0630 |
| period 1, 1 samples, 3 reads, pcr_free | D | 0.3189 | 0.1928 | 0.0951 | — | 0.0379 |
| period 1, 1 samples, 3 reads, pcr_free | D_working_F | 0.3189 | 0.1928 | 0.0951 | 0.4630 | 0.0379 |
| period 1, 1 samples, 30 reads, pcr | B_raw | 0.0599 | 0.0316 | 0.0196 | 0.1803 | 0.0567 |
| period 1, 1 samples, 30 reads, pcr | B_curves | 0.0617 | 0.0174 | 0.0112 | 0.1803 | 0.0567 |
| period 1, 1 samples, 30 reads, pcr | D | 0.0398 | 0.0089 | 0.0107 | — | 0.0444 |
| period 1, 1 samples, 30 reads, pcr | D_working_F | 0.0398 | 0.0089 | 0.0107 | 0.1620 | 0.0444 |
| period 1, 1 samples, 300 reads, pcr | A_raw | 0.0172 | 0.0100 | 0.0060 | 0.1596 | 0.0506 |
| period 1, 1 samples, 300 reads, pcr | A_curves | 0.0167 | 0.0058 | 0.0074 | 0.1596 | 0.0506 |
| period 1, 1 samples, 300 reads, pcr | B_raw | 0.0209 | 0.0082 | 0.0079 | 0.3068 | 0.0518 |
| period 1, 1 samples, 300 reads, pcr | B_curves | 0.0204 | 0.0065 | 0.0064 | 0.3068 | 0.0518 |
| period 1, 1 samples, 300 reads, pcr | D | 0.0195 | 0.0054 | 0.0060 | — | 0.0372 |
| period 1, 1 samples, 300 reads, pcr | D_working_F | 0.0195 | 0.0054 | 0.0060 | 0.1408 | 0.0372 |
| period 1, 1 samples, 300 reads, pcr_free | B_raw | 0.0753 | 0.0327 | 0.0224 | 0.3047 | 0.0583 |
| period 1, 1 samples, 300 reads, pcr_free | B_curves | 0.0766 | 0.0185 | 0.0114 | 0.3047 | 0.0583 |
| period 1, 1 samples, 300 reads, pcr_free | D | 0.0549 | 0.0195 | 0.0106 | — | 0.0336 |
| period 1, 1 samples, 300 reads, pcr_free | D_working_F | 0.0549 | 0.0195 | 0.0106 | 0.1489 | 0.0336 |
| period 1, 10 samples, 10 reads, pcr | B_raw | 0.0233 | 0.0139 | 0.0102 | 0.0523 | 0.0379 |
| period 1, 10 samples, 10 reads, pcr | B_curves | 0.0246 | 0.0093 | 0.0060 | 0.0523 | 0.0379 |
| period 1, 10 samples, 10 reads, pcr | D | 0.0210 | 0.0091 | 0.0058 | 0.0463 | 0.0307 |
| period 1, 10 samples, 3 reads, pcr | A_raw | 0.0447 | 0.0255 | 0.0238 | 0.0368 | 0.0936 |
| period 1, 10 samples, 3 reads, pcr | A_curves | 0.0552 | 0.0136 | 0.0133 | 0.0368 | 0.0936 |
| period 1, 10 samples, 3 reads, pcr | B_raw | 0.0728 | 0.0335 | 0.0293 | 0.0425 | 0.0398 |
| period 1, 10 samples, 3 reads, pcr | B_curves | 0.0651 | 0.0271 | 0.0229 | 0.0425 | 0.0398 |
| period 1, 10 samples, 3 reads, pcr | D | 0.0575 | 0.0256 | 0.0183 | 0.0314 | 0.0321 |
| period 1, 10 samples, 3 reads, pcr | D_tract_prior | 0.0561 | 0.0251 | 0.0180 | 0.0319 | 0.0328 |
| period 1, 10 samples, 30 reads, pcr | B_raw | 0.0192 | 0.0073 | 0.0100 | 0.0436 | 0.0398 |
| period 1, 10 samples, 30 reads, pcr | B_curves | 0.0217 | 0.0043 | 0.0069 | 0.0436 | 0.0398 |
| period 1, 10 samples, 30 reads, pcr | D | 0.0173 | 0.0044 | 0.0084 | 0.0352 | 0.0365 |
| period 1, 10 samples, 300 reads, pcr | B_raw | 0.0065 | 0.0028 | 0.0024 | 0.0479 | 0.0428 |
| period 1, 10 samples, 300 reads, pcr | B_curves | 0.0068 | 0.0020 | 0.0021 | 0.0479 | 0.0428 |
| period 1, 10 samples, 300 reads, pcr | D | 0.0063 | 0.0021 | 0.0020 | 0.0332 | 0.0411 |
| period 1, 2000 samples, 10 reads, pcr | B_raw | 0.0024 | 0.0019 | 0.0012 | 0.0255 | 0.0619 |
| period 1, 2000 samples, 10 reads, pcr | B_curves | 0.0024 | 0.0019 | 0.0014 | 0.0255 | 0.0619 |
| period 1, 2000 samples, 10 reads, pcr | D | 0.0019 | 0.0012 | 0.0012 | 0.0292 | 0.0567 |
| period 1, 2000 samples, 3 reads, pcr | B_raw | 0.0083 | 0.0042 | 0.0015 | 0.0246 | 0.0543 |
| period 1, 2000 samples, 3 reads, pcr | B_curves | 0.0083 | 0.0029 | 0.0020 | 0.0246 | 0.0543 |
| period 1, 2000 samples, 3 reads, pcr | D | 0.0063 | 0.0025 | 0.0015 | 0.0254 | 0.0552 |
| period 1, 2000 samples, 3 reads, pcr_free | B_raw | 0.0176 | 0.0064 | 0.0082 | 0.0265 | 0.0313 |
| period 1, 2000 samples, 3 reads, pcr_free | B_curves | 0.0167 | 0.0048 | 0.0076 | 0.0265 | 0.0313 |
| period 1, 2000 samples, 3 reads, pcr_free | D | 0.0104 | 0.0032 | 0.0061 | 0.0293 | 0.0293 |
| period 1, 500 samples, 10 reads, pcr | B_raw | 0.0061 | 0.0024 | 0.0025 | 0.0310 | 0.0447 |
| period 1, 500 samples, 10 reads, pcr | B_curves | 0.0061 | 0.0019 | 0.0026 | 0.0310 | 0.0447 |
| period 1, 500 samples, 10 reads, pcr | D | 0.0060 | 0.0013 | 0.0016 | 0.0327 | 0.0412 |
| period 1, 500 samples, 3 reads, pcr | B_raw | 0.0128 | 0.0066 | 0.0031 | 0.0340 | 0.0475 |
| period 1, 500 samples, 3 reads, pcr | B_curves | 0.0128 | 0.0048 | 0.0024 | 0.0340 | 0.0475 |
| period 1, 500 samples, 3 reads, pcr | D | 0.0118 | 0.0046 | 0.0020 | 0.0361 | 0.0501 |
| period 1, 500 samples, 3 reads, pcr | D_tract_prior | 0.0092 | 0.0048 | 0.0018 | 0.0330 | 0.0319 |
| period 1, 500 samples, 30 reads, pcr | B_raw | 0.0027 | 0.0019 | 0.0011 | 0.0517 | 0.0503 |
| period 1, 500 samples, 30 reads, pcr | B_curves | 0.0027 | 0.0016 | 0.0013 | 0.0517 | 0.0503 |
| period 1, 500 samples, 30 reads, pcr | D | 0.0027 | 0.0011 | 0.0011 | 0.0500 | 0.0581 |
| period 1, 500 samples, 300 reads, pcr | B_raw | 0.0020 | 0.0011 | 0.0010 | 0.0513 | 0.0388 |
| period 1, 500 samples, 300 reads, pcr | B_curves | 0.0020 | 0.0008 | 0.0010 | 0.0513 | 0.0388 |
| period 1, 500 samples, 300 reads, pcr | D | 0.0020 | 0.0011 | 0.0008 | 0.0534 | 0.0387 |
| period 1, 60 samples, 10 reads, pcr | B_raw | 0.0138 | 0.0063 | 0.0061 | 0.0389 | 0.0427 |
| period 1, 60 samples, 10 reads, pcr | B_curves | 0.0138 | 0.0036 | 0.0039 | 0.0389 | 0.0427 |
| period 1, 60 samples, 10 reads, pcr | D | 0.0132 | 0.0030 | 0.0036 | 0.0392 | 0.0390 |
| period 1, 60 samples, 10 reads, pcr | D_tract_prior | 0.0133 | 0.0031 | 0.0036 | 0.0375 | 0.0424 |
| period 1, 60 samples, 10 reads, pcr_free | B_raw | 0.0424 | 0.0208 | 0.0158 | 0.0381 | 0.0436 |
| period 1, 60 samples, 10 reads, pcr_free | B_curves | 0.0494 | 0.0125 | 0.0110 | 0.0381 | 0.0436 |
| period 1, 60 samples, 10 reads, pcr_free | D | 0.0400 | 0.0069 | 0.0099 | 0.0393 | 0.0392 |
| period 1, 60 samples, 3 reads, pcr | B_raw | 0.0271 | 0.0155 | 0.0107 | 0.0457 | 0.0524 |
| period 1, 60 samples, 3 reads, pcr | B_curves | 0.0185 | 0.0108 | 0.0072 | 0.0457 | 0.0524 |
| period 1, 60 samples, 3 reads, pcr | D | 0.0219 | 0.0095 | 0.0062 | 0.0494 | 0.0461 |
| period 1, 60 samples, 30 reads, pcr | B_raw | 0.0080 | 0.0040 | 0.0027 | 0.0406 | 0.0325 |
| period 1, 60 samples, 30 reads, pcr | B_curves | 0.0080 | 0.0029 | 0.0026 | 0.0406 | 0.0325 |
| period 1, 60 samples, 30 reads, pcr | D | 0.0076 | 0.0027 | 0.0022 | 0.0409 | 0.0362 |
| period 1, 60 samples, 300 reads, pcr | B_raw | 0.0033 | 0.0014 | 0.0015 | 0.0401 | 0.0408 |
| period 1, 60 samples, 300 reads, pcr | B_curves | 0.0033 | 0.0013 | 0.0018 | 0.0401 | 0.0408 |
| period 1, 60 samples, 300 reads, pcr | D | 0.0033 | 0.0010 | 0.0015 | 0.0402 | 0.0412 |
| period 2, 1 samples, 10 reads, pcr | B_raw | 0.1330 | 0.0584 | 0.0579 | 0.3641 | 0.0789 |
| period 2, 1 samples, 10 reads, pcr | B_curves | 0.0985 | 0.0347 | 0.0269 | 0.3641 | 0.0789 |
| period 2, 1 samples, 10 reads, pcr | D | 0.0687 | 0.0300 | 0.0225 | — | 0.0456 |
| period 2, 1 samples, 10 reads, pcr | D_working_F | 0.0687 | 0.0300 | 0.0225 | 0.2045 | 0.0456 |
| period 2, 1 samples, 3 reads, pcr | A_raw | 0.1918 | 0.1384 | 0.0698 | 0.6675 | 0.0585 |
| period 2, 1 samples, 3 reads, pcr | A_curves | 0.0946 | 0.0504 | 0.0250 | 0.6675 | 0.0585 |
| period 2, 1 samples, 3 reads, pcr | B_raw | 0.2187 | 0.1163 | 0.1157 | 0.4809 | 0.0873 |
| period 2, 1 samples, 3 reads, pcr | B_curves | 0.1527 | 0.0785 | 0.0422 | 0.4809 | 0.0873 |
| period 2, 1 samples, 3 reads, pcr | D | 0.1327 | 0.0538 | 0.0393 | — | 0.0485 |
| period 2, 1 samples, 3 reads, pcr | D_working_F | 0.1327 | 0.0538 | 0.0393 | 0.3802 | 0.0485 |
| period 2, 1 samples, 3 reads, pcr | D_tract_prior | 0.0934 | 0.0638 | 0.0284 | — | 0.0514 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_raw | 0.2716 | 0.1289 | 0.1030 | 0.0979 | 0.0895 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_curves | 0.1499 | 0.0759 | 0.0397 | 0.0979 | 0.0895 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D | 0.1029 | 0.0364 | 0.0343 | — | 0.0453 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D_working_F | 0.1029 | 0.0364 | 0.0343 | 0.4177 | 0.0453 |
| period 2, 1 samples, 3 reads, pcr_free | B_raw | 0.7556 | 0.3023 | 0.2497 | 0.5441 | 0.0714 |
| period 2, 1 samples, 3 reads, pcr_free | B_curves | 0.6303 | 0.2286 | 0.2542 | 0.5441 | 0.0714 |
| period 2, 1 samples, 3 reads, pcr_free | D | 0.3403 | 0.1690 | 0.1334 | — | 0.0435 |
| period 2, 1 samples, 3 reads, pcr_free | D_working_F | 0.3403 | 0.1690 | 0.1334 | 0.2226 | 0.0435 |
| period 2, 1 samples, 30 reads, pcr | B_raw | 0.0879 | 0.0331 | 0.0340 | 0.3493 | 0.0817 |
| period 2, 1 samples, 30 reads, pcr | B_curves | 0.0763 | 0.0166 | 0.0199 | 0.3493 | 0.0817 |
| period 2, 1 samples, 30 reads, pcr | D | 0.0687 | 0.0105 | 0.0116 | — | 0.0421 |
| period 2, 1 samples, 30 reads, pcr | D_working_F | 0.0687 | 0.0105 | 0.0116 | 0.1114 | 0.0421 |
| period 2, 1 samples, 300 reads, pcr | A_raw | 0.0124 | 0.0046 | 0.0075 | 0.1945 | 0.0509 |
| period 2, 1 samples, 300 reads, pcr | A_curves | 0.0162 | 0.0046 | 0.0036 | 0.1945 | 0.0509 |
| period 2, 1 samples, 300 reads, pcr | B_raw | 0.0223 | 0.0114 | 0.0100 | 0.3132 | 0.0834 |
| period 2, 1 samples, 300 reads, pcr | B_curves | 0.0220 | 0.0068 | 0.0056 | 0.3132 | 0.0834 |
| period 2, 1 samples, 300 reads, pcr | D | 0.0191 | 0.0056 | 0.0036 | — | 0.0426 |
| period 2, 1 samples, 300 reads, pcr | D_working_F | 0.0191 | 0.0056 | 0.0036 | 0.1071 | 0.0426 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_raw | 0.0217 | 0.0128 | 0.0108 | 0.2238 | 0.0699 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_curves | 0.0203 | 0.0066 | 0.0045 | 0.2238 | 0.0699 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D | 0.0196 | 0.0046 | 0.0030 | — | 0.0430 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D_working_F | 0.0196 | 0.0046 | 0.0030 | 0.1517 | 0.0430 |
| period 2, 1 samples, 300 reads, pcr_free | B_raw | 0.0770 | 0.0390 | 0.0274 | 0.3283 | 0.0542 |
| period 2, 1 samples, 300 reads, pcr_free | B_curves | 0.0664 | 0.0156 | 0.0186 | 0.3283 | 0.0542 |
| period 2, 1 samples, 300 reads, pcr_free | D | 0.0518 | 0.0116 | 0.0140 | — | 0.0399 |
| period 2, 1 samples, 300 reads, pcr_free | D_working_F | 0.0518 | 0.0116 | 0.0140 | 0.1193 | 0.0399 |
| period 2, 10 samples, 10 reads, pcr | B_raw | 0.0338 | 0.0173 | 0.0151 | 0.0590 | 0.0897 |
| period 2, 10 samples, 10 reads, pcr | B_curves | 0.0301 | 0.0115 | 0.0072 | 0.0590 | 0.0897 |
| period 2, 10 samples, 10 reads, pcr | D | 0.0303 | 0.0051 | 0.0074 | 0.0219 | 0.0654 |
| period 2, 10 samples, 3 reads, pcr | B_raw | 0.0828 | 0.0461 | 0.0360 | 0.0592 | 0.0643 |
| period 2, 10 samples, 3 reads, pcr | B_curves | 0.0782 | 0.0267 | 0.0159 | 0.0592 | 0.0643 |
| period 2, 10 samples, 3 reads, pcr | D | 0.0449 | 0.0253 | 0.0110 | 0.0224 | 0.0522 |
| period 2, 10 samples, 3 reads, pcr | D_tract_prior | 0.0469 | 0.0254 | 0.0103 | 0.0205 | 0.0544 |
| period 2, 10 samples, 30 reads, pcr | B_raw | 0.0252 | 0.0116 | 0.0111 | 0.0711 | 0.0824 |
| period 2, 10 samples, 30 reads, pcr | B_curves | 0.0219 | 0.0058 | 0.0051 | 0.0711 | 0.0824 |
| period 2, 10 samples, 30 reads, pcr | D | 0.0184 | 0.0056 | 0.0055 | 0.0467 | 0.0582 |
| period 2, 10 samples, 300 reads, pcr | B_raw | 0.0070 | 0.0036 | 0.0041 | 0.0522 | 0.0628 |
| period 2, 10 samples, 300 reads, pcr | B_curves | 0.0071 | 0.0030 | 0.0025 | 0.0522 | 0.0628 |
| period 2, 10 samples, 300 reads, pcr | D | 0.0078 | 0.0016 | 0.0021 | 0.0299 | 0.0601 |
| period 2, 2000 samples, 10 reads, pcr | B_raw | 0.0035 | 0.0013 | 0.0012 | 0.0596 | 0.0996 |
| period 2, 2000 samples, 10 reads, pcr | B_curves | 0.0035 | 0.0013 | 0.0013 | 0.0596 | 0.0996 |
| period 2, 2000 samples, 10 reads, pcr | D | 0.0030 | 0.0012 | 0.0012 | 0.0608 | 0.1034 |
| period 2, 2000 samples, 3 reads, pcr | B_raw | 0.0089 | 0.0040 | 0.0033 | 0.0539 | 0.0673 |
| period 2, 2000 samples, 3 reads, pcr | B_curves | 0.0087 | 0.0038 | 0.0030 | 0.0539 | 0.0673 |
| period 2, 2000 samples, 3 reads, pcr | D | 0.0042 | 0.0019 | 0.0018 | 0.0554 | 0.0767 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_raw | 0.0065 | 0.0018 | 0.0030 | 0.0245 | 0.0342 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_curves | 0.0062 | 0.0020 | 0.0026 | 0.0245 | 0.0342 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | D | 0.0029 | 0.0010 | 0.0019 | 0.0351 | 0.0330 |
| period 2, 2000 samples, 3 reads, pcr_free | B_raw | 0.0164 | 0.0038 | 0.0109 | 0.0536 | 0.0414 |
| period 2, 2000 samples, 3 reads, pcr_free | B_curves | 0.0140 | 0.0047 | 0.0076 | 0.0536 | 0.0414 |
| period 2, 2000 samples, 3 reads, pcr_free | D | 0.0110 | 0.0048 | 0.0081 | 0.0511 | 0.0526 |
| period 2, 500 samples, 10 reads, pcr | B_raw | 0.0063 | 0.0028 | 0.0021 | 0.0451 | 0.0696 |
| period 2, 500 samples, 10 reads, pcr | B_curves | 0.0062 | 0.0026 | 0.0016 | 0.0451 | 0.0696 |
| period 2, 500 samples, 10 reads, pcr | D | 0.0058 | 0.0013 | 0.0018 | 0.0528 | 0.0756 |
| period 2, 500 samples, 3 reads, pcr | B_raw | 0.0128 | 0.0049 | 0.0055 | 0.0455 | 0.0537 |
| period 2, 500 samples, 3 reads, pcr | B_curves | 0.0131 | 0.0044 | 0.0026 | 0.0455 | 0.0537 |
| period 2, 500 samples, 3 reads, pcr | D | 0.0097 | 0.0028 | 0.0020 | 0.0509 | 0.0508 |
| period 2, 500 samples, 3 reads, pcr | D_tract_prior | 0.0113 | 0.0011 | 0.0033 | 0.0473 | 0.0330 |
| period 2, 500 samples, 30 reads, pcr | B_raw | 0.0036 | 0.0015 | 0.0017 | 0.0477 | 0.0760 |
| period 2, 500 samples, 30 reads, pcr | B_curves | 0.0036 | 0.0015 | 0.0018 | 0.0477 | 0.0760 |
| period 2, 500 samples, 30 reads, pcr | D | 0.0035 | 0.0006 | 0.0013 | 0.0524 | 0.0804 |
| period 2, 500 samples, 300 reads, pcr | B_raw | 0.0015 | 0.0008 | 0.0008 | 0.0534 | 0.0587 |
| period 2, 500 samples, 300 reads, pcr | B_curves | 0.0015 | 0.0008 | 0.0007 | 0.0534 | 0.0587 |
| period 2, 500 samples, 300 reads, pcr | D | 0.0016 | 0.0007 | 0.0007 | 0.0574 | 0.0654 |
| period 2, 60 samples, 10 reads, pcr | B_raw | 0.0197 | 0.0073 | 0.0074 | 0.0542 | 0.0653 |
| period 2, 60 samples, 10 reads, pcr | B_curves | 0.0202 | 0.0040 | 0.0049 | 0.0542 | 0.0653 |
| period 2, 60 samples, 10 reads, pcr | D | 0.0132 | 0.0020 | 0.0039 | 0.0527 | 0.0616 |
| period 2, 60 samples, 10 reads, pcr | D_tract_prior | 0.0131 | 0.0021 | 0.0038 | 0.0524 | 0.0630 |
| period 2, 60 samples, 10 reads, pcr_F0.4 | B_raw | 0.0204 | 0.0078 | 0.0075 | 0.0495 | 0.0411 |
| period 2, 60 samples, 10 reads, pcr_F0.4 | B_curves | 0.0205 | 0.0039 | 0.0043 | 0.0495 | 0.0411 |
| period 2, 60 samples, 10 reads, pcr_F0.4 | D | 0.0157 | 0.0032 | 0.0037 | 0.0465 | 0.0401 |
| period 2, 60 samples, 10 reads, pcr_free | B_raw | 0.0487 | 0.0202 | 0.0233 | 0.0543 | 0.0665 |
| period 2, 60 samples, 10 reads, pcr_free | B_curves | 0.0406 | 0.0123 | 0.0119 | 0.0543 | 0.0665 |
| period 2, 60 samples, 10 reads, pcr_free | D | 0.0388 | 0.0084 | 0.0119 | 0.0522 | 0.0630 |
| period 2, 60 samples, 3 reads, pcr | B_raw | 0.0382 | 0.0156 | 0.0129 | 0.0710 | 0.0696 |
| period 2, 60 samples, 3 reads, pcr | B_curves | 0.0332 | 0.0093 | 0.0066 | 0.0710 | 0.0696 |
| period 2, 60 samples, 3 reads, pcr | D | 0.0248 | 0.0073 | 0.0048 | 0.0633 | 0.0734 |
| period 2, 60 samples, 30 reads, pcr | B_raw | 0.0106 | 0.0037 | 0.0046 | 0.0488 | 0.0550 |
| period 2, 60 samples, 30 reads, pcr | B_curves | 0.0104 | 0.0028 | 0.0029 | 0.0488 | 0.0550 |
| period 2, 60 samples, 30 reads, pcr | D | 0.0103 | 0.0017 | 0.0027 | 0.0450 | 0.0602 |
| period 2, 60 samples, 300 reads, pcr | B_raw | 0.0030 | 0.0015 | 0.0017 | 0.0460 | 0.0636 |
| period 2, 60 samples, 300 reads, pcr | B_curves | 0.0030 | 0.0015 | 0.0012 | 0.0460 | 0.0636 |
| period 2, 60 samples, 300 reads, pcr | D | 0.0031 | 0.0010 | 0.0014 | 0.0434 | 0.0583 |

### Error honesty, side by side

Share of estimates within one / two of their reported errors of the truth, pooled over classes and seeds (honest: about 0.68 / 0.95); in brackets, how many estimates had an error. Only the fits as they come carry errors (the curves replace them).

| cell | method | level | shorter_share | fall_off | between_locus_share | shares ≥ 0.05 | shares 0.01–0.05 | shares < 0.01 |
|---|---|---|---|---|---|---|---|---|
| period 1, 1 samples, 10 reads, pcr | B_raw | 0.68 / 0.94 (84) | 0.64 / 0.95 (83) | 0.67 / 0.95 (79) | 0.95 / 1.00 (22) | 0.68 / 0.93 (294) | 0.50 / 0.60 (402) | 0.23 / 0.33 (900) |
| period 1, 1 samples, 10 reads, pcr | B_curves | — | — | — | 0.95 / 1.00 (22) | 0.68 / 0.93 (294) | 0.50 / 0.60 (402) | 0.23 / 0.33 (900) |
| period 1, 1 samples, 10 reads, pcr | D | 0.73 / 0.98 (84) | 0.76 / 0.94 (84) | 0.36 / 0.85 (84) | — | 0.65 / 0.94 (294) | 0.69 / 0.89 (402) | 0.59 / 0.90 (900) |
| period 1, 1 samples, 10 reads, pcr | D_working_F | 0.73 / 0.98 (84) | 0.76 / 0.94 (84) | 0.36 / 0.85 (84) | — | 0.65 / 0.94 (294) | 0.69 / 0.89 (402) | 0.59 / 0.90 (900) |
| period 1, 1 samples, 3 reads, pcr | A_raw | 0.64 / 0.93 (14) | 0.58 / 0.92 (12) | 0.85 / 0.85 (13) | 0.00 / 0.50 (2) | 0.88 / 0.95 (41) | 0.79 / 0.97 (38) | 0.87 / 1.00 (23) |
| period 1, 1 samples, 3 reads, pcr | A_curves | — | — | — | 0.00 / 0.50 (2) | 0.88 / 0.95 (41) | 0.79 / 0.97 (38) | 0.87 / 1.00 (23) |
| period 1, 1 samples, 3 reads, pcr | B_raw | 0.61 / 0.92 (84) | 0.66 / 0.87 (77) | 0.54 / 0.86 (71) | 1.00 / 1.00 (3) | 0.67 / 0.89 (294) | 0.43 / 0.53 (402) | 0.26 / 0.35 (900) |
| period 1, 1 samples, 3 reads, pcr | B_curves | — | — | — | 1.00 / 1.00 (3) | 0.67 / 0.89 (294) | 0.43 / 0.53 (402) | 0.26 / 0.35 (900) |
| period 1, 1 samples, 3 reads, pcr | D | 0.69 / 0.98 (84) | 0.56 / 0.90 (84) | 0.58 / 0.93 (84) | — | 0.67 / 0.91 (294) | 0.70 / 0.89 (402) | 0.57 / 0.88 (900) |
| period 1, 1 samples, 3 reads, pcr | D_working_F | 0.69 / 0.98 (84) | 0.56 / 0.90 (84) | 0.58 / 0.93 (84) | — | 0.67 / 0.91 (294) | 0.70 / 0.89 (402) | 0.57 / 0.88 (900) |
| period 1, 1 samples, 3 reads, pcr | D_tract_prior | 0.69 / 0.95 (42) | 0.60 / 0.95 (42) | 0.69 / 0.93 (42) | — | 0.63 / 0.86 (147) | 0.75 / 0.95 (201) | 0.51 / 0.88 (450) |
| period 1, 1 samples, 3 reads, pcr_free | B_raw | 0.81 / 1.00 (31) | 0.68 / 0.95 (19) | 0.41 / 0.82 (17) | 1.00 / 1.00 (3) | 0.69 / 0.91 (143) | 0.47 / 0.57 (194) | 0.26 / 0.36 (442) |
| period 1, 1 samples, 3 reads, pcr_free | B_curves | — | — | — | 1.00 / 1.00 (3) | 0.69 / 0.91 (143) | 0.47 / 0.57 (194) | 0.26 / 0.36 (442) |
| period 1, 1 samples, 3 reads, pcr_free | D | 0.62 / 0.95 (42) | 0.55 / 0.88 (42) | 0.43 / 0.79 (42) | — | 0.72 / 0.93 (147) | 0.78 / 0.94 (201) | 0.66 / 0.90 (450) |
| period 1, 1 samples, 3 reads, pcr_free | D_working_F | 0.62 / 0.95 (42) | 0.55 / 0.88 (42) | 0.43 / 0.79 (42) | — | 0.72 / 0.93 (147) | 0.78 / 0.94 (201) | 0.66 / 0.90 (450) |
| period 1, 1 samples, 30 reads, pcr | B_raw | 0.61 / 0.95 (84) | 0.63 / 0.93 (84) | 0.76 / 0.98 (84) | 0.98 / 1.00 (43) | 0.64 / 0.91 (294) | 0.46 / 0.56 (402) | 0.23 / 0.32 (900) |
| period 1, 1 samples, 30 reads, pcr | B_curves | — | — | — | 0.98 / 1.00 (43) | 0.64 / 0.91 (294) | 0.46 / 0.56 (402) | 0.23 / 0.32 (900) |
| period 1, 1 samples, 30 reads, pcr | D | 0.62 / 0.98 (84) | 0.82 / 1.00 (84) | 0.67 / 0.94 (84) | — | 0.62 / 0.86 (294) | 0.59 / 0.86 (402) | 0.60 / 0.87 (900) |
| period 1, 1 samples, 30 reads, pcr | D_working_F | 0.62 / 0.98 (84) | 0.82 / 1.00 (84) | 0.67 / 0.94 (84) | — | 0.62 / 0.86 (294) | 0.59 / 0.86 (402) | 0.60 / 0.87 (900) |
| period 1, 1 samples, 300 reads, pcr | A_raw | 0.64 / 0.93 (14) | 0.50 / 1.00 (14) | 0.64 / 0.93 (14) | 0.75 / 1.00 (4) | 0.71 / 0.93 (45) | 0.83 / 0.94 (36) | 0.50 / 0.79 (28) |
| period 1, 1 samples, 300 reads, pcr | A_curves | — | — | — | 0.75 / 1.00 (4) | 0.71 / 0.93 (45) | 0.83 / 0.94 (36) | 0.50 / 0.79 (28) |
| period 1, 1 samples, 300 reads, pcr | B_raw | 0.63 / 0.90 (84) | 0.70 / 0.98 (84) | 0.62 / 0.93 (84) | 0.95 / 1.00 (37) | 0.65 / 0.93 (294) | 0.51 / 0.60 (402) | 0.24 / 0.32 (900) |
| period 1, 1 samples, 300 reads, pcr | B_curves | — | — | — | 0.95 / 1.00 (37) | 0.65 / 0.93 (294) | 0.51 / 0.60 (402) | 0.24 / 0.32 (900) |
| period 1, 1 samples, 300 reads, pcr | D | 0.67 / 0.92 (84) | 0.57 / 0.85 (84) | 0.51 / 0.90 (84) | — | 0.69 / 0.89 (294) | 0.67 / 0.93 (402) | 0.63 / 0.90 (900) |
| period 1, 1 samples, 300 reads, pcr | D_working_F | 0.67 / 0.92 (84) | 0.57 / 0.85 (84) | 0.51 / 0.90 (84) | — | 0.69 / 0.89 (294) | 0.67 / 0.93 (402) | 0.63 / 0.90 (900) |
| period 1, 1 samples, 300 reads, pcr_free | B_raw | 0.57 / 0.90 (42) | 0.69 / 0.86 (42) | 0.54 / 0.98 (41) | 0.94 / 1.00 (18) | 0.60 / 0.92 (147) | 0.49 / 0.60 (201) | 0.21 / 0.30 (450) |
| period 1, 1 samples, 300 reads, pcr_free | B_curves | — | — | — | 0.94 / 1.00 (18) | 0.60 / 0.92 (147) | 0.49 / 0.60 (201) | 0.21 / 0.30 (450) |
| period 1, 1 samples, 300 reads, pcr_free | D | 0.55 / 0.90 (42) | 0.40 / 0.67 (42) | 0.76 / 1.00 (42) | — | 0.76 / 0.92 (147) | 0.71 / 0.94 (201) | 0.60 / 0.87 (450) |
| period 1, 1 samples, 300 reads, pcr_free | D_working_F | 0.55 / 0.90 (42) | 0.40 / 0.67 (42) | 0.76 / 1.00 (42) | — | 0.76 / 0.92 (147) | 0.71 / 0.94 (201) | 0.60 / 0.87 (450) |
| period 1, 10 samples, 10 reads, pcr | B_raw | 0.86 / 0.95 (42) | 0.69 / 0.98 (42) | 0.76 / 1.00 (42) | 0.69 / 0.83 (42) | 0.64 / 0.90 (147) | 0.51 / 0.69 (201) | 0.35 / 0.46 (450) |
| period 1, 10 samples, 10 reads, pcr | B_curves | — | — | — | 0.69 / 0.83 (42) | 0.64 / 0.90 (147) | 0.51 / 0.69 (201) | 0.35 / 0.46 (450) |
| period 1, 10 samples, 10 reads, pcr | D | 0.76 / 0.98 (42) | 0.57 / 0.83 (42) | 0.67 / 1.00 (42) | 0.69 / 0.83 (42) | 0.63 / 0.90 (147) | 0.51 / 0.69 (201) | 0.44 / 0.55 (450) |
| period 1, 10 samples, 3 reads, pcr | A_raw | 0.64 / 1.00 (14) | 0.71 / 0.93 (14) | 0.62 / 0.92 (13) | 0.40 / 1.00 (5) | 0.37 / 0.51 (41) | 0.32 / 0.52 (50) | 0.36 / 0.52 (56) |
| period 1, 10 samples, 3 reads, pcr | A_curves | — | — | — | 0.40 / 1.00 (5) | 0.37 / 0.51 (41) | 0.32 / 0.52 (50) | 0.36 / 0.52 (56) |
| period 1, 10 samples, 3 reads, pcr | B_raw | 0.55 / 0.98 (42) | 0.76 / 0.98 (42) | 0.62 / 0.88 (40) | 0.69 / 0.95 (42) | 0.66 / 0.93 (147) | 0.54 / 0.66 (201) | 0.38 / 0.52 (450) |
| period 1, 10 samples, 3 reads, pcr | B_curves | — | — | — | 0.69 / 0.95 (42) | 0.66 / 0.93 (147) | 0.54 / 0.66 (201) | 0.38 / 0.52 (450) |
| period 1, 10 samples, 3 reads, pcr | D | 0.67 / 0.95 (42) | 0.43 / 0.76 (42) | 0.45 / 0.69 (42) | 0.76 / 0.95 (42) | 0.67 / 0.92 (147) | 0.55 / 0.72 (201) | 0.44 / 0.58 (450) |
| period 1, 10 samples, 3 reads, pcr | D_tract_prior | 0.71 / 0.95 (42) | 0.45 / 0.76 (42) | 0.48 / 0.69 (42) | 0.71 / 0.95 (42) | 0.67 / 0.93 (147) | 0.63 / 0.79 (201) | 0.63 / 0.85 (450) |
| period 1, 10 samples, 30 reads, pcr | B_raw | 0.67 / 0.90 (42) | 0.76 / 0.98 (42) | 0.50 / 0.88 (42) | 0.64 / 0.95 (42) | 0.65 / 0.89 (147) | 0.50 / 0.68 (201) | 0.38 / 0.51 (450) |
| period 1, 10 samples, 30 reads, pcr | B_curves | — | — | — | 0.64 / 0.95 (42) | 0.65 / 0.89 (147) | 0.50 / 0.68 (201) | 0.38 / 0.51 (450) |
| period 1, 10 samples, 30 reads, pcr | D | 0.69 / 0.90 (42) | 0.60 / 0.83 (42) | 0.57 / 0.79 (42) | 0.57 / 0.79 (42) | 0.68 / 0.88 (147) | 0.54 / 0.70 (201) | 0.45 / 0.60 (450) |
| period 1, 10 samples, 300 reads, pcr | B_raw | 0.62 / 0.95 (42) | 0.76 / 0.93 (42) | 0.62 / 0.95 (42) | 0.67 / 0.93 (42) | 0.61 / 0.90 (147) | 0.54 / 0.69 (201) | 0.29 / 0.41 (450) |
| period 1, 10 samples, 300 reads, pcr | B_curves | — | — | — | 0.67 / 0.93 (42) | 0.61 / 0.90 (147) | 0.54 / 0.69 (201) | 0.29 / 0.41 (450) |
| period 1, 10 samples, 300 reads, pcr | D | 0.60 / 0.93 (42) | 0.67 / 0.95 (42) | 0.69 / 0.95 (42) | 0.64 / 0.95 (42) | 0.61 / 0.89 (147) | 0.55 / 0.72 (201) | 0.34 / 0.44 (450) |
| period 1, 2000 samples, 10 reads, pcr | B_raw | 0.61 / 0.96 (28) | 0.57 / 0.89 (28) | 0.50 / 0.86 (28) | 0.89 / 1.00 (28) | 0.67 / 0.89 (98) | 0.50 / 0.68 (134) | 0.32 / 0.45 (300) |
| period 1, 2000 samples, 10 reads, pcr | B_curves | — | — | — | 0.89 / 1.00 (28) | 0.67 / 0.89 (98) | 0.50 / 0.68 (134) | 0.32 / 0.45 (300) |
| period 1, 2000 samples, 10 reads, pcr | D | 0.71 / 0.96 (28) | 0.61 / 0.93 (28) | 0.64 / 0.86 (28) | 0.82 / 1.00 (28) | 0.68 / 0.89 (98) | 0.51 / 0.67 (134) | 0.25 / 0.31 (300) |
| period 1, 2000 samples, 3 reads, pcr | B_raw | 0.43 / 0.68 (28) | 0.57 / 0.82 (28) | 0.50 / 0.82 (28) | 0.68 / 1.00 (28) | 0.66 / 0.92 (98) | 0.57 / 0.67 (134) | 0.36 / 0.48 (300) |
| period 1, 2000 samples, 3 reads, pcr | B_curves | — | — | — | 0.68 / 1.00 (28) | 0.66 / 0.92 (98) | 0.57 / 0.67 (134) | 0.36 / 0.48 (300) |
| period 1, 2000 samples, 3 reads, pcr | D | 0.61 / 0.86 (28) | 0.61 / 0.86 (28) | 0.68 / 0.93 (28) | 0.71 / 1.00 (28) | 0.62 / 0.92 (98) | 0.55 / 0.68 (134) | 0.27 / 0.32 (300) |
| period 1, 2000 samples, 3 reads, pcr_free | B_raw | 0.57 / 0.79 (14) | 0.79 / 0.93 (14) | 0.57 / 0.71 (14) | 0.57 / 1.00 (14) | 0.61 / 0.90 (49) | 0.55 / 0.66 (67) | 0.31 / 0.47 (150) |
| period 1, 2000 samples, 3 reads, pcr_free | B_curves | — | — | — | 0.57 / 1.00 (14) | 0.61 / 0.90 (49) | 0.55 / 0.66 (67) | 0.31 / 0.47 (150) |
| period 1, 2000 samples, 3 reads, pcr_free | D | 0.64 / 0.93 (14) | 0.64 / 1.00 (14) | 0.57 / 0.86 (14) | 0.64 / 1.00 (14) | 0.57 / 0.90 (49) | 0.54 / 0.67 (67) | 0.26 / 0.31 (150) |
| period 1, 500 samples, 10 reads, pcr | B_raw | 0.50 / 0.79 (42) | 0.64 / 0.98 (42) | 0.60 / 0.98 (42) | 0.69 / 0.95 (42) | 0.65 / 0.90 (147) | 0.52 / 0.68 (201) | 0.34 / 0.50 (450) |
| period 1, 500 samples, 10 reads, pcr | B_curves | — | — | — | 0.69 / 0.95 (42) | 0.65 / 0.90 (147) | 0.52 / 0.68 (201) | 0.34 / 0.50 (450) |
| period 1, 500 samples, 10 reads, pcr | D | 0.50 / 0.81 (42) | 0.67 / 0.95 (42) | 0.57 / 0.95 (42) | 0.69 / 0.93 (42) | 0.67 / 0.90 (147) | 0.53 / 0.67 (201) | 0.28 / 0.35 (450) |
| period 1, 500 samples, 3 reads, pcr | B_raw | 0.38 / 0.83 (42) | 0.69 / 0.93 (42) | 0.64 / 0.93 (42) | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.64 (201) | 0.37 / 0.47 (450) |
| period 1, 500 samples, 3 reads, pcr | B_curves | — | — | — | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.64 (201) | 0.37 / 0.47 (450) |
| period 1, 500 samples, 3 reads, pcr | D | 0.48 / 0.88 (42) | 0.38 / 0.62 (42) | 0.83 / 0.95 (42) | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.64 (201) | 0.28 / 0.32 (450) |
| period 1, 500 samples, 3 reads, pcr | D_tract_prior | 0.29 / 1.00 (14) | 0.21 / 0.50 (14) | 0.86 / 1.00 (14) | 0.57 / 1.00 (14) | 0.65 / 0.88 (49) | 0.49 / 0.60 (67) | 0.38 / 0.49 (150) |
| period 1, 500 samples, 30 reads, pcr | B_raw | 0.67 / 0.98 (42) | 0.57 / 0.93 (42) | 0.67 / 0.95 (42) | 0.50 / 0.88 (42) | 0.54 / 0.85 (147) | 0.48 / 0.66 (201) | 0.34 / 0.46 (450) |
| period 1, 500 samples, 30 reads, pcr | B_curves | — | — | — | 0.50 / 0.88 (42) | 0.54 / 0.85 (147) | 0.48 / 0.66 (201) | 0.34 / 0.46 (450) |
| period 1, 500 samples, 30 reads, pcr | D | 0.71 / 0.98 (42) | 0.69 / 0.98 (42) | 0.62 / 0.95 (42) | 0.52 / 0.90 (42) | 0.52 / 0.87 (147) | 0.49 / 0.66 (201) | 0.27 / 0.35 (450) |
| period 1, 500 samples, 300 reads, pcr | B_raw | 0.55 / 0.98 (42) | 0.48 / 0.95 (42) | 0.40 / 0.98 (42) | 0.67 / 0.90 (42) | 0.64 / 0.88 (147) | 0.52 / 0.67 (201) | 0.33 / 0.43 (450) |
| period 1, 500 samples, 300 reads, pcr | B_curves | — | — | — | 0.67 / 0.90 (42) | 0.64 / 0.88 (147) | 0.52 / 0.67 (201) | 0.33 / 0.43 (450) |
| period 1, 500 samples, 300 reads, pcr | D | 0.55 / 0.98 (42) | 0.50 / 0.90 (42) | 0.40 / 0.98 (42) | 0.57 / 0.88 (42) | 0.62 / 0.88 (147) | 0.53 / 0.66 (201) | 0.28 / 0.33 (450) |
| period 1, 60 samples, 10 reads, pcr | B_raw | 0.71 / 0.88 (42) | 0.74 / 1.00 (42) | 0.64 / 0.95 (42) | 0.64 / 0.98 (42) | 0.65 / 0.93 (147) | 0.54 / 0.69 (201) | 0.35 / 0.48 (450) |
| period 1, 60 samples, 10 reads, pcr | B_curves | — | — | — | 0.64 / 0.98 (42) | 0.65 / 0.93 (147) | 0.54 / 0.69 (201) | 0.35 / 0.48 (450) |
| period 1, 60 samples, 10 reads, pcr | D | 0.74 / 0.88 (42) | 0.76 / 1.00 (42) | 0.71 / 0.95 (42) | 0.62 / 0.98 (42) | 0.66 / 0.93 (147) | 0.53 / 0.68 (201) | 0.29 / 0.36 (450) |
| period 1, 60 samples, 10 reads, pcr | D_tract_prior | 0.74 / 0.88 (42) | 0.76 / 1.00 (42) | 0.71 / 0.95 (42) | 0.60 / 0.98 (42) | 0.67 / 0.93 (147) | 0.54 / 0.69 (201) | 0.37 / 0.52 (450) |
| period 1, 60 samples, 10 reads, pcr_free | B_raw | 0.64 / 0.93 (42) | 0.64 / 0.95 (42) | 0.74 / 0.95 (42) | 0.67 / 0.98 (42) | 0.65 / 0.93 (147) | 0.54 / 0.69 (201) | 0.33 / 0.50 (450) |
| period 1, 60 samples, 10 reads, pcr_free | B_curves | — | — | — | 0.67 / 0.98 (42) | 0.65 / 0.93 (147) | 0.54 / 0.69 (201) | 0.33 / 0.50 (450) |
| period 1, 60 samples, 10 reads, pcr_free | D | 0.64 / 0.95 (42) | 0.86 / 0.95 (42) | 0.69 / 0.93 (42) | 0.62 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.69 (201) | 0.28 / 0.37 (450) |
| period 1, 60 samples, 3 reads, pcr | B_raw | 0.79 / 0.95 (42) | 0.74 / 0.90 (42) | 0.67 / 0.95 (42) | 0.69 / 0.95 (42) | 0.61 / 0.86 (147) | 0.53 / 0.70 (201) | 0.36 / 0.51 (450) |
| period 1, 60 samples, 3 reads, pcr | B_curves | — | — | — | 0.69 / 0.95 (42) | 0.61 / 0.86 (147) | 0.53 / 0.70 (201) | 0.36 / 0.51 (450) |
| period 1, 60 samples, 3 reads, pcr | D | 0.74 / 0.98 (42) | 0.48 / 0.74 (42) | 0.62 / 0.93 (42) | 0.67 / 0.88 (42) | 0.65 / 0.87 (147) | 0.54 / 0.70 (201) | 0.26 / 0.35 (450) |
| period 1, 60 samples, 30 reads, pcr | B_raw | 0.76 / 0.90 (42) | 0.74 / 0.98 (42) | 0.76 / 0.90 (42) | 0.74 / 0.93 (42) | 0.65 / 0.90 (147) | 0.47 / 0.66 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 30 reads, pcr | B_curves | — | — | — | 0.74 / 0.93 (42) | 0.65 / 0.90 (147) | 0.47 / 0.66 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 30 reads, pcr | D | 0.76 / 0.90 (42) | 0.67 / 0.93 (42) | 0.74 / 0.95 (42) | 0.62 / 0.93 (42) | 0.65 / 0.89 (147) | 0.46 / 0.66 (201) | 0.27 / 0.32 (450) |
| period 1, 60 samples, 300 reads, pcr | B_raw | 0.64 / 0.93 (42) | 0.60 / 0.95 (42) | 0.48 / 0.88 (42) | 0.69 / 0.93 (42) | 0.61 / 0.91 (147) | 0.49 / 0.67 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 300 reads, pcr | B_curves | — | — | — | 0.69 / 0.93 (42) | 0.61 / 0.91 (147) | 0.49 / 0.67 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 300 reads, pcr | D | 0.64 / 0.93 (42) | 0.64 / 0.95 (42) | 0.50 / 0.88 (42) | 0.55 / 0.93 (42) | 0.61 / 0.91 (147) | 0.50 / 0.66 (201) | 0.31 / 0.35 (450) |
| period 2, 1 samples, 10 reads, pcr | B_raw | 0.71 / 0.92 (126) | 0.65 / 0.93 (124) | 0.65 / 0.87 (124) | 0.91 / 1.00 (11) | 0.65 / 0.89 (354) | 0.36 / 0.42 (744) | 0.19 / 0.28 (1296) |
| period 2, 1 samples, 10 reads, pcr | B_curves | — | — | — | 0.91 / 1.00 (11) | 0.65 / 0.89 (354) | 0.36 / 0.42 (744) | 0.19 / 0.28 (1296) |
| period 2, 1 samples, 10 reads, pcr | D | 0.68 / 0.96 (126) | 0.47 / 0.90 (126) | 0.75 / 0.97 (126) | — | 0.55 / 0.87 (354) | 0.71 / 0.91 (744) | 0.60 / 0.91 (1296) |
| period 2, 1 samples, 10 reads, pcr | D_working_F | 0.68 / 0.96 (126) | 0.47 / 0.90 (126) | 0.75 / 0.97 (126) | — | 0.55 / 0.87 (354) | 0.71 / 0.91 (744) | 0.60 / 0.91 (1296) |
| period 2, 1 samples, 3 reads, pcr | A_raw | 0.76 / 0.90 (21) | 0.61 / 0.89 (18) | 0.72 / 0.94 (18) | 0.00 / 0.17 (6) | 0.69 / 0.94 (49) | 0.86 / 0.98 (49) | 0.85 / 1.00 (27) |
| period 2, 1 samples, 3 reads, pcr | A_curves | — | — | — | 0.00 / 0.17 (6) | 0.69 / 0.94 (49) | 0.86 / 0.98 (49) | 0.85 / 1.00 (27) |
| period 2, 1 samples, 3 reads, pcr | B_raw | 0.66 / 0.92 (125) | 0.65 / 0.90 (109) | 0.65 / 0.88 (108) | — | 0.65 / 0.88 (354) | 0.36 / 0.42 (744) | 0.19 / 0.30 (1296) |
| period 2, 1 samples, 3 reads, pcr | B_curves | — | — | — | — | 0.65 / 0.88 (354) | 0.36 / 0.42 (744) | 0.19 / 0.30 (1296) |
| period 2, 1 samples, 3 reads, pcr | D | 0.68 / 0.90 (126) | 0.49 / 0.83 (126) | 0.60 / 0.86 (126) | — | 0.66 / 0.90 (354) | 0.71 / 0.93 (744) | 0.56 / 0.88 (1296) |
| period 2, 1 samples, 3 reads, pcr | D_working_F | 0.68 / 0.90 (126) | 0.49 / 0.83 (126) | 0.60 / 0.86 (126) | — | 0.66 / 0.90 (354) | 0.71 / 0.93 (744) | 0.56 / 0.88 (1296) |
| period 2, 1 samples, 3 reads, pcr | D_tract_prior | 0.76 / 0.89 (63) | 0.33 / 0.84 (63) | 0.83 / 0.92 (63) | — | 0.64 / 0.85 (177) | 0.73 / 0.95 (372) | 0.42 / 0.87 (648) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_raw | 0.66 / 0.90 (62) | 0.67 / 0.89 (55) | 0.60 / 0.84 (58) | 1.00 / 1.00 (1) | 0.60 / 0.89 (176) | 0.30 / 0.38 (371) | 0.19 / 0.30 (631) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_curves | — | — | — | 1.00 / 1.00 (1) | 0.60 / 0.89 (176) | 0.30 / 0.38 (371) | 0.19 / 0.30 (631) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D | 0.78 / 0.95 (63) | 0.67 / 1.00 (63) | 0.60 / 0.87 (63) | — | 0.59 / 0.87 (177) | 0.69 / 0.89 (372) | 0.59 / 0.93 (648) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D_working_F | 0.78 / 0.95 (63) | 0.67 / 1.00 (63) | 0.60 / 0.87 (63) | — | 0.59 / 0.87 (177) | 0.69 / 0.89 (372) | 0.59 / 0.93 (648) |
| period 2, 1 samples, 3 reads, pcr_free | B_raw | 0.90 / 1.00 (41) | 0.86 / 1.00 (7) | 0.17 / 0.17 (18) | 1.00 / 1.00 (1) | 0.67 / 0.90 (150) | 0.36 / 0.40 (320) | 0.21 / 0.30 (537) |
| period 2, 1 samples, 3 reads, pcr_free | B_curves | — | — | — | 1.00 / 1.00 (1) | 0.67 / 0.90 (150) | 0.36 / 0.40 (320) | 0.21 / 0.30 (537) |
| period 2, 1 samples, 3 reads, pcr_free | D | 0.44 / 0.78 (63) | 0.46 / 0.70 (63) | 0.44 / 0.76 (63) | — | 0.68 / 0.90 (177) | 0.67 / 0.90 (372) | 0.52 / 0.88 (648) |
| period 2, 1 samples, 3 reads, pcr_free | D_working_F | 0.44 / 0.78 (63) | 0.46 / 0.70 (63) | 0.44 / 0.76 (63) | — | 0.68 / 0.90 (177) | 0.67 / 0.90 (372) | 0.52 / 0.88 (648) |
| period 2, 1 samples, 30 reads, pcr | B_raw | 0.59 / 0.91 (126) | 0.62 / 0.94 (126) | 0.62 / 0.92 (125) | 0.96 / 1.00 (26) | 0.71 / 0.88 (354) | 0.35 / 0.43 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 30 reads, pcr | B_curves | — | — | — | 0.96 / 1.00 (26) | 0.71 / 0.88 (354) | 0.35 / 0.43 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 30 reads, pcr | D | 0.60 / 0.87 (126) | 0.58 / 0.93 (126) | 0.75 / 0.98 (126) | — | 0.60 / 0.87 (354) | 0.68 / 0.91 (744) | 0.62 / 0.93 (1296) |
| period 2, 1 samples, 30 reads, pcr | D_working_F | 0.60 / 0.87 (126) | 0.58 / 0.93 (126) | 0.75 / 0.98 (126) | — | 0.60 / 0.87 (354) | 0.68 / 0.91 (744) | 0.62 / 0.93 (1296) |
| period 2, 1 samples, 300 reads, pcr | A_raw | 0.76 / 0.95 (21) | 0.71 / 0.95 (21) | 0.62 / 0.95 (21) | 1.00 / 1.00 (1) | 0.82 / 0.94 (51) | 0.85 / 1.00 (46) | 0.67 / 0.90 (30) |
| period 2, 1 samples, 300 reads, pcr | A_curves | — | — | — | 1.00 / 1.00 (1) | 0.82 / 0.94 (51) | 0.85 / 1.00 (46) | 0.67 / 0.90 (30) |
| period 2, 1 samples, 300 reads, pcr | B_raw | 0.71 / 0.97 (126) | 0.67 / 0.95 (126) | 0.69 / 0.94 (126) | 1.00 / 1.00 (29) | 0.67 / 0.91 (354) | 0.35 / 0.41 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 300 reads, pcr | B_curves | — | — | — | 1.00 / 1.00 (29) | 0.67 / 0.91 (354) | 0.35 / 0.41 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 300 reads, pcr | D | 0.64 / 0.95 (126) | 0.65 / 0.96 (126) | 0.77 / 0.98 (126) | — | 0.60 / 0.92 (354) | 0.64 / 0.88 (744) | 0.58 / 0.88 (1296) |
| period 2, 1 samples, 300 reads, pcr | D_working_F | 0.64 / 0.95 (126) | 0.65 / 0.96 (126) | 0.77 / 0.98 (126) | — | 0.60 / 0.92 (354) | 0.64 / 0.88 (744) | 0.58 / 0.88 (1296) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_raw | 0.68 / 0.95 (63) | 0.54 / 0.87 (63) | 0.59 / 0.95 (63) | 0.91 / 1.00 (11) | 0.67 / 0.93 (177) | 0.35 / 0.41 (372) | 0.17 / 0.26 (648) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_curves | — | — | — | 0.91 / 1.00 (11) | 0.67 / 0.93 (177) | 0.35 / 0.41 (372) | 0.17 / 0.26 (648) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D | 0.65 / 0.95 (63) | 0.71 / 0.92 (63) | 0.87 / 1.00 (63) | — | 0.51 / 0.82 (177) | 0.71 / 0.91 (372) | 0.57 / 0.89 (648) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D_working_F | 0.65 / 0.95 (63) | 0.71 / 0.92 (63) | 0.87 / 1.00 (63) | — | 0.51 / 0.82 (177) | 0.71 / 0.91 (372) | 0.57 / 0.89 (648) |
| period 2, 1 samples, 300 reads, pcr_free | B_raw | 0.70 / 0.92 (63) | 0.60 / 0.89 (63) | 0.59 / 0.90 (63) | 1.00 / 1.00 (14) | 0.68 / 0.91 (177) | 0.34 / 0.40 (372) | 0.18 / 0.27 (648) |
| period 2, 1 samples, 300 reads, pcr_free | B_curves | — | — | — | 1.00 / 1.00 (14) | 0.68 / 0.91 (177) | 0.34 / 0.40 (372) | 0.18 / 0.27 (648) |
| period 2, 1 samples, 300 reads, pcr_free | D | 0.54 / 0.90 (63) | 0.70 / 0.95 (63) | 0.67 / 0.90 (63) | — | 0.53 / 0.90 (177) | 0.66 / 0.87 (372) | 0.60 / 0.90 (648) |
| period 2, 1 samples, 300 reads, pcr_free | D_working_F | 0.54 / 0.90 (63) | 0.70 / 0.95 (63) | 0.67 / 0.90 (63) | — | 0.53 / 0.90 (177) | 0.66 / 0.87 (372) | 0.60 / 0.90 (648) |
| period 2, 10 samples, 10 reads, pcr | B_raw | 0.70 / 0.97 (63) | 0.70 / 0.97 (63) | 0.68 / 0.94 (63) | 0.67 / 0.90 (63) | 0.66 / 0.89 (177) | 0.44 / 0.55 (372) | 0.26 / 0.38 (648) |
| period 2, 10 samples, 10 reads, pcr | B_curves | — | — | — | 0.67 / 0.90 (63) | 0.66 / 0.89 (177) | 0.44 / 0.55 (372) | 0.26 / 0.38 (648) |
| period 2, 10 samples, 10 reads, pcr | D | 0.70 / 0.98 (63) | 0.84 / 0.97 (63) | 0.68 / 0.97 (63) | 0.75 / 0.89 (63) | 0.68 / 0.89 (177) | 0.46 / 0.59 (372) | 0.35 / 0.47 (648) |
| period 2, 10 samples, 3 reads, pcr | B_raw | 0.63 / 0.95 (63) | 0.65 / 0.90 (62) | 0.59 / 0.87 (61) | 0.78 / 0.95 (63) | 0.63 / 0.84 (177) | 0.45 / 0.54 (372) | 0.29 / 0.41 (648) |
| period 2, 10 samples, 3 reads, pcr | B_curves | — | — | — | 0.78 / 0.95 (63) | 0.63 / 0.84 (177) | 0.45 / 0.54 (372) | 0.29 / 0.41 (648) |
| period 2, 10 samples, 3 reads, pcr | D | 0.59 / 0.90 (63) | 0.49 / 0.76 (63) | 0.79 / 1.00 (63) | 0.78 / 0.94 (63) | 0.63 / 0.86 (177) | 0.50 / 0.63 (372) | 0.37 / 0.50 (648) |
| period 2, 10 samples, 3 reads, pcr | D_tract_prior | 0.59 / 0.89 (63) | 0.43 / 0.75 (63) | 0.83 / 1.00 (63) | 0.59 / 1.00 (63) | 0.62 / 0.85 (177) | 0.60 / 0.75 (372) | 0.61 / 0.83 (648) |
| period 2, 10 samples, 30 reads, pcr | B_raw | 0.59 / 0.92 (63) | 0.67 / 0.92 (63) | 0.60 / 0.92 (63) | 0.56 / 0.92 (63) | 0.66 / 0.88 (177) | 0.41 / 0.52 (372) | 0.25 / 0.35 (648) |
| period 2, 10 samples, 30 reads, pcr | B_curves | — | — | — | 0.56 / 0.92 (63) | 0.66 / 0.88 (177) | 0.41 / 0.52 (372) | 0.25 / 0.35 (648) |
| period 2, 10 samples, 30 reads, pcr | D | 0.62 / 0.89 (63) | 0.60 / 0.98 (63) | 0.62 / 0.97 (63) | 0.43 / 0.78 (63) | 0.65 / 0.88 (177) | 0.42 / 0.56 (372) | 0.31 / 0.43 (648) |
| period 2, 10 samples, 300 reads, pcr | B_raw | 0.65 / 0.90 (63) | 0.67 / 0.90 (63) | 0.62 / 0.94 (63) | 0.71 / 0.94 (63) | 0.66 / 0.90 (177) | 0.45 / 0.54 (372) | 0.26 / 0.36 (648) |
| period 2, 10 samples, 300 reads, pcr | B_curves | — | — | — | 0.71 / 0.94 (63) | 0.66 / 0.90 (177) | 0.45 / 0.54 (372) | 0.26 / 0.36 (648) |
| period 2, 10 samples, 300 reads, pcr | D | 0.59 / 0.92 (63) | 0.68 / 0.98 (63) | 0.56 / 0.90 (63) | 0.83 / 0.98 (63) | 0.67 / 0.89 (177) | 0.46 / 0.56 (372) | 0.36 / 0.47 (648) |
| period 2, 2000 samples, 10 reads, pcr | B_raw | 0.60 / 0.95 (42) | 0.60 / 0.93 (42) | 0.67 / 0.93 (42) | 0.57 / 0.88 (42) | 0.56 / 0.78 (118) | 0.44 / 0.57 (248) | 0.30 / 0.42 (432) |
| period 2, 2000 samples, 10 reads, pcr | B_curves | — | — | — | 0.57 / 0.88 (42) | 0.56 / 0.78 (118) | 0.44 / 0.57 (248) | 0.30 / 0.42 (432) |
| period 2, 2000 samples, 10 reads, pcr | D | 0.74 / 1.00 (42) | 0.60 / 0.95 (42) | 0.55 / 0.95 (42) | 0.60 / 0.86 (42) | 0.55 / 0.79 (118) | 0.43 / 0.56 (248) | 0.22 / 0.25 (432) |
| period 2, 2000 samples, 3 reads, pcr | B_raw | 0.51 / 0.83 (63) | 0.62 / 0.92 (63) | 0.52 / 0.87 (63) | 0.56 / 0.89 (63) | 0.63 / 0.89 (177) | 0.44 / 0.57 (372) | 0.26 / 0.39 (648) |
| period 2, 2000 samples, 3 reads, pcr | B_curves | — | — | — | 0.56 / 0.89 (63) | 0.63 / 0.89 (177) | 0.44 / 0.57 (372) | 0.26 / 0.39 (648) |
| period 2, 2000 samples, 3 reads, pcr | D | 0.73 / 0.94 (63) | 0.67 / 0.90 (63) | 0.67 / 0.97 (63) | 0.52 / 0.86 (63) | 0.66 / 0.90 (177) | 0.41 / 0.55 (372) | 0.19 / 0.22 (648) |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_raw | 0.52 / 0.86 (21) | 0.67 / 1.00 (21) | 0.48 / 0.76 (21) | 0.62 / 0.90 (21) | 0.59 / 0.95 (59) | 0.47 / 0.61 (124) | 0.36 / 0.54 (216) |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_curves | — | — | — | 0.62 / 0.90 (21) | 0.59 / 0.95 (59) | 0.47 / 0.61 (124) | 0.36 / 0.54 (216) |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | D | 0.81 / 0.90 (21) | 0.76 / 1.00 (21) | 0.57 / 0.95 (21) | 0.62 / 0.76 (21) | 0.63 / 0.93 (59) | 0.45 / 0.60 (124) | 0.27 / 0.33 (216) |
| period 2, 2000 samples, 3 reads, pcr_free | B_raw | 0.81 / 0.95 (21) | 0.81 / 0.95 (21) | 0.38 / 0.90 (21) | 0.62 / 0.95 (21) | 0.66 / 0.88 (59) | 0.40 / 0.53 (124) | 0.31 / 0.45 (216) |
| period 2, 2000 samples, 3 reads, pcr_free | B_curves | — | — | — | 0.62 / 0.95 (21) | 0.66 / 0.88 (59) | 0.40 / 0.53 (124) | 0.31 / 0.45 (216) |
| period 2, 2000 samples, 3 reads, pcr_free | D | 0.71 / 0.95 (21) | 0.43 / 1.00 (21) | 0.29 / 0.81 (21) | 0.57 / 0.95 (21) | 0.66 / 0.90 (59) | 0.37 / 0.50 (124) | 0.21 / 0.25 (216) |
| period 2, 500 samples, 10 reads, pcr | B_raw | 0.73 / 0.98 (63) | 0.60 / 0.95 (63) | 0.65 / 0.95 (63) | 0.75 / 0.90 (63) | 0.61 / 0.87 (177) | 0.41 / 0.52 (372) | 0.26 / 0.38 (648) |
| period 2, 500 samples, 10 reads, pcr | B_curves | — | — | — | 0.75 / 0.90 (63) | 0.61 / 0.87 (177) | 0.41 / 0.52 (372) | 0.26 / 0.38 (648) |
| period 2, 500 samples, 10 reads, pcr | D | 0.67 / 0.97 (63) | 0.68 / 0.97 (63) | 0.67 / 0.92 (63) | 0.71 / 0.94 (63) | 0.60 / 0.88 (177) | 0.41 / 0.51 (372) | 0.18 / 0.21 (648) |
| period 2, 500 samples, 3 reads, pcr | B_raw | 0.56 / 0.89 (63) | 0.73 / 0.98 (63) | 0.63 / 1.00 (63) | 0.71 / 0.92 (63) | 0.64 / 0.93 (177) | 0.41 / 0.51 (372) | 0.27 / 0.40 (648) |
| period 2, 500 samples, 3 reads, pcr | B_curves | — | — | — | 0.71 / 0.92 (63) | 0.64 / 0.93 (177) | 0.41 / 0.51 (372) | 0.27 / 0.40 (648) |
| period 2, 500 samples, 3 reads, pcr | D | 0.67 / 0.97 (63) | 0.60 / 0.98 (63) | 0.79 / 0.97 (63) | 0.71 / 0.92 (63) | 0.62 / 0.93 (177) | 0.40 / 0.49 (372) | 0.18 / 0.22 (648) |
| period 2, 500 samples, 3 reads, pcr | D_tract_prior | 0.48 / 0.90 (21) | 0.95 / 1.00 (21) | 0.38 / 0.81 (21) | 0.52 / 0.90 (21) | 0.73 / 0.97 (59) | 0.33 / 0.44 (124) | 0.26 / 0.37 (216) |
| period 2, 500 samples, 30 reads, pcr | B_raw | 0.63 / 0.89 (63) | 0.57 / 0.94 (63) | 0.65 / 0.95 (63) | 0.67 / 0.90 (63) | 0.62 / 0.90 (177) | 0.38 / 0.50 (372) | 0.27 / 0.38 (648) |
| period 2, 500 samples, 30 reads, pcr | B_curves | — | — | — | 0.67 / 0.90 (63) | 0.62 / 0.90 (177) | 0.38 / 0.50 (372) | 0.27 / 0.38 (648) |
| period 2, 500 samples, 30 reads, pcr | D | 0.63 / 0.94 (63) | 0.81 / 0.95 (63) | 0.70 / 0.97 (63) | 0.67 / 0.89 (63) | 0.63 / 0.91 (177) | 0.37 / 0.49 (372) | 0.20 / 0.23 (648) |
| period 2, 500 samples, 300 reads, pcr | B_raw | 0.56 / 0.92 (63) | 0.57 / 0.94 (63) | 0.65 / 0.95 (63) | 0.52 / 0.92 (63) | 0.63 / 0.90 (177) | 0.41 / 0.53 (372) | 0.25 / 0.36 (648) |
| period 2, 500 samples, 300 reads, pcr | B_curves | — | — | — | 0.52 / 0.92 (63) | 0.63 / 0.90 (177) | 0.41 / 0.53 (372) | 0.25 / 0.36 (648) |
| period 2, 500 samples, 300 reads, pcr | D | 0.56 / 0.92 (63) | 0.65 / 0.97 (63) | 0.63 / 0.95 (63) | 0.51 / 0.89 (63) | 0.62 / 0.92 (177) | 0.41 / 0.52 (372) | 0.19 / 0.23 (648) |
| period 2, 60 samples, 10 reads, pcr | B_raw | 0.62 / 0.89 (63) | 0.71 / 0.92 (63) | 0.59 / 0.94 (63) | 0.68 / 0.95 (63) | 0.63 / 0.89 (177) | 0.43 / 0.52 (372) | 0.29 / 0.41 (648) |
| period 2, 60 samples, 10 reads, pcr | B_curves | — | — | — | 0.68 / 0.95 (63) | 0.63 / 0.89 (177) | 0.43 / 0.52 (372) | 0.29 / 0.41 (648) |
| period 2, 60 samples, 10 reads, pcr | D | 0.67 / 0.94 (63) | 0.90 / 0.97 (63) | 0.63 / 0.97 (63) | 0.60 / 0.95 (63) | 0.63 / 0.91 (177) | 0.41 / 0.51 (372) | 0.20 / 0.25 (648) |
| period 2, 60 samples, 10 reads, pcr | D_tract_prior | 0.68 / 0.90 (63) | 0.90 / 0.97 (63) | 0.65 / 0.97 (63) | 0.62 / 0.97 (63) | 0.66 / 0.90 (177) | 0.43 / 0.53 (372) | 0.35 / 0.51 (648) |
| period 2, 60 samples, 10 reads, pcr_F0.4 | B_raw | 0.57 / 0.86 (63) | 0.62 / 0.94 (63) | 0.57 / 0.87 (63) | 0.67 / 0.92 (63) | 0.64 / 0.92 (177) | 0.50 / 0.65 (372) | 0.37 / 0.53 (648) |
| period 2, 60 samples, 10 reads, pcr_F0.4 | B_curves | — | — | — | 0.67 / 0.92 (63) | 0.64 / 0.92 (177) | 0.50 / 0.65 (372) | 0.37 / 0.53 (648) |
| period 2, 60 samples, 10 reads, pcr_F0.4 | D | 0.59 / 0.87 (63) | 0.84 / 0.95 (63) | 0.63 / 0.90 (63) | 0.65 / 0.90 (63) | 0.62 / 0.94 (177) | 0.50 / 0.63 (372) | 0.29 / 0.40 (648) |
| period 2, 60 samples, 10 reads, pcr_free | B_raw | 0.76 / 0.92 (63) | 0.68 / 0.94 (63) | 0.59 / 0.95 (63) | 0.67 / 0.97 (63) | 0.63 / 0.89 (177) | 0.43 / 0.53 (372) | 0.28 / 0.42 (648) |
| period 2, 60 samples, 10 reads, pcr_free | B_curves | — | — | — | 0.67 / 0.97 (63) | 0.63 / 0.89 (177) | 0.43 / 0.53 (372) | 0.28 / 0.42 (648) |
| period 2, 60 samples, 10 reads, pcr_free | D | 0.67 / 0.92 (63) | 0.57 / 0.84 (63) | 0.60 / 0.98 (63) | 0.62 / 0.97 (63) | 0.63 / 0.91 (177) | 0.41 / 0.52 (372) | 0.20 / 0.25 (648) |
| period 2, 60 samples, 3 reads, pcr | B_raw | 0.67 / 0.87 (63) | 0.75 / 0.90 (63) | 0.65 / 0.95 (63) | 0.67 / 0.92 (63) | 0.68 / 0.88 (177) | 0.46 / 0.54 (372) | 0.26 / 0.40 (648) |
| period 2, 60 samples, 3 reads, pcr | B_curves | — | — | — | 0.67 / 0.92 (63) | 0.68 / 0.88 (177) | 0.46 / 0.54 (372) | 0.26 / 0.40 (648) |
| period 2, 60 samples, 3 reads, pcr | D | 0.54 / 0.83 (63) | 0.68 / 0.98 (63) | 0.78 / 0.98 (63) | 0.65 / 0.90 (63) | 0.67 / 0.90 (177) | 0.45 / 0.53 (372) | 0.18 / 0.24 (648) |
| period 2, 60 samples, 30 reads, pcr | B_raw | 0.62 / 0.90 (63) | 0.79 / 0.92 (63) | 0.63 / 0.94 (63) | 0.68 / 0.94 (63) | 0.59 / 0.91 (177) | 0.43 / 0.55 (372) | 0.26 / 0.37 (648) |
| period 2, 60 samples, 30 reads, pcr | B_curves | — | — | — | 0.68 / 0.94 (63) | 0.59 / 0.91 (177) | 0.43 / 0.55 (372) | 0.26 / 0.37 (648) |
| period 2, 60 samples, 30 reads, pcr | D | 0.57 / 0.92 (63) | 0.76 / 0.98 (63) | 0.59 / 0.89 (63) | 0.65 / 0.94 (63) | 0.65 / 0.90 (177) | 0.42 / 0.54 (372) | 0.19 / 0.23 (648) |
| period 2, 60 samples, 300 reads, pcr | B_raw | 0.67 / 0.95 (63) | 0.57 / 0.98 (63) | 0.57 / 0.90 (63) | 0.76 / 0.94 (63) | 0.58 / 0.88 (177) | 0.47 / 0.56 (372) | 0.26 / 0.36 (648) |
| period 2, 60 samples, 300 reads, pcr | B_curves | — | — | — | 0.76 / 0.94 (63) | 0.58 / 0.88 (177) | 0.47 / 0.56 (372) | 0.26 / 0.36 (648) |
| period 2, 60 samples, 300 reads, pcr | D | 0.68 / 0.94 (63) | 0.63 / 1.00 (63) | 0.57 / 0.97 (63) | 0.75 / 0.94 (63) | 0.57 / 0.88 (177) | 0.48 / 0.54 (372) | 0.20 / 0.25 (648) |

