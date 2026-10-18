# Outlier share, phase 2: can a fit recover the stray-read share without truth?

**Date:** 2026-10-18. **Plan:** [`str_outlier_share_plan.md`](str_outlier_share_plan.md) §4,
phases 2 and 3 (the estimator part). **Code:** `experiments/str_outlier_share/sim.py`; the
grid's output is `experiments/str_outlier_share/phase2_grid_seed1.txt`.

## What was asked

`estimate-parameters` never sees the truth. To fit the stray-read share it would have to separate
stray reads from stutter and from real second alleles using the reads alone, summing over each
sample's unknown genotype the way the stutter fit already does. This phase asks whether that is
possible across the range the caller has to cover — one sample to many, a few reads to hundreds —
and whether the noisy tracts phase 1 found (2–5 in 100, with a third or more of their reads at
another length) bias the answer.

## The simulation

One stratum (a homopolymer of about 15 repeats), drawn as follows:

- **Allele lengths** from a fitted homopolymer stratum's shares, 90% at the reference length.
- **Each sample's genotype** from a Dirichlet-multinomial with concentration 1 over those shares —
  the one-sample prior calling uses.
- **Reads** Poisson around the depth; each from one of the two copies at random. A read is a
  **stray** with the true share (uniform over offsets −12 to +12 repeats); otherwise it slips with
  probability 0.04, shorter two times in three, by k repeats with chance proportional to
  `0.23^(k−1)`.
- **Noisy tracts**, in half the runs: 3 tracts in 100 take 30% of their reads, in every sample,
  from one other length 1 to 4 repeats from the reference — phase 1's middle class.

**The fit** maximises, over the slip rate, its direction, its fall-off, the stray share and the 19
allele-length shares, the sum over (tract, sample) of the log of the genotype-prior-weighted
likelihood, where each read's probability is `(1 − share)·stutter + share/25`. No class for noisy
tracts is in the model. Each fit sees about 600,000 reads, the size of a stratum at the census cap,
except at one sample and 3x, where 2,000 tracts give 6,000.

## What was found

**The share comes back within about 10% from 10 samples at any depth, and from one sample at 30x and
above; one sample at 3x and 10x is loose.** Fitted against true share, one seed per cell:

| samples × depth | true 0.001 | true 0.003 | true 0.010 |
|---|---|---|---|
| 1 × 3 | 0.0002 | 0.0036 | 0.0079 |
| 1 × 10 | 0.0006 | 0.0026 | 0.0113 |
| 1 × 30 / 100 / 300 | 0.0009–0.0011 | 0.0027–0.0030 | 0.0097–0.0102 |
| 10 × 3 / 10 / 30 / 100 / 300 | 0.0010–0.0011 | 0.0030–0.0033 | 0.0097–0.0103 |
| 100 × 3 / 10 / 30 / 100 / 300 | 0.0009–0.0010 | 0.0029–0.0031 | 0.0098–0.0102 |

**Noisy tracts do not move the share.** With 3 tracts in 100 noisy, every cell from one sample at
10x upwards still lands within 7% of the truth (at most 0.0006 off a true 0.010). The fit explains a noisy tract's second
length as a heterozygote, which costs the share nothing. Where they leave a mark is the stutter: the
slip rate comes out 0.041–0.046 for a true 0.040 at one sample and low depth, and the fall-off 0.27
for a true 0.23 at 100 samples and 3x.

**One sample at low depth is loose**: at 3x, 0.0002 for a true 0.001, and 0.0010 for 0.003 with
noisy tracts; at 10x, 0.0006 for 0.001. At 3x that is 6,000 reads in this simulation; a real
genome at 3x carries far more tracts per stratum, so a fit pooled over strata (one share per
period class — whether the share varies with repeat count is not yet measured) would see much
more.
Whatever is fitted there needs a standard error and a fallback below some count, which the stutter
fit already has the machinery for.

## What this does not show

- **Calling.** Whether a fitted share near 0.003 calls as well as the stated 0.07 is phase 4. On
  HG002, the fitted model at 0.003 tied 0.03 at 5x and 10x and lost 27 tracts in 8,044 at 30x
  (phase 1 report); that is the question the program has not answered.
- **One seed per cell.** The spread between cells is the only estimate of noise here.
- **A simplified fit.** The prototype is a direct maximisation in Python with the one-sample prior,
  not the penalised fit `estimate-parameters` runs; the flank-deletion reads phase 1 found (0.5–2.4
  in 1,000) are not simulated.
- **Real data.** The truth-free fit has not yet been run on HG002's psps to see whether it lands on
  phase 1's truth-based 1–5 in 1,000.
