# Fit precision — how far the 256-point average over a tract's length frequencies is from the truth

**Date:** 2026-10-01. **Decided by the owner at checkpoint C:** measure this before step D (step 1 of the plan
proposed in [fit_precision_c3_stopped_2026-10-01.md](fit_precision_c3_stopped_2026-10-01.md)). **Branch:**
`fit-precision`. A measurement only: no production setting changed (`QUADRATURE_POINTS` and `dirichlet_points`
are as they were).

## 1. The question

At each repeat tract the repeat-tract fit does not know how the stratum's alleles are shared out among its
chromosomes. It averages the tract's likelihood over every way they could be, weighted by the stratum's length
spectrum and concentration (a Dirichlet distribution). That average is taken over a fixed set of 256 points: a
Halton sequence mapped onto the Dirichlet by stick-breaking (`dirichlet_points`). Step C2 found it about 3.7
log-likelihood units a tract too low at the truth on one stratum of thirteen classes. This measures how wrong it
is across the range: classes, concentration, depth and cohort size.

## 2. How it was measured

`the_average_over_a_tracts_frequencies_against_the_truth` (ignored test, `ssr_fit.rs`), 60 drawn tracts a cell, at
the numbers they were drawn with; logs `tmp/fit_precision/average/run1.log` and `run2_3reads.log`.

- **The reference** is the same average by importance sampling, with a standard error for each tract. The points
  are drawn from a mixture: a tenth from the stratum's own Dirichlet, the rest from Dirichlets moved all, half and
  a quarter of the way towards the allele copies the tract's reads say its samples carry. It used 20,000 points a
  tract, and 200,000 for the 3-read cells at 9 and 13 classes. A narrower first proposal gave a reference biased
  low at 3 reads, and was replaced. The reference now agrees with 65,536 fixed points to within 0.15 a tract in
  every 20-sample cell. The one cell where it is not firm is 63 samples × 3 reads, where 49 of 60 tracts have an
  error above 0.05.
- **The fixed points**: the fit's own `dirichlet_points` at 256 (production), 4,096 and 65,536 points; the table
  gives the mean, over the cell's tracts, of the average's log minus the reference's.
- **Placed by the reads**: the importance sampler itself at 256 points — a first look at the alternative of
  putting each tract's points where its reads say the frequencies are — with the spread of its error over tracts.

## 3. What was measured

The error a tract, in log-likelihood units; negative is too low.

| classes | concentration | samples × reads | reference, a tract (its largest error) | 256 points | 4,096 | 65,536 | placed by the reads, 256 (spread) |
|---|---|---|---|---|---|---|---|
| 3 | 0.1 | 20 × 3 | -15.2 (0.12) | -0.03 | -0.00 | -0.00 | -0.11 (0.62) |
| 3 | 0.5 | 20 × 3 | -24.7 (0.08) | +0.02 | +0.00 | +0.00 | -0.17 (0.42) |
| 3 | 2 | 20 × 3 | -38.0 (0.05) | +0.01 | +0.00 | +0.00 | -0.04 (0.19) |
| 3 | 10 | 20 × 3 | -45.1 (0.01) | +0.01 | +0.00 | -0.00 | -0.00 (0.04) |
| 5 | 0.1 | 20 × 3 | -18.4 (0.86) | -0.15 | +0.00 | -0.01 | -0.88 (1.90) |
| 5 | 0.5 | 20 × 3 | -35.3 (0.29) | -0.22 | -0.03 | -0.02 | -0.72 (1.05) |
| 5 | 2 | 20 × 3 | -50.6 (0.32) | -0.11 | -0.01 | -0.02 | -0.09 (0.45) |
| 5 | 10 | 20 × 3 | -63.7 (0.04) | -0.02 | -0.01 | -0.00 | -0.01 (0.06) |
| 9 | 0.1 | 20 × 3 | -24.8 (0.19) | -0.08 | -0.01 | -0.01 | -3.25 (4.39) |
| 9 | 0.5 | 20 × 3 | -36.6 (0.32) | -0.67 | -0.18 | +0.03 | -1.15 (1.43) |
| 9 | 2 | 20 × 3 | -61.3 (0.10) | -1.96 | -0.25 | -0.03 | -0.19 (0.41) |
| 9 | 10 | 20 × 3 | -75.8 (0.01) | -0.04 | -0.01 | -0.01 | +0.00 (0.09) |
| 13 | 0.1 | 20 × 3 | -22.4 (0.99) | -1.46 | -0.03 | +0.03 | -4.18 (6.40) |
| 13 | 0.5 | 20 × 3 | -38.4 (0.91) | -2.42 | -0.30 | +0.04 | -1.93 (2.06) |
| 13 | 2 | 20 × 3 | -66.4 (0.64) | -4.38 | -0.73 | -0.14 | -0.28 (0.70) |
| 13 | 10 | 20 × 3 | -78.8 (0.03) | -0.85 | -0.16 | +0.01 | -0.02 (0.19) |
| 13 | 0.5 | 4 × 3 | -11.1 (0.20) | -0.42 | -0.02 | -0.01 | -0.14 (0.42) |
| 13 | 0.5 | 63 × 3 | -108.8 (1.00) | -17.30 | -2.16 | -0.09 | -7.74 (5.05) |
| 3 | 0.1 | 20 × 30 | -121.1 (0.01) | +0.01 | -0.00 | +0.00 | +0.00 (0.03) |
| 3 | 0.5 | 20 × 30 | -158.0 (0.01) | -0.00 | -0.00 | -0.00 | +0.00 (0.04) |
| 3 | 2 | 20 × 30 | -217.7 (0.01) | +0.00 | +0.00 | +0.00 | -0.01 (0.04) |
| 3 | 10 | 20 × 30 | -267.2 (0.01) | -0.01 | -0.00 | -0.00 | +0.00 (0.03) |
| 5 | 0.1 | 20 × 30 | -145.7 (0.01) | -0.11 | -0.04 | -0.00 | +0.00 (0.04) |
| 5 | 0.5 | 20 × 30 | -221.2 (0.01) | -0.22 | -0.04 | -0.01 | -0.00 (0.04) |
| 5 | 2 | 20 × 30 | -273.9 (0.01) | -0.23 | +0.02 | +0.00 | -0.01 (0.07) |
| 5 | 10 | 20 × 30 | -337.2 (0.01) | +0.00 | +0.00 | +0.00 | -0.00 (0.05) |
| 9 | 0.1 | 20 × 30 | -171.7 (0.06) | -1.63 | +0.03 | +0.00 | +0.02 (0.08) |
| 9 | 0.5 | 20 × 30 | -236.3 (0.02) | -4.19 | -0.73 | -0.05 | -0.01 (0.05) |
| 9 | 2 | 20 × 30 | -319.5 (0.02) | -2.20 | -0.44 | +0.01 | -0.00 (0.08) |
| 9 | 10 | 20 × 30 | -393.2 (0.02) | -0.74 | -0.13 | -0.00 | -0.00 (0.10) |
| 13 | 0.1 | 20 × 30 | -189.4 (0.14) | -8.86 | -0.10 | -0.01 | -0.06 (0.48) |
| 13 | 0.5 | 20 × 30 | -227.4 (0.03) | -2.26 | -0.18 | -0.01 | +0.00 (0.04) |
| 13 | 2 | 20 × 30 | -348.8 (0.02) | -8.15 | -1.91 | -0.29 | -0.01 (0.07) |
| 13 | 10 | 20 × 30 | -399.6 (0.02) | -1.90 | -0.27 | -0.01 | +0.00 (0.09) |
| 13 | 0.5 | 4 × 30 | -51.6 (0.01) | -1.67 | -0.22 | +0.01 | -0.00 (0.04) |
| 13 | 0.5 | 63 × 30 | -771.9 (0.11) | -36.40 | -6.64 | -0.81 | +0.02 (0.20) |

- **At three classes the 256 points are right**: within 0.03 a tract in every cell. This is the regime the
  decision to use 256 points was checked at (`parameter_prepass_joint_fit.md` §4.2).
- **At nine and thirteen classes they are too low, by up to 8.9 a tract at 20 samples**. The worst cells are
  concentration 0.1 to 2. More points help slowly: 4,096 still leave up to 1.9 a tract at thirteen classes, and
  65,536 leave 0.29 at concentration 2 and 30 reads.
- **The error grows with the cohort**: at thirteen classes, concentration 0.5 and 30 reads, it is 1.7 a tract at
  4 samples, 2.3 at 20 and 36.4 at 63 (0.81 even at 65,536 points). Each sample multiplies a tract's likelihood
  by its own reads' term, so the product peaks more sharply where the tract's real frequencies are, and fewer of
  a fixed set of points land under that peak. A cohort of thousands is far worse than any cell here.
- **Points placed by the reads are right at 30 reads**: within 0.06 a tract on average in every cell, 63 samples
  included, at 256 points; the spread over tracts is at most 0.2 but in one cell (0.48 at thirteen classes and
  concentration 0.1). **At 3 reads they are not**: at low
  concentration they are worse than the fixed points — 1.9 to 4.2 too low at thirteen classes, the spread 2 to 6 —
  because three reads do not pin a sample's genotype, and a tract's reads say little about where its frequencies
  are.

## 4. What it means for the fit

Every repeat-tract stratum the fit draws its numbers from at thirteen classes is fitted to a likelihood that is
several units a tract too low, by an amount that depends on the numbers being fitted. That is what moved the
fitted spectrum and concentration in step C2. The errors C1 computes, step C3's check against a longer climb, step
D's precision target and the errors step E2 writes all stand on this likelihood.

## 5. Proposed next step (for the owner)

**Decided (owner, checkpoint C, 2026-10-02): the 256 points are precise enough.** The average that combines both
point sets is not measured, and `QUADRATURE_POINTS` and `dirichlet_points` stay as they are. The proposal below is
kept as it was made.

**Recommendation: measure an average that combines both point sets** — the fixed points and the points placed by
the reads, each tract's average taken over the two together with weights for the mixture they form (a
"defensive" importance sample). The fixed points carry the 3-read, low-concentration cells, where placing by the
reads fails, and the placed points carry the cells with many reads or many samples, where fixed points fail. The
same table will show whether 256 of each is enough, and the cost will be measured: placing points per tract builds
each tract's own point set at every evaluation, where the fixed set is built once an evaluation for every tract.

Raising the fixed count alone is not recommended: 4,096 points cost about ten times the time and still leave the
63-sample cell 6.6 a tract too low at 30 reads.
