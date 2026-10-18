# Outlier share, phase 1: where unexplained reads come from at repeat tracts

**Date:** 2026-10-18. **Plan:** [`str_outlier_share_plan.md`](str_outlier_share_plan.md) §4, phase 1.
**Sample:** HG002, GIAB tandem-repeat benchmark regions, at 5x, 10x and 30x, mapped with minimap2.
**Evidence:** every repeat tract's complete reads as the psp stores them, so tracts that calling
skips for having no reads off the reference are counted too. psps generated with the held
edge-deletion fix. **Truth:** the assembly-based benchmark VCF. Scripts: `tmp/p1/phase1.py`,
`tmp/p1/phase1b.py` in the `str-group-b` worktree.

## What was found

**Reads that one step of stutter cannot explain come from two different populations, and only one
of them is a share of reads.** At truth-homozygous tracts, counting reads two or more repeats from
the true length, a non-whole number of repeats off, or carrying a flank deletion:

- **Clean tracts — 95 to 98 in 100.** About 1 to 7 reads in 1,000 land beyond one repeat. This
  background is the same at every depth.
- **Noisy tracts — 2 to 5 in 100.** About a third of their reads land elsewhere (the middle
  class), or nearly all of them do (reads and truth disagree outright, about 1 tract in 100). The
  shares hold steady from 5x to 30x, so this is a property of the genome or the truth, not of the
  sampling.

Three tract classes, fitted per period class to the count of such reads at each tract:

| depth, period | clean share, rate | middle share, rate | outright-disagreement share, rate |
|---|---|---|---|
| 30x, 1 | 0.948, 0.0052 | 0.038, 0.31 | 0.015, 0.94 |
| 30x, 2 | 0.965, 0.0015 | 0.024, 0.24 | 0.011, 0.95 |
| 30x, 3–6 | 0.978, 0.0031 | 0.016, 0.42 | 0.007, 0.99 |
| 10x, 1 | 0.955, 0.0057 | 0.031, 0.33 | 0.015, 0.94 |
| 10x, 2 | 0.972, 0.0018 | 0.017, 0.32 | 0.011, 0.98 |
| 5x, 1 | 0.968, 0.0066 | 0.021, 0.48 | 0.012, 0.99 |

One rate is far from the data: at 30x homopolymers it predicts 0.4 tracts with six or more such
reads where 256 are observed. A beta-binomial (reads correlated within a tract) fits better than
one rate but over-predicts the long tail by half. The three classes fit best, except at one read,
where they still over-predict by about a fifth.

## What the clean background is made of

Share of all reads at clean tracts (posterior above 0.9 of being clean), 30x:

| | one repeat off | two repeats off | three or more | non-whole | flank deletion |
|---|---|---|---|---|---|
| period 1 | 0.0291 | 0.0012 | 0.0009 | — | 0.0024 |
| period 2 | 0.0125 | 0.0004 | 0.0002 | 0.0001 | 0.0005 |
| period 3–6 | 0.0071 | 0.0003 | 0.0013 | 0.0011 | 0.0003 |

- **One and two repeats off are stutter**, which the fitted stutter shape already prices.
- **Three or more repeats off and non-whole lengths are the stray reads proper:** about 1 read in
  1,000 at periods 1 and 2, and 2 in 1,000 at longer motifs. That agrees in order of size with the
  2026-09 count at homozygous-reference tracts (1 in 2,300 and 1 in 209).
- **Flank-deletion reads** are counted separately because the calling likelihood treats them
  separately: a read whose flank deletion differs from the allele's gets no stutter probability
  and falls wholly to the outlier term. At homopolymers they are the largest of the three, at
  2.4 in 1,000.

## What it means for the outlier share

- **The literal share a fit would find is about 0.001 to 0.005**, depending on period and on
  whether flank-deletion reads count. The stated value is 0.20, and the fitted model's best tuned
  value is 0.03.
- **Calling with a near-literal share costs nothing at 5x and 10x and 27 tracts in 8,044 at 30x.**
  Fitted stutter and prior, with the edge fix, correct tract genotypes over the tracts every arm
  built:

  | outlier share | 5x | 10x | 30x |
  |---|---|---|---|
  | 0.20 (shipped) | 2,811 | 4,622 | 7,216 |
  | 0.03 (best tuned) | 2,874 | 4,715 | 7,228 |
  | 0.01 | 2,876 | 4,710 | 7,213 |
  | 0.003 (about the literal value) | 2,876 | 4,710 | 7,201 |
  | `--defaults` with the edge fix | 2,870 | 4,706 | 7,243 |

- **The 30x difference is the noisy tracts.** A read-level floor bounds each read separately, so
  it bounds a tract's total evidence less tightly as depth grows. The noisy tracts carry 30% or
  more of their reads elsewhere, which at 30x is several reads a tract; a read-level share near 0.2
  turns some of them back into homozygotes and a share near 0.003 does not. That is the depth-
  dependent over-confidence the 2026-10-18 calibration study found at 30x.
- **The noisy class is the repeat-tract counterpart of the generic path's noisy sites**, at a
  similar size: 2–5 tracts in 100 here, 0.4–1.4 sites in 100 there. The 2026-08-19 decision not to
  carry such a class in the calling likelihood was made for the generic path; whether it carries
  over is the next decision (see the plan, §5).

## What this phase did not do

- No truth-free estimate yet: every number above uses the truth to say which tracts are
  homozygous. A fit cannot.
- Tomato (3x, a cohort) is not measured yet; there the noisy class's cohort signature — the same
  tract noisy in every sample — is testable.
- The `--defaults` arm's preference for 0.20 is not decomposed here. The 2026-09 study found 0.20
  still best when its flat-prior harness used measured slippage numbers instead of the shipped
  ones (`tract_genotype_investigation/readmodel.md` §3), which points at the prior (candidate E in
  the plan) rather than the stutter shape.
