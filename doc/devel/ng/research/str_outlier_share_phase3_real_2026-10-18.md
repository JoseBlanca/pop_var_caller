# Outlier share, phase 3: the truth-free fit on HG002's real reads

**Date:** 2026-10-18. **Plan:** [`str_outlier_share_plan.md`](str_outlier_share_plan.md) §4,
phase 3, run on real data before any calling test. **Sample:** HG002, GIAB tandem-repeat
benchmark regions, 5x / 10x / 30x, every tract's complete reads from the psp (phase 1's dumps).
**Code and output:** `experiments/str_outlier_share/real_fit.py`,
`experiments/str_outlier_share/phase3_real_fit_{5x,10x,30x}.tsv`.

## What was run

Phase 2's prototype fit, unchanged, on each stratum — one repeat-unit length and one reference
repeat count — holding at least 100 tracts (42 to 44 strata a depth, up to 3,000 tracts each). It
fits the slip rate, its direction and fall-off, the stray share and the allele-length shares,
summing over each tract's unknown genotype. It sees only reads a whole number of repeats from the
reference, as today's stutter fit does; reads a non-whole number off and reads with a flank deletion
are counted beside it. No truth is used.

## What was found

**Pooled over a period class, the share lands where phase 1's truth-based count put it, at every
depth.** Read-weighted over the class's strata:

| | 5x | 10x | 30x | phase 1, against truth (30x) |
|---|---|---|---|---|
| period 1 | 0.0031 | 0.0039 | 0.0027 | 0.0021 two or more repeats off, at clean tracts |
| period 2 | 0.0045 | 0.0044 | 0.0032 | 0.0006 |
| period 3–6 | 0.0002 | 0.0004 | 0.0002 | 0.0016 |

The phase 1 column counts reads two or more repeats off at tracts it classed as clean; the fit's
share is spread uniformly over 25 lengths, so the two are related but not the same quantity. The
fitted value sits above the count at periods 1 and 2 and below it at 3–6. The shortfall at 3–6 is
not explained here: reads three or more repeats off are 0.0013 of reads at clean tracts there, and
the fit sees them.

**Per stratum, the share is not one number.** At homopolymers it rises steadily with tract length
(30x): 0.0005 at 8 repeats, 0.0019 at 10–13, 0.004–0.008 at 16–20, 0.012–0.015 at 21–26. The
fall-off fitted alongside it is small, 0.04–0.25, so stutter's chance of slipping three or more
repeats is small, and the share is probably taking the long tail of stutter at long tracts rather
than reads from elsewhere. At dinucleotides it swings between neighbouring strata — 0.0000, 0.0202,
0.0017, 0.0107, 0.0000 from 10 to 14 repeats — so with 150–250 tracts a stratum the share and the
stutter tail are not separated.

**The reads the fit does not see are larger than the share it fits.** As shares of all reads at the
class's tracts:

| | non-whole number of repeats off | flank deletion |
|---|---|---|
| period 1 | — | 0.010–0.014 |
| period 2 | 0.0018–0.0022 | 0.0038–0.0062 |
| period 3–6 | 0.0039–0.0050 | 0.0018–0.0041 |

Flank-deletion reads include real alleles: since the held edge fix, a deletion allele that touches
the tract's edge is spelled with one. Phase 1 counted 0.0024 of reads at clean homozygous
homopolymer tracts carrying one; the 1–1.4% here is all tracts. At calling, a read whose flank
deletion differs from its allele's gets no stutter probability and falls entirely to the outlier
term, so a share fitted for calling has to account for these reads too.

## What it means

- **A fitted share near 2–4 reads in 1,000 at periods 1 and 2 is what the reads say**, at each of
  the three depths, and it agrees in size with the truth-based count.
- **Fitting it as one number per period class would be wrong at long homopolymers**, where it
  climbs to 15 in 1,000 together with a steep fitted stutter tail. Either the share is fitted per
  stratum and smoothed along repeat count, as the slip rate already is, or the stutter shape needs a
  heavier tail and the share stays flat. Which of the two is right is the next question.
- **The fit has to see the reads calling sends to the outlier term**: non-whole lengths and
  mismatched flank deletions.

## What this does not show

- Whether calling with these values beats `--defaults` at 0.07 (phase 4).
- Whether the rise with tract length is stutter's tail or real stray reads; the truth-based counts
  of phase 1 were not split by repeat count.
- One sample only; the cohort case is simulated (phase 2), not measured.
