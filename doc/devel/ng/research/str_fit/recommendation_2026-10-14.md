# The repeat-tract fit: what the research recommends

**2026-10-14, for the owner's decision.** This closes the research plan
[`../str_fit_research_plan.md`](../str_fit_research_plan.md). Every figure below comes from the
reports in this directory, named beside it. All accuracy figures are on simulated cohorts drawn
from a known model (the bench, `phase0_bench_2026-10-12.md`); real cohorts have not been scored yet.

## The recommendation

**Replace the fit in main with penalised EM, with a floor of one locus's worth on each allele
length.** For each motif length it fits, along the reference repeat count, smooth curves for the
stutter (how often a read slips, which way, how far), for the allele length variation (what
fraction of chromosomes carries each length) and for the between-locus share F (how much of that
variation lies between loci rather than among samples at one locus). It fits them directly from
the reads, in one loop: weigh each sample's possible genotypes, count with those weights, draw
smooth curves through the counts, repeat. How smooth each curve is, is estimated from the data.

**Why:**

- **It predicts loci it has not seen better than every alternative, in all 45 simulated cells
  tested** (1 to 2,000 samples, 3 to 300 reads a sample, homopolymers and dinucleotides, PCR,
  PCR-free, and dinucleotides at F = 0.4). The lead over the counting fit is large where data is
  thin and small but consistent in large cohorts (`compare_D_prior_2026-10-14.md`):

  | samples | lead over the counting fit, log units per hidden locus (homopolymers / dinucleotides) |
  |---|---|
  | 1 | +0.10 to +0.14 / +0.13 to +0.15 |
  | 10 | +0.08 to +0.10 / +0.11 to +0.15 |
  | 60 | +0.03 / +0.04 to +0.08 |
  | 500 | +0.01 to +0.03 / +0.03 to +0.04 |
  | 2,000, 3 reads | +0.009 / +0.051 |

  The fit in main predicted worse than both where it could be run (−0.27 to −0.46 a locus against
  the truth at 1 and 10 samples, where penalised EM is at −0.04 to −0.09).

- **It is thousands of times faster than the fit in main** (`cost_2026-10-14.md`). tomato2 (59
  accessions, 32 Mb): 5.8 seconds at the ±9 window, against 6.7 hours for main's fit at ±6.
  Estimated at kimura's scale (2,169 samples) on this Mac: 7–23 minutes at 3 reads a sample, up to
  28 at 30 reads, within the one-hour budget. Memory at kimura's scale: about 10 GiB.

- **Its uncertainties are honest for stutter, F and common allele lengths**: the truth falls
  within one reported error about 50–70 times in 100, as it should (`compare_D_prior_2026-10-14.md`).

## What it costs, stated plainly

1. **A slip rate 0.3–0.7% low at 3 reads a sample in cohorts of 500 or more**, largest in long,
   thin classes. It comes from the floor on the allele lengths, which buys the better prediction.
   Eight reshaped floors were tried; each removed the bias only as fast as it lost prediction
   (`shares_prior_2026-10-14.md`). Its uncertainties do not cover that bias. Whether a slip rate
   0.5% low changes any genotype call is not measured.
2. **Overconfident uncertainties on rare allele lengths** (under 1% of chromosomes) from 60 samples
   up: the truth falls within one error only 27–42 times in 100. They need a flag.
3. **An overconfident shorter share for homopolymers at 3 reads** (29–55 times in 100 within one
   error).
4. **Curves that run high at the ends of the repeat range at 1 sample and 3 reads** (+49% to +127%
   at the extreme class), because the end classes' own data there are biased and too few to show
   the curve bending (`track2_penalised_em_2026-10-12.md` §10.3).
5. **F cannot be measured at one sample**, and the fit reports it as absent.

## What its implementation spec must require

1. **Fit one motif length's classes so that memory stays bounded at deep reads.** Penalised EM
   holds every class's distinct rows of reads at once; at 2,000 samples × 30 reads it peaked at
   4.8 GB, which would not fit in 64 GiB at kimura's scale (`cost_2026-10-14.md` §9).
2. **Keep the two exact speed-ups**: each distinct row of reads scored once, and pairs of samples
   kept as four numbers. Every fitted number is unchanged by them to about 10⁻⁶.
3. **Report F as absent at one sample**, and its uncertainty with it.
4. **Flag uncertainties on rare allele lengths.**
5. **Write the between-locus share F, not the concentration, and "allele length variation" in
   place of "length spectrum"** (the owner's renaming, 2026-10-12; half done on `stutter-multiplier`).

## Decisions for the owner

1. **Adopt penalised EM with the one-locus floor as the repeat-tract fit**, and have its spec
   written. *Recommended.*
2. **Calling: keep a class's own allele length variation when its F is absent**, held at a stated
   default concentration, instead of falling to a flat shape. At one sample this makes held-out
   prediction about ten times closer to the truth (`track1_counting_fit_errors_2026-10-12.md` §13).
   *Recommended.* Without it, the one-sample gains above are mostly lost.

## Not yet measured

- **Real data.** Held-out prediction on tomato2 and kimura, and calling on HG002 and the GIAB trio,
  with the chosen method. This is the first step after the decision.
- **Calling accuracy** on simulated data with the fitted numbers.
- **The kimura shape run** (`kimura_prompt_shapes_span9.md`), handed to the owner.
- **Reading kimura's 2,169 census files**, which the one-hour budget also has to cover.
