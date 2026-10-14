# The prior on the allele-length shares: can its slip-rate bias be removed?

**Status:** research report, 2026-10-14, branch `str-fit-compare`. Follows
[`compare_D_prior_2026-10-14.md`](compare_D_prior_2026-10-14.md). That report found that penalised
counting with the counting fit's prior on the shares (**D+prior**, `D_tract_prior`) predicts held-out
loci better than the counting fit in every cell measured. It also found that the prior brings back a
small slip-rate bias at large cohorts: −0.3% to −0.7% at 500 and 2,000 samples × 3 reads, where
penalised counting without the prior (**D**) is within −0.1%. This report tries reshaped priors to
remove the bias and keep the prediction. All figures are on drawn data from the fitting bench (the
tomato-like PCR truth, at most 100 loci a class, the ±9 window); none is a cost measure.

---

## 1. The answer

**No reshaped prior removes the bias without giving up held-out prediction in proportion.** Across
the eight variants tried, the slip-rate bias and the held-out gain move together. Making the prior
smaller, or moving it to lengths far from the reference, shrinks the bias and loses the prediction in
step. Nothing measured has D+prior's prediction with less of its bias.

**Recommendation: keep D+prior as it is: one locus's worth of allele draws a class, spread evenly over
all 19 lengths** (the counting fit's own `SHARE_PRIOR_TRACTS = 1`). The trade-off it carries:

- **The cost:** at three reads and 500 to 2,000 samples, the slip rate comes out 0.34% to 0.53% low
  against D on the same draws (0.49% to 0.67% at the long classes, where the slip rate is highest). Its
  reported errors cover the truth less often there: within one error 0.38 to 0.57 times in 100,
  within two 0.76 to 0.94, against D's 0.48 to 0.73 and 0.88 to 0.97. **At 10 reads the bias is
  nearly gone** (−0.02% to −0.04% at 500 samples), so it is a low-depth effect.
- **What it buys:** held-out prediction better than D by 0.05 to 0.09 a locus at 10 samples, 0.15
  to 0.30 at 60, 0.52 to 0.96 at 500 and 0.62 to 1.31 at 2,000, and better than the counting fit in
  every cell (previous report).

For scale, the bias is about half of the slip rate's own error there. At 2,000 × 3, D+prior's slip
rate is off by 0.76% to 0.89% root-mean-square a class, of which the bias is −0.35% to −0.45%.

**One variant ties with D+prior: the same prior with the reference length's two neighbours left
out** (`D_prior_skip1`). It matches D+prior on held-out prediction within ±0.014 a locus in every
cell and removes only a twentieth to a tenth of the bias. It is not worth the extra rule.

**If the bias must go, the least costly option measured is the prior moved to lengths five or more
repeats from the reference** (`D_prior_skip4`). It halves the bias at 2,000 × 3 (−0.12% to −0.17%
against D, from −0.35% to −0.45%). It costs 0.15 to 0.48 a held-out locus at 500 to 2,000 samples,
which puts it behind the counting fit on prediction there. **Not recommended.**

---

## 2. What was tried

Every variant is a pseudo-count prior: a fixed number of extra allele draws added to each length of a
class before its shares are chosen. So every fit is still the maximum of a penalised objective, and
the sandwich and Laplace errors apply unchanged. D's own floor of one draw a class, spread over the 19
lengths, stays in every variant.

| method | the prior on each class's shares |
|---|---|
| `D` | none (only D's floor) |
| `D_tract_prior` (D+prior) | one locus's worth of draws (the class's draws ÷ its loci with reads), spread evenly over all 19 lengths |
| `D_prior_skip1` | the same total, spread evenly over the 16 lengths more than **1** repeat from the reference; none on −1, 0, +1 |
| `D_prior_skip2` | the same, over lengths more than **2** repeats away (14 lengths) |
| `D_prior_skip4` | over lengths more than **4** repeats away (10 lengths: ±5 to ±9) |
| `D_prior_skip6` | over lengths more than **6** repeats away (6 lengths: ±7 to ±9) |
| `D_prior_half` | half a locus's worth, spread over all 19 |
| `D_prior_fixed100` | **100 draws a class in all**, whatever the cohort, spread over all 19 |
| `D_prior_fixed30` | 30 draws a class in all, spread over all 19 |

- **"Lengths next to common alleles" is implemented as lengths next to the reference length.** In the
  bench's truth the reference is the most common length in every class (its share 0.50 to 0.96). The
  wider skips cover the common alleles one or two repeats from it as well. A prior that avoids each
  locus's own common alleles would need a prior per locus, which the class-level shares do not have;
  it was not built.
- **The two scalings with evidence.** A fixed total does not grow with the cohort. Half a locus's
  worth grows with it at half the rate.
- **The shape-based alternative** (a weak penalty pulling the shares towards a smooth falling tail
  rather than towards flat) was not built. It needs a new penalty term inside the shares' curve
  rather than a change to the counts. Section 3 also argues against it: the bias follows how much
  prior weight there is more than where it sits, and a falling tail concentrates weight next to the
  reference, where the skips show it costs most.

**Cells.** Tomato-like PCR truth, both motif lengths, seeds 1–3 (1–6 at one sample). Every variant
was run in the same process as D and D+prior on the same draws, so all comparisons are paired:

- **slip-rate bias:** 500 × 3, 500 × 10, 2,000 × 3;
- **held-out prediction:** 1 × 3 (both fallbacks), 10 × 3, 60 × 3, 60 × 10, 500 × 3, 500 × 10,
  2,000 × 3.

The held-out score uses 20 loci a class at 500 and 2,000 samples, every held-out locus elsewhere. The
counting fit's (B's) scores are from the earlier reports' runs on the same seeds.

## 3. Slip-rate bias

The slip rate's relative error, estimate minus D's estimate on the same draw and class, divided by the
truth, averaged over classes and seeds. Read it as "the bias the prior adds", with a standard error of
0.01% to 0.09%. "Long" means the classes of 16 to 21 repeats (one-base runs) and 14 to 27 (two-base),
where the slip rate is highest. Each cell is one-base / two-base runs. The last column is the share of
the slip rate's estimates whose truth lies within one / two reported errors at 2,000 × 3 (one-base ;
two-base); honest is 0.68 / 0.95.

| prior | 500 × 3 | 500 × 3, long | 2,000 × 3 | 2,000 × 3, long | 500 × 10 | coverage, 2,000 × 3 |
|---|---|---|---|---|---|---|
| D | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0.62 / 0.90 ; 0.73 / 0.94 |
| **D+prior** | −0.34% / −0.53% | −0.49% / −0.67% | −0.35% / −0.45% | −0.56% / −0.58% | −0.04% / −0.02% | 0.55 / 0.76 ; 0.46 / 0.84 |
| skip1 | −0.31% / −0.48% | −0.47% / −0.64% | −0.34% / −0.42% | −0.55% / −0.56% | −0.04% / −0.02% | 0.52 / 0.79 ; 0.43 / 0.89 |
| skip2 | −0.23% / −0.33% | −0.39% / −0.48% | −0.27% / −0.31% | −0.46% / −0.45% | −0.04% / −0.01% | 0.52 / 0.81 ; 0.51 / 0.92 |
| skip4 | −0.12% / −0.18% | −0.22% / −0.29% | −0.17% / −0.12% | −0.32% / −0.21% | −0.05% / −0.02% | 0.50 / 0.86 ; 0.62 / 0.92 |
| skip6 | −0.06% / −0.09% | −0.11% / −0.15% | −0.11% / −0.00% | −0.20% / −0.05% | −0.07% / −0.03% | 0.55 / 0.88 ; 0.68 / 0.94 |
| half | −0.24% / −0.36% | −0.33% / −0.46% | −0.24% / −0.32% | −0.38% / −0.41% | −0.02% / −0.02% | 0.52 / 0.81 ; 0.49 / 0.92 |
| fixed 100 | −0.14% / −0.21% | −0.19% / −0.27% | −0.06% / −0.07% | −0.09% / −0.10% | −0.01% / −0.01% | 0.60 / 0.88 ; 0.70 / 0.94 |
| fixed 30 | −0.06% / −0.09% | −0.08% / −0.11% | −0.02% / −0.03% | −0.03% / −0.04% | 0 / 0 | 0.62 / 0.90 ; 0.71 / 0.94 |

The full rows, with the absolute bias against the truth, are in the appendix.

**What it shows:**

- **The bias follows the prior's total weight.** Half the weight gives about two-thirds of the bias;
  a fixed 30 or 100 draws, a small fraction of one locus's worth at 500 and 2,000 samples, gives
  almost none.
- **Moving the weight away from the reference reduces the bias, but slowly.** Removing the prior from
  the reference and its two neighbours (skip1) cuts it by a twentieth to a tenth; from ±2 by a quarter
  to a third; from ±4 by a half to three-quarters; from ±6 by 70% to all of it. So slipped reads are absorbed by prior weight at every
  length a slip reaches, not only next to the reference. Slips of two to four repeats are common at
  the bench's fall-off of 0.1 to 0.3.
- **At 10 reads the bias is under 0.1% for every variant.** It is a three-read effect, as track 1
  found for the counting fit (§12.1, −0.19% at 30 reads).

## 4. Held-out prediction

Seed-paired difference against D+prior, a held-out locus (above zero, better than D+prior). Each cell
is one-base / two-base runs. Standard errors over seeds are 0.000 to 0.04 (appendix). At one sample
the classes keep their own shares at the stated F. With calling's fallback every D variant scores the
same (−1.628 and −1.828 a locus), because none has an F at one sample and the fallback replaces the
shares with a flat shape.

| prior | 1 × 3 | 10 × 3 | 60 × 3 | 60 × 10 | 500 × 3 | 500 × 10 | 2,000 × 3 |
|---|---|---|---|---|---|---|---|
| D | +0.002 / +0.006 | −0.055 / −0.087 | −0.151 / −0.250 | −0.184 / −0.298 | −0.572 / −0.865 | −0.524 / −0.961 | −0.620 / −1.305 |
| skip1 | −0.001 / −0.002 | +0.002 / +0.000 | +0.006 / +0.002 | +0.006 / +0.004 | +0.014 / −0.006 | +0.006 / +0.005 | +0.003 / −0.007 |
| skip2 | −0.001 / −0.004 | −0.003 / −0.007 | −0.001 / −0.022 | −0.007 / −0.030 | −0.030 / −0.097 | −0.016 / −0.105 | −0.029 / −0.109 |
| skip4 | −0.003 / −0.009 | −0.014 / −0.030 | −0.042 / −0.100 | −0.048 / −0.124 | −0.148 / −0.338 | −0.229 / −0.446 | −0.206 / −0.478 |
| skip6 | −0.006 / −0.012 | −0.030 / −0.061 | −0.096 / −0.175 | −0.104 / −0.219 | −0.381 / −0.541 | −0.341 / −0.676 | −0.435 / −0.809 |
| half | +0.001 / +0.004 | −0.015 / −0.022 | −0.031 / −0.042 | −0.036 / −0.051 | −0.073 / −0.103 | −0.077 / −0.127 | −0.056 / −0.096 |
| fixed 100 | −0.479 / −0.643 | −0.039 / −0.106 | +0.023 / +0.024 | +0.023 / +0.028 | −0.188 / −0.249 | −0.184 / −0.314 | −0.303 / −0.517 |
| fixed 30 | −0.175 / −0.281 | +0.012 / +0.005 | −0.031 / −0.040 | −0.037 / −0.053 | −0.314 / −0.438 | −0.306 / −0.534 | −0.436 / −0.730 |
| B (counting fit, after curves) | −0.102 / −0.133 | −0.100 / −0.149 | −0.034 / −0.039 | −0.027 / −0.040 | −0.010 / −0.034 | −0.026 / −0.042 | −0.009 / −0.051 |

- **Every variant that reduces the bias loses prediction at 500 and 2,000 samples**, and loses more
  the more bias it removes. Against the counting fit (the B row): half and skip2 stay ahead of B at 10
  and 60 samples and fall behind it at 500 and 2,000; skip4 and skip6 fall behind from 60 samples up.
- **A fixed total cannot fit every cohort size.** 100 draws a class swamps a one-sample class, which
  holds fewer than 200 allele draws (−0.48 / −0.64 a locus against D+prior). The same 100 draws is too
  little at 500 and 2,000 samples, where a class holds 100,000 or more. Only at 60 samples, where 100
  draws is about a locus's worth, does it beat D+prior (+0.023 / +0.024 to +0.028). The prior has to
  grow with the cohort to keep a held-out locus at an unseen length possible, and a prior that grows
  with the cohort is what brings the bias.
- **skip1 is D+prior's equal**: within ±0.014 a locus everywhere, above it in 9 of the 14 cell-period
  pairs, level in 1 and below in 4. It is ahead by more than two standard errors only on one-base
  runs at 60 × 3 and 60 × 10 (+0.006, standard error 0.001 or less) and 500 × 3 (+0.014 ± 0.004).

## 5. Small shares' error honesty

Shares whose truth is below 0.01: the share of estimates within one / two reported errors (one-base ;
two-base), pooled over classes and seeds.

| prior | 1 × 3 | 10 × 3 | 60 × 10 | 500 × 3 | 2,000 × 3 |
|---|---|---|---|---|---|
| D | 0.57 / 0.88 ; 0.56 / 0.88 | 0.44 / 0.58 ; 0.37 / 0.50 | 0.29 / 0.36 ; 0.20 / 0.25 | 0.28 / 0.32 ; 0.18 / 0.22 | 0.28 / 0.34 ; 0.19 / 0.22 |
| **D+prior** | 0.52 / 0.89 ; 0.42 / 0.87 | 0.63 / 0.85 ; 0.61 / 0.83 | 0.37 / 0.52 ; 0.35 / 0.51 | 0.37 / 0.49 ; 0.28 / 0.40 | 0.38 / 0.50 ; 0.27 / 0.39 |
| skip1 | 0.49 / 0.89 ; 0.39 / 0.85 | 0.62 / 0.87 ; 0.60 / 0.86 | 0.38 / 0.55 ; 0.38 / 0.57 | 0.35 / 0.50 ; 0.31 / 0.45 | 0.36 / 0.52 ; 0.29 / 0.45 |
| skip2 | 0.48 / 0.89 ; 0.36 / 0.83 | 0.62 / 0.86 ; 0.59 / 0.86 | 0.35 / 0.56 ; 0.39 / 0.58 | 0.35 / 0.50 ; 0.31 / 0.49 | 0.38 / 0.50 ; 0.31 / 0.48 |
| skip4 | 0.41 / 0.87 ; 0.29 / 0.77 | 0.54 / 0.81 ; 0.47 / 0.78 | 0.35 / 0.56 ; 0.37 / 0.61 | 0.31 / 0.48 ; 0.31 / 0.56 | 0.35 / 0.52 ; 0.33 / 0.56 |
| skip6 | 0.34 / 0.76 ; 0.34 / 0.67 | 0.42 / 0.68 ; 0.29 / 0.61 | 0.31 / 0.48 ; 0.17 / 0.33 | 0.28 / 0.46 ; 0.14 / 0.31 | 0.31 / 0.46 ; 0.14 / 0.30 |
| half | 0.55 / 0.90 ; 0.49 / 0.89 | 0.58 / 0.77 ; 0.54 / 0.74 | 0.33 / 0.45 ; 0.28 / 0.40 | 0.31 / 0.39 ; 0.22 / 0.31 | 0.35 / 0.42 ; 0.22 / 0.31 |
| fixed 100 | 0.00 / 0.00 ; 0.00 / 0.00 | 0.06 / 0.32 ; 0.01 / 0.09 | 0.36 / 0.58 ; 0.41 / 0.62 | 0.31 / 0.36 ; 0.22 / 0.29 | 0.29 / 0.35 ; 0.20 / 0.26 |
| fixed 30 | 0.00 / 0.00 ; 0.00 / 0.00 | 0.43 / 0.78 ; 0.35 / 0.76 | 0.34 / 0.44 ; 0.28 / 0.40 | 0.29 / 0.34 ; 0.21 / 0.27 | 0.29 / 0.34 ; 0.19 / 0.23 |

**No prior makes small shares honest from 60 samples up**: the best reach 0.31 to 0.41 within one
error. The skips do not help the small shares they leave without a prior, and skip6 makes them worse
(0.14 to 0.31 within one error). The fixed totals put every small share far above its truth at one
sample (coverage zero). D+prior and skip1 are the best at 10 samples (0.60 to 0.63 within one error)
and as good as any from 60 up.

## 6. What was not done

- **A prior that avoids each locus's common alleles** (it needs a prior per locus) and **the
  shape-based penalty** (a new term in the shares' curve). §2 gives the reasons.
- **A prior that scales with read depth.** The bias is a three-read effect; a prior that shrinks as
  depth grows would not remove it at three reads, where it matters.
- **The other two truths** (PCR-free, F 0.4) and depths above 10: the variants were run on the default
  truth only.
- **The effect on calling.** Whether a slip rate 0.5% low at long classes changes any genotype was not
  measured.

**How to rerun** one cell, every variant on the same draws:

```sh
cargo build --release --features bench-fixtures --example ng_str_fit_bench
target/release/examples/ng_str_fit_bench --period 2 --samples 500 --depth 3 --seed 1 --truth pcr \
    --loci-at-most 100 --heldout-loci 20 --both-fallbacks --out rows.tsv --methods \
    D,D_tract_prior,D_prior_skip1,D_prior_skip2,D_prior_skip4,D_prior_skip6,D_prior_half,D_prior_fixed100,D_prior_fixed30
```

The settings are `PenalisedFitConfig::{allele_prior_tracts, allele_prior_fixed_draws,
allele_prior_skips_within}`, all off by default. The rows are in this worktree's `tmp/grid_prior` and
`tmp/grid_prior2`, not committed. The draws are seeded and reproduce.

---

## Appendix A: slip-rate bias and coverage, every row

Mean relative bias against the truth (all classes, long classes), the same against D on the same
draws, and coverage of the slip rate's errors; ± is the standard error over class fits.

| cell | method | fits | mean relative bias, all classes | long classes | minus D, all classes | minus D, long classes | within 1 / 2 errors |
|---|---|---|---|---|---|---|---|
| p1_n500_d3 | D | 42 | +0.02% ± 0.24% | -0.62% ± 0.24% | +0.00% ± 0.00% | +0.00% ± 0.00% | 0.48 / 0.88 (42) |
| p1_n500_d3 | D_tract_prior | 42 | -0.32% ± 0.25% | -1.11% ± 0.24% | -0.34% ± 0.04% | -0.49% ± 0.07% | 0.38 / 0.83 (42) |
| p1_n500_d3 | D_prior_skip1 | 42 | -0.30% ± 0.26% | -1.08% ± 0.24% | -0.31% ± 0.04% | -0.47% ± 0.06% | 0.38 / 0.86 (42) |
| p1_n500_d3 | D_prior_skip2 | 42 | -0.21% ± 0.26% | -1.00% ± 0.24% | -0.23% ± 0.04% | -0.39% ± 0.05% | 0.40 / 0.86 (42) |
| p1_n500_d3 | D_prior_skip4 | 42 | -0.11% ± 0.26% | -0.84% ± 0.24% | -0.12% ± 0.04% | -0.22% ± 0.05% | 0.45 / 0.86 (42) |
| p1_n500_d3 | D_prior_skip6 | 42 | -0.04% ± 0.26% | -0.72% ± 0.24% | -0.06% ± 0.04% | -0.11% ± 0.06% | 0.48 / 0.88 (42) |
| p1_n500_d3 | D_prior_half | 42 | -0.22% ± 0.25% | -0.95% ± 0.23% | -0.24% ± 0.03% | -0.33% ± 0.05% | 0.43 / 0.86 (42) |
| p1_n500_d3 | D_prior_fixed100 | 42 | -0.12% ± 0.25% | -0.80% ± 0.23% | -0.14% ± 0.02% | -0.19% ± 0.03% | 0.45 / 0.86 (42) |
| p1_n500_d3 | D_prior_fixed30 | 42 | -0.04% ± 0.24% | -0.69% ± 0.24% | -0.06% ± 0.01% | -0.08% ± 0.01% | 0.50 / 0.86 (42) |
| p2_n500_d3 | D | 63 | -0.04% ± 0.18% | -0.08% ± 0.19% | +0.00% ± 0.00% | +0.00% ± 0.00% | 0.67 / 0.97 (63) |
| p2_n500_d3 | D_tract_prior | 63 | -0.57% ± 0.17% | -0.75% ± 0.18% | -0.53% ± 0.04% | -0.67% ± 0.04% | 0.57 / 0.94 (63) |
| p2_n500_d3 | D_prior_skip1 | 63 | -0.52% ± 0.18% | -0.72% ± 0.18% | -0.48% ± 0.04% | -0.64% ± 0.04% | 0.59 / 0.94 (63) |
| p2_n500_d3 | D_prior_skip2 | 63 | -0.38% ± 0.18% | -0.56% ± 0.19% | -0.33% ± 0.04% | -0.48% ± 0.05% | 0.62 / 0.94 (63) |
| p2_n500_d3 | D_prior_skip4 | 63 | -0.22% ± 0.18% | -0.37% ± 0.20% | -0.18% ± 0.04% | -0.29% ± 0.05% | 0.63 / 0.94 (63) |
| p2_n500_d3 | D_prior_skip6 | 63 | -0.13% ± 0.18% | -0.23% ± 0.20% | -0.09% ± 0.04% | -0.15% ± 0.06% | 0.63 / 0.94 (63) |
| p2_n500_d3 | D_prior_half | 63 | -0.41% ± 0.17% | -0.54% ± 0.18% | -0.36% ± 0.03% | -0.46% ± 0.03% | 0.63 / 0.95 (63) |
| p2_n500_d3 | D_prior_fixed100 | 63 | -0.26% ± 0.17% | -0.35% ± 0.18% | -0.21% ± 0.02% | -0.27% ± 0.02% | 0.62 / 0.98 (63) |
| p2_n500_d3 | D_prior_fixed30 | 63 | -0.13% ± 0.17% | -0.19% ± 0.18% | -0.09% ± 0.01% | -0.11% ± 0.01% | 0.62 / 0.97 (63) |
| p1_n500_d10 | D | 42 | +0.02% ± 0.19% | -0.57% ± 0.17% | +0.00% ± 0.00% | +0.00% ± 0.00% | 0.50 / 0.81 (42) |
| p1_n500_d10 | D_tract_prior | 42 | -0.02% ± 0.20% | -0.64% ± 0.19% | -0.04% ± 0.01% | -0.07% ± 0.03% | 0.50 / 0.79 (42) |
| p1_n500_d10 | D_prior_skip1 | 42 | -0.02% ± 0.20% | -0.65% ± 0.19% | -0.04% ± 0.01% | -0.08% ± 0.03% | 0.50 / 0.79 (42) |
| p1_n500_d10 | D_prior_skip2 | 42 | -0.02% ± 0.20% | -0.65% ± 0.19% | -0.04% ± 0.02% | -0.08% ± 0.03% | 0.50 / 0.79 (42) |
| p1_n500_d10 | D_prior_skip4 | 42 | -0.03% ± 0.20% | -0.66% ± 0.19% | -0.05% ± 0.02% | -0.09% ± 0.03% | 0.50 / 0.76 (42) |
| p1_n500_d10 | D_prior_skip6 | 42 | -0.05% ± 0.20% | -0.69% ± 0.20% | -0.07% ± 0.02% | -0.12% ± 0.04% | 0.50 / 0.79 (42) |
| p1_n500_d10 | D_prior_half | 42 | -0.00% ± 0.20% | -0.62% ± 0.18% | -0.02% ± 0.01% | -0.04% ± 0.02% | 0.48 / 0.79 (42) |
| p1_n500_d10 | D_prior_fixed100 | 42 | +0.01% ± 0.20% | -0.60% ± 0.18% | -0.01% ± 0.01% | -0.02% ± 0.01% | 0.48 / 0.81 (42) |
| p1_n500_d10 | D_prior_fixed30 | 42 | +0.02% ± 0.19% | -0.58% ± 0.17% | -0.00% ± 0.00% | -0.01% ± 0.00% | 0.48 / 0.81 (42) |
| p2_n500_d10 | D | 63 | -0.08% ± 0.09% | -0.05% ± 0.09% | +0.00% ± 0.00% | +0.00% ± 0.00% | 0.67 / 0.97 (63) |
| p2_n500_d10 | D_tract_prior | 63 | -0.10% ± 0.09% | -0.07% ± 0.10% | -0.02% ± 0.01% | -0.03% ± 0.02% | 0.65 / 0.97 (63) |
| p2_n500_d10 | D_prior_skip1 | 63 | -0.10% ± 0.09% | -0.07% ± 0.10% | -0.02% ± 0.01% | -0.02% ± 0.02% | 0.67 / 0.97 (63) |
| p2_n500_d10 | D_prior_skip2 | 63 | -0.09% ± 0.09% | -0.06% ± 0.10% | -0.01% ± 0.01% | -0.02% ± 0.02% | 0.67 / 0.97 (63) |
| p2_n500_d10 | D_prior_skip4 | 63 | -0.10% ± 0.09% | -0.07% ± 0.10% | -0.02% ± 0.01% | -0.02% ± 0.02% | 0.67 / 0.97 (63) |
| p2_n500_d10 | D_prior_skip6 | 63 | -0.11% ± 0.09% | -0.08% ± 0.10% | -0.03% ± 0.02% | -0.04% ± 0.03% | 0.67 / 0.97 (63) |
| p2_n500_d10 | D_prior_half | 63 | -0.10% ± 0.09% | -0.06% ± 0.09% | -0.02% ± 0.01% | -0.02% ± 0.01% | 0.63 / 0.97 (63) |
| p2_n500_d10 | D_prior_fixed100 | 63 | -0.09% ± 0.09% | -0.05% ± 0.09% | -0.01% ± 0.00% | -0.01% ± 0.01% | 0.67 / 0.97 (63) |
| p2_n500_d10 | D_prior_fixed30 | 63 | -0.08% ± 0.09% | -0.05% ± 0.09% | -0.00% ± 0.00% | -0.00% ± 0.00% | 0.67 / 0.97 (63) |
| p1_n2000_d3 | D | 42 | +0.00% ± 0.14% | -0.31% ± 0.12% | +0.00% ± 0.00% | +0.00% ± 0.00% | 0.62 / 0.90 (42) |
| p1_n2000_d3 | D_tract_prior | 42 | -0.35% ± 0.15% | -0.87% ± 0.17% | -0.35% ± 0.05% | -0.56% ± 0.09% | 0.55 / 0.76 (42) |
| p1_n2000_d3 | D_prior_skip1 | 42 | -0.33% ± 0.15% | -0.86% ± 0.17% | -0.34% ± 0.05% | -0.55% ± 0.08% | 0.52 / 0.79 (42) |
| p1_n2000_d3 | D_prior_skip2 | 42 | -0.26% ± 0.15% | -0.77% ± 0.16% | -0.27% ± 0.04% | -0.46% ± 0.07% | 0.52 / 0.81 (42) |
| p1_n2000_d3 | D_prior_skip4 | 42 | -0.17% ± 0.15% | -0.63% ± 0.16% | -0.17% ± 0.04% | -0.32% ± 0.06% | 0.50 / 0.86 (42) |
| p1_n2000_d3 | D_prior_skip6 | 42 | -0.11% ± 0.15% | -0.52% ± 0.16% | -0.11% ± 0.04% | -0.20% ± 0.06% | 0.55 / 0.88 (42) |
| p1_n2000_d3 | D_prior_half | 42 | -0.24% ± 0.15% | -0.69% ± 0.16% | -0.24% ± 0.03% | -0.38% ± 0.06% | 0.52 / 0.81 (42) |
| p1_n2000_d3 | D_prior_fixed100 | 42 | -0.05% ± 0.14% | -0.40% ± 0.13% | -0.06% ± 0.01% | -0.09% ± 0.02% | 0.60 / 0.88 (42) |
| p1_n2000_d3 | D_prior_fixed30 | 42 | -0.02% ± 0.14% | -0.35% ± 0.12% | -0.02% ± 0.00% | -0.03% ± 0.01% | 0.62 / 0.90 (42) |
| p2_n2000_d3 | D | 63 | -0.00% ± 0.10% | -0.11% ± 0.09% | +0.00% ± 0.00% | +0.00% ± 0.00% | 0.73 / 0.94 (63) |
| p2_n2000_d3 | D_tract_prior | 63 | -0.45% ± 0.11% | -0.69% ± 0.11% | -0.45% ± 0.04% | -0.58% ± 0.06% | 0.46 / 0.84 (63) |
| p2_n2000_d3 | D_prior_skip1 | 63 | -0.42% ± 0.11% | -0.67% ± 0.11% | -0.42% ± 0.04% | -0.56% ± 0.06% | 0.43 / 0.89 (63) |
| p2_n2000_d3 | D_prior_skip2 | 63 | -0.31% ± 0.11% | -0.56% ± 0.11% | -0.31% ± 0.04% | -0.45% ± 0.06% | 0.51 / 0.92 (63) |
| p2_n2000_d3 | D_prior_skip4 | 63 | -0.12% ± 0.11% | -0.32% ± 0.12% | -0.12% ± 0.04% | -0.21% ± 0.06% | 0.62 / 0.92 (63) |
| p2_n2000_d3 | D_prior_skip6 | 63 | -0.01% ± 0.11% | -0.16% ± 0.12% | -0.00% ± 0.05% | -0.05% ± 0.07% | 0.68 / 0.94 (63) |
| p2_n2000_d3 | D_prior_half | 63 | -0.32% ± 0.10% | -0.53% ± 0.10% | -0.32% ± 0.03% | -0.41% ± 0.04% | 0.49 / 0.92 (63) |
| p2_n2000_d3 | D_prior_fixed100 | 63 | -0.08% ± 0.10% | -0.21% ± 0.09% | -0.07% ± 0.01% | -0.10% ± 0.01% | 0.70 / 0.94 (63) |
| p2_n2000_d3 | D_prior_fixed30 | 63 | -0.03% ± 0.10% | -0.15% ± 0.09% | -0.03% ± 0.00% | -0.04% ± 0.01% | 0.71 / 0.94 (63) |

## Appendix B: held-out prediction with standard errors

Log-likelihood a held-out locus against the truth, mean over seeds ± standard error, and each
variant paired against D+prior (`D_tract_prior`). One-sample cells keep their own shares.

| cell | seeds | D | D_tract_prior | D_prior_skip1 | D_prior_skip2 | D_prior_skip4 | D_prior_skip6 | D_prior_half | D_prior_fixed100 | D_prior_fixed30 | D − D_tract_prior | D_prior_skip1 − D_tract_prior | D_prior_skip2 − D_tract_prior | D_prior_skip4 − D_tract_prior | D_prior_skip6 − D_tract_prior | D_prior_half − D_tract_prior | D_prior_fixed100 − D_tract_prior | D_prior_fixed30 − D_tract_prior |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| period 1, 1 samples, 3 reads, pcr | 6 | -0.037 ± 0.004 | -0.039 ± 0.004 | -0.040 ± 0.004 | -0.040 ± 0.004 | -0.042 ± 0.004 | -0.044 ± 0.004 | -0.038 ± 0.004 | -0.518 ± 0.006 | -0.214 ± 0.003 | +0.002 ± 0.001 | -0.001 ± 0.000 | -0.001 ± 0.000 | -0.003 ± 0.000 | -0.006 ± 0.001 | +0.001 ± 0.000 | -0.479 ± 0.009 | -0.175 ± 0.005 |
| period 1, 10 samples, 3 reads, pcr | 3 | -0.141 ± 0.007 | -0.087 ± 0.009 | -0.085 ± 0.009 | -0.090 ± 0.011 | -0.100 ± 0.009 | -0.117 ± 0.011 | -0.102 ± 0.009 | -0.126 ± 0.008 | -0.075 ± 0.009 | -0.055 ± 0.004 | +0.002 ± 0.000 | -0.003 ± 0.004 | -0.014 ± 0.003 | -0.030 ± 0.002 | -0.015 ± 0.001 | -0.039 ± 0.006 | +0.012 ± 0.001 |
| period 1, 60 samples, 3 reads, pcr | 3 | -0.344 ± 0.010 | -0.193 ± 0.002 | -0.187 ± 0.002 | -0.194 ± 0.001 | -0.235 ± 0.007 | -0.289 ± 0.007 | -0.224 ± 0.001 | -0.170 ± 0.003 | -0.224 ± 0.002 | -0.151 ± 0.012 | +0.006 ± 0.001 | -0.001 ± 0.002 | -0.042 ± 0.009 | -0.096 ± 0.008 | -0.031 ± 0.002 | +0.023 ± 0.002 | -0.031 ± 0.001 |
| period 1, 60 samples, 10 reads, pcr | 3 | -0.421 ± 0.011 | -0.237 ± 0.019 | -0.231 ± 0.019 | -0.244 ± 0.018 | -0.285 ± 0.021 | -0.340 ± 0.015 | -0.273 ± 0.018 | -0.214 ± 0.018 | -0.274 ± 0.018 | -0.184 ± 0.014 | +0.006 ± 0.000 | -0.007 ± 0.002 | -0.048 ± 0.010 | -0.104 ± 0.004 | -0.036 ± 0.001 | +0.023 ± 0.000 | -0.037 ± 0.001 |
| period 1, 500 samples, 3 reads, pcr | 3 | -0.989 ± 0.053 | -0.417 ± 0.017 | -0.403 ± 0.013 | -0.446 ± 0.008 | -0.564 ± 0.015 | -0.798 ± 0.031 | -0.490 ± 0.022 | -0.605 ± 0.032 | -0.730 ± 0.038 | -0.572 ± 0.040 | +0.014 ± 0.004 | -0.030 ± 0.013 | -0.148 ± 0.031 | -0.381 ± 0.022 | -0.073 ± 0.007 | -0.188 ± 0.017 | -0.314 ± 0.025 |
| period 1, 500 samples, 10 reads, pcr | 3 | -0.959 ± 0.045 | -0.435 ± 0.039 | -0.429 ± 0.034 | -0.451 ± 0.032 | -0.664 ± 0.032 | -0.776 ± 0.060 | -0.512 ± 0.040 | -0.619 ± 0.028 | -0.740 ± 0.031 | -0.524 ± 0.022 | +0.006 ± 0.005 | -0.016 ± 0.007 | -0.229 ± 0.024 | -0.341 ± 0.024 | -0.077 ± 0.005 | -0.184 ± 0.012 | -0.306 ± 0.012 |
| period 1, 2000 samples, 3 reads, pcr | 3 | -1.138 ± 0.025 | -0.518 ± 0.014 | -0.515 ± 0.020 | -0.547 ± 0.028 | -0.724 ± 0.034 | -0.953 ± 0.040 | -0.573 ± 0.009 | -0.821 ± 0.006 | -0.954 ± 0.041 | -0.620 ± 0.037 | +0.003 ± 0.006 | -0.029 ± 0.027 | -0.206 ± 0.030 | -0.435 ± 0.043 | -0.056 ± 0.006 | -0.303 ± 0.008 | -0.436 ± 0.039 |
| period 2, 1 samples, 3 reads, pcr | 6 | -0.047 ± 0.003 | -0.053 ± 0.002 | -0.055 ± 0.002 | -0.057 ± 0.002 | -0.062 ± 0.002 | -0.065 ± 0.003 | -0.049 ± 0.003 | -0.696 ± 0.006 | -0.334 ± 0.004 | +0.006 ± 0.001 | -0.002 ± 0.000 | -0.004 ± 0.001 | -0.009 ± 0.000 | -0.012 ± 0.000 | +0.004 ± 0.001 | -0.643 ± 0.007 | -0.281 ± 0.004 |
| period 2, 10 samples, 3 reads, pcr | 3 | -0.193 ± 0.016 | -0.106 ± 0.009 | -0.106 ± 0.009 | -0.113 ± 0.008 | -0.137 ± 0.007 | -0.167 ± 0.015 | -0.129 ± 0.012 | -0.212 ± 0.006 | -0.101 ± 0.006 | -0.087 ± 0.007 | +0.000 ± 0.001 | -0.007 ± 0.002 | -0.030 ± 0.003 | -0.061 ± 0.006 | -0.022 ± 0.003 | -0.106 ± 0.010 | +0.005 ± 0.005 |
| period 2, 60 samples, 3 reads, pcr | 3 | -0.502 ± 0.045 | -0.252 ± 0.016 | -0.250 ± 0.015 | -0.274 ± 0.015 | -0.352 ± 0.018 | -0.427 ± 0.042 | -0.294 ± 0.023 | -0.228 ± 0.010 | -0.292 ± 0.023 | -0.250 ± 0.029 | +0.002 ± 0.001 | -0.022 ± 0.006 | -0.100 ± 0.014 | -0.175 ± 0.026 | -0.042 ± 0.007 | +0.024 ± 0.006 | -0.040 ± 0.007 |
| period 2, 60 samples, 10 reads, pcr | 3 | -0.590 ± 0.027 | -0.292 ± 0.023 | -0.288 ± 0.024 | -0.322 ± 0.028 | -0.417 ± 0.031 | -0.511 ± 0.030 | -0.343 ± 0.025 | -0.264 ± 0.022 | -0.345 ± 0.025 | -0.298 ± 0.009 | +0.004 ± 0.001 | -0.030 ± 0.006 | -0.124 ± 0.010 | -0.219 ± 0.011 | -0.051 ± 0.002 | +0.028 ± 0.001 | -0.053 ± 0.002 |
| period 2, 500 samples, 3 reads, pcr | 3 | -1.374 ± 0.044 | -0.509 ± 0.011 | -0.515 ± 0.018 | -0.606 ± 0.029 | -0.847 ± 0.040 | -1.050 ± 0.028 | -0.612 ± 0.014 | -0.758 ± 0.018 | -0.947 ± 0.033 | -0.865 ± 0.046 | -0.006 ± 0.009 | -0.097 ± 0.021 | -0.338 ± 0.037 | -0.541 ± 0.025 | -0.103 ± 0.003 | -0.249 ± 0.009 | -0.438 ± 0.027 |
| period 2, 500 samples, 10 reads, pcr | 3 | -1.568 ± 0.054 | -0.608 ± 0.011 | -0.603 ± 0.006 | -0.712 ± 0.031 | -1.054 ± 0.024 | -1.284 ± 0.023 | -0.735 ± 0.014 | -0.922 ± 0.022 | -1.142 ± 0.033 | -0.961 ± 0.044 | +0.005 ± 0.009 | -0.105 ± 0.029 | -0.446 ± 0.031 | -0.676 ± 0.022 | -0.127 ± 0.005 | -0.314 ± 0.012 | -0.534 ± 0.023 |
| period 2, 2000 samples, 3 reads, pcr | 3 | -1.951 ± 0.073 | -0.645 ± 0.025 | -0.653 ± 0.033 | -0.755 ± 0.029 | -1.124 ± 0.025 | -1.454 ± 0.050 | -0.741 ± 0.032 | -1.163 ± 0.015 | -1.376 ± 0.004 | -1.305 ± 0.096 | -0.007 ± 0.009 | -0.109 ± 0.023 | -0.478 ± 0.012 | -0.809 ± 0.026 | -0.096 ± 0.007 | -0.517 ± 0.012 | -0.730 ± 0.021 |
