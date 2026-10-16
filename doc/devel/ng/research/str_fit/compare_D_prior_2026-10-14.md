# Penalised counting with the counting fit's prior on the shares, across the grid

**Status:** research report, 2026-10-14, branch `str-fit-compare`. Follows
[`compare_B_D_2026-10-13.md`](compare_B_D_2026-10-13.md), which found that penalised counting (D)
loses held-out prediction to the counting fit (B) from 60 samples up because of its floor under the
allele-length shares, and that D with B's prior on the shares (`D_tract_prior`) won the few cells it
was run in. This report runs that variant across the same grid. All figures are on drawn data from
the fitting bench; **no time in it is a cost measure** (the runs shared the machine with another
agent's profiling).

**The variant.** `D_tract_prior` is D unchanged — every class of a motif length fitted at once, each
number a penalised curve, F read off pairs of samples, its own floor of one allele draw a class —
with B's prior added: each round, every allele length gets the class's allele draws divided by its
loci with reads and by the 19 lengths, one locus's worth spread evenly (`allele_prior_tracts = 1`,
B's `SHARE_PRIOR_TRACTS`). Below it is called **D+prior**.

---

## 1. The answer

**D+prior predicts held-out loci better than B in every cell where both were scored: 45 cells, both
motif lengths, three truths, 1 to 2,000 samples, 3 to 300 reads.** The seed-paired difference, in
log-likelihood units a held-out locus (above zero, D+prior is closer to the truth), with its standard
error over seeds:

| samples | one-base runs | two-base runs |
|---|---|---|
| 1 (3 to 300 reads, own shares kept) | +0.10 to +0.14 (± 0.01) | +0.13 to +0.15 (± 0.01) |
| 10 (3 to 300 reads) | +0.08 to +0.10 (± 0.01) | +0.11 to +0.15 (± 0.01 to 0.02) |
| 60 (3 to 300 reads) | +0.026 to +0.034 (± 0.004) | +0.039 to +0.081 (± 0.005 to 0.03) |
| 500 (3 and 10 reads) | +0.010 ± 0.006, +0.026 ± 0.006 | +0.034 ± 0.008, +0.042 ± 0.010 |
| 2,000 (3 reads) | +0.009 ± 0.003 | +0.051 ± 0.008 |

**The margin shrinks as the cohort grows**: from more than half of B's whole loss at 1 and 10
samples, and a tenth to a fifth at 60, to 2% of it (one-base runs) and 7% (two-base) at 2,000
samples. On one-base runs at 500
samples and 3 reads it is 1.7 standard errors from zero, the least certain cell. The PCR-free truth
and two-base runs at F = 0.4 give the same sign in all 15 of their cells (§2.2).

**What the prior costs: the slip rate's small bias, the same as B's.** From 500 samples up D+prior's
slip rate comes out 0.3% to 0.7% low (median over classes), as B's does, where D without the prior
is within 0.1%. Its slip rate's error grows to B's (two-base runs at 2,000 × 3: relative error 0.0089
against D's 0.0042 and B's 0.0087), and its errors cover it less often (within one error 0.38 to 0.57
times in 100 at 500 and 2,000 samples and three reads, against D's 0.48 to 0.73). This is track 1
§12.1's trade between stutter bias and held-out prediction, now paid by D as well. At 1 to 60
samples D+prior recovers every number as well as D does.

**Error honesty** (§3):

- **The stutter numbers and F are honest** in most cells, with two exceptions below.
- **B's prior does not make small shares honest from 60 samples up; it brings D to B's level.**
  Shares below 0.01 lie within one error 0.27 to 0.42 times in 100 (D alone 0.18 to 0.31, B 0.25 to
  0.38). At 10 samples it does much better than either (0.46 to 0.63 against B's 0.25 to 0.38). At one
  sample D+prior's small shares are nearly honest within two errors (0.87 to 0.94) and somewhat low
  within one (0.42 to 0.55).
- **D's shorter share at three reads is still too narrow on one-base runs**: within one error 0.29 to
  0.55 times in 100 at 1, 10, 60, 500 and 2,000 samples. On two-base runs it is 0.43 at 10 samples
  and 0.65 to 0.76 from 60 up.

**So D+prior is the better predictor everywhere it was run, its errors are as honest as B's, and it
carries B's small slip-rate bias at large cohorts.** Whether a prior shaped to stay off the lengths
next to common alleles (track 1 §12.1's option 3) keeps the prediction and drops the bias is the next
question; it was not tried here.

---

## 2. Held-out prediction

### 2.1 The default truth

Seed-paired difference, D+prior minus B after the caller's curves, a held-out locus. Seeds: 6 at one
sample, 3 elsewhere. At one sample both keep their own shares; with calling's fallback both lose most
of what they fitted, as [`compare_B_D_2026-10-13.md`](compare_B_D_2026-10-13.md) §5 found for D.

| samples | reads | one-base runs | two-base runs |
|---|---|---|---|
| 1 | 3 | +0.102 ± 0.010 | +0.133 ± 0.008 |
| 1 | 10 | +0.144 ± 0.014 | +0.151 ± 0.013 |
| 1 | 30 | +0.106 ± 0.011 | +0.143 ± 0.009 |
| 1 | 300 | +0.115 ± 0.007 | +0.141 ± 0.011 |
| 10 | 3 | +0.100 ± 0.014 | +0.149 ± 0.012 |
| 10 | 10 | +0.088 ± 0.008 | +0.147 ± 0.016 |
| 10 | 30 | +0.095 ± 0.009 | +0.105 ± 0.012 |
| 10 | 300 | +0.082 ± 0.006 | +0.113 ± 0.018 |
| 60 | 3 | +0.034 ± 0.002 | +0.039 ± 0.005 |
| 60 | 10 | +0.027 ± 0.004 | +0.040 ± 0.006 |
| 60 | 30 | +0.026 ± 0.004 | +0.081 ± 0.031 |
| 60 | 300 | +0.031 ± 0.003 | +0.054 ± 0.010 |
| 500 | 3 | +0.010 ± 0.006 | +0.034 ± 0.008 |
| 500 | 10 | +0.026 ± 0.006 | +0.042 ± 0.010 |
| 2,000 | 3 | +0.009 ± 0.003 | +0.051 ± 0.008 |

For scale, B's own loss against the truth: 0.14 to 0.20 a locus at one sample (own shares), 0.19 to
0.30 at 10, 0.23 to 0.38 at 60, 0.43 to 0.65 at 500, 0.53 to 0.70 at 2,000 (appendix).

### 2.2 The other two truths

| truth | cell | seeds | one-base runs | two-base runs |
|---|---|---|---|---|
| PCR-free | 1 × 3 (own shares) | 3 | +0.107 ± 0.005 | +0.162 ± 0.004 |
| PCR-free | 1 × 300 (own shares) | 3 | +0.110 ± 0.002 | +0.131 ± 0.019 |
| PCR-free | 60 × 10 | 3 | +0.029 ± 0.008 | +0.040 ± 0.002 |
| PCR-free | 500 × 3 | 3 | +0.012 ± 0.003 | +0.029 ± 0.006 |
| PCR-free | 2,000 × 3 | 1 | +0.009 | +0.024 |
| F = 0.4 | 1 × 3 (own shares) | 3 | | +0.118 ± 0.009 |
| F = 0.4 | 1 × 300 (own shares) | 3 | | +0.128 ± 0.002 |
| F = 0.4 | 60 × 10 | 3 | | +0.055 ± 0.003 |
| F = 0.4 | 500 × 3 | 3 | | +0.044 ± 0.010 |
| F = 0.4 | 2,000 × 3 | 1 | | +0.040 |

"F = 0.4" is the default truth with its between-locus share replaced by 0.4 in every class (tomato2's
two-base classes give 0.29 to 0.50), two-base runs only.

## 3. Error honesty

Share of estimates whose truth lies within one / two reported errors, pooled over classes and seeds
(honest: about 0.68 / 0.95). Default truth; ranges over the cells named.

| number | cells | B | D | D+prior |
|---|---|---|---|---|
| shares < 0.01 | 1 sample | 0.18–0.26 / 0.26–0.35 | 0.56–0.63 / 0.87–0.93 | 0.42–0.55 / 0.87–0.94 |
| | 10 samples | 0.25–0.38 / 0.35–0.52 | 0.31–0.45 / 0.43–0.60 | **0.46–0.63 / 0.66–0.85** |
| | 60 samples | 0.26–0.36 / 0.36–0.51 | 0.18–0.31 / 0.23–0.36 | 0.31–0.42 / 0.46–0.58 |
| | 500 and 2,000 samples | 0.25–0.38 / 0.36–0.50 | 0.18–0.28 / 0.21–0.35 | 0.27–0.38 / 0.38–0.50 |
| shares 0.01–0.05 | 10 samples | 0.41–0.54 / 0.52–0.69 | 0.42–0.55 / 0.56–0.72 | 0.54–0.63 / 0.68–0.79 |
| | 60 to 2,000 samples | 0.38–0.55 / 0.50–0.70 | 0.37–0.55 / 0.49–0.70 | 0.42–0.55 / 0.51–0.70 |
| shares ≥ 0.05 | all | 0.54–0.71 / 0.78–0.93 | 0.52–0.69 / 0.79–0.94 | 0.50–0.69 / 0.82–0.95 |
| shorter share, one-base, 3 reads | 1, 10, 60, 500, 2,000 samples | 0.57–0.76 / 0.83–0.98 | 0.38–0.61 / 0.62–0.90 | 0.29–0.55 / 0.62–0.90 |
| shorter share, two-base, 3 reads | same | 0.62–0.75 / 0.90–0.98 | 0.49–0.68 / 0.76–0.98 | 0.43–0.76 / 0.75–0.98 |
| slip rate, 3 reads | 500 and 2,000 samples | 0.38–0.56 / 0.74–0.89 | 0.48–0.73 / 0.86–0.97 | 0.38–0.57 / 0.76–0.94 |
| F | 10 to 2,000 samples | 0.50–0.89 / 0.83–1.00 | 0.43–0.83 / 0.78–1.00 | 0.41–0.81 / 0.73–1.00 |

- **Small shares.** The prior lifts D's small shares back to B's coverage from 60 samples up, not to
  honest coverage. The reason is the one tracks 1 and 2 gave: at F = 0.7 a share of 0.005 at 100 loci
  is carried by a fraction of one locus, and an error read off the loci's spread cannot see a length
  no training locus carried. At 10 samples, where a locus's worth of draws is about 20, the prior
  holds such shares near their truth and the errors cover them in half to two-thirds of cases.
- **The shorter share at three reads** stays too narrow on one-base runs at every cohort size; the
  prior does not touch it. Its errors there are too small by roughly a third to a half.
- **The slip rate at 500 and 2,000 samples and three reads** loses the coverage D had and falls to
  B's, from the prior's bias (§1).
- F's coverage at 10 samples and 30 reads is the lowest of any method (0.41 / 0.73 two-base); the
  other F cells are between 0.56 and 0.81 within one error.

## 4. Recovery

Median over classes of each class's root-mean-square error over the seeds; the slip rate relative.

| cell | number | B after curves | D | D+prior |
|---|---|---|---|---|
| one-base, 1 × 3 | slip rate | 0.176 | 0.087 | 0.087 |
| two-base, 1 × 3 | slip rate | 0.153 | 0.133 | 0.132 |
| | reference-length share | 0.087 | 0.049 | 0.048 |
| two-base, 10 × 3 | slip rate | 0.078 | 0.045 | 0.047 |
| | F | 0.059 | 0.022 | 0.021 |
| one-base, 60 × 3 | slip rate | 0.019 | 0.022 | 0.026 |
| two-base, 60 × 10 | slip rate | 0.020 | 0.013 | 0.013 |
| two-base, 500 × 3 | slip rate | 0.013 | 0.010 | 0.012 |
| | reference-length share | 0.054 | 0.051 | 0.054 |
| two-base, 500 × 10 | reference-length share | 0.070 | 0.076 | 0.069 |
| one-base, 2,000 × 3 | slip rate | 0.0075 | 0.0063 | 0.0076 |
| two-base, 2,000 × 3 | slip rate | 0.0087 | 0.0042 | 0.0089 |
| | reference-length share | 0.067 | 0.077 | 0.067 |
| one-base, PCR-free, 1 × 3 | slip rate | 0.670 | 0.319 | 0.319 |
| two-base, F 0.4, 500 × 3 | slip rate | 0.016 | | 0.011 |

**From 1 to 60 samples D+prior recovers every number as D does**, the stutter well ahead of B at one
sample and on two-base runs. **From 500 samples up the prior moves the slip rate's error to B's and
the reference-length share's error back to B's**: it trades D's better stutter for B's shares, and
the held-out score (§2) says the shares matter more for prediction.

## 5. What was run

| cells | truth | seeds | B | D+prior |
|---|---|---|---|---|
| 1 × 3, 10, 30, 300 | default | 6 | from the earlier report | new |
| 10 and 60 × 3, 10, 30, 300 | default | 3 | from the earlier report | new |
| 500 × 3 and 10 | default | 3 | from the earlier report | new |
| 2,000 × 3 | default | 3 | one-base seeds 2–3 and two-base seed 3 new; the rest from the earlier report | new |
| 1 × 3, 1 × 300, 60 × 10 | PCR-free; F 0.4 | 3 | from the earlier report | new |
| 500 × 3 | PCR-free; F 0.4 | 3 | new | new |
| 2,000 × 3 | PCR-free; F 0.4 | 1 | new | new |

Held-out scores at 500 and 2,000 samples use 20 held-out loci a class; elsewhere every held-out
locus. **Not run:** D+prior's held-out score at 500 × 30 and 300 and 2,000 × 10, 30 and 300 (B has
none there either); seeds 2 and 3 of the other truths at 2,000 samples. Runs were held to four
processes at a time and 2,000-sample runs to two; the lowest free memory seen was 30.6 GiB.

**How to rerun** one cell:

```sh
target/release/examples/ng_str_fit_bench --period 2 --samples 500 --depth 3 --seed 1 --truth pcr \
    --loci-at-most 100 --methods B,D_tract_prior --both-fallbacks --heldout-loci 20 --out rows.tsv
uv run --no-project --with numpy python scripts/research/str_fit_summarise.py tmp/grid tmp/grid_b \
    tmp/grid_dtp --compare B_raw,B_curves,D,D_tract_prior
```

The rows are in this worktree's `tmp/grid_dtp` (with the earlier report's in `tmp/grid` and
`tmp/grid_b`), not committed; the draws are seeded and reproduce.

---

## Appendix A: one sample, own shares kept

Held-out score against the truth, a locus, with each class's own shares kept at the stated F; the
last column is seed-paired.

| cell | seeds | B after curves, own shares | D with B's prior, own shares | paired difference |
|---|---|---|---|---|
| period 1, 1 samples, 10 reads, pcr | 6 | -0.182 ± 0.015 | -0.038 ± 0.007 | +0.144 ± 0.014 |
| period 1, 1 samples, 3 reads, pcr | 6 | -0.141 ± 0.012 | -0.039 ± 0.004 | +0.102 ± 0.010 |
| period 1, 1 samples, 3 reads, pcr_free | 3 | -0.144 ± 0.006 | -0.036 ± 0.001 | +0.107 ± 0.005 |
| period 1, 1 samples, 30 reads, pcr | 6 | -0.139 ± 0.011 | -0.032 ± 0.003 | +0.106 ± 0.011 |
| period 1, 1 samples, 300 reads, pcr | 6 | -0.157 ± 0.006 | -0.042 ± 0.002 | +0.115 ± 0.007 |
| period 1, 1 samples, 300 reads, pcr_free | 3 | -0.148 ± 0.001 | -0.039 ± 0.002 | +0.110 ± 0.002 |
| period 2, 1 samples, 10 reads, pcr | 6 | -0.196 ± 0.017 | -0.045 ± 0.005 | +0.151 ± 0.013 |
| period 2, 1 samples, 3 reads, pcr | 6 | -0.186 ± 0.010 | -0.053 ± 0.002 | +0.133 ± 0.008 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | 3 | -0.169 ± 0.006 | -0.051 ± 0.003 | +0.118 ± 0.009 |
| period 2, 1 samples, 3 reads, pcr_free | 3 | -0.209 ± 0.009 | -0.048 ± 0.005 | +0.162 ± 0.004 |
| period 2, 1 samples, 30 reads, pcr | 6 | -0.192 ± 0.006 | -0.049 ± 0.003 | +0.143 ± 0.009 |
| period 2, 1 samples, 300 reads, pcr | 6 | -0.196 ± 0.011 | -0.055 ± 0.003 | +0.141 ± 0.011 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | 3 | -0.183 ± 0.003 | -0.055 ± 0.002 | +0.128 ± 0.002 |
| period 2, 1 samples, 300 reads, pcr_free | 3 | -0.176 ± 0.022 | -0.045 ± 0.003 | +0.131 ± 0.019 |

## Appendix B: the full tables

Generated by `scripts/research/str_fit_summarise.py tmp/grid tmp/grid_b tmp/grid_dtp --compare
B_raw,B_curves,D,D_tract_prior`. `—` where a method was not run in that cell or gave no estimate. In
the first table the one-sample rows are at calling's fallback; Appendix A has them with own shares.

### Held-out prediction, side by side

Log-likelihood a held-out locus under each fit minus under the truth (zero is as good as the truth; below zero is worse), at calling's fallback for a class without its own numbers. Mean over the seeds ± its standard error. The last columns are the paired difference, seed by seed, of D against B: above zero means D predicts better.

| cell | seeds | loci a seed | B_raw | B_curves | D | D_tract_prior | D − B_raw | D − B_curves | D_tract_prior − B_curves |
|---|---|---|---|---|---|---|---|---|---|
| period 1, 1 samples, 10 reads, pcr | 6 | 971 | -1.087 ± 0.088 | -1.092 ± 0.088 | -1.650 ± 0.015 | -1.650 ± 0.015 | -0.564 ± 0.082 | -0.559 ± 0.082 | -0.559 ± 0.082 |
| period 1, 1 samples, 3 reads, pcr | 6 | 921 | -1.537 ± 0.043 | -1.524 ± 0.042 | -1.628 ± 0.021 | -1.628 ± 0.021 | -0.091 ± 0.036 | -0.104 ± 0.037 | -0.104 ± 0.037 |
| period 1, 1 samples, 3 reads, pcr_free | 3 | 921 | -1.510 ± 0.109 | -1.489 ± 0.113 | -1.657 ± 0.047 | -1.657 ± 0.047 | -0.148 ± 0.103 | -0.168 ± 0.108 | -0.168 ± 0.108 |
| period 1, 1 samples, 30 reads, pcr | 6 | 971 | -0.887 ± 0.083 | -0.868 ± 0.082 | -1.695 ± 0.020 | -1.695 ± 0.020 | -0.807 ± 0.073 | -0.826 ± 0.073 | -0.826 ± 0.073 |
| period 1, 1 samples, 300 reads, pcr | 6 | 971 | -0.954 ± 0.059 | -0.942 ± 0.058 | -1.653 ± 0.027 | -1.653 ± 0.027 | -0.699 ± 0.072 | -0.712 ± 0.072 | -0.712 ± 0.072 |
| period 1, 1 samples, 300 reads, pcr_free | 3 | 971 | -1.125 ± 0.145 | -0.975 ± 0.125 | -1.688 ± 0.036 | -1.688 ± 0.036 | -0.564 ± 0.177 | -0.713 ± 0.156 | -0.713 ± 0.156 |
| period 1, 10 samples, 10 reads, pcr | 3 | 971 | -0.219 ± 0.005 | -0.207 ± 0.006 | -0.191 ± 0.017 | -0.119 ± 0.011 | +0.028 ± 0.014 | +0.016 ± 0.012 | +0.088 ± 0.008 |
| period 1, 10 samples, 3 reads, pcr | 3 | 971 | -0.196 ± 0.020 | -0.186 ± 0.022 | -0.141 ± 0.007 | -0.087 ± 0.009 | +0.055 ± 0.013 | +0.045 ± 0.015 | +0.100 ± 0.014 |
| period 1, 10 samples, 30 reads, pcr | 3 | 971 | -0.221 ± 0.015 | -0.211 ± 0.014 | -0.186 ± 0.018 | -0.116 ± 0.007 | +0.035 ± 0.009 | +0.025 ± 0.009 | +0.095 ± 0.009 |
| period 1, 10 samples, 300 reads, pcr | 3 | 971 | -0.210 ± 0.003 | -0.206 ± 0.004 | -0.191 ± 0.013 | -0.124 ± 0.006 | +0.019 ± 0.013 | +0.015 ± 0.014 | +0.082 ± 0.006 |
| period 1, 2000 samples, 3 reads, pcr | 3 | 276 | -0.549 ± 0.019 | -0.527 ± 0.012 | -1.157 | -0.518 ± 0.014 | -0.595 | -0.619 | +0.009 ± 0.003 |
| period 1, 2000 samples, 3 reads, pcr_free | 1 | 276 | -0.568 | -0.535 | — | -0.526 | — | — | +0.009 |
| period 1, 500 samples, 10 reads, pcr | 3 | 276 | -0.458 ± 0.036 | -0.461 ± 0.035 | -0.959 ± 0.045 | -0.435 ± 0.039 | -0.501 ± 0.013 | -0.498 ± 0.019 | +0.026 ± 0.006 |
| period 1, 500 samples, 3 reads, pcr | 3 | 276 | -0.441 ± 0.022 | -0.427 ± 0.016 | -0.989 ± 0.053 | -0.417 ± 0.017 | -0.548 ± 0.035 | -0.562 ± 0.038 | +0.010 ± 0.006 |
| period 1, 500 samples, 3 reads, pcr_free | 3 | 276 | -0.467 ± 0.019 | -0.445 ± 0.018 | — | -0.434 ± 0.017 | — | — | +0.012 ± 0.003 |
| period 1, 60 samples, 10 reads, pcr | 3 | 971 | -0.277 ± 0.025 | -0.264 ± 0.022 | -0.421 ± 0.011 | -0.237 ± 0.019 | -0.144 ± 0.020 | -0.157 ± 0.018 | +0.027 ± 0.004 |
| period 1, 60 samples, 10 reads, pcr_free | 3 | 971 | -0.279 ± 0.024 | -0.259 ± 0.024 | -0.399 ± 0.008 | -0.231 ± 0.019 | -0.120 ± 0.016 | -0.139 ± 0.016 | +0.029 ± 0.008 |
| period 1, 60 samples, 3 reads, pcr | 3 | 971 | -0.239 ± 0.004 | -0.227 ± 0.003 | -0.344 ± 0.010 | -0.193 ± 0.002 | -0.105 ± 0.014 | -0.117 ± 0.013 | +0.034 ± 0.002 |
| period 1, 60 samples, 30 reads, pcr | 3 | 971 | -0.250 ± 0.010 | -0.242 ± 0.010 | -0.391 ± 0.006 | -0.216 ± 0.007 | -0.141 ± 0.016 | -0.149 ± 0.016 | +0.026 ± 0.004 |
| period 1, 60 samples, 300 reads, pcr | 3 | 971 | -0.258 ± 0.038 | -0.256 ± 0.040 | -0.413 ± 0.047 | -0.225 ± 0.038 | -0.154 ± 0.010 | -0.157 ± 0.009 | +0.031 ± 0.003 |
| period 2, 1 samples, 10 reads, pcr | 6 | 878 | -1.614 ± 0.101 | -1.594 ± 0.101 | -1.903 ± 0.021 | -1.903 ± 0.021 | -0.289 ± 0.091 | -0.309 ± 0.090 | -0.309 ± 0.090 |
| period 2, 1 samples, 3 reads, pcr | 6 | 833 | -1.854 ± 0.025 | -1.802 ± 0.023 | -1.828 ± 0.020 | -1.828 ± 0.020 | +0.026 ± 0.022 | -0.026 ± 0.021 | -0.026 ± 0.021 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | 3 | 829 | -1.822 ± 0.059 | -1.778 ± 0.054 | -1.863 ± 0.020 | -1.863 ± 0.020 | -0.042 ± 0.076 | -0.085 ± 0.072 | -0.085 ± 0.072 |
| period 2, 1 samples, 3 reads, pcr_free | 3 | 830 | -1.980 ± 0.381 | -1.697 ± 0.110 | -1.880 ± 0.045 | -1.880 ± 0.046 | +0.100 ± 0.368 | -0.183 ± 0.085 | -0.183 ± 0.085 |
| period 2, 1 samples, 30 reads, pcr | 6 | 878 | -1.447 ± 0.114 | -1.412 ± 0.120 | -1.888 ± 0.026 | -1.888 ± 0.026 | -0.442 ± 0.121 | -0.476 ± 0.129 | -0.476 ± 0.129 |
| period 2, 1 samples, 300 reads, pcr | 6 | 878 | -1.269 ± 0.115 | -1.252 ± 0.116 | -1.885 ± 0.013 | -1.885 ± 0.013 | -0.615 ± 0.121 | -0.633 ± 0.122 | -0.633 ± 0.122 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | 3 | 878 | -1.231 ± 0.043 | -1.219 ± 0.041 | -1.912 ± 0.028 | -1.912 ± 0.028 | -0.681 ± 0.017 | -0.693 ± 0.014 | -0.693 ± 0.014 |
| period 2, 1 samples, 300 reads, pcr_free | 3 | 878 | -1.605 ± 0.356 | -1.273 ± 0.149 | -1.856 ± 0.015 | -1.856 ± 0.015 | -0.251 ± 0.372 | -0.584 ± 0.164 | -0.584 ± 0.164 |
| period 2, 10 samples, 10 reads, pcr | 3 | 878 | -0.320 ± 0.018 | -0.299 ± 0.021 | -0.274 ± 0.010 | -0.152 ± 0.005 | +0.046 ± 0.010 | +0.025 ± 0.013 | +0.147 ± 0.016 |
| period 2, 10 samples, 3 reads, pcr | 3 | 878 | -0.288 ± 0.021 | -0.255 ± 0.020 | -0.193 ± 0.016 | -0.106 ± 0.009 | +0.095 ± 0.013 | +0.061 ± 0.010 | +0.149 ± 0.012 |
| period 2, 10 samples, 30 reads, pcr | 3 | 878 | -0.279 ± 0.022 | -0.265 ± 0.021 | -0.263 ± 0.022 | -0.160 ± 0.017 | +0.015 ± 0.016 | +0.001 ± 0.019 | +0.105 ± 0.012 |
| period 2, 10 samples, 300 reads, pcr | 3 | 878 | -0.284 ± 0.026 | -0.271 ± 0.026 | -0.268 ± 0.012 | -0.158 ± 0.008 | +0.016 ± 0.019 | +0.003 ± 0.019 | +0.113 ± 0.018 |
| period 2, 2000 samples, 3 reads, pcr | 3 | 388 | -0.728 ± 0.033 | -0.697 ± 0.029 | -1.889 ± 0.067 | -0.645 ± 0.025 | -1.146 ± 0.117 | -1.169 ± 0.096 | +0.051 ± 0.008 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | 1 | 388 | -0.866 | -0.840 | — | -0.800 | — | — | +0.040 |
| period 2, 2000 samples, 3 reads, pcr_free | 1 | 388 | -0.698 | -0.651 | — | -0.627 | — | — | +0.024 |
| period 2, 500 samples, 10 reads, pcr | 3 | 388 | -0.661 ± 0.009 | -0.650 ± 0.010 | -1.568 ± 0.054 | -0.608 ± 0.011 | -0.907 ± 0.046 | -0.919 ± 0.048 | +0.042 ± 0.010 |
| period 2, 500 samples, 3 reads, pcr | 3 | 388 | -0.575 ± 0.022 | -0.543 ± 0.018 | -1.374 ± 0.044 | -0.509 ± 0.011 | -0.799 ± 0.047 | -0.831 ± 0.052 | +0.034 ± 0.008 |
| period 2, 500 samples, 3 reads, pcr_F0.4 | 3 | 388 | -0.722 ± 0.038 | -0.687 ± 0.029 | — | -0.643 ± 0.027 | — | — | +0.044 ± 0.010 |
| period 2, 500 samples, 3 reads, pcr_free | 3 | 388 | -0.627 ± 0.015 | -0.581 ± 0.018 | — | -0.553 ± 0.018 | — | — | +0.029 ± 0.006 |
| period 2, 60 samples, 10 reads, pcr | 3 | 878 | -0.348 ± 0.030 | -0.332 ± 0.029 | -0.590 ± 0.027 | -0.292 ± 0.023 | -0.242 ± 0.011 | -0.258 ± 0.010 | +0.040 ± 0.006 |
| period 2, 60 samples, 10 reads, pcr_F0.4 | 3 | 878 | -0.391 ± 0.005 | -0.375 ± 0.007 | -0.594 ± 0.006 | -0.320 ± 0.005 | -0.203 ± 0.010 | -0.219 ± 0.011 | +0.055 ± 0.003 |
| period 2, 60 samples, 10 reads, pcr_free | 3 | 878 | -0.352 ± 0.031 | -0.326 ± 0.029 | -0.582 ± 0.031 | -0.286 ± 0.027 | -0.231 ± 0.007 | -0.256 ± 0.009 | +0.040 ± 0.002 |
| period 2, 60 samples, 3 reads, pcr | 3 | 878 | -0.321 ± 0.015 | -0.291 ± 0.018 | -0.502 ± 0.045 | -0.252 ± 0.016 | -0.181 ± 0.032 | -0.211 ± 0.030 | +0.039 ± 0.005 |
| period 2, 60 samples, 30 reads, pcr | 3 | 878 | -0.390 ± 0.043 | -0.377 ± 0.040 | -0.605 ± 0.021 | -0.295 ± 0.014 | -0.215 ± 0.028 | -0.229 ± 0.026 | +0.081 ± 0.031 |
| period 2, 60 samples, 300 reads, pcr | 3 | 878 | -0.382 ± 0.029 | -0.380 ± 0.029 | -0.641 ± 0.030 | -0.326 ± 0.025 | -0.259 ± 0.017 | -0.261 ± 0.017 | +0.054 ± 0.010 |

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
| period 1, 1 samples, 10 reads, pcr | D_tract_prior | 0.0616 | 0.0198 | 0.0224 | — | 0.0361 |
| period 1, 1 samples, 3 reads, pcr | B_raw | 0.2175 | 0.1586 | 0.0962 | 0.2415 | 0.0602 |
| period 1, 1 samples, 3 reads, pcr | B_curves | 0.1757 | 0.0943 | 0.0602 | 0.2415 | 0.0602 |
| period 1, 1 samples, 3 reads, pcr | D | 0.0873 | 0.0752 | 0.0450 | — | 0.0420 |
| period 1, 1 samples, 3 reads, pcr | D_tract_prior | 0.0872 | 0.0738 | 0.0454 | — | 0.0430 |
| period 1, 1 samples, 3 reads, pcr_free | B_raw | 0.6567 | 0.3692 | 0.1973 | 0.4731 | 0.0630 |
| period 1, 1 samples, 3 reads, pcr_free | B_curves | 0.6696 | 0.3186 | 0.1282 | 0.4731 | 0.0630 |
| period 1, 1 samples, 3 reads, pcr_free | D | 0.3189 | 0.1928 | 0.0951 | — | 0.0379 |
| period 1, 1 samples, 3 reads, pcr_free | D_tract_prior | 0.3191 | 0.1926 | 0.0953 | — | 0.0422 |
| period 1, 1 samples, 30 reads, pcr | B_raw | 0.0599 | 0.0316 | 0.0196 | 0.1803 | 0.0567 |
| period 1, 1 samples, 30 reads, pcr | B_curves | 0.0617 | 0.0174 | 0.0112 | 0.1803 | 0.0567 |
| period 1, 1 samples, 30 reads, pcr | D | 0.0398 | 0.0089 | 0.0107 | — | 0.0444 |
| period 1, 1 samples, 30 reads, pcr | D_tract_prior | 0.0398 | 0.0089 | 0.0108 | — | 0.0451 |
| period 1, 1 samples, 300 reads, pcr | B_raw | 0.0209 | 0.0082 | 0.0079 | 0.3068 | 0.0518 |
| period 1, 1 samples, 300 reads, pcr | B_curves | 0.0204 | 0.0065 | 0.0064 | 0.3068 | 0.0518 |
| period 1, 1 samples, 300 reads, pcr | D | 0.0195 | 0.0054 | 0.0060 | — | 0.0372 |
| period 1, 1 samples, 300 reads, pcr | D_tract_prior | 0.0195 | 0.0054 | 0.0060 | — | 0.0408 |
| period 1, 1 samples, 300 reads, pcr_free | B_raw | 0.0753 | 0.0327 | 0.0224 | 0.3047 | 0.0583 |
| period 1, 1 samples, 300 reads, pcr_free | B_curves | 0.0766 | 0.0185 | 0.0114 | 0.3047 | 0.0583 |
| period 1, 1 samples, 300 reads, pcr_free | D | 0.0549 | 0.0195 | 0.0106 | — | 0.0336 |
| period 1, 1 samples, 300 reads, pcr_free | D_tract_prior | 0.0549 | 0.0194 | 0.0106 | — | 0.0397 |
| period 1, 10 samples, 10 reads, pcr | B_raw | 0.0233 | 0.0139 | 0.0102 | 0.0523 | 0.0379 |
| period 1, 10 samples, 10 reads, pcr | B_curves | 0.0246 | 0.0093 | 0.0060 | 0.0523 | 0.0379 |
| period 1, 10 samples, 10 reads, pcr | D | 0.0210 | 0.0091 | 0.0058 | 0.0463 | 0.0307 |
| period 1, 10 samples, 10 reads, pcr | D_tract_prior | 0.0211 | 0.0093 | 0.0059 | 0.0443 | 0.0307 |
| period 1, 10 samples, 3 reads, pcr | B_raw | 0.0728 | 0.0335 | 0.0293 | 0.0425 | 0.0398 |
| period 1, 10 samples, 3 reads, pcr | B_curves | 0.0651 | 0.0271 | 0.0229 | 0.0425 | 0.0398 |
| period 1, 10 samples, 3 reads, pcr | D | 0.0575 | 0.0256 | 0.0183 | 0.0314 | 0.0321 |
| period 1, 10 samples, 3 reads, pcr | D_tract_prior | 0.0561 | 0.0251 | 0.0180 | 0.0319 | 0.0328 |
| period 1, 10 samples, 30 reads, pcr | B_raw | 0.0192 | 0.0073 | 0.0100 | 0.0436 | 0.0398 |
| period 1, 10 samples, 30 reads, pcr | B_curves | 0.0217 | 0.0043 | 0.0069 | 0.0436 | 0.0398 |
| period 1, 10 samples, 30 reads, pcr | D | 0.0173 | 0.0044 | 0.0084 | 0.0352 | 0.0365 |
| period 1, 10 samples, 30 reads, pcr | D_tract_prior | 0.0173 | 0.0043 | 0.0084 | 0.0350 | 0.0321 |
| period 1, 10 samples, 300 reads, pcr | B_raw | 0.0065 | 0.0028 | 0.0024 | 0.0479 | 0.0428 |
| period 1, 10 samples, 300 reads, pcr | B_curves | 0.0068 | 0.0020 | 0.0021 | 0.0479 | 0.0428 |
| period 1, 10 samples, 300 reads, pcr | D | 0.0063 | 0.0021 | 0.0020 | 0.0332 | 0.0411 |
| period 1, 10 samples, 300 reads, pcr | D_tract_prior | 0.0063 | 0.0021 | 0.0020 | 0.0310 | 0.0376 |
| period 1, 2000 samples, 10 reads, pcr | B_raw | 0.0024 | 0.0019 | 0.0012 | 0.0255 | 0.0619 |
| period 1, 2000 samples, 10 reads, pcr | B_curves | 0.0024 | 0.0019 | 0.0014 | 0.0255 | 0.0619 |
| period 1, 2000 samples, 10 reads, pcr | D | 0.0019 | 0.0012 | 0.0012 | 0.0292 | 0.0567 |
| period 1, 2000 samples, 3 reads, pcr | B_raw | 0.0076 | 0.0040 | 0.0023 | 0.0357 | 0.0502 |
| period 1, 2000 samples, 3 reads, pcr | B_curves | 0.0075 | 0.0030 | 0.0024 | 0.0357 | 0.0502 |
| period 1, 2000 samples, 3 reads, pcr | D | 0.0063 | 0.0025 | 0.0015 | 0.0254 | 0.0552 |
| period 1, 2000 samples, 3 reads, pcr | D_tract_prior | 0.0076 | 0.0028 | 0.0019 | 0.0359 | 0.0502 |
| period 1, 2000 samples, 3 reads, pcr_free | B_raw | 0.0176 | 0.0064 | 0.0082 | 0.0265 | 0.0313 |
| period 1, 2000 samples, 3 reads, pcr_free | B_curves | 0.0167 | 0.0048 | 0.0076 | 0.0265 | 0.0313 |
| period 1, 2000 samples, 3 reads, pcr_free | D | 0.0104 | 0.0032 | 0.0061 | 0.0293 | 0.0293 |
| period 1, 2000 samples, 3 reads, pcr_free | D_tract_prior | 0.0169 | 0.0050 | 0.0062 | 0.0262 | 0.0313 |
| period 1, 500 samples, 10 reads, pcr | B_raw | 0.0061 | 0.0024 | 0.0025 | 0.0310 | 0.0447 |
| period 1, 500 samples, 10 reads, pcr | B_curves | 0.0061 | 0.0019 | 0.0026 | 0.0310 | 0.0447 |
| period 1, 500 samples, 10 reads, pcr | D | 0.0060 | 0.0013 | 0.0016 | 0.0327 | 0.0412 |
| period 1, 500 samples, 10 reads, pcr | D_tract_prior | 0.0061 | 0.0015 | 0.0016 | 0.0307 | 0.0448 |
| period 1, 500 samples, 3 reads, pcr | B_raw | 0.0128 | 0.0066 | 0.0031 | 0.0340 | 0.0475 |
| period 1, 500 samples, 3 reads, pcr | B_curves | 0.0128 | 0.0048 | 0.0024 | 0.0340 | 0.0475 |
| period 1, 500 samples, 3 reads, pcr | D | 0.0118 | 0.0046 | 0.0020 | 0.0361 | 0.0501 |
| period 1, 500 samples, 3 reads, pcr | D_tract_prior | 0.0122 | 0.0047 | 0.0022 | 0.0337 | 0.0473 |
| period 1, 500 samples, 3 reads, pcr_free | B_raw | 0.0274 | 0.0187 | 0.0097 | 0.0294 | 0.0476 |
| period 1, 500 samples, 3 reads, pcr_free | B_curves | 0.0272 | 0.0088 | 0.0071 | 0.0294 | 0.0476 |
| period 1, 500 samples, 3 reads, pcr_free | D_tract_prior | 0.0264 | 0.0070 | 0.0064 | 0.0290 | 0.0474 |
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
| period 1, 60 samples, 10 reads, pcr_free | D_tract_prior | 0.0408 | 0.0075 | 0.0100 | 0.0381 | 0.0425 |
| period 1, 60 samples, 3 reads, pcr | B_raw | 0.0271 | 0.0155 | 0.0107 | 0.0457 | 0.0524 |
| period 1, 60 samples, 3 reads, pcr | B_curves | 0.0185 | 0.0108 | 0.0072 | 0.0457 | 0.0524 |
| period 1, 60 samples, 3 reads, pcr | D | 0.0219 | 0.0095 | 0.0062 | 0.0494 | 0.0461 |
| period 1, 60 samples, 3 reads, pcr | D_tract_prior | 0.0263 | 0.0092 | 0.0068 | 0.0445 | 0.0526 |
| period 1, 60 samples, 30 reads, pcr | B_raw | 0.0080 | 0.0040 | 0.0027 | 0.0406 | 0.0325 |
| period 1, 60 samples, 30 reads, pcr | B_curves | 0.0080 | 0.0029 | 0.0026 | 0.0406 | 0.0325 |
| period 1, 60 samples, 30 reads, pcr | D | 0.0076 | 0.0027 | 0.0022 | 0.0409 | 0.0362 |
| period 1, 60 samples, 30 reads, pcr | D_tract_prior | 0.0074 | 0.0027 | 0.0022 | 0.0394 | 0.0313 |
| period 1, 60 samples, 300 reads, pcr | B_raw | 0.0033 | 0.0014 | 0.0015 | 0.0401 | 0.0408 |
| period 1, 60 samples, 300 reads, pcr | B_curves | 0.0033 | 0.0013 | 0.0018 | 0.0401 | 0.0408 |
| period 1, 60 samples, 300 reads, pcr | D | 0.0033 | 0.0010 | 0.0015 | 0.0402 | 0.0412 |
| period 1, 60 samples, 300 reads, pcr | D_tract_prior | 0.0033 | 0.0010 | 0.0015 | 0.0394 | 0.0385 |
| period 2, 1 samples, 10 reads, pcr | B_raw | 0.1330 | 0.0584 | 0.0579 | 0.3641 | 0.0789 |
| period 2, 1 samples, 10 reads, pcr | B_curves | 0.0985 | 0.0347 | 0.0269 | 0.3641 | 0.0789 |
| period 2, 1 samples, 10 reads, pcr | D | 0.0687 | 0.0300 | 0.0225 | — | 0.0456 |
| period 2, 1 samples, 10 reads, pcr | D_tract_prior | 0.0685 | 0.0297 | 0.0224 | — | 0.0504 |
| period 2, 1 samples, 3 reads, pcr | B_raw | 0.2187 | 0.1163 | 0.1157 | 0.4809 | 0.0873 |
| period 2, 1 samples, 3 reads, pcr | B_curves | 0.1527 | 0.0785 | 0.0422 | 0.4809 | 0.0873 |
| period 2, 1 samples, 3 reads, pcr | D | 0.1327 | 0.0538 | 0.0393 | — | 0.0485 |
| period 2, 1 samples, 3 reads, pcr | D_tract_prior | 0.1319 | 0.0540 | 0.0397 | — | 0.0477 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_raw | 0.2716 | 0.1289 | 0.1030 | 0.0979 | 0.0895 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_curves | 0.1499 | 0.0759 | 0.0397 | 0.0979 | 0.0895 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D | 0.1029 | 0.0364 | 0.0343 | — | 0.0453 |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D_tract_prior | 0.1049 | 0.0365 | 0.0329 | — | 0.0492 |
| period 2, 1 samples, 3 reads, pcr_free | B_raw | 0.7556 | 0.3023 | 0.2497 | 0.5441 | 0.0714 |
| period 2, 1 samples, 3 reads, pcr_free | B_curves | 0.6303 | 0.2286 | 0.2542 | 0.5441 | 0.0714 |
| period 2, 1 samples, 3 reads, pcr_free | D | 0.3403 | 0.1690 | 0.1334 | — | 0.0435 |
| period 2, 1 samples, 3 reads, pcr_free | D_tract_prior | 0.3385 | 0.1687 | 0.1333 | — | 0.0494 |
| period 2, 1 samples, 30 reads, pcr | B_raw | 0.0879 | 0.0331 | 0.0340 | 0.3493 | 0.0817 |
| period 2, 1 samples, 30 reads, pcr | B_curves | 0.0763 | 0.0166 | 0.0199 | 0.3493 | 0.0817 |
| period 2, 1 samples, 30 reads, pcr | D | 0.0687 | 0.0105 | 0.0116 | — | 0.0421 |
| period 2, 1 samples, 30 reads, pcr | D_tract_prior | 0.0686 | 0.0106 | 0.0116 | — | 0.0456 |
| period 2, 1 samples, 300 reads, pcr | B_raw | 0.0223 | 0.0114 | 0.0100 | 0.3132 | 0.0834 |
| period 2, 1 samples, 300 reads, pcr | B_curves | 0.0220 | 0.0068 | 0.0056 | 0.3132 | 0.0834 |
| period 2, 1 samples, 300 reads, pcr | D | 0.0191 | 0.0056 | 0.0036 | — | 0.0426 |
| period 2, 1 samples, 300 reads, pcr | D_tract_prior | 0.0191 | 0.0056 | 0.0036 | — | 0.0462 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_raw | 0.0217 | 0.0128 | 0.0108 | 0.2238 | 0.0699 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_curves | 0.0203 | 0.0066 | 0.0045 | 0.2238 | 0.0699 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D | 0.0196 | 0.0046 | 0.0030 | — | 0.0430 |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D_tract_prior | 0.0196 | 0.0046 | 0.0030 | — | 0.0404 |
| period 2, 1 samples, 300 reads, pcr_free | B_raw | 0.0770 | 0.0390 | 0.0274 | 0.3283 | 0.0542 |
| period 2, 1 samples, 300 reads, pcr_free | B_curves | 0.0664 | 0.0156 | 0.0186 | 0.3283 | 0.0542 |
| period 2, 1 samples, 300 reads, pcr_free | D | 0.0518 | 0.0116 | 0.0140 | — | 0.0399 |
| period 2, 1 samples, 300 reads, pcr_free | D_tract_prior | 0.0515 | 0.0116 | 0.0140 | — | 0.0523 |
| period 2, 10 samples, 10 reads, pcr | B_raw | 0.0338 | 0.0173 | 0.0151 | 0.0590 | 0.0897 |
| period 2, 10 samples, 10 reads, pcr | B_curves | 0.0301 | 0.0115 | 0.0072 | 0.0590 | 0.0897 |
| period 2, 10 samples, 10 reads, pcr | D | 0.0303 | 0.0051 | 0.0074 | 0.0219 | 0.0654 |
| period 2, 10 samples, 10 reads, pcr | D_tract_prior | 0.0305 | 0.0050 | 0.0073 | 0.0254 | 0.0582 |
| period 2, 10 samples, 3 reads, pcr | B_raw | 0.0828 | 0.0461 | 0.0360 | 0.0592 | 0.0643 |
| period 2, 10 samples, 3 reads, pcr | B_curves | 0.0782 | 0.0267 | 0.0159 | 0.0592 | 0.0643 |
| period 2, 10 samples, 3 reads, pcr | D | 0.0449 | 0.0253 | 0.0110 | 0.0224 | 0.0522 |
| period 2, 10 samples, 3 reads, pcr | D_tract_prior | 0.0469 | 0.0254 | 0.0103 | 0.0205 | 0.0544 |
| period 2, 10 samples, 30 reads, pcr | B_raw | 0.0252 | 0.0116 | 0.0111 | 0.0711 | 0.0824 |
| period 2, 10 samples, 30 reads, pcr | B_curves | 0.0219 | 0.0058 | 0.0051 | 0.0711 | 0.0824 |
| period 2, 10 samples, 30 reads, pcr | D | 0.0184 | 0.0056 | 0.0055 | 0.0467 | 0.0582 |
| period 2, 10 samples, 30 reads, pcr | D_tract_prior | 0.0185 | 0.0055 | 0.0055 | 0.0482 | 0.0564 |
| period 2, 10 samples, 300 reads, pcr | B_raw | 0.0070 | 0.0036 | 0.0041 | 0.0522 | 0.0628 |
| period 2, 10 samples, 300 reads, pcr | B_curves | 0.0071 | 0.0030 | 0.0025 | 0.0522 | 0.0628 |
| period 2, 10 samples, 300 reads, pcr | D | 0.0078 | 0.0016 | 0.0021 | 0.0299 | 0.0601 |
| period 2, 10 samples, 300 reads, pcr | D_tract_prior | 0.0078 | 0.0016 | 0.0021 | 0.0320 | 0.0494 |
| period 2, 2000 samples, 10 reads, pcr | B_raw | 0.0035 | 0.0013 | 0.0012 | 0.0596 | 0.0996 |
| period 2, 2000 samples, 10 reads, pcr | B_curves | 0.0035 | 0.0013 | 0.0013 | 0.0596 | 0.0996 |
| period 2, 2000 samples, 10 reads, pcr | D | 0.0030 | 0.0012 | 0.0012 | 0.0608 | 0.1034 |
| period 2, 2000 samples, 3 reads, pcr | B_raw | 0.0089 | 0.0040 | 0.0033 | 0.0539 | 0.0673 |
| period 2, 2000 samples, 3 reads, pcr | B_curves | 0.0087 | 0.0038 | 0.0030 | 0.0539 | 0.0673 |
| period 2, 2000 samples, 3 reads, pcr | D | 0.0042 | 0.0019 | 0.0018 | 0.0554 | 0.0767 |
| period 2, 2000 samples, 3 reads, pcr | D_tract_prior | 0.0089 | 0.0024 | 0.0030 | 0.0540 | 0.0672 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_raw | 0.0065 | 0.0018 | 0.0030 | 0.0245 | 0.0342 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_curves | 0.0062 | 0.0020 | 0.0026 | 0.0245 | 0.0342 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | D | 0.0029 | 0.0010 | 0.0019 | 0.0351 | 0.0330 |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | D_tract_prior | 0.0046 | 0.0018 | 0.0025 | 0.0243 | 0.0342 |
| period 2, 2000 samples, 3 reads, pcr_free | B_raw | 0.0164 | 0.0038 | 0.0109 | 0.0536 | 0.0414 |
| period 2, 2000 samples, 3 reads, pcr_free | B_curves | 0.0140 | 0.0047 | 0.0076 | 0.0536 | 0.0414 |
| period 2, 2000 samples, 3 reads, pcr_free | D | 0.0110 | 0.0048 | 0.0081 | 0.0511 | 0.0526 |
| period 2, 2000 samples, 3 reads, pcr_free | D_tract_prior | 0.0132 | 0.0042 | 0.0055 | 0.0537 | 0.0415 |
| period 2, 500 samples, 10 reads, pcr | B_raw | 0.0063 | 0.0028 | 0.0021 | 0.0451 | 0.0696 |
| period 2, 500 samples, 10 reads, pcr | B_curves | 0.0062 | 0.0026 | 0.0016 | 0.0451 | 0.0696 |
| period 2, 500 samples, 10 reads, pcr | D | 0.0058 | 0.0013 | 0.0018 | 0.0528 | 0.0756 |
| period 2, 500 samples, 10 reads, pcr | D_tract_prior | 0.0056 | 0.0013 | 0.0018 | 0.0444 | 0.0693 |
| period 2, 500 samples, 3 reads, pcr | B_raw | 0.0128 | 0.0049 | 0.0055 | 0.0455 | 0.0537 |
| period 2, 500 samples, 3 reads, pcr | B_curves | 0.0131 | 0.0044 | 0.0026 | 0.0455 | 0.0537 |
| period 2, 500 samples, 3 reads, pcr | D | 0.0097 | 0.0028 | 0.0020 | 0.0509 | 0.0508 |
| period 2, 500 samples, 3 reads, pcr | D_tract_prior | 0.0115 | 0.0029 | 0.0026 | 0.0452 | 0.0535 |
| period 2, 500 samples, 3 reads, pcr_F0.4 | B_raw | 0.0158 | 0.0044 | 0.0044 | 0.0527 | 0.0438 |
| period 2, 500 samples, 3 reads, pcr_F0.4 | B_curves | 0.0164 | 0.0031 | 0.0027 | 0.0527 | 0.0438 |
| period 2, 500 samples, 3 reads, pcr_F0.4 | D_tract_prior | 0.0126 | 0.0023 | 0.0032 | 0.0522 | 0.0437 |
| period 2, 500 samples, 3 reads, pcr_free | B_raw | 0.0385 | 0.0172 | 0.0107 | 0.0559 | 0.0537 |
| period 2, 500 samples, 3 reads, pcr_free | B_curves | 0.0427 | 0.0123 | 0.0066 | 0.0559 | 0.0537 |
| period 2, 500 samples, 3 reads, pcr_free | D_tract_prior | 0.0317 | 0.0079 | 0.0077 | 0.0553 | 0.0535 |
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
| period 2, 60 samples, 10 reads, pcr_F0.4 | D_tract_prior | 0.0152 | 0.0031 | 0.0037 | 0.0437 | 0.0400 |
| period 2, 60 samples, 10 reads, pcr_free | B_raw | 0.0487 | 0.0202 | 0.0233 | 0.0543 | 0.0665 |
| period 2, 60 samples, 10 reads, pcr_free | B_curves | 0.0406 | 0.0123 | 0.0119 | 0.0543 | 0.0665 |
| period 2, 60 samples, 10 reads, pcr_free | D | 0.0388 | 0.0084 | 0.0119 | 0.0522 | 0.0630 |
| period 2, 60 samples, 10 reads, pcr_free | D_tract_prior | 0.0380 | 0.0082 | 0.0113 | 0.0519 | 0.0638 |
| period 2, 60 samples, 3 reads, pcr | B_raw | 0.0382 | 0.0156 | 0.0129 | 0.0710 | 0.0696 |
| period 2, 60 samples, 3 reads, pcr | B_curves | 0.0332 | 0.0093 | 0.0066 | 0.0710 | 0.0696 |
| period 2, 60 samples, 3 reads, pcr | D | 0.0248 | 0.0073 | 0.0048 | 0.0633 | 0.0734 |
| period 2, 60 samples, 3 reads, pcr | D_tract_prior | 0.0251 | 0.0068 | 0.0047 | 0.0597 | 0.0652 |
| period 2, 60 samples, 30 reads, pcr | B_raw | 0.0106 | 0.0037 | 0.0046 | 0.0488 | 0.0550 |
| period 2, 60 samples, 30 reads, pcr | B_curves | 0.0104 | 0.0028 | 0.0029 | 0.0488 | 0.0550 |
| period 2, 60 samples, 30 reads, pcr | D | 0.0103 | 0.0017 | 0.0027 | 0.0450 | 0.0602 |
| period 2, 60 samples, 30 reads, pcr | D_tract_prior | 0.0103 | 0.0017 | 0.0027 | 0.0455 | 0.0549 |
| period 2, 60 samples, 300 reads, pcr | B_raw | 0.0030 | 0.0015 | 0.0017 | 0.0460 | 0.0636 |
| period 2, 60 samples, 300 reads, pcr | B_curves | 0.0030 | 0.0015 | 0.0012 | 0.0460 | 0.0636 |
| period 2, 60 samples, 300 reads, pcr | D | 0.0031 | 0.0010 | 0.0014 | 0.0434 | 0.0583 |
| period 2, 60 samples, 300 reads, pcr | D_tract_prior | 0.0031 | 0.0010 | 0.0014 | 0.0451 | 0.0599 |

### Error honesty, side by side

Share of estimates within one / two of their reported errors of the truth, pooled over classes and seeds (honest: about 0.68 / 0.95); in brackets, how many estimates had an error. Only the fits as they come carry errors (the curves replace them).

| cell | method | level | shorter_share | fall_off | between_locus_share | shares ≥ 0.05 | shares 0.01–0.05 | shares < 0.01 |
|---|---|---|---|---|---|---|---|---|
| period 1, 1 samples, 10 reads, pcr | B_raw | 0.68 / 0.94 (84) | 0.64 / 0.95 (83) | 0.67 / 0.95 (79) | 0.95 / 1.00 (22) | 0.68 / 0.93 (294) | 0.50 / 0.60 (402) | 0.23 / 0.33 (900) |
| period 1, 1 samples, 10 reads, pcr | B_curves | — | — | — | 0.95 / 1.00 (22) | 0.68 / 0.93 (294) | 0.50 / 0.60 (402) | 0.23 / 0.33 (900) |
| period 1, 1 samples, 10 reads, pcr | D | 0.73 / 0.98 (84) | 0.76 / 0.94 (84) | 0.36 / 0.85 (84) | — | 0.65 / 0.94 (294) | 0.69 / 0.89 (402) | 0.59 / 0.90 (900) |
| period 1, 1 samples, 10 reads, pcr | D_tract_prior | 0.73 / 0.98 (84) | 0.76 / 0.94 (84) | 0.37 / 0.85 (84) | — | 0.60 / 0.93 (294) | 0.69 / 0.91 (402) | 0.52 / 0.93 (900) |
| period 1, 1 samples, 3 reads, pcr | B_raw | 0.61 / 0.92 (84) | 0.66 / 0.87 (77) | 0.54 / 0.86 (71) | 1.00 / 1.00 (3) | 0.67 / 0.89 (294) | 0.43 / 0.53 (402) | 0.26 / 0.35 (900) |
| period 1, 1 samples, 3 reads, pcr | B_curves | — | — | — | 1.00 / 1.00 (3) | 0.67 / 0.89 (294) | 0.43 / 0.53 (402) | 0.26 / 0.35 (900) |
| period 1, 1 samples, 3 reads, pcr | D | 0.69 / 0.98 (84) | 0.56 / 0.90 (84) | 0.58 / 0.93 (84) | — | 0.67 / 0.91 (294) | 0.70 / 0.89 (402) | 0.57 / 0.88 (900) |
| period 1, 1 samples, 3 reads, pcr | D_tract_prior | 0.70 / 0.98 (84) | 0.55 / 0.90 (84) | 0.57 / 0.94 (84) | — | 0.62 / 0.87 (294) | 0.70 / 0.91 (402) | 0.52 / 0.89 (900) |
| period 1, 1 samples, 3 reads, pcr_free | B_raw | 0.81 / 1.00 (31) | 0.68 / 0.95 (19) | 0.41 / 0.82 (17) | 1.00 / 1.00 (3) | 0.69 / 0.91 (143) | 0.47 / 0.57 (194) | 0.26 / 0.36 (442) |
| period 1, 1 samples, 3 reads, pcr_free | B_curves | — | — | — | 1.00 / 1.00 (3) | 0.69 / 0.91 (143) | 0.47 / 0.57 (194) | 0.26 / 0.36 (442) |
| period 1, 1 samples, 3 reads, pcr_free | D | 0.62 / 0.95 (42) | 0.55 / 0.88 (42) | 0.43 / 0.79 (42) | — | 0.72 / 0.93 (147) | 0.78 / 0.94 (201) | 0.66 / 0.90 (450) |
| period 1, 1 samples, 3 reads, pcr_free | D_tract_prior | 0.62 / 0.95 (42) | 0.55 / 0.88 (42) | 0.43 / 0.79 (42) | — | 0.66 / 0.87 (147) | 0.77 / 0.94 (201) | 0.52 / 0.90 (450) |
| period 1, 1 samples, 30 reads, pcr | B_raw | 0.61 / 0.95 (84) | 0.63 / 0.93 (84) | 0.76 / 0.98 (84) | 0.98 / 1.00 (43) | 0.64 / 0.91 (294) | 0.46 / 0.56 (402) | 0.23 / 0.32 (900) |
| period 1, 1 samples, 30 reads, pcr | B_curves | — | — | — | 0.98 / 1.00 (43) | 0.64 / 0.91 (294) | 0.46 / 0.56 (402) | 0.23 / 0.32 (900) |
| period 1, 1 samples, 30 reads, pcr | D | 0.62 / 0.98 (84) | 0.82 / 1.00 (84) | 0.67 / 0.94 (84) | — | 0.62 / 0.86 (294) | 0.59 / 0.86 (402) | 0.60 / 0.87 (900) |
| period 1, 1 samples, 30 reads, pcr | D_tract_prior | 0.63 / 0.98 (84) | 0.83 / 1.00 (84) | 0.67 / 0.94 (84) | — | 0.61 / 0.85 (294) | 0.61 / 0.88 (402) | 0.53 / 0.91 (900) |
| period 1, 1 samples, 300 reads, pcr | B_raw | 0.63 / 0.90 (84) | 0.70 / 0.98 (84) | 0.62 / 0.93 (84) | 0.95 / 1.00 (37) | 0.65 / 0.93 (294) | 0.51 / 0.60 (402) | 0.24 / 0.32 (900) |
| period 1, 1 samples, 300 reads, pcr | B_curves | — | — | — | 0.95 / 1.00 (37) | 0.65 / 0.93 (294) | 0.51 / 0.60 (402) | 0.24 / 0.32 (900) |
| period 1, 1 samples, 300 reads, pcr | D | 0.67 / 0.92 (84) | 0.57 / 0.85 (84) | 0.51 / 0.90 (84) | — | 0.69 / 0.89 (294) | 0.67 / 0.93 (402) | 0.63 / 0.90 (900) |
| period 1, 1 samples, 300 reads, pcr | D_tract_prior | 0.67 / 0.92 (84) | 0.57 / 0.85 (84) | 0.51 / 0.90 (84) | — | 0.66 / 0.87 (294) | 0.69 / 0.95 (402) | 0.55 / 0.94 (900) |
| period 1, 1 samples, 300 reads, pcr_free | B_raw | 0.57 / 0.90 (42) | 0.69 / 0.86 (42) | 0.54 / 0.98 (41) | 0.94 / 1.00 (18) | 0.60 / 0.92 (147) | 0.49 / 0.60 (201) | 0.21 / 0.30 (450) |
| period 1, 1 samples, 300 reads, pcr_free | B_curves | — | — | — | 0.94 / 1.00 (18) | 0.60 / 0.92 (147) | 0.49 / 0.60 (201) | 0.21 / 0.30 (450) |
| period 1, 1 samples, 300 reads, pcr_free | D | 0.55 / 0.90 (42) | 0.40 / 0.67 (42) | 0.76 / 1.00 (42) | — | 0.76 / 0.92 (147) | 0.71 / 0.94 (201) | 0.60 / 0.87 (450) |
| period 1, 1 samples, 300 reads, pcr_free | D_tract_prior | 0.55 / 0.90 (42) | 0.40 / 0.69 (42) | 0.76 / 1.00 (42) | — | 0.72 / 0.90 (147) | 0.73 / 0.95 (201) | 0.54 / 0.91 (450) |
| period 1, 10 samples, 10 reads, pcr | B_raw | 0.86 / 0.95 (42) | 0.69 / 0.98 (42) | 0.76 / 1.00 (42) | 0.69 / 0.83 (42) | 0.64 / 0.90 (147) | 0.51 / 0.69 (201) | 0.35 / 0.46 (450) |
| period 1, 10 samples, 10 reads, pcr | B_curves | — | — | — | 0.69 / 0.83 (42) | 0.64 / 0.90 (147) | 0.51 / 0.69 (201) | 0.35 / 0.46 (450) |
| period 1, 10 samples, 10 reads, pcr | D | 0.76 / 0.98 (42) | 0.57 / 0.83 (42) | 0.67 / 1.00 (42) | 0.69 / 0.83 (42) | 0.63 / 0.90 (147) | 0.51 / 0.69 (201) | 0.44 / 0.55 (450) |
| period 1, 10 samples, 10 reads, pcr | D_tract_prior | 0.76 / 0.98 (42) | 0.57 / 0.83 (42) | 0.67 / 1.00 (42) | 0.71 / 0.81 (42) | 0.64 / 0.91 (147) | 0.56 / 0.74 (201) | 0.55 / 0.78 (450) |
| period 1, 10 samples, 3 reads, pcr | B_raw | 0.55 / 0.98 (42) | 0.76 / 0.98 (42) | 0.62 / 0.88 (40) | 0.69 / 0.95 (42) | 0.66 / 0.93 (147) | 0.54 / 0.66 (201) | 0.38 / 0.52 (450) |
| period 1, 10 samples, 3 reads, pcr | B_curves | — | — | — | 0.69 / 0.95 (42) | 0.66 / 0.93 (147) | 0.54 / 0.66 (201) | 0.38 / 0.52 (450) |
| period 1, 10 samples, 3 reads, pcr | D | 0.67 / 0.95 (42) | 0.43 / 0.76 (42) | 0.45 / 0.69 (42) | 0.76 / 0.95 (42) | 0.67 / 0.92 (147) | 0.55 / 0.72 (201) | 0.44 / 0.58 (450) |
| period 1, 10 samples, 3 reads, pcr | D_tract_prior | 0.71 / 0.95 (42) | 0.45 / 0.76 (42) | 0.48 / 0.69 (42) | 0.71 / 0.95 (42) | 0.67 / 0.93 (147) | 0.63 / 0.79 (201) | 0.63 / 0.85 (450) |
| period 1, 10 samples, 30 reads, pcr | B_raw | 0.67 / 0.90 (42) | 0.76 / 0.98 (42) | 0.50 / 0.88 (42) | 0.64 / 0.95 (42) | 0.65 / 0.89 (147) | 0.50 / 0.68 (201) | 0.38 / 0.51 (450) |
| period 1, 10 samples, 30 reads, pcr | B_curves | — | — | — | 0.64 / 0.95 (42) | 0.65 / 0.89 (147) | 0.50 / 0.68 (201) | 0.38 / 0.51 (450) |
| period 1, 10 samples, 30 reads, pcr | D | 0.69 / 0.90 (42) | 0.60 / 0.83 (42) | 0.57 / 0.79 (42) | 0.57 / 0.79 (42) | 0.68 / 0.88 (147) | 0.54 / 0.70 (201) | 0.45 / 0.60 (450) |
| period 1, 10 samples, 30 reads, pcr | D_tract_prior | 0.69 / 0.90 (42) | 0.60 / 0.81 (42) | 0.57 / 0.79 (42) | 0.60 / 0.79 (42) | 0.66 / 0.89 (147) | 0.59 / 0.75 (201) | 0.63 / 0.84 (450) |
| period 1, 10 samples, 300 reads, pcr | B_raw | 0.62 / 0.95 (42) | 0.76 / 0.93 (42) | 0.62 / 0.95 (42) | 0.67 / 0.93 (42) | 0.61 / 0.90 (147) | 0.54 / 0.69 (201) | 0.29 / 0.41 (450) |
| period 1, 10 samples, 300 reads, pcr | B_curves | — | — | — | 0.67 / 0.93 (42) | 0.61 / 0.90 (147) | 0.54 / 0.69 (201) | 0.29 / 0.41 (450) |
| period 1, 10 samples, 300 reads, pcr | D | 0.60 / 0.93 (42) | 0.67 / 0.95 (42) | 0.69 / 0.95 (42) | 0.64 / 0.95 (42) | 0.61 / 0.89 (147) | 0.55 / 0.72 (201) | 0.34 / 0.44 (450) |
| period 1, 10 samples, 300 reads, pcr | D_tract_prior | 0.60 / 0.93 (42) | 0.67 / 0.95 (42) | 0.69 / 0.95 (42) | 0.67 / 0.93 (42) | 0.62 / 0.90 (147) | 0.60 / 0.76 (201) | 0.46 / 0.66 (450) |
| period 1, 2000 samples, 10 reads, pcr | B_raw | 0.61 / 0.96 (28) | 0.57 / 0.89 (28) | 0.50 / 0.86 (28) | 0.89 / 1.00 (28) | 0.67 / 0.89 (98) | 0.50 / 0.68 (134) | 0.32 / 0.45 (300) |
| period 1, 2000 samples, 10 reads, pcr | B_curves | — | — | — | 0.89 / 1.00 (28) | 0.67 / 0.89 (98) | 0.50 / 0.68 (134) | 0.32 / 0.45 (300) |
| period 1, 2000 samples, 10 reads, pcr | D | 0.71 / 0.96 (28) | 0.61 / 0.93 (28) | 0.64 / 0.86 (28) | 0.82 / 1.00 (28) | 0.68 / 0.89 (98) | 0.51 / 0.67 (134) | 0.25 / 0.31 (300) |
| period 1, 2000 samples, 3 reads, pcr | B_raw | 0.45 / 0.74 (42) | 0.57 / 0.83 (42) | 0.57 / 0.81 (42) | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.55 / 0.67 (201) | 0.38 / 0.49 (450) |
| period 1, 2000 samples, 3 reads, pcr | B_curves | — | — | — | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.55 / 0.67 (201) | 0.38 / 0.49 (450) |
| period 1, 2000 samples, 3 reads, pcr | D | 0.61 / 0.86 (28) | 0.61 / 0.86 (28) | 0.68 / 0.93 (28) | 0.71 / 1.00 (28) | 0.62 / 0.92 (98) | 0.55 / 0.68 (134) | 0.27 / 0.32 (300) |
| period 1, 2000 samples, 3 reads, pcr | D_tract_prior | 0.55 / 0.76 (42) | 0.50 / 0.81 (42) | 0.57 / 0.76 (42) | 0.71 / 0.98 (42) | 0.69 / 0.93 (147) | 0.55 / 0.67 (201) | 0.38 / 0.50 (450) |
| period 1, 2000 samples, 3 reads, pcr_free | B_raw | 0.57 / 0.79 (14) | 0.79 / 0.93 (14) | 0.57 / 0.71 (14) | 0.57 / 1.00 (14) | 0.61 / 0.90 (49) | 0.55 / 0.66 (67) | 0.31 / 0.47 (150) |
| period 1, 2000 samples, 3 reads, pcr_free | B_curves | — | — | — | 0.57 / 1.00 (14) | 0.61 / 0.90 (49) | 0.55 / 0.66 (67) | 0.31 / 0.47 (150) |
| period 1, 2000 samples, 3 reads, pcr_free | D | 0.64 / 0.93 (14) | 0.64 / 1.00 (14) | 0.57 / 0.86 (14) | 0.64 / 1.00 (14) | 0.57 / 0.90 (49) | 0.54 / 0.67 (67) | 0.26 / 0.31 (150) |
| period 1, 2000 samples, 3 reads, pcr_free | D_tract_prior | 0.57 / 0.86 (14) | 0.43 / 1.00 (14) | 0.50 / 0.79 (14) | 0.57 / 1.00 (14) | 0.63 / 0.90 (49) | 0.55 / 0.66 (67) | 0.33 / 0.46 (150) |
| period 1, 500 samples, 10 reads, pcr | B_raw | 0.50 / 0.79 (42) | 0.64 / 0.98 (42) | 0.60 / 0.98 (42) | 0.69 / 0.95 (42) | 0.65 / 0.90 (147) | 0.52 / 0.68 (201) | 0.34 / 0.50 (450) |
| period 1, 500 samples, 10 reads, pcr | B_curves | — | — | — | 0.69 / 0.95 (42) | 0.65 / 0.90 (147) | 0.52 / 0.68 (201) | 0.34 / 0.50 (450) |
| period 1, 500 samples, 10 reads, pcr | D | 0.50 / 0.81 (42) | 0.67 / 0.95 (42) | 0.57 / 0.95 (42) | 0.69 / 0.93 (42) | 0.67 / 0.90 (147) | 0.53 / 0.67 (201) | 0.28 / 0.35 (450) |
| period 1, 500 samples, 10 reads, pcr | D_tract_prior | 0.50 / 0.79 (42) | 0.67 / 0.95 (42) | 0.55 / 0.95 (42) | 0.67 / 0.90 (42) | 0.65 / 0.90 (147) | 0.53 / 0.68 (201) | 0.36 / 0.50 (450) |
| period 1, 500 samples, 3 reads, pcr | B_raw | 0.38 / 0.83 (42) | 0.69 / 0.93 (42) | 0.64 / 0.93 (42) | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.64 (201) | 0.37 / 0.47 (450) |
| period 1, 500 samples, 3 reads, pcr | B_curves | — | — | — | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.64 (201) | 0.37 / 0.47 (450) |
| period 1, 500 samples, 3 reads, pcr | D | 0.48 / 0.88 (42) | 0.38 / 0.62 (42) | 0.83 / 0.95 (42) | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.64 (201) | 0.28 / 0.32 (450) |
| period 1, 500 samples, 3 reads, pcr | D_tract_prior | 0.38 / 0.83 (42) | 0.29 / 0.62 (42) | 0.81 / 0.95 (42) | 0.71 / 0.98 (42) | 0.67 / 0.93 (147) | 0.54 / 0.64 (201) | 0.37 / 0.49 (450) |
| period 1, 500 samples, 3 reads, pcr_free | B_raw | 0.64 / 0.93 (42) | 0.62 / 0.95 (42) | 0.76 / 0.93 (42) | 0.69 / 0.98 (42) | 0.67 / 0.93 (147) | 0.52 / 0.64 (201) | 0.36 / 0.48 (450) |
| period 1, 500 samples, 3 reads, pcr_free | B_curves | — | — | — | 0.69 / 0.98 (42) | 0.67 / 0.93 (147) | 0.52 / 0.64 (201) | 0.36 / 0.48 (450) |
| period 1, 500 samples, 3 reads, pcr_free | D_tract_prior | 0.71 / 0.93 (42) | 0.57 / 0.90 (42) | 0.79 / 1.00 (42) | 0.69 / 0.98 (42) | 0.67 / 0.93 (147) | 0.53 / 0.64 (201) | 0.37 / 0.48 (450) |
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
| period 1, 60 samples, 10 reads, pcr_free | D_tract_prior | 0.64 / 0.95 (42) | 0.86 / 0.93 (42) | 0.67 / 0.93 (42) | 0.62 / 1.00 (42) | 0.67 / 0.94 (147) | 0.55 / 0.68 (201) | 0.37 / 0.54 (450) |
| period 1, 60 samples, 3 reads, pcr | B_raw | 0.79 / 0.95 (42) | 0.74 / 0.90 (42) | 0.67 / 0.95 (42) | 0.69 / 0.95 (42) | 0.61 / 0.86 (147) | 0.53 / 0.70 (201) | 0.36 / 0.51 (450) |
| period 1, 60 samples, 3 reads, pcr | B_curves | — | — | — | 0.69 / 0.95 (42) | 0.61 / 0.86 (147) | 0.53 / 0.70 (201) | 0.36 / 0.51 (450) |
| period 1, 60 samples, 3 reads, pcr | D | 0.74 / 0.98 (42) | 0.48 / 0.74 (42) | 0.62 / 0.93 (42) | 0.67 / 0.88 (42) | 0.65 / 0.87 (147) | 0.54 / 0.70 (201) | 0.26 / 0.35 (450) |
| period 1, 60 samples, 3 reads, pcr | D_tract_prior | 0.64 / 0.95 (42) | 0.48 / 0.69 (42) | 0.64 / 0.90 (42) | 0.67 / 0.88 (42) | 0.63 / 0.87 (147) | 0.55 / 0.70 (201) | 0.42 / 0.58 (450) |
| period 1, 60 samples, 30 reads, pcr | B_raw | 0.76 / 0.90 (42) | 0.74 / 0.98 (42) | 0.76 / 0.90 (42) | 0.74 / 0.93 (42) | 0.65 / 0.90 (147) | 0.47 / 0.66 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 30 reads, pcr | B_curves | — | — | — | 0.74 / 0.93 (42) | 0.65 / 0.90 (147) | 0.47 / 0.66 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 30 reads, pcr | D | 0.76 / 0.90 (42) | 0.67 / 0.93 (42) | 0.74 / 0.95 (42) | 0.62 / 0.93 (42) | 0.65 / 0.89 (147) | 0.46 / 0.66 (201) | 0.27 / 0.32 (450) |
| period 1, 60 samples, 30 reads, pcr | D_tract_prior | 0.76 / 0.90 (42) | 0.64 / 0.93 (42) | 0.74 / 0.95 (42) | 0.71 / 0.93 (42) | 0.66 / 0.90 (147) | 0.47 / 0.66 (201) | 0.36 / 0.47 (450) |
| period 1, 60 samples, 300 reads, pcr | B_raw | 0.64 / 0.93 (42) | 0.60 / 0.95 (42) | 0.48 / 0.88 (42) | 0.69 / 0.93 (42) | 0.61 / 0.91 (147) | 0.49 / 0.67 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 300 reads, pcr | B_curves | — | — | — | 0.69 / 0.93 (42) | 0.61 / 0.91 (147) | 0.49 / 0.67 (201) | 0.32 / 0.42 (450) |
| period 1, 60 samples, 300 reads, pcr | D | 0.64 / 0.93 (42) | 0.64 / 0.95 (42) | 0.50 / 0.88 (42) | 0.55 / 0.93 (42) | 0.61 / 0.91 (147) | 0.50 / 0.66 (201) | 0.31 / 0.35 (450) |
| period 1, 60 samples, 300 reads, pcr | D_tract_prior | 0.64 / 0.93 (42) | 0.64 / 0.95 (42) | 0.50 / 0.88 (42) | 0.60 / 0.93 (42) | 0.62 / 0.93 (147) | 0.50 / 0.66 (201) | 0.36 / 0.46 (450) |
| period 2, 1 samples, 10 reads, pcr | B_raw | 0.71 / 0.92 (126) | 0.65 / 0.93 (124) | 0.65 / 0.87 (124) | 0.91 / 1.00 (11) | 0.65 / 0.89 (354) | 0.36 / 0.42 (744) | 0.19 / 0.28 (1296) |
| period 2, 1 samples, 10 reads, pcr | B_curves | — | — | — | 0.91 / 1.00 (11) | 0.65 / 0.89 (354) | 0.36 / 0.42 (744) | 0.19 / 0.28 (1296) |
| period 2, 1 samples, 10 reads, pcr | D | 0.68 / 0.96 (126) | 0.47 / 0.90 (126) | 0.75 / 0.97 (126) | — | 0.55 / 0.87 (354) | 0.71 / 0.91 (744) | 0.60 / 0.91 (1296) |
| period 2, 1 samples, 10 reads, pcr | D_tract_prior | 0.68 / 0.96 (126) | 0.47 / 0.90 (126) | 0.75 / 0.97 (126) | — | 0.50 / 0.82 (354) | 0.71 / 0.94 (744) | 0.46 / 0.88 (1296) |
| period 2, 1 samples, 3 reads, pcr | B_raw | 0.66 / 0.92 (125) | 0.65 / 0.90 (109) | 0.65 / 0.88 (108) | — | 0.65 / 0.88 (354) | 0.36 / 0.42 (744) | 0.19 / 0.30 (1296) |
| period 2, 1 samples, 3 reads, pcr | B_curves | — | — | — | — | 0.65 / 0.88 (354) | 0.36 / 0.42 (744) | 0.19 / 0.30 (1296) |
| period 2, 1 samples, 3 reads, pcr | D | 0.68 / 0.90 (126) | 0.49 / 0.83 (126) | 0.60 / 0.86 (126) | — | 0.66 / 0.90 (354) | 0.71 / 0.93 (744) | 0.56 / 0.88 (1296) |
| period 2, 1 samples, 3 reads, pcr | D_tract_prior | 0.68 / 0.90 (126) | 0.49 / 0.83 (126) | 0.60 / 0.86 (126) | — | 0.59 / 0.85 (354) | 0.74 / 0.95 (744) | 0.42 / 0.87 (1296) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_raw | 0.66 / 0.90 (62) | 0.67 / 0.89 (55) | 0.60 / 0.84 (58) | 1.00 / 1.00 (1) | 0.60 / 0.89 (176) | 0.30 / 0.38 (371) | 0.19 / 0.30 (631) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | B_curves | — | — | — | 1.00 / 1.00 (1) | 0.60 / 0.89 (176) | 0.30 / 0.38 (371) | 0.19 / 0.30 (631) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D | 0.78 / 0.95 (63) | 0.67 / 1.00 (63) | 0.60 / 0.87 (63) | — | 0.59 / 0.87 (177) | 0.69 / 0.89 (372) | 0.59 / 0.93 (648) |
| period 2, 1 samples, 3 reads, pcr_F0.4 | D_tract_prior | 0.78 / 0.95 (63) | 0.67 / 1.00 (63) | 0.62 / 0.92 (63) | — | 0.55 / 0.84 (177) | 0.70 / 0.92 (372) | 0.46 / 0.89 (648) |
| period 2, 1 samples, 3 reads, pcr_free | B_raw | 0.90 / 1.00 (41) | 0.86 / 1.00 (7) | 0.17 / 0.17 (18) | 1.00 / 1.00 (1) | 0.67 / 0.90 (150) | 0.36 / 0.40 (320) | 0.21 / 0.30 (537) |
| period 2, 1 samples, 3 reads, pcr_free | B_curves | — | — | — | 1.00 / 1.00 (1) | 0.67 / 0.90 (150) | 0.36 / 0.40 (320) | 0.21 / 0.30 (537) |
| period 2, 1 samples, 3 reads, pcr_free | D | 0.44 / 0.78 (63) | 0.46 / 0.70 (63) | 0.44 / 0.76 (63) | — | 0.68 / 0.90 (177) | 0.67 / 0.90 (372) | 0.52 / 0.88 (648) |
| period 2, 1 samples, 3 reads, pcr_free | D_tract_prior | 0.44 / 0.78 (63) | 0.46 / 0.70 (63) | 0.44 / 0.76 (63) | — | 0.61 / 0.84 (177) | 0.70 / 0.92 (372) | 0.40 / 0.88 (648) |
| period 2, 1 samples, 30 reads, pcr | B_raw | 0.59 / 0.91 (126) | 0.62 / 0.94 (126) | 0.62 / 0.92 (125) | 0.96 / 1.00 (26) | 0.71 / 0.88 (354) | 0.35 / 0.43 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 30 reads, pcr | B_curves | — | — | — | 0.96 / 1.00 (26) | 0.71 / 0.88 (354) | 0.35 / 0.43 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 30 reads, pcr | D | 0.60 / 0.87 (126) | 0.58 / 0.93 (126) | 0.75 / 0.98 (126) | — | 0.60 / 0.87 (354) | 0.68 / 0.91 (744) | 0.62 / 0.93 (1296) |
| period 2, 1 samples, 30 reads, pcr | D_tract_prior | 0.60 / 0.87 (126) | 0.57 / 0.93 (126) | 0.75 / 0.98 (126) | — | 0.50 / 0.83 (354) | 0.70 / 0.95 (744) | 0.49 / 0.92 (1296) |
| period 2, 1 samples, 300 reads, pcr | B_raw | 0.71 / 0.97 (126) | 0.67 / 0.95 (126) | 0.69 / 0.94 (126) | 1.00 / 1.00 (29) | 0.67 / 0.91 (354) | 0.35 / 0.41 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 300 reads, pcr | B_curves | — | — | — | 1.00 / 1.00 (29) | 0.67 / 0.91 (354) | 0.35 / 0.41 (744) | 0.18 / 0.26 (1296) |
| period 2, 1 samples, 300 reads, pcr | D | 0.64 / 0.95 (126) | 0.65 / 0.96 (126) | 0.77 / 0.98 (126) | — | 0.60 / 0.92 (354) | 0.64 / 0.88 (744) | 0.58 / 0.88 (1296) |
| period 2, 1 samples, 300 reads, pcr | D_tract_prior | 0.64 / 0.95 (126) | 0.65 / 0.96 (126) | 0.77 / 0.98 (126) | — | 0.57 / 0.85 (354) | 0.69 / 0.91 (744) | 0.45 / 0.89 (1296) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_raw | 0.68 / 0.95 (63) | 0.54 / 0.87 (63) | 0.59 / 0.95 (63) | 0.91 / 1.00 (11) | 0.67 / 0.93 (177) | 0.35 / 0.41 (372) | 0.17 / 0.26 (648) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | B_curves | — | — | — | 0.91 / 1.00 (11) | 0.67 / 0.93 (177) | 0.35 / 0.41 (372) | 0.17 / 0.26 (648) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D | 0.65 / 0.95 (63) | 0.71 / 0.92 (63) | 0.87 / 1.00 (63) | — | 0.51 / 0.82 (177) | 0.71 / 0.91 (372) | 0.57 / 0.89 (648) |
| period 2, 1 samples, 300 reads, pcr_F0.4 | D_tract_prior | 0.65 / 0.95 (63) | 0.71 / 0.92 (63) | 0.87 / 1.00 (63) | — | 0.45 / 0.82 (177) | 0.74 / 0.93 (372) | 0.46 / 0.91 (648) |
| period 2, 1 samples, 300 reads, pcr_free | B_raw | 0.70 / 0.92 (63) | 0.60 / 0.89 (63) | 0.59 / 0.90 (63) | 1.00 / 1.00 (14) | 0.68 / 0.91 (177) | 0.34 / 0.40 (372) | 0.18 / 0.27 (648) |
| period 2, 1 samples, 300 reads, pcr_free | B_curves | — | — | — | 1.00 / 1.00 (14) | 0.68 / 0.91 (177) | 0.34 / 0.40 (372) | 0.18 / 0.27 (648) |
| period 2, 1 samples, 300 reads, pcr_free | D | 0.54 / 0.90 (63) | 0.70 / 0.95 (63) | 0.67 / 0.90 (63) | — | 0.53 / 0.90 (177) | 0.66 / 0.87 (372) | 0.60 / 0.90 (648) |
| period 2, 1 samples, 300 reads, pcr_free | D_tract_prior | 0.54 / 0.90 (63) | 0.70 / 0.94 (63) | 0.67 / 0.90 (63) | — | 0.54 / 0.81 (177) | 0.72 / 0.91 (372) | 0.46 / 0.90 (648) |
| period 2, 10 samples, 10 reads, pcr | B_raw | 0.70 / 0.97 (63) | 0.70 / 0.97 (63) | 0.68 / 0.94 (63) | 0.67 / 0.90 (63) | 0.66 / 0.89 (177) | 0.44 / 0.55 (372) | 0.26 / 0.38 (648) |
| period 2, 10 samples, 10 reads, pcr | B_curves | — | — | — | 0.67 / 0.90 (63) | 0.66 / 0.89 (177) | 0.44 / 0.55 (372) | 0.26 / 0.38 (648) |
| period 2, 10 samples, 10 reads, pcr | D | 0.70 / 0.98 (63) | 0.84 / 0.97 (63) | 0.68 / 0.97 (63) | 0.75 / 0.89 (63) | 0.68 / 0.89 (177) | 0.46 / 0.59 (372) | 0.35 / 0.47 (648) |
| period 2, 10 samples, 10 reads, pcr | D_tract_prior | 0.70 / 0.98 (63) | 0.84 / 0.97 (63) | 0.71 / 0.97 (63) | 0.68 / 0.89 (63) | 0.68 / 0.90 (177) | 0.55 / 0.73 (372) | 0.56 / 0.79 (648) |
| period 2, 10 samples, 3 reads, pcr | B_raw | 0.63 / 0.95 (63) | 0.65 / 0.90 (62) | 0.59 / 0.87 (61) | 0.78 / 0.95 (63) | 0.63 / 0.84 (177) | 0.45 / 0.54 (372) | 0.29 / 0.41 (648) |
| period 2, 10 samples, 3 reads, pcr | B_curves | — | — | — | 0.78 / 0.95 (63) | 0.63 / 0.84 (177) | 0.45 / 0.54 (372) | 0.29 / 0.41 (648) |
| period 2, 10 samples, 3 reads, pcr | D | 0.59 / 0.90 (63) | 0.49 / 0.76 (63) | 0.79 / 1.00 (63) | 0.78 / 0.94 (63) | 0.63 / 0.86 (177) | 0.50 / 0.63 (372) | 0.37 / 0.50 (648) |
| period 2, 10 samples, 3 reads, pcr | D_tract_prior | 0.59 / 0.89 (63) | 0.43 / 0.75 (63) | 0.83 / 1.00 (63) | 0.59 / 1.00 (63) | 0.62 / 0.85 (177) | 0.60 / 0.75 (372) | 0.61 / 0.83 (648) |
| period 2, 10 samples, 30 reads, pcr | B_raw | 0.59 / 0.92 (63) | 0.67 / 0.92 (63) | 0.60 / 0.92 (63) | 0.56 / 0.92 (63) | 0.66 / 0.88 (177) | 0.41 / 0.52 (372) | 0.25 / 0.35 (648) |
| period 2, 10 samples, 30 reads, pcr | B_curves | — | — | — | 0.56 / 0.92 (63) | 0.66 / 0.88 (177) | 0.41 / 0.52 (372) | 0.25 / 0.35 (648) |
| period 2, 10 samples, 30 reads, pcr | D | 0.62 / 0.89 (63) | 0.60 / 0.98 (63) | 0.62 / 0.97 (63) | 0.43 / 0.78 (63) | 0.65 / 0.88 (177) | 0.42 / 0.56 (372) | 0.31 / 0.43 (648) |
| period 2, 10 samples, 30 reads, pcr | D_tract_prior | 0.62 / 0.89 (63) | 0.62 / 0.98 (63) | 0.62 / 0.97 (63) | 0.41 / 0.73 (63) | 0.66 / 0.89 (177) | 0.54 / 0.68 (372) | 0.56 / 0.79 (648) |
| period 2, 10 samples, 300 reads, pcr | B_raw | 0.65 / 0.90 (63) | 0.67 / 0.90 (63) | 0.62 / 0.94 (63) | 0.71 / 0.94 (63) | 0.66 / 0.90 (177) | 0.45 / 0.54 (372) | 0.26 / 0.36 (648) |
| period 2, 10 samples, 300 reads, pcr | B_curves | — | — | — | 0.71 / 0.94 (63) | 0.66 / 0.90 (177) | 0.45 / 0.54 (372) | 0.26 / 0.36 (648) |
| period 2, 10 samples, 300 reads, pcr | D | 0.59 / 0.92 (63) | 0.68 / 0.98 (63) | 0.56 / 0.90 (63) | 0.83 / 0.98 (63) | 0.67 / 0.89 (177) | 0.46 / 0.56 (372) | 0.36 / 0.47 (648) |
| period 2, 10 samples, 300 reads, pcr | D_tract_prior | 0.59 / 0.92 (63) | 0.70 / 0.98 (63) | 0.56 / 0.90 (63) | 0.81 / 0.97 (63) | 0.67 / 0.90 (177) | 0.54 / 0.69 (372) | 0.55 / 0.80 (648) |
| period 2, 2000 samples, 10 reads, pcr | B_raw | 0.60 / 0.95 (42) | 0.60 / 0.93 (42) | 0.67 / 0.93 (42) | 0.57 / 0.88 (42) | 0.56 / 0.78 (118) | 0.44 / 0.57 (248) | 0.30 / 0.42 (432) |
| period 2, 2000 samples, 10 reads, pcr | B_curves | — | — | — | 0.57 / 0.88 (42) | 0.56 / 0.78 (118) | 0.44 / 0.57 (248) | 0.30 / 0.42 (432) |
| period 2, 2000 samples, 10 reads, pcr | D | 0.74 / 1.00 (42) | 0.60 / 0.95 (42) | 0.55 / 0.95 (42) | 0.60 / 0.86 (42) | 0.55 / 0.79 (118) | 0.43 / 0.56 (248) | 0.22 / 0.25 (432) |
| period 2, 2000 samples, 3 reads, pcr | B_raw | 0.51 / 0.83 (63) | 0.62 / 0.92 (63) | 0.52 / 0.87 (63) | 0.56 / 0.89 (63) | 0.63 / 0.89 (177) | 0.44 / 0.57 (372) | 0.26 / 0.39 (648) |
| period 2, 2000 samples, 3 reads, pcr | B_curves | — | — | — | 0.56 / 0.89 (63) | 0.63 / 0.89 (177) | 0.44 / 0.57 (372) | 0.26 / 0.39 (648) |
| period 2, 2000 samples, 3 reads, pcr | D | 0.73 / 0.94 (63) | 0.67 / 0.90 (63) | 0.67 / 0.97 (63) | 0.52 / 0.86 (63) | 0.66 / 0.90 (177) | 0.41 / 0.55 (372) | 0.19 / 0.22 (648) |
| period 2, 2000 samples, 3 reads, pcr | D_tract_prior | 0.46 / 0.84 (63) | 0.57 / 0.86 (63) | 0.40 / 0.79 (63) | 0.56 / 0.89 (63) | 0.65 / 0.90 (177) | 0.45 / 0.57 (372) | 0.27 / 0.39 (648) |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_raw | 0.52 / 0.86 (21) | 0.67 / 1.00 (21) | 0.48 / 0.76 (21) | 0.62 / 0.90 (21) | 0.59 / 0.95 (59) | 0.47 / 0.61 (124) | 0.36 / 0.54 (216) |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | B_curves | — | — | — | 0.62 / 0.90 (21) | 0.59 / 0.95 (59) | 0.47 / 0.61 (124) | 0.36 / 0.54 (216) |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | D | 0.81 / 0.90 (21) | 0.76 / 1.00 (21) | 0.57 / 0.95 (21) | 0.62 / 0.76 (21) | 0.63 / 0.93 (59) | 0.45 / 0.60 (124) | 0.27 / 0.33 (216) |
| period 2, 2000 samples, 3 reads, pcr_F0.4 | D_tract_prior | 0.52 / 0.86 (21) | 0.67 / 1.00 (21) | 0.43 / 0.81 (21) | 0.67 / 0.86 (21) | 0.63 / 0.95 (59) | 0.48 / 0.63 (124) | 0.38 / 0.55 (216) |
| period 2, 2000 samples, 3 reads, pcr_free | B_raw | 0.81 / 0.95 (21) | 0.81 / 0.95 (21) | 0.38 / 0.90 (21) | 0.62 / 0.95 (21) | 0.66 / 0.88 (59) | 0.40 / 0.53 (124) | 0.31 / 0.45 (216) |
| period 2, 2000 samples, 3 reads, pcr_free | B_curves | — | — | — | 0.62 / 0.95 (21) | 0.66 / 0.88 (59) | 0.40 / 0.53 (124) | 0.31 / 0.45 (216) |
| period 2, 2000 samples, 3 reads, pcr_free | D | 0.71 / 0.95 (21) | 0.43 / 1.00 (21) | 0.29 / 0.81 (21) | 0.57 / 0.95 (21) | 0.66 / 0.90 (59) | 0.37 / 0.50 (124) | 0.21 / 0.25 (216) |
| period 2, 2000 samples, 3 reads, pcr_free | D_tract_prior | 0.81 / 1.00 (21) | 0.57 / 1.00 (21) | 0.62 / 0.81 (21) | 0.62 / 0.95 (21) | 0.66 / 0.88 (59) | 0.42 / 0.53 (124) | 0.32 / 0.44 (216) |
| period 2, 500 samples, 10 reads, pcr | B_raw | 0.73 / 0.98 (63) | 0.60 / 0.95 (63) | 0.65 / 0.95 (63) | 0.75 / 0.90 (63) | 0.61 / 0.87 (177) | 0.41 / 0.52 (372) | 0.26 / 0.38 (648) |
| period 2, 500 samples, 10 reads, pcr | B_curves | — | — | — | 0.75 / 0.90 (63) | 0.61 / 0.87 (177) | 0.41 / 0.52 (372) | 0.26 / 0.38 (648) |
| period 2, 500 samples, 10 reads, pcr | D | 0.67 / 0.97 (63) | 0.68 / 0.97 (63) | 0.67 / 0.92 (63) | 0.71 / 0.94 (63) | 0.60 / 0.88 (177) | 0.41 / 0.51 (372) | 0.18 / 0.21 (648) |
| period 2, 500 samples, 10 reads, pcr | D_tract_prior | 0.65 / 0.97 (63) | 0.67 / 0.95 (63) | 0.65 / 0.92 (63) | 0.75 / 0.92 (63) | 0.64 / 0.88 (177) | 0.43 / 0.52 (372) | 0.27 / 0.38 (648) |
| period 2, 500 samples, 3 reads, pcr | B_raw | 0.56 / 0.89 (63) | 0.73 / 0.98 (63) | 0.63 / 1.00 (63) | 0.71 / 0.92 (63) | 0.64 / 0.93 (177) | 0.41 / 0.51 (372) | 0.27 / 0.40 (648) |
| period 2, 500 samples, 3 reads, pcr | B_curves | — | — | — | 0.71 / 0.92 (63) | 0.64 / 0.93 (177) | 0.41 / 0.51 (372) | 0.27 / 0.40 (648) |
| period 2, 500 samples, 3 reads, pcr | D | 0.67 / 0.97 (63) | 0.60 / 0.98 (63) | 0.79 / 0.97 (63) | 0.71 / 0.92 (63) | 0.62 / 0.93 (177) | 0.40 / 0.49 (372) | 0.18 / 0.22 (648) |
| period 2, 500 samples, 3 reads, pcr | D_tract_prior | 0.57 / 0.94 (63) | 0.76 / 0.98 (63) | 0.65 / 0.92 (63) | 0.71 / 0.92 (63) | 0.68 / 0.94 (177) | 0.42 / 0.51 (372) | 0.28 / 0.40 (648) |
| period 2, 500 samples, 3 reads, pcr_F0.4 | B_raw | 0.51 / 0.89 (63) | 0.75 / 0.98 (63) | 0.65 / 0.94 (63) | 0.59 / 0.84 (63) | 0.64 / 0.92 (177) | 0.53 / 0.65 (372) | 0.31 / 0.51 (648) |
| period 2, 500 samples, 3 reads, pcr_F0.4 | B_curves | — | — | — | 0.59 / 0.84 (63) | 0.64 / 0.92 (177) | 0.53 / 0.65 (372) | 0.31 / 0.51 (648) |
| period 2, 500 samples, 3 reads, pcr_F0.4 | D_tract_prior | 0.56 / 0.94 (63) | 0.79 / 0.98 (63) | 0.60 / 0.83 (63) | 0.59 / 0.84 (63) | 0.67 / 0.93 (177) | 0.54 / 0.66 (372) | 0.33 / 0.53 (648) |
| period 2, 500 samples, 3 reads, pcr_free | B_raw | 0.65 / 0.97 (63) | 0.70 / 0.95 (63) | 0.71 / 0.98 (63) | 0.71 / 0.94 (63) | 0.64 / 0.93 (177) | 0.40 / 0.51 (372) | 0.27 / 0.41 (648) |
| period 2, 500 samples, 3 reads, pcr_free | B_curves | — | — | — | 0.71 / 0.94 (63) | 0.64 / 0.93 (177) | 0.40 / 0.51 (372) | 0.27 / 0.41 (648) |
| period 2, 500 samples, 3 reads, pcr_free | D_tract_prior | 0.60 / 0.89 (63) | 0.60 / 0.92 (63) | 0.60 / 0.98 (63) | 0.71 / 0.92 (63) | 0.67 / 0.93 (177) | 0.41 / 0.51 (372) | 0.29 / 0.41 (648) |
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
| period 2, 60 samples, 10 reads, pcr_F0.4 | D_tract_prior | 0.57 / 0.87 (63) | 0.84 / 0.95 (63) | 0.65 / 0.90 (63) | 0.70 / 0.94 (63) | 0.67 / 0.93 (177) | 0.51 / 0.66 (372) | 0.44 / 0.64 (648) |
| period 2, 60 samples, 10 reads, pcr_free | B_raw | 0.76 / 0.92 (63) | 0.68 / 0.94 (63) | 0.59 / 0.95 (63) | 0.67 / 0.97 (63) | 0.63 / 0.89 (177) | 0.43 / 0.53 (372) | 0.28 / 0.42 (648) |
| period 2, 60 samples, 10 reads, pcr_free | B_curves | — | — | — | 0.67 / 0.97 (63) | 0.63 / 0.89 (177) | 0.43 / 0.53 (372) | 0.28 / 0.42 (648) |
| period 2, 60 samples, 10 reads, pcr_free | D | 0.67 / 0.92 (63) | 0.57 / 0.84 (63) | 0.60 / 0.98 (63) | 0.62 / 0.97 (63) | 0.63 / 0.91 (177) | 0.41 / 0.52 (372) | 0.20 / 0.25 (648) |
| period 2, 60 samples, 10 reads, pcr_free | D_tract_prior | 0.68 / 0.92 (63) | 0.57 / 0.86 (63) | 0.62 / 0.98 (63) | 0.63 / 0.97 (63) | 0.66 / 0.90 (177) | 0.44 / 0.53 (372) | 0.35 / 0.50 (648) |
| period 2, 60 samples, 3 reads, pcr | B_raw | 0.67 / 0.87 (63) | 0.75 / 0.90 (63) | 0.65 / 0.95 (63) | 0.67 / 0.92 (63) | 0.68 / 0.88 (177) | 0.46 / 0.54 (372) | 0.26 / 0.40 (648) |
| period 2, 60 samples, 3 reads, pcr | B_curves | — | — | — | 0.67 / 0.92 (63) | 0.68 / 0.88 (177) | 0.46 / 0.54 (372) | 0.26 / 0.40 (648) |
| period 2, 60 samples, 3 reads, pcr | D | 0.54 / 0.83 (63) | 0.68 / 0.98 (63) | 0.78 / 0.98 (63) | 0.65 / 0.90 (63) | 0.67 / 0.90 (177) | 0.45 / 0.53 (372) | 0.18 / 0.24 (648) |
| period 2, 60 samples, 3 reads, pcr | D_tract_prior | 0.56 / 0.79 (63) | 0.65 / 0.98 (63) | 0.83 / 0.98 (63) | 0.68 / 0.90 (63) | 0.69 / 0.88 (177) | 0.47 / 0.55 (372) | 0.35 / 0.54 (648) |
| period 2, 60 samples, 30 reads, pcr | B_raw | 0.62 / 0.90 (63) | 0.79 / 0.92 (63) | 0.63 / 0.94 (63) | 0.68 / 0.94 (63) | 0.59 / 0.91 (177) | 0.43 / 0.55 (372) | 0.26 / 0.37 (648) |
| period 2, 60 samples, 30 reads, pcr | B_curves | — | — | — | 0.68 / 0.94 (63) | 0.59 / 0.91 (177) | 0.43 / 0.55 (372) | 0.26 / 0.37 (648) |
| period 2, 60 samples, 30 reads, pcr | D | 0.57 / 0.92 (63) | 0.76 / 0.98 (63) | 0.59 / 0.89 (63) | 0.65 / 0.94 (63) | 0.65 / 0.90 (177) | 0.42 / 0.54 (372) | 0.19 / 0.23 (648) |
| period 2, 60 samples, 30 reads, pcr | D_tract_prior | 0.57 / 0.92 (63) | 0.75 / 0.98 (63) | 0.59 / 0.89 (63) | 0.65 / 0.92 (63) | 0.63 / 0.92 (177) | 0.44 / 0.56 (372) | 0.31 / 0.46 (648) |
| period 2, 60 samples, 300 reads, pcr | B_raw | 0.67 / 0.95 (63) | 0.57 / 0.98 (63) | 0.57 / 0.90 (63) | 0.76 / 0.94 (63) | 0.58 / 0.88 (177) | 0.47 / 0.56 (372) | 0.26 / 0.36 (648) |
| period 2, 60 samples, 300 reads, pcr | B_curves | — | — | — | 0.76 / 0.94 (63) | 0.58 / 0.88 (177) | 0.47 / 0.56 (372) | 0.26 / 0.36 (648) |
| period 2, 60 samples, 300 reads, pcr | D | 0.68 / 0.94 (63) | 0.63 / 1.00 (63) | 0.57 / 0.97 (63) | 0.75 / 0.94 (63) | 0.57 / 0.88 (177) | 0.48 / 0.54 (372) | 0.20 / 0.25 (648) |
| period 2, 60 samples, 300 reads, pcr | D_tract_prior | 0.68 / 0.94 (63) | 0.63 / 1.00 (63) | 0.57 / 0.97 (63) | 0.76 / 0.94 (63) | 0.59 / 0.88 (177) | 0.48 / 0.56 (372) | 0.34 / 0.47 (648) |

