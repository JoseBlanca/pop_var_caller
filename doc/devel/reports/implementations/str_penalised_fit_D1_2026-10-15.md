# Step D1: the production fit on the bench, thin-class errors, and one-sample direction numbers

Measurement report, 2026-10-15 and 16, branch `str-penalised-fit`. Every run but the last 14 used
the bench built at 53a3e4f9 (after C10); the last 14 (§6) used it built at 8c6ce6f5, whose fitted
numbers are bit-identical to 53a3e4f9's (checked on the full reference fixture and end to end), so
the two sets pair. Plan step D1 (`doc/devel/ng/impl_plan/str_penalised_fit.md`, milestone D); spec
`doc/devel/ng/spec/str_penalised_fit.md` §1.2, §2.6, §2.10, §4.2. Every figure is on cohorts drawn
by the fitting bench (`examples/ng_str_fit_bench.rs`) from a known truth; **no time in it is a cost
measure**: runs up to 2026-10-15 shared the machine at four threads, the last 14 ran alone at 14.

**The fits compared.**

- **The production fit** (`D_tract_prior` on the bench, `PenalisedFitConfig::production()`): every
  repeat class of a motif length fitted at once, each number a penalised curve along repeat count,
  with the counting fit's prior on the allele-length shares and, since step C8, a prior of one read's
  worth a class on the three stutter numbers.
- **The counting fit** (`B` on the bench; `B_curves` after calling's stutter curves): each class
  fitted on its own, curves drawn afterwards. It is the baseline of
  `doc/devel/ng/research/str_fit/compare_D_prior_2026-10-14.md` (below, **the 2026-10-14 report**).

**Vocabulary.**

- **Held-out score.** The bench splits each class's loci into five groups; it fits on four and scores
  the fifth: the log-likelihood of the held-out loci's reads under the fitted numbers, minus the same
  under the true numbers. Zero is as good as the truth. The **paired difference** is production fit
  minus counting fit, seed by seed, in log-likelihood units a held-out locus; above zero, the
  production fit predicts better.
- **Coverage.** The share of estimates whose truth lies within one (or two) of the fit's own standard
  errors. Errors that are honest give about 68 in 100 within one and 95 in 100 within two.
- **Raw error.** The error the fit computes, before the rules that withhold it in the parameters file
  (`errors_of_a_class`, `penalised_fit.rs:982`): fewer than `LOCI_FOR_AN_ERROR` (8) loci with reads,
  a share carried by fewer than `LOCI_FOR_A_SHARE_ERROR` (3) loci, fewer slipped reads than the
  prior's weight, or an error above 3 on the logit scale.
- **Loci with reads.** The loci of a class at which at least one sample has a read spanning the
  repeat.

---

## 1. The answer

**Recommendation for checkpoint D: split `LOCI_FOR_AN_ERROR`, and lower the share rule.** All
figures are raw-error coverage on simulated classes (tomato-like truth, 1 to 2,000 samples, 3 and 30
reads, 3 seeds; §4), pooled over cohort sizes unless said otherwise.

- **Give the three stutter numbers (slip rate, shorter share, fall-off) an error from 3 loci with
  reads instead of 8.** At 3 to 7 loci the shorter share and the fall-off hold the truth within two
  errors 90 to 92 times in 100, as classes of 10 or more loci do (90 and 91). The slip rate holds it
  84 to 86 times in 100 against 92: it misses by more than two errors about 1 time in 7, where a thick
  class misses about 1 time in 12. At 1 or 2 loci the slip rate holds it 75 times in 100, and 59 at
  2,000 samples, so classes that thin keep no error.
- **Keep the line at 8 for F (the between-locus share) and for the allele-length shares.** F at 3 to
  7 loci holds the truth within two errors 72 to 80 times in 100 against 91 for thick classes, and 35
  at 1 or 2 loci. A thin class's shares are nearly all carried by under one locus, where they hold
  it 66 to 71 times in 100.
- **Move `LOCI_FOR_A_SHARE_ERROR` from 3 loci to 1.** In classes of 10 or more loci, a share carried
  by 1 to 3 loci (share × loci with reads) holds the truth within two errors 96 times in 100 (4,934
  shares), as often as shares carried by more. Shares carried by under one locus hold it 70 times in
  100 and stay withheld.
- **The recommendation does not depend on the error-matrix fix awaiting the owner** (below). In the
  runs that fix gives errors back to, the stutter numbers at 3 to 7 loci hold the truth within two
  errors 78 to 88 times in 100 (diagnosis report, its first table under option (c)).

**The results behind it, and the rest of D1:**

- **The production fit predicts held-out loci better than the counting fit in all 48 cells run**
  (1 to 2,000 samples, 3 to 300 reads, three truths), by +0.009 to +0.162 log-likelihood units a
  held-out locus; in 167 of the 168 seed-runs. The exception is one-base 500 × 30 seed 2, at −0.002.
  Below 2,000 samples every margin is the 2026-10-14 figure to within 0.002 (§2).
- **The cells spec §1.2 named as not run** (§2.2, §3): with three seeds at every 2,000 × 3 cell the
  margin is +0.009 ± 0.003 to +0.048 ± 0.002; at 500 × 30 it is +0.009 ± 0.006 (one-base) and
  +0.022 ± 0.008 (two-base); at two-base 2,000 × 10, +0.067 ± 0.018.
- **Two-base 2,000 × 3 seed 3**, where the production fit scores 0.015 a locus worse than the
  prototype did on the same draw, is within the seed-to-seed spread of the margin (+0.022, +0.053,
  +0.064 over seeds 3, 2, 1), and is the only such move: over the other 152 seed-runs that have a
  prototype score, none moved by more than 0.0032. Not investigated further (§2.2).
- **About one thin-design run in five gave no stutter error for any class; the cause is now known**
  (§4.1): the error matrix leaves out the curvature of the allele-length prior. A fix is awaiting the
  owner's decision (`str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`). It changes errors only;
  no held-out score in this report moves.
- **One-sample two-base direction numbers** (§5): C8's review counts are reproduced exactly (shorter
  share 31/63 and 21/126 with the prior against 18/62 and 12/126 without; fall-off 28/126 against
  18/126).
- **No code was changed.**

---

## 2. Held-out prediction at the 2026-10-14 report's cells

**Below 2,000 samples, the production fit's margin over the counting fit is what the 2026-10-14 report
measured, to within 0.002 log units a held-out locus in all 40 cells.** The stutter prior (C8) and the
review fixes since the prototype moved no paired difference by more than 0.002; the largest move is
+0.002 (one-base, PCR-free, 1 × 3). The production fit still predicts better than the counting fit in
every one of these cells.

Seed-paired difference, production fit minus counting fit after calling's curves, log-likelihood a
held-out locus, mean over seeds ± its standard error. At one sample both keep their own shares (the
`_own_shares` variants), as the 2026-10-14 report's §2 does; elsewhere both are scored as calling
would. "Default" is `CohortModel::tomato_like` with PCR stutter; "PCR-free" has a tenth of its slip
rate; "F = 0.4" sets every two-base class's between-locus share to 0.4.

| truth | motif length | samples | reads | seeds | now | 2026-10-14 | now − then |
|---|---|---|---|---|---|---|---|
| default | 1 | 1 | 3 | 6 | +0.102 ± 0.010 | +0.102 ± 0.010 | +0.000 |
| default | 2 | 1 | 3 | 6 | +0.133 ± 0.008 | +0.133 ± 0.008 | +0.001 |
| default | 1 | 1 | 10 | 6 | +0.144 ± 0.014 | +0.144 ± 0.014 | +0.000 |
| default | 2 | 1 | 10 | 6 | +0.151 ± 0.013 | +0.151 ± 0.013 | -0.000 |
| default | 1 | 1 | 30 | 6 | +0.106 ± 0.011 | +0.106 ± 0.011 | +0.000 |
| default | 2 | 1 | 30 | 6 | +0.143 ± 0.009 | +0.143 ± 0.009 | -0.000 |
| default | 1 | 1 | 300 | 6 | +0.115 ± 0.007 | +0.115 ± 0.007 | +0.000 |
| default | 2 | 1 | 300 | 6 | +0.141 ± 0.011 | +0.141 ± 0.011 | -0.000 |
| default | 1 | 10 | 3 | 3 | +0.100 ± 0.014 | +0.100 ± 0.014 | +0.000 |
| default | 2 | 10 | 3 | 3 | +0.148 ± 0.012 | +0.149 ± 0.012 | -0.001 |
| default | 1 | 10 | 10 | 3 | +0.088 ± 0.008 | +0.088 ± 0.008 | -0.000 |
| default | 2 | 10 | 10 | 3 | +0.147 ± 0.016 | +0.147 ± 0.016 | +0.000 |
| default | 1 | 10 | 30 | 3 | +0.095 ± 0.009 | +0.095 ± 0.009 | +0.000 |
| default | 2 | 10 | 30 | 3 | +0.105 ± 0.012 | +0.105 ± 0.012 | -0.000 |
| default | 1 | 10 | 300 | 3 | +0.082 ± 0.006 | +0.082 ± 0.006 | -0.000 |
| default | 2 | 10 | 300 | 3 | +0.113 ± 0.018 | +0.113 ± 0.018 | -0.000 |
| default | 1 | 60 | 3 | 3 | +0.034 ± 0.003 | +0.034 ± 0.002 | +0.000 |
| default | 2 | 60 | 3 | 3 | +0.039 ± 0.006 | +0.039 ± 0.005 | -0.000 |
| default | 1 | 60 | 10 | 3 | +0.027 ± 0.004 | +0.027 ± 0.004 | +0.000 |
| default | 2 | 60 | 10 | 3 | +0.040 ± 0.006 | +0.040 ± 0.006 | +0.000 |
| default | 1 | 60 | 30 | 3 | +0.026 ± 0.004 | +0.026 ± 0.004 | +0.000 |
| default | 2 | 60 | 30 | 3 | +0.081 ± 0.031 | +0.081 ± 0.031 | +0.000 |
| default | 1 | 60 | 300 | 3 | +0.031 ± 0.003 | +0.031 ± 0.003 | +0.000 |
| default | 2 | 60 | 300 | 3 | +0.054 ± 0.010 | +0.054 ± 0.010 | +0.000 |
| default | 1 | 500 | 3 | 3 | +0.010 ± 0.005 | +0.010 ± 0.006 | -0.001 |
| default | 2 | 500 | 3 | 3 | +0.034 ± 0.008 | +0.034 ± 0.008 | +0.000 |
| default | 1 | 500 | 10 | 3 | +0.026 ± 0.006 | +0.026 ± 0.006 | +0.000 |
| default | 2 | 500 | 10 | 3 | +0.042 ± 0.010 | +0.042 ± 0.010 | +0.000 |
| F = 0.4 | 2 | 1 | 3 | 3 | +0.118 ± 0.008 | +0.118 ± 0.009 | +0.000 |
| F = 0.4 | 2 | 1 | 300 | 3 | +0.128 ± 0.002 | +0.128 ± 0.002 | +0.000 |
| F = 0.4 | 2 | 60 | 10 | 3 | +0.055 ± 0.004 | +0.055 ± 0.003 | +0.000 |
| F = 0.4 | 2 | 500 | 3 | 3 | +0.043 ± 0.012 | +0.044 ± 0.010 | -0.001 |
| PCR-free | 1 | 1 | 3 | 3 | +0.109 ± 0.005 | +0.107 ± 0.005 | +0.002 |
| PCR-free | 2 | 1 | 3 | 3 | +0.162 ± 0.003 | +0.162 ± 0.004 | +0.000 |
| PCR-free | 1 | 1 | 300 | 3 | +0.110 ± 0.002 | +0.110 ± 0.002 | +0.000 |
| PCR-free | 2 | 1 | 300 | 3 | +0.131 ± 0.019 | +0.131 ± 0.019 | -0.000 |
| PCR-free | 1 | 60 | 10 | 3 | +0.029 ± 0.008 | +0.029 ± 0.008 | +0.000 |
| PCR-free | 2 | 60 | 10 | 3 | +0.040 ± 0.002 | +0.040 ± 0.002 | -0.000 |
| PCR-free | 1 | 500 | 3 | 3 | +0.012 ± 0.003 | +0.012 ± 0.003 | -0.000 |
| PCR-free | 2 | 500 | 3 | 3 | +0.028 ± 0.006 | +0.029 ± 0.006 | -0.000 |

Script `tmp/d1/heldout.py tmp/d1/grid <the 2026-10-14 run directories>`; full table with each fit's
loss against the truth in `tmp/d1/heldout_a.md`. For scale, the counting fit loses 0.14 to 0.20 a
locus against the truth at one sample, 0.19 to 0.30 at 10, 0.23 to 0.38 at 60 and 0.43 to 0.69 at 500.

### 2.1 The counting fit was not rerun where its 2026-10-14 scores exist

The counting fit's code has not changed since 2026-10-14, and the bench's draw and scorer give the
same numbers. Checked by running it again at four cells (two-base, seed 1: 1 × 3, 60 × 10, 500 × 3 and
2,000 × 3) and comparing every one of its 1,448 rows a cell with the 2026-10-14 rows
(`tmp/d1/bdiff.py`):

- at 1 × 3 every row is identical to the 8 printed decimals;
- elsewhere its numbers differ by at most 2×10⁻⁷ and its held-out score over all loci by at most
  2.5×10⁻⁵ log units a locus (the 2,000-sample cell, from the Monte Carlo sum's order), 400 times
  smaller than the smallest paired difference in the table;
- the truth's own held-out score, which every difference is taken against, agrees to 2×10⁻⁷.

So in every cell the 2026-10-14 runs scored the counting fit, this run fitted only the production
fit (`--methods D_tract_prior`) and paired it with the counting fit's 2026-10-14 scores, seed by seed.
That is every cell in the table above, five of the six default 2,000 × 3 runs (two-base seed 1 was
run with both fits, as the timing probe), and seed 1 of 2,000 × 3
PCR-free (both motif lengths) and F = 0.4. The counting fit was run fresh, beside the production fit,
in the cells it was never scored in: 2,000 × 3 PCR-free and F = 0.4 seeds 2 and 3, 500 × 30, and
2,000 × 10. The 2026-10-14 rows are in `../pop_var_caller-str-fit-compare/tmp/grid`, `grid_b` and
`grid_dtp`; the list of reused runs is `tmp/d1/old_b.txt`.

### 2.2 2,000 samples × 3 reads

**Three seeds at every 2,000 × 3 cell; the production fit is ahead at every seed, by +0.004 to
+0.064.** Same measure as the table above. "Then" is the prototype on 2026-10-14, which ran one seed
in three of these cells.

| truth | motif length | seeds now | now | 2026-10-14 | seeds then | counting fit's seeds reused | margin by seed (1, 2, 3) |
|---|---|---|---|---|---|---|---|
| default | 1 | 3 | +0.009 ± 0.003 | +0.009 ± 0.003 | 3 | 3 | +0.014, +0.010, +0.004 |
| default | 2 | 3 | +0.046 ± 0.013 | +0.051 ± 0.008 | 3 | 2 (seed 1 rerun) | +0.064, +0.053, +0.022 |
| PCR-free | 1 | 3 | +0.014 ± 0.003 | +0.009 | 1 | 1 (seeds 2–3 fitted fresh) | +0.008, +0.019, +0.015 |
| PCR-free | 2 | 3 | +0.027 ± 0.009 | +0.024 | 1 | 1 (seeds 2–3 fitted fresh) | +0.024, +0.013, +0.043 |
| F = 0.4 | 2 | 3 | +0.048 ± 0.002 | +0.040 | 1 | 1 (seeds 2–3 fitted fresh) | +0.043, +0.050, +0.050 |

Held-out loci: 20 a class (276 a seed one-base, 388 two-base). For scale, the counting fit loses
0.53 (one-base) to 0.84 (two-base, F = 0.4) a held-out locus against the truth here.

- **Where a 2026-10-14 seed exists, the production fit scores what the prototype scored on the same
  draw**, except at one seed. Over the 153 seed-runs at every cell size that have both scores, the
  production fit's score against the truth moved by at most 0.0032 a locus in 152.
- **That one is two-base default, seed 3**: −0.628 a locus against the prototype's −0.614, a move of
  0.015. It is within the seed-to-seed spread of the margin over the counting fit in this cell
  (+0.022 at seed 3 against +0.053 and +0.064 at seeds 1 and 2; the PCR-free cell spreads from +0.013
  to +0.043), and the production fit is still ahead at that seed (counting fit −0.650). Which change
  since the prototype moved it (the stutter prior, the review fixes, or C10's fixed sample order) was
  not traced, by instruction. The six new seeds here have no prototype score to compare with; their
  margins, +0.013 to +0.050, lie inside the range the seeds with one span.


## 3. The cells spec §1.2 names as not run

**The production fit is ahead in all three cells, by the most at two-base 2,000 × 10.** The reduced
set the lead approved: three seeds at every 2,000 × 3 cell (§2.2); 500 × 30, both motif lengths;
2,000 × 10, two-base only (one-base was dropped for time). Both fits were fitted in the same run, since
the 2026-10-14 runs never scored these cells. Default truth, 20 held-out loci a class.

| motif length | samples | reads | seeds | production − counting fit | counting fit − truth | production fit − truth | margin by seed (1, 2, 3) |
|---|---|---|---|---|---|---|---|
| 1 | 500 | 30 | 3 | +0.009 ± 0.006 | −0.527 ± 0.057 | −0.517 ± 0.061 | +0.016, −0.002, +0.015 |
| 2 | 500 | 30 | 3 | +0.022 ± 0.008 | −0.668 ± 0.028 | −0.646 ± 0.028 | +0.034, +0.006, +0.025 |
| 2 | 2,000 | 10 | 3 | +0.067 ± 0.018 | −0.833 ± 0.066 | −0.766 ± 0.048 | +0.034, +0.071, +0.097 |

- **One-base 500 × 30 seed 2 is the only seed-run of the 168 in D1 where the counting fit predicts
  better**, by 0.002 a held-out locus. One-base margins from 500 samples up are small throughout
  (+0.009 to +0.014 a cell; +0.001 at one-base 500 × 3 seed 2), so this sits at the low end of that
  spread rather than apart from it.
- **Depth does not move the two-base margin beyond seed noise.** At 2,000 samples it is +0.046 ± 0.013
  at 3 reads and +0.067 ± 0.018 at 10; at 500 samples +0.034 ± 0.008 at 3, +0.042 ± 0.010 at 10 and
  +0.022 ± 0.008 at 30. Each pair differs by about one standard error of the difference or less.
- The production fit's raw errors in these cells hold the truth as the thick classes of §4.4 do: within
  one / within two errors, slip rate 0.66 / 0.95 (500 × 30, 105 class estimates) and 0.62 / 0.97
  (2,000 × 10, 63); shorter share 0.79 / 0.96 and 0.60 / 0.97; fall-off 0.66 / 0.96 and 0.57 / 0.95;
  F 0.58 / 0.90 and 0.60 / 0.89 (`tmp/d1/final/new_cells_cov.py`).

---

## 4. Errors on classes with 1 to 7 loci

**The question** (spec §2.6). The parameters file withholds every error of a class with fewer than 8
loci with reads, because no measurement covered such classes. This section measures how often the
raw errors of thin classes cover the truth, for the owner's decision on `LOCI_FOR_AN_ERROR`, and the
same for the shares' rule `LOCI_FOR_A_SHARE_ERROR`.

### 4.1 How it was measured

Two designs of loci a class, drawn from the 2026-10-14 report's default PCR truth
(`CohortModel::tomato_like`), both motif lengths, 1, 10, 60, 500 and 2,000 samples, 3 and 30 reads, 3
seeds each:

- **Every class thin** (`--loci N`, N = 1 to 7): all 14 one-base or 21 two-base classes hold the same
  N loci. No class has neighbours with more reads to borrow from along the curve. 420 runs.
- **A thin tail** (`tmp/d1/mixed_p1.toml`, `mixed_p2.toml`): the short classes as tomato2 has them,
  capped at 100 loci (one-base 8–14 repeats, two-base 6–14), and the long classes thin, 1 to 7 loci
  each (one-base 15–21, two-base 15–27). This is the shape the kimura cohort has (107 of 141 classes
  under 8 loci, `kimura_shapes_span9_2026-10-14.md`). The files reproduce `tomato_like`'s truth
  exactly on every class (checked on the drawn truth tables). 60 runs.

Each run fits the production fit once on every locus (`--fixture-out`, which records each class's
loci with reads, numbers and raw errors at full precision) and writes the truth (`--truth-out`).
`tmp/d1/coverage.py` joins them into `tmp/d1/cov_rows.tsv`; `tmp/d1/cov_summary.py` makes the
tables (`tmp/d1/cov_summary.md`). Command for one run:

```sh
RAYON_NUM_THREADS=4 target/release/examples/ng_str_fit_bench --period 2 --samples 60 --depth 3 \
  --seed 1 --truth pcr --loci 4 --methods D_tract_prior \
  --truth-out tmp/d1/cov/uniform4_p2_n60_d3_s1.truth.tsv --fixture-out tmp/d1/cov/uniform4_p2_n60_d3_s1.json
```

**About one thin-design run in five returned no stutter errors at all, and the cause is now
known.** In 89 of the 480 runs every class's slip rate, shorter share, fall-off and allele-length
shares came back without an error; F kept its own (`UndeterminedErrors::EveryNumberButBetweenLocusShare`).
86 of the 89 are at 3 reads; by cohort size, 5 at 1 sample, 3 at 10, 23 at 60, 25 at 500 and 33 at
2,000. In the two-base thin-tail design at 60 × 3 and 2,000 × 3 all three seeds lost every error,
including those of the six classes with 100 loci. On the 2026-10-14 grid (no class under 11 loci) it
happened in 1 run of 153 (`p2_n60_d10_pcr_free` seed 2).

The diagnosis is `str_penalised_fit_thin_errors_diagnosis_2026-10-15.md`. In short: every error of a
motif length is read from one curvature matrix, and that matrix leaves out the curvature of the prior
the fit puts on each class's allele-length shares. At a thin class the reads alone then curve the
wrong way along the trade between the two outermost allele lengths on one side of the window (8 and 9
repeats shorter, mostly), the matrix cannot be inverted, and the fit withholds every error of the motif
length. Adding the prior's curvature made it invertible in all 384 runs of up to 500 samples (the 2,000-
sample runs were not re-run). **That fix is awaiting the owner's decision.** It changes errors only:
the fitted values, and so every held-out score in §2 and §3, stay as they are.

**The coverage tables below use the current errors**, without that fix, over the estimates that did
get an error (the counts in brackets). If the fix is adopted, the stutter numbers' coverage in the
runs that lost every error is the diagnosis report's first table under option (c): at 3 to 7 loci,
within two errors, slip rate 0.86 and 0.88 and fall-off 0.88 and 0.88, as in the runs that never
failed (0.86 to 0.89); shorter share 0.78 and 0.84, below them (0.91 and 0.89).

### 4.2 Stutter numbers and F

Within one / within two raw errors (estimates with an error), pooled over both motif lengths, both
depths and 3 seeds. "8+" in this table is only the thick classes of the thin-tail design; §4.4 adds
the 2026-10-14 grid.

**Every class thin:**

| number | samples | 1–2 loci | 3–4 loci | 5–7 loci |
|---|---|---|---|---|
| slip rate | 1 | 0.43 / 0.75 (368) | 0.50 / 0.79 (420) | 0.52 / 0.83 (600) |
| | 10 | 0.53 / 0.87 (392) | 0.63 / 0.91 (406) | 0.60 / 0.90 (630) |
| | 60 | 0.54 / 0.82 (308) | 0.60 / 0.90 (280) | 0.59 / 0.90 (532) |
| | 500 | 0.38 / 0.69 (315) | 0.44 / 0.79 (294) | 0.52 / 0.84 (518) |
| | 2,000 | 0.29 / 0.59 (273) | 0.48 / 0.79 (238) | 0.54 / 0.84 (441) |
| shorter share | 1 | 0.40 / 0.76 (368) | 0.52 / 0.86 (420) | 0.60 / 0.92 (600) |
| | 10 | 0.51 / 0.78 (392) | 0.51 / 0.92 (406) | 0.40 / 0.88 (630) |
| | 60 | 0.58 / 0.83 (308) | 0.57 / 0.91 (280) | 0.56 / 0.82 (532) |
| | 500 | 0.57 / 0.93 (315) | 0.71 / 0.96 (294) | 0.63 / 0.92 (518) |
| | 2,000 | 0.70 / 0.96 (273) | 0.66 / 0.96 (238) | 0.65 / 0.95 (441) |
| fall-off | 1 | 0.57 / 0.79 (368) | 0.60 / 0.86 (420) | 0.52 / 0.83 (600) |
| | 10 | 0.62 / 0.92 (392) | 0.60 / 0.93 (406) | 0.55 / 0.91 (630) |
| | 60 | 0.57 / 0.87 (308) | 0.73 / 0.95 (280) | 0.73 / 0.95 (532) |
| | 500 | 0.57 / 0.89 (315) | 0.55 / 0.88 (294) | 0.67 / 0.91 (518) |
| | 2,000 | 0.62 / 0.86 (273) | 0.65 / 0.92 (238) | 0.61 / 0.91 (441) |
| F (between-locus share) | 10 | 0.29 / 0.48 (420) | 0.53 / 0.75 (420) | 0.55 / 0.83 (630) |
| | 60 | 0.23 / 0.31 (420) | 0.54 / 0.72 (420) | 0.54 / 0.78 (630) |
| | 500 | 0.23 / 0.31 (385) | 0.53 / 0.71 (420) | 0.55 / 0.80 (630) |
| | 2,000 | 0.21 / 0.27 (385) | 0.54 / 0.72 (420) | 0.58 / 0.81 (630) |

**Thin tail** (fewer estimates: 7 thin one-base and 13 thin two-base classes a run):

| number | 1–2 loci | 3–4 loci | 5–7 loci | 8+ loci (the thick head) |
|---|---|---|---|---|
| slip rate | 0.68 / 0.89 (146) | 0.67 / 0.94 (144) | 0.60 / 0.93 (191) | 0.66 / 0.93 (389) |
| shorter share | 0.64 / 0.89 (146) | 0.68 / 0.94 (144) | 0.61 / 0.92 (191) | 0.60 / 0.92 (389) |
| fall-off | 0.70 / 0.86 (146) | 0.57 / 0.83 (144) | 0.61 / 0.84 (191) | 0.62 / 0.89 (389) |
| F | 0.39 / 0.56 (144) | 0.58 / 0.78 (144) | 0.56 / 0.84 (192) | 0.67 / 0.92 (384) |

The thin-tail figures are pooled over 1 to 2,000 samples; by cohort size they have 22 to 48 estimates
a cell (`tmp/d1/cov_summary.md`), too few to read alone.

What the tables say:

- **Stutter numbers on classes of 5 to 7 loci** cover the truth within two errors 0.82 to 0.95 of the
  time when every class is thin, and 0.84 to 0.93 in the thin tail, against 0.89 to 0.93 for the
  thick classes. Within one error they sit at 0.40 to 0.73, below the honest 0.68 in most cells; the
  thick classes are at 0.60 to 0.66.
- **Classes of 1 or 2 loci are worse where it matters most**: the slip rate at 500 and 2,000 samples
  covers within two errors only 0.69 and 0.59 of the time when every class is thin. In the thin tail,
  where the curve carries the class's number from its thick neighbours, the slip rate covers 0.89
  within two errors.
- **F needs more loci than the stutter numbers.** At 1 or 2 loci it covers within two errors 0.27 to
  0.48 of the time; at 3 to 7 loci 0.71 to 0.84; at 8+ 0.85 to 0.96. Its errors at 3 to 7 loci are too
  narrow by roughly a third at every cohort size.
- Requiring an error to pass the not-placed test (at most 3 on the logit scale, or on the log scale
  for F's concentration) changes no entry by more than 0.02 (`tmp/d1/cov_summary_placed.md`). The
  slipped-read rule (§2.10.4) could not be applied: the bench's fixture line does not record the
  class's expected slipped reads. Commit 7265efd0 has since changed that rule for the fall-off (its
  error is now withheld below one expected slipped read, as the owner ruled); the fixture lines record
  raw errors, so no table here moves.

### 4.3 The allele-length shares and `LOCI_FOR_A_SHARE_ERROR`

The rule withholds a share's error when the fitted share times the class's loci with reads, the
number of loci expected to carry that allele length, is below 3. Coverage by that number, both designs
pooled:

| class loci | samples | below 1 | 1 to 3 | 3 to 8 | 8 or more |
|---|---|---|---|---|---|
| 8+ | 1 | 0.44 / 0.85 (1,334) | 0.64 / 0.92 (183) | 0.63 / 0.92 (158) | 0.57 / 0.79 (149) |
| | 10 | 0.55 / 0.83 (1,302) | 0.69 / 0.93 (217) | 0.66 / 0.91 (141) | 0.70 / 0.96 (164) |
| | 60 | 0.31 / 0.43 (919) | 0.66 / 0.97 (157) | 0.58 / 0.97 (103) | 0.64 / 0.96 (132) |
| | 500 | 0.24 / 0.33 (909) | 0.71 / 0.97 (147) | 0.62 / 0.97 (94) | 0.64 / 0.96 (104) |
| | 2,000 | 0.28 / 0.35 (816) | 0.79 / 0.99 (136) | 0.75 / 0.99 (109) | 0.67 / 0.98 (117) |
| 5–7 | all | 0.37 / 0.66 (51,482) | 0.57 / 0.90 (2,445) | 0.56 / 0.84 (1,401) | — |
| 3–4 | all | 0.36 / 0.69 (32,253) | 0.42 / 0.84 (1,557) | 0.17 / 0.71 (48) | — |
| 1–2 | all | 0.40 / 0.71 (34,023) | 0.57 / 0.82 (215) | — | — |

- **In classes of 8 or more loci, a share carried by 1 to 3 loci is as honest as one carried by
  more**: 0.64 to 0.79 within one error and 0.92 to 0.99 within two, where shares carried by under one
  locus cover 0.24 to 0.55 within one. The rule's line at 3 could move to 1 there.
- In classes of 5 to 7 loci, shares carried by 1 or more loci cover 0.84 to 0.90 within two errors;
  below that size the shares are mostly carried by under one locus and cover poorly.

### 4.4 Reading for `LOCI_FOR_AN_ERROR`

**The reference: classes of 10 or more loci with reads.** The same raw errors on the 2026-10-14
report's default-truth cells (every class has 10 to 100 loci with reads; loci counted from the
committed fixture's record of the same draws), pooled over both motif lengths, the depths and the
seeds of §2 (`tmp/d1/grid_cov.py`, `tmp/d1/grid_cov_now.md`; 2,000 samples × 3 reads only):

| number | 1 sample | 10 | 60 | 500 | 2,000 | all |
|---|---|---|---|---|---|---|
| slip rate | 0.67 / 0.94 (840) | 0.65 / 0.93 (420) | 0.65 / 0.90 (420) | 0.54 / 0.90 (210) | 0.50 / 0.81 (105) | 0.64 / 0.92 (1,995) |
| shorter share | 0.58 / 0.90 (840) | 0.60 / 0.88 (420) | 0.70 / 0.95 (420) | 0.62 / 0.90 (210) | 0.54 / 0.84 (105) | 0.61 / 0.90 (1,995) |
| fall-off | 0.62 / 0.91 (840) | 0.65 / 0.92 (420) | 0.65 / 0.94 (420) | 0.66 / 0.93 (210) | 0.47 / 0.78 (105) | 0.63 / 0.91 (1,995) |
| F | — | 0.65 / 0.89 (420) | 0.66 / 0.93 (420) | 0.71 / 0.93 (210) | 0.62 / 0.92 (105) | 0.66 / 0.91 (1,155) |
| share, carried by < 1 locus | 0.55 / 0.90 | 0.54 / 0.74 | 0.34 / 0.45 | 0.30 / 0.40 | 0.31 / 0.41 | 0.47 / 0.70 (27,951) |
| share, 1 to 3 loci | 0.67 / 0.96 | 0.71 / 0.95 | 0.72 / 0.97 | 0.74 / 0.97 | 0.74 / 0.96 | 0.70 / 0.96 (4,934) |
| share, 3 to 8 loci | 0.58 / 0.86 | 0.69 / 0.95 | 0.66 / 0.96 | 0.69 / 0.96 | 0.68 / 0.96 | 0.64 / 0.92 (2,599) |
| share, 8+ loci | 0.47 / 0.80 | 0.70 / 0.96 | 0.67 / 0.96 | 0.67 / 0.95 | 0.67 / 0.96 | 0.60 / 0.89 (2,421) |

Pooled over every cohort size, against that reference (within one / within two; every class thin):

| number | 1–2 loci | 3–4 loci | 5–7 loci | 10+ loci (grid) |
|---|---|---|---|---|
| slip rate | 0.44 / 0.75 | 0.54 / 0.84 | 0.56 / 0.86 | 0.64 / 0.92 |
| shorter share | 0.54 / 0.84 | 0.58 / 0.92 | 0.56 / 0.90 | 0.61 / 0.90 |
| fall-off | 0.59 / 0.87 | 0.62 / 0.90 | 0.61 / 0.90 | 0.63 / 0.91 |
| F | 0.24 / 0.35 | 0.53 / 0.72 | 0.56 / 0.80 | 0.66 / 0.91 |

**My reading, for the owner:**

- **The stutter numbers' errors could be given from 3 loci with reads.** At 3 to 7 loci the shorter
  share and fall-off cover the truth within two errors 0.90 to 0.92 of the time, as the thick classes
  do (0.90, 0.91); the slip rate 0.84 to 0.86 against 0.92. Within one error all three sit 0.02 to
  0.11 below the thick classes. At 1 or 2 loci the slip rate drops to 0.75 within two errors, and to
  0.59 at 2,000 samples, so those stay withheld.
- **F's error should keep the line at 8**, or have its own: at 3 to 7 loci it covers 0.72 to 0.80
  within two errors, at 1 or 2 loci 0.35.
- **The share rule could move from 3 loci to 1**: in classes of 10 or more loci, shares carried by 1
  to 3 loci cover 0.70 / 0.96, as honest as any group of shares. In thin classes the rule matters
  less, since nearly every share there is carried by under one locus.
- **Lowering the line only helps in runs whose errors were computed.** In one thin-design run in five
  (86 of the 89 at 3 reads) no stutter error was computed for any class (§4.1), and in the thin tail
  that took the thick classes' errors too. For a kimura-shaped cohort that loss is the larger defect;
  `LOCI_FOR_AN_ERROR` does not touch it, and the error-matrix fix in the diagnosis report does.

The reference above is unchanged by the cells first run in D1: its loci with reads come from the
committed fixture, which does not hold 500 × 30 or 2,000 × 10. Their stutter numbers' and F's errors,
counted separately in §3, hold the truth within two errors 0.89 to 0.97 of the time.

---

## 5. One-sample two-base direction numbers with the stutter prior

**The question** (spec §2.6). At one sample the stutter prior pulls the two-base shorter share and
fall-off toward its centre (0.5 and 0.05), and the errors do not widen to cover the pull. C8's review
counted the classes whose number lies more than two of its own raw errors from the truth; this repeats
the count at HEAD, after C9 and C10. An honest error misses by more than two errors about 1 time in 22
(about 3 of 63 classes, 6 of 126).

Runs: `--period 2 --samples 1 --loci-at-most 100 --methods D_tract_prior --fixture-out …`, with
`--stutter-prior-reads 1` (production) and `0` (no prior); 21 classes a run. Script
`tmp/d1/direction.py`; outputs `tmp/d1/dir/w1/`, `tmp/d1/dir/w0/`.

Classes more than two errors from the truth, with the prior | without:

| cell (1 sample, two-base) | seeds | slip rate | shorter share | fall-off | C8's review |
|---|---|---|---|---|---|
| 3 reads, PCR-free | 1–3 | 14/63 \| 14/63 | **31/63 \| 18/62** | 12/63 \| 15/63 | shorter share 31/63 \| 18/62 |
| 3 reads, PCR-free | 1–6 | 33/126 \| 37/126 | 48/126 \| 35/125 | 18/126 \| 15/126 | not run |
| 3 reads, PCR | 1–3 | 7/63 \| 7/63 | 20/63 \| 10/63 | 7/63 \| 5/63 | |
| 3 reads, PCR | 1–6 | 11/126 \| 12/126 | 20/126 \| 21/126 | **28/126 \| 18/126** | fall-off 28/126 \| 18/126 |
| 10 reads, PCR | 1–3 | 4/63 \| 4/63 | 11/63 \| 4/63 | 0/63 \| 0/63 | |
| 10 reads, PCR | 1–6 | 5/126 \| 5/126 | **21/126 \| 12/126** | 5/126 \| 4/126 | shorter share 21/126 \| 12/126 |

- **Every count C8's review gave is reproduced exactly** (bold), so neither C9 nor C10 moved the
  one-sample fits.
- **With the prior, the shorter share misses by more than two errors 3 to 8 times as often as honest
  errors would**: 31 of 63 at 3 reads PCR-free (honest: about 3), 20 of 63 at 3 reads PCR seeds 1–3,
  21 of 126 at 10 reads (honest: about 6). Without the prior the same cells give 18, 10 and 12.
- **Seeds 4–6 at 3 reads PCR go the other way**: the shorter share misses in 0 of 63 with the prior
  and 11 of 63 without, which is why the six-seed count (20 against 21) is level.
- The slip rate is unaffected (11 against 12 of 126, 5 against 5).

---

## 6. What was run

The bench binary (`cargo build --release --features bench-fixtures --example ng_str_fit_bench`,
native on the Mac) was built at 53a3e4f9 for every run but the last 14, and at 8c6ce6f5 for those:
2,000 × 3 PCR-free seeds 2–3 (two-base) and seed 3 (one-base), F = 0.4 seeds 2–3, 500 × 30 and
2,000 × 10. The two builds' fitted numbers are bit-identical (the full reference fixture and an
end-to-end run). Every run went through `tmp/d1/drive.py`, one at a time, at `RAYON_NUM_THREADS=4` on
a shared machine until 2026-10-15 and at 14 threads, alone, for the last 14; the 2,000-sample runs
under the fixture script's memory guard (stop below 8 GiB available; lowest seen 37.6 GiB). Each
run's command, status, wall time and peak memory are one line of `tmp/d1/runs.tsv`; the job lists are
`tmp/d1/*.txt`, written by `tmp/d1/gen_jobs.py`. For planning only: at 14 threads with both fits, a
500 × 30 run took 66 to 81 s and 0.74 to 0.90 GB, a two-base 2,000 × 10 run 230 to 242 s and about
1.0 GB, a two-base 2,000 × 3 run 7 to 21 min and about 0.8 GB.

| what | runs | output | script |
|---|---|---|---|
| held-out comparison, every 2026-10-14 cell below 2,000 samples | 144 | `tmp/d1/grid/*.tsv` | `tmp/d1/heldout.py` |
| held-out comparison, 2,000 × 3 | 15 | `tmp/d1/grid/*n2000_d3*.tsv` | same |
| held-out comparison, 500 × 30 and 2,000 × 10 | 9 | `tmp/d1/grid/*n500_d30*.tsv`, `*n2000_d10*.tsv` | same; `tmp/d1/final/new_cells_cov.py` |
| thin-class coverage | 480 | `tmp/d1/cov/*.json`, `*.truth.tsv` | `tmp/d1/coverage.py`, `cov_summary.py`, `grid_cov.py` |
| one-sample direction numbers | 36 | `tmp/d1/dir/w1/`, `tmp/d1/dir/w0/` | `tmp/d1/direction.py` |

The final tables: `tmp/d1/heldout_final.md` (every held-out cell), `tmp/d1/grid_cov_final.md` (§4.4's
reference), and the seed-by-seed scores against the prototype's, `tmp/d1/final/perseed.py` and
`negatives.py`.

Held-out runs: `--loci-at-most 100 --both-fallbacks`, `--heldout-loci 20` from 500 samples up,
`--methods D_tract_prior` where the counting fit's 2026-10-14 scores are reused (§2.1) and
`--methods B,D_tract_prior` otherwise. For example:

```sh
RAYON_NUM_THREADS=4 target/release/examples/ng_str_fit_bench --period 2 --samples 2000 --depth 3 \
  --seed 2 --truth pcr --loci-at-most 100 --methods D_tract_prior --both-fallbacks --heldout-loci 20 \
  --truth-out tmp/d1/grid/p2_n2000_d3_pcr_s2.truth.tsv --out tmp/d1/grid/p2_n2000_d3_pcr_s2.tsv
```
