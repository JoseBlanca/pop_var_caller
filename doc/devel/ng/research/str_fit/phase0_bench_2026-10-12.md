# The repeat-tract fitting bench: phase 0

**Status:** built and run, 2026-10-12. Phase 0 of the research plan
(`doc/devel/ng/research/str_fit_research_plan.md` §2, §3.1, §5): a simulator that draws a cohort's
repeat-tract read counts from a known model, and one program that scores every fitting method on the
same draws. This report says how to run it, how a new method plugs in, and what the first run of the
two existing methods found at the plan's three required corners.

## 1. What was built

| piece | where |
|---|---|
| the simulator | `src/parameter_estimation/joint/ssr_fit/bench_fixtures/cohort_draw.rs` |
| the method interface, the curves, the held-out score | `src/parameter_estimation/joint/ssr_fit/bench_fixtures/fit_bench.rs` |
| the bench, one grid cell a run | `examples/ng_str_fit_bench.rs` |
| the grid driver | `scripts/research/str_fit_grid.sh` |
| the summariser | `scripts/research/str_fit_summarise.py` |

All of it is compiled only with the cargo feature `bench-fixtures`; nothing changes in the caller. One
library function changed visibility: the step that draws the curves across repeat count after every
class is fitted is now the public `ssr_fit::draw_and_apply_curves`, and `fit_strata` calls it, so the
bench applies the caller's own curves to any method's answers.

### 1.1 The simulator

It draws, for one motif length, every class (one reference repeat count) at once, as the census would
hand them to the fit: per locus, per sample with at least one read, the reads counted by how many
repeat units they differ from the reference, `−8 … +8`, the end buckets meaning "this far or further".

In the order it draws:

1. **The class's allele-length shares**, from `−14` to `+14` repeat units, so alleles reach past the
   fit's outermost classes (±9, "this far or further") and that meaning is tested.
2. **Each locus's allele frequencies**, a Dirichlet draw around those shares. Its concentration `c` is
   given as the between-locus share `F = 1 / (1 + c)`.
3. **Each sample's genotype**: with chance equal to the sample's homozygote excess its two copies are
   one draw, otherwise two.
4. **Each sample's reads**: a Poisson number with the cell's mean (spanning reads, not coverage), each
   from either copy with even odds and its length reported through the stutter of the sample's library.

Every number is a curve across the repeat count (straight lines between given points, flat past the
ends). The default truth, `CohortModel::tomato_like`, takes its numbers from the tomato cohorts:

| number | one-base runs | two-base runs | source |
|---|---|---|---|
| slippage level | 0.005 at 8 repeats → 0.27 at 16, flat after | 0.0075 at 6 → 0.17 at 14 → 0.21 at 18 | tomato2's counting fit at ±9, class by class |
| shorter share | 0.55 → 0.60 | 0.85 → 0.78 | the same run (0.53–0.62 and 0.74–0.90) |
| fall-off | 0.10 at 8 → 0.30 at 18 | 0.10 at 6 → 0.30 at 18 | the same run, rounded into 0.1–0.3 |
| between-locus share F | 0.7 | 0.7 | the plan; tomato2's two-base runs sit lower, 0.29–0.50 |
| reference share | 0.93 at 8 → 0.50 at 16 | 0.96 at 6 → 0.50 at 16 | the kimura refits of the six-number shape at ±6 |
| homozygote excess | uniform 0.85–0.99 a sample | | tomato accessions are largely inbred lines |

The rest of the shape (the split between shorter and longer alleles, the share one repeat away on each
side, how fast each tail falls) follows the same kimura refits; the doc comment of `tomato_like` gives
every point. A PCR-free variant divides every level by ten. Loci per class are tomato2's (8 to 4,947 a
class, copied into the code as `TOMATO2_LOCI`), or a fixed number, optionally capped.

Several libraries, each with its own stutter multiplier, can be drawn; the evidence either carries the
true multipliers or ones. The model has a place for two subpopulations (`PopulationStructure`), not
built, as the plan defers structure.

**Checked by tests**: over 4,000 loci the share of reads in each bucket is the class's shares mixed
through the stutter, for a library at multiplier one and one at three; across 100,000 draws a locus's
frequency of a length has the Dirichlet's mean and variance; the Poisson depths have their mean and
variance at 0.5, 3, 30 and 300; and the fit's old one-class drawer, `draw_stratum`, now draws through
the simulator and returns exactly what it returned before (a checksum taken before the change).

### 1.2 The interface a new method implements

```rust
pub trait PeriodFit: Sync {
    fn name(&self) -> String;
    fn fit_period(&self, classes: &[StratumEvidence], homozygote_excess: &[f64],
                  allele_span: i32) -> Vec<ClassEstimate>;
    fn draws_its_own_curves(&self) -> bool { false }
}
```

A method is handed **every class of one motif length at once** and returns one `ClassEstimate` per
class: the three stutter numbers, the 19 allele-length shares (`−9 … +9`), the concentration, and a
standard error for any of them it can give. Taking a whole motif length is what penalised EM (avenue D)
needs; a per-class method loops. **Any number may be absent** (`None`), which is the honest answer
where a method cannot estimate it. A method that does not draw its own curves is scored twice: as it
comes (`X_raw`) and after the caller's curves (`X_curves`). To add one, implement the trait in a module
of its own and add one line to `method_named` in the bench.

Two are wired in: **A**, the climb, run as the caller runs it (`fit_strata_on_sample_subsets`); and
**B**, the counting fit with the spread read from pairs of samples (`B_within`: from each sample's own
two copies).

### 1.3 What the bench measures

One run is one cell (motif length, samples, depth, truth, seed) and writes one row per method, class
and number: the truth, the estimate and the standard error.

- **Recovery**: estimate minus truth, summarised over seeds as bias and root mean square error,
  relative to the truth for the level and the concentration.
- **Error honesty**: of the estimates with a standard error, the share within one and within two
  errors of the truth.
- **Held-out prediction**: each class's loci are split into five folds, whole loci; the method is
  fitted on four and scored on the fifth, by the log of the probability of every sample's reads at the
  locus under the class's fitted numbers, the locus's own frequencies integrated out. The truth is
  scored on the same loci, so the reported figure is method minus truth (zero is as good as the
  truth). Per read and per locus.
- **Cost**: wall time of the fit on every locus and of the held-out fits, the process's CPU time and
  peak memory (`/usr/bin/time`, one method a process).

**How the held-out score integrates a locus's frequencies out.** It must not depend on any method's
own machinery, and the climb's 256-point average is several log-likelihood units a locus wrong at 13
allele lengths. So:

- up to two samples with reads at the locus, it is **summed exactly** (each sample's genotype
  probability is a polynomial in the frequencies, whose Dirichlet averages are known); checked against
  a brute-force sum over every allele drawn;
- beyond that, **Monte Carlo**: frequencies drawn from Dirichlets placed where the locus's reads put
  them, mixed with one draw in ten from the class's own Dirichlet; where that does not reach the error
  target within 8,192 draws, a **sequential sampler** that draws the samples' alleles one after
  another through the Pólya urn and so integrates the frequencies out exactly. Draws are added until
  the Monte Carlo error is below **0.05 log-likelihood units a locus**, or 65,536 draws. Both samplers
  agree with brute-force sums within their errors (tests on 2- and 5-sample loci, one of them a locus
  whose reads two alleles explain equally).

The second sampler exists because the first fails on such loci: with F = 0.7 a locus is nearly fixed
at one length, and when reads cannot tell −8 from −9 (both land in the −8 bucket) the frequencies have
two separate peaks. On 10 samples at three reads, 91 loci of 700 missed the target with the first
sampler alone, none with both.

## 2. How to run it

```sh
cargo build --release --features bench-fixtures --example ng_str_fit_bench

# one cell
target/release/examples/ng_str_fit_bench --period 1 --samples 60 --depth 3 --seed 1 \
    --truth pcr --methods A,B --loci-at-most 100 --out rows.tsv --truth-out truth.tsv

# a grid, one method a process, peak memory and CPU time recorded
OUT=tmp/grid METHODS="A B" SEEDS="1 2 3" PERIODS="1 2" CELLS="1:3 2000:3 1:300" \
    LOCI_AT_MOST=100 ./scripts/research/str_fit_grid.sh
# CELLS=full is the plan's grid: 1, 10, 60, 500, 2,000 samples × 3, 10, 30, 300 reads

uv run --no-project --with numpy python scripts/research/str_fit_summarise.py tmp/grid
```

`--truth` takes `pcr`, `pcr_free` or a TOML file whose keys are `CohortModel`'s fields. Other options:
`--loci tomato2|N`, `--repeats 8-16`, `--excess`, `--library-multipliers 1,1.5 [--multipliers-told]`,
`--folds` (0 switches the held-out score off), `--heldout-loci` (scored per class), `--target-se`,
`--max-draws`, `--sampler`.

## 3. The first run: the climb and the counting fit at the three corners

**What was run.** The default truth (PCR stutter), one-base runs (period 1: 14 classes, 8–21 repeats)
and two-base runs (period 2: 21 classes, 6–27 repeats), **at most 100 loci a class** (971 and 878
loci) so that both methods could be run on the same draws, seeds 1, 2 and 3. All on this Mac (18
cores, 64 GiB), with two to four runs sharing the machine, so wall times are inflated by up to about
two-fold; the CPU times are not.

| cell | the climb (A) | the counting fit (B) |
|---|---|---|
| 1 sample, 3 reads | 3 seeds; held-out score for seed 1, period 1 only | 3 seeds, all with held-out score |
| 1 sample, 300 reads | 3 seeds, no held-out score | 3 seeds, all with held-out score |
| 2,000 samples, 3 reads | **not run** | 3 seeds; held-out score for seed 1 |

**Where the climb could not go.** One climb of every class at one sample took 15 to 40 minutes
(880–2,410 s). Its held-out score needs five more fits; the one cell where it was run took 2.9 hours
(10,588 s wall, 22,267 CPU-seconds), so the other held-out cells were left out under the 30-minute
budget. At 2,000 samples it was not attempted: the climb fits a large cohort's classes on its first
256 samples, which at three reads is about 240 times as many sample-locus rows with reads as at one
sample (24,000 against 100 a class), and its work grows with those rows.

The rows are in `tmp/str_fit_grid_2026-10-12b/` (not committed; the commands in §2 reproduce them
exactly, the draws being seeded) and the summariser's full output is committed beside this report as
`phase0_corners_summary_2026-10-12.md`.

### 3.1 Recovery

Each figure is the median, over the classes, of that class's root mean square error over the three
seeds: relative to the truth for the level, absolute otherwise. "raw" is the method as it comes;
"curves" is after the caller's curves across repeat count, which change only the three stutter
numbers.

| cell | number | A raw | A curves | B raw | B curves |
|---|---|---|---|---|---|
| period 1, 1 sample, 3 reads | level (relative) | 0.25 | 0.17 | 0.25 | 0.17 |
| | shorter share | 0.19 | 0.09 | 0.20 | 0.09 |
| | reference-length share | 0.067 | 0.067 | 0.068 | 0.068 |
| | F (truth 0.7), median bias | −0.54 | −0.54 | −0.70 | −0.70 |
| period 2, 1 sample, 3 reads | level (relative) | 0.27 | 0.14 | 0.25 | 0.14 |
| | F, median bias | −0.62 | −0.62 | −0.70 | −0.70 |
| period 1, 1 sample, 300 reads | level (relative) | 0.017 | 0.019 | 0.017 | 0.019 |
| | shorter share | 0.007 | 0.004 | 0.007 | 0.004 |
| | F, median bias | −0.29 | −0.29 | −0.70 | −0.70 |
| period 2, 1 sample, 300 reads | level (relative) | 0.018 | 0.019 | 0.018 | 0.019 |
| | F, median bias | −0.34 | −0.34 | −0.70 | −0.70 |
| period 1, 2,000 samples, 3 reads | level (relative) | — | — | 0.006 | 0.006 |
| | shorter share | — | — | 0.004 | 0.003 |
| | reference-length share | — | — | 0.047 | 0.047 |
| | F, median bias (RMSE) | — | — | +0.008 (0.033) | |
| period 2, 2,000 samples, 3 reads | level (relative) | — | — | 0.005 | 0.005 |
| | F, median bias (RMSE) | — | — | −0.007 (0.055) | |

What the table says:

- **At one sample the two methods recover the stutter equally well** — the level to about 25% a
  class at three reads and 2% at 300 — and the curves cut the error at three reads by a third to a half
  (0.25 to 0.17 or 0.14), but not at 300 reads, where each class already holds its level to 2%.
- **At one sample neither method can estimate F, and neither says so.** The counting fit returns its
  search bound in every class (concentration 10,000, F = 0.0001, where the truth is 0.7): with one
  sample there is no pair of samples to read the spread from. The climb returns F between 0 and the
  truth, biased low by 0.29 to 0.62. The project's rule is that a number a method cannot estimate
  comes back absent, not as a confident value.
- **At 2,000 samples the counting fit recovers everything closely**: the level to 0.5–0.6% a class,
  F to within 0.03–0.06.

### 3.2 Error honesty (the climb only)

The counting fit gives no standard errors. The climb's, pooled over classes and seeds:

| cell | level | shorter share | allele-length shares | F |
|---|---|---|---|---|
| period 1, 1 sample, 3 reads | 0.64 / 0.90 | 0.55 / 0.87 | 0.78 / 0.96 | 0.60 / 0.80 (5 classes) |
| period 2, 1 sample, 3 reads | 0.68 / 0.95 | 0.74 / 0.93 | 0.83 / 0.97 | 0.00 / 0.38 (8) |
| period 1, 1 sample, 300 reads | 0.60 / 0.93 | 0.76 / 1.00 | 0.78 / 0.95 | 0.36 / 0.73 (11) |
| period 2, 1 sample, 300 reads | 0.76 / 0.98 | 0.63 / 0.95 | 0.83 / 0.96 | 0.56 / 0.78 (9) |

Each entry is the share of estimates within one / within two standard errors of the truth; honest
errors give 0.68 / 0.95. The stutter and share errors are about right at one sample. F's are too
small, on the few classes that have one at all.

### 3.3 Held-out prediction

The log-likelihood of the held-out loci's reads under the method's fit, minus the truth's on the
same loci, in log-likelihood units a locus. Zero is as good as the truth.

| cell | A raw | A curves | B raw | B curves |
|---|---|---|---|---|
| period 1, 1 sample, 3 reads, seed 1 (911 loci) | −0.32 | −0.27 | −1.24 | −0.46 |
| period 1, 1 sample, 3 reads, seeds 1–3 | | | −1.2, −10.4, −12.4 | −0.46, −0.68, −0.57 |
| period 2, 1 sample, 3 reads, seeds 1–3 | | | −5.5, −19.2, −34.9 | −1.0, −0.9, −0.8 |
| period 1, 1 sample, 300 reads, seeds 1–3 | | | −54, −43, −41 | the same |
| period 2, 1 sample, 300 reads, seeds 1–3 | | | −76, −87, −111 | the same |
| period 1, 2,000 samples, 3 reads, seed 1 (276 loci) | | | −385 | −386 |
| period 2, 2,000 samples, 3 reads, seed 1 (388 loci) | | | −438 | −438 |

The Monte Carlo error of these figures is zero at one sample (summed exactly) and 0.003 a locus at
2,000 samples.

- **On the one cell where both ran, the climb predicts held-out loci better than the counting fit**:
  0.27 against 0.46 log-likelihood units a locus below the truth after the curves, 0.32 against 1.24
  before them.
- **The counting fit's raw stutter at one sample and three reads is sometimes degenerate**: a class
  with its shorter share at 1.0 and fall-off at 0.0 (8 repeats, seed 1), so a read that slipped the
  other way is all but impossible. Those classes cost tens of units a locus (20 repeats, seed 1: −39
  against the truth's −4.3); the curves remove almost all of it.
- **At 300 reads and at 2,000 samples the counting fit loses 40 to 440 units a locus, and the curves
  do not help**, because the loss is not in the stutter. Two things in its answers can produce it:
  the concentration at its bound at one sample (every locus assumed to carry the class's shares,
  where the truth has most loci fixed at one length), and **allele-length shares of exactly zero** —
  74 of them in period 1 at 2,000 samples, seed 1, some where the truth is 0.0015. A held-out locus
  fixed at a length the fit gave no share costs thousands of units at 2,000 samples. With 80 loci
  a class to fit on and F = 0.7, a length of share 0.0015 is expected at 0.1 of a locus, so zero is
  what maximum likelihood says. At 2,000 samples the loss sits in 7 of 14 classes (81 to 1,836
  units a locus; the other seven lose under 1.5), which fits this explanation; it was not checked
  locus by locus. Free shares fitted on few loci are what the parametric shape (avenue C) would
  replace.

### 3.4 Cost

| cell | method | one fit, s | held-out fits, s | held-out scoring, s | peak memory |
|---|---|---|---|---|---|
| 1 sample, 3 or 300 reads | A | 880–2,410 | 9,556 (one cell) | 0.3 | 35–51 MiB |
| 1 sample, 3 or 300 reads | B | 1–3 | 5–13 | 0.1–0.2 | 10–12 MiB |
| period 1, 2,000 samples, 3 reads (1.85 M rows) | B | 403–1,096 | 2,103 | 1,181 | 3.0–3.8 GiB |
| period 2, 2,000 samples, 3 reads (1.67 M rows) | B | 850–1,237 | 4,651 | 1,897 | 1.7–2.6 GiB |

"Rows" are sample-locus pairs with at least one read. Peak memory is the whole process, the drawn
cohort included; the draw itself takes under a second at every cell.

**The held-out score at 2,000 samples is the bench's own large cost**: about 30 CPU-seconds a locus
(the truth's 276 loci took 466 s on 18 threads), most of it in the sequential sampler, which a
locus needed in more than half the cases (median 10,240 draws). After it, 11 of 276 and 4 of 388
loci were still above the 0.05 error target. Scoring 20 held-out loci a class (`--heldout-loci 20`)
kept a cell to under an hour. At one sample the score is exact and costs nothing.

## 4. Things found on the way

- **The counting fit is slower per row here than on kimura.** On kimura's evidence (86,688 tracts,
  2,169 samples, allele window ±6) the whole run took 931 s; here 971 loci of 2,000 samples at ±9
  took 403–1,096 s. The wider
  window doubles the genotypes (190 against 91); the rest is not explained. Worth a look by track 1.
- **A sampler for a locus's frequencies that ignores multimodality fails silently.** The first
  held-out score drew frequencies from one Dirichlet fitted to the locus; on loci whose reads two
  allele lengths explain equally it reported errors far above target, and fitting the Dirichlet to
  weighted draws made it worse. The sequential sampler fixed it; a Gaussian proposal in log-ratio
  coordinates fitted at the locus's peak would likely make the large cohorts cheaper.

## 5. Not done

- The rest of the plan's grid (10, 60 and 500 samples; 10 and 30 reads), the PCR-free truth, and the
  dinucleotide-realistic F (0.3–0.5): the bench runs them (`CELLS=full`, `TRUTH=pcr_free`, a model
  file), but they were not run.
- The climb's held-out score at three of its four cells, and the climb at 2,000 samples (§3).
- tomato2's full loci per class at 2,000 samples (up to 4,947 a class): the counting fit took 787 s
  for one period-2 fit at its full 1,328 loci, so a held-out cell is about an hour; the corners here
  were capped at 100 loci a class.
- Several libraries with different multipliers: the simulator draws them and a test checks their
  reads, but no method was scored on them.
