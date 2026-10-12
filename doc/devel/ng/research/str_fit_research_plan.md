# Fitting the repeat-tract model: a research plan

**Status:** plan, 2026-10-12, agreed in outline with the owner. Nothing in it is built yet. **No
change to the caller's fit is made until this plan has answered its two questions.** Written to be
picked up in a later conversation, so it states its own context.

The two questions:

1. **How good are the fits?** How close does each method come to the truth where the truth is
   known, how honest are its error bars, and how well does it predict data it was not fitted on?
2. **How much computing do they need?** Time and memory, and how both grow with samples, loci and
   read depth.

---

## 1. The problem, seen fresh

**For each motif length, the reads we observe are the product of hidden variation in several
dimensions plus a measurement error.** The owner's framing (2026-10-12):

- **Along the reference repeat count.** Every parameter below should change smoothly with the number
  of repeat units in the reference: PCR slippage and the germline mutation rate both rise with tract
  length.
- **Between loci with the same reference repeat count.** Each locus has its own allele frequencies
  over repeat counts. Across loci these follow a **typical allele-length shape**, a distribution
  that a few parameters may describe, and a **between-locus share** F that says how far individual
  loci depart from it. F is the Balding–Nichols F_ST with loci in place of populations: at 0 every
  locus carries the typical shape, at 1 each locus is fixed at one length drawn from it.
- **Between samples at one locus.** Samples draw their genotypes from the locus's frequencies, made
  more homozygous by each sample's inbreeding. Population structure among the samples also lives
  here.
- **Stutter, the measurement error.** A read misreports its allele's length by slipping repeats.
  Three numbers describe it: the slip rate, the share of slips that go shorter, and how fast longer
  slips become rarer. They too change smoothly with repeat count.

**What makes stutter separable from true variation is depth.** Stutter acts on each read and the
alleles on each chromosome, so several reads from one chromosome show which is which. At one read
per sample they cannot be told apart. That is why a single low-coverage sample is the hardest case
and must be measured, not assumed.

**In statistical terms this is a hierarchical latent-variable model with measurement error.** The
genotypes are hidden; the loci's frequencies are random effects around a typical shape; the
parameters are smooth functions of repeat count; stutter is a known-form error laid over the
alleles.

### 1.1 What exists today

- **The climb** (`parameter_estimation/joint/ssr_fit.rs`), the fit the caller uses. It fits each
  class (motif length × reference repeat count) on its own, averaging each locus's likelihood over
  its possible frequencies at 256 points, which is several log-likelihood units a locus wrong at 13
  allele lengths. Stutter is then smoothed across repeat count by curves fitted to the per-class
  answers (`spec/str_slippage_level_curve.md`). tomato2 (59 accessions, 32 Mb): 6.7 hours at an
  allele window of ±6, about 22 hours estimated at ±9.
- **The counting fit** (`ssr_fit/composite_fit.rs`), a prototype. Expectation maximisation per
  class, treating each sample at a locus as drawing from the class's distribution independently of
  the other samples there, with F read from pairs of samples. tomato2: 70 s at ±9. kimura (2,169
  accessions, at ±6): 15.5 min, 31.5 GiB peak. No standard errors yet.
- **A six-number allele-length shape** (`scripts/research/str_allele_length_shapes.py`), tested
  against the free shares on kimura at ±6: 2–6 changed calls per 10,000 for homopolymers of 8–11
  repeats, 37–123 at 17–23. Not yet measured at ±9.

---

## 2. What "good" is measured by

Every method is scored on the same five measures, on the same test beds (§3).

| measure | what it answers | where it can be measured |
|---|---|---|
| **Recovery** | how far each estimated parameter and curve is from the truth: bias and spread over repeated draws | simulation only |
| **Error honesty** | whether the reported standard errors match the actual spread of the estimates: the truth should fall within one error about 68 times in 100 | simulation only |
| **Held-out prediction** | how well the fitted model predicts the reads at loci it was not fitted on, in log-likelihood per read | simulation and real data |
| **Calling** | genotype accuracy when the caller uses the fitted parameters | simulation; HG002 and the GIAB trio against their truth sets; the trio's Mendelian consistency |
| **Cost** | wall time, peak memory, and how both grow with samples, loci and depth | everywhere |

**Held-out prediction is the one measure that works on real data without a truth set**, and it is
the fairest way to compare methods with different numbers of parameters. Loci are held out whole,
because samples at one locus share its alleles.

---

## 3. Test beds

### 3.1 A count-level simulator (to build first)

**Draws the census's read counts directly from a known model**, without simulating reads or
alignment, so a whole cohort is drawn in seconds and every parameter's truth is known. It extends
the fit's existing test helper `draw_stratum`, which draws one class, to a whole motif length:

- every parameter given as a curve across repeat count;
- loci per class as the real catalog has them, or chosen;
- per-sample inbreeding, and optionally two subpopulations to test structure;
- several libraries with their own stutter multipliers;
- alleles up to and beyond the recorded ±8.

**The range the caller must cover**, from `CLAUDE.md`, is the grid every method is run on:

| samples | 1 | 10 | 60 | 500 | 2,000 |
|---|---|---|---|---|---|
| **reads a position** | 3 | 10 | 30 | 300 | |

Not every cell is needed for every method; the corners (1 sample at 3 reads, 2,000 samples at 3
reads, 1 sample at 300 reads) are required for all.

### 3.2 The read-level simulator (exists)

`pop_var_caller-tract-ceiling/tmp/stutter_bench/first_pass/`: simulated reads, aligned, with known
stutter rates (`stutter_rates.tsv`) and library multipliers (`libraries.tsv`), PCR and PCR-free,
1 and 50 samples at 5–60×. Used once a method is chosen, to check it end to end through alignment
and the census.

### 3.3 Real cohorts

| cohort | samples | depth | truth |
|---|---|---|---|
| tomato2 | 59 | about 3× | none; held-out prediction only |
| kimura (tomato25) | 2,169 | about 3× | none; the scale test |
| HG002 (`benchmarks/ssr_hg002`) | 1 | about 300× | GIAB repeat-tract truth |
| GIAB trio | 3 | high | GIAB small-variant truth sets; Mendelian consistency |

---

## 4. The avenues

Each is a prototype: an example program over the library, no change to the caller.

### A. The climb, as it is (the reference)

Scored on every measure, so that every other avenue is compared with what ships. It is too slow for
the large cells of the grid; it is run where it finishes in a few hours.

### B. The counting fit per class, then the existing curves

The counting fit replaces the climb inside each class; the curves across repeat count are drawn
from its answers as today. **The question:** does the counting fit match or beat the climb on
recovery and held-out prediction, at a fraction of the cost? It is nearly built; what it lacks is
standard errors (avenue F).

### C. The counting fit with a parametric allele-length shape

The 19 free allele-length shares of each class are replaced by a shape of about 6 numbers, plus one
number for each edge class (which means "9 or more repeats away" and cannot sit on a falling
tail). **The question:** does the shape lose held-out prediction at large classes, and gain it at
thin ones? This includes rerunning the shape comparison at ±9.

### D. Penalised EM: the curves fitted during the fit

The E-step is the counting fit's: weigh each sample's genotypes and count expected reads and
alleles in each class. The M-step fits, for each motif length at once, every parameter as a
penalised spline across repeat count: the three stutter numbers, the shape's numbers, and F. The
penalty's strength is estimated from the data, by cross-validation over loci or by restricted
maximum likelihood, so how smooth each curve is gets measured, not chosen. This is the approach of
generalised additive models for every parameter of a distribution (GAMLSS), inside EM (Green,
1990). **The question:** does borrowing strength during the fit beat smoothing afterwards, at thin
classes and at 1 sample, without smoothing away real kinks such as HG002's 9-repeat homopolymers
(87% below the current curve on 3,875 slipped reads)? It would replace three separate pieces: the
per-class fit, the shape-choosing grid and the blend.

### E. A full Bayesian fit, as a yardstick

The same model written for Stan (Hamiltonian Monte Carlo), run only on small cells of the grid. It
gives the most honest uncertainty available and the exact likelihood, with no shortcut over samples
at a locus. **The question it answers is about the others:** how much do B, C and D lose by their
shortcuts? It is not a candidate for the caller.

### F. Standard errors for the EM avenues

The counting fit's objective treats samples at one locus as independent, so errors computed from it
directly come out too small (spec `str_counting_fit.md` §2.1). Two ways: a bootstrap that resamples
whole loci, and the sandwich formula for a composite likelihood. **The question:** which is honest
(measure 2 of §2) at what cost? Needed by B, C and D alike.

### Other ideas worth one paragraph each in the final report, not a prototype now

- **The exact likelihood by EM**, with each locus's frequencies as a further hidden quantity, in
  place of the counting fit's shortcut. Avenue E says whether the shortcut costs enough to matter.
- **Prior work.** popSTR is recalled as modelling stutter across loci by regression on locus
  features, refitted in turns with the genotypes. That needs checking against the paper before it
  is cited. HipSTR fits stutter by EM locus by locus.

---

## 5. Order, and what runs in parallel

**Phase 0, first and alone: the count-level simulator and one scoring harness** (§3.1, §2). Every
avenue is scored by the same code on the same draws, or their results cannot be compared. This is
small: the per-class drawing exists.

**Phase 1, in parallel**, each in its own worktree and branch, each producing one report:

| track | avenues | depends on |
|---|---|---|
| 1 | **B and F**: the counting fit with honest errors | phase 0 |
| 2 | **D**: penalised EM | phase 0; reuses B's E-step as it stands |
| 3 | **C**: the parametric shape, including the ±9 rerun on kimura | phase 0 for simulation; the kimura run can start now |
| 4 | **E**: the Stan yardstick on small cells | phase 0 |
| — | **A** is scored once, where it finishes, as part of phase 0's harness check | phase 0 |

**Phase 2: the survivors on real data.** Held-out prediction on tomato2 and kimura, calling on
HG002 and the trio, the read-level simulator end to end. Then one comparison report, and the owner
chooses.

**Why parallel works here:** after phase 0 the four tracks share no code they change; each adds its
own module and example program. The risk is four reports measured four ways, which phase 0
removes.

---

## 6. Decisions for the owner

1. **A computing budget.** For example: the fit at kimura's scale (2,169 samples) within an hour on
   one machine with 64 GiB. Without a target, "fast enough" cannot be judged.
2. **How many tracks at once.** Four is the natural split; two or three keeps each review
   manageable.
3. **Whether population structure (two subpopulations in the simulator) is in scope now**, or
   deferred until the single-population model is settled.
