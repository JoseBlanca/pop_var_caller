# Predicting a panel's spectrum: three ways, timed and compared

*Measurement report, 2026-08-22. Branch `ng-calling-prior`. Harness:
[`examples/ng_spectrum_projection_cost.rs`](../../../examples/ng_spectrum_projection_cost.rs).*

## Why this was measured

Step D of the genotype prior fits two numbers — the chromosomes' worth of prior belief on the
reference allele and on the alternatives — by searching for the pair whose **predicted** allele
frequency spectrum matches the one the parameter pre-pass measured. The prediction is the objective
of that search, so a fit pays it on the order of a hundred times.

The first version was cubic in the number of samples. At the panel sizes this project tests —
tomato's 63 accessions, GIAB's 3 — a whole fit was 40 milliseconds. At the top of the committed
cohort range, several thousand samples, it was **half an hour**. The question was whether that could
be bought back, and at what cost in accuracy.

## What was compared

- **term-by-term** — the straightforward sum: one exponential per term, nothing stepped, nothing
  skipped.
- **recurrence** — each term written as a beta-binomial weight times a hypergeometric one, both
  genuine probabilities, so the hypergeometric is stepped by an exact ratio instead of
  exponentiated. One exponential per `(branch split, draw count)` pair rather than one per term.
  **The model is unchanged**; this is arithmetic.
- **tail-trimmed** — the recurrence, plus dropping branch splits whose probability falls below a
  given fraction of the likeliest split's.

## The numbers

Release, tomato's fitted diversity of 6 in 10,000, inbreeding coefficient 0.8 — the corner the
projection is aimed at and the one where the branch sum is widest. One prediction:

| samples | term-by-term | recurrence | trimmed at `1e-18` | speed-up | worst class error |
|---|---|---|---|---|---|
| 400 | 43.8 ms | 17.4 ms | **5.8 ms** | 7.6× | 2e-13 |
| 800 | 339.6 ms | 121.8 ms | **29.9 ms** | 11.4× | 4e-13 |
| 1,600 | 2.1 s | 728.9 ms | **179.4 ms** | 11.7× | 4e-13 |
| 3,200 | 12.1 s | 3.7 s | **960.3 ms** | 12.6× | 1e-12 |

**Growth falls from `N^2.95` to `N^2.45`.** In fit terms, at about 160 objective evaluations: 3,200
samples goes from 32 minutes to **2.6 minutes**, 1,600 from 5.6 minutes to 29 seconds.

> **Corrected 2026-08-22, when the fit was built (step D2).** Both numbers in that sentence were
> assumptions and both were low. A fit runs **399** predictions, not 160, and a prediction
> **inside a fit** averages 1.78 s at 3,200 individuals against the 0.96 s measured here at the
> neutral pair — the search spends most of its predictions away from that pair, where the
> branch-tail trim drops fewer splits. Measured end to end: 3.8 s at 400 individuals, 22 s at 800,
> 2.2 minutes at 1,600 and **11.8 minutes at 3,200**. The per-prediction table below is unchanged
> and still correct for what it measures.

## What the accuracy costs, and why the answer is "nothing"

**At `1e-18` the trim is not measurable.** The worst class-by-class disagreement with the
term-by-term sum is 4e-13 at 800 and 1,600 samples — which is exactly the disagreement the
*untrimmed* recurrence has. That figure is floating-point accumulation, not the trim.

Loosening the tolerance does cost, and in proportion to it:

| tolerance | worst class error at 800 samples | time |
|---|---|---|
| `1e-18` | 4e-13 | 29.9 ms |
| `1e-12` | 3e-13 | 24.5 ms |
| `1e-8` | 2e-9 | 19.9 ms |
| `1e-6` | 2e-7 | 17.3 ms |

So the loose end of that range buys another 40% of speed for six orders of magnitude of accuracy.
Not worth it: `1e-18` is shipped.

**The trim's error is bounded, not estimated.** For a fixed number of inbred individuals the classes
are themselves a distribution summing to one, so dropping that split moves no class by more than
its own probability and the whole spectrum by no more than the dropped tail. That is why a
tolerance can be quoted as a guarantee.

## Two defects the measurement found, both in the fast version

**The multiplicative walk must start at its mode.** Started at the low end of its range, the first
weight underflows to zero long before the ones at the mode become small — and since each step
multiplies the last, the whole row then contributes nothing. Measured at 1,600 samples: one class
came back **5.7e-16 against its true 6.1e-7**, and the spectrum lost 3 parts in 10,000 of its mass,
with every entry still finite and non-negative. Nothing downstream could have told. It is pinned by
a release-only test, since the defect needs a panel of that size to appear.

**A `usize` subtraction that release was quietly wrapping.** One step ratio was written
`singles - draws + doubled + 1`, and Rust evaluates that left to right, so the subtraction
underflowed whenever the draw count exceeded the single-chromosome count. In release it wrapped and
the later additions brought it back to the right value; the harness therefore gave correct answers
for the whole first round of measurement. The debug build caught it the moment the code entered the
library. Reordered to `(singles + doubled + 1) - draws`, which cannot go negative.

**The term-by-term sum came out of this well**: at 1,600 samples its total is off by 1.3e-12. It was
never the weak link, and it is kept in the test module as the oracle the fast version is checked
against.

## What was not tried, and why

**Quadrature over the allele frequency.** Gauss–Jacobi with the Beta weight would be exact with
`N + 1` nodes and, with the class distribution obtained by transform rather than convolution, would
reach `N² log N`. It was not built: it needs Golub–Welsch nodes and a transform, the frequency
distribution here is sharply concentrated near zero at realistic diversities — `Beta(6e-4, 1)` puts
essentially all its mass below one in a thousand — and raising a generating function to a large
power at complex nodes is a known accuracy trap. The measured `N^2.45` already puts the whole
committed range inside a once-per-run cost, so the added machinery has nothing to buy yet. Recorded
here so the next person does not have to rediscover the option.

**Binning the classes.** It would have departed from spec §4.1's "over **all** classes including
monomorphic", and it turned out not to be needed.
