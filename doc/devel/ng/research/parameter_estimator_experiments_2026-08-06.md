# What the estimators of ng step 4 actually do — measurements, 2026-08-06

*Research note. Three harnesses, all in `examples/`, all seeded and all reproducible from
a clean checkout. They exist because the parameters ng step 4 emits have **no downstream
check**: nothing a caller does can tell a right error rate from a wrong one, so an error
here is a plausible number nobody notices. Everything below is measurement or derivation;
the design decisions it settled are in
[`../spec/parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md),
[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) and their
architecture companions, which point here rather than repeating the numbers.*

| what was asked | harness | what it settled |
|---|---|---|
| does the multi-library cell key give an unbiased error rate and heterozygosity? | [`ng_multilib_key_harness.rs`](../../../../examples/ng_multilib_key_harness.rs) | the key is fine; the **scoring rule** an earlier draft used was not. Spec §1 |
| is the inbreeding coefficient `F` well behaved? | [`ng_inbreeding_harness.rs`](../../../../examples/ng_inbreeding_harness.rs) | the estimator is sound; its **initialisation** is load-bearing and was unspecified. Spec §6 |
| what does binning the depth cost? | both | the ladder's edges are a **correctness** parameter, not only a memory one — but a good ladder costs 0.05 rungs. §4 |
| what does assuming a heterozygote is a half cost? | [`ng_multilib_key_harness.rs`](../../../../examples/ng_multilib_key_harness.rs) | under 0.5% of heterozygosity at 20 reads a site, up to 10% at three. §5 |
| can the STR stutter accumulator give an unbiased answer? | [`ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs) | **not as specified** — it pools reads across loci, which identifies nothing. Keyed by locus it is exact. §6 |

---

## 1. The method, and why it is not "simulate and compare"

Three properties are wanted of any estimator here: it should be **unbiased**, it should
**approach the truth as data accumulate**, and it should be **as precise as the data
allow**. A single simulated run compared against a single alternative separates none of
them — a difference between two answers mixes bias, sampling noise and grid resolution
into one number, and no amount of care in reading it recovers the three.

**Where the estimator maximises a sum over independent cells, bias can be computed
exactly, with no simulation at all.** Replace each cell's observed count with the cell's
*probability under a known truth* and maximise:

```text
Q(θ)  =  Σ over cells   P_true(cell) · ln L_rule(cell ; θ)
θ*    =  argmax Q(θ)               what the estimate converges to
bias  =  θ* − θ_true               a fixed number, with no sampling error in it
```

`θ*` is the value a misspecified maximum-likelihood fit converges to — the parameter whose
model sits closest in Kullback-Leibler divergence to the truth (White 1982). Because there
are no draws, `bias = 0` and `bias ≠ 0` are decided exactly rather than to within Monte
Carlo noise, and **`bias = 0` is precisely the statement that the estimator is
consistent**. It also runs in milliseconds, so the parameter space can be swept densely
instead of at a handful of hand-picked scenarios.

**A hidden Markov model does not admit this**, because a window's contribution to the
likelihood depends on every other window on its chromosome. For `F` the harness therefore
does the next best thing — a ladder in the number of windows, with repeats — which
separates *converges to the truth* from *converges to something else* from *is merely
noisy*, without pretending to be exact.

**Three habits earned their keep and are worth carrying to the next estimator.**

- **Fit from several starting points, and spread them over every axis the fit can get
  stuck on.** Both harnesses produced a confident wrong finding before this was done, and
  in both cases the wrong finding was the harness's, not the estimator's.
- **Check the algebra before running anything.** For the multi-library key, three
  identities — the scoring rule sums to one, no cell is charged a negative count of reads,
  and with the libraries' error rates equal the rule reproduces the exact likelihood —
  reject a broken rule in one line each, without a fit.
- **Report the *realised* truth of each draw, not the nominal one.** A finite genome does
  not have the `F` its transition rates imply, and comparing against the nominal value
  reads sampling as bias.

---

## 2. The multi-library cell key

### 2.1 The question

When two libraries with different error rates cover a site, the site's likelihood is one
sum over genotypes with both libraries' reads inside it. A cell key of total depth and
total alternative count has forgotten which library produced which read. The design keeps
the alternative reads' library attribution for sites with at most four of them; the
question was what that is worth, and how such a cell should be scored.

### 2.2 The pooled key is not biased — it is blind

Because `p_j`, the chance one read shows a non-reference base at `j` alternative copies, is
a straight line in `ε`, the share-weighted rate over libraries is

```text
Σ_g w_g · p_j(ε_g)  =  p_j(ε̄)        with   ε̄ = Σ_g w_g · ε_g
```

— a single mean. **The pooled key therefore contains `ε̄` and nothing else about the
individual rates**: the likelihood is exactly flat along every combination holding `ε̄`
fixed. The harness confirms it. With the key scored as a proper likelihood, the answer for
one library's rate **moves 23 to 38 rungs of the error-rate ladder depending only on where
the search starts**, while `ε̄` is right to three decimal places at every depth, ratio and
split tried.

### 2.3 The scoring rule an earlier draft used, and what it cost

That draft supplied the missing per-library depths by inventing them — each library given
its average share, `n̂_g = w_g·n`, clamped at zero where that went negative. Two
libraries at 3 reads a site, one erring four times as often as the other, splitting the
reads 90/10:

| what the key keeps, and how it is scored | `ε̄` | the two libraries' own rates | `Hobs` | `π_hom_alt` |
|---|---:|---|---:|---:|
| pooled, by average share | 0 | both return 0.0013 — the four-fold difference is gone | 0 | 0 |
| pooled, as a likelihood | 0 | not identified: drifts 38 rungs with the start | 0 | 0 |
| attribution, by average share | −3.3 | −6.9 / +2.8 rungs | **+64%** | **−74%** |
| **attribution, as a likelihood** | **0** | **0 / 0** | **0** | **0** |

Three readings.

**The average-share plug-in returns the same rate for every library.** It reports `ε̄` for
all of them, so a genuine four-fold chemistry difference comes back as none — the
read-group grain's entire purpose deleted by the scoring rule rather than by the key. An
earlier version of the spec read these numbers as *"the pooled key lands 30 rungs out"*,
and took their refusal to shrink with depth as evidence of how bad the key was. An
unidentified quantity produces exactly that signature: the same answer at 3, 6, 10, 20 and
60 reads.

**Attribution scored by average share is worse than the key it replaces**, and for a
reason unrelated to chemistry — it appears in full on two libraries with the *same* error
rate, where it reports heterozygosity 68% high, the homozygous-non-reference rate 78% low,
and a 2.3-fold difference between libraries that have none. The mechanism is that
`n̂_g = w_g·n` charges a library for reference reads it never had whenever its alternative
count exceeds its average share — at three reads, a library credited with 1.5 when the real
split was 3/0. That happens at **8 sites in 1,000**, and those 8 carry the whole
distortion. It fades with depth: at one rate and an even split, +66% at 3 reads, +15% at 6,
+1.3% at 10, and nothing by 20.

**`EXACT_BREAKDOWN_DEPTH = 4` was a patch for it.** Keeping the whole per-library breakdown
at sites of total depth four or less takes heterozygosity from 66% out to 0.03% at three
reads — but only because 8 sites in 10 at three reads have four reads or fewer and so skip
the plug-in entirely. It was set to the value that hides the defect at the one depth this
cohort sits at, and at six reads with a 90/10 split it still leaves an invented 0.4-rung
difference between two identical libraries.

### 2.4 The rule that works, and it is free

Sum over what the key forgot instead of guessing it. Each read independently picks a
library and then shows the alternative allele or the reference, so the cell — how many
alternative reads came from each library, and how many reads showed the reference in
total — is one multinomial over `G + 1` categories:

```text
                              n!                                                       n−k
L(n, k₁…k_G | θ)  =  Σ  π_j ────────────  ·  Π (w_g·p_j(ε_g))^{k_g}  ·  ( Σ w_g·(1 − p_j(ε_g)) )
                     j       Π k_g! (n−k)!    g                            g
```

Closed form, the same number of cells as the plug-in, and **bias exactly zero in all 31
worlds** — error-rate ratios of 1, 4 and 10, mean depths 3 to 60, even and 90/10 splits,
two libraries and four.

### 2.5 Precision

Ten seeds per rung, on two libraries at 3 reads a site with a 4× rate ratio and a 90/10
split. `sites×` is the extra genome a key needs to match the precision of scoring every
read against its own library.

| candidate | `Hobs` error at 10⁶ sites | at 10⁷ | scatter at 10⁷ | sites× |
|---|---:|---:|---:|---:|
| exact per-library (the ceiling) | −0.29% | −0.22% | 0.41% | 1.00× |
| pooled, as a likelihood | −0.25% | −0.21% | 0.46% | 1.10–1.23× |
| **attribution, by average share** | **+64.16%** | **+64.40%** | 0.52% | 1.62× |
| attribution, as a likelihood | −0.29% | −0.22% | 0.41% | 1.00× |
| attribution at a bound of two, as a likelihood | −0.29% | −0.22% | 0.41% | 1.00× |

The consistency signature is visible directly: scatter falls 12.5% → 5.2% → 0.9% → 0.4%
across 10⁴ to 10⁷ sites, about `1/√N`, and the consistent rules' mean error falls with it.
The plug-in's **flattens at +64.4%** while its scatter keeps shrinking — a floor, not
noise — and it matches the analytic prediction for that world (+64.45%) to 0.05 percentage
points, which validates the analytic calculation against simulation.

The last two rows matter for the design: the coarsened key reaches the exact oracle's
precision, and so does a bound of two alternative reads instead of four. **Once the scoring
is right, the attribution bound is a precision knob that is not currently buying
precision** — and at three reads a bound of two costs 28% fewer cells for the same answer.

### 2.6 The coupled fit — alternating between two tables is consistent

The error rates are read off the read-group table and the genotype frequencies off the
whole-sample one, and the two are coupled: a higher error rate explains the same
alternative reads as less real variation. The design resolves that by alternating. **That
is not a climb on any single objective** — it is a fixed point of two estimating
equations, one per table — so whether it even lands on the truth had to be checked rather
than assumed.

Both tables' cells are weighted by their exact probabilities under a known truth, so what
the loop converges to is what an infinite genome would give. Started deliberately away
from the answer — every error rate at three times the truth, every genotype frequency at
half of it — the fixed point is the truth in **all 25 worlds**: `0.000` rungs on every
error rate, `0.000%` on heterozygosity and on the homozygous-non-reference rate, across
ratios of 1, 4 and 10, depths 3 to 20, even and 90/10 splits, two libraries and four.

**Why that is not luck.** Each table's likelihood is correctly specified for its own
entries. A read-group entry holds one library's depth and alternative count at a site, and
its marginal is

```text
P(m, k)  =  PoissonTruncated(λ·w_g ; m) · Σ_j π_j · Binom(k ; m, p_j(ε_g))
```

— the genotype is still drawn once for the site and still enters through the same mixture.
So each block's score is an unbiased estimating equation at the truth, and the truth is a
fixed point of the pair.

**Splitting a site between two entries costs nothing but precision.** Fitting *everything*
— both the rates and the frequencies — from the read-group table alone is also exactly
unbiased in all 25 worlds. This matters because §1 of the spec rejects a *windowed*
histogram keyed per read group on the grounds that splitting a site lets its two entries
draw independent genotypes, and it is worth being clear that the objection is about the
windowed object rather than about splitting as such. In the standard vocabulary the
product over split entries is a **composite likelihood**: every factor is a true marginal,
so the estimator stays consistent and what the split throws away is the dependence between
a site's entries, which is precision.

**Practical note on the stopping rule.** Convergence to the fixed point is *linear*, as
alternating schemes are: many worlds ran the full 200 iterations without meeting a
movement tolerance of 10⁻¹², yet every one of them was already at the truth to better than
a thousandth of a ladder rung. A tolerance fine enough to be interesting is therefore far
finer than the answer needs. Stop on **the ladder rung ceasing to move**, cap the loop, and
keep the best-scoring iterate.

### 2.7 What is not measured here

Depth binning — every number in §2 is at an exact depth. It is measured in §4, which found
the scoring rule stays exactly unbiased under a good ladder and loses up to one rung of the
error-rate ladder under a bad one.

---

## 3. The inbreeding coefficient `F`

### 3.1 The setup

A tomato-shaped genome: 12 contigs, 8,004 windows of 100 kb, 100,000 covered sites a
window, 3 reads a site, error rate 0.001, one heterozygote per kilobase outside a run and a
homozygous-non-reference rate of 6 per kilobase. Runs are generated as the spec describes —
one allele draw doubled, so the homozygous-non-reference rate *rises* inside a run rather
than the heterozygotes merely being suppressed.

### 3.2 The estimator is sound

`F` recovers the **realised** autozygous fraction of each drawn genome to four decimal
places:

| true `F` | what the genome realised | `F` returned |
|---:|---:|---:|
| 0.05 | 0.0381 | 0.0381 |
| 0.15 | 0.1349 | 0.1349 |
| 0.30 | 0.3073 | 0.3073 |
| 0.60 | 0.6016 | 0.6016 |

### 3.3 The robustness claim holds — and had never been tested

Spec §6.2's second reason for choosing the runs estimator over the ratio one is that a
uniform floor of spurious heterozygotes — collapsed paralogs, mismapping — lifts both
states together and cancels out of the gap between them. Adding such a floor to both
states, in multiples of the real outside heterozygosity of 1 per kilobase:

| false heterozygotes added | realised `F` | `F` returned | fitted inside het | fitted outside het |
|---|---:|---:|---:|---:|
| none | 0.2817 | 0.2817 | 0.000052 | 0.000999 |
| 1 per kb | 0.3027 | 0.3028 | 0.001042 | 0.002001 |
| 3 per kb | 0.2629 | 0.2634 | 0.003055 | 0.003999 |
| 5 per kb | 0.3601 | 0.3592 | 0.005044 | 0.006009 |

Both states rise with the floor and `F` does not move. The claim is true. What the ratio
estimator would read on the same data — the whole-genome heterozygosity — rises from
0.000715 to 0.005715, an eight-fold inflation, so the contrast §6.2 draws is real.

### 3.4 The defect: the initialisation is load-bearing

**The table above is what you get from starting points that disagree about how far apart
the two states are.** From starting points that disagree only about `F`, all sharing one
guess at the separation, the same data gives:

| false heterozygotes added | realised `F` | starts spread over `F` only | starts also spread over the separation |
|---|---:|---:|---:|
| 2 per kb | 0.3255 | 0.3255 | 0.3255 |
| **3 per kb** | **0.2629** | **0.0000** | 0.2634 |
| 5 per kb | 0.3601 | **0.0000** | 0.3592 |

Once a floor lifts both states, the truth has them only 1.3-fold apart in heterozygosity.
A start guessing the inside state at a tenth of the outside rate then fits every window to
the outside state on the first pass, empties the inside one, and drives the rate of
entering a run to zero. **Every start made the same wrong guess, so "keep the best-scoring
fit" had nothing better to pick**, and the fit reported convergence.

The failure leaves a fingerprint: the inside state's fitted heterozygosity comes back as
**exactly its starting value** — 0.000400, 0.000450, 0.000500, 0.000550, 0.000600 across
the collapsing rows, precisely a tenth of each row's starting outside rate. Nothing was
ever assigned to that state, so nothing ever updated it.

**But `F = 0` is also a correct answer for an outcrossing species, and that case leaves the
same fingerprint.** So the two cannot be told apart from the fitted values alone. They can
be told apart from **how the search went** — whether starts that reached different `F`
scored measurably differently — which is why the fit must report the spread across its
starting points and not only the winner.

### 3.5 The chain does no work at 100 kb

A window's state is inferred from its own reads *and* from its neighbours, through the
transition rates. `undecided` is the share of windows whose posterior landed between 0.01
and 0.99 — the only ones where a neighbour could have mattered. `shuffle Δ` refits after
shuffling the window order within each contig, which destroys every run while preserving
every window's contents.

| sites per window | undecided | `F` change on shuffling |
|---:|---:|---:|
| 100,000 (100 kb at 3 reads) | **0.00%** | **0.0000** |
| 10,000 (10 kb) | 5.17% | −0.0014 |
| 3,000 | 19.34% | −0.0043 |
| 1,000 | 48.44% | +0.0149 |
| 300 | 83.95% | +0.2469 |

At the specified window size **not one of 8,004 windows is undecided**, and shuffling the
entire genome changes `F` by zero. Every window is classified on its own evidence, so the
model is a two-component mixture over windows wearing a hidden Markov model's clothes. Two
claims in the spec do not survive that: that forward–backward *"pools across windows... which
lets a lone quiet window between noisy ones read as a fluctuation"* — there are no such
windows — and that shuffling is a diagnostic, since at 100 kb it returns "no difference"
whether the runs are real or not. It becomes a diagnostic below about 1 kb.

Finer windows do not hurt accuracy. Holding the genome at 800 Mb and cutting it three ways,
with the 3-per-kilobase floor of §3.4 present:

| sites per window | windows | realised `F` | `F` returned |
|---:|---:|---:|---:|
| 100,000 | 8,004 | 0.2886 | 0.2884 |
| 10,000 | 80,040 | 0.3477 | 0.3463 |
| 3,000 | 266,796 | 0.2289 | 0.2278 |

So a finer grain is available if it is ever wanted, at ten to thirty times the windows —
which land on the most expensive accumulator in step 4.

### 3.6 The floor: what `F` returns on a genome with no runs

Not zero. A two-state model can always raise its score a little by calling some windows
*inside a run* on the strength of ordinary sampling wobble in their heterozygote counts.
Eight seeds each, on a genome generated with no runs at all:

| windows | mean `F` returned | worst seed |
|---:|---:|---:|
| 1,200 | 0.226 | **0.84** |
| 4,800 | 0.017 | 0.086 |
| 19,200 | 0.0040 | 0.021 |
| 76,800 | 0.0016 | 0.0046 |

It shrinks with windows, so it is imprecision rather than bias — but it is the estimator's
**resolution**, and it should be read as one: about **0.01 at tomato's 8,004 windows** and
about **0.003 at a human genome's 31,000**. An `F` below that means *nothing detected*, not
a small autozygous fraction. It also means a run over a few hundred windows — a development
fixture, a restricted region set — cannot estimate `F` at all: at 1,200 windows a genome
with no runs returned 0.84 on one seed of eight.

### 3.7 Two definitions of `F`, and they differ

The spec calls both of these `F`: the chain's stationary inside-probability
`tAZ/(tAZ + tHW)`, a property of the fitted *model*; and the coverage-weighted posterior
occupancy, a property of the *data*. For a finite genome they differ by ordinary sampling.

| true `F` | realised | posterior occupancy | stationary ratio |
|---:|---:|---:|---:|
| 0.05 | 0.0381 | **0.0381** | 0.0338 |
| 0.30 | 0.3073 | **0.3073** | 0.3181 |

11% apart in relative terms at `F` = 0.05 and 3.5% at `F` = 0.30. The posterior occupancy
recovers the realised value exactly, and it is also the quantity the caller's prior asks
for — whether *this* individual's two copies at *this* locus descend from one ancestral
copy, not what a hypothetical genome from the same pedigree would average.

### 3.8 One thing that is a proof rather than a measurement

If the two states carry identical genotype frequencies, every window's emission is the same
under either state, the observed sequence is independent draws from one distribution, and
`P(data | θ)` does not contain the transition rates at all — the likelihood is **exactly
flat in `F`**. The harness reproduces the flatness: nine starting points score identically
to eleven significant figures.

**The ridge is nonetheless unreachable.** A fit never sits at coincidence; it starts with
the states separated, every window fits one of them better, the other empties, and `F`
collapses to zero. Expectation-maximization walks away from the ridge rather than along
it, so this degenerate case does not produce an arbitrary `F` in practice.

### 3.9 What is not measured here

The noise floor of §3.6 at *fine* windows — many windows, little evidence in each, which is
exactly where a two-state model splits on noise. Everything said about finer grains in
§3.5 assumes it behaves, and that should be measured before the window is moved down.

Depth binning of the window emission, which §4.4 measures: `F` does not move at all.

**Every table in §3 is at 3 reads a site.** §4.4 and §6 add 20 and 60, and the second of
those turned up two defects in this harness that only fire above about 40 reads. The §3
tables are unaffected and were re-run against the fixes to confirm it.

---

## 4. What binning the depth costs

### 4.1 The question, and why it could not wait

The accumulator does not keep a site's exact depth. It keys a cell by a depth **bin** — one
bin per depth at the bottom, widening geometrically above, to a cap — and scores every cell
at the mean of the exact depths that fell in it
([`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §2.2). Three
claims rested on that being harmless and none was measured: the multi-library score of §2,
which is written for an exact `n`; the pooled score; and the window emission the runs model
of §3 reads.

It could not wait for the implementation because a bad answer changes a **type** rather than
a constant — either the attributed cell keeps an exact depth, or the exact-per-depth region
runs much further up, or the score needs the bin's whole depth distribution instead of its
mean, and those are three different accumulators.

A ladder is written **`exact≤E, B bins, cap C`**: exact integers to `E`, then `B − E − 1`
geometrically widening bins up to `C`. The **exact ladder** — one bin per depth — is the
control: under it every cell's mean depth is its own depth, so a binned fit must reproduce
the unbinned answer, and any bias it reported would be the harness's.

### 4.2 The checks that need no fit

At 20 reads a site, two libraries at a 4× rate ratio and a 90/10 split:

| ladder | bins | widest bin | a site's score moves, no alt read | …showing an alt read | mass scored below its own alt count |
|---|---:|---:|---:|---:|---:|
| exact (control) | 126 | 1 | **0** | **0** | 0.0000% |
| exact≤8, 16 bins, cap 124 | 16 | 40 | 1.9 × 10⁻⁷ | 8.7 × 10⁻³ | 0.0000% |
| exact≤8, 20 bins, cap 124 | 20 | 27 | 1.4 × 10⁻⁷ | 5.9 × 10⁻³ | 0.0000% |
| exact≤8, 16 bins, cap 300 | 16 | 121 | 2.2 × 10⁻⁷ | 1.1 × 10⁻² | 0.0000% |

"A site's score moves" is how far one site's contribution to the objective shifts when its
exact depth is replaced by its cell's mean, in nats, averaged over sites. It is taken on the
mixture over the three genotypes rather than on the components, because a component that
shifts by hundreds of nats while sitting 10⁻⁴⁰⁰ below its neighbours shifts nothing a fit can
see. Measured on the components instead — which was the first thing this column did — the
ladders above read about 15 nats a site at 60 reads, where the fits they produce differ from
the exact answer by 0.03 rungs. The number was measuring the deep, empty tail of the
homozygous-non-reference term and nothing else.

**Where binning lands is on the sites that show an alternative read**, by four to five orders
of magnitude, and those are one site in a few hundred. At 60 reads a site both columns fall
below 10⁻⁹: the genotype is certain at that depth whatever the exact `n`.

**Scored at its own cell's mean, no cell is ever charged a negative number of reference
reads** — zero truth mass, on every ladder and at every depth. §4.5 measures the variant
that is.

### 4.3 The bias, over twenty worlds

Error-rate ratios of 1 and 4, mean depths 3, 10, 20, 30 and 60, even and 90/10 read splits.
Each cell weighted by its exact probability under the truth, so these are what an infinite
genome returns, with no sampling noise in them. The worst of the twenty is shown, with the
world it came from; `Hobs` and `π_hom_alt` are relative errors, `ε̄` is in rungs of the
quarter-Phred error-rate ladder (5.9% each).

The last column is how many cells the 20-read world's own data occupies, which is what the
harness enumerates — not the accumulator's table size. A table sized for the whole ladder is
583 cells whatever the sample's depth, since a bin's row must be as wide as its deepest site's
alternative count.

| ladder | bins | widest bin | worst `ε̄` | worst `Hobs` | worst `π_hom_alt` | where the worst is | occupied at 20 reads |
|---|---:|---:|---:|---:|---:|---|---:|
| exact (the control) | 126 | 1 | **0.000** | **0.00%** | **0.00%** | — | 2,542 |
| exact≤4, 16 bins, cap 124 | 16 | 33 | 0.466 | 1.63% | 2.23% | 10 reads | 360 |
| exact≤8, 16 bins, cap 124 | 16 | 40 | 0.545 | 1.33% | 1.83% | 10 reads | 382 |
| exact≤16, 20 bins, cap 124 | 20 | 61 | 0.981 | 3.30% | 5.53% | 20 reads | 407 |
| exact≤12, 20 bins, cap 124 | 20 | 35 | 0.103 | 0.30% | 0.46% | 20 reads | 425 |
| **exact≤8, 20 bins, cap 124** | 20 | 27 | **0.054** | **0.23%** | **0.30%** | 10 reads | 495 |
| exact≤8, 24 bins, cap 124 | 24 | 21 | 0.019 | 0.11% | 0.14% | 10 reads | 609 |
| exact≤8, 16 bins, cap 300 | 16 | 121 | 1.038 | 4.72% | 7.96% | 10–20 reads | 283 |
| exact≤8, 20 bins, cap 300 | 20 | 84 | 0.190 | 0.64% | 0.88% | 10 reads | 419 |

The exact ladder returns `0.000` rungs and `0.00%` in all twenty worlds, so the binned
machinery reduces and the numbers above are the binning's.

**Four readings.**

**The edges are a correctness parameter, not only a memory one.** Holding the cap at 124 and
moving from 16 bins to 20 takes the error rate from 0.55 rungs to 0.05 and the
homozygous-non-reference rate from 1.8% to 0.3%, for 30% more cells. That is a tenfold change
in the answer from a choice the architecture doc had left `OPEN:` on the premise that it only
buys memory.

**The cap competes for bins, and 300 is not free.** At 16 bins, raising the cap from 124 to
300 doubles the error-rate bias (0.55 → 1.04 rungs) and quadruples the
homozygous-non-reference one (1.8% → 8.0%) — **on data where no site is deeper than 125**.
The extra reach is spent on depths nothing occupies, and it is paid for by the bins covering
the depths everything occupies.

**The band where it bites is the ordinary whole-genome one: 10 to 30 reads.** At tomato's 3
reads a site every ladder whose exact-per-depth region reaches 8 is inside 0.004 rungs,
because 97 sites in 100 sit at depth 6 or below and are never binned at all. At 60 reads
every ladder is inside 0.04 rungs, because the genotype is certain. In between, the bins are
already wide and the genotype is not yet decided, and that is where the two meet. *An
exact-per-depth region of 4 rather than 8 costs 0.15 to 0.21 rungs even at 3 reads*, because
the Poisson tail of a 3-read sample reaches depth 15.

**The homozygous-non-reference rate is the most sensitive of the three outputs**, running
about 1.5 times the heterozygosity error in every row — which matters because §5 of the spec
argues it is the output that carries the most for a landrace far from the reference.

**On this evidence the ladder to build is `exact≤8, 20 bins, cap 124`**: 0.054 rungs and
0.3% at its worst, on 495 cells at 20 reads against the exact table's 2,542.

### 4.4 The runs model's window emission

Spec §6.1 says of this emission that binning "makes it an approximation" and that "scoring
each cell at its own mean depth is what keeps the error small". The same drawn genome,
refitted with the same nine starting points, changing only how a window's sites are keyed
and scored — 3,600 windows of 100,000 sites, with window coverage alternating over 0.6×,
1.0× and 1.6× the stated mean so that a cell's per-window mean depth and its whole-sample
mean are genuinely different numbers:

| mean depth | ladder | cell mean taken | cells per window | `F` returned |
|---:|---|---|---:|---:|
| 3 | *realised truth* | | | **0.4125** |
| 3 | exact (control) | exact | 434 | 0.4128 |
| 3 | exact≤8, 16 bins, cap 124 | per window | 132 | 0.4128 |
| 3 | exact≤8, 16 bins, cap 124 | whole sample | 132 | 0.4128 |
| 20 | *realised truth* | | | **0.3408** |
| 20 | exact (control) | exact | 3,320 | 0.3408 |
| 20 | exact≤8, 16 bins, cap 124 | per window | 281 | 0.3408 |
| 20 | exact≤8, 16 bins, cap 300 | whole sample | 266 | 0.3408 |
| 60 | *realised truth* | | | **0.2236** |
| 60 | exact (control) | exact | 15,050 | 0.2236 |
| 60 | exact≤8, 16 bins, cap 124 | per window | 458 | 0.2236 |
| 60 | exact≤8, 16 bins, cap 300 | whole sample | 466 | 0.2236 |

**`F` does not move at all.** Every ladder, both mean grains, all three depths, flat coverage
and varied: **the same four decimal places as the unbinned fit**, which is itself within
0.0003 of the genome's realised autozygous fraction. At 60 reads that is 15,050 cells a window
reduced to 458 — a 33-fold saving — for an identical answer. The claim in spec §6.1 holds, and
it holds with room to spare.

**What does move is the two states' fitted heterozygote rates**, and only slightly: the
inside rate `h` moves from 0.000050 to 0.000054 and the outside rate `Hout` from 0.001001 to
0.000953 on the coarsest ladder tried — 8% and 5%. Those are reported diagnostics rather than
emitted parameters, and `F` reads only the gap between them, which is 20-fold and survives.

**Taking the mean per window rather than over the whole sample is worth about 8% of `h`**, and
nothing of `F`. That is the smaller of the two grains the architecture could have specified,
and it specifies the more accurate one, so nothing is owed here.

### 4.5 Scoring a cell at its bin's mean instead of its own — the doc overstates this

The architecture doc (§2.2) argues that a **bin** mean is not merely coarser but unbounded: a
cell whose alternative count exceeds its bin's mean depth is charged a negative number of
reference reads, its term `(ε/3)^(n−k)` "diverges as `ε → 0`", "the profile scan's objective
becomes unbounded and every fit rails to the ladder's floor".

**The mechanism is real and the consequence is not.** Fitting the same worlds with the bin
mean in place of the cell mean:

| | with the cell's own mean | with the bin's mean |
|---|---:|---:|
| truth mass scored below its own alt count | 0.0000% | 0.28% – 0.32% |
| `ε̄`, rungs | 0.00 – 0.55 | **−0.50 to −5.23** |
| `Hobs` | −1.3% to +0.4% | **+6.1% to +17.6%** |
| `π_hom_alt` | −0.1% to +1.8% | **−10.1% to −29.3%** |

The bin mean is badly wrong — it moves the homozygous-non-reference rate by up to 29% and
the error rate by 5.2 rungs, downward, which is the direction the doc predicts. But it does
**not** rail: 5.2 rungs below the truth is 0.74 times the true rate, where the ladder's floor
is 80 rungs below. The 0.3% of sites whose term grows as `ε` falls are outweighed by the
sites showing one or two alternative reads, whose terms fall faster, so the objective stays
bounded and the fit lands at a wrong but finite value.

**That makes the failure worse to live with, not better.** A railed fit announces itself —
`ScanResult::argmax_at_ladder_end` exists precisely to catch it. A rate 26% low with a
homozygous-non-reference rate 29% low announces nothing. The design decision the doc reaches
is right; the reason it gives would let a reviewer accept a rail flag as sufficient
protection, and it is not.

### 4.6 What is not measured

**Sites deeper than the cap.** No world here reaches one: at mean depth 60 the deepest site
the Poisson truncation keeps is 125, against caps of 124 and 300. So the cap enters these
numbers only through how far it stretches the geometric region. What a site *deeper* than the
cap costs is the subsampling rule
([`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §2.2's
hypergeometric draw), and this harness does not implement one — which leaves HG002's 300×
untested on the one mechanism that only fires there.

**Ladders between the ones tried.** The recommendation of §4.3 is the best of eight, not a
search. The gradient is clear enough — the widening ratio between bins is what the answer
tracks, 1.49 at 16 bins and 1.28 at 20 — but no ladder was tuned.

**Binning at more than two libraries, and at ploidy above two.** Both fits here are two
libraries and diploid.

---

## 5. What assuming a heterozygote is a half costs

### 5.1 The question

The model says a read at a heterozygote shows the alternative allele with probability
`½ + ε/3`. Reads carrying the alternative allele map slightly less often than
reference-carrying ones, so the truth sits nearer 0.47–0.49. Spec §8 leans toward adopting a
fitted per-read-group constant in place of the `½`; spec §11.3 makes the decision turn on
whether a fitted value departs from a half by more than its standard error on real data,
which needs a pipeline that does not exist.

**The prior question is answerable today, and it is the one that decides whether to build the
parameter at all:** generate at a true heterozygote allele balance `b`, fit with the model
that assumes a half, and read the exact bias. If the misspecification costs nothing, no
finding on real data would make the parameter worth having.

One library, so nothing about the multi-library key of §2 is in these numbers, and exact
depths throughout, so nothing about §4's binning is either. **`b = 0.50` is the control and
returns exactly zero in every row.** Every fit was run from three starting points — the
truth, three times it, and a third of it — and all three landed on the same answer to three
decimal places, so none of these zeros is a search that started on its own answer.

### 5.2 The numbers

At tomato's measured rates — 1 heterozygote and 6 homozygous-non-reference sites per
kilobase — `ε` in rungs of the quarter-Phred ladder and the two frequencies in relative
error:

| mean depth | `b` = 0.49 | 0.47 | 0.45 | 0.44 |
|---|---|---|---|---|
| **3** | 0.05 rungs, `Hobs` −1.4% | 0.16 rungs, **−4.5%** | 0.27 rungs, **−8.0%** | 0.32 rungs, **−9.8%** |
| **6** | 0.03 rungs, −0.9% | 0.08 rungs, −2.9% | 0.14 rungs, −5.2% | 0.17 rungs, −6.5% |
| **10** | 0.01 rungs, −0.4% | 0.03 rungs, −1.4% | 0.06 rungs, −2.7% | 0.08 rungs, −3.5% |
| **20** | 0.001 rungs, −0.04% | 0.002 rungs, −0.2% | 0.005 rungs, −0.4% | 0.007 rungs, −0.5% |
| **60** | 0.000 rungs, −0.00% | 0.001 rungs, −0.00% | 0.001 rungs, −0.00% | 0.001 rungs, −0.00% |

The homozygous-non-reference rate never moves more than 0.4% anywhere in that table.

**The error-rate bias scales with how heterozygous the genome is; the heterozygosity bias does
not.** At ten times the heterozygote rate — 10 per kilobase — the same worlds give, at 3
reads: 0.45 rungs at `b` = 0.49, 1.34 at 0.47, 2.23 at 0.45 and 2.67 at 0.44, while `Hobs`
stays at −1.3%, −4.2%, −7.4% and −9.2%. `ε` absorbs the misfit in proportion to how many
sites are misfitted; the heterozygote count loses a fixed *share* of itself.

### 5.3 The reading

**All of the cost is at low depth, and the mechanism is class confusion rather than
arithmetic.** At 3 reads a heterozygote usually shows one or two alternative reads of three;
lowering `b` moves mass from "two of three" to "one of three" and to "none of three", and a
site showing no alternative read is indistinguishable from a homozygous-reference one. The
heterozygotes that fall off the bottom are simply lost, so `Hobs` comes back low. At 20 reads
a heterozygote shows nine or ten of twenty and never zero, no site changes class, and the
misfit inside the heterozygote term is absorbed with no visible cost.

**Against the decision rule agreed in advance** — under one rung in `ε` and under 5% relative
in both frequencies, across `b` from 0.44 to 0.50:

- `π_hom_alt` passes everywhere: 0.4% at worst at tomato's heterozygote rate and 3.9% at ten
  times it, against a bar of 5%.
- `ε` passes at tomato's heterozygote rate — 0.32 rungs at worst, a third of the bar — and
  **fails at ten times it**, reaching 2.67 rungs at 3 reads and `b` = 0.44. It is inside one
  rung there for `b` ≥ 0.48.
- `Hobs` **fails at 3 reads for `b` ≤ 0.46** — 6.2% at 0.46, 8.0% at 0.45, 9.8% at 0.44 — and
  is just inside the bar at `b` = 0.47, at 4.5%. It passes at every depth of 6 and above, and
  the failure is the same size at both heterozygote rates.

So the rule does not close §11.3 as "no". But it does not open it as written either, because
what it identifies is narrower than the question §8 asks: **the parameter earns its keep only
on shallow samples, and there its whole effect is on heterozygosity.** A cohort at 20 reads a
site or more may adopt `½` with a clear conscience; tomato's 3 reads is exactly the case where
it cannot.

### 5.4 What is not measured

**Whether `b` really is 0.47–0.49 on our data.** That is §11.3's own question and it still
needs the pipeline. What is settled is that the answer matters at 3 reads and does not at 20.

**A model that also lowers the depth at homozygous-alternative sites.** Reference bias does
not only skew a heterozygote's reads; it lowers coverage wherever the sample is
non-reference. That is a richer misspecification than the one §8 proposes and it was
deliberately left out — these numbers are for the heterozygote term alone.

**Interaction with the multi-library key and with depth binning.** Both were held exact here,
on the reasoning that three separate approximations should be priced separately before being
priced together. Nothing says they add.

---

## 6. The STR stutter accumulator

*Harness: [`ng_str_stutter_harness.rs`](../../../../examples/ng_str_stutter_harness.rs). Same
method as §2 — every cell weighted by its exact probability under a known truth — applied to a
different accumulator. The design it measures is
[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.1.*

### 6.1 The question, and the truth every world is generated from

The STR path estimates, per (read group × stratum), **how often a read shows a different
number of repeats than the allele it came off, which way it moves, how far, and a per-base
substitution rate**. A stratum is one *(motif period, repeat count)*.

Every world below is generated from a stratum shaped like tomato's dinucleotides at six or
more repeats, as
[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §3 and §5 measure it:
**2.0% of reads slip; 17 in every 100 of those gain a repeat rather than lose one — a 4.9-fold
asymmetry; and 9 in every 100 slipped reads take a second step.** What varies between worlds is
the locus population — how many loci carry an allele other than the reference length, and how
many carry two different ones — and the depth.

**A correction the harness forced before any world was run: this path fits four numbers, not
three.** §1 of the spec, the summary table in
[`parameter_prepass.md`](../spec/parameter_prepass.md) §1, and §3.1's "scan its three noise
parameters" all count the substitution rate, the direction split and the fall-off — and omit
**how often a read slips at all**. That is the quantity §4 stratifies by, §1.1 holds monotonic
across strata, and §5 tabulates at 0.091%, 0.170% and 2.006% by repeat-count band. A model with
no level cannot express any of them.

### 6.2 The cell has to be a locus. Pooling reads across loci identifies nothing

**The accumulator as specified pools reads across loci**: one tally per (read group, period,
repeat count), holding how many *reads* landed at each whole-repeat offset (§4.1, and the
five-object table in [`parameter_prepass.md`](../spec/parameter_prepass.md) §5.1). A read
carries no genotype — it drew one allele and then slipped — so what such a tally holds is the
**allele spectrum convolved with the slippage kernel**, and recovering the kernel means undoing
a convolution with both halves unknown.

Each row fits the level from four starting points that disagree about all three slippage
parameters at once. `spread` is the ratio between the highest and lowest answer: the column
that says whether the fit found something or merely stopped.

| allele lengths present | loci off the reference length | allele spectrum | level | direction split | spread |
|---|---:|---|---:|---:|---:|
| reference only | 0% | supplied | **0.0%** | 0.0000 | 1.0× |
| reference only | 0% | fitted | **0.0%** | 0.0000 | 1.0× |
| −1 … +1 | 5% | supplied | 0.0% | 0.0000 | 1.0× |
| −1 … +1 | 5% | fitted | +1.8% | −0.0003 | 1.8× |
| −1 … +1 | 30% | supplied | 0.0% | 0.0000 | 1.0× |
| −1 … +1 | 30% | fitted | 0.0% | 0.0000 | 1.0× |
| −2 … +2 | 30% | supplied | 0.0% | 0.0000 | 1.0× |
| −2 … +2 | 30% | fitted | **−14.9%** | +0.0002 | 2.8× |
| −3 … +3 | 30% | supplied | 0.0% | 0.0000 | 1.0× |
| −3 … +3 | 30% | fitted | **+124%** | **+0.0687** | **333×** |

The first two rows are the control and read exactly zero: with every locus at the reference
length there is nothing for slippage to be confused with.

**Two readings, and the second is the one that decides a type.**

- **Every `supplied` row is exactly zero.** If the allele spectrum is known, a per-read tally
  is an unbiased estimator of all three slippage parameters, however wide the spectrum. The
  accumulator is not the problem; fitting the spectrum from the same tally is.
- **Fitted, the answer runs away as the allele spectrum widens.** Seven allele lengths against
  seven offset buckets leaves the fit more free parameters than the tally has independent
  numbers, and the level moves **333-fold depending only on where the search starts**. That is
  the signature §2.2 recorded for the pooled multi-library key — the same answer at every depth,
  moving with the start — and it means *not identified* rather than *badly estimated*. Seven
  allele lengths is not a stress case: it is one repeat unit either side of a tract at a stratum
  with ordinary variation.

**Keeping a locus's reads together removes it.** Every other measurement below keys the
accumulator by locus — one entry per locus, holding how many of that locus's reads fell at each
offset — and the same fit is then exactly unbiased in the control (0.000% on the level, 0.0000
on both shares, spread 1.000×).

### 6.3 Where the offsets are measured from

Spec §4.1 records offsets from **each locus's own modal observed length** and marks the choice
`OPEN`, with the answer to its own doubt being *"the origin is a binning choice and not a
genotype call: the fit marginalises over the genotype, so a heterozygous locus's second allele
is explained by the genotype term rather than charged to slippage."*

Three arms, all keyed by locus, differing only in where the origin sits and what the fit
believes about it. `hets` is how often the fitted genotype frequencies say a locus carries two
different alleles.

<!-- ORIGIN TABLE -->

### 6.4 Scoring a saturated end bucket

The offset range is small and its ends absorb everything beyond — "at least four repeats
short" is one bucket. Scoring that bucket as though every read in it sat exactly on the edge is
a plug-in; the marginal is the sum over everything the bucket takes in.

Loci keyed by locus, 6 reads apiece, 10 loci in 100 off the reference length.

| range | edge scoring | level | fall-off | level | fall-off |
|---|---|---:|---:|---:|---:|
| | | *9 in 100 take a second step* | | *30 in 100 do* | |
| ±1 | **marginal** | **0.00%** | **0.0000** | **0.00%** | **0.0000** |
| ±1 | plug at the edge | −6.5% | −0.0850 | −5.3% | −0.2980 |
| ±1 | plug, rescaled to sum to one | +8.4% | −0.0074 | **+32.9%** | −0.0466 |
| ±2 | **marginal** | **0.00%** | **0.0000** | **0.00%** | **0.0000** |
| ±2 | plug at the edge | +0.5% | −0.0121 | +3.6% | −0.0864 |
| ±2 | plug, rescaled to sum to one | +1.5% | +0.0086 | **+23.0%** | +0.1246 |
| ±3 | **marginal** | **0.00%** | **0.0000** | **0.00%** | **0.0000** |
| ±3 | plug at the edge | +0.04% | −0.0010 | +1.0% | −0.0241 |
| ±3 | plug, rescaled to sum to one | +0.08% | +0.0015 | +3.8% | +0.0426 |

**The marginal is exactly unbiased at every range**, including ±1, where the end buckets absorb
every read that moved more than one repeat. It costs nothing: a bucket's probability is a sum
over a handful of kernel terms.

**The plug-in's damage is a function of how much the buckets absorb, and it lands hardest on
the fall-off.** At ±1 with the measured fall-off, plugging at the edge returns a fall-off of
0.002 against a truth of 0.087 — it reports that a slipped read essentially never takes a
second step, which is the parameter's whole content. The rescaled plug-in protects the fall-off
and pays in the level instead: **+33% where 30 in 100 slipped reads take a second step**. That
is the regime long tracts sit in, so the error is largest exactly where slippage matters most.

**Neither plug-in survives the algebraic gate**, which is the cheaper way to reject it: the
un-rescaled rule sums to 0.9488 over the cell space at ±1, 0.9954 at ±2 and 0.9996 at ±3, so it
is not the likelihood of anything and no consistency result covers it.

### 6.5 The shape of the surface, and the size of the search

[`parameter_prepass.md`](../spec/parameter_prepass.md) §9.3 records that nobody has shown the
profile curve has a single hump, and §3.1 makes that the reason for stepping through the noise
parameters end to end. At one parameter the admission is cheap; at three, a flat scan is
4.2 million scores **per (read group × stratum)**.

The curve below profiles the level: at each value the other two slippage parameters *and* the
genotype frequencies are maximised out. 41 rungs from 0.0001 to 0.3.

| | interior local maxima | best rung | truth |
|---|---:|---:|---:|
| 3 reads a locus, 46 loci in 100 heterozygous | **1** | 0.0182 | 0.0200 |
| 8.5 reads a locus, 10 loci in 100 heterozygous | **1** | 0.0182 | 0.0200 |

The best rung is the ladder step nearest the truth — these rungs are 22% apart, against the
design's 6% — so the offset is resolution and not bias.

**One hump on both, which is worth exactly what it says: two worlds, one axis profiled.** It is
evidence against a second hump in the level and says nothing about the two shares. What it does
support is dropping the flat scan in favour of a search from several starts, which is what
§6.2's and §6.3's fits use throughout — coordinate ascent from four starts spread over all three
parameters — and which agrees with itself to 1.00× on every well-specified world in this
section.

**How sharp the peak is, on the other hand, differs sharply between the two.** At 3 reads with
46 loci in 100 heterozygous the profile is nearly flat over a three-fold range of the level
(0.002 nats between 0.0122 and 0.0272, either side of a peak at 0.0182); at 8.5 reads with 10 in
100 heterozygous the same span costs 0.017 nats. The level is weakly determined in exactly the
regime tomato sits in — which is a statement about precision, not bias, and the reason the
observation count beside each fit is load-bearing.

### 6.6 What merging two strata costs

Thin strata take their neighbours' value, and a fitted sequence that fails to rise with repeat
count is merged and refitted (§1.1, inherited from GATK). Both change the estimate and neither
had a bias measured against it. Two strata pooled and fitted as one, at 3 reads a locus:

| the two strata's true levels | first stratum's share of the loci | one level fitted for both | worst error carried by a stratum |
|---|---:|---:|---:|
| 2.0% and 2.0% | 50% | **2.000%** | **0%** |
| 2.0% and 2.0% | 80% | **2.000%** | **0%** |
| 2.0% and 3.0% | 50% | 2.50% | 25% |
| 2.0% and 3.0% | 80% | 2.20% | 27% |
| 2.0% and 4.0% | 50% | 2.98% | 49% |
| 2.0% and 4.0% | 80% | 2.39% | 40% |
| 2.0% and 8.0% | 50% | 4.82% | 141% |
| 2.0% and 8.0% | 80% | 3.08% | 61% |

The first two rows are the control and cost exactly nothing. **A merge returns close to the
loci-weighted mean of the two levels**, so what each stratum then carries is its own distance
from that mean: a 1.5-fold difference between neighbours costs about a quarter of the level, a
two-fold difference about half.

**What that is worth on real strata.** §5's measured levels run 0.091% below four repeats,
0.170% at four to five and 2.006% at six or more, and §4's dinucleotides reach 15.0% at 12–15
repeats — so slippage rises roughly 1.3-fold per repeat count over that range. Merging or
borrowing across **one** repeat count therefore costs on the order of 15 to 25% of the level.
That is the price of the noise it buys back, and it should be recorded per stratum rather than
inferred.

### 6.7 The substitution rate is a division, not a search

The composition channel keeps two running counts per stratum — bases compared and bases
mismatched — and each read is compared against the tract at **the length that read shows**, so
a mismatch is a substitution and not a slip.

Two consequences, both arithmetic:

- **The two counts are a sufficient statistic and the maximum-likelihood rate is mismatches
  over bases compared.** A read's mismatch count is binomial at `ε` whatever length it showed,
  so the length channel and the composition channel factorise exactly. Recovered 0.0030 from a
  truth of 0.0030 by search, which is the check that the closed form is the maximum and not
  merely a moment estimate.
- **Where a stratum holds reads of two different true rates, the pooled counters return the
  base-weighted mean** — 0.0025 from equal shares of 0.001 and 0.004, 0.0019 from 90% at 0.001
  and 10% at 0.010. That is the right answer for a model carrying one rate, and it is the same
  result §2.2 records for the generic path's shared `ε`.

**So `ε` is not an axis of the search**, and
[`parameter_prepass.md`](../spec/parameter_prepass.md) §3's "the STR path's three, scanned
together, are 4.2 million" prices the wrong three: 161³ is the right arithmetic for the wrong
parameters. What is searched is the slippage — how often, which way, how far.

### 6.8 What is not measured

**Depth above about 12 reads a locus.** The cell space is every way a locus's reads split
across the offset buckets, so it grows as the depth to the power of the bucket count; the worlds
here are exact to 12 reads. HG002 at 300× needs a coarsening of the per-locus cell that nothing
here has priced, and it is the one place the accumulator's memory could stop being kilobytes.

**Periods other than the one shaped like tomato's dinucleotides**, and ploidy above two.

**A non-whole-repeat outcome that depends on the allele.** The guard bucket is treated as an
independent per-read category, which makes it factorise out exactly; if a read is more likely to
produce a non-whole-repeat length at some allele lengths than others, it does not.

**The interaction between the origin choice and the end buckets.** Both are decisions about the
offset axis and both were priced alone. Centring on the mode shifts offsets before they
saturate, so it must change how often a read lands in an end bucket; nothing here varies the two
together.

---

## 7. Findings from the earlier review that these measurements overturned

Recorded because the same reasoning will be tempting again.

| claimed in review | what the measurement says |
|---|---|
| `F` is unidentifiable when the two states coincide, so Baum-Welch returns whatever it drifted to | the flat ridge is real but **unreachable** — every start collapses to `F` ≈ 0 (§3.8) |
| posteriors are strictly positive, so `F` has a floor and can never return exactly 0 | they underflow to exactly zero at 100 kb; `F = 0.0000` is returned routinely |
| the shuffle diagnostic would miss the outbred failure | it detects **nothing** at 100 kb, real runs included (§3.5) |
| the false-heterozygote robustness claim is untested and probably wrong | untested it was; wrong it is not (§3.3) |
| today's pooled key is 30 rungs out | the pooled key is **unbiased for `ε̄`** — the 30 rungs were the scoring rule (§2.2) |

| claimed in review | what the measurement says |
|---|---|
| a bin mean instead of a cell mean makes the objective unbounded, so every fit rails to the ladder's floor | the mechanism is real, the consequence is not: the fit lands 5.2 rungs low and 29% low on `π_hom_alt`, bounded and **silent** (§4.5) |
| the depth ladder's edges buy memory and not accuracy | 16 bins against 20, same cap, is 0.55 rungs against 0.05 (§4.3) |

And three produced by these harnesses themselves, each of which looked like a finding about
the estimator and was a defect in the measuring code:

- the multi-library sweep's **54% heterozygosity error**, which appeared identically on
  libraries with the *same* error rate and was a clamp in the measuring code;
- the **`F` collapse of §3.4**, first reported as an estimator failure when it was five
  starting points sharing one guess at the state separation;
- **`F` = 1.0000, converged, on a genome 32% covered by runs at 60 reads a site.** Two
  defects behind it, both found only because a control was run. The generator's binomial
  sampler computed `ln(1 − p)` as `(1.0 - p).ln()`, which for `p` below 2.2 × 10⁻¹⁶ rounds to
  `ln(1) = 0`, makes the gap to the next success `−∞`, saturates that to 0 on the cast to an
  unsigned integer, and so returns "every trial succeeded" where the answer is 0. Above a
  mean depth of about 40 the depth-1 cell falls under that threshold — `λ·e^(−λ)` is
  4 × 10⁻⁸ at λ = 20 and 5 × 10⁻²⁵ at λ = 60 — and every one of a window's 100,000 sites was
  assigned to "one read, no alternative": a genome with **no heterozygote anywhere**, which
  the runs model correctly read as one long run. Behind it, the fit's convergence test was
  relative to the whole genome's log-likelihood, which scales with sites × depth, so its
  tolerance was 151 nats at 3 reads a site and **2,013 nats at 60**; expectation-maximization
  stopped after four passes with the answer still moving and reported convergence. Both are
  fixed. **Every number this note published before today is at 3 reads a site**, where the
  depth-1 cell has probability 0.15 and nothing is near either cliff, and §3's tables are
  unchanged by the fixes.
