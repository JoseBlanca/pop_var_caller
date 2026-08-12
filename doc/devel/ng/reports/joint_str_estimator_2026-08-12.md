# The repeat-tract half of the joint fit: the model works, the way it was to be computed does not

*Research report, 2026-08-12. Covers `spec/parameter_prepass_joint_fit.md` §4.1 and §4.2, which were
written in one sitting by analogy with the ordinary-position path and had no program behind them.
**One program stands behind this**: `examples/ng_joint_str_harness.rs`. Raw output in `tmp/str/`.*

---

## 1. What was asked and what came back

§4 of the fit spec makes two decisions about a repeat tract, and this is the first measurement of
either.

- **§4.1** — a locus's length frequencies are a latent vector drawn from a Dirichlet whose mean is
  the stratum's own length spectrum and whose concentration `κ` **says how monomorphic loci are**. A
  large `κ` is meant to return the per-stratum model exactly, so that *per locus or per stratum*
  becomes a comparison of one fitted number rather than of two designs.
- **§4.2** — the integral over that vector is done by summing over the configurations a locus can be
  in: **fixed for one length, or segregating two**.

**Three answers.**

1. **§4.1 is right, and by more than it claimed.** `κ` is identified — fitted 0.487 against a truth
   of 0.500 — and the slippage numbers come back within about a percent where the per-stratum model,
   which is what the built route uses, puts the slippage level **70.9% too low at tomato's three
   reads a site**. Per-stratum really is the large-`κ` limit: its error runs from −37.1% to −0.9% as
   `κ` goes from 0.05 to 50.
2. **§4.2 breaks exactly the property §4.1 was adopted for.** Its support cannot hold a locus
   carrying three or more lengths, and the slippage level it returns tracks how many such loci there
   are: **+0.9% where there are none, +23.7% at 18% of loci, +722% at 99.9%**. The large-`κ` limit
   §4.1 needs is precisely the regime where every locus carries three or more lengths, so under §4.2
   the comparison between per-locus and per-stratum cannot be run at all.
3. **§4.2's stated reason for existing is answered by something cheaper.** It says a Dirichlet over
   thirteen length classes "cannot be integrated by quadrature". That is true of a **grid**: nested
   quantile quadrature needs `24^(classes − 1)` points — 576 at three classes, 331,776 at five,
   1.1 × 10¹¹ at the record's nine. A **fixed low-discrepancy point set** integrates the same
   Dirichlet in **256 points at any number of classes**, returns the same answers as the grid to
   within 0.3 percentage points, and at the record's nine classes is *smaller* than §4.2's own
   support, which needs 441.

**Recommendation: keep §4.1, replace §4.2's decision with the point set.** It costs less at the
class count that matters, it is unbiased where §4.2 is not, and it is the only one of the two under
which §4.1's central claim can be tested.

---

## 2. What was measured

`examples/ng_joint_str_harness.rs` draws a stratum, draws each locus's length frequencies, draws each
sample's genotype under a supplied inbreeding coefficient, draws reads, and fits. Four descriptions
of a locus, against one truth:

| candidate | what a locus's length frequencies are | points of support |
|---|---|---:|
| `per-stratum` | the stratum's spectrum, the same at every locus | 1 |
| `dirichlet` | §4.1, over a nested quantile grid | `24^(classes−1)` |
| `dirichlet-pts` | §4.1, over a fixed Halton point set — **not in the spec** | 256 |
| `faces` | §4.2's support: one length, or two with a frequency between | `classes + 12·C(classes,2)` |

**A read** reports one of the sample's two allele lengths, or slips: the three numbers
`parameter_prepass_ssr.md` §3 fits are how **often** a read slips, which **way** — reads showing a
shorter tract outnumber longer ones 4.9 to 1 at tomato dinucleotides — and how **far**,
geometrically. A step past the recorded range saturates into the end class, as the record does.

**§4.2 is fitted as its own model rather than as an approximation of the Dirichlet**, because a face
of the simplex carries no Dirichlet density to inherit: a mass on each *fixed for one length* corner,
the rest spread over the *segregating two* edges, with the frequency along an edge on a quadrature.
That is the same four-number shape §2.1.2 fits on the ordinary-position path, with one corner per
length class instead of two.

**The inbreeding coefficient is supplied**, as the spec has it arriving from the ordinary-position
path, so nothing here says whether it could be fitted jointly with the rest.

---

## 3. The model recovers what it was given, and the per-stratum route does not

20 samples, 1,500 loci, three length classes, `κ` = 0.5 — a stratum where 30.3% of loci carry one
length across the panel and 18.1% carry three or more.

| | slippage level | which way | how far | `κ` | seconds |
|---|---:|---:|---:|---:|---:|
| truth | 0.0800 | 0.830 | 0.250 | 0.500 | |
| `per-stratum` | **−28.9%** | +5.5% | −3.9% | — | 1.6 |
| `dirichlet` | +0.0% | −0.6% | +7.2% | **−1.3%** | 91.4 |
| `dirichlet-pts` | +0.0% | −0.7% | +7.0% | **−1.6%** | 47.6 |
| `faces` | **+23.7%** | −2.8% | −4.5% | — | 8.7 |

**The fitted score is above the score at the true values in every arm**, so none of these is an
optimiser that stopped early — they are the values each description converges to.

**At tomato's actual depth the gap is far larger.** Three reads a site, 6,000 loci:

| | slippage level | which way | how far | `κ` |
|---|---:|---:|---:|---:|
| truth | 0.0800 | 0.830 | 0.250 | 0.500 |
| `per-stratum` | **−70.9%** (0.0233) | +20.5% (pinned at 1.000) | −100% (collapsed to 0) | — |
| `dirichlet-pts` | **+0.3%** | −0.2% | +1.1% | −2.7% |
| `faces` | +25.8% | −3.2% | −5.6% | — |

**At three reads a site the per-stratum model does not merely lose accuracy, it loses the
parameters**: the direction split pins at 1.000 and the fall-off collapses to zero, which is the
shape of a quantity that is not identified rather than badly estimated. **This is the strongest
single argument for the joint route on this path**, and it is on the cohort the caller is aimed at.

---

## 4. Why §4.2's support fails, and where

The support cannot represent a locus carrying three or more lengths. So the fit's only way to explain
three lengths in a locus's reads is to say the reads slipped, and the slippage level absorbs it.

**The bias tracks that population and nothing else.** Every row is 20 samples and 1,200–1,500 loci:

| loci carrying 3+ lengths in the panel | `faces` slippage level |
|---:|---:|
| 0.0% (two length classes — the support is complete) | **+0.9%** |
| 4.5% | +5.4% |
| 18.1% | +23.7% |
| 31.8% | +39.4% |
| 40.5% | +63.0% |
| 91.9% | +217.1% |
| 99.9% | **+722.1%** |

**The first row is the control that makes this a statement about the support rather than about my
parametrisation of it.** With only two length classes there is no third length for the support to
miss, and there `faces` recovers the slippage level to within a percent — as `dirichlet` does.

**The last rows are where §4.1 needs it to work.** A large `κ` is the per-stratum model, and it is
also the regime where every locus carries every length. Under §4.2, `κ → ∞` returns a slippage level
of 0.658 against a truth of 0.080. So §4.2 does not merely lose accuracy in that regime — **it makes
§4.1's one-fitted-number comparison impossible to run.**

**The spec's data-dependent candidate set does not rescue it.** §4.2 narrows the candidate lengths to
those some read reported, one step either side; that makes the support *smaller*, not wider, and the
missing population is loci needing three lengths at once.

### 4.1 It is bias, not scatter — eight times the spread across five draws

The same truth redrawn five times, 20 samples and 1,500 loci each. Fitted slippage level, against a
truth of 0.0800:

| draw | 1 | 2 | 3 | 4 | 5 | spread | distance from the truth |
|---|---:|---:|---:|---:|---:|---:|---:|
| `faces` | 0.0978 | 0.1008 | 0.0981 | 0.0988 | 0.1001 | 0.0030 | **+0.0191** |
| `per-stratum` | 0.0570 | 0.0578 | 0.0562 | 0.0578 | 0.0580 | 0.0018 | **−0.0230** |

**Neither draw ever reaches the truth**, and both errors are six to twelve times the spread between
draws, so §3's percentages are the values these descriptions converge to rather than the luck of one
data set.

### 4.2 The comparison is not rigged in the Dirichlet's favour

Everything above draws the truth from a Dirichlet, which is the family `dirichlet-pts` fits — so the
run that matters is the other one. Drawing the truth from the **faces** model instead, where 71.8% of
loci are fixed for one length and none carries three:

| | slippage level | which way | how far |
|---|---:|---:|---:|
| truth | 0.0800 | 0.830 | 0.250 |
| `faces` — its own family | −0.7% | +0.0% | +1.2% |
| `dirichlet-pts` — **misspecified here** | **−1.0%** | +0.1% | +0.9% |
| `per-stratum` | −38.0% | +12.4% | −19.9% |

**The Dirichlet fit recovers the slippage numbers whether or not the truth is a Dirichlet; the faces
support recovers them only when the truth is a faces mixture.** Its fitted `κ` comes out 0.109 where
the Dirichlet truth had 0.500, which is the concentration doing its job — a truth that is 71.8%
monomorphic is a more monomorphic stratum, and `κ` says so. **`κ` is a description of the stratum,
not a quantity to check against a number from another family.**

---

## 5. The concentration is identified across a thousand-fold range, and per-stratum is its limit

Three length classes, 20 samples, six reads a locus, 1,200 loci.

| `κ` | loci at one length | loci at 3+ lengths | `per-stratum` | `dirichlet-pts` | `faces` |
|---:|---:|---:|---:|---:|---:|
| 0.05 | 87.1% | 0.4% | −37.1% | **−0.2%** | −0.0% |
| 0.20 | 59.5% | 4.5% | −32.8% | **+0.6%** | +5.4% |
| 1.00 | 11.7% | 40.5% | −22.5% | **+1.0%** | +63.0% |
| 5.00 | 0.1% | 91.9% | −9.9% | **−1.0%** | +217.1% |
| 50.0 | 0.0% | 99.9% | −0.9% | **+0.1%** | +722.1% |

Two things are in this table.

- **§4.1's claim is confirmed**: the per-stratum model's error runs from −37.1% to −0.9% as loci stop
  being monomorphic, so it really is the large-`κ` limit and the choice between the two really is one
  fitted number.
- **One estimator covers both regimes.** The Dirichlet fit is within one percent of the truth at
  every `κ` over a thousand-fold range, so nothing has to be decided in advance about how
  monomorphic a stratum is.

**Which regime real strata are in is not measured and cannot be until the STR records exist.** If
tomato's strata sit near `κ` = 0.05 then §4.2's support would have cost almost nothing on that data —
but it would still have made the comparison in §4.1 unrunnable, which is the reason to replace it
rather than to keep it and hope.

---

## 6. What the point set costs, and why the class count is the whole argument

The record holds offsets over ±4, which is **nine length classes**.

| classes | `dirichlet` grid | `faces` support | `dirichlet-pts` |
|---:|---:|---:|---:|
| 2 | 24 | 14 | 256 |
| 3 | 576 | 39 | 256 |
| 5 | 331,776 | 125 | 256 |
| 9 | 1.1 × 10¹¹ | 441 | **256** |

**At the record's own width the point set is the smallest of the three**, and it is the only one whose
cost does not grow with the class count at all. Measured at three classes it returns the grid's
answers to within 0.3 percentage points on `κ` and 0.2 on every slippage number, at half the time;
at five classes the grid cannot be run and the point set fits in 96 seconds.

**The points are fixed and the quantile map is continuous in `κ`**, so the objective is a smooth
function of the concentration rather than a jittery one — which is what stops the search chasing
sampling noise, and is why this is quadrature rather than Monte Carlo.

---

## 7. What this cannot say

- **The truth is drawn, not real.** §4.1 separates the bias from the scatter — the errors are six to
  twelve times the spread between draws — but nothing here is real reads.
- **No real reads.** Which `κ` a tomato stratum has, and how many of its loci carry three lengths, is
  the measurement that decides how much any of this matters in practice, and it needs the STR records
  (`parameter_prepass_joint_records.md` §3), which are typed but not yet filled by a walk.
- **The inbreeding coefficient is supplied.** Whether it can be fitted jointly on this path is
  untested.
- **One stratum at a time.** Borrowing across thin strata and the monotonicity constraint along the
  repeat-count axis (`parameter_prepass_ssr.md` §4.3) are not exercised.
- **Biallelic reads only in the sense that the record is.** The substitution channel — the difference
  list that carries interruptions — is not in this model at all; it is fitted by a division and does
  not interact with the slippage numbers.
