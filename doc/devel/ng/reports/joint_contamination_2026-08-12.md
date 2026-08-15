# Contamination: population structure hides it rather than inventing it, and the site budget lands on the threshold

*Research report, 2026-08-12. Covers `spec/parameter_prepass_joint_fit.md` §3.4 — the per-sample
contamination fraction, which had no measurement of any kind behind it and which sets the route's site
budget. **One program stands behind this**: `examples/ng_joint_contamination_harness.rs`. Raw output in
`tmp/contamination/`.*

---

## 1. What was asked and what came back

`α` is the fraction of a sample's reads that came from another individual — a second plant in the
tube, a neighbouring library on the same run. §3.4 fits one per sample, and §3.4.2 names the hazard
that decides whether the estimate is worth anything on a panel of landraces: **which allele frequency
the reads are scored against.**

**The test that matters is not that a contaminated sample returns its own fraction. It is that a
clean panel returns zero**, because a false positive here is a run the user is told to repeat.

**Three answers, and the first corrects the spec.**

1. **Structure breaks the estimate downward, not upward.** §3.4.2's closing sentence says a pooled
   spectrum would make structure "read as contamination". Measured, the opposite happens: at an `F_st`
   of 0.20 a sample genuinely contaminated at **3% comes back at 0.5%**, and one at 1% comes back at
   **exactly zero**. Contaminated samples pass as clean. **The paper §3.4.2 cites reports the same
   direction** — 2.9% returned for a true 10% — so the spec contradicted its own evidence in that
   sentence, and the section above it had the direction right.
2. **Scoring each sample against its own subpopulation's frequency is the fix, and it works at every
   divergence tested.** With the correct per-group frequency a clean panel gives **0 of 50 samples
   above 1%** at `F_st` 0.00 through 0.20, and a 3% contaminant comes back at 3.0% to 3.5%.
3. **But that frequency has to borrow strength across the panel, not partition it.** Estimating each
   group's frequency from its own twelve samples adds about **+0.015 to every sample's `α`** and flags
   **41 to 47 of 50 clean samples** at a 1% threshold. The danger is not which model produces the
   individual-specific frequency; it is how few samples stand behind it.

**And the budget lands exactly on the threshold.** A clean sample's worst fitted `α` falls from 1.85%
at 3,400 segregating markers to 0.86% at 13,700 and 0.53% at 55,000. The two-million-position
selection yields about 10,000 segregating markers on a fifty-sample panel, which puts the noise floor
at **about 1% — the flagging threshold itself.**

---

## 2. What was measured

`examples/ng_joint_contamination_harness.rs`. A panel of `N` samples in `G` subpopulations, drawn
Balding–Nichols: each locus has an ancestral allele frequency and each subpopulation's own frequency
is drawn around it with a spread set by `F_st`. A sample's genotype comes from its subpopulation's
frequency under a supplied inbreeding coefficient; each read comes from the sample's own genotype with
probability `1 − α` and from a contaminant's with probability `α`. The contaminant is another
individual of the same subpopulation — the second plant in the tube.

The fit is `verifyBamID`'s sequence-only two-genotype mixture (Jun et al. 2012): **both genotypes are
unknown and both are summed over**, and the search is restricted to `α ≤ ½` because that likelihood is
symmetric and cannot tell a 20% contaminated sample from an 80% one.

**Three arms, differing only in the frequency each sample's genotype is scored against.**

| arm | the frequency | what it stands for |
|---|---|---|
| `pooled` | one per locus, from the whole panel | a route that fits one spectrum |
| `by-group` | one per locus per subpopulation, from that subpopulation's members | an individual-specific frequency obtained by **partitioning** the panel |
| `true group` | the frequency the genotypes were actually drawn at | **no fit can have this** — the ceiling of the fix, and the arm that separates a frequency that is *wrong* from one that is *right and noisy* |

**The pooled frequency is estimated from the panel being scored, which includes the sample.** That is
a real difference from `verifyBamID`, whose frequencies come from an external reference panel, and it
is why the failure below is a loss of power rather than the collapse that paper reports.

---

## 3. A clean panel: what each frequency does to it

Fifty samples in four subpopulations, three reads a site, 20,000 loci, every sample's true `α` zero.
Flagged means fitted `α` at or above 1%.

| `F_st` | pooled: mean / max / flagged | by-group: mean / max / flagged | true group: mean / max / flagged |
|---:|---|---|---|
| 0.00 | 0.0040 / 0.0119 / **3 of 50** | 0.0149 / 0.0251 / **42 of 50** | 0.0016 / 0.0087 / **0 of 50** |
| 0.02 | 0.0019 / 0.0093 / 0 of 50 | 0.0147 / 0.0237 / 41 of 50 | 0.0016 / 0.0090 / 0 of 50 |
| 0.05 | 0.0003 / 0.0036 / 0 of 50 | 0.0149 / 0.0229 / 45 of 50 | 0.0011 / 0.0072 / 0 of 50 |
| 0.10 | 0.0000 / 0.0000 / 0 of 50 | 0.0166 / 0.0257 / 47 of 50 | 0.0017 / 0.0086 / 0 of 50 |
| 0.20 | 0.0000 / 0.0000 / 0 of 50 | 0.0162 / 0.0266 / 45 of 50 | 0.0014 / 0.0077 / 0 of 50 |

**Read alone, this table says the pooled frequency is the best of the three and gets better as the
panel diverges. That reading is wrong, and §4 is why**: its zeros at `F_st` 0.10 and 0.20 are an
estimator with no power, not one with good specificity. A table of a null test cannot tell those
apart, which is the reason the spike run exists.

**What the third column does say, and it is the load-bearing one:** a frequency that is *correct*
returns a clean panel clean at every divergence. So nothing about population structure makes this
parameter unmeasurable — the difficulty is entirely in getting the frequency.

**What the second column says is a warning about how.** `by-group` has the right frequency in
expectation and gets it from a twelfth of the panel; that noise alone puts 41 to 47 clean samples of
50 over a 1% threshold. **An individual-specific frequency that partitions the panel is worse than a
pooled one**, whatever it does about structure.

---

## 4. A contaminated sample: the direction of the failure

One sample of fifty genuinely contaminated, the other forty-nine clean. Fitted `α` for that sample,
and the largest fitted `α` among the forty-nine.

| `F_st` | true `α` | pooled | by-group | true group |
|---:|---:|---|---|---|
| 0.00 | 0.010 | **0.0103** / others 0.0081 | 0.0212 / 0.0206 | 0.0073 / 0.0047 |
| 0.00 | 0.030 | **0.0325** / 0.0081 | 0.0448 / 0.0206 | 0.0292 / 0.0047 |
| 0.00 | 0.100 | **0.1065** / 0.0081 | 0.1237 / 0.0206 | 0.1024 / 0.0047 |
| 0.10 | 0.010 | 0.0037 / 0.0000 | 0.0313 / 0.0248 | **0.0134** / 0.0094 |
| 0.10 | 0.030 | 0.0209 / 0.0000 | 0.0517 / 0.0248 | **0.0324** / 0.0094 |
| 0.10 | 0.100 | 0.0829 / 0.0000 | 0.1261 / 0.0248 | **0.1009** / 0.0094 |
| 0.20 | 0.010 | **0.0000** / 0.0000 | 0.0320 / 0.0272 | 0.0131 / 0.0078 |
| 0.20 | 0.030 | **0.0050** / 0.0000 | 0.0557 / 0.0272 | 0.0347 / 0.0078 |
| 0.20 | 0.100 | 0.0584 / 0.0000 | 0.1322 / 0.0272 | 0.1044 / 0.0078 |

**On an unstructured panel the pooled frequency is exactly right** — 0.0103, 0.0325 and 0.1065 for
truths of 0.010, 0.030 and 0.100. So the estimator itself works and the model is identified; there is
nothing wrong with §3.4.1.

**Structure takes its power away, in one direction.** The same estimator at `F_st` 0.20 returns
**0.0000 for a true 1%** and **0.0050 for a true 3%**. Both sit under any 1–3% flagging threshold, so
**a contaminated sample is reported as clean** — which is the failure mode that costs a study
something, since nobody re-examines a sample the pipeline passed.

**The correct per-group frequency restores it**: 0.0134, 0.0324, 0.1009 at `F_st` 0.10 and 0.0131,
0.0347, 0.1044 at 0.20 — the truth at every divergence.

**So §3.4.2's decision is right and its stated reason is not.** Individual-specific frequencies are
necessary, but not because structure would be *read as* contamination: because structure makes
contamination *invisible*. The consequence for the pipeline is the opposite of what a false-positive
story implies — the number to watch is not how many samples are flagged, it is whether a panel that
should have flagged something flagged nothing.

---

## 5. How many markers, and the budget lands on the threshold

The frequency each sample is scored against is its own subpopulation's, correct rather than estimated,
so what this prices is the estimator's own appetite and not the frequency's error. Fifty samples, four
subpopulations at `F_st` 0.10, three reads a site.

| loci | of them segregating | clean panel: mean / max `α` | one sample at 3%: fitted |
|---:|---:|---|---:|
| 5,000 | 3,434 | 0.0044 / **0.0185** | 0.0290 |
| 20,000 | 13,748 | 0.0017 / **0.0086** | 0.0300 |
| 80,000 | 54,735 | 0.0008 / 0.0053 | 0.0320 |
| 320,000 | 218,943 | 0.0003 / 0.0021 | 0.0308 |

**The contaminated sample's estimate is right at every budget** — 0.029 to 0.032 for a truth of 0.030,
from 3,400 markers up. **What more markers buy is the noise floor on the clean samples**, and that is
what a flagging threshold has to clear: 1.85% at 3,400 markers, 0.86% at 13,700, 0.21% at 219,000.

**So the budget is set by the threshold, not by the estimate.** To use a 1% threshold without false
positives the clean panel's worst `α` has to sit below 1%, which happens between 3,400 and 13,700
segregating markers.

**The two-million-position selection yields about 10,000 segregating markers on a fifty-sample panel**
([`joint_fit_estimator_2026-08-12.md`](joint_fit_estimator_2026-08-12.md) §5.6: 10,208 at 1.28 M
positions). Interpolating this table puts the noise floor there at **about 1%** — the flagging
threshold itself. **§3.4.4's arithmetic is confirmed, and it is confirmed to be tight rather than
comfortable.**

**Two courses, and the second is cheaper.** Raise the budget — 55,000 segregating markers would want
about 11 M positions, which is five times the census — or **report the noise floor beside `α` and let
the threshold be set from it**. The floor is measurable on each run without extra data: it is the
distribution of fitted `α` across the panel, and a sample is contaminated when it stands out from
that distribution rather than when it crosses a constant.

---

## 5a. Fitting the individual frequencies rather than being handed them — added 2026-08-13

§3's third column is an oracle. This section replaces it with the method: each sample gets coordinates
from a decomposition of the cohort's own genotype dosages, and at each locus the allele frequency is
fitted as a straight line in those coordinates across **all** samples. No groups, no assignment. Four
axes, and a **shrinkage** that pulls each locus's slopes towards zero by how much of that locus's
dosage spread the line actually explains — a locus whose slopes are indistinguishable from noise keeps
only its intercept, which is the pooled frequency.

### What it buys

Fifty samples in four equal subpopulations at `F_st` 0.20, three reads a site, 80,000 loci. One sample
contaminated; the largest estimate among the other forty-nine is the noise floor it has to clear.

| true `α` | pooled | 4 axes, fitted | 4 axes, shrunk | true group (oracle) |
|---:|---|---|---|---|
| 0.010 | **0.0000** / floor 0.0000 | 0.0262 / 0.0342 | **0.0146** / 0.0196 | 0.0112 / 0.0055 |
| 0.030 | **0.0022** / 0.0000 | 0.0443 / 0.0265 | **0.0320** / 0.0141 | 0.0305 / 0.0055 |
| 0.100 | 0.0570 / 0.0000 | 0.1188 / 0.0214 | **0.1018** / 0.0105 | 0.1011 / 0.0055 |

**It turns a blind estimate into a detectable one.** At a true 3% the pooled frequency returns 0.0022,
which is indistinguishable from its own clean samples; the fitted frequency returns 0.0320 against a
floor of 0.0141, which is 2.3 times it. The oracle returns 0.0305 against 0.0055, so **the fitted
version reaches the oracle's estimate and pays for it in floor**.

**Shrinkage is not optional.** Without it the estimate is 0.0443 for a true 0.030 — half again too
high — and on the unbalanced panel below it degenerates entirely.

**The floor is not a budget knob.** It falls from 0.0154 at 20,000 loci to 0.0141 at 80,000 while the
oracle's falls 0.0078 to 0.0055. So most of the floor is the cost of *fitting* the frequencies, and
more markers will not buy it back.

### The unbalanced panel, and the failure is the opposite of the one predicted

Subpopulations of 40, 5, 3 and 2 samples at `F_st` 0.20, 20,000 loci, **nobody contaminated**. The
worst spurious `α` in each group:

| group size | 40 | 5 | 3 | **2** |
|---|---:|---:|---:|---:|
| pooled | 0.0045 | 0.0000 | 0.0000 | 0.0000 |
| 4 axes, fitted | 0.0144 | 0.0389 | 0.0662 | **0.2346** |
| 4 axes, shrunk | 0.0136 | 0.0078 | 0.0133 | **0.0311** |
| true group (oracle) | 0.0069 | 0.0098 | 0.0037 | 0.0000 |

**A clean sample in the group of two is reported as 23% contaminated**, and with the contaminated
sample placed there the unshrunk fit runs to the search boundary at 0.5000 — no estimate at all.
Shrinkage brings 23% down to 3.1%, which is still three times what the group of forty gets.

**The mechanism is not the one I predicted, and the prediction was backwards.** I expected a small
group to fail to get an axis and fall back towards the panel average, which by §4 would *under*state
its contamination. The opposite happens: a small group sits at the **extreme** of an axis, which is
where a straight line is most sensitive to it, so its own noisy dosages bend the line towards
themselves. Its "expected" frequency is largely its own echo, and by §3's mechanism a noisy frequency
*manufactures* contamination.

### One number says whose estimate to trust, and it is free

How much of its own fitted frequency a sample supplies — its **leverage** — depends only on the
coordinates, not on any locus. It is therefore **one number per sample for the whole run, computable
before a single locus is fitted**, and it tracks the damage:

| group size | 40 | 5 | 3 | **2** |
|---|---:|---:|---:|---:|
| how much of its own frequency it supplies | 0.027 | 0.307 | 0.429 | **0.857** |
| spurious `α` without shrinkage | 0.0144 | 0.0389 | 0.0662 | 0.2346 |

A fair share, with four axes and fifty samples, would be 0.100. The lone pair supplies **0.857** — the
line at that sample's position is that sample — and collects 23% spurious contamination.

**So the design rule is the one ng uses everywhere else: refuse rather than guess.** Compute the
leverage, and where a sample supplies more than about half of its own frequency, emit its
contamination as *not identified* rather than as a number. It costs nothing, it needs no threshold on
the data, and it converts the failure this section measures from a silently wrong estimate into an
absent one.

---

## 6. What this cannot say

- **The drawn contaminant is from the sample's own subpopulation** — the second plant in the tube.
  A contaminant from a different subpopulation is a different and probably easier problem, and is not
  measured.
- **`by-group` is not `verifyBamID2`.** That method fits each individual's frequency as a smooth
  function of its principal-component coordinates, borrowing across the whole panel; `by-group`
  partitions. So §3's second column prices the *partitioning* failure. **§5a is the real method** and
  answers what partitioning could not.
- **The coordinates come from one pass, not from an iteration.** `PCAngsd` refines dosages and
  coordinates against each other until they settle; this fits dosages once under a pooled prior and
  decomposes them once. Whether iterating lowers §5a's noise floor is untested.
- **Four axes on four subpopulations is a generous case.** How many axes a real landrace panel needs
  is a judgement from its own scree plot, and eight axes was measurably worse here (0.0230 to 0.0369
  floor against four axes' 0.0141) — more axes mean more leverage for everybody.
- **No real reads, and no real structure.** Balding–Nichols is a model of divergence, not a tomato
  panel. What `F_st` the 63 tomato accessions actually have is unmeasured and would set which row of
  §4 they sit in.
- **The inbreeding coefficient and the error rate are supplied.** Fitting `α` jointly with them is
  what the route actually does, and the correlation between `α` and the error rate at three reads a
  site is untested here.
