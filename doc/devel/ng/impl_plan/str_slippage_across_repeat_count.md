# ng — how slippage varies with repeat count: the research, then the fix

**Status:** plan, 2026-08-19. Nothing built. **Written to be run in its own conversation**, so it
states its own context rather than assuming the reader was in the one that raised it.

This turns a question into build order. The question came out of
[`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §4.4 and it belongs to the parameter
fit, not to the caller: **the four slippage numbers are fitted separately for every
(motif period, repeat count) cell, and nothing makes neighbouring cells agree with each other.** The
owner's instruction is to find out what shape the numbers actually follow, decide how to fit that
shape, and then change the fitting module.

**Half the answer is already measured and §5 gives it.** What this plan adds is the other half: which
curve, fitted how, and per what.

---

## 1. What this decides

**Three decisions, in order, and the second cannot be taken before the first.**

1. **What shape does each slippage parameter follow as repeat count rises?** A straight line, an
   exponential, something that flattens at the top, or no assumed shape at all.
2. **How is that shape fitted?** What each cell is weighted by, whether the curve is fitted to the
   per-cell answers or straight to the reads, and what happens beyond the repeat counts the data
   reaches.
3. **Fitted per what?** The owner's guess is *per motif period, independently*. §6 argues that is
   the right default and the right thing to test first, and names the two ways it could be wrong.

**The output is a change to the fitting module** — `src/ng/parameter_estimation/joint/ssr_fit.rs`
and the spec that governs it, [`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md)
§4.3 — plus a report carrying the numbers that justified it. **The owner decides on the numbers**;
no threshold is fixed here.

---

## 2. Why this is worth doing

**Slippage rises steeply and smoothly with repeat count, and the fit currently treats each cell as
though the others did not exist.** On tomato, a read misreads a homopolymer of 8 repeats 2 times in
1,000 and one of 12 repeats 9.4 times in 1,000 — **4.7 times as often over four repeat counts**
([`../reports/str_fit_on_real_records_2026-08-13.md`](../reports/str_fit_on_real_records_2026-08-13.md)
§6.1). Nothing in the estimator knows that.

**Two devices stand in for the knowledge, and both are blunt.** A cell too thin to fit **borrows** a
neighbour's value outright, which costs 15 to 25% of the level per repeat count borrowed across. And
where the fitted sequence dips, two cells are **merged and refitted together**, which costs each of
them its distance from the pooled mean — about a quarter of the level at a 1.5-fold difference, half
at 2-fold, and up to 141% at 4-fold ([`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md)
§4.3). **How often the merge fires when the truth is monotone has never been measured**, and that
document records it as unmeasured.

**A curve replaces both.** It borrows from every cell rather than from one neighbour, weighted by
distance and by evidence; and it enforces monotonicity by its own shape, so the merge rule has
nothing left to do and its unmeasured failure mode goes with it.

**The thin end is where this pays.** Of the 71 cells the tomato bench run produced, **6 clear the
refusal floor of 50 tracts, and those 6 hold 88% of the tracts**. Everything else is currently either
borrowed or refused.

---

## 3. Scope

**In.**

- The four numbers fitted per (period, repeat count): the slippage level, the direction split, the
  fall-off, and the substitution rate.
- The shape each follows across repeat count, the fitting method, and the grouping.
- Replacing the borrow rule and the merge rule where the curve subsumes them.
- What the emitted parameters must carry so that a smoothed value is still distinguishable from a
  measured one.

**Out.**

- **The read likelihood.** It looks the numbers up by the candidate's own period and repeat count
  ([`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §4.4) and does not care whether the
  table behind that lookup came from a cell or a curve. Swapping a lookup for a function evaluation
  costs it nothing; §11 states the one interface consequence.
- **The generic path's error rate.** Different object, no repeat-count axis.
- **Per-locus adaptation of the stutter parameters inside the caller's own loop.** That is a separate
  question with its own home ([`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §6.1), and
  it is about a *locus* departing from its cell, not about cells relating to each other.
- **Changing the stratification axes.** §6.3 raises one candidate and hands it on rather than
  settling it here.

---

## 4. Principles — how the order below was chosen

1. **Fit the curve against cells that were fitted independently, or the exercise is circular.** A
   cell whose value was borrowed from its neighbour will of course lie on a smooth curve through its
   neighbours. **Every fit this plan learns from must run with borrowing and merging switched off**,
   which the 2026-08-13 run already demonstrates is possible (`SSR_BORROWING_FLOOR=0`).
2. **Read the data before the literature, and the literature before writing code.** The shape is an
   empirical question about these libraries; published work says what shapes are plausible and what
   others found, which is a prior and not an answer.
3. **The test of a smoothed value is a cell it never saw.** Curve-fitting quality on the cells used
   to fit it can always be improved by adding parameters. Held-out prediction cannot.
4. **A curve that fits beautifully and changes no genotype has not earned a change.** The last step
   is calls, not curves.
5. **Smoothing must not launder thinness.** After this change a value can come from a cell with two
   tracts behind it and look identical in the output to one with two thousand. What the provenance
   has to carry is §10.

---

## 5. What is already known, so that nobody re-measures it

**The tomato bench run, 63 accessions over 80 spans and 8 Mb, borrowing off**
([`../reports/str_fit_on_real_records_2026-08-13.md`](../reports/str_fit_on_real_records_2026-08-13.md)
§6.1). Five consecutive homopolymer cells, each fitted from its own tracts alone with nothing linking
them:

| repeats | tracts | reads crossing | level | shorter-share | fall-off |
|---:|---:|---:|---:|---:|---:|
| 8 | 2,082 | 937,557 | 0.0020 | 0.595 | 0.636 |
| 9 | 887 | 363,500 | 0.0027 | 0.637 | 0.631 |
| 10 | 350 | 121,244 | 0.0037 | 0.623 | 0.629 |
| 11 | 153 | 34,904 | 0.0059 | 0.713 | 0.540 |
| 12 | 83 | 15,976 | 0.0094 | 0.754 | 0.758 |

**Fitting a curve through those five, unweighted** — the first pass this plan is meant to replace
with something careful, reported here so the second pass has a number to beat:

| | span over 8→12 repeats | straight line in repeat count | in the **log** |
|---|---|---|---|
| **level** | 4.7-fold | worst error **43%**, and worst at 8 repeats | **R² 0.99**, 1.474× per repeat, worst error 10% |
| **shorter-share** | 1.27-fold | R² 0.88, worst error 0.041, dips at 10 | — |
| **fall-off** | 1.40-fold | **R² 0.10** — no trend, only noise | — |

**Three things this already suggests, each of which the research must confirm or overturn rather
than inherit:**

- **The level is exponential in repeat count, not linear.** Production fits a straight line
  (`baseline + slope · repeat_count`, [`em.rs:385`](../../../../src/ssr/cohort/em.rs)),
  and over a 4.7-fold rise a line misses by 43% — worst at 8 repeats, the bottom of the range, which
  is exactly where the copy-floor decision reads the level
  ([`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §5.1).
- **The three parameters do not behave alike.** Smoothing the fall-off on this evidence would
  manufacture a trend that is not there.
- **The five cells were monotone with no constraint applied**, so on this cohort the merge rule would
  never have fired. That is one data point on a question the spec calls unmeasured.

**Other numbers that bear on the design and should not be re-derived.** The direction split differs
by period far more than by repeat count: tomato homopolymers 1.4× and dinucleotides 4.9×; human
homopolymers 1.9× and dinucleotides 3.4× ([`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md)
§3). Precision differs by parameter: holding the direction split to 6% of itself takes about 1,400
slipped reads and the fall-off about 4,000, where the level reaches 5% of itself on 455 slipped reads
out of half a million (§4.5 there). And the whole run took 4,539 s of which the repeat-tract fit was
1,690 s — so a full refit is about half an hour, which is what §8's step budget is built around.

---

## 6. Fitted per what — the owner's guess, and where it could be wrong

**The proposal is one curve per motif period, fitted independently. That is the right default.** The
periods are not the same physics: the direction split runs 1.4× at tomato homopolymers against 4.9×
at dinucleotides, and the share of slippage that moves a whole repeat runs 98% at dinucleotides, 95%
at trinucleotides and 93% at hexamers. **One curve across periods would flatten a threefold
difference that is measured.**

**Two ways it could still be wrong, and the plan tests both.**

### 6.1 The thin periods are the ones that most need help, and independence denies it to them

Six periods, each running from its copy floor up to what a read can span — 338 cells by arithmetic,
of which the ones holding tracts are far fewer. Homopolymers and dinucleotides hold nearly all the
data; periods 5 and 6 may not have enough cells to fit a curve *at all*, and those are exactly the
cells that borrowing serves worst today.

**So test partial pooling as well as independence:** one shared slope across periods with a
per-period offset, against a free slope per period. **This project has a measured warning against
reflexive partitioning** and it should be read before the arm is designed: fitting allele frequencies
per subpopulation from twelve samples each did *worse* than pooling them, and the record is explicit
that this was "a defect of partitioning, not of sample size", so a larger cohort is not what makes
partitioning safe ([`../spec/calling_priors.md`](../spec/calling_priors.md) §10).

### 6.2 Repeat count may not be the only axis the curve should take

Candidates to test as a second variable, in the order they are cheap:

- **tract length in base pairs** (period × repeat count) against repeat count alone — these coincide
  at period 1 and diverge fast, and the spec's §4 chose repeat count on evidence that should be
  re-read rather than re-assumed;
- **purity** — production carries a per-locus purity factor on the level and ng's fit has no fitted
  source for one ([`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §4.6);
- **read group**, which is already the outer grain, but whose *slope* might or might not be shared
  across libraries.

### 6.3 A question this plan raises and does not settle: period, or motif?

**The stratification axis is the motif's length, not the motif.** A poly-A tract and a poly-G tract
of the same length are one cell. Sequencing chemistry is known to treat them differently, and if that
difference is large it is currently being averaged into a single homopolymer curve. **Measurable with
what this plan already builds** — split period 1 by base and fit both — but changing the axis is the
stratification spec's decision, not this plan's. **Report the size and hand it on.**

---

## 7. Instrument and discipline

- **Every fit that feeds a curve runs with borrowing and merging off.** *Since 2026-08-20 there is
  nothing to switch: the joint route has neither, and `SSR_BORROWING_FLOOR` is gone with them
  ([`../spec/str_slippage_level_curve.md`](../spec/str_slippage_level_curve.md) §5.1). The principle
  stands and is now structural rather than a run setting.* Formerly `SSR_BORROWING_FLOOR=0`, and
  the merge rule disabled the same way. Any run that does not is unusable for §4's first principle,
  and the report must say which setting produced each table.
- **Two cohorts, not one.** Tomato (63 accessions, ~3 reads a position, a selfing crop) and
  HG002 (one sample, high coverage, human). A shape that holds on one and not the other is a fact
  about that dataset; the specs say so about every other number and it applies here.
- **Every fitted cell carries its evidence into the analysis** — tracts, reads crossing, and slipped
  reads — because §8's weighting question cannot be answered without them and because §5's
  unweighted first pass is exactly what needs improving on.
- **Builds and runs go through `./scripts/dev.sh`** where a container runtime exists, and through
  `cargo` directly on machines that have none, per `CLAUDE.md`. Scratch output goes under the
  repository's own `tmp/`, never the system one.

---

## 8. The steps

### A — Establish the target

**A1. Refit both cohorts with borrowing and merging off, and dump every cell.** One table per cohort:
period, repeat count, tracts, reads crossing, slipped reads, and the four fitted numbers. About half
an hour per cohort on the tomato bench region set. **Produces:** the data every later step reads.
**Gate:** at least three periods with four or more consecutive populated cells, or the shape question
cannot be asked at that period and the report must say so rather than fitting two points.

**A2. Measure what the current rules do to that table.** Refit with borrowing and merging on, and
report per cell how far the value moved and which rule moved it. **This is the baseline any curve has
to beat**, and it also answers the spec's unmeasured question: how often the merge rule fires when
the unconstrained sequence was already monotone.

### B — Read

**B1. The published work on how slippage varies with repeat number.** The microsatellite
mutation-rate literature has measured this directly for decades and consistently reports a steep
rise with repeat number; **what the research must establish is the functional form the field fits and
whether it saturates**, not merely that the rate rises. Start from the direct human estimates and the
review literature, follow citations, and **record what each source actually fitted** — a rate per
generation is not a per-read slippage rate, and the two must not be conflated in the write-up.

**B2. What the neighbouring software does, from source we hold.** GATK's DRAGstr estimator is
vendored and fits per (period, repeat count) with its own borrowing and monotonicity devices — read
how it relates neighbouring cells, since ng's two rules were copied from it. HipSTR fits per locus by
expectation-maximization instead. GangSTR uses a per-locus model. **Licence discipline applies:
read to understand, implement from the publication** — HipSTR is GPL-2, GangSTR GPL-3-or-later, and
this project is neither ([`../spec/locus_generation_ssr.md`](../spec/locus_generation_ssr.md) §3).

**B3. One paragraph per source, in the report, saying what it fitted and what that implies here.**
A reading step with no written output is a step nobody can check.

### C — Choose the shape

**C1. Fit the candidate families to A1's cells, per period, weighted (see C2).** At least:

| candidate | why it is in the list |
|---|---|
| straight line in the level | production's current shape — the thing to beat |
| **log-linear** (exponential) | §5's first pass, R² 0.99 on five cells |
| log-linear with a ceiling | slippage cannot exceed 1, and long tracts may saturate; the exponential must eventually be wrong |
| power law in repeat count | the other standard two-parameter monotone family |
| **isotonic regression** | monotone by construction with **no shape assumed** — the principled replacement for the merge rule, and the control that says how much the parametric families are assuming |

**C2. Settle the weight before comparing the families, because it changes the ranking.** The cells
run from 2,082 tracts down to 83, so an unweighted fit lets the thinnest cells pull hardest — which
is what §5's first pass did. Candidates: tracts, reads crossing, and **slipped reads**, the last
being what actually sets the level's precision (§5). **State the choice and show one family fitted
both ways**, so the reader can see whether the ranking depended on it.

**C3. Choose on held-out cells, not on fit.** Leave out each populated cell in turn, fit the rest,
predict it, and report the error distribution per family. **This is the criterion that decides**, and
it is the one adding parameters cannot game.

**C4. Ask the same question of all four parameters separately, and expect different answers.** §5
already shows the level trending strongly, the direction split weakly, and the fall-off not at all.
**A parameter with no trend keeps its per-cell value and its borrowing rule** — say so explicitly
rather than smoothing everything because the machinery is there.

### D — Choose the fitting method

**D1. Two-stage against one-stage.** Two-stage fits each cell first and then a curve through the
point estimates — simple, and it treats noisy estimates as though they were data. One-stage fits the
curve's parameters directly against the reads, which is statistically the right thing and more work.
**Measure the gap on the held-out criterion** before deciding; if it is small, the simpler one wins
and the report says by how much.

**D2. Independence against partial pooling across periods** (§6.1), on the same criterion, with
periods 5 and 6 as the cases that decide it.

**D3. Beyond the observed range.** Long tracts a read can barely span have few or no cells. Say what
the curve does there — extrapolate, hold flat at the last fitted value, or refuse — and **report
which**, because a caller reading an extrapolated value cannot currently tell.

**D4. Where a curve and a well-measured cell disagree.** A cell with 2,000 tracts behind it is
better evidence than the curve through its neighbours. Decide whether the emitted value is the
curve, the cell, or a weighted blend, and on what.

### E — Implement

**E1. The curve, behind the existing fit's seam**, so A2's baseline can be re-run against it at any
time. Both must be selectable in one build.

**E2. Retire the two rules the curve subsumes**, and only those. If C4 leaves a parameter unsmoothed,
its borrowing rule stays and the spec says which parameters are on which route.

**E3. Provenance** (§10).

**E4. Amend [`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3** — the
borrowing and merging section is the one this replaces — and its architecture sibling.

### F — Validate

**F1. Held-out error against A2's baseline**, per cohort, per period, per parameter.

**F2. Genotype concordance.** The generator that chose the read-likelihood model builds a cohort with
known genotypes ([`sim.rs`](../../../../src/ssr/cohort/sim.rs)); score calls under the old parameters
and the new, at the committed depth range. **A curve that changes no call has not earned the change**
— report that outcome as readily as the other one.

**F3. HG002's tandem-repeat benchmark**, as a regression guard on real data. Its limits are on record
and should be quoted rather than discovered: it reaches the well-covered loci, where any parameter
set does well.

---

## 9. Verification summary

| step | what must be true before the next one starts |
|---|---|
| A1 | three periods with four or more consecutive populated cells, or the shape question is not asked at that period |
| A2 | every cell's movement attributed to borrowing or to merging, and the merge-firing rate reported |
| B3 | one written paragraph per source |
| C2 | the weight chosen, with one family shown fitted both ways |
| C3 | held-out error per family, per period, per parameter |
| C4 | an explicit keep-or-smooth verdict for each of the four parameters |
| D1–D4 | each answered with a number, not a preference |
| E2 | no rule retired that the curve does not subsume |
| F2 | calls compared, and the result reported whichever way it falls |

---

## 10. Provenance — what must not be lost

**After this change a value fitted from 2,000 tracts and one interpolated across a gap look the
same.** Today's `Provenance::Borrowed` marks the second kind and this plan removes the mechanism that
sets it, so it must be replaced rather than dropped. What the emitted parameters have to carry, per
cell:

- **how much evidence stood behind that cell itself** — tracts, reads crossing, slipped reads — which
  the library already records and the specification does not yet require of the emitted table;
- **whether the emitted value is the cell's own fit, the curve, or a blend**, and if a blend, its
  weight;
- **whether the cell was inside the fitted range or beyond it** (D3);
- **the curve's own goodness of fit** for that period, so a consumer can tell a well-determined curve
  from one fitted through three points.

---

## 11. The one consequence outside the fitting module

**The read likelihood looks these numbers up per candidate, not per locus**
([`../spec/read_likelihoods.md`](../spec/read_likelihoods.md) §4.4): a candidate of 6 repeats and one
of 12 at the same tract are drawn from different cells and slip at measurably different rates. If the
emitted object becomes a curve rather than a table, that lookup becomes an evaluation. **It costs
nothing and it is not optional to state**: the model must not be given a table that has silently been
resampled from a curve at cell centres, because then a candidate whose repeat count falls between
cells gets the wrong end of a rounding nobody documented.

---

## 12. What could make this exercise worthless

- **Fitting the curve to cells that were themselves borrowed or merged.** Circular, and it would look
  like a triumph. §4's first principle exists for this.
- **Judging on curve fit rather than on held-out cells.** Any family can be made to fit the cells it
  was fitted to.
- **One cohort and one period.** The level's exponential rise is currently five homopolymer cells on
  one crop. Two cohorts and at least three periods, or the finding is about tomato.
- **Smoothing all four parameters because the machinery exists.** The fall-off shows no trend on the
  evidence we have; fitting a line through it would invent one and hide the fact that it is borrowed.
- **Losing the record of thinness** (§10). The failure would be silent and would surface much later
  as a caller confidently using a number nobody measured.
- **Stopping at the curve.** The reason to do any of this is genotypes. F2 is the step that says
  whether it mattered.
