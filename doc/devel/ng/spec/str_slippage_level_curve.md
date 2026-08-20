# ng — the slippage level as a curve across repeat count

**Status:** design, 2026-08-20. **No code yet — this settles the design.** It amends the STR half
of the parameter pre-pass, [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3 and §4.5,
and is the first piece of the research plan
[`../impl_plan/str_slippage_across_repeat_count.md`](../impl_plan/str_slippage_across_repeat_count.md)
to reach a settled design. Build order lives in
[`../impl_plan/str_slippage_level_curve.md`](../impl_plan/str_slippage_level_curve.md). The
measurements it rests on are in
[`../reports/str_slippage_shape_2026-08-20.md`](../reports/str_slippage_shape_2026-08-20.md).

---

## 1. What this is

**Today the slippage level is fitted separately in every (motif period, repeat count) cell, and
nothing makes neighbouring cells agree.** A cell too thin to fit takes a neighbour's value whole,
which costs 15 to 25% of the level per repeat count borrowed across
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3). On the tomato bench run 65 of 71
cells got no answer from their own tracts at all.

**This replaces that, for one of the four numbers.** The slippage level becomes a curve in repeat
count, fitted once per motif period and evaluated at every cell — so a thin cell reads a line
drawn through all its neighbours weighted by their evidence, instead of copying whichever
neighbour the borrowing rule reached first.

### 1.1 Goals

1. **Every populated cell gets a level**, including the ones below today's refusal floor of 50
   tracts ([`ssr_fit.rs:438`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L438)),
   without any cell's answer being one neighbour's answer.
2. **The shape of the rise is fitted, not chosen.** The two cohorts we have prefer opposite
   shapes at the same repeat counts (§2), so a design that hard-codes either one is wrong on
   the other.
3. **How far a cell is pulled toward the curve is set by two measured precisions, not by a
   category.** There is no "well-measured" class of cell: the pull is the ratio of the curve's
   error to the cell's own, and both are numbers the fit already has (§7).
4. **A smoothed value stays distinguishable from a measured one** in what the pre-pass emits.

### 1.2 Non-goals, and what this does not do

- **It does not touch the other three numbers.** The direction split, the fall-off and the
  substitution rate keep their per-cell fit and today's borrowing rule. Whether they should be
  smoothed too is a separate question, and the measurements point both ways — the fall-off is
  flat across tomato's five cells and spans 3.14-fold across HG002's twenty-three
  ([`../reports/str_slippage_shape_2026-08-20.md`](../reports/str_slippage_shape_2026-08-20.md)
  §5). Its home is the research plan's step C4.
- **It does not fit the curve to the reads.** The curve is fitted to the per-cell answers — the
  two-stage route. Fitting the curve's own parameters directly against the reads is the research
  plan's D1 and is deferred with a home (§10).
- **It does not extrapolate.** Outside the repeat counts the curve was fitted over, the level is
  held at the nearest fitted end rather than continued (§6).
- **It does not change how the read likelihood looks a level up.** The consumer indexes by the
  candidate's own period and repeat count
  ([`read_likelihoods.md`](read_likelihoods.md) §4.4) and receives a table, exactly as now; what
  changes is where the numbers in that table came from.
- **It does not change the stratification axes.** Whether period should be motif is raised in the
  research plan §6.3 and stays there.
- **It does not adapt anything per locus.** That is the caller's own loop
  ([`read_likelihoods.md`](read_likelihoods.md) §6.1).

---

## 2. What the data says, and why it fixes the family

Two cohorts, both fitted with borrowing off so every cell speaks from its own tracts alone
(`SSR_BORROWING_FLOOR=0`). Numbers and method in
[`../reports/str_slippage_shape_2026-08-20.md`](../reports/str_slippage_shape_2026-08-20.md);
what matters here is the three facts that decide the design.

**First: the rise is steep at short tracts and flattens at long ones.** On HG002 — one sample at
about 300 reads a position, 23 consecutive homopolymer cells from 8 to 30 repeats — the level
runs from 3.7 reads slipping per 1,000 to about 120, a 37-fold rise. But the step-to-step ratio
falls from 2.45-fold at 9→10 repeats to 1.02-fold at 19→20. Predicting each cell from the other
22, a straight line in repeat count lands within **4.4%** of the held-out cell and an exponential
within **22.7%**.

**Second: the two cohorts prefer opposite shapes over the same repeat counts.** Restricted to
8–12 repeats, five cells each:

| | level, 8→12 repeats (reads slipping per 1,000) | held-out error, exponential | held-out error, straight line |
|---|---|---:|---:|
| tomato, 63 accessions at ~3 reads a position | 2.1, 2.9, 3.9, 6.6, 10.2 | **12.4%** | 33.6% |
| HG002, one sample at ~300 reads a position | 3.7, 6.7, 16.5, 22.0, 28.0 | 31.2% | **8.0%** |

So it is not a matter of one cohort seeing a wider window. **Neither shape can be hard-coded.**

**Third: both fixed shapes fail catastrophically outside the range they saw.** An exponential
fitted on HG002's 8–12 cells says the level at 30 repeats is **205** — as a probability — where
the cell at 30 repeats fits 0.120. A straight line fitted over 8–30 goes **negative below 7.4
repeats**, and the homopolymer copy floor puts real cells at 8
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5.1.1).

**One family covers both, with one number saying which.** Write `n` for the repeat count and
`level` for the slippage level:

```text
level ^ rise_shape  =  intercept + slope · n            for rise_shape > 0
log(level)          =  intercept + slope · n            for rise_shape = 0
```

**`rise_shape` says how the level compounds:** at **0**, each extra repeat *multiplies* the level
by a fixed factor; at **1**, each extra repeat *adds* a fixed amount. Fitted on a grid over
`[0, 1]`, it lands at 0.00 on tomato's homopolymers, 1.00 on HG002's, and 0.80 on HG002's
dinucleotides — where it beats both ends, 3.8% held-out against the straight line's 5.9% and the
exponential's 18.3%.

**`rise_shape = 1` is production's own shape.** Production fits `baseline + slope · units` and
clamps to `[0, 1]` ([`em.rs:385`](../../../../src/ssr/cohort/em.rs#L385)), so this family nests
what production does rather than replacing it with something unrelated.

### 2.1 Two families that were live and lost

- **A saturating exponential**, `ceiling · (1 − exp(−rate · (n − offset)))` — the shape the
  microsatellite mutation-rate literature settles on after finding the exponential too simple. It
  ties the straight line on HG002's homopolymers (4.3% against 4.4%) and loses on tomato's five
  cells (19.9% against the exponential's 12.4%), because three free numbers cannot be determined
  from five cells spanning four repeat counts. **Rejected for needing more cells than the thin
  periods have**, not for fitting worse where the cells exist.
- **GATK DRAGstr's per-base hazard**, `1 − (1 − q)^(period · n)` — the shape the neighbouring
  software assumes inside each cell
  ([`DragstrParametersEstimator.java`](../../../../gatk/src/main/java/org/broadinstitute/hellbender/tools/dragstr/DragstrParametersEstimator.java)
  builds `log10PCorrectByLength` as `lengthInBases × log10PCorrectPerPosition`). It predicts a
  held-out cell to 28.6% on HG002's homopolymers against the straight line's 4.4%. **Kept in the
  measurement as the control that says how much the free families are assuming**, and rejected as
  a family. This is also why DRAGstr must refit `q` freely in every cell.

---

## 3. The grain: what is fitted per read group, and what is shared

**Three numbers describe one period's curve, and they are not fitted at the same grain.**

| number | grain | why |
|---|---|---|
| `intercept`, `slope` | **per (slippage group, motif period)** | This is where a library's own chemistry lives: the two cohorts differ 1.8-fold in the level at 8 repeats and 2.7-fold at 12. Both are ordinary two-parameter fits and a period with four cells can carry them. |
| `rise_shape` | **per motif period, shared by every slippage group** | It is a curvature, and a curvature needs the whole span of repeat counts visible at once. On tomato one library at 3 reads a position puts about 12 slipped reads behind a cell — enough for a level, not for a shape. |

**A slippage group is the set of read groups a run fits one set of slippage numbers for**
([`ssr_fit.rs:291`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L291) — `slippage`
is one entry per group, `None` where that group put no read in the cell). The walk builds it
([`ng_joint_records_walk.rs:747`](../../../../examples/ng_joint_records_walk.rs#L747)): one group
per read group, or all of them pooled.

**⚠ What "per read group" can mean today is narrower than it sounds, and a coder must know
this before promising anything.** A read group is identified by its index *within one sample*
([`types.rs:210`](../../../../src/ng/types.rs#L210) — `ReadGroupId(pub u32)`), and the cohort's
list is the *union* of those indices across samples
([`census.rs:1670`](../../../../src/ng/parameter_estimation/joint/census.rs#L1670)). Every tomato
CRAM declares its own library — `LB:PRJNA454805_SRR7279481` and so on, 63 of them — and every one
of them is read group 0 inside its own sample, so the cohort sees **one** read group and the run
prints `1 read groups in 1 slippage group, pooled`. **So this design fits per slippage group as
the census defines it, and on a cohort of single-library samples that is one curve for the whole
cohort.** Making the read-group axis carry library identity across a cohort is a change to the
census with its own home (§10) — it is not smuggled in here, and no claim in this document
depends on it.

---

## 4. How the curve is fitted

**Two-stage: fit the cells, then fit a line through them.** Stage one is today's fit with
borrowing switched off, so every cell's level is its own tracts' answer and nothing has been
smoothed before the smoothing (the research plan §4's first principle). Stage two draws the
curve.

### 4.1 Which cells feed a curve

A cell feeds its period's curve when **its own fit returned a level** — that is, it cleared the
refusal floor on its own tracts. Cells below the floor consume the curve and do not contribute to
it: a cell with 8 tracts has nothing to say about the shape and would pull hardest if it were let
in unweighted.

**A period needs at least four contributing cells before a curve is drawn at all.** Below that,
the period keeps today's behaviour — per-cell fits with the borrowing rule — and says so in the
provenance. *Four is the smallest count at which leaving one cell out still leaves a line and a
spare; it is a floor on arithmetic, not a measured threshold, and it is soft.*

### 4.2 The weight

**Each cell is weighted by its slipped reads** — `level × spanning reads`, the count that actually
sets how precisely a cell determines its level
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.5).

**The choice barely matters and the measurement says so, which is why it is recorded as soft.**
On HG002's 23 homopolymer cells the winning family's held-out error moves from 5.13% unweighted to
4.77% by tracts, 4.65% by reads crossing and 4.39% by slipped reads, and the ranking of families
is identical under all four. Slipped reads is chosen because it is the one with a reason behind
it, not because it won.

**⚠ Do not use `StratumEvidence::reads_off_reference_length()` as the slipped-read count.** It
counts reads whose tract differs from the *reference* length, which at a polymorphic tract is
mostly genuine allele length: on HG002 at 30 repeats 60 reads in 100 sit off the reference length
and the fit attributes only 12 of them to slippage.

### 4.3 Choosing `rise_shape`

Over a grid of 21 rungs from 0.00 to 1.00 in steps of 0.05, for each rung: fit `intercept` and
`slope` for every slippage group at that period, then score the rung by leaving each contributing
cell out in turn, refitting that group's line without it, and predicting it. **The rung with the
lowest median relative error over all left-out cells of all groups wins.**

**Ties go to the larger `rise_shape`**, because the failure modes are not symmetric: a level that
is too small under-weights real slippage, where the exponential end can return a number above 1
and has to be clamped.

**The score continues the line where a deployed curve holds it flat (§6), and that difference is
deliberate — DECIDED 2026-08-20 (owner).** Leaving out the lowest or highest cell puts it outside
what remains, so scoring it the way a deployed curve behaves would predict it with its
neighbour's value, which says nothing about the shape — and the end cells are exactly where a
shape shows. Measured three ways over both cohorts, the chosen rung is the same at four of the
five periods that can be scored; holding flat instead moves HG002's dinucleotides from 0.80 to
0.70, and scoring only the cells that stay inside the range leaves 3 of tomato's 5 cells scored
and moves its answer from 0.00 to 0.15. **This does not reach a caller:** the score is discarded
once the rung is picked, and what a consumer reads is held flat.

### 4.4 What the curve records about itself

The held-out error the winning rung achieved is kept, because §7's blend reads it as the curve's
own precision. So is the repeat-count range fitted over, and how many cells stood behind it.

---

## 5. What replaces borrowing, and what borrowing still does

**The level stops borrowing. Nothing else does.**

- **The level** at every cell of a period with a curve comes from §7 — the curve, the cell's own
  fit, or a blend of the two. The nearest-neighbour ring
  ([`ssr_fit.rs:711`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L711)) no longer
  supplies it.
- **The direction split and the fall-off are copied from the nearest stratum that measured them
  well** — §5.1, brought into scope by the owner on 2026-08-20. Until then this route pooled a
  thin stratum's tracts with its neighbours' and **refitted** the pooled set against a floor
  counted in tracts, which is not what §4.5 of
  [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) has ever asked for.
- **The monotonicity merge is not retired here, because it does not exist here.** The joint route
  has borrowing only; `merge_until_monotone` lives in the per-sample route
  ([`ssr/mod.rs:1479`](../../../../src/ng/parameter_estimation/ssr/mod.rs#L1479)). What the curve
  does to that rule is the per-sample route's question.

**What the curve is worth against the rule it replaces, measured.** The dip the merge rule exists
to remove is real but rare and concentrated at the thin end: fitting HG002's cells with nothing
linking them, **6 of 50 steps between neighbouring cells run downhill** — 2 of 22 at period 1, 1 of
19 at period 2, 0 of 3 at period 3, and **3 of 6 at period 4**, where the deepest is a 2.12-fold
drop between two cells of 88 and 104 tracts. On tomato's four steps, none. A curve is monotone by
construction wherever `slope > 0`, so all six go away without a rule.

### 5.1 The shares' floor, and why copying replaces the pooled refit

**The level and the two shares starve at completely different rates, so one floor cannot protect
both.** A stratum of 100,000 tracts at 5 reads each — half a million reads — at a slippage level
of 0.091% has **455 slipped reads**. That measures the level to about 5% of itself, and the
fall-off to about 45%: the same stratum measures one of its numbers well and another barely at
all ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.5).

**So the shares get their own floor, counted in slipped reads rather than tracts.** At the values
§3 of that document measures, holding the direction split to 6% of itself takes about **1,400**
slipped reads and holding the fall-off to the same takes about **4,000**; the fall-off binds, and
[`../arch/parameter_prepass_ssr.md`](../arch/parameter_prepass_ssr.md) §4.1 fixes the floor at
`MIN_SLIPPED_READS_TO_FIT_SHARES = 4_000`. *Soft, and expected to be missed by every stratum at
the bottom of the repeat range — that is the rule working, since the alternative is a share fitted
on five reads and reported as measured.*

**Below the floor the two shares are copied from the nearest stratum at the same period that
clears it** — nearest by repeat count, and where two are equally near, the shorter tract wins
because there are more of them. Motif period is never crossed: the direction split runs 1.4× at
tomato homopolymers and 4.9× at its dinucleotides, so a copy across periods would carry a
threefold error.

**Copying replaces the pooled refit outright, and that is the point.** Once the level comes from
the curve and the shares are copied, **nothing needs a stratum's tracts pooled with its
neighbours' at all**. What that removes is the expensive arm: the perf review measures one run at
**1,036.8 s with pooled borrowing against 155.5 s without**
([`perf_ng-census-joint-fit_2026-08-15.md`](../../reports/reviews/perf_ng-census-joint-fit_2026-08-15.md)
§3), because a borrowing stratum refits about a thousand tracts where an independent one reads
only its own.

**A stratum with reads but no fit of its own is now emitted rather than refused**, which is spec
§1.1's first goal. Its level is the curve's and its shares are copied, and **none of it was
fitted**, so it carries no length spectrum, no concentration and no log-likelihood — there is
nothing to put there and a fitted-looking zero would be a lie. That is a different shape from a
stratum that was fitted, and §8 says how it is told apart.

**Two refusals survive**, and both mean there is genuinely no answer: a stratum no read crossed,
and a stratum whose period drew no curve *and* whose period has no stratum clearing the shares'
floor.

---

## 6. Beyond the repeat counts the curve saw

**Outside the fitted range the curve is held flat at the nearest fitted end, never continued.**
Below the lowest contributing cell the level is that cell's curve value; above the highest, that
one's.

The alternative was to extrapolate, and §2's third fact is why it lost: both fixed shapes produce
impossible numbers a few repeat counts outside their range, and `rise_shape` being fitted does not
protect against that — a period whose cells all sit in a narrow window will fit `rise_shape` near 0
and the curve will explode above it, exactly as HG002's 8–12 window does at 205.

**Holding flat is wrong in a known direction and that is the point:** the level genuinely keeps
rising, so a held value under-states slippage at a tract longer than anything fitted. It is
recorded as out-of-range in the provenance (§8) so a consumer can see it, and it cannot produce a
number that is not a probability.

---

## 7. Where the curve and a cell disagree

**The emitted level is an inverse-variance blend of the cell's own fit and the curve, on the log
scale.** Both quantities are relative errors, so they combine multiplicatively:

```text
cell's relative standard error     ≈ 1 / sqrt(slipped reads in that cell)
curve's relative standard error    =  the held-out error of §4.4

weight ∝ 1 / (relative standard error)²

log(level emitted) = w_cell · log(cell's level) + w_curve · log(curve's level)
```

**This *is* the protection against fitting each cell's noise, and it is not a switch between two
regimes.** The curve carries **93%** of the weight at a cell with 40 slipped reads behind it and
**6%** at one with 8,000, at HG002's homopolymer curve. Using the curve everywhere is the same
formula with the curve's weight pinned at 1; using each cell's own answer is it pinned at 0.
Neither end is a separate design, and neither has to be argued for separately — what has to be
argued is why the fitted middle beats the always-curve end, which is §7.1.

### 7.1 Why not use the curve everywhere

**Because where the curve misses, it misses systematically and by far more than the cell's own
noise, and it misses at the cells holding the most tracts.** Fitted over HG002's 23 homopolymer
cells, the winning curve sits within 0.5 to 12% of every cell from 10 repeats up — but at 8 and 9
repeats it is **27% and 55% high**, against those cells' own sampling errors of 1.8% and 1.7%.
Those two cells hold 4,194 and 2,608 tracts, more than any other, and 8 to 9 repeats is where most
homopolymer loci sit.

**No member of the family repairs it**, so this is not a bad choice of shape number: at every rung
from 0.00 to 1.00 the worst residual falls at 8 or 9 repeats, and the winning rung is the least
bad of them — 55% against 282% at the multiplying end. There is a knee between 9 and 10 repeats
that a two-parameter monotone curve cannot bend around. *That knee is also the most suspicious
step in the sequence — the level jumps 2.45-fold across it — and §11's first open question bears
directly on it.*

**Always-curve would therefore emit 10.4 reads slipping per 1,000 at a 9-repeat homopolymer where
that cell's own 3,520 slipped reads say 6.7.** The blend emits 7.1. That is the whole argument.

**Where the curve is as good as the cells, the blend costs nothing** and says so by itself: at
HG002's dinucleotides the curve's median distance from a cell is 3.5% against the cells' own 3.5%,
so the weights come out near even and the emitted values sit between two answers that agree.

### 7.2 One refinement: a cell may say the curve is wrong about it

**A disagreement far larger than either error explains is evidence about the curve, not about the
cell.** Scale the gap by the two errors combined:

```text
gap = |log(cell's level) − log(curve's level)| / sqrt(cell's error² + curve's error²)
```

and beyond a knee of 2.5 combined errors, divide the curve's weight by `(gap / 2.5)²`. At HG002's
9-repeat homopolymer the gap is **9.3 combined errors**, which no sampling noise produces, and the
curve is stood down.

**This is a refinement and the spec says so rather than dressing it up:** without it the blend
already emits 7.1 reads per 1,000 there against the cell's 6.7 — 5.8% high, not 55% — because a
cell with 3,520 slipped reads outweighs the curve on its own. With it, 6.8. It is worth having
because the bottom of the repeat range is what the copy-floor decision reads
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5.1), and it costs one comparison.

**Three cases fall out of the same formula and none needs a branch:**

- a cell with no fit of its own — `w_cell = 0`, the level is the curve's;
- a period with no curve — `w_curve = 0`, the level is the cell's, as today;
- everything else — a blend, with the curve's share recorded.

**⚠ A cell whose fitted level is exactly zero must not enter the logarithm.** The level is
bounded below by zero and a thin cell piles up against that boundary
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.5). Such a cell contributes nothing to
the curve (§4.1 excludes it — a zero level means zero slipped reads means zero weight) and takes
the curve's value whole.

---

## 8. What is emitted, and what the provenance must carry

**A cell that today is refused can now carry a level, so the emitted object has to be able to say
so.** Today's outcome is either a fit or a refusal
([`ssr_fit.rs:322`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L322)); after this
change a cell can hold a level with no shares, which is neither.

Per cell, beside the four numbers:

- **how much evidence stood behind the cell itself** — tracts, reads crossing, and slipped reads;
- **where the level came from** — the cell's own fit, the curve, or a blend, and for a blend the
  share the curve carried;
- **where the two shares came from** — the stratum's own fit, or the repeat count they were
  copied from (§5.1);
- **whether the cell sat inside the curve's fitted range or beyond it** (§6);
- **the curve's own held-out error and cell count** at that period, so a consumer can tell a
  curve through 23 cells from one through four.

**`Provenance::Borrowed` no longer marks what it used to for the level** — the mechanism that set
it is gone — and dropping it without replacement is the failure the research plan §10 exists to
prevent. The four items above are the replacement.

---

## 9. Cross-cutting concerns

**Cost: this run gets cheaper, not more expensive.** Stage one is the fit with borrowing off,
which on the perf review's own measurement is 155.5 s against 1,036.8 s for the same cohort with
borrowing on
([`perf_ng-census-joint-fit_2026-08-15.md`](../../reports/reviews/perf_ng-census-joint-fit_2026-08-15.md)
§3) — because a borrowing cell refits a pooled set of about a thousand tracts where an
independent cell reads only its own. Stage two is 21 rungs × a two-parameter least-squares fit
over at most a few hundred cells: microseconds beside a fit measured in minutes.

**⚠ Do not delete the pooled-set dedup in `fit_strata`
([`ssr_fit.rs:653`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L653)) as dead code.**
With borrowing off every pooled set is a singleton and the dedup does nothing, so it looks
removable; it is what makes the borrowing arm affordable, and that arm must stay selectable to
re-measure what this change bought.

**Errors.** A period that cannot draw a curve is not an error — it keeps today's behaviour and
says so. The one genuine failure is a curve whose `slope` comes out non-positive, which would
make the level fall with repeat count; that cell set is reported as no curve rather than emitted,
because a falling level contradicts every measurement this design rests on.

**One sample, and a thousand.** The curve is fitted from cells, and a cell pools every sample's
reads at every tract in it, so cohort size enters only through how much evidence each cell holds.
At one sample the curve matters *more*, not less: HG002 is one sample and 65 of its 137 cells
would otherwise have no level.

**Three reads a position, and three hundred.** Tomato at ~3 reads is the case where a cell's own
level is worst determined and the curve carries most of the weight — 3 of its 6 fitted cells sit
below §7's crossover. HG002 at ~300 is where cells win and the curve mostly stands aside. Both
are the blend working rather than two behaviours.

**Concurrency.** Stage two runs after every cell of a period is fitted, so it is a join across
cells that are already fitted independently. It adds no shared state to stage one.

---

## 10. Deferred, with a recommended home

- **One-stage fitting** — fitting the curve's parameters directly against the reads rather than
  through the per-cell answers. Home: the research plan
  [`str_slippage_across_repeat_count.md`](../impl_plan/str_slippage_across_repeat_count.md) D1,
  which asks for the gap to be measured on the held-out criterion before deciding.
- **Smoothing the other three numbers.** Home: the same plan's C4.
- **A read-group axis that carries library identity across a cohort** (§3). Home: the census
  specification [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) — it
  owns what a section is keyed by. Until then a multi-library cohort of single-library samples
  fits one curve.
- **Partial pooling of `intercept` and `slope` across periods** for periods 5 and 6, which on
  both cohorts hold too few cells to fit a line. Home: the research plan's D2.

---

## 11. Open questions

- **⚠ OPEN — is the flattening the polymerase or the census's recording window?** The census
  records a read's length offset only over ±4 repeats and folds everything beyond into an end
  bucket ([`census.rs:398`](../../../../src/ng/parameter_estimation/joint/census.rs#L398)). At 30
  repeats 60 reads in 100 that cross the tract report a length other than the reference's, so the
  end buckets carry a large share of the evidence exactly where the curve bends. **A flattening
  that begins where the buckets saturate is what a recording artefact looks like.** *Leaning: the
  design is unaffected either way — `rise_shape` is fitted, so if the bend is an artefact the
  rung simply moves toward 0 — but the fitted values reported before this is settled must not be
  quoted as chemistry.* **A second reason to run it, found while settling §7.1:** the one place
  no member of the family can follow the cells is the 9→10 step, where the level jumps 2.45-fold
  — and that is where the recording window would first bite. If the knee is an artefact, the
  curve may fit the whole range afterwards and §7.2's refinement stops earning its keep.
  **The measurement that settles it:** widen `RECORDED_OFFSET_RANGE` to 8, re-walk HG002, and
  compare the fitted `rise_shape` and the residual at 8 and 9 repeats. About an hour. **Confirm
  before the numbers are published; not a blocker on the code.**
- **⚠ OPEN — is four contributing cells enough to draw a curve?** §4.1's floor is arithmetic. *The
  measurement that settles it: draw cells at a known curve and at the cell counts real periods
  have, and report the held-out error against the number of contributing cells.* *Leaning: four
  is too few — HG002's period 3 has exactly four and its best rung predicts a held-out cell only
  to 31%, against 3.8% at period 2's twenty.* **Confirm before the floor is treated as measured.**
- **⚠ OPEN — should `rise_shape` be shared across periods as well as across slippage groups?** §3
  shares it across groups on an argument about evidence, and the same argument applies to a
  period with four cells. *Leaning: no — HG002's fitted rungs run 1.00, 0.80, 0.00 and 1.00
  across four periods, and the two extremes are the two periods with fewest cells, so sharing
  would either be harmless or would be letting the thin periods vote.* **The measurement that
  settles it:** score a shared rung against per-period rungs on the same held-out criterion, which
  is the research plan's D2 arm.
- **RESOLVED — which family.** One family with a fitted `rise_shape` (§2), because the two cohorts
  prefer opposite fixed shapes over the same repeat counts and both fixed shapes produce
  impossible values outside their range.
- **RESOLVED — the curve or the cell where they disagree, and why not the curve everywhere.** An
  inverse-variance blend (§7), rather than a rule with a threshold, because both precisions are
  estimable and the pull then falls out rather than being chosen. Always-curve was live and lost
  on a measurement: it would emit 10.4 reads slipping per 1,000 at HG002's 9-repeat homopolymer
  where that cell's own 3,520 slipped reads say 6.7 (§7.1).
- **RESOLVED — beyond the fitted range.** Held flat (§6), because a fitted `rise_shape` near 0
  extrapolates to numbers that are not probabilities.

---

## 12. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the per-cell fit that stage one runs | [`fit_stratum` / `fit_pooled`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L482) | called unchanged, with `borrowing_floor = 0` |
| the per-cell evidence the weight reads | [`StratumEvidence`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L169) | `spanning_reads()` and the fitted level give the slipped-read count |
| borrowing, for the two shares | [`borrow_up_to_the_floor`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs#L711) | kept, narrowed to the shares |
| the walk that produces a cell table | [`ng_joint_records_walk.rs`](../../../../examples/ng_joint_records_walk.rs) | the measurement vehicle and the parity oracle |

**The parity oracle** is the run this design was measured on: the same two cohorts, fitted with
borrowing off, must reproduce their per-cell levels to the precision the cell table prints when
the curve is switched off. A change that moves a cell's *own* fit is a defect in the plumbing,
because nothing here touches stage one.
