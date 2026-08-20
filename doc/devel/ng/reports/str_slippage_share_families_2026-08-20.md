# What shape the two slippage shares follow across repeat count — measured on two cohorts

**Status:** measurement, 2026-08-20. Produced by the first step of Milestone E of
[`../impl_plan/str_slippage_level_curve.md`](../impl_plan/str_slippage_level_curve.md), which
ships no code: its output is the list of curve families the next step implements. The design it
serves is [`../spec/str_slippage_level_curve.md`](../spec/str_slippage_level_curve.md) §5.1.

**In one sentence:** all three candidate shapes are needed — across ten cases a flat value wins
twice, a bending line three times and a curve that bends twice five times — so the shape is chosen
at run time by the same leave-one-cell-out score that already chooses the level's, and the
copy-from-a-neighbour rule it replaces is beaten in 8 of those 10 cases even when that rule is
given evidence it does not have in practice.

---

## Vocabulary

A **stratum** is one cell of the table the slippage parameters are fitted in: all tracts of one
motif period at one repeat count — say every 12-repeat dinucleotide in the genome. The **motif
period** is the length of the repeat unit: 1 for `AAAAAA`, 2 for `ATATAT`.

Each stratum carries three fitted numbers about how often the copying steps before sequencing add
or drop whole repeat units:

- the **level** — how many of the reads crossing the tract come back at a wrong length;
- the **direction split** — of those, the share that came back *shorter* rather than longer;
- the **fall-off** — how much rarer a two-unit slip is than a one-unit slip, and a three-unit than
  a two-unit.

**Slipped reads** is the level times the reads crossing the stratum's tracts: the count of reads
that actually carry information about the two shares. It is what sets how precisely a stratum
measures them.

The **leave-one-cell-out score** is how a shape is judged: fit the shape to every stratum of a
period except one, predict the one left out, and take the median of those errors across all the
strata in turn. It measures what the curve is for — supplying a number to a stratum that has none
of its own.

---

## 1. What was asked, and what was run

The level already gets a curve across repeat count. The two shares do not: a stratum either
measures its own on at least 4,000 slipped reads, or copies both whole from the nearest stratum
that did. That rule reaches 13 of HG002's 137 strata and none at all of tomato's
([the design](../spec/str_slippage_level_curve.md) §5.1.1), so it is being replaced by a curve per
period fitted from every stratum, weighted.

**A curve needs a family, and the level's cannot be reused** — it is monotone and refuses a falling
fit, and both shares fall as often as they rise. The design lists three candidates and says to pick
between them by measurement, per motif period and per parameter:

- a **constant**, the evidence-weighted mean of the period's strata — the answer when there is no
  trend to fit;
- a **line on the logit scale** — a share bending steadily one way as tracts get longer;
- a **quadratic on the logit scale** — a share that bends twice, falling then rising.

*The logit scale is `log(p / (1 − p))`. It is used because both shares are proportions: on that
scale a fitted curve can take any value at all and the share it maps back to still lies strictly
between 0 and 1, so nothing ever needs clamping.*

**The tables.** Both cohorts' per-stratum tables from the census's widened ±8 recording window,
with every stratum fitted from its own tracts and nothing copied from a neighbour:
`hg002_wide8_plain.csv` (a single human sample, GIAB HG002, about 300 reads a position) and
`tomato_d4_plain.csv` (63 tomato accessions at about 3 reads a position each). The script is
`tmp/slippage_curve/e1_share_families.py`.

**Reading a copied share as a measurement would be circular**, which is why the HG002 table used
here is the one made before the copying rule existed. Nothing was copied on tomato in any run —
no tomato stratum ever cleared 4,000 slipped reads, so there was never a lender.

**The script reproduces the earlier report where the two overlap.** Scored the same way — absolute
error, strata weighted by slipped reads — it returns
[`str_slippage_shape_2026-08-20.md`](str_slippage_shape_2026-08-20.md) §5's figures exactly: 0.060
for a logit-line on HG002's dinucleotide direction split against 0.240 for a flat mean, 0.043
against 0.122 for its homopolymer fall-off, 0.176 against 0.225 for its tetranucleotide fall-off.

---

## 2. How few strata the choice rests on

**A period can only be given a curve where at least four of its strata were fitted on their own
tracts** (the design's §4.1 floor: leaving one out still leaves a line and a spare). That is a
harder constraint than it sounds:

| | strata with reads crossing them | fitted on their own tracts | periods that reach four fitted strata |
|---|---:|---:|---|
| HG002, one sample at ~300 reads | 132 | 55 | 1, 2, 3, 4 — of six |
| tomato, 63 accessions at ~3 reads | 49 | 6 | 1 — of six |

Per period, the fitted strata that a curve would be drawn through:

| fitted strata at | period 1 | 2 | 3 | 4 | 5 | 6 |
|---|---:|---:|---:|---:|---:|---:|
| HG002 | 23 | 20 | 4 | 7 | 1 | 0 |
| tomato | 5 | 1 | 0 | 0 | 0 | 0 |

**So this measurement has ten cases to speak from** — four periods on HG002 and one on tomato,
each with two parameters — and two of the five periods are thin: HG002's trinucleotides have
exactly four strata and tomato's homopolymers five.

---

## 3. The result: which family wins, where

Leave-one-stratum-out median error, with strata weighted by the precision the design gives a share
(§5.1, `sqrt((1 − p) / (p · S))` on `S` slipped reads). **Bold is the winner.** The first column is
the rule being retired, scored the same way, and §4 says why it is scored generously.

| cohort | period | fitted strata | parameter | spans | copy the nearest | constant | logit-line | logit-quadratic |
|---|---:|---:|---|---:|---:|---:|---:|---:|
| HG002 | 1 | 23 | direction split | 1.60× | 5.5% | 5.4% | 6.1% | **4.4%** |
| HG002 | 1 | 23 | fall-off | 3.54× | 11.9% | 37.5% | 19.0% | **8.6%** |
| HG002 | 2 | 20 | direction split | 4.52× | **17.6%** | 40.3% | 27.0% | 17.9% |
| HG002 | 2 | 20 | fall-off | 5.17× | **20.5%** | 23.8% | 31.8% | 22.2% |
| HG002 | 3 | 4 | direction split | 2.71× | 61.6% | 76.7% | 121.3% | **5.2%** |
| HG002 | 3 | 4 | fall-off | 2.31× | 39.3% | **26.1%** | 36.8% | 37.5% |
| HG002 | 4 | 7 | direction split | 2.56× | 31.9% | 30.4% | **25.2%** | 52.8% |
| HG002 | 4 | 7 | fall-off | 1.92× | 23.6% | 36.6% | **8.1%** | 29.9% |
| tomato | 1 | 5 | direction split | 1.25× | 11.4% | 11.0% | **8.9%** | 13.7% |
| tomato | 1 | 5 | fall-off | 1.33× | 11.8% | **11.8%** | 15.6% | 27.9% |

**Three findings, and the first is the one the next step needs.**

**All three families earn their place, and no two periods want the same one.** Counting only the
three candidates — the first column is the rule being retired, not a candidate — the constant wins
two of the ten cases, the logit-line three and the logit-quadratic five. Dropping any one of them
costs a real amount somewhere: without the quadratic, HG002's homopolymer fall-off is
predicted to 19.0% instead of 8.6%; without the line, its tetranucleotide fall-off to 29.9% instead
of 8.1%; without the constant, tomato's homopolymer fall-off to 15.6% instead of 11.8%. **So the
family is not a constant of the code — it is fitted at run time, per period and per parameter, by
the same score that already picks the level's shape.**

**Sometimes there is nothing to fit, and the score says so plainly.** Two of the ten cases are won
by the constant, and in both the runners-up are no better than the mean: HG002's trinucleotide
fall-off (26.1% against 36.8% and 37.5%) and tomato's homopolymer fall-off (11.8% against 15.6%).
**This is the answer to the research plan's question of whether each parameter should be smoothed
at all** ([`../impl_plan/str_slippage_across_repeat_count.md`](../impl_plan/str_slippage_across_repeat_count.md),
step C4): sometimes yes and sometimes no, and which one is decided per period by measurement rather
than fixed in advance.

**How well the winner predicts varies more than which winner it is.** The best available error runs
from 4.4% (HG002's homopolymer direction split, fitted through 23 strata) to 26.1% (its
trinucleotide fall-off, fitted through 4). **A stratum with no fit of its own inherits that error**, so a consumer needs to
be told it — which is what the design already requires the emitted provenance to carry (§8: the
curve's own held-out error and how many strata stood behind it).

---

## 4. Against the rule being retired

**The copy-from-the-nearest-stratum rule is beaten in 8 of the 10 cases, and the two it wins it
wins by less than two points** — 17.6% against 17.9% and 20.5% against 22.2%, both on HG002's
dinucleotides. One of the eight is a dead heat: on tomato's homopolymer fall-off the constant
returns 11.76% against the rule's 11.80%. Where it loses it loses heavily — 61.6% against 5.2% on
HG002's trinucleotide direction split, and 23.6% against 8.1% on its tetranucleotide fall-off.

**And that comparison is generous to the rule by a wide margin.** As scored here it may copy from
*any* other fitted stratum. The rule as built may only copy from a stratum holding at least 4,000
slipped reads, and the best-measured stratum of a period holds that many at only one period of the
twelve measured across both cohorts (HG002's homopolymers, 8,840). At every other period the real
rule supplies nothing at all — the 13 strata it reached on HG002 are all homopolymers, and on
tomato it reached none.

**So the honest summary is not "the curve is somewhat better".** Only HG002's homopolymers hold a
stratum above 4,000 slipped reads, so at **8 of the 10 cases the rule being replaced has no number
to give at all** — and the curve gives one predicting a held-out stratum to between 5.2% and
26.1%.

---

## 5. What this delivers, and what it leaves for the refusal floor

Counting strata that end up with a complete parameter set — a level and both shares:

| | strata with reads crossing them | complete today | complete once the shares get a curve | still with nothing |
|---|---:|---:|---:|---:|
| HG002, one sample at ~300 reads | 132 | 68 | **113** | 19 |
| tomato, 63 accessions at ~3 reads | 49 | 6 | **16** | 33 |

**On a deep single sample the curve nearly finishes the job; on a shallow cohort it does not.**
HG002 ends with 19 strata unserved, all at periods 5 and 6, where 1 and 0 strata respectively were
fitted. Tomato ends with 33 unserved — two thirds of everything its reads touched — because only
its homopolymers reach four fitted strata.

**What is blocking tomato is not the curve but what is allowed to feed it.** Nothing below 50
tracts is fitted today, and tomato's dinucleotides have 23 strata carrying reads of which one
clears that floor. Counting how many strata each period would offer at a lower floor — an upper
bound, since it assumes every fit succeeds:

| strata holding at least this many tracts | 20 | 10 | 5 | 3 |
|---|---:|---:|---:|---:|
| tomato, period 2 (23 strata carry reads) | 3 | 6 | 11 | 15 |
| tomato, period 3 (8 strata carry reads) | 0 | 1 | 2 | 4 |
| HG002, period 5 (10 strata carry reads) | 2 | 3 | 5 | 7 |
| HG002, period 6 (10 strata carry reads) | 1 | 2 | 2 | 3 |

**A floor at five tracts would give tomato's dinucleotides a curve through up to 11 strata and
HG002's pentanucleotides one through up to 5**, taking tomato from 16 complete strata to about 39
of 49. **That is the step gated on the open question of whether a fit that thin is usable at all**
([the design](../spec/str_slippage_level_curve.md) §11), and these counts are the size of the prize,
not evidence that the fits work.

---

## 6. What the recording window changed about the shares

The census used to record a read's length offset over only ±4 repeat units and fold everything
beyond into an end bucket; it now records ±8, which agrees with ±12 to within 1.8%
([`str_slippage_shape_2026-08-20.md`](str_slippage_shape_2026-08-20.md) §7). Scoring the same ten
cases on the old tables and the new:

- **The winning family changes in 4 of the 10 cases** — HG002's trinucleotide direction split
  (constant → logit-quadratic), its tetranucleotide direction split and fall-off (both constant →
  logit-line), and tomato's homopolymer direction split (logit-quadratic → logit-line).
- **The measured shares themselves moved much further than that.** HG002's dinucleotide direction
  split ran from 0.828 to 0.950 across the range at ±4 — a 1.15-fold span, near enough flat — and
  runs from 0.217 to 0.983 at ±8, a 4.52-fold span that falls to a minimum at 10 repeats and climbs
  from there. **A share read from a ±4 table is not a noisy version of the right answer; it is a
  different shape.**

**No number in this report comes from a ±4 table**, and none in the design should.

---

## 7. Two decisions this measurement raises

### 7.1 Which precision weights a stratum — the design's formula is inconsistent with its own blend

**The design fixes both a weight and a scale, and they do not match.** §5.1 gives a share's
relative standard error as `sqrt((1 − p) / (p · S))` on `S` slipped reads, and §5.1 also says the
blend runs on the logit scale. The precision that belongs on the logit scale is a different
quantity: `1 / (S · p · (1 − p))` for the variance, so a weight of `S · p · (1 − p)`. The two
differ by a factor of `1 / (1 − p)²`, which is negligible for a share near zero and enormous for
one near one.

**What that does to real strata.** At HG002's dinucleotides, the 23-repeat stratum splits 0.983
shorter on 1,940 slipped reads and the 10-repeat stratum splits 0.223 on 2,262. The design's
formula gives the first stratum **173 times** the weight of the second; the logit-scale precision
gives it **one twelfth**. They do not merely differ in size — they rank the two strata in opposite
orders.

**Neither is uniformly better on the leave-one-out score**: the best error available at a case
improves under the logit-scale weight in 4 of the 10, worsens in 5 and ties in 1. The two largest
moves go opposite ways — the logit-quadratic on HG002's dinucleotide direction split improves from
17.9% to 11.2%, while the logit-line on its tetranucleotide fall-off worsens from 8.1% to 29.5%. *That score is itself computed as a relative error, which is the
log-scale measure, so it is not a neutral judge between the two.*

**Recommendation: fit and blend on the logit scale with the logit-scale weight**, and keep
`sqrt((1 − p) / (p · S))` in the design where it is doing its other job — explaining where the
figure of 4,000 slipped reads came from, which it reproduces to within 3% of the architecture's
own 1,400 and 4,000. The reason is not the score, which does not separate them; it is that a
weighted fit's weight must be the inverse variance of the quantity being fitted, and the quantity
being fitted is the logit. Under the design's formula a stratum counts as best-measured when its
share is nearest to certainty, which is a property of dividing by `p` and not of its evidence.

**RULED (owner, 2026-08-20): the logit-scale weight.** §9 records the rule as built and what it
changes about the table above.

### 7.2 A curve that bends twice, drawn through only four strata

HG002's trinucleotide direction split is won by the logit-quadratic at 5.2% against 76.7% for the
constant, and it has exactly four fitted strata. **With four strata and three coefficients, every
leave-one-out fit passes exactly through the three that remain** — there is no residual left for
the score to notice a bad shape with.

**The evidence that it is not merely interpolating:** all four held-out predictions land within 7%
of the stratum they predict, including the two at the ends of the range, which are extrapolations.
The underlying swing is far larger than sampling noise — the split runs 0.892, 0.330, 0.356, 0.895
across 6 to 9 repeats, where those four strata's own sampling errors are 2.9%, 8.0%, 10.1% and 2.6%
of themselves.

**Recommendation: allow it, with no extra floor beyond the design's existing four.** Two guards
already contain the failure it risks. Outside the fitted repeat range the curve is held flat and
never continued ([the design](../spec/str_slippage_level_curve.md) §6), so a quadratic cannot run
away at the ends. Inside the range, the winning quadratics stay close to the strata they were fitted
from: the furthest any of them travels beyond its own strata is 0.063, HG002's trinucleotide split
dipping to 0.267 where its lowest stratum sits at 0.330.

**What the alternative costs, if you would rather be conservative:** requiring five strata for the
quadratic — one more than its coefficients, the same reasoning that set the floor of four for a
line — changes exactly one of the ten cases, and it changes it from 5.2% to 76.7%.

---

## 8. What the next step takes from this

- **Three families, chosen at run time** by leave-one-stratum-out median error, per motif period
  and per parameter: constant, logit-line, logit-quadratic. None is droppable.
- **A period with fewer than four fitted strata gets no curve**, and its thin strata keep getting
  nothing until the refusal floor moves — 19 strata on HG002, 33 on tomato.
- **The provenance must carry the curve's own held-out error and how many strata it was fitted
  through**, because that error ranges from 4.4% to 26.1% across the ten cases measured here and a
  stratum with no fit of its own inherits it whole.
- **Two answers wanted before code is written:** the weight formula (§7.1) and whether a quadratic
  may be drawn through four strata (§7.2).

---

## 9. The rule as built, and what the approved weight changed

**Settled after §7 was written**, and the numbers here — not §3's — are what the code does. Four
decisions, one of them the owner's:

- **A stratum's weight is the inverse variance of its logit**, `slipped reads × p × (1 − p)`
  (§7.1, ruled by the owner).
- **The held-out score is in logit units too** — the median of `|logit(predicted) −
  logit(measured)|` over leaving each stratum out in turn. It has to be, because the blend weighs
  a curve's error against a stratum's own, and a stratum's own is now a logit-scale error.
- **The flat shape is the weighted mean of the logits**, not of the shares, so that all three
  shapes are least squares on one scale.
- **A tie goes to the simplest shape.**

What each period's strata then choose, asserted in `tests/share_curve_on_real_cells.rs`:

| cohort | period | fitted strata | parameter | shape | predicts a held-out stratum to |
|---|---:|---:|---|---|---:|
| tomato | 1 | 5 | direction split | sloping | 0.240 |
| tomato | 1 | 5 | fall-off | flat | 0.321 |
| HG002 | 1 | 23 | direction split | flat | 0.165 |
| HG002 | 1 | 23 | fall-off | sloping | 0.239 |
| HG002 | 2 | 20 | direction split | turning | 0.516 |
| HG002 | 2 | 20 | fall-off | flat | 0.527 |
| HG002 | 3 | 4 | direction split | turning | 0.210 |
| HG002 | 3 | 4 | fall-off | flat | 0.813 |
| HG002 | 4 | 7 | direction split | flat | 0.789 |
| HG002 | 4 | 7 | fall-off | flat | 0.840 |

*(logit units: 0.2 is about a tenth of the way from a share of 0.5 to one of 0.7.)*

**The headline survives the change and the detail does not.** All three shapes still win somewhere,
so none can be dropped — but under the approved weight the flat shape wins six of the ten rather
than two, and the two that changed hands did so by a hair: HG002's homopolymer direction split and
its dinucleotide fall-off are each within a few hundredths of a logit unit of the next shape up,
and the tie rule keeps the simpler one. **Read a shape as this period's best available
description, not as a fact about the chemistry.**

**What a curve is worth against a stratum that has its own answer.** At HG002's homopolymers the
median stratum holds its own direction split to 0.033 logit units where the period's curve
predicts one to 0.165 — five times better. So a stratum with a fit of its own will keep almost all
of its weight in the blend, and the curve's value is at the strata with no fit of their own —
58 of HG002's 132 populated strata and 10 of tomato's 49 sit at a period that has a curve, of
which 13 and 0 respectively get anything under the rule being retired.
