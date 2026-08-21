# What the slippage curves reach, on both cohorts

**Status:** measurement, 2026-08-21. The last step of Milestone E of
[`../impl_plan/str_slippage_level_curve.md`](../impl_plan/str_slippage_level_curve.md): what the
three curves deliver once every slippage number a stratum carries is its own answer blended with
its motif period's curve. Design: [`../spec/str_slippage_level_curve.md`](../spec/str_slippage_level_curve.md)
§5.1, §8 and §9.

**In one sentence, and it is not the sentence this milestone was expected to produce:** on a deep
single human sample the work's worth is a **wrong parameter at 31% of the repeat loci, corrected** —
the rule it replaced was overwriting the direction split and fall-off of strata that had measured
them perfectly well — and on a 63-accession cohort at three reads a position it is that **the
curves carry a third to a half of the weight** at strata too thin to measure themselves, where on
the deep sample they carry 3%. The coverage it adds is small at both: **2% and 7% of loci**, a
twentieth of what the count of strata suggests.

---

## Vocabulary

A **stratum** is one cell of the table the slippage parameters are fitted in: every repeat tract of
one motif period at one repeat count — say every 12-repeat dinucleotide in the analysed regions.
Each carries three numbers about how often the copying steps before sequencing add or drop whole
repeat units: the **level** (how often a read comes back at the wrong length), the **direction
split** (of the reads that slipped, the share showing a *shorter* tract) and the **fall-off** (how
much rarer a two-unit slip is than a one-unit one).

**A stratum is not a locus, and the difference is the whole of §4.** Strata differ in size by three
orders of magnitude: HG002's 8-repeat homopolymers hold 4,194 tracts and its 38-repeat homopolymers
hold 5. Counting strata weights those equally; counting the tracts in them does not.

Each number can come from three places, and the emitted table says which: the stratum's **own fit**,
its period's **curve**, or a **blend** of the two weighted by how precisely each determines it.

---

## 1. What was run

Both cohorts walked and fitted end to end at the census's ±8 recording window, each fitted twice —
once with the curves drawn and once with them switched off, which is the parity oracle:

- **HG002**, one human sample at about 300 reads a position, over 21.0 Mb of the GIAB tandem-repeat
  regions: 498,524 kept generic positions and 29,787 repeat tracts in 137 strata, of which 132 have
  a read crossing them and 27,399 tracts sit in those.
- **tomato**, 63 accessions at about three reads a position each, over 8.0 Mb: 1,999,404 kept
  generic positions and 4,164 tracts in 71 strata.

Raw output under `tmp/slippage_curve/{hg002,tomato}_e5*.csv` and their logs.

---

## 2. The oracle first: no stratum's own fit moved

**Every stratum fitted before this milestone is fitted to exactly the same three numbers now.**
HG002's no-curve arm against the same arm before the milestone: 55 strata fitted in both, three
numbers each, and **all 165 identical to the last digit**. On tomato, 6 strata in both and all 18
identical.

Nothing in this milestone touches how a stratum is fitted, so a stratum's own answer moving would
be a defect in the plumbing rather than a consequence of the design.

---

## 3. What the rule that was removed had been doing

**A stratum kept its own two shares only if its own slipped reads reached 4,000; below that both
were replaced by one neighbour's, whole.** On HG002 that overwrote **11 of the 55 strata fitted on
their own tracts, and they hold 8,363 of the 26,769 tracts in fitted strata — 31%.**

The two largest are where most repeat loci sit:

| stratum | tracts | its own reads said | it was reported as |
|---|---:|---:|---:|
| homopolymers, 8 repeats | 4,194 | 0.5095 | **0.6919** |
| homopolymers, 9 repeats | 2,608 | 0.5838 | **0.6919** |

**The 8-repeat stratum had measured its own direction split to within 1.6% of itself**, on 3,666
slipped reads, and lost it for a value 36% higher because 3,666 falls short of 4,000. The caller was
being told that 69 of every 100 slipped reads at an 8-repeat homopolymer come back short, where
that stratum's own reads say 51.

**It now keeps its own answer**: its period's curve carries 2% of the weight there and the emitted
split is 0.5124.

---

## 4. What moved, counted in loci

Over the 68 strata that carry parameters under both rules, holding 26,834 tracts:

| | median move | tracts in a stratum that moved more than a tenth |
|---|---:|---|
| the level | 0.4% | 65 of 26,834 — **none to speak of** |
| the direction split | 0.4% | **7,779 of 26,834 (29%)** |
| the fall-off | 0.6% | **7,827 of 26,834 (29%)** |

**The level barely moves because it never had a copy rule to remove** — it has had its curve since
the earlier milestones, and where a stratum measures its own level well the blend leaves it alone.

---

## 5. What coverage this adds, and why the stratum count flatters it

| HG002 | before | now |
|---|---:|---:|
| strata carrying a complete parameter set | 68 | **117** |
| tracts in those strata | 26,834 | **27,363** |
| tracts in strata carrying nothing | 565 | **36** |

**The 49 strata that gained a parameter set hold 529 tracts — 2% of the 27,399 in the run.** A thin
stratum is thin because few loci sit in it, so "68 strata to 117" reads as a large gain and is a
small one. The strata still carrying nothing are 15, at periods 5 and 6, holding 36 tracts between
them: their period has fewer than four strata fitted on its own tracts, so no *level* curve is
drawn. **The level's four-stratum floor is now the only thing that refuses a populated stratum** —
the two shares always have a curve to give.

---

## 6. What the curve is worth where a stratum already has an answer

**On a sample this deep, almost nothing, and that is the design working.** Over HG002's 79 strata
fitted on their own tracts, the curve carries a median of **3.9%** of the weight for the level,
**3.4%** for the direction split and **2.6%** for the fall-off, moving the emitted values by a
median of **0.4%, 0.2% and 0.3%**. The largest single move is 70%, at a stratum whose own fit and
its curve disagree by more than either error explains.

A stratum with thousands of slipped reads behind it should keep its own answer, and it does.

---

## 7. Where each number came from

HG002, over the 117 strata that carry numbers:

| | its own fit whole | a blend | the curve whole |
|---|---:|---:|---:|
| the level | 5 | 74 | 38 |
| the direction split | 0 | 79 | 38 |
| the fall-off | 0 | 79 | 38 |

The five levels taken whole are the strata at periods 5 and 6, whose period has no curve to blend
with — the formula's other end, and what every stratum did before the milestone.

**No stratum anywhere took a share from another motif period or from the built-in default.** The
shares' fallback ladder has four rungs; 112 strata took the top one — their own period's curve with
its shape chosen by measurement — and 5 took the second, a flat mean over a period with too few
strata to choose a shape between. Tomato is the same: 38 of its 39 took the top rung and 1 the
second. **So no stratum on either cohort is furnished from nothing**, which was the question this
step was asked to answer.

**`StratumOutcome::Derived` is used and stays** — 38 strata on HG002 and 22 on tomato are it.

---

## 8. tomato

**On the shallow cohort the curves do the opposite job, and this is where the smoothing itself
earns its place.** Of tomato's 49 strata with a read crossing them, holding 3,965 loci:

| | before | now |
|---|---:|---:|
| fitted on their own tracts | 6 | **17** |
| furnished from their period's curves | 0 | **22** |
| carrying nothing | 43 | **10** |
| loci in strata carrying parameters | 3,661 | **3,941 of 3,965** |

**Where a stratum does have its own fit, the curve carries a third to a half of the weight** —
a median of 34% for the level, 49% for the direction split and 23% for the fall-off, against
3.9%, 3.4% and 2.6% on HG002. At three reads a position a stratum barely measures itself, so its
period's line does much of the work; at 300 reads it does almost none. *That is the same formula
at both ends, and neither end was chosen.*

**Nothing that already carried parameters moved**: the 6 strata fitted under both rules shift by a
median of 0.4% to 1.4% and not one of their 3,661 loci moves by more than a tenth. **The copy rule
never fired on tomato at all** — no stratum ever reached 4,000 slipped reads, so there was never a
neighbour to copy from, which is why this cohort gains coverage where HG002 gained a correction.

**The 10 strata still carrying nothing hold 24 loci**, at periods 3, 4 and 6, and the reason is the
same as HG002's: fewer than four strata at those periods can measure themselves, so no *level* line
is drawn.

---

## 9. What the run costs

| fitting | before | now |
|---|---:|---:|
| HG002 | 366.1 s | **421.0 s** |
| tomato | 549.2 s | **665.0 s** |

**The extra time is the extra strata, at about 4 seconds each**: the floor for fitting a stratum at
all moved from 50 tracts to 8, which admitted 24 more on HG002. *Against the arm that pooled a thin
stratum's tracts with its neighbours' and refitted — deleted earlier in this milestone — the same
cohort took 1,036.8 s.*
