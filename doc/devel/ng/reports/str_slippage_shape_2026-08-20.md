# What shape slippage follows across repeat count — measured on two cohorts

**2026-08-20.** *Programs: `examples/ng_joint_records_walk.rs`, which walks real alignments,
fills the census and fits every (motif period, repeat count) cell; the per-cell table it now
writes under `SSR_CELL_TABLE`. Both runs used `SSR_BORROWING_FLOOR=0`, so every cell speaks from
its own tracts alone and nothing is smoothed before the smoothing is measured.*

This is step A of [`../impl_plan/str_slippage_across_repeat_count.md`](../impl_plan/str_slippage_across_repeat_count.md)
and the evidence behind [`../spec/str_slippage_level_curve.md`](../spec/str_slippage_level_curve.md).

---

## 1. What was run

| | tomato | HG002 |
|---|---|---|
| samples | 63 accessions | 1 |
| reads a position | ~3 | ~300 |
| analysed regions | `benchmarks/tomato1/regions.bed`, 80 spans, 8 Mb | the GIAB tandem-repeat benchmark's 50,000 Tier intervals, padded 150 bp and merged: 49,701 spans, 21 Mb |
| alignments | `benchmarks/tomato1/crams/*.bench.cram` | `benchmarks/ssr_hg002/bam/300x/HG002_TR_v1.0.1_Tier_300x.bam` |
| repeat tracts kept | 4,164 in 71 cells | 29,787 in 137 cells |
| reads crossing a tract | 1,587,705 | 5,221,537 |
| reads reaching one without crossing it | 1,538,184 | 4,040,816 |
| cells fitted on their own tracts | **6** | **55** |
| repeat-tract fit | 666 s | 492 s |

**Tomato does not meet the research plan's gate for asking the shape question, and no larger run
is available.** The plan asks for three periods with four or more consecutive populated cells;
tomato gives **one** period with five. Its CRAMs are sliced to those 80 spans — a read count of
zero in a region outside them, checked — so there is no more tomato genome to walk. **Every
cross-period statement below therefore rests on HG002**, and is labelled that way rather than
generalised.

**What the human region set is and is not.** It is chosen to *be* tandem repeats and to be the
loci an assembly-based truth set is confident about, which is what makes 137 populated cells
possible at all. It is not a sample of the genome, so the *mix* of cells is not representative;
the level within a cell is a property of the tracts in it.

---

## 2. The cells

**HG002, homopolymers — 23 consecutive cells, nothing linking them:**

| repeats | tracts | reads crossing | slipped reads | level | shorter share | fall-off |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 4,194 | 798,878 | 2,967 | 0.00371 | 0.630 | 0.426 |
| 9 | 2,608 | 522,983 | 3,520 | 0.00673 | 0.642 | 0.256 |
| 10 | 1,883 | 384,622 | 6,331 | 0.01646 | 0.726 | 0.167 |
| 12 | 1,337 | 281,382 | 7,874 | 0.02798 | 0.743 | 0.146 |
| 15 | 953 | 175,436 | 8,056 | 0.04592 | 0.732 | 0.195 |
| 20 | 373 | 54,384 | 3,671 | 0.06750 | 0.645 | 0.247 |
| 25 | 195 | 21,831 | 2,344 | 0.10737 | 0.583 | 0.307 |
| 30 | 68 | 6,250 | 752 | 0.12025 | 0.497 | 0.379 |

*Eight of the 23 shown; the full table is the run's own `SSR_CELL_TABLE` output.*

**Tomato, homopolymers — five cells, and they reproduce the 2026-08-13 run:**

| repeats | tracts | reads crossing | slipped reads | level | shorter share | fall-off |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 2,082 | 937,558 | 1,972 | 0.00210 | 0.598 | 0.635 |
| 9 | 887 | 363,500 | 1,037 | 0.00285 | 0.611 | 0.633 |
| 10 | 350 | 121,244 | 478 | 0.00394 | 0.602 | 0.608 |
| 11 | 153 | 34,904 | 232 | 0.00663 | 0.696 | 0.518 |
| 12 | 83 | 15,977 | 163 | 0.01021 | 0.736 | 0.735 |

**Slipped reads is `level × reads crossing`, not the count of reads sitting off the reference
length.** The second is much larger and is mostly genuine allele length: at HG002's 30-repeat
cell, 60 reads in 100 that cross the tract report a length other than the reference's and the fit
attributes 12 of those 60 to slippage.

---

## 3. The rise decelerates, and that is what decides the family

**On HG002 the level rises 37-fold over 8→30 repeats and the rise flattens as it goes.** Step by
step the ratio runs 1.81, 2.45, 1.33, 1.27, 1.22, 1.20, 1.13, 1.11, 1.14, 1.09, 1.04, 1.02, 1.04,
1.08, 1.10, 1.18, 1.08, 1.10, 0.94, 1.10, 1.14, 0.87.

**Every family was scored by leaving each cell out in turn, fitting the rest, and predicting the
one left out.** Fit quality on the cells a family saw is reported beside it and decides nothing.
Weighted by slipped reads, HG002 homopolymers, 23 cells:

| family | held-out median | held-out worst |
|---|---:|---:|
| straight line in repeat count | **4.4%** | 61% |
| saturating exponential | **4.3%** | 92% |
| isotonic (monotone, no shape assumed) | 11.3% | 81% |
| power law in repeat count | 17.7% | 201% |
| log-linear — the exponential | 22.7% | 345% |
| GATK DRAGstr's per-base hazard | 28.6% | 511% |
| flat (the mean) | 52.5% | 1,300% |

**The ranking does not depend on the weight.** The winner's held-out median is 5.13% unweighted,
4.77% weighted by tracts, 4.65% by reads crossing and 4.39% by slipped reads, and the order of the
families is identical under all four.

### 3.1 Where the winning curve misses, and by how much against the cells' own noise

**A held-out median hides where a family fails.** Fitted over all 23 HG002 homopolymer cells, the
winning curve's distance from each cell, beside that cell's own sampling error
(`1 / sqrt(slipped reads)`):

| repeats | tracts | slipped reads | cell's own error | curve's distance | ratio |
|---:|---:|---:|---:|---:|---:|
| 8 | 4,194 | 2,967 | 1.8% | **27%** | 15× |
| 9 | 2,608 | 3,520 | 1.7% | **55%** | 33× |
| 10–19 | — | 4,332–8,117 | 1.1–1.5% | 0.5–4.4% | 0.5–3.9× |
| 20–23 | — | 2,933–3,671 | 1.7–1.8% | 7.5–12.0% | 4.1–6.5× |
| 24–30 | — | 752–2,728 | 1.9–3.6% | 1.9–10.3% | 0.7–3.7× |

**The two worst cells are the two holding the most tracts**, and they sit where most homopolymer
loci are. The miss is one-directional — the curve says 4.7 and 10.4 reads slipping per 1,000 at 8
and 9 repeats where the cells fit 3.7 and 6.7.

**No rung of the family repairs it.** At every shape number from 0.00 to 1.00 the worst residual
falls at 8 or 9 repeats, and the winning rung is the least bad: 55% against 282% at the
multiplying end. There is a knee at the 9→10 step — the level jumps 2.45-fold across it, the
largest step in the sequence — that a two-parameter monotone curve cannot bend around.

**At the other two periods the curve is as accurate as the cells.** Median curve distance against
median cell error: HG002 dinucleotides 3.53% against 3.51%, tomato homopolymers 6.11% against
4.57%. So the failure above is one period's low end, not a property of the approach.

**The measuring machinery is not biased toward any of them.** Run against a table generated with
a known exponential rise and 12% scatter, it picks the exponential (held-out median 9–11%, which
is the scatter) and rejects the straight line at 70–88% and DRAGstr's form at 67–71%.

---

## 4. The two cohorts prefer opposite shapes over the same repeat counts

Restricted to 8–12 repeats — five cells each, the same window:

| | levels 8→12, reads slipping per 1,000 | exponential | straight line |
|---|---|---:|---:|
| tomato | 2.1, 2.9, 3.9, 6.6, 10.2 | **12.4%** | 33.6% |
| HG002 | 3.7, 6.7, 16.5, 22.0, 28.0 | 31.2% | **8.0%** |

**So the disagreement is not that one cohort saw a wider window.** At the same repeat counts the
human library rises faster and then bends; the tomato one compounds.

**What this does not establish.** The two runs differ in five things at once — library
preparation, species, region set, read depth and cohort size — so the shape difference cannot be
attributed to the chemistry. Testing that needs two libraries over the same loci, and the census
cannot express it today: a read group is an index within one sample, so tomato's 63 declared
libraries (`LB:PRJNA454805_SRR…`) all arrive as read group 0 and the run prints `1 read groups in
1 slippage group, pooled`.

### 4.1 One family covers both

Fitting `level ^ rise_shape = intercept + slope · repeat count`, with `rise_shape = 0` read as
the exponential, over 21 rungs from 0.00 to 1.00:

| | cells | repeats | best `rise_shape` | its held-out median | exponential | straight line |
|---|---:|---|---:|---:|---:|---:|
| tomato, period 1 | 5 | 8–12 | **0.00** | 12.4% | 12.4% | 33.6% |
| HG002, period 1 | 23 | 8–30 | **1.00** | 4.4% | 22.7% | 4.4% |
| HG002, period 2 | 20 | 6–25 | **0.80** | **3.8%** | 18.3% | 5.9% |
| HG002, period 3 | 4 | 6–9 | 0.00 | 31.1% | 31.1% | 53.2% |
| HG002, period 4 | 7 | 6–12 | 1.00 | 53.3% | 54.7% | 53.3% |

**At period 2 the fitted rung beats both ends** — 3.8% against 5.9% and 18.3%. Read that margin
as "the flexible family is not worse": the rung is chosen on the same held-out score it is scored
by.

**Period 4 is noise and should be read as such.** Its best rung predicts a held-out cell to 53%,
which is worse than any period-1 or period-2 family; seven cells at 52 to 332 tracts do not
determine a curve.

---

## 5. The other three numbers, and why this report does not settle them

**The fall-off trends on one cohort and not on the other.** Across tomato's five cells it spans
1.42-fold and a flat mean predicts a held-out cell to within 0.023 — as well as anything else.
Across HG002's 23 it spans **3.14-fold**, falling from 0.426 at 8 repeats to 0.146 at 12 and
rising again to 0.379 at 30, and the flat mean is the worst of four candidates (0.062 against
isotonic's 0.030 and a logit-linear's 0.028). **A U-shape is not what any monotone smoother
describes**, which is why this is left open rather than answered here.

**The direction split trends in opposite directions.** Tomato's shorter share rises 0.598 → 0.736
over 8→12 repeats. HG002's rises 0.630 → 0.743 by 13 repeats and then falls to 0.497 by 30. On
HG002 a logit-linear fit predicts a held-out cell to 0.024 against the flat mean's 0.059, so
there is a trend to fit; it is not the same trend in both cohorts.

**The substitution rate rises with repeat count on HG002** — 6.33-fold across the 23 cells, from
0.0012 to 0.0078 — where nothing in the specification expects it to depend on repeat count at all.
Part of this is mechanical: §4.1 of the STR pre-pass spec counts a mismatch against the motif
tiled to the read's length, so an interruption inside a long tract is charged here. It is recorded
and not pursued.

---

## 6. How often the fitted level dips — the question the specification records as unmeasured

[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3 asks how often the
monotonicity merge would fire when the truth is monotone. With every cell fitted independently
and no constraint anywhere:

| | steps between neighbouring cells | downward | deepest dip |
|---|---:|---:|---:|
| HG002, period 1 | 22 | 2 | 1.15-fold |
| HG002, period 2 | 19 | 1 | 1.12-fold |
| HG002, period 3 | 3 | 0 | — |
| HG002, period 4 | 6 | **3** | **2.12-fold** |
| tomato, period 1 | 4 | 0 | — |

**Six downward steps in fifty, and they concentrate where the cells are thinnest** — half of
period 4's steps, against 2 of 22 at period 1. The deepest is between cells of 88 and 104 tracts.

**On the joint route the merge rule does not exist to fire.** `joint/ssr_fit.rs` has borrowing
only; `merge_until_monotone` lives in the per-sample route
([`ssr/mod.rs:1479`](../../../../src/ng/parameter_estimation/ssr/mod.rs#L1479)). The count above
is what that rule *would* do if the joint route had it.

---

## 7. What would change these numbers

**The census records a read's length offset over ±4 repeats only**
([`census.rs:398`](../../../../src/ng/parameter_estimation/joint/census.rs#L398)) and folds
everything beyond into an end bucket. The share of crossing reads sitting off the reference length
rises from 3.5% at 8 repeats to 60% at 30, so the end buckets carry most of the evidence exactly
where §3's curve bends. **A flattening that begins where the buckets saturate is what a recording
artefact looks like, and this report cannot tell the two apart.** Widening that constant to 8 and
re-walking HG002 — about an hour — is the measurement; until it is run, the fitted `rise_shape`
values above should not be quoted as facts about chemistry.

---

## 8. What was changed to produce this

Two additions to the library and one to the walk, none of which moves a fitted number:

- `StratumEvidence` gained `bases_compared` and `mismatching_bases`, filled by `gather_strata`
  from the census sections that already record them, plus `substitution_rate()` and
  `reads_off_reference_length()`.
- The walk writes a per-cell CSV under `SSR_CELL_TABLE` — the evidence counts beside the fit —
  and, under `SSR_CELL_TABLE_BORROWED`, fits the same cells a second time with borrowing on and
  writes that table too, so "how far borrowing moved this cell" is a difference within one walk
  rather than a comparison across runs.

**The borrowed arm has not yet been run on either cohort.** On tomato it is the expensive one:
the perf review measures a borrowing run at 1,036.8 s against 155.5 s for the same cells fitted
independently, and at 63 samples that extrapolates to hours.
