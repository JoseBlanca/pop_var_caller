# What shape slippage follows across repeat count — measured on two cohorts

**2026-08-20, corrected the same day.** *Programs: `examples/ng_joint_records_walk.rs`, which walks
real alignments, fills the census and fits every (motif period, repeat count) cell; the per-cell
table it writes under `SSR_CELL_TABLE`. Every stratum speaks from its own tracts alone — nothing
is borrowed, pooled or smoothed before the smoothing is measured.*

> **⚠ CORRECTED. The first version of this report measured everything through a census that
> recorded a read's length offset over only ±4 repeats, and that window was losing real slippage
> at long tracts** — by a factor of **2.26 at 30-repeat homopolymers**, growing with tract length.
> It made the rise appear to flatten above about 20 repeats when it does not, and two of the
> report's headline claims were artefacts of it. The recording window is now ±8
> ([`census.rs`](../../../../src/ng/parameter_estimation/joint/census.rs)), verified converged
> against a ±12 arm that agrees within 1.8% at every repeat count. **Every number below is a ±8
> measurement.** §7 records what the correction changed, because two of the things it overturned
> were this report's own conclusions.

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
| 8 | 4,194 | 798,878 | 3,666 | 0.00459 | 0.51 | 0.54 |
| 10 | 1,883 | 384,622 | 6,574 | 0.01709 | — | — |
| 13 | 1,173 | 238,751 | 8,840 | 0.03703 | 0.65 | 0.39 |
| 16 | 783 | 136,712 | 7,916 | 0.05790 | — | — |
| 20 | 373 | 54,384 | 4,598 | 0.08455 | — | — |
| 25 | 195 | 21,831 | 3,241 | 0.14844 | — | — |
| 30 | 68 | 6,250 | 1,697 | 0.27157 | 0.82 | 0.74 |

*Seven of the 23 shown; the shares are quoted at the ends and the middle only, since §5 reads
their trend rather than their values.*

**Tomato, homopolymers — five cells:**

| repeats | tracts | reads crossing | slipped reads | level | shorter share | fall-off |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 2,082 | 937,558 | 2,338 | 0.00249 | 0.550 | 0.677 |
| 9 | 887 | 363,500 | 1,447 | 0.00398 | 0.648 | 0.757 |
| 10 | 350 | 121,244 | 485 | 0.00400 | 0.600 | 0.626 |
| 11 | 153 | 34,904 | 236 | 0.00676 | 0.678 | 0.569 |
| 12 | 83 | 15,977 | 172 | 0.01079 | 0.687 | 0.714 |

**Slipped reads is `level × reads crossing`, not the count of reads sitting off the reference
length.** The second is much larger and at a polymorphic tract is mostly genuine allele length.

---

## 3. The rise, and what shape fits it

**On HG002 the level rises 59-fold over 8→30 repeats** — 4.6 reads slipping per 1,000 at 8 repeats
to 272 at 30 — **and it does not flatten**. Step to step it is still rising 1.09 to 1.32-fold at
the top of the range.

**Every family was scored by leaving each cell out in turn, fitting the rest, and predicting the
one left out.** Fit quality on the cells a family saw decides nothing. Writing `shape` for the
number in `level ^ shape = intercept + slope · repeat count`, with 0 read as the exponential
(each repeat multiplies the level) and 1 as a straight line (each repeat adds to it):

| | cells | repeats | level spans | best shape | its held-out median | exponential | straight line |
|---|---:|---|---:|---:|---:|---:|---:|
| tomato, period 1 | 5 | 8–12 | 4.3× | **0.65** | **7.3%** | 14.7% | 26.2% |
| HG002, period 1 | 23 | 8–30 | 59.2× | **0.35** | **7.7%** | 18.8% | 21.8% |
| HG002, period 2 | 20 | 6–25 | 63.4× | **0.70** | **11.4%** | 15.5% | 21.9% |
| HG002, period 3 | 4 | 6–9 | 4.7× | 1.00 | 32.0% | 88.5% | 32.0% |
| HG002, period 4 | 7 | 6–12 | 3.0× | 1.00 | 21.3% | 28.0% | 21.3% |

**At every period with enough cells to mean anything, the fitted shape beats both ends** — and it
lands *between* them. Periods 3 and 4 predict a held-out cell to 32% and 21%; read them as noise,
not as answers.

**This is the finding that fixes the family.** Neither fixed shape is right: the exponential loses
at every period, and the straight line loses at the three that carry evidence. A family with a
fitted shape number covers all of them, and the shape is a number the data supplies rather than a
choice the design makes.

**The measuring machinery is not biased toward any family.** Run against a table generated with a
known exponential rise and 12% scatter, it picks the exponential (held-out median 9–11%, which is
the scatter) and rejects the straight line at 70–88%.

---

## 4. Do the two cohorts agree?

Restricted to 8–12 repeats — five cells each, the same window:

| | levels 8→12, reads slipping per 1,000 | best shape | exponential | straight line |
|---|---|---:|---:|---:|
| tomato | 2.5, 4.0, 4.0, 6.8, 10.8 | **0.65** | 14.7% | 26.2% |
| HG002 | 4.6, 7.4, 17.1, 23.1, 30.0 | **1.00** | 26.7% | 4.4% |

**They do not land on the same shape, and neither is at the exponential end.** Over its full range
HG002 fits 0.35 and tomato fits 0.65, so the two are nearer each other than either is to a fixed
family — but on the shared window they still differ, and the levels differ by about two-fold at
every repeat count.

**What this cannot establish.** The two runs differ in five things at once — library preparation,
species, region set, read depth and cohort size — so the difference cannot be attributed to the
chemistry. Testing that needs two libraries over the same loci, and the census cannot express it
today: a read group is an index within one sample, so tomato's 63 declared libraries all arrive as
read group 0.

---

## 5. The two shares, and what they ask of a curve

**Read from each stratum's own fit, with nothing copied.** Scored the same way — leave a cell out,
fit the rest, predict it — over a constant, a straight line, a logit-line and an isotonic fit:

| | spans | best family | its held-out median | a flat mean |
|---|---:|---|---:|---:|
| tomato p1, shorter share | 1.25× | isotonic | 0.047 | 0.085 |
| tomato p1, fall-off | 1.33× | **flat** | 0.074 | 0.074 |
| HG002 p1, shorter share | 1.60× | isotonic | 0.033 | 0.038 |
| HG002 p1, fall-off | 3.54× | logit-line | **0.043** | 0.122 |
| HG002 p2, shorter share | 4.52× | logit-line | **0.060** | 0.240 |
| HG002 p2, fall-off | 5.17× | **flat** | 0.119 | 0.119 |
| HG002 p4, fall-off | 1.92× | logit-line | 0.176 | 0.225 |

**Three things follow, and they are the reason the shares get a per-period family rather than one
rule.**

- **Sometimes there is a trend worth fitting and it is large.** HG002's dinucleotide direction
  split spans 4.52-fold and a logit-line predicts a held-out cell four times better than the mean.
- **Sometimes there is none.** Tomato's fall-off and HG002's dinucleotide fall-off are predicted
  as well by a flat mean as by anything else, which is the honest answer for them.
- **The trends are not the same direction at different periods.** HG002's homopolymer direction
  split rises 0.51 → 0.82; its dinucleotide split falls to 0.61 by the middle of the range and
  climbs to 0.97. **A single family imposed on all of them would be wrong somewhere.**

**The substitution rate rises with repeat count on HG002 — 6.33-fold across the 23 cells.** Part of
this is mechanical: the STR pre-pass counts a mismatch against the motif tiled to the read's
length, so an interruption inside a long tract is charged here. Recorded, not pursued.

---

## 6. How often the fitted level dips — the question the specification records as unmeasured

[`../spec/parameter_prepass_ssr.md`](../spec/parameter_prepass_ssr.md) §4.3 asks how often the
monotonicity merge would fire when the truth is monotone. With every cell fitted independently and
no constraint anywhere:

| | steps between neighbouring cells | downward | deepest dip |
|---|---:|---:|---:|
| HG002, period 1 | 22 | 2 | 1.08-fold |
| HG002, period 2 | 19 | **5** | 1.26-fold |
| HG002, period 3 | 3 | 0 | — |
| HG002, period 4 | 6 | **3** | 1.31-fold |
| tomato, period 1 | 4 | 0 | — |

**Ten downward steps in fifty, and they concentrate where the cells are thinnest** — half of period
4's steps and a quarter of period 2's, against 2 of 22 at period 1 and none on tomato. **No dip is
deeper than 1.31-fold**, which is what a curve fitted through all of a period's cells absorbs
without a rule.

**On the joint route the merge rule does not exist to fire.** `joint/ssr_fit.rs` never had it;
`merge_until_monotone` lives in the per-sample route
([`ssr/mod.rs`](../../../../src/ng/parameter_estimation/ssr/mod.rs)). The count above is what that
rule *would* do if this route had it.

---

## 7. What the recording window changed, including two of this report's own conclusions

The census records a read's length offset over a fixed window and folds everything beyond into an
end bucket, scored by its marginal. That marginal's justification was measured "at a recorded range
of ±1 on a stratum whose alleles reach three repeats either side" — nothing like a 30-repeat
homopolymer where 60 reads in 100 report a length other than the reference's.

**Widening the window from ±4 to ±8 raises the measured level, and the error grows with tract
length:**

| homopolymer repeats | at ±4 | at ±8 | at ±12 |
|---:|---:|---:|---:|
| 10 | 16.5 | 17.1 | 17.1 |
| 20 | 67.5 | 84.6 | 83.5 |
| 25 | 107.4 | 148.4 | 146.1 |
| 30 | 120.3 | 271.6 | 270.3 |

*(reads slipping per 1,000)*. **±12 agrees with ±8 to within 1.8% at every repeat count**, so the
widening has converged.

**Two of this report's first conclusions were artefacts of ±4 and are withdrawn:**

- **"The rise decelerates, and a straight line beats the exponential."** At ±4 the step-to-step
  ratio fell to 1.02 by 20 repeats and the fitted shape was 1.00. At ±8 it is still rising
  1.09–1.32 at the top and the shape is 0.35. **The flattening was the bucket, not the polymerase.**
- **"The two cohorts prefer opposite shapes."** At ±4 tomato fitted 0.00 and HG002 1.00 — the two
  ends. At ±8 they fit 0.65 and 0.35. They still differ; they are no longer opposite.

**A third was distorted rather than reversed.** HG002's homopolymer fall-off looked **U-shaped** at
±4 — 0.43 at 8 repeats, 0.15 at 12, back to 0.38 at 30. At ±8 it rises 0.54 → 0.74 with a dip in
the middle a logit-line fits well. Anything designed around that U specifically should be re-read.

**What the correction costs.** The offsets are a dense array a locus a read group, so exactly
`2·(2n+1)` bytes — **18.2 at ±4 against 34.2 at ±8**, measured; a HG002 sample's census goes 4.13 MB
to 4.61 MB. The owner accepted that on 2026-08-20. **A window scaled to the tract's own repeat
count would be cheaper and exact** — a 6-repeat tract cannot lose more than 6 — and needs a
variable-width record; deferred to the census design.

---

## 8. What was changed to produce this

- `StratumEvidence` gained `bases_compared` and `mismatching_bases`, filled by `gather_strata`
  from the census sections that already record them, plus `substitution_rate()` and
  `reads_off_reference_length()`.
- The walk writes a per-cell CSV under `SSR_CELL_TABLE`, with a second arm under
  `SSR_CELL_TABLE_NO_CURVE` that fits the same cells with no curve drawn — the parity oracle.
- `RECORDED_OFFSET_RANGE` moved 4 → 8 (§7).

**Which table each number here came from.** The levels and shapes are the `_plain` arm of the ±8
runs, where every stratum is fitted from its own tracts. The two shares are read from a run made
**before** the share-copying rule existed, or from tomato, where no period had a stratum clearing
the copy rule's floor and so nothing was copied — otherwise a copied share would have been read as
a measurement.

**The pooled-borrowing arm was never run and never will be**: pooling has since been deleted
([`ssr_fit.rs`](../../../../src/ng/parameter_estimation/joint/ssr_fit.rs)), on the reasoning that a
curve through all of a period's cells is a better answer than one neighbour's, and it removed the
run's expensive arm — the perf review measures 1,036.8 s against 155.5 s for the same cells fitted
independently.
