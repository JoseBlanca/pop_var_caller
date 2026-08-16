# The contamination fraction was half of the truth because of how the depth was written down

*2026-08-16. The share of a sample's reads that came from another individual has been coming back at
about a third of what it should be. The floor under it was fixed three days ago — clean samples now
read zero — but the value itself was still short, and the record said the shortfall was unexplained.
It is explained: the store of evidence the estimate is made from keeps each position's read count as
a code that stands for a range of depths, and the range is wider than the thing being measured. Every
depth now gets a code of its own, and a sample drawn at 3% contaminated reads 2.6% instead of 1.2%.*

*Program: `examples/ng_joint_contamination_control.rs`. Raw output: `tmp/contamination_depth/`.*

---

## 1. What the machine had to work with

Before any variant is called, one pass over each sample's reads writes down what that sample showed
at a fixed set of positions — the same positions in every sample. That store is the **census**, and
the contamination estimate is made from it and nothing else. At each position a sample's census holds
two things:

- **how many of its reads disagreed with the reference, exactly**, and which base they showed;
- **how deep the position was**, as a code on a ladder of depth bands rather than as a number.

The ladder was there for the *other* route to the same parameters, which keeps a table of counts per
band and whose memory grows with the square of the depth it resolves. Bands cost cells there, so the
ladder is exact to 8 reads and widens above: 10, 13, 17, 22, 28, 36, 46, 59, 75, 97, 124. The census
inherited it, and in the census a band costs nothing of the kind — one code per position, and the
only price of a finer ladder is the width of that code.

**At 30 reads a position that code said "somewhere between 29 and 36".** The count of disagreeing
reads beside it was exact. So the share of reads that disagreed — which is what a contamination
fraction is made of — carried a **±12% uncertainty that came from nothing but the record**, while a
3% contamination moves that share by about **1.5%**. The signal sat inside the resolution.

---

## 2. The measurement

A drawn panel, because a real one has no truth to check against: **40 samples in four subpopulations
diverged at F_st 0.20, 400,000 positions, 1 position in 30 planted as mismapped, and sample 0 drawn
genuinely 3% contaminated by a plant taken from the whole panel.** The same program draws the reads,
writes a census from them, fits every parameter and then estimates contamination, so the only thing
that differs between the rows below is what it was asked to do.

| reads a position | how the depth is written down | sample 0, drawn at 3% | worst of the 39 clean |
|---|---|---:|---:|
| 3 | exact — the ladder is exact to 8 | 0.0102 | 0.0003 |
| 8 | exact | 0.0198 | 0.0000 |
| 30 | a band: the code stands for 29–36 | 0.0120 | 0.0000 |
| **30** | **exact: one code per depth** | **0.0263** | **0.0000** |

**The value more than doubles and the floor does not move.** 0.0263 against a truth of 0.030 is 88%
of it, where the band gave 40%.

Two things this settles that were open:

- **It is not low coverage.** The reasonable guess was that at three reads a sample's own genotype is
  mostly guessed, so a heterozygote can absorb the stray reads. That would disappear at 30 reads and
  it does not: the band at 30 reads returns the same 40% as three reads does.
- **It is not the population structure.** The same panel drawn with no divergence at all (one
  subpopulation, F_st 0) and 30 reads a position returns 0.0108 — the same shortfall — so the fitted
  allele frequencies are not what is losing it.

**The three readings of a sample's own place on the panel's axes also stop mattering.** They were
0.0102, 0.0153 and 0.0183 at three reads a position, a trade in which everything that moved the value
towards the truth moved the clean samples up with it. At 30 reads with an exact ladder they are
0.0263, 0.0262 and 0.0263 — one answer. The apparatus built to work around a coarse depth is not
load-bearing once the depth is not coarse.

---

## 3. What was changed, and what it costs

**The census ladder now has one bin per depth to the cap of 124, then the same ten rungs above it**
that carried it to about 1,500 reads a position. Exact wherever the record itself is exact: above the
cap a position's allele counts are thinned to it, so up there the record is approximate however the
depth is written. The histogram route keeps its twenty-bin ladder, and every one of that ladder's
edges is still an edge of this one, so a census code still becomes a histogram bin by collapsing.

The price is the width of the code: **five bits to eight**, because 135 rungs plus the never-walked
sentinel do not fit in 32 values.

| | per sample, per read group, at 2 M positions | across 63 samples |
|---|---:|---:|
| depth array, five bits | 1.25 MB | 79 MB |
| depth array, eight bits | **2.00 MB** | **126 MB** |
| the same sample's list of disagreeing reads at 30× | 2.7 MB | — |

So the depth array grows by 60% and the generic census by something under a third. **Nothing changes
for a shallow cohort's numbers** — at three reads a position both ladders are exact through the whole
range the data reaches — but the memory is spent there too, since the ladder is one of the twelve
values every sample in a run must record under.

A census written by an earlier build cannot be read: its depth array is a different number of bytes
for the same position count and its codes index a different ladder. The file version says so —
`VERSION` is 2 — and the answer to it is to rebuild.

---

## 4. Leaving a sample out of its own allele frequency: built, measured, off

Each position's allele frequency is fitted as a straight line in the panel's axes of variation, using
**every** sample — including the one about to be judged against it. So a contaminated sample's stray
reads have already moved that frequency towards what its own reads show, by a share that depends only
on the panel size: `(components + 1) / samples`, an eighth on a panel of 40 and two fifths on a panel
of 12. Taking one sample back out of a least-squares fit has a closed form and the share to remove is
a number this module already computes for its own refusal, so it costs one multiplication a position.

The expectation was that it would matter most on a small panel. **Measured, it is worth nothing on
either.** Same drawn panel, 30 reads a position, sample 0 at 3%:

| panel | without | with |
|---|---:|---:|
| 40 samples, four subpopulations | 0.0263 | 0.0263 |
| 12 samples, two subpopulations | 0.0248 | 0.0248 |

**And at three reads a position it destroys the finding.** There it does move the value — 0.0102 to
0.0174, nearer the truth of 0.030 — and it moves the clean samples further: the worst of the 39 goes
from 0.0003 to **0.0192**, which is *above* the contaminated sample. A cohort at tomato's depth would
flag the wrong accession.

### Why it does nothing, when the argument for it was sound

The argument was: the sample's own reads move the frequency it is judged against, so take them out.
That is true of the *fit*. It is not true of the *likelihood*, and the difference is where the
frequency is used.

**A sample's own allele frequency enters this estimator in exactly one place: as the prior over that
sample's own genotype.** What a genotype predicts about the reads — a homozygote shows the allele at
the error rate, a heterozygote at a half — does not depend on any frequency. So the frequency only
matters to the extent that the *prior* still decides the genotype:

- **At 30 reads a position it decides nothing.** Thirty reads settle a genotype whatever the prior
  says, so shifting the prior by an eighth, or by two fifths, moves the likelihood by nothing that
  four decimal places can see. That is the whole of the 0.0263-against-0.0263 row.
- **At three reads the prior does decide** — and there the correction costs more in noise than it
  removes in bias. It divides the fitted value by `1 − leverage`, which inflates that value's own
  sampling error, and subtracts a dosage estimated from three reads, which is mostly noise itself. A
  noisy per-sample frequency manufactures contamination; this module already carries that mechanism,
  measured on panels split into groups, where the noise adds about 0.015 to *every* sample.

**The floor's behaviour is the signature.** The correction did not lift the contaminated sample
alone — the worst clean sample went 0.0003 → 0.0192 and the median 0.0000 → 0.0007. A bias being
removed lifts the contaminated sample; noise being added lifts everybody, and everybody is what rose.

That also settles what *leverage* has been measuring. On the unbalanced panel with nobody contaminated
at all, the groups' leverages were 0.027, 0.307, 0.429 and 0.857 and their spurious fractions rose
with them. That was read as self-influence biasing the estimate. It is better read as **noise**: a
sample supplying most of its own frequency has a frequency fitted from almost nothing else, and a
frequency fitted from almost nothing is noisy. Leave-one-out pushes every sample a little way along
that same axis, which is the wrong direction.

It is kept as `ContaminationConfig::leave_self_out`, defaulting to off, with the numbers on the field;
`LEAVE_SELF_OUT=1` runs the control from the other side. **What would confirm the account above**: the
same panel at three reads a position with 200 samples, where leverage is 0.025 rather than 0.125 — if
the story is right, the correction goes inert there too, because what it removes and what it adds both
scale with leverage. Not run; it would not change the switch either way.

---

## 5. On real alignments: the 63 tomato accessions

*Added the same day, after the rest of this report was written on drawn data alone. The accessions
were assumed to sit at three reads a position, where the change does nothing. **They do not** — the
cohort runs from 2.4× to 30.6×, so its deeper accessions had exactly the positions this change is
about.*

`examples/ng_joint_records_walk.rs` over the 63 benchmark CRAMs, 8 Mb of analysed regions, 1,999,404
kept positions, the same input as `joint_records_on_real_alignments_2026-08-13.md`. Raw output:
`tmp/contamination_depth/tomato_n63_exact_ladder.txt`.

**The encoding lands where it was priced, and nothing else moved.** The depth array is **1.999 MB at
8.00 bits a position**, against 1.25 MB before. The list of disagreeing reads is **89,241 entries for
the median accession — the same number the August 13th run recorded**, so what changed is the depth
encoding and not the walk. Every position below the cap now carries an exact depth; the two deepest
accessions hold **1 position in 1,000 above the cap**, which is the only place a range survives.

**The floor held and the top of the panel rose:**

| | the widening ladder (2026-08-13) | one bin per depth |
|---|---:|---:|
| median accession | 0.0000 | **0.0001** |
| highest accession | 0.0090 | **0.0398** |

**The median not moving is what says this is not a new floor** — a floor lifts everybody. Eight
accessions now read between 1% and 4%, and the highest supplies 0.024 of its own fitted frequency, so
it is not a sample reading its own echo. Beside the drawn control, where the exact ladder took a true
3% from 0.0120 to 0.0263 and left all 39 clean samples at 0.0000, the reading is that **a few of these
archive accessions are genuinely contaminated at a few percent and were being reported at less than
half of it**. There is no truth set on real reads; that is a reading and not a proof.

The estimate rests on 36,712 varying positions, one accession of the 63 was refused for supplying most
of its own frequency, and the rest of the fit is unremarkable: converged in 30 passes and 645 s, error
rate 0.00334 a base, 3.15% of positions judged mismapped, expected heterozygosity 4.155 per kilobase.

---

## 6. What this cannot say

- **One real cohort, and it is not deep.** The tomato run above reaches 30× in a handful of its
  accessions and sits far lower in most, so the depths where the change matters most are thinly
  sampled. A cohort that is deep throughout would say more, and this project does not have one.
- **A tenth of the value is still missing** at 30 reads: 0.0263 against 0.030. Neither the depth, nor
  the structure, nor the sample's own contribution to its frequency accounts for it, and no candidate
  here has been measured.
- **Only contamination was measured.** The same coarse depth is read by the heterozygosity and
  error-rate fits, where reading a band's midpoint rather than summing over it was already on record
  as costing 4 to 7% of heterozygosity at 30 reads. An exact ladder should help there too. Unmeasured.
