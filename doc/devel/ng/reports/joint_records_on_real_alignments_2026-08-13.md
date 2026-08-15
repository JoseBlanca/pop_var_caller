# The joint route's records, driven through real alignments — and what the cohort turned out to be

*2026-08-13. Everything measured for ng's second route to the calling parameters had, until now, been
measured against truths the measuring program made up itself. This is the first run of the real
writer — `parameter_estimation::joint::records::RecordWriter` — over real reads, on one human sample
and on sixty-three tomato accessions, and it answers three questions the specifications asserted
without evidence. Written for a reader who has read none of those specifications.*

*Programs: `examples/ng_joint_records_walk.rs`. Raw output: `tmp/records/tomato_cohort_n63.txt`,
`tmp/records/hg002_30x.txt`.*

---

## 1. What was run

Before a caller can genotype anything it has to know the things it will assume — how often the
sequencer misreads a base, how heterozygous each plant is, how diverse the population is. ng measures
those first, and it is building **two ways of doing it** so it can compare them. The one being built
here keeps raw evidence at a bounded set of positions — the same positions in every sample — and fits
everything once across the whole cohort.

What each sample writes down at those positions is a **record**, and the code that fills one had never
seen an alignment. Two runs now have:

| | reference | analysed stretch | kept positions | kept repeat tracts | samples |
|---|---|---|---|---:|---:|---:|
| tomato | SL4.00 | 8.0 Mb of the bench regions | 1,999,404 | 4,164 | 63 |
| GIAB HG002 at 30× | GRCh38 | 6.1 Mb over 50,000 tandem repeats | 1,999,148 | 20,204 | 1 |

Each walk fills the records **and**, in the same pass, the coverage-by-window summary the route needs
to tell a duplicated stretch of genome from a heterozygote. That summary existed only inside a
measuring program before this; it is now `parameter_estimation::joint::coverage`.

---

## 2. What the records weigh, measured rather than priced

The specification prices these objects by arithmetic and then says to measure them instead. Measured,
**two of the five rows are wrong, both in the same direction, and both because a cheap thing was
assumed to be cheaper than it is.**

### 2.1 The per-position depth array is exactly what was claimed

One five-bit code per kept position per read group: **1.250 MB at two million positions**, against
1.25 MB predicted. There is nothing to say about this row.

### 2.2 The list of disagreeing reads grows with depth, and the predicted range is a 3× range

Whenever a sample's reads at a kept position do not all show the reference base, the exceptions are
listed. The specification puts that list at **30–250 kB per sample**, saying it is driven by the
sequencer's error rate rather than by variants. **The mechanism is right and the range is a range for
one depth only** — the number of entries scales with depth, because a deeper sample gives error more
chances to show:

| sample | mean depth | entries | positions in 1,000 with an entry |
|---|---:|---:|---:|
| shallowest tomato accession | 2.4× | 9,213 | 5 |
| median tomato accession | 10.9× | 89,241 | 45 |
| deepest tomato accession | 30.6× | 331,036 | 135 |
| HG002 | 29.5× | 223,161 | 112 |

At four bytes an entry — the specification's own figure — that is **37 kB at 2.4× and 1.3 MB at
30.6×**. So the predicted range describes the shallow tomato archive, which is what it was written
for, and **understates a 30× sample by a factor of five**. Nothing needs redesigning; the row needs a
depth beside it.

*One number for whoever writes the byte format: an entry is 12 bytes in memory today (a position, an
allele, a count) against the 4 bytes the size claim assumes. The claim is reachable — a position needs
21 bits and an allele 3 — but only by packing, so it is a constraint on that format rather than a
property the code already has.*

### 2.3 The repeat-tract records cost 25 bytes a tract, not the 9 to 11 assumed

This is the row the specification itself flags as the one most likely to move, and it moved. Measured
identically on both genomes: **25.0 bytes per tract per read group**. The specification's figure was
4–5 MB across 462,701 tracts, which is 9–11 bytes each.

**What that does to the bill.** The tomato catalog holds 462,701 repeat tracts at the calling floors,
and the design keeps every one of them:

| | per sample | fifty samples |
|---|---:|---:|
| specification | 4–5 MB | 200–250 MB |
| measured | **11.6 MB** | **578 MB** |

The knob for it already exists and has never had a reason to fire: a cap on how many tracts each
class of repeat contributes. Nobody has to decide that today — the fit does not read the repeat
records yet — but the number to decide against is 578 MB and not 250 MB.

### 2.4 The list of mismatching bases inside a tract now exists, and it too scales with depth

This channel — which base of which read disagreed with the reference inside a repeat tract — was
specified and **was not being written at all**. It is written now, and it is what the repeat-tract
error rate has to be estimated from: a substitution that does not change a tract's length is invisible
to everything else the record holds.

Measured: **0.054 mismatching bases per tract at tomato's 2.4×, and 0.585 at HG002's 30×**. Scaled to
462,701 tracts that is 0.2 MB and 2.2 MB, against 0.3 MB predicted — right at three reads a site, ten
times low at thirty.

*One limitation, stated because it is deliberate.* A read whose tract has slipped a whole repeat unit
has no base-for-base correspondence with the reference, so this writer does not compare it: which of
its bases sits over which reference base is the aligner's answer, not the writer's. Such a read
contributes its length and **nothing to the denominator**, so the error rate stays a ratio of two
quantities counted over the same reads. At tomato roughly 95 reads in 100 are unslipped, so the
denominator loses about a twentieth of its support.

### 2.5 The coverage-by-window summary lands where it was priced

**1 byte per 500 bp window.** The runs cover only a benchmark's worth of genome, so the direct figures
are 16,066 windows on tomato and 61,759 on HG002; scaled to the whole reference that is **1.56 MB for
tomato and 6.2 MB for GRCh38**, against 1.6 MB and 6.2 MB predicted.

**The thing that had to be checked and was**: summing ten adjacent 500 bp windows gives exactly the
mean over the same 5 kb, weighting each window by how many of its positions are actually in scope. A
sum that treated every window as full agrees everywhere except at the end of a contig and the edge of
a region, which is precisely where it matters.

*One property of the stored byte, found while testing it: a window's depth is stored relative to the
sample's own median, in steps of 3%, reaching **eight times the median**. A window above that reads as
eight. The class this object exists to find sits at two, so the ceiling is not in the way — but a
tenfold pile-up is stored as eightfold.*

### 2.6 The whole cohort, held at once

**153 MB for 63 tomato accessions** at two million positions apiece, which is the state the estimator
runs in. That figure carries only the 4,164 repeat tracts inside the benchmark's 8 Mb; with the full
462,701 it becomes about 880 MB, against the specification's "roughly 360–420 MB for a fifty-sample
cohort".

---

## 3. A defect the walk exposed: two different things were being recorded as the same thing

The records are required to keep three states apart — a position the run **never opened** (which is a
bug), a position walked where **no read arrived** (which is data), and a position with reads. On the
first real walk, **93,150 of 1,999,404 kept positions on a 25× tomato accession came back as "never
opened"**, and every one of them was data.

The cause is not in the records. **The locus generator emits nothing at a position no read reached**,
so silence from the walk was indistinguishable from a region the run never opened. Fixed by having the
walk say which stretches it covered, whatever it found in them; a real depth is never overwritten by
that mark and the mark never overwrites a real depth. After the fix the same accession reports **0
never opened and 93,150 walked at zero depth**, which is the truth.

The share is not small and it varies more than depth alone explains: **1.4% to 35% of kept positions
across the 63 accessions**, and 0.5% on HG002. It is the denominator of every rate the route reports,
so a fifth of a shallow sample's heterozygosity estimate rests on positions where nothing was seen.

---

## 4. How far apart the 63 tomato accessions are — and it is the whole reason the contamination work matters

**This was unmeasured, and it decides how much of the contamination design is needed.**

Contamination is the fraction of a sample's reads that came from a different plant. It is identified
by finding reads carrying alleles the population does not expect — which means the estimator needs to
know what the population expects at each position. If the panel is genetically uniform, one pooled
answer serves every accession. If the panel is diverged, a pooled answer is wrong for everybody, and
the measured consequence is severe: at a divergence of `F_st` 0.20, **a genuinely 3%-contaminated
accession comes back at 0.5%** and passes as clean.

### What was measured

All 63 accessions' records were held at once — the state the estimator runs in — and 36,151 of the
two million kept positions turned out to segregate widely enough to say anything about who differs
from whom. Those positions were decomposed into axes of variation, the panel was split along the
leading axis, and the divergence across that split was measured.

**And the same pipeline was run on a control with no structure in it at all**: every sample's reads
redrawn at the position's own cohort-wide allele frequency and that sample's own depth, which destroys
structure and keeps the depths, the missing data and the frequencies exactly. Anything the control
returns is what the estimator manufactures on its own.

| | measured | control, no structure |
|---|---:|---:|
| share of variation on axis 1 | **22.4%** | 4.6% |
| share on axis 2 | 10.5% | 4.5% |
| divergence (`F_st`) across the axis-1 split | **0.44** | −0.01 |
| divergence across a split chosen at random | 0.003 | −0.01 |

**The panel divides 19 accessions against 44, and the two groups are as different from each other as
0.44.** Both controls agree that this is not an artefact: choosing the split from the data buys
nothing when there is nothing there (−0.01), and a split that ignores the structure on the *real* data
finds almost nothing (0.003).

**An internal check that the axis is genetic and not technical.** Six accessions appear twice in the
bench cohort, sequenced on separate runs. Every pair lands together on the axis — `SRS3394549` at
+0.114 and +0.112, `SRS3394688` at +0.109 and +0.108, `SRS3394712` at +0.062 and +0.070 — so the axis
is tracking the plant and not the run.

### What follows

**The tomato cohort sits at or past the worst row of the table this decision was to be read off.** The
pooled-frequency estimator, which is the simple one, would report a 3%-contaminated accession at
around half a percent on this panel. So the per-accession allele frequencies — fitting each position's
frequency as a straight line in the panel's own axes of variation, using every accession — are
**necessary rather than an improvement**, and the ceiling measurement that showed contamination is
recoverable when the frequency is right is the thing to build towards.

*One caveat about the size, not the direction.* `F_st` measured this way is biased upward in an inbred,
selfing panel because the correction it uses assumes random mating within each group. The random-split
figure of 0.003 is the honest floor including that bias, so the excess is real; whether the true value
is 0.44 or 0.35 does not change which row of the table the panel is in.

---

## 5. The estimator now exists for the ordinary-position half

`parameter_estimation::joint::fit` fits, in one place, from every sample's records:

- **two rates per read group** — how often a read misreads a base at an ordinary position, and at a
  mismapped one;
- **one cohort-wide share** of positions that are mismapped;
- **four numbers describing the population's allele frequencies** — how often a position carries only
  the reference base, how often only a non-reference one, and the shape of what is left;
- **one number per sample** saying how much less heterozygous it is than random mating would predict;
- and, derived from the converged answer rather than fitted, **each sample's heterozygosity and
  homozygous-non-reference rate**, with the count of positions that actually carried a read beside
  them.

### It recovers what it is given

A cohort of 10 samples at 8 reads a position, 3,000 positions, drawn at known values:

| | drawn | fitted |
|---|---:|---:|
| error rate at an ordinary position | 0.00200 | 0.00213 |
| error rate at a mismapped one | 0.0600 | 0.0534 |
| share of mismapped positions | 0.0200 | 0.0277 |
| share of positions carrying only the reference | 0.9000 | 0.9064 |
| share carrying only a non-reference base | 0.01000 | 0.00998 |
| homozygote excess | 0.200 | 0.186 |
| heterozygosity | 0.01700 | 0.01850 |

Twenty-seven passes, converged. **The two controls that say this measured something**: a cohort drawn
with no inbreeding at all comes back with none (every sample below 0.08 against a truth of zero,
started from a deliberately wrong shape), and a cohort of one sample still fits the error rate and
**marks the two parameters it cannot fit there as not fitted** rather than returning a plausible zero.

### One departure from the architecture, and it is why the program can be run

The architecture proposes a bounded search per parameter. **One evaluation of the likelihood is a pass
over two million positions times fifty samples times the frequency quadrature**, and a hundred
parameters searched that way is a thousand such passes — the program could not be run, which is
exactly the trap recorded in the previous session's notes. Instead there is **one pass per iteration**,
which accumulates the counts every parameter's own maximisation needs; each maximisation then runs
over a few dozen accumulated numbers with the reads untouched. Each iteration cannot lower the
likelihood, which a coordinate climb over this surface does not guarantee.

**A second departure, smaller:** the share of mismapped positions is one cohort-wide number, not one
per read group as the specification's table has it. A position is mismapped or it is not; it cannot be
mismapped for one library and clean for another at the same position. The per-read-group grain belongs
to the other route, which never sees more than one sample's libraries at a time.

### What is not in it

Contamination (it needs the per-accession frequencies of §4, which are specified and unbuilt), the
duplicated-stretch class (it needs the coverage summary wired into the likelihood), the repeat-tract
half, and any ploidy but two.

---

## 6. What I would do next, in order

1. **The comparison the route exists for** — both routes on the GIAB trio at 30× against its truth
   VCFs, a sweep in sample count, and the 63 accessions. The estimator now exists for the
   ordinary-position half, which was the thing blocking it.
2. **The per-accession allele frequencies**, because §4 says the tomato panel cannot be served without
   them.
3. **The repeat-tract half of the estimator**, whose model was measured and held.
4. **The byte-level record format**, where the 12-bytes-against-4 of §2.2 is the constraint.
