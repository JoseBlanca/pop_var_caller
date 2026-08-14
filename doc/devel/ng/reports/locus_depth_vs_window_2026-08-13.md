# The position's own read count cannot replace the coverage summary: at three reads a position it knows nothing, and calls 44 in every 100 positions it scores duplicated

*Research report, 2026-08-13. **Recommendation: keep the coverage-by-window summary.** One
program stands behind this, `examples/ng_joint_duplication_probe.rs`, extended with a second
and third discriminator; eight tomato accessions from 2.5× to 28.7× mean depth, two window
widths each. Raw output in `tmp/locusdepth/`. The window arm's numbers reproduce
[`duplicated_locus_probe_2026-08-12.md`](duplicated_locus_probe_2026-08-12.md) to the digit,
which is what says the two sets of numbers can be compared.*

---

## 1. The question, and the answer

Before calling anything, the caller estimates how heterozygous each sample is. **A stretch of
genome a plant carries twice while the reference carries it once corrupts that estimate**:
both copies' reads pile onto the one place the reference offers, and wherever the two copies
differ the position reads about half non-reference — which is what a heterozygote looks like.

Two ways to recognise such a position were on the table, and a subsystem hangs on which is
needed:

- **the coverage-by-window summary** — the sample's mean depth in every 500 bp window,
  divided by what that window's GC content predicts. Expensive: a genome-wide depth
  accumulator, a per-sample GC curve, denominators read off the reference, and a second full
  pass over every sample's stored pileup when the estimate runs from stored data;
- **the position's own read count** — already in every kept record, five bits, free.

**The read count is not a substitute, and the gap is largest exactly where the caller
lives.** Three findings, in the order that matters:

1. **At the depth the tomato cohort actually runs — 2.5 to 5.2 reads a position — the read
   count carries nothing.** Its enrichment is 0.63 to 1.14, and 1.00 is the value a
   discriminator unrelated to the artefact returns — one of the four shallow accessions comes
   out below it and one sits exactly on it. The window at 5 kb gives 6.6 to 14.4 on the same
   positions.
2. **It does not merely fail to enrich, it floods the class.** At 2.5 reads a position it
   calls **44 of every 100 positions it scores** "about two copies". The window calls 2.8.
3. **Where it does work it is half the window, and it stops working at 98 reads a
   position.** At 25.2× it reaches 10.8 against the window's 21.3, flagging 3.5 times as many
   positions to get there. From 76 reads a position the stored five-bit code no longer puts
   a doubled position 1.6 times above an ordinary one, and from 98 it gives the two the
   identical code. That ceiling is a property of the ladder, not of any sample — it is the
   same 98 for every run.

**What the read count *is* good for is adding to the window, from about ten reads a position
up.** Requiring both raises the enrichment from 9.6 to 16.7 at 9.9×, and from 21.3 to 25.7 at
25.2×. Below ten reads it makes things worse. §7 has the proposal and why it is not free.

---

## 2. What is being counted, and what it is not

**Near half** means a position whose reads disagree with the reference in a fraction between
0.35 and 0.65, with at least two disagreeing reads so that one stray read at low depth is not
a half. That is the artefact's signature.

**Enrichment** is how many positions a discriminator flags *and* that read near half, divided
by what the two happening independently would give. An enrichment of 1.0 means the
discriminator knows nothing; 20 means the flagged positions read near half twenty times as
often as the flagging rate and the near-half rate together predict.

**There is no truth set.** Nobody has a validated list of the stretches these accessions
carry twice, so nothing here is a detection rate and nothing here says a flagged position
*is* a duplication. Enrichment is a comparison between discriminators on the same positions,
and that is all it is. It is usable for this decision precisely because all three arms are
scored on the same positions and judged against the same near-half target.

**Every arm scores the same positions**: those at four reads or more whose window has a
relative coverage at all. Two runs per accession, one with 500 bp windows — the width the
records specify — and one with 5 kb, the width the earlier report found the shallow samples
need.

---

## 3. The three arms, per accession

Enrichment, positions at four reads or more. The read-count arm does not depend on the window
width; its column is from the 5 kb run.

| accession | mean depth | window, 500 bp | window, 5 kb | the position's own read count | both required |
|---|---:|---:|---:|---:|---:|
| SRR7279533 | 2.51× | 1.6× | **14.0×** | **1.14×** | 14.4× |
| SRR7279488 | 2.72× | 1.6× | **14.4×** | **1.01×** | 13.0× |
| SRR7279501 | 3.60× | 1.3× | **6.6×** | **0.63×** | 5.5× |
| SRR7279484 | 5.15× | 1.5× | **7.9×** | **1.00×** | 9.8× |
| SRR7279481 | 9.89× | 2.5× | 9.6× | 1.65× | **16.7×** |
| SRR7279483 | 13.32× | 7.7× | 16.0× | 3.21× | **20.5×** |
| SRR7279482 | 25.20× | 24.0× | 21.3× | 10.8× | **25.7×** |
| SRR7279540 | 28.69× | 24.9× | 23.4× | 10.8× | 24.0× |

**The read-count arm is flat at 1 across the whole bottom half of the cohort.** It reaches a
fifth of the window's enrichment at 13× and half of it at 25×. The four shallow accessions
are the ones this caller is aimed at.

**How many positions each arm has to flag to get that** — the share of scored positions it
calls "about two copies", which is the class the fit would then have to carry:

| accession | mean depth | window, 5 kb | the position's own read count |
|---|---:|---:|---:|
| SRR7279533 | 2.51× | 2.8% | **44.3%** |
| SRR7279488 | 2.72× | 3.2% | **46.5%** |
| SRR7279501 | 3.60× | 5.9% | **32.6%** |
| SRR7279484 | 5.15× | 3.0% | **19.5%** |
| SRR7279481 | 9.89× | 2.7% | 12.3% |
| SRR7279483 | 13.32× | 1.4% | 8.7% |
| SRR7279482 | 25.20× | 0.66% | 2.3% |
| SRR7279540 | 28.69× | 0.72% | 2.8% |

**At three reads a position the read-count arm puts nearly half the positions it scores in the
duplicated class**, and that is not a threshold wanting tuning. A position is only scored when
it has four reads or more, which at a mean of 2.8 is 31 positions in 100 to begin with; among
those, five and six reads are the two commonest values there are, and "about twice the median"
is exactly five and six reads. So the band that means *doubled* is the band that holds the
ordinary positions. Raising its lower edge trades the flood for an empty class: there is no
setting at which one Poisson draw at a mean of three separates one copy from two.

---

## 4. What the record's own encoding costs

The record does not hold a read count. It holds **a five-bit code standing for a range of
counts** — exact for zero to eight reads, then eleven widening rungs to a ceiling of 124 —
and the count is thinned to that ceiling before anything is written. So the arm was scored
twice: once reading the code as the fit would (the middle of the range it stands for), and
once reading the exact count the walk saw, which no stored record can offer.

| accession | mean depth | as the record stores it | as the walk sees it | what the encoding costs |
|---|---:|---:|---:|---:|
| SRR7279533 | 2.51× | 1.14× | 1.14× | nothing |
| SRR7279488 | 2.72× | 1.01× | 1.01× | nothing |
| SRR7279501 | 3.60× | 0.63× | 0.66× | 5% |
| SRR7279484 | 5.15× | 1.00× | 1.04× | 4% |
| SRR7279481 | 9.89× | 1.65× | 1.75× | 6% |
| SRR7279483 | 13.32× | 3.21× | 3.61× | 11% |
| SRR7279482 | 25.20× | 10.8× | 12.5× | 14% |
| SRR7279540 | 28.69× | 10.8× | 17.1× | **37%** |

**Below nine reads a position the ladder stores the count exactly, so the encoding costs
nothing — and there the arm has nothing to lose.** The cost arrives with the widening rungs,
and at 28.7× it removes more than a third of what the arm had. **This is the first
measurement of what the ladder costs a per-position discriminator**; the 0.054-rung figure in
the generic path's specification was measured on pooled histogram cells, where many positions
share a bin, and it does not describe a single record read alone.

---

## 5. Where the cap ends the arm, and it is inside this caller's range

The ceiling of 124 is not a tuning knob for this arm — it is the top of the ladder, and above
it a doubled position and an ordinary one are written identically. Two thresholds follow from
the ladder alone, so they are the same number for every sample:

- **from 76 reads a position**, a doubled position's stored code no longer reads 1.6 times an
  ordinary one's, which is the lower edge of the two-copy band — so the arm stops flagging
  duplications at all;
- **from 98 reads a position**, a doubled position is written as the very same code as an
  ordinary one. The arm is then blind, not merely quiet.

**On tomato the cap does not fire**: the deepest accession's median position holds 31 reads,
and the deepest single position any accession reaches is 234. Two accessions have positions
over the ceiling — 9,420 of 7.4 M on SRR7279483 and 1,799 of 7.5 M on SRR7279540, one in a
thousand and fewer. Those are single positions far above their sample's own level, not the
sample's depth reaching the ceiling.

**But 76 reads a position is inside the range this caller commits to.** A whole-genome human
run at 30× sits at a third of it; a targeted or exome run at 100× to 200× is past it
entirely. So the cheap arm's answer is not "works at high depth": it works over a band from
about ten to about seventy reads a position, and is blind on either side. The window has no
such ceiling — a doubled window reads 2.0 at any depth.

Two smaller facts, both checked rather than assumed:

- **The cap does not damage the near-half signal.** Above it the alternative reads are thinned
  by the same ratio as the depth, so the fraction survives. Across all eight accessions
  exactly **one position in 59 million** has its near-half status changed by that rounding.
- **The cap and the ladder act per read group**, and these accessions have one read group
  each, so a position's depth here is a read group's depth. A sample sequenced as several
  libraries reaches the ceiling later in total reads, and at the same point per library.

---

## 6. What the GC correction costs the cheap arm: nothing genome-wide

Depth tracks GC content by more than the doubling being looked for. Over single positions on
SRR7279482 the mean depth runs from **11.7 reads a position at 18% GC to 29.0 at 34%**, a
factor of 2.5 — larger than the two-against-one the arm is trying to see, and larger than the
1.79-fold swing the same sample shows across window means, because a window average dilutes
the extremes. So a read-count arm needs the same correction the window arm has, or the
comparison is rigged. **It does not need the window's machinery to get it**, and that is two
separate points:

- **The GC content itself comes off the reference**, not off any depth accumulator. Reading
  the GC of the 500 bp around a kept position costs one reference fetch per kept position.
- **The depth-against-GC curve can be fitted from the kept positions alone.** Fitting it from
  one walked position in 300 — about 25,000 positions a sample, the order of what the records
  keep — gives the same answer as fitting it from all 7.5 M: 10.72 against 10.77 at 25.2×,
  3.14 against 3.21 at 13.3×, 1.61 against 1.65 at 9.9×.

**Skipping the correction is what would cost**, and unevenly: dropping it takes 1.65 to 0.81
at 9.9× and 1.00 to 0.70 at 5.2×, but *raises* 10.8 to 13.9 at 28.7× and 3.21 to 3.82 at
13.3×. The window arm behaves the same way at high depth — the earlier report has it at 24.8
corrected against 32.6 uncorrected on SRR7279482. So on this data the correction removes
false positives at low depth and removes some true ones at high depth. **It is needed for the
comparison between arms to be fair, and it is not established that the fit wants it**; that
question is open for both arms and unchanged by this report.

**One thing this changes on the cost side of the window.** The summary as specified stores
the sample's depth-against-GC curve and each window's mean depth, and not the window's own GC
fraction — so nothing downstream can divide one by the other. The fix is one byte a window
and is now in the specification, but two things follow: the expensive arm is slightly more
expensive than it was priced at, and **no code has yet run the summary end to end**, since
anything that had would have hit this. The window arm's numbers here and in the earlier
report come from a correction the probe computes for itself during the walk. They are what
the summary *could* deliver once it stores GC, not what it delivers today.

---

## 7. Do the two add to each other? Above ten reads a position, yes

Enrichment says how good each discriminator is alone. It cannot say whether one carries
anything the other does not. The near-half rate in each cell of *window band × position band*
can, and it needs no threshold argument: read across a row to see what the position's count
adds at a fixed window verdict, and down a column to see what the window adds at a fixed
position verdict.

Positions at four reads or more, 5 kb windows. Each cell is the near-half rate and the number
of positions it is over.

**SRR7279533, 2.51× — the position's count adds nothing**

| window \ position | one copy | two copies |
|---|---|---|
| **one copy** | 0.0104% of 708,577 | 0.0190% of 817,891 |
| **two copies** | 0.3801% of 8,682 | 0.5281% of 29,160 |

Down the column the rate rises 37-fold (0.0104% → 0.3801%): the window turns positions into
near-half positions. Across the row it moves 1.8-fold, over 818,000 positions the read count
has called doubled and that read no more near-half than the rest.

**SRR7279540, 28.69× — both add**

| window \ position | one copy | two copies |
|---|---|---|
| **one copy** | 0.0117% of 6,296,517 | 0.3349% of 165,424 |
| **two copies** | 0.7118% of 8,008 | 1.1874% of 28,887 |

Here the row rises 29-fold and the column 61-fold. The two are finding partly different
positions, which is why requiring both beats either.

**Across all eight accessions**, how much each discriminator multiplies the near-half rate
once the other is held fixed. Above 1 means it is carrying something the other does not.

| accession | mean depth | read count, inside one-copy windows | read count, inside two-copy windows | window, inside one-copy read counts |
|---|---:|---:|---:|---:|
| SRR7279533 | 2.51× | 1.8× | 1.4× | 37× |
| SRR7279488 | 2.72× | 2.0× | 1.3× | 48× |
| SRR7279501 | 3.60× | 0.4× | 0.9× | 8.7× |
| SRR7279484 | 5.15× | 0.8× | 2.9× | 5.2× |
| SRR7279481 | 9.89× | 1.1× | **4.4×** | 8.7× |
| SRR7279483 | 13.32× | 2.9× | 1.8× | 32× |
| SRR7279482 | 25.20× | **12×** | **5.9×** | 8.4× |
| SRR7279540 | 28.69× | **29×** | 1.7× | 61× |

**The window carries something the read count does not at every depth** — 5.2-fold at the
worst accession and 61-fold at the best. **The read count carries something the window does
not only from about ten reads a position**, and at the shallow end its two columns straddle 1
in both directions, which is what noise looks like.

**What I would do with that, and it is not free.** A rule that requires both discriminators
above some depth and only the window below it is a depth-conditional knob in a fit that has
none. It buys 9.6 → 16.7 at 9.9× and 21.3 → 25.7 at 25.2×, against a loss of 14.4 → 13.0 and
6.6 → 5.5 at the shallow end. **I would not add it now.** Worth measuring inside the fit
first, as a change in fitted heterozygosity on a drawn panel — that is the only measure which
says whether a difference in enrichment reaches anything the caller emits.

---

## 8. What this cannot say

- **No truth set**, so every number here is a comparison between discriminators, never an
  accuracy. §2 says this and it bears repeating: an arm at 14× enrichment is not detecting 14
  times more duplications, it is flagging positions that read near half 14 times more often
  than chance would put them together.
- **8 Mb of tomato chromosomes 1 to 12**, randomly placed, one species, one aligner. Unplaced
  contigs are excluded and those are where an assembly's collapsed copies concentrate.
- **The near-half test needs at least two disagreeing reads**, so a position with more reads
  is mechanically a little likelier to pass it. That confound works *in the read-count arm's
  favour* — it is strongest at low depth, where the arm's flagged positions are exactly the
  deeper ones — and the arm still returns 1.0 there. At 25× and above, where the arm does
  score, both bands clear the two-read floor easily and the confound is small. So it does not
  explain the finding in either direction.
- **The read-count arm is scored on every generic position, not on the loci the records
  actually keep.** Keeping one position in a few hundred changes how many positions the fit
  has, not how well a single position's depth discriminates, so the enrichment carries over;
  the class's *size* on kept loci does not follow from anything here.
- **The window arm's numbers assume the summary can apply its own GC correction**, which as
  of today it cannot (§6).

---

## 9. What I would change in the specifications — not made here

These are proposals; nothing under `doc/devel/ng/spec/` or `doc/devel/ng/arch/` was edited.

| document | proposed change |
|---|---|
| `spec/parameter_prepass_joint_records.md` §4 | Record that the per-position read count was measured as an alternative to the summary and rejected, with the numbers of §3 — otherwise the same question returns every time someone notices the five bits are already there. State the window summary's justification in terms of the shallow end (14× at 2.5 reads a position against 1.1×), which is where it is decided, and not the 25× case. |
| `spec/parameter_prepass_joint_fit.md` §2.2 | The class's discriminator stays local relative coverage and never the position's own depth, which the section already says; add *why*, in one line with the 44-in-100 figure, since the reason is not obvious and the cheap alternative looks free. |
| `spec/parameter_prepass_generic.md` §4 (the depth ladder) | The ladder's cost is stated as 0.054 rungs of the error-rate ladder, measured on pooled histogram cells. Note that a *per-record* consumer pays differently — 14% and 37% of a discriminator's enrichment at 25× and 28.7× (§4) — and that reading a single stored code is a different use from summing many into a histogram. |
| `spec/parameter_prepass_joint_records.md` §4 | If the summary is ever to carry a per-position companion, the ceiling of 124 must be stated as a limit of that companion and not only of the histogram: a discriminator built on the stored code is blind from 98 reads a position (§5), which a 100× run reaches. |
