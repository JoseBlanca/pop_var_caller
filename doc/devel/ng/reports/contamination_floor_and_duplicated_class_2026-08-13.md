# The contamination floor on real reads had two causes, and the second one is the depth ladder

*Research report, 2026-08-13, later the same day as
[`joint_fit_against_truth_2026-08-13.md`](joint_fit_against_truth_2026-08-13.md), which reported the
floor and named one of the two causes. Written for a reader who has read none of the specifications.*

*Programs: `examples/ng_joint_contamination_control.rs` (the drawn control),
`examples/ng_joint_records_walk.rs` (the 63 tomato accessions),
`examples/ng_joint_duplicated_in_fit.rs` (the third class inside the estimator). Raw output under
`tmp/records/`.*

---

## 1. What was wrong, and what it is now

ng measures the numbers a caller will assume before any variant is called. One of them is
**contamination** — the share of a sample's reads that came from a different plant, a second seedling
in the tube or a neighbouring library on the same sequencing run. It is identified by *a small share
of reads carrying an allele the sample should not have*.

Run over 63 tomato accessions from a public archive, it came back with a **median accession 6.5%
contaminated and the highest 12.5%**. Sixty-three archive accessions are not all one part in fifteen
somebody else's plant. That was a floor and not a measurement, and this report is what was under it.

**Two causes, both now fixed, and the second had not been named by anyone.** On the same 63
accessions and the same reads:

| | median accession | highest |
|---|---:|---:|
| as reported this morning | 0.0684 | 0.1246 |
| positions that mismap excluded (§2) | 0.0091 | 0.0680 |
| …and the read depth read as the range it is, not as one number (§4) | **0.0000** | **0.0090** |

**The median accession now returns exactly zero and the worst 0.9%**, and a drawn panel with one
sample genuinely 3% contaminated still finds it at thirty times the worst clean sample (§3). The
first cause is the one the morning's report named. On its own it takes the median from 0.0684 to
0.0091; the depth ladder on its own takes it to 0.0300. **Neither alone reaches zero and both
together do.**

---

## 2. The first cause: a position that mismaps looks contaminated in everybody

**A position where two stretches of genome the reference holds once both pile reads up** — one
sequence mapped on top of another — has a few reads disagreeing with the reference **in every sample
at once**. That is the contamination signature exactly, and nothing excluded such positions from the
ones contamination was measured over.

Worse, those positions are *enriched* among the ones contamination uses. Contamination can only be
seen where the cohort varies, so the marker list is built from positions where some samples' reads
disagree — and a position that mismaps disagrees in everybody, so it looks like a position the
cohort varies at. **On tomato, 20,767 of the 52,525 markers were positions the fit says are
mismapped: two markers in five.**

### What was built

The fit already knows. Holding every sample at the same position is the whole reason this route
exists, and what it buys is exactly this: with one sample a few disagreeing reads are
indistinguishable from a rare heterozygote, and with sixty-three, a position that reads part
non-reference *in everybody* has nowhere else to go. The estimator computed that probability for
every position on every pass and then threw it away.

It is now kept — one four-byte number per position, 8 MB at the two-million-position budget — and
contamination drops every position more likely mismapped than not, weighting what survives by its own
probability of being ordinary.

---

## 3. The drawn control: the floor goes and the signal stays

The trap this project has fallen into before is that **an estimator with no power looks exactly like
a clean answer**. So the fix is measured on a drawn panel in three states, all at three reads a
position: 40 samples in 4 subpopulations diverged at `F_st` 0.20, 400,000 positions, with 1 position
in 30 planted as mismapped at a 2.4% disagreement rate — tomato's own numbers. One sample is drawn
genuinely 3% contaminated by a plant taken from the whole panel.

| the panel | sample 0, drawn at 3% | median of the other 39 | worst of them |
|---|---:|---:|---:|
| mismapped positions planted, **every position used** — this morning's estimator | 0.0251 | 0.0081 | 0.0144 |
| mismapped positions planted, **those positions dropped** | 0.0129 | 0.0011 | 0.0057 |
| mismapped positions planted, **dropped, and the depth summed over** (§4) | **0.0102** | **0.0000** | **0.0003** |
| **no mismapped positions at all**, dropped and depth summed over | 0.0115 | 0.0000 | 0.0008 |

And the same panel with **nobody contaminated at all**:

| | sample 0 | median | worst |
|---|---:|---:|---:|
| every position used | 0.0091 | 0.0086 | 0.0124 |
| those positions dropped | 0.0022 | 0.0016 | 0.0042 |
| dropped, and the depth summed over | **0.0000** | **0.0000** | **0.0014** |

**The floor is gone and the sample is still found.** Nobody contaminated returns zero for the median
accession and 0.0014 for the worst; the genuinely contaminated one returns 0.0102, which is **thirty
times the worst clean sample**. And the last row of the first table is the control that says the
exclusion is not simply blindness: handed a panel with no mismapped positions to exclude, the
estimator returns 0.0115 for the same sample — so excluding them costs the signal about a tenth of
itself and the floor all of itself.

The fit's own recovery of the class is what makes it work: it books 0.0315 of positions as mismapped
against 0.033 planted, and condemns 8,729 of the 400,000 individually.

---

## 4. The second cause, which nobody had named: the depth ladder manufactures contamination

Run the same drawn panel with **nobody contaminated and no mismapped positions at all** — nothing
whatever to find — and sweep the read depth:

| reads a position | median accession's fitted contamination | worst |
|---:|---:|---:|
| 3 | **0.0013** | 0.0061 |
| 10 | **0.0248** | 0.0283 |
| 30 | **0.0222** | 0.0268 |

**A panel with nothing in it returns two and a half percent at ten reads a position and one part in
eight hundred at three.** Nothing changed but the depth.

**It is not population structure.** The same sweep at `F_st` 0 and at 0.20 differs by less than a
fifth of a percentage point, and turning the structure model off entirely — one allele frequency for
the whole panel instead of one per sample — moves the median by 0.0001.

**It is how the depth is stored.** A position's read count is kept as one of twenty five-bit codes.
Below nine reads each code is one exact depth; above nine a code stands for a range that widens as it
climbs. **The count of disagreeing reads is exact and the depth is not**, so a heterozygote with 33
reads whose code says *26 to 35* is scored as having about 30, and its read share lands at 0.55
rather than 0.50 for a reason that has nothing to do with the sample. **A read share away from a half
is what a contamination fraction is made of**, so the estimator books the difference as
contamination — in every sample at once, which is why it reads as a floor.

The pattern matches the ladder exactly: three reads a position sits inside the exact region and
returns 0.0013; ten and thirty sit outside it and return 0.025 and 0.022.

### The fix is one sum, and it takes the floor to zero

Each sample's reads are now scored against **every depth its code could stand for**, weighted
equally, instead of against the middle of the range. Below nine reads the range is one value and this
is the plain binomial it always was. On the same drawn panel at ten reads a position, with nobody
contaminated and nothing planted:

| | median accession | worst |
|---|---:|---:|
| depth read as the middle of its range | 0.0405 | 0.0435 |
| **depth summed over its range** | **0.0000** | **0.0000** |

Every one of the forty samples returns exactly zero. It costs the width of the code's range in
arithmetic — four seconds over the whole tomato panel.

---

## 5. On the 63 tomato accessions

Two runs over the same reads, the same 1,999,404 kept positions and the same converged fit — 29
passes, a read misreading at 0.00333 at an ordinary position and 0.0239 at a mismapped one, 1 position
in 30 in the mismapped class.

| | markers | median accession | highest |
|---|---:|---:|---:|
| as reported this morning | 52,525 | 0.0684 | 0.1246 |
| depth summed over its range, mismapped positions kept | 52,525 | 0.0300 | 0.0905 |
| mismapped positions dropped, depth read as one number | 31,758 | 0.0091 | 0.0680 |
| **both** | 31,758 | **0.0000** | **0.0090** |

**Two markers in five were positions that mismap**: 20,767 of 52,525. Over the whole census the fit
condemns 62,728 of 1,999,404 positions — 1 in 32 — so the mismapped positions are enriched about
thirteen-fold among the markers, which is what §2 predicts and what makes them worth this much.

**Each fix on its own halves the median and the two together take it to zero**, which is what two
independent causes look like. Fifty-nine of the 62 accessions estimated now return exactly zero; the
three that do not are `SRS3394709` at 0.0090, `SRS3394633` at 0.0041 and `SRS3394663` at 0.0023.
`SRS3394702` is still refused, and for the reason it was before: it sits furthest out on the panel's
leading axis and supplies most of its own fitted allele frequency.

**`SRS3394709` is the one to look at, and it is not clear-cut.** It supplies 0.46 of its own fitted
frequency, just under the half at which a sample is refused outright, so part of its 0.9% is the same
effect that refuses `SRS3394702`. Against that, it is now ten times the median accession where this
morning it was 1.1 times it. **What the panel says is that it stands out; what it does not yet say is
by how much**, for the reason §6 gives.

---

## 6. What the three readings of a sample's own coordinates are worth

`verifyBamID2` maximises over the fraction **and the contaminated sample's own place on the panel's
axes of variation together**, and this route did not. The reason it matters: a stray read comes from
whoever else was on the plate, whose expected genotype is the panel average — which is the origin of
those axes — so a fraction `α` of stray reads drags the sample a fraction `α` of the way to the
origin, and the frequency the fit then predicts for it sits closer to the contaminant's than the
truth. The difference between the two is the entire signal.

Three readings are now built and measured. **As read**, unchanged. **Divided by `1 − α`**, which
undoes exactly that drag with no freedom of its own — at the true `α` the sample is put back where it
would have stood uncontaminated, and at zero nothing moves. And **each axis searched freely beside
`α`**, the literal reading of *maximise over both*.

On the drawn panel with a sample at 3%, no mismapped positions, depth summed over:

| where the sample is taken to stand | that sample | worst of the 39 clean | separation |
|---|---:|---:|---:|
| **as read** | 0.0115 | 0.0008 | **14×** |
| divided by `1 − α` | 0.0166 | 0.0046 | 3.6× |
| each axis free | 0.0190 | 0.0146 | 1.3× |

**The magnitude and the separation move together, so correcting the coordinates buys nothing.**
Every reading that gets the value nearer the truth of 0.030 raises the clean samples by about as
much. The reason is visible in the mechanism: at three reads a position a sample's genotypes are
mostly prior, so its coordinates are pulled towards the origin by the *prior* as well as by any
contamination, and nothing distinguishes the two — inflating them helps every sample's fit, not only
a contaminated one's.

**Recommendation: keep the coordinates as read, which is what ships.** What contamination has to do
is flag a sample, and 14× separation does that where 1.3× does not. **The bias in the value is
therefore still open**, and `α` should still be read as *this sample stands out from the panel* and
not as *this sample is 1.2% contaminated*. What has changed is that the panel it stands out from now
sits at zero rather than at 6.8%.

---

## 7. The duplicated-stretch class is now in the estimator, and it behaves

Where a plant carries **two copies of a stretch the reference holds once**, both copies' reads land
at the same place and, wherever the copies differ from each other, about half of them disagree with
the reference. That is what a heterozygote looks like, and until today the estimator had nowhere else
to put it.

[`duplicated_class_identification_2026-08-13.md`](duplicated_class_identification_2026-08-13.md)
measured what that costs and what identifies it, on a program written for the question. **This is the
same cohort shape put through the estimator itself**, which is a different claim: that the class
works where it has to live, beside the two noise classes, the frequency density and the homozygote
excess, all fitted at once.

Fifty samples, 60,000 positions, three reads a position, a selfing panel whose plants are 0.6 less
heterozygous than random mating predicts, one duplicated carrier position for every three genuinely
heterozygous ones:

| | heterozygosity | against the truth | homozygote excess | truth |
|---|---:|---:|---:|---:|
| the class off | 0.00100 | **+60.8%** | **0.4209** | 0.600 |
| identified from the cohort | 0.00062 | −1.2% | **0.5948** | 0.600 |
| …and from each sample's coverage too | 0.00063 | +0.9% | 0.5927 | 0.600 |

**And the control that decides whether it can be on by default: the same panel with no duplicated
positions in it at all.**

| samples | the class off | identified from the cohort | class weight it invents |
|---:|---:|---:|---:|
| 10 | −1.1% | −6.4% | 0.00098 |
| 25 | −0.4% | −1.6% | 0.00011 |
| 50 | −0.3% | −0.5% | 0.00003 |

**At twenty-five samples and above the class costs nothing when there is nothing to find** — it
shrinks its own weight to 3 positions in 100,000 and moves heterozygosity by half a percent. **At ten
it takes heterozygosity 6.4% low**, which is the same sample-count floor the other report measured
from the other side. Handed coverage readings as well, it returns exactly zero weight at every panel
size.

**So it ships on**, and a run of fewer than about twenty-five samples with no coverage summary should
turn it off rather than pay 6% of its heterozygosity for a class it cannot identify.

**One cost, and it is convergence.** A cohort with no duplicated positions has to shrink the class's
weight to nothing, and a weight that halves on every pass never stops moving by much *relative to
itself*, so the run reports it has not settled. Shares are now judged against a floor of one position
in ten thousand rather than against themselves, which fixes the reported state; the arithmetic still
takes more passes — several arms of the table above ran out at 120 where the same fit without the
class settles in 20.

---

## 8. The trio's heterozygosity excess is not the duplicated class

The benchmark trio's heterozygosity comes back **1.23 to 1.28 times its benchmark VCFs'**, and where
that comes from is open. The class above was the leading candidate: three samples is far below the
twenty-five the cohort pattern needs, but the trio's coverage summary exists and a coverage reading
works at any panel size — so the trio is the one place where a real truth set could say whether the
class is real rather than only arithmetically consistent.

**It is not.** With the class fitted and each sample's coverage reading supplied at every position:

| | heterozygosity, per kilobase | benchmark VCF | before |
|---|---:|---:|---:|
| HG002 | 0.806 | 0.639 | 0.806 |
| HG003 | 0.761 | 0.596 | 0.761 |
| HG004 | 0.806 | 0.654 | 0.806 |

**Unchanged to three decimal places, and the class fits a weight of exactly zero**: of 449,489
positions the fit places none in it. So whatever produces the excess, it is not a stretch these three
people carry twice.

**Two reasons to believe that rather than suspect the test, and one caveat.** The regions are GIAB's
own high-confidence set, which is constructed to exclude the collapsed and duplicated parts of the
genome — so a class of duplications is exactly what should be absent there. And the coverage reading
does see something: 1,687 of 449,489 positions read more like two copies than one in HG002, about 1
in 266, and the fit still declines to build a class out of them. **The caveat is that the coverage
reading has no GC correction** (§9, change 6), so it is copy number plus GC content rather than copy
number; on this benchmark at 330 reads a position that matters less than it would on tomato, but it
is not nothing.

**Where I would look next for the 1.26.** The fit judges 1.80 positions per kilobase to be
segregating within the trio where the truth is 1.219, and it puts 1 position in 108 in the mismapped
class at a 4.3% disagreement rate. Three samples is six chromosomes to estimate a population
frequency density from, and the fit still runs out of passes there — 200 without settling, where the
63-accession cohort settles in 29. **A frequency density that has not settled, fitted from six
chromosomes, is the remaining candidate**, and it is testable by fitting the trio inside a larger
human cohort rather than alone.

---

## 9. What I would change in the specifications

*Nothing under `spec/` or `arch/` was edited. These are the changes I would make, in the order I
would make them.*

**In [`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md):**

1. **§3.4, contamination — add the two exclusions, because neither is a tuning detail.** A position
   the fit judges more likely mismapped than not is not a marker, and a sample's reads are scored
   against every depth its stored code could stand for rather than against the middle of the range.
   Without the first the median tomato accession reads 6.8% contaminated; without the second a drawn
   panel with nothing in it reads 2.5% at ten reads a position and zero at three. Both belong in the
   specification as requirements with those numbers beside them.
2. **§3.4, the per-position probability of being mismapped is now an output of the fit, not an
   internal.** It is what having every sample at one position buys, and it has two consumers rather
   than one: contamination, and whatever calls variants afterwards. Four bytes a position.
3. **§3.4.4's noise floor of about 1% at 10,000 markers is superseded.** With those two exclusions the
   drawn floor is **zero** for the median accession and 0.0014 for the worst, at 8,879 markers — so
   the recommendation to raise the marker budget for the sake of the floor no longer has a floor to
   raise it against. What remains true is that the *value* is attenuated, so judging a sample against
   the panel's own spread stays right.
4. **§3.4, the joint maximisation over the sample's own coordinates — record it as measured and not
   adopted.** Three readings are built. Undoing the drag with `1 − α` moves a drawn 3% sample from
   0.0115 to 0.0166 against a truth of 0.030 and moves the worst clean sample from 0.0008 to 0.0046,
   so the separation falls from 14× to 3.6×; searching each axis freely is worse again. The attenuation
   and the floor move together and correcting the coordinates does not separate them.
5. **§2.2 and §5.1, the duplicated class — add what turning it off is worth and what it costs when
   there is nothing to find.** Measured inside the estimator: heterozygosity 60.8% high and the
   homozygote excess 0.4209 against a drawn 0.600 with the class off; −1.2% and 0.5948 with it on.
   And at ten samples with no duplications at all the class takes heterozygosity 6.4% low, which is
   the number a run under about twenty-five samples needs in order to decide to turn it off.

**In [`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md):**

6. **§4 — the coverage summary must keep each window's GC content, or it cannot be used as a
   copy-number discriminator at all.** It stores the sample's depth-against-GC curve and the window's
   mean depth, and not the window's GC fraction, so nothing downstream can divide one by the other.
   On tomato coverage runs from 16.2 reads a base at 20% GC to 29.0 at 36% — **a factor of 1.79, which
   is larger than the doubling being looked for** — so an uncorrected reading is mostly GC content.
   One byte a window fixes it, against the one the mean depth already costs.
7. **§2.2's depth ladder — record that a bin's width is a source of apparent contamination, and where
   it is paid.** The count of disagreeing reads is exact and the depth is a range, and any consumer
   that divides one by the other inherits an error of up to a sixth at thirty reads a position. The
   ladder is not the thing to change; what belongs in the document is the warning and the remedy —
   sum over the range — beside the description of the encoding.

---

## 10. The byte-level record format, as I would write it

*This is a proposal, not built. It is written against the five decisions of §6.1, §6.2 and the
architecture's §1.1a and §2.2, which is what constrains it.*

**One file per sample beside its pileup, and the file is a directory followed by sections.**

| | bytes | what it is |
|---|---:|---|
| magic and version | 8 | `ngrec\0` and two bytes of format version. A reader refusing an unknown version is the whole of forward compatibility here — the file is a cache and the answer to a version it cannot read is to rebuild. |
| the thirteen recording terms | ~200 | **A table of (field name, 16-byte digest)**, not the values. Their only use is equality with a name attached, and a digest table gives exactly that without a codec that has to track every field `StrRepeatCriteria` and `ScanParams` ever grow. |
| the kept-loci digest, per megabase | 16 × blocks | Outside the sections, because it is checked before anything large is decoded. |
| the pileup it was built from | 24 | A digest of that pileup's header and its record count. Never modification time. |
| the directory | 24 × sections | One entry per section: the key (kind, read group, stratum), an 8-byte offset and an 8-byte length. Sorted by the key's own order, so a reader binary-searches and a fit that iterates is deterministic. |
| the sections | the rest | Each one independently decodable, in the directory's order. |

**A generic section is two runs of bytes and nothing else**: the five-bit depth array exactly as
`PackedDepthCodes` already holds it, then the sparse list of non-reference observations as
delta-encoded position, allele and count. The array is already the wire form — 1.25 MB per read group
at two million positions — so this section is a copy rather than an encoding.

**An SSR section is one stratum's tracts for one read group**, with the stratum's first index in the
directory entry so that indices inside the section are stratum-local. That is the one storage change
§2.2 of the architecture says this decision forces, and the format is where it becomes visible.

**Three things I would put in the format that the documents do not yet require:**

- **Each section carries its own 4-byte checksum of its bytes.** A file that is a cache will be
  half-written by an interrupted run, and a reader that finds a torn section should rebuild rather
  than fit on rubbish. The directory alone cannot tell a truncated section from a short one.
- **The directory's last entry is a sentinel giving the file's total length**, so a truncated file is
  detected before any section is read rather than at the section that runs off the end.
- **Sections are written in the order the fit reads them** — every generic section, then strata in
  band order — so a sequential read of the file is also a valid fit, and the two-phase run over a
  thousand samples is a streaming read rather than a thousand seek patterns.

**The memory arithmetic, which is what the shape is for.** At the 5,000-tract cap a section is about
50 kB a sample, so one band of three strata across a thousand samples is 150 MB; the generic half is
1.25 MB a sample, so it is 1.25 GB across a thousand and **must not be held for every sample at once**
— which is the one place today's code is wrong for a large cohort. `fit_jointly` takes
`&[SampleRecords]` and holds every sample's every section for the whole run, and the contamination
step allocates about a gigabyte of scratch at 63 samples and two million positions. That is fine at
63 and wrong at a thousand, and the format above is what makes the fix possible rather than being the
fix.

---

## 11. What this cannot say

- **Every number outside §5 and §8 is against a made-up truth.** The drawn panels grade the
  arithmetic and the model's consistency, not whether the model describes tomato. The one place a real
  truth set speaks is §8, and what it says there is a negative.
- **The drawn panels' mismapped positions are drawn from the very model the fit assumes** — a
  position with a higher error rate in every sample. Real mismapping need not look like that, and a
  generator only reproduces the mismapping someone built into it. What §5 shows is that the exclusion
  finds two markers in five on real reads and takes the median accession from 6.8% to 0.9%; what it
  cannot show is whether the positions it condemns are the right ones.
- **The depth-ladder floor of §4 was found on drawn data and its fix is measured there.** Tomato's
  own contribution from it is §5's re-run and nothing more; the accessions run from 2.4 to 30.6 reads
  a position, so different accessions sit on different parts of the ladder and the effect is not one
  number for the panel.
- **The coverage discriminator is not usable on tomato as it stands**, because the stored summary
  cannot be GC-corrected (§9, change 6). Everything §7 says about it is from readings drawn with the
  scatter model they are scored against, which is generous to it in exactly the way the other
  report's coverage arm was.
- **A contaminated sample's fraction is still attenuated** — 0.0102 for a drawn 0.030 — and §6 says
  the coordinate correction does not fix it. What the two exclusions changed is the floor it is judged
  against, not the value.
- **One drawn panel per setting.** Differences of a percentage point are not separable from the draw;
  the differences carrying the conclusions here are ten to fifty times that.
