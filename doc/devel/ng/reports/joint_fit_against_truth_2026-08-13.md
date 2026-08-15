# The joint estimator against a real truth set, and what a benchmark trio said about it

*2026-08-13, later the same day as
[`joint_records_on_real_alignments_2026-08-13.md`](joint_records_on_real_alignments_2026-08-13.md),
which built the records this reads. Written for a reader who has read none of the specifications.*

*Programs: `examples/ng_joint_records_walk.rs` (walk and fit),
`examples/ng_joint_sample_count_sweep.rs` (the drawn curve). Raw output under `tmp/records/`.*

---

## 1. What this is

ng estimates the numbers a caller will assume — how often a sequencer misreads a base, how
heterozygous each plant is, how inbred, how diverse the population is — **before** any variant is
called. It is building two ways of doing that. The one measured here keeps raw evidence at a bounded
set of positions, the same positions in every sample, and fits everything once across the whole
cohort; that estimator now exists for ordinary positions and had, until this run, only ever been given
data it drew itself.

**Three things were run**, which is the comparison `parameter_prepass_joint_fit.md` §8's third
measurement asks for:

| | samples | positions | depth | is there a truth? |
|---|---:|---:|---:|---|
| the GIAB benchmark trio, on the 100 regions all three share | 3 | 449,489 | ~300 reads a position | **yes** — a benchmark VCF per sample |
| the tomato bench cohort | 63 | 1,999,404 | 2.4 to 30.6 | no |
| a drawn cohort, refitted at 2, 3, 5, 10, 25 and 50 of its samples | 2 → 50 | 200,000–400,000 | 3, 8, 30 | yes, by construction |

---

## 2. A defect the deep data exposed, and it moved every number

**Depth is not stored exactly.** A position's read count is kept as one of twenty bins in five bits,
and the top bin ends at 124 — because two samples' records have to be small enough to hold the whole
cohort at once. The design has always said that a position deeper than the cap is **thinned down to
it** before anything is recorded, and the cap is one of the thirteen values two samples must agree on
before they can be pooled.

**Nothing was doing the thinning.** On the trio at 300 reads a position, the consequence is not
subtle: the stored depth reads back as about 111 while the count of disagreeing reads is undiminished,
so a heterozygote showing 150 reads of one allele is recorded as having 150 alternative reads out of
111. The likelihood then charges it a negative number of reference reads, which is clamped to zero,
and the position reads as homozygous for the alternative allele.

Thinning was added — proportionally and without a random draw, so a genome split across parallel
workers still produces byte-identical records — and the top bin now answers the cap rather than its
own midpoint. What that one fix moved, on the same reads:

| | before | after | truth |
|---|---:|---:|---:|
| positions where all three carry only a non-reference base | 0.460/kb | **0.230/kb** | 0.256/kb |
| positions segregating within the trio | 2.540/kb | **1.80/kb** | 1.219/kb |
| HG002's homozygous-non-reference rate | 1.059/kb | **0.436/kb** | 0.441/kb |
| HG002's heterozygosity | 1.064/kb | **0.806/kb** | 0.639/kb |
| the shape of the allele-frequency density | pinned at its bound | Beta(4.7, 8.2) | — |
| how often a read misreads a base | 0.00583 | 0.00486 | — |

**It fires on samples above about 124 reads a position, which is this human benchmark and never
tomato**, so the tomato numbers below are unaffected. The reason it took real data to find is that
every synthetic cohort drawn for this route was drawn at 3 to 8 reads a position, where the cap can
never fire.

---

## 3. The trio against its benchmark VCFs

Three samples is the mechanism at its weakest useful strength, which is what makes this the arm that
cannot be argued with: a fit that needs fifty samples to work has to show *something* at three.

The truth is counted over exactly the positions the fit averages over — 449,489 of the 452,288 bases
the three benchmark region sets share — restricted to single-base substitutions, because the model
here is a substitution model and reads spanning an insertion or a deletion are held out of it.

| | truth het/kb | fitted het/kb | | truth hom-alt/kb | fitted hom-alt/kb | |
|---|---:|---:|---:|---:|---:|---:|
| HG002 | 0.639 | 0.806 | **1.26×** | 0.441 | 0.436 | 0.99× |
| HG003 | 0.596 | 0.761 | **1.28×** | 0.383 | 0.381 | 0.99× |
| HG004 | 0.654 | 0.806 | **1.23×** | 0.458 | 0.458 | 1.00× |

**The homozygous-non-reference rate is right to within one part in a hundred in all three samples.**
That is the harder of the two numbers to get by accident: it is the rate at which the sample differs
from the reference on both copies, and nothing in the fit is tuned to it.

**Heterozygosity comes back one and a quarter times the truth**, and the excess is one thing rather
than three: the fit judges **1.80 positions per kilobase to be segregating within the trio where the
truth is 1.219** — 1.48 times too many — and everything else follows. Two candidate causes, and this
run cannot separate them:

- **Real mismapping the model does not describe.** The fit puts 1 position in 108 in its mismapped
  class at a disagreement rate of 4.3%. Whether that absorbs what real mismapping does is exactly the
  open question §2.2 of the design was written around.
- **Three samples is very few for a population frequency.** The shape of the allele-frequency density
  is a cohort-level parameter, and with three samples the data that constrains it is six chromosomes.

**Two checks that say the fit is not simply confused.** Every one of the three comes back with a
homozygote excess of exactly zero — the fit's own statement that these people are not inbred, which is
true and which it was free to get wrong. And the shape it settles on is no longer at the edge of its
allowed range, which it was before §2's fix.

**One thing to be careful about**: the residual the histogram route reports for the same phenomenon —
1.09 times the benchmark's heterozygosity — was measured on **one** sample at 30 reads a position over
a different region set. It is not a like-for-like comparison with 1.26 at three samples and 300 reads
a position, and putting the two numbers in one sentence would suggest a ranking neither measurement
supports.

---

## 4. Sixty-three tomato accessions

No truth set, so what this arm shows is that the estimator behaves at the real shape: the depth, the
sample count and the mismapping a crop reference has. Two million positions, 63 accessions held at
once, **converged in 29 passes and 652 seconds**.

| | |
|---|---|
| a read misreads a base | 0.0033 at an ordinary position, 0.024 at a mismapped one |
| positions in the mismapped class | 1 in 30 |
| the population's expected heterozygosity | 4.89 per kilobase |
| each accession's own heterozygosity | 0.32 to 4.83 per kilobase, median **0.91** |
| each accession's homozygous-non-reference rate | 0.47 to 7.41 per kilobase, median 1.50 |
| how much less heterozygous than random mating predicts | 0.23 to 0.90, median **0.78** |

**The homozygote excess is the number to look at, and it is what a selfing crop should give.** A
median of 0.78 says a typical accession is heterozygous at about a fifth of the positions random
mating in this panel would predict — which is what repeated self-pollination does, and which nothing
in the fit was told to expect. The two numbers are consistent with each other and with the panel: the
population's 4.89 per kilobase times the 0.22 that survives the inbreeding is 1.08, against a measured
median of 0.91.

**Six accessions appear twice, sequenced on separate runs, and the pairs agree.** Their fitted
inbreeding differs by 0.011 to 0.081 on a scale of 0 to 1, and their heterozygosity by 3% to 30% —
so the estimator is reading the plant rather than the run, at least to that resolution.

*One external check: an earlier rough SNP caller measured this cohort's pooled observed heterozygosity
at 1.049 per kilobase, against the median of 0.91 here.*

---

## 5. The drawn curve in sample count — a null result, and it is the estimator's own model that makes it one

The design asks for the curve that says **where having many samples at one position starts paying**,
because that is the number a user needs and the trio alone cannot give it. A cohort is drawn with a
share of positions planted as mismapped, and the same drawn positions are refitted using the first 2,
3, 5, 10, 25 and 50 samples. **Every arm is run twice** — once as drawn and once with no mismapped
positions planted at all — so that an excess appearing in both would be exposed as a property of the
estimator rather than of the mechanism.

At tomato's shape (3 reads a position, 1.5 heterozygous positions per kilobase, 1 position in 100
mismapped at a 6% disagreement rate) and at a harder one (8 reads, mismapped positions disagreeing at
**40%**, which is what a duplicated stretch of genome looks like and is the case the design says the
two-class model cannot reach):

| samples | 2 | 3 | 5 | 10 | 25 | 50 |
|---|---:|---:|---:|---:|---:|---:|
| fitted ÷ drawn heterozygosity, 3 reads, mismapped at 6% | 1.02 | 1.02 | 0.98 | 1.00 | 1.00 | 1.00 |
| fitted ÷ drawn heterozygosity, 8 reads, mismapped at 40% | 1.04 | 1.03 | 1.03 | 1.03 | 1.03 | 1.03 |
| fitted ÷ drawn heterozygosity, 30 reads, mismapped at 40% | 1.12 | 1.12 | 1.11 | 1.09 | 1.09 | 1.10 |
| …the same, with **nothing** planted | 1.07 | 1.06 | 1.05 | 1.05 | 1.04 | 1.04 |

**The curve is flat, and that is not a result about sample count.** It is the framing the design
itself warns about: the mismapped positions here are drawn from the very two-class model the fit
assumes, so the fit recovers them at two samples as easily as at fifty — the share comes back at 0.0053
against a drawn 0.005 in every arm. **A generator only reproduces the mismapping someone built into
it.** So this measurement, as specified, cannot locate where the cohort starts paying; only real reads
carry mismapping the model does not describe, which makes §3's trio and §4's cohort the informative
arms and this one a control.

**What the control does say, and it is the one thing here that is a curve in sample count.** Handed a
cohort with **no mismapped positions in it at all**, the fit still books some — and how many depends
on how many samples it has:

| samples | 2 | 3 | 5 | 10 | 25 | 50 |
|---|---:|---:|---:|---:|---:|---:|
| share booked as mismapped, 3 reads a position (truth: none) | 0.0100 | 0.0081 | 0.0100 | 0.0100 | **0.0000** | **0.0004** |
| share booked as mismapped, 8 reads a position (truth: none) | 0.0047 | 0.0028 | 0.0105 | **0.0009** | **0.0009** | **0.0008** |
| share booked as mismapped, 30 reads a position (truth: none) | **0.0007** | 0.0008 | 0.0008 | 0.0009 | 0.0011 | 0.0013 |

**A shallow cohort of a few samples invents about 1 position in 100; one with enough samples finds
none.** That is the mechanism the whole route is built on, showing up in the class weight rather than
in heterozygosity — with several samples at one position the fit can see that no position behaves the
way a mismapped one does.

**Where "enough" falls is set by depth, and it runs the right way**: about twenty-five samples at 3
reads a position, about ten at 8, and two are already enough at 30. Deep reads settle the question at
one position by themselves; shallow reads need the cohort, which is the case this route was adopted
for.

It costs nothing in heterozygosity here — every arm at 3 and 8 reads is within 4% — which is why the
flat rows above are flat. On real reads it need not be free, and §3's trio is where it is not.

*One number in the table above is the estimator's own and not the cohort's:* at 30 reads a position it
returns heterozygosity 4 to 7% high **with nothing planted at all**, at every sample count. That is a
bias of the fit on its own model at high depth, it does not shrink with samples, and it is a
candidate — though not a sufficient one — for part of §3's 1.26.

---

## 6. Contamination, which the tomato panel's divergence made urgent

The first report measured this cohort's divergence at `F_st` 0.44 across its leading split. At that
divergence a single allele frequency for the whole panel does not merely lose precision: a sample
genuinely 3% contaminated comes back at half a percent and passes as clean. So contamination needs
**each sample's own allele frequency at each position**, and it is now built —
`parameter_estimation::joint::contamination`, fitting each position's frequency as a straight line in
the panel's own axes of variation, using every sample, with the slopes shrunk so that a position whose
structure is indistinguishable from noise keeps only the panel-wide frequency.

### What the implementation settled that the design had left implicit

**The two genotypes are drawn against two different frequencies.** The sample's own genotype is drawn
at *its* frequency, because that is a statement about its ancestry. The contaminant's is drawn at the
frequency of **whoever was sequenced beside it** — by default the whole panel — because a neighbouring
library on a plate is not chosen for ancestry. Scoring both against the sample's own frequency is the
obvious reading, and it is wrong in the expensive direction. Forty samples, four subpopulations at
`F_st` 0.20, three reads a position, one sample contaminated at 3%:

| the contaminant's genotype drawn at | that sample | worst of the 39 clean | a clean panel's mean |
|---|---:|---:|---:|
| **the panel's frequency** — correct | 0.0166 | 0.0032 | **0.0004** |
| the sample's own frequency | 0.0481 | 0.0195 | 0.0099 |

A contaminant from a different subpopulation carries alleles the *sample's own* frequency calls rare,
and rare alleles turning up is the contamination signature — so the wrong reading manufactures about
**1% of contamination in every clean sample**, which is the flagging threshold itself.

### It finds the sample, and it understates the fraction

| | 12,000 varying positions | 60,000 |
|---|---:|---:|
| the sample contaminated at 3% | 0.0166 | 0.0163 |
| the worst of the thirty-nine clean ones | 0.0032 | 0.0004 |

**The separation improves with the budget and the value does not.** Forty times the noise floor at
60,000 positions, which is what a threshold needs — and 0.0163 for a truth of 0.030, which does not
move, so it is a bias rather than noise.

**The cause is the one part of `verifyBamID2` not yet built.** It maximises over `α` *and the intended
sample's own coordinates together*; here those coordinates are estimated from the sample's own reads,
which the contamination has already pulled towards the panel average, so its fitted frequency sits
closer to the contaminant's than it should and the difference the estimator lives on shrinks. **Until
that lands, `α` says *this sample stands out from the panel*, not *this sample is 1.7% contaminated*.**

### The refusal

A sample sitting alone at the end of an axis has a fitted frequency that is mostly its own echo, and a
noisy frequency manufactures contamination. How much of its own frequency a sample supplies depends
only on the coordinates, so it is one number per sample for the whole run, computable before a single
position is fitted; above a half the sample is told *not identified* rather than given a number. On an
evenly filled panel of twenty the numbers run 0.05 at the middle of an axis to 0.19 at its ends, well
clear of the refusal — `(components + 1) / samples` is the panel's *mean*, not everyone's share.

### On the real tomato panel it returns a number that cannot be true, and the cause is nameable

Run on the 63 accessions over the 52,525 positions the cohort varies at, **the median accession comes
back at 6.5% contaminated**, and the highest at 12.5%. Sixty-three accessions from a public archive are
not all one part in fifteen someone else's plant. **That is a floor, not a measurement**, and reporting
it as contamination would be reporting the estimator's own noise as biology.

**The refusal did fire, once and on the right accession.** `SRS3394702` is told *not identified*, and it
is the accession sitting furthest out on the panel's leading axis — the same one at −0.426 in the first
report's structure measurement. So the leverage machinery works; it is not what is producing the floor,
because the floor is on everybody.

**What is producing it, on the evidence here.** Contamination is identified by *a small share of reads
carrying an allele the sample should not have* — and a mismapped position produces exactly that, in
every sample at once. The same fit puts **1 position in 30 in its mismapped class, at a disagreement
rate of 2.4%**, and nothing excludes those positions from the ones contamination is measured over. A
panel-wide floor of a few percent is what that would look like. The drawn panels of §6 above had **no
mismapped positions in them at all**, and their floor was 0.0004 — which is the same statement from the
other side.

**So the fix is named and it is this route's own mechanism.** The joint fit already computes, for every
position, the posterior that it belongs to the mismapped class — that is the thing having many samples
at one position buys, and §5's control shows it is recovered correctly once there are enough samples.
It is simply not surfaced: `expectation` accumulates it and throws it away. **Weighting the
contamination markers by it, or dropping the positions it condemns, is the next thing to build**, and
until it exists no contamination number from real reads should be read as a fraction.

---

## 7. What is not here

- **The other route is not run beside this one** — and per the owner, 2026-08-13, that comparison is
  not wanted: the per-sample histogram route is likely to be dropped because it cannot produce the
  population's diversity or a contamination fraction at all.
- **The duplicated-stretch class and the repeat-tract half** are not in the estimator.
- **The joint maximisation over the contaminated sample's own coordinates** — §6's bias.
- **The fit runs out of passes on the trio** — 200 without settling, where the tomato cohort settles in
  29. The best-scoring iterate is what is returned and it is marked as not converged, but a
  three-sample fit wandering is itself a symptom of the frequency density being weakly constrained
  there.
