# The benchmark trio's extra heterozygotes are 59 positions where a quarter of the reads disagree

*Research report, 2026-08-14, following
[`joint_fit_against_truth_2026-08-13.md`](joint_fit_against_truth_2026-08-13.md), which found the
excess, and
[`contamination_floor_and_duplicated_class_2026-08-13.md`](contamination_floor_and_duplicated_class_2026-08-13.md),
which eliminated the first explanation for it. Written for a reader who has read none of the
specifications.*

*Programs: `examples/ng_joint_records_walk.rs` (the trio and the tomato cohort),
`examples/ng_joint_sample_count_sweep.rs` (the drawn cohorts). Raw output under `tmp/records/`, the
scripts that read it under `tmp/`.*

---

## 1. The question, and the answer in one paragraph

Before calling any variant, ng estimates how heterozygous each sample is. Run on the three-sample
human benchmark trio — a mother, a father and their child, each with a published truth set — it
returns **1.23 to 1.28 times as many heterozygous positions as the benchmark says are there**, while
the rate at which a sample differs from the reference on *both* copies is right to one part in a
hundred. Three explanations were on the table; this report measured each.

**None of the three accounts for it.** The excess is a definite set of positions and they have been
looked at: for HG002, **76 positions' worth of heterozygosity out of 449,489, at positions the
benchmark calls the sample homozygous reference.** There **a quarter of the reads disagree with the
reference — a mean share of 0.24 against 0.49 at a real heterozygote** — and **59 of the 83 such
positions read that way in all three people at once**, which no parent-to-child inheritance explains and which
is what a stretch of genome the reference holds once, with two stretches piling reads onto it, looks
like. The estimator's model has two stories for a disagreeing read: a heterozygote, at half the
reads, or a misreading, at half a percent to four percent of them. **A quarter is nearer a half than
it is to four percent, so the fit chooses heterozygous, and it does so with a posterior of 1.000.**

| where the 1.26 came from | share of it |
|---|---:|
| the two counts being over different things | **none** — the fit places 0.0 of its heterozygosity at the benchmark's insertions, deletions and their surroundings |
| the fit stopping before it settled | **0.1%** — heterozygosity is 0.8032 per kilobase at pass 29 and 0.8039 at pass 200 |
| the read depth read as one number when it stands for a range | **0.2%** — and structurally bounded: 999 of every 1,000 trio positions carry a code that stands for a single depth |
| positions where a quarter of the reads disagree, in everybody | **the rest** |

**The depth fix was made anyway, because the trio is the one dataset where it could not matter.**
Ninety-nine point nine per cent of the trio's positions sit at the ladder's cap, which is one exact
depth; **on the 63-accession tomato cohort 53% of positions sit in a code that stands for a range**,
and that is where the defect lives. §4 is what the fix moved there and on drawn cohorts, and §5 is
what it cost.

---

## 2. What each explanation was, and what killed it

### 2.1 The two counts are over different things — no

A benchmark VCF's heterozygote count is defined over the regions the benchmark calls confident; the
fit's is defined over the positions it kept. **Here the two are the same set by construction**: the
analysed regions are the intersection of the three samples' own high-confidence region sets, so every
kept position is inside all three, and the truth was already counted over exactly the 449,489
positions the fit averages over.

What remained was subtler and is the reason this was worth checking. **The fit's model is a
substitution model** — reads carrying an insertion or a deletion are held out of it — so the truth
was counted over single-base substitutions only. A benchmark heterozygous *indel* is therefore a
position the truth count calls nothing and the fit might well call heterozygous, and there are enough
of them to matter: 71 positions in HG002 lie under a heterozygous indel's own reference span, and 873
more lie within ten bases of one. Counting those into the truth would take HG002's benchmark
heterozygosity from 0.639 to 0.797 per kilobase and the ratio from 1.26 to 1.01 — **an arithmetic
coincidence that would have closed the whole gap.**

**It is a coincidence, and the way to tell is to ask where the fit actually puts its
heterozygosity.** The fit now keeps, for every sample at every kept position, the posterior that the
sample is heterozygous there (`JointFitConfig::genotype_posteriors`). Summed inside each class of
position the benchmark names:

| the benchmark says | positions | the fit's heterozygosity there | its homozygous-non-reference |
|---|---:|---:|---:|
| a heterozygous substitution | 287 | **284.9** | 0.0 |
| a substitution with both copies non-reference | 198 | 0.0 | **196.0** |
| a heterozygous insertion or deletion | 71 | **0.0** | 0.0 |
| …the same with both copies non-reference | 20 | 0.0 | 0.0 |
| within ten bases of one | 873 | **0.0** | 0.0 |
| the sample is homozygous reference | 448,040 | **76.5** | 0.0 |
| | 449,489 | 361.4 | 196.0 |

*HG002; HG003 and HG004 give 76.7 and 71.7 in the last row and the same zeros above it.*

**The fit puts exactly none of its heterozygosity at the benchmark's indels or beside them**, which
is what holding those reads out of the model is supposed to achieve, and which means re-counting the
truth to include them would have been adding a number to one side of a comparison it does not belong
on. The two counts are over the same things.

*The same table says something the earlier report could not: the fit finds **284.9 of the benchmark's
287** heterozygous substitutions and **196 of its 198** non-reference homozygotes. Its sensitivity is
not in question; what it adds is.*

### 2.2 The fit had not settled — 0.1% of it

The trio's fit runs out of passes at 200 without reporting convergence, where the 63-accession tomato
cohort settles in 29, and a number taken from an unsettled fit is not the estimator's answer. The fit
now keeps one line per pass (`JointFitConfig::pass_trace`):

| pass | log-likelihood | largest parameter move | HG002's heterozygosity, per kilobase | the frequency density's shape |
|---:|---:|---:|---:|---:|
| 1 | −7,621,235 | — | 0.7098 | 0.500 |
| 5 | −6,430,717 | 1.00 | 0.7872 | 0.562 |
| 25 | −6,430,230 | 0.042 | 0.8030 | 0.885 |
| 29 | −6,430,226 | 0.039 | **0.8032** | 1.039 |
| 100 | −6,430,211 | 0.0057 | 0.8037 | 3.908 |
| 200 | −6,430,211 | 0.0012 | **0.8039** | 4.779 |

**Heterozygosity is settled by pass 25 and moves by 0.0009 per kilobase over the next 175 passes** —
one part in a thousand, against an excess of one part in four. **What is still moving is the shape of
the population's allele-frequency density**, climbing from 0.885 at pass 25 to 4.779 at pass 200 and
not finished. So the run is right to say it has not converged, and wrong to suggest that
heterozygosity is what has not converged.

*Why the shape wanders is worth one sentence, because it is the same fact as the sample count: the
density describes a population and the trio supplies six chromosomes to describe it with. A shape
that is barely constrained can drift a long way for a very small gain in likelihood — the last 175
passes buy 19 units of log-likelihood out of 6.4 million.*

### 2.3 The depth read as one number — 0.2% of it, and it could not have been more

A position's read count is stored as one of twenty five-bit codes: exact below nine reads, a widening
range above it, and a top bin ending at a cap of 124. **A position deeper than the cap is thinned
down to it — depth and read counts together — before it is recorded**, so the top bin is not a range
at all: every position in it sits exactly on 124.

The walk now reports where a run's positions sit on that ladder:

| | positions at an exact depth (0–8) | positions in a code standing for a range (9–97) | positions at the cap (124) |
|---|---:|---:|---:|
| the benchmark trio, ~300 reads a position | 0.0% | **0.1%** | **99.9%** |
| the tomato cohort, 2.4 to 30.6 reads a position | 45.0% | **53.4%** | 0.0% |

**One trio position in a thousand carries a code that stands for more than one depth.** The
hypothesis was that the trio, being the deep data, is where a point read off a wide bin does most
damage; the cap means the trio is past the widening region entirely, and it is tomato that sits in
it. Fitting the same reads both ways confirms the size directly:

| | the depth read as the middle of its range | the depth summed over its range | the benchmark |
|---|---:|---:|---:|
| HG002's heterozygosity, per kilobase | 0.806 | 0.804 | 0.639 |
| HG003 | 0.760 | 0.759 | 0.596 |
| HG004 | 0.804 | 0.803 | 0.654 |
| HG002's homozygous-non-reference rate | 0.436 | 0.436 | 0.441 |

---

## 3. What the excess positions are

Eighty-three positions of the 449,489 are called heterozygous in somebody by the fit and variant in
nobody by any of the three benchmark VCFs. They are not scattered one per sample:

| the fit calls it heterozygous in | positions | heterozygosity carried, summed over the three samples |
|---:|---:|---:|
| all three people | **59** | **176.5** of 224.4 |
| two of them | 19 | 38.9 |
| one of them | 5 | 5.0 |

**Seventy-nine per cent of the excess sits at positions all three read heterozygous at once.** HG003 and
HG004 are unrelated adults; a variant they both carry heterozygously, that their child also carries
heterozygously, is possible, and fifty-nine of them inside 449 kilobases where the benchmark calls
none is not.

**What the reads look like there.** Across the 177 sample-positions involved, the share of reads
disagreeing with the reference runs **0.15 to 0.50 with a mean of 0.24**. At the 287 positions HG002's
benchmark calls heterozygous substitutions, the same share has a mean of **0.49**. So the fit is
calling *a quarter of the reads* heterozygous, with a posterior of 1.000, and the reason it can is
that the model offers nothing else: within a class, a read either comes from a copy the sample
carries or is a misreading, and misreadings run at 0.49% at an ordinary position and 4.3% at a
mismapped one. Twenty-four per cent is five times too common to be a misreading and half as common as
it should be for a heterozygote, and the likelihood picks the nearer of the two by an enormous margin
— at 124 reads, 32 of them disagreeing, the heterozygous reading is about ten million times the
mismapped-homozygote reading.

**They come in runs, which is what a stretch of genome rather than a base looks like.** Six of the 59
lie between chr6:135,206,472 and chr6:135,206,504 — six positions inside 33 bases. Three more lie
between chr2:236,202,603 and chr2:236,202,628. A base that misreads is a property of a base; a run of
them 33 bases long is a property of a region.

**This is the class of position the design has always said the two-class model cannot reach**, and it
is the one the duplicated-stretch class was built for. That class was tested on this trio and fitted a
weight of exactly zero (`contamination_floor_and_duplicated_class_2026-08-13.md` §8), and §6 below is
why that is not surprising and what would be needed instead.

---

## 4. The depth fix, and what it moved

**What was wrong.** The likelihood was handed one depth for a code that stands for several, taken
from the middle of the range, while the count of reads disagreeing with the reference is exact. The
reference-read count is the difference between the two, so it carried the bin's width, and a
heterozygote's read share landed away from a half for a reason that has nothing to do with the
sample.

**What it is now.** The likelihood sums over every depth the code could stand for. Two things decide
how much weight each depth gets, and **the answer is short because they nearly cancel**: a deeper
position has more ways to have produced the reads that were seen, and a deeper position is rarer,
because a sample's read count at a position is a Poisson draw around its own coverage. Writing both
down, the factorials cancel and what is left is that the reference reads are themselves Poisson, at
the sample's coverage times the chance a read shows the reference base, cut to the depths the code
allows. Below nine reads the range is one value and this is the plain multinomial it always was.

**Both pulls are needed and neither may be dropped.** With the multinomial coefficient alone the sum
collapses onto the deepest depth in the range; with the Poisson alone onto the shallowest. Either way
the fix becomes a second point read, at an edge of the bin instead of its middle, which is worse than
what it replaced.

**And the expectation had to move with it.** The error rate is maximised over expected read counts,
and an expectation taken at the middle of the range while the likelihood sums over the range is not
the same statement twice. Left inconsistent, the two disagree in one direction — the likelihood
prefers a deeper position for a homozygous-reference sample that showed a couple of disagreeing
reads, the middle of the range books it fewer reference reads than that — and **the fitted error rate
came back 24% above the truth on a drawn cohort at eight reads a position** where the consistent pair
returns it to within 4%. That was found by the drawn cohorts in the module's own tests and is the
single largest thing this work got wrong before getting it right.

### 4.1 The drawn control, at thirty reads a position

A drawn cohort with nothing planted in it is the control that says whether an excess belongs to the
mechanism or to the estimator. At three and eight reads a position it showed almost nothing, because
the ladder is exact or nearly so there. **At thirty reads it showed the estimator returning
heterozygosity 4 to 7% above what the cohort was drawn with**, at every sample count, which is the
one number `joint_fit_against_truth_2026-08-13.md` §5 flagged as the estimator's own.

200,000 positions, drawn heterozygosity 1.491 per kilobase, every plant drawn 0.15 less heterozygous
than random mating predicts. **Fitted ÷ drawn heterozygosity:**

| samples | 2 | 3 | 5 | 10 | 25 | 50 |
|---|---:|---:|---:|---:|---:|---:|
| nothing planted, the depth read as the middle | 1.07 | 1.06 | 1.05 | 1.05 | 1.04 | 1.04 |
| nothing planted, **summed over the range** | **0.99** | **0.92** | **0.92** | **0.95** | **0.95** | **0.98** |
| 1 position in 200 planted mismapped, the middle | 1.12 | 1.12 | 1.11 | 1.09 | 1.09 | 1.10 |
| the same, **summed over the range** | **0.99** | **0.91** | **0.92** | **0.95** | **0.95** | **0.98** |

**Read the last two rows first, because they are the ones with something planted to find.** A
consistent 9 to 12% overshoot becomes a 2 to 9% undershoot. **This is an improvement and it is not a
cure**: the estimator now misses low by about as much as it used to miss high at three of the six
sample counts, and only at fifty samples is it within 2%.

**Two things make the after-numbers weaker evidence than the before-numbers.** Every arm above now
runs out of 200 passes, where before only the two-sample arms did — so those are values at pass 200
rather than settled values, and §5 is what that costs. And one drawn cohort per setting means
differences of a percentage point are not separable from the draw; the 9 to 12 points these rows move
are.

### 4.2 The same cohorts at eight reads a position

The module's own drawn cohort — ten samples, 3,000 positions, eight reads a position, every parameter
known — is where the fix can be graded against five truths at once rather than one. It converges in
28 passes both ways.

| | the middle of the range | summed over the range | drawn |
|---|---:|---:|---:|
| a read misreads at an ordinary position | 0.00213 | **0.00193** | 0.00200 |
| …at a mismapped one | 0.0534 | **0.0587** | 0.0600 |
| the share of positions mismapped | 0.0277 | **0.0226** | 0.0200 |
| each sample's observed heterozygosity | 0.01850 | **0.01794** | 0.01700 |
| how much less heterozygous than random mating predicts | 0.186 | 0.212 | 0.200 |

**Four of the five move toward the truth and the fifth crosses it**, from 0.014 below to 0.012 above.

### 4.3 Sixty-three tomato accessions, which is where the ladder actually widens

No truth set here, so what this arm shows is the size of the move on the data that sits in the band
where a code stands for several depths — **53% of tomato's positions against 0.1% of the trio's**.
Two million positions, 63 accessions, 2.4 to 30.6 reads a position, everything else held identical.

| | the depth read as the middle | summed over the range |
|---|---:|---:|
| passes, and how long | 29, 594 s | 30, **666 s** |
| a read misreads at an ordinary position | 0.00333 | 0.00338 |
| …at a mismapped one | 0.0239 | **0.0254** |
| the share of positions mismapped | 0.0335 | **0.0294** |
| the population's expected heterozygosity, per kilobase | 4.886 | 4.895 |
| each accession's own heterozygosity, median | **0.907** | **0.867** |
| …its range | 0.317 – 4.816 | 0.318 – 4.648 |
| each accession's homozygous-non-reference rate, median | 1.502 | 1.517 |
| how much less heterozygous than random mating predicts, median | 0.778 | 0.789 |

**Nothing here moves by more than a twentieth of itself.** The median accession's heterozygosity
falls 4.4%, the share of positions booked as mismapped falls 12%, and the rate at which a read
disagrees at a mismapped position rises 6%. **Convergence is not affected at these depths** — 30
passes against 29 — and the run costs 12% more time.

That is smaller than §4.1's drawn cohorts at thirty reads a position would suggest, and the reason is
in the same table as everything else: tomato's accessions run from 2.4 to 30.6 reads a position, so
about 45% of positions are at an exact depth and most of the rest are in the narrow codes just above
eight, where a range is two or three values wide. **The drawn cohort at a flat thirty reads sits in
the widest codes for every position and is the worst case, not the typical one.**

### 4.4 An incidental finding, and it is larger than anything else here

The tomato run above has the **duplicated-stretch class turned off**, matching the run it is compared
against. Turned on — which is how it ships — on the same reads:

| | the class off | the class on |
|---|---:|---:|
| each accession's heterozygosity, median, per kilobase | **0.867** | **0.064** |
| …its range | 0.318 – 4.648 | 0.009 – 3.448 |
| how much less heterozygous than random mating predicts, median | 0.789 | **0.983** |
| the share of positions a sample carries an extra copy of | — | 0.0079 |
| the population's expected heterozygosity, per kilobase | 4.895 | 4.153 |
| how long | 666 s | 1,311 s |

**The median accession's heterozygosity falls thirteenfold and the homozygote excess is pinned at
0.98**, which says these plants are very nearly free of heterozygous positions. An earlier rough SNP
caller measured this cohort's pooled observed heterozygosity at 1.049 per kilobase; 0.064 is not that
number.

**This is not the depth fix.** The two runs above differ only in the class, and the class was decided
on drawn cohorts, where at fifty samples with no duplications planted it invented a weight of
0.00003 and moved heterozygosity by half a per cent
(`contamination_floor_and_duplicated_class_2026-08-13.md` §7). On the real panel it books 0.0079 —
260 times that — and takes 93% of the heterozygosity with it. **Nobody had run it on real reads with
a cohort large enough to identify it**: the trio is three samples and fitted zero, and the tomato run
in the earlier report had the class off. It is reported here because this run is where it turned up,
and it needs its own investigation rather than a paragraph.

---

## 5. What the fix costs: passes

**Every arm of §4.1 now runs out of 200 passes where most of them used to settle**, and that is the
only place it happens. Real data is unaffected: the tomato cohort settles in 30 passes against 29
(§4.3), and the drawn cohort at eight reads a position in 28 against 27 (§4.2). What is peculiar to
§4.1 is that every position is drawn at a flat thirty reads, so **every position sits in one of the
ladder's widest codes at once** — which is the worst case for the sum and does not describe any real
sample, whose depth varies along the genome. The cause is visible in the arithmetic: the
reference-read term is no longer a single logarithm but a sum over the range, and its curvature in
the error rate is gentler, so each maximisation step moves less.

**What this means for a number read off a run that says RAN OUT OF PASSES.** The trio's trace (§2.2)
says heterozygosity is settled long before the run stops and the frequency density is not, so the
label is about the density. That was measured on the trio and does not transfer: **any run reporting
that it ran out should be re-read with its own trace before its numbers are quoted**, which is what
the trace was added for.

---

## 6. Where I would look next, and it is not another hypothesis about the trio

The excess is 59 positions where a quarter of the reads disagree in everybody. Two things follow.

**A share the model has no state for needs a state, not a better fit.** Between *a copy the sample
carries*, at half the reads, and *a misreading*, at a few per cent, there is nothing. The
duplicated-stretch class was added for exactly this shape and fits a weight of zero here for a reason
that is structural rather than accidental: **the class is identified across a cohort by the absence of
samples homozygous for the non-reference allele**, and three samples that are all carriers supply no
such absence. It needs about twenty-five samples, and the trio is three.

**So the trio can show the defect and cannot fix it.** What would test the fix is the same positions
inside a larger human cohort — several dozen people over the same regions — where a position that
reads a quarter in everybody has somewhere else to go. Failing that, a per-position statement that
does not need a cohort: these positions have a run structure (six inside 33 bases) and a coverage
signature that a single sample carries by itself.

---

## 7. What this cannot say

- **The benchmark VCF is a truth set for its confident regions and for nothing else**, and every
  number here is inside them. The 59 positions are places the benchmark calls the samples
  homozygous reference; that they are wrong in the fit rather than missing from the benchmark rests
  on the read shares (0.24 against 0.49) and on all three samples reading alike, not on the benchmark
  being taken as infallible.
- **Three samples is below the floor for two of this route's parameters.** The frequency density is
  fitted from six chromosomes and the duplicated class needs about twenty-five samples. The fit says
  so about itself: the homozygote excess comes back exactly zero in all three, which is correct for
  these people and is also what an unidentified parameter looks like.
- **The drawn cohorts grade the arithmetic, not the model.** Their mismapped positions are drawn from
  the very two-class model the fit assumes, so they cannot contain a position that reads a quarter in
  everybody — which is what §3 says the trio's excess is. That is why §4.1's control is evidence about
  the depth ladder and not about the trio.
- **One drawn cohort per setting**; differences of a percentage point are not separable from the draw.
- **The depth prior is one number per sample.** The Poisson the range is weighted with is centred on
  the sample's mean depth over all kept positions, and coverage is not one number — a position in a
  part of a genome that reads much shallower or deeper than the sample's average has its range
  weighted by the wrong Poisson. What that costs is bounded by the width of a bin. The per-window
  coverage summary the records already carry is where a per-position centre would come from.

---

## 8. What I would change in the specifications

*Nothing under `spec/` or `arch/` was edited. These are the changes I would make, in the order I
would make them.*

**In [`parameter_prepass_joint_fit.md`](../spec/parameter_prepass_joint_fit.md):**

1. **§2.2 — record what the two-class model does with a position where a quarter of the reads
   disagree, with the trio's number beside it.** It calls it heterozygous, with a posterior of 1.000,
   in every sample that shows the pattern; on the benchmark trio that is 59 positions in 449,489 and
   **79% of a heterozygosity 26% above a real truth set**. The document says the model cannot reach
   such positions; what it does not say is which way it fails, and *calls them heterozygous* is the
   answer a reader needs, because heterozygosity is one of the numbers this pass exists to produce.
2. **§3.2 — the per-sample rates now have a per-position form, and it is what a truth set is compared
   against.** `JointFit::genotype_posterior` keeps each sample's heterozygous and
   homozygous-non-reference posterior at each kept position, off by default at eight bytes a position
   a sample. Without it a disagreement with a benchmark is one ratio; with it the disagreement is a
   list of positions, and on this trio the two readings led to opposite conclusions — recounting the
   truth to include the benchmark's indels would have closed the gap arithmetically while the fit
   places **zero** heterozygosity there.
3. **§3.1 — the depth is summed over rather than read from, and the expectation must be summed the
   same way.** The likelihood sums over the depths a stored code stands for, weighting each by the
   Poisson its own coverage implies. The error rate's expected read counts must be taken under the
   same sum: taking them at the middle of the range instead puts the fitted error rate **24% above
   the truth** on a drawn cohort at eight reads a position, against 4% when the two agree. This is a
   requirement, not a tuning note — the two halves of an expectation-maximisation step have to be the
   same statement.
4. **§8 — add what the fix is worth and what it costs, at the depth where it acts.** On a drawn
   cohort at thirty reads a position with nothing planted, fitted ÷ drawn heterozygosity goes from
   1.04–1.07 to 0.92–0.99; with 1 position in 200 planted mismapped, from 1.09–1.12 to 0.91–0.99. On
   the 63 tomato accessions, which are 2.4 to 30.6 reads a position, the median accession's
   heterozygosity moves 4.4% and nothing else moves by more than a twentieth of itself. The cost is
   passes, and only on the flat-thirty-read drawn cohort, where every arm now runs out at 200:
   tomato settles in 30 against 29 and the eight-read drawn cohort in 28 against 27.
5. **§3.2 — say that a run reporting it ran out of passes must be read with its trace.** On the trio
   the label is about the frequency density's shape, which climbs from 0.885 at pass 25 to 4.779 at
   pass 200 and is not finished; heterozygosity is settled by pass 25 and moves one part in a thousand
   thereafter. Reporting *not converged* without saying which parameter is still moving invites a
   correct number to be discarded and a wandering one to be quoted.

**In [`parameter_prepass_generic.md`](../spec/parameter_prepass_generic.md):**

6. **§4 — the warning added on 2026-08-13 is now a built rule and should say so, and it should say
   what a run's own ladder occupancy is worth knowing.** *Sum over the depths the code stands for* is
   what the joint fit does. Beside it belongs the fact that decides whether it matters for a given
   run: **the trio has 0.1% of its positions in a code standing for a range and tomato has 53%**,
   because a cap at 124 puts a 300-read sample past the widening region altogether. A reader deciding
   whether the ladder can hurt their data needs that number, and it costs one pass over the codes.
7. **§4 — `representative_depth` is gone from `DepthBinEdges` and nothing replaces it.** The ladder
   now exposes `depth_range`, which is a bin's own definition and what a histogram's row width needs,
   and `recorded_depths`, which is what a stored code means for a recorded position and differs in
   one place: **the top bin is the cap and nothing else**, because a position deeper than the cap was
   thinned down to it. There is no longer any method that hands back a single depth from inside a
   bin, which is the point — the two defects this ladder has caused, a contamination floor of 2.5% on
   a drawn panel holding none and a copy-number discriminator losing 37% of its enrichment, were both
   a consumer taking one value from inside a range.

**In [`parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md):**

8. **§2.2 — the thinning is not only a way of fitting deep positions into five bits; it changes what
   the top bin means, and the document should say so.** Every other bin is a range a position could
   be anywhere in. The top bin is a single depth, because thinning puts every position deeper than
   the cap exactly on it. A consumer that treats the top bin as its nominal range 98–124 understates
   a 300-read sample's depth by a tenth while its read counts are undiminished — the defect §2 of
   `joint_fit_against_truth_2026-08-13.md` reported, in the form it would come back in.
