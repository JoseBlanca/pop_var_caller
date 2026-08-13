# ng — the joint parameters fit: what is recorded at each kept locus

*Design spec, 2026-08-10. **No code yet — this settles the design.** One of three documents covering
the **joint parameters fit**, ng step 4's second route to every parameter it emits; read
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) first — it says what the route is,
what it produces and why it exists. This one settles **what each sample records at a kept locus, and
how it is encoded**. Which loci are kept is
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md); the mathematics that reads these
records is the parameters fit document.*

*Types and interfaces: [`../arch/parameter_prepass_joint_records.md`](../arch/parameter_prepass_joint_records.md).*

***It changes one decision in [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)***
*§5 — the depth ladder is twenty bins and five bits, not sixteen and four (§2.2) — and that document
carries a note saying so.*

---

## 1. What this is, and why it is its own module

**Two words carry this document, and both are named here rather than assumed.**

- **The genome walk** — one pass over **one sample's** alignments, visiting each locus of the genome in turn
  and gathering what that sample's reads say there. It is where a record is filled. One sample, one
  genome walk, and each sample is walked separately.
- **The parameters fit** — the estimation of the parameters the caller will run on, done **once for the whole
  cohort and before any variant is called**: how often a read shows the wrong base, how heterozygous
  each sample is, how inbred, how diverse the population is, how often a repeat tract gains or loses a
  copy, and how much of each sample's DNA came from a different plant. It reads every sample's records
  at once. **It is not the variant caller** — it is what the variant caller is later handed.

**Everything the parameters fit learns about a sample's reads comes from these records**, so what is missing from
them cannot be recovered later. The genome walk visits every locus once; a field not written then would need
a second traversal of the reads, which is the one thing this whole step is built to avoid. *(The one
exception is the window summary of §4, and it proves the rule: it is not in the records, so a run that
did not build it during the genome walk has to go back over the sample's pileup for it.)*

**It is a module of its own because it defines the types that hold this information, and two other
modules stand on them.** The genome walk fills a record as it finishes each kept locus; the parameters fit reads every
sample's records and computes from them. Neither has to know what the other does: the genome walk
needs no likelihood, and the parameters fit needs no knowledge of how a record is packed — the same division the
code already draws between `parameter_estimation::generic`, which shapes the data, and
`parameter_estimation::fitting`, which does the mathematics. **Because the types are the whole of it,
it can be built before either user exists**, and tested by filling a record, writing it out, reading
it back and comparing.

**Two kinds of record, because the two paths observe different things.** At an ordinary position the
observation is *which base a read showed*; at a repeat tract it is *how long a tract a read showed*.
Five per-base buckets cannot express a length, and a length distribution cannot express which base was
substituted, so the two share the selection rule and nothing about their contents.

### 1.1 Goals

1. **Hold everything the parameters fit needs and nothing it does not**, with §2's and §3's lists being the
   contract.
2. **Be readable in parts, so the cohort's records never have to be resident all at once.** An
   earlier version of this goal said the opposite — *stay small enough to keep every sample's whole
   record set in memory* — which holds at fifty samples and fails at a thousand. What has to fit is
   **one section at a time**: the generic half, or one repeat-tract stratum. §6.2 says what that
   obliges the layout to do.
3. **Distinguish "no data" from "data saying nothing"** — a locus never walked, a locus walked with no
   coverage, and a locus with reads that all matched are three different things and only the first is a
   bug.

### 1.2 Non-goals

- **It does not choose loci** ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md))
  and does not fit anything ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md)).
- **It does not say how the records are framed on disk** — the byte layout of the file is the
  implementer's. Where they live is no longer open: §6.1 settles it. The requirement that drove that
  decision is a property rather than a mechanism — **the parameters fit must reach every sample's records without
  walking the reads again** — and one file per sample beside the pileup is the cheapest thing that has
  it.

---

## 2. The generic record: four allele counts, a spare, and a depth

**Per kept position, per sample, per read group:** how many reads supported A, C, G and T, plus one
bucket for anything else — an indel, a spanning deletion. The reference base is a property of the
position and is held **once for the cohort**, not once per sample.

### 2.1 Why per-allele counts and not a count of non-reference reads

Three reasons are already recorded
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §2): two samples'
"non-reference" may be different alleles, so a spectrum built from a bare count would credit them with
an allele they do not share; the reference is not always the major allele; and contamination is
identified by *which* allele a sample's stray reads carry.

**The joint parameters fit adds a fourth, and it is structural.** The quantity being fitted at a locus is **the
frequency of an allele** ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.1), so
the allele has to exist as a thing that can have a frequency. "Non-reference" is not one.

**And a fifth, which is a requirement on the *parameters fit* rather than on the record, stated here because
this is where the temptation is.** All four counts are stored, and the parameters fit **sums over which
non-reference base is the segregating one** rather than picking the largest
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.1.1). At a position where the
population carries only the reference base the non-reference reads are errors spread over three
bases, so choosing the observed largest and scoring it as one allele's evidence is conditioning on a
maximum — small per site, one-directional, and landing on exactly the rare-frequency classes
everything downstream reads. **A record that collapsed the four counts to "the alternative one and
the rest" would make that choice unavoidable**, which is the reason the four are kept apart.

### 2.2 The depth ladder: twenty bins, five bits

**Depth is stored binned, on the ladder [`parameter_prepass_generic.md`](parameter_prepass_generic.md)
§4 settled: exact integers up to 8, then geometrically widening bins to a cap of 124 — twenty bins.**
Fine at the bottom and wide at the top, because depth 1 and depth 5 are different kinds of observation
while depth 100 and depth 105 say almost the same thing, and at three reads a plant 97 sites in 100
sit at depth 6 or below.

**Change to [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5, which stores
sixteen bins in four bits.** That choice was made for the encoding; the ladders were later measured
against each other, and **sixteen bins cost ten times the bias of twenty** — 0.55 rungs of the
error-rate ladder against 0.054, and 1.8% in each genotype frequency against 0.3%. One bit is not
worth that.

**Five bits also has room for a state four bits does not.** Twenty bins plus a *never walked* sentinel
is 21 codes; four bits holds 16 and could not express the bins alone. That sentinel is what §1.1's
third goal needs, and
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §7 already requires the
distinction without having anywhere to put it.

**Allele counts are not binned at all — DECIDED 2026-08-13.** The
code stores an exact count.

**The reason they are exact is what the counts look like, and it holds at every depth this caller has
to run at.** The entries in the sparse list are overwhelmingly **counts of one, two or three reads** —
a miscall at a position the sample is homozygous for — and every ladder of the *"exact at the bottom,
widening at the top"* shape above keeps those exact. **The tail a ladder compresses is where an allele
count almost never is**, so binning saves close to nothing however deep the run.

*How the list grows with depth, since the shallow figure is the one this document quoted before and it
is not the general case.* A position needs an entry once any read miscalls, so the share of positions
carrying one is about `1 − (1 − ε)^depth`, and at a per-base error rate near 2 in 1,000:

| a sample's depth | 3× | 50× | 100× |
|---|---:|---:|---:|
| entries per two million positions | ~12,000 | ~190,000 | ~360,000 |
| the sparse list against the 1.25 MB depth array | a twentieth | a third | **larger** |

**So at 100× the exceptions are the bigger half of the generic record**, not the rounding error the
tomato archive makes them look like — and binning the count still does not help, because the extra
entries a deep sample brings are the shallow sample's counts of one and two, only more of them.
**What bounds the list at high depth is the per-position depth cap** (§5), which subsamples a
position's reads before anything is recorded; §7.8's 300× arm exists to measure this list rather than
the ladder, and says so.

**And it keeps a fourteenth value off §5's list.** Every bin width has to travel with the sample, or
two samples' rows mean different things. An exact count means the same thing everywhere and needs no
such guarantee — which is the second reason not to bin a field whose binning would save so little.

*Soft, and it is the one number carried across from a different arrangement of the same data:* the
0.55-against-0.054 measurement was made on **pooled** cells, where many sites share a bin. Here each
sample-locus record is binned separately. The argument transfers; the measurement does not, and
§7.3 checks it directly.

### 2.3 The encoding, which is what makes this cheap

**A dense array in position order, plus a sparse list beside it.**

| part | what it is | size at two million positions |
|---|---|---|
| depth array | entry *i* is the *i*-th kept position's binned depth, five bits, **no coordinates and no index** | **1.25 MB per read group** |
| non-reference observations | index, allele, count — about four bytes each | 30–250 kB, driven by *errors* rather than variants |

**The positions are never stored**, because they are reproducible from
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §2's rule: every sample and the
parameters fit derive the identical list. Storing coordinates instead would cost about five bytes each, 10 MB
before any data was recorded.

**The depth array alone reconstructs every quiet position exactly**: depth `n` with no sparse entry
means `n` reads on the reference base, and the reference base belongs to the position.

### 2.4 Why the read group and not the sample

**Two consumers want two grains, and one object serves both because the read group sits *below* the
position.** The cross-sample work — diversity, the spectrum, relatedness — wants the *sample's* counts
at a position; the error rate is chemistry and wants the *read group's*. Summing a position's read
groups is addition of raw counts at one place, exact and free; the reverse is not recoverable.

**What it costs is a multiplier equal to the read groups per sample: 1 for 1,550 of the 1,707 samples**
in the tomato archive survey, 2 or 3 for nearly all the rest.

**And it is what lets the joint parameters fit avoid a compromise the other route makes.**
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §1 keeps the per-library breakdown only
for sites with at most four alternative reads, because a histogram key holding it everywhere would be
large. A record holds it everywhere at every depth, so the parameters fit scores each read against its own
library's rate with no bound and no pooling arm.

---

## 3. The STR record: a length distribution, a guard, and two composition counts

**Per kept locus, per sample, per read group:**

- **how many reads showed each whole-repeat offset from the reference tract length**, over a recorded
  range of **±4** with saturating end buckets;
- **one guard bucket** for reads differing by something that is not a whole number of motif copies;
- **bases compared, and the mismatching bases as a difference list** — `(which read, offset from the
  tract start, which base)` per mismatch. Offsets record *length*, and a substitution that does not
  change a tract's length is invisible to them, so without this channel the error rate cannot be
  recovered at all. **A pair of counters would do for the rate and for nothing else**
  ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §2.1): it cannot separate
  a substitution inside the tract, which interrupts the motif and changes what a repeat unit is, from
  ordinary error in the flank, and it cannot say that two reads carried the *same* interruption, which
  is what makes an interruption an allele rather than an error. **The list costs about what the
  counters cost** — at 300 base comparisons a locus and an error rate of 0.002 a locus carries about
  0.6 mismatching bases — because it is driven by the error rate, exactly as §2.3's sparse list is.
  Storing the reads instead would be 75 bytes per locus per sample two-bit packed, 3.7 GB across fifty
  samples at a million loci, against the 60–110 MB of §6.

### 3.1 The origin is the reference tract length — settled, and measured

**Not the locus's own modal observed length**, which is the alternative and was the earlier leaning.
[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1 measured it: a modal origin moves with the
reads, so a fit treating it as a fixed property of the locus is misspecified. The slippage level came
out **+50% to +408%** high, and the direction split at **0.48 where the truth is 0.17** — a 1.1-fold
asymmetry where the truth is 4.9-fold.

Both routes therefore share one origin, which is also what makes the comparison between them a
comparison of the same quantity.

### 3.2 Two widths

- **The recorded range is ±4.** Narrow is fine.
- **The allele lengths the parameters fit may place mass on reach ±6.** That is the load-bearing one, because it
  is what lets an end bucket be attributed to a distant allele rather than to a far slip.

**Narrow is fine only because the end buckets are scored by their marginal**: "at least four repeats
short" gets the sum over every offset it absorbs, never the probability of sitting exactly on the
edge. Measured: at a recorded range of ±1, on a stratum whose alleles reach three repeats either side,
the marginal rule still returns the slippage level to within 0.05%; plugging in the edge instead costs
33% of the level.

### 3.3 The guard bucket is a diagnostic, not a parameter

A read differing by a non-whole number of copies is modelled as an independent per-read outcome, so the
likelihood splits exactly into *how many reads were non-whole-repeat* times *how the rest fell across
the offsets*. Nothing about the slippage parameters is estimated from it.

**It has a threshold**: one non-whole-repeat read in ten of the reads that differ from the reference
tract length. Above that, the locus is not something this noise model describes and the parameters fit should say
so rather than fit it.

**Record what it caught, not only how many**, by the same sparse mechanism as the difference list: a
non-whole-repeat read keeps its offset and its size. A partial unit at a tract edge is alignment
ambiguity; an indel in the flank is a different thing; a locus over the threshold for the first reason
and one over it for the second are not the same locus. A bare count can raise the threshold and never
explain it.

### 3.4 The read cap must match the other route's

[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1 caps how many of a locus's reads are
entered and subsamples uniformly down to it, seeding the draw from the locus's position so that a
region-sharded genome walk and a single-threaded one keep the same reads. **The record uses the same cap and
the same seeding**, because if the two routes cap differently the comparison between them confounds
the cap with the route. At tomato's three reads a locus neither cap fires; at HG002's 300× both would.

### 3.5 What is held once per locus rather than per sample

- **the reference tract length**, the analogue of §2's reference base;
- **the stratum** — the (motif period, reference repeat count) pair — a column of the catalog
  [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.3 makes an input, so it
  costs nothing per sample and is not re-derived here. It is what that document selects on and what
  the parameters fit fits within.

**Encoding follows §2.3's shape**: a dense binned count of spanning reads per locus, plus a sparse list
of the non-zero offsets. At three reads a locus nearly every read sits at offset 0, so the sparse list
is short.

**Four states, not two.** A sample at a kept STR locus may have had no read reach the locus at all;
reads that reached it but none crossing the whole tract, so none reports a length; or reads that
crossed it, whether they showed the reference length or another. (The fourth, *never walked*, is §6's
and is a bug rather than data.) The generic set's zero-depth-against-quiet distinction is the last
pair, not the first.

**A read that did not cross the tract is not nothing, and this record has nowhere to put it.** It
proves the tract is **at least** as long as the stretch it covered — a censored observation, which
[`locus_generation_ssr.md`](locus_generation_ssr.md) records deliberately and whose admission gate is
overlap rather than spanning for exactly that reason: 7,085 such reads on chromosome 1 of tomato
SRR7279503 alone. **The gap is worth stating because the censoring is not random.** A tract longer
than a read is never crossed, in every sample at every depth, so it runs along repeat count — the axis
the slippage numbers are fitted within — and a stratum unobservable with this read length must not
look like one that was merely unlucky with coverage. **What to record is open**; the cheapest form is
one count of covering-but-not-crossing reads per locus beside the offsets, which distinguishes the
states without asking the parameters fit to score a lower bound (which [`locus_generation_ssr.md`](locus_generation_ssr.md)
leaves to the caller).

---

## 4. The third object: coverage by window, per sample

**A position where half the reads disagree with the reference looks like a heterozygote, and some of
them are not.** If a plant carries two copies of a stretch of genome that the reference holds only
once, both copies' reads align to the same place. Wherever the two copies differ from each other, half
the reads show one base and half show the other — which is what a heterozygote looks like, at a
position where the plant is nothing of the kind. **The model's only home for such a position is
*heterozygous*, so each one inflates the heterozygosity this whole pass exists to measure.**
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2 gives them a class of their
own; this section is the evidence that class is recognised by.

**Nothing at the position itself tells the two apart.** Both show an even split of two bases, at
whatever depth the sample runs at. What differs is the stretch around it: a doubled stretch collects
two copies' reads, so **the read depth over it is about twice what the rest of that sample's genome
gives**. That is the whole of the signal, and it is why the quantity has to be a depth averaged over a
window rather than the position's own depth — at six reads a position, one copy against two is six
reads against twelve, and the scatter in a single position's count is larger than the gap. Tomato's
three reads a site is half that again. The window has to be wide enough to have collected about
12,000 aligned bases before the doubling is legible at all (§4.1).

**And the comparison has to be against the sample's own coverage, not an absolute depth.** Two
measurements on tomato say why. Across one panel the samples run from 2.5 to 28.7 reads a site, so no
fixed number means "doubled" in all of them. And within a single sample, median window depth runs from
16.2 reads a position at 20% GC content to 29.0 at 36% — a factor of **1.79**, which is larger than
the doubling being looked for, so a window at an extreme of GC content reads high for a reason that
has nothing to do with copy number. Hence *relative*: each window's depth is divided by what this
sample's own coverage, at this window's own GC content, would predict.

*Whether this object is required at all is open, and the question is not this document's.* The cohort
carries a second signal for the same loci — a real half-frequency variant leaves about a quarter of the
samples in each homozygous class, a duplication leaves every sample at a half — and if that suffices,
the summary becomes optional and a two-phase run stops having to read every pileup a second time to
rebuild it. [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §11 question 11 carries
it, with the measurement that would settle it. **Until then this object is specified as required**,
because the alternative is unmeasured and this one is not.

**None of that can come out of §2's records, which is why this is a third object and not a field.**
They hold one binned depth per kept position, and the kept positions are one in a few hundred — a
500 bp window holds one or two of them, which is the per-position measurement just ruled out. The
summary is over **every** position the genome walk visited, not over the kept ones.

**What it holds, per sample:**

- **mean depth in each fixed window of the reference** — 500 bp windows, which over tomato's 782 Mb is
  1.6 M windows and over GRCh38's 3.1 Gb is 6.2 M. One `f32` or a binned byte each;
- **the sample's own depth-against-GC curve**, a few hundred numbers, because coverage tracks GC
  content and an uncorrected window near an extreme of it reads high or low for a reason that has
  nothing to do with copy number. Production computes both in its Stage-1 pileup
  ([`../../specs/hidden_paralog_filter.md`](../../specs/hidden_paralog_filter.md) §2); `src/pileup/`
  is frozen, so ng builds its own.

**Decision, 2026-08-13 (owner): this summary is never written to disk.** Unlike the records, it is
recoverable from the per-sample pileup — every position's depth is in there — so storing it would be a
cached copy of something a later pass can rebuild, and the owner's call is not to keep that copy. It
is therefore built twice, by whichever code is holding the depths at the time:

- **in the run that goes straight from alignments to a fit**, the genome walk builds it as it goes, because
  the genome walk already has every position's depth in hand;
- **in the two-phase run** — alignments to a per-sample pileup file, then a cohort call some time
  later — the parameters fit builds it by reading the pileups. That read costs a full pass over each sample's
  pileup, which is the price of not storing it (§6.1 prices the same pass for the records, which *are*
  stored).

**One condition makes the rebuild possible, and it fails quietly if it is missed.** A window's mean
depth is the sum of depths divided by **the positions in scope**, and scope comes from the reference
and the analysed regions — never from the positions that happened to carry a read. A window whose
reads are missing genuinely has a low mean, and that is the truth about it; dividing by the covered
positions instead makes it read as normal, and it also shifts the sample's median, which every other
number here is expressed relative to. **So a rebuild from a pileup works only where the analysed
regions are recoverable from the pileup's own header.** Production's `.psp` header carries that
provenance; ng's must keep it.

**Size: 1.6 to 6.2 MB per sample at one byte a window**, which at fifty samples is 80 to 310 MB. That
is now a **resident** cost during the parameters fit and not a stored one, and it is still the same order as the
records themselves (§6) — which is the reason this section exists rather than a sentence in the parameters fit
spec. **§7.8 measures it beside the records rather than instead of them.**

**Two properties it inherits.** It is **per sample and needs no cohort**, which matters because this
caller must also run on one sample; and its window grid is a function of the reference alone, so two
samples' summaries are comparable by construction.

*The window size still travels in §5's identity, and the decision above weakens what that check is
for.* No summary is ever written, so the parameters fit builds every sample's at one width in one process and
there is nothing left to disagree about. The value costs nothing and stays; it is the one member of
§5's list whose subject has gone, and it should be retired rather than defended if a later revision
finds no path that can still pool two summaries built at different widths.

### 4.1 The stored window is 500 bp; the width the parameters fit reads at is the sample's own

**Measured, 2026-08-12** ([`../reports/duplicated_locus_probe_2026-08-12.md`](../reports/duplicated_locus_probe_2026-08-12.md)
§4). A window's mean depth separates one copy from two only once the window has collected about
**12,000 aligned bases**, which is depth times width. Enrichment of the joint cell — a two-copy
window *and* a near-half alternative fraction — over what independence predicts, at 500 bp:

| mean depth | 2.51× | 3.60× | 5.15× | 9.89× | 13.32× | 25.20× | 28.69× |
|---|---:|---:|---:|---:|---:|---:|---:|
| enrichment | 1.6× | **1.3×** | 1.5× | 2.5× | 7.7× | 24.0× | 24.9× |

**At 3.6× a 500 bp window's mean depth is scatter**: its two-copy band swells to 9.7% of positions,
eleven times a deep sample's share, and those positions are no likelier to read near half than any
others. Widening the window recovers it — the 2.5× sample goes 1.6× at 500 bp, 8.5× at 2 kb, 15.1×
at 10 kb — and 2.51× at 5 kb (12,550 bases a window) returns what 25.2× at 500 bp (12,600) does.
**The deep sample gains nothing above the floor**, so this is not an argument for a wider grid.

**So the width is two decisions, not one.**

- **Stored at 500 bp**, as this section prices it. Storing at each sample's own width would break the
  one property §4 leans on — that two samples' summaries are comparable by construction — and put a
  per-sample number into a value §5 requires every sample to agree on.
- **Read at whatever width the sample's depth requires**, by summing adjacent windows in the parameters fit.
  Summing binned means back to a wider mean needs the per-window position counts, which the summary
  already has as its denominators. **Summing is exact and free; unsumming is not possible**, which
  is the whole reason the fine grid is the stored one.

**The tomato archive is the cohort that makes this load-bearing.** At its three reads a site the
class is invisible at 500 bp and plain at 5 kb, so a fit that read the stored grid directly would
find no duplications in exactly the samples the class was adopted for.

### 4.2 What the gating measurement returned

The class is kept. On tomato SRR7279482 at 25×, **1 position in 8,600** sits in a window near two
copies and reads between 35% and 65% alternative; the near-half rate inside those windows is
**1.26%** against **0.033%** elsewhere. **The population is 150 to 590 positions per two million**
across the eight samples walked — about a third of the same sample's near-half positions in
ordinary-coverage windows, and thirty times smaller than
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2 assumed before it was
measured. That does not withdraw the object: it is a real population, concentrated 24 times where
coverage says it should be, and nothing else in the model can hold it.

**And the summary must be per sample rather than per locus**, which the same run checked: of 84
windows some sample reads near two copies, **40 are read that way by exactly one of the eight** and
11 by seven or eight. Both a shared component — the reference's own collapse — and copy number
segregating in the panel are present, and the segregating one is the larger.

**One thing it did not settle: whether the GC correction helps.** The depth-against-GC curve is real
on tomato — median window depth runs from 16.2 reads a position at 20% GC to 29.0 at 36%, a factor
of 1.79 — but correcting for it *lowered* the enrichment on SRR7279482 from 32.6× to 24.8×, by
adding two-copy windows that carry no near-half signal. The curve stays in the summary, because it
costs a few hundred numbers and the 1.79-fold range is larger than the signal it would otherwise
swamp; **whether the parameters fit divides by it is a flag, and the measurement that would settle it is the
same genome walk on a genome whose duplications are known.**

---

## 5. What travels once per sample, beside the records

**Thirteen values, and the parameters fit must refuse to pool samples that disagree on any of them.**

Seven identify the loci and are
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §5.1's: the seed, the reference
digest, the analysed-region-set digest, **the repeat catalog's build settings and scoring weights**,
**the STR routing criteria this run asked it for**, the generic target count, and the STR per-stratum
cap.

**An eighth is the only one that checks the answer rather than the question** — a digest of the loci
actually kept, computed as these records are filled and blocked per megabase so a mismatch names where
it happened ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §5.1,
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.2). The other seven all
agree when a hash function or a threshold's arithmetic has changed underneath them; this one does not.
**It must be produced by the code that writes the entries**, not by re-running the selection, or it
witnesses nothing about this file.

**Five more are this document's**, and every one of them is load-bearing rather than informational.
The first two say what was recorded; the last three say **in what units**, and an earlier version of
this section had none of those three — which left the whole check able to pass while two samples' rows
meant different things.

- **the per-stratum locus counts** — for every stratum, how many loci the analysed regions hold and how
  many were kept. Without them anything pooled across strata is biased and silently so
  ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.3);
- **the per-locus read cap** (§3.4), for the same reason as the other caps: two samples recorded under
  different caps did not record the same thing;
- **a digest of the depth ladder's edges** (§2.2). The generic record stores a five-bit *code*, not a
  depth. Two samples binned under different edges hold codes that mean different depths, **and all
  the other values agree** — the loci were the same, the seed was the same, the digest of the kept
  loci matches because the loci did match. A code is only a number until something says what ladder
  it indexes;
- **the per-position depth cap** — the depth above which a position's reads are subsampled down before
  anything is recorded. It is the generic path's twin of the STR read cap on the row above, it moves
  independently of the ladder's own top rung, and a sample recorded at a different one did not record
  the same evidence;
- **the coverage-by-window size** (§4), where a window summary exists at all. Windows of different
  widths are not comparable and a relative copy number computed across two grids is meaningless.
  **This is the weakest of the thirteen since 2026-08-13**, because the summary is no longer written
  anywhere (§4): the parameters fit builds every sample's itself, at one width, so the disagreement it guards
  against has no way to arise. It is kept because it costs one comparison.

*None of the five needs a new mechanism: they are equality comparisons beside the eight above, and
they fail in exactly the same silent way.*

---

## 6. Cross-cutting concerns

**Memory, and the STR half is now measured rather than guessed.** An earlier version of this section
priced the STR set at *"a couple of megabytes per read group"* without a locus count behind it. There
is one now: at the STR path's calling floors, tomato SL4.00 holds **462,701 STR loci** in 141 strata,
and every one of them is kept because the per-stratum cap never has to fire
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.5,
`examples/ng_joint_loci_probe.rs`).

| object | per sample, per read group | at fifty samples |
|---|---:|---:|
| generic depth array, 2 M positions at five bits | 1.25 MB | 63 MB |
| generic sparse non-reference entries | 30–250 kB | 1.5–13 MB |
| STR set at 462,701 loci — offsets, guard, censoring count, base-comparison denominator | ~4–5 MB | 200–250 MB |
| STR difference list, driven by the error rate | ~0.3 MB | 15 MB |
| **coverage by window** (§4), per sample and not per read group — **resident only, never stored** | 1.6 MB | 80 MB |

**So the STR records are the larger half of the bill, not the smaller** — the opposite of what the
earlier arithmetic implied, and it follows directly from keeping every STR locus rather than a
sample of them. Roughly **360–420 MB for a fifty-sample cohort**, held while the parameters fit runs. **These are
the only objects that accumulate across the cohort**; everything else in step 4 is dropped when its
sample finishes.

**That total is what a run holds when it holds everything, which §6.2 says it need not.** The column
above is the right one for a fifty-sample cohort taking the simple path; at a thousand samples the
number that matters is the largest single section, because the generic half is dropped before the
repeat-tract half is read and the strata are read one at a time.

*The row that would change most is the STR one, because it is the one whose per-locus size depends on
an encoding that is not written yet.* **§7.8 measures it rather than trusting this table**, and if the
number lands badly the knob is the per-stratum cap — which exists, costs nothing to use, and until
now had no reason to fire.

**Measured, 2026-08-13 — two rows are wrong and two are exact.**
The writer was driven through real alignments for the first time: 63 tomato accessions over 8 Mb of
SL4.00, and HG002 at 30× over 6.1 Mb of GRCh38
([`../reports/joint_records_on_real_alignments_2026-08-13.md`](../reports/joint_records_on_real_alignments_2026-08-13.md),
`examples/ng_joint_records_walk.rs`), each at two million kept positions.

| row | measured | against the row above |
|---|---|---|
| depth array | 1.250 MB at 2 M positions | **exact** |
| sparse non-reference entries | 9,213 at 2.4× to 331,036 at 30.6×, which is **37 kB to 1.3 MB** at four bytes each | the range describes a shallow sample; it understates 30× fivefold. The mechanism — driven by the error rate, not by variants — holds, and the row needs a depth beside it |
| STR set | **25.0 bytes a locus a read group**, identical on both genomes, so **11.6 MB** at 462,701 loci and **578 MB** at fifty samples | 2.3 to 2.9 times the row. This is the row the paragraph above predicted would move; the knob named there is the one to reach for |
| STR difference list | 0.054 mismatching bases a locus at 2.4× and 0.585 at 30×, so 0.2 MB and **2.2 MB** at 462,701 loci | right at three reads a site, ten times low at thirty |
| coverage by window | 1 byte a window, so **1.56 MB** over tomato and **6.2 MB** over GRCh38 | **exact** |

**One number for the byte format, which §1.2 leaves open.** A sparse entry is **12 bytes in memory**
against the four this table prices — a position, an allele and a count, none of them packed. Four is
reachable, since a position needs 21 bits and an allele 3, but only by packing: it is a constraint on
that format rather than a property the code already has.

**And the difference list was measurable only after it was written.** §3 specifies it and the writer
was not filling it; it is filled now, for the reads whose tract is the reference's length. A read that
slipped a whole unit has no base-for-base correspondence with the reference — that is the aligner's
answer, not the writer's — so it contributes its offset and **nothing to the base-comparison
denominator**, which keeps the STR error rate a ratio of two quantities counted over the same reads.
At tomato about 95 reads in 100 are unslipped.

**Concurrency.** A region-sharded genome walk fills the entries for the positions in its own region, and
merging is concatenation in position order. Nothing needs communication between shards.

**Errors.** Three states must be distinguishable after a write and a read: **never walked** (a bug),
**walked with zero depth** (data), and **reads present, none non-reference** (data). §2.2's fifth bit
is what makes the first expressible. **At an STR locus there is a fourth** — reads reached the locus
but none crossed the tract (§3) — and it is the one with no field today.

**And the first real genome walk collapsed the first two, 2026-08-13.** The generic locus generator emits
nothing at a position no read reached, so silence from the genome walk was indistinguishable from a region
the run never opened: **93,150 of 1,999,404 kept positions on a 25× tomato accession came back as
*never walked*, and every one of them was data.** A bit that can express a state is not the same as a
genome walk that produces it. **The genome walk must therefore say which stretches it covered, whatever it found in
them** — `RecordWriter::mark_walked`, called for each region handed to the generators; a real depth
never overwrites the mark and the mark never overwrites a real depth, so the order does not matter.
The share is **1.4% to 35% of kept positions across the 63 tomato accessions** and 0.5% on HG002 at
30×, and it is the denominator of every rate the parameters fit reports.

**Determinism.** The read subsample of §3.4 is seeded from the locus's position, so the same sample
walked in one region and in many produces byte-identical records.

### 6.1 Where the records live, and the two ways they are built — DECIDED 2026-08-13

**Decision (owner): one file per sample, written beside that sample's pileup, never inside it.**

**Why they are stored at all, given the pileup already holds the evidence.** They are a cache: every
number here can be recomputed from the pileup, so this is a question of cost rather than of
availability. The cost is that **rebuilding is a full pass, not a seek**. The kept positions are
scattered — one every 390 bases at two million on tomato — while a pileup block holds about a megabyte
of records covering thousands of consecutive positions, so every block contains several kept positions
and a rebuild decompresses the whole file. Caching is worth it because **a sample's records do not
depend on which other samples are in the cohort**: written once, they serve every future cohort call,
which is the same argument that justifies the pileup itself. Without the cache, each cohort call pays
one full pass over every sample's pileup before it can start.

**Why beside the pileup and not inside it.** Three reasons, and the first is this project's own
history:

- **A derived section inside a large file makes every change to it a rewrite of that file.** The
  per-sample summary section already lives inside production's `.psp`, and bumping its version meant
  regenerating every existing `.psp`. What is here will move more often than that summary did: the
  budget is a per-run knob ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md)
  §4.3), the selection carries a seed, the depth ladder can be re-cut, the repeat catalog can be
  rebuilt. Beside the pileup, "rebuild it" means deleting a small file.
- **More than one census can coexist** — a small one held in memory for the pooled rates and a much
  larger one streamed for contamination — as more than one file. Inside the pileup they compete for
  one section.
- **It is the clearer contract for whoever runs this.** A file that can be deleted and regenerated says
  what it is; a section inside the pileup means regenerating it costs a rewrite of the pileup.

**Two ways it is built, and they must be one builder.**

1. **During the genome walk**, in both kinds of run. The genome walk visits every position anyway and deciding
   whether a position is kept is a hash test, so building the records there is nearly free. In the run
   that goes straight to a fit, the records stay in memory; in the two-phase run, the same genome walk that
   writes the pileup also writes the records file.
2. **From an existing pileup**, which is the regeneration path rather than the normal one. It serves
   three cases: pileups written before this file existed, a records file lost or built at knobs that
   have since changed, and — the case that justifies writing the code — **wanting a larger census than
   the one on disk**, without going back to the alignments.

**Reading a pileup back immediately after writing it is not one of the two.** It costs a full
decompression pass over a file the genome walk has just finished producing, and buys nothing.

**The hazard the second path introduces, and how it is closed.** Once both exist they must produce the
same records from the same sample, or a cohort's parameters depend on which path happened to run, and
nothing in the output would show it. Two things close it: **the genome walk-time builder is fed from the same
stream of loci that the pileup writer consumes**, downstream of every read filter and depth cap the
pileup applies, so it cannot see reads the pileup will drop; and §7.12 asserts that one sample built
both ways gives byte-identical files.

**Staleness.** The records file names the pileup it was built from: a digest of that pileup's header —
reference, analysed regions, read filters, command line — together with its record count. **Never
modification time.** On a mismatch the parameters fit rebuilds silently when the pileup is available, and fails
naming the field that differs when it is not. That is the same shape as §5's refusal, pointed at a
different object.

**Size on disk.** About 6 MB per read group at a two-million-position census, and about 20–26 MB at
the twenty-seven-million-position one that contamination would use — roughly 1.5 GB across the
sixty-three-accession tomato cohort, read sequentially rather than held
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.4).

### 6.2 The file is read in sections, because the parameters fit is — DECIDED 2026-08-13

**Decision (owner): the file is laid out so that the generic half and each repeat-tract stratum can be
read on their own, and the reader decodes only what it is asked for.** Nothing about a records file
should require a program to hold a whole sample's evidence to use any of it.

**Two facts about the estimator make this possible, and neither is a hope.**

- **The two halves are fitted one after the other, not together.** The repeat-tract half consumes
  exactly one number per sample from the generic half — that sample's homozygote excess, which weights
  a genotype drawn from a locus's length frequencies
  ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §4.1) — and hands nothing back:
  contamination is fitted on the generic loci alone (§4.3 there), and so are the noise classes, the
  frequency density and the homozygote excess itself (§3.3 there). **So the generic records can be
  dropped before a single repeat-tract record is read.**
- **Within the repeat-tract half, a stratum is fitted on its own.** The four slippage numbers are per
  (read group × stratum), the concentration that says how monomorphic a stratum's loci are is per
  stratum, and the length spectrum a locus's frequencies are drawn from is that stratum's own. What
  crosses strata is sums of counts — the substitution denominator, the per-stratum weights — and a sum
  is accumulated as each stratum goes past rather than held.

**What the layout must therefore carry:** a directory at the head of the file giving the byte extent of
each section — the generic records per read group, and the repeat-tract records per (read group ×
stratum) — so a reader can take one section without decoding the rest. The thirteen identity values and
the kept-loci digest stay outside the sections, since they are checked before anything is read.

**And "one section at a time" has to be enforced by the interface rather than left to whoever writes
the fit.** A reader that *returns* a section lets its caller keep every section it ever asked for, at
which point a run reading a file has quietly reassembled the whole file in memory — the outcome this
section exists to prevent, arrived at without anybody deciding to. So a section is lent for the length
of a call and cannot be retained; the architecture document carries the shape
([`../arch/parameter_prepass_joint_records.md`](../arch/parameter_prepass_joint_records.md) §2.2), and
§7.16 is the test. **The unit lent is one band of strata across *every* sample**, because a tract's
length frequencies are fitted from every sample with reads there — per-sample access would be the wrong
grain — and a *band* rather than a single stratum because 68 of tomato's 141 strata hold fewer than a
hundred tracts and are fitted by borrowing from their neighbouring repeat counts
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.6).

**What it buys, and the measured limit on it.** Peak resident becomes the largest single section
rather than the sum: at two million positions the generic half is 1.25 MB per read group and the whole
repeat-tract set is 4–5 MB. **But one stratum — period 1 at 8 repeats — holds 217,812 of tomato's
462,701 kept tracts, 47% of them**
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.5), so reading a stratum at a
time cuts the repeat-tract peak by about half and not by the 141 strata there are. **The per-stratum
cap is what turns this into a bound rather than an improvement**, and that document called the memory
bill *"the first reason the cap has ever had to fire"*; this is the shape of the run where it fires.

**What the cap should be is not settled, and memory is no longer the thing that decides it.** A tract
costs about ten bytes a read group, so even a cap of 20,000 puts the largest section at 200 kB a
sample — 200 MB across a thousand samples, which is affordable, and every smaller cap more so. **So
the cap is set by what the estimator needs, and that is measured only from above**: the per-tract fit
recovered a slippage level of 0.0803 against a truth of 0.0800 at **6,000 tracts** in a stratum at
three reads a site, and the comparisons against the per-stratum model ran at **1,200 to 1,500 tracts**
with twenty samples ([`../reports/joint_str_estimator_2026-08-12.md`](../reports/joint_str_estimator_2026-08-12.md)).
Nobody has swept downwards, so **where it starts to hurt is unknown**, and a cap above a few thousand
buys nothing anyone has measured. [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md)
§6 question 1 carries the sweep that would settle it.

**Sectioning is one of three levers and bounds one of three axes.** It bounds *which object* is
resident. Reading every sample's file in genome order together bounds *how many loci* are resident,
and fitting on a subsample of samples bounds *how many samples* are
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §11, question 10). A cohort of
thousands needs all three; a cohort of fifty needs none of them, and can take the run that never
writes a file at all.

---

## 7. How we know it works

1. **Round-trip.** Write and read back a record set holding every corner: a never-walked position, a
   zero-depth position, a position with reads and none non-reference, a position at the depth cap, and
   a multi-allelic position. All must come back distinguishable and unchanged.
2. **The STR ends saturate and the guard catches.** A read at an offset beyond the recorded range must
   land in the end bucket rather than being dropped or wrapping, and a read whose tract differs by a
   non-whole number of motif copies must land in the guard bucket, with its offset and size. Both are
   cheap to assert and both are silent when wrong.
3. **The difference list places a mismatch where it happened** (§3). A substitution planted in the
   flank and one planted inside the tract must come back distinguishable, and two reads carrying the
   same interior substitution must come back as two entries at one offset — not as one entry, and not
   as a count of two. **That is the assertion an interrupted-repeat model would rest on**, and a
   read-blind encoding passes every other check in this list.
4. **The four states at an STR locus round-trip** (§3): no read reaching the locus, reads reaching it
   but none crossing the tract, reads crossing it, and the region never walked. Plant the second
   deliberately — it is the state with no field, so the test either shows a lower bound surviving or
   shows that it does not.
5. **The depth ladder costs what it was measured to cost.** §2.2 adopts the twenty-bin ladder on a
   measurement made against pooled cells, and here each sample-locus record is binned separately.
   Fit synthetic data at full depth resolution and at the ladder; the gap must be within the 0.054
   rungs that measurement reports. **This is the one number carried across from a different
   arrangement of the same data**, so it is the one to check rather than assume.
6. **Read groups fold exactly.** A sample's counts at a position must equal the sum of its read
   groups' — raw counts at one place, so the equality is exact and not approximate. Assert it on a
   sample carrying two read groups.
7. **Sharded recording is exact.** The same sample walked in one region and in many must give
   byte-identical records, including the STR read subsample.
8. **Size is measured, not assumed.** §6's figures are arithmetic. Measure the objects at rest on two
   runs that stress different axes: **HG002 at 300×**, where the sparse list grows with depth × error
   rate; and **the whole tomato cohort**, where what grows is the number of samples held at once.
   **Report the generic set, the STR set and the window summary separately** — §6 says the STR set is
   the larger half at 462,701 loci, and a single total would hide it being wrong. *At 300× the
   per-position depth cap fires long before the ladder's top rung does, so what that arm measures is
   the sparse list and not the ladder.*
9. **The five unit values refuse.** Write two record sets identical in every way except one of: the
   depth ladder's edges, the per-position depth cap, the STR read cap, the per-stratum counts, the
   coverage window size. Each must be refused, naming which. **This is the test that would have failed
   before this revision and passed after it**, because the eight identity values that existed then all
   agree in every one of those five cases (§5).
10. **The window summary is comparable across samples and is not derivable from the records** (§4).
    Two samples' summaries must be built on the same grid — the grid being a function of the reference
    — and a summary built at a different window size must be refused rather than resampled. **Assert
    also that a window's mean depth differs from the mean over the *kept* positions inside it**, on a
    fixture where the two diverge: that inequality is the reason the object exists, and an
    implementation that quietly derived one from the other would pass every other check here.
11. **Summing windows gives the mean over their union** (§4.1). Ten adjacent windows summed must equal
    one genome walk's mean over the same 5 kb, exactly — it is a ratio of two sums the summary already holds.
    **Plant a short window in the run**, one the analysed regions or the ambiguity mask cut down, and
    assert that weighting it by its own position count changes the answer: a sum that treated every
    window as full would agree with the direct mean everywhere else and be wrong only where it
    matters, which is at every contig end and every region edge.
12. **The two builders agree, byte for byte** (§6.1). Build one sample's records during the genome walk and
    build them again from the pileup that same genome walk wrote, and compare the files. **This is the test
    that decides whether the pileup really holds everything the records need**, and it fails on
    exactly the fields that turn out not to be recoverable rather than on all of them, so run it on a
    fixture carrying every corner §7.1 lists — including a repeat tract, whose per-read length is the
    field most likely to be missing.
13. **A stale records file is refused or rebuilt, and never used** (§6.1). Write records, then change
    the pileup they name — a different analysed-region set is the cheapest way — and assert that the
    parameters fit rebuilds when the pileup is there and fails naming the field when it is not. Assert also that
    touching the pileup's modification time alone changes nothing, since the check must not key on it.
14. **A smaller census is a subset of a larger one, and the parameters fit says so**
    ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3). Build records at the
    large budget, take the subset a smaller target selects, and assert it equals records built at that
    smaller target directly — same positions, same values. **This is what makes "write large once" safe**,
    and if it ever fails, every run that took a subset instead of rebuilding was reading a different set
    of loci than it thought.
15. **The rebuilt coverage summary equals the genome walk's** (§4). Build the window summary during a genome walk and
    again from the pileup that genome walk wrote, and assert they match. **Plant a window inside the analysed
    regions that no read reached**: it is the one the two paths disagree about if the rebuild takes its
    denominator from the covered positions rather than from the reference and the analysed regions, and
    it is the failure that also moves the sample's median and so every other window with it.
16. **Reading one section touches only that section** (§6.2). Ask a records file for one stratum's
    tracts and assert two things: the values match what a whole-file read gives for that stratum, and
    **the bytes actually read are the section's own**, which a counting reader makes checkable. The
    second half is the one worth writing — an implementation that decodes everything and then returns
    a slice passes the first half and delivers none of the memory this section exists for.
