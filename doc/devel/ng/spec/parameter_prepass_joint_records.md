# ng — the joint fit: what is recorded at each kept locus

*Design spec, 2026-08-10. **No code yet — this settles the design.** One of three documents covering
the **joint fit**, ng step 4's second route to every parameter it emits; read
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) first — it says what the route is,
what it produces and why it exists. This one settles **what each sample records at a kept locus, and
how it is encoded**. Which loci are kept is
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md); the mathematics that reads these
records is the fit document.*

*Types and interfaces: [`../arch/parameter_prepass_joint_records.md`](../arch/parameter_prepass_joint_records.md).*

***It changes one decision in [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)***
*§5 — the depth ladder is twenty bins and five bits, not sixteen and four (§2.2) — and that document
carries a note saying so.*

---

## 1. What this is, and why it is its own module

**The fit reads nothing but these records, so what is missing from them cannot be recovered later.**
The walk visits every locus once; a field not written then would need a second traversal of the reads,
which is the one thing this whole step is built to avoid.

**It is a module of its own because it defines the types that hold this information, and two other
modules stand on them.** The walk fills a record as it finishes each kept locus; the fit reads every
sample's records and computes from nothing else. Neither has to know what the other does: the walk
needs no likelihood, and the fit needs no knowledge of how a record is packed — the same division the
code already draws between `parameter_estimation::generic`, which shapes the data, and
`parameter_estimation::fitting`, which does the mathematics. **Because the types are the whole of it,
it can be built before either user exists**, and tested by filling a record, writing it out, reading
it back and comparing.

**Two kinds of record, because the two paths observe different things.** At an ordinary position the
observation is *which base a read showed*; at a repeat tract it is *how long a tract a read showed*.
Five per-base buckets cannot express a length, and a length distribution cannot express which base was
substituted, so the two share the selection rule and nothing about their contents.

### 1.1 Goals

1. **Hold everything the fit needs and nothing it does not**, with §2's and §3's lists being the
   contract.
2. **Stay small enough to keep for every sample at once**, since the whole cohort's records are
   resident during the fit.
3. **Distinguish "no data" from "data saying nothing"** — a locus never walked, a locus walked with no
   coverage, and a locus with reads that all matched are three different things and only the first is a
   bug.

### 1.2 Non-goals

- **It does not choose loci** ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md))
  and does not fit anything ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md)).
- **It does not decide where the records live between the walk and the fit** — in memory, one file per
  sample, or folded into whatever the pipeline already writes. What it requires is a property, not a
  mechanism: **the fit must reach every sample's records without walking the reads again.**

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

**The joint fit adds a fourth, and it is structural.** The quantity being fitted at a locus is **the
frequency of an allele** ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.1), so
the allele has to exist as a thing that can have a frequency. "Non-reference" is not one.

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

**Allele counts are not binned the same way.** The difference between 0, 1 and 2 reads supporting an
allele is most of the signal at low coverage, so small counts stay exact and only the tail is binned.

*Soft, and it is the one number carried across from a different arrangement of the same data:* the
0.55-against-0.054 measurement was made on **pooled** cells, where many sites share a bin. Here each
sample-locus record is binned separately. The argument transfers; the measurement does not, and
§6.3 checks it directly.

### 2.3 The encoding, which is what makes this cheap

**A dense array in position order, plus a sparse list beside it.**

| part | what it is | size at two million positions |
|---|---|---|
| depth array | entry *i* is the *i*-th kept position's binned depth, five bits, **no coordinates and no index** | **1.25 MB per read group** |
| non-reference observations | index, allele, count — about four bytes each | 30–250 kB, driven by *errors* rather than variants |

**The positions are never stored**, because they are reproducible from
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §2's rule: every sample and the
fit derive the identical list. Storing coordinates instead would cost about five bytes each, 10 MB
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

**And it is what lets the joint fit avoid a compromise the other route makes.**
[`parameter_prepass_generic.md`](parameter_prepass_generic.md) §1 keeps the per-library breakdown only
for sites with at most four alternative reads, because a histogram key holding it everywhere would be
large. A record holds it everywhere at every depth, so the fit scores each read against its own
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
  samples at a million loci, against the 60–110 MB of §5.

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
- **The allele lengths the fit may place mass on reach ±6.** That is the load-bearing one, because it
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
tract length. Above that, the locus is not something this noise model describes and the fit should say
so rather than fit it.

**Record what it caught, not only how many**, by the same sparse mechanism as the difference list: a
non-whole-repeat read keeps its offset and its size. A partial unit at a tract edge is alignment
ambiguity; an indel in the flank is a different thing; a locus over the threshold for the first reason
and one over it for the second are not the same locus. A bare count can raise the threshold and never
explain it.

### 3.4 The read cap must match the other route's

[`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1 caps how many of a locus's reads are
entered and subsamples uniformly down to it, seeding the draw from the locus's position so that a
region-sharded walk and a single-threaded one keep the same reads. **The record uses the same cap and
the same seeding**, because if the two routes cap differently the comparison between them confounds
the cap with the route. At tomato's three reads a locus neither cap fires; at HG002's 300× both would.

### 3.5 What is held once per locus rather than per sample

- **the reference tract length**, the analogue of §2's reference base;
- **the stratum** — the (motif period, reference repeat count) pair — a column of the catalog
  [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.3 makes an input, so it
  costs nothing per sample and is not re-derived here. It is what that document selects on and what
  the fit fits within.

**Encoding follows §2.3's shape**: a dense binned count of spanning reads per locus, plus a sparse list
of the non-zero offsets. At three reads a locus nearly every read sits at offset 0, so the sparse list
is short.

**Four states, not two.** A sample at a kept STR locus may have had no read reach the locus at all;
reads that reached it but none crossing the whole tract, so none reports a length; or reads that
crossed it, whether they showed the reference length or another. (The fourth, *never walked*, is §5's
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
states without asking the fit to score a lower bound (which [`locus_generation_ssr.md`](locus_generation_ssr.md)
leaves to the caller).

---

## 4. What travels once per sample, beside the records

**Ten values, and the fit must refuse to pool samples that disagree on any of them.**

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

**Two more are this document's**, and both are load-bearing rather than informational:

- **the per-stratum locus counts** — for every stratum, how many loci the analysed regions hold and how
  many were kept. Without them anything pooled across strata is biased and silently so
  ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.3);
- **the per-locus read cap** (§3.4), for the same reason as the other caps: two samples recorded under
  different caps did not record the same thing.

---

## 5. Cross-cutting concerns

**Memory.** About **1.25 MB per read group** for the generic set at two million positions plus
30–250 kB of sparse entries, and a couple of megabytes per read group for the STR set. Roughly
**60–110 MB for a fifty-sample cohort**, held for every sample during the fit. **These are the only
objects that accumulate across the cohort**; everything else in step 4 is dropped when its sample
finishes.

**Concurrency.** A region-sharded walk fills the entries for the positions in its own region, and
merging is concatenation in position order. Nothing needs communication between shards.

**Errors.** Three states must be distinguishable after a write and a read: **never walked** (a bug),
**walked with zero depth** (data), and **reads present, none non-reference** (data). §2.2's fifth bit
is what makes the first expressible. **At an STR locus there is a fourth** — reads reached the locus
but none crossed the tract (§3) — and it is the one with no field today.

**Determinism.** The read subsample of §3.4 is seeded from the locus's position, so the same sample
walked in one region and in many produces byte-identical records.

---

## 6. How we know it works

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
8. **Size is measured, not assumed.** §5's figures are arithmetic. Measure the records at rest on two
   runs that stress different axes: **HG002 at 300×**, where depth runs past the ladder's cap and the
   sparse list grows with depth × error rate; and **the whole tomato cohort**, where what grows is the
   number of samples held at once. Report the two sets separately.
