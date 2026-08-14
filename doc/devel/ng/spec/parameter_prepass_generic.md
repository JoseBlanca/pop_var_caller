# ng — the parameter pre-pass: the SNP/indel path

*Design spec, 2026-08-03. **No code yet — this settles the design.** One of five documents covering
ng step 4; the shared framing, the estimator and the map of accumulators are in
[`parameter_prepass.md`](parameter_prepass.md), which this assumes. **Scope: the generic path end to
end** — its two histograms, the per-base error rate, the sample's own rates, and the inbreeding
coefficient `F`. Everything fitted from those two histograms is here, because the document that
specifies an accumulator specifies what is fitted from it. `src/ssr/` and `src/pileup/` are frozen
production: everything said about them here is a record, not a change.*

---

## 1. What this path produces, and what it reads

**First, what this document is and is not.** The per-sample walk fills **five** accumulators
([`parameter_prepass.md`](parameter_prepass.md) §5.1): two generic histograms, the STR table, and the
two censuses. **This document covers the two histograms**, and "the generic path" is that narrower
thing throughout — not everything the SNP/indel side does. The censuses are filled by the same walk,
at the same time, over the same loci; they are specified in
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) because their consumers are
cross-sample.

**Four numbers, out of two accumulated objects.**

| parameter | grain | fitted from | section |
|---|---|---|---|
| per-base error rate `ε` | read group | the read-group histogram | §3 |
| observed heterozygosity `Hobs` | sample | the windowed histogram, summed over windows | §5 |
| homozygous-non-reference rate `π_hom_alt` | sample *and reference* | the same, summed over windows | §5 |
| inbreeding coefficient `F` | sample | the same, window by window | §6 |

**Two accumulated objects, not three.** The walk builds the read-group histogram and the windowed
histogram. A whole-sample histogram is wanted too — §5 fits from it — but it is **not accumulated**:
it is the windowed one summed over its windows. That fold is exact and free, because summing a site's
windows is addition and no site is in two windows. Building it separately would be a third object to
keep in step with the second.

**Why two objects and not one.** Both are `(depth, alt-count)` histograms differing only in the key:

- **the read-group histogram** enters a site **once per read group that covered it**. A site with 20
  reads from one library and 10 from another appears as two entries, at depth 20 and depth 10;
- **the windowed histogram** enters that site **once, at depth 30**, in the window it falls in.

`ε` is a rate per **read** and survives that splitting; heterozygosity is a rate per **site** and does
not. So `ε` comes from the first, and everything in §5 and §6 from the second.
[`parameter_prepass.md`](parameter_prepass.md) §5.1 states the rule in general — it governs all five
of step 4's accumulators, not only these two — and records the alternatives it rules out.

**The windowed histogram keeps one more thing, and it is measured rather than argued.** When two
libraries with different error rates cover a site, the site's likelihood is one sum over genotypes
with *both* libraries' reads inside it — the genotype belongs to the site, and both libraries are
reading the same one. A key of total depth and total alternative count has discarded which library
produced which read, so it cannot weigh each read against its own library's rate.

**The key therefore holds one extra thing, and only where it changes an answer: which libraries
the alternative reads came from, for sites with at most four of them.** Above that the site pools.
It is still one entry per site, and at one library the key reduces to today's exactly.

**What that buys is not precision — it is the ability to tell two libraries apart at all.** With
the library forgotten, each read shows a non-reference base at the share-weighted rate
`Σ_g w_g·p_j(ε_g)`, where `w_g` is library `g`'s share of the reads and `p_j` is
[`parameter_prepass.md`](parameter_prepass.md) §3's per-read probability at `j` alternative copies.
Because `p_j` is a straight line in `ε`, that weighted rate equals `p_j(ε̄)` at the single
share-weighted mean `ε̄ = Σ_g w_g·ε_g`. **So the pooled key sees `ε̄` and nothing else about the
individual rates** — no amount of genome separates them, because the likelihood is exactly flat
along every combination that holds `ε̄` fixed. Keeping which library each alternative read came
from is what breaks that flatness.

**And the key has to be *scored* as a likelihood, which is a separate decision from what it
keeps.** A key that has forgotten the per-library depths must still supply them somehow. The
tempting answer is to invent them — give each library its average share of the depth, `n̂_g = w_g·n`
— and it is wrong in a way that does not shrink with data. The right answer is to sum over what was
forgotten, which for this key has a closed form and costs the same:

```text
                              n!                                                       n−k
L(n, k₁…k_G | θ)  =  Σ  π_j ────────────  ·  Π (w_g·p_j(ε_g))^{k_g}  ·  ( Σ w_g·(1 − p_j(ε_g)) )
                     j       Π k_g! (n−k)!    g                            g
```

Each read independently picks a library and then shows the alternative allele or the reference, so
the cell — how many alternative reads came from each library, and how many reads showed the
reference in total — is one multinomial over `G + 1` categories. Nothing is approximated.

**Both halves are measured, and the measurement is exact rather than simulated** — each cell is
weighted by its probability under a known truth, so what the fit returns is what an infinite genome
would return and a departure from truth is bias with no sampling noise in it. The numbers, the
sweep behind them and the two failures they overturned are
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§2. What they settle:

- **Scored as a likelihood the key is exactly unbiased** — zero in all 31 worlds tried, spanning
  error-rate ratios of 1, 4 and 10, mean depths 3 to 60, even and 90/10 splits, two libraries and
  four. It costs nothing: the same number of cells as a plug-in, and a closed form.
- **Scored by average share it is not**, and the damage has nothing to do with chemistry — on two
  libraries with the *same* error rate it reports heterozygosity 68% high and the
  homozygous-non-reference rate 78% low at three reads, from 8 sites in a thousand.
- **The pooled key was never 30 rungs out.** Its `ε̄` is exact everywhere; the ±30-rung split an
  earlier version of this section attributed to it was the plug-in inventing an apportionment the
  data does not contain, and reporting the same wrong split at 3, 6, 10, 20 and 60 reads because
  the quantity is not identified rather than badly estimated.
- **The bound of four alternative reads is a precision choice, not a correctness one** — a bound of
  two is equally unbiased on 28% fewer cells, and at three reads neither loses measurable precision
  against scoring every read against its own library.

*Retired: keeping the whole per-library breakdown, depths included, at shallow sites.* An earlier
version carried a second arm for total depths of four or less. It existed to rescue the
average-share plug-in, and rescues it only because 8 sites in 10 at three reads have four reads or
fewer and skip the plug-in entirely — the threshold was the value that hides the defect at the one
depth this cohort happens to sit at. Scoring the key properly removes the need for it, and with it a
whole arm of the accumulator.

*Rejected: one windowed histogram per read group.* It looks like the obvious way to keep the
chemistry axis, and it is the failure this section's own rule forbids — it splits each site into one
entry per library, so the two entries each draw their own genotype independently. Our example site
would score the first library's entry as evidence for homozygous-reference and the second's as
evidence for homozygous-*alternative*, one alt read out of one, and multiply them. Jointly the same
data is decent evidence for a heterozygote.

**Most samples have one read group, and then the two objects coincide.** A sample sequenced from a
single library splits no site, so its read-group histogram is exactly its windowed histogram summed
over windows — the same cells, with the same counts. That is the common case rather than a corner:
1,550 of the 1,707 samples in the tomato archive survey carry one library, and so does every sample
in both cohorts this spec's numbers come from ([`parameter_prepass.md`](parameter_prepass.md) §5).

**Build both anyway, and assert the coincidence rather than exploiting it.** Skipping the read-group
histogram whenever a sample has one library would save 4.7 kB against the windowed histogram's 37 MB
(§9) — one part in eight thousand. What it would buy is a second accumulation path, taken by nine
samples in ten, in which the object the multi-library path builds separately is never built at all.
The coincidence is worth having as a test (§12.6): an exact equality, needing no simulated truth, on
the cohorts we hold today. **It is a plumbing test and not evidence about any of the four numbers** —
both histograms reduce the same locus through the same function and bin it by the same shared
object, so the equality is close to guaranteed by the types and fails only if a field is transposed.
The checks that carry real evidence are §12's external anchors.

---

## 2. The noise model — what `ε` actually contains

> **⚠ AMENDED 2026-08-10 by the second-class-of-site milestone**
> ([`../impl_plan/noise_model_extension.md`](../impl_plan/noise_model_extension.md)). This section
> said the noise model is *a per-base substitution rate and nothing else*. **It now has two classes
> of site**, and §2.1 below states the addition. Everything else in this section — what `ε`
> contains, why one sample cannot separate its causes, why base qualities are not modelled — is
> unchanged and is what the *clean* class means.

**The generic path's noise model is a per-base substitution rate, applied at a site drawn from one
of two classes** (§2.1). Within a class, the likelihood of
[`parameter_prepass.md`](parameter_prepass.md) §3 applies unchanged: a read over a reference copy
shows some other base with probability `ε`, and a read over an alternative copy reverts to the
reference with probability `ε/3`; there is no other way for a read to be wrong. (§3 explains the `3`: three
bases to go wrong into, one to come back to.) That is appropriate where the
alternatives to a reference base are three other bases. A repeat tract can also slip a whole copy,
which this model has nowhere to put — hence the STR path's separate one
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md)).

**What `ε` measures is broader than sequencing error, and one sample cannot narrow it.** `ε` is the
rate at which a read shows an allele the individual does not carry, whatever put it there — and the
sequencer misreading a base is only the most obvious contributor. Others: a polymerase
misincorporating during library amplification, two fragments joining into a chimera, a base damaged
before sequencing began, DNA from another individual (contamination) or from another library on the
same run (index hopping), a read from a paralogous locus mismapped here, and a read at the right
locus aligned through the wrong gap. **That list is not meant to be complete, and that is the
point** — the parameter is defined by its effect rather than by its causes. Nothing in a single
sample separates the contributors, so the fitted `ε` is their sum, and the read-group grain fits the
sum because most of them are properties of a library and a run, as chemistry is.

**Mostly this is the quantity the caller wants anyway**, since it asks *how often does a read show
an alternative allele when the individual is homozygous reference*, and every contributor answers
that question. The cost is that they are not spread the same way across the genome: misreads and
polymerase errors are roughly uniform, contamination concentrates at the sites where the contaminant
differs from the reference — which are the sites polymorphic in the population — and misalignment
and chimeras concentrate at repeats and indels. One flat rate therefore understates the
alternative-read rate wherever the non-uniform contributors pile up and overstates it everywhere
else.

### 2.1 A second class of site — ADDED 2026-08-10

**One rate per read group fits the body of the distribution and not its tail.** Measured on HG002's
GIAB confident regions, at the 550,976 loci the benchmark records no variant of any kind, **818
carry three or more alternative reads where one rate predicts 29**
([`../research/noise_model_overdispersion_2026-08-10.md`](../research/noise_model_overdispersion_2026-08-10.md)).
The three-genotype mixture has exactly one class that can explain such a locus, so the surplus
arrived as heterozygosity: **1.41 times the benchmark's count** on that sample.

**So a site is *clean* with probability `1 − w` and *noisy* with probability `w`, and the genotype
emission uses that site's own rate.** The noisy rate `ε_noisy` and the share `w` are **one pair per
sample**, where `ε` is one rate per read group; at a noisy site every library's reads disagree at
`ε_noisy`. Fitted on real alignments the pair comes out at **0.4% to 1.4% of sites at 4 to 10
percent**, against clean rates of 2 to 4 in a thousand, and it cuts HG002's fitted heterozygosity
from 1.41 to **1.085** times the benchmark. A beta-binomial was tried and loses by 425 nats for one
fewer parameter.

**What makes a site noisy** (owner, 2026-08-10) is not one thing, and the three known causes do not
belong to the same object. A **duplication the reference does not carry** collects two copies' reads
at one locus, so every position where the copies differ shows alternative reads at every depth —
a property of the *genome*, shared by every library made from it. **Contamination in a library**
raises the alternative count at exactly the loci where the contaminant differs from this sample — a
property of the *library*, and two libraries of one sample can differ in it. **Error-prone sequence
context and mismapping** are partly the library's too, since mapping difficulty depends on read
length and insert size. A noisy population is therefore expected in most samples rather than in
unlucky ones.

**The pair is fitted per sample all the same, and that is an assumption rather than a measurement.**
No data distinguishes a per-sample from a per-library share: every sample in both cohorts carries one
library. **Contamination is the case that would break it first**, because it can lift one library's
share while leaving its sibling's untouched. To revisit when a multi-library alignment exists.

**A sample that wants a noisier class than the model covers is refused and fitted with one
rate — it is an outlier, not a reason to widen the model** (owner, 2026-08-10). The error-rate
ladder runs Phred 10 to 50 because that is the range of *sequencing* noise (§3); an argmax on
its coarsest rung is the search asking to leave that range, and what such a sample holds is a
population of positions this model does not describe. Two of five real alignments do it —
tomato SRR7279482 and SRR7279483, at 0.42% and 0.49% of sites — and what they are asking for
fits a duplication the reference does not carry, where about **half** the reads disagree, five
times that rung.

*Rejected: widen the ladder for the noisy class.* As the noisy rate approaches a half, a noisy
site and a heterozygous site are the same distribution, so the class that exists to take mass
**away** from heterozygosity would begin taking real heterozygotes with it. **A model flexible
enough to absorb those two samples would serve every sample that does meet its assumptions
worse, and that is a regression.** *Rejected: take the next rung down instead* — an answer just
inside the range carrying none of the evidence that the sample is outside it.

**The refusal is reported rather than silent.** `site_noise` is then `None` for a reason quite
unlike the ordinary one, so `site_noise_off_the_ladder` says which happened. Such a sample gets
the one-rate answer it would have had before this milestone, and the caller can see that it did.

**What a sample emits stays one number: the share-weighted marginal** `(1 − w)·ε_clean + w·ε_noisy`
— the probability a read disagrees with the reference at a site drawn at random. That is the
quantity a model-free count measures, and measured against one it sits **3.1% high, half a rung of
the error-rate ladder**. Emitting `ε_clean` instead would report **16% below** that count, which
[`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §9 calls an
unambiguous bug. The pair travels beside it as a diagnostic, so a consumer that wants to score a
read against its own site's class can, without another fit. **It is folded into one number only at
the emitted surface** — the runs model of §6 is handed the pair, because averaging first puts the
tail misspecification back inside it: measured, the twenty-fold gap between the heterozygote rates
inside and outside a run collapses to 4.6-fold when the runs model is handed the marginal, and to
3.7-fold when it is handed the clean rate alone.

**Telling them apart needs the cohort**, which is why contamination is a cohort-gather parameter in
the interfaces ([`../arch/ng_step_interfaces.md:347-349`](../arch/ng_step_interfaces.md)) and not
fitted here: contaminant and index-hopped reads carry real segregating alleles and land
preferentially on polymorphic sites, and the error-like contributors do not. **The evidence for that
comparison is the census sites**
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)) — the same object, asked a
different question. **Deferred, with a home**
([`parameter_prepass.md`](parameter_prepass.md) §8).

**Decision: base qualities are not modelled at all, and reads are not filtered on them either. One
error rate per read group, estimated from the reads.** A base quality is the instrument's own claim
about itself, recalibrated by a tool that assumed a variant catalog. `ε` is estimated from the data
instead, and once it is, the qualities have nothing left to contribute. *Rejected: admit a read only
if its mean base quality clears a threshold.* Dropping the poorest reads selects on the very quantity
being estimated, so the fitted `ε` would describe the reads that survived rather than the data.

**`ε` is the error rate of *admitted* reads, and three selections upstream of this step are
correlated with it.** The same argument that rejects a base-quality threshold applies to
filters this step does not own: read filtering drops a read whose mismatch fraction exceeds 0.10
and one whose mapping quality is below 20
([`../arch/read_filtering.md`](../arch/read_filtering.md) §1). The first is a direct threshold on
the count of non-reference bases in a read — the numerator of `ε`. The second removes the mismapped
reads `ε` is defined above to contain. Both truncate the right tail of the very distribution being
fitted, and because they truncate it far out, the fit shows no strain: the estimate comes out clean,
self-consistent, and describing a read population the caller will only see if it is handed the same
filter. The mismatch filter is also reference-dependent, so it bites hardest where a sample is far
from the reference — `π_hom_alt`'s home ground.

**The third is the per-column depth cap, and it is a selection rather than a thinning.** When a
column carries more contributors than the cap allows, the walk keeps the ones admitted earliest —
that is, the reads whose alignment began furthest left, so the column sits nearer their 3' end,
where error rates are highest. It is not a random subsample and the reads it keeps are not
independent of what they show. The cap is 250 on indel-bearing columns and 8,000 elsewhere, so it
fires on the indel columns of a deep sample and never at tomato's 3 reads a site.
[`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §2.3 works through
what this step does about it — enter the locus, count it, and leave the fix upstream where the
truncation is.

**So: every emitted `ε` carries the admission policy it was fitted under**, and the caller must be
handed the same one. **Settled by:** fitting `ε` on one tomato sample at
`max_read_mismatch_fraction` of `None` and `0.10`, crossed with `min_mapq` of 0 and 20. If the four
numbers differ, `ε` is a property of the filter configuration and not only of the chemistry, and
that belongs in the parameter's provenance rather than in a footnote.

*Rejected: keep the qualities as a covariate — stratify the histogram into a few quality bins and fit
an error rate for each.* It is the obvious refinement, and an earlier draft of this section promised
it in passing without any accumulator that could deliver it. Two objections, and the second holds
however much memory we have.

- **It multiplies the accumulator that costs something.** A quality axis multiplies the cells in
  every generic histogram by the number of bins. On the read-group histogram that is nothing — four
  bins take it from 4.7 kB to about 19 kB. On the windowed histogram it is 37 MB to 150 MB per
  sample on tomato and 145 MB to 580 MB on human (§9), and that object is held once per sample in
  flight.
- **Nothing downstream could use the result.** Fitting a per-bin rate from the read-group histogram
  would itself be admissible: that object already splits a site's reads, and a rate per read survives
  splitting (§1). But a site enters every other likelihood here only through its depth and its
  alternative count, which is exactly what makes the `(depth, alt-count)` histogram a sufficient
  statistic (§4). Let `ε` vary by base quality and it stops being one — two sites with the same depth
  and the same alternative count but different quality profiles now have different likelihoods, so
  the cells no longer carry what the estimator needs. A site's reads span quality bins, so no
  re-binning of sites repairs it; the qualities would have to be kept per read, all the way to the
  caller.

The benefit against that is a rate varying across bins whose boundaries we would have to choose, and
the single-rate model is what §3's fit is written for. *An earlier version of this paragraph rejected
the covariate for want of evidence — "at 3 reads per plant there is not enough per bin" — which was
the wrong objection: that is an argument about per-site depth, and this rate is pooled over a whole
genome, where 800 Mb at 3 reads leaves hundreds of millions of base observations per bin.*

---

## 3. The read-group histogram, and the error rate

**What is accumulated: a `(depth, alt-count)` histogram per read group**, including every site with
**no** alternative reads. Those `k = 0` cells are the majority of the genome and the strongest
evidence there is about `ε`; production discards them at
[`het.rs:147-148`](../../../../src/sample_summary/het.rs) and that is
[`parameter_prepass.md`](parameter_prepass.md) §2.1's second finding.

**How `ε` is fitted — a scan on the outside, a climb on the inside.** The procedure is
[`parameter_prepass.md`](parameter_prepass.md) §3.1's, and it is not a plain grid over every
parameter at once: **step through the error rate, and at each step climb to the genotype frequencies
that fit best**, then keep the step that scored highest. Only the noise parameter is stepped through,
because with `ε` held fixed the climb over the frequencies provably cannot get stuck, so searching
them would be wasted resolution. Expectation-maximization is a reasonable default for that climb and
nothing here depends on it being EM. **On this path the scan is one-dimensional**, `ε` being the only
noise parameter — 161 steps at §3's quarter-Phred spacing, each a climb over a few hundred binned
cells. That other document carries the argument for the shape and the proof behind it.

**Only `ε` is fitted from this object.** The genotype frequencies appear in the likelihood — they
have to, since a site's alternative reads are explained jointly by error and by real variation — but
what comes out of *this* table is the error rate alone. §5 says where the genotype frequencies come
from and §5.1 how the two fits are reconciled.

---

## 4. The windowed histogram

**Decision: a `(depth, alt-count)` histogram keyed by sample × genomic window, the window fixed at
100 kb and the depth binned.** Both numbers are settled below. *Rejected:* a single genome-wide
tally — smaller, cheaper, and it forecloses the runs estimator of §6.

**A site enters this histogram once, whole** (§1) — because its genotype is one thing however many
read groups covered it. What that costs at accumulation time is nothing: summing a site's read groups
*before* binning it is exact rather than approximate, since these are raw counts at one position.

**A window's coverage is needed, and it is two different numbers.** A window with 10,000 covered
sites and one with 90,000 are not comparable as raw heterozygote counts, and coverage varies that
much across a genome. Both numbers are needed and they are not interchangeable:

- **How many loci the window holds** is the sum of its cell counts — every locus enters exactly one
  cell, including the ones with no alternative reads — so it is **derived, never stored**.
- **How many reference positions those loci covered** is *not* derivable, because a generic locus
  widened to an indel's reference span is one cell entry and several positions. It is accumulated
  alongside, as one counter per window
  ([`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §2.2).

**`F` is weighted by the second**, so that it is a fraction of the analysable genome rather than of
the locus list and a window dense in indels is not under-weighted (§6, §6.5). Neither number plays
any part in the runs model's *evidence* (§6.1) — that is the cells themselves.

**Trap: do not accumulate a count of heterozygous sites.** Counting requires calling, and calling is
the bias step 4 exists to remove ([`parameter_prepass.md`](parameter_prepass.md) §2). The heterozygote
count is derived afterwards, and softly: once the error rate and the genotype frequencies are fitted,
every site has a *probability* of being heterozygous rather than a verdict.

```text
                        π_het · (½+ε/3)^k (½−ε/3)^(n−k)
P(het | n, k) = ────────────────────────────────────────────────────────────────────────
                π_hom_ref·ε^k(1−ε)^(n−k) + π_het·(½+ε/3)^k(½−ε/3)^(n−k)
                                         + π_hom_alt·(1−ε/3)^k·(ε/3)^(n−k)

expected_hets(window) = Σ over cells   count(n, k) · P(het | n, k)
```

That is the three terms of [`parameter_prepass.md`](parameter_prepass.md) §3 divided by their sum
instead of multiplied across sites — and a window's expected heterozygote count is a sum over
**cells**, not over sites.

**This count is a readable summary, not the runs model's input.** It is what §12.3 checks against a
hard count where a hard count is safe, and it is the sensible thing to report per window. The runs
model of §6.1 does **not** consume it: it scores each state against the window's cells directly,
which keeps the homozygous-non-reference evidence that collapsing to a single heterozygote count
would discard.

Two properties earn this shape:

- **The histogram is a sufficient statistic — exactly so before binning, and near enough after.**
  `P(het | n, k)` depends on a site only through its
  depth and its alternative count, so the cell counts carry everything an estimator needs and the
  individual sites can be forgotten. Refit `ε`, change the model, swap the estimator entirely, and
  every window's heterozygote count recomputes without touching a read. **The reach of that is one
  sample's walk**, since the histogram is dropped when the sample finishes
  ([`parameter_prepass.md`](parameter_prepass.md) §1.3) — it is what lets the fits inside the walk be
  reordered or iterated freely (§5.1's alternating loop depends on it), not a promise about
  refitting a finished run.
- **An ambiguous site contributes its share and no more.** A 4-read site split 2/2 gives about half a
  heterozygote; a clean 40-read split gives essentially one; a 30-read site with no alternative gives
  essentially zero. Production's classifier discards the first, keeps the second, and never sees the
  third.

**Why the window key and not the runs themselves.** Accumulating windowed counts keeps accumulation
**associative**: a region-sharded walk merges by summing on the window key, and no shard needs to
know what its neighbour saw. Accumulating runs directly would not — a run crossing a shard boundary
would have to be stitched, which is the kind of seam that produces a bug found six months later. The
model is then a single sequential pass over some tens of thousands of windows, run once, after the
merge.

**Settled: bin the depth.** Under the rule in
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §4 — fine at the bottom,
widening geometrically, never equal-width. `P(het | n, k)` moves smoothly with depth up there, so
merging depth 100 with depth 105 shifts a window's expected heterozygote count by a fraction of one
out of hundreds, while merging depth 1 with depth 5 would discard most of what this cohort has.
Binning instead of full resolution cuts the memory nearly ninefold (§9) for an answer that moves by 0.054 rungs.

**Settled, and measured rather than argued: the ladder is exact integers to 8, then geometrically
widening bins to a cap of 124 — twenty bins in all.** An earlier version of this paragraph called
the bin count *soft*, "arithmetic rather than measurement". It is neither soft nor arithmetic: **the
edges are a correctness parameter.** Across twenty worlds the adopted ladder's asymptotic bias is
0.054 rungs of the error-rate ladder and 0.3% in each genotype frequency; the same cap at sixteen
bins is 0.55 rungs and 1.8%, and a cap of 300 at sixteen bins is 1.04 rungs and 8.0% — on data where
no site is deeper than 125, so the extra reach is spent on depths nothing occupies and paid for out
of the depths everything occupies
([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§4.3). **The bin count and the cap are one decision, not two.**

**That 0.054 is a pooled-cell figure, and a consumer reading one stored code alone pays differently —
ADDED 2026-08-13.** It was measured where many positions share a bin and the widths average out. Two
measurements since then are about the other arrangement, where a single sample's single position
carries one code. **A contamination estimator that divided an exact count of disagreeing reads by a
point read off the bin returned 2.5% contamination on a drawn panel holding none**
([`../reports/contamination_floor_and_duplicated_class_2026-08-13.md`](../reports/contamination_floor_and_duplicated_class_2026-08-13.md));
and **a copy-number discriminator reading the stored code instead of the exact count loses 11% of its
enrichment at 13× depth and 37% at 29×**
([`../reports/locus_depth_vs_window_2026-08-13.md`](../reports/locus_depth_vs_window_2026-08-13.md)).
Neither is an argument against the ladder — the histograms this section sizes are the pooled
arrangement, and there 0.054 stands. **It is a warning to whoever reads a single stored code: sum over
the depths the code stands for, rather than taking a value from inside the range.**

**Where a ladder can hurt is 10 to 30 reads a site.** At tomato's 3 reads, 97 sites in 100 sit at
depth 6 or below and are never binned at all; at 60 reads the genotype is certain whatever the exact
depth. So **check any change to the ladder in that band and against both consumers** — the runs of
§6 and the sample rates of §5 — because a check run at tomato's own depth would pass anything. The
homozygous-non-reference rate is the most sensitive of the three outputs, about 1.5 times the
heterozygosity error in every ladder tried, which matters most for the sample §5 was written about.

**Settled: the window is 100 kb, fixed, everywhere.** The window size is not a memory knob — it is
the resolution at which a run can be seen at all, and a run shorter than a few windows is invisible.
At 1 Mb, the research note's figure inherited from ROHan, the shortest resolvable run is several
megabases; at 100 kb it is a few hundred kilobases. Three things decided this:

- **It is affordable.** 37 MB per sample on tomato and 145 MB on human, against 3.7 MB and 14.5 MB at
  1 Mb (§9) — and this object lives only for the duration of one sample's walk
  ([`parameter_prepass.md`](parameter_prepass.md) §1.3).
- **The evidence per window is ample — in fact more than ample.** At tomato's roughly one
  heterozygote per kilobase a 100 kb window carries about 100 expected heterozygotes outside a run
  against near zero inside, and that separation is so wide that **every window is settled by its own
  reads and the chain adds nothing** (§6.1). The arithmetic gets thin around 10 kb.
- **One constant beats a knob.** *Rejected: expose the window size, with a per-organism default.* The
  organisms do differ — a selfing crop's runs reach tens of megabases while an outbred human genome
  carries much shorter autozygous segments — but 100 kb is adequate for both, and a window size is
  not a quantity a user is in a position to choose. An unsettable knob is worse than a constant,
  because it invites a wrong answer and offers no way to recognise one.

**What 100 kb misses, stated rather than left to be discovered.** A run shorter than a few windows
is invisible, so the shortest run this resolves is about 300 kb. That is far below anything a
selfing crop has — a tomato landrace is homozygous over tens of megabases — and far below a
consanguineous human's segments, which run from 5 to 50 Mb. **What it misses is the background
autozygosity of an ordinary outbred human**, most of which sits in segments of 0.1 to 2 Mb. Two
reasons that is the right trade rather than a gap to close:

- **What is missed contributes an `F` of order 0.01**, and a genotype prior mixing
  `F·π_i + (1−F)·π_i^ploidy` moves imperceptibly between `F` = 0.01 and `F` = 0.02. The segments
  100 kb cannot see are the ones the consumer cannot feel.
- **The estimator's own resolution is the same size.** On a genome with no runs at all it returns
  about 0.01 at tomato's 8,004 windows and 0.003 at a human genome's 31,000 (§6.1), so a human
  background `F` is within a factor of a few of the noise however fine the window. Cutting the
  window to 10 kb recovers `F` exactly and costs ten times the windows — 1.45 GB per human sample
  against 145 MB (§9), on the most expensive accumulator step 4 has. **Not worth paying for an `F`
  no caller can feel**, and the measurement is there if that judgement ever needs revisiting
  ([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md) §3.5).

**A check worth running anyway, though it no longer gates anything:** the length distribution of the
runs actually present in the tomato cohort, measurable today. If runs there turn out to be shorter
than a few hundred kilobases, 100 kb is the number to revisit.

**This object is built on every run, including one that never computes `F`.** It is the only
accumulator that keeps whole sites, and §5 needs whole sites. The window key is what inbreeding
costs; the object's existence is not. §6.4 works that through for the case where `F` is supplied on
the command line rather than fitted.

---

## 5. The sample's rates — heterozygosity and distance from the reference

Two numbers, both fitted from the windowed histogram summed over its windows: the **observed
heterozygosity** — how often an individual's two copies of a site differ, `Hobs` in the population
genetics and `π_het` in the likelihood — and the **homozygous-non-reference rate** `π_hom_alt`, how
often *both* copies differ from the reference.

**The second is not a leftover of the first; it measures something else entirely.** How often an
individual carries a non-reference allele at all is *heterozygosity + the homozygous-non-reference
rate*, and that quantity belongs to the **pair** (individual, reference), not to the individual:
swap in a different accession as the reference and it changes, while heterozygosity and inbreeding
do not. The two also come apart in the direction that matters here. A selfing landrace far from the
reference accession is **mostly homozygous and mostly non-reference at the same time** — low
heterozygosity, high homozygous-non-reference rate. A caller whose prior assumes "non-reference
implies rare" is wrong on exactly that sample, and tomato's reference is one cultivated accession.

**Two alternatives were rejected on the way to §1's split, and both are tempting enough to record.**

*Rejected — fit all three rates per read group, then average the genotype frequencies to the sample.*
This was an earlier decision here, and it is wrong twice: it fits a quantity the individual has only
one of once per library, from a table whose sites have been split, and then supplies no principled
weight for the average. The averaging step was the tell.

*Rejected — one fit per sample on the summed read-group histograms.* The sum §1 rules out, and it
also forces a single error rate across libraries whose chemistry differs, which is the whole reason
for the read-group grain.

### 5.1 The two fits are coupled — CLOSED: alternate, and it is consistent

A higher error rate explains the same alternative reads as less real variation, so `ε` and the
genotype frequencies trade off inside one likelihood — but they are read off two different tables.

**Decision: alternate.** Hold the genotype frequencies where they are and fit each read group's `ε`
from the read-group table; then hold the rates and climb to the frequencies from the whole-sample
table; repeat until the fitted rung stops moving.

**That is a fixed point of two estimating equations rather than a climb on one objective, so it
needed checking rather than assuming — and it checks out.** Weighting both tables by their exact
probabilities under a known truth and starting deliberately away from it — every error rate at three
times the truth, every frequency at half — the loop converges to the truth in all 25 worlds tried:
zero rungs on every rate and zero relative error on both frequencies, across error-rate ratios of 1,
4 and 10, depths 3 to 20, even and 90/10 splits, two libraries and four
([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md) §2.6).

**The reason it works is worth keeping, because it also settles a nearby worry.** A read-group
entry's own distribution is correctly specified — the genotype is still drawn once for the site and
still enters through the same mixture — so each block's score is unbiased at the truth and the truth
is a fixed point of the pair. The same fact means **splitting a site between two entries costs
precision and not correctness**: fitting *both* the rates and the frequencies from the read-group
table alone is also exactly unbiased. §1's objection to a windowed histogram keyed per read group is
about what that would do to a *per-window* statistic, not about splitting as such.

**Convergence is linear**, as alternating schemes are, so a tolerance fine enough to be interesting
is far finer than the answer needs — worlds that ran 200 iterations without meeting a movement
tolerance of 10⁻¹² were already at the truth to better than a thousandth of a rung. Stop on the
ladder rung ceasing to move, cap the loop, keep the best-scoring iterate
([`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §5.2).

**AMENDED 2026-08-10 — the second class of site sits outside this alternation, and every rung of
the `ε` scan is scored with it.** The alternation above is unchanged: two blocks, two tables. What
changed is that the pair `(w, ε_noisy)` of §2.1 is fitted around it, by holding `ε_noisy` at each
rung of the ladder in turn while the rates, the frequencies and the share settle, and keeping the
best-scoring rung. **Leaving the pair out of the `ε` scan was a defect and not a simplification**:
a candidate rate scored under the one-class rule is being asked to explain a table whose tail
belongs to the other class, so the scan returns the tail-inflated rate whatever pair sits beside
it — and since that scan is where this step's rate comes from, the rate then never moves. Measured
on a table generated at HG002's own parameters it came back **three rungs high**, with the
generating parameters scoring 351 nats better
([`../reports/implementations/ng_noise_model_extension_n5_fix_2026-08-10.md`](../../reports/implementations/ng_noise_model_extension_n5_fix_2026-08-10.md)).

**The two classes are not identified without an ordering constraint**, and it has the same shape as
§6.1's `h << Hout`. Swapping the labels of a mixture describes the same distribution, so a "noisy"
class *finer* than every library's clean rate and holding most of the genome scores identically to
the reading anyone means — and a sample emits the two rates weighted by the share, so accepting it
reports a rate an order of magnitude below what the reads support. The runs model resolves its
version by relabelling after the fit; that is not available here, because the clean rate is one per
read group and the noisy rate is one per sample, so there is no single rate to swap with. **A
candidate rate at or below the coarsest fitted clean rate is refused instead.**

**And the 25 worlds above were re-fitted through the extended model and are unchanged** — 0.000
rungs and 0.000% throughout. That is a regression check and not evidence the extension works: those
worlds carry one error rate, so there is no second class in them to find. **What is evidence** is a
further six worlds that do have one, each built as one library and as two and fitted through both
cell keys, which recover the share, both rates and both frequencies **exactly**
([`../reports/implementations/ng_noise_model_extension_n4_2026-08-10.md`](../../reports/implementations/ng_noise_model_extension_n4_2026-08-10.md)).

**Ask first how much this can matter.** Both cohorts in hand are single-library throughout
([`parameter_prepass.md`](parameter_prepass.md) §5), so the coupling bites only on multi-library
data — 157 of 1,707 samples in the tomato archive survey carry more than one library.

---

## 6. Inbreeding

**`F` is the probability that an individual's two copies of a stretch of DNA came from the same
ancestral copy** — one grandparent's chromosome reaching it down both parental lines. Where that
happens the individual has no choice but to be homozygous. The caller uses it as exactly that
probability: its genotype prior mixes the two cases, `F·π_i + (1−F)·π_i^ploidy`
([`../../specs/ssr_cohort_mark2.md`](../../specs/ssr_cohort_mark2.md) line 286).

**It is not chemistry.** Inbreeding belongs to the individual, so `F` is estimated **per sample**,
pooling that sample's read groups.

**`F` is a fraction of the analysable genome, not of the reference.** Each window's probability of
being inside a run is weighted by how many **reference positions** that window actually covered (§4)
— so a window where reads were scarce counts for less. That is the quantity the caller wants: it
applies `F` at the loci it is calling, and those are covered loci. Weighting every window equally
instead would let an unmappable stretch no caller will ever visit pull the estimate around. The
denominator therefore differs slightly between samples of different coverage, which is harmless
unless coverage is correlated with autozygosity — and the correlation that does exist runs the
protective way, since the under-covered regions are the repetitive ones where collapsed paralogs
manufacture false heterozygotes.

**`F` is computed here, not deferred.** The cohort's diversity divides by `1 − F`
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3), so without `F` there is no
diversity, no frequency-spectrum scale, and half of what
[`parameter_prepass.md`](parameter_prepass.md) §4.1's comparison is meant to measure cannot be
produced at all.

### 6.1 The two candidate estimators want different things

- **The ratio, `1 − Hobs/Hexp`.** Count how often this individual is heterozygous; compare it to how
  often random mating would, given the population's diversity. The observed half needs no counting —
  the heterozygosity of §5 *is* it. **No window key. Needs a cohort-level diversity.**
- **Runs of homozygosity.** Inbreeding does not arrive spread evenly along a chromosome. A shared
  great-grandparent hands down whole intact stretches, tens of megabases long, and inside one both
  copies are the same ancestral DNA — so heterozygous sites are nearly absent, while between them the
  genome looks ordinarily outbred. Walk the genome in windows and ask of each: **inside a run, or
  outside?** `F` is the fraction of the genome inside a run. **Needs a window key. Needs only this
  sample.**

**How a run is found.** A two-state hidden Markov model — *inside a run*, *outside a run*, the state
never observed and inferred only from the data. Each state carries its own **genotype frequencies**,
and what changes between the two states is worth being precise about, because it is not only the
heterozygotes.

Inside a run the two copies are the same ancestral copy, so the genotype is one allele draw rather
than two: heterozygotes are nearly absent, and the mass that would have been theirs sits in the
homozygous classes instead. **Heterozygotes do not simply vanish — their share moves to the
homozygous-non-reference class**, so those cells carry a second, independent signal for locating a
run, and a model tracking only a heterozygosity rate would discard it. That matters most for the
sample §5 was written about: a selfing landrace far from the reference is mostly inside runs, and its
high homozygous-non-reference rate is partly the run structure itself.

**Decision: each state carries its own three genotype frequencies, fitted freely.** Four free
numbers — outside `(1 − Hout − Aout, Hout, Aout)`, inside `(1 − h − Ain, h, Ain)` — with the
identifying constraint that the inside state's heterozygote share `h` sits well below the outside
state's `Hout`. **That ordering is what makes the model about autozygosity rather than about
bimodality**; nothing else ties the two states together.

*Rejected, and it is worth recording precisely because it looks right: tie both states to one
allele frequency `f`, outside `((1−f)², 2f(1−f), f²)` and inside `(1−f−h, h, f)`.* It is
`bcftools roh`'s parameterisation, and it is correct **there** because that tool reads `f` **per
site** from an allele-frequency tag. With one genome-wide `f` it is not merely approximate — it
forces the answer. Take the model's own identities, use `f = π_hom_alt + π_het/2` (the alternative
allele's frequency in a diploid), and `F` falls out with no run structure left in it:

```text
F  =  (π_hom_alt − f²) / (f − f²)
```

On GIAB HG002 — 2,427,887 heterozygous and 1,594,600 homozygous-non-reference sites over its
confident regions — that is **F = 0.57 for a non-consanguineous individual whose true value is
zero**, and 0.91 for a tomato landrace far from the reference. The cause is that real genomes carry
a frequency spectrum: heterozygotes come from variants at intermediate frequency, homozygous-
non-reference sites from near-fixed differences to the reference accession. One `f` fitted to the
heterozygotes gives `f² = 1.2e-6`, predicting **580** homozygous-non-reference sites where HG002 has
1,594,600 — short by a factor of 2,750. The only place the model can find them is the inside state,
whose weight is `f` and not `f²`, so it puts the genome inside runs. It is not a local optimum the
fit could escape: setting `F = 0` and raising `f` to `√π_hom_alt` predicts 48 heterozygotes per
kilobase against the 1 observed, and loses by an enormous margin.

**This is the same tie §11.4 rejects for the genome-wide fit**, and for the same reason — it presumes
the answer. §6.1 imposing it while §11.4 forbids it was a real contradiction, not a nuance.

*Also rejected: leave the inside heterozygote share at exactly zero.* Then a heterozygote inside a
run is impossible, and one false one — a collapsed paralog, which this project has measured on this
cohort — costs the inside state the whole ratio between the heterozygote and homozygous-reference
terms at that site: about 125 for a site showing one alternative read of three at `ε = 0.001`,
1.3 × 10⁵ for two of three, and past 10³⁵ for fifteen of thirty. The model breaks the run rather
than accept it, so `F` deflates in proportion to each sample's artifact density: the pathology the
tomato baseline recorded, reproduced by the estimator meant to diagnose it. `h` fitted rather than
fixed is what absorbs it, and
[`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§3.3 measures that it does.

**What freeing the states costs, stated rather than hidden.** Per-site allele frequencies would let
each site discriminate according to how informative it is — a heterozygote at a common variant is
strong evidence against a run, one at a private variant much weaker. Pooled into cells, every site
contributes the average. The estimate stays consistent, because a state's fitted frequency vector is
the correct average over the sites assigned to it; what is lost is **power, not correctness**, so
more windows are needed for the same confidence.

**And the residual risk runs the opposite way to the failure above, and it is now measured.** A
two-state model can always improve its likelihood by splitting on window-to-window noise, so a
genuinely outbred genome comes back with a small non-zero `F` rather than exactly zero. **That
number is the estimator's resolution, and it is a function of how many windows the genome holds**
([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§3.6): on a genome generated with no runs at all it comes back at about **0.01 at tomato's 8,004
windows** and **0.003 at a human genome's 31,000**, shrinking as windows accumulate. So an `F` below
that means *nothing detected*, not a small autozygous fraction, and §6.5 says what has to be emitted
for a reader to tell the two apart. It also means a run over a few hundred windows cannot estimate
`F` at all — at 1,200 windows a genome with no runs returned 0.84 on one seed of eight.

**The floor is also what makes the robustness claim below true rather than contradicted, and that
claim now has a measurement behind it.** With `h` fitted, a uniform floor of false heterozygotes
lifts `h` and the outside rate together and cancels out of the gap between them. Adding spurious
heterozygotes at up to five times the real rate of one per kilobase moves `F` not at all — both
fitted rates rise with the floor and the gap survives — while the whole-genome heterozygosity the
ratio estimator of §6.3 would read inflates eight-fold (research note §3.3). *Soft: whether one `h`
per sample suffices*, or whether it must vary with local mappability; the per-sample constant is the
cheap version and the one to measure first.

**Where the four frequencies come from.** All of them are fitted per sample from this same windowed
histogram, inside the per-sample walk — nothing is borrowed from the cohort, and nothing needs the
frequency spectrum, which is a gather output ([`parameter_prepass.md`](parameter_prepass.md) §1) and
not available here. Freeing the states is precisely what removes the need for it: a state absorbs
the spectrum's effect into its own fitted vector instead of trying to derive it from one number.

`bcftools roh` is built on exactly this contrast, and its emissions are the reference form to read
([`bcftools/vcfroh.c:476-499`](../../../../bcftools/vcfroh.c)): a mixture over the three genotypes
weighted by the state's frequencies, computed from genotype **likelihoods** and never from calls.

**What the model reads is the cells, not a heterozygote count.** Given the state, sites are
independent, and each site enters only through its depth and alternative count — so a window's
log-likelihood under state `s` is a sum over that window's cells:

```text
log P(window | s) = Σ over cells  count(n,k) · log [ π_hom_ref,s · ε^k (1−ε)^(n−k)
                                                   + π_het,s     · (½+ε/3)^k (½−ε/3)^(n−k)
                                                   + π_hom_alt,s · (1−ε/3)^k (ε/3)^(n−k) ]
```

which is [`parameter_prepass.md`](parameter_prepass.md) §3's likelihood with the state's frequencies
substituted. The factorisation is exact — the product over a window's sites is the product over its
cells raised to their counts — for cells that hold one exact depth. **Binning the depth makes it an
approximation**, because sites sharing a cell no longer share a likelihood; scoring each cell at its
own mean depth (§9) is what keeps the error small.

**How small is now measured, and the answer is that `F` does not move at all**
([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§4.4). The same drawn genome, refitted with only the keying and scoring changed, returns the same
`F` to four decimal places at 3, 20 and 60 reads a site, under every ladder tried, whether the cell
mean is taken per window or over the whole sample, and with window coverage flat or varying
threefold. At 60 reads that is 15,050 cells a window reduced to 458 for an identical answer. What
does move is the two states' fitted heterozygote rates, by up to 8% on the inside rate and 5% on the
outside one — and `F` reads only the gap between them, which is twentyfold and survives. **No
per-window site count appears here, and that
is deliberate** — dividing by one would make a half-covered window look autozygous merely for being
thinly read.

**Transitions are set per base, not per window — and there are two of them.** A two-state chain has
a rate for entering a run and a rate for leaving it; `bcftools roh` holds both as named parameters
([`bcftools/vcfroh.c:486-487`](../../../../bcftools/vcfroh.c), `tAZ = P(AZ|HW)` and
`tHW = P(HW|AZ)`). An expected run length fixes only the second. **Fixing both a priori would set
the chain's stationary inside-probability, `tAZ/(tAZ + tHW)`, and so assume the shape of the answer;
fixing one still shrinks every sample's `F` toward whatever the other implies wherever a window's
own evidence is thin.**

**That stationary probability is not `F`, and an earlier version of this section said it was.** It
is a property of the fitted *model* — what fraction of an infinitely long genome from a chain with
those rates would lie inside a run. `F` is a property of the *data*: what fraction of **this**
genome does. For a finite genome the two differ by ordinary sampling, measured at 11% relative at
`F` = 0.05 and 3.5% at `F` = 0.30, and only the second recovers the genome's realised autozygous
fraction (research note §3.7). The caller wants the second, because its prior asks whether this
individual's two copies at this locus descend from one ancestral copy — not what a hypothetical
genome from the same pedigree would average. §6.5 is where `F` is defined.

**Decision: fit both, per sample, inside the same loop that fits the state frequencies.** That makes
this Baum–Welch rather than a forward–backward pass at fixed transitions, so it inherits the caveat
of every expectation-maximization: it climbs to a stationary point, not provably to the maximum, and
the surface here is not concave. It therefore needs the same treatment as the coupled fit — an
iteration cap, the best-scoring iterate kept, and the termination reported
([`../arch/parameter_prepass_generic.md`](../arch/parameter_prepass_generic.md) §5.3). *Rejected:
hold both transitions at constants.* It is simpler and is what `bcftools roh` does by default — but
its defaults imply a stationary autozygous fraction of 0.93, and shipping an assumed `F` inside a
tool whose output *is* `F` is the circularity §6.3 exists to prevent, one level lower down.

**The per-base part still stands, and matters more now.** Whatever the fitted rates, they are held
per base and converted to a per-window probability at the window size; holding a per-window
probability instead would tie the model's expected run length to the window, so changing 100 kb to
1 Mb would grow the runs it expects tenfold with nothing to announce it. `bcftools roh` scales its
transitions by distance for the same reason
([`bcftools/vcfroh.c:452-472`](../../../../bcftools/vcfroh.c)).

**Windows with no data are not gaps in the chain.** The accumulator holds only windows that received
a site, so iterating it would place an unmappable 3 Mb block as a single window-to-window step and,
worse, follow the last window of one contig with the first of the next. The chain is walked over the
**full window range of each contig**, absent windows included as empty — an empty window emits
nothing and the transition still advances — and forward–backward restarts at every contig boundary.

The forward–backward algorithm then gives each window the probability it is inside a run, using every
window on the chromosome rather than each alone. `F` is those probabilities weighted by each window's
covered-site count (§4, and §6's opening).

**How much the neighbours contribute is a function of the window size, and at 100 kb the answer is
nothing** ([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
§3.5). A 100 kb window at three reads a site holds about 100,000 covered sites, and **not one of
8,004 such windows is left undecided** — every one is settled by its own evidence, and shuffling the
whole genome into a random window order changes `F` by zero. So the transitions, the contig
restarts and the empty-window rule above are **insurance rather than working machinery at this
grain**: they cost little, and they are what would carry an estimate through a window too thin to
settle itself. The share of undecided windows reaches 5% at 10 kb, 48% at 1 kb and 84% at 300 sites,
where shuffling moves `F` by 0.25. An earlier version of this paragraph claimed the pooling *"lets a
lone quiet window between noisy ones read as a fluctuation rather than a 100 kb run"*; at 100 kb
there are no such windows.

**Both states' frequencies are fitted by the model itself, and that is where the robustness comes
from.** A uniform floor of false heterozygotes — mismapping, collapsed paralogs — lifts the inside
rate and the outside rate together, and it is the *gap* between them that decides each window's
state. In the ratio estimator the same floor goes straight into `Hobs` with nothing to cancel it.

**One identity survives the untie, and it ties this section back to §5.** Averaged over the genome
the two states mix in proportion `F`, whatever each state's frequencies turn out to be:

```text
π_het      =  F·h    +  (1 − F)·Hout
π_hom_alt  =  F·Ain  +  (1 − F)·Aout
```

With `h` near zero the first reduces to `Hobs = Hout·(1 − F)`, which is §6.3's ratio estimator
recovered — so the two routes to `F` still have something to say to each other, and §6.3's
cross-check still works. **What no longer appears is `f`**: nothing here assumes the genome has one
allele frequency, which is what made the tied version return 0.57 on an outbred genome. **Worth
emitting as a residual**: the gap between these two identities at the fitted values and the
separately-fitted `π_het` and `π_hom_alt` of §5 is one division per sample, and it is the only thing
that ties the two fits together at all.

### 6.2 Decision: the runs estimator is the one a caller reads

Four reasons, in increasing weight:

- **It is the quantity the consumer asks for.** The genotype prior is a mixture over *whether the
  two alleles are one ancestral copy*, which is realized autozygosity. `1 − Hobs/Hexp` measures a
  deviation from Hardy-Weinberg proportions instead, and so absorbs population structure: a cohort
  that is really two subpopulations looks homozygote-excessive for reasons no individual's parents
  caused. Mark-2 warns about this and does not correct it (line 286).
- **It is robust to a uniform floor of false heterozygotes** (§6.1).
- **It separates artifact from biology** — whether a sample's excess heterozygosity is uniform (an
  artifact) or segmental (a real outcross), the pathology the tomato baseline recorded and that no
  read-quality heuristic could reach.
- **It is the only one that leaves the cohort's diversity estimable at all** (§6.3). This is the
  decisive one, because it is about correctness rather than preference.

### 6.3 The circularity, and the ratio as a labelled diagnostic

The cohort's diversity is estimated as `Hexp = mean over samples of Hobs / (1 − F)`
([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3). That is §6.1's ratio estimator
rearranged, which creates a trap:

> **Do not take `F` from the ratio estimator and then compute the cohort's diversity from it.** That
> is circular: the ratio estimator *needs* a diversity to produce `F`, so feeding its `F` back in
> returns whatever was assumed. The runs estimator has no such problem, because it reads `F` off the
> **genomic distribution** of heterozygosity and never needs a population expectation.

**This is the one place a downstream step constrains a choice made here**, and it is an argument
about correctness rather than robustness: the runs estimator is what makes the cohort's diversity
estimable at all.

**The ratio estimator is still worth computing, in the one order that is not circular.** Once the
runs estimator has produced `F` and the gather has turned it into a diversity, feed that diversity
back through `F = 1 − Hobs/Hexp` and see whether the same `F` comes out.

- **It costs nothing.** Both inputs already exist by then; the arithmetic is one division per sample.
- **It is a real check, not a formality.** The two routes read the data differently — one from the
  *distribution of heterozygosity along the genome*, the other from its *total* against a population
  expectation — so agreement is evidence the model holds, and disagreement localises the problem.
  A cohort with population structure should disagree in a predictable direction, since `1 − Hobs/Hexp`
  absorbs Wahlund effects that autozygosity does not (§11.2).
- **What it must not become is an input.** The moment its `F` is fed back into the diversity it
  helped produce, the loop closes and both numbers become whatever was assumed. Emit it as a
  diagnostic, clearly labelled, never as the `F` a caller consumes.

### 6.4 When `F` is supplied rather than fitted

**A run may be handed `F` instead of estimating it** — `F = 0` for an outcrossing species is the
obvious case, and the value would come from the command line. What that changes is smaller than it
looks, and the thing it does *not* change is worth stating first.

- **The window key is dropped.** No runs pass means no need to see the genome in order, so the
  windowed histogram collapses to one histogram per sample: 4.7 kB instead of 37 MB on tomato (§9).
- **The object is not dropped.** §5's two rates are counted per **site**, and the read-group
  histogram has split every site covered by two libraries into two shallower ones. A sample-keyed
  whole-site histogram is still required, window key or not. Only when the sample *also* has a single
  read group do the two coincide and the generic accumulator set collapse to one histogram (§1).

**A supplied `F` reopens §6.3's circularity from outside the program.** That trap is fed by
`Hexp = Hobs/(1 − F)`: an `F` that itself came from a ratio against some diversity returns that
diversity when fed back in. Within one run the order of operations prevents it; a number typed on a
command line carries no provenance. `--inbreeding 0` for an outcrosser is safe, an `F` copied from a
previous run's ratio diagnostic (§6.3) is not, and nothing in the program can tell them apart.

**So record a supplied `F` as supplied.** That is a fourth provenance state alongside fitted,
borrowed and defaulted ([`parameter_prepass.md`](parameter_prepass.md) §6) — more authoritative than
a default, less checkable than a fit, and the only one whose correctness rests on the person who
typed it.

### 6.5 What `F` is, how the fit must be started, and what has to be emitted beside it

**`F` is the coverage-weighted posterior occupancy**: each window's forward–backward probability of
being inside a run, weighted by how many reference positions that window covered, summed over the
genome. Not the transition rates' ratio (§6.1). It recovers a drawn genome's realised autozygous
fraction to four decimal places at every level tried, from 0.05 to 0.60
([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md) §3.2).

**The fit must be started from several points, and they must disagree about how far apart the two
states are — not only about `F`.** This is the one place the estimator produces a confident wrong
number, and it is worth being exact about how.

Baum–Welch climbs to a stationary point. A start that guesses the inside state's heterozygote rate
far below the truth fits every window to the outside state on its first pass, empties the inside
state, and drives the rate of entering a run to zero. **On a genome with 26% of its length in runs
and a floor of spurious heterozygotes three times the real rate, that returns `F` = 0.0000, reports
convergence, and gives no other signal** — because every start made the same wrong guess, so keeping
the best-scoring one had nothing better to pick. Starts spanning the separation return 0.2634 on the
same data, against a realised 0.2629 — which is where the 26% above comes from, so the two can be
checked against each other (research note §3.4). **At minimum: the inside heterozygote rate started at 1/20, 1/3 and
3/4 of the outside one, crossed with a few implied `F`.** Nine such starts cost seconds on 8,000
windows and are the difference between an answer and a plausible zero.

**Three things go out with `F`, and none of them is decoration.**

- **The fitted separation** — both states' three frequencies, and the ratio between their
  heterozygote rates. The failure above leaves the inside state's rate at *exactly* its starting
  value, because nothing was ever assigned to it.
- **The spread across starting points** — the best and second-best `F` and their scores. This is
  the only thing that separates a genuine `F` = 0 from a search that never found the second state:
  both leave an empty inside state, and only the scores say whether a better answer was looked for
  and rejected. **A run where every start returned the same `F` at the same score has not measured
  zero autozygosity; it has failed to find anything, and must say so.**
- **The resolution** — the noise floor at this run's window count (§6.1): about 0.01 at 8,000
  windows, 0.003 at 31,000. An `F` below the floor is *nothing detected*.

**Below a few thousand windows `F` is not estimable and must not be emitted.** At 1,200 windows a
genome generated with no runs at all returned `F` averaging 0.23, and 0.84 on one seed of eight. No
real run is that small — a tomato genome is 8,004 windows and a human 31,000 — but development
fixtures and region-restricted runs are, and a number produced there would look like any other.

---

## 7. Ploidy

[`parameter_prepass.md`](parameter_prepass.md) §3 makes the likelihood work at any ploidy, and both
histograms are ploidy-free, so nothing in the accumulation blocks a polyploid. Two definitions do,
and they belong to the estimators rather than to the accumulators:

- **Heterozygosity stops being a yes-or-no.** A tetraploid site can carry one, two or three
  alternative copies, and calling all of them "heterozygous" throws away the dosage. The natural
  replacement is **gene diversity** — the chance that two copies drawn at random from the individual
  differ — which reduces to heterozygosity at `P = 2` and is what the soft count of §4 should
  accumulate toward for `P > 2`.
- **There is no single inbreeding coefficient above diploidy.** Polyploids need several
  identity-by-descent coefficients, and autopolyploids add double reduction, so §6.1's two-state
  picture is a diploid simplification rather than a general model.

**Deferred, with a home (§10).** The accumulators serve whatever replaces them, which is the point of
keeping them model-free.

---

## 8. Two assumptions inside the `½`

**The het model's `½` is optimistic, in two independent ways.** The `(½)^n` term makes two
assumptions and both are wrong:

1. **That both alleles are sampled equally.** Reads carrying the alternative allele map slightly
   less often than reference-carrying ones, so a true het sits nearer 0.47–0.49 than 0.50. Bryc et
   al. fit exactly this term; it replaces `½` with one fitted per-read-group constant.

   **What the `½` costs is now measured, and it is a shallow-sample problem**
   ([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
   §5). Generating at a true balance `b` and fitting with the model that assumes a half, at
   tomato's one heterozygote per kilobase:

   | | at 3 reads a site | at 10 | at 20 and above |
   |---|---|---|---|
   | `b` = 0.47 | `Hobs` **−4.5%**, `ε` 0.16 rungs | −1.4%, 0.03 rungs | −0.2%, 0.002 rungs |
   | `b` = 0.44 | `Hobs` **−9.8%**, `ε` 0.32 rungs | −3.5%, 0.08 rungs | −0.5%, 0.007 rungs |

   The homozygous-non-reference rate never moves more than 0.4% anywhere. **The mechanism is class
   confusion, not arithmetic:** at 3 reads a heterozygote usually shows one or two of three, and
   lowering `b` pushes some of them to *none of three*, where they are indistinguishable from a
   homozygous-reference site and are simply lost. At 20 reads a heterozygote shows nine or ten of
   twenty and never zero, no site changes class, and the misfit costs nothing visible.

   *Leaning: adopt, and the reason has narrowed.* One parameter, and the grid is already there —
   but what it buys is heterozygosity on **shallow** cohorts specifically. A cohort at 20 reads a
   site or more may keep the `½` with a clear conscience. Tomato's 3 reads is exactly the case
   that cannot. Whether `b` really is 0.47–0.49 on our data is still §11.3's question.
2. **That errors are independent.** They are not, and the research note rates this the *larger* of
   the two effects: the false hets that matter — collapsed paralogs, mismapped repeat copies — are
   systematic, so 40 identical reads are not 40 independent observations. **Deferred, with a home
   (§10)**, and worth knowing that it is the bigger of the two while reading the smaller one's
   adoption above.

---

## 9. Memory

**The windowed histogram is the expensive accumulator of step 4.** Its key carries the genomic
window: 8,000 windows on tomato at 100 kb, about 31,000 on human, **per sample**. The sizes below
are **arithmetic, not measurement** — dense storage of every `(depth, alt-count)` cell up to a depth
cap, **eight bytes a cell**: four for the site count and four for the sum of the exact depths that
landed in it, which §4's scoring needs per cell rather than per bin.
[`parameter_prepass.md`](parameter_prepass.md) §10.6 replaces the arithmetic with a measurement.

| depth resolution | cells/window | per window | **per sample** (tomato / human) | if fifty were held at once |
|---|---:|---:|---:|---:|
| full, depth ≤ 100 | ~5,150 | 41 kB | 330 MB / 1.3 GB | 16 GB / 64 GB |
| **20 depth bins — adopted** | 583 | 4.7 kB | **37 MB / 145 MB** | 1.9 GB / 7.2 GB |

The adopted row is the ladder §4 fixes — exact integers to 8, then eleven widening bins to a cap of
124 — whose cell count is `Σ (bin's top depth + 1)` = 583, since a bin's row must be as wide as its
deepest site's alternative count. An earlier version of this table carried 465 cells and 30 MB, from
before the ladder was measured rather than assumed.

**Two things are not in that table, both small and both only on some runs.** A sample with two or
more libraries also carries the attributed cells — the cells that record which library each of a
site's four-or-fewer alternative reads came from (§1) — whose size is the number of *observed*
attributions rather than the number possible. At two libraries and tomato's 3 reads that is about
40% more cells than the pooled key alone, and at 20 reads about 15% more; a single-library sample
never builds them at all, since a lone library's attribution says nothing the pooled key does not.
And a genome with more than one ploidy carries one table per ploidy present, which multiplies
nothing in practice: a haploid sex chromosome is a few percent of the windows, not a second copy of
them.

**Read the per-sample column, not the last one.** This histogram is reduced to its parameters at the
end of its own sample's walk and then dropped
([`parameter_prepass.md`](parameter_prepass.md) §1.3), so a fifty-sample run holds one per sample
*in flight*, not fifty. The last column is what a design that kept them all would cost, and it is
here only because that design is the easy mistake to make. **How many are in flight is the multiplier
that decides peak memory, and nobody has chosen it yet** ([`parameter_prepass.md`](parameter_prepass.md)
§6): at eight concurrent samples the adopted row is 300 MB on tomato and 1.2 GB on human.

That ninefold gap is the argument for binning depth, and it is a ninefold gap per sample either
way. It is not a budget judgement — the full-resolution table costs nearly nine times as much for an
answer §4 shows moves by 0.054 rungs. **The unit is the sample, not the read group** — a consequence of §4
keying by sample so a site enters once — so a three-library sample costs the same here as a
one-library sample.

**What inbreeding costs is the window key, not the object.** §4 has already said the histogram is
built whether or not `F` is computed; the table above says what the *windows* add. Drop inbreeding —
or supply `F` outright (§6.4) — and it collapses to **one histogram per sample**, a single cell of
the "per window" column: 4.7 kB rather than 37 MB, a factor of 8,000. The windows are what `F` costs;
the cells are what step 4 costs.

**The read-group histogram is free by comparison**: a few hundred cells, about 4.7 kB per read
group — one part in eight thousand of the windowed object (§1).

**Concurrency.** Both accumulators are keyed and therefore associative: a region-sharded walk merges
by summing, with no communication between shards. The only sequential step is the runs pass itself,
which walks 8,000 windows of one tomato sample in order and is irrelevant to the walk's cost.

**Determinism.** The soft heterozygote count is a sum over cells; fix the summation order so it does
not vary with thread count.

---

## 10. Deferred, with a recommended home

- **Inbreeding and diversity above diploidy** (§7) — several identity-by-descent coefficients instead
  of one `F`, double reduction in autopolyploids, heterozygosity replaced by gene diversity.
  **Home:** a spec of its own, which should pick a definition that degrades to the diploid one at
  `P = 2` rather than maintaining two code paths. Nothing in the walk changes when it lands.
- **Overdispersion — a beta-binomial in place of the binomial** (§8, second assumption): one extra parameter for the
  correlated errors the independence assumption ignores. This project has been bitten by the same
  effect before, which is why it is recorded rather than dismissed. **Home:** the same spec that
  adopts or rejects the reference-bias term, since both change the same `½`.
- **Splitting contamination out of `ε`** (§2). **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §5.

---

## 11. Open questions

1. **How are the error rate and the genotype frequencies fitted, now that they come from two
   different tables?** — **CLOSED: alternate.** §5.1 carries the decision and the measurement that
   the fixed point is the truth. The follow-on this once carried — whether depth *binning* disturbs
   the multi-library score — is also closed: it does not, by 0.054 rungs under the ladder §4 now
   fixes (research note §4.3).
2. **Do the two `F` estimators agree on the tomato cohort, and where they disagree, why?** — OPEN,
   and it is now a measurement rather than a choice. *Leaning:* they will disagree on the
   het-inflated sample, and in the direction §6.2's first bullet predicts. **Settled by:** running
   both and comparing, which §6.3 makes cheap.
3. **Adopt the reference-bias term in place of `½`?** — OPEN, and **half of it is now closed**:
   the prior question, whether the misspecification costs anything worth a parameter, is answered
   yes for shallow cohorts and no for deep ones
   ([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
   §5). A decision rule agreed before the measurement — under one rung in `ε` and under 5% relative
   in both frequencies, over `b` from 0.44 to 0.50 — would have closed this as "no". It fails, and
   it fails in one place: **heterozygosity at 3 reads a site**, which comes back 6.2% low at
   `b` = 0.46 and 9.8% low at 0.44 (and 4.5% low at 0.47, just inside the bar). `π_hom_alt` stays
   under 4% everywhere, and everything passes at 6 reads and above. *Leaning:* yes; costs one parameter, and it earns its
   keep on tomato-depth cohorts rather than in general. **What remains to settle:** whether the
   fitted value departs from ½ by more than its standard error on real data — the same measurement
   as before, now with a stated size of what it would buy.
4. **Fit the homozygous-non-reference rate freely, or derive it from heterozygosity under
   Hardy-Weinberg?** — **CLOSED: fit it freely.** Hardy-Weinberg would save a parameter by tying the
   two together, and §6.1's identity shows what the tie assumes: `π_hom_alt = F·f + (1−F)·f²`, which
   equals the Hardy-Weinberg `f²` only at `F = 0`. Tying them therefore presumes the answer to the
   quantity §6 is estimating. §5's argument from the data — heterozygosity is within-individual, the
   homozygous-non-reference rate is a distance to the reference accession, and a selfing crop pulls
   them apart — points the same way, and is now the corroboration rather than the reason. **Worth
   checking anyway, as a signature rather than a decision:** fit both ways on the tomato cohort and
   confirm the free rate departs from the Hardy-Weinberg prediction most in the most inbred samples,
   which is what §6.1's identity predicts.

---

## 12. How we know it works

*These are this path's own. The tests that span both paths — recovering known parameters, the
histogram-versus-census comparison, the read-group grain, scan spacing, memory, determinism — are
[`parameter_prepass.md`](parameter_prepass.md) §10.*

1. **`F` is recovered from synthetic data at known inbreeding**, across at least two levels, since
   §6.3 shows the cohort's diversity divides by `1 − F` and is most sensitive where `F` is large.
   Simulate **runs** rather than a uniform excess of homozygosity — a genome with the right overall
   `F` but no segmental structure would let a broken runs model score well by accident. **Score
   against the genome's *realised* autozygous fraction, not the transition rates' nominal one**: a
   finite genome does not have the `F` its rates imply, and comparing against the nominal value
   reads sampling as bias. The harness is
   [`examples/ng_inbreeding_harness.rs`](../../../../examples/ng_inbreeding_harness.rs).
   1a. **The floor is measured at the run's own window count** and reported beside `F`: fit a genome
   generated with **no runs at all** and record what comes back (§6.1, §6.5). It is the estimator's
   resolution, it changes with genome size, and without it a small `F` cannot be read.
   1b. **A floor of false heterozygotes does not move `F`** — the property §6.2's second reason
   rests on, and untested until the research note measured it. Generate at up to five spurious
   heterozygotes per kilobase against a real rate of one, and assert `F` holds. **Run it from starts
   that disagree about the state separation**, since starts that do not are what turns this case
   into a silent `F` = 0 (§6.5).
2. **The two `F` estimators agree on synthetic data with no population structure**, where §11.2
   predicts they should, and the ratio is fed the diversity the runs estimator produced (§6.3's
   order). This is the check that the whole `Hobs = Hexp(1 − F)` relation is implemented consistently
   in both directions; disagreement here is a bug rather than biology.
3. **The soft heterozygote count matches a hard count where a hard count is safe.** On high-coverage
   data where genotypes are unambiguous, `expected_hets(window)` must agree with the number of
   heterozygous sites actually there. This checks the derivation of §4, not the estimator.
4. **The accumulator survives a change of model — within the sample's walk.** Refit with a different
   error rate and recompute every window's heterozygote count from the histogram alone, without
   re-reading a read. This is the sufficiency property §4 claims, and it is what makes the estimator
   swappable *while the histogram is still in hand*. **It is not a claim that the swap can happen
   later**: the histogram is dropped when its sample finishes
   ([`parameter_prepass.md`](parameter_prepass.md) §1.3), so revisiting §6.2 on a finished run means
   walking again, or having persisted the histogram on purpose.
5. **The sample's heterozygosity is right on a multi-library sample.** Generate two libraries from
   **one genome** with different error rates; the fit must return one heterozygosity. This is the
   test that would have caught the design §5 replaced, and it is specified in full in
   [`parameter_prepass.md`](parameter_prepass.md) §10.3 because it shares the read-group harness.
6. **The read-group histogram is the windowed one summed, on a single-library sample.** Fold the
   windowed histogram over its windows and compare it cell for cell against the read-group histogram,
   on a sample with one read group. They must be identical. This is the property §1 relies on to
   avoid a third object; it compares only objects the walk already builds, rather than accumulating a
   whole-sample histogram for the test alone; it needs no simulated truth; and it runs on the tomato
   CRAMs and on HG002 as they stand, since every sample in both carries one read group
   ([`parameter_prepass.md`](parameter_prepass.md) §5).
7. **The runs model uses the homozygous-non-reference signal, not only the heterozygote one.** In
   test 1's simulator, an autozygous stretch must be generated as §6.1 describes — one allele draw
   doubled, so the homozygous-non-reference rate *rises* inside a run — not merely as a stretch with
   its heterozygotes suppressed. A simulator that only suppressed heterozygosity would let a
   two-class emission pass, which is the model §6.1 replaced.
8. **The multi-library cell key is unbiased, and the check is arithmetic rather than a simulation**
   ([`examples/ng_multilib_key_harness.rs`](../../../../examples/ng_multilib_key_harness.rs)).
   Weight each cell by its exact probability under a known truth and maximise; what comes out is
   what an infinite genome would return, so any departure from truth is bias with no sampling noise
   to argue about. **Assert zero** on every one of the harness's worlds. Three algebraic checks run
   first, none needing a fit: the scoring rule sums to one over the cell space at any parameter
   values; no cell is ever charged a negative count of reference reads; and with every library's
   error rate set equal the rule reproduces the exact per-library likelihood to floating point,
   since there is then nothing to attribute. **A rule that fails any of the three cannot be
   unbiased**, and each of them is one line — this is where the average-share plug-in §1 retired
   would have been caught the day it was written. A sample-size ladder over 10⁴ to 10⁷ sites then
   confirms the implementation converges on the arithmetic's answer, and gives the precision the
   key costs against scoring every read against its own library.
9. **The depth ladder is checked in the band where a ladder can hurt, and against a control.** The
   sweep behind §4's ladder lives in the same harness (`--only=binning`), and two of its properties
   are what a replacement must reproduce rather than merely beat. **The control:** the exact
   ladder — one bin per depth — must return the unbinned answer to floating point, since a binned
   fit that moved there would be reporting the harness's arithmetic and not the binning's. **The
   band:** 10 to 30 reads a site. At tomato's 3 reads, 97 sites in 100 are never binned at all, and
   at 60 the genotype is certain whatever the depth, so a ladder checked only at either extreme
   passes whatever it does. The measured bias of the adopted ladder is 0.054 rungs and 0.3%
   ([`../research/parameter_estimator_experiments_2026-08-06.md`](../research/parameter_estimator_experiments_2026-08-06.md)
   §4.3).
10. **No cell is scored at a depth below its own alternative count.** Assert it directly on the
    accumulator, at every cell of every fit. It holds by construction when the mean is taken per
    cell and fails on 0.3% of sites when it is taken per bin, and what that costs is a fit 5.2 rungs
    low and 29% low on `π_hom_alt` **with nothing on the outside to show for it** — not the railed
    scan an earlier draft of the architecture doc predicted, which `argmax_at_ladder_end` would have
    caught (research note §4.5). The assertion is the only thing standing there.
