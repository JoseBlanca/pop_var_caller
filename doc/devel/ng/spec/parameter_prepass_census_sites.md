# ng — the parameter pre-pass: the census sites

*Design spec, 2026-08-03. **No code yet — this settles the design.** One of five documents covering
ng step 4; the shared framing is in
[`parameter_prepass.md`](parameter_prepass.md), which this assumes. **Scope: the two censuses — small
sets of loci at which every sample keeps its raw evidence instead of folding it into a histogram**,
one drawn from ordinary sites and one from STR tracts. Why they exist, how the loci are chosen, what
is stored, and what it costs. Their consumers are all in
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md). `src/ssr/` and `src/pileup/` are frozen
production: everything said about them here is a record, not a change.*

*The name is meant literally: a census asks **the same questions of everyone**. These are the loci at
which every individual is asked, so that answers can be compared — which is exactly what a histogram,
having forgotten which locus it saw, can never support.*

---

## 1. Why a second kind of object exists at all

**Every histogram in [`parameter_prepass.md`](parameter_prepass.md) is a per-sample marginal, and a
marginal cannot say which samples carried an allele together.** That correspondence is the whole
content of a frequency spectrum, of a contamination estimate, and of a relatedness matrix. It is not
recoverable from summaries at any grain, however fine the key — a histogram has already forgotten
which site each observation came from. So those parameters get an object of their own.

**What separates the two objects is not where the data comes from.** Every parameter in step 4 would
ideally use the whole genome. The split is whether the parameter can be computed **one sample at a
time**: an error rate can, because it is a property of one sample's reads; a frequency spectrum
cannot, because it is a comparison between samples at the same site.

**This object is a compromise, and it is worth being precise about what is scarce.** The *comparison*
is genuinely out of reach in a per-sample pass — it never holds two samples at once. The *data* is
not: the walk visits every position and could record all of it, at the encoding of §5, for about
400 MB per sample on tomato and 20 GB across fifty. **What is traded is storage, not availability.**
So the size of this object is a budget with a knob on it (§5), and a small budget suffices because
everything reading it is a distribution: the precision of a distribution estimated from a random
sample of sites depends on **how many sites it holds, not what fraction of the genome they are**.

**What it cannot do is anything that needs a site's neighbours.** Two million positions across 800 Mb
is one every four hundred bases, and they are spread out on purpose (§3) so that no two are close
enough to be inherited together. Linkage between variants, haplotypes, and anything else that reads a
stretch of genome rather than a pile of separate sites are out of reach by design.

### 1.1 Two sets, because STR loci are a different population

**Decision: build two censuses — one over ordinary sites, one over STR loci — sharing the
selection machinery and nothing else.** The generic set must exclude repeat tracts, because their
variability would distort every substitution statistic computed from it (§3). That same fact, read
the other way round, is why the STR set is worth having: **a repeat tract mutates orders of magnitude
faster than a base does, so the population's diversity at STR loci is a different quantity, not a
correction to the generic one.**

The STR path already depends on such a number and takes it on faith: `SFS_THETA = 0.01`
([`src/ssr/cohort/freebayes_emit.rs:42`](../../../../src/ssr/cohort/freebayes_emit.rs)), described in
its own comment as *"freebayes' default `-T`"* — a **SNP-scale** constant governing repeat tracts.
It is not incidental. The same comment explains that each distinct allele pays a factor of `θ`, so
this is precisely the number deciding how much read evidence a rare STR allele must show before it is
believed, and too small a value suppresses real STR variation. That is worth holding beside the
recorded gap between our STR emissions and HipSTR's — not as an explanation, since that gap has other
documented causes, but as a hypothesis nobody could test while the number was never measured.

**A second, less certain payoff.** Cross-sample evidence at STR loci is also what a per-locus stutter
estimator would need, which makes the comparison in
[`parameter_prepass.md`](parameter_prepass.md) §4.2 possible at all. Whether it improves the
slippage priors is genuinely unknown — the per-stratum histogram already pools across loci, and it is
not obvious that cross-sample data at one locus improves a *chemistry* parameter. **The diversity
payoff alone justifies the set; treat the slippage one as a question the set makes askable.**

**What the two sets share, and what they do not:**

| | generic set | STR set |
|---|---|---|
| unit | a position | a locus, as region typing delimits it |
| domain | analysed regions, **minus** STR loci | analysed regions, **restricted to** STR loci |
| selection rule | §3, identical | §3, identical |
| what is stored | reads supporting A/C/G/T + other (§2) | reads at each whole-repeat offset + non-whole-repeat (§2.1) |
| reference value held once per locus | the reference base | the reference tract length |
| consumers | diversity, spectrum, contamination, relatedness | STR diversity; possibly per-locus stutter |

---

## 2. What is stored at each site

**Decision: per-allele read counts, not a count of non-reference reads.** The histograms elsewhere
collapse everything that is not the reference into one number, which is right for a *rate* — an error
rate does not care which wrong base a read shows. It is wrong here, for three reasons that bite at
once:

- **Two samples' "non-reference" may be different alleles.** Where one sample carries T and another G
  against a reference C, an alternative count records both as non-reference, and a spectrum built
  from it would credit them with an allele they do not share.
- **The reference is not always the major allele.** Where the reference accession happens to carry
  the rare allele, "non-reference" *is* the common one, and any statistic treating the reference as
  the baseline has the site inside out.
- **Contamination is identified by *which* allele.** The test is whether a sample's low-level
  non-reference reads carry the alleles the panel carries at that site — a question an alternative
  count cannot answer at all.

**So each site records the reads supporting A, C, G and T, plus one bucket for anything else** (an
indel, a spanning deletion). Five numbers. The reference base is recorded once per site rather than
per sample, and no per-sample record assumes it is the common allele — deciding a site's allele set,
and which allele is major, is the gather's job because only the gather sees every sample.

**Store every selected position, including the empty ones.** A position where a sample has no
alternative read, or no coverage at all, is not missing data — it is the denominator, and dropping it
would leave the shape of a spectrum with no scale attached. Zero depth is recorded as zero depth,
which is different from monomorphic and the consumer must be able to tell.

**Decision: key by read group, and let the consumer sum.** Two consumers want two different grains.
Diversity, relatedness and the frequency spectrum are properties of individuals, so they want **the
sample's** counts at a position. The per-base error rate is chemistry, so it wants **the read
group's** ([`parameter_prepass.md`](parameter_prepass.md) §1.1) — and without that axis the censuses
could not stand as an alternative to the genome-wide histograms at all, since one of the four
parameters those produce is per read group.

**Keying at the finer grain serves both, because read group sits below the site.** Summing a
position's read groups is addition of raw counts at one place — 20 + 10 really is 30 — so the
sample's record is an exact fold of the read groups', not an approximation. The reverse is not
recoverable. *Rejected: key by sample and accept that the error rate must come from the histograms* —
it would make the comparison in [`parameter_prepass.md`](parameter_prepass.md) §4.1 unrunnable for
the one parameter it most wants to test.

**What it costs is a multiplier equal to the read groups per sample**, which is **1 for 1,550 of the
1,707 samples** in the archive survey ([`read_groups.md`](read_groups.md) §1), 2 or 3 for nearly all
the rest, and 42 for one outlier. §5's sizes are therefore per read group. **This holds for both
sets.**

### 2.1 What is stored at an STR locus

**The observation is a tract length, so the record is a length distribution.** Per-base allele counts
mean nothing here: the alleles are lengths, and two reads showing the same length agree whatever
bases they contain. Each sample records, per selected locus, **how many reads showed each whole-repeat
offset from the reference tract length**, plus one bucket for reads that differ by something that is
not a whole number of motif copies.

- **The reference tract length is a property of the locus**, held once for the cohort rather than
  once per sample — the analogue of the reference base in §2, and available from region typing's
  delimitation.
- **The offset range is bounded and small.** §3 of
  [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) measures where slippage mass lives: reads that move usually move one repeat, and the second step is
  already rare. A range of about −4 … +4 with saturating end buckets loses nothing that is fitted,
  and the "not a whole repeat" bucket is the same guard §2 uses — a locus where it is large is one
  the model cannot describe.
- **Bases compared and bases mismatched, as two counts.** The same composition channel the STR
  histogram carries ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1) and for the same
  reason: offsets record *length*, and a substitution that does not change a tract's length is
  invisible to them, so `ε` cannot be recovered from offsets alone. Without this the census route
  could not fit the STR path's error rate, and [`parameter_prepass.md`](parameter_prepass.md) §4.1's
  comparison would cover three of that path's four parameters instead of four.
- **Zero-depth and no-spanning-read are distinct, and both must be recorded.** A locus a sample
  simply did not span carries no information about its alleles, and is different from one it spanned
  and found unremarkable. The generic set makes the same distinction between zero depth and
  monomorphic.

**Why this cannot reuse §2's shape**: five per-base buckets cannot express a length, and a length
distribution cannot express which base was substituted. The two sets share the selection rule and the
binning rule, and nothing about their contents.

---

## 3. Choosing the positions

**The property to deliver:** *for a given genome and set of analysed regions, every sample selects
the identical positions, arriving at them independently.* Samples are walked separately — on
different machines, at different times, with no sample able to see what another chose — so the set
cannot be negotiated or handed round. It has to be computed from the run's inputs alone.

The rule is a pure function of the position and those inputs, never of the data:

> keep position `p` if `p ∈ analysed_regions` and `hash(contig, p, seed) < threshold`

- **The domain is the analysed regions, not the genome.** With no `--regions` BED that is the whole
  reference; with one, it is the reference **intersected with the BED**. Selecting genome-wide under
  a BED would leave nearly every chosen position unvisited, so a run over 1% of the genome would
  yield 1% of the intended sites and nobody would notice until the spectrum came out thin.
- **The count is set by the threshold, from the analysed length.** Choose
  `threshold / hash_range = n / analysed_length` for a target of `n` positions — `analysed_length`
  being the total length of the region set, not of the contig table. The realised count is binomial
  around `n`, so ±√n: about ±100 at `n` in the ten-thousands.
- **The run's identity is therefore three things: the seed, the reference and the region set.** All
  three travel with the sample's summary, as a seed value and two digests. Two samples analysed under
  different region sets share no loci, and a cross-sample estimate over mismatched sets is not noisy
  but meaningless, so the gather must be able to refuse rather than average
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2). *A BED is not an incidental
  detail of how a run was invoked; it is part of what defines the data produced.*
- **Never select on the data.** "Positions that looked variable" is ascertainment: variability is a
  function of depth and error, so selecting on it conditions on the quantity being measured — the
  bias [`parameter_prepass.md`](parameter_prepass.md) §2 exists to remove, one level up. The region
  set is *not* an exception: it is chosen before any read is examined, so restricting to it narrows
  the population being described without biasing the estimate within it.

**Scattered positions, not sampled regions.** Contiguous blocks would be cheaper to fetch and are the
wrong shape: sites within a block share a genealogy, so a block of *k* linked positions carries far
less independent information than *k* scattered ones. The estimators downstream assume sites are
independent, and scattering is what makes that nearly true.

**The two sets partition on STR loci, using the delimitation region typing has already produced
(step 3).** The generic set takes the analysed regions **minus** those loci; the STR set takes them.
A repeat tract mutates orders of magnitude faster than a base does, so mixing the two would have a
diversity estimate averaging a substitution rate together with a slippage rate, which describes
neither. The "other" bucket of §2 is the guard for what slips through the boundary: a site where it
is not small is not a clean substitution site, and the gather can drop it.

**The partition costs the selection rule nothing**, and that is worth checking rather than assuming.
Region typing delimits STR loci from the **reference alone**, never from the reads, so its output is
identical for every sample by the same argument the hash rule relies on. Both domains are therefore
still pure functions of the run's inputs, and both sets inherit the property of §3's opening. The
region-typing output does, however, join the seed, the reference and the region set as something the
gather must check agreement on — a run with different copy floors delimits different loci, so the
two sets would differ (§7's error model, and
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2).

**Applying the rule to loci rather than positions.** For the STR set, hash the locus's start
coordinate; the target count is a number of loci and the denominator is how many STR loci the
analysed regions contain, not their total length. Everything else — the seed, the ban on selecting
on the data, the independence of samples — carries over unchanged.

---

## 4. Binning

**Bin — but not into equal widths. This rule governs every depth binning in step 4.** Depth 1 and
depth 5 are different kinds of observation; depth 100 and depth 105 say almost the same thing. So the
bins must be **fine at the bottom and wide at the top** — exact integers up to about eight, widening
geometrically above, so that each bin costs a similar small *fraction* of what it holds rather than a
similar absolute slice of the depth axis. Two reasons, and the second decides it:

- **The likelihood moves fastest where depth is low.** Equal-width bins spend their resolution where
  it buys least.
- **At 3 reads per plant, nearly all the data is at the bottom.** On a Poisson model at mean depth 3,
  about **97 sites in 100 sit at depth 6 or below**. Equal-width bins coarse enough to be economical
  near depth 100 would collapse essentially the whole cohort into the first bin — the binning would
  destroy the dataset rather than compress it.

**Allele counts are not binned the same way.** The difference between 0, 1 and 2 reads supporting an
allele is most of the signal at low coverage, so small counts stay exact and only the tail is binned.
That applies to each of the five buckets, which is what keeps a per-allele record cheap: at three
reads a site nearly every bucket is 0 and at most two are non-zero.

---

## 5. How many positions, and what it costs

**The target is roughly ten thousand sites *variable across the cohort*** — ample, since under a
neutral shape that puts a couple of thousand in the singleton class and tens in the high-frequency
tail. How many positions must be selected to yield that many depends on how densely the panel
segregates, **which is the one number to measure before fixing the threshold**. At a segregating rate
near 1 in 200 bp it is about two million positions.

**Two million sounds expensive and is not. One step in the reason is easy to lose, and losing it
makes the object ten times bigger.**

**The positions are never stored.** They are reproducible from §3's rule, so every sample and the
gather derive the identical list. What a sample stores is a **dense array in position order** — entry
*i* is the *i*-th selected position — with no coordinates, no keys and no index. Storing coordinates
instead would cost about five bytes each, 10 MB before any data is recorded.

**Each entry is one binned depth, which fits in four bits.** Sixteen bins is exactly what §4's scheme
produces: integers 0–8, then seven geometric bins above. At three reads a site a sample's record is
nearly always *"n reads, all matching the reference"*, and that is the whole of it. Per-allele detail
is stored only where there is any, as a sparse list beside the array:

| part | size per sample | driven by |
|---|---:|---|
| depth array | 2 M × 4 bits = **1 MB** | the target count, nothing else |
| non-reference observations — index, allele, count, about 4 bytes each | **30–250 kB** | *errors*, not variants: 2 M × 3 reads × an error rate of 0.001–0.01 is 6,000–60,000 positions with a spurious non-reference read, against roughly 1,000 carrying a real one |

**The depth array alone reconstructs every empty position exactly**: depth `n` with no sparse entry
means `n` reads on the reference base, and the reference base belongs to the site, recorded once for
the cohort rather than once per sample.

So **about 1 MB per sample and roughly 50 MB for a fifty-sample cohort**, rising towards 1.3 MB if the
fitted error rate is at the high end. The same encoding over every position of the genome would be
400 MB per sample and 20 GB across fifty — a four-hundredfold saving, which is just 800 Mb ÷ 2 M
positions and does not depend on the bytes-per-site figure at all.

**The target count is a configuration knob, and should be presented as a memory budget.** It is the
one dial trading storage against the precision of everything in
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md); the relationship is transparent — bytes
scale linearly with sites, precision with their square root — and the sensible default is whatever
buys ten thousand variable sites on the cohort at hand. Two things move it: a panel that segregates
more sparsely than assumed needs a larger target for the same yield, and a study willing to spend
more storage gets proportionally better cross-sample estimates. Neither is a reason to rebuild
anything.

### 5.1 The STR set is a different regime, and may not need sampling at all

**The first thing to measure is how many STR loci there are**, because the arithmetic above may not
apply. Two differences pull in opposite directions:

- **There are far fewer STR loci than genome positions**, so the pool to sample from is smaller by
  orders of magnitude.
- **A far larger fraction of them vary in a cohort.** Repeat tracts are hypervariable, so where the
  generic set needs about two hundred positions to yield one that segregates, an STR set may need
  only two or three.

Together those mean the target of ~10,000 loci varying across the cohort could well be reached by
keeping **every** STR locus the analysed regions contain — in which case the sampling step is a
no-op and the threshold is set to keep everything. The cost stays small either way: at eight offset
buckets plus a spare, a locus is a couple of bytes per sample, so even a million loci is a couple of
megabytes per sample, comparable to the generic set.

**Decision: implement the selection anyway, even if the first cohort keeps everything.** Setting a
threshold that admits all loci is free; discovering later that a bigger genome needs sampling and
having no mechanism is not. **Soft:** both the locus count and the fraction that vary are unmeasured
here, and they are the two numbers that set the threshold (§9).

---

## 6. Could this replace the genome-wide histograms?

Worth asking, because two million sites is statistically ample for a *rate*: at three reads a site
that is six million read observations, which pins an error rate near 0.001 to about one part in
eighty. The answer differs per object.

- **The per-read-group generic histogram: plausibly yes on precision, but it saves nothing.** It is
  kilobytes per read group. Dropping it would simplify the code and shorten the walk, not reduce
  memory, so the case has to be made on simplicity rather than cost.
- **The windowed heterozygosity histogram
  ([`parameter_prepass_generic.md`](parameter_prepass_generic.md)): no, and this is clear-cut.**
  It estimates a *local* rate, window by window, not a global one. A 1 Mb window holds about a million
  sites but only ~2,500 selected ones, so at one heterozygote per kilobase it would carry about
  **2.5** expected heterozygotes instead of a thousand — far too thin to separate "inside a run" from
  "outside". Widening the window destroys the resolution the window exists for. **A subsample is the
  wrong instrument for a local quantity.**
- **The STR histogram: now askable, since the STR census (§1.1) supplies cross-sample evidence
  at repeat tracts.** Whether it can replace the per-stratum histogram is [`parameter_prepass.md`](parameter_prepass.md) §4.2, not this comparison.

**This is settled by measurement, not by leaning, and the plan is already set**
([`parameter_prepass.md`](parameter_prepass.md) §4.1): build both, fit the per-sample parameters both
ways against synthetic truth (§10.2 there), measure what each costs (§10.6 there), and decide
afterwards. Two things to hold on to while doing it:

- **Know what "good enough" means before reading the numbers, though the decision can still weigh
  them.** The histogram wins on precision by construction, so "which is more precise" is not the
  question. [`parameter_prepass.md`](parameter_prepass.md) §4.1 gives three criteria — is it *biased*
  (the one that cannot be waived, since bias is what a broken selection rule produces and more sites
  will not cure it), is the gap between the two routes wide enough to feel, and do any calls actually change — together
  with the precision this object is expected to reach, so that a surprise is recognisable as one.
- **Even a clean win here saves little memory.** The object it would displace is kilobytes per read
  group. The case for dropping it is simplicity — one accumulator instead of two — and it should be
  argued on that, honestly, rather than dressed up as a saving.
- **Expect the question to stay open, and plan for that.** Two of the three criteria need the
  downstream caller, which is not built, so only the bias check runs today. The question closes now
  **only** if the two routes agree to within rounding; short of that, both implementations stay and
  the numbers are recorded for whoever can finish the comparison
  ([`parameter_prepass.md`](parameter_prepass.md) §4.1).

---

## 7. Cross-cutting concerns

**Memory.** About 1 MB per **read group**, so 50 MB across a fifty-sample cohort where each sample
carries one (§5, §2) — the smallest of step 4's
accumulators despite being the only one kept per site. Held per sample during the walk; held for
every sample at once during the gather, where 50 MB is unremarkable.

**Concurrency.** The selection rule is a pure function of position, so a region-sharded walk needs no
communication: each shard fills the entries for the positions in its own region, and merging is
concatenation in position order.

**Errors.** A sample that fails to record a selected position — because the region was never walked,
say — must be distinguishable from one that recorded zero depth there. The first is a bug; the second
is data.

---

## 8. Deferred, with a recommended home

- **The estimators that read either set** — diversity, the frequency spectrum, contamination,
  relatedness, STR diversity, and any per-locus stutter model. This document settles what is
  accumulated and none of the arithmetic on top. **Home:**
  [`parameter_prepass_cohort.md`](parameter_prepass_cohort.md), except the per-locus stutter
  comparison, which is [`parameter_prepass.md`](parameter_prepass.md) §4.2.
- **Anything requiring linked sites** — a spectrum stratified by genomic context, linkage
  disequilibrium, haplotypes. The positions are scattered precisely so sites are near-independent
  (§3), which makes this the wrong object for those questions. **Home:** the cohort caller, which
  sees every site.

---

## 9. How we know it works

1. **Every sample selects the same positions.** Run §3's selection over samples with different
   coverage, different read lengths, and a region-sharded walk at several thread counts; the selected
   set must be **identical every time**, and identical to the set computed directly from the seed, the
   reference and the region set. **Include a `--regions` case**: the selected positions must all lie
   inside the BED, and the realised count must hit the target computed from the region set's length
   rather than the genome's — the arithmetic most likely to be written against the contig table by
   reflex.
2. **The selection is unbiased.** On synthetic data with a known frequency spectrum, the spectrum
   estimated from the selected sites must match the one computed from every site, within its own
   error. If the selection ever came to depend on the data, this is where it shows.
3. **Empty positions round-trip.** A position with zero depth, a position with reads but none
   non-reference, and a position never walked must be three distinguishable states after a write and
   read.
4. **Memory is measured, not assumed** — §5's figures are arithmetic. Covered by
   [`parameter_prepass.md`](parameter_prepass.md) §10.6, which reports each object separately.
5. **The two sets partition, and neither leaks.** No locus appears in both; every analysed position
   is in exactly one domain or excluded for a stated reason. Run it with region typing's copy floors
   moved and confirm both sets change together — the partition follows the delimitation rather than
   a second copy of its rules.
6. **The STR set records lengths, not bases, and its ends saturate.** A read at an offset beyond the
   stored range must land in the end bucket rather than being dropped or wrapping, and a read whose
   tract differs by a non-whole number of motif copies must land in the "other" bucket. Both are
   cheap to assert and both are silent when wrong.
7. **The two numbers that set the STR threshold are measured** (§5.1): how many STR loci the analysed
   regions contain, and what fraction of them vary across the cohort. Until they exist, the threshold
   is a guess, and §5.1's expectation that the set may need no sampling at all is untested.
