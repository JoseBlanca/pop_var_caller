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

**This object is a compromise, and what is scarce is worth being precise about.** The comparison
between samples is genuinely out of reach in a per-sample pass, which never holds two samples at once.
The *data* is not: the walk visits every position and could record all of it, at the encoding of §5,
for about 400 MB per sample on tomato and 20 GB across fifty. **What is traded is storage, not
availability.**
So the size of this object is a budget with a knob on it (§5), and a small budget suffices because
everything reading it is a distribution: the precision of a distribution estimated from a random
sample of sites depends on **how many sites it holds, not what fraction of the genome they are**.

**What it cannot do is anything that needs a site's neighbours.** Two million positions across 800 Mb
is one every four hundred bases **on average, and the average is all the rule guarantees**: a uniform
hash leaves the gaps geometrically distributed, with the mode at zero. About **442,000 of the kept
positions have their neighbour within 100 bases, and 49,000 within ten**. Nothing spreads them out,
and in a selfing panel like tomato's — where linkage reaches far past a kilobase — no budget this
object could afford would make the kept sites independent of one another.

**What that does and does not cost.** It does not touch the rates: a mean over correlated sites has
the same expectation as a mean over independent ones, and unbiasedness needs only that the choice
never looks at the data (§3). It does make **every precision figure computed from the site count**
optimistic, by a factor set by how far linkage reaches in the panel (§5). And it leaves linkage,
haplotypes and anything else reading a stretch of genome out of reach — **for which random thinning is
the worst possible shape**, since it discards exactly the close pairs that carry linkage information
and keeps unlimited distant ones, which carry none. The instrument for those is a small clustered
budget rather than a larger scattered one (§8).

### 1.1 Two sets, because STR loci are a different population

**Decision: build two censuses — one over ordinary sites, one over STR loci.** The generic set must exclude STR repeat tracts, because their
variability would distort every substitution statistic computed from it (§3). That same fact, read
the other way round, is why the STR set is worth having: **a repeat tract mutates orders of magnitude
faster than a base does, so the population's diversity at STR loci is a different quantity, not a
correction to the generic one.**

**What the two sets share, and what they do not:**

| | generic set | STR set |
|---|---|---|
| unit | a position | a locus, as region typing delimits it |
| domain | analysed regions, **minus** STR loci | analysed regions, **restricted to** STR loci |
| selection rule | §3, identical | §3, identical |
| what is stored | reads supporting A/C/G/T + other (§2) | reads at each whole-repeat offset, the non-whole-repeat guard, and the mismatching bases as a difference list (§2.1) |
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

**Keying at the finer grain serves both.** Summing a
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

**The alleles are lengths, so the allele record is a length distribution.** Each sample records, per
selected locus, **how many reads showed each whole-repeat offset from the reference tract length**,
plus one bucket for reads that differ by something that is not a whole number of motif copies.

**The bases are not thrown away, and keeping only a summary of them is a compromise to save space.**
Ideally a locus would keep the reads themselves: then every departure from the reference tract could
be *attributed* rather than counted — a substitution inside the tract, which interrupts the motif and
changes what a repeat unit is, against ordinary sequencing error in the flank; which base replaced
which; and whether two reads carried the same interruption, which is what makes an interruption an
allele rather than an error. Storing them is what costs: at about 100 bases compared per read and
three reads a locus, two-bit packed, a locus is 75 bytes per sample, so a million loci across fifty
samples is **3.7 GB** against the ~50 MB the whole record set is budgeted at (§5).

**So keep the differences instead, sparsely — mismatches are rare, which is what makes this cheap.**
At 300 base comparisons a locus and an error rate of 0.002 a locus carries about **0.6 mismatching
bases**, so a list of *(which read, offset from the tract start, which base)* costs about two bytes
per locus — roughly what a bare pair of counters costs, and the same trick §5 uses for the generic
set's non-reference observations.

- **The reference tract length is a property of the locus**, held once for the cohort rather than
  once per sample — the analogue of the reference base in §2, and available from region typing's
  delimitation.
- **The offset range is bounded and small.** §3 of
  [`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) measures where slippage mass lives: reads that move usually move one repeat, and the second step is
  already rare. A range of about −4 … +4 with saturating end buckets loses nothing that is fitted,
  and the "not a whole repeat" bucket is the same guard §2 uses — a locus where it is large is one
  the model cannot describe.
- **Bases compared, and the mismatches as a difference list.** The count of bases compared is the
  denominator the STR path's substitution error is fitted against; the differences are its numerator,
  and unlike a bare count of them they say *where*. Offsets record *length*, so a substitution that
  does not change a tract's length is invisible to them and `ε` cannot be recovered from offsets alone
  ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.1) — without this channel
  [`parameter_prepass.md`](parameter_prepass.md) §4.1's comparison would cover three of that path's
  four parameters instead of four. **The flank-against-tract split is what a pair of counters could
  not give**, and it is what any interrupted-repeat model would have to read.
- **The guard bucket should say what it caught**, by the same sparse mechanism: record a non-whole-repeat
  read with its offset and its size, so a partial unit at a tract edge — which is alignment ambiguity —
  is distinguishable from an indel in the flank. A bare count can only raise a threshold, never
  explain it.
- **Four states, not two, and the pair that needs saying is not the obvious one.** A sample at an STR
  locus may have had **no read reach the locus at all**; **reads that reached it but none crossing the
  whole tract**, so none reports a length; or **reads that crossed it**, whether they showed the
  reference length or another. (A fourth — *never walked* — is a bug rather than data.) The generic
  set's zero-depth-against-monomorphic distinction is the *last* pair, not the first.
- **And a read that did not cross the tract is not nothing.** It proves the tract is **at least** as
  long as the stretch it covered — a censored observation, which
  [`locus_generation_ssr.md`](locus_generation_ssr.md) records deliberately and whose admission gate is
  overlap rather than spanning for exactly that reason: 7,085 such reads on chromosome 1 of tomato
  SRR7279503 alone. **This record has nowhere to put one, which is a gap rather than a decision.** What
  makes the first pair worth keeping apart is that a tract longer than a read is never crossed, in
  every sample at every depth — censoring that runs along repeat count, the axis the slippage numbers
  are fitted within, so a stratum unobservable with this read length must not look like one that was
  unlucky with coverage.

**Why this cannot reuse §2's shape**: five per-base buckets cannot express a length, and a length
distribution cannot express which base was substituted — which is why the length distribution travels
with a difference list rather than alone. The two sets share the selection rule and the binning rule,
and nothing about their contents.

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
- **The run's identity starts with three things: the seed, the reference and the region set.** All
  three travel with the sample's summary, as a seed value and two digests. Two samples analysed under
  different region sets share no loci, and a cross-sample estimate over mismatched sets is not noisy
  but meaningless, so the gather must be able to refuse rather than average
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2). *A BED is not an incidental
  detail of how a run was invoked; it is part of what defines the data produced.* **The full list is
  longer than three** — the STR side adds the catalog's build settings, the routing criteria and the
  caps ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §5.1) — **and none of it
  is a substitute for the digest of the positions actually selected** (§5.2), which is the only value
  that witnesses the answer rather than the question.
- **Never select on the data.** "Positions that looked variable" is ascertainment: variability is a
  function of depth and error, so selecting on it conditions on the quantity being measured — the
  bias [`parameter_prepass.md`](parameter_prepass.md) §2 exists to remove, one level up. The region
  set is *not* an exception: it is chosen before any read is examined, so restricting to it narrows
  the population being described without biasing the estimate within it.

**Scattered positions, not sampled regions.** Contiguous blocks would be cheaper to fetch and carry
less per position: sites within a block share a genealogy, so *k* positions in one block say less
about a genome-wide rate than *k* scattered ones. **What scattering does not buy is independence** —
at one position in four hundred most kept sites still have a neighbour within a kilobase (§1), and in
a selfing panel linkage reaches further than that. The estimators downstream treat sites as
independent because it makes their arithmetic simple; what that costs is precision figures that are
optimistic (§5), not estimates that are wrong. Scattering buys the most information per byte stored,
and nothing more.

**Step 3 has already marked which stretches of the reference are repeat tracts, and both sets read
that marking rather than deciding again.** The generic set draws its positions from the analysed
regions *outside* those tracts; the STR set draws its loci from the tracts themselves, so every
analysed position belongs to exactly one of the two.
A repeat tract mutates orders of magnitude faster than a base does, so mixing them would have a
diversity estimate averaging a substitution rate together with a slippage rate, which describes
neither. The "other" bucket of §2 is the guard for what slips through the boundary: a site where it
is not small is not a clean substitution site, and the gather can drop it.

Region typing delimits STR loci from the **reference alone**, never from the reads, so its output is
identical for every sample by the same argument the hash rule relies on. Both domains are therefore
still pure functions of the run's inputs, and both sets inherit the property of §3's opening. The
region-typing output does, however, depend on the copy floors.

**Applying the rule to loci rather than positions.** For the STR set, hash the locus's start
coordinate; the target count is a number of loci and the denominator is how many STR loci the
analysed regions contain, not their total length. Everything else — the seed, the ban on selecting
on the data, the independence of samples — carries over unchanged.

> **Superseded for the STR set — equal per stratum, not proportional**
> ([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3). A uniform hash makes the
> retained set a scale model of the genome's, which gives the fewest loci to the strata where
> slippage is largest and least known — it varies twenty-two-fold across repeat counts. The
> replacement is a per-stratum cap filled by an even spread over the analysed regions, which keeps
> every property this section requires because a stratum is read off the reference tract. It has one
> price: the STR set stops being a random sample of STR loci, so anything pooled **across** strata —
> STR diversity — must be reweighted by each stratum's true locus count.

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

**The knob is a number of positions, and about two million of them is the default.** It has to be a
position count and not anything derived from the data: positions are what memory scales on, they are
what §3's rule can deliver without seeing a read, and they are the only figure a user can set before
the run.

### 5.1 Which count each parameter's precision actually rests on

**Not all of them rest on the same count, which is why an earlier version of this section stated the
target as "ten thousand sites variable across the cohort" and why that was wrong as a knob.** Three
groups:

- **The per-base error rate rests on read observations, and every position carries them.** Two million
  positions at three reads is six million observations whether the panel is diverse or clonal, which
  is what pins the rate to about one part in eighty (§6). Variable sites are irrelevant to it — they
  are a rounding error in the denominator.
- **Heterozygosity, the homozygous-non-reference rate and diversity rest on the position count too,
  for the error the caller consumes.** These are proportions: the standard error is `√(p(1−p)/n)`,
  which falls as `n` grows *and* as `p` falls, so a low-diversity panel is estimated more tightly in
  absolute terms, not less. It is only the **relative** error that goes as one over the square root of
  the variable-site count — and nothing downstream consumes a relative error, which
  [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.2 settles at length in
  rejecting a budget of the form `positions = target / min(Hobs)`. **A target stated in variable sites
  is that same rejected rule wearing a different name**, since it makes the budget rise as the panel
  gets less diverse.
- **The frequency spectrum, contamination and relatedness are the ones that genuinely need variable
  sites**, because a monomorphic position carries no information about any of them: the spectrum's
  classes are populated only by segregating sites, contamination is read from stray reads carrying
  alleles *the panel segregates*, and two samples' relatedness is invisible where everyone is
  identical. Ten thousand segregating sites is ample for these — under a neutral shape that puts a
  couple of thousand in the singleton class and tens in the high-frequency tail.

**So the two counts are linked by the panel's segregating rate, which is a property of the cohort and
not a knob.** At a rate near 1 in 200 bp, two million positions yield about ten thousand segregating
ones, which is where the default comes from. **The variable-site count is therefore an outcome to
report, not a target to hit** — the run should print it, and a run that returns a few hundred is one
whose spectrum, contamination and relatedness figures are thin while its error rate and its rates are
exactly as good as anywhere else. That is the honest statement, and it is per parameter.

**Two million sounds expensive and is not. One step in the reason is easy to lose, and losing it
makes the object ten times bigger.**

**The positions are never stored.** They are reproducible from §3's rule, so every sample and the
gather derive the identical list. What a sample stores is a **dense array in position order** — entry
*i* is the *i*-th selected position — with no coordinates, no keys and no index. **Coordinates would
cost about five bytes each, 10 MB — which is not much in itself, and is eight times the 1.25 MB of
data they would be indexing.** That ratio is the argument, not the absolute figure.

### 5.2 What proves that every sample really did select the same positions

**Agreeing on the inputs is an indirect check, and it is the weaker one.** The seed, the reference
digest, the region-set digest and the other identity values
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §5.1) say the two runs were
*asked* the same question. They cannot see a hash function that changed between versions, a threshold
computed in 64 bits on one machine and 128 on another, or a walk that filled its array from the wrong
end. Every one of those leaves the inputs matching and the sets different, and a set difference is not
a noisy estimate but a meaningless one.

**So carry a direct witness: a digest of the selected positions, computed by the walk itself.** Thirty-two
bytes per sample, against 1.25 MB of records. Two requirements make it worth its size:

- **It must be produced where the entries are, not re-derived.** A digest computed by running §3's rule
  a second time proves that the rule is deterministic, which nobody doubts; it says nothing about the
  array that was actually written. Hash each position **as its entry is filled**, so the digest
  witnesses the array's own order and length.
- **Block it, so a mismatch can be localised.** One digest per megabase alongside the whole-set one is
  800 entries of 8 bytes on tomato — **6.4 kB, half a percent of the record** — and it turns "these two
  samples disagree" into "they disagree in this megabase", without storing a single coordinate. That is
  the middle ground between 32 bytes that only says *no* and 10 MB that says everything.

**This is a different failure from the never-walked sentinel**
([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §2.2), which catches a
sample missing entries it should have. The digest catches a sample holding entries for *other*
positions — the same count, the same shape, and nothing anywhere else to notice.

**Each entry is one binned depth, which fits in four bits.** Sixteen bins is exactly what §4's scheme
produces: integers 0–8, then seven geometric bins above.

> **Superseded — five bits, twenty bins**
> ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §2.2). Four bits was chosen
> here for the encoding; [`parameter_prepass_generic.md`](parameter_prepass_generic.md) §4
> subsequently measured the two ladders and found sixteen bins cost ten times the bias of twenty —
> 0.55 rungs of the error-rate ladder against 0.054. Five bits also leaves room for the
> never-walked state §7 requires, which four bits does not. The sizes below rise by a quarter, to
> about 1.25 MB per read group. At three reads a site a sample's record is
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

**The position count is a configuration knob, and should be presented as a memory budget.** It is the
one dial trading storage against the precision of everything in
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md); the relationship is transparent — bytes
scale linearly with positions, precision with their square root — and **the two are the only things
worth reporting when the budget is swept** (§9): what it cost at rest, and what each parameter's
precision was, **parameter by parameter, because §5.1 says they do not degrade together.** A budget
cut ten-fold costs the error rate a factor of three and can cost the spectrum its rarest classes
outright.

**The square root is over *independent* sites, and the kept ones are not independent (§1).** For a
per-read quantity that changes nothing: errors are independent across reads whatever the sites do, so
§6's "one part in eighty" stands. For a per-site rate — heterozygosity, the spectrum, diversity — the
effective count is the number of independent stretches of genome the panel carries, not the number of
positions. If a 100 kb stretch behaves as one draw, two million positions carry 8,000 draws and every
such interval widens sixteen-fold; the true factor lies between about 3 and 16 and depends on how far
linkage reaches in the panel, **which is measurable from cohort calls that already exist and is not
measured here**.

One thing moves the default: a study willing to spend more storage gets proportionally better
cross-sample estimates. **A sparsely segregating panel is not a second reason** — that would be the
budget rising as diversity falls, which §5.1 and
[`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.2 both reject. What such a
panel gets is a thinner spectrum, reported as such.

### 5.3 The STR set is a different regime, and may not need sampling at all

**The first thing to measure is how many STR loci there are**, because the arithmetic above may not
apply. Two differences pull in opposite directions:

- **There are far fewer STR loci than genome positions**, so the pool to sample from is smaller by
  orders of magnitude.
- **A far larger fraction of them vary in a cohort.** Repeat tracts are hypervariable, so where the
  generic set needs about two hundred positions to yield one that segregates, an STR set may need
  only two or three.

Together those mean the ten thousand *segregating* loci the spectrum wants (§5.1) could well come out
of keeping **every** STR locus the analysed regions contain — in which case the sampling step is a
no-op and the threshold is set to keep everything. *Here the segregating count is a check on whether a
cap is needed at all, not a target the cap is solved for; the cap itself is still a locus count.* The cost stays small either way: at eight offset
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
- **Anything requiring linked sites** — linkage disequilibrium, how it decays with distance,
  haplotypes. **Random thinning is the wrong instrument, and not because it leaves too few sites**: it
  discards exactly the close pairs where linkage information lives and keeps unlimited distant ones,
  which carry none (§1). The instrument is a **clustered budget** — a few hundred short blocks with
  the positions inside them kept densely, chosen by the same hash rule applied to block starts. At 500
  blocks of 5 kb keeping one position in fifty that is 50,000 positions, about 2% of the generic
  budget, and it buys millions of pairs spanning 50 bp to 5 kb. **Not built, and nothing in step 4
  asks for a distance-dependent parameter today** — but it is also what would measure the linkage
  extent §5's precision figures depend on. **Home:** here, if a consumer appears; otherwise the cohort
  caller, which sees every site and needs none of it.

---

## 9. How we know it works

1. **Every sample selects the same positions.** Run §3's selection over samples with different
   coverage, different read lengths, and a region-sharded walk at several thread counts; the selected
   set must be **identical every time**, and identical to the set computed directly from the seed, the
   reference and the region set. **Include a `--regions` case**: the selected positions must all lie
   inside the BED, and the realised count must hit the target computed from the region set's length
   rather than the genome's — the arithmetic most likely to be written against the contig table by
   reflex.
2. **The digest witnesses the array, not the rule** (§5.2). Corrupt one entry's position in a walk
   that is otherwise correct — same count, same shape — and the whole-set digest must change and the
   megabase digest must name the block. **Then check the test can fail**: a digest re-derived from
   §3's rule instead of from the filled entries passes this unchanged, which is the whole reason the
   requirement is stated. And two samples whose digests differ must be **refused rather than pooled**.
3. **The selection is unbiased.** On synthetic data with a known frequency spectrum, the spectrum
   estimated from the selected sites must match the one computed from every site, within its own
   error. If the selection ever came to depend on the data, this is where it shows.
4. **Empty positions round-trip.** A position with zero depth, a position with reads but none
   non-reference, and a position never walked must be three distinguishable states after a write and
   read.
5. **Memory and precision are measured together, over a sweep of the position budget** — §5's figures
   are arithmetic, and these two are what the knob trades against each other (§5). Refit at several
   budgets and report, at each: **bytes at rest**, and **each parameter's error against a known truth,
   one row per parameter**. Pooling them into a single "precision" hides the finding §5.1 predicts —
   that the error rate and the rates degrade as the square root of the budget while the spectrum's
   rare classes empty out entirely, so there is no single budget at which everything is still good.
   **Report the segregating-site count each budget yielded** beside them, as the outcome it is. Memory
   alone is also covered by [`parameter_prepass.md`](parameter_prepass.md) §10.6, which reports each
   object separately.
6. **The two sets partition, and neither leaks.** No locus appears in both; every analysed position
   is in exactly one domain or excluded for a stated reason. Run it with region typing's copy floors
   moved and confirm both sets change together — the partition follows the delimitation rather than
   a second copy of its rules.
7. **The STR set's ends saturate, and its guard catches.** A read at an offset beyond the stored range
   must land in the end bucket rather than being dropped or wrapping, and a read whose tract differs
   by a non-whole number of motif copies must land in the guard bucket with its offset and size. Both
   are cheap to assert and both are silent when wrong.
8. **The difference list places a mismatch where it happened** (§2.1). A substitution planted in the
   flank and one planted inside the tract must come back distinguishable, and two reads carrying the
   same interior substitution must come back as two entries at one offset rather than as one entry or
   as a count of two. **That last one is what an interrupted-repeat model would read**, and a
   read-blind encoding passes every other check here.
9. **The four states at an STR locus round-trip** (§2.1): no read reaching the locus, reads reaching
   it but none crossing the tract, reads crossing it and showing the reference length, and the region
   never walked. The second is the one to plant deliberately — it is the state this record has no
   field for, so the test either shows a lower bound surviving or shows that it does not.
10. **The two numbers that set the STR threshold are measured** (§5.3): how many STR loci the analysed
   regions contain, and what fraction of them vary across the cohort. Until they exist, the threshold
   is a guess, and §5.3's expectation that the set may need no sampling at all is untested.
