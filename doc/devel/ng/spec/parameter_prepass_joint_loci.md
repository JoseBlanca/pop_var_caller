# ng — the joint parameters fit: which loci are kept

*Design spec, 2026-08-10, revised 2026-08-11. One of three documents covering the **joint parameters fit**, ng
step 4's second route to every parameter it emits; read
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) first — it says what the route is,
what it produces and why it exists. This one settles **which loci every sample keeps evidence at**,
and nothing else. What is recorded at each is
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md).*

*Types and interfaces: [`../arch/parameter_prepass_joint_loci.md`](../arch/parameter_prepass_joint_loci.md).*

***Both halves are built.*** *The repeat catalog ships the per-stratum sampler this document asked for
([`repeat_catalog.md`](repeat_catalog.md), `src/ng/repeat_catalog/`), so §3 records why the rule is
that one and what using it obliges a consumer to do, rather than proposing it. **The generic rule (§2)
landed 2026-08-12** in `src/ng/parameter_estimation/joint/loci.rs`, with the unambiguous-base mask §2
now requires and the measurements in §2 and §4.5 behind it
(`examples/ng_joint_loci_probe.rs`).*

***It changes one decision in [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)***
*§3 — how STR loci are chosen (§3 below) — and that document carries a note saying so.*

---

## 1. What this is, in one paragraph

**The joint parameters fit reads a small set of loci at which every sample kept its raw evidence, instead of the
per-sample summaries the other route builds.** *The **parameters fit** is the estimation of the parameters the
caller will run on — error rates, heterozygosity, inbreeding, the population's diversity, repeat
slippage, contamination — done once over the whole cohort **before any variant is called**; it is not
the variant caller. To **walk** a sample is to pass over its alignments once, visiting each locus in
turn.* For that to work at all, every sample has to keep
evidence at *the same* loci — and the samples are walked separately, on different machines, at
different times, with no sample able to see what any other chose. So the set cannot be negotiated or
handed round. **This document is the rule that lets each sample arrive at the identical set on its
own**, for the two kinds of locus the two paths care about: ordinary positions for the SNP/indel
path, and repeat tracts for the STR path.

**It is a module of its own because it is a pure function of the run's inputs** — the reference, the
analysed regions, the repeat catalog, a seed, and a couple of caps. It touches no read, so it can be
built and tested with no alignment file in sight: run it twice and compare the lists.

**Both sides are built.** The STR side is two calls on the repeat catalog
([`repeat_catalog.md`](repeat_catalog.md), `src/ng/repeat_catalog/`), which ships the per-stratum
sampler; the generic side is the hash rule §2 states in full, in
`src/ng/parameter_estimation/joint/loci.rs`. What this document holds is the policy — which caps,
which seed, what has to travel with a sample, and what the stratification obliges a consumer to do
afterwards (§3.5).

### 1.1 Goals

1. **The same loci in every sample**, arrived at independently.
2. **Represent what the estimates need**.
3. **Make the size a knob with a legible meaning**, since the number of loci is what this route's
   memory bill scales on (§4).

### 1.2 Non-goals, and what it does not do

- **It does not say what is stored at a chosen locus.** That is
  [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md).
- **It does not decide which loci exist**, or which are repeat tracts. Region typing (step 3) does,
  and this rule reads its output.
- **It never looks at the data.** Choosing "positions that looked variable" would condition the
  selection on the very quantity being estimated. The analysed-region set is not an exception to
  that: it is fixed before any read is examined, so restricting to it narrows the population being
  described without biasing the estimate within it.

---

## 2. The generic set: uniform, scattered, and reproducible

**Unchanged from [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §3**, which
this section summarises rather than restates so that the STR rule beside it can be read in one place.

> keep position `p` if `p` lies in the analysed regions, carries an unambiguous reference base, and
> `hash(contig, p, seed) < threshold`

- **The domain is the analysed regions, not the genome.** With no `--regions` BED that is the whole
  reference; with one, it is the reference intersected with the BED. Selecting genome-wide under a BED
  would leave nearly every chosen position unvisited.
- **And it is narrowed once more, by the reference's own bases — ADDED 2026-08-12, measured.** A
  position inside an assembly gap is kept by a rule that cannot see it and then covered by no read in
  any sample. It is worse than wasted budget: the parameters fit derives its per-sample rates as means of
  genotype posteriors over the kept loci
  ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.2), and a locus with no reads
  contributes its **prior** — the model's own prediction — rather than evidence. Every such position
  is a fixed weight on the answer.

  **Measured** (`examples/ng_joint_loci_probe.rs`), at a two-million-position budget:

  | reference | not `A`/`C`/`G`/`T` | kept positions landing there |
  |---|---:|---:|
  | tomato SL4.00 | 44,731 bases (0.01%) | **135** |
  | GRCh38 with the hs38d1 decoys | 165,046,090 bases (5.31%) | **106,423** (5.32%) |

  **So it is nothing on tomato and one kept position in nineteen on human**, which is why the rule
  states it rather than leaving it to whoever writes the BED. **It costs no extra read of the
  reference**: the mask comes off the same forward pass that computes the reference's digests, the
  seam [`repeat_catalog.md`](repeat_catalog.md) already uses. **The threshold's denominator is the
  masked length**, not the contig table's total — the same trap the BED sets, one level down, and the
  arch doc makes the two impossible to state separately.
- **The threshold sets the count**, from the **selectable** length — the analysed regions after the
  ambiguous bases are taken out, which is the domain the two points above define together. The hash is
  a 64-bit value, so it takes `2^64` values and a target of `n` positions out of `selectable_length`
  is `threshold = 2^64 · n / selectable_length`. **Nothing about that range is measured or discovered — it
  is the output width of the hash we chose**, and the one property it rests on is that the hash
  spreads uniformly across it. Use the same `xxh3` the STR sampler already uses
  (`src/ng/repeat_catalog/strata.rs`, `hash_locus`), so both halves of this document rest on one
  assumption rather than two. Two mechanics worth having in advance: compute the threshold in 128 bits,
  since `2^64 · n` overflows a `u64`; and where `n >= analysed_length` it saturates and every position
  is kept, which is the right answer rather than a case to guard. The realised count is binomial around
  `n`, so about `±√n` — roughly ±1,400 at two million.

  *The STR rule needs no such range* (§3.2): keeping the `cap` lowest hashes compares hashes against
  each other, so it is scale-free. The threshold rule needs the range only because it converts a hash
  into a keep-or-drop decision on its own, and it can do that because its denominator — the analysed
  length — is known without scanning anything. The STR denominator is not, which is the whole of §3.1.
- **Scattered positions, never contiguous blocks.** Sites within a block share a genealogy, so *k*
  positions in one block say less about a genome-wide rate than *k* scattered ones. **What scattering
  does not buy is independence**, and the census-sites document overstated this: at one position in
  three hundred and ninety-one the gaps are geometric, so — **measured on tomato, not predicted** —
  **453,445 of 2,002,470 kept positions have their neighbour within 100 bases** and 1,847,538 have one
  within a kilobase, while in a selfing panel linkage reaches past a kilobase
  ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §1). The rates are
  unaffected — a mean over correlated sites has the same expectation — but the standard errors of §4.2
  are optimistic, and anything distance-dependent needs a clustered budget rather than a bigger
  scattered one.
- **Repeat tracts are excluded.** Step 3 has already marked which stretches of the reference are
  repeat tracts and this rule reads that marking; their variability would distort every substitution
  statistic computed from the set.

**There is nothing to stratify on here**, which is worth one sentence because §3 does stratify. The
generic path's rates are not stratified at all — the error rate, the heterozygosity and the distance
from the reference are each fitted over all sites at once
([`parameter_prepass.md`](parameter_prepass.md) §3's glossary) — so there is no group whose
representation could come out wrong.

---

## 3. The STR set: equal per stratum, not proportional

**A stratum is one (motif period, reference repeat count) pair** — the group of repeat tracts that
share a fitted slippage behaviour ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4). The
STR path fits its four numbers **within** each stratum, one set per (read group × stratum).

**Change to [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §3, which applies
the generic hash rule to STR loci too.** A uniform rule makes the kept set a scale model of the
genome's, and that is the wrong shape when the parameter is fitted per stratum.

**Why.** Slippage varies twenty-two-fold across repeat counts — 9 reads in 10,000 below four repeats
against 2 in 100 at six or more ([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5) — and the
strata where it is largest are the rarest in the genome. A proportional sample gives them the fewest
loci. **The strata whose parameter matters most, and is least well known, would get the least
evidence.**

### 3.1 Three things are wanted at once, and one of them forces a pass over the reference

We want the kept loci to be **(a)** a random sample within each stratum, **(b)** exactly `cap` of
them, and **(c)** decided without anything precomputed. **All three together are impossible**, and
the reason is arithmetic rather than a shortcoming of any algorithm: *"exactly `k` of `n`, drawn at
random"* is not defined until `n` is known. A rule that sees one segment at a time and must decide
there and then does not know `n`.

So one of the three has to give, and which one is the decision:

| rule | random within the stratum? | exactly `cap`? | needs nothing precomputed? |
|---|---|---|---|
| keep the first `cap` | **no** — measured 32.6% low, below | yes | yes |
| a hash threshold per stratum, `keep if hash(locus) < t` | yes | **no** — binomial around `cap`, and `t` still needs the stratum's size to set | yes |
| an even spread over the stratum's segments | **no** — systematic, and it depends on the order segments arrive in | yes | **no** — needs each stratum's total |
| **the `cap` lowest hashes in the stratum** | **yes** | **yes** | **no** — needs the stratum enumerated first |

**"The first N" is not "N of them", and the difference is a quarter of the answer — measured.** Over
20 Mb of tomato chromosome 1, capping mononucleotides at 6 repeats to the **first** 2,953 of 28,125
segments returns a slippage rate of 0.0541% against 0.0803% on the full set — **32.6% low** — with
every uncapped stratum in the same run agreeing to the digit. Tracts near the start of a chromosome
stutter measurably less, so a prefix is a biased sample and the bias runs one way.

**The last two rows both need a pass over the reference, and that pass is the same object.** Once it
is being paid for, the fourth row is strictly better than the third — genuinely random rather than
evenly spaced, and independent of the order segments are seen in. So the decision is really one
decision, not two: **pay for the pass, and then take the best rule it buys.**

### 3.2 Decision: the `cap` lowest hashes per stratum — and it is already built

**Give every STR locus a value `hash(contig, start, seed)` and keep, within each stratum, the `cap`
loci whose value is smallest.** Everything is kept where a stratum holds fewer than `cap`.

**This is `sample_loci_per_stratum(criteria, region, cap, seed)` on the repeat catalog**
([`repeat_catalog.md`](repeat_catalog.md) §5.3, `src/ng/repeat_catalog/strata.rs`), which lands with
that spec rather than with this one. **So this document's consumer holds no selection logic of its
own**: it states a policy, a cap and a seed, and receives the loci. What follows is why the rule is
that one — the reasoning belongs here because it is this consumer's requirement that shaped it.

**The seed and the cap are this run's, not the file's**, which is why they travel with a sample's
records (§5.1) rather than living in the catalog's header.

Four properties, and the last two are what the even spread lacked:

- **It is a uniform random sample of the stratum.** The hash is uncorrelated with position, so the
  chromosome-start effect that made a prefix 32.6% low cannot reach it.
- **It keeps exactly `cap`.**
- **It does not depend on the order loci are seen in.** The kept set is a function of the *set* of
  loci and the seed, nothing else — so a region-sharded enumeration and a single-threaded one agree,
  and merging two shards is taking the lowest `cap` of the union. The even spread has no such
  property: it counts segments as they arrive, so a different traversal keeps different loci.
- **Its working state is `cap` values per stratum**, a bounded heap, so enumerating a genome costs
  nothing to hold — a few hundred strata at a cap in the thousands is a handful of megabytes, and the
  counts of §3.5 are a tally beside the heaps rather than a second traversal.

**Rejected: the even spread over each stratum's segments**
([`examples/ng_str_stutter_rate.rs`](../../../../examples/ng_str_stutter_rate.rs), where it is
written and in use). It was the right answer to *"the first N is biased"* and it fixed that; what it
does not fix is that a systematic sample is not a random one, and it buys its evenness with an
order dependence. **It costs the same pass over the reference**, so nothing is saved by keeping it.

**Rejected: a hash threshold per stratum, decided locally as each segment is seen.** This is the one
rule that needs no pass at all, and it is genuinely tempting. It fails on its own terms: the
threshold that yields about `cap` loci depends on how many loci the stratum holds, so the knowledge
we were trying to avoid is needed anyway — and what is bought with it is an *approximate* count, since
the realised number is binomial around the target. A stratum holding 40 loci would sometimes keep 25.

### 3.3 The pass the rule needs is the repeat catalog, and it exists

**What the rule needs is every STR locus in the analysed regions with its stratum** — the (motif
period, reference repeat count) pair. That is the **repeat catalog**
([`repeat_catalog.md`](repeat_catalog.md)): every tandem repeat in the reference, found once during
the pass that already streams the whole FASTA, and written beside it. It records **repeats rather than
loci** — period, span, score, motif, purity — so which of them is an STR is decided when the file is
read, and the strata at any copy floor come out without a second scan.

**This document is the reason it exists** and is one of its consumers; it is not built here.

**What this consumer asks of it — two calls, one pass** ([`repeat_catalog.md`](repeat_catalog.md) §5.3):

- `count_loci_per_stratum(criteria, region)` — the true locus count of every stratum, which §3.5 makes
  load-bearing;
- `sample_loci_per_stratum(criteria, region, cap, seed)` — the kept loci, §3.2's rule.

Both take `region`, so a run under a `--regions` BED counts and samples within the BED rather than
over the whole genome — which §4.1's property requires and which nothing here has to arrange.

**The bound this places on the run, and it is checked rather than hoped.** The catalog is built at
settings at least as permissive as any reader will ask for, and a reader more permissive on any bounded
axis is **refused, not served a short list** — one comparison, `built.serves(asked)`. The file's build
floors are `[5, 5, 4, 4, 4, 3]` copies for periods 1 to 6, over periods 1..=6, with 15 bp of sequence
beside each tract. **The STR path's calling floors are `[8, 6, 6, 6, 5, 4]`**
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §5.1.1), which clears the catalog's table at
every period — three repeats of headroom at period 1, one at period 6. So this route's selection is
serviceable, and a run that lowers its routing floor below the catalog's stops rather than
under-reporting.

**Not production's catalog.** `src/ssr/catalog/` builds one via the external `trf-mod` binary and is
frozen production; ng does not depend on it, and
[`typed_regions.md`](typed_regions.md) is explicit that it is a comparison oracle and never an ng
dependency.

**The failure mode a shared artefact introduces is a stale one**, and the mitigation is load-bearing
rather than advisory: because the kept set is a pure function of the catalog's contents, a catalog
built from a different reference selects a *different* set. The catalog is therefore opened with
`open_checking_against_reference`, which compares its per-contig MD5s against digests **recomputed
from the FASTA in this run** — nothing is trusted because it was written down — and names the contig
that differs rather than the genome. What travels with a sample's records is the file's build settings
(§5.1), so the parameters fit refuses to pool across two samples that read different catalogs.

### 3.4 It stays identical in every sample

§1's property survives stratification, and it is worth checking rather than assuming. Everything the
rule reads is fixed before any sample is walked:

- **the catalog** — a function of the reference and the region-typing parameters, computed once;
- **the stratum**, a column of that catalog, read off the reference tract and never from reads;
- **the seed and the cap**, run inputs.

Nothing depends on coverage, read length, traversal order or thread count. Two samples run against the
same catalog, seed and cap keep the identical loci.

### 3.5 The price: the kept set is a random sample *within* a stratum, not of STR loci as a whole

**This is the one real cost of stratifying, and it survives the move to §3.2's rule** — the sample is
now genuinely random inside each stratum, but the strata are still represented equally rather than in
proportion. It splits cleanly.

- **Anything fitted *within* a stratum is unaffected.** The four slippage numbers and the STR path's
  substitution error are estimated from the loci of one stratum, and how many loci other strata
  contributed cannot reach them. **These are the parameters the stratification exists to serve.**
- **Anything pooled *across* strata is biased.** Long-tract strata are rare in the genome and
  over-represented here, and they are also the most variable — so an unweighted statistic over the
  whole kept set comes out too high. **STR diversity is exactly such a statistic**
  ([`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §3).
- **So reweight, and store what makes reweighting possible.** For every stratum, record **how many
  loci the analysed regions hold** and **how many were kept**. STR diversity becomes a mean over
  strata weighted by the true counts. That is one line of arithmetic and **it is silent when
  omitted**, which is why the counts are a stored field
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §5) rather than a note
  here.

**The true counts are `count_loci_per_stratum`'s** (§3.3), answered over the same `region` and in the
same pass as the sample, so keeping them costs a tally beside the heaps rather than a second
traversal.

### 3.6 What selection cannot fix

**A stratum holding fewer loci than the cap keeps all of them and is still thin.** Nothing about
choosing loci helps a stratum with eleven of them; it borrows its value from its neighbouring repeat
counts, exactly as in the per-stratum route
([`parameter_prepass_ssr.md`](parameter_prepass_ssr.md) §4.3).

---

## 4. The size knobs, and what they buy

**Two knobs, one per path, and they have different shapes.**

| path | knob | total loci |
|---|---|---|
| generic | a target position count | the target |
| STR | a per-stratum cap | `cap × non-empty strata`, minus every stratum that holds fewer |

**Memory is linear in the total and precision goes as its square root.** Halving a budget costs about
1.4× the standard error; doubling it buys about 0.7×.

### 4.1 The accuracy does not plateau — usefulness does

**Expecting a plateau is the trap.** Precision keeps improving as `1/√loci` for as long as loci are
added; nothing in the statistics flattens out. **What plateaus is usefulness**: once a parameter is
pinned finer than the caller can feel, more loci buy nothing.

### 4.2 A low-heterozygosity sample does not need more sites — and the reason to say so is that the opposite is the natural guess

**Heterozygosity is a proportion, so measure its error on the scale it lives on: `[0, 1]`.** With `n`
loci and a true rate `p`, the standard error is `√(p(1−p)/n)`, and because `p` is small that is
essentially **`√(p/n)` — which *falls* as the sample gets less heterozygous.** The same budget
estimates a low rate more tightly in absolute terms, not less.

The tomato cohort's 63 samples span **0.149 to 4.483 heterozygotes per kilobase**, mean 1.050, median
0.865 (`benchmarks/tomato1/results/ours/cohort/het_baseline_flat_eps.txt`), so this is checkable
against the range we actually have:

| sample | heterozygous loci in 2 M | standard error | 95% interval |
|---|---:|---:|---|
| at the cohort median, 0.865 /kb | ~1,730 | 2.1 × 10⁻⁵ | 8.2 – 9.1 × 10⁻⁴ |
| at the cohort floor, 0.149 /kb | ~300 | 8.6 × 10⁻⁶ | 1.3 – 1.7 × 10⁻⁴ |

> **⚠ Both intervals assume the two million sites are independent draws, and they are not** (§2). The
> effective count is the number of independent stretches of genome the panel carries: if a 100 kb
> stretch behaves as one draw, two million positions carry 8,000 and both intervals widen sixteen-fold
> — which puts the floor sample's across zero, at ±1.4 × 10⁻⁴ around an estimate of 1.49 × 10⁻⁴. The
> true factor is between about 3 and 16 depending on how far linkage reaches in this panel, **which is
> measurable from the cohort calls that already exist and is not measured here.** What survives
> unchanged is the *comparison* the section is making, since both rows widen by the same factor.

**The inbred sample's interval is narrower**, and the same holds a further order of magnitude down.
Reading it the other way round — 2% of the estimate against 6% — is a *relative* error, and nothing
here consumes one: the caller multiplies a prior, the diversity divides, and both work on the number
itself.

**So: no rule of the form `positions = target / min(Hobs)`, and the budget does not move because a
panel is autogamous.** An earlier version of this section said it did. Two million positions is set by
the cross-sample statistics ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)
§5), and a Bayesian reading says the same thing the frequentist one does: a Beta posterior after 2 M
loci and 300 heterozygous ones is already concentrated, and a few hundred successes determine a small
proportion perfectly well.

**What *is* a real problem on such a sample is a different quantity entirely: how much artefact the
real heterozygotes are buried under.** Of the two million positions, at three reads and an error rate
of 0.001, about **6,000 show a spurious alternative read**; the noisy-locus class — one locus in 110
disagreeing with the reference at 5% rather than 0.19%
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2) — contributes about
**2,500** more. Against those, a sample at the cohort floor has about **260** heterozygous loci that
show an alternative read at all.

**A third contributor is real and, measured, it is the smallest of the three.** A duplication the
reference does not carry shows about half its reads disagreeing at every position where the copies
differ; the generic path's fit asks for such a population at **0.42% and 0.49% of sites** on two
tomato samples and is refused, because the class that would hold it cannot be widened without eating
real heterozygotes ([`parameter_prepass_generic.md`](parameter_prepass_generic.md) §2.1;
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2 is what this route does
instead).

**An earlier version of this paragraph called it the largest contributor at 1,700 to 8,400 of two
million positions, and a genome walk over eight tomato alignments makes it 150 to 590**
([`../reports/duplicated_locus_probe_2026-08-12.md`](../reports/duplicated_locus_probe_2026-08-12.md)).
The old figure read the fitted class weight as a count of positions showing an alternative read.
Those are two quantities: 0.6% to 3.2% of positions sit in a window near two copies, which is what
the weight measures, but **only 0.3% to 1.2% of those positions read near half**, because a
duplication is silent wherever its two copies agree.

**So the budget is about 6,000 error positions, 2,500 noisy ones, a few hundred duplicated ones, and
260 real heterozygotes at the cohort floor. At most 3 in every 100 positions carrying an alternative
read are really heterozygous, and the other 97 are noise** — the conclusion is unchanged, and it now
rests on the two terms that were always the large ones. That is the number to plan for, and it is not
a counting problem:

- **More sites do not fix it.** The signal and the background grow together, so the *ratio* is
  unchanged. Extra loci shrink the sampling error, which was never the binding constraint.
- **It is a bias problem, and bias does not shrink with `n`.** The estimate is what is left after the
  background is subtracted, so an error of a few percent in the background's fitted weight is a
  hundred-percent error in the heterozygosity. At the cohort median the same background error is a
  few percent of the answer; at the floor it is the whole of it.

**This is the case the route is built for, which is why it is stated here rather than treated as a
defeat.** The background is fitted from every sample and every locus at once, so it is pinned by the
whole panel rather than by the one inbred sample; and because a mismapped locus is noisy in *every*
sample, the offending loci are identifiable individually rather than merely averaged over (§4.4).
Whether that is enough to measure 0.149 per kilobase apart from the artefact floor is **unmeasured,
and it is the sharpest question these three documents carry** (§6.2).

**And it lands on inbreeding rather than on genotype calls.** A heterozygosity prior wrong by a factor
of two at 1.5 × 10⁻⁴ flips nothing — the posterior odds for a site showing one alternative read of
three stay near 10⁻³ either way. But `F_hom_excess` is `1 − Hobs/Hexp`
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §5), so where it is near 1 the
whole of `1 − F_hom_excess` is `Hobs/Hexp`: **the heterozygosity *is* the estimate**, and a bias in it
passes through at full size. What the caller's genotype prior multiplies is `1 − F_autozygosity`, a
different number — but the two run together on an autogamous panel, which is precisely the panel where
this bias is largest, so the exposure is real whichever of the two a caller is handed.

### 4.3 The sweep, run — and two million positions is a contamination budget, not an estimates budget

**Measured 2026-08-12** (`examples/ng_joint_fit_harness.rs`, fifty samples, three reads a site, a
fresh draw at each budget, three starting points, best taken). Errors against the drawn truth:

| positions | segregating | clean error rate | noisy error rate | noisy-locus share | `Hexp` |
|---:|---:|---:|---:|---:|---:|
| 5,000 | 39 | −1.1% | −2.0% | +15.9% | +14.0% |
| 20,000 | 149 | −1.4% | +2.8% | −13.3% | −3.5% |
| 80,000 | 627 | +0.2% | −0.4% | **−1.3%** | −7.0% |
| 320,000 | 2,614 | −0.4% | −1.3% | −1.0% | **+1.6%** |
| 1,280,000 | 10,208 | −0.0% | +0.2% | −0.9% | −0.9% |

**The three parameters degrade at completely different rates, which is why §4.3's older wording
insisted on one row per parameter.**

- **The error rates are finished at five thousand positions** — within about 2%, and 256 times the
  budget buys another 1%. They are measured from *reads*, and five thousand positions at fifty
  samples and three reads is three quarters of a million read observations. **Nothing about the site
  budget is set by them.**
- **The noisy-locus share needs about eighty thousand.** It is a property of *loci* rather than of
  reads, so it wants loci: wrong by a sixth at five thousand, within a hundredth at eighty thousand.
- **The diversity is the slowest, and it tracks the segregating count rather than the budget.** It
  scatters by several percent until a few thousand sites segregate, and settles near 1% at ten
  thousand.
- **Inbreeding needs about twenty thousand.** Against a drawn truth of 0.600 the parameters fit returns 0.563 at
  five thousand positions and lands within a percentage point at every budget above it; against a
  truth of zero it returns 0.052 at five thousand and 0.000 above. **The low-budget failure is an
  inbreeding coefficient invented out of scatter**, and it goes both ways.

*The whole table repeats at an inbreeding of 0.6 with every column within a percentage point of the
outbred one, so none of these budgets depends on how inbred the panel is.*

**So for everything in this table, three hundred and twenty thousand positions is enough at fifty
samples — a sixth of the two million.** What keeps the budget at two million is the parameter that
is *not* in this table: contamination wants about ten thousand segregating markers
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.4), and this panel yields
10,208 of them at 1.28 M positions. **Two million is a contamination budget.** A run that does not
need `α` to a couple of percent can have its records six times smaller, and that is now a knob with a
measured meaning rather than a number nobody had checked. **§4.3.4 says how that knob is set once and
turned down later**, and [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.4 says
why the contamination end of it is worth turning *up* now that the records are on disk.

### 4.3.1 A bigger panel does not buy fewer loci — MEASURED

The same sweep at **two hundred** samples, four times the panel:

| positions | noisy-locus share, 50 samples | at 200 | `Hexp`, 50 | at 200 | clean rate, 50 | at 200 |
|---:|---:|---:|---:|---:|---:|---:|
| 5,000 | +15.9% | **+16.0%** | +14.0% | −9.0% | −1.1% | +0.7% |
| 20,000 | −13.3% | **−12.4%** | −3.5% | +3.3% | −1.4% | −0.0% |
| 80,000 | −1.3% | +4.5% | −7.0% | −4.0% | +0.2% | +0.5% |
| 320,000 | −1.0% | −2.8% | +1.6% | +0.4% | −0.4% | −0.2% |

**The noisy-locus share is unmoved by the panel** — +16% at five thousand positions whether the
cohort holds fifty samples or two hundred, and settling only as the loci accumulate. That is what it
must do: `w` is *the share of loci that are noisy*, so it is counted in loci and **cannot be bought
with samples**. §4.4 says the same thing about the other direction, and this is its mirror: **the
budget buys what loci buy, and the panel buys what samples buy, and neither substitutes for the
other.**

**What the larger panel does buy** is the read-driven parameters at small budgets: the clean rate
goes from −1.1% to +0.7% at five thousand positions, and the inbreeding coefficient from 0.052 to
0.004 against a truth of zero. Both are fitted from reads, and four times the samples is four times
the reads.

**So the locus budget is set by two things, and a cohort of thousands relieves neither**: the
noisy-locus share, which needs about eighty thousand positions, and contamination, which needs its
segregating markers.

### 4.3.2 What a thousand samples does for the marker count: about a seventh — MEASURED

The counts above are sites segregating **in the population**, a property of the truth that does not
move with the panel. What contamination can actually use is markers segregating **in the panel** —
where the cohort's own chromosomes differ — and that does grow with the panel, because a bigger
sample reaches further down the rare tail. **How far it grows is a property of the population's
low-frequency tail, so it had to be measured rather than argued.**

| panel | of the population's segregating sites, the share the panel sees |
|---:|---:|
| 50 samples | 77–82% |
| 200 | 82–85% |
| 1,000 | **89–92%** |

**Twenty times the panel buys about fourteen percent more usable markers, not an order of magnitude.**
Under this truth's density — `Beta(0.3, 1.2)`, the rare-allele pile-up a neutral population has —
about a quarter of segregating sites sit below a frequency of one in a hundred, and most of those are
already visible to fifty samples. **So the contamination budget falls by roughly a seventh between
fifty samples and a thousand, and a cohort of thousands does not make the site budget collapse.**

*An earlier note in this section guessed the other way*, on the reasoning that a panel of two
thousand reaches alleles at one in four thousand where fifty reaches one in a hundred. That is true
and it is not the point: what matters is how much *mass* the density puts down there, and under a
neutral tail it is not much. **On a population with a heavier rare tail — a recent expansion, or a
much larger effective size — the same measurement would come out differently**, which is a reason to
re-run it on the real cohort rather than to carry this number as a constant.

*And one refinement of §4.3.1's claim.* At five thousand positions the noisy-locus share is +16% at
fifty, two hundred **and** a thousand samples, so there the panel buys nothing at all. At twenty
thousand, a thousand-sample panel does better than a fifty-sample one — −4.0% against −13.3% — so
samples are not entirely useless to it. **The budget still has to be counted in loci**; what is wrong
is only the stronger form of the claim.

*And it is at the cohort's own heterozygosity throughout; §4.2's low-heterozygosity case is the
separate question §6.2 carries.*

### 4.3.3 The experiment as it was specified, for the parts still to run

Refit at 2 M, 500 k, 100 k and 20 k positions on the same drawn data and find **the first budget at
which the fitted values move by more than a caller could feel**. Two million was chosen for the
cross-sample statistics and has never been checked against this use.

**Report two things at each budget and nothing else: what it cost at rest, and each parameter's error
against the drawn truth — one row per parameter.** They are the two quantities the knob trades, and
they must not be pooled: the error rate and the per-sample rates degrade as the square root of the
budget, while the spectrum's rare classes empty out altogether, so a single "precision" column would
report a budget as adequate that has already lost the spectrum
([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.1). **Report the
segregating-site count each budget yielded** beside them; it is an outcome of the panel, not a target.

**Run it at the low end of the heterozygosity distribution as well as at the median** — not because
the budget should differ (§4.2 says it should not) but because **the failure there is of a different
kind**: at the median a too-small budget shows up as scatter, and at the floor a mis-fitted background
shows up as a confident wrong number that more sites do not cure. A sweep run only at 1 per kilobase
would report the first and never see the second. **At minimum include a sample drawn at 0.149 per
kilobase** — tomato's measured floor — **and one an order of magnitude below it** for a selfing line.
§6 carries it.

### 4.3.4 The budget is a per-run knob, and a small census is contained in a large one — DECIDED 2026-08-13

**Decision (owner): where the records are written to a file, the genome walk writes the largest census that
will plausibly be wanted, and every smaller run takes a subset of it without rebuilding anything.**

**Both rules nest, which is what makes that safe.** The generic rule keeps a position when
`hash(contig, p, seed) < threshold` (§2), and a smaller target is a smaller threshold — so its set is
contained in the larger one's, position for position. The STR rule keeps the `cap` lowest hashes per
stratum (§3.2), and a smaller cap keeps a prefix of that same ordering. **Neither needs the file read
again**: a run wanting fewer loci recomputes the smaller threshold and skips the entries that fail it.

**Only one direction is free.** Shrinking costs nothing. Growing needs the census rebuilt, which
outside the direct run means one pass over every sample's pileup
([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1). That asymmetry is
the argument for erring high rather than low, and the depth array is five bits a position so erring
high is cheap: two million positions is 1.25 MB per read group.

*How far high is no longer a contamination question.* An earlier revision of this section argued for
about twenty-seven million positions because contamination's noise floor fell as markers were added.
**That floor was two defects and is fixed** ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md)
§3.4.4), so nothing now asks for a census larger than §4.3's own table does — three hundred and twenty
thousand positions for everything measured there, and two million for the margin.

**This corrects the reason §5.1 gives for carrying the generic target count**, though not the
requirement itself. That table says a different target is *"a different set, not a subset of the larger
one"*, and under §2's threshold rule it plainly is a subset. What still forces two samples to agree is
the **indexing**: the depth array stores no coordinates and entry *i* is the *i*-th kept position, so
records written at two targets have different entries at the same offset. The parameters fit must refuse, and the
value stays on the list.

**No reordering follows from any of this.** Records stay in genome order, which is what the genome walk
produces and what a region-sharded merge concatenates. The smaller set is not a contiguous prefix of
the file and does not need to be: it is found by recomputing the hash the selection already uses.

### 4.4 One thing the budget does not buy

The joint parameters fit's distinctive advantage — telling a mismapped locus from a heterozygous one, because a
mismapped locus is noisy in *every* sample
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2) — comes from having **many
samples at one locus**, not many loci. It does not improve as the budget grows and does not degrade as
it shrinks. **The budget buys precision for the pooled rates and nothing at all for that.**

### 4.5 The STR cap turns out to be no cap at all — MEASURED 2026-08-12

[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.3 suspected the whole STR
set would fit with no sampling: there are orders of magnitude fewer STR loci than genome positions.
**It does.** At the STR path's calling floors `[8, 6, 6, 6, 5, 4]` with a 30 bp flank and a 100 bp
satellite cap, tomato SL4.00 holds **462,701 STR loci in 141 strata**
(`examples/ng_joint_loci_probe.rs`) — under a quarter of the two million generic positions, and a cap
above the largest stratum keeps every one of them.

| cap | loci kept | strata capped, of 141 | what it does to the estimate (§6 question 1) |
|---:|---:|---:|---|
| 100 | 8,699 | 68 | below the floor at both depths |
| 500 | 27,698 | 35 | below the floor at both depths |
| 1,000 | 41,271 | 21 | the floor at six reads a site; too small at three |
| **5,000** | **86,688** | **8** | **the floor at three reads a site — the recommended cap** |
| 20,000 | 157,752 | 3 | above the floor; buys nothing measured |
| none | 462,701 | 0 | above the floor; buys nothing measured |

*The last column is measured, 2026-08-13*
([`../reports/str_stratum_size_sweep_2026-08-13.md`](../reports/str_stratum_size_sweep_2026-08-13.md)),
*and the cap is set by the shallow case because the number that fixes it — the concentration — does not
improve with depth. §6 question 1 says which of the five fitted numbers breaks first and why.*

**How many strata can be fitted at all is a property of the analysed regions, not of the reference —
MEASURED 2026-08-13**
([`../reports/str_fit_on_real_records_2026-08-13.md`](../reports/str_fit_on_real_records_2026-08-13.md)
§5). The whole tomato reference holds 462,701 tracts in 141 strata. A **452 kb region set** on the
human reference holds **216 tracts in 32 strata**, of which exactly one could be fitted on its own; and
on the 63-accession tomato run, **65 of 71 strata could say nothing from their own tracts**, the six
that could holding 88% of them. The consequence is not the count but what the fit then does: every
stratum below the floor borrows from its neighbours, and where nearly all of them are below it, they
pool each other and the repeat-count axis is averaged away
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §4).

**So a run reports, before fitting, how many strata clear the floor on their own.** It costs one pass
over counts the selection already holds, and it is the difference between a user reading a flat
repeat-count axis and a user knowing the axis was flattened.

**The distribution is what makes a cap look attractive and then unnecessary.** One stratum — period 1
at 8 repeats — holds 217,812 loci, 47% of the total, and the next two hold another 32%; **68 strata
hold fewer than a hundred loci each**. So a cap does almost all of its work on three strata, and
those three are the ones where the parameter is already best determined.

**Three things follow, and the third is the one that matters.**

- **No sampling means no reweighting.** §3.5's per-stratum weights are the correction for keeping
  equal numbers from unequal strata; keep everything and the strata are represented in proportion by
  construction. The counts stay a stored field (§3.5) because a *future* run may cap, and because
  they are free — but the arithmetic they enable does not have to fire.
- **The thin strata are unchanged by any of this.** Selection was never their problem (§3.6): 68
  strata below a hundred loci borrow from their neighbouring repeat counts whatever the cap is.
- **The memory bill moved.** Keeping 462,701 loci in every sample makes the STR records the *larger*
  half of what the cohort holds, not the smaller
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6). **That is the
  reason to keep the cap mechanism**, and the first reason it has ever had to fire.

---

## 5. Cross-cutting concerns

**Concurrency.** The generic rule is a pure function of position, so a region-sharded genome walk needs no
communication: each shard selects within its own region and merging is concatenation in position
order. **The STR rule is order-independent by construction** (§3.2), so a sharded enumeration of the
catalog needs none either: each shard keeps the lowest `cap` hashes it saw per stratum, and merging is
taking the lowest `cap` of the union. That the result equals the single-threaded one is a property of
the rule, not something the merge has to arrange.

**Memory.** Neither the generic positions nor the STR loci are written down as coordinates: the
generic set is reproducible from its hash rule, and the STR set is reproducible from the catalog, the
seed and the cap. Selecting the STR set holds `cap` values per stratum, a bounded heap — a few hundred
strata at a cap in the thousands is a handful of megabytes. What travels with a sample's data is the
*inputs* (§5.1).

**Compute.** The catalog is built once per (reference, scan settings) rather than once per sample —
measured at 30 s over tomato's 795 MB reference and 110 s over GRCh38, inside the digest pass. Reading
it to select is a hash and a heap push per surviving locus, in one forward pass.

**Errors.** A run whose catalog differs from another's has different loci in it, so the two samples did
not select the same set. That is not a degraded estimate but a meaningless one, and the parameters fit must be
able to refuse rather than average. **The catalog's own staleness is caught before that** — its
per-contig MD5s are checked against digests recomputed from the FASTA when it is opened (§3.3) — so
what §5.1 guards is the different failure: two samples that read compatible files and asked them
different questions.

### 5.1 What must travel with a sample so the parameters fit can check

Seven values identifying what was asked for, checked for agreement across every sample before anything
is pooled — and then an eighth, below, that checks what came back:

| value | why a mismatch is silent |
|---|---|
| the selection **seed** | a different seed selects a disjoint set; both look well-formed |
| the **reference** digest | coordinates mean different things, so "the same position" is not |
| the **analysed region set** digest | the likeliest to differ by accident, because a BED feels like a runtime convenience. It is not: it defines what population of sites the estimate describes |
| the catalog's **build settings and scoring weights** | `build_settings()` — the floors, period range and flank the file was built at, plus the two Ruzzo–Tompa weights. They decide which tracts exist to be sampled at all, and the weights are not a filter: a different weighting is a different set of tracts, not a subset of one |
| the **STR routing criteria** this run asked for | the copy floors, purity floor, satellite cap and bundle radius passed to `str_loci`. Two samples that filtered the same file differently hold different loci |
| the **generic target count** | a smaller target *is* a subset of a larger one under §2's threshold rule (§4.3.4) — what breaks is the indexing, since the depth array holds no coordinates and entry *i* is the *i*-th kept position, so two targets put different loci at the same offset |
| the **STR per-stratum cap** | a smaller cap *is* a subset of a larger one under §3.2's rule — but of a different size, so two samples still hold different loci |

**The catalog's settings replace "the region-typing parameters" that
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2 lists, and the pair above is strictly
stronger**: the build settings say what the file contains and the routing criteria say what this run
took out of it, and a mismatch in either changes the loci. **The catalog's own identity is not on this
list, and deliberately** — it is checked against MD5s recomputed from the FASTA when the file is
opened (§3.3), so by the time a sample has records it has already proved which reference it read. What
two samples can still disagree about is what they *asked* the file, which is what travels.

The first two rows are already required by
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2; the last four are this document's
addition and fail in exactly the same silent way.

**All seven check the question, and none of them checks the answer.** They say two runs were *asked*
for the same loci. They cannot see a hash function that changed between versions, a threshold computed
in 64 bits on one machine and 128 on another, a catalog read in a different order by a sampler whose
order-independence has regressed, or a genome walk that filled its array from the wrong end — every one of
those leaves all seven agreeing and the kept sets different.

**So an eighth value travels, and it is the only direct one: a digest of the loci actually kept,
computed as the records are filled rather than by re-deriving the rule.** Thirty-two bytes per sample,
plus one digest per megabase so a mismatch names the block it happened in — 6.4 kB on tomato, half a
percent of the record. A digest produced by running the selection a second time proves only that the
selection is deterministic, which nobody doubts; it must witness the array that was written.
[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.2 carries the reasoning
and the test that distinguishes the two.

---

## 6. Open questions

1. **How many STR strata are there, and how many loci does each hold?** — **CLOSED 2026-08-12: 141
   strata, 462,701 loci on tomato at the calling floors** (§4.5,
   `examples/ng_joint_loci_probe.rs`). Three decisions waited on it and all three now have an answer.

   - **What the cap should be**: none, for accuracy. A cap set above 217,812 keeps every locus, and
     the only reason to set one lower is the memory bill §4.5's last point names.

     **REOPENED and CLOSED 2026-08-13: the cap is 5,000 tracts**
     ([`../reports/str_stratum_size_sweep_2026-08-13.md`](../reports/str_stratum_size_sweep_2026-08-13.md)).
     That is the floor at tomato's three reads a site and not a margin above it — one step down, at
     1,000 tracts, one of the five fitted numbers has already lost its footing. At six reads a site
     the floor is 1,000, and **the cap is set by the shallow case** because the number that fixes it
     does not improve with depth (below). Tomato then keeps 86,688 of its 462,701 tracts, with 8 of
     its 141 strata capped at all, and the largest section of a records file is 50 kB a sample —
     50 MB across a thousand.

     **Which of the five sets it is the part worth carrying, because it is not the one anyone would
     guess.** The slippage level — the number the estimator exists to produce — is the most durable of
     the five, still within 3.9% at 250 tracts. What breaks first is *how fast two-repeat slips fall
     off against one-repeat slips*, which rests on the fifth of the slipped reads that slipped by two:
     roughly 225 reads at 250 tracts. What breaks second is *the concentration*, how monomorphic the
     stratum's tracts are, **and that is the one the cap exists to supply**: it is counted in tracts
     rather than in reads, so doubling the depth halves the scatter of the read-driven numbers and
     leaves it untouched — 14.3% against 14.2% at 100 tracts. A deeper cohort does not get to keep
     fewer tracts.

     *Why the question was reopened, kept because the reasoning still holds:* reading one stratum at a
     time bounds what the parameters fit holds, and the bound **is** the cap
     ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.2). Memory does not
     choose the number — a tract costs about ten bytes a read group, so at a thousand samples a cap of
     5,000 and one of 20,000 are 50 MB and 200 MB for the largest section, and both are affordable.
     What chooses it is the estimate, which is why the sweep was needed.

     **Two parts of that sweep were not run, and each could raise the floor.** It drew **one stratum
     shape** — three length classes at a concentration of 0.5 — and a more nearly monomorphic stratum
     carries less signal per tract for exactly the number that sets the cap. And it held the panel at
     **twenty samples** throughout, so how tracts and samples trade against one another is unmeasured.
     **Neither touches the thin strata**: the 68 below a hundred tracts borrow from their neighbouring
     repeat counts whatever the cap is (§3.6), and where that borrowing has to start is still nobody's
     measurement.

     *One consequence of capping that is easy to forget:* §3.5's per-stratum reweighting exists for
     runs that sample, and a capped run samples. The stored counts it needs are already there, but the
     arithmetic that has never had to fire would begin to.
   - **Whether any sampling is needed at all**: no, so §3's reweighting never has to fire — though
     its stored counts stay, being free and being what a capped run would need.
   - **Whether the comparison against the per-sample route is a comparison**: on this path, yes and
     fully. The joint route holds *the same STR loci* the per-stratum histogram holds and remembers
     which was which, so §8's second measurement is like-for-like here where it is not on the generic
     path.

   *What was measured before was a different number* — the catalog holds 6.4 M repeats over tomato's
   795 MB reference and 23.6 M over GRCh38, but those are rows at the build floors
   `[5, 5, 4, 4, 4, 3]`. The calling floors, the purity floor, the satellite cap and bundling take
   6.4 M down to 462,701: **a factor of about fourteen**, which is the number nothing had.

   **Still open, and it is the human half**: the same count on GRCh38, which the current catalog file
   cannot answer because it was written in an older header format. It matters for the GIAB arm of
   [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §8's third measurement and for
   nothing else.
2. **Can a sample at 0.149 heterozygotes per kilobase be told apart from the artefact floor?** —
   OPEN, and **the sharpest question in these three documents** (§4.2). At least 97 in every 100
   positions carrying an alternative read there are noise, so the answer turns on how well the panel
   pins the background and how many of the offending loci it can name individually (§4.4) — **not on
   the site budget, which changes the ratio not at all**. *Leaning:* yes at tomato's floor and unknown
   an order of magnitude below it; the failure mode is a **confident** wrong number rather than a
   noisy one, because a mis-fitted background is bias and bias does not shrink with more loci.
   **Settled by:** §4.3's sweep run at the low end, against drawn genomes whose true heterozygosity is
   known.

   **It is not an incidental output on such a panel** (§4.2): `1 − F_hom_excess` is `Hobs/Hexp` where
   that coefficient is near 1, so the heterozygosity is the inbreeding estimate under another name —
   and on an autogamous panel it tracks the `F_autozygosity` the caller's genotype prior multiplies.
3. **Where does the generic budget start to matter?** — **MEASURED 2026-08-12 at fifty samples: not
   until well below two million, and each parameter has its own answer** (§4.3). The error rates are
   within 2% at **five thousand** positions, the noisy-locus share within a hundredth at **eighty
   thousand**, and the diversity near 1% once a few thousand sites segregate, which takes about
   **three hundred and twenty thousand**. **The two-million budget is set by contamination alone**,
   which wants ten thousand segregating markers and gets 10,208 at 1.28 M positions.

   **Still open, and both are cheap**: the same sweep at a larger panel — more samples means more
   sites segregate somewhere, so the budget should fall further — and the same sweep on a sample at
   tomato's heterozygosity floor, which is question 2 and a different failure mode entirely.
4. **What fraction of STR loci vary across a cohort?** — OPEN, and it is the other half of question 1:
   it decides whether a cap set to keep everything actually delivers the ~10,000 varying loci the
   cross-sample statistics want ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.3).
   **The same low-heterozygosity concern applies here and is smaller**, because a repeat tract mutates
   orders of magnitude faster than a base does, so an autogamous panel that is nearly monomorphic for
   substitutions still segregates at its STR loci.

---

## 7. How we know it works

*The STR rule's own properties — exactly `cap` kept, spread rather than a prefix, order-independent,
refusal when a reader is too permissive — belong to
[`repeat_catalog.md`](repeat_catalog.md) §10, which ships with the sampler and is where they are
asserted. They are not restated here. What follows is this consumer's own.*

1. **Every sample selects the same loci.** Run the selection over samples with different coverage,
   different read lengths, and a region-sharded genome walk at several thread counts; the set must be
   identical every time, and identical to the set computed directly from the seven values of §5.1.
   **Include a `--regions` case**: the generic positions must all lie inside the BED and their count
   must come from the region set's length rather than the genome's — the arithmetic most likely to be
   written against the contig table by reflex — and the STR counts and sample must come from
   `count_loci_per_stratum` and `sample_loci_per_stratum` asked over that same region, not over the
   whole reference.
2. **The generic selection is unbiased.** On synthetic data with a known frequency spectrum, the
   spectrum over the chosen positions must match the one over every position, within its own error.
   If the selection ever came to depend on the data, this is where it shows.

   *Two properties of the rule itself are cheaper to check and are checked on real references
   instead of a fixture: that the realised count lands within `√target` of the target, and that the
   gaps between kept positions are geometric. Measured — tomato returns 2,002,505 positions for a
   target of 2,000,000 with gaps at 41 / 271 / 900 bp against a geometric prediction of 41 / 271 /
   900; GRCh38 returns 1,999,981 with 155 / 1,019 / 3,385 against 155 / 1,019 / 3,385. The kept
   positions fall across contigs as their lengths imply: chi-square 19.8 on 12 degrees of freedom on
   tomato, 127.9 on 127 on GRCh38.*
3. **No kept position sits on an ambiguous base** (§2). Build the selection over a reference holding
   an `N` run and require every kept position to be an `A`, `C`, `G` or `T`; **and require the
   realised count to come from the masked length rather than the contig's**, which is the same
   arithmetic error the `--regions` case sets and the one an implementation is most likely to make
   twice. On GRCh38 an unmasked rule puts 106,423 of 2,000,000 positions inside a gap, so a broken
   mask is visible without a fixture.
4. **The reweighting undoes the stratification.** STR diversity computed from the kept loci with the
   stored per-stratum weights must match the diversity computed from **every** STR locus, and must
   differ from the unweighted version. The second half is what makes the first half matter, and the
   whole thing is silent when the reweighting is omitted (§3.5).
5. **A larger cap contains a smaller one.** Selecting at `cap` and at `2·cap` must give nested sets.
   That is the sampler's property, asserted there; what this document needs from it is that the
   downward budget sweep of §4.1 is a sequence of subsets rather than unrelated draws, so **assert it
   at the two budgets the sweep actually uses**.
6. **The two sets partition, and neither leaks.** No generic position falls inside a locus the STR
   criteria keep. Run it with the routing floors moved and confirm both sets change together — the
   partition follows the criteria handed to the catalog rather than a second copy of its rules.
7. **A run more permissive than the catalog stops.** Ask for a copy floor below `[5, 5, 4, 4, 4, 3]`
   at any period, or a flank below 15 bp; the run must fail naming the axis and both values, **not
   proceed on a short list** — which would be a wrong per-stratum total that nothing downstream could
   notice (§3.3).
8. **A mismatch across samples is refused.** Two samples selected under different seeds, references,
   region sets, catalog build settings, routing criteria, target counts **or caps** must produce an
   error and not an average (§5.1).
9. **The kept-loci digest catches what the seven values cannot** (§5.1). Change the selection's answer
   while leaving all seven inputs identical — swap two kept loci, or drop one and add another — and the
   digest must change, the per-megabase digest must name the block, and the parameters fit must refuse. **Then
   check the check**: a digest re-derived by running the selection again passes this test unchanged,
   which is why §5.1 requires it to be computed where the records are written.
