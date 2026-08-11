# ng — the joint fit: which loci are kept

*Design spec, 2026-08-10, revised 2026-08-11. One of three documents covering the **joint fit**, ng
step 4's second route to every parameter it emits; read
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) first — it says what the route is,
what it produces and why it exists. This one settles **which loci every sample keeps evidence at**,
and nothing else. What is recorded at each is
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md).*

***The STR half of it is built.*** *The repeat catalog ships the per-stratum sampler this document
asked for ([`repeat_catalog.md`](repeat_catalog.md), `src/ng/repeat_catalog/`), so §3 now records why
the rule is that one and what using it obliges a consumer to do, rather than proposing it. The generic
rule (§2) has no code.*

***It changes one decision in [`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md)***
*§3 — how STR loci are chosen (§3 below) — and that document carries a note saying so.*

---

## 1. What this is, in one paragraph

**The joint fit reads a small set of loci at which every sample kept its raw evidence, instead of the
per-sample summaries the other route builds.** For that to work at all, every sample has to keep
evidence at *the same* loci — and the samples are walked separately, on different machines, at
different times, with no sample able to see what any other chose. So the set cannot be negotiated or
handed round. **This document is the rule that lets each sample arrive at the identical set on its
own**, for the two kinds of locus the two paths care about: ordinary positions for the SNP/indel
path, and repeat tracts for the STR path.

**It is a module of its own because it is a pure function of the run's inputs** — the reference, the
analysed regions, the repeat catalog, a seed, and a couple of caps. It touches no read, so it can be
built and tested with no alignment file in sight: run it twice and compare the lists.

**Half of it is already built.** The STR side is two calls on the repeat catalog
([`repeat_catalog.md`](repeat_catalog.md), `src/ng/repeat_catalog/`), which ships the per-stratum
sampler; the generic side is a hash rule this document states in full (§2). What is left here is the
policy — which caps, which seed, what has to travel with a sample, and what the stratification obliges
a consumer to do afterwards (§3.5).

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

> keep position `p` if `p` lies in the analysed regions and `hash(contig, p, seed) < threshold`

- **The domain is the analysed regions, not the genome.** With no `--regions` BED that is the whole
  reference; with one, it is the reference intersected with the BED. Selecting genome-wide under a BED
  would leave nearly every chosen position unvisited.
- **The threshold sets the count**, from the analysed length. The hash is a 64-bit value, so it takes
  `2^64` values and a target of `n` positions out of `analysed_length` is
  `threshold = 2^64 · n / analysed_length`. **Nothing about that range is measured or discovered — it
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
- **Scattered positions, never contiguous blocks.** Sites within a block share a genealogy, so a block
  of *k* linked positions carries far less independent information than *k* scattered ones. The
  estimators downstream treat sites as independent, and scattering is what makes that nearly true.
- **Repeat tracts are excluded**, using region typing's delimitation — their variability would distort
  every substitution statistic computed from the set.

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
(§5.1), so the fit refuses to pool across two samples that read different catalogs.

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
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §4) rather than a note
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

**So roughly 3 in every 100 positions carrying an alternative read are really heterozygous, and the
other 97 are noise.** That is the number to plan for, and it is not a counting problem:

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
three stay near 10⁻³ either way. But the homozygote-excess `F` is `1 − Hobs/Hexp`
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §5), so where `F` is near 1 the
whole of `1 − F` is `Hobs/Hexp` — and `1 − F` is what the caller's genotype prior multiplies. **On an
autogamous panel the per-sample heterozygosity is the inbreeding estimate wearing a different name**,
which is what makes a background-driven bias in it worth this much care.

### 4.3 The experiment: sweep downward, and at the low end

Refit at 2 M, 500 k, 100 k and 20 k positions on the same drawn data and find **the first budget at
which the fitted values move by more than a caller could feel**. Two million was chosen for the
cross-sample statistics and has never been checked against this use.

**Run it at the low end of the heterozygosity distribution as well as at the median** — not because
the budget should differ (§4.2 says it should not) but because **the failure there is of a different
kind**: at the median a too-small budget shows up as scatter, and at the floor a mis-fitted background
shows up as a confident wrong number that more sites do not cure. A sweep run only at 1 per kilobase
would report the first and never see the second. **At minimum include a sample drawn at 0.149 per
kilobase** — tomato's measured floor — **and one an order of magnitude below it** for a selfing line.
§6 carries it.

### 4.4 One thing the budget does not buy

The joint fit's distinctive advantage — telling a mismapped locus from a heterozygous one, because a
mismapped locus is noisy in *every* sample
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2) — comes from having **many
samples at one locus**, not many loci. It does not improve as the budget grows and does not degrade as
it shrinks. **The budget buys precision for the pooled rates and nothing at all for that.**

### 4.5 The STR cap may turn out to be no cap at all

[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.1 suspects the whole STR
set may fit with no sampling: there are orders of magnitude fewer STR loci than genome positions, and
a far larger fraction of them vary. **The mechanism costs nothing either way** — a cap high enough to
admit every locus keeps every locus, and `sample_loci_per_stratum` is the same call — so this is a
number to set from §6.1's measurement rather than a design to revisit.

---

## 5. Cross-cutting concerns

**Concurrency.** The generic rule is a pure function of position, so a region-sharded walk needs no
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
not select the same set. That is not a degraded estimate but a meaningless one, and the fit must be
able to refuse rather than average. **The catalog's own staleness is caught before that** — its
per-contig MD5s are checked against digests recomputed from the FASTA when it is opened (§3.3) — so
what §5.1 guards is the different failure: two samples that read compatible files and asked them
different questions.

### 5.1 What must travel with a sample so the fit can check

Seven values, checked for agreement across every sample before anything is pooled:

| value | why a mismatch is silent |
|---|---|
| the selection **seed** | a different seed selects a disjoint set; both look well-formed |
| the **reference** digest | coordinates mean different things, so "the same position" is not |
| the **analysed region set** digest | the likeliest to differ by accident, because a BED feels like a runtime convenience. It is not: it defines what population of sites the estimate describes |
| the catalog's **build settings and scoring weights** | `build_settings()` — the floors, period range and flank the file was built at, plus the two Ruzzo–Tompa weights. They decide which tracts exist to be sampled at all, and the weights are not a filter: a different weighting is a different set of tracts, not a subset of one |
| the **STR routing criteria** this run asked for | the copy floors, purity floor, satellite cap and bundle radius passed to `str_loci`. Two samples that filtered the same file differently hold different loci |
| the **generic target count** | a different target is a different set, not a subset of the larger one |
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

---

## 6. Open questions

1. **How many STR strata are there, and how many loci does each hold?** — OPEN, and it is the first
   thing to settle, because three decisions wait on it: what the cap should be (total loci is
   `cap × non-empty strata`), whether any sampling is needed at all (§4.3), and how much of the joint
   fit's comparison against the per-sample route is even a comparison — if every locus is kept, the
   joint route holds the same STR loci *and* remembers which was which.

   **It is now one call rather than an experiment**: `count_loci_per_stratum` at the STR path's
   calling floors, over the tomato reference. *What is measured already is a different number* — the
   catalog holds 6.4 M repeats over tomato's 795 MB reference and 23.6 M over GRCh38, but those are
   rows at the build floors `[5, 5, 4, 4, 4, 3]`, not loci at the calling floors `[8, 6, 6, 6, 5, 4]`
   after the purity floor, the satellite cap and bundling have been applied. The second number is
   smaller by an unknown factor and it is the one the cap depends on.

   **If the answer is "they all fit", most of §3 dissolves.** With no cap there is no sampling and no
   reweighting to remember — the strata are represented in proportion because everything is kept. §3
   is what happens when they do not fit, and it costs nothing to keep either way (§4.3), but the
   measurement decides whether it ever fires.
2. **Can a sample at 0.149 heterozygotes per kilobase be told apart from the artefact floor?** —
   OPEN, and **the sharpest question in these three documents** (§4.2). About 97 in every 100
   positions carrying an alternative read there are noise, so the answer turns on how well the panel
   pins the background and how many of the offending loci it can name individually (§4.4) — **not on
   the site budget, which changes the ratio not at all**. *Leaning:* yes at tomato's floor and unknown
   an order of magnitude below it; the failure mode is a **confident** wrong number rather than a
   noisy one, because a mis-fitted background is bias and bias does not shrink with more loci.
   **Settled by:** §4.3's sweep run at the low end, against drawn genomes whose true heterozygosity is
   known.

   **It is not an incidental output on such a panel** (§4.2): `1 − F` is `Hobs/Hexp` where `F` is near
   1, so this is the inbreeding estimate under another name, and inbreeding is what the caller's
   genotype prior multiplies.
3. **Where does the generic budget start to matter?** — OPEN. §4.3 gives the experiment. *Leaning:*
   below two million, since the error rate is already pinned to one part in eighty there and a low
   heterozygosity does not by itself demand more loci (§4.2). What could still force a per-run knob is
   question 2's answer rather than any counting argument.
4. **What fraction of STR loci vary across a cohort?** — OPEN, and it is the other half of question 1:
   it decides whether a cap set to keep everything actually delivers the ~10,000 varying loci the
   cross-sample statistics want ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.1).
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
   different read lengths, and a region-sharded walk at several thread counts; the set must be
   identical every time, and identical to the set computed directly from the seven values of §5.1.
   **Include a `--regions` case**: the generic positions must all lie inside the BED and their count
   must come from the region set's length rather than the genome's — the arithmetic most likely to be
   written against the contig table by reflex — and the STR counts and sample must come from
   `count_loci_per_stratum` and `sample_loci_per_stratum` asked over that same region, not over the
   whole reference.
2. **The generic selection is unbiased.** On synthetic data with a known frequency spectrum, the
   spectrum over the chosen positions must match the one over every position, within its own error.
   If the selection ever came to depend on the data, this is where it shows.
3. **The reweighting undoes the stratification.** STR diversity computed from the kept loci with the
   stored per-stratum weights must match the diversity computed from **every** STR locus, and must
   differ from the unweighted version. The second half is what makes the first half matter, and the
   whole thing is silent when the reweighting is omitted (§3.5).
4. **A larger cap contains a smaller one.** Selecting at `cap` and at `2·cap` must give nested sets.
   That is the sampler's property, asserted there; what this document needs from it is that the
   downward budget sweep of §4.1 is a sequence of subsets rather than unrelated draws, so **assert it
   at the two budgets the sweep actually uses**.
5. **The two sets partition, and neither leaks.** No generic position falls inside a locus the STR
   criteria keep. Run it with the routing floors moved and confirm both sets change together — the
   partition follows the criteria handed to the catalog rather than a second copy of its rules.
6. **A run more permissive than the catalog stops.** Ask for a copy floor below `[5, 5, 4, 4, 4, 3]`
   at any period, or a flank below 15 bp; the run must fail naming the axis and both values, **not
   proceed on a short list** — which would be a wrong per-stratum total that nothing downstream could
   notice (§3.3).
7. **A mismatch across samples is refused.** Two samples selected under different seeds, references,
   region sets, catalog build settings, routing criteria, target counts **or caps** must produce an
   error and not an average (§5.1).
