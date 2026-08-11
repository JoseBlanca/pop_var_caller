# ng — the joint fit: which loci are kept

*Design spec, 2026-08-10. **No code yet — this settles the design.** One of three documents covering
the **joint fit**, ng step 4's second route to every parameter it emits; read
[`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) first — it says what the route is,
what it produces and why it exists. This one settles **which loci every sample keeps evidence at**,
and nothing else. What is recorded at each is
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md).*

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
analysed regions, region typing's output, a seed, and a couple of caps. It touches no read, so it can
be built and tested with no alignment file in sight: run it twice and compare the lists.

### 1.1 Goals

1. **The same loci in every sample**, arrived at independently.
2. **Represent what the estimates need**, which is not the same as representing the genome — §3 is
   where those two come apart, and it is the substantive decision in this document.
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
- **The threshold sets the count**, from the analysed length: `threshold / hash_range = n /
  analysed_length` for a target of `n` positions. The realised count is binomial around `n`.
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

### 3.2 Decision: the `cap` lowest hashes per stratum, over a reference catalog

**Give every STR locus a value `hash(contig, start, seed)` and keep, within each stratum, the `cap`
loci whose value is smallest.** Everything is kept where a stratum holds fewer than `cap`.

Four properties, and the last two are what the even spread lacked:

- **It is a uniform random sample of the stratum.** The hash is uncorrelated with position, so the
  chromosome-start effect that made a prefix 32.6% low cannot reach it.
- **It keeps exactly `cap`.**
- **It does not depend on the order loci are seen in.** The kept set is a function of the *set* of
  loci and the seed, nothing else — so a region-sharded enumeration and a single-threaded one agree,
  and merging two shards is taking the lowest `cap` of the union. The even spread has no such
  property: it counts segments as they arrive, so a different traversal keeps different loci.
- **Its working state is `cap` values per stratum**, a bounded heap, so enumerating a genome costs
  nothing to hold.

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

### 3.3 The pass is the catalog, and it should be computed once

**What the rule needs is every STR locus in the analysed regions with its stratum** — the (motif
period, reference repeat count) pair — and nothing else. That is a **catalog**: a list of the loci a
reference contains, derived from the reference alone.

**That catalog is specified separately and built as independent work**:
[`repeat_catalog.md`](repeat_catalog.md). It is produced inside the pass that already streams the
whole FASTA, and it records **detections rather than loci** — period, span, score, motif and purity —
so the loci at any copy floor, and their strata, are derived from it without a re-scan. That is what
makes the count this rule needs available at whatever floor the run chooses.

**Decision: compute the catalog once per (reference, scan settings) and select from it.**
Three things it buys, in increasing weight:

- **The reference walk stops being per sample.** Region typing scans the whole reference with the
  tandem-repeat detector; a fifty-sample cohort does that fifty times today, for an answer that cannot
  differ between samples.
- **The per-stratum totals fall out of it**, and §3.5 needs them regardless of which selection rule is
  used — so the pass is not chargeable to this decision at all.
- **The identity check becomes a digest of the catalog rather than a list of parameters**, which is
  strictly stronger: it also catches a change in the detector that the parameters do not express.

**Not production's catalog.** `src/ssr/catalog/` builds one via the external `trf-mod` binary and is
frozen production; ng does not depend on it, and
[`typed_regions.md`](typed_regions.md) is explicit that it is a comparison oracle and never an ng
dependency.

**One thing to carry across from that spec, because it decides whether this rule can run at all**: the
catalog is built at settings **at least as permissive as any reader asks for**, and a run whose copy
floor is *lower* than the catalog's is refused rather than served a short list
([`repeat_catalog.md`](repeat_catalog.md) §4). A short list here is a wrong per-stratum total, and
nothing downstream would notice.

**What it costs is a stale-artefact failure mode that not having one cannot produce**, and the
mitigation has to be load-bearing rather than advisory: because the kept set is a pure function of the
catalog's contents, a catalog built from a different reference or different detector settings selects
a *different* set, silently. **So the catalog's digest travels with every sample's records and the fit
refuses to pool across a mismatch** (§5.1).

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
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §3) rather than a note
  here.

**The true counts now come from the catalog rather than from a pre-pass**, which is one of the three
things §3.3 lists it buying.

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

[`parameter_prepass.md`](parameter_prepass.md) §4.1 works out where two million generic positions
land — six million read observations at three reads a plant, pinning the error rate to about **one
part in eighty**, and about 2,000 heterozygous sites, pinning heterozygosity to about **3%**. Both are
far inside what a genotype likelihood can feel.

**So the experiment is to sweep downward, not upward.** Refit at 2 M, 500 k, 100 k and 20 k positions
on the same drawn data and find **the first budget at which the fitted values move by more than a
caller could feel**. That number is the real default; two million was chosen for the cross-sample
statistics ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5) and has never
been checked against this use. §6 carries it as the open question.

### 4.2 One thing the budget does not buy

The joint fit's distinctive advantage — telling a mismapped locus from a heterozygous one, because a
mismapped locus is noisy in *every* sample
([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §2.2) — comes from having **many
samples at one locus**, not many loci. It does not improve as the budget grows and does not degrade as
it shrinks. **The budget buys precision for the pooled rates and nothing at all for that.**

### 4.3 The STR cap may turn out to be no cap at all

[`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.1 suspects the whole STR
set may fit with no sampling: there are orders of magnitude fewer STR loci than genome positions, and
a far larger fraction of them vary. **Implement the cap anyway** — setting it high enough to admit
every locus is free, and discovering later that a bigger genome needs one and having no mechanism is
not.

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
seed and the cap. Selecting the STR set holds `cap` values per stratum, a bounded heap. What travels
with a sample's data is the *inputs* (§5.1) — a seed, three digests and two integers.

**Compute.** The catalog is one walk of the reference with the tandem-repeat detector, run **once per
(reference, region-typing parameters)** rather than once per sample. Selecting from it is a hash and a
heap push per locus.

**Errors.** A run whose catalog differs from another's has different loci in it, so the two samples did
not select the same set. That is not a degraded estimate but a meaningless one, and the fit must be
able to refuse rather than average. **A stale catalog is the way this happens in practice**, which is
why the check is a digest and not a promise (§3.3).

### 5.1 What must travel with a sample so the fit can check

Six values, checked for agreement across every sample before anything is pooled:

| value | why a mismatch is silent |
|---|---|
| the selection **seed** | a different seed selects a disjoint set; both look well-formed |
| the **reference** digest | coordinates mean different things, so "the same position" is not |
| the **analysed region set** digest | the likeliest to differ by accident, because a BED feels like a runtime convenience. It is not: it defines what population of sites the estimate describes |
| the **catalog** digest | it decides which loci are STR loci and what stratum each is in, so it sets both the boundary the two sets partition on and the pools the STR rule draws from. **A stale catalog is the realistic failure** (§3.3) |
| the **generic target count** | a different target is a different set, not a subset of the larger one |
| the **STR per-stratum cap** | a smaller cap *is* a subset of a larger one under §3.2's rule — but of a different size, so two samples still hold different loci |

**The catalog digest replaces "the region-typing parameters" that
[`parameter_prepass_cohort.md`](parameter_prepass_cohort.md) §2 lists, and is strictly stronger**: it
also catches a change in the detector itself, which the parameter list cannot express. The first three
are already required there; the last three are this document's addition and fail in exactly the same
silent way.

---

## 6. Open questions

1. **How many STR strata are there, and how many loci does each hold?** — OPEN, and it is the first
   thing to measure, because three decisions wait on it: what the cap should be (total loci is
   `cap × non-empty strata`), whether any sampling is needed at all (§4.3), and how much of the joint
   fit's comparison against the per-sample route is even a comparison — if every locus is kept, the
   joint route holds the same STR loci *and* remembers which was which.
   **Measurable today**, by running `type-regions` over the tomato reference and tabulating its rows:
   no reads, and it is the same artefact §3.3 makes an input.

   **If the answer is "they all fit", most of §3 dissolves.** With no cap there is no sampling, no
   selection rule to get right, and no reweighting to remember — the strata are represented in
   proportion because everything is kept. §3 is what happens when they do not fit, and building it is
   still right ([§4.3](#43-the-str-cap-may-turn-out-to-be-no-cap-at-all)), but the measurement decides
   whether it ever fires.
2. **Where does the generic budget start to matter?** — OPEN. §4.1 gives the experiment.
   *Leaning:* well below two million, since the error rate is already pinned to one part in eighty
   there.
3. **What fraction of STR loci vary across a cohort?** — OPEN, and it is the other half of question 1:
   it decides whether a cap set to keep everything actually delivers the ~10,000 varying loci the
   cross-sample statistics want ([`parameter_prepass_census_sites.md`](parameter_prepass_census_sites.md) §5.1).

---

## 7. How we know it works

1. **Every sample selects the same loci.** Run the selection over samples with different coverage,
   different read lengths, and a region-sharded walk at several thread counts; the set must be
   identical every time, and identical to the set computed directly from the six values of §5.1.
   **Include a `--regions` case**: the chosen positions must all lie inside the BED, and the realised
   count must hit the target computed from the region set's length rather than the genome's — the
   arithmetic most likely to be written against the contig table by reflex.
2. **The generic selection is unbiased.** On synthetic data with a known frequency spectrum, the
   spectrum over the chosen positions must match the one over every position, within its own error.
   If the selection ever came to depend on the data, this is where it shows.
3. **The stratified selection is random and exactly capped, not a prefix.** Every stratum with more
   loci than the cap keeps exactly `cap`. **Check the distribution, not only the count**, since a
   prefix passes a count check and is the failure measured at 32.6% low (§3.1): the kept loci's
   positions must be spread across the region set the way a uniform draw would be, and a fitted
   slippage rate over them must match the rate over the whole stratum.
4. **The selection does not depend on the order loci were enumerated in.** Shuffle the catalog, or
   shard the enumeration several ways at several thread counts; the kept set must be identical every
   time. **This is the property the even spread does not have** (§3.2), so it is the test that says
   which rule was implemented.
5. **A larger cap contains a smaller one.** Selecting at `cap` and at `2·cap` must give nested sets —
   a property of the lowest-hash rule that makes the downward budget sweep of §4.1 a sequence of
   subsets rather than a sequence of unrelated draws.
4. **The reweighting undoes the stratification.** STR diversity computed from the kept loci with the
   stored per-stratum weights must match the diversity computed from **every** STR locus, and must
   differ from the unweighted version. The second half is what makes the first half matter, and the
   whole thing is silent when the reweighting is omitted (§3.3).
5. **The two sets partition, and neither leaks.** No locus appears in both. Run it with region
   typing's copy floors moved and confirm both sets change together — the partition follows the
   delimitation rather than a second copy of its rules.
6. **A mismatch is refused.** Two samples selected under different seeds, references, region sets,
   region-typing parameters, target counts **or caps** must produce an error and not an average (§5.1).
