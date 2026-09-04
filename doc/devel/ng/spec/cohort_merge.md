# ng — the cohort merge: per-sample observations into cohort observations, in parallel

*Status: design spec draft, 2026-08-17. No code yet — this settles the design, and it
incorporates the owner's 2026-08-17 ruling: a cohort locus wider than `max_cohort_locus_span`
cannot be built — it is counted as a failed locus, never built, never sent downstream (§3.2).
This is the home [`run_streaming.md`](run_streaming.md) §10 promised for "the cohort merge's
reconciliation", and it also discharges the confirmation §10 demanded — that the merge frontier
stays bounded while spans are reconciled (§4.5 here). **That document's §1.2 and §10 entries were
repointed here on 2026-08-28**, having outlived this spec by eleven days. Companion architecture doc:
[`../arch/cohort_merge.md`](../arch/cohort_merge.md) (the types and interfaces). Reads on
[`locus_generation.md`](locus_generation.md) (what a per-sample observation is),
[`typed_regions.md`](typed_regions.md) (what a segment is), and
[`run_streaming.md`](run_streaming.md) (the run this component sits inside — especially its §3.2
streaming merge, §3.3 two-phase decode, §4.3 segment independence, §12 oracles).*

---

## 1. What this is

**The cohort merge turns k per-sample streams of locus observations, in coordinate order, into
one stream of cohort observations, in coordinate order, ready to be called — and it does that
work on several threads at once.** It is the stage that was the wall floor in production, where
the equivalent code — the cohort fold in
[`cohort_integration.rs`](../../../../src/var_calling/cohort_integration.rs) — advanced every
sample behind a single global frontier and therefore ran on one thread, however many the run had
(its module doc, [`:1-38`](../../../../src/var_calling/cohort_integration.rs), describes that
loop; §10 here carries what it teaches).

[`run_streaming.md`](run_streaming.md) already fixed the outer shape: the merge streams, keys on
coordinates, and its resident set is a frontier (§3.2 there); sources answer in two phases, a
cheap position summary for every position and a full build only for the positions the cohort kept
(§3.3 there). This document fixes what was deferred: **what a cohort observation is** — how
per-sample observations whose spans differ become one (§4) — and **how the work is divided** so
that builders can run concurrently along the genome without coordinating (§5, §6).

The design in one paragraph: the genome is dealt out in short adjacent **regions**, about twice
`max_cohort_locus_span` wide, one to each **builder**. A builder takes what its sources give it
over that region, groups positions into cohort loci by shared reference bases, discards the loci
no sample varied at,
**fails** the loci wider than `max_cohort_locus_span` — counted, never emitted (§3.2) — and fully
builds only
the survivors. It starts a locus only inside its own region but follows one past the end if a
deletion carries it there, so builders overlap. Everything goes to a single **organiser**, which
resolves those overlaps by one rule — the locus whose first position is earlier stands — and
releases loci in genome order. No builder talks to another.

**The first build is the direct path, and it has no position summaries.**
[`run_streaming.md`](run_streaming.md) §2 settles that direct mode is built first, and there the
walk has already minted every observation before this component sees it. So a builder groups and
filters the observations themselves, and assembling a survivor costs nothing beyond moving what it
already holds — nothing is deferred, because there is nothing to defer it from.

That split is what the psp path adds later, and only to avoid decompressing evidence it is about to
discard ([`run_streaming.md`](run_streaming.md) §3.3). **Everything else in this document is the
same in both paths** — the grouping, the two verdicts, the ownership rule, the organiser, the
determinism argument. Where the text below says a builder decides from a cheap per-position summary
and builds the survivors afterwards, read that as one step in the direct path and two in the psp
path; §4.4 says so at the point it matters.

### 1.1 Goals

1. **Builders share nothing but their results.** A builder reads observations from a cache it
   cannot modify (§6.4) and hands the loci it finished to the organiser. That hand-off is the only
   thing it does that another thread can see: no builder writes anything another builder reads, and there is no counter or cursor
   recording how far the cohort as a whole has got. Such a counter is the obvious way to write
   this — every builder advancing one shared position — and it is why the merge was serial in
   production: a position every builder must agree on is a lock every builder queues at (§6).
2. **The output is a function of the data.** The same cohort observations — and the same failed
   loci — in the same order, at any builder count, any look-ahead, any division into regions —
   [`run_streaming.md`](run_streaming.md) §12.2 requires the VCF byte-identical at any worker
   count, and this component is where that property is easiest to lose (§9).
3. **Bounded memory across the committed range** — one sample to several thousand, a few reads a
   position to several hundred (`CLAUDE.md`, *What this caller has to work on*). Every bound in
   §8 is a formula in the sample count and the look-ahead.
4. **The expensive work happens only for the loci that survive.** For each locus a builder first
   answers two cheap questions, using the position summaries alone: how wide is this locus, and did
   any sample record a non-reference read inside it? A locus that is wider than
   `max_cohort_locus_span`, or that
   no sample varied at, is dropped there and then — the builder never goes on to gather the
   samples' evidence for it (§3.2, §4.3). Gathering is the expensive step, and in the psp path it
   is a decompression as well, so what is dropped is never decompressed
   ([`run_streaming.md`](run_streaming.md) §3.3).

### 1.2 Non-goals, and what this document does not do

- **It does not decide which alleles are worth calling.** It unifies the alleles the samples
  actually showed into one table per locus (§4.2) — that is the point of a cohort observation. What
  it does not do is choose candidates from that table, normalise their representation for output,
  or weigh them: those are the calling steps' (§13).
- **It does not call.** No likelihood, no genotype, no QUAL. The output is evidence.
- **It does not define the psp encoding.** It adds two requirements to the encoding spec's
  inheritance — the position summary must carry the reference span, and the header must record the
  observation reach ceiling (§13) — and nothing else.
- **It does not schedule threads.** The pool, the look-ahead accounting and the yield point
  belong to the run objects ([`run_streaming.md`](run_streaming.md) §3.5); this component defines
  the work unit, the per-unit function, and the ordering invariant those objects enforce.
- **It does not decide where the calling function runs** relative to the reorder buffer — open
  question 5 (§14), with a leaning.

### 1.3 Vocabulary

- **observation** — one sample's evidence over one stretch of genome: a
  `SampleLocusObservations` ([`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)),
  with its `region` (1-based, inclusive — `GenomeRegion`,
  [`types.rs:79`](../../../../src/ng/types.rs)). Same use as [`run_streaming.md`](run_streaming.md) §1.3.
- **cohort observation** — every sample's observations over one cohort locus, grouped and
  collated by this component; the unit the caller consumes. §4 defines it.
- **cohort locus** — a group of positions chained by shared reference bases across the cohort
  (§4.1). A cohort locus either becomes a cohort observation, is dropped as non-variable (§4.3),
  or **fails** (below).
- **failed locus** — a cohort locus wider than `max_cohort_locus_span`: not built, nothing sent
  downstream, counted in the run summary (§3.2, §3.3). The owner's 2026-08-17 ruling.
- **span** — the count of reference bases an observation's `region` covers. An insertion's span
  is 1 (its anchor base); a deletion's span is the deleted run plus its anchor.
- **reach** — the last reference position an observation covers: `start + span − 1`. The word and
  the arithmetic are production's
  ([`cohort_integration.rs:46-48`](../../../../src/var_calling/cohort_integration.rs)).
- **`max_cohort_locus_span` / `max_cohort_locus_span`** — the *policy* constant: the widest cohort
  locus the
  caller undertakes to build. Default 50 bases. §3 is its section
  derivations only.
- **the reach ceiling** — the *physical* constant: the widest reference span any minted
  observation can have — the generic generator's `max_record_span`
  ([`pileup/generator.rs:93`](../../../../src/ng/locus_generation/pileup/generator.rs)). It bounds
  how far the observation cache must reach and nothing else (§5, §6.4).
- **position summary** — the cheap facts a builder needs about a position before it decides
  anything: what a record there covers, and whether any sample recorded something other than the
  reference. **In the direct path these are read straight off the observations the walk minted** —
  there is no separate object. In the psp path they are the file's light columns, decoded ahead of
  the heavy ones (production's `TwoPhaseSegment`,
  [`sample_reader.rs:698`](../../../../src/var_calling/sample_reader.rs)). The name is the role,
  not a data structure.
- **candidate position** — a position whose position summary, in at least one sample, records a
  non-reference observation count above zero. Candidates are what the variability filter keeps
  around (§4.3).
- **building region** — the stretch of genome assigned to one builder,
  `cohort_locus_builder_regions_len` bases, 20 by default. A builder starts a locus only inside its
  own building region and may finish one outside it (§6.1). Written just *region* below, where
  nothing else in this document is called one.
- **segment** — a different thing entirely, and never called a region here: the reference's own
  division into an STR tract, a bundle, a satellite, or the generic stretch between them, made by
  the typed-region generator before any read is examined
  ([`run_streaming.md`](run_streaming.md) §4.2). A building region is work; a segment is
  territory. Many building regions fit in one segment.
- **builder** — the worker that merges and builds the loci starting in one region.
- **organiser** — the single thread that collects the builders' output, resolves the overlaps
  between neighbouring regions, and releases loci in genome order (§6.3).

---

## 2. Where it sits

Both variant callers — `AlignedFilesVariantCaller`, which reads the alignment files, and
`PspVariantCaller`, which reads the psp files — drive `call_vars_in_segment(region, sources)` over a
loop of regions
([`run_streaming.md`](run_streaming.md) §3.1, §3.5). This component is the inside of that
arrangement, made precise:

- the **regions the loop deals out** are the builders' regions — about twice `max_cohort_locus_span`
  wide, with
  subdivision of long generic segments specced as the extension (§6.1). Subdivision refines
  [`run_streaming.md`](run_streaming.md) §4.3/§4.4 of the
  never-cut rule rather than a violation of it;
- the **merge inside one region is this component's builder** (§6.2): fold position summaries, group,
  filter, fail the over-wide, build survivors;
- the **in-order yield** the run objects already own ([`run_streaming.md`](run_streaming.md)
  §3.5) gains one invariant from this component: every region delivers exactly one result — its
  observations *and its failed-locus count* — even when both are empty (§6.3).

The walk stage (`SampleObservationGatherer`) does not use this component at all: it has one
sample and writes everything; there is no cohort to merge and no variability to filter.

---

## 3. `max_cohort_locus_span`

### 3.1 What it is and what it buys

**`max_cohort_locus_span` is the widest cohort observation this component will build, measured in
reference bases. A command-line parameter, default 50.** The default is the owner's number,
unmeasured — soft; open question 3 (§14) names the measurement. A run over long reads is expected
to set it higher, since the widest event worth merging into one locus grows with the reads.


What it buys:

- **a counted refusal instead of a silently wide locus.** A locus wider than `max_cohort_locus_span` is
  reported as a failure (§3.2) rather than built.
- **a hard ceiling on what this component emits, for generic loci.** A generic cohort observation
  is at most `max_cohort_locus_span` bases wide, so the assembly window (§4.5) and §8's memory
  prices are bounded by a
  number the operator chose rather than by the data. **An STR locus is not bounded by
  `max_cohort_locus_span`** (below),
  so the true ceiling on an emitted observation is the larger of `max_cohort_locus_span` and the
  widest STR tract the
  segmentation admits — a tract longer than `max_str_len` is a satellite and yields no locus at all
  ([`region_typing/mod.rs:223,235`](../../../../src/ng/region_typing/mod.rs), default 100 bases).
  At the two defaults that is 100 rather than 50, and both are the operator's to set.

**It does not make a cut safe.** Deciding where one builder's ground can end without splitting a
locus is a different question, and §5 answers it: none is looked for.

`max_cohort_locus_span` governs **generic** loci only. An STR observation's span is its reference
tract, which
the segmentation defines and which may exceed 50 bases; its segment's boundaries are cuts with no
check, and its width is a fact about the reference rather than a claim about the reads. The
same holds for bundles.
Satellites produce no observations at all ([`locus_generation.md`](locus_generation.md) §5).

An insertion never widens a locus past `max_cohort_locus_span` at any length: its reference span is
its anchor
base. `max_cohort_locus_span` constrains what a locus *covers on the reference*, not how many bases
a read
showed.

### 3.2 A locus wider than `max_cohort_locus_span` is a failed locus — the owner's ruling

> Owner, 2026-08-17: *"Then that locus can't be built, it has to be counted as a failed locus,
> not built and not sent downstream."*

**If a cohort locus comes out wider than 50 bases, it is not emitted.** Nothing goes downstream
over the ground it covers, and the run counts it as a failure.

**In every other respect it is an ordinary locus.** It is grouped like one, it owns its ground like
one, and it displaces the overlapping loci of neighbouring builders like one (§6.1). Only the two
things above are different. Treating it as a locus that happens not to be emitted, rather than as a
locus that does not exist, is what keeps the overlap rules free of special cases — and what stops
the ground it covers being emitted twice over by builders that never saw what opened there.

Two different things make a locus that wide, and the rule does not distinguish them. One read may
carry a deletion longer than 50 bases — ng mints observations up to `max_record_span`, 5,000 bases
by default ([`pileup/generator.rs:93,141`](../../../../src/ng/locus_generation/pileup/generator.rs)).
Or several short events in different samples may overlap each other in a chain, none of them wide
on its own, until the stretch they jointly cover passes 50. Either way the caller would have to
build one locus wider than it undertakes to call, so either way it declines. What matters is how
wide the locus ended up, not how it got there.

**The failure is per locus, not per sample — and that is the price of `max_cohort_locus_span`.** One
sample's
200-base deletion fails the locus for the whole cohort: a sample that had an ordinary SNP at a
position inside that deletion's span is chained into the same locus by the shared bases, and its
SNP is suppressed with it. `max_cohort_locus_span` buys a bounded, honest caller at the cost of
every bystander
variant inside a failed locus's span. That price is why the count must surface (§3.3) — it is
the only way an operator can see `max_cohort_locus_span` charging more than expected.

**Failure is decided on the fold, so it is cheap and deterministic.** Grouping and the span check
read only the folded position summaries (§4.1) — a failed locus is never assembled,
so a failure costs what a dropped non-variable locus costs. And because the fold is a pure
function of the per-sample streams (§9), the *set* of failed loci is too: the same input yields
the same failed loci at any builder count and any division into regions. A builder does not
"discover" an over-long event mid-build — the fold is complete for the region before any locus
is built — and it resumes at the next group, which the closure guarantees is disjoint from the
failed one (groups never share a base, §4.1). Nothing about failure involves timing.

**An observation wider than `max_cohort_locus_span` needs no special handling at all.** It makes a
locus wider than the bound, that locus fails, it is counted, and its span suppresses the loci that
overlap it (§6.1) — the ordinary path. Nothing in this component depends on observations staying
below any width (§5), so however wide one is, it is data rather than a defect.

### 3.3 What the caller sees, the count, and what a non-zero count means

Downstream of the merge, a failed locus is invisible in the stream: no cohort observation, no
variant, no VCF record over its ground — for any sample. In the VCF alone, that absence is
indistinguishable from "analysed and found nothing", which is exactly why the ruling says
*counted*, not merely *dropped*:

- **Where the count surfaces.** Each builder returns its region's failed-locus count with the
  region's observations (§6.3); the run object sums them at the end of iteration — the same
  finish-time tally path [`run_streaming.md`](run_streaming.md) §8 mandates for read-filter
  tallies, and with the same failure mode if it is skipped: a per-worker count that never gets
  summed under-reports by the worker count, silently. The total (with a per-contig breakdown left
  to the implementation) is part of the run summary every calling run reports. How failed loci
  surface beyond the count — spans in a log, a sidecar file — is open question 4 (§14).
- **What a reader does with a non-zero value.** The count is a calibration signal: it says the
  bound is smaller than what the data keeps presenting. The obvious case is long-read data under
  the 50-base default, where genuine mid-size deletions fail loci wholesale. The reader's move is
  to inspect the failed spans (open question 5's surface), and if they are real signal — lengths
  clustering just above `max_cohort_locus_span`, or a long-read run — raise `max_cohort_locus_span`
  and call
  again. Nothing about the observations changes (§3.1), so the second run starts from the same
  evidence.

**What is lost, stated flatly: every locus wider than `max_cohort_locus_span`, for the whole cohort
— the wide
event and every bystander variant chained into it.** At the default that means no deletion longer
than 50 bases is ever called, at any depth, any cohort size. If those events are ever wanted, the
home is a separate pass over the emitted records — the same home
[`run_streaming.md`](run_streaming.md) §4.3 gives cross-segment events — never a coupling between
in-flight regions (§13).

Alternatives the ruling supersedes, recorded because the earlier draft argued them:

- **Clamp at mint** (the earlier draft's choice: never mint the wide observation; demote its
  reads to `reads_without_observation`). Rejected by the ruling, and for cause: it silently calls
  the *rest* of the locus from a picture with the strongest evidence removed, it moves a calling
  policy into the making of observations so that changing `max_cohort_locus_span` means making them
  again, and
  it hides the event from the count.
- **Split the event** into sub-bound pieces — mints alleles no read showed.
- **Cap groups greedily and emit overlapping loci** — builds ground the caller does not
  undertake to determine, merely in smaller pieces; the ruling's refusal is the honest form.
- **Refuse the run** — one read event ending a several-thousand-sample run.

---

## 4. What a cohort observation is

### 4.1 Grouping: shared bases chain observations into cohort loci

The builder folds the k samples' evidence over its region into one cohort summary: the union of
positions, and per position the widest reference span any sample's observation covers and whether
any sample saw a non-reference read. This is production's `CohortSpanFold` — union on positions,
`max` on the aggregates, commutative and associative so the fold order cannot matter
([`cohort_integration.rs:64-149`](../../../../src/var_calling/cohort_integration.rs), proven
order-independent by its tests, [`:1302,:1960`](../../../../src/var_calling/cohort_integration.rs)).

**Only one of those two comes off the observation as a field.** The span is
`SampleLocusObservations::region`
([`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)). There is **no**
non-reference count: an observation carries its bases and its read count, and whether those bases
are the reference's is decided by comparing them — `*observation.bases == *locus.reference_bases`,
which is exactly what the census writer does today
([`census.rs:2084`](../../../../src/ng/parameter_estimation/joint/census.rs)), summing `num_obs`
over the observations that differ. So the fold derives that half rather than reading it, and a
builder pays one sequence comparison per observation to do so. It is the same derivation
production makes: its non-reference column is derived from the light columns rather than stored
([`sample_reader.rs:715-719`](../../../../src/var_calling/sample_reader.rs)).

*(In the psp path the same number may instead arrive precomputed, since a file is free to store
what a walk would have to recompute — that is the encoding spec's call, §13.)*

The builder then walks the region's positions from left to right, gathering them into loci. It
keeps one locus open at a time, and that locus covers the reference from its first position to the
furthest base any observation in it reaches.

**The next position joins the open locus if it falls within that reach**, because some observation
already in the locus covers that base — the two overlap, so they cannot be called apart. If it
falls beyond the reach, the open locus closes and the next position starts a new one. This is
production's rule, unchanged (`derive_is_kept`,
[`cohort_integration.rs:166-187`](../../../../src/var_calling/cohort_integration.rs)).

So a SNP at position 10 and a SNP at 11 are two loci, not one: neither observation covers the
other's base. A deletion at 10 spanning to 14 and a SNP at 12 are one locus, because the deletion
covers 12.

**Nothing caps how far this chains.** A locus grows as long as the data keeps overlapping, and each
closed locus is disjoint from its neighbours because consecutive ones never share a base. Whether
what came out is too wide is a separate question, asked next.

Then `max_cohort_locus_span` is applied per locus (§3.2): a locus at most that wide goes on to the
variability filter; a wider one is failed, counted, and done with.

**A locus may run past the end of the region its builder was given, and often does** — that is the
whole of §6.1's ownership rule, and it is what a deletion reaching forward makes necessary. What a
locus never crosses is a **segment** boundary: no observation crosses one
([`run_streaming.md`](run_streaming.md) §4.3), so no chain of overlapping observations can either.
Segments are the reference's own division; a builder's region is a slice of work. The first bounds
a locus, the second does not.

Because the chaining rule is a property of the observations and not of where the work was cut, the
set of loci — and the set of failures among them — is the same however the genome was divided
between builders (§9).

### 4.2 Collation: one allele table, and every sample's support for it

The observations that overlap the locus span are its **members** — and since loci are disjoint
(§4.1), every observation is a member of exactly one locus.

A sample whose reads all matched the reference over these positions still has members: its
reference observations are collected like everyone else's, so its depth is present when genotypes
are weighted. A sample with *no coverage* over the span has no members at all, which is a different
fact and stays one.

**The samples' observations are projected onto the locus span and unified into one allele table,
and this is the builder's work.** It is what makes a cohort observation a cohort's observation
rather than a bundle of per-sample pieces, and it is the reason the locus was grouped in the first
place.

Two steps, both mechanical:

- **Projection.** A sample's observation may be narrower than the locus — one sample has a SNP at
  a position another sample's deletion also covers — so each observation's sequence is widened to
  the full locus span, padded with the reference bases either side. Those bases travel on the
  observation already (`reference_bases`,
  [`locus_generation/mod.rs:46`](../../../../src/ng/locus_generation/mod.rs), whose doc calls them
  "what a wider-span projection needs when samples merge"), so nothing has to be fetched.
- **Unification.** Two projected sequences that are identical are the same allele, wherever they
  came from. The locus ends with one table of distinct alleles, the reference among them, and each
  sample's support expressed against that table — its per-read moments summed where two of its own
  observations projected onto the same allele **and came from the same read group**. The second
  half was added 2026-08-23 ([`calling_prerequisites.md`](../impl_plan/calling_prerequisites.md)
  B1): two reads showing the same bases from two lanes have different error rates, so a read
  likelihood may not fold them into one term ([`read_likelihoods.md`](read_likelihoods.md) §2.3).
  A sample with one read group is unaffected, which is most of them.

**Unification by exact match is only sound because indels were left-aligned upstream.** The same
deletion written at two placements would otherwise project to two sequences and become two alleles,
silently splitting one variant's evidence in half. Read preparation canonicalises that before the
generator mints anything (`LeftAlignPreparer`,
[`read/left_align.rs:92`](../../../../src/ng/read/left_align.rs)), so equal sequences here really
do mean the same allele. It is worth stating because the failure is invisible in the output: two
half-supported alleles look like a noisy site, not like a bug.

**A sample's own evidence is not merged across alleles, nor across read groups.** Its members'
counts and moments stay attached to the allele each projected onto and to the group the reads came
from, so the caller sees per-allele support per read group per sample — the first because a
genotype likelihood needs the alleles apart, the second because a read likelihood needs the groups
apart.

**Production does this in a merger too, and it is worth being exact about which one**, because it
has two and they do different jobs. `PerPositionMerger`
([`per_position_merger.rs:145`](../../../../src/var_calling/per_position_merger.rs)) is the k-way
merge over the per-sample files: it lines the streams up by position and says which samples had a
record there, and it unifies nothing. The allele work is `PerGroupMerger`
([`per_group_merger.rs:585`](../../../../src/var_calling/per_group_merger.rs)), which runs after
grouping and emits *"a unified allele set"* per group, having *"projected every sample-local allele
onto the group's reference span"* and *"deduplicated by byte equality"*
([`:1-20`](../../../../src/var_calling/per_group_merger.rs)). So unification is merging work done
after grouping — which is exactly where this component's builder does it. The one thing production
folds in there that ng keeps downstream is the likelihood: its Stage 5 emits both, and here the
builder emits evidence only (§1.2).

### 4.3 The variability filter: builders discard, nothing upstream does

**Most of the genome does not vary, and the builder throws those loci away itself.** A locus no
sample showed enough non-reference evidence at — which includes every locus where no sample varied
at all — is dropped by the builder: nothing is assembled, nothing goes to the organiser. Only what
is left waits to be emitted, which is why the organiser holds so little.

**The rule asks one sample at a time.** A locus is built when **some single sample** shows at
least

    max(min_alt_obs, ceil(min_alt_read_share × that sample's compared reads))

non-reference reads. *Compared reads* is the denominator the numerator was drawn from: the reads
whose whole sequence over the locus was compared against the reference
([`SampleLocusObservations::reads_compared_with_reference`](../../../../src/ng/locus_generation/mod.rs)),
which is neither read depth nor the reads that merely covered the ground. A partial read, whose
bases stop inside the locus, is in neither half.

The two parameters are the two ends of the depth axis, and neither alone would serve:

- `min_alt_obs`, default **2**, decides at low coverage, where a share of three reads rounds to
  nothing;
- `min_alt_read_share`, default **2 in 100**, decides at high coverage, where two reads out of
  three hundred is the sequencing error rate rather than an allele.

**Why one sample and not the cohort's sum.** ng summed every sample's non-reference reads and
compared that against `min_alt_obs` until 2026-08-19. The owner retired that rule on the
measurement below: because the bar was a fixed count and the sum grows with the cohort, the rule
stopped filtering as samples were added. On one 100 kb interval of the tomato panel it discarded
997 loci in every 1,000 at one accession, 878 at 16 and only 439 at 63 — where the merge went on
to assemble a cohort observation at 546 positions per kilobase, better than one every two bases,
and assembling them cost 170 ms of a 425 ms single-threaded merge.

Taking the share of the **cohort's** reads instead fixes the count and breaks something worse. That
bar grows with sample count while a rare allele's evidence does not: one carrier's handful of reads
is the same handful whether the run holds 10 samples or 10,000. At 63 tomato accessions it would
ask for about 14 non-reference reads, so an allele carried by one or two accessions at three reads
a position could not be seen; at 1,000 samples it would ask for about 200, near 100 carriers, and
nothing below about 5% frequency would survive. A study is run at that size to find exactly the
variation such a rule discards. The per-sample denominator has no such drift — it is the same
question at one sample and at a thousand.

**What each end of the range gets.** At **one sample** the rule is exactly the cohort-sum rule it
replaced, because a sum over one sample is that sample's own count. At **3 reads a position** the
share never reaches the floor and the floor is the whole rule — on the tomato panel it moved 12,034
kept loci to 12,033. At **313 reads a sample**, which is the GIAB trio's depth, the share does all
the work: 53,796 kept loci become 2,721 over 136 kb, and a true heterozygote there shows about 156
non-reference reads, a factor of 20 above the bar.

**What it costs, stated plainly:** a real variant no single sample showed twice is dropped, and no
downstream step can recover it, because nothing is emitted for that locus at all. At three reads a
position a heterozygous carrier shows a single read about four times in ten, so this is not an
empty class. It is the same trade production makes, whose rule is also per-sample (`derive_is_kept`
over `max_nonref_obs`,
[`cohort_integration.rs:166-187`](../../../../src/var_calling/cohort_integration.rs)), at the same
default of 2 ([`var_calling/mod.rs:72`](../../../../src/var_calling/mod.rs)). Lowering
`min_alt_obs` to 1 restores those variants and restores the cost.

**Neither number has been measured against a truth set.** Both were chosen against loci-kept counts
on two benchmarks a hundredfold apart in depth. What is still owed is a recall check: the tomato
panel's and GIAB's calls under the old rule and the new, so the variants this drops can be counted
rather than reasoned about.

**A kept locus is kept whole**, quiet positions and all — the threshold decides whether the locus
survives, never which of its positions do
([`cohort_integration.rs:1344`](../../../../src/var_calling/cohort_integration.rs)).

### 4.4 Assembling the survivors

For each locus that passed both verdicts, the builder collects every sample's observations at the
member positions into one cohort observation (§4.2).

**In the direct path this costs nothing extra**, and it is the path built first: the walk minted
those observations and is still holding them, so assembling a survivor moves what is already in
hand. There is no second read of anything.

**In the psp path it is where the heavy columns are decompressed** — the deferred inflation of
[`run_streaming.md`](run_streaming.md) §3.3, production's `set_variable_rows` shape
([`sample_reader.rs:789`](../../../../src/var_calling/sample_reader.rs), one column at a time, a
transient of about one column). This is the only step in the whole document whose *cost* differs
between the paths, which is why the discards above matter more there: what is discarded before this
point is never decompressed at all.

### 4.5 What is held while building

A builder works on one locus at a time and holds that locus's observations and no others. Since a
built locus is at most `max_cohort_locus_span` bases wide — wider ones fail before anything is
assembled (§3.2) — what one builder holds is one window of that width per sample, whatever its
region's length and whatever the region contains.

This is the confirmation [`run_streaming.md`](run_streaming.md) §10 asked for: reconciling the
samples cannot grow what is resident beyond one bounded window per sample, and
`max_cohort_locus_span` is what bounds it.

**The observations themselves are not the builder's to hold, though — they live in a cache the
organiser owns (§6.4), and that cache is this module's real memory.**

---

## 5. There are no safe places to cut, and none are looked for

**A builder is given a start position and works as though it were safe.** It may not be — a locus
may already be open there, opened by an observation that began earlier — and that is not the
builder's problem. It builds the loci it finds, starting from where it was told to start.

**The organiser is what makes this come out right.** It receives the loci from every builder,
looks for overlaps between them, and resolves each one the same way: **of two overlapping loci, the
one whose first position is earlier stands, and the other is removed** (§6.1, §6.3). A builder that
began inside a locus already open produced loci that overlap that locus; they lose, and they go.

**Why not look for safe cut points instead.** The only places known to be safe before any read is
examined are segment boundaries — no observation crosses one
([`run_streaming.md`](run_streaming.md) §4.3) — and they are far too far apart to divide work at
the 200-base grain this design needs (§6.1). Any finer cut would have to be justified from the data,
which means scanning backwards far enough that nothing could reach across, and that distance is
bounded only by how wide an observation may be: up to `max_record_span`, **5,000 bases by default**
and 65,535 at the ceiling
([`pileup/generator.rs:93,141,45`](../../../../src/ng/locus_generation/pileup/generator.rs);
production's default at [`walker/mod.rs:67`](../../../../src/pileup/walker/mod.rs)). Scanning five
kilobases backwards to hand out two hundred bases of work is not a trade worth making, and the scan
would still be a *check* whose failure the organiser has to catch anyway.

So the design does not check. It builds, overlaps, and resolves.

**What `max_record_span` still does here.** Nothing for correctness — no rule in this document
depends on it. It bounds one thing: how far past a builder's ground the observation cache may have
to reach, because an observation starting inside a region can cover up to that many bases beyond it
(§6.4, §8). It is a memory fact, not a safety one.

**What the segment boundaries still do.** A locus never crosses one, so the organiser will never
find an overlap that spans two segments; and a region boundary is placed at a segment's edge rather
than inside an STR or bundle segment, because ground inside such a segment holds no locus start and
a builder given it would have nothing to do (§6.1). Neither is a safety requirement now — both are
consequences of the segmentation that happen to be convenient.

## 6. The parallel scheme

### 6.1 How the work is handed out — decided 2026-08-17

**A segment cannot be the unit of work.** A segment can be enormous — in the worst case a whole
chromosome — and a builder that owns one must hold every variable locus it finds until the segment
is done, because nothing can be emitted before the ground in front of it is finished. That is the
memory this design exists to avoid, so segments are ruled out as the hand-out unit. *(This reverses
an earlier draft of this section, which made one region one segment. Recorded because the reversal
is the reason the rest of the section looks as it does.)*

**The genome is divided into short adjacent regions, each assigned to one builder —
`cohort_locus_builder_regions_len` bases, 20 by default.** Builders work their regions
concurrently. Because the regions are short, all the builders together cover only a few hundred
bases of genome at any moment, and that span is what the organiser must hold observations for
(§6.4). **That is what bounds this module's memory** — not the segment, not the locus, but how far
apart the builders are.

**A builder starts a locus only inside its own region, and may finish it outside.** Those two
halves are what make the scheme work:

- *Starting inside* is what assigns ownership. A cohort locus belongs to the builder whose region
  contains the locus's first position. That is a property of the data, so two builders can never
  disagree about who owns what, and no message between them is needed to decide it.
- *Finishing outside* is what keeps loci whole. A locus that grows past the end of its owner's
  region — a deletion carrying it forward — is followed to its end by its owner, reading into the
  next builder's ground. The next builder is not consulted and does not stop.

**A region boundary is never placed inside an STR or bundle segment.** Such a segment is one
locus whose span the reference fixes, and it can be wider than a region — up to `max_str_len`, 100
bases by default, against a region of about a hundred. Cutting inside one is not *incorrect*: the
builder owning the tract's start still follows it to its end and still owns it, by the rule above.
But it hands the next builder ground on which no locus can begin, so that builder has nothing to do
and the organiser discards whatever it produced there. So the divider snaps to the segment's end
instead, and a region containing an STR segment is as wide as it needs to be — the width is a
target, not a limit.

**So builders overlap, and the overlaps are cleaned up centrally.** A builder that starts inside
ground already covered by its predecessor's locus is working from a partial view: it never saw the
observations that opened before its region did. The loci it builds there are wrong. It does not
know that, and it is not asked to find out.

Instead each builder sends its surviving loci to a single organising thread, which holds the cache,
keeps the order, and prepares the emission (§6.3). That thread has every builder's output in front
of it, so it can see the overlaps that no individual builder can, and it removes the loci that lost:
**where a locus overlaps ground owned by an earlier locus, the earlier owner's locus stands and the
overlapping ones are dropped.**

**A failed locus takes part in all of this exactly like any other locus** (§3.2). It owns its
ground, it wins against the loci that overlap it, and the builders whose regions fall inside its
span lose their work there just as they would to a successful locus. Nothing about the overlap
rules has a special case for it. The two things that differ are at the end of the line: nothing is
emitted for it, and the run counts it. So the organiser is told a failed locus's span along with
the successful ones, and drops what overlaps it — otherwise the loci built inside its ground from a
partial view would survive into the output with nothing to displace them.

**What the overlaps cost.** A builder wastes the work it did inside a predecessor's locus. Since
roughly one position in a hundred is variable at the measured corner, most joins have no locus in
the overlap at all and cost nothing; the loss is bounded by how far a locus reaches past its own
region, and a locus that reaches further than `max_cohort_locus_span` is failed rather than built.
The exception is the ground of a failed locus, where the reach is bounded not by the 50-base
undertaking but by how wide an observation may be — up to 5,000 bases by default, spanning many regions,
all of whose loci are dropped. That is rare, it is counted, and a run where it is not rare is a run
whose bound is set wrong for its data (§3.3).

### 6.2 Builders: what one task does

One region, one builder, no shared mutable state. The builder:

1. pulls each sample's **position summary** for the region (phase 1);
2. **folds** them into the cohort summary and forms cohort loci by shared-base closure (§4.1);
3. **fails** every locus wider than `max_cohort_locus_span` — counted, not built (§3.2) — and
   **drops** every
   remaining locus with no candidate position (§4.3);
4. **assembles** each survivor: every sample's observations at its positions, collated into a
   `CohortObservation` (§4.2, §4.4);
5. delivers the region's result — its survivors, possibly none, and its failed-locus count,
   possibly zero — to the ordered emitter (§6.3).

**How a builder reaches its samples' data: it reads the organiser's observation cache, and
holds no reader of its own** (§6.4). **This reverses what this section said before, and the
reversal is the owner's, 2026-08-18.** The paragraph below is kept because its objection was
right and is now a cost this design accepts rather than one it avoids.

**What it said, and why it was overruled.** It said each region's task owns k cursors, one per
sample — the shape [`run_streaming.md`](run_streaming.md) §3.4 fixed ("a parallel loop gives each
in-flight region its own source per sample, over the sample's one open file") — and that a shared
per-sample window was considered and rejected, because *"it re-creates the global frontier this
design exists to remove: every builder's progress would couple through the window's trailing edge,
which is production's watermark (§10, lesson 5)."* §6.4, decided a day later, chose the shared
window anyway, for a reason this section had not weighed: a reader per builder per sample is
`builders × k` open cursors, where `k` is the cohort's sample count — at 3,000 samples and 16
builders, 48,000 of them, against 3,000.

**So the coupling is real and it is accepted, and this is its exact shape.** Because a builder
reads the cache while the organiser draws the readers forward, the two cannot run at once: the
organiser covers the ground for a **round** of regions, the round's builders run, and only when
every one of them has finished does the organiser resolve, release and evict. Every builder in a
round therefore waits for the slowest in that round. **That is a bounded frontier, not
production's**: it spans one round — `builders × cohort_locus_builder_regions_len` bases, 3,200 at
16 builders on 200-base regions — where production's watermark spanned the whole run and never
advanced past its slowest chunk. What it costs is the tail of each round; what it buys is a
reduction in open cursors and in decode work by a factor of the builder count, since two builders
whose regions land in the same psp block now share one decode instead of each paying for it.

**The memory product this section named still exists and is unchanged in kind.** It is
`samples × the ground the round spans`, which §6.4 states and §8 prices at both ends of the
range; the round replaces `look-ahead` as what sets it.

**What is not decided here** is whether the round is the final arrangement. A shape that lets one
builder start its next region while others are still working — an `RwLock` over the cache, or
windows handed out as owned copies — removes the round's tail at the cost of a lock or a copy per
region.

**The round's tail is not what costs, and it was the wrong thing to gate on.** This paragraph
used to say neither alternative should be reached for before the tail had been measured. What
was measured on 2026-08-18 instead was every phase of a round, with timers, at 200-base regions
and one region in flight per thread. On 8 threads the builders are 40–52% of the merge and the
organiser's own work — resolving overlaps, releasing loci, submitting outcomes — is **0.02 ms of
a 24 ms merge**. Everything else is the cache's two writers: drawing the readers forward 25–34%
and dropping what is past 16–32%. So **48–60% of the merge runs while every builder waits**, and
that is the round's *head*, not its tail. Fitting `time = serial + parallel ÷ threads` to the
one- and eight-thread times says the same from outside: 42–67%. It is why 4 threads buy only
1.4× and 8 buy 1.4–2.0× over one.

**Windows handed out as owned copies was then built and measured, and it does not pay.** The
driver took each round's window out of the cache as an owned vector — leaving behind a copy of
only the observations reaching past the round, which the next round chains through — so that the
builders read memory the organiser could not touch, and the cover for the next round ran beside
them. All 227 of the module's tests passed on it, output unchanged. What it did to the clock, on
8 threads against the arrangement above:

- **the overlap works**: covering the next round beside the builders took 12–14% off the merge
  at 1,000 and 3,000 samples on ground varying at about one position in a hundred;
- **the handover costs more than that**: 19–23% at the same two sizes, and 2–6% on ground
  twenty-five times denser;
- **net, it was 1–6% slower** everywhere but one cell, and it held **1.8 times the records** —
  the round the builders read plus the round the cache has drawn.

The handover cannot be made much cheaper as long as a sample's window is one vector: what it
copies is the observations that reach past the round's end, and those must be in the builders'
window *and* in what the cache keeps, so they are cloned once per sample per round. The version
that could pay is a different cache: per-sample **append-only storage whose addresses do not
move**, where a builder holds a range into memory the organiser is appending elsewhere in, so
the handover is a range and not a copy. That is a redesign of §6.4's cache, not a change to this
section, and nothing has been built or measured for it.

**The `RwLock` alternative is not a way round this.** A builder holds its window for the length
of its build and the organiser's draw would move that window, so the lock would be held for the
whole of one or the other — which is the round, arrived at through a lock.

**And the stage upstream is fourteen to twenty-three times the merge, which is what decides
whether any of this is worth building.** Producing a sample's observations and merging the
cohort's were timed side by side on the tomato benchmark
(`examples/ng_cohort_merge_real_cost.rs`, 2026-08-18): 16 accessions over 100 kb of SL4.0 cost
**2.21 s** in the generic locus generator against **95 ms** in the merge on 8 threads, and 63
accessions over 200 kb cost **12.09 s** against about **850 ms**. Both stages take threads —
the generator perfectly, since samples are independent files, the merge by 1.4× — so on 8
threads the generator is still about twice the merge. At **one** sample the gap is far wider:
HG002 over 76,530 bases cost 1.81 s to walk and 12.7 ms to merge, a factor of **142**.
Recovering the 12–14% the overlap is worth would move a whole run by roughly four parts in a
hundred at 16 samples and by nothing measurable at one. **Nothing in this section should be
built before the run it belongs to has been assembled and timed end to end.**

### 6.3 The organising thread: order, overlaps, emission

One thread receives every builder's output and does three things: it keeps the cache, it resolves
the overlaps of §6.1, and it emits in genome order.

Builders finish out of order; the consumer must see genome order. The reorder structure is
production's, carried whole: results keyed by region index in an ordered map, drained while the
head equals the next expected index — the `BTreeMap` the VCF writer drains on `next_expected`
([`vcf_writer.rs:168-176`](../../../../src/var_calling/vcf_writer.rs)).

**Resolving overlaps is why the drain waits for a gapless run of regions rather than emitting each
as it lands.** A locus can only be confirmed once the region before it has arrived, because that is
what says whether an earlier locus covers its start. So the organiser holds a region's loci until
its predecessor is in hand, applies the rule of §6.1 — earlier owner wins, overlapping loci and
everything inside a failed locus's span are dropped — and only then emits. The wait is one region,
not the whole convoy.

Only survivors are buffered — a builder fails and discards locally (§3.2, §4.3), so the buffer
holds about one locus per hundred positions of in-flight ground at the measured corner, not the
ground itself. But **every region delivers exactly one result — observations plus its
failed-locus count — even when both are empty**, because the drain advances only on a gapless
index sequence and the counts must all arrive to be summed (§3.3). Production distinguishes
exactly these things: variant-free spans emit nothing while the order counter stays gapless
([`cohort_integration.rs:1173-1179`](../../../../src/var_calling/cohort_integration.rs)), and the
gap is guarded at release level because a lost index silently truncates the output
(`MissingChunks`, [`vcf_writer.rs:152-158`](../../../../src/var_calling/vcf_writer.rs)). ng keeps
both: the empty-result rule and the release-level guard (arch doc §5).

Whether `call_vars_from_observation` runs inside the builder (buffering called variants) or after
the drain (buffering observations) is open question 5 — the choice commutes with everything above
because calling one cohort observation reads nothing outside it.

---

### 6.4 The observation cache — decided 2026-08-17

**Upstream produces observations in one forward pass; the builders want them in several places at
once. The organiser holds a cache between the two.**

A sample's observations arrive in coordinate order, from a walk that never goes backwards. But at
any moment several builders sit at different points of the genome, each wanting the observations
over its own region. Nothing can serve that by seeking, and giving every builder its own reader
would mean as many readers per sample as there are builders.

So there is **one reader per sample for the whole run**, advancing forward only, and the organiser
pulls from it into a cache:

- **Fill:** the cache holds every sample's observations covering the ground currently assigned to
  any builder. When the organiser hands out a region it draws the samples forward to cover it.
- **Evict:** when every builder that could own a locus over some ground has finished and the
  organiser has resolved and released it, the observations over that ground are dropped.
- **Builders read the cache and never write it.** They do not hold observations, do not seek, and
  do not talk to a reader. Goal 1 is unchanged: the only thing a builder produces is its result.

**This cache is where this module's memory goes**, and it is why the builders are kept close
together. What it holds is roughly

> `samples × (the ground the builders span + the reach of observations running past it)`

and the first factor is what the design controls: with `n` builders on regions of
`cohort_locus_builder_regions_len` bases, they span `n × len` bases of genome. At 16 builders and
200 bases that is 3,200 bases of ground; at 20 bases apiece it is 320, a tenth of the memory for
the same parallelism. The second factor is not controlled here — an observation may reach up to
`max_record_span` past where it starts (§5) — but wide observations are rare, so in practice the
first term is what moves.

**`cohort_locus_builder_regions_len` — a command-line parameter, default 200 bases.** The earlier
draft derived the width from `max_cohort_locus_span`, twice it; that was the wrong parent, because
the width's real cost is this cache and not anything about locus widths. Twenty was the owner's
starting value; two hundred replaced it on 2026-08-18, on the measurement §14 question 1 records.
What settled it is the other side of the trade, which the draft above did not have: a narrow
region makes the merge walk the whole cohort four times over — a cover, an eviction, a window and
a builder's set-up — whether or not the region holds a single record, and on ground where about
one position in a hundred varies a 20-base region holds a fifth of a locus.

## 7. Degradation at the edges

### 7.1 The dense region

The pathological input is ground where observations overlap wall to wall for a long distance —
multi-base spans chained end to end. Dense SNPs alone do not do it: a SNP covers one base, so
consecutive SNPs start separate loci (§4.1).

There, one locus can span many builders' regions. Its owner builds it, every builder inside its
span produces loci that lose the overlap, and their work is thrown away. The component degrades
toward serial behaviour over that ground and only that ground.

What bounds the damage:

- **Length: the segment.** A generic segment's end is always a cut, so a region never outgrows
  its segment, and the boundary scan always terminates there. This is why production's
  `StalledCut` has no counterpart here: production's cut could fail to advance and spin, so it
  guards the progress invariant with a release-level error
  ([`cohort_integration.rs:395-400`, fire site `:1145-1150`](../../../../src/var_calling/cohort_integration.rs));
  ng's progress guarantee is structural — the reference-defined segment — and production had no
  segments to lean on (§10, lesson 2).
- **Memory: the look-ahead.** At most `look-ahead` regions are in flight; when the slow region
  is the yield frontier, later builders idle rather than pull ahead
  ([`run_streaming.md`](run_streaming.md) §3.5). Peak memory does not grow with the stall; wall
  does, by the serial cost of the dense run.
- **Wasted work: the duplicated scans** (subdivision only). Builders at provisional points inside
  the dense run each scan position summaries to the same far cut. Position summaries only — cheap
  relative
  to one build — and bounded by the segment length times the in-flight count.

How long can a generic segment be? Unmeasured — [`run_streaming.md`](run_streaming.md) open
question 1 (kilobases at the routing floor is its leaning). This design inherits that number as
its worst-case serial unit; subdivision (§6.1) is its one mechanism, and it stops at the shape
with no cut, deliberately: going further would need cross-builder coupling inside an open locus,
which is the one thing this design refuses.

### 7.2 One sample

At k = 1 the fold of §4.1 is a copy of the sample's own position summary, and collation holds one
member list. Nothing else changes, deliberately: grouping, `max_cohort_locus_span` verdict, the
variability
filter and the assembly are calling-side work a single-sample run needs regardless — they are
not merge overhead, and there is no cheaper "do nothing" that still yields callable cohort
observations. The overhead unique to this component at k = 1 is the per-region result markers
(§6.3) and, under subdivision only, the boundary scans (§6.1) — per region, never per read.
Production's sharpest ordering bug was at exactly k = 1 across many intervals
([`cohort_integration.rs:1905-1928`](../../../../src/var_calling/cohort_integration.rs) — the
single-sample cohort was "what actually fails"), so k = 1 through the full parallel path is a
first-class oracle here, not a degenerate afterthought (§15).

At k = 0 there is nothing to iterate; construction of a caller with zero samples is refused at
the run level (production instead yields an empty stream —
[`:1997`](../../../../src/var_calling/cohort_integration.rs); ng prefers the refusal because a
zero-sample *calling run* is always a caller bug, where production's producer was also test
plumbing).

### 7.3 The large cohort, honestly

**This section used to warn that the keep rule would erode as the cohort grew, and it was
right — that was the cohort-sum rule, and it was replaced on 2026-08-19 (§4.3).** The
measurement it asked for was made on the tomato panel, one 100 kb interval, and it is the
reason the rule changed. Loci discarded out of every 1,000 the merge closed, under the old
rule and the new:

| accessions | closed | old rule (cohort sum ≥ 2) | new rule (some sample reaches it) |
|---|---|---|---|
| 1 | 96,082 | 997 in 1,000 | 997 in 1,000 |
| 16 | 98,726 | 878 | 960 |
| 63 | 97,408 | 439 | 876 |

The old rule's saving fell away as samples were added — at 63 accessions it built a cohort
observation at 546 positions per kilobase — because a fixed bar of 2 summed over the cohort is
reached by sequencing error somewhere as soon as the cohort is large. The new rule does not
drift: the question it asks of a sample does not change when a sample that carries nothing joins
the run, so the kept fraction tracks how much variation the cohort holds rather than how many
samples it has. Between 16 and 63 accessions it still fell, from 960 discarded in 1,000 to 876,
and that is real variation being found — three times the accessions carry more segregating
sites — not the rule weakening.

What remains true from the old warning is the shape of the cost if the kept fraction does climb:
wall and reorder-buffer occupancy, not correctness and not unbounded memory (§8's bounds do not
depend on the keep rate; the buffer holds at most `look-ahead` regions' survivors).

**Depth is the other axis and it is now covered too**, which the fixed bar never was: on the
GIAB trio at 313 compared reads a sample, three samples were enough for the old rule to build
57 loci in 100. The share brings that to 2 in 100 (§4.3). Neither number has been checked
against a truth set — see §4.3's last paragraph.

---

## 8. Where the memory goes

**One term dominates: the observation cache (§6.4).** Everything else is small or bounded by a
constant.

| what | size | source of the number |
|---|---|---|
| the observation cache | `samples × observations over the builders' span` — the span being `builders × cohort_locus_builder_regions_len`, 3,200 bases at 16 builders and the 200-base default, plus the tail of observations reaching past it | the formula is §6.4's; measured on fabricated ground with a record every hundred bases at 1,000 samples and 16 in flight: **33 records held per sample**, against 4 at the old 20-base default. What one record costs in bytes is still **unmeasured** |
| the current locus, per builder | one `max_cohort_locus_span`-wide window per sample, held while it is assembled (§4.5) | bounded by `max_cohort_locus_span`, never by the region |
| the cohort summary, per builder | ~24 B per position of its region (production's three-column layout, [`cohort_integration.rs:64-78`](../../../../src/var_calling/cohort_integration.rs)) | ~4.8 kB at a 200-base region; independent of the sample count |
| survivors awaiting release | about 1 locus per 100 positions of resolved ground (§4.3); a failed locus adds a span and a counter, no data (§3.2) | the 28,718 / 2.83 M measurement, 50 tomato samples at ~3× |

**The cache is why the builders are kept close together, and why the region width is bounded at
all.** Ten times the span is ten times the cache, for the same number of builders and the
same parallelism.

**What cannot be priced yet is what one observation costs**, which is the first factor in the only
term that matters. Nobody has measured it, so no total in this section is given: multiplying an
unmeasured per-observation size by a sample count would produce a figure that looks like a
measurement. Open question 1 (§14) sweeps the region width and should report the cache's peak
directly, which measures the product rather than either factor.

## 9. Determinism — why the output cannot depend on the builders

Three facts, each carrying part of the property, together the whole of it:

1. **The set of cohort loci — built, failed, and dropped alike — is a function of the per-sample
   streams and the configuration.** The fold is order-independent (union + `max`,
   production-tested with deliberate ties,
   [`cohort_integration.rs:1960`](../../../../src/var_calling/cohort_integration.rs)); the
   closure is a deterministic left-to-right walk of the folded summary (§4.1); both per-locus
   verdicts — `max_cohort_locus_span` and the variability filter — read only that summary (§3.2, §4.3);
   members are selected by span overlap (§4.2). No scheduling input appears anywhere in the
   chain, so the same input gives the same failed loci, and the same counts, at any builder
   count — the failure path is covered by the same argument as the happy one, because failure is
   decided on the same fold. Recovery is equally fixed: loci are disjoint (§4.1), so the locus
   after a failed one starts at the first folded position past the failed locus's furthest
   reach, a position the fold determines, not the builder.
2. **Where the work was cut changes no locus.** The chaining rule reads observations, not region
   boundaries (§4.1), so any
   partition into regions yields the same loci with the same members and the same verdicts —
   partition-invariance, production's own strongest producer test (identical output at chunk
   targets 1, 3, 17 and 100,000 —
   [`cohort_integration.rs:1577-1588`](../../../../src/var_calling/cohort_integration.rs)). The
   division into building regions is itself fixed (boundaries are
   functions of the data — §6.1), but the property would hold even if it were not.
3. **Emission order is region index, then position.** The ordered emitter (§6.3) makes the
   yield sequence the serial sequence regardless of completion order; within a region the
   builder emits in genome order; loci are disjoint, so position is a total order on them.

Together: the same bytes at 1 builder and at 16, at any look-ahead, at any `cohort_locus_builder_regions_len` —
the properties [`run_streaming.md`](run_streaming.md) §12.2 requires, plus partition-invariance,
which is this component's own and is tested directly (§15).

---

## 10. Lessons from production, by name

Production's producer is the closest thing to a prior implementation of this component. What
transfers, what does not, and why — the "why not" mattering because ng has reference-defined
segments and production had none.

1. **Cut where a bound proves safety; never speculate.** Production derived its independent
   intervals by gap-merging block ranges under `max_group_span` — two blocks join unless the gap
   between them exceeds it, "so every interval boundary is a safe gap no variant group can span"
   ([`merge_block_ranges`, `cohort_integration.rs:403-430`](../../../../src/var_calling/cohort_integration.rs))
   — and cut chunks only at clean group boundaries
   (`find_cut`, [`:266-287`](../../../../src/var_calling/cohort_integration.rs)). Nothing was ever
   built and revoked; safety was decided from cheap metadata before any real work. Carried
   **Not carried.** This design does the opposite — it builds, overlaps, and lets the organiser
   resolve (§5), because production's approach needs a backward scan bounded by how wide an
   observation can be, and at a 200-base region that scan is tens of times the work it
   protects. Production's own assumption behind the gap constant — that it is at least any
   record's span — lived only in a test-fixture comment
   ([`:1664`](../../../../src/var_calling/cohort_integration.rs)); ng promotes it to a named,
   header-recorded value (§5). The alternative — build speculatively, reconcile overlaps
   between workers after — puts coordination exactly where the data is hardest, and production's
   history contains no version that needed it.
2. **`StalledCut`: guard progress at release level — but ng's guarantee is structural.**
   Production's cut can fail to advance on a degenerate fold, and a `debug_assert!` was promoted
   to a release error because a spin is worse than a stop
   ([`:395-400`](../../../../src/var_calling/cohort_integration.rs)). Not carried as an error:
   ng's boundary scan terminates at the segment end by construction (§7.1). Carried as a test:
   the scan's bound is asserted, and the dense-region shape is a fixture, not a surprise.
   Production needed the guard precisely because, with no reference-defined structure, no
   position was *guaranteed* cuttable.
3. **Decide from the cheap facts, assemble only what survives.** Production reaches its verdicts
   from per-row summaries and touches the expensive evidence only for the rows it keeps
   ([`TwoPhaseSegment`, `sample_reader.rs:698-712`; `set_variable_rows`, `:789`](../../../../src/var_calling/sample_reader.rs)).
   That ordering is carried here — §4.3 and §4.4 are it, per builder, and §3.2 extends it to the
   failure verdict. **What it saves differs by path:** in the direct path the observations are
   already in hand, so what is saved is the assembling and everything after it; in the psp path
   the same ordering additionally means the discarded evidence is never decompressed. Also
   carried: production keeps an eager whole-segment decode alive purely as the byte-identity
   oracle its deferred path is checked against
   ([`sample_reader.rs:20-26`](../../../../src/var_calling/sample_reader.rs)) — §15 keeps that
   oracle shape.
4. **Reorder by index, drain on next-expected, and treat a gap as a release-level bug.** The
   writer's `BTreeMap` + `next_expected` cursor and the `MissingChunks` guard
   ([`vcf_writer.rs:152-176`](../../../../src/var_calling/vcf_writer.rs)). Carried whole, with
   the index moved from chunk order to region index and the "exactly one result per region,
   empty included" invariant stated (§6.3) — a result now carrying a count as well as
   observations.
5. **The watermark is the anti-pattern this design exists to remove.** Production advanced every
   sample behind `min` over per-sample coverage — the cohort watermark
   ([`cohort_integration.rs:920-939`](../../../../src/var_calling/cohort_integration.rs)) — and
   everything downstream of that single frontier was serial; its five-stage pipeline
   ([`pipeline.rs:1-30`](../../../../src/var_calling/pipeline.rs)), bounded queues
   ([`:98-111`](../../../../src/var_calling/pipeline.rs)) and measured thread split
   ([`resolve_split`, `:135-168`](../../../../src/var_calling/pipeline.rs)) are machinery for
   overlapping work *around* a barrier ng does not have. Not carried, and nothing like it may
   reappear: a shared per-sample read window (§6.2) would be this frontier under another name.
   What is kept from that history is one measured direction: finer work units balanced better
   ([`:86-97`](../../../../src/var_calling/pipeline.rs)), which informs `cohort_locus_builder_regions_len`'s sweep.
6. **The test shapes are the transferable specification of correctness.** §15's table maps each
   production test to the ng test it becomes; the two that found real bugs — partition-invariance
   and the k = 1 multi-interval ordering test
   ([`:1577`, `:1905`](../../../../src/var_calling/cohort_integration.rs)) — anchor the suite.

---

## 11. Traps — what will bite the coder

- **Merge the samples' observations by position before walking them; do not loop over samples at
  each position.** A deletion widens the locus while it is being closed, and the records the
  widening brings in sit at later positions in other samples. Walking one merged, position-ordered
  stream absorbs them as it goes. Looping over samples instead means every widening sends you back
  over samples already visited, repeatedly until the span stops growing — correct, but it turns one
  pass into a fixpoint, and the shape that needs it is the natural one to write.
- **The observation cache must hold what overlaps a region, not only what starts in it.** A
  builder given only the observations beginning inside its own ground cannot see a locus that
  opened earlier and reaches in, so it builds neighbouring loci from a partial picture — and
  the organiser cannot catch it, because those loci overlap nothing it was told about. The
  failing input needs a wide deletion beginning before a region and reaching into it; no test
  that varies only SNP fixtures will produce one (§6.4).
- **The boundary scan and the builder must read the same position summaries** (subdivision only).
  `boundary(g)` is computed independently by adjacent builders (§6.1); any divergence between the
  scan's view and the build's view — a filter applied in one and not the other — makes two
  builders resolve different boundaries, and a locus is built twice or not at all, silently. One
  code path for light-column access, used by both.
- **Reach arithmetic must saturate.** Production's `reach` saturates on both operations
  ([`cohort_integration.rs:46-48`](../../../../src/var_calling/cohort_integration.rs)); a
  non-saturating rewrite overflows at contig ends and the cut check quietly inverts.
- **A failed locus is a count; an over-ceiling observation is an error** (§3.2). Collapsing the
  second into the first hides file corruption inside a calibration counter; collapsing the first
  into the second lets one long deletion kill a run.
- **The per-region counts must be summed at the drain, not sampled.** The same failure shape as
  the read-filter tallies ([`run_streaming.md`](run_streaming.md) §8): per-builder counts that
  never get summed under-report by the worker count, and every number stays plausible.
- **Cursor order within a builder must stay monotonic.** Sources answer any order but pay a seek
  and a block decode going backward ([`run_streaming.md`](run_streaming.md) §8); a builder that
  interleaves the summary and assembly requests non-monotonically across samples turns every locus
  into k seeks.
- **The empty region's result is load-bearing** (§6.3). Dropping it stalls the drain at the
  first variant-free region — and now also loses a count. Production promoted exactly this class
  of invariant to release errors after silent truncation
  ([`vcf_writer.rs:152-158`](../../../../src/var_calling/vcf_writer.rs)).
- **The generator's region halo is not the merge's business.** The generic generator queries a
  halo past its region and clamps on the anchor
  ([`pileup/genome_walk.rs:200`](../../../../src/ng/locus_generation/pileup/genome_walk.rs));
  an observation belongs to the region containing its start position, and that
  is what the organiser's overlap rule rests on. Re-filtering members by any other rule double-counts at
  region edges.

---

## 12. Cross-cutting concerns

**Errors.** A builder's failure surfaces through the run's existing shape — an error naming the
sample and the span, ending iteration ([`run_streaming.md`](run_streaming.md) §9). This component
adds one error variant — an observation wider than its file's recorded reach ceiling (§3.2, names
sample and region) — plus the release-level gap guard on the drain (§6.3). A failed locus is
**not** an error: it is a counted outcome of a healthy run (§3.2).

**Concurrency.** Shared and read-only: the segmentation, the reference, the parameters, each
sample's open file. Owned per builder: k cursors, position summaries, the fold, the locus under
construction, the region's counts. Single-threaded by construction: the drain at the yield
point, where the counts are summed. No lock anywhere in the component; a lock is a defect (§6.2).

**Performance.** Two knobs reach this component: the run's look-ahead (memory) and — under
subdivision only — `cohort_locus_builder_regions_len` (balance). `max_cohort_locus_span` is policy, not tuning; the
reach ceiling is the generator's, not this component's. Any other tuning constant appearing in
the implementation is a defect ([`run_streaming.md`](run_streaming.md) §9 states the same rule
for the run).

---

## 13. Deferred, with a recommended home

- ~~**The psp header field for the observation reach ceiling**~~ — **landed 2026-09-04**
  ([`run_streaming.md`](run_streaming.md) §6.1's `observation_reach_ceiling_bp`, written by the
  walk stage and read by the calling stage at open, which takes the maximum over the cohort's
  files: `run_driver_psp_mode.md` steps A3 and E4). No refusal accompanies it, as this bullet
  said none would. **What a reader does with it is still nothing**: the ceiling bounds how far
  the observation cache may have to reach (§5, §6.4), and that cache grows on demand rather than
  taking a capacity — so the number is exposed for the psp-mode performance work rather than
  consumed by the merge.
- **The position summary's encoding, carrying the reference span** — to the psp encoding spec
  ([`run_streaming.md`](run_streaming.md) §10). Phase 1 must serve position, reference span and
  non-reference count per record (§1.3); the sketch in [`run_streaming.md`](run_streaming.md)
  §3.3 lists depth instead of span and should be corrected when the encoding spec lands.
- **Where the run summary lives, and the failed-locus count's exact surface in it** — to the
  emission step's document, which owns what a calling run reports; this spec fixes only that the
  count exists, is summed at the drain, and reaches the user (§3.3).
- **Choosing candidate alleles from the table, and their output representation** — to the calling
  steps' spec. This component unifies what the samples showed (§4.2); which of those alleles are
  worth calling, and how they are written, are decisions about calling. Production's variant
  grouping is the reuse candidate ([`run_streaming.md`](run_streaming.md) §10).
- **Calling deletions longer than `max_cohort_locus_span`** — if ever wanted, a pass over the
  emitted records,
  the home [`run_streaming.md`](run_streaming.md) §4.3 already names for cross-segment events.
  Never a coupling between in-flight regions.

---

## 14. Open questions

1. **What should `cohort_locus_builder_regions_len` be?** — **200 bases**, set by the owner on
   2026-08-18 (§6.1); the number is settled, the reason is measured on fabricated ground only.

   The draft of this question named one cost of a narrow region — more joins between builders,
   and so more overlapping work discarded — and it turned out not to be the one that decides.
   **A narrow region makes the merge walk the whole cohort four times over** — drawing the
   readers forward, dropping what is past, handing out a window, and the arrays a builder
   allocates before it reads anything — **whether or not the region holds a single record.** On
   ground where about one position in a hundred varies, four 20-base regions in five hold no
   record at all. Measured over 20,000 fabricated bases with a record every hundred, in a
   release build on 8 threads with one region in flight per thread: 63 samples took 12.2 ms at
   20 bases and 2.2 ms at 200, 1,000 samples 40.3 ms against 26.0, and 3,000 samples 177 ms
   against 113. Against that, the cache holds **33 records per sample rather than 4** at 1,000
   samples and 16 regions in flight — bounded by the round either way, and eight times as much
   of it.

   **The tomato cohort was then measured, and it changes the reason without changing the
   value.** `examples/ng_cohort_merge_real_cost.rs` walks the benchmark's 63 accessions through
   the generic locus generator over 100 kb of SL4.0 and merges what comes out. Two things it
   found:

   - **Observations arrive about one per covered base per sample, not one per varying
     position.** 96,605 records per sample over 100,000 bases — the generator emits at every
     position its reads cover, and the keep rule that discards the quiet ones runs inside the
     merge (§4.3), not before it. The fabricated ground this question was settled on had a
     record every hundred bases, so it was **a hundred times too thin**, and every ratio taken
     on it overstates what the width is worth.
   - **On real observations the width barely moves the merge on one thread and decides whether
     threads help at all.** One thread, 63 samples: 656 ms at 20 bases, 615 at 100, 616 at 200,
     624 at 500, 636 at 1,000 — 6% across a fiftyfold range. Eight threads at 16 samples: 173 ms
     at 20 bases against 130 ms on one thread, so **a 20-base region makes threads a
     pessimisation**, and 93 ms at 200, which is 1.4× one thread. The eight-thread optimum is
     100–200 bases at both cohort sizes measured.

   So 200 stands, and what it is for is the organiser: the per-region cost falls on the one
   thread that covers and evicts, so a narrow region starves the others.

   **The opposite corner says the same thing about density and moves the optimum.** HG002 at
   high coverage, one sample, over 20 of the GIAB benchmark's intervals (76,530 bases): again
   **one record per covered base** — 76,141 of them — and again **no building region empty at
   any width**, because high coverage leaves no gaps. The width optimum there is 500 bases at
   both 4 and 8 threads (11.7 and 12.7 ms, against 15.9 and 20.6 at 200), where at 63 samples
   it is 100–200. So the best width falls as the cohort grows, and 200 is within a tenth of the
   best at 16 and 63 samples and 1.6× off it at one — where the whole merge takes 20 ms.

   **That run also shows what this module is for, on real ground.** Twenty scattered intervals,
   one sample: the oracle takes **62.6 ms** and the cached driver **16.6**, because the oracle
   hands every analysed region the whole set of observations and closes the prefix in front of
   it each time (§6.4). Nearly four times, at one sample, on twenty intervals — the effect the
   C1 and C2 reviews measured on fabricated ground, now on real.

   **What is still owed** is the **discard rate** at the joins, which wider regions can only
   lower and which no measurement here has counted, and a cohort of *thousands*, which nothing
   real reaches: 63 accessions is the largest measured. The builder-idle profile this question
   originally asked for
   ([`pipeline.rs:86-97`](../../../../src/var_calling/pipeline.rs)) is superseded by the phase
   timings in §6.2, which say where the merge's time goes without one.
2. **When a sample has two separate observations inside one locus, is its allele the combination
   of both?** — OPEN, and it is the question projection raises that byte-equality does not answer.
   A sample with a SNP at one position and another three bases along, both inside one cohort locus,
   projects to a sequence carrying both changes — but only if the same reads carried both. If they
   sat on different haplotypes, that compound allele was never observed in any molecule.
   **Production answers it and ng's observations carry what is needed to do the same:** it admits a
   cross-record compound *"only when chain-id evidence inside at least one sample links the
   constituents"*, and where a sample's chain is broken at a compound another sample anchored, it
   falls back to treating the constituents independently and records that it did
   ([`per_group_merger.rs:1-20`](../../../../src/var_calling/per_group_merger.rs)); ng's
   `SequenceObservation` carries `chain_ids` for reads that disagree with the reference
   ([`locus_generation/mod.rs:157`](../../../../src/ng/locus_generation/mod.rs)). *Leaning:* follow
   production — compound only on chain evidence, constituents-independent otherwise, and the
   fallback recorded rather than silent. **Settled by:** deciding it before the builder is coded,
   because it changes what a projected allele *is*; it is not a tuning question.
3. **Is 50 the right default bound?** — OPEN; the value is the owner's, unmeasured — soft, and
   cheap to revisit because re-calling under a new bound needs no re-walk (§3.1). **Settled by:**
   counting, in one walk of a tomato sample and HG002, reads carrying deletions wider than 50
   bases (the events the caller gives up) and the length distribution just below — how much real
   signal sits within reach of a larger bound, against §8's assembly window growing with
   `max_cohort_locus_span`.
4. **Does the far end need a cohort-scaled keep threshold?** — **SETTLED 2026-08-19, and the
   answer is no in both directions.** The measurement this asked for was made (§7.3): the
   cohort-*sum* rule it was written about had the defect from the other side — its bar was fixed
   while the sum grew, so it stopped filtering as samples were added. The rule is now asked of
   one sample at a time, which neither scales with the cohort nor needs to. Scaling the bar with
   the cohort was also tried on paper and rejected for the reason the old leaning gave: at 1,000
   samples a 2% share of the cohort's reads asks for about 100 carriers, which drops exactly the
   rare variants a large cohort is run to find (§4.3). What is still owed is the recall check —
   §4.3's last paragraph, and it is a different question from this one.
5. **How do failed loci surface beyond the count?** — OPEN; leaning: log each failed locus's
   span at warning level up to a small cap, then count only — deterministic, since the drain
   emits in region order (§9). A machine-readable sidecar (a BED of failed spans) is the step up
   if operators need to intersect the refused ground with annotations. **Settled by:** the owner,
   when the run-summary surface is specced (§13).
6. **Does calling run inside the builder or after the drain?** — OPEN; leaning inside — the
   buffer then holds called variants, smaller than observations, and
   [`run_streaming.md`](run_streaming.md) §3.5's skeleton already collects per-region calls.
   Commutes either way (§6.3). **Settled by:** when the emission step fixes `Variant`'s shape,
   compare the two buffers' sizes at the §8 corners.

---

## 15. How we know it works

Production's tests around the cut and the intervals are the closest thing to a specification of
correctness that exists; each row names the ng test it becomes. All cited tests were read at the
lines given. The failed-locus rows have no production ancestor — the ruling is new — and are
listed after the table.

| production test ([`cohort_integration.rs`](../../../../src/var_calling/cohort_integration.rs) unless noted) | what it pins | the ng test it becomes |
|---|---|---|
| `streaming_matches_reference_across_chunk_sizes` `:1577` | identical output at chunk targets 1, 3, 17, 100,000 | **partition-invariance**: one-shot whole-segment build equals region-wise builds at any partition, and at 1–16 builders — observations *and* failed counts equal — this component's regression anchor |
| `staged_channel_path_matches_owned_across_intervals` `:1905` | k = 1 across many intervals equals the serial reference — the shape that found the stale-marker desync | k = 1 through the full parallel path over many regions and segments equals serial; run first, not last |
| `fold_unions_positions_and_maxes_aggregates` `:1302`, `merge_reduce_tree_is_order_independent_with_ties` `:1960` | fold commutes and associates, ties included | same properties on ng's fold, ties at shared positions included |
| `keep_threshold_one_is_variant_filter` `:1318` | threshold-1 keep is a variant filter | the variability filter keeps exactly the loci with a candidate position |
| `over_approximation_is_max_not_sum` `:1335` | production's keep is a `max` over samples, not a sum | **followed, not inverted, since 2026-08-19**: two samples with one non-reference read each at one position are **not** kept, because no single sample reached two (`one_non_reference_read_in_each_of_two_samples_does_not_reach_the_threshold`). ng adds a share of the sample's own reads on top, which production has no counterpart for (`the_share_decides_once_depth_makes_it_the_larger`) |
| `multi_position_group_kept_whole` `:1344` | a kept group keeps its quiet member positions | same: a kept locus carries its reference-only member positions |
| `interval_with_no_variants_yields_no_chunks` `:1640` | variant-free ground emits nothing, order stays gapless | an all-reference region emits no observation and still delivers its result; the drain's gap guard fires when it does not |
| the eager-decode oracle ([`sample_reader.rs:20-26`](../../../../src/var_calling/sample_reader.rs)) | two-phase equals eager, byte for byte | an eager whole-region build, tests only, as the oracle the two-phase builder is compared against |
| `produce_chunk_with_zero_samples_yields_none` `:1997` | zero samples does not spin | zero samples is refused at caller construction (§7.2) |

New tests the ruling requires, with no production ancestor:

- **The failed locus is counted and suppressed whole.** One sample's over-wide deletion in a
  cohort where another sample has a SNP inside its span: no observation over the locus's ground,
  the count is 1, and the loci on either side are unchanged (§3.2).
- **The chained case fails identically.** A chain of narrow observations whose closure exceeds
  `max_cohort_locus_span`: same verdict, same count, no single member over it (§3.2).
- **The failed set is scheduling-invariant.** The two fixtures above at 1, 2, 8 builders and at
  two partitions: identical failed counts and identical surviving streams (§9).
- **The width verdict comes before the variability one.** A reference-only chain wider than
  `max_cohort_locus_span` counts
  as failed, not silently dropped (§4.3).
## 16. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| reach arithmetic | `reach` ([`cohort_integration.rs:46-48`](../../../../src/var_calling/cohort_integration.rs)) | copied — saturating form and all (§11) |
| the cohort fold | `CohortSpanFold` ([`cohort_integration.rs:64-149`](../../../../src/var_calling/cohort_integration.rs)) | the layout and the union+`max` algebra, re-implemented over ng's position summaries; the order-independence tests come with it |
| grouping by shared bases | `derive_is_kept`'s group walk ([`:166-187`](../../../../src/var_calling/cohort_integration.rs)) | the chaining rule, unchanged; ng adds the per-locus bound verdict after it (§4.1) |
| the reach ceiling's existing owner | `max_record_span` ([`pileup/generator.rs:93,141,45`](../../../../src/ng/locus_generation/pileup/generator.rs)) | read, recorded in the psp header, never re-decided here |
| the psp path's deferred decode | `TwoPhaseSegment`, `set_variable_rows` ([`sample_reader.rs:698,789`](../../../../src/var_calling/sample_reader.rs)) | the shape the psp source must meet when that path is built; the direct path defers nothing, and the encoding spec owns ng's bytes |
| ordered emission | reorder `BTreeMap` + `next_expected` + `MissingChunks` ([`vcf_writer.rs:152-176`](../../../../src/var_calling/vcf_writer.rs)) | carried whole, keyed by region index; the result gains a count (§6.3) |
| per-sample stream heads | `MergedRegionReads` ([`sample_reads.md`](../arch/sample_reads.md) §4) | the argmin-over-heads layout for walking k position summaries in step during the fold |
| the flowing item | `SampleLocusObservations` ([`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | members of a cohort observation, unchanged |
| the segments and their boundaries | `TypedRegion`, `RegionKind` ([`region_typing/mod.rs:144,168`](../../../../src/ng/region_typing/mod.rs)) | the first build's whole partition (§6.1); their boundaries are the unconditional safe cuts (§5) |

**The parity oracle for the whole document is partition-invariance (§15, first row)** — the same
cohort observations and the same failed-locus counts from one builder over the whole genome and
from many builders over any safe partition — with the run-level mode-equivalence oracle
([`run_streaming.md`](run_streaming.md) §12.3) standing behind it.
