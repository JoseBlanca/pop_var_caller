# ng — direct mode: alignment files to called loci, in one process

**Status:** draft, 2026-08-28. The build order for **direct mode**: the object that takes every
sample's alignment files plus the run's parameters and yields called loci in genome order, and the
command that drives it. Design is settled in
[`../spec/run_streaming.md`](../spec/run_streaming.md) (the run's shape, §3 and §5.1) and
[`../spec/cohort_merge.md`](../spec/cohort_merge.md) (the merge it drives). This turns that design
into build order; it is **not** a place for new design.

**Why direct mode and not psp mode.** The spec chose this order and gave three reasons
([`run_streaming.md`](../spec/run_streaming.md) §2): it needs no file format, so it freezes nothing
while the psp encoding is still being measured; it is the shortest path to ng calling from real
reads; and once psp mode exists, direct mode is its oracle.

**Where this plan stops — ⛦ revised 2026-08-31, and it now goes further.** The draft stopped at
**called loci** because the emission step had no document. **It has one**: the format is settled in
[`../spec/vcf_output.md`](../spec/vcf_output.md) and the writer is coded, `src/ng/vcf/`, with
bcftools reading its output (owner's ruling, 2026-08-31). So **Milestone D goes through to VCF
records**, which is what [`../arch/run_streaming.md`](../arch/run_streaming.md) §3.4 already gives
both callers as their `Iterator::Item`. Expect to tweak that writer rather than to wait for another
plan.

---

## Scope

**In:** `src/ng/run/`'s two remaining pieces — a walker over one sample's alignment files behind the
merge's source interface, and the `AlignedFilesVariantCaller` that drives the merge and calls what
it yields; the wiring of `call_locus` into the merge's builder; the construction checks; and the
`call-from-alignments` subcommand.

**Out (later plans, or blocked):**

- ~~**The VCF writer and the `Variant` record**~~ — **no longer out, 2026-08-31.** They were listed
  as blocked on the emission step's spec; that spec is [`../spec/vcf_output.md`](../spec/vcf_output.md)
  and the writer is `src/ng/vcf/`. Milestone D reaches VCF records, and what it needs from the
  writer is a tweak rather than a build.
- **Everything psp** — the walk stage, the psp writer and reader, `PspVariantCaller`,
  `generate-psps`, `call-from-psps`, `generate-census`. Blocked on the encoding, which is under
  measurement ([`psp_encoding_experiments.md`](psp_encoding_experiments.md)).
- **Repeat tracts end to end** — candidate selection at a tract is specified and unbuilt
  ([`candidate_alleles_ssr.md`](candidate_alleles_ssr.md)). This plan's fixtures are the generic
  path; a tract in the analysed ground is routed and calls through the machinery that exists, but
  the plan does not claim the tract path is finished.
- **Splitting one sample's walk across workers** — goal 3 is unmet and §11 question 8 owns it.
- **The merge's own parallelism** — built, off by default, and
  [`../research/cohort_merge_parallel_cost_plan.md`](../research/cohort_merge_parallel_cost_plan.md)
  owns whether it is switched on.

## Principles (how the order was chosen)

- **Types first, then implementation**, within every milestone (project rule).
- **The heart before the plumbing, and here the heart is already built.** The merge and the calling
  loop exist and are tested; what is missing is the join. So the order is: make a walker look like a
  source (B), prove the merge reads it (C), join calling to the merge (D), then parallelise (E),
  then the command (F).
- **Each earlier stage is the next one's oracle.** The merge's single-threaded driver is the oracle
  for the walker-fed merge; the serial caller is the oracle for the pooled one. No stage is verified
  against itself.
- **Reuse over rewrite.** `SampleLocusObservationsIterator` already walks one sample's segments and
  yields observations
  ([`src/ng/locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs)); `call_locus`
  already turns one cohort locus into genotypes. Neither is re-derived — this plan adapts the first
  to the merge's interface and calls the second where the spec says the placement commutes.
- **Isolate the silent failures.** Two steps here produce wrong genotypes rather than a crash: the
  sample-order join (D2) and the concurrency invariance (E2). Own commits, oracle green before and
  after.
- **Incremental, with pauses.** One milestone, then stop for review.
- **Container builds.** All `cargo` through `./scripts/dev.sh` by absolute path (`CLAUDE.md`).

## Preconditions (already in place)

- **The cohort merge**, 13,285 lines under
  [`src/ng/run/cohort_merge/`](../../../../src/ng/run/cohort_merge/): the builder, the observation
  cache, the organiser's ordered release, and two single-threaded drivers that are the oracle for
  everything else.
- **The calling loop** — `call_locus` and everything under
  [`src/ng/calling/`](../../../../src/ng/calling/), through step E3b. One cohort locus in, genotypes
  out.
- **Candidate selection on the generic path** — `select_generic`, merged.
- **The per-sample walk** — `SampleLocusObservationsIterator`, which drives the typed-region
  generators over one sample's reads.
- **The parameters** — [`parameters_file.md`](parameters_file.md)'s Milestone F, so a run can be
  given its numbers or take the defaults. **This plan cannot start before that one reaches E2.**

**⚠ Two preconditions were not in place, and the executor must not work around them. One is now
done; the second still blocks F1.**

- ~~**The architecture document is stale.**~~ **Amended 2026-08-31**, and A1 is unblocked.
  [`../arch/run_streaming.md`](../arch/run_streaming.md) had specified `call_vars_in_segment`, a
  `LookAhead` knob, a segment pool and per-segment per-sample walkers — all retired by
  [`run_streaming.md`](../spec/run_streaming.md) §3.1 and §3.5. It now describes the built source
  trait, the serial merge feeding a pool of callers, and the two refusals A2 needs; its own
  amendment section lists what moved. **One question it records rather than answers**, and it is
  E1's rather than A1's: which thread the call runs on once several loci are called at a time —
  the merge's own region builders, which exist and are switched off, or a separate set of workers
  the merge feeds (its §8). It also surfaced A2's middle check, which the owner then ruled on;
  see that step.
- **The subcommand names are agreed and written nowhere.** `generate-psps`, `generate-census`,
  `call-from-psps`, `call-from-alignments` (owner, 2026-08-28). They belong with the rest of the
  command surface, which [`../spec/typed_regions_cli.md`](../spec/typed_regions_cli.md) owns.
  Record them there before F1.

---

## The steps

### Milestone A — the object, constructed but inert

✅ **A1. `AlignedFilesVariantCaller`'s type and construction.** Every sample's alignment files, the
segmentation's ingredients, the parameters in; the object holds one open `SampleReads` per sample
and the shared read-only state. No iteration yet.
*Depends:* the arch amendment above. *Source:* [`run_streaming.md`](../spec/run_streaming.md) §5.1.

✅ **A2. The construction checks.** Six refusals at construction, three of them before a single
file is opened:

- **a cohort of no alignment files** — a file pattern that matched nothing would otherwise open a
  caller over zero samples and die later inside the parameter assembly, a panic rather than a
  message;
- **parameters assembled for another cohort** — one inbreeding coefficient per sample of this run
  and one calibration per read group, each count checked with the other held fixed;
- **the file-descriptor headroom**, counted over *files* and not samples, so a run that would die
  at `EMFILE` refuses now, showing the arithmetic and the `ulimit` command;
- **the run's two views of its own reference agreeing** — the one the files were opened against
  and the one carrying the checksums are meant to be one reference at two moments, and the
  comparison downstream walks them in step;
- **the repeat catalog having been built on this run's reference** — *added from A2's own review,
  and the most consequential of the six*, because the catalog's own open compares digests only
  when the reference carries them and a `.fai`-read reference carries none. Without it a catalog
  from another build of the same assembly routes every repeat tract to the wrong position,
  genome-wide, with nothing to notice;
- **and each sample's contig checksums against the reference's** — the case names and lengths
  cannot catch.

**⛦ Two of the three checks this step was planned around turned out to be built already, and the
step changed shape twice. Both changes are the owner's, 2026-08-31.**

**⛦ First: the reference check.** It read "the analysed regions agreed across samples", which is a
psp fact — two files can record different ground, but a direct run computes one segmentation and
hands the same one to every sample, so the comparison cannot differ. The owner replaced it with the
contig agreement. **Then A2 found that check almost entirely built already, and building the rest
showed exactly what is left.** The open gate in `SampleReads::open`
([`open_bam.rs`](../../../../src/ng/read/input/open_bam.rs)) compares each file's `@SQ` list against
the reference's contig table — names, lengths, order **and the `M5` digests, whenever the reference
carries them**. A file aligned to another assembly, opened against a reference read from its FASTA,
never opens at all: the gate refuses it.

**So `check_assembly` covers one case and only one**, and it is the ordinary one: a run whose
reference was read from a `.fai`. That path hands back the contig table at once and verifies the
FASTA on a background thread, so the files open against a reference carrying *no* digests and the
gate has nothing to compare. The digests arrive only when the run joins that thread, and comparing
them then is what A2 builds. The caller therefore takes the verified reference and reports an
`AssemblyCheckOutcome` — because *no sample was aligned to a wrong assembly* and *no sample could be
checked* are different facts, and a run report has to tell them apart.

**⛦ Second: the sample-name match.** It cannot be done here. `RunParameters` carries no names — one
number per sample and one per read group, in the run's order — and a supplied file's names are
matched against the run's at that file's own door
(`ParametersFile::to_run_parameters_for`, which refuses naming the position where the two lists
diverge). What is left, and what nothing else prevents, is parameters assembled for one cohort
handed to a caller opened over another; A2 catches that on the counts.

*Depends:* A1. *Source:* §6.2, §7.1a; [`parameters_file.md`](../spec/parameters_file.md) §6;
[`../arch/run_streaming.md`](../arch/run_streaming.md) §5.

> **Checkpoint A: met 2026-08-31.** The refusals are covered by 40 tests across `callers.rs` and
> `segments.rs`. The 63-sample open is `examples/ng_open_cohort_descriptors.rs`, which opens the
> tomato accessions in `benchmarks/tomato1/crams/` — 63 files, 100,171 segments over the 80 BED
> regions, **819 of 819** contig checksums compared and agreeing — and counts what that costs in
> file descriptors. **It takes the catalog's path rather than finding it beside the reference**,
> because the reference is on the container's read-only `$HOME/genomes` mount and nothing can be
> written next to it. What it measured corrects the mechanism behind `DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS`'s
> constant without changing the constant; `callers.rs` carries the numbers and PROJECT_STATUS the
> correction against spec §7.1a.

### Milestone B — one sample's walk, behind the merge's interface

✅ **B1. The walker as an `ObservationSource`.** Adapt `SampleLocusObservationsIterator` to the
merge's trait: forward-only, one observation at a time, offering the spare record back for reuse.
One source per sample for the whole run — not one per worker, not one per segment.
*Depends:* A1. *Source:* §3.4; [`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs).

✅ **B2. One sample through the source equals one sample walked directly.** The observations a source
yields over the analysed ground are exactly those the iterator yields, in the same order. **Oracle:**
the existing iterator, driven directly.
*Depends:* B1. *Source:* §12 oracle 4.

**⛦ §12's fourth oracle is not literally true, in one field, and that is the owner's to settle.**
It asks that a segment walked alone emit *exactly* what the same span emits inside a whole walk.
Measured on the real generator, everything is equal but the **chain ids**: the type's own
documentation says "an id names a read within one walk", and the allocator counts up across a whole
walk and survives the per-chromosome reset — so a read is id 0 walked alone and id 4 walked fourth.
No implementation of a walk-scoped id can satisfy "exactly". B2 compares the **grouping** instead,
which still catches a read split in two, two reads merged, or a locus that lost its witnesses.
**The same question is owed to §12's first oracle**, byte-identical psps across worker counts,
which inherits the problem once chain ids are written to a file.

> **Checkpoint B: met 2026-08-31.** A walker is indistinguishable from any other source, proved
> against the machinery that existed before it: the real generic generator over a real indexed BAM,
> 62 loci across four generic segments, a satellite, a gap, a contig change and one analysed-but-
> empty stretch, compared whole. A 21-mutation pass on B1 killed 20; a 14-mutation pass on B2
> killed 13, six of them invisible to B1's tests. **Both survivors are the same one**: a walker
> that stashed every offered record and never freed one, which nothing can pin while there is no
> pool to count — recorded against G1. Pause for review.

### Milestone C — the merge, fed by walkers

☐ **C1. The merge driven over walker sources.** The single-threaded merge (`merge_cohort_through_cache`)
reading k walkers instead of the in-memory sources its tests use, yielding cohort loci in genome
order.
*Depends:* B2. *Source:* §3.2; [`cohort_merge.md`](../spec/cohort_merge.md).

☐ **C2. Cohort loci from reads match cohort loci from records.** Walk a small cohort, capture its
observations, feed the same observations to the merge from memory, compare. **Oracle:** the merge's
own in-memory driver, which is already the reference for its parallel one.
*Depends:* C1. *Source:* §12.

> **Checkpoint C:** alignment files to cohort loci, single-threaded, against an oracle.
> Pause for review.

### Milestone D — calling, joined to the merge

☐ **D1. `call_locus` in the builder.** The builder that assembles a cohort locus also calls it. The
spec says the placement commutes; this is the wiring
[`calling_loop.md`](calling_loop.md) lists as its own remaining work and no plan has claimed.
*Depends:* C2. *Source:* [`run_streaming.md`](../spec/run_streaming.md) §3.1;
[`calling_loop.md`](calling_loop.md), *Out of scope*.

☐ **D2. The sample-order join. Own commit, do not bundle.** The merge names samples by their index
in the run's order; the parameters name them by sample name; the scratch rows are the run's sample
order with the uncallable ones closed up. **Three numberings, and a mismatch produces wrong
genotypes rather than a crash** — the same accident that made six mutations survive the calling
loop's own tests. **Oracle:** a fixture where the three orders genuinely differ, so that swapping
any two changes a called genotype.
*Depends:* D1. *Source:* §5.1; [`calling_loop.md`](calling_loop.md) Milestone E1.

☐ **D3. The end-to-end fixture.** A handful of the tomato slices
(`benchmarks/tomato1/crams/`) over a small BED, alignments in, called genotypes out, on the generic
path. **Runs in the dev loop in minutes, not hours.** **It also reports where the time went** —
walking the reads, assembling loci, genotyping them — because Milestone E's shape is decided from
that split and this is the first run that can produce it (owner's ruling above; spec §11 question 7
asks the same question and says nobody has measured it).
*Depends:* D2. *Source:* §12 oracle 3 (its direct-mode half).

> **Checkpoint D:** ng calls genotypes from CRAM files. Pause for review — **this is the milestone
> the whole plan exists for.**

### Milestone E — the pool

**⛦ Owner's ruling, 2026-08-31: Milestones A to D build against the single-threaded merge, and what
this milestone becomes is decided from D3's measurement.** Two arrangements can genotype several
loci at once — the merge's own region batching switched on, so each thread assembles and genotypes
its own stretch of ground, or the merge left on one thread handing each finished locus to workers
that only genotype — and nothing measured says which. **The region batching may not survive either**;
whether it is worth keeping at all is the same measurement. So D3 comes first and reports where the
time goes, and this milestone is shaped after it.
*Background:* [`../arch/run_streaming.md`](../arch/run_streaming.md) §8; spec §11 question 7.

☐ **E1. Callers in flight.** The merge stays on one thread; each cohort locus goes to a free worker;
results are released in genome order. The bound is `callers in flight × one cohort locus`.
*Depends:* D3, and the ruling above. *Source:* §3.5, §5.1.

☐ **E2. Concurrency invariance. Own commit, do not bundle.** The same VCF-bound output at one caller
and at sixteen. **This is where a missed reset in `CallingScratch` shows** — the scratch is per
worker and reused across loci, the code already records that a dropped `clear()` is invisible in one
locus order, and under a pool the order is a scheduling artefact. **Oracle:** the serial caller of
D3, on a fixture whose loci differ in kind (ordinary sites and repeat tracts interleaved), run at
several worker counts.
*Depends:* E1. *Source:* §8 (the calling-scratch trap), §12 oracle 2.

> **Checkpoint E:** the answer does not depend on the worker count. Pause for review.

### Milestone F — the command

☐ **F1. `call-from-alignments`.** The subcommand: reference, catalog, alignment files, parameters
file or `--defaults`, analysed regions, output. Kebab-cased from its enum variant, like the three
that exist ([`cli.rs`](../../../../src/pop_var_caller_exp/cli.rs)).
*Depends:* E2, and the command-surface note above. *Source:* [`typed_regions_cli.md`](../spec/typed_regions_cli.md).

☐ **F2. The run writes the parameters it used, beside its output.**
*Depends:* F1. *Source:* [`parameters_file.md`](../spec/parameters_file.md) §7.

☐ **F3. The run report.** What the run refused and why: loci the merge would not build, samples with
no reads, every parameter that was defaulted rather than fitted.
*Depends:* F2. *Source:* §13 (`cohort_merge.md`'s refusal counts);
[`parameters_file.md`](../spec/parameters_file.md) §8.

> **Checkpoint F:** a person can run ng on a cohort of CRAMs from the command line. Pause for review.

### Milestone G — stop the merge freeing the records it was handed

**Last because it is an optimisation and everything above is correctness**, and because measuring it
needs a walker feeding a real merge — which is B1 and C1. Nothing here changes what a run answers.

**What it is.** The merge is the last owner of every per-sample observation record: it draws one,
walks it, evicts it, and drops it. On 63 tomato accessions over 100 kb that dropping is **25.9% of
the merge's CPU**, and it is expensive out of proportion to its size because a record is allocated by
whichever worker drew it and freed by whichever worker later evicts it — mimalloc takes a locked path
when those differ, and one atomic instruction inside `free` is 10.6% of the merge's whole CPU
([`../research/cohort_merge_parallel_cost_2026-08-28.md`](../research/cohort_merge_parallel_cost_2026-08-28.md) §2.2).

**The hook already exists and nothing fills it.** `ObservationSource::next_observation` takes the
spare record back, and the observation cache already offers one — B1's walker is the first thing that
could accept it. A walker that refills the returned record instead of allocating a new one removes
**92% of the merge's frees** (21.4 of 23.1 million on that ground, counted with dhat), and costs the
walker nothing it was not already doing: it fills a record either way.

☐ **G1. The walker refills the spare instead of allocating.** `next_observation(spare)` writes the
next observation into the returned record where its buffers are the right size, and allocates only
what it must grow. **Oracle:** B2 — the observations a source yields must not change, and a refilled
record that kept a stale field would show there.

**⛦ And it owes a test B1 could not write.** B1's suite pins that the spare does not come back out
as an observation; it cannot pin that the spare is *released*, because the dropping walker has no
pool to count. A walker that stashed every offered record for ever passed all fourteen of B1's
tests — the one survivor of a 21-mutation pass that killed the other twenty. **At 63 samples that
is unbounded growth of exactly the records this step exists to stop allocating**, so the step that
starts keeping records must bound how many it keeps and assert the bound.
*Depends:* B1, C1. *Source:* §3.4; [`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)'s `spare` list.

☐ **G2. What it was worth, on real reads.** The merge's wall time with the walker leasing against
minting, arms alternated inside one process. **Report the cohort size, the density, the machine and
the allocator**; and report it against the walk's own time, because the merge is 2.6–18% of
walk-plus-merge on the only cohort this has been measured on.
*Depends:* G1. *Source:* the research finding above, §5.5.

> **Checkpoint G:** the merge stops paying for records it did not make. Pause for review.

---

## Verification summary

| milestone | proven by |
|---|---|
| A | a 63-sample cohort opens; each construction refusal fires and names what differs, including the descriptor headroom |
| B | a source's observations equal the walk's, position for position (B2) |
| C | cohort loci from live reads equal cohort loci built from the same observations in memory — the merge's own oracle (C2) |
| D | genotypes on real tomato slices; and a fixture where the three sample numberings differ, so swapping any two changes a call (D2) |
| E | identical output at one caller and at sixteen, on a fixture mixing ordinary sites and repeat tracts (E2) |
| F | the command runs a cohort end to end and writes back the parameters it used |
| G | B2's oracle still holds with a leasing walker, and the merge's time is measured against the walk's |

## Out of scope (next plans)

- **The emission step's spec, then the VCF writer** — the next thing on the critical path, and
  writable once D3 produces called loci to describe.
- **psp mode's three stages** — their own plan, once the encoding settles.
- **Repeat tracts end to end** — [`candidate_alleles_ssr.md`](candidate_alleles_ssr.md) first.
- **Whether any of this needs more parallelism than E1 gives it** —
  [`run_streaming.md`](../spec/run_streaming.md) §11 question 7 and the merge's research plan.
