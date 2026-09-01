# ng — driving a run: types and interfaces

*Status: architecture draft (2026-08-16), **amended 2026-08-31**. Companion to the spec
[`../spec/run_streaming.md`](../spec/run_streaming.md) (the design and every *why*) and to the
shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). The merge these callers drive has its
own pair — [`../spec/cohort_merge.md`](../spec/cohort_merge.md) and
[`cohort_merge.md`](cohort_merge.md) — and is built; this document says what drives it. Naming
follows [`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for
types, verbs for functions, **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the
**contract** is the deliverable. This document does not re-argue a decision — the spec section is
cited instead.*

The public surface is the spec's three objects (spec §1): `AlignedFilesVariantCaller`,
`PspVariantCaller`, `SampleObservationGatherer` — each an iterator. Everything else in this
document is crate-private machinery inside one of them.

## What the amendment of 2026-08-31 changed

The 16 August draft ran the genome through a pool of workers, each owning one stretch. Spec §3.5
retired that arrangement and the merge was built to its replacement, so this document described
machinery that does not exist. Two later documents recorded the divergence on their own pages
instead of here — [`cohort_merge.md`](cohort_merge.md)'s note that it revises §3.2, and spec §3.2's
flag that its own text "describes the merge this document imagined, not the one that was built" —
and nobody came back to the original. **The amendment takes no design decision the spec or the
built code had not already taken.** Six things moved:

- **Where the parallelism is.** Stretches of genome are no longer dealt to a pool. One thread
  merges, and each cohort locus it produces goes to a free caller (spec §3.5). §3 is rewritten
  around that, and the `LookAhead` knob — segments in flight — is deleted. What survives is a count
  of **callers in flight**, and it belongs to the caller objects rather than to all three.
- **What a source is asked for.** `observations_in`, which returned an iterator over one segment,
  is replaced by the trait the merge built and uses: one observation at a time, with the spent
  record offered back for reuse (§2).
- **Where that trait lives.** In the merge, not in this module — recorded, not moved (§2).
- **The calling core.** `call_vars_in_segment` and `k_way_merge` are both gone: the merge is
  [`cohort_merge`](cohort_merge.md), and calling one cohort locus is the calling module's existing
  `LocusGenotyper::call_locus` reached through two shaping functions (§3.2).
- **Two refusals were missing** — the parameters' sample list against the run's, and the
  file-descriptor headroom (spec §6.2, §7.1a). Both are in `RunError` now (§5).
- **Every `file:line` in §7 was re-read.** The citations had drifted by up to 215 lines.

---

## Module home

`src/ng/run/`, beside `locus_generation/` and `parameter_estimation/`:

```
src/ng/run/
├── mod.rs         – re-exports the three objects + RunError
├── segments.rs    – Segmentation + SegmentationInputs: the run's segments and what they
│                    were computed from
├── walker.rs      – one sample's alignment files behind the merge's source trait
├── callers.rs     – AlignedFilesVariantCaller, and later PspVariantCaller
├── gatherer.rs    – SampleObservationGatherer
├── psp_header.rs  – PspHeader: the values every psp records
└── cohort_merge/  – built; its own arch doc (cohort_merge.md)
```

**`walker.rs` and `callers.rs` replace the draft's `source.rs` and `calling.rs`.** `source.rs` was
named for the trait it would hold, and the trait is in the merge (§2), so the file holds only the
alignment-backed implementation and is named for it. `calling.rs` was named for two free functions
that no longer exist; what goes in that file is the caller objects.

The psp **writer** and **reader** are deliberately not files here: their shapes belong to the
encoding spec (spec §10). This document defines what they plug into — the writer consumes the
gatherer's iterator (spec §5.2), the reader implements the merge's source trait (§2) — and the
header values they write and check (§4). The walk stage's per-sample loop (spec §5.2) is
composition the CLI owns; it needs no type of its own.

**This revises [`module_layout.md`](module_layout.md)'s `pipeline.rs` entry**, which named one
file holding "the `CallerRecipe` + the driver that runs it end-to-end". The recipe stays where
that document puts it; the drivers are the two caller objects here. Recorded rather than quietly
changed, because the tree in that document still shows the old placement.

---

## 1. The segments (`segments.rs`)

The segments every loop advances over are the typed-region generator's segments, held once and
shared read-only by every worker (spec §4.2, §9). No grouping routine exists: the loop unit for
*observation generation* is one segment (spec §4.4), so `Segmentation` is a list plus the record of
its inputs — the previous draft's `group_toward` is retired with the size choices it served.

**Built 2026-08-31**; the block below is the landed shape, and three things in it differ from
the sketch this document carried. It is neither `Clone` nor deriving `Debug` — a clone would
deep-copy the genome-sized list for a holder the design says should not exist, and a derived
`Debug` would print every segment. `build` takes the catalog's path, because four of the
catalog's own failures describe the file without naming it and a person with several catalogs
cannot act on that. And the analysed regions are held twice on purpose, in the two shapes their
two consumers take.

```rust
/// The run's segments, in genome order, plus the values they were computed from. A function of
/// the reference, the catalog, the repeat-tract criteria and the analysed regions — no
/// sample's reads — so it is identical in every sample of the run (spec §4.2).
pub struct Segmentation {
    inputs: SegmentationInputs,
    segments: Vec<TypedRegion>,
    /// What the merge takes, where `inputs.analysed_regions` is what the checks compare.
    analysed_regions: Vec<GenomeRegion>,
}

impl Segmentation {
    /// Consumes the typed-region generator's stream once. `analysed` is both recorded and
    /// used here, so what this reports and what the merge walks are one value; the catalog
    /// and the criteria are recorded as given, so reading the stream and building the record
    /// belong in one call site.
    pub fn build(
        segments: impl Iterator<Item = Result<TypedRegion, RepeatCatalogError>>,
        analysed: GenomeRegions,
        catalog: RepeatCatalogHeader,
        repeat_tract_criteria: StrRepeatCriteria,
        catalog_path: PathBuf,
    ) -> Result<Self, RunError>;

    pub fn inputs(&self) -> &SegmentationInputs;
    /// In genome order; a segment never crosses a contig and is never cut (spec §4.3).
    pub fn segments(&self) -> &[TypedRegion];
    /// The ground the merge advances over.
    pub fn analysed_regions(&self) -> &[GenomeRegion];
}
```

**What holding the segments whole costs is unmeasured, and it is the one term of a run that
grows with the genome rather than with the cohort.** How many segments a genome has at the
criteria a run asks with has never been counted — spec §11's question 1 is that measurement.
Beside 63 open alignment files it is small; at one sample it may not be.

`SegmentationInputs` is both the psp header's core (§4) and the operand of the file-against-run
check (spec §6.2). The parameters fit holds the equivalent object for its own compatibility
check: `RecordingTerms`
([`census.rs:1024`](../../../../src/ng/parameter_estimation/joint/census.rs)).

```rust
#[derive(Clone, PartialEq)]
pub struct SegmentationInputs {
    /// The catalog file's own header, **reused whole rather than restated**: it already carries
    /// the whole-reference MD5, the criteria the catalog was built under, the scan weights and
    /// the tool version (`repeat_catalog/mod.rs:283`).
    pub catalog: RepeatCatalogHeader,
    /// The criteria the *reader* asked with. Not the same value as `catalog.built_under` — the
    /// catalog is built below every floor a reader might ask with, so a reader filters rather
    /// than re-scans — and it is this one that decides where a segment ends.
    pub repeat_tract_criteria: StrRepeatCriteria,
    /// The regions the run was asked to analyse (`region_typing/mod.rs:77`). The field a user
    /// actually changes between runs, and the one compared across the cohort (spec §6.2).
    pub analysed_regions: GenomeRegions,
}

impl SegmentationInputs {
    /// The name of the first field that differs, for the error message; `None` when they
    /// agree. **A name rather than a `bool`**: "these two segmentations differ" leaves the
    /// user nothing to fix (spec §6.1). The names are noun phrases that read inside §5's
    /// "written under a different {field}" — *repeat catalog*, *set of repeat-tract
    /// criteria*, *set of analysed regions* — never field identifiers.
    ///
    /// **The order is the order a person should fix them in**, and it is checked: the catalog
    /// carries the reference's identity, so under a different catalog the other two
    /// comparisons are about different genomes.
    pub fn first_difference(&self, other: &Self) -> Option<&'static str>;
}
```

`PartialEq` and not `Eq`, because `StrRepeatCriteria`
([`repeat_catalog/criteria.rs:61`](../../../../src/ng/repeat_catalog/criteria.rs)) wraps
`SsrSegmentCriteria`, whose `min_purity` is an `f32`
([`segment_criteria.rs:502`](../../../../src/ng/region_typing/segment_criteria.rs)).

---

## 2. The source (`walker.rs`, and the trait in `cohort_merge/`)

**A source answers one question: what did this sample see next?** It hands back one observation at
a time, in coordinate order, never going backwards, for the whole run. One trait carries the entire
difference between the two callers (spec §3.3): behind it sits either a walk over alignment files
or a read from a psp, and nothing above it can tell which.

**The trait is already in the tree, and it is in the merge rather than here** —
[`ObservationSource`](../../../../src/ng/run/cohort_merge/observation_cache.rs), declared at
`observation_cache.rs:70`. The merge needed it before this module existed and defined it locally;
its own arch doc ([`cohort_merge.md`](cohort_merge.md) §2) says this document owns it. **Decided:
it stays where it is and this document describes it.** Moving a trait the whole merge already
implements against buys nothing a reader of either document cannot get from a cross-reference.

```rust
/// One sample's observations, one at a time, in coordinate order, forward only.
pub trait ObservationSource {
    /// What a failed read is. The merge adds nothing to it and passes it through, so it must
    /// name the sample it came from (§5).
    type Error;

    /// The next observation, or `None` once this sample is spent.
    ///
    /// `spare` is a record the merge will not read again, offered for reuse; `None` when it has
    /// none to hand back.
    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, Self::Error>>;
}
```

**Contract.**

- **Coordinate order, forward only.** Observations start at non-decreasing positions. Going
  backwards is not an error the source may report: it trips an assertion inside the merge's cache
  (`ObservationCache::cover`), because the cache has no error of its own to mint. That is right
  while observations come from this crate's generators, where backwards is a bug; `organise.rs`
  records that it owes the change to a returned `RunError` once observations are decoded from a
  file.
- **Exhaustion is final.** Once `next_observation` answers `None` it is never called again — the
  cache holds a flag of its own to guarantee it. A source that yielded `Some` after a `None` would
  be drawn in behind the cache window's right edge, and so silently out of coordinate order.
- **A failure leaves the source live.** `Some(Err(_))` ends the merge but does not end the source;
  it is `None` that ends it.
- **One source per sample for the whole run**, not one per worker and not one per segment (spec
  §3.4; [`cohort_merge.md`](cohort_merge.md) §6.4's "one reader per sample for the whole run"). The
  merge is the only consumer and it only moves forward, so each stretch is decoded once and the
  backward jump never happens. Both sources are still *capable* of a backward jump — asked for
  ground behind them they seek, and pay for it (spec §8).
- **Reuse is optional.** Every iterator of one sample's observations is already a source, through a
  blanket implementation (`observation_cache.rs:98`) that drops the spare and calls `next`. That is
  what lets fixtures and probes hand the merge a plain `Vec`'s iterator unchanged. A source opts
  into reuse by implementing the trait itself.
- **`Send + Sync` only for the parallel merge.** `merge_cohort_in_parallel`
  ([`parallel.rs:96`](../../../../src/ng/run/cohort_merge/parallel.rs)) shares the cache across
  rayon workers and so requires `S: Sync + Send`; the single-threaded driver does not. A walker
  built on `Rc` or `RefCell` would therefore work in direct mode and not under the parallel merge,
  which is off by default (spec §3.5).

**Two implementations:**

- **The walker** — **built 2026-08-31 as `AlignmentFilesWalker`**
  ([`run/walker.rs`](../../../../src/ng/run/walker.rs)). It is
  `SampleLocusObservationsIterator` ([`locus_generation/mod.rs:921`](../../../../src/ng/locus_generation/mod.rs))
  driven over the run's segments and made to answer the trait, and it is **a wrapper rather than
  a re-implementation**: everything the ownership shape this paragraph described — the
  `SampleCursor` (`SampleReads::cursor` takes `&self` and returns an owned, `Send` cursor —
  [`read/input/mod.rs:623`](../../../../src/ng/read/input/mod.rs), test
  [`:1441`](../../../../src/ng/read/input/mod.rs)), the reference accessor from the factory
  (`WindowedRefSeq` is `Send` and deliberately not `Sync` —
  [`read/input/mod.rs:606-611`](../../../../src/ng/read/input/mod.rs)), and the generator set whose
  drop order is load-bearing
  ([`locus_generation/mod.rs:1028-1040`](../../../../src/ng/locus_generation/mod.rs)) — is held by
  that iterator, so the walker inherits it including the `Drop` impl that guards the order. Spec
  §8 is the trap list this shape honours.

  **What the wrapper adds is the two things the trait asks for that a plain iterator cannot
  give.** A failure that names the sample and how far the walk had got (`RunError::SourceFailed`,
  §5), and the spare-record hook — taken, and **dropped**, which the trait's own contract permits;
  filling it is the plan's step G1 and this is what gives G1 something to change, since the
  blanket implementation is not editable for one source.

  **⚑ It does not honour one clause of the contract above, and the deviation is recorded rather
  than fixed.** *A failure leaves the source live* is what lets a cover be made again; the wrapped
  iterator latches `done` on an error, so a failed walk is **spent** and a consumer that swallowed
  the error and asked again would be told the sample is exhausted. Nothing does that today —
  `ObservationCache::draw_next` propagates without marking the source spent, and both drivers
  abandon the cache — so it is unreachable. It is written down because the failure it would cause
  is silent: cohort loci built without that sample, wrong genotypes rather than an error. **Any
  change that adds a retry has to fix this first.**

  **It is neither `Send` nor `Sync`, and not by its own choice.** `GeneratorSlot::Generator` holds
  a `Box<dyn LocusGenerator<S>>` with no auto-trait bound — deliberate, and recorded at that type.
  So a walker cannot go under `merge_cohort_in_parallel` without widening that trait object, and a
  walker and the merge drawing from it stay on one thread. E1 has to plan around it.

  **The walker is deliberately not an `Iterator`.** It could not be: the blanket implementation
  at [`observation_cache.rs:98`](../../../../src/ng/run/cohort_merge/observation_cache.rs) already
  makes every iterator of observations a source, so a type that was both would implement the trait
  twice and Rust refuses the overlap. The walk stays reachable as an iterator one level down,
  which is what the observations-equal-the-walk oracle (spec §12, plan step B2) drives.

  **Its region stream is `RunSegments`, an `Arc<Segmentation>` and an index**, whose item is
  `Result<TypedRegion, Infallible>`: a run's segments were read out of the catalog once, at
  `Segmentation::build`, so this stream has nothing left to fail at and says so in its type.
  `locus_generation` carries one `From<Infallible>` impl to admit it.

  **⛦ Shared ownership, and it was a borrow for a day** (owner's ruling, 2026-08-31). B1 shipped
  `RunSegments<'a>` holding a slice iterator, and B1's own review found what that costs: a run
  holds one walker per sample *and* the segmentation those walkers read, so with a borrow the run
  would be a struct whose walkers borrow its own field — self-referential, and safe Rust cannot
  express it. Cloning is not the escape, because `Segmentation` is deliberately not `Clone`, and
  neither is minting a walker per draw, which breaks the one-source-per-sample clause of §2's
  contract. So `AlignedFilesVariantCaller` holds `Arc<Segmentation>` and hands each walker a
  reference count; the genome-sized list is still stored once, and the walker type carries no
  lifetime, which is what lets a run store it. **`Arc` and not `Rc`**, even though a walker is
  `!Send` today: what makes it `!Send` is the generator set's unbounded trait object one layer
  down, and an `Rc` would add a second blocker to find and remove if that one is lifted.
- **The psp reader** (`src/ng/psp/`, built): a cursor over one open psp that decodes whichever
  blocks it needs and keeps the one it is in. Its resident state is the file's coarse index plus one
  decoded block; measured, that is **123 kB a cursor**, on top of **357 kB** for the open file
  itself on a human reference — almost all of the second being the reference's contig list (spec
  §7.2). No block is visible in its interface.

**The cheap question is not on this trait, and spec §10 owns that gap.** Deciding whether a
position is worth calling needs two small numbers per sample, and everything else only where the
cohort kept something (spec §3.3). This trait has one method and it returns a whole inflated
observation, so a psp reader behind it would decode everything everywhere. The shape of the fix is a
second method; it cannot be settled before the encoding is, and it costs a walker nothing either
way, so direct mode is unaffected.

---

## 3. The three objects (`callers.rs`, `gatherer.rs`)

### 3.1 The shape, and where the threads are

**The two stages parallelise along different axes** (spec §3.5), and the draft's single shared
skeleton — segments dealt to workers, several in flight — describes neither.

- **The two variant callers: a serial merge feeding a pool of callers.** One thread runs the merge
  and produces cohort loci in genome order; each locus goes to a free caller as it appears; results
  are released in genome order. The genome is not cut for calling, because a cohort locus is not a
  position — a deletion joins consecutive positions into one — so where loci begin and end is an
  output of the merge, not an input to it (spec §3.5).
- **The walk: one worker per sample.** Several samples are walked at once, one worker each, and
  each sample's walk is serial inside it (spec §5.2).

The knobs are newtypes over `NonZeroUsize` — a count whose zero is illegal:

```rust
/// How many cohort loci are being genotyped at the same time — one per worker, while the
/// merge that produced them runs on its own thread. The two caller objects' one concurrency
/// knob, and what their in-flight memory is a multiple of (spec §3.5, §7.1). The name is the
/// spec's *callers in flight*; what it counts is loci, not variant callers.
/// **No default is proposed** — spec §11 question 2 names the sweep that sets it.
pub struct CallersInFlight(pub NonZeroUsize);

/// Samples being walked at once in psp mode's walk stage. Each costs one open alignment file at
/// 11–15 MiB plus its census accumulator at about 6 MB per read group, so the read-group count
/// is part of the sizing and not a formality (spec §5.2). **No default is proposed** — spec §11
/// question 2 again.
pub struct SamplesInFlight(pub NonZeroUsize);

/// Segments of one sample's own walk in flight inside one gatherer. **Not on the default path**:
/// at one, a gatherer is the serial walk spec §5.2 describes. Whether raising it scales at all
/// is spec §11 question 3, and question 8 turns on the answer.
pub struct Workers(pub NonZeroUsize);
```

**`LookAhead` is deleted.** It counted segments in flight beyond the next to yield, which was the
memory knob of the pool spec §3.5 retired. Nothing in the built design has a segment in flight.

### 3.2 Calling one cohort locus

**There is no `call_vars_in_segment` and no `k_way_merge`.** The merge is the module
[`cohort_merge.md`](cohort_merge.md) documents, and calling one locus is three existing calls in the
calling module — the caller objects compose them and add nothing of their own:

```rust
// which alleles the locus is called over, and what each covering sample lost to the cut
let selection = select_generic(&observation, &selection_config, &mut selection_scratch);
// the merge's covering samples become one entry per sample of the run
let evidence = shape_generic_locus(
    &mut shaping, &observation, &selection, run_sample_count, &mut views,
);
// the call itself; the allele table leaves the selection by value (`LocusSelection::into_parts`)
let inference = genotyper.call_locus(&evidence, &frozen, alleles, &loop_config, &mut scratch);
```

- `select_generic` — [`allele_candidates/generic.rs:81`](../../../../src/ng/calling/allele_candidates/generic.rs).
  The repeat-tract path's equivalent is specified and unbuilt
  ([`../spec/candidate_alleles_ssr.md`](../spec/candidate_alleles_ssr.md)).
- `shape_generic_locus` / `shape_ssr_locus` —
  [`evidence_shaping.rs:403,443`](../../../../src/ng/calling/evidence_shaping.rs). **This is where
  the three sample numberings meet**, and that module's own header documents them: the merge's
  per-sample list holds only the samples that covered the locus, each naming its index in the run's
  order; the loop's list is one entry per sample of the run; and candidate selection's dropped-allele
  list is parallel to the merge's, not to the run's.
- `LocusGenotyper::call_locus` —
  [`inference/mod.rs:619,629`](../../../../src/ng/calling/inference/mod.rs). One cohort locus in,
  one `LocusInference` out.

**Where the call happens relative to the merge's builder is open, and the two documents lean
different ways.** [`../spec/cohort_merge.md`](../spec/cohort_merge.md) §6.3 leans to calling *inside*
the builder, because the buffer then holds called records rather than whole observations, and notes
that the choice commutes — calling one cohort locus reads nothing outside it. Spec §3.5 keeps the
merge on one thread and hands each finished locus *out* to a worker. The two agree on the answer and
not on which thread computes it. §8 records what is unsettled; a run that genotypes one locus at a
time calls in the same place either way, so nothing before that depends on it.

**The comparison that would settle it is now available.** `cohort_merge.md` §14 question 6 asked to
weigh the two buffers when the emission step fixed the record's shape; the record is
`VcfRecord` ([`vcf/mod.rs:84`](../../../../src/ng/vcf/mod.rs)) and the step is specified
([`../spec/vcf_output.md`](../spec/vcf_output.md)).

**Both built merge drivers accumulate, and a real run must not.**
`merge_cohort_through_cache` ([`serial.rs:146`](../../../../src/ng/run/cohort_merge/serial.rs)) and
`merge_cohort_in_parallel` ([`parallel.rs:96`](../../../../src/ng/run/cohort_merge/parallel.rs))
both return a `RegionOutcome` holding every cohort observation of the whole run in one `Vec`
([`build.rs:743`](../../../../src/ng/run/cohort_merge/build.rs)). That is what an oracle wants and
what a run cannot afford. The caller objects consume each locus where it is built and keep none, so
what they hold is the pool's loci and not the genome's.

### 3.3 `SampleObservationGatherer` — psp mode's walk

```rust
/// One sample's observations in genome order, census accumulated as they pass. psp mode's
/// walk stage is a loop of these, one per sample (spec §5.2).
pub struct SampleObservationGatherer { /* sources, CensusWriter, tallies */ }

impl SampleObservationGatherer {
    pub fn new(
        sample: &SampleInput,
        segmentation: &Segmentation,
        census: &CensusConfig,
        workers: Workers,
    ) -> Result<Self, RunError>;

    /// After the iterator is exhausted: the census the walk accumulated. Calling it earlier is
    /// an error — the census's per-stratum sums are complete only at the end (spec §5.2).
    pub fn finish(self) -> Result<SampleCensusEvidence, RunError>;
}

impl Iterator for SampleObservationGatherer {
    type Item = Result<SampleLocusObservations, RunError>;
}
```

**Contract.** Yields in genome order. At the yield point — single-threaded by construction — every
observation passes `CensusWriter::add_locus` and every completed segment passes `mark_walked`, empty
segments included
([`census.rs:2087,2106`](../../../../src/ng/parameter_estimation/joint/census.rs)); what the
iterator yields and what the census counted are the same stream, which is the whole of spec §5.2's
two closures. Read-filter tallies are per-cursor
([`read/input/mod.rs:620-622`](../../../../src/ng/read/input/mod.rs)) and summed at `finish` —
unsummed, drop rates under-report by the worker count (spec §8).

**A gatherer never sees more than one sample's files**, whichever way the stage is parallelised: the
pool is outside it, over samples; `workers` is inside it, over that sample's segments. At the
default of one there is no inside pool at all.

### 3.4 The two variant callers

**Direct mode's construction is built (2026-08-31) and the block below is the landed shape.**
Two things in it differ from the sketch this document carried, and both remove something. There
is **no `SampleInput` type**: a run builds one read-group table from its file paths and opens one
`SampleReads` per entry of `ReadGroups::read_groups_per_sample`, which is the rule
`SampleReads::open_only_sample` states for every tool that is not single-sample — and it is what
fixes the run's sample order in one place rather than two. And the constructor takes **no
concurrency knob**, under the ruling recorded in §8.

```rust
/// What a run opens and how it filters, grouped because the four travel together into every
/// sample's open. The read-group table's order is the run's sample order.
pub struct AlignmentInputs<'a> {
    pub read_groups: &'a ReadGroups,
    pub reference: &'a OpenReference,
    pub read_filters: ReadFilterConfig,
    /// **The five knobs the locus generator walks with** (2026-09-01) — the two per-column
    /// depth caps, the widest record footprint, the mate-lookup window and the ceiling on
    /// reads held open at once. Beside `read_filters` because the two answer one question:
    /// how this run turns bytes into evidence. **Checked at `open`**, with the other
    /// refusals, so `RunError::LocusGeneratorSettings` fires before a file is opened rather
    /// than at the first locus.
    pub locus_generator_settings: PileupGeneratorConfig,
    pub build_index_if_missing: bool,
    /// The reference **once its per-contig checksums are known** — what each sample's own
    /// contig checksums are compared against (§5's `SampleAlignedToAnotherReference`). **Not
    /// read at a sample's open**: the checksums to compare are captured as each file opens, so
    /// this is used after they all have — as are the generator settings above, which are read
    /// where the walkers are built. Only the caller can supply it: a
    /// `.fai`-only reference and one whose FASTA has not been read are the same value, so the
    /// run reports which it got rather than inferring it.
    pub reference_with_checksums: &'a ReferenceInfo,
}

/// What a run could learn about the assembly its samples were aligned to. **Two different
/// facts, and only one is reassuring**: *every sample agrees* is a check that ran, where
/// *nothing could be compared* is a check that did not, and a run report must tell them apart.
pub enum AssemblyCheckOutcome {
    /// The second number is what makes the first mean something: "1,386 of 1,512".
    EverySampleMatchedTheReference {
        alignment_files: usize,
        checksums_compared: usize,
        checksums_possible: usize,
    },
    /// Not one checksum could be compared. `because` names the side that had none — the
    /// reference (a `.fai`, an unread FASTA, or a trusted index) or the alignment files
    /// (ordinary: `@SQ M5` is optional).
    NothingCouldBeChecked { because: NoChecksums },
}

/// Direct mode (spec §5.1). Holds every sample's SampleReads open for the whole run —
/// 11–15 MiB **per open alignment file** — plus the shared read-only state, and from B1 one
/// walker per sample advancing at the merge frontier.
pub struct AlignedFilesVariantCaller { /* SampleReads per sample, read groups, reference,
                                          filters, Arc<Segmentation>, params, configs */ }
// The segmentation is behind an Arc so the walkers this run will own can hold it too; a
// walker that borrowed it could not be stored beside it (§2, and `shared_segmentation`).

impl AlignedFilesVariantCaller {
    /// Opens every sample's files — and, from A2, runs spec §6.2's and §7.1a's checks —
    /// before a read is decoded. Named `open` because that is what it does.
    ///
    /// **A fourth refusal came with the merge (2026-08-31)**: the reference is opened for
    /// walking here too, and one that holds no bases is refused before a file opens.
    pub fn open(
        alignments: AlignmentInputs<'_>,
        segmentation: Segmentation,
        parameters: RunParameters,
        calling_loop_config: RunnableCallingLoopConfig,
        candidate_selection: CandidateSelectionConfig,
        merge_parameters: MergeParameters,
    ) -> Result<Self, RunError>;
}

/// psp mode's calling stage (spec §5.3). `open` reads every header and runs both checks of
/// spec §6.2 before any block is decoded; the analysed regions come from the headers, and the
/// segmentation is rebuilt from the catalog and routing criteria the run was handed.
pub struct PspVariantCaller { /* open psps, segmentation, params, pool */ }

impl PspVariantCaller {
    pub fn open(
        psps: &[PathBuf],
        catalog: RepeatCatalogHeader,   // plus the open catalog + reference the rebuild needs
        routing: StrRepeatCriteria,
        parameters: RunParameters,
        loop_config: RunnableCallingLoopConfig,
        selection: CandidateSelectionConfig,
        merge: MergeParameters,
        callers_in_flight: CallersInFlight,
    ) -> Result<Self, RunError>;
}

impl Iterator for AlignedFilesVariantCaller { type Item = Result<VcfRecord, RunError>; }
impl Iterator for PspVariantCaller        { type Item = Result<VcfRecord, RunError>; }
```

**`MergeParameters` is a grouping, not a new decision.** The merge's five run parameters already
exist and are separate types
([`cohort_merge/mod.rs:269,314,367,465,532`](../../../../src/ng/run/cohort_merge/mod.rs)); passing
them one by one alongside four other arguments is what the grouping avoids. Whether it is worth a
struct is the constructor's to settle when it is coded.

**⛦ Direct mode calls its cohort, and calling happens in the builder (2026-09-01).**
`AlignedFilesVariantCaller::call_cohort(&genotyper)` drives the merge over one walker per sample
and genotypes each cohort locus **where it is built**, through the three calls §3.2 sketches. It
returns a `CalledCohort`: the called loci in genome order, the ground of the loci the width bound
refused, and what each sample's walk counted.

**How the call got inside the builder without the merge learning about calling.** `build_region`'s
locus walk moved into `build_region_handing_over`, which hands each surviving locus to a sink and
each refused locus's span to a vector; `build_region` is that function with `Vec::push` for a sink,
and `merge_cohort_through_cache` split the same way into `merge_cohort_handing_each_locus_over`.
The run supplies the calling sink, so `cohort_merge` imports nothing from `ng::calling`, the
ownership rule of spec §6.1 is still written once, and every existing oracle of the merge checks
both drivers at once. `merge_cohort` stays as the merge's oracle rather than the run's path.

**⛦ It still returns everything at once**, which is not spec §5.1's bound —
`callers in flight × one cohort locus` plus the frontier. What it *does* bound is the
observations: each is dropped as soon as its genotypes exist, which is what calling in the builder
buys. The refused-span list accumulates for the whole run as well. The pool milestone is where the
calls start being released singly, and it inherits a driver that no longer has to buffer the
observations to get there.

**⛦ What the run keeps from its walk (2026-09-01).** `ObservationCache::into_sources` hands the
per-sample readers back, so each walker's region accounting, its SNP/indel generator's counters and
the assembly-check outcome reach `CohortWalkTallies` before the walkers are dropped. **The
per-read-group read-filter tallies are still unreachable, and not for want of an accessor**: at
each contig boundary the retiring cursor's read-group counts are dropped rather than accumulated,
so by the end of a walk every contig but the last has lost them. Spec §8 requires a run to sum them
at the end; reaching them is a change to the locus generator, and it is F3's.

**Contract, both callers.** Records in genome order, identical at every number of callers in flight
(spec §12.2), and identical between the two callers on one cohort with fixed parameters (spec §12.3
— the regression anchor). Iteration ends at the first `Err`; direct mode leaves nothing to clean up,
and a psp without a valid trailer is refused at `open` (spec §9).

**What direct mode holds for the whole run:** one open `SampleReads` per sample, one walker per
sample, the segmentation, the parameters, the merge's observation cache, and the pool's loci. **The
open files are the whole bill** — 0.9 GB at 63 samples and 15 GB at a thousand, before a read is
decoded (spec §5.1) — which is what puts the descriptor check (spec §7.1a) and "where does direct
mode stop being usable" (spec §11, question 6) on the same axis.

---

## 4. What a psp header carries (`psp_header.rs`)

The bytes are the encoding spec's decision (spec §10); this file defines only the values.

```rust
/// The fields every psp header records (spec §6.1).
pub struct PspHeader {
    /// The values the segmentation was computed from (§1). `analysed` is the cohort check;
    /// catalog and routing are the file-against-run check; a refusal names the first
    /// differing field (spec §6.2).
    pub segmentation_inputs: SegmentationInputs,
    /// When the file was written. **Never compared** — the one field spec §12.1's
    /// byte-identity oracle excludes, skipped on purpose, not by omission.
    pub written_at: SystemTime,
}
```

**No `BlockBoundariesDigest`, no `writer_version` — dropped, spec §6.3, flagged there for the
owner.** No code path looks at two samples' blocks together, so boundary equality has no
consumer; format versioning belongs to the encoding spec. The one restriction the encoding spec
inherits: block cuts are a function of the observation stream alone, so the file is byte-identical
across worker counts (spec §6.3, §12.1).

---

## 5. Errors (`mod.rs`)

**The reason comes from the cause, not from the top line.** Every variant here names the sample
or the file the trouble is in; what is wrong with it arrives through the wrapped cause, so a
command reporting one of these renders the whole chain with `format_error_chain`
([`src/error_render.rs`](../../../../src/error_render.rs)) and never `Display` alone. A bare
`Display` says which sample would not open; the chain says its index is missing and where it was
looked for. **Two of the variants below are built (2026-08-31); the rest arrive with A2 and psp
mode.**

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    /// **Built.** The run's segments could not be read out of the repeat catalog. The path is
    /// carried here because the catalog's own error often does not have it: a digest mismatch,
    /// over-permissive criteria, differing scan weights and a differing tool version all
    /// describe the file without naming it.
    #[error("the run's segments could not be read from the repeat catalog {}", path.display())]
    Catalog { path: PathBuf, #[source] source: RepeatCatalogError },
    /// **Built.** One sample's alignment files could not be opened. The sample's name is a
    /// field of its own because the wrapped error does not carry it — an open failure knows
    /// which file it was, not which individual the file holds. Boxed to keep this type small.
    #[error("sample {sample}: its alignment files could not be opened")]
    OpeningSample { sample: String, source: Box<IngestError> },
    /// **Built 2026-08-31.** One sample's source failed. **Both the sample and where it had
    /// reached** — neither alone locates a failure in a run over thousands of samples (spec §9).
    ///
    /// **`reached` is an enum and not a bare position**, which is the one thing the landed
    /// shape changes: a source that fails on its very first draw has no position, and a run
    /// that said "failed at contig 0:1" when nothing had been read would send an operator to an
    /// innocent locus.
    ///
    /// **⚑ And this line carries no instruction, unlike the four refusals that do.** A cause is
    /// appended after a colon (`format_error_chain`), so an instruction here would land in the
    /// middle of the sentence, ahead of the thing it tells the reader to act on. Every variant
    /// that ends with one has no cause beneath it. Advice on a source failure belongs to
    /// whatever reports the run.
    #[error("sample {sample}: reading its observations failed; {reached}")]
    SourceFailed {
        sample: String,
        // NothingYet renders "it had produced no observations yet"; After(GenomePosition)
        // renders "its last complete observation ended at contig N position P". Both name
        // their own role, because the cause carries a SECOND genome coordinate — the region
        // that failed — and a bare position beside it reads as one fact said twice.
        reached: WalkProgress,
        source: Box<dyn std::error::Error + Send + Sync>,
    },
    /// **Built.** The run was given no alignment files, so there is no cohort to call.
    /// Refused rather than answered with an empty output, because assembling the parameters
    /// for a cohort of none panics inside the pre-pass — this refusal comes before that.
    #[error(
        "this run was given no alignment files, so it has no samples to call; \
         check the paths or the pattern it was given"
    )]
    NoAlignmentFiles,
    /// **Built 2026-08-31.** The run's reference was opened from a `.fai` alone, so it holds no
    /// bases and no locus can carry a reference allele. **Refused before a single alignment file
    /// is opened**, and it refuses nothing a real run does: every arm of
    /// `read_reference_verifying_or_creating_fai` keeps the FASTA's path beside the geometry it
    /// read from the index.
    #[error("this run's reference was opened from a `.fai` index alone, which holds no bases: …")]
    ReferenceHasNoBases,
    /// **Built 2026-08-31.** The `<reference>.fai` beside an otherwise-right FASTA is missing or
    /// damaged — `samtools faidx` rebuilds it. Named apart from the one above because the two ask
    /// for different things: a different argument, against a repair of the one already given.
    #[error("the index beside this run's reference {} could not be read", reference.display())]
    ReferenceIndexUnreadable { reference: PathBuf, #[source] source: std::io::Error },
    /// **Built 2026-08-31, and unreachable today**: nothing yet lets a run choose its locus
    /// generator's settings, so every run builds with the shipped defaults. It exists because the
    /// settings are an argument, and an argument nobody can pass today is one somebody passes
    /// tomorrow — spec §11 and the plan's Checkpoint C ask whether they should be run parameters.
    #[error("this run's locus generator would not accept its settings")]
    LocusGeneratorSettings { #[source] source: PileupGeneratorConfigError },
    /// **Built.** The parameters were not assembled for this cohort.
    ///
    /// **An arity check and not a match by name.** The assembled parameters carry no sample
    /// names — one number per sample and one per read group, in the run's order — so the run
    /// cannot compare names even in principle. A supplied file's names *are* matched against
    /// the run's, at that file's own door: `ParametersFile::to_run_parameters_for` refuses
    /// naming the position where the two lists diverge ([`parameters_file.md`](../spec/parameters_file.md)
    /// §6). What is left for a run to catch is parameters assembled for one cohort handed to a
    /// caller opened over another, which nothing else prevents.
    #[error(
        "the parameters were assembled for a different cohort: {counted} is \
         {in_the_parameters} in the parameters and {in_the_run} in this run; re-run the \
         parameter pre-pass for this cohort, or point the run at the file assembled for it"
    )]
    ParametersAreForAnotherCohort {
        counted: &'static str,   // "the number of samples", "the number of read groups …"
        in_the_parameters: usize,
        in_the_run: usize,
    },
    /// **Built.** The run needs more open files than the process is allowed. **Names both
    /// numbers and what to do**, because raising the limit is the operator's. The count is of
    /// *files*, not samples: a sample sequenced across four lanes is four files and eight
    /// descriptors (spec §7.1a).
    #[error(
        "this run needs {needed} open files for its {samples} samples and this process may \
         open {limit}; raise the limit (`ulimit -n`) or call fewer samples at once"
    )]
    NotEnoughFileDescriptors { samples: usize, needed: u64, limit: u64 },
    /// **Built.** The repeat catalog was built on a different reference from the one this run
    /// calls against. **Checked here because the catalog's own check cannot do it on the
    /// ordinary path**: opening a catalog compares digests only when the reference it was
    /// handed carries them, and one read from a `.fai` carries none, so there a catalog is
    /// admitted on contig names, lengths and order alone. The digests exist once the FASTA has
    /// been read. Silent and genome-wide otherwise — the catalog's coordinates are where the
    /// repeat tracts are, and every segment is drawn from it.
    #[error("the repeat catalog was built on a different reference: …")]
    CatalogIsForAnotherReference { reference: PathBuf, in_the_catalog: String, in_the_run: String },
    /// **Built.** The reference whose checksums the samples were checked against is not the one
    /// their files were opened against — a caller mistake, refused rather than trusted, because
    /// the comparison downstream walks the two contig lists in step.
    #[error("the reference the samples were checked against is not the one they were opened against: {difference}")]
    ReferenceCheckedAgainstAnotherGenome { difference: String },
    /// **Built.** A sample's reads are against a different build of the reference's assembly —
    /// every contig the right name and the right length, and different bases.
    ///
    /// **The open gate catches this itself whenever the reference carries digests**, so what
    /// remains for the run is the `.fai` path: there the contig table arrives at once and the
    /// digests only when the background verification is joined, so the files open against a
    /// reference with nothing to compare. That is why the caller takes the *verified* reference
    /// and reports an `AssemblyCheckOutcome` — *no sample was aligned to a wrong assembly* and
    /// *no sample could be checked* are different facts.
    #[error("sample {sample} was not aligned to this run's reference {}", reference.display())]
    SampleAlignedToAnotherReference {
        sample: String,
        /// Named because two inputs are in play and only one is wrong — the same reason
        /// `Catalog` carries its path. Which file and which contig come from the cause, so
        /// the two sentences say different things rather than one saying the other twice.
        reference: PathBuf,
        #[source] source: AssemblyMismatch,
    },
    /// Two samples were analysed over different segments, so they are not comparable — the
    /// cohort refusal (spec §6.2). **Reachable in psp mode only**, where each file records the
    /// ground it was written over. Spec §11's question 5 may later soften this to
    /// intersection-calling; until then it refuses.
    #[error("samples {left} and {right} were analysed over different segments")]
    AnalysedRegionsDiffer { left: String, right: String },
    /// One psp's recorded catalog or routing criteria differ from the run's own, so the
    /// segments this run loops over are not the segments the file's observations were minted
    /// inside — the file-against-run refusal (spec §6.2). `field` is
    /// `SegmentationInputs::first_difference`'s answer.
    #[error("psp for sample {sample} was written under a different {field}")]
    SegmentationInputsDiffer { sample: String, field: &'static str },
    /// Two files name the same sample: a duplicated argument, or a cohort that would call one
    /// individual twice and weight the allele frequencies by it (spec §6.2).
    #[error("sample {sample} appears twice: {first} and {second}")]
    SampleAppearsTwice { sample: String, first: PathBuf, second: PathBuf },
    /// A psp ended without a valid trailer: an interrupted walk, not a short sample (spec §9).
    #[error("psp for sample {sample} is incomplete")]
    IncompletePsp { sample: String },
    #[error("i/o while reading or writing a run's files")]
    Io(#[from] std::io::Error),
}
```

**`SourceFailed` replaces the draft's `WorkerFailed`, which named a segment.** The segment was the
unit of the pool §3.1 retired; under a serial merge what a caller can say is which sample failed and
where it had reached.

**The merge's `RunEndedShort` folds in here** when this type lands, together with its
`ObservationExceedsReachCeiling` — [`cohort_merge.md`](cohort_merge.md) §8 owes that, and three
other cleanups gated on the same moment. **They are not this document's to schedule**: they touch
the merge, and §8 below records them as still owed.

---

## 6. Design decisions — decided

- **Three public objects, each an iterator; everything else crate-private.** No caller of this
  module names a work unit, a block, or a range — spec §1, §3.
- **The genome is not cut for calling.** The merge runs on one thread and its loci go to a pool;
  no line is drawn through the genome in advance, because loci are what the merge produces —
  spec §3.5. This retires the draft's segment pool, `LookAhead`, and per-segment per-sample
  walkers.
- **The loop unit for *observation generation* is still one segment**, and no grouping routine
  exists — spec §4.4. The two statements are about different stages and do not conflict.
- **One source trait, two implementations; the merge consumes the trait and never invokes a
  walk.** The whole of "one calling function, whichever mode" — spec §3.1, §3.3.
- **The source trait stays in `cohort_merge/observation_cache.rs`**, where the merge built it,
  and this document describes it — §2. Moving it would change every implementation site to
  document ownership a cross-reference already carries.
- **One source per sample for the whole run**, in both callers — spec §3.4;
  [`cohort_merge.md`](cohort_merge.md) §6.4. This is what keeps every cursor forward-only.
- **Calling composes three existing calls and mints nothing** — selection, shaping, the call
  (§3.2). *Where* they run relative to the merge's builder is open (§8).
- **The callers' item is `VcfRecord`** ([`vcf/mod.rs:84`](../../../../src/ng/vcf/mod.rs)) — the
  emission step's document exists and named it, closing the draft's `OPEN: Variant's shape`.
- **What flows from a source is `SampleLocusObservations`, bare.** No span-carrying wrapper:
  walked-empty ground is recorded by the gatherer's `mark_walked` and readable back through the
  header plus trailer, so an empty container had no remaining job — spec §5.2, §8.
- **The census accumulator lives inside the gatherer, fed at the yield point.** One stream, two
  consumers impossible to desynchronise — spec §5.2.
- **The psp reader serves segments; blocks exist only inside the writer and reader** — spec §3.3;
  their sizing is wholly the encoding spec's — spec §6.3, §10.
- **The header carries no boundary digest and no writer version** — spec §6.3, **flagged for the
  owner as a reversal of an earlier draft**.
- **Eight refusal variants, eight axes**, and each one compares a different pair of things.
  `NoAlignmentFiles` asks whether the run has a cohort at all; `ParametersAreForAnotherCohort`
  compares the parameters to the run; `NotEnoughFileDescriptors` compares the run to the process;
  `SampleAlignedToAnotherReference` compares a sample's reads to the run's assembly;
  `CatalogIsForAnotherReference` compares the segments' catalog to it;
  `ReferenceCheckedAgainstAnotherGenome` compares the run's two views of its own reference;
  `AnalysedRegionsDiffer` compares two psps' recorded ground to each other; and
  `SegmentationInputsDiffer` compares a psp to the run — spec §6.2, §7.1a. **The first six are
  direct mode's and are built (2026-08-31); the last two only psp mode can reach.**
- **⚑ `CatalogIsForAnotherReference` exists because the layer that should hold it cannot, on the
  path that matters.** The catalog's own open guards every digest comparison on the reference
  having one, and the `.fai` path — the ordinary one — has none until its background read
  finishes. So the catalog is admitted on contig names, lengths and order, and nothing compares
  the digests afterwards. The run holds both values at construction and is the first place that
  can. Without it a catalog from another build of the same assembly routes every repeat tract to
  the wrong position, genome-wide, with nothing to notice.
- **`SampleAppearsTwice` is psp mode's alone**: in direct mode several files naming one
  individual are one sample by construction, since the read-group table groups them
  (`callers.rs`, built 2026-08-31), and a cohort that opens the same individual twice is a thing
  only two files can express.
- **Three of the six run before a single file is opened**, and the three that compare references
  after, because one of them reads the checksums each open captured. Each of the first three
  condemns the whole run, so opening a thousand files first would only make the message slower.
- **The refusals are refusals; the two built variants are failures.** `Catalog` and
  `OpeningSample` say a run could not read something it was handed, where the five above say a
  run was handed things that do not go together. They live in one enum because a caller has one
  `Result` either way, and they are told apart by what a person does next: fetch the file, or fix
  the arguments.
- **`CallersInFlight`, `SamplesInFlight`, `Workers` are newtypes over `NonZeroUsize`**, and
  **none of them is given a default here** — spec §11 questions 2 and 3 own the sweeps that set
  the first two, and the third is not on the default path.
- **`RepeatCatalogHeader` is reused whole as the reference-and-catalog identity.** No
  `ReferenceDigest` or `CatalogIdentity` is minted for the segmentation
  ([`repeat_catalog/mod.rs:283`](../../../../src/ng/repeat_catalog/mod.rs)). The parameters file
  has its own `ReferenceDigest` for a different job — proving the numbers were fitted against this
  assembly ([`parameters_file.md`](../spec/parameters_file.md) §6) — and the two are not unified.

---

## 7. Reconciliation with existing code

Every row read at the cited line on 2026-08-31.

| this doc | existing code | how they meet |
|---|---|---|
| the segments the loops advance over | `TypedRegion` [`region_typing/mod.rs:144`](../../../../src/ng/region_typing/mod.rs), `RegionKind` [`:168`](../../../../src/ng/region_typing/mod.rs) | consumed as-is; `Segmentation` holds the list; nothing re-classifies |
| the flowing item | `SampleLocusObservations` [`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs) | already owned and lifetime-free, which is what lets an observation outlive the worker that minted it |
| the source trait (§2) | `ObservationSource` [`cohort_merge/observation_cache.rs:70`](../../../../src/ng/run/cohort_merge/observation_cache.rs), blanket impl for iterators [`:98`](../../../../src/ng/run/cohort_merge/observation_cache.rs) | **built, and described here rather than moved**; the walker implements it |
| the walker behind the alignment source | `SampleLocusObservationsIterator` [`locus_generation/mod.rs:921`](../../../../src/ng/locus_generation/mod.rs) | **one per sample for the whole run**; drop order load-bearing [`:1028`](../../../../src/ng/locus_generation/mod.rs) |
| the per-sample read cursor | `SampleReads` [`read/input/mod.rs:398`](../../../../src/ng/read/input/mod.rs), `cursor` [`:623`](../../../../src/ng/read/input/mod.rs) | one `SampleReads` per sample, one owned cursor per walker (`Send` proven at [`:1441`](../../../../src/ng/read/input/mod.rs)) |
| the reference accessor factory | factory parameter of `cursor` [`read/input/mod.rs:606-611`](../../../../src/ng/read/input/mod.rs) | one accessor per walker; the factory exists because `WindowedRefSeq` is `Send`, not `Sync` |
| "correct in any order, fastest ascending" | per-segment fetch [`pileup/generator.rs:621`](../../../../src/ng/locus_generation/pileup/generator.rs); any-order cursor [`cursor.rs:92`](../../../../src/ng/read/input/cursor.rs), test [`:1207`](../../../../src/ng/read/input/cursor.rs) | the source ordering contract, walker side |
| read-filter tallies | in the cursor [`read/input/mod.rs:620-622`](../../../../src/ng/read/input/mod.rs) | summed at the gatherer's `finish` (spec §8) |
| the merge | `merge_cohort_through_cache` [`serial.rs:146`](../../../../src/ng/run/cohort_merge/serial.rs), `merge_cohort_in_parallel` [`parallel.rs:96`](../../../../src/ng/run/cohort_merge/parallel.rs) | **built**; both accumulate the whole run's loci in `RegionOutcome` [`build.rs:743`](../../../../src/ng/run/cohort_merge/build.rs), which is the oracle shape and not the run's |
| the cohort locus | `CohortObservation` [`build.rs:930`](../../../../src/ng/run/cohort_merge/build.rs) | what a source's observations become and what calling consumes; closes the draft's second `OPEN` |
| the streaming merge's model | `MergedCursors` [`read/input/sample_cursor.rs:178`](../../../../src/ng/read/input/sample_cursor.rs) | argmin over heads, keys beside heads. **The draft cited `MergedRegionReads`, a name from a superseded design that is not in the tree** (spec §3.2) |
| candidate selection | `select_generic` [`allele_candidates/generic.rs:81`](../../../../src/ng/calling/allele_candidates/generic.rs), `CandidateSelectionConfig` [`allele_candidates/mod.rs:46`](../../../../src/ng/calling/allele_candidates/mod.rs) | called per locus; the tract path is specified and unbuilt |
| the merge's loci → what the loop reads | `shape_generic_locus` [`evidence_shaping.rs:403`](../../../../src/ng/calling/evidence_shaping.rs), `shape_ssr_locus` [`:443`](../../../../src/ng/calling/evidence_shaping.rs) | data-shaping only; **where the three sample numberings are joined**, documented in that file's header |
| the call | `LocusGenotyper` [`inference/mod.rs:619`](../../../../src/ng/calling/inference/mod.rs), `call_locus` [`:629`](../../../../src/ng/calling/inference/mod.rs) → `LocusInference` [`calling/mod.rs:2942`](../../../../src/ng/calling/mod.rs) | one cohort locus in, genotypes out; the caller supplies `CallingScratch` [`calling/mod.rs:1241`](../../../../src/ng/calling/mod.rs) per worker |
| the fitted parameters | `RunParameters` [`calling/run_parameters.rs:98`](../../../../src/ng/calling/run_parameters.rs), borrowed as `FrozenParameters` [`calling/mod.rs:627`](../../../../src/ng/calling/mod.rs) | **the draft's `ModelParams` does not exist**; this is the type, and it holds only fitted numbers — the loop's own config is `RunnableCallingLoopConfig` [`inference/mod.rs:518`](../../../../src/ng/calling/inference/mod.rs) |
| the callers' item | `VcfRecord` [`vcf/mod.rs:84`](../../../../src/ng/vcf/mod.rs), `assemble_record` [`vcf/assemble.rs:117`](../../../../src/ng/vcf/assemble.rs), `VcfWriter` [`vcf/writer.rs:57`](../../../../src/ng/vcf/writer.rs) | built; the writer consumes the caller's stream in genome order |
| the census accumulator | `CensusWriter` [`census.rs:1928`](../../../../src/ng/parameter_estimation/joint/census.rs), `add_locus` [`:2087`](../../../../src/ng/parameter_estimation/joint/census.rs), `mark_walked` [`:2106`](../../../../src/ng/parameter_estimation/joint/census.rs), `finish` [`:2378`](../../../../src/ng/parameter_estimation/joint/census.rs) → `SampleCensusEvidence` [`:1378`](../../../../src/ng/parameter_estimation/joint/census.rs) | owned by the gatherer, fed at the yield point |
| the psp header's first consumer | `PileupIdentity::of_header` [`census_file.rs:96`](../../../../src/ng/parameter_estimation/joint/census_file.rs), `freshness` [`:131`](../../../../src/ng/parameter_estimation/joint/census_file.rs) | digests the psp header's reference, analysed regions, read filters and command line, plus the record count, so a psp header must carry all four (spec §6.1) |
| the census file | `write_census` [`census_file.rs:200`](../../../../src/ng/parameter_estimation/joint/census_file.rs), `open_census` [`:426`](../../../../src/ng/parameter_estimation/joint/census_file.rs) | the walk stage's loop writes `finish()`'s result through it |
| `SegmentationInputs`'s sibling in the fit | `RecordingTerms` [`census.rs:1024`](../../../../src/ng/parameter_estimation/joint/census.rs) | same shape, different stage; not unified — each stage refuses in its own vocabulary |
| `SegmentationInputs::catalog` | `RepeatCatalogHeader` [`repeat_catalog/mod.rs:283`](../../../../src/ng/repeat_catalog/mod.rs) — `reference_md5` [`:291`](../../../../src/ng/repeat_catalog/mod.rs), `built_under` [`:295`](../../../../src/ng/repeat_catalog/mod.rs) | **reused whole** — no identity type is minted |
| `SegmentationInputs::routing` | `StrRepeatCriteria` [`repeat_catalog/criteria.rs:61`](../../../../src/ng/repeat_catalog/criteria.rs); `min_purity: f32` [`segment_criteria.rs:502`](../../../../src/ng/region_typing/segment_criteria.rs) | stored, not restated; the `f32` is why the type is `PartialEq`, not `Eq` |
| `SegmentationInputs::analysed` | `GenomeRegions` [`region_typing/mod.rs:77`](../../../../src/ng/region_typing/mod.rs), `whole_contigs` [`:87`](../../../../src/ng/region_typing/mod.rs), `from_bed_path` [`:100`](../../../../src/ng/region_typing/mod.rs) | **reused whole** — "whole genome" is already the region set covering every contig |
| BED-edge clipping | `clips_at_a_bed_edge` [`region_typing/mod.rs:471`](../../../../src/ng/region_typing/mod.rs), emission rule [`:478-490`](../../../../src/ng/region_typing/mod.rs) | **consumed, not decided** — every loop sees finished segments |
| the resident index the psp reader must not copy | `BlockIndexEntry` [`src/psp/index.rs:42`](../../../../src/psp/index.rs), `decode_index` [`:110`](../../../../src/psp/index.rs) | **a model of what not to build** — 3.8 MB a file at 5 kb blocks, multiplied by the cohort (spec §7.2) |
| the open-file bill | [`examples/dhat_ng_open_files.rs`](../../../../examples/dhat_ng_open_files.rs) | where 11–15 MiB per open alignment file was measured (spec §5.1) |
| `GenomeRegion`, `GenomePosition`, `ContigId`, `Bp` | [`ng/types.rs:79`](../../../../src/ng/types.rs), [`:60`](../../../../src/ng/types.rs), [`:13`](../../../../src/ng/types.rs), [`:185`](../../../../src/ng/types.rs) | used as-is |

---

## 8. Open items

Genuinely open design questions:

- **OPEN: which plan wires a called locus to a `VcfRecord`.** This document's callers yield
  records; the mapper `assemble_record` is built; what is unsettled is whether the run driver's
  plan or the VCF module's plan does the joining, since both name it as their own work
  ([`../impl_plan/run_driver_direct_mode.md`](../impl_plan/run_driver_direct_mode.md) Milestone D,
  [`../impl_plan/vcf_output.md`](../impl_plan/vcf_output.md) Milestones D and E). **Owner's, and it
  is a sequencing question rather than a design one** — nothing about either shape changes with the
  answer.
- **OPEN: which threads do the genotype arithmetic, once several loci are genotyped at a time**
  (§3.2). Two readings, and they are the same run until the second locus starts before the first
  has finished. Calling *inside* the merge's builder means the threads are the merge's own region
  builders — they exist, in `merge_cohort_in_parallel`, and are off by default (spec §3.5).
  Calling *after* the builder means the merge stays on one thread and hands each finished locus to
  a separate set of workers, which is what spec §3.5 describes.
  **⛦ Owner's ruling, 2026-08-31: build against the single-threaded merge and settle this from a
  measurement.** A single-threaded run calls in the same place under either reading, so the run
  driver reaches genotypes from alignment files without answering it. **Whether the region
  batching is kept at all is part of the same question** — it is switched off today on a measured
  1.4× at eight threads, and nothing says it earns its place once genotyping is in the mix.
- **OPEN: the cheap question a source cannot be asked** — spec §10's second entry. It costs direct
  mode nothing and blocks nothing here.

Owed to [`cohort_merge.md`](cohort_merge.md) §8, all four gated on `RunError` landing in this
module, and **none of them is scheduled by this document**:

- the organiser should hold the cache, which it cannot while the cache is generic over its source's
  error type;
- `RunEndedShort` ([`organise.rs:59`](../../../../src/ng/run/cohort_merge/organise.rs)) should fold
  into `RunError`, together with `ObservationExceedsReachCeiling`;
- `ObservationCache::cover` and `evict_before` should become private once the organiser is their
  only caller;
- a source whose observations go backwards should return `RunError` rather than trip an assertion
  — which starts to matter when observations are decoded from a file rather than minted here.

Implementation-time confirmations:

- **Is `SampleReads` `Sync`?** §3.4 holds one per sample and the walkers advance on one thread, so
  direct mode does not need it; the parallel merge would. The suite proves the *cursor* `Send`
  ([`read/input/mod.rs:1441`](../../../../src/ng/read/input/mod.rs)) and nothing asserts `Sync` for
  `SampleReads`. Add the assertion beside that test when it is first needed.
- **Where the summed read-filter tallies ride at `finish`.** Beside the census in a small outcome
  struct, or a separate accessor — pin when coding; they must not be droppable silently (spec §8).
- **Whether `MergeParameters` is worth a struct** (§3.4).
- **Cost of comparing `SegmentationInputs` per sample at open.** A `RepeatCatalogHeader`
  comparison walks a per-contig vector; fine at hundreds of contigs and thousands of samples. If it
  is not, compare a digest — but keep the header stored, so a refusal can still name a field.
- **How the descriptor headroom is read.** `RLIMIT_NOFILE` on Unix; the count needed is two
  descriptors per sample for a CRAM and its index (spec §7.1a).
- **Which contig fields the reference comparison reads**, and whether an alignment header naming
  *extra* contigs the run does not analyse is a refusal or is ignored. The check itself is settled
  — §5's `SampleAlignedToAnotherReference`, owner's ruling 2026-08-31 — because spec §6.2's
  analysed-regions refusal is a psp fact and direct mode hands one segmentation to every sample, so
  that comparison cannot differ. What is left is how strict the comparison is.
- **`SampleInput` and `CensusConfig`** — concrete fields pinned when the constructors are coded.

---

## 9. Test shape

Unit tests beside each file; the run-level oracles are spec §12 and belong in `tests/`.

- `segments.rs`: `build` records the inputs it was given; `first_difference` names each field
  when that field alone is mutated.
- `walker.rs`: the observations a walker yields equal those the iterator yields driven directly,
  position for position — the merge's own fixtures cannot check this, because they are the
  in-memory sources. **Built 2026-08-31**, over the real generic generator and a real indexed BAM,
  and it carries **two things this document did not anticipate**.

  **The segment-independence oracle landed here rather than in `tests/`**, against the rule above.
  It is about one sample's walk, which is what this file owns, and an integration test would need
  the same three-read BAM to say the same thing. Recorded rather than quietly done; move it if the
  rule is meant to hold without exception.

  **And that oracle is not literally true, in one field.** Spec §12 asks that a segment walked
  alone emit *exactly* what the same span emits inside a whole walk. Measured, everything is equal
  but the **chain ids**: `SequenceObservation::chain_ids` says in its own documentation that "an id
  names a read within one walk", and the allocator counts up across a whole walk and survives the
  per-chromosome reset — so the `chr2` read is id 0 walked alone and id 4 walked fourth. No
  implementation of a walk-scoped id can satisfy "exactly". The test compares the **grouping** —
  every id replaced by the order of its first appearance — which still catches a read split in
  two, two reads merged, or a locus that lost its witnesses. **The spec sentence is the owner's to
  amend**, and the same question is owed to spec §12's first oracle, which asks for byte-identical
  psps across worker counts and would inherit the same problem once chain ids are written to a
  file.
- `callers.rs`: each construction refusal fires on its own, on a fixture that trips that one and no
  other; a cohort of the same loci called at one caller in flight and at sixteen yields the same
  records in the same order.
- `gatherer.rs`: everything yielded was counted by the census and everything completed was
  marked walked, empty segments included; tallies from several workers sum at `finish`.

**The regression anchor is spec §12.3** — the same cohort and parameters through
`AlignedFilesVariantCaller` and through the psp route give the same VCF; it is the only test that
can fail when the psp does not carry something the caller needs. The gatherer's own anchor is
spec §12.1: one sample gathered at 1, 2, 4, 8 and 16 workers gives byte-identical psps apart from
the header's timestamp.
