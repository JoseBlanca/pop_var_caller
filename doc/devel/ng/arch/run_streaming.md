# ng — driving a run: types and interfaces

*Status: architecture draft (2026-08-16), companion to the spec
[`../spec/run_streaming.md`](../spec/run_streaming.md) (the design and every *why*) and to the
shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for types,
verbs for functions, **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the
**contract** is the deliverable. This document does not re-argue a decision — the spec section is
cited instead.*

The public surface is the spec's three objects (spec §1): `AlignedFilesVariantCaller`,
`PspVariantCaller`, `SampleObservationGatherer` — each an iterator. Everything else in this
document is crate-private machinery inside one of them.

## Module home

`src/ng/run/`, a new folder beside `locus_generation/` and `parameter_estimation/`:

```
src/ng/run/
├── mod.rs        – re-exports the three objects + RunError
├── segments.rs   – Segmentation + SegmentationInputs: the run's segments and what they
│                   were computed from
├── source.rs     – the ObservationSource trait + the alignment-backed walker impl
├── gatherer.rs   – SampleObservationGatherer
├── calling.rs    – call_vars_in_segment, k_way_merge, and the two caller iterators
└── psp_header.rs – PspHeader: the values every psp records
```

The psp **writer** and **reader** are deliberately not files here: their shapes belong to the
encoding spec (spec §10). This document defines what they plug into — the writer consumes the
gatherer's iterator (spec §5.2), the reader implements `ObservationSource` (§2) — and the header
values they write and check (§4). The walk stage's per-sample loop (spec §5.2) is composition the
CLI owns; it needs no type of its own.

**This revises [`module_layout.md`](module_layout.md)'s `pipeline.rs` entry**, which named one
file holding "the `CallerRecipe` + the driver that runs it end-to-end". The recipe stays where
that document puts it; the drivers are the two caller objects here. Recorded rather than quietly
changed, because the tree in that document still shows the old placement.

---

## 1. The segments (`segments.rs`)

The segments every loop advances over are the typed-region generator's segments, held once and
shared read-only by every worker (spec §4.2, §9). No grouping routine exists: the loop unit is
one segment (spec §4.4), so `Segmentation` is a list plus the record of its inputs — the previous
draft's `group_toward` is retired with the size choices it served.

```rust
/// The run's segments, in genome order, plus the values they were computed from. A function of
/// the reference, the catalog, the routing criteria and the analysed regions — no sample's
/// reads — so it is identical in every sample of the run (spec §4.2).
pub struct Segmentation {
    inputs: SegmentationInputs,
    segments: Vec<TypedRegion>,
}

impl Segmentation {
    /// Consumes the typed-region generator's stream once. The inputs record is assembled here,
    /// so the values a segmentation was built from and the values it reports cannot disagree.
    pub fn build(
        segments: impl Iterator<Item = Result<TypedRegion, RepeatCatalogError>>,
        analysed: GenomeRegions,
        catalog: RepeatCatalogHeader,
        routing: StrRepeatCriteria,
    ) -> Result<Self, RunError>;

    pub fn inputs(&self) -> &SegmentationInputs;
    /// In genome order; a segment never crosses a contig and is never cut (spec §4.3).
    pub fn segments(&self) -> &[TypedRegion];
}
```

`SegmentationInputs` is both the psp header's core (§4) and the operand of the file-against-run
check (spec §6.2). The parameters fit holds the equivalent object for its own compatibility
check: `RecordingTerms`
([`census.rs:995`](../../../../src/ng/parameter_estimation/joint/census.rs)).

```rust
#[derive(Clone, PartialEq)]
pub struct SegmentationInputs {
    /// The catalog file's own header, **reused whole rather than restated**: it already carries
    /// the whole-reference MD5, the criteria the catalog was built under, the scan weights and
    /// the tool version (`repeat_catalog/mod.rs:283`).
    pub catalog: RepeatCatalogHeader,
    /// The criteria the *reader* asked with. Not the same value as `catalog.built_under` — the
    /// catalog is built below every routing floor so a reader filters rather than re-scans —
    /// and it is this one that decides where a segment ends.
    pub routing: StrRepeatCriteria,
    /// The regions the run was asked to analyse (`region_typing/mod.rs:77`). The field a user
    /// actually changes between runs, and the one compared across the cohort (spec §6.2).
    pub analysed: GenomeRegions,
}

impl SegmentationInputs {
    /// The name of the first field that differs, for the error message; `None` when they
    /// agree. **A name rather than a `bool`**: "these two segmentations differ" leaves the
    /// user nothing to fix (spec §6.1).
    pub fn first_difference(&self, other: &Self) -> Option<&'static str>;
}
```

`PartialEq` and not `Eq`, because `StrRepeatCriteria`
([`repeat_catalog/criteria.rs:61`](../../../../src/ng/repeat_catalog/criteria.rs)) wraps
`SsrSegmentCriteria`, whose `min_purity` is an `f32`
([`segment_criteria.rs:502`](../../../../src/ng/region_typing/segment_criteria.rs)).

---

## 2. The source (`source.rs`)

One trait carries the whole difference between the two callers (spec §3.3). The item is
`SampleLocusObservations` directly — no wrapper: an empty segment yields nothing, and
"analysed-but-empty" is carried by the psp header plus trailer and by the gatherer's own
`mark_walked`, not by an empty container (spec §5.2, §8).

```rust
/// Answers one question: this sample's observations in this segment, in coordinate order.
pub trait ObservationSource {
    fn observations_in(
        &mut self,
        segment: &GenomeRegion,
    ) -> impl Iterator<Item = Result<SampleLocusObservations, RunError>> + '_;
}
```

**Contract.** The yield is exactly the observations minted inside `segment`, in coordinate order,
produced lazily — the merge pulls, the source decodes or walks as pulled (spec §3.2). Segments may
be asked in any order and the answer is the same; ascending is the fast path, and a backward jump
costs a seek (spec §8). A source serves one consumer at a time (`&mut self`); a parallel loop
gives each in-flight segment its own source per sample over the sample's one shared open file
(spec §3.4, §7).

**Two implementations:**

- **The walker** (this file): owns a `SampleCursor` (`SampleReads::cursor` takes `&self` and
  returns an owned, `Send` cursor — [`read/input/mod.rs:623`](../../../../src/ng/read/input/mod.rs),
  test [`:1441`](../../../../src/ng/read/input/mod.rs)), a reference accessor from the factory
  (`WindowedRefSeq` is `Send` and deliberately not `Sync` —
  [`read/input/mod.rs:606-611`](../../../../src/ng/read/input/mod.rs)), and a generator set,
  whose drop order is load-bearing
  ([`locus_generation/mod.rs:707-737`](../../../../src/ng/locus_generation/mod.rs)). Spec §8 is
  the trap list this ownership shape honours.
- **The psp reader** (with the encoding, spec §10): a cursor over one open psp that decodes
  whichever blocks overlap the segment and keeps the one it is in. Its resident state is the
  file's coarse index plus one decoded block; the per-open-file share is bounded at tens of
  kilobytes (spec §7.2). No block is visible in its interface.

---

## 3. The three objects and the calling core (`gatherer.rs`, `calling.rs`)

All three objects share one internal skeleton (spec §3.4): segments dealt to `workers` threads
in genome order, at most `look_ahead` segments in flight, results drained at the yield point in
genome order. The knobs are newtypes over `NonZeroUsize` — a count whose zero is illegal:

```rust
pub struct Workers(pub NonZeroUsize);
/// Segments in flight beyond the next to yield — each object's one memory knob (spec §3.4,
/// §7.1). No default is proposed; spec §11 open question 2 names the sweep.
pub struct LookAhead(pub NonZeroUsize);
/// The walk stage's samples-at-once. Default 1 — one open alignment file, cohort-independent
/// peak (spec §5.2).
pub struct SamplesInFlight(pub NonZeroUsize);
```

### 3.1 `SampleObservationGatherer`

```rust
/// One sample's observations in genome order, census accumulated as they pass. psp mode's
/// walk stage is a loop of these, one per sample (spec §5.2).
pub struct SampleObservationGatherer { /* pool, sources, CensusWriter, tallies */ }

impl SampleObservationGatherer {
    pub fn new(
        sample: &SampleInput,
        segmentation: &Segmentation,
        census: &CensusConfig,
        workers: Workers,
        look_ahead: LookAhead,
    ) -> Result<Self, RunError>;

    /// After the iterator is exhausted: the census the walk accumulated. Calling it earlier is
    /// an error — the census's per-stratum sums are complete only at the end (spec §5.2).
    pub fn finish(self) -> Result<SampleCensusEvidence, RunError>;
}

impl Iterator for SampleObservationGatherer {
    type Item = Result<SampleLocusObservations, RunError>;
}
```

**Contract.** Yields in genome order. At the yield point — single-threaded by construction —
every observation passes `CensusWriter::add_locus` and every completed segment passes
`mark_walked`, empty segments included
([`census.rs:1965,1984`](../../../../src/ng/parameter_estimation/joint/census.rs)); what the
iterator yields and what the census counted are the same stream, which is the whole of spec
§5.2's two closures. Read-filter tallies are per-cursor
([`read/input/mod.rs:620-622`](../../../../src/ng/read/input/mod.rs)) and summed at `finish` —
unsummed, drop rates under-report by the worker count (spec §8).

### 3.2 The calling core

```rust
/// The one merge-and-call, over one segment. Nothing inside can tell a walker from a psp
/// reader (spec §3.1, goal 1). Serial, straight-line; the loop and its bookkeeping live in
/// the callers.
fn call_vars_in_segment(
    segment: &GenomeRegion,
    sources: &mut [impl ObservationSource],
    parameters: &ModelParams,
    variants: &mut Vec<Variant>,   // OPEN — the emission step's document owns Variant's shape
) -> Result<(), RunError>;

/// Streaming cohort merge: one head per source, keyed on coordinates, yielding one cohort
/// observation at a time. Modeled on MergedRegionReads — argmin over per-stream heads, keys
/// beside the heads (sample_reads.md §4). The reconciliation of differing spans inside a
/// cohort observation is the deferred merge spec's (spec §3.2, §10).
fn k_way_merge<'a>(
    per_sample: &'a mut [impl Iterator<Item = Result<SampleLocusObservations, RunError>>],
) -> impl Iterator<Item = Result<CohortObservation, RunError>> + 'a;
```

### 3.3 The two callers

```rust
/// Direct mode (spec §5.1). Holds every sample's SampleReads open for the whole run —
/// 11–15 MiB each — plus the shared read-only state; each in-flight segment's task owns its
/// per-sample walkers.
pub struct AlignedFilesVariantCaller { /* segmentation, SampleReads per sample, params, pool */ }

impl AlignedFilesVariantCaller {
    pub fn new(
        samples: &[SampleInput],
        segmentation: Segmentation,
        parameters: ModelParams,
        workers: Workers,
        look_ahead: LookAhead,
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
        parameters: ModelParams,
        workers: Workers,
        look_ahead: LookAhead,
    ) -> Result<Self, RunError>;
}

impl Iterator for AlignedFilesVariantCaller { type Item = Result<Variant, RunError>; }
impl Iterator for PspVariantCaller        { type Item = Result<Variant, RunError>; }
```

**Contract, both variant callers.** Variants in genome order, identical at every worker count and
look-ahead (spec §12.2), identical between the two callers on one cohort with fixed parameters
(spec §12.3 — the regression anchor). Iteration ends at the first `Err`; direct mode leaves
nothing to clean up, and a psp without a valid trailer is refused at `open` (spec §9).

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

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    /// A worker failed over a known stretch. **Both the sample and the span** — neither alone
    /// locates a failure in a run over thousands of samples (spec §9).
    #[error("sample {sample}: failed over {segment}")]
    WorkerFailed {
        sample: String,
        segment: GenomeRegion,
        source: Box<dyn std::error::Error + Send + Sync>,
    },
    /// Two psps were analysed over different segments, so the samples are not comparable —
    /// the cohort refusal (spec §6.2). Open question 5 may later soften this to
    /// intersection-calling; until then it refuses.
    #[error("samples {left} and {right} were analysed over different segments")]
    AnalysedRegionsDiffer { left: String, right: String },
    /// One psp's recorded catalog or routing criteria differ from the run's own, so the
    /// segments this run loops over are not the segments the file's observations were minted
    /// inside — the file-against-run refusal (spec §6.2). `field` is
    /// `SegmentationInputs::first_difference`'s answer.
    #[error("psp for sample {sample} was written under a different {field}")]
    SegmentationInputsDiffer { sample: String, field: &'static str },
    /// A psp ended without a valid trailer: an interrupted walk, not a short sample (spec §9).
    #[error("psp for sample {sample} is incomplete")]
    IncompletePsp { sample: String },
    #[error("i/o while reading or writing a run's files")]
    Io(#[from] std::io::Error),
}
```

---

## 6. Design decisions — decided

- **Three public objects, each an iterator; everything else crate-private.** No caller of this
  module names a work unit, a block, or a range — spec §1, §3.
- **The loop unit is one segment; no grouping routine exists.** `Segmentation` is a list;
  the previous draft's `group_toward`, `SampleSpanLoci`, `SampleLociSource`, `LociSink` and the
  three free entry points are retired with the size choices they served — spec §4.4.
- **One source trait, two implementations; the calling core consumes the trait and never
  invokes a walk.** The whole of "one calling function, whichever mode" — spec §3.1, §3.3.
- **`k_way_merge` streams and keys on coordinates**, modeled on `MergedRegionReads`
  ([`sample_reads.md`](sample_reads.md) §4); its residency is the frontier — spec §3.2.
- **What flows is `SampleLocusObservations`, bare.** No span-carrying wrapper: walked-empty
  ground is recorded by the gatherer's `mark_walked` and readable back through the header plus
  trailer, so an empty container had no remaining job — spec §5.2, §8.
- **The census accumulator lives inside the gatherer, fed at the yield point.** One stream, two
  consumers impossible to desynchronise — spec §5.2.
- **The psp reader serves segments; blocks exist only inside the writer and reader** — spec §3.3;
  their sizing is wholly the encoding spec's — spec §6.3, §10.
- **The header carries no boundary digest and no writer version** — spec §6.3, **flagged for the
  owner as a reversal of the previous draft**.
- **Two refusal variants, two axes.** `AnalysedRegionsDiffer` compares files to each other;
  `SegmentationInputsDiffer` compares a file to the run — spec §6.2.
- **`Workers`, `LookAhead`, `SamplesInFlight` are newtypes over `NonZeroUsize`** — zero is
  illegal for each; look-ahead is each object's one memory knob — spec §3.4, §7.1.
- **`RepeatCatalogHeader` is reused whole as the reference-and-catalog identity.** No
  `ReferenceDigest` or `CatalogIdentity` is minted
  ([`repeat_catalog/mod.rs:283`](../../../../src/ng/repeat_catalog/mod.rs)).

---

## 7. Reconciliation with existing code

Every row read at the cited line.

| this doc | existing code | how they meet |
|---|---|---|
| the segments the loops advance over | `TypedRegion` [`region_typing/mod.rs:144`](../../../../src/ng/region_typing/mod.rs), `RegionKind` [`:168`](../../../../src/ng/region_typing/mod.rs) | consumed as-is; `Segmentation` holds the list; nothing re-classifies |
| the flowing item | `SampleLocusObservations` [`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs) | already owned and lifetime-free, which is what lets an observation outlive the worker that minted it |
| the walker behind the alignment source | `SampleLocusObservationsIterator` [`locus_generation/mod.rs:706`](../../../../src/ng/locus_generation/mod.rs) | one per in-flight segment's task; drop order load-bearing [`:707-737`](../../../../src/ng/locus_generation/mod.rs) |
| the per-worker read cursor | `SampleReads` [`read/input/mod.rs:398`](../../../../src/ng/read/input/mod.rs), `cursor` [`:623`](../../../../src/ng/read/input/mod.rs) | one shared `SampleReads` per sample, one owned cursor per task (`Send` proven at [`:1441`](../../../../src/ng/read/input/mod.rs); `Sync` to confirm — §8) |
| the reference accessor factory | factory parameter of `cursor` [`read/input/mod.rs:606-611`](../../../../src/ng/read/input/mod.rs) | one accessor per task; the factory exists because `WindowedRefSeq` is `Send`, not `Sync` |
| "correct in any order, fastest ascending" | per-segment fetch [`pileup/generator.rs:621`](../../../../src/ng/locus_generation/pileup/generator.rs); any-order cursor [`cursor.rs:92-96`](../../../../src/ng/read/input/cursor.rs), test [`:1207`](../../../../src/ng/read/input/cursor.rs) | the `ObservationSource` ordering contract, walker side |
| read-filter tallies | in the cursor [`read/input/mod.rs:620-622`](../../../../src/ng/read/input/mod.rs) | summed at the gatherer's `finish` (spec §8) |
| the streaming merge's model | `MergedRegionReads` [`sample_reads.md`](sample_reads.md) §4 | argmin over heads, keys beside heads; `k_way_merge` copies the shape, swaps the item and the yield |
| the census accumulator | `CensusWriter` [`census.rs:1806`](../../../../src/ng/parameter_estimation/joint/census.rs), `add_locus` [`:1965`](../../../../src/ng/parameter_estimation/joint/census.rs), `mark_walked` [`:1984`](../../../../src/ng/parameter_estimation/joint/census.rs), `finish` [`:2252`](../../../../src/ng/parameter_estimation/joint/census.rs) → `SampleCensusEvidence` [`:1349`](../../../../src/ng/parameter_estimation/joint/census.rs) | owned by the gatherer, fed at the yield point; its doc comment already says it borrows the same locus stream its siblings see |
| the psp header's first consumer | `PileupIdentity::of_header` [`census_file.rs:91`](../../../../src/ng/parameter_estimation/joint/census_file.rs), `freshness` [`:126`](../../../../src/ng/parameter_estimation/joint/census_file.rs) | **already built on `main`** — it digests the psp header's reference, analysed regions, read filters and command line, plus the record count, so a psp header must carry all four (spec §6.1); which bytes exactly is the encoding spec's business |
| the census file | `write_census` [`census_file.rs:195`](../../../../src/ng/parameter_estimation/joint/census_file.rs), `open_census` [`:421`](../../../../src/ng/parameter_estimation/joint/census_file.rs) | the walk stage's loop writes `finish()`'s result through it |
| `SegmentationInputs`'s sibling in the fit | `RecordingTerms` [`census.rs:995`](../../../../src/ng/parameter_estimation/joint/census.rs) | same shape, different stage; not unified — each stage refuses in its own vocabulary |
| `SegmentationInputs::catalog` | `RepeatCatalogHeader` [`repeat_catalog/mod.rs:283`](../../../../src/ng/repeat_catalog/mod.rs) — `reference_md5` [`:291`](../../../../src/ng/repeat_catalog/mod.rs), `built_under` [`:295`](../../../../src/ng/repeat_catalog/mod.rs) | **reused whole** — no identity type is minted |
| `SegmentationInputs::routing` | `StrRepeatCriteria` [`repeat_catalog/criteria.rs:61`](../../../../src/ng/repeat_catalog/criteria.rs); `min_purity: f32` [`segment_criteria.rs:502`](../../../../src/ng/region_typing/segment_criteria.rs) | stored, not restated; the `f32` is why the type is `PartialEq`, not `Eq` |
| `SegmentationInputs::analysed` | `GenomeRegions` [`region_typing/mod.rs:77`](../../../../src/ng/region_typing/mod.rs), `whole_contigs` [`:87`](../../../../src/ng/region_typing/mod.rs), `from_bed_path` [`:100`](../../../../src/ng/region_typing/mod.rs) | **reused whole** — "whole genome" is already the region set covering every contig |
| BED-edge clipping | `clips_at_a_bed_edge` [`region_typing/mod.rs:471`](../../../../src/ng/region_typing/mod.rs), emission rule [`:482-488`](../../../../src/ng/region_typing/mod.rs) | **consumed, not decided** — every loop sees finished segments |
| the resident index the psp reader must not copy | `BlockIndexEntry` [`src/psp/index.rs:42`](../../../../src/psp/index.rs), `decode_index` [`:110`](../../../../src/psp/index.rs) | **a model of what not to build** — 3.8 MB a file at 5 kb blocks, multiplied by the cohort (spec §7.2) |
| `GenomeRegion`, `ContigId`, `Bp` | [`ng/types.rs:79`](../../../../src/ng/types.rs), [`:13`](../../../../src/ng/types.rs), [`:174`](../../../../src/ng/types.rs) | used as-is |

---

## 8. Open items

Genuinely open design questions:

- **OPEN: `Variant`'s shape.** Named in §3.2 only so the callers have an item type; the emission
  step's document owns it.
- **OPEN: `CohortObservation`'s shape** and the span reconciliation inside it — the deferred
  merge spec (spec §3.2, §10), which must also confirm the merge frontier stays bounded when
  spans overlap.

Implementation-time confirmations:

- **Is `SampleReads` `Sync`?** §2 shares it by reference across tasks; the suite proves the
  cursor `Send` ([`read/input/mod.rs:1441`](../../../../src/ng/read/input/mod.rs)) but nothing
  asserts `Sync` for `SampleReads`. Add an `assert_sync::<SampleReads>()` beside that test; if an
  interior non-`Sync` field turns up, the fix is per-task `SampleReads` plus a note in spec §7's
  formulas — not a lock.
- **Where the summed read-filter tallies ride at `finish`.** Beside the census in a small
  outcome struct, or a separate accessor — pin when coding; they must not be droppable silently
  (spec §8).
- **Worker reuse of cursors and generators across its consecutive segments.** An internal
  optimisation the design permits and does not require (spec §4.4); build the naive per-task
  shape first, measure, then reuse if per-segment overhead shows in a profile.
- **Cost of comparing `SegmentationInputs` per sample at open.** A `RepeatCatalogHeader`
  comparison walks a per-contig vector; fine at hundreds of contigs and thousands of samples. If
  it is not, compare a digest — but keep the header stored, so a refusal can still name a field.
- **`SampleInput` and `CensusConfig`** — concrete fields pinned when the constructors are coded.

---

## 9. Test shape

Unit tests beside each file; the run-level oracles are spec §12 and belong in `tests/`.

- `segments.rs`: `build` records the inputs it was given; `first_difference` names each field
  when that field alone is mutated.
- `source.rs`: a fake source proves the contract — the yield is exactly the segment's
  observations in coordinate order, an empty segment yields nothing, the same segment asked twice
  in different orders answers the same.
- `calling.rs`: `call_vars_in_segment` over fake sources — the merge keys on coordinates, not
  arrival order; a caller fed completions in every permutation yields in genome order; yields at
  look-ahead 1 equal yields at look-ahead 8.
- `gatherer.rs`: everything yielded was counted by the census and everything completed was
  marked walked, empty segments included; tallies from several workers sum at `finish`.

**The regression anchor is spec §12.3** — the same cohort and parameters through
`AlignedFilesVariantCaller` and through the psp route give the same VCF; it is the only test that
can fail when the psp does not carry something the caller needs. The gatherer's own anchor is
spec §12.1: one sample gathered at 1, 2, 4, 8 and 16 workers gives byte-identical psps apart from
the header's timestamp.
