# ng — the generic locus generator: types & interfaces

*Architecture draft (2026-07-27), the code-facing companion to
[`../spec/locus_generation_pileup.md`](../spec/locus_generation_pileup.md) — the second
[`LocusGenerator`](locus_generation.md), consuming the `Generic` region kind and producing one locus
per covered position. It **inherits** the locus type, the contract, the dispatch and the error model
from [`locus_generation.md`](locus_generation.md) (arch) — read that first; this doc adds only what
is generic-path-specific. Grounded in production
[`src/pileup/walker/`](../../../../src/pileup/walker/). The production-side defect this port fixes is
recorded separately in
[`pileup_partial_coverage_ref_fill_2026-07-27.md`](../../reports/research/pileup_partial_coverage_ref_fill_2026-07-27.md).
Naming: **STR** in prose, `ssr` in code. Signatures are illustrative; the **contract** is the
deliverable. Build order:
[prerequisites](../impl_plan/locus_generation_pileup_prerequisites.md) →
[the port](../impl_plan/locus_generation_pileup_port.md) →
[the generator](../impl_plan/locus_generation_pileup_generator.md).*

*⚠ **Branch state, 2026-07-27 — merged.** This worktree is now at `eb2857c`: `main` plus the
read-fetch perf work plus **`ng-read-groups`** (the read group as a first-class object — a new
`AlignedRead`, a new `read_groups.rs`, `src/ng/read/input/` substantially rewritten). Every
`file:line` below was re-verified against that tree after the merge. Two outcomes worth carrying at
the top: **the §2.2 prerequisite survives the rewrite intact**, and **the read type changed** —
`AlignedRead` replaces `MappedRead` as what the stream yields and what `ReadPreparer` consumes.*

## Module home

`src/ng/locus_generation/pileup/` — a **folder**, unlike its `ssr.rs` sibling, because it holds the
whole copied walker plus the generator that wraps it. **Production is not edited** (spec §3): ng
copies rather than reaches in, so no visibility lift and no field on a frozen type.

```
src/ng/locus_generation/pileup/
  mod.rs                – PileupGenerator, its config and counts, the RefSeq→fetcher shim
  genome_walk.rs        – the walk along genome coordinates: advance the position, admit and
                          expire reads, reconcile mates, drive the fold  (production's driver.rs)
  open_record.rs        – the open-record table + fold
  cigar_cursor.rs       – the per-op offset table + the adaptor mask
  decompose.rs          – ReadEvent + the indel base-quality proxy
  active_read_set.rs    – the active set
  chain_id_allocator.rs – chain ids + mate pairing
  errors.rs             – WalkerError
  tests.rs              – #[cfg(test)] production's own end-to-end walker suite, copied
                          verbatim: Milestone A's gate (spec §12). Dies with plan 3.
  copy_fidelity.rs      – #[cfg(test)] ng's own: the textual check that the eight copies
                          are still production's, from outside the files it checks
  parity.rs             – #[cfg(test)] the differential harness (spec §3)
```

**Eight of those are copies** — `genome_walk` through `errors`, plus `tests.rs` — landing verbatim
and changed only once the differential passes (spec §3). `mod.rs`, `copy_fidelity.rs` and `parity.rs`
are ng's own.

*(Inventory corrected 2026-07-29: `tests.rs` and `copy_fidelity.rs` were missing. `tests.rs` is not
optional — A4's gate is production's suite green **against the copy**, which is only meaningful if it
runs against ng's walker. `copy_fidelity.rs` cannot live inside `tests.rs`, because `tests.rs` is one
of the files it checks.)*
`src/ng/read/` gains the read type they all name:

```
src/ng/read/prepared_read.rs – PreparedRead (+ MateRole, ReadLengthError), copied from
                               pileup/walker/mod.rs and extended with read_group (§1.2)
```

**Still reused from production, unchanged and already `pub`:** `CigarOp`, `PileupRecord`,
`AlleleObservation`, `AlleleSupportStats`, `WalkerConfig`'s `DEFAULT_*` constants — none of which ng
modifies, so none of which it copies.

**One file is renamed on the way in: `driver.rs` → `genome_walk.rs`.** It is the only file in the
set named for a *role* rather than for what it owns — its siblings are `active_read_set`,
`cigar_cursor`, `open_record`, `chain_id_allocator` — and "driver" answers *driver of what?* with
nothing. `genome_walk` names the one job that file has and the others do not: advancing a position
cursor along genome coordinates over an active read set. It is deliberately **self-sufficient**,
readable in a stack trace or a grep hit without the folder for context — and ng has more than one
walk (region typing walks the reference too), so the qualifier earns its keep. *Two rejected:*
`pileup.rs`, because the folder is already `pileup/` so every file in it is the pileup and the name
would distinguish nothing; and bare `walk.rs`, which needs its parent to be understood. *One caveat,
recorded not hidden:* a walk covers **one region**, not the genome — `genome_walk` names the axis it
advances along, not the extent (§2.2). The **type** keeps production's `PileupWalker` so the
differential reads as a straight comparison; renaming it is a separable later call. Every other
copied file keeps its name.

## 1. The types

### 1.1 What the generator owns

```rust
/// The generic locus generator: a streaming pileup walk over one `Generic` region.
/// Holds its accessors as fields — the "a generator holds its own accessors" convention
/// (`locus_generation.md` §2). Two reference accessors, not one: `preparer` carries its own
/// (read preparation's rule), `reference` serves the walker's REF fetches (spec §2).
/// **Neither is rebuilt per segment** — that is the ~564k-opens trap (spec §8).
pub struct PileupGenerator<R: RefSeq, P: ReadPreparer> {
    reference: R,
    preparer: P,
    prep_scratch: P::Scratch,
    /// Lives across segments so `next_id` never repeats — but `reset()` is called at every
    /// region end, which clears `pending_mates` and `active_count` while preserving the
    /// counter. Carrying those two across regions cross-pairs mates between contigs and leaks
    /// `active_count` toward `ActiveReadsExhausted` (spec §8).
    chain_ids: ChainIdAllocator,
    config: PileupGeneratorConfig,
    counts: PileupGeneratorCounts,
    /* current region, current walk (§2.2), region clamp bounds */
}

/// This generator's knobs — owned, taken at construction (`locus_generation.md` §5). Every
/// value is **production's, inherited and never measured by ng**; that is the map of what is
/// safe to move (spec §7). Raw `u32`, not `Bp`: the copied walker speaks production's integer
/// widths and a port must not change types under itself (§3).
pub struct PileupGeneratorConfig {
    pub max_snp_column_depth: u32,
    pub max_indel_column_depth: u32,
    pub max_record_span: u32,
    pub mate_lookup_window: u32,
    pub max_active_reads: u32,
}
```

`Default` takes all five straight from production's own `pub const`s —
`DEFAULT_MAX_SNP_COLUMN_DEPTH` and siblings
([walker/mod.rs:67-89](../../../../src/pileup/walker/mod.rs#L67),
[chain_id_allocator.rs:41](../../../../src/pileup/walker/chain_id_allocator.rs#L41)) — **by name, not
by literal**, so there is one source of truth until ng deliberately diverges and the divergence shows
up as a diff. **There is no sixth knob** — the walk is bounded by the region, not by a window size
(§2.2), so nothing here tunes how far it reaches.

```rust
/// Run-level counts, alongside the shared `LocusCounts`. The first seven mirror production's
/// `RunSummary` field for field (spec §7); the last two are ng's.
pub struct PileupGeneratorCounts {
    pub reads_admitted: u64,
    pub record_widen_events: u64,
    pub mate_overlap_positions: u64,
    pub chain_allocations: u64,
    pub active_reads_high_water: u32,
    pub mate_lookup_evictions: u64,
    pub column_depth_truncations: u64,
    /// Reads silent over a whole record footprint — all-`N` or fully adaptor-masked, so never
    /// contributors and invisible to the per-locus tally (spec §6).
    pub reads_silent_over_footprint: u64,
    /// Records dropped by the region clamp. Observably zero-sum across neighbouring regions,
    /// which is how the gap-free tiling argument stays checkable (spec §7).
    pub records_outside_region: u64,
}
```

### 1.2 What changes inside the copied files

The port's changes, as a closed set a reviewer can check against. They reach the haplotype builder
and the walk, not just the fold — the port exists to stop fabricating, so the change belongs in the
function that fabricates (spec §4).

```rust
/// A stretch of reference positions, 1-based inclusive — matching
/// `PreparedRead::{alignment_start, alignment_end}`. Deliberately not `GenomeRegion`, which
/// carries a `ContigId` and a `u64` this fold has no use for.
#[derive(Clone, Copy)]
struct RefSpan { start: u32, end: u32 }

/// What one read currently contributes to one open record. **ng adds `witnessed` and
/// `read_group`.** `witnessed` is in *absolute reference coordinates*, not relative to the
/// footprint, which grows: coverage is resolved against the record's FINAL footprint at
/// `finalise()`, and the read may have expired long before that (`expire_passed` touches no
/// open record), so it cannot be recomputed from the active set (spec §4).
struct FoldedReadState {
    allele_index: usize,
    contribution: AlleleSupportStats,   // production's, minus `placed_start` only (§3)
    chain_id: ChainId,
    read_group: ReadGroupId,
    /// The positions this read actually witnessed inside the record — the union of its event
    /// footprints, **not** its alignment span. The span is blind to `N`, adaptor-masked,
    /// ref-skipped and dropped-indel positions, all of which production fills from the
    /// reference (spec §6).
    witnessed: RefSpan,
}

/// Build the haplotype a read presents, **emitting only what its events cover** — where
/// production fills every uncovered offset from the reference. Returns the extent covered, so
/// the caller can store it. Reads whose events tile the footprint produce byte-identical output
/// to production, which is what keeps the complete class parity-comparable (spec §4).
///
/// `ref_seq` is still needed, and not for filling: an `Insertion`/`Deletion` arm emits the
/// **anchor base** from the reference when no `Match` already emitted it
/// ([open_record.rs:546](../../../../src/pileup/walker/open_record.rs#L546), [:556](../../../../src/pileup/walker/open_record.rs#L556)).
/// That is one base, inside an event the read did witness — see spec §4 for the corner where
/// the read's own base there was dropped.
///
/// `None` = the witnessed positions are non-contiguous (an interior `N`, a ref-skip): one
/// `Observed` run cannot describe that honestly, so the read yields no observation and is
/// counted in `reads_without_observation` (spec §6).
fn apply_events_into(
    buf: &mut Vec<u8>,
    record_pos: u32,
    ref_seq: &[u8],
    events: &[ReadEvent],
) -> Option<RefSpan>;

/// How much of the finished record this read witnessed, in **locus positions**. Resolved once,
/// at `finalise()`, from `witnessed` against the final footprint — never during folding, when
/// the final footprint is not yet known (spec §4, §6).
fn coverage_of(witnessed: RefSpan, record_pos: u32, record_end_exclusive: u32) -> ReadCoverage;
```

**`widen` extends the REF bucket only.** Production appends the extra reference bases to *every*
bucket ([open_record.rs:348](../../../../src/pileup/walker/open_record.rs#L348)). `alleles[0]` is the
record's own reference sequence and genuinely grows; the other buckets hold what reads witnessed and
never grow. A live read re-folds against the wider window and lands wherever its bases put it; an
expired read keeps a bucket whose bases already say exactly what it saw. **This is what makes the
no-fabrication rule implementable at all** — production cannot express it, because the bases live on
the shared bucket and `FoldedReadState` holds none of its own.

**The read group arrives on the read, not beside it** (why it is carried: spec §6). ng's
`PreparedRead` carries `read_group: ReadGroupId` (module home, above), so `ReadContribution` copies
it off the active read like the fields around it and nothing here reconstructs anything.

**`finalise` returns `SampleLocusObservations`**, bucketing on
`(bases, read_coverage, read_group)` and computing `reads_without_observation` /
`reads_discarded_by_cap` per record. `AlleleSupportStats` is ng's own copy.

### 1.3 The shim

```rust
/// Adapts ng's `RefSeq` to the walker's `MultiChromRefFetcher`. **Semantically empty**: both
/// contracts are canonical uppercase `{A,C,G,T,N}`
/// ([fasta/mod.rs:117](../../../../src/fasta/mod.rs#L117)), so this moves bytes and decides
/// nothing (spec §3).
struct RefSeqFetcher<R: RefSeq>(R);
```

> **⚠ Superseded 2026-07-29 — this shim is temporary and plan 3 deletes it (owner).** It exists
> for one reason: the copies were transcribed *verbatim*, so their signatures are production's
> and say `MultiChromRefFetcher`. That is a consequence of the copy, not a design choice, and it
> stops being true the moment the two walkers diverge. **Plan 3's A0** — its first step, while the
> stage-1 differential can still prove the refactor free — has `open_record.rs` take a `RefSeq`
> directly, deletes `RefSeqFetcher` and its error translation, and switches to `fetch_into`, which
> is the allocation note in §4 below. ng then imports neither `MultiChromRefFetcher` nor
> `ChromRefFetchError`.
>
> *(A review also found the name was **deliberately retired** in this codebase for a different
> concept — `fasta/fetcher.rs:20-23`, a 2026-05-23 review — so a grep returns both the retirement
> note and a live type. Deleting it settles that too; a rename was the alternative and was not
> taken.)*

## 2. The interface

### 2.1 The generator

```rust
impl<R: RefSeq, P: ReadPreparer> LocusGenerator<()> for PileupGenerator<R, P> {
    fn begin_segment(&mut self, region: GenomeRegion);
    fn next_locus(&mut self, segment: &(), reads: &SampleReads)
        -> Result<Option<SampleLocusObservations>, LocusGenerationError>;
}
```

The segment payload is `()` because `RegionKind::Generic` carries none — the slot is already
`GeneratorSlot<()>` ([locus_generation/mod.rs:377](../../../../src/ng/locus_generation/mod.rs#L377)).

**Contract.** Lazy, one locus resident at a time, coordinate order. `begin_segment` **records the
region and nothing else** — it cannot fail, and opening a read query can, so the first `next_locus`
opens it and is where an `IngestError` surfaces. `next_locus` returns the next record whose anchor
falls inside the region, or `None` once the walk drains; a record anchored outside is dropped and
tallied, which is what makes neighbouring regions tile without duplicates or holes (spec §2).
Accessors, the preparer's scratch and the chain-id allocator persist **across** segments; only the
walk is per-segment. Errors are fatal and terminal — the walker's own
([errors.rs:12](../../../../src/pileup/walker/errors.rs#L12)) are latched by the walk, and every
variant reaches the caller wrapped (§3).

### 2.2 The walk, and the borrow that shapes it

`reads_in_region` returns `SampleRegionReads<'_, R>`, which **borrows the `SampleReads`**
([read/input/mod.rs:508](../../../../src/ng/read/input/mod.rs#L508)). The chain is
`SampleReads.files: Vec<AlignmentFile>` ([:330](../../../../src/ng/read/input/mod.rs#L330)) →
`&'a AlignmentFile` held by `RegionReads<'a>`
([open_bam.rs:215](../../../../src/ng/read/input/open_bam.rs#L215)), `BorrowedReader<'a>`
([:171](../../../../src/ng/read/input/open_bam.rs#L171)) and `RegionSource<'a>`
([region_query.rs:724](../../../../src/ng/read/input/region_query.rs#L724)) — the borrow exists
because a pooled reader has to be handed back to the file it came from. `LocusGenerator` lends
`reads: &SampleReads` **per call** and carries no lifetime parameter, so a generator **cannot hold
that stream between `next_locus` calls** — and the pileup, unlike the STR generator, yields many
loci per segment, so it must.

**Production does not solve this; it avoids needing to — and that is the useful finding.** Its
walker is generic over `I: IntoIterator<Item = PreparedRead>` and takes it by value, and the value
handed in *is itself a borrowing iterator*: `Box<dyn Iterator<Item = PreparedRead> + '_>`
([stage1_pipeline.rs:229](../../../../src/pop_var_caller/stage1_pipeline.rs#L229)). That is legal
there because the walker is a local driven to completion inside one function and handed to a
closure rather than returned ([:230-232](../../../../src/pop_var_caller/stage1_pipeline.rs#L230)):
the borrow never outlives the scope that created it. **ng's contract is resumable, and *resumable +
borrowed stream + no lifetime on `Self`* is not expressible.** One of the three has to give, and the
cheapest is the borrow.

**So the stream becomes owned, and the walk owns it exactly as production's walker owns its input:**

```rust
/// The walk over one region: **owns** its read stream, so nothing borrows `SampleReads` and
/// nothing has to be materialised. One `PileupWalker` per region, built on the first
/// `next_locus` and drained across the calls that follow (§2.1).
struct RegionWalk<R: RefSeq> { /* PileupWalker<PreparedReads<R>, RefSeqFetcher<R>>, clamp bounds */ }
```

This needs `reads_in_region` to hand back an owned stream, which is **two borrows, not a redesign** —
and the codebase already uses the remedy for both, one field away:

| borrow today | becomes | why it is not self-referential |
|---|---|---|
| `header: &'a sam::Header` in `BamRegionSource`/`CramRegionSource` ([region_query.rs:70](../../../../src/ng/read/input/region_query.rs#L70), [:388](../../../../src/ng/read/input/region_query.rs#L388)) | `Arc<sam::Header>`, cloned from `AlignmentFile` | an independent `Arc`, not a reference *into* the file. Precedent in the same struct: `entries: Arc<[crai::Record]>` ([:393](../../../../src/ng/read/input/region_query.rs#L393)) |
| `file: &'a AlignmentFile` in `BorrowedReader`/`RegionReads` ([open_bam.rs:171](../../../../src/ng/read/input/open_bam.rs#L171), [:215](../../../../src/ng/read/input/open_bam.rs#L215)) — held so `Drop` can return the pooled reader | `Arc<AlignmentFile>`, with `SampleReads.files: Vec<Arc<AlignmentFile>>` | the pool lives on `AlignmentFile`, so it stays shared per file; the clone is one atomic increment per query |

**Re-verified after the `ng-read-groups` merge (2026-07-27): all four survive it unchanged.**
`reads_in_region` still returns `SampleRegionReads<'_, R>`, the enum still carries `'a`,
`AlignmentFile` still owns its `header: sam::Header`, and `BorrowedReader<'a>` / `RegionReads<'a, R>`
still hold `&'a AlignmentFile` for the pool return. Only line numbers moved, and the citations above
are the post-merge ones. So this prerequisite is untouched by a rewrite of the very module it
targets — it is the *read type* that changed (§5), not the borrow.

Everything else in those sources is already owned. Precedent for the shape is not hypothetical:
`AlignmentFile.path` is an `Arc<Path>` for exactly this reason
([open_bam.rs:77](../../../../src/ng/read/input/open_bam.rs#L77) — *"`Arc` so the per-query order
guard can hold it for its error message without an allocation per query"*), and `4bc3ef9` applied
the same medicine to the reference accessor, closing `locus_generation.md` §8's "Arc gap". It is a
change to a *shipped ng module* rather than to frozen production, so it is a **prerequisite pass**
of its own (§4), on the `bundle_threshold` model.

## 3. Decisions — decided (why in the spec)

- **Copy the whole walker into ng; production is edited zero times** — ng needs its own read type
  (§1.2), and all four modules that looked reusable name `PreparedRead` in their signatures, so it
  reaches them anyway. *Rejected:* driving production's walker as a black box (cannot fill three
  fields of the locus type); and the earlier reuse/copy split, whose case evaporated with the read
  type (spec §3).
- **Transcribe first, change second.** Every copied file lands verbatim, still emitting
  `PileupRecord`, and is proven byte-identical before an edit — which now covers every line ng
  runs, not half of them (spec §3).
- **The copied files keep production's integer widths** — no fork: a port that retypes `u32` to `Bp`
  under itself cannot be shown byte-identical. `Bp` appears at the generator's own boundary or not at
  all (spec §3).
- **Nothing is written into an observation that its read did not witness** — the haplotype builder
  does not fill, and `widen` extends the REF bucket only. The witnessed extent comes from the
  **events**, not the alignment span, which is blind to `N`, adaptor-masked, ref-skipped and
  dropped-indel positions (spec §4, §6).
- **The witnessed extent is stored in absolute reference coordinates; coverage is resolved once at
  `finalise()`** — bases are fixed at fold time, the footprint is not, so the two are kept on
  different clocks (spec §4).
- **`ReadCoverage` is `Complete` + one `Observed` run in locus positions** — `Complete` is kept as
  the common case and as `complete_observations()`'s cheap test; a non-contiguous witness yields no
  observation and is counted (spec §6).
- **`placed_left` is carried; `placed_start` is not** — `placed_left` feeds the read-position-bias
  penalty in `vcf/qual_refine.rs` that production subtracts from QUAL, so dropping it would forfeit
  QUAL parity; `placed_start` is consumed by no model and is cheap to re-add later as a `finalise()`
  derivation from the span already stored (spec §6).
- **Chain ids: carried, never on the REF observation, as a type invariant** — absence *is* the
  encoding (the reference is the default haplotype and a default needs no tag), and coverage is
  encoded separately by locus existence and `reads_without_observation`, so nothing is conflated
  (spec §6).
- **A locus per covered position, REF-only included** — the per-sample/cohort split requires it
  (spec §2).
- **The chain-id allocator and both reference accessors live on the generator, across segments** —
  per-segment allocation collides ids between regions; a per-segment accessor re-pays the `.fai`
  parse. The allocator is nonetheless `reset()` at each region end, so only `next_id` crosses the
  boundary (spec §8).
- **The walk owns its read stream; `src/ng/read/input/` is changed to hand back an owned one
  (owner, 2026-07-27).** Production's walker owns its input too — ng only has to remove the borrow at
  its source, because ng's contract is *resumable* where production's is scoped (§2.2). *Rejected:*
  **bounded windows** inside this generator — it materialises a window's reads (contradicting spec
  §7's "no read buffer at all"), breaks chain continuity inside a long region, and adds a knob that
  exists only to satisfy a borrow rule; if the owned stream were somehow unreachable, the answer
  would be to revise read input's architecture, not to window the pileup. Also *rejected:* a lifetime
  on `LocusGenerator` (spreads `'a` to `GeneratorSlot`, `GeneratorSet` and the iterator —
  `locus_generation.md` §4 rejected exactly that), and production's own non-resumable shape (it means
  buffering a region's loci, which the one-at-a-time contract exists to prevent).
- **The read group joins the observation's identity, at `@RG` grain** — the bucket key becomes
  `(bases, read_coverage, read_group)`, so a per-group model gets the allele × group cross with its
  moments; free at one read group. **The near-term consumer is the STR path, not this one**: its
  stutter and `ε` are already fit per sample group, off groups it currently has to *infer* (spec §6).
  `ReadGroupId` is the finest grain, so library and experiment stay a downstream fold — exact,
  because every support field is additive and the merged cells share their key's other two
  components (spec §6, §11).
- **`reads_without_observation` is "considered minus folded", per record; the wholly-silent read is
  counted at run level** in `PileupGeneratorCounts::reads_silent_over_footprint` (§1.1). The
  per-locus value is therefore an honest **lower bound** — say so in its doc comment, since the
  shared type's own wording promises more (spec §11).
- **`LocusGenerationError` gains a `Walker` variant.** No fork: none of `TypedRegion` / `Reads` /
  `Reference` ([locus_generation/mod.rs:263](../../../../src/ng/locus_generation/mod.rs#L263))
  describes a malformed read or an exhausted chain-id space. ng's own `#[non_exhaustive]` enum, so
  the addition is source-compatible (spec §11).

## 4. Prerequisites and impl-time confirmations

*No open design questions remain; every one is resolved in spec §11.*

**One prerequisite pass, *before* this generator is coded**, on the model of the STR generator's
`bundle_threshold` rename: its own cargo-verified commit with its own review, not a drive-by inside
this work. It is not a design question — it is here because the plan must sequence it first.

- **P1 — `src/ng/read/input/` returns an owned region stream** (§2.2, §3). Two borrows become `Arc`s
  in a built, reviewed ng module. Regression anchor already exists: the BAM/CRAM parity oracle and
  the read-input suite must pass unchanged, since the change is representational and no read moves.
  *(A second prerequisite — a `read_group` field on production's frozen `PreparedRead`, preceded by
  relocating `ReadGroupId` to a crate-visible home — was drafted and is **withdrawn**: ng copies the
  type instead, so production is edited zero times.)*

- *Impl-time confirmations, not design items:*
  - **`MultiChromRefFetcher::fetch` allocates a `Vec<u8>` per call** where `RefSeq::fetch_into`
    writes into a caller buffer. It is called per `open_new`/`widen`, not per read, so the shim ships
    as-is; if it profiles, ng's copy of `open_record.rs` can take a `fetch_into`-shaped accessor,
    since that file is ng's to change.
  - **Where `ReadPreparer::prepare_read` sits in the chain** — as a `map` on the owned stream
    feeding `PileupWalker`, so preparation stays lazy and per read. It must not become a collect
    step; that is the property spec §7 names as the one a port can quietly destroy.
  - **No spec fold-in is owed here:** spec §2's "one read query per segment, lazily streamed" and
    §7's "no read buffer at all" are both literally true of §2.2's shape.

## 5. Reconciliation with existing code

Every row read at the cited line (2026-07-27). Convergence, not new types.

| what | existing code | ng action |
|---|---|---|
| the walk loop | `driver.rs` `run`/`PileupWalker` [:41](../../../../src/pileup/walker/driver.rs#L41), [:51](../../../../src/pileup/walker/driver.rs#L51) | **copy → `genome_walk.rs`**, then change (the region clamp; `I: IntoIterator` is kept, which is what lets `RegionWalk` own its stream) |
| the fold | `open_record::process_position` [:654](../../../../src/pileup/walker/open_record.rs#L654) | **copy**, then change (§1.2) |
| CIGAR cursor + adaptor mask | `CigarCursor` [cigar_cursor.rs:143](../../../../src/pileup/walker/cigar_cursor.rs#L143) | **copy** — `events_at`/`events_overlapping` take `&PreparedRead`, so ng's read type reaches it |
| per-position events | `ReadEvent` [decompose.rs:15](../../../../src/pileup/walker/decompose.rs#L15) | **copy** — `decompose` takes `&PreparedRead` [:79](../../../../src/pileup/walker/decompose.rs#L79); brings its `#[cfg(test)]` oracle, which the cursor is parity-tested against |
| active read set | `ActiveReads` [active_read_set.rs:32](../../../../src/pileup/walker/active_read_set.rs#L32) | **copy** — `admit` takes `PreparedRead` by value [:93](../../../../src/pileup/walker/active_read_set.rs#L93) |
| chain-id allocation + mate pairing | `ChainIdAllocator` [chain_id_allocator.rs:78](../../../../src/pileup/walker/chain_id_allocator.rs#L78) | **copy** — `allocate_for_read` takes `&PreparedRead` [:205](../../../../src/pileup/walker/chain_id_allocator.rs#L205); one instance per sample stream |
| walker config + errors + summary | `WalkerConfig` [walker/mod.rs:133](../../../../src/pileup/walker/mod.rs#L133), `WalkerError` [errors.rs:12](../../../../src/pileup/walker/errors.rs#L12), `RunSummary` [driver.rs:228](../../../../src/pileup/walker/driver.rs#L228) | ng mints its own (`PileupGeneratorConfig`/`Counts`, §1.1); production's `DEFAULT_*` constants are reused by name |
| `CigarOp` | [walker/mod.rs:43](../../../../src/pileup/walker/mod.rs#L43) | reuse as-is — already `pub`, unchanged by ng, already used across `src/ng/` |
| the stage-1 comparison type | `PileupRecord` [pileup_record.rs:174](../../../../src/pileup_record.rs#L174) | compare against; `PartialEq` [:208](../../../../src/pileup_record.rs#L208) is bitwise on the two `f32`s, which makes the `NaN` placeholders comparable |
| `ChainId` | [pileup_record.rs:30](../../../../src/pileup_record.rs#L30) (`= u64`) | reuse as-is |
| the locus type ng fills | `SampleLocusObservations` [locus_generation/mod.rs:34](../../../../src/ng/locus_generation/mod.rs#L34), `ObservedSequence` [:113](../../../../src/ng/locus_generation/mod.rs#L113) | fill; `kind` is `LocusKind::Generic` [:164](../../../../src/ng/locus_generation/mod.rs#L164) |
| `ReadCoverage` | [locus_generation/mod.rs:146](../../../../src/ng/locus_generation/mod.rs#L146) | **reshape** to `Complete` + `Observed { offset_in_locus, positions_covered }` (spec §6). Not `#[non_exhaustive]`; **six** exhaustive match sites — `num_obs_along_locus` [:81-83](../../../../src/ng/locus_generation/mod.rs#L81), the STR tally [ssr.rs:1015](../../../../src/ng/locus_generation/ssr.rs#L1015) and sort key [:1072](../../../../src/ng/locus_generation/ssr.rs#L1072), plus four dump tools — and four STR minting sites, two of which pass the variant **as a function value** ([ssr.rs:713](../../../../src/ng/locus_generation/ssr.rs#L713), `:716`) and so must be restructured |
| the read-group id | `ReadGroupId` [ng/types.rs:178](../../../../src/ng/types.rs#L178); the run's table maps it to `ReadGroup { library, experiment }` [read_groups.rs:43](../../../../src/ng/read/input/read_groups.rs#L43) | reuse as-is — it stays an ng type, which only works because ng owns the read type too |
| the prepared read | `PreparedRead` / `MateRole` / `ReadLengthError` [walker/mod.rs:236](../../../../src/pileup/walker/mod.rs#L236) | **copy into `src/ng/read/prepared_read.rs`**, extend with `read_group` (§1.2). Reverses [read_preparation.md](../spec/read_preparation.md) §3's reuse-as-is, which owes a fold-in |
| `ObservedSequence` | [locus_generation/mod.rs:113](../../../../src/ng/locus_generation/mod.rs#L113) | **extend** with the read group **and `placed_left`** (§3 — the type carries neither bias field today, so without it `finalise` has nowhere to put the value); a shared-type change, so the STR generator splits too — a fold-in for [locus_generation.md](locus_generation.md) and its arch. Shares one fixture rebaseline with the `ReadCoverage` reshape above |
| the generator contract + slot | `LocusGenerator` [:198](../../../../src/ng/locus_generation/mod.rs#L198), `GeneratorSet.generic: GeneratorSlot<()>` [:375](../../../../src/ng/locus_generation/mod.rs#L375) | implement; replaces `NoLoci { NotImplemented }` |
| read access | `SampleReads::reads_in_region` (`read/input/mod.rs`) | reuse — after the §2.2 prerequisite makes the stream owned. **`ng-read-groups` keeps the borrow**, so the prerequisite stands as written |
| the read type the stream yields | `AlignedRead` [read/aligned_read.rs:36](../../../../src/ng/read/aligned_read.rs#L36) — replaces `MappedRead` as `SampleRegionReads::Item` [mod.rs:590](../../../../src/ng/read/input/mod.rs#L590) | consume as-is; its `read_group: ReadGroupId` is what the preparer threads into ng's own `PreparedRead` (§1.2) |
| `AlignedRead` → `PreparedRead` | `ReadPreparer::prepare_read` [read/mod.rs:113](../../../../src/ng/read/mod.rs#L113) (takes `AlignedRead` since the merge; was `MappedRead`), `LeftAlignPreparer` [left_align.rs:87](../../../../src/ng/read/left_align.rs#L87) | **call** — this generator is step 2's only consumer |
| the reference | `RefSeq` [ng/ref_seq.rs:142](../../../../src/ng/ref_seq.rs#L142) → `MultiChromRefFetcher` [fasta/mod.rs:114](../../../../src/fasta/mod.rs#L114) | shim (§1.3); contracts already agree on canonical bytes |
| the region clamp | `drive_region_into_writer` [pileup_to_psp.rs:271](../../../../src/pileup/per_sample/pileup_to_psp.rs#L271) | reuse the **rule** (`(start..=end).contains(&record.pos)`), not the code — the rest of that file is `.psp` machinery |
| the parity harness shape | [alignment/delimit_parity.rs](../../../../src/ng/alignment/delimit_parity.rs), [read/left_align_parity.rs](../../../../src/ng/read/left_align_parity.rs) | model for `parity.rs` |
| walker test fixtures | `MockFasta`, `snp_read`, `paired_snp_reads` [walker/tests.rs:32-152](../../../../src/pileup/walker/tests.rs#L25) | `pub(crate)` under `#[cfg(test)]` — usable from ng's tests today |
| the windowed-statistics buffer | `SampleSummaryAccumulators` [pileup_to_psp.rs:57](../../../../src/pileup/per_sample/pileup_to_psp.rs#L57) | **do not port** — it exists to fill two `.psp`-only fields ng's type does not have |

## Test & bench shape

Tests live beside the code; `parity.rs` is `#[cfg(test)]`, the `delimit_parity` shape.

**Stage 1 is a gate, not a permanent test.** The verbatim copy is proven to emit `PileupRecord`
streams equal to production's on one shared `PreparedRead` stream, plus `RunSummary` — and that
harness dies when §1.2 lands, because the port then *deliberately* differs. It must therefore be
shown to discriminate before it is retired: mutate mate overlap, adaptor masking, widening, the
re-fold and the column cap in turn, and watch it fail (spec §12).

**What survives is a narrower permanent differential**: on loci where every folded read spans the
final footprint, ng and production must agree **forever** — that set is untouched by every change in
§1.2, so it is a real regression anchor rather than a snapshot. Loci outside it are the divergence,
enumerated and counted, and that count is the deliverable (spec §12).

**The dump tool** `examples/ng_generic_loci_dump.rs` is the definition of done, asserted on a
committed fixture — including a deletion long enough to widen a record past an already-expired read,
which is the only way to check §1.2's widening rule and the one thing production cannot pass by
construction (spec §12).
