# ng — the alignment cursor: types & interfaces

*Architecture draft, 2026-08-02. Code-facing companion to
[`../spec/alignment_cursor.md`](../spec/alignment_cursor.md) — every **why** points there and is
not re-argued here. Replaces the region-query half of
[`alignment_file.md`](alignment_file.md); that doc's §3.3 carries the supersession note. Under
[`module_layout.md`](module_layout.md). Build order:
[`../impl_plan/alignment_cursor.md`](../impl_plan/alignment_cursor.md).
Naming: **STR** in prose, `ssr` in code.*

---

## Module home

`src/ng/read/input/`, the module that already owns everything whose subject is one alignment
file. Two new files and one folder:

```
src/ng/read/input/
    open_bam.rs           AlignmentFile — loses the pool, gains cursor() and contigs()
    cursor.rs             AlignmentCursor, CursorError
    record_reader/
        mod.rs            enum RecordReader + the kept-reads contract both arms honour
        bam.rs            BamRecordReader
        cram.rs           CramRecordReader  (absorbs DecodedContainer)
        in_memory.rs      InMemoryRecordReader
    region_query.rs       deleted — its parts move to record_reader/ and cursor.rs
    reference.rs          unchanged (the per-chromosome registry is deferred, spec §12)
```

`record_reader/` is a folder rather than a file because the two real arms are large and share
only a contract, and because a third arm exists for tests. It is **not** a trait-per-step
bake-off: the set is closed and the arms do not compete (spec §5).

## 1. The types

### 1.1 The file — read-only once open

Nothing on `AlignmentFile` is written after construction, so it needs no locks and no atomics
(spec §3). It loses `readers_opened`, the reader pool, `ReaderHandle` and `BorrowedReader`.

```rust
pub struct AlignmentFile {
    path: Arc<Path>,
    header: Arc<sam::Header>,
    index: AlignmentIndex,
    resolution: ReadGroupResolution,
    /// CRAM only. **Deferred — see spec §12.** For now a cursor takes its
    /// chromosome's bases once at construction through the existing
    /// `OpenReference::bases_for_contig` (`reference.rs:238-252`); no new type.
}
```

### 1.2 The cursor — one consumer, one chromosome

```rust
/// A reader positioned in one chromosome of one file, holding the reads it has
/// recently decoded so a nearby region can be answered without unpacking again.
///
/// Not `Sync`: an open file position belongs to one consumer. Parallelism comes
/// from more cursors, never from sharing one (spec §3).
pub struct AlignmentCursor<R: RawRefSeq> {
    /// The whole chain below, owned: the filter holds `RegionRecords`, which
    /// holds the `RecordReader`. §2.3 shows the layering and why it is not a
    /// cycle. The reference accessor the mismatch filter needs sits inside,
    /// taken once here rather than per query (perf review L2).
    filter: ReadFilter<RegionRecords, R>,
    /// **Our** reads — decoded and filtered — kept for the next region, so a
    /// read is transformed once rather than once per region that returns it
    /// (spec §5). Drained in position order, extended by reading on.
    kept: VecDeque<AlignedRead>,
    /// The start of the last region served. This is the entire forget rule
    /// (spec §6): reuse `kept` when the next region starts at or after it.
    last_region_start: Option<u64>,
    contig: ContigId,
    /// CRAM only, taken once at construction; `fasta::Repository` is the
    /// existing type (`reference.rs:130`). Nothing new is minted (spec §10).
    bases: Option<fasta::Repository>,
}
```

### 1.3 Where records come from

One variant per format; the cursor interprets for both (spec §5). Each arm owns its own index
walk and its own unpacking. **Neither keeps anything across regions** — what is kept is *our*
reads, decoded and filtered, held one level up in the cursor, so a read is transformed once rather
than once per region that returns it (spec §5). A `RecordReader` therefore holds only its position
and a one-record pushback, for the record the sorted early stop consumes without yielding.

```rust
/// Finds records with the index, unpacks them, keeps the recent ones, and hands
/// over the next on demand.
enum RecordReader {
    Bam(BamRecordReader),
    Cram(CramRecordReader),
    /// A fixed list. Permanent, not test-only: it gives the tests and the
    /// differential harness a reader with no file behind it.
    InMemory(InMemoryRecordReader),
}
```

**The contract every arm honours**, stated once in `record_reader/mod.rs` so the two cannot
drift:

- `begin(region)` — position for a region on the cursor's chromosome. Reuses what is kept when
  the region's start is at or after the cut-off; otherwise drops and jumps.
- `next(&mut RecordBuf) -> io::Result<bool>` — the next record from the file, **in position
  order**. A block is unpacked only when the records before it have been handed over (spec §4).
- **Forgetting is one comparison, shared by both arms**: keep what is kept only when the new
  region starts at or after the last one served; evict a record once it ends before the current
  region's start. No index lookup — an earlier draft derived a byte cut-off from the index and was
  unsound three ways (spec §6).
- Records are yielded **raw**: no contig test, no overlap test, no read-group resolution. Those
  belong to the cursor.

### 1.4 Errors

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CursorError {
    /// A region on a chromosome this cursor does not cover. Make a cursor for
    /// that chromosome. A guard against a caller bug: correct code compares
    /// against `contig()` first and never sees it.
    #[error("cursor covers {cursor:?} but region is on {requested:?}")]
    WrongChromosome { cursor: ContigId, requested: ContigId },
    #[error("reading {path}")]
    Io { path: Arc<Path>, source: std::io::Error },
}
```

No ordering variant: within a chromosome there is no ordering rule (spec §4).

## 2. The interfaces

### 2.1 Making a cursor

```rust
impl AlignmentFile {
    /// The chromosomes this file may be read over — the reference's own list,
    /// proven equal to it by the open gate (spec §4).
    pub fn contigs(&self) -> &ContigList;

    /// One cursor for one chromosome. Opens its own descriptor and, for CRAM,
    /// takes a share of that chromosome's bases.
    ///
    /// Called once per chromosome per worker. Never per region — that is the
    /// whole point (spec §1).
    pub fn cursor<R: RawRefSeq>(
        self: &Arc<Self>,
        contig: ContigId,
        reference: R,
    ) -> Result<AlignmentCursor<R>, AlignmentFileError>;
}
```

**Contract.** Fallible only at construction: opening the descriptor, and for CRAM taking or
reading the bases. After it returns, a cursor cannot fail to exist.

### 2.2 Asking for reads

```rust
impl<R: RawRefSeq> AlignmentCursor<R> {
    /// The chromosome this cursor covers.
    pub fn contig(&self) -> ContigId;

    /// Point the cursor at `region`, and forget nothing that is still reachable.
    ///
    /// Any region on this chromosome, in any order — always correct. Free when
    /// the region already lies inside what the cursor holds; otherwise a jump.
    pub fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), CursorError>;

    /// The next read of the current region. `None` at the end of it.
    ///
    /// Owned, not filled into a caller's buffer: the reuse trick lives at the
    /// record seam below (`RecordSource::read_next`), not here — `decode`
    /// returns an owned `AlignedRead` (`filtering.rs:348`). Shaped as
    /// `Iterator::next` so `Iterator` stays available (spec §3).
    pub fn next_read(&mut self) -> Option<Result<AlignedRead, CursorError>>;
}
```

**Contract.**

- **`next_read` yields exactly the reads overlapping the current region**, after step-1 filtering,
  in position order. Choosing a wider region than you mean to walk is the caller's business, not
  this type's (spec §2).
- **Nothing is unpacked ahead of demand.** A block is decompressed only when a `next_read` needs a
  record inside it, so a caller that pulls one read and moves on has unpacked at most one block.
- **Memory is bounded by depth *plus the kept window*.** The walker's own bound — depth, never
  region size (`locus_generation_pileup.md:72-80`) — is about a walker that retains nothing, so it
  does not carry over unqualified. A cursor adds what it keeps: the overlap between two
  consecutive regions, ~5,000 bases' worth with this caller (spec §10, arithmetic not measured).
- **`move_to_region` leaves the cursor usable on any outcome.** `WrongChromosome` means "make a
  cursor for that chromosome"; the cursor itself is unharmed and still good for its own (spec §4).
- **Abandoning a region is free.** A caller that stops pulling and moves elsewhere leaves nothing
  to unwind: there is no stream object, so there is nothing to give back and nothing to forget.

**Why there is no stream type.** Earlier drafts returned one per region. It had to own or borrow
the cursor, and every difficulty — the lifetime, giving the cursor back, a caller forgetting to —
came from that. Removing it removed all of them (spec §3).

### 2.3 How the cursor holds the filter without a cycle

A cursor yields `AlignedRead`, so it must own a `ReadFilter`. A `ReadFilter` owning the cursor
would be a cycle — but it does not have to. The filter's source is the layer *below* the cursor,
not the cursor itself:

```text
RecordReader   raw records, from what is kept or from the file      (§1.3)
   ↓
RegionRecords  this region's records only: contig test, overlap test, sorted
               early stop, read-group resolution, per-read-group tallies
   ↓
ReadFilter     step-1 filtering                                      (filtering.rs)
   ↓
AlignedRead    what the caller gets
```

```text
kept: VecDeque<AlignedRead>     what a replay serves — above the filter, so a
   ↑                            replayed read skips both decode and filtering
ReadFilter<RegionRecords, R>    step-1 filtering; its counts() is the tally
   ↑
RegionRecords                   this region only; owns the RecordReader
   ↑
RecordReader                    finds and unpacks; keeps nothing across regions
```

`move_to_region` reaches through the filter to reposition the layer below:
`self.filter.source_mut().move_to(region)`.

**The tallies need no field and no hand-over.** `ReadFilter::counts()` is already a running tally,
readable at any point (`filtering.rs:797-799`), and the filter now lives as long as the cursor. So
a caller reads `cursor.counts()` whenever it likes — which is a better answer than spec §3's
"collect them when the cursor is retired", and shrinks the debt recorded there.

**One addition to an existing type:** `ReadFilter::source_mut(&mut self) -> &mut S`. Today the
source can only be reached by consuming the filter (`into_parts`, `filtering.rs:949`), which is
what forced a new filter per region. `OrderVerified` wraps the filter as it does now.

`RegionRecords` is today's `BamRegionSource`/`CramRegionSource` with the two format-specific
halves lifted out into `RecordReader` and the region narrowing kept — so this is a
re-layering of code that exists, not new machinery.

### 2.4 The sample layer, which is what callers hold

A sample is k files (k = 1 for most, more when one sample was sequenced in several runs), and
both generators ask a *sample* for reads. So the type a caller holds is a sample cursor:

```rust
/// One sample's reads over one chromosome. Holds a cursor per file and merges
/// them, which is what `SampleReads::reads_in_region` does per region today.
pub struct SampleCursor<R: RawRefSeq> {
    files: Vec<AlignmentCursor<R>>,
    /// The argmin merge's per-file head, as `MergedRegionReads` keeps today.
    /// Cleared by `move_to_region`, because the heads belong to a region.
    heads: Vec<Option<AlignedRead>>,
}

impl<R: RawRefSeq> SampleCursor<R> {
    pub fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), CursorError>;
    pub fn next_read(&mut self) -> Option<Result<AlignedRead, CursorError>>;
}
```

**The single-file case must stay free.** `SampleRegionReads` is deliberately an enum with a
`Single` arm rather than a merge of one, because dynamic dispatch on the hottest loop in the
module would stop the per-read chain being inlined (`mod.rs:560-566`). A sample cursor must keep
that property.

**Open descriptors.** A cursor holds its file open for as long as it lives, and there is no pool
bounding the count. The bound is `files × generators × workers` — with one file per sample, two
generators and 16 workers, 32; with ten files per sample, 320, which is close to macOS's default
soft limit of 256. This needs a stated ceiling before the fan-out lands, and it is a real cost
against the memory-efficiency thesis rather than a free consequence of removing the pool.

## 3. Design decisions — decided

Each records the code shape; the argument is the spec's.

1. **The caller keeps the cursor and pulls into its own buffer** — `move_to_region` then
   `next_read()`, no stream type. Two rejected: a callback per read, which would make the
   walker store every locus it produces; and a per-region stream object, which was the draft this
   replaced — spec §3, §13.1.
2. **No pool, no lock, no counters on the file.** Rejected: `&self` + an internally-shared pool,
   which existed to keep parallelism cheap to retrofit; per-worker cursors serve that with
   nothing shared — spec §3, and `alignment_file.md` §3.3's supersession note.
3. **Enum over trait for `RecordReader`.** The set is closed, and a trait adds a type parameter
   through four layers — spec §5.
4. **A cursor covers one chromosome; region order within it is unrestricted.** Two restrictions
   were weighed and only this one kept: a chromosome change costs a chromosome, a backward jump
   costs a block — spec §4, §13.4.
5. **The forget rule is a coordinate comparison**, shared by both arms, with no index lookup and
   nothing to tune. Rejected: deriving a byte cut-off from the index — `min_offset` collapses to
   byte 0 past the last populated window on both index kinds, and a byte range is not a record set
   — spec §6.
6. **The cursor owns its reference accessor**, taken at construction instead of per query. This
   is perf-review finding L2, and it **removes `MakeReference` from `PileupGenerator`'s type
   parameters** — three become two.
7. **`RecordReader` yields raw records; the cursor classifies.** One copy of the contig test,
   overlap test, early stop and read-group resolution for both formats, so a replayed record and
   a freshly decoded one cannot diverge — spec §5.

## 4. Reconciliation with existing code

Every line number below was checked by opening the file, against the code as committed at
`d95ce8b` with none of the experimental edits applied. No row is left marked "roughly right" —
that only moves the checking to whoever implements this.

| new name | existing code | what happens |
|---|---|---|
| `AlignmentCursor` | `ReaderHandle` (`open_bam.rs:166`), `BorrowedReader` (`:196`), `ReaderKind` (`:184`) | replaced; the pool and the borrow-and-return dance go |
| `RecordReader::Bam` | `BamRegionSource` (`region_query.rs:70`) | becomes it, minus the classification, plus kept records and the cut-off |
| `RecordReader::Cram` | `CramRegionSource` (`region_query.rs:570`), `DecodedContainer` (`:333`) | becomes it; the container cache stops being CRAM-only in *spirit* and stays CRAM-shaped in *fact* |
| `RecordReader` (enum) | `RegionSource` (`region_query.rs:950`) | same shape, new job: sourcing only |
| the cursor's classify step | steps 4–5 of `BamRegionSource::read_next` and the `RecordOwner` match in `CramRegionSource::read_next` | lifted out of both, written once |
| *(no replacement)* | `RegionReads<R>` (`open_bam.rs:248`) | **deleted**; its `Drop`-returns-to-pool goes with it |
| the filter seam | `trait RecordSource` (`filtering.rs:365`), `ReadFilter` (`:800`, `:949`) | **changes by one accessor** — see §2.3. An earlier row claimed "unchanged", which was wrong: `ReadFilter` owns its source, exposes no `source_mut` (zero occurrences), and returns it only by being consumed |
| the reference accessor | `reads_in_region<R>(…, reference: R)` (`open_bam.rs:489-493`) | moves to `cursor()`; `SampleReads::reads_in_region`'s `make_reference: F` factory (`mod.rs:541-549`) goes with it |
| *(none — deferred)* | `OpenReference` (`reference.rs:116-143`), `bases_for_contig` (`:238-252`) | **unchanged.** A cursor calls `bases_for_contig` once at construction and holds the handle; the existing one-contig bound is correct for every caller that exists today. The per-chromosome registry is deferred to spec §12 |
| `ContigList`, `ContigId`, `GenomeRegion` | `types.rs`, `reference_info.rs` | used as they are; nothing new minted |
| **the layer callers actually use** | `SampleReads::reads_in_region` (`mod.rs:541`), `SampleRegionReads` (`mod.rs:566-571`), `MergedRegionReads` (`merge.rs`) | **⚠ an earlier draft omitted this entirely.** Neither generator calls `AlignmentFile` — the generic one calls `SampleReads::reads_in_region` (`generator.rs:853`) and so does the STR one (`ssr.rs:375`). A sample spans k files, so a sample-level cursor is required: it holds k file cursors, re-points all of them on `move_to_region`, and merges their output. See §2.4 |

**Not a duplicate check:** `CursorError` is the only new error type, and it names two conditions
neither `AlignmentFileError` nor `IngestError` expresses.

## 5. Open items

- `OPEN:` **Must a replayed record be staged through a scratch buffer before decoding?** Worth
  +2.8 % of wall time. Decide after the first working cursor is measured end to end — spec §13.
- *Impl-time confirmation:* the exact `RecordReader` contract method names. The contract is fixed
  (§1.3); the spelling is the implementer's.

## 6. Test & bench shape

Tests live beside the code they cover, plus one cross-format oracle in `record_reader/mod.rs`.
The regression anchors are the ones that caught the first attempt, which **1,471 unit tests did
not** (spec §11):

- **`ng_generic_walk_probe`, real data, identical output.** chromosome 21 must print
  `loci=236081 observations=251786 reads_admitted=54709`; chromosome 1 `loci=1541788
  observations=1647161`.
- **The oracle, extended from one query to a sequence.** Today
  `t5_the_indexed_query_returns_exactly_what_a_linear_scan_returns` drives a *single* query;
  it must drive a run of regions through one cursor. This is the missing test.
- **Counters, asserted.** Records decoded must approach the true count (34,633 on chromosome 21),
  not a multiple of it — the saving measured rather than assumed.
- **Mutation.** Delete the cut-off, and the reuse check; each must fail something.

`benches/ng_generic_pileup_perf.rs` gains no new case: the cursor is measured through the walk it
exists to speed up, and the region-grain axis it already sweeps is the one that moves.
