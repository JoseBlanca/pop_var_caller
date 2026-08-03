# ng — read filtering in stages: types & interfaces

*Architecture draft, 2026-08-03. Code-facing companion to
[`../spec/read_filtering_stages.md`](../spec/read_filtering_stages.md) — every **why** points
there and is not re-argued here. Revises the single-type shape in
[`read_filtering.md`](read_filtering.md) §1; the filters that doc pins down are unchanged.
Under [`module_layout.md`](module_layout.md). Build order:
[`../impl_plan/read_filtering_stages.md`](../impl_plan/read_filtering_stages.md). Naming per
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md). Signatures are
illustrative; the **contract** is the deliverable.*

**This change adds no types.** It renames several, deletes three, and moves one loop. Read it as
a subtraction — spec §5.

**Both of the spec's questions are settled** (spec §9): the cursor holds the reference bases and
the buffer they are read into, and both filters stay plain functions; and the up-front contig
check becomes a comparison of two contig tables.

---

## 1. Module homes

| module | holds | change |
|---|---|---|
| `read/aligned_read.rs` | `RawAlignedRead`, `NoodlesRawAlignedRead`, `AlignedRead`, `decode_record` | **gains** the two raw types from `filtering.rs` — they are one thing in two states, and the conversion is already here |
| `read/filtering.rs` | `ReadFilterConfig`, `DropReason`, `FilterVerdict`, `ReadFilterCounts`, the two verdicts | **only the keep-or-drop rules and their thresholds** — loses the raw types, the loop, the source trait and `ReadFilter` |
| `read/input/aligned_reads_reader/` | the per-format readers | renamed from `record_reader/` |
| `read/input/region_raw_aligned_reads.rs` | `RegionRawAlignedReads` | renamed from `region_records.rs`; its trait impl becomes inherent methods |
| `read/input/cursor.rs` | `AlignmentCursor` | **gains** the loop, the reference, the buffer and the tally |

## 2. The renames

Mechanical, and the reasoning is spec §6. Listed here because every signature below uses the new
names.

| now | becomes |
|---|---|
| `RawRecord` (trait) | `RawAlignedRead` |
| `NoodlesRawRecord` | `NoodlesRawAlignedRead` |
| `RecordReader` | `AlignedReadsReader` |
| `BamRecordReader` / `CramRecordReader` / `InMemoryRecordReader` | `BamAlignedReadsReader` / `CramAlignedReadsReader` / `InMemoryAlignedReadsReader` |
| `record_reader/` | `aligned_reads_reader/` |
| `RegionRecords` | `RegionRawAlignedReads` |
| `DecodedContainer::fill_record` | `fill_raw_read` |
| `RecordIndex` (private) | `RawReadIndex` |

**The readers do not carry "raw" in their names**, so each reader's doc comment must state that
what it yields is undecoded — the name no longer says it.

## 3. The pieces

### 3.1 The item, and the conversion between its two states

```rust
// read/aligned_read.rs

/// One alignment record as it comes off the file, undecoded: a flag and a mapping
/// quality readable without unpacking anything else.
///
/// **An unmapped read is one of these.** Filter #5 rejects it, so it exists here
/// before it is rejected. `AlignedRead` has no such case — the conversion refuses a
/// record with no reference id or no position (spec §4).
pub trait RawAlignedRead {
    fn flag(&self) -> u16;
    fn mapq(&self) -> MapQual;
    /// Which read group this read belongs to, once its reader has resolved it.
    /// Read by the **tally**, not by any filter: a drop is charged to the read group
    /// it came from, and the first filter runs before any `AlignedRead` exists.
    fn read_group(&self) -> Option<ReadGroupId>;
    /// Convert into ng's own read. Fatal on failure, never a drop.
    fn decode(&self) -> io::Result<AlignedRead>;
}
```

Unchanged apart from the name and the module. `NoodlesRawAlignedRead` moves with it.

### 3.2 The two verdicts

Both are plain functions today — not methods on anything — and both stay that way (spec §9 Q1).
`filtering.rs` exports them; the cursor calls them.

```rust
// read/filtering.rs — the keep-or-drop rules; nothing here reads a file or converts.

/// Filters #1–#6, on the flag and the mapping quality. **Needs no reference**, so
/// anything wanting flag/MAPQ filtering without one can call it (spec §5).
pub(crate) fn verdict_on_raw_read(
    flag: u16, mapq: MapQual, config: &ReadFilterConfig,
) -> FilterVerdict;

/// Filters #7, #9 and #8, **in that order** — cheap checks first so a doomed read
/// never pays the reference fetch, and a read failing both #9 and #8 is charged to
/// the root cause. `Err` is a reference-fetch failure and is fatal, never a drop.
///
/// `ref_buf` is caller-owned scratch, reused across reads.
pub(crate) fn verdict_on_aligned_read(
    read: &AlignedRead, reference: &impl RawRefSeq,
    config: &ReadFilterConfig, ref_buf: &mut Vec<u8>,
) -> Result<FilterVerdict, RefSeqError>;
```

**Contract.** Both are pure functions of their arguments — no I/O, no state, no tally. Neither
knows what a BAM is. The first cannot fail.

Settled (spec §9 Q1): the cursor holds `reference` and `ref_buf` and passes them in, so neither
filter becomes a type.

### 3.3 The region narrowing

`RegionRawAlignedReads` loses its `RecordSource` impl and keeps the same four methods as
inherent ones. **No signature changes** beyond the rename.

```rust
// read/input/region_raw_aligned_reads.rs
impl RegionRawAlignedReads {
    pub(crate) fn header(&self) -> &sam::Header;
    /// Fill `buf` with this region's next raw aligned read. `Ok(false)` = the region
    /// is done — which is *not* the end of the file, and the caller knows which.
    pub(crate) fn read_next(&mut self, buf: &mut NoodlesRawAlignedRead) -> io::Result<bool>;
    pub(crate) fn other_sample_reads(&self) -> u64;
    pub(crate) fn jump_to(&mut self, region: GenomeRegion) -> io::Result<()>;
    pub(crate) fn continue_into(&mut self, region: GenomeRegion);
}
```

**Why the trait goes:** one production impl, and after §3.4 no generic consumer — spec §6.

### 3.4 The cursor owns the loop

```rust
// read/input/cursor.rs
pub struct AlignmentCursor<R: RawRefSeq> {
    reads: RegionRawAlignedReads,
    /// The single buffer reused across the whole walk.
    buffer: NoodlesRawAlignedRead,
    /// Held here because only the second filter needs it, and it is the cursor that
    /// keeps it alive for the chromosome (spec §9 Q1).
    reference: R,
    ref_buf: Vec<u8>,
    config: ReadFilterConfig,
    /// One entry per read group met, in first-seen order. Cumulative until
    /// `reset_counts` — spec §7.
    counts: Vec<ReadGroupCounts>,
    /// **Replaces `FilterState`.** The cursor causes region ends, so it never has to
    /// ask why reading stopped — which is what makes the three-way state
    /// unnecessary (spec §5).
    failed: bool,
    kept: VecDeque<AlignedRead>,
    /* …the retention fields, unchanged… */
}
```

The loop, in `next_read`, after the kept set is exhausted:

```
read_next → false: the region is done, return None
          → Err:   failed = true, return the error
verdict_on_raw_read       → Drop: charge it, continue
decode                    → Err:  failed = true, return the error
verdict_on_aligned_read   → Drop: charge it, continue
                          → Err:  failed = true, return the error
                          → Keep: tally, push onto `kept`, emit
```

**Contract, unchanged from today's composite.** Lazy; one raw read resident. A fatal error is
yielded once and then the cursor refuses every later region (`CursorError::AfterFailure`,
already present). The order guard still runs on emit.

**What the cursor gains publicly:**

```rust
impl<R: RawRefSeq> AlignmentCursor<R> {
    /// Step-1's tally, one entry per read group met. Cumulative — spec §7.
    pub fn read_group_counts(&self) -> Vec<ReadGroupCounts>;   // exists
    /// Start a fresh tally window. The caller chooses the window; the cursor does
    /// not reset on its own, and never per region — spec §7.
    pub fn reset_counts(&mut self);                            // new
}
```

## 4. Errors

**No new error type and no change of meaning.** `ReadFilterError`'s three variants already name
the three pieces:

| variant | raised by |
|---|---|
| `Source` | `RegionRawAlignedReads::read_next` |
| `Decode` | the conversion |
| `Reference` | the second filter's mismatch check |

`verdict_on_raw_read` cannot fail, which is why it returns a bare verdict.

## 5. Design decisions — decided

- **The two filters and the conversion sit below the cursor's kept set; the cursor owns the
  loop.** The rule that stops the pieces being rearranged — spec §3.
- **Two filters, not one** — spec §2.
- **A filter takes a borrow and returns a verdict, never a read** — spec §4.
- **The verdict carries the drop reason**, so an `Option` return is ruled out: the tally is keyed
  on reason and read group.
- **The read group is read off the raw aligned read** — it is already stamped there for exactly
  this ([`filtering.rs:350-357`](../../../../src/ng/read/filtering.rs#L350)).
- **The second filter's #7 → #9 → #8 order is kept**, not the numbering order — spec §4.
- **`FilterState` is deleted for a `failed` flag** — spec §5.
- **The source trait is deleted and the in-memory reader gains a scripted error** — spec §6.
- **The tally lives on the cursor, cumulative, with `reset_counts`; not on `AlignmentFile`** —
  spec §7.
- **The cursor holds the reference bases and the buffer; neither filter becomes a type** —
  spec §9 Q1.
- **The up-front contig check becomes a contig-table comparison**, not ~2,580 window fetches, and
  it proves more than the fetches did — spec §9 Q2.
- **No trait, no bake-off:** no competing implementations, so plain functions and concrete types
  — `module_layout.md` principle 1a.

## 6. Reconciliation with existing code

Every row read at the cited line, 2026-08-03.

| what | existing code | action |
|---|---|---|
| the six flag/MAPQ filters | `verdict_pre_decode` [`filtering.rs:210`](../../../../src/ng/read/filtering.rs#L210) | **rename** to `verdict_on_raw_read`; body unchanged |
| the three conversion-dependent filters | `verdict_post_decode` [`filtering.rs:269`](../../../../src/ng/read/filtering.rs#L269) | **rename** to `verdict_on_aligned_read`; body unchanged |
| the conversion | `RawRecord::decode` [`filtering.rs:349`](../../../../src/ng/read/filtering.rs#L349) → `decode_record` [`aligned_read.rs:67`](../../../../src/ng/read/aligned_read.rs#L67) | **reuse as-is** |
| the loop | `ReadFilter::next` [`filtering.rs:895`](../../../../src/ng/read/filtering.rs#L895) | **move** into `AlignmentCursor::next_read`; `ReadFilter` deleted |
| the up-front contig check | `ReadFilter::new`'s per-contig fetch loop [`filtering.rs:688`](../../../../src/ng/read/filtering.rs#L688) | **replace** with `self.contigs.first_disagreement(reference.contigs())` in `AlignmentFile::cursor`, the same comparison the open gate makes [`open_bam.rs:206`](../../../../src/ng/read/input/open_bam.rs#L206) — spec §9 Q2 |
| the tally and its fold | `ReadFilterCounts` [`filtering.rs:122`](../../../../src/ng/read/filtering.rs#L122), `ReadGroupCounts` [`:661`](../../../../src/ng/read/filtering.rs#L661), `tally_for_current_record` [`:846`](../../../../src/ng/read/filtering.rs#L846), `counts` [`:868`](../../../../src/ng/read/filtering.rs#L868) | **move** to the cursor, including the `other_sample` rider on the first entry |
| the errors | `ReadFilterError` [`filtering.rs:578`](../../../../src/ng/read/filtering.rs#L578) | **reuse as-is** |
| the raw read | `RawRecord` [`filtering.rs:334`](../../../../src/ng/read/filtering.rs#L334), `NoodlesRawRecord` [`:479`](../../../../src/ng/read/filtering.rs#L479) | **rename and move** to `aligned_read.rs` |
| the region narrowing | `RegionRecords`, [now `region_raw_aligned_reads.rs`](../../../../src/ng/read/input/region_raw_aligned_reads.rs) | **renamed at A3**; trait impl → inherent methods at C3 |
| the per-format readers | `RecordReader` and arms, [now `aligned_reads_reader/mod.rs`](../../../../src/ng/read/input/aligned_reads_reader/mod.rs) | **renamed at A2**; `InMemory` arm gains a scripted error at C1 |
| the three-way stop | `FilterState` [`filtering.rs:646`](../../../../src/ng/read/filtering.rs#L646), `restart_after_end_of_input` [`:778`](../../../../src/ng/read/filtering.rs#L778), `has_failed` [`:794`](../../../../src/ng/read/filtering.rs#L794), `source_mut` [`:816`](../../../../src/ng/read/filtering.rs#L816) | **delete**, all four |
| the source trait and its doubles | `RecordSource` [`filtering.rs:366`](../../../../src/ng/read/filtering.rs#L366), `FakeSource`, `ErroringSource` | **delete** |
| the probe-free constructor and lent buffers | `with_validated_contigs` [`filtering.rs:746`](../../../../src/ng/read/filtering.rs#L746), `ReadFilterBuffers` | **delete** — no caller once the cursor owns the loop |
| the cursor's existing filter field and its four call sites | [`input/cursor.rs:241`](../../../../src/ng/read/input/cursor.rs#L241), `:387`, `:445`, `:448` | **replace** with the fields in §3.4 |

## 7. Open items

**No open design questions.** Both of the spec's are settled — §9 there keeps the reasoning.

**Impl-time confirmations, not decisions:**

- Whether `verdict_on_raw_read` takes `(flag, mapq)` or `&impl RawAlignedRead`. The split form is
  what exists; the whole-read form reads better beside its sibling and changes no contract.
- What shape the in-memory reader's scripted error takes — an error at position *n*, or an arm
  that always fails. The three tests it has to serve are named in spec §8.
- Whether `region_records.rs` is renamed on disk or the type simply moves. The spec assumes the
  file follows the type.
- The `+ ContigTable` bound the contig comparison needs on `cursor`'s `R`. Every accessor in the
  tree implements it, but it propagates to `SampleReads::cursor` and to both generators'
  signatures — mechanical, and worth doing in one commit of its own.

## 8. Test & bench shape

Tests stay beside the code. `filtering.rs`'s 45 tests split: the ones about a rule stay, the loop ones
move to `cursor.rs`, the three test-double ones are replaced (spec §8).

**The regression anchors are output identity on real data, not the unit suite** — spec §8 has the
four dumps and the walk probe's figures. Three tests that do not exist yet are named there too,
and each covers something no output comparison can see.
