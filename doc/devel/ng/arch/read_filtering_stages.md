# ng — read filtering in stages: types & interfaces

*Architecture. **The whole plan is built** — A (the renames), B (the contig check and the
`fill_raw_read` signature), C (the loop moves into the cursor) and D (the two tests output
identity cannot see). §6's table marks what landed; ⚠ blocks record where this document was wrong
and building it corrected them, rather than the sentence being quietly rewritten. Code-facing
companion to
[`../spec/read_filtering_stages.md`](../spec/read_filtering_stages.md) — every **why** points
there and is not re-argued here. Revises the single-type shape in
[`read_filtering.md`](read_filtering.md) §1; the filters that doc pins down are unchanged.
Under [`module_layout.md`](module_layout.md). Build order:
[`../impl_plan/read_filtering_stages.md`](../impl_plan/read_filtering_stages.md). Naming per
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md). Signatures are
illustrative; the **contract** is the deliverable.*

**This change adds no types.** It renames several, deletes three, and moves one loop. Read it as
a subtraction — spec §5.

> ⚠ **Two additions, and neither is a type** (recorded for the same reason spec §1's ⚠ block
> exists — the line is close enough to be worth naming). `AlignmentFileError` gained a **variant**
> at B1 and `ReadFilterError` a fourth at Checkpoint C; `CursorCounts` gained a **field**,
> `reads_converted`, at D2. No new concept enters the design in any of the three. The rule did
> bite once, and held: D2's review proposed a zero-sized witness making the conversion's position
> a compile-time property, and it was **rejected on this sentence** — §5 has the decision.

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
| `RecordIndex` (private) | `PackedReadEntry` |

**The readers do not carry "raw" in their names**, so each reader's doc comment must state that
what it yields is undecoded — the name no longer says it.

**Two rows were revised by the owner at Checkpoint A (2026-08-03), after the review of the
step that landed them.**

`RecordIndex` first became `RawReadIndex`, which failed twice: "Index" named the container's
*field* rather than one entry — the code around it already calls a single one an `entry` — and
`RawRead` was a second, shorter name for `RawAlignedRead`, minted inside the very milestone
that existed to make that vocabulary consistent. `PackedReadEntry` says what the value is: one
read in the packed form this container stores, which is the whole reason the type exists.

`fill_raw_read` **kept its name and gained a signature change — done at B2.** It took
`&mut RecordBuf`, so it filled only the record half of a raw aligned read and its one caller
stamped the read group on the following line. It takes `&mut NoodlesRawAlignedRead` and sets
both halves now, which makes the name true and puts "on CRAM the read group is decided at
decode" in one place instead of two. `DecodedContainer::read_group(i)` went with the split it
existed for. Not a rename, so it could not travel with Milestone A.

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

`RegionRawAlignedReads` loses its `RecordSource` impl and keeps its methods as inherent ones.
**As built at C3** — this listing was corrected at Checkpoint C; the notes below say where it was
wrong.

```rust
// read/input/region_raw_aligned_reads.rs
impl RegionRawAlignedReads {
    /// Fill `buf` with this region's next raw aligned read. `Ok(false)` = the region
    /// is done — which is *not* the end of the file, and the caller knows which,
    /// because the caller is what set the region.
    pub(crate) fn read_next(
        &mut self, buf: &mut NoodlesRawAlignedRead,
    ) -> Result<bool, RegionReadError>;
    pub(crate) fn other_sample_records(&self) -> u64;
    pub(crate) fn jump_to(&mut self, region: GenomeRegion) -> io::Result<()>;
    pub(crate) fn continue_into(&mut self, region: GenomeRegion);
}
```

**Why the trait goes:** one production impl, and after §3.4 no generic consumer — spec §6.

> ⚠ **Three corrections, all made by building it.**
>
> **`header()` is not there, and should not be.** It went at B1 with the contig probe that was its
> only caller; §6 already carried the instruction *"re-add it at C3 only if a caller appears"*, and
> none has (grep-verified at C3). So there are **four** methods, not five.
>
> **The method is `other_sample_records`, not `other_sample_reads`.** The code's name is the older
> one, every call site uses it, and it is the more accurate: this layer counts *records* it stepped
> over, before anything became a read.
>
> **`read_next` does not return `io::Result`.** Two unrelated faults leave it — the reader failing,
> and a read group failing to resolve — and an `io::Result` cannot tell them apart, which is how
> an unresolvable `@RG` came to render as *"reading the next alignment record failed"*. It returns
> `Result<bool, RegionReadError>` and the cursor maps the two arms to `ReadFilterError::Source` and
> `::ReadGroup` (§4).

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
    /// `reset_read_group_counts` — spec §7.
    ///
    /// **Named `read_group_tally`, not `counts`**: the cursor already has a `counts` field of
    /// type `CursorCounts` — what the cursor *did* — and a `read_group_counts()` method whose
    /// value differs from this field's, because the method stamps the `other_sample` rider onto
    /// the first entry. A field and a method spelled alike and returning different things, both
    /// reachable in one function, is what C2's review filed as a Major.
    read_group_tally: Vec<ReadGroupCounts>,
    /// What the layer below had already skipped as another sample's when the current tally
    /// window opened. **Not in this sketch until Checkpoint C**, and it is what makes
    /// `reset_read_group_counts` honest: that count lives on the narrowing and is cumulative for
    /// the life of the cursor, so a reset that cleared only the tally would open a window with
    /// every foreign record the cursor had ever stepped over already in it.
    other_sample_at_window_start: u64,
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
convert_buffered_read     → Err:  failed = true, return the error
verdict_on_aligned_read   → Drop: charge it, continue
                          → Err:  failed = true, return the error
                          → Keep: tally, push onto `kept`, emit
```

> ⚠ **The conversion is a named method, not an inline `self.buffer.decode()`** — added at D2
> (2026-08-04), and it is the only production change Milestone D makes. `convert_buffered_read`
> increments `CursorCounts::reads_converted`, and that count is the **only** observable of the
> ordering this whole document is about: hoisting the conversion above the first filter keeps the
> same reads, charges the same drops to the same reasons and leaves the four acceptance dumps
> byte-identical, so before the counter existed the entire suite passed under it.
>
> **The increment has to be in the callee.** Beside the call it is left behind by a hoist that
> moves the call, and then reports the number the test expects rather than a wrong one — the
> dropped reads `continue` before reaching it — so the test keeps passing. §5 has the decision.

**Contract, unchanged from today's composite.** Lazy; one raw read resident. A fatal error is
yielded once and then the cursor refuses every later region (`CursorError::AfterFailure`,
already present). The order guard still runs on emit.

> ⚠ **One thing did change, and it is worth stating because nothing enforces it.** The old
> `ReadFilter` was a `FusedIterator` with an explicit guard: it stopped on a **clean end of
> input** as well as on a failure. The cursor's guard is `failed` alone, so re-asking a drained
> region re-enters `RegionRawAlignedReads::read_next` instead of short-circuiting.
>
> **Behaviour-preserving today, but only by the grace of two layers below.** All three reader arms
> latch their own "done" state, and the narrowing re-holds the record its early stop consumed, so
> a redundant call moves no counter and loses no read — verified by driving ten of them. That is a
> guarantee the cursor used to make for itself and now *assumes* of layers that never promised it.
> If it is worth keeping, it belongs in the `AlignedReadsReader` contract as a written rule; if
> not, the cursor should latch a clean stop again. **Not decided** — raised at Checkpoint C.
>
> A failed **reposition** also stops the cursor, not only a failed read: the reader's position is
> unknown afterwards, so no later region can be served from it, the reuse path included. The
> reposition therefore happens before any of the new region's state is committed, so a failed jump
> leaves the cursor exactly as it was.

**What the cursor gains publicly:**

```rust
impl<R: RawRefSeq> AlignmentCursor<R> {
    /// Step-1's tally, one entry per read group met. Cumulative — spec §7.
    pub fn read_group_counts(&self) -> Vec<ReadGroupCounts>;   // exists
    /// Start a fresh tally window. The caller chooses the window; the cursor does
    /// not reset on its own, and never per region — spec §7.
    pub fn reset_read_group_counts(&mut self);                 // new
}
```

## 4. Errors

**No new error type**, and no existing variant changes meaning. `ReadFilterError` gained a fourth
variant at Checkpoint C (owner, 2026-08-03) — see the correction below.

| variant | raised by |
|---|---|
| `Source` | the reader failing to hand over a record |
| `ReadGroup` | `resolve_read_group`, on a record that cannot be attributed |
| `Decode` | the conversion |
| `Reference` | the second filter's mismatch check |

`verdict_on_raw_read` cannot fail, which is why it returns a bare verdict.

> ⚠ **This table said "three variants name the three pieces" and it was wrong.** Two unrelated
> faults leave `RegionRawAlignedReads::read_next` — the reader failing, and a record's read group
> failing to resolve — and while that method returned an `io::Result` they were indistinguishable,
> so **both rendered as *"reading the next alignment record failed"***. An operator meeting that
> goes looking for a truncated file when what is wrong is the `@RG` header: a different fault, in a
> different file, wanting a different fix. Found by C2's review against a real BAM.
>
> The split is at the source: `read_next` returns `Result<bool, RegionReadError>`, whose two
> variants the cursor maps to `Source` and `ReadGroup`. Pinned by
> `an_unresolvable_read_group_is_fatal_and_charged_to_its_own_condition`, which fails if the two
> are re-conflated.
>
> **`ReadFilterError::Decode` is unreachable and unpinned**, and that is recorded on the variant
> itself rather than here. No input can reach it: the conversion refuses a record with no reference
> id, no alignment start, or no read group stamped, and the region narrowing guarantees all three
> before it yields. It is kept as defence in depth against that narrowing regressing.

## 5. Design decisions — decided

- **The two filters and the conversion sit below the cursor's kept set; the cursor owns the
  loop.** The rule that stops the pieces being rearranged — spec §3.
- **The conversion is a named private method, `convert_buffered_read`, not an inline
  `self.buffer.decode()`** — and this is the bullet that makes the rule above *enforceable*
  (D2, 2026-08-04). Rearranging the three pieces changes no output, so before D2 nothing in the
  project could tell. The method increments `CursorCounts::reads_converted`; two tests assert the
  identity **`reads_converted` = kept + second-filter drops**, one failing if the conversion is
  hoisted above the first filter and one if a second-filter check is hoisted below it. Inlining
  the method back removes the only pin on the ordering — both tests fail, which is the intended
  alarm, though the failure then reads as a filter having moved.
- **The counter ships; it is not `#[cfg(test)]`.** It lives on `CursorCounts`, the type that
  already answers "what did this cursor do", and is folded across a sample's files by that struct's
  exhaustive-destructure `AddAssign`. The first build used a `#[cfg(test)]` thread-local and it was
  the weaker shape: a reset protocol, thread-global rather than per-cursor, and the crate's first
  `#[cfg(test)]` static. Cost measured at 0.7 % on the walk probe, inside run-to-run noise.
- **A compile-time witness was considered and rejected.** A zero-sized token minted by the first
  filter's `Keep` arm and required by the conversion makes the hoist a build error rather than a
  test failure; it was built and it works. It is not adopted because §1 states this design adds no
  type, and because it pins only the first of the two directions — no type can forbid a second
  copy of the length rule being written against noodles' types and checked early. Reopening it is
  a §1 amendment, not an implementation choice.
- **Two filters, not one** — spec §2.
- **A filter takes a borrow and returns a verdict, never a read** — spec §4.
- **The verdict carries the drop reason**, so an `Option` return is ruled out: the tally is keyed
  on reason and read group.
- **The read group is read off the raw aligned read** — it is already stamped there for exactly
  this ([`filtering.rs:350-357`](../../../../src/ng/read/filtering.rs#L350)).
- **The second filter's #7 → #9 → #8 order is kept**, not the numbering order — spec §4.
- **`FilterState` is deleted for a `failed` flag** — spec §5.
- **The source trait is deleted and the in-memory reader gains a scripted error** — spec §6.
- **The tally lives on the cursor, cumulative, with `reset_read_group_counts`; not on `AlignmentFile`** —
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
| the six flag/MAPQ filters | `verdict_pre_decode` [`filtering.rs:215`](../../../../src/ng/read/filtering.rs#L215) | **rename** to `verdict_on_raw_read`; body unchanged |
| the three conversion-dependent filters | `verdict_post_decode` [`filtering.rs:274`](../../../../src/ng/read/filtering.rs#L274) | **rename** to `verdict_on_aligned_read`; body unchanged |
| the conversion | `RawAlignedRead::decode` [`aligned_read.rs:71`](../../../../src/ng/read/aligned_read.rs#L71) → `decode_record` [`aligned_read.rs:118`](../../../../src/ng/read/aligned_read.rs#L118) | **reuse as-is** |
| the loop | `ReadFilter::next` [`filtering.rs:776`](../../../../src/ng/read/filtering.rs#L776) | **move** into `AlignmentCursor::next_read`; `ReadFilter` deleted |
| the up-front contig check | `ReadFilter::new`'s per-contig fetch loop | ✅ **done at B1** — replaced by **two** checks in `AlignmentFile::cursor`: `self.contigs.first_disagreement(reference.contigs())`, the same comparison the open gate makes [`open_bam.rs:206`](../../../../src/ng/read/input/open_bam.rs#L206), **plus one zero-length fetch on this cursor's own contig**. The comparison alone loses a guarantee — `ResidentRefSeq`/`WindowedRefSeq` take their `ContigList` as a constructor argument unrelated to the bytes, so a matching table can front a FASTA that cannot serve the contig. Spec §9 Q2 |
| the tally and its fold | `ReadFilterCounts` [`filtering.rs:127`](../../../../src/ng/read/filtering.rs#L127), `ReadGroupCounts` [`:564`](../../../../src/ng/read/filtering.rs#L564), `tally_for_current_record` [`:727`](../../../../src/ng/read/filtering.rs#L727), `counts` [`:749`](../../../../src/ng/read/filtering.rs#L749) | **move** to the cursor, including the `other_sample` rider on the first entry |
| the errors | `ReadFilterError` [`filtering.rs:477`](../../../../src/ng/read/filtering.rs#L477) | **reuse as-is** |
| the raw read | ✅ **done at A1** — `RawAlignedRead` [`aligned_read.rs:56`](../../../../src/ng/read/aligned_read.rs#L56), `NoodlesRawAlignedRead` [`:206`](../../../../src/ng/read/aligned_read.rs#L206) | **renamed and moved** to `aligned_read.rs` |
| the region narrowing | `RegionRecords`, [now `region_raw_aligned_reads.rs`](../../../../src/ng/read/input/region_raw_aligned_reads.rs) | **renamed at A3**; trait impl → inherent methods at C3 |
| the per-format readers | `RecordReader` and arms, [now `aligned_reads_reader/mod.rs`](../../../../src/ng/read/input/aligned_reads_reader/mod.rs) | **renamed at A2**; `InMemory` arm gains a scripted error at C1 |
| the three-way stop | `FilterState` [`filtering.rs:549`](../../../../src/ng/read/filtering.rs#L549), `restart_after_end_of_input` [`:659`](../../../../src/ng/read/filtering.rs#L659), `has_failed` [`:675`](../../../../src/ng/read/filtering.rs#L675), `source_mut` [`:697`](../../../../src/ng/read/filtering.rs#L697) | **delete**, all four |
| the source trait and its doubles | `RecordSource` [`filtering.rs:338`](../../../../src/ng/read/filtering.rs#L338), `FakeSource`, `ErroringSource` | **delete**. Note: `RecordSource::header` **already went at B1** — the contig probe was its only caller — so `RegionRawAlignedReads` has no `header()` today. §3.3 lists one; re-add it at C3 only if a caller appears |
| the probe-free constructor and lent buffers | `with_validated_contigs` [`filtering.rs:627`](../../../../src/ng/read/filtering.rs#L627), `ReadFilterBuffers` | **delete** — no caller once the cursor owns the loop. Both are down to **one** caller each since B1 deleted `ReadFilter::new` |
| the cursor's existing filter field and its call sites | [`input/cursor.rs:245`](../../../../src/ng/read/input/cursor.rs#L245), `:355`, `:410`, `:468`, `:471` | **replace** with the fields in §3.4 |

*Line citations refreshed at Checkpoint B (2026-08-03), after Milestones A and B moved them by
up to 120 lines. **Milestone C executes against this table**, so check them before trusting one.*

## 7. Open items

**No open design questions.** Both of the spec's are settled — §9 there keeps the reasoning.

**Impl-time confirmations. All five are answered.** The two that were open for Milestone C:

- ✅ **`verdict_on_raw_read` keeps `(flag, mapq)`**, not `&impl RawAlignedRead`. The owner deferred
  the change at Checkpoint C, and D1 then gave the split form a reason of its own: the reference-free
  capability (spec §5) is pinned by a **function-pointer coercion**
  (`read/reference_free_first_filter.rs`), and a whole-read form would put a trait bound back into
  that signature — the coercion would still work, but it would no longer read as *"a flag, a mapping
  quality and the thresholds, and nothing else"*. A reviewer verified the whole-read form compiles
  and passes; it is available, and there is no longer a reason to take it.
- ✅ **The in-memory reader's scripted error is positional** — `with_failure_at_read(n)`, plus
  `with_failing_seek_at` for the second way a reader can break (C1). An always-failing arm would
  have made the two-fault distinction untestable.

Answered earlier:

- ✅ **`region_records.rs` was renamed on disk** (A3, `git mv` → `region_raw_aligned_reads.rs`),
  as the spec assumed.
- ✅ **The `+ ContigTable` bound landed in one commit of its own** (B1) and reached
  `SampleReads::cursor` and four `PileupGenerator` sites; the SSR generator already required it.
- ✅ **`RecordIndex` became `PackedReadEntry`, not `RawReadIndex`** — §2 has the reasoning
  (owner's call at Checkpoint A).

**One thing this document did not anticipate and Milestone C must know:** `AlignmentFile::cursor`
is now the *only* place that proves an accessor belongs to its file, and `AlignmentCursor::
over_records` — which C2 rewrites — is infallible and carries that as a **prose precondition**.
Three test call sites build cursors through it directly, bypassing the check;
`the_fixture_accessors_carry_the_same_contig_table_as_the_fixture_files` is the standing guard on
them. C2 must not reintroduce a path that reaches `over_records` from outside `cursor`.

**Carried past this plan, for whoever picks the work up next.**

- **The cursor does not latch a clean end of input, only a failure** — §3.4's ⚠ block. Still
  undecided: either it becomes a written rule in the `AlignedReadsReader` contract, or the cursor
  latches again.
- **`ReadFilterError::Decode` is unreachable and unpinned**, kept as defence in depth, recorded on
  the variant. D2's tests do not change this — they measure that the conversion *happens*, not
  that it fails.
- **Nothing calls `reset_read_group_counts`** — spec §7's stated capability, still awaiting its
  first caller.
- **`ResidentRefSeq::new` and `WindowedRefSeq::new` can build a lying accessor**; B1's third check
  contains the damage at the cursor. Closing it is spec §10's "the file owning its reference".
- **Spec §5's capability stops at `AlignedReadsReader`** — a reference-free *pass* is
  constructible, a reference-free *cursor* is not, because `R: RawRefSeq` is unconditional. Found
  at D1 and now written into spec §5.

## 8. Test & bench shape

Tests stay beside the code. `filtering.rs`'s 45 tests split: the ones about a rule stay, the loop ones
move to `cursor.rs`, the three test-double ones are replaced (spec §8).

**The regression anchors are output identity on real data, not the unit suite** — spec §8 has the
four dumps and the walk probe's figures. Three tests that do not exist yet are named there too,
and each covers something no output comparison can see.

✅ **All three exist** (C1, D1, D2), and two of them took a different shape from the one spec §8
sketched — the corrections are recorded there. What is worth carrying here is the general lesson,
because it cost this milestone two rebuilt steps: **a property that changes no output cannot be
tested by a fixture, only by a measurement or by the type system.** D1's first attempt tried to
pin a signature property with a *scope* and pinned nothing; D2's plan tried to pin a work property
with a *fixture* that could not exist. What worked was a signature coercion, an exhaustively-built
config, and a counter inside the step being counted.
