# ng — the alignment cursor, Milestone B: implementation report

*Plan: [impl_plan/alignment_cursor.md](../../ng/impl_plan/alignment_cursor.md).
Design: [spec](../../ng/spec/alignment_cursor.md), [arch](../../ng/arch/alignment_cursor.md).
Branch `ng-generic-perf`. One section per plan step.*

**Milestone B builds the forget rule** — the one part of this design that can lose reads
*silently*. A rule that drops a read it should have kept produces a wrong genotype, not a
crash, and a first attempt lost **3,830 of 236,081 loci while all 1,471 unit tests passed**.
So it is built against a scripted list with no file behind it, where what a region should
return is answerable by scanning the same list by hand.

---

## B1 — the cursor, without the rule

### Plan

Arch §1.2 and §2.2: `AlignmentCursor` over `RecordReader::InMemory` — `move_to_region`,
`next_read`, `contig()`, and the kept reads.

### Changes made

| file | change |
|---|---|
| [src/ng/read/input/cursor.rs](../../../../src/ng/read/input/cursor.rs) | `AlignmentCursor`, its three public methods and `kept_reads`, and 13 tests driven against a hand-written linear scan of the same script. |
| [src/ng/read/input/region_records.rs](../../../../src/ng/read/input/region_records.rs) | **new.** The layer between the record reader and the filter: contig test, overlap test, sorted early stop, read-group resolution. Six tests. |
| [src/ng/read/filtering.rs](../../../../src/ng/read/filtering.rs) | `done: bool` → a three-state `FilterState`, and `restart_after_end_of_input`. Three tests. |

### The deviation that matters: `RegionRecords` is a milestone early

The plan puts it at **C1**, "lifted out of `BamRegionSource`". It cannot wait, and the plan's
own dependency graph is what makes the stated order impossible: a cursor yields
`AlignedRead`, so it owns a `ReadFilter`, which needs a `RecordSource` underneath it — and C1
depends on B2. It is written here in the shape the in-memory arm needs; **C1's remaining job
shrinks to proving the BAM arm reuses this type rather than growing a second copy**.

To keep that from becoming a second copy immediately, the overlap rule is **borrowed rather
than rewritten**: `region_query::overlaps` — what the BAM and CRAM sources already apply — is
lifted to `pub(crate)` and called from here. The one thing that must never happen in this
design is two paths disagreeing about which reads a region contains, and sharing the function
makes them provably identical rather than identical-looking. Milestone F, which deletes
`region_query.rs`, moves the function here.

### The owner-approved change to `ReadFilter`

A region's end reaches the filter as an *ordinary* end of input — `Ok(false)` is the only end
a `RecordSource` can report — which fused it permanently, so the first region a cursor drained
silenced it for the whole chromosome. `done: bool` is now `Running` / `EndOfInput` / `Failed`,
and repositioning clears **only** the middle one: a filter that met a corrupt record or an I/O
fault is not resurrected to read the next region out of a file already known to be broken.

### What B1 deliberately does not do

**No forget rule.** `move_to_region` throws away every kept read, every time — correct, and no
faster than the code it replaces. The rule lands next, as a diff against machinery that is by
then actually exercised.

### Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ng::` — **1512 passed; 0 failed; 2 ignored** (1492 at the milestone's
  start).

---

## B2 — the forget rule

Arch spec §6. **Reuse what is held only when the new region begins at or after the last one
served; otherwise drop everything and jump.** Eviction is its mirror: drop a kept read once it
ends before the current region begins.

### Two things the rule needed that the plan did not anticipate

- **Reuse needs more than not-clearing.** `RegionRecords` repositioned unconditionally, so
  keeping the reads would have handed some out twice. It now has `jump_to` (repositions) and
  `continue_into` (does not), and the cursor picks by the rule.
- **The early stop was eating one read per region boundary.** The stop fires *on* a record,
  which has already been taken from the reader; a region that continues had nothing that would
  ever produce it again. `RegionRecords` holds it and hands it over first — **in that layer,
  not in the `RecordReader` where arch §1.3 puts it**, because the over-read happens where the
  region's end is known. A reader cannot hold back a record it was never told to stop at.

### Verification

The plan's oracle, over all five shapes it names — ascending, backward, overlapping, adjacent,
far-apart — through one cursor, against a linear scan at every step; plus a property test over
random scripts and random runs of regions.

## B3 — the counters

`reads_decoded`, `reads_replayed`, `regions_reusing` / `regions_jumping`, `reads_evicted`.
**The saving is invisible in the reads** — a cursor that keeps nothing and one that keeps
everything return identical output for identical regions, which is the correctness
requirement — so the only way to see the rule work is to count what it avoided.

## The review of B2 and B3, and what it changed

A differential harness over **3,000 random cases** — scripts with shared starts, deletions,
skips and zero-span CIGARs; runs of up to 13 regions including zero-width, duplicate,
adjacent, contig-edge and past-the-end — found **no case where the rule loses a read**. What it
found instead:

- **Blocker: `read_end`'s `.max(1)` had no test.** Deleting it left all 1,524 tests green. It
  is the line that keeps a replayed read's footprint equal to a freshly-read one, and without
  it an all-soft-clip read is yielded when read fresh and dropped when replayed — the answer
  changing with how the caller happened to walk. Now pinned by a test that drives the *same*
  record through both rules at five regions.
- **A cursor whose file had failed answered later regions short.** After a fatal read the
  filter is finished for good, but the cursor still holds reads — so it kept producing
  *plausible, truncated* answers rather than empty ones. Measured: region 60..=100 answered
  `[]` where a scan had two reads. `ReadFilter::has_failed` and `CursorError::AfterFailure` now
  make it refuse.
- **A latent trap for Milestone C, stated in the contract rather than left to be discovered.**
  `continue_into` depends on `begin_region` *positioning* and never *bounding*: reading on must
  yield every record to the end of the chromosome, not just the ones inside the region it was
  handed. noodles' `query()` returns an interval-bounded iterator, which satisfies every other
  word of the contract and breaks this one — and would lose every record past the previous
  region's end, silently, for every region after the first.
- `continue_into` and the pushback had no direct tests in their own module; they have four now.

### The mutation that survives, and the claim about it that was wrong

Turning the kept walk's early `return None` into `continue` fails no test, and B2's commit
said B3's counters would pin it. They cannot: kept reads never decrease in position, so
walking on only reaches a layer that early-stops for the same reason — same answer, same
counters. But the review refuted the *absolutism* of the follow-up claim that no test could
ever distinguish them: walking on consumes records off the reader that stopping leaves alone,
and once the reader is a real file that read can fail, turning a clean `None` into an error
and a dead filter. So it is a short-circuit that also narrows what can go wrong. The comment
now says that.

### Validation

- `cargo fmt --check`, `cargo clippy --all-targets --all-features -- -D warnings` — clean.
- `cargo test --lib ng::` — **1530 passed; 0 failed; 2 ignored**.
