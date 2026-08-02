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
