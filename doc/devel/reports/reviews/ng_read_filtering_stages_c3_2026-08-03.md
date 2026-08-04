# Code review — ng read filtering in stages, C3 (the source trait goes)

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `0b5f5e7` (C2)
**Impl report:** [`ng_read_filtering_stages_c3_2026-08-03.md`](../implementations/ng_read_filtering_stages_c3_2026-08-03.md).

---

## 1. Scope

The uncommitted diff for C3 — 6 files, ~38 insertions / ~225 deletions, of which four are
prose-only.

**Categories dispatched:** one `general-purpose` agent in its own git worktree, running
`refactor_safety` (C3 is a deletion, so the question is what the deletion took) plus
`module_structure` (the milestone's checkpoint condition is a module-boundary claim). One rather
than three, in proportion to a step that removes code and moves no logic.

## 2. Verdict

**Request changes** → all applied. **Two Blockers**, two Majors, several Minors.

The deletion itself was correct and the bodies provably unchanged. Every finding was about
something the deletion *exposed*: two lines that had never been tested, a claim in the impl report
that was false, and prose the step made stale or broke.

## 3. Execution status

The agent reproduced the gate: 2,856 passed, fmt clean, clippy clean, `cargo doc` at 12.

## 4. Findings

### Blocker

#### B1: `region_raw_aligned_reads.rs` — the held record's read-group replay has no test
**Confidence:** High (mutation-verified, both directions)

The early stop fires on a record the reader has already handed over, so it is put back — *with*
the read group it arrived with, because on CRAM the group is decided while the container is
decoded and travels attached to the record. **That line had zero coverage.**

- `buf.read_group = None` → **survives all 2,856**. On a multi-read-group CRAM, where the record
  carries no `RG` tag to fall back on, this kills an otherwise valid run at every region boundary.
- `buf.read_group = Some(ReadGroupId(0))` → **survives all 2,856**. This one is worse: a read
  **silently attributed to the wrong library**, which changes no output anyone looks at and
  destroys exactly the per-read-group signal spec §7 exists to protect.

C3's safety argument is "the bodies are byte-identical, so the tests still cover them". Correct in
principle — and for this line the tests never did.

#### B2: `region_raw_aligned_reads.rs` — the early stop's `on_this_contig &&` guard has no test
**Confidence:** High (mutation-verified)

A sorted file groups each contig's records together, so a record of a *later* contig always begins
"past" this region's end by the position comparison alone — while saying nothing about whether this
contig has more records to come. Stopping on it loses every remaining read of the region, silently.

Dropping the guard **passes all 2,856**. The existing fixture cannot see it: its foreign-contig
record sits *inside* the region, so the position test never fires on it. Demonstrated both
directions — `["inside"]` against `["inside", "also-inside"]`.

### Major

#### M1: the impl report's `decode_fails` claim is false, twice
**Confidence:** High

C3's report said the doubles took with them "the only construction of `ReadFilterError::Decode` in
the tree". Neither half holds. `FakeRecord` never constructed a `ReadFilterError` at all — it
returned an `io::Error` that the deleted filter wrapped — and by C3 `decode_fails: true` appeared
**nowhere**: both constructors set `false`, because the two tests that set it were deleted at
**C2**. The field was already dead code.

`ReadFilterError::Decode` did become unreachable *and* unpinned — at C2, not C3. Its doc also named
a cause it cannot have, "the unmapped flag clear yet no position", which is one of the shapes the
layer below discards before the conversion is reached.

#### M2: prose this patch introduced or made stale
**Confidence:** High

- `input/mod.rs` now read "the region narrowing **from `filtering.rs`**" — it is `input/`'s own.
- `read/mod.rs` read "which applies the filter consumes" — a dangling verb from the same sweep.
- `read/mod.rs` claimed `ReadFilterCounts` is `pub(crate)`; it is `pub` and must be, as part of
  `AlignmentCursor::read_group_counts`' return type.
- `filtering.rs` claimed the module "needs `RecordBuf` alone"; after C3 it names no noodles type.
- `cursor.rs` still said "C3 deletes the doubles… this has to exist first", in the past.

### Minor / cross-category

- **`AlignmentCursor`'s doc comment is severed** — pre-existing from **C2**, in an in-scope file.
  The insertion of the error and tally types landed *inside* the doc, so a 27-line block describing
  the cursor was attached to `ReadFilterError` and the cursor's own rustdoc was the orphan fragment
  *"is no stream object to give back."*

## 5. What the review confirmed

- **The bodies are identical.** A comment-stripped diff of the old `impl RecordSource for
  RegionRawAlignedReads` against the new inherent block is exactly three signature lines.
- **No `header()` caller appeared**, verified by grep — so omitting arch §3.3's fifth method was
  right.
- **On the naming clash the code is right:** `other_sample_records` matches the reader layer below
  and counts *records*, not reads. Arch §3.3's `other_sample_reads` should move.
- **The checkpoint condition is met.** `filtering.rs` holds exactly spec §6's six items.

## 6. Out of scope observations

None beyond the severed doc, which was fixed here because it is in an in-scope file and was this
milestone's own doing.
