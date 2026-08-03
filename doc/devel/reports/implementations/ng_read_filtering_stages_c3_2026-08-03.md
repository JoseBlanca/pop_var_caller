# ng — read filtering in stages, C3: the source trait goes

**Date:** 2026-08-03 · **Branch:** `ng-generic-perf` · **Base:** `0b5f5e7` (C2)
**Plan:** [`read_filtering_stages.md`](../../ng/impl_plan/read_filtering_stages.md) step **C3**.
**Design authority:** [spec](../../ng/spec/read_filtering_stages.md) §6, §11 ·
[arch](../../ng/arch/read_filtering_stages.md) §3.3, §6.

---

## 1. Plan

> Delete the `RecordSource` trait and its two doubles; `RegionRawAlignedReads`'s trait
> implementation becomes inherent methods. Nothing generic consumes it once C2 lands.

The trait existed so a filter living *apart* from the cursor could be driven by test doubles. C2
moved the loop into the cursor, leaving `RecordSource` with one production implementation and no
generic consumer at all.

## 2. Changes made

- **`RecordSource` deleted** from `filtering.rs`, with its associated type and its two
  required methods.
- **`FakeSource` and `FakeRecord` deleted**, along with `fake`, `mapped`'s fake-only uses and the
  one test whose subject was the seam.
- **`RegionRawAlignedReads`'s four methods are inherent**: `jump_to`, `continue_into`,
  `other_sample_records`, `read_next`. Bodies unchanged; `read_next` gains the doc the trait used
  to carry, including why `Ok(false)` is *this region is done* rather than *the file is finished*,
  and that the caller knows the difference because the caller is what set the region.
- **Prose swept** across five files that named the trait, plus `region_raw_aligned_reads.rs`'s
  module doc, which described a layering (`ReadFilter` between this and `AlignedRead`) that C2
  removed and a Milestone-B/C ordering argument that no longer applies to anything.

## 3. Assumptions and deviations, recorded

**Arch §3.3 lists a fifth method, `header()`, and it is deliberately absent.** `RecordSource::header`
went at B1 with the contig probe — its only caller — and arch §6 already records the instruction:
*"re-add it at C3 only if a caller appears."* None has, so it is not re-added.

**Arch §3.3 also spells the method `other_sample_reads`; the code has `other_sample_records`.** The
code's name is the one that has always been there and the one every call site uses, and it is the
more accurate of the two — the layer counts *records* it skipped before anything became a read.
Recorded as a doc-side discrepancy for the checkpoint rather than renamed.

## 4. Tests

**One test died: `fake_source_drives_the_seam`.** Its subject *was* the seam — it drove
`FakeSource` through `read_next` and then called the two verdicts by hand, mimicking a loop that no
longer exists in that shape. Every property it touched is covered on the real chain: the verdicts
have their own direct tests in this module, and the read/convert/filter order is what the cursor's
walk tests and the four acceptance dumps exercise.

**A claim this report first made about `FakeRecord::decode_fails` was wrong, and the review caught
it.** The first draft said the doubles took with them "the only construction of
`ReadFilterError::Decode` in the tree". Two things are false in that. `FakeRecord` never
constructed a `ReadFilterError` at all — it returned an `io::Error` that the deleted filter
wrapped — and by C3 `decode_fails: true` appeared **nowhere**: both constructors set `false`,
because the two tests that set it were deleted at **C2**. So the field was already dead code, and
C3 removed dead code rather than the last live construction.

What is true is that `ReadFilterError::Decode` became **unreachable *and* unpinned at C2**, not at
C3. That is now recorded on the variant itself — where it survives, since the note C1 wrote lived
in the test block C3 deleted. Its doc also named a cause it cannot have ("the unmapped flag clear
yet no position"), which is one of the shapes the layer below discards before the conversion is
reached; corrected.

### Two Blockers, both in the `read_next` this step rewrote

The review mutated the inherent method rather than trusting that identical bodies stay covered, and
found **two lines with no test at all**. C3's own safety argument — "the bodies are byte-identical,
so the tests still cover them" — is exactly right and exactly why this was worth checking: for these
two lines the tests never did.

| line | mutation | result before |
|---|---|---|
| the held record's read group replay | `buf.read_group = Some(ReadGroupId(0))` | **survived 2,856** — a read silently attributed to the wrong library |
| the same line | `buf.read_group = None` | **survived 2,856** — kills an otherwise valid multi-read-group CRAM run at every region boundary, since such a record carries no `RG` tag to fall back on |
| the early stop's `on_this_contig &&` guard | dropped | **survived 2,856** — a record of a later contig ends the walk and every remaining read of the region is lost |

Two tests added, each killing its own mutation alone:
`the_held_record_carries_its_own_read_group_into_the_next_region` and
`a_record_of_another_contig_does_not_trigger_the_early_stop`. The existing fixture could not see
the second because its foreign-contig record sits *inside* the region, so the position test never
fires on it.

### A defect this milestone introduced at C2, found here

`AlignmentCursor`'s doc comment was **severed**: C2's insertion of the error and tally types landed
inside it, so a 27-line block describing the cursor was attached to `ReadFilterError` and the
cursor's own rustdoc was the orphan fragment *"is no stream object to give back."* Repaired by
moving the inserted block above the doc rather than through it.

**Suite: 2,857 → 2,858 (+1)** — one deleted (`fake_source_drives_the_seam`), two added by the
review. Fully accounted.

## 5. The milestone's checkpoint condition, met

The plan's Checkpoint C requires *"`filtering.rs` holding only the keep-or-drop rules and their
thresholds"*. After C3 the file's entire top-level surface is:

| item | spec §6 lists it? |
|---|---|
| `ReadFilterConfig` (+ its `Default`) | ✅ |
| `FilterVerdict` | ✅ |
| `DropReason` | ✅ |
| `ReadFilterCounts` (+ `add`, `record_drop`) | ✅ |
| `verdict_on_raw_read` | ✅ the first filter |
| `verdict_on_aligned_read` | ✅ the second filter |

Nothing else. Its imports are production's reused predicates, the read it judges, the reference and
`types` — nothing from `read::input`, so the module cycle C2 closed stays closed.

**1,787 → 781 lines** across the milestone.

## 6. Validation

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo test --lib` | **2,858 passed**, 0 failed, 5 ignored |
| `cargo clippy --all-targets --all-features -- -D warnings` | 0 errors, 0 warnings |
| `cargo doc --no-deps` | **12** unresolved links — the pre-existing baseline |
| four acceptance dumps, `cmp` | **byte-identical** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` |

`cargo doc` reached 13 when `RawAlignedRead`'s import left with the trait, orphaning an intra-doc
link in the module header; the count caught it and it is back to 12.

## 7. Follow-ups

- **`ReadFilterError::Decode` is unreachable and unpinned** — §4, and recorded at the code.
- **Arch §3.3's method list is stale in two ways** (`header`, `other_sample_reads`) — §3, for the
  checkpoint.
