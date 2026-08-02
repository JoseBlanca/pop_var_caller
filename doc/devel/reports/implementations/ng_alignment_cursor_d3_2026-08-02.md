# ng — the alignment cursor, D3: the STR generator

*Plan: [impl_plan/alignment_cursor.md](../../ng/impl_plan/alignment_cursor.md).
Design: [spec](../../ng/spec/alignment_cursor.md), [arch](../../ng/arch/alignment_cursor.md).
Review: [ng_alignment_cursor_d3_2026-08-02.md](../reviews/ng_alignment_cursor_d3_2026-08-02.md).
Branch `ng-generic-perf`, on top of `5ec6a60`.*

D3 points the microsatellite generator at a cursor, the way the previous step pointed the
generic one. It is a smaller change than that one, because the STR generator does not walk:
it fetches the reads over each repeat once, into a list.

---

## What changed

| file | change |
|---|---|
| [ssr.rs](../../../../src/ng/locus_generation/ssr.rs) | `SsrGenerator` holds a `SampleCursor` for one chromosome and the identity of the sample it was opened for. `fetch_capped_reads` is pointed at a cursor instead of being handed a sample and a factory, and the region widening it needs comes with it. `cursor_counts`, and the sample guard. `MF` drops off the type parameters. Six tests. |
| [read/input/mod.rs](../../../../src/ng/read/input/mod.rs) | `SampleIdentity` and `SampleReads::identity` — see below. |
| [locus_generation/mod.rs](../../../../src/ng/locus_generation/mod.rs) | `LocusGenerationError::ForeignSample`. |
| [pileup/generator.rs](../../../../src/ng/locus_generation/pileup/generator.rs) | the same sample guard, and two tests for it. |
| [cursor.rs](../../../../src/ng/read/input/cursor.rs) | `AddAssign` for `CursorCounts`, replacing four hand-written folds. |
| [ng_ssr_cohort_stutter.rs](../../../../examples/ng_ssr_cohort_stutter.rs) | one generator per sample instead of one shared between them. |
| four research tools | type annotations only, following the dropped type parameter. Their baselines do not move. |

---

## The defect this step is really about

Review found it, and it is worth stating on its own because the code looked right.

**A generator opens a reader for one sample and keeps it for a whole chromosome — but the
`LocusGenerator` trait hands it a `&SampleReads` afresh on every call.** I keyed the kept
reader on the chromosome alone. The sample argument was then read only on the call that
opened the reader, and ignored on every call after it.

`ng_ssr_cohort_stutter` asks every sample about one repeat before moving to the next repeat.
So the first sample would open the reader and answer for all of them. Fifty-one plants would
have produced one plant's reads, fifty-one times, under fifty-one names — with no error, no
empty rows, and a table of exactly the right shape. The question that tool exists to answer is
whether samples differ.

Nothing already measured is affected: the recorded tomato stutter results predate this branch,
and that code opened a query per repeat, so it was never exposed.

### The fix, and why not the obvious one

The obvious fix — open a new reader whenever the sample changes — is worse than making no
change at all. The tool walks region-major, so the sample changes at *every* repeat: it would
open one reader per sample per repeat, against the old code's one query per repeat.

So, on the owner's decision: **one generator per sample, and the generator refuses a sample it
was not opened for.** That is what the architecture already assumes — it counts open files as
`files × generators × workers` — and it turns an unwritten precondition into something the code
checks.

`SampleReads::identity` is the new piece of interface that makes the check possible. It
compares samples by **the files they read**, not by name: a name is the wrong test in both
directions, since one individual sequenced twice can be opened as two samples under one name,
and one file set can be opened twice under names that differ only in spelling.

The same guard is on the generic generator. No caller loops samples through that one today, but
both generators are handed their reads by the same trait, so the mistake is equally available;
it was simply made on the STR one first.

---

## Deviations

**`make_reference` is boxed rather than deleted**, exactly as in the previous step, and for the
same reason: `SampleReads::cursor` needs a factory, because a sample's files each need their own
reader. The type parameter goes; the field stays. It is now called once per file per
chromosome, where it used to run at every repeat.

**It is not `Send`-bounded here, and the generic generator's is.** The previous review asked for
`Send` as free insurance for a future parallel run. On this generator it is not free:
`ng_ssr_cohort_stutter` hands every file a clone of one shared reference reader, and that reader
holds a `RefCell`, so it cannot be shared across threads and neither can a closure returning it.
Adding the bound broke that tool outright. Every caller of the *generic* generator builds a fresh
reader per call, so the bound costs nothing there and is kept.

---

## Verification

**Both dumps byte-identical** to binaries built from `ee0c94b` — `ng_generic_loci_dump`
(251,792 lines) and `ng_ssr_loci_dump` (4,406 lines), on chromosome 21 of HG002 at 30×.

**The real-data anchor unchanged:** `loci=236081 observations=251786 reads_admitted=54709`.

**`cargo test --lib ng::` 1,566 passed.** Clippy with warnings as errors, `cargo fmt --check`
and `cargo test --examples` all clean.

### The number

`ng_ssr_loci_dump` on chromosome 21, two runs each:

| | user CPU | wall |
|---|---:|---:|
| `ee0c94b`, a query per repeat | 17.32 s / 17.08 s | 11.00 s / 10.64 s |
| D3, one cursor per chromosome | 12.48 s / 12.45 s | 11.09 s / 11.08 s |

**−28 % of CPU time, and the output is byte-identical.** Wall time does not move, because this
tool's wall clock is set by the whole-genome reference checksum running on another thread, not
by the walk.

**Read it with its region shape.** Microsatellite repeats are far apart — the walk covers
102,659 of them on chromosome 21 — so consecutive repeats share far fewer reads than the generic
generator's tiling regions do, and the saving is correspondingly smaller. The audit at the
previous checkpoint measured the ratio collapsing with distance between regions (23× when they
tile, 3.8× at thirteen times that spacing, and a two-fold *loss* walking backwards), so a
measurement was owed here rather than an argument. This is it.

**Memory: +1.3 MB** peak (55.0 → 56.3 MB), which is one cursor's kept reads.
