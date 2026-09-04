# psp record head — H2: the owning documents follow the code

**Date:** 2026-09-04
**Plan step:** [psp_head_compared_reads.md](../../ng/impl_plan/psp_head_compared_reads.md) Milestone H, step H2
**Spec:** [psp_head_compared_reads.md](../../ng/spec/psp_head_compared_reads.md) §1, §2
**Branch:** `ng-psp-mode`

## Plan

H1 gave the record head two fields — the keep rule's denominator, and the locus kind's tag brought
forward from the body. This step makes the documents that own the format say so, in the same
milestone, so no reader of a spec meets a head that no longer exists.

Documents only. No code.

## Changes made

**[`psp_file_format.md`](../../ng/spec/psp_file_format.md), which owns the layout.**

- §2's vocabulary entry for *the record head* lists all six fields.
- §4.3's diagram carries them in wire order.
- §4.3's opening sentence said the head answers *the two questions a reader has*. It now names
  what the head answers — the ground covered, the kind of locus, whether enough reads varied to be
  worth building, and how far to skip — and says which two of those it answered when the section
  was settled. **⚠ [`psp_head_compared_reads.md`](../../ng/spec/psp_head_compared_reads.md) §3
  quotes the old phrase**; the paragraph it points at now narrates the change rather than carrying
  the sentence verbatim.
- The field-by-field list gains an entry for each new field, each with the reason it is there:
  the kind, because the width bound governs generic loci only and the never-mix assertion needs it
  before any evidence is assembled, and because it is a **move** — the tract's motif and flanks
  stay in the body, so nothing has two answers to check; and the denominator, with what it is
  counted over, why depth was the wrong choice, and the head-only refusal that comes with it.
- The measured-cost paragraph says its figures predate the sixth field, that the kind's move should
  cost nothing to first order, that the denominator is expected to compress worse than the
  numerator because it tracks depth where the numerator is almost always zero, and that H3
  re-measures. Its "the four head fields" is now "the head's four scalars as it then stood".

**[`psp_record_encoding.md`](../../ng/spec/psp_record_encoding.md) §2.3** — the head line updated.

**[`run_streaming.md`](../../ng/spec/run_streaming.md) §3.3** — the ⚑ is marked met and points at
the spec that met it. **And its three-step sketch is corrected**, which the plan did not ask for
but the ⚑ sits directly under: step 2 asked a sample for *how many reads in total*, which is depth
— the denominator this spec's own §3 rules out, because it raises the bar with reads that could
never clear it. It now asks for the record's reach, the locus kind, the reads compared whole
against the reference, and how many of those varied.

**[`cohort_merge.md`](../../ng/spec/cohort_merge.md)** — §1.3's *position summary* entry gains the
kind and the denominator and says what each is for; §13's deferred bullet on the position summary's
encoding is struck through as landed, naming the record head as the encoding and the two fields
that joined the three it originally asked for.

**[`arch/psp_file_format.md`](../../ng/arch/psp_file_format.md) §3.1** — not on the plan's list,
included because it carries a `RecordHead` sketch that no longer matched the type. Both fields
added, with a line saying why `LocusKindTag` exists.

## Validation results

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0;
  `cargo test --lib` **6,157 passed**, 0 failed, 14 ignored — unchanged from H1, as a
  documents-only step should leave them.
- Every line citation added here was opened and read in this tree before it was written.

## Left standing, deliberately

- **Dated log entries that spell the old head are not rewritten**: `PROJECT_STATUS.md`'s C2 entry
  and [`impl_plan/psp_file_format.md`](../../ng/impl_plan/psp_file_format.md)'s completed C2 step.
  They record what was built on a date, not what the format is now.
- **[`examples/ng_psp_head_encoding.rs`](../../../../examples/ng_psp_head_encoding.rs) still
  measures a four-scalar head** and so under-reports it. That is H3, the next step.
- **Three claims in `src/ng/psp/` that predate this branch**: the record module's doc and
  `BODY_FIELDS`' doc both say a record's chain ids are dropped and arrive at Milestone E, and
  `chain_ids.rs` says nothing calls its writer yet. All three were already false at the branch
  point, since the chain ids landed at E4; none is this milestone's to fix.
