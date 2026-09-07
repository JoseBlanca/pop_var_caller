# A3 — the repository goes

**Date:** 2026-09-07
**Branch:** `ng-cram-window`
**Plan:** [cram_reference_window.md](../../ng/impl_plan/cram_reference_window.md), Milestone A step 3
**Design:** [alignment_cursor.md](../../ng/spec/alignment_cursor.md) §10 point 6;
[alignment_file.md](../../ng/spec/alignment_file.md) §5

## Verdict

**Nothing in the read-input layer can hold a chromosome of the reference any more.** The shared
`fasta::Repository`, the one-contig bound that narrowed it, the `unbounded` escape hatch from that
bound, `bases_for_contig`, and the `OpenReference` handle each `AlignmentFile` kept so its cursors
could ask — all deleted. `OpenReference` is now a description plus one question asked at open.
6,280 library tests pass, the same count as before: three tests of the cache went and three of the
check arrived.

The two open-time refusals the plan named as non-negotiable both still fire, and **one of them was
untested until this commit**.

## What landed

`OpenReference` is `Arc<ReferenceInfo>` and one method,
`check_bases_can_be_read`: the reference names a FASTA, and that FASTA's `.fai` parses
(`WindowedRefSeq::read_index`). No sequence is read and nothing is cached — a GRCh38-shaped `.fai`
parses in about 137 µs (spec §10 point 1, measured there), so a thousand-sample cohort pays it once
per file at startup.

`AlignmentFile::open` asks that instead of building a repository, and stops keeping the reference
afterwards: the `reference: Option<OpenReference>` field is gone, because the only reader of it was
the cursor's `bases_for_contig` call.

`AlignmentFile::cursor`'s CRAM arm is now four lines shorter and has one source of bases — the
reader it mints and checks for the decode (A1), which the decode reads a slice's declared span
from (A2). `CramAlignedReadsReader` loses its `repository` field and `decode_container_at` its
argument; the empty repository noodles still requires by value is built inside
`decode_container_at`, named `no_repository`, so what it is for is readable at the call.

## Two deviations from the plan, both small

- **The plan did not mention the `AlignmentFile::reference` field.** It lists what `OpenReference`
  loses and what `CramAlignedReadsReader` trades, and this field is neither — but with
  `bases_for_contig` gone nothing reads it, and its own doc comment described a shared cache of
  bases that no longer exists. Deleted.
- **The `Build` error's source type changed** from `AlignmentInputError` to `std::io::Error`,
  because the operation it wraps changed from `build_fasta_repository` to
  `WindowedRefSeq::read_index`. Its meaning is unchanged; the review then renamed it
  `IndexUnreadable`, below.

## What the check no longer covers, and where that is caught instead

`build_fasta_repository` opened the indexed FASTA — the `.fai` **and** the FASTA. `read_index`
opens only the `.fai`. So a reference naming a FASTA that has been deleted, with its `.fai` left
beside it, now passes `open` and fails one level later: at `AlignmentFile::cursor`, where the
zero-length probe fetches through every reference reader it mints and cannot open the file. That
is spec §10 point 6's shape — *does this reference carry a FASTA, and does its `.fai` open* — and
the probe is the stronger of the two checks, because it proves the cursor's own contig can be
served rather than that the geometry parses.

## The tests

Three went with the repository: `every_caller_gets_one_repository`,
`moving_to_a_new_contig_drops_the_previous_contig_s_bases`, and
`an_unbounded_reference_keeps_every_contig`. Each asserted a property of a cache that no longer
exists.

Three arrived, and the count is unchanged at 305 in `ng::read`:

- `a_reference_with_a_fasta_and_its_index_can_be_read` — the passing case.
- `a_fasta_whose_index_is_missing_is_a_fault_naming_the_fasta` — the `.fai` half of the check, at
  the unit level, asserting that the error names the FASTA.
- `a_cram_whose_reference_has_no_fai_is_refused_at_open` — **the same fault through
  `AlignmentFile::open`, and nothing in this tree held it before.** The missing-`.fai` refusal was
  a side effect of building a repository, so no test had to name it; with the repository gone the
  ask is explicit and this is what fails if a later edit drops it.

`a_cram_against_a_fai_only_reference_is_refused_at_open` is unchanged and green, as the plan
required.

### Mutation-verified

Both new refusals were checked by breaking them on purpose, one at a time, each reverted and the
revert confirmed by `grep` before the next run.

| mutation | `cargo test --lib ng::read` | which tests failed |
|---|---|---|
| `check_bases_can_be_read` returns `Ok(())` whenever a FASTA is named — the `.fai` never read | 303 passed, 2 failed | `a_fasta_whose_index_is_missing_is_a_fault_naming_the_fasta`, `a_cram_whose_reference_has_no_fai_is_refused_at_open` |
| the open-time check never runs (`if false &&` on its guard) | 303 passed, 2 failed | `a_cram_against_a_fai_only_reference_is_refused_at_open`, `a_cram_whose_reference_has_no_fai_is_refused_at_open` |

And the review's added test, the same way:

| mutation | `cargo test --lib ng::read` | which test failed |
|---|---|---|
| the FASTA is left in place, so nothing about the reference is broken | 305 passed, 1 failed | `a_cursor_over_a_fasta_that_has_been_deleted_is_refused`, on its own assertion |

No mutation reached a commit: each file was restored from a copy taken before the edit and the
restore confirmed by `grep` for the mutation's own text before the next run.

## Validation

Run in the container on this tree.

- `cargo test --lib` — **6,280 pass** at the first commit, **6,281** after the review's added
  test; 0 failed, 15 ignored, in 67.67 s.
- `cargo test --lib ng::read` — **306 pass**, 0 failed, re-run as the last thing before `git add`
  so the tree tested is the tree committed.
- `cargo clippy --lib --tests --all-features -- -D warnings` — the three
  `explicit lifetimes could be elided` errors in `run/cohort_merge/{build.rs:820, build.rs:893,
  serial.rs:67}`, all `main`'s and none in a file this step touches.
- `cargo fmt --all -- --check` — nine files drift, all `main`'s: two under `examples/`,
  `ng/psp/block.rs`, five under `ng/run/cohort_merge/`, and `ng/run/psp_source.rs`.
  `open_bam.rs` was formatted here and is clean.

## One thing found while doing this, and what the owner ruled

**A CRAM slice whose records span several contigs would now panic, where before this branch it
decoded.** It is not in A3's diff — A2 introduced it — but the deletion here makes it permanent,
so it was raised at the checkpoint. The owner's ruling: support such slices rather than refuse
them, by decoding one twice. That is Milestone A′ in the plan, and it runs before Milestone B.

Spec §10 point 2 says a slice with no single reference "needs no external bases; it takes the
existing no-window call, whose repository argument is then an empty `fasta::Repository::default()`
and is never consulted". **The second half is false for the multi-contig case.** In
`vendor/noodles-cram/src/io/reader/container/slice.rs:291`, a slice whose reference context
`is_many()` resolves *each mapped record's own* contig through the repository —
`get_record_reference_sequence`, which ends in `.expect("invalid reference sequence name")`.
Against an empty repository that expectation fails, so the process aborts. Until 2026-09-07 the
run's shared repository answered it and the file read correctly, at the cost of holding whole
contigs. Unmapped slices are unaffected: their records are unmapped, so noodles asks for no bases.

### The first count of how rare this is was measured wrongly

**Corrected here rather than quietly.** The first survey searched each `.crai` for the reference
id `-2`, the value a multi-contig slice's *header* carries, and reported none. That check could not
have found one however many there were: htslib writes such a slice to the index as **one line per
contig, all sharing the container offset and the slice landmark**
(`htslib/cram/cram_index.c:715`, `cram_index_build_multiref`), and `-2` never reaches the index at
all.

Re-measured by grouping index lines on (offset, landmark) and counting distinct contigs per group,
the conclusion survives — **0 multi-contig slices in 179 CRAMs and 134,860 slices under
`benchmarks/`**, the whole-genome tomato file's 112,140 included — but the evidence for it is
different, and the wrong method is what let the next sentence be wrong too.

**These slices are not rare in general, which the first report implied they were.** `samtools`
1.16.1 writes one **by default**, with no options set, for a coordinate-sorted file with two short
contigs and five reads each: it merges under-full slices across contigs
(`htslib/cram/cram_encode.c:3964`). So a reference with many short contigs — a draft assembly, a
scaffold-level reference — produces them routinely, and coordinate-sorted input is no protection.
That is what turned the owner's ruling from "prepare for an edge case" into a correctness gap.

The benchmark files are untracked, in the main checkout at
`/Users/jose/devel/pop_var_caller/benchmarks/`, so the count cannot be reproduced from a clone.

## What is not proven here

**Nothing has been measured.** The saving this plan exists for — resident memory on a real
whole-genome CRAM — is B1's, on `benchmarks/tomato_big_cram/DRR000741.p1.cram`, together with the
wall time and a VCF compared line for line. Until then "the decode holds no chromosome" is a
statement about the code, proved by the deletions above and by A2's fixtures, not a measurement.

**A CRAM cohort's descriptor count is still arithmetic.** `DESCRIPTORS_A_CRAM_NEEDS` is 3 and
A3 does not move it: the repository it deleted was one descriptor for the whole *run*, not one per
file, and `DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES` carries 23 of slack. Counting a
real CRAM cohort with `examples/ng_open_cohort_descriptors.rs` is what would replace the
arithmetic, and it remains unrun.
