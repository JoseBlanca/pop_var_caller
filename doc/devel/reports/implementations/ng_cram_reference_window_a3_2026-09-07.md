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

## One thing found while doing this, and it is the owner's to rule on

**A CRAM slice whose records span several contigs would now panic, where before this branch it
decoded.** It is not in A3's diff — A2 introduced it — and no file this project has ever met
contains such a slice, but the deletion here makes it permanent, so it is recorded before the
checkpoint rather than after.

Spec §10 point 2 says a slice with no single reference "needs no external bases; it takes the
existing no-window call, whose repository argument is then an empty `fasta::Repository::default()`
and is never consulted". **The second half is false for the multi-contig case.** In
`vendor/noodles-cram/src/io/reader/container/slice.rs:291`, a slice whose reference context
`is_many()` resolves *each mapped record's own* contig through the repository —
`get_record_reference_sequence`, which ends in
`.expect("invalid reference sequence name")`. Against an empty repository that expectation fails,
so the process aborts. Until 2026-09-07 the run's shared repository answered it and the file read
correctly, at the cost of holding whole contigs. Unmapped slices are unaffected: their records are
unmapped, so noodles asks for no bases at all.

**Nothing in reach is affected, measured.** Every `.crai` under `benchmarks/` — 180 CRAMs, and the
whole-genome tomato file among them at 112,140 slices — was read for its reference ids: 1,876
slices with `-1` (unmapped) in that file and **not one `-2` (multi-reference) anywhere**. Those
files are untracked, in the main checkout at `/Users/jose/devel/pop_var_caller/benchmarks/`, so the
count cannot be reproduced from a clone; the command was
`gzip -dc <file>.crai | awk -F'\t' '{c[$1]++} END {...}'` over `find benchmarks -name '*.crai'`.
`samtools`
writes them when a file is name-sorted or unsorted, and for runs of small contigs; a fragmented
assembly is where one would first appear, which is the same case `open`'s comment on index-building
already warns about.

**Recommendation: refuse such a file at the decode, naming it, rather than either panicking or
building a window path for it.** One `if` where `reference_span()` returns `None` and the slice's
records are not all unmapped, returning an `io::Error` that says ng cannot decode a slice spanning
several contigs. It costs a few lines, it converts an abort into a message, and it does not commit
the project to a per-record window design for a file shape nobody has. The alternative — teaching
the fork to resolve a multi-contig slice record by record through the `RawRefSeq` — is real work in
the vendored crate and should wait for a file that needs it. Either way this is a change to spec
§10 point 2, so it is not the implementer's to make.

## What the review changed (commit 2)

The step's review found three things worth fixing and four nits; all seven were applied in a
follow-up commit rather than by amending, so the record shows what was found.

- **The missing-`.fai` message had got worse, and the test could not see it.** Reporting the
  fault as `AlignmentFileError::Open` kept only the inner `io::Error` — and `fai::fs::read` is a
  bare `File::open`, whose failure says `No such file or directory` with no path in it. The top
  line therefore read *"opening alignment file '<reference>.fa' failed"*: the wrong kind of file,
  a FASTA that is present and readable, and no mention of the index. It now has its own variant,
  `AlignmentFileError::CramReferenceIndexUnreadable`, naming the CRAM, the FASTA and the `.fai`
  and saying to run `samtools faidx`. Both tests now assert the message, not only the variant.
- **The check that moved out of `open` had no test where it landed.** A FASTA deleted with its
  `.fai` left behind is now caught by `cursor`'s zero-length probe, and nothing held that:
  the neighbouring `a_cram_cursor_checks_that_its_second_reader_can_serve_the_contig` fails at
  the index lookup, before any file is opened. `a_cursor_over_a_fasta_that_has_been_deleted_is_refused`
  holds it, and leaving the FASTA in place fails it — 305 passed / 1 failed, the one being this
  test on its own assertion.
- **Two comments in `container.rs` said opposite things about multi-contig slices.** The older
  one repeated spec §10's claim that such a slice consults no reference. Corrected, along with
  the same sentence in the fork's own `Slice::reference_span` doc, which is where the wrong
  belief came from.
- **`ReferenceBasesError::Build` was renamed `IndexUnreadable`.** Nothing is built any more. This
  contradicts the plan's "keeps its name and its meaning", and the plan is right that the
  *meaning* is unchanged — but the name described the deleted operation.
- The `open` comment numbered "5." sat sixteen lines above the code it described, with the
  `.crai` grouping comment butted against it; moved down to the check. The tombstone comment
  where the deleted `reference` field was is now one line pointing at that check.

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
