# A1 — the CRAM decode is handed a reference reader of its own

**Date:** 2026-09-07
**Branch:** `ng-cram-window`
**Plan:** [cram_reference_window.md](../../ng/impl_plan/cram_reference_window.md), Milestone A step 1
**Design:** [alignment_cursor.md](../../ng/spec/alignment_cursor.md) §10 points 1, 5 and 6;
[arch/alignment_cursor.md](../../ng/arch/alignment_cursor.md) §1.2, §1.3, §4

## Verdict

**Plumbing only, and it decodes exactly what it decoded before.** The CRAM decode now owns a
windowed reference reader and does not read it: `decode_container_at` takes it as `_reference`
and step A2 is where it starts fetching bases. 6,275 library tests pass and every integration
binary passes except the one that was already failing on `main`.

## What landed

| file | change |
|---|---|
| `read/input/open_bam.rs` | `AlignmentFile::cursor` takes `impl FnMut() -> R` in place of `reference: R`, and bounds `R: RawRefSeq + ContigTable + Send + 'static`. It mints the cursor's reader as before; on the CRAM arm it mints a second for the decode. **Both go through one `checked_reference` closure** that applies the contig-table comparison and the servability probe — the two checks that are about an accessor rather than about the argument — so minting a reader and checking it are one operation. |
| `read/input/mod.rs` | `SampleReads::cursor` forwards its factory to each file instead of calling it — how many readers a file needs is the format's question, and only `AlignmentFile::cursor` knows the format. |
| `aligned_reads_reader/cram.rs` | `CramAlignedReadsReader` gains `decode_reference: Box<dyn RawRefSeq + Send>`, owned for the reader's life and passed to every `decode_container_at`. |
| `aligned_reads_reader/container.rs` | `decode_container_at` gains the parameter, named `_reference` because A2 is what reads it. |
| `locus_generation/pileup/generator.rs`, `locus_generation/ssr.rs` | `Send + 'static` propagated onto `PileupGenerator` and `SsrGenerator`, which are generic over the accessor and call `SampleReads::cursor`. |

**Why the decode gets its own reader rather than sharing the cursor's** is the ruling recorded
in spec §10 point 1: a reference reader is an open file position plus a resident window, and
this project gives one to every consumer of reference bases rather than sharing — `SampleReads::cursor`
takes a factory for exactly that reason, and one worker already mints several
(`walker.rs:1612-1641`). What is shared instead is the parsed `.fai` and the contig table, by
`Arc`, which is what makes an extra reader cost about 18 µs rather than 189.

## Deviations from the plan, and one reversal

**The plan's first draft of this step was built, discarded, and rewritten.** It threaded the
cursor's own accessor down `read_next` as `&dyn RawRefSeq` — through `RegionRawAlignedReads` and
`AlignedReadsReader` to the CRAM arm. That compiled and its tests passed, and it was wrong:
[`read_filtering_stages.md`](../../ng/spec/read_filtering_stages.md) §5 carries a dated
correction (2026-08-04) stating that those two types carry no reference bound, so that a pass
filtering reads on flag and mapping quality with no reference in scope stays constructible —
and `src/ng/read/reference_free_first_filter.rs` exists to pin exactly that. The threaded
version forced that module to construct a reference accessor, which its own docs say must not
happen. Raised with the owner and reversed; spec §10 point 5 now records both shapes and why
this one won.

**The `Send + 'static` bound reaches further than spec §10 point 6 first said.** It sits on
`AlignmentFile::cursor` and `SampleReads::cursor`, and because both locus generators are generic
over the accessor and call the latter, it propagates to `PileupGenerator` and `SsrGenerator` — 8
bound sites in two files. Every accessor actually handed to a cursor satisfies it already:
`WindowedRefSeq` at the production sites and `InMemoryRefSeq` in the tests. The `Rc`-based
`SharedReference` in `pileup/parity.rs` is `!Send` and is untouched, because it is never handed
to a cursor — it serves the production walk's `MultiChromRefFetcher`.

**One defect was found in this step's own first draft, before it was committed, and it is the
reason `checked_reference` exists.** That draft checked the *contig table* of the cursor's reader
only, and gave the decode's reader nothing but the zero-length probe. The probe cannot catch the
case that matters: a **permuted** contig table carries every name and every length the file
declares, so a fetch at position 1 succeeds — those bases are readable, they are simply another
chromosome's — while `ContigId` is an index, so every read would be decoded against the wrong
chromosome. On the filter's reader that costs a wrong mismatch fraction; on the decode's it
corrupts the read sequences themselves, silently, because a CRAM stores a read as its differences
from the reference. Both readers are now minted and checked in one place.

`a_cram_cursor_checks_the_second_reference_reader_it_mints` pins it with a factory that is right
the first time and permuted the second. **Mutation-verified rather than argued:** replacing the
CRAM arm's `checked_reference(make_reference(), …)` with a bare `make_reference()` leaves
`cargo test --lib ng::read` at **300 passed, 1 failed** — only this test — and the mutation was
reverted and the revert confirmed by grep before anything else ran.

**Three test fixtures became factories rather than values.** `wrong_names`, `wrong_lengths` and
`short` in `open_bam.rs` build deliberately-bad accessors for the three refusal tests; each is
now a closure that builds one on demand rather than a binding, which is a smaller change than
deriving `Clone` on `InMemoryRefSeq` for the sake of three tests.

## Validation

Run in the container on this tree.

- `cargo fmt --all -- --check` — clean.
- `cargo test --lib --tests` — **6,276 library tests pass**, 0 failed, 15 ignored, in 57.12 s.
  Every integration binary passes except `ng_calling_loop_calls_genotypes`, where
  `a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` fails at line 1235.
- `cargo clippy --lib --tests --all-features -- -D warnings` — three errors remain, all
  `explicit lifetimes could be elided` in `run/cohort_merge/build.rs:820`, `:893` and
  `run/cohort_merge/serial.rs:67`.

**Both of those failures are `main`'s, not this step's, and both were checked rather than
assumed**: the test fails identically at the pre-merge commit, and the same three clippy errors
come back from `cargo clippy --lib --all-features -- -D warnings` run in the `main` worktree at
`52b7b787`. Nothing in the changed files is flagged.

## Open, for the review to rule on

**Should the two accessor errors say *which* reader failed?** `CursorAccessorContigTable` and
`Reference` are both keyed on the alignment file's path, which is the same for both readers, so
an operator cannot tell the read filter's reader from the decode's. The file already carries the
precedent for splitting — `CursorAccessorContigTable`'s own doc argues it for a different pair.
Not done here, on the reasoning that both readers now come from one factory and are checked
within microseconds of each other, so the realistic fault fails the *first* one; a second-only
failure means the factory is not deterministic, which is a different thing to name. A `role`
field on both variants would cost 11 call sites across two files.

## What this step deliberately does not do

Read the reader. `decode_container_at` still decodes with `records_discarding_tags` against the
whole-chromosome repository, so no decoded byte moves — which is the point: A2 changes what the
bases are, and an error in that is a silently wrong read rather than a crash, so it lands as a
diff of its own with its own oracles.
