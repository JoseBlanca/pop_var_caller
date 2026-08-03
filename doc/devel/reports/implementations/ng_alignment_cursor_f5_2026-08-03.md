# ng — the alignment cursor, F5: delete the old read path

*Implementation report, 2026-08-03. Branch `ng-generic-perf`, base `2a083c9`.*

Plan: [`alignment_cursor.md`](../../ng/impl_plan/alignment_cursor.md) step **F5**.
Design authority: [spec](../../ng/spec/alignment_cursor.md), [arch](../../ng/arch/alignment_cursor.md).

---

## What the step was

> **F5. Delete the old path — its own commit.** `SampleReads::reads_in_region`, `RegionReads`,
> `ReaderHandle`, `BorrowedReader`, the pool, `readers_opened`, and `region_query.rs` itself.
> Plus the test at `locus_generation/mod.rs:882` and the doc link at `ref_seq.rs:611`.
> **Verification is mechanical:** `cargo build --all-targets` and `cargo test` green, and
> `grep -rn "reads_in_region\|RegionReads\|readers_opened" src/ examples/ benches/` returning
> nothing.

**The verification is mechanical; the step is not.** `grep` proves the old API is gone. It
cannot say whether the *rules* the deleted tests pinned are still checked anywhere — and 22
tests lived in `region_query.rs`, plus 13 in `open_bam.rs` and 15 in `mod.rs` that reached
their reads through `reads_in_region`. So the work ran in three phases, in order: move what
the new path needs out of the doomed module, triage every test against the rule it names, and
only then delete.

---

## 1. What moved before anything was deleted

`record_reader/cram.rs` imported three things from `region_query.rs` that were never the
per-region query's:

| moved | to | note |
|---|---|---|
| `DecodedContainer` (+ `RecordIndex`, `Span`) | `record_reader/container.rs` (new) | the packed container representation |
| `decode_container_at` | same | the one body both CRAM paths shared; now the only one |
| `owner_of_cram_record` | same | read-group resolution on a borrowed CRAM record |

A new module rather than the body of `cram.rs`, because that keeps the existing encapsulation
exactly: `cram.rs` only ever touched `len`, `read_group`, `fill_record` and
`other_sample_records`, and the packing internals stay private to the type.

`pub(crate) use super::region_records::overlaps` was a re-export the deleted sources needed;
the function already lives in `region_records.rs` and did not move.

**Two fields died in the move, and both were dead only once the old path went.**

- **`Footprint`, and `RecordIndex::footprint` with it.** The per-region CRAM source filtered a
  held container's records against the region it was serving, so it wanted each record's
  reference span pre-computed — a CIGAR walk saved per record per re-filter. The cursor's
  reader **positions and never bounds**: it does no region test at all, and the layer above
  applies `region_records::overlaps` to the rebuilt record. So the field had no reader, and
  keeping it would have been a CIGAR walk per record at decode for nobody.
- **`DecodedContainer::offset`.** It was the reuse key: the old source compared it to decide
  whether the container it held was the one the walk had reached. The new reader tracks
  `last_decoded_offset` itself and drops its container on every `begin_region`, so nothing
  reads the field.

Both are behaviour-neutral removals with a small decode-time and residency saving. They are
recorded here rather than absorbed silently because a reviewer reading `container.rs` against
its ancestor will notice two fields missing.

**A trap worth naming:** `record_reader/mod.rs` carries a module-level `#![allow(dead_code)]`,
which covers every submodule under it — including `container.rs`. Neither field would have
warned. They were found by reading, not by the compiler.

---

## 2. The triage — every deleted test, and where its rule went

### `region_query.rs` — 22 tests

**Genuinely dead: the subject itself is gone.**

| test | why |
|---|---|
| `t5_the_indexed_query_returns_exactly_what_a_linear_scan_returns` | superseded by `a_run_of_regions_through_one_bam_cursor_matches_a_linear_scan`, which is the same oracle over a *run* of regions on a 200 kb contig |
| `the_oracle_is_not_vacuous` | superseded by `the_bam_cursor_oracle_is_not_vacuous` |
| `the_indexed_query_returns_reads_in_coordinate_order` | subsumed: both cursor oracles compare **ordered** name vectors against a linear scan |
| `the_scan_stops_once_it_passes_the_region` | the early stop moved to `RegionRecords` at C1 and is pinned there by `the_walk_stops_at_the_first_record_beginning_past_the_region` and `the_record_the_early_stop_took_is_handed_over_next` |
| `the_reader_survives_the_query_and_comes_back_usable` | the pool |
| `the_guard_stays_ended_even_if_what_it_wraps_resumes` | `OrderVerified`'s fuse; a cursor latches instead, pinned by `a_cursor_whose_file_failed_refuses_later_regions_instead_of_answering_short` |
| `a_filter_error_is_wrapped_and_fuses_the_guard` | same |
| `an_inverted_region_is_an_error`, `a_cram_planner_rejects_an_inverted_region` | the planners. **See the open question below** |
| `the_crai_walk_stops_once_it_passes_the_region` | **a rule the new arm deliberately does not have** — a record reader positions, it never bounds. Its death is the feature |
| `a_stale_container_is_not_reused` | **checked by mutation, not assumed.** Removing `self.container = None; self.served = 0;` from `CramRecordReader::begin_region` fails `a_run_of_regions_through_one_cram_cursor_matches_a_linear_scan`. The rule is covered |

**Moved, because the cursor still obeys them.**

| test | new home | new name |
|---|---|---|
| `t4a_a_position_regression_within_a_contig_is_an_error` | `cursor.rs` | `a_read_going_backwards_within_a_region_is_a_fatal_error` |
| `t4c_reads_sharing_a_start_position_are_not_a_regression` | `cursor.rs` | `reads_sharing_a_start_position_are_not_out_of_order` |
| `t4d_querying_a_later_region_then_an_earlier_one_is_not_a_regression` | `cursor.rs` | `an_earlier_region_after_a_later_one_is_not_out_of_order` |
| `records_without_a_footprint_never_surface_from_a_region_query` | `region_records.rs` | `a_record_with_no_footprint_never_surfaces` |
| `the_crai_walk_skips_containers_that_end_before_the_region` | `record_reader/cram.rs` | `containers_ending_before_the_region_are_skipped_rather_than_walked` |
| `the_crai_is_grouped_so_each_contig_sees_only_its_own_entries` | `open_bam.rs` | unchanged |
| `an_unplaced_entry_before_the_placed_ones_does_not_hide_a_contig` | `open_bam.rs` | unchanged |
| `a_contig_absent_from_the_crai_has_an_empty_walk` | `open_bam.rs` | unchanged |
| `a_region_naming_a_contig_the_file_does_not_have_is_an_error` + `a_cram_region_naming_a_contig_the_header_lacks_is_an_error` | `open_bam.rs`, merged | `a_cursor_for_a_contig_the_file_does_not_have_is_an_error` |

`t4b` (a contig-order regression) has **no cursor equivalent by construction**: a cursor covers
one chromosome and `move_to_region` refuses a foreign one before anything moves, which
`a_region_on_another_chromosome_is_refused_and_the_cursor_survives` already pins.

**`CursorError::OutOfOrderRead` had no test at all** before this step — the guard in
`AlignmentCursor::emit` was written, reviewed and shipped at B/D without one. That is the
single biggest thing the triage found, and it is exactly the failure mode the plan warns about:
a guard that never fires looks like a guard that works.

The three `.crai` grouping tests are stated against `group_crai_by_contig` directly rather than
through a planner. That was always their real subject; the planner was one `Vec` index.

### `open_bam.rs` — the pool tests and the composed chain

The whole "reader pool (C1)" section is gone: `sequential_borrows_open_the_file_once_and_reuse_the_reader`,
`concurrent_borrows_each_get_their_own_reader_and_all_return`,
`borrows_from_many_threads_each_get_a_reader_and_all_return`, `a_cram_file_pools_readers_the_same_way`,
`a_failed_open_does_not_count_as_an_opened_reader`, `a_taken_handle_is_not_also_returned_by_the_borrow`,
and with them `t13_many_region_queries_open_the_file_once`,
`concurrent_region_queries_each_get_a_reader_and_bank_every_tally`,
`a_failed_query_returns_no_reader_because_it_never_took_one`,
`abandoning_a_stream_returns_the_reader_and_banks_its_counts` and
`a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally`. Every one of them is
about lending and reclaiming a reader; there is no lending any more.

Two of their properties are not about the pool and were kept, restated:

- `cursors_on_one_file_read_the_same_thing_from_many_threads` — eight threads, one open file, a
  cursor each. The pool made concurrency a property of the *file*; it is now a property of one
  cursor per worker, and what can still fail is sharing that crept in below.
- `a_cursor_is_send_so_a_worker_can_own_one` — the compile-time half.

Converted to cursors: `t9_a_cursor_runs_the_reference_dependent_filter`,
`records_outside_the_region_are_dropped_without_being_counted` (now also covering the
footprint-less record, because this is where the *tally* is observable),
`t10_a_truncated_file_fails_once_and_then_refuses_later_regions` (which gained an assertion —
the cursor outlives the region, so it must refuse later ones rather than merely fuse),
`t8_a_cram_yields_the_same_ordered_reads_as_the_same_bam`, and
`t2b_the_assembly_check_runs_after_the_reads_have_flowed`.

`t4d_a_later_region_then_an_earlier_one_on_one_handle_is_not_a_regression` is subsumed by the
`cursor.rs` version, which is stronger: the old one needed two separate queries to state the
property, the new one states it on one cursor where forgetting to reset `last_emitted` is
fatal.

`a_multi_container_cram_walks_its_crai_and_stops_early` is **half subsumed and half
deliberately gone**, and an in-place comment says so: reaching later containers is covered by
the CRAM cursor oracle over a three-container fixture; stopping early is the rule the new arm
does not have.

### `merge.rs` — 7 tests, and **the first draft of this report did not triage them at all**

They were disposed of in one clause of §3 — *"`merge.rs` (`MergedRegionReads` —
`sample_cursor.rs`'s `MergedCursors` replaces it)"* — on the assumption that a type-for-type
replacement is a test-for-test one. The coverage review found it was not.

| test | successor in `sample_cursor.rs` |
|---|---|
| `t6_two_files_interleave_in_coordinate_order_with_their_file_tags` | `two_files_are_merged_in_position_order` |
| `t6_reads_at_one_position_break_to_the_lower_file_index` | `reads_at_the_same_position_break_to_the_first_file` |
| `t6_the_merged_order_is_identical_across_runs` | subsumed by the tie-break test above: run-to-run identity *is* the tie-break being deterministic |
| `t7_the_same_file_twice_errors_at_the_first_collision` | `the_same_read_from_two_files_is_a_hard_error` |
| `t7_reads_sharing_a_position_but_not_a_name_both_survive` | `different_reads_at_the_same_position_are_not_a_duplicate` |
| `an_exhausted_file_does_not_end_the_merge` | `a_file_with_no_reads_in_the_region_does_not_stall_the_others`, plus the mid-region case inside `two_files_are_merged_in_position_order` |
| `three_files_merge_in_coordinate_order` | **none — this was the gap.** Written for F5 (below) |

**The gap was a Major.** Every `MergedCursors` test in the tree used exactly two files, and with
two files an argmin that scans only the first two slots is indistinguishable from one that scans
all of them: mutating `MergedCursors::argmin`'s scan to `.take(2)` left **all 1,541 tests
green**. A sample of three or more files would silently lose every read of its third and later
files — which reads as a sample sequenced less deeply, not as a fault. The deleted test's own
doc comment said why it existed: *"Three files, to prove the argmin is not accidentally a
two-way compare."*

`three_files_merge_in_coordinate_order` is restored in `sample_cursor.rs`, with the files given
**out of coordinate order** so a merge that concatenated in file order fails too. Verified: the
`.take(2)` mutation now fails it and nothing else.

### `mod.rs` — the sample layer

All 15 converted to `SampleCursor`, keeping every assertion: `t11_the_single_file_arm_matches_the_merged_arm`,
`counts_are_reported_per_read_group_and_not_summed`, `a_cursor_outlives_the_sample_reads_it_was_made_from`,
`a_merged_cursor_outlives_the_sample_reads_it_was_made_from`, `a_cursor_can_be_stored_in_a_struct_without_a_lifetime`,
`a_second_cursor_does_not_disturb_one_already_held`, `a_sample_cursor_is_send_in_both_arms`,
`a_file_with_several_read_groups_resolves_each_record_by_its_tag`,
`an_untagged_record_in_a_multi_group_file_is_fatal`,
`a_file_shared_between_samples_serves_each_open_only_its_own_reads`,
`a_shared_cram_serves_each_open_only_its_own_reads`, `the_whole_stack_over_two_files_and_two_samples`,
`each_read_carries_the_read_group_it_came_from`, and the per-file error wrapping.

One deleted as subsumed: `a_cram_with_several_read_groups_resolves_each_record_by_its_tag` read
the same fixture through the per-region source that `a_cursor_keeps_every_read_group_of_its_sample_not_just_one`
reads through a cursor, asserting the identical `(qname, read group)` triple. Both arms settle
the read group inside `decode_container_at`, which is now the only place a CRAM record's owner
is decided — one rule under two names. An in-place comment records it.

---

## 3. What was deleted

`SampleReads::reads_in_region`, `SampleReads::counts`, `SampleRegionReads`,
`AlignmentFile::reads_in_region`, `AlignmentFile::counts`, `RegionReads`, `ReaderHandle`,
`BorrowedReader`, `ReaderKind`, the `readers` pool with `borrow_reader` / `return_handle` /
`open_reader` / `lock_pool` / `add_counts` / `readers_opened` / `pooled_readers`, `QueryPlan`,
`merge.rs` (`MergedRegionReads` — `sample_cursor.rs`'s `MergedCursors` replaces it),
`region_query.rs` itself, and `examples/ng_cursor_vs_query.rs`.

Plus three things that became **unconstructible** and would have been error cases nothing could
raise:

- `AlignmentFileError::OutOfOrderRead` — the cursor raises `CursorError::OutOfOrderRead`, which
  names the file. Its message test moved to `cursor.rs`.
- `AlignmentFileError::Filter` — the cursor raises `CursorError::ReadRecord`.
- `IngestError::DuplicateReadAcrossFiles` — the per-region merge raised it;
  `CursorError::DuplicateReadAcrossFiles` is the live one. The variant's own doc said
  *"F is where the two become one"*. Its message test moved to `cursor.rs`.

And `ReadFilter::into_parts`, which handed a pooled caller its reader, buffers and tally back.
It was the one deletion the compiler *forced*: with no caller it is dead code, and
`clippy -D warnings` is red until it goes. Its two tests
(`into_parts_returns_the_buffers_with_their_allocations_and_the_tally`,
`with_validated_contigs_adopts_the_lent_buffers`) drove the lend-and-reclaim protocol only, and
went with it.

`examples/ng_cursor_vs_query.rs` **cannot outlive the old path** — its whole job is running both
paths against each other. It produced the headline numbers (23× on BAM, 23.7× on CRAM); both are
recorded in the Milestone C and E reports and in `PROJECT_STATUS.md`, so deleting the harness
loses no evidence.

---

## 4. One thing added, and it is a replacement rather than new scope

**`SampleCursor::read_group_counts`.** Deleting `SampleReads::counts` would otherwise have
removed a capability — a sample's step-1 tally, per read group — with nothing on the cursor
side answering it, and taken four tests with it. The new method folds the k file cursors'
tallies together *by read group*, in first-seen order, exactly as `AlignmentFile::add_counts`
did. `AlignmentCursor::read_group_counts` already existed; this is the sample-level sibling.

**One semantic difference, and the test now says so.** `SampleReads::counts` was empty before
any read, because the file's tally vector was empty until a stream's `Drop` folded one in.
`ReadFilter::counts` always returns at least the `None`-keyed rider it carries for records
skipped as another sample's, so a fresh cursor reports one all-zero entry rather than nothing.
`counts_are_reported_per_read_group_and_not_summed` asserts "no *keyed* tally yet" instead of
"empty", and explains why.

---

## 5. Mutation checks

Every kept test was checked against a mutation of the thing it names. Each was applied, run, and
reverted, and the replacement was confirmed to have applied before running (the plan's own
warning: `cargo fmt` between edit and script turns a no-op replacement into a false "survived").

| mutation | outcome |
|---|---|
| `AlignmentCursor::emit`'s comparison forced false | `a_read_going_backwards_within_a_region_is_a_fatal_error` and `an_earlier_region_after_a_later_one_is_not_out_of_order` **both fail**; nothing else does |
| `move_to_region`'s `self.last_emitted = None` removed | 14 cursor tests fail, including the new one |
| `region_records::overlaps`'s `_ => false` arm → `_ => true` | **only** `a_record_with_no_footprint_never_surfaces` fails — it is the sole cover for that arm |
| `first_entry_reaching`'s binary search replaced by `0` | `containers_ending_before_the_region_are_skipped_rather_than_walked` and `a_container_beginning_earlier_but_reaching_the_region_is_stepped_back_to` **both fail** |
| `CramRecordReader::begin_region` keeps its container | `a_run_of_regions_through_one_cram_cursor_matches_a_linear_scan` fails — which is what retires `a_stale_container_is_not_reused` |
| `MergedCursors::argmin`'s scan → `.take(2)` | **before review: all 1,541 green.** Now fails `three_files_merge_in_coordinate_order`, alone |
| `BamRecordReader::chunks_from`'s `region.is_empty()` guard forced false | **before review: all 1,541 green.** Now fails `an_inverted_region_is_refused_rather_than_answered`, alone |

`first_entry_reaching` became a free function over `&[cram::crai::Record]` so the rule could be
stated against a hand-built index; a `CramRecordReader` needs a real file behind it, and the
property is about the index alone. `record_reader/cram.rs` had **no tests of its own** before
this step.

---

## 6. Verification

| check | result |
|---|---|
| `cargo test --lib ng::` | **1,543 passed**, 0 failed (1,573 at base) |
| `cargo test --lib` (whole crate) | 2,842 passed, 0 failed |
| `cargo test --bins --tests` | green |
| `cargo test --examples` | green |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo fmt` | applied |
| `grep -rn "reads_in_region\|RegionReads\|readers_opened" src/ examples/ benches/` | **empty** |
| `ng_generic_walk_probe` chr21 | `loci=236081 observations=251786 reads_admitted=54709` — exact |
| `ng_generic_loci_dump` chr21 (BAM) | 251,792 lines, **byte-identical** to the `2a083c9` build |
| `ng_ssr_loci_dump` chr21 (BAM) | 4,406 lines, **byte-identical** |
| `ng_generic_loci_dump` SL4.0ch01 (CRAM) | 1,718,914 lines, **byte-identical** |
| `ng_ssr_loci_dump` SL4.0ch01 (CRAM) | 11,945 lines, **byte-identical** |

The line counts match the figures recorded against the `ee0c94b` baseline, so the identity chain
runs all the way back.

The suite is 30 tests smaller: ~35 deleted (the pool, the planners, the order-guard wrapper, the
two lend/return tests) against ~12 added or moved in. Every deletion is accounted for above.

### What the review changed

Two agents, one worktree each, on the staged diff. **No Blockers.** One Major — the `merge.rs`
triage this report did not do, above — and six Minors, all applied:

- `SampleCursor::read_group_counts`'s comment claimed *declaration* order; it is first-**seen**
  order, and a file whose first read carries its second `@RG` reports that one first.
- A broken intra-doc link in `filtering.rs` to the `into_parts` this step deletes — the one
  newly-broken link in a crate that already has 12 pre-existing ones. `cargo doc` is back to 12.
- Four doc sites still describing `region_query.rs` in the present or future tense
  (`region_records.rs` ×2, `test_fixtures.rs`, `record_reader/cram.rs`), including *"when
  Milestone F deletes `region_query.rs`…"* — which is this patch.

The correctness reviewer independently confirmed the two things most worth confirming: the
`container.rs` move is byte-identical to its ancestor apart from the three intended deletions
(diffed with comments stripped, and re-checked with `record_reader/mod.rs`'s module-level
`#![allow(dead_code)]` removed so nothing could hide under it), and `SampleCursor::read_group_counts`
folds correctly — probed with both files carrying foreign records, and with a file that yields
the sample nothing so its rider is `None`-keyed. Neither loses nor double-counts.

**The gate is not the skill's.** `cargo test --release` is red on a clean tree (four tests assert
on `debug_assert!` messages that release compiles out); `cargo test --all-targets` aborts on a
pre-existing panic in `benches/psp_writer_perf.rs:386`; `cargo doc` is red with 12 pre-existing
unresolved links. All three predate this branch.

---

## 7. Deviations from the plan, absorbed and recorded

1. **`region_query.rs` was not deletable on its own** — three CRAM items had to move out first
   (§1). The plan lists only the deletion.
2. **`ReadFilter::into_parts` had to go too**, and it is in `filtering.rs`, which F5's inventory
   does not name. It is dead code the moment the pool goes, so `-D warnings` forces it.
3. **`SampleCursor::read_group_counts` was added** (§4) so the deletion does not remove a
   capability.
4. **`first_entry_reaching` became a free function** so its rule is testable.
5. **`ScriptedRegionReads` in `genome_walk.rs` was renamed** to `ScriptedRegionSource` — a
   test-local struct whose name matched the verification grep.
6. **Three unconstructible error variants were deleted** (§3), one of which its own doc had
   already scheduled for this step.
7. **`locus_generation/mod.rs:882` needed no work.** The plan's pointer is a stale line number:
   what sits there is `sample_reads_over_fixture`, which opens a `SampleReads` and never asks it
   for reads. The verification grep confirms nothing in that file referenced the old API.
8. **`ref_seq.rs`'s doc link is at :645, not :611** — same drift; fixed.

---

## 8. Open, for the owner

**The cursor accepts an inverted region and the planners did not.** `BamRegionSource::plan` and
`CramRegionSource::plan` rejected `region.is_empty()` (start > end) with
`AlignmentFileError::Region`; `AlignmentCursor::move_to_region` validates the **chromosome only**.
An inverted region is not answered dangerously — the overlap test needs a read spanning from
`end` back to `start`, so the answer is empty in practice — but it is answered silently rather
than refused. This is a consequence of `move_to_region`'s shape, settled and reviewed at
Milestones B–D, not a decision F5 made; flagged rather than changed, because adding a check to
`move_to_region` is a design edit.

**The inverted-region guard that survived is now pinned, and the inconsistency is not.**
`BamRecordReader::chunks_from` still refuses `region.is_empty()` — its own comment records that
the first version of the reader dropped the check and a region `80..=70` came back with a read
spanning it. That guard was covered only by the deleted planners' tests, so after F5 it could be
disabled with a green suite; `an_inverted_region_is_refused_rather_than_answered` pins it.
What remains is the disagreement described above: the BAM reader refuses on a *jump*, the CRAM
reader never checks, and a forward region reaches neither because `continue_into` does not
reposition. Recorded, not resolved.

**`ReadFilterBuffers` is now a vestigial seam.** It exists so a caller can lend a filter its two
reused buffers and take them back; the only caller was the pool. Nothing lends any more —
`ReadFilter::new` passes `Default::default()` and `with_validated_contigs` is reached only
through it. Both docs now say so. Folding the parameter away is a `filtering.rs` cleanup, not
this step's.

**`cursor.rs`'s module doc is stale from Milestone A** — it opens with "So far only its errors …
The cursor itself lands in Milestone B" and lists the BAM arm as landing at C and CRAM at E. All
of that is done. Pre-existing, not created here, and not touched because rewriting a module
header is not a deletion step's business.

**And the carried-over one, raised at three checkpoints now.** `arch/locus_generation_pileup.md`
reproduces a deleted three-parameter signature under an "As built" banner and says "one walker
per region" and "one read query per segment"; `spec/locus_generation_pileup.md`, both
`locus_generation_ssr.md` docs and one line of `arch/alignment_cursor.md` are in the same state.
After F5 these are not merely stale — they describe an API that no longer exists anywhere in the
tree. The plan-driven skill does not edit design docs, so this needs the owner.
