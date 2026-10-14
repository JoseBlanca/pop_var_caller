# Code Review: str_bundles_reads_psp
**Date:** 2026-10-11
**Reviewer:** rust-code-review skill (orchestrator)
**Scope:** step 3 of the STR bundle caller (reference runs, the bundle generator, the psp bundle record, format 2.0), uncommitted diff on `str-bundles` at 3db4e2fc
**Status:** Request-changes

---

## 1. Scope

- **What was reviewed:** the uncommitted diff (25 modified files) and three new files on branch
  `str-bundles`, HEAD `3db4e2fc`. It adds `ReferenceRun` and `find_reference_runs`
  (`src/region_typing/reference_runs.rs`). It adds the bundle generator
  (`src/locus_generation/ssr_bundle.rs`), which is the STR segment generator run over a hull with
  the flat-gap delimiter. It replaces the segment inside `SsrLocus` with `MeasuredRepeat`
  (`src/locus_generation/ssr.rs`) and fills the bundle slot (`src/run/walker.rs`). It stores
  `LocusKind::SsrBundle(SsrBundleDetail)` as a trailing byte string on every psp record, with
  format 2.0 (`src/psp/{record,header}.rs`). It also adds the scan route's `FlankClamped` count for
  contig-end bundles, the bundle cell of `pvc inspect-psp`, `examples/ng_str_bundle_runs.rs`, and
  spec edits.
- **In-scope files:** the files above; the `LocusKind` match arms in `src/calling/`,
  `src/run/{callers,explain,records}.rs`, `src/run/cohort_merge/close.rs`,
  `src/parameter_estimation/`, `src/window_coverage/depth.rs`; and
  `doc/devel/ng/spec/{str_bundle_loci,locus_generation_ssr,psp_file_format,run_ssr_observations}.md`.
  Direct callees read: `alignment/ssr_best_path_flat_gap.rs` (traceback),
  `alignment/ssr_unit_robust.rs` (`resolve_anchors`, `AnchorRule`),
  `psp/segmentation_section.rs` (the bundle-cap key).
- **Out of scope:** calling bundles (step 4), the simulator (reviewed separately), and the typing
  (reviewed in `str_bundles_typing_2026-10-11.md`).
- **Categories dispatched.** By the owner's rule, no worktrees. Two read-only reviewers worked in
  this worktree and neither built nor edited anything. The first covered *reliability + refactor
  safety*. The second covered *errors + defaults + naming + idiomatic + smells + module structure +
  extras*, including doc quality. I ran every build and measurement myself. `float_portability`
  does not apply: no float in the diff reaches output. `unsafe_concurrency` does not apply either:
  there is no new concurrency, and the generator's `Send` bound is checked by the compiler at
  `GeneratorSlot`. Their files are in `tmp/review_2026-10-11_str-bundles-reads-psp/`.

## 2. Verdict

Request-changes. Every STR segment is measured exactly as before. One reviewer checked the
per-edge escape check and the `MeasuredRepeat` conversion against HEAD case by case. The
non-bundle records were byte-identical with main's build on the stuttering benchmark. But two
defects decided before the review give wrong bundle observations:

- the runs miss motifs no member carries (M1);
- about 1 crossing read in 100 at a bundle is reported complete although it ends inside the hull
  (M2).

The psp reader also accepts a bundle run lying past the hull (M3), and the 2.0 compatibility
story is stated wrongly in three places (M4).

## 3. Execution status

All in the dev container, release build, before any fix
(`tmp/review_2026-10-11_str-bundles-reads-psp/baseline.sh`):

- `cargo fmt --check` → exit 0.
- `cargo clippy --release --all-targets -- -D warnings` → exit 101, on two lints in code this diff
  does not touch. One is `clippy::type_complexity` at `src/locus_generation/ssr.rs:2422`, which is
  the same line at HEAD. The other is `clippy::single_match` at
  `src/run/edge_deletion_gathering.rs:201`. Both are pre-existing and out of scope.
- `cargo test --release --no-fail-fast --lib -- locus_generation region_typing psp:: run:: calling:: window_coverage parameter_estimation cli:: repeat_catalog`
  → 3,787 passed, 0 failed, 11 ignored.
- `cargo test --release --test ng_repeat_catalog_differential` → 13 passed.
- `cargo test --release --lib -- cross_platform_digests` → 4 passed.
- Not run: `cargo doc`, `cargo audit` (no dependency changed).

**The implementer's quantitative claims.**

- HG002 runs per bundle (0:141, 1:922, 2:1,210, 3:124, 4:12, 5:1; 3,767 runs):
  **CHECKED-CORRECT** against `tmp/step3/hg002_runs.err`, the output of the committed example.
- The spelling table and the byte-identity of non-bundle records: **CHECKED-CORRECT** against the
  outputs in `tmp/step3/`. Both are re-measured after the fixes, in the fix report.
- *"`(CT)20`, `(AT)8`, `T8` is what the rule gives when `AT` is among the motifs; both outcomes are
  pinned by tests"*: **CHECKED-CORRECT**
  (`the_ct_at_t_bundle_comes_out_as_three_runs_…`, `a_motif_no_member_carries_is_not_looked_for`).
- *"35,116,701 to 35,309,906 bytes (0.55%)"*: the arithmetic is correct (193,205 bytes, 0.550%).
  The files were not re-measured.

Needs verification: 0. M7's premise (the test fixture reaches no dropped read) is settled by
rewriting the test so that it reaches dropped reads by construction.

## 4. Open questions and assumptions

1. **Should 1.x psps be refused up front (`OLDEST_READABLE_MAJOR = 2`)?** Today they pass the
   version check and are then refused for lacking `max-str-bundle-len-bp`, unless they came from
   this branch's own builds since 270da29f. The review fixes the documentation (M4) and leaves the
   behaviour for the owner. Affects M4.
2. **Should reads that end inside a bundle still count toward window coverage?** See Mi10. Today
   they count nowhere, so a bundle hull reports about half its read depth.

## 5. Top 3 priorities

1. **M2** — reads that end inside the hull are reported complete: about 1 crossing read in 100
   becomes a false short allele. This is wrong evidence in the psp.
2. **M1** — runs of motifs no member carries are never looked for: the spec's own example loses its
   16-base `AT` run to fixed sequence, with no stutter priced in it.
3. **M3** — the psp reader accepts runs past the hull: step 4 will index the hull by run offset, so
   a damaged file would panic instead of being refused.

## 6. Findings

### Blocker

None.

### Major

#### M1: src/region_typing/reference_runs.rs:100-120 — Runs are looked for only in the members' motifs
- **Confidence:** High
- **Problem:** Step 1 of `find_reference_runs` iterates over `members` and scans each member's
  motif. At GRCh38 chr1:8,007,443 the typing lists the members `CT` (all 64 bases) and `T`, and no
  `AT` member. So the hull `(CT)20 (AT)8 T8` comes out as `(CT)20` and `T9`, and the 15 bases
  between them become fixed sequence. That is a 16-base `AT` repeat in which no stutter will be
  priced. The test `a_motif_no_member_carries_is_not_looked_for` pins this outcome.
- **Why it matters:** The members overlap and nest (§5.1: 551 of HG002's 2,410 bundles hold one
  member inside another), so the members' motifs are an incomplete list of where slips happen.
- **Suggested fix (decided by the owner before the review):** look for maximal pure stretches of
  every motif of period 1 to 6. A motif and its rotations are one motif, and a stretch is counted
  under its smallest period, so `ATATAT` is period 2. Keep the copy floor and the longest-first
  rule. Exclude any unit holding a base other than `A`, `C`, `G` or `T`, so that no run of `N` is
  ever found. Update §6.1, §4 and §11.

#### M2: src/alignment/ssr_best_path_flat_gap.rs:766-769 — The flat-gap delimiter calls a flank anchored from the alignment's shape
- **Confidence:** High (measured by the implementer: `str_bundle_loci.md` §10 step 3)
- **Problem:** `left_anchored: left_flank_len > 0 && tract_start != 0` asks only whether the path
  crossed the junction. Take a read that ends inside the hull. Spending its last bases as
  mismatches against the right flank is cheaper than a 15-base deletion, so the path crosses,
  `classify` says `Between`, and `classify_read` returns a complete observation spelled over part
  of the hull. The segments' default delimiter has the evidence test that stops this
  (`ssr_unit_robust.rs`, `resolve_anchors` with `AnchorRule`). The bundle path never reaches its
  own `EndedInsideABundle` drop for these reads.
- **Measured:** on the stuttering benchmark (PCR, 30×), 7.4 crossing reads in 1,000 at bundles
  are spelled over part of the hull only, at 1 sample, and 7.7 at 50 samples. At segments the
  figures are 0.0 and 0.2.
- **Why it matters:** each such read is a complete observation of a shorter allele that does not
  exist. The candidate bar (two reads, a tenth of the sample's) usually keeps them out, but they
  enter the read likelihoods.
- **Suggested fix (decided):** give the flat-gap delimiter used for bundles the same anchor test,
  as an opt-in rule (`SsrFlatGapAligner::with_anchor_rule`). `SsrFlatGapAligner::new` stays
  unchanged, so the parity oracle, the bake-off examples and every other user keep algorithm 3's
  positional test. Segments use the unit-robust delimiter and are untouched. Reuse
  `resolve_anchors` and `FlankEvidence` from `ssr_unit_robust.rs` rather than copying them.

#### M3: src/psp/record.rs:1129-1180 — The psp reader accepts bundle runs that lie past the hull
- **Confidence:** High
- **Categories:** reliability, errors (convergent)
- **Problem:** `read_bundle_payload` refuses runs that overlap, that descend, or that hold no whole
  copy. It never compares a run's end with the record's `reference_bases`, which
  `decode_record_body` has already decoded (`record.rs:1246`). The round-trip fixture
  `every_locus_kind_round_trips` writes runs `(0,40)`, `(40,17)` and `(300,8)` on a 7-base hull,
  and the record decodes as valid.
- **Why it matters:** step 4 will index the hull by run offset (§6.2). A damaged file would then
  panic on a slice index instead of being refused by name.
- **Suggested fix:** pass `reference_bases.len()` into `read_bundle_payload` and refuse a run whose
  end lies past it. Make the fixtures' runs fit their own hulls, and add a refusal case.

#### M4: src/psp/header.rs:486-497, psp_file_format.md, str_bundle_loci.md §6.2 — The 2.0 compatibility story is stated wrongly
- **Confidence:** High
- **Categories:** errors, extras
- **Problem:** two separate errors.
  - **"A 2.0 reader still reads 1.x" is true only for files no released build wrote.**
    `segmentation_section.rs:394-399` refuses a psp without `max-str-bundle-len-bp`, and every psp
    written by main lacks it. The key arrived in 270da29f without a version bump. So the only 1.x
    files this reader opens are 1.2 files from this branch's own builds since that commit.
  - **The reason given for 2.0 is wrong.** `FORMAT_VERSION`'s doc says a 1.x reader would "take
    its payload for the next field". As built, the payload is the **last** body field. A 1.x
    reader would skip it as an unknown final field, which this diff's own test
    `a_manifest_that_renames_reorders_or_drops_a_field_is_refused` pins. The real hazards are the two that `psp_file_format.md` names: the reader
    would read every bundle as a tag with no runs, and it would ignore the bundle cap.
- **Why it matters:** someone holding an archive of main's psps reads "still reads 1.x" and expects
  not to regenerate. Every one of them has to be written again.
- **Suggested fix:** state what is true in all three places. A 1.x header passes the version check,
  but a 1.x psp is read only if it records the bundle cap, which no released build wrote, so in
  practice every 1.x psp must be written again. Give the real reason for the major bump. Whether
  to refuse 1.x up front is open question 1.

#### M5: src/cli/command_line.rs:49-50, src/cli/call_from_alignments.rs:35-37, src/run/callers.rs:1245-1246 — `--help` and docs still say bundles are "not built yet"
- **Confidence:** High
- **Problem:** the bundle slot is filled, so a bundle's ground is handled and
  `unhandled_not_implemented` stays at 0 for bundles. The diff's own test asserts this. Yet:
  - `--help` says a cluster of repeats is "charged to `not built yet` and reported at the end";
  - `call_from_alignments.rs` says it "has no generator and no model";
  - `bundles_set_aside`'s doc says "the bundle generator is deferred".
- **Why it matters:** `--help` points a user at a counter that now reads zero for this ground.
- **Suggested fix:** say that an STR bundle is measured and stored, not yet called, and counted as
  set aside.

#### M6: src/region_typing/reference_runs.rs:337 — No test tells the cut-back arithmetic from its two obvious wrong versions
- **Confidence:** High
- **Problem:** lines 135-139 apply two rules. A piece cut at its end goes back to whole copies; a
  piece cut only at its start keeps its natural part copy.
  `a_cut_end_goes_back_to_whole_copies_and_a_cut_start_keeps_the_part_copy` exercises neither:
  its start-cut piece is period 1, which has no part copy, and its end-cut piece `(AC)7` is
  already whole copies. A mutant that never cuts back passes every test, and so does one that
  always cuts back.
- **Why it matters:** a run's length decides which bases carry stutter and which are fixed.
- **Suggested fix:** add the reviewer's two hand-traced fixtures, `T (CAG)8 G30 T` and
  `T G30 (CAG)8 CA T`, which give 21 and 26 bases where the wrong rules give 23 and 24. Name the
  existing test for what it shows.

#### M7: src/locus_generation/ssr_bundle.rs:290 — The generator-level accounting test cannot show that reads ending inside are dropped and counted
- **Confidence:** Medium
- **Categories:** reliability, errors
- **Problem:** `a_bundle_keeps_no_partial_and_accounts_every_read` uses an all-`A` reference and
  all-`A` reads. It never asserts `ended_inside_a_bundle > 0` or that any observation exists. Its
  identity check (`supported + without == fetched`) holds by construction of `tally`. A bundle
  built with `KeptAsPartial`, or a tally arm that adds to `outside_tract` instead, could leave it
  green.
- **Suggested fix:** build a reference with flanks that continue nothing around a hull of two
  runs. Use reads copied from it, some crossing the hull and some ending in its middle. Assert the
  exact `ended_inside_a_bundle`, at least one complete observation, and no partial.

### Minor

#### Mi1: src/psp/record.rs:224-231 — `BODY_FIELDS`' doc describes the abandoned draft layout
- **Confidence:** High
- **Problem:** the doc names four `bundle-run-*` fields that do not exist. It says a bundle's flanks
  are in `repeat-*-flank`, but they are inside the payload. It counts twenty-one fields where there
  are twenty-four.
- **Suggested fix:** describe the two trailing fields every body carries.

#### Mi2: src/window_coverage/depth.rs:197-202, :620; src/locus_generation/mod.rs:1007-1008 — Comments say bundle records are never emitted
- **Confidence:** High
- **Suggested fix:** say that a bundle record covers its hull as a segment's does, and drop "Only
  one slot is filled today".

#### Mi3: src/psp/record.rs:3001-3012 — The refused-payload test accepts `Truncated` for every case
- **Confidence:** High
- **Problem:** each case has one right answer, `Malformed` on `str-bundle-payload` with its own
  reason. The `||` lets a regression to an unrelated truncation pass.
- **Suggested fix:** assert `Malformed` and a reason substring per case.

#### Mi4: src/psp/record.rs:1353 — No test puts a payload on a repeat-tract record
- **Confidence:** High
- **Suggested fix:** run the existing refusal test over `Generic` and `Ssr`.

#### Mi5: src/locus_generation/ssr.rs:2329 — The per-edge escape check is never tested where only the right edge escapes
- **Confidence:** High
- **Suggested fix:** add a right-flank-continues case with edges `(None, Some(3))` → `true`.

#### Mi6: src/locus_generation/ssr.rs:2317 — The flank-only test never reaches the ordering it claims to pin
- **Confidence:** High
- **Problem:** the read is settled by the CIGAR pre-check before any delimiter call. A mutant moving
  the `Dropped` check above `partial`'s outside-tract test passes the test.
- **Suggested fix:** use a read whose soft clip faces the hull (`8M6S`), which reaches the
  delimiter, and assert `OutsideTract`.

#### Mi7: src/locus_generation/ssr_bundle.rs:122-142 — `unknown_edges` claims "fixed sequence at both edges"
- **Confidence:** High
- **Problem:** `None` in `EdgeRunPeriods` means fixed sequence, not unknown. The fetch reads only the
  span, so this is harmless today, but the name contradicts the type.
- **Suggested fix:** rename it and say in a comment that the fetch reads only the span; or
  restructure in step 4.

#### Mi8: src/cli/inspect_psp.rs — The `motif` cell changed meaning for bundle rows and `COLUMNS_VERSION` stayed at 2
- **Confidence:** Medium
- **Problem:** the module's version policy, read literally, covers only adding, removing, renaming
  or moving a column.
- **Suggested fix:** an owner's call. Either bump to 3, or add a `reference_runs` column.

#### Mi9: src/psp/mod.rs:268-277 — `UnsupportedVersion` reads as "too new" when it refuses a too-old file
- **Confidence:** High
- **Suggested fix:** carry the oldest readable major in the error and print the range.

#### Mi10: src/window_coverage/depth.rs:203 — A bundle hull reports about half its read depth to the coverage windows
- **Confidence:** High
- **Problem:** the reads dropped for ending inside the hull are counted at no position. A segment's
  partials count at the positions they witnessed. On the stuttering benchmark (1 sample, 30×), a
  bundle cell has 22 crossing reads and 21 reads with no observation, against 20 and 6 at a
  segment. So the hull's positions report about half the depth reads give them.
- **Why it matters:** a coverage window (500 bases) that overlaps a 35-base bundle reads about 3%
  low there. That cuts slightly against the hidden-duplication filter's signal, which is excess
  depth.
- **Suggested fix:** an owner's call (open question 2). One option is to count a dropped read at
  the hull positions it witnessed, in the coverage depth only. Recorded, not fixed here.

#### Mi11: Other naming and structure points, recorded
- `MeasuredRepeat::tract_len` is a hull's length for a bundle; `repeat_len` would say what it is.
- The field `SsrBundleGenerator::hulls` holds a generator.
- `EdgeRunPeriods` and `bundle_hull` are `pub` with no caller outside `locus_generation`.
- `ReferenceRun` is interchange data (locus, psp, inspect, calling) defined in `region_typing`.
- The bundle generator's settings error is reported as `RunError::TractGeneratorSettings`.
- `SsrBundleDetail: Default` is documented as tests-only but is the decoder's placeholder.
- An empty member list silently becomes the delimiter motif `N`.
- `read_bundle_payload`'s inner errors carry offsets into the payload, not into the body.
- The 0-based conversion in the `FlankClamped` charge (`hull.start - 1`) is not pinned at a BED
  edge.

### Nits

The members' order cannot change the last tie-break (`motif.as_bytes()`), which never decides
anything. `LocusGenerator::<crate::…::SsrSegment>` could import `SsrSegment`. `resolve_features`
returns an unnamed tuple. The example's `walk.by_ref()` is not needed.

### The 67 `PARALOG_LR` changes — judged correct counting, not double counting

Region typing partitions the reference. A bundle's hull is one typed region, so no generic or
segment record covers its positions, and `depth.rs` reports a bundle record at every position of
its span, as for a segment. Before this change the hull had no record, so its positions were absent
from the coverage windows: not zero, absent (`depth.rs`: a position with no record is not a covered
position). Now they are present, once, at the depth of the reads that crossed the hull. Those reads
are real and were previously unseen. The only bias is the under-count of Mi10. So the 67 records
whose `PARALOG_LR` moved, out of 2,396, moved because their 500-base windows gained correctly
counted positions. No position is counted twice.

## 7. Out of scope observations

- `src/locus_generation/ssr.rs:2422` and `src/run/edge_deletion_gathering.rs:201` fail
  `clippy -D warnings` (pre-existing, newer lints).
- None of the STR generators' per-reason counts (`outside_tract`, `ended_inside_a_bundle`, …)
  reaches a run summary. This predates the diff.
- A cohort mixing 1.x and 2.0 psps would show the 1.x samples as uncovered at every bundle once
  step 4 calls bundles (open question 1).
- `generate_psps.rs:455-462` still describes `unhandled_not_implemented_bp` as bundle ground. No
  production slot is `NotImplemented` any more.

## 8. Missing tests to add now

- `find_reference_runs_cuts_an_end_cut_piece_back_to_whole_copies` and
  `find_reference_runs_keeps_the_part_copy_of_a_piece_cut_only_at_its_start` (M6).
- `find_reference_runs_finds_runs_of_motifs_no_member_carries`, the chr1:8,007,443 hull with
  members `CT` and `T` → `(CT)20 (AT)8 T8` (M1).
- `find_reference_runs_counts_a_stretch_under_its_smallest_period` and
  `find_reference_runs_finds_no_run_of_n` (M1).
- `the_flat_gap_delimiter_with_the_anchor_rule_demotes_a_read_ending_inside` and
  `the_flat_gap_delimiter_without_the_rule_is_unchanged` (M2).
- `decode_record_body_refuses_a_bundle_run_past_the_hull` (M3).
- `ssr_bundle_generator_drops_and_counts_reads_ending_inside_the_hull` (M7).
- `reference_run_escapes_the_tract_returns_true_when_only_the_right_run_continues` (Mi5).
- `classify_read_counts_a_flank_read_with_a_clip_facing_the_hull_as_outside` (Mi6).
- `decode_record_body_refuses_a_payload_on_a_repeat_tract_record` (Mi4).

## 9. What's good

- `MeasuredRepeat` names the four facts in which a hull differs from a segment, and segments
  convert into it with their old behaviour (`src/locus_generation/ssr.rs`, `From<&SsrSegment>`).
  Both a reviewer's case analysis and the byte-identical non-bundle records confirm the
  conversion.
- The bundle payload is a field every record carries, so a later reader can always step over it
  (`src/psp/record.rs`, `put_bundle_payload`). That keeps the 1.2 rule for adding fields.
- `ended_inside_a_bundle` is a fifth no-observation reason folded into the destructuring total, so
  a dropped read is counted, not lost (`src/locus_generation/ssr.rs`, `reads_without_observation`).
- The payload's bytes are pinned byte for byte (`a_bundle_records_payload_is_these_exact_bytes`),
  which catches a symmetric reordering that a round trip cannot.

## 10. Commands to re-verify

- `tmp/review_2026-10-11_str-bundles-reads-psp/baseline.sh` in the dev container.
- `cargo test --release --no-fail-fast --lib -- reference_runs ssr_bundle locus_generation::ssr psp::record alignment::ssr_best_path_flat_gap`
  after the fixes.
