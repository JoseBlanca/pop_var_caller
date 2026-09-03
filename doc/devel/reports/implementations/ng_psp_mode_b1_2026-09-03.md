# ng psp mode — B1: `SampleObservationGatherer` over the existing walker

**Date:** 2026-09-03
**Plan step:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone B, step B1
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §5.2; **Arch:** [run_streaming.md](../../ng/arch/run_streaming.md) §3.3
**Branch:** `ng-psp-mode`

## Plan

Build psp mode's walk-stage object: one sample's alignment files in, that sample's
observations out in genome order as an `Iterator`, and a driver that drains them into a
`PspWriter` record by record — serial within the sample, reusing the direct-mode chain
(`build_read_groups` → `SampleReads::open` → `generic_path_generators` →
`SampleLocusObservationsIterator` over `RunSegments`) rather than writing any new walking
code. The gatherer also builds the psp `Header` at construction — the first production
code that builds one — from spec §6.1's parts: the segmentation inputs whole, the
walk-local read-group table, the configured reach ceiling, and the read filters recorded
as provenance.

Two commits preceded this step in Milestone B, both owner-ruled at Checkpoint A:

- `cba10a0b` — a pre-existing standing failure fixed forward: the production
  `psp_writer_perf` bench panicked under `cargo test --all-targets` (documented in
  PROJECT_STATUS since before this branch), because both `psp_writer_phases` sub-benches
  assumed the byte cap was the writer's only block cut. Measured before fixing: the prime
  loop exhausted all 3,300,000 fixture records with 49 projected bytes left. Fixed with
  `new_with_block_layout` and a window no position can cross; an assert now names the
  premise.
- `a1fdab11` — `SegmentationInputs` lifted to `src/ng/segmentation_inputs.rs` (psp and run
  become mutually dependent this milestone; the lift resolves the direction), the `ng::run`
  re-export kept so no call site churned, and the Checkpoint A rulings recorded in the plan
  and spec §6.1.

## Assumptions (choices the plan left open)

- **The census configuration is not in the constructor yet.** Plan B1's text lists "a
  census configuration" among the constructed-from items, but the census wiring is
  Milestone G, and no census configuration type exists (`CensusWriter::new` takes eight
  positional arguments, with an explicit `clippy` note rejecting a config struct). A dead
  parameter pinned by no test would be recorded speculation, so the constructor grows the
  census inputs at G1 instead. Recorded in the type's doc comment.
- **`WriterProvenance` is split caller/gatherer.** The caller supplies what only it can
  know (tool, version, subcommand, command line, timestamp — a fixed timestamp is also what
  makes B3's byte-identity comparison possible); the gatherer overwrites what it knows
  better: the input basenames from the files it actually opens, and the read filters via
  `ReadFilterConfig::provenance_parameters()` through `record_parameters` (the A-review
  note for B1).
- **The manifest is `Manifest::as_this_build_writes_it()`**, no block-size knob. Nothing in
  Milestone B needs one; F2's "default block size" is exactly this value. A knob can be
  added at the command if a later step wants it.
- **The multi-sample refusal reuses `IngestError`'s shape.** The gatherer replicates
  `SampleReads::open_only_sample`'s three-way match (that call throws away the read-group
  table the header needs), and wraps the same `IngestError::SampleNameMismatch` /
  `IngestError::ReadGroups` in a new `RunError::NotOneSamplesFiles` — new because
  `OpeningSample` names one sample and this failure is that no one sample could be
  established.
- **Two `RunError` variants added**: `NotOneSamplesFiles` (above) and `PspNotWritten`
  (`create`/`finish` failures, which have no locus for `RecordNotWritten` to name — the
  path locates them). Both are walk-stage refusals arch §6 did not enumerate; its three
  named additions (`AnalysedRegionsDiffer`, `SegmentationInputsDiffer`,
  `SampleAppearsTwice`) are Milestone E's and remain unbuilt.
- **Door checks kept to what the object owns**: empty file list, generator-settings check,
  bases-holding reference, one-sample check. The cohort-scope refusals direct mode also
  runs (descriptor headroom, assembly check against the checksummed reference, catalog
  digest against the reference) stay with the command that loops samples (plan C1:
  "refusals shared with direct mode").

## Changes made

- `src/ng/run/gatherer.rs` (new): `SampleWalkInputs` (the per-sample counterpart of
  `AlignmentInputs`), `SampleObservationGatherer` with `open` / `header` / `sample_name` /
  `reached` / `counts` / `write_psp`, an `Iterator` impl with the walker's exact
  failure shape (`SourceFailed` naming the sample and its progress), and a manual `Debug`
  that prints no reads.
- `src/ng/run/mod.rs`: `pub mod gatherer;` + re-export; the two new `RunError` variants
  beside `RecordNotWritten`.

## Tests added

All in `gatherer.rs`; every asserted setting differs from its type's default (the
callers.rs rule), so a gatherer that recorded shipped constants fails.

- `open_fixes_the_header_a_calling_run_will_check` — two files, two read groups, two
  libraries: the walk-local table numbered from zero in file order; the reach ceiling is
  the configured 4,321 (asserted ≠ default); segmentation inputs recorded whole; contigs
  and both digests from the opened reference; caller provenance surviving beside the
  recorded filter values; basenames overwritten.
- `the_gatherer_yields_what_the_direct_walk_yields` — the same fixture through the bare
  direct-mode chain gives equal observations (guarded non-empty).
- `write_psp_round_trips_the_walk` — the file reads back record for record as the walk
  streamed; header equal to the gatherer's plus the store's own `zstd-compression-level`;
  `WriteStats.records` equals the walk's count. (The fixture-scale smoke of B2's oracle;
  B2 adds the shared run fixtures and a real CRAM slice.)
- `files_naming_two_samples_are_refused` — `NotOneSamplesFiles`, chain naming both samples.
- `an_empty_file_list_is_refused` — `NoAlignmentFiles`.
- `a_reference_without_bases_is_refused` — `ReferenceHasNoBases`.

## Validation results

In the container (`scripts/dev.sh`), on the exact tree reported:

- `cargo fmt` — applied; `cargo fmt --check` clean.
- `cargo clippy --all-targets --all-features -- -D warnings` — exit 0.
- `cargo test --lib 'ng::run::gatherer'` — 6 passed, 0 failed.
- Full-suite state is the branch's known one: red only on the three documented locus-dump
  behaviour failures that predate this branch (PROJECT_STATUS, "Three tests in the two
  locus dumps still fail"); the bench panic that used to sit beside them is fixed
  (`cba10a0b`).

## Tradeoffs and follow-ups

- The gatherer exposes `counts()` but not `generators()`; C2's per-sample report will say
  whether the per-generator counts are needed and can add the accessor then.
- Read-filter tallies (spec §8) live in cursors the generators own and are unreachable
  from here — the same standing gap the walker documents; nothing new is lost.
- The half-written file a failed `write_psp` leaves is refused by every reader as
  interrupted (format guarantee); the command-level exit, message and `--force` rerun
  behaviour are C3's to pin.
