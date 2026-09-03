# ng psp mode — B2+B3: the gatherer's oracle, and byte identity

**Date:** 2026-09-03
**Plan steps:** [run_driver_psp_mode.md](../../ng/impl_plan/run_driver_psp_mode.md) Milestone B, steps B2 and B3 — **bundled deliberately**: B3's oracle is one env flag on B2's harness and one test beside B2's, over the same fixtures; two loop iterations would review the same diff twice.
**Spec:** [run_streaming.md](../../ng/spec/run_streaming.md) §2, §12.1, §12.9
**Branch:** `ng-psp-mode`

## Plan

B2: prove the file is the walk — gather a sample to a psp, read it back, compare record for
record, field for field, against the same sample walked in memory, on the fixtures **and** on
one real CRAM slice. B3: the format-level byte-identity re-made through the gatherer on real
reads — the same sample gathered twice gives identical files.

## Assumptions

- **The real-data oracle's in-memory side is a second gatherer**, not the bare direct-mode
  chain: `generic_path_generators` and `WalkReference` are `pub(crate)`, so an example cannot
  build the bare chain without widening visibility, and that a gatherer *is* the bare chain
  is already pinned at fixture scale by `the_gatherer_yields_what_the_direct_walk_yields`.
  The equality chain — bare chain == gatherer (fixture test), gatherer's stream == file
  (this oracle) — covers the plan's claim without a new `pub`.
- **B3's "worker-count" is degenerate by ruling**: the walk runs at concurrency one (owner,
  2026-09-03), so the schedule-invariance oracle is gather-twice identity. With the caller's
  provenance timestamp fixed, identity is asserted over the whole file — stronger than
  §12.1's "identical but the timestamp".
- **The real slice lives in `tmp/tomato_slice/`** (one CRAM + index + the repeat catalog,
  83 MB, copied from the main checkout's `benchmarks/tomato1/crams/` — the worktree's
  container cannot see the main checkout, and the benchmark data is not in git).

## Changes made

- `examples/ng_psp_gather_oracle.rs` (new): reference + catalog + BED assembled exactly as
  `ng_call_cohort_end_to_end` does, then per sample: gather → file, walk again in memory,
  compare header and every record; `NG_TWICE=1` gathers again and compares whole files as
  bytes. Non-zero exit on the first disagreement, naming the record index and field group.
- `src/ng/run/gatherer.rs`, three tests:
  - `analysed_but_empty_ground_round_trips` (§12.9) — ground analysed over chr2, reads only
    on chr1: zero records written, the file opens whole, and its header still names the
    analysed ground.
  - `tract_bearing_ground_round_trips_the_walk` — the review's carry-forward: a segmentation
    with a declared AT tract and a read crossing it; the guard asserts tract-kind records
    are actually on the walk (they are), and the file reads back equal.
  - `gathering_the_same_sample_twice_gives_identical_bytes` — B3 at fixture scale, whole
    files compared.

## Validation results

All in the container, on the committed tree:

- `cargo fmt --check` clean; `cargo clippy --all-targets --all-features -- -D warnings` exit 0.
- `cargo test --lib 'ng::run::gatherer'` — **16 passed** (was 13).
- The real-data oracle, release build, one tomato accession (`SRR7279481.p1.bench.cram`,
  read group sample `SRS3394712`) over the first two intervals of
  `benchmarks/tomato1/regions.bed` (200 kb of SL4.0, about three reads a position), quoted
  from the harness's own output:
  - 183,807 records, 4 blocks, 948,689 bytes, gathered in 0.32 s (setup 2.75 s);
  - **1,217 of the walk's records are repeat tracts** — the tract path is on the oracle;
  - 2,713 of 2,991 segments walked by a filled generator (the rest: bundles/satellites,
    refused by kind);
  - **the file IS the walk**: all 183,807 records equal field for field, header equal;
  - `NG_TWICE=1`: the second gather **byte-identical, all 948,689 bytes**.

## Tradeoffs and follow-ups

- The run-level shared test-fixture module (B1 review's carry-forward "lands with B2") is
  **still owed**: B2's two tests reuse gatherer.rs's local fixtures rather than hoisting
  them, because the hoist rewrites callers.rs's and walker.rs's test modules and belongs in
  its own commit. Carried to the Checkpoint B note.
- The oracle holds every record of the in-memory walk at once (183,807 records ≈ the walk's
  working set); fine for a 200 kb probe, not a whole-genome tool — it is an oracle, not a
  pipeline.
