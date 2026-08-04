# ng — the alignment cursor, Milestone D1+D2: fixes applied

*Against [the review](../reviews/ng_alignment_cursor_d_2026-08-02.md). One row per finding,
including the ones deliberately not fixed.*

## Applied

| finding | fix |
|---|---|
| **Blocker — D2 could be switched off with the suite green** | `PileupGenerator::cursor_counts()`, summing retired chromosomes' cursors and asking the live one, reached through a new `PileupWalker::reads()`. Two tests: `the_cursor_is_kept_across_regions_rather_than_minted_for_each` (asserts one jump and two reuses over three ascending regions, and that reads come back replayed) and `cursor_tallies_are_taken_from_a_chromosome_before_it_is_retired`. Kept **off** `PileupGeneratorCounts`, which the dump tools print verbatim against byte-identical baselines. |
| **Major — abandoned region leaked a record** | `a_region_abandoned_half_walked_leaks_no_record_into_the_next_one`: take one locus, walk away, and assert the next region emits nothing anchored before its own start. |
| **Major — failure blamed on the chromosome's first region** | `a_failure_names_the_region_it_happened_in_not_the_chromosomes_first`: a clean region, then one whose read fails preparation, asserting `error.region()`. |
| **Major — `RegionReadSource` did not require replay** | The requirement is now in the trait method's contract, with what breaks without it and why the walker depends on it. |
| **Major — `fold_region_walk`'s stale justification** | Rewritten to say what actually makes the plain sum safe, and what swapping `begin_region` for `reset` would do. Same for `silent_exits`'s field doc and the delta fold's `debug_assert` message. |
| Minor — `IngestError::Cursor`'s "every `CursorError` names the path" | Corrected, and the `DuplicateReadAcrossFiles` collision with the enum's own top-level variant recorded as F's to resolve. |
| Minor — the box made the generator unconditionally `!Send` | `Box<dyn FnMut() -> R + Send>`. Costs nothing; the fan-out gives each worker its own generator. |
| Minor — `AfterFailure` reaches `OpenReadQuery` | Documented on the variant: it is the other category by the split's own logic, and unreachable through this generator because the `failed` latch fires first. |
| Minor — `enter_chromosome`'s redundant `chain_ids.reset()` | Kept, comment corrected: `begin_region` does it one call later on the ordinary path, so what this covers is the **failure** path, where the allocator sits in `self.chain_ids` and no `begin_region` is coming. |
| Minor — the ordering rationale in `move_to_region` was untrue | Corrected: both steps must precede the re-anchoring peek, which of them leads is not load-bearing, and the reposition leads because a source that cannot move leaves the walker untouched. |
| Minor — a chromosome test discarded its middle walk | Now asserted non-empty. |
| Minor — the fallible move set state before it could fail | `region` and `done` are adopted **after** the cursor move succeeds. |
| Minor — eleven stale doc sites | Fixed in `pileup/mod.rs`, `locus_generation/mod.rs`, `read/input/mod.rs` (×3), `benches/ng_generic_pileup_perf.rs` (×3), `examples/ng_generic_loci_dump.rs`, `examples/ng_generic_walk_probe.rs`. |
| Minor — `genome_walk.rs`'s "changed in one respect" header | Replaced with the divergence list, A0 through D1. |

## Measured rather than fixed

**The reference-accessor Blocker.** Measured on a synthetic contiguously-covered 20 Mb contig at
30×, four spans, fixture build isolated in its own process. The delta does not scale with contig
length; +3.4 MB on a 25.6 MB baseline at the largest span, and 1.32× faster. Refuted on
magnitude — table and reasoning in the implementation report.

What the measurement found instead — a ~1 byte-per-base term on **both** sides, from
run-lifetime accessors nothing can evict — is pre-existing and is raised for the owner at
Checkpoint D.

## Not fixed, and why

- **Surviving mutations on defensive lines** (`pending.clear()`, the three `begin_region` field
  resets, `next_read_id`/`by_read_id`, the dead error arm). Each was probed; none has a
  reachable failure that could be written without a contrived fixture. Listed in the review under
  *Known and accepted* so the next reader finds them named rather than discovering them again.
- **The spec/arch fold-in for `locus_generation_pileup.md`.** This skill does not edit design
  documents, and Milestone D's plan does not include a fold-in. Raised at Checkpoint D.
- **The generic path's unreachable per-read-group drop tallies.** Milestone F's inventory.

## Validation after the fixes

`cargo test --lib ng::` **1,559 passed**. `cargo clippy --all-targets --all-features -- -D
warnings` clean. `cargo fmt --check` clean. `cargo test --examples` green (34 suites).

Anchors re-run after every fix: `ng_generic_walk_probe` chr21 prints `loci=236081
observations=251786 reads_admitted=54709`; `ng_generic_loci_dump` and `ng_ssr_loci_dump` are
byte-identical to `ee0c94b`'s binaries.

Six mutations re-run and each still killed, plus three new ones for the new tests: minting a
cursor per region, dropping the stream's region reassignment, and discarding retired cursors'
tallies.
