# window coverage — A1: production's sliding window, transcribed into ng

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone A, step A1
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.3, §3.4, §3.5, §3.6, §7
**Review:** [ng_window_coverage_a1_2026-09-06.md](../reviews/ng_window_coverage_a1_2026-09-06.md)
**Branch:** `ng-window-coverage`

## The answer

**The centred sliding window that production computes now exists inside ng, and a differential
against the original says the copy computes the same numbers.** Fed 200 pseudo-random streams —
jumping and repeating positions, `N` bases in both cases, contig changes, window widths from one
base to 600 over contigs of up to 2,000 covered positions — the two accumulators emit the same
windows at the same centres, with **447,581 mean depths and 447,581 GC fractions equal bit for
bit**, and identical histogram cells. Nothing else in the run changed: this step adds a module
and one `mod` line, and no existing code calls it yet.

## Plan

The plan's A1 asks for two files under `src/ng/window_coverage/` holding the configuration, the
emitted pair, the histogram type and the accumulator, with production's `sliding_*` tests
transcribed and green, and without the fixed-tile accumulator or the heterozygosity fields.
Delivered as three files — the third is the differential described above, which the plan does
not ask for.

## Assumptions and deviations, all minor and all recorded

The plan names the items to transcribe by production's identifiers; the spec's §3.6 gives the
shapes ng is to end with. Where the two differ, the spec won, because A2 and A3 extend exactly
these types and a rename in the middle of that would be churn:

1. **`SlidingWindowCoverageAccumulator` → `WindowCoverageAccumulator`, `CoverageBinScheme` →
   `WindowCoverageConfig`** (spec §3.6's names). The "Sliding" prefix distinguished production's
   two accumulators; ng brings over only one, so it distinguishes nothing here.
2. **The emitted pair carries no coordinate.** Production's `WindowCoverage` holds `chrom_id`
   and `pos`; spec §3.5's holds only `gc_fraction` and `mean_depth`, and §3.6 has `pop_ready`
   return `(GenomePosition, WindowCoverage)`. Followed, because the cache keys the pair by its
   centre anyway and a locus would otherwise repeat a coordinate every consumer already has.
3. **`PartialEq` on the pair compares bit patterns, not values.** No `NaN` is produced until A2
   adds the floor, so this is written ahead of its need: the type is the one that will carry
   `NaN` as *absent*, and a derived `PartialEq` would make an absent window unequal to itself
   the moment the floor lands (spec §6 trap 4). Production's `LocusWindowCoverage`
   ([types.rs:233](../../../../src/var_calling/types.rs)) does the same.
4. **Two histogram fields dropped.** Spec §7's reuse map names `callable_positions` (and
   "heterozygosity", which is not a field of this type); `n_skipped_tiles` is dropped on the same
   argument without the spec naming it, because the sliding model can only ever write `0` there.
   Neither is read by the coverage model fit — checked: `src/paralog/coverage_model.rs` reads
   `window_bp`, `gc_bins`, `depth_bins`, `depth_bin_width` and `counts` and nothing else.
   `callable_positions` has exactly one in-crate reader, `src/var_calling/diversity.rs:223-235`,
   which uses it as a cohort diversity denominator — a number ng obtains elsewhere. Production's
   `n_positions` becomes `windows_folded`, which is what it counts.
5. **`finish` returns the histogram, not `Option<CoverageByGcHistogram>`.** Spec §3.6's `Option`
   exists because a sample that finalises no window has no median to set its depth bin width
   from — that is A3's business, and `finish`'s doc says the shape changes there. At A1 the width
   is a configuration constant and a histogram always exists.
6. **ng's coordinate types replace production's two `u32`s**: `ContigId` and the 1-based
   `Position`, which is a `u64`. The window arithmetic is unchanged; only the saturation comment
   moves, since a centre within half a window of `u64::MAX` is unreachable for a different
   reason than it was at `u32::MAX`.

**Beyond the plan, and deliberate:** the differential in `production_parity.rs`. A1's whole
claim is that the copy is faithful, and the transcribed unit tests check it against
hand-computed means on streams they were written for. The differential checks it against the
thing it was copied from, on streams neither was written for. It is `#[cfg(test)]`, so nothing
shipped depends on `src/sample_summary/`; ng reads production as a test oracle in two other
places already (`ng/scanner_parity.rs`, `calling::genotype_table_parity`).

**Added during the review's fix pass, and each recorded there with its finding:** two settled
constants (`WINDOW_BP`, `GC_BINS`) so that C3's look-ahead and C2's configuration derive their
numbers from one place; `saturating_add` on the histogram cell counter, which above about 4.3 Gbp
of reference is the difference between a cell that under-reports and one that reports a
near-empty cell where the single-copy peak is; an exhaustive destructure of the configuration in
`finish`, so A2's and A3's new fields cannot silently fail to reach the histogram; `#[must_use]`
on `finish`; three field renames (`half` → `half_window_bp`, `buf` → `positions`, `count` →
`summed_positions`); and two narrowed visibilities.

## Changes made

- **[`src/ng/window_coverage/mod.rs`](../../../../src/ng/window_coverage/mod.rs)** — the module
  doc (one computation, two consumers), `WINDOW_BP` and `GC_BINS`, `WindowCoverageConfig` with
  its validity assertion and cell-index arithmetic, `WindowCoverage`, `CoverageByGcHistogram`.
- **[`src/ng/window_coverage/accumulator.rs`](../../../../src/ng/window_coverage/accumulator.rs)**
  — `WindowCoverageAccumulator`: the two-pointer buffer, the finalisation frontier, the contig
  reset, the histogram fold, and 25 tests.
- **[`src/ng/window_coverage/production_parity.rs`](../../../../src/ng/window_coverage/production_parity.rs)**
  — the differential against `src/sample_summary/coverage.rs`.
- **[`src/ng/mod.rs`](../../../../src/ng/mod.rs)** — one `pub mod window_coverage;` line.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched; the freeze
holds, checked with `git diff` over those three trees.

## Tests

**26 in the module: eleven transcribed from production, fourteen new unit tests, and the
differential.**

The eleven keep production's fixtures and assertions, with the tuple return and ng's coordinate
types substituted, and two small changes: assertions on the two dropped histogram fields become
`windows_folded`, and `sliding_empty_stream_is_empty` gains an all-zero-counts assertion
production does not have. One is renamed —
`sliding_n_position_excluded_but_advances_frontier` → `sliding_n_position_is_excluded_from_every_window`
— because its fixture cannot check the second half of the old name (every covered depth in it is
10, so deferring a finalisation changes nothing); the frontier now has a test of its own.

| test | what it pins |
|---|---|
| `sliding_ramp_means_match_hand_computed` | ten centred means over a linear ramp, against hand arithmetic |
| `sliding_uniform_all_one_cell` | a flat stream puts all 20 windows in one histogram cell |
| `sliding_n_position_is_excluded_from_every_window` | an `N` at depth 999 enters no mean |
| `sliding_window_does_not_span_contigs` | two contigs at depth 10 and 20; no window averages across the boundary |
| `sliding_is_deterministic` | kept for parity; cannot fail until a clock or a map order is introduced |
| `sliding_odd_window_uses_floor_half` | `window_bp = 5` gives `half = 2`, the documented truncation |
| `sliding_empty_stream_is_empty` | no windows, an all-zero histogram |
| `sliding_window_bp_one_emits_self_only_window` | `half = 0`, where the two-pointer extents collapse |
| `sliding_large_gaps_yield_singleton_windows` | positions 10, 100, 1000: no stale sum crosses a gap |
| `sliding_single_covered_position_emits_one_self_window` | the divisor-is-one path |
| `sliding_finish_returns_the_undrained_tail` | a caller that never drains still gets every window |
| `finish_echoes_the_configured_bin_scheme` | the four fields a consumer reads a cell by |
| `mean_depth_on_the_top_bin_edge_lands_in_the_overflow_column` | the fit's own rejection guard |
| `gc_fraction_of_one_saturates_into_the_last_gc_bin` | the GC clamp |
| `lowercase_reference_bases_are_read_the_same_as_uppercase` | soft-masked reference, GC and `N` |
| `a_repeated_position_emits_one_window_per_observation` | the boundary of the non-decreasing order rule |
| `a_centre_closes_when_the_frontier_reaches_its_right_edge` | the `<=` frontier, visible only under a repeat |
| `every_folded_window_lands_in_exactly_one_histogram_cell` | cell counts sum to `windows_folded` |
| four `new_panics_on_*` | each of the configuration's four asserts |
| two `observe_panics_in_debug_on_*` | the coordinate-order guard, both axes |
| `an_absent_window_equals_itself_and_signed_zeroes_differ` | the bit-pattern equality, ahead of A2 |
| `the_transcription_matches_production_on_streams_neither_test_was_written_for` | the differential |

**Fourteen of the fifteen new tests exist because a mutation survived the twelve the step
started with.** Nineteen mutations were run in the review; seven survived, none of them because
the mutation changed nothing. Seven mutations were then re-run on the fixed tree **with the
differential excluded**, so that the unit tests have to carry them: the histogram's
`depth_bin_width` echo forced to a constant, the overflow column merged into the last regular
bin, the frontier's `<=` narrowed to `<`, case-sensitive GC, case-sensitive `N`, a derived
`PartialEq`, and the `window_bp` assert removed. **Each fails exactly one test** — 24 passed, 1
failed, every time.

The differential's own discrimination was established the same way: pulling the window's left
edge **in** by one base — `centre - half_window_bp` to `centre - (half_window_bp - 1)`, which
narrows the window rather than widening it — makes it fail on a value at the first seed, and
fails three of the eleven transcribed tests. (At `window_bp = 1` the half-window is zero and that
literal mutation underflows rather than producing a wrong window.) The count of windows compared,
447,581, is asserted inside the test, so a change that quietly stops the streams producing
windows cannot leave it passing on nothing.

## Validation results

Run in the container (`./scripts/dev.sh`), on this worktree:

- `cargo test --lib --all-features` — **6,301 passed, 0 failed, 15 ignored**. The 26 above are
  the increase over `main`'s 6,275, measured by commenting out the `mod` line and re-running.
- `cargo test --lib --all-features ng::window_coverage` — 26 passed, 0 failed.
- `cargo clippy --lib --all-features` — **no diagnostic in `src/ng/window_coverage/`**. Three
  `needless_lifetimes` warnings remain elsewhere in the crate; all three predate this branch
  (`cohort_merge/build.rs:820`, `build.rs:893`, `cohort_merge/serial.rs:67`), measured on a tree
  with this module removed.
- `rustfmt --check` clean on all three module files.

**Four checks are red on `main` before this work and are left alone**, because three of them sit
in `src/ng/run/cohort_merge/`, where another branch is working: nine files fail
`cargo fmt --check`; the three clippy lints above; the integration test
`a_contaminants_reads_at_a_tract_are_not_called_as_a_second_allele` fails; and the example
`ng_candidate_selection_probe` does not compile against the current `CohortObservation::over`.
Each was measured with this module taken back out of the tree.

**The standing calling oracle is green and was not re-run for this step.** Six tomato accessions
over the first two 100 kb intervals of `benchmarks/tomato1/regions.bed`, called from the CRAMs
and from stored files: 2,311 records, the two VCFs identical apart from `##commandline` (sha256
`84ad19c2…`), and the parameters file beside each identical too. This step cannot have moved it —
nothing calls the new module — and the plan runs it either side of A3, which lands alone.

**One hazard worth naming for whoever measures next.** Restoring a file with `mv file.bak file`
carries the backup's older modification time back with it, and cargo then reuses a stale build
and reports the *previous* tree's test count with no compilation and no warning. It happened once
here, showing 6,275 where 6,301 was true. Restore with `cp`, or `touch` afterwards.

## Tradeoffs and follow-ups

- **The differential's histogram half has a limited life.** A3 fits the depth bin width per
  sample, at which point production's fixed-width histogram is no longer the same object; the
  window means and GC fractions stay comparable and are what the test is really for. Four
  behaviours it used to be the sole guard of now have unit tests, so the narrowing costs nothing.
- **`WindowCoverageConfig` still carries `depth_bin_width`.** It leaves the configuration at A3,
  when the width becomes a fitted per-sample number; `min_window_positions` (A2) and
  `depth_scale_windows` (A3) join it there.
- **Four items deferred with a home**, each recorded in the review's §7: promoting the
  coordinate-order guard to a release assert (C2, where the caller exists); a named `absent()`
  value (A2, the step that produces one); grouping the six per-contig fields into a sub-struct
  (A3, which reshapes the struct anyway); and widening the histogram's cell counter to `u64`
  (D3, priced against measured memory).
- **A number in the spec does not add up, and it is not a code defect.** Spec §5 budgets "one
  12-byte entry per held position per sample" for the ready deque; the entry spec §3.6 actually
  specifies — `(GenomePosition, WindowCoverage)`, with a `u64` position — is 24 bytes. At the
  spec's own 3,000 samples and half a 500-base window held per sample that is about 18 MB rather
  than 9. Raised with the owner at Checkpoint A; plan step D3 measures it.
- **No caller yet.** The module is dead code until C2 puts the accumulator in the cache's
  per-sample state. That is the plan's order — the algorithmic heart before the plumbing.
