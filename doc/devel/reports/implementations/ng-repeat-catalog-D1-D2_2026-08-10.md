# ng repeat catalog — D1 + D2: the segmentation, and the differential that guards it

*Implementation report, 2026-08-10. Plan:
[`impl_plan/repeat_catalog.md`](../../ng/impl_plan/repeat_catalog.md) steps **D1** and **D2**.
Design: [`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) §5.1, §5.3, §10.1 and
[`arch/repeat_catalog.md`](../../ng/arch/repeat_catalog.md) §2.3.*

## Plan

`genome_segments` derives the genome's segmentation from the file, and the differential asserts it
equals what a live scan produces. One commit, because the guard is the step.

## Assumptions and deviations

1. **Step 3's own functions do the work.** `prefilter`, `split_bundles`, `is_compound`,
   `coverage_runs`, `resolve_features` and `fill_generic_gaps` became `pub(crate)` and are called
   as-is; the only thing rewritten is `finish_locus`, and only because its three base-reading values
   are stored. That is what makes "identical, not close" reachable.
2. **`genome_segments` streams one contig at a time.** A contig's rows and its segments are in
   memory; the genome's never are.
3. **The criteria check is eager**, at the call, before a row is read — a refusal is a fact about the
   file and the policy, not about a row.
4. **Scoring-weight equality is a separate call, `check_scored_with`.** A reader that only reads has
   no `ScanParams` to offer; one that also scans must not mix the two, and calls it.

## Changes made

- [`src/ng/repeat_catalog/segments.rs`](../../../../src/ng/repeat_catalog/segments.rs) —
  `segments_of_contig`, `admit` (classify's gates over stored fields) and `finish_from_row`.
- [`src/ng/repeat_catalog/reader.rs`](../../../../src/ng/repeat_catalog/reader.rs) —
  `genome_segments`, the `GenomeSegments` iterator, `check_serves`, `check_scored_with`.
- [`src/ng/region_typing/`](../../../../src/ng/region_typing/) — six items widened to `pub(crate)`;
  no logic changed.

## Tests added

Four integration tests in
[`tests/ng_repeat_catalog_differential.rs`](../../../../tests/ng_repeat_catalog_differential.rs),
over a six-contig fixture carrying clean tracts of three periods, two tracts inside the bundle
radius, a 1.2 kb array, a 180 bp array, an interrupted tract and tracts hard against both contig ends.

- **The derived segmentation equals the scanned one at every policy** — six policies, including ones
  differing on every bounded axis and both unbounded ones.
- **The derived loci are the scanned loci**, compared by coordinates, motif and purity, with a guard
  that the fixture actually produces loci.
- **A more permissive reader is refused eagerly**, naming the axis and the value; and the mirror
  case, that the unbounded axes are served.
- **The fixture drives the purity gate and the satellite cap** — see below.

## What mutation testing found, and it mattered

The first version of all three comparisons passed. **Removing the purity floor from the derived path
left them all green**, because no locus in that fixture had a purity between the base policy's 0.8
and the strict policy's 0.95 — the interrupted tract scored 0.96, so the "stricter purity" policy
discriminated nothing. The fixture now uses `(CAG)3 CTG` three times over, which scores 0.92, and a
fourth test asserts the discrimination directly: some locus must sit between the two floors, the
strict policy must drop at least one, and the small satellite cap must produce a satellite the large
one does not. With that fixture the purity mutation fails two tests.

A second mutation — dropping `.max(1)` from the left flank clamp — stays green, and that is honest
rather than a gap: the file's own 15 bp flank floor excludes every tract nearer the contig start than
the default bundle radius, so the clamp cannot fire for a catalog row. The code keeps it to mirror
`finish_locus` exactly, and says so.

A third — an off-by-one in the detected-span conversion — fails two tests, which is the silent
failure this step was isolated for.

## Validation

In the dev container: `cargo fmt` clean; `cargo clippy --lib --tests --all-features -- -D warnings`
clean; `cargo test --lib` **2,924 passed**; the differential **4 passed**.

**Pre-existing and untouched:** `cargo clippy --all-targets` still fails on
`examples/ng_inbreeding_harness.rs` (`076cb5e9`).

## Tradeoffs and follow-ups

- `admit` matches pre-screened intervals back to their rows with a linear `find` per interval, which
  is quadratic in a contig's row count. Correct but not the shape for a genome-scale contig; a sorted
  merge is the fix, and E2's measurement is what should decide whether it matters.
- The comparison excludes regions within the flank floor of a contig's ends — the one stated place
  the file holds less than a scan. The `chr_edges` contig exists so that exclusion is exercised
  rather than hypothetical.
