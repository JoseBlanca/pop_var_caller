# ng repeat catalog — A3 + A4: the shared arithmetic, and one interval becomes one row

*Implementation report, 2026-08-10. Plan:
[`impl_plan/repeat_catalog.md`](../../ng/impl_plan/repeat_catalog.md) steps **A3** and **A4**.
Design: [`spec/repeat_catalog.md`](../../ng/spec/repeat_catalog.md) §1, §3.1, §3.2 and
[`arch/repeat_catalog.md`](../../ng/arch/repeat_catalog.md) §2.2, §5.*

## Plan

**A3** shares `classify`'s tract arithmetic with the catalog builder; **A4** uses it to turn one
detected interval into one row. Two commits, as the plan asks: A4 is the module's one silent-failure
site, so it lands alone with its tests.

## Assumptions and deviations

1. **A3 lifted visibility, not code.** The plan says the trim, the motif slice and the purity
   recomputation are "lifted into `pub(crate)` helpers so the builder and `finish_locus` share one
   implementation". `minimal_trim`, `recompute_purity` and `upper` were **already** exactly those
   helpers, called in that order by `finish_locus`; the only thing missing was that the builder could
   not see them. So A3 changed three signatures to `pub(crate)` and documented why each is shared. No
   body moved, which is the strongest available evidence that behaviour is unchanged — and the 113
   `region_typing` tests are green before and after.
2. **A4 computes the purity but not the compound-motif gate.** `classify` rejects a compound motif
   (`ATAT` = `(AT)²`) *before* bundling, and the check reads only the stored motif, so the reader
   re-runs it. Nothing is lost and the builder stays criteria-free.
3. **The flank is measured against the cut tract**, falling back to the detected span when there is no
   clean cut. The spec says "15 bp of sequence beside it"; the locus is the cut tract, so that is what
   a read would be anchored around.
4. **`row_for_interval` takes the whole `StrRepeatCriteria`**, though it reads only the copy floors and
   the flank. Passing the two fields separately would let a caller mix a floor from one policy with a
   flank from another.

## Changes made

- [`src/ng/region_typing/segment_criteria.rs`](../../../../src/ng/region_typing/segment_criteria.rs) —
  `upper`, `minimal_trim` and `recompute_purity` are `pub(crate)`, each documenting that it is shared
  with the catalog builder and that it reads only the tract and its motif, never the criteria. **That
  last property is what makes storing the results sound**: a different copy floor or purity floor
  changes which tracts survive, never where the cut falls.
- [`src/ng/repeat_catalog/row.rs`](../../../../src/ng/repeat_catalog/row.rs) — `row_for_interval`,
  `clears_detected_copy_floor`, and `RowRejection` (four admission outcomes, each documenting when it
  fires).

`row_for_interval` runs `classify`'s own order: motif → copy floor **on the detected span** → cut →
purity over the cut → flank. It applies no purity floor, no satellite cap and no bundling.

## Tests added

10 new tests in `row.rs` (23 in the module, 2,898 in the lib, all green).

- **The conversion, checkable by hand**: a `(CAG)6` tract at offsets 40..58 comes out as 41..=58,
  18 bases, 6 copies, motif `CAG`.
- **A perfect tract trims to itself** at purity 1; **a ragged tract cuts back to whole copies**, the
  cut lies inside the detection, and its length is a whole multiple of the period.
- **A tract with no clean cut is still a row**, with `trimmed`, `purity` and `stratum()` all absent
  and its detected span intact — the rule that keeps a neighbour's bundling identical to a live scan.
- **A tract below the floor after trimming is kept**: `(CAG)3 CAT` is 4 detected copies (the period-3
  floor) and 3 trimmed ones (below it), and the row stands. This is the measurement-site rule stated
  as a test, and the fixture asserts it really does fall below the floor, so it cannot quietly stop
  proving anything.
- **A tract below the floor on the detected span is rejected.**
- **15 bp of flank is kept and 14 is not**, on both contig ends — the pair that pins the boundary
  rather than merely exercising it.
- **A soft-masked tract yields an upper-case motif** (tomato SL4.0 carries 227,170 lowercase bases).
- **The pre-screen agrees with the row builder** about the copy floor, in both directions.

## Validation

In the dev container:

- `cargo fmt` — clean; `cargo clippy --lib --tests --all-features -- -D warnings` — clean.
- `cargo test --lib` — **2,898 passed, 0 failed, 5 ignored**.

One failure found and fixed during the run, in a test rather than the code: the first
"below the floor after trimming" fixture was an `AT` repeat whose detected span began mid-copy, so
the motif came out as `AG` and no clean cut existed at all. Replaced with `(CAG)3 CAT`, which loses a
copy to the trim as intended.

**Pre-existing and untouched:** `cargo clippy --all-targets` still fails on
`examples/ng_inbreeding_harness.rs` (three `dead_code` errors, committed in `076cb5e9`).

## Tradeoffs and follow-ups

- `clears_detected_copy_floor` has no caller outside its test yet; C2's builder is its consumer. It
  exists now because the property it names — the floor is measured where `prefilter` measures it — is
  what keeps bundling identical, and a named function is what a test can hold.
- Rejections are counted, not logged, and the tally surfaces in C2's builder and E1's CLI output.
