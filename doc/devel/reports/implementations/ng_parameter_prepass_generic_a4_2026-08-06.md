# ng step 4, generic path — A4: the depth ladder

**Date:** 2026-08-06
**Plan:** [parameter_prepass_generic.md](../../ng/impl_plan/parameter_prepass_generic.md), Milestone A step A4 — **own commit, do not bundle**
**Design authority:** [arch](../../ng/arch/parameter_prepass_generic.md) §2.2 · [spec](../../ng/spec/parameter_prepass_generic.md) §4 · [research note](../../ng/research/parameter_estimator_experiments_2026-08-06.md) §4.3

---

## 1. Plan

`DepthBin` and `DepthBinEdges` — which depths share a bin, and where each bin's row
starts in a histogram's flat cell vector. `bin_for`, `row_start`, `depth_range`,
`cell_count`, `bin_count`, `max_depth`.

**Why this is its own commit.** The plan isolates six steps whose failure is silent,
and this is the first. The ladder's edges are a **correctness** parameter, not the
memory-only knob an earlier draft of the architecture assumed: at the same cap of 124,
sixteen bins bias the fitted error rate by 0.55 rungs and the homozygous-non-reference
rate by 1.8%, against 0.05 rungs and 0.3% for twenty. A cap of 300 on sixteen bins
costs 1.04 rungs and 8.0%. Every one of those is a plausible number, and nothing
downstream would show any of them — so if a parameter moves later, `git bisect` has to
be able to land on one commit.

## 2. Assumptions

1. **`bin_for` is total rather than panicking above the cap.** The architecture does
   not say what a depth above `max_depth()` answers. It returns the last bin. The cap
   is enforced by the subsampling that runs *before* a site reaches a histogram
   (Milestone C2), so a panic here would be a second guard on the same invariant, and
   the one that fired would be the less informative of the two. Pinned by test.
2. **The ladder is generated, not tabulated.** `bin_tops` is computed from
   `EXACT_DEPTH_LIMIT`, `MAX_BINNED_DEPTH` and `DEPTH_BIN_COUNT` by the same rule the
   research harness used — `EXACT_DEPTH_LIMIT · r^k` rounded, with `r` the ratio
   reaching the cap in the bins left over — rather than being written out as eleven
   literals. The eleven literals are in the *test*, so the generator is checked against
   the measured ladder rather than being the only statement of it.
3. **`DepthBinEdges` owns `row_starts`.** The architecture says the offsets are "a pure
   function of the edges, so they are computed once per run and shared". They are a
   field, computed in `new()`, rather than a method that recomputes.

## 3. Changes made

**New — [src/ng/parameter_estimation/generic/depth_bins.rs](../../../../src/ng/parameter_estimation/generic/depth_bins.rs)**

Three constants stating the ladder's shape — `EXACT_DEPTH_LIMIT = 8`,
`MAX_BINNED_DEPTH = 124`, `DEPTH_BIN_COUNT = 20` — each with what the research note
measured for the alternatives. `DepthBin(u16)`, and `DepthBinEdges` holding `bin_tops`
and `row_starts` with the six accessors.

`depth_range` returns a `RangeInclusive<u32>` rather than a pair, because its sole
consumer is a row width where an off-by-one silently mis-sizes the table — inclusivity
belongs in the type.

**A deviation from the plan's file list, recorded.** The plan's A1 names four files
under `generic/`; this is a fifth, `depth_bins.rs`. The alternative was to put the
ladder in `histogram.rs`, which the architecture's own module table describes as "the
cell table". The binning rule is not the table — it is shared by every table in the
run, by `Arc`, precisely so that two of them cannot be binned differently. Putting the
shared rule inside one of its consumers would have made that relationship unreadable.

**Changed — [generic/mod.rs](../../../../src/ng/parameter_estimation/generic/mod.rs)** — `pub mod depth_bins;`.

## 4. Tests added

Seven, and the first three are the plan's stated oracle:

- `the_widening_bins_top_out_at_the_measured_depths` — the eleven tops are exactly
  `10, 13, 17, 22, 28, 36, 46, 59, 75, 97, 124`.
- `the_ladder_holds_583_cells` — 45 in the nine exact bins, 538 in the eleven widening
  ones. This is the number spec §9 prices the accumulator's memory against.
- `bin_for_is_monotone_and_total_over_every_depth_the_ladder_covers` — a deeper site
  never lands in an earlier bin, every bin is reachable, and the bin a depth answers is
  the bin whose range contains it. A ladder that skipped a bin would leave a histogram
  row permanently empty, which no fit would report.

And four more that the plan does not name but the shape needs:

- `the_bottom_of_the_ladder_is_one_bin_per_depth` — depths 0–8 are never merged.
- `the_bin_ranges_partition_the_depths_from_zero_to_the_cap` — consecutive ranges, no
  gap, no overlap. Stated separately from `bin_for` because a rule that agrees with
  itself is not the same as a rule that covers everything.
- `each_row_is_as_wide_as_its_bins_deepest_alternative_count` — the rows sit end to
  end; an off-by-one would have one bin's cells reading into the next bin's row.
- `a_depth_above_the_cap_answers_the_last_bin` — assumption 1, pinned.

## 5. Validation results

All in the container.

| command | result |
|---|---|
| `cargo fmt --check` | clean |
| `cargo clippy --all-targets --all-features -- -D warnings` | clean |
| `cargo doc --no-deps --lib` | 12 unresolved links, all pre-existing; none in this file |
| `cargo test --lib ng::parameter_estimation` | 13 passed, 0 failed |

**The oracle was shown to bite, against the three ladders the research note measured
as worse.** Each was applied as a one-constant mutation and the suite re-run:

| mutation | what the note measured for it | tests failing |
|---|---|---|
| `DEPTH_BIN_COUNT = 16` | 0.55 rungs, 1.8% | 2 |
| `MAX_BINNED_DEPTH = 300` | 0.190 rungs, 0.88% (at 20 bins) | 2 |
| `EXACT_DEPTH_LIMIT = 4` | 0.15–0.21 rungs even at 3 reads | 2 |

**Corrected after review, twice over.** This table first read 2 / 3 / 3; those came from
a `grep -c FAILED` that also matched the `test result: FAILED.` summary line. Two
reviewers independently re-ran all three and measured **2 each, and the same two every
time** — `the_widening_bins_top_out_at_the_measured_depths` and
`the_ladder_holds_583_cells`. The `MAX_BINNED_DEPTH = 300` row also quoted the
sixteen-bin cost against a twenty-bin ladder.

**The correction changes what the table means, which matters more than the number.**
The other five tests are written *in terms of* the three shape constants, so their
expectations move with a mutation and they cannot detect a constant edit at all. They
check that whatever ladder is configured is internally consistent — a different
question from whether it is the ladder that was measured. **The oracle against a moved
parameter is two tests**, the two carrying hard literals. That is now said in the doc
comment of the first of them, where the next reader meets it, and a third test gained a
literal `max_depth()` assertion so it joins them.

The commit message of `54378b9` carries the original 2/3/3 figures and cannot be
corrected without rewriting history.

## 5a. Review outcome

Three agents, seven categories, each in its own worktree: 0 Blocker, **3 Major**, 10
Minor, 4 Nit; 16 applied, 1 deferred, 0 disputed. See the
[review](../reviews/ng_parameter_prepass_generic_a4_2026-08-06.md) and the
[fixes applied](../reviews/fixes_applied_2026-08-06_v2.md).

**Two of the three Majors were real defects in this file, and only mutation found
them.** The guard against a degenerate ladder fired zero times at the adopted constants
— deleting it left all seven tests green — and its two clamps were ordered so that when
it *did* fire it produced the empty bins it claimed to prevent. And `row_start` answered
`583` for a bin index one past the end, where `depth_range` panicked on the same input;
583 is `cell_count()`, a number shaped exactly like a row offset.

Both are fixed, and the fix's own first attempt reproduced the defect it was fixing:
reordering the clamps makes strict increase provable by construction, so the assertion
proposed alongside it was itself unreachable. The tests caught that.

After fixes the module holds 20 tests (was 7) and the suite is 2,920 passed / 1 failed
(pre-existing) / 5 ignored.

## 6. Tradeoffs and follow-ups

- **The generator is checked against the ladder, not the reverse.** If a future change
  wants different edges, the eleven literals in
  `the_widening_bins_top_out_at_the_measured_depths` are what must be re-measured and
  re-stated — not adjusted to match whatever the generator then produces. The test's
  doc says so.
- **Nothing consumes the edges yet.** `DepthAltHistogram` takes them by `Arc` in B2,
  and `merge` proves two tables share one rule by `Arc::ptr_eq` in B4.
- **`bin_for`'s partition-point search is O(log 11).** The exact region — which is
  where tomato's cohort lives, 97 sites in 100 — short-circuits before it.
