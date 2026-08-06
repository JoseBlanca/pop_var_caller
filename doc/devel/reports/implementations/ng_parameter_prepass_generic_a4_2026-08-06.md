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
| `MAX_BINNED_DEPTH = 300` | 1.04 rungs, 8.0% (at 16 bins) | 3 |
| `EXACT_DEPTH_LIMIT = 4` | 0.15–0.21 rungs even at 3 reads | 3 |

All three restored, 7 passed.

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
