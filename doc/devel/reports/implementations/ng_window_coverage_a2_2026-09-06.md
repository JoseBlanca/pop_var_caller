# window coverage — A2: the floor, so a window built from a handful of positions says nothing

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone A, step A2
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.3
**Review:** [ng_window_coverage_a2_2026-09-06.md](../reviews/ng_window_coverage_a2_2026-09-06.md)
**Branch:** `ng-window-coverage`
**Builds on:** [A1](ng_window_coverage_a1_2026-09-06.md)

## The answer

**A window is now allowed to refuse.** Production's window reports a number however few covered
positions it was built from — one is enough, and the resulting figure looks exactly as confident
as one averaged over five hundred. ng's holes are not scattered noise: a repeat-tract region that
emits no loci, an analysed-region edge, a stretch no read reached. So a window finalised over
fewer than `min_window_positions` **distinct** covered positions comes back **absent** — both
numbers `NaN` — and the filter skips that sample there, which its scorer already does for an
absent sample.

*Distinct* is the review's correction and not a wording choice: the accumulator's input rule is
non-decreasing, so one base can be observed more than once, and the first version counted records.
Probed, one base observed five times cleared a floor of five and was folded into the histogram —
the exact window the floor exists to silence. The divisor of the two means still counts every
observation, as production's does, so the differential is unaffected.

The value is **50 of 500 and it is provisional**, marked so in the code: plan step D1 measures the
distribution of covered positions per window on the tomato slice and on HG002 and sets it.

## Assumptions, and one the spec does not settle

**An absent window is emitted but not folded into the histogram.** The spec says a window under
the floor "emits an absent value" (§3.3) and separately that "every finalised window is folded
into a per-sample count matrix" (§3.4); read together those two sentences say to fold a `NaN`,
and that is not a thing the histogram can hold. Both casts of `NaN` come back `0`, so folding
would pile every window too sparse to speak into **one cell — the lowest GC bin's lowest depth
bin** — in exactly the samples and regions where coverage is thinnest. The yardstick the filter
fits would then carry a spike of windows that had no depth to report.

So the pair is emitted at its centre, absent, and the histogram never sees it. `windows_folded`
consequently stops equalling the sample's covered positions, and its doc says so. This is the one
choice A2 makes that the spec leaves open; it is cheap to reverse, and it is raised at
Checkpoint A.

## Changes made, after the review's fixes

- **[`mod.rs`](../../../../src/ng/window_coverage/mod.rs)** — `MIN_WINDOW_POSITIONS = 50`, whose
  doc says in its first sentence that the value is soft and D1 sets it; the configuration type's
  doc now says which two of the three constants are production's and inherited and which one
  nobody has yet defended with data. `WindowCoverageConfig::min_window_positions`, with two
  `assert_valid` guards: at least 1, and no more than the window is wide enough to hold — a floor
  above that would silence every window of every sample and say nothing about it.
  `WindowCoverage::absent()` and `is_absent()`, so the absent state is one named thing rather
  than a pair of `f32::NAN`s written at each site; `is_absent` answers "absent" if **either**
  number is missing, because both fields are `pub` and a real depth over a missing GC fraction is
  not a usable measurement.
- **[`accumulator.rs`](../../../../src/ng/window_coverage/accumulator.rs)** — a second counter
  beside the divisor, holding how many *different* coordinates the running sums sit on, which is
  what the floor is judged against; `finalise_centre` decides in one branch and pushes and
  advances its cursor once, on both paths, because `finalise_all` loops on that cursor and an
  exit that forgot it would hang rather than fail; `finish`'s exhaustive destructure names the
  new configuration field and ignores it explicitly, because a floor is not part of the bin
  scheme a consumer reads a cell by.
- **[`production_parity.rs`](../../../../src/ng/window_coverage/production_parity.rs)** — the
  differential runs with the floor at 1, **and asserts that no window it compares is absent**, so
  "the floor is switched off in this comparison" is checked rather than commented.

**The floor is at 1 in every test that predates this step and in the differential** — which is
stronger than switching it off, because a window always holds at least its own centre, so a floor
of 1 sits exactly on the boundary and those tests are themselves what catch a `<` widened to
`<=`. The differential still compares **447,581 window means
and GC fractions bit for bit**, the same count as at A1, which is what says the floor is an
addition rather than a change.

## Tests added

Eleven — five written with the step, six more the review required — and the boundary is tested by
moving the floor across a fixed stream rather than by moving the stream, so that nothing but the
floor can explain the difference:

| test | what it pins |
|---|---|
| `a_window_at_the_floor_speaks_and_one_position_short_of_it_does_not` | five positions inside one window: at a floor of 5 every window reports 10, at a floor of 6 every window is absent |
| `an_absent_window_is_absent_in_both_fields_and_compares_by_bits` | both numbers go absent together, and the pair equals `WindowCoverage::absent()` |
| `a_window_under_the_floor_is_emitted_but_not_folded` | the centre is still emitted; no cell takes it and `windows_folded` stays 0 — and the same stream with the floor off folds both, so the claim is about the floor and not the fixture |
| `an_n_position_does_not_count_toward_the_floor` | the floor counts covered positions, not reference span |
| `a_stream_holding_both_a_silent_window_and_a_speaking_one_is_scored_per_window` | one isolated position and a run of four in one stream: the verdict is per window, not per stream |
| `a_repeated_position_does_not_clear_the_floor_on_its_own` | one base observed five times is one base |
| `the_floor_counts_the_positions_of_its_own_contig` | five per contig meets a floor of 5 and not one of 6, so the count does not leak across the boundary |
| `is_absent_is_true_when_only_one_field_is_missing` | fail-closed on either number |
| `is_absent_accepts_any_nan_payload_but_equality_does_not` | why `is_absent` exists rather than `== absent()` |
| `new_panics_on_a_zero_floor`, `new_panics_on_a_floor_no_window_could_reach` | both ends of the configuration guard |

**Seven mutations, each caught, run in the container on this tree.** Three against the step as
first written, four against it after the review's fixes; each of the last four fails exactly one
test:

| mutation | result |
|---|---|
| the floor comparison `<` → `<=` (a floor of *n* silences a window holding *n*) | 8 tests fail, `a_window_at_the_floor_speaks_and_one_position_short_of_it_does_not` among them |
| an absent window folded anyway, into cell 0 | `a_window_under_the_floor_is_emitted_but_not_folded` |
| only `mean_depth` goes absent, `gc_fraction` stays a real number | `an_absent_window_is_absent_in_both_fields_and_compares_by_bits` |
| the floor judged against records rather than distinct coordinates | `a_repeated_position_does_not_clear_the_floor_on_its_own` |
| absence made sticky — once one window is silenced, so is every later one | `a_stream_holding_both_a_silent_window_and_a_speaking_one_is_scored_per_window` |
| the distinct count left un-reset at a contig change | `the_floor_counts_the_positions_of_its_own_contig` |
| the unreachable-floor assertion removed | `new_panics_on_a_floor_no_window_could_reach` |

The review's own reliability pass ran 13 further mutations against the step as first written: 11
caught, 2 survived, and one of those two was indistinguishable on any value the module can
produce.

## Validation results

In the container (`./scripts/dev.sh`):

- `cargo test --lib --all-features` — **6,312 passed, 0 failed, 15 ignored** (6,306 before the
  review's fixes; 6,301 at A1; 6,275 without this module).
- `cargo test --lib --all-features ng::window_coverage` — 37 passed, 0 failed.
- `cargo clippy --lib --all-features` — no diagnostic in `src/ng/window_coverage/`; the three
  `needless_lifetimes` warnings elsewhere predate this branch.
- `rustfmt --check` clean on all three module files.

**The calling oracle was not re-run and does not need to be**: nothing calls the module, so no
VCF byte can move. The plan runs it either side of A3, which lands alone.

## Tradeoffs and follow-ups

- **The floor's value is not measured yet**, and until D1 measures it nothing here says whether
  50 of 500 is right. What A2 fixes is that the mechanism exists and is tested; what D1 supplies
  is the number.
- **`is_absent` is fail-closed rather than enforced.** Both fields are `pub`, so a consumer can
  still build a half-absent pair; what the predicate guarantees is that such a pair is treated as
  absent rather than scored. The stronger form — a private field and a constructor — is a
  follow-up for plan step C2, when the first builder outside this module exists; today a
  two-`f32` positional constructor would buy a transposition hazard for a lock nothing needs.
- **Three things are for the owner at Checkpoint A**, none of them the implementer's to change:
  spec §3.4's "Every finalised window is folded" now contradicts the code; a
  `windows_under_the_floor` counter on the histogram is recommended but changes the type the
  filter plan reads; and spec §5's 12-byte deque entry is 24 bytes, carried over from A1.
- **A window under the floor still costs what a window costs.** The floor decides what a window
  *says*, not whether it is computed; the two-pointer work is unchanged.
