# ng cohort merge — E1: ordered release

*Implementation report, 2026-08-18. Step E1 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §6.3 and [arch](../../ng/arch/cohort_merge.md) §4, §5.*

> **This is the first draft, and [the review](../reviews/ng_cohort_merge_e1_2026-08-18.md)
> changed three things in it.** A gap at the *tail* of a run finished `Ok` — the step's own
> contract, missed; the reorder map was written inside an `assert!` condition; and
> `MissingRegionResults` became the three-variant `RunEndedShort`. What landed is in
> [the fix report](ng_cohort_merge_e1_fixes_2026-08-18.md); §4 below describes the draft.

## 1. Plan

Add the organiser's reorder buffer: take one outcome per building region, keyed by the index
the run handed the region out under, and let the loci out along an unbroken run of indexes.
**Resolve no overlaps** — that is the next step — and **hold no cache** — that is the step
after.

## 2. What it is

Three items in
[`organise.rs`](../../../../src/ng/run/cohort_merge/organise.rs), beside the observation cache:

- **`RegionIndex(u64)`** — where a building region falls in the order the run hands them out,
  0 for the run's first, counting on across every analysed region. It is the run's numbering
  and not a coordinate, because the drain has to tell *"region 7 has not arrived yet"* from
  *"region 7 found nothing"*, which two genome positions cannot say.
- **`Organiser`** — `submit`, `drain_ready`, `failed_locus_count`, `is_finished`, `finish`.
  A `BTreeMap<RegionIndex, RegionOutcome>` of what arrived early, a `next_expected` cursor, a
  `VecDeque` of released loci, and the failed-locus total. This is production's reorder map
  carried whole: `VcfWriter` holds the same `BTreeMap` and drains it on `next_expected`
  (`var_calling/vcf_writer.rs:168-176`, drained at `:246`).
- **`MissingRegionResults`** — the refusal at the end of a run.

## 3. The three things E1 had to get right

**Exactly one outcome per region, empty ones included** (spec §6.3). A region that found
nothing must still submit, or every region behind it is held for ever — the drain advances
only along an unbroken index run. Two tests state the two halves:
`an_empty_region_still_lets_the_region_behind_it_out`, and
`a_region_that_never_submits_holds_back_every_region_after_it`.

**A gap is an error, not a truncation.** `finish` consumes the organiser and returns the
failed-locus total, or refuses and names what would have been lost. This is production's
`WriterError::MissingChunks` (`var_calling/vcf_writer.rs:152-158`), which arch §5 will fold
into `RunError::MissingRegionResults`.

**The count is of regions released, not of regions submitted** (spec §3.3). A region held
behind a gap contributes nothing to the total until it is let out, which is what lets E2 drop
a failed locus that an earlier one displaced without the total having counted it already.
Pinned by `the_failed_locus_count_ignores_a_region_still_held_behind_a_gap`.

## 4. Assumptions and deviations, both small

**Superseded by the review — see the note at the top.** **`MissingRegionResults` carries two
counts where arch §5 wrote one.** Arch's variant is
`{ count: usize }` — regions never emitted — because production emits from inside its own
`submit` and has no second step to forget. Here taking the released loci is the caller's own
call, so a run can end short a second way: loci released and never drained. Both truncate the
output, so the struct carries `regions_never_released` and `loci_never_drained`, and the
message names both. `finish_refuses_a_run_with_loci_released_and_never_drained` is the case
arch's single count could not express.

**`drain_ready` takes one locus at a time from the front** rather than emptying the buffer
into the returned iterator. Arch's signature (`impl Iterator<Item = CohortObservation> + '_`)
is unchanged; the difference is what a caller who stops halfway loses, which is nothing —
`a_drain_stopped_halfway_leaves_the_rest_for_the_next_call` pins it, and it is the one thing
`self.released.drain(..)` would get wrong.

**The release runs inside `submit`, not inside `drain_ready`.** A region becomes releasable
the moment the index in front of it arrives, which has nothing to do with when the caller next
asks for loci; production does the same. It is also what makes `failed_locus_count` meaningful
before anything has been drained.

**Two panics rather than errors.** A region submitted twice, and a region submitted after it
was released, are bugs in whoever hands the regions out — not facts about the data — so they
stay panics when the observations arrive from a psp file. That is the opposite of the three
assertions in `build.rs`, which are about the producer's data and are recorded as owed
`RunError` variants.

**Not in this step, and named where they will land:** the organiser does not hold the
observation cache (arch §4 gives it one) and resolves no overlaps. The cache is drawn forward
by whoever hands regions out, which is the parallel arrangement's shape to settle (E3), and
the overlap rule is E2's.

## 5. Tests

Sixteen, all in `organise.rs`'s test module, which went from 32 tests to 48. The whole
`ng::run::cohort_merge` module went from 168 to 184 (`cargo test --lib ng::run::cohort_merge`:
`184 passed`). **The review added five more, and the counts after it are in
[the fix report](ng_cohort_merge_e1_fixes_2026-08-18.md).**

| test | what it pins |
|---|---|
| `a_region_arriving_in_its_turn_releases_at_once` | the plain case |
| `a_region_that_arrives_early_waits_for_the_one_before_it` | the buffer's whole job |
| `a_convoy_held_behind_a_gap_releases_in_region_order_when_the_gap_closes` | four regions arriving 3,1,2,0 come out 0,1,2,3 |
| `an_empty_region_still_lets_the_region_behind_it_out` | spec §6.3's empty-result rule |
| `a_region_that_never_submits_holds_back_every_region_after_it` | its other half, and the refusal |
| `the_loci_of_one_region_come_out_in_the_order_the_builder_gave_them` | within-region order |
| `the_failed_locus_count_sums_every_released_region` | spec §3.3's total |
| `the_failed_locus_count_ignores_a_region_still_held_behind_a_gap` | the count is of released regions |
| `a_region_is_released_on_arrival_rather_than_at_the_next_drain` | the release point, without draining |
| `a_drain_stopped_halfway_leaves_the_rest_for_the_next_call` | the front-at-a-time drain |
| `an_organiser_is_finished_once_every_region_is_released_and_drained` | `is_finished` covers both buffers |
| `finish_refuses_a_run_with_loci_released_and_never_drained` | the second count |
| `finish_returns_the_failed_locus_total_of_a_run_that_ended_clean` | the healthy end |
| `an_organiser_that_was_given_nothing_finishes_clean` | a run over no regions |
| `a_region_submitted_twice_is_refused` | one outcome per region |
| `a_region_submitted_after_it_was_released_is_refused` | the cursor is respected |

**Nine mutations, nine killed** (`tmp/mutate_e1.sh`, `tmp/mutate_e1b.sh`; each restores a
pristine copy after every run, so it is safe in the main worktree). The mutation, then the
tests that failed:

| mutation | killed by |
|---|---|
| release moved from `submit` to `drain_ready` | 13 tests, `a_region_is_released_on_arrival_rather_than_at_the_next_drain` among them |
| release ignores the order (`pop_first` instead of `remove(&next_expected)`) | 7, incl. `a_region_that_arrives_early_waits_for_the_one_before_it` |
| failed spans not counted | 4, incl. both count tests |
| `drain_ready` uses `drain(..)` | 1 — `a_drain_stopped_halfway_leaves_the_rest_for_the_next_call` |
| `is_finished` ignores undrained loci | 2, incl. `finish_refuses_a_run_with_loci_released_and_never_drained` |
| `finish` never refuses | 2 |
| a second outcome per region accepted | 1 — `a_region_submitted_twice_is_refused` |
| a released index can be resubmitted | 1 — `a_region_submitted_after_it_was_released_is_refused` |
| `next_expected` never advances | 9 |

Two mutations are killed by exactly one test each, which is the point of those two tests
existing; the rest are killed broadly, as a reorder buffer's mutations should be.

## 6. Validation

Run in the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean. (`--all-targets` is red on this
  branch and was before milestone A: 49 pre-existing errors in `examples/`, `benches/` and
  other modules' test code.)
- `cargo test --lib ng::run::cohort_merge` — `184 passed; 0 failed`.
- `cargo test --lib` — `3807 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out; finished in 701.00s`.

## 7. Follow-ups

- **The organiser does not own the cache yet**, so the open item D2 recorded stands: a failed
  merge leaves the cache advanced, and making that unrepresentable waits on E3 deciding how
  builders and the organiser share it.
- **`MissingRegionResults` is a struct of this module's own**, to be folded into `RunError`
  with arch §5's other variants when the run's error type lands. Its two counts go with it.
