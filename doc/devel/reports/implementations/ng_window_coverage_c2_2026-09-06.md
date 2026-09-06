# window coverage — C2: the accumulator in the merge's per-sample window

**Date:** 2026-09-06
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone C, step C2
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.3, §3.5
**Review:** [ng_window_coverage_c2_2026-09-06.md](../reviews/ng_window_coverage_c2_2026-09-06.md)
**Branch:** `ng-window-coverage`
**Builds on:** [B1](ng_window_coverage_b1_2026-09-06.md), [C1](ng_window_coverage_c1_2026-09-06.md)

## The answer

**Every sample now carries a window-coverage accumulator inside the merge's cache, is fed every
record it holds exactly once, and a builder can read the window at any position.** On the tomato
slice the measurement runs end to end through both modes — the accumulator sees every record of
every sample, in coordinate order, against the bases the cover read — and **no VCF byte moved**:
2,311 records and sha256 `84ad19c22dd14de583cd85805dcd2e5169e799d7a63691c979b7fa43d400590d` on
both routes, as before A1.

Nothing reads the windows yet beyond this step's own tests; the locus carries them at C4.

## The thing implementation found, and it changed C1's fetch

**A cover crosses contigs, and a buffer is one contig's.** Drawing stops at the first record past
the reach, and a reach on a later contig is past every position of an earlier one — so the cover
that first reaches contig *n* also draws whatever the sample still had on contig *n − 1*. C1 read
one contig's ground per cover, so those records had no bases, and the observation pass met them at
once: seven existing tests failed the moment the accumulator was wired in, every one of them a
fixture that crosses a contig.

They are real covered positions of that sample and their bases are on the contig they sit on, so
**the cover now reads each contig it added records on, in ascending order.** In practice that is
one fetch per cover and two at a contig boundary.

**This paragraph said "ending on the region's own", and the correctness review measured that
false**: ascending puts the region's contig last only while no contig sorts after it, and the
record a sample is drawn one *past* the reach is routinely the first of the next contig. C3's own
commit adds the second fetch that puts the region's ground back in the buffer.

This is a correction to code C1 committed, made here rather than amended into it, because it was
the accumulator that showed the gap.

## Assumptions and deviations, all minor and all recorded

1. **The observation happens after the cover's fetch, not inside `draw_to`.** The plan's C2 text
   says `draw_to` observes each drawn record; it cannot, because the bases are not in hand until
   the drawing has stopped (C1's assumption 1). Every record is still observed exactly once, in
   coordinate order, against bases this cover read — which is what the sentence is for.
2. **A per-sample cursor is what makes it exactly once.** A record is held across every cover it
   reaches into, and each of those covers reads ground containing it. The cursor is the start of
   the last record the accumulator saw; a sample's records are disjoint and ascending, so a start
   is enough to say what has been seen.
3. **The finalised windows are a `Vec` drained from the front, not a `VecDeque`.** Spec §3.3 says
   deque; a consumer is handed a contiguous slice, and a deque's two halves are not one — the same
   reason `held_observations` is a `Vec`.
4. **Eviction drops windows by their own coordinate**, not by the records'. A window is centred on
   a position and the record that position came from may already be gone.
5. **B1's rule gained a second entry point.** The cache keeps a summary and its evidence side by
   side rather than the `Drawn` it drew, so `for_each_reported_depth_of` takes the evidence
   borrowed (`EvidenceForOneRecord`) and the `Drawn` form calls it after checking the summary is
   the draw's own. One rule, two ways in, and the tests that were written against the `Drawn` form
   are unchanged.
6. **A position with no base panics.** The ground a cover reads is computed from the very records
   being observed, so a missing base is a defect in that computation and not an absent value.

## Changes made

- **[`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)** —
  `SampleWindow` gains the accumulator, the finalised windows and the cursor; `over` spells the
  configuration, which is the one site that does; `read_the_ground_and_measure_coverage_over_it`
  walks the contigs a cover added records on, fetching each and observing the records there;
  `measure_coverage_on` is the per-sample half; both evictors drop the windows behind the
  evicted position; `WindowedCohort` gains the third view and `window_at`.
- **[`depth.rs`](../../../../src/ng/window_coverage/depth.rs)** — `for_each_reported_depth`,
  the rule over borrowed evidence, and `EvidenceForOneRecord`, the borrowed counterpart of
  `Drawn`.
- **[`build.rs`](../../../../src/ng/run/cohort_merge/build.rs)** — the two `WindowedCohort`
  literals name the new field.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched.

## Tests added

Five, in the cache's own module. The fixture they share is 600 one-base records at three reads
each, which is what makes a window finalise at all: a centre closes only once the stream has
passed 250 bases beyond it, so a fixture shorter than that emits nothing.

| test | what it pins |
|---|---|
| `a_builder_reads_the_window_a_sample_has_at_a_position` | the measurement end to end: the window centred on 300 is a mean depth of **exactly 3** over the covered positions 50 to 550, **251 of whose 501 bases are `G` or `C`** — both numbers predicted from the fixture, not read off the run |
| `a_record_held_across_two_covers_is_observed_by_one_of_them` | every finalised centre is distinct and ascending, and covering the same ground in two goes gives the **same windows** as covering it in one |
| `a_cover_that_crosses_a_contig_observes_the_records_left_on_the_one_it_leaves` | the discovery above: a cover that moves to contig 1 draws contig 0's leftovers, and both contigs' records reach the accumulator |
| `eviction_drops_the_windows_behind_it_and_keeps_the_rest` | the deque does not grow with the contig, and nothing before the evicted position survives |
| `a_sample_with_no_record_at_a_position_has_no_window_there` | a centre the stream has not passed by half a window is not finalised, and reading it gives `None` rather than a fabricated pair |

**This paragraph claimed the seven existing tests that failed when the accumulator was first wired
in were the strongest evidence here, and the correctness review measured that false.** Reverting
the multi-contig read on the committed tree fails **no** pre-existing test — only
`a_cover_that_crosses_a_contig_observes_the_records_left_on_the_one_it_leaves`, which this step
wrote. The seven were real of an intermediate tree that had the accumulator without the contig
`break` in the observation pass; with that break a record on another contig is silently skipped
rather than panicking on a missing base. One test guards the multi-contig read, and it is this
step's own.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,358 passed, 0 failed, 15 ignored** (6,353 before this
  step). The five above are the whole increase.
- `cargo test --lib --all-features cohort_merge::observation_cache` — 47 passed, 0 failed.
- `cargo clippy --lib --all-features` — 3 warnings, all `needless_lifetimes` in
  `src/ng/run/cohort_merge/`, all predating this branch.
- `rustfmt --check` — every touched file back at exactly the hunk count it had before this step
  (`observation_cache.rs` **4**, one fewer than before this step because the evictor extraction
  rewrote a hunk out of existence; the rest unchanged); `depth.rs` clean.
- **The standing oracle after this step**: 2,311 records, sha256 `84ad19c2…0590d` on both routes.
  **This is the first step where the measurement runs on real data**, so it is the first oracle run
  that could have moved and did not.

## What the review changed

**A design review ran on this step and its fixes are in this commit. The correctness review was
still running when the session ended and wrote nothing**; it was re-run against this commit a
session later, found four Major defects and six wrong figures in this report, and is written up
at [ng_window_coverage_c2_2026-09-06.md](../reviews/ng_window_coverage_c2_2026-09-06.md). Its
fixes are in their own commit after C3, and the corrections it forced are marked in place above
rather than edited away.

- **`ground_this_cover_holds_on` scanned every sample's whole held window, once per contig, per
  cover.** Held summaries ascend in `(contig, position)`, so one contig's are a contiguous run
  and its ends are two binary searches. Against `cover`'s own measured shape — 3,000 samples,
  200 held each — the scan was **600,000 summary reads a contig** where the ends are 6,000, and
  a cover at a contig boundary reads two.
- **The measurement's three fields left `SampleWindow` for a `WindowCoverageInProgress` of their
  own.** They shared nothing with the other eight, and keeping them apart is what collapses the
  cursor's partition, which existed twice verbatim, and what makes C5's move-out one field.
- **The two evictors are one method now.** They had grown to share five statements, and the
  standing oracle for this module is that the serial and parallel paths do not drift.
- **The base lookup existed twice.** C2 wrote the offset arithmetic a second time inside the
  per-position closure, and that copy dropped the contig check the accessor makes. One free
  function, `base_in`, is what both read — the same argument the depth rule makes one level
  down about two derivations of "depth".
- **`measure_coverage_on` took the contig twice**, once directly and once inside the
  buffer's origin, and had to be trusted to get both the same. It takes the origin alone now.
- **`Drawn::evidence` became `impl From<&Drawn> for EvidenceForOneRecord`**, in the module that
  owns the borrowed type: written the other way, the merge named a `window_coverage` type in a
  method it does not itself call.
- **B1's `Drawn`-shaped entry point is deleted.** C2 supplied a caller for the *borrowed* form,
  and the compiler then said the other had none outside its own tests. One rule, one way in;
  its contract prose moved with it, and the one test whose subject was the deleted check now
  pins the `From` mapping instead.
- **`window_at` became `window_coverage_at`**, and `WindowedCohort::coverage` became
  `finalised_windows`, because "coverage" already means read depth in this file and "window" a
  stretch of held records.

**One finding is left open on purpose, and it is C4's to settle.** `window_coverage_at` answers
`None` both for "this sample has no window here" — an ordinary answer — and for "this window
carries no measurement at all", which is a gap. Three live construction sites pass no
measurement, and one of them is `merge_cohort_serially`, which is **the oracle C4's plan names
as its green criterion**: comparing absent against absent would pass. C4 must either make the
two answers different at the type level or have its comparison refuse an unmeasured window
rather than compare it.

## Tradeoffs and follow-ups

- **The accumulator's memory is now per sample and unbounded until `finish`.** The histogram is
  80.2 kB a sample once its depth axis is fitted, and until 10,000 windows have been finalised the
  buffered pairs are held instead. Plan step D3 measures the milestone's memory; this is the step
  that adds most of it.
- **The finalised windows are held until eviction**, one 24-byte entry per covered position per
  sample. The organiser's eviction is what bounds it, exactly as it bounds the held records.
- **A cover at a contig boundary makes two fetches**, one per contig. The alternative — a buffer
  spanning contigs — would need a second origin and a second bound for one cover in a genome's
  worth of them.
- **Nothing reads `window_at` yet.** C4 puts the pair on the cohort locus and beside the record at
  the sink; C5 takes the histograms out of the cache.
