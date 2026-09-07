# window coverage — C5: the histograms out of the cache

**Date:** 2026-09-07
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone C, step C5
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.4, §3.5, §10
**Review:** [ng_window_coverage_c5_2026-09-07.md](../reviews/ng_window_coverage_c5_2026-09-07.md)
**Branch:** `ng-window-coverage`
**Builds on:** [C2](ng_window_coverage_c2_2026-09-06.md), [C4](ng_window_coverage_c4_2026-09-07.md)

## The answer

**Every sample's coverage histogram now leaves the merge when it returns its sources, and it is
the one a straight walk of the same store makes.** On the six tomato accessions over two 100 kb
intervals:

- **All six samples' histograms are byte-identical to the whole-store recomputation's.**
- **The two modes' are byte-identical to each other**: 6 samples, 6 of them fitted, which the
  mode-equivalence oracle now compares itself.
- The VCF and the per-record windows are unmoved: 2,311 records, sha256 `84ad19c2…`, 13,866 window
  rows of which 13,589 carry a measurement.

**`windows_folded` accounts for every covered position.** On this slice the floor silenced none,
so each sample's folded count *is* its covered-position count — 194,942, 198,710, 193,547,
187,762, 198,681 and 198,600, each equal to the positions the recomputation observed. The review
reproduced the same identity on a different slice where the floor did silence some: 300,943 folded
plus 377 silenced against 301,320 covered.

## Why finishing is the step, not a formality

The histogram is the **yardstick** half of the hidden-duplication filter's input: the per-locus
pairs (step C4) say what a sample's depth was around a locus, and this says what one copy's depth
looks like in that sample, so the filter divides one by the other. Both halves have to come from
one accumulator or they are measurements on different scales — which is what spec §1 exists to
arrange.

`WindowCoverageAccumulator::finish` does two things nothing else does:

- **It closes the centres no cover can.** A window centred at `p` is finalised once the sample's
  stream reaches `p + 250`; the last half-window of a sample's own records has no later position,
  so no cover reaches it however far the look-ahead draws. Measured by the review on three
  accessions over 400 kb: 250, 250 and 151 centres, **about 6 windows in 10,000**. Without it every
  histogram is short by its tail and the byte-identity above fails outright.
- **It fits the depth axis for a sample whose pass ended early** — one that finalised fewer than
  `depth_scale_windows` (10,000) windows never reaches the point where the width is set, and comes
  back with no histogram at all rather than an unfitted one. **That half is carried by unit tests
  only**: every sample of both slices passed 10,000 long before its stream ended, so the case
  arises below about 10 kb of covered ground, or at a sample sparse enough that the caller's own
  low end is what would meet it.

**The tail windows `finish` hands back are dropped.** A window is read at the locus a builder is
building, and by the time this runs every builder has run. They are the same centres the
whole-store recomputation keeps and a run does not — which is why the *windows* comparison shows a
handful of loci with nothing on the run's side while the *histograms* compare equal.

## Changes made

- **[`observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)** —
  `into_sources_and_histograms` finishes every accumulator and hands the histograms back beside the
  sources. **`into_sources` is gone**: it would now be a shorter name that does that work and drops
  the result, and a caller reaching for it would lose the one measurement nothing downstream can
  recompute without a second pass over the reads. Two tests.
- **[`callers.rs`](../../../../src/ng/run/callers.rs)**,
  **[`psp_caller.rs`](../../../../src/ng/run/psp_caller.rs)** — `CalledCohort`, `WrittenCohort` and
  `StoredCohortTallies` each carry the list, one entry per sample, unread until the filter's plan.
  psp mode asserts the pairing against its file count, as it already does for the read counts.
- **[`recorded_windows.rs`](../../../../src/ng/run/cohort_merge/recorded_windows.rs)** —
  `record_the_histograms` writes them to the file `histograms_beside` names, one line a sample;
  `write_the_histogram` is `pub` so the recomputation renders its own side through it. Four tests.
- **[`ng_window_coverage_probe.rs`](../../../../examples/ng_window_coverage_probe.rs)** — the
  comparison, per sample.
- **[`ng_mode_equivalence_oracle.sh`](../../../../scripts/ng_mode_equivalence_oracle.sh)** — the
  two modes' histogram files compared, with three guards.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched.

## Assumptions and deviations, all minor and all recorded

1. **The plan says `into_sources` is "joined by" a form that returns the histograms; it is
   *replaced* by one.** Joined, the shorter name finishes every accumulator and throws the result
   away, and it had no callers the moment the two real ones moved. A caller with no use for the
   histograms writes `.0`, which shows a reader what is being dropped.
2. **The histograms are written as text and compared as text**, through one library function both
   sides render with. Spec §10 asks for byte-identity and the shell oracle can only compare bytes;
   a structural comparison would need a second, weaker instrument for the same claim. What text
   gives up is a readable failure, which is what the two "name the differing bin" fixes below buy
   back.
3. **psp mode carries the list on `StoredCohortTallies`** rather than beside the per-file counts.
   A `StoredSample` is built per file as the psps are paired back with what they gave; the
   histograms come out of the cache in one piece at the end of the pass.

## What the review changed

Two agents; the full report is [beside this one](../reviews/ng_window_coverage_c5_2026-09-07.md).
**5 Major, 8 Minor**, all applied. The measurement itself survived every attack on it — `finish` is
called exactly once per sample, only after both error returns, and the run's histograms equal the
walk's on two independent slices. What needed fixing was everything around it:

- **Both failure reports truncated away the difference they exist to show.** A fitted line is eight
  header fields and 20,050 cells; the probe printed the first 160 characters of each side and the
  oracle diffed the two files cut to 120 columns, so a cell difference reported itself as two
  identical prefixes — or, in the shell, as nothing at all. Both now name what differs: a header
  field by its name, a cell by its GC bin and depth bin.
- **A new function landed inside another function's doc comment**, so one function was documented
  as doing what the other does and the other had no comment. The same defect this project's review
  history has recorded once before.
- **The `.histograms` path was derived three times**, twice as the same four lines in different
  crate targets — in a module whose own doc argues that both halves of a format must live together.
  One `histograms_beside` now, with a test over the odd paths.
- **A fitted line's cell count was never checked against its own scheme.** Both sides render
  through one function, so a `counts` disagreeing with `gc_bins × (depth_bins + 1)` produces a line
  no reader can index — and two such lines compare *equal* when both are wrong the same way, which
  is a pass in both comparisons this file exists for. Prevention, not a live defect: the
  accumulator allocates exactly that and never resizes, but every field is `pub`.

**Two mutations survive the library suite and are worth stating.** A wrong sample index in
`record_the_histograms`, and a probe that compared only the first sample, are both caught only by
running the probe on a real store — and the *mode-equivalence* oracle cannot catch the first,
because both modes render the same wrong index. That is the same standing the row half has had
since C3: the oracle proves the two modes agree, the probe proves what they agree on is right, and
the probe is run by hand.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,378 passed, 0 failed, 15 ignored** (6,373 before this
  step); `cargo check --lib --tests --all-features` clean; the probe's own 7 passed.
- `cargo clippy --lib --all-features --bins --example ng_window_coverage_probe` — 3 warnings, all
  `needless_lifetimes` in `src/ng/run/cohort_merge/`, all predating this branch.
- `rustfmt --check` — every file this step touched at the hunk count it had before it.
- **The standing oracle**: 2,311 records, sha256 `84ad19c2…0590d` on both routes; the two modes'
  window coverage identical over 13,866 rows, 13,589 of them a measurement; **and their coverage
  histograms identical, 6 samples, 6 of them fitted**.
- **The whole-store recomputation**: all six samples' histograms agree with the run's byte for
  byte; 13,589 of 13,866 sample-loci agree on their window, 0 disagree.

## Tradeoffs and follow-ups

- **80.2 kB a sample, now held past the merge** — 80 MB at a thousand samples, 240 MB at three
  thousand, which is spec §3.4's own figure and what plan step D3 prices against
  [`run_streaming.md`](../../ng/spec/run_streaming.md) §7.2's budget. Finishing *lowers* live heap
  at the moment it runs: the accumulator held the same counts vector plus its position buffer.
- **Nothing reads the histograms.** The hidden-duplication filter is their consumer and is built on
  its own branch.
- **The early-pass half of `finish` is untested on real data**, above.
