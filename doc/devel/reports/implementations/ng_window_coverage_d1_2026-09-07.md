# window coverage — D1: how many positions a window holds, and what the floor should be

**Date:** 2026-09-07
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone D, step D1
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.3, §9
**Reviews:** [correctness](../reviews/ng_window_coverage_d1_correctness_2026-09-07.md),
[naming and structure](../reviews/ng_window_coverage_d1_design_2026-09-07.md)
**Branch:** `ng-window-coverage`
**Builds on:** [C3](ng_window_coverage_c3_2026-09-06.md), [C5](ng_window_coverage_c5_2026-09-07.md)

## The answer

**50 of 500 stands, and the floor turns out to be a rule about how long the run's analysed
intervals are rather than about how deep the sample is.** Measured over five stores spanning
5.1 to 301 reads compared with the reference a position and analysed intervals from 122 bases to
100 kb. Windows the floor would silence, per 10,000:

| store | intervals asked for | reads a position | windows | under 50 | under 100 | under 200 |
|---|---|---|---|---|---|---|
| tomato, 6 accessions | 2 × 100 kb | 14.4 | 1,172,242 | **0.00** | 0.08 | 8.10 |
| tomato, 1 accession | 80 × 100 kb | 10.3 | 7,667,288 | **1.13** | 1.69 | 15.39 |
| HG002, GIAB high-confidence | 1,000 × 5.1 kb | 301.4 | 5,046,746 | **0.06** | 6.68 | 16.18 |
| HG002, tandem-repeat tiers | 50,000 × 122 bp | 30.3 | 5,910,300 | **360.13** | 4,393.07 | 5,848.64 |
| HG002, tandem-repeat tiers | 50,000 × 122 bp | 5.1 | 5,833,898 | **381.69** | 4,432.50 | 5,913.10 |

*Reads a position is the mean `reads_compared_with_reference` at a record covering one base — the
number the depth rule reads at 99 records in 100.*

**The last two rows are the same ground at a sixth of the depth, and the share silenced at 50
barely moves: 360 windows in every 10,000 against 382.** Between the third row and the fourth the
intervals shorten 42-fold and that share goes from 28 windows of 5,046,746 to 212,850 of
5,910,300 — **6,500 times as many**. So depth is not what decides how often the floor bites;
interval length is.

**Why, and it is arithmetic rather than a property of these data.** A window is the covered
positions within 250 bases either side of a centre. An interval shorter than 250 bases therefore
lies inside *every one* of its own windows, so on a run whose intervals are 122 bases every window
holds about 122 positions and the floor is asking how short an interval is too short. Checked
against the request rather than the answer: **2.9 in 100 of the tandem-repeat tiers' requested
bases lie in intervals under 50 bases and 42.8 in 100 in intervals under 100** (from the BED
alone, `awk` over `HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed`), against the **3.6 and
43.9 in 100** of windows the floor silences there.

Two things could have made that check weaker than it looks, and neither does. **It could have
begged the question** — it does not, because the two BED figures are computed from the region file
alone and never touch the store, the walk or the accumulator, so a shared misunderstanding of what
a window spans cannot cancel out. And **it assumes a window never reaches into a neighbouring
interval within 250 bases**, which would lift a count above its own interval's length; computing
each requested base's whole 501-base neighbourhood against the whole BED moves the prediction only
to 2.89 and 42.48 in 100, so the assumption costs at most a third of a point. The remaining excess
over the measured 3.6 and 43.9 is the positions no read reached: 5,910,300 of the 6,089,411
requested bases were covered at 30 reads a position, so 2.9 in 100 were not.

## What settles the number

**Against raising it to 100.** On short-interval ground that would silence 44 windows in 100, and
nothing measured here says such a window is wrong: at 122 positions and 5.1 reads a position a
window is still a mean over about 630 reads, where spec §4's reason for windowing at all is that a
*single* position at 3 reads cannot separate one copy from two. Whether a thin window's copy
number is wrong is the hidden-duplication filter's own question, measured on its own branch —
**this step measures the cost of the floor and not its benefit**, and says so rather than reading
the cost as a verdict.

**Against lowering it to 10.** The one-accession tomato store holds **261 windows with fewer than
10 covered positions** and 867 with fewer than 50, on ground made of 100 kb intervals. Those are
the isolated covered positions the floor exists for, and a floor of 10 would let 606 of them
speak.

**On the ordinary case it costs nothing.** Across the three stores whose intervals are 5 kb or
longer, the worst is 1.13 in 10,000 — **1 window in 8,850**.

## What the six-accession slice could not say, and what was built to say it

**On the slice this branch has used at every step, the floor silences nothing** — 0 windows of
1,172,242, in all six samples. A default confirmed there would have been confirmed by data that
never tested it. Three stores were built for this step so that it was:

- **tomato, one accession over all 80 benchmark intervals** (`tmp/b2_wide_psp`, from plan step
  B2): the same interval length, forty times the ground, and the first store where the shipped
  floor silences anything at all;
- **HG002 over the `human_genome_bottle` benchmark's 1,000 intervals of about 5 kb**, from
  `benchmarks/human_genome_bottle/crams/HG002_reads_selected_1000_rg.cram` — the top of the depth
  range this caller commits to, at 301 reads a position;
- **HG002 over the 50,000 tandem-repeat Tier intervals at 30× and at 5×**, from
  `benchmarks/ssr_hg002/bam/{30x,5x}/`, which are the same ground at two depths and the only pair
  here that isolates depth from interval length.

**No whole-genome store of either species exists on this machine and none can be built from these
files** — every alignment file here is cut to its benchmark's intervals, so a run over whole
contigs would see the same islands with empty stretches between them. What that costs is the one
case this measurement does not have: a sample whose *coverage* is patchy over continuous analysed
ground, as opposed to one whose analysed ground is itself in pieces. The 5× row is the nearest
thing to it, and it moved the share silenced from 360 windows in every 10,000 to 382.

## How the distribution is measured

**Read off the shipped accumulator, at one instance per candidate floor, rather than recomputed.**
The count of covered positions in a window is private to `WindowCoverageAccumulator` and no method
hands one back; what it does report is how many windows a given floor silenced, and *that* number
at floor `F` is exactly the number of windows holding fewer than `F` positions. So
`ng_window_coverage_probe --covered-positions-per-window` runs 25 accumulators — one per candidate
floor, each called an *arm* — over the positions the whole-store walk is already feeding one,
identical in every setting but the floor. Their answers, read in floor order, are the cumulative
distribution.

Nothing in the probe reimplements the window. That matters because the alternative — a second
sweep over the same positions counting them by hand — would be one copy of the rule checked
against another, and a shared misunderstanding of what a window spans would cancel out.

**Six things are asserted before any of it is printed**, because a row of zeroes at every floor
is what both a correct measurement on contiguous ground and a completely broken one look like:

- every arm finalised the same number of windows as the walk beside it — the floor decides whether
  a window speaks, never whether it exists, so a mismatch means the two were not fed the same
  positions;
- **the arm at the floor the run actually ships silenced exactly what the walk itself counted
  absent** — the one place the sweep and the recomputation answer the same question. It has teeth
  on four of the five stores, where that number is 28, 867, 212,850 and 222,673 rather than zero;
  on the six-accession slice it compares 0 with 0, and the printed `sweep-tied-to-the-walk` row
  says which of the two a store gave;
- **the arm at a floor of 1 silenced nothing**, since a centre always lies in its own window —
  the only one of the six that names a number rather than comparing two quantities that can both
  be zero;
- arms exist at both floors those two checks are written against, since each is written as "if
  this is the row at floor `F`" and a hand-edited floor list that lost `F` would run neither;
- the walk finalised at least one window;
- a higher floor never silenced fewer windows than a lower one.

The shipped floor is folded into the sweep's floor list at runtime rather than written into it, so
that the second check cannot quietly stop applying if the default moves.

**What none of the six can see is the arms' own configuration**, and the correctness review
constructed the defect that exploits it: build the arms with a 1,000-base window instead of the
run's 500. A window is finalised per covered position *whatever the configuration says*, so the
window-count check is blind to it; on the six-accession slice the shipped-floor check compares 0
with 0; monotonicity still holds — and the run prints a full distribution of a window the caller
never builds. A *narrowing* drift is caught, but by the library rather than by the sweep
(`WindowCoverageConfig::assert_valid` refuses a floor wider than the window). What pins the width
is the unit tests, which fix the counts a 500-base window produces over 600 consecutive positions;
the practical guard is the one this step honoured — walk at least one store whose own windows the
floor silences something in.

## Assumptions and deviations, all minor and all recorded

1. **The plan names "the tomato slice and HG002"; five stores were walked, not two.** The slice
   silences nothing, so it settles nothing on its own, and "HG002" turned out to name two
   different benchmarks in this tree with intervals 42-fold apart — which is the axis the answer
   depends on. Both are reported with their interval length beside them.
2. **The floor list stops at 501, not higher.** `WindowCoverageConfig::assert_valid` refuses a floor
   above the widest a window can be (a 500-base window spans 501 positions), so 501 is the
   highest arm and it counts the windows whose ground is *not* completely covered.
3. **The sweep is behind a flag.** It costs one more pass of the window arithmetic per arm, and
   the probe is run at every step of this plan for the C3/C4/C5 comparisons, which do not need it.

## Changes made

- **[`examples/ng_window_coverage_probe.rs`](../../../../examples/ng_window_coverage_probe.rs)** —
  `FloorSweep` and its arms, the `WindowsUnderEachFloor` answer a closed sweep hands back, the
  `--covered-positions-per-window` flag, the per-store and total rows, and
  `the_configuration_a_run_uses` extracted so that the walk's accumulator and every arm cannot
  drift apart. Eight tests.
- **[`src/ng/window_coverage/mod.rs`](../../../../src/ng/window_coverage/mod.rs)** —
  `MIN_WINDOW_POSITIONS` keeps its value of 50 and stops being marked soft; its doc carries the
  table above and the reason.
- **[`doc/devel/ng/spec/window_coverage.md`](../../ng/spec/window_coverage.md)** — §3.3's "until
  measured" and §9's OPEN both closed with the measurement.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched, and no
library code changed behaviour: the only edit under `src/` is a constant's doc comment.

## Tests added

Eight, in the probe.

| test | what it pins |
|---|---|
| `an_arms_silenced_count_is_the_windows_holding_fewer_positions_than_its_floor` | that an arm's silenced count is the windows holding fewer positions than its floor, on ten positions each inside every one of the ten windows, where the answer is countable by hand |
| `the_sweep_traces_a_distribution_whose_windows_hold_different_counts` | 600 consecutive positions, where the counts run 251 to 501 by arithmetic: three floors, below the distribution, through it, and at its top. **It is also the only test that reaches the mid-stream drain**, since ten positions never pass a centre's right edge |
| `two_sweeps_sum_floor_by_floor` | two identical stores summed double every count (1,200 windows against 600), and the total's floors are the same floors in the same order |
| `a_sweep_whose_windows_differ_from_the_walks_is_refused` | a sweep fed different positions from the walk beside it stops rather than reporting a distribution of a stream nothing else saw |
| `a_sweep_disagreeing_with_the_walk_at_the_shipped_floor_is_refused` | the one check that ties the sweep to the measurement it is about |
| `an_arm_at_a_floor_of_one_silencing_a_window_is_refused` | the one check that names a number rather than comparing two that can both be zero |
| `a_sweep_that_lost_the_floor_a_check_is_written_against_is_refused` | a floor list that no longer holds a floor a check is written against fails, rather than running no check |
| `a_sweep_over_an_empty_store_is_refused` | a row of zeroes at every floor from an empty store is refused, not printed |

## What else the five walks said

The same walks answer plan step B2's question, because the probe asks it on every store it opens.
**The head count equals the evidence's own sum at every one of 24,793,279 records covering one
base**, across all five stores — 0 disagreements, 0 records carrying a witness that stops short,
and no one-base record that was a repeat tract. B2 measured this on 8,784,182 records of tomato
only; it now also holds on human data at 5.1, 30.3 and 301.4 reads a position. **The cheap branch
decides 99.2 and 99.4 positions in 100 on the two tomato stores, 98.9 in 100 on the GIAB
high-confidence store, and 93.8 in 100 on the tandem-repeat one**, where more of the ground is
repeat tract.

## What the reviews changed

Two agents, each detached at the step in its own worktree; the full reports are beside this one,
[one on correctness, errors and every number this prose claims](../reviews/ng_window_coverage_d1_correctness_2026-09-07.md)
and [one on naming, idiom, module structure and refactor safety](../reviews/ng_window_coverage_d1_design_2026-09-07.md).
**6 Major, 14 Minor.** All applied except four the correctness review recorded rather than fixed,
listed at the end of this section. **Neither agent could break the measurement**: the correctness
one read the accumulator against the claim case by case — repeated positions, `N` bases, contig
changes, the `finish` tail, contig-edge truncation, the depth-bin-width fit — and re-derived every
hand-computed test expectation and every figure in the table above from the raw output and from
the region files. Its mutation table lists seven deliberate defects, six of them caught.

The four that were defects rather than clarity:

- **The four checks reduced to one, and that one is silent on the slice.** Both agents arrived at
  the same missing check independently — the arm at a floor of 1 must silence nothing — and the
  correctness agent then built the defect that survives even so: a bank configured with a
  1,000-base window. Two further checks and a printed row now bound it, and what none of them can
  see is stated above rather than left for the next reader to discover.
- **Two numbers in this report were wrong, both flattering the measurement.** "6,000-fold" came
  from dividing two *rounded* printed shares and is 6,491; the shipped-floor check has teeth on
  four stores, not three, because the GIAB store's 28 is also not zero. Both corrected, with the
  counts beside them.
- **The counting rule was written out twice** — once in the mid-stream drain, once for the tail —
  and the mutation that inverts it in the first was caught by only one of the six tests, because
  ten positions never pass a centre's right edge and so never reach that drain at all. One rule
  now, called from both.
- **`FloorSweep::empty()` forced four other things**: an `Option` around every arm's accumulator,
  two `expect`s, a destructure arm explaining why an accumulator is not summed, and a comment
  describing a clone that did not happen. Splitting the live sweep from its answer removes all
  five — the answer is a `WindowsUnderEachFloor` a closed sweep hands back.

And the prose: the module's own "one walk, several measurements" section described a file with two
measurements after this step added a third; `WindowCoverageConfig::min_window_positions` still
called the value provisional 45 lines from the constant this step un-provisioned; and "moves by
6 in 100" carried a *relative change in a rate* in three documents where every other "in 100" is a
*share of windows* — a reader carrying the established unit across would have read the depth effect
as about 600 times its size.

**Four findings recorded rather than fixed**, all from the correctness review:

- **The surviving mutation.** Draining `pop_ready` once per observed position instead of looping
  survives all fifteen tests and all six checks, because `finish` hands back the undrained queue
  and the counts come out identical. What it costs is memory — about 142 MB an arm over 5.9 M
  windows, times 25 arms — and nothing here measures memory, so nothing here can catch it.
- **"about 630 reads" is 626** — 122 positions × 5.13 reads a position. Left as written, since it
  says "about" and 626 rounds to 630 at two figures.
- **"Those are the isolated covered positions the floor exists for" is inferred, not measured.**
  The probe reports counts, not where those 261 windows sit or what surrounds them. The inference
  is safe — a window holding under 10 of a possible 501 positions *is* near-empty ground — but the
  word does work this measurement does not supply.
- **The depth column is a mean over one-base records only**, which its footnote says exactly and
  spec §9's "spanning 5.1 to 301 reads compared with the reference a position" drops. On the
  tandem-repeat store the wide tract records cover 359,907 bases — 6.1 in 100 of the positions —
  whose depth is not in that mean.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — 6,378 passed, 0 failed, 15 ignored, unchanged: this step adds
  no library test and changes no library code but two doc comments.
- `cargo test --all-features --example ng_window_coverage_probe` — **15 passed**, 0 failed
  (7 before this step).
- `cargo check --lib --tests --all-features` — clean.
- `cargo clippy --lib --all-features --bins --example ng_window_coverage_probe` — 3 warnings, all
  `needless_lifetimes` in `src/ng/run/cohort_merge/`, all predating this branch.
- `rustfmt --check` — every file this step touched at the hunk count it had before it.
- **The standing oracle**: 2,311 records, sha256 `84ad19c2…0590d` on both routes; the two modes'
  windows identical over 13,866 rows, 13,589 of them a measurement; their histograms identical,
  6 samples, 6 fitted.

## Tradeoffs and follow-ups

- **The cost of the floor is measured; its benefit is not.** How wrong a thin window's copy number
  is can only be asked once something divides one by the other, which is
  [`hidden_paralog_filter.md`](../../ng/spec/hidden_paralog_filter.md)'s branch. If that
  validation finds thin windows biased, this default is the knob it should reach for, and the
  table above says what each setting costs.
- **The case this measurement does not have** is a sample whose coverage is patchy over
  continuous analysed ground. Nothing in this tree can build one.
- **The sweep runs by hand.** Like the rest of the probe's real-data checks, nothing in CI walks a
  store, so the six unit tests are what guard the arithmetic between runs.
