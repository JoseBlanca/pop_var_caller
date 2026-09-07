# window coverage — D1 review: correctness, reliability, and every number the step claims

**Date:** 2026-09-07
**Reviewing:** commit `444672b0` ("wip: D1 for review") against `cdfd7821`
**Step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone D, step D1
**Report under review:** [ng_window_coverage_d1_2026-09-07.md](../implementations/ng_window_coverage_d1_2026-09-07.md)
**Remit:** reliability, correctness, error handling, and every self-claimed number. Naming and
structure are another reviewer's.

## The short version

**The measurement is sound and the ruling it supports stands.** The claim the whole step rests on
— that an arm's count of silenced windows at floor `F` is exactly the number of windows holding
fewer than `F` covered positions — is true of this accumulator in every case I could construct,
and I re-derived the three hand-computed test expectations and every figure in the report's table
from the raw probe output and from the BED files. Two numbers in the prose are wrong, both fixed.

**One reliability finding is worth the reader's attention.** The four assertions `FloorSweep::close`
makes reduce, in practice, to *one* — and that one compares two zeroes on the six-accession slice,
which is the store this branch used at every earlier step. I construct below a defect that passes
all four and prints a whole distribution measured over a window the caller never builds. Three
cheap checks are added that close most of it, and the rest is named in the report so the next
reader is not misled about what the assertions cover.

**Nine of the ten deliberate defects I wrote into `FloorSweep` were caught by the six tests.** The
one survivor is not a numeric defect at all, and I explain why below rather than dressing it up.

Everything compiles and passes: `cargo test --all-features --example ng_window_coverage_probe` —
15 passed (13 before this review's two tests); `cargo test --lib --all-features` — 6,378 passed, 0
failed, 15 ignored; `cargo check --lib --tests --all-features` clean; `cargo clippy` unchanged at
3 pre-existing `needless_lifetimes`; `cargo fmt --check` clean on every file this review touched.

---

## 1. The central claim, checked against the accumulator

**Claim.** An arm at floor `F`, fed exactly the positions the walk beside it is fed, silences
exactly the windows that hold fewer than `F` distinct covered positions — so the bank's answers in
floor order are the cumulative distribution.

**It holds.** `src/ng/window_coverage/accumulator.rs` reads `min_window_positions` in exactly one
place, `finalise_centre`:

- the *absent* branch is the only producer of `WindowCoverage::absent()`, and it is entered on
  `distinct_positions_summed < min_window_positions` and nothing else;
- the other branch divides `sum_gc` and `sum_depth` by `summed_positions`, which the same function
  documents and the shape guarantees is `>= 1` (the centre lies in its own window), so that branch
  cannot produce a `NaN` and `is_absent()` cannot be true of a window that spoke;
- `ready.push_back(...)` and `centre_offset += 1` happen on **both** branches, so *how many*
  windows a stream finalises is one per buffered covered position whatever the floor is.

So across a bank of arms differing only in the floor, the sequence of centres and each centre's
`distinct_positions_summed` are bit-identical, and only the emit-or-silence decision differs. The
cases the brief asks about:

| case | why it does not break the claim |
|---|---|
| repeated positions | `distinct_positions_summed` counts coordinates, not records: incremented only when an entry differs from its predecessor, decremented only when the last entry at a coordinate leaves the buffer. Floor-independent. |
| `N` reference bases | not buffered, so they become no centre and join no window's sums; they still advance the frontier. Every arm sees the identical `(contig, position, base, depth)` call, so every arm skips the same ones. |
| contig change | `finalise_all()` then `reset_contig()`, driven by the position stream alone. |
| the tail closed by `finish` | `finish` calls `finalise_all` and hands back the **whole** `ready` queue, drained or not; `FloorSweep::close` counts all of it. |
| windows truncated at a contig edge | a truncated window simply holds fewer positions, which is the quantity being measured. Nothing special is needed. |
| the depth-bin-width fit | touches only `fold_or_hold_back`, the held-back list and `counts`. It never reaches the emitted `WindowCoverage` and never reaches `windows_under_the_floor`. A high-floor arm can end `EveryWindowUnderTheFloor` or `Unfittable`; the sweep discards each arm's histogram (`let (tail, _) = accumulator.finish()`), so that is invisible. |

**Sound.** No change needed.

## 2. Blocker / Major / Minor findings

### Major 1 — the shipped-floor check has teeth on four stores, not three (CONFIRMED, fixed)

The report says the check "has teeth on three of the five stores, where that number is 867,
212,850 and 222,673 rather than zero". The `human_genome_bottle` store's walk counts **28** absent
windows (`tmp/d1/sweep_hg002_bottle.txt`, `windows-under-the-floor 28`), which is also non-zero. The
check is vacuous on exactly one store — the six-accession slice, where every one of the six samples
reports 0.

Corrected in the report to "four of the five stores ... 28, 867, 212,850 and 222,673", with the
slice named as the one where it compares 0 with 0.

### Major 2 — the four assertions reduce to one, and that one is silent on the slice (CONFIRMED by construction, partly fixed)

The brief asks whether a defect can survive all four. **It can.** Build the arms with a wider
window than the run uses — `window_bp: 1000` in `FloorSweep::new`'s config literal — and:

- **assertion 1** (every arm finalised the same number of windows as the walk) passes. A window is
  finalised for every buffered covered position *whatever the configuration says*, so this
  assertion is blind to `window_bp`, to `gc_bins`, to `depth_bins` — to everything but a difference
  in the position stream itself;
- **assertion 2** (the arm at the shipped floor silenced what the walk counted absent) passes on
  the six-accession slice, where the walk silences 0 and a *wider* window silences 0 too;
- **assertion 3** (the walk finalised at least one window) passes;
- **assertion 4** (monotone in floor) passes — every arm shares the wider window, so the counts
  are still monotone.

The run then prints a full distribution of a window the caller never builds, over the one store
this branch has used at every earlier step. I verified the two load-bearing pieces rather than
reasoning alone: a *narrowing* drift is caught by the library, not by the sweep — I built the arms
with `window_bp: 250` and got `min_window_positions 300 is more positions than a 250-base window
can hold (251)` out of `WindowCoverageConfig::assert_valid`; a widening drift passes that check,
because `assert_valid` only bounds the floor from above.

**What was applied.** Three things, none of which fully closes it and all of which help:

1. **A fifth check with a value in it, not a comparison:** the arm at a floor of 1 must silence
   nothing, because a centre lies in its own window and no window holds zero covered positions.
   `CANDIDATE_FLOORS`'s own doc comment already stated this invariant; nothing asserted it. It is
   the only check in `close` that names a number rather than comparing two quantities that can both
   be zero.
2. **A sixth check that the floor-keyed checks can find their arm.** Both are written as "if this
   arm is the one at floor `F`", which on a hand-edited floor list that had stopped holding `F`
   would be no check at all rather than a failure. Now `close` asserts that arms at 1 and at
   `MIN_WINDOW_POSITIONS` exist.
3. **A printed row, `sweep-tied-to-the-walk`,** giving the shipped floor and what the tie to the
   walk actually compared on this store. A `0` there says the distribution below it is confirmed by
   nothing outside the bank — which is the honest reading of the slice's own output file.

**What is left open, and is now written into the report:** no runtime check can see the arms'
configuration. The unit tests are what pin it — `the_sweep_traces_a_distribution_whose_windows_hold_different_counts`
fixes the counts a 500-base window produces over 600 consecutive positions, and any width drift
moves them. And the practical guard is the one D1 in fact honoured: run the sweep over at least
one store whose walk silences something.

### Minor 3 — "6,000-fold" is 6,491-fold (CONFIRMED, fixed in three places)

The report, spec §3.3 and `MIN_WINDOW_POSITIONS`'s doc all say the answer moves 6,000-fold between
the 5.1 kb-interval store and the 122-base-interval store. That figure comes from dividing the two
*rounded* printed shares, 360.13 / 0.06 = 6,002. From the counts:

```
28 / 5,046,746  = 5.548 in 1,000,000
212,850 / 5,910,300 = 36,013 in 1,000,000
ratio = 6,491
```

Fixed to "6,500-fold", with both counts given beside it so the next reader does not have to
recompute it either.

### Minor 4 — spec §9 quotes the milder of the two short-interval rows (CONFIRMED, fixed)

§9 reads "It silences **at most** 1 window in 8,850 on intervals of 5 kb and longer, and 1 in 28 on
a run whose intervals average 122 bases". The first half is a worst case across three stores; the
second is the 30× row (360.13 in 10,000 → 1 in 27.8) while the same ground at 5× gives 381.69 in
10,000 → **1 in 26.2**. Fixed to "1 in 26 to 1 in 28".

### Minor 5 — `WindowCoverageConfig::validate` does not exist (CONFIRMED, fixed)

Deviation 2 names it; the method is `assert_valid`. Fixed.

### Minor 6 — "99.2 positions in 100 on tomato" is one of two tomato stores (CONFIRMED, fixed)

`positions-from-the-head-share` is 0.9920 on the slice and 0.9940 on the one-accession wide store,
and the 5 kb-interval human store's 0.9893 is not reported at all. Fixed to give all four.

### Minor 7 — `WindowTotals::sweep`'s comment describes the summation backwards (CONFIRMED, fixed)

"The first store's arms are cloned into by summing an empty bank into them" — the code is
`self.sweep.get_or_insert_with(FloorSweep::empty).add(store)`, which sums each store, the first
included, into a bank that starts empty. That is the better design and the reason given for it is
right; the sentence says the opposite of what happens. Rewritten.

### Minor 8 — a mutation that survives everything, and why it is not a defect (CONFIRMED)

Changing `while let Some(...) = accumulator.pop_ready()` to `if let` in `FloorSweep::observe`
survives all 13 tests and all four assertions. **This is a real gap in what the checks can see, but
it is not a wrong number:** `finish` returns the whole undrained `ready` queue as part of the tail,
which `close` counts, so `windows` and `silenced` come out identical. What it costs is memory —
each arm's queue would grow to the whole store, about 142 MB per arm at 24 bytes an entry over
5.9M windows, times 25 arms. Nothing here measures memory, so nothing here can catch it. Recorded
rather than fixed.

### Minor 9 — "about 630 reads" is 626 (not fixed)

122 positions × 5.13 reads a position = 626. The report says "about", and 626 rounds to 630 at two
significant figures. Left as written.

### Minor 10 — "isolated" is inferred, not measured (not fixed)

Of the 261 tomato windows holding fewer than 10 positions, the report says "Those are the isolated
covered positions the floor exists for". The probe reports counts, not where those windows sit or
what surrounds them. The inference is safe — a window holding under 10 of a possible 501 positions
*is* near-empty ground — but the word does work the measurement does not supply. The report is
otherwise scrupulous about this line ("this step measures the cost of the floor and not its
benefit", said explicitly), so this is a lone slip.

### Minor 11 — the depth column is a mean over one-base records only (not fixed)

The footnote says so exactly ("the mean `reads_compared_with_reference` at a record covering one
base"), so the report is honest. Spec §9's "five stores spanning 5.1 to 301 reads compared with the
reference a position" drops the qualifier. On the tandem-repeat store the wide tract records cover
359,907 bases — 6.1 in 100 of the positions — whose depth is not in that mean.

## 3. The six tests, put under deliberate defects

Each defect was written into the shipped file, built and run through
`./scripts/dev.sh cargo test --all-features --example ng_window_coverage_probe`, then restored with
`cp` from a copy taken before any edit.

| defect written into `FloorSweep` | caught? | by what |
|---|---|---|
| an arm's own floor off by one (`min_window_positions: floor - 1`) | **yes** | `the_sweep_traces_a_distribution...` (96 against 98) and `two_sweeps_sum...` (192 against 196) |
| silenced counted on `!window.is_absent()` in the `observe` drain | **yes** | the shipped-floor assertion, via `the_sweep_traces_a_distribution...` (350 against 0) |
| the tail not counted in `close` | **yes** | four tests; the window-count assertion fires first (350 against 600) |
| the shipped-floor check comparing `arm.silenced` with itself | **yes** | `a_sweep_disagreeing_with_the_walk_at_the_shipped_floor_is_refused` — "did not panic as expected" |
| `add` zipping the other bank's arms in reverse | **yes** | `add`'s own floor-order assertion, via `two_sweeps_sum...` (1 against 501) |
| `pop_ready` drained once per `observe` rather than looped | **no** | — see Minor 8; the counts are unchanged, only memory is |
| an arm built with `window_bp: 250` | **yes** | `WindowCoverageConfig::assert_valid` rejects floor 300 in a 250-base window |

**One gap the mutations exposed and the fix closes.** The wrong-condition defect in the `observe`
drain was *not* caught by `an_arms_silenced_count_is_the_windows_holding_fewer_positions_than_its_floor`,
the test the report calls "the claim the whole measurement rests on". Ten consecutive positions
never reach any centre's right edge, so every window there comes back through `finish` and the
`observe` drain is never exercised at all; only the 600-position test reaches it, where 350 of the
600 centres close before the stream ends. The counting rule is now written once
(`count_one_window`) and called from both drains, so it cannot be right on one path and inverted on
the other, and the flagship test's doc says which path it exercises.

## 4. The hand-computed expectations, re-derived

`half_window_bp = 500 / 2 = 250`, so the window centred at `p` spans `[p − 250, p + 250]` — 501
positions when nothing clips it. Over 600 consecutive covered positions the count at centre `p` is
`min(600, p+250) − max(1, p−250) + 1`:

| centres | count |
|---|---|
| `p` in 1..=250 | `p + 250`, so 251 up to 500 |
| `p` in 251..=350 | 501 |
| `p` in 351..=600 | `851 − p`, so 500 down to 251 |

- **floor 250 → 0.** The thinnest window is 251 (at `p = 1` and at `p = 600`), which clears 250.
  The test's comment, "the one centred at base 1, holding bases 1 to 251", is right.
- **floor 300 → 98.** `p + 250 < 300` for `p` in 1..=49, and `851 − p < 300` for `p` in 552..=600.
  Forty-nine at each end.
- **floor 501 → 500.** Only the 100 centres at 251..=350 hold the full 501.

**All three correct.** The ten-position test likewise: every window over 10 consecutive positions
holds all ten, so floors 1..10 silence none and floors 15, 20 and 501 silence all ten, and "ten
positions, ten centres" is right.

## 5. Every number in the report and the doc comment

Reproduced from `/Users/jose/devel/pop_var_caller-window-coverage/tmp/d1/sweep_*.txt` and from the
BED files under `benchmarks/`. Everything below is **as claimed** unless the row says otherwise.

**The table** (all five rows, all four columns, both in the report and in `MIN_WINDOW_POSITIONS`'s
doc):

| store | windows | reads a position | under 50 | under 100 | under 200 |
|---|---|---|---|---|---|
| tomato slice | 1,172,242 ✓ | 14.37 → 14.4 ✓ | 0.00 ✓ | 0.08 ✓ | 8.10 ✓ |
| tomato wide | 7,667,288 ✓ | 10.29 → 10.3 ✓ | 1.13 ✓ | 1.69 ✓ | 15.39 ✓ |
| HG002 bottle | 5,046,746 ✓ | 301.36 → 301.4 ✓ | 0.06 ✓ | 6.68 ✓ | 16.18 ✓ |
| HG002 TR 30× | 5,910,300 ✓ | 30.30 → 30.3 ✓ | 360.13 ✓ | 4,393.07 ✓ | 5,848.64 ✓ |
| HG002 TR 5× | 5,833,898 ✓ | 5.13 → 5.1 ✓ | 381.69 ✓ | 4,432.50 ✓ | 5,913.10 ✓ |

**Interval lengths, from the region files the build scripts name:**

- `benchmarks/ssr_hg002/regions/HG002_GRCh38_TandemRepeats_v1.0.1_Tier_50000.bed` — 50,000
  intervals, 6,089,411 bases, **mean 121.79** → "50,000 × 122 bp" ✓. Its first three columns,
  sorted, are identical to `tmp/d1/hg002_all.bed`, which is the file the store was actually built
  over — so the report's naming of the BED it checked against is right.
- `tmp/d1/hg002_bottle.bed` — 1,000 intervals, 5,078,968 bases, **mean 5,079** → "1,000 × 5.1 kb" ✓.
- `benchmarks/tomato1/regions.bed` — 80 intervals, 8,000,000 bases, mean 100,000 → "80 × 100 kb" ✓.
- The six-accession slice — `tmp/d_baseline/psps.log` records "2 intervals ... 200000 bases" →
  "2 × 100 kb" ✓, and the same log's "records written: 2311" matches the oracle's 2,311 ✓.
- "the intervals shorten 42-fold" — 5,079 / 121.79 = 41.7 ✓.

**Derived figures:**

- "**6 in 100**", the 5×-against-30× move: 381.69 / 360.13 = 1.0599 ✓ as a relative change. The
  sentence did not say relative; it now gives both figures.
- "**6,000-fold**" — **wrong, it is 6,491**. Finding Minor 3.
- "**1 window in 8,850**" — 10,000 / 1.13 = 8,850 ✓ (from the counts, 7,667,288 / 867 = 8,844).
- "**261**" windows under 10 positions on the tomato wide store ✓.
- "**606** would speak at a floor of 10" — 867 − 261 = 606 ✓.
- "**44 in 100**" at floor 100 on short-interval ground — 43.93 and 44.33 ✓.
- "**3.6 and 43.9 in 100**" of windows silenced there — 360.13 and 4,393.07 per 10,000 ✓.
- "**2.9 in 100** of requested bases in intervals under 50 and **42.8 in 100** under 100" — from
  the BED: 178,102 / 6,089,411 = 2.9248% and 2,606,472 / 6,089,411 = 42.8034% ✓.
- "**24,793,279** records covering one base" — 1,162,872 + 7,621,310 + 4,992,697 + 5,544,484 +
  5,471,916 = 24,793,279 ✓. The tomato pair alone sums to 8,784,182, which is exactly the count B2
  reported, consistent with B2 having walked the same two stores ✓.
- "0 disagreements, 0 partial witnesses, no one-base tract record" — zero in all five files ✓.
- "**93.8 in 100** on the tandem-repeat store" — 0.9381 and 0.9380 ✓. "99.2 on tomato" — finding
  Minor 6.
- "**630 reads**" — 626. Finding Minor 9.
- "0 windows of 1,172,242, in all six samples" ✓ — each of the six per-sample blocks reports 0.
- "**501** is the highest arm ... a 500-base window spans 501 positions" ✓ —
  `assert_valid` computes `2 * (window_bp / 2) + 1`.
- "a bank of **25** accumulators" ✓ — `CANDIDATE_FLOORS` has 25 entries and already contains 50, so
  the merge with the shipped floor dedups back to 25.

**Validation results**, all re-run in the container on this worktree:

- `cargo test --lib --all-features` — 6,378 passed, 0 failed, 15 ignored ✓ exactly as claimed.
- `cargo clippy --lib --all-features` — 3 warnings, all `needless_lifetimes`, at
  `src/ng/run/cohort_merge/build.rs:821`, `build.rs:895` and `serial.rs:69` ✓.
- `cargo test --all-features --example ng_window_coverage_probe` — 13 passed as the step was
  written ✓ (15 with this review's two tests).
- The standing oracle's sha256 and its 13,866/13,589 window rows were **not** re-run; they need the
  reference and a calling run, and nothing in the diff can plausibly move them (the only `src/`
  edit is a doc comment). Its record count is corroborated: `tmp/d_baseline/psps.log` reports
  "records written: 2311".

## 6. The mechanism claim — is the arithmetic capable, and does the check beg the question?

**The reason is arithmetically capable of the observed sizes.** An interval of length `L ≤ 251`
containing a centre `p` lies wholly inside `[p − 250, p + 250]`, because every position of it is at
most `L − 1 ≤ 250` from `p`. So on 122-base intervals a window holds about 122 positions, and the
floor is asking how short an interval is too short. If every window held exactly its own interval's
length, floor 50 would silence the 2.92 in 100 of bases in intervals under 50 and floor 100 the
42.80 in 100 under 100 — against 3.60 and 43.93 measured. **Right size, right direction.**

**The check does not beg the question.** The BED figures are computed from the region file alone
and never touch the store, the walk or the accumulator, so a shared misunderstanding of what a
window spans cannot cancel out.

**It was under-specified in one way, which I measured and added to the report.** Comparing
"interval shorter than `F`" with "window silenced at `F`" assumes a window never reaches into a
*neighbouring* interval within 250 bases, which would lift a count above its own interval's length.
Computing each requested base's full 501-base neighbourhood against the whole BED — still from the
BED alone — moves the prediction only from 2.92 and 42.80 to **2.89 and 42.48 in 100**, so the
assumption costs at most a third of a point and the check stands. The remaining excess over the
measured 3.60 and 43.93 is the report's stated cause, the positions no read reached, and it is the
right size: 5,910,300 of the 6,089,411 requested bases were covered at 30 reads a position, so 2.9
in 100 were not. The same reading explains the top of the distribution — the BED predicts 87.35 in
100 of windows holding fewer than 501 positions, and 91.86 was measured.

## 7. The ruling

**Supported by what was measured, and honest about what was not.**

- "on intervals of 5 kb and longer, at most 1 window in 8,850" — measured, three stores.
- "raising it to 100 would silence 44 in 100 on short-interval ground" — measured, two stores.
- "nothing measured here says such a window is wrong" — correctly declared *not* measured, twice,
  with the benefit assigned to the hidden-duplication filter's own branch. This is the place a
  report of this kind usually overreaches and this one does not.
- "lowering it to 10 would let 606 windows speak" — measured, one store; the word "isolated"
  attached to them is inferred (Minor 10).
- The one arm the measurement lacks — patchy coverage over *continuous* analysed ground — is named,
  with the reason nothing in this tree can build it and the nearest substitute identified.

Nothing else in the ruling is asserted as measured that was inferred.

## 8. `the_configuration_a_run_uses()`

**The property a previous review praised is preserved exactly.** All six fields are spelled and
there is no `..` in the function, so a field added to `WindowCoverageConfig` still has to be
answered for in one place before the file compiles.

**The sweep's `..the_configuration_a_run_uses()` does not reopen that hole**, because a new field
gets its value in that same function. The residual risk is a different one, and it is Major 2's
mechanism: the `..` also means that a *hand-added override inside the arm literal* — someone
narrowing the window while debugging and forgetting to take it out — is visible nowhere. No runtime
check in `close` can see it. The tests can, and now say so.

## 9. What I changed, so it can be lifted

All four files are in this review's worktree at commit `444672b0`; `git diff HEAD` shows 171
insertions and 30 deletions.

**`examples/ng_window_coverage_probe.rs`**

1. `count_one_window(windows, silenced, window)` — the counting rule written once, called from both
   `FloorSweep::observe`'s drain and `FloorSweep::close`'s tail. A free function taking the two
   counters rather than a method on `OneFloor`, because `observe` is holding the same arm's
   accumulator across the drain; `observe` now destructures the arm exhaustively so the borrows
   split.
2. `FloorSweep::close` gains two checks: **the arm at a floor of 1 silenced nothing**, and **arms
   exist at both floors the floor-keyed checks are written against**. Its doc says six things
   rather than four.
3. `FloorSweep::print` gains a `sweep-tied-to-the-walk` row, giving the shipped floor and what the
   tie to the walk actually compared on this store — `0` says the rows below it are confirmed by
   nothing outside the bank.
4. Two tests: `an_arm_at_a_floor_of_one_silencing_a_window_is_refused` and
   `a_sweep_that_lost_the_arm_a_check_is_written_against_is_refused`.
5. `an_arms_silenced_count_...`'s doc now says which drain it exercises and which test covers the
   other.
6. `WindowTotals::add`'s comment rewritten to describe the direction the summation actually runs.

**`doc/devel/reports/implementations/ng_window_coverage_d1_2026-09-07.md`** — "three of the five
stores" → four, with 28 added and the vacuous store named; "6,000-fold" → 6,500-fold with both
counts; "6 in 100" given its two figures; the neighbourhood-aware BED prediction and the
uncovered-position arithmetic added to the mechanism paragraph; `validate` → `assert_valid`; "99.2
on tomato" → all four head shares; the four-assertion list rewritten as six, with a new paragraph
saying plainly what none of them can see; test count six → eight and the two new rows; probe test
count 13 → 15.

**`doc/devel/ng/spec/window_coverage.md`** — §3.3's "6,000-fold" → 6,500-fold with the two shares;
§9's "1 in 28" → "1 in 26 to 1 in 28".

**`src/ng/window_coverage/mod.rs`** — `MIN_WINDOW_POSITIONS`'s doc: "6,000-fold" → 6,500-fold with
the counts, and "6 in 100" given its two figures. **Still a doc comment only** — the step's claim
that no library code changed behaviour survives this review's edits.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` was touched.

## 10. What I looked at and found sound

- `WindowCoverageAccumulator::finalise_centre`, `observe`, `pop_ready`, `finish`, `finalise_all`,
  `reset_contig`, `fold_or_hold_back` and `fit_depth_bin_width_and_fold_held_back` — read in full
  against the claim in §1 above, including the two-pointer invariants. The shrink-left step's
  `centre_offset -= 1` cannot underflow: a front is popped only when it sits strictly left of the
  current centre, so it is a different, earlier entry and `centre_offset >= 1`. The
  distinct-coordinate decrement fires only for the last entry sharing a coordinate, which is right
  for a buffer that is sorted and whose co-located entries are always summed together.
- **Ordering in `main`.** `windows.close()` runs before `windows.print(...)` and before
  `WindowTotals::add`, so the totals include every arm's tail, and a failed assertion stops the run
  before any sweep row is printed.
- **`FloorSweep::add`** — the per-arm floor assertion makes a mis-zip loud, and the exhaustive
  destructure of `OneFloor` in the loop means a field added later has to be answered for.
- **`FloorSweep::empty`** — arms with no accumulator; `close` on such a bank would panic at its
  `expect`, and nothing calls it. Its floors come from the same `the_floors_asked_about()` as
  `new`'s, so a total's rows are the store rows in the same order.
- **`print`'s denominator** is the first arm's window count, which assertion 1 has made equal to
  every other arm's on every store summed in.
- **Argument handling** — `--covered-positions-per-window` without `--reference` is refused with a
  message that says why, matching the existing `--windows-from-the-run` guard.
- **The floor list** — ascending, deduped, 25 entries, top at 501 which is exactly the widest
  `assert_valid` allows; the shipped floor is merged in at runtime, so moving the default cannot
  quietly drop the check written against it.
- **Report structure** — the deviations are all recorded, the plan's "tomato slice and HG002" is
  correctly reported as five stores with the reason, and the arm the measurement lacks is named
  rather than glossed.
