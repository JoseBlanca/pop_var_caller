# window coverage — D2 review: reliability, correctness, and every number recomputed

**Date:** 2026-09-07
**Reviewing:** commit `f3228dff` (`wip: D2 for review`) against `38e06a3f` (D1)
**Remit:** reliability, correctness, error handling, and a recomputation of every figure the step
claims. Naming and structure are a second reviewer's.
**Raw output checked against:** `/Users/jose/devel/pop_var_caller-window-coverage/tmp/d2/scheme_*.txt`
(five files, ten sample-stores)

## The short of it

**No blocker, and the step's conclusion survives intact.** All three constants are kept for
reasons the raw output supports: at the shipped range the worst of ten sample-stores overflows
1,972 windows of 7,666,421, and every arithmetic step from that number to "780 times under the
fit's guard" checks out.

Three things were wrong or unpinned, and all three are now fixed in this worktree:

- **one wrong count in the prose, repeated in three files** — the prefix median comes out below
  the whole-store median on **eight** of the ten sample-stores, not nine (Major 1, CONFIRMED);
- **the overflow fraction's denominator was not the fit's, and nothing would have noticed** — the
  probe divided by the accumulator's `windows_folded` counter rather than by the cells, and two
  deliberate defects in that denominator survived every test (Major 2, CONFIRMED by mutation);
- **the fit's guard is a strict `>` and the probe's `>` was untested** — flipping it to `>=`
  survived all 21 tests (Minor 4, CONFIRMED by mutation).

The step now has 22 tests, and the two mutations above are killed by the one added.

## 1. Findings

### Major 1 — "nine of the ten" is eight of the ten (CONFIRMED, fixed in three files)

The report, spec §3.4 and `DEPTH_SCALE_WINDOWS`'s doc all say the median fitted from the first
10,000 windows "came out below the median over every window of the same store in nine of the ten".
Recomputing the ratio for every sample-store from `scheme_*.txt` — the `as-the-run-fits-it` arm's
median over the `scale-sample=every-window` arm's — gives eight below 1 and two above:

| store | sample | prefix median | whole-store median | ratio |
|---|---|---|---|---|
| tomato, 80 × 100 kb | SRS3394712 | 6.4092 | 9.7725 | 0.656 |
| HG002, GIAB high-confidence | HG002 | 246.2080 | 302.6613 | 0.814 |
| tomato slice | SRS3394712 | 6.4092 | 7.8463 | 0.817 |
| HG002, tandem-repeat tiers, 5× | HG002 | 4.0318 | 4.9301 | 0.818 |
| tomato slice | SRS3394713 | 10.1238 | 12.3234 | 0.822 |
| HG002, tandem-repeat tiers, 30× | HG002 | 25.1182 | 30.1525 | 0.833 |
| tomato slice | SRS3394711 | 13.5868 | 14.5988 | 0.931 |
| tomato slice | SRS3394712_SRR7279484 | 3.9721 | 4.2395 | 0.937 |
| tomato slice | SRS3394606 | 27.6108 | 27.0160 | **1.022** |
| tomato slice | SRS3394714 | 21.6228 | 20.0200 | **1.080** |

The report's own table prints those two ratios as 1.02 and 1.08, so the table and the sentence
above it contradict each other. The "13% on average" that follows is right and is a mean over all
**ten** ratios (mean 0.8728, so 12.7% shallow); over the eight that are below it would be 17.2%.

**Why it matters beyond the count.** The sentence is the premise of the step's load-bearing
argument — that the prefix is systematically shallow, so the range's headroom is what makes it
safe. Eight of ten still carries that argument; nine of ten overstates a systematic effect that
two stores contradict, and those two are the *shallow-interval* tomato accessions where a reader
would most want to know the bias can go the other way.

**Fixed** in the report, in spec §3.4 and §9, and in `DEPTH_SCALE_WINDOWS`'s doc comment, with the
two exceptions named and sized (2% and 8% above).

### Major 2 — the overflow fraction was divided by the accumulator's counter, not by the cells (CONFIRMED by mutation, fixed)

The report says the fraction is computed "**the way `coverage_model.rs` computes it** — the
overflow column summed over GC rows, over the regular bins plus that column". The code did not do
that. `FittedAxis::overflow_fraction` was

```rust
self.overflowed as f64 / self.folded as f64          // folded == hist.windows_folded
```

where production (`coverage_model.rs:244-249`) is

```rust
let regular_total: u64 = layout.depth_marginal(hist).iter().sum();
let overflow_total = layout.overflow_total(hist);
let overflow_fraction = overflow_total as f64 / (regular_total + overflow_total) as f64;
```

The two numbers **do** agree on all ten sample-stores — I checked the identity against the raw
output: `windows-finalised − windows-under-the-floor` equals every arm's printed `folded` on all
five files (7,667,288 − 867 = 7,666,421; 5,910,300 − 212,850 = 5,697,450; 5,833,898 − 222,673 =
5,611,225; 5,046,746 − 28 = 5,046,718; and each slice sample likewise). **So no published figure
changes.** They can part in exactly one way, which the module's own doc names: `fold` writes cells
with `saturating_add` but increments `windows_folded` unconditionally, so a cell that saturated at
`u32::MAX` would leave the counter high and every overflow fraction reading **low** — the unsafe
direction, since low is what makes a sample look acceptable. That needs one GC-and-depth cell to
hold 4.3 billion windows, so a reference above about 4.3 Gbp: larger than tomato or human, not
larger than the plant genomes this caller commits to taking.

**The reliability finding is not the saturation case — it is that nothing was checking.** Two
deliberate denominators survived the whole suite (see the mutation table): dividing by
`folded + overflowed`, and folding `windows_under_the_floor` into `folded`. The second is the one
that would bite: on the tandem-repeat store 212,850 of 5,910,300 windows are silenced by the floor,
so counting them would divide every fraction by 1.037 and understate the overflow by 3.6% on
exactly the store where short intervals make windows thin.

**Fixed.** `FittedAxis` now carries `in_a_regular_bin` — summed out of the cells, production's
`regular_total` — beside `overflowed`, and `folded()` is their sum. The identity with
`windows_folded` is asserted where the axis is read, with a message naming saturation as the only
way it can fail, so the one case that could silently lower a fraction stops the run instead.

### Minor 3 — "16 MB at a million windows" is derived differently from the "262 kB at 10,000" it is compared with (CONFIRMED, fixed)

The report and `DEPTH_SCALE_WINDOWS`'s doc say holding a million windows back is "16 MB of `f64`
pairs a sample — 16 GB across a thousand — against the 262 kB at 10,000 that spec §3.4 already
charges". The 262 kB is **capacity**: spec §3.4 and `accumulator.rs`'s type doc both say so
explicitly — 10,000 pairs is 160 kB of live data, the `Vec` grows by doubling, capacity reaches
16,384, and 16,384 × 16 = 262,144 bytes. The 16 MB is **length**: 1,000,000 × 16. Derived the same
way, a million windows held back reaches capacity 1,048,576, which is **16.8 MB** a sample and
16.8 GB across a thousand.

The argument is unaffected — the honest figure is worse, not better — but the sentence compares two
numbers taken on different bases in the same breath. Fixed in both places, with the basis stated.

### Minor 4 — the fit's guard is a strict `>` and nothing pinned it (CONFIRMED by mutation, fixed)

`the_fit_would_reject_it` uses `> THE_FITS_OVERFLOW_GUARD`, which is right: production rejects on
`overflow_fraction > cfg.max_overflow_fraction`, so a sample sitting exactly on 0.20 is **taken**.
Changing the probe's `>` to `>=` survived all 21 tests. The three arms the tests exercise sit at
0.0, 0.3 and 0.657 — nothing near the boundary.

**Fixed** by the added test, which constructs a stream whose share is exactly 200 in 1,000 and
asserts the sample is taken.

### Minor 5 — `WindowCoverageConfig` still calls the constants provisional after D1 and D2 settled them (CONFIRMED, fixed)

The step's report says the three constants "stop being marked soft". They do, at each `const`. But
`WindowCoverageConfig`'s own doc still read "`MIN_WINDOW_POSITIONS`, `DEPTH_BINS`,
`DEPTH_SCALE_WINDOWS` and `DEPTH_RANGE_IN_MEDIANS` are ng's own and are provisional: a reader who
finds a run's behaviour turning on one of them is looking at a value nobody has yet defended with
data, and plan steps D1 and D2 are where that happens" — with D1 and D2 both landed — and three of
its fields still said "Provisionally". A reader arriving at the type rather than at the constants
would be told the opposite of what the step concluded. Fixed in all four places.

### Minor 6 — "4.0 to 302.7 reads a position" mixes two different quantities (CONFIRMED, fixed)

Assumption 1 describes the ten sample-stores as spanning "4.0 to 302.7 reads a position". 302.66 is
a whole-store median; 4.03 is a *prefix* median (HG002 at 5×), whose whole-store median is 4.93.
Stated on one basis it is 3.97 to 246.21 (prefix medians — the ones the shipped setting produces,
and the pair the report quotes as "62-fold" elsewhere) or 4.24 to 302.66 (whole-store medians).
Fixed to the whole-store pair, with the prefix pair named beside it.

### Minor 7 — "780 times the headroom it needs" is a ratio of overflow shares, not of headroom (CONFIRMED, fixed)

"A prefix 34% shallow makes the axis 34% short, and a range with 780 times the headroom it needs
absorbs that without the fit noticing." The 780 is `0.20 / 0.000257` — the guard over the *measured
overflow share*, not a margin in depth. The claim it supports is sound and is a measurement, not an
inference (see §4 below), but the sentence asks the reader to read a share ratio as a depth margin.
Rewritten to say what was measured: the whole cost of being 34% short on that store is 1,972
windows, 2.6 in 10,000 against a guard at 2,000 in 10,000.

### Minor 8 — "Tests added: Five" over a table of six (CONFIRMED, fixed)

The "Changes made" section says six tests and the validation line says 21 passed against 15 before,
which is six. Only the "Tests added" heading says five. Now seven, with the one this review added.

### Minor 9 — the `every-window` arm's held-back list is unbounded, and nothing says so (CONFIRMED, not fixed)

`BinSchemeName::ScaleSampleOf(u32::MAX)` sets `depth_scale_windows` to a number no store reaches,
so that arm holds **every** window back until `finish`. At 16 bytes a window with a doubling `Vec`
that is about 134 MB on the 7.67-million-window tomato store, plus 16.8 MB for the
`ScaleSampleOf(1_000_000)` arm — around 151 MB of transient beside the walk. Stores are walked one
at a time, so it does not multiply by the six-sample slice, and it plainly ran. Left as is: this is
the arm's whole purpose, the probe is not production, and the run completed. Recorded so that a
later reader pointing this sweep at a whole human genome knows what the last two arms cost.

## 2. The mutation table

Each defect was written into the shipped file, built and run through
`<root>/scripts/dev.sh cargo test --all-features --example ng_window_coverage_probe`, then restored
with `cp` from a copy taken before any edit. Baseline: 21 passed.

| defect written into `BinSchemeSweep` / `OneBinSchemeAnswer` | caught, before the fixes? | by what |
|---|---|---|
| overflow column read one index low (`+ depth_bins - 1`) | **yes** | `the_range_decides...` (0 against 300) and `a_scale_sample_that_closes...` (0 against 19,700) |
| cell index transposed to column-major (`depth_bins * gc_bins + gc_bin`) | **yes** | the same two tests, same numbers |
| median recovered with the run's range instead of the arm's | **yes** | `the_range_decides...` — the `RangeOf(2.5)` arm reports a median of 0.25 |
| arm-to-config mapping shifted by one | **yes** | three tests; `a_scale_sample_that_closes...` reports a median of 0.25 |
| the walk-tie assertion comparing an answer with itself | **yes** | `a_bin_scheme_sweep_disagreeing_with_the_walk_is_refused` — "did not panic as expected" |
| `the_fit_would_reject_it` using `>=` rather than `>` | **NO — survivor** | nothing; Minor 4. Now caught by the added test |
| the fraction divided by `folded + overflowed` | **NO — survivor** | nothing; Major 2. Now caught by the added test |
| `folded` set to `windows_folded + windows_under_the_floor` | **NO — survivor** | nothing; Major 2. Now structurally impossible — `folded()` is derived from the cells — and the added test pins it |

After the fixes, the three survivors were re-run: `>=` fails
`the_overflow_share_is_over_the_windows_in_the_cells_and_the_guard_is_not_reached_at_it`; dividing
by `in_a_regular_bin` alone fails the same test; and the third cannot be expressed any more without
tripping the new cells-against-counter assertion.

**What the survivors have in common** is the answer to check 3, below: every one of them is a
misreading of a histogram that both sides of the walk-tie share.

## 3. The three assertions in `BinSchemeSweep::close`

The tie is **stronger than D1's** in the way D1's review asked about, and blind in a different way.

D1's Major 2 was that its four assertions could not see the arms' *window width*: an arm built with
`window_bp: 1000` finalises the same number of windows, so every check passed. D2's tie compares
`as_a_row()`, which carries the fitted median, the bin width, the top of the range, the folded and
overflowed counts and the fraction — so a width drift in the run-configured arm changes the mean
depths it folds, changes the fitted median, and shows up. And unlike D1's, it is never a comparison
of two zeroes: all ten sample-stores fit a real axis (every `histogram-depth-bin-width` line in the
raw output is non-zero).

**What it cannot see.** Both sides of the comparison are read by the *same* function,
`OneBinSchemeAnswer::of`, so any defect in reading a histogram applies identically to both and
cancels. That is exactly what the first mutation in the table demonstrates: reading the overflow
column at the wrong index left the tie green and was caught only by the unit tests. The tie's real
subject is narrower than the prose suggests — it checks that the sweep was fed the same positions
as the walk, and nothing else. That is worth having and is what the sweep's correctness most
depends on, but the report's "the one place this sweep and the walk answer the same question"
should be read as "the one place it is checked that they saw the same stream".

**On the formatting.** `{:.4}` on the median, `{:.6}` on the width and `{:.5}` on the fraction
could hide a difference below about 5e-7 in the width. In practice they cannot hide anything that
matters, because `folded` and `overflowed` are printed exactly, and any difference in the position
stream — the only difference the tie can see at all — moves `folded` by a whole window.

**The other two assertions.** "Some arm is `AsTheRunFitsIt`" guards a check written as "if this is
that arm", which would otherwise become vacuous on an edited list; that is D1's Minor lesson
applied, and it is right. "Some arm folded a window" refuses an empty store. Neither is
weakened by anything I found.

## 4. The load-bearing argument, checked

**"A prefix 34% shallow makes the axis 34% short."** Arithmetically exact and not an approximation:
the width is `median × range / bins` and the top of the range is `width × bins`, so the top is
`median × range` — linear in the median. 64.09 / 97.72 = 0.656, the same 0.656 as the median ratio.
✔

**"The entire overflow of the worst store — those 1,972 windows — is attributable to the short
scale sample: with the width fitted from every window instead, that store overflows nothing."**
CONFIRMED from `scheme_tomato_wide.txt`: `scale-sample=every-window` and `scale-sample=1000000`
both print `overflowed=0`, against `overflowed=1972` at the run's setting. ✔ Worth adding, and now
added to the report: the effect is **not monotone in the scale sample** — 2,995 at 100 windows,
1,945 at 1,000, 1,972 at 10,000, 2,017 at 100,000, then 0 at a million. The step's stronger claim,
that no *affordable* prefix fixes it, is what the 2,017 at 100,000 windows actually supports, and
that is the claim the report makes.

**"The two constants cannot be moved independently."** Follows from the two above and is a
judgement about a future step, not a measurement. Sound as stated.

## 5. Every number, recomputed

All from `scheme_*.txt`. **Right unless marked.**

**The range table.** Worst overflow share across all ten sample-stores at each setting:
2.5 → 0.11026 (tomato_wide) = 1,102.6 in 10,000 → **1,103** ✔, and 0.20/0.11026 = 1.81 → **1.8** ✔.
5 → 0.00349 = 34.9 → **35** ✔, 0.20/0.00349 = 57.3 → **57** ✔.
10 → 1,972/7,666,421 = 0.00025723 = 2.57 → **2.6** ✔, 0.20/0.00025723 = 777.5 → **780** ✔.
20 and 40 → every arm `overflowed=0` ✔.
Top-of-axis column, read off the 5.1× store (`scheme_hg002_5x.txt`): 10.08, 20.16, 40.32, 80.64,
161.27 → **10.1 / 20.2 / 40.3 / 80.6 / 161.3** ✔.

**"1,972 windows of 7,666,421"** ✔ verbatim. **"nine of the ten sample-stores overflow nothing at
all"** at the shipped setting ✔ — six slice samples, both tandem-repeat stores and the bottle store
all print `overflowed=0`; only tomato_wide does not. (This "nine" is right; the wrong one is Major
1's, about the median table.) **"a guard that fires at 2,000 in 10,000"** ✔ — 0.20.

**"614 in 10,000"** ✔ — the slice's SRS3394712 at range 2.5, 0.06144. **"Both are tomato at 6 to 10
reads a position"** ✔ — 6.41 prefix / 7.85 and 9.77 whole-store.

**"thirteenfold"** ✔ — 780/57 = 13.7.

**The median table**: all ten prefix medians, whole-store medians and ratios reproduce to the digits
printed (table in Major 1). **"34% at worst"** ✔ (1 − 0.656). **"13% on average"** ✔ over all ten
(1 − 0.8728 = 12.7%). **"nine of the ten"** ✘ — Major 1.

**"5.25, 6.50, 6.41, 6.25 … 9.43"** ✔ — tomato_wide at scale samples 100 / 1,000 / 10,000 (the
run's own arm) / 100,000, then 1,000,000. **"100,000 windows is 100 kb, which is one of that store's
eighty intervals; a million is ten of them"** ✔ — that store finalises 7,667,288 windows over 80
intervals, so a window is a position and 100,000 is about one interval.

**"16 MB of f64 pairs a sample", "16 GB across a thousand"** ✘ — Minor 3; 16.8 MB and 16.8 GB on
the basis the 262 kB uses. **"262 kB at 10,000"** ✔ — 16,384 × 16 bytes, and spec §3.4 and
`accumulator.rs` both state the doubling.

**"3.97 to 246.21 … 62-fold"** ✔ — 246.208/3.9721 = 61.99. **"bin 7"** ✔ — 3.9721/0.5 = 7.94,
floor 7. **"bin 492 of 2,000"** ✔ — 246.208/0.5 = 492.4, floor 492, and production's scheme is
2,000 bins of 0.5 to 1,000×. **"three quarters of the axis above"** ✔ — 492/2,000 = 24.6%.
**"bin 40 of 400 by construction"** ✔ — `median / (median × 10 / 400)` = 40 identically.
**"50 × 401 × 4 bytes, 80.2 kB"** ✔ = 80,200. **"400 kB a sample"** for production ✔ =
50 × 2,001 × 4 = 400,200.

**"ten sample-stores", "the same five stores D1 used"** ✔ — 6 + 1 + 1 + 1 + 1. **"three interval
lengths"** ✔ — 100 kb, 5.1 kb, 122 bp. **"4.0 to 302.7 reads a position"** ✘ — Minor 6.

**"Ten arms — the run's own configuration, five scale samples … and four ranges"** ✔ = 1 + 5 + 4.

**Validation section.** `cargo test --lib --all-features` — I ran it: **6,378 passed, 0 failed, 15
ignored** ✔ exactly as claimed. `--example ng_window_coverage_probe` — **21 passed** on the
unmodified tree ✔, and 15 tests at `38e06a3f` ✔. Clippy — **3 warnings, all `needless_lifetimes`
in `src/ng/run/cohort_merge/`** ✔. `rustfmt --check` — the two files this step touched are clean
✔ (nine other files in the tree are not, all predating this step).

**"Tests added: Five"** ✘ over a table of six — Minor 8.

## 6. What I looked at and found sound

- **`median_depth` is recovered, not stored, and the recovery inverts the fit.** The accumulator
  computes `width = median × range / bins`; the probe computes `median = width × bins / range`.
  Four roundings separate them, so exactness is not free — I measured it over 400,000 random
  medians in 0.5 to 400 reads a position and the worst relative error is **1.8e-16, one ulp**. At
  the deepest store that is 4e-14 reads against a value printed to four decimals. The tests' 1e-9
  tolerances are the right shape and are not hiding anything.
- **It uses the arm's own range everywhere.** `named.configuration().depth_range_in_medians`
  returns the moved field for a `RangeOf` arm and the run's for the other nine, which is correct
  for all ten; and the walk's histogram is read as `AsTheRunFitsIt`, which is the configuration the
  walk in fact used. The mutation that swaps in the run's range for the arm's is caught.
- **The sweep is fed exactly the walk's positions.** One loop, one `observe` per arm per position,
  the same four values in the same order, with an exhaustive destructure of
  `WindowRecomputation`'s fields so a new field cannot go silently unfed. Closed after the walk's
  own histogram is taken, which is what the tie needs.
- **The six new tests' exact counts and medians re-derive correctly** from the accumulator's real
  behaviour, contig edges included. A contig of *N* ≥ 251 consecutive covered positions finalises
  *N* windows, the thinnest of them (at the first and last position) holding 251 covered positions
  — five times the floor of 50 — so nothing near a contig edge is silenced and every window's mean
  depth is exactly its contig's. Hence 600, 700 + 300 = 1,000, and 10,300 + 19,700 = 30,000, all
  as asserted. The medians follow from `median_depth`'s upper-middle rule, which takes the element
  at index `len / 2` of the sorted depths: at 1,000 windows that is index 500, which falls among
  the 700 windows at depth 1, so the median is 1 ✔; at 30,000 windows it is index 15,000, which
  falls among the 19,700 at depth 100 ✔; and the run's arm fits at the 10,000th held-back window,
  all of which
  are still on the first contig at depth 1 ✔ — the last 250 windows of that contig do not finalise
  until the contig changes, which does not disturb the count.
- **Clamping cannot lose a window.** `cell_index` clamps GC into the last GC bin and depth into the
  overflow column, so a clamped window still reaches a cell and is still counted; and an absent
  window never reaches `fold_or_hold_back` at all, so no `NaN` can enter the histogram.
- **`SampleHistogram`'s three silences are handled**, and an arm that is silent while another
  fits — reachable, if a store's first 10,000 windows have a zero median but the whole store does
  not — prints as a named silence rather than as a fitted zero.
- **The `--bin-scheme` flag is opt-in** (`bin_scheme.then(BinSchemeSweep::new)`), so a walk that
  did not ask for the sweep pays none of its ten accumulators.

## 7. Exactly what I changed

Four files. No library behaviour changed: the only `src/` edits are doc comments. Nothing under
`src/sample_summary/`, `src/paralog/` or `src/var_calling/` was touched.

**`examples/ng_window_coverage_probe.rs`**

1. `FittedAxis::folded` (a copy of `hist.windows_folded`) replaced by `in_a_regular_bin`, summed
   out of the cells over every GC row's regular bins — production's `regular_total` — with
   `folded()` returning `in_a_regular_bin + overflowed`. `overflow_fraction` now divides one by
   the other, which is `coverage_model.rs`'s formula rather than a quantity that happens to equal
   it.
2. In `OneBinSchemeAnswer::of`, an assertion that the cells and `windows_folded` agree, whose
   message names cell saturation as the only way they can part and says which direction the
   fraction would err.
3. A test,
   `the_overflow_share_is_over_the_windows_in_the_cells_and_the_guard_is_not_reached_at_it`:
   800 windows at 1 read a position, 200 at 20, and 30 on a contig too short to clear the floor.
   It pins the share at 200 in 1,000 rather than 200 in 1,030, and pins that a sample sitting
   exactly on the guard is **taken** — the two survivors of the mutation run.

**`src/ng/window_coverage/mod.rs`** (doc comments only)

4. `DEPTH_SCALE_WINDOWS`: "nine of the ten" → "eight of the ten", with the two exceptions sized;
   "13% on average" said to be over all ten; "16 MB" → "16.8 MB", with the doubling `Vec`'s
   capacity named as the basis both it and the 262 kB are taken on.
5. `WindowCoverageConfig`'s doc and its `depth_bins` / `depth_scale_windows` /
   `depth_range_in_medians` fields: no longer "provisional", now pointing at the measurement that
   settled each, in the shape `min_window_positions` already used after D1.

**`doc/devel/ng/spec/window_coverage.md`**

6. §3.4 and §9: "nine of the ten" → "eight of the ten", and §9's "reads about 13% shallow" given
   its basis and its worst case.

**`doc/devel/reports/implementations/ng_window_coverage_d2_2026-09-07.md`**

7. The median-table sentence (Major 1), the memory sentence (Minor 3), the depth range in
   assumption 1 (Minor 6), the "780 times the headroom" sentence (Minor 7), "Five" → "Seven"
   (Minor 8) and "Six tests" → "Seven".
8. "How it is measured" now says what the fraction is and is not divided by, and that the two
   readings of the total are asserted equal.
9. The load-bearing paragraph now carries the non-monotone series (2,995 / 1,945 / 1,972 / 2,017 /
   0), which is what actually supports "no affordable prefix fixes it".
10. The tests table gains the new row; the validation line goes to 22 passed and says where the
    twenty-second came from; "Changes made" records the `WindowCoverageConfig` edit.

## 8. Validation, on the tree as I am leaving it

All three run in the container through this worktree's `scripts/dev.sh`.

- `cargo test --all-features --example ng_window_coverage_probe` — **22 passed; 0 failed; 0
  ignored**.
- `cargo check --lib --tests --all-features` — clean, `Finished dev profile`.
- `rustfmt --check --edition 2024 examples/ng_window_coverage_probe.rs
  src/ng/window_coverage/mod.rs` (rustfmt 1.9.0-stable, from the toolchain in the image) — **exit
  0, no hunks**.
- Also re-run: `cargo test --lib --all-features` — 6,378 passed, 0 failed, 15 ignored, unchanged;
  `cargo clippy --lib --all-features --bins --example ng_window_coverage_probe` — the same 3
  pre-existing `needless_lifetimes` warnings and nothing new.

**Not re-run: the five walks.** They need the `.psp` stores under the branch worktree's `tmp/`,
which is outside this worktree. Nothing I changed can move a published figure: the fix replaces a
denominator with one I verified equal to it on all ten sample-stores from the raw output, and the
new assertion fires only where they are unequal.
