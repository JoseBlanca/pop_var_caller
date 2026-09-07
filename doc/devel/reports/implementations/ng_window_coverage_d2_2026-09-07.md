# window coverage — D2: where the depth axis is cut, and what falls off the end of it

**Date:** 2026-09-07
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone D, step D2
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.4, §9
**Branch:** `ng-window-coverage`
**Builds on:** [C5](ng_window_coverage_c5_2026-09-07.md), [D1](ng_window_coverage_d1_2026-09-07.md)

## The answer

**All three constants stand — 400 bins spanning ten times the sample's median, fitted from its
first 10,000 windows — and the reason is that at the shipped range the worst sample-store overflows
780 times less than the fit tolerates.** Measured on the same five stores the floor measurement
walked, ten sample-stores in all.

The model fit that reads this histogram rejects a sample outright once more than **a fifth** of its
windows are in the overflow column
([`DEFAULT_MAX_OVERFLOW_FRACTION`](../../../../src/paralog/coverage_model.rs)), because at that
point the sample's own single-copy peak has left the range and the regular bins hold noise. At the
shipped setting, **nine of the ten sample-stores overflow nothing at all**, and the tenth overflows
**1,972 windows of 7,666,421 — 2.6 in every 10,000**, against a guard that fires at 2,000 in
10,000.

## What the factor of ten buys, and what narrowing it would cost

The worst overflow across all ten sample-stores, at each setting of the range:

| range, in medians | top of the axis, on the 5.1× store | worst overflow anywhere | against the fit's guard |
|---|---|---|---|
| 2.5 | 10.1 | **1,103 in 10,000** | within a factor of **1.8** |
| 5 | 20.2 | 35 in 10,000 | 57 times under |
| **10 (shipped)** | **40.3** | **2.6 in 10,000** | **780 times under** |
| 20 | 80.6 | 0 | — |
| 40 | 161.3 | 0 | — |

**No setting in that table made the fit reject a sample, 2.5 included**, so what the table ranks is
margin and not pass-or-fail. Ten was kept for the size of its margin: halving to 5 would halve the
bin width — twice the resolution — and cut the margin under the guard from 780 times to 57, and
nothing measured asks for that resolution. Widening past 20 buys nothing the guard can see and
costs resolution everywhere. The upper end of the table is where the settings tried stop, not where
a cost appears.

**2.5 is the setting that would start to matter.** Three of the ten sample-stores overflow anything
at all there, and all three are tomato: the one-accession store at 1,103 in 10,000, one accession
of the six-sample slice at 614, and a second accession of the same slice at 191. Their fitted
medians are 6.41, 6.41 and 3.97 reads a position — the shallow end of the ten — but not in that
order, so what these rows show is that the tomato stores overflow first, not that overflow falls
off with depth.

## The scale sample stands, and it is not because it is representative

**The median fitted from the first 10,000 windows is below the whole store's median on eight of the
ten sample-stores** — by 34% at worst and, averaged over all ten, 13%. The two that come out above
do so by 2% and 8%:

| store | sample | median from the first 10,000 windows | median from every window | ratio |
|---|---|---|---|---|
| tomato, 80 × 100 kb | SRS3394712 | 6.41 | 9.77 | **0.66** |
| HG002, GIAB high-confidence | HG002 | 246.21 | 302.66 | 0.81 |
| tomato slice | SRS3394712 | 6.41 | 7.85 | 0.82 |
| tomato slice | SRS3394713 | 10.12 | 12.32 | 0.82 |
| HG002, tandem-repeat tiers, 5× | HG002 | 4.03 | 4.93 | 0.82 |
| HG002, tandem-repeat tiers, 30× | HG002 | 25.12 | 30.15 | 0.83 |
| tomato slice | SRS3394711 | 13.59 | 14.60 | 0.93 |
| tomato slice | SRS3394712_SRR7279484 | 3.97 | 4.24 | 0.94 |
| tomato slice | SRS3394606 | 27.61 | 27.02 | 1.02 |
| tomato slice | SRS3394714 | 21.62 | 20.02 | 1.08 |

**It is not that 10,000 windows is too few.** On the one-accession tomato store the fitted median
barely moves between 100 windows and 100,000 — 5.25, 6.50, 6.41, 6.25 — and then jumps to 9.43 at
a million. A hundred thousand windows is 100 kb, which is one of that store's eighty intervals; a
million is ten of them. **The scale sample is a prefix of the run's own ground, and the first
interval is not the genome.** No affordable prefix fixes that: holding a million windows back
before folding any is **16.8 MB a sample — 16.8 GB across a thousand** — against the 262 kB at
10,000 that spec §3.4 already charges to the per-sample budget. Both figures are the capacity a
doubling `Vec` of `f64` pairs reaches (1,048,576 and 16,384 entries at 16 bytes), which is the
basis spec §3.4 states its own number on.

**The effect is not monotone in the length of the prefix**, which is what makes "a longer prefix"
the wrong fix rather than merely an expensive one. On the worst store the overflow at the shipped
range runs 2,995 windows at a scale sample of 100, 1,945 at 1,000, 1,972 at 10,000, 2,017 at
100,000, and only then 0 at a million. What the 2,017 at 100,000 says is that a prefix ten times
longer than the shipped one is no better; a hundred times longer is, and costs 16.8 MB a sample.

**What makes the bias harmless is the range, and the two are not independent.** The entire overflow
of the worst store — those 1,972 windows — is attributable to the short scale sample: with the
width fitted from every window instead, that store overflows nothing. A prefix 34% shallow makes
the axis 34% short, and at a range of ten what that costs is 2.6 windows in 10,000 over the top,
which the guard does not notice. That cost is what a narrower range would multiply: at a range of
5 the same prefix on the same store puts 35 in 10,000 over. **If the range is ever narrowed, this
is the term that decides how far it can go**, not the overflow measured at ten.

## 400 bins stands, and what it is for

The medians measured here span **3.97 to 246.21 reads a position — 62-fold**. A fixed bin cannot
serve both ends of that: production's 0.5× bin puts the shallowest sample's median at its bin 7
and the deepest at bin 492 of 2,000, so the shallow sample gets seven bins to resolve everything
below its single-copy peak and the deep one spends three quarters of the axis above a carrier it
will never see. Fitting the width to the sample puts every median at bin 40 of 400 **by
construction** — that is the formula, not a measurement — and the measurement is the 62-fold spread
that makes it worth doing.

Nothing here forces 400 either way. What it decides is memory — 50 × 401 × 4 bytes, 80.2 kB a
sample — and that is plan step D3's to price.

## What this step does not measure

**The cost of each setting, not its benefit**, exactly as D1. Whether a coarser or finer axis makes
the coverage model's fit better or worse is a question about the fit, and the fit is built on the
hidden-duplication filter's own branch. What is settled here is that at the shipped setting no
sample of any store on this machine comes close to the one failure the fit can see for itself, and
that the margin is large enough to absorb a scale-sample bias of a third.

## How it is measured

**The same device the floor measurement used**: one shipped `WindowCoverageAccumulator` per
candidate setting, fed exactly the positions the whole-store walk is already feeding one,
differing from it in a single field. Ten arms — the run's own configuration, five scale samples from 100 windows to
every window the store has, and four ranges from 2.5 to 40 medians. Each reports the axis it fitted
and how much of the sample fell off the end of it. Nothing in the probe recomputes a histogram.

The overflow fraction is computed **the way `coverage_model.rs` computes it** — the overflow column
summed over GC rows, over the regular bins plus that column — so the number reported is the one
that would be compared against the guard rather than a near relative of it.

**Three things are asserted before any of it is printed:**

- **the arm configured exactly as the run reproduces the walk's own histogram** — the axis read off
  each, compared as the row it prints, so the tolerance is the row's own precision. It is the one
  place this sweep and the recomputation answer the same question, and unlike D1's tie to the walk
  it is never a comparison of two zeroes, because every store here fits a real axis;
- an arm running the run's own configuration exists at all: the check above finds that arm and
  fails if there is none, rather than skipping silently;
- some arm folded a window, so a store that folded nothing is refused rather than printing ten
  rows of zeroes.

## Changes made

- **[`examples/ng_window_coverage_probe.rs`](../../../../examples/ng_window_coverage_probe.rs)** —
  `BinSchemeSweep` and its arms, the answer each one hands back, the `--bin-scheme` flag, and the
  per-sample rows. Six tests.
- **[`src/ng/window_coverage/mod.rs`](../../../../src/ng/window_coverage/mod.rs)** —
  `DEPTH_BINS`, `DEPTH_SCALE_WINDOWS` and `DEPTH_RANGE_IN_MEDIANS` keep their values and stop being
  marked soft; each doc carries what the measurement found.
- **[`doc/devel/ng/spec/window_coverage.md`](../../ng/spec/window_coverage.md)** — §3.4's "soft,
  and marked so" and §9's second OPEN closed with the measurement.

Nothing under `src/sample_summary/`, `src/paralog/` or `src/var_calling/` is touched, and no
library code changed behaviour: the only edits under `src/` are doc comments.

## Assumptions and deviations, all minor and all recorded

1. **The plan says "on both benchmarks"; the same five stores D1 walked were used**, which is ten
   sample-stores across two species, three interval lengths and — as median window depths over
   every window of each store — **4.24 to 302.66 reads a position**. The reason is D1's: the
   six-accession slice alone settles nothing, and the two benchmarks named turn out to be four.
   *(The medians the shipped setting actually fits, from each store's first 10,000 windows, run
   3.97 to 246.21; the two pairs are different quantities and the 62-fold spread quoted above is
   the second.)*
2. **The fit's guard is read out of production rather than copied.** The probe's
   `THE_OVERFLOW_SHARE_THE_FIT_ALLOWS` is bound to
   `paralog::coverage_model::DEFAULT_MAX_OVERFLOW_FRACTION`, which is a `pub` constant of a `pub`
   module. The local name is there to say in the probe's own words what the number is; the value
   cannot drift from the one that would judge the sample. Reading a public constant is not an edit
   to frozen production.
3. **The median each arm fitted is recovered from its width**, since the accumulator keeps the
   median nowhere: `median = width × depth_bins / depth_range_in_medians`, using the arm's own
   range and not the run's.

## Tests added

Seven, in the probe — 22 pass where 15 did before.

| test | what it pins |
|---|---|
| `a_stream_of_one_depth_fits_an_axis_this_file_can_recover_that_depth_from` | 600 windows all at 8 reads: median 8, width 0.2, axis to 80, nothing over — the width the accumulator fits and the median this file recovers from it are two directions of one formula |
| `the_range_decides_whether_the_deep_windows_overflow_and_the_fit_refuses_the_sample` | 700 windows at 1 and 300 at 20, on two contigs so neither mixes: at ranges of 2.5 and 10 the axis stops below 20 and the fit rejects the sample; at 40 it does not |
| `a_scale_sample_that_closes_before_the_depth_changes_costs_the_sample` | 10,300 windows at 1 then 19,700 at 100: the run's setting fits from a prefix that is all at 1, overflows 19,700, and the fit refuses the sample; holding every window back fits 100 and overflows nothing |
| `an_arm_that_holds_every_window_back_answers_as_the_run_does_on_a_short_stream` | the control for the above — under 10,000 windows both fit at `finish` from everything they have, so the difference there is the scale sample and not a mis-wired arm |
| `a_bin_scheme_sweep_disagreeing_with_the_walk_is_refused` | the arm configured as the run must reproduce the walk's histogram; closed against a walk that saw twice the depth, it stops |
| `a_bin_scheme_sweep_over_an_empty_store_is_refused` | a store that folded nothing is refused, not printed as ten rows of zeroes |
| `the_overflow_share_is_over_the_windows_in_the_cells_and_the_guard_is_not_reached_at_it` | 800 windows at 1, 200 at 20, and 30 the floor silences: the share is exactly the fit's guard, and the sample is taken — the denominator is the cells and the comparison is a strict `>`, the two places this file could disagree with the fit while every other test stayed green |

## What the reviews changed

Two agents, each detached at the step in its own worktree; the full reports are beside this one,
[one on correctness and every number](../reviews/ng_window_coverage_d2_correctness_2026-09-07.md)
and [one on naming, structure and prose](../reviews/ng_window_coverage_d2_design_2026-09-07.md).
**9 Major, 13 Minor**, no Blocker, all applied but one recorded below. **Neither agent could break
the measurement**: the median recovery inverts the fit to one part in 10^16 over 400,000 random
medians, the sweep is fed exactly the walk's positions, and all six original tests' exact counts
and medians re-derive correctly including at contig edges.

Four that were defects rather than clarity:

- **"Nine of the ten" is eight of the ten**, and this report's own table printed the two
  exceptions (ratios 1.02 and 1.08) directly beneath the sentence that denied them. **Both agents
  found it independently.** Eight of ten still carries the argument that the prefix reads shallow;
  nine overstated a systematic effect that two stores contradict, and those two are the
  shallow-interval tomato accessions where a reader would most want to know the bias can run the
  other way. The 13% average is right and is over all ten.
- **The overflow fraction was not computed the way this report said it was.** The prose claimed
  production's denominator — the regular bins plus the overflow column, out of the cells — and the
  code divided by the accumulator's `windows_folded` counter. **The two agree on all ten
  sample-stores**, which was checked, so no figure here changes; what was wrong was that nothing
  checked. Two deliberate denominators survived all 21 tests, and one of them — folding the
  floor-silenced windows in — would understate overflow by 3.6% on the tandem-repeat store, where
  212,850 of 5,910,300 windows are silenced. Both totals now come out of the cells, with the
  identity against the counter asserted and cell saturation named as the one way it can fail.
- **The fit's guard is a strict `>` and nothing pinned it.** Flipping the probe's comparison to
  `>=` survived all 21 tests: the arms the tests exercise sit at 0.0, 0.3 and 0.657, nowhere near
  the boundary. The added test puts a sample exactly on the guard and asserts it is taken.
- **Three of D1's four prose rulings were re-broken in the same forms** — the module doc still
  said "three questions" after this step added a fourth; "The fourth question" still called the
  three constants provisional in the commit that measured them; and "780 times the headroom"
  carried a ratio of overflow *shares* in the language of range *width*, which is D1's "6 in 100"
  defect one step later.

And the smaller ones, all mine: "16 MB at a million windows" was a length compared against a
capacity-based 262 kB in the same breath (16.8 MB on one basis); `WindowCoverageConfig`'s own doc
still called all four constants provisional after D1 and D2 settled them; "4.0 to 302.7 reads a
position" mixed a prefix median with a whole-store one; and "the safe band is 5 to 40" read as
measured at both ends when no setting tried, 2.5 included, actually made the fit reject a sample.

**One recorded rather than fixed**: the arm that holds every window back is about 134 MB of
transient on the 7.67-million-window tomato store, plus 16.8 MB for the million-window arm. That
is the arm's whole purpose and the probe is not production, but a later reader pointing this sweep
at a whole human genome should know what the last two arms cost.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — 6,378 passed, 0 failed, 15 ignored, unchanged: this step
  adds no library test and changes no library code but doc comments.
- `cargo test --all-features --example ng_window_coverage_probe` — **22 passed**, 0 failed
  (15 before this step).
- `cargo check --lib --tests --all-features` — clean.
- `cargo clippy --lib --all-features --bins --example ng_window_coverage_probe` — 3 warnings, all
  `needless_lifetimes` in `src/ng/run/cohort_merge/`, all predating this branch.
- `rustfmt --check` — both files this step touched at 0 hunks.
- **The standing oracle**: 2,311 records, sha256 `84ad19c2…0590d` on both routes; the two modes'
  windows identical over 13,866 rows, 13,589 of them a measurement; their histograms identical,
  6 samples, 6 fitted.
- **The five walks re-run on the tree as it ships**, and every figure quoted above reproduces.

## Tradeoffs and follow-ups

- **The range and the scale sample trade against each other**, and the report above says by how
  much: a prefix a third shallow costs a third of the axis, and only the range's margin makes that
  free. A later step that narrows the range has to re-check the prefix bias first.
- **The prefix bias has a fix that is not a longer prefix** — sampling the scale windows across
  the run's ground rather than taking the first of them. It is not needed at the shipped range and
  is not built.
- **The bin count is a memory knob**, priced at plan step D3.
