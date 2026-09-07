# window coverage — D3: what this plan costs a run, per sample and in total

**Date:** 2026-09-07
**Plan step:** [window_coverage.md](../../ng/impl_plan/window_coverage.md) Milestone D, step D3
**Spec:** [window_coverage.md](../../ng/spec/window_coverage.md) §3.4, §4, §5;
[run_streaming.md](../../ng/spec/run_streaming.md) §7.2
**Branch:** `ng-window-coverage`
**Builds on:** [C5](ng_window_coverage_c5_2026-09-07.md), [D2](ng_window_coverage_d2_2026-09-07.md)

## The answer

**106.4 kB a sample for the whole pass, and 368 kB a sample while the depth axis is being
fitted.** The figure is *added up from the code* — the struct sizes and allocations a run makes,
each pinned by a library test — rather than read off a run, because a term this size is far below
what a whole-run memory measurement can resolve. The run was made anyway, and what it gives is a
ceiling: the plan costs less than 1.8 MB a sample, seventeen times looser than the figure above.

| what | bytes a sample, psp mode | when it is live |
|---|---|---|
| the histogram, 50 × 401 cells of `u32` | **80,200** | allocated whole in `new`, held to the pass's end |
| the sliding window's buffer, 512 covered positions of 16 bytes at capacity | **8,192** | the whole pass |
| the look-ahead's extra retained positions, half a window at one record a base — 72 bytes each: a record's summary and range at 48, and the window its centre finalised at 24 | **18,000** | between one eviction and the next |
| **held for the whole pass** | **106,392 — 106.4 kB** | |
| the held-back windows the axis is fitted from, 16,384 pairs of `f64` | 262,144 | until the sample's ten-thousandth window, then never again |
| **peak, while the axis is being fitted** | **368,536 — 368 kB** | |

**Against [`run_streaming.md`](../../ng/spec/run_streaming.md) §7.2's 500 kB an open sample: 21%
for the whole pass, and 74% at the peak.** On tomato, where §7.2 measures an open psp at 108 kB, a
sample goes to 214 kB and touches 477 kB while its axis is fitted — inside the budget at both. On
a human reference that open psp costs 480 kB, so the pass alone is 586 kB, which is over, and
stays over until the run's readers share one contig list — the psp path's own step D1, which takes
the 480 to 123 and so the pass to 229 kB and the peak to 492 kB. **This plan does not create that
problem and does not fix it**: §7.2 measures a human open psp at 480 kB, of which 357 kB is the
header and almost all of *that* is the reference's contig list, kept once per open sample.

**One term is unpriced and is this report's own open item.** The 24 bytes a finalised window is
charged over the *whole* stretch the merge retains, not only the half-window the look-ahead added
— and how long that stretch runs is set by the merge's eviction rhythm, which nothing here
measures. The 18,000 above is the look-ahead's marginal share of it. A high-water mark on
`WindowCoverageInProgress::finalised` would settle it; it is not built.

**Direct mode is dearer and is not priced here.** It holds the built record beside the summary —
`SampleLocusObservations` is 120 bytes inline plus its evidence on the heap — so a retained
position there costs at least 192 bytes rather than 72. Every figure in this report is psp mode's.

At scale the whole-pass term is **106 MB across a thousand samples and 319 MB across three
thousand**, against spec §3.4's own estimate of 80 MB and 240 MB — which counted the histogram
alone.

**What is asserted, and what is assumed.** Each term's per-unit size is pinned by a library test
(`one_samples_accumulator_is_88_kb_for_the_pass_and_262_kb_more_until_the_axis_is_fitted` in
`window_coverage/accumulator.rs`,
`a_retained_position_in_psp_mode_costs_a_summary_a_range_and_a_window_at_72_bytes` in
`cohort_merge/observation_cache.rs`), so a struct that grows a field fails a test rather than
quietly making this report wrong. **What no test pins is how many records a sample holds**: the
18 kB is half a window at *one record a base*, and the tomato slice's samples average 0.99. It is
also a psp-mode figure, which is the mode measured here — direct mode holds the built record
beside the summary, 120 bytes of it inline plus its evidence on the heap, so its look-ahead term
is several times larger and is not priced.

## What the whole-run measurement says, and what it cannot

**Peak resident of a psp-mode calling run on the tomato slice, before and after this plan**, over
cohort sizes 1, 2, 3, 4, 6, 10, 20, 40 and 63 with three repeats each — 27 runs an arm, 54 in all,
the same 63 stores in the same order at every size so that a step from *N* to *M* adds samples and
changes nothing else:

| | peak resident, MB |
|---|---|
| **slope** | **39.5 a sample before, 40.3 after** — a difference of 0.82, standard error 0.49 |
| **intercept** | −11.8 before, +14.5 after |
| **residual spread about the line** | 38.7 before, 34.3 after |

**The slope did not move, and could not have been seen to.** What this plan adds is 0.10 MB a
sample. What the fit resolves is a difference of 0.82 MB a sample with a standard error of 0.49
— a 95% interval running from −0.1 to +1.8, which spans zero and is eighteen times wider than
the whole effect. Two runs of the *same* binary at 63 samples differ by as much as 134 MB.

**So what the run gives is a ceiling, not a value**: it says the plan costs less than about
1.8 MB a sample, seventeen times looser than the 106.4 kB the code adds up to. Nothing above that
ceiling appeared.

The three cohort sizes the plan names:

| samples | before, MB | after, MB |
|---|---|---|
| 1 | 25.9 | 59.1 |
| 6 | 216.9 | 240.3 |
| 63 | 2,447.9 | 2,539.0 |

**The 39.5 MB a sample in both arms is not the budget's quantity** and must not be read against it.
§7.2's 500 kB is the resident state of *one open psp* — its header, block index, footer and cursor,
measured at 108 kB on tomato. The 39.5 MB is everything a calling run holds per sample, which on
this slice is dominated by the merge's held records and by the psp source's arena, which
[`cohort_merge_psp_path.md`](../../ng/spec/cohort_merge_psp_path.md) §3.4 records as never
released. **That term is not this plan's** and this plan neither worsened nor fixed it.

## The one term the measurement can see, and what is known about it

**A fixed term of tens of megabytes that follows the ground a run walks, not the cohort.** The
intercepts differ by 26.3 MB with a standard error of 12.9 — the only difference the regression
resolves at all. Whether it really is cohort-independent was settled by changing the ground
instead of the cohort, one sample either way:

| ground walked, one sample | before, MB | after, MB | difference |
|---|---|---|---|
| 10 kb (one interval) | 13.3 | 18.8 | **5.5** |
| 200 kb (two intervals) | 25.9 | 59.1 | **33.2** |

Twenty times the ground, six times the difference — **so it follows the ground**. What this table
cannot settle is whether it is *also* per sample, because both its rows are at one sample, where a
per-run term and a per-sample term are the same number. *(Five repeats at 10 kb and three at
200 kb. The two arms do not overlap at either ground — before spans 5.3 to 17.8 MB at 10 kb and
after 18.6 to 19.7 — though at 10 kb the before arm's own spread, 12.5 MB, is wider than the
5.5 MB difference between them.)*

**What rules out per-sample is the cohort sweep, and it does so decisively.** A 26 MB term charged
per sample would appear as a 26 MB difference in *slope*; the slopes differ by 0.82 with a standard
error of 0.49, which puts a per-sample reading **52 standard errors** from what was measured. Said
in the raw numbers: at 200 kb the two arms differ by 33.2 MB at one sample and by 91.1 MB at 63,
where per-sample would have put about two gigabytes between them. So it is not paid once per
sample, which is what rules out the
two terms that are: the extra records the look-ahead makes each sample hold, and the arena each
sample's psp source never releases
([`cohort_merge_psp_path.md`](../../ng/spec/cohort_merge_psp_path.md) §3.4).

**Which allocation it actually is was not measured, and the obvious candidate does not fit the
size.** The reference the cache now reads once per cover ([the step that gave the cache its
reference accessor](ng_window_coverage_c1_2026-09-06.md)) is charged to the run and grows with
ground, which is the right shape — but 200 kb of ground is 200 kB of reference bases, against a
term of 33.2 MB there, some 166 bytes for every base walked. Naming the allocation needs a heap
profile of a one-sample run at the two ground sizes. Nothing in this plan turns on the answer:
what a per-sample budget needs is that the term does not multiply by the cohort, and that is
what the two tables show.

**And the term's own size is not robust, because the instrument reads low at exactly the runs it
rests on.** `peak_rss.sh` samples the kernel's high-water mark every 20 ms, so it can only
*under*-report a peak, and the shortest runs — the one- and two-sample ones — are where a sample
is most likely to be missed. Those are also the rows an intercept is fitted from. Refitting with
the one-sample rows dropped takes the intercept difference from 2.05 standard errors to 1.66, and
dropping two and three sizes takes it to 1.59 and 0.80. **So "tens of megabytes, following the
ground" is what this measurement supports; 26.3 MB is not a figure to quote.** The one- against
63-sample comparison above does not depend on the fit and is unaffected.

## Where the numbers came from

- **The stores**: the 63 tomato accessions of `benchmarks/tomato1/crams`, walked once into psps
  over the two 100 kb intervals every step of this plan has used, so that no run re-reads a CRAM.
- **The two binaries**: this branch's head, and one built from `a33ada0f`, the commit this branch
  starts at. **The psps are the same files for both** — this plan changes no psp byte (spec §1.2),
  and both binaries produce the same VCF record counts from them at all three sizes (21, 245 and
  5,492), which is what makes the comparison a comparison.
- **Peak resident** through [`scripts/peak_rss.sh`](../../../../scripts/peak_rss.sh), which
  samples `/proc/<pid>/status`'s `VmHWM` — the kernel's own high-water mark — because the dev
  container ships no `time(1)`.
- **The harness names its binary rather than searching for the newer one**, which the repository's
  other memory scripts do: this measurement compares two builds, and "take the newer" would
  silently measure one of them twice.

## Assumptions and deviations, all minor and all recorded

1. **The per-sample cost is added up from the code, not read off the run.** The plan asks for
   both, and the run cannot resolve 100 kB against a 34 MB residual. Both are reported, and the
   run is reported as what it is: a ceiling, at 1.8 MB a sample.
2. **Nine cohort sizes rather than the three the plan names**, because three points cannot separate
   a slope from an intercept when the scatter is this wide, and the intercept turned out to be the
   only term the measurement resolves.
3. **The two footprint figures are asserted in library tests**, which is a small addition to
   `src/` in a step the plan describes as a measurement. The alternative was arithmetic in prose
   over private struct layouts, which is exactly the claim that goes stale silently.
4. **`before` was built by detaching this worktree to `a33ada0f`**, building, and returning. No
   other worktree was touched.
5. **Tomato only, and psp mode only**, where the plan's verification summary asks Milestone D's
   measurements for both benchmarks. What a human reference does to the *pass* is arithmetic
   above (§7.2's 480 kB header, this plan's 100 kB beside it) and needs no run; what it would do
   to the whole-run terms is unmeasured. Direct mode's look-ahead is dearer per held record —
   the built record beside the summary — and is unmeasured too.

## Changes made

- **Two library tests**, below. No behaviour changes: `src/` gains assertions and nothing else.
- **Spec [`window_coverage.md`](../../ng/spec/window_coverage.md) §4** — the three-thousand-sample
  bullet carries the per-sample breakdown and the total, and reconciles them with §3.4's 240 MB,
  which counted the histogram alone.
- **Spec §5** — the memory bullet names the terms that make the 106.4 kB, the transient that
  makes the peak, and the run-charged term the whole-run measurement resolved.

## Tests added

Two, both pinning a number this report quotes.

| test | what it pins |
|---|---|
| `window_coverage::accumulator::tests::one_samples_accumulator_is_88_kb_for_the_pass_and_262_kb_more_until_the_axis_is_fitted` | the histogram's 80,200 bytes, allocated whole at `new` and never resized; the sliding buffer's 8,192 and the held-back transient's 262,144, both read off a fed accumulator's own capacity rather than multiplied out; and the struct sizes they rest on |
| `cohort_merge::observation_cache::tests::a_retained_position_in_psp_mode_costs_a_summary_a_range_and_a_window_at_72_bytes` | a retained position's 72 bytes in psp mode — 48 for the record's summary and range, 24 for the window its centre finalised — and that half a window of them at one record a base is 18 kB, with the half-window taken from `WINDOW_BP` rather than retyped |

Each assertion's message names the figure it guards and the documents that quote it, so a struct
that grows a field tells its author which paragraphs to revisit.

## What the reviews changed

Two agents, each detached at the step in its own worktree; the full reports are beside this one,
[one on correctness and every number](../reviews/ng_window_coverage_d3_correctness_2026-09-07.md)
and [one on naming, structure and prose](../reviews/ng_window_coverage_d3_design_2026-09-07.md).
**1 Blocker, 9 Major, 11 Minor.** Every statistic in the run reproduced to the digit — both
slopes, both intercepts, both residual standard errors, and the standard errors of both
differences. What was wrong was almost everything around them.

- **Blocker: the answer compared the wrong quantity against the budget.** The first draft said
  100.2 kB a sample and put the 262 kB transient below the total line as if it were not charged,
  because it had the histogram's lifetime wrong — "from the depth axis being fitted". But
  `counts: vec![0; cell_count()]` runs in `new`, so all four terms are live together and the peak
  is 368 kB, 74% of the budget rather than a fifth. **The step's own test proved it and the prose
  read it backwards**: the test asserts the counts vector's length immediately after `new`.
- **A per-sample term was missing altogether.** Each sample's finalised windows are retained in
  the cache and evicted with the records — 24 bytes a retained position, on top of the 48 the
  record costs. The look-ahead's share goes from 12,000 bytes to 18,000 and the per-sample total
  from 100.2 to 106.4 kB. **The rest of that stretch is still unpriced**, and is named above as
  this report's open item.
- **Two of the four "pinned" figures were not pinned.** The sliding buffer was asserted as
  `501 × 16` and the transient as `16_384 × 16` with 16,384 a bare literal — neither ran the
  accumulator or read a capacity. A change to reserve the held-back list up front would have
  taken the transient from 262 kB to 160 kB with the test still green. Both are now read off a
  fed accumulator's own capacities, which also corrected the buffer: a `VecDeque` allocates its
  capacity, so it is 8,192 bytes and not 8,016.
- **An argument did not follow.** "So the term is not per sample" was drawn from a table whose
  two rows are both at one sample, where a per-run term and a per-sample term are the same
  number. The valid argument was a section away and is now in its place: a per-sample reading
  sits 52 standard errors from the measured slope difference. And "the difference is larger than
  either spread" was simply false — 5.5 MB against a before-arm spread of 12.5.

And a methodological point neither the plan nor I had considered: **`peak_rss.sh` samples the
kernel's high-water mark every 20 ms, so it can only read low**, and the shortest runs are where
a peak is most likely missed — which are exactly the rows an intercept is fitted from. Dropping
them takes the intercept difference from 2.05 standard errors to 0.80, so its size is not a
figure to quote, only its shape.

Smaller, and all mine: 48 bytes a held record is psp mode's and direct mode's is 152, which the
draft never said; "nineteen times wider" used a 95% interval where the ± is a standard error;
"the 480 kB header" inverted §7.2, where 480 kB is the whole open sample and 357 kB the header;
the two tests were named for a report rather than for what they pin, and eleven of their thirteen
assertions carried no failure message; and spec §3.3 still said the look-ahead costs "some tens of
kilobytes at a thousand samples" against this step's own 18 MB.

## Validation results

In the container, on this worktree:

- `cargo test --lib --all-features` — **6,380 passed**, 0 failed, 15 ignored (6,378 before this
  step; the two are its own).
- `cargo check --lib --tests --all-features` — clean.
- `cargo clippy --lib --all-features --bins --example ng_window_coverage_probe` — 3 warnings, all
  `needless_lifetimes` in `src/ng/run/cohort_merge/`, all predating this branch.
- `rustfmt --check --edition 2024` — `accumulator.rs` 0 hunks; `observation_cache.rs` 4, the same
  four that predate this branch and none of them in the region this step touched.
- **The standing oracle**: 2,311 records, sha256 `84ad19c2…0590d` on both routes; the two modes'
  windows identical over 13,866 rows, 13,589 of them a measurement; their histograms identical,
  6 samples, 6 fitted.

## Tradeoffs and follow-ups

- **A calling run holds 39.5 MB a sample on this slice, 79 times §7.2's per-open-psp budget, and
  it is not this plan's.** Both arms carry it and it was there before this branch started. The
  psp source's unreleased arena
  ([`cohort_merge_psp_path.md`](../../ng/spec/cohort_merge_psp_path.md) §3.4) is the named
  candidate; nothing here measures which term it actually is. Chasing it needs a heap profile —
  which allocation sites hold what — of a run at two or three cohort sizes on this same slice;
  whole-run peak resident, the instrument used here, cannot attribute a byte to anything. It is
  about four hundred times what this plan adds a sample, and it belongs to whoever next works on
  the psp source.
- **The 262 kB transient could be removed** by fitting the depth axis from a running quantile
  rather than a held-back list, at the cost of a different median. Nothing asks for it yet, but
  it is over half the 500 kB budget while it is live, and it is what puts a sample's peak at
  368 kB — so it is the first thing to cut if a human-reference run finds the budget binding.
- **The histogram is the term that scales**, and the bin-scheme measurement one step earlier
  ([its report](ng_window_coverage_d2_2026-09-07.md)) found nothing forcing its 400 bins.
  Halving them halves 80.2 kB to 40.2 and doubles the bin width; if the per-sample budget ever
  binds, that is the knob, and the filter's own validation is what says whether the resolution
  can be spared.
