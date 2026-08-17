# ng cohort merge — D1: the observation cache

*Implementation report, 2026-08-17. Step D1 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §6.4 and [arch](../../ng/arch/cohort_merge.md) §4.
Revised after [the D1 review](../reviews/ng_cohort_merge_d1_2026-08-17.md), which found four
wrong claims in the first draft of this report — each is corrected in place and named in §8.*

## 1. Plan

One forward reader per sample, and a window over what the builders currently span: `cover`,
`with_observations`, `evict_before`. Nothing reads the cache yet — D2 is where the serial
driver goes through it.

## 2. Why it exists, in one number

A builder closes loci from the beginning of whatever observations it is handed and discards
those that opened before its own ground (`build_region`), so it pays for its prefix: **3.3 µs
per prefix base at 63 samples, 40 µs at 250** — the C1 review, release build. The same effect
end to end, measured **on one sample** by the C2 review: **20,000 observations cost 5.4 ms
merged as one analysed region and 184 ms as a thousand.** (The two figures come from different
cohorts and are quoted with their own; the second is not a 63-sample number.) Short building
regions are what the parallel arrangement needs and what the memory bound rests on
(spec §6.4), and they are only affordable if each builder is handed a window over its own
ground. That window is this file.

## 3. What it is

New file [`organise.rs`](../../../../src/ng/run/cohort_merge/organise.rs) —
`ObservationCache`, one `SampleWindow` per sample, each a source and the observations drawn
from it that have not been evicted. The organiser proper (ordered release, overlap resolution)
joins it at milestone E, which is why the file is named for the whole job; keeping the cache
there is what will let `cover` and `evict_before` become **private to the file** when the
organiser arrives, so that `build.rs` cannot reach them however it holds the cache.

- **`cover(region)`** draws every sample forward until the region is covered *and* until the
  chain a locus opening inside it can reach has ended.
- **`with_observations(span, f)`** hands each sample's window out for the length of the call —
  `&self`, so a builder cannot disturb it.
- **`evict_before(position)`** drops what ends before `position`.

## 4. The one piece of real algorithm: the chain's reach is a fixpoint across samples

The reach starts at the region's last base and grows with every observation that begins at or
before it — the same `<=` the closer chains on (spec §4.1). **It cannot be computed in one
sweep of the samples**, because one sample's deletion is what makes another sample's later
observation part of the locus, and that other sample may already have been swept. So the
samples are swept repeatedly until a whole sweep moves the reach no further.

**The sweep count is a property of the data.** Two sweeps is ordinary; three samples can need
three (`a_chain_that_needs_a_third_sweep_is_drawn_whole`), and the worst case is one sweep per
sample, when a chain runs through the cohort in decreasing sample order.

**Each sweep re-reads a sample's window from its first observation rather than resuming where
the last sweep stopped.** Re-reading one already inside the reach cannot move it, so the scan
is idempotent — and keeping no resume mark means there is no mark for eviction to correct or
for the next region to inherit: a cover is a function of the window, the source and the region
asked for. The review made the alternative and measured what it costs: a mark carried across
covers skips a survivor of an eviction (drawing one observation too far), and points **past
the end** of a window that lost two entries at once. Both cases are now tests.

## 5. Deviations from the architecture's sketch, and why

Three, all recorded at the code:

- **The source's error type is the source's, not `RunError`.** Arch §4 writes
  `cover(&mut self, region) -> Result<(), RunError>` and names `run_streaming.md` arch §2's
  `ObservationSource` as what the cache is handed. **Neither exists in the code** — both belong
  to the run's own document, which this plan puts out of scope — so the cache is generic over
  an iterator of `Result<SampleLocusObservations, E>` and passes failures through untouched.
  When `ObservationSource` lands, `observations_in` is exactly this shape and `RunError` is
  exactly this `E`. The cache's doc now also states the requirement the passthrough leaves on
  the source: **`E` must name its own sample**, because this cache knows which reader failed
  and adds nothing.
- **`with_observations`'s `span` bounds the left edge and checks the right, rather than
  trimming it.** The arch says "every sample's observations overlapping `span`"; trimming the
  right would cut a locus that opens inside the span and reaches past it, which is the deletion
  the whole ownership rule exists to keep whole (spec §6.1). So the start selects the left edge
  and the **end is checked against the ground `cover` actually reached** — which is what keeps
  the parameter honest: three review categories independently found that an ignored `span.end`
  was a field no behaviour depended on.
- **`evict_before` and the left trim share one predicate** — keep an observation whose last
  base is at or after the position — rather than "starts after". That is what makes both safe:
  the observation that chains a locus across the edge is exactly the one that reaches over it,
  and it is kept. A locus whose first position is at or after the edge has every member at or
  after it too, so nothing it needs is ever dropped.

## 6. What it costs to be forward-only

**The window overshoots by at most one observation per sample.** The only way to learn that the
next observation begins beyond the reach is to draw it, and once drawn it is held rather than
thrown away. `drawing_stops_one_observation_past_the_reach` pins the count exactly (two draws,
not one and not three) so the overshoot is a stated property rather than an accident.

**The overshoot is also what makes the fixpoint hard to test**, and it is where the first draft
of this report was wrong: an observation one place beyond the reach is in the window whether or
not the sweep that would have folded it ran. A fixture that shows the difference needs a
*second* observation behind the first.

## 7. What a cover costs, measured

The D1 review measured the loop on a synthetic cohort in a release build, over 20-base regions:

| what | cost |
|---|---|
| one whole cover at 3,000 samples | 2.87 ms |
| one extra sweep at 3,000 samples | 3.1 µs |
| the worst case (a chain running backwards through the cohort), one 11-base region at 3,000 samples | 28 ms |
| a cover at 1,000 samples with 4 observations held per sample | 616 µs |
| the same with 200 held per sample | 1,028 µs |

The last two are why the doc no longer says the window is "short by construction" without
naming the construction: it stays short **only while the organiser evicts at the pace it
releases ground**, and the organiser is milestone E. The `held` term is what the re-read
multiplies.

## 8. Tests — 27, and the mutations they kill

Twenty-seven tests in `organise.rs`, and the module's suite is 142. Sixteen mutations were
written against the source and the module's tests re-run in the container for each; one of the
sixteen does not compile, so **fifteen were run and fourteen fail at least one test.** The
driver is `tmp/mutate_d1.sh` (scratch, not tracked).

| mutation | killed by |
|---|---|
| sweep the samples once instead of to a fixpoint | 4 tests, including `the_chain_reach_follows_a_widening_in_a_later_sample` |
| stop the fixpoint after two sweeps | `a_chain_that_needs_a_third_sweep_is_drawn_whole` and the differential |
| stop drawing *at* the reach rather than past it | `an_observation_beginning_on_the_reach_is_drawn_and_widens_it` and the differential |
| trim the window at observations reaching *past* the span rather than *to* it | 3 tests |
| evict by where an observation starts rather than where it ends | 5 tests |
| evict from the first sample only | `eviction_drops_from_every_sample_not_only_the_first` |
| hand out the whole window instead of trimming | 4 tests |
| start the reach at the region's first base | 21 tests |
| never grow the reach | 6 tests |
| drop the coordinate-order check | 2 tests |
| compare coordinate order only *within* a contig | `a_source_that_goes_back_a_contig_is_refused` |
| never check the window against the covered ground | `a_window_over_ground_no_cover_reached_is_refused` |
| check the covered ground against the span's *start* rather than its end | `a_window_reaching_past_the_covered_ground_is_refused` |
| read an inverted span's left edge as its `start` field | `a_region_given_end_first_still_covers_its_ground` |
| drop the inverted-region defence in `cover` | the same test |

**The one survivor is a survivor by design.** Replacing eviction's prefix drain with a filter
over the whole window passes every test *and* the differential. It cannot be otherwise: a
sample's records are disjoint and ascending (`build_region` asserts it), so reach is monotone
across one sample's window and the first survivor is the last non-survivor's successor. The
prefix form is chosen for its cost, not for a difference in what it keeps, and the doc now says
so rather than promising a property no fixture can reach.

**The strongest test is the one the review wrote**:
`a_builder_fed_from_the_cache_closes_the_loci_a_whole_stretch_would` drives a real
`build_region` through the cache over 200 seeded random layouts — 1 to 8 samples, 1 or 2
contigs, one observation in ten a deletion up to 150 bases wide, building regions 1 to 12 bases,
with gaps in the ground — and compares each region's whole `RegionOutcome` against the same
builder handed the entire stretch. It agrees everywhere as the code stands (600 of 600 layouts
in the reviewer's run, 200 of 200 here), and disagrees on 410 of 600 under a cover capped at two
sweeps. That is the property this file exists for, and until the review nothing asserted it.

**Four claims in this report's first draft were wrong, all of them about my own fixtures**, and
each is corrected above: that the two-sample widening fixture pinned the fixpoint (it did not —
the overshoot keeps its far observation either way); that the same fixture killed the
single-sweep and never-grow mutations (the boundary test killed both); that the boundary case
needed three samples to be visible (two suffice, if the second sample carries a second
observation); and a mutation table of eleven rows under a sentence claiming twelve. The pattern
is the one the plan's own skill names: figures quoted from the design documents were right, and
every wrong one was my own claim about my own test.

## 9. Validation

`cargo fmt --check` clean; `cargo clippy --lib --all-features -- -D warnings` clean;
`cargo test --lib ng::run::cohort_merge` → `142 passed; 0 failed`. `--all-targets` clippy is red
on this branch for 49 pre-existing reasons in other modules — the standing item.
