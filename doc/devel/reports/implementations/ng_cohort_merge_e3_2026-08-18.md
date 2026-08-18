# ng cohort merge — E3: the builders run at the same time

*Implementation report, 2026-08-18. Step E3 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §6.1, §6.2 and §6.4, with §15's oracle as the thing to
reproduce. **§6.2 was amended by the owner as part of this step** — see §2.*

> **This is the first draft, and [the review](../reviews/ng_cohort_merge_e3_2026-08-18.md)
> changed what the tests can see.** The concurrency held under attack — 400 random layouts, eight
> pool sizes, 200 repeats, no disagreement — but **the round itself was invisible to the suite**:
> its size was pinned by an inequality a doubled round also passes, and `.par_iter()` could be
> replaced by `.iter()` with everything still green. One number below is wrong too: a short last
> round costs 2 idle builders on the fixture named, not 15. What landed is in
> [the fix report](ng_cohort_merge_e3_fixes_2026-08-18.md).

## 1. Plan

Hand the building regions out to several builders at once, reading the shared observation
cache, with the organiser collecting — and change nothing about the answer.

## 2. The design question this step had to settle, and the owner's ruling

`ObservationCache::with_observations` takes `&self` while `cover` and `evict_before` take
`&mut self`, so **the borrow checker forbids drawing the readers forward while any builder holds
a window.** That is the right refusal on one thread; several builders reading while the
organiser advances needs a shape, and the spec named none.

Worse, **the spec contradicted itself**. §6.2 said each region's task owns k cursors of its own
and explicitly *rejected* a shared window — "it re-creates the global frontier this design
exists to remove: every builder's progress would couple through the window's trailing edge."
§6.4, decided a day later, chose the shared window anyway. D1 and D2 built §6.4's shape, and the
owner accepted it at Checkpoint D, so the code had already chosen; §6.2's objection had been
overwritten rather than answered.

**The owner's ruling, 2026-08-18: build the round, and amend §6.2 to say the coupling was
accepted and why.** §6.2 now records that the reversal is deliberate, gives the reason §6.2 had
not weighed — a reader per builder per sample is `n × k` open cursors, 48,000 at 3,000 samples
and 16 builders against 3,000 — states the coupling's exact shape and bound, and leaves the two
alternatives (an `RwLock`, or windows handed out as owned copies) named and deferred until the
round's tail has been measured.

## 3. What it is

[`parallel.rs`](../../../../src/ng/run/cohort_merge/parallel.rs), beside the two serial drivers
it must reproduce. Per round:

1. evict at the round's first region's first base;
2. cover every region of the round, on the driver's own thread;
3. run the round's builders concurrently over `&ObservationCache`, through rayon;
4. submit the outcomes to the organiser in region order, and drain what that releases.

**`builders` is how many regions are in flight, not how many threads run them.** Threads are
rayon's; what this number sets is the ground the cache must hold —
`builders × cohort_locus_builder_regions_len` bases, 320 at 16 builders on 20-base regions
(spec §6.4).

**A round never crosses an analysed region.** The run's intervals are disjoint and may sit on
different contigs, and a round opens by evicting at its first region's first base — a round
spanning two intervals would evict on the strength of ground the second does not contain. The
cost is a short last round per interval.

**The failed spans are gathered by the driver and the organiser's count is the cross-check.**
The organiser resolves overlaps and counts the failures that survived, but does not hand the
spans back (arch §4: "failed spans participate in that resolution and are never released"), so
the driver keeps the spans of every outcome it submitted and asserts the two agree.

## 4. Deviations, and what is still owed

**The organiser does not own the cache**, which arch §4 gives it. The cache is generic over its
source's error type because the run's `ObservationSource` and `RunError` do not exist yet, and
making `Organiser` generic over the same parameters to hold it would push that genericity into
a type that has no other use for it. The driver owns both instead, exactly as
`merge_cohort_through_cache` owns the cache today. Arch §4's `Organiser::cache()` is still owed,
and is a change to make when the run's own types land rather than before.

**`builders` is a bare `NonZeroUsize`, not a newtype with a default** like the module's other
three run parameters. It is not a psp-recorded or spec-defaulted value yet — spec §6.4 names the
region width as the command-line parameter and prices the builder count only as `n`. Recorded as
a choice rather than made silently.

**Two assertions in the driver are safety nets no test reaches**, and the doc says so: both can
only fire if the organiser displaced something, and nothing a builder produces can make it.
Mutation testing agrees from the other side — disabling either leaves the whole suite green,
because every defect that would trip them is already caught by the output comparison.

**One tidy taken while here:** the `member` fixture moved from `serial.rs`'s test module into
the module-wide `#[cfg(test)] mod fixtures`, because three test modules now want it and a
fixture that differs between the files comparing their outputs is a difference nobody would look
for. That is one of the two items `mod.rs` recorded as owed.

## 5. Tests

Eleven, in `parallel.rs`. The module went from 209 tests to **219**
(`cargo test --lib ng::run::cohort_merge`).

| test | what it pins |
|---|---|
| `one_builder_gives_the_serial_drivers_answer` | the base case: a round of one is the cached driver |
| `the_answer_is_the_same_at_every_builder_count` | 1, 2, 4, 8, 16 against the oracle |
| `the_parallel_driver_matches_the_cached_one_at_every_width` | five widths × five builder counts |
| `the_cache_holds_a_whole_round_and_not_one_region` | **the structural pin** — see below |
| `a_round_stops_at_the_end_of_its_analysed_region` | two intervals on two contigs, a round of 8 over 4 regions |
| `the_failed_spans_come_through_the_parallel_driver_unchanged` | a refused locus at three builder counts |
| `a_failing_source_ends_the_parallel_merge` | the source's own error, untouched |
| `merging_no_analysed_regions_in_parallel_yields_nothing` | the empty run |
| `merging_no_samples_in_parallel_yields_nothing` | the empty cohort |
| `an_analysed_region_with_inverted_ends_is_refused_in_parallel` | the guard both serial drivers keep |

**The structural pin is the one that matters**, because every output here is identical by
construction: a driver that ignored the builder count entirely would pass every comparison
above. `the_cache_holds_a_whole_round_and_not_one_region` reads what the cache still holds when
the merge returns, on a record every ten bases at 20-base regions: **2 records at one builder,
4 at two and at four, 12 at eight, 29 at sixteen.** A driver evicting per region rather than per
round would hold the same 2 at every count. The assertion is the inequality rather than the
table, because the exact records depend on where the last round's ground falls.

**Eleven mutations, nine killed and two proved redundant** (`tmp/mutate_e3.sh`): the round
forced to one region, and to the whole analysed region — both killed only by the structural pin;
eviction at the round's *last* region; only the round's first region covered; outcomes submitted
in reverse; failed spans not gathered; the drain never taken; rounds allowed to cross analysed
regions; the malformed-region check dropped. The two that survive are the driver's own end
assertions, which §4 explains.

## 6. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `219 passed; 0 failed`.
- `cargo test --lib` — `3842 passed; 0 failed; 11 ignored` at this draft; `3847` after the
  review's fixes.
- `tmp/mutate_e3.sh` — 11 mutations, 9 killed, 2 recorded as safety nets.

## 7. Follow-ups

- **The round's tail is unmeasured.** Spec §6.2 now says the two alternatives should not be
  reached for before it is, and nothing yet measures it — that wants a real cohort and belongs
  with the sweep spec §14 question 1 already names.
- **Arch §4's `Organiser::cache()`**, owed until the run's own source and error types land.
