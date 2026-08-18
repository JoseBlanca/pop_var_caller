# ng cohort merge — E4: asserting the milestone

*Implementation report, 2026-08-18. Step E4 of
[the plan](../../ng/impl_plan/cohort_merge.md); design authority
[spec](../../ng/spec/cohort_merge.md) §15 — the partition-invariance oracle.
Revised after [the E4 review](../reviews/ng_cohort_merge_e4_2026-08-18.md), which found the
sweep's own fixture-shape assertions unable to fail; the fixes are in place and named in §7.*

## 1. Plan

State milestone E's claim once, where a reader looking for it will find it: **the same cohort
observations and the same failed spans as the oracle, at 1, 2, 4, 8 and 16 builders and at
several building-region widths.**

## 2. What the earlier steps already covered, and what was missing

E3 landed two sweeps of its own, and both compare against the **cached** serial driver rather
than the oracle. That is sound — D2 proved the cached driver byte-identical to the oracle — but
it makes the milestone's claim transitive where the plan asks for it directly, and neither sweep
runs on a fixture carrying a refused locus.

Two things were missing, and E4 adds both.

**The cross-product against the oracle itself**, in
[`parallel.rs`](../../../../src/ng/run/cohort_merge/parallel.rs):
`the_parallel_merge_is_the_oracles_at_every_width_and_count` — five counts of regions in flight
(1, 2, 4, 8, 16) by five widths (1, 3, 20, 47, 600), twenty-five comparisons of the whole
rendering, which is the observations and the failed spans together.

**The fixture has to reach both shapes or the sweep proves nothing**, so it carries a deletion
at 305–330, which chains two samples into one locus, and a 91-base record the 50-base default
bound refuses, which puts a failed span across many regions at the narrow widths. **Both are
asserted present in the oracle before the sweep runs**, so the twenty-five comparisons cannot
pass by comparing nothing. The deletion reaches across a region boundary at **four of the five
widths** — at 1, 3, 20 and 47 the region holding base 305 ends before 330; at 600 the whole
stretch is one region and there is no boundary. That row earns its place differently: it is
where one round stands fifteen-sixteenths idle at 16 in flight. At the other end, the width of 1
makes 600 regions, so every round but the last is full — at 16 in flight the last holds 8.

**The parallel driver over the whole existing fixture corpus.** The claim is worth more on
ground nobody wrote for it, so `serial.rs`'s shared driver-agreement helper now also merges in
parallel and refuses any difference from the oracle — at 1, 4 and 16 regions in flight, since
its busiest caller runs it two hundred times and a full sweep there would not pay. That reaches **the two hundred random layouts** and,
more valuable, **the locus built from observations the generic generator actually minted**: two
samples' reads on disk, a substitution and a five-base deletion, chained into one cohort locus.
Nothing in `parallel.rs`'s own fixtures is real data.

## 3. That the wiring is discriminating, measured

Breaking the parallel driver — evicting at the round's *last* region instead of its first, one
line — now fails **eleven** tests, and three of them are in `serial.rs`:
`several_contigs_come_through_the_cache_unchanged`, `the_two_drivers_agree_on_random_layouts`
and `the_cache_changes_nothing_at_every_building_region_width`. Before E4 the same break failed
only `parallel.rs`'s own tests.

## 4. Tests

Three changes, two tests added. The module went from 224 tests to **226**; most of the extra
coverage is in the count of comparisons rather than the count of tests, which is the point of
routing the parallel driver through the shared helper.

**The sweep now asserts its own division.** The review's Major: the two fixture-shape assertions
are made against the oracle, which divides nothing, so neither could fail if the width list were
edited to widths that stop dividing the fixture — the sweep would then compare five undivided
merges against an undivided oracle and report the milestone proved. Each width now asserts how
many building regions the deletion touches — 26, 9, 2, 2 and 1 at widths 1, 3, 20, 47 and 600 —
straight from the divider the driver uses. The oracle's observation count is pinned at 50 for
the same reason: the fixture is shared by three files, and an edit dropping most of its loci
would otherwise leave the sweep agreeing on a nearly empty answer.

**Sixteen regions in flight joined the helper's sample**, because the review measured what 1 and
4 miss: a defect confined to a round's fifth region or later fails five tests and none of them is
in `serial.rs`, since the helper's rounds never held more than four.

## 5. What the review corrected

Besides the Major above: the report claimed the width of 1 keeps every round full, where at 16
in flight the last of thirty-eight rounds holds 8; and `serial.rs`'s helper claimed four in
flight puts several builders in one round at every width the file uses, where at width 600 the
analysed stretch is one region and every count gives the same one-builder round. Both are
corrected in place. `render` — which is what "the same answer" means for every comparison in the
module — read `RegionOutcome`'s fields rather than destructuring it, so a field added to that
type would have dropped silently out of all of them; it now destructures, as the driver beside
it does.

## 6. Validation

In the container (`./scripts/dev.sh`):

- `cargo fmt --check` — clean.
- `cargo clippy --lib --all-features -- -D warnings` — clean.
- `cargo test --lib ng::run::cohort_merge` — `226 passed; 0 failed`.
- `cargo test --lib` — `3849 passed; 0 failed; 11 ignored; 0 measured; 0 filtered out;
  finished in 862.10s` (3,847 before this step).

## 7. Follow-ups

- **The round's tail is unmeasured** — spec §6.2 says the `RwLock` and owned-window
  alternatives should not be reached for before it is.
- **That the builders occupy several threads is untested**; the type-level guard on
  `in_region_order` stands in for it.
- **Nothing in CI guards the merge's cost.** The repo runs ten criterion benches and neither the
  merge walk nor the parallel driver is among them. A bench is a project decision, raised rather
  than taken.
