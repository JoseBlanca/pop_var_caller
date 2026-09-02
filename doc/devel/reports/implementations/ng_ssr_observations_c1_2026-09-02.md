# ng STR observations — C1: the tract generator says why its reader dropped reads

*2026-09-02. Step C1 of
[`run_ssr_observations.md`](../../ng/impl_plan/run_ssr_observations.md), realizing
[spec §3.2](../../ng/spec/run_ssr_observations.md)'s first accounting debt. Branch
`ng-ssr-observations`.*

## Plan

`LocusGenerator::read_filter_counts` defaults to an empty list, and `SsrGenerator` did not
override it. The 2026-09-01 run-report work recorded that as the thing that *"will bite whoever
fills a tract slot"*: with the slot filled and no override, every run's per-read-group drop
rates would count only the SNP/indel generator's reader and report them as the sample's.

The generator already keeps its own cursor, one per chromosome, and `SampleCursor` already
tallies per read group. What was missing is the same two pieces the SNP/indel generator has: a
map of the retired chromosomes' tallies, and an accessor that adds the live cursor's to it.

## Assumptions

None the spec left open. The shape is `PileupGenerator`'s, deliberately — same field name, same
`BTreeMap<Option<ReadGroupId>, ReadFilterCounts>`, same "retired are folded in as they go, the
live one is asked directly" so the answer is current at any moment.

## Changes made

`src/ng/locus_generation/ssr.rs`: `SsrGenerator::retired_read_group_counts`, folded in where the
cursor is retired at a chromosome boundary — beside the aggregate `retired_cursor_counts`, which
is taken at the same line and for the same reason; `SsrGenerator::read_filter_counts`; and the
`LocusGenerator` override that makes it reachable once the generator is boxed.

## Tests added

- `the_generator_reports_its_readers_drops_per_read_group_across_chromosomes` — two kept reads
  and two duplicates, one pair on each chromosome, read after the first chromosome and again
  after the second: (1 kept, 1 dropped) then (2, 2). A generator that lost the retiring cursor
  answers (1, 1) the second time. It then asserts the same answer comes back **through the
  trait**, which is the only way a run can ask.
- `the_tract_generators_per_group_tallies_and_its_aggregate_count_the_same_reads` — the two
  harvests of one cursor are read by different callers and nothing compared them.

## Validation

In the dev container. `cargo fmt --check` and `cargo clippy --all-targets --all-features -D
warnings` clean; `ng::locus_generation` 378 passed, 0 failed.

**Two mutations, both killed**, and both are the silent failure this step exists for:

| mutation | what a run would have reported | tests failed |
|---|---|---|
| the trait override answers the empty default | no tract drops at all, every other number right | 1 |
| the retirement fold is dropped | the last chromosome's drops as the whole walk's | 2 |

Restored from a backup, and the mutations' absence checked by grep before committing.

## Tradeoffs and follow-ups

**Nothing reads this yet.** The slot is still unfilled — that is C2 — so the override is
correct and inert until then. The run report's own summing is already written for it: it walks
whatever generators the set holds.
