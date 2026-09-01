# G2 — what leasing the merge's records would be worth, measured before it is built

**Date:** 2026-09-01. **Plan:** [`run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md),
Milestone G. **Ruling this follows:** the owner's, 2026-09-01 — *do G2 before G1*, the same
order that made E1 come out right. **Design:** [`cohort_merge.md`](../../ng/spec/cohort_merge.md)
§5, §6.4; [`run_streaming.md`](../../ng/spec/run_streaming.md) §3.4, §5.5.
**Source under measurement:** [`cohort_merge_parallel_cost_2026-08-28.md`](../../ng/research/cohort_merge_parallel_cost_2026-08-28.md)
§2.2, §3, §5.5.

---

## The answer

**On a calling run, the per-sample observation records are 24% of what the calling phase
allocates — not the 92% the milestone was written around.** The share is flat across the whole
cohort range the tomato benchmark offers.

| samples | blocks allocated, calling phase | records drawn | observations in them | leasing could remove |
|---|---|---|---|---|
| 3 | 10,848,402 | 549,244 | 571,859 | 2,242,206 — **20.7%** |
| 12 | 38,817,290 | 2,190,217 | 2,269,125 | 8,918,684 — **23.0%** |
| 24 | 74,846,956 | 4,376,026 | 4,517,272 | 17,786,596 — **23.8%** |
| 63 | 194,871,747 | 11,405,016 | 11,844,580 | 46,499,192 — **23.9%** |

Tomato accessions from `benchmarks/tomato1/crams/` over the **first two intervals** of
`benchmarks/tomato1/regions.bed` — 200 kb of SL4.0 at about three reads a position — in the
development container, release build, `NG_COVER=serial`:

```
./scripts/dev.sh cargo run --release --example ng_call_cohort_end_to_end \
    --no-default-features --features dhat-heap,merge-timing -- \
    $HOME/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
    benchmarks/tomato1/crams/S_lycopersicum_chromosomes.4.00.repeats.parquet \
    benchmarks/tomato1/regions.bed benchmarks/tomato1/crams
```

**The allocator is dhat's, not the one a run ships with.** The project forbids `unsafe`, so it
cannot install a counting allocator of its own; `dhat::Alloc` is the route every other
`examples/dhat_*.rs` takes. Mimalloc is what ships. **So the block counts transfer and the times
do not** — which is why this report quotes no time.

## Why it is a count and not an attribution

The share is arithmetic on the record's own shape rather than a profiler's guess at whose
allocation is whose. An emitted `SampleLocusObservations` owns exactly:

- one allocation for `reference_bases`,
- one for the `observations` vector,
- and, per observation, one for its `bases` and one for its `chain_ids`

— so `2 × records + 2 × observations` is its whole heap footprint, at both mint sites
([`fast_column.rs:373`](../../../src/ng/locus_generation/pileup/fast_column.rs),
[`open_record.rs:1115`](../../../src/ng/locus_generation/pileup/open_record.rs)). On the fast
lane the chain-id vector is `mem::take`n out of the column scratch, which leaves that slot empty
to regrow, so it is an allocation either way.

**It is an upper bound on what leasing removes**, twice over: reuse only avoids an allocation
where the returned buffer is already big enough, and two of the four classes are `Box<[u8]>`,
which has no spare capacity to grow into.

**The counts are taken where every record enters the merge** —
`ObservationCache::draw_next` — rather than at the mint, because the generator also produces
records the region clamp discards, and the merge never owns those, so it never frees them.

## Why this differs from the 92% the plan quotes

Both numbers are right about different denominators. The research note counted **21.4 of 23.1
million frees, 92%** ([§5.5](../../ng/research/cohort_merge_parallel_cost_2026-08-28.md)) on a
probe that handed the merge records **made before its clock started**. Its denominator therefore
held the merge's own allocations and nothing else.

A calling run decodes reads. Every read becomes a `MappedRead` owning its sequence, its
qualities and its CIGAR, and that is where three quarters of the allocation goes. The 92% was
never wrong; it answers *what fraction of the merge's own traffic is records*, and the question
Milestone G exists to serve is *what fraction of a run's*.

## The other two ceilings, restated because the plan's copy of them is stale

- **⚑ The plan's G2 asks for the merge to be reported against "2.6–18% of walk-plus-merge".
  Those two numbers were retracted by their own source.** The research note's §3 says in terms:
  *"Those replace the 2.6% and 18% this section carried before the record-making was
  measured"* — the honest pair is **1.4%** with the walk serial and **10%** with the walk on
  eight threads. The plan's sentence is the owner's and is not edited here.
- **Where the frees actually sit in a calling run.** E1 measured `call_cohort` at 63 accessions
  as 88.1% drawing the readers, 0.8% evicting, 5.5% assembling, 5.3% genotyping. The record
  frees are inside the 88%, not the 0.8%: eviction moves a record to the cache's spare list, and
  it is freed later when the walker declines the offer
  ([`walker.rs`](../../../src/ng/run/walker.rs), `drop(spare)`), which happens inside the draw.
  So no existing timing column isolates them, and none of the three is the prize on its own.

## The ruling

**Milestone G is dropped — owner, 2026-09-01**, on the measurement above: *"drop G. No problem.
We'll work on the performance in future sessions. The critical objective right now is to have a
first working variant caller to be able to improve upon it in future sessions."*

The merge goes on freeing the records it was handed. **What replaces the checkpoint is the
number**, so the milestone can be re-opened on evidence rather than on the retracted 92%.

## What this leaves Milestone G

**G1 is a large change in the hottest code in the caller.** The walker cannot refill anything:
a record is minted four layers below it — walker → `SampleLocusObservationsIterator` →
`GeneratorSet` → `PileupGenerator` → the chromosome walk → one of the two mint sites — so
leasing needs a new method on the `LocusGenerator` trait (7 implementors), a spare threaded down
through `PileupWalker::fill_pending`, and changes at both mint sites. It also has to overturn a
decision already written down: `finalise_recycling`'s own comment says the emitted bytes are
"the one part of a record that genuinely has to be new each time".

**Against that, every ceiling is now measured rather than assumed:** at most 24% of the calling
phase's allocations, sitting inside a decode that is 88% of `call_cohort`, in a merge that is
1.4–10% of walk-plus-merge.

**What is still not measured, and what it would cost to measure.** None of this is wall time. A
share of allocator *calls* bounds a share of allocator *work*, and neither is a share of the
run. Getting the time would take either building G1 and timing it, or a sampling profiler —
and this machine cannot run one: `perf_event_open` is gated by the host sysctl on the Linux box
and on macOS the container targets the VM kernel rather than the host (`CLAUDE.md`).

## What landed, so the measurement can be repeated

Two counters and a heap report. Both are no-ops in a build without their feature, so nothing
here is carried by a shipped run.

- **`timing::RECORDS_DRAWN` and `timing::OBSERVATIONS_DRAWN`**
  ([`cohort_merge/timing.rs`](../../../src/ng/run/cohort_merge/timing.rs)), incremented in
  `ObservationCache::draw_next`, cleared by `timing::reset` with the rest.
- **`report_the_heap`** on `examples/ng_call_cohort_end_to_end.rs`, printing the calling phase's
  blocks and bytes and the record share, under `--features dhat-heap`. Without the feature it
  prints the two counts and says which build gives the share, so the share cannot be guessed at
  from the counts alone.

## What was measured

| check | result |
|---|---|
| `cargo test --lib` | **5,908 passed, 0 failed, 14 ignored** — unchanged by this step |
| `cargo test --lib ng::run` | **438 passed** — unchanged by this step |
| `cargo fmt --check` | exit 0 |
| `cargo clippy --all-targets --all-features -- -D warnings` | exit 0 |
| `cargo clippy --all-targets --no-default-features --features dhat-heap,merge-timing` | exit 0 |
| `cargo doc --no-deps` | 26 `error: unresolved link`, 23 `warning: redundant explicit link target` — the standing baseline |
