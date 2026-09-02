# G1 — record leasing, built and measured: it does not pay

**Date:** 2026-09-02. **Plan:** [`run_driver_direct_mode.md`](../../ng/impl_plan/run_driver_direct_mode.md),
Milestone G step G1 — deferred on 2026-09-01, re-opened on the owner's ruling of 2026-09-02
("I meant to defer it more than to drop it"). **Branch:** `ng-record-lease`, where the
experiment is commits `732f1068`, `f7c561b4`, `af92c896`, `d57c234d`, reverted in `849d283e`.

---

## The answer first

**G1 works and buys nothing.** The walk fills 48% of its records into ones the merge handed
back, the run gets no faster at any thread count, and it holds about 13% more memory. The
code is reverted; this report is what the branch is for.

**The reason is not a defect in the implementation, and it is worth more than the
milestone was.** G was written on a measurement of a cross-thread `free` taking a *locked*
path. That lock was the macOS system allocator's, and the binary was using the system
allocator only because `src/main_exp.rs` never declared a `#[global_allocator]` — a defect
found and fixed earlier the same day (`1f40c833`). With mimalloc a cross-thread free goes
to a per-thread heap without a lock, and a small allocation there costs about what the
refill's copying and bookkeeping cost. **The allocator fix already took the win leasing was
aimed at.**

---

## What was built

Exactly what the plan specifies, in three steps.

1. **`732f1068` — the record's byte buffers become growable.**
   `SampleLocusObservations::reference_bases` and `SequenceObservation::bases` were
   `Box<[u8]>`, which has no spare capacity to refill into. Both became `Vec<u8>`. The psp
   file's bytes are unchanged — the writer length-prefixes the same slice either way — and
   the whole cost is eight bytes per observation and per record while they are alive. 34
   files, all mechanical.

2. **`f7c561b4` — the ordinary-column lane fills a pooled record.** A `RecordPool` on
   `WalkerState`, beside the walk's other scratch, bounded at two records. The lane takes a
   record, overwrites every field, and reuses the buffers.

   Two shapes were chosen so a later field cannot arrive stale. The fill rebuilds each
   observation as an **exhaustive struct literal** whose two vectors are the old buffers
   taken out and cleared, so a field added to `SequenceObservation` fails to compile at the
   fill site rather than silently keeping the previous locus's value. And the empty record a
   cold pool hands out is an `empty_shell()` constructor rather than a `Default` impl, so an
   observation of nothing is not a value anything can build by accident.

   **The bound is asserted, which is the test the milestone owed.** B1's suite could pin
   that a returned record does not come back out as an observation but not that it is
   *released*, because a walker that drops it has nothing to count — a walker that stashed
   every offered record passed all fourteen of those tests, the one survivor of a
   twenty-one-mutation pass. `RecordPool::put` refuses past its bound and a test offers it a
   thousand records and asserts the count.

3. **`af92c896` — the seam.** `ObservationSource::next_observation` had taken the spare back
   since Milestone B and the walker dropped it on its first line, so eviction never avoided
   a free: it deferred one into the draw, and the draw allocated a fresh record four layers
   below. The walker now offers it to its generators; the trait method is defaulted to
   dropping so the six generators that mint no records need no change; the generic slot
   routes it to its chromosome walk's pool.

4. **`d57c234d` — surplus observation slots go back to the pool, not to the allocator.** A
   first version shortened a record's observation list by dropping the surplus, which frees
   two buffers and allocates them again at the next locus that needs the slot — churn on
   exactly the fluctuation that is most common at three reads a position. Fixed; it changed
   neither the time nor the memory.

---

## What it does, counted

The deterministic gate, because the machine could not give a trustworthy clock. A temporary
counter in the pool, tomato benchmark, 63 accessions, the first four BED regions of
`benchmarks/tomato1/regions.bed`:

| | count |
|---|---:|
| records the lane asked the pool for | 38,384,881 |
| — filled into a record handed back | **18,424,572 (48%)** |
| — allocated fresh | 19,960,309 |
| records the merge offered back | 21,621,112 |
| — accepted | 18,424,572 (85% of offers) |
| — refused at the pool's bound | 3,196,540 |

The draw count is exactly the lane's own emitted-column census (`FAST_COLUMNS`,
38,384,881), which is the check that the two agree.

**The ceiling is the cache, not the pool.** The merge offers a spare only when its own
cache has one — 21.6 M offers against 38.4 M draws, 56%. Raising the pool's bound would buy
at most the gap between 48% and 56%; the rest is a question about
`ObservationCache::evict_before`, not about leasing.

---

## What it costs

Arms alternated in one loop against `main` at `e270ce14`, same ground, same cohort, native
host release build. **Two virtual machines were taking about five of this host's eighteen
cores throughout**, so minima are quoted beside medians and no single run should be read on
its own.

**One thread**, five rounds:

| | wall, five rounds | min | median | peak resident |
|---|---|---:|---:|---:|
| `main` | 22.15, 21.90, 21.76, 25.90, 22.24 | 21.76 s | 22.24 s | 966–989 MB |
| with G1 | 22.53, 26.35, 23.37, 23.10, 25.79 | 22.53 s | 23.37 s | 1,118–1,124 MB |

**Eighteen threads**, five rounds:

| | wall, five rounds | min | median | peak resident |
|---|---|---:|---:|---:|
| `main` | 9.18, 9.37, 11.01, 10.10, 9.98 | 9.18 s | 9.98 s | 1,192–1,224 MB |
| with G1 | 9.23, 9.76, 11.41, 8.95, 12.15 | 8.95 s | 9.76 s | 1,323–1,397 MB |

**The time is a wash and the memory is not.** At eighteen threads the two overlap
completely; at one thread G1 is slower on both statistics. Peak resident is up about 13% in
every round of every arm, which is the one figure the background load does not move.

An earlier five-round pair at one thread had `main` at 22.4–23.9 s and G1 clustered at
33.0–33.9 s — a 40% gap that did not reproduce in the two later experiments. It is recorded
because it is the reason this report quotes five rounds and two thread counts rather than
one number: **at this host's load, a single alternated pair cannot tell 4% from 40%.**

### Where the memory goes, and it is inherent rather than a defect

Keeping capacity is what a pool is for. Here the capacity kept is the wrong size: a refilled
record's `Vec<SequenceObservation>` holds the capacity it **doubled** to — four — where
`.collect()` previously sized it exactly to the one or two observations a locus at three
reads a position carries.

At 63 accessions the round width the run chooses is 7,936 bases, so about **499,968 records
are live in the merge's window at once**. Three surplus slots at 112 bytes over that many
records is **168 MB predicted against about 150 MB measured**. Shrinking the vector on each
refill would give the memory back and undo the saving in the same line.

---

## What this changes about the plan

**Milestone G is closed by measurement rather than by assumption**, which is what its own G2
step asked for and could not get: G2's report says in terms *"None of this is wall time …
Getting the time would take either building G1 and timing it, or a sampling profiler — and
this machine cannot run one."* Both have now happened.

Two of the numbers G was written around should not be quoted again without the allocator
attached:

- *"one atomic instruction inside `free` is 10.6% of the merge's whole CPU"*
  (`cohort_merge_parallel_cost_2026-08-28.md` §2.2) — measured under the system allocator.
- *"removes 92% of the merge's frees"* — true, and the frees turned out not to be worth 92%
  of anything once they stopped taking a lock.

**What is still open, and it is not this.** The profile of the shipped binary puts 65% of
CPU in ng's own code with the ordinary-column lane's own arithmetic the largest single item
at 13.4%; the allocator, after mimalloc, is 7% of CPU in total. Leasing was aimed at a
fraction of that 7%. Anything worth another step change is in the lane's arithmetic or in
the decode/calling overlap, not in allocation counts.

---

## Provenance

Every number above was produced by the scripts in
`tmp/perf_review_2026-09-02_ng-calling/` — `run_ng.sh` for the command, `ab.sh` for the
alternated arms — against binaries built from the commits named at the top. The VCF is
byte-identical to `main`'s on the same cohort and ground in every arm, and `cargo test
--lib` was 5,914 passed, 0 failed with the feature in (three added by it) and is 5,911 with
it reverted.
