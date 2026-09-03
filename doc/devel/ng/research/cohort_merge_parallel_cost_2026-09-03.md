# Finding — the whole calling run at 8 threads, and an instrument regression that made the merge look barrier-bound again

**Date:** 2026-09-03. **Status:** measurement finished; one code change adopted (the timing
feature's two per-record counters sharded). **Builds on the 2026-08-28 finding**
([`cohort_merge_parallel_cost_2026-08-28.md`](cohort_merge_parallel_cost_2026-08-28.md)), which
settled the merge-only questions; this session widens the measurement to the **whole calling
run** (the owner's ask) and repairs the instrument both findings depend on. Together they answer
the direct-mode half of [`run_streaming.md`](../spec/run_streaming.md) §11 question 7.

**The one-paragraph answer.** At 8 threads a 63-sample calling run over 200 kb of tomato ground
spends **84.5% of its wall decoding reads, 8.0% genotyping, 6.2% assembling loci** — the pool
belongs to the decode, and per-locus calling parallelism stays refuted at this corner. The
decode is already spread across samples, and its missing speed-up now has sizes: **the decode
work itself swells ~2.6-fold in total CPU when eight copies run at once** (the same phenomenon
the 2026-08-28 profile attributed largely to the allocator on the merge's ground), and the
per-sweep **wait for the slowest of the 63 samples costs a further ×1.48**. Two candidate
explanations died today: swapping mimalloc for glibc changes nothing, and the cover's fixpoint
re-sweeps are 14 extra sweeps in 414 — both refuted as costs. Separately, the **per-record
counters added to the timing feature on 2026-09-01 (G2) had been throttling every instrumented
parallel merge since**: at 8 threads they turned a 2.25× merge into a measured 1.40×. They are
sharded now, and the instrumented build again matches a build with no instrumentation (545 ms
against 552 at 8 threads).

---

## 1. Ground, machine, and trust

- **Data:** the tomato benchmark's 63 accessions (`benchmarks/tomato1/crams/`, sliced CRAMs,
  about three reads a position) over the first two intervals of `regions.bed` — **200 kb of
  SL4.0**; 16 accessions for one cohort-size contrast. A sample yields ~194,000 observation
  records over the 200,000 bases — about one per base — and a 63-sample run draws 11.5 M
  records carrying 12.1 M sequence observations. Facts about this corner, not the caller.
- **Machine:** the development Mac (Apple M5 Pro, 18 cores, 6 of them performance) running
  the Apple-container VM `dev.sh` gives 8 vCPUs and 16 GB; mimalloc unless said otherwise;
  release builds; medians. The host was shared with other work in some sittings — the serial
  oracle swung 17.6 → 26.2 s across sittings of one unchanged binary — so **only
  within-sitting comparisons decide anything below**, and cross-sitting reads are labelled.
- **Instruments:** `examples/ng_call_cohort_end_to_end.rs` (whole run; serial oracle against
  the run's record path; this session added lines printing the drawing's summed busy time
  beside its wall) and `examples/ng_cohort_merge_real_cost.rs` (the merge alone on real
  observations; **its runs here use 200-base building regions, the probe's default — the
  shipped default has been 500 since 2026-08-28**, so its scaling figures are the
  conservative width). Raw outputs under the worktree's `tmp/parallel_cost/`, not committed.

---

## 2. The whole run: where 8 threads go

`NG_COVER=parallel NG_SAMPLES=63 NG_REGIONS=2`, one process per pool size, one sitting,
timing feature on (§4 shows its run-level distortion is inside drift):

| rayon threads | calling wall | speed-up |
|---|---|---|
| 1 | 17.56 s | — |
| 2 | 11.39 s | 1.54× |
| 4 | 9.53 s | 1.84× |
| 8 | 8.77 s | 2.00× |

A bare build in a noisier sitting reads 17.73 / 11.28 / 9.60 / 9.60 s — the same shape, so
this curve is not the instrument's. Three alternated serial–parallel pairs, first sitting:
serial 17.62–18.16 s, parallel 8.62–8.74 s — consistent with Milestone E1's recorded 1.8×.

**Inside the 8.77 s at 8 threads** (the merge's own stopwatches):

| part | ms | share |
|---|---|---|
| drawing the readers forward (the decode, spread across samples) | 7,411 | 84.5% |
| genotyping the loci | 705 | 8.0% |
| assembling the loci | 548 | 6.2% |
| evicting what the run has passed | 83 | 1.0% |

Even with the decode parallel across samples it is five-sixths of the wall. **A genotyping
pool's ceiling here is 8%** — the arrangement Milestone E was once sketched around stays
unbuilt for the third measurement in a row. (Its share grows with cohort — 2.2% at 3 samples,
10.8% at 63, 2026-09-01 — and where that flattens is unknown; nothing turns on it until the
decode is fixed.)

## 3. The decode's missing 6×, split in two

New lines printing the drawing's two sides, 63 samples, one sitting
(`drawing_sides.out`):

| arm | drawing wall | the samples' own drawing, summed over threads |
|---|---|---|
| record path, 1 thread | 17.49 s | 17.49 s |
| record path, 8 threads | 8.25 / 8.29 s (twice) | **44.63 / 44.88 s** |

Two multiplicative losses against the one-thread work spread perfectly (2.19 s):

- **Work inflation, ×2.55–2.57**: the decode that costs 17.5 s of CPU alone costs 44.6–44.9 s
  of CPU when eight copies run together. Reproduced twice here and echoed in the merge-only
  probe (×3.36 there, §5) and in the 2026-08-28 finding (×2.08 at 8 threads on its ground,
  where a `samply` profile put **43.5% of the merge's CPU inside the allocator**, 86% of its
  hottest atomic instruction under `free`, freeing per-sample records).
- **Spread loss, ×1.48**: even the inflated work spread perfectly over 8 threads would be
  5.58 s of wall; the measured wall is 8.25–8.29 s. The cover advances all 63 samples to a
  ~500-base frontier and waits for the slowest, ~400 times a run — real decode is skewed
  (CRAM container boundaries land unevenly), where the 2026-08-28 minted-record ground spread
  to 98% of perfect.

**Refuted today, each by a measurement:**

- **The choice of allocator.** The identical run on glibc: 18.25 s at 1 thread, 9.41/9.82 s
  at 8 — indistinguishable from mimalloc's 17.73 and 9.60. Neither of the two allocators on
  offer is better; this does *not* clear the allocation **traffic** itself, which the
  2026-08-28 profile implicates and which only Milestone-G-style leasing (fewer
  allocations, not a different allocator) would test at run level.
- **The fixpoint's re-sweeps.** 414 cover sweeps for 400 working windows — the chain-closing
  iterations E1's report guessed at are 3% extra sweeps, nothing.
- **The frees, in isolation.** The merge probe's two hand-fed arms differ only in whether the
  merge's released records are dropped inside the clock or hoarded: 254 against 281 ms at 8
  threads — inside noise, agreeing with the 2026-08-28 §5.5 verdict that the free's
  wall-clock price cannot be shown this way.

**What remains, with the arithmetic a mechanism owes.** Scheduling 8 VM threads against 6
performance cores explains at most ×1.1–1.3 of summed busy; ×2.56 needs most of its size from
**shared cache/memory contention in the decode and its allocation traffic**. The inflation
curve (one noisier sitting): ×1.08 at 2 threads, ×2.34 at 4, ×2.65 at 6, ×2.85 at 8 —
near-nothing at 2 threads, steep from 4. **Splitting cache contention from allocation traffic
from VM scheduling needs `perf` on real cores — the Linux box — and that is the one question
this report hands on.**

## 4. The instrument regression: found, fixed, verified

G2 (2026-09-01, `f6d38ea9`) added two counters bumped **once per record drawn** —
`RECORDS_DRAWN`, `OBSERVATIONS_DRAWN`, two adjacent global atomics — to price Milestone G's
leasing. Under the parallel cover every worker bumps them at ~100 ns intervals and the shared
cache line ping-pongs; **every instrumented parallel-merge figure taken between 2026-09-01 and
today carries it.** The merge-only probe, 63 samples, 200-base regions, 10 rounds, medians:

| build | 1 thread | 8 threads | scaling |
|---|---|---|---|
| timing on, counters as G2 left them | 1,105 ms | 791 ms | **1.40×** |
| no timing feature at all | 1,244 ms | 552 ms | **2.25×** |
| timing on, the two counters sharded (adopted) | 1,181 ms | **545 ms** | 2.17× |

That 1.40× is numerically the figure §11 question 7 used to quote from a 16-sample sweep —
a coincidence that made the merge look barrier-bound again after the 2026-08-28 finding had
already shown 3.1× (at 1,000-base regions, clean: those counters did not exist yet, and that
session measured its instrument's whole overhead at 2.9%). The fix is in `timing.rs`: a
`ShardedCounter` of 16 cache-line-padded cells indexed by rayon worker, summed on read;
every other counter is touched at most once per sample per sweep and stays a plain atomic.
**The whole-run measurements never needed retracting** — there a record arrives every ~4 µs
per thread, the line does not ping-pong, and instrumented walls match bare walls within
drift, which is why E1's 1.8×/1.5× stands.

## 5. The merge today, at the conservative width

With the counters sharded, the trustworthy 8-thread split (63 samples, 200-base regions,
summed over 10 rounds): whole merge 5,438 ms = cover wall 3,903 (71.8%; its busy summed
27,996 — the ×3.36 inflation again, on minted records) + builders 1,206 (22.2%) + evict 318
(5.8%) + **slowest-builder wait 252 (4.6%) + launching 178 (3.3%) + ordered release 10
(0.2%)**. The scheduling terms the plan's §5a/§5b would remove total ~8% of a stage that is
~10% of walk-plus-merge here — the 2026-08-28 refutation of the overlap driver and the
owner's dropping of the sliding window both stand.

**The record-supply arms** (bare build, one process, 63 samples, `supply_arms.out`):

| records supplied by | 1 thread | 8 threads | scaling |
|---|---|---|---|
| minted — a fresh record per draw (what a run does) | 1,200 ms | 554 ms | 2.2× |
| leased — the returned record refilled | 1,110 ms | 490 ms | 2.3× |
| handed — made before the clock; frees inside | 954 ms | **254 ms** | **3.75×** |
| hoarded — made before; nothing freed inside | 904 ms | 281 ms | 3.2× |

The merge's machinery — cover bookkeeping, windows, builders, ordered release — parallelises
at 3.75×; **what stops the shipped arm at 2.2× is constructing records inside the clock**,
the same class of work as the run's decode. Leasing reads 12% at 8 threads here; the
2026-08-28 caveat applies unchanged — the probe's refill stands in for the generator's fill,
so the honest test of leasing remains a producer that leases, at run level.

**Small cohorts parallelise worse:** at 16 samples the bare merge reads 229 / 142 / 167 /
175 ms at 1 / 2 / 4 / 8 threads — best at **two** threads and degrading beyond, because
per-region fixed costs and scheduling swamp 16 samples' work. A thread default should not
scale past the cohort's ability to feed it.

## 6. What this leaves for the parallelisation plan

Ranked by measured size at this corner (16–63 samples, 200 kb, ~1 record/base/sample —
the 1-to-10,000-sample range this caller commits to extrapolates only where marked):

1. **The decode's ×2.6 self-slowdown is the whole game.** Until it is understood, threads buy
   2× and then stop. Next step: `perf` on the Linux box — real cores, no VM, and the
   allocator's share measurable again. If the inflation shrinks there, much of this was the
   Mac VM and the lab servers never had the problem; what persists is cache/allocation
   traffic, and the remedies are fatter per-sweep decode work and fewer allocations
   (leasing), not more threads.
2. **The per-sweep straggler wait (×1.48) is the second lever** and is a scheduling shape:
   samples advancing to a shared frontier in lock-step. Letting each sample draw ahead into
   its own window would decouple them; re-measure after (1), which also shrinks the skew.
3. **Buffer reuse / leasing** stays the adoptable small win (12% on the merge at 8 threads
   here; 83–92% of the merge's allocator traffic by the 2026-08-28 dhat count) and is the
   only run-level test of the allocation-traffic hypothesis on the table.
4. **Genotyping pools stay unbuilt** — 8.0% at 63 samples; re-measure the share past ~500
   samples before revisiting.
5. **The barrier/ordered-release rework stays refused** — ~8% of ~10%.

**Not measured, said out loud:** more than 8 threads (the VM's ceiling); cohorts past 63 on
real reads; the whole 8 Mb ground (host shared; E1's 1.5× stands as recorded); repeat-tract
ground (no catalog on this path — density, not record shape, drives these costs).
