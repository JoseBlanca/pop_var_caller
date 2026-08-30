# What the cohort merge's parallel driver costs, and what recovers it

**Date:** 2026-08-28. **Status:** finding. **Plan:** [`cohort_merge_parallel_cost_plan.md`](cohort_merge_parallel_cost_plan.md).
**Branch:** `ng-merge-parallel-cost`.

The parallel merge gives back much less speed than its thread count suggests. This is what its
time is made of, which of the candidates were worth anything, and which were built and set aside.

**Where every number comes from, and it is one corner.** 63 tomato accessions
(`benchmarks/tomato1/crams`) over one 100 kb interval of SL4.0, walked through the generic locus
generator: 6.09 million per-sample records, about one per covered base per sample, 1.04
observations per record, 12,029 cohort loci built and 87.6% of closed loci discarded as too
quiet. The allocator is mimalloc, which is this project's default. The machine is an Apple M5
Pro (18 cores) running the dev container, which gets 8 CPUs and 16 GB. **A figure here is a fact
about 63 samples at that density on that machine, not a property of the merge** — §6 gives what
changes with cohort size.

---

## 1. The answer in one paragraph

**The round barrier is not what costs.** Waiting for each round's slowest builder is 3.9% of the
merge and launching the round's builders is another 0.5%. What caps the speed-up is that **eight
threads do 2.08 times the total work one thread does for the identical answer** — 4.90 seconds of
CPU against 10.18, summed over eight merges — and the extra is almost entirely the allocator, because a record is allocated
by whichever worker drew it and freed by whichever worker later evicts it, and mimalloc takes a
locked path when those differ. Eight threads divided by 2.08 times the work is 3.8× at best; the
measured speed-up is 3.1×, and the gap is the part that ran on one thread while seven were idle.

---

## 2. Where the merge's time goes

Measured by stopwatches inside the merge itself, behind the `merge-timing` cargo feature — every
counter is a zero-sized type when the feature is off, so the drivers' source reads the same in
both builds. Overhead when it is on: **2.9%** (238 ms plain against 245 ms timed).

At 1,000-base building regions, 8 threads, 8 merges summed:

| part | share | does it spread over threads? |
|---|---|---|
| drawing every sample's reader forward | 63.7% | yes, and nearly perfectly |
| evicting what the round has passed | 14.9% | **no — one thread** (until §5.3) |
| the builders' own work | 16.8% | yes |
| waiting for the round's slowest builder | 3.9% | — |
| launching and collecting the builders | 0.5% | — |
| releasing loci in region order | 0.1% | — |

**The drawing phase spreads almost perfectly and still does not scale.** Its wall time is its
summed work divided by the thread count to within 7% at every thread count measured — so the
spreading works. What grows is the work:

| threads | one merge | drawing: wall | drawing: work summed | builders: work summed |
|---|---|---|---|---|
| 1 | 650 ms | 3,465 ms | 3,463 ms | 1,434 ms |
| 2 | 367 ms | 1,798 ms | 3,478 ms | 1,661 ms |
| 4 | 250 ms | 1,249 ms | 4,790 ms | 1,845 ms |
| 8 | 210 ms | 1,062 ms | 7,930 ms | 2,246 ms |

Taken in one process, arms in that order with the 1-thread arm repeated last, because the first
arm of a run is up to 1.6× slower cold and **this machine has drifted by a factor of two between
two runs of one unchanged binary**.

### 2.1 Where the extra work goes

From a sampled profile of the merge alone (`samply`, inside the container with
`perf_event_paranoid` lowered on the VM's own kernel):

- **43.5%** of the merge's CPU is inside the allocator.
- **10.6%** of *all* of it is one atomic compare-and-swap instruction. 86% of those samples are
  called from the allocator's `free`, and 94% of `free`'s own time is dropping per-sample
  observation records.

### 2.2 Which records — and whose cost they are

| | share of the merge's CPU | whose work is it? |
|---|---|---|
| making a per-sample record | 25.5% | **not the merge's** — the generator or psp reader makes it once, upstream |
| destroying a per-sample record | 25.9% | **the merge's** — it is the last owner; it evicts and drops |
| assembling the cohort loci | 21.0% | the merge's |
| the cache's bookkeeping, eviction's walk, the organiser | the rest | the merge's |

**The cohort loci are not where the allocator time is**: dropping the built loci is 1.6% of
`free`'s time against 93.7% for the per-sample records, and 96.8% of the allocator's fast path is
entered from making a per-sample record. The reason is arithmetic on the counts — 6.09 million
records in for 12,029 loci out, 506 to 1, because the keep rule discards 87.6% of the positions
it closes.

**A caveat that changes what may be concluded, and it has now been measured.** The probe *creates*
each record inside the merge's clock, by cloning a template, standing in for a generator that mints
one from reads. So the 25.5% row above is charged to the merge in every other wall-clock figure in
this document and should not be. Giving the probe a source that hands over records made **before**
the clock starts prices it — 500-base regions, 8 threads, 10 merges an arm, alternated in one
process:

| records | pass 1 | pass 2 |
|---|---|---|
| made inside the merge's clock | 203.5 ms | 202.9 ms |
| made before it, handed over | 114.0 ms | 114.2 ms |

**So 44% of what this document calls "the merge" is making the records, and the merge's own work at
63 accessions is 114 ms rather than 203.** Every other merge figure here is therefore about 1.8
times what the merge itself costs; the *shares* in the tables are unaffected, because they divide
one measured merge by another part of the same one.

The 25.9% row — the freeing — is genuinely the merge's, because nothing owns those records
afterwards. **It could not be isolated in wall clock, and §5.5 says why.**

---

## 3. Is the merge worth optimising at all?

Nobody had ever timed it against another stage. Producing these observations took **7.81 s** for
the 63 samples one after another on one thread; the merge's own work, with the records made before
its clock starts, is **0.114 s** (§2.2).

- With the walk as it stands, serial: the merge is **1.4%** of walk-plus-merge.
- With the walk given all eight threads: **10%**.

**Those replace the 2.6% and 18% this section carried before the record-making was measured**, and
they are the honest pair: a run's generator makes each record once, and the merge should not be
charged for it twice.

That bracket is the ceiling on everything in this document, and it does not cover the rest of a
run — the parameters fit, the calling loop and the VCF writing have never been timed either.

---

## 4. What each candidate was worth

| candidate | verdict |
|---|---|
| widen the building region, 200 → 500 bases | **adopt** — faster in every sitting, +3% peak resident |
| eviction on the pool instead of one thread | **adopt, small** — 4–5%, and less at larger cohorts |
| give evicted records back to the producer | **the largest lever, not built** — 25.9% of the merge's CPU; needs a producer that leases, so it is milestone G of the run driver's plan (§5.5) |
| overlap the reader advance with the building | **refuted** — 2–4% slower, +52% peak resident |
| fold the held window in by bisection | **refuted** — no change; the window is 0.9 records |
| drop the rounds for a sliding window | **not built** — the owner dropped it once the barrier priced at 3.9% |

---

## 5. Each one, with its numbers

### 5.1 The building region's width

At 8 threads, one sitting, medians of 5:

| width | 20 | 100 | 200 | 500 | 1,000 |
|---|---|---|---|---|---|
| merge | 579 ms | 406 ms | 260 ms | 220 ms | 219 ms |

**The recommendation rests on the ordering, not on those milliseconds.** Five sweeps were taken
under different host loads and their absolute numbers span a factor of two; every one of them puts
200 slower than 500, and 500 level with or slightly slower than 1,000. Two sweeps taken while the
host was indexing a fresh build tree are excluded outright (§8).

Peak resident rises monotonically and cheaply: 3.148 GB at 200, 3.244 at 500 (+3%), 3.343 at
1,000 (+6%), 3.590 at 2,000 (+14%). **500 bases buys the speed for 3% more memory**; 1,000 is not
reliably better than 500 and 2,000 is worse.

The width is also the plan's diagnostic, and it points the way the profile does: at 20 bases,
launching the round's builders costs 12.0% of the merge and waiting for the slowest costs 5.0%; at
1,000 bases those are 0.5% and 3.9%. So the loss at narrow regions is per-round overhead, not
stragglers.

**This overturns the current default's own documentation**, which says the eight-thread optimum is
100–200 bases (`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`). That note's figures were taken at 16
samples and across separate runs.

### 5.2 What the fabricated fixture got wrong about real ground

- **No building region was empty.** 0 of 5,000 at every width swept. On the fabricated fixture the
  skip for an empty region was worth about a third of the merge, because there every sample carries
  a record at the *same* positions.
- **The fixed cost per building region is small**: assembling each sample's window plus opening the
  walk is 1.0% of the merge at 1,000 bases and 3.4% at 20.

### 5.3 Eviction on the pool

Eviction walks every sample, shares nothing between them, and ran on the organiser's thread
between rounds while the other seven sat idle. Putting it on the pool, measured before and after,
alternated twice on a quiet machine:

| | eviction | whole merge (8 merges) |
|---|---|---|
| one thread | 212.1 / 237.9 ms | 1539.9 / 1759.7 ms |
| on the pool | 102.8 / 107.6 ms | 1525.4 / 1628.6 ms |

Eviction halves — 225 ms mean to 105 — but about a third of that comes back as a slower cover
(1,076 → 1,113 ms mean), because records are now freed by whichever worker drained that sample and
the recycled ones cross threads more often on the next draw. **Net 4.4%.**

### 5.4 Two changes built, measured and set aside

**Overlapping the reader advance with the building** (the plan's §5a). A round's builders and the
next round's cover run in one `rayon::join`; the builders read what the last cover released while
the cover appends to buffers no builder can see. It does what it was meant to — drawing falls from
65% of the merge to 4% — and it is **2–4% slower**:

| driver | pass 1 | pass 2 |
|---|---|---|
| round at a time | 197.8 ms | 206.3 ms |
| overlapping | 205.9 ms | 210.2 ms |

**The phase it hides was never idle time.** Inside the overlapped round the builders' summed work
is 2,076 ms and the cover's is 7,179 ms; over 8 threads that is 1,157 ms perfectly spread, and the
round's wall is 1,181 ms — 98% of perfect. There is no idle capacity to hide anything in. What
remains is new cost: releasing what the cover drew, 6.6% of the merge, and eviction up 31% because
the window it drains holds two rounds. Peak resident goes 3.34 GB → 5.07 GB.

**Folding the already-held window in by bisection instead of walking it.** A sample's records are
disjoint and ascending, so the reach is monotone and the records a given reach admits are a prefix —
findable by search rather than by walking. It changes nothing, because **the window a cover starts
from is 0.9 records**: eviction runs immediately before each cover, so by then there is nothing to
walk. It is load-bearing only in the overlapping driver, which cannot evict while builders read and
whose covers start from **7,145** records — which is why that driver came in 2–4% slow rather than
far worse, and so is why its refutation is a fair one.

### 5.5 The lever that was mis-measured, and is the largest

**Give the merge's evicted records back to the producer instead of freeing them.** The interface
hook exists: `ObservationSource::next_observation` takes a spare record back and the cache already
offers one. Nothing upstream fills it.

Counted with dhat over the same ground, leasing removes **21.4 of the 25.8 million heap blocks a
merge allocates and 21.4 of the 23.1 million it frees** — 83% and 92%.

**The first verdict on this was wrong and is retracted.** Timed with the probe, leasing bought
nothing on eight threads (188 / 198 ms against 190 / 194) and cost 40% on one (846 / 856 against
621 / 587), and it was reported as refuted. **The probe charges the merge for filling the record in
both arms** — minting clones a template, leasing refills the returned buffer, and both happen inside
the merge's clock — so what was compared was "clone then free" against "refill", which cost about
the same. A real run pays neither: the generator fills the record either way, and leasing removes
the merge's *free*, which §2.2 measures at 25.9% of the merge's CPU.

**Three attempts to settle it, and none did.** Re-run on the current driver with the arms
alternated four times, leasing and minting are indistinguishable — 181 / 232 ms against 191 / 234,
a 2.6% mean difference inside a run whose own second pass was 20% slower than its first. And the
profile cannot separate the fill from the rest, because the probe's `refill` is inlined and has no
symbol left to attribute.

And the third: a source that hands over records made before the clock and **keeps** the ones the
merge gives back, against one that lets the merge drop them — the two differing in the freeing and
nothing else. That fails for a reason worth recording: **the device that stops the merge freeing must
hoard**, which grows the process from 5.2 GB to 6.4 GB, and the memory pressure costs more than the
free saves. Hoarding measured 109.2 and 134.4 ms against handing's 114.8 and 112.3 — one arm either
side, on a spread of 94 to 165.

**So this is what stands.** The merge's *freeing* of per-sample records is 25.9% of its CPU, by
attribution from the minting profile (§2.2), and leasing removes 92% of those frees by count. What
is **not** established is the wall-clock saving: it cannot be had by removing the free from this
probe, because nothing here can accept a returned record without either refilling it — which is the
generator's work, charged to the merge — or holding it, which costs more than it saves.

**Settling it needs a real leasing producer**, one that fills a returned record instead of
allocating a new one — which is a change to how the generator or the psp reader gives the merge its
records, not to the merge. That is where the measurement belongs, and it is now written down as the
last milestone of the run driver's plan
([`../impl_plan/run_driver_direct_mode.md`](../impl_plan/run_driver_direct_mode.md), milestone G).

---

## 6. What changes with cohort size

The tomato benchmark holds 63 accessions, so the range above 63 is bracketed on fabricated ground —
a record every four bases, every sample at the same positions, which is far sparser per region than
real data. **What transfers is the direction, not the percentages.**

| cohort | drawing readers forward | eviction | building the loci |
|---|---|---|---|
| 63 | 56.6% | 12.2% | 30.7% |
| 1,000 | 48.8% | 4.2% | 46.9% |
| 3,000 | 37.2% | 3.0% | 59.7% |

Two consequences, and the first corrects an argument made while this work was running:

1. **Eviction does not become more important at large cohorts** — its share falls, because
   assembling loci grows faster than throwing records away. The case for putting it on the pool is
   the 4.4% at 63 samples and nothing more.
2. **Parallelism pays better as the cohort grows**, because the phase that grows fastest is the one
   that spreads cleanly. Against the honest one-thread baseline that uses no pool at all: 1.6× at 63
   samples, 2.1× at 1,000, 2.9× at 3,000.

---

## 7. How the answer was kept identical

**The one hard constraint is that the output cannot change**, and the fixture suite is not enough
for a schedule change. `ng_cohort_merge_real_cost` now compares each driver against
`merge_cohort_serially` — which holds every sample's observations at once and divides the ground not
at all — **on the real tomato observations**, at every width swept and at 1, 4 and 16 regions in
flight. Both drivers give the oracle's answer over all 12,033 lines of output. The 258 cohort-merge
tests stay green.

---

## 8. Measurement hygiene, because most of the effects here are smaller than the noise

- This machine drifted by a factor of two between two runs of one unchanged binary, and by 1.6×
  between the first and last arm of one run. **Every comparison that decides something in this
  document was taken inside one process**, which is why the probe's thread count, width, record
  supply and driver are all lists.
- Two width sweeps taken while Spotlight was indexing a fresh build tree — host load 24 on 18 cores —
  disagreed with each other and are excluded.
- Sweeping a knob in ascending order confounds it with warm-up; the arms are run as a palindrome.

---

## 9. The answer owed to `run_streaming.md` §11, question 7

That question asks **which stage of the psp-to-VCF path is worth a pool at all**, and names three
ideas for the merge half of it. This is the merge half's answer; the decode and calling halves still
need the run-level profile the question asks for, which nothing in this work touches.

**The merge is already worth a pool, and it is not where a run's time is.** Eight threads give
**3.1×** on 63 tomato accessions at 1,000-base regions — not the 1.4× that question currently
quotes, which came from 200-base regions at a cohort size where the per-round overhead is a fifth of
the merge. The merge is **2.6% to 18%** of walking-plus-merging depending on whether the walk is
given the machine, so tripling it moves a few per cent of those two stages.

The three ideas the question names, settled:

1. **Sweep the building-region width** — done, and it is the one to act on. 500 bases against the
   shipped 200, for 3% more peak resident (§5.1). Setting the shipped default remains the owner's.
2. **Overlap the reader advance with the building** — built, measured, **refuted**: 2–4% slower and
   52% more peak resident, because the phase it hides was already saturating the pool (§5.4).
3. **Drop the rounds for a sliding window** — **not built, and dropped by the owner on the strength
   of the profile**: the barrier it removes is 3.7% of the merge and the launching another 0.6%,
   while the machinery it needs cost the overlapping experiment 6.6% before removing any barrier.

**What the question did not name, and what the profile says to do instead:** stop the merge freeing
the per-sample records. That is 25.9% of its CPU, the interface hook already exists, and it needs
one measurement this work did not manage to take honestly (§5.5).

**The concurrency defaults this document can speak to.** `cohort_locus_builder_regions_len`: 500.
`cohort_locus_builder_regions_in_flight`: one per worker thread is what a run already takes with no
value given, and nothing measured here argues against it — what it sets is memory, and the width is
the knob that moves time.
