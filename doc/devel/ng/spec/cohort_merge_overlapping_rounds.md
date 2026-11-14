# Overlapping the cohort merge's rounds

**Status: built for `call-from-psps`, 2026-11-13, on the owner's go-ahead** (`stream.rs`,
`observation_cache/round_view.rs`). §0 says what was built and what it gave; §1 to §5 are the
proposal as it was reviewed, kept because its measurements are what the decision rested on, and
§6 says where the build departed from it and why.

## 0. What was built, and what it gave

**A stream of regions.** Each round's records are sealed once read into a read-only view its
builders read (`ObservationCache::seal_round`); a psp source keeps its bodies in sealed stretches
of its arena that the views share (`ObservationSource::seal`); up to six rounds are held, each a
sixth of the run's ground; the next rounds are read on a thread and a pool of their own while
earlier rounds' regions are built; and regions go to the pool as their round is sealed.

Measured on the tomato25 cohort's stored files, 2,176 samples at about three reads a position,
48 threads, 24,000 bases of ground, refinement off, the VCF byte-identical to the round driver's:

| ground called | round driver | stream | peak resident, round driver → stream | CPU |
|---|---|---|---|---|
| 300 kb of chromosome 0 | 57–59 s | 43 s, 43 s | 21.2 → 21.3 GB | 41–42 → 42 min |
| 1 Mb of chromosome 0 | 3 min 17 s | 2 min 25 s | 25.8 → 33.3 GB | 1 h 57 → 2 h 00 |

So about a quarter off calling, for up to a third more memory. The 3× of §2 did not
materialise because two cheaper changes made before it (the fragment keeper taking records in
batches, and regions of 125 bases, four a thread) had already recovered most of the round's idle
tail: with them, the round driver kept 38–42 of the 48 threads running.

## 1. The problem, measured (as proposed)

The cohort merge calls the genome in **rounds**. A round:

1. drops the observations behind it (*evict*);
2. reads every sample's observations over its stretch of genome from the psp files (*cover*);
3. builds and calls the loci of its regions on the thread pool, one region per thread;
4. waits until every region has finished, and only then starts the next round.

Step 4 is the cost. A thread whose region is done has nothing to do until the round's slowest
region finishes, because the next round's observations are not read yet, and they cannot be
read while builders are still reading the cache (§6.2: the organiser's draw moves the memory a
builder is reading).

Measured on the tomato25 cohort — 2,176 samples at about three reads a position — over the first
300 kb of chromosome 0 with 48 threads, refinement off, on 2026-11-13:

| | |
|---|---|
| calling | 1 min 45 s, 20,293 records |
| time inside rounds | about the whole of calling |
| a round's regions, average | 3.6 s each |
| a round's slowest region | about 6 s |
| threads running, sampled during calling | about 29 of 48 |
| reading the round's stretch (step 2) | about 9% of calling |

So **about a third of the pool sits idle at the end of every round**, and step 2 runs while
nothing else does.

**Making regions smaller or more numerous does not fix it.** The 500-base region is already the
floor (`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`), and giving a round more regions than threads
— two or four per thread, over a correspondingly longer stretch — left calling no faster (1 min
47 s at one per thread, 1 min 48 s at two, 2 min 03 s at four) while peak memory rose from
22 GB to 34 and 58 GB, because a round holds every sample's observations over its whole stretch.

**Why §6.2's earlier verdict does not carry over.** §6.2 measured the round's overlap on 2026-08-18
and judged it not worth building, for two reasons. Both depended on the setting then measured:

- *The handover cost more than the overlap saved* (19–23% against 12–14%) — because the
  window was handed to the builders as an owned copy. The copy is what the proposal below removes.
- *The stage upstream of the merge was 14 to 23 times the merge*, so the merge's own speed hardly
  moved a run. That was alignment mode at 16 to 63 samples, where walking the alignments
  dominates. **In psp mode there is no upstream stage left in the calling pass**: the psps were
  written beforehand, and the merge with its calling is the whole pass. At 2,176 samples that
  pass is projected at about 64 hours for the 783 Mb genome (300 kb in 88 s, scaled).

## 2. The proposal

**Let the next round's regions start while the current round's last regions finish.** Three
changes, each needed for the next:

1. **A per-sample store whose records do not move.** Today each sample's held observations are
   one `Vec`; evicting shifts what survives to the front, and drawing may reallocate, so a
   builder holding a slice and the organiser drawing cannot coexist. Instead, each sample's
   records go into fixed-size **chunks** that are appended and never moved; eviction drops whole
   chunks that every unfinished region is past. A builder holds ranges into chunks; the organiser
   appends further on. This is the "append-only storage whose addresses do not move" that §6.2
   named as the version that could pay.
2. **Cover beside the builders.** With records that do not move, the organiser reads round
   *n + 1*'s stretch while round *n*'s regions are being built — step 2 stops being serial.
3. **Hand regions out as a stream, not as rounds.** Regions go to the pool in genome order as
   soon as their stretch is covered; a thread that finishes takes the next. Eviction follows the
   oldest unfinished region instead of the round's start. The organiser still puts results back
   in genome order before handing them over, as it does now across rounds.

**What it should buy — replayed, not guessed (§4).** On 1 Mb of the same cohort the rounds took
371 s. Replaying every region's measured build time through a simulated 48-thread pool:

| arrangement | seconds |
|---|---|
| today's rounds, as measured | 371 |
| today's rounds, replayed (cover, then regions, then a barrier) | 277 |
| a stream of regions, cover overlapped | 118 |
| every region's work spread perfectly over 48 threads | 116 |

The average region took 2.75 s and each round's slowest 5.76 s, so a round keeps the pool about
42% busy; a stream comes within 2% of perfect balance. The 94 s between the first two rows is the
hand-over of the previous round's records (the fragment keeper and record preparation, which run
on the pool since `4111e859`): its tasks queue behind the round's regions, run in the round's tail
and spill past its end. In a stream they are more tasks for the pool. With their CPU added (about
230 s over these 58,355 records) the stream predicts about **125 s against 371 s — calling about
three times faster** at this cohort size, so roughly 64 hours to 21 for the genome.

**What it costs.**

- *Memory.* The cache holds the stretch of every region in flight plus the next one being
  covered — roughly 1.5 to 2 rounds of ground instead of one. Doubling a round's ground took
  peak memory from 22 GB to 34 GB at 2,176 samples (§1), so the estimate is about 35 GB, under
  the 80 GB the machine allows. It must be measured, not assumed.
- *Code.* `observation_cache.rs` (4,100 lines) changes its storage; `rounds.rs` (500) becomes a
  stream of regions; the builder's interface to the cache (`with_observations`, `build_at`,
  `window_starts_into`) keeps its shape but returns ranges into chunks. The organiser's
  ordering and overlap rules (§6.1, §6.3) do not change. The existing driver equivalence tests
  — the round driver writes what the streaming driver writes — are the main safety net.

## 3. How it behaves across the range

- **One sample.** A region is cheap and the pool has little to wait for; the stream behaves like
  today's rounds. The chunked store adds a little indexing and nothing else.
- **Tens of samples.** Regions are fast and uniform; the gain is the overlapped cover, a few
  percent.
- **Thousands of samples.** The case this is for: regions are slow and uneven, and the idle tail
  is a third of the pool.
- **Deep coverage.** Regions cost more and vary more, so the tail grows and so does the gain. The
  memory term grows with depth too, which is why §4 measures memory as well as time.
- **Alignment mode.** Walking the alignments still dominates there (§6.2's finding stands), so
  the gain is small; the change must not make it slower, which the measurement checks.

## 4. Before building all of it: a measurement that can refute it

Done on 2026-11-13 over 1 Mb (`tmp/replay/replay.py` replaying `REGION`/`COVER`/`ROUND` timings
from an instrumented run); the result is the table in §2. The method, for repeating it on other
ground or another cohort: **Log each region's build time from a real run, then replay the log through a
simulated pool of 48 threads that takes regions in genome order with no barrier**, the cover
charged as it would be when overlapped. The replay predicts what streaming would give on that
ground; comparing it with the measured round time splits the idle third into the part only the
stream recovers (the slowest region's tail) and the part the overlapped cover recovers.

If the replay predicts under 15% on calling, stop: the chunked store is not worth its code. If
the cover's share is most of the gain, build steps 1 and 2 and leave the rounds. On tomato25 the
cover was 35 s of the 277 and the tail most of the rest, so the stream is what pays.

**One caveat the replay cannot settle.** Region times were measured with the pool 42% busy; with
it near 100% busy, memory bandwidth is shared by more threads at once, and each region may take
longer. The replay is therefore an upper bound on the gain, and the built version must be timed.

## 5. Decisions for the owner

1. **Build it?** Recommended: yes, in the order of §2 — chunked store, then cover beside the
   builders (measured on its own), then the stream — stopping after any step whose measurement
   says the next will not pay.
2. **The memory allowance.** Recommended: the cache may hold up to two rounds' ground; a run that
   would exceed the machine's limit says so up front, as the open-files check does today.
3. **§6.2.** If built, §6.2's paragraph on the round is rewritten to say the round was replaced
   and why, keeping the 2026-08-18 measurement as history.

## 6. Where the build departed from the proposal

- **No chunked store.** §2's step 1 proposed per-sample storage whose addresses do not move. What
  was built instead moves each round's records into a sealed view and copies back only what
  reaches past the round — the 2026-08-18 handover `cohort_merge.md` §6.2 rejected — with the one
  change that makes it pay: the evidence behind a psp record stays in sealed stretches of the
  source's arena, shared by every view that needs it, so the copy is summaries and coverage
  windows only. The 4,100-line cache kept its storage.
- **Many small rounds, not two.** With two rounds held the next read could start only when the
  older round was built, and the pool ran out of work at every read: over 1 Mb, 2 × 12,000 bases
  called in 2 min 53 s at 37.3 GB. Six rounds of 4,000 called in 2 min 24 s at 32.4 GB; eight and
  twelve gained nothing more (2 min 28 s, 2 min 34 s).
- **The reads on a pool of their own.** As a task on the regions' pool a read waited behind
  region builds at each of its parallel steps and took 4 to 6 s a round, against 0.6 s for the
  round driver's; on 16 threads of its own it takes about 1 s and is off the critical path.
- **Regions dispatched on a thread of their own.** Done on the calling thread, which also hands
  finished regions over through the record pipeline, the pool went without regions while the
  pipeline ran: the stream was slower than the rounds (1 min 27 s against 57 s over 300 kb) until
  the two were split.
- **Memory came out the same to a third above the round driver's** for the same ground — the
  same over 300 kb, 33.3 GB against 25.8 over 1 Mb — not the 1.5 to 2 rounds §2 estimated,
  because the ground is cut into six rounds rather than held twice.
