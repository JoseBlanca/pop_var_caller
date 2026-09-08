# `generate-psps`: each CRAM container read once, decoded a step early, and the psp compressed off the walking thread

**Date:** 2026-09-08
**Branch:** `main`, commits `a2a77780`, `aa598441`, `c9018a70`, `a6457dbb`, `5548c5fb` and
`9969cc65` on `8203d218`; then `2c582c68` and `02ffa51c` (§7a, §7b)
**Acting on:** the BAM/CRAM → psp performance review of 2026-09-08

**One tomato accession over 10 Mb of `SL4.0ch01`, out of a 49 GB whole-genome CRAM, now takes
28.17 s where it took 57.03 s** — **2.02× the throughput**, with the psp and the census
byte-identical outside their own timestamp, and **30 MB more resident**, about 295 MB against
264.6, which is the one number that went the wrong way. Six things carry it: the walk's two
readers were decoding every container of the file separately; the walk was compressing every psp
block itself; it was inflating each container in the middle of building a locus rather than a step
ahead, on a thread of its own; every admitted read was scanning the whole pending-mates map to
find the stale entries in it; a deletion anywhere along a read sent every base that read covered
down the general path; and every covered base sorted its column's chain ids that the active set
could have kept in order.

**The first four are §1 to §6c and take it to 35.59 s. The last two are §7a and §7b, and take it
from 34.25 s to 28.17 s** — the two figures differ because the same binary was re-measured in a
later sitting.

**The fixture is a 103.5× sample**, measured with `samtools coverage` on 100 kb at
SL4.0ch01:1.0–1.1 Mb — the high end of the depth range this caller commits to, not the tomato
cohort's 3×. Every share quoted below is a fact about that corner.

The command:

    ./target/release/pop_var_caller_exp generate-psps \
        --reference ~/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa \
        --catalog benchmarks/tomato1/crams/S_lycopersicum_chromosomes.4.00.repeats.parquet \
        --alignment benchmarks/tomato_big_cram/DRR000741.p1.cram \
        --regions tmp/psp_perf_chr1_10mb.bed --output-dir <dir> --force

Host: macOS 15, Apple Silicon, 18 logical cores, 64 GB. Native host release build — fat LTO, one
codegen unit, `-C target-cpu=apple-m1`, mimalloc.

---

## 1. What was doubled, and why nobody had seen it

A walk drives two locus generators over one sample: the SNP/indel one and the repeat-tract one.
Each holds its own [`SampleCursor`](../../../../src/ng/read/input/sample_cursor.rs), and
therefore its own [`CramAlignedReadsReader`](../../../../src/ng/read/input/aligned_reads_reader/cram.rs),
on purpose — [`walker.rs`](../../../../src/ng/run/walker.rs)'s `generic_path_generators` says why:
*"their regions interleave, and sharing one would tie their lifetimes together"*.

**Interleaving is exactly the case where both readers reach the same container**, and each was
inflating it, decoding its ten thousand records and fetching its reference window separately. The
review counted it with an instrument in `decode_container_at` that recorded which reader pulled
which offset:

| | 2 Mb of `SL4.0ch01` | 10 Mb of `SL4.0ch01` |
|---|---:|---:|
| container decodes | 498 | 2,534 |
| distinct containers | **249** | **1,267** |
| containers decoded by both readers | **249 — all of them** | **1,267 — all of them** |
| decodes per reader | 249 and 249 | 1,267 and 1,267 |
| seconds inside `decode_container_at` | 3.23 s of 14.13 s | 16.54 s of 56.89 s |

Neither reader ever decoded a container twice on its own. **The whole doubling was across the two
cursors**, and it was invisible from either reader's own accounting, which is why
`cram_read_path_2026-09-04.md` §7 could say ng *"already keeps exactly one decoded container per
open file"* and be wrong: one per *cursor*, and there are two.

## 2. The fix: the decoded containers belong to the open file

[`DecodedContainerCache`](../../../../src/ng/read/input/aligned_reads_reader/container.rs) sits on
`AlignmentFile` and both of that file's readers take from it.

**Two slots, and the number came from the measurement rather than from taste.** Of the 249 repeats
over 2 Mb, 245 were the very next decode and all 249 fell within two; over 10 Mb, 1,246 of 1,267
fell within two, and the worst gap on either fixture was four. So two slots catch 98% of the
repeats — and they hold **what the two readers were holding anyway**, one decoded container each.
Four slots would catch the last 2% for about 7 MB more per open file. The cheap option is also the
one that costs no memory, and for the stretch where both readers sit on one container it holds one
where they held two.

**`Arc<Mutex<..>>` rather than `Rc<RefCell<..>>`, decided on evidence.** The review flagged this as
the question to answer rather than assume. `AlignmentFile` is shared behind an `Arc` and its
cursors are minted from several threads at once by
`cursors_on_one_file_read_the_same_thing_from_many_threads`
([`open_bam.rs`](../../../../src/ng/read/input/open_bam.rs)), so anything on it must be `Sync` — a
`RefCell` cannot go there whatever the walk does today on one thread. The lock is taken twice per
container, to look and to store, and **never across the decode itself**, so two readers that ever
did run at once would at worst both decode one container rather than wait on each other.

**One cache per open file is one cache per walk**, which is the scope the sharing wants.
`SampleReads` is not `Clone` and opens its own `AlignmentFile` per sample, so nothing shares this
between two samples or two workers; a walk split across workers (§11 question 3) would give each
its own, which is what a cache of *the container a reader is on* needs.

**Two things the change could have lost quietly**, both named by the review and both kept: a cache
hit is still charged to the second reader's `other_sample_records`, which is a per-reader tally
that would otherwise silently halve; and `last_decoded_offset`, the guard against serving a
multi-slice container's records twice, stays per reader.

**The reference window goes with it.** The window each container is decoded against is fetched
inside `decode_container_at`, so a cache hit skips that fetch too — the review asked whether the
change removed the doubled fetch, and it does.

## 3. What it is worth

`/usr/bin/time -l`, `RAYON_NUM_THREADS=1`, arms interleaved in one sitting as one binary each
copied to a fixed path and writing into one output directory, so the command line is
byte-identical between arms.

| | | before | after step 1 | and after the four smaller changes |
|---|---|---:|---:|---:|
| **2 Mb**, 4 pairs | instructions | 309,072 G | 286,257 G | **276,706 G** |
| | cycles | 62,646 G | 55,573 G | **54,170 G** |
| | wall | 14.43 s | 12.85 s | **12.26 s** |
| | peak resident | 204.70 MB | 196.28 MB | **190.92 MB** |
| **10 Mb**, 3 pairs | instructions | 1,205,099 G | 1,088,779 G | **1,037,554 G** |
| | wall | 58.27 s | 49.97 s | **47.28 s** |
| | peak resident | 272.30 MB | 233.96 MB | **241.12 MB** |

The baseline's own instruction spread over four runs is 12.8 M on 309 G — one part in 24,000 — so
the container-sharing effect is 1,800 times the noise.

**And with the rayon pool left alone, which is how the command actually runs**, three interleaved
pairs at 10 Mb: **57.03 s → 47.20 s** and **264.6 MB → 245.7 MB** for these two commits, then
**43.99 s and 249.0 MB** with §6a's compression offload, then **36.75 s and 292.4 MB** with §6b's
container prefetcher. The pinned figures above are what the arms are *compared* on; these are what
the command costs.

The 10 Mb peak rises 7 MB between the second and third columns where the 2 Mb peak falls 5 MB.
Both are mimalloc holding different segments for a different allocation mix, both reproduce across
three runs, and both sit well inside the 272.30 MB the command started at.

## 4. The four smaller changes, and the two that were refused

Each had been built and measured by the review and each was too small to lead with. Together they
are −3.34% of retired instructions on top of the container cache. Judged by whether they also
improve the code:

- **The per-base mismatch-fraction loop** (`src/bam/alignment_input.rs`, −2.08%) indexed `seq`,
  `ref_seq` and `qual` through `get(read_pos + k).unwrap_or(..)` — one bounds test and one branch
  each, per aligned base, three times over. The three slices are now taken once per CIGAR run and
  iterated together, branchlessly. **Kept**: it also makes the two base-set tests symmetric, where
  the old loop tested the read against four bytes and the reference against eight.
- **A third of every CRAM slice's blocks were inflated for a caller that reads no tags**
  (`vendor/noodles-cram`, −1.14%) — 780 of 2,254 blocks and 14.8 MB of 111.1 MB per 60 containers.
  FORK.md gains change 7. **Kept**: `tags_can_be_left_unread` becomes a wrapper over the function
  that returns the set, so there is one predicate where there could have been two that drift.
- **Two guards that state what the data is** (−0.37% and −0.34%): `canonicalise_runs` returns a
  single run unsorted, which is what a DNA-seq read gives at almost every position;
  `find_overlapping` answers an empty table without two `partition_point` searches, and the
  ordinary-column lane asks it at every covered base. **Kept**: both say something true about the
  data that the code did not say.
- **A census locus the selection kept nothing in no longer builds a depth array.** Every write in
  `add_generic` sits behind `generic_index`, which answers `None` at every position of such a
  locus, so the array was filled from every observation and thrown away — 4 loci in 5 here, and
  399 in 400 on a whole-genome walk, the census budget being the same two million positions over
  four hundred times the ground. **Kept**: the predicate is exact rather than conservative, and it
  is one binary search.

**Refused, both from the same review:**

- **The DP row hoist** in `ssr_unit_robust.rs` earns −0.13% for about thirty-five lines of
  `split_at_mut` plus `split_first_mut` plus an index-remapping closure, in place of
  `rows[prev_slot][column]`. That is a little speed for code that is harder to read, and the
  review's own note says the cost that remains in that loop is arithmetic and not addressing.
- **Turning the census's `BTreeMap<ReadGroupId, Vec<u32>>` into a `Vec<Vec<u32>>` searched
  linearly.** The early exit above already removes 399 of every 400 of those maps, so what the
  rework buys is allocations that no longer happen, and what it costs is a slot indirection and
  two linear scans per observation.

## 5. A bug in one of the patches, found by the test written to check its argument

The mismatch-fraction rewrite rests on two arguments rather than on a transcription: that a
position past the end of the read or of the reference reaches neither counter, so the loop may
stop there; and that `& 0xdf` accepts exactly the eight reference bytes the two-case list did. An
argument is what a randomised comparison is for, so
`the_rewritten_mismatch_loop_agrees_with_the_one_it_replaced` runs 20,000 random reads against the
loop as it was, kept as the oracle.

**It panicked.** With the length clamped but not the start, a malformed
record whose CIGAR claims more read bases than it carries leaves `read_pos` past the end of `seq`,
and Rust rejects `&seq[read_pos..read_pos]` there — an empty range whose *start* is out of bounds
is still out of bounds. So a record the old loop merely counted nothing for would have stopped the
run. The starts are clamped too, which is an identity wherever the usable length is non-zero.

This matters more than the usual because **the function is shared with the frozen production
caller** (`src/pileup/per_sample/read_processor.rs`, `src/ng/read/filtering.rs`). The generator
aims at what the arguments turn on: sequences, qualities and reference deliberately allowed to run
short of what the CIGAR claims; bases drawn from `ACGTNacgtn` so both cases and the excluded byte
appear; and a quality floor of 0 among the four, so the short-quality tail's *does an absent score
clear the floor* is exercised in both directions.

## 6. The BAM arm doubles too, and this does not fix it

Nobody had looked. Instrumented over the same 2 Mb converted to BAM with `samtools view -b`, the
two readers served **2,485,552 and 2,483,893 records for a region holding 2,485,551** — each of
them reads every record of the region.

Sized from a 10 s `sample` of that run, out of 5,962 busy samples: bgzf inflation 17.0%, BAM
record decoding 4.5%, the mismatch filter 2.1%. **So the second pass is about 12% of the walking
thread there**, which is the same size as the CRAM arm's was at 2 Mb.

**There is no equivalent fix.** The CRAM cache reuses an object that already exists — a decoded
container — and a BAM has no such object: records are decoded one at a time straight into the
caller's buffer. Reuse would mean either a bgzf block cache inside noodles, which is a 64 kB
granularity instead of 3.6 MB and would leave the record decoding doubled anyway, or one reader
serving both cursors, which is the arrangement `walker.rs` declines for a reason. Recorded here so
it is not re-derived.

## 6a. Compressing the psp's blocks on a thread of their own

`push` closed a block, decoded its head, compressed it at zstd level 9 and wrote it, all inline —
and **nothing downstream of a block waits for its compressed bytes**. No record after it depends
on them, and its index entry is not built until the write returns. The only thing tying
compression to the walk is the *order* blocks reach the file in, which one compressing thread with
FIFO queues keeps by construction; and `BlockCompressor::compress` is a pure function of one
payload, so which thread ran it cannot be read out of the file.

**Measured on wall, not on retired instructions, and that is the point.** The change moves work
rather than removing it, so the instruction counter — which counts the whole process — cannot see
it. Instructions move +0.02% and user time +0.2%, which is the copy the offload adds; wall is what
it buys.

| | | before | after |
|---|---|---:|---:|
| **2 Mb**, 5 pairs, pinned | wall | 12.57 s | **11.92 s** (−5.2%, won 5/5 by 0.65–0.73 s) |
| | instructions | 276,719 G | 276,748 G (+0.01%) |
| | user | 12.31 s | 12.34 s (+0.24%) |
| | peak resident | 190.92 MB | 206.06 MB (+15.1 MB) |
| **10 Mb**, 3 pairs, pinned | wall | 48.10 s | **44.67 s** (−7.1%, won 3/3 by 3.20–3.43 s) |
| | peak resident | 241.12 MB | 235.47 MB (−5.6 MB) |
| **10 Mb**, 3 pairs, pool left alone | wall | 47.44 s | **43.99 s** (−7.3%, won 3/3 by 2.84–3.45 s) |
| | peak resident | 245.7 MB | 249.0 MB (+3.3 MB) |

**The memory cost is small and does not keep a sign**: +15 MB at 2 Mb and within a few megabytes
either way at 10 Mb. It is two block payloads in flight, and a payload grows with read depth —
which this fixture already exercises near the top of the range at 103.5×.

**A profile says directly that the compression left the walking thread.** Taken 6 s into a 10 Mb
walk: **47 ZSTD frames on the `psp-block-compression` thread and none on the main one**, where
before they were 1,631 of that thread's 22,480 busy samples. The walking thread parks zero times,
so it never waits on the compressor.

### No knob, which is a change from what the review recommended

It proposed one because a cohort is parallelised by running invocations and a saturated machine
cannot use the extra thread. Measured rather than assumed — **eight concurrent invocations over
disjoint 2 Mb windows, four pairs**:

| | wall for the whole set | sum of the eight peaks |
|---|---|---:|
| before | 15.65, 16.71, 15.82, 16.30 s | 1,701 MB |
| after | 16.19, 16.07, 15.50, 16.74 s | 1,796 MB |

**The two arms cannot be separated on wall** — the variant is ahead in two pairs of four, and each
arm's own spread is seven times the difference of the medians. What it costs on a saturated
machine is memory, and it is **about +12 MB a sample on a 240 MB footprint**. A flag, its
documentation and a second code path are not worth 5% of one sample's memory when the wall is
unchanged; the thread is blocked on a channel rather than spinning, and it is created and joined
per psp, so a cohort walked in one process accumulates none.

### Why it went before splitting the walk

Under a walk split across workers (§11 question 3) the workers feed **one serial merger that owns
the psp writer and the census** — one writer, because the block cut follows the coordinate grid
rather than the workers, which is exactly what
`one_sample_gathered_at_any_worker_count_gives_byte_identical_files` guarantees. So writing is that
design's serial floor:

| the writer's share of the walking thread | 4 workers | 8 workers | however many |
|---|---:|---:|---:|
| 14.2% — compressing 7.3%, encoding 5.9%, census 1.0% | 2.8× | 4.0× | 7.0× |
| **7.4% — encoding 6.4%, census 1.1%** | **3.3×** | **5.3×** | **13.5×** |

**This roughly doubles what splitting the walk is worth**, which is the reason it went first
rather than sitting with the other small wins. (§7 restates the same cap in seconds, which is the
form that survives the walking thread getting shorter again.) What is left on that thread is
encoding records into block payloads — order-dependent *within* a block and independent *between*
blocks — so if four workers is not enough, building a whole block off-thread is the next place to
look, and it is an encoding question rather than a scheduling one.

### Two invariants, and one behaviour that changed

The result channel is **unbounded**, or the two sides can wait on each other; and the line is
**joined on drop**, so no thread outlives its psp.
`writers_abandoned_with_blocks_in_flight_neither_hang_nor_pile_up` covers the path nothing else
reached — the existing drop tests push one record, so no block is ever in flight — by abandoning
twenty writers in turn, each with ten blocks closed, in one process. Mutation-verified: a drop that
joins without first closing the channel hangs it, and 90 seconds was not enough for it to finish.

**A compression or write failure is now raised a few records after the record that caused it**,
because handing a block over is not writing it. Nothing a reader can see changes — the file stops
in the same place and `finish` still refuses to seal it — and the line hands back the sentence
`spent` records, so the message still says which of the two things went wrong.

## 6b. Decoding the containers a step ahead of the walk

The decode does not sit *beside* locus generation, it sits **inside** it. The generator asks for
the next read; that call goes down through the cursor to the CRAM reader, and when the reader has
run out of records it inflates and decodes a whole container — 3.6 MB, ten thousand records —
before returning the one read that was asked for. The call tree, as samples on the walking thread:

```
16,260  GeneratorSet::next_locus              "give me the next position"
 └ 4,236  PreparedSampleReads::next             … which needs another read
   └ 3,745  AlignmentCursor::next_read
     └ 2,070  RegionRawAlignedReads::read_next
       └ 2,048  Slice::decode_blocks_inner      ← a whole container, decoded here
```

Nothing about building a position needs the *next* container, and which container is next is
knowable — the `.crai` lists them in file order and the walk goes forward through it. So a thread
now decodes it into the cache the two readers already share.

**A hint, not a queue.** The prefetcher holds one wanted offset. A reader that takes a container
overwrites it with the next offset it expects to want, and a hint the thread never picked up is
simply replaced. Work the walk has already passed is worthless, and one cell can hold nothing
stale.

**One decode per container even when two threads want it at once, which is the part that had to
be right.** The prefetcher makes real a race that was only theoretical: the walk can reach a
container the prefetch thread is midway through, and left alone both would decode it — §1's
doubling arriving by a new route, worst exactly when the prefetcher is behind. So a slot is
claimed before its decode starts, a second thread waits on a condition variable rather than
starting its own, and a decode that fails clears its claim and wakes the waiters.

| | | before | after |
|---|---|---:|---:|
| **2 Mb**, 3 pairs, pinned | wall | 11.85 s | **10.46 s** (−11.7%, won 3/3) |
| | instructions | 276,764 G | 275,561 G (−0.43%) |
| | peak resident | 206.1 MB | 243.0 MB (+36.9 MB) |
| **10 Mb**, 3 pairs, pinned | wall | 44.53 s | **37.47 s** (−15.9%, won 3/3 by 6.98–7.24 s) |
| | instructions | 1,037,821 G | 1,030,059 G (−0.75%) |
| | peak resident | 235.4 MB | 282.1 MB (+46.7 MB) |
| **10 Mb**, 3 pairs, pool left alone | wall | 44.23 s | **36.75 s** (−16.9%, won 3/3 by 7.00–7.53 s) |
| | peak resident | 249.0 MB | 292.4 MB (+43.4 MB) |

**Instructions fall rather than rise, which is the evidence that the sharing survived.** A
prefetcher that raced the readers would decode containers twice and the count would go up; it goes
down 0.75%, because the thread also catches containers a reader would have re-decoded after a
region jump.

### The thread is the whole win, and the extra cache slots are none of it

The change bundles two things — four cache slots instead of two, and the thread — so one binary
settled the split by an env-var ablation rather than by argument:

| | wall, 10 Mb | peak resident |
|---|---:|---:|
| two slots, no thread | 44.52 s | 235.4 MB |
| **four slots, no thread** | **45.15 s** | **251.8 MB** |
| four slots and the thread | 37.84 s | 283.0 MB |

**Four slots on their own buy nothing and cost 16 MB.** They are not there to be faster; they are
there so a prefetched container does not give up one a reader is about to want. Three slots was
tried and is not cheaper: the same wall, 0.15% more instructions, and a peak that swings wider.

### It wins on a saturated machine too, which the compression offload did not

Eight concurrent invocations over disjoint 2 Mb windows, three pairs of one binary:

| | wall for the whole set | sum of the eight peaks |
|---|---|---:|
| without the prefetcher | 15.45, 15.92, 16.10 s | 1,829 MB |
| with it | 13.83, 14.11, 14.07 s | 1,982 MB |

**−11.6%, ahead in all three pairs**, for about +20 MB a sample. The prefetch thread is mostly
asleep, so eight of them do not cost eight cores.

### How well it keeps up, and what that kills

A profile taken 6 s into a 10 Mb walk, by thread:

| thread | busy |
|---|---:|
| the walk | **99.8%** — it waits 47 samples of 20,502 |
| `cram-container-prefetch` | 22.6% |
| `psp-block-compression` | 10.2% |

**The walking thread is essentially never blocked on the prefetcher**, which is the question the
review said to ask and the answer that decides the rest. Because it is not, the expensive half of
this candidate is dead: moving the read *preparation* off as well — two prefetchers, two queues of
prepared reads, a prefetcher that must be stoppable and rewindable — is now worth **2.3% of the
walking thread, about 0.85 s of a 37 s run**, because what was left of that share went with the
decode.

### What it costs

**+43 MB on a single-sample run**: two more cache slots (7.2 MB), a third CRAM descriptor and a
third reference reader — both are file positions rather than values, so neither can be shared with
a reader that is seeking elsewhere — and up to three decodes in flight, which the claim rule
refuses to evict. Against the 264.6 MB this command cost before any of today's work it is now
292.4 MB, and that is the one number that has gone the wrong way.

## 6c. Making the walk legible, and the scan that was hiding in it

**The largest single leaf in the profile was a trait forwarder.** `generator.rs:1421` is
`PileupGenerator::next_locus(self, reads)` — one call, no body — and it carried **44.6% of the
walking thread**. Fat LTO with one codegen unit had inlined the whole generator into it, and
`sample` aggregates that node across every address inside it, so 15 seconds of a 37-second run
could not be attributed to anything.

**Temporary `#[inline(never)] // PROFILING SCRATCH` markers on the walker's phases broke it
open**, in the same build configuration so nothing else moved. The marked build ran at 36.69 s
against the unmarked 37.47 s and 1.2% fewer instructions, so the barriers cost nothing and the
profile is representative. They are a measuring instrument and are not in the tree.

What it showed, as shares of the walking thread:

| | before the fix below | after |
|---|---:|---:|
| `WalkerState::process_position` | 16.3% | **22.0%** |
| `fast_column::try_ordinary_column` | 13.6% | 13.8% |
| `ChainIdAllocator::allocate_for_read` | **11.8%** | below 0.5% |
| sorting chain ids (`Vec<u64>`) | 11.2% | 11.5% |
| `ssr::classify::delimit`, the tract aligner | 9.9% | 10.0% |
| `apply_events_into` | 3.1% | 5.9% |
| `finalise_recycling` | 3.8% | 4.2% |

**`allocate_for_read` was the surprise, and a second round of markers said why**: 2,333 of its
2,424 samples were in `evict_stale_pending` and 60 in its own body. That function ran an
`AHashMap::retain` over the **whole** pending-mates map on every admitted read, to drop the first
mates the walk had passed the 10 kb lookup window of. The map holds every first mate still inside
that window — thousands of entries on a 104× sample — so it was a full scan of the map about
twelve million times over 10 Mb.

**The fix is a queue in registration order**, so eviction pops instead of scanning. Reads arrive
sorted by alignment start (`admit_read` refuses one that does not) and an entry's `seen_at` *is*
its read's alignment start, so the queue is non-decreasing by construction: everything stale is at
the front, and the first entry inside the window ends the loop. The set it drops is exactly the
set the scan dropped, which matters beyond tidiness — `mate_lookup_evictions` is asserted against
production's walker in `parity.rs`.

**A pending entry now carries a serial, and it is not bookkeeping.** The queue names entries by
qname, and a file may register one qname twice — a completed pair, then another read with the same
name. Without the serial the older slot would find the *newer* map entry under that name and drop
it; its second mate would mint a fresh chain id instead of pairing, and the psp would differ, on
input that is merely malformed rather than unreadable. The scan could not get this wrong, because
it read each entry's own `seen_at`. `a_qname_registered_twice_keeps_the_second_registration` pins
it, mutation-verified.

| | | before | after |
|---|---|---:|---:|
| **10 Mb**, 3 pairs, pinned | wall | 36.95 s | **35.70 s** (−3.4%, won 3/3) |
| | instructions | 1,030,041 G | **967,422 G** (−6.1%) |
| | cycles | 208,908 G | 202,888 G (−2.9%) |
| **10 Mb**, 3 pairs, pool left alone | wall | 37.30 s | **35.59 s** (−4.6%, won 3/3) |

**Instructions fall twice as far as wall**, which is the shape a linear sweep of a few thousand map
entries has — high instructions per cycle, few stalls. And what it removes is not only the 1.7 s
but the *scaling*: the scan cost reads × pending entries, so it grew with read depth and with the
mate-lookup window together. At 300× it would have been proportionally worse.

## 7. Where the wall goes now, and what is left

**Two more commits landed after this section was first written** — `2c582c68` and `02ffa51c`,
§7a and §7b below — and they take the same 10 Mb run from 34.25 s to **28.17 s**. The table
immediately below is the split *before* those two, and is kept because §7a and §7b are written
against it; the split after them is in §7c.

**The walk's one thread is locus generation and little else.** Profiled after §1 to §6c, with
the inline barriers of §6c in place so the parts are separable:

| | share of the walking thread | ~seconds |
|---|---:|---:|
| `WalkerState::process_position` | **22.0%** | 7.8 |
| `fast_column::try_ordinary_column` | 13.8% | 4.9 |
| sorting, all sites — chain ids are 11.5% of it | **14.9%** | 5.3 |
| `ssr::classify::delimit`, the tract aligner | 10.0% | 3.6 |
| `apply_events_into` | 5.9% | 2.1 |
| encoding psp records | 7.3% | 2.6 |
| `finalise_recycling` | 4.2% | 1.5 |
| `memmove` | 4.6% | 1.6 |
| the census | 1.3% | 0.5 |
| reading and filtering reads | 2.3% | 0.8 |
| decoding containers, compressing psp blocks | on threads of their own | |

Splitting the walk across k workers (§11 question 3) is what is left, and **its cap is best stated
in seconds, because a share moves when the thread it is a share of gets shorter**. A 10 Mb run of
this sample decomposes as:

| | seconds | threads |
|---|---:|---|
| reading and MD5-ing the reference, the catalog, the segmentation | **2.6** | serial, before the walk begins |
| locus generation, the read cursor and its filters | 31.2 | k ways |
| encoding psp records, and the census | **3.0** | serial, on the merger |
| decoding containers, compressing psp blocks | — | already off the walking thread |

which is `2.6 + 31.2/k + 3.0`: **13.4 s at four workers, 9.5 s at eight, and 5.6 s however many**
— against the 36.75 s the run cost at that point. (§7c restates it against 28.17 s.) The emulation of that split as concurrent processes over disjoint BEDs
reached 2.97× at four ways and 3.88× at eight, before any of today's changes; it is an upper bound
on the worker side, since processes share nothing and each wrote its own psp.

**The two serial terms are nearly equal and neither of them is the walk.** Past four workers the
thing to attack is one of them, and they need different fixes: the setup is re-reading and
re-digesting one unchanging reference on every invocation — paid 63 times across a cohort for a
file that did not change — and the merger's 3.0 s is encoding records into block payloads, which
is order-dependent *within* a block and independent *between* blocks, so a whole block could be
built off-thread the way its compression already is.

### The largest thing left, measured rather than guessed

`process_position`'s 22.0% splits, under a second round of inline barriers, into the **general
fold at 16.0%**, its own body at 4.9% and `CigarCursor::events_at` at 3.4%;
`may_have_mate_overlap_at` is below 0.45%. The fold is what runs at a column the ordinary-column
lane hands back — and it handles **about two columns in ten while costing 16% of the thread**.
Per column the general path is about **8.7×** the fast lane.

So the lever is fast-lane coverage, and a counter probe over the whole 10 Mb walk says exactly
what stops it. Of **10,641,693 columns**:

| | columns | share |
|---|---:|---:|
| taken by the ordinary-column lane | 7,603,904 | **71.5%** |
| handed back — **some active read's CIGAR contains an `I` or a `D`** | 2,031,990 | **19.1%** |
| handed back — two contributors share a chain id (mate overlap) | 970,978 | 9.1% |
| handed back — a record is already open over this base | 25,275 | 0.2% |
| handed back — depth over the column cap, or no active read | 9,085 | 0.1% |
| handed back — no contributor | 461 | 0.0% |

**The dominant reason is a property of the read and not of the column.**
`CigarCursor::matches_only` asks whether the read's *whole* CIGAR is free of indels, so one
indel-carrying read poisons every column it covers — about a read length of them, for every
column where it is active. Those 2.03 M columns are 67% of the general path's traffic and so
roughly **18% of the walking thread, about 6.5 s**; taken by the fast lane instead they would
cost about 3.4%, so the prize is around **5 s of a 35.6 s run**.

**Asked per-base instead, and the answer came back split (2026-09-08, `9969cc65`).** The
whole-read test bought two things: a fact about this base — no read here is doing anything but
showing a letter — and a fact about history — no event from an earlier base is still reaching in.
The first is now asked at the base, through `CigarCursor::plain_match_at`. The second was left a
whole-read question, because relaxing it too moved the psp from byte 33,403,518 of 157,822,123
and the case that diverged had not been identified. **It has been, and §7a is what it was.**

**Mate overlap at 9.1% is the second reason and it will grow with depth**, which the module's own
note says: at 300× a pair is present at most columns and the skip stops firing. This fixture is
104×.

**What was measured and refused**, in one place:

- **Moving the read preparation off the walking thread too** — two prefetchers, two queues of
  prepared reads, a prefetcher that must be stoppable and rewindable. Its ceiling fell three times
  before anything was built: 29% of the walking thread when the review priced it, 19.7% once the
  containers were shared, 17–18.5% after the smaller changes and the compression offload. Then
  §6b took the containers off, and what is left is **2.3%, about 0.85 s of a 37 s run**. Not worth
  its complexity.
- **Four cache slots without the prefetch thread** — no gain and +16 MB (§6b).
- **Three cache slots with it** — same wall, 0.15% more instructions, a wider peak.
- **The DP row hoist and the census map rework** — §4.

## 7a. The column that made the deletion relaxation unsafe, and it is one in 10.6 million

`9969cc65` shipped the insertion half of the fast lane's guard and left the deletion half alone,
because relaxing it moved the psp and nobody could say why. **The case is the first base of a
generic region that begins immediately after repeat ground**, and it is the only base of the walk
whose record is built from a window that does not start at it.

`window_of_read` widens that record's event query one base to the left, so that an insertion
anchored on the repeat's last base can be claimed by this region — and it then keeps **every
deletion the widened window returns, whatever its anchor**. So a read whose deleted run stops on
the repeat's last base contributes its deletion's quality proxy to `min_bq_for_read` at the
region's first base, and the ordinary-column lane, which asks only about this base, mints that
read's `ln ε` from the base's own quality instead.

**Found by diffing the loci the two arms emit, not the psp's bytes.** A temporary dump of every
locus the walk yields, from one binary with the guard and one without, over the same 15 kb: **one
line of 15,786 differs, and only in `q_sum`**. Over the whole 10 Mb it is **one locus of
10,641,693 columns** — `SL4.0ch01:2,185,561`, the first base after a tract, where 22 of the 76
reads carry a 10-base deletion ending on 2,185,560 and the general path's window is
`[2,185,560, 2,185,562)`.

**How often the general path does this at all, counted rather than inferred.** Of the 10,641,693
columns walked, **8,685 are a region's first base beside repeat ground** — 8 in 10,000 — and they
carry 796,886 read-into-record folds. **31** of those folds keep a deletion whose footprint stopped
before the record, at three columns; in **21** of them it lowers the read's minted error, by 1 to
10 Phred points. Only one of the three shows in the emitted psp when the lane's junction test is
removed, because at the other two the lane refuses the column on one of its other tests regardless
— checked by dumping every locus in a window around each and diffing the arms: 1 differing line of
15,786, then 0 of 10,360 and 0 of 13,541.

**The general path's behaviour there is left alone (owner, 2026-09-08)**, and `window_of_read`'s
own note now states the case rather than resting on the ordinary window's promise, which does not
cover it. Narrowing it would buy no speed: the lane must refuse that column regardless, because the
same widened window can hand the fold a claimed junction insertion.

So the lane refuses that column outright — one base per region — and asks about deletions per base
everywhere else. Refusing it also closes the **claimed junction insertion**, which the per-base
insertion test of `9969cc65` would otherwise have let through: a read carrying one shows a plain
letter at the region's first base, and the lane has no allele for the inserted bases the fold
spells before it.

| | before | after |
|---|---:|---:|
| instructions, 4 pairs pinned | 936,326 G | **861,122 G** (−8.03%) |
| its own spread | 95 G (0.010%) | 73 G (0.008%) |
| wall, 4 pairs pinned | 34.79 s | **30.68 s** (−11.8%, won 4/4 by 3.78–4.29 s) |
| wall, 3 pairs, pool left alone | 34.25 s | **30.23 s** (−11.7%, won 3/3) |
| peak resident, pool left alone | 294.8 MB | 293.9 MB |

## 7b. The chain ids come out ordered because the active set keeps them so

Every covered base built a fresh vector of the column's ~87 chain ids and sorted and deduped it,
because the emitted observation has to carry them ascending.

**Priced before anything was changed, and the pricing needed a control.** Two arms, one running an
*extra* copy-and-sort of the same list and one running the copy alone, both behind the same
`black_box`: **the sort is +5.72% of the run's instructions and +2.83 s of its wall**. The control
matters — the same barrier placed in `admit`/`expire_passed` was at first read as a 6% cost of the
index below, and it was the barrier.

**The gate the review set, because this could have cost more than it saved.** The 2026-09-02
review measured one extra linear pass over the same 87-element column at +0.43% instructions and
+3.7% wall, which is the same order as the sort. So the chain-id-ordered index was built in
`ActiveReads` **with the sort left in place and unused**, and measured against a control carrying
the same barriers in the same places: **+0.67% of instructions**, against a 2% stop threshold. It
cleared, so it was wired.

`ActiveReads` now keeps a second view of its own reads — `(chain_id, read_id)` ascending, one
entry per live read. `read_id` is in the key because two mates of a pair share one chain id; the
entry is **inserted** rather than appended, because a second mate takes the id its first mate was
given, minted up to a mate-lookup window earlier. The lane fills each observation's list by
walking that order: a subsequence of an ascending sequence is ascending.

Finding each read's observation without a search needs one dense table per column, indexed by
`read_id` less the column's smallest. The live set's ids are all but contiguous — reads leave in
very nearly the order they arrive — and on this fixture the span never reached twice the column's
depth. A span past 4,096, or a live set larger than the table's `u16` sentinel, falls back to the
sort; neither was reached.

| | before | after |
|---|---:|---:|
| instructions, 3 pairs pinned | 861,000 G | **841,608 G** (−2.25%) |
| its own spread | 120 G (0.014%) | 20 G (0.002%) |
| wall, 3 pairs pinned | 30.71 s | **28.85 s** (−6.1%, won 3/3 by 1.80–1.86 s) |
| wall, 3 pairs, pool left alone | 30.10 s | **28.17 s** (−6.4%, won 3/3 by 1.73–2.06 s) |
| peak resident, pool left alone | 295.9 MB | 295.9 MB |

**Instructions fall a quarter as far as wall**, which is the difference between a sort — data
dependent branches, few of them predictable — and two linear passes.

## 7c. What the walking thread holds after those two

Re-profiled the same way, temporary inline barriers in the same configuration. The marked build
ran 29.23 s against the unmarked 28.85 s and 2.5% more instructions, so the barriers are visible
but small and the profile is representative.

| | share of the walking thread |
|---|---:|
| `WalkerState::process_position`, whole | **48.8%** |
| — of which `fast_column::try_ordinary_column` | 31.6% |
| — of which `open_record::process_position`, the general fold | 8.4% |
| `ssr::classify::delimit`, the tract aligner | 12.4% |
| sorting, all remaining sites | 8.1% |
| `memmove` | 7.8% |
| `close_aged_records_into` (of which `finalise` 4.0%) | 4.4% |
| `expire_passed` | 4.2% |
| `admit_read` (of which `ActiveReads::admit` 2.3%) | 2.6% |

**The shares are not comparable term by term with the table at the top of §7**, which was taken
against a 35.6-second thread with different barriers; what is comparable is that sorting fell from
14.9% to 8.1% *while the thread itself got shorter*, and that none of the 8.1% is a chain-id sort
in the fast lane any more. What remains of it is the contributor and mate-overlap sorts under
`process_position` (3.1%), the lane's own `read_id` sort of its contribution buffer (1.8%), the
general path's `finalise` (1.3%), and the tract aligner and psp encoder (1.9%).

**The fast lane now takes the columns it was refusing.** `try_ordinary_column` is 31.6% of the
thread and the general fold 8.4%, where before §7a the fold was 16.0% against the lane's 13.8%.

**What is left, and it is question 3.** Neither commit touches the setup or the merger, so the
decomposition of §7 stands with its shareable term shortened: **2.6 s of serial setup, 22.6 s that
k workers can share, 3.0 s of encoding and census on one merger** — `2.6 + 22.6/k + 3.0`, so about
11.3 s at four workers, 8.4 s at eight, and 5.6 s however many, against 28.17 s today. The two
serial terms are now more than a fifth of the run between them, so the ratio the split can reach
has fallen with every one of today's changes; the seconds it can remove have not.

## 8. How it was checked

- **The psp and the census, with the command line held byte-identical.** The psp header records
  the process's own argv, so two runs from different binary paths or into different output
  directories differ from the header onward for a reason that has nothing to do with the change.
  Held identical, at both fixture sizes the psp is **the same length** and differs only in the 3
  to 5 bytes of its `created` timestamp; the census differs in the 16 bytes at offsets 1472–1487
  that two runs of one *unchanged* binary already differ in. That is the floor, and every change
  here sits on it.
- **`examples/ng_cram_decode_layers`**, which hashes every field of every decoded record, still
  gives `c0bdbb0e464a920b966ad487fc3ca678` over 60 containers and 600,000 records.
- **`cargo test --lib`: 6,662 passing**, eleven of them new. noodles-cram's own 223 and 75 pass.
  `diff -r` against the registry copy prints exactly FORK.md's list.
- **Mutation-verified.** A `take` that ignores the offset fails
  `two_cram_cursors_over_one_file_each_see_every_read`; one cache slot instead of two fails the
  policy test; the unclamped slice start fails the mismatch-loop comparison; a drop that joins the
  compressing thread without first closing its channel hangs
  `writers_abandoned_with_blocks_in_flight_neither_hang_nor_pile_up`; and one that joins the
  prefetch thread without telling it to stop hangs
  `a_cram_file_that_is_closed_leaves_no_prefetch_thread_behind` — 100 seconds was not enough for
  either to finish.
- **The CRAM decode was not re-run for §6a**, which touches only the psp writer. It was for §6b,
  and the digest is unmoved.
- **Every CRAM cursor test now runs with a prefetcher**, since the first CRAM cursor on a file
  starts one. `a_run_of_regions_through_one_cram_cursor_matches_a_linear_scan` and
  `two_cram_cursors_over_one_file_each_see_every_read` therefore compare a walk *with a thread
  running ahead of it* against a linear scan of the file.

`cargo test` on its own still fails for three examples that have not compiled since before this
work and one integration test that fails on `main` and did before — neither touched here.
