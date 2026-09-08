# `generate-psps`: each CRAM container read once, and the psp compressed off the walking thread

**Date:** 2026-09-08
**Branch:** `main`, commits `a2a77780`, `aa598441` and `c9018a70` on `8203d218`
**Acting on:** the BAM/CRAM → psp performance review of 2026-09-08

**One tomato accession over 10 Mb of `SL4.0ch01`, out of a 49 GB whole-genome CRAM, now takes
43.99 s and 249.0 MB where it took 57.03 s and 264.6 MB** — 1.30× the throughput for 16 MB less
resident memory, with the psp and the census byte-identical outside their own timestamp. Two
changes carry it: the walk's two readers were decoding every container of the file separately,
and the walk was compressing every psp block itself.

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
pairs at 10 Mb: **57.03 s → 47.20 s** and **264.6 MB → 245.7 MB** for these two commits, and
**43.99 s and 249.0 MB** with §6a's compression offload on top. The pinned figures above are what
the arms are *compared* on; these are what the command costs.

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
rather than sitting with the other small wins. What is left on that thread is encoding records
into block payloads — order-dependent *within* a block and independent *between* blocks — so if
3.3× at four workers is not enough, building a whole block off-thread is the next place to look,
and it is an encoding question rather than a scheduling one.

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

## 7. What was measured and stopped

**The second candidate under §11 question 8 — a prefetcher moving the decode and the read
preparation off the walking thread — is not recommended, and the reason is that its ceiling moved
under it.** The review priced it from a profile taken before any of this: 20% decode plus 9% read
preparation, so 29% and about 1.4×. Re-measured the same way after each change:

| | CRAM decode | reading and filtering the reads, and the reference bases the decode needs | together | ceiling |
|---|---:|---:|---:|---:|
| before | 20% | 9% | **29%** | 1.41× |
| with the containers shared | 12.0% | 7.7% | **19.7%** | 1.25× |
| and with the four smaller changes | 11.0% | 6.1% | **17.1%** | 1.21× |
| and with the psp compression on its own thread | 11.9% | 6.6% | **18.5%** | 1.23× |

macOS `sample`, 25–30 s into a 10 Mb walk, `RAYON_NUM_THREADS=1`, parked and waiting threads
excluded so the denominator is the walking thread alone. **The last row goes back up because the
walking thread got smaller, not because the decode got bigger**: the decode's own samples are flat
at about 2,460 while the thread it sits on fell from 22,480 to 20,643.

**1.23× is the ceiling and not the gain**: a bounded queue never hides a stage completely, and the
design adds a copy of every prepared read through a channel — the same review's compression
offload, the one stage-offload anyone here has actually built, bought its wall at the cost of 1.1%
more user time for exactly that reason.

What it costs is two prefetch threads, two bounded queues and a prefetcher that has to be
stoppable and rewindable — two, because the two locus generators hold separate cursors whose
regions interleave. And it wins nothing on a cohort run as one invocation a sample, where the
machine is already saturated, so it needs a knob as well.

**Question 3's split of the walk across segments is the axis worth building instead.** Emulated as
concurrent processes over disjoint BEDs — an upper bound, since processes share nothing — it
reached 2.97× at four ways and 3.88× at eight on this same 10 Mb fixture. It needs the same knob,
and a sharded worker holds its own cursors and decodes on its own thread anyway, so **the 1.23× is
largely inside the 2.97×**. §11 question 8 now says so.

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
- **`cargo test --lib`: 6,658 passing**, six of them new. noodles-cram's own 223 and 75 pass.
  `diff -r` against the registry copy prints exactly FORK.md's list.
- **Mutation-verified.** A `take` that ignores the offset fails
  `two_cram_cursors_over_one_file_each_see_every_read`; one cache slot instead of two fails the
  policy test; the unclamped slice start fails the mismatch-loop comparison; and a drop that
  joins the compressing thread without first closing its channel hangs
  `writers_abandoned_with_blocks_in_flight_neither_hang_nor_pile_up`.
- **The CRAM decode was not re-run for §6a**, which touches only the psp writer.

`cargo test` on its own still fails for three examples that have not compiled since before this
work and one integration test that fails on `main` and did before — neither touched here.
