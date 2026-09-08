# `generate-psps` reads each CRAM container once instead of twice

**Date:** 2026-09-08
**Branch:** `main`, commits `a2a77780` and `aa598441` on `8203d218`
**Acting on:** the BAM/CRAM → psp performance review of 2026-09-08

**One tomato accession over 10 Mb of `SL4.0ch01`, out of a 49 GB whole-genome CRAM, now takes
47.20 s and 245.7 MB where it took 57.03 s and 264.6 MB** — 1.21× the throughput for 19 MB less
resident memory, with the psp and the census byte-identical outside their own timestamp. Most of
it is one change: the walk's two readers were decoding every container of the file separately.

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
pairs at 10 Mb: **57.03 s → 47.20 s** and **264.6 MB → 245.7 MB**. The pinned figures above are
what the arms are *compared* on; these are what the command costs.

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

macOS `sample`, 30 s into a 10 Mb walk, `RAYON_NUM_THREADS=1`, parked threads excluded from the
busy count.

**1.21× is the ceiling and not the gain**: a bounded queue never hides a stage completely, and the
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
and a sharded worker holds its own cursors and decodes on its own thread anyway, so **the 1.21× is
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
- **`cargo test --lib`: 6,657 passing**, five of them new. noodles-cram's own 223 and 75 pass.
  `diff -r` against the registry copy prints exactly FORK.md's list.
- **Mutation-verified.** A `take` that ignores the offset fails
  `two_cram_cursors_over_one_file_each_see_every_read`; one cache slot instead of two fails the
  policy test; the unclamped slice start fails the mismatch-loop comparison.

`cargo test` on its own still fails for three examples that have not compiled since before this
work and one integration test that fails on `main` and did before — neither touched here.
