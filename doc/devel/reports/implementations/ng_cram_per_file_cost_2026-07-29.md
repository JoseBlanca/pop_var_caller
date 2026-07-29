# ng cohort CRAM ingest — what an open file costs, and why it was three times too much

**Date:** 2026-07-29 · **Branch:** `ng-cram-repository-sharing` ·
**Follows:** [ng_cram_reference_sharing_2026-07-29.md](ng_cram_reference_sharing_2026-07-29.md),
which removed the reference from the per-file cost and left this behind.

Once the reference was shared and bounded to one contig, a cohort's peak memory was ~166 MiB fixed
plus **14.6 MiB for every file opened** — at 51 samples, 82% of the total, and ~30 GiB at the
2,085-file archive this is heading for. This is what that 14.6 MiB was and what it is now.

## The answer

Each open file retained one fully decoded CRAM container: **12.7 MiB of live heap in ~72,000
allocations**. It is `DecodedContainer` (`region_query.rs`), the cache from `b918fb6` that took the
per-locus cost from ~38 ms to ~0.6 ms. It travels with the pooled reader, so there is one per open
file, held for the life of the run.

Attribution, from diffing DHAT profiles at 16 and 32 files (per file):

| MiB/file | site | what it is |
| -------: | ---- | ---------- |
| 6.10 | `insert` — noodles-sam `record_buf/data.rs:228` | the auxiliary-tag map on each record |
| 2.75 | `decode_container_at` — `region_query.rs` | the `Vec<RecordBuf>` itself |
| 2.83 | `try_clone_from_alignment_record` ×3 | sequence, qualities, CIGAR copied out |
| 0.71 | `bstr::from` | read names |
| 0.31 | `collect<… -> Footprint>` | the footprint table |
| **12.71** | | ~1.3 KB live per record, ~10⁴ records a container |

## Method, including the part that went wrong

**DHAT's `At t-gmax` nearly produced the opposite conclusion.** It is live heap at *one instant*,
the global peak. Below ~16 files that instant falls during reference loading, so `t-gmax` came out
identical at 1 and 8 files (219.33 vs 219.50 MiB) and appeared to prove that nothing is retained
per file. It proved no such thing: the per-file retention had not yet grown past the reference
spike. **Comparing a peak across configurations is only valid when the peak falls in the same
phase.** An explicit `HeapStats::get()` probe at a known point — files open, walk finished — is
what settled it.

Two false trails were closed on the way, both worth the time they cost:

- **Not the allocator.** All growth is anonymous (`RssFile` flat at 2.3 MiB), but glibc's
  dynamic-mmap-threshold fix moved it 5% and mimalloc moved the slope 12.0 → 10.8 MiB/file. Freed
  pages were not the story.
- **Not merely holding files open.** Opening 32 files costs 3.8 MiB *total* (38 KiB each).
  Querying **one** span through each jumps to 507 MiB, and ten spans or a hundred add nothing. The
  whole cost arrives on the first query and stays — which is what named the container cache.

The harness is [dhat_ng_open_files.rs](../../../../examples/dhat_ng_open_files.rs); it runs the
input layer alone, so what scales with file count there scales with it in the cohort dump and
nothing else is in the way. `tmp/dhat_diff.py` diffs two profiles by allocation site.

## What was fixed

### 1. The tags are dropped, and the read group resolved at decode

Nothing in ng reads an auxiliary tag except `resolve_read_group`, and only the `RG` one:
`AlignedRead`, the only thing a record ever becomes, has no tag field, so no other consumer can be
reading one. The whole map was decoded, stored, cloned per query and never read.

`DecodedContainer` now resolves each record's read group **while the tags are still in hand**, and
then drops every tag. Resolving at decode is what generalises it: the tags go on both the
single-read-group and the several-read-groups arms, not just the first. It also lifts a per-record
tag lookup out of the query path, which ran on every record of every query for an answer that
cannot change once the record is decoded.

**`clear()` is not enough, and this is the trap.** `Data` is a `Vec<(Tag, Value)>`; `clear` drops
the values but keeps the vector's capacity, which is the larger half. Measured: `clear` saved
0.13 MiB per file where replacing the whole `Data` saved 6.2 MiB.

**One behaviour change, visible only on a malformed file:** a record with an unreadable `RG` now
fails when its container is decoded rather than when that record is served — earlier, and possibly
for records outside the queried region. The condition was already fatal either way.

### 2. The container is stored packed, not as records

A `RecordBuf` per record was ~680 bytes across seven allocations, against ~270 bytes of read that
anything consumes. The rest was `Vec` headers, capacity slack and per-allocation overhead.

The bytes now go into two flat buffers — one of names, sequences and qualities, one of CIGAR
operations — and each record becomes a fixed-size index entry naming its slices of them. A query
collects **indices**; a record is rebuilt only when it is actually served, straight into the
caller's reused buffer.

That deletes the per-query `record.clone()` as a side effect: **the CRAM query path now allocates
nothing per read**, the buffer-reuse property the BAM path already had and CRAM had given up by
moving whole records in.

`shrink_to_fit` on the three buffers after decode is worth its own line — 1.3 MiB per file was pure
doubling slack, held for the life of the run, for one copy per container against thousands of
reads.

**A hazard found while writing it:** the index's `u32` offsets would silently wrap on a container
holding more than 4 GiB of read data and serve another record's bytes — wrong reads, not a crash.
Nothing near that exists in practice (10⁴ records a container), but the CRAM spec's record count is
a 32-bit field, so it now refuses rather than wraps.

## Results

Live heap per open file, measured with the harness at 1, 8 and 32 files:

| | MiB/file |
| --- | -------: |
| before | 12.71 |
| tags dropped | 6.48 |
| packed container | 5.35 |
| + `shrink_to_fit` | **4.04** |

Peak RSS of the cohort dump, 556-span BED — the last column is this work, the rest is
[the reference report](ng_cram_reference_sharing_2026-07-29.md):

| files | original | one repository | + one contig | + this work |
| ----: | -------: | -------------: | -----------: | ----------: |
|     1 |  835.1 MiB |    840.1 MiB |    180.4 MiB |   180.4 MiB |
|     8 | 6105.8 MiB |    946.9 MiB |    280.8 MiB |   192.8 MiB |
|    16 |12125.8 MiB |   1063.7 MiB |    395.3 MiB |   235.9 MiB |
|    32 |    *OOM*   |   1307.0 MiB |    638.7 MiB |   322.5 MiB |
|    51 |    *OOM*   |   1590.6 MiB |    908.9 MiB |   422.1 MiB |

The 51-sample cohort that started this — OOM-killed at ~80 s against a 16 GB cap, extrapolating to
37.6 GiB — now peaks at **422 MiB** and emits its 305,846 rows in 154 s.

**Wall time improves too**: 49.4 s → 42.4 s at 8 files, 83.9 s → 65.7 s at 16, and 51 samples run
in 154 s against 209.6 s for the shared-reference build. Decode pays one extra memcpy per record
(89 containers a file — nothing); every query loses N allocations.

**Verified:** output byte-identical to the pre-change binary at 8 samples (47,462 rows) after each
of the two changes and again on the final code; `cargo test` 2656 lib tests plus every integration
target, 0 failures; `cargo clippy --all-targets` clean.

## What is left

The per-file term is 4.04 MiB, of which 2.60 MiB is the reads themselves — name, sequence,
qualities, CIGAR — and ~1.4 MiB is the index and CIGAR buffers. Going below that means storing less
than the reads, which is a different question from this one.

Two things were considered and rejected, with reasons worth keeping:

- **Cache the decompressed blocks instead of the records** (1.93 MiB a container, *smaller* than
  the payload it encodes, because CRAM stores sequences as differences from the reference). Blocked
  by CRAM's structure: `cram::Record` borrows the decompressed buffers, so caching both is
  self-referential, and the data series are sequentially coded — record 4,000 cannot be decoded
  without the 3,999 before it. noodles exposes `slice.records()` (all of them) and nothing finer,
  so the only way to use a block cache is to re-materialise every record per query, which hands
  back `b918fb6` in full.
- **Decode one slice of a container at a time.** These files have exactly one slice per container
  (89 of 89, checked in the `.crai`), which is samtools' default. There is nothing to skip.
