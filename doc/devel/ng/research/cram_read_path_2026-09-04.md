# Reading a CRAM: where the time goes, and what a fork of noodles buys

**Date:** 2026-09-04
**Branch:** `ng-cram-perf`
**Question asked:** the CRAM decode looked like most of a calling run's wall time. Could
modifying noodles make it faster, and could requiring sorted files and reading them
sequentially cut the memory?

**Answer in one line.** Three changes, all measured and all leaving every decoded read
byte-identical, take the CRAM read path from 0.847 s to 0.454 s per 600,000 whole-genome
reads — **1.87×** — and take a whole calling run from 79.0 s to 58.1 s — **1.36×**. A fourth,
built and verified but not yet wired into ng, cuts what the decode holds resident by
**6.4×**: 198 MB to 31 MB on a tomato chromosome, 493 MB to 77 MB on a human one. None of
these is the change the previous review pointed at, and the thing that review pointed at
turns out to be an artefact of this repository's own benchmark files.

---

## 1. The instrument, because there was none

`perf_ng-calling_2026-09-02.md` §3 recorded the gap in terms: *"No benchmark or harness in
the repository reads a CRAM."* So the first thing built was one —
[`examples/ng_cram_decode_layers.rs`](../../../../examples/ng_cram_decode_layers.rs).

It walks a fixed run of containers four times, each pass doing one layer more than the last,
so each layer's own cost is a difference between two passes:

| pass | what it adds |
|---|---|
| `read` | seek, and pull each container's bytes into memory. No decompression. |
| `blocks` | parse the compression header and inflate every block. |
| `records` | decode every record out of the inflated blocks — and this is where noodles checks the reference digest. |
| `convert` | what ng does with a record: resolve its read group, and copy its fields into ng's flat buffers. |

It also hashes every field of every decoded record — flags, position, mapping quality, name,
bases, quality scores and CIGAR — so a change to the decoder has an oracle that is not a
timing. Every number below comes with that digest unchanged at `a127ecf5083f535a002a5f461f150ede`.

**The two files it was run on, and why both.**

- `benchmarks/tomato_big_cram/DRR000741.p1.cram` — one tomato accession, whole genome,
  49 GB, coordinate-sorted, CRAM 3.0 written by samtools 1.21, 112,140 slices of about
  7,075 bases each. **This is what a real input looks like.**
- `benchmarks/tomato1/crams/SRR7279481.p1.bench.cram` — the benchmark fixture: the same
  species cut down to 80 scattered 100 kb windows, 65 slices.

Hardware: macOS 15 on Apple Silicon, 18 logical cores, native host release build (fat LTO,
one codegen unit, `-C target-cpu=apple-m1`). Every figure is the fastest of seven passes with
the machine otherwise idle.

---

## 2. The reference digest is a fixture artefact, not a caller cost

The previous review's largest single finding was that noodles re-hashes the reference span of
every CRAM slice and that `md5::compress` was 11.2 % of busy single-threaded CPU. It also
said, correctly, that the benchmark files might be inflating it. They are, and by more than
the review supposed.

Measured on the two files, hashing the same spans the same way:

| file | slices | reference bases those slices span | digest as a share of the CRAM read path |
|---|---:|---:|---:|
| whole-genome CRAM, 60 containers | 60 | 0.5 Mb | **0.1 %** |
| benchmark fixture, 7 containers | 7 | 82.6 Mb | **47.6 %** |

A CRAM slice holds a fixed number of *records*. On a file cut down to scattered windows one
slice's records run from one window to the next, and the slice's span covers the megabases of
empty reference in between — 82.6 Mb of reference spanned to carry 61,427 reads, against
0.5 Mb to carry 600,000. **So on the benchmark the digest is half the read path and on a real
file it is one part in a thousand.**

**What follows for the project.** Nothing should be optimised against `benchmarks/tomato1`'s
CRAMs without checking the same thing on a whole-genome file; on the CRAM read path those
fixtures are not a small version of the real input, they are a different workload. The
`[patch.crates-io]` fork the previous review proposed for the digest is not worth carrying.

---

## 3. Where the time actually goes

600,000 reads of the whole-genome CRAM, one thread, seconds:

| layer | stock noodles | with §4 | with §4–5 | with §4–5 and §5a |
|---|---:|---:|---:|---:|
| pull the bytes off disk | 0.001 | 0.001 | 0.001 | 0.001 |
| inflate the blocks | 0.346 | **0.266** | 0.266 | 0.266 |
| decode the records | 0.145 | 0.144 | 0.144 | **0.049** |
| build ng's read from each record | 0.355 | 0.355 | **0.139** | 0.139 |
| **total** | **0.847** | **0.767** | **0.551** | **0.454** |

Reading bytes off disk is 1 part in 800. Two things are large and they are roughly equal:
**inflating the compressed blocks (41 %)** and **turning a decoded record into ng's own read
(42 %)**. The record decode itself — unpacking the data series, resolving mates, rebuilding
bases against the reference — is 17 %.

**For scale against the reference implementation.** `samtools view -O sam` on the same 600,254
reads takes 0.42 s, and that includes formatting them as 210 MB of text and writing it. Stock
noodles takes 0.492 s to reach decoded records and does not format anything. So htslib decodes
and formats in less time than noodles decodes — the gap is real, and §4 and §5 close about
half of it.

---

## 4. noodles builds a megabyte of lookup table for every order-1 rANS block, mostly for
symbols that never occur

**What the codecs cost.** Instrumenting the vendored copy's block decode, over the same 60
containers:

| compression method | blocks | compressed | inflated | seconds |
|---|---:|---:|---:|---:|
| none | 374 | 0.0 MB | 0.0 MB | 0.000 |
| gzip | 529 | 7.0 MB | 35.5 MB | 0.030 |
| **rANS 4x8** | **1,351** | **17.3 MB** | **74.5 MB** | **0.323** |

**rANS is 91 % of decompression, and it inflates at 231 MB/s where gzip does 1,183 MB/s** —
one fifth the throughput of the general-purpose codec beside it, on the same file. gzip here
is `flate2` on the `zlib-rs` backend; rANS is noodles' own Rust.

**Where that goes.** Splitting the rANS time into building the decode tables and decoding
symbols: 1,051 of the blocks are order-0 over 10.1 MB, 300 are order-1 over 64.5 MB, and
**table building is 0.161 s of the 0.397 s** — 41 % of the codec, before a symbol is decoded.

The order-1 model keys on the previous symbol, so it declares 256 contexts and
`build_cumulative_frequencies_symbols_table` fills a 4,096-entry lookup for every one of them:
a megabyte of table per block. A block of DNA bases uses about five of those contexts and a
block of quality scores a few dozen; the rest have an all-zero frequency row that no state can
ever land in.

**The change is to skip them** — `noodles-cram/src/codecs/rans_4x8/decode/order_1.rs`, about
fifteen lines. A row is skipped only when every frequency in it is zero, and the decoder
indexes row *i* only after emitting symbol *i*, which requires a non-zero frequency for *i*
somewhere; so a decoder reaching a skipped row would already have decoded a symbol the
frequency table says cannot occur.

**Measured: block inflation 0.346 s → 0.266 s, −23 %; the whole read path −9.4 %.** All 218 of
noodles-cram's own tests pass, and the record digest is unchanged.

**What this does not cover.** These files are CRAM 3.0, so only gzip and rANS 4x8 appear.
CRAM 3.1 files — `samtools --output-fmt-option version=3.1` — use rANS Nx16, fqzcomp and the
name tokenizer instead, and reading noodles' implementations of those suggests they are worse
still: Nx16 does a linear scan of up to 256 cumulative frequencies *per decoded byte* where
4x8 uses a lookup table, and fqzcomp constructs 65,536 adaptive models, each holding two
`Vec`s, before decoding a block. Neither is measured here because no file in this repository
uses them. **If ng is ever handed 3.1 input, this section has to be re-run before anything is
concluded from it.**

---

## 5. Every record is built twice: once as noodles' `RecordBuf`, once as ng's read

`decode_container_at` called `RecordBuf::try_clone_from_alignment_record` and then copied the
result into ng's flat buffers. The `RecordBuf` is a whole second copy of the read, and getting
the fields out of the CRAM record to build it goes through
`sam::alignment::Record`, whose accessors each return a `Box<dyn …>` — **eight heap
allocations per record before a byte is copied**, and the tags are decoded into it although
ng drops them on the next line.

The change is a small direct API on the vendored `cram::Record`
(`vendor/noodles-cram/src/record/direct.rs`): `write_bases_into`, `write_quality_scores_into`,
`write_cigar_into`, `name_bytes`, `read_group_index` — each writing into a buffer the caller
owns and reuses for a whole container. `decode_container_at` then fills its flat buffers
straight from the CRAM record and `RecordBuf` disappears from the path.

**Measured: building ng's read from each record 0.355 s → 0.139 s, 2.6×; the whole read path
−28 %.** The harness hashes every field both ways in the same run and the two digests agree,
so this is not "the tests still pass" — it is the same 600,000 reads, field for field.

---

## 5a. Four fifths of decoding a record was reading tags ng throws away

`Records::read_data` decoded every auxiliary tag of every record and there was no way to ask
it not to. ng reads none of them: the one thing it wants from a record's tags is the read
group, and **a CRAM does not store that as a tag** — it stores a number, an index into the
header's `@RG` list, which `data().get(&READ_GROUP)` answers from that number without
touching a stored field.

**Skipping them is not simply a matter of not asking.** A CRAM's tags and its data series are
decoded out of the same set of streams, and a stream is a cursor: a read skipped when
something else depended on it having happened leaves every value after it wrong and says
nothing. So the fork decides the question from the compression header before any record is
read — `tag_streams::tags_can_be_left_unread` — and refuses unless every tag in the container
reads only external blocks that no data series reads. When it refuses, the tags are read as
before and only their values are discarded, which is smaller and always safe.

One case earns its own arm, and finding it is the difference between 3 % and 20 %: **a
one-symbol Huffman alphabet reads no bits**, because it encodes a constant, and that is how a
writer stores the length of a fixed-width tag. Refusing it refused every real file here.

Measured, each on records that hash the same with the tags and without:

| file | before | after |
|---|---:|---:|
| whole-genome tomato, 600,000 reads | 0.570 s | **0.454 s** |
| GIAB human, 400,000 reads | 0.551 s | 0.505 s |
| the region-subset tomato fixture, 61,427 reads | 0.205 s | 0.190 s |

**A correction to a number earlier in this note's own working.** An instrumented build put the
tags at 0.232 s of 0.287 s of record decoding. That was mostly the instrument — an `Instant`
pair per record over 1.2 million records. Merely discarding the decoded values, the safe
fallback above, is worth 3.4 %, not 40 %; the 20 % comes from not reading the streams at all,
which is a different change and needed the check to be safe.

---

## 6. End to end, on a real call

`call-from-alignments`, one sample, 10 Mb of tomato chromosome 1 out of the whole-genome CRAM,
`--defaults`, three repeats each with the arms alternated:

| | wall (18 threads) | wall (1 thread) | peak resident |
|---|---:|---:|---:|
| `main` as it stands | 81.5, 80.4 s | **79.0 s** | 298–322 MB |
| §4, §5 and §5a applied | 60.6, 59.6 s | **58.1 s** | 322 MB |

**1.36×, and the VCF is byte-identical** — 50,203 records, compared line for line except
`##commandline` and `##parametersFile`. Peak resident does not move. The single-threaded pair
is the one quoted because the 18-thread figures swing by 10 s between repeats of the same
binary on this machine while the one-thread pair does not, and a single-sample run has little
to parallelise anyway. The 63-accession
benchmark cohort over 20 of its regions is also byte-identical and shows no wall-time change,
which is what §2 predicts: on those fixtures the read path's cost is the reference digest,
which neither change touches.

The library's 6,189 tests pass, and noodles-cram's own 218.

---

## 7. Memory: the sliding-window reference is built, and it is 6.4×

The suggestion was that requiring sorted files and reading them sequentially would let ng hold
a window of the reference rather than a whole chromosome. **It works, it is built in the fork,
and it saves more than the first estimate in this note said** — because noodles was holding
not one copy of the chromosome but, at its peak, two.

**Why it was possible.** A CRAM may embed its own reference, and for those noodles already
rebases coordinates against the embedded block's start. A caller-supplied window is the same
problem — bases whose first byte is not position 1 — so the fork adds a `Window` variant
beside `Embedded` and the two indexing sites handle them together. Everything a window needs
is known before the fetch: the slice header states `alignment_start` and `alignment_span`, and
is parsed before any block is decoded. `Slice::reference_span` now says so and
`Slice::records_over_window` takes the bases.

**Measured**, decoding the same containers into records that hash the same, `/usr/bin/time -l`
on the harness with one pass running so the process holds only what that pass needs:

| | whole contig | per-slice window | |
|---|---:|---:|---|
| tomato chromosome 1 (90.9 Mb) | 198.2 MB | **31.0 MB** | 0.453 s → 0.454 s |
| human chromosome 1 (249.0 Mb) | 492.8 MB | **76.5 MB** | 0.500 s → 0.530 s |

**6.4× less resident, for the same records and the same time.** 198 MB for a 90.9 Mb
chromosome is about two bytes a base, and that is the finding the earlier estimate missed:
`Repository::get` fetches the contig *and clones it into a cache that never evicts*, so the
peak holds both. The human file's 6 % of extra time is its 2.5 Mb slice spans — its windows
are large and re-read once per slice.

**A window that does not cover the slice is refused and the message names both spans.** The
failure that prevents is the silent one: a short window decodes every read against bases that
are not under it. Checked by running the harness with the window one base short.

**What is not built: ng's side.** ng holds one `fasta::Repository` for the whole run, shared
by every open file, and switching to windows means one sliding window shared by the cohort
instead — which has to be safe for the samples decoding in parallel and has to advance
without ever going backwards. That is a design change in `read/input/reference.rs`, not a
call-site swap, and it is its own piece of work.

**And the proportion has not changed even though the number has.** The reference is shared by
every file in a run, so this is a *fixed* cost: on the one-sample whole-genome run whose peak
is 271 MB it is most of what is left after the changes above, and on the 63-accession
benchmark's 2,257 MB it is under a tenth. **It helps most where memory is least of a
problem** — one sample against a large genome — and a cohort's memory is elsewhere.

**One cohort-scaling piece is already taken.** Every open CRAM held its `.crai` twice, once
flat and once bucketed by contig, and on the CRAM arm the flat copy was never read again:
56 bytes a slice, and this repository's whole-genome tomato CRAM has 112,140 slices, so 6.3 MB
a file and 396 MB across 63. It is now `Option<AlignmentIndex>` and the grouping returns one
or the other, so both cannot be held. **Peak resident on a one-file run does not fall** — it
rises 2.8 MB, reproducibly, because the peak falls during reference loading and freeing 6.3 MB
just before it leaves the allocator holding a segment the contig does not fit into. It is a
live-heap change taken for the cohort case this machine cannot run.

**One thing sequential reading would not fix.** ng already reads containers in file order
through the `.crai`, and already keeps exactly one decoded container per open file rather than
caching. What requiring sorted input would additionally buy is the index itself — a purely
sequential walk needs none — which is the 6.3 MB a file above, now taken by a narrower change
that does not require the guarantee.

---

## 8. What was checked and is not worth doing

- **The reference digest fork** (previous review's largest finding). One part in a thousand on
  a real file — §2.
- **Buffering the CRAM reader.** noodles' `Reader<R>` does not wrap `R`, so ng's
  `Reader<File>` issues about twenty one-to-eight-byte `read(2)` calls per container header.
  It is real and it is 2 parts in 1,000 of the read path on a warm page cache; the whole
  `read` layer is 0.001 s of 0.847 s.
- **A faster MD5.** `md-5` 0.11 has no assembly backend on aarch64, and §2 makes the question
  moot.
- **Merely discarding the decoded tag values** rather than not reading them: 3.4 %, against
  20 % for not reading them — §5a.

## 9. What is next, in the order it pays

1. **Wire ng to the windowed reference** — §7. The noodles half is built and verified; the ng
   half is one sliding window shared by the cohort in place of a `Repository`, which is a
   design change in `read/input/reference.rs`. Worth 167 MB on a tomato chromosome and 416 MB
   on a human one, fixed per run.
2. **rANS symbol decoding is now the largest single item left** — 0.196 s of the 0.454 s the
   read path costs, 43 %. A sampling profile of the harness puts `Block::decode_inner` at 8,815
   self-samples of 25 seconds, four times the next leaf. What was *not* tried is making the
   inner loop cheaper: it does three table lookups a byte (symbol, frequency, cumulative
   frequency) where htslib packs all three into one word, and it renormalises a byte at a time
   through an `io::Result`. Both are plausible and neither is measured, so neither is a
   finding yet.
3. **Re-measure everything on a CRAM 3.1 file before ever claiming it applies there** — §4.

**What was measured and left alone.** After the tag skip, `RandomState::hash_one` was the
third-largest leaf in the profile at 1,013 samples — noodles keys tag content ids and its
external data readers in `HashMap`s with the default SipHash. Most of that is gone for a
caller that skips tags, which is why it is recorded rather than pursued.

**The fork is the open question and it is the owner's.** §4 is a change to noodles that should
go upstream: it is a clear defect, it is fifteen lines, and it needs no new API. §5, §5a and
§7's decode are not — each adds an interface noodles has no other caller for, so they either
live in a fork this project carries and rebases at every noodles bump, or they are proposed
upstream and wait. All of them live in `vendor/noodles-cram`, pinned by a `[patch.crates-io]`
entry against the same 0.93.0 the manifest already names, and `FORK.md` beside them lists
every one with what it costs and whether it belongs upstream.

---

## 10. Where each piece is

The branch is `ng-cram-perf`, and each change is its own commit so that the two that touch
noodles can be judged — and if the owner decides so, offered upstream — one at a time.

| commit | what it is | upstream? |
|---|---|---|
| `build(ng): vendor noodles-cram 0.93.0 unchanged` | the copy, byte for byte the registry's, and `[patch.crates-io]` pointing the build at it. `vendor/noodles-cram/FORK.md` is the list of what this copy has that crates.io does not, checkable with one `diff -r`. | — |
| `test(ng): a harness that times the CRAM read path layer by layer` | `examples/ng_cram_decode_layers.rs`, and the per-codec counters behind the `cram-perf-counters` feature. | no — an instrument |
| `perf(cram): rANS builds a megabyte of lookup table per block…` | §4. | **yes** — a defect, no API change |
| `perf(cram): read a decoded record's fields into the caller's buffers…` | §5's noodles half. | not as it stands — a new API with one caller |
| `perf(ng): a CRAM record goes straight into the container's buffers` | §5's ng half. | — |
| `perf(cram): a caller that reads no auxiliary tags need not decode them…` | §5a's noodles half, with the check that says when it is safe. | not as it stands |
| `perf(ng): the CRAM decode stops reading tags it was already throwing away` | §5a's ng half, and the 1.36× of §6. | — |
| `perf(ng): a CRAM's `.crai` is grouped by contig at open and the flat copy is then let go` | §7's cohort-scaling piece. | — |
| `perf(cram): decode a slice against a window of the reference…` | §7's windowed decode. **ng does not use it yet.** | the shape may be worth proposing |
| this report | | |
