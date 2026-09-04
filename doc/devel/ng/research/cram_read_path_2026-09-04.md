# Reading a CRAM: where the time goes, and what a fork of noodles buys

**Date:** 2026-09-04
**Branch:** `ng-cram-perf`
**Question asked:** the CRAM decode looked like most of a calling run's wall time. Could
modifying noodles make it faster, and could requiring sorted files and reading them
sequentially cut the memory?

**Answer in one line.** Two changes, both measured and both leaving every decoded read
byte-identical, take the CRAM read path from 0.847 s to 0.551 s per 600,000 whole-genome
reads — **1.54×** — and take a whole calling run from 79.4 s to 66.2 s — **1.20×**. Neither
is the change the previous review pointed at, and the thing that review pointed at turns out
to be an artefact of this repository's own benchmark files.

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

| layer | stock noodles | with §4 | with §4 and §5 |
|---|---:|---:|---:|
| pull the bytes off disk | 0.001 | 0.001 | 0.001 |
| inflate the blocks | 0.346 | **0.266** | 0.266 |
| decode the records | 0.145 | 0.144 | 0.144 |
| build ng's read from each record | 0.355 | 0.355 | **0.139** |
| **total** | **0.847** | **0.767** | **0.551** |

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

## 6. End to end, on a real call

`call-from-alignments`, one sample, 10 Mb of tomato chromosome 1 out of the whole-genome CRAM,
`--defaults`, three repeats each with the arms alternated:

| | wall (18 threads) | wall (1 thread) | peak resident |
|---|---:|---:|---:|
| `main` as it stands | 79.4, 79.2, 77.9, 80.2, 79.1 s | 78.1 s | 298–322 MB |
| §4 and §5 applied | 66.2, 65.7, 66.3, 69.2, 67.0 s | 64.9 s | 322–352 MB |

**1.20×, and the VCF is byte-identical** — 50,203 records, compared line for line except
`##commandline` and `##parametersFile`. Peak resident does not move. The 63-accession
benchmark cohort over 20 of its regions is also byte-identical and shows no wall-time change,
which is what §2 predicts: on those fixtures the read path's cost is the reference digest,
which neither change touches.

The library's 6,189 tests pass, and noodles-cram's own 218.

---

## 7. Memory: what a sliding-window reference would buy, and what it would not

The suggestion was that requiring sorted files and reading them sequentially would let ng hold
a window of the reference rather than a whole chromosome, as `WindowedRefSeq` already does for
the walk's own view of the bases.

**It would work, and the mechanism is already in noodles for a different case.** A CRAM
decodes a read's bases against the reference by absolute coordinate: the record decoder is
handed `(**sequence).as_ref()` — the whole contig — plus an absolute alignment start, and
indexes it at `position - 1`. But noodles already carries a *rebased* variant beside it, for
CRAMs that embed their own reference: `ReferenceSequence::Embedded` holds a `reference_start`
and the record decoder subtracts it. Giving the external variant the same field, and passing
the window's first position through the digest check, is the whole of the change. Everything a
window needs is known before the fetch: the slice header carries `reference_sequence_id`,
`alignment_start` and `alignment_span`, parsed twenty lines before the repository is asked.

**What it is worth, measured.** Calling 200 kb out of the whole-genome CRAM, so that almost
nothing but the fixed cost is resident:

| contig held | length | peak resident |
|---|---:|---:|
| `SL4.0ch00` | 9.6 Mb | 150 MB |
| `SL4.0ch01` | 90.9 Mb | 271 MB |

121 MB for 81.2 Mb of difference — **about 1.5 bytes of resident memory per base of the
contig in hand**, which is one copy in noodles' repository plus ng's own window. So on tomato
a windowed reference is worth roughly 130 MB, and on a human chromosome 1 (248.9 Mb) roughly
370 MB.

**And here is the part that decides whether to build it: that cost is fixed, not per sample.**
One reference is shared by every file in a run (`read/input/reference.rs`), so it is 48 % of a
one-sample run's 271 MB peak and 4 % of the 63-sample benchmark's 2,316 MB. **A windowed
reference helps most exactly where memory is least of a problem** — a single sample against a
large genome — and barely moves a cohort. It is the right change for one human sample on a
small machine and the wrong place to look for a thousand-sample run.

**What does scale with the cohort, and is nearly free to fix.** Every open CRAM holds its
`.crai` twice: once as the flat `AlignmentIndex::Crai` and once grouped by contig
(`read/input/open_bam.rs`). On the CRAM arm the flat copy is **never read again after open** —
`self.index` is used only when building a BAM cursor. A `crai::Record` is 56 bytes and this
file's index has 112,140 of them, so that is 6.3 MB per open file held for nothing: 400 MB
across 63 whole-genome samples. Dropping it is a field's type, not an algorithm.

**One thing sequential reading would not fix.** ng already reads containers in file order
through the `.crai`, and already keeps exactly one decoded container per open file rather than
caching. The gain from requiring sorted input is the index itself — a purely sequential walk
does not need one — not the walk.

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

## 9. What is next, in the order it pays

1. **Land §4 and §5.** Both are measured, both leave the VCF byte-identical, and together they
   are 1.20× on a whole-genome single-sample call. They need a decision about how the fork of
   noodles-cram is carried — see below.
2. **Drop the unread flat `.crai`** — §7, a few lines, ~400 MB on a 63-sample whole-genome
   cohort.
3. **The windowed reference** — §7, if and only if the target is one sample against a large
   genome.
4. **Re-measure everything on a CRAM 3.1 file before ever claiming it applies there** — §4.

**The fork is the open question and it is the owner's.** §4 is a change to noodles that should
go upstream: it is a clear defect, it is fifteen lines, and it needs no new API. §5 is not —
it adds an API noodles has no other caller for, so it either lives in a fork this project
carries and rebases at every noodles bump, or it is proposed upstream as a
`Record::write_*_into` family and waits. Today both live in `vendor/noodles-cram`, pinned by a
`[patch.crates-io]` entry against the same 0.93.0 the manifest already names.

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
| `perf(ng): a CRAM record goes straight into the container's buffers` | §5's ng half, and the end-to-end numbers of §6. | — |
| this report | | |
