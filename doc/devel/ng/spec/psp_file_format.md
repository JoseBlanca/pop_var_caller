# ng — the psp file format: the container

*Design spec, 2026-08-25. **No code for the format itself**; a working writer and reader for
the compression scheme exist as a measuring program (`examples/psp_row_stream_roundtrip.rs`) and
every number here was taken from it or from production.*

*This document owns the **container**: the sections a file is made of, how a reader finds them,
what a reader is allowed to hold, and what the writer and reader offer their callers. What goes
**inside** a record — which fields, how each is encoded, what may be approximated — is
[`psp_record_encoding.md`](psp_record_encoding.md), and the read names are
[`psp_chain_id_encoding.md`](psp_chain_id_encoding.md). The requirement both serve is a
per-open-sample memory budget, stated in §1.1.*

*Companion measurements: [the memory review](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md).*

---

## 1. What this is, in one page

A **psp** holds everything one sample's reads showed at every reference position a run analysed —
one record per covered position, at three reads a position and at three hundred, for a cohort of
one sample and of several thousand. A caller opens one per sample and holds them all open for the
whole run.

Anything a psp costs to have open is multiplied by the cohort size, so at a thousand samples the
number that decides whether a run fits on a machine is what one open file costs. Measured on production's format: **2.6 MB per open sample**, which is 2.6 GB at a thousand
and 7.7 GB at three thousand.

**The cause is that production ties two things to one number.** A `.psp` groups records into
blocks and compresses each block whole, so a block is simultaneously *how far back the compressor
may look for a repeated pattern* and *how much a reader must inflate before it can hand out its
first record*. Making blocks smaller to save memory therefore costs compression, and multiplies the
block index besides. Measured across that trade on a 50-accession cohort: 80 kb blocks cost 2,511 MB
on disk and 1,666 MB of memory; 1 kb blocks cost 3,930 MB on disk and 192 MB of memory. **Every
setting is a payment.**

**This format unties them.** The block stays large — for its compression ratio and for a small
index — while a separate, declared number caps how much a reader must hold, and the reader decodes
a block incrementally instead of inflating it. Measured against production's reader on the same
records: **0.34 MB per open sample instead of 2.6 MB, with the file 35 % smaller, the index 11×
smaller and the read 1.8× faster.** Nothing is traded away, which is what makes this different from
any point on the curve above.

### 1.1 Goals

1. **An open file costs no more than 500 kB of resident memory**, and does not grow with the block
   size, the depth or the length of the genome. *(The owner's working budget, 2026-08-25: 1.5 GB
   across three thousand samples. It supersedes the "tens of kilobytes" of
   [`run_streaming.md`](run_streaming.md) §7.2, which should be corrected when that document is next
   touched.)* **The independence from depth is measured, not assumed — §5.2.**
2. **A reader can start at any block** without reading what comes before it.
3. **A file that is not complete is refused, not read short.** A run killed part-way must not look
   like a sample with less of the genome.
4. **The file carries its own parameters.** Block size, compression settings and field encodings are
   the writer's choices, recorded in the file; a reader is driven by what it finds and assumes
   nothing.
5. **The same sample gathered at any worker count produces the same bytes.** This is what makes a
   byte-identity check on the writer mean anything.

### 1.2 Non-goals, and what this does not do

- **It does not fix which fields a record has, or how each is encoded.** That is
  [`psp_record_encoding.md`](psp_record_encoding.md). This document fixes only that the file
  *declares* those choices and where the declaration lives.
- **It does not choose a block size, a compression level, or a quantisation step.** Those are the
  user's, per the owner's ruling of 2026-08-25. §4.4 says how they travel.
- **It does not change production's `.psp`.** ng replaces it; production stays as it is. Several
  findings here apply to it unchanged and that is not a reason to touch it.
- **It does not define the cohort reader's scheduling** — which blocks it fetches, in what order,
  with how much look-ahead. That is [`run_streaming.md`](run_streaming.md)'s.
- **It does not specify a random-access API within a block.** A reader seeks to a block and streams
  from its start; reaching a record inside a block means decoding from that block's beginning. **This
  is a deliberately cheap corner** — the owner's ruling of 2026-08-25 is that the common case is
  reading a file from beginning to end, so the format is not shaped around fast arbitrary seeks.

---

## 2. Vocabulary

Production's code and this document use the word *trailer* for two different things, and a coder
moving between them will get it wrong.

| this document | what it is | production's name for it |
|---|---|---|
| **header** | plain text, at the start; what is known before any record is written | header |
| **block** | a run of consecutive records, compressed as one unit | block |
| **block index** | one entry per block, so a reader can seek | block index |
| **trailer** | a payload the writer supplies when it closes the file; may be empty | **metadata section** |
| **footer** | fixed-size, at the very end; offsets and magic. Its presence means the file is complete | **trailer** — 32 bytes there ([`src/psp/trailer.rs:28`](../../../../src/psp/trailer.rs)); this one is wider, §3.3 |

Two more terms this document needs:

- **look-back window** — when a compressor finds a repeated byte sequence it stores a
  back-reference: *copy 40 bytes from 3,000 bytes ago*. The look-back window caps how far back
  "ago" may be. **It is the reason a reader needs memory at all**: to resolve a back-reference the
  decoder must still be holding those bytes, so it must retain the last *window* bytes of what it
  has already decompressed, and nothing more. zstd takes it as a power of two, so a file declares
  the exponent.
- **genomic block size** — how much reference one block covers, and the rule that cuts one.
  **Default 100 kb** (§4.1). *Naming: there are now two "block sizes" — this one in base pairs and a
  ceiling in bytes — so the field is `genomic_block_size_bp` and the ceiling is
  `block_byte_ceiling`. Neither is called "block size" unqualified, and the unit is in the name
  because production's own `block_window_bp` sets that precedent and because the two will otherwise
  be read for each other.*
- **stream** — a block may be cut into more than one separately-compressed piece, so that a reader
  wanting only one cheap number per record does not have to decompress whole records. Each piece a
  reader has open costs it a decoder. §5.2 measures that.

---

## 3. The file

```
+---------------------------+
| header                    |  plain text: magic, length, TOML body, sentinel
+---------------------------+
| block 0                   |  one or more compressed streams
| block 1                   |
| ...                       |
+---------------------------+
| block index               |  one entry per block
+---------------------------+
| trailer                   |  the writer's closing payload; may be empty
+---------------------------+
| footer                    |  fixed size: offsets, counts, checksum, magic
+---------------------------+
```

**The order of the last three is a decision, not an accident** (the owner asked, 2026-08-25): the
index sits *before* the trailer so that rewriting the trailer means truncating at the trailer's
offset and writing forward from there, leaving the index untouched. With the trailer first, every
trailer rewrite would move the index and force it to be rewritten too. It also makes appending
cheap: the footer's index offset is where the blocks end, so a writer reopening a file truncates
there and carries on.

### 3.1 The header — plain text, and it stays that way

Production's shape is kept as it is ([`src/psp/header.rs:1-38`](../../../../src/psp/header.rs)):
a 4-byte magic, an 8-byte little-endian length for the TOML body, the body, and a sentinel line.
The length is authoritative and the sentinel is a cross-check.

**It stays plain text for one reason worth stating: `head` and a TOML parser tell you what a file
is without a special tool.** That has repeatedly been worth more than the bytes it costs, and the
header is written once per file rather than once per record.

The header carries what is known *before* any record is written:

- **the format version**, as a plain string. **This one field must be readable without knowing the
  version**, which is why the header cannot become binary: it would create a chicken-and-egg
  problem the first time the binary layout changed.
- the sample name, the reference it was called against, the contig list with lengths and checksums;
- the writer's provenance and the parameters it ran with;
- **the manifest** (§4.4) — the block cut rule, the look-back window, the stream layout, and each
  field's encoding with its parameters.

### 3.2 Blocks

A block is a run of consecutive records from one chromosome, compressed as one unit, with the
look-back window capped at the value the header declares.

**Every block is self-contained.** It opens with its chromosome, its first position and its record
count, and **every running difference inside it restarts** — the position difference, the coverage
difference, the read-name difference. A block never crosses a chromosome. So the restart points are
the block boundaries and there is no separate seek mechanism to build.

**This is the property goal 2 rests on, and the one most easily broken by accident.** A block that
fails to reset a running difference reads back wrong from its first record — and *plausibly* wrong,
because coverage is smooth and a position difference that is slightly off still parses.

### 3.3 The block index and the footer

**One entry per block**, in genomic order: the chromosome, the first position, the byte offset, and
a summary of the one number a cohort scan reads (the largest non-reference support anywhere in the
block), so a reader can decide whether to touch a block **without decompressing it**.

**A flat vector, decoded whole at open.** This is production's shape
([`BlockIndexEntry`, `src/psp/index.rs:42`](../../../../src/psp/index.rs)) and it was going to have
to change: at production's 5 kb blocks a whole-genome sample has about 156,000 entries and a 3.8 MB
index per open file, which [`run_streaming.md`](run_streaming.md) §7.2 rejects. **A large block
removes the problem rather than solving it.** Measured on a tomato accession: 154 blocks at 1 MB
blocks against 1,674 in the `.psp`. Scaled to a whole genome — arithmetic from that measured count —
roughly 14,000 entries and a few hundred kilobytes. **So the coarse-index-and-chain scheme that
document asks for must not be built.**

The **footer** is fixed-size, at the very end, and carries the offsets and lengths of the index and
the trailer, the block count, a checksum over the index, and a magic placed last so a four-byte read
at end-of-file rejects a truncated or foreign file before anything else is touched. Production's
32-byte layout ([`src/psp/trailer.rs:26-45`](../../../../src/psp/trailer.rs)) is the model; it needs
two more offsets for the trailer, so it will be wider.

**A file with no valid footer is refused** — goal 3. It is the only signal that distinguishes a
completed file from a run that was killed, and there is no safe way to read one short: a caller
would silently get a sample that stops in the middle of a chromosome.

### 3.4 The trailer

**An arbitrary payload the writer supplies when it closes the file.** It may be empty.

**Why this exists as a separate section rather than more header:** the header is what you know
*before* writing; the trailer is what you only know *after*. The per-sample summary is the case that
forces it — a coverage-against-GC histogram and heterozygosity counts accumulated *from the records
as they are written* — so it cannot be in the header, and it is read by the genotype prior, by
production's hidden-paralog filter, and by ng's duplication filter when that is built.

**The trailer's payload is binary, not text**, per the owner's ruling of 2026-08-25 and for a
measured reason: production stores that histogram as TOML, which costs **1.05 MB of text per open
sample to yield a 0.52 MB parsed histogram**, and its reader holds both for the length of the run.
[`psp_record_encoding.md`](psp_record_encoding.md) §4.1 owns that decision.

**The container does not interpret the trailer.** It stores bytes and hands them back. What is in
them is the writer's business, which is what lets ng's summary change shape without a container
version bump.

---

## 4. What the format fixes, and what the file carries

**The organising rule, from the owner (2026-08-25): the format fixes the grammar; the file carries
the values.** Block size and quantisation are the user's choices, so they cannot be constants in
this document — but a reader must still be able to read any file, so every such choice is declared.

### 4.1 The genomic block size — the rule that cuts a block

**The genomic block size is how much reference a block covers**, and it is the rule that cuts one —
not a count of bytes. **Default 100 kb** (the owner, 2026-08-25). Production has the same idea under
another name (`block_window_bp`, default 5,000,
[`src/psp/writer.rs:92,119`](../../../../src/psp/writer.rs)).

**It is a grid on the coordinate, not a running total.** A block ends when a position crosses into
the next multiple of the genomic block size. That distinction is the whole point: a grid makes every
sample cut at the *same* coordinates, and a running count does not.

Keeping the rule genomic rather than byte-based buys three things:

- **Blocks align across samples.** Every sample cuts at the same coordinate, so a cohort reader
  stepping across a region touches one aligned block per sample instead of one in some samples and
  two in others. A byte-based cut loses this, and losing it was going to be an accident rather than
  a decision.
- **Goal 2 comes free.** A span cut *is* a bound on restart granularity. A byte cut needs an extra
  rule — *also cut every so many kilobases* — to stop a block spanning an enormous stretch on a
  sparse sample. That rule is not needed here, and one fewer rule is one fewer thing to get wrong.
- **Goal 5 comes free.** A cut that depends only on the coordinate cannot depend on how the writer
  was scheduled.

**One consequence to design for rather than discover: a span-cut block has a variable size in
bytes.** At three hundred reads a position a 5 kb span is a great deal of data; on a sample with
almost no coverage it is a handful of records, and a very small block compresses badly because the
compressor starts cold each time. **So the writer may declare a secondary rule** — close a block
early if it exceeds a byte ceiling, and keep accumulating across empty spans rather than emitting
near-empty blocks. Both are the user's choices and both are recorded.

#### What the size actually costs — measured 2026-08-25

Bytes a record and blocks per sample, byte ceiling off unless stated:

| cut rule | tomato 3× | HG002 279× |
|---|---|---|
| byte count only, 1 MiB | 148 blocks · 4.629 | 3 blocks · 16.255 |
| genomic 5 kb | 1,674 · 4.579 | 629 · 18.242 |
| genomic 20 kb | 480 · 4.629 | 540 · 18.084 |
| **genomic 100 kb — the default** | **160 · 4.627** | **281 · 17.557** |
| genomic 1,000 kb | 90 · 4.626 | 34 · 16.444 |

**On contiguously covered data the size barely matters — about 1 % of the file across a two-hundred
fold range.** That is the capped window doing what §4.2 says it does: a match cannot reach past
32 kB, so a larger block gives the match finder nothing extra, and the entropy tables are already
amortised at the small end.

**On patchy coverage it costs 10 %, and that is what the secondary rule is for.** The human sample
here is 74,623 covered positions scattered over 644 small regions, so a 5 kb grid gives blocks of
about 119 records and the compressor's cold start dominates. **Its patchiness is an artefact of the
region-restricted pileup that produced it, not a property of a 279× sample** — but patchy coverage is
real for exomes, panels and thin samples, so a writer that accumulates across empty spans rather than
emitting near-empty blocks earns its place.

**Why 100 kb and not something else.** It is a round number that satisfies goal 2 directly, sits in
the flat part of the tomato curve, and recovers most of the patchy-data penalty (17.557 against
18.242 at 5 kb). **It is not an optimum** — 1,000 kb is smaller on both samples — and it is a
starting value, not a derived one — and **the argument that was holding it down has been
withdrawn.** I had reasoned that larger blocks cost seek time, since a reader starting mid-block
decodes from that block's beginning. **The owner's ruling of 2026-08-25: seeking is not the common
case — a run reads a file from the beginning to the end.** So seek cost should not set this number.

**What still argues against a very large block is skipping, which is not the same thing.** The index
carries a per-block summary so a scan can decide whether to touch a block at all (§3.3), and that
decision is only as fine as the block. At 100 kb a scan skips in 100 kb units; at 1,000 kb it cannot
skip anything smaller. **How much that costs is unmeasured**, and it is the same measurement that
decides §12 question 1 — how much of a cohort scan the index summary alone can serve.

**So 1,000 kb is live and may well be better.** It was smaller on both samples measured — 4.626
against 4.627 on tomato, which is nothing, and **16.444 against 17.557 on the patchy human one,
which is 6 %**. The case for moving rests entirely on that second figure, and on the skip measurement
coming back saying blocks are rarely skipped.

### 4.2 The look-back window — declared, and the reader honours it

**The window is what a reader must hold, so it is a per-file number the reader reads before it
allocates.** zstd also records it in each frame's own header, so a reader *could* learn it there —
but by then it is about to allocate, and a mismatch between our declaration and the frame's is a
corruption worth detecting rather than tolerating.

**A reader configures its decoder from the declared value and refuses a file that exceeds its
budget with a clear error.** This is not defensive padding: zstd refuses a frame whose window
exceeds the decoder's configured maximum, so a reader that assumes 32 kB will reject a perfectly
good file written with more, and the error it produces says nothing useful. *"This file needs a
512 kB window and this reader is configured for 32 kB"* is an actionable message; a zstd error code
is not.

### 4.3 How many streams a block is cut into

**Declared, and it is the one parameter with a hard ceiling — see §5.2 for the measurement.** A
block may be one compressed stream carrying whole records, or two — the records, plus a small one
carrying the single number a cohort scan reads per record so a scan need not decompress records at
all.

**A reader opens only the streams it needs**, which is what keeps the cost a function of what the
reader is doing rather than of the file. A scanning reader opens the small one; a materialising
reader opens the records one; a reader doing both at once pays for both.

### 4.4 The manifest

The header carries, for each field of a record: its name, its cardinality, its encoding, and that
encoding's parameters. Encodings come from a **fixed, named menu** — plain variable-length integer,
zig-zag variable-length integer, fixed-point-with-a-declared-step, difference-from-previous-with-
reset, raw bytes. Production has the shape of this already (its `[[column]]` array,
[`src/psp/header.rs:24`](../../../../src/psp/header.rs), typed by
[`src/psp/registry.rs`](../../../../src/psp/registry.rs), which already distinguishes a column's
cardinality from its element type).

**A fixed menu rather than an open-ended scheme, deliberately.** An arbitrary plug-in encoding would
buy flexibility nobody has asked for and cost speed on the hot path — decoding runs at about
20 million records a second, and every field goes through it.

**What this buys beyond the owner's rule: a reader can skip a field it does not recognise**,
provided the manifest says how to measure its length. That is a better compatibility story than a
version number alone, and it comes free.

### 4.5 How records are serialised — by hand, not through `serde`

**Decided: the record encoder and decoder are written by hand against the manifest.** Three reasons,
in order of weight:

1. **A `serde` derive fixes the layout at compile time, and §4.4 requires it to be a per-file
   choice.** Expressing *"this value is a float stored as a variable-length integer count of a step
   declared in the header"* needs a custom serialiser — at which point the encoder has been written
   by hand and a `serde` layer added on top of it.
2. **Speed.** Decoding is the hot loop, measured at about 20 million records a second in the
   prototype, and a generic serialiser gives up borrowing directly from the decompression buffer.
3. **Compatibility.** `bincode`-style formats have no field identity, so a reader meeting an
   unfamiliar field misparses silently. The manifest handles this properly (§4.4).

**`serde` is right for the header** — it already is, through TOML — **and for the trailer's
payload**, which is small, written once and read once. Spending the convenience there costs nothing
on any hot path.

---

## 5. What a reader may hold, and what it costs

### 5.1 The reader's contract — two conditions, and the second is the one that gets lost

**A reader must not hold an amount that depends on the block size.** Two things are required and
only the first is the compressor's doing:

1. **Do not inflate the whole block.** Pull decompressed bytes out incrementally, bounded by the
   declared window.
2. **Do not accumulate what you inflated.** Parse a record out of the rolling buffer, hand it to the
   caller, keep nothing.

**Satisfying only the first moves the memory rather than removing it.** In production's cohort run
the assembled per-sample columns are the largest single mass of the heap — 29.6 %, larger than the
decompression buffers they came from
([the memory review](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md) §2). A reader that
streams a block in perfectly and then gathers every record into an array before returning has
achieved nothing.

So one open sample holds: the declared window, one read buffer, one rolling decompressed buffer, the
decoder's own state, and the record being built. **None of those is a function of the block size.**

### 5.2 What it costs, measured — and the ceiling on stream count

62 tomato accessions opened at once and advanced one record each per round, the shape a cohort merge
reads in. `rolling` and `read chunk` are the reader's two buffers, which are its own choice and not
the file's:

| reader buffers | 1 stream | 2 streams |
|---|---:|---:|
| 64 kB each | 346 kB | 691 kB |
| 16 kB each | 250 kB | 501 kB |
| **4 kB each** | **227 kB** | **453 kB** |

**And it barely moves with depth**, which goal 1 asserts and this is the evidence for. The same walk
over 62 copies of a sample at **279 reads a position** — a hundred times the tomato cohort's depth:

| | per open sample | R² |
|---|---:|---:|
| tomato, 3 reads a position | 346.3 kB | 1.0000 |
| HG002, 279 reads a position | 367.1 kB | 1.0000 |

**6 % more memory for 93 times the depth**, both fits straight to four digits. That is what the
design predicts and the reason it predicts it: what a reader holds is its decoder and its two
buffers, none of which is a function of the data. The 6 % is the record being built, which grows
from 4.6 to 16.3 bytes.

*Both are one sample's characteristics replicated, so this measures the reader's cost and not the
variety of a real cohort.*

Three things follow:

- **The cost is per live decoder and it is linear.** The second decoder costs the same as the first,
  as it should for the same object.
- **About 190 kB of it is zstd's own context state**, which no buffer choice reaches. That is the
  floor for one stream, and §5.3 explains where it comes from and why we cannot currently move it.
- **Against the 500 kB budget: one stream fits at any buffer size, two fit only with small buffers,
  and four do not fit at all.** So **one compressed stream per field — which is what a columnar
  layout amounts to — is ruled out on memory grounds**, not on preference. That is worth stating
  plainly because the columnar shape is what production uses and what a port would reach for.

*Measured with `--streams N`, which opens the same store N times per sample. A stream's buffers are
the same size whatever it carries, so this measures the multiplier without needing the streams to
differ — which is the quantity in question. It does **not** measure how much a real second stream
would add to the file, which is [`psp_record_encoding.md`](psp_record_encoding.md) §2.3's.*

### 5.3 Where the 190 kB floor comes from, and why it is not currently reachable

**It is not the look-back window.** A 32 kB window explains 32 kB of it. The rest is zstd's own
**internal block** — its private subdivision of a compressed frame, at most 128 KiB
(`ZSTD_BLOCKSIZE_MAX`, defined as `1<<17` in the vendored `zstd-sys` C headers) — which the decoder sizes
its working space from. 32 kB of window plus 128 KiB of block staging plus the entropy tables is
about the 190 kB measured, which is why the floor sits where it does.

**zstd 1.5.7 can be told to make that unit smaller**, on both sides: `ZSTD_c_maxBlockSize` when
compressing and `ZSTD_d_maxBlockSize` when decompressing, the second being the one that actually
reduces what a decoder allocates.

**We cannot reach either from Rust as the dependency currently stands, and this was checked rather
than assumed.** `zstd-safe` exposes `CParameter::MaxBlockSize` only behind its `experimental`
feature, which this project does not enable; and it does not wrap the decompression-side parameter
at all — `DParameter` carries `WindowLogMax` and three other experimental entries, and none of them
is this one.

**So the 190 kB is a property of the bindings, not of the format.** Recorded because it decides
§12 question 1: if the floor could be brought down to, say, 60 kB, two parts would fit the budget
comfortably and the light/heavy split would stop being a trade. **Settled by:** enabling
`zstd-safe`'s experimental feature and calling `ZSTD_DCtx_setParameter` through `zstd-sys` for the
decompression side, then re-running the sweep in §5.2. Neither is difficult; both are a dependency
change and unsafe FFI for a question that is not blocking.

### 5.4 Against production, on the same records

| | production `.psp` | this format | |
|---|---:|---:|---|
| per open sample | 2.614 MB | 0.338 MB | **7.7× less** |
| bytes a record | 8.188 | 5.356 | 35 % smaller |
| cohort on disk (62 samples) | 3.52 GB | 2.38 GB | 32 % smaller |
| blocks per sample | 1,674 | 154 | index 11× smaller |
| 62-sample walk | 42.4 s | 23.1 s | 1.8× faster |
| records read | 471,520,156 | 471,520,156 | identical, same checksum |

Both memory figures are the slope of a straight-line fit over cohorts of 1 to 62, R² 0.9998 and
1.0000. *Extrapolated to three thousand samples that is 7.66 GB against 0.99 GB — arithmetic from a
fit that stops at 62, not a measurement.*

---

## 6. The interface

**What a user of this format wants to do**, and what the container offers. The types are the
implementation's; this fixes the operations and their guarantees.

### 6.1 Reading

- **Open a file.** Reads the footer, then the index, then parses the header. No block is touched.
  Fails if the footer is absent or the version is unknown.
- **Get the header** — available immediately after open, including the manifest.
- **Get the trailer** — one seek and one read.
- **Walk the blocks without decoding them** — chromosome, span, and the scan summary, straight from
  the index. **This is what lets a cohort skip regions where no sample varies**, and it is most of
  why reading psp files is cheaper than re-reading alignments.
- **Stream records from a chosen block, the first by default.**
- **Stream records from a genomic position** — the index turns a coordinate into a block with one
  lookup. Users think in coordinates, so a position-based entry point belongs here even though a
  block-based one is what it is built on.

### 6.2 Writing

- **Create a file, given the header content.** The manifest is fixed at this moment; a writer cannot
  change a field's encoding half-way through a file.
- **Push records**, in coordinate order. **The writer rejects an out-of-order record** — coordinate
  order is what the index and every seek depend on, and a writer that accepts a stray record
  produces a file that seeks wrongly rather than one that fails.
- **Close, with a trailer payload that may be empty.** Closing writes the index, the trailer and the
  footer, in that order. **Before the close there is no footer, so any reader correctly refuses the
  file** — which is exactly what should happen to a killed run.
- **Reopen an existing complete file to add more records.** The footer's index offset says where the
  blocks end; the writer truncates there, discarding the old index, trailer and footer, and carries
  on. The header — and therefore the manifest — is not rewritten, so the appended records must use
  the encodings already declared.

**Closing must be durable and this is easy to get wrong.** Flush the format, then surface the
buffered writer's errors, then sync. A `BufWriter` dropped without that can swallow a failed flush,
and a truncated footer on a billions-of-records file looks exactly like an interrupted run —
production's writer carries this warning in its own doc comment
([`src/psp/writer.rs:694-702`](../../../../src/psp/writer.rs)) because it has bitten before.

---

## 7. Cross-cutting concerns

**Memory.** The subject of the document; §5 measures it. One number to carry: **an open file is 227
to 346 kB depending on the reader's buffer choices, and the budget is 500 kB.**

**Errors.** Four classes, and they want to be distinguishable because they mean different things to
whoever sees them:

| what happened | what the user has to do |
|---|---|
| no valid footer | the run was interrupted — rebuild the file |
| unknown format version | upgrade the reader |
| declared window exceeds the reader's budget | raise the budget, or rewrite the file |
| a block fails to decompress, or a record runs past its block | the file is corrupt |

None of these may reach a caller as a half-built record.

**Concurrency.** Blocks are independently decodable, so nothing here serialises. Each open sample has
its own decoder state; two readers of one file share nothing but the bytes on disk.

**Byte-identity.** The block cut is a function of the coordinate and the observation stream alone
(§4.1), so the same sample gathered at any worker count gives the same file apart from any timestamp
in the header.

---

## 8. Traps — what will bite the coder

- **The word *trailer* means two different things** in this document and in production's code (§2).
  Production's `trailer.rs` is this document's *footer*; production's *metadata section* is this
  document's *trailer*.
- **⚠ Streaming the block is half the job; not accumulating the records is the other half**, and it
  is the half that gets lost (§5.1). A reader that streams perfectly and then gathers records into
  an array has moved the memory, not removed it.
- **A record can straddle the rolling buffer, and the parse must be restartable.** The buffer holds
  whatever has been decompressed so far, which may be the first half of a record. Running out of
  bytes has to be an answer the parser can give — not a panic, not a short read — and the retry must
  resume from the record's *start* with the running position, coverage and read-name bases restored.
  **A parse that half-advances that state before failing corrupts every record after it, plausibly.**
- **A single record can exceed the rolling buffer** — many alleles, many read names. The buffer must
  grow; a fixed maximum record size is not a safe assumption.
- **Every running difference resets at a block boundary** (§3.2), and a block that forgets one reads
  back wrong from its first record without failing.
- **A missing footer means refuse, not read short** (§3.3). This is the only thing standing between a
  killed pileup and a silently truncated sample.
- **Do not size any buffer from the block.** A block header that carries its uncompressed length is a
  temptation to allocate it; the whole design is that a reader never does.
- **`--target-dir` on any heap-profiled run.** Without it the instrumented binary overwrites the
  release one at the same path and every later "release" run silently re-executes it five to six
  times slower. That has already cost one performance review its whole measurement set.

---

## 9. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| header framing | [`src/psp/header.rs`](../../../../src/psp/header.rs) | the pattern as-is — magic, length prefix, TOML body, sentinel — so `head` still works on an ng psp |
| the field manifest | the `[[column]]` array, [`src/psp/header.rs:24`](../../../../src/psp/header.rs) + [`src/psp/registry.rs`](../../../../src/psp/registry.rs) | the shape; §4.4 widens it to carry encoding parameters |
| footer | [`src/psp/trailer.rs`](../../../../src/psp/trailer.rs) | the layout and the tail-magic-last trick; widened for the trailer's offsets |
| block index | [`src/psp/index.rs`](../../../../src/psp/index.rs) | the flat vector, unchanged — §3.3 says why it no longer needs replacing |
| variable-length integers | [`src/psp/varint.rs`](../../../../src/psp/varint.rs) | as-is: LEB128 and zig-zag LEB128, specified and tested |
| the compression seam | `new_column_compressor`, [`src/psp/block.rs:718`](../../../../src/psp/block.rs) | the shape — one long-lived compressor per writer, frame checksums on. The window cap is new |
| the streaming reader and writer | `examples/psp_row_stream_roundtrip.rs` | the working prototype every number here was measured on; the parity oracle below is its `verify-streaming` |

**Parity oracle.** The prototype's `verify-streaming` walks the new store and the `.psp` it was
written from in lockstep and fails on the first record that disagrees — 7,687,686 records, every
integer field, allele sequence and read-name list compared exactly, the approximated fields against
their own step. That is the model: a deliberately simple decoder used only by tests, against which
the real reader is compared record for record.

---

## 10. How we know it works

1. **Round-trip, with the right strictness per field.** Integer fields, allele sequences and
   read-name lists identical; approximated fields inside their own step. A round-trip that compares
   whole records with a blanket tolerance will pass while a read-name list is being corrupted.
2. **Restart equals sequential.** Reading from an arbitrary block gives exactly the records a full
   sequential read gives from that point. This is the test that catches a running difference that
   was not reset (§3.2), and it cannot be skipped because the failure is plausible rather than loud.
3. **The per-open-file budget, measured rather than argued.** N samples open and walked in lockstep,
   peak resident reported, against the 500 kB of goal 1. The prototype does this at 1, 8, 31 and 62.
4. **Worker-count invariance.** One sample gathered at 1, 2, 4, 8, 16 workers gives byte-identical
   files apart from the header's timestamp.
5. **An interrupted write is refused.** Kill a writer before close; the reader must reject the file,
   not read the blocks that made it to disk.
6. **The cohort scan still saves what it saved.** Scanning the index over a segment must skip the
   blocks it should skip; a file that round-trips but forces every block to be inflated has passed
   test 1 and failed the design.

### 10.1 ⚠ The end-to-end oracle cannot be a byte-identical VCF

**Measured, not feared.** Changing only how records are grouped into blocks changes the order in
which per-sample evidence is summed, and floating-point addition is not associative. On the
50-accession tomato cohort, rewriting the same records at a 20 kb and an 80 kb block span left the
set of sites and alleles **identical to the line** and returned **1,194 records of 180,366 with a
different QUAL** — median difference 0.010 Phred, largest 4.48. No genotype, GQ, DP, AD or INFO field
differed anywhere.

**One site in 180,366 crossed the emission gate as a result**: `SL4.0ch10:58265030 T>A` sits at
QUAL 30.075 against a minimum of 30, and is emitted at 20 kb and 80 kb but not at 5 kb. Seventeen of
the differing records sit below QUAL 40, within reach of the gate.

**So the oracle is: the site list, genotypes, GQ, DP and AD exactly; QUAL within a tolerance; and an
explicit statement about what happens at the gate.** A test demanding byte-identity fails on a
correct implementation. A test using a tolerance everywhere passes while a read-name list is being
corrupted.

---

## 11. Deferred, with a recommended home

- **Which fields a record has and how each is encoded** —
  [`psp_record_encoding.md`](psp_record_encoding.md), which this document's manifest carries the
  declarations for.
- **The read names' encoding** — [`psp_chain_id_encoding.md`](psp_chain_id_encoding.md).
- **The cohort reader's scheduling** — which blocks to fetch, in what order, with how much
  look-ahead — [`run_streaming.md`](run_streaming.md).
- **Correcting [`run_streaming.md`](run_streaming.md) §7.2**, whose "tens of kilobytes" is
  superseded by the 500 kB budget. That document's owner should make the change; it is the sentence
  that made a columnar shape look impossible, and leaving it will send the next reader down that
  path.
- **Whether a record carries its own reference bases.** Today it does, so every sample carries a copy
  of the reference over its footprint. The leaning is to drop it and re-fetch; nobody has timed the
  re-fetch. [`run_streaming.md`](run_streaming.md) §11 question 4.

---

## 12. Open questions

1. **How many streams does a block have — one or two?** — OPEN, and §5.2 gives it a hard ceiling
   rather than an answer. One stream is 227–346 kB per open sample and two are 453–691 kB, so two
   fit the budget only with small buffers. What is not measured is what the second stream *buys*: how
   often a cohort scan can skip a block using only the index summary (§3.3), which may make a
   per-record scan stream unnecessary. *Leaning:* start with one stream and the index summary, and
   add the second only if the skip rate is measured to be poor. **§5.3 is the other way this could
   move:** the per-decoder floor is set by a zstd parameter the Rust bindings do not currently
   expose, and lowering it would make two parts cheap. **Settled by:** counting, on a real
   cohort segment, how many blocks the index summary alone lets a scan skip.
2. **What byte ceiling, if any, should a writer put on a block?** — OPEN (§4.1), and it costs
   nothing measurable so far: a 100 kb grid with a 1 MiB ceiling gives 4.628 bytes a record against
   4.627 without, on tomato. *Leaning:* offer it, default it off, and let the first whole-genome
   deep-coverage run set it — at 279 reads a position a fully covered 100 kb block is about 1.6 MB,
   which is a large thing to hold while writing. **Settled by:** the block-size distribution on a
   whole-genome deep-coverage sample, which nothing here has produced. *Not* by seek time: the owner
   ruled on 2026-08-25 that reading a file start to end is the common case.
3. **How badly do near-empty blocks compress on a patchy sample?** — **ANSWERED 2026-08-25 (§4.1):
   about 10 % of the file.** A 5 kb grid on a sample with 74,623 covered positions scattered over 644
   regions gives blocks of about 119 records and costs 18.242 bytes a record against 16.444 at
   1,000 kb. **So the rule that accumulates across empty spans ships**, and the question that
   replaces it is what threshold it uses. *Still unmeasured:* a genuinely thin whole-genome sample —
   1× rather than a region-restricted one — where the gaps are between positions rather than between
   regions.
4. **What does the reader's buffer pair cost in speed as it shrinks?** — PARTLY MEASURED. Smaller
   was slightly *faster* at one stream (31 s against 35 s over 471 M records), which is the opposite
   of the expected direction and is not understood. *Leaning:* default to 16 kB, which is inside
   budget at two streams and near the fast end. **Settled by:** repeating the sweep on a machine that
   is not sharing its cores, which this one was.
