# ng — the psp file format: the container

*Design spec, 2026-08-25. **No code for the format itself**; a working writer and reader for
the compression scheme exist as a measuring program (`examples/psp_row_stream_roundtrip.rs`) and
every number here was taken from it or from production.*

*This document owns the **container**: the sections a file is made of, how a reader finds them,
what a reader is allowed to hold, and what the module offers its callers (§6). What goes
**inside** a record — which fields, how each is encoded, what may be approximated — is
[`psp_record_encoding.md`](psp_record_encoding.md), and the chain ids are
[`psp_chain_id_encoding.md`](psp_chain_id_encoding.md). The requirement both serve is a
per-open-sample memory budget, stated in §1.1.*

*Downstream: [`../arch/psp_file_format.md`](../arch/psp_file_format.md) is the code shape,
[`../impl_plan/psp_file_format.md`](../impl_plan/psp_file_format.md) the build order. Companion
measurements: [the memory review](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md).*

---

## 1. What this is

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
   touched.)* **The independence from depth is measured, not assumed — §5.2.** The one thing that
   *can* push a reader past this is a single record larger than its buffer, which §8 forbids fixing
   a limit on in the format; §4.4 gives the reader its own ceiling instead, set so that three
   thousand readers at it come to the budget rather than a multiple of it.
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
  user's, per the owner's ruling of 2026-08-25. §4.5 says how they travel.
- **It does not change production's `.psp`.** ng replaces it; production stays as it is. Several
  findings here apply to it unchanged and that is not a reason to touch it.
- **It does not define the cohort reader's scheduling** — which blocks it fetches, in what order,
  with how much look-ahead. That is [`run_streaming.md`](run_streaming.md)'s.
- **It does not fix the module's Rust types or names.** §6 fixes the operations a caller can
  perform and what each guarantees and refuses; the signatures there are indicative and the
  architecture doc may name things differently.
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
| **psp block** | a run of consecutive records over one span of reference, compressed as one unit | block |
| **block index** | one entry per psp block, so a reader can seek | block index |
| **trailer** | a payload the writer supplies when it closes the file; may be empty | **metadata section** |
| **footer** | fixed-size, at the very end; offsets and magic. Its presence means the file is complete | **trailer** — 32 bytes there ([`src/psp/trailer.rs:28`](../../../../src/psp/trailer.rs)); this one is wider, §3.3 |

**⚠ zstd has its own "block", and it is not ours.** Three things end up called a block or a frame,
so this document uses all three names in full wherever either could be meant:

| term | what it is | who chooses its size |
|---|---|---|
| **psp block** | a span of reference and the records in it | **us** — the `genomic_block_size_bp` of §4.1 |
| **zstd frame** | one independently decompressable compressed unit. **One per psp block** | follows the psp block |
| **zstd block** | zstd's own subdivision of a frame, at most 128 KiB | **nobody here** — zstd makes them as it works |

The zstd block is not a knob and appears in this document exactly once as anything but a caution:
it is where most of the per-decoder memory floor comes from (§5.3). Everywhere else, an unqualified
"block" in this document means a psp block.

Two more terms this document needs:

- **⚠ "reach" is not used in this document for anything to do with compression.**
  [`cohort_merge.md`](cohort_merge.md) defines it as *the last reference position an observation
  covers*, `start + span − 1`, and the word and the arithmetic are production's. An earlier draft
  here used it for the compressor's look-back distance; that was a collision and it is gone.
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
- A **chain id** identifies the DNA fragment a piece of evidence came from: **one `u64` per read
  *pair*, not per read**, mates collapsed onto a single id, allocated in order and never reused
  ([`chain_id_allocator.rs`](../../../../src/ng/locus_generation/pileup/chain_id_allocator.rs)).
- **the record head** — the fixed fields at the front of every record that let a reader decide
  whether it wants the record without building it: the position offset, the reference span, the
  non-reference read count, and the body's length (§4.3).

---

## 3. The file

```
+---------------------------+
| header                    |  plain text: magic, length, TOML body, sentinel
+---------------------------+
| psp block 0               |  records, one compressed stream, each record
| psp block 1               |    opening with the head of §4.3
| ...                       |
+---------------------------+
| block index               |  one entry per psp block
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
- **the manifest** (§4.5) — the block cut rule, the look-back window, the stream layout, and each
  field's encoding with its parameters.
**One further field is wanted by the cohort merge and is not this document's to add.**
[`cohort_merge.md`](cohort_merge.md) asks the header to record the **observation reach ceiling** —
the widest reference span any observation in the file can have, the generator's `max_record_span`
([`src/ng/locus_generation/pileup/generator.rs:93`](../../../../src/ng/locus_generation/pileup/generator.rs))
— and routes it to [`run_streaming.md`](run_streaming.md) §6.1, which owns the header's contents. It
is flagged there and not yet written.

**It is worth knowing what it is and is not for, because it is easy to over-read.** It is *not* an
artefact of the columnar layout: it is a fact about how wide a record can be, which is true whatever
the layout. And it is *not* needed for correctness — that document says so in as many words: *"no
refusal accompanies it, and none is needed"*. What it buys is that a reader can **size** its
observation cache up front, taking the maximum over the cohort's files, instead of growing it. A
forward reader learns each record's span from the position summary as it goes and never needs the
ceiling at all.

### 3.2 Blocks

A block is a run of consecutive records from one chromosome, compressed as one unit, with the
look-back window capped at the value the header declares.

**Every block is self-contained.** It opens with its chromosome, its first position and its record
count, and **every running difference inside it restarts** — the position difference, the coverage
difference, the chain-id difference. A block never crosses a chromosome. So the restart points are
the block boundaries and there is no separate seek mechanism to build.

**This is the property goal 2 rests on, and the one most easily broken by accident.** A block that
fails to reset a running difference reads back wrong from its first record — and *plausibly* wrong,
because coverage is smooth and a position difference that is slightly off still parses.

### 3.3 The block index and the footer

**One entry per block**, in genomic order: the chromosome, the first position, and the byte offset.
**And nothing else.** An earlier draft of this section had each entry also carry the largest
non-reference support anywhere in the block, so a reader could decide whether to touch a block
without decompressing it. At a 100 kb block with roughly one varying position in a hundred, every
block holds about a thousand of them, so that field never says *skip*: it cost index bytes and
writer bookkeeping for a decision that never fires, and it is gone (the owner, 2026-08-25 —
[`psp_record_encoding.md`](psp_record_encoding.md) §2.4 carries the reasoning). *This paragraph
described it as present until 2026-08-26; the removal had been agreed and was not finished here.*

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
them is the writer's business.

**And what is in them is deliberately not settled** (the owner, 2026-08-25). The per-sample summary
is the case that exists today; the coverage-against-GC histogram and the census are expected to join
it. **What the set finally is depends on the statistical work, which is not finished**, and the
container is built so that it does not have to be: because the payload is opaque here, adding to it
or reshaping it is a writer-side change and **not a container version bump**. That is the property
this section exists to protect, and it is worth more than any particular list would be.

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
compressor starts cold each time. **So the writer may declare one secondary rule** — close a block
early once it exceeds a byte ceiling. It is the user's choice and it is recorded in the manifest.

**A grid cell holding no records produces no block.** The rule cuts where a block *ends*; it never
asks for one per 100 kb of reference. A sample covering two cells ninety apart writes two blocks, so
a thin sample pays no index entry and no compressed frame for reference it did not cover.

**⚠ Blocks that are merely *small* are not merged, and this section used to say they were.**
An earlier version offered a second secondary rule — accumulate across empty spans so a patchy
sample gets one large block instead of several thin ones — and §12 question 3 recorded it as
shipping, with only its threshold left open. **The owner ruled against it on 2026-08-27: merging
would complicate the alignment between samples**, which is the first thing this section lists as
what the grid buys. Merge, and one sample's block may begin ninety cells before its neighbour's, so
which block holds a given position differs from sample to sample and a cohort reader wanting one
position decodes from far behind it.

**What merging would have saved was never measured.** The figures below compare *block sizes* on the
same non-merging writer; no writer that accumulates across empty spans has been built. If someone
later wants the price, they have to build it to get it.

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

**On patchy coverage a 5 kb grid costs 10 % against a 1,000 kb one**, because blocks of about 119
records make the compressor start cold each time. The human sample here is 74,623 covered positions
scattered over small regions. **Its patchiness is an artefact of the region-restricted pileup that
produced it, not a property of a 279× sample** — but patchy coverage is real for exomes, panels and
thin samples. **The answer is the block size, not merging** (owner, 2026-08-27): the 100 kb default
already recovers most of that penalty, 17.557 bytes a record against 18.242 at 5 kb, and 1,000 kb
recovers the rest — which is the live question below, not a second cut rule.

*⚠ The "644 regions" in the sentence above is not re-derivable from the file: the corpus gives 1,217
maximal runs of consecutive covered positions in 281 occupied 100 kb cells. The figure is kept as
written because the measurements in the table were taken against it; the count itself carries no
weight in the argument.*

**Why 100 kb and not something else.** It is a round number that satisfies goal 2 directly, sits in
the flat part of the tomato curve, and recovers most of the patchy-data penalty (17.557 against
18.242 at 5 kb). **It is not an optimum** — 1,000 kb is smaller on both samples — and it is a
starting value, not a derived one — and **the argument that was holding it down has been
withdrawn.** I had reasoned that larger blocks cost seek time, since a reader starting mid-block
decodes from that block's beginning. **The owner's ruling of 2026-08-25: seeking is not the common
case — a run reads a file from the beginning to the end.** So seek cost should not set this number.

**What argues against a very large block is where a reader may start, and nothing else.** A reader
seeks to a block and streams from its beginning (§1.2), so the block size *is* the granularity of
restarting: at 100 kb a reader that wants one position decodes at most 100 kb of records to reach it,
and at 1,000 kb, ten times that.

*An earlier version of this paragraph argued instead from block skipping — that the index's per-block
summary lets a scan decide whether to touch a block at all, and that the decision is only as fine as
the block. **That argument is void: the summary was removed** (§3.3), because at 100 kb essentially
every block contains something and the field never fires. Skipping now happens a record at a time,
through the record head of §4.3, and a record head does not care how large its block is.*

**So 1,000 kb is live and may well be better.** It was smaller on both samples measured — 4.626
against 4.627 on tomato, which is nothing, and **16.444 against 17.557 on the patchy human one,
which is 6 %**. **The case for moving now rests on that second figure alone**, against a tenfold
coarser restart granularity; the skip measurement the earlier version made it conditional on does
not exist to be taken, since there is nothing left to skip at block grain.

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

### 4.3 The record head — how a reader judges a record without building it

**Settled 2026-08-25, and it replaces a two-stream design that was adopted and then measured
away.** A psp block is **one** compressed stream. Each record in it opens with a fixed head that
answers, cheaply, the two questions a reader has before it decides whether it wants the record:

```
record = position_offset | reference_span | non_reference_reads | record_length | body
         └────────────────── the head ──────────────────────┘   └── skip this ──┘
```

A reader takes the head, decides, and either builds the body or advances `record_length` bytes past
it. **Nothing else in the block has to be touched to make that decision.**

#### Why a head and not a second stream

**Because the alternative does not save what it appears to.** The obvious design — and the one this
document carried until it was measured — puts those numbers in their own separately-compressed
stream so a scan can read them without touching records at all. It is what production's two-phase
decode does with its columns.

**A reader still has to walk the records.** It cannot seek to record *N* inside a psp block: the
block is one zstd frame, it comes out sequentially from its start, and the decompressed bytes carry
no separators — so finding where record *N* begins means decoding every variable-length integer in
the records before it. A second stream therefore adds its own walk **on top of** that one rather
than replacing it, and it costs a second decompressor.

Timed on a tomato accession at three reads a position, 7.69 M records, one stream:

| what a walk does | time |
|---|---:|
| decompress only, nothing else | 0.104 s |
| + walk each record's bytes to find where it ends | 0.163 s |
| + build the record objects — a full walk | 0.30 s |

**So the three costs are 0.104 s of decompression, 0.059 s of finding record ends, and 0.137 s of
building records.** Only the last is avoided by knowing a record is unwanted. The middle one is
avoided *only* by a record saying how long it is, which is what `record_length` is for.

| design | a walk keeping one record in a hundred | file | memory at 5,000 samples |
|---|---:|---:|---:|
| no head — build every record | 0.29 s | 4.628 | 1.14 GB |
| **head with `record_length`** | **0.141 s** | **5.056 (+9.2 %)** | **1.14 GB** |
| a second stream | ~0.204 s *(composed)* | +1.7 % | 2.27 GB |

**All but the second-stream row are now measured on a reader that exists.** The skipping walk is
**2.06× faster** than a full one; reading heads and skipping bodies costs 0.027 s over 7.69 M
records, and the rest of its 0.141 s is decompression, which nothing avoids. **The second stream is
last on both axes**, which is why it is gone.

**⚠ The head costs 9.2 % of the file at three reads a position and 5.8 % at 279 — not the 1.4 % an
earlier draft quoted**, and the difference is the point of the next paragraph.

#### Why a skippable body costs more than a length field

**The body has to stand on its own, and that is most of the price.** A record's coverage and its read
names are normally coded as differences from the *previous record* — which is cheap, because
neighbouring positions are alike. **A reader that skips a body never sees those differences, so it
loses the base both are measured from and every record after it decodes wrong.** So with a head,
both restart at each record: the coverage is absolute, and the chain-id difference is measured
from zero rather than from the record before.

That is what turns a 1.4 % length field into a 9.2 % head. **It is a real cost of skippability, not
an implementation detail**, and it is why the earlier figure was wrong rather than merely imprecise.

*Not measured: a head that carries the two differences itself — the coverage step and the advance in
chain-id numbering — so the body can keep its cross-record coding and a skipping reader still track
both. It would move those fields rather than duplicate them, so it should recover much of the 9.2 %,
at the cost of two more head fields and a reader that maintains state while skipping.*

#### What the head carries, and why each field

- **`position_offset`** — the distance from the previous record. Every running difference restarts at
  a block boundary (§3.2).
- **`reference_span`** — how many reference bases the record covers. **Required rather than chosen**:
  [`cohort_merge.md`](cohort_merge.md) names it as one of two things it asks of this document, because
  a record widened by a deletion covers more than one position, so a reader indexed by position
  cannot work out what a record reaches from its start alone.
- **`non_reference_reads`** — the reads at this position that supported something other than the
  reference. **The owner's correction of his own first suggestion, 2026-08-25**: a count of
  alternative *alleles* answers *does anything vary here* identically, since an allele exists only
  because reads showed it — but the read count also lets a reader apply a threshold.
- **`record_length`** — the body's length in bytes, so an unwanted record is skipped rather than
  decoded. Measured cost: **1.4 % of the file** at three reads a position, 3.3 % at 279.

Together these are what [`cohort_merge.md`](cohort_merge.md) calls the **position summary** —
*"the cheap facts a builder needs about a position before it decides anything"* — and the name is
that document's, not this one's.

**Fixed-width or variable-length is the manifest's to say** (§4.5), not this section's. A fixed
width is quicker to read, and costs less than it looks after compression because a column of small,
repetitive values collapses — the four head fields together compressed to 0.077 bytes a record when
measured on their own. *Unmeasured: the two encodings against each other in place.*

### 4.4 The reader's two buffers — 16 kB each

A reader holds a buffer of compressed bytes read from the file and a rolling buffer of decompressed
bytes it parses records out of. **Both are the reader's choice, not the file's.**

**Bigger is not faster, which is the opposite of what one would assume.** Re-measured 2026-08-26 on a
machine with nothing else running, median of three runs — a walk over 7.69 M records keeping one in a
hundred, and the same store held open 62 times:

| buffers | walk | per open sample | at 5,000 samples |
|---|---:|---:|---:|
| 4 kB | 0.149 s | 233 kB | 1.11 GB |
| **16 kB** | **0.143 s** | **257 kB** | **1.23 GB** |
| 64 kB | 0.161 s | 353 kB | 1.69 GB |
| 256 kB | 0.200 s | — | — |

**There is an optimum rather than a trend, and it is near 16 kB.** Going up from there costs both
memory and time: 64 kB is 13 % slower than 16 kB and 256 kB is 40 % slower. Going down costs a
little time and saves a little memory.

**So 16 kB, and there is no trade to weigh.** *Why the curve turns is not established. Buffers
falling out of cache is the obvious guess and it is a guess; nothing here measured it.*

**Against the 500 kB budget an open sample is 257 kB**, which leaves room to spend later if something
justifies it.

#### The rolling buffer's ceiling — the third number, and the only one a corrupt file can move

The rolling buffer starts at 16 kB and **grows for a record that does not fit**, because §8 says a
maximum record size is not safe to assume. On a well-formed file that growth is bounded by the
largest record the caller can produce, and it shrinks back at the next block. **On a corrupt file
nothing bounds it**: no record parses, so nothing ever fits, and the buffer doubles until the frame
runs out — and a block's *decompressed* size is not bounded by its size on disk. Measured, a
4,132-byte block drove a reader with no ceiling to hold **67,125,248 bytes**.

**So the reader carries a ceiling on how far the rolling buffer may grow for one record, and it is
the reader's budget rather than a maximum record size the format fixes.** Two numbers set the
default of **512 kB** and both are measured:

- **It is ten times the largest record this caller's own depth cap can produce.** At three hundred
  reads a position, one observation each, a record encodes to 18,292 bytes over a 50-base span and
  48,693 over 150. *(Milestone E's chain ids are not in those figures — the encoder drops them
  today — and the number wants re-measuring when they arrive.)*
- **It keeps the worst case inside goal 1's budget.** Three thousand open samples every one of which
  met a damaged block at the same moment is 1,572,864,000 bytes, against the 1.5 GB §1.1 gives three
  thousand — the same number to within 5 %. At 1 MiB, the value first proposed, it would have been
  3.07 GB.

**It is a per-reader setting rather than a constant**, so a run that genuinely needs more can raise
it without a format change, and the refusal names it (§7). A ceiling in the *format* would make a
legitimate file unreadable everywhere at once, which is what §8 rules out; a ceiling in the *reader*
is the same shape as §4.2's look-back window budget.

*An earlier draft of this section recommended 4 kB on a sweep taken while other work was running on
the same machine, and that sweep's timings were not usable. The memory figures were.*

### 4.5 The manifest

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

### 4.6 How records are serialised — by hand, not through `serde`

**Decided: the record encoder and decoder are written by hand against the manifest.** Three reasons,
in order of weight:

1. **A `serde` derive fixes the layout at compile time, and §4.5 requires it to be a per-file
   choice.** Expressing *"this value is a float stored as a variable-length integer count of a step
   declared in the header"* needs a custom serialiser — at which point the encoder has been written
   by hand and a `serde` layer added on top of it.
2. **Speed.** Decoding is the hot loop, measured at about 20 million records a second in the
   prototype, and a generic serialiser gives up borrowing directly from the decompression buffer.
3. **Compatibility.** `bincode`-style formats have no field identity, so a reader meeting an
   unfamiliar field misparses silently. The manifest handles this properly (§4.5).

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

### 5.2 What it costs, measured

62 tomato accessions opened at once and advanced one record each per round, the shape a cohort merge
reads in. `rolling` and `read chunk` are the reader's two buffers, which are its own choice and not
the file's:

| reader buffers | per open sample |
|---|---:|
| 64 kB each | 353 kB |
| **16 kB each — settled, §4.4** | **257 kB** |
| 4 kB each | 233 kB |

*A second column here once priced a second compressed stream per block, at roughly double each of
these. That design is gone (§4.3), and what it would have cost is recorded there rather than here.*

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
  floor for a stream, and §5.3 explains where it comes from and why we cannot currently move it.
- **Against the 500 kB budget, one stream fits at any buffer size and four do not fit at all.** So
  **one compressed stream per field — which is what a columnar layout amounts to — is ruled out on
  memory grounds**, not on preference. That is worth stating plainly because the columnar shape is
  what production uses and what a port would reach for. It is also why §4.3 puts the cheap fields in
  a record's head rather than in a stream of their own.

*Measured with `--streams N`, which opens the same store N times per sample. This is why a second
stream is priced the way §4.3 prices it. A stream's buffers are
the same size whatever it carries, so this measures the multiplier without needing the streams to
differ — which is the quantity in question.*

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

## 6. Using the module — what a caller can do

**This section fixes the operations, what each guarantees, and what each refuses.** The Rust is
indicative: the exact types are the implementation's, and the names are a suggestion the arch doc may
overrule. What is not negotiable is the set of operations and the guarantees attached to them.

### 6.1 The five things a caller does

| | operation | what it needs |
|---|---|---|
| **read** | open a finished file and get its header, its trailer, its block list or its records | a path |
| **write** | create a file from nothing and fill it | a path and the header content |
| **append** | reopen a finished file and add more records to it | a path |
| **re-trailer** | reopen a finished file and replace only its trailer | a path and a payload |
| **inspect** | read the header of a file without opening it as a reader | a path |

### 6.2 Reading

```rust
let psp = PspReader::open(path)?;          // footer, then index, then header
psp.header();                              // the manifest included
psp.trailer();                             // one seek and one read; may be empty
psp.blocks();                              // index entries: chromosome, first position, offset
psp.records();                             // every record, from the first block
psp.records_from(chromosome, position)?;   // from the block containing that position
psp.records_from_block(n)?;                // the block-level entry point the above is built on
```

**Opening touches no block.** It reads the footer for the offsets, the index, and the header. That is
the cost a cohort pays per open sample before it reads anything, and §5 is about keeping it small.

**`records_from` exists because callers think in coordinates.** The index turns a coordinate into a
block with one lookup, and the reader then streams from that block's start — it cannot start
mid-block (§1.2), so the position asked for is where the *records* start, not where the reading
starts.

**`blocks()` is the cheap survey.** Chromosome, first position and offset without decompressing
anything. What it does *not* carry is a per-block summary of variation:
[`psp_record_encoding.md`](psp_record_encoding.md) §2.4 says why that field was removed rather than
kept, and §3.3 here says the same.

**Skipping records is the reader's decision, not a separate call.** Every record opens with the head
of §4.3, so a walk hands the caller the head and the caller says whether it wants the body:

```rust
psp.records_where(|head| head.non_reference_reads > 0)   // builds only what the predicate wants
```

This is the shape the cohort's first pass uses, and it is why a walk that keeps one record in a
hundred is about twice as fast as one that keeps all of them (§4.3).

### 6.3 Writing a new file

```rust
let mut out = PspWriter::create(path, header)?;   // the manifest is fixed here
out.push(&record)?;                               // coordinate order, enforced
out.finish(trailer_payload)?;                     // index, trailer, footer — then durable
```

**The manifest is fixed at `create` and cannot change afterwards.** A field's encoding, the genomic
block size and the look-back window are all decided before the first record and recorded in the
header. A writer that could change a field's encoding half-way through would produce a file no reader
could interpret without re-reading the header per block.

**`push` rejects an out-of-order record**, rather than accepting it and producing a file that seeks
wrongly. Coordinate order is what the index and every seek depend on.

**`finish` is what makes the file readable at all.** It writes the index, then the trailer, then the
footer. **Before it, there is no footer, so every reader refuses the file** — which is exactly what
should happen to a killed run, and is goal 3.

**⚠ `finish` must be durable, and this is easy to get wrong.** Flush the format, *then* surface the
buffered writer's errors, *then* sync. A `BufWriter` dropped without that can swallow a failed flush,
and a truncated footer on a billions-of-records file looks exactly like an interrupted run.
Production's writer carries this warning in its own doc comment
([`src/psp/writer.rs:694-702`](../../../../src/psp/writer.rs)) because it has bitten before.

### 6.4 Appending to a finished file

```rust
let mut out = PspWriter::append(path)?;   // truncates at the index; header and manifest kept
out.push(&record)?;
out.finish(trailer_payload)?;
```

**The footer says where the blocks end**, so appending means truncating at the index offset —
discarding the old index, trailer and footer — and carrying on. **The header is not rewritten**, so
the appended records must use the encodings already declared, and `append` fails on a file whose
manifest the writer cannot honour.

**It also inherits the coordinate-order rule across the seam**: the first appended record must not
precede the last one already in the file.

**A file being appended to has no footer while the writer holds it open**, exactly like a new one. An
append interrupted half-way leaves a file that every reader refuses — which is right, but note that
it has *lost* the trailer the file had before, so an append is not a safe in-place edit of a file
whose trailer matters. **Write to a new path and rename if that matters to the caller.**

### 6.5 Replacing the trailer

```rust
psp::replace_trailer(path, payload)?;   // no writer, no records touched
```

**This is why the index sits before the trailer** (§3): the trailer's offset is where the rewrite
starts, so the blocks and the index are untouched and only the trailer and the footer are written.
It is the cheap operation, and it exists because the trailer is where things computed *after* the
records land — the per-sample summary today, whatever the statistical work adds later (§3.4).

### 6.6 Inspecting without opening

```rust
psp::read_header(path)?;   // header only; no footer, no index
```

**A file with no footer still has a header**, so this works on an interrupted file where
`PspReader::open` correctly fails. That is its point: a tool that reports what a half-written file
*was going to be* needs it, and so does anything that checks a cohort's files agree on a reference
before committing to a run.

### 6.7 What the caller sees when something is wrong

The five error classes of §7, and which operation raises each:

| | raised by |
|---|---|
| no valid footer — the run was interrupted | `open`, `append`, `replace_trailer` |
| unknown format version | `open`, `read_header`, `append`, `replace_trailer` |
| the file's look-back window exceeds the reader's budget | `open` |
| a record needs more of the reader's buffer than it allows one to hold | any record walk |
| a block fails to decompress, or a record runs past its block | any record walk |
| a record out of coordinate order | `push` |

**None of these may reach a caller as a half-built record**, and none is a panic: a corrupt file is
an input, not a bug.

**⚠ `replace_trailer` joined the second row on 2026-08-28 (the owner), and it is a correction the
implementation earned.** This section originally gave that operation two refusals — no valid
footer, and the file's own bytes — on the reasoning that the trailer is opaque and the cheap
operation should read nothing it does not need. **It has to read the header, and a file was
destroyed proving it**: the footer says where the trailer starts and *nothing in the footer bounds
that below*, because the only thing that knows where the blocks begin is the header's length. A
file whose footer claimed the trailer began at byte 4 passed every check the 48 bytes can make
about themselves, and the rewrite put a trailer over the header — a 3,742-byte psp reduced to 56
bytes, reported as success. Reading the header to bound the offset means a file written by a newer
format is now refused rather than rewritten, which is the safe answer and is this row.

---

## 7. Cross-cutting concerns

**Memory.** The subject of the document; §5 measures it. One number to carry: **an open file is 227
to 346 kB depending on the reader's buffer choices, and the budget is 500 kB.**

**Errors.** Five classes, and they want to be distinguishable because they mean different things to
whoever sees them:

| what happened | what the user has to do |
|---|---|
| no valid footer | the run was interrupted — rebuild the file |
| unknown format version | upgrade the reader |
| declared window exceeds the reader's budget | raise the budget, or rewrite the file |
| a record needs more of the buffer than the reader allows one to hold | raise the reader's buffer ceiling, or the block is corrupt |
| a block fails to decompress, or a record runs past its block | the file is corrupt |

**The fourth is the newest and the one whose two readings differ most**, which is why it is its own
class rather than damage: a genuine record can be larger than any fixed budget — §8 says a maximum
record size is not safe to assume — but a corrupt block that never parses grows the buffer until the
frame runs out, and both arrive at the same line. Measured: a 4,132-byte block drove a reader with no
ceiling to hold 67,125,248 bytes. The refusal names the ceiling so the first reading has an action;
the class exists so the second is not silently reported as the first.

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
  resume from the record's *start* with the running position, coverage and chain-id bases restored.
  **A parse that half-advances that state before failing corrupts every record after it, plausibly.**
- **A single record can exceed the rolling buffer** — many alleles, many chain ids. The buffer must
  grow; a fixed maximum record size is not a safe assumption. **⚠ And the growth still needs a
  bound, because on a corrupt block nothing else provides one**: no record parses, so nothing ever
  fits, and the buffer doubles until the frame runs out. The bound belongs to the *reader*, not the
  format (§4.4) — a limit in the format would make a legitimate file unreadable everywhere at once,
  which is what this trap is about. Getting the two confused is the trap inside the trap: refusing a
  well-formed large record as damage sends the operator to rebuild a file that is fine.
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
| the field manifest | the `[[column]]` array, [`src/psp/header.rs:24`](../../../../src/psp/header.rs) + [`src/psp/registry.rs`](../../../../src/psp/registry.rs) | the shape; §4.5 widens it to carry encoding parameters |
| footer | [`src/psp/trailer.rs`](../../../../src/psp/trailer.rs) | the layout and the tail-magic-last trick; widened for the trailer's offsets |
| block index | [`src/psp/index.rs`](../../../../src/psp/index.rs) | the flat vector, unchanged — §3.3 says why it no longer needs replacing |
| variable-length integers | [`src/psp/varint.rs`](../../../../src/psp/varint.rs) | as-is: LEB128 and zig-zag LEB128, specified and tested |
| the compression seam | `new_column_compressor`, [`src/psp/block.rs:718`](../../../../src/psp/block.rs) | the shape — one long-lived compressor per writer, frame checksums on. The window cap is new |
| the streaming reader and writer | `examples/psp_row_stream_roundtrip.rs` | the working prototype every number here was measured on; the parity oracle below is its `verify-streaming` |

**Parity oracle.** The prototype's `verify-streaming` walks the new store and the `.psp` it was
written from in lockstep and fails on the first record that disagrees — 7,687,686 records, every
integer field, allele sequence and chain-id list compared exactly, the approximated fields against
their own step. That is the model: a deliberately simple decoder used only by tests, against which
the real reader is compared record for record.

---

## 10. How we know it works

1. **Round-trip, with the right strictness per field.** Integer fields, allele sequences and
   chain-id lists identical; approximated fields inside their own step. A round-trip that compares
   whole records with a blanket tolerance will pass while a chain-id list is being corrupted.
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
correct implementation. A test using a tolerance everywhere passes while a chain-id list is being
corrupted.

---

## 11. Deferred, with a recommended home

- **Which fields a record has and how each is encoded** —
  [`psp_record_encoding.md`](psp_record_encoding.md), which this document's manifest carries the
  declarations for.
- **The chain ids' encoding** — [`psp_chain_id_encoding.md`](psp_chain_id_encoding.md).
- **The cohort reader's scheduling** — which blocks to fetch, in what order, with how much
  look-ahead — [`run_streaming.md`](run_streaming.md).
- **The header field for the observation reach ceiling** —
  [`run_streaming.md`](run_streaming.md) §6.1, which owns the header's contents.
  [`cohort_merge.md`](cohort_merge.md) §13 routed it there and it is not yet written. Needed only so
  a reader can size its observation cache up front; a forward reader learns each record's span from
  the position summary and never needs it.
- **Correcting [`run_streaming.md`](run_streaming.md) §7.2**, whose "tens of kilobytes" is
  superseded by the 500 kB budget. That document's owner should make the change; it is the sentence
  that made a columnar shape look impossible, and leaving it will send the next reader down that
  path.
- **Whether a record carries its own reference bases.** Today it does, so every sample carries a copy
  of the reference over its footprint. The leaning is to drop it and re-fetch; nobody has timed the
  re-fetch. [`run_streaming.md`](run_streaming.md) §11 question 4.

---

## 12. Open questions

1. **Does a psp block ever need a second compressed stream?** — **NO, settled 2026-08-25** (§4.3).
   It was adopted on 2026-08-25 and reversed the same day when the walk was timed in three parts: a
   second stream adds its own pass on top of the walk over the records rather than replacing it,
   because a reader cannot seek inside a zstd frame. It is slower than a record head *and* costs a
   second decompressor — 2.27 GB against 1.14 GB at 5,000 samples. **The reversal is mine to own: I
   priced the split as though a separate stream of cheap fields removed the walk over the records, and it does
   not.**

2. **What byte ceiling, if any, should a writer put on a block?** — OPEN (§4.1), and it costs
   nothing measurable so far: a 100 kb grid with a 1 MiB ceiling gives 4.628 bytes a record against
   4.627 without, on tomato. *Leaning:* offer it, default it off, and let the first whole-genome
   deep-coverage run set it — at 279 reads a position a fully covered 100 kb block is about 1.6 MB,
   which is a large thing to hold while writing. **Settled by:** the block-size distribution on a
   whole-genome deep-coverage sample, which nothing here has produced. *Not* by seek time: the owner
   ruled on 2026-08-25 that reading a file start to end is the common case.
3. **How badly do near-empty blocks compress on a patchy sample?** — **ANSWERED 2026-08-25 (§4.1):
   a 5 kb grid costs about 10 % against a 1,000 kb one.** A sample with 74,623 covered positions
   scattered over small regions gives blocks of about 119 records and costs 18.242 bytes a record
   against 16.444 at 1,000 kb. **⚠ This entry used to conclude "so the rule that accumulates across
   empty spans ships", leaving only its threshold open. That is reversed: the owner ruled against
   merging on 2026-08-27** — it would break the cross-sample block alignment the grid exists for
   (§4.1). **The lever is the block size instead**, which is question 2's neighbour and live.
   *Still unmeasured:* a genuinely thin whole-genome sample — 1× rather than a region-restricted one
   — where the gaps are between positions rather than between regions.
4. **What does the reader's buffer pair cost in speed as it shrinks?** — **CLOSED 2026-08-26**
   (§4.4). Re-measured on a quiet machine with repeats: there is an optimum near 16 kB rather than a
   trend, and both larger and smaller are slower. 64 kB is 13 % slower and 256 kB is 40 % slower, so
   the earlier draft's puzzle — smaller coming out faster — was half of a curve, not a paradox. *Why
   it turns is still not established.*
