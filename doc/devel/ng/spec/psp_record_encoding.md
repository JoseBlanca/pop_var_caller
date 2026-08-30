# ng — how the psp stores one sample's observations

*Status: design spec, 2026-08-19; **revised 2026-08-25, when the shape was built and
measured.** This settles **what a record holds and in what form** — which fields, how each is
encoded, which may be stored approximately and at what step, and what a reader is required to hold
while reading one back. It does not fix a byte layout. **What that costs in disk and in memory
follows from those choices and from the settings a user picks; the figures throughout are measured
consequences, not quantities this document fixes.** Every one was taken on real files with the
programs named in §14.*

*What the revision changed, so a reader of the old version knows what moved: the block is
now **large** and the compressor's **look-back window** is capped separately, where before both were one
small block (§2, §2.2); the **dictionary is gone**, because a large block warms its own
context (§2.2); the **index question is closed**, because a large block has few entries
(§2.4, §13 Q4); and the byte-identity oracle is replaced by an exact-plus-tolerance one,
because block boundaries move QUAL in its last digits (§12). The measured result: **the cost
of an open sample falls 7.7-fold, 2,677 kB to 346 kB, while the file gets 35 % smaller, the
index 11× smaller and the read 1.8× faster** (§2.1).*

*Downstream: [`../arch/psp_file_format.md`](../arch/psp_file_format.md) is the code shape for this
and for the container, and [`../impl_plan/psp_file_format.md`](../impl_plan/psp_file_format.md) is the
build order.*

*This is the document [`run_streaming.md`](run_streaming.md) §10 defers to when it says
"the psp file's encoding — byte layout, compression, block sizing …, checksums, format
versioning, the index, the trailer. Its own spec beside this one". It inherits that
document's header fields (§6.1), its reader contract (§2.3), its per-open-file budget
(§7.2) and its worker-count-invariance restriction (§12.1). The chain ids are a field
big enough to have their own document —
[`psp_chain_id_encoding.md`](psp_chain_id_encoding.md) — and §6 here reports the
measurement that document's experiment was waiting for.*

---

## 1. What this is

A psp holds, for one sample, everything its reads showed at every position the run
analysed — one record per covered reference position, at three reads a position and at
three hundred, for a cohort of one sample and of several thousand. **This document settles what one
of those records holds and in what form**: which fields, how each is encoded, which may be stored
approximately and at what step, and what a reader must hold while reading one back.

**It does not settle what a record costs, and cannot.** The disk it takes and the memory a reader
spends are *consequences* of those choices together with settings that are deliberately the user's —
the block size and the approximation steps are declared per file, not fixed here
([`psp_file_format.md`](psp_file_format.md) §4). The same design gives 4.6 or 8.8 bytes a record
depending on how it is set. Every figure below is a measurement of one such setting, reported so the
choices can be made on numbers; none of them is a property the format guarantees.

**In production's design those two costs are one question**, because a block is both how far back the
compressor may look and how much a reader must inflate. Separating them is what the shape in §2 is
for.

### 1.1 The problem

Production's `.psp` groups records into blocks of a target uncompressed size
([`TARGET_BLOCK_BYTES`, `src/psp/writer.rs:72`](../../../../src/psp/writer.rs) — 1 MiB),
stores each block **columnar** — all the records' positions in one buffer, all their depths in
the next, and so on — and compresses each of those buffers as a separate zstd frame at level 9
([`ZSTD_COMPRESSION_LEVEL`, `src/psp/block.rs:709`](../../../../src/psp/block.rs)). The
block is therefore **both** the furthest back the compressor may look for a repeated
pattern **and** the amount a reader must decode before it can hand out its first record.
Shrinking it to save memory also costs bytes, which is the trade ng inherits if it ports
the shape — and ng cannot afford it, because
[`run_streaming.md`](run_streaming.md) §7.2 requires an open psp to cost *a few hundred
kilobytes, not megabytes*, at three thousand open samples.

### 1.2 Goals

1. **An open psp costs a few hundred kilobytes at most.** The owner set the working budget
   at **500 kB resident per open sample** on 2026-08-25 — 1.5 GB across three thousand — and
   [`run_streaming.md`](run_streaming.md) §7.2 was corrected to it on 2026-08-30. **Measured,
   the shape in §2 costs 346 kB** in the prototype, and **ng's own store costs 480 kB an open
   sample on a human reference**, of which 123 kB is the reader and the rest the header's copy
   of the reference's contig list ([`psp_file_format.md`](psp_file_format.md) §5.2). The goal is
   met rather than argued, with 4 % to spare.
2. **The file is smaller than production's at a comparable reader memory**, not only where
   production is allowed to spend more.
3. **A reader can start part-way through**, without reading what comes before, at the grain the
   **genomic block size** sets — default 100 kb
   ([`psp_file_format.md`](psp_file_format.md) §4.1).
4. **The cheap first pass survives.** [`run_streaming.md`](run_streaming.md) §3.3 makes
   psp mode's saving come from scanning one cheap number per position and building the
   full record only where some sample might vary — about one position in a hundred. An
   encoding that forces a reader to inflate every field to read one throws that away.
5. **Nothing about the writer's scheduling reaches the bytes.** Inherited from
   [`run_streaming.md`](run_streaming.md) §12.1: a psp block cut is a function of the
   observation stream alone, so the same sample gathered at any worker count gives the
   same file.

### 1.3 Non-goals, and what this document does not do

- **It does not change production's `.psp`.** ng replaces it; production stays as it is.
  Several findings below would apply there unchanged, and that is not a reason to touch it.
- **It does not fix a byte layout** — field order within a record, the exact framing
  integers, the trailer's bytes, the version tag. It fixes what is stored, in what form,
  in what unit, and what a reader must hold; the layout is the implementation's, guided
  by this.
- **It does not decide what the header contains.** [`run_streaming.md`](run_streaming.md)
  §6.1 already does, and its list is binding because the census names a psp by digesting
  exactly those fields
  ([`census_file.rs:89-98`](../../../../src/ng/parameter_estimation/joint/census_file.rs)).
- **It does not settle the chain ids' encoding.** That is
  [`psp_chain_id_encoding.md`](psp_chain_id_encoding.md)'s experiment; §6 here supplies
  numbers it was explicitly waiting for and states a leaning, nothing more.
- **It does not specify compression of anything but the observation stream** — the
  header stays plain text so that `head` and a TOML parser can read it, as production's
  does.

### 1.4 Vocabulary, defined once

- A **record** is one sample's observations at one covered reference position:
  [`SampleLocusObservations`](../../../../src/ng/locus_generation/mod.rs), which holds the
  region, the reference bases over it, a list of observed sequences with their support,
  and two counts of reads that produced no observation.
- A **chain id** identifies the DNA fragment a piece of evidence came from. **One `u64` per read
  *pair*, not per read** — mates are collapsed onto a single id — allocated in order from zero and
  never reused
  ([`chain_id_allocator.rs`](../../../../src/ng/locus_generation/pileup/chain_id_allocator.rs)).
  **The merge never uses its value; it only asks whether two records name the same fragment**
  ([`psp_chain_id_encoding.md`](psp_chain_id_encoding.md) §1.1), which is what leaves the stored form
  free to be anything that preserves equality. *An earlier draft of this document called these "read
  names", which is wrong twice: they are identifiers rather than names, and they are per pair rather
  than per read.*
- **psp block**, **zstd frame**, **zstd block** — three different things, and this document says
  which it means wherever either could be read. A psp block is a span of reference and the records
  in it; it becomes exactly one zstd frame; and a zstd block is zstd's own subdivision of a frame,
  at most 128 KiB, which nothing here chooses.
  [`psp_file_format.md`](psp_file_format.md) §2 has the table.
- **Columnar** is how production lays a block out, and the word is the format's own: a psp header
  carries a `[[column]]` array with one entry per field
  ([`src/psp/header.rs:24`](../../../../src/psp/header.rs)). Rather than writing record 1's fields,
  then record 2's fields, and so on, it gathers *all* the records' positions into one buffer — one
  **column** — all their depths into a second, all their allele sequences into a third, and
  compresses each column on its own.
- **Record-major** is the alternative this document proposes: each record's fields written together,
  then the next record's, with the bytes cut into blocks and each block compressed on its own. It is
  the layout a reader can walk without having to inflate anything it does not need.
- A **dictionary** is a few tens of kilobytes of representative bytes stored once in the
  file and handed to the compressor before every psp block, so that a small block is not
  compressed from a cold start. It is a standard zstd facility, not something we build.
- **Reader memory** is the quantity every choice here is measured against: the bytes one
  open sample forces a reader to hold before it can produce a record.

---

## 2. The proposal

*Revised 2026-08-25 after the shape below was built and measured. The earlier proposal —
small independent blocks compressed against a shared dictionary — is superseded, and §2.2
records why, because the reasoning that led to it was sound against a question that turned
out to be the wrong one.*

**Write the records record-major, not columnar — and make the psp block large while capping what a
reader has to hold.** Each record's fields go one after another; the bytes are cut into large **psp
blocks**; each psp block becomes **one zstd frame**, compressed with the **look-back window** capped
at write time. A reader never inflates a psp block: it pulls decompressed bytes into a small rolling
buffer, parses one record out of it, hands that record over and keeps nothing. Beside the records
runs a second, tiny stream carrying the one number the merge scans. A tail index names every psp
block.

```
psp file
  header                  plain-text, the fields run_streaming.md §6.1 fixes
  psp block 0             the records, one compressed stream, each record opening
                          with the head that lets a reader skip it (file format §4.3)
  psp block 1
  ...
  index                   one entry per block: chromosome, first position,
                          and the byte offset — nothing else (§2.4)
  trailer                 locates the index; its absence means "interrupted"
```

**The one idea the whole shape rests on: the block and the reader's memory are different
numbers.** In production's `.psp` they are the same one — a block is both how far back the
compressor may look for a repeat and how much must be inflated before the first record
exists — so every setting is a trade. Capping the window separates them: the block can be a
megabyte, for its ratio and for a small index, while a reader holds about a hundred kilobytes
— measured at 123 kB in ng ([`psp_file_format.md`](psp_file_format.md) §5.2).

**Two conditions have to hold together, and only the first is the compressor's doing:**

1. **Do not inflate the whole block.** The capped look-back window and an incremental decoder.
2. **Do not accumulate what you inflated.** The reader hands each record over and retains
   nothing. Satisfy only the first and the memory reappears in the caller's arrays — which
   is where the largest single mass of a cohort run's heap was found
   ([the memory review](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md) §2).

Everything below is a consequence of that shape or a number that justifies it.

### 2.1 Why record-major rather than columnar — priced

Measured like for like, same field encodings and same compressor, on a tomato accession at
about three reads a position:

| shape | reader memory | bytes a record |
|---|---|---|
| columnar | 1024 KiB | 6.27 |
| columnar | 256 KiB | 6.60 |
| columnar | 64 KiB | 7.33 |
| **record-major, 32 KiB blocks + dictionary** | **36 KiB** | **7.39** |
| record-major, 8 KiB blocks + dictionary | 12 KiB | 7.57 |

**Going record-major costs about 18 % of the bytes and buys about one
twenty-eighth of the memory.** That is the whole trade, and goal 1 decides it: three
thousand open samples at a megabyte each is 3 GB before a single record is read, at 36 KiB
each it is 108 MB.

The second property is what makes the small blocks work at all: **the record stream is
nearly flat in block size** — 132 KiB of reader memory down to 12 KiB costs 3 %, where the
columnar shape over the same range costs a third of its size. A design whose memory knob
is nearly free is a different kind of object from one whose memory knob is its ratio.

Against the file production writes today (11.85 bytes a record on that sample, 11.86 on
HG002), the record stream is 37 % and 43 % smaller — but most of that is §6's, not this
section's.

**Those rows all price memory against bytes, and the shape in §2 does not have to pay
that.** Built and walked against production's reader on the same 62-accession cohort, one
record at a time, the way a merge reads — 1 MB psp blocks, a 32 KiB look-back window, no dictionary:

| samples open | production `.psp` | this shape |
|---:|---:|---:|
| 1 | 5.0 MB | 2.5 MB |
| 8 | 24.3 MB | 4.8 MB |
| 31 | 85.8 MB | 12.7 MB |
| 62 | 164.5 MB | 23.1 MB |

| | fixed | **per open sample** | R² |
|---|---:|---:|---:|
| production `.psp` | 3.3 MB | **2.614 MB** | 0.9998 |
| this shape | 2.1 MB | **0.338 MB** | 1.0000 |

**The per-open-sample cost falls 7.7-fold, from 2,677 kB to 346 kB**, and both lines are
straight to four digits, which is what a per-open-sample cost looks like. Extrapolated over
the committed range — arithmetic from a fit over 1 to 62 samples, not a measurement — three
thousand samples goes from 7.66 GB to 0.99 GB.

**And nothing is traded away for it**, which is what makes this different from every row in
the table above:

| | production `.psp` | this shape | |
|---|---:|---:|---|
| bytes a record | 8.188 | 5.356 | 35 % smaller |
| cohort on disk | 3.52 GB | 2.38 GB | 32 % smaller |
| blocks per sample | 1,674 | 154 | index 11× smaller |
| 62-sample walk | 42.4 s | 23.1 s | 1.8× faster |
| records read | 471,520,156 | 471,520,156 | identical, same checksum |

*Measured on `examples/psp_row_stream_roundtrip.rs`, phases `encode-streaming`,
`verify-streaming` and `many-ngs`; tomato at about three reads a position. This is the
per-open-sample term only — a whole cohort run also carries the merger's projections and
the genotype fit, and neither is sized by the block.*

### 2.2 Why the dictionary was proposed, and why it is no longer needed — a superseded conclusion, kept

**This section used to argue for small psp blocks with a shared dictionary, against one long
stream with a capped look-back window. The measurement behind it was correct and the conclusion was
wrong, and the reason is worth keeping** — it is a clean example of sweeping a grid that
does not contain the answer.

What was measured: on HG002, a continuous stream with a **512 KiB window** gave 6.54 bytes a
record and 32 KiB blocks with a dictionary gave 6.52 — the same size for fourteen times the
memory, so the window looked like something you pay for and do not get back.

**What was never swept is the combination §2 proposes: a *large block* with a *small
window*.** Every streaming arm in that grid tied the window to the psp block, or used an 8 MiB
window for both. The two were separate knobs in the measuring program and were never varied
independently, so the one configuration that wins was not in the table. **The window is not
what costs memory when a reader decodes incrementally — the inflated block is** — and that
distinction is invisible in a sweep where the two always move together.

Two consequences:

- **The dictionary is no longer part of the design.** It existed because a 32 KiB block is
  compressed from a cold start; a 1 MB block warms its own context within the first few
  kilobytes. The measured store in §2.1 carries **no dictionary at all** and is still 35 %
  smaller a record than production's file.
- **The trap it created goes with it.** A dictionary is held per open reader unless
  deliberately shared, and 112 KiB across three thousand samples is about 330 MB — larger
  than the blocks the whole design existed to shrink. That risk is now simply absent rather
  than mitigated.

*Kept rather than deleted so that nobody re-derives the dictionary from the same table.*

### 2.3 The record head, so the cohort builds only the records that might vary

**Record-major has one real loss: a reader cannot read one field without decoding all of them.** The
cohort's first pass is exactly that shape — before anything can be called at a position, the caller
asks whether *any* sample shows something other than the reference there, and at about 99 positions
in 100 none does. Production measured that as 28,718 positions worth calling out of 2.83 million.

**So every record opens with a fixed head carrying what that first pass reads**, and the head ends
with the body's length so an unwanted record is skipped rather than decoded:

```
position_offset | reference_span | non_reference_reads | record_length | body
```

[`psp_file_format.md`](psp_file_format.md) §4.3 owns the layout and the field-by-field reasons. What
belongs here is the cost, because it is a cost this document's record pays.

**Three levels of work, timed on a tomato accession at three reads a position, 7.69 M records:**

| what a walk does | time |
|---|---:|
| decompress only, nothing else | 0.104 s |
| + walk each record's bytes to find where it ends | 0.163 s |
| + build the record objects — a full walk | 0.30 s |

**Decompression is 0.104 s, finding record ends is 0.059 s, building records is 0.137 s.** A reader
that wants one record in a hundred pays the first, skips the third, and the head's `record_length`
is what lets it skip most of the second: **about 0.126 s against 0.30 s.** The length costs a
measured 1.4 % of the file at three reads a position and 3.3 % at 279.

**Why not the columns production uses for this.** Its two-phase decode leaves the heavy columns
compressed while it reads the light ones
([`TwoPhaseSegment`, `src/var_calling/sample_reader.rs:698-712`](../../../../src/var_calling/sample_reader.rs),
then [`:789`](../../../../src/var_calling/sample_reader.rs)) — and *there* that saves memory, because
a block is inflated whole. **Here it would cost memory** and save nothing, because a reader never
inflates a block and cannot seek inside one: a separate stream of cheap fields adds its own pass on
top of the walk over the records rather than replacing it. §4.3 of the container spec has the
comparison; the short version is 2.27 GB against 1.14 GB at 5,000 samples, and slower.

### 2.4 The index, and where a reader may start

**Every psp block is self-contained.** It opens with its chromosome, its first position and its
record count, and every running difference inside it — the position difference, the coverage
difference, the chain-id difference — restarts at zero. A psp block never crosses a chromosome. So
the restart points are the psp block boundaries, and there is no separate seek mechanism to build.

**The index has one entry per psp block**, in genomic order: the chromosome, the first position, the
and its byte offset. **And nothing else**: an earlier draft had each entry carry the largest
non-reference support in the block, so a reader could skip a whole block without decompressing it.
At a 100 kb block with roughly one varying position in a hundred, every block holds about a thousand
of them, so that field never says *skip*. It cost index bytes and writer bookkeeping for a decision
that never fires, and it is gone (the owner, 2026-08-25).

#### How big the index is, at the settled block size

**A genomic grid makes the block count arithmetic rather than a measurement**: one block per 100 kb
of reference per chromosome. For tomato's 782 Mb across 13 chromosomes that is about **7,840
blocks**, and at production's 24 bytes an entry
([`BlockIndexEntry`, `src/psp/index.rs:42`](../../../../src/psp/index.rs)) an index of about
**190 kB per open sample** — under 40 kB if the entries are variable-length integers.

That is the whole of the problem this section used to have. **The earlier design cut blocks by byte
count at 32 KiB, which on a whole-genome sample gives roughly 500,000 blocks and a 5 MB index per
open file** — worse than the 3.8 MB [`run_streaming.md`](run_streaming.md) §7.2 rejects when it
measures production's flat vector, decoded whole at open
([`decode_index`, `src/psp/index.rs:110`](../../../../src/psp/index.rs)). Small blocks were forced
then because the block was also what a reader had to inflate; once the look-back window bounds that
instead (§2), the block is sized for compression and for the index, and both want it large.

**So the coarse-index-and-chain scheme this section once owed is not needed and must not be built.**
A flat vector of one entry per block, decoded at open, is a couple of hundred kilobytes. §13's open
question 4 is closed.

**And the extra cut rule it owed is gone too.** A byte-cut block can span an unbounded stretch of
reference on a sparse sample, so that design needed a second rule — *also cut every 100 kb* — to keep
goal 3's restart guarantee. **A block cut on a 100 kb genomic grid satisfies goal 3 by
construction**, so there is no second rule. What survives is a byte *ceiling*, and it is there for
the opposite case: at three hundred reads a position a fully covered 100 kb block is about 1.6 MB
([`psp_file_format.md`](psp_file_format.md) §4.1).

---

## 3. Settings, and what each is worth

Every number below was measured on both corners; the spread is what makes each of these a
knob rather than a constant.

**The psp block size and the look-back window are two settings, not one, and that is the point of
the design.** Confusing them is what made the earlier proposal (§2.2) reach for a dictionary.

| setting | proposed | what it moves |
|---|---|---|
| **block size** | **1 MB of records** | the compression ratio and the number of index entries. **Not the reader's memory.** |
| **look-back window** | **32 KiB** | **the reader's memory**, and almost nothing else once the psp block is large |
| rolling buffer | 64 KiB | how often the reader refills; one record must fit, and it grows if one does not |
| zstd level | 9 | 16 % of the file (tomato 8.15 / 7.44 / 6.86); write time 1.7 s / 4.0 s / 16.2 s |
| dictionary | **none** | superseded — a 1 MB block warms its own context (§2.2) |

**Level 9 is inherited, not derived** — it is what production uses
([`src/psp/block.rs:709`](../../../../src/psp/block.rs)) and the sweep gives no reason to
move: level 19 costs four times the write for 8 % of the file.

**Neither the psp block size nor the look-back window is a format change**, provided the window is written
into the file and the reader takes what it finds there rather than assuming a value. A
reader that assumes will either allocate more than it needs or refuse a legitimate file.

**What has not been swept: the block size itself.** 1 MB was chosen as comfortably large
and measured once; nothing here says whether 256 KB gives the same file and the same index
or whether 16 MB gives a better one. The window was capped at 32 KiB for the same reason.
Both are cheap to sweep now that a working writer and reader exist, and neither can move
the per-open-sample figure much — the window is what that is made of.

---

## 4. What is *not* stored

[`run_streaming.md`](run_streaming.md) §11 question 4 asks whether an observation's
reference bases should go in the file at all —
[`SampleLocusObservations::reference_bases`](../../../../src/ng/locus_generation/mod.rs) is
a `Box<[u8]>` per record, so written out it is a per-sample copy of the reference over the
analysed ground.

**Leaning: do not store them, and re-fetch when a psp block is decoded.** The measurement
supporting it is indirect but firm: the allele sequences, which are the same kind of
content, compressed to 0.32–0.34 bytes a record on both corners, and the reference bases
are longer than the allele on every record where a deletion widened the footprint. The
calling stage holds the reference already. **Confirm before code** by writing one sample
both ways and timing the re-fetch inside a segment walk; the question is
[`run_streaming.md`](run_streaming.md)'s, not this document's, and this is a leaning
offered to it rather than a decision taken from it.

### 4.1 The per-sample summary is not stored as text — settled 2026-08-25

**A psp carries one more thing besides records: a summary of the whole sample.** Today that is a
coverage-against-GC histogram and a count of heterozygous sites, and its consumers are the genotype
prior, production's hidden-paralog filter, and — when it is built — **ng's own duplication filter,
which does not exist yet but is coming** (the owner, 2026-08-25). The section is therefore permanent
and its encoding is worth fixing rather than inheriting.

**Production writes it as TOML — text — inside the metadata section, and that is what this decision
reverses.** A histogram is a matrix of counts; written as text each count becomes digits, separators
and indentation, and the parser must materialise the text before it can produce the matrix. Measured
on a 50-accession tomato run at about three reads a position: **1.05 MB of text per open sample
yielding a 0.52 MB parsed histogram** — and both are resident at once, because the reader keeps the
decompressed text alive for the life of the reader while the parsed structure sits beside it.

That cost is per open sample and nothing else moves it — not the block size, not the record layout,
not the chain ids, not the depth. At 50 samples the two together are 23.1% of the whole run's live
heap, against 12.7% for the compressed block decode this document spends most of its pages on.

**So, for ng:**

- **The sample summary is a binary section, not TOML.** Counts as variable-length integers, the bin
  schemes as a small fixed head. The plain-text header stays text, because `head` and a TOML parser
  reading a psp's identity is worth keeping (§1.3); the summary is not part of that and gains
  nothing from being readable in a pager.
- **The reader does not hold the encoded bytes once the summary is decoded.** Holding the source
  representation beside the decoded one doubles a per-open-sample cost for nobody. On production this
  was measured at 1.03 MB per sample of peak — 13.7% of the per-sample cost of a whole cohort run.
- **Both properties are per-open-sample costs and belong in the same budget as the blocks** (§2.1).
  A design that gets a reader down to 36 kB and then spends 1.5 MB on a summary has not made an open
  sample cheap.

**Open, and small:** whether the decoded summary is itself held per sample for the whole run or
folded into whatever the consumers need at startup. At 0.52 MB per sample it is 1.5 GB at three
thousand samples, which is inside the working budget but not free. Measured by whoever builds the
duplication filter, since the answer depends on what that filter reads and when.

*Numbers from [`../../reports/reviews/psp_memory_milestone_z_2026-08-25.md`](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md).*

---

## 5. The three quantities that arrive as floating point

**They are 70 % of the compressed file (§11) and none of them needs the precision it is
stored at.** Each is written instead as an integer count of steps, as a variable-length
integer, so a value needing one byte takes one and a value needing three takes three in the
same file. **There is no field width to choose** — the question that looks like "8 bits or
16?" is really "how big is the step?", and the compressor then charges for how much the
numbers move rather than for a declared width.

Swept one field at a time, at 32 KiB blocks and level 9, in bytes a record for the whole
file (so a row's distance from the baseline is that field's own contribution):

| quantity | step | HG002 | tomato |
|---|---|---|---|
| GC fraction of the window | 1/100,000 | 6.948 | 7.610 |
| | 1/10,000 | 6.727 | 7.443 |
| | 1/1,000 | 6.658 | 7.394 |
| | **1/100 — settled** | **6.026** | **6.792** |
| | 1/20 | 5.767 | 6.567 |
| mean coverage of the window | 1/256 read | 7.085 | 7.808 |
| | 1/16 read | 6.727 | 7.443 |
| | **1/4 read — settled** | **6.566** | **7.288** |
| | 1 read | 6.480 | 7.198 |
| summed log-error per allele | 1/4,096 ln | 7.047 | 7.774 |
| | 1/256 ln | 6.727 | 7.443 |
| | 1/16 ln | 6.415 | 7.091 |
| | 1/4 ln | 6.041 | 6.732 |

### 5.1 Re-measured 2026-08-25, on the streaming store and at 279 reads a position

**The sweep below was taken at 32 KiB blocks and on samples at about 3 and 30 reads a position.
This one is on the shape the container spec settled — 1 MB blocks, a 32 kB window — and it adds the
deep corner that was missing.** The deep sample is HG002 over chr21's tandem-repeat regions,
**279 reads a position** by its own coverage histogram, 74,623 covered positions.

Each field swept alone, the other two held at the prototype's settings, in bytes a record for the
whole file:

| quantity | step | tomato 3× | HG002 279× |
|---|---|---:|---:|
| GC fraction of the window | 1/100 | **4.659** | **16.104** |
| | 1/1,000 | 5.241 | 16.401 |
| | 1/10,000 | 5.356 | 16.630 |
| | 1/100,000 | 5.459 | 16.729 |
| mean coverage of the window | 1 read | 5.110 | 16.471 |
| | 1/4 read | **5.179** | **16.532** |
| | 1/16 read | 5.356 | 16.630 |
| | 1/64 read | 5.508 | 16.775 |
| summed log-error per allele | 1/64 ln | 5.296 | 16.383 |
| | 1/256 ln | 5.356 | 16.630 |
| | **1/1,024 ln** | **5.496** | **16.836** |
| | 1/4,096 ln | 5.793 | 17.273 |

**Three things this settles.**

**The GC fraction is the expensive one, and it was set far finer than anything reads it.** Taking it
from 1/10,000 to 1/100 is worth **13 % of the tomato file** and 3.2 % of the human one — more than
the other two together. Its consumer bins its input, so 1 % of GC and 0.01 % of GC are the same
number by the time anything uses it.

**A finer summed log-error is cheap, and cheaper at depth.** 1/256 → 1/1,024 costs **2.6 % of the
file at 3 reads a position and 1.2 % at 279** — a quarter of the error in the likelihood term
(0.4 % → 0.1 %) for one part in forty of the file. *This corrects an expectation recorded here
earlier:* the field's magnitude does grow with depth, so it costs more **bytes** at 279× (+0.206 a
record against +0.140), but the record it sits in grows faster, so its **share** falls.

**Combined, and the savings do not add — measured together:**

| GC · coverage · log-error | tomato 3× | HG002 279× |
|---|---:|---:|
| 1/10,000 · 1/16 · 1/256 (the prototype) | 5.356 | 16.630 |
| 1/10,000 · 1/16 · 1/1,024 | 5.496 (+2.6 %) | 16.836 (+1.2 %) |
| 1/100 · 1/4 · 1/1,024 | 4.629 (−13.6 %) | 16.255 (−2.3 %) |
| **1/100 · 1/4 · 1/4,096 — settled** | **4.907 (−8.4 %)** | **16.689 (+0.4 %)** |

**So the finer log-error is free and then some**: coarsening the two fields whose consumers cannot
tell the difference pays for it four times over, and the file still ends up 13.6 % smaller than the
prototype at low depth.

**Round-tripped at the settled steps on both samples** — 7,687,686 tomato records and 74,623 human
ones — with every integer field, allele sequence and chain-id list identical and each approximated
field inside half its own step.

*Two caveats. The human sample is 74,623 positions of tandem repeat, chosen because it is where the
300× reads are; it is small and it is not a random slice of a genome, so treat its per-record figures
as the right size rather than as precise. And a tomato row is not monotone at the finest end —
1/16,384 came back at 5.782 against 5.793 for 1/4,096 — which is unexplained and too small to chase.*

### 5.1.1 The bit-identity cost, and why the answer is to quantise *upstream*

**Approximating a float in the psp alone breaks the oracle the whole psp path is checked against.**
[`run_streaming.md`](run_streaming.md) §1.2 requirement 4: *"the same VCF from direct mode and from
psp mode with the parameters held fixed. This is the oracle for everything below."* Direct mode is
built first precisely to be that oracle. If the psp rounds a value that direct mode uses whole, the
two routes see different numbers and the check degrades from *identical* to *within a tolerance* —
which the block-boundary QUAL measurement in
[`psp_file_format.md`](psp_file_format.md) §10.1 shows is a much weaker test, one that can pass while
a chain-id list is being corrupted.

**So it matters what quantisation is worth.** Measured with the writer storing raw IEEE bytes as the
alternative — verified bit-exact, worst error 0.000000 on all three fields over 7,687,686 records —
starting from the GC fraction already being an integer:

| coverage · summed log-error | tomato 3× | HG002 279× |
|---|---:|---:|
| raw · raw | 8.842 | 22.043 |
| **1/4 read** · raw | 5.513 (−37.7 %) | 20.841 (−5.5 %) |
| raw · **1/1,024 ln** | 7.869 (−11.0 %) | 17.477 (−20.7 %) |
| **1/4 read · 1/1,024 ln** | **4.629 (−47.6 %)** | **16.255 (−26.3 %)** |

**Nearly half the file at three reads a position, and a quarter at 279.**

**Which of the two dominates flips with depth, and the reason is structural rather than
statistical.** The window's mean coverage is **one value per record**; the summed log-error is **one
value per allele** — `FrameWriter::push` writes it inside the per-allele loop. Depth multiplies the
second and leaves the first alone, so at 279 reads a position the log-error is the expensive one and
at three it is the coverage.

#### The summed log-error on its own — measured 2026-08-25

**The tables above move several fields at once, which makes the one field that carries a modelling
risk hard to read.** This one moves only it, in the settled configuration — record head, 100 kb psp
blocks, GC at 1/100, coverage at 1/4 of a read:

| the summed log-error stored as | tomato 3× | HG002 279× |
|---|---:|---:|
| a raw `f64` — bit-exact, no approximation | 6.021 | 23.440 |
| a count of 1/4,096 of a natural-log unit | 5.331 (−11.5 %) | 19.094 (−18.5 %) |
| **a count of 1/1,024** | **5.056 (−16.0 %)** | **18.570 (−20.8 %)** |
| a count of 1/256 | 4.927 (−18.2 %) | 18.270 (−22.1 %) |

**Two things to read off it, and the second is the useful one.**

**Storing this field as an integer at all is worth 16 % of the file at three reads a position and
21 % at 279.** That is the size of the decision.

**Which step barely matters.** From 1/256 to 1/1,024 — a quarter of the error — costs 2.6 % of the
file on tomato and 1.6 % on HG002, against the 16–21 % that separates any of them from a raw `f64`.
**So the precision is nearly free once the field is an integer**, and there is no reason to choose a
coarse step to save bytes. Even 1/4,096, a sixteenth of 1/256's error, still keeps three quarters of
the saving.

*What the error means: the value is a log, so rounding it by δ is a relative error of about δ in the
probability it stands for. 1/256 is 0.4 %, 1/1,024 is 0.1 %, 1/4,096 is 0.024 %. It is one rounding
of an already-summed quantity, not one per read, so it does not compound across the reads that went
into it.*

#### Settled: the rounding happens in the type, at ng level — the owner, 2026-08-25

**The three quantities are rounded where they are computed, not where they are written**, and the
step for the summed log-error is **1/4,096 of a natural-log unit**.

**This dissolves the bit-identity problem rather than trading against it.** If the rounding is
upstream of both routes, direct mode and psp mode compute the same rounded value; the file stores an
integer because the value *is* an integer by then; and there is nothing left to diverge. It is
better than storing a raw `f64` on the oracle's own terms, not merely equal to it: two routes that
sum in a different order differ in the last bits of an `f64`, and rounding to 1/4,096 **absorbs**
that difference and makes them agree where full precision would not.

**Why 1/4,096 rather than the 1/1,024 this document proposed.** Once the field is an integer at all,
the step is nearly free: 1/4,096 is a sixteenth of 1/256's error and still keeps three quarters of
the saving — 5.331 bytes a record against 5.056 at 1/1,024 and 6.021 raw. **The owner took the
precision** because the cost of taking it is 5 % of the file and the cost of being wrong about a
likelihood term is not measured in bytes.

**What this means for this document, precisely.** The psp does not approximate anything. It stores
an integer it was handed, and the header records the step so a reader can interpret it — **as a
property inherited from the type, not as a psp setting a writer chooses.** The step is therefore not
one of the knobs of §3, and a psp cannot be written with a different one than the types produce.

**And it is not this document's change to make.** Rounding at the type touches every ng module that
computes or carries these quantities, so it is ng-wide work. **It goes into the psp implementation
plan when that is written** (the owner, 2026-08-25) — recorded here so it reaches it, not owned
here. The two remaining
consequences for the encoding are the ones above: an integer field, and a header that records the
step it means.

*The window's GC fraction and its mean coverage take the same treatment, and were never in doubt:
they are terminal per-window statistics — computed once, read by a curve that bins them, never added
to again. The owner's own proposal for the GC fraction is a type holding it as an integer from 0 to
100.*

---

**A record at 279 reads a position costs 16.3 bytes against 4.6 at three reads** — three and a half
times, for a hundred times the depth. That is the chain ids, and it is §6's subject.

---

**Two of the three are free to coarsen and one is not.** The GC fraction feeds a
coverage-against-GC curve that bins its input, and the mean coverage feeds a ratio of
observed to expected depth; neither consumer can tell 1 % of GC from 0.01 %, or a quarter of
a read from a sixteenth. Taking those two to the proposed steps and leaving the log-error
alone gives **5.843 bytes a record on HG002 and 6.621 on tomato — 11–13 % below the
baseline, for no accuracy anything downstream consumes.**

**The summed log-error is the one that carried a modelling risk, and §5.1.1 records how it was
settled: 1/4,096 of a natural-log unit, rounded in the type rather than in the file.** It goes
straight into a likelihood: a step of 1/16 of a natural-log unit is a 6 % error in that
term, where 1/256 is 0.4 % and the settled 1/4,096 is 0.024 %. It is also the field whose magnitude grows with depth, which is
the second reason a fixed 8- or 16-bit width is the wrong shape for it — at three hundred
reads a position the value needs more range than sixteen bits hold, while at three reads it
needs six. **Proposed: 1/256 of a natural-log unit until the modelling side rules
otherwise** (§13, open question 1); coarsening it to 1/16 would buy a further 5 %.

**The three steps are written into the file.** A reader never has to be told them, so
changing one is a writer decision and not a format change. Encoding of a missing value: a
window that does not exist — an `N` reference position — is a real state and not a zero, so
the code 0 is reserved for it and every present value is shifted by one.

**Round-tripping this is lossy by construction and must be checked as such**, which is
§12's first oracle: every integer field, allele sequence and chain-id list comes back
identical, and these three come back inside their own step. Measured over 7.59 M tomato
records at the proposed steps: worst error 0.005 in GC fraction, 0.125 of a read in
coverage, 0.002 natural-log units in the summed log-error.

---

## 6. The chain ids, and why they decide the file's size at depth

**ng gives every read pair an identifier and records it at every position the pair covers.** This is the owner's ruling of
2026-08-17 and it is already in the code: the fast single-base path pushes a chain id for
every read with no reference test
([`fast_column.rs:315`](../../../../src/ng/locus_generation/pileup/fast_column.rs)), and the
general path's field carries the ruling in its own doc comment
([`open_record.rs:494`](../../../../src/ng/locus_generation/pileup/open_record.rs)). The
cohort merge needs it to answer *was this read here at all* — when a cohort locus spans
several of a sample's records, a read that agreed with the reference at one of them has to
be distinguishable from a read that never reached it.

Production names about 3.4 % of the reads it folds
([`src/pileup/walker/open_record.rs:150-160`](../../../../src/pileup/walker/open_record.rs)),
so **its files say nothing about what this costs.** Measured instead from the alignments
themselves — one identifier per read pair, allocated in order, the reference stretches taken
from each read's CIGAR, the resulting per-position live sets written three ways and
compressed identically at 32 KiB blocks with a block cut every 1,500 positions:

**The three ways, defined.** A position's chain ids are the identifiers of the read pairs covering
it, so at 300 reads a position that is about 300 identifiers, and the same pair reappears at every
one of the ~150 positions it covers.

- **The whole list, as raw identifiers.** At each position, the count and then each identifier as a
  fixed 8-byte integer. Identifiers run into the millions, so each costs its 8 bytes.
- **The whole list, as differences.** Sort the identifiers, store the first, and store each one after
  it as *the distance from the one before*. Because identifiers are allocated in order, pairs
  covering the same position have nearby ones, so those distances are small — often 1 or 2 — and a
  **variable-length integer** stores a small number in a single byte. The list is still written out at
  every position; each entry is just cheaper.
- **Only the changes.** Do not write the list at all. At each position write which pairs *started*
  covering it and which *stopped*, and let a reader carry the set forward. A pair covering 150
  positions then costs about two entries in total instead of appearing in 150 lists.

| | tomato slice, 11.4 reads a position | HG002 slice, 293 reads a position |
|---|---|---|
| whole list per position, raw 8-byte identifiers | 1.020 bytes a position | 43.78 |
| whole list, each identifier as its distance from the one before | 0.668 | 11.72 |
| **only the changes — who arrived, who left** | **0.432** | **6.42** |

*The deep corner is the HG002 benchmark slice, whose reads are concentrated into 1,000 small
regions; the depth is real and inside the committed range, but produced by that selection
rather than by a 300× library.*

Set against the roughly 5.4 bytes a record everything else costs, **the chain ids in
production's shape are 16 % of ng's file at eleven reads a position and 89 % of it at three
hundred.** As changes they are 7 % and 54 %. The reason to act is not the shallow corner: the
naïve form grows faster than depth — 25.7 times the depth cost 43 times the bytes — while
the differential form grows slower, 14.9 times.

Four things this settles that
[`psp_chain_id_encoding.md`](psp_chain_id_encoding.md) §10 listed as open:

- **The saving survives zstd** (its question 1). zstd is already very good at the repeated
  lists — 679 MB of raw identifiers on tomato became 7.4 MB, ninety-two fold — which is why
  the field looks cheap at low depth. At three hundred reads a position the same collapse
  only reaches thirty-six fold, and that is where the gap opens.
- **Delta-varints alone are worth having and may be enough** (its question 2). They capture
  **60 % of the available saving at eleven reads a position and 86 % at three hundred**, with
  no reader state, no residual arithmetic and no new error class.
- **An identifier goes live more than once for most reads** (its question 4), and not
  marginally: **83 % of identifiers on HG002 and 91 % on tomato** cover two stretches with a
  gap between them, because a pair's mates rarely overlap. An arrivals-and-departures stream
  that assumes one stretch per read loses the second mate of nine reads in ten — silently,
  because the merge would simply see a read that was not there. **A re-entry form is part of
  the first version, not a later fix.**
- **Restating the live set at every block is affordable.** Cutting blocks every 1,500
  positions rather than by byte count — which is what §2.4 does — costs the differential form
  12 % of its own bytes on tomato (0.385 → 0.432) and leaves it far ahead of both
  alternatives.

**Settled 2026-08-25 by the owner: changes-only ships, and distances are not built.** The
intermediate form was proposed here as a first version because it is a few lines and has no
interactions; the owner has taken the destination directly.

#### What that means for the record head, because the two do collide

**Changes-only works by carrying a set forward**: a reader knows which read pairs are live only
because it applied every arrival and departure since the block began. **A reader that skips a
record's body would never see that record's changes, so its set would go stale and every later
record it wanted would be wrong.** Restating the live set at every record is not a fix — that *is*
writing the whole list at every position, which is the form we are leaving.

**The resolution falls out of the shape the encoding already has.**
[`psp_chain_id_encoding.md`](psp_chain_id_encoding.md) §4 splits the chain ids into two parts, and
they sit on opposite sides of the skip:

- **the live-set changes — who arrived, who left — go in the head.** They carry the state, so every
  reader decodes them whether or not it wants the record.
- **the exception lists — the ids of every observation except the residual one — stay in the
  skippable body.** They carry no state and are only needed by a reader that is building the record.
  They are the ~3.4 % of ids production stored before the 2026-08-17 ruling, so this half is small.

A reader that wants a record then has everything: the live set from the heads it has been reading
all along, and the exceptions from the body it just decoded. **The residual observation's ids are
the live set minus the others**, which is what makes the whole scheme cheap and is also where it
fails silently (§7).

**What this costs, and it moves the wrong way with depth.** The live-set changes are 0.432 bytes a
position at 11.4 reads and **6.42 at 293**. So at low depth the head stays small and a skipping
reader avoids most of the record; **at high depth the head carries most of the bytes and the skip
saves much less.** The chain-id saving and the skip saving therefore pull in opposite directions
across the committed range, and **how much of the 2.06× survives at depth is unmeasured** — it is
the first thing to measure once a writer exists.

**⚠ A stale paragraph stood here and has been cut, 2026-08-28.** It repeated the four settled
questions above verbatim and then argued *"Distances ship. Changes-only is now in doubt"* — that
the two forms are alternatives at depth, and that putting the changes in the head and the rest in
the skippable body was *"a third shape nobody has priced"*. That third shape **is** what §13
question 2 settles, what this section's own resolution above describes, and what the architecture's
decisions list records; it was built as Milestone E of
[`../impl_plan/psp_file_format.md`](../impl_plan/psp_file_format.md) on 2026-08-27. The paragraph
was the only text in three documents saying otherwise, and it survived long enough to be read as an
open question by someone implementing from it.

**The one thing it was right about is the paragraph above**, which stands: how much of §2.3's
2.06× survives at depth is unmeasured, and **nothing can measure it until Milestone F opens a
file**. Every figure in this section is from a prototype over alignments, not from the writer that
now exists.

---

## 7. Traps — what will bite the coder

- **[`run_streaming.md`](run_streaming.md) §8 says chain ids are already omitted for reads
  that agree with the reference, "recorded so nobody re-adds them when the psp writer is
  built". That bullet is stale and following it produces a file the cohort merge cannot
  use.** The code names every read
  ([`fast_column.rs:315`](../../../../src/ng/locus_generation/pileup/fast_column.rs),
  [`open_record.rs:494`](../../../../src/ng/locus_generation/pileup/open_record.rs)), by the
  2026-08-17 ruling that post-dates the bullet. Fix the bullet when that document is next
  touched.
- **⚠ Streaming the block is half the job; not accumulating the records is the other half,
  and it is the half that gets lost.** A reader can decode incrementally and still gather
  every record into per-sample arrays before handing anything on, at which point the memory
  is exactly where it was. In production's cohort run those assembled per-sample columns are
  the largest single mass of the heap, larger than the decompression buffer they came from
  ([the memory review](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md) §2). The
  reader must hand each record over and retain nothing.
- **A record can straddle the rolling buffer, and the parse has to be restartable.** The
  buffer holds whatever has been decompressed so far, which may be the first half of a
  record. Running out of bytes has to be an answer the parser can give — not a panic and not
  a short read — and the retry must resume from the record's *start* with the running
  position, coverage and chain-id bases restored to what they were. A parse that half-
  advances that state before failing corrupts every record after it, plausibly.
- **A single record can exceed the rolling buffer.** Many alleles, many chain ids. The
  buffer has to grow rather than fail, and a fixed maximum record size is not a safe
  assumption to bake in.
- **The look-back window must be written into the file and honoured on read.** zstd will refuse a zstd frame
  whose window exceeds the decoder's configured maximum, so a reader that assumes 32 KiB will
  reject a legitimate file written with more. Read it from the header and configure the
  decoder from that.
- **A dictionary is no longer part of this design (§2.2)**, and if one is ever reintroduced
  it is held *per open reader* unless deliberately shared: 112 KiB across three thousand
  samples is about 330 MB, larger than the buffers the whole design exists to shrink.
- **A dictionary trained on the blocks it is then measured against reports a saving no
  reader ever gets.** This is easy to do by accident and the number it produces is
  spectacular — in an early run of our own probe, a hundredfold. Train on one half of the
  blocks and measure on the other.
- **A psp block's span in reference is fixed by the grid, but its size in bytes is not.**
  §2.4's 100 kb cut rule exists for that; without it goal 3 quietly fails on exactly the
  sparse samples that need it.
- **The three quantised fields are lossy and the integer fields are not.** A round-trip test
  that compares whole records with a tolerance will pass while the chain-id list is being
  corrupted. Compare the integer fields, the sequences and the name lists **exactly**, and
  only the three quantised fields against their step.
- **The window's mean coverage is stored as a difference from the previous record**, so a
  psp block that does not reset that difference reads back wrong from its first record — and
  plausibly, because coverage is smooth. The same applies to the position and the chain-id
  bases. Every running base resets at a psp block boundary; this is the property §2.4's restart
  guarantee rests on.
- **Do not decode the whole index at open** — §2.4, and the 3.8 MB per file
  [`run_streaming.md`](run_streaming.md) §7.2 measures on production's.

---

## 8. Cross-cutting concerns

**Memory.** The point of the design. One open sample holds, per stream: a buffer of compressed
bytes, a rolling buffer of decompressed ones, the decompressor's own state, and one record. The
whole-cohort figure was measured — 63 tomato accessions open at once and advanced one record
each per round, the shape the merge reads in: **22 MB peak against 170 MB for the same walk
over the `.psp`s**, and the psp side is flattered there because those particular files
average 258 KiB of block rather than the writer's 1 MiB default.

**Speed.** Decoding was measured at 20.7 M records a second against 12.6 M for the `.psp`
reader on the same sample, and the cohort walk finished in 22 s against 38 s. The stream
never builds a field array it must then walk again to assemble records. This is a
single-threaded sequential-read measurement and says nothing about the seek path.

**Errors.** A zstd frame that decompresses to the wrong length, a record that runs off the end of
its psp block, a chain-id list longer than the observation's read count: all are corrupt-input
failures belonging to the psp reader's error type, and none may reach the merge as a
half-built record. The trailer's absence means the writer was interrupted and the file is
refused rather than read as a short sample —
[`run_streaming.md`](run_streaming.md) §9's rule, unchanged.

**Concurrency.** psp blocks are independently decodable, so nothing here serialises. The block
cut must depend only on the observation stream, never on scheduling, or the byte-identity
oracle of [`run_streaming.md`](run_streaming.md) §12.1 breaks — a byte target counted over
the records as they arrive satisfies this; a flush driven by a timer or a queue depth does
not.

---

## 9. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| variable-length integer codec | `src/psp/varint.rs` | as-is: LEB128 and zig-zag LEB128, already specified and tested |
| plain-text header framing | `src/psp/header.rs` | the pattern — magic, length prefix, TOML, sentinel — so `head` still works on an ng psp |
| the zstd seam | `new_column_compressor`, `zstd_compress_into` ([`src/psp/block.rs:718,730`](../../../../src/psp/block.rs)) | the shape: one long-lived compressor per writer, frame checksums on. The dictionary is new |
| the two-phase read | `TwoPhaseSegment`, `set_variable_rows` ([`src/var_calling/sample_reader.rs:698,789`](../../../../src/var_calling/sample_reader.rs)) | **the idea, not the shape.** Its light/heavy split becomes a head on each record rather than separate columns, for the reason in §2.3 |
| the eager whole-segment decode | [`sample_reader.rs:20-26`](../../../../src/var_calling/sample_reader.rs) | the parity oracle's model: a simple decoder used only by tests, against which the real one is byte-compared |
| the record | `SampleLocusObservations` ([`src/ng/locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | what is written and what must come back |

**Parity oracle:** the same shape production uses — a deliberately simple whole-file decoder
used only by tests, and the real reader compared against it record for record.

---

## 10. Deferred, with a recommended home

- **The byte layout itself** — field order inside a record, the framing integers' widths, the
  trailer's bytes, the format version tag. To the implementation, guided by this document;
  it is small enough not to need its own spec once the shape is fixed.
- **The chain ids' final encoding** — [`psp_chain_id_encoding.md`](psp_chain_id_encoding.md),
  whose experiment §6 here feeds rather than replaces.
- **Whether the reference bases are stored** — [`run_streaming.md`](run_streaming.md) §11
  question 4; §5 offers a leaning and the measurement that would close it.
- **Multi-library samples.** Both corners measured here are one read group. A sample with
  several will have more distinct values per psp block and compress somewhat worse; nobody has
  measured how much. To the first run on a multi-library cohort.

---

## 11. What a record costs today — the measurement the design was chosen against

**This is where §2's shape came from.** Production's file, broken down by field, per record,
compressed, at 512 KiB blocks — GIAB HG002 at about thirty reads a position:

| field | today, as written | as a fixed-point integer |
|---|---|---|
| mean coverage of the window | 2.762 bytes | 0.120 |
| summed log-error per allele | 2.817 | 1.967 |
| GC fraction of the window | 0.806 | 0.597 |
| *the eleven other fields together* | 2.795 | 2.754 |
| **total** | **9.180** | **5.438** |

A tomato accession at about three reads a position gives the same shape: 2.788 + 2.905 +
0.613 of 9.909 bytes a record, falling to 0.121 + 2.067 + 0.512 of 6.269.

**Three floating-point fields are 70 % of the compressed file.** They are stored at full
IEEE precision — two `f32`s and an `f64`
([`src/psp/registry.rs:366,382,406`](../../../../src/psp/registry.rs)) — and the bottom
bits of a window's mean coverage or a sum of log error probabilities are arithmetic noise,
which is the one thing a compressor cannot shrink. §6 is about that, and it is worth more
than every layout decision in this document put together.

**One field measured here says nothing about ng, and it is the one that matters most at
depth.** Production stores the names of about 3.4 % of the reads it folds
([`src/pileup/walker/open_record.rs:150-160`](../../../../src/pileup/walker/open_record.rs)),
ng stores all of them, and §7 measures that field separately.

---

## 12. How we know it works

1. **Round-trip, with the right strictness per field.** Every integer field, every allele
   sequence and every chain-id list identical; the three quantised fields inside their own
   step. Already demonstrated on 7.59 M tomato records and 0.60 M HG002 records with the
   probe of §14.
2. **The per-open-file budget, measured rather than argued.** N samples open and walked in
   lockstep, peak resident reported, against the 500 kB budget of
   [`run_streaming.md`](run_streaming.md) §7.2. The probe does this at 8, 32 and 63 samples
   today; ng's own store was measured this way at 1 to 5,000 samples on 2026-08-30
   ([`psp_file_format.md`](psp_file_format.md) §5.2).
3. **Worker-count invariance**, inherited from [`run_streaming.md`](run_streaming.md) §12.1:
   one sample gathered at 1, 2, 4, 8, 16 workers gives byte-identical files apart from the
   header's timestamp. This is what §8's block-cut rule exists to preserve.
4. **The cheap first pass is still cheap.** Walking a segment reading only record heads must
   build about one record in a hundred, the ratio
   [`run_streaming.md`](run_streaming.md) §3.3 measured. A file that round-trips but forces
   every record to be built has failed goal 4 while passing oracle 1.
5. **Restart equals sequential.** Reading from an arbitrary psp block gives exactly the records a
   full sequential read gives from that point — the test that catches a running base that was
   not reset (§7).
6. **Mode equivalence**, [`run_streaming.md`](run_streaming.md) §12 oracle 3: the same cohort
   called through the direct route and through the psp route gives the same VCF. It is the
   sufficiency test for everything this document chose not to store.

**⚠ Oracle 6 cannot be a byte-identical VCF, and this was measured rather than feared.**
Changing only how records are grouped into blocks changes the order in which per-sample
evidence is summed, and floating-point addition is not associative. On the 50-accession
tomato cohort, rewriting the same records at a 20 kb and an 80 kb block window left the set
of sites and alleles **identical to the line** and returned **1,194 records of 180,366 with a
different QUAL** — median difference 0.010 Phred, largest 4.48. Nothing else moved: no
genotype, no GQ, no DP, no AD, no INFO field differed anywhere in the cohort.

**One site in 180,366 crossed the emission gate as a result** — `SL4.0ch10:58265030 T>A`
sits at QUAL 30.075 against a minimum of 30, and it is emitted at 20 kb and 80 kb and not at
5 kb. Seventeen differing records sit below QUAL 40, so within reach of the gate.

**So the oracle is: the site list, genotypes, GQ, DP and AD exactly; QUAL within a
tolerance; and an explicit statement about the gate.** A test that demands byte-identity
will fail on a correct implementation, and a test that compares everything with a tolerance
will pass while a chain-id list is being corrupted.

---

## 13. Open questions

1. **What step may the summed log-error be stored at?** — **SETTLED 2026-08-25: 1/4,096 of a
   natural-log unit, rounded in the type rather than in the file** (§5.1.1). A 0.024 % error in the
   term, against 0.4 % at 1/256. The step was taken finer than this document proposed because it is
   nearly free once the field is an integer at all: 1/4,096 costs 5 % of the file against 1/1,024 and
   still keeps three quarters of the 16–21 % that separates any integer step from a raw `f64`.

   **Two things follow and neither is this document's to build.** The rounding happens where the
   value is computed, so it is ng-wide work in whichever implementation plan owns those types. And
   because it happens upstream of both routes, direct mode and psp mode see the same value — which is
   the point: approximating in the psp alone would have broken the oracle
   ([`run_streaming.md`](run_streaming.md) §1.2) that the whole psp path is checked against.
2. **Does the changes-only chain-id encoding ship, or only
   distances only?** — **SETTLED 2026-08-25: changes-only, and distances are not built**, and
   **BUILT 2026-08-27** as Milestone E of
   [`../impl_plan/psp_file_format.md`](../impl_plan/psp_file_format.md). The two do collide with the
   skippable records of §2.3, because a set carried forward cannot survive a skipped record — and
   §6 resolves it by putting the live-set changes in the head and leaving the exception lists in
   the body, which is what was built. **What is not settled is what that costs at depth**: the
   changes are 0.432 bytes a position at 11.4 reads and 6.42 at 293, so the head grows with depth
   and the skip's value shrinks. **Settled by:** measuring the skipping walk on a head that carries
   the changes — which needs a file, so Milestone F.
3. **Are a block's parts interleaved or separated into two regions of the file?** — **MOOT**: a psp
   block is one stream (§2.3), so there is nothing to interleave. *Kept because the reasoning applies
   again if the head's live-set changes are ever split off into a part of their own: interleaving
   keeps a segment's parts adjacent, so serving one segment is one seek; separating makes a whole-file
   scan of one part
   one sequential read. *Leaning:* interleave, because
   [`run_streaming.md`](run_streaming.md) §3.4's reader serves segments rather than files.
   **Settled by:** timing a cohort merge over one chromosome both ways, once a writer exists.
4. **How does the index stay small when a file has 100,000 blocks?** — **CLOSED 2026-08-25:
   the file does not have 100,000 blocks.** The question only existed because a psp block had to
   be small to bound a reader's memory; once the look-back window does that instead, the psp block is sized
   for compression and for the index, and both want it large. Measured: 154 blocks in a tomato
   accession at 1 MB blocks against 1,674 in production's `.psp`, so a flat index of one entry
   per block is a few hundred kilobytes on a whole-genome sample (§2.4). **The coarse-index-
   and-chain scheme must not be built.**
5. **How much do the numbers move on a multi-library sample?** — OPEN. Both corners here are
   one read group, and the read group joins an observation's identity
   ([`src/ng/locus_generation/mod.rs:245`](../../../../src/ng/locus_generation/mod.rs)), so a
   multi-library sample has more, smaller observations per record. *Leaning:* somewhat worse,
   not structurally different. **Settled by:** running the probe on a multi-library sample.

---

## 14. Where the numbers came from

Two throwaway programs, both in the tree, both runnable:

- **`examples/psp_record_stream_compression.rs`** — sweeps the grid of §2.1 and §2.2: field
  layout (record-major or columnar), field width (fixed-width, variable-length, fixed-point),
  batch size, and framing (independent blocks, with or without a dictionary; one continuous
  stream with a chosen window, flushed or not). Reports bytes a record against the reader
  memory each combination implies, and each field's own compressed bytes.
- **`examples/psp_row_stream_roundtrip.rs`** — a working encoder, decoder and verifier for the
  proposed shape, plus the many-sample walk of §8. `verify` walks the new store and the `.psp`
  in lockstep and fails on the first record that disagrees.
- **`examples/ng_chain_id_column_cost.rs`** — §6's measurement, taken from a CRAM: read pairs
  given one identifier each, their reference stretches from the CIGAR, the per-position live
  sets written three ways.

The full report, with the sweeps these tables summarise, is
[`../research/per_sample_record_store_compression_2026-08-19.md`](../research/per_sample_record_store_compression_2026-08-19.md).

**What was measured on, and what that covers.** A tomato accession at about three reads a
position (`SRR7279540`, 7.59 M records) and GIAB HG002 at about thirty (0.60 M records) for
everything but §6; for §6, a tomato benchmark slice at 11.4 reads a position and an HG002
benchmark slice at 293. Cohort figures are the 63-accession tomato panel. **One sample and one
cohort of 63 is not the committed range** — nothing here was measured at one sample of a
thousand-sample cohort, and the memory claims at three thousand samples are arithmetic from the
per-sample figure, not measurements.
