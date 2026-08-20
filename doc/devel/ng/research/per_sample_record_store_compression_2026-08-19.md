# A per-sample store that compresses a stream of records, not a grid of columns

**Status:** research complete, no production code changed. **Date:** 2026-08-19.
**Measured on:** two real per-sample pileup files written by the shipped
`pileup` subcommand — one tomato accession at about three reads a position
(`SRR7279540`, 7.59 M records) and the GIAB human sample HG002 at about thirty
(`HG002.30x`, 0.60 M records) — plus all 63 tomato accessions together for the
memory measurement.
**Programs:** `examples/psp_record_stream_compression.rs` (the layout sweep) and
`examples/psp_row_stream_roundtrip.rs` (a working encoder, decoder and verifier).

---

## The answer, in one paragraph

**A per-sample file that just serialises records one after another and cuts them
into small independently-compressed frames is not a compromise: measured against
today's files it is 37 % smaller on tomato and 43 % smaller on HG002, decodes 1.6
to 1.8 times faster, and needs 32 KiB of reader memory instead of a whole block.
With all 63 tomato accessions open at once and walked in lockstep, peak memory
falls from 170 MB to 22 MB and the walk finishes in 22 s instead of 38 s.** But
the credit does not belong where it looks like it belongs. Two thirds of the size
win comes from a change that has nothing to do with layout — three floating-point
fields are stored at full IEEE precision and are 70 % of the compressed file
between them, and replacing them with fixed-point integers costs nothing anyone
downstream can measure. Dropping the columns *by itself* costs about **18 % more
bytes**, in exchange for holding about **1/28 of the memory**. Both changes are
available independently, and the float one is available to the current format
without touching its design at all.

---

## Vocabulary, defined once

- A **record** is what one covered reference position of one sample holds: the
  position, its alleles, and per allele five support counts, a summed log-error,
  two mapping-quality moments and the ids of the reads folded into it.
- The **store** is the per-sample file these records are written to — today the
  `.psp`, in ng whatever replaces it.
- **Transposed** (the current design): the store cuts records into **blocks** of
  a target uncompressed size, and inside a block every field becomes its own
  buffer — all the positions together, then all the allele counts together — each
  compressed as a separate zstd frame.
- **A record stream** (the alternative): each record's fields are written one
  after another and the bytes are cut into **frames** of a target uncompressed
  size, each compressed on its own.
- A **dictionary** is a few tens of kilobytes of representative bytes stored once
  in the file and fed to the compressor before every frame, so a small frame is
  not compressed from a cold start. This is a standard zstd feature.
- **Reader memory** is the number this whole report is against: the bytes one
  open sample forces a reader to hold before it can hand out its first record.
  For the transposed store that is the block, twice over — the decompressed bytes
  and the decoded field arrays. For a record stream it is one frame.

---

## 1. Why reader memory is the axis, and not file size

The two designs are not on a size ladder with the transposed one above. They are
both trading the *same* resource, and the current design ties two uses of it
together that need not be tied.

In the transposed store the block is simultaneously (a) the furthest back the
compressor may look for a repeated byte pattern and (b) the amount a reader must
decode before the first record exists. Making it smaller to save memory therefore
also costs ratio, which is exactly the trade the format spec records: a sweep on
tomato with 18 samples on 4 threads moved peak memory from 2.5 GB at 16 MiB blocks
down to 59 MB at 64 KiB blocks, and paid for it in bytes.

A record stream unties them. A frame can be 32 KiB — the reader's whole cost —
while a stored dictionary supplies the compressor with the context that a 32 KiB
window would otherwise not have. So the honest comparison is not "which is
smaller" but **"which is smaller at the same reader memory"**, and that is what
the tables below plot.

## 2. Where the bytes actually are today

Before any layout question: three of the fourteen fields are floating-point, and
they are most of the file. Per record, compressed, on HG002 with 512 KiB blocks:

| field | today (IEEE float) | as fixed-point integers |
|---|---|---|
| mean coverage of the window | 2.762 bytes | 0.120 |
| summed log-error per allele | 2.817 | 1.967 |
| GC fraction of the window | 0.806 | 0.597 |
| *all eleven other fields together* | 2.795 | 2.754 |
| **total** | **9.180** | **5.438** |

Tomato is the same shape: 2.788 + 2.905 + 0.613 of 9.909 bytes a record, falling
to 0.121 + 2.067 + 0.512 of 6.269.

The reason is mechanical. A window's mean coverage and a sum of log error
probabilities are smooth physical quantities whose bottom mantissa bits are
arithmetic noise, and noise is the one thing a compressor cannot shrink. Storing
the GC fraction to one part in 10,000, the mean coverage to 1/16 of a read as a
difference from the position before, and the summed log-error to 1/256 of a
natural-log unit **cuts the file by 37 % on tomato and 41 % on HG002 with no
layout change whatever**. The worst error any record suffered, checked over all
7.59 M tomato records: 0.00005 in GC fraction, 1/32 of a read in coverage, and
0.002 natural-log units in the summed log-error. Coarser quantisation keeps
paying — 1/4 of a log unit instead of 1/256 takes the log-error field from 1.967
to 1.507 bytes — so there is a knob here, not a cliff.

**This is worth stating separately because it is independent of every other
decision here** — it is a change to how three numbers are written down, not to
how the file is organised, and either layout gets all of it. (It is *not*
proposed for production's `.psp`, which stays as it is: ng replaces it.)

## 3. Size against reader memory

Bytes per record, at zstd level 9. `reader KiB` is the memory one open sample
forces. The first block of rows is today's design; the second is today's design
with the two encoding fixes (varints and fixed-point floats) and nothing else;
the third is the record stream.

**Tomato accession SRR7279540, ≈3 reads a position, 4.0 M records
(today's file: 11.85 bytes a record):**

| design | reader KiB | bytes/record |
|---|---|---|
| transposed, fixed-width fields (today) | 4096 | 10.24 |
| " | 1024 | 11.10 |
| " | 256 | 12.71 |
| " | 64 | 15.41 |
| transposed + varints + fixed-point floats | 1024 | 6.27 |
| " | 256 | 6.60 |
| " | 64 | 7.33 |
| **record stream + dictionary** | **132** | **7.32** |
| **record stream + dictionary** | **36** | **7.39** |
| record stream + dictionary | 12 | 7.57 |
| record stream, no dictionary | 36 | 7.87 |

**GIAB HG002, ≈30 reads a position, 0.60 M records
(today's file: 11.86 bytes a record):**

| design | reader KiB | bytes/record |
|---|---|---|
| transposed, fixed-width fields (today) | 4096 | 9.28 |
| " | 1024 | 10.17 |
| " | 256 | 11.81 |
| " | 64 | 14.41 |
| transposed + varints + fixed-point floats | 1024 | 5.44 |
| " | 256 | 5.77 |
| " | 64 | 6.47 |
| **record stream + dictionary** | **132** | **6.48** |
| **record stream + dictionary** | **36** | **6.52** |
| record stream + dictionary | 12 | 6.70 |
| record stream, no dictionary | 36 | 7.12 |

Three things to read off these.

**The record stream is nearly flat in frame size.** Going from 132 KiB of reader
memory down to 12 KiB costs 3 % on both samples. That flatness is what the
dictionary buys, and it is the property that makes the design work: memory can be
dialled to almost nothing without watching the file grow. The transposed store is
the opposite — with the identical field encoding, moving it from 256 KiB of reader
memory to 16 KiB takes tomato from 6.60 to 8.80 bytes a record, a third larger.

**Dropping the columns costs 18 %, not more.** Compare like with like at the same
encoding: transposed at 1024 KiB gives 6.27 bytes a record on tomato, the record
stream at 36 KiB gives 7.39. That is the whole price of the layout change, and it
buys a 28-fold cut in reader memory.

**Against today's file, the record stream wins outright.** 7.39 against 11.10 on
tomato, 6.52 against 10.17 on HG002, at 1/28 of the memory.

## 4. A store that actually works, not a projection

The numbers above come from an encoder that throws its output away, which proves
nothing about whether the bytes can be read back. So the second program writes a
real file — records serialised end to end, cut into 32 KiB frames, each
compressed against a dictionary trained on that file's own records and stored in
its header — and reads it back.

`verify` walks the new store and the `.psp` in lockstep and fails on the first
record that disagrees. Over all 7.59 M tomato records and all 0.60 M HG002
records, every integer field, every allele sequence and every read-id list came
back identical, and the three quantised fields came back inside the tolerance
their quantisation implies.

| | tomato SRR7279540 | GIAB HG002 |
|---|---|---|
| today's `.psp`, bytes a record | 11.85 | 11.86 |
| record store, bytes a record (dictionary included) | **7.44** | **6.73** |
| decode, records a second | 20.7 M | 21.8 M |
| today's `.psp`, records a second | 12.6 M | 11.8 M |
| decode peak memory, one open sample | 2.4 MB | 2.4 MB |
| today's `.psp`, same | 5.3 MB | 10.2 MB |

The store is smaller **and** faster to read. The speed is not a surprise in
hindsight: the record stream never builds a field array it then walks a second
time to assemble records, and reading a varint is cheaper than the transposed
path's slab bookkeeping.

## 5. The cohort case, which is what the memory question was about

All 63 tomato accessions opened at once and advanced one record each per round —
the shape a cohort merge reads in — on one thread:

| samples open | `.psp` peak memory | record store | `.psp` seconds | record store |
|---|---|---|---|---|
| 8 | 24.9 MB | 4.6 MB | 4.30 | 2.75 |
| 32 | 90.5 MB | 12.3 MB | 18.61 | 11.36 |
| 63 | 170.4 MB | 22.0 MB | 38.25 | 22.37 |

Both walks yielded the same 476,758,046 records. On disk the cohort is 3,665 MB
of `.psp` against 2,260 MB of record store.

Two caveats on the memory column, in opposite directions. **In the `.psp`'s
favour:** these particular files average about 258 KiB of uncompressed block, not
the 1 MiB the writer defaults to, so a cohort written with the default would sit
roughly four times higher than the 170 MB shown. **Against the record store:**
its 312 KiB per open sample is mostly the 112 KiB dictionary, held once per
reader; zstd lets one prepared dictionary be shared by every reader, and a 16 KiB
dictionary costs only 3 % more bytes, so this figure has a lot of slack in it.

## 6. What a record stream gives up, and what it costs to keep

One thing, and it matters here more than in most projects: **a reader can no
longer read one field without decoding all of them.**

The cohort merge depends on that. The locus-stream-shape experiment of 2026-08-06
found that when the merge reads a sample it should first scan one small number per
position and build the full record only where some sample might carry a variant —
about one position in a hundred — and that nine tenths of that experiment's
6.9-fold saving came from the cheap column existing at all. A record stream
destroys that unless something is done.

The fix is cheap and does not bring the columns back. Carry **that one scalar in
its own parallel stream**, frame for frame with the record stream: two streams,
not fourteen. Its price is visible in the field table — the comparable per-record
scalar compresses to 0.12 bytes a record on tomato — and splitting a batch into
more separately-compressed pieces was measured to *help*, not hurt (on HG002 at
512 KiB, fourteen separate field frames give 5.44 bytes a record against 5.82 for one
combined frame), because grouping like values together is most of what
transposition was doing. A per-frame summary in the tail index would additionally
let whole frames be skipped without decompressing them at all.

## 7. The knobs, with their measured spread

| knob | range measured | effect on size | effect on other things |
|---|---|---|---|
| frame size | 8 – 512 KiB | 4 % (tomato: 7.66 → 7.36) | it *is* the reader memory |
| dictionary size | 0 – 112 KiB | 6 % (tomato: 7.90 → 7.44) | held once per open reader unless shared |
| zstd level | 3 / 9 / 19 | 16 % (tomato: 8.15 / 7.44 / 6.86) | encode 1.7 s / 4.0 s / 16.2 s |
| log-error precision | 1/256 – 1/4 log unit | that field: 1.97 → 1.51 bytes | accuracy of a likelihood term |

Level 9 and a 32 KiB frame are the settings every number in this report was taken
at, and nothing in the sweep argues for moving them: level 19 costs four times
the write for 8 % of the file, and the frame size is flat.

## 8. What else was considered

- **One continuous compressed stream with a bounded sliding window and periodic
  flush points.** This is the shape the question proposed, and it buys nothing
  over small dictionary-primed frames. On HG002, same record encoding throughout:
  a stream with a 512 KiB window gives 6.54 bytes a record, and 32 KiB frames with
  a dictionary give 6.52 — the same size for fourteen times the reader memory.
  Only an 8 MiB window pulls ahead, to 6.25, which is 4 % of the file for 230
  times the memory. The reason is that a window is paid for in memory by every
  open reader on every read, while a dictionary is paid for once, is smaller, and
  can be shared. The stream also gives up free random access between its flush
  points, and flushing every 8 KiB to get it back costs 3 % more bytes.
- **Keeping the transposed layout but shrinking the blocks.** Measured above: it
  is the curve the record stream flattens, and at 64 KiB it is already worse than
  the record stream at 36 KiB on both samples.
- **Row records in gzip blocks — what BAM does.** The same architecture as the
  recommendation with a weaker codec and no dictionary; zstd at level 3 already
  beats gzip on both ratio and speed, so there is no argument for it.
- **Column stores with row groups (Parquet, ORC, Arrow IPC, CRAM).** All make the
  same block-is-the-memory-unit trade the current `.psp` makes; the repeat catalog
  already uses Parquet where that trade is the right one, because its heavy
  readers touch three of seven columns. The per-sample store's readers touch all
  of them.
- **Faster or stronger codecs (lz4, xz).** Not measured. lz4 trades roughly a
  third of the ratio for speed the decode is not short of, and xz trades a large
  multiple of the encode time and far more decoder memory for a modest ratio gain.

## 9. What this does not tell you

- **Decode was measured on one thread, one reader at a time, front to back.** The
  region-seek path — jumping to a coordinate through an index — was not
  implemented or measured. Nothing in the design obstructs it (a frame carries its
  own chromosome and first position, and never crosses a chromosome), but it is
  untested.
- **The record shape measured in §§1–7 is production's, not ng's**, and on one
  field they are not comparable at all: production stores the names of about 3.4 %
  of the reads it folds and ng stores all of them. §12 measures that field on its
  own, from the alignments. The rest of the fields — the summed log-error, the
  mapping-quality moments, the sequences — are the same ones.
- **§12's read-name measurement is derived from alignment geometry, not from ng's
  generator**, which does not yet write a file. It uses the reference stretch each
  read covers and gives each read pair one identifier, which is what ng's own
  allocator does; what it cannot see is how ng splits those names across several
  observations of one record.
- **Both samples are one library, one chemistry.** A cohort of mixed read groups
  will have more distinct values per frame and compress somewhat worse.
- **The quantisation is a lossy change and needs a modelling ruling, not a
  measurement.** The tolerances are stated above; whether 1/256 of a natural-log
  unit in a summed log-error is acceptable is a question about the likelihood, not
  about the file format.

## 10. What to build

1. **Quantise the three floating-point fields**, at the steps §11 settles on. It
   is 44–48 % of the file and independent of every other decision here.
2. **Store the read names as changes, not as a list per position** (§12). At a few
   reads a position it saves about 0.6 bytes a position; at a few hundred it saves
   37. If only one change is affordable, delta-varints alone (§12) capture 60 % of
   that saving at eleven reads a position and 86 % at three hundred, for a few
   lines and no reader state.
3. **Make ng's per-sample store a record stream**: records serialised end to end,
   32 KiB frames, one dictionary in the file header, a tail index of (chromosome,
   first position, byte offset) per frame. No transposition, no column manifest,
   no per-block schema.
4. **Keep one scalar per record in a parallel stream** so the merge can still scan
   cheaply and materialise about one position in a hundred, and put a summary of
   it per frame in the index so whole frames can be skipped.
5. **Share one prepared dictionary across the cohort's open readers**, or shrink
   it to 16 KiB, so the per-reader cost is the frame and not the dictionary.

---

## 11. How finely the three quantities need to be stored (2026-08-19, follow-up)

The question was whether 16 bits would buy anything over 8. **The format has no
field width to choose:** each quantity is stored as an integer count of steps and
written as a variable-length integer, so a value that needs one byte takes one and
a value that needs three takes three, in the same file. What there is to choose is
the **step**, and the compressor then charges for how much the numbers actually
move rather than for a declared width. Swept on the two samples, at 32 KiB frames,
zstd level 9, changing one field at a time:

| quantity | step | HG002 | tomato |
|---|---|---|---|
| GC fraction of the window | 1/100,000 | 6.948 | 7.610 |
| | 1/10,000 (baseline) | 6.727 | 7.443 |
| | 1/1,000 | 6.658 | 7.394 |
| | **1/100** | **6.026** | **6.792** |
| | 1/20 | 5.767 | 6.567 |
| mean coverage of the window | 1/256 read | 7.085 | 7.808 |
| | 1/16 read (baseline) | 6.727 | 7.443 |
| | **1/4 read** | **6.566** | **7.288** |
| | 1 read | 6.480 | 7.198 |
| | 4 reads | 6.456 | 7.160 |
| summed log-error per allele | 1/4,096 ln | 7.047 | 7.774 |
| | 1/256 ln (baseline) | 6.727 | 7.443 |
| | 1/16 ln | 6.415 | 7.091 |
| | 1/4 ln | 6.041 | 6.732 |
| | 1 ln | 5.638 | 6.355 |

Bytes per record for the whole file, so a row's distance from the baseline is what
that field alone contributes.

**Two of the three are free to coarsen and one is not.** The GC fraction feeds a
coverage-against-GC curve that bins its input anyway, and the mean coverage feeds a
ratio of observed depth to expected depth; neither can tell 1 % from 0.01 % of GC,
or a quarter of a read from a sixteenth. Taking the GC fraction to 1 % and the
coverage to a quarter of a read, and **leaving the log-error where it is**, gives
5.843 bytes a record on HG002 and 6.621 on tomato — 11–13 % below the baseline,
for no accuracy anyone downstream consumes.

The summed log-error is the one to leave alone without a modelling ruling: it goes
straight into a likelihood, so a step of 1/16 of a natural-log unit is a 6 % error
in that term where 1/256 is 0.4 %. It is also the field whose magnitude grows with
depth, which is the other reason a fixed 8- or 16-bit width is the wrong frame for
it — at three hundred reads a position the value needs more range than 16 bits
hold, while at three reads it needs six bits.

**Recommended steps:** GC fraction 1/100, mean coverage 1/4 of a read, summed
log-error 1/256 of a natural-log unit until someone rules otherwise. The three
scales are written into the file's own header in the probe, so changing one is a
writer decision and no reader has to be told.

## 12. Storing the read names as changes rather than as a list (2026-08-19, follow-up)

**This is the field that decides ng's file size at depth, and nothing measured in
§§1–7 says about it.** Production names the reads that disagreed with the
reference and drops the rest — about 3.4 % of them — so its read-name field was
0.33 bytes a record on tomato and 0.18 on HG002, a rounding error. ng names every
read at every position it covers, and that changes the field's size by two orders
of magnitude.

Measured directly from the alignments (`examples/ng_chain_id_column_cost.rs`): each
read pair gets one identifier allocated in order, the reference stretches it covers
are worked out from its CIGAR, and the resulting per-position live sets are written
three ways and compressed identically — 32 KiB frames, zstd level 9, a frame cut
every 1,500 covered positions so the restart points match what the record stream
gives.

| | tomato bench slice, 11.4 reads a position | HG002 bench slice, 293 reads a position |
|---|---|---|
| whole list, raw 8-byte identifiers (production's shape) | 1.020 bytes a position | 43.78 |
| whole list, each identifier as its distance from the one before | 0.668 | 11.72 |
| **changes only — who arrived, who left** | **0.432** | **6.42** |

Set against the roughly 5.4 bytes a record the other thirteen fields cost, the
read names in production's shape would be **16 % of ng's file at eleven reads a
position and 89 % of it at three hundred**. Stored as changes they are 7 % and
54 %. So the answer is yes, and the reason to do it is not the tomato corner — it
is that the naïve field grows faster than depth (25.7 times the depth cost 43
times the bytes) while the differential form grows slower than depth (14.9 times).

Four things the measurement settles that the encoding spec listed as open:

- **The saving survives zstd** — the spec's first open question. zstd is already
  extremely good at the repeated lists (679 MB of raw identifiers on tomato became
  7.4 MB, ninety-two fold), which is why the field looks cheap at low depth. It is
  not good enough at three hundred reads a position, where the same collapse only
  reaches thirty-six fold.
- **Delta-varints alone are worth having, and may be enough** — the spec's second
  open question. They capture 57 % of the available saving at eleven reads a
  position and **86 % at three hundred**, with no reader state, no residual
  arithmetic and no new error class. If the differential form is deferred, this is
  the change to make meanwhile.
- **An identifier goes live more than once for most reads** — the spec's fourth
  open question, and the answer is not marginal: **83 % of identifiers on HG002 and
  91 % on tomato** cover two stretches with a gap between them, because a pair's
  mates rarely overlap. An arrivals-and-departures stream that assumes one stretch
  per identifier would lose the second mate of nine reads in ten. It needs a
  re-entry form, and that is not an edge case to handle later.
- **Restating the live set at every frame is affordable.** Cutting a frame every
  1,500 positions instead of by byte count — the shape a shared file has, where the
  record stream decides the cut — costs the differential form 12 % of its own bytes
  on tomato (0.385 → 0.432) and leaves it far ahead of both alternatives.

A dictionary adds little here (0.432 → 0.334 on tomato, nothing on HG002): the
repetition this field carries is between neighbouring positions, which the frame
already contains, not between distant frames.

## 13. Restart points and the index (2026-08-19, follow-up)

The requirement is to be able to start reading part-way through a sample rather
than only at its beginning. **The record stream already gives that at a much finer
grain than asked for, and it costs almost nothing to index.**

Every frame is self-contained by construction: it opens with its chromosome, its
first position and its record count, and every running base — the position, the
coverage difference, the read-identifier difference — restarts at zero inside it. A
frame never crosses a chromosome. So the restart points are the frame boundaries,
and at 32 KiB there are 5,085 of them in a 56 MB tomato sample, one about every
1,500 covered positions. Asking for one every 100 kb would be asking for fewer.

**What has to go in the index at each point:** the chromosome, the first position,
and the byte offset — about ten bytes a frame as variable-length integers, so
roughly 50 KB for that tomato sample, one part in a thousand of the file. The index
sits at the tail, as the `.psp`'s does, so the writer stays single-pass.

**What has to go in the index beyond that, and why:**

- **A summary for skipping.** The one scalar the merge scans (§6) summarised per
  frame — the largest non-reference support in it — lets a reader decide from the
  index alone whether to decompress the frame at all.
- **Any state a frame does not restate itself.** In the store as measured there is
  none. If the differential read-name form of §12 lands, the live read set is
  exactly such state, and §12 measures the cost of restating it at every frame
  rather than carrying it: 12 % of that field on tomato. **This is the one
  interaction between the two proposals**, and it is why the frame boundary and the
  restart point should be the same thing rather than two mechanisms.

**One rule to add for sparse samples.** A frame holds a fixed number of *bytes*, so
on a sample whose coverage is patchy a frame can span a long stretch of reference.
If a guaranteed restart every 100 kb of *reference* is wanted, cut a frame at those
boundaries too. On the tomato sample this would add frames only where coverage is
thin, so its cost is a fraction of the 0.1 % the index already spends.
