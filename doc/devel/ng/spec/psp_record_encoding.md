# ng — how the psp stores one sample's observations

*Status: design spec, 2026-08-19. **No code.** This settles what a psp record costs and
what shape the file has; it does not fix a byte layout. Every number in it was measured
on real files with three throwaway programs, named in §14.*

*This is the document [`run_streaming.md`](run_streaming.md) §10 defers to when it says
"the psp file's encoding — byte layout, compression, block sizing …, checksums, format
versioning, the index, the trailer. Its own spec beside this one". It inherits that
document's header fields (§6.1), its reader contract (§3.3), its per-open-file budget
(§7.2) and its worker-count-invariance restriction (§12.1). The read names are a field
big enough to have their own document —
[`psp_chain_id_encoding.md`](psp_chain_id_encoding.md) — and §7 here reports the
measurement that document's experiment was waiting for.*

---

## 1. What this is

A psp holds, for one sample, everything its reads showed at every position the run
analysed — one record per covered reference position, at three reads a position and at
three hundred, for a cohort of one sample and of several thousand. **This document
settles how much of the disk one of those records costs and how much memory a reader
has to spend to get it back.** Those two are the same question asked twice, and the
current production format ties them together in a way ng should not inherit.

### 1.1 The problem, in one paragraph

Production's `.psp` groups records into blocks of a target uncompressed size
([`TARGET_BLOCK_BYTES`, `src/psp/writer.rs:72`](../../../../src/psp/writer.rs) — 1 MiB),
transposes each block so that every field becomes its own buffer, and compresses each of
those buffers as a separate zstd frame at level 9
([`ZSTD_COMPRESSION_LEVEL`, `src/psp/block.rs:709`](../../../../src/psp/block.rs)). The
block is therefore **both** the furthest back the compressor may look for a repeated
pattern **and** the amount a reader must decode before it can hand out its first record.
Shrinking it to save memory also costs bytes, which is the trade ng inherits if it ports
the shape — and ng cannot afford it, because
[`run_streaming.md`](run_streaming.md) §7.2 requires an open psp to cost *tens of
kilobytes, not megabytes*, at three thousand open samples.

### 1.2 Goals

1. **An open psp costs tens of kilobytes.** Inherited from
   [`run_streaming.md`](run_streaming.md) §7.2, and the constraint the whole shape is
   chosen against.
2. **The file is smaller than production's at every point on that curve**, not merely at
   the memory-hungry end.
3. **A reader can start part-way through**, at a grain no coarser than 100 kb of
   reference, without reading what comes before.
4. **The two-phase decode survives.** [`run_streaming.md`](run_streaming.md) §3.3 makes
   psp mode's saving come from scanning one cheap number per position and building the
   full record only where some sample might vary — about one position in a hundred. An
   encoding that forces a reader to inflate every field to read one throws that away.
5. **Nothing about the writer's scheduling reaches the bytes.** Inherited from
   [`run_streaming.md`](run_streaming.md) §12.1: a frame cut is a function of the
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
- **It does not settle the read names' encoding.** That is
  [`psp_chain_id_encoding.md`](psp_chain_id_encoding.md)'s experiment; §7 here supplies
  numbers it was explicitly waiting for and states a leaning, nothing more.
- **It does not specify compression of anything but the observation stream** — the
  header stays plain text so that `head` and a TOML parser can read it, as production's
  does.

### 1.4 Vocabulary, defined once

- A **record** is one sample's observations at one covered reference position:
  [`SampleLocusObservations`](../../../../src/ng/locus_generation/mod.rs), which holds the
  region, the reference bases over it, a list of observed sequences with their support,
  and two counts of reads that produced no observation.
- **Transposed** describes production's shape: inside a block, every field of every record
  is gathered into its own buffer and compressed separately.
- A **record stream** is the alternative this document proposes: each record's fields
  written one after another, the bytes cut into **frames**, each frame compressed on its
  own.
- A **dictionary** is a few tens of kilobytes of representative bytes stored once in the
  file and handed to the compressor before every frame, so that a small frame is not
  compressed from a cold start. It is a standard zstd facility, not something we build.
- **Reader memory** is the quantity every choice here is measured against: the bytes one
  open sample forces a reader to hold before it can produce a record.

---

## 2. What a record costs today, and which fields it is

Before any layout question, the fields are not equal, and the two largest are not the ones
anybody would guess. Per record, compressed, measured on a real per-sample file at 512 KiB
blocks — GIAB HG002 at about thirty reads a position:

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
which is the one thing a compressor cannot shrink. §7 is about that, and it is worth more
than every layout decision in this document put together.

**One field measured here says nothing about ng, and it is the one that matters most at
depth.** Production stores the names of about 3.4 % of the reads it folds
([`src/pileup/walker/open_record.rs:150-160`](../../../../src/pileup/walker/open_record.rs)),
ng stores all of them, and §8 measures that field separately.

---

## 3. The proposal, in one page

**Write the records as a stream, not as a grid.** Each record's fields go one after
another; the bytes are cut into small frames; each frame is compressed on its own against
a dictionary stored in the file. Beside it runs a second, tiny stream carrying the one
number the merge scans. A tail index names every frame.

```
psp file
  header                  plain-text, the fields run_streaming.md §6.1 fixes
  dictionary              one, ~16-112 KiB, trained on this file's own frames
  frame 0   light         the scan scalar for each record in the frame
            heavy         the records themselves
  frame 1   light
            heavy
  ...
  index                   one entry per frame: chromosome, first position,
                          the two byte offsets, and the frame's scan summary
  trailer                 locates the index; its absence means "interrupted"
```

Everything below is a consequence of that shape or a number that justifies it.

### 3.1 Why a stream rather than the transposed blocks — priced

Measured like for like, same field encodings and same compressor, on a tomato accession at
about three reads a position:

| shape | reader memory | bytes a record |
|---|---|---|
| transposed, one frame per field | 1024 KiB | 6.27 |
| transposed, one frame per field | 256 KiB | 6.60 |
| transposed, one frame per field | 64 KiB | 7.33 |
| **record stream, 32 KiB frames + dictionary** | **36 KiB** | **7.39** |
| record stream, 8 KiB frames + dictionary | 12 KiB | 7.57 |

**Dropping the transposition costs about 18 % of the bytes and buys about one
twenty-eighth of the memory.** That is the whole trade, and goal 1 decides it: three
thousand open samples at a megabyte each is 3 GB before a single record is read, at 36 KiB
each it is 108 MB.

The second property is what makes the small frames work at all: **the record stream is
nearly flat in frame size** — 132 KiB of reader memory down to 12 KiB costs 3 %, where the
transposed shape over the same range costs a third of its size. A design whose memory knob
is nearly free is a different kind of object from one whose memory knob is its ratio.

Against the file production writes today (11.85 bytes a record on that sample, 11.86 on
HG002), the record stream is 37 % and 43 % smaller — but most of that is §7's, not this
section's.

### 3.2 Frames with a dictionary, not one long stream with a sliding window

The obvious alternative to independent frames is one continuous compressed stream whose
sliding window is capped, flushed periodically so a reader can resume. **Measured, it buys
nothing.** On HG002 with the identical record encoding, a stream with a 512 KiB window
gives 6.54 bytes a record and 32 KiB frames with a dictionary give 6.52 — the same size for
fourteen times the reader memory. Only an 8 MiB window pulls ahead, to 6.25, which is 4 % of
the file for 230 times the memory.

The reason is structural rather than incidental: **a window is paid for in memory by every
open reader on every read; a dictionary is paid for once, is smaller, and can be shared
between readers.** The stream also gives up free random access between its flush points,
and flushing every 8 KiB to get it back costs 3 % more bytes.

The dictionary itself is a knob, not a fixture: at 32 KiB frames it is worth 6 % of the
file (tomato 7.90 → 7.44 bytes a record with a 112 KiB dictionary), and a 16 KiB
dictionary keeps most of that (7.66) at a seventh of the per-reader cost. **A file with no
dictionary at all still beats production**, so this is a tuning decision and not a
dependency.

### 3.3 The light stream — how the two-phase decode survives

A record stream's one real loss is that a reader cannot read one field without decoding all
of them, and [`run_streaming.md`](run_streaming.md) §3.3 depends on exactly that: psp mode's
saving is scanning a cheap number per position across the cohort and inflating only the
positions some sample varied at — production's `TwoPhaseSegment` decodes the light columns
for every row and leaves the heavy ones compressed
([`src/var_calling/sample_reader.rs:698-712`](../../../../src/var_calling/sample_reader.rs)),
then inflates the kept rows
([`:789`](../../../../src/var_calling/sample_reader.rs)). An earlier measurement of that
saving on the cohort reader materialised 28,718 loci instead of 2.83 million.

**So the file has two streams, not fourteen and not one.** The light stream carries, per
record, what phase one needs — the position, the record's reference span, and the summed
non-reference support — and nothing else. The heavy stream carries the record.

This costs almost nothing and is not a compromise between the two shapes: **splitting a
frame into more separately-compressed pieces was measured to help, not hurt** (on HG002 at
512 KiB, fourteen field-frames give 5.44 bytes a record against 5.82 for one combined
frame), because putting like values next to each other is most of what transposition was
doing. The nearest field measured on its own, a small per-record varint, compressed to
0.12 bytes a record on tomato; the light stream is a few of those.

**The unit of the light stream is the record, and its frames are cut with the heavy
stream's** — one index entry names both. A reader that scans a segment reads only light
frames; a reader that then needs three records decompresses the three heavy frames holding
them.

### 3.4 The index, and where a reader may start

**Every frame is self-contained by construction.** It opens with its chromosome, its first
position and its record count, and every running base inside it — the position difference,
the coverage difference, the read-name difference — restarts at zero. A frame never crosses
a chromosome. So the restart points are the frame boundaries, and there is no separate
mechanism to build.

At 32 KiB frames the measured tomato sample — 7.59 M covered positions, 56 MB — has 5,085
of them, one about every 1,500 covered positions, **finer than the 100 kb the requirement
asks for.** An entry naming a frame's chromosome, first position and two byte offsets is
about ten bytes as variable-length integers: roughly 50 KB for that sample, one part in a
thousand of the file.

Two rules the index has to carry beyond that:

- **A per-frame summary of the scan scalar** — the largest non-reference support anywhere in
  the frame — reachable **without decompressing the frame**, so a reader can decide whether
  to touch it at all. Whether that lives in the index or in the frame's own uncompressed head
  depends on how the index is sized, which is open (question 4 below).
- **A frame is cut at every 100 kb of reference as well as at its byte target.** A frame
  holds a fixed number of *bytes*, so on a sample whose coverage is patchy it can span a
  long stretch of reference and the guarantee in goal 3 would not hold. This costs frames
  only where coverage is thin.

**That 50 KB is a sliced benchmark sample, and the whole-genome number is the one that
matters.** The measured file covers about one per cent of the tomato genome. A whole-genome
sample at three reads a position covers most of 800 Mb, so at 1,500 positions a frame it
holds roughly **500,000 frames — an index of about 5 MB per open file**, which is *worse*
than the 3.8 MB [`run_streaming.md`](run_streaming.md) §7.2 rejects when it measures
production's flat vector of one 24-byte entry per block
([`BlockIndexEntry`, `src/psp/index.rs:42`](../../../../src/psp/index.rs), decoded whole at
open, [`decode_index`, `:110`](../../../../src/psp/index.rs)).

**So finer frames make the index problem worse, and this design does not solve it — it
inherits §7.2's shape and owes the sizing.** Index at a coarse grain and chain the frames
within it, each frame carrying enough to reach the next, so a reader seeks once and then
walks. **This is the one requirement in this document with no measurement behind it**, and
it should be settled before the writer is built rather than after (§13, open question 4).

---

## 4. Settings, and what each is worth

Every number below was measured on both corners; the spread is what makes each of these a
knob rather than a constant.

| setting | proposed | range measured | what it moves |
|---|---|---|---|
| frame size | 32 KiB | 8 – 512 KiB | 4 % of the file (tomato 7.66 → 7.36); it **is** the reader memory |
| dictionary | 16 KiB, shared across readers | 0 – 112 KiB | 6 % of the file (tomato 7.90 → 7.44); held per reader unless shared |
| zstd level | 9 | 3 / 9 / 19 | 16 % of the file (tomato 8.15 / 7.44 / 6.86); write time 1.7 s / 4.0 s / 16.2 s |

**Level 9 is inherited, not derived** — it is what production uses
([`src/psp/block.rs:709`](../../../../src/psp/block.rs)) and the sweep gives no reason to
move: level 19 costs four times the write for 8 % of the file. **The frame size is a
genuine choice and 32 KiB is the middle of a flat region**, not an optimum. Both are safe
to move later; neither is a format change if the reader takes what the file says.

---

## 5. What is *not* stored

[`run_streaming.md`](run_streaming.md) §11 question 4 asks whether an observation's
reference bases should go in the file at all —
[`SampleLocusObservations::reference_bases`](../../../../src/ng/locus_generation/mod.rs) is
a `Box<[u8]>` per record, so written out it is a per-sample copy of the reference over the
analysed ground.

**Leaning: do not store them, and re-fetch when a frame is decoded.** The measurement
supporting it is indirect but firm: the allele sequences, which are the same kind of
content, compressed to 0.32–0.34 bytes a record on both corners, and the reference bases
are longer than the allele on every record where a deletion widened the footprint. The
calling stage holds the reference already. **Confirm before code** by writing one sample
both ways and timing the re-fetch inside a segment walk; the question is
[`run_streaming.md`](run_streaming.md)'s, not this document's, and this is a leaning
offered to it rather than a decision taken from it.

---

## 6. The three quantities that arrive as floating point

**They are 70 % of the compressed file (§2) and none of them needs the precision it is
stored at.** Each is written instead as an integer count of steps, as a variable-length
integer, so a value needing one byte takes one and a value needing three takes three in the
same file. **There is no field width to choose** — the question that looks like "8 bits or
16?" is really "how big is the step?", and the compressor then charges for how much the
numbers move rather than for a declared width.

Swept one field at a time, at 32 KiB frames and level 9, in bytes a record for the whole
file (so a row's distance from the baseline is that field's own contribution):

| quantity | step | HG002 | tomato |
|---|---|---|---|
| GC fraction of the window | 1/100,000 | 6.948 | 7.610 |
| | 1/10,000 | 6.727 | 7.443 |
| | 1/1,000 | 6.658 | 7.394 |
| | **1/100 — proposed** | **6.026** | **6.792** |
| | 1/20 | 5.767 | 6.567 |
| mean coverage of the window | 1/256 read | 7.085 | 7.808 |
| | 1/16 read | 6.727 | 7.443 |
| | **1/4 read — proposed** | **6.566** | **7.288** |
| | 1 read | 6.480 | 7.198 |
| summed log-error per allele | 1/4,096 ln | 7.047 | 7.774 |
| | **1/256 ln — proposed** | **6.727** | **7.443** |
| | 1/16 ln | 6.415 | 7.091 |
| | 1/4 ln | 6.041 | 6.732 |

**Two of the three are free to coarsen and one is not.** The GC fraction feeds a
coverage-against-GC curve that bins its input, and the mean coverage feeds a ratio of
observed to expected depth; neither consumer can tell 1 % of GC from 0.01 %, or a quarter of
a read from a sixteenth. Taking those two to the proposed steps and leaving the log-error
alone gives **5.843 bytes a record on HG002 and 6.621 on tomato — 11–13 % below the
baseline, for no accuracy anything downstream consumes.**

**The summed log-error is the one that needs a ruling and does not have one.** It goes
straight into a likelihood: a step of 1/16 of a natural-log unit is a 6 % error in that
term where 1/256 is 0.4 %. It is also the field whose magnitude grows with depth, which is
the second reason a fixed 8- or 16-bit width is the wrong frame for it — at three hundred
reads a position the value needs more range than sixteen bits hold, while at three reads it
needs six. **Proposed: 1/256 of a natural-log unit until the modelling side rules
otherwise** (§13, open question 1); coarsening it to 1/16 would buy a further 5 %.

**The three steps are written into the file.** A reader never has to be told them, so
changing one is a writer decision and not a format change. Encoding of a missing value: a
window that does not exist — an `N` reference position — is a real state and not a zero, so
the code 0 is reserved for it and every present value is shifted by one.

**Round-tripping this is lossy by construction and must be checked as such**, which is
§12's first oracle: every integer field, allele sequence and read-name list comes back
identical, and these three come back inside their own step. Measured over 7.59 M tomato
records at the proposed steps: worst error 0.005 in GC fraction, 0.125 of a read in
coverage, 0.002 natural-log units in the summed log-error.

---

## 7. The read names, and why they decide the file's size at depth

**ng names every read at every position it covers.** This is the owner's ruling of
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
compressed identically at 32 KiB frames with a frame cut every 1,500 positions:

| | tomato slice, 11.4 reads a position | HG002 slice, 293 reads a position |
|---|---|---|
| whole list per position, raw 8-byte identifiers | 1.020 bytes a position | 43.78 |
| whole list, each identifier as its distance from the one before | 0.668 | 11.72 |
| **only the changes — who arrived, who left** | **0.432** | **6.42** |

*The deep corner is the HG002 benchmark slice, whose reads are concentrated into 1,000 small
regions; the depth is real and inside the committed range, but produced by that selection
rather than by a 300× library.*

Set against the roughly 5.4 bytes a record everything else costs, **the read names in
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
- **Restating the live set at every frame is affordable.** Cutting frames every 1,500
  positions rather than by byte count — which is what §3.4 does — costs the differential form
  12 % of its own bytes on tomato (0.385 → 0.432) and leaves it far ahead of both
  alternatives.

**Proposal: delta-varints in the first version, the differential form second.** The first is
a few lines inside the record encoder, has no interaction with anything else in this
document, and takes the deep corner from 43.8 to 11.7 bytes a position. The second is worth
another 5.3 bytes a position there and needs the live set restated at every frame, which is
the one place these two documents' designs touch (§3.4).

---

## 8. Traps — what will bite the coder

- **[`run_streaming.md`](run_streaming.md) §8 says chain ids are already omitted for reads
  that agree with the reference, "recorded so nobody re-adds them when the psp writer is
  built". That bullet is stale and following it produces a file the cohort merge cannot
  use.** The code names every read
  ([`fast_column.rs:315`](../../../../src/ng/locus_generation/pileup/fast_column.rs),
  [`open_record.rs:494`](../../../../src/ng/locus_generation/pileup/open_record.rs)), by the
  2026-08-17 ruling that post-dates the bullet. Fix the bullet when that document is next
  touched.
- **A dictionary is held per open reader unless it is shared.** 112 KiB across three
  thousand samples is about 330 MB — larger than the frames the whole design exists to
  shrink. zstd can prepare one dictionary and lend it to every reader
  (`Decompressor::with_prepared_dictionary`); use that, or a 16 KiB dictionary, or both.
  Measured: 63 open samples cost 22 MB, about 312 KiB each, of which the 112 KiB dictionary
  is a third and the rest is the read buffer, the decompressed frame and the decompressor's
  own workspace.
- **A dictionary trained on the frames it is then measured against reports a saving no
  reader ever gets.** This is easy to do by accident and the number it produces is
  spectacular — in an early run of our own probe, a hundredfold. Train on one half of the
  frames and measure on the other.
- **A frame's size is in bytes, so its span in reference is whatever coverage makes it.**
  §3.4's 100 kb cut rule exists for that; without it goal 3 quietly fails on exactly the
  sparse samples that need it.
- **The three quantised fields are lossy and the integer fields are not.** A round-trip test
  that compares whole records with a tolerance will pass while the read-name list is being
  corrupted. Compare the integer fields, the sequences and the name lists **exactly**, and
  only the three quantised fields against their step.
- **The window's mean coverage is stored as a difference from the previous record**, so a
  frame that does not reset that difference reads back wrong from its first record — and
  plausibly, because coverage is smooth. The same applies to the position and the read-name
  bases. Every running base resets at a frame boundary; this is the property §3.4's restart
  guarantee rests on.
- **Do not decode the whole index at open** — §3.4, and the 3.8 MB per file
  [`run_streaming.md`](run_streaming.md) §7.2 measures on production's.

---

## 9. Cross-cutting concerns

**Memory.** The point of the design. One open sample holds: the compressed frame it is
reading, the decompressed frame (32 KiB), one record, and its share of the dictionary. The
whole-cohort figure was measured — 63 tomato accessions open at once and advanced one record
each per round, the shape the merge reads in: **22 MB peak against 170 MB for the same walk
over the `.psp`s**, and the psp side is flattered there because those particular files
average 258 KiB of block rather than the writer's 1 MiB default.

**Speed.** Decoding was measured at 20.7 M records a second against 12.6 M for the `.psp`
reader on the same sample, and the cohort walk finished in 22 s against 38 s. The stream
never builds a field array it must then walk again to assemble records. This is a
single-threaded sequential-read measurement and says nothing about the seek path.

**Errors.** A frame that decompresses to the wrong length, a record that runs off the end of
its frame, a read-name list longer than the observation's read count: all are corrupt-input
failures belonging to the psp reader's error type, and none may reach the merge as a
half-built record. The trailer's absence means the writer was interrupted and the file is
refused rather than read as a short sample —
[`run_streaming.md`](run_streaming.md) §9's rule, unchanged.

**Concurrency.** Frames are independently decodable, so nothing here serialises. The frame
cut must depend only on the observation stream, never on scheduling, or the byte-identity
oracle of [`run_streaming.md`](run_streaming.md) §12.1 breaks — a byte target counted over
the records as they arrive satisfies this; a flush driven by a timer or a queue depth does
not.

---

## 10. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| variable-length integer codec | `src/psp/varint.rs` | as-is: LEB128 and zig-zag LEB128, already specified and tested |
| plain-text header framing | `src/psp/header.rs` | the pattern — magic, length prefix, TOML, sentinel — so `head` still works on an ng psp |
| the zstd seam | `new_column_compressor`, `zstd_compress_into` ([`src/psp/block.rs:718,730`](../../../../src/psp/block.rs)) | the shape: one long-lived compressor per writer, frame checksums on. The dictionary is new |
| the two-phase read | `TwoPhaseSegment`, `set_variable_rows` ([`src/var_calling/sample_reader.rs:698,789`](../../../../src/var_calling/sample_reader.rs)) | the light/heavy split, reduced from fourteen columns to two streams (§3.3) |
| the eager whole-segment decode | [`sample_reader.rs:20-26`](../../../../src/var_calling/sample_reader.rs) | the parity oracle's model: a simple decoder used only by tests, against which the real one is byte-compared |
| the record | `SampleLocusObservations` ([`src/ng/locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | what is written and what must come back |

**Parity oracle:** the same shape production uses — a deliberately simple whole-file decoder
used only by tests, and the real reader compared against it record for record.

---

## 11. Deferred, with a recommended home

- **The byte layout itself** — field order inside a record, the framing integers' widths, the
  trailer's bytes, the format version tag. To the implementation, guided by this document;
  it is small enough not to need its own spec once the shape is fixed.
- **How the index stays small at 100,000 frames a file** — the coarse-index-plus-chaining
  shape [`run_streaming.md`](run_streaming.md) §7.2 asks for. Named as open question 4 below
  rather than deferred to another document, because nothing else can settle it.
- **The read names' final encoding** — [`psp_chain_id_encoding.md`](psp_chain_id_encoding.md),
  whose experiment §7 here feeds rather than replaces.
- **Whether the reference bases are stored** — [`run_streaming.md`](run_streaming.md) §11
  question 4; §5 offers a leaning and the measurement that would close it.
- **Multi-library samples.** Both corners measured here are one read group. A sample with
  several will have more distinct values per frame and compress somewhat worse; nobody has
  measured how much. To the first run on a multi-library cohort.

---

## 12. How we know it works

1. **Round-trip, with the right strictness per field.** Every integer field, every allele
   sequence and every read-name list identical; the three quantised fields inside their own
   step. Already demonstrated on 7.59 M tomato records and 0.60 M HG002 records with the
   probe of §14.
2. **The per-open-file budget, measured rather than argued.** N samples open and walked in
   lockstep, peak resident reported, against [`run_streaming.md`](run_streaming.md) §7.2's
   "tens of kilobytes". The probe does this at 8, 32 and 63 samples today.
3. **Worker-count invariance**, inherited from [`run_streaming.md`](run_streaming.md) §12.1:
   one sample gathered at 1, 2, 4, 8, 16 workers gives byte-identical files apart from the
   header's timestamp. This is what §9's frame-cut rule exists to preserve.
4. **The two-phase saving is still there.** Scanning the light stream over a cohort segment
   must materialise about one record in a hundred, the ratio
   [`run_streaming.md`](run_streaming.md) §3.3 measured. A file that round-trips but forces
   every record to inflate has failed goal 4 while passing oracle 1.
5. **Restart equals sequential.** Reading from an arbitrary frame gives exactly the records a
   full sequential read gives from that point — the test that catches a running base that was
   not reset (§8).
6. **Mode equivalence**, [`run_streaming.md`](run_streaming.md) §12 oracle 3: the same cohort
   called through the direct route and through the psp route gives the same VCF. It is the
   sufficiency test for everything this document chose not to store.

---

## 13. Open questions

1. **What step may the summed log-error be stored at?** — OPEN, and it is a modelling
   question, not an encoding one. The term enters a likelihood; 1/256 of a natural-log unit is
   a 0.4 % error in it and 1/16 is 6 %. *Leaning:* 1/256 — a 5 % file saving is not worth an
   unquantified change to a likelihood. **Settled by:** whoever owns the emission model saying
   what error in that term is acceptable.
2. **Does the differential read-name encoding ship in the first version, or only
   delta-varints?** — OPEN. §7 measures both. *Leaning:* delta-varints first, the differential
   form second, because the first is a few lines with no interaction and takes the deep corner
   from 43.8 to 11.7 bytes a position. **Settled by:**
   [`psp_chain_id_encoding.md`](psp_chain_id_encoding.md) §7's experiment, which now has its
   sizes and needs its decode and merge times.
3. **Are the light and heavy streams interleaved frame by frame, or separated into two
   regions of the file?** — OPEN, and not measured. Interleaving keeps a segment's light and
   heavy adjacent, so serving one segment is one seek; separating makes a whole-file light scan
   one sequential read. *Leaning:* interleave, because
   [`run_streaming.md`](run_streaming.md) §3.4's reader serves segments rather than files.
   **Settled by:** timing a cohort merge over one chromosome both ways, once a writer exists.
4. **How does the index stay small when a file has 100,000 frames?** — OPEN, and the one
   requirement in this document with no measurement behind it. Ten bytes an entry is 1 MB per
   open sample, which fails
   [`run_streaming.md`](run_streaming.md) §7.2 at three thousand samples. *Leaning:* index at a
   coarse grain — one entry per megabase, say — and let each frame carry the offset of the
   next, so a reader seeks once and then walks; the per-frame scan summary then has to live in
   the frame's own head rather than in the index, which costs a decode per skipped frame unless
   the summaries are gathered into the coarse entry. **Settled by:** pricing both against a
   real segment-serving read, and it should be settled before the writer is built, not after.
5. **How much do the numbers move on a multi-library sample?** — OPEN. Both corners here are
   one read group, and the read group joins an observation's identity
   ([`src/ng/locus_generation/mod.rs:245`](../../../../src/ng/locus_generation/mod.rs)), so a
   multi-library sample has more, smaller observations per record. *Leaning:* somewhat worse,
   not structurally different. **Settled by:** running the probe on a multi-library sample.

---

## 14. Where the numbers came from

Two throwaway programs, both in the tree, both runnable:

- **`examples/psp_record_stream_compression.rs`** — sweeps the grid of §3.1 and §3.2: field
  order (records or transposed), field width (fixed-width, variable-length, fixed-point),
  batch size, and framing (independent frames, with or without a dictionary; one continuous
  stream with a chosen window, flushed or not). Reports bytes a record against the reader
  memory each combination implies, and each field's own compressed bytes.
- **`examples/psp_row_stream_roundtrip.rs`** — a working encoder, decoder and verifier for the
  proposed shape, plus the many-sample walk of §9. `verify` walks the new store and the `.psp`
  in lockstep and fails on the first record that disagrees.
- **`examples/ng_chain_id_column_cost.rs`** — §7's measurement, taken from a CRAM: read pairs
  given one identifier each, their reference stretches from the CIGAR, the per-position live
  sets written three ways.

The full report, with the sweeps these tables summarise, is
[`../research/per_sample_record_store_compression_2026-08-19.md`](../research/per_sample_record_store_compression_2026-08-19.md).

**What was measured on, and what that covers.** A tomato accession at about three reads a
position (`SRR7279540`, 7.59 M records) and GIAB HG002 at about thirty (0.60 M records) for
everything but §7; for §7, a tomato benchmark slice at 11.4 reads a position and an HG002
benchmark slice at 293. Cohort figures are the 63-accession tomato panel. **One sample and one
cohort of 63 is not the committed range** — nothing here was measured at one sample of a
thousand-sample cohort, and the memory claims at three thousand samples are arithmetic from the
per-sample figure, not measurements.
