# ng — cutting a cohort run's memory by redesigning the per-sample store: experiment plan

*Draft, 2026-08-25. **Experiments only — no writer, no format, no production code.** This plan
runs three measurements and returns one recommended configuration with the curve behind it. The
designs being measured are [`../spec/psp_record_encoding.md`](../spec/psp_record_encoding.md) and
[`../spec/psp_chain_id_encoding.md`](../spec/psp_chain_id_encoding.md); the requirement they both
serve is a per-open-sample memory budget, which the owner reset on 2026-08-25 and which is stated
below rather than taken from [`../spec/run_streaming.md`](../spec/run_streaming.md) §7.2. A question this plan cannot
answer from those three goes back to them.*

---

## What this is for

**ng has to store one sample's observations for a caller that runs from one read a position to
three hundred, on any species.** A per-sample store is a single object, but the thing that
dominates its size changes completely across that range: at low depth almost all of the file is
three quantities that arrive as floating-point numbers, and at high depth almost all of it is the
list of which reads were present. **A configuration chosen at one depth is therefore not a
configuration — it is a point on a curve nobody has drawn.**

So this plan does not measure "low coverage and high coverage". It measures a ladder, and the
deliverable is the shape of the curve over that ladder, on two species.

### The three things being tried

1. **Large blocks with streaming decompression.** Today's production `.psp` makes a reader
   decompress a whole block before it can hand out its first record, so the block is both how far
   back the compressor may look for a repeat *and* how much memory a reader spends. Decompressing
   a block incrementally breaks that link: the memory becomes the compressor's *reach* — a number
   we set at write time and can cap at, say, 32 kB — while the block itself can be a megabyte.

   **The target to aim at, stated by the owner on 2026-08-25: a reader that holds one observation
   per open file at a time.** Nothing reaches that exactly — a compressor needs its reach and a
   decoder needs somewhere to put its output — but it is the direction every choice here is judged
   against, and it is why the number that matters is what a reader *holds*, not what the file
   weighs.
2. **Approximate the three floating-point quantities.** The window's mean coverage, the window's
   GC fraction and the summed log-error per allele are stored at full precision, and their bottom
   bits are arithmetic noise. Stored instead as an integer count of a chosen step, they are much
   smaller — and one of the three feeds a likelihood, so how coarse it may be is a modelling
   question, not an encoding one.
3. **Encode the read names as changes.** ng names every read at every position it covers. A read
   covering 150 positions is therefore named 150 times. Storing instead who arrived and who left
   costs about two entries for the whole read.

### What each measurement must return

For every rung of the depth ladder, on both species, and for every configuration:

| quantity | how it is taken | why this one |
|---|---|---|
| **peak resident of a cohort run at N samples** | the whole caller run, N samples, peak RSS | **the objective** — everything else is how it was bought |
| **compressed bytes per record** | the whole file divided by its records | what the memory is paid for in |
| **compressed bytes per covered position** | the whole file divided by positions | diverges from the above as depth rises, because a record gets fatter rather than more numerous |
| **compressed bytes per field** | each field's own frames, compressed alone | which field to attack next; the answer changes with depth |
| **measured peak resident per open sample** | N samples opened and walked in lockstep, peak RSS, divided by N | the budget below is about this and nothing else |
| **decode throughput** | records a second, sequential; and the cost of a seek | a store that is small and slow is not a win |
| **CPU-seconds per million records decoded** | total processor time, not wall time | wall time can be bought with cores; this cannot, and it is what a cohort run pays |
| **thread scaling of decode** | wall time at 1, 2, 4, 8, 16 threads | says how much of the wall-time cost cores can actually buy back |
| **the growth exponent** | slope of log(bytes per record) against log(depth) | the one number that says whether a design survives to 300× |

**The growth exponent is the headline and the reason the ladder exists.** It is already known to
separate the read-name encodings: measured between 11.4 and 293 reads a position, storing whole
lists of raw identifiers cost 43 times the bytes for 25.7 times the depth, while storing only the
changes cost 14.9 times. A design whose cost grows faster than depth fails at the top of the range
however good it looks at the bottom.

### The objective, and the one thing that could cap it

**The objective is to reduce the peak memory a cohort run needs, and file size is the constraint
rather than the prize** (the owner, 2026-08-25). Everything below is measured in bytes as well,
because bytes are what a smaller memory footprint is usually paid for in — but the number this
programme is judged on is **peak resident of a cohort run at N samples**, not bytes a record.

**The design being tried aims squarely at what is already known to be production's largest
memory knob.** Of every producer-side setting swept on a 50-sample tomato run, the psp block's
reference window is the only one that moves peak resident, and it moves it a long way — 80 kb
blocks cost 3,237 MB and 5 kb blocks 1,062 MB. The default was lowered from 20 kb to 5 kb for
exactly that reason: **58% less memory, for 6.8% more wall time and 13% more disk.** That trade —
memory bought with bytes — is the tie every experiment here exists to cut, because a block is
simultaneously how far back the compressor may look and how much a reader must hold.

**But part of production's peak is not the store, and how much decides the ceiling on all of this.**
Attributed at 50 samples: per-sample chunks held for the whole cohort about 53% of live heap (the
read names alone 33.5%), the per-group merger's projections about 30%, the posterior fit about 8%.
None of those three is a decode buffer. **If the part a store redesign controls is a third of the
peak, then a store costing a tenth of the memory reduces a cohort run's peak by at most 30%** — a
real result, but not the one anybody would assume from a table of bytes a record.

So **Milestone Z below runs before anything else in this plan**, and takes half a day.

### The memory budget this is measured against

**Working budget: 500 kB resident per open sample** — 1.5 GB across three thousand samples — set by
the owner on 2026-08-25, who priced 450 MB across a thousand samples as comfortable. This is much
looser than the "tens of kilobytes, not megabytes" that
[`run_streaming.md`](../spec/run_streaming.md) §7.2 asks for, and **the difference reopens a design
choice that document treated as settled** (§B0 below). It is a working figure, not a ruling: every
sweep reports the whole curve of bytes against memory, so the point on it can be moved without
re-running anything.

### One thing that is settled and must not be reopened

**ng names every read at every position it covers, and that is a requirement of the cohort merge,
not an encoding preference.** A cohort locus can span several of one sample's records, so a read
that covered a position and agreed with the reference has to be distinguishable from a read that
never reached it; unnamed, those two are the same absence. Every experiment here changes how that
column is *stored* and none of them may change what the merge is given back.

---

## Scope

**In — the probe programs and the data ladder they run on:**

- `benchmarks/ssr_hg002/src/subsample_coverages.sh` — extended to reach 1× and 3×.
- A second ladder for tomato, built the same way from `benchmarks/tomato_big_cram/DRR000741.p1.cram`.
- `examples/psp_record_stream_compression.rs` — the configuration sweep. Gains: records sourced
  from ng's own locus generator rather than a production `.psp`; block size and compressor reach
  as independent knobs; per-rung reporting.
- `examples/psp_row_stream_roundtrip.rs` — the working encoder, decoder and verifier, and the
  many-samples-open walk that produces the measured memory figure.
- `examples/ng_chain_id_column_cost.rs` — the read-name forms, swept against block size.
- `examples/blocksize_rewrite.rs`, `examples/psp_rechunk.rs`, `examples/psp_block_stats.rs` — already
  written, and what puts production's columnar blocks into the grid as an arm rather than a footnote.

**Out — later plans, and nothing here may quietly become them:**

- **The psp writer itself**, its byte layout, its version tag and its trailer. This plan chooses a
  configuration; building it is the next plan's job.
- **The index's shape.** Milestone B measures how many blocks a file has at each block size, which
  is what decides whether the index is a problem at all. If large blocks make it small enough, the
  coarse-index-and-chain scheme [`run_streaming.md`](../spec/run_streaming.md) §7.2 asks for is not
  needed and should not be built on speculation.
- **Whether the reference bases are stored** — [`run_streaming.md`](../spec/run_streaming.md) §11
  question 4. Independent of all three experiments here.
- **Any change to production's `.psp`.** Two findings here would apply to it unchanged. That is not
  a reason to touch it.

## Where this runs, and the portability that requires

**Two machines, and the split is set by memory.** Encoding sweeps, probe development and everything
at 63 samples or fewer run on the development Mac, whose container is capped at 16 GB and 8 CPUs.
**Every measurement at 500 samples or more runs on a Linux server with 32 threads and 128 GB**,
which the owner has made available for this — a 1,000-sample run was measured at about 30 GB, so it
fits with room to spare.

**That server is offered on the condition that the experiment is easy to port, so portability is a
requirement of every step here and not an afterthought:**

- **⚠ The Mac has six fast cores, not eight.** The rest are energy-saving cores, and the container
  wrapper hands the runtime `--cpus 8` by default, which oversubscribes the fast ones and makes any
  thread-scaling curve taken here meaningless. **Cap Mac runs at four threads, and take every
  thread-scaling measurement on the 32-thread Linux server.** A wall-time figure from this machine
  carries its thread count and the caveat, or it is not reported.
- **One script per experiment, taking paths as arguments** — no hard-coded locations, no dependence
  on the developer's home directory.
- **Assume no container runtime.** The dev container wrapper refuses to run where neither podman nor
  Apple's `container` exists, and a Linux box may well have neither. Every script must build with
  plain `cargo` as an alternative, **and must look for its binary in both `target/release/` and
  `target-container/release/`, taking the newer** — the container build points its output elsewhere,
  so a script that checks one tree silently runs a stale binary.
- **Name the data it needs and check for it up front**, failing with the missing path rather than
  part-way through a long run. The alignment files are gigabytes, are not in git, and have to be
  staged on the server by hand.
- **Write results as a table to a named output path**, not to standard output, so a run that takes
  hours is not lost to a closed terminal.
- **Record the machine, the thread count and the resolved data paths in the results header.** A
  memory figure without the machine on it is not a measurement, and this plan produces figures from
  two machines that will end up in the same table.

## Principles — how the order was chosen

- **One grid, not three studies.** The three savings interact: at 1× the floating-point fields are
  most of the file and the read names are nearly free, and at 300× it is the other way round.
  Measuring each alone and adding the savings would overstate the total. Every configuration is
  measured on the same records, at every rung.
- **The multiplier before the sweep it multiplies.** How much a streaming decompressor costs to hold
  is what decides how many separately-compressed streams a block can afford, so B1 is measured before
  B0's grid is run, even though B0 states the question. It is a half-day measurement, not a milestone.
- **Measure after compression, and measure memory rather than computing it.** Every existing
  memory figure for the streaming shape is arithmetic — the compressor's reach plus the record
  being built. A real decoder's resident cost is what the requirement is about.
- **Separate the field encoding from the layout, always.** Most of the difference between what
  production writes today and anything proposed here is the fields, not the shape — 44% of the file
  on tomato. A comparison that bundles them credits a layout change with a saving a two-line field
  change also delivers (B0b).
- **Both species at every rung.** Read length, repeat content and GC differ, and all three change
  what compresses. A result on one species is a result about that species.

---

## Milestone Z — how much of production's peak the store actually owns

**This runs before everything else and sets the budget for the rest of the plan's effort.** Half a
day, on code that already exists.

Run production's cohort caller at **50 and at 1,000 samples**, and split peak resident into the part
that follows the psp block's reference window and the part that does not. The window is the handle —
sweep it at 5 kb, 20 kb and 80 kb and read off how much of the peak moves with it. Everything that
does not move is the floor, and the floor is what a store redesign cannot touch.

**Report one number at each sample count: the largest reduction in a cohort run's peak that a
perfect store could deliver.** Then compare it against the attribution already on record at 50
samples — per-sample chunks about 53% of live heap, the merger's projections about 30%, the
posterior fit about 8% — and say whether it agrees.

### Z1 — the 53% is the target, and a smaller file does not by itself shrink it

**The owner named that term as what this programme is for (2026-08-25), and it needs one distinction
made before any encoding is swept: those are decoded values held in memory, not compressed bytes on
disk.** The mass is the per-sample columns assembled for a whole cohort chunk and held for every
sample at once — at 50 samples, the read names alone were 1.66 GB, a third of the peak, as one
64-bit integer per read per position. **A file encoding that reconstructs exactly those lists saves
disk and decode time and changes that 1.66 GB not at all.**

Both specs say so deliberately:
[`psp_chain_id_encoding.md`](../spec/psp_chain_id_encoding.md) §2 makes "changing what the merge
consumes" a non-goal and §9 defers the in-memory shape elsewhere. **So the term the owner wants back
is, as things stand, out of scope of the very documents this plan was written to test.** This step
brings it in.

**⚠ And ng's own ruling makes that term bigger, not smaller.** Production stores names for about
3.4% of the reads it folds — it drops the reference allele's, which were 96.6% of them and are never
read, a change that alone cut a 50-sample run's peak by about a fifth. **ng names every read at
every position, by the 2026-08-17 ruling the cohort merge needs.** Held in the same shape,
that column is not a third of ng's peak; it is far more. **Measure it: fold ng-shaped records at
each rung of the ladder and report what the names cost in memory, held, against production's.**

So there are two ways to get the 53% back and this plan must price both, because only the second is
an encoding question:

- **Hold less per sample.** Do not assemble a whole chunk across all N samples before working on it.
  That is [`run_streaming.md`](../spec/run_streaming.md)'s architecture rather than a store format,
  and ng does not inherit production's chunk shape — so part of this saving may already be designed
  in and merely unmeasured. **Confirm that before crediting it to anything here.**
- **Hold the same information in a narrower form.** A name that is an index within a block rather
  than a genome-wide 64-bit integer; or keeping the names differential in memory and materialising a
  list only for the locus in hand.

  **⛔ Not first, and only if measurement forces it (the owner, 2026-08-25):** *"that would
  complicate the code quite a bit. So first let's try simpler approaches, let's evaluate them and
  let's increase the complexity only if needed."* Production already tried a block-local name and
  rejected it as considerably more complex and slower, for the file rather than for what is held.

  **So the order is fixed: measure what the simple things give first** — holding less per sample
  (above), and the plain encoding changes of Milestones B, C and D — and open this only if the
  measured gap to the objective still needs it. **Report the gap explicitly when that point is
  reached**, so the decision to add complexity is taken against a number rather than against a
  hunch.

**Report, per rung: bytes held per sample per position, for each candidate in-memory shape**, beside
the file figures the rest of the plan produces. A design that halves the file and doubles what is
held is a loss against this objective, and nothing else in this plan would notice.

*A caveat that has already caught this measurement once: heap profiling reports a larger live heap
than the process's resident size — 4.94 GB against about 3.0 GB on one 50-sample run — because
serialised allocation deepens the read-ahead and holds more chunks at once. Use the attribution for
proportions and the resident figure for totals; do not mix them in one sentence.*

## Milestone A — one instrument and one ladder

Nothing can be compared until every configuration is measured on the same records at a known
depth. This milestone is most of the work in the plan and none of the interest.

### A1 — the depth ladder, with honest depths

**HG002.** `benchmarks/ssr_hg002/bam/` already holds 5, 10, 15, 20, 30, 50 and 300 reads a
position, subsampled from the 300× source with a fixed seed so the ladders are nested — the 5×
reads are a subset of the 10× reads. Extend `subsample_coverages.sh`'s target list with **1, 3 and
100**. The script already computes fractions from the *measured* mean depth over the benchmark
regions rather than from the nominal label, so the rungs are honest; keep that.

**Tomato.** Build the same ladder from `benchmarks/tomato_big_cram/DRR000741.p1.cram`. Its depth
has to be measured first — the file is 49 GB and nobody has written down what depth that is. The
top rung is whatever the file gives; if it is below 100×, say so and let the human ladder carry the
deep end alone.

**Report per rung, in the results table and not only in a log:** measured mean depth over the
region set, number of covered positions, number of reads, mean read length.

**⚠ The HG002 ladder is reads restricted to 50,000 tandem-repeat intervals, not a whole genome.**
Reads sit in dense islands with empty ground between them, and that changes two things this plan
measures directly: how many reference positions a block of a given byte size spans, and therefore
how many blocks a file has and how big its index is. **The tomato ladder is the whole-genome shape
and must carry every index and block-count claim.** The human ladder carries depth.

### A2 — records as ng makes them, not as production wrote them

**The configuration sweep today re-encodes a production `.psp`, and production names about 3.4% of
the reads it folds while ng names all of them.** So every file-size number that sweep has produced
understates ng's at depth, and understates it worst exactly where the read names dominate. This is
the single largest source of error in the existing measurements and it has to be closed before any
of the three experiments means anything.

Source the records from ng's own generic locus generator instead. The walk already exists and is
already driven from a reference and an alignment file by two probes —
`examples/ng_generic_walk_probe.rs` (which retains nothing, and is the one to copy) and
`examples/ng_generic_loci_dump.rs` (which buffers every row and is the wrong shape for a memory
number). Feed the encoder from that walk.

**Done when** the record count and the per-position read counts from the encoder's input match
`ng_generic_loci_dump`'s output on one committed fixture, exactly.

### A3 — the sweep program reports per rung

One program, one grid, one table. Behind one interface: record source, field encoding, block size,
compressor reach, framing, read-name form. Every cell reports the six quantities above.

**Done when** one run over the whole ladder produces one table, and re-running it produces the same
table.

---

## Milestone B — how a block is cut, decoded, and how big it is

### B0 — how many separate streams a block is cut into, which the budget reopened

**With a budget of 500 kB rather than 30 kB per open sample, the choice of writing records
end-to-end instead of one buffer per field is no longer forced, and on the measurements we have it
is the more expensive choice.** The sweep behind
[`psp_record_encoding.md`](../spec/psp_record_encoding.md) measured both shapes at identical field
encodings, and at comparable reader memory the field-per-buffer shape is the smaller file:

| shape | reader memory | tomato, ≈3 reads a position | HG002, ≈30 |
|---|---|---|---|
| one buffer per field | 1024 kB | 6.27 bytes a record | 5.44 |
| one buffer per field | 256 kB | 6.60 | 5.77 |
| one buffer per field | 64 kB | 7.33 | 6.47 |
| records end to end | 132 kB | 7.32 | 6.48 |
| records end to end | 36 kB | 7.39 | 6.52 |

Read the middle rows against each other: **records end to end at 132 kB costs what one buffer per
field costs at 64 kB** — the same bytes for twice the memory — and at 256 kB the field-per-buffer
shape is about 10% smaller than anything the record stream reaches. At 30 kB per sample that
comparison was irrelevant, because the field-per-buffer shape cannot go there: dialled from 256 kB
down to 16 kB it grows by a third, where the record stream over the same range grows 3%. **At
500 kB per sample it is not irrelevant at all.**

So the axis to sweep is not two shapes but one number: **how many separately-compressed streams a
block is cut into.** One stream is records end to end; fourteen is a buffer per field; two is the
split that keeps the cheap per-position scan readable without inflating the rest. Grouping like
values together is worth having — at 512 kB on HG002, fourteen field buffers give 5.44 bytes a
record against 5.82 for one combined buffer, so about 7% of the file — and it is a knob, not a
shape.

**Sweep 1, 2, 4 and 14 streams**, at every block size and rung. Four is fields grouped by how they
behave — the position and the scan scalar; the three approximated quantities; the read names;
everything else.

**⚠ Streaming decompression multiplies by the stream count, and that is the whole tension.** A
reader must hold one decompressor per stream to assemble a record, so memory is roughly *streams ×
(reach + zstd's working buffers)*. At 150 kB a stream, fourteen streams is 2 MB a sample and four
times over budget; two streams is 300 kB and inside it. **So the budget does not simply favour the
field-per-buffer shape — it opens a middle that nobody has measured**, and finding where that curve
crosses the budget line is what this milestone is for.

### B0b — production's shape as an arm, and the three choices it bundles

**What production does is one point in this grid and it belongs in every table, at every rung, not
as a footnote at the end.** But it is not one decision — it is four made together, and comparing
against it as a lump would credit the new design with savings that a much smaller change to the old
one also delivers. The four, separated:

| choice | what production does | why it must be its own column |
|---|---|---|
| **field encoding** | fixed-width fields, floats at full precision | **this is most of the difference and it is not a layout question at all** |
| **streams per block** | one buffer per field, fourteen of them | B0's axis |
| **decode mode** | the whole block is decompressed before the first record | the thing streaming changes |
| **where a block is cut** | on a reference-coordinate grid, `DEFAULT_BLOCK_WINDOW_BP` = 5,000 bp ([`src/psp/writer.rs:119`](../../../../src/psp/writer.rs)), so blocks start at the same coordinate in every sample | **a property the new design does not have — see below** |

**The field encoding is worth more than the layout, and that has already been measured.** On tomato
at about three reads a position, production's own layout costs 11.10 bytes a record as it writes
them today and **6.27 with variable-length integers and the approximated floats and nothing else
changed** — 44% of the file, from a change to the fields rather than to the shape. On HG002 it is
10.17 against 5.44, 47%. **So the baseline that matters is production's layout with modern field
encoding, not production as shipped**; report both, and the difference between them is Milestone C's
result appearing in a different table.

**Blocks aligned across samples is the one thing production has that the new design does not, and it
was never priced.** Because every sample cuts at the same reference coordinate, a cohort reader
advancing over a segment touches one block per sample and they line up; a block cut by byte count
does not, so a segment can straddle two blocks in one sample and sit inside one in another. This is
also production's real memory lever — narrowing that window from 20 kb to 5 kb cut peak resident by
58% on a 50-sample run at four threads. **Sweep it as its own axis** — cut by bytes, and cut on a
reference grid at 5 kb, 20 kb and 100 kb — and report what alignment costs in bytes and what it buys
in the cohort walk of B4.

Two constraints it must satisfy to be admissible at all: a grid cut is a function of the coordinate
and the observation stream alone, so it keeps the property that the same sample gathered at any
worker count gives the same file; and a grid cut answers the "a reader may start every 100 kb of
reference" requirement directly, where a byte cut needs the extra rule of B3.

**Tools that already exist for this arm**, so none of it needs writing: `examples/blocksize_rewrite.rs`
rewrites a directory of `.psp` files at a chosen block target, `examples/psp_rechunk.rs` rewrites
them at a chosen reference-grid window, and `examples/psp_block_stats.rs` reports block count, span
and bytes per block from the index without decoding any payload.

### B1 — what a streaming decoder actually costs to hold

A reader that decompresses a whole 32 kB block holds 32 kB. A reader that streams through a 1 MB
block holds the compressor's reach — which we cap — *plus* the decompressor's internal working
buffers, and those are not something we choose: zstd's own internal unit runs to 128 kB.

**This is no longer a gate on the idea** — 150 kB a sample is inside the budget above — **but it is
the multiplier in B0's arithmetic and it decides how many streams we can afford.** Measure it two
ways and report both: what zstd itself predicts for a given file without decoding it, and the peak
resident of N readers actually open and streaming, at 1, 2, 4 and 14 streams.

### B2 — block size and compressor reach, swept independently

The existing sweep never separated them: every streaming configuration in it tied the compressor's
reach to the block, or used an 8 MB reach for both. The two knobs are already separate fields in
that program; the grid is what was missing.

Sweep block size **32 kB → 16 MB** against reach **8 kB → 1 MB**, at every rung, on both species.

**What to expect, and what would be a surprise.** Capping the reach means a large block cannot
compress better than a small one by finding more distant repeats — the only savings a large block
can produce are the ones it avoids paying at each restart: the framing, the running position and
coverage differences that reset, the live read set that must be restated, and the cold start a
small block suffers before the compressor has learned anything. If a large block wins by *more*
than those, something in the measurement is wrong.

**A consequence to price, not to assume away:** restart points get coarser, so serving a segment
that begins mid-block costs decoding from the block's start. Report the cost of a seek to a random
position at each block size — it is time, not memory, but it is the thing being traded.

### B3 — the block count, the index, and the coarse-restart rule

**At each block size, on the whole-genome tomato ladder, report the number of blocks in the file
and what an index over them would cost per open sample.** This is the measurement that decides
whether an index problem exists: at 32 kB blocks a whole-genome sample at three reads a position
has on the order of 500,000 blocks and an index around 5 MB, which is *worse* than the 3.8 MB
production index [`run_streaming.md`](../spec/run_streaming.md) §7.2 rejects. At 1 MB blocks it
should be tens of thousands and hundreds of kilobytes. **If large blocks bring it under the budget,
the coarse-index-and-chain design is not needed and must not be built.**

Two rules to check while here:

- **A block's size is in bytes, so its span in reference is whatever the coverage makes it.** At 1×
  on a sparse sample a block can cover a very long stretch, and
  [`psp_record_encoding.md`](../spec/psp_record_encoding.md) goal 3 asks that a reader be able to
  start no coarser than every 100 kb of reference. Measure the widest reference span any block
  reaches at each block size and each rung, and report whether the extra "cut at 100 kb" rule binds
  at all once blocks are large. At 300× it will never bind; at 1× it may be what sets the block
  count.
- **A large block may remove the need for the dictionary.** The dictionary — a few tens of
  kilobytes of representative bytes handed to the compressor before each block so a small one is
  not compressed from a cold start — exists only because 32 kB blocks are small. A 1 MB block warms
  its own context. Measure with and without: if a large block loses nothing without one, that also
  removes the trap that a dictionary is held *per open reader*, which at 112 kB across three
  thousand samples is 330 MB.

### B4 — memory against the number of open samples

The requirement is about a cohort, and the committed range reaches several thousand samples.
Measure peak resident with **1, 8, 64, 500 and 3,000 samples open and walked in lockstep**, at the
chosen block size, at a low rung and a high one. Files can be replicated with patched sample names
to reach the high counts; a precedent exists for that.

**Done when** the per-open-sample cost is a measured number at 3,000 open files rather than
arithmetic from one file, which is what every existing claim at that end is.

---

## Milestone C — the three approximated quantities

### C1 — the step sweep, per rung

Sweep each of the three one at a time and then together, at every rung, reporting each field's own
compressed bytes and the whole file's.

The steps to start from, and what each is for:

| quantity | what it feeds | proposed step | why that is thought to be safe |
|---|---|---|---|
| GC fraction of the window | a coverage-against-GC curve that bins its input | 1/100 | no consumer can tell 1% of GC from 0.01% |
| mean coverage of the window | a ratio of observed to expected depth | 1/4 of a read | no consumer can tell a quarter of a read from a sixteenth |
| summed log-error per allele | a likelihood, directly | 1/256 of a natural-log unit | the owner's ruling of 2026-08-25 — the **coarsest** allowed, not the proposal; see C2 |

**Report the share of the file each field holds at each rung.** This is the number that makes the
whole plan coherent: the reason to approximate them is expected to fade as depth rises, and how
fast it fades has never been measured. At about thirty reads a position the three together are 70%
of the compressed file. At 300× they are expected to be a small minority. **Nobody has measured
either end of that.**

### C2 — the summed log-error: the floor is ruled, the gain is not

**Ruling, 2026-08-25 (the owner): 1/256 of a natural-log unit is the coarsest step this field may
take.** That is a 0.4% error in a term that enters a likelihood directly. The step may be finer if
finer turns out to be nearly free; it may not be coarser, so the 1/16 and 1/4 steps that earlier
sweeps included are out of the sweep, not merely disfavoured.

What is left to measure is **how much the file grows on the way from 1/256 to full precision, at
each rung**: 1/256, 1/1024, 1/4096, then unapproximated. If the whole span is a few per cent at
every depth, the field stops being interesting and 1/256 ships because it is the cheapest allowed.

**The step interacts with depth, and this is the finding to watch for.** This is a *sum* of log
error probabilities over the reads supporting an allele, so its magnitude grows with depth — at
three reads a position the value needs about six bits and at three hundred it needs more than
sixteen. Stored as a variable-length integer count of a fixed step, a value three hundred times
larger costs more bytes, so **1/256 is not the same price at 300× as at 3×.** Report that field's
own compressed bytes at every rung.

**If it turns out expensive at the deep end, one question goes back to the owner and this plan does
not answer it:** whether the 0.4% tolerance is a fraction of the term (in which case the step may
grow with the value, and the cost stays flat across the ladder) or an absolute quantity of
natural-log units (in which case it may not). The two are the same thing at three reads a position
and are not at three hundred.

### C3 — round-trip, with the right strictness per field

**Three fields are lossy and every other field is not.** A round-trip check that compares whole
records with a tolerance will pass while the read-name list is being corrupted. Compare the integer
fields, the allele sequences and the read-name lists **exactly**; compare only these three against
their own step.

Also: a window that does not exist — an `N` reference position — is a real state and not a zero, so
one code is reserved for it and every present value shifted. A round-trip that never sees an `N`
region has not tested this; include one deliberately.

---

## Milestone D — the read names

### D1 — three forms across the ladder

Three ways to store which reads were present, measured at every rung on both species:

- **the whole list, raw** — a count then one 8-byte identifier per live read at every position,
  which is what a straight port of production's column gives;
- **the whole list, as distances** — the same list, each identifier stored as its distance from the
  one before, as a variable-length integer;
- **the changes only** — who arrived and who left, arrivals as distances and departures as their
  place in the live set.

Between 11.4 and 293 reads a position these cost 1.02 → 43.78, 0.67 → 11.72 and 0.43 → 6.42 bytes a
position. **The ladder is what turns those four points into two curves**, and the growth exponent
of each is what decides the design. Note also what those numbers already say: the distances form,
which is a few lines and needs no reader state at all, captures 86% of the available saving at the
deep end.

*Those figures were taken on benchmark slices — tomato at 11.4 reads a position and HG002 at 293 —
and the deep one comes from reads concentrated into 1,000 small regions rather than from a 300×
library. A1's ladder is what replaces them.*

### D2 — a read covers two stretches, not one

**A name is allocated per read *pair*, and the two mates rarely overlap, so most names go live,
stop, and go live again.** Measured on the two corners: 83% of names on HG002 and 91% on tomato
cover two stretches with a gap between them. **A changes-only stream that assumes one stretch per
name loses the second mate of nine reads in ten — and loses it silently**, because the cohort merge
would simply see a read that was not there.

So a re-entry form is part of the first version, not a later fix, and this step's job is to prove
it: the reconstructed live set at every position must equal the one the walk produced, exactly, at
every rung. **Re-measure the two-stretch fraction across the ladder** — subsampling changes which
mates survive, so it is not a constant.

### D3 — what restating the live set costs against block size

**This is the step the block-size question opened, and it is why Milestones B and D are not
independent.** Every restart point has to restate the whole live set, and at 300× that set is 300
names rather than 11.

What is known: on tomato at 11.4 reads a position, cutting a block every 1,500 positions costs
0.432 bytes a position against 0.385 when cut by byte count — so the entire restatement overhead
there is 0.047 bytes a position, 12% of the read-name column and under 1% of the file. **The deep
end is unmeasured and is where it should matter**, because the set being restated is twenty-six
times bigger.

Sweep the block size against the rung and report the restatement's share of the read-name column at
each. The program already takes both a byte target and a positions target for the cut.

**One constraint the cut rule may not break:** where a block ends must be a function of the
observation stream alone, never of how the writer was scheduled, or the same sample gathered at
different worker counts gives different files —
[`run_streaming.md`](../spec/run_streaming.md) §12.1. A byte target counted over records as they
arrive satisfies this; a flush driven by a queue depth or a timer does not.

**Two classes of read are counted and never named, and must not leak in:** reads that produced no
observation, and reads a depth cap discarded. A read the cap dropped is in no observation, so if it
is in the live set the derived reference list gains a read nobody folded. The check that catches it:
a derived list's length must be at most the observation's read count and at least half of it, since
at most two mates share a name.

---

## Milestone E — the combination, and what goes to the owner

### E1 — one configuration, measured end to end

Take the winner of each of B, C and D, run them **together** over the whole ladder on both species,
and report the six quantities. **The savings do not add** — the point of the grid — so this is a
measurement and not a sum.

### E1b — what the memory is being paid for in time, and how much of it cores buy back

**Every choice in this plan may be trading processor time for memory, and cores can pay some of
that back but not all of it.** Two quantities, and conflating them is the failure to avoid:

- **CPU-seconds per million records decoded.** This is what a cohort run actually spends and no
  number of cores reduces it. Report it for every configuration on the ladder.
- **Wall time at 1, 2, 4, 8 and 16 threads.** This is what an operator feels, and it is what cores
  buy back.

**Decoding should parallelise almost perfectly and the experiment has to check that rather than
assume it**, because it is the argument for accepting a slower configuration. Two independent
levers exist: every open sample is a separate file, and every block is independently decodable
within one file. So a cohort run has thousands of independent units of work and should scale until
something else is the bottleneck.

**What to watch for, and it is the reason this step exists.** Independent blocks scale; shared
state does not. A dictionary shared across readers, a single input buffer pool, or a merge that
must consume samples in lockstep will each put a floor under the scaling curve, and the curve is
the only thing that shows it. Report the curve, not a speedup figure at one thread count.

**Report a wall-time-against-memory frontier at the end**, at a low rung and a high one: for each
configuration, its per-open-sample memory and its wall time at the thread count the machine has.
That frontier, with the file sizes beside it, is what the final choice is made on.

### E2 — against production's columnar blocks, on the same ladder

**Production's shape is a first-class arm swept across the whole grid (B0b), and this step is where
it is reported against the winner rather than inside the sweep.** Four rows, on the same records at
every rung:

1. **production as it writes files today** — fixed-width fields, fourteen buffers, whole-block
   decode, a 1 MiB target cut on a 5 kb reference grid;
2. **production's layout with the field encoding of Milestone C** — the honest layout baseline, and
   the one the design argument actually turns on;
3. **production's layout at its best block size for the budget**, found with `blocksize_rewrite.rs`;
4. **the recommended configuration** from E1.

Rows 1 and 2 differ by 44% of the file on tomato and 47% on HG002 at the depths measured so far.
**If most of the remaining gap between rows 2 and 4 is small, the honest recommendation may be to
keep the columnar blocks and change the fields** — and this plan has to be able to produce that
answer, or it is not an experiment.

**⚠ Two things make the comparison not like-for-like, and both have to be reported apart from the
totals.** Production names 3.4% of the reads it folds and ng names all of them, so ng's records are
bigger for a reason that has nothing to do with encoding: report every row a second time with the
read-name column excluded, so *"ng stores more"* and *"ng stores it worse"* are separable. And
production's records are not ng's — no partial read witness, no read group in the observation
identity — so a per-record byte figure compares two different records. **Report bytes per covered
position alongside**, which is the unit both shapes share.

### E3 — the recommendation

One configuration, with: the curve of bytes per record over 1× to 300× on two species, the measured
per-open-sample memory at 3,000 open files, the wall-time-against-memory frontier from E1b, and the
growth exponent of each candidate read-name form.

**And one number the owner should be given rather than asked for: where the 500 kB working budget
sits on the bytes-against-memory curve.** If halving it costs 2% of the file the budget is not
binding and should be halved; if doubling it saves 10% the budget is the wrong place to be. That is
a reading off a curve this plan produces, not a separate study.

---

## Traps — what will bite whoever runs this

- **A dictionary trained on the blocks it is then measured against reports a saving no reader ever
  gets.** Easy to do by accident and the number is spectacular: an early run of the existing probe
  reported a hundredfold. Train on one half of the blocks, measure on the other.
- **Heap profiling overwrites the plain binary.** A `dhat` run without `--target-dir` replaces
  `target/release/examples/<probe>`, and the next plain run silently re-executes the instrumented
  build at five to six times slower with nothing in the output saying so. This has already
  destroyed one performance review's measurement set.
- **Two build trees hold different binaries.** The container build points its target directory at
  `target-container/`; a direct `cargo` run on a machine with no container runtime writes
  `target/`. A script that looks for a built binary must check both and take the newer.
- **Every running difference resets at a block boundary.** The position, the coverage difference and
  the read-name base all restart, and a block that does not reset one reads back wrong from its
  first record — *plausibly* wrong, because coverage is smooth. The test that catches it: reading
  from an arbitrary block gives exactly what a full sequential read gives from that point.
- **Compare after compression, never before.** The changes-only form looks far better than it is on
  raw bytes, because the whole-list form is highly compressible already: 679 MB of raw identifiers
  on tomato became 7.4 MB, a ninety-two-fold collapse. At 293 reads a position the same collapse
  only reaches thirty-six-fold, and that gap is the whole finding.
- **The depth label is not the depth at a position.** A "300×" run is a mean; across a benchmark
  fixture a sample's reads at a locus have been seen to run from 8 to 428. Report the distribution
  at each rung, not only its mean, or a claim about behaviour "at 300×" is a claim about an average
  that few positions have.
- **Scratch goes in this repository's `tmp/`**, never in the system temp directory or the
  assistant's own scratch directory: a path outside the project mount is invisible inside the dev
  container, so nothing written there can be handed to `cargo`.

---

## Settled on 2026-08-25, and not to be re-derived

- **Every read is named at every position it covers.** A requirement of the cohort merge, not an
  encoding preference.
- **500 kB resident per open sample is affordable**, so the memory constraint is not what chooses
  the design. This reopened B0.
- **The summed log-error's step is 1/4,096 of a natural-log unit** (the owner, 2026-08-25),
  superseding the earlier "1/256 is the coarsest tolerable". Measured: once the field is an integer
  at all it is worth 16 % of the file at three reads a position and 21 % at 279, and which step is
  chosen barely matters — 1/4,096 costs 5 % against 1/1,024 and buys a sixteenth of 1/256's error.
- **⚠ The rounding happens in the type, not in the psp writer, and that is ng-wide work this plan
  does not own.** Approximating in the psp alone would make direct mode and psp mode see different
  numbers and break the oracle the whole psp path is checked against
  ([`../spec/run_streaming.md`](../spec/run_streaming.md) §1.2). Rounding where the value is computed
  makes both routes agree — and agree *better* than full precision would, since it absorbs the
  last-bit differences of two summation orders. **Carried here so it reaches the psp implementation plan, which is
  where the owner put it on 2026-08-25**; the same treatment applies to the window's GC fraction and its mean
  coverage.
- **The per-sample summary does not go in ng's psp as TOML text.** Settled after Milestone Z
  measured it at 1.05 MB of text per open sample for a 0.52 MB parsed histogram, both resident at
  once. The section is permanent — ng's duplication filter is not built yet but is coming — so it is
  encoded as binary and the reader does not keep the encoded bytes after decoding.
  [`../spec/psp_record_encoding.md`](../spec/psp_record_encoding.md) §4.1 owns it.
- **Wall time may be traded for memory, within reason, because cores can pay some of it back** —
  which is why E1b measures processor time and the thread-scaling curve separately, rather than
  reporting one wall-time figure.

## Open questions this plan does not close

1. **Is the summed log-error's tolerance a fraction of the term, or an absolute number of
   natural-log units?** — the owner's, and only if C2 finds the fixed step expensive at 300×. The
   two are indistinguishable at three reads a position.
2. **Are the separately-compressed streams within a block interleaved or held in two regions of the
   file?** — unmeasured and untouched here. It is a seek-pattern question that only a real segment-serving
   reader can answer, and no writer exists yet.
3. **Multi-library samples.** Every rung in this plan is one read group, and the read group joins an
   observation's identity, so a sample with several has more and smaller observations per record.
   Expected somewhat worse, not structurally different. Left to the first run on a multi-library
   cohort.

## Where the existing numbers came from

Every figure quoted above without a new measurement behind it comes from the sweeps in
[`../research/per_sample_record_store_compression_2026-08-19.md`](../research/per_sample_record_store_compression_2026-08-19.md),
taken on a tomato accession at about three reads a position, GIAB HG002 at about thirty, and two
benchmark slices at 11.4 and 293. **One sample and one cohort of 63 is not the committed range**,
which is the reason for this plan.
