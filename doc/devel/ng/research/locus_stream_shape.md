# How the calling pipeline should pass loci from one stage to the next

**Date:** 2026-08-06. **Status:** research finding, settled by measurement; one decision
deliberately left open.

---

## What this document answers

The pileup produces one locus for every covered base of the genome. Five stages downstream
consume those loci. This document answers a question that had to be settled before any of them
is written, because it is expensive to change later:

> Should a locus travel through the pipeline as **one self-contained object per locus**, or
> should many loci travel together as **one array per field**?

The answer is that a locus travels as a self-contained object everywhere in memory, and the
arrays are used only inside the file. This document explains why, because the reasoning is not
obvious and the opposite arrangement is what most people expect.

Four experiments settled it. They built both arrangements at every stage, ran them on real
tomato and human alignments, and checked that both produced identical output before comparing
their cost. The numbers are in §7; the reasoning is in §4 and §5.

---

## 1. Vocabulary

A few software terms carry the whole argument, so they are defined here rather than assumed.

**A locus.** One covered base of the reference genome, together with what the reads said about
it: which sequences were seen, how many reads showed each, their base qualities and mapping
qualities. The pileup emits one per covered base — about 2.9 billion for a human sample at 30×
coverage.

**Record layout.** Each locus is a self-contained object. It owns its own reference sequence,
its own list of observations, and inside each observation its own sequence of bases. This is
what the pileup produces today.

**Columnar layout.** Many loci are held together, with one array per field: all the depths in
one array, all the quality sums in another, all the allele sequences end to end in a third.
Because different loci have different numbers of observations, the arrays that hold them need a
companion array of positions saying where each locus's part begins and ends. Fields shaped like
that — variable length per locus — are the awkward ones, and they matter later in the argument.

**A block.** One batch of loci in columnar layout — the unit that is filled, handed over, or
written to disk together. The `.psp` file is a sequence of blocks.

**An allocation.** A request to the operating system's memory manager for a piece of memory,
and later a matching release. Each one costs roughly 340 machine instructions on this project's
hardware, measured. A locus in record layout needs about thirteen of them: one for its reference
sequence, one for its list of observations, and one for each observation's sequence and chain
identifiers.

**Cache.** A small, very fast memory close to the processor. This machine's fastest level holds
128 kilobytes per fast core. Data that fits there is reached in a few cycles; data that does not
costs hundreds of cycles to fetch from main memory. Much of the usual argument for columnar
layouts is about this distinction.

**Instructions retired.** The count of machine instructions a program actually executed. This
project measures with it rather than with elapsed time, because elapsed time on this machine
varies by more than the effects being measured, while the instruction count is reproducible to
one part in ten thousand.

---

## 2. The recommended architecture

**Loci are self-contained objects everywhere in memory. Arrays of fields are used only inside
the file.**

Stage by stage, and what each hands the next:

1. **The pileup emits one locus at a time**, as a self-contained object, exactly as it does
   today.
2. **The parameter pre-pass takes them one at a time.** It builds its histograms, estimates its
   parameters, and copies out the hundred thousand loci it keeps for later. It never holds a
   reference into anything the pileup owns.
3. **The file stores arrays of fields**, in blocks of about a megabyte. This is the one place
   the columnar layout is used, and it is where it pays — see §5.
4. **The merge reads each sample's file and scans one small number per locus first**: how much
   non-reference evidence that sample carries there. It builds the full locus only where some
   sample might have a variant — about one position in a hundred. For those positions it works
   with self-contained objects, one locus of lookahead per sample.
5. **The calling step takes self-contained objects.**
6. **The VCF writer takes self-contained objects.**

Only one of these six stages sees arrays of fields, and it is the one that writes and reads a
file.

---

## 3. Why the opposite arrangement is the natural expectation

Columnar layouts have a real and well-known advantage, and it is worth stating fairly before
explaining why it does not apply here.

They come from analytical databases, where the same very large table is swept many times, and
each sweep reads two or three fields out of forty. Storing each field in its own array means a
sweep touches only the fields it needs. Everything else stays on disk or out of cache. The
saving is enormous and it is why every analytical data format is built this way.

The expectation that this pipeline should work the same way is therefore reasonable. It is also
wrong, for four separate reasons that happen to point the same direction.

---

## 4. Why arrays of fields do not pay between the stages

### 4.1 The data is swept once, not many times

Each locus is produced once, read once to build histograms, written once, read once by the
merge, and consumed once by the calling step. That is five touches, and each stage wants a
different part of the locus.

Only two of the five skip anything worth skipping. The pre-pass reads three or four of the ten
or so fields. The merge's first pass reads exactly one. **Those two are where the
recommendation keeps the columnar layout** — inside the file, which both of them read from.
Every other stage wants nearly the whole locus, so there is nothing to skip and nothing for the
layout to win.

### 4.2 Nothing is waiting on memory

The other classical advantage of columnar storage is that it brings in only the bytes a stage
needs, instead of whole cache lines full of fields it will not read. That advantage only exists
when the data does not fit in cache.

It fits here, with room to spare. Everything one covered base needs — the active reads, their
sequences and qualities, the accumulating locus — is **10 kilobytes at 30× coverage against a
128-kilobyte cache**. Across 2.5 million covered bases of a human chromosome, **not one exceeded
it**. The data is already as close to the processor as it can be put. Rearranging it cannot
bring it closer.

### 4.3 The saving inside the pileup is about allocation, not arrangement

Filling arrays inside the pileup does save real work: **1,344 instructions per covered base**,
which is 7.8 % of what a covered base costs at 30× coverage. That saving is genuine and it was
measured four ways.

But the mechanism is that thirteen small allocations per covered base stop happening — not that
the data ends up arranged by field. A pool of reusable objects, or one large buffer that many
loci are carved out of, would recover most of the same saving with no change of layout at all.
The columnar arrangement is one way to stop allocating. It is not the source of the benefit, and
attributing the benefit to the arrangement is what makes the arrangement look more valuable than
it is.

### 4.4 Every boundary between the two shapes costs a conversion

This is the reason that reverses the intuition rather than merely weakening it.

Variable-length fields are the difficulty. A locus has a variable number of observations, and
each observation has a variable-length sequence. You cannot take a range of loci out of a block
without rewriting every position marker after it. So a stage that receives arrays and needs
objects pays to build them, and a stage that receives arrays and needs a *different* range of
arrays pays to copy them.

This project has already paid for that twice, in the production caller:

- A per-locus object sitting between two columnar stages was identified as pure overhead and
  removed.
- What replaced it — slicing ranges of loci out of one set of arrays into another — then became
  **the largest single consumer of memory in the whole cohort path: 49 % of peak, with 20.8
  gigabytes moved.** That is more than double what the objects it replaced had cost.

The same pattern appeared again in the new experiments, twice. At the calling boundary, the
columnar version ended up **one allocation per locus** away from the object version, because the
VCF writer needs an object anyway and the data had to be copied back out. At the merge, building
objects out of a block cost **nearly four times** what receiving objects from the pileup had
cost.

**The layout moves the copying. It does not remove it.**

### 4.5 Where the real work is, the layout is irrelevant

The calling step spends **999,224 instructions per called locus** on its arithmetic. Building
that locus as a self-contained object costs **5,999** — the arithmetic is 167 times the data
handling. Whether the calling step is handed an object or a set of array slices changes its
total cost by **0.10 %**, while the object version is half the source code.

Nothing about the calculation resisted the columnar form. An entry point taking array slices
already existed and was used, so both versions ran identical arithmetic. The columnar version
simply had nothing to win.

---

## 5. Where arrays of fields do pay, and why

Two places, and they share a mechanism: both **skip** data rather than touch it.

**In the file.** Values of the same kind sit next to each other, which is what lets a
general-purpose compressor work well on them. A reader that does not need a field never
decompresses it. And the block written to disk is close to the block held in memory, so writing
is close to a copy per field. None of this was ever in question, and the file format already
works this way.

**In the merge's first pass.** Reading one small number per locus across every sample, and
building the full locus only where some sample might carry a variant, is the single largest
saving found anywhere outside the pileup: **the merge does 2.0 times less work, makes 414 times
fewer allocations, and moves 252 megabytes instead of 4,638.**

The important detail is what makes that possible. It is not that the whole locus is laid out by
field. It is that **one summary number per locus is available without touching anything else**.
That is a much weaker requirement, and it is why the merge in §2 stays object-shaped: it reads
that one number, then works with self-contained objects for the one position in a hundred that
survives.

This was tested directly. A merge that reads the summary number but is otherwise written the
conventional way, with one locus of lookahead per sample, captures **90 %** of the benefit of a
fully columnar merge — without any of the machinery a columnar merge needs to keep its
references valid.

---

## 6. The rule this leaves behind

**The cost is in how many objects get built, not in how the data is laid out.**

Every saving that mattered came from *not building* something:

- not building thirteen small pieces of memory per covered base, inside the pileup;
- not building 99 loci in 100 at the merge, because a summary number said no sample could have a
  variant there.

Every attempt to win by changing the arrangement alone moved a copy somewhere else instead. That
rule decides the next case without another round of experiments: when a stage is expensive, ask
what it is building that nobody needs, before asking how its data is arranged.

---

## 7. What was measured

Explanation ends here; this section is the record.

**Instrument.** Instructions retired, taken from `/usr/bin/time -l`, with start-up cost
subtracted by differencing two runs of the same input at two different locus counts. Minimum of
three to five runs per arrangement, with the binaries alternated inside one script. Elapsed time
was never used and no sampling profile was used as the source of any number: other work was
loading the machine, and both of those measurements are distorted by it while instruction counts
are not. Peak memory came from the same command; allocation counts from a heap profiler.

**Fixtures.** Real alignments throughout, never synthetic data. Both of this project's largest
past savings came from skew in real data — about 4 positions in 100 are variable, and 96.6 % of
the haplotype identifiers stored were for the reference allele and were discarded unread —
and synthetic data has neither property, so it would have reported that the columnar layout buys
nothing.

- One tomato sample, 49 GB, about 130× coverage.
- One human sample at 10×, 30× and 300× coverage.
- One sparse tomato sample at about 5× coverage.
- Ten tomato samples walked in step, for the merge.
- Fifty tomato samples as built `.psp` files, for the merge and the calling step.

**Correctness.** Every arrangement had to produce identical output before its cost was compared.
The merge and the calling step matched to the exact bit, with no tolerance allowed — 0 of 23,552
loci differed in genotype, quality or allele frequency. The pileup arrangements matched on a
digest covering every field of every locus.

**Headline results.**

| what was compared | result |
|---|---|
| filling arrays inside the pileup, against one object per locus | −1,344 instructions per covered base: **−7.8 % at 30×**, −12.1 % at 5×, −2.0 % at 130×; 45 % fewer allocations |
| the pre-pass reading arrays, against reading one locus at a time | **1.5 % apart** — the consumer's shape does not matter |
| the merge reading a file, scanning a summary number first | **2.0× less work**, 414× fewer allocations, 252 MB moved instead of 4,638 MB |
| the same, but the merge otherwise written conventionally | captures **90 %** of that benefit |
| the merge reading the pileup directly, arrays against objects | **indistinguishable**; the merge is 1.6 % of that path |
| the calling step, array slices against one object | **0.10 %**, and the object version is half the source |
| peak memory, columnar against objects | columnar **0.2 % to 4.8 % higher**, never lower |

**Two findings unrelated to this question, both worth acting on.** The largest single source of
allocations in the whole pileup is **copying each read's name** — 224,000 of 1.5 million, and
nothing in this study touches it. And the calling step allocates 13.85 pieces of memory per
locus, 1.8 times what building the locus costs, because a shape check builds and sorts a nested
list purely to compare its length and then discards it; the table it needs is already cached
elsewhere.

---

## 8. The decision left open

**The pileup could fill a block internally and still hand out one locus at a time**, invisible
to every stage downstream. That is worth about 8 % of what a covered base costs at 30×
coverage.

It is not recommended yet, for three reasons:

- It costs **830 lines** in the three files that carry the walk's hottest code, which have just
  absorbed three rounds of careful optimisation.
- Two of its new invariants are the kind that are either correct or silently wrong, with nothing
  in between, and **one of them cannot be caught by any test this project has**, because getting
  it wrong costs memory rather than changing any answer.
- Running the pileup on more than one core has never been attempted and is worth several times
  more.

**If it is taken later, take it in this form**: blocks inside the pileup, a four-byte summary
number per locus stored alongside them, and every stage downstream still receiving one
self-contained locus at a time. That captures nine tenths of the available benefit and leaves
every consumer in the shape this document recommends. The summary number is not optional — a
block without it makes the merge **worse** than today, because building objects out of arrays
costs more than receiving objects that were never taken apart.

---

## 9. Where the detail lives

The four experiment reports, with full tables and the code as diffs, are in
[experiments/locus_stream_shape/](experiments/locus_stream_shape/); its
[README](experiments/README.md) says what each one measured. The plan they were run under, and
a condensed results section, is
[locus_stream_shape_experiments.md](../impl_plan/locus_stream_shape_experiments.md).

Raw profiles and heap dumps were not kept. They were hundreds of megabytes and are reproducible
from the commands each report records.
