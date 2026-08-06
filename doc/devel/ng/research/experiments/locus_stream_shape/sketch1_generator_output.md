# Sketch 1 — what the generator emits, and what the parameter pre-pass reads

**Date:** 2026-08-06. **Base commit:** `1e5ffa8`. **Status:** throwaway sketch, to be deleted
after the decision. Plan: `doc/devel/ng/impl_plan/locus_stream_shape_experiments.md` §4.

---

## The answer, in one sentence

Emitting blocks of parallel arrays instead of one owned locus per covered base removes
**1,344 instructions of the 4,533 a covered base costs regardless of its depth** — about a
third — which is **7.8 % of a covered base at 30× and 2.0 % at 130×**; the record-shaped
consumer that refills its own buffer (arm C) gives back a fifth of that and lands **1.5 %
behind** the columnar consumer, so **no consumer ever needs to see columns**, and the real
question is whether 8 % at the target depth is worth **830 lines of producer** that a
record-shaped pipeline does not need.

---

## Vocabulary, defined once

- A **covered base** is one reference position the walk turns over and delivers a locus for.
  On the generic path a locus is one covered base 992 times in 1,000 (cited, `price.md`).
- A **read at a position** is one `(covered base × read)` pair the fold visits — what the
  2026-08-05 price table calls a *contributor visit*. Measured here as `Σ num_obs` over
  every emitted observation, which counts the same pairs.
- **Depth** is reads-at-a-position per covered base, which is what the walk pays for, not
  the sample's nominal coverage.
- A **record** is today's `SampleLocusObservations`: one owned object per covered base,
  carrying a boxed slice of reference bases, a `Vec` of observations, and per observation a
  boxed slice of allele bytes and a `Vec` of chain ids.
- A **block** is many loci with each field held as one array — a flat byte buffer plus an
  offsets array for each ragged field, and an index saying where each locus's observations
  begin.

The four arms:

| arm | the walk emits | the pre-pass reads |
|---|---|---|
| **A** | one owned record per covered base — **today's shipped generator** | the record, record-shaped |
| **B** | a block the caller owns and reuses; no locus is ever materialised | the block's columns directly |
| **C** | the same block | one locus at a time, through a view over a buffer **the pre-pass owns** and refills per locus |
| **Cb** | the same block | the same record-shaped code, but over slices **borrowed straight out of the block** |

**Cb is not a candidate.** It is the shape production built and made test-only
(`BlockColumns<'a>`, `src/psp/reader.rs`), because a consumer holding a locus across calls is
self-referential. It is measured here for one purpose: `C − Cb` is exactly what refilling the
consumer's own buffer costs, and without it that cost is a guess.

---

## §7's table

Instructions and peak resident memory are **measured in this worktree**, on the three
fixtures below, all in whole-contig mode (`PVC_PROBE_WHOLE_CONTIG=1`). Per-base figures come
from differencing two locus counts of one fixture, which cancels start-up and reference
loading exactly. Every figure is the **minimum of three runs**, with the four arms alternated
inside one script.

| arm | instructions per covered base, 30× BAM | per read at a position, BAM | fitted depth-independent term | peak RSS, 30× / 130× | allocations per base, 30× | bytes copied at the handoff, per base | lines of source | what it cannot do |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| **A** records | **19,757** | 663 | **6,035** | 149.2 MB / 1,069 MB | **7.72** | 0 copied, 7.72 allocated | **61** consumer, **0** producer | cannot stop allocating; hands over one locus at a time and nothing else |
| **B** columns | **18,216** (−7.8 %) | 654 | **4,691** (−1,344) | 149.8 MB / 1,071 MB | **4.25** (−45 %) | **83** | **48** consumer, **830** producer | every stage downstream must be written columnar, including the ones nobody knows how to write that way |
| **C** scratch view | **18,482** (−6.5 %) | 653 | **4,979** (−1,056) | 149.8 MB / 1,073 MB | **4.25** | **130** (83 + 47 refill) | **131** consumer, **830** producer | the view is a different type from `SampleLocusObservations`, so downstream code is re-typed even though it is not re-shaped |
| **Cb** borrowed view | **18,199** (−7.9 %) | 654 | **4,665** (−1,370) | 149.7 MB / 1,073 MB | **4.25** | **83** | **66** consumer, **830** producer | a consumer cannot hold a locus across `next()` — self-referential, and production already retired this shape |

Notes on each column, so none of it is read as more than it is.

- **Instructions per covered base** — measured, HG002 30× BAM `chr1`, whole-contig mode,
  770,000 → 1,540,000 loci. Whole-contig, not the shipped region stream, because that fixture
  is tandem-repeat-targeted: through the region stream its `chr1` types 165 reference bases
  for every covered base it delivers (cited), and half the run is region setup rather than
  walking.
- **Per read at a position** — derived from the fitted depth-independent term, exactly as the
  2026-08-05 price table derives it. Arm A's **663** against the shipped baseline's cited
  **659** is the check that arm A is still the shipped walk: 0.6 %.
- **Fitted depth-independent term** — from two tomato CRAMs on the same reference and the
  same contig at depths **4.87** and **97.97**, solving `per_base + per_visit × depth` for
  each arm. Arm A's 6,035 is higher than the shipped 4,533 (cited) because arm A carries the
  pre-pass; the arm-to-arm **difference** is what this sketch is about.
- **Peak RSS** — measured, `/usr/bin/time -l`, at the upper endpoint of each fixture. The
  columnar arms are 0.2–0.4 % **higher**, not lower.
- **Allocations per base** — measured with `--features dhat-heap --target-dir target-dhat`,
  200,000-locus prefix. Whole-process, so it includes the read decoder, which is why arm A
  reads 7.72 and not the thirteen a locus by itself carries.
- **Bytes copied at the handoff** — measured. Arm A copies nothing: the walk hands over an
  object it built and the consumer takes ownership. Arms B/C/Cb write the whole payload into
  the block once (83 bytes per covered base at 30×, from 490 blocks × 256 KiB over 1,541,702
  loci); arm C copies it a second time into its own buffer (47 bytes per base, counted
  directly by the sketch).
- **Lines of source** — code lines only: blank lines and comments stripped. Consumer counts
  are the pre-pass's own reduction plus its buffer machinery. The producer count is the whole
  cost of teaching the walk to fill a block: a new 367-line module plus 516 lines added and
  53 removed across five existing files.

### The same arms at three depths

Instructions per covered base, measured, minimum of three:

| fixture | depth | A | B | C | Cb | A → B |
|---|---:|---:|---:|---:|---:|---:|
| tomato big CRAM, `SL4.0ch01` | 97.97 | 114,648 | 112,510 | 112,248 | 112,608 | **−2.0 %** |
| HG002 30× BAM, `chr1` | 20.69 | 19,757 | 18,216 | 18,482 | 18,199 | **−7.8 %** |
| tomato sparse CRAM, `SL4.0ch01` | 4.87 | 11,436 | 10,054 | 10,316 | 10,035 | **−12.1 %** |

**The saving is almost all a fixed number of instructions per covered base, not a fraction of
one.** In absolute terms A → B is 2,138 instructions per covered base at depth 98, 1,541 at
depth 20.7 and 1,382 at depth 4.9, and the fitted split says why: **1,344 of it is the
per-base term and about 8 more come off each read at a position**. So what changes with depth
is mostly the denominator, not the saving. At 130× a covered base costs 114,648 instructions
and the saving is a fiftieth of it; at 4.9× it costs 11,436 and the saving is an eighth.
**The stated target is one human sample at 30× whole-genome, where it is 7.8 %.**

The per-read term barely moves — 1,108 → 1,101 on a CRAM, 663 → 654 on a BAM, under 1.5 % —
which is what should happen: the change is to what a *column* costs, not to what a *read*
costs.

### What that is worth over a genome

At the target — one human sample, 30×, about 3.1 × 10⁹ covered bases — the measured saving of
1,541 instructions per covered base is **4.8 × 10¹² instructions**, against the cited
1.4 × 10¹³ of column overhead the plan identified and about 5.7 × 10¹³ for the whole walk.

---

## What the pre-pass itself costs, which sizes the whole question

The shipped probe (`ng_generic_walk_probe`) walks the identical pipeline and drops every
locus, so arm A minus the probe is the pre-pass stand-in's own price. Measured, same
endpoints, minimum of three:

| fixture | bare walk, per covered base | arm A | the pre-pass costs | arm B |
|---|---:|---:|---:|---:|
| HG002 30× BAM `chr1` | 18,406 | 19,757 | **1,351** (7.3 % of arm A) | **18,216** |
| tomato sparse CRAM | 10,052 | 11,436 | **1,384** (12.1 % of arm A) | **10,054** |

**The columnar walk with a full pre-pass on top costs what today's walk costs with no
pre-pass at all.** On the sparse tomato that is 10,054 against 10,052 — two instructions in
ten thousand. That is two independent numbers of similar size landing on top of each other,
not a law, but it is the cleanest statement of the scale: **what the block layout removes
from the producer is about what the entire parameter pre-pass adds.**

It also answers the second half of the plan's question directly. The pre-pass is **7 % of a
covered base at 30×**. Arguing about its shape is arguing about 7 %, and the arms differ
inside that by 1.5 points.

---

## The block-size sweep

Block size is the known memory lever in this project — production measured 16 MiB blocks
against 1 MiB at **2,501 MB against 261 MB** of peak resident memory with wall time flat
(cited). **It does not behave that way here.**

HG002 30× `chr1`, arm B, whole-contig, minimum of three, normalised by each run's own locus
count (the budget is checked per block, so a larger block overshoots the locus cap):

| block budget | blocks | instructions per covered base | peak RSS |
|---:|---:|---:|---:|
| 4 KiB | 31,027 | 19,303 | 149.1 MB |
| 16 KiB | 7,814 | 19,300 | 149.0 MB |
| 64 KiB | 1,957 | 19,294 | 148.9 MB |
| 256 KiB | 490 | 19,295 | 149.5 MB |
| 1 MiB | 123 | 19,265 | 151.2 MB |
| 4 MiB | 31 | 19,250 | 156.9 MB |
| 16 MiB | 8 | 19,229 | 176.2 MB |

Tomato big CRAM `SL4.0ch01`, arm B, 1 M loci:

| block budget | blocks | instructions per covered base | peak RSS |
|---:|---:|---:|---:|
| 16 KiB | 5,255 | 114,335 | 725.9 MB |
| 256 KiB | 330 | 115,113 | 733.0 MB |
| 4 MiB | 21 | 114,304 | 735.4 MB |
| 16 MiB | 6 | 113,485 | 787.4 MB |

Two readings, and both are findings.

**Instructions are flat across a four-thousand-fold range of block size** — 19,303 down to
19,229 on the BAM, a 0.4 % decline, and within 1.4 % on the CRAM with no trend. There is no
knee. Whatever a bigger block saves in loop overhead it gives back in cache misses.

**Peak resident memory grows by roughly twice the block budget, and by nothing else.** A
16 MiB block costs 27 MB over a 4 KiB one on the BAM and 61 MB over a 16 KiB one on the CRAM
— the payload plus the capacity slack of twenty-two separately-growing arrays. Against the
walk's own 149 MB and 726 MB that is +18 % and +8 %.

**Why this differs from production's 2,501 MB → 261 MB.** There the block *was* the peak: the
`.psp` decoder held whole decoded blocks and nothing comparable beside them. Here the walk's
peak is the reference window, the read decoder and the active read set, and the block is a
rounding error beside them. **Block size is a real memory knob and a small one, and it is not
an instruction knob at all.** Anything from 16 KiB to 1 MiB is free.

---

## What the random-locus sampler cost

The sampler keeps loci by a **hash threshold on the position** — `mix64(contig, pos) % n`, a
pure function of where the site is — not by reservoir sampling. The spec asks for exactly
this (`parameter_prepass_census_sites.md` §3: *"keep position p if … hash(contig, p, seed) <
threshold"*), for a reason the sketch also depends on: every arm then selects the identical
set, which is what makes a bit-exact cross-arm gate possible at all. A fixed rate was
sufficient; nothing in the statistics needed a reservoir.

Priced by sweeping the rate on HG002 30× `chr1`, 1.54 M loci, and differencing against the
sampler switched off. Measured, minimum of three:

| arm | sampler off | 1 in 15 (102,554 kept) | per kept locus |
|---|---:|---:|---:|
| A | 32.1508 G | 32.1975 G | **456** |
| B | 29.7261 G | 29.9468 G | **2,152** |
| C | 30.1650 G | 30.3829 G | **2,125** |

At the real rate — 100,000 loci out of a human genome at 30×, which is **1 site in about
29,000** — copying kept loci out of blocks costs **0.22 G instructions** against roughly
5.5 × 10¹³ for the run: **one part in 250,000**. It is as small as the plan expected. Arm A
pays 456 instructions rather than 2,155 because it does not copy at all; it moves an object
it already owns.

The kept set's **memory** is the larger half of the story, and it favours the block: keeping
102,554 loci raises peak RSS by **48.6 MB in arm A and 37.0 MB in arms B and C**, because
`copy_out` allocates each vector at exactly its length while the walk's own `collect` leaves
slack.

One honest deflation. The spec's census does **not** store loci — it stores *"a dense array
in position order … one binned depth in four bits"* plus a sparse list of non-reference
observations, about 1 MB for two million sites. Keeping whole loci is what the plan asked
this sketch to price, and it is an **upper bound** on what the real census will copy.

---

## Correctness — the arms do the same work

**The shipped walk did not change behaviour.** All four acceptance dumps, `cmp` against the
stored copies in `tmp/perf_review_2026-08-04_ng-generic-walk/`:

```
ng_ssr_loci_dump     chr21       byte-identical
ng_ssr_loci_dump     SL4.0ch01   byte-identical
ng_generic_loci_dump chr21       one line: record_widen_events=423 → 425   (accepted)
ng_generic_loci_dump SL4.0ch01   one line: record_widen_events=622 → 621   (accepted)
```

Those are exactly the two accepted divergences and nothing else. (The stored tomato baselines
were taken on the **sparse bench CRAM**, not the 49 GB one — their header reads
`reads_admitted=105894`. My first attempt ran them against the big CRAM, produced a 5 GB dump
and a meaningless diff, and cost half an hour; noted here because the file names do not say
which fixture they came from.)

`cargo test --lib`: **2,893 passed, 1 failed**, and the failure is
`parity::every_divergence_from_production_is_one_of_the_six_named_classes` — the accepted
clean-tree divergence. That is the clean-tree baseline exactly.

**The four arms agree bit for bit.** Run to exhaustion, with the payload digest on — which
hashes *every field of every locus*, including `q_sum`'s bit pattern, not merely what the
pre-pass reads:

```
chr21, shipped typed-region stream
  A  loci=236081 covered=236505 visits=4616340 observations=251786
     rg_digest=790d3e89c23840d3 win_digest=cd20cdc40ff177a2
     kept=498 kept_digest=793fd7a48a735dd9 payload_digest=8d25867e6967cb51
     records_outside_region=127231  mate_overlap_positions=39312
     fit rg=0 phred=24.50 het=0.00316 hom_alt=0.00063 lnL=-125245.798173
  B, C, Cb — identical, every field

tomato sparse CRAM SL4.0ch01, whole contig
  A  loci=1851982 covered=1852975 visits=8030465 observations=1859657
     rg_digest=05f61fe6e3f6b09f win_digest=b452a623b792907e
     kept=3670 kept_digest=13e577e30baca06b payload_digest=dadb3728df3b1c2f
     fit rg=0 phred=30.25 het=0.00040 hom_alt=0.00063 lnL=-75269.841244
  B, C, Cb — identical, every field
```

`loci=236081 observations=251786` on chr21 reproduces the shipped probe's own gate figures
(cited).

**No floating-point tolerance was needed anywhere.** The accumulators are integer; the fit's
`f64` sums run over a dense cell table in a fixed index order, and the fitted log-likelihood
agrees to all seventeen digits printed. That was designed for rather than discovered: all
four arms call **one** `add_site`, and only the ten to fifteen lines that walk a locus
differ.

**The pre-pass stand-in accumulates what the spec says**, so it is not measuring plumbing:
two objects (a `(depth, alt-reads)` histogram keyed by read group × ploidy, and one keyed by
contig × 100 kb window × ploidy), a ragged cell table with depth bins exact to eight and
widening geometrically to 124, `depth_sums` **per cell** rather than per bin, stochastic
downscaling of over-deep sites keyed on the position, complete witnesses only, and no base or
mapping quality anywhere. At the end it runs the real fit — a profile scan over the spec's
161-rung Phred ladder against a grid of genotype frequencies — so that nothing accumulates
into an object no one ever reads.

---

## What I changed about the design while building it, and why

**1. A runtime enum, not a type parameter.** Making the walk generic over its sink would
ripple a parameter through `PileupGenerator`, `GeneratorSlot`, `GeneratorSet` and
`SampleLocusObservationsIterator` — including the STR generator, which this sketch does not
touch. `LocusSink` is one enum with a `wants_records()` predicate instead. The cost is one
perfectly-predicted branch per locus against a per-locus cost in the thousands; the benefit
is that arm A's code path is unchanged and the acceptance dumps still pass.

**2. The fast lane's held locus had to be re-solved.** The ordinary-column path
(`fast_column.rs`, 78 columns in 100, cited) does not emit at the base it walks: it holds the
locus one step, because that is where the general path drains the one-base record it would
have left open. A general record ending at `p−1` is still in the table when the fast lane
fires at `p`, so **releasing the fast locus immediately would put it ahead of that record**.
Given the cited split — 27 columns in 100 reach the general path — the interleave that makes
this matter is roughly one column in five. The columnar path cannot hold a *record*, so it
holds the accumulating scratch instead and swaps it with a second buffer (a pointer swap, no
copy), and it releases at the top of `process_position` rather than inside the fast branch,
which is order-identical because nothing between those two points writes to the sink. I found
this by reading the existing comment on `sealed`, which states the invariant plainly; had I
not, the payload digest would have caught it, but only after a wrong measurement had been
taken.

**3. `Vec::clear` is the enemy of the thing being measured.** Both the fast lane's observation
buffer and the general fold's keyed-observation buffer hold `Vec`s inside their entries.
Clearing the outer vector **drops the inner ones**, which frees exactly the allocations the
block layout exists to keep. Both had to become pools with a live-prefix counter, and the
fast lane's reset has to clear *every* entry rather than the live prefix, because the buffer
arrives by a swap and which of its entries are dirty is the other buffer's bookkeeping. This
is the subtlest thing in the sketch and it is invisible in the output — it costs allocations,
not correctness, so no test can fail on it. It is the reason the allocation count is a
first-class result and not a diagnostic.

**4. The general fold shares one function between the arms rather than having two.**
`keyed_observations_counting` became a thin wrapper over `keyed_observations_into`, which
takes the caller's buffer. Handed a fresh `Vec` — arm A — every add pushes and the behaviour
is what it always was. Handed a pool it reuses in place. That kept one copy of the fold's
fifty lines. **The close did not get the same treatment**: `finalise_into_columns` duplicates
everything from the witness tally to the canonical sort and has to be kept in step with
`finalise_recycling` by hand. That is the shape of the cost production reported as *"the
range-append is written twice"*, arrived at independently and within a day.

**5. The block outlives the walker.** A chromosome boundary mints a new walker; the block must
not be reset there, because it is the consumer's handover unit and knows nothing about
chromosomes. The sink is taken off the retiring walker and installed on the new one. Three
lines, and without them a multi-chromosome run silently loses whatever the block held.

**6. The region clamp moved.** The record path drops out-of-region loci in
`PileupGenerator::next_locus`, after the walker has yielded. The columnar path has no such
moment, so the sink is told the region bound and refuses before writing. Both count
`records_outside_region`, and the gate confirms they agree — 127,231 on chr21, in all four
arms.

---

## What building each arm was actually like

The plan asked for this as directly as it asked for the numbers, so here it is plainly.

**Arm A was free**, because it already exists. Reading a locus is thirty lines and every one
of them says what it means.

**Arm B's consumer was the shortest to write and is the one whose types protect it least.**
It is 48 code lines against arm A's 61, and it is genuinely shorter — indexing
`obs_num_obs[j]` beats reaching through two levels of ownership. But **every field access is
an unchecked index into a different array, and the two index spaces are both `usize`**:
indexing a per-observation array with a locus index compiles, runs, and returns a plausible
number. Nothing in arm B is protected by a type; what protects it is the payload digest, and
I would not have trusted any measurement taken without one. The arms agreed on their first
run — but that is a statement about how carefully the code was written, not about how hard it
is to get wrong.

**Arm C's consumer was the easiest to trust and the most tedious to build.** The loop is arm
A's loop, and I could have copied it — but it needs 62 lines of buffer machinery underneath,
and `refill_from` is thirteen field copies that must be kept in step with the block layout by
hand. The failure modes are opposite in the two arms: add a field to an observation and arm B
stops compiling at the one place the block is written, while arm C's `refill_from` compiles
unchanged and silently carries whatever was in the buffer before. **That is the reverse of
what I expected.** The record-shaped arm is safer where the *reader* is and more dangerous
where the *copier* is.

**Arm Cb took twenty minutes and is not on the table.** It reads exactly like arm C, costs
what arm B costs, and is the shape production retired. Building it was worth it only to learn
that arm C's copy is 1.5 % — a number I would otherwise have had to guess.

**The producer is where the work actually was**, and I did not expect that either. The
sketch's consumers are 48 to 131 lines each; teaching the walk to fill a block is **830 code
lines** across six files, three of which — `fast_column.rs`, `open_record.rs`,
`genome_walk.rs` — are the walk's hottest and most carefully commented. Two of the six design
decisions above (the held-locus ordering and the `Vec::clear` pooling) are the kind that are
correct or silently wrong with nothing in between, and one of those two cannot be caught by
any test the project has. **If the decision is made on the consumer's ergonomics it is being
made on the small end of the problem.**

---

## Instrument and its limits

- **`instructions retired` from `/usr/bin/time -l`**, floor-subtracted by differencing two
  locus counts of one fixture, minimum of three runs, arms alternated inside one script.
  Floors were also measured per binary per fixture for reference — HG002 `chr1` 0.170 G,
  tomato sparse 1.141 G, tomato big 1.296 G, arm A — and they cancel in every figure quoted.
- **`PVC_TRUST_REFERENCE_INDEX=1` throughout**, echoed as `reference_check=trusted_unverified`.
  No run here is compared against one that verified the FASTA.
- **No wall clock is quoted anywhere and none was recorded.** This host has 6
  high-performance and 12 low-energy cores, two other sketches were measuring concurrently,
  and a further agent was loading the machine for part of the session.
- **No sampling profile was taken and none is quoted.** Every attribution here is either an
  instruction-count ablation (the sampler; the pre-pass against the bare probe) or a static
  allocation tally from dhat. Both are per-process and unaffected by host contention. The one
  place I would otherwise have reached for a profile — *where the remaining allocations are* —
  is answered below from dhat's block counts instead.
- **Peak resident memory** from the same `/usr/bin/time -l`, per-process and unaffected.
- Every number above says whether it was measured here or cited. The cited ones are: the
  4,533-instruction column term, the 659-instruction visit, the 78-in-100 fast-lane share and
  27-in-100 general-path share, the 992-in-1,000 one-base record, the 165-to-1 typing ratio,
  and production's 2,501 MB → 261 MB block sweep.

### Where the allocations that remain are

From dhat block counts — a static tally, not a profile — HG002 30× `chr1`, 200,000 loci. The
sites that vanish between arm A and arm B are exactly the locus-building ones:

| site | arm A | arm B |
|---|---:|---:|
| `try_ordinary_column`'s per-observation `SequenceObservation` | 147,718 | **0** |
| the `Vec<SequenceObservation>` each locus collects into | 139,581 | **0** |
| `grow_one<KeyedObservation>` — the general fold's buffer | 120,874 | **0** |
| the reference-fetch buffer (`ref_seq.rs`) | 120,843 | 61,643 |
| **cloning each read's name (`String::clone`)** | **223,710** | **226,980** |
| **total, whole process** | **1,543,681** | **862,120** |

**The largest single allocation site in the walk is not the locus at all — it is cloning each
read's name**, and nothing in this sketch touches it. It is 223,710 of arm A's 1,543,681
blocks and 226,980 of arm B's 862,120, so it rises from 14 % of the allocations to 26 % of
them simply by standing still. If allocation count is the target, that is the next thing to
look at and it has nothing to do with the shape of the locus stream.

---

## What I would tell the owner

**The default survives, and if any arm displaces it the arm is C.**

- B and C land **1.5 % apart** on instructions and **identical** on allocations, peak memory
  and emitted bytes. The plan said *"if B and C land close together, that is a finding, and a
  large one: it would mean no consumer ever needs to see columns."* They did, and it does.
- The prize is **7.8 % of a covered base at the target depth**, falling to 2.0 % at 130× and
  rising to 12.1 % at 4.9×. It is 1,344 instructions per covered base plus about 8 per read
  at a position, so it is worth most exactly where breadth dominates — 30× whole-genome.
- The price is **830 lines of producer** in the walk's hottest three files, two of whose
  design decisions are silently-wrong-or-right, plus a second `finalise` kept in step with the
  first by hand.
- **Peak memory is not an argument either way.** The columnar arms are 0.2–0.4 % *higher*, and
  block size buys nothing between 16 KiB and 1 MiB.
- **The pre-pass is not what decides it**, exactly as the plan predicted. It is 7 % of a
  covered base, and the arms differ inside that by a fifth.

The one number I did not expect: **arm B with a full pre-pass costs what today's walk costs
with no pre-pass at all.** If the pre-pass has to come out of the walk's current budget, the
block layout is how it gets paid for. If it does not, the 8 % is a straightforward trade
against 830 lines and the stated preference for records everywhere.

---

## Files

Beside this report:

- `sketch1_producer.diff` — the whole producer change against `1e5ffa8`, six files.
- `sketch1_src_block.rs` — the new `LocusBlock` / `LocusSink` module.
- `sketch1_example_ng_prepass_sketch.rs` — the four arms and the pre-pass stand-in.
- `sketch1_raw/` — every measurement in run order (`hg30.txt`, `tomsparse.txt`, `tom130.txt`,
  `floors.txt`, `sampler.txt`, `sweep_hg30.txt`, `sweep_tom130.txt`, `probe_ref.txt`), the
  drivers that produced them, and `gate/` with the four dump diffs.

To reproduce any arm:

```
PVC_TRUST_REFERENCE_INDEX=1 PVC_PROBE_WHOLE_CONTIG=1 \
PVC_SKETCH_ARM=B PVC_SKETCH_BLOCK_KB=256 PVC_PROBE_MAX_LOCI=1540000 \
  /usr/bin/time -l target/release/examples/ng_prepass_sketch \
    $HREF benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam chr1
```
