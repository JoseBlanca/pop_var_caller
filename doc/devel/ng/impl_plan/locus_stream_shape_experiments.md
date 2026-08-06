# What shape the locus stream should have — three experiments

**Status:** plan, nothing built.
**Date:** 2026-08-06.
**Decides:** what ng's generic locus generator hands the stages after it, and where — if
anywhere — that shape should stop being one record per locus.

---

## 1. The pipeline this is about

1. **Produce the loci.** Built and shipped: one locus per covered reference base, single
   sample.
2. **Read every locus** to build histograms, estimate parameters, and keep about 100,000 loci
   drawn at random. Read-only; it alters nothing.
3. **Write to a `.psp` file and read it back.** Optional, but usually done, because it saves
   everything above.
4. **Merge the evidence of all samples**, k-way, position by position.
5. **Call the SNPs** (the expectation-maximisation work).
6. **Write the VCF.**

Steps 3 to 6 are close to what production already does. Step 1 is new and shipped; step 2 is
new and in flight on another branch.

## 2. The default architecture, and what an experiment has to beat

> **Records everywhere. Columns only in the file.**

That is the owner's position and it is the starting point, not a straw man. It is simpler,
every stage can be written the obvious way, and **some of the algorithms are not known in a
columnar form at all** — which is not a detail to be designed away later.

So the burden of proof sits on any departure from it. An experiment does not ask *"is
columnar faster here"* — it will usually be, in isolation. It asks **"is it enough faster,
here, to be worth writing this stage in a shape we find harder to reason about"**. There is no
threshold set in advance; the numbers get discussed.

**What "record" and "column" mean below.** A *record* is one owned object per locus, carrying
its own reference bases, its observations, and per observation an allele string and a list of
chain ids — what the generator emits today, at about thirteen heap allocations per covered
base. A *column layout* is a block of many loci where each field is one array shared across
the block, with a flat byte buffer plus an offsets array for the ragged fields, and an index
saying where each locus's part of each array begins.

## 3. What is already measured — do not re-derive any of it

**ng's generator, as shipped at `1e5ffa8`** (measured 2026-08-05, `tmp/perf_2026-08-05_shipped_reprice/price.md`):

| | |
|---|---:|
| per read at a position, BAM | 659 instructions |
| per read at a position, CRAM | 1,103 instructions |
| per covered base, independent of depth | **4,533 instructions** |
| one covered base at 30× | 18,116 instructions |
| hand-count of the work one read-at-a-position performs | 300–600 instructions |

The per-read cost is within a factor of two of the work it does; the per-base cost is not
where the reads are, and is roughly thirteen allocations at the measured ~340 instructions
each. **That is the number these experiments exist to attack.** At 30× it is a quarter of what
a covered base costs, and about 10¹³ instructions over a human genome.

**Production's cohort path** (measured, cited from `perf_var_calling_cohort_2026-06-06.md`,
`re_architecture_p6_measurement.md`, `src/psp/`, `src/var_calling/`):

- A row object between two columnar ends was called *"pure overhead: columnar → row →
  columnar"*, and removing it was the point of that redesign.
- Its replacement — slicing row ranges out of columns into new columns — then became the
  **largest heap item at 49 % of peak, 20.8 GB churned**, against 19 % for the per-locus row
  materialisation it replaced. **The layout moved the copy; it did not remove it.**
- **65 % of the in-flight payload was chain ids, and 96.6 % of chain ids were reference-allele
  ids the merger discards.** Storing less beat every layout change.
- Only about **4 positions in 100** are variable, and folding one cheap column across samples
  to find them means only those ever become objects.
- Block size is the memory knob: 16 MiB blocks → 1 MiB took peak resident memory from
  **2,501 MB to 261 MB** with wall time flat.
- **A borrowed per-locus view over the decoder's block lost to lifetimes.** `BlockColumns<'a>`
  exists and is test-only. A cursor holding one locus per sample across calls to `next()` is
  self-referential; both mergers own their lookahead for that reason, and both say so in the
  source.
- The offset arithmetic for the ragged columns is re-derived at **five call sites**, and the
  range-append is written twice.
- Blocks whose boundaries did not line up with where consumers wanted to cut needed a second
  mechanism: shared segments, a decode-once cache, two compaction paths, and a requirement
  that both produce identical bytes.

## 4. The three experiments

Each is a **sketch**: a few hundred lines, throwaway, in its own branch, deleted after the
decision. None of them commits an interface. Each carries enough of the real computation that
it is not measuring plumbing alone.

### Sketch 1 — what the generator emits, and what the parameter pre-pass reads

**The question.** Does emitting blocks of arrays instead of one record per covered base remove
enough of the 4,533-instruction per-base cost to be worth it — and does the pre-pass get
faster or only more awkward?

**Built.** Three arms over the same walk:

- **A (the default):** today's shipped generator, one record per locus. The baseline exists
  and is measured to the instruction.
- **B:** the walk fills a block of arrays that the caller owns and reuses; the pre-pass reads
  columns; loci are materialised nowhere.
- **C:** the walk fills the same block, but the pre-pass is handed one locus at a time through
  a view borrowing from a buffer **it** owns and refills per locus — the pattern that survived
  in production, rather than the borrow into the block that did not.

C exists because it is the honest middle: it keeps every consumer's code record-shaped while
removing the allocation, and if it lands close to B then the columnar consumer buys nothing
and the answer is C.

**Measured.** Instructions retired per covered base and per read-at-a-position, split the same
way as the existing baseline so the numbers are directly comparable; peak resident memory;
allocations per covered base; and, for the sampler that keeps 100,000 loci, the cost of
copying them out.

**What I expect, stated so it can be wrong.** B and C both remove most of the per-base cost;
the gap between them is small; the pre-pass is not the thing that decides it.

### Sketch 2 — the k-way merge

**The question.** Does the merge want columns, and what does the per-sample cursor cost when
it cannot borrow?

**Built.** Two arms over a real multi-sample cohort:

- **A (the default):** merge records directly — one owned locus of lookahead per sample, an
  O(N) scan over the heads at each output position. This is production's
  `per_position_merger` shape.
- **B:** fold one cheap column across samples first to find which positions are variable, then
  materialise only those. This is production's measured design and the sketch is to find out
  whether it still wins when the producer is ng rather than a `.psp` file.

**Measured.** Instructions retired end to end, peak resident memory, **bytes copied at the
handoff** — the number production's redesign got wrong twice — and how many positions ever
become objects.

**The specific thing to watch.** Production's cursors own their lookahead because a borrowing
one is self-referential. Whatever B does about that (shared blocks behind a reference count, or
a per-sample owned buffer the view borrows from) is the part to report honestly, because it is
the part that will still be there in five years.

### Sketch 3 — what the calling maths is handed

**The question.** For the expectation-maximisation step, does the input shape matter at all?

**Built.** The same calculation, twice: once taking a per-locus record, once taking per-locus
slices over a buffer the caller owns and reuses.

**Measured.** Instructions retired for the calling step alone, and the share of it that is data
movement rather than arithmetic.

**Why this one may end quickly, and that is a good outcome.** If the arithmetic dwarfs the
data movement, the shape at this boundary is free and the right answer is whichever reads
better — which, by section 2, is records. That result would remove a question rather than
answer it, and it is worth the day it costs. **This is also the sketch that meets the
"some algorithms we don't know how to do in columns" problem head on**: if a step cannot be
written columnar, the sketch is where that is discovered cheaply.

## 5. Fixtures — real data, and why synthetic would mislead

**Both of production's largest wins came from skew**: only 4 positions in 100 are variable,
and 96.6 % of chain ids are reference-allele ids. Synthetic data has neither property, so a
sketch driven by it would report that columns buy nothing, and it would be wrong. Every
measurement uses real alignments.

| fixture | what it is | used by |
|---|---|---|
| `benchmarks/tomato_big_cram/DRR000741.p1.cram` | 49 GB CRAM, ~130×, one sample | sketches 1 and 3 (deep) |
| `benchmarks/ssr_hg002/bam/30x/HG002_TR_v1.0.1_Tier_30x.bam` | human 30×, BAM, one sample | sketch 1 (the target depth) |
| `benchmarks/ssr_hg002/bam/{10x,300x}/…` | the same sample at two other depths | sketch 1, if depth changes the answer |
| `benchmarks/ssr_tomato1/crams/*.bench.cram` | 63 tomato CRAMs, one sample each | sketch 2 (the cohort) |
| `benchmarks/tomato1/results/ours/cohort/psp/*.psp` | 56 built `.psp` files | sketch 2, so the merge can be measured without re-walking |
| `/Users/jose/devel/pop_var_caller/tmp/aligned_psp/` | 50 aligned `.psp` files | sketch 2 at cohort scale |

**⚠ The human fixture is tandem-repeat-targeted, not whole-genome.** Walking its `chr1` types
165 reference bases for every covered base it delivers, and about half that run's instructions
are region setup rather than walking. Any per-base number from it must say which mode it came
from, and the tomato CRAM is the one that looks like real whole-genome work.

## 6. Instrument and discipline

Carried unchanged from the three performance rounds that produced the numbers in §3, because
they are what made those numbers trustworthy:

- **`instructions retired` from `/usr/bin/time -l`** is the instrument, floor-subtracted with a
  `PVC_PROBE_MAX_LOCI=1` run **measured per binary per fixture**, minimum of three runs a side,
  binaries alternated within one script.
- **Wall-clock comparisons are not admissible on this host** — 6 high-performance and 12
  low-energy cores — and criterion's `change:` line is not admissible across commits: two runs
  of an *identical* binary reported four of six points "statistically significant", spanning
  −9.2 % to +15.9 %.
- **Peak resident memory matters as much as instructions here**, because the block is the unit
  that must be live at once and block size is the known memory knob.
- **Allocations and bytes copied at each handoff are first-class results**, not diagnostics.
  They are what production's redesign got wrong twice.
- Sketches build natively on the host; `--features dhat-heap` always with
  `--target-dir target-dhat`.

## 7. What each sketch reports

One table per sketch, same columns, so three sets of numbers can be laid side by side:

| arm | instructions per covered base | per read at a position | peak RSS | allocations per base | bytes copied at the handoff | lines of source | what it cannot do |
|---|---|---|---|---|---|---|---|

**The last two columns are not decoration.** A sketch that wins on instructions and doubles the
source it takes to express the stage has not obviously won, and the place production paid for
this is on record: offset arithmetic re-derived at five call sites, a range-append written
twice, and a whole second mechanism for blocks that did not cut where consumers wanted to cut.

## 8. Non-goals

- **Not building a pipeline.** No sketch is a step toward an implementation; each is deleted
  after the decision.
- **Not re-deriving §3.** Production's cohort measurements stand. Sketch 2 tests whether they
  survive a different producer, not whether they were right.
- **Not choosing one shape for everything.** The plausible outcome is records everywhere with
  one exception, or records everywhere full stop.
- **No threshold decided in advance.** The numbers get discussed.

## 9. What could make this exercise worthless, so it can be watched for

- **A sketch that omits the real computation** measures plumbing that real work would hide.
  Each arm carries the actual arithmetic or a proxy of the same size.
- **Synthetic data**, for the reason in §5.
- **Per-stage microbenchmarks**, which favour columns everywhere by construction and hide the
  conversions — which is where the cost turned out to be. Every number is taken over a whole
  path, even when a sketch only varies one stage of it.
- **Measuring a stage nobody was going to change.** Step 6, writing the VCF, is record-shaped
  by nature and is not in any sketch.

---

## 10. Results — all three sketches, 2026-08-06

Full reports, raw measurements and sketch code in gitignored
`tmp/locus_stream_experiments_2026-08-06/`. Every number below was measured on real
alignments; no wall clock was recorded and no sampling profile was quoted, because a further
agent was loading the host — every attribution is an instruction-count ablation or a static
allocation tally, both per-process.

### The default survives. Two of the three questions are closed.

**Sketch 3 — the calling boundary: records, and it is not close.** The
expectation-maximisation step costs **999,224 instructions per called locus**; materialising
that locus as an owned record costs **5,999**. The arithmetic is **167 times** the data
movement. Handing the step a record instead of borrowed slices costs **0.59 %** of the calling
step, and **0.10 %** once it must also produce the record the VCF writer consumes. Records are
85 lines against 162 for the same stage, and hold one locus resident instead of a 242 MB
block. Calls bit-identical, 0 of 23,552 loci diverging, no tolerance allowed.

*Nothing in the arithmetic resisted columns* — a borrowed-slice entry point already existed —
but the **output** side did: pruning unsupported alleles changes one locus's ragged widths,
which shifts every offset after it in a block. And because the VCF writer wants a record, the
columnar arm copies the passthrough columns back out: 32.44 allocations per locus against the
record arm's 33.44. **One allocation apart.** The ratio also moves the right way: 167 at 50
samples, 35 at 10, where the record penalty rises to 2.6 %.

**Sketch 2 — the merge: it depends entirely on what is on the other side.** Reading `.psp`,
folding one cheap column across samples before materialising anything is **2.00× fewer
instructions**, 414-fold fewer allocations, 252 MB copied at the handoff against 4,638 MB, and
5.5 M objects materialised against 52.2 M. Reading **ng's generator directly**, the merge is
**1.6 % of the whole path** and the two shapes are indistinguishable (26.542 against
26.565 × 10⁹ instructions). The mechanism: the generator's only emitter returns a fully owned
record, so there is no cheap summary to fold and "materialise late" saves a copy the record
arm was never going to make.

**So the merge is not a reason to make ng's in-memory stream columnar.** It is a reason to
fold-then-materialise when reading the file — which is about `.psp`, whose columnar shape was
never in question. The fold arm is **3.9× the source** of production's merger and its author
found it materially harder to write correctly; both bugs hit were in cursor and window
bookkeeping that exists only because a view cannot outlive a block. Peak heap went the *wrong*
way by 4.8 %.

The borrowing constraint was solved without owning or reference-counting: N distinct readers
can all be borrowed *shared* at once, so the loop is load-mutably → fold and materialise from
shared views → drop views → advance mutably. Zero bytes copied, 3.87 MB of scratch, and one
permanent bound — the merge can only advance to the minimum across samples of the last
position in each currently-loaded block.

### The one question left open, and it is inside the generator

**Sketch 1 — blocks remove 1,344 instructions of the 4,533 a covered base costs regardless of
depth.** Four arms, all agreeing bit-for-bit on a digest over every field of every locus:

| arm | instr / covered base, 30× BAM | allocations / base | peak RSS 30× | lines: consumer / producer |
|---|---:|---:|---:|---|
| **A** records (today) | 19,757 | 7.72 | 149.2 MB | 61 / 0 |
| **B** columns end to end | **18,216 (−7.8 %)** | 4.25 (−45 %) | 149.8 MB | 48 / 830 |
| **C** record view over the consumer's own buffer | **18,482 (−6.5 %)** | 4.25 | 149.8 MB | 131 / 830 |
| **Cb** view borrowed into the block (production retired this) | 18,199 (−7.9 %) | 4.25 | 149.7 MB | 66 / 830 |

Against depth: **−12.1 % at 4.9×, −7.8 % at 20.7×, −2.0 % at 98×.** The saving is a fixed
1,344 instructions per covered base plus about 8 per read at a position, so what changes with
depth is the denominator. At the stated target — one human sample at 30×, ~3.1 × 10⁹ covered
bases — it is **4.8 × 10¹² instructions**.

**B and C land 1.5 % apart, which is the finding §4 flagged as the large one: no consumer ever
needs to see columns.** The block can be a producer-side implementation detail.

**Peak memory is not an argument either way** — the columnar arms are 0.2–0.4 % *higher* — and
**block size is not the lever here that it is in the `.psp` reader**: instructions are flat
from 4 KiB to 16 MiB. The random-locus sampler is free: copying 100,000 kept loci out of
blocks is one part in 250,000 of the run.

**The price is 830 code lines of producer** — a new 367-line module plus 516 added and 53
removed across five files, three of them (`fast_column.rs`, `open_record.rs`,
`genome_walk.rs`) the walk's hottest. Two of its six design decisions are correct or silently
wrong with nothing in between, and one of those — pooling the observation buffers, because
`Vec::clear` drops the inner `Vec`s and frees exactly what the block layout exists to keep —
**cannot be caught by any test the project has**, since it costs allocations rather than
correctness. The close had to be written twice (`finalise_into_columns` beside
`finalise_recycling`, kept in step by hand), which is production's *"the range-append is
written twice"* arrived at independently within a day.

**The number that frames the decision: arm B with a full parameter pre-pass costs what today's
walk costs with no pre-pass at all** — 10,054 against 10,052 instructions per covered base on
the sparse tomato fixture. What the layout removes is about what the whole pre-pass adds. So
the question is not "is 7.8 % worth 830 lines" in isolation; it is **whether the pre-pass has
to fit inside the walk's current budget**. If it does, blocks are how it is paid for, and arm C
is the shape, because every consumer keeps its record-shaped code. If it does not, the trade is
7.8 % against 830 lines in the three files that most need to stay correct.

### Sketch 4 — the combination the first three could not measure

A block-filling generator feeding a fold-first merge in memory. 10 tomato CRAMs in lockstep,
300 regions, 2,830,932 loci, **1 position in 101 variable**, depth 4.55×. Instructions
retired, minimum of five runs, all modes alternated in one script, spread 0.07–0.20 %:

| state | producer | merge | whole path | merge alone (ablation) | peak RSS | loci materialised |
|---|---|---|---:|---:|---:|---:|
| **A** | records | records | 27.034 G | 0.448 G | 260.1 MB | 2,830,932 |
| **B** | records | fold | 27.007 G | 0.420 G | 268.5 MB | 2,830,932 |
| **C** | blocks | records | 23.932 G | **1.642 G** | 262.3 MB | 2,830,932 |
| **E** | blocks | records + keep column | 22.820 G | 0.384 G | 262.5 MB | 28,718 |
| **D** | blocks | fold | **22.674 G** | **0.239 G** | 262.4 MB | **0** |

**The design point §2 of this sketch predicted is real: the block as sketch 1 built it cannot
be folded.** It carries no cheap per-locus summary, and deriving one means walking the
heaviest array. A `locus_nonref_obs` column was added — 29 lines, all inside the block module,
no hot walk file touched, **4 bytes per locus and 51 instructions per locus** in a deliberately
naive form. **The merge's needs shape the generator's output layout**, and that is a standing
constraint rather than a one-off.

**⚠ Arm C is a trap, and it is the obvious thing to build.** Blocks with a merge that
materialises records is **worse at the merge than today** — 1.642 G against arm A's 0.448 G —
because the block-to-record copy is work the record producer never had to do. A block producer
without the keep column gives back a third of what it saves.

**Arm E is the finding.** It is record-shaped — the same head scan as arm A — and simply reads
the keep column, materialising 28,718 loci instead of 2.83 M. It banks **90 % of the fold's
saving** (C→E is 1.258 G of 1.403 G; the dense-window fold adds the remaining 0.146 G) with
**no round structure and no shared-borrow phase**. Sketch 1's *"no consumer ever needs to see
columns"* arrives a second time by a different route.

**The merge is one twentieth of the story.** A→D is 4.360 G, of which **95.2 % is the producer
and 4.8 % is the merge**. The fold's prize does **not** grow with sample count — flat from 3 to
40 samples, and the ratio between arms C and D actually falls from 7.8× to 6.1×. Sketch 2's
memory penalty for the fold (4.8 %, 10 MB) **vanishes** with a block producer, because the
block already is the buffer.

**⚠ This fixture is 4.55× deep, which is where sketch 1's depth curve says the producer wins
most.** The 15.6 % producer saving here is *not* the 30× number; sketch 1 measures −12.1 % at
4.9×, −7.8 % at 20.7× and −2.0 % at 98×. A clean-tree control also shows sketch 1's runtime
sink enum taxing the record path by 1.81 %, so the honest A→D range is **14.6–16.1 %** on this
fixture.

All five states and both floors agree bit-for-bit on sketch 2's own digest
(`0xf46f8f924ee468aa`, `loglik_acc = -56602.839011`), reached by routes sketch 2 did not have,
and invariant across block budgets from 4 KiB to 1 MiB.

**So the merge does not change the standing decision.** At the 30× target the producer is worth
about 7.8 % and the merge under 1 %, against 830 lines in the walk's hottest three files plus 29
in the block module. If the 7.8 % was not worth it before, the merge does not make it worth it —
**but if it is taken, take it as arm E**: blocks inside the producer, a four-byte keep column,
and a record-shaped merge that reads it. That is the cheapest 90 %, and it leaves every
consumer record-shaped.

### Two findings independent of the whole question

- **The largest single allocation site in the walk is cloning each read's name** — 223,710
  blocks of arm A's 1,543,681, untouched by any arm here, and it *rises* from 14 % of
  allocations to 26 % under arm B simply by standing still.
- **The calling step allocates 13.85 blocks per locus in both arms**, 1.8× what materialising
  a record costs: `validate_record_shape` calls `genotype_order`, which builds and sorts a
  nested vector purely to compare its length and then discards it. `shape_for` already caches
  that table.
