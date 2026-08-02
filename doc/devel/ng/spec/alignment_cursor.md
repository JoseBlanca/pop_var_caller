# ng — reading a sorted BAM or CRAM forwards, once

**Date:** 2026-08-01 · **Status:** design agreed with the owner, **no code yet**.
**Companions:** [alignment_file.md](alignment_file.md) (the current design this replaces),
[sample_reads.md](sample_reads.md), [ref_seq.md](ref_seq.md) (the same split, for the reference),
[locus_generation_pileup.md](locus_generation_pileup.md) (the caller).
**Code-facing companion:** [`../arch/alignment_cursor.md`](../arch/alignment_cursor.md) (types &
interfaces). **Build order:** [`../impl_plan/alignment_cursor.md`](../impl_plan/alignment_cursor.md).
**Evidence:** [perf_ng-generic-pileup_2026-07-31.md](../../reports/reviews/perf_ng-generic-pileup_2026-07-31.md).

---

## Words used in this document

A few file-format terms appear throughout. They are all about *how the bytes are stored*, not
about biology.

- **Block-compressed file (BGZF).** A BAM file is not one compressed stream. It is a chain of
  small blocks, each compressed on its own, each about 64 KB when unpacked. To read one read you
  must unpack the whole block it sits in.
- **Index** (the `.bai`, `.csi` or `.crai` file beside the alignment file). A small sidecar that
  answers one question: *which byte ranges of this file could hold reads overlapping this stretch
  of the chromosome?*
- **Byte range.** One such answer — "between byte X and byte Y". The index usually returns
  several per question. (noodles calls these *chunks*.)
- **Cursor.** The new thing this document describes. In Python terms: `AlignmentFile` is the
  *iterable* and a cursor is an *iterator* over it — it remembers where it got to.

---

## 1. The problem

ng reads alignment files almost entirely forwards. The region walk produces regions in order
along the chromosome and asks, for each one, "give me the reads here".

Today each of those questions is answered from scratch. We ask the index where to look, jump to
that byte, unpack the block, and read forward. Then the next region arrives and we do all of it
again — usually landing in **the same block we just unpacked**.

Measured on HG002 at 30× depth, chromosome 21. **Two caveats that travel with every number in
this document.** The alignment file is tandem-repeat-targeted — 0.64 % of chromosome 1's positions
are covered — while the workload this is for is 30× whole-genome, so *ratios transfer and
microseconds do not*
([perf review §1](../../reports/reviews/perf_ng-generic-pileup_2026-07-31.md)). And the
region-grain sweep that bounds the available saving was taken with the probe's whole-contig mode:
its region-size knob only *splits* generic regions
(`examples/ng_generic_walk_probe.rs:171-195`), so a 10 kb or whole-contig row cannot be
reproduced against the real typed-region stream, whose generic regions average 392 bp. That
ceiling is a ceiling on a synthetic mode — which is also why caller-side coalescing is not
available (§13, N2).

| | |
|---|---:|
| jumps that land in the block already unpacked | **82 %** (27,731 of 33,671) |
| share of the walk spent on those jumps | **30 %** |
| times the same 35,228 reads are unpacked and decoded | **1,067,729 — about 30 times each** |

The waste is built into the shape, not into a bad line of code, and it comes from **how our
caller asks**. The pileup generator does not request the region it is walking; it requests a
window 5,000 bases wider, because a locus anchored inside the region can be stretched by a
deletion until the reads supporting it lie outside the region entirely
(`generator.rs:820-850`). Regions average about 390 bases, so each request overlaps the one
before it by most of its width: **roughly 93 % of what a request decodes was already decoded by
the request before it.**

That widening is the caller's business, not this type's — see the responsibility split in
section 2 — but it is what makes remembering decoded reads worth doing at all.

The unpacking is not something ng can skip. The library we use, noodles, re-unpacks
unconditionally: its `seek` does a file jump, a read, and a full block decompression with a
checksum, and it never checks whether the block it already holds is the one being asked for
(noodles-bgzf 0.47.0, `src/io/reader.rs:175-186`). There is also no way to ask it to move to a
different read inside the block it is already holding.

Caching the compressed blocks would not help. A cache of compressed bytes still has to
decompress, and decompression is what the profile charges for — a cache saves the `lseek` and the
`read`, not the inflate and the checksum. (An earlier draft argued this from a cold-jump timing of
103 µs against 68.6 µs warm; that is a single sample, and the structural reason above does not
need it.)

**So the fix is not to jump.** Read forward, and remember the reads already decoded so the next
region can be answered from memory.

### Two changes, and only one of them is what the measurement bought

This document proposes two things. They are separable, they are worth very different amounts, and
conflating them would let an architecture preference borrow the authority of a number it did not
earn.

| | what it is | what it is worth |
|---|---|---|
| **Keeping decoded reads** | remember what was decoded; reuse it when the next region can use it | **the whole saving** — the walk falls from 5.18 s to 2.69 s on chromosome 21 |
| **Changing who owns the reader** | a long-lived per-worker cursor instead of a per-region borrow from a shared pool | **noise today** — the pool's locks are 3 samples in 25,446, and the *entire* per-region query setup is 0.58 µs, 0.9 % of the walk |

**Keeping decoded reads does not require the ownership change.** The reuse rule has to be total
anyway (§6), so it does not care whether the reader was owned or borrowed: a pooled reader holding
some other region's records simply fails the check and jumps, which is today's behaviour. This is
not a hypothesis — **the CRAM path already retains across queries through the pool**, keeping its
last decoded container on the pooled handle (`open_bam.rs:169-172`) and skipping the decode when
it is asked for the same one again (`region_query.rs:711-718`).

So the ownership change stands on its own argument, not on the measurement: it makes the retention
rule live in one place instead of being threaded through a handle, it removes the per-region
reference-accessor factory (perf review L2, which drops a type parameter from `PileupGenerator`),
and it is the shape a per-worker fan-out needs. Those are worth having. None of them is 48 %.

**Both are being built together (owner, 2026-08-02.)** The alternative — keep the current
per-region API, put the retained records and the §6 comparison on the pooled reader handle, and
leave the cursor for later — was put by adversarial review and declined. It would have bought the
whole measured saving while touching only `region_query.rs` and `open_bam.rs`, leaving the walker,
the generator and the filter untouched; against that, it keeps the reuse rule threaded through a
pooled handle, keeps the per-region reference-accessor factory, and does not move towards the
per-worker shape the fan-out needs. The reason for the decision is the owner's and is not recorded
here.

What the split above is still for: **when the result is measured, the number belongs to
retention.** Do not report the walker rewrite as though the profile justified it.

### Why this was not done earlier

[alignment_file.md §3.3](alignment_file.md) chose the current design and explicitly left the door
open for this one: *"A sorted batch of loci could later be served by one forward sweep instead of
N seeks… the interface must not preclude it; but it is **not built now** — the container cache
recovers most of the CRAM cost and **BAM seeks are cheap**."*

Only that last clause turned out to be wrong, and the numbers above are what showed it. Everything
else in that section still holds.

## 2. What we are building

A **cursor**: a reader used by one part of the program at a time, which remembers where it got to
in the file and keeps the reads it has recently decoded, so a nearby region can be answered
without unpacking anything again.

**Goals.**

1. **Consecutive requests do not unpack the same data twice.** When the next region needs blocks
   the last one already unpacked, they are reused rather than decompressed and decoded again.
   This is deliberately *not* a promise that a read is decoded only once in a run: two cursors on
   one file each keep their own copy (section 7), and a long jump throws everything away.
2. **One rule per file format, short enough to say in a sentence.** The rule that decides "can I
   answer this from what I hold, or must I jump?" lives in one file and is tested there. It is a
   check the cursor makes, never a promise the caller has to keep.
3. **The cursor knows exactly which reads it may forget**, from one comparison of two region
   starts — no window size, no capacity, not even a question to the index. Section 6 gives the
   rule, the argument that it is sufficient, and why the obvious alternative is unsound.
4. **When in doubt, do what the code does today.** If carrying on is ever not provably safe, the
   cursor jumps and re-reads. That is the current behaviour, so the fallback carries no new risk.

**Non-goals.**

- **This is not a parallel design.** One cursor serves one consumer. More parallelism means more
  cursors, which is why nothing inside a cursor is shared or locked.
- **This is not a cache with a policy.** No capacity, no eviction strategy, no knob to tune. What
  is kept is what the index says is still reachable.
- **This does not change what we call.** Regions, the walk, the loci emitted — all unchanged. The
  test for "did we get it right" is that the output is identical.

**It does not** reorder regions, merge them, widen them, serve two consumers at once, or decide
what a region means.

**The responsibility split, stated because it is easy to get wrong.** The cursor returns the reads
overlapping the region it was given — nothing more and nothing less. Deciding *which* region to
ask for, including asking for a wider one than you mean to walk, belongs to the caller and already
does: the pileup generator widens the region itself before calling. A read layer that second-
guessed the region would be deciding something it cannot know, since only the caller knows what it
intends to do with the reads.

## 3. The shape: one shared file, one reader each

Two types, and the split is the one the reference reader (`RefSeq`) already uses.

**`AlignmentFile` is the iterable.** Everything settled when the file is opened and never changed
after: the path, the header, the parsed index, the chromosome list, the read-group table, and —
for CRAM only — a **way to fetch reference bases**, which CRAM needs because it decodes against
the reference.

That last one is a way to get bases, not the sequence itself. `OpenReference` holds a lazily-built
accessor (`bases: OnceLock<fasta::Repository>`, `reference.rs:116-143`), so a run over BAM never
opens the FASTA at all. A run over CRAM keeps resident only the chromosomes that a cursor is
currently reading — normally one, briefly two while workers cross a boundary — and frees each when
the last cursor on it goes away. Section 10 gives the rule and why the two modes that exist today
do not fit. Holding a whole reference would cost more memory than the entire walk.

**`AlignmentCursor` is the iterator.** The open file descriptor, the position in the file, the
reads being kept, the last region served. None of this can be shared, ever — it is one consumer's
place in the file.

| | shared between threads? | how |
|---|---|---|
| path, header, index, chromosome list, read groups | yes — never changes | one shared copy |
| the CRAM reference bases | yes — one copy per chromosome being read (section 10) | shared while a cursor needs it, freed when none does |
| file descriptor, position, kept reads, last region served, tallies | **no** | owned by one cursor |

**Nothing on the shared side is written to, and the file keeps no counters.** It keeps two kinds
today and neither survives. `readers_opened` exists only to catch a reader leaked back to the
pool (`open_bam.rs:704`), and there is no pool. It is read at **ten sites across nine tests**
plus `Debug` (`:742`, `:1087`, `:1106`, `:1123`, `:1134`, `:1162`, `:1238`, `:1254`, `:1428`,
`:1554`, `:1614`) — so deleting it is a test edit, not a one-line removal. The per-read-group
drop tallies stay where they already are — `ReadFilter::counts()`, a running tally readable at any
point (`filtering.rs:797-799`) — and the filter now lives as long as the cursor, so a caller reads
`cursor.counts()` whenever it likes. No atomic, no hand-over.

**A worker takes a cursor once and keeps it** — it does not borrow one for each region and give
it back. Three things follow, all of them wanted.

*The lock on the read path disappears.* Today a shared pool hands out readers, and taking one
costs three lock operations per region — about 1.8 million lock pairs on chromosome 1. Profiling
put that at 3 samples out of 25,446, so this is not a speed argument. An entire mechanism simply
stops existing.

*The drop tallies stay readable, and no hand-over appears.* Today the per-region stream folds them
into the file when it is dropped (`open_bam.rs:263-282`). A cursor does not need to: the tally
already lives on `ReadFilter` as a running count readable at any point (`filtering.rs:797-799`),
and the filter lives as long as the cursor, so a caller reads `cursor.counts()` whenever it wants.

Worth stating because an earlier draft got this wrong in both directions. It claimed a drop-order
defect here that does not exist — `locus_generation/mod.rs:700-701` declares `generators` before
`reads`, and `:776-800` adds an explicit `Drop` so the property does not rest on field order — and
then invented a hand-over obligation to replace it. Neither is real.

*"One reader per worker" stops being a comment.* `locus_generation.md` §9 already requires it.
Nothing enforces it today; here the types do.

### The caller keeps the cursor; there is no stream object

A caller holds the cursor as an ordinary value, moves it to a region, and pulls reads out one at
a time into a buffer it owns:

```text
cursor = file.cursor(chromosome, reference)

cursor.move_to_region(region)
while read = cursor.next_read():  use read

cursor.move_to_region(next_region)      # same cursor, still holding what it kept
while read = cursor.next_read():  use read
```

**Deleting the stream object removes the difficulties it created** — the lifetime, giving the
cursor back, a caller forgetting to. It does **not** leave nothing to get wrong: the generator
already performs exactly this take-and-give-back dance for the run-lifetime chain-id allocator
(`generator.rs:869-873`, taken on open with an `expect`, returned at `end_walk`), and moving the
cursor into the walker below adds a second one.

One smaller property falls out: two readers of one cursor are impossible, because pulling a read
needs exclusive use of the cursor and there is no second object to hold one.

**`next_read` hands back an owned read; it does not fill a caller's buffer.** An earlier draft did,
justified as "no read costs an allocation". That justification was wrong twice. The reuse trick
exists at the *record* seam — `RecordSource::read_next(&mut buf)`, *"reusing its allocations"*
(`filtering.rs:376-379`) — and that stays, one layer below the cursor. It does **not** exist at
the *read* seam: `decode` returns an owned `AlignedRead` (`filtering.rs:348`) and builds four
vectors doing it (`aligned_read.rs:71`, `:102`, `:103`, `:109`). Promising otherwise meant
rewriting that into a fill-a-buffer form across the trait and its implementors — and this
project's own measurements say that class of change buys nothing: removing 36 % of all
allocations moved wall time by −0.5 %, and swapping the allocator outright by +2.4 %
([perf review §2](../../reports/reviews/perf_ng-generic-pileup_2026-07-31.md)). An owned read
also leaves `Iterator` available to whatever consumes the cursor.

**⚠ What this costs elsewhere, costed honestly — it is the largest single edit in the plan, and
an earlier draft of this paragraph dismissed it in one sentence.**

- **The walker.** It owns a `Peekable` of reads and peeks at construction
  (`genome_walk.rs:103-113`) — the peek H4 measured at 50.6 % of the walk. It cannot borrow a
  cursor the generator also holds. So either the cursor moves *into* the walker and
  `move_to_region` is forwarded through it, which makes the walker long-lived and requires a
  per-region reset of `WalkerState`, `pending`, `done` and `stop_after`; or the walker stops being
  an `Iterator` and takes `&mut cursor` per step. Either is a substantial edit to a 1,276-line
  file, still shadowed by a 4,233-line differential harness.
- **The filter seam — solved, at the price of one accessor.** `ReadFilter` owns its source,
  exposes no mutable access, and returns it only by being consumed (`filtering.rs:800`, `:949`).
  A cursor yielding `AlignedRead` must own a `ReadFilter`; a `ReadFilter` owning the *cursor*
  would be a cycle, but it does not have to — the filter's source is the layer below the cursor.
  Adding `ReadFilter::source_mut` lets `move_to_region` reposition through it. See
  [the arch doc §2.3](../arch/alignment_cursor.md) for the layering.

The walker cost is real and is not bought by the measurement: §1 shows the saving comes from
keeping decoded reads, which needs none of it.

## 4. What the cursor promises

**A cursor covers one chromosome. Within it, ask for any region in any order and the answer is
always right.** Whether it is *fast* depends on how close the region is to the last one.

There is no ordering rule for the caller to keep. On each `move_to_region` the cursor decides
whether what it already holds can answer part of the question:

- **partly held** — the region starts inside what is kept and runs past it. **This is the common
  case and the one this type exists for**: a walk moving forward hits it at every step. Hand over
  the kept reads first, then carry on reading from the file. No jump, and nothing unpacked twice.
- **entirely held** — a region wholly inside what is kept, which happens when the caller looks
  back a little. Answered with no file access at all.
- **not held** — a big jump forward past everything kept, or **any** region starting before the
  last one served. Drop what is held, jump, read. This is exactly what the code does today for
  *every* region, so it is neither new nor risky.

**Reuse is partial.** The condition is only that the new region starts at or after the last one
served (§6) — *not* that the whole region fits in what is kept. A region running past the kept
reads is the normal case, and it is served by handing over the kept ones and then reading on.

**Blocks are unpacked only when the reads before them have been handed over.** A region spanning
two blocks yields everything from the first before the second is touched, and if the caller stops
early the second is never unpacked at all. The stream stays lazy, which
[locus_generation_pileup.md](locus_generation_pileup.md) §7 requires — and the order is right for
free, because reads kept from an earlier read are in file order, which in a sorted file is
position order, and reading forward continues it.

**The check must be total, and this is the one thing that can go wrong.** "Any region is allowed"
must never become "assume what I hold is usable". The first attempt at this (section 6) failed for
exactly this reason — it assumed where it should have checked — and no rule imposed on callers
would have saved it.

**Why no ordering rule, but a chromosome rule?** Because the two restrictions are not comparable,
and only one pays for itself.

*Requiring regions to move forward was rejected.* It buys one thing — a caller that walks
backwards by mistake gets a loud error instead of quiet slowness — and the per-cursor tally of
jumps versus reuses shows that anyway. Meanwhile a backward jump costs a seek and a block, which
is what **every** request costs today. Restricting it protects against nothing much and puts error
handling at every call site.

*Requiring a new cursor to change chromosome was taken.* Three things separate it from the above.
The cost it prevents is chromosome-sized, not block-sized: on CRAM a chromosome change means
re-reading the reference bases, hundreds of megabytes. The caller already has that boundary — the
region walk goes chromosome by chromosome — so it is not a new state machine. And **nothing in a
cursor survives a chromosome change anyway**: the kept reads are useless, the bases are useless,
and the only thing that could carry over is the file descriptor, at about 19 µs to reopen. The
restriction states what is already true instead of imposing something new, and it stops the
interface inviting an operation that is catastrophic on CRAM.

A cursor therefore publishes the chromosome it covers, so a caller compares and mints a new cursor
at the boundary. The error exists as a guard against a bug, not as a step in normal control flow.

## 5. Two formats, one set of decisions

BAM and CRAM store their data differently but we ask the same questions of every read once it is
decoded. So the split is: **each format finds and decodes its own reads; the shared part decides
what to do with them.**

**The shared part, written once for both formats:** is this read on the chromosome we asked for,
does it overlap the region, does it start past the end (in which case we can stop, since the file
is sorted), which read group is it, and does it belong to another sample. A read replayed from
memory and a read freshly decoded go through *the same lines*, which is what makes them behave
identically — rather than that being something we claim.

**What differs by format is only finding and unpacking:**

| | BAM | CRAM |
|---|---|---|
| index gives | byte ranges | container entries |
| unpacks | one read at a time from a block | a whole container at a time |

**What is kept is the same on both arms: our reads, already decoded and filtered** (owner,
2026-08-02). A read is returned by about a dozen consecutive regions — the caller widens each one
by 5,000 bases and regions average 390 — so keeping raw records would mean turning that same read
into an `AlignedRead` a dozen times. Keeping our reads means doing it once.

That also removes a problem this section used to have. Keeping *raw* records on the CRAM arm means
keeping a container, and the reader holds exactly one (`container: Option<DecodedContainer>`,
`region_query.rs:597`); a region straddling a container boundary would re-decode the one it
dropped, and real retention would have needed a *set* of containers, each on the order of 10⁴
records. Keeping our reads makes that moot — CRAM's single-container cache stays what it is today,
a within-query optimisation.

**The cost, unmeasured on both sides.** An `AlignedRead` owns its name, CIGAR, sequence and
qualities, so it is larger than the record it came from. Nobody has numbers. **Revisit if memory
bites**; the fallback is to keep raw records and pay the repeated transform.

**One type with two variants, not a trait.** The set of formats is closed, so an enum fits. A
trait would add a fourth type parameter to `PileupGenerator`, which already carries three, and it
would spread to the walker, the generator and the dispatcher; the dynamic-dispatch alternative
would pay a virtual call per read, about a million per run. A third variant holding a fixed list
of reads in memory serves tests and the differential harness. (`ReaderKind`'s doc comment at `open_bam.rs:175-183` is sometimes cited here; it is about a
*different* open question — one pool holding an enum versus production's split `BamFile`/`CramFile`
— not enum-versus-trait for the record source. It is also stale: it says the container cache "will
sit beside this rather than inside it", and the cache shipped **inside** `ReaderHandle` (`:172`).)

## 6. Knowing which reads can be forgotten

**This is the section to read.** Keeping reads is easy; knowing when one may be dropped is the
entire correctness problem. A first attempt got it wrong and **lost 3,830 of 236,081 loci — 1.6 %
— while all 1,471 unit tests passed.**

That attempt used two rules, each sensible on its own, written in two different places:

- *drop reads that end before the current region starts* — they cannot overlap a later region, and
  without this the kept set grows to the whole chromosome;
- *skip a part of the file the reader has already gone past* — its reads are in memory, so there
  is no need to jump back.

The second rule assumed the first had not run. And the assumption underneath both was false:
**moving forward along the chromosome does not mean moving forward through the file.** The index
is a hierarchy of nested bins, so its answer for a *later* stretch of the chromosome can point at
an *earlier* part of the file. By the time the cursor is asked for it, the reader has gone past
and the first rule has already dropped those reads. Nothing can recover them.

**The fix is a comparison of two region starts, and nothing else.**

> **Reuse what is kept only when the new region starts at or after the last one served.
> Otherwise drop everything and jump.**

That single test is sufficient, and the argument is short enough to check. Split the records the
new region needs by where they start:

- **starting at or before the last region's end.** Such a record reaches forward to the new
  region's start, and the new start is at or after the last start — so the record overlaps the
  *last* region too. It was therefore inside the byte ranges the index gave for that region, and
  the scan decoded it before its early stop. **It is kept.**
- **starting after the last region's end.** The scan stopped at the first of these, so the reader
  is sitting on it. **It is read forward**, with no jump.

Every record is in one of those two groups, so nothing is missing. Eviction is the mirror image:
drop a kept record once it **ends before the current region's start**, because every later region
starts at or after this one and cannot overlap it.

**There is no number to tune and no index lookup.** Not a window size, not a capacity, not a
per-format derivation — one comparison of two coordinates.

### Why the obvious alternative — asking the index — is wrong

An earlier draft of this section derived the cut-off from the index instead: for each position the
index records the earliest byte at which a read overlapping it can begin, so keep everything from
there to the reader. It is a natural idea and it fails three ways, all found in review.

1. **The index answer is not monotone; it collapses to zero.** `LinearIndex::min_offset` is
   `self.get(i).copied().unwrap_or_default()`
   (noodles-csi-0.56.0 `…/index/linear_index.rs:15-18`), so past the last populated 16 kb window —
   the tail of every chromosome — it returns byte 0. CSI's `BinnedIndex::min_offset` does the same
   when no ancestor bin carries an offset (`…/index/binned_index.rs:11-27`). "Keep from the
   cut-off" then means "keep everything decoded since the file began", which breaks the memory
   bound this design rests on.
2. **A byte range is not a record set.** The scan seeks from one byte range to the next and skips
   the gaps between them, so bytes inside `[cut-off, reader]` were never decoded. "Between the
   cut-off and here" therefore does not mean "held".
3. **Those two combine into silent loss.** A backward region landing in the same 16 kb window
   gets an equal cut-off, passes the test, and is served without the gap records it needs — the
   first attempt's failure reached by a different road, and §11's backward-request test would
   catch it only after the fact.

The rule above has none of these: it never consults the index, so it cannot inherit the index's
zero; and it never reasons about byte ranges, so it cannot mistake a range for a record set.

**The mistake this also avoids, from elsewhere in ng.** `RawChromReader` widens its buffer
whenever the next request is within 64 KB of the last (`raw_chrom_reader.rs:62`, `:348-349`) and
only an explicit call ever shrinks it (`:374`). Walking a whole chromosome through one shared
reader grows it to about 250 MB on chromosome 1, against a 20 MB baseline. That is the tuned-guess
approach; a coordinate comparison has nothing to tune.

## 7. Cases that must be handled, because they are where it breaks

- **A region abandoned half-way.** The caller stops pulling and moves to another region. The
  cursor is simply repositioned; there is no stream to unwind and nothing to give back. Still
  worth a test, because it leaves the cursor mid-region, which is next to where the first attempt
  lost reads.
- **Two generators, two cursors.** The repeat and the general generators take turns over one
  sample's regions, so each holds its own cursor. The general one's kept reads must therefore span
  the repeat regions sitting between its own — short, against requests that already overlap by
  5,000 bases. Memory doubles:
  about 0.5 MB per cursor at 30× depth, about 5 MB at 300×.
- **A jump forward beyond what is kept.** Drop everything, jump, carry on.
- **A region on another chromosome.** Rejected — make a cursor for that chromosome. Correct code
  never sees this, because it compares against `contig()` first.

## 8. The types

```rust
/// The iterable: fixed at open, shared freely.
/// Read-only once open. No counters, nothing written after construction.
pub struct AlignmentFile { /* path, header, index, read groups, CRAM base registry */ }

impl AlignmentFile {
    /// The chromosomes this file can be read over — the file's own `@SQ` list,
    /// which the checks at open proved agrees with the reference's on names,
    /// lengths and order. Not *identical* to it: an absent digest is a wildcard
    /// on either side, so the digests published here are the file's.
    pub fn contigs(&self) -> &ContigList;

    /// Make a cursor for one chromosome. Opens its own file descriptor and,
    /// for CRAM, takes a share of that chromosome's bases.
    /// Called once per chromosome per worker — never per region.
    pub fn cursor(self: &Arc<Self>, contig: ContigId)
        -> Result<AlignmentCursor, AlignmentFileError>;
}

/// The iterator: one consumer, one chromosome, everything mutable is its own.
/// Everything mutable in the read path lives here. Plain fields — no atomics,
/// no locks — because exactly one consumer touches them.
pub struct AlignmentCursor {
    /* the filter chain (which owns the record reader), the reads kept for the
       next region, the last region's start, the chromosome, and for CRAM the
       reference-base handle. Field-by-field in the arch doc §1.2. */
}

impl AlignmentCursor {
    /// The chromosome this cursor covers. Compare against it before asking, and
    /// the error below never fires.
    pub fn contig(&self) -> ContigId;

    /// Point the cursor at `region`. Any region on this cursor's chromosome, in
    /// any order — always correct.
    ///
    /// Costs nothing when the region is already inside what the cursor holds;
    /// otherwise it jumps and re-reads, as the code does today for every region.
    pub fn move_to_region(&mut self, region: GenomeRegion) -> Result<(), CursorError>;

    /// The next read of the current region, in position order. `None` once the
    /// region is exhausted; `move_to_region` starts the next one.
    ///
    /// Shaped as `Iterator::next` so the cursor can implement `Iterator` if a
    /// consumer wants adapters — the walker peeks, for instance.
    pub fn next_read(&mut self) -> Option<Result<AlignedRead, CursorError>>;
}

/// Where records come from: finds them with the index, unpacks them, keeps the
/// recent ones, and hands over the next one on demand. One variant per format,
/// a closed set.
enum RecordReader {
    Bam(BamRecordReader),
    Cram(CramRecordReader),
    /// A fixed list of records. Permanent, not test-only: it gives the tests and
    /// the differential harness a reader with no file behind it.
    InMemory(Vec<RecordBuf>),
}

/// No ordering errors — within a chromosome there is no ordering rule.
pub enum CursorError {
    /// A region on another chromosome. Make a cursor for it; this one is
    /// unaffected and still good for its own. A guard against a caller bug,
    /// not a step in normal control flow.
    WrongChromosome { cursor: ContigId, requested: ContigId },
    Io(std::io::Error),
}
```

There is no stream type. `RegionReads`, and the borrow-and-return dance around it, are deleted.

## 9. What existing code becomes

| what | today | after |
|---|---|---|
| deciding what to do with a decoded read | `BamRegionSource::read_next` steps 4–5 | moved into the cursor, one copy for both formats |
| asking the BAM index where to look | `BamRegionSource::plan` | into the BAM half, unchanged — measured at 0.43 µs per region, so not a cost |
| decoding and keeping a CRAM container | `DecodedContainer`, `CramRegionSource` | becomes the CRAM half; the keeping stops being CRAM-only |
| the pool of readers and its buffers | `ReaderHandle`, `BorrowedReader` | **deleted** |
| checking the file's chromosomes against the reference | the open-time checks | unchanged; now also published as `contigs()` |
| the correctness oracle | `t5_the_indexed_query_returns_exactly_what_a_linear_scan_returns` | **extended** — see section 11 |

## 10. Memory, errors, threads

**Memory — arithmetic, not measured.** No RSS figure exists for the failed implementation, so
every number here is derived and should be read as such. A cursor keeps about as many reads as two
consecutive requests overlap by. Under §6's rule that span is the caller's widening, ~5,000 bases,
because eviction is by coordinate; an earlier draft derived the same figure while using an
index-rounded cut-off, where the true span would have been a 16 kb tile plus the widening plus a
read length — about four times larger. The number is right now for a different reason than it was
first given. That gives roughly 0.5 MB at 30× depth, 5 MB at
300×, against a whole-genome walk that peaks at 30 MB. What is kept is bounded by how far apart
two consecutive regions are, never by chromosome length. Keeping reads was measured on its own in the first attempt at
**+2.8 % of wall time** — the cost of copying each read — against the 43 % that decompression
occupies.

**Errors.** Two new cases, both meaning "make a new cursor". Input/output errors keep their
current shape. The cursor is consumed when it rejects a region, so a half-valid cursor cannot
exist.

**Threads.** Nothing inside a cursor is locked, because nothing is shared. `AlignmentFile` can be
used from many threads; a cursor cannot be, in the same way and for the same reason as the
windowed reference reader — an open file position belongs to one consumer.

### The CRAM reference bases — unchanged, and here is why the redesign is deferred

A cursor covers one chromosome, so it calls `bases_for_contig` **once at construction** and holds
the handle for its life. For every caller that exists today — all single-threaded — the accessor's
existing one-chromosome bound is then exactly right: one chromosome resident, cleared when the
walk moves on.

**An earlier draft replaced that with a registry of weak per-chromosome references, and it is now
deferred (§12).** The design is sound and the motivation is real — workers on different
chromosomes evicting each other's bases — but the motivation is the fan-out, *which does not
exist*. Meanwhile every measurement behind this document is BAM, and **a BAM run never opens the
FASTA at all** (perf review L4 prices this at "zero cost today"). Rebuilding a piece of the
cohort-memory design on the strength of a workload nobody has run is the kind of unmeasured work
this document is otherwise arguing against.

## 11. How we will know it works

**The unit test suite is not the bar.** All 1,471 of them passed while the first attempt lost
3,830 loci, because every existing read-path test drives a *single* region query. Consecutive
queries through one cursor are the untested surface, and they are the entire feature.

**Every threshold below is an absolute count from one tandem-repeat-targeted fixture** (§1). They
are regression anchors against themselves, not properties of the generator: a different alignment
file changes all of them.

1. **Identical output on real data.** `ng_generic_walk_probe` on HG002 30×, chromosome 21, prints
   exactly `loci=236081 observations=251786 reads_admitted=54709`; chromosome 1 prints
   `loci=1541788 observations=1647161`. This is the check that caught the defect.
2. **Identical output at every region size** — 400 bases, 10 kb, 100 kb, and one region for the
   whole chromosome. The first attempt first diverged at 100 kb.
3. **The oracle, extended to a sequence.** Today it compares one indexed query against a plain
   scan of the file. It must compare *a run of ascending queries through one cursor* against the
   same scan. This is the missing test that let the defect through.
4. **One test per case in section 7** — abandoned stream, jump beyond what is kept, two cursors
   interleaved, a backward request, and a request on another chromosome. The last two must return
   the same reads a plain scan of the file would, not an error.
5. **The saving is asserted, not assumed.** Count reads decoded, replayed from memory, and jumps.
   On chromosome 21 the number decoded must approach the true read count — **and the threshold
   must name its mode**, because the two figures in circulation are from different runs: 34,633 is
   the *typed-region* walk (the failed implementation's own counters), 35,228 the *whole-contig*
   probe mode (perf review H4). Assert against the mode the test runs in.
6. **Break it on purpose.** Remove the start comparison, the eviction test, and the
   already-passed-it skip in turn; each must fail a test. This project keeps finding tests that cannot fail, and a cache is where
   one hides best.

## 12. Left for later, with a home

- **Answering a batch of sorted regions in one sweep.** `alignment_file.md` §3.3 raised it; this
  design is what makes it possible, not the thing itself. Home: a follow-up once the cursor works.
- **One cursor shared by both generators.** Their regions interleave but still move forward
  overall, so one cursor could serve both and halve the memory. Not now — it ties two generators'
  lifetimes together. Home: revisit when the parallel fan-out lands.
- **The per-region reference-accessor factory** (perf review finding L2). The cursor is the
  natural owner of the reference accessor its mismatch filter needs. Home: fold in during
  implementation if it is free; otherwise its own change.

## 13. Decisions made, questions open

**Agreed with the owner, 2026-08-01.**

1. **The caller keeps the cursor and pulls reads into a buffer it owns** — `move_to_region` then
   `next_read()`. Two alternatives lost. *The cursor calling the caller back for each
   read* would force the walker to store every locus it produces, which
   [locus_generation_pileup.md](locus_generation_pileup.md) §7 forbids. *Returning a stream object
   per region* was the draft this document carried for several revisions; it owned or borrowed the
   cursor, and every difficulty about lifetimes, giving the cursor back, and callers forgetting to
   came from that object existing. Deleting it removed all of them — spec §3.
2. **No pool, no lock; a worker takes a cursor once.** This reverses `alignment_file.md` §3.3,
   whose stated reason was future parallelism. Per-worker cursors serve that reason better, with
   nothing shared at all.
3. **One type with a variant per format**, rather than a trait — the format set is closed, and a
   trait would add a fourth type parameter to an already heavily-parameterised generator.
4. **Any region within the cursor's chromosome, in any order; a new cursor to change chromosome.**
   Two restrictions were weighed separately and only one was kept, which is the useful part of
   the record.

   *Ordering was not restricted.* A backward jump costs a seek and a block — what every request
   costs today — so forbidding it would protect against almost nothing while putting error
   handling at every call site.

   *Changing chromosome was restricted.* The cost is chromosome-sized rather than block-sized: on
   CRAM it means re-reading the reference bases, hundreds of megabytes. The caller already has
   that boundary, since the region walk goes chromosome by chromosome. And nothing in a cursor
   survives the change anyway — kept reads and bases are both useless — so the rule states what is
   already true rather than imposing something. It also removes the swap logic and the
   intermediate state that free movement would need inside the cursor.
5. **The reuse check must be total.** This is what the first attempt got wrong (section 6): it
   assumed rather than checked. No rule imposed on callers would have saved it, which is why the
   rule now lives in the cursor.
6. **Replace the current region readers outright**, with no parallel path behind a flag.

**Resolved by reading the code, 2026-08-02.**

7. **The forget rule is one comparison of region starts, and both formats share it** — no index
   lookup on either side. This replaces an earlier draft that derived a byte cut-off from the
   index per format; §6 records the three ways that failed, all found in adversarial review. The
   consequence worth noting: a question this document once carried — whether CRAM needed a cut-off
   derivation of its own — no longer exists, and neither does the unwritten CSI arm the index
   version would have required (`AlignmentIndex::BamCsi` is a supported input,
   `region_query.rs:130`, and CSI has no linear index).

   What each format owns is *what* it keeps — and on the CRAM arm that is not settled: the reader
   holds one container, not a set, so it cannot retain across a container boundary today (§5).
   Whether retention belongs per file at all is the open question below.

**Also rejected, recorded so it is not re-proposed.**

- **The caller keeps the reads it already pulled, and asks only for the new suffix.** Because the
  widening is one-sided — `[region.start, region.end + max_record_span]`,
  `generator.rs:839-848` — consecutive queries are nested forward, so a caller could retain the
  reads whose end reaches the next region and query only what is new. It is **provably complete,
  needs no index reasoning and changes nothing in the read layer** — the one caller-side shape
  that is actually correct. It loses on cost, not on soundness: the walker takes `AlignedRead` by
  value, so each retained read is cloned once per region it is replayed into — on the order of
  10⁶ clones of the pipeline's largest object, against the cursor's ~35 k record copies. And the
  per-region drop tallies would change meaning, each read counted once instead of ~13 times,
  which breaks "identical output" as the acceptance test.

**Two things the one-sided widening banks** (`generator.rs:839-848`): because starts *and* ends
both ascend, §6's comparison is as cheap at the call site as inside; and N1 already closed
narrowing the widening — bit-identical output for 7 % of wall — so there is no tension between
the two knobs to resolve.

**Still open.**

- **Do not assume "kept from the cut-off" and "kept what was scanned" name the same range.**
  `optimize_chunks` drops chunks ending at or before the cut-off but does not truncate the start
  of the one that survives (noodles-csi-0.56.0 `binning_index.rs:134-140`), so a scan can begin
  before it. Harmless in itself; it would make an assertion that the two are equal fail.

- **Must a replayed record be staged through a scratch buffer before it is decoded?** Reads reach
  the caller owned — as `AlignedRead`, which owns its bases, qualities and CIGAR — and that does
  not change. But the replay path copies a kept record into a scratch buffer first, only because
  `read_next` fills a buffer the caller supplies. Decoding straight from the kept record would
  remove that copy, measured at **+2.8 % of wall time**, with the borrow living inside a single
  `next()` call and never reaching a struct field.

  Leaning: keep the copy for the first version. Settle it by measuring the finished cursor end to
  end; against the ~40 % the first attempt showed, 2.8 % does not justify reshaping the record
  source before there is a working baseline.

  **Handing borrowed reads to the caller is *not* the alternative, and should not be revisited.**
  Rust would make it safe — a live borrowed read makes the cursor immutably borrowed, so the next
  region cannot be requested and the cache cannot be wiped underneath it. That is the guarantee,
  and it is too strong: the walker holds reads while advancing across positions, and a stream with
  a lifetime cannot be stored in a field beside the cursor that lends it. `Arc<RecordBuf>` avoids
  the lifetime but pays an allocation and an atomic per read, which is worse than the copy.
