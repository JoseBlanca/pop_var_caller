# ng — how a run is driven: three iterators, two modes

*Status: design spec, 2026-08-16. No code yet — this settles the design. It replaces the
2026-08-14 draft, which drove a run through three public size choices — the walk unit, the psp
block, the calling range — and the machinery that kept them apart. The public shape is now three
objects, each an iterator, and every division of work is internal to one of them. One consequence
of the new shape is argued in §6.3 and **flagged there for the owner to rule on: the psp header's
block-boundary digest is dropped**, and its `writer_version` field with it. Companion architecture
doc: [`../arch/run_streaming.md`](../arch/run_streaming.md) (the types and interfaces). Reads on
[`locus_generation.md`](locus_generation.md) (what an observation is),
[`typed_regions.md`](typed_regions.md) (what a segment is), and
[`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) (the census — the
per-sample evidence the parameters fit reads).*

---

## 1. What this is

ng turns alignment files into a VCF. Between those two points sit a walk over the reads, a fit of
the model parameters, and the calling itself. This document settles **how a run is driven**: the
three objects a run is made of, the one calling function two of them share, how the loop inside
each object is parallelised, and what bounds the memory held at once.

The whole public shape is **three objects, each an iterator**:

- **`AlignedFilesVariantCaller`** — every sample's alignment files plus the model parameters in;
  variants out, in genome order, ready for the VCF writer. This is **direct mode**, entire.
- **`PspVariantCaller`** — every sample's psp file plus the parameters in; the same variants out.
  This is **psp mode's calling stage**.
- **`SampleObservationGatherer`** — *one* sample's alignment files in; that sample's observations
  out, in genome order, ready to be written to a psp — and, when the iterator is exhausted, the
  sample's census (`finish()`). psp mode's walk stage is a loop over samples with one of these
  each.

Everything else — how the genome is divided among workers, where a psp's compressed blocks fall,
how out-of-order workers become an in-order output — is internal to one of the three. No caller of
these objects ever names a work unit, a block, or a range.

### 1.1 Goals

1. **One calling function, whichever mode the run is in.** The code that merges samples and calls
   genotypes is written once; the two caller objects differ only in where each sample's
   observations come from. This is a goal rather than an observation because it is easy to lose:
   production grew a second, separate driver for its no-file path and deleted it on 2026-06-01
   rather than keep maintaining two.
2. **Degrade across the committed range** — one sample to several thousand, a few reads a position
   to several hundred (`CLAUDE.md`, *What this caller has to work on*). Every memory bound in §7
   is a formula in the sample count, so the large end has an answer rather than an omission.
3. **Parallelism inside one sample as well as across samples.** Production's single-sample pileup
   reaches 1.81× at four threads and gets *worse* at eight (59.2 s → 32.7 s → 34.4 s), because
   only its BAQ stage threads while reading, walking and writing stay serial — about 40% of that
   run cannot be sped up
   ([`pileup_thread_scaling_2026-06-11.md`](../../reports/pileup_thread_scaling_2026-06-11.md)).
   A single sample must be able to use the machine.
4. **The output must not depend on how the work was divided.** The same VCF at one thread and at
   sixteen; the same VCF from direct mode and from psp mode with the parameters held fixed. This
   is the oracle for everything below (§12).

### 1.2 Non-goals, and what this document does not do

- **It does not define the psp file's encoding** — byte layout, compression, block sizing,
  checksums, format versioning, the index. It fixes only the contract the file's reader must
  answer (§3.3), what the header must record (§6), and the resident-state budget an open file
  must meet (§7.2). The encoding is deferred with a home (§10).
- **It does not define the cohort merge's reconciliation** — how observations whose spans differ
  between samples, because each sample's generic observations are cut from its own data, become
  one cohort observation. Deferred with a home (§10). This document fixes that the merge keys on
  coordinates and that it streams (§3.2).
- **It does not define the caller** (steps 6–13 of the proposal) or the parameters file's format.
  Both sit downstream of the objects named here; the parameters file is deferred with a home and
  is on direct mode's critical path (§10).
- **It does not decide what a BED does to a segment.** The typed-region generator already decided
  that (§4.2); this design takes the segments as given.
- **It does not re-open where the census lives.** One file per sample, beside its psp, never
  inside it — decided 2026-08-13
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1). This
  document only says when that file is written (§5.2).

### 1.3 Vocabulary

- **analysed regions** — the set of genome stretches the run was asked to analyse. Without a BED
  it is one region per contig; with a BED it is the BED's intervals (§4.1).
- **segment** — one stretch of the reference as the typed-region generator cuts it: an STR tract,
  a bundle, a satellite, or the generic stretch between them (§4.2). Nothing downstream ever cuts
  one (§4.3).
- **observation** — a stretch of genome plus what one sample's reads showed there: a
  `SampleLocusObservations` ([`locus_generation.md`](locus_generation.md) §3). The locus-generation
  documents call this a *locus*; this document says *observation* throughout, because what matters
  here is the item flowing between objects, and every observation is minted inside one segment.
- **cohort observation** — every sample's observations over one stretch of genome, brought
  together by the merge. How overlapping spans are reconciled inside one is the deferred merge
  spec's question (§10); that it is the unit the caller consumes is this document's.
- **census** — a sample's evidence at a scattered subset of the analysed positions, chosen in
  advance and identical in every sample, which is all the parameters fit reads. About one
  analysed position in 390 on tomato.
- **psp** — the per-sample file holding all of one sample's observations, in genome order.
- **block** — the psp's unit of compression: a run of observations compressed as one payload and
  decoded as one payload. Named here only so this document can say where it lives: **inside the
  psp writer and reader, and nowhere above them** (§3.3, §6.3).
- **source** — the object that answers one question: *one sample's observations in one segment, as
  an iterator, in coordinate order.* Two exist — a walker over alignment files, a psp reader
  (§3.3).
- **look-ahead** — how many segments an object may hold in flight beyond the one it must yield
  next (§3.5). The memory knob.

---

## 2. The two modes, and the build order

**A run is in one of two modes, and what selects them is where the evidence comes from** — the
reads in the alignment files, or the psp files an earlier run wrote from them.

**Direct mode — evidence from the alignment files, nothing stored.** One object,
`AlignedFilesVariantCaller`, walks every sample together and yields variants; no intermediate file
is written.

**It requires the model parameters up front, and that is a precondition rather than what defines
the mode.** The parameters have to be fitted from the whole cohort's evidence before a single
variant can be called, so a run that has to fit them cannot call anything on its first pass over
the reads — it would need a second pass, and decoding every read twice is what the psp exists to
avoid. So *evidence from the alignments* and *parameters in hand* go together: without the second,
the first is not available, and the run belongs in psp mode.

Where the parameters come from is not ng's business — a previous fit, a published set, a service
whose chemistry does not change from week to week — which is how freebayes and GATK are normally
run. The consequence for us: the parameters file is a **user-facing input format** — documented,
readable, hand-editable — and its spec (§10) is on this mode's critical path.

**psp mode — evidence from the psp files, in three uncoupled stages.** Each may be a separate
invocation:

1. **the walk** — each sample is walked on its own by a `SampleObservationGatherer`; the
   alignments are read once; the stage writes two files per sample, its psp and its census file.
2. **the fit** — the model parameters are fitted from the cohort's census files. The fit never
   opens a psp or an alignment file.
3. **the calling** — a `PspVariantCaller` reads the psps, merges the samples, calls, and yields
   the variants.

The stages couple only through the two files, so a sample can be added to a cohort later without
re-walking the others, and a failed sample is one sample to re-run. The alignments are read
exactly once, which is the reason the psp exists: re-walking costs a full decode of every read,
where reading back a stretch of psp costs a seek and a decompression.

**Build order — decided: direct mode first.** Three reasons. It needs no file format, so it
freezes nothing while the caller's own steps are unwritten. It is the shortest path to ng emitting
a VCF. And once psp mode exists, direct mode is its oracle: the same cohort with the same
parameters must give the same VCF from both modes, which is also the test that the psp carries
everything calling needs (§12).

---

## 3. The calling function, the merge, and the shape all three objects share

### 3.1 One segment at a time, through one function

Both variant callers — the one reading alignment files and the one reading psp files — drive one
function. Given one segment and one source per sample, it merges the
samples' observations and calls variants over that segment:

```
# a "source" answers one question: this sample's observations in this segment,
# as an iterator, in coordinate order.
#   the walker      – reads the sample's alignment files
#   the psp reader  – decodes whichever blocks overlap, keeping the one it is in

call_vars_in_segment(segment, sources):
    per_sample = [source.observations_in(segment) for source in sources]
    for cohort_observation in k_way_merge(per_sample):
        yield each variant of call_vars_from_observation(cohort_observation)
```

and each caller's `next()` is, in its serial skeleton:

```
next():
    while pending is empty:
        segment = segments.next()          # no more segments → iteration ends
        pending = collect(call_vars_in_segment(segment, sources))
    return pending.pop_front()
```

`call_vars_in_segment` takes **one** segment, deliberately. The loop over segments is what gets
parallelised, so it belongs where the in-flight bookkeeping is — inside the iterator — and not
inside the calling function, which stays a straight-line pipeline anyone can read or test on one
segment. The two callers keep their own `next()` — they differ in what they hold and how they
parallelise (§5) — while `k_way_merge` and `call_vars_from_observation` are written once. That is
the whole of goal 1: nothing inside the merge or the caller can tell a walker from a file reader.

The serial skeleton above is the **contract** — what the iterator yields and in what order. §3.5
says how the same yield sequence is produced in parallel.

### 3.2 The merge is streaming

`k_way_merge` consumes the per-sample iterators and yields one cohort observation at a time. Its
resident set is its **frontier** — the head observation of each sample's iterator, so about one
observation per sample — never a whole segment's observations per sample. Nothing accumulates a
segment before merging it; the merge pulls, the sources produce lazily, and a cohort observation is
dropped the moment the caller has consumed it.

The repo already has this shape in its read layer: `MergedRegionReads` is an argmin k-way merge
over per-file read streams, holding one head per stream with the sort keys beside the heads —
[`sample_reads.md`](../arch/sample_reads.md) §4 describes it, and it is the model for
`k_way_merge`. The differences are the item (observations, not reads) and the yield (a cohort
observation groups every sample's entries at one stretch of genome rather than interleaving
them); where two samples' generic observations overlap without coinciding, the frontier may
briefly hold a few observations of one sample while spans are reconciled — the deferred merge
spec (§10) owns that reconciliation and must confirm the frontier stays bounded by it.

**The consequence: segment size is not a memory knob.** A segment's cost while being called is its
merge frontier plus the sources' own working state, not its length. Segments can therefore be
exactly the segments (§4.4), with no machinery to group them toward a size target — the previous
draft's walk units and calling ranges existed to hit size targets that no longer control anything.

**The merge keys on coordinates — decided, carried forward.** Generic observations are split out
of each sample's own data, so they do not line up between samples. A merge keyed on anything
positional-within-a-file — a block number, an index within a decoded payload — would be wrong;
coordinates are the only key all producers share.

### 3.3 Almost every position is never built, and the saving lands in two different places

**A psp holds every locus, because no sample knows which positions another sample varies at.** What
is skipped is not the storing and not the visiting — it is the *building*. A position where every
sample's reads matched the reference produces no cohort observation and no call at all, and an
earlier measurement of exactly this question on the cohort reader materialised **28,718 loci
instead of 2.83 million**, about one position in a hundred.

**The decision is a cohort one and cannot be taken per sample.** A position where *this* sample saw
only reference reads may be a variant site in another, and there this sample's evidence is still
needed — its depth is what separates a confident homozygous reference from no coverage, and it is
part of the allele frequency every genotype at that locus is weighted by. So no sample's quiet
position is dropped on its own account; what is dropped is a position no sample varied at.

**Where the saving is taken differs by mode, because the expensive thing differs.**

- **Direct mode: at the merge.** The walker has already read the reads, so its observations exist
  whatever happens next. The saving is downstream — where the cohort's non-reference evidence at a
  position is nothing, `k_way_merge` builds no cohort observation and `call_vars_from_observation`
  is never entered.
- **psp mode: during decompression**, which is the production shape and worth copying rather than
  re-inventing. `TwoPhaseSegment`
  ([`var_calling/sample_reader.rs:698-712`](../../../../src/var_calling/sample_reader.rs)) decodes
  a block's **light** columns for every row — the positions, the reference span, and the sum of
  non-reference observations — while leaving the heavy columns compressed. The cohort fold sums the
  light column across samples and decides which rows are variable; only then does
  `set_variable_rows` ([`:789`](../../../../src/var_calling/sample_reader.rs)) inflate the heavy
  columns of the kept rows, in every sample. A quiet sample at a kept position is therefore fully
  present, and a position nobody varied at costs one integer per sample.

**This makes a source's answer two-phase rather than one call**, and the same two phases in both
modes:

```
# phase 1 — cheap, every position in the segment
source.light_column(segment)     -> per position: depth, and whether any read was non-reference

#   the merge folds phase 1 across the samples and decides which positions are variable

# phase 2 — expensive, only at the positions the fold kept
source.build(kept_positions)    -> the full observations there
```

The walker's phase 1 costs a walk and its phase 2 is free, because it is holding the reads already;
the psp reader's phase 1 is a light decode and its phase 2 is the deferred inflation. **One
interface, two cost profiles** — and the merge, which cannot tell them apart, is where the decision
is made in both.

*The production path also ships its own oracle, which is worth carrying: an eager whole-segment
decode used only by tests, as the byte-identity check the two-phase path is measured against
([`sample_reader.rs:20-26`](../../../../src/var_calling/sample_reader.rs)).*

### 3.4 The two sources

A source answers *one sample's observations in one segment, in coordinate order*. It is the single
place the two callers differ.

- **The walker** wraps one sample's open alignment files. Asked for a segment, it points a read
  cursor at the segment's span, runs the locus generators over it, and yields the observations as
  they are minted. It is the same machinery the gatherer uses (§5.2); direct mode is that
  machinery driven by the calling loop instead of by a file-writing loop.
- **The psp reader** wraps one sample's open psp. **It serves segments, not blocks**: asked for a
  segment, it decides internally which blocks overlap, decodes them one at a time, and yields the
  observations inside the segment. It keeps the block it is currently in, so successive segments
  inside the same block cost no further decode. The block never surfaces above the reader — no
  caller mentions one, the two of them loop over the same segments, and there is no "line the segment
  up with a block" problem: any segment is servable, and the segments they loop over are
  simply the common case its cache is warmest for.

Both sources are cursors: asked for ascending segments they stream forward; a backward jump is
legal and costs a seek (§8). A source serves one consumer at a time, so a parallel loop gives
each in-flight segment its own source per sample, over the sample's one open file (§5, §7).

### 3.5 Out-of-order work, in-order output: bounded look-ahead

All three objects share one internal skeleton. The segments — the segments, in genome order — are
handed to a pool of workers; workers finish out of order; the iterator yields in genome order, so
a finished segment waits until every earlier segment has been yielded. The number of segments in
flight — being worked on, or finished and waiting — is bounded: the **look-ahead**. When the
segment at the yield frontier is slow, up to `look-ahead` later segments may complete and be held;
then the pool idles rather than pull further ahead.

With a streaming merge, the look-ahead is the objects' only real memory knob: peak resident is
`look-ahead × one segment's in-flight cost`, plus per-sample open-file state (§7.1). One knob per
object, spent the same way in all three; at look-ahead 1 each object degrades to its serial
skeleton and its minimum memory.

What differs per object is only what a segment's in-flight cost *is* — a walking working set, a
decoded block per sample, or one segment's observations — and §5 prices each.

---

## 4. The ground the loop runs over

### 4.1 A run is given a set of segments to analyse

A run analyses a **set of segments**, held in `GenomeRegions`
([`src/ng/region_typing/mod.rs:77`](../../../../src/ng/region_typing/mod.rs)), which wraps the
same `RegionSet` production's `--regions` runs on. Its own doc states the rule — *"'Whole genome'
is not a special case — it is the region set whose every region covers an entire contig"* — and
`whole_contigs` ([`:87`](../../../../src/ng/region_typing/mod.rs)) is the default rather than a
bypass; a BED loads through `from_bed_path`
([`:100`](../../../../src/ng/region_typing/mod.rs)). No code in this design tests "is this the
whole genome".

### 4.2 Segments, and what a BED edge does to one

The typed-region generator cuts the reference into segments — `TypedRegion`
([`src/ng/region_typing/mod.rs:144`](../../../../src/ng/region_typing/mod.rs)) with a
`RegionKind` ([`:168`](../../../../src/ng/region_typing/mod.rs)): an STR tract, a bundle, a
satellite, or the generic stretch between them. The segmentation is a function of the reference,
the repeat catalog, the routing criteria and the analysed regions — **no sample's reads** — so it
is identical in every sample of a run.

What a BED edge does to a segment is already the generator's decision, not this design's: a
finding (STR tract, bundle, satellite) is emitted whole even where it crosses the edge, and only
a `Generic` stretch is clipped to the requested span (`clips_at_a_bed_edge`,
[`src/ng/region_typing/mod.rs:471`](../../../../src/ng/region_typing/mod.rs), emission rule
[`:482-488`](../../../../src/ng/region_typing/mod.rs);
[`typed_regions.md`](typed_regions.md) §2.5: the BED chooses what you are shown, never what
things are). So a BED edge is a segment boundary before this design sees it.

### 4.3 No boundary splits a segment

**Owner's rule, 2026-08-09: a segment is never cut.** Every observation is minted inside one
segment, so a loop whose segments are whole segments never splits an observation — which is the
whole correctness argument for calling segments independently: nothing a worker needs sits on the
far side of its segment's edge.

The rule is measured, not aesthetic. A test that chopped each generic stretch into thirds on
`benchmarks/tomato1/crams/SRR7279481.p1.bench.cram` lost 17 positions out of 7,429,336: one read
carried a 91-base deletion spanning `SL4.0ch01:32,931,518–32,931,608`, a cut landed 74 bases
inside it, and the part of the deletion past the cut was emitted by no segment at all.

**What the independence argument excludes:** treating evidence in two adjacent segments as one
variant — say, a long deletion whose reads touch a generic stretch on each side of a repeat
tract. If that is ever wanted, it must be a pass over the emitted records, never a coupling
between in-flight segments. The residue is measured and small: three sites on one tomato sample
where a read's deletion starts in one generic segment, jumps a repeat tract, and ends in a later
one — a few tens of positions in 7.4 million.

### 4.4 The loop unit is one segment — decided

Each caller's loop, and the gatherer's, advances one segment per iteration. No grouping of
segments into larger hand-out units exists in this design.

The previous draft argued a lone segment was too small a unit of work, from a measured average of
391 bases (613,682 generic segments on human chromosome 1). **Owner's ruling: that average does
not apply, and the mechanism built on it goes.** The 391-base figure was taken at the *catalog's*
admission floor, `CATALOG_MIN_COPIES = [5, 5, 4, 4, 4, 3]`
([`src/ng/repeat_catalog/criteria.rs:16`](../../../../src/ng/repeat_catalog/criteria.rs)), which
admits a five-base homopolymer — in random sequence a run of five or more identical bases begins
about one position in 341, which matches the average. At the floor a caller actually routes on —
the copy number at which a repeat starts to stutter, `[8, 6, 6, 6, 5, 4]`, measured over 2,457
tomato libraries
([`src/ng/region_typing/segment_criteria.rs:402-414`](../../../../src/ng/region_typing/segment_criteria.rs))
— the same arithmetic gives one homopolymer run every ~22,000 bases, so a generic segment at the
routing floor is probably kilobases long. Probably: it is unmeasured, and open question 1 (§11)
names the measurement — a filter over the stored catalog file, not a genome scan.

If segments at the routing floor turn out short enough that per-segment overhead shows up in a
profile, the fix is internal to the loop — for instance, a worker keeping its cursors and
generators across the consecutive segments it processes (§5.1) — and changeable without touching
any interface, because no interface names a unit larger than a segment.

---

## 5. The three objects

### 5.1 `AlignedFilesVariantCaller` — direct mode

**Constructed from:** every sample's alignment files, the segmentation's ingredients (reference,
catalog, routing criteria, analysed regions), and the model parameters. **Yields:** variants, in
genome order.

**What it holds for the whole run:** one open `SampleReads` per sample — about 11 to 15 MiB of
live heap per open alignment file (measured slope 12.0 MiB per file, 10.8 under a different
allocator — [`examples/dhat_ng_open_files.rs`](../../../../examples/dhat_ng_open_files.rs)) —
plus the shared read-only state (§9).

**How it parallelises:** the segment loop of §3.5. One in-flight segment is one task: every
sample's walker over that segment, the merge, the calls. Within a segment the merge is one
consumer, so the parallelism is across segments — which serves both ends of the cohort range: at
one sample only the segment axis exists, and at a thousand samples the segments are still millions.
Each in-flight segment's task owns its own per-sample cursors and generators (the traps of §8 make
sharing them wrong); a worker that processes consecutive segments keeps each cursor's movement
forward-only, the fast path (§8).

**What bounds it:** `look-ahead × samples × one segment's walking set` +
`samples × 11–15 MiB` of open files. A segment's walking set is the active reads and generator
state at the merge frontier, not the segment's length (§3.2).

**Where the mode stops fitting.** Direct mode is for the run that has the memory, already has the
parameters, and is not coming back to these samples in another cohort. Needing every sample open
at once is what calling *is*, not a defect of the mode — but the open files alone price it: 0.9 GB
at 63 samples, 15 GB at a thousand, 44 GB at three thousand, before a read is decoded. When that
bill does not fit, the answer is psp mode, and the failure must be a message saying so, not an
allocator kill. Where exactly the crossover sits is open question 6 (§11).

### 5.2 `SampleObservationGatherer` — psp mode's walk

**Constructed from:** one sample's alignment files, the segmentation's ingredients, and the
census configuration. **Yields:** the sample's observations, in genome order, ready to be written
to a psp. **`finish()`** — called after the iterator is exhausted — hands over the census the
gatherer accumulated: concretely, `CensusWriter::finish()` returning the sample's
`SampleCensusEvidence`
([`src/ng/parameter_estimation/joint/census.rs:1806,2252,1349`](../../../../src/ng/parameter_estimation/joint/census.rs)).

**The walk stage is a loop over samples**, one gatherer each, default one at a time:

```
for each sample:
    gatherer = SampleObservationGatherer::new(sample's files, ...)
    psp writer consumes the gatherer            # blocks, header, trailer: the writer's business
    census file written from gatherer.finish()  # write_census — census_file.rs:195
```

One sample at a time keeps this stage's peak memory independent of the cohort size — the property
psp mode exists for — and holds one alignment file open. Samples-in-flight is a knob, raised only
when one sample's segments cannot fill the pool (a small genome, a very high worker count); each
extra sample in flight costs another open file at 11–15 MiB plus its workers' working sets. Never
a thread per sample, which at a thousand samples is a thousand threads.

**The census is fed from inside the gatherer, on the same stream it yields.** Every observation
the gatherer yields passes the census accumulator (`CensusWriter::add_locus`,
[`census.rs:1965`](../../../../src/ng/parameter_estimation/joint/census.rs)) at the ordered yield
point, and every segment the loop completes is marked walked (`mark_walked`,
[`census.rs:1984`](../../../../src/ng/parameter_estimation/joint/census.rs)) whether or not it
produced observations. Two consequences, each closing a defect the previous draft closed by
contract:

- **The psp and the census cannot see different evidence.** There is one stream: what the
  iterator yields is what the census counted. No delivery discipline between two sinks is needed,
  because there are no longer two sinks — the second consumer lives inside the producer.
- **Quiet ground cannot be lost.** A walked stretch with no reads yields nothing, and to the
  census a segment never visited and a segment with no coverage would look the same — only the
  first is a bug. The loop that knows a segment is finished and the census that must be told are
  now the same object, so no worker can forget to report an empty segment; there is no call to
  forget.

The census accumulator is held for the whole sample — about 6 MB per read group at a
two-million-position census — because two of its counts are per-stratum sums the walk accumulates
as it goes; that is why the census file is written at the end. The psp and the census are two
files rather than two parts of one because the census is a rebuildable cache of the psp: its
positions depend on a per-run budget and a seed, so "rebuild it" must mean deleting a small file,
not rewriting a large one — decided 2026-08-13
([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §6.1).

**How it parallelises: over the segments of the one sample it was constructed from**, never across
samples — the sample loop above is serial by default and a gatherer only ever sees one sample's
files. Within that sample it is the same skeleton (§3.5) — its segments to workers, each with its own
cursor, reference accessor and generators, observations drained in genome order. Read-filter
tallies live in each worker's cursor and are summed when the gatherer finishes, or drop rates
under-report by the worker count (§8). Whether a within-sample walk scales is the open question
that can invalidate this schedule — production's equivalent tops out at 1.81× on four threads
(goal 3); ng starts from a better position because its per-segment work is independent, but that
is a leaning, not a measurement (open question 3, §11). If it fails, the fallback is several
samples at once with few workers each; nothing else in this document changes.

**What bounds it:** `look-ahead × one segment's observations and working set` + the census
accumulator. Nothing crosses samples.

### 5.3 `PspVariantCaller` — psp mode's calling

**Constructed from:** every sample's psp, the segmentation's ingredients, and the parameters.
Construction opens every file, reads every header, and runs the checks of §6.2 before any block
is decoded; the analysed regions come from the headers, not from a flag — the files know what
ground they cover. **Yields:** variants, in genome order.

**How it parallelises:** the segment loop of §3.5, one source per sample per in-flight segment.
Each source is a psp reader's cursor over the sample's one open file; the open file's resident
state — required to be tens of kilobytes (§7.2) — is paid once per sample, and each cursor holds
its own current decoded block. Successive segments in one worker land in the block its cursor
already holds far more often than not, which is what makes segment-grained segments cheap here;
two workers whose segments share a block each decode it, a duplication bounded by the look-ahead
and paid in time, not correctness.

**What bounds it:** `look-ahead × samples × one decoded block` + `samples × per-open-file state`.
Both multiplicands are priced in §7. At look-ahead 1 this collapses to `samples × one decoded
block` paid once — which is the low-memory shape the previous draft deferred as a separate
"lockstep fallback" driver. It is no longer a separate anything: it is this object at its
smallest setting.

---

## 6. What a psp header records — **changed from the previous draft**

### 6.1 The fields

- **The analysed regions.** The one field compared across the cohort (§6.2). Recording the set
  buys two more things: a reader can tell an *unanalysed* position from an analysed position that
  had *no coverage*, which are different facts (the same distinction the census draws per
  position — [`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §1.1);
  and the census can be rebuilt from the psp, because census positions are chosen from the
  analysed stretch, so a rebuild that did not know that stretch would pick different positions.
- **The other values the segmentation was computed from** — the repeat catalog's identity
  (which carries the reference's identity) and the routing criteria. Compared against the calling
  run's own (§6.2), and recorded so a refusal can name the field that differs instead of
  reporting only that two things disagree — the same shape as the parameters fit's own
  compatibility check
  ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §5).
- **The read filters the walk applied, and the command line that produced the file.** Not compared
  by anything in this document — they are here because **the psp header already has a consumer in
  the code, and that consumer fixes part of its contents.** A census file names the psp it was
  built from by digesting exactly these things: *"a digest of the pileup's header — its reference,
  its analysed regions, its read filters and the command line that produced it"*, together with the
  psp's record count
  ([`census_file.rs:76-98`](../../../../src/ng/parameter_estimation/joint/census_file.rs)). A psp
  whose header omits them cannot be named, and a census built from it cannot be told from a census
  built from different reads — which is the whole point of naming it
  ([`Freshness`, `:106-127`](../../../../src/ng/parameter_estimation/joint/census_file.rs) turns the
  comparison into *use it*, *rebuild it*, or *refuse*).
- **How many records the file holds.** The second half of that identity, for the same reason.
- **When the file was written.** Never compared — it is what makes §12's byte-identity oracle
  read *"identical apart from the timestamp"* rather than *"identical"*, so it is a field the
  comparison deliberately skips rather than an omission.

**The census is the header's first consumer, and it is already built.** `PileupIdentity::of_header`
takes the psp header's bytes and a record count and says only this: two psps with the same header
and the same record count get the same identity and no others do — *"which bytes exactly is the
pileup writer's business"*
([`census_file.rs:91-98`](../../../../src/ng/parameter_estimation/joint/census_file.rs)). So the
encoding spec (§10) is free to lay the header out as it likes, and is not free to leave any of the
four values out.

### 6.2 The two refusals

- **Across the cohort: the analysed regions must be equal.** Two samples analysed over different
  segments are not comparable outside the ground they share — a sample has no records over ground
  it never looked at, and a caller that read that absence as homozygous reference would invent
  genotypes. This is the one cohort check, and the one mismatch that would produce a wrong answer.
  The refusal names both samples. (Whether the shared ground could be called instead of refusing
  is open question 5, §11.)
- **Each file against the run: the segmentation inputs must match.** The caller's segments are the
  segments it computes from the reference, catalog and routing criteria it was handed (§5.3). The
  observations in a psp were minted inside the segments of the *writing* run's segmentation, and
  "no observation crosses a segment's edge" (§4.3) is only true when the two segmentations are the
  same — under a different catalog or different routing floors, an STR observation could straddle
  a calling segment, and the independence argument of §4.3 fails. So a psp whose recorded catalog
  identity or routing criteria differ from the run's is refused, and the refusal names the first
  field that differs.

### 6.3 Dropped: the block-boundary digest, and its `writer_version` field

> **Flagged for the owner — this reverses the previous draft, which the owner has not ruled on.**
> The 2026-08-14 draft's §3.4 made every psp header carry a digest of the file's block boundaries,
> checked for equality across the cohort at open, plus a `writer_version` field so that the one
> mismatch the digest could not otherwise explain — equal inputs, unequal boundaries, therefore
> different writer code — could be named in the refusal.

**Decision: both fields are dropped.** The digest existed for one consumer: a calling stage that
merged the cohort in lockstep, block *k* of every sample covering the same stretch, so that no
sample's data ever waited for another's. In this design that consumer does not exist. No code
path ever looks at two samples' blocks together: each sample's psp reader independently serves
segments, decoding whichever of its own blocks overlap (§3.3), and the cohort merge is a streaming
merge keyed on coordinates whose resident set is its frontier (§3.2). Whether two files' blocks
line up changes nothing the run can observe — each reader holds one decoded block per cursor
either way. A check that guards nothing is a refusal waiting to fire on a harmless difference —
two writer versions in one cohort, say — and the previous draft itself classified the mismatch as
a performance matter, refused only "as a simplicity choice". The simplicity now runs the other
way: no digest, no comparison, no refusal to explain.

Two knock-on consequences, each following from the same removal:

- **`writer_version` goes with it.** Its stated purpose was to let a boundary-mismatch refusal
  name its cause; with no boundary check there is no such refusal. Versioning the psp *format* —
  so a reader can refuse bytes it does not understand — is a different job the encoding spec
  (§10) owns as a matter of course; nothing at the run level requires a writer-version field.
- **"Nothing about a sample's own data may move a block boundary" lapses.** The previous draft's
  rule, and its argument was cross-sample boundary equality alone: a data-dependent flush would
  cut a deep sample's blocks where its cohort's are not cut, and the digest would then refuse
  exactly the sample that most needed the run. With no equality to preserve, the rule has no
  remaining justification, and block sizing becomes wholly the encoding spec's trade — including
  a cap on a block's decoded size, which production's 1 MiB mid-window force-flush
  ([`src/psp/writer.rs:72,289-296`](../../../../src/psp/writer.rs)) shows the shape of and which
  the previous draft had to forbid. The pathological-block question (the densest stretch of a
  300× sample decoding to an unbounded payload) thereby gets its natural fix back. One property
  survives the hand-over because an oracle depends on it: block cuts must be a function of the
  observation stream alone — a stream that is identical at every worker count — so that the psp
  stays byte-identical across worker counts (§12.1). A cut that depended on scheduling would
  break that; a cut that depends on the observations does not.

If the digest is ever wanted again, the argument for it will have to name a consumer of boundary
equality, and this design has none.

---

## 7. Where the memory goes

### 7.1 The three bounds

| object | peak resident |
|---|---|
| `SampleObservationGatherer` | `look-ahead × one segment's observations + working set` + census accumulator (~6 MB per read group) |
| `PspVariantCaller` | `look-ahead × samples × one decoded block` + `samples × per-open-file state` |
| `AlignedFilesVariantCaller` | `look-ahead × samples × one segment's walking set` + `samples × 11–15 MiB` open alignment files |

Only the first is independent of the cohort size, which is the reason psp mode exists.

Pricing the `PspVariantCaller` row at the far end of the committed range — three thousand
samples, look-ahead 8, and per sample a 16 kB read buffer plus roughly 10 kB of decoded
observations per block (both estimates, not measurements): `8 × 3,000 × 26 kB ≈ 620 MB`.
Workable — but only because a sample costs 26 kB. At a megabyte per sample the same run needs
24 GB. That arithmetic is what turns "keep per-sample state small" from a preference into the
requirement below, and re-running it on real files once the encoding exists is part of open
question 2 (§11).

### 7.2 Requirement: an open psp costs tens of kilobytes, not megabytes

Everything in the calling stage multiplies by the sample count, and the per-open-file state is
the easiest cost to get wrong because it looks like bookkeeping. Production's psp index is the
counter-example ng would otherwise copy: a flat vector of one 24-byte entry per block
(`BlockIndexEntry`, [`src/psp/index.rs:42`](../../../../src/psp/index.rs)), decoded whole at open
(`decode_index`, [`:110`](../../../../src/psp/index.rs)). At a 5 kb block over an 800 Mb genome
that is 160,000 entries — **3.8 MB per open file, 11.5 GB across three thousand samples** —
before any data is read.

The shape of the fix: index at a much coarser grain and chain blocks within it — each block
carrying enough to reach the next — so a reader seeks once and then streams. A few hundred to a
few thousand index entries per file is kilobytes. This constrains the psp encoding without
specifying it; the encoding spec (§10) inherits the requirement, together with the reader
contract of §3.3.

---

## 8. Traps — what will bite the coder

Each is a property of code that exists today, and each produces a wrong answer or a silent
under-count rather than a crash.

- **A locus generator holds state across segments, so it cannot be shared between workers.** The
  iterator that owns the generators documents a load-bearing drop order
  ([`src/ng/locus_generation/mod.rs:706-737`](../../../../src/ng/locus_generation/mod.rs)). Each
  in-flight segment's task builds its own generator set — or a worker reuses one set across *its
  own* consecutive segments, never across workers.
- **The reference accessor is `Send` but deliberately not `Sync`.** `WindowedRefSeq` holds an
  open per-contig reader; the input layer takes a *factory* rather than a shared accessor for
  exactly this reason
  ([`src/ng/read/input/mod.rs:606-611`](../../../../src/ng/read/input/mod.rs)). One accessor per
  worker per file; wrapping a shared one in a lock serialises the walk's hottest path.
- **Whether `SampleReads` is `Sync` is unverified.** Its cursor is owned and `Send` (test
  `a_sample_cursor_is_send_in_both_arms`,
  [`src/ng/read/input/mod.rs:1441`](../../../../src/ng/read/input/mod.rs)), so the intended
  shape — one shared `SampleReads` per sample, one cursor per worker via `cursor(&self, ...)`
  ([`:623`](../../../../src/ng/read/input/mod.rs)) — needs `SampleReads: Sync`. Confirm at
  implementation time before building the pool on it.
- **Read-filter tallies live in the cursor, not the file.** They belong to a cursor from the
  moment it is made ([`src/ng/read/input/mod.rs:620-622`](../../../../src/ng/read/input/mod.rs)),
  so per-worker cursors give per-worker tallies. Sum them when the gatherer or caller finishes,
  or drop rates under-report by a factor of the worker count — silently, since every number stays
  plausible.
- **Dropping the walk's pieces in the wrong order loses the tallies too.**
  `SampleLocusObservationsIterator` releases its generators before its reads, and the comment
  records that no test can fail if that breaks
  ([`src/ng/locus_generation/mod.rs:707-737`](../../../../src/ng/locus_generation/mod.rs)). A
  task that dismantles the pieces itself must keep the order.
- **A backward jump is legal and costs a seek plus a block decode.** The read cursor answers
  segments in any order, and a test asserts backwards-walked segments return what a linear scan
  returns ([`src/ng/read/input/cursor.rs:92-96`](../../../../src/ng/read/input/cursor.rs), test
  at [`:1207`](../../../../src/ng/read/input/cursor.rs)); the psp reader's coarse index (§7.2)
  has the same character. So a work-stealing pool is correct however it steals, and slow unless
  each worker's own sequence of segments stays monotonic — schedule for that.
- **Chain ids are already omitted for reads that agree with the reference**, in both generic
  paths
  ([`src/ng/locus_generation/pileup/fast_column.rs:312`](../../../../src/ng/locus_generation/pileup/fast_column.rs),
  [`src/ng/locus_generation/pileup/open_record.rs:472`](../../../../src/ng/locus_generation/pileup/open_record.rs)).
  Recorded so nobody re-adds them when the psp writer is built: production's reference-side
  equivalents were about 31% of its peak live heap.
- **An analysed stretch with no reads must read back as analysed-and-empty, not as unknown.**
  Inside a run the gatherer closes this itself (§5.2). Across the files it survives through the
  header's analysed regions plus the trailer: a psp with a valid trailer covers everything its
  header says was analysed, so absence inside that ground means *no observations*, and a file
  without a valid trailer must be refused as interrupted rather than read as a short sample (§9).

---

## 9. Cross-cutting concerns

**Errors.** The iterators yield `Result` items: a worker's failure surfaces as an error naming
the sample and the genome span, and ends the iteration. A gatherer that dies part-way leaves a
psp without a valid trailer, and the reader must refuse such a file rather than read it as a
short one — an interrupted sample must not pass for a sample with fewer observations. In direct
mode there is nothing to clean up: the VCF is truncated at the last yielded segment and the run
reports where it stopped.

**Concurrency.** Shared and read-only: the reference bases, the repeat catalog, the segments,
the read-group table, the model parameters, and each sample's open file (`SampleReads` or psp —
the former's `Sync` to be confirmed, §8). Mutable and per worker: cursors, reference accessors,
generator sets, decoded blocks. Single-threaded by construction: the yield point — each object
is an iterator, so its consumer (the VCF writer, the psp writer, the census accumulator) runs on
the consuming thread and needs no lock.

**Performance.** The knobs are the look-ahead and the worker count, per object, plus
samples-in-flight in the walk stage (§5.2). Block size is the encoding spec's trade (§6.3) —
it reaches this design only as "one decoded block", the multiplicand in §7.1. Any other tuning
constant that appears in the implementation is a defect, not a lever.

---

## 10. Deferred, with a recommended home

- **The psp file's encoding** — byte layout, compression, block sizing and any cap on a block's
  decoded size (a sample-data-dependent cut is now legal — §6.3), checksums, format versioning,
  the index, the trailer. Its own spec beside this one; it inherits the header values (§6.1), the
  reader contract (§3.3 — segments served, blocks internal, observations in coordinate order),
  the per-open-file budget (§7.2), and one restriction: block cuts are a function of the
  observation stream alone, so the file is worker-count invariant (§12.1). One inherited sizing
  rule: whatever target the writer cuts toward must count the segments' own bases, not the
  chromosome span — under a BED the gaps belong to no segment, and a span-counted target would
  make blocks over sparse ground nearly empty payloads that cost their framing and carry little.
- **The cohort merge's reconciliation** — how per-sample observations whose spans differ become
  one cohort observation, and the confirmation that the merge frontier stays bounded while doing
  it (§3.2). Its own spec; production's cohort variant grouping and genotype joining are the
  reuse candidates; `MergedRegionReads` ([`sample_reads.md`](../arch/sample_reads.md) §4) is the
  streaming shape to keep.
- **The parameters file's format** — what the user supplies in direct mode and the fit writes in
  psp mode. Its own spec, **on direct mode's critical path** (§2): the mode cannot run without
  it.
- **The VCF writer** — consumes a caller's iterator; the variants arrive in genome order, so it
  writes as it reads. Its shape, and the `Variant` record's, belong to the emission step's
  document.

---

## 11. Open questions

1. **How long is a segment at the routing floor?** — OPEN. Measured only at the catalog floor:
   391 bases average (§4.4). *Leaning:* kilobases at the routing floor. **Settled by:** counting
   segments and their length distribution over the existing catalog file at both floors, tomato
   and human — a filter over a stored file, not a genome scan.
2. **The look-ahead's default, per object.** — OPEN; no value proposed. The gatherer's and the
   callers' cost per unit of look-ahead differ by a factor of the sample count, so one default
   will not serve all three; the `PspVariantCaller`'s should come from the cohort size, not the
   core count. **Settled by:** sweeping look-ahead on one tomato sample (gatherer) and on the
   tomato cohort, in both modes — wall time and peak resident — and, for the psp caller,
   re-running §7.1's arithmetic with the two per-sample costs measured on real files.
3. **Does a within-sample segment loop scale, and what does per-segment fetching cost when
   segments are short?** — OPEN, and the first half can invalidate §5.2's schedule (production:
   1.81× at four threads, worse at eight). **Settled by:** driving the gatherer at 1, 2, 4, 8, 16
   workers on a tomato sample and HG002 — wall time, peak resident, observations identical to
   serial. The second half needs its own run: how often a segment's fetch lands in the
   already-decoded block, on a whole-genome walk and on a fragmented BED over a full, un-sliced
   CRAM. The nearest measurement — production's `--regions` costs a flat ~14% (24.7 s → 28.2 s)
   and does not grow from 80 to 4,000 intervals (29.0 / 28.8 / 28.5 s) — does not cover that
   shape.
4. **Should an observation's reference bases be stored in the psp?** — OPEN.
   `SampleLocusObservations::reference_bases` is a `Box<[u8]>` per observation
   ([`src/ng/locus_generation/mod.rs:44-46`](../../../../src/ng/locus_generation/mod.rs));
   written out, it is a per-sample copy of the reference — megabytes per sample at 7.4 million
   generic observations, multiplied by the cohort. *Leaning:* do not store it; re-fetch when a
   block is decoded, since the calling stage holds the reference anyway. **Settled by:** measuring
   the file-size difference on one sample and the re-fetch's cost in calling wall time.
5. **Refuse a cohort analysed over different segments, or call the shared ground?** — OPEN.
   Refusal is what §6.2 specifies and the safe first behaviour, but the case is legitimate — a
   panel run and a whole-genome run in one cohort, or calling one chromosome from whole-genome
   psps. *Leaning:* call the intersection and emit no-call — never homozygous reference — for a
   sample over ground it did not analyse, the rule the census already applies per position; that
   needs per-sample analysed regions carried into the emission step, so it is not free.
   **Settled by:** deciding when the emission step is specified; until then, refuse and name both
   files.
6. **Where does direct mode stop being usable?** — OPEN. The open-file bill is measured (§5.1);
   the in-flight walking sets are not, and they depend on depth. *Leaning:* low hundreds of
   samples at ordinary depth. **Settled by:** running direct mode at 1, 10, 50 and 200 samples on
   the tomato cohort, peak resident against sample count.

---

## 12. How we know it works

Each oracle is a property of the run, not of one type — which is why they live here.

1. **Worker-count invariance of the psp.** One sample gathered at 1, 2, 4, 8, 16 workers gives
   byte-identical psps apart from the header's timestamp. Production already holds this for its
   own pileup output across thread counts, so it is a reachable bar — and it is what §6.3's
   restriction on block cuts preserves.
2. **Worker-count and look-ahead invariance of the VCF**, from each of the two callers.
3. **Mode equivalence — the oracle that justifies the design.** The same cohort and the same
   parameters, run through `AlignedFilesVariantCaller` and through the psp route, give the same
   VCF. This is simultaneously the proof that the calling function is mode-blind (goal 1) and the
   sufficiency test for the psp: anything the file fails to carry surfaces here, where a
   write-read round-trip test would pass.
4. **Segment independence of the observations.** One segment walked alone emits exactly the
   observations the same span emits inside a whole-genome, single-threaded walk. This asserts
   §4.3 is honoured; the thirds-chopping test that lost 17 positions is the failure shape it
   catches.
5. **The two refusals stay two.** Unequal analysed regions refuse naming both samples; a psp
   whose segmentation inputs differ from the run's refuses naming the first differing field.
   Both fire at construction, before any block is decoded (§6.2).
6. **The census built during the walk equals the census built from the psp.** Specified already
   as a byte-for-byte comparison
   ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §7.12) and named
   here because this document is what puts the two producers in one run: it is the sharpest test
   that a psp holds everything a census needs.
7. **Analysed-but-empty survives a round trip.** A stretch a sample analysed and found no reads
   in reads back as analysed and empty — distinguishable, via the header's analysed regions, from
   a stretch the sample never looked at (§8, last trap).

---

## 13. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the segments the loops run over | `TypedRegion`, `RegionKind` ([`src/ng/region_typing/mod.rs:144,168`](../../../../src/ng/region_typing/mod.rs)) | consumed as-is; a segment handed to `call_vars_in_segment` is one of them |
| the analysed regions | `GenomeRegions` ([`src/ng/region_typing/mod.rs:77,87,100`](../../../../src/ng/region_typing/mod.rs)) | reused whole — it wraps production's `RegionSet`, so ng and production agree on what a BED means; its value is recorded in the psp header |
| what a BED edge does to a segment | `clips_at_a_bed_edge` and the emission rule ([`src/ng/region_typing/mod.rs:471,482-488`](../../../../src/ng/region_typing/mod.rs)) | taken as given — findings whole, generic clipped; nothing here re-decides it |
| one sample's observations | `SampleLocusObservations` ([`src/ng/locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | the item of every stream in §3, unchanged |
| the walker behind a source | `SampleLocusObservationsIterator` ([`src/ng/locus_generation/mod.rs:706`](../../../../src/ng/locus_generation/mod.rs)) | one per task, fed that segment's segments |
| per-segment reads | `SampleReads` and `cursor()` ([`src/ng/read/input/mod.rs:398,623`](../../../../src/ng/read/input/mod.rs)) | one shared `SampleReads` per sample, one owned cursor per worker (`Send` proven; `Sync` to confirm — §8) |
| the streaming merge's shape | `MergedRegionReads` ([`sample_reads.md`](../arch/sample_reads.md) §4) | the model for `k_way_merge`: argmin over per-stream heads, keys beside the heads, frontier-sized residency |
| the census accumulator | `CensusWriter::add_locus`, `mark_walked`, `finish` ([`src/ng/parameter_estimation/joint/census.rs:1965,1984,2252`](../../../../src/ng/parameter_estimation/joint/census.rs)) | fed inside the gatherer, at the ordered yield point (§5.2) |
| the census file | `write_census`, `open_census` ([`src/ng/parameter_estimation/joint/census_file.rs:195,421`](../../../../src/ng/parameter_estimation/joint/census_file.rs)) | written by the walk stage's per-sample loop from `finish()`'s result |
| psp block index | `BlockIndexEntry`, `decode_index` ([`src/psp/index.rs:42,110`](../../../../src/psp/index.rs)) | **a model of what not to build** — §7.2 rejects the flat per-block index at ng's sample counts |
| block cutting | `PspWriter`'s grid and force-flush ([`src/psp/writer.rs:297-301,72,289-296`](../../../../src/psp/writer.rs)) | neither carries as a rule: block sizing is wholly the encoding spec's, and the force-flush is now a legal shape for capping a block's decoded size (§6.3) |

**The parity oracle for the whole document is §12.3** — direct mode against psp mode, one cohort,
parameters held fixed.
