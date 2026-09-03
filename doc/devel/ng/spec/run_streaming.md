# ng — how a run is driven: three iterators, two modes

*Status: design spec, 2026-08-16, **amended 2026-08-28**. No code yet — this settles the design.
The 2026-08-28 amendment changes where the parallelism is and nothing else: the walk goes several
samples at a time with each sample serial inside (§5.2), and the two callers stop cutting the
genome into segments for workers and instead call the merge's loci from a pool (§3.1, §3.5). §7.1's
memory bounds and §11's questions follow from that; **every concurrency default is now unset on
purpose**, because the psp format's cost, the calling loop's cost and both of their memory
footprints are unmeasured. The public shape is three objects, each an iterator, and every division
of work is internal to one of them. **One thing waits on the owner: §6.3 drops the psp header's
block-boundary digest**, and its `writer_version` field with it. **⚠ The companion architecture doc has not been amended and now specifies the deleted design**:
[`../arch/run_streaming.md`](../arch/run_streaming.md) still gives `call_vars_in_segment`, a
`LookAhead` knob, a segment pool, and per-segment per-sample walkers, all of which §3.1 and §3.5
retire. Read it for the types it names and not for the shape; amending it is owed. Reads on
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
   to several hundred (`CLAUDE.md`, *What this caller has to work on*). Every memory bound in §7 is
   a formula whose terms are stated, so the large end has an answer rather than an omission. **The
   walk stage's bound is in *samples in flight* rather than in cohort size** (§7.1), which is what
   makes it the one stage whose peak does not grow with the cohort.
3. **A single sample must be able to use the machine.** Production's single-sample pileup reaches
   1.81× at four threads and gets *worse* at eight (59.2 s → 32.7 s → 34.4 s), because only its
   BAQ stage threads while reading, walking and writing stay serial — about 40% of that run cannot
   be sped up
   ([`pileup_thread_scaling_2026-06-11.md`](../../reports/pileup_thread_scaling_2026-06-11.md)).
   **This goal is not met today.** The walk runs several samples at once with each sample serial
   inside (§5.2), which is right for a cohort and gives a lone sample nothing. **Question 8 (§11)
   owns what to do about it** and is the only place this document argues it.
4. **The output must not depend on how the work was divided.** The same VCF at one thread and at
   sixteen; the same VCF from direct mode and from psp mode with the parameters held fixed. This
   is the oracle for everything below (§12).

### 1.2 Non-goals, and what this document does not do

- **It does not define the psp file's encoding** — byte layout, compression, block sizing,
  checksums, format versioning, the index. It fixes only the contract the file's reader must
  answer (§3.3), what the header must record (§6), and the resident-state budget an open file
  must meet (§7.2). The encoding is deferred with a home (§10).
- **It does not define the cohort merge itself** — how the samples' streams become one stream of
  cohort loci, including how observations whose spans differ between samples become one
  cohort locus. That is [`cohort_merge.md`](cohort_merge.md), written the day after this
  document and built since. What this document fixes about the merge is only what the run needs
  from it: that it keys on coordinates and that it streams (§3.2).
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
- **cohort locus** — every sample's observations over one stretch of genome, brought together by
  the merge; the unit the caller consumes. How overlapping spans are reconciled inside one is
  [`cohort_merge.md`](cohort_merge.md)'s. **This is the only name this document uses for it**, and
  it is the merge's spec's and the code's.
- **census** — a sample's evidence at a scattered subset of the analysed positions, chosen in
  advance and identical in every sample, which is all the parameters fit reads. About one
  analysed position in 390 on tomato.
- **psp** — the per-sample file holding all of one sample's observations, in genome order.
- **block** — how far back the psp's compressor is allowed to look for a repeat. **A reader can
  decode less than a whole one before handing out a record**, and the
  encoding work since 2026-08-19 is trying to make it much larger than that. Decompressing a block
  incrementally separates *how far back the compressor may look* from *how much memory a reader
  spends*, so a block may be a megabyte while a reader holds far less
  ([`../impl_plan/psp_encoding_experiments.md`](../impl_plan/psp_encoding_experiments.md), the
  first of the three things it measures). Production couples the two, and that coupling is what
  ng is trying not to inherit ([`psp_record_encoding.md`](psp_record_encoding.md) §1.1).
- **a reader's working set** — what one open psp costs while it is being read. **This, not the
  block size, is the quantity every memory bound in this document multiplies by the sample
  count.** How the encoding keeps it small is the encoding's business; that it must fit §7.2's
  budget is the only thing this document asks of any shape it chooses.
- **source** — one sample's whole run of observations, in coordinate order, as an iterator that
  only moves forward. It is *asked for* one segment at a time (§3.4), but it is not rebuilt per
  segment and not held per worker: a caller has exactly one per sample for the whole run (§3.5).
  Two exist — a walker over the sample's alignment files, and a reader over its psp (§3.3).
- **samples in flight** — how many samples the walk stage walks at once (§5.2). The walk's memory
  knob.
- **callers in flight** — how many cohort loci the two callers may have out with workers at once
  (§3.5). Small, and the only concurrency setting on the calling side.

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

### 3.1 One locus at a time, through one function

Both variant callers — the one reading alignment files and the one reading psp files — drive one
merge and one calling function:

```
# a "source" answers one question: this sample's observations, as an iterator,
# in coordinate order, advancing forward only, for the whole run.
#   the walker      – reads the sample's alignment files
#   the psp reader  – decodes the block it is in, and the next when it passes the end

next():
    while pending is empty:
        cohort_observation = merge.next()   # merge exhausted → iteration ends
        pending = collect(call_vars_from_observation(cohort_observation, parameters))
    return pending.pop_front()
```

where `merge` is the k-way merge over the sources (§3.2), yielding one cohort locus at a
time in genome order.

**The two callers differ only in what a source is.** `k_way_merge` and
`call_vars_from_observation` are written once and neither can tell a walker from a file reader,
which is the whole of goal 1.

**The skeleton above is the contract** — what the iterator yields, and in what order. §3.5 says
how the same yield sequence is produced with a pool of callers, and the sequence is required to be
identical either way (§12).

### 3.2 The merge is streaming

`k_way_merge` consumes the per-sample iterators and yields one cohort locus at a time. Its
resident set is its **frontier** — the head observation of each sample's iterator, so about one
observation per sample — never a whole segment's observations per sample. Nothing accumulates a
segment before merging it; the merge pulls, the sources produce lazily, and a cohort locus is
dropped the moment the caller has consumed it.

The repo already has this shape in its read layer: `MergedCursors`
([`src/ng/read/input/sample_cursor.rs:178`](../../../../src/ng/read/input/sample_cursor.rs)) is an
argmin k-way merge over a sample's per-file cursors, holding one head per stream with the sort keys
beside the heads, and it is the model for `k_way_merge`. *(The architecture doc calls this shape
`MergedRegionReads` ([`sample_reads.md`](../arch/sample_reads.md) §4); that name is from a
superseded per-region design and is not in the tree.)* The differences are the item (observations, not reads) and the yield (a cohort locus
groups every sample's entries at one stretch of genome rather than interleaving
them); where two samples' generic observations overlap without coinciding, the frontier may
briefly hold a few observations of one sample while spans are reconciled —
[`cohort_merge.md`](cohort_merge.md) §4.2 owns that reconciliation and its §4.5 confirms the
frontier stays bounded through it.

**⚑ "About one observation per sample" describes the merge this document imagined, not the one that
was built.** The built merge reads through an observation cache spanning the ground its builders
cover, and `cohort_merge.md` §8 measures **33 records held per sample** at 1,000 samples on the
shipped defaults. §7.1 carries the consequence.

**The consequence: segment size is not a memory knob.** A segment's cost while being called is its
merge frontier plus the sources' own working state, not its length. Segments can therefore be
exactly the segments (§4.4), with no machinery to group them toward a size target — the previous
draft's walk units and calling ranges existed to hit size targets that no longer control anything.

**The merge keys on coordinates — decided, carried forward.** Generic observations are split out
of each sample's own data, so they do not line up between samples. A merge keyed on anything
positional-within-a-file — a block number, an index within a decoded payload — would be wrong;
coordinates are the only key all producers share.

### 3.3 Most positions produce no cohort locus at all, and the saving lands in two different places

**A psp holds every position, because no sample knows which positions another sample varies at.**
What is skipped is not the storing and not the visiting — it is the *building*. A position where
every sample's reads matched the reference produces no cohort locus and no call. Measured on 10
tomato bench samples over 300 BED regions: **3,060 of 309,018 covered positions varied — one
position in 101**
([`../research/experiments/locus_stream_shape/sketch4_columnar_producer_plus_fold.md`](../research/experiments/locus_stream_shape/sketch4_columnar_producer_plus_fold.md)
§2).

**The decision is a cohort one and cannot be taken per sample.** A position where *this* sample saw
only reference reads may be a variant site in another, and there this sample's evidence is still
needed — its depth is what separates a confident homozygous reference from no coverage, and it is
part of the allele frequency every genotype at that locus is weighted by. So no sample's quiet
position is dropped on its own account; what is dropped is a position no sample varied at.

**Where the saving is taken differs by mode, because the expensive thing differs.**

- **Direct mode: at the merge.** The walker has already read the reads, so its observations exist
  whatever happens next. The saving is downstream — where the cohort's non-reference evidence at a
  position is nothing, `k_way_merge` builds no cohort locus and `call_vars_from_observation`
  is never entered.
- **psp mode: while reading the file.** The evidence a sample stored is not all equally expensive
  to get at. Deciding *whether a position is worth calling* needs two small numbers from each
  sample — how many of its reads were non-reference there, and how many of its reads were
  **compared against the reference** over the whole locus. That second one is neither read depth
  nor the count of reads that covered the ground: a read whose bases stop inside the locus is in
  neither number ([`cohort_merge.md`](cohort_merge.md) §4.3). A file that cheaply stored *depth*
  would not answer the rule.
  Reconstructing *what that sample actually saw* needs everything else: the sequences, the
  qualities, the read names. **So the file is written so the small numbers can be read without
  reconstructing the rest**, and a run reads them for every position but reconstructs the rest only
  where the cohort decided to call.

  **This is production's one good idea here, and its layout is not being copied.** Production does
  the same trick over a transposed block — the cheap fields decoded for every row, the expensive
  ones left compressed, then inflated for the rows that were kept (`TwoPhaseSegment` and
  `set_variable_rows`,
  [`var_calling/sample_reader.rs:698-712,789`](../../../../src/var_calling/sample_reader.rs)).
  **How ng's file is cut up so the same saving is available is not settled** — the encoding work is
  sweeping it, because the cheaper the cut is in bytes the more a reader must hold to use it, and
  the two pull against each other ([`../impl_plan/psp_encoding_experiments.md`](../impl_plan/psp_encoding_experiments.md),
  its milestone B0). **What this document requires is only that the saving survives whatever wins**,
  and that the two small numbers above are among what a run can read cheaply. That second half is a
  requirement the encoding did not have written down: the rule that admits a position asks each
  sample for a *share* of its own reads as well as a flat count
  ([`cohort_merge.md`](cohort_merge.md) §4.3), so the read count is not optional.

**And there is a cheaper step in front of that one.** For each stretch of genome
the file stores together, it also records **the most non-reference reads any position in that
stretch carried, readable without decompressing anything**. A reader can compare that against the
bar a position must clear to be worth calling and, when no position in the stretch could have
cleared it, move on without reading the stretch at all.

**The test is the admission bar and not "any non-reference read at all" — the difference decides
whether this is worth having.** Sequencing error alone puts a non-reference read somewhere in any
substantial stretch: at 3 reads a position over 5,000 positions with Q30 bases, about fifteen
erroneous bases are expected, so a test of *any* would answer "no" almost never and the step would
save nothing. The bar — two non-reference reads, or two in a hundred of that sample's reads,
whichever is larger ([`cohort_merge.md`](cohort_merge.md) §4.3) — is cleared by error alone far more
rarely, and **the share half is what keeps that true at depth**: at 300 reads a position two
non-reference reads is noise, six in three hundred is not. **A stored summary that cannot be
compared against the share half makes this step inert at exactly the depths where a heavy record is
largest** — the same requirement the ⚑ below states for step 2, and one number serves both.

**Passing over a sample is not the same as leaving it out, and confusing the two is the trap.** Two
different questions get asked at a genomic position, and they need different things:

- **Is this position worth calling at all?** One sample showing enough non-reference reads admits
  it for everybody ([`cohort_merge.md`](cohort_merge.md) §4.3). A sample that showed no
  non-reference read anywhere in the stretch cannot be the one that admits it, so passing over that
  sample cannot change this answer.
- **What is each sample's genotype at that position?** This needs *every* sample's evidence, the
  quiet ones included. Thirty reads that all match the reference is a confident homozygous
  reference; no reads at all is a no-call. A sample left out entirely makes those two
  indistinguishable.

**So the shortcut may excuse a sample from the first question and never from the second.** A reader
that passed over a stretch must still be asked for it, in full, at every position the cohort ended
up keeping — and must answer.

**Passing over is one sample's decision, never the cohort's.** Each sample's file is cut into
stretches by its own data, and those cuts do not coincide between samples (§6.3). One sample can
skip ground another is reading in full. There is no cohort-wide skip and none is needed: the saving
is already taken sample by sample.

**So reading one sample's evidence happens in three steps.** The first exists only when the
evidence comes from a psp; a walker over alignment files starts at the second, because it is
holding the reads already.

```
# step 1 — psp only, and it reads nothing but a stored summary
did this sample show a non-reference read anywhere in this stretch?
#   "no" excuses the sample from step 2 over that stretch. It is still asked in
#   step 3 at every position the cohort kept, and must answer there in full.

# step 2 — cheap, at every position of the segment
how many non-reference reads, and how many reads in total, did this sample have here?
#   the merge puts every sample's answers together and decides which positions are
#   worth calling

# step 3 — expensive, only at the positions the merge kept
everything the sample saw there — sequences, qualities, read names
```

**Step 2 costs the two sources completely different things, and the merge cannot tell.** A walker
has the reads in hand, so step 2 costs it a walk and step 3 is free; a psp reader pays a cheap read
for step 2 and a dear one for step 3. **One interface, two cost profiles** — and the merge, which
sees only the answers, is where the decision is taken either way.

**⚑ The read count in step 2 is a requirement the encoding did not have written down.** The rule
that admits a position asks each sample for a *share* of its own reads as well as a flat count —
`max(floor reads, share × that sample's reads compared against the reference over the locus)`
([`src/ng/run/cohort_merge/mod.rs`](../../../../src/ng/run/cohort_merge/mod.rs)) — so a cheap read
that returns only the non-reference count cannot answer it. Named here because this document is
where the requirement comes from; the encoding is where it has to be met.

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
  segment, it works out internally where in the file that ground lies, decodes what it needs, and
  yields the observations inside the segment. It keeps its place in the file, so a segment
  following the one it just served costs no fresh seek and no re-decode. The block never surfaces
  above the reader — no caller mentions one, and there is no "line the segment up with a block"
  problem: any segment is servable, and consecutive segments are simply the case the reader is
  warmest for. **How much it holds while doing this is §7.2's budget, not a block's size** (§1.3,
  *a reader's working set*).

Both sources are cursors: asked for ascending segments they stream forward; a backward jump is
legal and costs a seek (§8). A source serves one consumer at a time. **In a caller there is one
source per sample for the whole run** and it only ever moves forward, because the single merge is
its only consumer (§3.5) — so each stretch of the file is decoded once, by one cursor, and the
backward jump never happens. In a gatherer splitting one sample across workers, each worker needs its own source over
that sample's one open file (§5.2).

### 3.5 Where the parallelism is

**Everything in this section is provisional until §11's measurements exist** — the psp format's
cost, the calling loop's cost and the memory of both are unknown, and no arrangement here should be
improved on before they are known.

**The two stages parallelise along different axes, because their work has different shapes.**

**The walk: one worker per sample.** Samples are independent — nothing one sample's walk computes
is read by another's — so the pool is over samples and each sample's walk is serial inside. §5.2
gives the reasoning and what it costs.

**The two callers: a serial merge feeding a pool of loci.** The merge produces cohort loci one at
a time in genome order on a single thread; each locus is handed to a free worker as it appears;
the called results are released in genome order. Nothing is decided in advance and no line is
drawn through the genome.

**Why the genome is not cut for calling, and this is the part that is easy to get wrong.** Two
reasons. **The first carries the argument on its own**; the second is a guess until §11's question 1
measures it, and is here because it says what would go wrong if the first were somehow solved:

- **The loci are not known until the merge has made them.** A cohort locus is not a position: a
  deletion joins consecutive positions into one locus, so where one begins and ends is an output
  of the merge and not an input to it. Anything that hands out loci in advance has to have merged
  them first.
- **The one genome cut that respects loci — the segment — is wildly uneven.** An STR tract is tens
  of bases; a generic stretch between tracts is unmeasured and probably far longer (§4.4 gives the
  arithmetic and says why it is a guess). Give segments to a pool and one enormous generic segment both
  starves the pool, since only one worker is inside it, and forces every other worker's finished
  output to wait behind it — which is the memory the ordered release then has to hold.

**What is in flight, and it is small.** `callers in flight × one cohort locus`, plus whatever the
sources hold at the merge frontier (§3.2) and the merge's own observation cache (§7.1). A locus is the same size whether it came from a kilobase of
quiet generic ground or from a tract, so the bound does not depend on how the genome is shaped —
which is exactly what the segment pool could not promise.

**The merge has its own parallelism, and it is off unless measured to pay.** `merge_cohort_in_parallel`
([`src/ng/run/cohort_merge/parallel.rs`](../../../../src/ng/run/cohort_merge/parallel.rs)) cuts
each analysed region into fixed-width **building regions** — 200 bases by default — works a round
of them concurrently, and releases their loci by region index. It is byte-identical to the
single-threaded merge at any region width and any number in flight. **Its measured gain is poor,
and its own design says why:** the organiser draws the readers forward while no builder runs, then
a round of builders runs, and nothing is released until every builder in that round has finished —
so a barrier fires every `regions_in_flight × regions_len` bases — 3,200 on a sixteen-thread
machine at the 200-base default width, because regions-in-flight has **no constant default** and
takes one per worker thread
([`cohort_merge/mod.rs:592`](../../../../src/ng/run/cohort_merge/mod.rs)) — with a serial phase
between rounds. It stays available and off; §11 question 7 names what would turn it
on and what to try first.

**⚑ This arrangement puts the largest known cost on one thread, and that is the risk §11's
question 7 exists to check.** Everything each sample's source does — opening its data, decoding it,
minting observations — happens where the merge pulls, which is one thread; the pool gets only the
genotype arithmetic. The merge's own assembly is on that thread too, and it is not small:
`cohort_merge.md` §4.3 measures projecting and assembling a locus at **170 ms of a 425 ms
single-threaded merge** over 100 kb at 63 samples. **Nothing here is defended as optimal** — it is
the simplest arrangement that is obviously correct, chosen so that the measurement decides the rest
rather than an argument.

**Which of the two stages is even worth a pool is not yet known.** Production's cohort
variant-calling profile put its EM at about 3% against roughly 30% for the producer's decode
([`cohort_varcalling_perf_2026-07-03.md`](../../reports/cohort_varcalling_perf_2026-07-03.md)) — CPU
shares, and that report warns in the same passage that a 30% CPU share was "real but not
wall-relevant". So it is a reason to suspect that in the psp-to-VCF path the work is mostly
*decoding k psp files* rather than calling loci, and not a reason to believe it. ng's calling loop is heavier than production's — several passes, a
genotype table per locus — so that is a reason to measure the split before building either pool,
not a reason to assume it carries over (§11, question 7).

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

### 4.4 The segment is the unit of *observation generation* — decided

Locus generation advances one segment at a time: the gatherer's loop, and each source inside a
caller, produce a segment's observations before moving to the next. No grouping of segments into
larger hand-out units exists in this design.

**It is not the unit of parallel work** (§3.5). The argument below is about a segment being too
small a *step*, which holds whether the loop is serial or not.

**A measured average of 391 bases (613,682 generic segments on human chromosome 1) does not apply
here, and the owner ruled so.** That figure was taken at the *catalog's*
admission floor, `CATALOG_MIN_COPIES = [5, 5, 4, 4, 4, 3]`
([`src/ng/repeat_catalog/criteria.rs:16`](../../../../src/ng/repeat_catalog/criteria.rs)), which
admits a five-base homopolymer — in random sequence a run of five or more identical bases begins
about one position in 341, against a measured 391. **That is agreement in the wrong direction and
should not be read as confirmation**: the calculation counts homopolymers only, while the catalog
admits periods 2 to 6, all of which also cut generic ground — so it should over-predict cuts and
under-predict segment length, not the reverse. At the floor a caller actually routes on —
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

**How it parallelises:** the merge runs on one thread and its cohort loci go to a pool of callers
(§3.5). Every sample's walker advances at the merge frontier, in one place; the walkers are not a
pool here, because a walker only produces what the merge is about to consume and running them
ahead would buffer observations nobody has asked for. Each walker keeps its cursors and
generators for the whole run, which is what keeps every cursor's movement forward-only — the fast
path (§8), and the traps there are why cursors are never shared.

**What bounds it:** `callers in flight × one cohort locus` + the merge frontier (about one
observation per sample, §3.2) + `samples × 11–15 MiB` of open files. **The open files are the
whole bill**: everything else is small and does not grow with the genome's shape.

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
([`src/ng/parameter_estimation/joint/census.rs:1928,2378,1378`](../../../../src/ng/parameter_estimation/joint/census.rs)).

**The walk stage is a loop over samples**, one gatherer each. **Several samples are walked at
once, one worker each, and each sample's own walk is serial** (owner's decision). **And the
first implementation runs that loop at concurrency one** (owner's ruling, 2026-09-03): samples
are processed one at a time, in the order given, and a cohort is parallelised by running
invocations — typically one sample each — because each sample's generation is independent of
every other's. That independence is the critical difference from direct mode, which must hold
every sample open at one shared frontier; the walk stage never has to. The in-process fan-out
below remains what the arrangement permits, and its knob is §11 question 2's:

```
for each sample:
    gatherer = SampleObservationGatherer::new(sample's files, ...)
    psp writer consumes the gatherer            # blocks, header, trailer: the writer's business
    census file written from gatherer.finish()  # write_census — census_file.rs:200
```

**Why across samples rather than within one.** One worker per sample buys three things that
splitting a sample would cost back. Each psp is written from one deterministic stream, so the
writer never has to reassemble out-of-order work and §12.1's worker-count invariance comes for free
instead of being designed for. The census accumulator is fed from one thread. And the read-filter
tallies cannot under-report by the worker count (§8). Against those: within-sample scaling in ng is
unmeasured, and the only measurement anybody has is production's, which gets worse past four
threads (goal 3).

**Samples in flight is bounded and small, and it is not a thread per sample** — at a thousand
samples that would be a thousand open alignment files. Each sample in flight costs one open file
at 11–15 MiB plus its census accumulator at about 6 MB per read group. **At one read group a
sample** — every sample of both benchmark cohorts — six at once is under 150 MB. **The read-group
count is a real multiplier, not a formality**: the 300-reads-a-position end of the committed range
is reached by sequencing one sample over many lanes, and each lane is a read group, so the same six
samples at sixteen read groups each is closer to 700 MB. Whoever sets this knob must be told the
read-group count, not only the sample count. **Peak memory is a function of samples *in flight*, not of cohort size**, which is
the property psp mode exists for and which survives this change.

**A run of one sample gets one worker here, which is goal 3 unmet** — §11's question 8 owns it.
The one consequence to carry meanwhile: **the psp writer must keep honouring §12.1's worker-count
invariance**, because splitting one sample's walk is still on the table, and it is the only thing
that would ever make a sample's observations arrive out of order.

**The census is fed from inside the gatherer, on the same stream it yields.** Every observation
the gatherer yields passes the census accumulator (`CensusWriter::add_locus`,
[`census.rs:2087`](../../../../src/ng/parameter_estimation/joint/census.rs)) at the ordered yield
point, and every segment the loop completes is marked walked (`mark_walked`,
[`census.rs:2106`](../../../../src/ng/parameter_estimation/joint/census.rs)) whether or not it
produced observations. Two consequences:

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

**A gatherer never sees more than one sample's files**, whichever way the stage is parallelised.
The pool is outside it, over samples; the knob is inside it, over that sample's segments — its
segments to workers, each with its own cursor, reference accessor and generators, observations
drained in genome order. Read-filter tallies then live in each worker's cursor and are summed when
the gatherer finishes, or drop rates under-report by the worker count (§8).

**What bounds one gatherer:** its own segment concurrency × one segment's observations and working
set, plus the census accumulator. At the default that concurrency is one. Nothing crosses samples
at any setting.

### 5.3 `PspVariantCaller` — psp mode's calling

**Constructed from:** every sample's psp, the segmentation's ingredients, and the parameters.
Construction opens every file, reads every header, and runs the checks of §6.2 before any block
is decoded; the analysed regions come from the headers, not from a flag — the files know what
ground they cover. **Yields:** variants, in genome order.

**How it parallelises:** the merge on one thread, its cohort loci to a pool of callers (§3.5).
There is **one source per sample for the whole run**, not one per sample per in-flight unit: a psp
reader's cursor over that sample's one open file, advancing only forward, at the merge frontier.
So each stretch of the file is decoded once, by one cursor — where a source per in-flight segment
would have two workers sharing a block each decode it.

**What that costs is now measured** (§7.2): the open file's resident state is **357 kB on a human
reference and 7 kB on tomato**, almost all of it the reference's contig list, and the cursor
walking it adds **123 kB**. Both are paid once per sample, because under this design there is one
of each.

**What bounds it:** `samples × one reader's working set` + `callers in flight × one cohort locus`.
The first term is the run's floor and is paid whatever the concurrency; §7 prices it, and §7.2's
requirement is what keeps it affordable at three thousand samples. The second is small and is the
only part a concurrency setting moves.

**If this stage needs more parallelism than a pool of callers gives it, the next axis is the
decode** — every sample's cursor inflating its next block at the same time, still feeding one
merge — and not a second cut through the genome. Which axis is worth building is §11's question 7.

---

## 6. What a psp header records

### 6.1 The fields

- **The analysed regions.** Compared across the cohort (§6.2). Recording the set
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
- **The sample's name.** Without it a calling run has only the order its psp files were named on
  the command line, and the per-sample parameters — the inbreeding coefficient above all — are
  keyed by name ([`parameters_file.md`](parameters_file.md) §3.5). A reordered file list would then
  give every sample another sample's coefficient and put the wrong name on every VCF column, with
  nothing failing. The calling stage matches on this name and refuses a set it cannot match (§6.2).
- **The sample's read groups: each one's `@RG ID`, its library, and the identifier this walk gave
  it.** **Without this, no cohort can be assembled from separately-walked samples at all.** A
  gatherer sees one sample's files, so it numbers that sample's read groups from zero, and every
  sample's first read group comes back as identifier `0`
  ([`read_groups.md`](read_groups.md) §4). The parameters fit keys its per-library error rates on
  that identifier, so it refuses a cohort whose samples collide — which is every psp-mode cohort of
  two or more.

  **The fix is to renumber at the calling stage rather than at the walk.** Each psp carries its own
  table; the calling run reads them all at open and builds the run-wide numbering by merging them.
  *The alternative — numbering across the whole cohort before any sample is walked — was rejected
  because it destroys the property psp mode exists for*: the numbering would depend on the file
  list, so adding one sample later would renumber every psp already written and force a re-walk
  (§2).
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
- ~~How many records the file holds~~ — **dropped from the header (owner's ruling,
  2026-09-03).** The identity's second half travels *beside* the header rather than inside it:
  `PileupIdentity::of_header` already takes the count as its own argument, and the writer's
  `WriteStats` supplies it at the one moment a census is built from a walk. A header written
  before the first record cannot carry a count without being rewritten, and rewriting is what
  §6.3's invariance argument rules out.
- **When the file was written.** Never compared — it is what makes §12's byte-identity oracle
  read *"identical apart from the timestamp"* rather than *"identical"*, so it is a field the
  comparison deliberately skips rather than an omission.

**The census is the header's first consumer, and it is already built.** `PileupIdentity::of_header`
takes the psp header's bytes and a record count and says only this: two psps with the same header
and the same record count get the same identity and no others do — *"which bytes exactly is the
pileup writer's business"*
([`census_file.rs:91-98`](../../../../src/ng/parameter_estimation/joint/census_file.rs)). So the
encoding spec (§10) is free to lay the header out as it likes, and is not free to leave out any of
what §6.1 lists: the analysed regions, the catalog identity, the routing criteria, the read filters
and the command line.

### 6.2 The refusals

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
- **Across the cohort: two psps may not name the same sample.** Two files for one sample is either
  a duplicated argument or a cohort that would call one individual twice and weight the allele
  frequencies by it. The refusal names the sample and both files.
- **Each file against the parameters: every sample the parameters name must be present, and the
  reverse.** The per-sample values are keyed by name, not by position
  ([`parameters_file.md`](parameters_file.md) §6), so the match is by name and a gap either way is
  refused naming the samples. **This is what makes the file order irrelevant**, which is the point
  of recording the name at all.

**The read-group tables are merged rather than compared.** Two samples that numbered their read
groups from zero are the normal case, not an error (§6.1), so the calling stage builds one run-wide
numbering from the tables it read and every sample's identifiers are remapped into it. What *is*
refused is a table that cannot be merged — two read groups in one sample sharing an `@RG ID`, or a
sample whose table is absent — because neither can be renumbered without guessing.

### 6.3 The header does not describe where the psp's blocks fall — **ruled, 2026-09-03**

**Decision: the header carries no digest of the file's block boundaries, and no `writer_version`
field.** Confirmed by the owner 2026-09-03, with the reason stated stronger than this section
had it: a cohort whose files cut their blocks in different places **must not be refused**. In
the production caller, psps cut at the same points were more memory-efficient; whether that
holds in ng is unmeasured, and either way alignment is a property files may have — the
coordinate-grid cut (`psp_file_format.md` §4.1) gives it to same-block-size files by
construction — never one a run demands.

**What such a digest would be for.** If the calling stage merged the cohort in lockstep — block *k*
of every sample covering the same stretch of genome, so no sample's data ever waited for another's
— then two psps whose blocks fell in different places could not be merged, and a check at open
would catch it. **This design has no such consumer.** Each sample's reader independently serves the
segments it is asked for (§3.3), and the merge is keyed on coordinates with its own frontier
(§3.2); no code path ever looks at two samples' blocks together. Whether two files' blocks line up
changes nothing a run can observe.

**A check that guards nothing is a refusal waiting to fire on a harmless difference** — two writer
versions in one cohort, say — so it is not kept for safety's sake.

Three things follow:

- **`writer_version` goes with it.** Its only stated job was to explain a boundary-mismatch
  refusal. Versioning the psp *format*, so a reader can refuse bytes it does not understand, is a
  different job and the encoding spec (§10) owns it.
- **Block sizing becomes wholly the encoding spec's trade** — including a cap on how much a block
  decodes to, which matters at the top of the depth range: the densest stretch of a 300× sample
  would otherwise decode to an unbounded payload. Production's 1 MiB mid-window force-flush
  ([`src/psp/writer.rs:72,289-296`](../../../../src/psp/writer.rs)) shows the shape of that cap.
- **One property survives the hand-over, because an oracle depends on it:** where a block is cut
  must be a function of the observation stream alone — a stream that is identical at every worker
  count — so that the psp stays byte-identical across worker counts (§12.1). A cut that depended on
  scheduling would break that; a cut that depends on the sample's own data does not.

**If a digest is ever wanted back, the case for it has to name a consumer of boundary equality**,
and this design has none.

---

## 7. Where the memory goes

### 7.1 The three bounds

*One term is now measured — a reader's working set, at 480 kB on a human reference (§7.2). The
rest is arithmetic over estimates, and §11's question 2 is the measurement that replaces them.*

| object | peak resident |
|---|---|
| the walk stage | `samples in flight × (one open alignment file at 11–15 MiB + census accumulator at ~6 MB per read group + one sample's walking set)` |
| `PspVariantCaller` | `samples × one reader's working set` + **the merge's observation cache** + `callers in flight × one cohort locus` |
| `AlignedFilesVariantCaller` | `samples × 11–15 MiB` open alignment files + the merge frontier + `callers in flight × one cohort locus` |

**Only the first is independent of the cohort size**, which is the reason psp mode exists — and it
now depends on *samples in flight* rather than on walking one sample at a time, so the property
survives the change of default but is bought with a knob rather than for free.

**⚑ The observation cache is the term this table cannot price, and §3.2 understates it.** That
section says the merge's resident set is "about one observation per sample". The built merge holds
more: the cache spans the ground the builders in flight cover, and
[`cohort_merge.md`](cohort_merge.md) §8 measures **33 records held per sample** at 1,000 samples and
16 regions in flight on the 200-base default — against 4 at the old 20-base one. That document
declines to give a total, and gives its reason: **what one record costs in bytes is unmeasured**,
and at 300 reads a position a record is at its largest. So this row is a lower bound until that
number exists, and §11's question 2 owes it.

Pricing the `PspVariantCaller` row at the far end of the committed range — three thousand samples
on a human reference, **with both costs now measured on files ng writes** (§7.2) rather than
estimated. One open file plus the one cursor walking it is **480 kB**, so the cohort's floor is
**1.44 GB**, plus a caller pool's loci, which are negligible beside it. That is inside §7.2's
budget of 500 kB a sample, with 4% to spare.

**The look-ahead does not multiply that**, and this is where the change of unit (§3.1, §5.3) pays:
one source per sample for the whole run means the look-ahead multiplies the caller pool's loci,
not the readers. *(The measurement's own write-up prices a look-ahead of 8 at 4.0 GB, adding one
123 kB cursor per sample per unit. That arithmetic belongs to the source-per-in-flight-segment
design this document replaced; the per-cursor figure it rests on is kept, the multiplication is
not.)*

**The largest lever is not the look-ahead's**: three quarters of the 480 kB is a copy of the
reference's contig list, one per open sample, and §7.2 says what to do about it. It replaces the
estimate this paragraph used to carry — a 16 kB read buffer plus roughly 10 kB of decoded
observations, 26 kB a sample — which priced the same run at 620 MB, low by a factor of eighteen.

### 7.1a Memory is not the only per-open-file resource

**Three thousand open files is a file-descriptor limit before it is a memory bill**, and nothing
else in this document says so. Linux commonly ships `RLIMIT_NOFILE` at 1,024 and macOS at 256, so
the calling stage hits `EMFILE` around the thousandth sample whatever §7.2's budget says. Direct
mode reaches it sooner: a CRAM and its index are two descriptors each, so three thousand samples is
six thousand.

**What a run must do about it is refuse with a message that names the limit and the sample count**,
in the same shape §5.1 asks for when the memory bill does not fit — not die at file 1,020 with an
operating-system error. Raising the limit is the operator's to do, and the message should say so.
**This constrains nothing about the encoding**; it is a check at construction, beside the header
checks of §6.2.

### 7.2 Requirement: 500 kB resident per open psp — **reset 2026-08-25**

*This section asked for "tens of kilobytes, not megabytes" until the owner reset it on 2026-08-25.
The change is recorded rather than edited away because the encoding work was already running
against the old figure and one of its design choices turns on which number applies
([`../impl_plan/psp_encoding_experiments.md`](../impl_plan/psp_encoding_experiments.md), its
milestone B0).*

**The budget: 500 kB resident per open sample** — 500 MB across a thousand, 1.5 GB across three
thousand. The encoding spec sets the same number as its first goal
([`psp_file_format.md`](psp_file_format.md) §1.1), against that 1.5 GB. *(The plan that records the
reset says the owner priced "450 MB across a thousand samples" as comfortable, which is 450 kB a
sample; the 10% gap between the two figures is unexplained and neither document resolves it. Treat
500 kB as the budget and 450 MB as the comfort the owner actually expressed.)* **It was a working
figure rather than a ruling** when it was set, because the encoding's sweeps report the whole curve
of file size against reader memory — but the store has since been built and measured against it,
and the answer is below: **480 kB, inside the budget with 4% to spare**.

**Why there is a budget at all.** Everything in the calling stage multiplies by the sample count,
and the per-open-file state is the easiest cost to get wrong because it looks like bookkeeping.
Production's psp index is the counter-example ng would otherwise copy: a flat vector of one 24-byte
entry per block (`BlockIndexEntry`, [`src/psp/index.rs:42`](../../../../src/psp/index.rs)), decoded
whole at open (`decode_index`, [`:110`](../../../../src/psp/index.rs)). At a 5 kb block over an
800 Mb genome that is 160,000 entries — **3.8 MB per open file, 11.5 GB across three thousand
samples** — before any data is read. That is over budget by sevenfold even at the looser figure.

**The shape of the fix is now contingent, not settled.** A coarse index with the blocks chained
within it — each carrying enough to reach the next — is what this section used to prescribe. The
encoding plan declines to build it on speculation: how many blocks a file has at each block size is
being measured first, and *"if large blocks make it small enough, the coarse-index-and-chain scheme
`run_streaming.md` §7.2 asks for is not needed and should not be built"*
([`../impl_plan/psp_encoding_experiments.md`](../impl_plan/psp_encoding_experiments.md)). So this
section asks for the budget and no longer names the mechanism.

**And the measurement below has since answered it: the block index is not what costs.** Of an open
sample's 480 kB, the header, block index and footer together are 357 kB, and almost all of *that*
is the reference's contig list. The coarse-index-and-chain scheme this section used to prescribe
would be aimed at the wrong term.

**What the reset changed is not this requirement but a choice underneath it.** At tens of kilobytes,
writing each record's fields together was the only shape that fitted; at 500 kB, gathering like
fields together is affordable again and is smaller on the measurements taken so far. Which shape
wins is the encoding's to settle and this document does not care, provided the budget holds and the
cheap read of §3.3 survives.

#### What it actually costs, measured

**480 kB an open sample on a human reference, and 108 kB on tomato.** One store was opened 1, 2,
4, 8, 16, 32, 62, 125, 250, 500, 1,000 and 5,000 times over, every reader walked a record a round
in lockstep, and peak resident taken against the sample count — least squares, R² = 0.99999
([the measurement](../../reports/implementations/ng_psp_h4_2026-08-30.md)). The budget is met
with 19.7 kB, or 4 %, to spare. It is two parts, and only one of them is the reader:

| | human, 2,580 contigs | tomato, 13 contigs |
|---|---:|---:|
| the open file, before a block is touched — header, block index, footer | 357 kB | 7 kB |
| the cursor walking it — two 16 kB buffers, the decoder, the record being built | 123 kB | 101 kB |
| **per open sample** | **480 kB** | **108 kB** |

**The cursor costs near enough the same on both** — 123 kB against 101 kB — on corpora whose read
depth differs by a factor of 27 (280.0 reads a record against 10.3) and whose contig count differs
by a factor of 200. That is the *does not grow with the depth* half of the requirement, measured
rather than assumed.

**The header is the whole of the difference, and almost all of the header is the reference's
contig list**: about 138 bytes a contig, so 2,580 contigs cost 357 kB where 13 cost 7 kB. **Three
quarters of the human figure is a list that is identical in every sample of the cohort**, kept
once per open sample. A reference of about 3,700 contigs would spend the entire 500 kB on the
header before a record was read.

Giving a run's readers one copy of that list instead of one each takes the human figure from
480 kB to 123 kB, and it is the largest memory lever the store has. **The owner ruled on
2026-08-30 that it is not a question for the psp format**: the header goes on carrying the list,
because a psp has to be interpretable on its own, and the copy a reader works from can come from
the code that already handles the fasta reference. It is this document's to arrange — §10 carries
it — because this document owns the run objects.

---

## 8. Traps — what will bite the coder

Each is a property of code that exists today, and each produces a wrong answer or a silent
under-count rather than a crash.

- **The calling scratch is per worker and reused across loci, and a pool of per-locus callers makes
  a missed reset non-deterministic.** `CallingScratch`
  ([`src/ng/calling/mod.rs`](../../../../src/ng/calling/mod.rs)) is every buffer a locus's calling
  fills, allocated once per worker and cleared between loci. The code already records that a
  dropped `clear()` is invisible to the existing tests because it only shows in one locus order
  ([`inference/repeat_tract_parameters.rs`](../../../../src/ng/calling/inference/repeat_tract_parameters.rs),
  [`inference/summarise_condition.rs`](../../../../src/ng/calling/inference/summarise_condition.rs)).
  **With a pool of per-locus callers the order a given scratch sees is a scheduling artefact**, so
  a missed reset makes the VCF depend on the worker count — and only sometimes, which is worse than
  always. §12.2's invariance oracle is what
  must catch this, and it must therefore be run at more than one worker count on a fixture whose
  loci differ in kind.
- **A locus generator holds state across segments, so it cannot be shared between workers.** The
  iterator that owns the generators documents a load-bearing drop order
  ([`src/ng/locus_generation/mod.rs:917-942`](../../../../src/ng/locus_generation/mod.rs)). One
  generator set per source, and a source belongs to one worker: in a caller that is one set per
  sample for the run, and in a gatherer splitting one sample across workers, one set per worker,
  reused across *its own* consecutive segments and never shared.
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
  ([`src/ng/locus_generation/mod.rs:917-942`](../../../../src/ng/locus_generation/mod.rs)). A
  task that dismantles the pieces itself must keep the order.
- **A backward jump is legal and costs a seek plus a block decode.** The read cursor answers
  segments in any order, and a test asserts backwards-walked segments return what a linear scan
  returns ([`src/ng/read/input/cursor.rs:92-96`](../../../../src/ng/read/input/cursor.rs), test
  at [`:1207`](../../../../src/ng/read/input/cursor.rs)); the psp reader's coarse index (§7.2)
  has the same character. So a work-stealing pool is correct however it steals, and slow unless
  each worker's own sequence of segments stays monotonic — schedule for that.
- **Every read is named at every position it covers, reference-matching reads included — and the
  psp writer must not "optimise" that away.** The opposite is what a coder would reach for: production's reference-side read names were
  about 31% of its peak live heap, so dropping them looks like free money. **They are not
  droppable.** Both generic paths record the ids of every read folded at a position, whether it
  departed from the reference or agreed with it — the owner's ruling of 2026-08-17
  ([`open_record.rs:472-473`](../../../../src/ng/locus_generation/pileup/open_record.rs),
  [`fast_column.rs:120-126`](../../../../src/ng/locus_generation/pileup/fast_column.rs)) — because
  a cohort locus can span several of one sample's records and the merge needs to know which reads
  are the same read across them. The psp encoding plan lists this under *"one thing that is settled
  and must not be reopened"*
  ([`../impl_plan/psp_encoding_experiments.md`](../impl_plan/psp_encoding_experiments.md)); what it
  does instead is encode the names as arrivals and departures rather than repeating them.
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
generator sets, and each reader's own working set. Single-threaded by construction: the yield point — each object
is an iterator, so its consumer (the VCF writer, the psp writer, the census accumulator) runs on
the consuming thread and needs no lock.

**Performance.** Three knobs at this level: **samples in flight** in the walk stage (§5.2),
**callers in flight** on the calling side (§3.5), and the within-sample walk's worker count, which
is off by default (§5.2). **The cohort merge owns two more of its own** — the building region's
width and how many are in flight
([`cohort_merge/mod.rs:531,592`](../../../../src/ng/run/cohort_merge/mod.rs)) — which §3.5 and
§12.2 both treat as settings, so five in total and not three. Block size is the encoding spec's trade (§6.3) — it reaches this design
only as one reader's working set, the multiplicand in §7.1. Any other tuning constant that appears in
the implementation is a defect, not a lever. **None of the three has a proposed default**, because
nothing about the psp format's cost, the calling loop's cost, or either one's memory has been
measured yet (§11, questions 2 and 7).

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
- **A way to ask a source the cheap question.** §3.3's whole saving rests on reading two small
  numbers per position and inflating the rest only where the cohort kept something. **The built
  merge cannot express that**: `ObservationSource` has one method, which returns a whole inflated
  observation ([`cohort_merge/observation_cache.rs`](../../../../src/ng/run/cohort_merge/observation_cache.rs)),
  so a psp reader behind it would decode everything for every position. The shape of the fix is a
  second method — *the cheap numbers over this range* — with the existing one becoming *now inflate
  these positions*; the alternative is a scan pass ahead of the build pass, which costs a second
  traversal. **It belongs with the psp reader**, since a walker's cheap answer is free and only a
  file reader has anything to gain, and it cannot be settled before the encoding is.
- ~~**The cohort merge's reconciliation**~~ — **written, 2026-08-17:
  [`cohort_merge.md`](cohort_merge.md), and built since** (`src/ng/run/cohort_merge/`). It settles
  how per-sample observations whose spans differ become one cohort locus — every member
  projected onto the locus span and unified into one allele table (its §4.2) — and discharges the
  confirmation this entry asked for, that the merge frontier stays bounded while spans are
  reconciled (its §4.5). It adds one rule this document should know about: **a cohort locus wider
  than `max_cohort_locus_span` is never built**, and is counted as a failed locus rather than sent
  downstream.
- ~~**The parameters file's format**~~ — **written, 2026-08-28:
  [`parameters_file.md`](parameters_file.md).** A TOML file carrying every number calling runs on,
  each with its warrant. Two things in it change this document: **every run writes one beside its
  VCF**, whatever the numbers came from (that spec §7), and **the defaults live in the binary
  rather than in a shipped file**, so "run without a fit" is a flag and not a path (§8 there).
- **One contig list for the run, not one per open sample** — on a human reference 357 kB of an
  open sample's 480 kB is a copy of the reference's contigs, identical in every sample (§7.2), and
  at three thousand samples that is 1.07 GB of the same list. Its home is here, in the object that
  opens the files: `PspVariantCaller` construction already reads every header, so it can check the
  lists agree and then hand every reader one shared list — sourced from the code that handles the
  fasta reference, which the run holds anyway. The psp format does not change: a psp still carries
  its own contig list, or it stops being interpretable on its own (owner's ruling, 2026-08-30).
- **The VCF writer** — consumes a caller's iterator; the variants arrive in genome order, so it
  writes as it reads. Its shape, and the `Variant` record's, belong to the emission step's
  document.

---

## 11. Open questions

1. **How long is a segment at the routing floor?** — OPEN. Measured only at the catalog floor:
   391 bases average (§4.4). *Leaning:* kilobases at the routing floor. **Settled by:** counting
   segments and their length distribution over the existing catalog file at both floors, tomato
   and human — a filter over a stored file, not a genome scan.
2. **The two concurrency defaults: samples in flight for the walk, callers in flight for the two
   callers.** — the walk half is **ANSWERED** (owner, 2026-09-03): `generate-psps` processes its
   samples one at a time, and a cohort is parallelised by running invocations — no in-process
   fan-out and no default owed (§5.2). The callers-in-flight half stays OPEN; no value proposed —
   a caller in flight costs one locus. **Settled by:** sweeping it on the tomato slices and on
   HG002, wall time and peak resident, with the output required identical at every setting.
3. **Does splitting one sample's walk across workers scale?** — OPEN. It is not on the default
   path, and it is the measurement question 8 turns on. **Settled by:** driving one gatherer at 1, 2, 4, 8, 16 workers on a tomato slice
   and HG002 — wall time, peak resident, observations identical to serial. Production's is 1.81×
   at four threads and worse at eight.
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
7. **Which stage of the psp-to-VCF path is worth a pool at all?** — OPEN, and it comes before any
   of the improvements below. The candidates are decoding the psps, merging, and calling. **Settled
   by:** a profile of one cohort run over a handful of tomato slices, wall time split three ways.
   Production's shape — posterior engine about 3%, producer decode about 30% — suggests decode
   rather than calling, but production's engine is not ng's loop.

   **The direct-mode half is answered, in two findings.** The merge half
   ([`../research/cohort_merge_parallel_cost_2026-08-28.md`](../research/cohort_merge_parallel_cost_2026-08-28.md)):
   eight threads give **3.1×** on 63 tomato accessions at 1,000-base regions — not the 1.4×
   this question used to quote, which came from 200-base regions at 16 samples — and the merge
   is 1.4–10% of walking-plus-merging, so it is already worth its pool and is not where a
   run's time is. The run half
   ([`../research/cohort_merge_parallel_cost_2026-09-03.md`](../research/cohort_merge_parallel_cost_2026-09-03.md)):
   at 8 threads a 63-sample run over 200 kb is **84.5% decoding reads, 8.0% genotyping, 6.2%
   assembling** — the pool belongs to the decode and nothing else earns one yet. The decode is
   already spread across samples; what caps it at 2× is the work itself slowing ~2.6-fold when
   eight copies run at once (the choice of allocator, the barrier, the fixpoint's re-sweeps
   and the frees are each refuted by measurement there) plus a ×1.5 per-sweep wait for the
   slowest sample. Splitting cache contention from allocation traffic from the Mac VM's
   scheduling needs `perf` on the Linux box, and that is where this question now lives. The
   psp half stays OPEN until the psp format exists.

   **⚠ An instrument regression sat between those two findings** (the 09-03 finding §4): the
   two per-record counters G2 added to the `merge-timing` feature on 2026-09-01 shared one
   cache line across every worker, and an instrumented parallel merge measured between then
   and 2026-09-03 reads ~1.4× where the true figure is ~2.2× (200-base regions). They are
   sharded now and the instrumented build matches an uninstrumented one again.

   **The merge half's plan**
   ([`../research/cohort_merge_parallel_cost_plan.md`](../research/cohort_merge_parallel_cost_plan.md))
   is settled by the two findings above. The sketch below is kept for the record; of it, the
   width extension was done (500 is the default now), the overlap driver was built and
   refuted, and the owner dropped the sliding window once the barrier priced at ~4%:

   - **Extend the building-region width sweep** past eight threads and toward the far end of the
     cohort range; it is a run parameter, so it costs no code, and it has already been swept at one
     and eight threads on 63 samples, where 100–200 bases was the optimum and 200 is the shipped
     default. What is untested is whether that optimum holds past eight threads, and whether it
     moves as the cohort grows — the organiser's per-region work walks the whole cohort, so its
     share should grow with sample count.
   - **Overlap the reader advance with the building.** The organiser draws the readers forward
     while no builder runs, and that phase is I/O and decode — exactly the work worth hiding
     behind compute. Filling the next round's ground into a second buffer while the current
     round's builders run removes the serial phase; the cost is holding two rounds of observations
     instead of one, and the released order is unchanged because release is still by region index.
   - **Drop the rounds for a sliding window.** The barrier exists because builders read the
     observation cache while the organiser writes it, not because of ordering — the organiser
     already releases along a gapless run of region indexes, which supports arbitrary out-of-order
     completion. A cache that is extended ahead of the claim frontier and evicted behind the
     release frontier, with a published "covered to here" mark that builders never read past,
     removes the barrier entirely. It is the largest of the three and should not be attempted
     before the first two have said whether it is needed.
8. **How does a run of one sample use the machine?** — OPEN, and it is goal 3 unmet. Walking
   several samples at once gives a lone sample one worker. The only candidate anybody has is
   question 3's — split that sample's walk across its segments — and if that does not scale there
   is no second idea on the table, so the honest state is *goal accepted, mechanism unknown*.
   **Settled by:** question 3's sweep, which decides between the candidate working and the
   question being genuinely open.

   **Two things follow from the answer, which is why this is not idle.** If the split works, the
   psp writer must go on honouring §12.1's worker-count invariance, because observations then
   arrive out of order within one sample. If it does not, that requirement has no remaining
   consumer and the writer gets simpler. **Nothing about the writer should be simplified on the
   strength of the current default alone** — the default is a schedule, and the requirement is
   about what schedules must remain possible.

---

## 12. How we know it works

Each oracle is a property of the run, not of one type — which is why they live here.

1. **Worker-count invariance of the psp.** One sample gathered at 1, 2, 4, 8, 16 workers gives
   byte-identical psps apart from the header's timestamp. **Production already holds this** — its
   `.psp` bodies are byte-identical across every thread count, the `created` timestamp being the
   only difference
   ([`pileup_thread_scaling_2026-06-11.md`](../../reports/pileup_thread_scaling_2026-06-11.md)) —
   so it is a reachable bar — and it is what §6.3's
   restriction on block cuts preserves.
2. **Concurrency invariance of the VCF**, from each of the two callers: the same VCF at any number
   of callers in flight, and — where the merge's own region parallelism is switched on — at any
   building-region width and any number of regions in flight. The merge already holds the second
   half against its single-threaded oracle
   ([`parallel.rs`](../../../../src/ng/run/cohort_merge/parallel.rs)).
3. **Mode equivalence — the oracle that justifies the design.** The same cohort and the same
   parameters, run through `AlignedFilesVariantCaller` and through the psp route, give the same
   VCF. This is simultaneously the proof that the calling function is mode-blind (goal 1) and the
   sufficiency test for the psp: anything the file fails to carry surfaces here, where a
   write-read round-trip test would pass.
4. **Segment independence of the observations.** One segment walked alone emits exactly the
   observations the same span emits inside a whole-genome, single-threaded walk. This asserts
   §4.3 is honoured; the thirds-chopping test that lost 17 positions is the failure shape it
   catches.
5. **Every refusal fires at construction, and names what differs.** Unequal analysed regions name
   both samples; differing segmentation inputs name the first differing field; two psps for one
   sample name the sample and both files; a parameters file whose sample list does not match names
   the samples missing from each side (§6.2). None of them may wait for a locus.
6. **The file order does not matter.** The same cohort called with its psp arguments in any order
   gives the same VCF, sample for sample — the test that the name-keyed join of §6.2 is real and
   that nothing joins by position.
7. **A cohort of separately-walked samples is callable.** Every sample walked in its own
   invocation, each numbering its read groups from zero, merges into one run-wide numbering and
   calls (§6.1, §6.2). Without this the ordinary psp-mode run does not work at all, and the failure
   is a refusal at the fit rather than anything visible in the walk.
8. **The census built during the walk equals the census built from the psp.** Specified already
   as a byte-for-byte comparison
   ([`parameter_prepass_joint_records.md`](parameter_prepass_joint_records.md) §7.12) and named
   here because this document is what puts the two producers in one run: it is the sharpest test
   that a psp holds everything a census needs.
9. **Analysed-but-empty survives a round trip.** A stretch a sample analysed and found no reads
   in reads back as analysed and empty — distinguishable, via the header's analysed regions, from
   a stretch the sample never looked at (§8, last trap).

---

## 13. Reuse map

| what | existing code | how it is reused |
|---|---|---|
| the segments the loops run over | `TypedRegion`, `RegionKind` ([`src/ng/region_typing/mod.rs:144,168`](../../../../src/ng/region_typing/mod.rs)) | consumed as-is; the segments a source is asked for, one at a time, are these |
| the analysed regions | `GenomeRegions` ([`src/ng/region_typing/mod.rs:77,87,100`](../../../../src/ng/region_typing/mod.rs)) | reused whole — it wraps production's `RegionSet`, so ng and production agree on what a BED means; its value is recorded in the psp header |
| what a BED edge does to a segment | `clips_at_a_bed_edge` and the emission rule ([`src/ng/region_typing/mod.rs:471,482-488`](../../../../src/ng/region_typing/mod.rs)) | taken as given — findings whole, generic clipped; nothing here re-decides it |
| one sample's observations | `SampleLocusObservations` ([`src/ng/locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | the item of every stream in §3, unchanged |
| the walker behind a source | `SampleLocusObservationsIterator` ([`src/ng/locus_generation/mod.rs:915`](../../../../src/ng/locus_generation/mod.rs)) | one per task, fed that segment's segments |
| per-segment reads | `SampleReads` and `cursor()` ([`src/ng/read/input/mod.rs:398,623`](../../../../src/ng/read/input/mod.rs)) | one shared `SampleReads` per sample, one owned cursor per worker (`Send` proven; `Sync` to confirm — §8) |
| the streaming merge's shape | `MergedCursors` ([`src/ng/read/input/sample_cursor.rs:178`](../../../../src/ng/read/input/sample_cursor.rs)) | the model for `k_way_merge`: argmin over per-stream heads, keys beside the heads, frontier-sized residency. *(`arch/sample_reads.md` §4 calls this `MergedRegionReads`, a name from a superseded design that is not in the tree.)* |
| the census accumulator | `CensusWriter::add_locus`, `mark_walked`, `finish` ([`src/ng/parameter_estimation/joint/census.rs:2087,2106,2378`](../../../../src/ng/parameter_estimation/joint/census.rs)) | fed inside the gatherer, at the ordered yield point (§5.2) |
| the census file | `write_census`, `open_census` ([`src/ng/parameter_estimation/joint/census_file.rs:200,426`](../../../../src/ng/parameter_estimation/joint/census_file.rs)) | written by the walk stage's per-sample loop from `finish()`'s result |
| psp block index | `BlockIndexEntry`, `decode_index` ([`src/psp/index.rs:42,110`](../../../../src/psp/index.rs)) | **a model of what not to build** — §7.2 rejects the flat per-block index at ng's sample counts |
| block cutting | `PspWriter`'s grid and force-flush ([`src/psp/writer.rs:297-301,72,289-296`](../../../../src/psp/writer.rs)) | neither carries as a rule: block sizing is wholly the encoding spec's, and the force-flush is now a legal shape for capping a block's decoded size (§6.3) |

**The parity oracle for the whole document is §12.3** — direct mode against psp mode, one cohort,
parameters held fixed.
