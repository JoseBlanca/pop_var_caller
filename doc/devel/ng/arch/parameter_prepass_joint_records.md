# ng — the joint parameters fit, what is recorded at each kept locus: types & interfaces

*Status: architecture draft (2026-08-12), companion to the spec
[`../spec/parameter_prepass_joint_records.md`](../spec/parameter_prepass_joint_records.md) (the
design and its rationale) and to the shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md)
(vocabulary) and [`module_layout.md`](module_layout.md) (the `src/ng/` tree). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for types, verbs
for functions, **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the **contract** is the
deliverable. Every "why" lives in the spec — this doc does not re-argue one.*

*Which loci these are is [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md); the
mathematics that reads them is [`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md).*

## Module home

`src/ng/parameter_estimation/joint/records.rs`, beside `loci.rs` and `fit.rs`. **No trait**: one
record shape per path, no competing candidate, so a codec trait would be ceremony over two structs
and a writer.

---

## 1. Types

### 1.1 `SampleRecords` — one sample's evidence, however it is being held

**One type for both kinds of run.** A run that goes straight from alignments to a fit holds every
section in memory; a run reading a file fills a section when it is asked for and drops it after. That
is the *only* difference between them, so it is a property of the value rather than a second type —
and the parameters fit is written once, against calls that behave the same either way.

```rust
pub struct SampleRecords {
    pub sample: SampleName,
    /// The thirteen values the parameters fit refuses to pool across (spec §5), and the
    /// digest of the loci actually kept. **Outside the sections**: they are compared
    /// across every sample before anything large is decoded.
    pub identity: RecordIdentity,
    /// Resident, or backed by a file and loaded on demand. Not public: what a caller may
    /// do with a section is §2.2's scoped access, and a field would let it keep one.
    sections: Sections,
}
```

**Contract.** A sample's counts at a position are the sum of its read groups' — exact, so the
parameters fit may fold freely. Sections are enumerated in a fixed order (read groups by id, strata by
their stratum key), so a fit that iterates is deterministic.

**What was here before, and why it is gone.** An earlier version of this section held
`generic: BTreeMap<ReadGroupId, GenericRecords>` and `ssr: BTreeMap<ReadGroupId, SsrRecords>` as public
fields, and a companion `SampleRecordsFile` for the run that could not hold them. Two shapes for one
object made the fit two code paths, and public maps made *"hold one section at a time"* a convention
rather than a property — nothing stopped a caller keeping every stratum it had ever asked for. **The
storage stays exactly as §1.2 and §1.4 describe it**; what changed is that no caller reaches it
directly.

*Rejected, and it is worth saying why because it is the natural first idea:* one flat collection of
per-locus records, `Vec<CensusRecord>` with an enum for the two kinds. **A Rust enum is as wide as its
widest variant**, and the repeat-tract variant carries nine two-byte offset buckets before its other
fields — so each of the two million ordinary positions would pay about 20 bytes where five bits is what
it needs, some 40 MB a read group against the packed array's 1.25 MB. Boxing the variants trades the
width for a pointer and an allocation per locus. And one collection in genome order interleaves the two
halves, so reading one stratum would mean walking all of it. **`CensusRecord` survives one level up**,
as the item the genome walk emits and the fit iterates (§2.2.1), where the uniformity is worth having and
the storage is not paying for it.

### 1.2 `GenericRecords` — a dense array and a sparse list

Two parts, because at three reads a site nearly every position is *"n reads, all matching the
reference"* and only a few thousand in a million are anything else. The dense array holds the depth of
every kept position in order; the sparse list holds the exceptions, each naming the position it
belongs to by its **index into the selection** rather than by a coordinate.

```rust
pub struct GenericRecords {
    /// Entry `i` is the `i`-th kept position's depth code, five bits, packed.
    depth: PackedDepthCodes,
    /// Only where a read was not on the reference base. **Sorted by `index`**, so the
    /// parameters fit walks both halves in one pass.
    non_reference: Vec<AlleleObservation>,
}

#[derive(Copy, Clone)]
pub struct AlleleObservation {
    /// Index into `KeptLoci::generic` — the position's only identity.
    pub index: u32,
    pub allele: ObservedAllele,
    /// **Exact, not binned** (spec §2.2): nearly every entry is a miscall of one to three
    /// reads, so a ladder's compressed tail stays empty at any depth, and binning would
    /// cost a fourteenth bin width in `RecordIdentity` for nothing.
    pub reads: u32,
}

/// A, C, G, T, or anything else — an indel or a spanning deletion (spec §2.1).
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum ObservedAllele { A, C, G, T, Other }
```

**Contract.** The dense half reconstructs a quiet position exactly: a code with no sparse entry means
"that depth, all reads on the reference base", and the reference base belongs to the position.

### 1.3 `DepthCode` — the ladder plus one state it cannot express

A depth bin alone cannot say *this position was never visited*, and that state has to be
distinguishable from *visited and empty* because only the first is a bug. So the stored code is the
ladder's bins plus one sentinel — 21 codes, which is why five bits rather than four.

```rust
#[derive(Copy, Clone, PartialEq, Eq)]
pub enum DepthCode {
    /// The region holding this position was never walked — a bug, not data (spec §1.1).
    NeverWalked,
    /// A rung of the shared ladder, bin 0 being zero depth.
    Binned(DepthBin),
}
```

**Contract.** Five bits per entry, asserted at build time (`const _: () = assert!(…)`) so a ladder
grown past 31 codes fails to compile rather than truncating. **The ladder is not defined here** — it
is `DepthBinEdges`, shared with the histogram route so the two cannot bin differently, which is what
makes their comparison a comparison.

### 1.4 `SsrRecords` — lengths, a guard, and the differences

At a repeat tract the observation is a tract *length*, so the allele record is a distribution over
offsets from the reference tract length. Three things travel beside it, and each answers a question
the offsets cannot: how many reads reached the locus without crossing it (a lower bound on the
length), how many differed by something that is not a whole number of copies (the guard), and which
bases mismatched and where (the substitution channel — the error rate cannot be recovered from
lengths).

```rust
pub struct SsrRecords {
    /// Spanning reads at each whole-repeat offset from the reference tract length,
    /// ±4 with saturating ends.
    offsets: Vec<OffsetCounts>,
    /// Reads that reached a locus and crossed no whole tract — a censored lower bound
    /// (spec §3).
    covering_not_crossing: Vec<u16>,
    /// Whether the genome walk reached this locus at all. **The STR half's answer to the
    /// generic half's never-walked sentinel**, and it has to be its own field: the
    /// other four vectors are all zero both for a locus never walked and for one
    /// walked with no read, and only the first is a bug (spec §6).
    walked: Vec<bool>,
    /// Reads differing by a non-whole number of copies, with what they were.
    guard: Vec<GuardObservation>,
    /// Per locus, the denominator the STR error rate is fitted against.
    bases_compared: Vec<u32>,
    differences: Vec<TractDifference>,
}

#[derive(Copy, Clone)]
pub struct TractDifference {
    pub locus: u32,
    /// Which of this locus's reads carried it, in the locus's own read order.
    /// **Two entries at one offset on one read is a different observation from the same
    /// two on two reads** — a read-blind encoding passes every other check (spec §7.3).
    pub read: u8,
    /// Signed offset from the tract start: negative in the left flank, `0..len` inside,
    /// beyond `len` in the right flank.
    pub offset: i16,
    pub base: ObservedAllele,
}
```

**Contract.** All the per-locus vectors are indexed by the locus's position in `KeptLoci`, as the
generic half is indexed by position. An offset beyond ±4 saturates into the end bucket rather than
being dropped or wrapping.

**`SsrLocusState` is what a consumer asks**, rather than reading the vectors and re-deriving the four
states — never walked, walked with no read, reached but not crossed, crossed — each time. The guard's
one-in-ten threshold is a query on the same type (`guard_is_over_threshold`) and **not** a write-time
error: a locus over it is well-formed data the parameters fit should decline to fit, and refusing the write would
make a property of the sample look like a property of the file.

### 1.5 `RecordIdentity` — the thirteen values

What travels with a sample so the parameters fit can refuse to pool two runs that did not record the same thing.
Seven say which loci were *asked for*, one says which came back, and **five say in what units the
evidence was written down** — the group an earlier version of this section did not have, which left
the whole check able to pass while two samples' rows meant different things (spec §5).

```rust
#[derive(Clone, PartialEq)]
pub struct RecordIdentity {
    /// The seven that say which loci were asked for.
    pub selection: SelectionIdentity,
    /// The one that says which loci came back (`loci.rs` §1.3).
    pub kept_loci: KeptLociDigest,

    // ---- the five that say in what units ----
    /// Per stratum, held against kept: anything pooled across strata is biased without
    /// it and silently so (spec §5).
    pub ssr_stratum_counts: StratumCounts,
    /// Must match the other route's, or the comparison confounds the cap with the
    /// route (spec §3.4).
    pub read_cap: ReadCap,
    /// A digest of `DepthBinEdges`. **The generic record stores a code, not a depth**, and
    /// two samples binned under different edges hold codes meaning different things while
    /// every other value here agrees — the loci did match (spec §5).
    pub depth_ladder: DepthLadderDigest,
    /// The depth above which a position's reads are subsampled before anything is
    /// recorded — the generic twin of `read_cap`, and it moves independently of the
    /// ladder's top rung.
    pub depth_cap: DepthCap,
    /// The coverage-by-window grid, where a summary exists (§1.6). Windows of different
    /// widths are not comparable and a relative copy number across two grids is
    /// meaningless.
    pub coverage_window: Option<Bp>,
}
```

**`PartialEq` and not `Eq`**, forced by `SelectionIdentity` holding float floors — compare those by
bit pattern (`loci.rs` §1.2).

### 1.6 `CoverageByWindow` — the third object, and the one that is not a record

**The parameters fit's duplicated-site class is conditioned on the sample's local relative coverage, and that
cannot be derived from §1.2's records**: those hold one binned depth per kept position, and the kept
positions are one in a few hundred, so a 500 bp window holds one or two of them — the per-base
measurement the class's own constraint rules out (fit spec §2.2, records spec §4).

```rust
/// One sample's depth over fixed windows of the reference, plus the GC curve that
/// corrects it. **Over every position the genome walk visited, not over the kept ones.**
pub struct CoverageByWindow {
    /// The grid, a function of the reference alone — so two samples' summaries are
    /// comparable by construction. Travels in `RecordIdentity` (§1.5). **Always the
    /// stored 500 bp**; the width the parameters fit reads at is its own decision (spec §4.1).
    window_bp: Bp,
    /// The sample's median window depth, in reads a position. The one number that
    /// makes `depth`'s bytes mean something, and it is per sample deliberately: what
    /// the class reads is coverage *relative to this sample*.
    median_depth: f32,
    /// Mean depth per window, in reference order, as `round(32 × mean / median_depth)`
    /// saturating at 255. Resolution is 3% of the sample's own median and the range
    /// reaches 8 times it — a byte is enough only because the quantity is a ratio.
    /// **Storing the mean itself in a byte would not do**: at three reads a position
    /// the difference between one copy and two is the difference between 3 and 6.
    /// 1.6 M windows on tomato, 6.2 M on GRCh38 — 1.6 to 6.2 MB per sample (spec §6).
    depth: Vec<u8>,
    /// Windows holding fewer than `window_bp` walked positions, as
    /// `(window index, positions)` — the ends of contigs, and anything the analysed
    /// regions or the ambiguity mask cut into. **The parameters fit needs these to sum adjacent
    /// windows** (spec §4.1): a wider mean is `Σ depth × positions / Σ positions`, and
    /// a short window weighted as if it were full pulls the wider mean towards it.
    /// Sparse because nearly every window is full, exactly as §1.2's non-reference
    /// entries are.
    short_windows: Vec<(u32, u16)>,
    /// Depth against GC content, a few hundred numbers. Coverage tracks GC, and an
    /// uncorrected window at an extreme of it reads high for a reason that has nothing
    /// to do with copy number: on tomato the median window depth runs from 16.2 reads
    /// a position at 20% GC to 29.0 at 36%, a factor of 1.79.
    gc_curve: Vec<f32>,
}
```

**Contract.** Per sample and **not** per read group — copy number is a property of the individual, not
of a library, and eight tomato samples on one grid confirm it: of 84 windows some sample reads near
two copies, 40 are read that way by exactly one of them (spec §4.2). Needs no cohort, which is what
lets the caller run on one sample. **A window's mean depth is not the mean over the kept positions
inside it**, and spec §7.10 asserts that inequality because an implementation that quietly derived one
from the other would pass every other check.

**The parameters fit sums, the genome walk does not.** Adjacent windows are summed at read time up to the width the
sample's depth requires — about 12,000 aligned bases, so 500 bp at 25× and 5 kb at 2.5× (spec §4.1).
The type therefore exposes the denominators rather than only the means, and offers no resampling: a
summary built at a different `window_bp` is refused, not converted (spec §7.10).

**Settled, 2026-08-12** — the object is built. Fit spec §2.2's gating measurement returned 1 position
in 8,600 in a two-copy window reading near half, a near-half rate of 1.26% inside those windows
against 0.033% outside, and 24.8 times what independence would give
(`../reports/duplicated_locus_probe_2026-08-12.md`). `coverage_window` is `None` where a run chose not
to build the summary and **where a `SampleRecords` has just been read from a file** (§2.2, spec §4:
the summary is never serialized) — never because the class was not there.

---

## 2. Interfaces

### 2.1 Filling records during the genome walk

The writer is handed the same locus stream the histogram accumulators are handed, so one genome walk fills
both routes and the comparison between them is over identical evidence. It knows which loci are kept
and ignores the rest.

**Built, 2026-08-12.**

```rust
pub struct RecordWriter { /* … */ }

impl RecordWriter {
    pub fn new(
        sample: String,
        loci: &KeptLoci,
        read_groups: Vec<ReadGroupId>,
        contig_of: &dyn Fn(&str) -> Option<ContigId>,
        identity: SelectionIdentity,
        edges: DepthBinEdges,
        read_cap: ReadCap,
        depth_cap: DepthCap,
    ) -> Self;

    /// Record this locus if it is a kept one. **Borrows and does not take**, so the genome walk
    /// passes the locus on to the histogram accumulators untouched and one pass fills
    /// both routes.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations);

    pub fn finish(self, coverage: Option<CoverageByWindow>) -> SampleRecords;
}
```

**The read groups are declared, not discovered.** A group that put no read at a position must still
get its zero there, because that entry is the denominator its own error rate is fitted against;
discovering the groups from the observations would start each group's record at its first read and
leave every position before it indistinguishable from never walked.

**Depth is the read group's own and never the sample's.** `SampleLocusObservations::num_obs_along_locus`
pools, which is the wrong grain for a record keyed by read group — and it is the grain the error rate
is fitted at, so pooling would score every read against a rate fitted from a depth its own library
never had. The writer derives per-group depth from the observations itself, honouring the same
`ReadWitness` runs that method does.

**`contig_of` because a kept STR locus is an `SsrSegment`**, which names its contig by string where
a locus names it by index. One closure at construction resolves the whole selection once.

**No `merge` yet.** A region-sharded genome walk fills disjoint index ranges and merging is concatenation,
but nothing in this build shards, and a merge with no caller is a merge with no test that could fail.

**Contract.** Every kept locus gets an entry whether or not a read reached it — the entry is the
denominator. A locus in a region never visited keeps `DepthCode::NeverWalked`, so the three states
survive a write and a read. **`add_locus` is the only place the kept-loci digest is fed**, so a
record set and its digest cannot disagree.

### 2.2 Getting records to the parameters fit

**Decided 2026-08-13 (spec §6.1): one file per sample, written beside that sample's pileup and never
inside it.** The requirement is unchanged — the parameters fit reaches every sample's records without walking
reads again — but where they live is no longer open, so this section carries the shape rather than
deferring it.

```rust
pub fn write_records(records: &SampleRecords, out: impl Write) -> Result<(), RecordError>;
pub fn read_records(input: impl Read) -> Result<SampleRecords, RecordError>;
```

**Reading is by section, and a section is borrowed for the length of a call — 2026-08-13 (spec §6.2).**
The parameters fit holds the generic half and the repeat-tract half at different times, and one band of
strata at a time within the second. **Returning a section would make that a convention**: a caller
could keep every one it had asked for, and a file-backed value would grow into the whole file it was
supposed to avoid. So access is scoped, and holding two sections is not expressible.

**The unit is one stratum across every sample, not one sample's stratum.** A tract's length frequencies
are fitted from every sample with reads there (fit spec §4.1), so the fit needs sample 1 to *N* at
stratum *k* together and none of them afterwards. That puts the scoped call on the cohort rather than
on one sample:

```rust
/// Every sample's records, however each one is held.
pub struct CohortRecords { /* Vec<SampleRecords>, in sample-name order */ }

impl CohortRecords {
    /// Opens or adopts each sample's records and checks the thirteen identity values
    /// across all of them **before any section is decoded** (spec §5).
    pub fn new(samples: Vec<SampleRecords>) -> Result<Self, RecordError>;

    pub fn read_groups(&self) -> &[ReadGroupId];
    pub fn strata(&self) -> &[Stratum];

    /// Every sample's generic records for one read group, for the length of the call.
    pub fn with_generic<R>(
        &mut self,
        group: ReadGroupId,
        f: impl FnOnce(&[&GenericRecords]) -> R,
    ) -> Result<R, RecordError>;

    /// Every sample's tracts for a **band** of strata, for the length of the call.
    pub fn with_strata<R>(
        &mut self,
        group: ReadGroupId,
        strata: &[Stratum],
        f: impl FnOnce(&[&SsrRecords]) -> R,
    ) -> Result<R, RecordError>;
}
```

**Why a band and not one stratum.** Of tomato's 141 strata, 68 hold fewer than a hundred tracts each
and are fitted by borrowing from their neighbouring repeat counts rather than alone
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §3.6), so the fit will sometimes
want a stratum and its neighbours together. Those are the thin ones — a band of three costs almost
nothing — but a signature taking a single stratum would have to be widened once code depended on it.

**Contract.** A call decodes **its own bytes and no others**; spec §7.16 asserts the byte count and not
only the values, because an implementation that decodes the whole file and hands back a slice satisfies
every value comparison and delivers none of the memory this exists for. A file-backed sample drops what
it loaded when the call returns; a resident one lends what it already holds and drops nothing. **Peak
memory is samples × the largest band**: at a per-stratum cap of 5,000 tracts and about ten bytes a
tract, 50 kB a sample, so 50 MB across a thousand.

**`SsrRecords` becomes per (read group × stratum) rather than per read group.** Today it holds vectors
indexed by a locus's position in the whole kept set; a section holds one stratum's slice of that, so
the index is stratum-local and the stratum's first index travels in the file's directory. **That is the
one storage-type change this decision forces**, and it is the reason it is recorded here rather than
left to implementation.

### 2.2.1 `CensusRecord` — the item both paths iterate

The two halves store different things and are read through different calls, but a genome walk emits one
locus at a time and a fit consumes one locus at a time, and at *that* level they are the same shape.

```rust
/// One locus's evidence for one read group, as the walk emits it and the fit reads it.
pub enum CensusRecord<'a> {
    Generic { index: u32, depth: DepthCode, non_reference: &'a [AlleleObservation] },
    Ssr { index: u32, stratum: Stratum, tract: &'a SsrLocus },
}
```

**Contract.** This is an **iteration item, never a stored one** — it borrows from a section that is
already decoded, so a two-million-entry sequence of them exists only as a cursor. §1.1 says what
storing them instead would cost.

**The coverage-by-window summary is not in that stream.** It is recoverable from the pileup and the
owner's decision is not to keep a copy (spec §4), so `write_records` omits it and `read_records`
returns a `SampleRecords` whose `coverage_window` is `None`. Whoever runs the parameters fit fills it: from the
walk in the direct run, and from a pass over the pileups in the two-phase one. **The type keeps the
field** — it is what the parameters fit reads — and only the serialized form drops it.

**Two producers, one builder.** The genome walk-time producer is `RecordWriter` (§2.1). The second reads an
existing pileup and drives the same `RecordWriter` through the same locus stream, so there is one
implementation of what a record means and two sources of loci. It exists for pileups written before
this file did, for a records file lost or built at knobs since changed, and above all for **growing** a
census, which is the one direction that cannot be served by subsetting an existing file
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3.4).

**What the file names, so a stale one cannot be used.** Beside the thirteen identity values of §1.5 it
carries the identity of the pileup it was built from: a digest of that pileup's header — reference,
analysed regions, read filters, command line — and its record count. **Not modification time.** On a
mismatch the caller rebuilds when the pileup is reachable and refuses naming the field when it is not.

**CLI surface.** The pileup run writes the records file as a side effect; a subcommand builds one from
an existing pileup. Nothing else needs a knob: the census size travels in the identity, and a run
wanting fewer loci subsets what is there rather than asking for a different file.

**NOT BUILT — and what is built instead.** The encoding that carries the *content* is in the types:
`PackedDepthCodes` is the five-bit array and its round trip is asserted at every bit offset, the
sparse lists are plain vectors, and `CoverageByWindow` stores the byte it means to store. What is
missing is only the framing that puts them in a file, plus a wire form for `SelectionIdentity` —
which holds `StrRepeatCriteria` and `ScanParams`, so its codec reaches into two other modules'
private fields. **The likely shape when it is written is a per-field digest table rather than the
values**: the identity's only use is equality with a name attached, and a table of
`(field name, digest)` gives exactly that without the codec having to track every field either type
grows.

### 2.3 Errors

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RecordError {
    /// The stream is not a record set, or is a version this build does not read.
    #[error("record stream is malformed or of an unknown version")]
    Malformed,
    #[error("i/o while reading or writing records")]
    Io(#[from] std::io::Error),
}
```

**The guard threshold is not an error here.** A locus over one non-whole-repeat read in ten is
well-formed data the *parameters fit* should decline to fit (spec §3.3); encoding it as a write failure would
make a property of the sample look like a property of the file.

---

## 3. Design decisions — decided

- **`DepthCode` wraps `DepthBin` and adds one sentinel; the ladder stays in `generic/`.** Five bits
  for twenty bins plus *never walked* — spec §2.2. A second ladder here would silently uncouple the
  two routes' binning.
- **The difference list carries `read`.** Without it an interruption cannot be told from two errors,
  and no other assertion in the suite notices — spec §3, §6.3.
- **`covering_not_crossing` is a field, not an inference.** The state has no field today and the
  censoring runs along repeat count — spec §3.
- **`ReadCap` is a newtype, not a `usize`.** It travels in `RecordIdentity` and is compared for
  equality; a bare integer beside the other counts there is transposable. **`DepthCap` is a second
  newtype for the same reason** and is not the same number.
- **Five values say in what *units* the evidence was recorded, not only which loci** — §1.5. The
  ladder digest is the one that matters most, because a depth code without its ladder is a number.
- **The coverage-by-window summary is a third object, not a field of a record** — §1.6; spec §4.
  Per sample rather than per read group, and it accumulates across the cohort like the records do,
  which is why spec §6 sizes it in the same table. **It is never serialized** (2026-08-13): it is
  rebuildable from the pileup, so it is built by the genome walk in a direct run and by a pass over the
  pileups otherwise — resident cost, not stored cost.
- **One `SampleRecords` for both kinds of run, differing only in whether its sections are resident** —
  §1.1. A second type for the file case would make the parameters fit two code paths for one object.
- **Sections are borrowed for the length of a call, never returned** — §2.2. The estimator holds the
  generic half and the repeat-tract half at different times and one band of strata within the second;
  an accessor that handed a section back would make that a convention a caller could ignore, and a
  file-backed value would grow into the file it exists to avoid.
- **The scoped call is on the cohort, and its unit is one stratum across every sample** — §2.2. A
  tract's length frequencies are fitted from every sample with reads there, so per-sample access would
  be the wrong grain; and it takes a **band** of strata, because 68 of tomato's 141 are fitted by
  borrowing from their neighbours.
- **`SsrRecords` is keyed by stratum** — §2.2; spec §6.2. The generic half and the repeat-tract half
  are never resident together, because the second consumes one number per sample from the first and
  returns nothing.
- **`CensusRecord` is an iteration item, not a stored one** — §1.1, §2.2.1. Stored, its enum width
  would put about 20 bytes on each of two million ordinary positions where five bits is what they need.
- **The records are a cache of the pileup, and the pileup is the source of truth** — spec §6.1.
  Everything in a records file can be recomputed from the sample's pileup; it is kept because
  recomputing means a full decompression pass, and because the file serves every future cohort call
  rather than one.
- **OPEN:** the recorded offset range (±4) and the span the parameters fit may place allele mass on (±6) are
  **two constants, not one** — spec §3.2. `OffsetCounts`'s width must not be read as the parameters fit's allele
  span.

---

## 4. Reconciliation with existing code

| this doc | existing code | how they meet |
|---|---|---|
| `DepthBin`, the ladder behind `DepthCode` | [`generic/depth_bins.rs:106,141`](../../../../src/ng/parameter_estimation/generic/depth_bins.rs) | reused unchanged; **do not mint a second ladder** |
| `RecordWriter::add_locus` | `GenericAccumulators::add_locus` ([`generic/accumulators.rs:278`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)) | same signature and same borrow, so one genome walk feeds both routes |
| `RecordWriter::merge` | `GenericAccumulators::merge` ([`generic/accumulators.rs:392`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)) | same contract: shard, then fold |
| what the writer is fed | `SampleLocusObservations` ([`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | taken as-is: `region`, `observations`, and the no-observation scalar |
| `ReadGroupId` | [`src/ng/types.rs:199`](../../../../src/ng/types.rs) | used as-is |
| `StratumCounts` | [`repeat_catalog/strata.rs:15`](../../../../src/ng/repeat_catalog/strata.rs) | stored, not restated |
| `SelectionIdentity`, `KeptLociDigest` | [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §1.2–1.3 | held, not redefined |

---

## 5. Open items

- **Impl-time:** `PackedDepthCodes`'s representation (5-bit packing over `Vec<u8>`, or `bitvec`). The
  contract is five bits per entry, index-addressable.
- ~~**Impl-time:** where `AlleleCount` stops being exact and starts binning.~~ **CLOSED 2026-08-13:
  it never starts.** The count is exact (spec §2.2). Almost every entry in the sparse list is one,
  two or three reads — a miscall — so the tail a ladder would compress is empty at any depth, and
  binning would buy a fourteenth bin width in `RecordIdentity` for nothing. `AlleleCount` was a name
  in this document that the code never needed; the field is a plain count.
- **Impl-time:** the wire format of `write_records`. **Not Parquet by reflex** — the catalog's reasons
  for it were columnar range queries, which nothing here does. It must be readable **sequentially in
  genome order across samples at once**, since that is how contamination reads a large census
  ([`parameter_prepass_joint_fit.md`](parameter_prepass_joint_fit.md) §3.4.4) — one locus's evidence
  from every sample, never a whole sample resident.
- **Impl-time:** how a records file names the pileup it came from (§2.2). The contract is a header
  digest plus a record count, compared for equality; which digest is an implementation choice.

---

## 6. Test shape

Tests live in `joint/records.rs` and need no alignment file: fill, write, read back, compare. Four
assertions carry the weight, each with a mutation that must fail it — the round-trip over every corner
state; the difference list separating a flank substitution from an interior one **and** two reads from
one read twice; the four STR states including the one with no field; and read groups folding exactly
on a two-group sample. **Sizes are measured rather than asserted**, on HG002 at 300× and on the whole
tomato cohort, reported separately (spec §7.8).

**One test does need an alignment file, and it is the one that validates the whole two-phase design**:
build a sample's records during a genome walk and again from the pileup that genome walk wrote, and compare the
files byte for byte (spec §7.12). It is what shows the pileup really holds everything a record needs,
and it fails on precisely the fields that do not survive the round trip.
