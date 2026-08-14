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

`src/ng/parameter_estimation/joint/census.rs`, beside `loci.rs` and `fit.rs`. **No trait**: one
record shape per path, no competing candidate, so a codec trait would be ceremony over two structs
and a writer.

---

## 1. Types

### 1.1 `SampleCensusEvidence` — one sample's evidence, however it is being held

**One type for both kinds of run.** A run that goes straight from alignments to a fit holds every
section in memory; a run reading a file fills a section when it is asked for and drops it after. That
is the *only* difference between them, so it is a property of the value rather than a second type —
and the parameters fit is written once, against calls that behave the same either way.

```rust
pub struct SampleCensusEvidence {
    pub sample: SampleName,
    /// The twelve values the parameters fit refuses to pool across (spec §5), and the
    /// digest of the loci actually kept. **Outside the sections**: they are compared
    /// across every sample before anything large is decoded.
    pub terms: RecordingTerms,
    /// Resident, or backed by a file and read on demand — §1.1a. Not public: what a
    /// caller may do with a section is §2.2's scoped access, and a field would let it
    /// keep one.
    sections: Sections,
}
```

**Contract.** A sample's counts at a position are the sum of its read groups' — exact, so the
parameters fit may fold freely. Sections are enumerated in a fixed order (read groups by id, strata by
their stratum key), so a fit that iterates is deterministic.

**What was here before, and why it is gone.** An earlier version of this section held the two halves as
**public maps keyed by read group**, and beside them a **second type** for the run that could not hold
them all. Two shapes for one object made the fit two code paths, and public fields made *"hold one
section at a time"* a convention rather than a property — nothing stopped a caller keeping every
stratum it had ever asked for. **The
storage stays exactly as §1.2 and §1.4 describe it**; what changed is that no caller reaches it
directly.

*Rejected, and it is worth saying why because it is the natural first idea:* one flat collection of
per-locus records, `Vec<LocusEvidence>` with an enum for the two kinds. **A Rust enum is as wide as its
widest variant**, and the repeat-tract variant carries nine two-byte offset buckets before its other
fields — so each of the two million ordinary positions would pay about 20 bytes where five bits is what
it needs, some 40 MB a read group against the packed array's 1.25 MB. Boxing the variants trades the
width for a pointer and an allocation per locus. And one collection in genome order interleaves the two
halves, so reading one stratum would mean walking all of it. **`LocusEvidence` survives one level up**,
as the item the genome walk emits and the fit iterates (§2.2.1), where the uniformity is worth having and
the storage is not paying for it.

### 1.1a `SectionKey` and `Sections` — what a section is, and where it is

**A section is the smallest piece of a sample's evidence anything ever asks for**, and there are
exactly two kinds: one read group's ordinary-position records, and one read group's tracts for one
stratum. That division is not an encoding detail — it is the estimator's own shape (spec §6.2), which
finishes the ordinary positions before reading a tract and reads one band of strata at a time.

```rust
/// Which section. **One name for three things that must not drift**: the entries of a
/// file's directory, the unit a call in §2.2 lends, and the unit spec §7.16 counts bytes
/// against.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum SectionKey {
    /// One read group's ordinary-position records (§1.2).
    Generic(ReadGroupId),
    /// One read group's tracts for one stratum (§1.4).
    Ssr(ReadGroupId, Stratum),
}

/// The decoded contents of one section.
pub enum Section {
    Generic(GenericEvidence),
    Ssr(SsrEvidence),
}

/// Where a sample's sections are. **Two states and no more**: a genome walk produces the
/// first, opening a file produces the second, and §1.1 exists so that no caller can tell
/// which it has.
enum Sections {
    /// Built by the walk. Every section decoded and held — this is the run that never
    /// writes a file, at fifty samples or at one.
    Resident(BTreeMap<SectionKey, Section>),
    /// Opened and checked, the directory read, **nothing else decoded**. Holds a reader
    /// and where each section sits in it, and no section at all between calls.
    Backed {
        reader: Box<dyn ReadSeek>,
        directory: BTreeMap<SectionKey, ByteExtent>,
    },
}

/// Where a section sits in a file: an offset and a length.
pub struct ByteExtent { offset: u64, len: u64 }
```

**Contract, and it is what makes §2.2's scoped access work for both states.** A scoped call decodes
the sections it was asked for **into a local**, hands the closure borrows of them, and drops them when
the closure returns; in the `Resident` state it borrows what is already there and drops nothing. The
closure sees `&GenericEvidence` or `&SsrEvidence` either way, so the parameters fit is one code path.
**`Sections` therefore retains nothing between calls in the `Backed` state** — not as a policy the
implementation follows, but because there is no field it could put a decoded section in.

**`SectionKey`'s ordering is the enumeration order** §1.1's contract promises: read groups by id, then
strata by their stratum key, with the ordinary-position sections before the tracts. A fit that
iterates sections is deterministic because this `Ord` is, and not because anything sorts.

**One read per section, then decode from a slice — and this is a requirement, not a preference.**
`Box<dyn ReadSeek>` is a deliberate choice: a type parameter for the reader would propagate through
`SampleCensusEvidence`, `CohortCensusEvidence` and `fit_jointly` to describe something only the
`Backed` state has, and would force every sample in a cohort to be backed the same way. Its cost is
one indirect call wherever the reader is touched, which is nothing **at the granularity this design
uses it**: the directory gives an offset and a length, so a call seeks once, fills a buffer once, and
decodes from the slice. A 63-sample cohort with 141 strata is about nine thousand section reads a
pass, moving 50 kB to 1.25 MB each.

**An implementation that read fields through the reader instead would pay that indirection millions of
times and lose inlining across the byte-level loop** — and if the reader were an unbuffered `File`,
each of those would be a syscall, which costs a thousand times what the indirection does. Decoding
from a slice also leaves the door open to memory-mapping the file, where a section is a subslice and
there is no read at all.

### 1.2 `GenericEvidence` — a dense array and a sparse list

Two parts, because at three reads a site nearly every position is *"n reads, all matching the
reference"* and only a few thousand in a million are anything else. The dense array holds the depth of
every kept position in order; the sparse list holds the exceptions, each naming the position it
belongs to by its **index into the selection** rather than by a coordinate.

```rust
pub struct GenericEvidence {
    /// Entry `i` is the `i`-th kept position's depth code, five bits, packed.
    depth: PackedDepthCodes,
    /// Only where a read was not on the reference base. **Sorted by `index`**, so the
    /// parameters fit walks both halves in one pass.
    non_reference: Vec<AlleleObservation>,
}

#[derive(Copy, Clone)]
pub struct AlleleObservation {
    /// Index into `CensusLoci::generic` — the position's only identity.
    pub index: u32,
    pub allele: ObservedAllele,
    /// **Exact, not binned, and one byte** (spec §2.2). Nearly every entry is a miscall of
    /// one to three reads, so a ladder's compressed tail stays empty at any depth and
    /// binning would cost a fourteenth bin width in `RecordingTerms` for nothing. The
    /// width is safe because the genome walk thins a position's counts by the same ratio
    /// it thins the depth, so `DepthCap` bounds this exactly — and `DepthCap` refuses a
    /// value a byte cannot hold.
    pub reads: u8,
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

### 1.4 `SsrEvidence` — lengths, a guard, and the differences

At a repeat tract the observation is a tract *length*, so the allele record is a distribution over
offsets from the reference tract length. Three things travel beside it, and each answers a question
the offsets cannot: how many reads reached the locus without crossing it (a lower bound on the
length), how many differed by something that is not a whole number of copies (the guard), and which
bases mismatched and where (the substitution channel — the error rate cannot be recovered from
lengths).

```rust
pub struct SsrEvidence {
    /// Spanning reads at each whole-repeat offset from the reference tract length,
    /// ±4 with saturating ends.
    offsets: Vec<OffsetCounts>,
    /// Reads that reached a tract in this stratum and crossed none of it — **one count for
    /// the whole stratum, not one per locus** (spec §3). The estimator summed the per-locus
    /// version the moment it read it, and the loss runs along repeat count, which is what a
    /// stratum is.
    covering_not_crossing: u32,
    /// Whether the genome walk reached this locus at all. **The STR half's answer to the
    /// generic half's never-walked sentinel**, and it has to be its own field: the
    /// other four vectors are all zero both for a locus never walked and for one
    /// walked with no read, and only the first is a bug (spec §6).
    /// A bit per locus, not a byte: `Vec<bool>` spent 0.46 MB a read group on tomato to
    /// carry 58 kB of information.
    walked: BitVec,
    /// Reads differing by a non-whole number of copies, with what they were.
    guard: Vec<GuardObservation>,
    /// The denominator the STR substitution rate is fitted against — **per stratum**, which
    /// is the grain that rate is fitted at (spec §3). Per locus it was 1.85 MB a read group
    /// on tomato that nothing read.
    bases_compared: u64,
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

**Contract.** All the per-locus vectors are indexed by the locus's position in `CensusLoci`, as the
generic half is indexed by position. An offset beyond ±4 saturates into the end bucket rather than
being dropped or wrapping.

**`SsrLocusState` is what a consumer asks**, rather than reading the vectors and re-deriving the four
states — never walked, walked with no read, reached but not crossed, crossed — each time. The guard's
one-in-ten threshold is a query on the same type (`guard_is_over_threshold`) and **not** a write-time
error: a locus over it is well-formed data the parameters fit should decline to fit, and refusing the write would
make a property of the sample look like a property of the file.

### 1.5 `RecordingTerms`

**The terms under which one sample's evidence was recorded**, travelling with it so the parameters fit
can refuse to pool two runs that did not record the same thing. Seven say which loci were *asked for*,
one says which came back, and **four say in what units the evidence was written down** — the group an
earlier version of this section did not have, which left the whole check able to pass while two
samples' rows meant different things (spec §5).

*Renamed 2026-08-13 from `RecordIdentity`, and its sibling in `loci.rs` from `SelectionIdentity` to
`SelectionTerms`.* The pair now says what each one is the terms **of**: `SelectionTerms` is how the
loci were chosen, `RecordingTerms` is how the evidence was recorded at them, and the second contains
the first. *"Record"* named no particular record and *"identity"* named an abstraction; the terms are
the thing itself, and the sentence a mismatch produces reads in plain English — **these two samples
recorded under different terms.**

```rust
#[derive(Clone, PartialEq)]
pub struct RecordingTerms {
    /// The seven that say which loci were asked for.
    pub selection: SelectionTerms,
    /// The one that says which loci came back (`loci.rs` §1.3).
    pub kept_loci: CensusLociDigest,

    // ---- the four that say in what units ----
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
}
```

**`PartialEq` and not `Eq`**, forced by `SelectionTerms` holding float floors — compare those by
bit pattern (`loci.rs` §1.2).

**A thirteenth value was here and is gone — 2026-08-13.** It was the coverage-by-window grid width,
guarding against two samples' window summaries being built on different grids. **That summary is never
written to a file** (spec §4), so it is not in the artifact these terms describe, and whoever runs the
parameters fit builds every sample's itself, in one process, at one width — the disagreement had no way
to arise. A term describing an object the file does not contain is worse than no term: a reader takes
each one as guarding something in the file. If summaries are ever stored, the width comes back
**attached to the summary**, not to evidence that does not hold one.

### 1.6 `CoverageByWindow` — REMOVED 2026-08-14

**The type is gone and so is the object.** It held one sample's mean read depth over fixed 500 bp
windows, the sample's depth-against-GC curve and each window's GC fraction, and it existed so that the
duplicated-copy class could be conditioned on the coverage around a position rather than on the
position itself. Spec §4 says why it went and what replaced it: the cohort's genotype composition above
about twenty-five samples, and **the position's own depth against the sample's median** above about
twenty-five reads a position — the second of which is already in `GenericEvidence` and needed only
§2.2's ladder extended to reach past 76 reads a position.

**Nothing else in this document referred to it**, which is the cleanest evidence that it was a
subsystem rather than a field: no other type held one, no interface passed one, and the only signature
that mentioned it was `CensusWriter::finish`.

**`coverage.rs` is now dead code** — `src/ng/parameter_estimation/joint/coverage.rs`, about 560 lines,
built during the walk and read by nothing. Delete it with this change rather than leaving it to be
found.

---

## 2. Interfaces

### 2.1 Filling records during the genome walk

The writer is handed the same locus stream the histogram accumulators are handed, so one genome walk fills
both routes and the comparison between them is over identical evidence. It knows which loci are kept
and ignores the rest.

**Built, 2026-08-12.**

```rust
pub struct CensusWriter { /* … */ }

impl CensusWriter {
    pub fn new(
        sample: String,
        loci: &CensusLoci,
        read_groups: Vec<ReadGroupId>,
        contig_of: &dyn Fn(&str) -> Option<ContigId>,
        terms: SelectionTerms,
        edges: DepthBinEdges,
        read_cap: ReadCap,
        depth_cap: DepthCap,
    ) -> Self;

    /// Record this locus if it is a kept one. **Borrows and does not take**, so the genome walk
    /// passes the locus on to the histogram accumulators untouched and one pass fills
    /// both routes.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations);

    pub fn finish(self) -> SampleCensusEvidence;
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
pub fn write_records(records: &SampleCensusEvidence, out: impl Write) -> Result<(), CensusError>;
pub fn read_records(input: impl Read) -> Result<SampleCensusEvidence, CensusError>;
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
pub struct CohortCensusEvidence { /* Vec<SampleCensusEvidence>, in sample-name order */ }

impl CohortCensusEvidence {
    /// Opens or adopts each sample's records and checks the twelve recording terms
    /// across all of them **before any section is decoded** (spec §5).
    pub fn new(samples: Vec<SampleCensusEvidence>) -> Result<Self, CensusError>;

    pub fn read_groups(&self) -> &[ReadGroupId];
    pub fn strata(&self) -> &[Stratum];

    /// Every sample's generic records for one read group, for the length of the call.
    pub fn with_generic<R>(
        &mut self,
        group: ReadGroupId,
        f: impl FnOnce(&[&GenericEvidence]) -> R,
    ) -> Result<R, CensusError>;

    /// Every sample's tracts for a **band** of strata, for the length of the call.
    pub fn with_strata<R>(
        &mut self,
        group: ReadGroupId,
        strata: &[Stratum],
        f: impl FnOnce(&[&SsrEvidence]) -> R,
    ) -> Result<R, CensusError>;
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
memory is samples × the largest band**, and **what bounds one section is the per-stratum cap** — the
memory guarantee the whole by-section shape exists for, so it belongs in this contract rather than
being left to the loci document. At the cap measured on 2026-08-13, **5,000 tracts**, a tract costs
about ten bytes a read group: **50 kB a sample, 50 MB across a thousand**
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §6 question 1). Uncapped, the
largest stratum tomato holds is 217,812 tracts, so the same section would be 2.2 MB a sample and
2.2 GB across a thousand — which is what the cap is for.

**`SsrEvidence` becomes per (read group × stratum) rather than per read group.** Today it holds vectors
indexed by a locus's position in the whole kept set; a section holds one stratum's slice of that, so
the index is stratum-local and the stratum's first index travels in the file's directory. **That is the
one storage-type change this decision forces**, and it is the reason it is recorded here rather than
left to implementation.

### 2.2.1 `LocusEvidence` — the item both paths iterate

The two halves store different things and are read through different calls, but a genome walk emits one
locus at a time and a fit consumes one locus at a time, and at *that* level they are the same shape.

```rust
/// One locus's evidence for one read group, as the walk emits it and the fit reads it.
pub enum LocusEvidence<'a> {
    Generic { index: u32, depth: DepthCode, non_reference: &'a [AlleleObservation] },
    Ssr { index: u32, stratum: Stratum, tract: &'a SsrLocus },
}
```

**Contract.** This is an **iteration item, never a stored one** — it borrows from a section that is
already decoded, so a two-million-entry sequence of them exists only as a cursor. §1.1 says what
storing them instead would cost.

**Two producers, one builder.** The genome walk-time producer is `CensusWriter` (§2.1). The second reads an
existing pileup and drives the same `CensusWriter` through the same locus stream, so there is one
implementation of what a record means and two sources of loci. It exists for pileups written before
this file did, for a records file lost or built at knobs since changed, and above all for **growing** a
census, which is the one direction that cannot be served by subsetting an existing file
([`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §4.3.4).

**The directory is `SectionKey` → `ByteExtent`** (§1.1a), so the file's layout and the calls above are
described by one type. A reader that cannot find the key it was asked for fails naming it, rather than
returning an empty section — a stratum with no tracts in this sample is still a section, holding zero.

**What the file names, so a stale one cannot be used.** Beside the twelve recording terms of §1.5 it
carries the identity of the pileup it was built from: a digest of that pileup's header — reference,
analysed regions, read filters, command line — and its record count. **Not modification time.** On a
mismatch the caller rebuilds when the pileup is reachable and refuses naming the field when it is not.

**CLI surface.** The pileup run writes the records file as a side effect; a subcommand builds one from
an existing pileup. Nothing else needs a knob: the census size travels in the terms, and a run
wanting fewer loci subsets what is there rather than asking for a different file.

**NOT BUILT — and what is built instead.** The encoding that carries the *content* is in the types:
`PackedDepthCodes` is the five-bit array and its round trip is asserted at every bit offset, the
sparse lists are plain vectors. What is
missing is only the framing that puts them in a file, plus a wire form for `SelectionTerms` —
which holds `StrRepeatCriteria` and `ScanParams`, so its codec reaches into two other modules'
private fields. **The likely shape when it is written is a per-field digest table rather than the
values**: these terms' only use is equality with a name attached, and a table of
`(field name, digest)` gives exactly that without the codec having to track every field either type
grows.

### 2.3 Errors

```rust
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CensusError {
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
- **`ReadCap` is a newtype, not a `usize`.** It travels in `RecordingTerms` and is compared for
  equality; a bare integer beside the other counts there is transposable. **`DepthCap` is a second
  newtype for the same reason** and is not the same number. **It also refuses a value above
  `u8::MAX`**, because it is the bound on `AlleleObservation::reads` (§1.2): the walk thins a
  position's counts by the same ratio it thins the depth, so a cap a byte cannot hold would saturate
  the counts silently while the depth field said otherwise (spec §2.2).
- **Five values say in what *units* the evidence was recorded, not only which loci** — §1.5. The
  ladder digest is the one that matters most, because a depth code without its ladder is a number.
- **There is no third object** — §1.6; spec §4. A per-sample coverage-by-window summary was specified
  until 2026-08-14 and is removed: the duplicated-copy class is conditioned on the cohort's genotype
  composition and on the position's own depth, both of which the records already carry.
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
| `CensusWriter::add_locus` | `GenericAccumulators::add_locus` ([`generic/accumulators.rs:278`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)) | same signature and same borrow, so one genome walk feeds both routes |
| `CensusWriter::merge` | `GenericAccumulators::merge` ([`generic/accumulators.rs:392`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)) | same contract: shard, then fold |
| what the writer is fed | `SampleLocusObservations` ([`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | taken as-is: `region`, `observations`, and the no-observation scalar |
| `ReadGroupId` | [`src/ng/types.rs:199`](../../../../src/ng/types.rs) | used as-is |
| `StratumCounts` | [`repeat_catalog/strata.rs:15`](../../../../src/ng/repeat_catalog/strata.rs) | stored, not restated |
| `SelectionTerms`, `CensusLociDigest` | [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §1.2–1.3 | held, not redefined |

**The code has not been renamed yet, and this document leads it — 2026-08-13.** The module is still
`records.rs` and the types still carry their old names. The sweep is mechanical:

| this document | the code today |
|---|---|
| `census.rs` | `records.rs` |
| `CensusLoci`, `CensusLociDigest` | `KeptLoci`, `KeptLociDigest` |
| `SampleCensusEvidence` | `SampleRecords` |
| `CohortCensusEvidence` | *(new — no code yet)* |
| `GenericEvidence`, `SsrEvidence` | `GenericRecords`, `SsrRecords` |
| `CensusWriter` | `RecordWriter` |
| `LocusEvidence` | *(new — no code yet)* |
| `CensusError` | `RecordError` |
| `RecordingTerms`, `SelectionTerms` | `RecordIdentity`, `SelectionIdentity` |

**It waits rather than colliding.** The sweep touches `records.rs`, `loci.rs`, `fit.rs`,
`contamination.rs` and five examples, and the last two are being changed elsewhere as this is written.
Do it in one pass when that work is committed — a partial rename leaves the module saying both, which
is the state this rename exists to end.

*Why "census" at all.* It is the word for what this route does that "records" does not carry: **the
same questions asked of every sample**. A record is any written-down thing; a census is a set of
questions put identically to a whole population, which is the one property the entire route rests on.

---

## 5. Open items

- **Impl-time:** `PackedDepthCodes`'s representation (5-bit packing over `Vec<u8>`, or `bitvec`). The
  contract is five bits per entry, index-addressable.
- ~~**Impl-time:** where `AlleleCount` stops being exact and starts binning.~~ **CLOSED 2026-08-13:
  it never starts.** The count is exact (spec §2.2). Almost every entry in the sparse list is one,
  two or three reads — a miscall — so the tail a ladder would compress is empty at any depth, and
  binning would buy a fourteenth bin width in `RecordingTerms` for nothing. `AlleleCount` was a name
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

Tests live in `joint/census.rs` and need no alignment file: fill, write, read back, compare. Four
assertions carry the weight, each with a mutation that must fail it — the round-trip over every corner
state; the difference list separating a flank substitution from an interior one **and** two reads from
one read twice; the four STR states including the one with no field; and read groups folding exactly
on a two-group sample. **Sizes are measured rather than asserted**, on HG002 at 300× and on the whole
tomato cohort, reported separately (spec §7.8).

**One test does need an alignment file, and it is the one that validates the whole two-phase design**:
build a sample's records during a genome walk and again from the pileup that genome walk wrote, and compare the
files byte for byte (spec §7.12). It is what shows the pileup really holds everything a record needs,
and it fails on precisely the fields that do not survive the round trip.
