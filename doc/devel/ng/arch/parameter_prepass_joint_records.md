# ng — the joint fit, what is recorded at each kept locus: types & interfaces

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

### 1.1 `SampleRecords` — the whole input to the fit

One sample's evidence at the kept loci, plus the values the fit checks before pooling. The two paths
are separate maps because they hold different things, and both are keyed by read group because that
is the grain the walk fills them at.

```rust
pub struct SampleRecords {
    pub sample: SampleName,
    /// Indexed by position in `KeptLoci::generic`; carries no coordinates (spec §2.3).
    pub generic: BTreeMap<ReadGroupId, GenericRecords>,
    pub ssr: BTreeMap<ReadGroupId, SsrRecords>,
    /// The ten values the fit refuses to pool across (spec §4).
    pub identity: RecordIdentity,
}
```

**Contract.** A sample's counts at a position are the sum of its read groups' — exact, so the fit may
fold freely. `BTreeMap` and not `HashMap`, so a fit that iterates is deterministic.

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
    /// fit walks both halves in one pass.
    non_reference: Vec<AlleleObservation>,
}

#[derive(Copy, Clone)]
pub struct AlleleObservation {
    /// Index into `KeptLoci::generic` — the position's only identity.
    pub index: u32,
    pub allele: ObservedAllele,
    pub reads: AlleleCount,
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
    /// two on two reads** — a read-blind encoding passes every other check (spec §6.3).
    pub read: u8,
    /// Signed offset from the tract start: negative in the left flank, `0..len` inside,
    /// beyond `len` in the right flank.
    pub offset: i16,
    pub base: ObservedAllele,
}
```

**Contract.** All five vectors are indexed by the locus's position in `KeptLoci`, as the generic half
is indexed by position. An offset beyond ±4 saturates into the end bucket rather than being dropped
or wrapping.

### 1.5 `RecordIdentity` — the ten values

What travels with a sample so the fit can refuse to pool two runs that did not record the same thing.
Nine say what was *asked for* and one — the digest — says what came back; the last two are this
module's own, being properties of the recording rather than of the selection.

```rust
#[derive(Clone, PartialEq, Eq)]
pub struct RecordIdentity {
    /// The seven that say which loci were asked for.
    pub selection: SelectionIdentity,
    /// The one that says which loci came back (`loci.rs` §1.3).
    pub kept_loci: KeptLociDigest,
    /// Per stratum, held against kept: anything pooled across strata is biased without
    /// it and silently so (spec §4).
    pub ssr_stratum_counts: StratumCounts,
    /// Must match the other route's, or the comparison confounds the cap with the
    /// route (spec §3.4).
    pub read_cap: ReadCap,
}
```

---

## 2. Interfaces

### 2.1 Filling records during the walk

The writer is handed the same locus stream the histogram accumulators are handed, so one walk fills
both routes and the comparison between them is over identical evidence. It knows which loci are kept
and ignores the rest.

```rust
pub struct RecordWriter { /* … */ }

impl RecordWriter {
    pub fn new(loci: &KeptLoci, identity: SelectionIdentity, read_cap: ReadCap) -> Self;

    /// Record this locus if it is a kept one. **Borrows and does not take**, so the walk
    /// passes the locus on to the histogram accumulators untouched and one pass fills
    /// both routes.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations);

    /// Fold a shard's writer in. Concatenation in position order.
    pub fn merge(&mut self, other: Self);

    pub fn finish(self) -> SampleRecords;
}
```

**Contract.** Every kept locus gets an entry whether or not a read reached it — the entry is the
denominator. A locus in a region never visited keeps `DepthCode::NeverWalked`, so the three states
survive a write and a read. **`add_locus` is the only place the kept-loci digest is fed**, so a
record set and its digest cannot disagree.

### 2.2 Getting records to the fit

**Where they live is a non-goal** (spec §1.2); the requirement is that the fit reaches every sample's
without walking reads again. So the type is serializable and says nothing about where it is put.

```rust
pub fn write_records(records: &SampleRecords, out: impl Write) -> Result<(), RecordError>;
pub fn read_records(input: impl Read) -> Result<SampleRecords, RecordError>;
```

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
well-formed data the *fit* should decline to fit (spec §3.3); encoding it as a write failure would
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
  equality; a bare integer beside the other counts there is transposable.
- **OPEN:** the recorded offset range (±4) and the span the fit may place allele mass on (±6) are
  **two constants, not one** — spec §3.2. `OffsetCounts`'s width must not be read as the fit's allele
  span.

---

## 4. Reconciliation with existing code

| this doc | existing code | how they meet |
|---|---|---|
| `DepthBin`, the ladder behind `DepthCode` | [`generic/depth_bins.rs:106,141`](../../../../src/ng/parameter_estimation/generic/depth_bins.rs) | reused unchanged; **do not mint a second ladder** |
| `RecordWriter::add_locus` | `GenericAccumulators::add_locus` ([`generic/accumulators.rs:278`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)) | same signature and same borrow, so one walk feeds both routes |
| `RecordWriter::merge` | `GenericAccumulators::merge` ([`generic/accumulators.rs:392`](../../../../src/ng/parameter_estimation/generic/accumulators.rs)) | same contract: shard, then fold |
| what the writer is fed | `SampleLocusObservations` ([`locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)) | taken as-is: `region`, `observations`, and the no-observation scalar |
| `ReadGroupId` | [`src/ng/types.rs:199`](../../../../src/ng/types.rs) | used as-is |
| `StratumCounts` | [`repeat_catalog/strata.rs:15`](../../../../src/ng/repeat_catalog/strata.rs) | stored, not restated |
| `SelectionIdentity`, `KeptLociDigest` | [`parameter_prepass_joint_loci.md`](parameter_prepass_joint_loci.md) §1.2–1.3 | held, not redefined |

---

## 5. Open items

- **Impl-time:** `PackedDepthCodes`'s representation (5-bit packing over `Vec<u8>`, or `bitvec`). The
  contract is five bits per entry, index-addressable.
- **Impl-time:** where `AlleleCount` stops being exact and starts binning — a `pub const` with units
  and source in its doc comment (spec §2.2).
- **Impl-time:** the wire format of `write_records`. **Not Parquet by reflex** — the catalog's reasons
  for it were columnar range queries, which nothing here does.

---

## 6. Test shape

Tests live in `joint/records.rs` and need no alignment file: fill, write, read back, compare. Four
assertions carry the weight, each with a mutation that must fail it — the round-trip over every corner
state; the difference list separating a flank substitution from an interior one **and** two reads from
one read twice; the four STR states including the one with no field; and read groups folding exactly
on a two-group sample. **Sizes are measured rather than asserted**, on HG002 at 300× and on the whole
tomato cohort, reported separately (spec §6.8).
