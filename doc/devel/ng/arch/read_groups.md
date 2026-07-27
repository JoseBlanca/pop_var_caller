# ng — read groups: types & interfaces

*Status: architecture draft (2026-07-27), companion to the spec
[`../spec/read_groups.md`](../spec/read_groups.md) (the design and every "why"). Stands on
[`alignment_file.md`](alignment_file.md) and [`sample_reads.md`](sample_reads.md), and **changes
both** — the removals are in §5. Shared vocabulary: [`ng_step_interfaces.md`](ng_step_interfaces.md)
§1; module rules: [`module_layout.md`](module_layout.md). Naming per
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md). Signatures are
illustrative; the **contract** is the deliverable.*

## Module home

```
src/ng/types.rs           – ReadGroupId (shared vocabulary, §1.1)
src/ng/read/
├── aligned_read.rs       – NEW: AlignedRead, the ng-owned decoded read (§1.4)
├── filtering.rs          – RecordSource/RawRecord: carry and resolve the read group (§3.3)
└── input/
    ├── read_groups.rs    – NEW: ReadGroup, ReadGroups, SampleReadGroups, NameWithOrigin,
    │                       ReadGroupResolution, build_read_groups, ReadGroupError (§1–§3)
    ├── open_bam.rs       – AlignmentFile: loses sample_name, gains its resolution
    ├── region_query.rs   – the record sources: carry the resolution per query
    ├── mod.rs            – SampleReads: opened from read groups, not paths
    └── merge.rs          – unchanged
```

`read_groups.rs` sits in `input/` because the only `@RG` parser in ng is already there
([`open_bam.rs:856`](../../../../src/ng/read/input/open_bam.rs#L856)) and because the plan it
produces is what `SampleReads` is built from. The ng read type is one level up, in `read/`, because
filtering, preparation and locus generation all speak it.

## 1. Types

### 1.1 `ReadGroupId` — shared vocabulary, not this module's

"Which read group a read came from" is spoken by reads, by observations, and later by the parameter
pre-pass. Like `GenomePosition` before it ([`sample_reads.md`](sample_reads.md) §1.1) it goes in
`ng::types` and is added to `ng_step_interfaces.md` §1; this module is its first consumer, not its
owner. Unconstrained, so a public field and `get()`, matching `ContigId`
([`types.rs:13`](../../../../src/ng/types.rs#L13)).

```rust
/// Which read group a read belongs to: an index into the run's [`ReadGroups`] table.
/// Unconstrained — any `u32` is a legal index at the type level, and an out-of-range id
/// is caught at lookup — so the field is public and there is no checked constructor.
///
/// Minted only by `build_read_groups`, in input-file order then header order, so the same input
/// list always yields the same ids whatever order the files are opened in (spec §4).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct ReadGroupId(pub u32);
```

### 1.2 The read group and the table

One record per `@RG`. The atoms (`file`, `id`, `sample`) are read from the input; the two grouping
names may have been synthesized, and say so (spec §3, §6).

```rust
/// One `@RG` record: the file that declared it, what it says, and the grouping names —
/// declared or synthesized.
pub struct ReadGroup {
    /// `Arc` because a file's k read groups share one path, and because
    /// `AlignmentFile` already holds its path that way (`open_bam.rs:75`).
    pub file: Arc<Path>,
    /// `@RG ID`, verbatim. Unique within its file only — never an identity (spec §4).
    pub id: Box<str>,
    /// `@RG SM`. Always present: absence is a hard error at the pre-pass.
    pub sample: Box<str>,
    /// `@RG LB`, or synthesized from sample + id + file stem.
    pub library: NameWithOrigin,
    /// The sequencing experiment; falls back to the library.
    pub experiment: NameWithOrigin,
    /// `@RG PL`. Reports only — nothing keys on it (spec §6).
    pub platform: Option<Box<str>>,
}

/// A name used for grouping, plus where it came from. The origin cannot be recovered
/// later, and a chemistry-group report has to be able to say "we made this one up".
pub struct NameWithOrigin {
    pub value: Box<str>,
    pub origin: NameOrigin,
}

pub enum NameOrigin { Declared, Synthesized }

/// Every read group in the run, in two views of the **same** set: by identifier, and
/// grouped by the sample each names. Built once by `build_read_groups`, then read-only
/// and shared.
pub struct ReadGroups { /* Vec<ReadGroup> + Vec<SampleReadGroups> */ }

impl ReadGroups {
    pub fn get(&self, id: ReadGroupId) -> &ReadGroup;
    pub fn len(&self) -> usize;
    pub fn iter(&self) -> impl Iterator<Item = (ReadGroupId, &ReadGroup)>;
    /// One entry per sample, in first-seen order — what each open is built from.
    pub fn read_groups_per_sample(&self) -> &[SampleReadGroups];
}

/// One sample and the read groups that name it. The unit `SampleReads::open` takes.
pub struct SampleReadGroups {
    pub sample: Box<str>,
    pub read_groups: Vec<ReadGroupId>,
}
```

The by-sample view is not a second collection — it is the same read groups grouped by
`ReadGroup::sample`, so it belongs on this type rather than in a wrapper holding both.

### 1.3 How a file's records are read

Two questions have to be answered about every record that comes out of a file:

1. **Which read group does it belong to?**
2. **Does it belong to the sample this open is for?** A file may declare read groups for several
   samples, and an open always serves exactly one (spec §9).

`ReadGroupResolution` answers both, and it is built once — when the file is opened — then consulted
for every record of that file.

**Why once, and not per record.** *How* a file's records must be read depends only on how many
`@RG` its header declares, which cannot change while the file is open (spec §7). A file declaring
one read group needs no per-record work at all.

**Who builds it: `SampleReads`, once per open.** Question 2 is about the sample being opened, not
about the file, so the same file opened for two samples gets two different resolutions — each
marking a different subset of the file's read groups as its own.

It is an enum rather than an `Option<ReadGroupId>`; the reason is in §4.

```rust
/// How this open reads this file's records — decided once, when the file is opened.
pub enum ReadGroupResolution {
    /// The header declares exactly one read group; every record is that one and the
    /// record's `RG` is not read. Such a file is single-sample by construction.
    Sole(ReadGroupId),
    /// The header declares several; each record's `RG` names which. A record with no
    /// `RG`, or naming none of these, is fatal.
    PerRecord(Box<[(Box<str>, RecordOwner)]>),
}

/// What a declared read group means to *this* open.
pub enum RecordOwner {
    /// This open's sample. The read is yielded, carrying this id.
    Mine(ReadGroupId),
    /// Declared in the file but naming another sample: the read is skipped, and tallied
    /// apart from the drop categories — it is not a quality drop (spec §9).
    OtherSample,
}
```

`PerRecord` is a small array, not a `HashMap`: k is the number of `@RG` in one file. Linear or
sorted is a code-review question, not this doc's.

### 1.4 The ng read

ng owns its decoded read (spec §8), so the identifier can ride on it and production stays unaware
of ng. Same fields as production's `MappedRead`
([`alignment_input.rs:78`](../../../../src/bam/alignment_input.rs#L78)) with `source_file_index`
replaced — the read group knows its file.

```rust
/// One decoded, filtered read. ng's own, modelled on `bam::alignment_input::MappedRead`;
/// the difference is the last field.
pub struct AlignedRead {
    pub qname: Vec<u8>,
    pub flag: u16,
    pub ref_id: usize,
    pub pos: u64,
    pub mapq: u8,
    pub cigar: Vec<CigarOp>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub mate_ref_id: Option<usize>,
    pub mate_pos: Option<u64>,
    pub adaptor_boundary: Option<u32>,
    /// Which read group this read came from. Resolved at decode (§3.3).
    pub read_group: ReadGroupId,
}
```

**Not** the place to revisit the per-read allocations (`qname`, `seq`, `qual`, `cigar`) — that is
parked as measure-first in [`../impl_plan/read_filtering.md`](../impl_plan/read_filtering.md) and
this change only unblocks it.

## 2. Errors

Two surfaces, because they fire at two times. Header failures are the pre-pass's and happen before
any read; a record failure travels in the read stream.

```rust
/// Header-level failures, all raised by `build_read_groups` before any read flows (spec §6).
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReadGroupError {
    /// The file declares no `@RG` at all. Message names the file and says it must be
    /// re-headered.
    NoReadGroups { path: PathBuf },
    /// An `@RG` carries no `SM`. Distinct from the above: different fix, so different
    /// variant (spec §6).
    MissingSampleName { path: PathBuf, read_group_id: String },
    /// Two read groups synthesized the same library name — two input files with the
    /// same stem in different directories. Names both full paths (spec §6).
    DuplicateSynthesizedLibrary { library: String, paths: (PathBuf, PathBuf) },
}

/// What a record source can fail with. Replaces the bare `io::Error` in
/// `ReadFilterError::Source`, so an unresolvable read group can name the read.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum RecordSourceError {
    Io(#[source] io::Error),
    /// A record with no `RG` in a file that declares several — nothing to assign it to.
    RecordWithoutReadGroup { path: Arc<Path>, qname: Vec<u8>, position: GenomePosition },
    /// A record naming a read group the header never declared.
    UndeclaredReadGroup { path: Arc<Path>, qname: Vec<u8>, named: Box<str> },
}
```

**Why a new source error rather than folding these into `io::Error`.** `ReadFilterError::Source`
wraps `io::Error` today ([`filtering.rs:586`](../../../../src/ng/read/filtering.rs#L586)); stuffing
a read-group failure in there would render it as a string and lose the read's identity, which is the
one thing the message needs (spec §7).

## 3. Interfaces

### 3.1 The pre-pass

```rust
/// Read every input file's header — headers only — parse its `@RG` records, apply the
/// §6 rules, assign ids, and group by sample.
///
/// **Contract.** Total: every returned `ReadGroup` has a sample, a library and an
/// experiment. Deterministic: ids follow `paths` order, then header order, whatever
/// order the files were read in. Fails on any of the three header errors; opens no
/// index and reads no record.
pub fn build_read_groups(paths: &[PathBuf]) -> Result<ReadGroups, ReadGroupError>;
```

### 3.2 What changes on `AlignmentFile` and `SampleReads`

```rust
impl AlignmentFile {
    /// Takes the resolution `SampleReads` built for this open instead of deriving a
    /// sample name from the header. `sample_name()` is gone (§5).
    pub fn open_as(
        path: &Path,
        reference: &ReferenceInfo,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
        resolution: ReadGroupResolution,
    ) -> Result<Self, AlignmentFileError>;

    /// Per read group, not per file: a drop rate is a read group's property, and a file
    /// may now hold several (spec §8).
    pub fn counts(&self) -> Vec<(ReadGroupId, ReadFilterCounts)>;
}

impl SampleReads {
    /// Open the files this sample's read groups live in, then check they all name one
    /// sample — the single-sample-per-open invariant (spec §4), enforced here even
    /// though a plan-built call never trips it.
    pub fn open(
        sample: &SampleReadGroups,
        read_groups: &ReadGroups,
        reference: &ReferenceInfo,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
    ) -> Result<Self, IngestError>;
}
```

**Contract, unchanged where it matters.** `reads_in_region` still yields one sample's reads in
coordinate order, lazily, fused. What is added: every read carries its `ReadGroupId`, and a read
belonging to another sample's read group in a shared file is not yielded and is **not** counted as a
drop — it gets its own tally (spec §9).

### 3.3 Resolving a record

The resolution rides on the reused record buffer, refreshed once per query alongside the header the
sources already borrow ([`region_query.rs:67`](../../../../src/ng/read/input/region_query.rs#L67)),
and is applied in `decode` — **not** in `read_next`. Reads dropped by the pre-decode gate (unmapped,
secondary, duplicate, low MAPQ) never reach `decode`, so resolving there costs nothing for them; the
`Sole` arm is a match and a copy either way.

```rust
pub trait RawRecord {
    /// Decode, resolving the read group. Fatal on failure, as before — the variants are
    /// `RecordSourceError`'s.
    fn decode(&self) -> Result<AlignedRead, RecordSourceError>;
}
```

## 4. Design decisions — decided

- **`ReadGroupId` lives in `ng::types`, not here — decided.** It is cross-step vocabulary, and
  naming it after one consumer invites the next to mint a duplicate (§1.1; the `GenomePosition`
  precedent, [`sample_reads.md`](sample_reads.md) §1.1).
- **The identifier is a run-wide index, not a `(file, @RG)` pair or a string — decided.** One number
  per read, one lookup; a per-open numbering would give one physical read group two ids when a file
  is opened for two samples. Why: **spec §4**.
- **`ReadGroupResolution` is an enum, not `Option<ReadGroupId>` — decided.** `None` would be a
  sentinel for "use the map" (§1.3); spec §7.
- **The other-sample filter lives inside the resolution, not beside it — decided** (§1.3). A
  separate "wanted ids" set passed alongside would be a second structure that has to agree with the
  first; folding it in means one lookup, and disagreement is unrepresentable.
- **Resolution happens at `decode`, not `read_next` — decided** (§3.3). The rejected order costs an
  auxiliary-tag lookup on every read the pre-decode gate is about to drop.
- **A new `RecordSourceError` replaces `io::Error` in `ReadFilterError::Source` — decided** (§2).
- **`AlignedRead` is ng's own type — decided.** The live alternative was adding a field to
  production's `MappedRead`, which is a small edit but makes production carry an ng concept; spec
  §8.
- **`counts()` is keyed by read group — decided** (§3.2); spec §8.
- **The by-sample view lives on `ReadGroups`, not in a wrapper holding both — decided** (§1.2). It
  is the same read groups grouped by `ReadGroup::sample`, so a type pairing the two would be one
  collection wearing two names.
- **`build_read_groups` is a free function, not a builder — decided.** It runs once and its output
  is immutable; there is no partial state worth naming.
- **No trait, no bake-off.** One implementation of everything here. This is a single module, not a
  swappable step.

## 5. Reconciliation with existing code

Every row read at the cited line. The first four are **removals** — code that exists and goes.

| ng name | existing code | action |
|---|---|---|
| — | `SampleNames` enum + `sample_names()` [`open_bam.rs:835`, `:856`](../../../../src/ng/read/input/open_bam.rs#L835) | **replace**: parse the whole `@RG` record, not just `SM` |
| — | `AlignmentFileError::MultipleSampleNames` [`mod.rs:132`](../../../../src/ng/read/input/mod.rs#L132) | **delete**: several samples per file is a normal input (spec §8) |
| `ReadGroupError::{NoReadGroups, MissingSampleName}` | `AlignmentFileError::MissingSampleName { read_group: Option<String> }` [`mod.rs:151`](../../../../src/ng/read/input/mod.rs#L151) — one variant for both states | **split** into two, moved to the pre-pass |
| — | `AlignmentFile::sample_name()` [`open_bam.rs:417`](../../../../src/ng/read/input/open_bam.rs#L417) | **delete**: the sample is the plan's, not a file's |
| `build_read_groups` | no counterpart; [`ng_ssr_cohort_stutter.rs:175`](../../../../examples/ng_ssr_cohort_stutter.rs#L175) groups paths by `SM` by hand, opening each file once to probe | **new**; that ad-hoc grouping becomes a call |
| the single-sample check | `agreed_sample_name` [`mod.rs:550`](../../../../src/ng/read/input/mod.rs#L550) | **keep**, re-keyed on read groups (§3.2) |
| `ReadGroupId` | no counterpart — reads carry `source_file_index: usize` [`alignment_input.rs:102`](../../../../src/bam/alignment_input.rs#L102) | **new**, in `ng::types` |
| `AlignedRead` | `MappedRead` [`alignment_input.rs:78`](../../../../src/bam/alignment_input.rs#L78) | **model, not reuse**; copy `record_buf_to_mapped_read` [`:803`](../../../../src/bam/alignment_input.rs#L803) |
| decode helpers | `compute_adaptor_boundary` is `pub(crate)` [`:883`](../../../../src/bam/alignment_input.rs#L883); `cigar_to_ops` is module-private [`:1105`](../../../../src/bam/alignment_input.rs#L1105) | reuse the first as-is; the second needs `pub(crate)` or a copy |
| the per-record stamp | `NoodlesRawRecord.source_file_index` [`filtering.rs:382`](../../../../src/ng/read/filtering.rs#L382), stamped in `read_next` [`:438`](../../../../src/ng/read/filtering.rs#L438) | same slot, new payload: the file's `ReadGroupResolution` |
| the region sources | `BamRegionSource` [`region_query.rs:64`](../../../../src/ng/read/input/region_query.rs#L64), which borrows `&'a sam::Header` [`:67`](../../../../src/ng/read/input/region_query.rs#L67) | carry the resolution the same way |
| `ReadFilterCounts` | [`filtering.rs:117`](../../../../src/ng/read/filtering.rs#L117) | reuse the type; change its key (§3.2) |
| `@RG` tag constants | noodles `read_group::tag::{SAMPLE, LIBRARY, PLATFORM}` (`SM`/`LB`/`PL`) | use as-is; there is **no** `SRX` constant — a non-standard tag read through `other_fields()`, as `SAMPLE` already is [`open_bam.rs:861`](../../../../src/ng/read/input/open_bam.rs#L861) |
| `DuplicateReadAcrossFiles` | [`mod.rs:622`](../../../../src/ng/read/input/mod.rs#L622) | unchanged — still scoped to one sample's files (spec §9) |

## 6. Open items

- **`OPEN:` which tag spells the experiment** — the spec's one unresolved design question (spec
  §13). Until it closes, `ReadGroup::experiment` is always `Synthesized` from the library, which is
  a working default, not the answer.
- *Impl-time.* How the reused record buffer holds the resolution and still satisfies its `Default`
  bound ([`filtering.rs:346`](../../../../src/ng/read/filtering.rs#L346)) — an `Arc` refreshed per
  query is the obvious route; confirm when coding.
- *Impl-time.* `PerRecord` as a linear array or a sorted one (§1.3). Measure only if a real
  multi-`@RG` file appears; every file in the surveyed archive declares one.
- *Impl-time.* Whether `counts()` returns pairs or a table indexed by the file's read groups — the
  same shape question `AlignmentFile::counts` already carries ([`alignment_file.md`](alignment_file.md) §7).

## Test & bench shape

Unit tests beside `read_groups.rs` for everything header-shaped: the three hard errors, the
synthesized library name and its collision, the experiment fallback, and identifier determinism
across a shuffled open order. The full obligation list is **spec §11**.

Two that need a fixture rather than a header literal, and belong beside `mod.rs`: a file declaring
several `@RG` naming two samples, opened twice, each open seeing only its own reads with the foreign
ones outside the drop tallies; and the same-`@RG`-file-untagged-records case, which is the universal
path and must read normally.

**No new bench.** Nothing here is on a per-read path in the `Sole` arm, and the merge's existing
per-read budget bench ([`sample_reads.md`](sample_reads.md), T14) already guards the loop these reads
travel; it should keep passing unchanged, which is the regression anchor.
