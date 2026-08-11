# ng — the reference's tandem-repeat catalog: types & interfaces

*Status: architecture draft (2026-08-10), companion to the spec
[`../spec/repeat_catalog.md`](../spec/repeat_catalog.md) (the design and its rationale) and to the
shared arch docs [`ng_step_interfaces.md`](ng_step_interfaces.md) (vocabulary) and
[`module_layout.md`](module_layout.md) (the `src/ng/` tree). Naming follows
[`naming.md`](../../../../ai/skills/rust-code-review/code_review/naming.md): nouns for types, verbs
for functions, **STR** in prose ↔ `ssr` in code. Signatures are illustrative; the **contract** is the
deliverable. Every "why" lives in the spec — this doc does not re-argue one.*

## Module home

`src/ng/repeat_catalog/` — a folder, not a file, because it owns three separable concerns with
different reasons to change (the builder that rides the reference pass, the reader that answers
queries, the Parquet wire format). It is **not a pipeline step**: like `reference_info.rs` and
`ref_seq.rs` it is reference-side infrastructure, so it sits beside them rather than in the numbered
sequence (`module_layout.md`, *The tree*).

```
src/ng/repeat_catalog/
├── mod.rs        – RepeatCatalog (the reader), StrRepeatCriteria, FoundRepeat, RepeatCatalogError
├── criteria.rs   – StrRepeatCriteria + the serves() comparison and its refusal type
├── builder.rs    – RepeatCatalogBuilder: the ReferenceBasesObserver impl (spec §2.3)
└── parquet_file.rs – the schema, the writer, the row-group reader (spec §3.5)
```

Two edits land outside it: the observer seam in `src/ng/reference_info.rs` (§2.1), and lifting
`classify`'s motif/trim/purity arithmetic into callable helpers in
`src/ng/region_typing/segment_criteria.rs` (§5, row *motif/trim/purity*). The CLI subcommand lands in
`src/pop_var_caller_exp/repeat_catalog.rs`, mirroring `typed_regions.rs` (`src/pop_var_caller_exp/mod.rs:8-15`).

---

## 1. Types

### 1.1 `StrRepeatCriteria` — the one policy value

Everything a caller states about which tracts are STRs. It **wraps** `SsrSegmentCriteria` rather than
restating its five fields, so classification's own gates keep one home and cannot drift; the two extra
fields are the ones step 3 has no concept of (spec §5.2).

```rust
/// Which tandem repeats count as STR loci: how short they may be, how long, and how much
/// room they need beside them. One value, used on both sides of the file — the builder is
/// given one and records it, every reader passes one in (spec §5.2).
#[derive(Debug, Clone, PartialEq)]
pub struct StrRepeatCriteria {
    /// Periods, copy floors, purity floor, score floor, bundle radius — step 3's own
    /// admission rules, unchanged and not duplicated.
    pub classification: SsrSegmentCriteria,
    /// Sequence required on each side of a tract, to the contig's end. **Not**
    /// `SsrSegmentCriteria::bundle_threshold`, which drops a locus only when its flank
    /// clamps to zero (`segment_criteria.rs:1126`); this is a floor, and the file is
    /// built at 15 bp (spec §1, §4.1).
    pub min_flank_bp: Bp,
    /// A tract longer than this is a satellite, not a locus. A pure read-time filter over
    /// stored spans — the file caps no length (spec §4.2).
    pub max_str_len_bp: Bp,
}
```

**Defaults are the catalog's, not step 3's** (spec §4.1): `periods` 1..=6, `min_copies`
`[5, 5, 4, 4, 4, 3]`, `min_flank_bp` 15, `max_str_len_bp` 500, with the purity/score/bundle fields
left at `SsrSegmentCriteria::default()`. Each is a named `pub const` with its units and the spec
section that fixed it — never a literal at a use site.

```rust
/// Whether a catalog built under `self` can answer a reader asking `wanted` (spec §4.3).
/// `Ok(())` or the first axis that fails, carrying both values so the error names them.
pub fn serves(&self, wanted: &StrRepeatCriteria) -> Result<(), CriteriaRefusal>;

/// Why a reader cannot be served. One variant per **bounded** axis (spec §4.1); the
/// unbounded ones (purity, score, satellite cap, bundle radius) cannot appear here.
#[non_exhaustive]
pub enum CriteriaRefusal {
    CopyFloor { period: u8, built: u32, wanted: u32 },
    MinFlank  { built: Bp, wanted: Bp },
    PeriodRange { built: PeriodRange, wanted: PeriodRange },
}
```

**Contract.** `serves` compares only what the builder honoured: every period's copy floor, the flank
floor, and the period range. It is a total function with no I/O, and it is the *only* place the
permissiveness rule is expressed — the read methods call it, they do not re-implement it.

### 1.2 `FoundRepeat` — one row

```rust
/// One tandem repeat as the file holds it: what the scanner found, where the whole-motif
/// cut falls, and the two values that needed the bases (spec §3.1, §3.2).
///
/// Coordinates are **1-based inclusive genomic**, converted from the detector's 0-based
/// half-open `RepeatInterval` at the builder's edge — the one site an off-by-one could
/// live (spec §3.1).
#[derive(Debug, Clone, PartialEq)]
pub struct FoundRepeat {
    pub contig: ContigId,
    /// The span the scanner reported. Bundling and the pre-screen use this one.
    pub detected: SpanBp,
    /// The same tract cut back to whole motif copies; `None` when no clean cut exists
    /// (`minimal_trim` fails — `RejectionReason::NoCleanTrim`). **This is the locus**, and
    /// the copy floor, the purity floor, the flank test and the stratum count all read it.
    pub trimmed: Option<SpanBp>,
    pub period: u8,
    pub score: i32,
    pub motif: Motif,
    /// Fraction of the **trimmed** tract matching a perfect motif tiling; `None` exactly
    /// when `trimmed` is.
    pub purity: Option<f32>,
}

/// A 1-based inclusive span on one contig. Ng has `GenomeRegion` for a named region; this
/// is the bare pair a row carries, and it exists so the two spans above cannot be swapped.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SpanBp { pub start: Position, pub end: Position }
```

`trimmed` and `purity` are `Option`, never a sentinel span or `-1` purity: a tract with no clean trim
is a real thing the file records, because it can still be a **bundle member** for its neighbour
(`segment_criteria.rs:973-977` bundles before `finish_locus` runs).

### 1.3 The header

```rust
/// What the file says about itself, read without touching a row (spec §3.4).
#[derive(Debug, Clone, PartialEq)]
pub struct RepeatCatalogHeader {
    /// One entry per contig, in reference order: name, length, MD5 — the `@SQ M5`
    /// (`reference_info.rs:51-67`). Reused verbatim; geometry fields carry through unused.
    pub contigs: Vec<ContigInfo>,
    /// The whole-reference digest (`reference_info.rs:71-80`).
    pub reference_md5: [u8; 16],
    /// What the builder honoured: periods, copy floors, min flank. The other four fields of
    /// `StrRepeatCriteria` are read-time filters and are recorded for provenance only.
    pub built_under: StrRepeatCriteria,
    /// The scoring weights the tracts came out of. **Equal or refuse** (spec §4.2).
    pub scan: ScanParams,
    /// Crate version; a detector change invalidates the file even at identical settings.
    pub tool_version: String,
}
```

### 1.4 Errors

```rust
/// Everything that can stop a catalog being read or built. Fatal by design: none of these
/// is a warning, a fallback or a silent rebuild (spec §4.3, §6).
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RepeatCatalogError {
    /// No file at that path — the caller can act on this one, so it names the command
    /// that writes it (spec §5.3).
    #[error("no repeat catalog at {path}; build one with `pop_var_caller_exp repeat-catalog`")]
    NotFound { path: PathBuf },
    /// A contig's MD5 differs from the one this run computed from the FASTA (spec §4.3).
    #[error("catalog does not describe this reference: contig {contig} digests {found}, catalog says {expected}")]
    ContigDigestMismatch { contig: String, found: String, expected: String },
    /// A contig present in one and absent from the other, or the order differs.
    #[error("catalog contig table does not match the reference: {detail}")]
    ContigTableMismatch { detail: String },
    /// The reader asked about tracts the file does not hold (spec §4.1, §4.3).
    #[error("catalog is not permissive enough: {0}")]
    CriteriaTooPermissive(#[from] CriteriaRefusal),
    /// Different scoring weights are a different set of tracts, not a subset (spec §4.2).
    #[error("catalog was scored with {built:?}, reader asked for {wanted:?}")]
    ScoringWeightsDiffer { built: ScanParams, wanted: ScanParams },
    /// Built by a different version of the detector.
    #[error("catalog was written by tool version {built}, this is {running}")]
    ToolVersionDiffers { built: String, running: String },
    /// The file is unreadable, has no footer (a build that died mid-write), or its schema
    /// is not this one.
    #[error("catalog file {path} is unreadable or truncated")]
    Unreadable { path: PathBuf, source: Box<dyn std::error::Error + Send + Sync> },
    #[error("writing the catalog to {path}")]
    Write { path: PathBuf, source: std::io::Error },
}
```

**Contract.** Every variant names the axis or the contig; none is constructed without both values.
`Unreadable` is what a truncated Parquet file produces on open — the missing footer is the check
(spec §3.5), so there is no "short but valid" state to detect separately.

---

## 2. Interfaces

### 2.1 The observer seam, in `reference_info.rs`

The pass keeps its leaf property (`reference_info.rs:12-15` — it knows nothing about ng): it hands
uppercased bases forward and learns nothing about repeats.

```rust
/// Something that wants each contig's bases as the reference pass streams them by
/// (spec §2.2). Bases arrive **uppercased, in coordinate order, in whatever chunks the
/// pass batched** — never a whole contig, so an observer that needs one accumulates it.
///
/// Deliberately **infallible**: the pass is a leaf and cannot own an observer's error type.
/// An observer that fails records it and no-ops for the rest of the pass, then reports at
/// its own `finish` (see `RepeatCatalogBuilder`).
pub trait ReferenceBasesObserver {
    /// A new contig has started; `index` is its position in reference order.
    fn contig_started(&mut self, name: &str, index: usize);
    /// The next bases of the current contig, uppercased.
    fn bases(&mut self, upper: &[u8]);
    /// The current contig is complete, with the geometry and digest the pass reconstructed.
    fn contig_finished(&mut self, info: &ContigInfo);
}

/// `read_reference_info` with an observer attached (spec §2.2). The existing entry point
/// becomes this one with a no-op observer, so no caller changes.
pub fn read_reference_info_observing(
    source: ReferenceSource,
    observer: &mut dyn ReferenceBasesObserver,
) -> Result<ReferenceInfo, ReferenceInfoError>;
```

**Contract.** Calls are ordered `contig_started → bases* → contig_finished`, once per contig, in file
order. `bases` is called only for `ReferenceSource::Fasta` (a `.fai` has no bases); a `.fai`-only read
calls none of the three. The observer sees the same bytes the digests see — `FastaPass::flush_md5`
(`reference_info.rs:664-671`) is the single site that feeds both.

### 2.2 The builder

```rust
/// Accumulates one contig, scans it whole, and writes its rows (spec §2.3, §2.4).
pub struct RepeatCatalogBuilder { /* criteria, scan params, writer, per-contig buffer */ }

impl RepeatCatalogBuilder {
    pub fn new(out: &Path, criteria: StrRepeatCriteria, scan: ScanParams) -> Result<Self, RepeatCatalogError>;
    /// Flush the last row group, write the footer metadata, rename into place. Returns the
    /// per-period row tally the CLI prints (spec §2.6).
    pub fn finish(self, info: &ReferenceInfo) -> Result<RowsByPeriod, RepeatCatalogError>;
}

impl ReferenceBasesObserver for RepeatCatalogBuilder { /* … */ }
```

**Contract.** `contig_finished` is where the work happens: `find_tandem_repeats` over the accumulated
contig (`tandem_repeat.rs:483`), then per interval the copy floor **on the detected span** — as
`prefilter` measures it (spec §3.1) — the trim, the motif, the purity, and the flank floor, then rows
appended in `(start, period, end)` order. The buffer is cleared before the next contig, so nothing
accumulates across them. With `--threads N` up to N contigs are in flight and completed contigs are
written **in reference order**, never in completion order.

### 2.3 The reader

```rust
/// An opened catalog: its header, validated against this run's reference, plus the handle
/// the row methods read through (spec §5.3).
pub struct RepeatCatalog { /* header, file handle */ }

impl RepeatCatalog {
    pub fn open_checking_against_reference(path: &Path, reference: &ReferenceInfo)
        -> Result<Self, RepeatCatalogError>;

    pub fn header(&self) -> &RepeatCatalogHeader;

    /// Every repeat the file holds there, as stored, no criteria applied.
    pub fn repeats_in_region(&self, region: Option<&GenomeRegions>)
        -> impl Iterator<Item = Result<FoundRepeat, RepeatCatalogError>> + '_;

    /// The genome's segments in coordinate order — STR segments and the generic spans
    /// between them (spec §5.1).
    pub fn genome_segments(&self, criteria: &StrRepeatCriteria, region: Option<&GenomeRegions>)
        -> Result<impl Iterator<Item = Result<TypedRegion, RepeatCatalogError>> + '_, RepeatCatalogError>;

    /// The surviving STR loci alone, without the generic spans between them.
    pub fn str_loci(&self, criteria: &StrRepeatCriteria, region: Option<&GenomeRegions>)
        -> Result<impl Iterator<Item = Result<SsrSegment, RepeatCatalogError>> + '_, RepeatCatalogError>;

    /// How many loci in each (period, repeat count) stratum.
    pub fn count_loci_per_stratum(&self, criteria: &StrRepeatCriteria, region: Option<&GenomeRegions>)
        -> Result<StratumCounts, RepeatCatalogError>;

    /// Up to `cap` loci per stratum — the ones whose `hash(contig, start, seed)` is lowest —
    /// and the full counts, from **one** pass (spec §5.3).
    pub fn sample_loci_per_stratum(
        &self, criteria: &StrRepeatCriteria, region: Option<&GenomeRegions>, cap: u32, seed: u64,
    ) -> Result<(StratumCounts, StratumSample), RepeatCatalogError>;
}
```

**Contract, and it is the same for all four criteria-taking methods.** The criteria check runs
**once, eagerly**, at the call — hence `Result<impl Iterator, _>` rather than an iterator of results
that refuses on its first item: a refusal is about the file and the policy, not about a row. After
that the iterators are lazy and stream a row group at a time, so peak memory is one row group plus
the caller's own state. `region = None` means the whole reference. Rows arrive in coordinate order,
and `genome_segments` covers its region with no gap — the property `partition_resident` has
(`region_typing/mod.rs:380`) and the differential checks.

**`sample_loci_per_stratum` returns the counts too**, because the pre-pass needs both and the tally
is a counter beside the heaps rather than a second traversal (spec §5.3). Its working state is
`cap` values per stratum in a bounded max-heap, and the result is order-independent: merging two
shards is taking the lowest `cap` of the union.

---

## 3. The file (spec §3.5)

One row group per contig; the header of §1.3 in the footer's key-value metadata, JSON-encoded under
one key so a `duckdb` user can read it too.

| column | Arrow type | encoding |
|---|---|---|
| `contig` | `Dictionary(UInt16, Utf8)` | dictionary — tens of distinct values per genome |
| `detected_start`, `detected_end` | `UInt64` | delta; ascending within a row group |
| `trimmed_start`, `trimmed_end` | `UInt64`, nullable | delta |
| `period` | `UInt8` | RLE |
| `score` | `Int32` | plain |
| `motif` | `Dictionary(UInt16, Utf8)` | dictionary — at most 5,460 distinct primitive motifs of period 1..6 (`MAX_MOTIF_LEN = 6`, `types.rs:286`) |
| `purity` | `Float32`, nullable | plain |

**Fixed by us, not left to the crate's defaults** (spec §6, determinism): the compression codec and
level, one row group per contig, and the writer-version string. A default that moves under a crate
upgrade would change the bytes without changing the content, and §10.5 would start failing for a
reason that is not a bug in this code.

**Arrow types do not leave `parquet_file.rs`.** The reader hands back `FoundRepeat`, `SsrSegment` and
`TypedRegion`; nothing in `mod.rs`'s signatures mentions Arrow, which is what keeps a format change
one module deep.

---

## 4. Design decisions — decided

- **`StrRepeatCriteria` wraps `SsrSegmentCriteria`; it does not restate its fields.** Five of the
  seven values are step 3's own admission rules, and duplicating them is how two copies of a policy
  drift. *Rejected:* a flat seven-field struct plus a `to_segment_criteria()` — it reads better at a
  call site and it is the drift. — spec §5.2
- **The builder's copy floor is measured on the detected span, as `prefilter` does.** Measuring it on
  the trimmed span would drop rows that a reader's pre-screen keeps as bundle members, and bundling
  would differ from a live scan — silently, and only near clusters. — spec §3.1
- **Both spans are stored.** Bundling needs the detected one, everything downstream needs the trimmed
  one, and a row carrying one would send the reader back to the FASTA for the other. — spec §3.1
- **The observer trait is infallible; the builder holds its own first error.** A fallible trait would
  put an error type from `repeat_catalog` into `reference_info`, which is a leaf and must stay one.
  — spec §2.2
- **The criteria check is eager, the row iteration is lazy.** A refusal is a fact about the file, not
  about a row; discovering it on `next()` would let a caller start a loop that cannot run. — spec §5.3
- **`open_checking_against_reference` takes `&ReferenceInfo`, not a path.** The digests it checks
  against are the ones this run computed from the FASTA; taking a path would invite a second read and
  a second source of truth. — spec §4.3
- **Contig identity is `ContigId` (the header's index), resolved on open.** Rows store the index, not
  the name; the contig table is what maps it, and it is validated before any row is read. — spec §3.4
- **No trait, no bake-off.** One builder, one reader, one file format. `RepeatCatalog` is a concrete
  type. — `module_layout.md` principle 1a

---

## 5. Reconciliation with existing code

Every row was read at the cited line.

| this doc | existing code | how they meet |
|---|---|---|
| detection | `find_tandem_repeats` ([`tandem_repeat.rs:483`](../../../../src/ng/tandem_repeat.rs)) | called as-is, once per whole contig; `scan_windowed` (`:1003`) is the fallback of spec §2.3, not used |
| the detector's row | `RepeatInterval` ([`tandem_repeat.rs:208-219`](../../../../src/ng/tandem_repeat.rs)) | 0-based half-open; converted to `FoundRepeat`'s 1-based inclusive spans at the builder's edge |
| scoring weights | `ScanParams` ([`tandem_repeat.rs:123-135`](../../../../src/ng/tandem_repeat.rs)) | stored in the header verbatim; `min_copies` there stays the scanner's permissive floor (the table's minimum, 3) |
| period range | `PeriodRange` ([`tandem_repeat.rs:55`](../../../../src/ng/tandem_repeat.rs)) | reused; the catalog builds at 1..=6 |
| copy floors | `MinCopies` ([`segment_criteria.rs:355-362`](../../../../src/ng/region_typing/segment_criteria.rs)) | reused; the catalog's table is `[5, 5, 4, 4, 4, 3]` against `MinCopies::default()`'s `[8, 6, 6, 6, 5, 4]` (`:444`) |
| classification policy | `SsrSegmentCriteria` ([`segment_criteria.rs:478-540`](../../../../src/ng/region_typing/segment_criteria.rs)) | wrapped by `StrRepeatCriteria` (§1.1) |
| pre-screen, overlap resolution | `prefilter` ([`segment_criteria.rs:677`](../../../../src/ng/region_typing/segment_criteria.rs)) | called by the reader on rows from the file instead of on a live scan |
| admission | `classify` ([`segment_criteria.rs:849`](../../../../src/ng/region_typing/segment_criteria.rs)) → `Classified` (`:717-733`) | the reader's segmentation runs it; the builder does **not** — it needs only the criteria-free parts |
| motif / trim / purity | inside `finish_locus` ([`segment_criteria.rs:1011-1097`](../../../../src/ng/region_typing/segment_criteria.rs)): the `minimal_trim` call (`:1031`, defined `:1217`), `Motif::new` (`:1083`), `recompute_purity` (`:1097`) | **lifted into `pub(crate)` helpers** so the builder and `finish_locus` share one implementation; nothing is re-derived here |
| the locus type | `SsrSegment` ([`segment_criteria.rs:141-153`](../../../../src/ng/region_typing/segment_criteria.rs)) | `str_loci` returns these, built through `SsrSegment::new` so its coordinate invariant still guards them |
| the segment type | `TypedRegion` / `RegionKind` ([`region_typing/mod.rs:146`, `:170`](../../../../src/ng/region_typing/mod.rs)) | `genome_segments` returns these — the same four kinds, from a file instead of a scan |
| the parity oracle | `partition_resident` ([`region_typing/mod.rs:380`](../../../../src/ng/region_typing/mod.rs)) | §10.1's differential runs it over the same reference at the same policy |
| the satellite cap's other home | `TypedRegionConfig::max_str_len` ([`region_typing/mod.rs:209-232`](../../../../src/ng/region_typing/mod.rs)) | there it is the cap **and** the window margin; here there are no windows, so `StrRepeatCriteria::max_str_len_bp` is the cap alone |
| the reference pass | `read_reference_info` ([`reference_info.rs:270`](../../../../src/ng/reference_info.rs)), `FastaPass` (`:516`), `flush_md5` (`:664`) | gains the observer seam of §2.1; the digests and geometry are untouched |
| contig table | `ContigInfo` ([`reference_info.rs:51-67`](../../../../src/ng/reference_info.rs)), `ReferenceInfo` (`:71-93`) | stored in the header verbatim — name, length, MD5 (`@SQ M5`) |
| sibling-file pattern | `sibling_fai_path` ([`reference_info.rs:796`](../../../../src/ng/reference_info.rs)), `write_fai` (`:816`) | the same shape: derived file beside the reference, written atomically, created on request |
| error style | `RefSeqError` ([`ref_seq.rs:41`](../../../../src/ng/ref_seq.rs)), `ReferenceInfoError` ([`reference_info.rs:183`](../../../../src/ng/reference_info.rs)) | `#[non_exhaustive]` `thiserror` enum, one doc comment per variant saying when it fires |
| region restriction | `GenomeRegions` ([`region_typing/mod.rs:78-80`](../../../../src/ng/region_typing/mod.rs)), wrapping `RegionSet` ([`regions.rs:74`](../../../../src/regions.rs)) | the `region` argument; `None` = whole reference |
| scalars | `ContigId`, `Position`, `Bp`, `Motif` ([`types.rs:13`, `:34`, `:174`, `:338`](../../../../src/ng/types.rs)) | reused; no new scalar newtype is minted except `SpanBp` (§1.2) |
| the CLI | `PopVarCallerExpCommand` ([`pop_var_caller_exp/cli.rs:22`](../../../../src/pop_var_caller_exp/cli.rs)), `typed_regions.rs` | a second variant `RepeatCatalog(RepeatCatalogArgs)`, one module beside `typed_regions` |

---

## 6. Open items

- `OPEN:` **whether `classify`'s coordinate arithmetic is criteria-independent** — spec §9.4. The
  stored trim, motif and purity are only sound if it is. Settled by the differential (§7) run at
  several policies, not by reading the code again.
- **Impl-time confirmation, not a decision:** the exact Parquet writer knobs (codec, level, dictionary
  page size) that make two builds byte-identical — pinned when the writer is written, asserted by the
  thread-count test.
- **Impl-time confirmation:** whether `minimal_trim`, `recompute_purity` and the motif slice lift
  cleanly out of `finish_locus` or need a small struct to carry the trimmed slice with them.

---

## 7. Test and bench shape

Tests live in per-module `#[cfg(test)]` blocks; the differential and the CLI test are integration
tests under `tests/`, since they need a fixture reference on disk.

**The regression anchor is `partition_resident`** (spec §10.1): the segmentation derived from a
catalog must equal a live scan's, run at several policies including ones differing from the build
settings on every bounded axis, and including a fixture with overlapping detections, a 2 kb tract, and
repeats at both contig ends. No `bench/` — the build runs once per reference and the only numbers
wanted (file size, row count by period, wall clock, peak memory) come out of the CLI itself (spec
§10.4).
