# ng — the psp store: types and interfaces

*Architecture draft, 2026-08-26. Companion to [`../spec/psp_file_format.md`](../spec/psp_file_format.md)
(the container) and [`../spec/psp_record_encoding.md`](../spec/psp_record_encoding.md) (what a record
holds). Every **why** lives there; this doc adds the code shape and nothing else.*

*Module-tree rules: [`module_layout.md`](module_layout.md). **That document's §"Crate boundary" says
ng "does not adopt the `.psp` file yet — the artifact starts in memory and gains a serialization when
memory forces it". Memory forced it** ([the memory review](../../reports/reviews/psp_memory_milestone_z_2026-08-25.md)),
so this doc is the serialization it deferred; that paragraph should be corrected when it is next
touched.*

---

## 1. Module home

```
src/ng/psp/
├── mod.rs        – PspReader, PspWriter, the free functions, the error enums, re-exports
├── header.rs     – Header, Manifest, FieldEncoding: build, encode, parse, validate
├── block.rs      – the psp block: cutting it, compressing it, streaming it back
├── record.rs     – one SampleLocusObservations to bytes and back; the RecordHead
├── index.rs      – BlockIndex, BlockIndexEntry
├── footer.rs     – the fixed tail: offsets, checksum, magic
└── chain_ids.rs  – the live set and its changes (spec psp_chain_id_encoding.md)
```

**`psp/`, mirroring production's `src/psp/`** (the owner, 2026-08-26). *An earlier draft called it
`store/`, reasoning that a module should be named for what it owns rather than for a file
extension. That was wrong here twice over: `psp` **is** the project's name for this artifact — it
titles all three specs, it is the extension and it is the type prefix — and the parallel with
`src/psp/` tells a reader moving between the two what the relationship is, which `store/` would
have hidden.*

**A folder, not a file**, because the container splits along seams a reader crosses one at a time:
opening touches `footer` then `index` then `header` and no block; a walk touches `block` and
`record` and neither of the others. Those are different concerns with different tests.

**It is not a pipeline step**, so there is no `LocusGenerator`-style trait and no bake-off shape. It
is infrastructure beside `ref_seq.rs`, used by the locus generator (writing) and by the cohort merge
(reading).

**`chain_ids.rs` is its own file because it is its own problem.** It is the only field with a spec of
its own, it holds state across records within a block, and it is where the silent failure lives
(§5.3).

---

## 2. The record, and what differs from the prototype

**⚠ The measuring prototype encodes production's `PileupRecord`. ng's record is
`SampleLocusObservations` ([`src/ng/locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs)),
and it is not the same shape.** Every byte figure in the specs was taken on the first; the second
carries three things the first has no equivalent of, so **expect the measured bytes-a-record numbers
to move, and do not treat them as predictions of ng's file.**

| | production `PileupRecord` | ng `SampleLocusObservations` |
|---|---|---|
| extent | one anchor `pos` | a `GenomeRegion` — **the span is a field, not derived** |
| reference bases | absent | `reference_bases: Box<[u8]>` — spec §4 leans to dropping it |
| per observation | allele sequence + support | **+ `read_witness`, + `read_group`** |
| reads that showed nothing | absent | `reads_without_observation`, `reads_discarded_by_cap` |

**`read_witness` is the one that bites.** It is `Complete`, or a set of runs
([`src/ng/locus_generation/witness.rs:259`](../../../../src/ng/locus_generation/witness.rs)) — so an
observation's identity carries a variable-length field the prototype never encoded, and it must be in
the record body.

**The span being a field is a simplification, not a cost.** The record head needs it (spec §4.3) and
`region` already holds it, where the prototype had to take the reference allele's length as a
stand-in.

---

## 3. The types

### 3.1 The record head

The fixed fields at the front of every record that let a reader decide whether it wants the body
without decoding it. Spec [`psp_file_format.md`](../spec/psp_file_format.md) §4.3 for why each is
here.

```rust
/// What a reader learns about a record before deciding to build it.
///
/// Every field is read from the record's head; none requires touching the body.
/// `body_bytes` is what makes skipping possible at all — the encoded bytes carry no
/// separators, so without it a reader must decode every integer in a record to find
/// where the next one starts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordHead {
    /// Where the record sits. Absolute, rebuilt by the reader from the block's first
    /// position and the encoded difference.
    pub region: GenomeRegion,
    /// Reads supporting anything other than the reference at this locus, summed over
    /// the non-residual observations. Zero and "no alternative allele here" are the
    /// same condition.
    pub non_reference_reads: u32,
    /// Length of the body that follows, in bytes.
    pub body_bytes: u32,
}
```

**`GenomeRegion` rather than a position and a span** — the newtype exists
([`src/ng/types.rs:79`](../../../../src/ng/types.rs)) and holds exactly this pair, so a head that
carried two integers would be minting a duplicate.

### 3.2 The header and its manifest

```rust
/// A file's own account of how it was written. Everything a reader needs before it
/// touches a block, and the only part of the file that is plain text.
pub struct Header {
    /// `MAJOR.MINOR`. **Parsed before anything else and never behind a binary
    /// encoding** — a reader must be able to learn the version of a file it cannot
    /// otherwise read.
    pub format_version: (u16, u16),
    pub sample: String,
    pub reference: ReferenceInfo,
    pub contigs: Vec<ContigEntry>,
    pub writer: WriterProvenance,
    pub manifest: Manifest,
}

/// How this file encodes what it holds. Every value here is the writer's choice,
/// recorded so a reader is driven by the file rather than by an assumption.
pub struct Manifest {
    /// Reference bases per psp block. The cut is a grid on the coordinate: a block ends
    /// when a position crosses into the next multiple.
    pub genomic_block_size_bp: Bp,
    /// A secondary cut, when one span holds too much at depth. `None` = no ceiling.
    pub block_byte_ceiling: Option<u32>,
    /// zstd's look-back window as its exponent. **A reader configures its decoder from
    /// this**; assuming a value makes it reject legitimate files.
    pub window_log: u8,
    /// One entry per field, in encoding order.
    pub fields: Vec<FieldSpec>,
}

/// A field, and how to read it. The set of encodings is closed (§3.3) so a reader stays
/// a match rather than a plugin host.
pub struct FieldSpec {
    pub name: FieldName,
    pub encoding: FieldEncoding,
}
```

**`Bp` and `ContigId` are ng's existing newtypes** ([`src/ng/types.rs:185,13`](../../../../src/ng/types.rs));
this doc mints neither.

### 3.3 Field encodings — a closed set

```rust
/// How one field's values are laid down. **Closed deliberately**: an open-ended scheme
/// would buy flexibility nobody asked for and cost speed on a path that runs at about
/// twenty million records a second (spec psp_file_format.md §4.5).
#[non_exhaustive]
pub enum FieldEncoding {
    /// LEB128.
    Varint,
    /// Zig-zag LEB128, for values that go negative.
    SignedVarint,
    /// A fixed-width little-endian integer of `bytes` width.
    Fixed { bytes: u8 },
    /// Raw IEEE bytes. **The escape hatch, not the default** — a float stored raw is a
    /// float two callers can disagree about in its last bits.
    Ieee { bytes: u8 },
    /// A count of steps of `1 / scale`. **The step is inherited from the type that
    /// produced the value, not chosen by the writer** (spec psp_record_encoding.md
    /// §5.1.1): the psp records it so a reader can interpret the integer, and cannot
    /// write a file with a step the types did not produce.
    FixedPoint { scale: u32 },
    /// Bytes with a varint length in front.
    LengthPrefixedBytes,
}
```

**`FixedPoint` carries a scale it does not own**, and that is the whole point of the ruling: the
rounding happens where the value is computed, so direct mode and psp mode see the same number. The
field exists to *interpret*, not to *decide*.

### 3.4 The block index and the footer

```rust
/// One entry per psp block, in genomic order.
///
/// **And nothing more.** An earlier draft carried the largest non-reference support in
/// the block so a reader could skip whole blocks; at a 100 kb block essentially every
/// block contains something, so the field never fired (spec psp_record_encoding.md §2.4).
pub struct BlockIndexEntry {
    pub contig: ContigId,
    pub first_position: GenomePosition,
    pub offset: u64,
}

/// The fixed tail. **Its presence is what says the file is complete**; a reader that
/// cannot find it refuses the file rather than reading the blocks that reached disk.
pub struct Footer {
    pub index_offset: u64,
    pub index_bytes: u64,
    pub trailer_offset: u64,
    pub trailer_bytes: u64,
    pub n_blocks: u64,
    pub index_checksum: u32,
    // magic last, so a four-byte read at end-of-file rejects a foreign or truncated
    // file before anything else is touched
}
```

---

## 4. The interfaces

### 4.1 Reading

```rust
impl PspReader {
    /// Footer, then index, then header. **No block is touched**, so this is the cost a
    /// cohort pays per open sample before reading anything.
    pub fn open(path: &Path) -> Result<Self, PspReadError>;

    pub fn header(&self) -> &Header;
    /// The writer's closing payload. Opaque here; the caller interprets it.
    pub fn trailer(&mut self) -> Result<&[u8], PspReadError>;
    pub fn blocks(&self) -> &[BlockIndexEntry];

    /// Every record, from the first block.
    pub fn records(&mut self) -> RecordIter<'_>;
    /// From the block holding `position`. **Reading starts at that block's first
    /// record**, not at `position` — a reader cannot start mid-block.
    pub fn records_from(&mut self, contig: ContigId, position: GenomePosition)
        -> Result<RecordIter<'_>, PspReadError>;

    /// The cohort's first pass: `want` sees each record's head and says whether to
    /// build the body. A record the predicate declines costs the head plus a pointer
    /// advance.
    pub fn records_where<F>(&mut self, want: F) -> SelectiveIter<'_, F>
    where
        F: FnMut(&RecordHead) -> bool;
}
```

**Contract.** Both iterators are lazy and borrow the reader; neither retains a record it has
yielded. A reader holds, per open file: the compressed read buffer, the rolling decompressed buffer
(16 kB each — spec §4.4), the decompressor's state, and the record being built. **Nothing it holds is
a function of the block size**, which is goal 1 and the reason the design exists. An iterator that
fails yields `Err` once and then `None`; it never yields a half-built record.

**`records_where` is the only way to skip.** A caller that filters the output of `records()` has
already paid to build every record, which is the cost the head exists to avoid.

### 4.2 Writing

```rust
impl PspWriter {
    /// The manifest is fixed here and cannot change for the life of the file.
    pub fn create(path: &Path, header: Header) -> Result<Self, PspWriteError>;

    /// **Coordinate order is enforced**: a record that starts before the previous one
    /// is `OutOfOrder`, not a file that seeks wrongly.
    pub fn push(&mut self, record: &SampleLocusObservations) -> Result<(), PspWriteError>;

    /// Index, then trailer, then footer — then flush, surface the buffered writer's
    /// errors, and sync. **Consumes the writer**: there is no way to hold one that has
    /// finished, and no way to produce a readable file without calling this.
    pub fn finish(self, trailer: &[u8]) -> Result<WriteStats, PspWriteError>;

    /// Reopen a finished file and add records. Truncates at the index offset, keeping
    /// the header and therefore the manifest — so appended records must use the
    /// encodings already declared, and **the old trailer is discarded with the index**.
    pub fn append(path: &Path) -> Result<Self, PspWriteError>;
}

/// Replace a finished file's trailer, touching neither blocks nor index. Cheap because
/// the index sits before the trailer (spec §3).
pub fn replace_trailer(path: &Path, trailer: &[u8]) -> Result<(), PspWriteError>;

/// The header of a file that may have no footer. Works where `PspReader::open` correctly
/// fails, which is its purpose.
pub fn read_header(path: &Path) -> Result<Header, PspReadError>;
```

**Contract.** A `PspWriter` dropped without `finish` leaves a file with no footer — refused by every
reader, which is the intended outcome for a killed run and **must not be softened into a partial
read**. `finish` consuming `self` is what makes that unambiguous in the type system.

### 4.3 Errors

```rust
/// **Every variant is an input problem, not a bug** — a corrupt file is data, so none of
/// these is a panic.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum PspReadError {
    /// No valid footer: the writer was interrupted. Rebuild the file.
    #[error("{path} has no valid footer — the writer did not finish")]
    Incomplete { path: PathBuf },
    /// Written by a newer format. Upgrade the reader.
    #[error("{path} is format {found:?}; this reader understands up to {supported:?}")]
    UnsupportedVersion { path: PathBuf, found: (u16, u16), supported: (u16, u16) },
    /// The file's look-back window is larger than this reader budgeted for. **A distinct
    /// variant because the fix is a knob, not a rebuild** — and because zstd's own error
    /// here says nothing a user can act on.
    #[error("{path} needs a {needed_bytes}-byte window; this reader allows {allowed_bytes}")]
    WindowTooLarge { path: PathBuf, needed_bytes: u64, allowed_bytes: u64 },
    /// A block failed to decompress, or a record ran past its block.
    #[error("{path}: block {block} is corrupt")]
    CorruptBlock { path: PathBuf, block: u64, source: std::io::Error },
    #[error(transparent)]
    Io(#[from] std::io::Error),
}
```

**The four classes are separate because the instructions differ** — rebuild the file, upgrade the
reader, raise a limit, the data is damaged. Collapsing them into one loses that.

---

## 5. Design decisions — decided

- **The module is `src/ng/psp/`, a folder, not a pipeline step.** No trait, no bake-off: there is
  one implementation. §1.
- **`RecordHead` carries `GenomeRegion`, not a position and a span** — the newtype exists and holds
  that pair. §3.1.
- **The field-encoding set is closed.** Speed on the hot path over flexibility nobody asked for —
  spec [`psp_file_format.md`](../spec/psp_file_format.md) §4.5.
- **`FixedPoint`'s scale is recorded, not chosen.** The rounding is upstream in the type, which is
  what keeps direct mode and psp mode bit-identical — spec
  [`psp_record_encoding.md`](../spec/psp_record_encoding.md) §5.1.1.
- **`finish` consumes the writer.** A file with no footer is refused, and the type system should say
  so rather than a doc comment. §4.2.
- **`records_where` is the skip, and filtering `records()` is not.** A caller that filters afterwards
  has already paid to build everything — spec [`psp_file_format.md`](../spec/psp_file_format.md) §6.2.
- **The chain ids' live-set changes ride in the record head; the exception lists stay in the body.**
  The changes carry state across records and a skipped body would strand them; the exception lists
  carry none — spec [`psp_record_encoding.md`](../spec/psp_record_encoding.md) §6.
- **Records are serialised by hand, not through `serde`.** A derive fixes the layout at compile time
  and the manifest requires it per-file — spec
  [`psp_file_format.md`](../spec/psp_file_format.md) §4.6. `serde` stays for the header's TOML and
  for the trailer's payload, neither on a hot path.

---

## 6. Reconciliation with existing code

Every row read before it was written down.

| what this doc names | existing code | how it is used |
|---|---|---|
| the header's framing | [`HEAD_MAGIC`, `src/psp/header.rs:53`](../../../../src/psp/header.rs) | **the pattern, re-implemented**: magic, `u64` length, TOML body, sentinel — so `head` still works on an ng file |
| `Header` / `Manifest` | `WriterHeader` / `ParsedHeader`, [`src/psp/header.rs:88,151`](../../../../src/psp/header.rs) | the build-side / parse-side split, kept. ng's adds `Manifest`; production's `[[column]]` array is the ancestor |
| `FieldEncoding` | [`src/psp/registry.rs`](../../../../src/psp/registry.rs)'s `Cardinality` and element types | **the idea, widened.** Production's registry names a type; this names a type *and its parameters* |
| `Footer` | `TRAILER_BYTES = 32`, `TRAILER_MAGIC = PSPE`, [`src/psp/trailer.rs:28,33`](../../../../src/psp/trailer.rs) | the layout and the magic-last trick. **Wider than 32**: ng's needs the trailer's offset and length too |
| `BlockIndexEntry` | [`src/psp/index.rs:42`](../../../../src/psp/index.rs) | same shape, same flat vector decoded at open. ng's drops nothing and adds nothing |
| varints | `encode_u64_leb128` / `decode_u64_leb128` / `encode_i64_svarint`, [`src/psp/varint.rs:46,83,119`](../../../../src/psp/varint.rs) | **called as-is.** Specified and tested; ng mints no second implementation |
| the compression seam | `ZSTD_COMPRESSION_LEVEL = 9`, [`src/psp/block.rs:709`](../../../../src/psp/block.rs) | the level and the long-lived-compressor shape. **The window cap is new** — production never set one |
| `ng::psp::PspReadError` / `PspWriteError` | production's same-named pair, [`src/psp/errors.rs:204,596`](../../../../src/psp/errors.rs) | the `#[non_exhaustive]` `thiserror` shape and the per-variant doc comment. **Same names, different types, different modules** — ng's surface is its own, and the two must not be confused at a `use` site |
| the record | `SampleLocusObservations`, [`src/ng/locus_generation/mod.rs:40`](../../../../src/ng/locus_generation/mod.rs) | **what is written and what must come back.** Not production's `PileupRecord` — §2 |
| `ReadWitness` | [`src/ng/locus_generation/witness.rs:259`](../../../../src/ng/locus_generation/witness.rs) | encoded as-is; the variable-length field the prototype never carried |
| `GenomeRegion`, `ContigId`, `Bp`, `ReadGroupId` | [`src/ng/types.rs:79,13,185,210`](../../../../src/ng/types.rs) | used, not re-minted |
| `ChainId` | `pub type ChainId = u64`, [`src/pileup_record.rs:30`](../../../../src/pileup_record.rs) | **an alias, not a newtype, and it is production's.** OPEN below |
| the working prototype | `examples/psp_row_stream_roundtrip.rs` | the parity oracle (§8) and the source of every measurement the specs quote |

---

## 7. Open items

- **`OPEN:` should `ChainId` become an ng newtype?** It is `pub type ChainId = u64` in production
  ([`src/pileup_record.rs:30`](../../../../src/pileup_record.rs)), so it is an alias with no type
  safety, and ng imports it. The store's chain-id path mixes ids with *positions in the live set*
  (spec [`psp_chain_id_encoding.md`](../spec/psp_chain_id_encoding.md) §4) — two `u64`s that must not
  be transposed, which is precisely the case for a newtype. **Not this module's to decide**: it
  belongs with whoever owns ng's chain-id minting.
- **Impl-time: the exact width of `Fixed`/`Ieee`, and whether the head's fields are fixed or
  varint.** Spec [`psp_file_format.md`](../spec/psp_file_format.md) §4.3 leaves it to the manifest
  and records that a fixed width costs less than it looks after compression, unmeasured in place.
  Pin it with a measurement at implementation time, not by argument.
- **Impl-time: whether `reference_bases` is stored.** Spec
  [`psp_record_encoding.md`](../spec/psp_record_encoding.md) §4 leans to dropping and re-fetching;
  nobody has timed the re-fetch.
- **`OPEN:` how much of the record head's speed-up survives at 300 reads a position.** The chain
  ids' live-set changes ride in the head, and they grow with depth — 0.432 bytes a position at 11.4
  reads, 6.42 at 293 — so the head grows while the skip's value shrinks. Spec
  [`psp_record_encoding.md`](../spec/psp_record_encoding.md) §6. **A measurement, not a design
  question**, and the first one to take once a writer exists.

---

## 8. Test and bench shape

**Unit tests live beside their module** (`#[cfg(test)] mod tests` per file), as the rest of ng does.

**The parity oracle is `examples/psp_row_stream_roundtrip.rs`'s `verify`**, promoted to a test: it
walks a written store and its source in lockstep and fails on the first record that disagrees, with
the strictness the fields require — every integer, sequence, witness and chain-id list **exactly**,
and the fixed-point fields inside their own step. **A blanket tolerance would pass while a chain-id
list was being corrupted**, which is the failure mode the strictness exists for.

**Three regression anchors**, in order of what they catch:

| anchor | catches |
|---|---|
| restart equals sequential — read from an arbitrary block, get what a full read gives from there | a running difference that was not reset at a block boundary: silent, and plausible |
| an interrupted write is refused — kill a writer before `finish`, the reader must reject | a partial read being treated as a short sample |
| worker-count invariance — one sample at 1, 2, 4, 8, 16 workers is byte-identical | a block cut that depends on scheduling |

**The bench shape is the prototype's**: peak resident with N samples open and walked in lockstep,
against the 500 kB per-open-sample budget, and a walk keeping one record in a hundred against one
keeping all.
