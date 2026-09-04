//! The **psp store**: everything one sample's reads showed at every reference position a
//! run analysed, written once by the locus generator and read back by the cohort merge.
//!
//! A caller opens one file per sample and holds them all open for the whole run, so what
//! one open file costs is multiplied by the cohort size. That is the requirement the
//! format is shaped around: **an open file costs no more than 500 kB of resident memory,
//! and that figure does not grow with the block size, the read depth or the length of the
//! genome** (`doc/devel/ng/spec/psp_file_format.md` §1.1).
//!
//! The way it gets there is to untie two numbers production's `.psp` ties together. There,
//! a block is at once *how far back the compressor may look for a repeated pattern* and
//! *how much a reader must inflate before it can hand out its first record*, so shrinking a
//! block to save memory costs compression and multiplies the index. Here the block stays
//! large — for its ratio and for a small index — while a separate declared number caps the
//! compressor's look-back, and the reader decodes incrementally instead of inflating a
//! whole block.
//!
//! ```text
//! +---------------------------+
//! | header                    |  plain text: magic, length, TOML body, sentinel
//! +---------------------------+
//! | psp block 0               |  records, one compressed stream, each record
//! | psp block 1               |    opening with the head that lets a reader skip it
//! | ...                       |
//! +---------------------------+
//! | block index               |  one entry per psp block
//! +---------------------------+
//! | trailer                   |  the writer's closing payload; may be empty
//! +---------------------------+
//! | footer                    |  fixed size: offsets, counts, checksum, magic last
//! +---------------------------+
//! ```
//!
//! **Words this module uses differently from production's `src/psp/`**, and a coder moving
//! between the two will get them wrong: what production calls its *metadata section* is the
//! **trailer** here, what production calls its *trailer* is the **footer** here, and *block*
//! here always means a **psp block** — a span of reference and the records in it — never one of
//! zstd's own internal subdivisions.
//!
//! **⚠ That collision reaches the method names, and it is sharpest on one pair.** Production has
//! a `PspReader` too, and `PspReader::trailer()` means *the fixed tail* there and *the writer's
//! closing payload* here — where the fixed tail is [`PspReader::footer`]. Two types of the same
//! name with a method of the same name returning different sections of a file is the worst shape
//! this vocabulary can take, and it is kept because spec §6.2 fixes both words. `seek_to` is the
//! milder case: a genomic coordinate in production, a byte offset here.
//!
//! **This module is not a pipeline step.** It is infrastructure beside
//! [`crate::ng::ref_seq`], used by the locus generator to write and by the cohort merge to
//! read, so there is no trait and no competing implementation.
//!
//! Design authority: `doc/devel/ng/spec/psp_file_format.md` (the container),
//! `doc/devel/ng/spec/psp_record_encoding.md` (what a record holds),
//! `doc/devel/ng/spec/psp_chain_id_encoding.md` (the chain ids), and
//! `doc/devel/ng/arch/psp_file_format.md` (the code shape).

use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use crate::ng::types::GenomeRegion;

// **`pub(crate)`, with the surface re-exported below.** The split into files is how the code
// is arranged — a reader crosses the container's seams one at a time — not a second set of
// paths for callers. Left `pub`, every type here would have two public spellings, rustdoc
// would render each twice, and the first example to reach for `ng::psp::header::Header` would
// pin the path the re-export list exists to make canonical. Production does the same for two
// of its own submodules.
pub(crate) mod block;
pub(crate) mod chain_ids;
pub(crate) mod footer;
pub(crate) mod header;
pub(crate) mod index;
pub(crate) mod reader;
pub(crate) mod record;
pub(crate) mod segmentation_section;
pub(crate) mod trailer;
pub(crate) mod walk;
pub(crate) mod writer;

pub use block::{
    BlockBuilder, BlockCompressError, BlockCompressor, BlockCutRuleError, BlockHead,
    BlockHeadDecodeError, BlockReadError, BlockRecords, BlockStream, BlockWriteError,
    COMPRESSED_BLOCK_LENGTH_BYTES, CompressedBlockAt, DecodedBlockHead, READ_CHUNK_BYTES,
    ROLLING_BUFFER_CEILING_BYTES, ROLLING_BYTES, StreamedRecord, ZSTD_COMPRESSION_LEVEL,
    compressed_block_at,
};
pub use chain_ids::{LiveSet, LiveSetChanges, LiveSetReader, LiveSetWriter};
pub use footer::{
    FOOTER_BYTES, FOOTER_MAGIC, FileSection, Footer, FooterDecodeError, decode_footer,
    encode_footer,
};
pub use header::{
    ContigIdentity, DEFAULT_BLOCK_BYTE_CEILING, DEFAULT_GENOMIC_BLOCK_SIZE_BP,
    DEFAULT_LOOK_BACK_WINDOW_LOG, FIXED_INTEGER_WIDTHS_BYTES, FORMAT_VERSION, FieldEncoding,
    FieldName, FieldShape, FieldSpec, HEAD_MAGIC, HEAD_SENTINEL, HEADER_FRAMING_BYTES, Header,
    IEEE_FLOAT_WIDTHS_BYTES, MAX_HEADER_BODY_BYTES, MAX_LOOK_BACK_WINDOW_LOG,
    MIN_LOOK_BACK_WINDOW_LOG, Manifest, ParameterValue, ReadGroupIdentity, ReferenceIdentity,
    WriterProvenance,
};
pub use index::{
    BlockIndexEntry, IndexDecodeError, IndexEntryField, checksum_index, decode_index, encode_index,
};
pub use reader::{DEFAULT_LOOK_BACK_WINDOW_BUDGET_BYTES, PspReader};
pub use record::{
    BodyByteCount, DecodedRecord, DecodedRecordBody, LocatedRecord, OffsetBase, RecordDecodeError,
    RecordEncodeError, RecordEncoder, RecordHead, RecordLayout, RecordLayoutError, decode_record,
    decode_record_body, decode_the_body_of, encode_record_body, read_record_head, record_fields,
};
pub use trailer::{FileAfterAFailedReplacement, TrailerReplacementFailure, replace_trailer};
pub use walk::{RecordIter, SelectiveRecordIter};
pub use writer::{PspWriter, WriteStats, ZSTD_COMPRESSION_LEVEL_KEY};

// ---------------------------------------------------------------------
// Inspecting a file without opening it as a reader
// ---------------------------------------------------------------------

/// The header of a psp, read from the file and nothing else — no footer, no index, no block.
///
/// **It works where opening the file as a reader correctly fails**, and that is its purpose
/// (spec §6.6). A file with no footer is one a run was killed part-way through writing, and
/// every reader refuses it; but it still *has* a header, so a tool that reports what a
/// half-written file was going to be needs this, and so does anything checking that a
/// cohort's files agree on a reference before a run commits to them.
///
/// **A newer format comes back as [`PspReadError::UnsupportedVersion`], not as a parse
/// failure** — *upgrade the reader*, not *rebuild the file*. That is the whole reason the
/// header is plain text and the version is read before anything else in it.
///
/// **Nothing is allocated from a length this function has not checked.** The declared body
/// length is read from a fixed 12-byte prefix and bounded against
/// [`MAX_HEADER_BODY_BYTES`] before a buffer for it exists, so a corrupt or hostile file
/// cannot ask for memory on its own say-so.
pub fn read_header(path: &Path) -> Result<Header, PspReadError> {
    read_header_and_its_length(path).map(|(header, _)| header)
}

/// The same, and **how many bytes it occupies** — which is where the file's blocks begin.
///
/// `PspReader::open` needs that to say whether a block offset from the index lands in the blocks
/// or inside the header; without it an index entry pointing at byte 0 opened, and the refusal was
/// deferred to a corrupt-block error at read time, which is the wrong instruction arriving after
/// a cohort has committed to the sample.
pub(crate) fn read_header_and_its_length(path: &Path) -> Result<(Header, usize), PspReadError> {
    let file = File::open(path).map_err(|source| PspReadError::Io {
        path: path.to_path_buf(),
        while_doing: "opening the file",
        source,
    })?;
    read_header_from(file, path)
}

/// The same, from a handle already open on the file.
///
/// **`PspReader::open` holds one**, and opening the file a second time to read its header is a
/// second `open(2)` per sample in a cohort that opens thousands. It rewinds first, so it does not
/// depend on where the caller left the cursor.
pub(crate) fn read_header_from(
    mut file: File,
    path: &Path,
) -> Result<(Header, usize), PspReadError> {
    use std::io::Seek;
    file.rewind().map_err(|source| PspReadError::Io {
        path: path.to_path_buf(),
        while_doing: "rewinding to the header",
        source,
    })?;
    let io = |while_doing: &'static str| {
        move |source: std::io::Error| PspReadError::Io {
            path: path.to_path_buf(),
            while_doing,
            source,
        }
    };

    // The magic and the declared length are a fixed-size prefix. Reading exactly those first
    // is what keeps the length from sizing a buffer before it has been believed.
    let mut framing = [0u8; HEAD_MAGIC.len() + 8];
    file.read_exact(&mut framing)
        .map_err(|source| match source.kind() {
            std::io::ErrorKind::UnexpectedEof => PspReadError::MalformedHeader {
                path: path.to_path_buf(),
                reason: "is too short to hold a header's magic and length".to_string(),
                source: None,
            },
            _ => io("reading the header's magic and length")(source),
        })?;

    if framing[..HEAD_MAGIC.len()] != HEAD_MAGIC {
        let mut found = [0u8; 4];
        found.copy_from_slice(&framing[..HEAD_MAGIC.len()]);
        return Err(PspReadError::NotAnNgPsp {
            path: path.to_path_buf(),
            found,
            expected: HEAD_MAGIC,
        });
    }

    let body_bytes = u64::from_le_bytes(
        framing[HEAD_MAGIC.len()..]
            .try_into()
            .expect("the remainder of a 12-byte array after 4 bytes is 8 bytes"),
    );
    if body_bytes == 0 || body_bytes > MAX_HEADER_BODY_BYTES {
        return Err(PspReadError::MalformedHeader {
            path: path.to_path_buf(),
            reason: format!(
                "declares a {body_bytes}-byte header body; this reader allows 1 to \
                 {MAX_HEADER_BODY_BYTES}"
            ),
            source: None,
        });
    }

    // Bounded above, so this cannot be driven past a megabyte however damaged the file is.
    let mut whole_header = framing.to_vec();
    whole_header.resize(HEADER_FRAMING_BYTES + body_bytes as usize, 0);
    file.read_exact(&mut whole_header[framing.len()..])
        .map_err(|source| match source.kind() {
            std::io::ErrorKind::UnexpectedEof => PspReadError::MalformedHeader {
                path: path.to_path_buf(),
                reason: format!("declares a {body_bytes}-byte header body and ends before it does"),
                source: None,
            },
            _ => io("reading the header body")(source),
        })?;

    Header::decode(&whole_header, path)
}

// ---------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------

/// Everything that can go wrong reading a psp.
///
/// **Every variant but one is an input problem, not a bug.** A corrupt or truncated file is
/// data a run was handed, so none of these is a panic and none of them may reach a caller as a
/// half-built record (spec §6.7). The exception is [`NoSuchBlock`](Self::NoSuchBlock), which
/// reports a caller asking for a block the file does not have, and its own doc says why that is
/// an error here rather than the panic a slice index would give.
///
/// The variants are separate because **the instruction to whoever sees them differs**:
/// rebuild the file, upgrade the reader, raise a limit, the data is damaged. Collapsing
/// them into one error loses exactly that.
///
/// **Same name as production's `psp::errors::PspReadError`, a different type in a
/// different module.** ng's surface is its own; the two must not be confused at a `use`
/// site.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum PspReadError {
    /// No valid footer, so the writer never finished: the run was interrupted and the file
    /// holds however many blocks reached disk. **Refusing it is the point** — read short,
    /// it is indistinguishable from a sample covering less of the genome.
    #[error("{} has no valid footer — the writer did not finish", path.display())]
    Incomplete { path: PathBuf },

    /// Written by a newer format than this reader understands. Upgrade the reader; the
    /// file is fine.
    ///
    /// **Reached by parsing the version and nothing else**, which is why the header stays
    /// plain text: a reader must be able to learn the version of a file it cannot
    /// otherwise read (spec §3.1).
    #[error(
        "{} is format {}.{}; this reader understands up to {}.{}",
        path.display(), found.0, found.1, supported.0, supported.1
    )]
    UnsupportedVersion {
        path: PathBuf,
        found: (u16, u16),
        supported: (u16, u16),
    },

    /// The file declares a compressor look-back window larger than this reader budgeted
    /// for, so its decoder would have to hold more than it is allowed to.
    ///
    /// **A variant of its own because the fix is a knob, not a rebuild**, and because
    /// zstd's own error here names a code rather than a number anyone can act on.
    #[error(
        "{} needs a {needed_bytes}-byte look-back window; this reader allows {allowed_bytes}",
        path.display()
    )]
    WindowTooLarge {
        path: PathBuf,
        needed_bytes: u64,
        allowed_bytes: u64,
    },

    /// The file is not an ng psp at all: its first four bytes are not this format's magic.
    ///
    /// **Its own class because the instruction is different from every other one here** —
    /// not *rebuild this file* but *you handed me the wrong file*. The everyday case is a
    /// production `.psp`, which shares the extension and nothing else; a BAM, a gzip or a
    /// text file lands here too, and `found` is what tells them apart.
    #[error(
        "{} does not begin with {}; it is not an ng psp (its first bytes are {found:02x?})",
        path.display(),
        String::from_utf8_lossy(expected.as_slice()).escape_debug()
    )]
    NotAnNgPsp {
        path: PathBuf,
        found: [u8; 4],
        expected: [u8; 4],
    },

    /// The header is an ng psp's, and this reader cannot make sense of it: the declared body
    /// length disagrees with the sentinel, the TOML does not parse, or a field in it breaks a
    /// rule the format requires — an empty contig list, an MD5 that is not 32 hex characters,
    /// a field encoding with no width.
    ///
    /// **One class rather than the dozen production distinguishes**, because they share an
    /// instruction: the file is damaged, rebuild it. `reason` carries which rule broke, and
    /// `source` carries the parser's own account of it where there is one.
    ///
    /// **Not the same thing as [`UnsupportedVersion`](Self::UnsupportedVersion)**, and the
    /// order the two are checked in is the point: a file written by a newer format must say
    /// so, not come back as unparseable TOML (spec §3.1).
    #[error("{}: {reason}", path.display())]
    MalformedHeader {
        path: PathBuf,
        reason: String,
        /// The TOML or UTF-8 error underneath, when the fault was a parser's rather than a
        /// rule's. Kept as a source rather than flattened into `reason` so that a caller
        /// walking the chain reaches the parser's own span.
        #[source]
        source: Option<Box<dyn std::error::Error + Send + Sync>>,
    },

    /// The file is damaged somewhere outside its blocks: its footer's offsets do not describe
    /// the file it sits in, its block index does not match the checksum the footer carries, or
    /// the index names a block that is not in the blocks.
    ///
    /// **One class rather than a variant per rule**, because they share an instruction — the
    /// file is damaged, rebuild it — and `reason` carries which rule broke. It is deliberately
    /// *not* [`Incomplete`](Self::Incomplete): that one means a writer was killed and the file
    /// can be rebuilt by re-running, where this one means the bytes on disk disagree with each
    /// other.
    #[error("{}: {reason}", path.display())]
    Damaged {
        path: PathBuf,
        /// Which section refused, in a sentence — *the footer does not decode*. The detail is
        /// the cause's.
        reason: String,
        /// The codec's own account underneath, when one refused. `None` when the rule that
        /// broke is one this reader checks itself, having the file's length, which no codec can
        /// (spec §6.2's three checks at `open`).
        ///
        /// **Kept as a cause rather than flattened into `reason`**, so that a caller walking
        /// the chain reaches which field of which entry the index stopped at. Four different
        /// causes reached this variant as one string until Milestone G.
        ///
        /// **Typed rather than `Box<dyn Error>`**, so that a caller deciding whether a re-run
        /// would help can tell the footer from the block index by matching rather than by
        /// downcasting or by reading `reason`.
        #[source]
        source: Option<DamageFound>,
    },

    /// A block failed to decompress, or a record ran past the end of its block. The file
    /// is damaged.
    ///
    /// **The cause is the walk's own error, not a string.** A record's
    /// [`Truncated`](RecordDecodeError::Truncated) means *the body stopped early*, which the
    /// streaming reader answers by decompressing more and only turns into damage once the block
    /// has no more to give; its [`Unsupported`](RecordDecodeError::Unsupported) means *upgrade
    /// the reader*, which is [`UnsupportedRecordEncoding`](Self::UnsupportedRecordEncoding)
    /// here and must not arrive as damage. [`BlockReadError`] is where those distinctions are
    /// already drawn, so this variant carries one rather than restating it (the C1 review, M13;
    /// the D1 review, which asked that they not be flattened).
    #[error("{}: block {block} is corrupt", path.display())]
    CorruptBlock {
        path: PathBuf,
        /// Which block, counted from the file's first. Every fault *inside* a block names the
        /// block itself, and is an ordinal into [`PspReader::block_index`].
        ///
        /// **⚠ One case is one past the last entry the index has, and so is not an index into
        /// it**: the four bytes that introduce a block are read before the block is counted, so
        /// a fault in them names a block that never began. Look it up with
        /// `block_index().get(block)`, never `block_index()[block]`.
        block: u64,
        #[source]
        source: BlockReadError,
    },

    /// One record needs more of the reader's rolling buffer than it allows a single record to
    /// hold.
    ///
    /// **A class of its own because the fix is a knob, not a rebuild** (spec §7): a genuine
    /// record can be larger than any fixed budget — §8 refuses to fix a maximum record size in
    /// the format — while a corrupt block that never parses grows the buffer until the frame
    /// runs out, and the two arrive at the same line. The refusal names the ceiling and the
    /// knob that raises it, which is what spec §7 asks of it.
    ///
    /// ⚠ **The message named the ceiling without naming a way to raise it** until the G1
    /// review pointed out that the sentence underneath said *raise the ceiling to read it*
    /// while no `PspReader` had one. [`PspReader::with_a_record_buffer_ceiling`] is that knob.
    #[error(
        "{}: a record in block {block} needs more than the {allowed_bytes} bytes this reader \
         allows one record to hold; raise it with PspReader::with_a_record_buffer_ceiling",
        path.display()
    )]
    RecordLargerThanTheReaderAllows {
        path: PathBuf,
        /// Which block, counted as [`CorruptBlock`](Self::CorruptBlock) counts it.
        block: u64,
        allowed_bytes: usize,
        #[source]
        source: BlockReadError,
    },

    /// A record, or the record layout the file's manifest declares, names something this reader
    /// does not know. **Upgrade the reader**; the file is not damaged.
    ///
    /// **Not [`UnsupportedVersion`](Self::UnsupportedVersion)**, which is the whole file's
    /// format number and is answered before a block is touched. This one is a field encoding
    /// inside a record, and it is the same instruction about a smaller subject.
    ///
    /// **The instruction is in the message and not only in this doc**, because the everyday
    /// caller prints `Display` and never walks the chain — which is where the field that
    /// disagrees is named.
    #[error(
        "{}: this reader cannot read the records in this file — upgrade the reader",
        path.display()
    )]
    UnsupportedRecordEncoding {
        path: PathBuf,
        #[source]
        source: BlockReadError,
    },

    /// A record-buffer ceiling at or under the buffer it is a ceiling on, which would turn away
    /// records the buffer already holds without growing at all.
    ///
    /// **The second of the two variants here that report the caller's mistake rather than the
    /// file's** — see [`NoSuchBlock`](Self::NoSuchBlock) — and it is refused at
    /// [`PspReader::with_a_record_buffer_ceiling`] rather than at the record that would have
    /// tripped it, so a setting that could never work fails where it was made.
    #[error(
        "{}: a record-buffer ceiling of {ceiling} bytes, which is not above the \
         {buffer_bytes}-byte buffer it is a ceiling on",
        path.display()
    )]
    RecordBufferCeilingTooSmall {
        path: PathBuf,
        ceiling: usize,
        buffer_bytes: usize,
    },

    /// A block was asked for by an ordinal the file does not have.
    ///
    /// **One of the two variants here that report the caller's mistake rather than the file's**
    /// — the other is [`RecordBufferCeilingTooSmall`](Self::RecordBufferCeilingTooSmall) — and
    /// it is an error rather than a panic because the caller is a cohort holding thousands of
    /// these open: an ordinal computed wrongly must fail one sample, not abort the run.
    /// [`PspReader::records_from`] cannot raise it — it derives the ordinal from the index —
    /// so it reaches only a caller that indexed [`PspReader::block_index`] itself.
    #[error(
        "{} has {blocks_in_the_file} blocks; block {ordinal_asked_for} was asked for",
        path.display()
    )]
    NoSuchBlock {
        path: PathBuf,
        /// **An ordinal, not a coordinate**: the index into [`PspReader::block_index`].
        ordinal_asked_for: u64,
        /// A count, not the blocks themselves.
        blocks_in_the_file: u64,
    },

    /// Reading the file's bytes failed.
    ///
    /// **It names the file and what was being done to it**, because a bare
    /// `No such file or directory (os error 2)` from a cohort opening a thousand samples
    /// says nothing anyone can act on. Not `#[from]`: every site that raises this knows the
    /// path, and a blanket conversion is how that gets lost.
    #[error("{} could not be read while {while_doing}", path.display())]
    Io {
        path: PathBuf,
        /// What the reader was doing, as a phrase that follows "while": `reading the
        /// header`, `seeking to block 12`.
        while_doing: &'static str,
        #[source]
        source: std::io::Error,
    },
}

/// Which of a psp's own decoders refused, under a [`PspReadError::Damaged`].
///
/// **Two inhabitants, and they are the two sections a reader decodes before it touches a
/// block**: the fixed tail and the block index. A caller deciding whether re-running the
/// pileup would help can tell them apart by matching on this rather than by downcasting a
/// `Box<dyn Error>` or by reading a sentence.
///
/// **Transparent**, so a chain reads *the footer does not decode* from the parent's `reason`
/// and then the decoder's own account — rather than the same sentence twice.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum DamageFound {
    /// The fixed 48-byte tail did not decode.
    #[error(transparent)]
    Footer(#[from] FooterDecodeError),
    /// The block index did not decode, though its bytes matched the footer's checksum.
    #[error(transparent)]
    BlockIndex(#[from] IndexDecodeError),
}

/// Which of a psp's own decoders refused a structure the **writer** read back before believing
/// it, under [`PspWriteError::WouldNotBeReadable`].
///
/// **`finish` decodes the index and the footer it is about to write with the very functions a
/// reader will use** (spec §6.3), so this is the read side's own account arriving on the write
/// side. Transparent, so a chain reads *which structure* from the parent and *what was wrong with
/// it* from here.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum NotReadable {
    /// The block index the writer built does not decode.
    #[error(transparent)]
    BlockIndex(#[from] IndexDecodeError),
    /// The footer the writer built does not decode.
    #[error(transparent)]
    Footer(#[from] FooterDecodeError),
    /// A block's opening fields do not decode out of the block the writer just built.
    #[error(transparent)]
    BlockHead(#[from] BlockHeadDecodeError),
}

/// Why a writer cannot honour what a file's header declares, under
/// [`PspWriteError::UnsupportedHeader`].
///
/// **Two things are built from the header before a byte is written** — the rule that cuts a block
/// and the compressor — and they fail for different reasons, which a caller reading only a
/// sentence could not tell apart.
///
/// **They also come from two different halves of the header**, which is why neither this type nor
/// the error above says *manifest* any more: the cut rule comes from the manifest, and the
/// compression level from the writer's own recorded parameters. *This type was `ManifestRefusal`
/// and the error `UnsupportedManifest` until 2026-08-30, when the second of those was found being
/// reported as a fault in the first.*
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum HeaderRefusal {
    /// The genomic block size or the byte ceiling is one no cut rule can be built from.
    #[error(transparent)]
    CutRule(#[from] BlockCutRuleError),
    /// The header records a compression level in a shape this writer cannot read — a string, a
    /// float, a boolean.
    ///
    /// **A refusal rather than a fallback**, and the distinction is the point: *absent* is what a
    /// file written before the parameter existed looks like and legitimately falls back to this
    /// build's level, where *present and unreadable* is **a setting recorded and then ignored**,
    /// which is worse than one not recorded (spec goal 4). Until the G4 review they were the same
    /// arm, and a file recording `"1"` was appended to at level 9 with nothing said.
    #[error("the header records a compression level of {recorded}, which is not an integer")]
    UnreadableLevel { recorded: String },

    /// The header records a compression level no `i32` can hold.
    ///
    /// **Its own variant so that the refusal names the level the file records.** It used to
    /// saturate to `i32::MAX` before the range check saw it, so a file recording −5,000,000,000
    /// was refused for recording 2,147,483,647 — a number that is not in the file, which is the
    /// clamp `BlockCompressor::with_level` exists to prevent, moved one function upstream.
    #[error(
        "the header records a compression level of {recorded}, which is past what a level can be"
    )]
    LevelPastAnyLevel { recorded: i64 },

    /// The look-back window or the compression level is one no compressor can be built from.
    ///
    /// **Reached by `append` on a file whose recorded compression level zstd will not take** —
    /// 99, say. It is *not* reachable through `create`, whose header is encoded first and whose
    /// window log `header.rs`'s own validation refuses as
    /// [`PspWriteError::InvalidHeaderField`] before the compressor sees it.
    ///
    /// ⚠ **The G3 review found it unreachable and this doc said so; G4 made it reachable and the
    /// doc was not updated until G4's own review caught that.** A variant documented as
    /// unreachable, quietly reachable, is a path nothing protects.
    #[error(transparent)]
    Compressor(#[from] BlockCompressError),
}

impl PspReadError {
    /// A rule this reader checks itself broke, and there is no error underneath it.
    ///
    /// **The rules with no cause are the ones that need the file's length** — that the sections
    /// end where the footer begins, that the index starts after the header, that every block
    /// offset lands in the blocks. No codec in the file knows how long the file is, so none of
    /// them can have refused first.
    pub(crate) fn damaged(path: &Path, reason: String) -> Self {
        Self::Damaged {
            path: path.to_path_buf(),
            reason,
            source: None,
        }
    }

    /// A decoder underneath refused, and **its own account is kept as the cause** rather than
    /// pressed into `reason`.
    ///
    /// `reason` names the section; the cause says what was wrong with it. **Neither repeats the
    /// other**: a chain that prints the same sentence twice reads as a fault in the printer.
    pub(crate) fn damaged_by(path: &Path, reason: &str, source: DamageFound) -> Self {
        Self::Damaged {
            path: path.to_path_buf(),
            reason: reason.to_string(),
            source: Some(source),
        }
    }
}

/// Everything that can go wrong writing a psp.
///
/// **Same name as production's `psp::errors::PspWriteError`, a different type in a
/// different module** — see [`PspReadError`].
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum PspWriteError {
    /// A record was pushed that starts before the one before it. **Refused rather than
    /// accepted**: coordinate order is what the block index and every seek rest on, and a
    /// file that breaks it seeks wrongly rather than failing (spec §6.3).
    #[error(
        "{}: record at {offered} starts before the previous record at {previous}",
        path.display()
    )]
    OutOfOrder {
        /// Which file was being written. **The one write error raised per record**, so in a
        /// cohort gathering sixty samples at once it is the only thing that says which.
        path: PathBuf,
        previous: GenomeRegion,
        offered: GenomeRegion,
    },

    /// A header field the writer was handed cannot be written: an empty contig list, a
    /// non-finite parameter value, a version this writer does not produce.
    #[error("header field {field} is not writable: {reason}")]
    InvalidHeaderField { field: String, reason: String },

    /// An append was asked to extend a file whose header declares something this writer cannot
    /// honour. **Appending does not rewrite the header**, so the added records must use the
    /// encodings the file already declares (spec §6.4).
    ///
    /// **The sentence says *header* rather than *manifest*, and the difference is where an
    /// operator looks.** The cut rule is built from the manifest, but the compression level is
    /// read from the writer's recorded parameters — the header's provenance half — so a level
    /// this build cannot read was arriving under a sentence pointing at the wrong section.
    #[error("{} declares something in its header this writer cannot honour", path.display())]
    UnsupportedHeader {
        path: PathBuf,
        /// **Typed, so a caller can tell the cut rule from the compressor** without reading a
        /// sentence — and so `Display` and the chain do not print the same words twice.
        #[source]
        source: HeaderRefusal,
    },

    /// Reopening a finished file failed. **Append and trailer replacement read before they
    /// write** — the footer for the offsets, the header for the manifest — so both surface
    /// the read-side classes (spec §6.7) rather than restating them here.
    #[error("{}: the file could not be reopened", path.display())]
    Reopen {
        path: PathBuf,
        #[source]
        source: PspReadError,
    },

    /// A record the block builder or the record codec refused: an empty region, a region at
    /// the coordinate ceiling, a body longer than a head can describe, or a contig already
    /// written and revisited.
    ///
    /// **Not [`OutOfOrder`](Self::OutOfOrder)**, which is the one refusal a caller can act on by
    /// reordering its input; these say the record itself cannot be written.
    #[error("{}: a record could not be written", path.display())]
    RecordRefused {
        path: PathBuf,
        #[source]
        source: BlockWriteError,
    },

    /// A block the compressor refused — its own configuration, or a frame longer than the
    /// four-byte length in front of it can describe.
    #[error("{}: a block could not be compressed", path.display())]
    BlockRefused {
        path: PathBuf,
        #[source]
        source: BlockCompressError,
    },

    /// **The writer was about to produce a file its own reader would refuse**, and stopped.
    ///
    /// Raised by `finish`, which decodes the index and the footer it is about to write with the
    /// very functions that will read them back. A file that reaches disk and is then rejected is
    /// a walk thrown away; this turns that into an error before the bytes land.
    #[error("{}: this writer would have produced a file it cannot read — {reason}", path.display())]
    WouldNotBeReadable {
        path: PathBuf,
        /// Which structure would not read back, in a sentence. The detail is the cause's.
        reason: String,
        /// The decoder's own account, when a decoder refused. `None` when the rule that broke is
        /// one the writer checks itself — the index's checksum against the footer's copy of it,
        /// which no decoder can see.
        #[source]
        source: Option<NotReadable>,
    },

    /// Writing the file's bytes failed. Named and sourced the way
    /// [`PspReadError::Io`] is, and for the same reason.
    #[error("{} could not be written while {while_doing}", path.display())]
    Io {
        path: PathBuf,
        /// What the writer was doing, as a phrase that follows "while": `writing the
        /// header`, `flushing block 12`, `syncing`.
        while_doing: &'static str,
        #[source]
        source: std::io::Error,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position};

    /// **A half-written file still has a header, at every cut of a whole psp** — spec §6.6's
    /// claim, which nothing held before this step.
    ///
    /// Its neighbour [`a_file_cut_short_inside_its_header_is_refused_at_every_cut`] sweeps a bare
    /// encoded header and only requires *some* refusal; this sweeps a finished file and requires
    /// the class and the boundary.
    ///
    /// [`read_header`] exists so that a tool can report what a killed run *was going to be*, on a
    /// file every reader correctly refuses. That is only worth anything if it holds wherever the
    /// kill landed, so this cuts a finished psp at every byte and requires exactly two outcomes,
    /// with the boundary between them at the header's own length:
    ///
    /// - **at or past the header, it reads** — and gives back the same header the whole file
    ///   does, not merely *a* header;
    /// - **below it, it refuses without panicking**, and always as
    ///   [`MalformedHeader`](PspReadError::MalformedHeader).
    ///
    /// ⚠ **Never as [`NotAnNgPsp`](PspReadError::NotAnNgPsp).** The magic is compared only
    /// *after* a twelve-byte read of the magic and the declared body
    /// length has succeeded, so a file cut inside its own magic is refused for being too short
    /// before its first four bytes are ever looked at. **A truncated ng psp is therefore never
    /// reported as the wrong kind of file** — which is the right answer, since it is an ng psp,
    /// and *you handed me the wrong file* would send its owner looking for a file that is not
    /// missing. This test asserted the opposite when it was written, and the sweep is what said
    /// so: 0 of 3,136 cuts below the header gave `NotAnNgPsp`.
    ///
    /// ⚠ **The boundary is asserted, not assumed.** A `read_header` that returned `Ok` on a file
    /// cut one byte short of its header would pass a test that only checked *some* cut reads.
    #[test]
    fn a_half_written_file_still_has_a_header_at_every_cut() {
        use crate::ng::psp::writer::tests_support::{a_finished_psp, bytes_of, rewrite};

        let (_dir, path) = a_finished_psp();
        let whole = bytes_of(&path);
        let (want, header_bytes) =
            read_header_and_its_length(&path).expect("the whole file's header reads");

        assert!(
            header_bytes > HEAD_MAGIC.len(),
            "the fixture's header has to be longer than its magic, or the cuts inside the magic \
             — the ones that would tempt a `NotAnNgPsp` — are never taken"
        );

        let mut cuts_that_read = 0usize;
        let mut cuts_refused = 0usize;
        for cut in 0..whole.len() {
            rewrite(&path, &whole[..cut]);
            match read_header(&path) {
                Ok(found) => {
                    assert!(
                        cut >= header_bytes,
                        "a cut at {cut} is inside a {header_bytes}-byte header and read anyway"
                    );
                    assert_eq!(
                        found, want,
                        "a cut at {cut} gave a different header from the whole file's"
                    );
                    cuts_that_read += 1;
                }
                Err(PspReadError::MalformedHeader { .. }) => {
                    assert!(
                        cut < header_bytes,
                        "a cut at {cut} holds a whole {header_bytes}-byte header and was refused"
                    );
                    cuts_refused += 1;
                }
                Err(other) => panic!("a cut at {cut} gave {other}"),
            }
        }
        assert_eq!(
            cuts_that_read,
            whole.len() - header_bytes,
            "every cut at or past the header must read, and no other"
        );
        assert_eq!(
            cuts_refused, header_bytes,
            "and every cut below it must be refused as a malformed header"
        );
    }

    /// The messages are the contract spec §6.7 tabulates, so they are pinned rather than
    /// left to whatever `thiserror` last rendered. Each has to say what happened *and*
    /// carry the number whoever sees it must act on.
    #[test]
    fn read_errors_name_the_file_and_the_limit_that_was_exceeded() {
        let incomplete = PspReadError::Incomplete {
            path: PathBuf::from("SRR7279481.psp"),
        };
        assert_eq!(
            incomplete.to_string(),
            "SRR7279481.psp has no valid footer — the writer did not finish"
        );

        let version = PspReadError::UnsupportedVersion {
            path: PathBuf::from("SRR7279481.psp"),
            found: (2, 0),
            supported: (1, 0),
        };
        assert_eq!(
            version.to_string(),
            "SRR7279481.psp is format 2.0; this reader understands up to 1.0"
        );

        let window = PspReadError::WindowTooLarge {
            path: PathBuf::from("SRR7279481.psp"),
            needed_bytes: 524_288,
            allowed_bytes: 32_768,
        };
        assert_eq!(
            window.to_string(),
            "SRR7279481.psp needs a 524288-byte look-back window; this reader allows 32768"
        );
    }

    /// An out-of-order push names the file and both regions. It is the one write error
    /// raised per record, so in a cohort gathering sixty samples at once the file's name is
    /// the only thing that says which of them stopped.
    #[test]
    fn an_out_of_order_push_names_the_file_and_both_regions() {
        let at = |start: u64, end: u64| GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        };
        let refused = PspWriteError::OutOfOrder {
            path: PathBuf::from("SRR7279481.psp"),
            previous: at(1_000, 1_000),
            offered: at(900, 900),
        };
        assert_eq!(
            refused.to_string(),
            "SRR7279481.psp: record at contig 0:900-900 starts before the previous record at \
             contig 0:1000-1000"
        );
    }

    /// Reopening for an append or a trailer rewrite reads the footer first, so a file with no
    /// footer surfaces as the read-side class. **The wrapper says a write was in progress and
    /// the chain still reaches the read error** — a transparent wrapper would have said
    /// neither.
    #[test]
    fn a_reopen_says_a_write_was_in_progress_and_still_carries_its_cause() {
        let refused = PspWriteError::Reopen {
            path: PathBuf::from("half-written.psp"),
            source: PspReadError::Incomplete {
                path: PathBuf::from("half-written.psp"),
            },
        };
        assert_eq!(
            refused.to_string(),
            "half-written.psp: the file could not be reopened"
        );

        let cause = std::error::Error::source(&refused).expect("the read error is the cause");
        assert_eq!(
            cause.to_string(),
            "half-written.psp has no valid footer — the writer did not finish"
        );
    }

    /// Handing this reader a production `.psp` — the same extension, a different format — is
    /// not the same instruction as handing it a damaged ng one, so it is not the same error.
    /// The message has to carry what the bytes actually were, or a BAM and a gzip are
    /// indistinguishable in a log.
    #[test]
    fn a_foreign_file_is_its_own_class_and_the_message_names_the_bytes_found() {
        let refused = PspReadError::NotAnNgPsp {
            path: PathBuf::from("SRR7279481.psp"),
            found: *b"PSP\n",
            expected: header::HEAD_MAGIC,
        };
        let said = refused.to_string();
        assert!(said.contains("it is not an ng psp"), "got {said}");
        assert!(
            said.contains("50, 53, 50, 0a"),
            "the message must show the bytes that were there; got {said}"
        );
    }

    /// An I/O failure from a cohort opening a thousand samples has to say which file and what
    /// was being done to it. `No such file or directory (os error 2)` on its own does not.
    #[test]
    fn an_io_failure_names_the_file_and_what_was_being_done_to_it() {
        let refused = PspReadError::Io {
            path: PathBuf::from("SRR7279481.psp"),
            while_doing: "reading the header",
            source: std::io::Error::from(std::io::ErrorKind::NotFound),
        };
        assert_eq!(
            refused.to_string(),
            "SRR7279481.psp could not be read while reading the header"
        );
        assert!(std::error::Error::source(&refused).is_some());
    }

    // -----------------------------------------------------------------
    // A3 — read_header, from a file on disk
    // -----------------------------------------------------------------

    /// The same tomato-shaped header the encoder's own tests use, trimmed to what a file
    /// needs.
    fn a_header() -> Header {
        let contigs = vec![ContigIdentity {
            name: "SL4.0ch01".to_string(),
            length: 90_863_682,
            md5: Some([0x1b; 16]),
        }];
        Header {
            format_version: FORMAT_VERSION,
            sample: "SRR7279481".to_string(),
            reference: ReferenceIdentity {
                name: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                md5: Some([0x0a; 16]),
            },
            segmentation_inputs:
                crate::ng::psp::segmentation_section::segmentation_inputs_for_tests(&contigs),
            contigs,
            read_groups: crate::ng::psp::header::read_groups_for_tests(),
            observation_reach_ceiling_bp: crate::ng::types::Bp(4_000),
            writer: WriterProvenance {
                tool: "ng".to_string(),
                version: "0.1.0".to_string(),
                subcommand: "pileup".to_string(),
                input_alignments: vec!["SRR7279481.cram".to_string()],
                input_reference: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                command_line: "ng pileup --sample SRR7279481".to_string(),
                parameters: std::collections::BTreeMap::from([(
                    "depth-cap".to_string(),
                    ParameterValue::Integer(300),
                )]),
                created: "2026-08-26T09:15:00Z".parse().expect("an RFC 3339 stamp"),
            },
            manifest: Manifest {
                genomic_block_size_bp: DEFAULT_GENOMIC_BLOCK_SIZE_BP,
                block_byte_ceiling: None,
                look_back_window_log: DEFAULT_LOOK_BACK_WINDOW_LOG,
                fields: vec![FieldSpec {
                    name: FieldName("position-offset".to_string()),
                    encoding: FieldEncoding::Varint,
                }],
            },
        }
    }

    /// Put `bytes` in a file and hand back the directory that owns it, because dropping the
    /// directory deletes the file.
    fn a_file_holding(bytes: &[u8]) -> (tempfile::TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = dir.path().join("SRR7279481.psp");
        std::fs::write(&path, bytes).expect("the fixture is written");
        (dir, path)
    }

    /// **The case this function exists for.** A run killed part-way leaves a file with a
    /// header, some blocks and no footer; every reader refuses it, and a tool reporting what
    /// it was going to be still has to read the header. So the fixture is a header followed
    /// by bytes that are not a footer.
    #[test]
    fn a_header_reads_from_a_file_that_has_no_footer() {
        let written = a_header();
        let mut on_disk = written.encode().expect("a valid header encodes");
        on_disk.extend_from_slice(&[0x5a; 8_192]);
        let (_dir, path) = a_file_holding(&on_disk);

        assert_eq!(read_header(&path).expect("the header reads"), written);
    }

    /// And it reads a file that holds nothing but a header, which is what a writer that was
    /// killed before its first block leaves behind.
    #[test]
    fn a_header_reads_from_a_file_that_holds_nothing_else() {
        let written = a_header();
        let (_dir, path) = a_file_holding(&written.encode().expect("encodes"));
        assert_eq!(read_header(&path).expect("the header reads"), written);
    }

    /// **A newer format is `UnsupportedVersion`, not a parse failure** — spec §6.7's
    /// distinction, and the reason the header is plain text. The body here carries a key no
    /// 1.0 reader knows, which is what a real 2.0 file would do.
    #[test]
    fn a_file_from_a_newer_format_says_so_rather_than_failing_to_parse() {
        let body = "format-version = \"2.0\"\n\
                    whatever-version-two-added = { shape = \"nothing here knows\" }\n";
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&(body.len() as u64).to_le_bytes());
        bytes.extend_from_slice(body.as_bytes());
        bytes.extend_from_slice(HEAD_SENTINEL);
        let (_dir, path) = a_file_holding(&bytes);

        match read_header(&path) {
            Err(PspReadError::UnsupportedVersion {
                found, supported, ..
            }) => {
                assert_eq!(found, (2, 0));
                assert_eq!(supported, FORMAT_VERSION);
            }
            other => panic!("expected UnsupportedVersion, got {other:?}"),
        }
    }

    /// Both formats use the extension `.psp`. Handed production's, this says which it was
    /// given rather than reporting a damaged ng file.
    #[test]
    fn a_production_psp_on_disk_is_refused_as_a_foreign_file() {
        let mut productions = Vec::new();
        productions.extend_from_slice(&crate::psp::header::HEAD_MAGIC);
        productions.extend_from_slice(&(1u64).to_le_bytes());
        productions.push(b'x');
        let (_dir, path) = a_file_holding(&productions);

        match read_header(&path) {
            Err(PspReadError::NotAnNgPsp { found, .. }) => {
                assert_eq!(found, crate::psp::header::HEAD_MAGIC);
            }
            other => panic!("expected NotAnNgPsp, got {other:?}"),
        }
    }

    /// A file cut short at any point inside its header is refused, at every cut. **The
    /// interesting half is the twelve-byte prefix**: below it there is not even a length to
    /// believe, and above it the file promises a body it does not have.
    #[test]
    fn a_file_cut_short_inside_its_header_is_refused_at_every_cut() {
        let whole = a_header().encode().expect("encodes");
        for cut in 0..whole.len() {
            let (_dir, path) = a_file_holding(&whole[..cut]);
            let refused = read_header(&path);
            assert!(
                refused.is_err(),
                "a file of {cut} bytes must not yield a header"
            );
        }
        let (_dir, path) = a_file_holding(&whole);
        assert!(read_header(&path).is_ok(), "and the whole file still reads");
    }

    /// The declared length is bounded **before a buffer for it exists**, so a corrupt file
    /// cannot ask this reader for four exabytes of memory. Each fixture is twelve bytes long:
    /// if the bound were checked after the allocation, this test would not return.
    ///
    /// **Both ends of the bound, and a zero is one of them.** Dropping the zero check alone
    /// survived every other test in this module — a zero-length body is caught further in,
    /// where the file is refused for ending early rather than for the length it declared, and
    /// those are different things to tell whoever is looking at the file.
    #[test]
    fn a_declared_length_outside_the_readers_bound_is_refused_without_allocating_for_it() {
        for declared in [0u64, MAX_HEADER_BODY_BYTES + 1, u64::MAX] {
            let mut bytes = Vec::new();
            bytes.extend_from_slice(&HEAD_MAGIC);
            bytes.extend_from_slice(&declared.to_le_bytes());
            let (_dir, path) = a_file_holding(&bytes);

            let refused = match read_header(&path) {
                Err(refused) => refused,
                Ok(_) => panic!("a {declared}-byte body must be refused"),
            };
            assert!(
                refused
                    .to_string()
                    .contains(&format!("this reader allows 1 to {MAX_HEADER_BODY_BYTES}")),
                "for a declared length of {declared}, got {refused}"
            );
        }
    }

    /// A file that is not there is an I/O failure that names the file and what was being
    /// done, not a bare operating-system message.
    #[test]
    fn a_missing_file_names_itself_and_the_operation() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let path = dir.path().join("never-written.psp");

        match read_header(&path) {
            Err(PspReadError::Io {
                while_doing,
                source,
                ..
            }) => {
                assert_eq!(while_doing, "opening the file");
                assert_eq!(source.kind(), std::io::ErrorKind::NotFound);
            }
            other => panic!("expected Io, got {other:?}"),
        }
    }

    /// An empty file has no magic to check, so it is refused for being too short rather than
    /// for being some other format — there are no bytes to say which.
    #[test]
    fn an_empty_file_is_refused_for_being_too_short() {
        let (_dir, path) = a_file_holding(&[]);
        let refused = read_header(&path).expect_err("an empty file holds no header");
        assert!(refused.to_string().contains("too short"), "got {refused}");
    }
}
