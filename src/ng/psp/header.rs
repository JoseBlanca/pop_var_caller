//! The header: a file's own account of what it holds and how it was written.
//!
//! It sits at the start of every psp and **stays plain text**, framed exactly as
//! production's `.psp` frames its own (`src/psp/header.rs`): a 4-byte magic, an 8-byte
//! little-endian length for the TOML body, the body, and a sentinel line. The length is
//! authoritative and the sentinel is a cross-check.
//!
//! **Why plain text is worth its bytes**: `head` and a TOML parser tell you what a file is
//! without a special tool, and the format version in particular *must* be readable without
//! knowing the version — a binary header would create a chicken-and-egg problem the first
//! time its layout changed (spec §3.1). It is written once per file rather than once per
//! record, so the cost does not scale with anything.
//!
//! The header carries what is known **before** any record is written. What is only known
//! afterwards — the per-sample summary, the coverage-against-GC histogram — goes in the
//! trailer instead (spec §3.4).

use std::collections::BTreeMap;
use std::path::Path;

use serde::{Deserialize, Serialize};

use super::{PspReadError, PspWriteError};
use crate::ng::types::Bp;

/// A file's own account of how it was written: everything a reader needs before it touches
/// a block, and the only part of the file that is plain text.
///
/// **Fixed at `create` and never rewritten.** An append reuses the header it finds, so the
/// appended records must use the encodings already declared (spec §6.3, §6.4).
#[derive(Debug, Clone, PartialEq)]
pub struct Header {
    /// `MAJOR.MINOR`. **Parsed before anything else and never behind a binary encoding** —
    /// a reader must be able to learn the version of a file it cannot otherwise read.
    pub format_version: (u16, u16),
    /// The sample whose reads these records describe: one psp, one sample.
    pub sample: String,
    /// Which reference the reads were called against, in the form two files can be
    /// compared on before a cohort run commits to them.
    pub reference: ReferenceIdentity,
    /// Every contig, **in the order that defines [`ContigId`]** — `ContigId(i)` is
    /// `contigs[i]`, the `@SQ` / `.fai` order the reference was read in. This is the
    /// coordinate space every record in the file is written against.
    ///
    /// [`ContigId`]: crate::ng::types::ContigId
    pub contigs: Vec<ContigIdentity>,
    /// What produced the file and what it ran with.
    pub writer: WriterProvenance,
    /// How this file encodes what it holds.
    pub manifest: Manifest,
}

/// Which reference a psp was called against: a name to print and a digest to compare.
///
/// **The pair is what a cohort-agreement check needs** — anything that verifies a run's
/// samples were called against the same assembly before committing to them (spec §6.6),
/// which is one of the two reasons `read_header` exists.
///
/// **Deliberately not [`crate::ng::reference_info::ReferenceInfo`]**, which is the whole
/// description of a reference on disk: it carries each contig's `.fai` geometry — the byte
/// offset of its first base, its line width — and the absolute path of the FASTA. Neither
/// survives a round trip through a psp honestly (a parsed header would have to invent
/// geometry it never stored), and a file that recorded the producer's directory layout
/// would leak it to everyone the file is shared with.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReferenceIdentity {
    /// The reference FASTA's **basename**, no directory component — what to print when
    /// telling a user which assembly a file belongs to.
    pub name: String,
    /// MD5 of every contig's uppercased bases concatenated in file order — the
    /// whole-assembly digest [`ReferenceInfo::md5`] carries, **not** a SAM `@SQ M5`.
    ///
    /// `None` when the reference was read from a `.fai` alone, which holds no bases. An
    /// absent digest is a check that cannot be made, not a check that passed.
    ///
    /// [`ReferenceInfo::md5`]: crate::ng::reference_info::ReferenceInfo::md5
    pub md5: Option<[u8; 16]>,
}

/// One contig as the header records it: its name, its length, and its own MD5.
///
/// **The same name as [`crate::fasta::ContigEntry`] and a different type**, because that
/// one compares an absent MD5 as a wildcard — `Some(a) == None` is `true` there, so a
/// round-trip test that silently dropped every digest would still pass. Equality here is
/// plain field equality, which is what a test of the header's encoding needs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContigIdentity {
    /// The contig's name: the first word of the FASTA `>` header, the same string as a SAM
    /// `@SQ SN`.
    pub name: String,
    /// Total bases — the same number as `@SQ LN`.
    pub length: u64,
    /// MD5 of this contig's uppercased bases: the SAM `@SQ M5`. `None` from a `.fai`-only
    /// read, which holds no bases.
    pub md5: Option<[u8; 16]>,
}

/// What produced the file, in a form that reproduces a run on any host without recording
/// the producer's directory layout or username.
#[derive(Debug, Clone, PartialEq)]
pub struct WriterProvenance {
    /// The program: `ng`.
    pub tool: String,
    /// Its version string.
    pub version: String,
    /// Which subcommand of it wrote the file.
    pub subcommand: String,
    /// The alignment files the records were gathered from, **basenames only** — no
    /// directory component.
    pub input_alignments: Vec<String>,
    /// The reference FASTA, **basename only**.
    pub input_reference: String,
    /// The full invoking command line, space-joined, so the file records exactly how it
    /// was produced.
    pub command_line: String,
    /// One entry per exposed knob. A `BTreeMap` because the order has to be the same on
    /// every run: goal 5 is that the same sample gathered at any worker count gives the
    /// same bytes, and a map with a run-dependent order would break it in the header.
    pub parameters: BTreeMap<String, ParameterValue>,
    /// When the writer ran.
    ///
    /// **The one field that legitimately differs between two otherwise identical runs**,
    /// which is what the worker-count invariance check has to exempt (spec §7).
    pub created: toml::value::Datetime,
}

/// One parameter's value, in the four shapes TOML has for a scalar.
///
/// `untagged`, so a parameter appears in the header as the value itself — `depth-cap = 300`,
/// not `depth-cap = { integer = 300 }`. Reading the header has to be worth doing by eye.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(untagged)]
pub enum ParameterValue {
    /// `i64` to match TOML's signed integer width.
    Integer(i64),
    /// Must be finite — a NaN or an infinity has no TOML spelling to round-trip through.
    Float(f64),
    Boolean(bool),
    String(String),
}

/// How this file encodes what it holds. **Every value here is the writer's choice,
/// recorded so a reader is driven by the file rather than by an assumption** — which is
/// goal 4, and the reason none of these is a constant in the code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Manifest {
    /// Reference bases per psp block. **A grid on the coordinate, not a running total**: a
    /// block ends when a position crosses into the next multiple of this, so every sample
    /// cuts at the same coordinates and a cohort reader stepping across a region touches
    /// one aligned block per sample. Default 100 kb (spec §4.1).
    pub genomic_block_size_bp: Bp,
    /// A secondary cut for when one span holds too much: close a block early once it
    /// exceeds this many bytes. `None` = no ceiling.
    ///
    /// It exists because a span-cut block has a variable size in bytes — at 300 reads a
    /// position a 5 kb span is a great deal of data, and on a sparse sample it is a handful
    /// of records that compress badly because the compressor starts cold (spec §4.1).
    pub block_byte_ceiling: Option<u32>,
    /// The compressor's look-back window, as the exponent of two zstd takes.
    ///
    /// **A reader configures its decoder from this**; assuming a value makes it reject
    /// legitimate files with an error that names a zstd code rather than a number anyone
    /// can act on (spec §4.2).
    pub look_back_window_log: u8,
    /// One entry per field of a record, **in encoding order**.
    pub fields: Vec<FieldSpec>,
}

impl Manifest {
    /// What a writer that was told nothing records in the file it writes.
    ///
    /// **The one place every default is assembled**, so that a caller does not have to know
    /// which four values a manifest holds and F3's writer is not the fifth hand-built literal
    /// to get one of them wrong: [`DEFAULT_GENOMIC_BLOCK_SIZE_BP`],
    /// [`DEFAULT_BLOCK_BYTE_CEILING`], [`DEFAULT_LOOK_BACK_WINDOW_LOG`], and the fields this
    /// build of the record codec encodes.
    ///
    /// **Every one of them still goes into the file.** A reader is driven by what it finds
    /// there and never by these (spec §4.5); they are only what a writer starts from.
    pub fn as_this_build_writes_it() -> Self {
        Self {
            genomic_block_size_bp: DEFAULT_GENOMIC_BLOCK_SIZE_BP,
            block_byte_ceiling: DEFAULT_BLOCK_BYTE_CEILING,
            look_back_window_log: DEFAULT_LOOK_BACK_WINDOW_LOG,
            fields: crate::ng::psp::record::record_fields(),
        }
    }
}

/// A field of a record, and how to read it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldSpec {
    pub name: FieldName,
    pub encoding: FieldEncoding,
}

/// A record field's name, as the file spells it.
///
/// **A string rather than a closed set of known fields**, and deliberately: a reader has
/// to be able to carry a name it does not recognise, because that is what lets it skip an
/// unfamiliar field instead of misparsing the record around it (spec §4.5). The newtype is
/// here so a field's name cannot be passed where a sample's or a contig's is wanted.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct FieldName(pub String);

/// How one field's values are laid down.
///
/// **Closed deliberately.** An open-ended scheme would buy flexibility nobody has asked
/// for and cost speed on a path that decodes about twenty million records a second, where
/// every field goes through it (spec §4.5). A closed set keeps the reader a `match` rather
/// than a plugin host.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldEncoding {
    /// LEB128: a value needing one byte takes one and a value needing three takes three,
    /// in the same file.
    Varint,
    /// Zig-zag LEB128, for values that go negative — a difference from the previous
    /// record, typically.
    SignedVarint,
    /// A fixed-width little-endian integer, `width_bytes` wide.
    FixedWidthInteger { width_bytes: u8 },
    /// Raw IEEE bytes. **The escape hatch, not the default** — a float stored raw is a
    /// float two callers can disagree about in its last bits, and the three quantities
    /// that arrive as floating point are 70 % of the compressed file when stored this way
    /// (spec psp_record_encoding.md §5).
    IeeeFloat { width_bytes: u8 },
    /// A count of steps of `1 / steps_per_unit`, written as a **zig-zag** variable-length
    /// integer. The name says which way round the arithmetic goes: 4,096 steps to one
    /// natural log, so the stored integer is the value multiplied by it.
    ///
    /// Signed, because the one field carrying this today is a sum of log error
    /// probabilities and every probability's log is at most zero. Zig-zag frames its bytes
    /// exactly as the unsigned form does, so a reader walking past a field of this kind that
    /// it does not recognise measures it without knowing the sign convention.
    ///
    /// **The step is inherited from the type that produced the value, not chosen by the
    /// writer.** The rounding happens where the value is computed, so a run reading its
    /// observations straight from memory and a run reading them back from a psp see the
    /// same number — which is the oracle the whole psp path is checked against. This field
    /// exists so a reader can *interpret* the integer; it cannot make a file with a step
    /// the types did not produce (spec psp_record_encoding.md §5.1.1).
    FixedPoint { steps_per_unit: u32 },
    /// Bytes with a varint length in front — an allele's sequence, a chain-id list.
    LengthPrefixedBytes,
}

// ---------------------------------------------------------------------
// Framing — the bytes around the TOML body
// ---------------------------------------------------------------------

/// The 4-byte head magic. Printable ASCII with the newline last, so `head` on a psp opens
/// the TOML body on a following line.
///
/// **`NGP\n` and not production's `PSP\n`, deliberately.** Both formats use the extension
/// `.psp` and both will sit on the same disks while ng is being built, so the first four
/// bytes are what tells a reader — or a person running `file` — which one it is holding.
/// Sharing the magic would leave the two to disagree later, over the TOML body, with an
/// error that names a missing key rather than the wrong format.
pub const HEAD_MAGIC: [u8; 4] = *b"NGP\n";

/// The line that closes the header. **A cross-check, not a boundary**: the declared body
/// length is authoritative, and this catches a length that disagrees with the bytes.
pub const HEAD_SENTINEL: &[u8; 17] = b"---END-HEADER---\n";

/// Magic, the 8-byte body length, and the sentinel — everything in the header that is not
/// the TOML body.
pub const HEADER_FRAMING_BYTES: usize = HEAD_MAGIC.len() + 8 + HEAD_SENTINEL.len();

/// The largest TOML body this reader will read: 16 MiB less the framing.
///
/// **Checked before anything is allocated**, so a corrupt or hostile length field cannot
/// drive a large allocation on its own say-so.
///
/// **It is really a limit on how many contigs a reference may have**, and that is why it is
/// 16 MiB rather than production's 1 MiB (the owner, 2026-08-26). Every contig costs a
/// `[[contig]]` table in the body, so the cap measured out at about **11,300 contigs with
/// their MD5s and about 35,000 without** — and a draft assembly with tens of thousands of
/// scaffolds is ordinary in plant genomics, which is this caller's own domain. The header is
/// written once and read once per file, so nothing on any hot path notices the room.
pub const MAX_HEADER_BODY_BYTES: u64 = (16 * 1024 * 1024) - HEADER_FRAMING_BYTES as u64;

/// The format this writer produces and this reader understands. A file whose **major**
/// differs is refused as [`PspReadError::UnsupportedVersion`], not read as damaged.
pub const FORMAT_VERSION: (u16, u16) = (1, 0);

/// The largest integer a TOML value can carry: TOML's integer is signed 64-bit.
///
/// **A limit of the file's own syntax, not a choice.** Two of the header's numbers are `u64`
/// in Rust — a contig's length and the genomic block size — and the `toml` crate will
/// *serialise* a value above this that its own parser then refuses. So the rule belongs in
/// [`check_rules`] with the others, or the writer produces a file its own reader rejects.
const MAX_TOML_INTEGER: u64 = i64::MAX as u64;

/// The smallest look-back window zstd will accept, as its exponent: 2^10 bytes.
pub const MIN_LOOK_BACK_WINDOW_LOG: u8 = 10;

/// The largest look-back window zstd will accept, as its exponent: 2^31 bytes.
///
/// **This is what zstd allows, not what a reader should budget for.** A file may legally
/// declare a 2 GiB window and a cohort reader holding one file per sample cannot afford it —
/// that refusal is [`PspReadError::WindowTooLarge`], a reader's budget rather than a format
/// rule, and it belongs to the reader that Milestone D builds.
pub const MAX_LOOK_BACK_WINDOW_LOG: u8 = 31;

/// The look-back window a writer takes when nothing else says otherwise: 2^15 = 32 kB.
///
/// **This is what the memory measurements were taken at** (spec §4.2, §5.2), and it is the
/// number that makes an open file cost 0.34 MB rather than production's 2.6 MB. It is a
/// starting value recorded in every file it writes, not a property of the format — a reader
/// is driven by what the file declares.
pub const DEFAULT_LOOK_BACK_WINDOW_LOG: u8 = 15;

/// The genomic block size a writer takes when nothing else says otherwise: 100 kb of
/// reference per psp block (the owner, 2026-08-25; spec §4.1).
///
/// **Not an optimum, and the spec says so.** Measured on a tomato accession at three reads a
/// position, bytes a record barely moves across a two-hundred-fold range of this number —
/// 4.629 at 20 kb against 4.626 at 1,000 kb — because the capped look-back window means a
/// larger block gives the match finder nothing extra. 100 kb is a round number in the flat
/// part of that curve, and spec §4.1 records that 1,000 kb is live and may be better.
pub const DEFAULT_GENOMIC_BLOCK_SIZE_BP: Bp = Bp(100_000);

/// The byte ceiling a writer takes when nothing else says otherwise: **none at all**.
///
/// **Off because the value is not known yet**, not because none is wanted. Spec §12 question 2
/// is open, and what settles it is the block-size distribution on a whole-genome
/// deep-coverage sample, which nothing has produced. Leaving it off costs nothing measurable
/// on the data there is: a 100 kb grid with a 1 MiB ceiling gives 4.628 bytes a record against
/// 4.627 without, on a tomato accession at three reads a position (spec §4.1). What a ceiling
/// would be *for* is the other end of the range — at 279 reads a position a fully covered
/// 100 kb block is about 1.6 MB, which is a large thing to hold while writing.
///
/// **Named rather than spelled `None` at each site**, so that when the question closes the
/// value has somewhere to land with the reasoning beside it.
pub const DEFAULT_BLOCK_BYTE_CEILING: Option<u32> = None;

// ---------------------------------------------------------------------
// Build, encode, parse, validate
// ---------------------------------------------------------------------

/// A rule the header broke, before it is dressed as a writer's error or a reader's.
///
/// The rules are written once and checked on both sides: the writer refuses to produce a
/// file that breaks one, and the reader refuses to believe a file that does. Two copies of
/// one rule are two things that can disagree.
struct BrokenRule {
    field: String,
    reason: String,
}

impl BrokenRule {
    fn new(field: impl Into<String>, reason: impl Into<String>) -> Self {
        BrokenRule {
            field: field.into(),
            reason: reason.into(),
        }
    }
}

impl Header {
    /// The header's bytes: magic, the body's length, the TOML body, the sentinel.
    ///
    /// Refuses a header that breaks any rule the format requires rather than writing a file
    /// no reader would accept.
    pub fn encode(&self) -> Result<Vec<u8>, PspWriteError> {
        if self.format_version != FORMAT_VERSION {
            return Err(PspWriteError::InvalidHeaderField {
                field: "format-version".to_string(),
                reason: format!(
                    "this writer produces {}.{} only",
                    FORMAT_VERSION.0, FORMAT_VERSION.1
                ),
            });
        }
        if let Err(broken) = check_rules(self) {
            return Err(PspWriteError::InvalidHeaderField {
                field: broken.field,
                reason: broken.reason,
            });
        }

        let body = toml::to_string_pretty(&WireHeader::from(self)).map_err(|e| {
            PspWriteError::InvalidHeaderField {
                field: "header".to_string(),
                reason: format!("could not be written as TOML: {e}"),
            }
        })?;
        let body = body.as_bytes();
        if body.len() as u64 > MAX_HEADER_BODY_BYTES {
            return Err(PspWriteError::InvalidHeaderField {
                field: "header".to_string(),
                reason: format!(
                    "the TOML body is {} bytes; the format allows {MAX_HEADER_BODY_BYTES}",
                    body.len()
                ),
            });
        }

        let mut bytes = Vec::with_capacity(HEADER_FRAMING_BYTES + body.len());
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&(body.len() as u64).to_le_bytes());
        bytes.extend_from_slice(body);
        bytes.extend_from_slice(HEAD_SENTINEL);
        Ok(bytes)
    }

    /// The header at the start of `bytes`, and how many bytes it occupied — which is where
    /// the first psp block begins.
    ///
    /// `path` names the file in any error; it is not read.
    ///
    /// **The version is checked before the rest of the body is interpreted.** A file written
    /// by a newer format has to come back as
    /// [`UnsupportedVersion`](PspReadError::UnsupportedVersion) — *upgrade the reader* —
    /// rather than as a header full of keys this reader does not know, which would send
    /// whoever sees it looking for a damaged file.
    pub fn decode(bytes: &[u8], path: &Path) -> Result<(Self, usize), PspReadError> {
        let malformed = |reason: String| PspReadError::MalformedHeader {
            path: path.to_path_buf(),
            reason,
            source: None,
        };
        let malformed_by = |reason: &str, cause: Box<dyn std::error::Error + Send + Sync>| {
            PspReadError::MalformedHeader {
                path: path.to_path_buf(),
                reason: format!("{reason}: {cause}"),
                source: Some(cause),
            }
        };

        // The magic and the length are read from a fixed-size prefix, so this is the one
        // place the length of `bytes` has to be checked before anything is indexed.
        let Some((magic, after_magic)) = bytes.split_at_checked(HEAD_MAGIC.len()) else {
            return Err(malformed(format!(
                "the file is {} bytes, too short to hold a header's magic",
                bytes.len()
            )));
        };
        if magic != HEAD_MAGIC {
            let mut found = [0u8; 4];
            found.copy_from_slice(magic);
            return Err(PspReadError::NotAnNgPsp {
                path: path.to_path_buf(),
                found,
                expected: HEAD_MAGIC,
            });
        }
        let Some((declared_length, _)) = after_magic.split_first_chunk::<8>() else {
            return Err(malformed(format!(
                "the file is {} bytes, too short to hold the header body's length",
                bytes.len()
            )));
        };

        let body_bytes = u64::from_le_bytes(*declared_length);
        if body_bytes == 0 || body_bytes > MAX_HEADER_BODY_BYTES {
            return Err(malformed(format!(
                "declares a {body_bytes}-byte header body; this reader allows 1 to \
                 {MAX_HEADER_BODY_BYTES}"
            )));
        }

        let body_at = HEAD_MAGIC.len() + 8;
        // `body_bytes` is at most `MAX_HEADER_BODY_BYTES`, checked above, so neither sum can
        // overflow a `usize` on any target this builds for.
        let sentinel_at = body_at + body_bytes as usize;
        let header_bytes = sentinel_at + HEAD_SENTINEL.len();
        if bytes.len() < header_bytes {
            return Err(malformed(format!(
                "declares a {body_bytes}-byte header body but only {} bytes follow the length",
                bytes.len() - body_at
            )));
        }
        if &bytes[sentinel_at..header_bytes] != HEAD_SENTINEL.as_slice() {
            return Err(malformed(
                "the header's declared length does not reach its closing line".to_string(),
            ));
        }

        let body = std::str::from_utf8(&bytes[body_at..sentinel_at])
            .map_err(|e| malformed_by("the header body is not valid UTF-8", Box::new(e)))?;

        let table: toml::Table = body
            .parse()
            .map_err(|e| malformed_by("the header body is not valid TOML", Box::new(e)))?;

        let format_version = version_in(&table, path)?;
        if format_version.0 != FORMAT_VERSION.0 {
            return Err(PspReadError::UnsupportedVersion {
                path: path.to_path_buf(),
                found: format_version,
                supported: FORMAT_VERSION,
            });
        }

        // **The body is parsed twice, and it has to be.** The version must be read before the
        // body is interpreted as this version's shape, so it comes from the bare table above.
        // Feeding that table to the wire types instead of re-reading the text loses the one
        // thing a `toml::Table` cannot carry back out: `writer.created` arrives as a string
        // rather than as a TOML datetime, and every header fails to parse. A header is about a
        // kilobyte and is read once per file open, so the second pass costs nothing that
        // matters — and the alternative was measured by trying it.
        let wire: WireHeader = toml::from_str(body).map_err(|e| {
            malformed_by(
                "the header body is not a header this reader can read",
                Box::new(e),
            )
        })?;
        let header = wire
            .into_header(format_version)
            .map_err(|broken| malformed(format!("{}: {}", broken.field, broken.reason)))?;
        check_rules(&header)
            .map_err(|broken| malformed(format!("{}: {}", broken.field, broken.reason)))?;

        Ok((header, header_bytes))
    }
}

/// The format version alone, read from the body's bare TOML table before anything else in it
/// is interpreted.
///
/// **This is what keeps the header plain text.** A reader has to be able to learn the version
/// of a file it cannot otherwise read, so the version is taken from the table before the body
/// is deserialised into types this version's reader knows (spec §3.1).
fn version_in(table: &toml::Table, path: &Path) -> Result<(u16, u16), PspReadError> {
    let malformed = |reason: String| PspReadError::MalformedHeader {
        path: path.to_path_buf(),
        reason,
        source: None,
    };

    let Some(spelled) = table.get("format-version") else {
        return Err(malformed("the header has no format-version".to_string()));
    };
    let Some(spelled) = spelled.as_str() else {
        return Err(malformed(format!(
            "format-version is {spelled}, which is not a MAJOR.MINOR string"
        )));
    };
    let Some((major, minor)) = spelled.split_once('.') else {
        return Err(malformed(format!(
            "format-version {spelled:?} is not MAJOR.MINOR"
        )));
    };
    match (major.parse::<u16>(), minor.parse::<u16>()) {
        (Ok(major), Ok(minor)) => Ok((major, minor)),
        _ => Err(malformed(format!(
            "format-version {spelled:?} is not MAJOR.MINOR, each part a number below \
             {}",
            u16::MAX
        ))),
    }
}

/// Every rule a header must satisfy, in one place, checked on both sides.
fn check_rules(header: &Header) -> Result<(), BrokenRule> {
    if header.sample.trim().is_empty() {
        return Err(BrokenRule::new("sample", "is empty"));
    }
    check_basename("reference.name", &header.reference.name)?;

    if header.contigs.is_empty() {
        return Err(BrokenRule::new(
            "contig",
            "is empty; a psp's coordinates mean nothing without the contig list they index",
        ));
    }
    let mut seen = std::collections::HashSet::with_capacity(header.contigs.len());
    for contig in &header.contigs {
        if contig.name.trim().is_empty() {
            return Err(BrokenRule::new("contig.name", "is empty"));
        }
        if contig.name.chars().any(char::is_whitespace) {
            return Err(BrokenRule::new(
                "contig.name",
                format!("{:?} holds whitespace; a SAM @SQ SN cannot", contig.name),
            ));
        }
        if contig.length == 0 {
            return Err(BrokenRule::new(
                "contig.length",
                format!("{} is zero bases long", contig.name),
            ));
        }
        if contig.length > MAX_TOML_INTEGER {
            return Err(BrokenRule::new(
                "contig.length",
                format!(
                    "{} is {} bases; a TOML integer is signed, so a header cannot carry more \
                     than {MAX_TOML_INTEGER}",
                    contig.name, contig.length
                ),
            ));
        }
        if !seen.insert(contig.name.as_str()) {
            return Err(BrokenRule::new(
                "contig.name",
                format!(
                    "{:?} appears twice; a ContigId must name one contig",
                    contig.name
                ),
            ));
        }
    }

    check_basename("writer.input-reference", &header.writer.input_reference)?;
    for input in &header.writer.input_alignments {
        check_basename("writer.input-alignments", input)?;
    }
    for (name, value) in &header.writer.parameters {
        if let ParameterValue::Float(number) = value
            && !number.is_finite()
        {
            return Err(BrokenRule::new(
                "writer.parameters",
                format!("{name} is {number}, which has no TOML spelling"),
            ));
        }
    }

    check_manifest(&header.manifest)
}

/// A path recorded in the header must be a basename.
///
/// **Provenance, not tidiness**: a file that recorded the producer's directory layout would
/// carry it to everyone the file is shared with, and a run reproduced on another host cannot
/// use the path anyway.
fn check_basename(field: &str, spelled: &str) -> Result<(), BrokenRule> {
    if spelled.trim().is_empty() {
        return Err(BrokenRule::new(field, "is empty"));
    }
    // **One ordinary component, not merely one component.** Counting components alone lets
    // `/`, `.` and `..` through, because each of those is a single component of its own kind —
    // so the rule read as enforced while `..` recorded a directory the same way a path would.
    let mut components = Path::new(spelled).components();
    let sole_ordinary_component = matches!(
        (components.next(), components.next()),
        (Some(std::path::Component::Normal(_)), None)
    );
    if !sole_ordinary_component {
        return Err(BrokenRule::new(
            field,
            format!("{spelled:?} holds a directory component; only a basename is recorded"),
        ));
    }
    Ok(())
}

fn check_manifest(manifest: &Manifest) -> Result<(), BrokenRule> {
    if manifest.genomic_block_size_bp.get() == 0 {
        return Err(BrokenRule::new(
            "manifest.genomic-block-size-bp",
            "is zero; the block cut is a grid on the coordinate and a zero grid has no cells",
        ));
    }
    if manifest.genomic_block_size_bp.get() > MAX_TOML_INTEGER {
        return Err(BrokenRule::new(
            "manifest.genomic-block-size-bp",
            format!(
                "is {}; a TOML integer is signed, so a header cannot carry more than \
                 {MAX_TOML_INTEGER}",
                manifest.genomic_block_size_bp.get()
            ),
        ));
    }
    if manifest.block_byte_ceiling == Some(0) {
        return Err(BrokenRule::new(
            "manifest.block-byte-ceiling",
            "is zero; a ceiling no block can stay under gives every record a block of its own",
        ));
    }
    if !(MIN_LOOK_BACK_WINDOW_LOG..=MAX_LOOK_BACK_WINDOW_LOG)
        .contains(&manifest.look_back_window_log)
    {
        return Err(BrokenRule::new(
            "manifest.look-back-window-log",
            format!(
                "is {}; zstd takes a look-back window between 2^{MIN_LOOK_BACK_WINDOW_LOG} and \
                 2^{MAX_LOOK_BACK_WINDOW_LOG} bytes",
                manifest.look_back_window_log
            ),
        ));
    }

    if manifest.fields.is_empty() {
        return Err(BrokenRule::new(
            "manifest.field",
            "is empty; a reader is driven by the file's declared encodings and there are none",
        ));
    }
    let mut seen = std::collections::HashSet::with_capacity(manifest.fields.len());
    for field in &manifest.fields {
        if field.name.0.trim().is_empty() {
            return Err(BrokenRule::new("manifest.field.name", "is empty"));
        }
        if !seen.insert(field.name.0.as_str()) {
            return Err(BrokenRule::new(
                "manifest.field.name",
                format!(
                    "{:?} appears twice; the fields are in encoding order and a repeat leaves \
                     the reader no way to know which is which",
                    field.name.0
                ),
            ));
        }
        check_encoding(&field.name, field.encoding)?;
    }
    Ok(())
}

/// The widths a fixed-width integer field may declare.
///
/// **Named rather than spelled once into the check and again into the message.** The two
/// drifted apart under review: widening the bound in one and leaving the other alone passed
/// every test while admitting a file no decoder can read.
pub const FIXED_INTEGER_WIDTHS_BYTES: [u8; 4] = [1, 2, 4, 8];

/// The widths a raw IEEE float field may declare — `f32` and `f64`. No other width has an
/// IEEE meaning.
pub const IEEE_FLOAT_WIDTHS_BYTES: [u8; 2] = [4, 8];

fn check_encoding(name: &FieldName, encoding: FieldEncoding) -> Result<(), BrokenRule> {
    let field = "manifest.field.encoding";
    match encoding {
        FieldEncoding::FixedWidthInteger { width_bytes } => {
            if !FIXED_INTEGER_WIDTHS_BYTES.contains(&width_bytes) {
                return Err(BrokenRule::new(
                    field,
                    format!(
                        "{:?} is a {width_bytes}-byte fixed-width integer; the widths are {:?}",
                        name.0, FIXED_INTEGER_WIDTHS_BYTES
                    ),
                ));
            }
        }
        FieldEncoding::IeeeFloat { width_bytes } => {
            if !IEEE_FLOAT_WIDTHS_BYTES.contains(&width_bytes) {
                return Err(BrokenRule::new(
                    field,
                    format!(
                        "{:?} is a {width_bytes}-byte IEEE float; the widths are {:?}",
                        name.0, IEEE_FLOAT_WIDTHS_BYTES
                    ),
                ));
            }
        }
        FieldEncoding::FixedPoint { steps_per_unit } => {
            if steps_per_unit == 0 {
                return Err(BrokenRule::new(
                    field,
                    format!("{:?} counts steps of 1/0", name.0),
                ));
            }
        }
        FieldEncoding::Varint
        | FieldEncoding::SignedVarint
        | FieldEncoding::LengthPrefixedBytes => {}
    }
    Ok(())
}

// ---------------------------------------------------------------------
// The wire types — what the TOML body actually is
// ---------------------------------------------------------------------
//
// Kept apart from the public types on purpose, which is production's split
// (`src/psp/header.rs`): the public types carry strong types — a `(u16, u16)` version, a
// `Bp`, a `FieldEncoding` — and these carry what TOML has, which is strings, integers and
// tables. Every field is a value before it is a table, because TOML requires it.

/// Every MD5 travels as 32 lowercase hex characters, which is what a SAM `@SQ M5` is.
fn hex_of(digest: [u8; 16]) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn digest_of(field: &str, spelled: &str) -> Result<[u8; 16], BrokenRule> {
    let wrong = || {
        BrokenRule::new(
            field,
            format!("{spelled:?} is not 32 lowercase hex characters"),
        )
    };
    if spelled.len() != 32 || !spelled.bytes().all(|b| b.is_ascii_hexdigit()) {
        return Err(wrong());
    }
    if spelled.bytes().any(|b| b.is_ascii_uppercase()) {
        return Err(wrong());
    }
    let mut digest = [0u8; 16];
    for (byte, pair) in digest.iter_mut().zip(spelled.as_bytes().chunks_exact(2)) {
        let pair = std::str::from_utf8(pair).map_err(|_| wrong())?;
        *byte = u8::from_str_radix(pair, 16).map_err(|_| wrong())?;
    }
    Ok(digest)
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireHeader {
    format_version: String,
    sample: String,
    reference: WireReference,
    #[serde(default)]
    contig: Vec<WireContig>,
    writer: WireWriter,
    manifest: WireManifest,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireReference {
    name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    md5: Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireContig {
    name: String,
    length: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    md5: Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireWriter {
    tool: String,
    version: String,
    subcommand: String,
    input_alignments: Vec<String>,
    input_reference: String,
    command_line: String,
    created: toml::value::Datetime,
    #[serde(default)]
    parameters: BTreeMap<String, ParameterValue>,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireManifest {
    genomic_block_size_bp: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_byte_ceiling: Option<u32>,
    look_back_window_log: u8,
    #[serde(default)]
    field: Vec<WireFieldSpec>,
}

/// One field's declaration, **flat rather than nested**: `encoding` names the scheme and
/// `width-bytes` or `steps-per-unit` carries its one parameter. A nested table per field
/// would read as `[manifest.field.encoding]` inside an array of tables, which is legal TOML
/// and hard to read in `head` — and readability is the reason the header is text at all.
///
/// The two parameters are named for what they are rather than for their Rust types: a
/// header carries `fixed-width-integer` four lines from `fixed-point`, and *width* against
/// *step* is what tells a person reading it which is which.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireFieldSpec {
    name: String,
    encoding: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    width_bytes: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    steps_per_unit: Option<u32>,
}

impl From<&Header> for WireHeader {
    fn from(header: &Header) -> Self {
        WireHeader {
            format_version: format!("{}.{}", header.format_version.0, header.format_version.1),
            sample: header.sample.clone(),
            reference: WireReference {
                name: header.reference.name.clone(),
                md5: header.reference.md5.map(hex_of),
            },
            contig: header
                .contigs
                .iter()
                .map(|contig| WireContig {
                    name: contig.name.clone(),
                    length: contig.length,
                    md5: contig.md5.map(hex_of),
                })
                .collect(),
            writer: WireWriter {
                tool: header.writer.tool.clone(),
                version: header.writer.version.clone(),
                subcommand: header.writer.subcommand.clone(),
                input_alignments: header.writer.input_alignments.clone(),
                input_reference: header.writer.input_reference.clone(),
                command_line: header.writer.command_line.clone(),
                created: header.writer.created,
                parameters: header.writer.parameters.clone(),
            },
            manifest: WireManifest {
                genomic_block_size_bp: header.manifest.genomic_block_size_bp.get(),
                block_byte_ceiling: header.manifest.block_byte_ceiling,
                look_back_window_log: header.manifest.look_back_window_log,
                field: header
                    .manifest
                    .fields
                    .iter()
                    .map(|field| {
                        let (encoding, width_bytes, steps_per_unit) = field.encoding.spelled();
                        WireFieldSpec {
                            name: field.name.0.clone(),
                            encoding: encoding.to_string(),
                            width_bytes,
                            steps_per_unit,
                        }
                    })
                    .collect(),
            },
        }
    }
}

impl WireHeader {
    fn into_header(self, format_version: (u16, u16)) -> Result<Header, BrokenRule> {
        let reference = ReferenceIdentity {
            name: self.reference.name,
            md5: self
                .reference
                .md5
                .map(|spelled| digest_of("reference.md5", &spelled))
                .transpose()?,
        };
        let contigs = self
            .contig
            .into_iter()
            .map(|contig| {
                Ok(ContigIdentity {
                    name: contig.name,
                    length: contig.length,
                    md5: contig
                        .md5
                        .map(|spelled| digest_of("contig.md5", &spelled))
                        .transpose()?,
                })
            })
            .collect::<Result<Vec<_>, BrokenRule>>()?;
        let fields = self
            .manifest
            .field
            .into_iter()
            .map(|field| {
                Ok(FieldSpec {
                    encoding: encoding_of(&field)?,
                    name: FieldName(field.name),
                })
            })
            .collect::<Result<Vec<_>, BrokenRule>>()?;

        Ok(Header {
            format_version,
            sample: self.sample,
            reference,
            contigs,
            writer: WriterProvenance {
                tool: self.writer.tool,
                version: self.writer.version,
                subcommand: self.writer.subcommand,
                input_alignments: self.writer.input_alignments,
                input_reference: self.writer.input_reference,
                command_line: self.writer.command_line,
                parameters: self.writer.parameters,
                created: self.writer.created,
            },
            manifest: Manifest {
                genomic_block_size_bp: Bp(self.manifest.genomic_block_size_bp),
                block_byte_ceiling: self.manifest.block_byte_ceiling,
                look_back_window_log: self.manifest.look_back_window_log,
                fields,
            },
        })
    }
}

/// Every scheme the format has, listed once.
///
/// **This is the list, and both sides read it.** The writer takes each scheme's spelling
/// from [`FieldEncoding::spelled`] and the reader recognises exactly what this array
/// produces, so a seventh scheme cannot reach one side without reaching the other. When they
/// were two independent lists a scheme added to the writer's alone made
/// [`Header::encode`] write a file that [`Header::decode`] refused — the exact failure the
/// two-sided rule set exists to prevent, and it passed the whole suite.
///
/// The parameters here are placeholders: the array names the schemes, and each field's own
/// parameter is read from the header beside it.
const ALL_ENCODINGS: [FieldEncoding; 6] = [
    FieldEncoding::Varint,
    FieldEncoding::SignedVarint,
    FieldEncoding::FixedWidthInteger { width_bytes: 1 },
    FieldEncoding::IeeeFloat { width_bytes: 4 },
    FieldEncoding::FixedPoint { steps_per_unit: 1 },
    FieldEncoding::LengthPrefixedBytes,
];

impl FieldEncoding {
    /// What this scheme is called in the header, and the one parameter it carries there.
    ///
    /// **The only place a scheme's spelling is written.** Exhaustive on purpose: adding a
    /// scheme is a compile error here, which is what makes [`ALL_ENCODINGS`] the whole list
    /// rather than most of it.
    fn spelled(self) -> (&'static str, Option<u8>, Option<u32>) {
        match self {
            FieldEncoding::Varint => ("varint", None, None),
            FieldEncoding::SignedVarint => ("signed-varint", None, None),
            FieldEncoding::FixedWidthInteger { width_bytes } => {
                ("fixed-width-integer", Some(width_bytes), None)
            }
            FieldEncoding::IeeeFloat { width_bytes } => ("ieee-float", Some(width_bytes), None),
            FieldEncoding::FixedPoint { steps_per_unit } => {
                ("fixed-point", None, Some(steps_per_unit))
            }
            FieldEncoding::LengthPrefixedBytes => ("length-prefixed-bytes", None, None),
        }
    }
}

/// The declared scheme and its one parameter, together.
///
/// **A parameter that belongs to another scheme is refused rather than ignored.** A file
/// saying `encoding = "varint"` beside `steps-per-unit = 4096` was written by something that
/// meant one of the two, and reading it as a plain varint would silently multiply every
/// value in that field by 4,096.
fn encoding_of(field: &WireFieldSpec) -> Result<FieldEncoding, BrokenRule> {
    let named = "manifest.field.encoding";
    let wrong_parameter = |wanted: &str| {
        BrokenRule::new(
            named,
            format!(
                "{:?} is {:?} and carries a {wanted} it has no use for",
                field.name, field.encoding
            ),
        )
    };
    let missing = |wanted: &str| {
        BrokenRule::new(
            named,
            format!(
                "{:?} is {:?} and carries no {wanted}",
                field.name, field.encoding
            ),
        )
    };

    let scheme = ALL_ENCODINGS
        .iter()
        .find(|candidate| candidate.spelled().0 == field.encoding)
        .ok_or_else(|| {
            let known: Vec<&str> = ALL_ENCODINGS
                .iter()
                .map(|candidate| candidate.spelled().0)
                .collect();
            BrokenRule::new(
                named,
                format!(
                    "{:?} is {:?}, which is not one of {}",
                    field.name,
                    field.encoding,
                    known.join(", ")
                ),
            )
        })?;

    // Exhaustive, no wildcard: a scheme added to `ALL_ENCODINGS` has to say here which of the
    // two parameters it reads, rather than inheriting "takes none" in silence.
    let encoding = match scheme {
        FieldEncoding::Varint => FieldEncoding::Varint,
        FieldEncoding::SignedVarint => FieldEncoding::SignedVarint,
        FieldEncoding::LengthPrefixedBytes => FieldEncoding::LengthPrefixedBytes,
        FieldEncoding::FixedWidthInteger { .. } => FieldEncoding::FixedWidthInteger {
            width_bytes: field.width_bytes.ok_or_else(|| missing("width"))?,
        },
        FieldEncoding::IeeeFloat { .. } => FieldEncoding::IeeeFloat {
            width_bytes: field.width_bytes.ok_or_else(|| missing("width"))?,
        },
        FieldEncoding::FixedPoint { .. } => FieldEncoding::FixedPoint {
            steps_per_unit: field.steps_per_unit.ok_or_else(|| missing("step"))?,
        },
    };

    // Which parameters this scheme carries is the same answer the writer gives, so a
    // parameter it did not put there is one the file should not have.
    let (_, carries_width, carries_step) = encoding.spelled();
    if carries_width.is_none() && field.width_bytes.is_some() {
        return Err(wrong_parameter("width"));
    }
    if carries_step.is_none() && field.steps_per_unit.is_some() {
        return Err(wrong_parameter("step"));
    }
    Ok(encoding)
}
#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------
    // Fixtures
    // -----------------------------------------------------------------

    /// A header that says what a real tomato run says, minus the eleven contigs the shape of
    /// these tests does not need.
    fn a_written_header() -> Header {
        Header {
            format_version: FORMAT_VERSION,
            sample: "SRR7279481".to_string(),
            reference: ReferenceIdentity {
                name: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                md5: Some([0x0a; 16]),
            },
            contigs: vec![
                ContigIdentity {
                    name: "SL4.0ch00".to_string(),
                    length: 9_643_250,
                    md5: Some([0x1b; 16]),
                },
                ContigIdentity {
                    name: "SL4.0ch01".to_string(),
                    length: 90_863_682,
                    md5: None,
                },
            ],
            writer: WriterProvenance {
                tool: "ng".to_string(),
                version: "0.1.0".to_string(),
                subcommand: "pileup".to_string(),
                input_alignments: vec!["SRR7279481.cram".to_string()],
                input_reference: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                command_line: "ng pileup --sample SRR7279481".to_string(),
                parameters: BTreeMap::from([
                    ("depth-cap".to_string(), ParameterValue::Integer(300)),
                    ("min-base-quality".to_string(), ParameterValue::Integer(20)),
                    ("support-share".to_string(), ParameterValue::Float(0.1)),
                    ("realign".to_string(), ParameterValue::Boolean(true)),
                ]),
                created: "2026-08-26T09:15:00Z"
                    .parse()
                    .expect("a valid RFC 3339 stamp"),
            },
            manifest: a_manifest(),
        }
    }

    /// One field per encoding the format has, so anything that walks the manifest meets all
    /// six rather than the two a minimal fixture would carry.
    fn a_manifest() -> Manifest {
        Manifest {
            genomic_block_size_bp: DEFAULT_GENOMIC_BLOCK_SIZE_BP,
            block_byte_ceiling: Some(1_048_576),
            look_back_window_log: DEFAULT_LOOK_BACK_WINDOW_LOG,
            fields: vec![
                FieldSpec {
                    name: FieldName("position-offset".to_string()),
                    encoding: FieldEncoding::Varint,
                },
                FieldSpec {
                    name: FieldName("coverage-step".to_string()),
                    encoding: FieldEncoding::SignedVarint,
                },
                FieldSpec {
                    name: FieldName("body-bytes".to_string()),
                    encoding: FieldEncoding::FixedWidthInteger { width_bytes: 4 },
                },
                FieldSpec {
                    name: FieldName("allele-bases".to_string()),
                    encoding: FieldEncoding::LengthPrefixedBytes,
                },
                FieldSpec {
                    name: FieldName("window-mean-coverage".to_string()),
                    encoding: FieldEncoding::FixedPoint { steps_per_unit: 4 },
                },
                FieldSpec {
                    name: FieldName("summed-log-error".to_string()),
                    encoding: FieldEncoding::FixedPoint {
                        steps_per_unit: 4_096,
                    },
                },
                FieldSpec {
                    name: FieldName("raw-escape-hatch".to_string()),
                    encoding: FieldEncoding::IeeeFloat { width_bytes: 8 },
                },
            ],
        }
    }

    /// The error a call had to produce. `expect_err` would do, but it needs `Debug` on the
    /// success type and gives no room for a formatted message naming which row failed.
    fn refusal<T, E>(outcome: Result<T, E>, what: &str) -> E {
        match outcome {
            Err(refused) => refused,
            Ok(_) => panic!("{what}"),
        }
    }

    fn decoded(bytes: &[u8]) -> Result<(Header, usize), PspReadError> {
        Header::decode(bytes, Path::new("SRR7279481.psp"))
    }

    /// The TOML body alone, sliced out by the declared length.
    ///
    /// **Never `String::from_utf8` over the whole framed header.** The eight bytes of length
    /// between the magic and the body are binary, and whether they happen to be valid UTF-8
    /// depends on how long the body is — so a test that decodes the whole frame passes until
    /// somebody adds a header key, then panics inside the helper instead of failing its own
    /// assertion. It happened during this module's review.
    fn body_of(bytes: &[u8]) -> &str {
        let declared = u64::from_le_bytes(
            bytes[HEAD_MAGIC.len()..HEAD_MAGIC.len() + 8]
                .try_into()
                .expect("eight bytes"),
        ) as usize;
        let body_at = HEAD_MAGIC.len() + 8;
        std::str::from_utf8(&bytes[body_at..body_at + declared]).expect("the body is UTF-8")
    }

    /// Frame a TOML body into a header's bytes **without running the writer's rules over
    /// it**, so a file holding a value the writer would never produce can be handed to the
    /// reader.
    fn framed(body: &str) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&(body.len() as u64).to_le_bytes());
        bytes.extend_from_slice(body.as_bytes());
        bytes.extend_from_slice(HEAD_SENTINEL);
        bytes
    }

    /// Frame a header that breaks a rule, so the reader meets it in a file.
    ///
    /// A non-finite parameter has no TOML spelling at all, so it is written as the string
    /// `nan` and the reader meets it as a bad value rather than a bad float; every other
    /// broken header here serialises as itself.
    fn smuggle(header: &Header) -> Vec<u8> {
        let body = toml::to_string_pretty(&WireHeader::from(header)).unwrap_or_else(|_| {
            let mut without_the_unspellable = header.clone();
            without_the_unspellable.writer.parameters.insert(
                "share".to_string(),
                ParameterValue::String("nan".to_string()),
            );
            toml::to_string_pretty(&WireHeader::from(&without_the_unspellable))
                .expect("everything else has a TOML spelling")
        });
        framed(&body)
    }

    /// Frame a header whose manifest holds exactly the one field declaration given, so a
    /// declaration no `FieldEncoding` can represent can still be put in a file.
    fn one_field_declared_as(declaration: &str) -> Vec<u8> {
        let whole = a_written_header();
        let body = toml::to_string_pretty(&WireHeader::from(&whole)).expect("encodes");
        let up_to_the_fields = body
            .split("[[manifest.field]]")
            .next()
            .expect("the manifest declares fields")
            .to_string();
        framed(&format!(
            "{up_to_the_fields}[[manifest.field]]\n{declaration}"
        ))
    }

    // -----------------------------------------------------------------
    // The round trip
    // -----------------------------------------------------------------

    /// Everything in the header comes back, field for field. **Equality here is strict** —
    /// the contig row is this module's own type precisely so that a dropped MD5 is a
    /// difference and not a wildcard match.
    #[test]
    fn a_header_round_trips_field_for_field() {
        let written = a_written_header();
        let bytes = written.encode().expect("a valid header encodes");
        let (read_back, header_bytes) = decoded(&bytes).expect("its own bytes parse");
        assert_eq!(read_back, written);
        assert_eq!(header_bytes, bytes.len());
    }

    /// The second half of `decode`'s return is **where block 0 begins**, and on a real file
    /// that is not the end of the buffer. A header alone cannot tell the two apart: the
    /// round-trip fixture's correct answer and its buffer length are the same number, so
    /// returning the buffer's length instead passed every test in this module.
    #[test]
    fn decode_reports_the_header_length_and_not_the_buffer_length() {
        let header = a_written_header();
        let encoded = header.encode().expect("a valid header encodes");
        let header_only = encoded.len();

        let mut with_a_block = encoded;
        with_a_block.extend_from_slice(&[0xAA; 4_096]);

        let (read_back, first_block_at) = decoded(&with_a_block).expect("the header parses");
        assert_eq!(read_back, header);
        assert_eq!(
            first_block_at, header_only,
            "block 0 begins where the header ends, not where the buffer does"
        );
        assert_eq!(
            with_a_block[first_block_at], 0xAA,
            "and that byte is the block's first"
        );
    }

    /// The two fields that would round-trip through a wrong step without complaining: a step
    /// of a quarter of a read and one of 1/4,096 of a natural log sit in the same manifest,
    /// so a decoder that read the step from the wrong field would still produce integers.
    #[test]
    fn two_fixed_point_steps_in_one_manifest_stay_apart() {
        let written = a_written_header();
        let bytes = written.encode().expect("a valid header encodes");
        let (read_back, _) = decoded(&bytes).expect("its own bytes parse");
        let steps: Vec<_> = read_back
            .manifest
            .fields
            .iter()
            .filter_map(|field| match field.encoding {
                FieldEncoding::FixedPoint { steps_per_unit } => {
                    Some((field.name.0.as_str(), steps_per_unit))
                }
                _ => None,
            })
            .collect();
        assert_eq!(
            steps,
            vec![("window-mean-coverage", 4), ("summed-log-error", 4_096)]
        );
    }

    /// Every encoding the format has survives the trip, at every width the rules allow.
    /// Without this, five of the seven legal widths were never written to a file at all —
    /// dropping `8` from the accepted set passed the whole suite.
    #[test]
    fn every_encoding_at_every_legal_width_round_trips() {
        let mut fields = vec![
            FieldSpec {
                name: FieldName("varint".to_string()),
                encoding: FieldEncoding::Varint,
            },
            FieldSpec {
                name: FieldName("signed".to_string()),
                encoding: FieldEncoding::SignedVarint,
            },
            FieldSpec {
                name: FieldName("bytes".to_string()),
                encoding: FieldEncoding::LengthPrefixedBytes,
            },
            FieldSpec {
                name: FieldName("finest-step".to_string()),
                encoding: FieldEncoding::FixedPoint { steps_per_unit: 1 },
            },
            FieldSpec {
                name: FieldName("widest-step".to_string()),
                encoding: FieldEncoding::FixedPoint {
                    steps_per_unit: u32::MAX,
                },
            },
        ];
        for width in FIXED_INTEGER_WIDTHS_BYTES {
            fields.push(FieldSpec {
                name: FieldName(format!("fixed-{width}")),
                encoding: FieldEncoding::FixedWidthInteger { width_bytes: width },
            });
        }
        for width in IEEE_FLOAT_WIDTHS_BYTES {
            fields.push(FieldSpec {
                name: FieldName(format!("ieee-{width}")),
                encoding: FieldEncoding::IeeeFloat { width_bytes: width },
            });
        }

        let mut written = a_written_header();
        written.manifest.fields = fields;
        let bytes = written.encode().expect("every legal width encodes");
        let (read_back, _) = decoded(&bytes).expect("and decodes");
        assert_eq!(read_back.manifest.fields, written.manifest.fields);
    }

    /// A contig as long as a TOML integer goes, and one base longer. **The wire format sets
    /// this limit, not the rule set**, which is how it was missed: `encode` wrote a `u64`
    /// above `i64::MAX` that its own reader then refused.
    #[test]
    fn the_widest_number_toml_can_carry_round_trips_and_one_more_is_refused() {
        let mut at_the_edge = a_written_header();
        at_the_edge.contigs[0].length = MAX_TOML_INTEGER;
        at_the_edge.manifest.genomic_block_size_bp = Bp(MAX_TOML_INTEGER);
        let bytes = at_the_edge
            .encode()
            .expect("the widest TOML integer encodes");
        let (read_back, _) = decoded(&bytes).expect("and decodes");
        assert_eq!(read_back, at_the_edge);

        for (what, break_it) in [
            (
                "a contig one base too long",
                Box::new(|header: &mut Header| header.contigs[0].length = MAX_TOML_INTEGER + 1)
                    as Box<dyn Fn(&mut Header)>,
            ),
            (
                "a block grid one base too wide",
                Box::new(|header: &mut Header| {
                    header.manifest.genomic_block_size_bp = Bp(MAX_TOML_INTEGER + 1)
                }),
            ),
        ] {
            let mut header = a_written_header();
            break_it(&mut header);
            let refused = refusal(header.encode(), &format!("the writer must refuse {what}"));
            let _ = refused;
        }
    }

    // -----------------------------------------------------------------
    // What the file looks like
    // -----------------------------------------------------------------

    /// The reason the header is text: `head` on a psp tells you what it is. The magic ends in
    /// a newline so the body starts on its own line, and the body is TOML a person can read.
    #[test]
    fn the_body_is_readable_toml_after_a_newline_terminated_magic() {
        let bytes = a_written_header().encode().expect("a valid header encodes");
        assert_eq!(&bytes[..4], b"NGP\n");
        let body = body_of(&bytes);
        for expected in [
            "format-version = \"1.0\"",
            "sample = \"SRR7279481\"",
            "genomic-block-size-bp = 100000",
            "look-back-window-log = 15",
            "encoding = \"fixed-point\"",
            "steps-per-unit = 4096",
            "encoding = \"fixed-width-integer\"",
            "width-bytes = 4",
        ] {
            assert!(body.contains(expected), "{expected:?} missing from: {body}");
        }
        assert_eq!(&bytes[bytes.len() - HEAD_SENTINEL.len()..], HEAD_SENTINEL);
    }

    /// A parameter is written as the value itself, not as a tagged pair, because the header
    /// is meant to be read by eye.
    #[test]
    fn a_parameter_is_written_as_its_bare_value() {
        let bytes = a_written_header().encode().expect("a valid header encodes");
        let body = body_of(&bytes);
        assert!(body.contains("depth-cap = 300"), "body was: {body}");
        assert!(body.contains("realign = true"), "body was: {body}");
    }

    /// Goal 5 is that the same sample gathered at any worker count gives the same bytes, so
    /// nothing in the header may depend on the order a map happened to be filled in. The
    /// parameters are built here in reverse alphabetical order and inserted one at a time,
    /// which is what a `HashMap` would reorder and a `BTreeMap` must not.
    #[test]
    fn the_same_header_encodes_to_the_same_bytes_whatever_order_it_was_built_in() {
        let canonical = a_written_header().encode().expect("encodes");

        let mut backwards = a_written_header();
        let mut reversed = BTreeMap::new();
        let mut entries: Vec<_> = backwards.writer.parameters.iter().collect();
        entries.sort_by(|left, right| right.0.cmp(left.0));
        let entries: Vec<_> = entries
            .into_iter()
            .map(|(key, value)| (key.clone(), value.clone()))
            .collect();
        for (key, value) in entries {
            reversed.insert(key, value);
        }
        backwards.writer.parameters = reversed;

        assert_eq!(backwards.encode().expect("encodes"), canonical);
    }

    // -----------------------------------------------------------------
    // Versions
    // -----------------------------------------------------------------

    /// A file from a later format has to say so. The version is read from a bare TOML table
    /// before the body is deserialised, so a body full of keys this reader has never seen
    /// still yields the right answer instead of a parse failure.
    #[test]
    fn a_newer_major_version_is_refused_as_a_version_and_not_as_damage() {
        let bytes = framed(
            "format-version = \"2.0\"\n\
             whatever-version-two-added = { shape = \"nothing here knows\" }\n",
        );
        match decoded(&bytes) {
            Err(PspReadError::UnsupportedVersion {
                found, supported, ..
            }) => {
                assert_eq!(found, (2, 0));
                assert_eq!(supported, FORMAT_VERSION);
            }
            other => panic!("expected UnsupportedVersion, got {other:?}"),
        }
    }

    /// A later *minor* of the same major is this reader's to read — that is what splitting
    /// the version into major and minor is **for**, and it only means something if a minor
    /// that actually added something still reads.
    ///
    /// The key added here is the one spec §3.1 says the header will gain: the observation
    /// reach ceiling. Before this test the module refused it as a damaged file, because the
    /// wire types were `deny_unknown_fields` and the only later-minor fixture added nothing.
    #[test]
    fn a_later_minor_that_added_a_key_still_reads() {
        let mut written = a_written_header();
        let bytes = written.encode().expect("a valid header encodes");
        let body = body_of(&bytes)
            .replace("format-version = \"1.0\"", "format-version = \"1.4\"")
            .replace(
                "sample = ",
                "observation-reach-ceiling-bp = 4000\nsample = ",
            );

        written.format_version = (1, 4);
        let (read_back, _) = decoded(&framed(&body)).expect("a later minor of this major reads");
        assert_eq!(
            read_back, written,
            "everything this reader knows comes back; the key it does not know is ignored"
        );
    }

    /// This writer produces one version. A caller that asks for another is told so, rather
    /// than being handed a file stamped with a version whose layout the writer does not
    /// actually produce.
    #[test]
    fn the_writer_refuses_to_stamp_a_version_it_does_not_produce() {
        let mut header = a_written_header();
        header.format_version = (1, 7);
        let refused = header.encode().expect_err("this writer produces 1.0");
        assert!(
            refused.to_string().contains("produces 1.0"),
            "got {refused}"
        );
    }

    // -----------------------------------------------------------------
    // Framing
    // -----------------------------------------------------------------

    /// Both formats use the extension `.psp`, so the first four bytes are what tells them
    /// apart — and being handed the wrong file is a different instruction from being handed a
    /// damaged one, so it is a different error.
    #[test]
    fn a_production_psp_is_refused_as_a_foreign_file_not_as_a_damaged_one() {
        let mut productions = Vec::new();
        productions.extend_from_slice(&crate::psp::header::HEAD_MAGIC);
        productions.extend_from_slice(&(1u64).to_le_bytes());
        productions.push(b'x');
        productions.extend_from_slice(HEAD_SENTINEL);

        match decoded(&productions) {
            Err(PspReadError::NotAnNgPsp { found, .. }) => {
                assert_eq!(found, crate::psp::header::HEAD_MAGIC);
            }
            other => panic!("expected NotAnNgPsp, got {other:?}"),
        }
    }

    /// The declared length is authoritative and the sentinel is the cross-check on it. A
    /// length one byte short leaves the sentinel where the reader does not look for it.
    #[test]
    fn a_length_that_disagrees_with_the_sentinel_is_refused() {
        let mut bytes = a_written_header().encode().expect("a valid header encodes");
        let body_bytes = u64::from_le_bytes(bytes[4..12].try_into().expect("eight bytes"));
        bytes[4..12].copy_from_slice(&(body_bytes - 1).to_le_bytes());

        let refused = decoded(&bytes).expect_err("a short length must not be believed");
        assert!(
            refused.to_string().contains("closing line"),
            "got {refused}"
        );
    }

    /// The length field is read before anything is allocated, so a corrupt one cannot ask for
    /// a large buffer on its own say-so — and the ceiling is this reader's, so moving it has
    /// to fail a test rather than pass quietly.
    #[test]
    fn a_declared_length_beyond_the_readers_ceiling_is_refused() {
        let good = a_written_header().encode().expect("a valid header encodes");
        for declared in [0u64, MAX_HEADER_BODY_BYTES + 1, u64::MAX] {
            let mut bytes = good.clone();
            bytes[4..12].copy_from_slice(&declared.to_le_bytes());
            let refused = decoded(&bytes).expect_err("a header cannot be that size");
            assert!(
                refused
                    .to_string()
                    .contains(&format!("this reader allows 1 to {MAX_HEADER_BODY_BYTES}")),
                "got {refused}"
            );
        }
    }

    /// **The header's size cap is really a cap on how many contigs a reference may have**, so
    /// what it has to be tested against is a fragmented assembly rather than a size.
    ///
    /// 30,000 scaffolds with their MD5s is an ordinary draft plant genome and is well past
    /// what the original 1 MiB cap allowed — measured in review at about 11,300 contigs with
    /// digests. Raised to 16 MiB by the owner, 2026-08-26.
    #[test]
    fn a_reference_of_thirty_thousand_scaffolds_still_fits_the_header() {
        let mut fragmented = a_written_header();
        fragmented.contigs = (0u32..30_000)
            .map(|scaffold| ContigIdentity {
                name: format!("scaffold_{scaffold}"),
                length: 1_000 + u64::from(scaffold),
                md5: Some([(scaffold % 256) as u8; 16]),
            })
            .collect();

        let bytes = fragmented.encode().expect("a fragmented assembly encodes");
        assert!(
            bytes.len() > 1024 * 1024,
            "the fixture must exceed the old 1 MiB cap to be testing anything; it is {} bytes",
            bytes.len()
        );
        let (read_back, _) = decoded(&bytes).expect("and decodes");
        assert_eq!(read_back.contigs.len(), 30_000);
        assert_eq!(read_back, fragmented);
    }

    /// A buffer too short to hold even the framing is refused rather than indexed into.
    #[test]
    fn a_buffer_shorter_than_the_framing_is_refused_at_every_length() {
        let good = a_written_header().encode().expect("a valid header encodes");
        for cut in 0..HEADER_FRAMING_BYTES {
            let refused = decoded(&good[..cut]);
            assert!(refused.is_err(), "a {cut}-byte buffer must not parse");
        }
    }

    // -----------------------------------------------------------------
    // The rules — every one, from both sides
    // -----------------------------------------------------------------

    /// One broken header: the rule it breaks, the change that breaks it, and the words the
    /// refusal has to carry. Named because the triple is otherwise complex enough that clippy
    /// asks for a name, and it reads better with one.
    type BrokenHeaderCase = (&'static str, Box<dyn Fn(&mut Header)>, &'static str);

    /// Each rule is written once and checked on both sides, so the writer refuses to make a
    /// file the reader would refuse to read. The table names the rule, a header that breaks
    /// it, and the words the message has to carry.
    ///
    /// **Every fixture sits one step over its rule's boundary**, not far beyond it. Three of
    /// them did not, and the mutations that widened those bounds survived the whole suite: a
    /// five-component absolute path let a one-directory relative path through, a look-back
    /// window of 40 let the floor drop from 10 to 0, and a three-character digest let the
    /// lowercase clause be deleted.
    #[test]
    fn every_rule_is_refused_by_the_writer_and_by_the_reader_alike() {
        let broken: Vec<BrokenHeaderCase> = vec![
            (
                "an empty sample name",
                Box::new(|header| header.sample = "  ".to_string()),
                "sample",
            ),
            (
                "a reference recorded with one directory above it",
                Box::new(|header| header.reference.name = "genomes/tomato.fa".to_string()),
                "directory component",
            ),
            (
                "a reference recorded as a parent directory",
                Box::new(|header| header.reference.name = "..".to_string()),
                "directory component",
            ),
            (
                "an empty reference name",
                Box::new(|header| header.reference.name = " ".to_string()),
                "reference.name",
            ),
            (
                "no contigs",
                Box::new(|header| header.contigs.clear()),
                "contig",
            ),
            (
                "an empty contig name",
                Box::new(|header| header.contigs[0].name = String::new()),
                "contig.name",
            ),
            (
                "a contig name holding a space",
                Box::new(|header| header.contigs[0].name = "SL4.0 ch00".to_string()),
                "holds whitespace",
            ),
            (
                "one contig named twice",
                Box::new(|header| header.contigs[1].name = header.contigs[0].name.clone()),
                "appears twice",
            ),
            (
                "a contig of no length",
                Box::new(|header| header.contigs[0].length = 0),
                "zero bases long",
            ),
            (
                "a contig longer than a TOML integer",
                Box::new(|header| header.contigs[0].length = MAX_TOML_INTEGER + 1),
                "a TOML integer is signed",
            ),
            (
                "the input reference recorded with its directory",
                Box::new(|header| header.writer.input_reference = "genomes/tomato.fa".to_string()),
                "writer.input-reference",
            ),
            (
                "an input alignment recorded with its directory",
                Box::new(|header| {
                    header.writer.input_alignments = vec!["runs/SRR7279481.cram".to_string()]
                }),
                "writer.input-alignments",
            ),
            (
                "a parameter that is not a number",
                Box::new(|header| {
                    header
                        .writer
                        .parameters
                        .insert("share".to_string(), ParameterValue::Float(f64::NAN));
                }),
                "no TOML spelling",
            ),
            (
                "a block grid of no width",
                Box::new(|header| header.manifest.genomic_block_size_bp = Bp(0)),
                "grid has no cells",
            ),
            (
                "a block grid wider than a TOML integer",
                Box::new(|header| header.manifest.genomic_block_size_bp = Bp(MAX_TOML_INTEGER + 1)),
                "a TOML integer is signed",
            ),
            (
                "a byte ceiling no block can stay under",
                Box::new(|header| header.manifest.block_byte_ceiling = Some(0)),
                "a block of its own",
            ),
            (
                "a look-back window one step below zstd's floor",
                Box::new(|header| {
                    header.manifest.look_back_window_log = MIN_LOOK_BACK_WINDOW_LOG - 1
                }),
                "2^10 and 2^31",
            ),
            (
                "a look-back window one step above zstd's ceiling",
                Box::new(|header| {
                    header.manifest.look_back_window_log = MAX_LOOK_BACK_WINDOW_LOG + 1
                }),
                "2^10 and 2^31",
            ),
            (
                "no declared fields",
                Box::new(|header| header.manifest.fields.clear()),
                "there are none",
            ),
            (
                "an empty field name",
                Box::new(|header| header.manifest.fields[0].name = FieldName(" ".to_string())),
                "manifest.field.name",
            ),
            (
                "one field declared twice",
                Box::new(|header| {
                    header.manifest.fields[1].name = header.manifest.fields[0].name.clone()
                }),
                "appears twice",
            ),
            (
                "a fixed-width integer one byte off a legal width",
                Box::new(|header| {
                    header.manifest.fields[2].encoding =
                        FieldEncoding::FixedWidthInteger { width_bytes: 3 }
                }),
                "the widths are [1, 2, 4, 8]",
            ),
            (
                "an IEEE float one byte off a legal width",
                Box::new(|header| {
                    header.manifest.fields[6].encoding = FieldEncoding::IeeeFloat { width_bytes: 5 }
                }),
                "the widths are [4, 8]",
            ),
            (
                "a fixed-point field counting steps of nothing",
                Box::new(|header| {
                    header.manifest.fields[4].encoding =
                        FieldEncoding::FixedPoint { steps_per_unit: 0 }
                }),
                "1/0",
            ),
        ];

        for (what, break_it, expected) in broken {
            let mut header = a_written_header();
            break_it(&mut header);

            let refused = refusal(header.encode(), &format!("the writer must refuse {what}"));
            assert!(
                refused.to_string().contains(expected),
                "the writer's message for {what} was {refused:?}, which does not say \
                 {expected:?}"
            );

            // The same rule met from the other side: a file that already holds the broken
            // value has to be refused when it is read, not only when it is written.
            //
            // **Two rows are refused earlier than the rule**, and that is the point of them:
            // a number above `i64::MAX` has no TOML spelling a parser will take back, so the
            // reader stops at the syntax. The rule exists so the *writer* stops first, rather
            // than producing a file only its own parser discovers is unreadable.
            let stopped_at_the_syntax = expected == "a TOML integer is signed";
            let expected_of_the_reader = if stopped_at_the_syntax {
                "is not valid TOML"
            } else {
                expected
            };
            let smuggled = smuggle(&header);
            let refused = refusal(
                decoded(&smuggled),
                &format!("the reader must refuse {what}"),
            );
            assert!(
                refused.to_string().contains(expected_of_the_reader),
                "the reader's message for {what} was {refused:?}, which does not say \
                 {expected_of_the_reader:?}"
            );
        }
    }

    /// An MD5 is 32 **lowercase** hex characters, which is what a SAM `@SQ M5` is. The
    /// uppercase clause is not decoration: a digest that came back in another case would
    /// re-encode to different bytes, and goal 5's byte identity would be gone.
    #[test]
    fn a_digest_is_refused_for_its_length_and_for_its_case_alike() {
        let lower = hex_of([0x1b; 16]);
        for (what, spelled) in [
            ("too short", "1b1b1b".to_string()),
            ("too long", format!("{lower}ab")),
            ("uppercase", lower.to_uppercase()),
            ("not hex at all", "z".repeat(32)),
        ] {
            let whole = a_written_header();
            let body = toml::to_string_pretty(&WireHeader::from(&whole))
                .expect("encodes")
                .replace(&lower, &spelled);
            let refused = refusal(
                decoded(&framed(&body)),
                &format!("a digest that is {what} must be refused"),
            );
            assert!(
                refused.to_string().contains("32 lowercase hex"),
                "for {what}, got {refused}"
            );
        }
    }

    // -----------------------------------------------------------------
    // The manifest's encodings
    // -----------------------------------------------------------------

    /// A field carrying a parameter that belongs to another scheme was written by something
    /// that meant one of the two. Reading it as the scheme it names would silently multiply
    /// every value in that field by the step it was not supposed to have.
    ///
    /// **All three arms**, because only one was covered and the check had a catch-all.
    #[test]
    fn an_encoding_carrying_the_wrong_parameter_is_refused() {
        for (declaration, what) in [
            (
                "name = \"x\"\nencoding = \"varint\"\nsteps-per-unit = 4096\n",
                "a varint with a step",
            ),
            (
                "name = \"x\"\nencoding = \"varint\"\nwidth-bytes = 4\n",
                "a varint with a width",
            ),
            (
                "name = \"x\"\nencoding = \"fixed-width-integer\"\nwidth-bytes = 4\nsteps-per-unit = 4096\n",
                "a fixed-width integer with a step",
            ),
            (
                "name = \"x\"\nencoding = \"fixed-point\"\nsteps-per-unit = 4096\nwidth-bytes = 4\n",
                "a fixed-point field with a width",
            ),
        ] {
            let refused = refusal(
                decoded(&one_field_declared_as(declaration)),
                &format!("{what} must be refused"),
            );
            assert!(
                refused.to_string().contains("has no use for"),
                "for {what}, got {refused}"
            );
        }
    }

    /// A width or a step that is simply missing is the same class of error as one that does
    /// not belong: the file does not say how to read the field.
    #[test]
    fn a_scheme_without_its_parameter_is_refused() {
        for (declaration, wanted) in [
            (
                "name = \"body-bytes\"\nencoding = \"fixed-width-integer\"\n",
                "carries no width",
            ),
            (
                "name = \"raw\"\nencoding = \"ieee-float\"\n",
                "carries no width",
            ),
            (
                "name = \"q-sum\"\nencoding = \"fixed-point\"\n",
                "carries no step",
            ),
        ] {
            let refused = refusal(
                decoded(&one_field_declared_as(declaration)),
                &format!("{declaration:?} must be refused"),
            );
            assert!(refused.to_string().contains(wanted), "got {refused}");
        }
    }

    /// An encoding this reader has never heard of is refused rather than guessed at, and the
    /// message lists what it does know — which is what tells whoever sees it whether to
    /// upgrade the reader or rebuild the file.
    #[test]
    fn an_unknown_encoding_is_refused_and_the_message_lists_the_known_ones() {
        let refused = decoded(&one_field_declared_as(
            "name = \"chain-ids\"\nencoding = \"roaring-bitmap\"\n",
        ))
        .expect_err("that is not one of the six");
        let said = refused.to_string();
        assert!(said.contains("roaring-bitmap"), "got {said}");
        for known in ALL_ENCODINGS {
            assert!(
                said.contains(known.spelled().0),
                "the message must list {:?}; got {said}",
                known.spelled().0
            );
        }
    }

    /// The writer's spellings and the reader's are one list. When they were two, a scheme
    /// added to the writer's alone made `encode` produce a file `decode` refused — and the
    /// whole suite stayed green.
    #[test]
    fn every_scheme_the_writer_can_spell_is_one_the_reader_recognises() {
        for scheme in ALL_ENCODINGS {
            let mut header = a_written_header();
            header.manifest.fields = vec![FieldSpec {
                name: FieldName("the-one-field".to_string()),
                encoding: scheme,
            }];
            let bytes = header
                .encode()
                .unwrap_or_else(|e| panic!("{scheme:?} must encode: {e}"));
            let (read_back, _) =
                decoded(&bytes).unwrap_or_else(|e| panic!("{scheme:?} must decode: {e}"));
            assert_eq!(read_back.manifest.fields[0].encoding, scheme);
        }
    }

    /// `ALL_ENCODINGS` has to be all of them. The exhaustive match is what makes adding a
    /// seventh scheme a compile error here rather than a file one side cannot read.
    #[test]
    fn the_encoding_list_holds_every_scheme_the_type_has() {
        let sample = FieldEncoding::Varint;
        let named = match sample {
            FieldEncoding::Varint
            | FieldEncoding::SignedVarint
            | FieldEncoding::FixedWidthInteger { .. }
            | FieldEncoding::IeeeFloat { .. }
            | FieldEncoding::FixedPoint { .. }
            | FieldEncoding::LengthPrefixedBytes => 6,
        };
        assert_eq!(
            ALL_ENCODINGS.len(),
            named,
            "a scheme was added to the type without being added to the list both sides read"
        );
        let mut spellings: Vec<&str> = ALL_ENCODINGS.iter().map(|e| e.spelled().0).collect();
        spellings.sort_unstable();
        spellings.dedup();
        assert_eq!(
            spellings.len(),
            ALL_ENCODINGS.len(),
            "two schemes share a name"
        );
    }

    // -----------------------------------------------------------------
    // The parser, over bytes it did not write
    // -----------------------------------------------------------------

    /// `decode` reads a file. **It must refuse damage rather than panic on it**, and nothing
    /// else in this module puts a byte in front of it that the writer did not produce.
    ///
    /// Deterministic on purpose — a fixed seed, so a failure is reproducible from the test
    /// name alone.
    #[test]
    fn decode_refuses_damaged_bytes_without_panicking() {
        let good = a_written_header().encode().expect("a valid header encodes");
        let mut seed = 0x2026_0826_u64;
        let mut next = move || {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            seed
        };

        let mut accepted = 0usize;
        for _ in 0..20_000 {
            let mut damaged = good.clone();
            for _ in 0..1 + (next() % 4) {
                let at = (next() as usize) % damaged.len();
                damaged[at] ^= 1u8 << (next() % 8);
            }
            if next() % 3 == 0 {
                damaged.truncate((next() as usize) % good.len());
            }
            if decoded(&damaged).is_ok() {
                accepted += 1;
            }
        }
        // Some flips land inside a sample name or a command line and leave a header that is
        // still well-formed — that is the format working, not a hole. What matters is that
        // every input either parsed or was refused, and none of them panicked.
        assert!(
            accepted < 20_000,
            "damage that changes the framing must not parse"
        );

        for length in 0..400usize {
            let mut arbitrary = Vec::with_capacity(length);
            for _ in 0..length {
                arbitrary.push((next() % 256) as u8);
            }
            assert!(
                decoded(&arbitrary).is_err(),
                "random bytes must not parse as a header"
            );
        }
    }

    /// A header that came off disk re-encodes to the bytes it came from. This is what lets an
    /// append reuse the header it finds rather than rewriting it, and it is the property a
    /// silently-dropped field would break.
    #[test]
    fn a_decoded_header_re_encodes_to_the_bytes_it_came_from() {
        let bytes = a_written_header().encode().expect("a valid header encodes");
        let (read_back, _) = decoded(&bytes).expect("its own bytes parse");
        assert_eq!(read_back.encode().expect("and re-encode"), bytes);
    }

    // -----------------------------------------------------------------
    // The types themselves
    // -----------------------------------------------------------------

    /// The header's contig list is what pins ng's coordinate space, so its equality has to be
    /// plain field equality. `crate::fasta::ContigEntry` — three fields the same — treats an
    /// absent MD5 as a wildcard, and a round-trip test written against *that* would pass while
    /// the encoder dropped every digest in the file.
    #[test]
    fn a_dropped_contig_md5_is_a_difference_here_and_a_wildcard_in_the_fasta_type() {
        let with_digest = ContigIdentity {
            name: "SL4.0ch01".to_string(),
            length: 90_863_682,
            md5: Some([7u8; 16]),
        };
        let without = ContigIdentity {
            md5: None,
            ..with_digest.clone()
        };
        assert_ne!(with_digest, without);

        let fasta_with_digest = crate::fasta::ContigEntry {
            name: "SL4.0ch01".to_string(),
            length: 90_863_682,
            md5: Some([7u8; 16]),
        };
        let fasta_without = crate::fasta::ContigEntry {
            md5: None,
            ..fasta_with_digest.clone()
        };
        assert_eq!(
            fasta_with_digest, fasta_without,
            "the fasta type's wildcard MD5 is why the header mints its own contig row"
        );
    }

    /// A field's encoding is part of the file's identity: two files whose manifests differ
    /// only in a fixed-point step hold different numbers behind the same integers, so the
    /// step has to take part in equality rather than being carried alongside it.
    #[test]
    fn a_fixed_point_step_is_part_of_the_encoding_not_a_note_beside_it() {
        let quarter_read = FieldEncoding::FixedPoint { steps_per_unit: 4 };
        let sixteenth_read = FieldEncoding::FixedPoint { steps_per_unit: 16 };
        assert_ne!(quarter_read, sixteenth_read);
        assert_ne!(quarter_read, FieldEncoding::Varint);
    }
}
