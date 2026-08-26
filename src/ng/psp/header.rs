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
    pub contigs: Vec<ContigEntry>,
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
pub struct ContigEntry {
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
    pub window_log: u8,
    /// One entry per field of a record, **in encoding order**.
    pub fields: Vec<FieldSpec>,
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
    /// A fixed-width little-endian integer of `bytes` width.
    Fixed { bytes: u8 },
    /// Raw IEEE bytes. **The escape hatch, not the default** — a float stored raw is a
    /// float two callers can disagree about in its last bits, and the three quantities
    /// that arrive as floating point are 70 % of the compressed file when stored this way
    /// (spec psp_record_encoding.md §5).
    Ieee { bytes: u8 },
    /// A count of steps of `1 / scale`, written as a varint.
    ///
    /// **The step is inherited from the type that produced the value, not chosen by the
    /// writer.** The rounding happens where the value is computed, so a run reading its
    /// observations straight from memory and a run reading them back from a psp see the
    /// same number — which is the oracle the whole psp path is checked against. This field
    /// exists so a reader can *interpret* the integer; it cannot make a file with a step
    /// the types did not produce (spec psp_record_encoding.md §5.1.1).
    FixedPoint { scale: u32 },
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

/// The largest TOML body this reader will read: 1 MiB less the framing.
///
/// **Checked before anything is allocated**, so a corrupt or hostile length field cannot
/// drive a large allocation on its own say-so.
pub const MAX_HEADER_BODY_BYTES: u64 = (1024 * 1024) - HEADER_FRAMING_BYTES as u64;

/// The format this writer produces and this reader understands. A file whose **major**
/// differs is refused as [`PspReadError::UnsupportedVersion`], not read as damaged.
pub const FORMAT_VERSION: (u16, u16) = (1, 0);

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
        };

        if bytes.len() < HEAD_MAGIC.len() + 8 {
            return Err(malformed(format!(
                "the file is {} bytes, too short to hold a header's magic and length",
                bytes.len()
            )));
        }
        if bytes[..HEAD_MAGIC.len()] != HEAD_MAGIC {
            return Err(malformed(format!(
                "does not start with {:?}; it is not an ng psp",
                String::from_utf8_lossy(&HEAD_MAGIC).trim_end()
            )));
        }

        let length_at = HEAD_MAGIC.len();
        let body_bytes = u64::from_le_bytes(
            bytes[length_at..length_at + 8]
                .try_into()
                .expect("the slice is eight bytes long"),
        );
        if body_bytes == 0 || body_bytes > MAX_HEADER_BODY_BYTES {
            return Err(malformed(format!(
                "declares a {body_bytes}-byte header body; the format allows 1 to \
                 {MAX_HEADER_BODY_BYTES}"
            )));
        }

        let body_at = length_at + 8;
        let sentinel_at = body_at + body_bytes as usize;
        let header_bytes = sentinel_at + HEAD_SENTINEL.len();
        if bytes.len() < header_bytes {
            return Err(malformed(format!(
                "declares a {body_bytes}-byte header body but only {} bytes follow the \
                 length",
                bytes.len() - body_at
            )));
        }
        if &bytes[sentinel_at..header_bytes] != HEAD_SENTINEL.as_slice() {
            return Err(malformed(
                "the header's declared length does not reach its closing line".to_string(),
            ));
        }

        let body = std::str::from_utf8(&bytes[body_at..sentinel_at])
            .map_err(|e| malformed(format!("the header body is not valid UTF-8: {e}")))?;

        let format_version = version_of(body, path)?;
        if format_version.0 != FORMAT_VERSION.0 {
            return Err(PspReadError::UnsupportedVersion {
                path: path.to_path_buf(),
                found: format_version,
                supported: FORMAT_VERSION,
            });
        }

        let wire: WireHeader = toml::from_str(body)
            .map_err(|e| malformed(format!("the header body is not valid TOML: {e}")))?;
        let header = wire
            .into_header(format_version)
            .map_err(|broken| malformed(format!("{}: {}", broken.field, broken.reason)))?;
        check_rules(&header)
            .map_err(|broken| malformed(format!("{}: {}", broken.field, broken.reason)))?;

        Ok((header, header_bytes))
    }
}

/// The format version alone, read without interpreting anything else in the body.
///
/// **This is what keeps the header plain text.** A reader has to be able to learn the
/// version of a file it cannot otherwise read, so the version is taken from a bare TOML
/// table before the body is deserialised into types this version's reader knows
/// (spec §3.1).
fn version_of(body: &str, path: &Path) -> Result<(u16, u16), PspReadError> {
    let malformed = |reason: String| PspReadError::MalformedHeader {
        path: path.to_path_buf(),
        reason,
    };

    let table: toml::Table = body
        .parse()
        .map_err(|e| malformed(format!("the header body is not valid TOML: {e}")))?;
    let spelled = table
        .get("format-version")
        .and_then(toml::Value::as_str)
        .ok_or_else(|| malformed("the header has no format-version".to_string()))?;
    let (major, minor) = spelled
        .split_once('.')
        .ok_or_else(|| malformed(format!("format-version {spelled:?} is not MAJOR.MINOR")))?;
    let parsed = major.parse::<u16>().ok().zip(minor.parse::<u16>().ok());
    parsed.ok_or_else(|| malformed(format!("format-version {spelled:?} is not MAJOR.MINOR")))
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
    if Path::new(spelled).components().count() > 1 {
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
    if manifest.block_byte_ceiling == Some(0) {
        return Err(BrokenRule::new(
            "manifest.block-byte-ceiling",
            "is zero; a ceiling no block can stay under closes every block empty",
        ));
    }
    // zstd's own bounds on `windowLog`: 2^10 is the smallest window it will take and 2^31
    // the largest. A file outside them cannot be decompressed by any reader, ours included.
    if !(10..=31).contains(&manifest.window_log) {
        return Err(BrokenRule::new(
            "manifest.window-log",
            format!(
                "is {}; zstd takes a look-back window between 2^10 and 2^31 bytes",
                manifest.window_log
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

fn check_encoding(name: &FieldName, encoding: FieldEncoding) -> Result<(), BrokenRule> {
    let field = "manifest.field.encoding";
    match encoding {
        FieldEncoding::Fixed { bytes } => {
            if !matches!(bytes, 1 | 2 | 4 | 8) {
                return Err(BrokenRule::new(
                    field,
                    format!(
                        "{:?} is a {bytes}-byte fixed integer; the widths are 1, 2, 4 and 8",
                        name.0
                    ),
                ));
            }
        }
        FieldEncoding::Ieee { bytes } => {
            if !matches!(bytes, 4 | 8) {
                return Err(BrokenRule::new(
                    field,
                    format!(
                        "{:?} is a {bytes}-byte IEEE float; the widths are 4 and 8",
                        name.0
                    ),
                ));
            }
        }
        FieldEncoding::FixedPoint { scale } => {
            if scale == 0 {
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
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
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
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
struct WireReference {
    name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    md5: Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
struct WireContig {
    name: String,
    length: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    md5: Option<String>,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
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
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
struct WireManifest {
    genomic_block_size_bp: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    block_byte_ceiling: Option<u32>,
    window_log: u8,
    #[serde(default)]
    field: Vec<WireFieldSpec>,
}

/// One field's declaration, **flat rather than nested**: `encoding` names the scheme and
/// `bytes` or `scale` carries its one parameter. A nested table per field would read as
/// `[manifest.field.encoding]` inside an array of tables, which is legal TOML and hard to
/// read in `head` — and readability is the reason the header is text at all.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case", deny_unknown_fields)]
struct WireFieldSpec {
    name: String,
    encoding: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    bytes: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    scale: Option<u32>,
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
                window_log: header.manifest.window_log,
                field: header
                    .manifest
                    .fields
                    .iter()
                    .map(|field| {
                        let (encoding, bytes, scale) = match field.encoding {
                            FieldEncoding::Varint => ("varint", None, None),
                            FieldEncoding::SignedVarint => ("signed-varint", None, None),
                            FieldEncoding::Fixed { bytes } => ("fixed", Some(bytes), None),
                            FieldEncoding::Ieee { bytes } => ("ieee", Some(bytes), None),
                            FieldEncoding::FixedPoint { scale } => {
                                ("fixed-point", None, Some(scale))
                            }
                            FieldEncoding::LengthPrefixedBytes => {
                                ("length-prefixed-bytes", None, None)
                            }
                        };
                        WireFieldSpec {
                            name: field.name.0.clone(),
                            encoding: encoding.to_string(),
                            bytes,
                            scale,
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
                Ok(ContigEntry {
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
                window_log: self.manifest.window_log,
                fields,
            },
        })
    }
}

/// The declared scheme and its one parameter, together.
///
/// **A parameter that belongs to another scheme is refused rather than ignored.** A file
/// saying `encoding = "varint"` beside `scale = 4096` was written by something that meant
/// one of the two, and reading it as a plain varint would silently divide every value in
/// that field by 4,096.
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

    let encoding = match field.encoding.as_str() {
        "varint" => FieldEncoding::Varint,
        "signed-varint" => FieldEncoding::SignedVarint,
        "length-prefixed-bytes" => FieldEncoding::LengthPrefixedBytes,
        "fixed" => FieldEncoding::Fixed {
            bytes: field.bytes.ok_or_else(|| missing("width"))?,
        },
        "ieee" => FieldEncoding::Ieee {
            bytes: field.bytes.ok_or_else(|| missing("width"))?,
        },
        "fixed-point" => FieldEncoding::FixedPoint {
            scale: field.scale.ok_or_else(|| missing("step"))?,
        },
        unknown => {
            return Err(BrokenRule::new(
                named,
                format!(
                    "{:?} is {unknown:?}, which is not one of varint, signed-varint, fixed, \
                     ieee, fixed-point, length-prefixed-bytes",
                    field.name
                ),
            ));
        }
    };
    match encoding {
        FieldEncoding::Fixed { .. } | FieldEncoding::Ieee { .. } => {
            if field.scale.is_some() {
                return Err(wrong_parameter("step"));
            }
        }
        FieldEncoding::FixedPoint { .. } => {
            if field.bytes.is_some() {
                return Err(wrong_parameter("width"));
            }
        }
        _ => {
            if field.bytes.is_some() {
                return Err(wrong_parameter("width"));
            }
            if field.scale.is_some() {
                return Err(wrong_parameter("step"));
            }
        }
    }
    Ok(encoding)
}
#[cfg(test)]
mod tests {
    use super::*;

    /// The header's contig list is what pins ng's coordinate space, so its equality has to
    /// be plain field equality. `crate::fasta::ContigEntry` — the same name, one module
    /// over — treats an absent MD5 as a wildcard, and a round-trip test written against
    /// *that* would pass while the encoder dropped every digest in the file.
    #[test]
    fn a_dropped_contig_md5_is_a_difference_here_and_a_wildcard_in_the_fasta_type() {
        let with_digest = ContigEntry {
            name: "SL4.0ch01".to_string(),
            length: 90_863_682,
            md5: Some([7u8; 16]),
        };
        let without = ContigEntry {
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
        let quarter_read = FieldEncoding::FixedPoint { scale: 4 };
        let sixteenth_read = FieldEncoding::FixedPoint { scale: 16 };
        assert_ne!(quarter_read, sixteenth_read);
        assert_ne!(quarter_read, FieldEncoding::Varint);
    }

    /// A header that says what a real tomato run says, minus the 11 contigs the shape of
    /// the test does not need.
    fn a_written_header() -> Header {
        Header {
            format_version: FORMAT_VERSION,
            sample: "SRR7279481".to_string(),
            reference: ReferenceIdentity {
                name: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                md5: Some([0x0a; 16]),
            },
            contigs: vec![
                ContigEntry {
                    name: "SL4.0ch00".to_string(),
                    length: 9_643_250,
                    md5: Some([0x1b; 16]),
                },
                ContigEntry {
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

    fn a_manifest() -> Manifest {
        Manifest {
            genomic_block_size_bp: Bp(100_000),
            block_byte_ceiling: Some(1_048_576),
            window_log: 15,
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
                    encoding: FieldEncoding::Fixed { bytes: 4 },
                },
                FieldSpec {
                    name: FieldName("allele-bases".to_string()),
                    encoding: FieldEncoding::LengthPrefixedBytes,
                },
                FieldSpec {
                    name: FieldName("window-mean-coverage".to_string()),
                    encoding: FieldEncoding::FixedPoint { scale: 4 },
                },
                FieldSpec {
                    name: FieldName("summed-log-error".to_string()),
                    encoding: FieldEncoding::FixedPoint { scale: 4_096 },
                },
            ],
        }
    }

    fn decoded(bytes: &[u8]) -> Result<(Header, usize), PspReadError> {
        Header::decode(bytes, Path::new("SRR7279481.psp"))
    }

    /// Everything in the header comes back, field for field. **Equality here is strict** —
    /// the contig row is this module's own type precisely so that a dropped MD5 is a
    /// difference and not a wildcard match.
    #[test]
    fn a_header_round_trips_field_for_field() {
        let written = a_written_header();
        let bytes = written.encode().expect("a valid header encodes");
        let (read_back, header_bytes) = decoded(&bytes).expect("its own bytes parse");
        assert_eq!(read_back, written);
        assert_eq!(
            header_bytes,
            bytes.len(),
            "the reported length is where the first block begins"
        );
    }

    /// The two fields that would round-trip through a wrong scale without complaining: a
    /// step of 1/4 of a read and one of 1/4,096 of a natural log sit in the same manifest,
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
                FieldEncoding::FixedPoint { scale } => Some((field.name.0.as_str(), scale)),
                _ => None,
            })
            .collect();
        assert_eq!(
            steps,
            vec![("window-mean-coverage", 4), ("summed-log-error", 4_096)]
        );
    }

    /// The reason the header is text: `head` on a psp tells you what it is. The magic ends
    /// in a newline so the body starts on its own line, and the body is TOML a person can
    /// read.
    #[test]
    fn the_body_is_readable_toml_after_a_newline_terminated_magic() {
        let bytes = a_written_header().encode().expect("a valid header encodes");
        assert_eq!(&bytes[..4], b"NGP\n");
        let body_bytes = u64::from_le_bytes(bytes[4..12].try_into().expect("eight bytes")) as usize;
        let body = std::str::from_utf8(&bytes[12..12 + body_bytes]).expect("UTF-8");
        assert!(
            body.contains("format-version = \"1.0\""),
            "body was: {body}"
        );
        assert!(body.contains("sample = \"SRR7279481\""), "body was: {body}");
        assert!(
            body.contains("genomic-block-size-bp = 100000"),
            "body was: {body}"
        );
        assert!(
            body.contains("encoding = \"fixed-point\""),
            "body was: {body}"
        );
        assert!(body.contains("scale = 4096"), "body was: {body}");
        assert_eq!(&bytes[12 + body_bytes..], b"---END-HEADER---\n");
    }

    /// A parameter is written as the value itself, not as a tagged pair, because the header
    /// is meant to be read by eye.
    #[test]
    fn a_parameter_is_written_as_its_bare_value() {
        let bytes = a_written_header().encode().expect("a valid header encodes");
        let body = String::from_utf8(bytes).expect("UTF-8");
        assert!(body.contains("depth-cap = 300"), "body was: {body}");
        assert!(body.contains("realign = true"), "body was: {body}");
    }

    /// Goal 5 is that the same sample gathered at any worker count gives the same bytes, so
    /// nothing in the header may depend on the order a map happened to be filled in.
    #[test]
    fn the_same_header_encodes_to_the_same_bytes_whatever_order_it_was_built_in() {
        let mut backwards = a_written_header();
        let parameters: Vec<_> = backwards.writer.parameters.clone().into_iter().collect();
        backwards.writer.parameters = parameters.into_iter().rev().collect();
        assert_eq!(
            backwards.encode().expect("encodes"),
            a_written_header().encode().expect("encodes")
        );
    }

    /// A file from a later format has to say so. The version is read from a bare TOML table
    /// before the body is deserialised, so a body full of keys this reader has never seen
    /// still yields the right answer instead of a parse failure.
    #[test]
    fn a_newer_major_version_is_refused_as_a_version_and_not_as_damage() {
        let body = "format-version = \"2.0\"\n\
                    whatever-version-two-added = { shape = \"nothing here knows\" }\n";
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&(body.len() as u64).to_le_bytes());
        bytes.extend_from_slice(body.as_bytes());
        bytes.extend_from_slice(HEAD_SENTINEL);

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

    /// A later *minor* of the same major is this reader's to read: that is what the split
    /// into major and minor is for.
    #[test]
    fn a_newer_minor_version_is_not_refused_by_the_version_check() {
        let mut written = a_written_header();
        let bytes = written.encode().expect("a valid header encodes");
        let body_bytes = u64::from_le_bytes(bytes[4..12].try_into().expect("eight bytes")) as usize;
        let body = std::str::from_utf8(&bytes[12..12 + body_bytes])
            .expect("UTF-8")
            .replace("format-version = \"1.0\"", "format-version = \"1.7\"");
        let mut later_minor = Vec::new();
        later_minor.extend_from_slice(&HEAD_MAGIC);
        later_minor.extend_from_slice(&(body.len() as u64).to_le_bytes());
        later_minor.extend_from_slice(body.as_bytes());
        later_minor.extend_from_slice(HEAD_SENTINEL);

        written.format_version = (1, 7);
        let (read_back, _) = decoded(&later_minor).expect("a later minor of this major reads");
        assert_eq!(read_back, written);
    }

    /// Both formats use the extension `.psp`, so the first four bytes are what tells them
    /// apart. A production `.psp` handed to this reader must be refused here rather than
    /// misread further in.
    #[test]
    fn a_production_psp_is_refused_at_the_magic() {
        let mut productions = Vec::new();
        productions.extend_from_slice(&crate::psp::header::HEAD_MAGIC);
        productions.extend_from_slice(&(1u64).to_le_bytes());
        productions.push(b'x');
        productions.extend_from_slice(HEAD_SENTINEL);

        let refused = decoded(&productions).expect_err("production's magic is not ours");
        assert!(
            matches!(refused, PspReadError::MalformedHeader { .. }),
            "got {refused:?}"
        );
        assert!(
            refused.to_string().contains("not an ng psp"),
            "got {refused}"
        );
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

    /// The length field is read before anything is allocated, so a corrupt one cannot ask
    /// for a large buffer on its own say-so.
    #[test]
    fn an_enormous_declared_length_is_refused_before_anything_is_allocated() {
        let mut bytes = a_written_header().encode().expect("a valid header encodes");
        bytes[4..12].copy_from_slice(&u64::MAX.to_le_bytes());

        let refused = decoded(&bytes).expect_err("a header cannot be that big");
        assert!(
            refused.to_string().contains("the format allows 1 to"),
            "got {refused}"
        );
    }

    /// Each rule is written once and checked on both sides, so the writer refuses to make a
    /// file the reader would refuse to read. The table names the rule, a header that breaks
    /// it, and the words the message has to carry.
    #[test]
    fn every_rule_is_refused_by_the_writer_and_by_the_reader_alike() {
        let broken: Vec<(&str, Box<dyn Fn(&mut Header)>, &str)> = vec![
            (
                "an empty sample name",
                Box::new(|header| header.sample = "  ".to_string()),
                "sample",
            ),
            (
                "a reference recorded with its directory",
                Box::new(|header| {
                    header.reference.name = "/home/jose/genomes/tomato.fa".to_string()
                }),
                "directory component",
            ),
            (
                "no contigs",
                Box::new(|header| header.contigs.clear()),
                "contig",
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
                "a byte ceiling no block can stay under",
                Box::new(|header| header.manifest.block_byte_ceiling = Some(0)),
                "closes every block empty",
            ),
            (
                "a look-back window zstd cannot take",
                Box::new(|header| header.manifest.window_log = 40),
                "2^10 and 2^31",
            ),
            (
                "no declared fields",
                Box::new(|header| header.manifest.fields.clear()),
                "there are none",
            ),
            (
                "one field declared twice",
                Box::new(|header| {
                    header.manifest.fields[1].name = header.manifest.fields[0].name.clone()
                }),
                "appears twice",
            ),
            (
                "a fixed integer of a width that is not a power of two",
                Box::new(|header| {
                    header.manifest.fields[2].encoding = FieldEncoding::Fixed { bytes: 3 }
                }),
                "1, 2, 4 and 8",
            ),
            (
                "an IEEE float of no legal width",
                Box::new(|header| {
                    header.manifest.fields[2].encoding = FieldEncoding::Ieee { bytes: 2 }
                }),
                "4 and 8",
            ),
            (
                "a fixed-point field counting steps of nothing",
                Box::new(|header| {
                    header.manifest.fields[4].encoding = FieldEncoding::FixedPoint { scale: 0 }
                }),
                "1/0",
            ),
        ];

        for (what, break_it, expected) in broken {
            let mut header = a_written_header();
            break_it(&mut header);

            let refused = header
                .encode()
                .expect_err(&format!("the writer must refuse {what}"));
            assert!(
                refused.to_string().contains(expected),
                "the writer's message for {what} was {refused:?}, which does not say \
                 {expected:?}"
            );

            // The same rule, met from the other side: a file that already holds the broken
            // value has to be refused when it is read, not only when it is written.
            let smuggled = smuggle(&header);
            let refused = decoded(&smuggled).expect_err(&format!("the reader must refuse {what}"));
            assert!(
                refused.to_string().contains(expected),
                "the reader's message for {what} was {refused:?}, which does not say \
                 {expected:?}"
            );
        }
    }

    /// Frame a header's TOML **without** running the writer's rules over it, so a file
    /// holding a value the writer would never produce can be handed to the reader.
    ///
    /// A NaN parameter has no TOML spelling at all, so it is written as the string `nan`
    /// and the reader meets it as a bad value rather than a bad float; every other broken
    /// header here serialises as itself.
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
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&(body.len() as u64).to_le_bytes());
        bytes.extend_from_slice(body.as_bytes());
        bytes.extend_from_slice(HEAD_SENTINEL);
        bytes
    }

    /// A field carrying a parameter that belongs to another scheme was written by something
    /// that meant one of the two. Reading it as the scheme it names would silently divide
    /// every value in that field by the step it was not supposed to have.
    #[test]
    fn an_encoding_carrying_the_wrong_parameter_is_refused() {
        let body = one_field_declared_as(
            "name = \"summed-log-error\"\nencoding = \"varint\"\nscale = 4096\n",
        );
        let refused = decoded(&body).expect_err("a varint has no step");
        assert!(
            refused.to_string().contains("has no use for"),
            "got {refused}"
        );
    }

    /// A width or a step that is simply missing is the same class of error as one that does
    /// not belong: the file does not say how to read the field.
    #[test]
    fn a_fixed_width_field_without_its_width_is_refused() {
        let body = one_field_declared_as("name = \"body-bytes\"\nencoding = \"fixed\"\n");
        let refused = decoded(&body).expect_err("a fixed integer needs a width");
        assert!(
            refused.to_string().contains("carries no width"),
            "got {refused}"
        );
    }

    /// An encoding this reader has never heard of is refused rather than guessed at. The
    /// message lists what it does know, because that is what tells the reader whether to
    /// upgrade or to rebuild.
    #[test]
    fn an_unknown_encoding_is_refused_and_the_message_lists_the_known_ones() {
        let body = one_field_declared_as("name = \"chain-ids\"\nencoding = \"roaring-bitmap\"\n");
        let refused = decoded(&body).expect_err("that is not one of the six");
        assert!(
            refused.to_string().contains("roaring-bitmap"),
            "got {refused}"
        );
        assert!(
            refused.to_string().contains("length-prefixed-bytes"),
            "got {refused}"
        );
    }

    /// Frame a header whose manifest holds exactly the one field declaration given, so a
    /// declaration that no `FieldEncoding` can represent can still be put in a file.
    fn one_field_declared_as(declaration: &str) -> Vec<u8> {
        let whole = a_written_header();
        let body = toml::to_string_pretty(&WireHeader::from(&whole)).expect("encodes");
        let up_to_the_fields = body
            .split("[[manifest.field]]")
            .next()
            .expect("the manifest declares fields")
            .to_string();
        let body = format!("{up_to_the_fields}[[manifest.field]]\n{declaration}");

        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&(body.len() as u64).to_le_bytes());
        bytes.extend_from_slice(body.as_bytes());
        bytes.extend_from_slice(HEAD_SENTINEL);
        bytes
    }

    /// An MD5 is 32 lowercase hex characters, which is what a SAM `@SQ M5` is. A digest of
    /// the wrong length is not a shorter digest; it is a different file's.
    #[test]
    fn a_contig_digest_that_is_not_thirty_two_hex_characters_is_refused() {
        let whole = a_written_header();
        let body = toml::to_string_pretty(&WireHeader::from(&whole))
            .expect("encodes")
            .replace(&hex_of([0x1b; 16]), "1b1b1b");
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&HEAD_MAGIC);
        bytes.extend_from_slice(&(body.len() as u64).to_le_bytes());
        bytes.extend_from_slice(body.as_bytes());
        bytes.extend_from_slice(HEAD_SENTINEL);

        let refused = decoded(&bytes).expect_err("that is not a digest");
        assert!(
            refused.to_string().contains("32 lowercase hex"),
            "got {refused}"
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
}
