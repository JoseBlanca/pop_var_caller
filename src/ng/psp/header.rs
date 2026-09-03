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
use crate::ng::psp::segmentation_section::{self, WireSegmentation};
use crate::ng::segmentation_inputs::SegmentationInputs;
use crate::ng::types::{Bp, ReadGroupId};

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
    /// The sample's read groups, **in the order that defines the file's
    /// [`ReadGroupId`]s** — entry `i` is the read group this walk numbered `i`, and every
    /// observation in the file carries these walk-local numbers.
    ///
    /// **Without this table, no cohort can be assembled from separately-walked samples
    /// at all** (spec `run_streaming.md` §6.1): a gatherer sees one sample's files, so it
    /// numbers that sample's read groups from zero, and every sample's first read group
    /// comes back as identifier 0. The calling stage reads every file's table at open and
    /// merges them into one run-wide numbering (§6.2) — which is what lets a sample be
    /// walked once and joined to any cohort later.
    pub read_groups: Vec<ReadGroupIdentity>,
    /// The widest reference span any observation in this file can have, in bases — the
    /// locus generator's own cap
    /// ([`PileupGeneratorConfig::max_record_span`](crate::ng::locus_generation::pileup::PileupGeneratorConfig::max_record_span)),
    /// known before the first record because it is a setting, not a measurement.
    ///
    /// **A sizing fact, not a correctness fact** (`psp_file_format.md` §3.1): a cohort
    /// reader can size its observation cache up front, taking the maximum over its
    /// files, instead of growing it — and a forward reader never needs it at all, so
    /// nothing refuses on it. `cohort_merge.md` §13 is the consumer; it reads this at
    /// open (plan step E4).
    pub observation_reach_ceiling_bp: Bp,
    /// What produced the file and what it ran with.
    pub writer: WriterProvenance,
    /// The ground the sample was analysed over, and what the segmentation that shaped
    /// its observations was computed from — the operands of a calling run's
    /// compatibility checks (spec `run_streaming.md` §6.1, §6.2). The wire shape is
    /// [`crate::ng::psp::segmentation_section`]'s.
    pub segmentation_inputs: SegmentationInputs,
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

/// One read group as the header records it: what the alignment file called it, which
/// library its reads came from, and the number this walk gave it.
///
/// The three things spec `run_streaming.md` §6.1 asks the table to carry, and no more:
/// the sample is the header's own `sample` field (one psp, one sample), and the file the
/// group was declared in is a path, which the header records only as provenance
/// basenames.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroupIdentity {
    /// The `@RG ID`, verbatim. **A label, never an identity** — the SAM specification
    /// makes it unique within its file and says nothing across files, so two entries of
    /// this table may share one `id` when a sample was sequenced across files
    /// (`read_groups.md` §4). Identity is [`walk_local_id`](Self::walk_local_id).
    pub id: String,
    /// The library the group's reads came from — `@RG LB`, or the name the walk
    /// synthesized when the file declared none. The parameters fit keys its per-library
    /// error rates on groupings of this.
    pub library: String,
    /// The number this walk gave the group — **walk-local**, as the name says: every
    /// psp numbers its own read groups from zero, so these collide across files by
    /// construction and the calling stage renumbers them at open (spec §6.2).
    ///
    /// **Also this entry's position in the table** — entry `i` carries number `i`, a
    /// rule checked on both sides — so the number a person reads beside an `@RG ID`
    /// and the number the code derives from order cannot disagree.
    pub walk_local_id: ReadGroupId,
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

impl WriterProvenance {
    /// Put a batch of settings into [`parameters`](Self::parameters), each under its own
    /// key — the seam through which a walk records what it ran with (spec
    /// `run_streaming.md` §6.1: recorded, never compared).
    ///
    /// **Insert-only**: a key already present is overwritten, a key not in `entries` is
    /// left standing. A caller replacing an earlier recording whose key set may have
    /// shrunk — the read filters' conditional floor, say — clears its own keys first
    /// (the producer of the entries names them, e.g.
    /// [`READ_FILTER_PROVENANCE_KEYS`](crate::ng::read::READ_FILTER_PROVENANCE_KEYS)).
    pub fn record_parameters(
        &mut self,
        entries: impl IntoIterator<Item = (String, ParameterValue)>,
    ) {
        self.parameters.extend(entries);
    }
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

/// A field of a record, and how to read it: what it is called, what one appearance of it
/// looks like, and how each one is laid down.
///
/// **The three things spec §4.5 asks the manifest to carry for every field.** The middle one
/// is [`FieldShape`], and it is not a field of this struct: an encoding fixes it, so it is
/// derived rather than stored beside the thing that decides it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldSpec {
    pub name: FieldName,
    pub encoding: FieldEncoding,
}

impl FieldSpec {
    /// What one appearance of this field looks like on the wire.
    ///
    /// **Derived from the encoding rather than stored beside it**, which is what keeps the
    /// manifest's two accounts of a field's shape from being able to disagree: in memory
    /// there is one account, and a *file's* own `shape` key is checked against it as the file
    /// is read ([`check_declared_shape`]). So a `FieldSpec` that came out of
    /// [`Header::decode`] says what its file said, and one this build assembled cannot say
    /// anything else.
    pub fn shape(&self) -> FieldShape {
        self.encoding.shape()
    }
}

/// Generates [`FieldShape`], the list of every shape, and each one's spelling **from one
/// source**, so a shape cannot be added to the type without reaching the list and the
/// spelling too.
///
/// This is the drift `ALL_ENCODINGS`' own doc comment records having been paid for once
/// already: when the writer's spellings and the reader's were two hand-maintained lists, a
/// scheme added to one alone made [`Header::encode`] write a file [`Header::decode`] refused,
/// and the whole suite stayed green. A hand-written array cannot prevent that, because its
/// length is a number in the source — a new variant leaves it short and still compiling.
macro_rules! field_shapes {
    (
        $(#[$enum_meta:meta])*
        pub enum $name:ident {
            $( $(#[$variant_meta:meta])* $variant:ident => $spelling:literal ),+ $(,)?
        }
    ) => {
        $(#[$enum_meta])*
        #[non_exhaustive]
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub enum $name {
            $( $(#[$variant_meta])* $variant ),+
        }

        /// Every shape, listed once — **generated beside the enum**, so the writer spells from
        /// the list the reader recognises from and neither can go short of the other.
        const ALL_SHAPES: &[$name] = &[ $( $name::$variant ),+ ];

        impl $name {
            /// What this shape is called in the header.
            fn spelled(self) -> &'static str {
                match self {
                    $( $name::$variant => $spelling ),+
                }
            }
        }
    };
}

field_shapes! {
    /// What **one appearance of a field** looks like on the wire: a single value, or a count
    /// followed by that many values.
    ///
    /// **What it buys, and it is why spec §4.5 asks the manifest to carry it**: a reader
    /// meeting a field name it has never heard of can still say where that field ends and step
    /// over it, so a psp written by a later version of ng stays readable by an older one. A
    /// scalar ends after one value; a list ends once the count in front of it has been
    /// honoured.
    ///
    /// **⚠ It does not say how often the field appears in a record**, which is the other half
    /// of the idea and the half production's store calls *cardinality*
    /// (`Cardinality { PerRecord, PerAllele }`, `src/psp/registry.rs`). The two must not be
    /// confused, and the wire key here is `shape` rather than `cardinality` for exactly that
    /// reason: both formats use the extension `.psp` and will sit on the same disks, so the
    /// same key name meaning two different things is a trap for whoever reads one with `head`.
    /// This type is production's `Shape`, and it is named and spelled to match it.
    ///
    /// Concretely: `mapq-sum` is a plain integer, so its shape is [`Scalar`](Self::Scalar) —
    /// yet a record with five observations holds five of them, because
    /// `encode_record_body` writes it once per observation. Ten of the body's twenty-two
    /// fields repeat that way and two more repeat once per witness run; **nothing in the
    /// manifest says so**, and a record's own counts are what a reader uses. So a later writer
    /// that appends a field repeating per observation writes a file this reader cannot step
    /// over, and must raise the format version instead. The two fields waiting to be
    /// added — a window's GC fraction and its mean coverage — appear once per record.
    ///
    /// **Every encoding fixes exactly one of these**, and [`FieldEncoding::shape`] is where
    /// that is written down; a file declaring anything else is refused
    /// ([`check_declared_shape`]).
    pub enum FieldShape {
        /// One value, whose bytes the encoding alone measures.
        ///
        /// **A length-prefixed byte string is one value**, not a list: its prefix counts its
        /// *bytes*, not how many appearances of the field there are. An allele's sequence is
        /// one sequence however long it is.
        Scalar => "scalar",
        /// A count, then that many values of the field's element type — an observation's read
        /// identifiers, or the arrivals and departures of a record's live set.
        List => "list",
    }
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
    /// **Which reads started covering this record and which stopped**: a count of departures,
    /// each a position in the live set as the previous record left it; then a count of arrivals,
    /// each an identifier. Both runs are strictly ascending and written as their gaps
    /// (spec `psp_chain_id_encoding.md` §4).
    ///
    /// **The one composite in an otherwise scalar set, and it earns the exception.** The others
    /// exist so a reader can measure a field it does not recognise and walk past it. This one can
    /// be measured the same way — two counted runs of varints — but **walking past the copy in a
    /// record's head would be wrong**, because that copy carries the state every later record is
    /// decoded against, and a reader that stepped over it would build all of them against a stale
    /// set. So the head's copy is a field this reader knows by name and refuses a file that moves
    /// or renames; what the measure-and-step-over path exists for is a *later* writer putting
    /// another one at the end of a body, which nothing today does.
    ChainIdChanges,
    /// **One observation's own reads**: a count, then the identifiers as ascending gaps.
    ///
    /// A different shape from [`ChainIdChanges`](Self::ChainIdChanges), which is *two* counted
    /// runs — so the two cannot share a scheme even though both are runs of identifiers. ⚠ They
    /// briefly did, and a reader stepping over a field of one under the other's rule measures the
    /// wrong number of bytes: `[5, 9]` is three bytes as a list and seven as a set of changes.
    ///
    /// **This one may be stepped over**, unlike the changes, because it carries no state: it is
    /// the half of the chain-id column that lives in a record's skippable body.
    ChainIdList,
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
pub(crate) const MAX_TOML_INTEGER: u64 = i64::MAX as u64;

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
///
/// `pub(crate)` because the segmentation section of the header lives in its own file
/// ([`crate::ng::psp::segmentation_section`]) and reports its broken rules in the same
/// shape.
#[derive(Debug)]
pub(crate) struct BrokenRule {
    pub(crate) field: String,
    pub(crate) reason: String,
}

impl BrokenRule {
    pub(crate) fn new(field: impl Into<String>, reason: impl Into<String>) -> Self {
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
    if header.observation_reach_ceiling_bp.get() == 0 {
        return Err(BrokenRule::new(
            "observation-reach-ceiling-bp",
            "is zero; an observation covers at least one base, so no record could exist \
             under this ceiling",
        ));
    }
    if header.observation_reach_ceiling_bp.get() > MAX_TOML_INTEGER {
        return Err(BrokenRule::new(
            "observation-reach-ceiling-bp",
            format!(
                "is {}; a TOML integer is signed, so a header cannot carry more than \
                 {MAX_TOML_INTEGER}",
                header.observation_reach_ceiling_bp.get()
            ),
        ));
    }
    check_contigs(&header.contigs)?;
    check_read_groups(&header.read_groups)?;
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

    segmentation_section::check_segmentation(&header.segmentation_inputs, &header.contigs)?;
    check_manifest(&header.manifest)
}

/// The contig-list rules, split out of [`check_rules`] because the reader needs them
/// **before** the rest: the segmentation section anchors its analysed spans to this
/// list by name, so a duplicated or broken contig has to be refused as what it is,
/// not as the span-resolution failure it would cause two lines later.
fn check_contigs(contigs: &[ContigIdentity]) -> Result<(), BrokenRule> {
    if contigs.is_empty() {
        return Err(BrokenRule::new(
            "contig",
            "is empty; a psp's coordinates mean nothing without the contig list they index",
        ));
    }
    let mut seen = std::collections::HashSet::with_capacity(contigs.len());
    for contig in contigs {
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
    Ok(())
}

/// The read-group table's rules: the table exists, its identifiers are the walk's own
/// numbering from zero in order, and its strings cannot forge a line in the header's
/// text.
///
/// **A duplicated `@RG ID` is legal here**, deliberately: SAM makes the id unique only
/// within one file, so a sample sequenced across files may carry two entries with one
/// `id` and different libraries. What *calling* refuses is a table it cannot merge —
/// that refusal is the calling stage's (spec §6.2), because it is about assembling a
/// cohort, not about this file being well-formed.
fn check_read_groups(read_groups: &[ReadGroupIdentity]) -> Result<(), BrokenRule> {
    if read_groups.is_empty() {
        return Err(BrokenRule::new(
            "read-group",
            "is empty; a sample's reads come from at least one read group, and a table \
             with none cannot be renumbered at calling",
        ));
    }
    for (position, group) in read_groups.iter().enumerate() {
        // A table of more than u32::MAX read groups cannot exist: each entry costs
        // header bytes and the body ceiling is reached long before.
        if group.walk_local_id.get() != position as u32 {
            return Err(BrokenRule::new(
                "read-group.walk-local-id",
                format!(
                    "entry {position} carries identifier {}; the identifiers are the \
                     walk's own numbering from zero, in table order",
                    group.walk_local_id.get()
                ),
            ));
        }
        if group.id.is_empty() {
            return Err(BrokenRule::new("read-group.id", "is empty"));
        }
        // Control characters only — not all whitespace, because SAM allows a space in
        // an @RG ID and a library name, and refusing one would refuse real archives.
        // A space cannot forge a header line; a newline or control character can.
        check_no_control_characters("read-group.id", &group.id)?;
        if group.library.is_empty() {
            return Err(BrokenRule::new(
                "read-group.library",
                format!("of read group {:?} is empty", group.id),
            ));
        }
        check_no_control_characters("read-group.library", &group.library)?;
    }
    Ok(())
}

/// A string that may legitimately hold spaces (SAM allows them in `@RG` values) but
/// must not hold what could rewrite the header's text: a newline lands in the body as
/// further lines, and any control character is nothing an alignment header can carry.
fn check_no_control_characters(field: &str, value: &str) -> Result<(), BrokenRule> {
    if value.chars().any(char::is_control) {
        return Err(BrokenRule::new(
            field.to_string(),
            format!(
                "{value:?} holds a control character; an alignment header cannot carry \
                 one, and a newline would land in this header's text as lines no field \
                 declares"
            ),
        ));
    }
    Ok(())
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
        // **A name carrying a newline rewrites what the header appears to say.** TOML writes
        // such a name as a multi-line string, so its own bytes land in the body as further
        // lines — and a name holding `\ncardinality = "..."` shows a reader running `head` a
        // key that no field declared, while the file still parses and round-trips. Readability
        // is the reason this header is text at all, so a name that can forge a line in it is
        // refused. The contig-name rule a hundred lines above says the same thing for the same
        // reason; field names had no such rule until a header was built that exploited it.
        if field
            .name
            .0
            .chars()
            .any(|character| character.is_whitespace() || character.is_control())
        {
            return Err(BrokenRule::new(
                "manifest.field.name",
                format!(
                    "{:?} holds whitespace or a control character; such a name is written as a \
                     multi-line string and can make the header's text show a key no field \
                     declares",
                    field.name.0
                ),
            ));
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
        | FieldEncoding::LengthPrefixedBytes
        | FieldEncoding::ChainIdChanges
        | FieldEncoding::ChainIdList => {}
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

/// An `f32` widened for TOML through its own shortest decimal, so the header shows
/// `0.93` rather than `0.9300000071525574`.
///
/// Exact both ways: `Display` on an `f32` prints the shortest decimal that reads back
/// to the same `f32`, and that decimal's nearest `f64` narrows back to it.
///
// PANIC-FREE: `f64`'s parser accepts every string `f32`'s `Display` produces — `NaN`
// and `inf` included — so the expect cannot fire on any input, even the test-only path
// that serialises a rule-breaking header.
pub(crate) fn wire_float_of(value: f32) -> f64 {
    format!("{value}")
        .parse()
        .expect("a float's own Display re-parses")
}

/// Every MD5 travels as 32 lowercase hex characters, which is what a SAM `@SQ M5` is.
pub(crate) fn hex_of(digest: [u8; 16]) -> String {
    digest.iter().map(|byte| format!("{byte:02x}")).collect()
}

pub(crate) fn digest_of(field: &str, spelled: &str) -> Result<[u8; 16], BrokenRule> {
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
    // A bare value, so it must sit with the other top-level scalars, before the first
    // table opens.
    observation_reach_ceiling_bp: u64,
    reference: WireReference,
    #[serde(default)]
    contig: Vec<WireContig>,
    #[serde(default)]
    read_group: Vec<WireReadGroup>,
    writer: WireWriter,
    // Before the manifest, deliberately: the manifest's field declarations close the
    // body, so tests (and people) can cut a body at `[[manifest.field]]` and keep every
    // other section intact.
    segmentation: WireSegmentation,
    manifest: WireManifest,
}

/// One `[[read-group]]` row: the `@RG ID`, the library, and the walk-local number.
///
/// `walk-local-id` is also the row's position — redundancy checked on both sides, like
/// the header's declared length against its sentinel — so the number a person reads
/// beside an id and the number the code derives from order cannot disagree.
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
struct WireReadGroup {
    id: String,
    library: String,
    walk_local_id: u32,
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

/// One field's declaration, **flat rather than nested**: `shape` says whether one appearance
/// of the field is a single value or a counted run of them, `encoding` names the scheme and
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
    shape: String,
    encoding: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    width_bytes: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    steps_per_unit: Option<u32>,
}

impl From<&Header> for WireHeader {
    fn from(header: &Header) -> Self {
        // Exhaustive, so a field added to `Header` in a later step fails to compile
        // here, on the encode side — not only in the decode literal, where a default
        // could be filled in without the writer ever recording the field.
        let Header {
            format_version,
            sample,
            reference,
            contigs,
            read_groups,
            observation_reach_ceiling_bp,
            segmentation_inputs,
            writer,
            manifest,
        } = header;
        WireHeader {
            format_version: format!("{}.{}", format_version.0, format_version.1),
            sample: sample.clone(),
            observation_reach_ceiling_bp: observation_reach_ceiling_bp.get(),
            reference: WireReference {
                name: reference.name.clone(),
                md5: reference.md5.map(hex_of),
            },
            contig: contigs
                .iter()
                .map(|contig| WireContig {
                    name: contig.name.clone(),
                    length: contig.length,
                    md5: contig.md5.map(hex_of),
                })
                .collect(),
            read_group: read_groups
                .iter()
                .map(|group| WireReadGroup {
                    id: group.id.clone(),
                    library: group.library.clone(),
                    walk_local_id: group.walk_local_id.get(),
                })
                .collect(),
            writer: WireWriter {
                tool: writer.tool.clone(),
                version: writer.version.clone(),
                subcommand: writer.subcommand.clone(),
                input_alignments: writer.input_alignments.clone(),
                input_reference: writer.input_reference.clone(),
                command_line: writer.command_line.clone(),
                created: writer.created,
                parameters: writer.parameters.clone(),
            },
            manifest: WireManifest {
                genomic_block_size_bp: manifest.genomic_block_size_bp.get(),
                block_byte_ceiling: manifest.block_byte_ceiling,
                look_back_window_log: manifest.look_back_window_log,
                field: manifest
                    .fields
                    .iter()
                    .map(|field| {
                        let (encoding, width_bytes, steps_per_unit) = field.encoding.spelled();
                        WireFieldSpec {
                            name: field.name.0.clone(),
                            shape: field.shape().spelled().to_string(),
                            encoding: encoding.to_string(),
                            width_bytes,
                            steps_per_unit,
                        }
                    })
                    .collect(),
            },
            segmentation: WireSegmentation::from_inputs(segmentation_inputs, contigs),
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
                let encoding = encoding_of(&field)?;
                check_declared_shape(&field, encoding)?;
                Ok(FieldSpec {
                    encoding,
                    name: FieldName(field.name),
                })
            })
            .collect::<Result<Vec<_>, BrokenRule>>()?;
        let read_groups = self
            .read_group
            .into_iter()
            .map(|group| ReadGroupIdentity {
                id: group.id,
                library: group.library,
                walk_local_id: ReadGroupId(group.walk_local_id),
            })
            .collect();
        // The contig rules run before the segmentation section is resolved, because the
        // section anchors each analysed span to this list by name — resolving against a
        // duplicated or zero-length contig would report the span as broken when the
        // contig is. `check_rules` runs the same checks again afterwards; one function,
        // called twice, cannot disagree with itself.
        check_contigs(&contigs)?;
        let segmentation_inputs = self.segmentation.into_inputs(&contigs)?;

        Ok(Header {
            format_version,
            sample: self.sample,
            observation_reach_ceiling_bp: Bp(self.observation_reach_ceiling_bp),
            reference,
            contigs,
            read_groups,
            segmentation_inputs,
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
const ALL_ENCODINGS: [FieldEncoding; 8] = [
    FieldEncoding::Varint,
    FieldEncoding::SignedVarint,
    FieldEncoding::FixedWidthInteger { width_bytes: 1 },
    FieldEncoding::IeeeFloat { width_bytes: 4 },
    FieldEncoding::FixedPoint { steps_per_unit: 1 },
    FieldEncoding::LengthPrefixedBytes,
    FieldEncoding::ChainIdChanges,
    FieldEncoding::ChainIdList,
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
            FieldEncoding::ChainIdChanges => ("chain-id-changes", None, None),
            FieldEncoding::ChainIdList => ("chain-id-list", None, None),
        }
    }

    /// What one appearance of a field written this way looks like.
    ///
    /// **The one place the pair is decided**, which is what keeps the manifest's two accounts
    /// of a field's shape from being able to disagree: a file declaring anything else is
    /// refused (see [`FieldShape`]).
    ///
    /// **Exhaustive with no wildcard**, so a scheme added to the closed set has to say here
    /// what its bytes look like rather than inheriting *one value* in silence — and inheriting
    /// *one value* for a counted run of identifiers is exactly the mistake that would let a
    /// reader measure the wrong number of bytes.
    ///
    /// **⚠ Nothing else in this module can tell you whether an answer here is right.** Every
    /// round-trip test asks this function what a scheme's shape is and then checks the file
    /// agrees, so a wrong answer applied consistently passes all of them;
    /// `every_encoding_lays_down_the_shape_its_bytes_have` is the one test that states each
    /// answer independently, and it exists because three of these eight schemes are used by no
    /// record field today and nothing else reaches them.
    pub fn shape(self) -> FieldShape {
        match self {
            // A number, a fixed-width integer, a float, a count of fixed-point steps: one value
            // each. **A length-prefixed byte string is one value too** — its prefix counts its
            // bytes, not how many appearances of the field there are.
            FieldEncoding::Varint
            | FieldEncoding::SignedVarint
            | FieldEncoding::FixedWidthInteger { .. }
            | FieldEncoding::IeeeFloat { .. }
            | FieldEncoding::FixedPoint { .. }
            | FieldEncoding::LengthPrefixedBytes => FieldShape::Scalar,
            // A count and then that many identifiers; the changes are two such runs, which is
            // still a counted run rather than one value.
            FieldEncoding::ChainIdChanges | FieldEncoding::ChainIdList => FieldShape::List,
        }
    }
}

/// The declared shape, as a name from the closed set.
///
/// **An unrecognised spelling is refused rather than guessed at**, and the match is exact —
/// no trimming, no case folding. A reader that fell back on the encoding's own shape would
/// accept `shape = "per-observation"` from a later writer, which names a field it cannot in
/// fact step over, and would then read the record after it from a position that is not a field
/// boundary.
fn shape_of(field: &WireFieldSpec) -> Result<FieldShape, BrokenRule> {
    ALL_SHAPES
        .iter()
        .copied()
        .find(|candidate| candidate.spelled() == field.shape)
        .ok_or_else(|| {
            let known: Vec<&str> = ALL_SHAPES
                .iter()
                .map(|candidate| candidate.spelled())
                .collect();
            BrokenRule::new(
                "manifest.field.shape",
                format!(
                    "field {:?} declares shape {:?}, which is not one of {}",
                    field.name,
                    field.shape,
                    known.join(", ")
                ),
            )
        })
}

/// A file's declared shape must be the one its encoding lays down.
///
/// **The rule that stops the manifest's two accounts of a field's shape from disagreeing**,
/// and it lives on the reading side because that is the only side where two accounts exist: in
/// memory a [`FieldSpec`] derives its shape from its encoding, so this build cannot write a
/// disagreement in the first place.
///
/// A file saying `shape = "scalar"` beside `encoding = "chain-id-list"` was written by
/// something that meant one of the two, and a reader that believed the shape would step over a
/// count and stop inside a run of identifiers — landing in the middle of the next field rather
/// than failing.
fn check_declared_shape(field: &WireFieldSpec, encoding: FieldEncoding) -> Result<(), BrokenRule> {
    let declared = shape_of(field)?;
    let laid_down = encoding.shape();
    if declared != laid_down {
        return Err(BrokenRule::new(
            "manifest.field.shape",
            format!(
                "field {:?} declares shape {:?}, but its encoding {:?} lays down {:?}",
                field.name,
                declared.spelled(),
                encoding.spelled().0,
                laid_down.spelled(),
            ),
        ));
    }
    Ok(())
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
        FieldEncoding::ChainIdChanges => FieldEncoding::ChainIdChanges,
        FieldEncoding::ChainIdList => FieldEncoding::ChainIdList,
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
/// A read-group table for tests: two groups, two libraries, in walk order.
///
/// Two rather than one, so a round trip that dropped or reordered rows fails; the
/// values are non-defaults with the id and the library visibly different strings, so a
/// decode that read one column into the other cannot pass.
#[cfg(test)]
pub(crate) fn read_groups_for_tests() -> Vec<ReadGroupIdentity> {
    vec![
        ReadGroupIdentity {
            id: "SRR7279481".to_string(),
            library: "tomato-pe-1".to_string(),
            walk_local_id: ReadGroupId(0),
        },
        ReadGroupIdentity {
            id: "SRR7279481.L2".to_string(),
            library: "tomato-pe-2".to_string(),
            walk_local_id: ReadGroupId(1),
        },
    ]
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
        let contigs = vec![
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
        ];
        let mut header = Header {
            format_version: FORMAT_VERSION,
            sample: "SRR7279481".to_string(),
            reference: ReferenceIdentity {
                name: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                md5: Some([0x0a; 16]),
            },
            segmentation_inputs: segmentation_section::segmentation_inputs_for_tests(&contigs),
            contigs,
            read_groups: read_groups_for_tests(),
            observation_reach_ceiling_bp: Bp(4_000),
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
        };
        // Through the real A4 path, so the recorded filters in every test's header are
        // what `provenance_parameters` writes, not a hand-kept copy of it.
        header.writer.record_parameters(
            crate::ng::read::ReadFilterConfig::default().provenance_parameters(),
        );
        header
    }

    /// **One field per encoding the format has — all eight**, so anything that walks the
    /// manifest meets every one rather than the handful a minimal fixture would carry.
    ///
    /// ⚠ It carried six for a while, and the two it omitted were exactly the two list-shaped
    /// encodings — so `"list"` never once appeared in a header text any test read, and the
    /// whole of that half of the wire vocabulary was pinned only as a side effect of an
    /// error-message assertion. A fixture that omits a case cannot fail on it.
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
                FieldSpec {
                    name: FieldName("chain-id-changes".to_string()),
                    encoding: FieldEncoding::ChainIdChanges,
                },
                FieldSpec {
                    name: FieldName("observation-reads".to_string()),
                    encoding: FieldEncoding::ChainIdList,
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

    /// Two entries sharing one `@RG ID` are **accepted** — deliberately, and pinned so
    /// a later "obvious" uniqueness tightening cannot land quietly. SAM makes the id
    /// unique only within one file; a sample sequenced across two files may declare
    /// `ID:1` in each, and refusing that would refuse real archives the walk handles
    /// fine. The table calling cannot *merge* is the calling stage's refusal (spec
    /// §6.2), not this file's.
    #[test]
    fn two_read_groups_sharing_an_rg_id_round_trip() {
        let mut written = a_written_header();
        written.read_groups[1].id = written.read_groups[0].id.clone();
        assert_ne!(
            written.read_groups[0].library, written.read_groups[1].library,
            "the fixture's two libraries stay distinct, so the rows stay distinguishable"
        );
        let bytes = written.encode().expect("a shared @RG ID encodes");
        let (read_back, _) = decoded(&bytes).expect("and decodes");
        assert_eq!(read_back, written);
    }

    /// A space in an `@RG ID` or a library name is legal — SAM allows it, and real
    /// archives carry it — so the control-character rule must not widen into the
    /// whitespace rule the contig names use. A widened rule would refuse real files.
    #[test]
    fn a_read_group_id_with_a_space_round_trips() {
        let mut written = a_written_header();
        written.read_groups[0].id = "HiSeq 2000 lane 3".to_string();
        written.read_groups[0].library = "prep A".to_string();
        let bytes = written.encode().expect("a spaced @RG ID encodes");
        let (read_back, _) = decoded(&bytes).expect("and decodes");
        assert_eq!(read_back, written);
    }

    /// The default run's shape — analysed regions covering whole contigs — encodes and
    /// decodes through the full header path. The shared fixture deliberately uses proper
    /// sub-spans, so without this a `>` → `>=` regression in the span-end rule would
    /// refuse every whole-genome run while the suite stayed green.
    #[test]
    fn a_whole_genome_header_round_trips() {
        let mut written = a_written_header();
        let bounds = segmentation_section::contig_bounds_of(&written.contigs);
        written.segmentation_inputs.analysed_regions =
            crate::ng::region_typing::GenomeRegions::whole_contigs(&bounds);
        let bytes = written.encode().expect("a whole-genome ground encodes");
        let (read_back, _) = decoded(&bytes).expect("and decodes");
        assert_eq!(read_back, written);
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
        at_the_edge.observation_reach_ceiling_bp = Bp(MAX_TOML_INTEGER);
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
            // The observation reach ceiling, value and all (the fixture's 4,000).
            "observation-reach-ceiling-bp = 4000",
            // The read-group table: the row marker and each row's three keys, with the
            // first row's values pinned whole.
            "[[read-group]]",
            "id = \"SRR7279481\"",
            "library = \"tomato-pe-1\"",
            "walk-local-id = 0",
            "walk-local-id = 1",
            // The read filters, as the fixture's provenance records them.
            "read-filter-min-mapq = ",
            "read-filter-min-read-length-bp = ",
            "read-filter-drop-qc-fail = ",
            "read-filter-drop-duplicates = ",
            "read-filter-max-read-mismatch-fraction = ",
            "read-filter-mismatch-bq-floor = ",
            // The segmentation section: an analysed span, a routing criterion spelled
            // as the short decimal a person would write, and the catalog's identity.
            "[[segmentation.analysed-region]]",
            "min-purity = 0.93",
            "[segmentation.catalog.scan]",
            "[segmentation.catalog.built-under]",
            "[[segmentation.catalog.contig]]",
            "tool-version = \"trf-port-0.9.1\"",
            // Every key the section writes, pinned by name: a round trip cannot catch a
            // serde field rename, and a renamed key orphans every file already written.
            "contig = ",
            "start = ",
            "end = ",
            "period-min = ",
            "period-max = ",
            "min-copies-by-period = ",
            "min-copies-for-wider-periods = ",
            "min-score = ",
            "bundle-threshold-bp = ",
            "min-flank-bp = ",
            "max-str-len-bp = ",
            "reference-md5 = ",
            "match-reward = ",
            "mismatch-penalty = ",
            "min-copies = ",
            "length = ",
            "offset = ",
            "line-bases = ",
            "line-width = ",
            "longest-tract-bp = ",
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
    /// (The key this test first used — the observation reach ceiling — has since been
    /// added for real, so the unknown key is now an invented one.)
    #[test]
    fn a_later_minor_that_added_a_key_still_reads() {
        let mut written = a_written_header();
        let bytes = written.encode().expect("a valid header encodes");
        let body = body_of(&bytes)
            .replace("format-version = \"1.0\"", "format-version = \"1.4\"")
            .replace("sample = ", "a-key-a-later-minor-added = 4000\nsample = ");

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
        // The segmentation section grows with the assembly too — an analysed span and a
        // catalog row per scaffold — so rebuilding it here makes this the whole
        // header's worst case, not just the contig list's.
        fragmented.segmentation_inputs =
            segmentation_section::segmentation_inputs_for_tests(&fragmented.contigs);

        let bytes = fragmented.encode().expect("a fragmented assembly encodes");
        assert!(
            bytes.len() > 1024 * 1024,
            "the fixture must exceed the old 1 MiB cap to be testing anything; it is {} bytes",
            bytes.len()
        );
        eprintln!(
            "30,000 scaffolds with digests, segmentation included: {} bytes of a {} ceiling",
            bytes.len(),
            MAX_HEADER_BODY_BYTES
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
            (
                "an observation reach ceiling of zero",
                Box::new(|header| header.observation_reach_ceiling_bp = Bp(0)),
                "no record could exist",
            ),
            (
                "an observation reach ceiling wider than a TOML integer",
                Box::new(|header| header.observation_reach_ceiling_bp = Bp(MAX_TOML_INTEGER + 1)),
                "a TOML integer is signed",
            ),
            (
                "a read-group table with no entries",
                Box::new(|header| header.read_groups.clear()),
                "renumbered at calling",
            ),
            (
                "a read-group identifier out of walk order",
                Box::new(|header| header.read_groups[1].walk_local_id = ReadGroupId(7)),
                "the walk's own numbering from zero",
            ),
            (
                // The boundary case the row above cannot see: entry 1 repeating
                // identifier 0 is what two zero-numbered walks pasted together look
                // like, and a check loosened from `!=` to `>` would accept it.
                "a read-group identifier repeating zero",
                Box::new(|header| header.read_groups[1].walk_local_id = ReadGroupId(0)),
                "the walk's own numbering from zero",
            ),
            (
                "an empty @RG ID",
                Box::new(|header| header.read_groups[0].id = String::new()),
                "read-group.id",
            ),
            (
                "an @RG ID holding a newline",
                Box::new(|header| header.read_groups[0].id = "rg\nforged = 1".to_string()),
                "holds a control character",
            ),
            (
                // A second control character besides the newline, so the rule cannot
                // quietly narrow to newline-only.
                "an @RG ID holding a tab",
                Box::new(|header| header.read_groups[0].id = "rg\tlane".to_string()),
                "holds a control character",
            ),
            (
                "a library name holding a newline",
                Box::new(|header| header.read_groups[0].library = "lib\nforged = 1".to_string()),
                "holds a control character",
            ),
            (
                "a segmentation recording no analysed ground at all",
                Box::new(|header| {
                    header.segmentation_inputs.analysed_regions =
                        crate::ng::region_typing::GenomeRegions::whole_contigs(&[])
                }),
                "records the ground",
            ),
            (
                "a catalog contig name holding a newline",
                Box::new(|header| {
                    header.segmentation_inputs.catalog.contigs[0].name =
                        "evil\nfabricated-key = \"x\"".to_string()
                }),
                "holds whitespace or a control character",
            ),
            (
                "a catalog tool version holding a newline",
                Box::new(|header| {
                    header.segmentation_inputs.catalog.tool_version = "trf\nport".to_string()
                }),
                "holds whitespace or a control character",
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

    /// The segmentation section's rules fire through the whole header path on both
    /// sides — the writer refuses to produce the file and the reader refuses to believe
    /// it. The section's own tests cover each rule; this pins the wiring.
    #[test]
    fn a_broken_segmentation_section_is_refused_by_writer_and_reader_alike() {
        // Writer side: analysed spans built against a longer contig list than the
        // header declares.
        let mut broken = a_written_header();
        let a_third_contig = ContigIdentity {
            name: "SL4.0ch02".to_string(),
            length: 1_000_000,
            md5: None,
        };
        let mut three = broken.contigs.clone();
        three.push(a_third_contig);
        broken.segmentation_inputs = segmentation_section::segmentation_inputs_for_tests(&three);
        let refused = refusal(
            broken.encode(),
            "the writer must refuse a span outside the file's contig list",
        );
        assert!(
            refused.to_string().contains("segmentation.analysed-region"),
            "got {refused}"
        );

        // Reader side: the same file-level fact, met in a file — a span naming a contig
        // the header does not declare.
        let body = toml::to_string_pretty(&WireHeader::from(&a_written_header()))
            .expect("encodes")
            .replace("contig = \"SL4.0ch00\"", "contig = \"SL4.0ch09\"");
        let refused = refusal(
            decoded(&framed(&body)),
            "the reader must refuse a span naming an undeclared contig",
        );
        assert!(refused.to_string().contains("SL4.0ch09"), "got {refused}");
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
                "name = \"x\"\nshape = \"scalar\"\nencoding = \"varint\"\nsteps-per-unit = 4096\n",
                "a varint with a step",
            ),
            (
                "name = \"x\"\nshape = \"scalar\"\nencoding = \"varint\"\nwidth-bytes = 4\n",
                "a varint with a width",
            ),
            (
                "name = \"x\"\nshape = \"scalar\"\nencoding = \"fixed-width-integer\"\nwidth-bytes = 4\nsteps-per-unit = 4096\n",
                "a fixed-width integer with a step",
            ),
            (
                "name = \"x\"\nshape = \"scalar\"\nencoding = \"fixed-point\"\nsteps-per-unit = 4096\nwidth-bytes = 4\n",
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
                "name = \"body-bytes\"\nshape = \"scalar\"\nencoding = \"fixed-width-integer\"\n",
                "carries no width",
            ),
            (
                "name = \"raw\"\nshape = \"scalar\"\nencoding = \"ieee-float\"\n",
                "carries no width",
            ),
            (
                "name = \"q-sum\"\nshape = \"scalar\"\nencoding = \"fixed-point\"\n",
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
            "name = \"chain-ids\"\nshape = \"scalar\"\nencoding = \"roaring-bitmap\"\n",
        ))
        .expect_err("that is not one of the eight");
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

    // -----------------------------------------------------------------
    // What one appearance of a field looks like
    // -----------------------------------------------------------------

    /// **Each encoding's shape, named one at a time against a literal.**
    ///
    /// ⚠ This is the only test that says what the answers *are*. Every other test of the shape
    /// column asks [`FieldEncoding::shape`] what a scheme lays down and then checks the file
    /// agrees, so a wrong answer applied consistently passes all of them — measured: moving
    /// `SignedVarint`, `FixedWidthInteger` and `IeeeFloat` into the list arm left the suite
    /// green while the header wrote `shape = "list"` beside a 4-byte integer and an 8-byte
    /// float. Those three are reached by no record field today, and they are exactly what the
    /// two queued fields — a window's GC fraction and its mean coverage — will use.
    #[test]
    fn every_encoding_lays_down_the_shape_its_bytes_have() {
        for scheme in [
            FieldEncoding::Varint,
            FieldEncoding::SignedVarint,
            FieldEncoding::LengthPrefixedBytes,
            FieldEncoding::FixedWidthInteger { width_bytes: 4 },
            FieldEncoding::IeeeFloat { width_bytes: 8 },
            FieldEncoding::FixedPoint { steps_per_unit: 4 },
        ] {
            assert_eq!(
                scheme.shape(),
                FieldShape::Scalar,
                "{scheme:?} lays down one value"
            );
        }
        for scheme in [FieldEncoding::ChainIdChanges, FieldEncoding::ChainIdList] {
            assert_eq!(
                scheme.shape(),
                FieldShape::List,
                "{scheme:?} lays down a counted run"
            );
        }
    }

    /// **Every scheme, both ways round**: a file declaring the shape its encoding lays down
    /// round-trips, and one declaring the other shape is refused.
    ///
    /// The disagreeing half is built as **text**, because a `FieldSpec` can no longer hold a
    /// disagreement — which is the point of deriving the shape rather than storing it.
    #[test]
    fn every_scheme_round_trips_with_its_own_shape_and_is_refused_with_the_other() {
        for scheme in ALL_ENCODINGS {
            let laid_down = scheme.shape();
            let mut agreeing = a_written_header();
            agreeing.manifest.fields = vec![FieldSpec {
                name: FieldName("the-one-field".to_string()),
                encoding: scheme,
            }];
            let bytes = agreeing
                .encode()
                .unwrap_or_else(|e| panic!("{scheme:?} must encode: {e}"));
            let (read_back, _) =
                decoded(&bytes).unwrap_or_else(|e| panic!("{scheme:?} must decode: {e}"));
            assert_eq!(read_back.manifest.fields[0].shape(), laid_down);

            let other = match laid_down {
                FieldShape::Scalar => FieldShape::List,
                FieldShape::List => FieldShape::Scalar,
            };
            let (spelling, width_bytes, steps_per_unit) = scheme.spelled();
            let mut declaration = format!(
                "name = \"the-one-field\"\nshape = \"{}\"\nencoding = \"{spelling}\"\n",
                other.spelled()
            );
            if let Some(width) = width_bytes {
                declaration.push_str(&format!("width-bytes = {width}\n"));
            }
            if let Some(step) = steps_per_unit {
                declaration.push_str(&format!("steps-per-unit = {step}\n"));
            }
            let refused = refusal(
                decoded(&one_field_declared_as(&declaration)),
                &format!("{scheme:?} declared {:?} must be refused", other.spelled()),
            );
            assert_eq!(
                refused.to_string(),
                format!(
                    "SRR7279481.psp: manifest.field.shape: field \"the-one-field\" declares \
                     shape {:?}, but its encoding {:?} lays down {:?}",
                    other.spelled(),
                    spelling,
                    laid_down.spelled()
                ),
                "the message must say which declaration was which"
            );
        }
    }

    /// A file whose shape disagrees with its encoding is refused, **and the message names both
    /// declarations the right way round**.
    ///
    /// Asserted as the whole sentence rather than as `contains` of its ingredients: swapping
    /// the two operands leaves every ingredient present, and produces a message claiming
    /// `chain-id-list` lays down a scalar.
    #[test]
    fn a_file_whose_shape_disagrees_with_its_encoding_is_refused_by_the_reader() {
        let refused = refusal(
            decoded(&one_field_declared_as(
                "name = \"observation-reads\"\nshape = \"scalar\"\n\
                 encoding = \"chain-id-list\"\n",
            )),
            "a counted run declared as one value must be refused",
        );
        assert_eq!(
            refused.to_string(),
            "SRR7279481.psp: manifest.field.shape: field \"observation-reads\" declares shape \
             \"scalar\", but its encoding \"chain-id-list\" lays down \"list\"",
        );
    }

    /// A shape this reader has never heard of is refused rather than fallen back from, and the
    /// message lists the two it knows.
    ///
    /// **The near-miss spellings are the point.** A reader that trimmed, folded case or
    /// matched a prefix would accept `"Scalar"`, `" list "` or `"s"` and reach the very
    /// fallback the exact match exists to forbid — and a fixture whose only bad value looks
    /// like neither known spelling cannot tell an exact match from a sloppy one. The field is
    /// deliberately **not** named after the bad value: a fixture named `per-observation-thing`
    /// makes `contains("per-observation")` true even when the message has dropped the value.
    #[test]
    fn an_unknown_shape_is_refused_and_the_message_lists_the_known_ones() {
        for bad in [
            "per-observation",
            "Scalar",
            "LIST",
            " list ",
            "list ",
            "s",
            "scalars",
            "",
        ] {
            let refused = refusal(
                decoded(&one_field_declared_as(&format!(
                    "name = \"a-field\"\nshape = \"{bad}\"\nencoding = \"varint\"\n"
                ))),
                &format!("{bad:?} is not one of the two"),
            );
            let said = refused.to_string();
            assert_eq!(
                said,
                format!(
                    "SRR7279481.psp: manifest.field.shape: field \"a-field\" declares shape \
                     {bad:?}, which is not one of scalar, list"
                ),
                "the message must name the value it refused"
            );
            for known in ALL_SHAPES {
                assert!(
                    said.contains(known.spelled()),
                    "the message must list {:?}; got {said}",
                    known.spelled()
                );
            }
        }
    }

    /// A field entry that declares no shape at all is refused. **Not defaulted**: the key is
    /// one of the three spec §4.5 requires, and no psp exists that omits it.
    ///
    /// The refusal must come from the **parser**, naming the missing key and the line it was
    /// missing from. Adding `#[serde(default)]` to the wire field would still refuse the file —
    /// as an empty string that is not a known spelling — but the message would then claim the
    /// file declared `""` and would lose the line number, so asserting merely that the word
    /// `shape` appears does not hold this in place.
    #[test]
    fn a_field_that_declares_no_shape_is_refused_by_the_parser() {
        let refused = refusal(
            decoded(&one_field_declared_as(
                "name = \"position-offset\"\nencoding = \"varint\"\n",
            )),
            "a field entry with no shape must be refused",
        );
        let said = refused.to_string();
        assert!(said.contains("missing field `shape`"), "got {said}");
        assert!(
            said.contains("[[manifest.field]]"),
            "the message must point at the entry it faulted; got {said}"
        );
    }

    /// The header a reader sees says it in words, **for both spellings**.
    ///
    /// Read off the file's own text, because a round-trip through this module's types would
    /// pass just as well if the key never reached the bytes. Two traps this test has already
    /// fallen into and now covers:
    ///
    /// - **Anchored at the start of a line.** As a bare substring it passed with the key
    ///   renamed to `field-cardinality`, because that name ends in the one being looked for.
    /// - **The fixture must carry a list.** `a_manifest()` declared six of the eight encodings
    ///   and omitted exactly the two list-shaped ones, so `"list"` never appeared in any
    ///   header text any test read.
    #[test]
    fn the_header_text_names_each_fields_shape() {
        let written = a_written_header();
        let bytes = written.encode().expect("a valid header encodes");
        let text = body_of(&bytes);
        assert!(text.contains("\nshape = \"scalar\""), "got {text}");
        assert!(
            text.contains("\nshape = \"list\""),
            "the fixture must declare a list-shaped field, or this test cannot see half the \
             vocabulary; got {text}"
        );
        assert_eq!(
            text.lines()
                .filter(|line| *line == "shape = \"scalar\"")
                .count()
                + text
                    .lines()
                    .filter(|line| *line == "shape = \"list\"")
                    .count(),
            written.manifest.fields.len(),
            "one line per declared field"
        );
    }

    /// **A field name may not carry whitespace or a control character**, because TOML writes
    /// such a name as a multi-line string and its own bytes then land in the header body as
    /// further lines.
    ///
    /// Measured before the rule existed: a field named `evil"\nshape = "list"\nx = "` produced
    /// a header whose text showed a `shape` line for a field that declared another, nine such
    /// lines for eight fields — while the file still round-tripped, because the TOML parser is
    /// not fooled. It is the person running `head` who is, and readability is the whole reason
    /// this header is text.
    #[test]
    fn a_field_name_that_could_forge_a_line_in_the_header_is_refused() {
        for forged in [
            "evil\"\nshape = \"list\"\nx = \"",
            "two words",
            "trailing\t",
            "null\u{0}byte",
        ] {
            let mut header = a_written_header();
            header.manifest.fields = vec![FieldSpec {
                name: FieldName(forged.to_string()),
                encoding: FieldEncoding::Varint,
            }];
            let refused = refusal(
                header.encode(),
                &format!("{forged:?} must not be writable as a field name"),
            );
            let said = refused.to_string();
            assert!(
                said.contains("whitespace or a control character"),
                "got {said}"
            );
        }
    }

    /// **The record's own fields, and which two of the twenty-seven are lists.** A count and
    /// the names, so a field changing shape has to be changed here on purpose — the two halves
    /// of the chain-id column are the only counted runs a record carries.
    #[test]
    fn exactly_two_of_the_records_fields_are_lists() {
        let fields = crate::ng::psp::record::record_fields();
        assert_eq!(fields.len(), 27);
        let lists: Vec<&str> = fields
            .iter()
            .filter(|field| field.shape() == FieldShape::List)
            .map(|field| field.name.0.as_str())
            .collect();
        assert_eq!(lists, ["chain-id-changes", "observation-reads"]);
    }

    /// Every shape the type has is in [`ALL_SHAPES`], spelled exactly once.
    ///
    /// **The list is generated beside the enum**, so a variant cannot be added to one without
    /// reaching the other — which is what this test can *not* establish on its own, and why the
    /// generation exists. What it does hold is that no two shapes share a spelling: two
    /// variants spelled alike would make `shape_of` return the first for both, and the
    /// disagreement check would then pass on a file it should refuse.
    #[test]
    fn every_shape_is_spelled_once_and_differently() {
        for shape in ALL_SHAPES {
            assert_eq!(
                ALL_SHAPES
                    .iter()
                    .filter(|other| other.spelled() == shape.spelled())
                    .count(),
                1,
                "{:?} shares its spelling with another shape",
                shape
            );
            assert!(!shape.spelled().is_empty());
        }
        // Round-tripping every spelling through the parser is what ties the list to the reader:
        // a shape in the list that `shape_of` cannot find would be one the writer can spell and
        // the reader cannot read.
        for shape in ALL_SHAPES {
            let field = WireFieldSpec {
                name: "a-field".to_string(),
                shape: shape.spelled().to_string(),
                encoding: "varint".to_string(),
                width_bytes: None,
                steps_per_unit: None,
            };
            match shape_of(&field) {
                Ok(parsed) => assert_eq!(parsed, *shape),
                Err(broken) => panic!(
                    "{:?} is in the list and the reader will not parse it: {}",
                    shape.spelled(),
                    broken.reason
                ),
            }
        }
    }

    /// `ALL_ENCODINGS` has to be all of them. The exhaustive match is what makes adding an
    /// eighth scheme a compile error here rather than a file one side cannot read.
    #[test]
    fn the_encoding_list_holds_every_scheme_the_type_has() {
        let sample = FieldEncoding::Varint;
        let named = match sample {
            FieldEncoding::Varint
            | FieldEncoding::SignedVarint
            | FieldEncoding::FixedWidthInteger { .. }
            | FieldEncoding::IeeeFloat { .. }
            | FieldEncoding::FixedPoint { .. }
            | FieldEncoding::LengthPrefixedBytes
            | FieldEncoding::ChainIdChanges
            | FieldEncoding::ChainIdList => 8,
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
