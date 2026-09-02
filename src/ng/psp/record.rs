//! One record on the wire: the head a reader uses to decide whether it wants the record, and
//! the body it decodes once it has.
//!
//! # The head
//!
//! The fixed fields at the front of every record, which let a reader decide without building
//! anything. [`RecordEncoder::encode_record`] writes one and hands it back; [`read_record_head`]
//! reads one and bounds the body behind it, which is the whole of the skip — a reader that does
//! not want the record advances past it and touches no byte of its body. [`decode_record`] does
//! both halves and checks that the head's declared length and the body's real shape agree.
//!
//! ```text
//! record = position-offset | reference-span | non-reference-reads
//!          | record-body-byte-count | body
//!          └──────────── the head, as the manifest names it ─────┘   └─ skipped ─┘
//! ```
//!
//! **The chain-id live-set changes belong in that head and are not written yet.** They join it
//! at Milestone E3, between the count and the body; nothing above shows them because no code
//! writes one.
//!
//! **[`RecordHead`] is the fixed part of the head, and only the fixed part.** The chain-id
//! live-set changes — which reads arrived at this position and which left — are in the head
//! too, because they carry state a skipping reader must keep up to date or the merge composes
//! an allele for a read that was never there (spec psp_record_encoding.md §6). They are a
//! variable-length list, 6.42 bytes a position at 293 reads a position, so they are handed
//! straight to the reader's live set rather than stored in a `Copy` struct. Milestone E3 is
//! where they land; nothing in this type has to change for them.
//!
//! A reader takes the head, decides, and either builds the body or advances past it by the byte
//! count the head declared; nothing else in the block has to be touched to make that decision.
//! **Measured
//! on a tomato accession at three reads a position, 7.69 M records: a walk keeping one
//! record in a hundred takes 0.141 s against 0.29 s for one that builds every record —
//! 2.06× faster** (spec §4.3).
//!
//! **The head is not free, and most of its price is not the length field.** It costs 9.2 %
//! of the file at three reads a position and 5.8 % at 279, of which the length field alone
//! is 1.4 % and 3.3 %. The rest is what skippability forces on the body: a record's
//! coverage and its chain ids would otherwise be coded as differences from the previous
//! record, and a reader that skips a body never sees those differences — so both restart at
//! every record instead (spec §4.3).
//!
//! # The body
//!
//! [`encode_record_body`] and [`decode_record_body`] are the other half of this file: one
//! [`SampleLocusObservations`] to bytes and back, exactly, with no compression, no file and
//! nothing read from outside the bytes themselves. The fields and their encodings are
//! [`BODY_FIELDS`], which is what a writer declares in the header's manifest and what a reader
//! refuses a file for disagreeing with — **a fingerprint of the layout rather than a recipe a
//! third party could parse a body from**, since nothing in it says which fields repeat.
//!
//! **Two things are deliberately not in it yet.** A record's chain ids are dropped, because
//! they hold state across records and that is Milestone E of the plan. And the record's
//! coordinate is not in the body at all: it rides in the head, so [`decode_record_body`] is
//! handed the region rather than reading one.

use crate::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation, SsrDetail,
    WitnessedLocusPositions,
};
use crate::ng::psp::chain_ids::{
    LiveSet, LiveSetReader, LiveSetWriter, as_a_read_set, decode_read_list, encode_read_list,
    residual_reads,
};
use crate::ng::psp::header::{FieldEncoding, FieldName, FieldSpec, Manifest};
use crate::ng::types::{ContigId, GenomeRegion, Motif, Position, ReadGroupId, SummedLogError};
use crate::pileup_record::ChainId;
use crate::psp::errors::VarintError;
use crate::psp::varint::{
    decode_i64_svarint, decode_u64_leb128, encode_i64_svarint, encode_u64_leb128,
};

/// What a reader learns about a record before deciding to build it.
///
/// Every field is read from the record's head; none requires touching the body.
/// **`body_bytes` is what makes skipping possible at all** — the encoded bytes carry no
/// separators, so without it a reader that wants to reach the next record must decode
/// every variable-length integer in this one to find where it ends.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RecordHead {
    /// Where the record sits. **Absolute, and rebuilt by the reader** from the block's
    /// first position and the difference the head encodes — the difference restarts at
    /// every block boundary (spec §3.2).
    ///
    /// It is a region and not a position because a record widened by a deletion covers
    /// more than one base, so a reader indexed by position cannot work out what a record
    /// reaches from its start alone. The cohort merge names this as one of two things it
    /// asks of the format.
    pub region: GenomeRegion,
    /// Reads at this locus that supported something other than the reference, summed over the
    /// observations **that spanned the whole locus**.
    ///
    /// **An observation from a read whose witness stopped inside the locus is counted by
    /// neither side.** Its bases end where the read's evidence ended, so there is nothing to
    /// compare them against — which means a locus whose only varying evidence came from partial
    /// observations reports zero here.
    ///
    /// **Zero and "nothing varies here" are the same condition** for the reads that could be
    /// compared, and that is what the cohort's first pass filters on. A count of reads rather
    /// than of alternative alleles: it answers *does anything vary here* just as well, and it
    /// also lets a reader apply a threshold.
    pub non_reference_reads: u32,
    /// Length of the body that follows, in bytes.
    pub body_bytes: BodyByteCount,
}

// ---------------------------------------------------------------------
// The head: the fixed fields a reader judges a record by
// ---------------------------------------------------------------------

/// Every field a record's **head** carries, in encoding order — the fixed part in front of the
/// body, which a reader decodes for every record whether or not it wants the record.
///
/// **Variable-length integers, and that is an implementation choice with a measurement owed.**
/// Spec `psp_file_format.md` §4.3 leaves the width to the manifest and records that a fixed width
/// is quicker to read and costs less than it looks after compression — the four scalar head fields
/// together compressed to 0.077 bytes a record when measured on their own — but that the two
/// encodings have never been compared *in place*. That comparison needs a compressor, which
/// arrives at Milestone D2, so the choice here is the one that composes with every other field and
/// **that comparison belongs to Milestone D2**.
///
/// **The fifth field is not a scalar and is here for the one reason that outranks tidiness.**
/// The chain ids' live-set changes carry state: a reader knows which reads are live only because
/// it applied every arrival and departure since the block began. **A reader that skipped a
/// record's body would never see changes kept there, so its set would go stale and every later
/// record it did want would be wrong** — silently, because a stale set is still a plausible set.
/// So the changes go in front of the body-byte count's reach, where every reader decodes them
/// whether or not it wants the record, and the *exception lists* — the ids of every observation
/// except the residual one — stay in the skippable body, because they carry no state
/// (spec `psp_record_encoding.md` §6, arch §5).
///
/// ⚠ **Switching them then is a format change, not just a manifest one.** This reader is driven
/// by [`RECORD_HEAD_FIELDS`], and [`RecordLayout::from_manifest`] refuses a file declaring
/// anything else — so a build that switched would reject every psp written before it, and a
/// build before it would reject every psp written after. It costs nothing today because no psp
/// exists; from Milestone F it costs a version.
const RECORD_HEAD_FIELDS: [(&str, FieldEncoding); 5] = [
    ("position-offset", FieldEncoding::Varint),
    ("reference-span", FieldEncoding::Varint),
    ("non-reference-reads", FieldEncoding::Varint),
    ("record-body-byte-count", FieldEncoding::Varint),
    ("chain-id-changes", FieldEncoding::ChainIdChanges),
];

const fn record_head_field(position: usize) -> &'static str {
    RECORD_HEAD_FIELDS[position].0
}

const POSITION_OFFSET: &str = record_head_field(0);
const REFERENCE_SPAN: &str = record_head_field(1);
const NON_REFERENCE_READS: &str = record_head_field(2);
const RECORD_BODY_BYTE_COUNT: &str = record_head_field(3);

// ---------------------------------------------------------------------
// The body: one record's fields to bytes and back
// ---------------------------------------------------------------------

/// Every field a record's **body** carries, in encoding order.
///
/// **A field list, not a count of values per record.** Twelve of these twenty-two are written
/// more than once: ten once per *observation* — one of those ten for every observation but the
/// residual — and two more once per *witness run* inside an observation. The last three are
/// written only when the locus kind is a repeat tract. The manifest says
/// what one appearance of a field looks like — a single value, or a counted run of them
/// ([`FieldShape`]); it does not say how often a field appears. So what this list gives a
/// reader is the order, the encodings and each field's own shape; how many of each a body
/// holds is read from the counts the body itself carries.
///
/// [`FieldShape`]: crate::ng::psp::header::FieldShape
///
/// **What it is for: a file's fingerprint of its own layout, and a reader's check against it.**
/// A writer declares this array after [`RECORD_HEAD_FIELDS`] ([`record_fields`]) and a reader refuses any
/// file declaring something else ([`RecordLayout::from_manifest`]), so a file from another version
/// is rejected rather than decoded into plausible nonsense.
///
/// **It does not drive the codec, and the compiler cannot make it.** [`encode_record_body`] and
/// [`decode_record_body`] write and read these fields in this order by hand — which they must,
/// since nothing here says which of them repeat. What holds the two in step is
/// `the_fixture_encodes_to_these_exact_bytes`: reorder the codec and that test fails; reorder
/// this array and the manifest tests fail. **Change both together and you have changed the
/// format — raise the version.**
///
/// **The head's five fields are not here.** They are the part in front of the body (see
/// the module doc) and they arrive with C2; this array is what a reader meets *after* deciding
/// it wants the record.
///
/// **The chain ids are not here either, and that is Milestone E.** A record's chain ids are
/// dropped by [`encode_record_body`] and come back empty from [`decode_record_body`], which is
/// stated on both and pinned by a test rather than left for a reader to discover.
const BODY_FIELDS: [(&str, FieldEncoding); 22] = [
    ("reference-bases", FieldEncoding::LengthPrefixedBytes),
    ("observation-count", FieldEncoding::Varint),
    // Which observation's reads are derived rather than stored, or the observation count for
    // "none of them" — see `residual_observation_of`.
    ("residual-observation", FieldEncoding::Varint),
    // How many reads that observation names. **Written only when there is one**, and it is what
    // turns the reader's guard from an inequality with slack into an equality — see
    // `check_a_derived_read_list`.
    ("residual-read-count", FieldEncoding::Varint),
    // The next ten, once per observation:
    ("observation-bases", FieldEncoding::LengthPrefixedBytes),
    ("witness-run-count", FieldEncoding::Varint),
    // The next two, once per witness run:
    ("witness-run-start", FieldEncoding::Varint),
    ("witness-run-length", FieldEncoding::Varint),
    ("read-group", FieldEncoding::Varint),
    ("reads-showing-the-sequence", FieldEncoding::Varint),
    ("reads-on-the-forward-strand", FieldEncoding::Varint),
    ("summed-log-error", SUMMED_LOG_ERROR_ENCODING),
    ("mapq-sum", FieldEncoding::Varint),
    ("mapq-sum-of-squares", FieldEncoding::Varint),
    ("reads-starting-left-of-the-locus", FieldEncoding::Varint),
    // Every observation except the residual one — the half of the chain-id column that carries
    // no state, so it stays in the skippable body (spec psp_record_encoding.md §6).
    ("observation-reads", FieldEncoding::ChainIdList),
    // Once per record again:
    ("reads-without-observation", FieldEncoding::Varint),
    ("reads-discarded-by-the-depth-cap", FieldEncoding::Varint),
    ("locus-kind", FieldEncoding::Varint),
    // The last three only when the locus kind is a repeat tract:
    ("repeat-motif", FieldEncoding::LengthPrefixedBytes),
    ("repeat-left-flank", FieldEncoding::LengthPrefixedBytes),
    ("repeat-right-flank", FieldEncoding::LengthPrefixedBytes),
];

/// The declared name of one field, by its position in [`BODY_FIELDS`].
///
/// **The one place a field's name is written**, so that the name a header carries and the name
/// a decode error reports cannot drift apart. The decoder's primitives take one of these rather
/// than a hand-written phrase: before this existed there were two vocabularies for the same
/// nineteen fields, and someone reading `the record's repeat motif is unreadable` had to guess
/// that the header's `repeat-motif` was the same field.
const fn body_field(position: usize) -> &'static str {
    BODY_FIELDS[position].0
}

const REFERENCE_BASES: &str = body_field(0);
const OBSERVATION_COUNT: &str = body_field(1);
const RESIDUAL_OBSERVATION: &str = body_field(2);
const RESIDUAL_READ_COUNT: &str = body_field(3);
const OBSERVATION_BASES: &str = body_field(4);
const WITNESS_RUN_COUNT: &str = body_field(5);
const WITNESS_RUN_START: &str = body_field(6);
const WITNESS_RUN_LENGTH: &str = body_field(7);
const READ_GROUP: &str = body_field(8);
const READS_SHOWING_THE_SEQUENCE: &str = body_field(9);
const READS_ON_THE_FORWARD_STRAND: &str = body_field(10);
const SUMMED_LOG_ERROR: &str = body_field(11);
const MAPQ_SUM: &str = body_field(12);
const MAPQ_SUM_OF_SQUARES: &str = body_field(13);
const READS_STARTING_LEFT: &str = body_field(14);
const OBSERVATION_READS: &str = body_field(15);
const READS_WITHOUT_OBSERVATION: &str = body_field(16);
const READS_DISCARDED_BY_THE_DEPTH_CAP: &str = body_field(17);
const LOCUS_KIND: &str = body_field(18);
const REPEAT_MOTIF: &str = body_field(19);
const REPEAT_LEFT_FLANK: &str = body_field(20);
const REPEAT_RIGHT_FLANK: &str = body_field(21);

/// The witness as a whole, for the two faults that are about the set rather than about one of
/// its numbers. Not a declared field: the three `witness-run-*` entries are.
const WITNESS: &str = "witness";

/// How the summed log-error is declared, and it is the type's step rather than a choice.
///
/// [`SummedLogError`] rounds to whole steps of 1/4,096 of a natural log **where the value is
/// computed**, so both routes into the caller — observations read straight from memory and
/// observations read back from a psp — land on the same integer. The file records the step so
/// a reader can turn it back into a quantity; it does not get to pick one (spec
/// psp_record_encoding.md §5.1.1). A file declaring any other step is refused rather than
/// rescaled.
const SUMMED_LOG_ERROR_ENCODING: FieldEncoding = FieldEncoding::FixedPoint {
    steps_per_unit: SummedLogError::STEPS_PER_NAT as u32,
};

const _: () = assert!(
    SummedLogError::STEPS_PER_NAT > 0 && SummedLogError::STEPS_PER_NAT <= u32::MAX as i64,
    "the summed log-error's step has to fit the u32 the manifest declares it in"
);

/// How many bytes a record body may be, as a type.
///
/// **The head declares the body's length in this, and everything that bounds a body reads the
/// width from here** rather than repeating the number. Three places used to spell it
/// independently — the head's field, the decoder's ceiling and the encoder's guard — and a
/// review changed all three one at a time with the whole suite green
/// (`ng_psp_c2_2026-08-26.md`, the refactor-safety pass).
pub type BodyByteCount = u32;

/// The locus-kind tags, as the bytes spell them.
///
/// **A tag is on disk in every file already written, so adding a kind adds a tag and never
/// renumbers one.** `the_locus_kind_tags_are_the_numbers_the_files_carry` is what makes a
/// renumbering a test failure rather than a silent reinterpretation of every tract record in
/// the field.
const KIND_GENERIC: u64 = 0;
/// See [`KIND_GENERIC`].
const KIND_SSR: u64 = 1;
/// See [`KIND_GENERIC`].
const KIND_SSR_BUNDLE: u64 = 2;

/// The fewest bytes one observation can occupy — nine single-byte variable-length integers,
/// one per field, with an empty sequence and a complete witness. Used to bound what a
/// declared observation count may make this reader reserve, never to decide a length.
/// `the_least_an_observation_and_a_run_can_cost_is_what_the_bounds_say` measures it.
const LEAST_BYTES_PER_OBSERVATION: usize = 9;

/// The fewest bytes one witnessed run can occupy: its start and its length, one byte each.
const LEAST_BYTES_PER_WITNESS_RUN: usize = 2;

const _: () = assert!(
    LEAST_BYTES_PER_OBSERVATION > 0 && LEAST_BYTES_PER_WITNESS_RUN > 0,
    "a per-entry byte floor is a divisor, so it cannot be zero"
);

/// Never reserve room for more entries than this, whatever a body's size.
///
/// **An absolute ceiling beside the relative one.** Reserving in proportion to the bytes
/// actually present already stops a declared count of 2⁶⁴−1 from reaching the allocator, but a
/// body is itself bounded only by the head's `u32`, so "proportional to the input" is still
/// gigabytes from one corrupt record — and a caller holds one open file per sample, at up to
/// several thousand. An honest record with more entries than this pays a handful of
/// reallocations, which is nothing beside the case this closes.
const MOST_ENTRIES_RESERVED: usize = 4_096;

/// The longest byte string or list a body may declare.
///
/// **A record body's length is a `u32` in the head** (spec `psp_file_format.md` §4.3), so
/// nothing inside one can honestly declare more bytes than that. A larger value is a byte
/// sequence no writer produced — damage — rather than a buffer that stopped early, and telling
/// the two apart is what Milestone D's restartable read is built on: it grows its buffer for
/// the second and must not for the first.
const MOST_BYTES_A_BODY_CAN_DECLARE: u64 = BodyByteCount::MAX as u64;

/// What a writer declares in the header's manifest: every field of a record, **head then
/// body**, in encoding order, with the encoding this version writes it in.
///
/// **The manifest is how a reader is driven by the file rather than by an assumption**
/// (spec §4.5). It is a fingerprint rather than a parsing recipe — see [`BODY_FIELDS`] for what
/// it does and does not promise.
pub fn record_fields() -> Vec<FieldSpec> {
    fields_this_build_knows().collect()
}

/// The head's five fields and then the body's twenty-two — **this build's**, not a file's, which
/// is the distinction that matters where `manifest.fields` is in scope beside it.
fn fields_this_build_knows() -> impl Iterator<Item = FieldSpec> {
    RECORD_HEAD_FIELDS
        .iter()
        .chain(BODY_FIELDS.iter())
        .map(|(name, encoding)| FieldSpec {
            name: FieldName((*name).to_string()),
            encoding: *encoding,
        })
}

/// How many fields this build knows, before anything a later writer added: the head's five and
/// the body's twenty-two. The fields past it were declared too — by that later writer — which is
/// why the name says whose set it is.
const KNOWN_FIELD_COUNT: usize = RECORD_HEAD_FIELDS.len() + BODY_FIELDS.len();

/// How this reader must read the records in one particular file.
///
/// Built once per file from its manifest and then used for every record, because checking a
/// twenty-seven-field declaration per record on a path that decodes about twenty million records a
/// second (spec §4.5) is work with one possible answer.
///
/// **What it carries is the part that differs between files: the fields this reader does not
/// know.** A later version of the writer may add a per-record scalar — the window's GC fraction
/// and its mean coverage are the two waiting to be computed — and it adds them *after* every
/// field named in [`BODY_FIELDS`]. Each encoding in the closed set measures its own length, so
/// this reader walks past such a field without knowing anything about it: a variable-length
/// integer ends at its own last byte, a fixed-width one is its declared width, and a byte
/// string carries its length in front.
///
/// **⚠ That works for a field written once per record, and this reader cannot tell such a field
/// from one written once per observation.** The manifest says what one *appearance* of a field
/// looks like ([`FieldShape`]) and it does *not* say how many appearances a record holds — a
/// record's own counts do. So a file whose extra field repeats per observation is *accepted*
/// and decoded into plausible nonsense from the second observation onwards; it is not refused,
/// and nothing here can refuse it. A later writer adding a per-observation field must raise the
/// format version, which a reader does refuse (`header.rs`, `UnsupportedVersion`).
///
/// **Nothing here checks a file's declared shape**, and it needs none: `Header::decode` has
/// already refused any file whose shape disagrees with its encoding, and this walk steps over
/// an unknown field by its encoding — the same source the shape is derived from.
///
/// [`FieldShape`]: crate::ng::psp::header::FieldShape
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordLayout {
    /// The encodings of the fields this reader does not recognise, in the order they appear at
    /// the end of every body. Empty for a file this version wrote.
    ///
    /// The declared *names* are not kept: this module's errors name a field with a
    /// `&'static str` and a later writer's names are not static. A fault in one of these
    /// carries the byte offset it happened at, and they are walked in declared order, which is
    /// what identifies it.
    unknown_final_encodings: Vec<FieldEncoding>,
}

impl RecordLayout {
    /// The layout this build of the code writes: every field it knows, and nothing after them.
    ///
    /// # When this is right
    ///
    /// For bytes this process encoded itself. It skips every check
    /// [`from_manifest`](Self::from_manifest) makes — a renamed field, a reordered pair, a
    /// different quantisation step — because there is no file to check against.
    pub fn as_this_build_writes_it() -> Self {
        Self {
            unknown_final_encodings: Vec::new(),
        }
    }

    /// The layout a file declares, checked against what this reader knows.
    ///
    /// **The known fields must come first, in this order, with these encodings.** A file that
    /// renames one, drops one, reorders two or declares one differently is refused — those
    /// are the shapes that would otherwise decode into plausible values rather than failing.
    /// Whatever the manifest lists *after* them is carried as something to walk past
    /// (see the type's own documentation, including what that cannot do).
    pub fn from_manifest(manifest: &Manifest) -> Result<Self, RecordLayoutError> {
        for (position, expected) in fields_this_build_knows().enumerate() {
            let Some(declared) = manifest.fields.get(position) else {
                return Err(RecordLayoutError::MissingField {
                    position,
                    expected: expected.name,
                });
            };
            if declared.name != expected.name {
                return Err(RecordLayoutError::UnexpectedField {
                    position,
                    expected: expected.name,
                    found: declared.name.clone(),
                });
            }
            if declared.encoding != expected.encoding {
                return Err(RecordLayoutError::WrongEncoding {
                    field: declared.name.clone(),
                    expected: expected.encoding,
                    found: declared.encoding,
                });
            }
        }
        Ok(Self {
            unknown_final_encodings: manifest
                .fields
                .iter()
                .skip(KNOWN_FIELD_COUNT)
                .map(|field| field.encoding)
                .collect(),
        })
    }

    /// How many fields this file carries that this reader does not know. Zero for a file
    /// this version wrote.
    pub fn unknown_field_count(&self) -> usize {
        self.unknown_final_encodings.len()
    }
}

/// Why a file's manifest cannot drive this reader.
///
/// **Every variant is an input problem, not a bug** — the same rule the rest of the module
/// keeps, and the reason none of them is a panic. Raised once per file, at open, never per
/// record.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RecordLayoutError {
    #[error("the manifest ends before field {position}, which must be {:?}", expected.0)]
    MissingField {
        position: usize,
        expected: FieldName,
    },
    #[error(
        "the manifest's field {position} is {:?}; this reader reads {:?} there",
        found.0, expected.0
    )]
    UnexpectedField {
        position: usize,
        expected: FieldName,
        found: FieldName,
    },
    #[error(
        "the manifest declares {:?} as {found:?}; this reader reads {expected:?}",
        field.0
    )]
    WrongEncoding {
        field: FieldName,
        expected: FieldEncoding,
        found: FieldEncoding,
    },
}

/// Why a record's bytes could not be read back.
///
/// **Every variant is an input problem, not a bug**, and there are three of them because a
/// caller should do three different things. Either the body **stops** before the record does,
/// and a reader meeting a record straddling its buffer reads more bytes and tries again rather
/// than treating it as damage (Milestone D); or the bytes **say something that cannot be
/// true**, and the file is corrupt; or the bytes are well formed and **name something this
/// reader does not know**, and the instruction is to upgrade the reader.
///
/// ⚠ **The first two are what Milestone D's restartable read branches on**, so a fault put in
/// the wrong class there makes a streaming reader either reject a good record or retry forever
/// on a bad one. `a_body_cut_short_is_truncated_at_every_cut_and_never_malformed` and
/// `a_length_no_body_could_hold_is_malformed_and_never_truncated` are the two tests that hold
/// the line.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RecordDecodeError {
    /// The bytes ran out while this field was being read, and what it declared was possible.
    #[error("the record's {field} runs past the bytes it was given, {bytes_in} bytes in")]
    Truncated {
        field: &'static str,
        /// How far in the reader had got. **A record at 300 reads a position holds hundreds of
        /// observations**, so a field name alone does not say where to look.
        ///
        /// **Measured from the record's first byte** by [`read_record_head`] and
        /// [`decode_record`], which is where a caller holding a block's bytes needs it. It is
        /// measured from the *body's* first byte only by [`decode_record_body`] called on its
        /// own, which has no head in front of it.
        bytes_in: usize,
    },
    /// The bytes were there and cannot mean what they say.
    #[error("the record's {field}, {bytes_in} bytes in, is unreadable: {reason}")]
    Malformed {
        field: &'static str,
        bytes_in: usize,
        reason: String,
    },
    /// The bytes are well formed and name something a later writer added.
    ///
    /// **Not damage.** Only the *major* format version is enforced when a file is opened, and
    /// the unknown-field machinery above exists so that a file from a later *minor* version is
    /// readable — so a locus kind this reader has never heard of is *upgrade the reader*, the
    /// same instruction `header.rs` gives for a newer format.
    #[error("the record's {field}, {bytes_in} bytes in, is {tag}, which this reader does not know")]
    Unsupported {
        field: &'static str,
        bytes_in: usize,
        tag: u64,
    },
}

impl RecordDecodeError {
    /// Move a fault's offset forward by the bytes in front of the parser that reported it.
    ///
    /// **The class is untouched**, unlike [`inside_a_bounded_body`](Self::inside_a_bounded_body):
    /// a sub-parser handed the *rest* of the record rather than a bounded slice still means
    /// "fetch more bytes" when it runs out, because more bytes is exactly what would help.
    pub(super) fn further_in(self, bytes_before: usize) -> Self {
        match self {
            Self::Truncated { field, bytes_in } => Self::Truncated {
                field,
                bytes_in: bytes_before + bytes_in,
            },
            Self::Malformed {
                field,
                bytes_in,
                reason,
            } => Self::Malformed {
                field,
                bytes_in: bytes_before + bytes_in,
                reason,
            },
            Self::Unsupported {
                field,
                bytes_in,
                tag,
            } => Self::Unsupported {
                field,
                bytes_in: bytes_before + bytes_in,
                tag,
            },
        }
    }

    /// Re-express a fault the body reported as a fault in the record that contains it.
    ///
    /// **Two things change, and both are wrong without it.** The offset becomes record-relative,
    /// because the body decoder counted from the body's first byte and a caller holding a block
    /// counts from the record's. And a `Truncated` becomes damage: the body was handed exactly
    /// the byte count the head declared, so a field running off *that* end is not a buffer that
    /// stopped early — no quantity of further bytes changes the answer, and leaving it in the
    /// class Milestone D fetches more bytes for is a retry that never ends.
    fn inside_a_bounded_body(self, head_bytes: usize, body_bytes: usize) -> Self {
        match self {
            Self::Truncated { field, bytes_in } => Self::Malformed {
                field,
                bytes_in: head_bytes + bytes_in,
                reason: format!(
                    "it runs past the {body_bytes} bytes the head declared for the body"
                ),
            },
            Self::Malformed {
                field,
                bytes_in,
                reason,
            } => Self::Malformed {
                field,
                bytes_in: head_bytes + bytes_in,
                reason,
            },
            Self::Unsupported {
                field,
                bytes_in,
                tag,
            } => Self::Unsupported {
                field,
                bytes_in: head_bytes + bytes_in,
                tag,
            },
        }
    }
}

/// One record's body, read back, and how much of the buffer it took.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct DecodedRecordBody {
    pub record: SampleLocusObservations,
    /// How many bytes of the buffer the body occupied. **What the head's `body_bytes` must
    /// equal**: the two disagreeing means the file and this reader disagree about the record's
    /// shape, which is what a version mismatch or a corrupt block looks like. A named field
    /// rather than a tuple element, because `let (record, _) = …` is shorter than the check and
    /// this is the check.
    pub bytes_read: usize,
}

/// **The guard on the one arithmetic in this format that fails silently.**
///
/// A derived list is the live set minus what the other observations named, and getting it wrong
/// by one gives the reference allele a read that does not exist — which the cohort merge composes
/// an allele for without complaint (spec `psp_chain_id_encoding.md` §5). **There is no second copy
/// of the list to compare against**, so the record carries its *length*: one varint against the
/// several bytes the list itself would cost, and it makes the check an equality.
///
/// **⚠ An inequality alone is not enough, and this was measured.** Spec §5 proposes bounding the
/// derived count by the observation's read count — at most `num_obs`, at least half of it, since
/// an identifier names one read or two. That window is `num_obs / 2` wide, and its slack is
/// exactly the number of read pairs whose two mates both cover this record — which is the shape
/// paired-end data has. Measured on a record whose residual names two reads with `num_obs = 4`: a
/// live set carrying **two identifiers no observation named** derives a list of four, passes
/// `4 <= 4` and `8 >= 4`, and the reference allele silently gains two reads. The declared length
/// closes that.
///
/// The inequality is kept beside it, because it is a second, independent statement: it says the
/// *declared* count could describe those reads at all. ⚠ It is the same inequality the walk's own
/// differential against production asserts **on one side only** — `parity.rs` asserts the lower
/// bound; the upper is this reader's own.
fn check_a_derived_read_list(
    body: &FieldReader<'_>,
    derived: &[ChainId],
    declared: u64,
    num_obs: u32,
) -> Result<(), RecordDecodeError> {
    if derived.len() as u64 != declared {
        return Err(body.malformed(
            OBSERVATION_READS,
            format!(
                "a derived list of {} reads where the record says {declared}; the live set and \
                 the lists stored beside it do not agree",
                derived.len()
            ),
        ));
    }
    if !check_a_read_list_against_its_read_count(derived.len(), num_obs) {
        let names = derived.len();
        let reads = u64::from(num_obs);
        return Err(body.malformed(
            OBSERVATION_READS,
            format!(
                "a derived list of {names} reads for an observation showing {reads}; an \
                 identifier names one read or two, so it has to be between {} and {reads}",
                reads.div_ceil(2)
            ),
        ));
    }
    Ok(())
}

/// **Which observation's reads this record derives rather than stores**, or `observations.len()`
/// for none of them.
///
/// **The largest list, because that is the largest saving.** Spec `psp_chain_id_encoding.md` §4
/// says the residual is *named rather than inferred* — "removes any guessing about which 'the
/// reference one' is when observations split by witness and by read group" — and naming it is
/// what lets the choice be the cheapest one rather than a guess at which allele is the
/// reference's. Ties go to the lowest index, so the choice is a function of the record alone and
/// the file stays the same at any worker count (spec `run_streaming.md` §12.1).
///
/// **⚠ And a record whose derivation would not be exact stores every list instead.** A chain id
/// names a read *pair*, and if both mates cover this record showing different sequences the same
/// id is in two observations — the derivation, which is the live set minus the others, would then
/// drop it from the residual. The spec's definition assumes that away; Checkpoint E asks for
/// chain ids that "round-trip **exactly**". So the writer checks, and falls back to naming
/// everything when the two disagree. The cost is one record's residual list; the alternative is a
/// read silently missing from the reference allele, which is exactly the failure §5 names.
fn residual_observation_of(
    record: &SampleLocusObservations,
    scratch: &mut ResidualScratch,
) -> usize {
    let none = record.observations.len();
    // **The one-observation record, which is about ninety-nine positions in a hundred at the
    // tomato panel's depth.** Its only observation names every read the record names, so
    // `named_elsewhere` below is empty, `every_read` is that observation's own list sorted and
    // deduplicated, and the derivation reproduces the list exactly precisely when the list was
    // already a strictly ascending set. That is one pass over the identifiers, where the general
    // path builds two lists, sorts both, derives, and compares.
    if none == 1 {
        let only = &record.observations[0];
        if only.chain_ids.is_empty() {
            return none;
        }
        let already_a_set = only.chain_ids.windows(2).all(|pair| pair[0] < pair[1]);
        if already_a_set
            && check_a_read_list_against_its_read_count(only.chain_ids.len(), only.num_obs)
        {
            return 0;
        }
        return none;
    }
    let Some(at) = (0..none).max_by_key(|at| {
        (
            record.observations[*at].chain_ids.len(),
            std::cmp::Reverse(*at),
        )
    }) else {
        return none;
    };
    if record.observations[at].chain_ids.is_empty() {
        return none;
    }

    // **Three lists a record, and none of them allocates after the first record.** They were
    // three `collect()`s and a `Vec::new` when this was written — three allocations a record at
    // any depth, plus the growth chain on each.
    let ResidualScratch {
        named_elsewhere,
        every_read,
        derived,
    } = scratch;

    named_elsewhere.clear();
    for (other, observation) in record.observations.iter().enumerate() {
        if other != at {
            named_elsewhere.extend_from_slice(&observation.chain_ids);
        }
    }
    as_a_read_set(named_elsewhere);

    every_read.clear();
    for observation in &record.observations {
        every_read.extend_from_slice(&observation.chain_ids);
    }
    as_a_read_set(every_read);

    residual_reads(every_read, named_elsewhere, derived);
    // **Two conditions, and the second is what keeps the reader's guard a guard.** The derivation
    // has to reproduce the list exactly, and it has to satisfy the inequality the reader checks it
    // against — otherwise this writer would produce a file its own reader refuses. Real evidence
    // always satisfies it (an identifier names one read or two), so the second condition only
    // fires on a record whose counts do not describe reads; storing that record's lists costs its
    // residual's bytes and nothing else.
    let exact = derived.as_slice() == record.observations[at].chain_ids.as_slice();
    let within_the_inequality =
        check_a_read_list_against_its_read_count(derived.len(), record.observations[at].num_obs);
    if exact && within_the_inequality {
        at
    } else {
        none
    }
}

/// Whether a list of `names` identifiers can describe an observation showing `num_obs` reads.
///
/// **At most `num_obs`**, because an identifier names at least one read and the list is
/// deduplicated; **at least half of it**, because at most two mates share one.
fn check_a_read_list_against_its_read_count(names: usize, num_obs: u32) -> bool {
    let names = names as u64;
    let reads = u64::from(num_obs);
    names <= reads && names * 2 >= reads
}

/// Append `record`'s body to `out`, and never fail.
///
/// **The chain ids are dropped.** [`SequenceObservation::chain_ids`] is Milestone E of
/// `doc/devel/ng/impl_plan/psp_file_format.md`, so a record encoded here and read back has
/// empty lists where it had ids, and everything else identical. Encoding a record that
/// carries them is not an error today and produces a file the cohort merge could not use;
/// the writer that produces real files arrives after E does.
///
/// **Appends rather than returning a buffer**, so a writer keeps one and clears it between
/// records instead of allocating per record; the same path decodes at about twenty million
/// records a second in the measuring prototype (spec §4.5).
///
/// # Nothing here is a difference from another record, and the signature is the guarantee
///
/// **A body stands on its own.** Every count it carries is written absolutely — the reads behind
/// each observation, the reads that showed nothing, the reads the depth cap dropped — and each
/// observation's own list of reads is written from zero rather than from the previous record's.
/// **That is why a skipped body costs nothing**: a reader that never saw it has missed nothing
/// the next record needs. The head is the half that carries state (spec §4.3), which is why the
/// chain ids' live-set *changes* go there and their exception lists stay here.
///
/// **⚠ With one exception, and it is the whole of Milestone E4**: the decoder takes a `&LiveSet`,
/// because one observation's reads are not in the body at all — they are the live set minus every
/// other observation's. That is a difference from the records before this one, and it is the one
/// thing here a skipped body *would* have stranded had it been stored rather than derived. It is
/// safe because the set is built from the record *heads*, which a skipping reader decodes anyway.
///
/// **⚠ The live set the decoder is given must be the one this record's head produced** — after
/// its changes are applied, not before. Anything else is accepted quietly whenever the derived
/// count still matches what the record declares.
///
/// **Getting it wrong would be silent**: a body carrying a difference still round-trips under a
/// walk that builds every record, and misreads only the records that follow a skipped one. The
/// tests under this module's `C3` banner are what see it; the injected-defect run that checked
/// they can is in `doc/devel/reports/implementations/ng_psp_c3_2026-08-26.md`.
///
/// It cannot fail. Every field of a [`SampleLocusObservations`] has a representation here and
/// the variable-length integers are unbounded; the one subtraction, a witnessed run's length,
/// cannot underflow because [`WitnessedLocusPositions`] keeps its runs private and its
/// constructors reject a run that covers no position.
pub fn encode_record_body(record: &SampleLocusObservations, out: &mut Vec<u8>) {
    encode_record_body_reusing(record, out, &mut ResidualScratch::default());
}

/// The lists the residual derivation needs, kept by whoever encodes more than one record.
///
/// **Every one of them is cleared and refilled a record**, so what a writer holds is the widest
/// record it has met — three lists of one identifier a read, about 7.2 kB at three hundred reads
/// a position. **A writer is one per sample being written, not one per sample held open for a
/// run**, so this is not spent against spec §1.1's 500 kB an open sample. A reader has no
/// equivalent: it derives its residual straight into the record it hands over.
#[derive(Debug, Default)]
pub struct ResidualScratch {
    /// The identifiers every observation but the residual one names, as a set.
    named_elsewhere: Vec<ChainId>,
    /// Every identifier the record names, as a set — the live set the residual is derived
    /// against.
    every_read: Vec<ChainId>,
    /// What the derivation produced, to be compared against the list it would replace.
    derived: Vec<ChainId>,
}

/// Append `record`'s body to `out`, reusing `scratch` rather than allocating the residual
/// derivation's three lists afresh.
///
/// See [`encode_record_body`], which is this with a scratch of its own and is the whole of the
/// difference between them.
pub fn encode_record_body_reusing(
    record: &SampleLocusObservations,
    out: &mut Vec<u8>,
    scratch: &mut ResidualScratch,
) {
    // Destructured with no `..`: **a field added to either type is a compile error here**, not
    // a field silently left out of every psp written from that day. The compiler already stops
    // at `decode_record_body`'s struct literals, but naming a field there is not encoding it.
    let SampleLocusObservations {
        region: _, // rides in the head, not the body
        reference_bases,
        observations,
        reads_without_observation,
        reads_discarded_by_cap,
        kind,
    } = record;

    put_bytes(out, reference_bases);
    put_varint(out, observations.len() as u64);
    let residual_at = residual_observation_of(record, scratch);
    put_varint(out, residual_at as u64);
    if let Some(residual) = observations.get(residual_at) {
        put_varint(out, residual.chain_ids.len() as u64);
    }
    for (at, observation) in observations.iter().enumerate() {
        let SequenceObservation {
            bases,
            read_witness,
            read_group,
            num_obs,
            num_fwd,
            q_sum,
            mapq_sum,
            mapq_sum_sq,
            placed_left,
            chain_ids,
        } = observation;
        put_bytes(out, bases);
        put_witness(out, read_witness);
        put_varint(out, u64::from(read_group.get()));
        put_varint(out, u64::from(*num_obs));
        put_varint(out, u64::from(*num_fwd));
        put_signed_varint(out, q_sum.steps());
        put_varint(out, u64::from(*mapq_sum));
        put_varint(out, *mapq_sum_sq);
        put_varint(out, u64::from(*placed_left));
        // **Every observation except the residual one.** Its reads are the live set minus these,
        // which is where most of this column's saving is (spec `psp_chain_id_encoding.md` §4).
        //
        // **`encode_read_list` makes the list a set itself**, because the bytes are gaps and
        // `SequenceObservation::chain_ids` says what the ids mean rather than what order they are
        // in. A precondition here instead would be one a `debug_assert` cannot hold in the
        // profile that writes real files.
        if at != residual_at {
            encode_read_list(chain_ids, out);
        }
    }
    put_varint(out, u64::from(*reads_without_observation));
    put_varint(out, u64::from(*reads_discarded_by_cap));
    put_kind(out, kind);
}

/// Read one record's body back, and say how many bytes it took.
///
/// `region` comes from the record's head, which is where the format keeps it — a body carries
/// no coordinate of its own (spec §4.3), and nothing here checks the region against the
/// reference bases the body holds, because a record's reference bases are not required to
/// cover its region. `layout` comes from the file's manifest, once, at open
/// ([`RecordLayout::from_manifest`]); [`RecordLayout::as_this_build_writes_it`] is for
/// bytes this process encoded itself and skips every check that function makes.
///
/// **The chain-id lists come back empty** — see [`encode_record_body`], which also states the
/// property this half has to keep: nothing this reads out of `bytes` depends on any record read
/// before it, which is what lets a reader skip a body and still trust the next one. (`region`
/// does come from every record before it, through the head's running offset — but it is a
/// parameter, not something the body's bytes encode.)
///
/// **Every sequence is copied out.** [`SampleLocusObservations`] owns its bases as `Box<[u8]>`,
/// so a record costs one allocation for the reference bases and one per observation. The
/// reading primitives hand back slices borrowed from `bytes`, and spec §4.6's second reason for
/// rejecting `serde` is exactly that borrow — but nothing consumes it yet, because the record
/// type is owned by design (arch §2). Whether a borrowed record view is worth building is
/// Milestone D3's question, with its own measurement.
///
/// **The decode is not injective.** A witness whose runs touch is normalised on the way in, and
/// LEB128 admits a non-canonical encoding of a small number, so two different bodies can
/// produce one record. Compare records, never bytes, when checking a reader against a writer:
/// `encode(decode(b)) == b` is not a property this format has.
///
/// **⚠ `live_reads` must be the set this record's head produced.** See
/// [`decode_the_body_of`], which is the call a reader that has read a head should be making.
pub fn decode_record_body(
    bytes: &[u8],
    region: GenomeRegion,
    live_reads: &LiveSet,
    layout: &RecordLayout,
) -> Result<DecodedRecordBody, RecordDecodeError> {
    let mut body = FieldReader::new(bytes);

    let reference_bases: Vec<u8> = body.read_length_prefixed(REFERENCE_BASES)?.into();

    let declared_observations = body.read_count(OBSERVATION_COUNT, LEAST_BYTES_PER_OBSERVATION)?;
    let residual_at = body.read_varint(RESIDUAL_OBSERVATION)?;
    if residual_at > declared_observations {
        return Err(body.malformed(
            RESIDUAL_OBSERVATION,
            format!(
                "observation {residual_at} of a record holding {declared_observations}; the count \
                 itself means none of them are derived"
            ),
        ));
    }
    let residual_reads_declared = if residual_at < declared_observations {
        Some(body.read_count(RESIDUAL_READ_COUNT, LEAST_BYTES_PER_OBSERVATION)?)
    } else {
        None
    };
    let mut observations = Vec::with_capacity(entries_to_reserve(
        declared_observations,
        LEAST_BYTES_PER_OBSERVATION,
        body.bytes_left(),
    ));
    let mut named_elsewhere: Vec<ChainId> = Vec::new();
    let mut one_list: Vec<ChainId> = Vec::new();
    for at in 0..declared_observations {
        let bases: Vec<u8> = body.read_length_prefixed(OBSERVATION_BASES)?.into();
        let read_witness = body.read_witness()?;
        let read_group = ReadGroupId(body.read_u32(READ_GROUP)?);
        let num_obs = body.read_u32(READS_SHOWING_THE_SEQUENCE)?;
        let num_fwd = body.read_u32(READS_ON_THE_FORWARD_STRAND)?;
        let q_sum = SummedLogError::from_steps(body.read_signed_varint(SUMMED_LOG_ERROR)?);
        let mapq_sum = body.read_u32(MAPQ_SUM)?;
        let mapq_sum_sq = body.read_varint(MAPQ_SUM_OF_SQUARES)?;
        let placed_left = body.read_u32(READS_STARTING_LEFT)?;
        let chain_ids = if at == residual_at {
            // Filled after the loop, from the live set minus every list above.
            Vec::new()
        } else {
            decode_read_list(&mut body, OBSERVATION_READS, &mut one_list)?;
            named_elsewhere.extend_from_slice(&one_list);
            std::mem::take(&mut one_list)
        };
        observations.push(SequenceObservation {
            bases,
            read_witness,
            read_group,
            num_obs,
            num_fwd,
            q_sum,
            mapq_sum,
            mapq_sum_sq,
            placed_left,
            chain_ids,
        });
    }

    if let Ok(at) = usize::try_from(residual_at)
        && at < observations.len()
    {
        as_a_read_set(&mut named_elsewhere);
        let mut derived = Vec::new();
        residual_reads(live_reads.ids(), &named_elsewhere, &mut derived);
        check_a_derived_read_list(
            &body,
            &derived,
            residual_reads_declared.unwrap_or_default(),
            observations[at].num_obs,
        )?;
        observations[at].chain_ids = derived;
    }

    let reads_without_observation = body.read_u32(READS_WITHOUT_OBSERVATION)?;
    let reads_discarded_by_cap = body.read_u32(READS_DISCARDED_BY_THE_DEPTH_CAP)?;
    let kind = body.read_locus_kind()?;

    for encoding in &layout.unknown_final_encodings {
        body.skip_unknown_field(*encoding)?;
    }

    Ok(DecodedRecordBody {
        record: SampleLocusObservations {
            region,
            reference_bases,
            observations,
            reads_without_observation,
            reads_discarded_by_cap,
            kind,
        },
        bytes_read: body.bytes_read(),
    })
}

// ---------------------------------------------------------------------
// Writing primitives
// ---------------------------------------------------------------------

pub(super) fn put_varint(out: &mut Vec<u8>, value: u64) {
    encode_u64_leb128(value, out);
}

fn put_signed_varint(out: &mut Vec<u8>, value: i64) {
    encode_i64_svarint(value, out);
}

fn put_bytes(out: &mut Vec<u8>, bytes: &[u8]) {
    put_varint(out, bytes.len() as u64);
    out.extend_from_slice(bytes);
}

/// A witness as a count of runs and then the runs themselves — **zero runs meaning the read
/// reached both borders**, which is [`ReadWitness::Complete`] and by far the common case, so
/// it costs one byte.
///
/// **Each run is written as its own start and its own length, not as a step from the run
/// before.** A step would be one byte smaller on a witness with holes in it and it would
/// depend on the runs arriving sorted and non-overlapping — an invariant
/// [`WitnessedLocusPositions`] does keep, but one this encoder would then be silently wrong
/// about if it ever stopped keeping it. Locus positions are at most 65,535, so an absolute
/// start costs at most three bytes and on a real locus costs one.
fn put_witness(out: &mut Vec<u8>, witness: &ReadWitness) {
    match witness {
        ReadWitness::Complete => put_varint(out, 0),
        ReadWitness::Partial { positions } => {
            put_varint(out, positions.runs().len() as u64);
            for (start, end) in positions.runs() {
                put_varint(out, u64::from(start));
                put_varint(out, u64::from(end - start));
            }
        }
    }
}

/// The kind tag, and whatever that kind carries.
///
/// **Exhaustive with no wildcard**, though [`LocusKind`] is `#[non_exhaustive]`: that
/// attribute binds other crates, not this one, so a kind added later is a compile error here
/// rather than a record written without its payload. The decoder cannot be made exhaustive the
/// same way — it matches a tag read from a file — so `every_locus_kind_round_trips` is what
/// stops a kind gaining a write side and no read side.
fn put_kind(out: &mut Vec<u8>, kind: &LocusKind) {
    match kind {
        LocusKind::Generic => put_varint(out, KIND_GENERIC),
        LocusKind::Ssr(detail) => {
            put_varint(out, KIND_SSR);
            put_bytes(out, detail.motif.as_bytes());
            put_bytes(out, &detail.left_flank);
            put_bytes(out, &detail.right_flank);
        }
        LocusKind::SsrBundle => put_varint(out, KIND_SSR_BUNDLE),
    }
}

// ---------------------------------------------------------------------
// Reading primitives
// ---------------------------------------------------------------------

/// A record's bytes being read, and how far into them the reader has got — the head's fields
/// and the body's alike.
///
/// **It never indexes past what it holds**: every method that advances checks the bytes are
/// there first, so a truncated or hostile record produces an error rather than a panic. The
/// invariant every method keeps is `bytes_read <= bytes.len()`, which is what makes
/// [`bytes_left`](Self::bytes_left) and the slicing below total.
///
/// Every method here advances the cursor, which is why each is named `read_` or `skip_` — a
/// column of them in [`decode_record_body`] has to say that it consumes as it goes.
///
/// **`pub(super)` so that [`chain_ids`](super::chain_ids) reads its bytes through this and not
/// through a second spelling.** The line this type draws — a field that ran out of bytes is
/// `Truncated` and a field that cannot mean what it says is `Malformed` — is what Milestone D's
/// restartable reader branches on, and a chain-id stream that drew it slightly differently would
/// make a streaming reader either reject a good record or retry for ever on a bad one.
pub(super) struct FieldReader<'a> {
    bytes: &'a [u8],
    bytes_read: usize,
}

impl<'a> FieldReader<'a> {
    pub(super) fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            bytes_read: 0,
        }
    }

    pub(super) fn bytes_read(&self) -> usize {
        self.bytes_read
    }

    pub(super) fn bytes_left(&self) -> usize {
        // Saturating rather than a bare subtraction: the invariant holds, and stating it here
        // costs nothing where relying on it silently would cost a panic.
        self.bytes.len().saturating_sub(self.bytes_read)
    }

    fn truncated(&self, field: &'static str) -> RecordDecodeError {
        RecordDecodeError::Truncated {
            field,
            bytes_in: self.bytes_read,
        }
    }

    pub(super) fn malformed(&self, field: &'static str, reason: String) -> RecordDecodeError {
        RecordDecodeError::Malformed {
            field,
            bytes_in: self.bytes_read,
            reason,
        }
    }

    /// One variable-length integer, read through production's codec.
    /// Advance past bytes another parser has already accounted for.
    ///
    /// **Bounded, unlike a plain addition**, because this is the one advancing method whose
    /// argument comes from somewhere else. The type's invariant is `bytes_read <= bytes.len()`,
    /// and `read_varint` slices on it unconditionally — so an overshoot here is a panic in a
    /// decoder whose whole contract is that corrupt input gives an error. Today's one caller
    /// cannot overshoot; the next one has no compiler help, which is why the bound is here
    /// rather than in a comment about the caller.
    fn skip(&mut self, bytes: usize) {
        debug_assert!(
            self.bytes_read + bytes <= self.bytes.len(),
            "a skip past the end would make every later slice a panic"
        );
        self.bytes_read = self.bytes_read.saturating_add(bytes).min(self.bytes.len());
    }

    pub(super) fn read_varint(&mut self, field: &'static str) -> Result<u64, RecordDecodeError> {
        self.accept(decode_u64_leb128(&self.bytes[self.bytes_read..]), field)
    }

    /// One zig-zag variable-length integer — what a value that goes negative is written as.
    fn read_signed_varint(&mut self, field: &'static str) -> Result<i64, RecordDecodeError> {
        self.accept(decode_i64_svarint(&self.bytes[self.bytes_read..]), field)
    }

    /// Advance by what the varint codec consumed, and turn its fault into this module's.
    ///
    /// **Both integer readers meet the two classes here**, so the line between "this record
    /// stopped early" — which Milestone D reads more bytes for — and "no writer produced these
    /// bytes" is drawn once.
    fn accept<T>(
        &mut self,
        decoded: Result<(T, usize), VarintError>,
        field: &'static str,
    ) -> Result<T, RecordDecodeError> {
        match decoded {
            Ok((value, used)) => {
                self.bytes_read += used;
                Ok(value)
            }
            Err(VarintError::Truncated) => Err(self.truncated(field)),
            Err(VarintError::Overflow) => Err(self.malformed(
                field,
                "a variable-length integer longer than any 64-bit value needs".to_string(),
            )),
        }
    }

    fn read_u32(&mut self, field: &'static str) -> Result<u32, RecordDecodeError> {
        let value = self.read_varint(field)?;
        u32::try_from(value).map_err(|_| {
            self.malformed(
                field,
                format!("{value}, which is past the {} this field holds", u32::MAX),
            )
        })
    }

    fn read_u16(&mut self, field: &'static str) -> Result<u16, RecordDecodeError> {
        let value = self.read_varint(field)?;
        u16::try_from(value).map_err(|_| {
            self.malformed(
                field,
                format!("{value}, which is past the {} this field holds", u16::MAX),
            )
        })
    }

    /// A count of entries the body says follow, refused when no body could hold that many.
    ///
    /// **Bounding it here is what keeps `Truncated` meaning what it says.** A count of 2⁴⁰
    /// observations is not a buffer that stopped early — no record body is that long — so
    /// reporting it as a short read would ask Milestone D's reader to grow its buffer to a
    /// terabyte instead of reporting damage.
    pub(super) fn read_count(
        &mut self,
        field: &'static str,
        least_bytes_each: usize,
    ) -> Result<u64, RecordDecodeError> {
        let declared = self.read_varint(field)?;
        // **A multiply where this was a divide**, and the two say the same thing: for positive
        // integers `declared > M / L` and `declared * L > M` are the same statement, because the
        // truncation in the division is exactly the slack the multiplication keeps. This runs
        // about five times a record — once for the observation count, once for the residual read
        // count, once a witness, once a stored read list, twice from the live set's changes — and
        // an aarch64 64-bit divide is tens of cycles to compare against a constant. The quotient
        // is still computed for the message, inside the branch that refuses.
        let most_possible = MOST_BYTES_A_BODY_CAN_DECLARE / least_bytes_each as u64;
        if declared > most_possible {
            return Err(self.malformed(
                field,
                format!(
                    "{declared} entries, more than the {most_possible} a record body's \
                     {MOST_BYTES_A_BODY_CAN_DECLARE} bytes could hold"
                ),
            ));
        }
        Ok(declared)
    }

    fn take(&mut self, count: usize, field: &'static str) -> Result<&'a [u8], RecordDecodeError> {
        let to = self
            .bytes_read
            .checked_add(count)
            .filter(|to| *to <= self.bytes.len())
            .ok_or_else(|| self.truncated(field))?;
        let taken = &self.bytes[self.bytes_read..to];
        self.bytes_read = to;
        Ok(taken)
    }

    /// A byte string, refused when it declares more bytes than a record body can hold.
    ///
    /// See [`read_count`](Self::read_count) for why the ceiling matters rather than being
    /// tidiness: past it, the declared length is damage rather than a short buffer.
    fn read_length_prefixed(&mut self, field: &'static str) -> Result<&'a [u8], RecordDecodeError> {
        let declared = self.read_varint(field)?;
        if declared > MOST_BYTES_A_BODY_CAN_DECLARE {
            return Err(self.malformed(
                field,
                format!(
                    "{declared} bytes, more than the {MOST_BYTES_A_BODY_CAN_DECLARE} a record \
                     body can be"
                ),
            ));
        }
        // Bounded above by a `u32`, so the cast cannot lose anything.
        self.take(declared as usize, field)
    }

    /// The positions of the locus a read witnessed — no runs meaning it reached both borders.
    ///
    /// **Runs are required to ascend and never touch, rather than being sorted and merged into
    /// whatever they describe.** [`WitnessedLocusPositions::from_half_open_runs`] would happily
    /// normalise `(0,5),(3,9)` into the single run `(0,9)` — a witness the encoder could not
    /// have written, decoding into a *different valid record* with the byte count still right,
    /// so even C2's length check would not see it. That is the one class this decoder would
    /// otherwise repair instead of reporting, and depth is raised over a witness, so a merged
    /// run is a read credited with positions it never saw.
    fn read_witness(&mut self) -> Result<ReadWitness, RecordDecodeError> {
        let declared_runs = self.read_count(WITNESS_RUN_COUNT, LEAST_BYTES_PER_WITNESS_RUN)?;
        if declared_runs == 0 {
            return Ok(ReadWitness::Complete);
        }
        let mut runs = Vec::with_capacity(entries_to_reserve(
            declared_runs,
            LEAST_BYTES_PER_WITNESS_RUN,
            self.bytes_left(),
        ));
        let mut previous_end: Option<u16> = None;
        for _ in 0..declared_runs {
            let start = self.read_u16(WITNESS_RUN_START)?;
            let length = self.read_u16(WITNESS_RUN_LENGTH)?;
            if length == 0 {
                return Err(self.malformed(
                    WITNESS_RUN_LENGTH,
                    format!("a run at {start} covering no position"),
                ));
            }
            if previous_end.is_some_and(|previous_end| start <= previous_end) {
                let previous_end = previous_end.unwrap_or_default();
                return Err(self.malformed(
                    WITNESS,
                    format!(
                        "a run starting at {start} after one ending at {previous_end}; \
                         a witness's runs ascend and never touch"
                    ),
                ));
            }
            let end = start.checked_add(length).ok_or_else(|| {
                self.malformed(
                    WITNESS,
                    format!(
                        "a run from {start} covering {length} positions ends past the last \
                         locus position a witness can name"
                    ),
                )
            })?;
            previous_end = Some(end);
            runs.push((start, end));
        }
        let positions = WitnessedLocusPositions::from_half_open_runs(runs)
            .ok_or_else(|| self.malformed(WITNESS, "runs that describe no position".to_string()))?;
        Ok(ReadWitness::Partial { positions })
    }

    fn read_locus_kind(&mut self) -> Result<LocusKind, RecordDecodeError> {
        let at = self.bytes_read;
        let tag = self.read_varint(LOCUS_KIND)?;
        match tag {
            KIND_GENERIC => Ok(LocusKind::Generic),
            KIND_SSR => {
                let motif_bases = self.read_length_prefixed(REPEAT_MOTIF)?;
                let motif = Motif::new(motif_bases)
                    .map_err(|source| self.malformed(REPEAT_MOTIF, source.to_string()))?;
                let left_flank: Box<[u8]> = self.read_length_prefixed(REPEAT_LEFT_FLANK)?.into();
                let right_flank: Box<[u8]> = self.read_length_prefixed(REPEAT_RIGHT_FLANK)?.into();
                Ok(LocusKind::Ssr(SsrDetail {
                    motif,
                    left_flank,
                    right_flank,
                }))
            }
            KIND_SSR_BUNDLE => Ok(LocusKind::SsrBundle),
            tag => Err(RecordDecodeError::Unsupported {
                field: LOCUS_KIND,
                bytes_in: at,
                tag,
            }),
        }
    }

    /// Walk past one field of a declared encoding without interpreting it — how a reader
    /// meets a field a later writer added (see [`RecordLayout`]).
    ///
    /// **Exhaustive with no wildcard**, so an encoding added to the closed set has to say
    /// here how long one of its values is, rather than inheriting a guess.
    fn skip_unknown_field(&mut self, encoding: FieldEncoding) -> Result<(), RecordDecodeError> {
        let field = "field a later writer added";
        match encoding {
            FieldEncoding::Varint
            | FieldEncoding::SignedVarint
            | FieldEncoding::FixedPoint { .. } => {
                self.read_varint(field)?;
            }
            FieldEncoding::FixedWidthInteger { width_bytes }
            | FieldEncoding::IeeeFloat { width_bytes } => {
                self.take(usize::from(width_bytes), field)?;
            }
            FieldEncoding::LengthPrefixedBytes => {
                self.read_length_prefixed(field)?;
            }
            // **Measurable, which is not the same as usable.** Two counted runs of varints, so
            // its length is a function of its own bytes and a reader can step over it. But this
            // encoding carries state every later record is decoded against, so a reader that
            // stepped over one it *needed* would build every record after it against a stale
            // live set. That is why the head's copy is a field this reader knows by name
            // (`RecordLayout::from_manifest` refuses a file that renames or moves it) rather
            // than something reached through this arm. What this arm handles is a *later*
            // writer putting another one at the end of a body, which nothing here uses.
            FieldEncoding::ChainIdChanges => {
                let departures = self.read_count(field, 1)?;
                for _ in 0..departures {
                    self.read_varint(field)?;
                }
                let arrivals = self.read_count(field, 1)?;
                for _ in 0..arrivals {
                    self.read_varint(field)?;
                }
            }
            // One counted run, where the changes above are two. ⚠ The two shared a scheme
            // briefly, and stepping over a list under the changes' rule measures the wrong number
            // of bytes: `[5, 9]` is three bytes as a list and seven as a set of changes.
            FieldEncoding::ChainIdList => {
                let reads = self.read_count(field, 1)?;
                for _ in 0..reads {
                    self.read_varint(field)?;
                }
            }
        }
        Ok(())
    }
}

/// How many entries to reserve room for, when a body declares how many follow.
///
/// **Never the declared count on its own.** A corrupt or hostile body can say it holds a
/// million observations in eleven bytes, and a reader that believes it asks the allocator for
/// the number before reading a single one. Each entry costs at least a byte or two, so what is
/// left of the body bounds what can really be there; and [`MOST_ENTRIES_RESERVED`] bounds that
/// in turn, because a body is itself bounded only by a `u32`. The loop still reads the declared
/// count and still fails when the bytes run out, which is where a wrong length is reported.
pub(super) fn entries_to_reserve(
    declared: u64,
    least_bytes_each: usize,
    bytes_left: usize,
) -> usize {
    let could_be_there = (bytes_left / least_bytes_each).min(MOST_ENTRIES_RESERVED);
    declared.min(could_be_there as u64) as usize
}

// ---------------------------------------------------------------------
// The head, and the skip it exists for
// ---------------------------------------------------------------------

/// Why a record could not be laid down.
///
/// **Every variant is a record the writer was handed that the format cannot hold**, not an
/// internal fault. The body encoder cannot fail at all; these are the four things the *head*
/// has to be able to say no to, because each of them would otherwise be written as a number
/// that means something else.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RecordEncodeError {
    /// A record starting before the position its offset is measured from. The head stores a
    /// distance forwards, and a distance backwards has no representation.
    ///
    /// **Refused rather than accepted**: coordinate order is what the block index and every
    /// seek rest on, and a file that breaks it seeks wrongly rather than failing. The writer
    /// of Milestone F3 turns this into `PspWriteError::OutOfOrder`, which names the file.
    ///
    /// The message says *the position its offset is measured from* rather than *the previous
    /// record*, because for a block's first record there is no previous record — the base is
    /// the block's own first position.
    #[error(
        "a record starting at {} starts before {}, the position its offset is measured from",
        offered_start.get(), previous_position.get()
    )]
    StartsBeforeThePreviousRecord {
        previous_position: Position,
        offered_start: Position,
    },

    /// A body longer than the head can describe. The head's byte count is what a reader
    /// advances to skip a record, so a length it cannot hold is a record nothing could skip.
    #[error(
        "a record body of {body_bytes} bytes, longer than the {} a head can describe",
        BodyByteCount::MAX
    )]
    BodyTooLong { body_bytes: usize },

    /// A region covering no bases, which is `end` before `start`. `GenomeRegion` has public
    /// fields and no constructor, so this is reachable and is a caller's mistake rather than a
    /// state the format has a spelling for.
    #[error("a record over {region}, which covers no reference base")]
    EmptyRegion { region: GenomeRegion },

    /// A region whose last base is the last coordinate a `u64` can name.
    ///
    /// **Refused rather than written, and the reader refuses it too**, so that a writer cannot
    /// produce a record its own reader would reject. The reason is not this module's:
    /// [`GenomeRegion::len`] is documented as wrong at exactly that coordinate — it overflows
    /// in a debug build and reports the region empty in a release one — so a head handing such
    /// a region to a consumer hands out a value that detonates when anything asks its width.
    /// It is not reachable from real contig coordinates, which is why refusing it costs
    /// nothing.
    #[error("a record over {region}, whose last base is the last coordinate there is")]
    EndsAtTheCoordinateCeiling { region: GenomeRegion },
}

/// The coordinate a record's position offset is measured from: the block's own first position
/// for a block's first record, and the previous record's start after that
/// (spec `psp_file_format.md` §3.2).
///
/// **A type rather than a bare [`Position`], because every `Position` type-checks in that
/// argument** — including the record's own start, which is the slip that reads most naturally
/// and puts every record in a block at the same coordinate. A wrong base does not fail: the
/// head carries a difference and nothing else, so a reader given one that is stale by one
/// record rebuilds every later record at the wrong coordinate and reports success.
/// `a_base_stale_by_one_record_moves_every_record_after_it` is what prices that.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct OffsetBase(Position);

impl OffsetBase {
    /// A block's first record is measured from the block's own first position — which is that
    /// record's start, so its offset is zero.
    pub fn at_block_start(first_position: Position) -> Self {
        Self(first_position)
    }

    /// Every record after it is measured from the record before it, taken from the head that
    /// was just written or just read so that a caller cannot invent one.
    #[must_use]
    pub fn after(head: &RecordHead) -> Self {
        Self(head.region.start)
    }

    /// The coordinate itself, for a message or a comparison.
    pub fn position(self) -> Position {
        self.0
    }
}

/// Lays records down, and keeps what one block's worth of writing needs.
///
/// **A record's head cannot be written until its body has been**, because the head ends with
/// the body's length in bytes — that field is the whole reason a reader can skip a record
/// rather than decoding every variable-length integer in it to find where it ends. So the body
/// goes into a scratch buffer first and is copied out behind its head. **The scratch is kept
/// across records**, which is what stops a writer allocating one per record; the same path
/// decodes at about twenty million records a second in the measuring prototype (spec §4.5), and
/// nobody has measured the writer's own rate.
///
/// **It also holds the base every position offset is measured from, and advances it itself.**
/// That was a parameter until a review priced the alternative: a caller that forgets to advance
/// it writes a file at coordinates nobody asked for, without an error anywhere. Milestone D1
/// gives the block its own reader-side state; this is the writer's half of it.
#[derive(Debug)]
pub struct RecordEncoder {
    /// Kept across blocks: reused so a writer does not allocate one per record.
    body_scratch: Vec<u8>,
    /// The residual derivation's three lists, kept for the same reason and cleared a record.
    residual_scratch: ResidualScratch,
    /// Which reads are live, and the changes it writes to move between records.
    ///
    /// **Its own per-block state is inside it**, reset by
    /// [`encode_record_starting_a_block`](Self::encode_record_starting_a_block) — which is why
    /// that is the only way to start a block. It sits here rather than in [`PerBlockState`]
    /// because it owns scratch buffers that must survive a boundary; what must not survive one
    /// is the set, and that is [`LiveSetWriter::start_block`]'s business.
    live_reads: LiveSetWriter,
    /// Replaced whole at every block start.
    block: PerBlockState,
}

/// Everything of the encoder's that **restarts when a block does** (spec §3.2).
///
/// **A struct rather than loose fields**, because a field added beside a running difference and
/// initialised once is one that silently never resets: measured at Milestone D, a field added to
/// [`RecordEncoder`] and initialised in its constructor alone builds clean, passes every test,
/// and writes a file that is wrong from each block's first record — plausibly wrong, because
/// coverage is smooth.
///
/// **⚠ The second running difference is *not* in here, and that is deliberate.** Milestone E's
/// chain-id live set lives in [`RecordEncoder::live_reads`] beside this struct, because
/// [`LiveSetWriter`] owns scratch buffers that must survive a block boundary — only the *set*
/// must not, and `LiveSetWriter::start_block` is what empties it. What keeps that from being the
/// hazard this doc describes is that there is one way to open a block,
/// [`RecordEncoder::encode_record_starting_a_block`], and it resets both.
///
/// **What the split buys, exactly, because it is less than it looks.** A field added *here* is
/// a compile error at `at_block_start`, which is the only constructor, so it cannot be left out
/// of the reset — measured, `error[E0063]: missing field … in initializer of PerBlockState`.
/// A field added to [`RecordEncoder`] *beside* this one and initialised in `for_block` still
/// compiles and still is never reset — also measured, 148 tests green. So this does not force
/// the choice; it makes the reset automatic once the choice is made, and names the two
/// lifetimes so that making it wrong is visible rather than invisible.
///
/// **⚠ Not `Copy`, deliberately** — the same reason `BlockCursor` on the reading side is not.
/// Milestone E's chain-id live set is a collection, and a `Copy` struct cannot hold one: adding
/// it here under a `Copy` derive gives `error[E0204]: the trait `Copy` cannot be implemented for
/// this type` *before* the `E0063` that is the point, and the cheapest way past `E0204` is to put
/// the field on [`RecordEncoder`] instead, where it compiles and is never reset. Dropping `Copy`
/// leaves only the error that sends the coder to the right place.
#[derive(Debug, Clone)]
struct PerBlockState {
    measured_from: OffsetBase,
}

impl PerBlockState {
    fn at_block_start(first_position: Position) -> Self {
        Self {
            measured_from: OffsetBase::at_block_start(first_position),
        }
    }
}

impl RecordEncoder {
    /// A writer for one block, whose first record sits at `first_position`.
    ///
    /// **⚠ A second way in, and the weaker one.** It sets the coordinate base from its argument
    /// and leaves the live set empty, so a caller that reaches
    /// [`encode_record`](Self::encode_record) without going through
    /// [`encode_record_starting_a_block`](Self::encode_record_starting_a_block) first gets a
    /// block measured from whatever it passed — a file that parses perfectly and is wrong from
    /// its first record if that was a guess. `BlockBuilder` passes a placeholder and always
    /// opens through the other one; this exists for the codec's own tests, which lay records down
    /// with no builder above them.
    pub fn for_block(first_position: Position) -> Self {
        Self {
            body_scratch: Vec::new(),
            residual_scratch: ResidualScratch::default(),
            live_reads: LiveSetWriter::new(),
            block: PerBlockState::at_block_start(first_position),
        }
    }

    /// **Write the record that opens a block**, resetting everything a block restarts first.
    ///
    /// This is the only way to *re*start a block, and that is the point. Spec §3.2 requires every
    /// running difference to restart at a block boundary, and there are two of them now — the
    /// coordinate the next offset is measured from, and the set of reads live. A separate
    /// `start_block` call that a caller had to remember beside `encode_record` would be one more
    /// thing to forget, and forgetting it writes a file that parses perfectly and is wrong from
    /// this record onward.
    ///
    /// **Nothing is reset until the record is known to be writable.** The refusals below depend
    /// on the record alone, not on either running difference, so they are made first — which is
    /// what lets a caller meet one at a block cut and carry on with the block it already had
    /// open, rather than needing to put back a live set that a reset had already thrown away.
    pub fn encode_record_starting_a_block(
        &mut self,
        record: &SampleLocusObservations,
        out: &mut Vec<u8>,
    ) -> Result<RecordHead, RecordEncodeError> {
        let span = record_span(record.region)?;
        self.body_scratch.clear();
        encode_record_body_reusing(record, &mut self.body_scratch, &mut self.residual_scratch);
        let body_bytes = declared_body_bytes(self.body_scratch.len())?;

        self.block = PerBlockState::at_block_start(record.region.start);
        self.live_reads.start_block();
        self.write_a_record(record, 0, span, body_bytes, out)
    }

    /// What the next record's offset will be measured from — the block's first position until a
    /// record has been written, and the last written record's start after that.
    pub fn measured_from(&self) -> OffsetBase {
        self.block.measured_from
    }

    /// The reads live at the record last written.
    pub fn live_reads(&self) -> &LiveSet {
        self.live_reads.live()
    }

    /// Append `record` to `out` as a head and then a body, and hand back the head that was
    /// written.
    ///
    /// **`out` is untouched unless the record is written.** Every refusal happens before the
    /// first byte, so a caller that meets one and carries on writes a block with no trace of
    /// the rejected record rather than a block whose next record starts inside a fragment.
    ///
    /// The head's `non_reference_reads` is **derived here, not supplied**: it is the reads at
    /// this locus that showed something other than the reference, which the record already
    /// knows and which a caller could otherwise get wrong.
    pub fn encode_record(
        &mut self,
        record: &SampleLocusObservations,
        out: &mut Vec<u8>,
    ) -> Result<RecordHead, RecordEncodeError> {
        let span = record_span(record.region)?;
        let offset = record
            .region
            .start
            .get()
            .checked_sub(self.block.measured_from.position().get())
            .ok_or(RecordEncodeError::StartsBeforeThePreviousRecord {
                previous_position: self.block.measured_from.position(),
                offered_start: record.region.start,
            })?;

        self.body_scratch.clear();
        encode_record_body_reusing(record, &mut self.body_scratch, &mut self.residual_scratch);
        let body_bytes = declared_body_bytes(self.body_scratch.len())?;

        self.write_a_record(record, offset, span, body_bytes, out)
    }

    /// The bytes, once every refusal has been made and the body is in `body_scratch`.
    ///
    /// **The chain-id changes go between the body's byte count and the body itself**, which is
    /// what puts them in the head: `body_bytes` does not reach them, so a reader that skips the
    /// body has already decoded them. That is the whole of spec §6's resolution — the half of
    /// the column that carries state is in front of the skip, and the half that does not stays
    /// behind it.
    fn write_a_record(
        &mut self,
        record: &SampleLocusObservations,
        offset: u64,
        span: u64,
        body_bytes: BodyByteCount,
        out: &mut Vec<u8>,
    ) -> Result<RecordHead, RecordEncodeError> {
        let (non_reference_reads, _) = record.non_reference_and_compared_reads();

        put_varint(out, offset);
        put_varint(out, span);
        put_varint(out, u64::from(non_reference_reads));
        put_varint(out, u64::from(body_bytes));
        self.live_reads.write_changes(
            record
                .observations
                .iter()
                .flat_map(|observation| observation.chain_ids.iter().copied()),
            out,
        );
        out.extend_from_slice(&self.body_scratch);

        let head = RecordHead {
            region: record.region,
            non_reference_reads,
            body_bytes,
        };
        self.block.measured_from = OffsetBase::after(&head);
        Ok(head)
    }
}

/// How many reference bases `region` covers, or the refusal if the format cannot say.
///
/// **Not [`GenomeRegion::len`]**, which its own documentation calls wrong when the region's last
/// base is the last coordinate a `u64` can name: it adds one to `end` before subtracting, so it
/// overflows in a debug build and reports the region empty in a release one. Deriving the width
/// first and adding one afterwards moves the only overflow to the whole-axis region, and both
/// that and the ceiling are refused here — a record the writer cannot describe is refused
/// rather than written as a number meaning something else.
fn record_span(region: GenomeRegion) -> Result<u64, RecordEncodeError> {
    if region.is_empty() {
        return Err(RecordEncodeError::EmptyRegion { region });
    }
    if region.end.get() == u64::MAX {
        return Err(RecordEncodeError::EndsAtTheCoordinateCeiling { region });
    }
    region
        .end
        .get()
        .checked_sub(region.start.get())
        .and_then(|width| width.checked_add(1))
        .ok_or(RecordEncodeError::EmptyRegion { region })
}

/// The head's byte count for a body of `len` bytes, or the refusal if no head can describe it.
///
/// **Lifted out of [`RecordEncoder::encode_record`] so that it has a test**: reaching it through
/// the encoder needs a body over four gibibytes, which no test may allocate.
fn declared_body_bytes(len: usize) -> Result<BodyByteCount, RecordEncodeError> {
    BodyByteCount::try_from(len).map_err(|_| RecordEncodeError::BodyTooLong { body_bytes: len })
}

/// One record found at the front of a buffer: what its head says, where its body sits, and how
/// many bytes the whole record takes.
///
/// **This is the skip.** A reader that does not want the record advances
/// [`record_bytes`](Self::record_bytes) and touches nothing else; a reader that does hands
/// [`body`](Self::body) to [`decode_record_body`]. What that is worth is in the module's own
/// documentation, with the sample it was measured on.
///
/// *Named for having been located rather than for the buffer it was found in: the format gives
/// "buffer" to the reader's own two (spec §4.4), and this is in neither of them particularly.*
#[derive(Debug, Clone, PartialEq, Eq)]
#[must_use]
pub struct LocatedRecord<'a> {
    /// What the head said — enough to decide whether the body is wanted at all.
    pub head: RecordHead,
    /// **Exactly the byte count the head declared, and bounded here.** A body handed a slice of
    /// its own record cannot read into the next one however damaged it is, so a length field
    /// that disagrees with the body's real shape is caught rather than absorbed.
    pub body: &'a [u8],
    /// The head and the body together — what a reader advances to reach the next record.
    ///
    /// Deliberately not `bytes_read`: on this type the body's bytes were *not* read, and
    /// [`DecodedRecordBody::bytes_read`] means the different thing it says.
    pub record_bytes: usize,
}

/// Read one record's head, and bound its body.
///
/// `contig` comes from the block, which never crosses one (spec §3.2), and `measured_from` is
/// the base the head's offset is a distance from — see [`OffsetBase`], and note that a base
/// which is merely *wrong* is not an error the format can detect.
///
/// **Nothing in the body is touched**, which is the point: this is what the cohort's first pass
/// runs, and on the corner it was measured on — the tomato panel at about three reads a
/// position — it is all that runs at roughly ninety-nine positions in a hundred
/// (`cohort_merge.md`, "roughly one position in a hundred is variable at the measured corner").
///
/// **⚠ This parses the chain ids' changes and leaves them waiting; the caller applies them.**
/// The head carries the arrivals and departures, and the live set moves only when the caller
/// calls [`LiveSetReader::apply_the_changes_just_parsed`] — which it must do for **every** record
/// it advances past, wanted or not, or every record after this one is built against a stale set.
/// `live_reads` is the walk's own reader — one per block, emptied by
/// [`LiveSetReader::start_block`] at every boundary — not a fresh one per record.
///
/// **Why the caller and not this function**: between reading a head and committing to the record
/// there is one thing that can still fail or refuse — the walk's predicate, which is a caller's
/// closure and may panic. Applying here would move the set for a record the walk then re-parses
/// from its first byte, and the second parse meets an arrival for a read already live: a sound
/// file refused as damaged. A caller that has read a head reaches the body through
/// [`decode_the_body_of`], never by parsing the same bytes again.
///
/// # Errors
///
/// - [`RecordDecodeError::Truncated`] when the bytes stop inside the record: **fetch more and
///   call this again from the record's first byte.** `live_reads` is exactly where it was, which
///   is what makes that a retry rather than a second application of this record's changes.
/// - [`RecordDecodeError::Malformed`] when the head cannot mean what it says — a span of no
///   bases, a coordinate off the axis, or a chain-id stream inconsistent with the live set.
pub fn read_record_head<'a>(
    bytes: &'a [u8],
    contig: ContigId,
    measured_from: OffsetBase,
    live_reads: &mut LiveSetReader,
) -> Result<LocatedRecord<'a>, RecordDecodeError> {
    let mut reader = FieldReader::new(bytes);

    let offset = reader.read_varint(POSITION_OFFSET)?;
    let start = measured_from
        .position()
        .get()
        .checked_add(offset)
        .ok_or_else(|| {
            reader.malformed(
                POSITION_OFFSET,
                format!(
                    "{offset} past {}, which is off the coordinate axis",
                    measured_from.position().get()
                ),
            )
        })?;

    let span = reader.read_varint(REFERENCE_SPAN)?;
    if span == 0 {
        return Err(reader.malformed(
            REFERENCE_SPAN,
            "a record covering no reference base; every record covers at least one".to_string(),
        ));
    }
    // `end` is the last base covered, so a one-base record has `end == start`. The ceiling is
    // refused for the reason `RecordEncodeError::EndsAtTheCoordinateCeiling` gives: a region
    // ending on the last coordinate reports its own width wrongly, so handing one out puts the
    // fault in whichever consumer asks — the writer refuses the same regions, so no file this
    // module wrote can reach here.
    let end = start
        .checked_add(span - 1)
        .filter(|end| *end < u64::MAX)
        .ok_or_else(|| {
            reader.malformed(
                REFERENCE_SPAN,
                format!("{span} bases from {start}, which runs off the coordinate axis"),
            )
        })?;

    let non_reference_reads = reader.read_u32(NON_REFERENCE_READS)?;
    let body_bytes = reader.read_u32(RECORD_BODY_BYTE_COUNT)?;

    // **The chain ids' changes, and every reader decodes them.** They are read here rather than
    // in the body because they carry the state later records are decoded against: a reader that
    // skipped them would build every record after this one against a stale live set — silently,
    // because a stale set is still a plausible set (spec `psp_record_encoding.md` §6).
    //
    // ⚠ **Parsed here and applied below, with the body's bounding in between.** A body that stops
    // early is a `Truncated`, whose contract is *fetch more bytes and re-parse this record from
    // its first byte* — so the set may not have moved by then. Measured on the version that
    // applied them here: a well-formed file of 1,999 records read a byte at a time was refused at
    // record 149 with *"id 150, which is already live"*, and a record that only departs reads
    // retried to `Ok` with a read silently gone from the set for the rest of the block.
    //
    // ⚠ **And the applying is the caller's, one call per record advanced past.** It is deferred
    // that far because the walk's predicate runs in between and is a caller's closure: a set moved
    // before a closure that panics is a set moved for a record the walk then re-parses, and the
    // second parse refuses a sound file. [`decode_record`] and the block walk both apply the
    // moment they commit to the record, and both reach the body through the head this returns
    // rather than parsing the head again — which is what keeps it exactly once.
    let changes_bytes = live_reads
        .parse_changes(&bytes[reader.bytes_read()..])
        .map_err(|fault| fault.further_in(reader.bytes_read()))?;
    reader.skip(changes_bytes);

    let body = reader.take(body_bytes as usize, RECORD_BODY_BYTE_COUNT)?;

    Ok(LocatedRecord {
        head: RecordHead {
            region: GenomeRegion {
                contig,
                start: Position(start),
                end: Position(end),
            },
            non_reference_reads,
            body_bytes,
        },
        body,
        record_bytes: reader.bytes_read(),
    })
}

/// One whole record, read back.
#[derive(Debug, Clone, PartialEq)]
#[must_use]
pub struct DecodedRecord {
    pub record: SampleLocusObservations,
    /// What the head said before the body was built — including
    /// [`non_reference_reads`](RecordHead::non_reference_reads), which the body does not carry.
    pub head: RecordHead,
    /// The head and the body together: what a reader advances to reach the next record.
    pub record_bytes: usize,
}

/// Read one whole record — head, then body — and check that the two agree.
///
/// **The check is the reason this exists rather than being two calls.** The head declares how
/// long the body is and the body says how much of itself it used; the two disagreeing means the
/// file and this reader disagree about the record's shape, which is what a version mismatch or a
/// corrupt block looks like. Split across two calls it is a comparison a caller can forget to
/// make, and forgetting it is silent.
///
/// **The head's non-reference read count is checked too, against the body just built.** It is
/// derivable, and the disagreement that matters reads *low*: a varying position whose head says
/// zero is one the cohort's first pass skips, and nothing downstream looks at it again. The
/// check costs one more pass over the observations already in hand, which is small beside
/// building them.
///
/// **⚠ And it is two calls now, which this doc used to argue against.** Reading a head moves the
/// live set, so the body must be reached through [`decode_the_body_of`] rather than by parsing
/// the same bytes twice — that is what this function does internally, and what a caller that has
/// already read a head must do instead of calling this.
///
/// **A fault inside the body is damage here, never a short read**, and that is the third thing
/// this function does that its two halves cannot. [`read_record_head`] hands the body a slice of
/// exactly the length the head declared, so a field running off the end of *that* is not a
/// buffer that stopped early — no quantity of further bytes changes the answer. Left in the
/// `Truncated` class it would ask Milestone D's reader to fetch more and try again, for ever, on
/// one damaged record.
pub fn decode_record(
    bytes: &[u8],
    contig: ContigId,
    measured_from: OffsetBase,
    live_reads: &mut LiveSetReader,
    layout: &RecordLayout,
) -> Result<DecodedRecord, RecordDecodeError> {
    let found = read_record_head(bytes, contig, measured_from, live_reads)?;
    // **The head is committed to here**: nothing between this and the body can refuse the record,
    // so this is where the live set moves. `decode_the_body_of` reads the set below.
    live_reads.apply_the_changes_just_parsed();
    decode_the_body_of(&found, live_reads.live(), layout)
}

/// Build the record a head has already located.
///
/// **A caller that has read the head passes it here rather than the bytes**, because reading a
/// head applies that record's chain-id changes to the live set and doing it twice would apply
/// them twice — a read arriving into a set it is already in, refused as damage, on a file that is
/// perfectly good.
///
/// **⚠ `live_reads` must be the set this record's head produced** — after
/// [`read_record_head`] applied its changes, not before, and not a fresh one. One observation's
/// reads are not in the body: they are that set minus every other observation's. A different set
/// is refused only when the count it derives disagrees with the one the record declares, which is
/// most wrong sets and not all of them.
pub fn decode_the_body_of(
    found: &LocatedRecord<'_>,
    live_reads: &LiveSet,
    layout: &RecordLayout,
) -> Result<DecodedRecord, RecordDecodeError> {
    let head_bytes = found.record_bytes - found.body.len();
    let decoded_body = decode_record_body(found.body, found.head.region, live_reads, layout)
        .map_err(|fault| fault.inside_a_bounded_body(head_bytes, found.body.len()))?;
    let (non_reference_reads, _) = decoded_body.record.non_reference_and_compared_reads();
    if non_reference_reads != found.head.non_reference_reads {
        return Err(RecordDecodeError::Malformed {
            field: NON_REFERENCE_READS,
            bytes_in: head_bytes,
            reason: format!(
                "a head reading {} over a body whose reads show {non_reference_reads}",
                found.head.non_reference_reads
            ),
        });
    }
    if decoded_body.bytes_read != found.body.len() {
        return Err(RecordDecodeError::Malformed {
            field: RECORD_BODY_BYTE_COUNT,
            bytes_in: head_bytes + decoded_body.bytes_read,
            reason: format!(
                "a head declaring {} body bytes over a body that used {}",
                found.body.len(),
                decoded_body.bytes_read
            ),
        });
    }
    Ok(DecodedRecord {
        record: decoded_body.record,
        head: found.head,
        record_bytes: found.record_bytes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A reader for one block's worth of records.
    ///
    /// **One of these per walk, not one per record.** A record's head carries which reads started
    /// covering it and which stopped, and a departure is a position in the set as the previous
    /// record left it — so a fresh reader meets the second record's departures against nothing.
    fn a_live_set_reader() -> LiveSetReader {
        let mut reader = LiveSetReader::new();
        reader.start_block();
        reader
    }

    /// [`decode_record`], with the walk's reader carried in.
    fn decode_a_record(
        live_reads: &mut LiveSetReader,
        bytes: &[u8],
        contig: ContigId,
        measured_from: OffsetBase,
        layout: &RecordLayout,
    ) -> Result<DecodedRecord, RecordDecodeError> {
        decode_record(bytes, contig, measured_from, live_reads, layout)
    }

    /// [`read_record_head`] and the apply that goes with it, with the walk's reader carried in.
    ///
    /// **The apply is here because these tests read a head and then ask what the live set holds**,
    /// which is the state a committed record leaves. A test about the *deferral* — that a head
    /// read and not committed to leaves the set where it was — calls the two separately and says
    /// so.
    fn read_a_record_head<'a>(
        live_reads: &mut LiveSetReader,
        bytes: &'a [u8],
        contig: ContigId,
        measured_from: OffsetBase,
    ) -> Result<LocatedRecord<'a>, RecordDecodeError> {
        let found = read_record_head(bytes, contig, measured_from, live_reads)?;
        live_reads.apply_the_changes_just_parsed();
        Ok(found)
    }

    /// A head is read once per record on a path that runs at about twenty million records
    /// a second, and a reader that had to clone one — or that paid a pointer chase to
    /// reach it — would be paying that on every record it skips. `Copy` is the property
    /// that says it does not.
    ///
    /// **This is a bound on the type, not on the head on disk.** The chain-id changes are
    /// part of the record's head and are not part of this struct (see the module doc), so a
    /// growing wire head does not have to move this number.
    #[test]
    fn a_head_is_copied_not_cloned_and_stays_small() {
        let head = RecordHead {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1_000),
                end: Position(1_002),
            },
            non_reference_reads: 3,
            body_bytes: 47,
        };
        let copied = head;
        assert_eq!(copied, head, "a head is Copy, so this is not a move");
        assert!(
            std::mem::size_of::<RecordHead>() <= 32,
            "the head is {} bytes; it is read once per record and holds no allocation",
            std::mem::size_of::<RecordHead>()
        );
    }

    /// The span is a field of the region rather than something a reader derives from the
    /// next record's start: a record widened by a deletion reaches past the record that
    /// follows it, so the distance between two starts is not a span.
    #[test]
    fn a_head_carries_a_span_wider_than_one_base() {
        let deletion = RecordHead {
            region: GenomeRegion {
                contig: ContigId(1),
                start: Position(90_667_287),
                end: Position(90_667_293),
            },
            non_reference_reads: 162,
            body_bytes: 231,
        };
        assert_eq!(deletion.region.len(), 7);
    }

    // -----------------------------------------------------------------
    // C1 — the record body, encoded and decoded
    // -----------------------------------------------------------------

    use crate::ng::psp::header::{
        ContigIdentity, DEFAULT_GENOMIC_BLOCK_SIZE_BP, DEFAULT_LOOK_BACK_WINDOW_LOG,
        FORMAT_VERSION, Header, Manifest, ReferenceIdentity, WriterProvenance,
    };
    use proptest::prelude::*;

    fn a_region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(3),
            start: Position(start),
            end: Position(end),
        }
    }

    /// A witness with a hole in it — the shape the generic fold mints for a read spliced
    /// across a widened record, and the one a two-number witness could not describe.
    fn a_witness_with_a_hole() -> ReadWitness {
        ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 2), (5, 7)])
                .expect("two runs with a gap are a canonical set"),
        }
    }

    /// A record with every kind of thing a body carries: three observations, a complete
    /// witness and two partial ones, two read groups, a negative error sum, a sequence longer
    /// than the reference's and one that is empty.
    ///
    /// **Its seven reference bases cover its seven-base region**, which is what a producer
    /// emits. The codec does not check that — a record's reference bases are not required to
    /// cover its region, and one test below turns on their not doing so — but every later
    /// milestone builds on this fixture, and a fixture that cannot occur teaches the wrong
    /// shape.
    ///
    /// **Every count in it is non-zero**, and `the_fixture_leaves_no_field_at_a_value_a_decoder_
    /// could_invent` is what keeps it that way: a field added to the record and filled in here
    /// with the value `decode_record_body` invents for it would round-trip by coincidence.
    fn a_rich_record() -> SampleLocusObservations {
        SampleLocusObservations {
            region: a_region(90_667_287, 90_667_293),
            reference_bases: b"ACGTACG".to_vec(),
            observations: vec![
                SequenceObservation {
                    bases: b"ACGTACG".to_vec(),
                    read_witness: ReadWitness::Complete,
                    read_group: ReadGroupId(1),
                    num_obs: 137,
                    num_fwd: 61,
                    q_sum: SummedLogError::from_nats(-42.5),
                    mapq_sum: 8_220,
                    mapq_sum_sq: 493_200,
                    placed_left: 44,
                    chain_ids: Vec::new(),
                },
                SequenceObservation {
                    bases: b"ACGTACGTACGTAC".to_vec(),
                    read_witness: a_witness_with_a_hole(),
                    read_group: ReadGroupId(7),
                    num_obs: 3,
                    num_fwd: 2,
                    q_sum: SummedLogError::from_nats(-0.25),
                    mapq_sum: 180,
                    mapq_sum_sq: 10_800,
                    placed_left: 1,
                    chain_ids: Vec::new(),
                },
                SequenceObservation {
                    bases: b"".to_vec(),
                    read_witness: ReadWitness::Partial {
                        positions: WitnessedLocusPositions::one_run_from_offset_and_length(4, 2)
                            .expect("one run of two positions"),
                    },
                    read_group: ReadGroupId(7),
                    num_obs: 1,
                    num_fwd: 1,
                    q_sum: SummedLogError::from_steps(-7),
                    mapq_sum: 60,
                    mapq_sum_sq: 3_600,
                    placed_left: 2,
                    chain_ids: Vec::new(),
                },
            ],
            reads_without_observation: 12,
            reads_discarded_by_cap: 305,
            kind: LocusKind::Generic,
        }
    }

    /// The reads a record names, as the writer's live set would hold them.
    fn the_live_reads_of(record: &SampleLocusObservations) -> LiveSet {
        let mut ids: Vec<ChainId> = record
            .observations
            .iter()
            .flat_map(|observation| observation.chain_ids.iter().copied())
            .collect();
        as_a_read_set(&mut ids);
        LiveSet::from_sorted_ids(ids)
    }

    /// Encode, decode, and hand back what came out together with the bytes that went in.
    fn round_trip(record: &SampleLocusObservations) -> (DecodedRecordBody, Vec<u8>) {
        let mut bytes = Vec::new();
        encode_record_body(record, &mut bytes);
        let decoded = decode_record_body(
            &bytes,
            record.region,
            &the_live_reads_of(record),
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect("what this encoder wrote, this decoder reads");
        (decoded, bytes)
    }

    /// Decode with the layout this build writes — the form every hostile-input test wants.
    ///
    /// **With no reads live**, which is what a reader meets when the bytes are not a record this
    /// writer produced.
    fn decoded(bytes: &[u8]) -> Result<DecodedRecordBody, RecordDecodeError> {
        decode_record_body(
            bytes,
            a_region(1, 1),
            &LiveSet::new(),
            &RecordLayout::as_this_build_writes_it(),
        )
    }

    /// A manifest that declares `fields` and is otherwise the default writer's.
    fn a_manifest_declaring(fields: Vec<FieldSpec>) -> Manifest {
        Manifest {
            genomic_block_size_bp: DEFAULT_GENOMIC_BLOCK_SIZE_BP,
            block_byte_ceiling: None,
            look_back_window_log: DEFAULT_LOOK_BACK_WINDOW_LOG,
            fields,
        }
    }

    // -----------------------------------------------------------------
    // The round trip
    // -----------------------------------------------------------------

    /// The whole of C1: what goes in comes back, field for field, and the decoder's own
    /// count of bytes read is the length of what the encoder wrote — which is the number
    /// C2's head will be checked against.
    #[test]
    fn a_generic_record_round_trips_field_for_field() {
        let written = a_rich_record();
        let (decoded, bytes) = round_trip(&written);
        assert_eq!(decoded.record, written);
        assert_eq!(
            decoded.bytes_read,
            bytes.len(),
            "the decoder stopped where the encoder did"
        );
    }

    /// **The bytes, spelled out.** Every other test here goes through the encoder and the
    /// decoder together, so a field moved in both stays green while the manifest keeps
    /// declaring the old order and every psp already written is misread. This is the one
    /// assertion in the module that does not come from the code it is testing.
    ///
    /// **Never regenerate it to make a test pass.** A change here is a change to the format:
    /// the manifest in [`BODY_FIELDS`] has to move with it and the version has to rise.
    ///
    /// ⚠ **It changed at Milestone E4**, deliberately: two body fields arrived — which
    /// observation's reads are derived rather than stored, and each stored observation's list —
    /// and [`BODY_FIELDS`] moved with them. The version did not rise, for the reason
    /// [`RECORD_HEAD_FIELDS`] already gives about its own additions: **no psp file exists yet**,
    /// so there is nothing to be incompatible with. From Milestone F, which writes files that
    /// outlive a build, a change here costs a version.
    #[test]
    fn the_fixture_encodes_to_these_exact_bytes() {
        let mut bytes = Vec::new();
        encode_record_body(&a_rich_record(), &mut bytes);
        assert_eq!(
            bytes,
            vec![
                // reference-bases: seven bases over a seven-base region
                7, b'A', b'C', b'G', b'T', b'A', b'C', b'G', //
                3,    // observation-count
                3,    // residual-observation: the count itself, so no list is derived — this
                // fixture's observations name no reads at all
                // — observation 1 —
                7, b'A', b'C', b'G', b'T', b'A', b'C', b'G', // observation-bases
                0,    // witness-run-count: a complete witness
                1,    // read-group
                137, 1,  // reads-showing-the-sequence
                61, // reads-on-the-forward-strand
                255, 159, 21, // summed-log-error: zig-zag of −174,080 steps
                156, 64, // mapq-sum
                144, 141, 30, // mapq-sum-of-squares
                44, // reads-starting-left-of-the-locus
                0,  // observation-reads: this observation names none
                // — observation 2 —
                14, b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T', b'A', b'C', b'G', b'T', b'A',
                b'C', // observation-bases
                2,    // witness-run-count
                0, 2, // run at 0, two positions
                5, 2, // run at 5, two positions
                7, // read-group
                3, // reads-showing-the-sequence
                2, // reads-on-the-forward-strand
                255, 15, // summed-log-error: zig-zag of −1,024 steps
                180, 1, // mapq-sum
                176, 84, // mapq-sum-of-squares
                1,  // reads-starting-left-of-the-locus
                0,  // observation-reads: none
                // — observation 3 —
                0, // observation-bases: empty
                1, // witness-run-count
                4, 2,  // run at 4, two positions
                7,  // read-group
                1,  // reads-showing-the-sequence
                1,  // reads-on-the-forward-strand
                13, // summed-log-error: zig-zag of −7 steps
                60, // mapq-sum
                144, 28, // mapq-sum-of-squares
                2,  // reads-starting-left-of-the-locus
                0,  // observation-reads: none
                // — the record's own tail —
                12, // reads-without-observation
                177, 2, // reads-discarded-by-the-depth-cap
                0, // locus-kind: Generic
            ]
        );
    }

    /// A record's reference bases are in the body, so decoding one needs no reference on hand
    /// — four bases over a one-base region, which nothing about the region implies. It is also
    /// the only shape that catches an observation count written conditionally: with
    /// observations present, an omitted zero cannot show.
    ///
    /// *Why they are stored at all is spec `psp_record_encoding.md` §4 and the C1 report §2.1.*
    #[test]
    fn a_record_with_no_observations_round_trips_and_needs_no_reference() {
        let written = SampleLocusObservations {
            region: a_region(1_000, 1_000),
            reference_bases: b"ACGT".to_vec(),
            observations: Vec::new(),
            reads_without_observation: 4,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        let (decoded, bytes) = round_trip(&written);
        assert_eq!(&*decoded.record.reference_bases, b"ACGT");
        assert_eq!(decoded.record, written);
        assert_eq!(decoded.bytes_read, bytes.len());
    }

    /// **Every kind the encoder can write, this decoder reads back as itself.** The inner
    /// `match` is what binds the list to the enum: a kind added to [`LocusKind`] is a compile
    /// error here as well as in `put_kind`, so a tag cannot be given a write side and left
    /// without a read side — which compiles and passes every round-trip test otherwise,
    /// because each of those names one kind.
    #[test]
    fn every_locus_kind_round_trips() {
        let kinds = vec![
            LocusKind::Generic,
            LocusKind::Ssr(SsrDetail {
                motif: Motif::new(b"AT").expect("a dinucleotide is a motif"),
                // The flanks differ in length because a left flank is clamped at a contig's
                // start, so a record can carry a short one and a long one.
                left_flank: b"GG".to_vec().into_boxed_slice(),
                right_flank: b"CCCCCCCCCC".to_vec().into_boxed_slice(),
            }),
            LocusKind::SsrBundle,
        ];
        for kind in &kinds {
            match kind {
                // Exhaustive on purpose: a kind added to `LocusKind` is a compile error until
                // it is in the list above.
                LocusKind::Generic | LocusKind::Ssr(_) | LocusKind::SsrBundle => {}
            }
        }
        for kind in kinds {
            let mut written = a_rich_record();
            written.kind = kind;
            let (decoded, _) = round_trip(&written);
            assert_eq!(decoded.record, written);
        }
    }

    /// **The tags as the bytes spell them.** A tag is on disk in every file already written, so
    /// this is what makes renumbering one a failure rather than a silent reinterpretation of
    /// every tract record in the field. The round-trip tests cannot see it: they are symmetric
    /// and pass under any consistent renumbering.
    #[test]
    fn the_locus_kind_tags_are_the_numbers_the_files_carry() {
        let tract = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide is a motif"),
            left_flank: Box::from(&b""[..]),
            right_flank: Box::from(&b""[..]),
        });
        for (kind, tag) in [
            (LocusKind::Generic, 0u8),
            (tract, 1),
            (LocusKind::SsrBundle, 2),
        ] {
            let mut out = Vec::new();
            put_kind(&mut out, &kind);
            assert_eq!(out[0], tag, "{kind:?}");
        }
    }

    /// **The chain ids round-trip exactly, which is what Checkpoint E asks for.**
    ///
    /// ⚠ This test was `chain_ids_come_back_empty_and_nothing_else_changes` from Milestone C1 —
    /// *"the day E lands, this test is what fails if the ids are still being dropped"*. E4 is
    /// that day, and this is what it became.
    ///
    /// The record below is the interesting shape: three observations naming reads, so the writer
    /// derives the largest list and stores the other two.
    #[test]
    fn the_chain_ids_round_trip_exactly() {
        let mut written = a_rich_record();
        written.observations[0].chain_ids = vec![4, 17, 900_001];
        written.observations[1].chain_ids = vec![41];
        written.observations[2].chain_ids = vec![70_000];
        for observation in &mut written.observations {
            observation.num_obs = observation.chain_ids.len() as u32;
        }

        let (decoded, bytes) = round_trip(&written);
        assert_eq!(decoded.record, written, "every list, exactly");

        // **And the saving is real**, measured against the same record with its residual made
        // underivable. ⚠ Two assertions stood here that a writer deriving *nothing* satisfied
        // just as well: one compared two direct `encode_read_list` calls and never touched the
        // writer's choice, the other only said a record with fewer identifiers is shorter.
        assert_eq!(
            residual_observation_of(&written, &mut ResidualScratch::default()),
            0,
            "the largest list is the one derived"
        );
        let mut nothing_derived = written.clone();
        // A read named by two observations: the derivation cannot reproduce it, so every list is
        // stored — the same record, the same identifiers, one list more on disk.
        nothing_derived.observations[1].chain_ids.push(4);
        nothing_derived.observations[1].chain_ids.sort_unstable();
        nothing_derived.observations[1].num_obs += 1;
        assert_eq!(
            residual_observation_of(&nothing_derived, &mut ResidualScratch::default()),
            nothing_derived.observations.len(),
            "no observation can be derived from this one"
        );
        let mut stored = Vec::new();
        encode_record_body(&nothing_derived, &mut stored);

        let mut the_residuals_list = Vec::new();
        encode_read_list(&written.observations[0].chain_ids, &mut the_residuals_list);
        assert!(
            the_residuals_list.len() > 1,
            "the residual has to name reads for its absence to be worth measuring"
        );
        assert!(
            bytes.len() + the_residuals_list.len() <= stored.len() + 2,
            "the derived body is about the residual's list shorter: {} against {}, and that \
             list is {} bytes",
            bytes.len(),
            stored.len(),
            the_residuals_list.len()
        );
    }

    /// The negative sums are the real ones — a sum of log error probabilities is at most
    /// zero — and the extremes are what a zig-zag encoding gets wrong if it is not one.
    #[test]
    fn a_summed_log_error_round_trips_negative_and_at_its_extremes() {
        for steps in [0, -1, 1, -4_096, i64::MIN, i64::MAX, -123_456_789] {
            let mut written = a_rich_record();
            written.observations[0].q_sum = SummedLogError::from_steps(steps);
            let (decoded, _) = round_trip(&written);
            assert_eq!(
                decoded.record.observations[0].q_sum.steps(),
                steps,
                "a summed log-error of {steps} steps"
            );
        }
    }

    /// The squared sum is the record's one 64-bit count, and no other fixture exceeds what a
    /// `u32` holds — so nothing else would notice the width it is read at shrinking, and a
    /// sample at several hundred reads a position would stop decoding.
    #[test]
    fn a_squared_mapping_quality_sum_past_u32_round_trips() {
        let mut written = a_rich_record();
        written.observations[0].mapq_sum_sq = u64::from(u32::MAX) + 1_000;
        let (decoded, _) = round_trip(&written);
        assert_eq!(decoded.record, written);
    }

    /// **Three runs, further apart than the fixture's two.** The fixture already carries a
    /// holed witness, so any round-trip test notices runs merged into their span; what this
    /// one adds is a run count above two and a start needing more than one byte, checked run
    /// by run rather than through whole-record equality.
    #[test]
    fn a_partial_witness_with_three_runs_round_trips_run_for_run() {
        let mut written = a_rich_record();
        written.observations[0].read_witness = ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 1), (3, 4), (9, 40)])
                .expect("three separated runs"),
        };
        let (decoded, _) = round_trip(&written);
        match &decoded.record.observations[0].read_witness {
            ReadWitness::Partial { positions } => {
                assert_eq!(
                    positions.runs().collect::<Vec<_>>(),
                    vec![(0, 1), (3, 4), (9, 40)]
                );
            }
            other => panic!("expected a partial witness, got {other:?}"),
        }
    }

    /// A record widened by a long deletion carries a long reference and long sequences, so its
    /// byte strings need a length prefix of more than one byte — and its bases are bytes rather
    /// than an alphabet, which no other fixture says.
    #[test]
    fn a_long_non_acgt_sequence_round_trips_through_a_multi_byte_length_prefix() {
        let mut written = a_rich_record();
        written.reference_bases = vec![b'N'; 300];
        written.observations[0].bases = (0u8..=255).collect::<Vec<_>>();
        let (decoded, bytes) = round_trip(&written);
        assert_eq!(decoded.record, written);
        assert_eq!(decoded.bytes_read, bytes.len());
    }

    /// **Every scalar in the fixture is non-zero, and that is what makes the round trip a
    /// test.** A field added to the record and filled in here with the value
    /// `decode_record_body` invents for it — which is what the compiler's own error message
    /// invites — would round-trip by coincidence and be dropped from every psp written from
    /// that day. This fails first, at the fixture.
    #[test]
    fn the_fixture_leaves_no_field_at_a_value_a_decoder_could_invent() {
        let record = a_rich_record();
        assert!(!record.reference_bases.is_empty());
        assert_ne!(record.reads_without_observation, 0);
        assert_ne!(record.reads_discarded_by_cap, 0);
        for (at, observation) in record.observations.iter().enumerate() {
            assert_ne!(observation.read_group.get(), 0, "observation {at}");
            assert_ne!(observation.num_obs, 0, "observation {at}");
            assert_ne!(observation.num_fwd, 0, "observation {at}");
            assert_ne!(observation.q_sum, SummedLogError::NONE, "observation {at}");
            assert_ne!(observation.mapq_sum, 0, "observation {at}");
            assert_ne!(observation.mapq_sum_sq, 0, "observation {at}");
            assert_ne!(observation.placed_left, 0, "observation {at}");
        }
    }

    // -----------------------------------------------------------------
    // The manifest a file carries, and this reader's check against it
    // -----------------------------------------------------------------

    /// **The names a file carries on disk, spelled out.** Every other manifest test compares
    /// [`BODY_FIELDS`] with itself; this one is the copy that has to be edited on purpose, so a
    /// rename that would change what a written file declares cannot pass unnoticed.
    #[test]
    fn record_fields_declares_the_names_a_written_file_carries() {
        let names: Vec<String> = record_fields()
            .into_iter()
            .map(|field| field.name.0)
            .collect();
        assert_eq!(
            names,
            vec![
                "position-offset",
                "reference-span",
                "non-reference-reads",
                "record-body-byte-count",
                "chain-id-changes",
                "reference-bases",
                "observation-count",
                "residual-observation",
                "residual-read-count",
                "observation-bases",
                "witness-run-count",
                "witness-run-start",
                "witness-run-length",
                "read-group",
                "reads-showing-the-sequence",
                "reads-on-the-forward-strand",
                "summed-log-error",
                "mapq-sum",
                "mapq-sum-of-squares",
                "reads-starting-left-of-the-locus",
                "observation-reads",
                "reads-without-observation",
                "reads-discarded-by-the-depth-cap",
                "locus-kind",
                "repeat-motif",
                "repeat-left-flank",
                "repeat-right-flank",
            ]
        );

        // **And the two chain-id fields are declared under the schemes that describe their
        // bytes**, which is what a *later* reader stepping over one is driven by. ⚠ They shared a
        // scheme when the second landed: `observation-reads` is one counted run of gaps and
        // `chain-id-changes` is two, so a reader measuring a list under the changes' rule reads
        // its identifiers as a departure count and then asks for an arrival count that is not
        // there. Nothing round-trips a manifest here, so this is the assertion that holds it.
        let declared: Vec<(String, FieldEncoding)> = record_fields()
            .into_iter()
            .map(|field| (field.name.0, field.encoding))
            .collect();
        assert!(
            declared.contains(&(
                "chain-id-changes".to_string(),
                FieldEncoding::ChainIdChanges
            )),
            "the head's live-set changes are two counted runs"
        );
        assert!(
            declared.contains(&("observation-reads".to_string(), FieldEncoding::ChainIdList)),
            "an observation's own reads are one counted run"
        );
    }

    /// What a writer puts in the header is what this reader checks it against, so a file this
    /// version wrote has nothing in it this version does not know.
    #[test]
    fn the_manifest_a_writer_declares_is_the_one_this_reader_checks() {
        let layout = RecordLayout::from_manifest(&a_manifest_declaring(record_fields()))
            .expect("the fields this writer declares are the fields this reader reads");
        assert_eq!(layout, RecordLayout::as_this_build_writes_it());
        assert_eq!(layout.unknown_field_count(), 0);
    }

    /// The manifest a writer declares has to be one the header can carry. The module has been
    /// bitten here before: the A1+A2 review found `encode` writing headers its own reader
    /// refused, twice over.
    #[test]
    fn the_fields_a_writer_declares_survive_a_header_round_trip() {
        let header = Header {
            format_version: FORMAT_VERSION,
            sample: "SRR7279481".to_string(),
            reference: ReferenceIdentity {
                name: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                md5: None,
            },
            contigs: vec![ContigIdentity {
                name: "SL4.0ch01".to_string(),
                length: 90_863_682,
                md5: None,
            }],
            writer: WriterProvenance {
                tool: "ng".to_string(),
                version: "0.1.0".to_string(),
                subcommand: "pileup".to_string(),
                input_alignments: vec!["SRR7279481.cram".to_string()],
                input_reference: "S_lycopersicum_chromosomes.4.00.fa".to_string(),
                command_line: "ng pileup --sample SRR7279481".to_string(),
                parameters: std::collections::BTreeMap::new(),
                created: "2026-08-26T09:15:00Z".parse().expect("an RFC 3339 stamp"),
            },
            manifest: a_manifest_declaring(record_fields()),
        };

        let bytes = header
            .encode()
            .expect("a writer's own manifest must be a legal header");
        let (back, _) = Header::decode(&bytes, std::path::Path::new("SRR7279481.psp"))
            .expect("and its own reader must read it");
        assert_eq!(back.manifest.fields, record_fields());
        assert_eq!(
            RecordLayout::from_manifest(&back.manifest).expect("and drive this reader"),
            RecordLayout::as_this_build_writes_it()
        );
    }

    /// **The step is the type's, not the file's, and this is the seam between them.** The
    /// declaration exists so a reader can turn the stored integer back into a quantity; a
    /// writer declaring anything else would put every error sum in the file out by the ratio.
    #[test]
    fn the_declared_step_is_the_types_step_and_is_four_thousand_and_ninety_six() {
        assert_eq!(SummedLogError::STEPS_PER_NAT, 4_096);
        let summed_log_error = record_fields()
            .into_iter()
            .find(|field| field.name.0 == "summed-log-error")
            .expect("the field is declared");
        assert_eq!(
            summed_log_error.encoding,
            FieldEncoding::FixedPoint {
                steps_per_unit: 4_096
            }
        );
    }

    /// A file declaring a coarser step holds integers that mean something else, and reading
    /// them as this type's would divide every error sum in the file by four — so the file is
    /// refused rather than rescaled. The wrong value is derived from the right one, so this
    /// cannot be satisfied by a coincidence of literals.
    #[test]
    fn a_manifest_that_declares_a_different_step_for_the_summed_log_error_is_refused() {
        let quarter_of_the_real_step = SummedLogError::STEPS_PER_NAT as u32 / 4;
        let mut fields = record_fields();
        let position = fields
            .iter()
            .position(|field| field.name.0 == "summed-log-error")
            .expect("the field is declared");
        fields[position].encoding = FieldEncoding::FixedPoint {
            steps_per_unit: quarter_of_the_real_step,
        };

        match RecordLayout::from_manifest(&a_manifest_declaring(fields)) {
            Err(RecordLayoutError::WrongEncoding { field, found, .. }) => {
                assert_eq!(field.0, "summed-log-error");
                assert_eq!(
                    found,
                    FieldEncoding::FixedPoint {
                        steps_per_unit: quarter_of_the_real_step
                    }
                );
            }
            other => panic!("expected WrongEncoding, got {other:?}"),
        }
    }

    /// Every field, one at a time: renamed, dropped, or swapped with the field next to it. All
    /// three are files whose records would otherwise decode into plausible values rather than
    /// failing — two adjacent variable-length integer fields swapped is two numbers exchanged
    /// and no error at all. The rename is asserted by variant, because which refusal it is is
    /// what tells whoever holds the file which of its fields is wrong.
    #[test]
    fn a_manifest_that_renames_reorders_or_drops_a_field_is_refused() {
        let declared = record_fields();
        for position in 0..declared.len() {
            let mut renamed = declared.clone();
            renamed[position].name = FieldName("something-else".to_string());
            match RecordLayout::from_manifest(&a_manifest_declaring(renamed)) {
                Err(RecordLayoutError::UnexpectedField {
                    position: at,
                    expected,
                    found,
                }) => {
                    assert_eq!(at, position);
                    assert_eq!(expected.0, declared[position].name.0);
                    assert_eq!(found.0, "something-else");
                }
                other => panic!("field {position} renamed: expected UnexpectedField, {other:?}"),
            }

            let mut dropped = declared.clone();
            dropped.remove(position);
            assert!(
                RecordLayout::from_manifest(&a_manifest_declaring(dropped)).is_err(),
                "field {position} dropped"
            );

            if position + 1 < declared.len() {
                let mut swapped = declared.clone();
                swapped.swap(position, position + 1);
                assert!(
                    RecordLayout::from_manifest(&a_manifest_declaring(swapped)).is_err(),
                    "fields {position} and {} swapped",
                    position + 1
                );
            }
        }
    }

    /// A manifest that ends early names the field it ended before, so whoever reads the
    /// message knows what the file is missing rather than that it is short.
    #[test]
    fn a_manifest_that_ends_early_names_the_field_it_stopped_before() {
        let mut fields = record_fields();
        fields.truncate(2);
        match RecordLayout::from_manifest(&a_manifest_declaring(fields)) {
            Err(RecordLayoutError::MissingField { position, expected }) => {
                assert_eq!(position, 2);
                assert_eq!(expected.0, record_fields()[2].name.0);
            }
            other => panic!("expected MissingField, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------
    // A field a later writer added
    // -----------------------------------------------------------------

    /// **What the owner's ruling of 2026-08-26 asks of this codec, tested.** The window's GC
    /// fraction and its mean coverage are not computed anywhere in ng yet, and adding either
    /// later must not break a reader already in the field. A later writer declares the new
    /// field after the ones here and writes it after the body's last byte; this reader walks
    /// past it on the encoding's own self-measurement, gets the record it always got, and
    /// reports the whole body as read.
    ///
    /// Every encoding in the closed set, because each measures itself differently.
    #[test]
    fn a_field_this_reader_does_not_know_is_walked_past_when_it_comes_last() {
        let record = a_rich_record();
        let mut body = Vec::new();
        encode_record_body(&record, &mut body);

        let later_versions = [
            (FieldEncoding::Varint, vec![0xac, 0x02]),
            (FieldEncoding::SignedVarint, vec![0xd7, 0x04]),
            (
                FieldEncoding::FixedWidthInteger { width_bytes: 4 },
                vec![1, 2, 3, 4],
            ),
            (FieldEncoding::IeeeFloat { width_bytes: 8 }, vec![9; 8]),
            // Multi-byte on purpose: a one-byte value cannot tell a varint skip from a
            // one-byte skip, and this is the only encoding whose framing is not obvious.
            (
                FieldEncoding::FixedPoint {
                    steps_per_unit: 100,
                },
                vec![0xc0, 0x1f],
            ),
            (
                FieldEncoding::LengthPrefixedBytes,
                vec![3, b'a', b'b', b'c'],
            ),
        ];

        for (encoding, written) in later_versions {
            let mut fields = record_fields();
            fields.push(FieldSpec {
                name: FieldName("a-field-from-a-later-writer".to_string()),
                encoding,
            });
            let layout = RecordLayout::from_manifest(&a_manifest_declaring(fields))
                .expect("an unknown field at the end is not a reason to refuse the file");
            assert_eq!(layout.unknown_field_count(), 1);

            let mut newer_body = body.clone();
            newer_body.extend_from_slice(&written);

            let decoded = decode_record_body(
                &newer_body,
                record.region,
                &the_live_reads_of(&record),
                &layout,
            )
            .unwrap_or_else(|refused| panic!("{encoding:?} was not walked past: {refused}"));
            assert_eq!(decoded.record, record, "{encoding:?}");
            assert_eq!(
                decoded.bytes_read,
                newer_body.len(),
                "{encoding:?} — the whole body has to be accounted for, or the next record \
                 starts in the middle of this one"
            );
        }
    }

    /// Two later fields, of different encodings, so the order they are walked in matters — with
    /// one unknown field every ordering is the same ordering, and the two the design actually
    /// anticipates are two.
    #[test]
    fn two_fields_this_reader_does_not_know_are_walked_past_in_the_order_declared() {
        let record = a_rich_record();
        let mut body = Vec::new();
        encode_record_body(&record, &mut body);

        let mut fields = record_fields();
        fields.push(FieldSpec {
            name: FieldName("window-gc-percent".to_string()),
            encoding: FieldEncoding::FixedWidthInteger { width_bytes: 4 },
        });
        fields.push(FieldSpec {
            name: FieldName("window-mean-coverage".to_string()),
            encoding: FieldEncoding::LengthPrefixedBytes,
        });
        let layout = RecordLayout::from_manifest(&a_manifest_declaring(fields))
            .expect("two unknown fields at the end are fine");
        assert_eq!(layout.unknown_field_count(), 2);

        body.extend_from_slice(&[7, 7, 7, 7, 2, b'h', b'i']);
        let decoded =
            decode_record_body(&body, record.region, &the_live_reads_of(&record), &layout)
                .expect("both are walked past");
        assert_eq!(decoded.record, record);
        assert_eq!(decoded.bytes_read, body.len());
    }

    /// A later writer's trailing field cut off by a block boundary is a short body — the class
    /// Milestone D reads more bytes for — and not a panic.
    #[test]
    fn an_unknown_trailing_field_cut_short_is_truncated_without_panicking() {
        let record = a_rich_record();
        let mut body = Vec::new();
        encode_record_body(&record, &mut body);
        let mut fields = record_fields();
        fields.push(FieldSpec {
            name: FieldName("window-mean-coverage".to_string()),
            encoding: FieldEncoding::IeeeFloat { width_bytes: 8 },
        });
        let layout = RecordLayout::from_manifest(&a_manifest_declaring(fields))
            .expect("an unknown field at the end is fine");

        let whole = body.len();
        body.extend_from_slice(&[1, 2, 3]);
        for cut in whole..body.len() {
            match decode_record_body(
                &body[..cut],
                record.region,
                &the_live_reads_of(&record),
                &layout,
            ) {
                Err(RecordDecodeError::Truncated { .. }) => {}
                other => panic!("{cut} bytes: expected Truncated, got {other:?}"),
            }
        }
    }

    // -----------------------------------------------------------------
    // Damaged and hostile bodies
    // -----------------------------------------------------------------

    /// A body cut at any point is refused at that point, **as `Truncated` and never as
    /// damage** — the split Milestone D's restartable read branches on, and one no other test
    /// in this module asserts. Nothing in the decoder indexes past what it holds.
    #[test]
    fn a_body_cut_short_is_truncated_at_every_cut_and_never_malformed() {
        let record = a_rich_record();
        let mut whole = Vec::new();
        encode_record_body(&record, &mut whole);

        for cut in 0..whole.len() {
            match decode_record_body(
                &whole[..cut],
                record.region,
                &the_live_reads_of(&record),
                &RecordLayout::as_this_build_writes_it(),
            ) {
                Err(RecordDecodeError::Truncated { .. }) => {}
                other => panic!(
                    "{cut} bytes of a record is a short buffer, which Milestone D reads more \
                     bytes for — not damage; got {other:?}"
                ),
            }
        }
        assert!(
            decode_record_body(
                &whole,
                record.region,
                &the_live_reads_of(&record),
                &RecordLayout::as_this_build_writes_it()
            )
            .is_ok(),
            "and the whole body still reads"
        );
    }

    /// A length or a count no record body could hold is **damage and never a short buffer**:
    /// reported as truncation it would ask Milestone D's reader to grow its buffer to a
    /// terabyte instead of refusing the file.
    #[test]
    fn a_length_no_body_could_hold_is_malformed_and_never_truncated() {
        let a_length_no_body_can_have = MOST_BYTES_A_BODY_CAN_DECLARE + 1;
        let cases: [(&str, Vec<u8>); 3] = [
            ("a reference longer than a body", {
                let mut body = Vec::new();
                encode_u64_leb128(a_length_no_body_can_have, &mut body);
                body
            }),
            ("more observations than a body could hold", {
                let mut body = vec![0u8];
                encode_u64_leb128(u64::MAX, &mut body);
                body
            }),
            ("more witness runs than a body could hold", {
                let mut body = vec![0u8, 1, 1, 0];
                encode_u64_leb128(u64::MAX, &mut body);
                body
            }),
        ];
        for (what, body) in cases {
            match decoded(&body) {
                Err(RecordDecodeError::Malformed { .. }) => {}
                other => panic!("{what}: expected Malformed, got {other:?}"),
            }
        }
    }

    /// A length prefix past what an address space holds, read at a non-zero offset — the case
    /// where a cursor that added without checking would land *behind* itself and index a slice
    /// backwards.
    #[test]
    fn a_length_prefix_past_the_address_space_is_refused_not_indexed() {
        let mut huge = Vec::new();
        encode_u64_leb128(u64::MAX, &mut huge);
        for prefix in [vec![], vec![3u8, b'A', b'C', b'G', 1]] {
            let mut body = prefix;
            body.extend_from_slice(&huge);
            body.extend_from_slice(&[0u8; 8]);
            assert!(
                decoded(&body).is_err(),
                "a length of 2^64 − 1 must be refused, not indexed"
            );
        }
    }

    /// **A count the file declares never sizes an allocation on its own.** Each fixture is a
    /// handful of bytes claiming more observations — or more witnessed runs — than there are
    /// bytes in the universe. Without the bound the reserve is asked for before a single entry
    /// is read: this test would fail with a capacity overflow rather than pass, and a merely
    /// large count would quietly reserve gigabytes.
    #[test]
    fn a_declared_count_larger_than_the_body_never_reaches_the_allocator() {
        let mut observations_claimed = vec![0u8];
        encode_u64_leb128(u64::MAX, &mut observations_claimed);
        assert!(decoded(&observations_claimed).is_err());

        let mut runs_claimed = vec![0u8, 1, 1, 0];
        encode_u64_leb128(u64::MAX, &mut runs_claimed);
        assert!(decoded(&runs_claimed).is_err());
    }

    /// **The exact count `read_count` accepts and the first it refuses**, both sides of the
    /// boundary, on every `least_bytes_each` the module uses.
    ///
    /// **It exists because the boundary was never pinned, and a review proposed moving it.** The
    /// candidate was to replace the division with `declared * L > M`, which is the same statement
    /// for positive integers — the truncation in the division is exactly the slack the
    /// multiplication keeps. **Measured and not taken:** callgrind on the deep full walk, the same
    /// binary twice each way, put the multiplication at 6,762,039,538 instructions against the
    /// division's 6,762,760,740 — **721,200 saved out of 6.76 billion, 0.011 %**, because LLVM
    /// already folds the division wherever `least_bytes_each` is a constant, which is every call
    /// site in this file. A non-obvious arithmetic identity inside a guard is not worth 0.011 %.
    ///
    /// The test stays, because nothing pinned this boundary before and the next person to have
    /// the idea should find it held rather than have to re-derive it.
    #[test]
    fn the_largest_count_a_body_could_hold_is_accepted_and_the_next_one_is_not() {
        for least_bytes_each in [1, LEAST_BYTES_PER_WITNESS_RUN, LEAST_BYTES_PER_OBSERVATION] {
            let most = MOST_BYTES_A_BODY_CAN_DECLARE / least_bytes_each as u64;
            for (declared, accepted) in [(most, true), (most + 1, false)] {
                let mut bytes = Vec::new();
                put_varint(&mut bytes, declared);
                let mut reader = FieldReader::new(&bytes);
                let got = reader.read_count(OBSERVATION_COUNT, least_bytes_each);
                assert_eq!(
                    got.is_ok(),
                    accepted,
                    "{declared} entries at {least_bytes_each} bytes each: the ceiling is {most}"
                );
                if accepted {
                    assert_eq!(got.expect("accepted"), declared);
                }
            }
        }
    }

    /// The reserve a declared count is allowed to make is what the bytes left could hold,
    /// capped absolutely — and "not the declared count" is not the same as "bounded". A guard
    /// that ignored the bytes and capped at some large constant would pass the hostile-body
    /// test above and fail here.
    #[test]
    fn entries_to_reserve_is_bounded_by_the_bytes_left_and_by_a_ceiling() {
        assert_eq!(
            entries_to_reserve(u64::MAX, LEAST_BYTES_PER_OBSERVATION, 0),
            0
        );
        assert_eq!(
            entries_to_reserve(u64::MAX, LEAST_BYTES_PER_OBSERVATION, 20),
            2
        );
        assert_eq!(
            entries_to_reserve(1_000_000, LEAST_BYTES_PER_OBSERVATION, 11),
            1
        );
        assert_eq!(
            entries_to_reserve(3, LEAST_BYTES_PER_OBSERVATION, 900),
            3,
            "an honest count is kept"
        );
        assert_eq!(
            entries_to_reserve(u64::MAX, LEAST_BYTES_PER_WITNESS_RUN, 1 << 20),
            MOST_ENTRIES_RESERVED,
            "a megabyte of body still reserves no more than the ceiling"
        );
    }

    /// **The bound is measured, not recalled.** It is what a declared count is trusted against,
    /// so a field added to an observation moves it and this is what says so.
    #[test]
    fn the_least_an_observation_and_a_run_can_cost_is_what_the_bounds_say() {
        let empty = SampleLocusObservations {
            region: a_region(1, 1),
            reference_bases: b"".to_vec(),
            observations: Vec::new(),
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        let cheapest_observation = SequenceObservation {
            bases: b"".to_vec(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: 0,
            num_fwd: 0,
            q_sum: SummedLogError::NONE,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        };

        let mut with_one = empty.clone();
        with_one.observations.push(cheapest_observation.clone());
        let mut with_a_run = empty.clone();
        with_a_run.observations.push(SequenceObservation {
            read_witness: ReadWitness::Partial {
                positions: WitnessedLocusPositions::one_run_from_offset_and_length(0, 1)
                    .expect("one run of one position"),
            },
            ..cheapest_observation.clone()
        });

        let bytes_of = |record: &SampleLocusObservations| {
            let mut bytes = Vec::new();
            encode_record_body(record, &mut bytes);
            bytes.len()
        };
        // **An observation costs the constant when its own list is derived rather than stored**,
        // and one byte more when it is stored — which is the case for a record whose reads the
        // writer cannot derive, and for the one below, whose observation names none. The constant
        // is a lower bound over both, because that is what bounding a declared count needs.
        let mut with_one_that_is_derived = empty.clone();
        with_one_that_is_derived
            .observations
            .push(SequenceObservation {
                num_obs: 1,
                chain_ids: vec![5],
                ..cheapest_observation.clone()
            });
        assert_eq!(
            bytes_of(&with_one_that_is_derived) - bytes_of(&empty),
            LEAST_BYTES_PER_OBSERVATION + 1,
            "a derived observation costs the constant plus the record's residual-read-count, \
             which is per record and only written when something is derived"
        );
        assert_eq!(
            bytes_of(&with_one) - bytes_of(&empty),
            LEAST_BYTES_PER_OBSERVATION + 1,
            "an observation whose list is stored costs its count byte too"
        );
        assert_eq!(
            bytes_of(&with_a_run) - bytes_of(&with_one),
            LEAST_BYTES_PER_WITNESS_RUN
        );
    }

    /// Bytes after the record are not the record's business: a decoder stops where the record
    /// ends and says so, because in a block the next record starts there.
    #[test]
    fn a_body_followed_by_more_bytes_stops_where_the_record_ends() {
        let record = a_rich_record();
        let mut stream = Vec::new();
        encode_record_body(&record, &mut stream);
        let ends_at = stream.len();
        stream.extend_from_slice(&[0xff; 64]);

        let decoded = decode_record_body(
            &stream,
            record.region,
            &the_live_reads_of(&record),
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect("the record is complete; what follows it is another record's");
        assert_eq!(decoded.record, record);
        assert_eq!(decoded.bytes_read, ends_at);
    }

    /// A run covering no position is a set [`WitnessedLocusPositions`] refuses to hold, and a
    /// witness with no runs at all already means [`ReadWitness::Complete`] — so a body saying
    /// both is damaged, and is refused rather than turned into a complete witness.
    #[test]
    fn a_witness_run_covering_no_position_is_refused() {
        // An empty reference; one observation; empty bases; one run; start 4, length 0.
        let body = vec![0u8, 1, 1, 0, 1, 4, 0];
        match decoded(&body) {
            Err(RecordDecodeError::Malformed { field, .. }) => {
                assert_eq!(field, "witness-run-length");
            }
            other => panic!("expected a malformed witness run, got {other:?}"),
        }
    }

    /// **Runs this encoder could not have written are refused, not repaired.**
    /// `WitnessedLocusPositions::from_half_open_runs` sorts and merges whatever it is handed,
    /// so an unsorted or overlapping pair would otherwise decode into a *different valid
    /// record* — with the byte count still right, so even C2's length check would not see it.
    /// A witness is what read depth is raised over, so a merged run credits a read with
    /// positions it never saw.
    #[test]
    fn witness_runs_that_do_not_ascend_or_that_touch_are_refused() {
        // An empty reference, one observation, empty bases, two runs, then the rest as zeroes.
        let out_of_order = vec![0u8, 1, 1, 0, 2, 5, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let overlapping = vec![0u8, 1, 1, 0, 2, 0, 5, 3, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let touching = vec![0u8, 1, 1, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        for (what, body) in [
            ("runs out of order", out_of_order),
            ("runs that overlap", overlapping),
            ("runs that touch", touching),
        ] {
            match decoded(&body) {
                Err(RecordDecodeError::Malformed { field, reason, .. }) => {
                    assert_eq!(field, "witness", "{what}");
                    assert!(
                        reason.contains("ascend and never touch"),
                        "{what}: got {reason}"
                    );
                }
                other => panic!("{what}: expected a refusal, got {other:?}"),
            }
        }
    }

    /// A witness coordinate past what a locus position holds is damage, not a coordinate to
    /// narrow: 65,536 must not come back as 0, or a read's evidence lands at the wrong offsets
    /// inside the record with no error anywhere.
    #[test]
    fn a_witness_coordinate_past_a_locus_position_is_refused_rather_than_narrowed() {
        for (start, length, expected) in [
            (65_536u64, 4u64, "witness-run-start"),
            (4, 65_536, "witness-run-length"),
        ] {
            let mut body = vec![0u8, 1, 1, 0, 1];
            encode_u64_leb128(start, &mut body);
            encode_u64_leb128(length, &mut body);
            body.extend_from_slice(&[0u8; 8]);
            match decoded(&body) {
                Err(RecordDecodeError::Malformed { field, reason, .. }) => {
                    assert_eq!(field, expected);
                    assert!(reason.contains("65536"), "got {reason}");
                }
                other => panic!("expected {expected} refused, got {other:?}"),
            }
        }
    }

    /// A run that starts inside the range and ends outside it — the overflow the sum can have
    /// and neither of its parts can.
    #[test]
    fn a_witness_run_ending_past_the_last_locus_position_is_refused() {
        let mut body = vec![0u8, 1, 1, 0, 1];
        encode_u64_leb128(65_535, &mut body);
        encode_u64_leb128(1, &mut body);
        body.extend_from_slice(&[0u8; 8]);
        match decoded(&body) {
            Err(RecordDecodeError::Malformed { field, reason, .. }) => {
                assert_eq!(field, "witness");
                assert!(reason.contains("ends past"), "got {reason}");
            }
            other => panic!("expected an overrunning run refused, got {other:?}"),
        }
    }

    /// **A kind tag this reader does not know is `Unsupported`, not damage.** Only the major
    /// format version is enforced when a file is opened, and the unknown-field machinery exists
    /// so a file from a later minor version is readable — so the instruction here is *upgrade
    /// the reader*, the same one `header.rs` gives for a newer format. It must not be read as
    /// the kind whose tag is next to it either, or a later writer's repeat tracts become SNP
    /// candidates.
    #[test]
    fn a_locus_kind_from_a_later_writer_says_upgrade_the_reader() {
        let record = a_rich_record();
        let mut body = Vec::new();
        encode_record_body(&record, &mut body);
        let last = body.len() - 1;
        assert_eq!(
            body[last], 0,
            "the fixture's kind is Generic, which is tag 0"
        );
        body[last] = 99;

        match decode_record_body(
            &body,
            record.region,
            &the_live_reads_of(&record),
            &RecordLayout::as_this_build_writes_it(),
        ) {
            Err(RecordDecodeError::Unsupported { field, tag, .. }) => {
                assert_eq!(field, "locus-kind");
                assert_eq!(tag, 99);
            }
            other => panic!("expected Unsupported, got {other:?}"),
        }
    }

    /// A motif of no bases, or of more than a repeat unit can hold, is refused where the
    /// motif is built rather than stored as a motif that could not have been minted.
    #[test]
    fn a_stored_motif_outside_the_repeat_period_range_is_refused() {
        for motif_bases in [&b""[..], &b"ATATATA"[..]] {
            // An empty reference, no observations, no observation derived, no reads either
            // way, kind 1, then the motif and two empty flanks.
            let mut body = vec![0u8, 0, 0, 0, 0, 1];
            body.push(motif_bases.len() as u8);
            body.extend_from_slice(motif_bases);
            body.extend_from_slice(&[0, 0]);

            match decoded(&body) {
                Err(RecordDecodeError::Malformed { field, .. }) => {
                    assert_eq!(field, "repeat-motif", "for {motif_bases:?}");
                }
                other => panic!("expected a refused motif for {motif_bases:?}, got {other:?}"),
            }
        }
    }

    /// A count that does not fit the field it is read into is damage, and the message says
    /// which field, what the number was, and how far into the body it sat — not a silently
    /// truncated depth.
    #[test]
    fn a_count_too_large_for_its_field_is_refused_rather_than_narrowed() {
        // An empty reference, one observation, empty bases, a complete witness, read group 0,
        // then a read count of 2^32.
        let mut body = vec![0u8, 1, 1, 0, 0, 0];
        let at = body.len();
        encode_u64_leb128(u64::from(u32::MAX) + 1, &mut body);

        match decoded(&body) {
            Err(RecordDecodeError::Malformed {
                field,
                bytes_in,
                reason,
            }) => {
                assert_eq!(field, "reads-showing-the-sequence");
                assert!(reason.contains("4294967296"), "got {reason}");
                assert!(
                    bytes_in > at,
                    "the offset has to be past where the field began, {at}; got {bytes_in}"
                );
            }
            other => panic!("expected a refused count, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------
    // Over the whole input domain, rather than the shapes anyone thought of
    // -----------------------------------------------------------------

    fn an_arbitrary_observation() -> impl Strategy<Value = SequenceObservation> {
        (
            prop::collection::vec(any::<u8>(), 0..40),
            prop::option::of((0u16..50, 1u16..20)),
            any::<u32>(),
            any::<u32>(),
            any::<u32>(),
            any::<i64>(),
            any::<u32>(),
            any::<u64>(),
            any::<u32>(),
        )
            .prop_map(
                |(bases, run, group, obs, forward, steps, mapq, mapq_squared, left)| {
                    SequenceObservation {
                        bases: bases,
                        read_witness: match run.and_then(|(start, length)| {
                            WitnessedLocusPositions::one_run_from_offset_and_length(start, length)
                        }) {
                            Some(positions) => ReadWitness::Partial { positions },
                            None => ReadWitness::Complete,
                        },
                        read_group: ReadGroupId(group),
                        num_obs: obs,
                        num_fwd: forward,
                        q_sum: SummedLogError::from_steps(steps),
                        mapq_sum: mapq,
                        mapq_sum_sq: mapq_squared,
                        placed_left: left,
                        chain_ids: Vec::new(),
                    }
                },
            )
    }

    proptest! {
        /// The round trip is the step's whole contract, over the field combinations no hand
        /// fixture enumerates.
        #[test]
        fn any_record_round_trips_exactly(
            reference in prop::collection::vec(any::<u8>(), 0..40),
            observations in prop::collection::vec(an_arbitrary_observation(), 0..6),
            without in any::<u32>(),
            capped in any::<u32>(),
        ) {
            let written = SampleLocusObservations {
                region: a_region(1_000, 1_001),
                reference_bases: reference,
                observations,
                reads_without_observation: without,
                reads_discarded_by_cap: capped,
                kind: LocusKind::Generic,
            };
            let mut bytes = Vec::new();
            encode_record_body(&written, &mut bytes);
            let decoded = decode_record_body(
                &bytes,
                written.region,
                &the_live_reads_of(&written),
                &RecordLayout::as_this_build_writes_it(),
            )
            .expect("what this encoder wrote, this decoder reads");
            prop_assert_eq!(decoded.record, written);
            prop_assert_eq!(decoded.bytes_read, bytes.len());
        }

        /// Arbitrary bytes are a typed error or a record, never a panic — a psp on disk may be
        /// anything, and the shapes a fixture resembles are not the shapes damage takes.
        #[test]
        fn arbitrary_bytes_decode_or_are_refused_but_never_panic(
            body in prop::collection::vec(any::<u8>(), 0..80)
        ) {
            if let Ok(decoded) = decoded(&body) {
                prop_assert!(decoded.bytes_read <= body.len());
            }
        }
    }

    // -----------------------------------------------------------------
    // -----------------------------------------------------------------
    // C2 — the record head, and the skip it exists for
    // -----------------------------------------------------------------

    /// The contig every fixture here sits on. A block never crosses one, so a reader takes it
    /// from the block rather than from each record.
    const A_CONTIG: ContigId = ContigId(3);

    /// Lay `records` down one after another, the way a block holds them, and hand back the
    /// bytes together with the position the block starts at.
    fn a_run_of_records(records: &[SampleLocusObservations]) -> (Vec<u8>, Position) {
        let block_starts_at = records
            .first()
            .map(|first| first.region.start)
            .expect("a run holds at least one record");
        let mut encoder = RecordEncoder::for_block(block_starts_at);
        let mut bytes = Vec::new();
        for record in records {
            encoder
                .encode_record(record, &mut bytes)
                .expect("the fixtures are in coordinate order");
        }
        (bytes, block_starts_at)
    }

    /// Three records over one contig, at increasing positions, of three different shapes: the
    /// rich one, a bare covered position where every read agreed, and a tract.
    ///
    /// **The gaps between them are deliberately not the spans.** The first record covers seven
    /// bases and the next starts thirteen later, so a reader that mistook the distance to the
    /// next record's start for the record's own width fails — which is the property
    /// `positions_are_rebuilt_from_the_blocks_first_position_and_the_offsets_since` exists for,
    /// and which an earlier version of this fixture did not have, because its gap and its span
    /// were both seven.
    fn three_records_in_order() -> Vec<SampleLocusObservations> {
        let mut a_position_where_every_read_agreed = SampleLocusObservations {
            region: GenomeRegion {
                contig: A_CONTIG,
                start: Position(90_667_300),
                end: Position(90_667_300),
            },
            reference_bases: b"A".to_vec(),
            observations: Vec::new(),
            reads_without_observation: 9,
            reads_discarded_by_cap: 3,
            kind: LocusKind::Generic,
        };
        a_position_where_every_read_agreed
            .observations
            .push(SequenceObservation {
                bases: b"A".to_vec(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(1),
                num_obs: 31,
                num_fwd: 15,
                q_sum: SummedLogError::from_steps(-900),
                mapq_sum: 1_860,
                mapq_sum_sq: 111_600,
                placed_left: 12,
                chain_ids: Vec::new(),
            });

        let mut tract = a_rich_record();
        tract.region = GenomeRegion {
            contig: A_CONTIG,
            start: Position(90_670_000),
            end: Position(90_670_011),
        };
        tract.reference_bases = b"ATATATATATAT".to_vec();
        tract.kind = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide is a motif"),
            left_flank: b"GGCC".to_vec().into_boxed_slice(),
            right_flank: b"TTAA".to_vec().into_boxed_slice(),
        });

        vec![a_rich_record(), a_position_where_every_read_agreed, tract]
    }

    fn record_body_length(record: &SampleLocusObservations) -> usize {
        let mut body = Vec::new();
        encode_record_body(record, &mut body);
        body.len()
    }

    /// The whole of C2: a record laid down as a head and a body comes back as the same record,
    /// with the head's own fields intact and the byte count covering both halves.
    #[test]
    fn a_record_round_trips_through_its_head_and_its_body() {
        let mut live_reads = a_live_set_reader();
        let written = a_rich_record();
        let (bytes, block_starts_at) = a_run_of_records(std::slice::from_ref(&written));

        let decoded = decode_a_record(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect("what this encoder wrote, this decoder reads");

        assert_eq!(decoded.record, written);
        assert_eq!(decoded.head.region, written.region);
        assert_eq!(decoded.record_bytes, bytes.len());
    }

    /// **What the encoder hands back is what it wrote.** The block index and the position
    /// summary of Milestones D and F are built from this return value rather than by re-reading
    /// the bytes, so a head that disagrees with its own record would seek wrongly rather than
    /// fail. The fixture discriminates because its three records differ in all three fields.
    #[test]
    fn the_head_the_encoder_hands_back_is_the_head_it_wrote() {
        let mut live_reads = a_live_set_reader();
        for record in three_records_in_order() {
            let mut bytes = Vec::new();
            let handed = RecordEncoder::for_block(record.region.start)
                .encode_record(&record, &mut bytes)
                .expect("the fixture is at its own base");
            let read = read_a_record_head(
                &mut live_reads,
                &bytes,
                A_CONTIG,
                OffsetBase::at_block_start(record.region.start),
            )
            .expect("and reads back");
            assert_eq!(handed, read.head, "over {}", record.region);
        }
    }

    /// **The head's non-reference read count is derived, not carried by the body**, and it
    /// counts only observations that spanned the whole locus — a partial one's bases stop where
    /// its read's witness stopped, so there is nothing to compare them against.
    ///
    /// The rich fixture is the zero case: its one complete observation shows the reference's own
    /// bases, and its two partial ones are scored by neither half. So a second record is built
    /// here with two complete observations, one of them a variant.
    #[test]
    fn the_head_carries_the_reads_that_showed_something_other_than_the_reference() {
        let mut live_reads = a_live_set_reader();
        let no_reads_varied = a_rich_record();
        assert_eq!(
            no_reads_varied.non_reference_and_compared_reads(),
            (0, 137),
            "every read that could be compared showed the reference"
        );

        let mut nineteen_reads_varied = a_rich_record();
        nineteen_reads_varied
            .observations
            .push(SequenceObservation {
                bases: b"ACGTACT".to_vec(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(1),
                num_obs: 19,
                num_fwd: 8,
                q_sum: SummedLogError::from_steps(-400),
                mapq_sum: 1_140,
                mapq_sum_sq: 68_400,
                placed_left: 5,
                chain_ids: Vec::new(),
            });
        assert_eq!(
            nineteen_reads_varied.non_reference_and_compared_reads(),
            (19, 156)
        );

        for (record, expected) in [(no_reads_varied, 0u32), (nineteen_reads_varied, 19)] {
            let (bytes, block_starts_at) = a_run_of_records(std::slice::from_ref(&record));
            let found = read_a_record_head(
                &mut live_reads,
                &bytes,
                A_CONTIG,
                OffsetBase::at_block_start(block_starts_at),
            )
            .expect("the head reads");
            assert_eq!(found.head.non_reference_reads, expected);
        }
    }

    /// **A head whose non-reference count disagrees with its own body is refused.** It is
    /// derivable from the body the reader just built, and the disagreement that matters is the
    /// one that reads *low*: a varying position declaring zero is a position the cohort's first
    /// pass skips, and nothing downstream ever looks at it again.
    #[test]
    fn a_head_whose_non_reference_count_disagrees_with_its_body_is_refused() {
        let mut live_reads = a_live_set_reader();
        let mut varying = a_rich_record();
        varying.observations.push(SequenceObservation {
            bases: b"ACGTACT".to_vec(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(1),
            num_obs: 19,
            num_fwd: 8,
            q_sum: SummedLogError::from_steps(-400),
            mapq_sum: 1_140,
            mapq_sum_sq: 68_400,
            placed_left: 5,
            chain_ids: Vec::new(),
        });
        let (mut bytes, block_starts_at) = a_run_of_records(std::slice::from_ref(&varying));

        // The head is four one-byte fields; the third is the count, and 19 fits in one byte.
        assert_eq!(bytes[2], 19);
        bytes[2] = 0;

        match decode_a_record(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
            &RecordLayout::as_this_build_writes_it(),
        ) {
            Err(RecordDecodeError::Malformed { field, reason, .. }) => {
                assert_eq!(field, "non-reference-reads");
                assert!(reason.contains("19"), "got {reason}");
            }
            other => panic!("expected a refused count, got {other:?}"),
        }
    }

    /// **The head, spelled out**, for the same reason the body has one: the encoder and the
    /// decoder move together, so only a byte string says where a field went. A change here is a
    /// format change.
    #[test]
    fn the_head_this_version_writes_is_these_exact_bytes() {
        let record = SampleLocusObservations {
            region: GenomeRegion {
                contig: A_CONTIG,
                start: Position(1_040),
                end: Position(1_046),
            },
            reference_bases: b"ACGTACG".to_vec(),
            observations: Vec::new(),
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        let mut bytes = Vec::new();
        RecordEncoder::for_block(Position(1_000))
            .encode_record(&record, &mut bytes)
            .expect("a record forty bases past the base it is measured from");

        let body_bytes = record_body_length(&record);
        assert_eq!(
            bytes[..6],
            [
                40, // position-offset: 1,040 − 1,000
                7,  // reference-span: seven bases, inclusive
                0,  // non-reference-reads: no observation showed anything
                body_bytes as u8,
                0, // chain-id departures: this record names no reads, so none stopped
                0, // chain-id arrivals: and none started
            ]
        );
        assert_eq!(bytes.len(), 6 + body_bytes);
    }

    // -----------------------------------------------------------------
    // The skip
    // -----------------------------------------------------------------

    /// **A reader that does not want a record advances past it and touches no byte of its
    /// body.**
    ///
    /// ⚠ Its evidence is the comparison against `records`, which is built independently of the
    /// codec — **not** the agreement between the two walks. `decode_record` returns
    /// `read_record_head`'s head and byte count verbatim, so those two agreeing is arithmetic.
    #[test]
    fn a_walk_over_heads_alone_reaches_the_same_records_as_a_full_decode() {
        let mut live_reads = a_live_set_reader();
        let records = three_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);
        let layout = RecordLayout::as_this_build_writes_it();

        let mut at = 0usize;
        let mut measured_from = OffsetBase::at_block_start(block_starts_at);
        for (index, expected) in records.iter().enumerate() {
            let found = read_a_record_head(&mut live_reads, &bytes[at..], A_CONTIG, measured_from)
                .expect("the head reads");
            assert_eq!(found.head.region, expected.region, "record {index}");

            let decoded = decode_a_record(
                &mut live_reads,
                &bytes[at..],
                A_CONTIG,
                measured_from,
                &layout,
            )
            .expect("and the record builds");
            assert_eq!(&decoded.record, expected, "record {index}");
            assert_eq!(decoded.record_bytes, found.record_bytes, "record {index}");

            at += found.record_bytes;
            measured_from = OffsetBase::after(&found.head);
        }
        assert_eq!(at, bytes.len(), "the walk consumed the run exactly");
    }

    /// **A walk that builds only every other record reaches the same records as one that builds
    /// them all** — a skipped body strands nothing. This is Milestone C3's oracle in miniature,
    /// over one block's worth of records rather than a file's.
    #[test]
    fn a_walk_that_skips_every_other_record_matches_a_full_decode_on_the_ones_it_keeps() {
        let mut live_reads = a_live_set_reader();
        let records = three_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);
        let layout = RecordLayout::as_this_build_writes_it();

        let mut at = 0usize;
        let mut measured_from = OffsetBase::at_block_start(block_starts_at);
        let mut kept = 0usize;
        for (index, expected) in records.iter().enumerate() {
            if index % 2 == 0 {
                let decoded = decode_a_record(
                    &mut live_reads,
                    &bytes[at..],
                    A_CONTIG,
                    measured_from,
                    &layout,
                )
                .expect("a kept record builds");
                assert_eq!(&decoded.record, expected, "record {index}");
                at += decoded.record_bytes;
                measured_from = OffsetBase::after(&decoded.head);
                kept += 1;
            } else {
                let found =
                    read_a_record_head(&mut live_reads, &bytes[at..], A_CONTIG, measured_from)
                        .expect("a skipped record's head reads");
                assert_eq!(found.head.region, expected.region, "record {index}");
                at += found.record_bytes;
                measured_from = OffsetBase::after(&found.head);
            }
        }
        assert_eq!(kept, 2, "the fixture has records on both sides of the skip");
        assert_eq!(at, bytes.len());
    }

    /// A record's position is rebuilt from the block's first position and the differences since.
    /// **The fixture's gaps are not its spans** — the first record covers seven bases and the
    /// next starts thirteen later — so a reader that took the distance to the next record's
    /// start for this record's width fails here.
    #[test]
    fn positions_are_rebuilt_from_the_blocks_first_position_and_the_offsets_since() {
        let mut live_reads = a_live_set_reader();
        let records = three_records_in_order();
        assert_ne!(
            records[0].region.len(),
            records[1].region.start.get() - records[0].region.start.get(),
            "the fixture's first gap has to differ from its first span, or this proves nothing"
        );
        let (bytes, block_starts_at) = a_run_of_records(&records);

        let mut at = 0usize;
        let mut measured_from = OffsetBase::at_block_start(block_starts_at);
        for expected in &records {
            let found = read_a_record_head(&mut live_reads, &bytes[at..], A_CONTIG, measured_from)
                .expect("the head reads");
            assert_eq!(found.head.region, expected.region);
            at += found.record_bytes;
            measured_from = OffsetBase::after(&found.head);
        }
        assert_eq!(at, bytes.len());
    }

    /// **A base stale by one record is silent, and this is what it costs.** The head carries a
    /// difference and nothing else, so a reader that does not advance its base rebuilds every
    /// later record at the wrong coordinate and reports no error — the walk even consumes the
    /// run exactly, so it looks healthy. It is why the base is a type with two constructors and
    /// why the encoder holds its own.
    #[test]
    fn a_base_stale_by_one_record_moves_every_record_after_it() {
        let mut live_reads = a_live_set_reader();
        let records = three_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);

        let mut at = 0usize;
        let mut rebuilt = Vec::new();
        while at < bytes.len() {
            // Never advanced: the fault this test exists to price.
            let found = read_a_record_head(
                &mut live_reads,
                &bytes[at..],
                A_CONTIG,
                OffsetBase::at_block_start(block_starts_at),
            )
            .expect("a stale base reads without error, which is the point");
            at += found.record_bytes;
            rebuilt.push(found.head.region.start);
        }
        assert_eq!(
            at,
            bytes.len(),
            "and the walk consumes the run, so it looks healthy"
        );

        let truth: Vec<Position> = records.iter().map(|record| record.region.start).collect();
        assert_eq!(
            rebuilt[0], truth[0],
            "the first record is measured from the block"
        );
        assert_eq!(
            rebuilt[1], truth[1],
            "and so is the second, whose base is the same"
        );
        assert_ne!(
            rebuilt[2], truth[2],
            "the third is not, and nothing said so"
        );
        assert_eq!(
            truth[2].get() - rebuilt[2].get(),
            truth[1].get() - truth[0].get(),
            "it lands early by exactly the gap the base never absorbed"
        );
    }

    /// A body is handed a slice of exactly its own record, so however damaged it is it cannot
    /// read into the record after it.
    #[test]
    fn a_records_body_is_bounded_by_what_its_head_declared() {
        let mut live_reads = a_live_set_reader();
        let records = three_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);

        let found = read_a_record_head(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
        )
        .expect("the head reads");
        assert_eq!(found.body.len(), found.head.body_bytes as usize);
        assert!(
            found.record_bytes < bytes.len(),
            "the first record is not the whole run"
        );
    }

    // -----------------------------------------------------------------
    // Records the format cannot hold, and heads it cannot believe
    // -----------------------------------------------------------------

    /// A record starting before the position its offset is measured from has no representation —
    /// the head stores a distance forwards — so it is refused, and nothing is written first.
    #[test]
    fn a_record_that_starts_before_the_previous_one_is_refused() {
        let mut record = a_rich_record();
        record.region.start = Position(900);
        record.region.end = Position(900);

        let mut bytes = Vec::new();
        match RecordEncoder::for_block(Position(1_000)).encode_record(&record, &mut bytes) {
            Err(RecordEncodeError::StartsBeforeThePreviousRecord {
                previous_position,
                offered_start,
            }) => {
                assert_eq!(previous_position, Position(1_000));
                assert_eq!(offered_start, Position(900));
            }
            other => panic!("expected a refusal, got {other:?}"),
        }
        assert!(bytes.is_empty(), "and nothing was written before refusing");
    }

    /// **An offset of zero is legal and cannot be refused here.** A block's first record is
    /// measured from its own start, so zero is what that record writes. Refusing a *repeated*
    /// position belongs to the block writer of Milestone D1, which knows which record it is on.
    #[test]
    fn a_record_at_the_base_it_is_measured_from_is_written_with_a_zero_offset() {
        let record = a_rich_record();
        let mut bytes = Vec::new();
        let head = RecordEncoder::for_block(record.region.start)
            .encode_record(&record, &mut bytes)
            .expect("zero forward is not backwards");
        assert_eq!(bytes[0], 0);
        assert_eq!(head.region, record.region);
    }

    /// A region covering no reference base is a caller's mistake rather than a state the format
    /// spells, and `GenomeRegion` has public fields, so it is reachable.
    #[test]
    fn a_record_over_no_reference_base_is_refused() {
        let mut record = a_rich_record();
        record.region.start = Position(1_010);
        record.region.end = Position(1_000);

        let mut bytes = Vec::new();
        assert!(matches!(
            RecordEncoder::for_block(Position(1_000)).encode_record(&record, &mut bytes),
            Err(RecordEncodeError::EmptyRegion { .. })
        ));
        assert!(bytes.is_empty(), "and nothing was written before refusing");
    }

    /// **A record whose last base is the last coordinate there is, refused on both sides.**
    /// `GenomeRegion::len` is documented as wrong exactly there — it overflows in a debug build
    /// and reports the region empty in a release one — so a head carrying such a region hands a
    /// consumer a value that detonates when anything asks its width. The writer refuses it so
    /// that no file it writes can contain one, and the reader refuses it so that no other
    /// writer's file can either.
    #[test]
    fn a_record_at_the_coordinate_ceiling_is_refused_by_the_writer_and_by_the_reader() {
        let mut live_reads = a_live_set_reader();
        let mut record = a_rich_record();
        record.region.start = Position(u64::MAX);
        record.region.end = Position(u64::MAX);

        let mut bytes = Vec::new();
        assert!(matches!(
            RecordEncoder::for_block(Position(u64::MAX)).encode_record(&record, &mut bytes),
            Err(RecordEncodeError::EndsAtTheCoordinateCeiling { .. })
        ));
        assert!(bytes.is_empty());

        // offset 1 from u64::MAX − 1 gives a start of u64::MAX; span 1 ends there too.
        let from_a_file = vec![1u8, 1, 0, 0];
        match read_a_record_head(
            &mut live_reads,
            &from_a_file,
            A_CONTIG,
            OffsetBase::at_block_start(Position(u64::MAX - 1)),
        ) {
            Err(RecordDecodeError::Malformed { field, .. }) => {
                assert_eq!(field, "reference-span");
            }
            other => panic!("expected a refused span, got {other:?}"),
        }
    }

    /// The same on the way back: a head declaring a span of zero is damage, not a record over
    /// nothing.
    #[test]
    fn a_head_declaring_no_reference_span_is_refused() {
        let mut live_reads = a_live_set_reader();
        let bytes = vec![0u8, 0, 0, 0];
        match read_a_record_head(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(Position(1_000)),
        ) {
            Err(RecordDecodeError::Malformed { field, .. }) => {
                assert_eq!(field, "reference-span");
            }
            other => panic!("expected a refused span, got {other:?}"),
        }
    }

    /// **A head whose declared length disagrees with the body's real shape is refused**, and the
    /// offset it reports is measured from the record's first byte, not the body's.
    #[test]
    fn a_head_that_declares_more_body_than_the_body_uses_is_refused() {
        let mut live_reads = a_live_set_reader();
        let record = a_rich_record();
        let (mut bytes, block_starts_at) = a_run_of_records(std::slice::from_ref(&record));

        let body_bytes = record_body_length(&record);
        let head_bytes = bytes.len() - body_bytes;
        assert_eq!(
            usize::from(bytes[head_bytes - 3]),
            body_bytes,
            "the body's length, three bytes before the body: the two chain-id counts follow it"
        );
        bytes[head_bytes - 3] = (body_bytes + 1) as u8;
        bytes.push(0);

        // ⚠ The extra byte goes at the *end*, so the body reads one byte past what it holds. The
        // head's chain-id changes are untouched — a byte inserted before them would be read as a
        // departure count instead, which is a different fault in a different field.
        match decode_a_record(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
            &RecordLayout::as_this_build_writes_it(),
        ) {
            Err(RecordDecodeError::Malformed {
                field,
                reason,
                bytes_in,
            }) => {
                assert_eq!(field, "record-body-byte-count");
                assert!(
                    reason.contains(&(body_bytes + 1).to_string()),
                    "got {reason}"
                );
                assert!(reason.contains(&body_bytes.to_string()), "got {reason}");
                assert_eq!(
                    bytes_in,
                    head_bytes + body_bytes,
                    "the offset is where the body stopped, counted from the record's first byte"
                );
            }
            other => panic!("expected a refused body length, got {other:?}"),
        }
    }

    /// The mirror: a head under-declaring its body. The body is bounded by what the head said,
    /// so it stops inside its own slice and says so rather than reaching into the record after
    /// it — which is the property `a_records_body_is_bounded_by_what_its_head_declared` claims.
    #[test]
    fn a_head_that_declares_less_body_than_the_body_uses_is_refused() {
        let mut live_reads = a_live_set_reader();
        let records = three_records_in_order();
        let (mut bytes, block_starts_at) = a_run_of_records(&records);

        let body_bytes = record_body_length(&records[0]);
        // The first record's four scalar head fields are one byte each, so the fourth is the
        // body's length. (The chain-id changes follow it, and this record names no reads.)
        let at_the_body_length = 3;
        assert_eq!(usize::from(bytes[at_the_body_length]), body_bytes);
        bytes[at_the_body_length] = (body_bytes - 1) as u8;

        assert!(
            decode_a_record(
                &mut live_reads,
                &bytes,
                A_CONTIG,
                OffsetBase::at_block_start(block_starts_at),
                &RecordLayout::as_this_build_writes_it(),
            )
            .is_err(),
            "a body that stops short of its own slice cannot be a record"
        );
    }

    /// **A fault inside a body the head already bounded is damage, and never a short read.**
    /// `Truncated` is the class Milestone D fetches more bytes for, and more bytes cannot help
    /// here: the body is re-bounded to the same declared length however large the buffer grows.
    /// Left in that class, one damaged record is a retry that never ends.
    #[test]
    fn a_fault_inside_a_bounded_body_is_damage_and_more_bytes_do_not_help() {
        let mut live_reads = a_live_set_reader();
        let record = a_rich_record();
        let (whole, block_starts_at) = a_run_of_records(std::slice::from_ref(&record));
        let layout = RecordLayout::as_this_build_writes_it();

        // Truncate the body's own first length prefix's worth of content by claiming far more
        // reference bases than the bounded body holds.
        let mut damaged = whole.clone();
        let head_bytes = damaged.len() - record_body_length(&record);
        damaged[head_bytes] = 120;

        let exactly = decode_a_record(
            &mut live_reads,
            &damaged,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
            &layout,
        );
        let mut grown = damaged.clone();
        grown.extend_from_slice(&[0u8; 4_096]);
        let with_more = decode_a_record(
            &mut live_reads,
            &grown,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
            &layout,
        );

        match &exactly {
            Err(RecordDecodeError::Malformed { field, reason, .. }) => {
                assert_eq!(*field, "reference-bases");
                assert!(reason.contains("the head declared"), "got {reason}");
            }
            other => panic!("expected damage, got {other:?}"),
        }
        assert_eq!(
            exactly, with_more,
            "more bytes change nothing, which is why this cannot be a retry"
        );
    }

    /// **A fault a sub-parser reports is re-based to the record's first byte, and keeps its
    /// class.**
    ///
    /// The chain-id changes are parsed by their own reader, which counts from its own first
    /// byte; a caller holding the record needs the offset measured from the record's. ⚠ And the
    /// class must *not* change the way [`RecordDecodeError::inside_a_bounded_body`] changes it:
    /// that one converts `Truncated` to `Malformed` because a body is handed a slice of exactly
    /// its declared length, so running off *that* end can never be fixed by more bytes. The
    /// chain-id parser is handed the rest of the record, where more bytes is exactly what helps.
    ///
    /// Nothing tested this until a review pointed out that both mutations — dropping the
    /// re-basing, and converting the class — left every test green.
    #[test]
    fn a_sub_parsers_fault_is_re_based_and_keeps_its_class() {
        let truncated = RecordDecodeError::Truncated {
            field: "chain-id arrival",
            bytes_in: 3,
        }
        .further_in(7);
        assert!(
            matches!(
                truncated,
                RecordDecodeError::Truncated {
                    field: "chain-id arrival",
                    bytes_in: 10
                }
            ),
            "got {truncated:?}"
        );

        let malformed = RecordDecodeError::Malformed {
            field: "chain-id departure position",
            bytes_in: 2,
            reason: "position 9 of a live set holding 2 reads".to_string(),
        }
        .further_in(5);
        let RecordDecodeError::Malformed { bytes_in, .. } = malformed else {
            panic!("the class must not change: {malformed:?}")
        };
        assert_eq!(bytes_in, 7);

        let unsupported = RecordDecodeError::Unsupported {
            field: "locus-kind",
            bytes_in: 1,
            tag: 9,
        }
        .further_in(4);
        let RecordDecodeError::Unsupported { bytes_in, .. } = unsupported else {
            panic!("the class must not change: {unsupported:?}")
        };
        assert_eq!(bytes_in, 5);

        // And through a real head: a record whose chain-id stream is cut reports the fault at the
        // record's own offset, not the stream's.
        let mut record = a_rich_record();
        record.observations[0].chain_ids = vec![5, 6, 7];
        let (bytes, block_starts_at) = a_run_of_records(std::slice::from_ref(&record));
        let mut live_reads = a_live_set_reader();
        let cut_in_the_changes = 6;
        match read_a_record_head(
            &mut live_reads,
            &bytes[..cut_in_the_changes],
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
        ) {
            Err(RecordDecodeError::Truncated { field, bytes_in }) => {
                assert!(field.starts_with("chain-id "), "got {field}");
                assert!(
                    bytes_in >= 4,
                    "the offset counts from the record's first byte, not the chain-id stream's: \
                     {bytes_in}"
                );
            }
            other => panic!("a cut inside the changes is a short read, got {other:?}"),
        }
    }

    /// **A later writer's chain-id-changes field at the end of a body is measured and stepped
    /// over**, whatever it holds.
    ///
    /// ⚠ This is the one arm of `skip_unknown_field` the suite could not see: replacing its
    /// arrival-count read with a constant left all 4,770 library tests green. It is reached only
    /// by a file from a *later* writer, so nothing this build writes exercises it.
    #[test]
    fn a_later_writers_chain_id_changes_field_is_measured_and_stepped_over() {
        let record = a_rich_record();
        let mut fields = record_fields();
        fields.push(FieldSpec {
            name: FieldName("some-later-live-set".to_string()),
            encoding: FieldEncoding::ChainIdChanges,
        });
        let layout = RecordLayout::from_manifest(&a_manifest_declaring(fields))
            .expect("one unknown field at the end is not a reason to refuse the file");

        // Two departures and three arrivals, as a later writer would lay them down.
        let mut trailing = Vec::new();
        put_varint(&mut trailing, 2);
        put_varint(&mut trailing, 0);
        put_varint(&mut trailing, 4);
        put_varint(&mut trailing, 3);
        put_varint(&mut trailing, 900);
        put_varint(&mut trailing, 0);
        put_varint(&mut trailing, 12);

        let mut body = Vec::new();
        encode_record_body(&record, &mut body);
        body.extend_from_slice(&trailing);

        let decoded =
            decode_record_body(&body, record.region, &the_live_reads_of(&record), &layout)
                .expect("a field this reader does not know is measured, not refused");
        assert_eq!(decoded.record, record);
        assert_eq!(
            decoded.bytes_read,
            body.len(),
            "the whole of the later writer's field was stepped over, and nothing else"
        );
    }

    // -----------------------------------------------------------------
    // The residual: the reads one observation names are derived, not stored
    // -----------------------------------------------------------------

    /// A record whose observations name reads, with `num_obs` consistent with them.
    fn a_record_naming_reads(lists: &[&[ChainId]]) -> SampleLocusObservations {
        let mut record = a_rich_record();
        record.observations.truncate(lists.len());
        while record.observations.len() < lists.len() {
            let mut more = record.observations[0].clone();
            more.bases = vec![b'T'; record.observations.len()];
            record.observations.push(more);
        }
        for (observation, ids) in record.observations.iter_mut().zip(lists) {
            observation.chain_ids = ids.to_vec();
            observation.num_obs = ids.len() as u32;
            observation.num_fwd = 0;
            observation.placed_left = 0;
        }
        record
    }

    /// **The largest list is derived and the others are stored**, and the record reads back
    /// exactly.
    ///
    /// That is where the column's saving is: the live set in the record's head already carries
    /// the union of every observation's reads, so storing the biggest list again is the one
    /// duplication worth removing (spec `psp_chain_id_encoding.md` §4).
    #[test]
    fn the_largest_list_is_derived_and_the_record_reads_back_exactly() {
        let record = a_record_naming_reads(&[&[3, 9], &[1, 4, 40, 41, 900], &[7]]);
        assert_eq!(
            residual_observation_of(&record, &mut ResidualScratch::default()),
            1,
            "the largest list, and ties would go to the lowest index"
        );

        let (decoded, _) = round_trip(&record);
        assert_eq!(decoded.record, record);
    }

    /// **A record the derivation cannot reproduce stores every list instead.**
    ///
    /// A chain id names a read *pair*. If both mates cover this record and show different
    /// sequences, the same identifier is in two observations — and the derivation, which is the
    /// live set minus the others, would then drop it from the residual. Checkpoint E asks for
    /// chain ids that round-trip **exactly**, so the writer checks and falls back.
    #[test]
    fn a_read_named_by_two_observations_makes_the_writer_store_every_list() {
        // Identifier 4 is in both: a pair whose two mates showed different sequences here.
        let record = a_record_naming_reads(&[&[4, 8], &[4, 30, 31]]);
        assert_eq!(
            residual_observation_of(&record, &mut ResidualScratch::default()),
            record.observations.len(),
            "no observation can be derived, so none is"
        );

        let (decoded, _) = round_trip(&record);
        assert_eq!(
            decoded.record, record,
            "and the record still reads back exactly, which is the point of falling back"
        );
    }

    /// **The reads counted and never named do not leak into the derived list.**
    ///
    /// Spec `psp_chain_id_encoding.md` §5: `reads_without_observation` and
    /// `reads_discarded_by_cap` are read counts with no identifiers behind them — a read the
    /// depth cap discarded is in no observation, so if it reached the live set the residual would
    /// gain a read nobody folded. They cannot, because the live set is built from the
    /// observations' own lists and from nothing else; this is the test that says so.
    #[test]
    fn reads_counted_but_never_named_stay_out_of_the_derived_list() {
        let mut record = a_record_naming_reads(&[&[2, 3], &[10]]);
        record.reads_without_observation = 17;
        record.reads_discarded_by_cap = 41;

        let mut bytes = Vec::new();
        let mut encoder = RecordEncoder::for_block(record.region.start);
        encoder
            .encode_record_starting_a_block(&record, &mut bytes)
            .expect("a record at its block's first position");
        assert_eq!(
            encoder.live_reads().ids(),
            [2, 3, 10],
            "the live set is the observations' reads and nothing else"
        );

        let (decoded, _) = round_trip(&record);
        assert_eq!(decoded.record, record);
        assert_eq!(decoded.record.observations[0].chain_ids, [2, 3]);
        assert_eq!(decoded.record.reads_discarded_by_cap, 41, "still counted");
    }

    /// **A derived list that cannot describe its observation's reads is refused.**
    ///
    /// There is no second copy to compare a derived list against, so the guard is an inequality
    /// against a number the record already carries: an identifier names one read or two, so the
    /// list is at most `num_obs` and at least half of it. It is the same inequality the walk's
    /// own differential against production asserts.
    ///
    /// Both directions, because a derivation can go wrong either way: too many when the
    /// subtraction is lost, too few when the live set is short.
    #[test]
    fn a_derived_list_that_cannot_describe_its_reads_is_refused() {
        let record = a_record_naming_reads(&[&[2, 3, 4], &[10]]);
        let mut bytes = Vec::new();
        encode_record_body(&record, &mut bytes);
        let live = the_live_reads_of(&record);

        // As written it reads: three derived reads for an observation showing three.
        let _ = decode_record_body(
            &bytes,
            record.region,
            &live,
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect("the record this writer wrote");

        // **The writer never produces either shape**, because it refuses to derive a list its own
        // reader would refuse — so the check is put to the decoder by damaging one byte of a
        // record the writer *did* derive: the read count of the derived observation.
        //
        // reference-bases(1 + 7) + observation-count(1) + residual-observation(1)
        // + residual-read-count(1) + observation-bases(1 + 7) + witness-run-count(1)
        // + read-group(1) = 21.
        const AT_THE_DERIVED_OBSERVATIONS_READ_COUNT: usize = 21;
        assert_eq!(
            bytes[AT_THE_DERIVED_OBSERVATIONS_READ_COUNT], 3,
            "the byte damaged below has to be the read count this fixture wrote"
        );
        for (what, num_obs) in [("too many", 2u8), ("too few", 7u8)] {
            let mut spliced = bytes.clone();
            spliced[AT_THE_DERIVED_OBSERVATIONS_READ_COUNT] = num_obs;
            let refused = decode_record_body(
                &spliced,
                record.region,
                &live,
                &RecordLayout::as_this_build_writes_it(),
            )
            .expect_err("a derived list that cannot describe its reads");
            let RecordDecodeError::Malformed { field, reason, .. } = refused else {
                panic!("{what}: got {refused:?}")
            };
            assert_eq!(field, "observation-reads", "{what}");
            assert!(reason.contains("derived list"), "{what}: {reason}");
        }
    }

    /// **A list that is not a set is written as one, rather than as different reads.**
    ///
    /// An observation's reads go on the wire as ascending gaps, so a repeated or out-of-order
    /// identifier has no honest spelling. ⚠ Measured on the version that trusted its caller: an
    /// observation handed `[3, 3]` wrote bytes that read back as the reads `[3, 4]` — gaining
    /// read 4, which nothing folded, with neither side reporting anything. That is spec §5's
    /// failure reached from the writer.
    ///
    /// `SequenceObservation::chain_ids` is a `pub` field documenting what the ids *mean* and not
    /// what order they are in, and `encode_record_body` is a `pub` function. Both of ng's pileup
    /// paths do sort and deduplicate; the codec no longer depends on it.
    #[test]
    fn a_list_that_is_not_a_set_is_written_as_the_set_it_names() {
        for (what, given, wanted) in [
            ("a repeat", vec![3u64, 3], vec![3u64]),
            ("out of order", vec![7u64, 3], vec![3u64, 7]),
            ("both", vec![5u64, 5, 1, 5], vec![1u64, 5]),
        ] {
            let mut record = a_record_naming_reads(&[&[], &[]]);
            record.observations[0].chain_ids = given.clone();
            record.observations[0].num_obs = given.len() as u32;
            // The other observation names more, so the first one's list is the stored one.
            record.observations[1].chain_ids = vec![100, 101, 102, 103];
            record.observations[1].num_obs = 4;

            let mut bytes = Vec::new();
            encode_record_body(&record, &mut bytes);
            let live = LiveSet::from_sorted_ids({
                let mut every = wanted.clone();
                every.extend_from_slice(&record.observations[1].chain_ids);
                as_a_read_set(&mut every);
                every
            });
            let decoded = decode_record_body(
                &bytes,
                record.region,
                &live,
                &RecordLayout::as_this_build_writes_it(),
            )
            .unwrap_or_else(|refused| panic!("{what}: {refused}"));
            assert_eq!(
                decoded.record.observations[0].chain_ids, wanted,
                "{what}: {given:?} names the reads {wanted:?} and no others"
            );
        }
    }

    /// **A later writer's chain-id *list* is measured as a list**, which is one counted run where
    /// the changes are two.
    ///
    /// ⚠ The two shared a scheme when this field first landed, and nothing noticed: stepping over
    /// a list under the changes' rule reads its identifiers as a departure count and then asks
    /// for an arrival count that is not there. `[5, 9]` is three bytes as a list and seven as a
    /// set of changes.
    #[test]
    fn a_later_writers_chain_id_list_is_measured_as_a_list() {
        let record = a_rich_record();
        let mut fields = record_fields();
        fields.push(FieldSpec {
            name: FieldName("some-later-read-list".to_string()),
            encoding: FieldEncoding::ChainIdList,
        });
        let layout = RecordLayout::from_manifest(&a_manifest_declaring(fields))
            .expect("one unknown field at the end is not a reason to refuse the file");

        let mut trailing = Vec::new();
        encode_read_list(&[5, 9], &mut trailing);
        assert_eq!(trailing.len(), 3, "a count and two gaps");

        let mut body = Vec::new();
        encode_record_body(&record, &mut body);
        body.extend_from_slice(&trailing);

        let decoded =
            decode_record_body(&body, record.region, &the_live_reads_of(&record), &layout)
                .expect("a field this reader does not know is measured, not refused");
        assert_eq!(decoded.record, record);
        assert_eq!(
            decoded.bytes_read,
            body.len(),
            "the whole of the later writer's list was stepped over, and nothing else"
        );
    }

    /// **A residual index past the record's own observations is damage.**
    ///
    /// ⚠ The bound had no test: deleting it left all 241 `ng::psp` tests green while a body
    /// claiming observation 200 of a one-observation record flipped from refused to accepted.
    /// The index is a `usize` the reader would then compare against a `Vec`'s length — which is
    /// safe by luck rather than by check, and reports nothing.
    #[test]
    fn a_residual_index_past_the_records_observations_is_damage() {
        // An empty reference, one observation, residual index 200.
        let body = vec![0u8, 1, 200, 1];
        let refused = decode_record_body(
            &body,
            a_region(1, 1),
            &LiveSet::new(),
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect_err("observation 200 of a record holding one");
        let RecordDecodeError::Malformed { field, reason, .. } = refused else {
            panic!("got {refused:?}")
        };
        assert_eq!(field, "residual-observation");
        assert!(reason.contains("200"), "{reason}");

        // And the count itself is the sentinel for "none", which must still be accepted. Built
        // through the encoder rather than by hand, so the tail of the body is whatever this
        // build writes.
        let names_nothing = a_record_naming_reads(&[&[]]);
        assert_eq!(
            residual_observation_of(&names_nothing, &mut ResidualScratch::default()),
            names_nothing.observations.len(),
            "an observation naming no reads is not worth deriving, so nothing is"
        );
        let mut none = Vec::new();
        encode_record_body(&names_nothing, &mut none);
        let decoded = decode_record_body(
            &none,
            names_nothing.region,
            &LiveSet::new(),
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect("the observation count means no observation is derived");
        assert_eq!(decoded.record, names_nothing);
    }

    /// **A live set carrying reads no observation named is refused**, which an inequality alone
    /// could not do.
    ///
    /// Spec §5 proposes bounding the derived count by the observation's read count. That window
    /// is `num_obs / 2` wide, and its slack is exactly the number of read pairs whose two mates
    /// both cover this record — the shape paired-end data has. Measured by the review: a residual
    /// naming two reads with `num_obs = 4`, against a live set carrying two identifiers no
    /// observation named, derives a list of four, passes `4 <= 4` and `8 >= 4`, and the reference
    /// allele silently gains two reads.
    ///
    /// The record carries the residual's *length* now, so the check is an equality and that
    /// record is refused.
    #[test]
    fn a_live_set_carrying_reads_no_observation_named_is_refused() {
        let mut record = a_record_naming_reads(&[&[10, 11], &[20]]);
        // Both mates of each pair cover here, so the read count is twice the identifier count —
        // the regime where the inequality has slack.
        record.observations[0].num_obs = 4;
        record.observations[1].num_obs = 2;

        let mut bytes = Vec::new();
        encode_record_body(&record, &mut bytes);

        // What the writer meant: reads 10, 11 and 20 are live.
        let honest = LiveSet::from_sorted_ids(vec![10, 11, 20]);
        let decoded = decode_record_body(
            &bytes,
            record.region,
            &honest,
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect("the record this writer wrote");
        assert_eq!(decoded.record.observations[0].chain_ids, [10, 11]);

        // And a live set with two reads nobody named: four derived where the record says two.
        let phantom = LiveSet::from_sorted_ids(vec![10, 11, 12, 13, 20]);
        let refused = decode_record_body(
            &bytes,
            record.region,
            &phantom,
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect_err("two reads no observation named");
        let RecordDecodeError::Malformed { field, reason, .. } = refused else {
            panic!("got {refused:?}")
        };
        assert_eq!(field, "observation-reads");
        assert!(
            reason.contains("where the record says"),
            "the equality is what catches this, not the inequality: {reason}"
        );
    }

    /// **A tie for the largest list goes to the lower index**, which is what makes the choice a
    /// function of the record alone — and so the file the same at any worker count.
    ///
    /// ⚠ Reversing the tie-break left all 4,778 library tests green.
    #[test]
    fn a_tie_for_the_largest_list_goes_to_the_lower_index() {
        let record = a_record_naming_reads(&[&[3, 4], &[10, 11], &[20]]);
        assert_eq!(
            residual_observation_of(&record, &mut ResidualScratch::default()),
            0
        );

        let the_other_way = a_record_naming_reads(&[&[20], &[10, 11], &[3, 4]]);
        assert_eq!(
            residual_observation_of(&the_other_way, &mut ResidualScratch::default()),
            1
        );
    }

    /// **The reads a record names are the union over all of its observations**, not the first
    /// one's.
    ///
    /// A record's observations are split by allele, by witness and by read group, so a locus with
    /// two alleles from two lanes is four observations and the reads are spread across them. ⚠
    /// Measured: a writer that named only the first observation's reads passed all 4,770 library
    /// tests, because no fixture had ever put ids on two observations.
    #[test]
    fn the_reads_a_record_names_are_the_union_over_its_observations() {
        let mut record = a_rich_record();
        record.observations[0].chain_ids = vec![3, 7];
        record.observations[1].chain_ids = vec![3, 900];
        record.observations[2].chain_ids = vec![41];

        let mut bytes = Vec::new();
        let mut encoder = RecordEncoder::for_block(record.region.start);
        encoder
            .encode_record_starting_a_block(&record, &mut bytes)
            .expect("a record at its own block's first position");
        assert_eq!(
            encoder.live_reads().ids(),
            [3, 7, 41, 900],
            "sorted, deduplicated, and from every observation"
        );

        let mut live_reads = a_live_set_reader();
        let _ = read_a_record_head(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(record.region.start),
        )
        .expect("the head reads");
        assert_eq!(live_reads.live().ids(), [3, 7, 41, 900]);
    }

    /// A whole record cut anywhere is `Truncated` — including the cut where the head is complete
    /// and the body is one byte short, which is what a record straddling Milestone D's rolling
    /// buffer looks like.
    ///
    /// **And the live set does not move at any of those cuts**, which is the half this test was
    /// blind to. ⚠ Its record had no chain ids, so every cut met the two bytes `0, 0` and a set
    /// moved by them is a set unchanged — the fixture names reads now. Measured on the version
    /// that applied the changes before bounding the body: the retry after a cut in the body met
    /// an arrival already live and refused a perfectly good record.
    #[test]
    fn a_record_cut_short_is_truncated_at_every_cut() {
        let mut live_reads = a_live_set_reader();
        let mut record = a_rich_record();
        record.observations[0].chain_ids = vec![10, 11, 12];
        record.observations[1].chain_ids = vec![12, 40];
        let (bytes, block_starts_at) = a_run_of_records(std::slice::from_ref(&record));
        let layout = RecordLayout::as_this_build_writes_it();
        let head_bytes = bytes.len() - record_body_length(&record);

        for cut in 0..bytes.len() {
            match read_a_record_head(
                &mut live_reads,
                &bytes[..cut],
                A_CONTIG,
                OffsetBase::at_block_start(block_starts_at),
            ) {
                Err(RecordDecodeError::Truncated { field, .. }) if cut >= head_bytes => {
                    // The head is whole and the body is short, so it is the body's own bounding
                    // that ran out — the one field a cut here may name.
                    assert_eq!(
                        field, "record-body-byte-count",
                        "{cut} bytes stops in the body"
                    );
                }
                Err(RecordDecodeError::Truncated { field, .. }) => {
                    // Inside the head, the field named is whichever one the bytes ran out on.
                    // The chain-id changes are a head field whose *parts* have their own names,
                    // because a fault in one of them has to say which run it was in.
                    assert!(
                        RECORD_HEAD_FIELDS.iter().any(|(name, _)| *name == field)
                            || field.starts_with("chain-id "),
                        "{cut} bytes stops in the head, and {field} is not one of its fields"
                    );
                }
                other => panic!("{cut} bytes of a record must be Truncated, got {other:?}"),
            }
            assert!(
                live_reads.live().is_empty(),
                "a record cut at {cut} moved the live set before it failed, so the retry the \
                 fault instructs would meet its own changes a second time"
            );
            assert!(
                matches!(
                    decode_a_record(
                        &mut live_reads,
                        &bytes[..cut],
                        A_CONTIG,
                        OffsetBase::at_block_start(block_starts_at),
                        &layout
                    ),
                    Err(RecordDecodeError::Truncated { .. })
                ),
                "{cut} bytes through decode_record must be Truncated too"
            );
        }
        assert!(
            read_a_record_head(
                &mut live_reads,
                &bytes,
                A_CONTIG,
                OffsetBase::at_block_start(block_starts_at)
            )
            .is_ok()
        );
    }

    /// A position offset that would run off the coordinate axis is damage rather than a wrapped
    /// coordinate, and so is a span that runs off it from a legal start.
    #[test]
    fn a_head_that_runs_off_the_coordinate_axis_is_refused() {
        let mut live_reads = a_live_set_reader();
        let mut runaway_offset = Vec::new();
        encode_u64_leb128(u64::MAX, &mut runaway_offset);
        runaway_offset.extend_from_slice(&[1, 0, 0]);
        match read_a_record_head(
            &mut live_reads,
            &runaway_offset,
            A_CONTIG,
            OffsetBase::at_block_start(Position(1_000)),
        ) {
            Err(RecordDecodeError::Malformed { field, .. }) => {
                assert_eq!(field, "position-offset");
            }
            other => panic!("expected a refused offset, got {other:?}"),
        }

        let mut runaway_span = vec![0u8];
        encode_u64_leb128(u64::MAX, &mut runaway_span);
        runaway_span.extend_from_slice(&[0, 0]);
        match read_a_record_head(
            &mut live_reads,
            &runaway_span,
            A_CONTIG,
            OffsetBase::at_block_start(Position(1_000)),
        ) {
            Err(RecordDecodeError::Malformed { field, .. }) => {
                assert_eq!(field, "reference-span");
            }
            other => panic!("expected a refused span, got {other:?}"),
        }
    }

    /// A head count wider than the `u32` it is read into is damage, not a narrowed number — the
    /// same rule the body's counts keep. A truncated body length does not fail: it re-bases
    /// every record after it in the block, plausibly, which is the silent failure the head
    /// exists to make impossible.
    #[test]
    fn a_head_count_too_large_for_its_field_is_refused_rather_than_narrowed() {
        let mut live_reads = a_live_set_reader();
        for (position, expected) in [
            (2usize, "non-reference-reads"),
            (3, "record-body-byte-count"),
        ] {
            let mut bytes = vec![0u8, 1]; // offset 0, span 1
            for at in 2..=3 {
                if at == position {
                    encode_u64_leb128(u64::from(u32::MAX) + 6, &mut bytes);
                } else {
                    bytes.push(0);
                }
            }
            bytes.extend_from_slice(&[0u8; 8]);
            match read_a_record_head(
                &mut live_reads,
                &bytes,
                A_CONTIG,
                OffsetBase::at_block_start(Position(1_000)),
            ) {
                Err(RecordDecodeError::Malformed { field, reason, .. }) => {
                    assert_eq!(field, expected);
                    assert!(reason.contains("4294967301"), "got {reason}");
                }
                other => panic!("{expected} past a u32: expected Malformed, got {other:?}"),
            }
        }
    }

    /// **A record body is at most 4,294,967,295 bytes, and the head's field, the decoder's
    /// ceiling and the encoder's guard all have to say so.** Each is a separate spelling of one
    /// rule, and a review moved all three one at a time with the whole suite green.
    #[test]
    fn a_body_is_at_most_four_gibibytes_wherever_that_is_written() {
        assert_eq!(MOST_BYTES_A_BODY_CAN_DECLARE, 4_294_967_295);
        assert_eq!(std::mem::size_of::<BodyByteCount>(), 4);
        assert_eq!(
            declared_body_bytes(BodyByteCount::MAX as usize)
                .expect("the largest body a head can describe"),
            BodyByteCount::MAX
        );
        assert!(matches!(
            declared_body_bytes(BodyByteCount::MAX as usize + 1),
            Err(RecordEncodeError::BodyTooLong { body_bytes }) if body_bytes == BodyByteCount::MAX as usize + 1
        ));
    }

    /// **A head declaring a body larger than any buffer holds is still `Truncated`,
    /// deliberately.** Spec §8 refuses a fixed maximum record size — "a single record can exceed
    /// the rolling buffer … the buffer has to grow rather than fail" — so there is no length at
    /// which this reader may call a declared body damage. What bounds it is the head's own
    /// `u32`, which the test above pins, and Milestone D's block.
    #[test]
    fn a_head_declaring_a_body_no_buffer_holds_is_truncated_not_malformed() {
        let mut live_reads = a_live_set_reader();
        let mut bytes = vec![0u8, 1, 0];
        encode_u64_leb128(u64::from(BodyByteCount::MAX), &mut bytes);
        // No departures and no arrivals, so the head is whole and what runs out is the body.
        bytes.extend_from_slice(&[0, 0]);
        match read_a_record_head(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(Position(1_000)),
        ) {
            Err(RecordDecodeError::Truncated { field, .. }) => {
                assert_eq!(field, "record-body-byte-count");
            }
            other => panic!("expected Truncated, got {other:?}"),
        }
    }

    /// **The top of the depth range, where a head's fields stop fitting in one byte.** Every
    /// other fixture here has a single-byte body length and a single-byte read count; at three
    /// hundred reads a position a record carries hundreds of observations and both need several.
    /// That is the regime the hand-written head/body offsets elsewhere in these tests assume
    /// away, and the caller is committed to it.
    #[test]
    fn a_record_at_the_depth_the_caller_must_survive_round_trips_through_its_head() {
        let mut live_reads = a_live_set_reader();
        let deep = a_record_at_three_hundred_reads();
        let (bytes, block_starts_at) = a_run_of_records(std::slice::from_ref(&deep));

        let found = read_a_record_head(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
        )
        .expect("the head reads");
        assert!(
            found.head.body_bytes > 255,
            "the body length needs more than one byte; it is {}",
            found.head.body_bytes
        );
        assert!(
            found.head.non_reference_reads > 255,
            "and so does the read count; it is {}",
            found.head.non_reference_reads
        );

        let decoded = decode_a_record(
            &mut live_reads,
            &bytes,
            A_CONTIG,
            OffsetBase::at_block_start(block_starts_at),
            &RecordLayout::as_this_build_writes_it(),
        )
        .expect("and the record builds");
        assert_eq!(decoded.record, deep);
        assert_eq!(decoded.record_bytes, bytes.len());
        assert_eq!(decoded.head, found.head);
    }

    proptest! {
        /// Any head this encoder writes reads back as the head it wrote, from any base — the
        /// offset is the one running quantity in the head and the only thing a wrong base moves.
        #[test]
        fn any_head_round_trips_from_any_base(
            base in 0u64..1_000_000_000,
            forward in 0u64..1_000_000,
            span in 1u64..500,
        ) {
            let mut live_reads = a_live_set_reader();
            let record = SampleLocusObservations {
                region: GenomeRegion {
                    contig: A_CONTIG,
                    start: Position(base + forward),
                    end: Position(base + forward + span - 1),
                },
                reference_bases: vec![b'A'; span as usize],
                observations: Vec::new(),
                reads_without_observation: 1,
                reads_discarded_by_cap: 0,
                kind: LocusKind::Generic,
            };
            let mut bytes = Vec::new();
            RecordEncoder::for_block(Position(base))
                .encode_record(&record, &mut bytes)
                .expect("a record at or past its base");
            let found = read_a_record_head(&mut live_reads, &bytes, A_CONTIG, OffsetBase::at_block_start(Position(base)))
                .expect("what this encoder wrote, this decoder reads");
            prop_assert_eq!(found.head.region, record.region);
            prop_assert_eq!(found.record_bytes, bytes.len());
        }

        /// **Arbitrary bytes are read or refused, never a panic, and never a region a consumer
        /// cannot ask about.** These two are the module's outermost untrusted-input entry
        /// points, and `region.len()` is called deliberately: a head this reader accepts is a
        /// head a block reader will measure, and the coordinate ceiling is where measuring one
        /// used to overflow.
        #[test]
        fn arbitrary_bytes_read_as_a_record_or_are_refused_but_never_panic(
            bytes in prop::collection::vec(any::<u8>(), 0..256),
            base in prop::sample::select(vec![0u64, 1, 1_000, u64::MAX / 2, u64::MAX - 1, u64::MAX]),
        ) {
            let mut live_reads = a_live_set_reader();
            let measured_from = OffsetBase::at_block_start(Position(base));
            if let Ok(found) = read_a_record_head(&mut live_reads, &bytes, A_CONTIG, measured_from) {
                prop_assert_eq!(found.body.len(), found.head.body_bytes as usize);
                prop_assert!(found.record_bytes <= bytes.len());
                // `len()` on purpose, not `is_empty()`: it is the accessor that used
                // to overflow at the coordinate ceiling, so calling it is the assertion.
                let width = found.head.region.len();
                prop_assert!(width >= 1);
            }
            if let Ok(decoded) = decode_a_record(
            &mut live_reads,
                &bytes,
                A_CONTIG,
                measured_from,
                &RecordLayout::as_this_build_writes_it(),
            ) {
                prop_assert!(decoded.record_bytes <= bytes.len());
                let width = decoded.head.region.len();
                prop_assert!(width >= 1);
            }
        }
    }

    // -----------------------------------------------------------------
    // C3 — the body stands on its own
    // -----------------------------------------------------------------

    /// How many records [`twelve_records_in_order`] holds.
    ///
    /// **Twelve, and the number is load-bearing three times over**: the span list has that many
    /// entries, the "only the last" skip pattern names the last by ordinal, and the exhaustive
    /// mask test walks 2^12 patterns. Three records have no *inside* — every pattern over them is
    /// the first, the last or the middle — and twelve gives patterns that skip runs.
    const RECORDS_IN_THE_RUN: usize = 12;

    // The exhaustive skip test indexes a mask by record, so the run has to fit one.
    const _: () = assert!(RECORDS_IN_THE_RUN < 32);

    /// A record at the top of the depth range — **three hundred reads a position**, which is
    /// where a head's body length and its read count stop fitting in one byte each.
    ///
    /// Shared by the head test that asserts those two widths and by the run below, so the depth
    /// this caller is committed to is written down once.
    fn a_record_at_three_hundred_reads() -> SampleLocusObservations {
        let mut deep = a_rich_record();
        deep.observations.clear();
        for read in 0..300u32 {
            deep.observations.push(SequenceObservation {
                // Not the reference's own bases, so every one of these reads is non-reference.
                bases: b"ACGTACT".to_vec(),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(read % 7),
                num_obs: 1,
                num_fwd: read % 2,
                q_sum: SummedLogError::from_steps(-(i64::from(read) + 1)),
                mapq_sum: 60,
                mapq_sum_sq: 3_600,
                placed_left: 1,
                chain_ids: Vec::new(),
            });
        }
        deep
    }

    /// Twelve records over one contig in coordinate order, of every shape this codec can be
    /// handed.
    ///
    /// **Not called a block**, though a later milestone's blocks will hold runs like it: a psp
    /// block is a span of reference compressed as one unit and cut on a 100 kb grid
    /// (spec §2, §4.1), and this is 232 bases of records in memory. Blocks arrive at D1.
    ///
    /// What varies, and every clause of this is asserted by
    /// `the_run_is_the_run_its_doc_describes` rather than left in prose:
    ///
    /// - **spans of one to twelve bases**, and **gaps that are never the record's own span**, so
    ///   a reader that took the distance to the next record for this record's width fails;
    /// - all three locus kinds, with **a different motif and different flanks on every tract
    ///   record**, so a payload taken from the previous tract record rather than from this one's
    ///   own bytes fails;
    /// - **reference bases that differ between records of the same length**, for the same reason;
    /// - records with no observations at all, records where every read agreed with the reference,
    ///   and three at three hundred reads a position;
    /// - **chain ids on some observations.** C1's encoder drops them, so the comparisons go
    ///   through [`as_a_decode_can_return_it`] — and the day Milestone E writes them, that helper
    ///   becomes the identity and this oracle covers the exception lists with no test touched.
    fn twelve_records_in_order() -> Vec<SampleLocusObservations> {
        let spans = [1u64, 7, 1, 3, 1, 12, 1, 1, 4, 9, 1, 2];
        let mut records = Vec::new();
        // The block's first position. Inert: every offset in the file is relative, so starting
        // at zero would give byte-identical output.
        let mut at = 1_000u64;

        for (index, span) in spans.into_iter().enumerate() {
            let mut record = match index % 4 {
                0 => {
                    // A covered position where every read agreed with the reference — most of a
                    // real file. Its one observation's bases are set below, once the reference
                    // this record carries is known.
                    let mut every_read_agreed = a_rich_record();
                    every_read_agreed.observations.truncate(1);
                    every_read_agreed.observations[0].read_witness = ReadWitness::Complete;
                    every_read_agreed.reads_without_observation = 4 + index as u32;
                    every_read_agreed
                }
                1 => {
                    // The rich fixture, whose partial witnesses name locus offsets up to seven —
                    // so the spans on this arm are seven bases or more, or the record could not
                    // occur.
                    assert!(
                        span >= 7,
                        "record {index} inherits a witness reaching offset 7 and has {span} bases"
                    );
                    a_rich_record()
                }
                2 => {
                    let mut no_read_showed_anything = a_rich_record();
                    no_read_showed_anything.observations.clear();
                    no_read_showed_anything.reads_without_observation = 11 + index as u32;
                    no_read_showed_anything.reads_discarded_by_cap = 1 + index as u32;
                    no_read_showed_anything
                }
                _ => a_record_at_three_hundred_reads(),
            };

            record.region = GenomeRegion {
                contig: A_CONTIG,
                start: Position(at),
                end: Position(at + span - 1),
            };
            // Different bases in records of the same length, so a reference taken from another
            // record is visible here and not only in a property test that happens to differ.
            record.reference_bases = (0..span)
                .map(|offset| b"ACGT"[((index as u64 + offset) % 4) as usize])
                .collect::<Vec<_>>();
            if index % 4 == 0 {
                record.observations[0].bases = record.reference_bases.clone();
            }
            record.kind = match index % 3 {
                0 => LocusKind::Generic,
                1 => LocusKind::Ssr(SsrDetail {
                    // A different motif and different flanks on every tract record.
                    motif: Motif::new(match index {
                        1 => &b"AT"[..],
                        4 => &b"CG"[..],
                        7 => &b"AAT"[..],
                        _ => &b"GCCT"[..],
                    })
                    .expect("each is a motif of one to six bases"),
                    left_flank: vec![b'G'; index].into_boxed_slice(),
                    right_flank: vec![b'T'; 1 + index % 3].into_boxed_slice(),
                }),
                _ => LocusKind::SsrBundle,
            };
            // Chain ids on some observations, distinct per record so ids taken from another
            // record's list would be visible.
            //
            // ⚠ **`num_obs` moves with the list, and that is what makes these records exercise
            // the residual at all.** The writer refuses to derive a list whose length cannot
            // describe its observation's read count, so while these observations kept
            // `a_rich_record`'s 137 reads against two identifiers, **not one of the twelve
            // records ever derived anything** — every one fell back to storing every list, and
            // the whole record-level oracle below ran past the step it was extended for.
            // Measured by the review that found it: a probe that panics on the derive path fired
            // in 8 tests of 241 before this line, and 16 after.
            for (which, observation) in record.observations.iter_mut().enumerate().take(2) {
                if index % 3 == 1 {
                    observation.chain_ids = vec![
                        (index * 100 + which) as u64,
                        (index * 100 + which + 7) as u64,
                    ];
                    observation.num_obs = observation.chain_ids.len() as u32;
                    observation.num_fwd = observation.num_fwd.min(observation.num_obs);
                }
            }

            // The next record starts past this one's end, by something that is not its span.
            at += span + 13 + index as u64;
            records.push(record);
        }
        assert_eq!(records.len(), RECORDS_IN_THE_RUN);
        records
    }

    /// The record a decode can currently return for `record`: everything it carries except the
    /// chain ids, which C1's encoder drops
    /// (`chain_ids_come_back_empty_and_nothing_else_changes`).
    ///
    /// **When Milestone E starts writing them this becomes the identity and can be deleted**, and
    /// the oracle covers the exception lists from that moment without a test being touched.
    fn as_a_decode_can_return_it(record: &SampleLocusObservations) -> SampleLocusObservations {
        record.clone()
    }

    /// One record as a walk met it: what its head said, how far it advanced, and the record
    /// itself where the walk was asked to build one.
    #[derive(Debug, Clone, PartialEq)]
    struct WalkedRecord {
        head: RecordHead,
        record_bytes: usize,
        built: Option<SampleLocusObservations>,
    }

    /// Walk `bytes`, building the records `is_kept` says to build and reading only the head of
    /// the rest. **Every head is read** — a skipping reader must, because the head carries the
    /// running position and, from Milestone E3, the chain-id changes.
    ///
    /// `what` names the walk, so a refusal from inside says which of a test's several patterns
    /// was running.
    ///
    /// **It asserts that it built exactly what it was asked to.** Without that, a harness that
    /// quietly built nothing would leave every test below green and proving nothing — which is
    /// what a review measured before this assertion existed.
    fn walk_records(
        what: &str,
        bytes: &[u8],
        block_starts_at: Position,
        layout: &RecordLayout,
        is_kept: impl Fn(usize) -> bool,
    ) -> Vec<WalkedRecord> {
        let mut live_reads = a_live_set_reader();
        let mut walked = Vec::new();
        let mut at = 0usize;
        let mut measured_from = OffsetBase::at_block_start(block_starts_at);
        while at < bytes.len() {
            let index = walked.len();
            let found = read_a_record_head(&mut live_reads, &bytes[at..], A_CONTIG, measured_from)
                .unwrap_or_else(|refused| panic!("{what}: record {index}'s head reads: {refused}"));
            let built = if is_kept(index) {
                // **From the head just read, not from the bytes again.** Reading a head applies
                // that record's chain-id changes, so parsing it twice applies them twice — and
                // the second time an arriving read is already live, which is damage.
                let decoded = decode_the_body_of(&found, live_reads.live(), layout)
                    .unwrap_or_else(|refused| panic!("{what}: record {index} builds: {refused}"));
                assert_eq!(
                    decoded.head, found.head,
                    "{what}: record {index}'s head read the same both ways"
                );
                assert_eq!(
                    decoded.record_bytes, found.record_bytes,
                    "{what}: record {index}"
                );
                Some(decoded.record)
            } else {
                None
            };
            at += found.record_bytes;
            measured_from = OffsetBase::after(&found.head);
            walked.push(WalkedRecord {
                head: found.head,
                record_bytes: found.record_bytes,
                built,
            });
        }
        assert_eq!(at, bytes.len(), "{what}: the walk consumed the run exactly");
        for (index, met) in walked.iter().enumerate() {
            assert_eq!(
                met.built.is_some(),
                is_kept(index),
                "{what}: record {index} was built exactly when the walk was asked to build it"
            );
        }
        walked
    }

    /// Every record the walk met is where it was written, and every one it built is exactly the
    /// record it was written from.
    ///
    /// **The evidence is the comparison against `records`, which is built independently of the
    /// codec** — not the agreement between two walks, which would be arithmetic.
    fn walk_matches(what: &str, walked: &[WalkedRecord], records: &[SampleLocusObservations]) {
        assert_eq!(walked.len(), records.len(), "{what}: every record was met");
        for (index, met) in walked.iter().enumerate() {
            assert_eq!(
                met.head.region, records[index].region,
                "{what}: record {index} landed where it was written"
            );
            if let Some(built) = &met.built {
                assert_eq!(
                    built,
                    &as_a_decode_can_return_it(&records[index]),
                    "{what}: record {index}"
                );
            }
        }
    }

    /// One skip pattern: what it is called, and which records it builds.
    type SkipPattern = (&'static str, fn(usize) -> bool);

    /// One way of damaging a skipped body: what it is called, and what to fill it with.
    type BodyDamage = (&'static str, fn(usize) -> Vec<u8>);

    /// The six skip patterns the hand-written tests below run: the edges, and the two halves.
    fn the_skip_patterns() -> [SkipPattern; 6] {
        [
            ("every record", |_| true),
            ("every even record", |index| index.is_multiple_of(2)),
            ("every odd record", |index| !index.is_multiple_of(2)),
            ("one record in four", |index| index.is_multiple_of(4)),
            ("only the last", |index| index + 1 == RECORDS_IN_THE_RUN),
            ("none at all", |_| false),
        ]
    }

    /// **What `twelve_records_in_order` claims about itself.** Every number in that fixture is
    /// chosen for a reason its doc gives, and this is where the reasons are checked — so that an
    /// edit made for something else cannot quietly narrow what the oracle covers. A review
    /// measured the previous fixture collapsing to twelve identical one-base records with the
    /// whole suite still green.
    #[test]
    fn the_run_is_the_run_its_doc_describes() {
        let records = twelve_records_in_order();
        assert_eq!(records.len(), RECORDS_IN_THE_RUN);

        let spans: Vec<u64> = records.iter().map(|record| record.region.len()).collect();
        assert_eq!(spans.iter().copied().min(), Some(1), "a one-base record");
        assert_eq!(spans.iter().copied().max(), Some(12), "and a widened one");
        for pair in records.windows(2) {
            let gap = pair[1].region.start.get() - pair[0].region.start.get();
            assert_ne!(
                gap,
                pair[0].region.len(),
                "a gap equal to its own record's span proves nothing about rebuilding positions"
            );
        }

        assert!(records.iter().any(|record| record.observations.is_empty()));
        assert!(
            records
                .iter()
                .any(|record| matches!(record.kind, LocusKind::Generic))
        );
        assert!(
            records
                .iter()
                .any(|record| matches!(record.kind, LocusKind::SsrBundle))
        );
        let tracts: Vec<&SsrDetail> = records
            .iter()
            .filter_map(|record| match &record.kind {
                LocusKind::Ssr(detail) => Some(detail),
                _ => None,
            })
            .collect();
        assert!(tracts.len() >= 2, "several tract records");
        for pair in tracts.windows(2) {
            assert_ne!(
                pair[0], pair[1],
                "two tract records with the same payload cannot show one taken from the other"
            );
        }
        assert!(
            records
                .iter()
                .any(|record| record.observations.iter().any(|o| !o.chain_ids.is_empty())),
            "chain ids, so Milestone E is covered the day it writes them"
        );

        // Two records of the same length whose reference bases differ.
        let mut by_length: std::collections::BTreeMap<usize, Vec<&[u8]>> =
            std::collections::BTreeMap::new();
        for record in &records {
            by_length
                .entry(record.reference_bases.len())
                .or_default()
                .push(&record.reference_bases);
        }
        assert!(
            by_length
                .values()
                .any(|bases| bases.windows(2).any(|pair| pair[0] != pair[1])),
            "reference bases taken from another record of the same length must be visible"
        );

        // And at least one record whose head fields each need more than a byte.
        let (bytes, block_starts_at) = a_run_of_records(&records);
        let walked = walk_records(
            "the fixture's own heads",
            &bytes,
            block_starts_at,
            &RecordLayout::as_this_build_writes_it(),
            |_| false,
        );
        assert!(
            walked
                .iter()
                .any(|met| met.head.body_bytes > 255 && met.head.non_reference_reads > 255),
            "no record at the top of the depth range"
        );
    }

    /// **The oracle C3 exists for.** A walk that builds only some of a run's records gets, for
    /// each one it does build, exactly the record a walk that built them all gets — under every
    /// pattern, including the ones that skip the first record, the last, and all but one.
    ///
    /// The failure this catches is silent by construction: if any field of a body were coded as
    /// a difference from an earlier record, the records *after* a skipped one would still decode,
    /// into plausible and wrong values.
    #[test]
    fn a_walk_that_skips_records_builds_the_kept_ones_exactly() {
        let records = twelve_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);
        let layout = RecordLayout::as_this_build_writes_it();

        for (what, is_kept) in the_skip_patterns() {
            let walked = walk_records(what, &bytes, block_starts_at, &layout, is_kept);
            walk_matches(what, &walked, &records);
        }
    }

    /// **Every** skip pattern a twelve-record run has — all 4,096 of them, enumerated rather than
    /// sampled, so the test is deterministic and a defect that fires only under some particular
    /// pattern cannot hide. It costs a few hundredths of a second.
    #[test]
    fn every_skip_pattern_builds_exactly_the_records_it_keeps() {
        let records = twelve_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);
        let layout = RecordLayout::as_this_build_writes_it();

        for mask in 0u32..(1 << RECORDS_IN_THE_RUN) {
            let what = format!("mask {mask:#05x}");
            let walked = walk_records(&what, &bytes, block_starts_at, &layout, |index| {
                mask & (1u32 << index) != 0
            });
            walk_matches(&what, &walked, &records);
        }
    }

    /// Replace every skipped record's body with `filling`, leaving the heads alone, and hand back
    /// the damaged bytes together with how many bodies were overwritten.
    fn with_skipped_bodies_replaced(
        bytes: &[u8],
        block_starts_at: Position,
        is_kept: impl Fn(usize) -> bool,
        filling: impl Fn(usize) -> Vec<u8>,
    ) -> (Vec<u8>, usize) {
        let mut live_reads = a_live_set_reader();
        let mut damaged = bytes.to_vec();
        let mut at = 0usize;
        let mut measured_from = OffsetBase::at_block_start(block_starts_at);
        let mut index = 0usize;
        let mut replaced = 0usize;
        while at < bytes.len() {
            let found = read_a_record_head(&mut live_reads, &bytes[at..], A_CONTIG, measured_from)
                .expect("every head reads");
            if !is_kept(index) && !found.body.is_empty() {
                let body_at = at + found.record_bytes - found.body.len();
                let written = filling(found.body.len());
                assert_eq!(
                    written.len(),
                    found.body.len(),
                    "a replacement body fits exactly"
                );
                damaged[body_at..body_at + found.body.len()].copy_from_slice(&written);
                replaced += 1;
            }
            at += found.record_bytes;
            measured_from = OffsetBase::after(&found.head);
            index += 1;
        }
        (damaged, replaced)
    }

    /// A body of exactly `bytes` bytes that decodes: some reference bases, no observations, no
    /// reads either way, and the generic kind tag.
    fn a_decodable_body_of(bytes: usize) -> Option<Vec<u8>> {
        for prefix_bytes in 1..=5usize {
            let reference_bases = bytes.checked_sub(prefix_bytes + 4)?;
            let mut body = Vec::new();
            encode_u64_leb128(reference_bases as u64, &mut body);
            if body.len() != prefix_bytes {
                continue;
            }
            body.extend(std::iter::repeat_n(b'C', reference_bases));
            body.extend_from_slice(&[0, 0, 0, 0]);
            if body.len() == bytes {
                return Some(body);
            }
        }
        None
    }

    /// **A skipped record's body can be anything at all, and the records after it are
    /// unchanged.** This is the strongest form of "the body stands on its own", and it runs under
    /// every skip pattern and two kinds of damage.
    ///
    /// The two fillings catch different things and are cheap together. `0xff` is not a decodable
    /// body — every byte carries a variable-length integer's continuation bit, so such a body can
    /// never terminate — which turns *a reader touched this body* into a hard failure. A filling
    /// that **is** a valid body catches the shape the first cannot: a reader that touches a
    /// skipped body, decodes it successfully, and carries something out of it.
    ///
    /// If a body ever carried a difference from an earlier record, this is what would fail, and
    /// the walk would still consume the run exactly — because the *head* is what says how far to
    /// advance.
    #[test]
    fn a_skipped_records_body_can_be_anything_and_the_kept_records_are_unchanged() {
        let records = twelve_records_in_order();
        let (whole_bytes, block_starts_at) = a_run_of_records(&records);
        let layout = RecordLayout::as_this_build_writes_it();

        let fillings: [BodyDamage; 2] = [
            ("filled with a byte no body can end on", |len| {
                vec![0xff; len]
            }),
            ("filled with another record's body", |len| {
                a_decodable_body_of(len).unwrap_or_else(|| vec![0xff; len])
            }),
        ];

        for (what, is_kept) in the_skip_patterns() {
            for (damage, filling) in fillings {
                let (damaged_bytes, replaced) =
                    with_skipped_bodies_replaced(&whole_bytes, block_starts_at, is_kept, filling);
                let skipped = (0..records.len()).filter(|index| !is_kept(*index)).count();
                assert_eq!(
                    replaced, skipped,
                    "{what}, {damage}: every skipped body was replaced"
                );
                if replaced > 0 {
                    assert_ne!(damaged_bytes, whole_bytes, "{what}, {damage}");
                }

                let label = format!("{what}, {damage}");
                let walked =
                    walk_records(&label, &damaged_bytes, block_starts_at, &layout, is_kept);
                walk_matches(&label, &walked, &records);
            }
        }
    }

    /// **A body decodes the same alone as it does after every record before it.** Handed nothing
    /// but its own slice and its own region — no earlier record read, in reverse order so that no
    /// prior state could have been built up even by accident — each body produces the record a
    /// full forward walk produces.
    #[test]
    fn a_body_decodes_the_same_alone_as_it_does_after_every_record_before_it() {
        let mut live_reads = a_live_set_reader();
        let records = twelve_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);
        let layout = RecordLayout::as_this_build_writes_it();

        // Locate every record first, so each body can then be decoded with nothing else in hand.
        let mut located = Vec::new();
        let mut at = 0usize;
        let mut measured_from = OffsetBase::at_block_start(block_starts_at);
        while at < bytes.len() {
            let found = read_a_record_head(&mut live_reads, &bytes[at..], A_CONTIG, measured_from)
                .expect("the head reads");
            let body_at = at + found.record_bytes - found.body.len();
            located.push((body_at..body_at + found.body.len(), found.head));
            at += found.record_bytes;
            measured_from = OffsetBase::after(&found.head);
        }
        assert_eq!(located.len(), records.len(), "every record was located");

        for (index, (body, head)) in located.iter().enumerate().rev() {
            let alone = decode_record_body(
                &bytes[body.clone()],
                head.region,
                &the_live_reads_of(&records[index]),
                &layout,
            )
            .unwrap_or_else(|refused| {
                panic!("record {index}'s body decodes with no record before it: {refused}")
            });
            assert_eq!(
                alone.record,
                as_a_decode_can_return_it(&records[index]),
                "record {index}"
            );
            assert_eq!(alone.bytes_read, body.len(), "record {index}");
        }
    }

    /// **The skip the cohort's first pass actually makes.** It decides from the head's
    /// `non_reference_reads` rather than from a record's ordinal, and never reads the body of a
    /// position where every read matched the reference — about ninety-nine positions in a hundred
    /// on the corner that was measured, the tomato panel at three reads a position.
    #[test]
    fn a_walk_that_skips_by_what_the_head_says_builds_the_kept_ones_exactly() {
        let records = twelve_records_in_order();
        let (bytes, block_starts_at) = a_run_of_records(&records);
        let layout = RecordLayout::as_this_build_writes_it();

        // What the head will say, worked out from the records rather than from the file.
        let varies: Vec<bool> = records
            .iter()
            .map(|record| record.non_reference_and_compared_reads().0 > 100)
            .collect();
        assert!(
            varies.iter().any(|it| *it) && varies.iter().any(|it| !it),
            "the head has to separate the run, or this proves nothing"
        );

        let walked = walk_records(
            "skipping by what the head says",
            &bytes,
            block_starts_at,
            &layout,
            |index| varies[index],
        );
        for (index, met) in walked.iter().enumerate() {
            let (expected, _) = records[index].non_reference_and_compared_reads();
            assert_eq!(
                met.head.non_reference_reads, expected,
                "record {index}: the head a skipping reader filters on"
            );
        }
        walk_matches("skipping by what the head says", &walked, &records);
    }

    /// The same walk under a file written by a **later version** — one whose manifest declares
    /// fields this reader does not know, at the end of every body.
    ///
    /// **That is the configuration the next per-record field arrives in**, and it is the one
    /// spec psp_record_encoding.md §7 names as a trap: production stores a window's mean coverage
    /// as a difference from the previous record. A field walked past with a length that depended
    /// on an earlier record would fail here and nowhere else.
    #[test]
    fn a_walk_under_a_layout_with_unknown_trailing_fields_still_skips_cleanly() {
        let records = twelve_records_in_order();
        let mut fields = record_fields();
        fields.push(FieldSpec {
            name: FieldName("window-gc-percent".to_string()),
            encoding: FieldEncoding::FixedWidthInteger { width_bytes: 4 },
        });
        fields.push(FieldSpec {
            name: FieldName("window-mean-coverage".to_string()),
            encoding: FieldEncoding::LengthPrefixedBytes,
        });
        let layout = RecordLayout::from_manifest(&a_manifest_declaring(fields))
            .expect("two unknown fields at the end are not a reason to refuse the file");
        assert_eq!(layout.unknown_field_count(), 2);

        // Write the run as that later version would: each body followed by the two extra fields.
        let block_starts_at = records[0].region.start;
        let mut encoder = RecordEncoder::for_block(block_starts_at);
        let mut bytes = Vec::new();
        for (index, record) in records.iter().enumerate() {
            let mut one = Vec::new();
            encoder
                .encode_record(record, &mut one)
                .expect("the fixture is in coordinate order");
            let body_at = one.len() - record_body_length(record);
            let mut widened = one[..body_at].to_vec();
            let mut body = one[body_at..].to_vec();
            body.extend_from_slice(&[7, 7, 7, index as u8]);
            body.extend_from_slice(&[2, b'h', b'i']);
            // The head's declared body length has to grow with the body it describes, so the
            // head is rewritten rather than reused — **and the chain-id changes behind the count
            // are copied through verbatim**, because they are the head's too and this test is
            // pretending to be a later writer, not a different codec.
            let mut head = FieldReader::new(&widened);
            let offset = head.read_varint("position-offset").expect("the head reads");
            let span = head.read_varint("reference-span").expect("the head reads");
            let non_reference = head
                .read_varint("non-reference-reads")
                .expect("the head reads");
            head.read_varint("record-body-byte-count")
                .expect("the head reads");
            let changes = widened[head.bytes_read()..].to_vec();
            widened.clear();
            put_varint(&mut widened, offset);
            put_varint(&mut widened, span);
            put_varint(&mut widened, non_reference);
            put_varint(&mut widened, body.len() as u64);
            widened.extend_from_slice(&changes);
            widened.extend_from_slice(&body);
            bytes.extend_from_slice(&widened);
        }

        for (what, is_kept) in the_skip_patterns() {
            let walked = walk_records(what, &bytes, block_starts_at, &layout, is_kept);
            walk_matches(what, &walked, &records);
        }
    }
}
