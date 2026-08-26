//! One record on the wire: the head a reader uses to decide whether it wants the record, and
//! the body it decodes once it has.
//!
//! # The head
//!
//! The fixed fields at the front of every record, which let a reader decide without building
//! anything.
//!
//! ```text
//! record = position_offset | reference_span | non_reference_reads | body_bytes
//!          | chain_id_changes | body
//!          └──────────────────────── the head ─────────────────────┘   └─ skip ─┘
//! ```
//!
//! **[`RecordHead`] is the fixed part of that head, and only the fixed part.** The chain-id
//! live-set changes — which reads arrived at this position and which left — are in the head
//! too, because they carry state a skipping reader must keep up to date or the merge composes
//! an allele for a read that was never there (spec psp_record_encoding.md §6). They are a
//! variable-length list, 6.42 bytes a position at 293 reads a position, so they are handed
//! straight to the reader's live set rather than stored in a `Copy` struct. Milestone E3 is
//! where they land; nothing in this type has to change for them.
//!
//! A reader takes the head, decides, and either builds the body or advances `body_bytes`
//! past it; nothing else in the block has to be touched to make that decision. **Measured
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
use crate::ng::psp::header::{FieldEncoding, FieldName, FieldSpec, Manifest};
use crate::ng::types::{GenomeRegion, Motif, ReadGroupId, SummedLogError};
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
    /// Reads at this locus that supported something other than the reference, summed over
    /// the observations. **Zero and "nothing varies here" are the same condition**, which
    /// is what the cohort's first pass filters on.
    ///
    /// A count of reads rather than of alternative alleles: it answers *does anything vary
    /// here* just as well, and it also lets a reader apply a threshold.
    pub non_reference_reads: u32,
    /// Length of the body that follows, in bytes.
    pub body_bytes: u32,
}

// ---------------------------------------------------------------------
// The body: one record's fields to bytes and back
// ---------------------------------------------------------------------

/// Every field a record's **body** carries, in encoding order.
///
/// **A field list, not a count of values per record.** Nine of these nineteen are written once
/// per *observation*, two of those nine once per *witness run* inside it, and the last three
/// only when the locus kind is a repeat tract. The manifest carries no cardinality
/// ([the A1+A2 report](../../../doc/devel/reports/implementations/ng_psp_a1_a2_2026-08-26.md)
/// §2.8), so what this list gives a reader is the order and the encodings; how many of each a
/// body holds is read from the counts the body itself carries.
///
/// **What it is for: a file's fingerprint of its own layout, and a reader's check against it.**
/// A writer declares this array ([`record_body_fields`]) and a reader refuses any file
/// declaring something else ([`RecordBodyLayout::from_manifest`]), so a file from another
/// version is rejected rather than decoded into plausible nonsense.
///
/// **It does not drive the codec, and the compiler cannot make it.** [`encode_record_body`] and
/// [`decode_record_body`] write and read these fields in this order by hand — which they must,
/// since nothing here says which of them repeat. What holds the two in step is
/// `the_fixture_encodes_to_these_exact_bytes`: reorder the codec and that test fails; reorder
/// this array and the manifest tests fail. **Change both together and you have changed the
/// format — raise the version.**
///
/// **The head's four fields are not here.** They are the fixed part in front of the body (see
/// the module doc) and they arrive with C2; this array is what a reader meets *after* deciding
/// it wants the record.
///
/// **The chain ids are not here either, and that is Milestone E.** A record's chain ids are
/// dropped by [`encode_record_body`] and come back empty from [`decode_record_body`], which is
/// stated on both and pinned by a test rather than left for a reader to discover.
const BODY_FIELDS: [(&str, FieldEncoding); 19] = [
    ("reference-bases", FieldEncoding::LengthPrefixedBytes),
    ("observation-count", FieldEncoding::Varint),
    // The next nine, once per observation:
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
const fn field(position: usize) -> &'static str {
    BODY_FIELDS[position].0
}

const REFERENCE_BASES: &str = field(0);
const OBSERVATION_COUNT: &str = field(1);
const OBSERVATION_BASES: &str = field(2);
const WITNESS_RUN_COUNT: &str = field(3);
const WITNESS_RUN_START: &str = field(4);
const WITNESS_RUN_LENGTH: &str = field(5);
const READ_GROUP: &str = field(6);
const READS_SHOWING_THE_SEQUENCE: &str = field(7);
const READS_ON_THE_FORWARD_STRAND: &str = field(8);
const SUMMED_LOG_ERROR: &str = field(9);
const MAPQ_SUM: &str = field(10);
const MAPQ_SUM_OF_SQUARES: &str = field(11);
const READS_STARTING_LEFT: &str = field(12);
const READS_WITHOUT_OBSERVATION: &str = field(13);
const READS_DISCARDED_BY_THE_DEPTH_CAP: &str = field(14);
const LOCUS_KIND: &str = field(15);
const REPEAT_MOTIF: &str = field(16);
const REPEAT_LEFT_FLANK: &str = field(17);
const REPEAT_RIGHT_FLANK: &str = field(18);

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
const MOST_BYTES_A_BODY_CAN_DECLARE: u64 = u32::MAX as u64;

/// What a writer declares in the header's manifest for a record's body: every field, in
/// encoding order, with the encoding this version writes it in.
///
/// **The manifest is how a reader is driven by the file rather than by an assumption**
/// (spec §4.5). It is a fingerprint rather than a parsing recipe — see [`BODY_FIELDS`] for what
/// it does and does not promise.
pub fn record_body_fields() -> Vec<FieldSpec> {
    BODY_FIELDS
        .iter()
        .map(|(name, encoding)| FieldSpec {
            name: FieldName((*name).to_string()),
            encoding: *encoding,
        })
        .collect()
}

/// How this reader must read the bodies in one particular file.
///
/// Built once per file from its manifest and then used for every record, because checking a
/// nineteen-field declaration per record on a path that decodes about twenty million records a
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
/// from one written once per observation.** The manifest carries no cardinality, so a file
/// whose extra field repeats per observation is *accepted* and decoded into plausible nonsense
/// from the second observation onwards — it is not refused, and nothing here can refuse it. A
/// later writer adding a per-observation field must raise the format version, which a reader
/// does refuse (`header.rs`, `UnsupportedVersion`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordBodyLayout {
    /// The encodings of the fields this reader does not recognise, in the order they appear at
    /// the end of every body. Empty for a file this version wrote.
    ///
    /// The declared *names* are not kept: this module's errors name a field with a
    /// `&'static str` and a later writer's names are not static. A fault in one of these
    /// carries the byte offset it happened at, and they are walked in declared order, which is
    /// what identifies it.
    unknown_final_encodings: Vec<FieldEncoding>,
}

impl RecordBodyLayout {
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
    pub fn from_manifest(manifest: &Manifest) -> Result<Self, RecordBodyLayoutError> {
        for (position, (name, encoding)) in BODY_FIELDS.iter().enumerate() {
            let Some(declared) = manifest.fields.get(position) else {
                return Err(RecordBodyLayoutError::MissingField {
                    position,
                    expected: FieldName((*name).to_string()),
                });
            };
            if declared.name.0 != *name {
                return Err(RecordBodyLayoutError::UnexpectedField {
                    position,
                    expected: FieldName((*name).to_string()),
                    found: declared.name.clone(),
                });
            }
            if declared.encoding != *encoding {
                return Err(RecordBodyLayoutError::WrongEncoding {
                    field: declared.name.clone(),
                    expected: *encoding,
                    found: declared.encoding,
                });
            }
        }
        Ok(Self {
            unknown_final_encodings: manifest
                .fields
                .iter()
                .skip(BODY_FIELDS.len())
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
pub enum RecordBodyLayoutError {
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
    /// The body ran out while this field was being read, and what it declared was possible.
    #[error("the record's {field} runs past the end of its body, {bytes_in} bytes in")]
    Truncated {
        field: &'static str,
        /// How far into the body the reader had got. **A record at 300 reads a position holds
        /// hundreds of observations**, so a field name alone does not say where to look.
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

/// Append `record`'s body to `out`, and never fail.
///
/// **The chain ids are dropped.** [`SequenceObservation::chain_ids`] is Milestone E of
/// `doc/devel/ng/impl_plan/psp_file_format.md`, so a record encoded here and read back has
/// empty lists where it had ids, and everything else identical. Encoding a record that
/// carries them is not an error today and produces a file the cohort merge could not use;
/// the writer that produces real files arrives after E does.
///
/// **Appends rather than returning a buffer**, so a writer keeps one and clears it between
/// records instead of allocating per record on a path that runs at about twenty million
/// records a second (spec §4.5).
///
/// It cannot fail. Every field of a [`SampleLocusObservations`] has a representation here and
/// the variable-length integers are unbounded; the one subtraction, a witnessed run's length,
/// cannot underflow because [`WitnessedLocusPositions`] keeps its runs private and its
/// constructors reject a run that covers no position.
pub fn encode_record_body(record: &SampleLocusObservations, out: &mut Vec<u8>) {
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
    for observation in observations {
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
            chain_ids: _, // Milestone E; named rather than swept up by `..`, so a rename shows
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
/// ([`RecordBodyLayout::from_manifest`]); [`RecordBodyLayout::as_this_build_writes_it`] is for
/// bytes this process encoded itself and skips every check that function makes.
///
/// **The chain-id lists come back empty** — see [`encode_record_body`].
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
pub fn decode_record_body(
    bytes: &[u8],
    region: GenomeRegion,
    layout: &RecordBodyLayout,
) -> Result<DecodedRecordBody, RecordDecodeError> {
    let mut body = BodyReader::new(bytes);

    let reference_bases: Box<[u8]> = body.read_length_prefixed(REFERENCE_BASES)?.into();

    let declared_observations = body.read_count(OBSERVATION_COUNT, LEAST_BYTES_PER_OBSERVATION)?;
    let mut observations = Vec::with_capacity(entries_to_reserve(
        declared_observations,
        LEAST_BYTES_PER_OBSERVATION,
        body.bytes_left(),
    ));
    for _ in 0..declared_observations {
        let bases: Box<[u8]> = body.read_length_prefixed(OBSERVATION_BASES)?.into();
        let read_witness = body.read_witness()?;
        let read_group = ReadGroupId(body.read_u32(READ_GROUP)?);
        let num_obs = body.read_u32(READS_SHOWING_THE_SEQUENCE)?;
        let num_fwd = body.read_u32(READS_ON_THE_FORWARD_STRAND)?;
        let q_sum = SummedLogError::from_steps(body.read_signed_varint(SUMMED_LOG_ERROR)?);
        let mapq_sum = body.read_u32(MAPQ_SUM)?;
        let mapq_sum_sq = body.read_varint(MAPQ_SUM_OF_SQUARES)?;
        let placed_left = body.read_u32(READS_STARTING_LEFT)?;
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
            chain_ids: Vec::new(),
        });
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

fn put_varint(out: &mut Vec<u8>, value: u64) {
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

/// A record body being read, and how far into it the reader has got.
///
/// **It never indexes past what it holds**: every method that advances checks the bytes are
/// there first, so a truncated or hostile body produces an error rather than a panic. The
/// invariant every method keeps is `bytes_read <= bytes.len()`, which is what makes
/// [`bytes_left`](Self::bytes_left) and the slicing below total.
///
/// Every method here advances the cursor, which is why each is named `read_` or `skip_` — a
/// column of them in [`decode_record_body`] has to say that it consumes as it goes.
struct BodyReader<'a> {
    bytes: &'a [u8],
    bytes_read: usize,
}

impl<'a> BodyReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self {
            bytes,
            bytes_read: 0,
        }
    }

    fn bytes_read(&self) -> usize {
        self.bytes_read
    }

    fn bytes_left(&self) -> usize {
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

    fn malformed(&self, field: &'static str, reason: String) -> RecordDecodeError {
        RecordDecodeError::Malformed {
            field,
            bytes_in: self.bytes_read,
            reason,
        }
    }

    /// One variable-length integer, read through production's codec.
    fn read_varint(&mut self, field: &'static str) -> Result<u64, RecordDecodeError> {
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
    fn read_count(
        &mut self,
        field: &'static str,
        least_bytes_each: usize,
    ) -> Result<u64, RecordDecodeError> {
        let declared = self.read_varint(field)?;
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
    /// meets a field a later writer added (see [`RecordBodyLayout`]).
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
fn entries_to_reserve(declared: u64, least_bytes_each: usize, bytes_left: usize) -> usize {
    let could_be_there = (bytes_left / least_bytes_each).min(MOST_ENTRIES_RESERVED);
    declared.min(could_be_there as u64) as usize
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position};

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
            reference_bases: b"ACGTACG".to_vec().into_boxed_slice(),
            observations: vec![
                SequenceObservation {
                    bases: b"ACGTACG".to_vec().into_boxed_slice(),
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
                    bases: b"ACGTACGTACGTAC".to_vec().into_boxed_slice(),
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
                    bases: Box::from(&b""[..]),
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

    /// Encode, decode, and hand back what came out together with the bytes that went in.
    fn round_trip(record: &SampleLocusObservations) -> (DecodedRecordBody, Vec<u8>) {
        let mut bytes = Vec::new();
        encode_record_body(record, &mut bytes);
        let decoded = decode_record_body(
            &bytes,
            record.region,
            &RecordBodyLayout::as_this_build_writes_it(),
        )
        .expect("what this encoder wrote, this decoder reads");
        (decoded, bytes)
    }

    /// Decode with the layout this build writes — the form every hostile-input test wants.
    fn decoded(bytes: &[u8]) -> Result<DecodedRecordBody, RecordDecodeError> {
        decode_record_body(
            bytes,
            a_region(1, 1),
            &RecordBodyLayout::as_this_build_writes_it(),
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
            reference_bases: b"ACGT".to_vec().into_boxed_slice(),
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

    /// **The chain ids are Milestone E and this is the test that says so out loud.** A record
    /// carrying them encodes without complaint and comes back with empty lists and everything
    /// else identical — so the day E lands, this test is what fails if the ids are still
    /// being dropped.
    #[test]
    fn chain_ids_come_back_empty_and_nothing_else_changes() {
        let mut written = a_rich_record();
        written.observations[0].chain_ids = vec![4, 17, 900_001];
        written.observations[1].chain_ids = vec![17];

        let (decoded, _) = round_trip(&written);

        assert!(
            decoded
                .record
                .observations
                .iter()
                .all(|obs| obs.chain_ids.is_empty()),
            "C1 does not write chain ids; E1 to E4 are where they arrive"
        );
        assert_ne!(
            decoded.record, written,
            "so the round trip is not yet exact"
        );

        let mut without_ids = written.clone();
        for observation in &mut without_ids.observations {
            observation.chain_ids.clear();
        }
        assert_eq!(decoded.record, without_ids, "and nothing else is touched");
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
        written.reference_bases = vec![b'N'; 300].into_boxed_slice();
        written.observations[0].bases = (0u8..=255).collect::<Vec<_>>().into_boxed_slice();
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
    fn record_body_fields_declares_the_names_a_written_file_carries() {
        let names: Vec<String> = record_body_fields()
            .into_iter()
            .map(|field| field.name.0)
            .collect();
        assert_eq!(
            names,
            vec![
                "reference-bases",
                "observation-count",
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
                "reads-without-observation",
                "reads-discarded-by-the-depth-cap",
                "locus-kind",
                "repeat-motif",
                "repeat-left-flank",
                "repeat-right-flank",
            ]
        );
    }

    /// What a writer puts in the header is what this reader checks it against, so a file this
    /// version wrote has nothing in it this version does not know.
    #[test]
    fn the_manifest_a_writer_declares_is_the_one_this_reader_checks() {
        let layout = RecordBodyLayout::from_manifest(&a_manifest_declaring(record_body_fields()))
            .expect("the fields this writer declares are the fields this reader reads");
        assert_eq!(layout, RecordBodyLayout::as_this_build_writes_it());
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
            manifest: a_manifest_declaring(record_body_fields()),
        };

        let bytes = header
            .encode()
            .expect("a writer's own manifest must be a legal header");
        let (back, _) = Header::decode(&bytes, std::path::Path::new("SRR7279481.psp"))
            .expect("and its own reader must read it");
        assert_eq!(back.manifest.fields, record_body_fields());
        assert_eq!(
            RecordBodyLayout::from_manifest(&back.manifest).expect("and drive this reader"),
            RecordBodyLayout::as_this_build_writes_it()
        );
    }

    /// **The step is the type's, not the file's, and this is the seam between them.** The
    /// declaration exists so a reader can turn the stored integer back into a quantity; a
    /// writer declaring anything else would put every error sum in the file out by the ratio.
    #[test]
    fn the_declared_step_is_the_types_step_and_is_four_thousand_and_ninety_six() {
        assert_eq!(SummedLogError::STEPS_PER_NAT, 4_096);
        assert_eq!(
            record_body_fields()[9].encoding,
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
        let mut fields = record_body_fields();
        let position = fields
            .iter()
            .position(|field| field.name.0 == "summed-log-error")
            .expect("the field is declared");
        fields[position].encoding = FieldEncoding::FixedPoint {
            steps_per_unit: quarter_of_the_real_step,
        };

        match RecordBodyLayout::from_manifest(&a_manifest_declaring(fields)) {
            Err(RecordBodyLayoutError::WrongEncoding { field, found, .. }) => {
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
        let declared = record_body_fields();
        for position in 0..declared.len() {
            let mut renamed = declared.clone();
            renamed[position].name = FieldName("something-else".to_string());
            match RecordBodyLayout::from_manifest(&a_manifest_declaring(renamed)) {
                Err(RecordBodyLayoutError::UnexpectedField {
                    position: at,
                    expected,
                    found,
                }) => {
                    assert_eq!(at, position);
                    assert_eq!(expected.0, BODY_FIELDS[position].0);
                    assert_eq!(found.0, "something-else");
                }
                other => panic!("field {position} renamed: expected UnexpectedField, {other:?}"),
            }

            let mut dropped = declared.clone();
            dropped.remove(position);
            assert!(
                RecordBodyLayout::from_manifest(&a_manifest_declaring(dropped)).is_err(),
                "field {position} dropped"
            );

            if position + 1 < declared.len() {
                let mut swapped = declared.clone();
                swapped.swap(position, position + 1);
                assert!(
                    RecordBodyLayout::from_manifest(&a_manifest_declaring(swapped)).is_err(),
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
        let mut fields = record_body_fields();
        fields.truncate(2);
        match RecordBodyLayout::from_manifest(&a_manifest_declaring(fields)) {
            Err(RecordBodyLayoutError::MissingField { position, expected }) => {
                assert_eq!(position, 2);
                assert_eq!(expected.0, BODY_FIELDS[2].0);
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
            let mut fields = record_body_fields();
            fields.push(FieldSpec {
                name: FieldName("a-field-from-a-later-writer".to_string()),
                encoding,
            });
            let layout = RecordBodyLayout::from_manifest(&a_manifest_declaring(fields))
                .expect("an unknown field at the end is not a reason to refuse the file");
            assert_eq!(layout.unknown_field_count(), 1);

            let mut newer_body = body.clone();
            newer_body.extend_from_slice(&written);

            let decoded = decode_record_body(&newer_body, record.region, &layout)
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

        let mut fields = record_body_fields();
        fields.push(FieldSpec {
            name: FieldName("window-gc-percent".to_string()),
            encoding: FieldEncoding::FixedWidthInteger { width_bytes: 4 },
        });
        fields.push(FieldSpec {
            name: FieldName("window-mean-coverage".to_string()),
            encoding: FieldEncoding::LengthPrefixedBytes,
        });
        let layout = RecordBodyLayout::from_manifest(&a_manifest_declaring(fields))
            .expect("two unknown fields at the end are fine");
        assert_eq!(layout.unknown_field_count(), 2);

        body.extend_from_slice(&[7, 7, 7, 7, 2, b'h', b'i']);
        let decoded =
            decode_record_body(&body, record.region, &layout).expect("both are walked past");
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
        let mut fields = record_body_fields();
        fields.push(FieldSpec {
            name: FieldName("window-mean-coverage".to_string()),
            encoding: FieldEncoding::IeeeFloat { width_bytes: 8 },
        });
        let layout = RecordBodyLayout::from_manifest(&a_manifest_declaring(fields))
            .expect("an unknown field at the end is fine");

        let whole = body.len();
        body.extend_from_slice(&[1, 2, 3]);
        for cut in whole..body.len() {
            match decode_record_body(&body[..cut], record.region, &layout) {
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
                &RecordBodyLayout::as_this_build_writes_it(),
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
                &RecordBodyLayout::as_this_build_writes_it()
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
                let mut body = vec![0u8, 1, 0];
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

        let mut runs_claimed = vec![0u8, 1, 0];
        encode_u64_leb128(u64::MAX, &mut runs_claimed);
        assert!(decoded(&runs_claimed).is_err());
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
            reference_bases: Box::from(&b""[..]),
            observations: Vec::new(),
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        let cheapest_observation = SequenceObservation {
            bases: Box::from(&b""[..]),
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
            ..cheapest_observation
        });

        let bytes_of = |record: &SampleLocusObservations| {
            let mut bytes = Vec::new();
            encode_record_body(record, &mut bytes);
            bytes.len()
        };
        assert_eq!(
            bytes_of(&with_one) - bytes_of(&empty),
            LEAST_BYTES_PER_OBSERVATION
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
            &RecordBodyLayout::as_this_build_writes_it(),
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
        let body = vec![0u8, 1, 0, 1, 4, 0];
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
        let out_of_order = vec![0u8, 1, 0, 2, 5, 2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let overlapping = vec![0u8, 1, 0, 2, 0, 5, 3, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let touching = vec![0u8, 1, 0, 2, 0, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
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
            let mut body = vec![0u8, 1, 0, 1];
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
        let mut body = vec![0u8, 1, 0, 1];
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
            &RecordBodyLayout::as_this_build_writes_it(),
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
            // An empty reference, no observations, no reads either way, kind 1, then the
            // motif and two empty flanks.
            let mut body = vec![0u8, 0, 0, 0, 1];
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
        let mut body = vec![0u8, 1, 0, 0, 0];
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
                        bases: bases.into_boxed_slice(),
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
                reference_bases: reference.into_boxed_slice(),
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
                &RecordBodyLayout::as_this_build_writes_it(),
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
}
