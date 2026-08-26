//! One record's head: the fixed fields at the front of every record that let a reader
//! decide whether it wants the record without building it.
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
//! [`encode_body`] and [`decode_body`] are the other half of this file: one
//! [`SampleLocusObservations`] to bytes and back, exactly, with no compression, no file and
//! nothing read from outside the bytes themselves. The fields and their encodings are
//! [`BODY_FIELDS`], which is also what a writer declares in the header's manifest — so the
//! file's account of itself and the code that writes it are one list rather than two.
//!
//! **Two things are deliberately not in it yet.** A record's chain ids are dropped, because
//! they hold state across records and that is Milestone E of the plan. And the record's
//! coordinate is not in the body at all: it rides in the head, so [`decode_body`] is handed
//! the region rather than reading one.

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

/// Every field a record's **body** carries, in the order the bytes carry them.
///
/// **One list, from which both the manifest and the codec are written.** A header that
/// declares one field order while the encoder writes another is a file that decodes into
/// plausible nonsense, and nothing in the types would catch it — so the declaration a writer
/// puts in the header ([`record_body_fields`]) and the check a reader makes against the
/// header it found ([`RecordBodyLayout::from_manifest`]) both read this array, and the codec
/// below writes the fields in this order.
///
/// **The head's four fields are not here.** They are the fixed part in front of the body
/// (see the module doc) and they arrive with C2; this array is what a reader meets *after*
/// deciding it wants the record.
///
/// **The chain ids are not here either, and that is Milestone E.** A record's chain ids are
/// dropped by [`encode_body`] and come back empty from [`decode_body`], which is stated on
/// both and pinned by a test rather than left for a reader to discover.
const BODY_FIELDS: [(&str, FieldEncoding); 19] = [
    ("reference-bases", FieldEncoding::LengthPrefixedBytes),
    ("observation-count", FieldEncoding::Varint),
    ("observation-bases", FieldEncoding::LengthPrefixedBytes),
    ("witness-run-count", FieldEncoding::Varint),
    ("witness-run-start", FieldEncoding::Varint),
    ("witness-run-length", FieldEncoding::Varint),
    ("read-group", FieldEncoding::Varint),
    ("reads-showing-the-sequence", FieldEncoding::Varint),
    ("reads-on-the-forward-strand", FieldEncoding::Varint),
    ("summed-log-error", SUMMED_LOG_ERROR_ENCODING),
    ("mapq-sum", FieldEncoding::Varint),
    ("mapq-sum-of-squares", FieldEncoding::Varint),
    ("reads-placed-left", FieldEncoding::Varint),
    ("reads-without-observation", FieldEncoding::Varint),
    ("reads-discarded-by-cap", FieldEncoding::Varint),
    ("locus-kind", FieldEncoding::Varint),
    ("ssr-motif", FieldEncoding::LengthPrefixedBytes),
    ("ssr-left-flank", FieldEncoding::LengthPrefixedBytes),
    ("ssr-right-flank", FieldEncoding::LengthPrefixedBytes),
];

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

/// The locus-kind tags, as the bytes spell them. Adding a kind adds a tag and never
/// renumbers one, because a tag is on disk in every file already written.
const KIND_GENERIC: u64 = 0;
const KIND_SSR: u64 = 1;
const KIND_SSR_BUNDLE: u64 = 2;

/// The fewest bytes one observation can occupy — nine single-byte variable-length integers,
/// one per field, with an empty sequence and a complete witness. Used to bound what a
/// declared observation count may make this reader allocate, never to decide a length.
const LEAST_BYTES_PER_OBSERVATION: usize = 9;

/// The fewest bytes one witnessed run can occupy: its start and its length, one byte each.
const LEAST_BYTES_PER_WITNESS_RUN: usize = 2;

/// What a writer declares in the header's manifest for a record's body: every field, in
/// encoding order, with the encoding this version writes it in.
///
/// **The manifest is how a reader is driven by the file rather than by an assumption**
/// (spec §4.5), so this is what makes the two halves one decision instead of two.
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
/// nineteen-field declaration per record on a path that decodes about twenty million records
/// a second is work with one possible answer.
///
/// **What it carries is the part that differs between files: the fields this reader does not
/// know.** A later version of the writer may add a per-record scalar — the window's GC
/// fraction and its mean coverage are the two waiting to be computed — and it adds them
/// *after* every field named here. Each encoding in the closed set measures its own length,
/// so this reader walks past such a field without knowing anything about it: a
/// variable-length integer ends at its own last byte, a fixed-width one is its declared
/// width, and a byte string carries its length in front. **That is only true for a field
/// that appears once per record**; the manifest carries no cardinality
/// (`doc/devel/reports/implementations/ng_psp_a1_a2_2026-08-26.md` §2.8), so a repeated field
/// added later cannot be skipped and this reader refuses the file instead of guessing.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct RecordBodyLayout {
    /// The encodings of the trailing fields this reader does not recognise, in the order
    /// they appear at the end of every body. Empty for a file this version wrote.
    unknown_trailing: Vec<FieldEncoding>,
}

impl RecordBodyLayout {
    /// The layout this version of the code writes: every field it knows, and nothing after
    /// them.
    pub fn current() -> Self {
        Self {
            unknown_trailing: Vec::new(),
        }
    }

    /// The layout a file declares, checked against what this reader knows.
    ///
    /// **The known fields must come first, in this order, with these encodings.** A file that
    /// renames one, drops one, reorders two or declares one differently is refused — those
    /// are the shapes that would otherwise decode into plausible values rather than failing.
    /// Whatever the manifest lists *after* them is carried as something to walk past
    /// (see the type's own documentation).
    pub fn from_manifest(manifest: &Manifest) -> Result<Self, RecordLayoutError> {
        for (position, (name, encoding)) in BODY_FIELDS.iter().enumerate() {
            let expected = FieldName((*name).to_string());
            let Some(declared) = manifest.fields.get(position) else {
                return Err(RecordLayoutError::MissingField { position, expected });
            };
            if declared.name != expected {
                return Err(RecordLayoutError::UnexpectedField {
                    position,
                    expected,
                    found: declared.name.clone(),
                });
            }
            if declared.encoding != *encoding {
                return Err(RecordLayoutError::WrongEncoding {
                    field: expected,
                    expected: *encoding,
                    found: declared.encoding,
                });
            }
        }
        Ok(Self {
            unknown_trailing: manifest.fields[BODY_FIELDS.len()..]
                .iter()
                .map(|field| field.encoding)
                .collect(),
        })
    }

    /// How many fields this file carries that this reader does not know. Zero for a file
    /// this version wrote.
    pub fn unknown_field_count(&self) -> usize {
        self.unknown_trailing.len()
    }
}

/// Why a file's manifest cannot drive this reader.
///
/// **Every variant is an input problem, not a bug** — the same rule the rest of the module
/// keeps. Raised once per file, at open, never per record.
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
/// Two classes, because there are two instructions: the body **stops** before the record does
/// — which a reader meeting a straddled buffer must be able to retry rather than treat as
/// damage (Milestone D) — or the bytes **say something that cannot be true**, which is
/// corruption.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RecordDecodeError {
    /// The body ran out while this field was being read.
    #[error("the record's {field} runs past the end of its body")]
    Truncated { field: &'static str },
    /// The bytes were there and cannot mean what they say.
    #[error("the record's {field} is unreadable: {reason}")]
    Malformed { field: &'static str, reason: String },
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
/// records a second.
///
/// It cannot fail: every field of a [`SampleLocusObservations`] has a representation here,
/// and the variable-length integers are unbounded.
pub fn encode_body(record: &SampleLocusObservations, out: &mut Vec<u8>) {
    put_bytes(out, &record.reference_bases);
    put_varint(out, record.observations.len() as u64);
    for observation in &record.observations {
        put_bytes(out, &observation.bases);
        put_witness(out, &observation.read_witness);
        put_varint(out, u64::from(observation.read_group.get()));
        put_varint(out, u64::from(observation.num_obs));
        put_varint(out, u64::from(observation.num_fwd));
        encode_i64_svarint(observation.q_sum.steps(), out);
        put_varint(out, u64::from(observation.mapq_sum));
        put_varint(out, observation.mapq_sum_sq);
        put_varint(out, u64::from(observation.placed_left));
    }
    put_varint(out, u64::from(record.reads_without_observation));
    put_varint(out, u64::from(record.reads_discarded_by_cap));
    put_kind(out, &record.kind);
}

/// Read one record's body back, and say how many bytes it took.
///
/// `region` comes from the record's head, which is where the format keeps it — a body carries
/// no coordinate of its own (spec §4.3). `layout` comes from the file's manifest, once, at
/// open.
///
/// **The count of bytes consumed is worth checking against the head's `body_bytes`**: the two
/// disagreeing means the file and this reader disagree about the record's shape, which is
/// exactly what a version mismatch or a corrupt block looks like. C2 makes that check; here
/// the number is returned so it can.
///
/// **The chain-id lists come back empty** — see [`encode_body`].
pub fn decode_body(
    bytes: &[u8],
    region: GenomeRegion,
    layout: &RecordBodyLayout,
) -> Result<(SampleLocusObservations, usize), RecordDecodeError> {
    let mut body = BodyBytes::new(bytes);

    let reference_bases: Box<[u8]> = body.length_prefixed("reference bases")?.into();

    let declared_observations = body.varint("observation count")?;
    let mut observations = Vec::with_capacity(room_for(
        declared_observations,
        LEAST_BYTES_PER_OBSERVATION,
        body.bytes_left(),
    ));
    for _ in 0..declared_observations {
        let bases: Box<[u8]> = body.length_prefixed("observed bases")?.into();
        let read_witness = body.witness()?;
        let read_group = ReadGroupId(body.u32("read group")?);
        let num_obs = body.u32("count of reads showing the sequence")?;
        let num_fwd = body.u32("count of reads on the forward strand")?;
        let q_sum = SummedLogError::from_steps(body.signed_varint("summed log-error")?);
        let mapq_sum = body.u32("sum of mapping qualities")?;
        let mapq_sum_sq = body.varint("sum of squared mapping qualities")?;
        let placed_left = body.u32("count of reads placed left")?;
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

    let reads_without_observation = body.u32("count of reads that showed nothing")?;
    let reads_discarded_by_cap = body.u32("count of reads the depth cap discarded")?;
    let kind = body.kind()?;

    for encoding in &layout.unknown_trailing {
        body.skip(*encoding)?;
    }

    Ok((
        SampleLocusObservations {
            region,
            reference_bases,
            observations,
            reads_without_observation,
            reads_discarded_by_cap,
            kind,
        },
        body.read_so_far(),
    ))
}

// ---------------------------------------------------------------------
// Writing primitives
// ---------------------------------------------------------------------

fn put_varint(out: &mut Vec<u8>, value: u64) {
    encode_u64_leb128(value, out);
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
/// rather than a record written without its payload.
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
/// there first, so a truncated or hostile body produces a
/// [`Truncated`](RecordDecodeError::Truncated) rather than a panic.
struct BodyBytes<'a> {
    bytes: &'a [u8],
    read: usize,
}

impl<'a> BodyBytes<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, read: 0 }
    }

    fn read_so_far(&self) -> usize {
        self.read
    }

    fn bytes_left(&self) -> usize {
        self.bytes.len() - self.read
    }

    /// One variable-length integer, read through production's codec.
    ///
    /// **`decoded` is where both integer readers meet it**, so the two spellings of "the
    /// bytes ran out here" cannot drift apart: a `Truncated` from the codec is this record
    /// stopping early, which Milestone D's reader retries with more bytes, and an `Overflow`
    /// is a byte sequence no 64-bit value can have produced, which is damage.
    fn varint(&mut self, field: &'static str) -> Result<u64, RecordDecodeError> {
        self.decoded(decode_u64_leb128(&self.bytes[self.read..]), field)
    }

    /// One zig-zag variable-length integer — what a value that goes negative is written as.
    fn signed_varint(&mut self, field: &'static str) -> Result<i64, RecordDecodeError> {
        self.decoded(decode_i64_svarint(&self.bytes[self.read..]), field)
    }

    fn decoded<T>(
        &mut self,
        decoded: Result<(T, usize), VarintError>,
        field: &'static str,
    ) -> Result<T, RecordDecodeError> {
        match decoded {
            Ok((value, used)) => {
                self.read += used;
                Ok(value)
            }
            Err(VarintError::Truncated) => Err(RecordDecodeError::Truncated { field }),
            Err(VarintError::Overflow) => Err(RecordDecodeError::Malformed {
                field,
                reason: "a variable-length integer longer than any 64-bit value needs".to_string(),
            }),
        }
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, RecordDecodeError> {
        let value = self.varint(field)?;
        u32::try_from(value).map_err(|_| RecordDecodeError::Malformed {
            field,
            reason: format!("{value}, which is past the {} this field holds", u32::MAX),
        })
    }

    fn u16(&mut self, field: &'static str) -> Result<u16, RecordDecodeError> {
        let value = self.varint(field)?;
        u16::try_from(value).map_err(|_| RecordDecodeError::Malformed {
            field,
            reason: format!("{value}, which is past the {} this field holds", u16::MAX),
        })
    }

    fn take(&mut self, count: usize, field: &'static str) -> Result<&'a [u8], RecordDecodeError> {
        let to = self
            .read
            .checked_add(count)
            .filter(|to| *to <= self.bytes.len())
            .ok_or(RecordDecodeError::Truncated { field })?;
        let taken = &self.bytes[self.read..to];
        self.read = to;
        Ok(taken)
    }

    fn length_prefixed(&mut self, field: &'static str) -> Result<&'a [u8], RecordDecodeError> {
        let declared = self.varint(field)?;
        // The cast cannot lose anything that `take` would then accept: a length beyond
        // `usize` is beyond the buffer too, and `take` refuses it.
        let count = usize::try_from(declared).unwrap_or(usize::MAX);
        self.take(count, field)
    }

    fn witness(&mut self) -> Result<ReadWitness, RecordDecodeError> {
        let field = "witness";
        let declared_runs = self.varint("witness run count")?;
        if declared_runs == 0 {
            return Ok(ReadWitness::Complete);
        }
        let mut runs = Vec::with_capacity(room_for(
            declared_runs,
            LEAST_BYTES_PER_WITNESS_RUN,
            self.bytes_left(),
        ));
        for _ in 0..declared_runs {
            let start = self.u16("witness run start")?;
            let length = self.u16("witness run length")?;
            let end = start
                .checked_add(length)
                .ok_or_else(|| RecordDecodeError::Malformed {
                    field,
                    reason: format!(
                        "a run from {start} covering {length} positions ends past the \
                                     last locus position a witness can name"
                    ),
                })?;
            runs.push((start, end));
        }
        let positions = WitnessedLocusPositions::from_half_open_runs(runs).ok_or_else(|| {
            RecordDecodeError::Malformed {
                field,
                reason: "a partial witness that covers no position; only `Complete` is written \
                         with no runs"
                    .to_string(),
            }
        })?;
        Ok(ReadWitness::Partial { positions })
    }

    fn kind(&mut self) -> Result<LocusKind, RecordDecodeError> {
        let field = "locus kind";
        let tag = self.varint(field)?;
        match tag {
            KIND_GENERIC => Ok(LocusKind::Generic),
            KIND_SSR => {
                let motif_bases = self.length_prefixed("repeat motif")?;
                let motif =
                    Motif::new(motif_bases).map_err(|source| RecordDecodeError::Malformed {
                        field: "repeat motif",
                        reason: source.to_string(),
                    })?;
                let left_flank: Box<[u8]> = self.length_prefixed("left flank")?.into();
                let right_flank: Box<[u8]> = self.length_prefixed("right flank")?.into();
                Ok(LocusKind::Ssr(SsrDetail {
                    motif,
                    left_flank,
                    right_flank,
                }))
            }
            KIND_SSR_BUNDLE => Ok(LocusKind::SsrBundle),
            unknown => Err(RecordDecodeError::Malformed {
                field,
                reason: format!("kind {unknown}, which this reader does not know"),
            }),
        }
    }

    /// Walk past one field of a declared encoding without interpreting it — how a reader
    /// meets a field a later writer added (see [`RecordBodyLayout`]).
    ///
    /// **Exhaustive with no wildcard**, so an encoding added to the closed set has to say
    /// here how long one of its values is, rather than inheriting a guess.
    fn skip(&mut self, encoding: FieldEncoding) -> Result<(), RecordDecodeError> {
        let field = "a field this reader does not know";
        match encoding {
            FieldEncoding::Varint
            | FieldEncoding::SignedVarint
            | FieldEncoding::FixedPoint { .. } => {
                self.varint(field)?;
            }
            FieldEncoding::FixedWidthInteger { width_bytes }
            | FieldEncoding::IeeeFloat { width_bytes } => {
                self.take(usize::from(width_bytes), field)?;
            }
            FieldEncoding::LengthPrefixedBytes => {
                self.length_prefixed(field)?;
            }
        }
        Ok(())
    }
}

/// How much to reserve for a list whose length the file declares.
///
/// **Never the declared length on its own.** A corrupt or hostile body can say it holds four
/// billion observations in nine bytes, and a reader that believes it asks the allocator for
/// the number before reading a single one. Each entry costs at least a byte or two, so what
/// is left of the body bounds what can really be there; the loop still reads the declared
/// count and still fails when the bytes run out, which is where the wrong length is reported.
fn room_for(declared: u64, least_bytes_each: usize, bytes_left: usize) -> usize {
    let could_be_there = bytes_left / least_bytes_each;
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

    use crate::ng::locus_generation::{LocusKind, ReadWitness, SsrDetail, WitnessedLocusPositions};
    use crate::ng::psp::header::{
        DEFAULT_GENOMIC_BLOCK_SIZE_BP, DEFAULT_LOOK_BACK_WINDOW_LOG, Manifest,
    };
    use crate::ng::types::{Motif, ReadGroupId, SummedLogError};

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

    /// A record with every kind of thing a body carries: several observations, a complete
    /// witness and a partial one with a hole, two read groups, a negative error sum, a
    /// sequence longer than the reference's and one that is empty.
    fn a_rich_record() -> SampleLocusObservations {
        SampleLocusObservations {
            region: a_region(90_667_287, 90_667_293),
            reference_bases: b"ACGTAC".to_vec().into_boxed_slice(),
            observations: vec![
                SequenceObservation {
                    bases: b"ACGTAC".to_vec().into_boxed_slice(),
                    read_witness: ReadWitness::Complete,
                    read_group: ReadGroupId(0),
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
                    num_fwd: 0,
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
                    q_sum: SummedLogError::NONE,
                    mapq_sum: 60,
                    mapq_sum_sq: 3_600,
                    placed_left: 0,
                    chain_ids: Vec::new(),
                },
            ],
            reads_without_observation: 12,
            reads_discarded_by_cap: 305,
            kind: LocusKind::Generic,
        }
    }

    /// Encode, decode, and hand back both the record that came out and how many bytes it
    /// said it took.
    fn round_tripped(
        record: &SampleLocusObservations,
    ) -> (SampleLocusObservations, usize, Vec<u8>) {
        let mut bytes = Vec::new();
        encode_body(record, &mut bytes);
        let (back, read) = decode_body(&bytes, record.region, &RecordBodyLayout::current())
            .expect("what this encoder wrote, this decoder reads");
        (back, read, bytes)
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

    /// The whole of C1: what goes in comes back, field for field, and the decoder's own
    /// count of bytes read is the length of what the encoder wrote — which is the number
    /// C2's head will be checked against.
    #[test]
    fn a_generic_record_round_trips_field_for_field() {
        let written = a_rich_record();
        let (back, read, bytes) = round_tripped(&written);
        assert_eq!(back, written);
        assert_eq!(
            read,
            bytes.len(),
            "the decoder stopped where the encoder did"
        );
    }

    /// A record's reference bases are in the body, so decoding one needs no reference on
    /// hand. **This is the choice C1 took**: the record encoding spec leans to dropping them
    /// and re-fetching, on a measurement nobody can take until a writer and a reader exist,
    /// so they are stored and declared, and dropping them later is a manifest change rather
    /// than a format break. The fixture's reference is four bases over a one-base region —
    /// nothing about the region says what they are.
    #[test]
    fn the_reference_bases_are_in_the_body_so_a_decode_needs_no_reference() {
        let written = SampleLocusObservations {
            region: a_region(1_000, 1_000),
            reference_bases: b"ACGT".to_vec().into_boxed_slice(),
            observations: Vec::new(),
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        let (back, ..) = round_tripped(&written);
        assert_eq!(&*back.reference_bases, b"ACGT");
        assert_eq!(back, written);
    }

    /// A tract carries its motif and both flanks, and the flanks differ in length because a
    /// left flank is clamped at the contig's start.
    #[test]
    fn an_str_record_round_trips_with_its_motif_and_its_flanks() {
        let mut written = a_rich_record();
        written.kind = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide is a motif"),
            left_flank: b"GG".to_vec().into_boxed_slice(),
            right_flank: b"CCCCCCCCCC".to_vec().into_boxed_slice(),
        });
        let (back, ..) = round_tripped(&written);
        assert_eq!(back, written);
    }

    /// The third kind, which carries nothing, still has to come back as itself rather than
    /// as the kind whose tag is next to it.
    #[test]
    fn a_repeat_bundle_record_round_trips() {
        let mut written = a_rich_record();
        written.kind = LocusKind::SsrBundle;
        let (back, ..) = round_tripped(&written);
        assert_eq!(back.kind, LocusKind::SsrBundle);
        assert_eq!(back, written);
    }

    /// A covered position where every read agreed and nothing was recorded as a sequence is
    /// a real record, and it is most of a file.
    #[test]
    fn a_record_with_no_observations_round_trips() {
        let written = SampleLocusObservations {
            region: a_region(500, 500),
            reference_bases: b"A".to_vec().into_boxed_slice(),
            observations: Vec::new(),
            reads_without_observation: 4,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        let (back, read, bytes) = round_tripped(&written);
        assert_eq!(back, written);
        assert_eq!(read, bytes.len());
    }

    /// **The chain ids are Milestone E and this is the test that says so out loud.** A record
    /// carrying them encodes without complaint and comes back with empty lists and everything
    /// else identical — so the day E lands, this test is what fails if the ids are still
    /// being dropped.
    #[test]
    fn chain_ids_are_dropped_until_milestone_e_and_come_back_empty() {
        let mut written = a_rich_record();
        written.observations[0].chain_ids = vec![4, 17, 900_001];
        written.observations[1].chain_ids = vec![17];

        let (back, ..) = round_tripped(&written);

        assert!(
            back.observations.iter().all(|obs| obs.chain_ids.is_empty()),
            "C1 does not write chain ids; E1 to E4 are where they arrive"
        );
        assert_ne!(back, written, "so the round trip is not yet exact for them");

        let mut without_ids = written.clone();
        for observation in &mut without_ids.observations {
            observation.chain_ids.clear();
        }
        assert_eq!(back, without_ids, "and nothing else is touched");
    }

    /// The negative sums are the real ones — a sum of log error probabilities is at most
    /// zero — and the extremes are what a zig-zag encoding gets wrong if it is not one.
    #[test]
    fn a_summed_log_error_round_trips_negative_and_at_its_extremes() {
        for steps in [0, -1, 1, -4_096, i64::MIN, i64::MAX, -123_456_789] {
            let mut written = a_rich_record();
            written.observations[0].q_sum = SummedLogError::from_steps(steps);
            let (back, ..) = round_tripped(&written);
            assert_eq!(
                back.observations[0].q_sum.steps(),
                steps,
                "a summed log-error of {steps} steps"
            );
        }
    }

    /// A witness with a hole is the shape two numbers could not describe, and the runs have
    /// to come back in the same places rather than merged into their span.
    #[test]
    fn a_partial_witness_with_holes_round_trips() {
        let mut written = a_rich_record();
        written.observations[0].read_witness = ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 1), (3, 4), (9, 40)])
                .expect("three separated runs"),
        };
        let (back, ..) = round_tripped(&written);
        match &back.observations[0].read_witness {
            ReadWitness::Partial { positions } => {
                assert_eq!(
                    positions.runs().collect::<Vec<_>>(),
                    vec![(0, 1), (3, 4), (9, 40)]
                );
            }
            other => panic!("expected a partial witness, got {other:?}"),
        }
        assert_eq!(back, written);
    }

    // -----------------------------------------------------------------
    // The manifest and the codec are one decision
    // -----------------------------------------------------------------

    /// What a writer puts in the header is what this reader checks it against, so a file this
    /// version wrote has nothing in it this version does not know.
    #[test]
    fn the_manifest_a_writer_declares_is_the_one_this_reader_checks() {
        let layout = RecordBodyLayout::from_manifest(&a_manifest_declaring(record_body_fields()))
            .expect("the fields this writer declares are the fields this reader reads");
        assert_eq!(layout, RecordBodyLayout::current());
        assert_eq!(layout.unknown_field_count(), 0);
    }

    /// **The step is the type's, not the file's.** A file declaring a coarser step holds
    /// integers that mean something else, and reading them as if they were this type's would
    /// multiply every error sum in the file by four — so the file is refused rather than
    /// rescaled.
    #[test]
    fn a_manifest_that_declares_a_different_step_for_the_summed_log_error_is_refused() {
        let mut fields = record_body_fields();
        let position = fields
            .iter()
            .position(|field| field.name.0 == "summed-log-error")
            .expect("the field is declared");
        fields[position].encoding = FieldEncoding::FixedPoint {
            steps_per_unit: 1_024,
        };

        match RecordBodyLayout::from_manifest(&a_manifest_declaring(fields)) {
            Err(RecordLayoutError::WrongEncoding { field, found, .. }) => {
                assert_eq!(field.0, "summed-log-error");
                assert_eq!(
                    found,
                    FieldEncoding::FixedPoint {
                        steps_per_unit: 1_024
                    }
                );
            }
            other => panic!("expected WrongEncoding, got {other:?}"),
        }
    }

    /// Every field, one at a time: renamed, dropped, or given the encoding of the field next
    /// to it. All three are files whose records would otherwise decode into plausible values
    /// rather than failing.
    #[test]
    fn a_manifest_that_renames_reorders_or_drops_a_field_is_refused() {
        let declared = record_body_fields();
        for position in 0..declared.len() {
            let mut renamed = declared.clone();
            renamed[position].name = FieldName("something-else".to_string());
            assert!(
                RecordBodyLayout::from_manifest(&a_manifest_declaring(renamed)).is_err(),
                "field {position} renamed"
            );

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
            Err(RecordLayoutError::MissingField { position, expected }) => {
                assert_eq!(position, 2);
                assert_eq!(expected.0, BODY_FIELDS[2].0);
            }
            other => panic!("expected MissingField, got {other:?}"),
        }
    }

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
        encode_body(&record, &mut body);

        let later_versions = [
            (FieldEncoding::Varint, vec![0xac, 0x02]),
            (FieldEncoding::SignedVarint, vec![0xd7, 0x04]),
            (
                FieldEncoding::FixedWidthInteger { width_bytes: 4 },
                vec![1, 2, 3, 4],
            ),
            (FieldEncoding::IeeeFloat { width_bytes: 8 }, vec![9; 8]),
            (
                FieldEncoding::FixedPoint {
                    steps_per_unit: 100,
                },
                vec![0x40],
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

            let (back, read) = decode_body(&newer_body, record.region, &layout)
                .unwrap_or_else(|refused| panic!("{encoding:?} was not walked past: {refused}"));
            assert_eq!(back, record, "{encoding:?}");
            assert_eq!(
                read,
                newer_body.len(),
                "{encoding:?} — the whole body has to be accounted for, or the next record \
                 starts in the middle of this one"
            );
        }
    }

    // -----------------------------------------------------------------
    // Damaged and hostile bodies
    // -----------------------------------------------------------------

    /// A body cut at any point is refused at that point, and nothing in the decoder indexes
    /// past what it holds. **The cut is what a straddled buffer looks like**, so this is also
    /// the shape Milestone D's restartable parse will retry rather than treat as damage.
    #[test]
    fn a_body_cut_short_is_refused_at_every_cut_and_never_panics() {
        let record = a_rich_record();
        let mut whole = Vec::new();
        encode_body(&record, &mut whole);

        for cut in 0..whole.len() {
            let refused = decode_body(&whole[..cut], record.region, &RecordBodyLayout::current());
            assert!(refused.is_err(), "{cut} bytes of a record must not decode");
        }
        assert!(
            decode_body(&whole, record.region, &RecordBodyLayout::current()).is_ok(),
            "and the whole body still reads"
        );
    }

    /// Bytes after the record are not the record's business: a decoder stops where the record
    /// ends and says so, because in a block the next record starts there.
    #[test]
    fn a_body_followed_by_more_bytes_stops_where_the_record_ends() {
        let record = a_rich_record();
        let mut stream = Vec::new();
        encode_body(&record, &mut stream);
        let ends_at = stream.len();
        stream.extend_from_slice(&[0xff; 64]);

        let (back, read) = decode_body(&stream, record.region, &RecordBodyLayout::current())
            .expect("the record is complete; what follows it is another record's");
        assert_eq!(back, record);
        assert_eq!(read, ends_at);
    }

    /// **A count the file declares never sizes an allocation on its own.** Each fixture is a
    /// handful of bytes claiming to hold more observations — or more witnessed runs — than
    /// there are bytes in the universe. If the declared count reached the allocator, this
    /// test would not return.
    #[test]
    fn a_declared_count_larger_than_the_body_never_reaches_the_allocator() {
        let huge = {
            let mut bytes = Vec::new();
            encode_u64_leb128(u64::MAX, &mut bytes);
            bytes
        };

        // An empty reference, then an observation count of 2^64 − 1, then nothing.
        let mut observations_claimed = vec![0u8];
        observations_claimed.extend_from_slice(&huge);
        assert!(
            decode_body(
                &observations_claimed,
                a_region(1, 1),
                &RecordBodyLayout::current()
            )
            .is_err()
        );

        // An empty reference, one observation, an empty sequence, then that many runs.
        let mut runs_claimed = vec![0u8, 1, 0];
        runs_claimed.extend_from_slice(&huge);
        assert!(decode_body(&runs_claimed, a_region(1, 1), &RecordBodyLayout::current()).is_err());
    }

    /// A run covering no position is a set [`WitnessedLocusPositions`] refuses to hold, and a
    /// witness with no runs at all already means [`ReadWitness::Complete`] — so a body saying
    /// both is damaged, and is refused rather than turned into a complete witness.
    #[test]
    fn a_witness_run_covering_no_position_is_refused() {
        // An empty reference; one observation; empty bases; one run; start 4, length 0.
        let body = vec![0u8, 1, 0, 1, 4, 0];
        match decode_body(&body, a_region(1, 1), &RecordBodyLayout::current()) {
            Err(RecordDecodeError::Malformed { field, .. }) => assert_eq!(field, "witness"),
            other => panic!("expected a malformed witness, got {other:?}"),
        }
    }

    /// A kind tag this reader does not know is refused, not read as the kind before it — a
    /// file from a later writer must not have its repeat tracts silently turned into SNP
    /// candidates.
    #[test]
    fn an_unknown_locus_kind_is_refused() {
        let record = a_rich_record();
        let mut body = Vec::new();
        encode_body(&record, &mut body);
        let last = body.len() - 1;
        assert_eq!(
            body[last], 0,
            "the fixture's kind is Generic, which is tag 0"
        );
        body[last] = 99;

        match decode_body(&body, record.region, &RecordBodyLayout::current()) {
            Err(RecordDecodeError::Malformed { field, reason }) => {
                assert_eq!(field, "locus kind");
                assert!(reason.contains("99"), "got {reason}");
            }
            other => panic!("expected an unknown kind, got {other:?}"),
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

            match decode_body(&body, a_region(1, 1), &RecordBodyLayout::current()) {
                Err(RecordDecodeError::Malformed { field, .. }) => {
                    assert_eq!(field, "repeat motif", "for {motif_bases:?}");
                }
                other => panic!("expected a refused motif for {motif_bases:?}, got {other:?}"),
            }
        }
    }

    /// A count that does not fit the field it is read into is damage, and the message says
    /// which field and what the number was — not a silently truncated depth.
    #[test]
    fn a_count_too_large_for_its_field_is_refused_rather_than_narrowed() {
        // An empty reference, one observation, empty bases, a complete witness, read group 0,
        // then a read count of 2^32.
        let mut body = vec![0u8, 1, 0, 0, 0];
        encode_u64_leb128(u64::from(u32::MAX) + 1, &mut body);

        match decode_body(&body, a_region(1, 1), &RecordBodyLayout::current()) {
            Err(RecordDecodeError::Malformed { field, reason }) => {
                assert_eq!(field, "count of reads showing the sequence");
                assert!(reason.contains("4294967296"), "got {reason}");
            }
            other => panic!("expected a refused count, got {other:?}"),
        }
    }
}
