//! The census on disk: a header, a directory, and the sections.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_records.md` §6.1 and §6.2; types in
//! `doc/devel/ng/arch/parameter_prepass_joint_records.md` §1.1a and §2.2.
//!
//! **One file per sample, beside that sample's pileup and never inside it** (spec §6.1). It is a
//! cache: everything in it can be recomputed from the pileup, and what it saves is a full
//! decompression pass over that pileup every time a cohort is fitted.
//!
//! # Why there is a directory
//!
//! The parameters fit never needs a whole sample's evidence at once. It finishes the ordinary
//! positions before it reads a tract, and within the tracts it fits one band of strata at a
//! time — so the smallest piece anything asks for is *one read group's ordinary positions* or
//! *one read group's tracts for one stratum*. The directory gives each of those an offset and a
//! length, so a reader can take one and decode nothing else.
//!
//! **The header is outside the sections on purpose.** The twelve recording terms and the
//! kept-loci digest are compared across every sample *before* anything large is decoded (spec
//! §5), so a cohort of a thousand samples that cannot be pooled says so after reading a few
//! hundred bytes each.
//!
//! # The layout
//!
//! Every integer is little-endian. Offsets are from the start of the file.
//!
//! ```text
//!   magic     8 bytes   "NGCENSUS"
//!   version   u16       1
//!   header              the sample's name, the twelve terms, and the pileup it was built from
//!   directory u32 n, then n × (section key, offset u64, length u64)
//!   sections            the bytes each directory entry points at
//! ```
//!
//! **The directory is written before the sections and holds absolute offsets**, so a reader
//! seeks once per section. That means the writer has to know each section's length before it
//! writes the directory, which it does by encoding the sections first and then placing them.
//!
//! # What is not here yet
//!
//! **Reading one section on its own.** This module writes the file and reads it whole; the
//! seeking reader that fills one section at a time is the next unit of work, and the whole-file
//! read here is its parity oracle.

use md5::{Digest, Md5};

use std::io::{Read, Write};
use std::path::Path;

use crate::ng::parameter_estimation::joint::loci::{BlockDigest, CensusLociDigest};
use crate::ng::repeat_catalog::StratumCounts;
use crate::ng::types::{ContigId, ReadGroupId};
use std::collections::BTreeMap;

use super::census::{
    AlleleObservation, ByteExtent, CensusError, DEPTH_CODE_BITS, DepthCap, DepthLadderDigest,
    GenericEvidence, GuardObservation, NamedReadGroup, OFFSET_BUCKETS, ObservedAllele,
    OffsetCounts, PackedDepthCodes, ReadCap, RecordingTerms, SampleCensusEvidence, Section,
    SectionKey, SelectionTermsDigest, SsrEvidence, Stratum, TractDifference, WalkedBits,
};

/// What every census file starts with, so a file that is not one is refused rather than decoded.
const MAGIC: &[u8; 8] = b"NGCENSUS";

/// The layout this build writes and the only one it reads.
///
/// **A version and not a feature flag.** A census is a cache with a pileup behind it, so the
/// answer to a version this build does not know is to rebuild rather than to interpret.
///
/// **2 since 2026-08-16**, when the depth code went from five bits on a widening ladder to
/// eight on one with a bin for every depth to the cap. A version-1 file's depth array is a
/// different number of bytes for the same position count and its codes index a different
/// ladder, so nothing about it can be salvaged by reading it more carefully.
const VERSION: u16 = 3;

/// Which pileup a census was built from.
///
/// **A digest and a count, never a modification time** (spec §6.1). A modification time changes
/// when a file is copied and does not change when its contents are rewritten in place, so it
/// answers a question nobody asked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PileupIdentity {
    /// A digest of the pileup's header — its reference, its analysed regions, its read filters
    /// and the command line that produced it.
    pub header: [u8; 16],
    /// How many records that pileup holds.
    pub records: u64,
}

impl PileupIdentity {
    /// The identity of a pileup whose header is `header` and which holds `records` records.
    ///
    /// **The bytes are the pileup's own header** — its reference, its analysed regions, its read
    /// filters and the command line that produced it (spec §6.1). Which bytes exactly is the
    /// pileup writer's business; what this promises is that two pileups with the same header and
    /// the same record count get the same identity and no others do.
    pub fn of_header(header: &[u8], records: u64) -> Self {
        let mut hasher = Md5::new();
        hasher.update(header);
        Self {
            header: hasher.finalize().into(),
            records,
        }
    }
}

/// What a run should do with a census file, given the pileup it means to fit from.
///
/// **A verdict and not an action.** Rebuilding a census from a pileup is the second producer
/// (milestone C); until that exists, this says what should happen and the caller decides.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Freshness {
    /// The census was built from this pileup. Use it.
    Fresh,
    /// It was built from a different pileup, and the pileup in hand can rebuild it. The value
    /// names what differs, so that a run can say why it spent the time.
    Rebuild(&'static str),
    /// It was built from a different pileup and there is none here to rebuild from. **Refuse**:
    /// a census whose evidence came from other reads is not this run's evidence, and the whole
    /// point of naming the pileup is that the difference is otherwise invisible.
    Refused(&'static str),
}

/// Whether a census may be used as it stands, rebuilt, or refused (spec §6.1).
///
/// `named` is what the census file says it was built from and `in_hand` is the pileup this run
/// has, `None` when it cannot be reached.
///
/// **Nothing here reads a modification time**, which is the point: a modification time changes
/// when a file is copied and does not change when its contents are rewritten in place, so it
/// answers a question nobody asked.
pub fn freshness(named: Option<PileupIdentity>, in_hand: Option<PileupIdentity>) -> Freshness {
    match (named, in_hand) {
        (Some(named), Some(here)) if named == here => Freshness::Fresh,
        (Some(named), Some(here)) => {
            let field = if named.header == here.header {
                "the pileup's record count"
            } else {
                "the pileup's header"
            };
            Freshness::Rebuild(field)
        }
        // A census that names a pileup, with none here: nothing can rebuild it and nothing can
        // check it, so it is refused rather than trusted.
        (Some(_), None) => Freshness::Refused("the pileup it was built from, which is not here"),
        // A census that names no pileup at all — written by a run that had none. It cannot be
        // checked, so it is rebuilt where that is possible and refused where it is not.
        (None, Some(_)) => {
            Freshness::Rebuild("the pileup it was built from, which it does not name")
        }
        (None, None) => Freshness::Refused("the pileup it was built from, which it does not name"),
    }
}

/// **Whether a census may be used, judged on its psp's header alone** — which is all a fit can
/// afford to check.
///
/// [`freshness`] compares the whole identity: the digest of the psp's header *and* how many
/// records it holds. A fit can get the first for the price of one short read
/// ([`psp::header_digest`](crate::ng::psp::header_digest)). **It cannot cheaply get the second.**
/// A psp's footer carries its block count and its byte offsets and no record count, and the
/// per-block counts sit inside each block's compressed stream — so obtaining one means
/// decompressing the psp, at every fit, for every sample. That is the cost psp mode exists to
/// avoid.
///
/// **What that leaves unchecked, precisely**: a psp whose header is unchanged and whose records
/// are not. Only [`PspWriter::append`](crate::ng::psp::PspWriter::append) can produce one — it
/// reopens a finished file and carries the header forward — and nothing in the shipped commands
/// calls it. **The fix is to put a record count in the psp's footer**, which makes both halves
/// one cheap read; that is a psp format change and belongs with the walk stage rather than here.
///
/// `named` is what the census says it was built from and `header_in_hand` is the digest of the
/// psp beside it, `None` when there is no psp there.
#[must_use]
pub fn freshness_by_header(
    named: Option<PileupIdentity>,
    header_in_hand: Option<[u8; 16]>,
) -> Freshness {
    match (named, header_in_hand) {
        (Some(named), Some(here)) if named.header == here => Freshness::Fresh,
        (Some(_), Some(_)) => Freshness::Rebuild("the pileup's header"),
        (Some(_), None) => Freshness::Refused("the pileup it was built from, which is not here"),
        (None, Some(_)) => {
            Freshness::Rebuild("the pileup it was built from, which it does not name")
        }
        (None, None) => Freshness::Refused("the pileup it was built from, which it does not name"),
    }
}

thread_local! {
    /// Bytes this thread has read out of census files, section by section.
    static BYTES_READ: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
}

/// Count `bytes` against this thread's total — called once a section read.
pub(super) fn count_bytes_read(bytes: u64) {
    BYTES_READ.with(|counter| counter.set(counter.get() + bytes));
}

/// How many bytes of census file this thread has read since [`reset_bytes_read`].
///
/// **This is spec §7.15's counting reader, and it is the half worth having.** An
/// implementation that decoded a whole file and handed back a slice would match every value a
/// section-by-section reader gives and deliver none of the memory the by-section design exists
/// for; only the byte count tells them apart.
///
/// **Per thread, because a read happens on the thread that asked for it** — so a test measuring
/// its own calls is not measuring another test's.
pub fn bytes_read() -> u64 {
    BYTES_READ.with(std::cell::Cell::get)
}

/// Start counting again from zero.
pub fn reset_bytes_read() {
    BYTES_READ.with(|counter| counter.set(0));
}

/// A census file, read whole.
#[derive(Debug)]
pub struct CensusFile {
    pub census: SampleCensusEvidence,
    /// Absent where the census was built during a walk that wrote no pileup.
    pub pileup: Option<PileupIdentity>,
}

// ---------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------

/// Write one sample's census, with the identity of the pileup it was built from.
///
/// # Errors
///
/// [`CensusError::Io`] and nothing else: every value here is already well-formed, since the
/// types that carry them refuse to hold anything that is not.
pub fn write_census(
    census: &SampleCensusEvidence,
    pileup: Option<PileupIdentity>,
    out: &mut impl Write,
) -> Result<(), CensusError> {
    // **The sections are encoded before the directory is written**, because the directory holds
    // each one's length and its absolute offset, and neither is known until the bytes exist.
    let mut sections: Vec<(SectionKey, Vec<u8>)> = Vec::new();
    for (group, records) in census.generic_sections() {
        sections.push((SectionKey::Generic(group), encode_generic(records)));
    }
    for (group, stratum, records) in census.ssr_sections() {
        sections.push((SectionKey::Ssr(group, stratum), encode_ssr(records)));
    }

    let mut head = Vec::new();
    head.extend_from_slice(MAGIC);
    put_u16(&mut head, VERSION);
    encode_header(&mut head, census, pileup);

    // The directory's own size depends only on its entries' keys, so it can be laid out at its
    // final size before the offsets are known.
    let mut directory = Vec::new();
    put_u32(&mut directory, sections.len() as u32);
    for (key, _) in &sections {
        encode_key(&mut directory, *key);
        put_u64(&mut directory, 0);
        put_u64(&mut directory, 0);
    }
    let first_section = (head.len() + directory.len()) as u64;

    // Now the same directory again, with each section placed end to end after it.
    directory.clear();
    put_u32(&mut directory, sections.len() as u32);
    let mut at = first_section;
    for (key, bytes) in &sections {
        encode_key(&mut directory, *key);
        put_u64(&mut directory, at);
        put_u64(&mut directory, bytes.len() as u64);
        at += bytes.len() as u64;
    }

    out.write_all(&head)?;
    out.write_all(&directory)?;
    for (_, bytes) in &sections {
        out.write_all(bytes)?;
    }
    Ok(())
}

fn encode_header(out: &mut Vec<u8>, census: &SampleCensusEvidence, pileup: Option<PileupIdentity>) {
    put_str(out, &census.sample);

    // **Who the read groups are, beside who the sample is**, and in the sample's own numbering:
    // entry `i` names the group its sections are keyed under as `i`. A cohort of censuses is
    // merged on these — every census numbers its groups from zero, because a walk sees one
    // sample, so the numbers collide by construction and only the names can tell two libraries
    // apart.
    let declared = census.declared_read_groups();
    put_u32(out, declared.len() as u32);
    for (id, group) in declared {
        put_u32(out, id.get());
        put_str(out, &group.declared_id);
        put_str(out, &group.library);
    }

    let terms = &census.terms;

    put_u16(out, terms.selection.fields().len() as u16);
    for (name, digest) in terms.selection.fields() {
        put_str(out, name);
        out.extend_from_slice(digest);
    }

    out.extend_from_slice(&terms.kept_loci.whole());
    let blocks = terms.kept_loci.blocks();
    put_u32(out, blocks.len() as u32);
    for block in blocks {
        put_u32(out, block.contig.get());
        put_u32(out, block.megabase);
        put_u64(out, block.digest);
    }

    let counts = terms.ssr_stratum_counts.iter_sorted();
    put_u32(out, counts.len() as u32);
    for ((period, repeats), loci) in counts {
        out.push(period);
        put_u64(out, repeats);
        put_u64(out, loci);
    }

    put_u32(out, terms.read_cap.0);
    out.extend_from_slice(&terms.depth_ladder.0);
    put_u32(out, terms.depth_cap.get());

    match pileup {
        Some(identity) => {
            out.push(1);
            out.extend_from_slice(&identity.header);
            put_u64(out, identity.records);
        }
        None => out.push(0),
    }
}

fn encode_key(out: &mut Vec<u8>, key: SectionKey) {
    match key {
        SectionKey::Generic(group) => {
            out.push(0);
            put_u32(out, group.get());
        }
        SectionKey::Ssr(group, stratum) => {
            out.push(1);
            put_u32(out, group.get());
            out.push(stratum.period);
            put_u64(out, stratum.reference_repeats);
        }
    }
}

fn encode_generic(records: &GenericEvidence) -> Vec<u8> {
    let mut out = Vec::new();
    put_u64(&mut out, records.depth().len() as u64);
    let bits = records.depth().as_bytes();
    put_u64(&mut out, bits.len() as u64);
    out.extend_from_slice(bits);
    put_u64(&mut out, records.non_reference().len() as u64);
    for entry in records.non_reference() {
        put_u32(&mut out, entry.index);
        out.push(entry.allele.code());
        out.push(entry.reads);
    }
    out
}

fn encode_ssr(records: &SsrEvidence) -> Vec<u8> {
    let mut out = Vec::new();
    put_u64(&mut out, records.len() as u64);
    for locus in 0..records.len() {
        for count in records.offsets(locus).counts() {
            put_u16(&mut out, *count);
        }
    }
    put_u64(&mut out, records.covering_not_crossing());
    let bits = records.walked_bits().as_bytes();
    put_u64(&mut out, bits.len() as u64);
    out.extend_from_slice(bits);
    put_u64(&mut out, records.bases_compared());
    put_u64(&mut out, records.guard().len() as u64);
    for entry in records.guard() {
        put_u32(&mut out, entry.locus);
        put_u32(&mut out, entry.length_difference as u32);
        put_u16(&mut out, entry.reads);
    }
    put_u64(&mut out, records.differences().len() as u64);
    for entry in records.differences() {
        put_u32(&mut out, entry.locus);
        put_u16(&mut out, entry.read);
        put_u16(&mut out, entry.offset as u16);
        out.push(entry.base.code());
    }
    out
}

fn put_u16(out: &mut Vec<u8>, value: u16) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn put_str(out: &mut Vec<u8>, value: &str) {
    put_u32(out, value.len() as u32);
    out.extend_from_slice(value.as_bytes());
}

// ---------------------------------------------------------------------
// Reading
// ---------------------------------------------------------------------

/// Read a whole census file.
///
/// **This is the parity oracle for the seeking reader**, not the shape a large run uses: it
/// decodes every section. What a run reading one stratum at a time will use is the directory
/// this returns alongside.
///
/// # Errors
///
/// [`CensusError::Malformed`] when the stream is not a census, is of another version, or ends
/// inside a value; [`CensusError::Io`] when the stream itself fails.
pub fn read_census(input: &mut impl Read) -> Result<CensusFile, CensusError> {
    let mut bytes = Vec::new();
    input.read_to_end(&mut bytes)?;
    decode_census(&bytes)
}

/// The same, from bytes already in hand.
pub fn decode_census(bytes: &[u8]) -> Result<CensusFile, CensusError> {
    let mut cursor = Cursor::new(bytes);
    if cursor.take(MAGIC.len())? != MAGIC {
        return Err(CensusError::Malformed);
    }
    if cursor.u16()? != VERSION {
        return Err(CensusError::Malformed);
    }
    let (sample, declared, terms, pileup) = decode_header(&mut cursor)?;
    let directory = decode_directory(&mut cursor)?;

    let mut sections = std::collections::BTreeMap::new();
    for (key, extent) in directory {
        let at = usize::try_from(extent.offset()).map_err(|_| CensusError::Malformed)?;
        let len = usize::try_from(extent.len()).map_err(|_| CensusError::Malformed)?;
        let slice = bytes
            .get(at..at.checked_add(len).ok_or(CensusError::Malformed)?)
            .ok_or(CensusError::Malformed)?;
        sections.insert(key, decode_section(key, slice)?);
    }

    Ok(CensusFile {
        census: SampleCensusEvidence::resident(sample, terms, declared, sections),
        pileup,
    })
}

/// Open a census file **without decoding a section**: the header is checked and the directory is
/// read, and nothing else is touched until a scoped call asks for something.
///
/// **This is what a run over a cohort too large to hold uses.** What it costs is the header and
/// the directory — a few hundred bytes plus sixteen a section — where reading the file whole
/// costs the file.
///
/// # Errors
///
/// [`CensusError::Io`] when the file will not open or read, [`CensusError::Malformed`] when it is
/// not a census this build reads.
pub fn open_census(
    path: &Path,
) -> Result<(SampleCensusEvidence, Option<PileupIdentity>), CensusError> {
    // The header and the directory sit at the front, so only their bytes are read. How many that
    // is is not known until they are decoded, so the front of the file is read in one go and the
    // rest is never touched.
    let mut file = std::fs::File::open(path)?;
    let mut head = vec![0_u8; HEAD_READ_BYTES];
    let filled = read_as_much_as_there_is(&mut file, &mut head)?;
    head.truncate(filled);

    let mut cursor = Cursor::new(&head);
    if cursor.take(MAGIC.len())? != MAGIC || cursor.u16()? != VERSION {
        return Err(CensusError::Malformed);
    }
    let (sample, declared, terms, pileup) = decode_header(&mut cursor)?;
    let directory = decode_directory(&mut cursor)?;

    Ok((
        SampleCensusEvidence::backed(
            sample,
            terms,
            declared,
            path.to_path_buf(),
            directory.into_iter().collect(),
        ),
        pileup,
    ))
}

/// How much of a census file's front is read to find its header and its directory.
///
/// **One read rather than a walk of growing reads.** A census of tomato's 141 strata over a
/// handful of read groups has a directory of a few thousand bytes and a header of a few hundred;
/// a megabyte covers a directory of about 40,000 sections, and reading it costs one seekless read
/// of a file that is megabytes long anyway. A file shorter than this is read whole.
const HEAD_READ_BYTES: usize = 1 << 20;

/// Fill as much of `into` as the stream has, and say how much that was — a short file is not an
/// error here, since the header may be all there is to read.
fn read_as_much_as_there_is(from: &mut impl Read, into: &mut [u8]) -> Result<usize, CensusError> {
    let mut filled = 0;
    while filled < into.len() {
        match from.read(&mut into[filled..])? {
            0 => break,
            n => filled += n,
        }
    }
    Ok(filled)
}

/// Where each section sits, without decoding one — what the seeking reader will open with.
///
/// # Errors
///
/// As [`decode_census`].
pub fn decode_directory_of(bytes: &[u8]) -> Result<Vec<(SectionKey, ByteExtent)>, CensusError> {
    let mut cursor = Cursor::new(bytes);
    if cursor.take(MAGIC.len())? != MAGIC || cursor.u16()? != VERSION {
        return Err(CensusError::Malformed);
    }
    decode_header(&mut cursor)?;
    decode_directory(&mut cursor)
}

type Header = (
    String,
    BTreeMap<ReadGroupId, NamedReadGroup>,
    RecordingTerms,
    Option<PileupIdentity>,
);

fn decode_header(cursor: &mut Cursor<'_>) -> Result<Header, CensusError> {
    let sample = cursor.string()?;

    let declared_count = cursor.u32()? as usize;
    let mut declared = BTreeMap::new();
    for _ in 0..declared_count {
        let id = ReadGroupId(cursor.u32()?);
        let named = NamedReadGroup {
            declared_id: cursor.string()?,
            library: cursor.string()?,
        };
        // **Two entries under one identifier is a malformed file, not a last-one-wins.** They
        // would name one section's read group two ways, and only one of the two could reach a
        // cohort's merge.
        if declared.insert(id, named).is_some() {
            return Err(CensusError::Malformed);
        }
    }

    let field_count = cursor.u16()? as usize;
    let mut fields = Vec::with_capacity(field_count);
    for _ in 0..field_count {
        let name = cursor.string()?;
        let digest = cursor.digest()?;
        // **The names are this build's own or the file is not comparable.** A table whose fields
        // are not these would be compared entry by entry against fields that mean something
        // else, which is a silent wrong answer rather than a refusal.
        let known = super::census::SELECTION_FIELDS
            .into_iter()
            .find(|known| *known == name)
            .ok_or(CensusError::Malformed)?;
        fields.push((known, digest));
    }
    if fields.len() != super::census::SELECTION_FIELDS.len()
        || fields
            .iter()
            .zip(super::census::SELECTION_FIELDS)
            .any(|((read, _), mine)| *read != mine)
    {
        return Err(CensusError::Malformed);
    }

    let whole = cursor.digest()?;
    let block_count = cursor.u32()? as usize;
    let mut blocks = Vec::with_capacity(block_count.min(1 << 20));
    for _ in 0..block_count {
        blocks.push(BlockDigest {
            contig: ContigId(cursor.u32()?),
            megabase: cursor.u32()?,
            digest: cursor.u64()?,
        });
    }

    let stratum_count = cursor.u32()? as usize;
    let mut counts = Vec::with_capacity(stratum_count.min(1 << 20));
    for _ in 0..stratum_count {
        let period = cursor.u8()?;
        let repeats = cursor.u64()?;
        counts.push(((period, repeats), cursor.u64()?));
    }

    let read_cap = ReadCap(cursor.u32()?);
    let depth_ladder = DepthLadderDigest(cursor.digest()?);
    let depth_cap = cursor.u32()?;
    if depth_cap > u32::from(u8::MAX) {
        return Err(CensusError::Malformed);
    }

    let pileup = match cursor.u8()? {
        0 => None,
        1 => Some(PileupIdentity {
            header: cursor.digest()?,
            records: cursor.u64()?,
        }),
        _ => return Err(CensusError::Malformed),
    };

    Ok((
        sample,
        declared,
        RecordingTerms {
            selection: SelectionTermsDigest::from_fields(fields),
            kept_loci: CensusLociDigest::from_parts(whole, blocks),
            ssr_stratum_counts: StratumCounts::from_counted(counts),
            read_cap,
            depth_ladder,
            depth_cap: DepthCap::new(depth_cap),
        },
        pileup,
    ))
}

fn decode_directory(cursor: &mut Cursor<'_>) -> Result<Vec<(SectionKey, ByteExtent)>, CensusError> {
    let count = cursor.u32()? as usize;
    let mut entries = Vec::with_capacity(count.min(1 << 20));
    for _ in 0..count {
        let key = decode_key(cursor)?;
        let offset = cursor.u64()?;
        let len = cursor.u64()?;
        entries.push((key, ByteExtent::new(offset, len)));
    }
    // **No two sections may share a byte**, since a call reading one seeks to its own offset and
    // takes its own length; overlapping extents would hand a caller another section's bytes.
    for (index, (_, mine)) in entries.iter().enumerate() {
        if entries[index + 1..]
            .iter()
            .any(|(_, theirs)| mine.overlaps(*theirs))
        {
            return Err(CensusError::Malformed);
        }
    }
    Ok(entries)
}

fn decode_key(cursor: &mut Cursor<'_>) -> Result<SectionKey, CensusError> {
    match cursor.u8()? {
        0 => Ok(SectionKey::Generic(ReadGroupId(cursor.u32()?))),
        1 => {
            let group = ReadGroupId(cursor.u32()?);
            let period = cursor.u8()?;
            Ok(SectionKey::Ssr(
                group,
                Stratum {
                    period,
                    reference_repeats: cursor.u64()?,
                },
            ))
        }
        _ => Err(CensusError::Malformed),
    }
}

/// One section, from exactly its own bytes — **what a seeking reader calls**, once it has read
/// the extent the directory gave it.
///
/// # Errors
///
/// [`CensusError::Malformed`] when the bytes are not a section of the kind the key names, or
/// when they end inside a value.
pub(super) fn decode_section(key: SectionKey, bytes: &[u8]) -> Result<Section, CensusError> {
    let mut cursor = Cursor::new(bytes);
    let section = match key {
        SectionKey::Generic(_) => Section::Generic(decode_generic(&mut cursor)?),
        SectionKey::Ssr(_, _) => Section::Ssr(decode_ssr(&mut cursor)?),
    };
    // **A section that does not use all its own bytes is not this section.** The extent came
    // from the directory, so bytes left over mean the two disagree about what is stored there.
    if cursor.at != bytes.len() {
        return Err(CensusError::Malformed);
    }
    Ok(section)
}

fn decode_generic(cursor: &mut Cursor<'_>) -> Result<GenericEvidence, CensusError> {
    let positions = usize::try_from(cursor.u64()?).map_err(|_| CensusError::Malformed)?;
    let bits = cursor.bytes()?.to_vec();
    // **The width comes from the encoding and is not written again here.** It stood as a bare
    // `5` until the depth code widened to eight bits, and a second statement of a number the
    // writer takes from `DEPTH_CODE_BITS` is exactly the bug this check exists to catch.
    if bits.len()
        != positions
            .saturating_mul(DEPTH_CODE_BITS as usize)
            .div_ceil(8)
    {
        return Err(CensusError::Malformed);
    }
    let entries = usize::try_from(cursor.u64()?).map_err(|_| CensusError::Malformed)?;
    let mut non_reference = Vec::with_capacity(entries.min(1 << 22));
    for _ in 0..entries {
        let index = cursor.u32()?;
        let allele = ObservedAllele::of_code(cursor.u8()?).ok_or(CensusError::Malformed)?;
        non_reference.push(AlleleObservation {
            index,
            allele,
            reads: cursor.u8()?,
        });
    }
    // The sparse list arrives in position order and names positions the dense array has, or
    // `from_parts` panics — a file is not a caller, so the check is a refusal here.
    if non_reference
        .windows(2)
        .any(|pair| pair[0].index > pair[1].index)
        || non_reference
            .last()
            .is_some_and(|entry| entry.index as usize >= positions)
    {
        return Err(CensusError::Malformed);
    }
    Ok(GenericEvidence::from_parts(
        PackedDepthCodes::from_bytes(bits, positions),
        non_reference,
    ))
}

fn decode_ssr(cursor: &mut Cursor<'_>) -> Result<SsrEvidence, CensusError> {
    let loci = usize::try_from(cursor.u64()?).map_err(|_| CensusError::Malformed)?;
    let mut offsets = Vec::with_capacity(loci.min(1 << 22));
    for _ in 0..loci {
        let mut counts = [0_u16; OFFSET_BUCKETS];
        for count in &mut counts {
            *count = cursor.u16()?;
        }
        offsets.push(OffsetCounts::from_counts(counts));
    }
    let covering_not_crossing = cursor.u64()?;
    let bits = cursor.bytes()?.to_vec();
    if bits.len() != loci.div_ceil(8) {
        return Err(CensusError::Malformed);
    }
    let bases_compared = cursor.u64()?;

    let guard_count = usize::try_from(cursor.u64()?).map_err(|_| CensusError::Malformed)?;
    let mut guard = Vec::with_capacity(guard_count.min(1 << 22));
    for _ in 0..guard_count {
        guard.push(GuardObservation {
            locus: cursor.u32()?,
            length_difference: cursor.u32()? as i32,
            reads: cursor.u16()?,
        });
    }

    let difference_count = usize::try_from(cursor.u64()?).map_err(|_| CensusError::Malformed)?;
    let mut differences = Vec::with_capacity(difference_count.min(1 << 22));
    for _ in 0..difference_count {
        differences.push(TractDifference {
            locus: cursor.u32()?,
            read: cursor.u16()?,
            offset: cursor.u16()? as i16,
            base: ObservedAllele::of_code(cursor.u8()?).ok_or(CensusError::Malformed)?,
        });
    }

    let out_of_range = loci as u32;
    if guard.iter().any(|entry| entry.locus >= out_of_range)
        || differences.iter().any(|entry| entry.locus >= out_of_range)
    {
        return Err(CensusError::Malformed);
    }

    Ok(SsrEvidence::from_parts(
        offsets,
        covering_not_crossing,
        WalkedBits::from_bytes(bits, loci),
        bases_compared,
        guard,
        differences,
    ))
}

/// A position in a byte slice that refuses to read past the end.
///
/// **Every short read is [`CensusError::Malformed`] and never a panic**, because the bytes come
/// from a file: a truncated census is a thing that happens, and a decoder that panics on one
/// takes the whole run with it.
struct Cursor<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Cursor<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, at: 0 }
    }

    fn take(&mut self, len: usize) -> Result<&'a [u8], CensusError> {
        let end = self.at.checked_add(len).ok_or(CensusError::Malformed)?;
        let slice = self.bytes.get(self.at..end).ok_or(CensusError::Malformed)?;
        self.at = end;
        Ok(slice)
    }

    fn u8(&mut self) -> Result<u8, CensusError> {
        Ok(self.take(1)?[0])
    }

    fn u16(&mut self) -> Result<u16, CensusError> {
        Ok(u16::from_le_bytes(
            self.take(2)?.try_into().expect("two bytes"),
        ))
    }

    fn u32(&mut self) -> Result<u32, CensusError> {
        Ok(u32::from_le_bytes(
            self.take(4)?.try_into().expect("four bytes"),
        ))
    }

    fn u64(&mut self) -> Result<u64, CensusError> {
        Ok(u64::from_le_bytes(
            self.take(8)?.try_into().expect("eight bytes"),
        ))
    }

    fn digest(&mut self) -> Result<[u8; 16], CensusError> {
        Ok(self.take(16)?.try_into().expect("sixteen bytes"))
    }

    /// A `u64` length and that many bytes.
    fn bytes(&mut self) -> Result<&'a [u8], CensusError> {
        let len = usize::try_from(self.u64()?).map_err(|_| CensusError::Malformed)?;
        self.take(len)
    }

    fn string(&mut self) -> Result<String, CensusError> {
        let len = self.u32()? as usize;
        std::str::from_utf8(self.take(len)?)
            .map(str::to_string)
            .map_err(|_| CensusError::Malformed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;

    use std::collections::BTreeMap;

    use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
    use crate::ng::parameter_estimation::joint::census::{
        CohortCensusEvidence, DepthCode, SsrLocusState,
    };
    use crate::ng::parameter_estimation::joint::loci::CensusLociDigester;
    use crate::ng::types::{GenomePosition, Position};

    const AT_SIX_REPEATS: Stratum = Stratum {
        period: 2,
        reference_repeats: 6,
    };
    const ATG_FOUR_REPEATS: Stratum = Stratum {
        period: 3,
        reference_repeats: 4,
    };

    fn selection_terms() -> crate::ng::parameter_estimation::joint::loci::SelectionTerms {
        use crate::ng::parameter_estimation::joint::loci::{
            CatalogBuildSettings, ReferenceDigest, RegionSetDigest, SelectionTerms,
        };
        use crate::ng::repeat_catalog::StrRepeatCriteria;
        use crate::ng::tandem_repeat::ScanParams;
        SelectionTerms {
            seed: 42,
            reference: ReferenceDigest([7; 16]),
            analysed_regions: RegionSetDigest([9; 16]),
            catalog_built_under: CatalogBuildSettings {
                criteria: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "0.1.0".to_string(),
            },
            ssr_criteria: StrRepeatCriteria::default(),
            generic_target: 2_000_000,
            ssr_cap: 1_000,
        }
    }

    /// Terms whose kept-loci digest witnesses two real positions in two different megabases, so
    /// the block list is not empty and its round trip is asserted rather than assumed.
    fn terms() -> RecordingTerms {
        let mut digester = CensusLociDigester::new();
        for (index, position) in [7_u64, 4_000_003].into_iter().enumerate() {
            digester.observe(
                index,
                GenomePosition {
                    contig: ContigId(0),
                    position: Position(position),
                },
            );
        }
        let mut counts = StratumCounts::default();
        counts.count(2, 6);
        counts.count(2, 6);
        counts.count(3, 4);
        RecordingTerms {
            selection: SelectionTermsDigest::of(&selection_terms()),
            kept_loci: digester.finish(),
            ssr_stratum_counts: counts,
            read_cap: ReadCap(1_000),
            depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
            depth_cap: DepthCap::new(124),
        }
    }

    /// **Every corner spec §7.1 lists, in one sample.** Five ordinary positions — never walked,
    /// walked at zero depth, reads with none non-reference, one non-reference allele, and two at
    /// one position — plus two strata of tracts carrying a saturating offset, a guard entry, a
    /// difference, and a locus the walk never reached.
    fn every_corner() -> SampleCensusEvidence {
        let edges = DepthBinEdges::for_census();
        let mut depth = PackedDepthCodes::never_walked(5);
        depth.set(1, DepthCode::Binned(edges.bin_for(0)));
        depth.set(2, DepthCode::Binned(edges.bin_for(6)));
        depth.set(3, DepthCode::Binned(edges.bin_for(124)));
        depth.set(4, DepthCode::Binned(edges.bin_for(300)));
        let generic = GenericEvidence::from_parts(
            depth,
            vec![
                AlleleObservation {
                    index: 3,
                    allele: ObservedAllele::G,
                    reads: u8::MAX,
                },
                AlleleObservation {
                    index: 4,
                    allele: ObservedAllele::C,
                    reads: 1,
                },
                AlleleObservation {
                    index: 4,
                    allele: ObservedAllele::Other,
                    reads: 3,
                },
            ],
        );

        let mut offsets = vec![OffsetCounts::default(); 2];
        offsets[0].add(0, 5);
        // Both saturate into an end bucket at any recorded range, which is the corner this
        // fixture is here to carry through the file and back.
        offsets[0].add(-RECORDED_OFFSET_RANGE - 1, 2);
        offsets[0].add(RECORDED_OFFSET_RANGE + 3, 3);
        let mut walked = WalkedBits::none_of(2);
        walked.set(0); // locus 1 is never walked, which is the state with no other field
        let tracts = SsrEvidence::from_parts(
            offsets,
            17,
            walked,
            60,
            vec![GuardObservation {
                locus: 0,
                length_difference: -3,
                reads: 2,
            }],
            vec![TractDifference {
                locus: 0,
                read: 299,
                offset: 8,
                base: ObservedAllele::T,
            }],
        );

        SampleCensusEvidence::resident(
            "corners".to_string(),
            terms(),
            NamedReadGroup::drawn_for("corners", [ReadGroupId(0)]),
            BTreeMap::from([
                (
                    SectionKey::Generic(ReadGroupId(0)),
                    Section::Generic(generic),
                ),
                (
                    SectionKey::Generic(ReadGroupId(1)),
                    Section::Generic(GenericEvidence::never_walked(5)),
                ),
                (
                    SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS),
                    Section::Ssr(tracts),
                ),
                (
                    SectionKey::Ssr(ReadGroupId(0), ATG_FOUR_REPEATS),
                    Section::Ssr(SsrEvidence::never_walked(1)),
                ),
            ]),
        )
    }

    fn round_trip(census: &SampleCensusEvidence, pileup: Option<PileupIdentity>) -> CensusFile {
        let mut bytes = Vec::new();
        write_census(census, pileup, &mut bytes).expect("a vector accepts every write");
        decode_census(&bytes).expect("what this build wrote, this build reads")
    }

    /// **The assertion the whole step rests on**, and the reason B1 is its own commit: a codec
    /// that reads a field at the wrong offset produces a plausible number rather than a crash.
    #[test]
    fn every_corner_state_survives_a_round_trip() {
        let census = every_corner();
        let read = round_trip(&census, None);
        assert_eq!(read.census, census, "the whole sample, field for field");
        assert_eq!(read.pileup, None);
    }

    /// The same value written twice is the same bytes — what §7.12's byte-for-byte comparison
    /// between the two builders will rest on, and what a directory whose offsets depended on a
    /// map's iteration order would break.
    #[test]
    fn writing_the_same_census_twice_gives_the_same_bytes() {
        let census = every_corner();
        let (mut first, mut second) = (Vec::new(), Vec::new());
        write_census(&census, None, &mut first).expect("a vector accepts every write");
        write_census(&census, None, &mut second).expect("a vector accepts every write");
        assert_eq!(first, second);
    }

    /// The three states an ordinary position can be in, read back off the file rather than off
    /// the value that was written — the distinction that has no field of its own.
    #[test]
    fn the_three_states_at_an_ordinary_position_come_back_apart() {
        let read = round_trip(&every_corner(), None);
        let mut census = read.census;
        let groups = census.read_groups();
        census
            .with_generic(&groups, |sections| {
                let edges = DepthBinEdges::for_census();
                let records = sections[0];
                assert_eq!(records.at(0).0, DepthCode::NeverWalked, "a bug");
                assert_eq!(
                    records.at(1).0,
                    DepthCode::Binned(edges.bin_for(0)),
                    "walked and empty is data, not a bug"
                );
                assert_eq!(records.at(2).0, DepthCode::Binned(edges.bin_for(6)));
                assert!(
                    records.at(2).1.is_empty(),
                    "six reads and no sparse entry means six reads on the reference base"
                );
                let (depth, alleles) = records.at(4);
                assert_eq!(depth, DepthCode::Binned(edges.bin_for(300)));
                assert_eq!(alleles.len(), 2, "a multi-allelic position keeps both");
            })
            // **Not discarded.** If lending the sections failed, the closure above would never
            // run and every assertion in it would be skipped — a test that passes by not
            // looking.
            .expect("a decoded census is resident and has no file to fail on");
    }

    /// The tract half's own corners: the saturating end buckets, the guard, the difference and
    /// its read number, the two per-stratum counts, and the locus the walk never reached.
    #[test]
    fn a_tracts_offsets_guard_and_difference_come_back_unchanged() {
        let read = round_trip(&every_corner(), None);
        let mut census = read.census;
        census
            .with_strata(ReadGroupId(0), &[AT_SIX_REPEATS], |sections| {
                let tracts = sections[0];
                assert_eq!(tracts.offsets(0).at(0), 5);
                assert_eq!(
                    tracts.offsets(0).at(-RECORDED_OFFSET_RANGE),
                    2,
                    "one past the short end saturates into it"
                );
                assert_eq!(
                    tracts.offsets(0).at(RECORDED_OFFSET_RANGE),
                    3,
                    "and three past the long end into that one"
                );
                assert_eq!(tracts.covering_not_crossing(), 17);
                assert_eq!(tracts.bases_compared(), 60);
                assert_eq!(tracts.guard()[0].length_difference, -3);
                assert_eq!(tracts.differences()[0].read, 299);
                assert_eq!(tracts.differences()[0].offset, 8);
                assert_eq!(tracts.differences()[0].base, ObservedAllele::T);
                assert_eq!(tracts.state(0), SsrLocusState::Crossed);
                assert_eq!(
                    tracts.state(1),
                    SsrLocusState::NeverWalked,
                    "the one state every other field reads the same as walked-and-empty"
                );
            })
            // **Not discarded**, for the reason the test above gives: a lending failure would
            // skip every assertion in the closure and the test would still pass.
            .expect("a decoded census is resident and has no file to fail on");
    }

    /// The directory says where a section is without decoding one — which is what the seeking
    /// reader will open with, and what a counting reader will measure against.
    #[test]
    fn the_directory_places_every_section_end_to_end_and_none_overlap() {
        let census = every_corner();
        let mut bytes = Vec::new();
        write_census(&census, None, &mut bytes).expect("a vector accepts every write");

        let directory = decode_directory_of(&bytes).expect("this build's own file");
        assert_eq!(
            directory.iter().map(|(key, _)| *key).collect::<Vec<_>>(),
            vec![
                SectionKey::Generic(ReadGroupId(0)),
                SectionKey::Generic(ReadGroupId(1)),
                SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS),
                SectionKey::Ssr(ReadGroupId(0), ATG_FOUR_REPEATS),
            ],
            "the enumeration order is the key's own"
        );
        for pair in directory.windows(2) {
            assert!(!pair[0].1.overlaps(pair[1].1));
            assert_eq!(
                pair[0].1.offset() + pair[0].1.len(),
                pair[1].1.offset(),
                "sections abut, so the file is its header plus its sections and nothing else"
            );
        }
        let last = directory.last().expect("four sections");
        assert_eq!(
            last.1.offset() + last.1.len(),
            bytes.len() as u64,
            "and the last one ends at the end of the file"
        );
    }

    /// The pileup a census was built from travels with it, so B3's staleness check has something
    /// to compare. **A census built during a walk that wrote no pileup carries none**, and the
    /// two cases have to be told apart rather than one standing for the other.
    #[test]
    fn the_pileup_a_census_was_built_from_survives_and_its_absence_is_not_a_zero() {
        let identity = PileupIdentity {
            header: [3; 16],
            records: 12_345,
        };
        assert_eq!(
            round_trip(&every_corner(), Some(identity)).pileup,
            Some(identity)
        );
        assert_eq!(round_trip(&every_corner(), None).pileup, None);
    }

    /// **A file this build did not write is refused rather than decoded.** Every one of these
    /// would otherwise produce a plausible census: a truncated file reads a length off the end
    /// of a buffer, and a wrong version reads this build's fields at another build's offsets.
    #[test]
    fn a_stream_that_is_not_this_builds_census_is_refused() {
        let census = every_corner();
        let mut bytes = Vec::new();
        write_census(&census, None, &mut bytes).expect("a vector accepts every write");

        assert!(matches!(decode_census(&[]), Err(CensusError::Malformed)));
        let mut wrong_magic = bytes.clone();
        wrong_magic[0] = b'X';
        assert!(matches!(
            decode_census(&wrong_magic),
            Err(CensusError::Malformed)
        ));
        let mut wrong_version = bytes.clone();
        wrong_version[8] = VERSION as u8 + 1;
        assert!(matches!(
            decode_census(&wrong_version),
            Err(CensusError::Malformed)
        ));
        // Every truncation, not one: a length read off the end of the buffer must be a refusal
        // wherever it falls, and a decoder that checks only at the top of the file is not one.
        for cut in 0..bytes.len() {
            assert!(
                matches!(decode_census(&bytes[..cut]), Err(CensusError::Malformed)),
                "a census cut to {cut} of {} bytes decoded",
                bytes.len()
            );
        }
    }

    /// **The oracle for the file-backed reader: it must answer what the resident value answers.**
    /// Every section of the corner fixture, asked for through the same scoped calls, off a file
    /// on disk against the value in memory.
    #[test]
    fn a_census_read_from_a_file_answers_what_the_one_in_memory_answers() {
        let dir = tempfile::tempdir().expect("a scratch directory");
        let path = dir.path().join("corners.census");
        let mut resident = every_corner();
        write_census(
            &resident,
            None,
            &mut std::fs::File::create(&path).expect("a new file"),
        )
        .expect("a file accepts every write");

        let (mut backed, pileup) = open_census(&path).expect("this build's own file");
        assert_eq!(pileup, None);
        assert_eq!(backed.sample, resident.sample);
        assert_eq!(
            backed.terms, resident.terms,
            "the terms are read at the door"
        );
        assert_eq!(backed.read_groups(), resident.read_groups());
        assert_eq!(backed.strata(), resident.strata());

        let groups = resident.read_groups();
        let from_memory = resident
            .with_generic(&groups, |sections| {
                sections.iter().map(|g| (*g).clone()).collect::<Vec<_>>()
            })
            .expect("a resident census has no file to fail on");
        let from_file = backed
            .with_generic(&groups, |sections| {
                sections.iter().map(|g| (*g).clone()).collect::<Vec<_>>()
            })
            .expect("this build's own file");
        assert_eq!(from_file, from_memory, "every ordinary-position section");

        // Read group 1 recorded ordinary positions and no tracts in this fixture, so the tract
        // half is asked of the group that has one — a real census gives every declared group
        // every stratum, and asking for a section that was never recorded is a panic by
        // contract.
        let strata = resident.strata();
        // One read group, because `every_corner` declares one — a loop over a single value
        // read as if the count were open, and it is not.
        let group = ReadGroupId(0);
        {
            let from_memory = resident
                .with_strata(group, &strata, |sections| {
                    sections.iter().map(|s| (*s).clone()).collect::<Vec<_>>()
                })
                .expect("a resident census has no file to fail on");
            let from_file = backed
                .with_strata(group, &strata, |sections| {
                    sections.iter().map(|s| (*s).clone()).collect::<Vec<_>>()
                })
                .expect("this build's own file");
            assert_eq!(from_file, from_memory, "every tract section of {group:?}");
        }
    }

    /// **The bytes a call reads are the section's own — spec §7.15's second half.** The first
    /// half is the values, which the test above pins; this is the one that tells a reader
    /// seeking to one section from one that decodes the file and hands back a slice, because
    /// both give the same values and only one of them delivers the memory.
    #[test]
    fn asking_for_one_section_reads_that_section_and_no_other_byte() {
        let dir = tempfile::tempdir().expect("a scratch directory");
        let path = dir.path().join("counted.census");
        let census = every_corner();
        write_census(
            &census,
            None,
            &mut std::fs::File::create(&path).expect("a new file"),
        )
        .expect("a file accepts every write");
        let whole_file = std::fs::metadata(&path).expect("the file exists").len();

        let bytes = std::fs::read(&path).expect("the file reads");
        let directory = decode_directory_of(&bytes).expect("this build's own file");
        let extent_of = |wanted: SectionKey| {
            directory
                .iter()
                .find(|(key, _)| *key == wanted)
                .map(|(_, extent)| extent.len())
                .expect("the section is in the directory")
        };

        // Opening reads the head of the file, which is not a section — so the count starts
        // after it, where the sections do.
        let (mut backed, _) = open_census(&path).expect("this build's own file");
        reset_bytes_read();
        backed
            .with_strata(ReadGroupId(0), &[AT_SIX_REPEATS], |sections| {
                sections[0].len()
            })
            .expect("this build's own file");
        let one_stratum = extent_of(SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS));
        assert_eq!(
            bytes_read(),
            one_stratum,
            "one stratum's read is that stratum's bytes and nothing else"
        );
        assert!(
            one_stratum < whole_file,
            "and that is less than the file: {one_stratum} of {whole_file}"
        );

        // A band of two reads exactly the two, and a second call for the same band reads them
        // again — nothing was retained between the calls.
        reset_bytes_read();
        backed
            .with_generic(&[ReadGroupId(0), ReadGroupId(1)], |sections| sections.len())
            .expect("this build's own file");
        let both = extent_of(SectionKey::Generic(ReadGroupId(0)))
            + extent_of(SectionKey::Generic(ReadGroupId(1)));
        assert_eq!(bytes_read(), both);
        backed
            .with_generic(&[ReadGroupId(0), ReadGroupId(1)], |sections| sections.len())
            .expect("this build's own file");
        assert_eq!(
            bytes_read(),
            2 * both,
            "a second call reads them again, because the first kept nothing"
        );
    }

    /// **What milestone B is for: the same cohort fitted from memory and from files gives the
    /// same parameters** (spec §7.15's first half, at the grain a user sees it). Three samples
    /// over 400 positions, drawn so that some carry a non-reference allele and some do not.
    #[test]
    fn a_cohort_fitted_from_files_gives_the_parameters_it_gives_from_memory() {
        use crate::ng::parameter_estimation::joint::fit::{JointFitConfig, fit_jointly};

        let dir = tempfile::tempdir().expect("a scratch directory");
        let samples = drawn_cohort();
        let mut from_memory =
            CohortCensusEvidence::new(samples.clone()).expect("every sample recorded one way");

        let mut backed = Vec::new();
        for sample in &samples {
            let path = dir.path().join(format!("{}.census", sample.sample));
            write_census(
                sample,
                None,
                &mut std::fs::File::create(&path).expect("a new file"),
            )
            .expect("a file accepts every write");
            backed.push(open_census(&path).expect("this build's own file").0);
        }
        let mut from_files =
            CohortCensusEvidence::new(backed).expect("every sample recorded one way");

        let config = JointFitConfig::default();
        let memory = fit_jointly(&mut from_memory, &config).expect("a drawn cohort pools");
        let files = fit_jointly(&mut from_files, &config).expect("a drawn cohort pools");

        assert_eq!(
            files.log_likelihood, memory.log_likelihood,
            "the likelihood is the same number, not a near one"
        );
        assert_eq!(files.density.value.a, memory.density.value.a);
        assert_eq!(files.density.value.b, memory.density.value.b);
        assert_eq!(files.noisy_share, memory.noisy_share);
        for sample in &samples {
            assert_eq!(
                files.rates[&sample.sample].value.heterozygous,
                memory.rates[&sample.sample].value.heterozygous,
                "{}'s heterozygosity",
                sample.sample
            );
            assert_eq!(
                files.hom_excess[&sample.sample].value.get(),
                memory.hom_excess[&sample.sample].value.get()
            );
        }
        for group in memory.noise.keys() {
            assert_eq!(
                files.noise[group].value.clean,
                memory.noise[group].value.clean
            );
            assert_eq!(
                files.noise[group].value.noisy,
                memory.noise[group].value.noisy
            );
        }
    }

    /// Three samples over 400 ordinary positions, drawn from one reproducible stream: most
    /// positions quiet, one in twenty carrying a non-reference allele in some samples.
    fn drawn_cohort() -> Vec<SampleCensusEvidence> {
        let edges = DepthBinEdges::for_census();
        let terms = terms();
        let mut state = 0x9E37_79B9_7F4A_7C15_u64;
        let mut next = move || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        (0..3)
            .map(|s| {
                let mut depth = PackedDepthCodes::never_walked(400);
                let mut sparse = Vec::new();
                for index in 0..400 {
                    let reads = 2 + next() % 6;
                    depth.set(index, DepthCode::Binned(edges.bin_for(reads)));
                    if next() % 20 == 0 {
                        sparse.push(AlleleObservation {
                            index: index as u32,
                            allele: ObservedAllele::C,
                            reads: (1 + next() % reads) as u8,
                        });
                    }
                }
                SampleCensusEvidence::resident(
                    format!("s{s}"),
                    terms.clone(),
                    NamedReadGroup::drawn_for(&format!("s{s}"), [ReadGroupId(s)]),
                    BTreeMap::from([(
                        // **One read group a sample**: a library is one plant's DNA
                        // preparation, so a cohort's samples never share one, and the cohort's
                        // door refuses a set that does.
                        SectionKey::Generic(ReadGroupId(s)),
                        Section::Generic(GenericEvidence::from_parts(depth, sparse)),
                    )]),
                )
            })
            .collect()
    }

    /// A resident census reads no bytes at all, which is the other half of the same property.
    #[test]
    fn a_resident_census_reads_nothing() {
        let mut census = every_corner();
        reset_bytes_read();
        let groups = census.read_groups();
        census
            .with_generic(&groups, |sections| sections.len())
            .expect("a resident census has no file to fail on");
        assert_eq!(bytes_read(), 0);
    }

    /// **A call asks for one section and reads one section's bytes.** The band here is one
    /// stratum of two, and what it must not do is decode the other — which a reader that read
    /// the file and handed back a slice would.
    #[test]
    fn a_call_for_one_stratum_reads_that_stratum_and_leaves_the_rest_alone() {
        let dir = tempfile::tempdir().expect("a scratch directory");
        let path = dir.path().join("one.census");
        let census = every_corner();
        write_census(
            &census,
            None,
            &mut std::fs::File::create(&path).expect("a new file"),
        )
        .expect("a file accepts every write");
        let file_len = std::fs::metadata(&path).expect("the file exists").len();

        let (mut backed, _) = open_census(&path).expect("this build's own file");
        let asked = backed
            .with_strata(ReadGroupId(0), &[AT_SIX_REPEATS], |sections| {
                assert_eq!(sections.len(), 1, "one stratum was asked for");
                sections[0].len()
            })
            .expect("this build's own file");
        assert_eq!(asked, 2, "the two tracts that stratum holds");

        // The section's own extent, which is what the call above may read and no more.
        let bytes = std::fs::read(&path).expect("the file reads");
        let directory = decode_directory_of(&bytes).expect("this build's own file");
        let extent = directory
            .iter()
            .find(|(key, _)| *key == SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS))
            .map(|(_, extent)| *extent)
            .expect("the stratum is in the directory");
        assert!(
            extent.len() < file_len,
            "a section is a part of the file, not the whole of it"
        );
    }

    // ---- the pileup a census names, and what a run does about it ---------------------

    /// **Spec §7.13.** A census built from one pileup and fitted against another must not be
    /// used as it stands: it is rebuilt where the pileup is there, and refused where it is not.
    /// The two cases are told apart by *which* value differs, so a run says why.
    #[test]
    fn a_census_naming_another_pileup_is_rebuilt_where_it_can_be_and_refused_where_it_cannot() {
        let built_from = PileupIdentity::of_header(b"reference=A regions=1-100 filters=q20", 4_000);
        assert_eq!(
            freshness(Some(built_from), Some(built_from)),
            Freshness::Fresh
        );

        // The cheapest thing to change, and the likeliest to change by accident: a different
        // analysed-region set, which is part of the header.
        let other_regions =
            PileupIdentity::of_header(b"reference=A regions=1-200 filters=q20", 4_000);
        assert_ne!(built_from.header, other_regions.header);
        assert_eq!(
            freshness(Some(built_from), Some(other_regions)),
            Freshness::Rebuild("the pileup's header")
        );
        assert_eq!(
            freshness(Some(built_from), None),
            Freshness::Refused("the pileup it was built from, which is not here")
        );

        // The same pileup with records added under it — the header is unchanged and the count
        // is not, which is exactly what the count is carried for.
        let grown = PileupIdentity {
            records: built_from.records + 1,
            ..built_from
        };
        assert_eq!(
            freshness(Some(built_from), Some(grown)),
            Freshness::Rebuild("the pileup's record count")
        );
    }

    /// **The check must not key on a modification time** (spec §6.1, §7.13), and this is what
    /// says it does not: the census file is touched, re-opened, and names the same pileup — so
    /// a run that copied its files, or one whose file system rewrote a timestamp, gets the same
    /// answer it got before.
    #[test]
    fn touching_a_census_changes_nothing_about_which_pileup_it_names() {
        let dir = tempfile::tempdir().expect("a scratch directory");
        let path = dir.path().join("touched.census");
        let identity = PileupIdentity::of_header(b"reference=A regions=1-100 filters=q20", 4_000);
        write_census(
            &every_corner(),
            Some(identity),
            &mut std::fs::File::create(&path).expect("a new file"),
        )
        .expect("a file accepts every write");

        let before = std::fs::metadata(&path).expect("the file exists");
        let (_, named) = open_census(&path).expect("this build's own file");

        // Touch it: the same bytes, written again, so only the timestamp moves.
        let bytes = std::fs::read(&path).expect("the file reads");
        std::fs::write(&path, &bytes).expect("the file is writable");
        let after = std::fs::metadata(&path).expect("the file exists");
        assert_eq!(before.len(), after.len(), "the same bytes");

        let (_, named_again) = open_census(&path).expect("this build's own file");
        assert_eq!(named_again, named);
        assert_eq!(named, Some(identity));
        assert_eq!(freshness(named_again, Some(identity)), Freshness::Fresh);
    }

    /// A census with no tracts at all, and one with no ordinary positions — the two ends of the
    /// range this caller works over, where a directory with one entry or an empty section could
    /// each be a special case nobody wrote code for.
    #[test]
    fn a_census_with_one_half_empty_round_trips() {
        let generic_only = SampleCensusEvidence::resident(
            "generic".to_string(),
            terms(),
            NamedReadGroup::drawn_for("generic", [ReadGroupId(0)]),
            BTreeMap::from([(
                SectionKey::Generic(ReadGroupId(0)),
                Section::Generic(GenericEvidence::never_walked(3)),
            )]),
        );
        assert_eq!(round_trip(&generic_only, None).census, generic_only);

        let tracts_only = SampleCensusEvidence::resident(
            "tracts".to_string(),
            terms(),
            NamedReadGroup::drawn_for("tracts", [ReadGroupId(0)]),
            BTreeMap::from([(
                SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS),
                Section::Ssr(SsrEvidence::never_walked(0)),
            )]),
        );
        assert_eq!(round_trip(&tracts_only, None).census, tracts_only);
    }

    /// **The read groups' names survive the round trip**, which is what a cohort of censuses is
    /// merged on.
    ///
    /// `every_corner_state_survives_a_round_trip` compares the whole value and so covers this
    /// too; this one exists because that comparison would go on passing if both sides recorded
    /// *no* names, and a census that names none is one no cohort can be built from.
    #[test]
    fn the_read_groups_names_survive_the_round_trip() {
        let census = SampleCensusEvidence::resident(
            "named".to_string(),
            terms(),
            BTreeMap::from([
                (
                    ReadGroupId(0),
                    NamedReadGroup {
                        declared_id: "HK5N7.1".to_string(),
                        library: "lib-A".to_string(),
                    },
                ),
                (
                    ReadGroupId(3),
                    NamedReadGroup {
                        declared_id: "HK5N7.2".to_string(),
                        library: "lib-B".to_string(),
                    },
                ),
            ]),
            BTreeMap::from([(
                SectionKey::Generic(ReadGroupId(0)),
                Section::Generic(GenericEvidence::never_walked(3)),
            )]),
        );

        let read = round_trip(&census, None);

        let named = read.census.declared_read_groups();
        assert_eq!(named.len(), 2, "both entries come back");
        assert_eq!(named[&ReadGroupId(0)].declared_id, "HK5N7.1");
        assert_eq!(named[&ReadGroupId(0)].library, "lib-A");
        assert_eq!(
            named[&ReadGroupId(3)].declared_id,
            "HK5N7.2",
            "an identifier that is not the entry's position comes back under its own number",
        );
        assert_eq!(named[&ReadGroupId(3)].library, "lib-B");
    }

    /// **A census naming one read group twice is malformed**, not a last-one-wins.
    ///
    /// Two entries under one identifier would name a section's read group two ways, and only one
    /// of the two could reach a cohort's merge — so the file is refused rather than half-read.
    ///
    /// **The duplicate is made by encoding the header's own bytes twice**, rather than by poking
    /// an offset this test computed: an offset would be a second copy of the layout, and would
    /// break on a change to the layout that this refusal does not care about.
    #[test]
    fn a_census_that_names_one_read_group_twice_is_refused() {
        let census = SampleCensusEvidence::resident(
            "named".to_string(),
            terms(),
            NamedReadGroup::drawn_for("named", [ReadGroupId(0)]),
            BTreeMap::from([(
                SectionKey::Generic(ReadGroupId(0)),
                Section::Generic(GenericEvidence::never_walked(1)),
            )]),
        );
        let mut bytes = Vec::new();
        write_census(&census, None, &mut bytes).expect("a vector accepts every write");

        // What one entry looks like on the wire, and where the count of them sits: both built
        // by the encoder rather than restated here.
        let mut one_entry = Vec::new();
        put_u32(&mut one_entry, 0);
        put_str(&mut one_entry, "named:rg0");
        put_str(&mut one_entry, "named:lib0");
        let mut count_at = Vec::new();
        count_at.extend_from_slice(MAGIC);
        put_u16(&mut count_at, VERSION);
        put_str(&mut count_at, "named");
        let at = count_at.len();

        assert_eq!(
            &bytes[at + 4..at + 4 + one_entry.len()],
            one_entry.as_slice(),
            "the one entry is where the encoder puts it",
        );
        bytes[at..at + 4].copy_from_slice(&2_u32.to_le_bytes());
        bytes.splice(at + 4..at + 4, one_entry);

        assert!(
            matches!(decode_census(&bytes), Err(CensusError::Malformed)),
            "a second entry under one identifier is a malformed file",
        );
    }
}
