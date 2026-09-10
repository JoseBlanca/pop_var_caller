//! The census on disk: a header, a directory, and the sections.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_records.md` §6.1 and §6.2; types in
//! `doc/devel/ng/arch/parameter_prepass_joint_records.md` §1.1a and §2.2.
//!
//! **One census per sample, and since `psp_census_pair.md` §3 it lives *inside* that sample's
//! psp**, as the file's closing payload. (It was a file of its own beside the psp until then;
//! nothing writes one now, and this codec is the same either way — a census is bytes, and where
//! they sit is the caller's.) It is a cache: everything in it
//! can be recomputed from the psp, and what it saves is a full decompression pass over that psp
//! every time a cohort is fitted.
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
//! Every integer is little-endian. **Offsets are from the census's own first byte**, not from the
//! start of whatever file it sits in — the directory is written before anyone knows where the
//! census will be put, so a reader that finds one at an offset adds that offset to every seek
//! ([`open_census_within`]).
//!
//! ```text
//!   magic     8 bytes   "NGCENSUS"
//!   version   u16       1
//!   header              the sample's name, its read groups, their minted read-error totals,
//!                       the twelve terms it was recorded under, and one byte that is always
//!                       zero (`encode_header`)
//!   directory u32 n, then n × (section key, offset u64, length u64)
//!   sections            the bytes each directory entry points at
//! ```
//!
//! **The directory is written before the sections and holds each one's offset**, so a reader
//! seeks once per section. That means the writer has to know each section's length before it
//! writes the directory, which it does by encoding the sections first and then placing them.
//!
//! # What is not here yet
//!
//! **Reading one section on its own.** This module writes the file and reads it whole; the
//! seeking reader that fills one section at a time is the next unit of work, and the whole-file
//! read here is its parity oracle.

use std::io::{Read, Seek, Write};
use std::path::Path;

use crate::ng::parameter_estimation::generic::calibration::MintedReadErrors;
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
/// **A version and not a feature flag.** A census can always be rebuilt from the psp that holds
/// it, so the answer to a version this build does not know is to rebuild rather than to
/// interpret.
///
/// **Three bumps so far, each one a census this build cannot read at all.** 2 on 2026-08-16,
/// when the depth code went from five bits on a widening ladder to eight on one with a bin for
/// every depth to the cap — a version-1 file's depth array is a different number of bytes for
/// the same position count and its codes index a different ladder, so nothing about it can be
/// salvaged by reading it more carefully. 3 on 2026-09-05, when the census began recording who
/// its read groups are; 4 the same day, when it began carrying what the base qualities claimed.
///
/// **Public because a psp's census is judged before it is decoded** (`psp_census_pair.md` §4.2):
/// [`CensusVerdict`](crate::ng::run::CensusVerdict) compares the version word at a trailer's
/// front against this, and names both numbers when they differ.
pub const VERSION: u16 = 4;

/// **The version, written out as a literal, so that changing it is a deliberate edit in two
/// places.**
///
/// A bump makes every census already written unreadable — every psp on disk needs
/// `regenerate-census` — so the one thing worth a test here is that the number did not move by
/// accident. Every other test compares the file's version word against the constant, so they all
/// pass whichever value it holds.
#[cfg(test)]
const THE_VERSION_THIS_BUILD_WRITES: u16 = 4;

/// **How many bytes of a census say which format it is** — the magic above and the version word
/// behind it, which is all [`version_word_of`] reads.
///
/// **Beside the two values it is made of**, so that widening either is a line away from the
/// number that says how far a reader must read to find them.
pub const BYTES_THAT_NAME_THE_VERSION: usize = MAGIC.len() + size_of::<u16>();

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
}

// ---------------------------------------------------------------------
// Writing
// ---------------------------------------------------------------------

/// Write one sample's census.
///
/// # Errors
///
/// [`CensusError::Io`] and nothing else: every value here is already well-formed, since the
/// types that carry them refuse to hold anything that is not.
pub fn write_census(
    census: &SampleCensusEvidence,
    out: &mut impl Write,
) -> Result<(), CensusError> {
    // **The sections are encoded before the directory is written**, because the directory holds
    // each one's length and its offset within the census, and neither is known until the bytes
    // exist.
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
    encode_header(&mut head, census);

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

fn encode_header(out: &mut Vec<u8>, census: &SampleCensusEvidence) {
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

    // **What each library's own base qualities claimed**, in a table of its own rather than a
    // field beside the names. A calibration is fitted from this and the run's measured error
    // rate together, so a census that carried the rate's evidence and not this one could produce
    // no calibration at all.
    //
    // **Its own length, so that "no entry" and "an entry that saw no reads" stay different on
    // the wire.** Writing one entry per declared group would turn the first into the second on
    // every round trip, and they are different claims: a group nothing was accumulated for
    // against a group whose reads all had a quality of zero.
    let minted = census.minted_read_errors();
    put_u32(out, minted.len() as u32);
    for (id, totals) in minted {
        put_u32(out, id.get());
        out.extend_from_slice(&totals.log_error_sum_scaled().to_le_bytes());
        put_u64(out, totals.reads());
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

    // **One byte that is always zero, and it stays.** It used to say whether the census named
    // the psp it was built from — a digest and a record count, so that a census kept in a file of
    // its own could be checked against that file. A census that *is* its psp's trailer cannot come
    // apart from it, so no shipped writer has set it since the census moved inside
    // (`psp_census_pair.md` §3), and none can now: there is nothing left to name a psp with.
    //
    // **Keeping it is what leaves every census this build has written byte for byte what it was.**
    // Dropping it would move every field of the directory and cost a format version, which makes
    // every psp already on disk unreadable, to save one byte a sample.
    out.push(0);
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
    let (sample, declared, minted, terms) = decode_header(&mut cursor)?;
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
        census: SampleCensusEvidence::resident(sample, terms, declared, minted, sections),
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
pub fn open_census(path: &Path) -> Result<SampleCensusEvidence, CensusError> {
    let whole = std::fs::metadata(path)?.len();
    open_census_within(path, ByteExtent::new(0, whole))
}

/// **The same, for a census that is a stretch of a larger file** — a psp's trailer
/// (`psp_census_pair.md` §5).
///
/// `census` says where in `path` the census's first byte is and how many bytes it occupies. Every
/// section this value later reads seeks to that offset plus the section's own, because the
/// directory's offsets are from the census's own front: it is written before anyone knows where
/// the census will be put.
///
/// **This is what a cohort too large to hold is opened with**, so it reads the census's head and
/// nothing more: a thousand samples' trailers read whole would be tens of gigabytes.
///
/// **The length is not there to make the read safe — it is there so the directory can be
/// checked.** Reading a few bytes past a census's end would be harmless on its own (the header
/// and the directory are decoded from the front of the buffer and whatever follows is never
/// looked at, and a short read at the end of a file is not an error). What the length buys is the
/// one check below: **every section has to end inside the census**, so a directory that points
/// outside it is refused here rather than turning into a seek into someone else's bytes and a
/// `resize` to a length the file supplied.
///
/// # Errors
///
/// [`CensusError::Io`] when the file will not open or read, [`CensusError::Malformed`] when the
/// bytes at the extent's offset are not a census this build reads, or when the directory places a
/// section outside the extent.
pub fn open_census_within(
    path: &Path,
    census: ByteExtent,
) -> Result<SampleCensusEvidence, CensusError> {
    // The header and the directory sit at the census's front, so only their bytes are read. How
    // many that is is not known until they are decoded, so the front is read in one go and the
    // rest is never touched. **Capped at the census's own length**, which matters only for a
    // census smaller than the buffer — above it the two are the same read.
    let mut file = std::fs::File::open(path)?;
    if census.offset() > 0 {
        file.seek(std::io::SeekFrom::Start(census.offset()))?;
    }
    let head_bytes = HEAD_READ_BYTES.min(usize::try_from(census.len()).unwrap_or(HEAD_READ_BYTES));
    let mut head = vec![0_u8; head_bytes];
    let filled = read_as_much_as_there_is(&mut file, &mut head)?;
    head.truncate(filled);

    let mut cursor = Cursor::new(&head);
    if cursor.take(MAGIC.len())? != MAGIC || cursor.u16()? != VERSION {
        return Err(CensusError::Malformed);
    }
    let (sample, declared, minted, terms) = decode_header(&mut cursor)?;
    let directory = decode_directory(&mut cursor)?;

    // **Every section ends inside the census, checked once here rather than at each read.** A
    // section read is a seek and a `resize` to a length the file supplied, so a directory naming
    // an extent the census does not hold would read another part of the file — or ask for an
    // allocation that aborts the process rather than returning an error.
    for (_, extent) in &directory {
        let ends_at = extent
            .offset()
            .checked_add(extent.len())
            .ok_or(CensusError::Malformed)?;
        if ends_at > census.len() {
            return Err(CensusError::Malformed);
        }
    }

    Ok(SampleCensusEvidence::backed(
        sample,
        terms,
        declared,
        minted,
        path.to_path_buf(),
        census.offset(),
        directory.into_iter().collect(),
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

/// **The version word of the census `head` begins with**, or `None` when those bytes are not a
/// census: they do not start with its magic, or they end before the word.
///
/// **This is what lets a census of another format be named as one rather than as damage**
/// (`psp_census_pair.md` §4.2). [`decode_census`] refuses both with the same
/// [`CensusError::Malformed`], which sends a user looking for a corrupted file when what happened
/// is that this build changed what a census holds. A judgement that reads the word first can tell
/// the two apart, and it costs [`BYTES_THAT_NAME_THE_VERSION`] bytes rather than a decode.
///
/// **It does not say the census is whole.** Everything past the version word — the header, the
/// directory, the sections — is unread here, so a truncated census of this version answers with
/// this build's own version, and the failure surfaces wherever something decodes it.
#[must_use]
pub fn version_word_of(head: &[u8]) -> Option<u16> {
    let (magic, rest) = head.split_at_checked(MAGIC.len())?;
    if magic != MAGIC {
        return None;
    }
    let word = rest.get(..size_of::<u16>())?;
    Some(u16::from_le_bytes([word[0], word[1]]))
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
    BTreeMap<ReadGroupId, MintedReadErrors>,
    RecordingTerms,
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

    let minted_count = cursor.u32()? as usize;
    let mut minted = BTreeMap::new();
    for _ in 0..minted_count {
        let id = ReadGroupId(cursor.u32()?);
        let sum = i128::from_le_bytes(
            cursor
                .take(16)?
                .try_into()
                .map_err(|_| CensusError::Malformed)?,
        );
        if minted
            .insert(id, MintedReadErrors::from_parts(sum, cursor.u64()?))
            .is_some()
        {
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

    // **The byte the pileup identity used to occupy** — see `encode_header`. It must be zero.
    // A census with it set was written before the census moved into the psp, names a psp by a
    // digest and a record count, and **this build does not read one**: the naming is what E2
    // deleted, so a reader that stepped over those 24 bytes would be accepting a file it can no
    // longer say anything true about. Nothing in the shipped commands can produce one — the fit
    // and the repair both read a psp's trailer, and the walk has written this byte zero since
    // Milestone A — so this refusal is reachable only from a census file kept from an older
    // build.
    if cursor.u8()? != 0 {
        return Err(CensusError::Malformed);
    }

    Ok((
        sample,
        declared,
        minted,
        RecordingTerms {
            selection: SelectionTermsDigest::from_fields(fields),
            kept_loci: CensusLociDigest::from_parts(whole, blocks),
            ssr_stratum_counts: StratumCounts::from_counted(counts),
            read_cap,
            depth_ladder,
            depth_cap: DepthCap::new(depth_cap),
        },
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

/// Fixtures shared by this module's tests and by those of the judgement that reads a psp's
/// trailer ([`ng::run::census_freshness`](crate::ng::run::census_freshness)).
///
/// **Here rather than duplicated**, for the reason `psp::writer`'s own `tests_support` gives:
/// a test that hand-wrote what a census begins with would keep passing after the magic or the
/// version word moved, and a moved version word is exactly what the judgement exists to
/// report.
#[cfg(test)]
pub(crate) mod tests_support {
    use super::*;

    use std::collections::BTreeMap;

    use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
    use crate::ng::parameter_estimation::joint::census::{DepthCode, RECORDED_OFFSET_RANGE};
    use crate::ng::parameter_estimation::joint::loci::CensusLociDigester;
    use crate::ng::types::{GenomePosition, Position};

    pub(super) const AT_SIX_REPEATS: Stratum = Stratum {
        period: 2,
        reference_repeats: 6,
    };
    pub(super) const ATG_FOUR_REPEATS: Stratum = Stratum {
        period: 3,
        reference_repeats: 4,
    };

    pub(super) fn selection_terms() -> crate::ng::parameter_estimation::joint::loci::SelectionTerms
    {
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
    pub(super) fn terms() -> RecordingTerms {
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
    pub(super) fn every_corner() -> SampleCensusEvidence {
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
            BTreeMap::new(),
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

    /// **The bytes of a census this build writes**, for a test elsewhere that needs a psp trailer
    /// this build would accept.
    ///
    /// **Written by [`write_census`] rather than spelled out**, because a test that hand-wrote a
    /// census's first bytes would keep passing after the magic or the version moved, and it is
    /// exactly a moved version word that the judgement it feeds
    /// ([`ng::run::census_freshness`](crate::ng::run::census_freshness)) exists to report.
    ///
    /// It is `every_corner`'s census, so it is also far longer than the ten bytes that judgement
    /// reads — which is the point of a trailer fixture for it.
    pub(crate) fn a_census_this_build_wrote() -> Vec<u8> {
        let mut bytes = Vec::new();
        write_census(&every_corner(), &mut bytes).expect("a vector accepts every write");
        bytes
    }
}

#[cfg(test)]
mod tests {
    use super::tests_support::*;
    use super::*;
    use crate::ng::parameter_estimation::joint::census::RECORDED_OFFSET_RANGE;

    use std::collections::BTreeMap;

    use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
    use crate::ng::parameter_estimation::joint::census::{
        CohortCensusEvidence, DepthCode, SsrLocusState,
    };

    fn round_trip(census: &SampleCensusEvidence) -> CensusFile {
        let mut bytes = Vec::new();
        write_census(census, &mut bytes).expect("a vector accepts every write");
        decode_census(&bytes).expect("what this build wrote, this build reads")
    }

    /// **The version has not moved.**
    ///
    /// Every other test here compares a file's version word against [`VERSION`], so they all pass
    /// whatever it holds. A bump is a real event — every census already written stops being
    /// readable and every psp on disk needs `regenerate-census` — and this is what makes it a
    /// deliberate edit in two places rather than a side effect of one.
    #[test]
    fn the_version_this_build_writes_has_not_moved() {
        assert_eq!(VERSION, THE_VERSION_THIS_BUILD_WRITES);
    }

    /// **A census that names the psp it was built from is refused, not read past.**
    ///
    /// The byte that said so is still in the header and this build always writes it zero
    /// (`encode_header`), so no census it writes can reach this. What can is a census file kept
    /// from a build before the census moved into the psp: it carries a digest and a record count
    /// that named a pairing this build no longer checks, and reading past them would be accepting
    /// a file it can say nothing true about.
    ///
    /// **The byte is found rather than guessed.** `encode_header` writes it last, so the header's
    /// own length places it, and the test cannot drift from the layout it is poking at.
    #[test]
    fn a_census_that_names_a_pileup_is_refused() {
        let census = every_corner();
        let mut bytes = a_census_this_build_wrote();

        let mut header = Vec::new();
        encode_header(&mut header, &census);
        let flag = MAGIC.len() + size_of::<u16>() + header.len() - 1;
        assert_eq!(bytes[flag], 0, "this build writes the byte absent");

        assert!(
            decode_census(&bytes).is_ok(),
            "the fixture decodes before the byte is touched, or this test proves nothing",
        );
        bytes[flag] = 1;
        assert!(
            matches!(decode_census(&bytes), Err(CensusError::Malformed)),
            "a census naming a pileup is refused",
        );
    }

    /// **The version word is read out of a census's first bytes, and out of nothing else.**
    /// Whether a psp's trailer is a census of this build's format, of another, or not a census
    /// at all is a judgement made on ten bytes (`psp_census_pair.md` §4.2), and each of the
    /// three answers is reachable.
    #[test]
    fn the_version_word_is_read_from_a_censuss_first_bytes_and_from_nothing_else() {
        let census = a_census_this_build_wrote();
        assert_eq!(version_word_of(&census), Some(VERSION));
        assert_eq!(
            version_word_of(&census[..BYTES_THAT_NAME_THE_VERSION]),
            Some(VERSION),
            "the word is inside the first {BYTES_THAT_NAME_THE_VERSION} bytes",
        );
        assert_eq!(
            version_word_of(&census[..BYTES_THAT_NAME_THE_VERSION - 1]),
            None,
            "bytes that end inside the version word do not name a version",
        );
        assert_eq!(
            version_word_of(b"a per-sample summary"),
            None,
            "a trailer holding something else is not a census",
        );
        assert_eq!(version_word_of(&[]), None);

        // **A word that is not this build's, read as the number it is.** Every assertion above
        // is satisfied by a function that answers `VERSION` whenever the magic matches, and
        // what a run does with an old census depends on which old version it names.
        let mut of_version_seven = census.clone();
        of_version_seven[MAGIC.len()..BYTES_THAT_NAME_THE_VERSION]
            .copy_from_slice(&7_u16.to_le_bytes());
        assert_eq!(version_word_of(&of_version_seven), Some(7));
    }

    /// **The assertion the whole step rests on**, and the reason B1 is its own commit: a codec
    /// that reads a field at the wrong offset produces a plausible number rather than a crash.
    #[test]
    fn every_corner_state_survives_a_round_trip() {
        let census = every_corner();
        let read = round_trip(&census);
        assert_eq!(read.census, census, "the whole sample, field for field");
    }

    /// **Decoding a census and encoding it again gives the bytes it was decoded from.**
    ///
    /// Every other round-trip test here compares decoded *values*, which is the right check for
    /// a codec on its own. This one is a check on the codec's bytes, and it exists because a
    /// caller leaned on it: the parity oracle used to decode the census *file* a rebuild had
    /// written and re-encode it before comparing it with the psp's trailer. **No caller decodes
    /// and re-encodes a census today** — since plan step D1 the rebuild writes into the trailer
    /// and the comparison is over those bytes directly. What this keeps guarding is the codec
    /// itself: a decode made lossy or normalising — a dropped empty section, a reordered
    /// directory — would be invisible to every round-trip test that compares decoded *values*.
    #[test]
    fn write_census_after_decode_census_returns_the_bytes_it_was_given() {
        let mut first = Vec::new();
        write_census(&every_corner(), &mut first).expect("a vector accepts every write");
        let read = decode_census(&first).expect("what this build wrote, this build reads");
        let mut again = Vec::new();
        write_census(&read.census, &mut again).expect("a vector accepts every write");
        assert_eq!(
            first, again,
            "decoding a census and encoding it again is the identity on bytes",
        );
    }

    /// The same value written twice is the same bytes — what §7.12's byte-for-byte comparison
    /// between the two builders will rest on, and what a directory whose offsets depended on a
    /// map's iteration order would break.
    #[test]
    fn writing_the_same_census_twice_gives_the_same_bytes() {
        let census = every_corner();
        let (mut first, mut second) = (Vec::new(), Vec::new());
        write_census(&census, &mut first).expect("a vector accepts every write");
        write_census(&census, &mut second).expect("a vector accepts every write");
        assert_eq!(first, second);
    }

    /// The three states an ordinary position can be in, read back off the file rather than off
    /// the value that was written — the distinction that has no field of its own.
    #[test]
    fn the_three_states_at_an_ordinary_position_come_back_apart() {
        let read = round_trip(&every_corner());
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
        let read = round_trip(&every_corner());
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
        write_census(&census, &mut bytes).expect("a vector accepts every write");

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

    /// **A file this build did not write is refused rather than decoded.** Every one of these
    /// would otherwise produce a plausible census: a truncated file reads a length off the end
    /// of a buffer, and a wrong version reads this build's fields at another build's offsets.
    #[test]
    fn a_stream_that_is_not_this_builds_census_is_refused() {
        let census = every_corner();
        let mut bytes = Vec::new();
        write_census(&census, &mut bytes).expect("a vector accepts every write");

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
            &mut std::fs::File::create(&path).expect("a new file"),
        )
        .expect("a file accepts every write");

        let mut backed = open_census(&path).expect("this build's own file");
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
        let mut backed = open_census(&path).expect("this build's own file");
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

    /// **A census read from the middle of a bigger file is the census it was written as**, and
    /// each section read still costs that section's bytes and no others.
    ///
    /// This is the shape a psp's trailer has (`psp_census_pair.md` §5): the census sits at an
    /// offset, with records before it and a footer after it, and its directory's offsets are
    /// relative to its own front because they were written before anyone knew where it would go.
    /// So every seek is the directory's offset plus where the census starts, and **a reader that
    /// forgot to add the offset would seek into the records and decode whatever is there** — a
    /// wrong section, or a malformed one, with nothing about the file being wrong.
    ///
    /// **The padding before the census is what does the work here**: without it a dropped offset
    /// would still land on the right section and this would pass. The padding after it makes the
    /// file the shape a psp is — something follows the census — and proves nothing on its own.
    #[test]
    fn a_census_at_an_offset_reads_the_same_as_one_at_the_front() {
        let dir = tempfile::tempdir().expect("a scratch directory");
        let path = dir.path().join("census-inside-something-else");
        let census = every_corner();

        let mut encoded = Vec::new();
        write_census(&census, &mut encoded).expect("a vector accepts every write");
        let before = vec![0x5a_u8; 4_097];
        let after = vec![0xa5_u8; 1_024];
        let mut whole = before.clone();
        whole.extend_from_slice(&encoded);
        whole.extend_from_slice(&after);
        std::fs::write(&path, &whole).expect("the scratch dir is ours");

        let at = before.len() as u64;
        let mut backed = open_census_within(&path, ByteExtent::new(at, encoded.len() as u64))
            .expect("this build wrote it");
        assert_eq!(backed.sample, census.sample, "the header is read at `at`");
        assert_eq!(backed.terms, census.terms);

        let groups = census.read_groups();
        let resident = decode_census(&encoded).expect("this build's own bytes");
        let mut from_memory = resident.census;
        assert_eq!(
            backed
                .with_generic(&groups, |sections| sections
                    .iter()
                    .map(|it| (*it).clone())
                    .collect::<Vec<_>>())
                .expect("this build wrote it"),
            from_memory
                .with_generic(&groups, |sections| sections
                    .iter()
                    .map(|it| (*it).clone())
                    .collect::<Vec<_>>())
                .expect("a resident census has no file to fail on"),
            "every ordinary-position section, read at an offset",
        );

        let strata = census.strata();
        assert_eq!(
            backed
                .with_strata(ReadGroupId(0), &strata, |sections| sections
                    .iter()
                    .map(|it| (*it).clone())
                    .collect::<Vec<_>>())
                .expect("this build wrote it"),
            from_memory
                .with_strata(ReadGroupId(0), &strata, |sections| sections
                    .iter()
                    .map(|it| (*it).clone())
                    .collect::<Vec<_>>())
                .expect("a resident census has no file to fail on"),
            "and every tract section",
        );

        // **And still only what was asked for.** The offset changes where a section is, not how
        // much of the file a call touches — which is the property that lets a cohort too large to
        // hold be opened at all.
        let directory = decode_directory_of(&encoded).expect("this build's own bytes");
        let one = directory
            .iter()
            .find(|(key, _)| *key == SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS))
            .map(|(_, extent)| extent.len())
            .expect("the section is in the directory");
        reset_bytes_read();
        backed
            .with_strata(ReadGroupId(0), &[AT_SIX_REPEATS], |sections| {
                sections[0].len()
            })
            .expect("this build wrote it");
        assert_eq!(
            bytes_read(),
            one,
            "one stratum's read is that stratum's bytes and no more — the counter covers the \
             section reads, not the head read that opening does",
        );
    }

    /// **A directory that places a section outside the census is refused when it is opened**,
    /// rather than seeking there when a fit asks for it.
    ///
    /// A section read is a seek and a `resize` to a length the file supplied. Left unchecked, a
    /// census whose directory outgrew it would read a stretch of whatever it sits in — a psp's
    /// records — and, for a large enough length, ask the allocator for it, which aborts the
    /// process instead of returning an error. **The length the caller passes is what makes the
    /// check possible**, and this is what it is for.
    #[test]
    fn a_section_that_ends_outside_the_census_is_refused_at_the_door() {
        let dir = tempfile::tempdir().expect("a scratch directory");
        let path = dir.path().join("truncated.census");
        let mut encoded = Vec::new();
        write_census(&every_corner(), &mut encoded).expect("a vector accepts every write");
        std::fs::write(&path, &encoded).expect("the scratch dir is ours");

        // The same file, opened as though it were one byte shorter than it is — so the last
        // section the directory names ends past the end.
        let refused = open_census_within(&path, ByteExtent::new(0, encoded.len() as u64 - 1))
            .expect_err("the last section ends outside a census this length");
        assert!(
            matches!(refused, CensusError::Malformed),
            "and got: {refused:?}",
        );

        open_census_within(&path, ByteExtent::new(0, encoded.len() as u64))
            .expect("at its real length it opens");
    }

    /// **A census assembled into a cohort keeps where it is in its file.**
    ///
    /// Every sample of a cohort goes through this: a census is written under its own walk's
    /// read-group identifiers, and building a cohort relabels each one onto run-wide identifiers,
    /// which rebuilds its directory. **An offset dropped there would send every seek to byte zero
    /// of the psp** — into the header and the records, which decode as something or as nothing,
    /// with no file being wrong.
    ///
    /// It is asserted by reading a section *after* the cohort is built, because that is the only
    /// thing that touches the file. The other tests here read a census that was never relabelled,
    /// so none of them can see this.
    #[test]
    fn a_census_in_a_cohort_keeps_where_it_is_in_its_file() {
        let dir = tempfile::tempdir().expect("a scratch directory");
        let samples = drawn_cohort();

        // **Each sample at a different offset**, so a reader that dropped the offset could not
        // pass by accident on a cohort that happened to agree at zero.
        let mut backed = Vec::new();
        for (which, sample) in samples.iter().enumerate() {
            let path = dir.path().join(format!("{}.inside", sample.sample));
            let mut encoded = Vec::new();
            write_census(sample, &mut encoded).expect("a vector accepts every write");
            let before = vec![0x5a_u8; 1_024 + which * 997];
            let mut whole = before.clone();
            whole.extend_from_slice(&encoded);
            std::fs::write(&path, &whole).expect("the scratch dir is ours");
            backed.push(
                open_census_within(
                    &path,
                    ByteExtent::new(before.len() as u64, encoded.len() as u64),
                )
                .expect("this build wrote it"),
            );
        }
        let mut cohort = CohortCensusEvidence::new(backed).expect("every sample recorded one way");
        let mut from_memory =
            CohortCensusEvidence::new(samples).expect("and so did the same ones in memory");

        let groups: Vec<ReadGroupId> = cohort.read_groups().to_vec();
        assert!(!groups.is_empty(), "the fixture declares read groups");
        // **The depth arrays, sample by sample** — the cohort lends its samples' sections rather
        // than handing them over, so what is compared is what a fit reads through it.
        let depths_of = |evidence: &mut CohortCensusEvidence| {
            evidence
                .with_generic(&groups, |samples| {
                    samples
                        .iter()
                        .map(|sample| {
                            sample
                                .iter()
                                .map(|(group, section)| {
                                    (*group, section.depth().as_bytes().to_vec())
                                })
                                .collect::<Vec<_>>()
                        })
                        .collect::<Vec<_>>()
                })
                .expect("this build wrote it")
        };
        assert_eq!(
            depths_of(&mut cohort),
            depths_of(&mut from_memory),
            "a census read out of the middle of a file gives the cohort the sections it holds",
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
                &mut std::fs::File::create(&path).expect("a new file"),
            )
            .expect("a file accepts every write");
            backed.push(open_census(&path).expect("this build's own file"));
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
                    BTreeMap::new(),
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
            &mut std::fs::File::create(&path).expect("a new file"),
        )
        .expect("a file accepts every write");
        let file_len = std::fs::metadata(&path).expect("the file exists").len();

        let mut backed = open_census(&path).expect("this build's own file");
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

    /// A census with no tracts at all, and one with no ordinary positions — the two ends of the
    /// range this caller works over, where a directory with one entry or an empty section could
    /// each be a special case nobody wrote code for.
    #[test]
    fn a_census_with_one_half_empty_round_trips() {
        let generic_only = SampleCensusEvidence::resident(
            "generic".to_string(),
            terms(),
            NamedReadGroup::drawn_for("generic", [ReadGroupId(0)]),
            BTreeMap::new(),
            BTreeMap::from([(
                SectionKey::Generic(ReadGroupId(0)),
                Section::Generic(GenericEvidence::never_walked(3)),
            )]),
        );
        assert_eq!(round_trip(&generic_only).census, generic_only);

        let tracts_only = SampleCensusEvidence::resident(
            "tracts".to_string(),
            terms(),
            NamedReadGroup::drawn_for("tracts", [ReadGroupId(0)]),
            BTreeMap::new(),
            BTreeMap::from([(
                SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS),
                Section::Ssr(SsrEvidence::never_walked(0)),
            )]),
        );
        assert_eq!(round_trip(&tracts_only).census, tracts_only);
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
            // What each library's own base qualities claimed, which rides beside the names.
            BTreeMap::from([(ReadGroupId(3), MintedReadErrors::from_parts(-4_096, 7))]),
            BTreeMap::from([(
                SectionKey::Generic(ReadGroupId(0)),
                Section::Generic(GenericEvidence::never_walked(3)),
            )]),
        );

        let read = round_trip(&census);

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

        let minted = read.census.minted_read_errors();
        assert_eq!(
            minted[&ReadGroupId(3)],
            MintedReadErrors::from_parts(-4_096, 7),
            "the sum comes back as the scaled integer it went in as, not through a float",
        );
        assert!(
            !minted.contains_key(&ReadGroupId(0)),
            "a group the walk accumulated nothing for has no entry at all, which is a different \
             claim from an entry whose sum is zero over some reads — and the table carries its \
             own length so the two stay apart on the wire",
        );
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
            BTreeMap::new(),
            BTreeMap::from([(
                SectionKey::Generic(ReadGroupId(0)),
                Section::Generic(GenericEvidence::never_walked(1)),
            )]),
        );
        let mut bytes = Vec::new();
        write_census(&census, &mut bytes).expect("a vector accepts every write");

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
