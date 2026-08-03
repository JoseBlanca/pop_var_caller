//! `AlignmentFile::open` — the validate-on-open gate, and the handle it
//! produces.
//!
//! Every check lives inside the function that opens the file, so a file is
//! either opened *and* validated or it is an `Err`: there is no window in which
//! an unvalidated handle exists. The checks run fail-fast in this order —
//! `@HD SO`, `@SQ`↔reference, the index, `@RG SM`
//! (`doc/devel/ng/spec/alignment_file.md` §3.1).
//!
//! **The invariant this establishes** is what every later layer leans on: once
//! the gate passes, `ref_id == ContigId` holds by construction for every record
//! the file can yield. That is what makes it sound for the merge one layer up
//! to compare positions across files without re-checking anything.
//!
//! ng reads `@HD SO`, `@SQ` and `@RG SM` off the noodles `sam::Header` itself
//! rather than reusing production's extractors: those are module-private, and
//! making them visible would be a production edit the ng freeze forbids
//! (`doc/devel/ng/arch/alignment_file.md` §5).
//!
//! The name says BAM but the module opens CRAM too — "BAM" in the everyday
//! sense of "the alignment file" (spec §6).

use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use noodles_bam as bam;
use noodles_cram as cram;
use noodles_sam as sam;

use crate::bam::index_preflight::{
    AlignmentFileKind, AlignmentIndex, load_alignment_index, preflight_alignment_indexes,
};
use crate::fasta::{ContigEntry, ContigList};
use crate::ng::read::filtering::ReadFilterConfig;
use crate::ng::read::input::aligned_reads_reader::{
    AlignedReadsReader, BamAlignedReadsReader, CramAlignedReadsReader,
};
use crate::ng::read::input::cursor::AlignmentCursor;
use crate::ng::read::input::read_groups::ReadGroupResolution;
use crate::ng::read::input::reference::{OpenReference, ReferenceBasesError};
use crate::ng::ref_seq::RawRefSeq;
use crate::ng::types::ContigId;

use super::AlignmentFileError;

/// One opened, **validated** alignment file.
///
/// The gate of [`AlignmentFile::open`] has passed, so `ref_id == ContigId`
/// holds for every record this file will ever yield — the invariant the region
/// query, the order guard and the cross-file merge all lean on without
/// re-checking. There is no way to hold one of these unvalidated: `open` either
/// returns a validated handle or an error.
///
/// The parsed index is owned here and lives for the whole run: minting a cursor
/// is an in-memory lookup plus one file open, never an index parse (spec §3.3).
/// Not `Clone` — it owns an index, and two copies would be two parses of it.
///
/// **Sharing is by `Arc`, not by reference**, because a cursor outlives the call
/// that made it and carries no lifetime, so the file it was minted from has to be
/// shareable rather than borrowed ([`cursor`](Self::cursor)).
pub struct AlignmentFile {
    /// `Arc` so the per-query order guard can hold it for its error message
    /// without an allocation per query.
    path: Arc<Path>,
    /// Kept for the region query, which resolves a contig to a `ref_id` and
    /// hands the header to the record source. Read from C2 onwards.
    ///
    /// `Arc` so a region source can **own** its header rather than borrow this
    /// one: a borrowed header ties the returned stream to the file's lifetime,
    /// and a resumable generator has to hold that stream across calls
    /// (`doc/devel/ng/arch/locus_generation_pileup.md` §2.2). An independent
    /// `Arc`, cloned out per query — *not* a reference into an `Arc`'d file,
    /// which is what would make the source self-referential. Parsed once at
    /// open either way; the clone is one atomic increment per query.
    header: Arc<sam::Header>,
    /// Parsed once, at open — never re-read per query (spec §3.3). Queried from
    /// C2 onwards, which is the guarantee the whole per-query cost model rests
    /// on: a query is an in-memory lookup plus a seek.
    index: AlignmentIndex,
    /// How this file's records are assigned to read groups. Settled at open and
    /// never recomputed; the record sources consult it per record only when it
    /// says they must.
    ///
    /// A region source **owns** what it consults, so the stream it serves does
    /// not borrow this file — but the sharing lives inside
    /// [`ReadGroupResolution`] rather than around it, so this field is a plain
    /// value and a per-query copy costs no atomic at all in the common
    /// single-read-group case. The header, which has no such cheap-clone form,
    /// is the one that needs an `Arc` of its own.
    resolution: ReadGroupResolution,
    /// **This file's own `@SQ` list**, which the open gate reconciled against
    /// the reference's: same names, same lengths, same order, and the same
    /// digests wherever both sides carry one (`alignment_file.md` §3.1,
    /// check 2). So `entries[i]` is `ContigId(i)`, and a caller asking which
    /// chromosomes it may make a cursor for gets the answer from the file it is
    /// about to read rather than from a reference handle it would have to carry
    /// alongside.
    ///
    /// **Reconciled is not identical, and the difference is the digests.**
    /// `first_disagreement` treats an absent `M5` as a wildcard
    /// (`fasta/mod.rs`), so a file declaring digests passes against a
    /// `.fai`-only reference that has none — and what is stored here is then
    /// the *file's* claim, not the reference's. Names, lengths and order are
    /// the part that is proved equal, and they are the part a cursor needs.
    contigs: ContigList,
    /// The `@SQ M5` tags, indexed by `ContigId`. `None` where the file carries
    /// no usable `M5`.
    ///
    /// **A projection of `contigs` above, and it should not have to exist.**
    /// It is stored because `sq_md5s()` hands out a *slice*, which a value built
    /// per call could not outlive the call to lend, and because
    /// `SampleReads::assembly_inputs` (`mod.rs`) lends one from every open file
    /// at once for `check_assembly` — which is `pub` and takes
    /// `&[Option<[u8; 16]>]`. Deleting the field and having `check_assembly`
    /// take a `&ContigList` was tried during review and works (the suite stays
    /// green, ~18 lines shorter), but it changes a public signature and eight
    /// call sites, which is more than the step that introduced `contigs` should
    /// carry. Recorded as a follow-up rather than left as an unexplained pair.
    ///
    /// Captured at open and compared much later, by `check_assembly` (D1),
    /// once the caller has joined `reference_info`'s background verification.
    sq_md5s: Vec<Option<[u8; 16]>>,
    /// Handed to each cursor's `ReadFilter`. Held on the file rather than passed
    /// per cursor because the filtering policy is the file's for the whole run.
    filter_config: ReadFilterConfig,
    /// This CRAM's `.crai` entries, grouped by contig — `crai_by_contig[i]`
    /// holds contig `i`'s entries in file order. Empty for a BAM.
    crai_by_contig: Vec<Arc<[cram::crai::Record]>>,
    /// **The run's reference**, held so each cursor can ask it for the bases
    /// narrowed to that cursor's contig — not a repository of this file's own.
    ///
    /// A handle, not a copy: it is an `Arc` inside, so every file in a run
    /// points at one cache of bases. That is what keeps a cohort's resident
    /// reference at one contig rather than `files × genome` (see
    /// `reference.rs`, and the note at the open site).
    ///
    /// `None` for a BAM, which stores its own sequences and needs no
    /// reference to decode.
    reference: Option<OpenReference>,
}

impl AlignmentFile {
    /// Open one indexed BAM/CRAM and **validate it or fail**.
    ///
    /// Four checks, fail-fast **in this order** (spec §3.1): `@HD SO` is
    /// `coordinate`; the `@SQ` list equals `reference.contig_list()` exactly,
    /// order included; the index loads; the `@RG` records name exactly one
    /// sample. The order is cheapest-and-most-fundamental first, so a file that
    /// is wrong in several ways reports the most basic fault rather than
    /// whichever check happens to run.
    ///
    /// A CRAM then needs one more thing before it can be read at all — the
    /// reference *bases*, to decode against — so the repository is built after
    /// those four, and a reference that cannot supply one is rejected here
    /// rather than at the first query.
    ///
    /// **The `@SQ` check is the permutation fix.** Comparing *in order* against
    /// the reference is what distinguishes this from the resolves-only check it
    /// replaces: a file whose contig list is a re-ordering of the reference's
    /// resolves every index and then fetches the wrong contig for every read.
    ///
    /// With `build_index_if_missing`, an absent index is built next to the file
    /// rather than rejected; that is a caller policy, not this module's.
    ///
    /// **Hands back an `Arc`, because that is the only shape the type is
    /// useful in.** Minting a cursor takes `self: &Arc<Self>`
    /// ([`cursor`](Self::cursor)), so a bare `AlignmentFile` could be asked for
    /// its path, its digests and its contig table and nothing else. Putting the
    /// share in the constructor states that invariant once, where a caller
    /// cannot skip it, instead of leaving every call site to wrap the result
    /// itself.
    pub fn open(
        path: &Path,
        reference: &OpenReference,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
        resolution: ReadGroupResolution,
    ) -> Result<Arc<Self>, AlignmentFileError> {
        let header = read_header(path)?;

        // 1. @HD SO.
        let observed_sort_order = sort_order(&header);
        if observed_sort_order.as_deref() != Some("coordinate") {
            return Err(AlignmentFileError::NotCoordinateSorted {
                path: path.to_path_buf(),
                sort_order: observed_sort_order,
            });
        }

        // 2. @SQ vs the reference's contig table — names, lengths, order, and
        //    digests where both sides carry one. `first_disagreement` applies
        //    the absent-digest-is-a-wildcard rule itself (`fasta/mod.rs`, its
        //    `if let (Some, Some)` arm) rather than going through
        //    `ContigEntry`'s `PartialEq`, which encodes the same rule
        //    separately — so a `.fai`-only reference compares on name and
        //    length alone with no branch needed here.
        //
        //    The file is the left operand so the message reads "file value vs
        //    reference value", the direction a user needs: the reference is the
        //    authority and the file is the thing that is wrong.
        let file_contigs =
            contig_list(&header).map_err(|bad| AlignmentFileError::MalformedMd5 {
                path: path.to_path_buf(),
                contig: bad.contig,
                detail: bad.detail,
            })?;
        file_contigs
            .first_disagreement(&reference.info().contig_list())
            .map_err(|detail| AlignmentFileError::ContigReconcile {
                path: path.to_path_buf(),
                detail,
            })?;

        // 3. The index — parsed here, once, and held for the life of the run.
        //
        // Building a missing index reaches `noodles_cram::fs::index`, which
        // decodes *multi-reference* slices against an empty
        // `fasta::Repository` (marked `// TODO` in noodles 0.93) and therefore
        // **panics** on a CRAM whose reads are stored as differences from the
        // reference and whose slices span contigs — plausible on a fragmented
        // assembly. Building the index outside ng avoids it. Noted here
        // because this is the only ng call that can reach that path.
        if build_index_if_missing {
            preflight_alignment_indexes(&[path.to_path_buf()], true).map_err(|source| {
                AlignmentFileError::Index {
                    path: path.to_path_buf(),
                    source,
                }
            })?;
        }
        let index = load_alignment_index(path).map_err(|source| AlignmentFileError::Index {
            path: path.to_path_buf(),
            source,
        })?;

        // 4. The `@RG` records were read, validated and identified before any
        //    file was opened, by the read-group pre-pass — so there is no check
        //    left here. What used to be checked at this point (exactly one
        //    `@RG SM`) has moved twice over: a file may now declare several read
        //    groups, and naming one sample is a property of the *open*, enforced
        //    by `SampleReads` (spec §4, §8).
        //
        // 5. A CRAM needs the reference *bases* to decode at all, so the
        //    repository is taken here — from the **run's** `OpenReference`, not
        //    built per file — and a reference that cannot supply one is a hard
        //    error now rather than a mystery at the first query.
        //
        //    Asking `OpenReference` rather than building here is the whole of
        //    the memory fix. A `fasta::Repository` memoises whole contigs and
        //    never evicts, so a per-file repository costs `files × genome`:
        //    measured at ~752 MiB per open file against the 746 MiB tomato
        //    reference, which is a 51-sample cohort dying at 38 GiB. Sharing
        //    the run's one makes that one contig, once (`reference.rs`).
        //
        //    Nothing is *read* here — this only opens the FASTA and proves it
        //    can be, so "this CRAM has no bases to decode against" is a fault
        //    at open rather than a mystery at the first query. The bases
        //    themselves arrive per contig, at query time. A BAM never asks, so
        //    a BAM-only run still never touches the FASTA.
        let crai_by_contig = match &index {
            AlignmentIndex::Crai(crai) => group_crai_by_contig(crai, file_contigs.entries.len()),
            _ => Vec::new(),
        };

        let file_reference = match AlignmentFileKind::from_path(path) {
            Some(AlignmentFileKind::Cram) => {
                reference.bases().map_err(|source| match source {
                    ReferenceBasesError::NoFasta => AlignmentFileError::CramNeedsReferenceFasta {
                        path: path.to_path_buf(),
                    },
                    ReferenceBasesError::Build { fasta, source } => AlignmentFileError::Open {
                        path: fasta,
                        source: std::io::Error::other(source),
                    },
                })?;
                Some(reference.clone())
            }
            _ => None,
        };

        // Free, and only sound now: check 2 proved this order is the
        // reference's, so position i really is `ContigId(i)`.
        let sq_md5s = file_contigs.entries.iter().map(|entry| entry.md5).collect();

        Ok(Arc::new(Self {
            path: Arc::from(path),
            header: Arc::new(header),
            index,
            resolution,
            contigs: file_contigs,
            sq_md5s,
            filter_config,
            crai_by_contig,
            reference: file_reference,
        }))
    }

    /// How this file's records are assigned to read groups (see the field).
    pub fn read_group_resolution(&self) -> &ReadGroupResolution {
        &self.resolution
    }

    /// The chromosomes this file may be read over.
    ///
    /// **The file's own `@SQ` list**, which the open gate proved agrees with the reference's
    /// on names, lengths and order, and on digests wherever both sides carry one
    /// (`alignment_file.md` §3.1, check 2). The digests are therefore the *file's*: a file may
    /// declare an `M5` where a `.fai`-only reference has none, and an absent digest is a
    /// wildcard on either side. Use [`sq_md5s`](Self::sq_md5s) and `check_assembly` to compare
    /// digests against a verified reference; use this to know which contigs exist and what
    /// `ContigId(i)` means.
    ///
    /// Position `i` is `ContigId(i)`, so a caller minting one cursor per chromosome can walk
    /// this without holding the reference open beside it.
    pub fn contigs(&self) -> &ContigList {
        &self.contigs
    }

    /// One cursor for one chromosome of this file.
    ///
    /// **Called once per chromosome per worker, never per region — that is the whole point.**
    /// A cursor opens its own descriptor and keeps it for as long as it lives, so it can stay
    /// positioned between regions and hand back reads it has already decoded and filtered
    /// rather than seeking to a block it is usually already sitting in.
    ///
    /// Fallible only here. Opening the descriptor is the one thing that can fail at
    /// construction; after this returns, a cursor cannot fail to exist.
    ///
    /// **CRAM is not served yet** (Milestone E). It is refused rather than silently mis-read,
    /// because a CRAM opened as a BAM would fail deep inside a decode with an error naming
    /// neither the format nor this decision.
    pub fn cursor<R: RawRefSeq>(
        self: &Arc<Self>,
        contig: ContigId,
        reference: R,
    ) -> Result<AlignmentCursor<R>, AlignmentFileError> {
        let open_error = |source: std::io::Error| AlignmentFileError::Open {
            path: self.path.to_path_buf(),
            source,
        };

        // A contig this file does not declare has no reads and no chunks, and every layer
        // below would answer that as "nothing here" — which is indistinguishable from a
        // chromosome that is genuinely empty. Refused where the caller can still act on it.
        if usize::try_from(contig.get())
            .ok()
            .is_none_or(|id| id >= self.contigs.entries.len())
        {
            return Err(AlignmentFileError::CursorContigNotInFile {
                path: self.path.to_path_buf(),
                contig,
                contigs_in_file: self.contigs.entries.len(),
            });
        }

        let aligned_reads_reader = match AlignmentFileKind::from_path(&self.path) {
            Some(AlignmentFileKind::Bam) => {
                let mut reader = bam::io::reader::Builder
                    .build_from_path(&self.path)
                    .map_err(open_error)?;
                reader.read_header().map_err(open_error)?;
                AlignedReadsReader::Bam(BamAlignedReadsReader::new(
                    reader,
                    Arc::clone(&self.header),
                    self.index.clone(),
                    Arc::clone(&self.path),
                ))
            }
            Some(AlignmentFileKind::Cram) => {
                let mut reader = cram::io::reader::Builder::default()
                    .build_from_path(&self.path)
                    .map_err(open_error)?;
                reader.read_header().map_err(open_error)?;
                // **The bases are taken once, for this cursor's chromosome, and held for its
                // life.** A CRAM decodes against the reference, and a cursor covers one
                // chromosome — so asking here rather than per region is the whole of what the
                // cursor changes about the reference (spec §10). Asking the *run's*
                // `OpenReference`, rather than holding a repository of our own, is what lets
                // the run drop a chromosome's bases when every cursor on it is gone.
                let repository = self
                    .reference
                    .as_ref()
                    .and_then(|reference| reference.bases_for_contig(contig))
                    .ok_or_else(|| AlignmentFileError::Open {
                        path: self.path.to_path_buf(),
                        source: std::io::Error::other(
                            "a CRAM was opened without reference bases to decode against",
                        ),
                    })?;
                let entries = self
                    .crai_by_contig
                    .get(usize::try_from(contig.get()).unwrap_or(usize::MAX))
                    .cloned()
                    .unwrap_or_else(|| Vec::new().into());
                AlignedReadsReader::Cram(CramAlignedReadsReader::new(
                    reader,
                    Arc::clone(&self.header),
                    repository,
                    entries,
                    self.resolution.clone(),
                    Arc::clone(&self.path),
                ))
            }
            _ => {
                return Err(AlignmentFileError::CursorFormatUnsupported {
                    path: self.path.to_path_buf(),
                    kind: "this format",
                });
            }
        };

        AlignmentCursor::over_records(
            aligned_reads_reader,
            contig,
            self.resolution.clone(),
            reference,
            self.filter_config,
            Arc::clone(&self.path),
        )
        .map_err(|source| AlignmentFileError::Reference {
            path: self.path.to_path_buf(),
            source,
        })
    }

    /// The `@SQ M5` tags, indexed by `ContigId`, for the deferred assembly
    /// check. `None` where the file carried no usable digest.
    pub fn sq_md5s(&self) -> &[Option<[u8; 16]>] {
        &self.sq_md5s
    }

    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Hand-written so the parsed index does not have to be `Debug`, and so the
/// output says what identifies the file rather than dumping one.
impl std::fmt::Debug for AlignmentFile {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Destructured exhaustively so a field that is added, removed or
        // retyped is a compile error *here* rather than a silent omission from
        // the output. What is printed is a deliberate subset — the `_` bindings
        // are the ones deliberately left out — and only an explicit list can
        // say which is which.
        let Self {
            path,
            header: _,
            index: _,
            resolution,
            contigs,
            sq_md5s: _,
            filter_config: _,
            crai_by_contig: _,
            reference: _,
        } = self;

        f.debug_struct("AlignmentFile")
            .field("path", path)
            .field("read_groups", resolution)
            // Counted off `contigs` now rather than off `sq_md5s`. The two are the same
            // length by construction, but one of them *is* the contig list and the other is
            // a projection of it, and the field is labelled "contigs".
            .field("contigs", &contigs.entries.len())
            .finish_non_exhaustive()
    }
}

/// Bucket a `.crai`'s entries by contig, keeping each contig's entries in file
/// order. One O(n) pass, at open.
///
/// **This replaces a binary search, and the reason is correctness rather than
/// speed.** Searching assumes the `.crai` is ordered by contig with unplaced
/// entries last — and that assumption is false for an index noodles itself
/// writes: within a multi-reference slice it emits entries sorted by
/// `Option<usize>`, and `None < Some(0)` in Rust, so the unplaced entry comes
/// *first*. A `partition_point` over an unpartitioned slice returns an
/// unspecified index, and the walk would then find a foreign entry and report
/// end-of-input — **losing every read of the region, with no error**.
///
/// Grouping assumes nothing about the order between contigs. It also makes the
/// per-query lookup O(1) rather than O(log n), which matters at one query per
/// STR locus. Within a contig the file order is kept, because the
/// container-level early stop does rely on containers ascending by start —
/// which follows from the file being coordinate-sorted, the same assumption the
/// BAM early stop makes.
///
/// Unplaced entries are dropped: they can never overlap a region.
pub(crate) fn group_crai_by_contig(
    index: &cram::crai::Index,
    contig_count: usize,
) -> Vec<Arc<[cram::crai::Record]>> {
    let mut by_contig: Vec<Vec<cram::crai::Record>> = vec![Vec::new(); contig_count];

    for record in index.iter() {
        if let Some(contig) = record.reference_sequence_id()
            && contig < contig_count
        {
            by_contig[contig].push(record.clone());
        }
    }

    by_contig.into_iter().map(Arc::from).collect()
}

/// Read just the SAM header, dispatching on the file's extension.
///
/// The reader is opened and dropped here: the gate needs the header, and the
/// readers a query will use come from the pool (step C1).
///
/// `pub(crate)` for the read-group pre-pass, which reads every input's header —
/// and only its header — before any file is opened for reading.
pub(crate) fn read_header(path: &Path) -> Result<sam::Header, AlignmentFileError> {
    let open_error = |source: std::io::Error| AlignmentFileError::Open {
        path: path.to_path_buf(),
        source,
    };

    let file = File::open(path).map_err(open_error)?;

    match AlignmentFileKind::from_path(path) {
        Some(AlignmentFileKind::Bam) => bam::io::Reader::new(file).read_header(),
        Some(AlignmentFileKind::Cram) => cram::io::Reader::new(file).read_header(),
        // `load_alignment_index` rejects this too, but the header read comes
        // first, so the extension has to be understood here.
        None => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "expected a '.bam' or '.cram' extension",
        )),
    }
    .map_err(open_error)
}

/// The `@HD SO` value, or `None` when the header carries no `@HD` line or no
/// sort-order tag at all.
///
/// The value is returned rather than compared here so the caller owns the
/// policy: "missing" and "queryname" are the same rejection, but they are not
/// the same *diagnosis*, and only the caller knows whether its error can say
/// which it was.
pub(crate) fn sort_order(header: &sam::Header) -> Option<String> {
    use sam::header::record::value::map::header::tag::SORT_ORDER;

    let raw = header.header()?.other_fields().get(&SORT_ORDER)?;
    Some(String::from_utf8_lossy(raw.as_ref()).into_owned())
}

/// The header's `@SQ` list as a [`ContigList`] — name, length, and the `M5`
/// digest where the file carries one.
///
/// The order is the header's own, which is the whole point: the gate compares
/// this list against the reference's **including order**, and that is what
/// catches a permutation. Building it in header order is therefore load-bearing
/// rather than incidental.
///
/// **A missing `M5` is fine; a malformed one is fatal.** Spec §3.1 settles the
/// missing case: never an error, never a warning, because the digest check is
/// opportunistic and refusing a file for not offering a bonus would punish the
/// common case. A tag that is *present* but not 32 hex characters is a
/// different thing — a header its writer got wrong — and reading it as absent
/// would pass an error silently. The caller turns this into
/// `MalformedMd5`, naming the contig.
pub(crate) fn contig_list(header: &sam::Header) -> Result<ContigList, MalformedMd5Tag> {
    let entries = header
        .reference_sequences()
        .iter()
        .map(|(name, reference_sequence)| {
            let name = String::from_utf8_lossy(name.as_ref()).into_owned();
            Ok(ContigEntry {
                length: usize::from(reference_sequence.length()) as u64,
                md5: md5_tag(reference_sequence).map_err(|detail| MalformedMd5Tag {
                    contig: name.clone(),
                    detail,
                })?,
                name,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(ContigList { entries })
}

/// An `@SQ` line whose `M5` tag could not be read as a digest. Carries the
/// facts; the path and the phrasing belong to the caller's error.
#[derive(Debug)]
pub(crate) struct MalformedMd5Tag {
    pub(crate) contig: String,
    pub(crate) detail: String,
}

/// The `M5` tag of one `@SQ` entry, decoded from hex.
///
/// `Ok(None)` = no tag, which is ordinary. `Err(detail)` = a tag that is there
/// but unreadable, which is not.
fn md5_tag(
    reference_sequence: &sam::header::record::value::Map<
        sam::header::record::value::map::ReferenceSequence,
    >,
) -> Result<Option<[u8; 16]>, String> {
    use sam::header::record::value::map::reference_sequence::tag::MD5_CHECKSUM;

    let Some(hex) = reference_sequence.other_fields().get(&MD5_CHECKSUM) else {
        return Ok(None);
    };
    let hex = hex.as_ref();

    decode_md5_hex(hex).map(Some).ok_or_else(|| {
        if hex.len() != 32 {
            format!("expected 32 hex characters, got {}", hex.len())
        } else {
            "contains a non-hex character".to_string()
        }
    })
}

/// Decode 32 hex characters into the 16 digest bytes they spell. `None` for any
/// other length or a non-hex character.
fn decode_md5_hex(hex: &[u8]) -> Option<[u8; 16]> {
    let hex: &[u8; 32] = hex.try_into().ok()?;

    let mut digest = [0u8; 16];
    for (byte, pair) in digest.iter_mut().zip(hex.chunks_exact(2)) {
        *byte = (hex_nibble(pair[0])? << 4) | hex_nibble(pair[1])?;
    }
    Some(digest)
}

fn hex_nibble(character: u8) -> Option<u8> {
    match character {
        b'0'..=b'9' => Some(character - b'0'),
        b'a'..=b'f' => Some(character - b'a' + 10),
        b'A'..=b'F' => Some(character - b'A' + 10),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::test_fixtures::{
        BIG_FIXTURE_CONTIG, FIXTURE_CONTIGS, bam_header, big_contig_specs, big_fixture_reference,
        big_spread_of_reads, fixture_read_group, fixture_reference, header, indexed_bam,
        matching_contigs, one_read, only_tally, unindexed_bam,
    };
    // `GenomeRegion` left `super::*` when the per-region query did: the file layer's own API
    // no longer names one, and only the tests below build them.
    use crate::ng::types::GenomeRegion;

    // --- @HD SO ---

    #[test]
    fn sort_order_reads_the_tag_when_present() {
        assert_eq!(
            sort_order(&header(Some("coordinate"), &[], &[])).as_deref(),
            Some("coordinate")
        );
        assert_eq!(
            sort_order(&header(Some("queryname"), &[], &[])).as_deref(),
            Some("queryname"),
            "a wrong value is reported, not silently normalised — the gate \
             wants to name what it found"
        );
    }

    #[test]
    fn sort_order_is_none_when_the_tag_or_the_hd_line_is_absent() {
        assert_eq!(sort_order(&header(None, &[], &[])), None, "@HD with no SO");
        assert_eq!(
            sort_order(&sam::Header::default()),
            None,
            "no @HD line at all"
        );
    }

    // --- @SQ ---

    #[test]
    fn contig_list_preserves_header_order_with_names_and_lengths() {
        let contigs = contig_list(&header(
            Some("coordinate"),
            &[("chr2", 200, None), ("chr1", 100, None)],
            &[],
        ))
        .expect("no M5 tags to be malformed");

        let names: Vec<&str> = contigs.entries.iter().map(|e| e.name.as_str()).collect();
        assert_eq!(
            names,
            vec!["chr2", "chr1"],
            "header order is kept as-is — comparing it against the reference \
             *in order* is what catches a permutation"
        );
        assert_eq!(contigs.entries[0].length, 200);
        assert_eq!(contigs.entries[1].length, 100);
    }

    #[test]
    fn contig_list_decodes_an_m5_tag_and_tolerates_its_absence() {
        let contigs = contig_list(&header(
            Some("coordinate"),
            &[
                ("chr1", 100, Some("0123456789abcdef0123456789ABCDEF")),
                ("chr2", 200, None),
            ],
            &[],
        ))
        .expect("both tags are well formed");

        assert_eq!(
            contigs.entries[0].md5,
            Some([
                0x01, 0x23, 0x45, 0x67, 0x89, 0xab, 0xcd, 0xef, 0x01, 0x23, 0x45, 0x67, 0x89, 0xab,
                0xcd, 0xef
            ]),
            "decoded from hex, and upper case is accepted alongside lower"
        );
        assert_eq!(
            contigs.entries[1].md5, None,
            "a contig without M5 is ordinary input, not a fault"
        );
    }

    /// A malformed digest is a **hard error**, not a silently-absent one.
    /// Reading it as absent would pass a real defect silently — and a file
    /// whose tags were all malformed would look untagged, switching
    /// wrong-assembly detection off with nobody told.
    #[test]
    fn a_malformed_m5_is_rejected_naming_the_contig_and_the_reason() {
        let wrong_length = contig_list(&header(
            Some("coordinate"),
            &[("chr1", 100, Some("abc"))],
            &[],
        ))
        .expect_err("a three-character digest is malformed");
        assert_eq!(wrong_length.contig, "chr1");
        assert_eq!(wrong_length.detail, "expected 32 hex characters, got 3");

        let not_hex = contig_list(&header(
            Some("coordinate"),
            &[("chr2", 200, Some("zzzz56789abcdef0123456789abcdef0"))],
            &[],
        ))
        .expect_err("a non-hex character is malformed");
        assert_eq!(not_hex.contig, "chr2");
        assert_eq!(not_hex.detail, "contains a non-hex character");

        for bad in ["", &"a".repeat(33)] {
            assert!(
                contig_list(&header(
                    Some("coordinate"),
                    &[("chr1", 100, Some(bad))],
                    &[]
                ))
                .is_err(),
                "malformed M5 {bad:?} must be rejected"
            );
        }
    }

    /// The same fault through the gate, so the path and the phrasing are
    /// covered too — and so it is unmistakable that this rejects the *file*,
    /// not merely the tag.
    #[test]
    fn the_gate_rejects_a_file_with_a_malformed_m5() {
        let contigs = vec![("chr1", 100, Some("not-a-digest")), ("chr2", 200, None)];

        let error = open_fixture(&contigs, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::MalformedMd5 { contig, detail, .. } => {
                assert_eq!(contig, "chr1");
                assert_eq!(detail, "expected 32 hex characters, got 12");
            }
            other => panic!("expected MalformedMd5, got {other:?}"),
        }
        assert!(
            error
                .to_string()
                .contains("@SQ 'chr1' has a malformed M5 tag"),
            "{error}"
        );
    }

    #[test]
    fn contig_list_is_empty_for_a_header_with_no_reference_sequences() {
        // B2 compares this against the reference, so an empty list is rejected
        // on length rather than on anything subtler.
        assert!(
            contig_list(&header(Some("coordinate"), &[], &[]))
                .expect("no tags at all")
                .entries
                .is_empty()
        );
    }

    // -----------------------------------------------------------------
    // The chromosomes a cursor may be made for (alignment cursor, A2)
    // -----------------------------------------------------------------

    /// The list carries the reference's names, lengths and **order** — the part the open
    /// gate proves equal, and the part that makes `ContigId(i)` mean the same thing to the
    /// file and to the reference. That is what `contigs()` exists to offer: a caller minting
    /// one cursor per chromosome walks this without holding the reference open beside it.
    ///
    /// **Compared field by field, not with `assert_eq!` on the lists.** `ContigEntry`'s
    /// `PartialEq` treats an absent `M5` as a wildcard, so comparing whole lists passes
    /// whether the field holds the file's or the reference's — a review mutation swapping one
    /// for the other left exactly that assertion green. Naming the three fields that really
    /// are proved equal is what makes this test able to fail.
    #[test]
    fn the_contig_list_has_the_reference_names_and_lengths_in_the_reference_order() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("opens");

        let expected: Vec<(&str, u64)> = matching_contigs()
            .iter()
            .map(|(name, length, _)| (*name, *length as u64))
            .collect();
        assert!(
            !expected.is_empty(),
            "an empty list would make the comparison below vacuous",
        );

        let observed: Vec<(&str, u64)> = file
            .contigs()
            .entries
            .iter()
            .map(|entry| (entry.name.as_str(), entry.length))
            .collect();
        assert_eq!(observed, expected);
    }

    /// `contigs()` and `sq_md5s()` are the list and a projection of it, so they cannot
    /// disagree — not about how many contigs there are, and not about any one contig's
    /// digest. Both are indexed by `ContigId`, and a divergence would silently pair a file's
    /// contig with a *different* reference contig inside `check_assembly`.
    ///
    /// A **distinct digest per contig**, so a projection that returned the same digest for
    /// every entry, or the entries in another order, fails here. The duplication itself is
    /// recorded as a follow-up on the `sq_md5s` field; while it stands, this is what stops
    /// the two drifting.
    #[test]
    fn the_contig_list_and_the_digest_list_index_alike() {
        let digests = distinct_digests();
        let fixture = open_fixture(&contigs_declaring(&digests), false);
        let file = fixture
            .file
            .as_ref()
            .expect("opens against a .fai reference");

        assert_eq!(file.contigs().entries.len(), file.sq_md5s().len());
        for (entry, digest) in file.contigs().entries.iter().zip(file.sq_md5s()) {
            assert_eq!(&entry.md5, digest);
        }
        assert!(
            file.sq_md5s().iter().all(Option::is_some),
            "the fixture declared an M5 for every contig, so a `None` here is the projection \
             losing one",
        );
    }

    /// **What `contigs()` holds is the file's claim, not the reference's** — and this is the
    /// test that can tell the two apart. The reference is the `.fai` arm, which carries no
    /// digests at all, while the file declares one per contig; if the field were populated
    /// from the reference, as this module's doc comments claimed until the A2 review, every
    /// digest here would be `None`.
    #[test]
    fn the_contig_list_carries_the_digests_the_file_declared_not_the_reference_s() {
        let (reference_dir, reference) = fixture_reference(false);
        assert!(
            reference
                .info()
                .contig_list()
                .entries
                .iter()
                .all(|entry| entry.md5.is_none()),
            "the .fai reference arm is supposed to carry no digests; without that this test \
             proves nothing",
        );

        let digests = distinct_digests();
        let fixture = open_fixture(&contigs_declaring(&digests), false);
        let file = fixture
            .file
            .as_ref()
            .expect("opens against a .fai reference");

        assert!(
            file.contigs()
                .entries
                .iter()
                .all(|entry| entry.md5.is_some()),
            "the file declared an M5 for every contig, so `contigs()` cannot be holding the \
             reference's list",
        );
        drop(reference_dir);
    }

    /// The `Debug` line is what a panic message or a log carries, and its `contigs` field is
    /// a **count** — the one shape that looks plausible whatever number it holds. A review
    /// mutation adding 999 to it survived every other test here.
    #[test]
    fn the_debug_line_counts_the_contigs_the_file_actually_has() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("opens");

        let rendered = format!("{file:?}");
        assert!(
            rendered.contains(&format!("contigs: {}", file.contigs().entries.len())),
            "the Debug line does not report this file's {} contigs:\n{rendered}",
            file.contigs().entries.len(),
        );
    }

    /// One distinct 32-hex-digit `M5` per fixture contig, so a test can tell the digests
    /// apart rather than only telling present from absent.
    fn distinct_digests() -> Vec<String> {
        (1..=matching_contigs().len())
            .map(|nth| format!("{nth:032x}"))
            .collect()
    }

    /// The fixture's `@SQ` shape with a digest attached to each contig.
    fn contigs_declaring(digests: &[String]) -> Vec<(&str, usize, Option<&str>)> {
        matching_contigs()
            .into_iter()
            .zip(digests)
            .map(|((name, length, _), digest)| (name, length, Some(digest.as_str())))
            .collect()
    }

    // -----------------------------------------------------------------
    // The cursor over a real BAM (alignment cursor, C3)
    // -----------------------------------------------------------------

    /// **The cursor tests run on a 200,000-base contig, and the reason is a defect this
    /// milestone shipped and a review caught.**
    ///
    /// The first version of these tests used the 100-base fixture contig. BAI's finest bins
    /// are 16 kb and a BGZF block is 64 kB, so *every* region there resolves to the same
    /// single chunk — and the run-of-regions oracle below, whose entire job is to catch a
    /// reader that stops at the previous region's end, **passed with exactly that defect in
    /// place**. An oracle that cannot fail for the thing the milestone is about is worse than
    /// no oracle, because it is counted as coverage.
    ///
    /// Reads are 30 bases because the default `min_read_length` silently drops anything
    /// shorter — a fixture that ignores it produces an empty walk, which reads as a broken
    /// cursor rather than as a badly-chosen read length.
    const BIG_STRIDE: usize = 90;

    fn big_spread() -> Vec<RecordBuf> {
        big_spread_of_reads(BIG_STRIDE)
    }

    /// **The oracle: a whole-file linear scan, filtered to the region.** It never seeks and
    /// never consults the index, so nothing about chunks, positioning or the early stop is
    /// shared with the thing it checks.
    fn names_by_linear_scan(path: &Path, region: GenomeRegion) -> Vec<String> {
        let mut reader = bam::io::reader::Builder
            .build_from_path(path)
            .expect("open bam");
        let header = reader.read_header().expect("read header");
        let mut names = Vec::new();
        for record in reader.record_bufs(&header) {
            let record = record.expect("a fixture record");
            let (Some(first), Some(last)) = (record.alignment_start(), record.alignment_end())
            else {
                continue;
            };
            if record.reference_sequence_id() == Some(region.contig.get() as usize)
                && usize::from(first) as u64 <= region.end.get()
                && usize::from(last) as u64 >= region.start.get()
            {
                names.push(
                    String::from_utf8_lossy(record.name().expect("named").as_ref()).into_owned(),
                );
            }
        }
        names
    }

    /// **The test the first attempt at this feature did not have.**
    ///
    /// Every other read-path test in this crate drives a *single* region query. That is
    /// precisely why 1,471 of them passed while 3,830 of 236,081 loci went missing:
    /// consecutive queries through one reader are the untested surface, and they are the
    /// entire feature. So a run of regions is asked of **one cursor**, in order, and each
    /// answer is compared against a scan of the whole file.
    #[test]
    fn a_run_of_regions_through_one_bam_cursor_matches_a_linear_scan() {
        let (reference_dir, reference) = big_fixture_reference();
        let (_bam_dir, path) = indexed_bam(&bam_header(&big_contig_specs()), &big_spread());
        let file = AlignmentFile::open(
            &path,
            &reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("the fixture opens");

        let mut cursor = file
            .cursor(ContigId(0), big_reference_bases())
            .expect("a cursor for contig 0");

        // Ascending, adjacent, overlapping, far apart, backward, repeated, and empty — the
        // shapes the plan names, all through the one cursor, and spanning enough of the
        // contig that the index has to resolve more than one chunk.
        for (start, end) in [
            (1u64, 5_000u64),
            (5_001, 10_000),
            (9_000, 20_000),
            (9_001, 20_001),
            (150_000, 160_000),
            (1, 8_000),
            (1, 8_000),
            (199_999, 200_000),
            (1, 200_000),
        ] {
            let asked = GenomeRegion {
                contig: ContigId(0),
                start: Position(start),
                end: Position(end),
            };
            cursor.move_to_region(asked).expect("on this chromosome");
            let mut actual = Vec::new();
            while let Some(read) = cursor.next_read() {
                let read = read.expect("the fixture reads decode and filter");
                actual.push(String::from_utf8_lossy(&read.qname).into_owned());
            }

            assert_eq!(
                actual,
                names_by_linear_scan(&path, asked),
                "the cursor disagreed with a linear scan at [{start}, {end}] — and it had \
                 already served every region before it",
            );
        }
        drop(reference_dir);
    }

    // -----------------------------------------------------------------
    // The cursor over a real CRAM (alignment cursor, E2)
    // -----------------------------------------------------------------

    /// A raw reference accessor over the CRAM fixture's FASTA, for the cursor's mismatch
    /// filter — the same shape the BAM cursor tests use, over the fixture the CRAM needs.
    fn cram_cursor_reference_bases(fasta: &Path) -> crate::ng::ref_seq::WindowedRefSeq {
        let info = read_reference_info(ReferenceSource::Fasta {
            fasta: fasta.to_path_buf(),
            fai: None,
        })
        .expect("read reference");
        crate::ng::ref_seq::WindowedRefSeq::new(fasta.to_path_buf(), info.contig_list())
    }

    /// The CRAM cursor fixture: a long contig with **more records than fit in one container**.
    ///
    /// noodles writes 10,240 records per container, so a fixture under that produces one
    /// container and one `.crai` entry — which exercises the container decode and **none** of
    /// the walk: not the multi-entry loop, not the positioning, not the carry-on past a
    /// container boundary. That is the CRAM shape of the trap the BAM tests above describe,
    /// and it would let the oracle pass with a reader that stops at the first container.
    const CRAM_CURSOR_CONTIG_LENGTH: usize = 400_000;
    const CRAM_CURSOR_READS: usize = 25_000;

    /// A whole-file linear scan of a CRAM, filtered to the region — the same oracle shape the
    /// BAM tests use, sharing nothing with the thing it checks: it never seeks and never opens
    /// the `.crai`.
    fn cram_names_by_linear_scan(path: &Path, fasta: &Path, region: GenomeRegion) -> Vec<String> {
        use noodles_fasta as fasta;

        let repository = fasta::Repository::new(
            fasta::io::indexed_reader::Builder::default()
                .build_from_path(fasta)
                .map(fasta::repository::adapters::IndexedReader::new)
                .expect("the fixture FASTA is indexed"),
        );
        let mut reader = cram::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(path)
            .expect("open cram");
        let header = reader.read_header().expect("read header");
        let mut names = Vec::new();
        for record in reader.records(&header) {
            let record = record.expect("a fixture record");
            let record = RecordBuf::try_from_alignment_record(&header, &record)
                .expect("a fixture record converts");
            let (Some(first), Some(last)) = (record.alignment_start(), record.alignment_end())
            else {
                continue;
            };
            if record.reference_sequence_id() == Some(region.contig.get() as usize)
                && usize::from(first) as u64 <= region.end.get()
                && usize::from(last) as u64 >= region.start.get()
            {
                names.push(
                    String::from_utf8_lossy(record.name().expect("named").as_ref()).into_owned(),
                );
            }
        }
        names
    }

    /// **The parity oracle, now through cursors — E2.**
    ///
    /// A run of regions is asked of **one CRAM cursor**, in order, and each answer is compared
    /// against a scan of the whole file. This is the CRAM half of what
    /// [`a_run_of_regions_through_one_bam_cursor_matches_a_linear_scan`] does, and it is the
    /// test the milestone exists to pass: the two formats share every line above the record
    /// reader, so a disagreement here is a container-walk bug.
    ///
    /// The regions deliberately cross container boundaries, go backwards, repeat, and end with
    /// the whole contig — because a reader that positioned correctly but *bounded* at the
    /// region's end would answer the first region right and every later one short.
    #[test]
    fn a_run_of_regions_through_one_cram_cursor_matches_a_linear_scan() {
        let (_cram_dir, path, _fasta_dir, fasta) =
            multi_container_cram(CRAM_CURSOR_CONTIG_LENGTH, CRAM_CURSOR_READS);
        let reference = OpenReference::from(
            read_reference_info(ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            })
            .expect("read reference"),
        );
        let file = AlignmentFile::open(
            &path,
            &reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("the fixture opens");

        let mut cursor = file
            .cursor(ContigId(0), cram_cursor_reference_bases(&fasta))
            .expect("a cursor for contig 0");

        // Ascending, adjacent, overlapping, far apart, backward, repeated, and the whole
        // contig — plus **a region just inside each container's first base**, taken from the
        // `.crai` itself rather than hard-coded so it stays a boundary whatever the writer's
        // container size turns out to be.
        //
        // Those boundary regions are what makes this able to fail for a reader that positions
        // at the container *starting* at or after the region and forgets the one before it: a
        // read beginning a few bases earlier reaches into the region and lives in the previous
        // container. Without them the fixture's regions all sit comfortably inside a container
        // and a reader that never stepped back would pass.
        let mut regions: Vec<(u64, u64)> = vec![
            (1, 20_000),
            (20_001, 40_000),
            (30_000, 90_000),
            (30_001, 90_001),
            (300_000, 340_000),
            (1, 25_000),
            (1, 25_000),
            (399_000, 400_000),
        ];
        let crai = cram::crai::fs::read(format!("{}.crai", path.display()))
            .expect("the fixture writes a .crai");
        for entry in &crai {
            let Some(first) = entry.alignment_start() else {
                continue;
            };
            let first = usize::from(first) as u64;
            if first > 1 {
                regions.push((first + 4, first + 60));
                regions.push((first, first + 60));
            }
        }
        regions.push((1, 400_000));

        for (start, end) in regions {
            let asked = GenomeRegion {
                contig: ContigId(0),
                start: Position(start),
                end: Position(end),
            };
            cursor.move_to_region(asked).expect("on this chromosome");
            let mut actual = Vec::new();
            while let Some(read) = cursor.next_read() {
                let read = read.expect("the fixture reads decode and filter");
                actual.push(String::from_utf8_lossy(&read.qname).into_owned());
            }

            assert_eq!(
                actual,
                cram_names_by_linear_scan(&path, &fasta, asked),
                "the CRAM cursor disagreed with a linear scan at [{start}, {end}] — and it \
                 had already served every region before it",
            );
        }
    }

    /// The CRAM oracle must be able to fail, and on a fixture that reaches past one container.
    #[test]
    fn the_cram_cursor_oracle_is_not_vacuous() {
        let (_cram_dir, path, _fasta_dir, fasta) =
            multi_container_cram(CRAM_CURSOR_CONTIG_LENGTH, CRAM_CURSOR_READS);

        let covered = cram_names_by_linear_scan(
            &path,
            &fasta,
            GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(CRAM_CURSOR_CONTIG_LENGTH as u64),
            },
        );
        assert_eq!(
            covered.len(),
            CRAM_CURSOR_READS,
            "the whole contig holds every fixture read",
        );

        let crai = cram::crai::fs::read(format!("{}.crai", path.display()))
            .expect("the fixture writes a .crai");
        assert!(
            crai.len() > 1,
            "the fixture must span more than one container, or the walk this checks is one \
             decode and the oracle cannot fail for a reader that stops there: {} entries",
            crai.len(),
        );
    }

    /// The oracle must be able to fail: a region the fixture covers has to return reads, or
    /// the test above could be comparing two empty vectors and calling it agreement.
    #[test]
    fn the_bam_cursor_oracle_is_not_vacuous() {
        let (_bam_dir, path) = indexed_bam(&bam_header(&big_contig_specs()), &big_spread());

        let covered = names_by_linear_scan(
            &path,
            GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(200_000),
            },
        );
        assert_eq!(
            covered.len(),
            big_spread().len(),
            "the whole contig holds every fixture read",
        );

        let narrow = names_by_linear_scan(
            &path,
            GenomeRegion {
                contig: ContigId(0),
                // Inside the read starting at 91, which is 30 bases long. The reads are
                // 90 apart, so most single bases fall in a gap.
                start: Position(100),
                end: Position(101),
            },
        );
        assert!(
            !narrow.is_empty() && narrow.len() < covered.len(),
            "and the regions must not all return the same thing: {narrow:?}",
        );
    }

    /// **A forward walk over a real BAM decodes each read once**, which is the claim the whole
    /// feature rests on and the number the perf review says has to come down.
    #[test]
    fn a_forward_walk_over_a_bam_decodes_each_read_once() {
        let (reference_dir, reference) = big_fixture_reference();
        let (_bam_dir, path) = indexed_bam(&bam_header(&big_contig_specs()), &big_spread());
        let file = AlignmentFile::open(
            &path,
            &reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("the fixture opens");

        let mut cursor = file
            .cursor(ContigId(0), big_reference_bases())
            .expect("a cursor for contig 0");

        // Overlapping, ascending — what a real region walk looks like.
        for (start, end) in [
            (1u64, 60_000u64),
            (50_000, 110_000),
            (100_000, 160_000),
            (150_000, 200_000),
        ] {
            let asked = GenomeRegion {
                contig: ContigId(0),
                start: Position(start),
                end: Position(end),
            };
            cursor.move_to_region(asked).expect("on this chromosome");
            while let Some(read) = cursor.next_read() {
                read.expect("the fixture reads decode and filter");
            }
        }

        let counts = cursor.counts();
        assert_eq!(
            counts.reads_decoded,
            big_spread().len() as u64,
            "every read in the file, and a forward walk must decode each of them once",
        );
        assert!(
            counts.reads_replayed > 0,
            "reads shared by consecutive regions must be served from what is held",
        );
        assert_eq!(
            counts.regions_jumping, 1,
            "only the first region has nothing to reuse"
        );
        drop(reference_dir);
    }

    // -----------------------------------------------------------------
    // The composed chain, through a cursor — T9, T10, and what the reader
    // pool's tests said before Milestone F deleted the pool
    // -----------------------------------------------------------------

    use noodles_sam::alignment::RecordBuf;

    use crate::ng::read::input::test_fixtures::read_named_with_length;
    use crate::ng::ref_seq::InMemoryRefSeq;
    use crate::ng::types::Position;

    /// All-`A` bases over the big fixture contig, so nothing is dropped for mismatching.
    fn big_reference_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![vec![b'A'; BIG_FIXTURE_CONTIG.1]])
    }

    /// An all-`A` reference matching the fixture contigs, so a read of `A`s
    /// matches perfectly and a read of `C`s mismatches at every base.
    fn reference_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        )
    }

    /// A read whose bases all mismatch the all-`A` reference — filter #8's
    /// business, and only reachable if the *whole* filter is composed in.
    fn all_mismatching_read(qname: &str, start: usize) -> RecordBuf {
        let mut record = read_named_with_length(qname, 0, start, 30);
        record.sequence_mut().as_mut().fill(b'C');
        record
    }

    /// The file comes back behind an `Arc` because that is what querying one
    /// now takes: a region stream shares ownership of the file it reads
    /// through, so it can hand the pooled reader back on `Drop`.
    fn opened_over(records: &[RecordBuf]) -> (TempDir, TempDir, Arc<AlignmentFile>) {
        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, path) = indexed_bam(&bam_header(&matching_contigs()), records);
        let file = AlignmentFile::open(
            &path,
            &reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("the fixture matches");
        (reference_dir, bam_dir, file)
    }

    fn whole_first_contig() -> GenomeRegion {
        GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(1),
            end: crate::ng::types::Position(FIXTURE_CONTIGS[0].1 as u64),
        }
    }

    /// Every read a cursor yields for one region of one file, by name — the whole real chain
    /// (aligned-reads reader → region narrowing → step-1 filter → order guard).
    fn cursor_names(file: &Arc<AlignmentFile>, region: GenomeRegion) -> Vec<String> {
        let mut cursor = file
            .cursor(region.contig, reference_bases())
            .expect("a cursor for this contig");
        cursor.move_to_region(region).expect("on this chromosome");
        let mut names = Vec::new();
        while let Some(read) = cursor.next_read() {
            let read = read.expect("no fatal read error");
            names.push(String::from_utf8_lossy(&read.qname).into_owned());
        }
        names
    }

    /// **T9 — the full step-1 filter runs, not just the cheap subset.**
    ///
    /// A read that mismatches the reference at every base is filter #8's business, and #8 is
    /// the *reference-dependent* filter — the one a reader applying only flag/MAPQ checks
    /// would miss. Its drop being charged to `high_mismatch_fraction` is what proves
    /// `ReadFilter` is composed into the chain rather than a subset of it. This is the
    /// property that decided the rebuild over reusing production's reader.
    ///
    /// Converted from the per-region query at Milestone F. The tally is read off the cursor,
    /// which is where it lives now: the filter lasts as long as the cursor, so nothing has to
    /// be folded back into the file when a region ends.
    #[test]
    fn t9_a_cursor_runs_the_reference_dependent_filter() {
        let (_reference_dir, _bam_dir, file) = opened_over(&[
            read_named_with_length("clean", 0, 1, 30),
            all_mismatching_read("mismatching", 40),
        ]);

        let mut cursor = file
            .cursor(ContigId(0), reference_bases())
            .expect("a cursor for contig 0");
        cursor
            .move_to_region(whole_first_contig())
            .expect("on this chromosome");
        let mut reads = Vec::new();
        while let Some(read) = cursor.next_read() {
            reads.push(read.expect("no fatal error"));
        }

        assert_eq!(reads.len(), 1, "only the clean read survives");
        assert_eq!(reads[0].qname, b"clean");

        let counts = only_tally(&cursor.read_group_counts());
        assert_eq!(
            counts.high_mismatch_fraction, 1,
            "the drop is charged to filter #8, so the whole cascade ran"
        );
        assert_eq!(counts.kept, 1);
    }

    /// **What the reader over-returned is dropped uncounted.** A read outside the region, and
    /// a record with no footprint at all, are not reads the filter rejected — they are reads
    /// the reader was never asked about. Charging them to a `DropReason` would make the tally
    /// mean something different for a narrowed read than for a whole-file one.
    ///
    /// The footprint-less record is here rather than only in `region_raw_aligned_reads` because
    /// this is where the *tally* is observable: that layer can say the record never surfaced,
    /// and only this one can say it was charged to nothing.
    #[test]
    fn records_outside_the_region_are_dropped_without_being_counted() {
        let (_reference_dir, _bam_dir, file) = opened_over(&[
            read_named_with_length("inside", 0, 1, 30),
            read_named_with_length("outside", 0, 60, 30),
        ]);

        let region = GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(35),
        };
        let mut cursor = file
            .cursor(ContigId(0), reference_bases())
            .expect("a cursor for contig 0");
        cursor.move_to_region(region).expect("on this chromosome");
        let mut reads = Vec::new();
        while let Some(read) = cursor.next_read() {
            reads.push(read.expect("no fatal error"));
        }

        assert_eq!(reads.len(), 1);
        let counts = only_tally(&cursor.read_group_counts());
        assert_eq!(counts.kept, 1);
        assert_eq!(
            counts.duplicate
                + counts.low_mapq
                + counts.supplementary
                + counts.secondary
                + counts.unmapped
                + counts.qc_fail
                + counts.too_short
                + counts.high_mismatch_fraction
                + counts.bad_cigar,
            0,
            "the out-of-region read is charged to no drop reason at all"
        );
    }

    /// **T10 — a fatal mid-stream error is yielded once, and then the cursor refuses.** A
    /// truncated file must not look like a short region.
    ///
    /// The per-region query fused its iterator and the next query started afresh; a cursor
    /// goes further, because it outlives the region — it **refuses every later region**
    /// rather than answering them from what it happens to be holding, which would be a
    /// plausible, silently short answer for the rest of the chromosome.
    #[test]
    fn t10_a_truncated_file_fails_once_and_then_refuses_later_regions() {
        let (_reference_dir, reference) = fixture_reference(false);
        // Enough records to span several BGZF blocks: truncating a single-block file would
        // destroy the *header*, which is a different fault (the gate's) from the mid-stream
        // one this test is about.
        let mut records: Vec<RecordBuf> = Vec::new();
        for start in 1..=60 {
            for copy in 0..120 {
                records.push(read_named_with_length(
                    &format!("r{start}_{copy}"),
                    0,
                    start,
                    30,
                ));
            }
        }
        let (_bam_dir, path) = indexed_bam(&bam_header(&matching_contigs()), &records);

        let file = AlignmentFile::open(
            &path,
            &reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("opens");

        // The cursor is minted before the truncation, as a real run's would be: the gate and
        // the index are intact and the fault can only appear mid-stream.
        let mut cursor = file
            .cursor(ContigId(0), reference_bases())
            .expect("a cursor for contig 0");

        let full = std::fs::metadata(&path).expect("stat").len();
        std::fs::OpenOptions::new()
            .write(true)
            .open(&path)
            .expect("reopen")
            .set_len(full * 3 / 4)
            .expect("truncate");

        cursor
            .move_to_region(whole_first_contig())
            .expect("the move itself still succeeds");

        let mut reads_before_error = 0;
        let mut errors = 0;
        let mut after_error = 0;
        while let Some(item) = cursor.next_read() {
            match item {
                Ok(_) => reads_before_error += 1,
                Err(_) => {
                    errors += 1;
                    while cursor.next_read().is_some() {
                        after_error += 1;
                    }
                    break;
                }
            }
        }

        // Without this the test would pass against a chain that yielded *nothing* — it would
        // still reach the truncation and still stop. That is exactly the shape of the fixture
        // bug this module already hit once.
        assert!(
            reads_before_error > 0,
            "reads must flow before the truncation is reached"
        );
        assert_eq!(errors, 1, "the truncation surfaced as a fatal error");
        assert_eq!(after_error, 0, "and the walk stopped rather than resuming");

        assert!(
            matches!(
                cursor.move_to_region(whole_first_contig()),
                Err(crate::ng::read::input::cursor::CursorError::AfterFailure { .. })
            ),
            "a cursor whose file failed must refuse later regions, not answer them short",
        );
    }

    /// **A cursor is `Send`, so a worker can own one.** The parallel fan-out is one cursor per
    /// worker sharing nothing, so an `Rc` or a `RefCell` creeping into the chain would break
    /// that — and would surface at the first parallel call site, a plan or two away, rather
    /// than here.
    #[test]
    fn a_cursor_is_send_so_a_worker_can_own_one() {
        fn assert_send<T: Send>() {}
        assert_send::<AlignmentCursor<InMemoryRefSeq>>();
        assert_send::<&AlignmentFile>();
    }

    /// And the same thing dynamically: eight threads, one open file, a cursor each, every one
    /// of them reading the whole contig correctly.
    ///
    /// The reader pool this replaces made concurrency a property of the *file*; it is now a
    /// property of there being one cursor per worker, each with its own descriptor. What can
    /// still fail is any sharing that crept in below — which is what this can catch.
    #[test]
    fn cursors_on_one_file_read_the_same_thing_from_many_threads() {
        let (_reference_dir, _bam_dir, file) = opened_over(&[
            read_named_with_length("a", 0, 1, 30),
            read_named_with_length("b", 0, 40, 30),
        ]);

        std::thread::scope(|scope| {
            for _ in 0..8 {
                scope.spawn(|| {
                    assert_eq!(
                        cursor_names(&file, whole_first_contig()),
                        vec!["a".to_string(), "b".to_string()],
                    );
                });
            }
        });
    }

    /// **A cursor for a chromosome the file does not have is refused, and refused here.**
    ///
    /// Every layer below would answer an unknown contig as "nothing here", which is
    /// indistinguishable from a chromosome that is genuinely empty — so a caller looping over
    /// a stale contig list would silently produce no calls for it. Inherited from the region
    /// planners the per-region query had, which raised the same refusal per query.
    #[test]
    fn a_cursor_for_a_contig_the_file_does_not_have_is_an_error() {
        let (_reference_dir, _bam_dir, file) =
            opened_over(&[read_named_with_length("r", 0, 1, 30)]);

        let error = file
            .cursor(ContigId(9), reference_bases())
            .err()
            .expect("a contig the file does not declare must be refused");
        match error {
            AlignmentFileError::CursorContigNotInFile {
                contig,
                contigs_in_file,
                ..
            } => {
                assert_eq!(contig, ContigId(9));
                assert_eq!(contigs_in_file, FIXTURE_CONTIGS.len());
            }
            other => panic!("expected CursorContigNotInFile, got {other:?}"),
        }

        // And the contigs it *does* have are fine, so the check is a check and not a refusal
        // of everything.
        assert!(file.cursor(ContigId(0), reference_bases()).is_ok());
    }

    // -----------------------------------------------------------------
    // BAM/CRAM parity (C5) — T8
    // -----------------------------------------------------------------

    use crate::ng::read::input::test_fixtures::{indexed_cram, multi_container_cram};

    /// A spread of reads with pile-ups, so parity means more than "both
    /// returned nothing".
    ///
    /// One contig, because the CRAM fixture cannot hold more — see
    /// `test_fixtures::indexed_cram` for why (a noodles indexing limitation,
    /// not a choice).
    fn parity_records() -> Vec<RecordBuf> {
        let mut records = Vec::new();
        let mut start = 1;
        while start + 30 < FIXTURE_CONTIGS[0].1 {
            records.push(read_named_with_length(&format!("r{start}"), 0, start, 30));
            if start % 20 == 1 {
                records.push(read_named_with_length(&format!("r{start}_b"), 0, start, 30));
            }
            start += 7;
        }
        records
    }

    /// One region's read names through a cursor of its own — the shape the parity comparison
    /// needs, since the point is that two *files* agree rather than that one cursor's run of
    /// regions does (which `a_run_of_regions_through_one_bam_cursor_matches_a_linear_scan`
    /// covers).
    fn names_from(file: &Arc<AlignmentFile>, region: GenomeRegion) -> Vec<String> {
        cursor_names(file, region)
    }

    /// **T8 — the same reads written as BAM and as CRAM produce the same
    /// ordered stream.**
    ///
    /// The two containers share nothing below the aligned-reads-reader seam: BAM reads one record at
    /// a time from bgzf chunks, CRAM decodes whole containers against the reference and walks
    /// a `.crai`. Everything above — the region narrowing, the filter, the order guard — is
    /// the same code. So a disagreement here is a aligned-reads-reader bug, which is exactly what
    /// this is looking for, and BAM is the oracle because it is the simpler reader and was
    /// verified first against a linear scan.
    ///
    /// **Through cursors since Milestone F**, which is the point of converting it rather than
    /// deleting it: the two arms it compares are the two that survive.
    #[test]
    fn t8_a_cram_yields_the_same_ordered_reads_as_the_same_bam() {
        let records = parity_records();

        let (_reference_dir, bam_reference) = fixture_reference(false);
        let (_bam_dir, bam_path) = indexed_bam(&bam_header(&matching_contigs()), &records);
        let bam = AlignmentFile::open(
            &bam_path,
            &bam_reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("the BAM opens");

        let (_cram_dir, cram_path, _fasta_dir, fasta) = indexed_cram(&records);
        let cram_reference = OpenReference::from(
            read_reference_info(ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            })
            .expect("read reference"),
        );
        let cram = AlignmentFile::open(
            &cram_path,
            &cram_reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("the CRAM opens");

        let regions = [
            (0u32, 1u64, 100u64),
            (0, 1, 30),
            (0, 45, 55),
            (0, 60, 62),
            (0, 90, 100),
        ];

        let mut total = 0;
        for (contig, start, end) in regions {
            let region = GenomeRegion {
                contig: crate::ng::types::ContigId(contig),
                start: crate::ng::types::Position(start),
                end: crate::ng::types::Position(end),
            };

            let from_bam = names_from(&bam, region);
            let from_cram = names_from(&cram, region);

            assert_eq!(
                from_cram, from_bam,
                "BAM and CRAM disagreed for contig {contig} [{start}, {end}]"
            );
            total += from_bam.len();
        }

        assert!(
            total > 20,
            "the fixture must actually cover these regions — {total} reads is \
             too few for the comparison to mean anything"
        );
    }

    /// **The owner's rule, enforced at open.** A CRAM stores sequences as
    /// differences from the reference, so it cannot be decoded from a `.fai`,
    /// which holds only geometry. That is a hard error the moment the file is
    /// opened — not a mystery at the first query, and not a silently empty
    /// stream.
    #[test]
    fn a_cram_against_a_fai_only_reference_is_refused_at_open() {
        let (_cram_dir, cram_path, _fasta_dir, fasta) = indexed_cram(&one_read());

        let fai_only = OpenReference::from(
            read_reference_info(ReferenceSource::Fai(
                crate::ng::reference_info::sibling_fai_path(&fasta),
            ))
            .expect("read reference"),
        );

        let error = AlignmentFile::open(
            &cram_path,
            &fai_only,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect_err("a CRAM needs the bases");
        assert!(matches!(
            error,
            AlignmentFileError::CramNeedsReferenceFasta { .. }
        ));
        assert!(
            error.to_string().contains("supply the reference FASTA"),
            "the message must say what to do about it: {error}"
        );
    }

    /// A BAM is unaffected — it stores its own sequences, so a `.fai`-only
    /// reference is perfectly ordinary input for it.
    #[test]
    fn a_bam_against_a_fai_only_reference_opens_normally() {
        let (_reference_dir, _bam_dir, file) =
            opened_over(&[read_named_with_length("r", 0, 1, 30)]);
        assert_eq!(file.read_group_resolution(), &fixture_read_group());
    }

    // **`a_multi_container_cram_walks_its_crai_and_stops_early` died with the per-region
    // query, and both halves of it are accounted for.** It asserted two things: that a walk
    // reaches later containers rather than only ever decoding the first, and that it *stops
    // early* rather than decoding the contig to the end.
    //
    // The first is subsumed — `a_run_of_regions_through_one_cram_cursor_matches_a_linear_scan`
    // runs a whole sequence of regions, including late ones and the whole contig, over a
    // three-container fixture, and compares every answer against a scan of the file.
    //
    // The second is a rule the cursor **deliberately does not have**: an aligned-reads reader
    // positions, it never bounds (`aligned_reads_reader/mod.rs`), because the cursor serves a
    // forward region by not repositioning at all. Where the walk *starts* is still asserted —
    // by `containers_ending_before_the_region_are_skipped_rather_than_walked`, in
    // `aligned_reads_reader::cram` — and where it stops is now the layer above's business.

    /// **T2b's sequencing half.** The digest check is deferred *on purpose*:
    /// the reference a file is opened against usually carries no digests yet
    /// (they arrive from a background genome pass), so waiting for them would
    /// block startup on a whole-genome read.
    ///
    /// The file here is aligned to the **wrong assembly** — its `@SQ M5` tags
    /// disagree with the reference's real digests — and the whole point is that
    /// it opens anyway, streams its reads anyway, and is caught only at the end.
    ///
    /// Note what this can and cannot pin. Moving the comparison *into* `open`
    /// would not fail this test — and that is not a gap, it is the argument
    /// for deferring: `open` receives the `.fai`-arm reference, which carries
    /// no digests, so a check there would compare nothing and let the wrong
    /// assembly through regardless of where it sits. What is pinned is the
    /// sequence — the file opens, reads flow, and the fault is still caught —
    /// and that the catching is real (blinding `check_assembly`'s comparison
    /// fails this test).
    #[test]
    fn t2b_the_assembly_check_runs_after_the_reads_have_flowed() {
        use crate::ng::read::input::check_assembly;

        let wrong_digest = "ffffffffffffffffffffffffffffffff";
        let contigs: Vec<(&str, usize, Option<&str>)> = FIXTURE_CONTIGS
            .iter()
            .map(|(name, length)| (*name, *length, Some(wrong_digest)))
            .collect();

        // Opened against a `.fai`-only reference, which carries no digests — so
        // the gate's own comparison is a wildcard and the wrong assembly sails
        // through.
        let (_reference_dir, reference) = fixture_reference(false);
        let (_bam_dir, path) = indexed_bam(
            &bam_header(&contigs),
            &[read_named_with_length("r", 0, 1, 30)],
        );
        let file = AlignmentFile::open(
            &path,
            &reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        )
        .expect("a wrong M5 is a wildcard against a .fai-only reference");

        assert_eq!(
            cursor_names(&file, whole_first_contig()),
            vec!["r".to_string()],
            "the reads flowed first — startup never blocked"
        );

        // Only now, against a reference read through the FASTA arm and so
        // carrying real digests, does the fault surface.
        let (_verified_dir, verified) = fixture_reference(true);
        let error = check_assembly(file.path(), file.sq_md5s(), verified.info())
            .expect_err("the wrong assembly is caught, after the fact");
        assert_eq!(error.contig, "chr1");
        assert_eq!(error.observed, [0xff; 16]);
    }

    // -----------------------------------------------------------------
    // Grouping the `.crai` by contig — moved here at Milestone F
    // -----------------------------------------------------------------
    //
    // These three tests came from the per-region CRAM source, which reached the grouping
    // through a planner. The planner is gone; the grouping is not — `cursor` looks a
    // chromosome's entries up in it directly — so the tests are stated against
    // `group_crai_by_contig` itself, which was always their real subject.

    /// A `.crai` entry for `contig` covering `[start, start + span - 1]`, at `offset` — which
    /// the tests use to tell entries apart.
    fn crai_entry(
        contig: Option<usize>,
        start: usize,
        span: usize,
        offset: u64,
    ) -> cram::crai::Record {
        cram::crai::Record::new(
            contig,
            noodles_core::Position::new(start),
            span,
            offset,
            0,
            0,
        )
    }

    /// **The fix for production's O(n) prefix rescan** — and, more importantly, for an
    /// ordering assumption that does not hold. A query on a late contig must reach that
    /// contig's entries without walking (or mis-searching) the earlier ones.
    ///
    /// Hand-built rather than read from a file, because the file fixture cannot hold more than
    /// one contig (`test_fixtures::indexed_cram`).
    #[test]
    fn the_crai_is_grouped_so_each_contig_sees_only_its_own_entries() {
        let index = vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(Some(0), 51, 50, 200),
            crai_entry(Some(1), 1, 100, 300),
            crai_entry(Some(1), 101, 100, 400),
            crai_entry(Some(2), 1, 300, 500),
        ];
        let grouped = group_crai_by_contig(&index, 3);

        for (contig, expected_offsets) in [
            (0usize, vec![100u64, 200]),
            (1, vec![300, 400]),
            (2, vec![500]),
        ] {
            let offsets: Vec<u64> = grouped[contig].iter().map(|e| e.offset()).collect();
            assert_eq!(offsets, expected_offsets, "contig {contig}");
        }
    }

    /// **The ordering assumption a binary search would have made, violated by noodles
    /// itself.** Within a multi-reference slice, `fs::index` sorts by `Option<usize>` — and
    /// `None < Some(0)` — so the *unplaced* entry is emitted first. A `partition_point` over
    /// that is unspecified, and a walk would then meet a foreign entry and report
    /// end-of-input, losing every read of the region with no error at all.
    ///
    /// Grouping assumes nothing about order, so an interleaved index is fine.
    #[test]
    fn an_unplaced_entry_before_the_placed_ones_does_not_hide_a_contig() {
        let index = vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(None, 1, 0, 150),
            crai_entry(Some(1), 1, 100, 300),
        ];
        let grouped = group_crai_by_contig(&index, 3);

        let offsets: Vec<u64> = grouped[1].iter().map(|e| e.offset()).collect();
        assert_eq!(
            offsets,
            vec![300],
            "contig 1's entry must be found even though an unplaced entry precedes it — the \
             case that silently returned nothing before"
        );

        // And the unplaced entry belongs to no contig at all.
        for (contig, entries) in grouped.iter().enumerate() {
            assert!(
                entries.iter().all(|e| e.reference_sequence_id().is_some()),
                "an unplaced entry leaked into contig {contig}"
            );
        }
    }

    /// A contig with no entries yields an empty walk rather than someone else's entries.
    #[test]
    fn a_contig_absent_from_the_crai_has_an_empty_walk() {
        let index = vec![
            crai_entry(Some(0), 1, 50, 100),
            crai_entry(Some(2), 1, 300, 500),
        ];
        let grouped = group_crai_by_contig(&index, 3);

        assert!(grouped[1].is_empty());
    }

    // --- the hex decoder's edges ---

    /// The four malformed cases above put the bad character early or vary the
    /// length; none proves the **last** pair is inspected. A 32-character string
    /// that is valid hex except for its final character catches a decoder that
    /// stops one pair short — the off-by-one the `zip`/`chunks_exact` pairing
    /// could plausibly hide.
    #[test]
    fn decode_md5_hex_rejects_a_non_hex_character_in_the_final_position() {
        let mut hex = [b'0'; 32];
        assert!(decode_md5_hex(&hex).is_some(), "all-zeros is valid hex");

        hex[31] = b'z';
        assert!(decode_md5_hex(&hex).is_none(), "last nibble is inspected");

        let mut hex = [b'0'; 32];
        hex[30] = b'z';
        assert!(
            decode_md5_hex(&hex).is_none(),
            "high nibble of the last byte is inspected too"
        );
    }

    /// Range patterns fail at their edges, and a `z`-only test cannot see it.
    /// These six characters each sit one step outside a valid range.
    // -----------------------------------------------------------------
    // The gate — `AlignmentFile::open` (T1, T2, T3, T12a)
    // -----------------------------------------------------------------
    use tempfile::TempDir;

    use crate::ng::reference_info::{ReferenceSource, read_reference_info};

    /// An opened fixture **plus the temp dirs its files live in**.
    ///
    /// The dirs are returned rather than dropped at helper exit because they
    /// own the BAM and the reference on disk: drop them and the handle points
    /// at deleted files. Today's assertions would survive that (the index is
    /// already parsed in memory), but the first test that actually *queries*
    /// the file would fail somewhere far from the cause.
    struct OpenedFixture {
        file: Result<Arc<AlignmentFile>, AlignmentFileError>,
        _reference_dir: TempDir,
        _bam_dir: TempDir,
    }

    fn open_fixture(
        contigs: &[(&str, usize, Option<&str>)],
        reference_has_digests: bool,
    ) -> OpenedFixture {
        open_fixture_with_header(&bam_header(contigs), reference_has_digests)
    }

    fn open_fixture_with_header(
        header: &sam::Header,
        reference_has_digests: bool,
    ) -> OpenedFixture {
        let (reference_dir, reference) = fixture_reference(reference_has_digests);
        let (bam_dir, path) = indexed_bam(header, &one_read());
        OpenedFixture {
            file: AlignmentFile::open(
                &path,
                &reference,
                ReadFilterConfig::default(),
                false,
                fixture_read_group(),
            ),
            _reference_dir: reference_dir,
            _bam_dir: bam_dir,
        }
    }

    /// The `ContigReconcile` detail for a header that should fail check 2.
    fn reconcile_detail(
        contigs: &[(&str, usize, Option<&str>)],
        reference_has_digests: bool,
    ) -> String {
        match open_fixture(contigs, reference_has_digests)
            .file
            .expect_err("must not open")
        {
            AlignmentFileError::ContigReconcile { detail, .. } => detail,
            other => panic!("expected ContigReconcile, got {other:?}"),
        }
    }

    #[test]
    fn a_matching_file_opens_and_exposes_its_read_groups_and_digests() {
        let fixture = open_fixture(&matching_contigs(), false);
        let file = fixture.file.as_ref().expect("the file matches");

        // The sample is no longer the file's to know: the read-group table owns
        // it. What the file carries is how to read its records.
        assert_eq!(file.read_group_resolution(), &fixture_read_group());
        assert_eq!(
            file.sq_md5s().len(),
            2,
            "one slot per contig, indexed by ContigId"
        );
        assert!(
            file.sq_md5s().iter().all(Option::is_none),
            "this fixture's @SQ lines carry no M5, which is ordinary input"
        );
    }

    /// **T1 — the permutation hole, closed.**
    ///
    /// `@SQ` lists the reference's contigs in the wrong order. Every `ref_id`
    /// still *resolves*, which is why the check this replaces let it through
    /// and then fetched the wrong contig for every read. Order-aware equality
    /// catches it on the first transposed name.
    ///
    /// Mutation-verified below by
    /// `a_resolves_only_check_accepts_the_permutation_this_gate_rejects`.
    #[test]
    fn t1_a_permuted_sq_list_is_rejected_naming_the_first_transposed_contig() {
        let permuted = vec![("chr2", 200, None), ("chr1", 100, None)];

        let detail = reconcile_detail(&permuted, false);

        assert_eq!(
            detail, "name disagreement at index 0 ('chr2' vs 'chr1')",
            "the first transposed position, file value before reference value"
        );
    }

    /// **Why T1 needs an ordered comparison — the superseded probe, recorded.**
    ///
    /// `ReadFilter::new`'s check asks only "does every `@SQ` index resolve in
    /// the reference?", fetching each contig and accepting if the fetch
    /// succeeds. Run here against the permuted header, it **passes** — order is
    /// never consulted, so the file goes on to fetch the wrong contig for every
    /// read.
    ///
    /// **This is documentation, not the mutation guard.** It does not call
    /// `AlignmentFile::open`, so gutting the gate's comparison would leave it
    /// green; T1 is what catches that. The real mutation was performed against
    /// the gate itself during B2 — the `first_disagreement` call was replaced
    /// with a resolves-only length check and T1 duly failed — and this test is
    /// the standing record of *what* the old check did, so the reason the gate
    /// is stricter does not have to be taken on trust.
    #[test]
    fn the_superseded_resolves_only_probe_cannot_see_order() {
        use crate::ng::ref_seq::{InMemoryRefSeq, RefSeq};
        use crate::ng::types::ContigId;

        let (_reference_dir, reference) = fixture_reference(false);
        let permuted = bam_header(&[("chr2", 200, None), ("chr1", 100, None)]);

        // The superseded check: fetch every contig the file's `@SQ` list can
        // name, and accept if each one resolves. Order is never consulted.
        let reference_bases = InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        );
        let mut probe = Vec::new();
        let every_contig_resolves = (0..permuted.reference_sequences().len()).all(|index| {
            reference_bases
                .fetch_into(ContigId(index as u32), 1, 0, &mut probe)
                .is_ok()
        });
        assert!(
            every_contig_resolves,
            "every @SQ index of the permuted file resolves — so the old probe \
             accepts it, and then every read fetches the wrong contig"
        );

        // The gate's own comparison rejects that very same header.
        assert!(
            contig_list(&permuted)
                .expect("no M5 tags")
                .first_disagreement(&reference.info().contig_list())
                .is_err(),
            "the ordered comparison must reject what the probe accepted"
        );
    }

    /// **T2 — name, length and count mismatches**, each naming the right field
    /// and index.
    #[test]
    fn t2_name_length_and_count_mismatches_are_each_named() {
        let cases = [
            (
                vec![("chrX", 100, None), ("chr2", 200, None)],
                "name disagreement at index 0",
            ),
            (
                vec![("chr1", 100, None), ("chr2", 999, None)],
                "length disagreement at index 1",
            ),
            (
                vec![("chr1", 100, None)],
                "@SQ list length differs (1 vs 2)",
            ),
        ];

        for (contigs, expected) in cases {
            let detail = reconcile_detail(&contigs, false);
            assert!(
                detail.contains(expected),
                "expected {expected:?}, got {detail:?}"
            );
        }
    }

    /// **T2, the digest half.** A wrong `@SQ M5` is caught at open **only**
    /// when the `ReferenceInfo` in hand carries digests — the `Fasta` arm. Under
    /// a `.fai`-only table the comparison is a wildcard and the same file opens
    /// cleanly, which is the behaviour production ships and what makes the
    /// deferred `check_assembly` (D1) necessary rather than redundant.
    #[test]
    fn t2_a_wrong_m5_is_caught_only_when_the_reference_carries_digests() {
        let wrong_digest = "ffffffffffffffffffffffffffffffff";
        let contigs = vec![
            ("chr1", 100, Some(wrong_digest)),
            ("chr2", 200, Some(wrong_digest)),
        ];

        let detail = reconcile_detail(&contigs, true);
        assert!(
            detail.contains("md5 disagreement at index 0"),
            "got: {detail}"
        );

        assert!(
            open_fixture(&contigs, false).file.is_ok(),
            "against a .fai-only reference the digest is a wildcard, so the \
             same file opens — the gap check_assembly closes later"
        );
    }

    /// **T3 — `SO` wrong or missing**, rejected before anything else, and the
    /// message says which it was.
    #[test]
    fn t3_a_file_that_is_not_coordinate_sorted_is_rejected_at_open() {
        let contigs = matching_contigs();

        let queryname = header(Some("queryname"), &contigs, &[("rg1", Some("NA12878"))]);
        let error = open_fixture_with_header(&queryname, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::NotCoordinateSorted { sort_order, .. } => assert_eq!(
                sort_order.as_deref(),
                Some("queryname"),
                "carries the observed value as a fact, not pre-quoted prose"
            ),
            other => panic!("expected NotCoordinateSorted, got {other:?}"),
        }
        assert!(
            error.to_string().contains("@HD SO is 'queryname'"),
            "and the rendered message quotes it: {error}"
        );

        let no_sort_order = header(None, &contigs, &[("rg1", Some("NA12878"))]);
        let error = open_fixture_with_header(&no_sort_order, false)
            .file
            .expect_err("must not open");
        match &error {
            AlignmentFileError::NotCoordinateSorted { sort_order, .. } => {
                assert_eq!(sort_order.as_deref(), None)
            }
            other => panic!("expected NotCoordinateSorted, got {other:?}"),
        }
        assert!(error.to_string().contains("@HD SO is missing"), "{error}");
    }

    /// The gate's fail-fast order is spec §3.1's — `SO` → `@SQ` → index → `SM`
    /// — and it is only observable through a file that fails two checks at
    /// once. Each boundary gets its own case, so re-ordering any adjacent pair
    /// fails a test rather than quietly changing which fault a user is told
    /// about.
    #[test]
    fn the_four_checks_run_in_the_specified_order() {
        // SO before @SQ.
        let bad_sort_and_bad_contigs = header(
            Some("queryname"),
            &[("chrX", 1, None)],
            &[("rg1", Some("NA12878"))],
        );
        assert!(matches!(
            open_fixture_with_header(&bad_sort_and_bad_contigs, false).file,
            Err(AlignmentFileError::NotCoordinateSorted { .. })
        ));

        // @SQ before the index: bad contigs on a BAM with no index at all.
        let bad_contigs = bam_header(&[("chrX", 1, None)]);
        assert!(matches!(
            open_unindexed_with_header(&bad_contigs).1,
            Err(AlignmentFileError::ContigReconcile { .. })
        ));

        // The index before @RG SM: no index and no sample.
        let no_sample = header(Some("coordinate"), &matching_contigs(), &[]);
        assert!(matches!(
            open_unindexed_with_header(&no_sample).1,
            Err(AlignmentFileError::Index { .. })
        ));
    }

    /// Write a BAM with **no** index beside it and try to open it. Returns the
    /// temp dirs so the file outlives the call.
    fn open_unindexed_with_header(
        header: &sam::Header,
    ) -> (
        (TempDir, TempDir),
        Result<Arc<AlignmentFile>, AlignmentFileError>,
    ) {
        let (reference_dir, reference) = fixture_reference(false);
        let records = if header.reference_sequences().is_empty() {
            Vec::new()
        } else {
            one_read()
        };
        let (dir, path) = unindexed_bam(header, &records);

        let opened = AlignmentFile::open(
            &path,
            &reference,
            ReadFilterConfig::default(),
            false,
            fixture_read_group(),
        );
        ((reference_dir, dir), opened)
    }

    // The gate's fourth check — "exactly one `@RG SM`" — and its three tests are
    // gone. Two of the states it rejected are now rejected earlier and better,
    // by the read-group pre-pass, which names the file *and* the remedy and
    // tells "no `@RG` at all" apart from "an `@RG` with no `SM`" instead of
    // folding both into one variant. The third is no longer a fault at all: a
    // file whose read groups name two samples is ordinary input, and naming one
    // sample is a property of the open, enforced by `SampleReads`. The
    // replacements live in `read_groups.rs`.

    /// A missing index is an error, not a silent whole-file scan — and with the
    /// build flag it is repaired instead.
    #[test]
    fn a_missing_index_is_rejected_unless_the_caller_asks_for_one() {
        let header = bam_header(&matching_contigs());
        let (dirs, opened) = open_unindexed_with_header(&header);

        assert!(
            matches!(opened, Err(AlignmentFileError::Index { .. })),
            "an unindexed file is an error, never a silent whole-file scan"
        );

        // The same file in the same place — only the caller's policy differs.
        let (_reference_dir, bam_dir) = &dirs;
        let (_fresh_reference_dir, reference) = fixture_reference(false);
        let path = bam_dir.path().join("sample.bam");
        assert!(
            AlignmentFile::open(
                &path,
                &reference,
                ReadFilterConfig::default(),
                true,
                fixture_read_group(),
            )
            .is_ok(),
            "with build_index_if_missing the index is created next to the file"
        );
    }

    #[test]
    fn hex_nibble_rejects_the_characters_bordering_each_valid_range() {
        for character in [b'/', b':', b'`', b'g', b'@', b'G'] {
            assert_eq!(
                hex_nibble(character),
                None,
                "{:?} borders a valid range but is not hex",
                character as char
            );
        }
        assert_eq!(hex_nibble(b'0'), Some(0));
        assert_eq!(hex_nibble(b'9'), Some(9));
        assert_eq!(hex_nibble(b'a'), Some(10));
        assert_eq!(hex_nibble(b'f'), Some(15));
        assert_eq!(hex_nibble(b'A'), Some(10));
        assert_eq!(hex_nibble(b'F'), Some(15));
    }
}
