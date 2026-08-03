//! Step 1's **input edge**: turning alignment files on disk into ordered,
//! filtered streams of reads for a region.
//!
//! Two layers, split by subject. [`open_bam`], [`record_reader`],
//! [`region_records`] and [`cursor`] own the reader logic whose subject is
//! **one alignment file** — the validate-on-open gate, the per-format record
//! readers, the region narrowing and the long-lived cursor
//! (`doc/devel/ng/spec/alignment_file.md`, `alignment_cursor.md`).
//! [`sample_cursor`] and this module's `SampleReads` own what is true only of
//! **the sample** — the cross-file checks and the k-way merge
//! (`doc/devel/ng/spec/sample_reads.md`). Both layers' error types live here,
//! where the arch docs put them.
//!
//! **There is exactly one way to read a BAM or a CRAM here**, and that is
//! Milestone F's whole result: the per-region query, its reader pool and its
//! order-guard wrapper are gone, and a caller holds a cursor for a chromosome
//! and points it at each region.
//!
//! The layering is strict: the sample layer consumes what the file layer
//! produces and contains no reader logic of its own. It is also where the
//! module's two guarantees come from — every file's `ref_id` was proved equal
//! to its `ContigId` at open, and every cursor proves the reads it hands out
//! are coordinate-monotonic while streaming — which is what makes it sound for
//! the merge to compare positions *across* files without re-checking anything.
//!
//! This is step 1's input edge rather than a new step, so it lives under
//! `read/` beside `filtering.rs` and reuses that module's
//! [`RecordSource`](super::RecordSource)/[`RawRecord`](super::RawRecord) seam
//! (`doc/devel/ng/arch/module_layout.md` principle 1, note b).

pub mod cursor;
pub mod open_bam;
pub mod read_groups;
pub(crate) mod record_reader;
pub mod reference;
pub(crate) mod region_records;
pub mod sample_cursor;

use cursor::CursorError;
use sample_cursor::SampleCursor;

#[cfg(test)]
pub(crate) mod test_fixtures;

use std::io;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::bam::errors::AlignmentIndexError;
use crate::ng::read::filtering::ReadFilterConfig;
use crate::ng::read::input::read_groups::{
    ReadGroup, ReadGroupError, ReadGroupResolution, ReadGroups, RecordOwner, SampleReadGroups,
    build_read_groups,
};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::types::{ContigId, GenomeRegion, ReadGroupId};
use crate::pop_var_caller::common::format_md5_hex;

use crate::ng::ref_seq::RawRefSeq;

use open_bam::AlignmentFile;

/// Everything that can go wrong for **one** alignment file.
///
/// **Every variant fires before a read flows** — at `AlignmentFile::open`, or
/// when a cursor is minted. A handle never exists in an unvalidated state, and
/// nothing here can reach a caller mid-stream: once reads are moving, a failure
/// is a [`CursorError`], which names the file it came from.
///
/// Design: `doc/devel/ng/spec/alignment_file.md` §4,
/// `doc/devel/ng/arch/alignment_file.md` §2.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum AlignmentFileError {
    /// The file could not be opened, or its header could not be parsed.
    #[error("opening alignment file '{path}' failed")]
    Open {
        path: PathBuf,
        #[source]
        source: io::Error,
    },

    /// `@HD SO` is absent or is not `coordinate` (spec §3.1, check 1).
    ///
    /// Cheap and checked first, but **not** a substitute for the streaming
    /// order guard: this catches a file that does not *claim* to be sorted,
    /// while [`CursorError::OutOfOrderRead`] catches the file that claims it
    /// and lies.
    ///
    /// `sort_order` carries what the header actually said — `queryname` and a
    /// missing tag are the same rejection but not the same diagnosis, and a
    /// user whose file says `SO:queryname` should be told that rather than
    /// left to guess. `None` means the tag (or the whole `@HD` line) was
    /// absent. (Arch §2 sketched this variant with only `path`; signatures
    /// there are illustrative, and production's equivalent carries the observed
    /// value too.)
    #[error(
        "alignment file '{path}' is not coordinate-sorted (@HD SO is {})",
        match sort_order {
            Some(value) => format!("'{value}'"),
            None => "missing".to_string(),
        }
    )]
    NotCoordinateSorted {
        path: PathBuf,
        sort_order: Option<String>,
    },

    /// The file's `@SQ` list is not the reference's contig table (spec §3.1,
    /// check 2). `detail` comes from `ContigList::first_disagreement` and names
    /// the first differing field and index.
    ///
    /// This is the variant that closes the **permutation hole**: a file whose
    /// `@SQ` list is a re-ordering of the reference's passes a resolves-only
    /// check and then fetches the wrong contig for every read.
    #[error("alignment file '{path}' does not match the reference contig table: {detail}")]
    ContigReconcile { path: PathBuf, detail: String },

    /// An `@SQ` line carries an `M5` tag that is not a 32-character hex digest.
    ///
    /// A *missing* `M5` is ordinary input and never an error (spec §3.1) — the
    /// digest check is opportunistic, and refusing a file for not offering a
    /// bonus would punish the common case. A **malformed** one is different: it
    /// is a header its writer got wrong, and reading it as absent would hide a
    /// real defect. Worse, it would hide it *silently and completely* — a file
    /// whose digests are all malformed would report as untagged, and
    /// wrong-assembly detection would switch off with nobody told. Errors are
    /// not passed silently, so this is fatal at open.
    #[error("alignment file '{path}': @SQ '{contig}' has a malformed M5 tag ({detail})")]
    MalformedMd5 {
        path: PathBuf,
        contig: String,
        detail: String,
    },

    /// The index is missing or unparseable (spec §3.1, check 3).
    #[error("loading the index for alignment file '{path}' failed")]
    Index {
        path: PathBuf,
        #[source]
        source: AlignmentIndexError,
    },

    // The two `@RG SM` variants are gone. A file that names no sample and a
    // file that names several are both the read-group pre-pass's business now,
    // and it reports them better: `ReadGroupError` tells "no `@RG` at all" apart
    // from "an `@RG` with no `SM`" as separate variants with separate remedies,
    // where the pair here folded both into one variant carrying an `Option`. A
    // file naming *several* samples has stopped being a fault entirely.
    //
    // `OutOfOrderRead` and `Filter` are gone too, and they left with the reader
    // that raised them. Both were yielded *in the item stream* by the per-region
    // query's order guard, which Milestone F deleted; the rules they carried did
    // not go with it — a cursor raises `CursorError::OutOfOrderRead` for the
    // first and `CursorError::ReadRecord` for the second, each naming the file
    // it came from, which the per-region pair could only do because a fresh
    // guard was built for every query. Kept as unconstructible variants they
    // would have been two error cases nothing could ever produce.
    /// The requested region is invalid, or names a contig absent from the
    /// reference.
    #[error(
        "invalid region: contig {} [{}, {}]",
        region.contig.get(), region.start.get(), region.end.get()
    )]
    Region { region: GenomeRegion },

    /// A cursor was asked for over a chromosome this file does not have.
    ///
    /// Every layer below would answer this as "nothing here", which is indistinguishable from
    /// a chromosome that is genuinely empty — so it is refused where the caller can still act
    /// on it rather than reported as an absence of reads.
    #[error(
        "alignment file '{path}' has {contigs_in_file} contigs and no contig {}",
        contig.get()
    )]
    CursorContigNotInFile {
        path: PathBuf,
        contig: crate::ng::types::ContigId,
        contigs_in_file: usize,
    },

    /// A cursor was asked for over a format that cannot serve one yet.
    ///
    /// Its own variant rather than an `Open` failure carrying a sentence, because nothing
    /// failed to open: the file is fine and the feature is absent. An operator who reads
    /// "opening … failed" goes looking at the file.
    #[error(
        "reading '{path}' through a cursor is not implemented for {kind} yet; \
         it arrives with the CRAM arm"
    )]
    CursorFormatUnsupported { path: PathBuf, kind: &'static str },

    /// The reference could not answer for a contig this file declares, at the point a cursor
    /// was made for it.
    ///
    /// Surfaced here rather than at the first read for the reason the open gate exists: a
    /// mismatched reference is a fault in what the run was given, not in the region a caller
    /// happened to ask for first.
    #[error("the reference cannot serve alignment file '{path}'")]
    Reference {
        path: PathBuf,
        #[source]
        source: crate::ng::ref_seq::RefSeqError,
    },

    /// A CRAM was opened against a reference that carries no FASTA — a
    /// `.fai`-only `ReferenceInfo`.
    ///
    /// CRAM stores its sequences as differences from the reference, so decoding
    /// a record *requires the bases themselves*. A `.fai` describes a genome's
    /// geometry and holds no bases, so there is nothing to decode against. BAM
    /// is unaffected: it stores its own sequences.
    ///
    /// Raised at open rather than at the first query, so the run fails while
    /// the cause is still obvious (owner, 2026-07-20).
    #[error(
        "'{path}' is a CRAM, which can only be decoded against the reference \
         sequence itself, but the reference was read from a .fai index (which \
         holds no bases) — supply the reference FASTA"
    )]
    CramNeedsReferenceFasta { path: PathBuf },

    /// Querying the parsed index for a region's chunks failed.
    ///
    /// Distinct from [`Self::Region`]: the region has already been checked
    /// against the contig table by the time this can fire, so a failure here is
    /// the *index* being unusable — corrupt or truncated — not the caller
    /// asking for something silly. Reporting it as a bad region would name the
    /// wrong culprit and throw away the cause.
    #[error(
        "querying the index of '{path}' for contig {} [{}, {}] failed",
        region.contig.get(), region.start.get(), region.end.get()
    )]
    IndexQuery {
        path: PathBuf,
        region: GenomeRegion,
        #[source]
        source: io::Error,
    },
}

/// What [`check_assembly`] found: how many contigs could be compared, out of
/// how many there are.
///
/// Returned as **information, not a complaint**. `compared == 0` is ordinary —
/// it means no contig carried a digest on both sides — and a caller that
/// ignores this is behaving correctly. It exists so a caller that wants to say
/// "assembly verified for 18 of 24 contigs" can, rather than having to imply a
/// guarantee it does not have.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AssemblyCheck {
    /// Contigs where **both** sides carried a digest and the two agreed.
    pub compared: usize,
    /// Contigs in the reference — the denominator, so `compared` can be read
    /// as a fraction of the genome actually verified.
    pub total: usize,
}

/// A contig whose `@SQ M5` disagrees with the reference's digest: the file is
/// aligned to a **different assembly** with the same contig names and lengths.
///
/// A separate type rather than an [`AlignmentFileError`] variant, because
/// [`check_assembly`] runs long after the streams are gone — so this can never
/// appear as a stream item, and a variant of that enum would be unreachable
/// there (spec §4).
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error(
    "alignment file '{path}' is aligned to a different assembly: contig \
     '{contig}' has M5 {} but the reference's is {}",
    format_md5_hex(*observed),
    format_md5_hex(*expected)
)]
pub struct AssemblyMismatch {
    pub path: PathBuf,
    pub contig: String,
    pub expected: [u8; 16],
    pub observed: [u8; 16],
}

/// Compare the `@SQ M5` tags captured at open against a **verified**
/// `ReferenceInfo` — the one carrying real per-contig digests, available once
/// the caller has joined `reference_info`'s background verification.
///
/// **Why this is deferred rather than done at open.** With a `.fai` present,
/// `reference_info` hands back the cheap table immediately and verifies the
/// genome on a background thread, so the `ReferenceInfo` a file is opened
/// against usually carries no digests at all. Waiting for them would mean
/// blocking startup on a whole-genome read. Instead the tags are captured for
/// free at open and compared here, at the point the caller joins that
/// verification — before anything is published. A wrong-assembly run therefore
/// wastes its work before aborting, which is the right trade for something
/// catastrophic and rare (spec §3.1).
///
/// **A contig is compared only when both sides carry a digest**, and a contig
/// that carries none on either side is skipped, not failed. Many BAM/CRAMs have
/// no `M5` at all, and a `.fai`-only reference has none either; refusing or
/// nagging would punish the common case for missing an opportunistic extra
/// check. (A *malformed* tag is a different matter and was already rejected at
/// open, so everything reaching here is absent or valid.)
///
/// `observed` is indexed by `ContigId`, parallel to `verified.contigs` — which
/// the open gate guarantees, having proved the file's `@SQ` list *is* the
/// reference's contig table.
pub fn check_assembly(
    path: &Path,
    observed: &[Option<[u8; 16]>],
    verified: &ReferenceInfo,
) -> Result<AssemblyCheck, AssemblyMismatch> {
    // The open gate proves the file's `@SQ` list *is* the reference's contig
    // table, so these are parallel. Checked because this is `pub` and takes two
    // loose slices: if they ever diverged, pair `i` would compare the file's
    // contig against a *different* reference contig — missing a real mismatch,
    // or naming the wrong contig in the error.
    debug_assert_eq!(
        observed.len(),
        verified.contigs.len(),
        "observed digests must be parallel to the reference's contigs"
    );

    let mut compared = 0;

    for (observed_digest, contig) in observed.iter().zip(verified.contigs.iter()) {
        let (Some(observed_digest), Some(expected)) = (observed_digest, contig.md5) else {
            continue;
        };

        if *observed_digest != expected {
            return Err(AssemblyMismatch {
                path: path.to_path_buf(),
                contig: contig.name.clone(),
                expected,
                observed: *observed_digest,
            });
        }
        compared += 1;
    }

    Ok(AssemblyCheck {
        compared,
        // The reference's count, not the number of pairs walked: if the two
        // slices disagree in length — which the open gate rules out — the
        // shortfall shows up as `compared < total` rather than vanishing.
        total: verified.contigs.len(),
    })
}

/// One sample's files, opened and cross-validated. Built once, queried per
/// region.
///
/// **A sample is usually several files, sometimes one.** The dominant reason
/// for several is that the same sample was *sequenced in several experiments* —
/// separate runs, libraries, often separate projects — rather than a per-lane
/// or per-chromosome split. That shapes everything here: per-experiment files
/// span the same coordinate range, so the merge interleaves at essentially
/// every position and ties are routine; they carry different read groups and
/// chemistries, so "are these really all the same sample?" is a live question;
/// and "how did each file behave?" is worth reporting separately.
///
/// This layer does exactly three things — check what can only be checked
/// *across* files, merge k ordered streams into one, and present an entry point
/// whose one-file and k-file arms are indistinguishable. It contains **no
/// reader logic**: opening, validating, seeking and filtering are
/// [`AlignmentFile`](open_bam::AlignmentFile)'s.
///
/// Not `Clone` — cloning it would hand out two views onto one set of reader
/// pools and one set of tallies. The files themselves are shared by `Arc`, so
/// the `Vec` no longer refuses a `Clone` derive on its own: the rule is a
/// decision now, not a compile error.
#[derive(Debug)]
pub struct SampleReads {
    /// `Arc` because a reader **owns** the file it reads through, and outlives the
    /// borrow of this `SampleReads` that made it: a generator is lent `&SampleReads`
    /// per `next_locus` call while holding its reader across many of them (arch
    /// `locus_generation_pileup.md` §2.2). That is a region stream for the STR
    /// generator, which hands its pooled reader back on `Drop`, and since D1 a
    /// **chromosome-lived cursor** for the generic one, which has no pool to hand
    /// anything back to. One pool per file either way — the `Arc` is shared, not
    /// cloned deeply.
    files: Vec<Arc<AlignmentFile>>,
    /// The one name all k files agree on.
    sample_name: String,
}

/// Which sample a reader was opened for — see [`SampleReads::identity`].
///
/// **Two samples are the same one when they read the same files**, not when they carry the
/// same name. Names are the wrong test in both directions: a cohort can hold two `SampleReads`
/// built from different files under one name (the same individual sequenced twice, opened as
/// two samples), and one file set can be opened twice under names that differ only by how the
/// caller spelled them. The files are what the reads actually come from, so the files are what
/// this compares.
///
/// It holds shared handles rather than raw addresses, so a freed sample's address cannot be
/// reused by a later one and quietly compare equal.
#[derive(Clone)]
pub struct SampleIdentity {
    files: Vec<Arc<AlignmentFile>>,
}

impl PartialEq for SampleIdentity {
    fn eq(&self, other: &Self) -> bool {
        self.files.len() == other.files.len()
            && std::iter::zip(&self.files, &other.files).all(|(a, b)| Arc::ptr_eq(a, b))
    }
}

impl Eq for SampleIdentity {}

impl std::fmt::Debug for SampleIdentity {
    /// The paths, because that is what an operator can act on. A pointer would say only that
    /// two things differ, which the comparison already said.
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_list()
            .entries(
                self.files
                    .iter()
                    .map(|file| file.path().display().to_string()),
            )
            .finish()
    }
}

impl SampleReads {
    /// Open every file, validating each, then check they all name one sample.
    /// Both happen before any read flows.
    ///
    /// Every read a file yields carries its **read group**, which is
    /// **downstream-meaningful, not just provenance**: a read group is one
    /// library preparation, and that is the grouping a per-chemistry error model
    /// keys on, because PCR stutter and per-base error are properties of the
    /// preparation rather than of the individual. This is where that tag is
    /// settled — a positional file index, which is what reads used to carry,
    /// could not answer it, since a sample's files need not be one library each.
    ///
    /// A per-file failure is wrapped with the position of the file that raised
    /// it, so a message can always name the culprit.
    ///
    /// **`reference` is the run's, not this sample's**, and that is what makes
    /// a cohort affordable: a CRAM decodes against the reference bases, and
    /// [`OpenReference`] holds the one copy of them that every sample's every
    /// file shares (`reference.rs`). A per-sample reference would put the whole
    /// genome in memory once per sample.
    pub fn open(
        sample: &SampleReadGroups,
        read_groups: &ReadGroups,
        reference: &OpenReference,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
    ) -> Result<Self, IngestError> {
        if sample.read_groups.is_empty() {
            return Err(IngestError::NoFiles);
        }

        check_one_sample(sample, read_groups)?;

        // The files this sample's read groups live in, each once, in the order
        // its read groups name them — which is the order the run's input list
        // gave, because that is the order the identifiers were minted in.
        let mut paths: Vec<Arc<Path>> = Vec::new();
        for id in &sample.read_groups {
            let file = &read_groups.get(*id).file;
            if !paths.iter().any(|seen| seen == file) {
                paths.push(Arc::clone(file));
            }
        }

        let files = paths
            .iter()
            .enumerate()
            .map(|(source_file_index, path)| {
                let resolution = resolution_for(path, sample, read_groups)?;
                AlignmentFile::open(
                    path,
                    reference,
                    filter_config,
                    build_index_if_missing,
                    resolution,
                )
                .map_err(|source| IngestError::File {
                    source_file_index,
                    source,
                })
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            files,
            sample_name: sample.sample.to_string(),
        })
    }

    /// Open the **one** sample these files hold, building the run's read-group
    /// table on the way.
    ///
    /// The convenience for a single-sample tool: it takes paths, as everything
    /// did before read groups existed, and it fails the same way if they turn
    /// out to name more than one sample. Eight of this repo's tools want exactly
    /// this, so the assertion lives here once rather than eight times.
    ///
    /// A tool that handles a cohort must not use it: it calls
    /// [`build_read_groups`] itself and opens one `SampleReads` per entry of
    /// [`ReadGroups::read_groups_per_sample`], which is the whole point of the
    /// table being run-wide.
    pub fn open_only_sample(
        paths: &[PathBuf],
        reference: &OpenReference,
        filter_config: ReadFilterConfig,
        build_index_if_missing: bool,
    ) -> Result<Self, IngestError> {
        let read_groups = build_read_groups(paths).map_err(IngestError::ReadGroups)?;

        match read_groups.read_groups_per_sample() {
            [only] => Self::open(
                only,
                &read_groups,
                reference,
                filter_config,
                build_index_if_missing,
            ),
            [] => Err(IngestError::NoFiles),
            // Listed per **read group**, not per distinct sample: the two
            // vectors have to stay parallel so each file is paired with what it
            // claims, which is the only way to see which one is the stray. A
            // file whose read groups name two samples appears twice, correctly.
            _ => Err(IngestError::SampleNameMismatch {
                files: read_groups
                    .iter()
                    .map(|(_, group)| group.file.to_path_buf())
                    .collect(),
                names: read_groups
                    .iter()
                    .map(|(_, group)| group.sample.to_string())
                    .collect(),
            }),
        }
    }

    /// The one sample all this sample's files name.
    pub fn sample_name(&self) -> &str {
        &self.sample_name
    }

    /// Which sample this is, as a token a caller can keep and compare later.
    ///
    /// **What it is for.** A locus generator opens a cursor for one sample and keeps it for a
    /// whole chromosome, but it is handed a `&SampleReads` afresh on every call. Without
    /// something to compare, a generator handed a *second* sample would keep answering out of
    /// the first sample's files — every sample in a cohort reporting the first one's reads,
    /// with no error and output of exactly the right shape. Holding this token lets the
    /// generator notice and refuse.
    ///
    /// Cheap: one atomic increment per file, and a sample has one file in almost every run.
    pub fn identity(&self) -> SampleIdentity {
        SampleIdentity {
            files: self.files.clone(),
        }
    }

    /// How many files this sample has.
    ///
    /// **Not** the space a read's `read_group` lives in: that identifier is
    /// run-wide, so it indexes the whole run's read-group table, not this
    /// sample's files. The two coincide only when a run holds one single-file
    /// sample, which is true of most fixtures and of nothing else.
    pub fn file_count(&self) -> usize {
        self.files.len()
    }

    /// One cursor for this sample over one chromosome — **the layer a caller holds, and
    /// since Milestone F the only way to read a file here.**
    ///
    /// Made **once per chromosome**, it answers every region on that chromosome from readers
    /// that stay where they are. What it replaced returned a fresh stream per region and
    /// reopened the question every time — index query, seek, decode, filter — for regions
    /// that usually sat in the block the reader was already in.
    ///
    /// **One reference accessor per file, taken here once** for the cursor's whole life
    /// rather than rebuilt per region (perf review L2). `RawRefSeq` impls are stateful
    /// readers, so the k files cannot share one: `WindowedRefSeq` holds an open per-contig
    /// reader and is `Send` but not `Sync`, precisely because each consumer is meant to own
    /// one. A factory gives each file cursor its own, which is what the type requires; the
    /// caller writes one closure.
    ///
    /// A read's step-1 tally is the cursor's from here on
    /// ([`SampleCursor::read_group_counts`]) — it needs no hand-back, because the filter now
    /// lives as long as the cursor rather than as long as one region.
    pub fn cursor<R, F>(
        &self,
        contig: ContigId,
        mut make_reference: F,
    ) -> Result<SampleCursor<R>, IngestError>
    where
        R: RawRefSeq,
        F: FnMut() -> R,
    {
        let mut cursors = Vec::with_capacity(self.files.len());
        for (source_file_index, file) in self.files.iter().enumerate() {
            cursors.push(file.cursor(contig, make_reference()).map_err(|source| {
                IngestError::File {
                    source_file_index,
                    source,
                }
            })?);
        }
        Ok(SampleCursor::new(cursors))
    }

    /// What the deferred assembly check needs, per file: the path, and the
    /// `@SQ M5` tags captured at its open.
    ///
    /// Forwarded rather than merged, because a disagreement has to name *which
    /// file* is aligned to the wrong assembly. The caller feeds each pair to
    /// [`check_assembly`] once it has joined `reference_info`'s background
    /// verification.
    ///
    /// Shaped as pairs rather than as spec §3.4's bare `sq_md5s`, because
    /// `check_assembly` needs both halves and this is its only consumer —
    /// and because handing out the `AlignmentFile`s themselves would let a
    /// caller query one file directly, skipping the merge and the
    /// same-file-twice check that are this type's whole reason to exist.
    pub fn assembly_inputs(&self) -> impl Iterator<Item = (&Path, &[Option<[u8; 16]>])> {
        self.files.iter().map(|file| (file.path(), file.sq_md5s()))
    }
}

/// Enforce the invariant that gives this type its meaning: **one open serves
/// exactly one sample.**
///
/// A file may declare read groups for several samples, so the check is over the
/// read groups this open was handed, not over the files they live in. A caller
/// that took its argument from `ReadGroups::read_groups_per_sample` can never
/// trip it — the grouping is by sample. It earns its keep for the caller that
/// assembles an open by hand, which is how tools and tests do it, and it is what
/// stops a foreign read group being read as part of a sample (spec §4).
fn check_one_sample(
    sample: &SampleReadGroups,
    read_groups: &ReadGroups,
) -> Result<(), IngestError> {
    let foreign: Vec<&ReadGroup> = sample
        .read_groups
        .iter()
        .map(|id| read_groups.get(*id))
        .filter(|group| group.sample != sample.sample)
        .collect();

    if foreign.is_empty() {
        return Ok(());
    }

    // Every read group is listed, not just the offenders, and the two vectors
    // are built from the same iteration so they stay parallel: the user has to
    // see which file claims which name to know which one is the stray.
    Err(IngestError::SampleNameMismatch {
        files: sample
            .read_groups
            .iter()
            .map(|id| read_groups.get(*id).file.to_path_buf())
            .collect(),
        names: sample
            .read_groups
            .iter()
            .map(|id| read_groups.get(*id).sample.to_string())
            .collect(),
    })
}

/// How this open reads one file's records.
///
/// A file declaring one read group resolves to [`Sole`](ReadGroupResolution::Sole)
/// and its records' `RG` tags are never read — the common case, and the one that
/// lets a file re-headered without rewriting its records be read at all. A file
/// declaring several resolves per record.
///
/// A file whose read groups name **several samples** resolves too: this open's
/// read groups become [`Mine`](RecordOwner::Mine) and the rest
/// [`OtherSample`](RecordOwner::OtherSample), and the record sources leave the
/// foreign records out — counted apart from every drop, because a read belonging
/// to someone else says nothing about how this read group behaved.
fn resolution_for(
    path: &Path,
    sample: &SampleReadGroups,
    read_groups: &ReadGroups,
) -> Result<ReadGroupResolution, IngestError> {
    let declared: Vec<(ReadGroupId, &ReadGroup)> = read_groups
        .iter()
        .filter(|(_, group)| *group.file == *path)
        .collect();

    Ok(match declared.as_slice() {
        // No read groups at all cannot happen — a file with no `@RG` fails when
        // the table is built, long before an open — but saying so beats letting
        // it fall into the arm below, which would build an empty `PerRecord`
        // that fails every record at read time instead of this file at open.
        [] => {
            return Err(IngestError::NoReadGroupsForFile {
                path: path.to_path_buf(),
            });
        }
        // One declared read group and it is this sample's: the records' own tags
        // need never be read. A lone *foreign* read group cannot reach here —
        // the file would hold no read group of this sample, so this open would
        // not have been given it.
        [(id, group)] if group.sample == sample.sample => ReadGroupResolution::Sole(*id),
        several => ReadGroupResolution::PerRecord(
            several
                .iter()
                .map(|(id, group)| {
                    let owner = if group.sample == sample.sample {
                        RecordOwner::Mine(*id)
                    } else {
                        RecordOwner::OtherSample
                    };
                    (group.id.clone(), owner)
                })
                .collect(),
        ),
    })
}

/// Sample-level failures — the layer above [`AlignmentFileError`].
///
/// The two own variants are the checks that can only be made *across* files.
/// Everything else arrives from one file and is **wrapped** with the index of
/// the file that raised it, rather than flattened: flattening would duplicate
/// every [`AlignmentFileError`] variant and lose which file was at fault
/// (`doc/devel/ng/arch/sample_reads.md` §2).
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum IngestError {
    /// The k files do not name the same sample. Fires at `SampleReads::open`,
    /// before any read (spec `sample_reads.md` §3.1).
    ///
    /// This guard earns its keep because a sample's files are usually separate
    /// *experiments* gathered by hand from different projects, where picking up
    /// a foreign file is the realistic failure mode.
    #[error(
        "the files of one sample name different samples: {}",
        files.iter().zip(names.iter())
            .map(|(file, name)| format!("'{}' names '{name}'", file.display()))
            .collect::<Vec<_>>().join(", ")
    )]
    SampleNameMismatch {
        files: Vec<PathBuf>,
        names: Vec<String>,
    },

    /// The run's read groups could not be read at all — a file with no `@RG`, an
    /// `@RG` with no `SM`, a library-name collision, or a header that would not
    /// parse. Raised before any file is opened.
    #[error("reading the input files' read groups failed")]
    ReadGroups(#[source] ReadGroupError),

    /// The run's read-group table holds no read group for a file this open was
    /// asked to read.
    ///
    /// Unreachable through the ordinary path: a file with no `@RG` fails when the
    /// table is built. It exists so that a caller assembling an open by hand from
    /// a table that never saw the file gets told, rather than getting a resolution
    /// that rejects every record it later reads.
    #[error("no read group in the run's table belongs to alignment file '{}'", path.display())]
    NoReadGroupsForFile { path: PathBuf },

    /// `SampleReads::open` was given no files.
    ///
    /// A sample with no files has no reads and no name, so there is nothing to
    /// build — and reporting it as a name *mismatch*, which is where zero files
    /// would otherwise land, would blame the wrong thing entirely.
    #[error("a sample must have at least one alignment file")]
    NoFiles,

    // The top-level `DuplicateReadAcrossFiles` is gone, and Milestone F is where its doc
    // said it would go. It was raised by the per-region merge, which no longer exists; the
    // live version of the rule is
    // [`CursorError::DuplicateReadAcrossFiles`](CursorError::DuplicateReadAcrossFiles),
    // reached through [`Cursor`](Self::Cursor). Two variants for one condition was the
    // transitional state, not the destination.
    /// Anything one file raised, at open or mid-stream, tagged with which one.
    #[error("alignment file #{source_file_index} failed")]
    File {
        source_file_index: usize,
        #[source]
        source: AlignmentFileError,
    },

    /// A cursor could not be pointed at a region, or failed reading one.
    ///
    /// **Not tagged with a file index, unlike [`File`](Self::File)**, because a
    /// [`CursorError`] carries its own provenance: all but one variant names the path it
    /// came from, which is what a run holding one cursor per chromosome per generator per
    /// worker actually needs — a bare index would not say which of them spoke.
    ///
    /// **The exception is
    /// [`DuplicateReadAcrossFiles`](CursorError::DuplicateReadAcrossFiles)**, which names two
    /// cursor *indexes* and no path. It used to have a twin at this level, raised by the
    /// per-region merge; Milestone F deleted that path, so there is now exactly one variant
    /// for the condition and a caller matching on it cannot miss half the cases.
    #[error("this sample's cursor failed")]
    Cursor(#[from] CursorError),
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::ng::read::aligned_read::AlignedRead;
    use crate::ng::read::filtering::ReadFilterCounts;
    use crate::ng::types::{ContigId, Position};

    // -----------------------------------------------------------------
    // SampleReads::open (E1) — T12b
    // -----------------------------------------------------------------

    use tempfile::TempDir;

    use crate::bam::index_preflight::preflight_alignment_indexes;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, header, indexed_bam, indexed_named_bam, matching_contigs, named_bam,
        one_read, only_tally, read_named_with_length,
    };

    /// An indexed BAM whose one read group names `sample`.
    /// A one-read BAM whose read group names `sample`, under `file_name`.
    ///
    /// **Every file a test opens together must have a different name.** These
    /// fixtures declare no `LB`, so each read group's library name is
    /// synthesized from its sample, its `@RG ID` and its file's name — and two
    /// files sharing all three would be indistinguishable, which the pre-pass
    /// rejects. Real inputs opened together have different names for the same
    /// reason, so this is the fixture matching reality, not working around it.
    fn bam_naming(sample: &str, file_name: &str) -> (TempDir, PathBuf) {
        bam_naming_with_reads(sample, 1, file_name)
    }

    fn bam_naming_with_reads(sample: &str, reads: usize, file_name: &str) -> (TempDir, PathBuf) {
        // Read names carry the file's, because two files of one sample opened
        // together must not repeat a name at a position: that is the
        // same-file-twice check, and it fires on identical `(qname, flag,
        // position)` across files whether or not the files really are the same.
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records: Vec<_> = (0..reads)
            .map(|i| read_named_with_length(&format!("{stem}-r{i}"), 0, 1 + i, 30))
            .collect();
        let header = header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg1", Some(sample))],
        );
        let (dir, path) = named_bam(&header, &records, file_name);
        preflight_alignment_indexes(std::slice::from_ref(&path), true).expect("build index");
        (dir, path)
    }

    /// **T12b — two files naming different samples, rejected at `open`, before
    /// any read.**
    ///
    /// The single-file half is `alignment_file`'s T12a. This is the check that
    /// cannot live one layer down, because agreement is not a property of any
    /// one file — and it is the one that catches the realistic mistake:
    /// gathering a sample's experiments by hand from different projects and
    /// picking up a stray.
    #[test]
    fn t12b_files_naming_different_samples_are_rejected_at_open() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_first_dir, first) = bam_naming("NA12878", "first.bam");
        let (_second_dir, second) = bam_naming("NA12892", "second.bam");

        let error = SampleReads::open_only_sample(
            &[first.clone(), second.clone()],
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect_err("two samples is not one sample");

        match &error {
            IngestError::SampleNameMismatch { files, names } => {
                assert_eq!(files, &[first, second]);
                assert_eq!(names, &["NA12878".to_string(), "NA12892".to_string()]);
            }
            other => panic!("expected SampleNameMismatch, got {other:?}"),
        }

        // The message has to pair each path with what it claims, or the user
        // cannot tell which file is the stray.
        let message = error.to_string();
        assert!(message.contains("NA12878"), "{message}");
        assert!(message.contains("NA12892"), "{message}");
    }

    #[test]
    fn files_agreeing_on_one_sample_open_and_expose_it() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_first_dir, first) = bam_naming("NA12878", "first.bam");
        let (_second_dir, second) = bam_naming("NA12878", "second.bam");

        let sample = SampleReads::open_only_sample(
            &[first, second],
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("both name the same sample");

        assert_eq!(sample.sample_name(), "NA12878");
        assert_eq!(sample.file_count(), 2);
    }

    /// A per-file failure is wrapped with **which** file failed — with k files
    /// gathered from k projects, "one of them is not coordinate-sorted" is not
    /// a usable message.
    #[test]
    fn a_per_file_failure_names_the_file_that_raised_it() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_good_dir, good) = bam_naming("NA12878", "good.bam");
        let (_bad_dir, bad) = indexed_bam(
            &header(
                Some("queryname"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            &[read_named_with_length("r", 0, 1, 30)],
        );

        let error = SampleReads::open_only_sample(
            &[good, bad],
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect_err("the second file is not coordinate-sorted");

        match error {
            IngestError::File {
                source_file_index,
                source: AlignmentFileError::NotCoordinateSorted { .. },
            } => assert_eq!(source_file_index, 1, "the second file, not the first"),
            other => panic!("expected File{{1, NotCoordinateSorted}}, got {other:?}"),
        }
    }

    /// Counts stay **per read group**, never summed: a bad run shows up as one library
    /// with an anomalous drop rate, and summing erases exactly that signal.
    ///
    /// Read off the sample's cursor since Milestone F. The rule is unchanged; what moved is
    /// where the number lives — the per-region streams folded theirs back into the file as
    /// each one ended, and a filter that lasts as long as a cursor simply has it.
    ///
    /// The files are given *different* read counts on purpose. Asserting only
    /// that there are k tallies would pass against an implementation returning
    /// the summed total k times — which is precisely the shape spec §3.3 argues
    /// against.
    #[test]
    fn counts_are_reported_per_read_group_and_not_summed() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_first_dir, first) = bam_naming_with_reads("NA12878", 3, "first.bam");
        let (_second_dir, second) = bam_naming_with_reads("NA12878", 1, "second.bam");

        let sample = SampleReads::open_only_sample(
            &[first, second],
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("opens");

        // A read group appears only once this sample's reads have met it, so before any read
        // there is no *keyed* tally — not two zeroed ones. What is there is the `None`-keyed
        // rider `ReadFilter::counts` always carries for records skipped as another sample's,
        // which belongs to no read group and is zero here.
        let mut cursor = sample
            .cursor(ContigId(0), reference_bases)
            .expect("a cursor for contig 0");
        assert!(
            cursor
                .read_group_counts()
                .iter()
                .all(|(id, tally)| id.is_none() && *tally == ReadFilterCounts::default()),
            "no read group has been met yet: {:?}",
            cursor.read_group_counts(),
        );

        cursor
            .move_to_region(whole_first_contig())
            .expect("on this chromosome");
        while let Some(read) = cursor.next_read() {
            read.expect("resolves");
        }

        // Each fixture file declares its own `@RG`, so the run's table minted
        // two identifiers and the reads split between them — 3 and 1, never
        // summed to 4. Summing is what erases one bad run among several.
        let counts = cursor.read_group_counts();
        assert_eq!(counts.len(), 2, "one tally per read group met");
        let kept: Vec<u64> = counts.iter().map(|(_, tally)| tally.kept).collect();
        assert_eq!(kept, vec![3, 1]);
        assert_eq!(
            counts
                .iter()
                .map(|(id, _)| id.expect("stamped").get())
                .collect::<Vec<_>>(),
            vec![0, 1],
            "keyed by the run-wide identifier, not by the file's position"
        );

        assert_eq!(
            sample.assembly_inputs().count(),
            2,
            "and one (path, digests) pair per file"
        );
    }

    /// **Three files, the stray in the middle.** A two-file test cannot tell a
    /// correctly-ordered `names` vector from a reversed one; this can, and it
    /// is what pins `files` and `names` being parallel.
    #[test]
    fn the_mismatch_lists_every_file_in_order_so_the_stray_is_identifiable() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_naming("NA12878", "a.bam");
        let (_stray_dir, stray) = bam_naming("NA12892", "stray.bam");
        let (_b_dir, b) = bam_naming("NA12878", "b.bam");

        let error = SampleReads::open_only_sample(
            &[a.clone(), stray.clone(), b.clone()],
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect_err("one stray is enough");

        match &error {
            IngestError::SampleNameMismatch { files, names } => {
                assert_eq!(files, &[a, stray, b]);
                assert_eq!(
                    names,
                    &[
                        "NA12878".to_string(),
                        "NA12892".to_string(),
                        "NA12878".to_string()
                    ],
                    "names must be parallel to files, not the distinct set"
                );
            }
            other => panic!("expected SampleNameMismatch, got {other:?}"),
        }
    }

    /// A one-file sample is ordinary — the arm `SampleReads` exists to make
    /// indistinguishable from the k-file one.
    #[test]
    fn a_single_file_sample_opens() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_dir, only) = bam_naming("NA12878", "only.bam");

        let sample =
            SampleReads::open_only_sample(&[only], &reference, ReadFilterConfig::default(), false)
                .expect("one file is a sample");
        assert_eq!(sample.file_count(), 1);
        assert_eq!(sample.sample_name(), "NA12878");
    }

    /// Zero files is not a sample. Reported as such, rather than falling
    /// through to a name *mismatch* with two empty lists — which would blame
    /// the wrong thing and render as a truncated sentence.
    #[test]
    fn a_sample_with_no_files_is_rejected_as_such() {
        let (_reference_dir, reference) = fixture_reference(false);

        let error =
            SampleReads::open_only_sample(&[], &reference, ReadFilterConfig::default(), false)
                .expect_err("no files is not a sample");
        assert!(matches!(error, IngestError::NoFiles));
        assert!(error.to_string().contains("at least one"), "{error}");
    }

    // -----------------------------------------------------------------
    // SampleCursor — T11 and the resumable shape
    // -----------------------------------------------------------------

    use crate::ng::read::input::test_fixtures::{FIXTURE_CONTIGS, bam_header};
    use crate::ng::ref_seq::InMemoryRefSeq;

    fn reference_bases() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(
            FIXTURE_CONTIGS
                .iter()
                .map(|(_, length)| vec![b'A'; *length])
                .collect(),
        )
    }

    fn whole_first_contig() -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(FIXTURE_CONTIGS[0].1 as u64),
        }
    }

    fn sample_over(paths: &[PathBuf]) -> (TempDir, SampleReads) {
        let (reference_dir, reference) = fixture_reference(false);
        let sample =
            SampleReads::open_only_sample(paths, &reference, ReadFilterConfig::default(), false)
                .expect("opens");
        (reference_dir, sample)
    }

    /// `(qname, read group)` per read of the region a cursor is pointed at — shared by
    /// `drain` and by the detached-cursor tests, so "the same reads" means literally the same
    /// code path on both sides of their comparison.
    fn collect_reads(cursor: &mut SampleCursor<InMemoryRefSeq>) -> Vec<(String, usize)> {
        let mut reads = Vec::new();
        while let Some(item) = cursor.next_read() {
            let read = item.expect("no fatal error");
            reads.push((
                String::from_utf8_lossy(&read.qname).into_owned(),
                read.read_group.get() as usize,
            ));
        }
        reads
    }

    /// A cursor over the whole first contig, already pointed at it.
    fn cursor_over_first_contig(sample: &SampleReads) -> SampleCursor<InMemoryRefSeq> {
        let mut cursor = sample
            .cursor(ContigId(0), reference_bases)
            .expect("a cursor for contig 0");
        cursor
            .move_to_region(whole_first_contig())
            .expect("on this chromosome");
        cursor
    }

    fn drain(sample: &SampleReads) -> Vec<(String, usize)> {
        collect_reads(&mut cursor_over_first_contig(sample))
    }

    /// **T11 — the caller cannot tell the arms apart.**
    ///
    /// A one-file sample yields exactly what the same file yields as one of two inputs whose
    /// other file is empty over the region — asserted through the same [`SampleCursor`] type,
    /// so the merge-free arm is provably equivalent rather than merely believed to be.
    #[test]
    fn t11_the_single_file_arm_matches_the_merged_arm() {
        let reads = [
            read_named_with_length("r1", 0, 1, 30),
            read_named_with_length("r2", 0, 30, 30),
            read_named_with_length("r3", 0, 60, 30),
        ];
        let (_only_dir, only) =
            indexed_named_bam(&bam_header(&matching_contigs()), &reads, "only.bam");
        // Empty over contig 0: its only read is on the *second* contig.
        let (_empty_dir, empty) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &[read_named_with_length("elsewhere", 1, 1, 30)],
            "empty.bam",
        );

        let (_a_dir, single) = sample_over(std::slice::from_ref(&only));
        let (_b_dir, merged) = sample_over(&[only, empty]);

        assert!(
            matches!(
                single.cursor(ContigId(0), reference_bases),
                Ok(SampleCursor::Single(_))
            ),
            "one file must take the merge-free arm"
        );
        assert!(
            matches!(
                merged.cursor(ContigId(0), reference_bases),
                Ok(SampleCursor::Merged(_))
            ),
            "two files must take the merge"
        );

        assert_eq!(
            drain(&single),
            drain(&merged),
            "the arms are indistinguishable to the caller"
        );
        assert_eq!(drain(&single).len(), 3);
    }

    /// **The cursor outlives the borrow that made it, and the sample too.**
    ///
    /// `cursor` takes `&self`, so the borrow ends when the call returns; the cursor it hands
    /// back carries no lifetime, holding its own open descriptor per file. Dropping the
    /// `SampleReads` before a single read has been pulled is the sharpest statement of that:
    /// **it would not compile at all** if the cursor still borrowed (`cannot move out of
    /// sample because it is borrowed`).
    ///
    /// **What this does not do is fail at run time**, and saying so is the honest version:
    /// once it compiles, safe Rust makes `after == expected` unfalsifiable. It is a
    /// compile-time anchor, and it is the shape both generators actually need — a resumable
    /// generator holds a cursor across many `next_locus` calls while being lent
    /// `&SampleReads` per call.
    ///
    /// The reads are compared against a control drain rather than merely counted, because
    /// "the cursor survives" and "the cursor survives and still yields the right reads" are
    /// different claims and only the second one is worth anything.
    #[test]
    fn a_cursor_outlives_the_sample_reads_it_was_made_from() {
        let reads = [
            read_named_with_length("r1", 0, 1, 30),
            read_named_with_length("r2", 0, 30, 30),
            read_named_with_length("r3", 0, 60, 30),
        ];
        let (_bam_dir, path) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &reads,
            "outlives_the_borrow.bam",
        );

        let (_reference_dir, sample) = sample_over(std::slice::from_ref(&path));
        let expected = drain(&sample);
        assert_eq!(expected.len(), 3, "the fixture must actually carry reads");

        let mut cursor = cursor_over_first_contig(&sample);

        // The load-bearing line: the sample — and with it the `Vec` of files —
        // is gone before the first read is pulled.
        drop(sample);

        assert_eq!(
            collect_reads(&mut cursor),
            expected,
            "a cursor that outlived its sample must still yield the same reads"
        );
    }

    /// **The k-file arm of the same property**, which is the production shape: a sample is
    /// usually several experiment files, so `Merged` — not `Single` — is what a generator
    /// will normally be holding.
    ///
    /// Worth its own test because the merged arm holds k cursors, hence k open descriptors
    /// and k tallies, past the sample's death; the interleave has to survive that, not merely
    /// compile.
    #[test]
    fn a_merged_cursor_outlives_the_sample_reads_it_was_made_from() {
        let (_first_dir, first) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &[
                read_named_with_length("a1", 0, 1, 30),
                read_named_with_length("a2", 0, 40, 30),
            ],
            "merged_outlives_a.bam",
        );
        let (_second_dir, second) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &[
                read_named_with_length("b1", 0, 20, 30),
                read_named_with_length("b2", 0, 60, 30),
            ],
            "merged_outlives_b.bam",
        );

        let (_reference_dir, sample) = sample_over(&[first, second]);
        let expected = drain(&sample);
        assert_eq!(
            expected.len(),
            4,
            "the fixture must carry both files' reads"
        );

        let mut cursor = cursor_over_first_contig(&sample);
        assert!(
            matches!(cursor, SampleCursor::Merged(_)),
            "two files must take the merge, or this tests the wrong arm"
        );

        drop(sample);

        assert_eq!(
            collect_reads(&mut cursor),
            expected,
            "a merged cursor that outlived its sample must still interleave the same reads"
        );
    }

    /// The same property in the shape that needs it: a **lifetime-free struct**
    /// holding one region's stream across several calls.
    ///
    /// That is `LocusGenerator::next_locus(&mut self, segment, reads:
    /// &SampleReads)` in miniature — the contract that cannot carry a lifetime
    /// (arch `locus_generation_pileup.md` §2.2) — and a locus generator yields many
    /// loci per segment, so it must hold its reader between calls. **This is now literally
    /// what both generators do**: each holds a [`SampleCursor`] for a chromosome and points
    /// it at each region. A **compile-time anchor**: with a borrowed cursor the `Generator`
    /// struct below needs a lifetime parameter and the test does not build. Named for what
    /// the body proves, since the run-time assertion repeats `drain`; the property that a
    /// *second* cursor cannot disturb a held one is next door, where it is falsifiable.
    #[test]
    fn a_cursor_can_be_stored_in_a_struct_without_a_lifetime() {
        struct Generator {
            cursor: Option<SampleCursor<InMemoryRefSeq>>,
        }

        impl Generator {
            /// `begin_segment`: mints the cursor, points it at the region, keeps it.
            fn begin(&mut self, reads: &SampleReads, region: GenomeRegion) {
                let mut cursor = reads
                    .cursor(region.contig, reference_bases)
                    .expect("a cursor for this contig");
                cursor.move_to_region(region).expect("on this chromosome");
                self.cursor = Some(cursor);
            }

            /// `next_locus`: lent the sample again, but pulls from the cursor
            /// it is already holding.
            ///
            /// `_reads` is unused **on purpose** — it stands in for
            /// `next_locus`'s per-call lend, and it is the parameter whose
            /// presence, not whose use, is the point. Do not "clean it up".
            fn next_qname(&mut self, _reads: &SampleReads) -> Option<String> {
                let read = self.cursor.as_mut()?.next_read()?.expect("no fatal error");
                Some(String::from_utf8_lossy(&read.qname).into_owned())
            }
        }

        let reads = [
            read_named_with_length("r1", 0, 1, 30),
            read_named_with_length("r2", 0, 30, 30),
            read_named_with_length("r3", 0, 60, 30),
        ];
        let (_bam_dir, path) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &reads,
            "held_across_calls.bam",
        );
        let (_reference_dir, sample) = sample_over(std::slice::from_ref(&path));

        let mut generator = Generator { cursor: None };
        generator.begin(&sample, whole_first_contig());

        let mut qnames = Vec::new();
        while let Some(qname) = generator.next_qname(&sample) {
            qnames.push(qname);
        }

        assert_eq!(
            qnames,
            vec!["r1".to_string(), "r2".to_string(), "r3".to_string()],
            "a held cursor must resume where the previous call left it"
        );
    }

    /// **A held cursor and a second one over the same sample are independent** — the
    /// run-time half of the resumable shape, and unlike the compile-time anchors above this
    /// one can actually fail.
    ///
    /// It falsifies any sharing of a file position, a record buffer or a reference reader
    /// between two live cursors on one file: each opens its own descriptor and owns its own
    /// accessor, which is exactly what the parallel fan-out will rely on. The threaded case
    /// is covered at the file level
    /// (`open_bam::tests::cursors_on_one_file_read_the_same_thing_from_many_threads`); this
    /// is the single-threaded interleaving.
    #[test]
    fn a_second_cursor_does_not_disturb_one_already_held() {
        let reads = [
            read_named_with_length("r1", 0, 1, 30),
            read_named_with_length("r2", 0, 30, 30),
            read_named_with_length("r3", 0, 60, 30),
        ];
        let (_bam_dir, path) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &reads,
            "two_live_streams.bam",
        );
        let (_reference_dir, sample) = sample_over(std::slice::from_ref(&path));

        let name = |item: Option<Result<AlignedRead, CursorError>>| {
            item.map(|read| {
                String::from_utf8_lossy(&read.expect("no fatal error").qname).into_owned()
            })
        };

        let mut held = cursor_over_first_contig(&sample);
        assert_eq!(name(held.next_read()).as_deref(), Some("r1"));

        // A whole second cursor, minted and drained while the first is parked.
        let second = collect_reads(&mut cursor_over_first_contig(&sample));
        assert_eq!(
            second
                .iter()
                .map(|(qname, _)| qname.as_str())
                .collect::<Vec<_>>(),
            vec!["r1", "r2", "r3"],
            "the second cursor is a full walk of its own"
        );

        assert_eq!(
            name(held.next_read()).as_deref(),
            Some("r2"),
            "the held cursor resumed where it was, undisturbed"
        );
        assert_eq!(name(held.next_read()).as_deref(), Some("r3"));
        assert!(held.next_read().is_none());
    }

    /// The cursor a generator stores has to be able to cross a thread boundary — **the
    /// merged arm included**, since that is the k-file production shape and the one a worker
    /// pool would move.
    ///
    /// `open_bam`'s sibling assertion covers the per-file `AlignmentCursor`; these are the
    /// two sample-level types, and an `Rc` or a non-`Send` field creeping into the merge
    /// would otherwise surface at the first parallel call site, a plan or two away.
    #[test]
    fn a_sample_cursor_is_send_in_both_arms() {
        fn assert_send<T: Send>() {}
        assert_send::<SampleCursor<InMemoryRefSeq>>();
        assert_send::<sample_cursor::MergedCursors<InMemoryRefSeq>>();
    }

    /// **The capability this step adds.** A file declaring several read groups
    /// is read, and each record goes to the read group its own `RG` tag names —
    /// which is the first point in the pipeline where a read can be attributed
    /// to the wrong library, and doing so would be silent.
    #[test]
    fn a_file_with_several_read_groups_resolves_each_record_by_its_tag() {
        use crate::ng::read::input::test_fixtures::{
            FixtureReadGroup, header_with_read_groups, indexed_named_bam,
            read_named_with_length_in_read_group,
        };

        let (_reference_dir, reference) = fixture_reference(false);
        let header = header_with_read_groups(
            Some("coordinate"),
            &matching_contigs(),
            &[
                FixtureReadGroup::new("rg1", Some("NA12878")).with_library("lib-A"),
                FixtureReadGroup::new("rg2", Some("NA12878")).with_library("lib-B"),
            ],
        );
        let (_dir, path) = indexed_named_bam(
            &header,
            &[
                read_named_with_length_in_read_group("a", 0, 1, 30, "rg1"),
                read_named_with_length_in_read_group("b", 0, 20, 30, "rg2"),
                read_named_with_length_in_read_group("c", 0, 40, 30, "rg1"),
            ],
            "two-groups.bam",
        );

        let read_groups = build_read_groups(std::slice::from_ref(&path)).expect("one sample");
        let [sample] = read_groups.read_groups_per_sample() else {
            panic!("both read groups name NA12878, so there is one sample");
        };
        let reads = SampleReads::open(
            sample,
            &read_groups,
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("a file with two read groups of one sample opens");

        let observed: Vec<(String, u32)> = collect_reads(&mut cursor_over_first_contig(&reads))
            .into_iter()
            .map(|(qname, read_group)| (qname, read_group as u32))
            .collect();

        // rg1 was declared first, so the table minted it 0 and rg2 got 1.
        assert_eq!(
            observed,
            vec![
                ("a".to_string(), 0),
                ("b".to_string(), 1),
                ("c".to_string(), 0),
            ],
            "each read carries the group its own tag names, not the file's first"
        );
    }

    /// An untagged record in such a file cannot be assigned to either candidate.
    /// Guessing would put a read in the wrong library, so it ends the run.
    #[test]
    fn an_untagged_record_in_a_multi_group_file_is_fatal() {
        use crate::ng::read::input::test_fixtures::{
            FixtureReadGroup, header_with_read_groups, indexed_named_bam,
            read_named_with_length_in_read_group,
        };

        let (_reference_dir, reference) = fixture_reference(false);
        let header = header_with_read_groups(
            Some("coordinate"),
            &matching_contigs(),
            &[
                FixtureReadGroup::new("rg1", Some("NA12878")),
                FixtureReadGroup::new("rg2", Some("NA12878")),
            ],
        );
        let (_dir, path) = indexed_named_bam(
            &header,
            &[
                read_named_with_length_in_read_group("tagged", 0, 1, 30, "rg1"),
                read_named_with_length("untagged", 0, 20, 30),
            ],
            "mixed.bam",
        );

        let read_groups = build_read_groups(std::slice::from_ref(&path)).expect("one sample");
        let [sample] = read_groups.read_groups_per_sample() else {
            panic!("one sample");
        };
        let reads = SampleReads::open(
            sample,
            &read_groups,
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("opens");

        let mut cursor = cursor_over_first_contig(&reads);
        let mut outcome = Vec::new();
        while let Some(item) = cursor.next_read() {
            outcome.push(item);
        }

        let error = outcome
            .iter()
            .find_map(|read| read.as_ref().err())
            .expect("the untagged record is fatal");
        // The name is in the innermost cause, so walk the chain rather than
        // rendering only the top, which says just "read filtering failed".
        let mut message = error.to_string();
        let mut cause: Option<&dyn std::error::Error> = std::error::Error::source(error);
        while let Some(source) = cause {
            message.push_str(&format!(": {source}"));
            cause = std::error::Error::source(source);
        }
        assert!(
            message.contains("untagged"),
            "the message names the read: {message}"
        );
    }

    // **`a_cram_with_several_read_groups_resolves_each_record_by_its_tag` died with the
    // per-region query, and it is subsumed rather than merely gone.** It read the same
    // three-record, two-read-group CRAM through the per-region source; the test immediately
    // below reads exactly that fixture through a cursor and asserts exactly the same
    // `(qname, read group)` triple. Both arms decide the read group inside
    // `decode_container_at`, which is now the only place a CRAM record's owner is settled, so
    // there was one rule under two names.

    /// **The same, through a cursor** — because a sample's libraries are the case a cursor
    /// could most plausibly get wrong, and the test above drives the path the cursor replaces.
    ///
    /// **A sample is not one read group.** One individual sequenced as several libraries
    /// declares an `@RG` for each, and every one of them is ours: the reads must all be kept,
    /// each carrying the library it came from, because that is the grouping a per-chemistry
    /// error model is fitted on. "Belongs to another sample" is the only thing that makes a
    /// record foreign.
    ///
    /// That matters more on the CRAM arm than anywhere else, because it is the one arm that
    /// decides a record's read group **itself**, while decoding a container — so a mistake
    /// here would drop a whole library's reads before any other layer could see them, and the
    /// output would look like a sample that was simply sequenced less deeply.
    #[test]
    fn a_cursor_keeps_every_read_group_of_its_sample_not_just_one() {
        use crate::ng::read::input::test_fixtures::{
            indexed_cram_declaring, read_named_with_length_in_read_group,
        };
        use crate::ng::reference_info::{ReferenceSource, read_reference_info};

        let (_cram_dir, cram_path, _fasta_dir, fasta) = indexed_cram_declaring(
            &[
                read_named_with_length_in_read_group("a", 0, 1, 30, "rg1"),
                read_named_with_length_in_read_group("b", 0, 20, 30, "rg2"),
                read_named_with_length_in_read_group("c", 0, 40, 30, "rg1"),
            ],
            // Two libraries, **one sample**. Not two samples — that is the other case, and it
            // is what `a_shared_cram_serves_each_open_only_its_own_reads` covers.
            &[("rg1", Some("NA12878")), ("rg2", Some("NA12878"))],
        );
        let reference = OpenReference::from(
            read_reference_info(ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            })
            .expect("read reference"),
        );

        let read_groups = build_read_groups(std::slice::from_ref(&cram_path)).expect("one sample");
        let [sample] = read_groups.read_groups_per_sample() else {
            panic!("both read groups name NA12878");
        };
        let reads = SampleReads::open(
            sample,
            &read_groups,
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("opens");

        let mut cursor = reads
            .cursor(ContigId(0), reference_bases)
            .expect("a cursor for contig 0");
        cursor
            .move_to_region(whole_first_contig())
            .expect("on this chromosome");

        let mut observed = Vec::new();
        while let Some(read) = cursor.next_read() {
            let read = read.expect("every record resolves");
            observed.push((
                String::from_utf8_lossy(&read.qname).into_owned(),
                read.read_group.get(),
            ));
        }

        assert_eq!(
            observed,
            vec![
                ("a".to_string(), 0),
                ("b".to_string(), 1),
                ("c".to_string(), 0),
            ],
            "every library of this sample is kept, each read carrying its own — a cursor that \
             kept only the first read group would drop 'b' and look like shallower sequencing",
        );
    }

    // --- how an open resolves one file's records ---

    /// `resolution_for` is where a file's read groups become this open's view of
    /// them, and every per-record decision follows from what it returns. It had
    /// been covered only sideways, through fixtures that happened to be untagged.
    #[test]
    fn one_own_read_group_resolves_without_consulting_any_record() {
        use crate::ng::read::input::test_fixtures::{FixtureReadGroup, header_with_read_groups};

        let (_dir, path) = named_bam(
            &header_with_read_groups(
                Some("coordinate"),
                &matching_contigs(),
                &[FixtureReadGroup::new("rg1", Some("NA12878"))],
            ),
            &one_read(),
            "only.bam",
        );
        let read_groups = build_read_groups(std::slice::from_ref(&path)).expect("valid");
        let sample = &read_groups.read_groups_per_sample()[0];

        let resolution = resolution_for(&path, sample, &read_groups).expect("resolves");

        assert_eq!(resolution, ReadGroupResolution::Sole(ReadGroupId(0)));
    }

    /// Several read groups, all this sample's: each record's own tag decides,
    /// and all of them are ours.
    #[test]
    fn several_own_read_groups_resolve_per_record_and_all_are_mine() {
        use crate::ng::read::input::test_fixtures::{FixtureReadGroup, header_with_read_groups};

        let (_dir, path) = named_bam(
            &header_with_read_groups(
                Some("coordinate"),
                &matching_contigs(),
                &[
                    FixtureReadGroup::new("rg1", Some("NA12878")),
                    FixtureReadGroup::new("rg2", Some("NA12878")),
                ],
            ),
            &one_read(),
            "two.bam",
        );
        let read_groups = build_read_groups(std::slice::from_ref(&path)).expect("valid");
        let sample = &read_groups.read_groups_per_sample()[0];

        let resolution = resolution_for(&path, sample, &read_groups).expect("resolves");

        assert_eq!(
            resolution,
            ReadGroupResolution::PerRecord(
                vec![
                    ("rg1".into(), RecordOwner::Mine(ReadGroupId(0))),
                    ("rg2".into(), RecordOwner::Mine(ReadGroupId(1))),
                ]
                .into()
            )
        );
    }

    /// A shared file: this open's read groups are `Mine`, the rest are not — and
    /// getting that backwards would yield another sample's reads as ours,
    /// carrying their read group.
    #[test]
    fn a_shared_file_resolves_only_this_samples_read_groups_as_mine() {
        use crate::ng::read::input::test_fixtures::{FixtureReadGroup, header_with_read_groups};

        let (_dir, path) = named_bam(
            &header_with_read_groups(
                Some("coordinate"),
                &matching_contigs(),
                &[
                    FixtureReadGroup::new("rg1", Some("NA12878")),
                    FixtureReadGroup::new("rg2", Some("HG002")),
                ],
            ),
            &one_read(),
            "shared-resolution.bam",
        );
        let read_groups = build_read_groups(std::slice::from_ref(&path)).expect("valid");

        for (index, expected) in [
            (
                0,
                [RecordOwner::Mine(ReadGroupId(0)), RecordOwner::OtherSample],
            ),
            (
                1,
                [RecordOwner::OtherSample, RecordOwner::Mine(ReadGroupId(1))],
            ),
        ] {
            let sample = &read_groups.read_groups_per_sample()[index];
            let resolution = resolution_for(&path, sample, &read_groups).expect("resolves");

            let ReadGroupResolution::PerRecord(entries) = resolution else {
                panic!("a shared file cannot resolve to a single read group");
            };
            let owners: Vec<RecordOwner> = entries.iter().map(|(_, owner)| *owner).collect();
            assert_eq!(owners, expected, "sample {}", sample.sample);
        }
    }

    /// A file the run's table never saw. Unreachable through the ordinary path —
    /// a file with no `@RG` fails when the table is built — but a caller
    /// assembling an open by hand gets told, rather than a resolution that
    /// rejects every record it later reads.
    #[test]
    fn a_file_absent_from_the_table_is_refused_at_open() {
        use crate::ng::read::input::test_fixtures::{FixtureReadGroup, header_with_read_groups};

        let (_dir, known) = named_bam(
            &header_with_read_groups(
                Some("coordinate"),
                &matching_contigs(),
                &[FixtureReadGroup::new("rg1", Some("NA12878"))],
            ),
            &one_read(),
            "known.bam",
        );
        let read_groups = build_read_groups(std::slice::from_ref(&known)).expect("valid");
        let sample = &read_groups.read_groups_per_sample()[0];

        let error = resolution_for(Path::new("/elsewhere/unknown.bam"), sample, &read_groups)
            .expect_err("the table holds no read group for that file");

        assert!(matches!(error, IngestError::NoReadGroupsForFile { .. }));
    }

    /// **The capability this step adds.** One file, two samples, opened twice:
    /// each open yields only its own reads, the foreign ones are counted apart
    /// from every drop, and the two opens' reads reunited are the file.
    ///
    /// Counting a foreign read as a drop is the failure this guards: it would
    /// make a shared file look like a low-quality one, and a drop rate is what a
    /// per-read-group diagnostic is *for*.
    #[test]
    fn a_file_shared_between_samples_serves_each_open_only_its_own_reads() {
        use crate::ng::read::input::test_fixtures::{
            FixtureReadGroup, header_with_read_groups, indexed_named_bam,
            read_named_with_length_in_read_group,
        };

        let (_reference_dir, reference) = fixture_reference(false);
        let header = header_with_read_groups(
            Some("coordinate"),
            &matching_contigs(),
            &[
                FixtureReadGroup::new("rg1", Some("NA12878")),
                FixtureReadGroup::new("rg2", Some("HG002")),
            ],
        );
        let (_dir, path) = indexed_named_bam(
            &header,
            &[
                read_named_with_length_in_read_group("a", 0, 1, 30, "rg1"),
                read_named_with_length_in_read_group("b", 0, 20, 30, "rg2"),
                read_named_with_length_in_read_group("c", 0, 40, 30, "rg1"),
                read_named_with_length_in_read_group("d", 0, 60, 30, "rg2"),
            ],
            "shared.bam",
        );

        let read_groups = build_read_groups(std::slice::from_ref(&path)).expect("two samples");
        assert_eq!(read_groups.read_groups_per_sample().len(), 2);

        let read_names_of = |index: usize| -> (Vec<String>, u64, u64) {
            let sample = &read_groups.read_groups_per_sample()[index];
            let reads = SampleReads::open(
                sample,
                &read_groups,
                &reference,
                ReadFilterConfig::default(),
                false,
            )
            .expect("a shared file opens for each of its samples");

            let mut cursor = cursor_over_first_contig(&reads);
            let names: Vec<String> = collect_reads(&mut cursor)
                .into_iter()
                .map(|(qname, _)| qname)
                .collect();
            let all = cursor.read_group_counts();
            let counts = &only_tally(&all);
            let dropped = counts.duplicate
                + counts.low_mapq
                + counts.supplementary
                + counts.secondary
                + counts.unmapped
                + counts.qc_fail
                + counts.too_short
                + counts.high_mismatch_fraction
                + counts.bad_cigar;
            (names, counts.other_sample, dropped)
        };

        let (first, first_foreign, first_dropped) = read_names_of(0);
        let (second, second_foreign, second_dropped) = read_names_of(1);

        assert_eq!(first, vec!["a".to_string(), "c".to_string()]);
        assert_eq!(second, vec!["b".to_string(), "d".to_string()]);
        assert_eq!(
            (first_foreign, second_foreign),
            (2, 2),
            "each open skipped the other sample's two reads"
        );
        assert_eq!(
            (first_dropped, second_dropped),
            (0, 0),
            "and counted none of them as a quality drop"
        );
    }

    /// The same shared-file behaviour over a **CRAM**, which reaches the other record reader
    /// — and with it the enum's CRAM arm, whose forwarding of the skipped-record count was
    /// added blind.
    #[test]
    fn a_shared_cram_serves_each_open_only_its_own_reads() {
        use crate::ng::read::input::test_fixtures::{
            indexed_cram_declaring, read_named_with_length_in_read_group,
        };
        use crate::ng::reference_info::{ReferenceSource, read_reference_info};

        let (_cram_dir, cram_path, _fasta_dir, fasta) = indexed_cram_declaring(
            &[
                read_named_with_length_in_read_group("a", 0, 1, 30, "rg1"),
                read_named_with_length_in_read_group("b", 0, 20, 30, "rg2"),
                read_named_with_length_in_read_group("c", 0, 40, 30, "rg1"),
            ],
            &[("rg1", Some("NA12878")), ("rg2", Some("HG002"))],
        );
        let reference = OpenReference::from(
            read_reference_info(ReferenceSource::Fasta {
                fasta: fasta.clone(),
                fai: None,
            })
            .expect("read reference"),
        );

        let read_groups = build_read_groups(std::slice::from_ref(&cram_path)).expect("two samples");
        assert_eq!(read_groups.read_groups_per_sample().len(), 2);

        let sample = &read_groups.read_groups_per_sample()[0];
        let reads = SampleReads::open(
            sample,
            &read_groups,
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("a shared CRAM opens for one of its samples");

        let mut cursor = cursor_over_first_contig(&reads);
        let names: Vec<String> = collect_reads(&mut cursor)
            .into_iter()
            .map(|(qname, _)| qname)
            .collect();

        assert_eq!(names, vec!["a".to_string(), "c".to_string()]);
        assert_eq!(
            only_tally(&cursor.read_group_counts()).other_sample,
            1,
            "the other sample's single read is reported through the record reader's CRAM arm"
        );
    }

    /// **The whole stack, on the shape none of the narrower tests has.** Two
    /// files and two samples at once, where one sample spans both files and the
    /// other lives inside a file it shares — so every seam this work touched is
    /// crossed in one pass: the table built from several paths, an open that
    /// merges two files, another that reads one sample out of a shared file,
    /// per-record resolution, the foreign-record skip, and the tally keyed by
    /// read group.
    ///
    /// It lives here rather than in `tests/` because the fixtures that write a
    /// BAM and its index are `pub(crate)`; making them public to host one test
    /// would widen the crate's surface for a test's convenience.
    ///
    /// Each of the three silent failures this milestone hit would fail here:
    /// resolving against a stale buffer (reads carry a neighbour's group),
    /// losing a skipped record from the tally (`other_sample` short), and
    /// charging a foreign read to a drop counter (`kept + drops` disagree).
    #[test]
    fn the_whole_stack_over_two_files_and_two_samples() {
        use crate::ng::read::input::test_fixtures::{
            FixtureReadGroup, header_with_read_groups, indexed_named_bam,
            read_named_with_length_in_read_group,
        };

        let (_reference_dir, reference) = fixture_reference(false);

        // File one is shared: NA12878's rg1 and HG002's rg2.
        let (_first_dir, first) = indexed_named_bam(
            &header_with_read_groups(
                Some("coordinate"),
                &matching_contigs(),
                &[
                    FixtureReadGroup::new("rg1", Some("NA12878")).with_library("lib-A"),
                    FixtureReadGroup::new("rg2", Some("HG002")).with_library("lib-B"),
                ],
            ),
            &[
                read_named_with_length_in_read_group("shared-a", 0, 1, 30, "rg1"),
                read_named_with_length_in_read_group("shared-b", 0, 20, 30, "rg2"),
                read_named_with_length_in_read_group("shared-c", 0, 40, 30, "rg1"),
            ],
            "first.bam",
        );
        // File two is NA12878's alone — a second library of the same individual,
        // which is the case the whole design exists for.
        let (_second_dir, second) = indexed_named_bam(
            &header_with_read_groups(
                Some("coordinate"),
                &matching_contigs(),
                &[FixtureReadGroup::new("rg1", Some("NA12878")).with_library("lib-C")],
            ),
            &[read_named_with_length_in_read_group(
                "second-a", 0, 10, 30, "rg1",
            )],
            "second.bam",
        );

        let read_groups =
            build_read_groups(&[first.clone(), second.clone()]).expect("both headers are valid");

        // Three read groups over two files, two samples, three libraries — and
        // the identifiers follow the input order, not the samples'.
        assert_eq!(read_groups.len(), 3);
        assert_eq!(
            read_groups
                .iter()
                .map(|(id, group)| (id.get(), &*group.sample, &*group.library.value))
                .collect::<Vec<_>>(),
            vec![
                (0, "NA12878", "lib-A"),
                (1, "HG002", "lib-B"),
                (2, "NA12878", "lib-C"),
            ]
        );

        let by_sample = read_groups.read_groups_per_sample();
        assert_eq!(
            by_sample
                .iter()
                .map(|entry| (&*entry.sample, entry.read_groups.len()))
                .collect::<Vec<_>>(),
            vec![("NA12878", 2), ("HG002", 1)],
            "one sample spans two files; the other is inside a file it shares"
        );

        let drain = |index: usize| -> (Vec<String>, Vec<(u32, u64)>, u64) {
            let sample = &by_sample[index];
            let reads = SampleReads::open(
                sample,
                &read_groups,
                &reference,
                ReadFilterConfig::default(),
                false,
            )
            .expect("opens");

            let mut cursor = cursor_over_first_contig(&reads);
            let names: Vec<String> = collect_reads(&mut cursor)
                .into_iter()
                .map(|(qname, _)| qname)
                .collect();

            let counts = cursor.read_group_counts();
            let kept: Vec<(u32, u64)> = counts
                .iter()
                .filter_map(|(id, tally)| Some((id.as_ref()?.get(), tally.kept)))
                .collect();
            let foreign = counts.iter().map(|(_, tally)| tally.other_sample).sum();
            (names, kept, foreign)
        };

        let (first_names, first_kept, first_foreign) = drain(0);
        let (second_names, second_kept, second_foreign) = drain(1);

        // NA12878's reads come from both files, merged in coordinate order.
        assert_eq!(
            first_names,
            vec![
                "shared-a".to_string(),
                "second-a".to_string(),
                "shared-c".to_string(),
            ],
        );
        assert_eq!(
            first_kept,
            vec![(0, 2), (2, 1)],
            "two libraries, tallied apart — never summed to one count of three"
        );
        assert_eq!(first_foreign, 1, "HG002's read, skipped in the shared file");

        assert_eq!(second_names, vec!["shared-b".to_string()]);
        assert_eq!(second_kept, vec![(1, 1)]);
        assert_eq!(
            second_foreign, 2,
            "NA12878's two reads in the shared file, skipped"
        );

        // Nothing was charged to a drop counter: every record either belonged to
        // the open or was skipped as another sample's.
        for sample in by_sample {
            let reads = SampleReads::open(
                sample,
                &read_groups,
                &reference,
                ReadFilterConfig::default(),
                false,
            )
            .expect("opens");
            let mut cursor = cursor_over_first_contig(&reads);
            while let Some(read) = cursor.next_read() {
                read.expect("resolves");
            }
            for (_, tally) in cursor.read_group_counts() {
                assert_eq!(
                    tally.duplicate
                        + tally.low_mapq
                        + tally.supplementary
                        + tally.secondary
                        + tally.unmapped
                        + tally.qc_fail
                        + tally.too_short
                        + tally.high_mismatch_fraction
                        + tally.bad_cigar,
                    0,
                    "a foreign read is never a quality drop"
                );
            }
        }
    }

    /// A read must carry the read group it came from. It is load-bearing rather
    /// than provenance — a read group is one library preparation, the grouping a
    /// per-chemistry error model keys on — and nothing proved it reached the
    /// reads until the reads could be read.
    ///
    /// The fixtures make each file its own read group, so "which read group" and
    /// "which file" coincide here. They do not in general: the identifier is
    /// run-wide, and one file can hold several read groups.
    #[test]
    fn each_read_carries_the_read_group_it_came_from() {
        let (_first_dir, first) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &[
                read_named_with_length("a1", 0, 1, 30),
                read_named_with_length("a2", 0, 50, 30),
            ],
            "first.bam",
        );
        let (_second_dir, second) = indexed_named_bam(
            &bam_header(&matching_contigs()),
            &[read_named_with_length("b1", 0, 25, 30)],
            "second.bam",
        );

        let (_dir, sample) = sample_over(&[first, second]);
        let reads = drain(&sample);

        assert_eq!(
            reads,
            vec![
                ("a1".to_string(), 0),
                ("b1".to_string(), 1),
                ("a2".to_string(), 0),
            ],
            "every read tagged with its own file, through the merge"
        );
    }

    /// A per-file failure is wrapped with the index of the file that raised it, so a caller
    /// matching on `IngestError::File` learns *which* file was at fault. With one file that
    /// index is 0 — trivially, which is the point: the wrapping is unconditional rather than
    /// something the merged arm does and the single one skips.
    #[test]
    fn the_single_file_arm_wraps_its_errors_like_the_merged_one() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_dir, path) = indexed_bam(
            &bam_header(&matching_contigs()),
            &[read_named_with_length("r", 0, 1, 30)],
        );
        let sample =
            SampleReads::open_only_sample(&[path], &reference, ReadFilterConfig::default(), false)
                .expect("opens");

        // A chromosome the file does not have.
        let error = sample
            .cursor(ContigId(9), reference_bases)
            .err()
            .expect("a contig the file lacks must not yield a cursor");
        match error {
            IngestError::File {
                source_file_index, ..
            } => assert_eq!(source_file_index, 0),
            other => panic!("expected a wrapped per-file error, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------
    // check_assembly (D1) — T2b
    // -----------------------------------------------------------------
    use crate::ng::reference_info::ContigInfo;

    fn contig(name: &str, md5: Option<[u8; 16]>) -> ContigInfo {
        ContigInfo {
            name: name.to_string(),
            length: 100,
            offset: 0,
            line_bases: 60,
            line_width: 61,
            md5,
        }
    }

    fn verified_reference(digests: &[Option<[u8; 16]>]) -> ReferenceInfo {
        ReferenceInfo {
            md5: None,
            contigs: digests
                .iter()
                .enumerate()
                .map(|(i, md5)| contig(&format!("chr{}", i + 1), *md5))
                .collect(),
            fasta_path: None,
        }
    }

    const A: [u8; 16] = [0xaa; 16];
    const B: [u8; 16] = [0xbb; 16];

    /// **T2b, the fault it exists to catch.** A file aligned to a different
    /// assembly with the same contig names and lengths passes every check at
    /// open — that is exactly why this check is worth having, and why
    /// production, which trusts the `@SQ M5`, cannot make it.
    #[test]
    fn t2b_a_disagreeing_digest_names_the_contig_and_both_values() {
        let error = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(B)],
            &verified_reference(&[Some(A), Some(A)]),
        )
        .expect_err("the second contig disagrees");

        assert_eq!(error.contig, "chr2");
        assert_eq!(error.expected, A);
        assert_eq!(error.observed, B);

        // Asserted whole, because the mistake worth catching is transposing
        // `expected` and `observed` in the message: that inverts the diagnosis
        // and sends the user after the wrong artefact, while leaving every
        // "contains" assertion green.
        assert_eq!(
            error.to_string(),
            format!(
                "alignment file '/data/sample.bam' is aligned to a different \
                 assembly: contig 'chr2' has M5 {} but the reference's is {}",
                "bb".repeat(16),
                "aa".repeat(16)
            )
        );
    }

    /// **T2b, the ordinary case.** A file with no `M5` tags at all is not a
    /// fault — it simply cannot be checked. `compared == 0` says so plainly
    /// rather than implying a guarantee that was never obtained.
    #[test]
    fn t2b_a_file_without_digests_is_not_an_error_and_reports_nothing_compared() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[None, None, None],
            &verified_reference(&[Some(A), Some(A), Some(A)]),
        )
        .expect("no tags is ordinary input");

        assert_eq!(
            check,
            AssemblyCheck {
                compared: 0,
                total: 3
            }
        );
    }

    /// **T2b, the partial case.** `compared` counts exactly the contigs that
    /// could be checked — so "verified 2 of 3" is a statement a caller can make
    /// truthfully.
    #[test]
    fn t2b_compared_counts_exactly_the_contigs_both_sides_carry() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), None, Some(A)],
            &verified_reference(&[Some(A), Some(A), Some(A)]),
        )
        .expect("the tagged contigs agree");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 2,
                total: 3
            }
        );

        // And symmetrically: a digest on the file's side is useless if the
        // reference has none, which is the whole `.fai`-only situation.
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(A), Some(A)],
            &verified_reference(&[Some(A), None, None]),
        )
        .expect("only the first contig is comparable");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 1,
                total: 3
            }
        );
    }

    /// The headline case a caller wants to report: everything tagged, all
    /// agreeing, so `compared == total`.
    #[test]
    fn a_fully_tagged_agreeing_file_reports_every_contig_verified() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(A), Some(A)],
            &verified_reference(&[Some(A), Some(A), Some(A)]),
        )
        .expect("all agree");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 3,
                total: 3
            }
        );
    }

    /// A `.fai`-only reference carries no digests at all, so nothing is
    /// comparable and nothing is claimed — the case that makes this check
    /// *deferred* rather than redundant with the one at open.
    #[test]
    fn nothing_is_comparable_against_a_fai_only_reference() {
        let check = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(A), Some(B)],
            &verified_reference(&[None, None]),
        )
        .expect("a .fai has no digests to disagree with");
        assert_eq!(
            check,
            AssemblyCheck {
                compared: 0,
                total: 2
            }
        );
    }

    /// The first disagreement wins, so the message names the earliest contig
    /// rather than an arbitrary one.
    #[test]
    fn the_first_disagreeing_contig_is_the_one_reported() {
        let error = check_assembly(
            Path::new("/data/sample.bam"),
            &[Some(B), Some(B)],
            &verified_reference(&[Some(A), Some(A)]),
        )
        .expect_err("both disagree");
        assert_eq!(error.contig, "chr1");
    }

    /// The variant carries `files` as well as `names`, and the whole point of
    /// the check is that someone picked up a stray file — so the message must
    /// name the *path*, not just repeat the sample names the user already
    /// knows disagree.
    #[test]
    fn sample_name_mismatch_message_names_the_offending_paths() {
        let message = IngestError::SampleNameMismatch {
            files: vec![PathBuf::from("/data/a.bam"), PathBuf::from("/data/b.bam")],
            names: vec!["NA12878".to_string(), "NA12892".to_string()],
        }
        .to_string();

        assert_eq!(
            message,
            "the files of one sample name different samples: \
             '/data/a.bam' names 'NA12878', '/data/b.bam' names 'NA12892'"
        );
    }

    /// **`SampleReads::cursor` must build one cursor per file, and nothing tested that.**
    ///
    /// Changing its loop to `.take(1)` — building only the first file's cursor — left the
    /// whole suite green. A sample sequenced across several runs would then be genotyped from
    /// one of its BAMs, silently: the reads it returned would be real, and the ones it never
    /// looked for invisible.
    #[test]
    fn a_sample_cursor_covers_every_file_of_the_sample() {
        use crate::ng::read::input::test_fixtures::{
            bam_header, fixture_reference, indexed_named_bam, matching_contigs,
            read_named_with_length,
        };
        use crate::ng::ref_seq::InMemoryRefSeq;
        use crate::ng::types::{ContigId, Position};

        let (_reference_dir, reference) = fixture_reference(false);
        let header = bam_header(&matching_contigs());
        // Two files, disjoint reads, so a cursor that skipped one would still answer
        // plausibly — which is the point.
        let (_dir_a, path_a) = indexed_named_bam(
            &header,
            &[read_named_with_length("from-a", 0, 1, 30)],
            "a.bam",
        );
        let (_dir_b, path_b) = indexed_named_bam(
            &header,
            &[read_named_with_length("from-b", 0, 41, 30)],
            "b.bam",
        );

        let sample = SampleReads::open_only_sample(
            &[path_a, path_b],
            &reference,
            ReadFilterConfig::default(),
            true,
        )
        .expect("two files naming one sample open together");
        assert_eq!(sample.file_count(), 2);

        let bases = || {
            InMemoryRefSeq::from_contigs(
                crate::ng::read::input::test_fixtures::FIXTURE_CONTIGS
                    .iter()
                    .map(|(_, length)| vec![b'A'; *length])
                    .collect(),
            )
        };
        let mut cursor = sample
            .cursor(ContigId(0), bases)
            .expect("a cursor per file of this sample");
        cursor
            .move_to_region(GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(100),
            })
            .expect("on this chromosome");

        let mut names = Vec::new();
        while let Some(read) = cursor.next_read() {
            names.push(
                String::from_utf8_lossy(&read.expect("the fixture reads decode").qname)
                    .into_owned(),
            );
        }
        assert_eq!(
            names,
            ["from-a", "from-b"],
            "a sample's cursor must reach every file it was opened with",
        );
    }
}
