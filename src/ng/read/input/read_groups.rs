//! Read groups — the unit that carries the sequencing chemistry.
//!
//! A **read group** is one SAM `@RG` record: in practice one library
//! preparation sequenced in one run. It names the sample the DNA came from
//! (`SM`), the library it was prepared into (`LB`), and the platform. The
//! pipeline's unit has been the *sample*, but PCR stutter and per-base error are
//! properties of the library preparation and of the DNA's condition, not of the
//! individual — so the read group is what a per-chemistry error model keys on,
//! and it has to survive the read path for that model to be estimable at all.
//!
//! This module owns the read groups of a whole run: parsing them out of each
//! file's header, giving each one a run-wide [`ReadGroupId`], filling in the
//! names a file did not declare, and grouping them by sample so each open knows
//! which reads are its own.
//!
//! Design: `doc/devel/ng/spec/read_groups.md` (what and why),
//! `doc/devel/ng/arch/read_groups.md` (types and interfaces).

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ng::types::ReadGroupId;

// ---------------------------------------------------------------------
// The read group
// ---------------------------------------------------------------------

/// One `@RG` record: the file that declared it, what it says, and the grouping
/// names — declared or synthesized.
///
/// The **atoms** (`file`, `id`, `sample`) are read from the input and never
/// modified. `library` and `experiment` are grouping keys, computed from those
/// atoms when the file does not supply them. Because the atoms survive beside
/// them, no naming choice made here can destroy information: a consumer that
/// dislikes the grouping can always regroup on `(sample, id)` or on the file
/// (spec §3).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroup {
    /// The file that declared it. `Arc` because a file's k read groups share one
    /// path, and because `AlignmentFile` already holds its path that way.
    pub file: Arc<Path>,
    /// `@RG ID`, verbatim.
    ///
    /// A **label, never an identity**: the SAM specification makes it unique
    /// within its file and says nothing across files, so two files may each
    /// declare `ID:1`. Identity is the [`ReadGroupId`] this table mints
    /// (spec §4).
    pub id: Box<str>,
    /// `@RG SM` — the individual the DNA came from. Always present: a read group
    /// without it is [`ReadGroupError::MissingSampleName`].
    pub sample: Box<str>,
    /// `@RG LB`, or synthesized when the file declares none.
    pub library: NameWithOrigin,
    /// The sequencing experiment — one library preparation and its sequencing
    /// configuration. Falls back to the library, which is why several read
    /// groups sharing a declared library (the lanes of one preparation) share
    /// one experiment.
    pub experiment: NameWithOrigin,
    /// `@RG PL`. Carried for reports only — **nothing keys on it.** Real archives
    /// carry misspellings (`ilumina`) and values repurposed to mean something
    /// else, so it is not trustworthy enough to group by (spec §6).
    pub platform: Option<Box<str>>,
}

/// A name used for grouping, plus where it came from.
///
/// The origin cannot be recovered once the name exists, and any later report
/// about chemistry has to be able to say "this grouping is ours, not the
/// file's".
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NameWithOrigin {
    pub value: Box<str>,
    pub origin: NameOrigin,
}

/// Whether a grouping name came from the file or from us.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum NameOrigin {
    /// The file's own tag, verbatim.
    Declared,
    /// Absent from the file, so this module built one (spec §6).
    Synthesized,
}

// ---------------------------------------------------------------------
// The table
// ---------------------------------------------------------------------

/// Every read group in the run, in two views of the **same** set: by identifier,
/// and grouped by the sample each names.
///
/// Built once, before any file is opened for reading, then read-only and shared.
/// The by-sample view is not a second collection — it is these read groups
/// grouped by [`ReadGroup::sample`] — which is why it lives here rather than in
/// a wrapper holding both (arch §1.2).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReadGroups {
    /// Indexed by [`ReadGroupId`]: `read_groups[i]` has id `i`.
    read_groups: Vec<ReadGroup>,
    /// One entry per sample, in first-seen order.
    per_sample: Vec<SampleReadGroups>,
}

impl ReadGroups {
    /// One read group, by its identifier. Panics on an id this table never
    /// minted — the out-of-range check the unconstrained newtype defers to its
    /// lookup.
    pub fn get(&self, id: ReadGroupId) -> &ReadGroup {
        &self.read_groups[id.get() as usize]
    }

    /// How many read groups the run has.
    pub fn len(&self) -> usize {
        self.read_groups.len()
    }

    /// Whether the run has no read groups at all. Only reachable from an empty
    /// input list: a file with no `@RG` is an error, not an empty contribution.
    pub fn is_empty(&self) -> bool {
        self.read_groups.is_empty()
    }

    /// Every read group with its identifier, in table order.
    pub fn iter(&self) -> impl Iterator<Item = (ReadGroupId, &ReadGroup)> {
        self.read_groups.iter().enumerate().map(|(index, group)| {
            let id = u32::try_from(index).expect("a read-group table fits in u32");
            (ReadGroupId(id), group)
        })
    }

    /// The same read groups grouped by sample, one entry per sample in
    /// first-seen order — what each `SampleReads` open is built from.
    pub fn read_groups_per_sample(&self) -> &[SampleReadGroups] {
        &self.per_sample
    }
}

/// One sample and the read groups that name it.
///
/// The unit an open takes: a `SampleReads` serves exactly one sample, whichever
/// files its read groups live in and whatever else those files contain
/// (spec §4).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleReadGroups {
    pub sample: Box<str>,
    pub read_groups: Vec<ReadGroupId>,
}

// ---------------------------------------------------------------------
// Resolving a file's records
// ---------------------------------------------------------------------

/// How one open reads one file's records — decided when the file is opened,
/// never per record.
///
/// Two questions have to be answered about every record: which read group it
/// belongs to, and whether it belongs to the sample this open serves. This
/// answers both. It is built per open, not per file, because the second question
/// is about the sample (spec §7, §9).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadGroupResolution {
    /// The header declares exactly one read group, so every record is that one
    /// and the record's `RG` is **not read**. There is nothing to check: a read
    /// in such a file has one group it could belong to whatever its tag says —
    /// which is also what lets a file re-headered without rewriting its records
    /// be read at all. Such a file is single-sample by construction.
    Sole(ReadGroupId),
    /// The header declares several, so each record's `RG` names which. A record
    /// with no `RG`, or naming none of these, is fatal: with several candidates
    /// there is no way to assign it, and guessing would attribute a read to the
    /// wrong library.
    PerRecord(Box<[(Box<str>, RecordOwner)]>),
}

/// What a declared read group means to *this* open.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RecordOwner {
    /// This open's sample: the read is yielded, carrying this identifier.
    Mine(ReadGroupId),
    /// Declared in the file but naming another sample: the read is skipped and
    /// tallied **apart from the drop categories**, because a foreign read says
    /// nothing about how this read group behaved (spec §9).
    OtherSample,
}

// ---------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------

/// Header-level failures, all raised before any read flows.
///
/// Each names the file and says what to do about it, because every one of them
/// is fixed by re-headering an input rather than by changing anything here.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum ReadGroupError {
    /// The file declares no `@RG` record at all, so nothing says which sample or
    /// library its reads came from.
    #[error(
        "alignment file '{path}' declares no @RG record, so nothing says which sample \
         its reads came from; add one (e.g. `samtools addreplacerg`) and try again",
        path = path.display()
    )]
    NoReadGroups { path: PathBuf },

    /// An `@RG` record carries no `SM` tag. Distinct from
    /// [`NoReadGroups`](Self::NoReadGroups): the file has read groups, one of
    /// them just does not name its sample, and the fix is a different one.
    #[error(
        "@RG '{read_group_id}' in alignment file '{path}' has no SM tag, so it does not \
         name the sample its reads came from",
        path = path.display()
    )]
    MissingSampleName {
        path: PathBuf,
        read_group_id: String,
    },

    /// Two read groups without an `LB` tag would take the same synthesized
    /// library name. Only reachable when two input files share a file name — the
    /// same name in different directories — since `@RG ID` is unique within a
    /// file (spec §6).
    #[error(
        "two read groups carry no LB tag and would both be called library '{library}': \
         '{first}' and '{second}' have the same file name, so their read groups cannot \
         be told apart; rename one input or give the read groups LB tags",
        first = paths.0.display(),
        second = paths.1.display()
    )]
    DuplicateSynthesizedLibrary {
        library: String,
        paths: (PathBuf, PathBuf),
    },
}
