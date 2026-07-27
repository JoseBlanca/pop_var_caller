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

use noodles_sam as sam;

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
// Reading one file's header
// ---------------------------------------------------------------------

/// What one `@RG` record actually said, before any missing name is filled in.
///
/// Private and intermediate, because reading a header (what the file *says*) and
/// filling in the names it left out are separate jobs — and only the second one
/// needs to know the file's name.
///
/// **No experiment tag is read.** `SRX` is not a tag the SAM specification
/// defines, so which one to look for is still open (spec §13); until it closes,
/// every experiment name is synthesized from the library.
#[derive(Debug, Clone, PartialEq, Eq)]
struct DeclaredReadGroup {
    id: Box<str>,
    sample: Box<str>,
    library: Option<Box<str>>,
    platform: Option<Box<str>>,
}

/// Every `@RG` record of one file, as its header declares them.
///
/// Both hard errors live here because both are properties of a single header
/// (spec §6). `path` is carried only so a message can name the file — a user who
/// hits either of these has to go and re-header something, and the path is what
/// tells them which.
///
/// **The first fault in header order is the one reported**, the rule the code
/// this replaces already followed: both are fatal and both are fixed by
/// correcting the header, so ranking them would add a rule to explain without
/// changing what the user does about it.
///
/// Several read groups are normal, and so are several *samples* among them: a
/// file naming two samples is no longer an error here. Keeping one sample per
/// open is a check on the open, not on the file (spec §4, §8).
// Only this module's tests call it until `build_read_groups` does. Two things
// keep the attribute from becoming a licence: `expect` (not `allow`) turns into
// an error the moment the function *is* called from the library, and `not(test)`
// scopes it to the build where it really is uncalled. Silencing it here also
// makes it a live root, so `DeclaredReadGroup` and `owned_str` need nothing.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "called only by this module's tests until build_read_groups consumes it"
    )
)]
fn declared_read_groups(
    header: &sam::Header,
    path: &Path,
) -> Result<Vec<DeclaredReadGroup>, ReadGroupError> {
    use sam::header::record::value::map::read_group::tag::{LIBRARY, PLATFORM, SAMPLE};

    if header.read_groups().is_empty() {
        return Err(ReadGroupError::NoReadGroups {
            path: path.to_path_buf(),
        });
    }

    let mut declared = Vec::with_capacity(header.read_groups().len());

    for (raw_id, read_group) in header.read_groups() {
        let id = String::from_utf8_lossy(raw_id.as_ref()).into_owned();
        let tag = |tag| {
            read_group
                .other_fields()
                .get(&tag)
                .map(|raw| owned_str(raw.as_ref()))
        };

        let Some(sample) = tag(SAMPLE) else {
            return Err(ReadGroupError::MissingSampleName {
                path: path.to_path_buf(),
                read_group_id: id,
            });
        };

        declared.push(DeclaredReadGroup {
            id: id.into_boxed_str(),
            sample,
            library: tag(LIBRARY),
            platform: tag(PLATFORM),
        });
    }

    Ok(declared)
}

/// A header tag's bytes as an owned string. Lossy, deliberately: a header tag
/// with invalid UTF-8 is a broken header, but it is a *label*, and refusing to
/// read the file over it would be a worse outcome than carrying the replacement
/// character into a report.
fn owned_str(raw: &[u8]) -> Box<str> {
    String::from_utf8_lossy(raw).into_owned().into_boxed_str()
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::input::test_fixtures::{
        FixtureReadGroup, header, header_with_read_groups, matching_contigs,
    };

    fn path() -> PathBuf {
        PathBuf::from("/data/project/SRR1.cram")
    }

    // --- reading one header's @RG records (A4) ---

    #[test]
    fn every_read_group_is_returned_in_header_order_with_its_tags() {
        let header = header_with_read_groups(
            Some("coordinate"),
            &matching_contigs(),
            &[
                FixtureReadGroup::new("rg1", Some("NA12878"))
                    .with_library("lib-A")
                    .with_platform("ILLUMINA"),
                FixtureReadGroup::new("rg2", Some("NA12878")).with_library("lib-B"),
            ],
        );

        let declared = declared_read_groups(&header, &path()).expect("both read groups are valid");

        assert_eq!(declared.len(), 2);
        assert_eq!(&*declared[0].id, "rg1", "header order, not sorted");
        assert_eq!(&*declared[0].sample, "NA12878");
        assert_eq!(declared[0].library.as_deref(), Some("lib-A"));
        assert_eq!(declared[0].platform.as_deref(), Some("ILLUMINA"));
        assert_eq!(&*declared[1].id, "rg2");
        assert_eq!(declared[1].library.as_deref(), Some("lib-B"));
        assert_eq!(
            declared[1].platform, None,
            "an absent tag is absent, not empty"
        );
    }

    #[test]
    fn a_file_with_no_read_groups_is_rejected_naming_the_file() {
        let header = header(Some("coordinate"), &matching_contigs(), &[]);

        let error = declared_read_groups(&header, &path()).expect_err("no @RG is fatal");

        assert!(matches!(error, ReadGroupError::NoReadGroups { .. }));
        let message = error.to_string();
        assert!(
            message.contains("/data/project/SRR1.cram"),
            "the message names the file: {message}"
        );
        assert!(
            message.contains("addreplacerg"),
            "the message says what to do about it: {message}"
        );
    }

    #[test]
    fn a_read_group_without_a_sample_is_rejected_naming_that_read_group() {
        let header = header(Some("coordinate"), &matching_contigs(), &[("rg1", None)]);

        let error = declared_read_groups(&header, &path()).expect_err("no SM is fatal");

        match &error {
            ReadGroupError::MissingSampleName { read_group_id, .. } => {
                assert_eq!(read_group_id, "rg1")
            }
            other => panic!("expected MissingSampleName, got {other:?}"),
        }
        let message = error.to_string();
        assert!(
            message.contains("@RG 'rg1'") && message.contains("/data/project/SRR1.cram"),
            "the message names both the read group and the file: {message}"
        );
    }

    /// **The first fault in header order wins.** A valid read group ahead of the
    /// offender does not shield it, and a later valid one does not rescue it.
    #[test]
    fn the_missing_sample_is_reported_for_whichever_read_group_lacks_it() {
        let second_is_bad = header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg1", Some("NA12878")), ("rg2", None)],
        );

        let error = declared_read_groups(&second_is_bad, &path()).expect_err("rg2 has no SM");

        match error {
            ReadGroupError::MissingSampleName { read_group_id, .. } => {
                assert_eq!(read_group_id, "rg2", "the offender, not the first record")
            }
            other => panic!("expected MissingSampleName, got {other:?}"),
        }
    }

    /// The behaviour change of spec §8: a file whose read groups name two
    /// samples used to be rejected at open. It is now an ordinary input — one
    /// sample per *open* is a separate check, made where the open is.
    #[test]
    fn a_file_naming_two_samples_is_read_without_complaint() {
        let two_samples = header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg1", Some("NA12878")), ("rg2", Some("HG002"))],
        );

        let declared = declared_read_groups(&two_samples, &path()).expect("no longer an error");

        assert_eq!(
            declared.iter().map(|rg| &*rg.sample).collect::<Vec<_>>(),
            vec!["NA12878", "HG002"]
        );
    }
}
