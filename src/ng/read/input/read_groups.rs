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

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use noodles_sam as sam;

use crate::ng::read::input::AlignmentFileError;
use crate::ng::read::input::open_bam::read_header;
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

// ---------------------------------------------------------------------
// Filling in the names a file did not declare
// ---------------------------------------------------------------------

/// Separates the parts of a synthesized library name. Distinct enough to read
/// back — `NA12878:rg1:SRR1` — and rare enough inside accessions and specimen
/// names that the parts stay legible.
const SYNTHESIZED_NAME_SEPARATOR: char = ':';

/// One declared `@RG` record, completed into a [`ReadGroup`]: whatever the file
/// left out is filled in here, and marked as ours.
///
/// The two fallbacks differ in what they protect against, so they are not
/// symmetric (spec §6):
///
/// - a missing **library** gets a name built from the sample, the `@RG ID` and
///   the file's name, because those three are what distinguishes one undeclared
///   read group from another;
/// - a missing **experiment** becomes the library, coarser rather than finer.
///   Several read groups legitimately share one declared library — that is a
///   library sequenced across lanes, one preparation, one experiment — and
///   falling back to `(library, id)` would split them. Falling back coarse is
///   safe here only because `id` survives on the record, so anything wanting the
///   finer split can still take it.
fn into_read_group(declared: DeclaredReadGroup, file: Arc<Path>) -> ReadGroup {
    let library = match declared.library {
        Some(value) => NameWithOrigin {
            value,
            origin: NameOrigin::Declared,
        },
        None => NameWithOrigin {
            value: synthesized_library(&declared.sample, &declared.id, &file),
            origin: NameOrigin::Synthesized,
        },
    };

    // Always synthesized: no experiment tag is read yet (spec §13). The origin
    // describes *this* name, so it is `Synthesized` even when the library it
    // copies was declared — we chose to call the experiment that, the file
    // did not.
    let experiment = NameWithOrigin {
        value: library.value.clone(),
        origin: NameOrigin::Synthesized,
    };

    ReadGroup {
        file,
        id: declared.id,
        sample: declared.sample,
        library,
        experiment,
        platform: declared.platform,
    }
}

/// A library name for a read group whose file declared none.
///
/// All three parts are needed. Without the file's name, two files that each
/// declare `ID:1` for the same sample would be given the same library and
/// silently fused; with it, the name is unique as long as no two input files
/// share a name, since `@RG ID` is unique within a file. Two files that *do*
/// share a name is the residual case, and it is a hard error raised where the
/// files meet, not papered over here (spec §6).
fn synthesized_library(sample: &str, read_group_id: &str, file: &Path) -> Box<str> {
    let stem = file_name_without_extensions(file);
    format!(
        "{sample}{sep}{read_group_id}{sep}{stem}",
        sep = SYNTHESIZED_NAME_SEPARATOR
    )
    .into_boxed_str()
}

/// The file's name with every extension removed: `SRR1.p1.cram` → `SRR1`.
///
/// Cut at the **first** dot rather than the last, because alignment files
/// routinely carry compound suffixes (`.p1.cram`, `.sorted.bam`) and it is the
/// part before the first dot that is the accession or specimen name a reader
/// recognises. Falls back to the whole file name, and then to the whole path,
/// rather than yielding an empty string for a dotfile or a path with no final
/// component — an empty part would make two names collide that are not the same.
fn file_name_without_extensions(file: &Path) -> String {
    let Some(name) = file.file_name().map(|name| name.to_string_lossy()) else {
        return file.to_string_lossy().into_owned();
    };
    match name.split_once('.') {
        Some((stem, _)) if !stem.is_empty() => stem.to_owned(),
        _ => name.into_owned(),
    }
}

/// A header tag's bytes as an owned string. Lossy, deliberately: a header tag
/// with invalid UTF-8 is a broken header, but it is a *label*, and refusing to
/// read the file over it would be a worse outcome than carrying the replacement
/// character into a report.
fn owned_str(raw: &[u8]) -> Box<str> {
    String::from_utf8_lossy(raw).into_owned().into_boxed_str()
}

// ---------------------------------------------------------------------
// Building the run's table
// ---------------------------------------------------------------------

/// Read every input file's header — **headers only** — and build the run's
/// read-group table.
///
/// This is not a preliminary step bolted on for an edge case: it is what decides
/// what the samples *are*. Nothing can build a per-sample file list without
/// reading the headers first, and doing it once here is also what keeps the
/// identifiers unique — a file opened later for two different samples must not
/// mint its read groups twice (spec §4, §5).
///
/// **Deterministic by construction.** Identifiers follow `paths` order and then
/// header order within a file, so the same input list always yields the same
/// identifiers. Nothing here depends on how fast a file was read, which is what
/// lets an identifier reach a report or an output file.
///
/// Opens no index and reads no record. An empty `paths` yields an empty table;
/// a file with no `@RG` is an error, not an empty contribution.
pub fn build_read_groups(paths: &[PathBuf]) -> Result<ReadGroups, ReadGroupError> {
    let mut read_groups: Vec<ReadGroup> = Vec::new();

    for path in paths {
        let header = read_header(path).map_err(|source| ReadGroupError::HeaderRead { source })?;
        let file: Arc<Path> = Arc::from(path.as_path());
        for declared in declared_read_groups(&header, path)? {
            read_groups.push(into_read_group(declared, Arc::clone(&file)));
        }
    }

    reject_colliding_synthesized_libraries(&read_groups)?;
    let per_sample = group_by_sample(&read_groups);

    Ok(ReadGroups {
        read_groups,
        per_sample,
    })
}

/// Fail if two read groups that had no `LB` ended up with the same library name.
///
/// **Only synthesized names are checked.** A *declared* library repeating is
/// ordinary and meaningful — that is one preparation sequenced across several
/// lanes or files, which is exactly the grouping we want to keep. Only names we
/// invented can collide by accident, and when they do the two read groups have
/// become indistinguishable, which no downstream step could undo (spec §6).
fn reject_colliding_synthesized_libraries(read_groups: &[ReadGroup]) -> Result<(), ReadGroupError> {
    let mut seen: HashMap<&str, &ReadGroup> = HashMap::new();

    for group in read_groups {
        if group.library.origin != NameOrigin::Synthesized {
            continue;
        }
        if let Some(first) = seen.insert(&group.library.value, group) {
            return Err(ReadGroupError::DuplicateSynthesizedLibrary {
                library: group.library.value.to_string(),
                paths: (first.file.to_path_buf(), group.file.to_path_buf()),
            });
        }
    }

    Ok(())
}

/// The same read groups, grouped by the sample they name, samples in first-seen
/// order and each sample's read groups in identifier order.
///
/// A linear scan over the samples found so far rather than a map, because that
/// *is* first-seen order — with a map the order would have to be reconstructed,
/// and the number of samples in a run is small enough that the scan costs
/// nothing.
fn group_by_sample(read_groups: &[ReadGroup]) -> Vec<SampleReadGroups> {
    let mut per_sample: Vec<SampleReadGroups> = Vec::new();

    for (index, group) in read_groups.iter().enumerate() {
        let id = ReadGroupId(u32::try_from(index).expect("a read-group table fits in u32"));
        match per_sample
            .iter_mut()
            .find(|entry| entry.sample == group.sample)
        {
            Some(entry) => entry.read_groups.push(id),
            None => per_sample.push(SampleReadGroups {
                sample: group.sample.clone(),
                read_groups: vec![id],
            }),
        }
    }

    per_sample
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
///
/// **Cloning is cheap, deliberately.** A region source takes a copy per query
/// and there are ~10⁶ queries in a run, so the sharing lives *inside* the one
/// variant that owns anything — `PerRecord`'s table, behind an `Arc` — rather
/// than around the whole enum. `Sole` is then a tag and a `u32`, which is the
/// overwhelmingly common case and costs no atomic at all; `PerRecord` costs one
/// increment. Wrapping the enum instead would have put a pointer chase and an
/// atomic on every query including the `Sole` ones, and an `Arc` in a struct
/// built a million times is worth avoiding (owner, 2026-07-28).
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
    ///
    /// `Arc<[_]>` rather than `Box<[_]>`: the table is settled at open and read
    /// but never written, and this is what makes a per-query clone one atomic
    /// increment instead of an allocation (see the type's doc).
    PerRecord(Arc<[(Box<str>, RecordOwner)]>),
}

impl ReadGroupResolution {
    /// Which read group a record belongs to, given the `RG` tag it carries (or
    /// does not).
    ///
    /// Pure, and takes the tag's bytes rather than a record, so the decision can
    /// be tested on its own — it is the point at which a read starts being
    /// attributed to a library, and a wrong answer here is silent.
    ///
    /// The `Sole` arm ignores `tag` entirely, which is the whole reason a file
    /// re-headered without rewriting its records can be read at all.
    pub fn owner_of(&self, tag: Option<&[u8]>) -> Result<RecordOwner, UnresolvedReadGroup> {
        match self {
            Self::Sole(id) => Ok(RecordOwner::Mine(*id)),
            Self::PerRecord(declared) => {
                let Some(tag) = tag else {
                    return Err(UnresolvedReadGroup::NoReadGroupTag);
                };
                declared
                    .iter()
                    .find(|(name, _)| name.as_bytes() == tag)
                    .map(|(_, owner)| *owner)
                    .ok_or_else(|| UnresolvedReadGroup::UndeclaredReadGroup {
                        named: String::from_utf8_lossy(tag).into_owned(),
                    })
            }
        }
    }
}

/// Why a record could not be assigned to a read group.
///
/// Only reachable in a file declaring several: with one, every record belongs to
/// it whatever its tag says. Both are fatal — with several candidates there is no
/// way to assign the record, and guessing would attribute a read to the wrong
/// library, which nothing downstream could detect.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum UnresolvedReadGroup {
    #[error("the record carries no RG tag, and its file declares several read groups")]
    NoReadGroupTag,
    #[error("the record names read group '{named}', which its file's header does not declare")]
    UndeclaredReadGroup { named: String },
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

    /// A file could not be opened, or its header could not be parsed.
    ///
    /// Wrapped rather than flattened, and the message stays short because the
    /// source already names the file and says how it failed. It exists so that
    /// "this input is unreadable" and "this input's read groups are wrong" stay
    /// distinguishable: they are different problems with different fixes.
    #[error("reading a file's read groups failed")]
    HeaderRead {
        #[source]
        source: AlignmentFileError,
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

    // --- filling in the names a file did not declare (A5) ---

    fn declared(id: &str, sample: &str, library: Option<&str>) -> DeclaredReadGroup {
        DeclaredReadGroup {
            id: id.into(),
            sample: sample.into(),
            library: library.map(Into::into),
            platform: None,
        }
    }

    fn file() -> Arc<Path> {
        Arc::from(path().as_path())
    }

    #[test]
    fn a_declared_library_is_carried_verbatim_and_marked_declared() {
        let group = into_read_group(declared("rg1", "NA12878", Some("lib-A")), file());

        assert_eq!(&*group.library.value, "lib-A");
        assert_eq!(group.library.origin, NameOrigin::Declared);
    }

    #[test]
    fn a_missing_library_is_built_from_sample_read_group_and_file_name() {
        let group = into_read_group(declared("rg1", "NA12878", None), file());

        assert_eq!(&*group.library.value, "NA12878:rg1:SRR1");
        assert_eq!(group.library.origin, NameOrigin::Synthesized);
    }

    /// Every extension goes, not just the last: `.p1.cram` is one compound
    /// suffix, and the accession before the first dot is what a reader knows the
    /// file by.
    #[test]
    fn the_file_name_loses_every_extension() {
        assert_eq!(
            file_name_without_extensions(Path::new("/data/SRR1.p1.cram")),
            "SRR1"
        );
        assert_eq!(
            file_name_without_extensions(Path::new("/data/HG002.sorted.bam")),
            "HG002"
        );
        assert_eq!(
            file_name_without_extensions(Path::new("/data/plain")),
            "plain"
        );
    }

    /// A dotfile has no part before its first dot. Falling back to the whole
    /// name matters: an empty part would make two names collide that are not the
    /// same read group at all.
    #[test]
    fn a_name_with_no_stem_falls_back_to_the_whole_name() {
        assert_eq!(
            file_name_without_extensions(Path::new("/data/.hidden.bam")),
            ".hidden.bam"
        );
    }

    /// The lanes of one preparation: several read groups, one declared library.
    /// They must stay one library — and one experiment — because that is what
    /// they are.
    #[test]
    fn read_groups_sharing_a_declared_library_share_one_library_and_experiment() {
        let lane_one = into_read_group(declared("rg1", "NA12878", Some("lib-A")), file());
        let lane_two = into_read_group(declared("rg2", "NA12878", Some("lib-A")), file());

        assert_eq!(lane_one.library.value, lane_two.library.value);
        assert_eq!(lane_one.experiment.value, lane_two.experiment.value);
        assert_ne!(lane_one.id, lane_two.id, "the finer split is still there");
    }

    /// Within one file the `@RG ID` is what separates undeclared libraries, so
    /// two of them cannot be fused.
    #[test]
    fn two_undeclared_libraries_in_one_file_get_different_names() {
        let first = into_read_group(declared("rg1", "NA12878", None), file());
        let second = into_read_group(declared("rg2", "NA12878", None), file());

        assert_ne!(first.library.value, second.library.value);
    }

    /// The experiment is always ours, even when it copies a declared library:
    /// the origin describes this name, not the one it was taken from.
    #[test]
    fn the_experiment_falls_back_to_the_library_and_says_it_is_synthesized() {
        let group = into_read_group(declared("rg1", "NA12878", Some("lib-A")), file());

        assert_eq!(group.experiment.value, group.library.value);
        assert_eq!(group.experiment.origin, NameOrigin::Synthesized);
        assert_eq!(
            group.library.origin,
            NameOrigin::Declared,
            "copying it does not make the library ours"
        );
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

    // --- resolving one record's read group ---

    fn per_record() -> ReadGroupResolution {
        ReadGroupResolution::PerRecord(
            vec![
                ("rg1".into(), RecordOwner::Mine(ReadGroupId(4))),
                ("rg2".into(), RecordOwner::Mine(ReadGroupId(5))),
            ]
            .into(),
        )
    }

    /// The universal path: one declared read group, and the tag is never
    /// consulted — which is what lets a file re-headered without rewriting its
    /// records be read at all.
    #[test]
    fn a_sole_read_group_claims_every_record_whatever_its_tag_says() {
        let sole = ReadGroupResolution::Sole(ReadGroupId(9));

        for tag in [None, Some(&b"rg1"[..]), Some(&b"nonsense"[..])] {
            assert_eq!(
                sole.owner_of(tag),
                Ok(RecordOwner::Mine(ReadGroupId(9))),
                "the tag is not read"
            );
        }
    }

    #[test]
    fn a_tagged_record_resolves_to_the_read_group_it_names() {
        assert_eq!(
            per_record().owner_of(Some(b"rg1")),
            Ok(RecordOwner::Mine(ReadGroupId(4)))
        );
        assert_eq!(
            per_record().owner_of(Some(b"rg2")),
            Ok(RecordOwner::Mine(ReadGroupId(5))),
            "not merely the first entry"
        );
    }

    /// With several candidates there is no way to assign an untagged record, and
    /// guessing would attribute a read to the wrong library — which nothing
    /// downstream could detect.
    #[test]
    fn an_untagged_record_is_unresolvable_when_several_are_declared() {
        assert_eq!(
            per_record().owner_of(None),
            Err(UnresolvedReadGroup::NoReadGroupTag)
        );
    }

    #[test]
    fn a_record_naming_an_undeclared_read_group_is_unresolvable() {
        let error = per_record()
            .owner_of(Some(b"rg7"))
            .expect_err("rg7 is not declared");

        assert_eq!(
            error,
            UnresolvedReadGroup::UndeclaredReadGroup {
                named: "rg7".to_string()
            }
        );
        assert!(
            error.to_string().contains("rg7"),
            "the message names it: {error}"
        );
    }

    /// A read group of another sample resolves — it is declared, after all — but
    /// resolves to *not ours*. Skipping it is a separate change; what matters
    /// here is that it is distinguishable from our own rather than silently
    /// counted as one.
    #[test]
    fn a_read_group_belonging_to_another_sample_resolves_as_not_ours() {
        let shared = ReadGroupResolution::PerRecord(
            vec![
                ("mine".into(), RecordOwner::Mine(ReadGroupId(1))),
                ("theirs".into(), RecordOwner::OtherSample),
            ]
            .into(),
        );

        assert_eq!(
            shared.owner_of(Some(b"theirs")),
            Ok(RecordOwner::OtherSample)
        );
    }
}

#[cfg(test)]
mod build_tests {
    use std::error::Error as _;

    use super::*;
    use crate::ng::read::input::test_fixtures::{
        FixtureReadGroup, header, header_with_read_groups, matching_contigs, named_bam, one_read,
    };
    use tempfile::TempDir;

    /// A BAM on disk whose header declares `read_groups`, under a chosen file
    /// name. The `TempDir` must be kept alive: dropping it removes the file.
    fn file_declaring(read_groups: &[FixtureReadGroup<'_>], file_name: &str) -> (TempDir, PathBuf) {
        let header = header_with_read_groups(Some("coordinate"), &matching_contigs(), read_groups);
        named_bam(&header, &one_read(), file_name)
    }

    #[test]
    fn identifiers_follow_input_file_order_then_header_order() {
        let (_first_dir, first) = file_declaring(
            &[
                FixtureReadGroup::new("rg1", Some("NA12878")),
                FixtureReadGroup::new("rg2", Some("NA12878")),
            ],
            "first.bam",
        );
        let (_second_dir, second) =
            file_declaring(&[FixtureReadGroup::new("rg1", Some("HG002"))], "second.bam");

        let table = build_read_groups(&[first, second]).expect("both headers are valid");

        assert_eq!(table.len(), 3);
        let named: Vec<(u32, &str, &str)> = table
            .iter()
            .map(|(id, group)| (id.get(), &*group.id, &*group.sample))
            .collect();
        assert_eq!(
            named,
            vec![
                (0, "rg1", "NA12878"),
                (1, "rg2", "NA12878"),
                (2, "rg1", "HG002"),
            ],
            "file order first, then header order within a file"
        );
    }

    #[test]
    fn read_groups_are_grouped_by_the_sample_they_name() {
        let (_first_dir, first) = file_declaring(
            &[
                FixtureReadGroup::new("rg1", Some("NA12878")),
                FixtureReadGroup::new("rg2", Some("HG002")),
            ],
            "first.bam",
        );
        let (_second_dir, second) = file_declaring(
            &[FixtureReadGroup::new("rg1", Some("NA12878"))],
            "second.bam",
        );

        let table = build_read_groups(&[first, second]).expect("both headers are valid");

        let grouped: Vec<(&str, Vec<u32>)> = table
            .read_groups_per_sample()
            .iter()
            .map(|entry| {
                (
                    &*entry.sample,
                    entry.read_groups.iter().map(|id| id.get()).collect(),
                )
            })
            .collect();
        assert_eq!(
            grouped,
            vec![("NA12878", vec![0, 2]), ("HG002", vec![1])],
            "samples in first-seen order; one sample's read groups may span files"
        );
    }

    /// The residual collision the file name cannot prevent: two inputs with the
    /// same name in different directories, whose read groups agree on sample and
    /// `@RG ID` and declare no library. They would become indistinguishable, so
    /// this is fatal rather than papered over.
    #[test]
    fn two_files_with_one_name_collide_on_a_synthesized_library() {
        let (_first_dir, first) =
            file_declaring(&[FixtureReadGroup::new("rg1", Some("NA12878"))], "run.bam");
        let (_second_dir, second) =
            file_declaring(&[FixtureReadGroup::new("rg1", Some("NA12878"))], "run.bam");

        let error = build_read_groups(&[first.clone(), second.clone()])
            .expect_err("the two synthesized libraries collide");

        match &error {
            ReadGroupError::DuplicateSynthesizedLibrary { library, paths } => {
                assert_eq!(library, "NA12878:rg1:run");
                assert_eq!(paths, &(first, second), "both paths, in input order");
            }
            other => panic!("expected DuplicateSynthesizedLibrary, got {other:?}"),
        }
        assert!(
            error.to_string().contains("rename one input"),
            "the message says what to do: {error}"
        );
    }

    /// A *declared* library repeating is the thing we want to keep, not a
    /// collision: one preparation sequenced into two files.
    #[test]
    fn a_declared_library_may_repeat_across_files() {
        let (_first_dir, first) = file_declaring(
            &[FixtureReadGroup::new("rg1", Some("NA12878")).with_library("lib-A")],
            "run.bam",
        );
        let (_second_dir, second) = file_declaring(
            &[FixtureReadGroup::new("rg1", Some("NA12878")).with_library("lib-A")],
            "run.bam",
        );

        let table = build_read_groups(&[first, second]).expect("a shared declared library is fine");

        assert_eq!(table.len(), 2);
        assert_eq!(
            table.get(ReadGroupId(0)).library.value,
            table.get(ReadGroupId(1)).library.value
        );
    }

    /// Identifiers are positions in the input list, so reordering the inputs
    /// reorders them — deliberately. What must not change is *which* read groups
    /// the table holds.
    #[test]
    fn reordering_the_inputs_reorders_the_identifiers_but_not_the_set() {
        let (_first_dir, first) = file_declaring(
            &[FixtureReadGroup::new("rg1", Some("NA12878"))],
            "first.bam",
        );
        let (_second_dir, second) =
            file_declaring(&[FixtureReadGroup::new("rg1", Some("HG002"))], "second.bam");

        let forward = build_read_groups(&[first.clone(), second.clone()]).expect("valid");
        let reversed = build_read_groups(&[second, first]).expect("valid");

        assert_eq!(&*forward.get(ReadGroupId(0)).sample, "NA12878");
        assert_eq!(&*reversed.get(ReadGroupId(0)).sample, "HG002");

        let samples = |table: &ReadGroups| {
            let mut names: Vec<String> = table
                .iter()
                .map(|(_, group)| group.sample.to_string())
                .collect();
            names.sort();
            names
        };
        assert_eq!(
            samples(&forward),
            samples(&reversed),
            "the same read groups"
        );
    }

    /// Building the table twice over the same list gives the same table — the
    /// property anything printing an identifier depends on.
    #[test]
    fn the_same_input_list_yields_the_same_table() {
        let (_dir, path) = file_declaring(
            &[
                FixtureReadGroup::new("rg1", Some("NA12878")).with_library("lib-A"),
                FixtureReadGroup::new("rg2", Some("HG002")),
            ],
            "run.bam",
        );

        let once = build_read_groups(std::slice::from_ref(&path)).expect("valid");
        let twice = build_read_groups(std::slice::from_ref(&path)).expect("valid");

        assert_eq!(once, twice);
    }

    #[test]
    fn a_header_that_cannot_be_read_is_its_own_failure() {
        let missing = PathBuf::from("/nonexistent/does-not-exist.bam");

        let error = build_read_groups(&[missing]).expect_err("the file is not there");

        assert!(matches!(error, ReadGroupError::HeaderRead { .. }));
        assert!(
            error.source().is_some(),
            "the underlying failure names the file"
        );
    }

    #[test]
    fn an_empty_input_list_yields_an_empty_table() {
        let table = build_read_groups(&[]).expect("nothing to read is not a failure");

        assert!(table.is_empty());
        assert!(table.read_groups_per_sample().is_empty());
    }

    /// The per-file failures reach the caller unchanged: a file with no `@RG`
    /// stops the whole pre-pass, because there is no partial table worth having.
    #[test]
    fn a_file_with_no_read_groups_stops_the_pre_pass() {
        let (_good_dir, good) =
            file_declaring(&[FixtureReadGroup::new("rg1", Some("NA12878"))], "good.bam");
        let empty_header = header(Some("coordinate"), &matching_contigs(), &[]);
        let (_bad_dir, bad) = named_bam(&empty_header, &one_read(), "bad.bam");

        let error = build_read_groups(&[good, bad]).expect_err("the second file has no @RG");

        assert!(matches!(error, ReadGroupError::NoReadGroups { .. }));
    }
}
