//! psp mode's calling stage: a cohort of stored samples, opened and checked before a block is
//! decoded.
//!
//! **The other half of spec §2's bargain.** The walk stage
//! ([`gatherer`](super::gatherer)) reads a sample's alignment files once and stores what it
//! saw; this is what reads a cohort of those files back and calls them. Everything between —
//! the merge, the genotyping, the records — is the body direct mode drives
//! ([`call_cohort_from_sources_handing_each_record_over`](super::callers::call_cohort_from_sources_handing_each_record_over)),
//! and the only thing that differs is where an observation comes from (spec §3.1).
//!
//! **Opening is two moves, not one, because the run's ground comes out of the files.** A
//! calling run over alignment files is told what to analyse and then opens its samples; a run
//! over psps is told nothing — the files know what ground they cover, and spec §5.3 says the
//! analysed regions come from the headers rather than from a flag. So:
//!
//! 1. [`OpenPspCohort::open`] opens every file, reads every header, and settles what the
//!    cohort *is*: the ground they agree on, one run-wide read-group numbering, and the
//!    refusals that need nothing but the headers (spec §6.2's cohort checks).
//! 2. Its caller builds the run's segmentation over that ground — through the same
//!    `run_ground::segments_over` both other subcommands use, which is what makes the check in
//!    move 3 mean anything.
//! 3. [`PspVariantCaller::open`] runs the refusals that compare each file *against the run*.
//!
//! **The split is what keeps the ground assembled in one place.** Building the segmentation
//! inside move 1 would mean a second copy of that assembly inside `ng::run`, where the first
//! copy — the one `call-from-alignments` and `generate-psps` share — lives a layer up, in the
//! command module. Reaching down for it from here would compile; what it would break is the
//! direction the dependencies run in, a pipeline stage calling the commands that drive it. Two
//! copies is exactly the drift psp mode's whole compatibility check exists to catch, so the
//! run object takes the segmentation as an argument, as direct mode's does.
//!
//! **No block is decoded by any of this.** Opening a psp reads its footer, its index and its
//! header and touches no block ([`PspReader::open`]), which is what lets a cohort of thousands
//! be opened and refused before any of it is spent.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ng::calling::inference::RunnableCallingLoopConfig;
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::psp::{ContigIdentity, Header, PspReader};
use crate::ng::read::input::read_groups::{
    NameOrigin, NameWithOrigin, ReadGroup, ReadGroups, SampleReadGroups,
};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::region_typing::GenomeRegions;
use crate::ng::types::{Bp, ReadGroupId};

use crate::ng::calling::allele_candidates::CandidateSelectionConfig;

use super::callers::{
    A_REFERENCE_WITH_NO_PATH, DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES, MergeParameters,
    refuse_a_catalog_built_on_another_reference, refuse_parameters_assembled_for_another_cohort,
    refuse_two_references_that_are_not_one,
};
use super::walker::WalkReference;
use super::{RunError, Segmentation};

/// Every psp of a cohort, open, with what its headers agree on.
///
/// **What this settles is what the cohort is**, before anything is built over it: which
/// individuals are in it and in what order, the ground they were all walked over, one run-wide
/// read-group numbering that the files' walk-local numberings merge into, and the widest span
/// any observation in any of them can have.
///
/// **It holds the open readers**, so the files are opened once for the run rather than once to
/// be checked and again to be read. What it does not hold is anything a block would have to be
/// decoded to learn.
#[derive(Debug)]
pub struct OpenPspCohort {
    /// One per sample, **in the order the paths were given** — which is the run's sample order
    /// and the order every per-sample list is indexed by.
    psps: Vec<PspReader>,
    /// Where each came from, kept beside the readers because a refusal names the file and a
    /// [`PspReader`] does not hand its path out.
    paths: Vec<PathBuf>,
    /// The ground every file agrees it was walked over. **Taken from the first file and then
    /// required of the rest** — the cohort refusal (spec §6.2).
    analysed_regions: GenomeRegions,
    /// Every read group of every sample, numbered run-wide. **This is also the remap** — see
    /// [`read_group_remap`](Self::read_group_remap).
    read_groups: ReadGroups,
    /// The widest reference span any observation in any of these files can have — the maximum
    /// over the files' own ceilings, which is what a cohort reader must size for.
    observation_reach_ceiling: Bp,
}

impl OpenPspCohort {
    /// Open every psp, read every header, and refuse a cohort that cannot be called as one.
    ///
    /// **The order of `psps` is the run's sample order** and reaches the VCF's sample columns,
    /// so it is the caller's to fix and this preserves it.
    ///
    /// # Errors
    ///
    /// [`RunError::NoPsps`] for an empty list; [`RunError::NotEnoughFileDescriptorsForPsps`]
    /// for a cohort this process may not hold open; [`RunError::PspNotRead`] naming the file
    /// for one that will not open — including a walk that was interrupted, which the cause says;
    /// [`RunError::SampleAppearsTwice`] for two files naming one individual;
    /// [`RunError::AnalysedRegionsDiffer`] for two walked over different ground; and
    /// [`RunError::PspReadGroupsCannotBeMerged`] for a table that cannot be renumbered.
    pub fn open(paths: &[PathBuf]) -> Result<Self, RunError> {
        // `current` is the soft limit — what this process may open now.
        let limit = rustix::process::getrlimit(rustix::process::Resource::Nofile).current;
        Self::open_within_a_descriptor_limit(paths, limit)
    }

    /// The same, with the descriptor limit passed in.
    ///
    /// **Split from the syscall so that the refusal's *place* has a test and not only its
    /// message.** The limit a machine reports is far above any cohort a test can build, so a
    /// check that read it inline could be moved below the open loop — or deleted — with
    /// nothing failing, and where it sits is the whole point of it (spec §7.1a).
    pub(crate) fn open_within_a_descriptor_limit(
        paths: &[PathBuf],
        limit: Option<u64>,
    ) -> Result<Self, RunError> {
        if paths.is_empty() {
            return Err(RunError::NoPsps);
        }
        // **Before the first `open(2)`, not after the last.** A cohort of thousands would
        // otherwise meet the operating system's own limit part-way through and report it
        // against whichever file happened to be next — an innocent path, and a message saying
        // nothing about the limit or how many samples the run has (spec §7.1a).
        refuse_if_more_descriptors_are_needed_than_allowed_for_psps(paths.len(), limit)?;

        let mut psps = Vec::with_capacity(paths.len());
        for path in paths {
            let psp = PspReader::open(path).map_err(|source| RunError::PspNotRead {
                path: path.clone(),
                source: Box::new(source),
            })?;
            psps.push(psp);
        }

        let headers: Vec<&Header> = psps.iter().map(PspReader::header).collect();
        refuse_a_sample_named_twice(&headers, paths)?;
        let analysed_regions = the_ground_every_file_agrees_on(&headers)?;
        let read_groups = merge_the_read_group_tables(&headers, paths)?;
        // **The maximum, not the first file's**: a cohort reader sizes for the widest
        // observation it can meet, and the files were written by separate invocations that may
        // have been given different caps.
        let observation_reach_ceiling = headers
            .iter()
            .map(|header| header.observation_reach_ceiling_bp)
            .max()
            .expect("the empty cohort was refused above");

        Ok(Self {
            psps,
            paths: paths.to_vec(),
            analysed_regions,
            read_groups,
            observation_reach_ceiling,
        })
    }

    /// How many samples this cohort holds.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.psps.len()
    }

    /// The individuals, in the run's sample order.
    pub fn sample_names(&self) -> impl Iterator<Item = &str> {
        self.psps.iter().map(|psp| psp.header().sample.as_str())
    }

    /// **The ground this run analyses** — read out of the files rather than asked for (spec
    /// §5.3), and the same in every one of them because a cohort that disagreed was refused.
    #[must_use]
    pub fn analysed_regions(&self) -> &GenomeRegions {
        &self.analysed_regions
    }

    /// Every read group of every sample, numbered run-wide.
    #[must_use]
    pub fn read_groups(&self) -> &ReadGroups {
        &self.read_groups
    }

    /// The widest reference span an observation in this cohort can have — what the merge sizes
    /// its cache for (`cohort_merge.md` §13).
    #[must_use]
    pub fn observation_reach_ceiling(&self) -> Bp {
        self.observation_reach_ceiling
    }

    /// The file each sample was read from, in the run's sample order.
    #[must_use]
    pub fn paths(&self) -> &[PathBuf] {
        &self.paths
    }

    /// One sample's map from the numbers its own walk gave its read groups to the numbers this
    /// run gives them, by the sample's position in the run's order.
    ///
    /// **The merged table is already the map, and keeping a second copy would be the hazard
    /// [`ReadGroups::of_merged_tables`] panics about** — a third per-sample vector that could
    /// disagree with the other two and put one sample's reads under another's calibration. It
    /// is the map because two rules meet: the run-wide identifiers of one sample's groups are
    /// pushed in the file's own header order, and a psp's header order *is* its walk-local
    /// numbering — entry `i` carries number `i`, which the format checks on both sides
    /// ([`ReadGroupIdentity`](crate::ng::psp::ReadGroupIdentity)).
    #[must_use]
    pub fn read_group_remap(&self, sample: usize) -> Option<&[ReadGroupId]> {
        self.read_groups
            .read_groups_per_sample()
            .get(sample)
            .map(|sample| sample.read_groups.as_slice())
    }
}

/// The reference a stored cohort is called against, in the two forms a run holds it.
///
/// **Two views of one genome, and both are needed.** The first is what the run opened for
/// reading bases — a record with an empty allele is padded from it. The second is the same
/// reference once its per-contig checksums are known, which is what every psp's contig table
/// is compared against; only the caller can supply it, because a `.fai`-only reference and one
/// whose FASTA has not been read are the same value and a run must report which it got rather
/// than infer it. [`AlignmentInputs`](super::AlignmentInputs) carries the same pair for the
/// same reasons.
pub struct StoredCohortInputs<'a> {
    /// What the run opened for reading bases.
    pub reference: &'a OpenReference,
    /// The same reference once its per-contig checksums are known.
    pub reference_with_checksums: &'a ReferenceInfo,
}

/// psp mode's calling stage (spec §5.3).
///
/// Holds every sample's psp open for the whole run — **357 kB an open file on a human
/// reference and 7 kB on tomato, measured** (spec §7.2), almost all of it the reference's
/// contig list — plus the shared read-only state a call needs. A cursor walking a file costs
/// 123 kB more, and this stage holds none: the sources the calling loop draws through are
/// built where the walk starts.
/// The observations come out of the files; everything above them is the body direct mode drives.
#[derive(Debug)]
pub struct PspVariantCaller {
    cohort: OpenPspCohort,
    segmentation: Arc<Segmentation>,
    parameters: RunParameters,
    /// **Checked and kept here, read at step E3.** These four are what the calling loop is
    /// handed once this caller drives it, and they are settled at `open` for the reason every
    /// other refusal is: a run whose reference cannot serve a padding base should be told
    /// before a block is decoded rather than at its first locus.
    ///
    /// **`expect` rather than `allow`, so the first real reader turns the line into a compile
    /// error** rather than leaving behind a suppression nobody removes — the rule this crate
    /// already follows where it defers a field.
    #[expect(dead_code, reason = "the calling loop reads it at plan step E3")]
    walk_reference: WalkReference,
    #[expect(dead_code, reason = "the calling loop reads it at plan step E3")]
    calling_loop_config: RunnableCallingLoopConfig,
    #[expect(dead_code, reason = "the calling loop reads it at plan step E3")]
    candidate_selection: CandidateSelectionConfig,
    #[expect(dead_code, reason = "the calling loop reads it at plan step E3")]
    merge_parameters: MergeParameters,
}

impl PspVariantCaller {
    /// Check an opened cohort against the run it is about to be called by, and keep it.
    ///
    /// **Every refusal here compares one file with the run**, where
    /// [`OpenPspCohort::open`]'s compare the files with each other.
    ///
    /// **`segmentation` must have been built over [`OpenPspCohort::analysed_regions`], and that
    /// does get checked** — not directly, but as a consequence: every psp's analysed regions
    /// were just forced equal to the cohort's, and `first_difference` compares that field, so a
    /// segmentation built over other ground is refused as
    /// [`RunError::SegmentationInputsDiffer`] naming the analysed regions. **It reads as though
    /// the file were at fault when it is the caller**, which is why the two moves are
    /// documented together at the top of this module. What is genuinely unchecked is narrower:
    /// nothing compares the segmentation's *segments* against the ground its own record claims,
    /// only the record.
    ///
    /// # Errors
    ///
    /// [`RunError::ParametersAreForAnotherCohort`] when the parameters were assembled for a
    /// different cohort; [`RunError::ReferenceHasNoBases`] or
    /// [`RunError::ReferenceIndexUnreadable`] for a reference a record cannot be padded from;
    /// [`RunError::ReferenceCheckedAgainstAnotherGenome`] when the two views of the reference
    /// are not one genome; [`RunError::PspAgainstAnotherReference`] for a file whose contig
    /// table is not this run's; and [`RunError::SegmentationInputsDiffer`] for one walked under
    /// a different catalog or different repeat-tract criteria.
    pub fn open(
        cohort: OpenPspCohort,
        reference: StoredCohortInputs<'_>,
        segmentation: Segmentation,
        parameters: RunParameters,
        calling_loop_config: RunnableCallingLoopConfig,
        candidate_selection: CandidateSelectionConfig,
        merge_parameters: MergeParameters,
    ) -> Result<Self, RunError> {
        let StoredCohortInputs {
            reference,
            reference_with_checksums,
        } = reference;
        // **A count, where spec §6.2 asks for a match by name.** It cannot be made here:
        // `RunParameters` is assembled per sample by position and carries no names at all, so
        // the by-name match belongs where the parameters *file* is read against this cohort's
        // sample list — `ParametersFile::to_run_parameters_for`, which the subcommand calls
        // (plan step F1). What is left for a run to catch is parameters assembled for one
        // cohort handed to a caller opened over another, which is what this counts.
        refuse_parameters_assembled_for_another_cohort(&parameters, cohort.read_groups())?;
        refuse_two_references_that_are_not_one(reference.info(), reference_with_checksums)?;
        // **The reference is opened for reading bases here**, with direct mode's own refusal:
        // a record with an empty allele is padded from it, so a `.fai` alone cannot call.
        let walk_reference = WalkReference::of(reference)?;

        // **Direct mode's own catalog check, which psp mode needs for the same reason and had
        // no counterpart for.** The catalog's coordinates are where the repeat tracts are and
        // every segment this run loops over is drawn from it, so a catalog built on another
        // build of the assembly puts every tract at the wrong position — silently, and
        // genome-wide.
        let reference_path = reference_with_checksums
            .fasta_path
            .as_deref()
            .or(reference.info().fasta_path.as_deref())
            .unwrap_or_else(|| Path::new(A_REFERENCE_WITH_NO_PATH));
        refuse_a_catalog_built_on_another_reference(
            &segmentation,
            reference_with_checksums,
            reference_path,
        )?;

        for psp in &cohort.psps {
            let header = psp.header();
            refuse_a_file_against_another_reference(header, reference_with_checksums)?;
            if let Some(field) = header
                .segmentation_inputs
                .first_difference(segmentation.inputs())
            {
                return Err(RunError::SegmentationInputsDiffer {
                    sample: header.sample.clone(),
                    field,
                });
            }
        }

        Ok(Self {
            cohort,
            segmentation: Arc::new(segmentation),
            walk_reference,
            parameters,
            calling_loop_config,
            candidate_selection,
            merge_parameters,
        })
    }

    /// How many samples this run calls.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.cohort.sample_count()
    }

    /// The sample names, in the run's sample order.
    pub fn sample_names(&self) -> impl Iterator<Item = &str> {
        self.cohort.sample_names()
    }

    /// Every read group of every sample, numbered run-wide.
    #[must_use]
    pub fn read_groups(&self) -> &ReadGroups {
        self.cohort.read_groups()
    }

    /// The ground this run analyses, and the record of what it was computed from.
    #[must_use]
    pub fn segmentation(&self) -> &Segmentation {
        &self.segmentation
    }

    /// The cohort's open files and what their headers agreed on.
    #[must_use]
    pub fn cohort(&self) -> &OpenPspCohort {
        &self.cohort
    }

    /// The numbers every call is scored against.
    #[must_use]
    pub fn parameters(&self) -> &RunParameters {
        &self.parameters
    }
}

/// Refuse a cohort holding one individual twice (spec §6.2).
///
/// **A linear scan and not a map**, for `group_by_sample`'s reason: it keeps the refusal
/// naming the *first* pair in the order the files were given, which is the order a person
/// typed them in, and a cohort is small enough that the scan costs nothing beside opening the
/// files.
fn refuse_a_sample_named_twice(headers: &[&Header], paths: &[PathBuf]) -> Result<(), RunError> {
    for (later, header) in headers.iter().enumerate() {
        if let Some(earlier) = headers[..later]
            .iter()
            .position(|seen| seen.sample == header.sample)
        {
            return Err(RunError::SampleAppearsTwice {
                sample: header.sample.clone(),
                first: paths[earlier].clone(),
                second: paths[later].clone(),
            });
        }
    }
    Ok(())
}

/// The ground the cohort was walked over, or the refusal naming the first two files that
/// disagree (spec §6.2).
///
/// **The first file's ground is the run's**, and every other file is required to match it. That
/// makes the refusal name a pair rather than describe a set, which is what a person can act on
/// — and it is why [`RunError::AnalysedRegionsDiffer`]'s two fields are named for their roles
/// rather than symmetrically.
fn the_ground_every_file_agrees_on(headers: &[&Header]) -> Result<GenomeRegions, RunError> {
    let first = headers.first().expect("the empty cohort was refused above");
    let ground = &first.segmentation_inputs.analysed_regions;
    for header in &headers[1..] {
        if header.segmentation_inputs.analysed_regions != *ground {
            return Err(RunError::AnalysedRegionsDiffer {
                left: first.sample.clone(),
                right: header.sample.clone(),
            });
        }
    }
    Ok(ground.clone())
}

/// Refuse a file written against a different reference from the run's (spec §6.2).
///
/// **Two things are compared and they fail differently.** The contig *table* — names, lengths
/// and order — is what every record's coordinates are written against, so a file whose table
/// is not the run's has every observation on the wrong chromosome. The per-contig checksums are
/// compared only where both sides carry one: a psp written from a `.fai`-only reference has
/// none, and an absent digest is a check that could not be made rather than one that passed.
fn refuse_a_file_against_another_reference(
    header: &Header,
    reference: &ReferenceInfo,
) -> Result<(), RunError> {
    let refuse = |difference: String| RunError::PspAgainstAnotherReference {
        sample: header.sample.clone(),
        difference,
    };
    // **The whole-assembly digest first, because it is the one comparison that is exact.** It
    // is what the header's `reference` field exists for, and where both sides carry one it
    // settles the question before any contig is walked; a `.fai`-only read on either side
    // carries none, and an absent digest is a check that could not be made.
    if let (Some(stored), Some(run)) = (header.reference.md5, reference.md5)
        && stored != run
    {
        return Err(refuse(format!(
            "it was walked against the assembly whose checksum is {} and this run's reference \
             is {}",
            crate::pop_var_caller::common::format_md5_hex(stored),
            crate::pop_var_caller::common::format_md5_hex(run),
        )));
    }
    if header.contigs.len() != reference.contigs.len() {
        return Err(refuse(format!(
            "the psp describes {} contigs and this run's reference {}",
            header.contigs.len(),
            reference.contigs.len(),
        )));
    }
    for (at, (stored, run)) in header.contigs.iter().zip(&reference.contigs).enumerate() {
        let ContigIdentity {
            name,
            length,
            md5: stored_md5,
        } = stored;
        if *name != run.name {
            return Err(refuse(format!(
                "contig {at} is '{name}' in the psp and '{}' in this run's reference",
                run.name,
            )));
        }
        if *length != run.length {
            return Err(refuse(format!(
                "contig '{name}' is {length} bases in the psp and {} in this run's reference",
                run.length,
            )));
        }
        if let (Some(stored_md5), Some(run_md5)) = (stored_md5, run.md5)
            && *stored_md5 != run_md5
        {
            return Err(refuse(format!(
                "contig '{name}' has checksum {} in the psp and {} in this run's reference",
                crate::pop_var_caller::common::format_md5_hex(*stored_md5),
                crate::pop_var_caller::common::format_md5_hex(run_md5),
            )));
        }
    }
    Ok(())
}

/// Build the run-wide read-group table out of the cohort's walk-local ones, and the map from
/// each file's numbering into it (spec §6.2).
///
/// **Merged, never compared.** Every psp numbers its own read groups from zero, so every
/// sample's first group is identifier 0 in its own file; that collision is the normal case and
/// the whole reason the table is recorded. The run's numbering is *first file first, header
/// order within a file*, which is `build_read_groups`' rule for alignment files — so a cohort
/// walked and a cohort stored produce the same numbering from the same sample order, which is
/// what spec §12.3's mode-equivalence oracle needs — it compares the two modes' VCF bytes, and
/// a read group's number reaches the file.
///
/// **What a psp cannot carry, and what is done about it.** A `ReadGroup` records where its
/// library name came from and which experiment it belongs to; the header records the `@RG ID`,
/// the library and the walk-local number, and no more (spec §6.1). So the library is taken as
/// **declared** — it is the file's own tag, verbatim, which is exactly what `Declared` means at
/// this layer — and the experiment falls back to the library, which is the same fallback
/// `read_groups` applies to a file that declares no experiment. The platform is dropped: it is
/// carried for reports and nothing keys on it.
fn merge_the_read_group_tables(
    headers: &[&Header],
    paths: &[PathBuf],
) -> Result<ReadGroups, RunError> {
    let mut groups: Vec<ReadGroup> = Vec::new();
    let mut per_sample: Vec<SampleReadGroups> = Vec::new();

    for (header, path) in headers.iter().zip(paths) {
        refuse_a_table_that_cannot_be_renumbered(header)?;
        let file: Arc<Path> = Arc::from(path.as_path());
        let mut mine = Vec::with_capacity(header.read_groups.len());
        for stored in &header.read_groups {
            let id =
                ReadGroupId(u32::try_from(groups.len()).expect("a read-group table fits in u32"));
            groups.push(ReadGroup {
                file: Arc::clone(&file),
                id: stored.id.clone().into_boxed_str(),
                sample: header.sample.clone().into_boxed_str(),
                // **`Synthesized`, and it is the weaker of the two claims on purpose.** A psp
                // records the library the walk *resolved* — `@RG LB`, or the name the walk
                // invented when the file declared none — and not which of the two it was (spec
                // §6.1). The field exists to tell "the file's" from "ours", and `run_report`
                // states the rule: *a synthesized name reported as a declared one is a claim
                // about the run that nobody made*. So this says the thing that cannot be
                // false; recording the origin in the header would let it say the true one.
                library: NameWithOrigin {
                    value: stored.library.clone().into_boxed_str(),
                    origin: NameOrigin::Synthesized,
                },
                // **Always synthesized, which is direct mode's own rule**: the origin
                // describes *this* name, and the experiment is the library copied because
                // nothing reads an experiment tag yet — we chose to call it that, the file
                // did not (`read_groups.rs`, `into_read_group`).
                experiment: NameWithOrigin {
                    value: stored.library.clone().into_boxed_str(),
                    origin: NameOrigin::Synthesized,
                },
                platform: None,
            });
            mine.push(id);
        }
        per_sample.push(SampleReadGroups {
            sample: header.sample.clone().into_boxed_str(),
            read_groups: mine,
        });
    }

    Ok(ReadGroups::of_merged_tables(groups, per_sample))
}

/// Refuse a psp's read-group table that cannot be renumbered (spec §6.2).
///
/// **One way it cannot be: a sample with no table at all**, because then nothing says which
/// group a record's reads came from and no renumbering can be invented.
///
/// **⚑ Spec §6.2 names a second — two entries of one sample sharing an `@RG ID` — and it is
/// deliberately not refused here.** The reason the spec gives is that such a table "cannot be
/// renumbered without guessing", and that is not this format's situation: a psp's identity is
/// the walk-local *number*, which is the entry's own position, checked on both sides
/// ([`ReadGroupIdentity`](crate::ng::psp::ReadGroupIdentity)) — nothing in the merge reads the
/// id at all. And the format's own validator declares the case legal in as many words: a **psp
/// holds one sample, not one alignment file**, and a sample sequenced across several files may
/// carry two entries with one `@RG ID` and different libraries (`psp/header.rs`). Direct mode
/// calls that cohort without complaint, so refusing it here would break spec §1.1's goal 1 for
/// every multi-lane sample whose lanes reuse an id — common in real archives — and cost a
/// re-walk nobody could avoid. **Raised for the owner at Checkpoint E: §6.2's clause should say
/// what it means, which is the empty table.**
fn refuse_a_table_that_cannot_be_renumbered(header: &Header) -> Result<(), RunError> {
    if header.read_groups.is_empty() {
        return Err(RunError::PspReadGroupsCannotBeMerged {
            sample: header.sample.clone(),
            problem: "it names no read groups at all, so nothing says which group a record's \
                      reads came from"
                .to_string(),
        });
    }
    Ok(())
}

/// Refuse a cohort this process may not hold open (spec §7.1a).
///
/// **The arithmetic is written out because the message asks an operator to act on it**: an open
/// psp costs one descriptor, held for the whole run (`PspReader` keeps one `File` and opens
/// nothing else), against an alignment file's two — its own reader, and the reference accessor
/// its cursor mints.
///
/// `None` means the platform reports no limit at all, and then there is nothing to refuse
/// against.
fn refuse_if_more_descriptors_are_needed_than_allowed_for_psps(
    samples: usize,
    limit: Option<u64>,
) -> Result<(), RunError> {
    let needed = samples as u64 + DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES;
    let Some(limit) = limit else {
        return Ok(());
    };
    if needed > limit {
        return Err(RunError::NotEnoughFileDescriptorsForPsps {
            samples,
            allowance: DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES,
            needed,
            limit,
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::psp::writer::tests_support::a_header;
    use crate::ng::psp::{PspWriter, ReadGroupIdentity};
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, fixture_reference_from_its_index,
    };
    use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion};
    use crate::ng::run::test_fixtures::{
        build_segmentation_under, catalog_header, catalog_header_built_on,
    };
    use crate::ng::types::{ContigId, GenomeRegion, Ploidy, Position};
    use crate::regions::{ContigBounds, Region};
    use tempfile::TempDir;

    /// The two contigs every fixture here shares — the same table
    /// [`fixture_reference_from_its_index`] builds, so a psp written over these and a run
    /// opened over that reference agree by construction and a test that wants them to disagree
    /// has to say so.
    fn fixture_contigs() -> Vec<ContigIdentity> {
        crate::ng::read::input::test_fixtures::FIXTURE_CONTIGS
            .iter()
            .map(|(name, length)| ContigIdentity {
                name: (*name).to_string(),
                length: *length as u64,
                md5: None,
            })
            .collect()
    }

    fn bounds() -> Vec<ContigBounds<'static>> {
        crate::ng::read::input::test_fixtures::FIXTURE_CONTIGS
            .iter()
            .map(|(name, length)| ContigBounds {
                name,
                length: *length as u32,
            })
            .collect()
    }

    /// One span on each contig, or a narrower one when `narrow` — the two grounds a
    /// disagreeing cohort is built from.
    fn ground(narrow: bool) -> GenomeRegions {
        let end = if narrow { 40 } else { 80 };
        GenomeRegions::from_normalized_spans(
            vec![
                Region {
                    chrom_id: 0,
                    start: 10,
                    end,
                },
                Region {
                    chrom_id: 1,
                    start: 10,
                    end: 150,
                },
            ],
            &bounds(),
        )
        .expect("the fixture's spans are normalized")
    }

    /// A run over `ground`, cut into one segment per contig, under the fixture catalog.
    fn a_segmentation(over: GenomeRegions) -> Arc<Segmentation> {
        a_segmentation_under(over, catalog_header())
    }

    /// The same, with the catalog named — for the tests that are about the catalog's identity
    /// or that need it to agree with a reference whose digests are real.
    fn a_segmentation_under(
        over: GenomeRegions,
        catalog: crate::ng::repeat_catalog::RepeatCatalogHeader,
    ) -> Arc<Segmentation> {
        build_segmentation_under(
            vec![
                TypedRegion {
                    region: GenomeRegion {
                        contig: ContigId(0),
                        start: Position(10),
                        end: Position(40),
                    },
                    kind: RegionKind::Generic,
                },
                TypedRegion {
                    region: GenomeRegion {
                        contig: ContigId(1),
                        start: Position(10),
                        end: Position(150),
                    },
                    kind: RegionKind::Generic,
                },
            ],
            over,
            catalog,
        )
    }

    /// One sample's psp: the fixture contig table, the run's own segmentation inputs, and no
    /// records at all — **E1 decodes no block, so an empty file exercises everything it does.**
    fn a_psp(
        dir: &TempDir,
        sample: &str,
        segmentation: &Segmentation,
        edit: impl FnOnce(&mut Header),
    ) -> PathBuf {
        let mut header = a_header(1_000);
        header.sample = sample.to_string();
        header.contigs = fixture_contigs();
        header.segmentation_inputs = segmentation.inputs().clone();
        edit(&mut header);
        let path = dir.path().join(format!("{sample}.psp"));
        let writer = PspWriter::create(&path, header).expect("the header writes");
        let _ = writer.finish(b"").expect("the file seals");
        path
    }

    /// The ordinary two-sample cohort, agreeing about everything.
    fn a_cohort() -> (TempDir, Arc<Segmentation>, Vec<PathBuf>) {
        a_cohort_whose_second_file(|_| {})
    }

    /// The same cohort with `edit` applied to the **second** file's header.
    ///
    /// **Second and not first, and that is the point of the helper.** A refusal provoked on a
    /// one-file cohort, or on the first file of two, cannot tell a loop over every psp from one
    /// that checks the first and stops — and a mistake confined to a later file (one accession
    /// re-walked against another build) is the ordinary shape of the fault at 63 samples or at
    /// 3,000. Measured: with every per-file refusal pinned on a one-file cohort, a mutant that
    /// checked only `psps[0]` passed all 490 tests of `ng::run`.
    fn a_cohort_whose_second_file(
        edit: impl FnOnce(&mut Header),
    ) -> (TempDir, Arc<Segmentation>, Vec<PathBuf>) {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let segmentation = a_segmentation(ground(false));
        let paths = vec![
            a_psp(&dir, "alpha", &segmentation, |_| {}),
            a_psp(&dir, "beta", &segmentation, edit),
        ];
        (dir, segmentation, paths)
    }

    /// Parameters assembled for exactly the cohort a table describes.
    fn parameters_for(read_groups: &ReadGroups) -> RunParameters {
        RunParameters::of_defaults(
            read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        )
    }

    fn open_the_caller(
        cohort: OpenPspCohort,
        reference: &OpenReference,
        segmentation: Segmentation,
        parameters: RunParameters,
    ) -> Result<PspVariantCaller, RunError> {
        PspVariantCaller::open(
            cohort,
            StoredCohortInputs {
                reference,
                reference_with_checksums: reference.info(),
            },
            segmentation,
            parameters,
            RunnableCallingLoopConfig::default(),
            CandidateSelectionConfig::default(),
            MergeParameters::DEFAULT,
        )
    }

    // -----------------------------------------------------------------
    // What the headers settle (spec §6.2's cohort checks)
    // -----------------------------------------------------------------

    /// **The run's ground comes out of the files** (spec §5.3), and the run-wide read-group
    /// numbering out of their tables.
    #[test]
    fn a_cohort_of_agreeing_psps_opens_and_says_what_it_is() {
        let (_dir, segmentation, paths) = a_cohort();
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree");

        assert_eq!(cohort.sample_count(), 2);
        assert_eq!(
            cohort.sample_names().collect::<Vec<_>>(),
            vec!["alpha", "beta"],
            "the run's sample order is the order the paths were given",
        );
        assert_eq!(
            cohort.analysed_regions(),
            &segmentation.inputs().analysed_regions,
            "the ground is read out of the headers, not asked for",
        );
        assert_eq!(cohort.paths(), paths.as_slice());
    }

    /// **Every psp numbers its read groups from zero, so the run must renumber them.** Both
    /// files here carry the same two walk-local numbers; a run that did not renumber would put
    /// beta's reads under alpha's calibration.
    #[test]
    fn the_read_group_tables_merge_into_one_run_wide_numbering() {
        let (_dir, _segmentation, paths) = a_cohort();
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree");

        assert_eq!(cohort.read_groups().len(), 4, "two read groups a sample");
        assert_eq!(
            cohort.read_group_remap(0),
            Some([ReadGroupId(0), ReadGroupId(1)].as_slice()),
            "the first file's groups keep their numbers",
        );
        assert_eq!(
            cohort.read_group_remap(1),
            Some([ReadGroupId(2), ReadGroupId(3)].as_slice()),
            "the second file's start where the first's ended, which is what renumbering means",
        );

        // **Which entry got which number, not only how many.** Both of a sample's groups carry
        // the same sample name, so a merge that walked a file's table backwards would still
        // hand back `[0, 1]` — and every observation would then land on the wrong read group's
        // calibration with nothing failing.
        let ids: Vec<&str> = cohort
            .read_groups()
            .iter()
            .map(|(_, group)| &*group.id)
            .collect();
        assert_eq!(
            ids,
            vec!["SRR7279481", "SRR7279481.L2", "SRR7279481", "SRR7279481.L2"],
            "first file first, header order within a file — `build_read_groups`' own rule",
        );
        assert_eq!(
            &*cohort
                .read_groups()
                .get(cohort.read_group_remap(1).expect("the second sample")[1])
                .library
                .value,
            "tomato-pe-2",
            "the remap's entry `i` is the run number of that file's walk-local group `i`",
        );

        let per_sample = cohort.read_groups().read_groups_per_sample();
        assert_eq!(per_sample.len(), 2);
        assert_eq!(&*per_sample[0].sample, "alpha");
        assert_eq!(&*per_sample[1].sample, "beta");
        for (id, group) in cohort.read_groups().iter() {
            let sample = if id.get() < 2 { "alpha" } else { "beta" };
            assert_eq!(&*group.sample, sample, "read group {id:?}");
        }
    }

    /// **The maximum, not the first file's**: the files are written by separate invocations,
    /// which may have been given different caps, and a cohort reader must size for the widest.
    #[test]
    fn the_reach_ceiling_is_the_widest_any_file_declares() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let segmentation = a_segmentation(ground(false));
        let paths = vec![
            a_psp(&dir, "alpha", &segmentation, |header| {
                header.observation_reach_ceiling_bp = Bp(300);
            }),
            a_psp(&dir, "beta", &segmentation, |header| {
                header.observation_reach_ceiling_bp = Bp(1_200);
            }),
        ];
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree");
        assert_eq!(cohort.observation_reach_ceiling(), Bp(1_200));
    }

    #[test]
    fn a_run_given_no_psps_is_refused() {
        assert!(matches!(OpenPspCohort::open(&[]), Err(RunError::NoPsps)));
    }

    /// **An interrupted walk is refused as interrupted, not read as a short sample** (spec §9).
    /// The file has no name to be reported under — the sample lives in the header, which comes
    /// after the footer that is missing — so the refusal names the path.
    #[test]
    fn a_psp_whose_walk_was_interrupted_is_refused_naming_the_file() {
        let (_dir, _segmentation, paths) = a_cohort();
        let whole = std::fs::read(&paths[1]).expect("the file reads");
        std::fs::write(&paths[1], &whole[..whole.len() - 8]).expect("the truncated copy writes");

        let refused = OpenPspCohort::open(&paths).expect_err("a file with no footer is refused");
        let RunError::PspNotRead { path, source } = refused else {
            panic!("a file that will not read must name itself: {refused:?}");
        };
        assert_eq!(path, paths[1]);
        assert!(
            source.to_string().contains("the writer did not finish"),
            "the cause tells an interrupted walk from a damaged file: {source}",
        );
    }

    /// **A cohort holds each individual once** — two files for one sample would weight every
    /// allele frequency by it.
    #[test]
    fn two_psps_naming_one_sample_are_refused_naming_both_files() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let segmentation = a_segmentation(ground(false));
        let first = a_psp(&dir, "alpha", &segmentation, |_| {});
        let again = dir.path().join("alpha-again.psp");
        std::fs::copy(&first, &again).expect("the second copy writes");
        let paths = vec![first.clone(), again.clone()];

        let refused = OpenPspCohort::open(&paths).expect_err("one individual, two files");
        let RunError::SampleAppearsTwice {
            sample,
            first: named_first,
            second,
        } = refused
        else {
            panic!("the duplicate must be named: {refused:?}");
        };
        assert_eq!(sample, "alpha");
        assert_eq!(named_first, first);
        assert_eq!(second, again);
    }

    /// **The one cohort mismatch that would give a wrong answer rather than a missing one.**
    /// A sample has no records over ground it never walked, and that absence reads exactly like
    /// *no variant here*.
    #[test]
    fn two_psps_walked_over_different_ground_are_refused_naming_both_samples() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let wide = a_segmentation(ground(false));
        let narrow = a_segmentation(ground(true));
        let paths = vec![
            a_psp(&dir, "alpha", &wide, |_| {}),
            a_psp(&dir, "beta", &narrow, |_| {}),
        ];

        let refused = OpenPspCohort::open(&paths).expect_err("two grounds, one cohort");
        let RunError::AnalysedRegionsDiffer { left, right } = refused else {
            panic!("both samples must be named: {refused:?}");
        };
        assert_eq!((left.as_str(), right.as_str()), ("alpha", "beta"));
    }

    /// A table with nothing in it cannot say which group a record's reads came from.
    ///
    /// **Provoked at the function rather than through a file, because no file this store wrote
    /// can carry such a table**: `PspWriter::create` refuses the header outright — *"is empty;
    /// a sample's reads come from at least one read group"*. The check stands for a file
    /// written by something that is not this store, and this is how a check whose input this
    /// build cannot produce is still pinned.
    #[test]
    fn a_table_naming_no_read_groups_cannot_be_renumbered() {
        let mut header = a_header(1_000);
        header.sample = "alpha".to_string();
        header.read_groups.clear();

        let refused = refuse_a_table_that_cannot_be_renumbered(&header)
            .expect_err("an empty table cannot be merged");
        let RunError::PspReadGroupsCannotBeMerged { sample, problem } = refused else {
            panic!("the sample must be named: {refused:?}");
        };
        assert_eq!(sample, "alpha");
        assert!(problem.contains("no read groups"), "{problem}");
    }

    /// **A sample sequenced across two files whose lanes reuse an `@RG ID` calls.**
    ///
    /// A psp holds one *sample*, not one alignment file, and SAM makes an id unique only within
    /// one file — so this table is what the format's own validator declares legal and what
    /// direct mode calls without complaint. Nothing in the merge reads the id: identity is the
    /// walk-local number, which is the entry's position. Spec §6.2 names this as a refusal and
    /// this is the fixture that says why it must not be one — refusing it would cost a re-walk
    /// for a cohort direct mode calls.
    #[test]
    fn a_sample_walked_across_two_files_sharing_a_read_group_id_calls() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let segmentation = a_segmentation(ground(false));
        let paths = vec![a_psp(&dir, "alpha", &segmentation, |header| {
            header.read_groups = vec![
                ReadGroupIdentity {
                    id: "L1".to_string(),
                    library: "lib-a".to_string(),
                    walk_local_id: ReadGroupId(0),
                },
                ReadGroupIdentity {
                    id: "L1".to_string(),
                    library: "lib-b".to_string(),
                    walk_local_id: ReadGroupId(1),
                },
            ];
        })];

        let cohort = OpenPspCohort::open(&paths).expect("two lanes of one sample are a cohort");
        assert_eq!(cohort.read_groups().len(), 2);
        let libraries: Vec<&str> = cohort
            .read_groups()
            .iter()
            .map(|(_, group)| &*group.library.value)
            .collect();
        assert_eq!(
            libraries,
            vec!["lib-a", "lib-b"],
            "the two lanes stay apart, which is what the walk-local number keeps",
        );
    }

    /// **One sample is a cohort**, and it is the end of the range this caller commits to.
    #[test]
    fn a_cohort_of_one_sample_opens_and_calls() {
        let dir = tempfile::tempdir().expect("a temporary directory");
        let segmentation = a_segmentation(ground(false));
        let paths = vec![a_psp(&dir, "alpha", &segmentation, |_| {})];
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let cohort = OpenPspCohort::open(&paths).expect("one file is a cohort");
        assert_eq!(cohort.sample_count(), 1);
        let parameters = parameters_for(cohort.read_groups());

        let caller = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(a_segmentation(ground(false))).expect("one handle"),
            parameters,
        )
        .expect("one sample matches the run");
        assert_eq!(caller.sample_names().collect::<Vec<_>>(), vec!["alpha"]);
    }

    /// **The refusal spec §7.1a asks for by name**, and it has to fire before the first file is
    /// opened: a cohort of thousands otherwise meets the operating system's limit part-way
    /// through and blames whichever psp happened to be next.
    ///
    /// **Split from the syscall so the message has a reader.** The limit a machine reports is
    /// far above any cohort a test can build.
    #[test]
    fn a_cohort_this_process_could_not_hold_open_is_refused_naming_the_limit() {
        let refused = refuse_if_more_descriptors_are_needed_than_allowed_for_psps(3_000, Some(256))
            .expect_err("3,000 psps do not fit in 256 descriptors");
        let RunError::NotEnoughFileDescriptorsForPsps {
            samples,
            allowance,
            needed,
            limit,
        } = refused
        else {
            panic!("the limit and the count must be named: {refused:?}");
        };
        assert_eq!((samples, limit), (3_000, 256));
        assert_eq!(
            needed,
            3_000 + allowance,
            "one descriptor a psp, plus what the run needs besides them",
        );
        assert!(
            refuse_if_more_descriptors_are_needed_than_allowed_for_psps(3_000, None).is_ok(),
            "a platform that reports no limit has nothing to refuse against",
        );
        assert!(
            refuse_if_more_descriptors_are_needed_than_allowed_for_psps(64, Some(4_096)).is_ok(),
            "an ordinary cohort is not refused",
        );
    }

    /// **The refusal fires before the first file is opened**, which is the whole of what it
    /// buys — the paths here do not exist, so a run that opened first would come back
    /// [`RunError::PspNotRead`] naming one of them.
    #[test]
    fn the_descriptor_refusal_comes_before_any_psp_is_opened() {
        let paths: Vec<PathBuf> = (0..8)
            .map(|at| PathBuf::from(format!("/nowhere/{at}.psp")))
            .collect();
        let refused = OpenPspCohort::open_within_a_descriptor_limit(&paths, Some(4))
            .expect_err("eight psps do not fit in four descriptors");
        assert!(
            matches!(refused, RunError::NotEnoughFileDescriptorsForPsps { .. }),
            "the descriptors are refused before the paths are touched: {refused:?}",
        );
    }

    /// **The two views of the reference must be one genome**, because every comparison
    /// downstream walks them in step — a run checked against one reference and opened against
    /// another would pair a psp's chromosome with something else's.
    #[test]
    fn a_run_checked_against_another_genome_is_refused() {
        let (_dir, _segmentation, paths) = a_cohort();
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let mut elsewhere = reference.info().clone();
        elsewhere.contigs[1].name = "chrX".to_string();
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree");
        let parameters = parameters_for(cohort.read_groups());

        let refused = PspVariantCaller::open(
            cohort,
            StoredCohortInputs {
                reference: &reference,
                reference_with_checksums: &elsewhere,
            },
            Arc::try_unwrap(a_segmentation(ground(false))).expect("one handle"),
            parameters,
            RunnableCallingLoopConfig::default(),
            CandidateSelectionConfig::default(),
            MergeParameters::DEFAULT,
        )
        .expect_err("the two references are not one genome");
        let RunError::ReferenceCheckedAgainstAnotherGenome { difference } = refused else {
            panic!("the difference must be named: {refused:?}");
        };
        assert!(difference.contains("chrX"), "{difference}");
    }

    // -----------------------------------------------------------------
    // What the run checks each file against (spec §6.2's file-against-run checks)
    // -----------------------------------------------------------------

    #[test]
    fn a_checked_cohort_opens_as_a_caller() {
        let (_dir, segmentation, paths) = a_cohort();
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree");
        let parameters = parameters_for(cohort.read_groups());
        let over = a_segmentation(ground(false));

        let caller = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(over).expect("one handle"),
            parameters,
        )
        .expect("the cohort matches the run");

        assert_eq!(caller.sample_count(), 2);
        assert_eq!(
            caller.sample_names().collect::<Vec<_>>(),
            vec!["alpha", "beta"]
        );
        assert_eq!(
            caller.segmentation().inputs(),
            segmentation.inputs(),
            "the run loops over the segmentation it was handed",
        );
    }

    /// **The independence the calling loop rests on is what this protects.** Under a different
    /// catalog a stored repeat-tract observation can straddle a calling segment's edge, and
    /// *no observation crosses a segment's edge* stops being true.
    #[test]
    fn a_psp_walked_under_another_catalog_is_refused_naming_the_field() {
        let (_dir, _segmentation, paths) = a_cohort_whose_second_file(|header| {
            header.segmentation_inputs.catalog.reference_md5 = [9; 16];
        });
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree about the ground");
        let parameters = parameters_for(cohort.read_groups());

        let refused = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(a_segmentation(ground(false))).expect("one handle"),
            parameters,
        )
        .expect_err("the catalog differs");
        let RunError::SegmentationInputsDiffer { sample, field } = refused else {
            panic!("the field must be named: {refused:?}");
        };
        assert_eq!(
            (sample.as_str(), field),
            ("beta", "repeat catalog"),
            "the check reaches the second file, not only the first",
        );
    }

    /// **Every record is written against the contig table its header carries**, in that order,
    /// so a file whose table is not the run's puts every observation on the wrong chromosome.
    #[test]
    fn a_psp_written_against_another_contig_table_is_refused_naming_the_contig() {
        for (what, edit, expected) in [
            (
                "a contig count",
                // **Appended rather than removed**: the header's own consistency check refuses
                // analysed regions that index a contig the file does not declare, so a shorter
                // table cannot be written beside this ground.
                Box::new(|header: &mut Header| {
                    header.contigs.push(ContigIdentity {
                        name: "chr3".to_string(),
                        length: 300,
                        md5: None,
                    });
                }) as Box<dyn FnOnce(&mut Header)>,
                "describes 3 contigs",
            ),
            (
                "a contig name",
                Box::new(|header: &mut Header| {
                    header.contigs[1].name = "chrX".to_string();
                }),
                "contig 1 is 'chrX'",
            ),
            (
                "a contig length",
                Box::new(|header: &mut Header| {
                    header.contigs[1].length = 999;
                }),
                "is 999 bases in the psp",
            ),
        ] {
            let (_dir, _segmentation, paths) = a_cohort_whose_second_file(edit);
            let (_reference_dir, reference) = fixture_reference_from_its_index();
            let cohort = OpenPspCohort::open(&paths).expect("the two files agree about the ground");
            let parameters = parameters_for(cohort.read_groups());

            let refused = open_the_caller(
                cohort,
                &reference,
                Arc::try_unwrap(a_segmentation(ground(false))).expect("one handle"),
                parameters,
            )
            .expect_err("the contig table differs");
            let RunError::PspAgainstAnotherReference { sample, difference } = refused else {
                panic!("{what}: the difference must be named: {refused:?}");
            };
            assert_eq!(sample, "beta", "{what}: the check reaches the second file");
            assert!(difference.contains(expected), "{what}: {difference}");
        }
    }

    /// **An absent digest is a check that could not be made, not one that passed**, so the
    /// comparison happens only where both sides carry one — which needs the reference read
    /// from its FASTA rather than from its index.
    #[test]
    fn a_psp_written_against_other_bases_is_refused_naming_the_checksum() {
        let (_reference_dir, reference) = fixture_reference(true);
        let real = reference.info().contigs[1]
            .md5
            .expect("the fasta arm carries digests");
        // The catalog has to claim this reference, or the catalog check refuses first.
        let catalog = catalog_header_built_on(reference.info().md5.expect("a whole-genome digest"));
        let dir = tempfile::tempdir().expect("a temporary directory");
        let segmentation = a_segmentation_under(ground(false), catalog.clone());
        let paths = vec![
            a_psp(&dir, "alpha", &segmentation, |header| {
                header.reference.md5 = reference.info().md5;
                header.contigs[0].md5 = reference.info().contigs[0].md5;
                header.contigs[1].md5 = reference.info().contigs[1].md5;
            }),
            // The assembly digest agrees, so the per-contig comparison is what has to catch
            // this — which is the case that arises when two builds share a whole-genome
            // digest's worth of ground and differ inside one chromosome.
            a_psp(&dir, "beta", &segmentation, |header| {
                header.reference.md5 = reference.info().md5;
                header.contigs[0].md5 = reference.info().contigs[0].md5;
                header.contigs[1].md5 = Some([0xab; 16]);
            }),
        ];
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree about the ground");
        let parameters = parameters_for(cohort.read_groups());

        let refused = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(a_segmentation_under(ground(false), catalog)).expect("one handle"),
            parameters,
        )
        .expect_err("the bases differ");
        let RunError::PspAgainstAnotherReference { sample, difference } = refused else {
            panic!("the checksum must be named: {refused:?}");
        };
        assert_eq!(sample, "beta", "the check reaches the second file");
        assert!(
            difference.contains("has checksum")
                && difference.contains(&crate::pop_var_caller::common::format_md5_hex(real)),
            "the refusal names both checksums: {difference}",
        );
    }

    /// **The header's `reference` field is what this check exists for**, and it is the one
    /// comparison that is exact: a whole-assembly digest either matches or it does not, where
    /// contig names and lengths agree across builds of one assembly.
    #[test]
    fn a_psp_walked_against_another_assembly_is_refused_naming_both_digests() {
        let (_reference_dir, reference) = fixture_reference(true);
        let catalog = catalog_header_built_on(reference.info().md5.expect("a whole-genome digest"));
        let dir = tempfile::tempdir().expect("a temporary directory");
        let segmentation = a_segmentation_under(ground(false), catalog.clone());
        let paths = vec![
            a_psp(&dir, "alpha", &segmentation, |header| {
                header.reference.md5 = reference.info().md5;
            }),
            a_psp(&dir, "beta", &segmentation, |header| {
                header.reference.md5 = Some([0xcd; 16]);
            }),
        ];
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree about the ground");
        let parameters = parameters_for(cohort.read_groups());

        let refused = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(a_segmentation_under(ground(false), catalog)).expect("one handle"),
            parameters,
        )
        .expect_err("the assemblies differ");
        let RunError::PspAgainstAnotherReference { sample, difference } = refused else {
            panic!("the assembly must be named: {refused:?}");
        };
        assert_eq!(sample, "beta");
        assert!(
            difference.contains("walked against the assembly") && difference.contains("cdcdcdcd"),
            "the refusal names the assembly the file was walked against: {difference}",
        );
    }

    /// **Direct mode's own refusal, which psp mode needs for the same reason.** A catalog built
    /// on another build of the assembly puts every repeat tract at the wrong position, and
    /// every segment this run loops over is drawn from it.
    #[test]
    fn a_catalog_built_on_another_reference_is_refused() {
        let (_reference_dir, reference) = fixture_reference(true);
        let dir = tempfile::tempdir().expect("a temporary directory");
        // The fixture catalog claims `[7; 16]`, which is not this reference's digest.
        let segmentation = a_segmentation(ground(false));
        let paths = vec![a_psp(&dir, "alpha", &segmentation, |header| {
            header.reference.md5 = reference.info().md5;
            header.contigs[0].md5 = reference.info().contigs[0].md5;
            header.contigs[1].md5 = reference.info().contigs[1].md5;
        })];
        let cohort = OpenPspCohort::open(&paths).expect("one file agrees with itself");
        let parameters = parameters_for(cohort.read_groups());

        let refused = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(a_segmentation(ground(false))).expect("one handle"),
            parameters,
        )
        .expect_err("the catalog was built on something else");
        assert!(
            matches!(refused, RunError::CatalogIsForAnotherReference { .. }),
            "{refused:?}",
        );
    }

    /// The parameters carry no sample names, so what a run can catch is a count — the same
    /// check direct mode makes, on a table this run built from the headers.
    #[test]
    fn parameters_assembled_for_another_cohort_are_refused() {
        let (_dir, _segmentation, paths) = a_cohort();
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree");
        // Parameters for the first sample alone.
        let one = OpenPspCohort::open(&paths[..1]).expect("one file agrees with itself");
        let parameters = parameters_for(one.read_groups());

        let refused = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(a_segmentation(ground(false))).expect("one handle"),
            parameters,
        )
        .expect_err("the parameters are for one sample and the run has two");
        let RunError::ParametersAreForAnotherCohort {
            counted,
            in_the_parameters,
            in_the_run,
        } = refused
        else {
            panic!("the counts must be named: {refused:?}");
        };
        assert_eq!(counted, "the number of samples");
        assert_eq!((in_the_parameters, in_the_run), (1, 2));
    }

    /// **A `.fai` holds no bases**, and a record with an empty allele is padded from the
    /// reference — so a run that could not write such a record is refused at the door rather
    /// than at the locus that needs one.
    #[test]
    fn a_reference_that_holds_no_bases_is_refused() {
        let (_dir, _segmentation, paths) = a_cohort();
        let (_reference_dir, reference) = fixture_reference(false);
        let cohort = OpenPspCohort::open(&paths).expect("the two files agree");
        let parameters = parameters_for(cohort.read_groups());

        let refused = open_the_caller(
            cohort,
            &reference,
            Arc::try_unwrap(a_segmentation(ground(false))).expect("one handle"),
            parameters,
        )
        .expect_err("a .fai cannot serve a padding base");
        assert!(
            matches!(refused, RunError::ReferenceHasNoBases),
            "{refused:?}"
        );
    }
}
