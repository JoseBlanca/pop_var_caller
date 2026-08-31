//! The object a direct-mode run is: every sample's alignment files open at once, the ground
//! to analyse, and the numbers to call with.
//!
//! **Direct mode** reads the evidence out of the alignment files and writes no intermediate
//! file (`doc/devel/ng/spec/run_streaming.md` §2). It needs the model parameters up front,
//! because they are fitted from the whole cohort and nothing can be called before they exist;
//! a run that still has to fit them belongs in **psp mode**, which walks each sample on its
//! own, writes what it saw to a per-sample file, fits the parameters from those, and calls
//! from the files.
//!
//! What this file holds today is the object and its construction. Iteration — the merge, the
//! call, the records — lands in the milestones after it
//! (`doc/devel/ng/impl_plan/run_driver_direct_mode.md`).
//!
//! **`pub`, though the architecture calls the accessors crate-private machinery.** The steps
//! that consume them — the per-sample walker, the construction checks — do not exist yet, so
//! `pub(crate)` items here would have no consumer and the crate's `-D warnings` gate would
//! reject them as dead code. The intent is the architecture's; narrow this when the walker
//! lands. `cohort_merge` carries the same note for the same reason.

use crate::ng::calling::allele_candidates::CandidateSelectionConfig;
use crate::ng::calling::inference::RunnableCallingLoopConfig;
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::read::filtering::ReadFilterConfig;
use crate::ng::read::input::SampleReads;
use crate::ng::read::input::read_groups::{ReadGroups, SampleReadGroups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::read::input::{AssemblyMismatch, check_assembly};
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::run::cohort_merge::{CohortLocusBuilderRegionsLen, MaxCohortLocusSpan, MinAltReads};
use crate::pop_var_caller::common::format_md5_hex;

use super::RunError;
use super::segments::Segmentation;

/// The alignment files a run reads, and how it reads them.
///
/// **Grouped rather than passed one by one**, because four of the five are answers to one
/// question — what this run's reads are — and travel together into every sample's open. The
/// fifth is read after those opens finish, and says so at its own field.
///
/// **The read-group table is the run's rather than each sample's**, because a read group is
/// one library preparation and that is the grouping the error model keys on: one sample's
/// read group 0 and another's read group 0 are different preparations, so the numbering has to
/// be run-wide or the two would share a fitted error rate.
///
/// **The reference is the run's for a different reason.** One open copy serves every file of
/// every sample; a per-sample copy would hold the genome in memory once per sample.
pub struct AlignmentInputs<'a> {
    /// Every read group of every sample, numbered run-wide. **The order of
    /// [`ReadGroups::read_groups_per_sample`] is the run's sample order**, which every
    /// per-sample structure downstream is indexed by.
    pub read_groups: &'a ReadGroups,
    /// The one copy of the reference every file decodes against.
    pub reference: &'a OpenReference,
    /// Which reads are admitted, applied per file as they are read.
    pub read_filters: ReadFilterConfig,
    /// Build a missing alignment index, **writing it beside the alignment file**, rather than
    /// refusing the file. Fails when that directory is not writable, which is the ordinary
    /// case for a read-only archive mount.
    pub build_index_if_missing: bool,
    /// The reference **once its per-contig checksums are known** — what each sample's own
    /// contig checksums are compared against.
    ///
    /// **The one field here that is not read at any sample's open**: it is used after every
    /// file is open, because the checksums to compare are captured as each one opens.
    ///
    /// **Only the caller can supply it, and only once the background read of the FASTA has
    /// finished.** A reference read from a `.fai` alone and one whose FASTA has not been read
    /// yet are the same value — no checksums anywhere — so nothing here can tell them apart.
    /// Handing over the second makes the check compare nothing; it does not fail, it reports
    /// that it had nothing to compare ([`AlignedFilesVariantCaller::assembly_check`]), which
    /// is what a run report has to say out loud.
    pub reference_with_checksums: &'a ReferenceInfo,
}

/// The merge's knobs that a single-threaded merge takes.
///
/// **Three values covering four of the merge's five run parameters**, because [`MinAltReads`]
/// is itself a floor and a share. The one left out is how many building regions are worked at
/// once, which only means anything once the merge is threaded; this is built against the merge
/// that runs on one thread (owner's ruling, 2026-08-31).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MergeParameters {
    /// How many reference bases one builder's region covers.
    pub cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen,
    /// The widest a cohort locus may be, in reference bases, before it is refused rather than
    /// built.
    pub max_cohort_locus_span: MaxCohortLocusSpan,
    /// How much non-reference evidence a position needs from some sample to be worth calling.
    pub min_alt_reads: MinAltReads,
}

impl MergeParameters {
    /// The shipped defaults, each taken from its own type rather than restated here.
    pub const DEFAULT: Self = Self {
        cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen::DEFAULT,
        max_cohort_locus_span: MaxCohortLocusSpan::DEFAULT,
        min_alt_reads: MinAltReads::DEFAULT,
    };
}

impl Default for MergeParameters {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// A calling run over open alignment files: **reads in, called variants out in genome order**.
///
/// **What it holds for the whole run** is one open [`SampleReads`] per sample, plus the shared
/// read-only state every sample's walk and every call reads (spec §5.1).
///
/// **The open files are the memory bill.** A measured 11 to 15 MiB of live heap **per open
/// alignment file** (slope 12.0 MiB a file, `examples/dhat_ng_open_files.rs`), so a sample
/// sequenced across four lanes costs four times that. Spec §5.1's figures — 0.9 GB at 63
/// samples, 15 GB at a thousand — are at one file a sample; multiply by the mean file count.
/// When that bill does not fit, the answer is psp mode and the failure must say so.
///
/// **Not `Clone`**: two of these over one cohort would open every file twice.
pub struct AlignedFilesVariantCaller {
    /// One per sample, in the run's sample order — the order
    /// [`ReadGroups::read_groups_per_sample`] gave.
    samples: Vec<SampleReads>,
    /// Every read group of every sample, numbered run-wide.
    ///
    /// **Kept, not borrowed and dropped.** The fitted parameters are keyed by read-group
    /// identifier and hold no table of their own, so writing the parameters file beside the
    /// output and reporting what the run used both need this table back. Re-deriving it there
    /// would mint a second numbering, which is the accident the run's one sample order exists
    /// to prevent.
    read_groups: ReadGroups,
    /// The one copy of the reference the whole run reads against.
    ///
    /// **Kept rather than borrowed at open**, because every walk needs it again: a cursor is
    /// built with a reference accessor taken from here. Cloning shares the cache rather than
    /// copying it.
    reference: OpenReference,
    /// Which reads are admitted. Kept for the same reason as the reference: each sample's walk
    /// opens its cursors with it.
    read_filters: ReadFilterConfig,
    /// The ground, computed once and shared by every sample's walk.
    segmentation: Segmentation,
    /// Every number the pre-pass fitted, frozen for the run.
    parameters: RunParameters,
    /// How the calling loop runs — already validated, because
    /// [`RunnableCallingLoopConfig`] is the only shape a checked configuration takes.
    calling_loop_config: RunnableCallingLoopConfig,
    /// Which alleles a locus is called over.
    candidate_selection: CandidateSelectionConfig,
    /// What the merge admits and how wide a locus it will build.
    merge_parameters: MergeParameters,
    /// What the assembly check could do at construction.
    assembly_check: AssemblyCheckOutcome,
}

impl AlignedFilesVariantCaller {
    /// Open every sample's alignment files and hold them for the run.
    ///
    /// **Named `open` rather than `new` because that is what it does** — every file of every
    /// sample is opened, validated and index-checked here, before a single read flows. A
    /// cohort that cannot be opened fails at once, naming the sample.
    ///
    /// **One `SampleReads` per entry of [`ReadGroups::read_groups_per_sample`]**, in that
    /// order, which is the run's sample order. That is what a cohort tool is required to do
    /// with a run-wide read-group table rather than opening each sample's paths on their own —
    /// the rule [`SampleReads::open_only_sample`] states for every tool that is not
    /// single-sample.
    pub fn open(
        alignments: AlignmentInputs<'_>,
        segmentation: Segmentation,
        parameters: RunParameters,
        calling_loop_config: RunnableCallingLoopConfig,
        candidate_selection: CandidateSelectionConfig,
        merge_parameters: MergeParameters,
    ) -> Result<Self, RunError> {
        let per_sample = alignments.read_groups.read_groups_per_sample();

        // **Three refusals before a single file is opened**, because each of them condemns the
        // whole run and opening a thousand files first would only make the message slower.
        refuse_an_empty_cohort(per_sample)?;
        refuse_parameters_assembled_for_another_cohort(&parameters, alignments.read_groups)?;
        refuse_without_descriptor_headroom(alignments.read_groups)?;

        let mut samples = Vec::with_capacity(per_sample.len());
        for sample in per_sample {
            let reads = SampleReads::open(
                sample,
                alignments.read_groups,
                alignments.reference,
                alignments.read_filters,
                alignments.build_index_if_missing,
            )
            .map_err(|source| RunError::OpeningSample {
                sample: sample.sample.to_string(),
                source: Box::new(source),
            })?;
            samples.push(reads);
        }

        // **And one after**, because it reads the `@SQ M5` tags each open captured.
        let reference_path = alignments
            .reference_with_checksums
            .fasta_path
            .as_deref()
            .or(alignments.reference.info().fasta_path.as_deref())
            .unwrap_or_else(|| std::path::Path::new(A_REFERENCE_WITH_NO_PATH));
        refuse_two_references_that_are_not_one(
            alignments.reference.info(),
            alignments.reference_with_checksums,
        )?;
        refuse_a_catalog_built_on_another_reference(
            &segmentation,
            alignments.reference_with_checksums,
            reference_path,
        )?;
        let assembly_check = check_every_sample_against_the_reference(
            &samples,
            per_sample,
            alignments.reference_with_checksums,
            reference_path,
        )?;

        Ok(Self {
            assembly_check,
            samples,
            read_groups: alignments.read_groups.clone(),
            reference: alignments.reference.clone(),
            read_filters: alignments.read_filters,
            segmentation,
            parameters,
            calling_loop_config,
            candidate_selection,
            merge_parameters,
        })
    }

    /// How many samples this run calls — the length every per-sample row downstream has.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    /// Every sample's open files, **in the run's sample order**.
    pub fn samples(&self) -> impl Iterator<Item = &SampleReads> {
        self.samples.iter()
    }

    /// The sample names, in the run's sample order: index `i` here is the sample every
    /// per-sample structure of this run means by `i`.
    pub fn sample_names(&self) -> impl Iterator<Item = &str> {
        self.samples().map(SampleReads::sample_name)
    }

    /// One sample's open files, by its position in the run's sample order.
    #[must_use]
    pub fn sample_reads(&self, sample_index: usize) -> Option<&SampleReads> {
        self.samples.get(sample_index)
    }

    /// Every read group of every sample, numbered run-wide.
    #[must_use]
    pub fn read_groups(&self) -> &ReadGroups {
        &self.read_groups
    }

    /// The ground this run analyses, and the record of what it was computed from.
    #[must_use]
    pub fn segmentation(&self) -> &Segmentation {
        &self.segmentation
    }

    /// The numbers the pre-pass fitted, which every call is scored against.
    #[must_use]
    pub fn parameters(&self) -> &RunParameters {
        &self.parameters
    }

    /// How the calling loop is configured for this run.
    #[must_use]
    pub fn calling_loop_config(&self) -> &RunnableCallingLoopConfig {
        &self.calling_loop_config
    }

    /// Which alleles a locus of this run is called over.
    #[must_use]
    pub fn candidate_selection(&self) -> &CandidateSelectionConfig {
        &self.candidate_selection
    }

    /// What the merge of this run admits and bounds.
    #[must_use]
    pub fn merge_parameters(&self) -> MergeParameters {
        self.merge_parameters
    }

    /// The reference every sample's walk reads against.
    #[must_use]
    pub fn reference(&self) -> &OpenReference {
        &self.reference
    }

    /// Which reads this run admits.
    #[must_use]
    pub fn read_filters(&self) -> ReadFilterConfig {
        self.read_filters
    }

    /// Whether this run could check what assembly its samples were aligned to, and over how
    /// much.
    ///
    /// **A run report has to say this out loud.** "No sample was aligned to a wrong assembly"
    /// and "no sample could be checked" are different facts, and only one of them is
    /// reassuring.
    #[must_use]
    pub fn assembly_check(&self) -> AssemblyCheckOutcome {
        self.assembly_check
    }
}

/// **The sample names and the sizes, not the contents.** A derived `Debug` would print every
/// open file, every segment and every fitted number — megabytes for a real cohort, in a
/// message someone is reading to find out which run this is.
impl std::fmt::Debug for AlignedFilesVariantCaller {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("AlignedFilesVariantCaller")
            .field("samples", &self.sample_names().collect::<Vec<_>>())
            .field("segments", &self.segmentation.segments().len())
            .finish_non_exhaustive()
    }
}

/// What a run could learn about the assembly its samples were aligned to.
///
/// **Two different facts, and only one is reassuring**: *every sample agreed with the reference*
/// is a check that ran, where *nothing could be checked* is a check that did not. A run report
/// that printed one for the other would be telling an operator their cohort is sound when
/// nothing looked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblyCheckOutcome {
    /// Every checksum that could be compared was compared, and all agreed.
    EverySampleMatchedTheReference {
        /// How many alignment files were looked at.
        alignment_files: usize,
        /// How many contig checksums were actually compared, summed over those files.
        checksums_compared: usize,
        /// How many could have been, had every contig carried a checksum on both sides:
        /// the files times the reference's contig count. **The denominator is what makes the
        /// first number mean something** — "1,386 of 1,512" is a fact, "1,386" alone is not.
        checksums_possible: usize,
    },
    /// Not one checksum could be compared, so no sample's assembly was checked at all.
    NothingCouldBeChecked {
        /// Which side had none — the side an operator would have to change.
        because: NoChecksums,
    },
}

/// Why no contig checksum could be compared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NoChecksums {
    /// The reference carries none. Three ways to get here, and an operator can act on all
    /// three: it was read from a `.fai` alone, which describes a genome's geometry and holds
    /// no bases to checksum; or its FASTA was read but the run did not wait for that read to
    /// finish; or the run asked for the index to be trusted without checking it
    /// (`ReferenceCheck::TrustIndexWithoutChecking`).
    TheReferenceCarriesNone,
    /// The reference has checksums and no alignment file offers one to compare against.
    /// **Ordinary**: `@SQ M5` is optional and plenty of BAM and CRAM files omit it.
    TheAlignmentFilesCarryNone,
}

/// A run with no alignment files is refused rather than answered with an empty output.
fn refuse_an_empty_cohort(per_sample: &[SampleReadGroups]) -> Result<(), RunError> {
    if per_sample.is_empty() {
        return Err(RunError::NoAlignmentFiles);
    }
    Ok(())
}

/// The parameters must have been assembled for **this** cohort.
///
/// **A count, not a match by name**, because [`RunParameters`] carries no names: one inbreeding
/// coefficient per sample of the run and one calibration per read group, both in the run's own
/// order. A supplied file's names are matched against the run's at that file's own door. What is
/// left, and what nothing else prevents, is parameters assembled for one cohort handed to a
/// caller opened over another.
fn refuse_parameters_assembled_for_another_cohort(
    parameters: &RunParameters,
    read_groups: &ReadGroups,
) -> Result<(), RunError> {
    let samples = read_groups.read_groups_per_sample().len();
    let coefficients = parameters.inbreeding_coefficient_by_sample().len();
    if coefficients != samples {
        return Err(RunError::ParametersAreForAnotherCohort {
            counted: "the number of samples",
            in_the_parameters: coefficients,
            in_the_run: samples,
        });
    }

    let of_the_run = read_groups.len();
    let of_the_parameters = parameters.read_group_count();
    if of_the_parameters != of_the_run {
        return Err(RunError::ParametersAreForAnotherCohort {
            counted: "the number of read groups with an error-model calibration",
            in_the_parameters: of_the_parameters,
            in_the_run: of_the_run,
        });
    }
    Ok(())
}

/// Descriptors one alignment file needs while a run holds it open.
///
/// **Measured, on the 63-accession tomato cohort** — `examples/ng_open_cohort_descriptors.rs`,
/// which opens that cohort through [`AlignedFilesVariantCaller::open`] and counts
/// `/proc/self/fd` at three points. The number is right and spec §7.1a's reason for it — "a CRAM
/// and its index are two descriptors each" — is not, so what was counted is written here instead:
///
/// - **63 open alignment files cost one descriptor between them** (3 → 4). An index is parsed
///   into memory at open and an open
///   [`AlignmentFile`](crate::ng::read::input::open_bam::AlignmentFile) keeps no handle, so
///   neither half of the spec's sentence is where the cost is.
/// - **A cursor costs 2 a file** (4 → 130 for 63 cursors, 2.00 a file): one for the file's own
///   reader, and one for the reference accessor [`SampleReads::cursor`] mints per file, which
///   opens the FASTA. A run holds one cursor per file for the whole walk (spec §5.1), so this is
///   the shape the refusal has to size for.
///
/// **Milestone E can move it**: spec §11's question 2 puts several callers in flight, and nobody
/// has counted what that opens. Re-run the probe there.
const DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS: u64 = 2;

/// Descriptors a run needs for everything that is not an alignment file: the three standard
/// streams, the reference and its index, the repeat catalog, and the output and its index —
/// eight — plus 24 of slack for whatever the runtime holds open on its own.
const DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES: u64 = 32;

/// A run that would run out of file descriptors refuses now, naming the arithmetic.
///
/// **The count is of files, not of samples**, because a sample sequenced across four lanes is
/// four files. At two a file, a Linux soft limit of 1,024 is reached at 496 files and macOS's
/// 256 at 112 — while the memory bill spec §7.2 budgets is 500 kB a sample, so 1,000 samples is
/// 500 MB. The descriptors bind first, by a wide margin.
fn refuse_without_descriptor_headroom(read_groups: &ReadGroups) -> Result<(), RunError> {
    // `current` is the soft limit — what this process may open now.
    let limit = rustix::process::getrlimit(rustix::process::Resource::Nofile).current;
    refuse_if_more_descriptors_are_needed_than_allowed(read_groups, limit)
}

/// The decision, with the limit passed in.
///
/// **Split from the syscall so the refusal itself has a test.** The limit a machine reports is
/// far above any cohort a test can build, so a check that read it inline would have a branch no
/// fixture could reach — and an unreachable refusal is one nobody has read the message of.
///
/// `None` means the platform reports no limit at all, and then there is nothing to refuse
/// against.
fn refuse_if_more_descriptors_are_needed_than_allowed(
    read_groups: &ReadGroups,
    limit: Option<u64>,
) -> Result<(), RunError> {
    let files: std::collections::BTreeSet<&std::path::Path> = read_groups
        .iter()
        .map(|(_, read_group)| read_group.file.as_ref())
        .collect();
    let alignment_files = files.len();
    let needed = descriptors_needed_for(alignment_files);

    let Some(limit) = limit else {
        return Ok(());
    };

    if needed > limit {
        return Err(RunError::NotEnoughFileDescriptors {
            samples: read_groups.read_groups_per_sample().len(),
            alignment_files,
            per_file: DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS,
            allowance: DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES,
            needed,
            limit,
        });
    }
    Ok(())
}

/// How many descriptors a run over `alignment_files` files needs, allowance included.
fn descriptors_needed_for(alignment_files: usize) -> u64 {
    alignment_files as u64 * DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS
        + DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES
}

/// The repeat catalog the segments came from must have been built on this run's reference.
///
/// **The catalog's own open cannot check this on the ordinary path.** It compares digests only
/// where the reference it was handed has them, and a reference read from a `.fai` has none — so
/// there the catalog is admitted on contig names, lengths and order alone. Once the FASTA has
/// been read the comparison is possible, and both values are in hand here: the catalog's
/// whole-reference checksum is not optional, and the reference now has one too.
///
/// **What it prevents is silent and genome-wide.** A catalog built on another build of the same
/// assembly puts every repeat tract at the wrong position, and every segment this run walks is
/// drawn from it.
fn refuse_a_catalog_built_on_another_reference(
    segmentation: &Segmentation,
    with_checksums: &ReferenceInfo,
    reference_path: &std::path::Path,
) -> Result<(), RunError> {
    let Some(of_the_run) = with_checksums.md5 else {
        return Ok(());
    };
    let of_the_catalog = segmentation.inputs().catalog.reference_md5;
    if of_the_catalog != of_the_run {
        return Err(RunError::CatalogIsForAnotherReference {
            reference: reference_path.to_path_buf(),
            in_the_catalog: format_md5_hex(of_the_catalog),
            in_the_run: format_md5_hex(of_the_run),
        });
    }
    Ok(())
}

/// The two references must be one reference at two moments.
///
/// **Checked because the comparison downstream walks them in step.** `check_assembly` pairs a
/// file's contig `i` — numbered by the reference the file was *opened* against — with contig
/// `i` of the reference carrying the checksums. It asserts that parity in debug and zips in
/// release, so two genuinely different genomes would pair a file's chromosome against something
/// else's, blaming the wrong contig or missing a real mismatch past the end of the shorter
/// list.
fn refuse_two_references_that_are_not_one(
    opened_against: &ReferenceInfo,
    with_checksums: &ReferenceInfo,
) -> Result<(), RunError> {
    if opened_against.contigs.len() != with_checksums.contigs.len() {
        return Err(RunError::ReferenceCheckedAgainstAnotherGenome {
            difference: format!(
                "one describes {} contigs and the other {}",
                opened_against.contigs.len(),
                with_checksums.contigs.len(),
            ),
        });
    }
    for (opened, checked) in opened_against.contigs.iter().zip(&with_checksums.contigs) {
        if opened.name != checked.name {
            return Err(RunError::ReferenceCheckedAgainstAnotherGenome {
                difference: format!(
                    "contig {} is '{}' in one and '{}' in the other",
                    opened.name, opened.name, checked.name,
                ),
            });
        }
        if opened.length != checked.length {
            return Err(RunError::ReferenceCheckedAgainstAnotherGenome {
                difference: format!(
                    "contig '{}' is {} bases in one and {} in the other",
                    opened.name, opened.length, checked.length,
                ),
            });
        }
    }
    Ok(())
}

/// What a refusal names when the reference has checksums but no path to point at — unreachable
/// while checksums only come from reading a FASTA, and readable rather than blank if it ever is.
const A_REFERENCE_WITH_NO_PATH: &str = "(a reference read from its index alone)";

/// Compare every sample's contig checksums against the reference's.
///
/// **What this catches is a sample aligned to a different build of the same assembly** — every
/// contig the right name and the right length, and different bases. Calling it beside the others
/// would compare genotypes against different sequence.
///
/// **It covers one case, and it is the ordinary one: a run whose reference was read from a
/// `.fai`.** The open gate already compares these same checksums — but only when the reference
/// it was handed carries them, and on the `.fai` path it does not: that read hands back the
/// contig table at once and reads the FASTA on a background thread, so the files open against a
/// reference with nothing to compare. The checksums exist only once that read has finished,
/// which is why this is deferred to here and why it needs a reference that has them rather than
/// the one the files were opened against.
///
/// A run that read its reference from the FASTA directly is already covered at the gate — a
/// wrong-assembly file never opens — and this walks it again, one pass over every file's contig
/// list, to say how much was covered.
fn check_every_sample_against_the_reference(
    samples: &[SampleReads],
    per_sample: &[SampleReadGroups],
    with_checksums: &ReferenceInfo,
    reference_path: &std::path::Path,
) -> Result<AssemblyCheckOutcome, RunError> {
    debug_assert_eq!(
        samples.len(),
        per_sample.len(),
        "one open per sample of the run, in the run's order"
    );

    if with_checksums
        .contigs
        .iter()
        .all(|contig| contig.md5.is_none())
    {
        return Ok(AssemblyCheckOutcome::NothingCouldBeChecked {
            because: NoChecksums::TheReferenceCarriesNone,
        });
    }

    let mut alignment_files = 0;
    let mut checksums_compared = 0;
    for (reads, sample) in samples.iter().zip(per_sample) {
        for (path, observed) in reads.assembly_inputs() {
            let checked = check_assembly(path, observed, with_checksums).map_err(
                |source: AssemblyMismatch| RunError::SampleAlignedToAnotherReference {
                    sample: sample.sample.to_string(),
                    reference: reference_path.to_path_buf(),
                    source,
                },
            )?;
            alignment_files += 1;
            checksums_compared += checked.compared;
        }
    }

    if checksums_compared == 0 {
        return Ok(AssemblyCheckOutcome::NothingCouldBeChecked {
            because: NoChecksums::TheAlignmentFilesCarryNone,
        });
    }

    Ok(AssemblyCheckOutcome::EverySampleMatchedTheReference {
        alignment_files,
        checksums_compared,
        checksums_possible: alignment_files * with_checksums.contigs.len(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::calling::inference::CallingLoopConfig;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, header, matching_contigs, named_bam, read_named_with_length,
        read_named_with_length_in_read_group,
    };
    use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion};
    use crate::ng::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};
    use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare};
    use crate::ng::tandem_repeat::ScanParams;
    use crate::ng::types::{ContigId, GenomeRegion, MapQual, Ploidy, Position};
    use crate::regions::ContigBounds;
    use std::num::NonZeroU32;
    use std::path::PathBuf;
    use tempfile::TempDir;

    /// A one-read indexed BAM whose one read group names `sample`, under `file_name`.
    ///
    /// **Every file opened together must have a different name**: these fixtures declare no
    /// `LB`, so a read group's library is synthesized from its sample, its `@RG ID` and its
    /// file's name, and two files sharing all three are indistinguishable. Real inputs differ
    /// for the same reason, so this matches reality rather than working around a check.
    pub(super) fn bam_for(sample: &str, file_name: &str) -> (TempDir, PathBuf) {
        let (dir, path) = unindexed_bam_for(sample, file_name);
        index(&path);
        (dir, path)
    }

    /// The same file with no index beside it — what a run refuses when it was not asked to
    /// build one.
    pub(super) fn unindexed_bam_for(sample: &str, file_name: &str) -> (TempDir, PathBuf) {
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records = vec![read_named_with_length(&format!("{stem}-r0"), 0, 1, 30)];
        let header = header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg1", Some(sample))],
        );
        named_bam(&header, &records, file_name)
    }

    /// **One indexed BAM holding two samples**, first `first` and then `second`.
    ///
    /// This is the fixture that makes "the run's sample order is the read-group table's" a
    /// claim a test can fail. With one sample per file, first-seen order and the order the
    /// paths were listed are the same sequence, so nothing could tell them apart. Here the
    /// order is decided *inside* one file, by which read group the header declares first.
    pub(super) fn bam_holding_two_samples(
        first: &str,
        second: &str,
        file_name: &str,
    ) -> (TempDir, PathBuf) {
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records = vec![
            read_named_with_length_in_read_group(&format!("{stem}-a"), 0, 1, 30, "rg1"),
            read_named_with_length_in_read_group(&format!("{stem}-b"), 0, 2, 30, "rg2"),
        ];
        let header = header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg1", Some(first)), ("rg2", Some(second))],
        );
        let (dir, path) = named_bam(&header, &records, file_name);
        index(&path);
        (dir, path)
    }

    fn index(path: &PathBuf) {
        crate::bam::index_preflight::preflight_alignment_indexes(std::slice::from_ref(path), true)
            .expect("build index");
    }

    pub(super) fn catalog_header() -> RepeatCatalogHeader {
        RepeatCatalogHeader {
            contigs: Vec::new(),
            reference_md5: [7; 16],
            built_under: StrRepeatCriteria::default(),
            scan: ScanParams::default(),
            tool_version: "test".to_string(),
            longest_tract_bp: Vec::new(),
        }
    }

    /// A segmentation over one short contig, with one ordinary stretch in it, whose catalog
    /// claims it was built on `reference_md5`.
    ///
    /// **The digest matters**: a run whose reference carries a checksum refuses a catalog built
    /// on a different one, so a test that hands over a checksummed reference has to hand over a
    /// catalog that agrees with it.
    pub(super) fn segmentation_built_on(reference_md5: [u8; 16]) -> Segmentation {
        let bounds = [ContigBounds {
            name: "chr1",
            length: 100,
        }];
        let segments = vec![TypedRegion {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(100),
            },
            kind: RegionKind::Generic,
        }];
        Segmentation::build(
            segments.into_iter().map(Ok),
            GenomeRegions::whole_contigs(&bounds),
            RepeatCatalogHeader {
                reference_md5,
                ..catalog_header()
            },
            StrRepeatCriteria::default(),
            PathBuf::from("/genomes/test.catalog.parquet"),
        )
        .expect("a clean stream builds")
    }

    /// The segmentation the tests that never reach the catalog check use.
    pub(super) fn segmentation() -> Segmentation {
        segmentation_built_on([7; 16])
    }

    // -----------------------------------------------------------------
    // Settings that are NOT their type's default
    //
    // A test that hands in a default and asserts the default back cannot tell "held what it
    // was given" from "replaced with the default" — the shape that let four mutations survive
    // an earlier draft of this suite. Everything below differs from what ships.
    // -----------------------------------------------------------------

    pub(super) fn unusual_read_filters() -> ReadFilterConfig {
        ReadFilterConfig {
            min_mapq: Some(MapQual(37)),
            ..ReadFilterConfig::default()
        }
    }

    fn unusual_candidate_selection() -> CandidateSelectionConfig {
        CandidateSelectionConfig {
            min_allele_support: MinAltReads {
                floor: MinAltObs(NonZeroU32::new(7).expect("not zero")),
                share: MinAltReadShare::new_or_panic(0.11),
            },
            ..CandidateSelectionConfig::DEFAULT
        }
    }

    fn unusual_merge_parameters() -> MergeParameters {
        MergeParameters {
            max_cohort_locus_span: MaxCohortLocusSpan(NonZeroU32::new(37).expect("not zero")),
            ..MergeParameters::DEFAULT
        }
    }

    fn unusual_calling_loop_config() -> RunnableCallingLoopConfig {
        CallingLoopConfig {
            max_passes: NonZeroU32::new(11).expect("not zero"),
            ..CallingLoopConfig::DEFAULT
        }
        .validate()
        .expect("eleven passes is a runnable setting")
    }

    /// Open a caller over the given BAM paths, with every setting deliberately unlike its
    /// default.
    fn open_over(
        paths: &[PathBuf],
        reference: &OpenReference,
    ) -> Result<AlignedFilesVariantCaller, RunError> {
        let read_groups = build_read_groups(paths).expect("the fixtures declare read groups");
        let parameters = RunParameters::of_defaults(
            &read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        );
        AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                reference_with_checksums: reference.info(),
            },
            segmentation(),
            parameters,
            unusual_calling_loop_config(),
            unusual_candidate_selection(),
            unusual_merge_parameters(),
        )
    }

    /// A cohort opens, and the object knows how many samples it is calling and what they are
    /// called. **The names are not in alphabetical order**, so a caller that sorted them would
    /// fail here rather than pass by coincidence.
    #[test]
    fn a_cohort_of_three_opens_and_names_its_samples() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("zeta", "a.bam");
        let (_b_dir, b) = bam_for("alpha", "b.bam");
        let (_c_dir, c) = bam_for("mu", "c.bam");

        let caller = open_over(&[a, b, c], &reference).expect("three readable samples open");

        assert_eq!(caller.sample_count(), 3);
        assert_eq!(
            caller.sample_names().collect::<Vec<_>>(),
            vec!["zeta", "alpha", "mu"],
        );
    }

    /// **The run's sample order is the read-group table's, and it is decided inside a file.**
    ///
    /// One BAM declares `zeta` before `alpha`, so first-seen order is neither alphabetical nor
    /// the order of the path list — there is only one path. A caller that sorted, or that
    /// re-derived an order of its own, cannot pass this.
    ///
    /// It matters because three different sample numberings meet in the calling loop, and a
    /// mismatch between them produces wrong genotypes rather than a crash.
    #[test]
    fn the_sample_order_is_decided_by_the_read_group_table_not_by_the_paths() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_dir, both) = bam_holding_two_samples("zeta", "alpha", "both.bam");

        let caller = open_over(std::slice::from_ref(&both), &reference).expect("both open");

        assert_eq!(
            caller.sample_names().collect::<Vec<_>>(),
            vec!["zeta", "alpha"]
        );
        assert_eq!(
            caller.sample_reads(0).expect("in range").sample_name(),
            "zeta",
        );
        assert_eq!(
            caller.sample_reads(1).expect("in range").sample_name(),
            "alpha",
        );
    }

    /// One sample per individual, whatever the file count: two files of one sample are one
    /// open, not two.
    #[test]
    fn two_files_of_one_sample_are_one_sample() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_first_dir, first) = bam_for("NA12878", "first.bam");
        let (_second_dir, second) = bam_for("NA12878", "second.bam");

        let caller = open_over(&[first, second], &reference).expect("one sample, two files");

        assert_eq!(caller.sample_count(), 1);
        assert_eq!(caller.sample_names().collect::<Vec<_>>(), vec!["NA12878"]);
    }

    /// **A file that cannot be opened fails at construction, and the whole reason reaches the
    /// person.** The wrapper names the sample; the cause names the file and says the index is
    /// missing. A run over a thousand samples that failed with an operating-system error
    /// part-way through the genome would leave nobody anywhere to look.
    #[test]
    fn a_sample_whose_index_is_missing_is_refused_naming_the_sample_and_the_file() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_good_dir, good) = bam_for("NA12878", "good.bam");
        let (_bad_dir, bad) = unindexed_bam_for("NA12891", "bad.bam");

        let error =
            open_over(&[good, bad.clone()], &reference).expect_err("the second has no index");

        match &error {
            RunError::OpeningSample { sample, .. } => {
                assert_eq!(sample, "NA12891", "the failing sample is named");
            }
            other => panic!("expected OpeningSample, got {other:?}"),
        }

        let rendered = format_error_chain(&error);
        assert!(rendered.contains("NA12891"), "names the sample: {rendered}");
        assert!(
            rendered.contains(&bad.display().to_string()),
            "names the file: {rendered}",
        );
        assert!(
            rendered.contains("index"),
            "says what is wrong with it: {rendered}",
        );
    }

    /// **Every setting comes back as it was handed in, and none of them is its default.**
    ///
    /// A run that quietly substituted a default here would call every position under the
    /// shipped thresholds while a person read their own numbers off the command line — wrong
    /// genotypes, no failure.
    #[test]
    fn every_setting_comes_back_as_it_was_given() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");

        let caller = open_over(std::slice::from_ref(&a), &reference).expect("opens");

        assert_eq!(caller.read_filters(), unusual_read_filters());
        assert_eq!(*caller.candidate_selection(), unusual_candidate_selection());
        assert_eq!(caller.merge_parameters(), unusual_merge_parameters());
        assert_eq!(*caller.calling_loop_config(), unusual_calling_loop_config());

        // And each of those really is unlike what ships, or the assertions above would hold
        // for a caller that ignored its arguments entirely.
        assert_ne!(unusual_read_filters(), ReadFilterConfig::default());
        assert_ne!(
            unusual_candidate_selection(),
            CandidateSelectionConfig::DEFAULT
        );
        assert_ne!(unusual_merge_parameters(), MergeParameters::DEFAULT);
        assert_ne!(
            unusual_calling_loop_config(),
            RunnableCallingLoopConfig::default(),
        );
    }

    /// The ground and the read-group table come back too.
    #[test]
    fn the_ground_and_the_read_group_table_come_back() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");

        let paths = [a];
        let expected_table = build_read_groups(&paths).expect("read groups");
        let caller = open_over(&paths, &reference).expect("opens");

        assert_eq!(caller.segmentation().segments().len(), 1);
        assert_eq!(caller.segmentation().inputs().catalog, catalog_header());
        assert_eq!(caller.segmentation().analysed_regions().len(), 1);
        assert_eq!(*caller.read_groups(), expected_table);
    }

    /// **The shipped defaults are each type's own**, not numbers restated here. Restating one
    /// would let the merge's default and the caller's drift apart in a later sweep.
    #[test]
    fn the_default_merge_parameters_are_each_types_default() {
        assert_eq!(
            MergeParameters::DEFAULT.cohort_locus_builder_regions_len,
            CohortLocusBuilderRegionsLen::DEFAULT,
        );
        assert_eq!(
            MergeParameters::DEFAULT.max_cohort_locus_span,
            MaxCohortLocusSpan::DEFAULT,
        );
        assert_eq!(MergeParameters::DEFAULT.min_alt_reads, MinAltReads::DEFAULT);
        assert_eq!(MergeParameters::default(), MergeParameters::DEFAULT);
    }

    /// The debug rendering names the samples and counts the segments. **It must not print
    /// their contents**: a real cohort's would be megabytes, in a line someone is reading to
    /// find out which run this is.
    #[test]
    fn the_debug_rendering_is_names_and_sizes() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");

        let caller = open_over(std::slice::from_ref(&a), &reference).expect("opens");

        let rendered = format!("{caller:?}");
        assert!(rendered.contains("NA12878"), "{rendered}");
        assert!(rendered.contains("segments: 1"), "{rendered}");
        assert!(
            !rendered.contains("Generic"),
            "the segments themselves are not printed: {rendered}",
        );
    }
}

#[cfg(test)]
mod construction_checks {
    use super::tests::{bam_for, segmentation_built_on, unusual_read_filters};
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, header, matching_contigs, named_bam, read_named_with_length,
    };
    use crate::ng::types::Ploidy;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn parameters_for(read_groups: &ReadGroups) -> RunParameters {
        RunParameters::of_defaults(
            read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        )
    }

    /// Open, with the parameters and the verified reference chosen separately — which is what
    /// lets a test hand over a mismatched pair.
    fn open_with(
        paths: &[PathBuf],
        reference: &OpenReference,
        verified: &ReferenceInfo,
        parameters: RunParameters,
    ) -> Result<AlignedFilesVariantCaller, RunError> {
        let read_groups = build_read_groups(paths).expect("read groups");
        AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                reference_with_checksums: verified,
            },
            // **The catalog claims the reference it is actually given**, so these tests reach
            // the checks they are about rather than stopping at the catalog one. The catalog
            // check has its own fixtures below.
            segmentation_built_on(verified.md5.unwrap_or([7; 16])),
            parameters,
            RunnableCallingLoopConfig::default(),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
    }

    /// **A run given no alignment files is refused, not answered with an empty output.** The
    /// likeliest way to reach it is a file pattern that matched nothing, and a VCF with no
    /// samples in it looks like a finished run.
    #[test]
    fn a_cohort_with_no_alignment_files_is_refused() {
        let (_reference_dir, reference) = fixture_reference(false);

        // **The parameters are another cohort's, and they have to be.** Assembling parameters
        // for a cohort of none is not possible — it panics inside the pre-pass's own assembly,
        // which is the failure this refusal exists to come before. So the fixture shows two
        // things at once: the empty cohort is refused, and it is refused *first*.
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let of_one = build_read_groups(std::slice::from_ref(&a)).expect("read groups");

        let error = open_with(&[], &reference, reference.info(), parameters_for(&of_one))
            .expect_err("no files is no cohort");

        assert!(matches!(error, RunError::NoAlignmentFiles), "{error:?}");
        assert!(
            error.to_string().contains("no alignment files"),
            "the message says what was missing: {error}",
        );
    }

    /// **Parameters assembled for another cohort are refused, with both counts.** Nothing else
    /// prevents it: the parameters and the read-group table reach `open` as separate arguments.
    #[test]
    fn parameters_assembled_for_another_cohort_are_refused_with_both_counts() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let (_b_dir, b) = bam_for("NA12891", "b.bam");

        // Parameters for a cohort of two, handed to a run over one.
        let two = build_read_groups(&[a.clone(), b]).expect("read groups");
        let error = open_with(
            std::slice::from_ref(&a),
            &reference,
            reference.info(),
            parameters_for(&two),
        )
        .expect_err("two samples' parameters do not fit a run of one");

        match &error {
            RunError::ParametersAreForAnotherCohort {
                counted,
                in_the_parameters,
                in_the_run,
            } => {
                assert_eq!(*counted, "the number of samples");
                assert_eq!(*in_the_parameters, 2);
                assert_eq!(*in_the_run, 1);
            }
            other => panic!("expected ParametersAreForAnotherCohort, got {other:?}"),
        }

        // **The rendered message, not just the fields.** Dropping a count from the format
        // string would leave the fields right and the person with nothing, and this step's
        // whole product is what the person reads.
        let message = error.to_string();
        assert!(
            message.contains("the number of samples is 2 in the parameters and 1 in this run"),
            "both counts, in a sentence that reads at one sample: {message}",
        );
        assert!(
            message.contains("re-run the parameter pre-pass"),
            "says what to do next: {message}",
        );
    }

    /// The same cohort's own parameters are accepted — the check must not refuse what it is
    /// meant to let through.
    #[test]
    fn a_cohorts_own_parameters_are_accepted() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");

        let paths = [a];
        let read_groups = build_read_groups(&paths).expect("read groups");
        open_with(
            &paths,
            &reference,
            reference.info(),
            parameters_for(&read_groups),
        )
        .expect("its own parameters fit");
    }

    // -----------------------------------------------------------------
    // The descriptor headroom
    // -----------------------------------------------------------------

    /// **The refusal names the limit, the count and what to do about it.** Raising the limit is
    /// the operator's to do, so a message that does not say so leaves them nowhere.
    #[test]
    fn a_run_that_would_run_out_of_descriptors_is_refused_naming_both_numbers() {
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let read_groups = build_read_groups(std::slice::from_ref(&a)).expect("read groups");

        let error = refuse_if_more_descriptors_are_needed_than_allowed(&read_groups, Some(4))
            .expect_err("one file needs more than four descriptors, allowance included");

        match &error {
            RunError::NotEnoughFileDescriptors {
                samples,
                alignment_files,
                needed,
                limit,
                ..
            } => {
                assert_eq!(*samples, 1);
                assert_eq!(*alignment_files, 1);
                assert_eq!(*needed, descriptors_needed_for(1));
                assert_eq!(*limit, 4);
            }
            other => panic!("expected NotEnoughFileDescriptors, got {other:?}"),
        }

        // **The rendered message, not the fields.** This step's whole product is what a person
        // reads, so the arithmetic they are asked to check has to be in the string: what it
        // needs, what it may open, and the two numbers those come from.
        let message = error.to_string();
        assert!(
            message.contains("needs 34 open files"),
            "names what it needs: {message}",
        );
        assert!(
            message.contains("may open 4"),
            "names the limit, not just a digit of it: {message}",
        );
        assert!(
            message.contains("1 alignment files at 2 each"),
            "shows the arithmetic behind the total: {message}",
        );
        assert!(
            message.contains("plus 32 for the reference"),
            "and the part that is not the files: {message}",
        );
        assert!(
            message.contains("ulimit -n 34"),
            "gives the command with the number in it: {message}",
        );
    }

    /// A limit above what the run needs lets it through, and a platform reporting no limit at
    /// all is not a refusal.
    #[test]
    fn headroom_above_what_the_run_needs_is_no_refusal() {
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let read_groups = build_read_groups(std::slice::from_ref(&a)).expect("read groups");

        refuse_if_more_descriptors_are_needed_than_allowed(
            &read_groups,
            Some(descriptors_needed_for(1)),
        )
        .expect("exactly enough is enough");
        refuse_if_more_descriptors_are_needed_than_allowed(&read_groups, None)
            .expect("no limit reported is no refusal");
    }

    /// **The count is of files and not of samples**, because a sample sequenced across four
    /// lanes is four files and eight descriptors.
    #[test]
    fn the_descriptor_count_grows_with_files_not_with_samples() {
        assert_eq!(
            descriptors_needed_for(1) + 2 * DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS,
            descriptors_needed_for(3),
            "each further file costs the file and its index",
        );
        assert_eq!(
            descriptors_needed_for(0),
            DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES,
            "with no alignment files, only the run's own allowance is needed",
        );
    }

    // -----------------------------------------------------------------
    // The assembly check
    // -----------------------------------------------------------------

    /// A one-read indexed BAM declaring `checksums[i]` as contig `i`'s `@SQ M5`.
    pub(super) fn bam_declaring_checksums(
        sample: &str,
        file_name: &str,
        checksums: &[String],
    ) -> (TempDir, PathBuf) {
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records = vec![read_named_with_length(&format!("{stem}-r0"), 0, 1, 30)];
        let contigs: Vec<(&str, usize, Option<&str>)> = matching_contigs()
            .into_iter()
            .zip(checksums)
            .map(|((name, length, _), checksum)| (name, length, Some(checksum.as_str())))
            .collect();
        let header = header(Some("coordinate"), &contigs, &[("rg1", Some(sample))]);
        let (dir, path) = named_bam(&header, &records, file_name);
        crate::bam::index_preflight::preflight_alignment_indexes(std::slice::from_ref(&path), true)
            .expect("build index");
        (dir, path)
    }

    /// The reference's own per-contig checksums, as `@SQ M5` hex.
    pub(super) fn checksums_of(verified: &ReferenceInfo) -> Vec<String> {
        verified
            .contigs
            .iter()
            .map(|contig| {
                format_md5_hex(
                    contig
                        .md5
                        .expect("the fasta arm carries per-contig checksums"),
                )
            })
            .collect()
    }

    /// The same checksums with one contig's replaced — a file aligned to a different build of the
    /// same assembly: right names, right lengths, different bases.
    pub(super) fn checksums_with_one_wrong(verified: &ReferenceInfo) -> Vec<String> {
        let mut checksums = checksums_of(verified);
        checksums[1] = "ffffffffffffffffffffffffffffffff".to_string();
        checksums
    }

    /// **The only case `check_assembly` can catch, and the reason it is deferred.**
    ///
    /// A run reading a reference with a `.fai` beside it gets the contig table at once and the
    /// checksums on a background thread, so the files open against a reference that carries
    /// none and the open gate cannot compare them. They exist only once that read has
    /// finished, and only then can a sample's own be checked. Opened against the FASTA
    /// instead, the gate does this itself and the file never opens — a different test, one
    /// layer down.
    #[test]
    fn a_sample_aligned_to_another_assembly_is_refused_naming_the_sample() {
        let (_open_dir, opened_against) = fixture_reference(false);
        let (_verified_dir, verified) = fixture_reference(true);
        let (_bad_dir, bad) = bam_declaring_checksums(
            "NA12891",
            "bad.bam",
            &checksums_with_one_wrong(verified.info()),
        );

        let read_groups = build_read_groups(std::slice::from_ref(&bad)).expect("read groups");
        let error = open_with(
            std::slice::from_ref(&bad),
            &opened_against,
            verified.info(),
            parameters_for(&read_groups),
        )
        .expect_err("its checksums are not this reference's");

        let named_reference = match &error {
            RunError::SampleAlignedToAnotherReference {
                sample, reference, ..
            } => {
                assert_eq!(sample, "NA12891");
                reference.clone()
            }
            other => panic!("expected SampleAlignedToAnotherReference, got {other:?}"),
        };

        // **Two references are in play and only one is wrong**, so the refusal has to say which
        // one the run was calling against — the same reason the catalog failure carries its
        // path.
        assert_eq!(
            Some(named_reference.as_path()),
            verified.info().fasta_path.as_deref(),
            "the reference named is the one whose checksums it was compared against",
        );

        let rendered = format_error_chain(&error);
        assert!(rendered.contains("NA12891"), "names the sample: {rendered}");
        assert!(
            rendered.contains(&named_reference.display().to_string()),
            "names the reference: {rendered}",
        );
        assert!(
            rendered.contains(&bad.display().to_string()),
            "names the file: {rendered}",
        );
        assert!(
            rendered.contains("ffffffffffffffffffffffffffffffff"),
            "shows the checksum that did not match: {rendered}",
        );
        assert!(
            rendered.contains("chr2"),
            "names the contig that differs, not just the file: {rendered}",
        );
        // **The chain must not say one thing twice.** The wrapper names the sample and the
        // reference; the cause names the file and the contigs. A wrapper that repeated the
        // cause's own sentence would make the message longer without making it say more.
        assert_eq!(
            rendered.matches("aligned to a different assembly").count(),
            1,
            "the diagnosis appears once, in the cause that can name the contig: {rendered}",
        );
    }

    /// A sample whose checksums are the reference's own passes, and the run records how much
    /// was actually compared — out of how much could have been.
    #[test]
    fn a_sample_carrying_the_references_checksums_passes_and_the_run_says_what_it_compared() {
        let (_open_dir, opened_against) = fixture_reference(false);
        let (_verified_dir, verified) = fixture_reference(true);
        let (_good_dir, good) =
            bam_declaring_checksums("NA12878", "good.bam", &checksums_of(verified.info()));

        let read_groups = build_read_groups(std::slice::from_ref(&good)).expect("read groups");
        let caller = open_with(
            std::slice::from_ref(&good),
            &opened_against,
            verified.info(),
            parameters_for(&read_groups),
        )
        .expect("its checksums are the reference's");

        let contigs = verified.info().contigs.len();
        assert_eq!(
            caller.assembly_check(),
            AssemblyCheckOutcome::EverySampleMatchedTheReference {
                alignment_files: 1,
                checksums_compared: contigs,
                checksums_possible: contigs,
            },
            "every contig carried a checksum on both sides, so every one was compared",
        );
    }

    /// **A reference with no checksums reports that it checked nothing, rather than passing
    /// silently** — and says which side had none, because that is the side an operator can
    /// change. "No sample was aligned to a wrong assembly" and "no sample could be checked" are
    /// different facts.
    #[test]
    fn a_reference_without_checksums_reports_that_nothing_could_be_checked() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");

        let paths = [a];
        let read_groups = build_read_groups(&paths).expect("read groups");
        let caller = open_with(
            &paths,
            &reference,
            reference.info(),
            parameters_for(&read_groups),
        )
        .expect("a fai-only reference is not a failure");

        assert_eq!(
            caller.assembly_check(),
            AssemblyCheckOutcome::NothingCouldBeChecked {
                because: NoChecksums::TheReferenceCarriesNone,
            },
        );
    }

    /// **And when the reference has checksums but no alignment file does, the run says that
    /// instead** — the ordinary case, since `@SQ M5` is optional and plenty of files omit it.
    /// Reporting it as "compared" over some number of files would be the substitution this
    /// type exists to prevent: a check that looked at nothing, printed as a check that passed.
    #[test]
    fn files_without_checksums_report_that_nothing_could_be_checked() {
        let (_open_dir, opened_against) = fixture_reference(false);
        let (_verified_dir, verified) = fixture_reference(true);
        // `bam_for` builds its header from `matching_contigs()`, which declares no `M5`.
        let (_a_dir, a) = bam_for("NA12878", "a.bam");

        let paths = [a];
        let read_groups = build_read_groups(&paths).expect("read groups");
        let caller = open_with(
            &paths,
            &opened_against,
            verified.info(),
            parameters_for(&read_groups),
        )
        .expect("files without checksums are ordinary, not a failure");

        assert_eq!(
            caller.assembly_check(),
            AssemblyCheckOutcome::NothingCouldBeChecked {
                because: NoChecksums::TheAlignmentFilesCarryNone,
            },
        );
    }

    /// **The cheap refusals run before any file is opened**, so a cohort that fails one of them
    /// does not pay to open a thousand files first.
    ///
    /// Shown where the two outcomes genuinely differ: the run's one file has no index, so
    /// opening it *would* fail — and the parameters are for another cohort, so the check that
    /// runs first is the one whose error comes back. A caller that opened first would report
    /// the open failure instead.
    #[test]
    fn the_cheap_refusals_come_before_the_files_are_opened() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_bad_dir, bad) = super::tests::unindexed_bam_for("NA12878", "bad.bam");
        let (_b_dir, b) = bam_for("NA12891", "b.bam");

        // Parameters for two samples, a run of one, and that one file unopenable.
        let two = build_read_groups(&[bad.clone(), b]).expect("read groups");
        let error = open_with(
            std::slice::from_ref(&bad),
            &reference,
            reference.info(),
            parameters_for(&two),
        )
        .expect_err("the parameters do not fit, and the file would not open either");

        assert!(
            matches!(error, RunError::ParametersAreForAnotherCohort { .. }),
            "the parameter mismatch is reported, so no file was opened: {error:?}",
        );
    }
}

#[cfg(test)]
mod checks_that_needed_their_own_fixtures {
    //! **Every test here exists because a deliberate defect survived without it.**
    //!
    //! The first draft's fixtures were one sample, one file, one read group — a shape in which
    //! "samples", "files" and "read groups" are the same number, both arity checks fail
    //! together, and reversing a two-element pairing is the identity. Seven mutations lived in
    //! that blind spot. Each test below names the one it kills.

    use super::tests::{bam_for, catalog_header, segmentation_built_on, unusual_read_filters};
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, header, matching_contigs, named_bam,
        read_named_with_length_in_read_group,
    };
    use crate::ng::repeat_catalog::RepeatCatalogHeader;
    use crate::ng::types::Ploidy;
    use crate::pop_var_caller::common::format_md5_hex;
    use std::path::PathBuf;
    use tempfile::TempDir;

    fn parameters_for(read_groups: &ReadGroups) -> RunParameters {
        RunParameters::of_defaults(
            read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        )
    }

    /// Open over a table already built, with the two references chosen separately.
    fn open_over(
        read_groups: &ReadGroups,
        opened_against: &OpenReference,
        with_checksums: &ReferenceInfo,
        parameters: RunParameters,
    ) -> Result<AlignedFilesVariantCaller, RunError> {
        AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups,
                reference: opened_against,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                reference_with_checksums: with_checksums,
            },
            segmentation_built_on(with_checksums.md5.unwrap_or([7; 16])),
            parameters,
            RunnableCallingLoopConfig::default(),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
    }

    /// One indexed BAM whose two read groups both name `sample` — one sample, two read groups,
    /// one file. The shape that separates the two arity counts.
    fn bam_with_two_read_groups(sample: &str, file_name: &str) -> (TempDir, PathBuf) {
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records = vec![
            read_named_with_length_in_read_group(&format!("{stem}-a"), 0, 1, 30, "rg1"),
            read_named_with_length_in_read_group(&format!("{stem}-b"), 0, 2, 30, "rg2"),
        ];
        let header = header(
            Some("coordinate"),
            &matching_contigs(),
            &[("rg1", Some(sample)), ("rg2", Some(sample))],
        );
        let (dir, path) = named_bam(&header, &records, file_name);
        crate::bam::index_preflight::preflight_alignment_indexes(std::slice::from_ref(&path), true)
            .expect("build index");
        (dir, path)
    }

    // -----------------------------------------------------------------
    // The two arity counts, each with the other held fixed
    // -----------------------------------------------------------------

    /// **The sample count is checked on its own.** Two samples with one read group each,
    /// against parameters for one sample with two read groups: the read-group counts agree, so
    /// only the sample check can refuse.
    ///
    /// Kills: disabling the inbreeding-coefficient comparison.
    #[test]
    fn the_sample_count_is_checked_with_the_read_group_count_agreeing() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let (_b_dir, b) = bam_for("NA12891", "b.bam");
        let (_one_dir, one_sample_two_groups) = bam_with_two_read_groups("NA12892", "both.bam");

        let run = build_read_groups(&[a, b]).expect("read groups");
        let other =
            build_read_groups(std::slice::from_ref(&one_sample_two_groups)).expect("read groups");
        assert_eq!(run.len(), other.len(), "the read-group counts must agree");
        assert_ne!(
            run.read_groups_per_sample().len(),
            other.read_groups_per_sample().len(),
            "and the sample counts must not",
        );

        let error = open_over(&run, &reference, reference.info(), parameters_for(&other))
            .expect_err("two samples' run, one sample's parameters");

        match &error {
            RunError::ParametersAreForAnotherCohort {
                counted,
                in_the_parameters,
                in_the_run,
            } => {
                assert_eq!(*counted, "the number of samples");
                assert_eq!(*in_the_parameters, 1);
                assert_eq!(*in_the_run, 2);
            }
            other => panic!("expected the sample count to be named, got {other:?}"),
        }
    }

    /// **The read-group count is checked on its own.** One sample with two read groups, against
    /// parameters for one sample with one: the sample counts agree, so only the read-group
    /// check can refuse.
    ///
    /// Kills: disabling the calibration comparison.
    #[test]
    fn the_read_group_count_is_checked_with_the_sample_count_agreeing() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_one_dir, two_groups) = bam_with_two_read_groups("NA12878", "two.bam");
        let (_a_dir, one_group) = bam_for("NA12891", "one.bam");

        let run = build_read_groups(std::slice::from_ref(&two_groups)).expect("read groups");
        let other = build_read_groups(std::slice::from_ref(&one_group)).expect("read groups");
        assert_eq!(
            run.read_groups_per_sample().len(),
            other.read_groups_per_sample().len(),
            "the sample counts must agree",
        );
        assert_ne!(run.len(), other.len(), "and the read-group counts must not");

        let error = open_over(&run, &reference, reference.info(), parameters_for(&other))
            .expect_err("two read groups' run, one read group's parameters");

        match &error {
            RunError::ParametersAreForAnotherCohort {
                counted,
                in_the_parameters,
                in_the_run,
            } => {
                assert!(
                    counted.contains("read groups"),
                    "the read-group count is named: {counted}",
                );
                assert_eq!(*in_the_parameters, 1);
                assert_eq!(*in_the_run, 2);
            }
            other => panic!("expected the read-group count to be named, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------
    // The descriptor count is over files
    // -----------------------------------------------------------------

    /// **One sample across two files needs two files' worth of descriptors.**
    ///
    /// Kills: counting samples instead of files. With one sample and two files the two answers
    /// differ, which no earlier fixture could show.
    #[test]
    fn one_sample_across_two_files_is_counted_as_two_files() {
        let (_first_dir, first) = bam_for("NA12878", "first.bam");
        let (_second_dir, second) = bam_for("NA12878", "second.bam");
        let read_groups = build_read_groups(&[first, second]).expect("read groups");
        assert_eq!(read_groups.read_groups_per_sample().len(), 1, "one sample");

        let error = refuse_if_more_descriptors_are_needed_than_allowed(&read_groups, Some(4))
            .expect_err("two files need more than four descriptors");

        match &error {
            RunError::NotEnoughFileDescriptors {
                samples,
                alignment_files,
                needed,
                ..
            } => {
                assert_eq!(*samples, 1);
                assert_eq!(*alignment_files, 2);
                assert_eq!(
                    *needed,
                    descriptors_needed_for(2),
                    "two files, not one sample",
                );
            }
            other => panic!("expected NotEnoughFileDescriptors, got {other:?}"),
        }
    }

    /// **And two samples sharing one file are counted as one file.** The other direction, so
    /// neither count can stand in for the other.
    #[test]
    fn two_samples_in_one_file_are_counted_as_one_file() {
        let (_dir, both) = super::tests::bam_holding_two_samples("zeta", "alpha", "both.bam");
        let read_groups = build_read_groups(std::slice::from_ref(&both)).expect("read groups");
        assert_eq!(read_groups.read_groups_per_sample().len(), 2, "two samples");

        let error = refuse_if_more_descriptors_are_needed_than_allowed(&read_groups, Some(4))
            .expect_err("one file still needs more than four descriptors");

        match &error {
            RunError::NotEnoughFileDescriptors {
                samples,
                alignment_files,
                needed,
                ..
            } => {
                assert_eq!(*samples, 2);
                assert_eq!(*alignment_files, 1);
                assert_eq!(
                    *needed,
                    descriptors_needed_for(1),
                    "one file, not two samples",
                );
            }
            other => panic!("expected NotEnoughFileDescriptors, got {other:?}"),
        }
    }

    // -----------------------------------------------------------------
    // The assembly loop: which sample is blamed, and how much was counted
    // -----------------------------------------------------------------

    /// **The refusal names the sample the bad file belongs to, and not another.**
    ///
    /// Kills: pairing the opens against the sample list reversed. With one sample a reversal is
    /// the identity; with three, and the wrong assembly on the *second*, it is not.
    #[test]
    fn the_refusal_names_the_sample_whose_file_is_wrong_not_another() {
        let (_open_dir, opened_against) = fixture_reference(false);
        let (_verified_dir, verified) = fixture_reference(true);
        let right = super::construction_checks::checksums_of(verified.info());
        let wrong = super::construction_checks::checksums_with_one_wrong(verified.info());

        let (_a_dir, a) =
            super::construction_checks::bam_declaring_checksums("first", "a.bam", &right);
        let (_b_dir, b) =
            super::construction_checks::bam_declaring_checksums("second", "b.bam", &wrong);
        let (_c_dir, c) =
            super::construction_checks::bam_declaring_checksums("third", "c.bam", &right);

        let read_groups = build_read_groups(&[a, b, c]).expect("read groups");
        assert_eq!(
            read_groups
                .read_groups_per_sample()
                .iter()
                .map(|sample| sample.sample.as_ref())
                .collect::<Vec<&str>>(),
            vec!["first", "second", "third"],
            "the middle sample is the one with the wrong checksums",
        );

        let error = open_over(
            &read_groups,
            &opened_against,
            verified.info(),
            parameters_for(&read_groups),
        )
        .expect_err("the second sample is on another assembly");

        match &error {
            RunError::SampleAlignedToAnotherReference { sample, .. } => {
                assert_eq!(
                    sample, "second",
                    "the sample blamed is the one whose file is wrong",
                );
            }
            other => panic!("expected SampleAlignedToAnotherReference, got {other:?}"),
        }
    }

    /// **A sample across two files is two files' worth of comparisons.**
    ///
    /// Kills: pinning the file counter to one rather than accumulating it.
    #[test]
    fn a_sample_across_two_files_reports_both_of_them() {
        let (_open_dir, opened_against) = fixture_reference(false);
        let (_verified_dir, verified) = fixture_reference(true);
        let right = super::construction_checks::checksums_of(verified.info());
        let (_first_dir, first) =
            super::construction_checks::bam_declaring_checksums("NA12878", "one.bam", &right);
        let (_second_dir, second) =
            super::construction_checks::bam_declaring_checksums("NA12878", "two.bam", &right);

        let read_groups = build_read_groups(&[first, second]).expect("read groups");
        let caller = open_over(
            &read_groups,
            &opened_against,
            verified.info(),
            parameters_for(&read_groups),
        )
        .expect("both files carry the reference's checksums");

        let contigs = verified.info().contigs.len();
        assert_eq!(
            caller.assembly_check(),
            AssemblyCheckOutcome::EverySampleMatchedTheReference {
                alignment_files: 2,
                checksums_compared: 2 * contigs,
                checksums_possible: 2 * contigs,
            },
            "one sample, two files, every contig compared in each",
        );
    }

    // -----------------------------------------------------------------
    // The catalog against the run's reference
    // -----------------------------------------------------------------

    /// **A repeat catalog built on another reference is refused**, and it is refused here
    /// because nothing else can: opening the catalog compares digests only when the reference
    /// it was handed carries them, and on the `.fai` path it carries none.
    ///
    /// Without this the segments — every repeat tract's coordinates — would be drawn from a
    /// different build of the genome and applied silently, genome-wide.
    #[test]
    fn a_catalog_built_on_another_reference_is_refused() {
        let (_open_dir, opened_against) = fixture_reference(false);
        let (_verified_dir, verified) = fixture_reference(true);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let read_groups = build_read_groups(std::slice::from_ref(&a)).expect("read groups");

        let of_the_run = verified.info().md5.expect("the fasta arm has one");
        let of_another = [0xab; 16];
        assert_ne!(of_the_run, of_another);

        let error = AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &opened_against,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                reference_with_checksums: verified.info(),
            },
            segmentation_built_on(of_another),
            parameters_for(&read_groups),
            RunnableCallingLoopConfig::default(),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
        .expect_err("the catalog was built on a different genome");

        match &error {
            RunError::CatalogIsForAnotherReference {
                in_the_catalog,
                in_the_run,
                ..
            } => {
                assert_eq!(*in_the_catalog, format_md5_hex(of_another));
                assert_eq!(*in_the_run, format_md5_hex(of_the_run));
            }
            other => panic!("expected CatalogIsForAnotherReference, got {other:?}"),
        }

        let rendered = format_error_chain(&error);
        assert!(
            rendered.contains("built on a different reference"),
            "says what is wrong: {rendered}",
        );
    }

    /// **A reference with no checksums cannot condemn a catalog**, so the check stands aside
    /// rather than refusing what it cannot judge.
    #[test]
    fn a_reference_without_checksums_does_not_condemn_the_catalog() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let read_groups = build_read_groups(std::slice::from_ref(&a)).expect("read groups");

        AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &reference,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                reference_with_checksums: reference.info(),
            },
            segmentation_built_on([0xab; 16]),
            parameters_for(&read_groups),
            RunnableCallingLoopConfig::default(),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
        .expect("nothing to compare against is not a disagreement");
    }

    /// **Two different genomes are refused rather than walked in step.** The comparison
    /// downstream pairs contig `i` with contig `i`; two genomes would pair a file's chromosome
    /// against something else's.
    #[test]
    fn a_checked_reference_that_is_not_the_opened_one_is_refused() {
        let (_open_dir, opened_against) = fixture_reference(false);
        let (_big_dir, another_genome) =
            crate::ng::read::input::test_fixtures::big_fixture_reference();
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let read_groups = build_read_groups(std::slice::from_ref(&a)).expect("read groups");

        let error = AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &opened_against,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                reference_with_checksums: another_genome.info(),
            },
            segmentation_built_on([7; 16]),
            parameters_for(&read_groups),
            RunnableCallingLoopConfig::default(),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
        .expect_err("those are two different genomes");

        assert!(
            matches!(error, RunError::ReferenceCheckedAgainstAnotherGenome { .. }),
            "{error:?}",
        );
    }

    /// The catalog header fixture is still what the other modules build on.
    #[test]
    fn the_catalog_fixture_carries_the_digest_it_is_given() {
        let built = segmentation_built_on([9; 16]);
        assert_eq!(built.inputs().catalog.reference_md5, [9; 16]);
        assert_eq!(
            catalog_header().built_under,
            built.inputs().catalog.built_under,
            "only the digest differs from the shared header",
        );
        let _: RepeatCatalogHeader = catalog_header();
    }
}
