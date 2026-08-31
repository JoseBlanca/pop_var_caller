//! ng's calling run — the machinery the two variant callers drive.
//!
//! A calling run reads every sample's locus observations in coordinate order and emits called
//! variants in coordinate order. `doc/devel/ng/spec/run_streaming.md` owns that outer shape
//! (the caller objects, the sources, the VCF writing) and
//! `doc/devel/ng/arch/run_streaming.md` the types; this module is where its parts land as
//! they are built.
//!
//! **Landed so far:** [`cohort_merge`], which turns k samples' observations into one stream of
//! cohort observations; [`segments`]'s [`Segmentation`], the ground every sample of a run
//! walks; and [`callers`]'s [`AlignedFilesVariantCaller`], constructed and checked but not yet
//! iterating.

pub mod callers;
pub mod cohort_merge;
pub mod segments;

pub use callers::{
    AlignedFilesVariantCaller, AlignmentInputs, AssemblyCheckOutcome, MergeParameters,
};
pub use segments::{Segmentation, SegmentationInputs};

use std::path::PathBuf;

use crate::ng::read::input::{AssemblyMismatch, IngestError};
use crate::ng::repeat_catalog::RepeatCatalogError;

/// What can go wrong driving a run.
///
/// **Every variant names the sample, the file or the counts a person can act on.** A run over
/// thousands of samples that says only "it failed" leaves nobody anywhere to look (spec §9).
///
/// **Some of the reason comes from the cause rather than the top line**, so a command reporting
/// one of these must render the whole chain with
/// [`format_error_chain`](crate::error_render::format_error_chain), never `Display` alone. A
/// bare `Display` says which sample would not open; the chain says its index is missing and
/// where it was looked for.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    /// The run's segments could not be read out of the repeat catalog.
    ///
    /// **The path is carried here because the catalog's own error often does not have it**:
    /// a digest mismatch, over-permissive criteria, differing scan weights and a differing
    /// tool version all describe the file without naming it, and a person with several
    /// catalogs on disk cannot act on that.
    #[error("the run's segments could not be read from the repeat catalog {}", path.display())]
    Catalog {
        /// The catalog the run read.
        path: PathBuf,
        /// What reading it hit.
        #[source]
        source: RepeatCatalogError,
    },

    /// One sample's alignment files could not be opened.
    ///
    /// **The sample's name is a field of its own because the wrapped error does not carry
    /// it**: an open failure knows which file it was, not which individual the file holds.
    /// The cause is boxed to keep this type small, since it travels in every `Result` a run
    /// returns.
    #[error("sample {sample}: its alignment files could not be opened")]
    OpeningSample {
        /// The sample as its read groups name it.
        sample: String,
        /// What opening it hit — the file, and why.
        source: Box<IngestError>,
    },

    /// The run was given no alignment files, so it has no cohort to call.
    ///
    /// **Refused rather than answered with an empty output.** Assembling the parameters for a
    /// cohort of none panics inside the pre-pass, so this refusal has to come first; and a VCF
    /// with no samples in it would otherwise look like a finished run.
    #[error(
        "this run was given no alignment files, so it has no samples to call; \
         check the paths or the pattern it was given"
    )]
    NoAlignmentFiles,

    /// The parameters were assembled for a different cohort from the one this run opened.
    ///
    /// **A count, not a name.** The assembled parameters carry no sample names — one number
    /// per sample and one per read group, in the run's order — so a run cannot compare names
    /// even in principle. A supplied file's names *are* matched against the run's, at that
    /// file's own door: `ParametersFile::to_run_parameters_for` refuses naming the position
    /// where the two lists diverge ([`parameters_file.md`](../spec/parameters_file.md) §6).
    /// What is left for a run to catch is parameters assembled for one cohort handed to a
    /// caller opened over another, which nothing else prevents.
    #[error(
        "the parameters were assembled for a different cohort: {counted} is {in_the_parameters} \
         in the parameters and {in_the_run} in this run; re-run the parameter pre-pass for this \
         cohort, or point the run at the file assembled for it"
    )]
    ParametersAreForAnotherCohort {
        /// What was counted, as a noun phrase that reads inside the sentence — "the number of
        /// samples", "the number of read groups with an error-model calibration".
        counted: &'static str,
        /// How many the parameters hold.
        in_the_parameters: usize,
        /// How many this run has.
        in_the_run: usize,
    },

    /// A sample's reads are against a different build of the reference's assembly.
    ///
    /// **Its contigs already have the right names and lengths** — a file whose header did not
    /// agree with the reference never opened. What differs is the bases, and calling this
    /// sample beside the others would compare genotypes against different sequence.
    ///
    /// **The reference is named because two inputs are in play and only one is wrong**, the
    /// same reason [`Catalog`](Self::Catalog) carries its path. Which file, which contig and
    /// which two checksums come from the cause, so the two sentences say different things
    /// rather than one saying the other twice.
    #[error(
        "sample {sample} was not aligned to this run's reference {}",
        reference.display()
    )]
    SampleAlignedToAnotherReference {
        /// The sample as its read groups name it.
        sample: String,
        /// The reference this run is calling against.
        reference: PathBuf,
        /// Which file, which contig, and the two checksums.
        #[source]
        source: AssemblyMismatch,
    },

    /// The run needs more open files than this process is allowed.
    ///
    /// **Named at construction rather than met as `EMFILE` part-way through a genome.**
    /// Raising the limit is the operator's to do, so the message carries the arithmetic it is
    /// asking them to act on: how many alignment files, what each costs, what the run needs
    /// besides them, and the command that raises it (spec §7.1a).
    #[error(
        "this run needs {needed} open files and this process may open {limit}: \
         {alignment_files} alignment files at {per_file} each, plus {allowance} for the \
         reference, the repeat catalog and the output ({samples} samples). \
         Raise the limit with: ulimit -n {needed}, or call fewer samples at once"
    )]
    NotEnoughFileDescriptors {
        /// How many samples the run holds — what an operator counts their cohort in.
        samples: usize,
        /// How many alignment files those samples are spread over. **This is what the
        /// arithmetic is over**: a sample sequenced across four lanes is four files.
        alignment_files: usize,
        /// Descriptors one alignment file needs.
        per_file: u64,
        /// Descriptors the run needs besides its alignment files.
        allowance: u64,
        /// The two above, combined.
        needed: u64,
        /// What this process may open — the soft `RLIMIT_NOFILE`.
        limit: u64,
    },

    /// The repeat catalog was built on a different reference from the one this run calls
    /// against.
    ///
    /// **Checked here because the catalog's own check cannot do it on the ordinary path.**
    /// Opening a catalog compares its digests against the reference's only when the reference
    /// it was handed carries them, and one read from a `.fai` carries none — so on that path a
    /// catalog is admitted on contig names, lengths and order alone. The digests exist once the
    /// FASTA has been read, and this is where they are compared.
    ///
    /// **Silent and genome-wide otherwise.** The catalog's coordinates are where the repeat
    /// tracts are; a catalog from another build of the same assembly routes every tract to the
    /// wrong position, and every one of this run's segments is drawn from it.
    #[error(
        "the repeat catalog was built on a different reference: the catalog's whole-reference \
         checksum is {in_the_catalog} and {}'s is {in_the_run}",
        reference.display()
    )]
    CatalogIsForAnotherReference {
        /// The reference this run calls against.
        reference: PathBuf,
        /// What the catalog says it was built on.
        in_the_catalog: String,
        /// What this run's reference actually is.
        in_the_run: String,
    },

    /// The reference whose checksums the samples were checked against is not the one their
    /// files were opened against.
    ///
    /// **A caller mistake rather than a user's, and refused rather than trusted.** The two are
    /// meant to be one reference at two moments — before and after its FASTA was read — and the
    /// comparison walks them in step, contig by contig. Two genuinely different genomes would
    /// pair contig *i* of a file against contig *i* of something else, reporting a mismatch on
    /// the wrong chromosome or missing one past the end of the shorter list.
    #[error(
        "the reference the samples were checked against is not the one they were opened against: {difference}"
    )]
    ReferenceCheckedAgainstAnotherGenome {
        /// The first contig that differs, and how.
        difference: String,
    },
}
