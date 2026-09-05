//! ng's runs — the machinery both calling modes drive, and psp mode's walk stage.
//!
//! A calling run reads every sample's locus observations in coordinate order and emits called
//! variants in coordinate order; psp mode's walk stage runs the same per-sample machinery and
//! stores the observations instead. `doc/devel/ng/spec/run_streaming.md` owns that outer shape
//! (the run objects, the sources, the VCF writing) and
//! `doc/devel/ng/arch/run_streaming.md` the types; this module is where its parts land as
//! they are built.
//!
//! **Landed so far:** [`cohort_merge`], which turns k samples' observations into one stream of
//! cohort observations; [`segments`]'s [`Segmentation`], the ground every sample of a run
//! walks; [`walker`]'s [`AlignmentFilesWalker`], one sample's alignment files behind the merge's
//! source interface; [`callers`]'s [`AlignedFilesVariantCaller`], which drives that merge
//! over one walker per sample and genotypes each cohort locus where it is built;
//! [`records`], which turns each called locus into what a VCF record states; [`report`]'s
//! [`RunReport`], what a finished run says about itself; [`gatherer`]'s
//! [`SampleObservationGatherer`], psp mode's walk stage — one sample's observations drained
//! into a psp on disk (spec §5.2); and [`psp_source`]'s [`PspObservationSource`], the same
//! merge interface served from a psp instead of from alignment files — the one place the two
//! calling modes differ (spec §3.1).
//!
//! **Three ways out, and only two are a real run's.**
//! [`AlignedFilesVariantCaller::call_cohort`] collects every called locus and hands them back at
//! once, which is what an oracle wants and what no real run can afford;
//! [`AlignedFilesVariantCaller::call_cohort_handing_each_record_over`] hands each record over as
//! it is finished and keeps none, which is the path
//! [`call-from-alignments`](crate::pop_var_caller_exp::call_from_alignments) takes; and
//! [`SampleObservationGatherer::write_psp`] stores one sample's walk as a psp, calling
//! nothing — psp mode's walk half, whose calling half reads a cohort of those files back.

pub mod callers;
pub mod census_from_psp;
pub mod cohort_merge;
pub mod gatherer;
pub mod psp_caller;
pub mod psp_source;
pub mod records;
pub mod report;
pub mod segments;
#[cfg(test)]
pub(crate) mod test_fixtures;
#[cfg(test)]
mod tract_junction_ownership;
pub mod walker;

pub use callers::{
    AlignedFilesVariantCaller, AlignmentInputs, AssemblyCheckOutcome, CalledCohort,
    CohortWalkTallies, MergeParameters, SampleWalkTallies, WrittenCohort,
};
pub use census_from_psp::{CensusFromPspError, CensusOfStoredPileup, CensusTally, census_from_psp};
pub use gatherer::{CensusPlan, CensusSelection, SampleObservationGatherer, SampleWalkInputs};
pub use psp_caller::{
    OpenPspCohort, PspVariantCaller, StoredCohortInputs, StoredCohortTallies, StoredSample,
};
pub use psp_source::{PspObservationSource, PspSourceError, StoredSampleTallies};
pub use report::RunReport;
pub use segments::{Segmentation, SegmentationInputs};
pub use walker::{AlignmentFilesWalker, RunSegments};

use std::path::PathBuf;

use crate::ng::read::input::{AssemblyMismatch, IngestError};
use crate::ng::repeat_catalog::RepeatCatalogError;
use crate::ng::types::{GenomePosition, GenomeRegion};

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

    /// One sample's source could not produce its next observation.
    ///
    /// **Which sample, and how far it had got.** Neither alone locates a failure in a run over
    /// thousands of samples: the sample without the position says nothing about where to look,
    /// and the position without the sample says nothing about which file to look in (spec §9).
    /// What went wrong arrives through the cause — a read query that would not open, a
    /// reference fetch past a contig's end — which names the region it failed on in every case
    /// but one, a failure of the region *stream* itself having no region to name
    /// ([`LocusGenerationError::region`](crate::ng::locus_generation::LocusGenerationError::region)).
    ///
    /// **The cause is a boxed trait object because the two modes fail differently**: a walker
    /// fails at reading alignment files
    /// ([`LocusGenerationError`](crate::ng::locus_generation::LocusGenerationError)) and a psp
    /// reader at decoding a file. The merge adds nothing to either and passes it through, so
    /// this variant is where they meet.
    ///
    /// **Two genome coordinates appear in the rendered chain and they mean different things**,
    /// so each says which it is: this line's is how far the sample *succeeded*, and the cause's
    /// is the region that *failed*. Without the distinction a reader sees "…ended at contig 0
    /// position 13: reference fetch over contig 0:10-20 failed" and reasonably takes it for one
    /// fact said twice.
    ///
    /// **It carries no "what to do next", unlike its neighbours here**, and that is the error
    /// chain's shape rather than an omission: [`format_error_chain`](crate::error_render::format_error_chain)
    /// appends each cause after a colon, so an instruction on this line would be buried in the
    /// middle of the sentence, ahead of the thing it is telling the reader to act on. The
    /// variants that end with an instruction — [`NoAlignmentFiles`](Self::NoAlignmentFiles),
    /// [`NotEnoughFileDescriptors`](Self::NotEnoughFileDescriptors) — all have no cause beneath
    /// them. Whatever reports a run is where the advice belongs.
    #[error("sample {sample}: reading its observations failed; {reached}")]
    SourceFailed {
        /// The sample as its read groups name it.
        sample: String,
        /// How far that sample's walk had got before it failed.
        reached: WalkProgress,
        /// What the source hit.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    /// The run's reference describes a genome's geometry but holds none of its bases.
    ///
    /// **A `.fai` is a table of contig names, lengths and offsets, and nothing else.** A run
    /// needs the sequence itself: every locus the walk emits carries the reference bases over its
    /// span, and the read preparer left-aligns against them. So a reference opened from an index
    /// alone can be checked against a cohort's headers and cannot be called against.
    ///
    /// **Refused before a single file is opened**, because it condemns the whole run and nothing
    /// about it improves by discovering it at the first locus, a genome's worth of setup later.
    #[error(
        "this run's reference was opened from a `.fai` index alone, which holds no bases: a run \
         needs the FASTA itself, both for the reference allele at every locus and to left-align \
         the reads. Point the run at the `.fa` beside that index"
    )]
    ReferenceHasNoBases,

    /// The reference's `.fai` index could not be read.
    ///
    /// **Named apart from the reference having no bases at all**: that one is a wrong argument — a
    /// `.fai` handed in where a FASTA was needed — while this is a missing or damaged
    /// `<reference>.fai` beside a FASTA that is otherwise right, which `samtools faidx
    /// <reference>` rebuilds.
    ///
    /// **The instruction is not on this line, and that is the chain's shape rather than an
    /// omission** — the same reason [`SourceFailed`](Self::SourceFailed) carries none: a cause is
    /// appended after a colon, so advice here would sit in front of the thing it is telling the
    /// reader to act on.
    #[error("the index beside this run's reference {} could not be read", reference.display())]
    ReferenceIndexUnreadable {
        /// The FASTA whose `.fai` would not read.
        reference: PathBuf,
        /// What reading it hit.
        #[source]
        source: std::io::Error,
    },

    /// A locus generator would not accept the settings it was built with.
    ///
    /// **A user's mistake, since 2026-09-01**: `AlignmentInputs::locus_generator_settings` is
    /// how a run says how deep to fold each position and how many reads to hold open, so these
    /// are numbers somebody typed. Checked at `AlignedFilesVariantCaller::open`, before a file
    /// is opened, so a cohort of a thousand samples is not opened to be told at its first
    /// locus.
    ///
    /// **The message is the cause's alone.** Wrapping it would put a sentence about locus
    /// generators in front of the one sentence that names the setting and the limit, and the
    /// reader needs the second.
    #[error(transparent)]
    LocusGeneratorSettings {
        /// Which setting, and why it was refused.
        source: crate::ng::locus_generation::pileup::PileupGeneratorConfigError,
    },

    /// The repeat-tract generator's settings contradict the ground's classification.
    ///
    /// One cross-config rule, and it is the only way this fires: the flank the generator
    /// fetches beside a tract must fit inside the radius within which the classification
    /// bundles two tracts together. A flank wider than that reaches into a neighbour the
    /// classification had already decided was far enough away.
    ///
    /// **The message is the cause's alone**, for the reason above.
    #[error(transparent)]
    TractGeneratorSettings {
        /// Which setting, and why it was refused.
        source: crate::ng::locus_generation::ssr::SsrGeneratorConfigError,
    },

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
         {alignment_files} alignment files at {per_file} each, {samples} samples at \
         {per_sample} more each for the reference bases their walks read, and {allowance} for \
         the reference, the repeat catalog and the output. \
         Raise the limit with: ulimit -n {needed}, or call fewer samples at once"
    )]
    NotEnoughFileDescriptors {
        /// How many samples the run holds — what an operator counts their cohort in.
        samples: usize,
        /// How many alignment files those samples are spread over. **This is what the
        /// arithmetic is over**: a sample sequenced across four lanes is four files.
        alignment_files: usize,
        /// Descriptors one alignment file needs: its reader, and the reference accessor the
        /// cursor over it holds.
        per_file: u64,
        /// Descriptors one **sample** needs on top of its files: the two reference accessors
        /// its locus generator holds for the run.
        per_sample: u64,
        /// Descriptors the run needs besides those two terms.
        allowance: u64,
        /// The three above, combined.
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

    /// The reference base a record with an empty allele has to be padded with could not be read.
    ///
    /// **VCF cannot spell an empty allele**, so an insertion's or a deletion's record is written
    /// by prefixing every allele with the reference base beside its span
    /// (`doc/devel/ng/spec/vcf_output.md` §5). The run reads that one base from its own
    /// reference, and where the read fails the record cannot be written at all — production's
    /// repeat-tract writer invents the letter `N` in the one case it cannot read, and this
    /// format deliberately does not port that.
    ///
    /// **The ordinary reachable cause is the FASTA becoming unreadable part-way through a run**,
    /// since the position asked for is inside a contig this run has already been walking.
    #[error("the reference base beside the locus at {locus} could not be read")]
    PaddingBaseUnreadable {
        /// The locus whose record needed the base.
        locus: GenomeRegion,
        /// What the reference fetch hit.
        #[source]
        source: crate::ng::ref_seq::RefSeqError,
    },

    /// Whatever the run was handing its records to would not take one.
    ///
    /// **The locus is named because a run writes hundreds of thousands of them** and the cause —
    /// a full disk, a directory that went away — says nothing about where the file stopped. A
    /// consumer holding the two knows how much of its output is complete.
    ///
    /// **It carries no "what to do next"**, for [`SourceFailed`](Self::SourceFailed)'s reason:
    /// the cause is appended after a colon, so an instruction here would sit in front of the
    /// thing it is telling the reader to act on.
    #[error("the record for the locus at {locus} could not be written")]
    RecordNotWritten {
        /// The locus whose record was refused.
        locus: GenomeRegion,
        /// What refused it.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    /// The files handed to one sample's walk could not be read as one sample's.
    ///
    /// psp mode's walk writes one psp per sample, so a gatherer takes one sample's alignment
    /// files and builds its read-group table over exactly those (spec `run_streaming.md`
    /// §5.2, §6.1 — the table's numbering starts at zero *because* it sees one sample).
    /// Two ways the files can fail that: their read-group headers cannot be read at all, or
    /// they read fine and name more than one sample. Both arrive here, and the source —
    /// [`IngestError`] — lists every file beside the sample it claims, so the stray is
    /// visible.
    ///
    /// **Deliberately not [`OpeningSample`](Self::OpeningSample)**, which names the sample
    /// that failed to open: this failure is that no one sample could be established, so a
    /// variant naming one would name it wrongly.
    #[error("a psp holds one sample, and these alignment files could not be read as one sample's")]
    FilesNotFromOneSample {
        /// What the read-group table build found.
        #[source]
        source: Box<IngestError>,
    },

    /// The sample's psp file itself could not be produced — created at the start, or sealed
    /// at the end.
    ///
    /// The record-by-record failures in between carry
    /// [`RecordNotWritten`](Self::RecordNotWritten), which names the locus; these two
    /// moments have no locus, so what locates them is the path. A file that fails at the
    /// seal may be left half-written at that path — the format guarantees a reader refuses
    /// it as interrupted rather than reading it as whole (`psp_file_format.md` §10).
    #[error("the psp at {} could not be written", path.display())]
    PspNotWritten {
        /// The file that could not be produced.
        path: PathBuf,
        /// What the store refused.
        #[source]
        source: Box<crate::ng::psp::PspWriteError>,
    },

    /// The sample's census file could not be produced, though its psp was.
    ///
    /// **This fails the sample's walk rather than being reported and passed over** (spec §2,
    /// plan step G2). The two files are one product: the census is what a parameters fit reads,
    /// and a psp without one forces the sample to be walked again — which is the single thing
    /// psp mode exists to avoid. A run that left a finished-looking psp beside a missing census
    /// would be storing that re-walk for somebody to discover later.
    ///
    /// **Three different failures reach it and the source says which**: the psp's header would
    /// not encode, the file would not be created, or the census encoder refused. The path is
    /// what locates all three.
    /// This run's census positions could not be chosen, so no sample can build one.
    ///
    /// **Refused before the first sample is walked**, because the selection is the run's and not
    /// the sample's: a run that discovered this at its third sample would have spent two walks.
    #[error("this run's census positions could not be chosen")]
    CensusNotPlanned {
        /// What the selection refused.
        #[source]
        source: Box<crate::ng::parameter_estimation::joint::loci::SelectionError>,
    },

    #[error("the census at {} could not be written", path.display())]
    CensusNotWritten {
        /// The file that could not be produced.
        path: PathBuf,
        /// What refused it.
        #[source]
        source: Box<dyn std::error::Error + Send + Sync>,
    },

    // ---------------------------------------------------------------------
    // psp mode's calling stage: what a cohort of stored files can be wrong about
    // (spec `run_streaming.md` §6.2)
    // ---------------------------------------------------------------------
    /// The run was given no psps, so it has no cohort to call.
    ///
    /// **[`NoAlignmentFiles`](Self::NoAlignmentFiles)'s sibling, not a reuse of it**: the two
    /// name different arguments, and a run told to look in a directory of psps that held none
    /// must be told about psps. The reason for refusing rather than writing an empty file is
    /// the same — assembling the parameters for a cohort of none panics inside the pre-pass.
    #[error(
        "this run was given no psps, so it has no samples to call; check the paths or the \
         directory it was given"
    )]
    NoPsps,

    /// One of the run's psps could not be read.
    ///
    /// **The path and not the sample**, and that is forced rather than chosen: a psp's sample
    /// name lives in its header, and a reader reaches the header third — footer, then index,
    /// then header, the order the layout forces — so a file that was cut short, or whose
    /// footer will not parse, has no name to be reported under. What went wrong comes from the cause, which tells an interrupted walk
    /// (`the writer did not finish`) apart from a damaged one.
    #[error("the psp at {} could not be read", path.display())]
    PspNotRead {
        /// The file that would not read.
        path: PathBuf,
        /// What the store refused.
        #[source]
        source: Box<crate::ng::psp::PspReadError>,
    },

    /// Two psps name the same sample.
    ///
    /// Either a duplicated argument, or a cohort that would call one individual twice and
    /// weight every allele frequency by it (spec §6.2).
    #[error(
        "sample {sample} appears twice, in {} and {}; a cohort holds each individual once, so \
         drop one of the two",
        first.display(),
        second.display()
    )]
    SampleAppearsTwice {
        /// The individual named by both files.
        sample: String,
        /// The first psp naming it.
        first: PathBuf,
        /// The second.
        second: PathBuf,
    },

    /// Two of the run's psps were walked over different ground.
    ///
    /// **The one cohort-wide check, and the one mismatch that would produce a wrong answer**
    /// rather than a missing one (spec §6.2). A sample has no records over ground it never
    /// looked at, and that absence is indistinguishable from *no variant here* — so a cohort
    /// mixing two grounds would call the ground only one sample walked and read the other's
    /// silence as homozygous reference.
    ///
    /// Whether the shared ground could be called instead of refusing is spec §11's question 5.
    #[error(
        "samples {left} and {right} were walked over different ground, so they cannot be \
         called as one cohort; re-walk one of them over the other's regions, or call each over \
         the ground it has"
    )]
    AnalysedRegionsDiffer {
        /// The sample whose ground the run took as the cohort's.
        left: String,
        /// The first sample whose ground differed from it.
        right: String,
    },

    /// A psp was written under a different segmentation from the one this run loops over.
    ///
    /// **Why this is not pedantry**: the observations in a psp were minted inside the segments
    /// of the walking run's segmentation, and *no observation crosses a segment's edge*
    /// (spec §4.3) is true only while the two segmentations are the same. Under a different
    /// catalog or different repeat-tract criteria, a stored repeat-tract observation can
    /// straddle a calling segment's edge, and the independence the calling loop rests on is
    /// gone.
    ///
    /// `field` is [`SegmentationInputs::first_difference`]'s answer, written to read inside
    /// this sentence.
    #[error(
        "the psp for sample {sample} was written under a different {field} from this run's; \
         call it with the catalog and repeat-tract settings it was walked with, or walk it \
         again with this run's"
    )]
    SegmentationInputsDiffer {
        /// The sample whose file disagrees with the run.
        sample: String,
        /// The first field of the segmentation's inputs that differs.
        field: &'static str,
    },

    /// A psp's coordinate space is not the run's reference's.
    ///
    /// **Every record in a psp is written against the contig table its header carries**, in
    /// that order — `ContigId(i)` *is* `contigs[i]` — so a file whose table is not this run's
    /// puts every one of its observations on the wrong chromosome. `difference` names the
    /// first contig that disagrees and how, because a run over a whole genome cannot be
    /// checked by hand.
    ///
    /// **Named apart from [`SampleAlignedToAnotherReference`](Self::SampleAlignedToAnotherReference)**,
    /// which is direct mode's and says something narrower: there the file's contigs are known
    /// to be the reference's and only the *bases* differ. Here the table itself may differ.
    #[error(
        "the psp for sample {sample} was written against a different reference from this \
         run's: {difference}"
    )]
    PspAgainstAnotherReference {
        /// The sample whose file disagrees.
        sample: String,
        /// The first contig that differs, and how.
        difference: String,
    },

    /// A run over stored files needs more open files than this process is allowed.
    ///
    /// **Named at construction rather than met as `EMFILE` at the two-hundred-and-fiftieth
    /// psp**, which is where a macOS default limit lands and where the operating system's own
    /// message would blame whichever file happened to be next (spec §7.1a). The arithmetic is
    /// stated because it is what the operator is being asked to act on, and it is simpler than
    /// direct mode's: a psp costs one descriptor, held for the whole run.
    #[error(
        "this run needs {needed} open files for its {samples} psps and this process may open \
         {limit}: one descriptor a psp, held open for the whole run, and {allowance} for the \
         reference, the repeat catalog and the output. \
         Raise the limit with: ulimit -n {needed}, or call fewer samples at once"
    )]
    NotEnoughFileDescriptorsForPsps {
        /// How many psps the run was given — one per individual.
        samples: usize,
        /// Descriptors the run needs besides its psps.
        allowance: u64,
        /// The two above, combined.
        needed: u64,
        /// What this process may open — the soft `RLIMIT_NOFILE`.
        limit: u64,
    },

    /// A psp's read-group table cannot be renumbered into the run's.
    ///
    /// **Tables are merged, not compared** (spec §6.2): every psp numbers its own read groups
    /// from zero, so identifiers colliding across files is the normal case. What cannot be
    /// merged is a table that does not identify its own groups at all.
    #[error(
        "the read groups of sample {sample}'s psp cannot be renumbered into this run's: {problem}"
    )]
    PspReadGroupsCannotBeMerged {
        /// The sample whose table will not merge.
        sample: String,
        /// What is wrong with it, as a clause that reads inside the sentence.
        problem: String,
    },
}

/// How far a sample's source had got when something went wrong.
///
/// **The second half of locating a failure** (spec §9), and an enum rather than a position
/// because a source that fails on its very first draw has no position to report and must not
/// invent one. A run that said "failed at contig 0:1" when nothing had been read would send an
/// operator to a locus that is innocent.
///
/// The position is the **last base the last observation covers** — `SampleLocusObservations::reach_position`,
/// which is the same value the merge orders on — and not where the source was decoding: a source
/// is ahead of what it has yielded, and how far ahead is the generators' business and nobody
/// else's.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WalkProgress {
    /// Nothing had been yielded yet, so the failure is on the first draw.
    NothingYet,
    /// The last observation yielded ended here.
    After(GenomePosition),
}

/// **Written as a clause that names its own role**, because [`RunError::SourceFailed`] renders
/// it beside a second genome coordinate — the region the cause failed on — and a bare "after
/// contig 0:13" beside "over contig 0:10-20" reads as one fact in two notations. Saying *last
/// complete observation* is what tells the two apart.
///
/// **The contig is printed by index, not by name, and that is a real cost to the reader**:
/// someone whose genome is `SL4.0ch01`…`SL4.0ch13` has to count from zero in the `.fai` to use
/// it. It is not that the names are unreachable — a run's [`Segmentation`] carries the catalog's
/// contig table, [`ContigInfo::name`](crate::ng::reference_info::ContigInfo) and all — but that
/// this type is a position with no reference beside it, so it cannot spend them alone. The place
/// to spend them is wherever a [`RunError`] is rendered for a person, which has the run's
/// reference in hand.
///
/// The spelling is `contig {n} position {p}`, which is what ng's other genome *positions* print
/// ([`vcf::writer`](crate::ng::vcf::writer), [`psp::index`](crate::ng::psp::index)).
/// [`GenomeRegion`]'s `contig {n}:{start}-{end}` is the spelling
/// for a *region*, and keeping the two apart is what stops this line and the cause's from
/// looking like the same kind of thing.
impl std::fmt::Display for WalkProgress {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NothingYet => write!(formatter, "it had produced no observations yet"),
            Self::After(position) => write!(
                formatter,
                "its last complete observation ended at contig {} position {}",
                position.contig.get(),
                position.position.get(),
            ),
        }
    }
}
