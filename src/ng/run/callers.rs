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
//! What this file holds today is the object, its construction, the merge driven over one
//! walker per sample, and the calling joined to it — reads in, one called locus per locus the
//! merge kept. What is still missing is the shape and the emission: the object returns every
//! called locus in one go rather than yielding VCF records as it goes, which is the rest of
//! this milestone and the pool's (`doc/devel/ng/impl_plan/run_driver_direct_mode.md`).
//!
//! **`pub`, though the architecture calls all of this crate-private** (arch §6: three public
//! objects, each an iterator, and nothing else). **The intent is real and acting on it is still
//! blocked, and the block is mechanical**: neither `merge_cohort` nor `call_cohort` has a
//! consumer outside tests until the
//! subcommand lands, so narrowing it makes it dead code and the crate's `-D warnings` gate
//! rejects the build — measured, not assumed, when this step tried it. What *did* narrow is
//! everything reachable only through it: `WalkReference` and `generic_path_generators` are
//! `pub(crate)`. Narrow the rest when the command exists to call it. `cohort_merge` carries the
//! same note for the same reason.

use std::num::NonZeroU32;
use std::sync::Arc;

use crate::ng::calling::allele_candidates::generic::select_generic;
use crate::ng::calling::allele_candidates::ssr::{
    SsrLocusSelection, SsrSelectionConfig, select_ssr,
};
use crate::ng::calling::allele_candidates::{
    AlleleRemap, CandidateSelectionConfig, SelectionVerdict, UnmatchedSupport,
};
use crate::ng::calling::evidence_shaping::{
    GenericEvidenceScratch, SsrEvidenceScratch, shape_generic_locus, shape_ssr_locus,
};
use crate::ng::calling::inference::{LocusGenotyper, RunnableCallingLoopConfig};
use crate::ng::calling::run_parameters::RunParameters;
use crate::ng::calling::{CallingScratch, FrozenParameters, LocusInference};
use crate::ng::locus_generation::pileup::{PileupGeneratorConfig, PileupGeneratorCounts};
use crate::ng::locus_generation::{
    GeneratorCounts, LocusCounts, LocusKind, SequenceObservation, SsrDetail,
};
use crate::ng::read::filtering::ReadFilterConfig;
use crate::ng::read::filtering::ReadFilterCounts;
use crate::ng::read::input::SampleReads;
use crate::ng::read::input::read_groups::{ReadGroups, SampleReadGroups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::read::input::{AssemblyMismatch, check_assembly};
use crate::ng::reference_info::ReferenceInfo;
use crate::ng::run::cohort_merge::build::{CohortObservation, RegionOutcome};
use crate::ng::run::cohort_merge::observation_cache::ObservationCache;
use crate::ng::run::cohort_merge::serial::{
    merge_cohort_handing_each_locus_over,
    merge_cohort_handing_each_locus_over_covering_samples_in_parallel, merge_cohort_through_cache,
};
use crate::ng::run::cohort_merge::{CohortLocusBuilderRegionsLen, MaxCohortLocusSpan, MinAltReads};
use crate::ng::types::{GenomeRegion, ReadGroupId};
use crate::ng::vcf::VcfRecord;
use crate::ng::vcf::assemble::assemble_record;
use crate::pop_var_caller::common::format_md5_hex;

use super::RunError;
use super::records::{
    a_written_genotype_carries_an_alternative, evidence_for_output, padding_base_beside,
};
use super::segments::Segmentation;
use super::walker::{AlignmentFilesWalker, RunSegments, WalkReference, generic_path_generators};

/// The alignment files a run reads, and how it reads them.
///
/// **Grouped rather than passed one by one**, because all six answer one question — how this
/// run turns its files into evidence. Four of them travel together into every sample's open;
/// the other two are read afterwards, and each says so at its own field.
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
    /// **The five knobs the locus generator walks with** — the two per-column depth caps, the
    /// widest record footprint, the mate-lookup window and the ceiling on reads held open at
    /// once.
    ///
    /// **Here beside the read filters because the two answer one question**: how this run
    /// turns bytes into evidence. The filters decide which reads are admitted; these decide
    /// how many of the admitted ones a position is scored on.
    ///
    /// **Three of the five are the depth axis, and it is not a formality.** The ceiling on
    /// reads held open was 4,096 until 2026-08-05, and at that value one ~130× tomato
    /// chromosome silently refused 19,725 reads
    /// ([`PileupGeneratorConfig::max_active_reads`]). A cohort deeper than the constants were
    /// set for needs a run that can raise them, and one shallower gains nothing by holding
    /// them high.
    ///
    /// **Checked at [`AlignedFilesVariantCaller::open`]**, with the rest of the refusals,
    /// rather than at the first locus a generator is built for.
    pub locus_generator_settings: PileupGeneratorConfig,
    /// Build a missing alignment index, **writing it beside the alignment file**, rather than
    /// refusing the file. Fails when that directory is not writable, which is the ordinary
    /// case for a read-only archive mount.
    pub build_index_if_missing: bool,
    /// The reference **once its per-contig checksums are known** — what each sample's own
    /// contig checksums are compared against.
    ///
    /// **Not read at any sample's open**: it is used after every file is open, because the
    /// checksums to compare are captured as each one opens. The settings above are the other
    /// field of this struct no sample's open reads — they are read where the walkers are
    /// built.
    ///
    /// **Only the caller can supply it, and only once the background read of the FASTA has
    /// finished.** A reference read from a `.fai` alone and one whose FASTA has not been read
    /// yet are the same value — no checksums anywhere — so nothing here can tell them apart.
    /// Handing over the second makes the check compare nothing; it does not fail, it reports
    /// that it had nothing to compare ([`AlignedFilesVariantCaller::assembly_check`]), which
    /// is what a run report has to say out loud.
    pub reference_with_checksums: &'a ReferenceInfo,
}

/// The merge's knobs that a run's own drivers take.
///
/// **Three values covering four of the merge's five run parameters**, because [`MinAltReads`]
/// is itself a floor and a share. The one left out is how many building regions are worked at
/// once, which only means anything under the merge's own region batching — which no calling
/// run reaches: since E1 the cover threads across samples, and the building itself stays on
/// one thread (the 2026-08-31 ruling to build against the single-threaded merge, unchanged by
/// where the cover's sweep runs).
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
    /// The reference the walk fetches its bases from, with its index parsed once for the run.
    ///
    /// **Opened at construction, not at the first locus**, so a reference that holds no bases is
    /// a refusal rather than a failure a genome's worth of setup later. Every walker takes three
    /// accessors from it and shares nothing but the index and the contig table.
    walk_reference: WalkReference,
    /// The ground, computed once and shared by every sample's walk.
    ///
    /// **Held behind an `Arc` because the walkers this run will own read it.** A walker keeps a
    /// handle on the segmentation for the whole run (`walker.rs`), so a run that held the
    /// segmentation by value and its walkers beside it would be a struct whose fields borrow one
    /// another — which safe Rust cannot express, and which cloning cannot escape, since
    /// [`Segmentation`] is deliberately not `Clone`. The list is still stored once: what each
    /// walker takes is a reference count, not a copy (owner's ruling, 2026-08-31).
    segmentation: Arc<Segmentation>,
    /// Every number the pre-pass fitted, frozen for the run.
    parameters: RunParameters,
    /// How the calling loop runs — already validated, because
    /// [`RunnableCallingLoopConfig`] is the only shape a checked configuration takes.
    calling_loop_config: RunnableCallingLoopConfig,
    /// Which alleles a locus is called over.
    candidate_selection: CandidateSelectionConfig,
    /// What the merge admits and how wide a locus it will build.
    merge_parameters: MergeParameters,
    /// The settings every one of this run's locus generators is built with — checked at
    /// `open`, so building a generator from them cannot fail later.
    locus_generator_settings: PileupGeneratorConfig,
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

        // **Five refusals before a single alignment file is opened**, because each of them condemns the
        // whole run and opening a thousand files first would only make the message slower.
        refuse_an_empty_cohort(per_sample)?;
        refuse_parameters_assembled_for_another_cohort(&parameters, alignments.read_groups)?;
        refuse_without_descriptor_headroom(alignments.read_groups)?;
        // **The locus generator's settings are checked before anything is built with them**,
        // for the reason the three above are checked before a file opens: a run whose depth
        // caps are impossible is wrong at the door, and a refusal at the first locus would
        // arrive after every file of a thousand-sample cohort had been opened.
        alignments
            .locus_generator_settings
            .check()
            .map_err(|source| RunError::LocusGeneratorSettings { source })?;
        // **The fifth is the walk's own precondition, checked at the door.** Opening the
        // reference for walking parses its index and refuses a reference that carries no bases at
        // all; doing it here rather than at the first locus is the same rule as the other three,
        // and it is the only one whose product the run keeps.
        let walk_reference = WalkReference::of(alignments.reference)?;

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
            walk_reference,
            segmentation: Arc::new(segmentation),
            parameters,
            calling_loop_config,
            candidate_selection,
            merge_parameters,
            locus_generator_settings: alignments.locus_generator_settings,
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

    /// A handle on that ground, for a walker to keep for the whole run.
    ///
    /// **A second accessor rather than a wider return on the first**, because the two answer
    /// different questions: [`segmentation`](Self::segmentation) is for reading it here, and this
    /// is for holding it elsewhere. Handing out the `Arc` from one accessor would make every
    /// reader take shared ownership of the genome-sized list to ask how many segments it has.
    #[must_use]
    pub fn shared_segmentation(&self) -> Arc<Segmentation> {
        Arc::clone(&self.segmentation)
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

    /// The settings every locus generator of this run is built with.
    #[must_use]
    pub fn locus_generator_settings(&self) -> PileupGeneratorConfig {
        self.locus_generator_settings
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

    /// One walker per sample, in the run's sample order, each over the whole segmentation.
    ///
    /// **Consumes the run's open files**, because a walker owns the `SampleReads` it reads and
    /// `SampleReads` is deliberately not `Clone` — one sample's files are opened once and walked
    /// once. Everything else a walker needs is shared: the segments by reference count, and the
    /// reference's index and contig table by one handle each — the accessors built over them are
    /// the sample's own, never shared.
    ///
    /// **One generator set per walker, never one shared between walkers**, because a locus
    /// generator carries state
    /// across segments (spec §8). Each is the generic path filled and both repeat-tract slots
    /// refused as unbuilt, which is what [`generic_path_generators`] documents.
    fn walkers(self) -> Result<RunReadyToWalk, RunError> {
        let mut walkers = Vec::with_capacity(self.samples.len());
        for reads in self.samples {
            let generators = generic_path_generators(
                &self.walk_reference,
                self.locus_generator_settings,
                // The generators derive the tract generator's bundle radius from here —
                // the criteria the ground was actually cut with, by construction.
                self.segmentation.inputs(),
            )?;
            walkers.push(AlignmentFilesWalker::over(
                Arc::clone(&self.segmentation),
                reads,
                generators,
            ));
        }
        Ok(RunReadyToWalk {
            segmentation: self.segmentation,
            merge_parameters: self.merge_parameters,
            walkers,
            parameters: self.parameters,
            calling_loop_config: self.calling_loop_config,
            candidate_selection: self.candidate_selection,
            assembly_check: self.assembly_check,
        })
    }

    /// **Read this run's cohort into cohort loci, in genome order, on one thread.**
    ///
    /// Every sample's walker advances at the merge frontier, in one place; the merge draws each
    /// forward only as far as the ground it is building, so each stretch of each file is decoded
    /// once and no walker runs ahead of what has been asked for (spec §5.1).
    ///
    /// **This is the merge's oracle, not the path a run takes.** It proves the join — real
    /// reads through real walkers into the merge that was until now fed from memory — and it
    /// stops at the evidence, so a test can compare cohort loci without a genotyper's answers
    /// in the way. [`call_cohort`](Self::call_cohort) is what a run drives, and it is this
    /// same driver with the calls made where each locus is built.
    ///
    /// **Its memory is the whole cohort's surviving loci**, which no real run can afford
    /// (spec §5.1's `callers in flight × one cohort locus` plus the frontier). That is what an
    /// oracle wants and what `call_cohort` no longer does.
    ///
    /// **It keeps none of what the walk counted**, and that is deliberate rather than a gap:
    /// this is the merge's own oracle, not the path a run takes.
    /// [`call_cohort`](Self::call_cohort) is the run's, and it hands the tallies back.
    pub fn merge_cohort(self) -> Result<RegionOutcome, RunError> {
        let pieces = self.walkers()?;
        let merge = pieces.merge_parameters;
        let segmentation = Arc::clone(&pieces.segmentation);
        let mut cache = ObservationCache::over(pieces.walkers);
        merge_cohort_through_cache(
            segmentation.analysed_regions(),
            &mut cache,
            merge.cohort_locus_builder_regions_len,
            merge.max_cohort_locus_span,
            merge.min_alt_reads,
        )
    }

    /// **Call this run's cohort: reads in, one called locus per locus the merge kept, in
    /// genome order.**
    ///
    /// **A locus is not a position** — a deletion joins consecutive positions into one, so how
    /// many loci a stretch of genome yields is the merge's answer rather than its length
    /// (arch §3.1).
    ///
    /// The whole of direct mode's join in one call — every sample's files walked at the merge
    /// frontier, the cohort's loci assembled over the analysed ground, and each locus
    /// genotyped **where it is built**, before the next one is closed.
    ///
    /// # Calling happens inside the builder, and the spec says it may
    ///
    /// A call reads nothing outside its own locus, so where it runs relative to the merge is
    /// free (`run_streaming.md` §3.1, `cohort_merge.md` §6.3). Calling inside the builder is
    /// what lets the cohort observation be dropped as soon as its genotypes exist: a cohort
    /// observation carries every covering sample's reads folded onto the locus's alleles,
    /// where a called locus carries one genotype per sample. The proof that it is free is
    /// `calling_inside_the_builder_gives_what_calling_after_the_merge_gives`, which calls the
    /// same cohort both ways and refuses any difference.
    ///
    /// # What it does not do
    ///
    /// **It still returns everything at once**, and **its cover is the serial one** — this is
    /// the oracle the parallel-covered record path is compared against (Milestone E), so it
    /// deliberately takes no part in the parallel cover: one driver has to keep the schedule
    /// the fixtures were reasoned about under. Spec §5.1 bounds a run at `callers in flight ×
    /// one cohort locus`; what this bounds is the *observations*, not the calls.
    /// **`loci_too_wide_to_assemble` accumulates for the whole run too**, and nothing
    /// bounds it either.
    ///
    /// **Every locus goes down the generic path.** Repeat-tract candidate selection is
    /// specified and unbuilt, and both tract generator slots are refused as such — so a run
    /// over ground with tracts in it is short rather than wrong, and
    /// [`CohortWalkTallies`] says by how much.
    ///
    /// # Errors
    ///
    /// The first sample whose walk fails ends the run, naming the sample and how far it got
    /// ([`RunError::SourceFailed`]). Calling itself cannot fail: a locus whose loop did not
    /// settle comes back with `converged` false rather than as an error.
    pub fn call_cohort<S, G>(self, genotyper: &G) -> Result<CalledCohort, RunError>
    where
        G: LocusGenotyper<S>,
        S: Default,
    {
        let run_sample_count = self.samples.len();
        let sample_names: Vec<String> = self.sample_names().map(str::to_owned).collect();
        let pieces = self.walkers()?;
        let RunReadyToWalk {
            segmentation,
            merge_parameters,
            walkers,
            parameters,
            calling_loop_config,
            candidate_selection,
            assembly_check,
        } = pieces;

        let frozen = parameters.view();
        // **Per worker, not per locus** — there is one worker here. The shaping buffers are
        // cleared and refilled at every locus (`evidence_shaping`'s module note); the calling
        // scratch is resized and refilled by `CallingScratch::prepare_for_locus`, which
        // `call_locus` runs before it reads anything.
        let mut shaping = GenericEvidenceScratch::default();
        let mut tract_shaping = SsrEvidenceScratch::default();
        let tract_selection = SsrSelectionConfig::at_ploidy(frozen.ploidy());
        let mut scratch: CallingScratch<S> = CallingScratch::default();
        let mut called_loci = Vec::new();
        let mut loci_too_wide_to_assemble = Vec::new();
        let mut loci_with_nobody_to_call = Vec::new();
        // Counted rather than collected: a run over tract-rich ground meets millions of these,
        // and what the report owes is how many, not where each one was.
        let mut tracts = TractOutcomes::default();

        let mut cache = ObservationCache::over(walkers);
        merge_cohort_handing_each_locus_over(
            segmentation.analysed_regions(),
            &mut cache,
            merge_parameters.cohort_locus_builder_regions_len,
            merge_parameters.max_cohort_locus_span,
            merge_parameters.min_alt_reads,
            &mut |observation| {
                let region = observation.region;
                match call_one_cohort_locus(
                    genotyper,
                    &observation,
                    &frozen,
                    &candidate_selection,
                    &tract_selection,
                    &calling_loop_config,
                    run_sample_count,
                    &mut shaping,
                    &mut tract_shaping,
                    &mut tracts,
                    &mut scratch,
                    // This entry point wants the call and nothing beside it, so the observation's
                    // remapping and leftover are dropped where they were built.
                    |inference, _remap, _unmatched, _verdict| inference,
                ) {
                    LocusOutcome::Called(called) => called_loci.push(called),
                    LocusOutcome::NobodyToCall => loci_with_nobody_to_call.push(region),
                    // Both counted inside the dispatch, where the verdict that decided them is.
                    LocusOutcome::BundleSetAside | LocusOutcome::TractWithoutWholeRepeats => {}
                }
            },
            &mut loci_too_wide_to_assemble,
        )?;

        Ok(CalledCohort {
            called_loci,
            loci_too_wide_to_assemble,
            loci_with_nobody_to_call,
            tracts,
            walk: CohortWalkTallies::of(sample_names, cache.into_sources(), assembly_check),
        })
    }

    /// **Call this run's cohort and hand every record over as it is finished, keeping none.**
    ///
    /// The path a command takes, where [`call_cohort`](Self::call_cohort) is the oracle:
    /// identical calling, and the answers leave one at a time instead of accumulating. What a
    /// run holds is then its open files, the merge's frontier and one record — spec §5.1's
    /// bound.
    ///
    /// # Where this run's parallelism is (Milestone E, decided 2026-09-01)
    ///
    /// **Each cover's samples are drawn forward concurrently; everything else is one
    /// thread.** Measured on 63 tomato accessions, decoding reads is 88% of `call_cohort`
    /// and genotyping 5–6%, so the parallelism went to the decode
    /// ([`merge_cohort_handing_each_locus_over_covering_samples_in_parallel`]'s own note
    /// carries the numbers) and no pool of genotyping workers was built. The output is
    /// identical at every thread count — the cover reaches the same fixpoint by any
    /// schedule, and assembly and calling stay on this thread in genome order — which is
    /// spec §12.2's oracle: pinned at the driver by
    /// `the_parallel_cover_gives_the_serial_drivers_answer`, and to be pinned end to end by
    /// Milestone E2's concurrency-invariance fixture, which is not built yet.
    ///
    /// `hand_over` is given each record in genome order and may refuse. **The first refusal is
    /// the run's answer**, naming the locus it was writing, and no record is handed over after
    /// it — but the merge is not stopped, because its sink cannot say so (see the note beside
    /// `stopped` in the body). So the walk finishes the analysed ground before the error is
    /// returned.
    ///
    /// # Not every called locus becomes a record
    ///
    /// **A locus no written genotype carries an alternative at is not written** — spec §9's
    /// rule, which is the whole of why this file has no gVCF and no reference blocks: the
    /// record's absence says *nothing here*. The count of those is in the answer, because
    /// *called* and *written* differing is a fact about a run and not an accident of it.
    ///
    /// # The reference is read once more, and only where a record needs it
    ///
    /// A record with an empty allele — an insertion's or a deletion's nature — is written by
    /// padding every allele with the reference base beside the span, which the locus does not
    /// carry. So this holds one reference accessor of its own for the whole run, over the index
    /// and contig table the walkers already share, and reads one base per such record
    /// ([`padding_base_beside`]). It is one more open file, and what it spends is the slack
    /// inside [`DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES`] rather than a raise.
    ///
    /// **On today's path that base is never fetched**, because the generic mint anchors its
    /// indels and so no generic allele is ever empty — see [`records`](super::records)'s own
    /// note. The accessor is still opened, and the answer is still computed where a record needs
    /// one, because [`VcfRecord::new`] asserts a padding base is carried exactly when some
    /// allele is empty.
    ///
    /// # Errors
    ///
    /// The first sample whose walk fails ends the run ([`RunError::SourceFailed`]); so does the
    /// first record `hand_over` refuses ([`RunError::RecordNotWritten`]) and a padding base the
    /// reference will not serve ([`RunError::PaddingBaseUnreadable`]). Calling itself cannot
    /// fail: a locus whose loop did not settle comes back with `converged` false and is written
    /// on the `EMNoConv` filter.
    pub fn call_cohort_handing_each_record_over<S, G, E>(
        self,
        genotyper: &G,
        hand_over: &mut impl FnMut(&VcfRecord) -> Result<(), E>,
    ) -> Result<WrittenCohort, RunError>
    where
        G: LocusGenotyper<S>,
        S: Default,
        E: std::error::Error + Send + Sync + 'static,
    {
        let run_sample_count = self.samples.len();
        let sample_names: Vec<String> = self.sample_names().map(str::to_owned).collect();
        // **Minted before the walkers take the run apart**, because `walkers` consumes the
        // caller. One accessor for the whole run, never shared: it walks forward with the merge
        // and releases what it has passed.
        let padding_reference = self.walk_reference.accessor();
        let pieces = self.walkers()?;
        let RunReadyToWalk {
            segmentation,
            merge_parameters,
            walkers,
            parameters,
            calling_loop_config,
            candidate_selection,
            assembly_check,
        } = pieces;

        let frozen = parameters.view();
        let mut shaping = GenericEvidenceScratch::default();
        let mut tract_shaping = SsrEvidenceScratch::default();
        let tract_selection = SsrSelectionConfig::at_ploidy(frozen.ploidy());
        let mut scratch: CallingScratch<S> = CallingScratch::default();
        let mut padding_scratch = Vec::new();
        let mut records_written = 0_u64;
        let mut loci_called_but_not_written = 0_u64;
        let mut loci_too_wide_to_assemble = Vec::new();
        let mut loci_with_nobody_to_call = Vec::new();
        // Counted rather than collected — see the other driver.
        let mut tracts = TractOutcomes::default();
        // **The merge's sink cannot fail**, so the first failure is stashed and the sink does
        // nothing after it. Reported ahead of whatever the merge itself then returns, because it
        // happened first.
        //
        // **⚑ The walk does not stop, and that is a cost worth naming rather than a choice.**
        // `merge_cohort_handing_each_locus_over` takes a sink that returns nothing
        // (`cohort_merge/serial.rs`), so nothing this side can do ends the merge — it runs to the
        // end of the analysed ground, decoding every remaining read of every sample, before the
        // error is returned. On a fixture that is invisible; on a cohort whose disk fills at the
        // first chromosome it is the rest of the genome decoded for nothing. Fixing it means the
        // merge's sink saying *stop* — one `ControlFlow` through both drivers and the region
        // builder — which is a change to the merge's interface and not this step's.
        let mut stopped: Option<RunError> = None;

        let mut cache = ObservationCache::over(walkers);
        let merged = merge_cohort_handing_each_locus_over_covering_samples_in_parallel(
            segmentation.analysed_regions(),
            &mut cache,
            merge_parameters.cohort_locus_builder_regions_len,
            merge_parameters.max_cohort_locus_span,
            merge_parameters.min_alt_reads,
            &mut |observation| {
                if stopped.is_some() {
                    return;
                }
                let region = observation.region;
                let built = call_one_cohort_locus(
                    genotyper,
                    &observation,
                    &frozen,
                    &candidate_selection,
                    &tract_selection,
                    &calling_loop_config,
                    run_sample_count,
                    &mut shaping,
                    &mut tract_shaping,
                    &mut tracts,
                    &mut scratch,
                    |inference, remap, unmatched, verdict| {
                        // **Asked before the reference is read**, so a locus that establishes
                        // no variant costs no fetch and no evidence gathering.
                        if !a_written_genotype_carries_an_alternative(&inference) {
                            return Ok(None);
                        }
                        let alleles = inference.alleles();
                        let padding = padding_base_beside(
                            &padding_reference,
                            inference.region,
                            alleles,
                            &mut padding_scratch,
                        )
                        .map_err(|source| {
                            RunError::PaddingBaseUnreadable {
                                locus: inference.region,
                                source,
                            }
                        })?;
                        let evidence = evidence_for_output(
                            &inference,
                            &observation,
                            remap,
                            unmatched,
                            verdict,
                            padding,
                        );
                        Ok(Some(assemble_record(&inference, evidence)))
                    },
                );
                match built {
                    LocusOutcome::NobodyToCall => loci_with_nobody_to_call.push(region),
                    // Both counted inside the dispatch, where the verdict that decided them is.
                    LocusOutcome::BundleSetAside | LocusOutcome::TractWithoutWholeRepeats => {}
                    LocusOutcome::Called(Err(error)) => stopped = Some(error),
                    LocusOutcome::Called(Ok(None)) => loci_called_but_not_written += 1,
                    LocusOutcome::Called(Ok(Some(record))) => match hand_over(&record) {
                        Ok(()) => records_written += 1,
                        Err(source) => {
                            stopped = Some(RunError::RecordNotWritten {
                                locus: region,
                                source: Box::new(source),
                            });
                        }
                    },
                }
            },
            &mut loci_too_wide_to_assemble,
        );
        if let Some(stopped) = stopped {
            return Err(stopped);
        }
        merged?;

        Ok(WrittenCohort {
            records_written,
            loci_called_but_not_written,
            loci_too_wide_to_assemble,
            loci_with_nobody_to_call,
            tracts,
            walk: CohortWalkTallies::of(sample_names, cache.into_sources(), assembly_check),
        })
    }
}

/// **What became of this run's repeat tracts** — the partition the run report prints, the way
/// the SNP/indel path's own outcomes already partition.
///
/// **The five are disjoint and they sum to every tract-kind locus the merge built.** Two of them
/// were never scored; three were, and two of those three carry a `FILTER` the file states. The
/// partition is what makes the report checkable: a run whose tract ground is charged to *called*
/// while a third of its tracts were refused would be saying something false, and the sum is what
/// says it is not.
///
/// **The refusals have to be counted because the file cannot state them.** A tract refused as
/// not periodic is called over the reference tract alone, so every sample is homozygous
/// reference and no record is written — in the file it is indistinguishable from a tract nobody
/// varied at (`doc/devel/ng/spec/vcf_output.md` §9). The count is the only place it appears.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct TractOutcomes {
    /// Scored, and carrying no repeat-tract filter.
    pub called: u64,
    /// Scored, and refused as not varying in whole motif units — `notPeriodic`.
    pub not_periodic: u64,
    /// Scored, and carrying more candidate sequences than the cap admits — `tooManyAlleles`.
    /// The locus is still called over the ones the cap kept.
    pub too_many_alleles: u64,
    /// **Not scored: a candidate carrying no whole motif copy**, which the stutter ladder has no
    /// rung for. See [`repeat_counts_the_tract_model_can_take`] for how often this fires.
    pub without_whole_repeats: u64,
    /// **Not scored: a repeat cluster with no clean flanks**, which nothing in the run builds a
    /// caller for — the bundle generator is deferred and has no model.
    pub bundles_set_aside: u64,
}

impl TractOutcomes {
    /// Every tract-kind locus the merge built — the sum of the five.
    #[must_use]
    pub fn built(&self) -> u64 {
        self.called
            + self.not_periodic
            + self.too_many_alleles
            + self.without_whole_repeats
            + self.bundles_set_aside
    }

    /// Those that were scored and then refused by a filter of their own.
    #[must_use]
    pub fn refused_by_a_filter(&self) -> u64 {
        self.not_periodic + self.too_many_alleles
    }
}

/// **What became of one cohort locus** — the four ends a driver has to tell apart.
///
/// Three of them produce no record and are not the same fact: a locus nobody can be called at,
/// a repeat cluster nothing builds a caller for, and a repeat tract the tract model cannot
/// describe are three different things to report (`calling_loop_ssr.md` §3.2).
enum LocusOutcome<R> {
    /// Called, and this is what the caller's own `finish` made of the answer.
    Called(R),
    /// **No sample of the run can be called here** — every one of them lost the alleles it
    /// earned, so there are no rows to genotype. Counted, never fatal (owner's ruling,
    /// 2026-09-01).
    NobodyToCall,
    /// **A repeat cluster with no clean flanks**, which nothing in the run builds a caller for
    /// yet: the bundle generator is deferred (`candidate_alleles_ssr.md` §11) and the tract
    /// model wants a single tract with two flanks. Set aside and counted.
    BundleSetAside,
    /// **A repeat tract holding a candidate that is not a whole number of motif copies** —
    /// a sequence shorter than one copy of the unit, which the stutter model is not written
    /// on: its ladder is in whole repeats and this candidate sits below the bottom rung.
    ///
    /// **Rare, and measured rather than assumed**: over HG002's 50,000-region Tier set through
    /// ng's own catalog, 1 of 17,315 kept candidates at 30× and none at all at 50× or 300×
    /// (`ng_ssr_selection_e2_2026-09-02.md`'s runs). Refused and counted rather than scored
    /// under a model that does not apply, and rather than widening the evidence's repeat count
    /// to admit zero — which is the change to make if that count ever stops being one in ten
    /// thousand.
    TractWithoutWholeRepeats,
}

/// **One cohort locus, down whichever path its kind names** — the dispatch the STR calling
/// loop's design asks for (`calling_loop_ssr.md` §3.2).
///
/// The two paths share their genotyper, their parameters, their calling-loop configuration and
/// their scratch; what differs is which selector narrows the locus and which shaper builds its
/// evidence. **The kind is read off the observation and never inferred from which fields are
/// populated**, so a third kind has to be decided here rather than falling into whichever arm
/// the compiler picks for it.
///
/// **A repeat tract's reads are not narrowed to its candidates**, which is the one shape
/// difference worth naming beside the two calls. The SNP/indel path hands the loop each
/// sample's reads folded onto the surviving alleles; the tract path hands it every observation
/// the sample showed, because the stutter model scores a read against a candidate rather than
/// matching it to one (`doc/devel/ng/spec/read_likelihoods.md` §4). Selection still decides the
/// candidate list, and the leftover it produces is still what says a sample must be emitted as
/// missing.
#[expect(
    clippy::too_many_arguments,
    reason = "the same nine `call_one_generic_locus` takes, plus the tract path's own selection \
              configuration and shaping buffers — a struct grouping them would exist for this \
              signature alone, and both paths' arguments are already the run's own values"
)]
fn call_one_cohort_locus<S, G, R>(
    genotyper: &G,
    observation: &CohortObservation,
    parameters: &FrozenParameters<'_>,
    candidate_selection: &CandidateSelectionConfig,
    tract_selection: &SsrSelectionConfig,
    calling_loop_config: &RunnableCallingLoopConfig,
    run_sample_count: usize,
    shaping: &mut GenericEvidenceScratch,
    tract_shaping: &mut SsrEvidenceScratch,
    tracts: &mut TractOutcomes,
    scratch: &mut CallingScratch<S>,
    finish: impl FnOnce(LocusInference, &AlleleRemap, &[UnmatchedSupport], SelectionVerdict) -> R,
) -> LocusOutcome<R>
where
    G: LocusGenotyper<S>,
{
    match &observation.kind {
        LocusKind::Generic => match call_one_generic_locus(
            genotyper,
            observation,
            parameters,
            candidate_selection,
            calling_loop_config,
            run_sample_count,
            shaping,
            scratch,
            finish,
        ) {
            Some(called) => LocusOutcome::Called(called),
            None => LocusOutcome::NobodyToCall,
        },
        LocusKind::Ssr(detail) => call_one_ssr_locus(
            genotyper,
            observation,
            detail,
            parameters,
            tract_selection,
            calling_loop_config,
            run_sample_count,
            tract_shaping,
            tracts,
            scratch,
            finish,
        ),
        LocusKind::SsrBundle => {
            tracts.bundles_set_aside += 1;
            LocusOutcome::BundleSetAside
        }
    }
}

/// **Which of the report's tract outcomes this locus is**, from selection's own verdict.
///
/// **Counted here rather than read off the written record**, because two of the three leave no
/// record to read. A tract refused as `notPeriodic` is called over the reference tract alone, so
/// every sample is homozygous reference and the locus establishes no variant — it is left out of
/// the file, where it is indistinguishable from a tract nobody varied at
/// (`doc/devel/ng/spec/vcf_output.md` §9). This count is the only place it appears.
///
/// **A truncated tract is still called** over the sequences the cap kept, and it is counted as
/// `tooManyAlleles` rather than as called so that the two are not summed into one number a
/// reader would take for clean calls.
fn count_this_tract(verdict: SelectionVerdict, tracts: &mut TractOutcomes) {
    match verdict {
        SelectionVerdict::NotPeriodic => tracts.not_periodic += 1,
        SelectionVerdict::Truncated { .. } => tracts.too_many_alleles += 1,
        _ => tracts.called += 1,
    }
}

/// **Selection's repeat counts as the tract's evidence takes them, or nothing** — the one place
/// a candidate carrying no whole motif copy stops the locus.
///
/// Selection counts a candidate's whole repeats by flooring its length by the motif's period,
/// so a sequence shorter than one copy of the unit comes back as zero; the evidence's counts are
/// [`NonZeroU32`], because the stutter ladder is written in whole repeats and a candidate below
/// its bottom rung has no rung. The conversion is therefore the refusal, and it is written once
/// here rather than as a `map` inside the tract arm, so that its own test can reach it.
///
/// **Its frequency is measured, not assumed**: over HG002's 50,000-region Tier set through ng's
/// own repeat catalog, 1 of 17,315 kept candidates at 30× carries no whole repeat, and none at
/// all at 50× or 300×
/// (`doc/devel/reports/implementations/ng_ssr_selection_e2_2026-09-02.md`'s runs). If that ever
/// stops being one in ten thousand, the change to make is to the evidence's own count type
/// rather than to this refusal.
fn repeat_counts_the_tract_model_can_take(counts: &[u32]) -> Option<Vec<NonZeroU32>> {
    counts.iter().copied().map(NonZeroU32::new).collect()
}

/// One repeat tract from evidence to genotypes — [`call_one_cohort_locus`]'s tract arm.
///
/// Three calls, the same three the SNP/indel path makes: narrow the locus
/// ([`select_ssr`]), shape its evidence ([`shape_ssr_locus`]), and genotype it with the run's
/// own genotyper. The differences are all in the first two.
///
/// # Why there is no "nobody to call" here
///
/// The SNP/indel path can leave a sample with no rows at all — the cap removes an allele the
/// sample earned and it is set aside for the rest of the locus — and a locus where that
/// happens to every sample has nothing to genotype. **A tract sets no sample aside**: a
/// discovery round can put back a length the cap cut, so nobody is locked out
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1), and a sample that showed nothing carries an
/// empty observation list whose sum is zero, which is the right answer rather than a special
/// case. So the tract arm's only refusals are the two the outcome names.
#[expect(
    clippy::too_many_arguments,
    reason = "the tract path's own spelling of `call_one_generic_locus`'s nine — see there"
)]
fn call_one_ssr_locus<S, G, R>(
    genotyper: &G,
    observation: &CohortObservation,
    detail: &SsrDetail,
    parameters: &FrozenParameters<'_>,
    tract_selection: &SsrSelectionConfig,
    calling_loop_config: &RunnableCallingLoopConfig,
    run_sample_count: usize,
    tract_shaping: &mut SsrEvidenceScratch,
    tracts: &mut TractOutcomes,
    scratch: &mut CallingScratch<S>,
    finish: impl FnOnce(LocusInference, &AlleleRemap, &[UnmatchedSupport], SelectionVerdict) -> R,
) -> LocusOutcome<R>
where
    G: LocusGenotyper<S>,
{
    let narrowed = select_ssr(
        observation,
        tract_selection,
        scratch.candidate_selection_mut(),
    );
    // **Converted before anything is shaped**, so a candidate the tract model cannot describe
    // stops the locus rather than reaching the emission as a one-repeat allele.
    let Some(repeat_counts) = repeat_counts_the_tract_model_can_take(&narrowed.repeat_counts)
    else {
        tracts.without_whole_repeats += 1;
        return LocusOutcome::TractWithoutWholeRepeats;
    };
    let SsrLocusSelection { selection, .. } = narrowed;
    let (alleles, verdict, unmatched, remap) = selection.into_parts();
    count_this_tract(verdict, tracts);

    let observations_of_each_run_sample = tract_shaping.rebuild(observation, run_sample_count);
    // **Two per-locus allocations, and both are the borrow checker's price rather than a
    // choice.** The slice-of-slices cannot live in the scratch beside the buffers it borrows,
    // and `views` cannot outlive the locus for the reason `shape_generic_locus` documents: a
    // `Vec` is invariant in its element type, so one held across two loci would hold the
    // first's borrow open into the second.
    let per_sample: Vec<&[SequenceObservation]> = observations_of_each_run_sample
        .iter()
        .map(Vec::as_slice)
        .collect();
    let mut views = Vec::new();
    let evidence = shape_ssr_locus(
        observation.region,
        &per_sample,
        detail,
        &repeat_counts,
        &mut views,
    );
    let inference =
        genotyper.call_locus(&evidence, parameters, alleles, calling_loop_config, scratch);
    LocusOutcome::Called(finish(inference, &remap, &unmatched, verdict))
}

/// One cohort locus from evidence to genotypes: **which alleles it is called over, whose reads
/// say what, and what each sample's genotype is.**
///
/// Three existing calls and nothing of its own (arch §3.2). It is a free function rather than
/// a method because it belongs to no run in particular — and if a pool of genotyping workers
/// is ever built (measured out at 63 samples, conditional on a large-cohort measurement), a
/// worker runs it over a scratch that worker owns.
///
/// **`views` is declared here and nowhere higher up, and the type says why.** It borrows the
/// rows `shaping` holds, and a `Vec` is invariant in its element type, so a list kept across
/// two loci would hold the first locus's borrow open into the second — an `E0499`
/// ([`shape_generic_locus`]). It is the shaping step's one per-locus allocation; the call
/// itself makes more, since a [`LocusInference`] owns its per-sample calls and its expected
/// allele copies.
///
/// # `finish` is what the caller does with the call *before the locus is gone*
///
/// **The call is not everything a run wants from a locus.** A record also needs what each
/// sample's reads showed and which of them no written allele explains, and those live in the
/// cohort observation and in candidate selection's leftover — both of which end here
/// ([`records`](super::records)). So the answer is handed to a closure rather than returned:
/// inside it the observation, the allele remapping and the leftover are all still in scope, and
/// after it they are dropped. [`call_cohort`](AlignedFilesVariantCaller::call_cohort), which
/// wants the call alone, passes a closure that returns it unchanged.
#[expect(
    clippy::too_many_arguments,
    reason = "the nine are the five things a locus is called against — the model, the \
              evidence, the parameters, and the selection's and the loop's configurations — \
              plus the run's sample count, the two scratches and what to do with the answer. \
              Grouping any of them would make a struct whose only \
              purpose is this signature — and a genotyping pool, if a large-cohort \
              measurement ever justifies one, hands the same nine to a worker."
)]
fn call_one_generic_locus<S, G, R>(
    genotyper: &G,
    observation: &CohortObservation,
    parameters: &FrozenParameters<'_>,
    candidate_selection: &CandidateSelectionConfig,
    calling_loop_config: &RunnableCallingLoopConfig,
    run_sample_count: usize,
    shaping: &mut GenericEvidenceScratch,
    scratch: &mut CallingScratch<S>,
    finish: impl FnOnce(LocusInference, &AlleleRemap, &[UnmatchedSupport], SelectionVerdict) -> R,
) -> Option<R>
where
    G: LocusGenotyper<S>,
{
    // Selection's buffers live on the calling scratch and are sized by selection itself, so
    // this borrow ends before the locus's shape is known and the one below begins.
    let selection = select_generic(
        observation,
        candidate_selection,
        scratch.candidate_selection_mut(),
    );
    let mut views = Vec::new();
    let evidence = shape_generic_locus(
        shaping,
        observation,
        &selection,
        run_sample_count,
        &mut views,
    );
    // The allele table leaves the selection by value: a discovery round appends to it and the
    // final prune shrinks it, so the loop owns it and hands it back inside the inference.
    // **A locus no sample of the run can be called at is counted, not called, and never
    // fatal** (owner's ruling, 2026-09-01). The genotyper's precondition is that there is
    // somebody to call — its scratch cannot be prepared for no rows — so the question is asked
    // here rather than left to fail there.
    if evidence.callable_sample_count() == 0 {
        return None;
    }
    // The allele table leaves the selection by value: a discovery round appends to it and the
    // final prune shrinks it, so the loop owns it and hands it back inside the inference.
    // **The other two parts stay here for `finish`**: the remapping says which of the merge's
    // alleles the loop was given, and the leftover says what each covering sample showed that it
    // was not — neither of which the inference carries and neither of which outlives this call.
    let (alleles, verdict, unmatched, remap) = selection.into_parts();
    let inference =
        genotyper.call_locus(&evidence, parameters, alleles, calling_loop_config, scratch);
    Some(finish(inference, &remap, &unmatched, verdict))
}

/// **A run whose walkers are built and which has not yet read a byte** — one walker per
/// sample, plus everything each locus will be called against.
///
/// **Private, and a struct rather than a tuple**: seven positional fields, four of which are
/// this run's configuration, is a signature a reader has to count their way through. The run's
/// two entry points destructure it and take what each needs.
struct RunReadyToWalk {
    segmentation: Arc<Segmentation>,
    merge_parameters: MergeParameters,
    /// One per sample, in the run's sample order.
    walkers: Vec<RunWalker>,
    parameters: RunParameters,
    calling_loop_config: RunnableCallingLoopConfig,
    candidate_selection: CandidateSelectionConfig,
    assembly_check: AssemblyCheckOutcome,
}

/// What a run's walker is, spelled once: the alignment-file source over the run's own segments.
///
/// A type alias rather than a wrapper, because it adds nothing — it is
/// [`AlignmentFilesWalker`] with its region stream named. It buys the one signature that names
/// it the right words, and it is the name any later concurrency work will want.
pub type RunWalker = AlignmentFilesWalker<RunSegments>;

/// **What a calling run produced: the genotypes, the two kinds of ground it produced no
/// genotypes for, and what the walk counted on the way.**
///
/// The four are one value because they are one run, and a report that quoted the genotypes
/// alone would be describing a cohort rather than a run. **`called_loci` is not every locus the
/// merge assembled**: add [`Self::loci_with_nobody_to_call`] for that. And neither of the two
/// refusal lists is the same fact as a locus nobody varied at, which is counted nowhere by
/// design (`cohort_merge.md` §3.3).
#[derive(Debug)]
pub struct CalledCohort {
    /// One per surviving locus, in genome order.
    pub called_loci: Vec<LocusInference>,
    /// **The ground of the loci the merge declined to assemble for being wider than
    /// `max_cohort_locus_span`**, in genome order.
    ///
    /// **Not a failure, and the name says so** — the merge's own field for these is called
    /// `failed_locus_spans`, and a run report that echoed that word would tell an operator
    /// their caller failed N times when what it did was decline ground it was configured not
    /// to assemble.
    ///
    /// Nor is it the ground nobody varied at: a locus too quiet to build is ground the caller
    /// examined and found matching the reference. Only this one is a setting worth reporting
    /// (`cohort_merge.md` §3.3).
    pub loci_too_wide_to_assemble: Vec<GenomeRegion>,
    /// **The ground of the loci no sample of the run could be called at**, in genome order —
    /// loci the merge assembled and that produced no genotypes.
    ///
    /// The allele cap cuts a sequence rather than refusing a locus
    /// (`doc/devel/ng/spec/candidate_alleles.md` §4.1), and a sample that had reads on the cut
    /// sequence is ruled uncallable. This is the case where **no sample of the run** is left —
    /// which needs every sample to have covered the locus, since one that covered nothing is
    /// callable and scored by the prior alone
    /// ([`LocusEvidence::callable_sample_count`](crate::ng::calling::LocusEvidence::callable_sample_count)).
    ///
    /// **It is not an error** — one hard locus must not end a cohort's run (owner's ruling,
    /// 2026-09-01). **And it is a third fact**: not a locus the width bound refused, which was
    /// never assembled; not a locus nobody varied at, which is counted nowhere; not a sample
    /// set aside at a locus other samples were called at, which is
    /// [`SampleGenotypeCall::Missing`](crate::ng::calling::SampleGenotypeCall::Missing).
    ///
    /// **A non-empty list is worth acting on**: raising `max_candidate_alleles` keeps more of
    /// what those loci vary over.
    pub loci_with_nobody_to_call: Vec<GenomeRegion>,
    /// **What became of this run's repeat tracts**, partitioned five ways — called, refused by
    /// each of the two tract filters, and the two kinds of tract nothing scored.
    ///
    /// **Not the same fact as the walk's `unhandled_not_implemented`**, which counts *regions*
    /// whose generator slot is unfilled. This counts *loci* that were built and merged across
    /// the cohort.
    pub tracts: TractOutcomes,
    /// What each sample's walk saw, and what the run could check about the assembly.
    pub walk: CohortWalkTallies,
}

/// **What a run that wrote its calls out produced** — everything [`CalledCohort`] carries
/// except the loci themselves, which such a run hands over one at a time rather than keeping.
///
/// **Written and called are different counts, and the difference is a fact about the run.** A
/// locus no written genotype carries an alternative at establishes no variant and is left out
/// of the file (`doc/devel/ng/spec/vcf_output.md` §9); there is no gVCF and no reference block,
/// so its absence is the file saying *nothing here*. A run whose two counts are far apart
/// called a great deal of ground that came back matching the reference, which is ordinary at
/// low depth and worth being able to see.
#[derive(Debug)]
pub struct WrittenCohort {
    /// Records handed over, which is the file's record count.
    pub records_written: u64,
    /// **Loci called where no written genotype carried an alternative**, and so left out of the
    /// file. Add this to [`Self::records_written`] for the loci that were called.
    pub loci_called_but_not_written: u64,
    /// **The ground of the loci the merge declined to assemble for being wider than
    /// `max_cohort_locus_span`**, in genome order — [`CalledCohort::loci_too_wide_to_assemble`].
    pub loci_too_wide_to_assemble: Vec<GenomeRegion>,
    /// **The ground of the loci no sample of the run could be called at**, in genome order —
    /// [`CalledCohort::loci_with_nobody_to_call`].
    pub loci_with_nobody_to_call: Vec<GenomeRegion>,
    /// **What became of this run's repeat tracts** — [`CalledCohort::tracts`].
    pub tracts: TractOutcomes,
    /// What each sample's walk saw, and what the run could check about the assembly.
    pub walk: CohortWalkTallies,
}

impl WrittenCohort {
    /// How many loci were called — written, plus those that established no variant.
    #[inline]
    #[must_use]
    pub fn loci_called(&self) -> u64 {
        self.records_written
            .saturating_add(self.loci_called_but_not_written)
    }
}

/// **What every sample's walk counted, kept past the walk.**
///
/// **The run does not own its walkers once the merge has them** — the observation cache does,
/// for the merge's whole duration, and hands them back spent
/// ([`ObservationCache::into_sources`]) — so these are copied out at that point rather than
/// read from a walker later. What is here is what a run report has to be able to state: for
/// each sample, how much of the analysed ground its walk handled, how much it could not and
/// why, and what the SNP/indel generator counted while doing it.
///
/// **⛦ The per-read-group read-filter tallies are here now (2026-09-01), and what they needed
/// was the generator rather than an accessor.** They belong to a cursor from the moment it is
/// made, and a cursor is rebuilt at every chromosome change — so a walk had already lost every
/// contig but its last, and reading them off the live cursor would have reported one
/// chromosome's drops as the run's. The generator now takes a retiring cursor's counts at the
/// boundary, the way it already took the aggregate ones, and sums the live cursor in when asked.
/// Spec §8's finish-time tally, and the failure it names if it is skipped: drop rates
/// under-report "silently, since every number stays plausible".
///
/// **One thing a run report also wants is still not here, and it is not an oversight.** The
/// repeat-tract slots' counts are absent because both slots are unfilled — a tract's ground is
/// charged to `unhandled_not_implemented` in [`SampleWalkTallies::regions`], which is where a
/// reader looks to see how short this caller's coverage is.
#[derive(Debug)]
pub struct CohortWalkTallies {
    /// One per sample, in the run's sample order.
    pub per_sample: Vec<SampleWalkTallies>,
    /// Whether the run could check what assembly its samples were aligned to, and over how
    /// much. **Carried here because it is computed at construction and the run is consumed**,
    /// so a report assembled after the run would otherwise have nowhere to read it.
    pub assembly_check: AssemblyCheckOutcome,
}

impl CohortWalkTallies {
    /// Copy what each walker counted out of it, in the order the walkers were handed over —
    /// which is the run's sample order, since that is the order the cache holds them in.
    ///
    /// # Panics
    ///
    /// If the names and the walkers are different lengths, which would pair one sample's name
    /// with another's counts — a wrong run report rather than a crash, and the accident the
    /// run's single sample order exists to prevent.
    fn of(
        sample_names: Vec<String>,
        walkers: Vec<RunWalker>,
        assembly_check: AssemblyCheckOutcome,
    ) -> Self {
        assert_eq!(
            sample_names.len(),
            walkers.len(),
            "the run holds {} samples and {} per-sample readers came back, so a report built \
             from these would put one sample's counts under another sample's name. This is a \
             defect in ng rather than anything about the data.",
            sample_names.len(),
            walkers.len()
        );
        Self {
            per_sample: sample_names
                .into_iter()
                .zip(walkers)
                .map(|(sample_name, walker)| SampleWalkTallies {
                    regions: walker.counts().clone(),
                    // **The slot holds a pileup generator or nothing**, which is what
                    // `generic_path_generators` builds, so the two arms below are the two
                    // real cases and the third is unreachable rather than merged in: a tally
                    // of another kind there would be a generator set built for a different
                    // path, and this reports it as counting nothing — which is what a
                    // generator set nobody mis-wired also reports. Nothing in a run can
                    // produce it; if `generic_path_generators` ever fills the slot from a
                    // caller's choice, this becomes two facts and needs two answers.
                    read_filters: walker.generators().read_filter_counts(),
                    snp_indel: match walker.generators().generic_counts() {
                        Some(GeneratorCounts::Pileup(counts)) => Some(*counts),
                        Some(_) | None => None,
                    },
                    sample_name,
                })
                .collect(),
            assembly_check,
        }
    }
}

/// **What one sample's walk counted** — its share of [`CohortWalkTallies`].
#[derive(Debug, Clone)]
pub struct SampleWalkTallies {
    /// The sample, by the name the run's read-group table gave it.
    pub sample_name: String,
    /// How this sample's share of the analysed ground was accounted for: regions in, regions
    /// handled, loci emitted, and the two kinds of region this caller produced nothing for —
    /// a gap it has not filled yet, and ground it will never call.
    ///
    /// **`regions_handled`, `unhandled_not_implemented` and `unhandled_out_of_scope` sum
    /// exactly to `regions_in`**, which is what makes "how much did this run not look at"
    /// answerable rather than an estimate.
    ///
    /// **In regions, not in bases**, for the handled share: the two unhandled classes each
    /// carry their base count and the handled one does not, so *what fraction of the genome
    /// did this run call* cannot be worked out from here. Regions differ in length by orders
    /// of magnitude, so a ratio of region counts is not that fraction.
    pub regions: LocusCounts,
    /// **Why this sample's reads were dropped, per read group, over the whole walk** — one
    /// entry a read group the sample's files declared, `None` for reads that named none.
    ///
    /// **Summed over every chromosome, which is what makes it the walk's** (2026-09-01). These
    /// counts belong to a cursor from the moment it is made, and a cursor is rebuilt at every
    /// chromosome change, so a run that read them off the live cursor reported its last
    /// chromosome's drops as the whole walk's. Spec §8 asks for the finish-time sum and names
    /// what skipping it costs: drop rates under-report, "silently, since every number stays
    /// plausible".
    pub read_filters: Vec<(Option<ReadGroupId>, ReadFilterCounts)>,
    /// What the SNP/indel generator counted while walking this sample — **`None` where it
    /// counted nothing at all**, which is a sample whose ground held no such region.
    ///
    /// Named for what it counts rather than for the code's word for it: this project calls
    /// the SNP/indel path the *generic* path, against the repeat-tract path, and *generic* to
    /// a reader of a run report means *unspecific*.
    ///
    /// # Two depth numbers, and reading one for the other is reading the opposite
    ///
    /// **`positions_short_of_cap` answers one question: did the read-hold ceiling cost this
    /// sample coverage?** Zero means no position was scored on fewer reads than
    /// `max_snp_column_depth` allows *because the ceiling had already given reads up*;
    /// `short_of_cap_deficit` says how many reads those positions were missing altogether.
    ///
    /// **It says nothing about the two per-position caps**, which is
    /// `column_depth_truncations` — the positions where `max_snp_column_depth` or
    /// `max_indel_column_depth` cut contributors the walk was holding. A run can have
    /// `positions_short_of_cap` at zero and `column_depth_truncations` in the millions: the
    /// ceiling kept every read and the caps then declined to score on all of them. So "did my
    /// depth settings shape the evidence" is answered by **both**, and by neither alone.
    pub snp_indel: Option<PileupGeneratorCounts>,
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

/// Descriptors one sample's **locus generator** holds, on top of what its files cost.
///
/// **Measured after the merge was driven over walkers, and the count changed** —
/// `examples/ng_open_cohort_descriptors.rs` again, on the same 63 tomato accessions. A generator
/// keeps two reference accessors of its own: one for the walk's REF fetches and one for the read
/// preparer. Each is a `WindowedRefSeq`, each opens a reader on the FASTA at its first fetch, and
/// each holds it for the run. They are per **sample**, not per file, so this is a second term in
/// the arithmetic rather than a correction to the first.
///
/// **This was missing and the refusal was wrong in the unsafe direction**: 63 samples over 63
/// files were budgeted 158 descriptors and the walking run held **253**, so a run could pass the
/// check and then die at `EMFILE` — the exact failure the check exists to prevent. Counted:
/// 3 before any sample opened, 4 with every file open, 130 with a cursor on each, **256** once
/// the generators' accessors were fetched through.
const DESCRIPTORS_A_SAMPLE_NEEDS_BESIDES_ITS_FILES: u64 = 2;

/// Descriptors a run needs for everything that is not per file or per sample: the three standard
/// streams, the reference and its index, the repeat catalog, the output and its index, and — from
/// 2026-09-01 — the reference accessor a writing run reads its padding bases from
/// ([`call_cohort_handing_each_record_over`](AlignedFilesVariantCaller::call_cohort_handing_each_record_over))
/// — **nine** — plus 23 of slack for whatever the runtime holds open on its own. The constant
/// did not move; what the ninth spent is slack.
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
    let samples = read_groups.read_groups_per_sample().len();
    let needed = descriptors_needed_for(alignment_files, samples);

    let Some(limit) = limit else {
        return Ok(());
    };

    if needed > limit {
        return Err(RunError::NotEnoughFileDescriptors {
            samples,
            alignment_files,
            per_file: DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS,
            per_sample: DESCRIPTORS_A_SAMPLE_NEEDS_BESIDES_ITS_FILES,
            allowance: DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES,
            needed,
            limit,
        });
    }
    Ok(())
}

/// How many descriptors a walking run needs: **two terms and an allowance**, because two of the
/// four a sample costs are per file and two are per sample.
fn descriptors_needed_for(alignment_files: usize, samples: usize) -> u64 {
    alignment_files as u64 * DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS
        + samples as u64 * DESCRIPTORS_A_SAMPLE_NEEDS_BESIDES_ITS_FILES
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
    use crate::ng::locus_generation::pileup::MAX_RECORD_SPAN_CEILING;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference_from_its_index, header, matching_contigs, named_bam,
        read_named_with_length, read_named_with_length_in_read_group,
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

    /// **Locus-generator settings no run would arrive at by default**, so a run that dropped
    /// them on the floor and built its generators with the shipped constants is visible.
    ///
    /// The knob moved is the hold ceiling, because it is the one whose default has already
    /// been wrong once in a way nothing noticed: at 4,096 it silently refused 19,725 reads on
    /// one ~130× tomato chromosome.
    pub(super) fn unusual_locus_generator_settings() -> PileupGeneratorConfig {
        PileupGeneratorConfig {
            max_active_reads: 4_096,
            ..PileupGeneratorConfig::default()
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
                locus_generator_settings: unusual_locus_generator_settings(),
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_a_dir, a) = bam_for("NA12878", "a.bam");

        let caller = open_over(std::slice::from_ref(&a), &reference).expect("opens");

        assert_eq!(caller.read_filters(), unusual_read_filters());
        assert_eq!(
            caller.locus_generator_settings(),
            unusual_locus_generator_settings()
        );
        assert_eq!(*caller.candidate_selection(), unusual_candidate_selection());
        assert_eq!(caller.merge_parameters(), unusual_merge_parameters());
        assert_eq!(*caller.calling_loop_config(), unusual_calling_loop_config());

        // And each of those really is unlike what ships, or the assertions above would hold
        // for a caller that ignored its arguments entirely.
        assert_ne!(unusual_read_filters(), ReadFilterConfig::default());
        assert_ne!(
            unusual_locus_generator_settings(),
            PileupGeneratorConfig::default()
        );
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

    /// **Locus-generator settings the walk could not use are refused at the door, before a
    /// file is opened.**
    ///
    /// The run is given both an impossible record-span ceiling *and* a BAM with no index
    /// beside it, and must report the settings — the refusal that costs nothing — rather than
    /// the open failure. That ordering is the whole of the check: a run whose depth caps are
    /// unusable would otherwise learn so at its first locus, after every file of a
    /// thousand-sample cohort had been opened.
    #[test]
    fn locus_generator_settings_the_walk_cannot_use_are_refused_before_a_file_is_opened() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_bad_dir, unindexed) = unindexed_bam_for("zeta", "zeta.bam");
        let read_groups = build_read_groups(&[unindexed]).expect("the fixture declares a group");
        let parameters = RunParameters::of_defaults(
            &read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        );

        let refused = AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &reference,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                locus_generator_settings: PileupGeneratorConfig {
                    // One past what a witness run can describe, which is the one knob whose
                    // ceiling this caller sets rather than inherits from production.
                    max_record_span: MAX_RECORD_SPAN_CEILING + 1,
                    ..PileupGeneratorConfig::default()
                },
                reference_with_checksums: reference.info(),
            },
            segmentation(),
            parameters,
            unusual_calling_loop_config(),
            unusual_candidate_selection(),
            unusual_merge_parameters(),
        )
        .expect_err("both are wrong, and the cheap one is reported");

        assert!(
            matches!(refused, RunError::LocusGeneratorSettings { .. }),
            "the settings, not the unopenable file: {refused:?}"
        );
        // **Asserted on the setting's name and on both numbers**, not on either alone: a
        // message carrying the values and not the knob leaves a reader with nothing to change,
        // and this is the one refusal of the six that a person reaches by typing a number.
        let rendered = format_error_chain(&refused);
        for expected in ["max_record_span", "65536", "65535"] {
            assert!(
                rendered.contains(expected),
                "the refusal must name {expected}: {rendered}",
            );
        }
    }

    /// The ground and the read-group table come back too.
    #[test]
    fn the_ground_and_the_read_group_table_come_back() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
    use super::tests::{
        bam_for, segmentation_built_on, unusual_locus_generator_settings, unusual_read_filters,
    };
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, fixture_reference_from_its_index, header, matching_contigs, named_bam,
        read_named_with_length,
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
                locus_generator_settings: unusual_locus_generator_settings(),
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();

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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
                assert_eq!(*needed, descriptors_needed_for(1, 1));
                assert_eq!(*limit, 4);
            }
            other => panic!("expected NotEnoughFileDescriptors, got {other:?}"),
        }

        // **The rendered message, not the fields.** This step's whole product is what a person
        // reads, so the arithmetic they are asked to check has to be in the string: what it
        // needs, what it may open, and the two numbers those come from.
        let message = error.to_string();
        assert!(
            message.contains("needs 36 open files"),
            "names what it needs: {message}",
        );
        assert!(
            message.contains("may open 4"),
            "names the limit, not just a digit of it: {message}",
        );
        assert!(
            message.contains("1 alignment files at 2 each"),
            "shows the per-file term: {message}",
        );
        assert!(
            message.contains("1 samples at 2 more each"),
            "and the per-sample term, which is the one a file count cannot stand in for: \
             {message}",
        );
        assert!(
            message.contains("and 32 for the reference"),
            "and the part that is neither: {message}",
        );
        assert!(
            message.contains("ulimit -n 36"),
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
            Some(descriptors_needed_for(1, 1)),
        )
        .expect("exactly enough is enough");
        refuse_if_more_descriptors_are_needed_than_allowed(&read_groups, None)
            .expect("no limit reported is no refusal");
    }

    /// **The count has two terms and they move independently**, because two of the four
    /// descriptors a one-file sample costs belong to the file and two to the sample.
    ///
    /// **It was one term until the merge was driven over walkers**, and that was wrong in the
    /// unsafe direction: 63 samples over 63 files were budgeted 158 and the walking run held
    /// 253, so a run could pass this check and die at `EMFILE` — the failure it exists to
    /// prevent (`examples/ng_open_cohort_descriptors.rs`).
    #[test]
    fn the_descriptor_count_grows_with_files_and_with_samples_separately() {
        assert_eq!(
            descriptors_needed_for(1, 1) + 2 * DESCRIPTORS_AN_ALIGNMENT_FILE_NEEDS,
            descriptors_needed_for(3, 1),
            "a sample sequenced across three lanes pays for three files and one walk",
        );
        assert_eq!(
            descriptors_needed_for(1, 1) + 2 * DESCRIPTORS_A_SAMPLE_NEEDS_BESIDES_ITS_FILES,
            descriptors_needed_for(1, 3),
            "three samples sharing one file pay for one file and three walks",
        );
        assert_eq!(
            descriptors_needed_for(0, 0),
            DESCRIPTORS_A_RUN_NEEDS_BESIDES_ITS_ALIGNMENT_FILES,
            "with no cohort at all, only the run's own allowance is needed",
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
        let (_open_dir, opened_against) = fixture_reference_from_its_index();
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
        let (_open_dir, opened_against) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_open_dir, opened_against) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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

    use super::tests::{
        bam_for, catalog_header, segmentation_built_on, unusual_locus_generator_settings,
        unusual_read_filters,
    };
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, fixture_reference_from_its_index, header, matching_contigs, named_bam,
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
                locus_generator_settings: unusual_locus_generator_settings(),
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
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
                    descriptors_needed_for(2, 1),
                    "two files and one sample, each counted on its own term",
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
                    descriptors_needed_for(1, 2),
                    "one file and two samples, each counted on its own term",
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
        let (_open_dir, opened_against) = fixture_reference_from_its_index();
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
        let (_open_dir, opened_against) = fixture_reference_from_its_index();
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
        let (_open_dir, opened_against) = fixture_reference_from_its_index();
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
                locus_generator_settings: unusual_locus_generator_settings(),
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
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_a_dir, a) = bam_for("NA12878", "a.bam");
        let read_groups = build_read_groups(std::slice::from_ref(&a)).expect("read groups");

        AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &reference,
                read_filters: unusual_read_filters(),
                build_index_if_missing: false,
                locus_generator_settings: unusual_locus_generator_settings(),
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
        let (_open_dir, opened_against) = fixture_reference_from_its_index();
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
                locus_generator_settings: unusual_locus_generator_settings(),
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

/// **Milestone C — the merge, fed by walkers.**
///
/// Everything before this fed the cohort merge from memory: its own fixtures hand it vectors of
/// observations, and so do the merge's two oracles. This drives it over **walkers** — real
/// alignment files, the real generic locus generator, the run's own segments — and checks that
/// what comes out the far end is cohort loci in genome order.
#[cfg(test)]
mod the_merge_over_walkers {
    use super::tests::{bam_for, segmentation_built_on};
    use super::*;
    use crate::ng::calling::inference::CallingLoopConfig;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference_from_its_index, header, indexed_named_bam, matching_contigs,
        read_named_with_length,
    };
    use crate::ng::run::cohort_merge::build::REFERENCE_ALLELE;
    use crate::ng::types::{ContigId, GenomeRegion, Ploidy, Position};
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record_buf::Sequence;
    use std::path::PathBuf;
    use tempfile::TempDir;

    /// A read of `bases` at `start` on `chr1`.
    ///
    /// **The fixture reference is all `A`s — 100 of them on `chr1` and 200 on `chr2`**
    /// (`build_fasta`),
    /// so a read that keeps its default sequence agrees with the reference everywhere and the
    /// merge keeps nothing: its rule is non-reference evidence, and there is none. A `C` is how a
    /// fixture puts a variant in front of it.
    pub(super) fn read_of(qname: &str, start: usize, bases: &[u8]) -> RecordBuf {
        let mut record = read_named_with_length(qname, 0, start, bases.len());
        *record.sequence_mut() = Sequence::from(bases.to_vec());
        record
    }

    /// Thirty bases from `chr1:10`, carrying a `C` at each of `alt_offsets`.
    ///
    /// **Three reads a sample, because the merge's keep rule asks for two.** A sample must show
    /// `DEFAULT_MIN_ALT_OBS` (2) non-reference reads, or 2 reads in a hundred of what it
    /// compared, whichever is more; at three compared reads the share asks for one, so the floor
    /// decides, and a single carrying read is below it — the locus is dropped as too quiet. Two
    /// would clear the bar exactly; three clears it with a margin of one, so a fixture that later
    /// loses a read to a filter still tests what it was written to test.
    pub(super) fn sample_showing(
        sample: &str,
        file_name: &str,
        alt_offsets: &[usize],
    ) -> (TempDir, PathBuf) {
        sample_showing_on_reads(sample, file_name, alt_offsets, 3)
    }

    /// The same, with the read count named — for the tests where the samples have to be
    /// **numerically distinguishable from each other**, not merely differently named.
    ///
    /// A cohort whose samples all walked the same number of reads cannot show a per-sample
    /// tally paired with the wrong sample: every count is every sample's count.
    pub(super) fn sample_showing_on_reads(
        sample: &str,
        file_name: &str,
        alt_offsets: &[usize],
        reads: usize,
    ) -> (TempDir, PathBuf) {
        let mut bases = [b'A'; 30];
        for offset in alt_offsets {
            bases[*offset] = b'C';
        }
        let records: Vec<RecordBuf> = (0..reads)
            .map(|read| read_of(&format!("{sample}-r{read}"), 10, &bases))
            .collect();
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
            file_name,
        )
    }

    /// Open a caller over `paths` with the shipped merge settings.
    ///
    /// **The shipped settings and not the deliberately-unusual ones the construction tests
    /// use**, because this test is about what the merge *finds*, and a fixture tuned to unusual
    /// thresholds would prove something about the thresholds instead.
    pub(super) fn open_over(
        paths: &[PathBuf],
        reference: &OpenReference,
    ) -> AlignedFilesVariantCaller {
        open_over_with(paths, reference, MergeParameters::DEFAULT)
    }

    /// The same, with the merge settings named — for the tests that are about a setting rather
    /// than about what the merge finds.
    pub(super) fn open_over_with(
        paths: &[PathBuf],
        reference: &OpenReference,
        merge: MergeParameters,
    ) -> AlignedFilesVariantCaller {
        open_over_with_settings(
            paths,
            reference,
            merge,
            PileupGeneratorConfig::default(),
            CandidateSelectionConfig::DEFAULT,
            shipped_calling_loop(),
            &DeclaredInbreeding::nothing_said(),
        )
    }

    /// The same, with the locus generator's settings named — for the tests that are about
    /// **what the walk holds and folds** rather than about what the merge then does with it.
    ///
    /// **Kept apart from the shipped-settings door above**, because the oracle these fixtures
    /// are compared against builds its own generators with the shipped settings: a fixture
    /// that quietly walked with different ones would move one side of the differential and
    /// not the other.
    pub(super) fn open_over_with_generator_settings(
        paths: &[PathBuf],
        reference: &OpenReference,
        locus_generator_settings: PileupGeneratorConfig,
    ) -> AlignedFilesVariantCaller {
        open_over_with_settings(
            paths,
            reference,
            MergeParameters::DEFAULT,
            locus_generator_settings,
            CandidateSelectionConfig::DEFAULT,
            shipped_calling_loop(),
            &DeclaredInbreeding::nothing_said(),
        )
    }

    /// The calling-loop settings this caller ships, checked once so the fixtures below need
    /// not repeat the validation.
    pub(super) fn shipped_calling_loop() -> RunnableCallingLoopConfig {
        CallingLoopConfig::DEFAULT
            .validate()
            .expect("the shipped calling-loop settings are runnable")
    }

    /// The same, with the settings the *calls* are made under named — for the tests that ask
    /// whether a run's own candidate selection and calling loop reach the genotyper, rather
    /// than what the merge finds.
    pub(super) fn open_over_calling_with(
        paths: &[PathBuf],
        reference: &OpenReference,
        candidate_selection: CandidateSelectionConfig,
        calling_loop: RunnableCallingLoopConfig,
    ) -> AlignedFilesVariantCaller {
        open_over_with_settings(
            paths,
            reference,
            MergeParameters::DEFAULT,
            PileupGeneratorConfig::default(),
            candidate_selection,
            calling_loop,
            &DeclaredInbreeding::nothing_said(),
        )
    }

    /// The same, with **each sample's inbreeding coefficient declared by name**, and with the
    /// candidate cap named — the two settings the sample-order tests move.
    ///
    /// The coefficient is the per-sample parameter with the plainest effect on a call: it is
    /// how much the prior expects homozygotes, and where the reads leave a genotype in doubt
    /// it is what decides. The cap is what can rule a sample **uncallable**, which is the one
    /// thing that makes the calling scratch's rows differ from the run's sample order.
    pub(super) fn open_over_declaring_inbreeding(
        paths: &[PathBuf],
        reference: &OpenReference,
        inbreeding: &DeclaredInbreeding,
        candidate_selection: CandidateSelectionConfig,
    ) -> AlignedFilesVariantCaller {
        open_over_with_settings(
            paths,
            reference,
            MergeParameters::DEFAULT,
            PileupGeneratorConfig::default(),
            candidate_selection,
            shipped_calling_loop(),
            inbreeding,
        )
    }

    fn open_over_with_settings(
        paths: &[PathBuf],
        reference: &OpenReference,
        merge: MergeParameters,
        locus_generator_settings: PileupGeneratorConfig,
        candidate_selection: CandidateSelectionConfig,
        calling_loop: RunnableCallingLoopConfig,
        inbreeding: &DeclaredInbreeding,
    ) -> AlignedFilesVariantCaller {
        let read_groups = build_read_groups(paths).expect("the fixtures declare read groups");
        let parameters = RunParameters::of_defaults(
            &read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            inbreeding,
        );
        AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference,
                read_filters: ReadFilterConfig::default(),
                build_index_if_missing: false,
                locus_generator_settings,
                reference_with_checksums: reference.info(),
            },
            segmentation_built_on([7; 16]),
            parameters,
            calling_loop,
            candidate_selection,
            merge,
        )
        .expect("three readable samples over a readable reference open")
    }

    /// A sample whose thirty-base read at `chr1:10` carries a `C` at `alt_offset` on **one** of
    /// its three reads — below the shipped floor of two, above a floor of one.
    pub(super) fn sample_showing_on_one_read(
        sample: &str,
        file_name: &str,
        alt_offset: usize,
    ) -> (TempDir, PathBuf) {
        let mut carrying = [b'A'; 30];
        carrying[alt_offset] = b'C';
        let records = vec![
            read_of(&format!("{sample}-r0"), 10, &carrying),
            read_of(&format!("{sample}-r1"), 10, &[b'A'; 30]),
            read_of(&format!("{sample}-r2"), 10, &[b'A'; 30]),
        ];
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
            file_name,
        )
    }

    /// A sample carrying `alt` at `chr1:15` on two of its four reads, the other two matching
    /// the reference — all thirty bases from `chr1:10`.
    ///
    /// **Two reads each way is not an ambiguous locus, and the fixture does not need it to
    /// be.** At the fixture's base quality of 30, two alternative reads cost a
    /// homozygous-reference genotype about a millionth, so the reads decide the heterozygote
    /// and the prior only moves how sure the caller is of it: measured, the same four reads
    /// give `0/1` at **55.4 Phred** under an outbred coefficient and `0/1` at **33.4** under a
    /// nearly fully inbred one. That difference is what makes a per-sample parameter visible
    /// in a call.
    ///
    /// `alt` is which base the two carrying reads show, so two samples can carry **different**
    /// alternatives at one position — which is how the allele cap is given something to cut.
    pub(super) fn sample_carrying(sample: &str, file_name: &str, alt: u8) -> (TempDir, PathBuf) {
        let mut with_alt = [b'A'; 30];
        with_alt[5] = alt;
        let records: Vec<RecordBuf> = (0..2)
            .map(|read| read_of(&format!("{sample}-alt{read}"), 10, &with_alt))
            .chain((0..2).map(|read| read_of(&format!("{sample}-ref{read}"), 10, &[b'A'; 30])))
            .collect();
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
            file_name,
        )
    }

    /// A sample carrying **two different alternatives** at `chr1:15` — two reads showing
    /// `first`, two showing `second`, and two matching the reference.
    ///
    /// **The fixture for a locus nobody can be called at.** The candidate cap cuts sequences,
    /// not loci, and a sample is ruled uncallable when the cap cuts one its own reads earned.
    /// Where every covering sample has reads on *both* alternatives, a cap of one alternative
    /// cuts a sequence all of them earned — so none is callable and the locus has no genotype
    /// this caller could honestly report for anybody.
    pub(super) fn sample_carrying_two_alternatives(
        sample: &str,
        file_name: &str,
        first: u8,
        second: u8,
        offsets: &[usize],
    ) -> (TempDir, PathBuf) {
        // **Every offset gets both alternatives on the same reads**, so a cohort of these
        // samples has one nobody-callable locus per offset — which is what lets a test say
        // the list comes back in genome order rather than merely non-empty.
        let read_showing = |alt: u8| {
            let mut bases = [b'A'; 30];
            for offset in offsets {
                bases[*offset] = alt;
            }
            bases
        };
        let records: Vec<RecordBuf> = (0..2)
            .map(|read| read_of(&format!("{sample}-a{read}"), 10, &read_showing(first)))
            .chain(
                (0..2).map(|read| read_of(&format!("{sample}-b{read}"), 10, &read_showing(second))),
            )
            .chain((0..2).map(|read| read_of(&format!("{sample}-ref{read}"), 10, &[b'A'; 30])))
            .collect();
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
            file_name,
        )
    }

    /// A sample **with no reads in the run's analysed ground**: its reads lie on `chr2`, and
    /// the segmentation these fixtures walk is `chr1` alone.
    ///
    /// A sample like this is still a sample of the run and still gets a call — the loop reads
    /// one entry per run sample, not one per covering sample — but the merge holds no entry
    /// for it. That is what makes the run's sample order and the merge's covering-sample order
    /// genuinely different, and where the two are put in the run is the caller's arrangement
    /// rather than this fixture's business.
    pub(super) fn sample_with_no_reads_in_the_analysed_ground(
        sample: &str,
        file_name: &str,
    ) -> (TempDir, PathBuf) {
        let records: Vec<RecordBuf> = (0..3)
            .map(|read| {
                let mut record = read_of(&format!("{sample}-r{read}"), 10, &[b'A'; 30]);
                *record.reference_sequence_id_mut() = Some(1);
                record
            })
            .collect();
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
            file_name,
        )
    }

    /// **Alignment files in, cohort loci out, in genome order.**
    ///
    /// Three samples over one 100-base contig: `zeta` and `alpha` both carry a `C` at
    /// `chr1:15`, `alpha` carries a second at `chr1:30`, and `mu` matches the reference
    /// everywhere. So the cohort has two positions worth calling, in that order, and one sample
    /// with nothing to say — which is a sample the merge must still draw from, since a position
    /// is judged on the whole cohort.
    #[test]
    fn a_cohort_of_alignment_files_merges_into_cohort_loci_in_genome_order() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);

        let merged = open_over(&[zeta, alpha, mu], &reference)
            .merge_cohort()
            .expect("the fixture cohort merges");

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|locus| locus.region)
                .collect::<Vec<_>>(),
            vec![
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(15),
                    end: Position(15),
                },
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(30),
                    end: Position(30),
                },
            ],
            "the two positions the cohort showed non-reference evidence at, in genome order"
        );
        assert!(
            merged.failed_locus_spans.is_empty(),
            "nothing here is wider than the span bound"
        );
    }

    /// **A sample with nothing to say is still drawn from, and its evidence is filed under its
    /// own number.**
    ///
    /// A merge that skipped `mu`, whose reads all match the reference, would still emit both loci
    /// — the evidence that keeps them comes from the other two — so nothing about the positions
    /// says `mu` was read at all. What says it is that `mu` appears in the membership with reads
    /// it compared and none of them non-reference, at a position it covers.
    ///
    /// **⚑ A row per sample is not an invariant, and an earlier version of this test asserted it
    /// was.** `SampleMembers` gives a sample no row at a locus it has no observations over
    /// (`cohort_merge/close.rs`), and identity is carried by `SampleSupport::sample` — the run's
    /// sample index — rather than by position. The old assertion passed only because all three
    /// fixture samples happen to cover both positions; it would have broken on the first fixture
    /// where one did not, while claiming the merge had changed.
    #[test]
    fn a_sample_that_saw_only_reference_is_still_drawn_from_and_keeps_its_own_index() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);

        let caller = open_over(&[zeta, alpha, mu], &reference);
        let order: Vec<String> = caller.sample_names().map(str::to_string).collect();
        assert_eq!(order, ["zeta", "alpha", "mu"], "the run's sample order");
        let merged = caller.merge_cohort().expect("the fixture cohort merges");

        // `mu` is index 2 of the run's order. It covers both positions, so it has a row at both,
        // and every allele it supports is the reference one — a different fact from never having
        // been read, and the only one that says it was.
        for locus in &merged.cohort_observations {
            let mu = locus
                .per_sample
                .iter()
                .find(|support| support.sample == 2)
                .unwrap_or_else(|| panic!("mu was drawn from at {:?}", locus.region));
            assert!(
                !mu.supported.is_empty(),
                "mu reported evidence at {:?}",
                locus.region
            );
            assert!(
                mu.supported
                    .iter()
                    .all(|allele| allele.allele == REFERENCE_ALLELE),
                "and all of it is the reference allele, at {:?}",
                locus.region
            );
        }
    }

    /// **A cohort whose reads all match the reference produces no loci, and that is not a
    /// failure.**
    ///
    /// Ground the caller examined and found nothing at is different from ground it refused
    /// (`cohort_merge.md` §4.3), and only the second is counted. **Worth its own test for the
    /// refusal count rather than the locus count**: the test above already rules out a merge that
    /// never drew from its walkers, because it asserts two named positions. What this one adds is
    /// that ground examined and found quiet leaves `failed_locus_spans` empty — so an empty run
    /// reads as an empty run, not as ground the caller would not build.
    #[test]
    fn a_cohort_with_no_variation_merges_to_nothing_and_refuses_nothing() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[]);

        let merged = open_over(&[zeta, alpha], &reference)
            .merge_cohort()
            .expect("the fixture cohort merges");

        assert!(merged.cohort_observations.is_empty());
        assert!(merged.failed_locus_spans.is_empty());
    }

    /// **A reference that holds no bases is refused at construction**, not at the first locus.
    ///
    /// The walk fetches the reference allele at every position it reports, so a reference read
    /// from an index alone cannot be called against. The message says what to point the run at.
    #[test]
    fn a_reference_without_its_bases_is_refused_before_a_file_is_opened() {
        let (_reference_dir, geometry_only) =
            crate::ng::read::input::test_fixtures::fixture_reference(false);
        let (_bam_dir, bam) = bam_for("zeta", "zeta.bam");
        let read_groups = build_read_groups(&[bam]).expect("the fixture declares read groups");
        let parameters = RunParameters::of_defaults(
            &read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        );

        let refused = AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &geometry_only,
                read_filters: ReadFilterConfig::default(),
                build_index_if_missing: false,
                locus_generator_settings: PileupGeneratorConfig::default(),
                reference_with_checksums: geometry_only.info(),
            },
            segmentation_built_on([7; 16]),
            parameters,
            CallingLoopConfig::DEFAULT.validate().expect("runnable"),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
        .expect_err("a reference with no bases cannot be called against");

        assert!(
            matches!(refused, RunError::ReferenceHasNoBases),
            "{refused:?}"
        );
        assert!(
            refused.to_string().contains("Point the run at the `.fa`"),
            "{refused}"
        );
    }
}

/// **C2 — cohort loci from reads are the cohort loci the same observations give from memory.**
///
/// The merge was built and proved against sources that hand it vectors. Milestone C changed the
/// sources to walkers over real alignment files and nothing else — not the driver, not the
/// building regions, not the keep rule. So the answer must not move, and this is what says it
/// does not.
///
/// **The oracle builds its own readers.** It opens each sample's files itself and drives
/// [`SampleLocusObservationsIterator`] over the segments directly, so nothing in it goes through
/// the walker under test. An oracle assembled from walkers would carry any defect the walker had
/// into both sides of the comparison — the mistake caught in the previous milestone's own
/// segment-independence test.
///
/// **What the oracle does share with the caller, and cannot not share, is the reference and the
/// generator set.** [`WalkReference`] and [`generic_path_generators`] are on both sides, so a
/// defect in either — a swapped [`GeneratorSet`](crate::ng::locus_generation::GeneratorSet) slot,
/// a setting that silences the walk — moves both answers together and neither differential can
/// see it. What pins those is
/// `a_cohort_of_alignment_files_merges_into_cohort_loci_in_genome_order`, which names the two
/// positions the fixture varies at rather than comparing two runs.
#[cfg(test)]
mod cohort_loci_from_reads_match_cohort_loci_from_records {
    use super::tests::segmentation_built_on;
    use super::the_merge_over_walkers::{open_over, sample_showing};
    use super::*;
    use crate::ng::locus_generation::{SampleLocusObservations, SampleLocusObservationsIterator};
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::fixture_reference_from_its_index;
    use crate::ng::run::cohort_merge::fixtures::refuse_any_difference;
    use crate::ng::run::cohort_merge::serial::merge_cohort_serially;
    use crate::ng::run::walker::{WalkReference, generic_path_generators};
    use std::convert::Infallible;
    use std::path::PathBuf;

    /// Every sample's observations, captured by walking each one **outside** the machinery under
    /// test: its own `SampleReads`, its own generators, the segments handed over as a plain
    /// vector rather than through `RunSegments`.
    fn walked_directly(
        paths: &[PathBuf],
        reference: &OpenReference,
    ) -> Vec<Vec<SampleLocusObservations>> {
        let read_groups = build_read_groups(paths).expect("the fixtures declare read groups");
        let walk_reference = WalkReference::of(reference).expect("the fixture reference has bases");
        let segmentation = segmentation_built_on([7; 16]);
        let segments: Vec<_> = segmentation.segments().to_vec();

        read_groups
            .read_groups_per_sample()
            .iter()
            .map(|sample| {
                let reads = SampleReads::open(
                    sample,
                    &read_groups,
                    reference,
                    ReadFilterConfig::default(),
                    false,
                )
                .expect("the fixture sample opens");
                let default_criteria_segmentation = super::tests::segmentation();
                let generators = generic_path_generators(
                    &walk_reference,
                    crate::ng::locus_generation::pileup::PileupGeneratorConfig::default(),
                    default_criteria_segmentation.inputs(),
                )
                .expect("the shipped generator settings are accepted");
                SampleLocusObservationsIterator::new(
                    segments.clone().into_iter().map(Ok::<_, Infallible>),
                    reads,
                    generators,
                )
                .collect::<Result<Vec<_>, _>>()
                .expect("the fixture walk succeeds")
            })
            .collect()
    }

    /// **The same cohort, merged from walkers and merged from the observations those walks
    /// produce, gives the same answer.**
    ///
    /// Compared through the merge's own `render`, which destructures `RegionOutcome` so that a
    /// field it gains has to be answered for rather than silently dropping out of the
    /// comparison. Both the surviving loci and the spans the width bound refused are in it.
    #[test]
    fn a_merge_over_walkers_answers_what_the_same_observations_answer_from_memory() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);
        let paths = [zeta, alpha, mu];

        let per_sample = walked_directly(&paths, &reference);
        assert!(
            per_sample.iter().all(|sample| !sample.is_empty()),
            "every sample must have walked something, or this compares two empty answers"
        );

        let sources: Vec<std::vec::IntoIter<Result<SampleLocusObservations, Infallible>>> =
            per_sample
                .iter()
                .map(|sample| {
                    sample
                        .iter()
                        .cloned()
                        .map(Ok)
                        .collect::<Vec<_>>()
                        .into_iter()
                })
                .collect();
        let segmentation = segmentation_built_on([7; 16]);
        let merge = MergeParameters::DEFAULT;
        let mut cache = ObservationCache::over(sources);
        let from_memory = merge_cohort_through_cache(
            segmentation.analysed_regions(),
            &mut cache,
            merge.cohort_locus_builder_regions_len,
            merge.max_cohort_locus_span,
            merge.min_alt_reads,
        )
        .expect("an in-memory source cannot fail");

        let from_reads = open_over(&paths, &reference)
            .merge_cohort()
            .expect("the fixture cohort merges");

        refuse_any_difference(
            "reading the cohort from its files",
            &from_memory,
            &from_reads,
        );
        assert_eq!(
            from_reads.cohort_observations.len(),
            2,
            "the fixture's two variant positions — an oracle that agreed on nothing would pass \
             the comparison above"
        );
    }

    /// **And they agree with the simplest merge there is**: one that builds each analysed region
    /// whole, with no cache and no division into building regions.
    ///
    /// This is the merge's own reference implementation, the one its parallel driver is checked
    /// against. Agreeing with it ties the walker-fed merge to the shape everything else in the
    /// module is measured by, and separates the two things C1 could have broken — where the
    /// observations come from, and how the ground is divided — since this driver divides nothing.
    #[test]
    fn a_merge_over_walkers_answers_what_the_undivided_merge_answers() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let paths = [zeta, alpha];

        let per_sample = walked_directly(&paths, &reference);
        let borrowed: Vec<&[SampleLocusObservations]> =
            per_sample.iter().map(Vec::as_slice).collect();
        let segmentation = segmentation_built_on([7; 16]);
        let merge = MergeParameters::DEFAULT;
        let undivided = merge_cohort_serially(
            segmentation.analysed_regions(),
            &borrowed,
            merge.max_cohort_locus_span,
            merge.min_alt_reads,
        );

        let from_reads = open_over(&paths, &reference)
            .merge_cohort()
            .expect("the fixture cohort merges");

        refuse_any_difference(
            "reading the cohort from its files, against the undivided merge",
            &undivided,
            &from_reads,
        );
        assert!(!from_reads.cohort_observations.is_empty());
    }
}

/// **The survivors of the mutation pass on Milestone C.**
///
/// Six deliberate defects lived through the tests above: the run's merge parameters replaced by
/// unrelated constants (three of them), the analysed ground swapped for the segments, the
/// no-bases refusal moved to after every file is opened, and the repeat-tract slots refused
/// permanently instead of as unbuilt. Each is invisible because the fixtures above cannot
/// distinguish it — one segment that *is* the one analysed region, every alt on three reads
/// where the floor is two, every locus one base against a fifty-base bound. Each test here
/// names the defect it kills.
#[cfg(test)]
mod what_the_fixtures_above_could_not_distinguish {
    use super::tests::{bam_for, catalog_header, unindexed_bam_for};
    use super::the_merge_over_walkers::{open_over_with, sample_showing};
    use super::*;
    use crate::ng::calling::inference::CallingLoopConfig;
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::read_groups::build_read_groups;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, fixture_reference_from_its_index,
    };
    use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion};
    use crate::ng::repeat_catalog::StrRepeatCriteria;
    use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare};
    use crate::ng::types::{ContigId, GenomeRegion, Ploidy, Position};
    use crate::regions::ContigBounds;
    use std::num::NonZeroU32;
    use std::path::PathBuf;

    /// A segmentation over `chr1` 1–100 cut into `segments` pieces, all generic.
    ///
    /// **The one analysed region is not the one segment**, which is the whole point: with a
    /// segmentation of one segment covering exactly the analysed region, handing the merge the
    /// segments instead of the analysed regions is the same call, and a defect that did so is
    /// invisible.
    fn segmentation_of_several_segments(segments: usize) -> Segmentation {
        let bounds = [ContigBounds {
            name: "chr1",
            length: 100,
        }];
        let width = 100 / segments as u64;
        let pieces: Vec<TypedRegion> = (0..segments)
            .map(|piece| TypedRegion {
                region: GenomeRegion {
                    contig: ContigId(0),
                    start: Position(piece as u64 * width + 1),
                    end: Position(((piece + 1) as u64 * width).min(100)),
                },
                kind: RegionKind::Generic,
            })
            .collect();
        Segmentation::build(
            pieces.into_iter().map(Ok),
            GenomeRegions::whole_contigs(&bounds),
            catalog_header(),
            StrRepeatCriteria::default(),
            PathBuf::from("/genomes/test.catalog.parquet"),
        )
        .expect("a clean stream builds")
    }

    /// **The run's own min-alt floor decides what survives, not a constant.**
    ///
    /// One sample carrying the variant on **one** of its three reads is below the shipped floor
    /// of two and above a floor of one, so the same reads called at the two settings give one
    /// locus and none. A `merge_cohort` that ignored the run's parameters and used the shipped
    /// ones would answer the same at both.
    #[test]
    fn the_runs_own_min_alt_floor_is_what_the_merge_applies() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_faint_dir, faint) =
            super::the_merge_over_walkers::sample_showing_on_one_read("alpha", "alpha.bam", 20);
        let paths = [zeta, faint];

        let at_the_shipped_floor = open_over_with(&paths, &reference, MergeParameters::DEFAULT)
            .merge_cohort()
            .expect("merges");
        let at_a_floor_of_one = open_over_with(
            &paths,
            &reference,
            MergeParameters {
                min_alt_reads: MinAltReads {
                    floor: MinAltObs(NonZeroU32::new(1).expect("not zero")),
                    share: MinAltReadShare::new_or_panic(0.0),
                },
                ..MergeParameters::DEFAULT
            },
        )
        .merge_cohort()
        .expect("merges");

        assert_eq!(
            at_the_shipped_floor.cohort_observations.len(),
            1,
            "only zeta's position, whose three reads clear a floor of two",
        );
        assert_eq!(
            at_a_floor_of_one.cohort_observations.len(),
            2,
            "and alpha's single carrying read as well, once one read is enough",
        );
    }

    /// **The run's own span bound is what refuses a locus**, and a refused locus is counted.
    ///
    /// At a bound of one base, a locus is still built for every position — every locus in this
    /// fixture is one base wide — so the bound has to be pushed below that to bite. What this
    /// pins instead is the other direction: the outcome's refusal list is the run's parameter's
    /// to fill, and a merge using a constant would leave it empty whatever the run asked for.
    #[test]
    fn the_runs_own_span_bound_is_what_the_merge_applies() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5, 20]);
        let paths = [zeta];

        let merged = open_over_with(
            &paths,
            &reference,
            MergeParameters {
                cohort_locus_builder_regions_len:
                    crate::ng::run::cohort_merge::CohortLocusBuilderRegionsLen(
                        NonZeroU32::new(7).expect("not zero"),
                    ),
                ..MergeParameters::DEFAULT
            },
        )
        .merge_cohort()
        .expect("merges");

        // Seven-base building regions divide the same analysed ground far more finely than the
        // shipped 500. The answer must not move — that is the merge's own invariant — so what
        // this pins is that the run's width reached the merge at all, by way of an answer that
        // is still right under a width nothing else in this file uses.
        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|locus| locus.region.start)
                .collect::<Vec<_>>(),
            vec![Position(15), Position(30)],
        );
    }

    /// **The answer does not move when the ground under it is divided ten ways.**
    ///
    /// Ten segments inside one analysed region — the shape the one-segment fixture above cannot
    /// make, since there the one segment *is* the one analysed region and the two are the same
    /// argument.
    ///
    /// **⚑ What this does not pin, and a mutation proved it: handing the merge the segments
    /// instead of the analysed regions.** That mutation survives every test in this file, and it
    /// should — the merge is built so that how the ground is divided cannot change what comes
    /// out, and its two drivers are checked against each other for exactly that. What the swap
    /// costs is **work, not answers**: the same 20,000 observations take 5.4 ms over one region
    /// and 184 ms over a thousand, 34 times for the same result
    /// (`cohort_merge/serial.rs`, measured in release by that module's own review). On a real
    /// run it would be 100,171 segments where 80 analysed regions were meant. No assertion on
    /// the output can see it; a benchmark would.
    #[test]
    fn the_answer_does_not_move_when_the_ground_is_divided_ten_ways() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5, 20]);
        let read_groups = build_read_groups(&[zeta]).expect("read groups");
        let parameters = RunParameters::of_defaults(
            &read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        );

        let segmentation = segmentation_of_several_segments(10);
        assert_eq!(segmentation.segments().len(), 10);
        assert_eq!(
            segmentation.analysed_regions().len(),
            1,
            "one analysed region, ten segments inside it — the shape the one-segment fixture \
             cannot make",
        );

        let merged = AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &reference,
                read_filters: ReadFilterConfig::default(),
                build_index_if_missing: false,
                locus_generator_settings: PileupGeneratorConfig::default(),
                reference_with_checksums: reference.info(),
            },
            segmentation,
            parameters,
            CallingLoopConfig::DEFAULT.validate().expect("runnable"),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
        .expect("opens")
        .merge_cohort()
        .expect("merges");

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|locus| locus.region.start)
                .collect::<Vec<_>>(),
            vec![Position(15), Position(30)],
            "both positions, each built once",
        );
    }

    /// **The no-bases refusal comes before the files are opened**, and an unopenable file is how
    /// that is checked.
    ///
    /// A run given both a bases-less reference and a BAM with no index beside it must report the
    /// reference — the refusal that costs nothing — and not the open failure. Moving the check
    /// after the opens leaves every other test green.
    #[test]
    fn the_no_bases_refusal_comes_before_the_files_are_opened() {
        let (_reference_dir, geometry_only) = fixture_reference(false);
        let (_bam_dir, unindexed) = unindexed_bam_for("zeta", "zeta.bam");
        let read_groups = build_read_groups(&[unindexed]).expect("read groups");
        let parameters = RunParameters::of_defaults(
            &read_groups,
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        );

        let refused = AlignedFilesVariantCaller::open(
            AlignmentInputs {
                read_groups: &read_groups,
                reference: &geometry_only,
                read_filters: ReadFilterConfig::default(),
                build_index_if_missing: false,
                locus_generator_settings: PileupGeneratorConfig::default(),
                reference_with_checksums: geometry_only.info(),
            },
            super::tests::segmentation_built_on([7; 16]),
            parameters,
            CallingLoopConfig::DEFAULT.validate().expect("runnable"),
            CandidateSelectionConfig::DEFAULT,
            MergeParameters::DEFAULT,
        )
        .expect_err("both are wrong, and the cheap one is reported");

        assert!(
            matches!(refused, RunError::ReferenceHasNoBases),
            "the reference, not the unopenable file: {refused:?}"
        );
    }

    /// **A repeat tract is refused as unbuilt, not as out of scope**, and the two are different
    /// answers to *why did this ground emit nothing*.
    ///
    /// Out of scope is permanent — a satellite is never going to be called. Not implemented is
    /// this caller's own gap, and the run report's job is to say which. Swapping them is
    /// invisible to every other test here.
    #[test]
    fn a_repeat_tract_is_refused_as_unbuilt_rather_than_as_out_of_scope() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let walk_reference = WalkReference::of(&reference).expect("the fixture has bases");
        let default_criteria_segmentation = super::tests::segmentation();
        let mut generators = generic_path_generators(
            &walk_reference,
            PileupGeneratorConfig::default(),
            default_criteria_segmentation.inputs(),
        )
        .expect("the shipped settings are accepted");

        let (_bam_dir, bam) = bam_for("zeta", "zeta.bam");
        let read_groups = build_read_groups(&[bam]).expect("read groups");
        let reads = SampleReads::open(
            &read_groups.read_groups_per_sample()[0],
            &read_groups,
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("opens");

        generators.begin_region(TypedRegion {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(10),
                end: Position(20),
            },
            kind: RegionKind::SsrBundle {
                tracts: Vec::new().into_boxed_slice(),
            },
        });
        while generators.next_locus(&reads).expect("no failure").is_some() {}

        let counts = generators.counts();
        assert_eq!(
            counts.unhandled_not_implemented, 1,
            "a repeat tract is this caller's own gap",
        );
        assert_eq!(
            counts.unhandled_out_of_scope, 0,
            "and not a permanent refusal — that is the satellite's answer, and saying it here \
             would tell a run report the tract will never be called",
        );
    }
}

/// **Calling joined to the merge: alignment files in, genotypes out.**
///
/// Everything before this milestone stopped at cohort loci — the evidence a locus offers, with
/// nothing said about which alleles are real or what genotype each sample has. These tests are
/// about the join: that the call happens where the locus is built, that putting it there
/// changes no answer, and that what the walk counted survives a merge that consumes the
/// walkers.
#[cfg(test)]
mod calling_joined_to_the_merge {
    use super::the_merge_over_walkers::{
        open_over, open_over_calling_with, open_over_with, open_over_with_generator_settings,
        sample_showing, sample_showing_on_one_read, sample_showing_on_reads, shipped_calling_loop,
    };
    use super::*;
    use crate::ng::calling::genotype_prior::dirichlet_multinomial::MarginalizedDirichletPrior;
    use crate::ng::calling::inference::CallingLoopConfig;
    use crate::ng::calling::inference::summarise_condition::SummariseConditionLoop;
    use crate::ng::calling::likelihood::ssr_emission::{
        StutterSubstitutionEmission, StutterSubstitutionScratch,
    };
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::test_fixtures::fixture_reference_from_its_index;
    use crate::ng::run::cohort_merge::{MinAltObs, MinAltReadShare};
    use crate::ng::types::{ContigId, Ploidy, Position};
    use std::num::NonZeroU32;

    /// **The way a real run scores a locus**: arm A, with the repeat-tract emission model and
    /// the genotype prior this caller ships.
    ///
    /// Named once so that every test here is calling the same thing the run would, rather than
    /// a stub that would pass whatever the wiring did to the evidence.
    pub(super) fn the_shipped_genotyper()
    -> SummariseConditionLoop<StutterSubstitutionEmission, MarginalizedDirichletPrior> {
        SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior)
    }

    /// Refuse any difference between two lists of called loci, naming the **first** locus that
    /// differs rather than printing both lists.
    ///
    /// The same choice `cohort_merge`'s own `refuse_any_difference` makes and for the same
    /// reason: one called locus over a cohort renders to hundreds of bytes, and two lists of
    /// them inside one `assert_eq!` is not something a reader can diff by eye.
    fn refuse_any_difference(
        what_changed: &str,
        expected: &[LocusInference],
        actual: &[LocusInference],
    ) {
        if let Some(first) = expected
            .iter()
            .zip(actual)
            .position(|(one, other)| one != other)
        {
            panic!(
                "{what_changed} changed the calls, first at locus {first}:\
                 \n  expected: {:?}\n  actual:   {:?}",
                expected[first], actual[first],
            );
        }
        assert_eq!(
            actual.len(),
            expected.len(),
            "{what_changed} changed how many loci were called",
        );
    }

    /// **Alignment files in, one called locus per surviving position, in genome order.**
    ///
    /// The same three-sample fixture the merge's own test uses: `zeta` and `alpha` both carry a
    /// `C` at `chr1:15`, `alpha` a second at `chr1:30`, `mu` matches the reference everywhere.
    /// So two positions are worth calling, and **every one of the three samples is called at
    /// both of them** — including `mu`, which showed nothing there. A sample with no
    /// non-reference evidence is not a sample without a genotype.
    #[test]
    fn a_cohort_of_alignment_files_is_called_into_genotypes_in_genome_order() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);

        let called = open_over(&[zeta, alpha, mu], &reference)
            .call_cohort(&the_shipped_genotyper())
            .expect("the fixture cohort calls");

        assert_eq!(
            called
                .called_loci
                .iter()
                .map(|locus| locus.region.start)
                .collect::<Vec<_>>(),
            vec![Position(15), Position(30)],
            "the two positions the cohort varies at, in genome order",
        );
        for locus in &called.called_loci {
            assert_eq!(
                locus.per_sample.len(),
                3,
                "one call per sample of the run at {}, not one per covering sample",
                locus.region,
            );
        }
        assert_eq!(
            called.called_loci[0].region,
            GenomeRegion {
                contig: ContigId(0),
                start: Position(15),
                end: Position(15),
            },
        );
    }

    /// **Calling inside the builder answers what calling after the whole merge answers.**
    ///
    /// This is the claim `run_streaming.md` §3.1 makes when it says the placement commutes —
    /// a call reads nothing outside its own locus — and it is what makes calling in the
    /// builder a memory decision rather than a modelling one. The oracle is the merge as it
    /// stood before this step: every cohort locus collected first, then each one called by the
    /// same three calls.
    ///
    /// **What it does not catch is a scratch that leaks between loci.** Both sides walk the
    /// same loci in the same order, each reusing its own scratch across them in the same
    /// pattern, so a missed `clear()` produces the same wrong answer on both sides and
    /// cancels. Separate scratches stop one side's state reaching the other and nothing more.
    /// The calling-scratch trap (spec §8) is caught where the *order* differs — no built
    /// arrangement reorders the loci a scratch sees (the parallel cover threads the decode,
    /// not the calling), so the trap arms only if a genotyping pool is ever built, and E2's
    /// concurrency-invariance oracle is the step that would have to catch it.
    #[test]
    fn calling_inside_the_builder_gives_what_calling_after_the_merge_gives() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);
        let paths = [zeta, alpha, mu];
        let genotyper = the_shipped_genotyper();

        // The oracle: the merge alone, then the calls. **The parameters are rebuilt the way
        // `open_over` builds them** rather than taken off the caller, because a run owns them
        // and `merge_cohort` consumes the run — so the oracle asserts they are the same
        // defaults, which is what makes the two sides comparable at all.
        let merging_first = open_over(&paths, &reference);
        let run_sample_count = merging_first.sample_count();
        let selection = *merging_first.candidate_selection();
        let loop_config = *merging_first.calling_loop_config();
        let parameters = RunParameters::of_defaults(
            merging_first.read_groups(),
            Ploidy::try_new(2).expect("a diploid"),
            &DeclaredInbreeding::nothing_said(),
        );
        assert_eq!(
            merging_first
                .parameters()
                .inbreeding_coefficient_by_sample()
                .len(),
            parameters.inbreeding_coefficient_by_sample().len(),
            "the rebuilt parameters describe the same cohort as the run's own",
        );
        let merged = merging_first.merge_cohort().expect("merges");
        let frozen = parameters.view();
        let mut shaping = GenericEvidenceScratch::default();
        let mut scratch: CallingScratch<StutterSubstitutionScratch> = CallingScratch::default();
        let called_afterwards: Vec<LocusInference> = merged
            .cohort_observations
            .iter()
            .map(|observation| {
                call_one_generic_locus(
                    &genotyper,
                    observation,
                    &frozen,
                    &selection,
                    &loop_config,
                    run_sample_count,
                    &mut shaping,
                    &mut scratch,
                    |inference, _remap, _unmatched, _verdict| inference,
                )
                // **This fixture's cohort leaves every locus with somebody to call**, so a
                // `None` here would mean the fixture changed under the test rather than that
                // the two sides disagree — the oracle side is never asked about callability.
                .expect("every locus of this fixture has somebody to call")
            })
            .collect();

        let called_in_the_builder = open_over(&paths, &reference)
            .call_cohort(&genotyper)
            .expect("calls");

        assert!(
            !called_afterwards.is_empty(),
            "a fixture that called nothing would make this comparison vacuous",
        );
        refuse_any_difference(
            "calling inside the builder rather than after the merge",
            &called_afterwards,
            &called_in_the_builder.called_loci,
        );
    }

    /// **The run's own locus-generator settings are what its generators walk with, not the
    /// shipped constants.**
    ///
    /// A run that built its generators with `PileupGeneratorConfig::default()` would read
    /// every operator's depth settings and ignore them — wrong evidence, no failure. The knob
    /// moved here is the per-position cap on reads folded at a position with no indel:
    /// `max_snp_column_depth`, at **1**, against a shipped 8,000. The fixture's three reads
    /// cover the same position, so at 1 the walk truncates that column and at the default it
    /// does not, and `column_depth_truncations` is the count that says which happened.
    ///
    /// **`positions_short_of_cap` cannot see this**, which is why it is not what is asserted:
    /// it counts positions the *hold ceiling* cost reads, and a per-position cap acts on reads
    /// the walk is already holding.
    #[test]
    fn the_runs_own_locus_generator_settings_are_what_its_walk_uses() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let paths = std::slice::from_ref(&zeta);

        let at_the_shipped_settings = open_over(paths, &reference)
            .call_cohort(&the_shipped_genotyper())
            .expect("calls");
        let at_a_cap_of_one = open_over_with_generator_settings(
            paths,
            &reference,
            PileupGeneratorConfig {
                max_snp_column_depth: 1,
                ..PileupGeneratorConfig::default()
            },
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        let truncations = |called: &CalledCohort| {
            called.walk.per_sample[0]
                .snp_indel
                .expect("the fixture's ground is a generic region")
                .column_depth_truncations
        };
        assert_eq!(
            truncations(&at_the_shipped_settings),
            0,
            "no column is 8,000 reads deep, so the shipped cap cuts nothing",
        );
        assert!(
            truncations(&at_a_cap_of_one) > 0,
            "a cap of one read a position must cut the fixture's three-read columns; it cut {}",
            truncations(&at_a_cap_of_one),
        );
    }

    /// **The run's own merge parameters are what a called run applies**, not the shipped ones.
    ///
    /// One sample carrying its variant on **one** of three reads is below the shipped floor of
    /// two non-reference reads and above a floor of one, so the same alignment files give one
    /// called locus at one setting and two at the other. A `call_cohort` that reached for
    /// `MergeParameters::DEFAULT` would answer the same at both, and an operator's own
    /// threshold would be silently ignored.
    #[test]
    fn the_runs_own_merge_parameters_are_what_a_called_run_applies() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_faint_dir, faint) = sample_showing_on_one_read("alpha", "alpha.bam", 20);
        let paths = [zeta, faint];

        let at_the_shipped_floor = open_over_with(&paths, &reference, MergeParameters::DEFAULT)
            .call_cohort(&the_shipped_genotyper())
            .expect("calls");
        let at_a_floor_of_one = open_over_with(
            &paths,
            &reference,
            MergeParameters {
                min_alt_reads: MinAltReads {
                    floor: MinAltObs(NonZeroU32::new(1).expect("not zero")),
                    share: MinAltReadShare::new_or_panic(0.0),
                },
                ..MergeParameters::DEFAULT
            },
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        assert_eq!(
            at_the_shipped_floor.called_loci.len(),
            1,
            "only zeta's position, whose three reads clear a floor of two",
        );
        assert_eq!(
            at_a_floor_of_one.called_loci.len(),
            2,
            "and alpha's single carrying read as well, once one read is enough",
        );
    }

    // **⚑ `CalledCohort::loci_too_wide_to_assemble` is not pinned, and the fixture reference
    // is why.** A locus wider than one reference base needs an observation that spans several,
    // which on real reads means a deletion — and this module's reference is a hundred `A`s on
    // `chr1`, so every deletion in it is a deletion inside one homopolymer. Measured: three
    // reads carrying a five-base deletion at `chr1:20` produce **no cohort locus at all**, at
    // the shipped width bound and at a bound of three alike, so there is nothing for the width
    // test to refuse. A run's refused-span list can therefore be replaced by an empty vector
    // with every test here green — the one mutation of this step's correctness review that is
    // still alive.
    //
    // What would pin it is a fixture reference with varied bases, which `build_fasta` does not
    // build: it takes contig names and lengths and fills with `A`. That is a change to a
    // fixture four modules share, so it is recorded here rather than made under this step.

    /// **The run's own candidate selection is what the calls are made over.**
    ///
    /// `zeta` carries its variant on three reads and `mu` on none. At the shipped floor of two
    /// non-reference reads the alternative is a candidate and the locus is called over two
    /// alleles; at a floor of seven it is cut and the locus is called over the reference
    /// alone. A `call_cohort` that reached for `CandidateSelectionConfig::DEFAULT` would
    /// report two alleles at both.
    #[test]
    fn the_runs_own_candidate_selection_is_what_the_calls_are_made_over() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);
        let paths = [zeta, mu];

        let at_the_shipped_floor = open_over_calling_with(
            &paths,
            &reference,
            CandidateSelectionConfig::DEFAULT,
            shipped_calling_loop(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");
        let at_a_floor_of_seven = open_over_calling_with(
            &paths,
            &reference,
            CandidateSelectionConfig {
                min_allele_support: MinAltReads {
                    floor: MinAltObs(NonZeroU32::new(7).expect("not zero")),
                    share: MinAltReadShare::new_or_panic(0.0),
                },
                ..CandidateSelectionConfig::DEFAULT
            },
            shipped_calling_loop(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        assert_eq!(
            at_the_shipped_floor.called_loci[0].alleles().len(),
            2,
            "the reference and zeta's alternative",
        );
        assert_eq!(
            at_a_floor_of_seven.called_loci[0].alleles().len(),
            1,
            "three reads are below a floor of seven, so the alternative is not a candidate",
        );
    }

    /// **The run's own calling-loop settings are what the loop runs under.**
    ///
    /// The pass cap is the one loop setting whose effect a fixture can see without a locus
    /// contrived not to settle: at a cap of one the loop reports one pass, and at the shipped
    /// cap of fifty it reports however many it took. A `call_cohort` that reached for the
    /// shipped configuration would report the same number at both.
    #[test]
    fn the_runs_own_calling_loop_settings_are_what_the_loop_runs_under() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);
        let paths = [zeta, mu];

        let at_the_shipped_cap = open_over_calling_with(
            &paths,
            &reference,
            CandidateSelectionConfig::DEFAULT,
            shipped_calling_loop(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");
        let at_a_cap_of_one = open_over_calling_with(
            &paths,
            &reference,
            CandidateSelectionConfig::DEFAULT,
            CallingLoopConfig {
                max_passes: NonZeroU32::new(1).expect("not zero"),
                ..CallingLoopConfig::DEFAULT
            }
            .validate()
            .expect("one pass is a runnable setting"),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        assert_eq!(
            at_a_cap_of_one.called_loci[0].passes, 1,
            "a cap of one pass stops the loop after one",
        );
        assert!(
            at_the_shipped_cap.called_loci[0].passes > 1,
            "and the shipped cap of fifty lets this locus take the {} it needs — a fixture \
             that settled in one pass could not tell the two caps apart",
            at_the_shipped_cap.called_loci[0].passes,
        );
    }

    /// **What each sample's walk counted survives the merge that consumed its walker.**
    ///
    /// The observation cache owns the walkers for the whole merge, so without a route back
    /// every one of these numbers is dropped where the merge returns — and a run report has
    /// nothing to say about a sample beyond its genotypes.
    ///
    /// The three samples are checked **by name and in the run's order**, because tallies
    /// carrying the right numbers under the wrong sample is exactly the failure that reads as
    /// correct.
    #[test]
    fn what_each_walk_counted_comes_back_from_the_merge() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing_on_reads("zeta", "zeta.bam", &[5], 3);
        let (_alpha_dir, alpha) = sample_showing_on_reads("alpha", "alpha.bam", &[5, 20], 5);
        let (_mu_dir, mu) = sample_showing_on_reads("mu", "mu.bam", &[], 7);

        let called = open_over(&[zeta, alpha, mu], &reference)
            .call_cohort(&the_shipped_genotyper())
            .expect("calls");

        assert_eq!(
            called
                .walk
                .per_sample
                .iter()
                .map(|walk| walk.sample_name.as_str())
                .collect::<Vec<_>>(),
            vec!["zeta", "alpha", "mu"],
            "one entry per sample, in the run's own sample order",
        );
        // **Each sample walked a different number of reads**, so the counts identify the
        // sample they came from and a permutation is visible. With three samples of three
        // reads each — which is every other fixture in this file — reversing the walkers on
        // the way out of the merge pairs zeta's name with mu's counts and nothing can see it.
        assert_eq!(
            called
                .walk
                .per_sample
                .iter()
                .map(|walk| {
                    walk.snp_indel
                        .expect("the fixture's ground is a generic region")
                        .reads_admitted
                })
                .collect::<Vec<_>>(),
            vec![3, 5, 7],
            "each sample's own read count, under its own name",
        );
        for walk in &called.walk.per_sample {
            assert!(
                walk.regions.regions_in > 0,
                "{}'s walk was handed regions",
                walk.sample_name,
            );
            assert_eq!(
                walk.regions.regions_handled
                    + walk.regions.unhandled_not_implemented
                    + walk.regions.unhandled_out_of_scope,
                walk.regions.regions_in,
                "{}'s regions partition exactly",
                walk.sample_name,
            );
            let counted = walk
                .snp_indel
                .expect("the fixture's ground is a generic region, so that generator counted");
            assert!(
                counted.reads_admitted > 0,
                "{}'s three reads were admitted",
                walk.sample_name,
            );
            assert_eq!(
                counted.positions_short_of_cap, 0,
                "the read-hold ceiling cost {} no coverage",
                walk.sample_name,
            );
        }
        assert!(
            matches!(
                called.walk.assembly_check,
                AssemblyCheckOutcome::NothingCouldBeChecked { .. }
            ),
            "the fixture's reference carries no checksums, and the run says so rather than \
             claiming every sample agreed: {:?}",
            called.walk.assembly_check,
        );
    }
}

/// **The sample-order join: three numberings meet at every locus, and a mismatch is a wrong
/// genotype rather than a crash.**
///
/// The three, in the words of the modules that own them:
///
/// - **the merge's**, which holds only the samples that covered the locus, each carrying its
///   own index in the run's order (`CohortObservation::per_sample`);
/// - **the run's**, which is `ReadGroups::read_groups_per_sample`'s first-seen order. Every
///   per-sample list the calling loop is *given* is in it — the evidence, and the model
///   parameters;
/// - **the calling scratch's rows**, which are the run's samples with the uncallable ones
///   closed up. The loop's own working buffers are indexed by these, so a row index is
///   neither of the other two, and the scratch is what maps between them
///   (`CallingScratch::claim_row_for`).
///
/// **A per-sample parameter is what makes a swap visible.** Nothing about a call says which
/// sample it belongs to, so a permutation of the run's samples produces a well-formed answer
/// for every one of them. What it changes is *which* answer, and only where the samples are
/// scored under something of their own: here, each sample's inbreeding coefficient, which is
/// how much the prior expects homozygotes.
///
/// **What the coefficient moves at this depth is the confidence, not the genotype, and the
/// tests compare the whole call for that reason.** At the fixture's base quality of 30, two
/// alternative reads cost a homozygous-reference genotype about a millionth, so the reads
/// decide the heterozygote and the prior only moves how sure the caller is: measured, the same
/// four reads give `0/1` at **55.4 Phred** under an outbred coefficient and `0/1` at **33.4**
/// under a nearly fully inbred one. A test comparing genotypes alone would have been blind to
/// a parameters list joined by the merge's entry — measured, that defect leaves both samples
/// at `0/1` and 55.450. It would **not** have been blind to a wrongly joined *evidence* list,
/// which moves the genotypes; the genotype assertions here are what guard that.
#[cfg(test)]
mod the_sample_order_join {
    use super::calling_joined_to_the_merge::the_shipped_genotyper;
    use super::the_merge_over_walkers::{
        open_over_declaring_inbreeding, sample_carrying,
        sample_with_no_reads_in_the_analysed_ground,
    };
    use crate::ng::calling::SampleGenotypeCall;
    use crate::ng::calling::allele_candidates::{CandidateSelectionConfig, MaxCandidateAlleles};
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::test_fixtures::fixture_reference_from_its_index;
    use crate::ng::types::{ContigId, GenomeRegion, Genotype, InbreedingF, Phred, Position};
    use std::path::PathBuf;
    use tempfile::TempDir;

    /// A coefficient at the outbred end — the prior expects Hardy–Weinberg proportions. It is
    /// also the shipped default, so a test that means to *move* a sample must not use it.
    fn outbred() -> InbreedingF {
        InbreedingF::try_new(0.0).expect("zero is a coefficient")
    }

    /// A coefficient at the inbred end — the prior expects almost no heterozygotes. Below one,
    /// which [`InbreedingF`] excludes.
    fn nearly_fully_inbred() -> InbreedingF {
        InbreedingF::try_new(0.99).expect("0.99 is a coefficient")
    }

    /// The one position this module's cohort varies at.
    fn the_locus() -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(15),
            end: Position(15),
        }
    }

    /// **Four samples whose run order, merge order and scratch rows are three different
    /// numberings**, returned with their temporary directories so the files outlive the call.
    ///
    /// | run index | sample | at `chr1:15` | merge entry | scratch row |
    /// |---|---|---|---|---|
    /// | 0 | `zeta` | two reads show `C` | 0 | 0 |
    /// | 1 | `nu` | two reads show `G` | 1 | — *uncallable* |
    /// | 2 | `alpha` | no reads in the analysed ground | — | 1 |
    /// | 3 | `mu` | two reads show `C` | 2 | 2 |
    ///
    /// **`nu` is what makes the third column differ from the second.** Its `G` is the cohort's
    /// lower-ranked alternative — two reads against the four behind `C` — so at a candidate cap
    /// of two, the reference and one alternative, selection cuts it and rules `nu` uncallable
    /// for having earned a sequence the cap removed. The scratch then holds a row for every
    /// sample but `nu`, so `alpha` and `mu` sit one row below their run index.
    ///
    /// **`alpha` is what makes the second column differ from the first**, by covering nothing
    /// where the others cover.
    ///
    /// **`zeta` and `mu` bring identical reads**, so any difference between their calls is a
    /// difference in what they were scored under and nothing else. **The names are not in
    /// alphabetical order**, and the hazard that guards against is real rather than
    /// theoretical: `DeclaredInbreeding` holds its per-sample values in a `BTreeMap`, so a
    /// defect that zipped that map's key order onto the run's samples would hand `zeta` what
    /// `mu` was declared.
    fn four_samples_that_do_not_line_up() -> (Vec<TempDir>, Vec<PathBuf>) {
        let (zeta_dir, zeta) = sample_carrying("zeta", "zeta.bam", b'C');
        let (nu_dir, nu) = sample_carrying("nu", "nu.bam", b'G');
        let (alpha_dir, alpha) = sample_with_no_reads_in_the_analysed_ground("alpha", "alpha.bam");
        let (mu_dir, mu) = sample_carrying("mu", "mu.bam", b'C');
        (
            vec![zeta_dir, nu_dir, alpha_dir, mu_dir],
            vec![zeta, nu, alpha, mu],
        )
    }

    /// The candidate cap that cuts the cohort's lower-ranked alternative: the reference and one
    /// alternative, which is the smallest cap the type allows.
    fn a_cap_of_one_alternative() -> CandidateSelectionConfig {
        CandidateSelectionConfig {
            max_candidate_alleles: MaxCandidateAlleles::new_or_panic(2),
            ..CandidateSelectionConfig::DEFAULT
        }
    }

    /// **The fixture really is three different numberings**, checked rather than described.
    ///
    /// **Without it every other test here could be passing on the identity permutation.** If
    /// `alpha` covered the locus, the merge's samples would be the run's own order; if `nu`
    /// stayed callable, the scratch's rows would be too. Either way a run that indexed one
    /// list by another would look correct in every fixture in this module. What this asserts
    /// is the shape that makes a swap visible: the merge holds **three** entries naming run
    /// samples **0, 1 and 3**, and exactly one sample comes back with no genotype, so the
    /// scratch's rows are **three** where the run has four.
    #[test]
    fn the_three_numberings_are_three_different_numberings() {
        let (_dirs, paths) = four_samples_that_do_not_line_up();
        let (_reference_dir, reference) = fixture_reference_from_its_index();

        let merged = open_over_declaring_inbreeding(
            &paths,
            &reference,
            &DeclaredInbreeding::nothing_said(),
            a_cap_of_one_alternative(),
        )
        .merge_cohort()
        .expect("merges");
        assert_eq!(merged.cohort_observations.len(), 1);
        assert_eq!(merged.cohort_observations[0].region, the_locus());
        assert_eq!(
            merged.cohort_observations[0]
                .per_sample
                .iter()
                .map(|support| support.sample)
                .collect::<Vec<_>>(),
            vec![0, 1, 3],
            "zeta, nu and mu covered the locus and alpha did not, so the merge's three entries \
             name the run's first, second and fourth samples",
        );

        let called = open_over_declaring_inbreeding(
            &paths,
            &reference,
            &DeclaredInbreeding::nothing_said(),
            a_cap_of_one_alternative(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");
        let calls = &called.called_loci[0].per_sample;
        assert_eq!(calls.len(), 4, "one call per sample of the run");
        assert_eq!(
            calls
                .iter()
                .map(SampleGenotypeCall::is_missing)
                .collect::<Vec<_>>(),
            vec![false, true, false, false],
            "nu earned the G the cap cut, so it alone is set aside — which is what makes the \
             scratch's three rows a different numbering from the run's four samples",
        );
    }

    /// **Each sample is scored under its own coefficient, and the fixture is built so that
    /// scoring it under a neighbour's would show.**
    ///
    /// `zeta` and `mu` bring identical reads to the locus and are declared at opposite ends of
    /// the coefficient's range, so their calls must differ. Between them in the run's order
    /// sit `nu`, which the cap sets aside, and `alpha`, which covers nothing — so `mu` is the
    /// run's sample 3, the merge's entry 2 and the scratch's row 2, and a run that reached for
    /// its coefficient by either of the other two numberings would score it under something
    /// else.
    #[test]
    fn each_samples_own_inbreeding_coefficient_is_what_scores_it() {
        let (_dirs, paths) = four_samples_that_do_not_line_up();
        let (_reference_dir, reference) = fixture_reference_from_its_index();

        let called = open_over_declaring_inbreeding(
            &paths,
            &reference,
            &DeclaredInbreeding::nothing_said()
                .and_this_sample("zeta", outbred())
                .and_this_sample("mu", nearly_fully_inbred()),
            a_cap_of_one_alternative(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        assert_eq!(called.called_loci.len(), 1);
        assert_eq!(called.called_loci[0].region, the_locus());
        let calls = &called.called_loci[0].per_sample;
        assert_ne!(
            calls[0], calls[3],
            "zeta and mu brought the same two-and-two reads and were declared at opposite ends \
             of the coefficient's range, so the only thing that can separate their calls is \
             the coefficient — and it must",
        );
        // **The direction is checked as well as the difference.** A join that handed each
        // sample *some* coefficient other than its own would also make these two calls
        // differ; what says the right one arrived is which of them is the less confident, and
        // the sample the prior pushes away from a heterozygote is the inbred one.
        let (zeta, mu) = (quality_of(&calls[0]), quality_of(&calls[3]));
        assert!(
            mu < zeta,
            "mu is declared nearly fully inbred and called a heterozygote, so it must be the \
             less confident of the two: zeta {zeta:?} Phred, mu {mu:?} Phred",
        );
        // **This one is not a restatement of "the genotype does not move" — it is the guard
        // on the *evidence* join.** Measured: a sample scored with no reads of its own comes
        // back `0/1` at 2.2 Phred, and one scored on somebody else's four reads comes back at
        // 55.4, so if the shaping put the merge's entries on the wrong run rows the two
        // genotypes would part company here. Deleting it because the qualities already differ
        // loses the only assertion that sees that defect.
        assert_eq!(
            genotype_of(&calls[0]),
            genotype_of(&calls[3]),
            "zeta and mu must be called the same genotype — four reads at Q30 decide the \
             heterozygote at both ends of the coefficient's range, so a difference here is \
             the evidence reaching the wrong sample rather than the prior doing its work",
        );
    }

    /// **Swapping the two coefficients swaps the two calls** — the oracle the plan asks for.
    ///
    /// Nothing else about the run moves: the same four files, the same reads, the same order,
    /// the same cap. Only which sample each declared coefficient names. If the run joined the
    /// parameters to the samples by anything but the run's own order — the merge's entry, the
    /// scratch's row, the order the coefficients were declared in, a sort of the names — the
    /// two calls would not exchange.
    #[test]
    fn swapping_two_samples_coefficients_swaps_their_calls() {
        let (_dirs, paths) = four_samples_that_do_not_line_up();
        let (_reference_dir, reference) = fixture_reference_from_its_index();

        let called_with = |zeta_coefficient, mu_coefficient| {
            open_over_declaring_inbreeding(
                &paths,
                &reference,
                &DeclaredInbreeding::nothing_said()
                    .and_this_sample("zeta", zeta_coefficient)
                    .and_this_sample("mu", mu_coefficient),
                a_cap_of_one_alternative(),
            )
            .call_cohort(&the_shipped_genotyper())
            .expect("calls")
        };
        let one_way = called_with(outbred(), nearly_fully_inbred());
        let the_other = called_with(nearly_fully_inbred(), outbred());

        let (first, second) = (
            &one_way.called_loci[0].per_sample,
            &the_other.called_loci[0].per_sample,
        );
        assert_eq!(
            first[0], second[3],
            "zeta declared outbred must call exactly as mu does when mu is declared outbred",
        );
        assert_eq!(first[3], second[0], "and the inbred end likewise",);
        assert_ne!(
            first[0], first[3],
            "a fixture whose two calls were equal would satisfy the two assertions above \
             without the coefficients reaching anybody",
        );
        // The two samples between them were named by neither declaration, so the swap must
        // leave both alone — including the one the cap set aside, which has no call to move.
        assert_eq!(first[1], second[1], "nu is set aside under both");
        assert_eq!(
            first[2], second[2],
            "and alpha was not named by either declaration",
        );
    }

    /// **The sample that covers nothing is still called, in its own place in the run's order,
    /// by the prior alone.**
    ///
    /// `alpha` has no reads in the analysed ground, so the merge holds no entry for it — and
    /// the loop reads one entry per sample of the *run*. An empty sum is zero, so every
    /// genotype scores alike and the prior decides on its own, which is the right answer
    /// rather than a special case. Measured: it comes back `0/1` at **2.2 Phred**, a call the
    /// caller is almost entirely unsure of, which is what "the prior alone" produces and what
    /// distinguishes it from a sample the candidate step set aside.
    #[test]
    fn a_sample_that_covered_nothing_is_called_by_the_prior_alone() {
        let (_dirs, paths) = four_samples_that_do_not_line_up();
        let (_reference_dir, reference) = fixture_reference_from_its_index();

        let called = open_over_declaring_inbreeding(
            &paths,
            &reference,
            &DeclaredInbreeding::nothing_said(),
            a_cap_of_one_alternative(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        let alpha = &called.called_loci[0].per_sample[2];
        assert!(
            !alpha.is_missing(),
            "a sample with no reads here is scored, not set aside — set aside is what the \
             candidate step does to a sample whose own allele the cap cut, which is nu",
        );
        let confidence = quality_of(alpha).expect("a scored sample carries a quality");
        assert!(
            confidence < 5.0,
            "with no reads of its own the prior decides alone and the call is barely held: \
             {confidence} Phred",
        );
        // **And the call says so itself**, which is what emission needs and what nothing
        // downstream could work out: the likelihoods live in scratch the next locus overwrites.
        // A sample with no reads has a flat row; one with reads does not.
        assert_eq!(
            called.called_loci[0]
                .per_sample
                .iter()
                .map(reads_said_nothing)
                .collect::<Vec<_>>(),
            vec![Some(false), None, Some(true), Some(false)],
            "alpha's reads said nothing because it has none; zeta and mu covered the locus and \
             theirs did; nu was set aside by the cap and so has no call to carry the fact at \
             all — which is the second, different route to a `./.`",
        );
        // **The whole list, not the one entry.** Index 2 of four is a fixed point under
        // several ways of getting the order wrong — reversal about the middle, and a zip that
        // drops a sample and leaves the rest shifted — so asserting one name passes on a list
        // that is wrong in both its length and its pairing.
        assert_eq!(
            called
                .walk
                .per_sample
                .iter()
                .map(|walk| walk.sample_name.as_str())
                .collect::<Vec<_>>(),
            vec!["zeta", "nu", "alpha", "mu"],
            "one entry per sample, in the run's own sample order",
        );
    }

    /// A call's genotype, or `None` where the sample was set aside.
    fn genotype_of(call: &SampleGenotypeCall) -> Option<&Genotype> {
        call.genotype()
    }

    /// Whether the loop found this sample's own reads said nothing about its genotype, or
    /// `None` where the sample was set aside — which is a different fact and a different `./.`.
    fn reads_said_nothing(call: &SampleGenotypeCall) -> Option<bool> {
        match call {
            SampleGenotypeCall::Called {
                reads_were_uninformative,
                ..
            } => Some(*reads_were_uninformative),
            SampleGenotypeCall::Missing => None,
        }
    }

    /// How sure the caller is, or `None` where the sample was set aside.
    ///
    /// **Compared as an `f32` and not as an `Option`**, wherever the comparison is what a test
    /// rests on: `None < Some(_)`, so a sample that went missing would slip past an ordering
    /// assertion written on the options and print `None` where the message says Phred.
    fn quality_of(call: &SampleGenotypeCall) -> Option<f32> {
        call.score_best_genotype().map(Phred::get)
    }
}

/// **A locus nobody can be called at is counted and reported, and never ends the run.**
///
/// The candidate step cuts an allele rather than refusing a locus, on the ground that most
/// samples stay callable (`doc/devel/ng/arch/candidate_alleles.md` §4.1). Where the cap cuts a
/// sequence that **every** covering sample had reads on, none of them is callable, and there is
/// no genotype this caller could honestly report for anybody there.
///
/// **That was a panic until 2026-09-01**, which meant one hard locus ended a cohort's run after
/// however many hours of walking. The owner's ruling is that such a locus is counted and
/// reported at the end; these tests are what says it is.
#[cfg(test)]
mod a_locus_nobody_can_be_called_at {
    use super::calling_joined_to_the_merge::the_shipped_genotyper;
    use super::the_merge_over_walkers::{
        open_over_declaring_inbreeding, sample_carrying, sample_carrying_two_alternatives,
    };
    use crate::ng::calling::allele_candidates::{CandidateSelectionConfig, MaxCandidateAlleles};
    use crate::ng::calling::parameters_file::DeclaredInbreeding;
    use crate::ng::read::input::test_fixtures::fixture_reference_from_its_index;
    use crate::ng::types::{ContigId, GenomeRegion, Position};

    /// The cap that keeps the reference and one alternative — the smallest the type allows, and
    /// what makes a second alternative something a sample can have earned and lost.
    fn a_cap_of_one_alternative() -> CandidateSelectionConfig {
        CandidateSelectionConfig {
            max_candidate_alleles: MaxCandidateAlleles::new_or_panic(2),
            ..CandidateSelectionConfig::DEFAULT
        }
    }

    /// **Two samples that both carry both alternatives: the run finishes, and says where it
    /// could call nobody.**
    ///
    /// Each sample shows `C` on two reads, `G` on two and the reference on two, so both
    /// alternatives clear the merge's floor and both are earned by both samples. At a cap of
    /// one alternative the lower-ranked is cut, every covering sample has lost a sequence its
    /// own reads earned, and the locus has nobody to call.
    #[test]
    fn a_locus_where_the_cap_cut_everybodys_allele_is_counted_and_the_run_finishes() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) =
            sample_carrying_two_alternatives("zeta", "zeta.bam", b'C', b'G', &[5]);
        let (_mu_dir, mu) = sample_carrying_two_alternatives("mu", "mu.bam", b'C', b'G', &[5]);

        let called = open_over_declaring_inbreeding(
            &[zeta, mu],
            &reference,
            &DeclaredInbreeding::nothing_said(),
            a_cap_of_one_alternative(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("a locus nobody can be called at does not end the run");

        assert_eq!(
            called.loci_with_nobody_to_call,
            vec![GenomeRegion {
                contig: ContigId(0),
                start: Position(15),
                end: Position(15),
            }],
            "the one locus, reported with its ground so a person can go and look at it",
        );
        assert!(
            called.called_loci.is_empty(),
            "and no record is made for it: {:?}",
            called.called_loci,
        );
        // **The three lists are three facts, which the type's own documentation claims and
        // nothing else here checks.** A locus counted in two of them would be reported twice
        // and would make the loci-assembled total wrong in both directions.
        assert!(
            called.loci_too_wide_to_assemble.is_empty(),
            "this locus was assembled — the width bound refused nothing: {:?}",
            called.loci_too_wide_to_assemble,
        );
    }

    /// **Two such loci come back in genome order**, which is what the field claims.
    ///
    /// The same two samples, each carrying both alternatives at **two** positions, so the cap
    /// leaves nobody callable at both. One locus could not tell an ordered list from a
    /// reversed one.
    #[test]
    fn the_loci_nobody_can_be_called_at_come_back_in_genome_order() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) =
            sample_carrying_two_alternatives("zeta", "zeta.bam", b'C', b'G', &[5, 20]);
        let (_mu_dir, mu) = sample_carrying_two_alternatives("mu", "mu.bam", b'C', b'G', &[5, 20]);

        let called = open_over_declaring_inbreeding(
            &[zeta, mu],
            &reference,
            &DeclaredInbreeding::nothing_said(),
            a_cap_of_one_alternative(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("neither locus ends the run");

        assert_eq!(
            called
                .loci_with_nobody_to_call
                .iter()
                .map(|region| region.start)
                .collect::<Vec<_>>(),
            vec![Position(15), Position(30)],
            "both loci, in genome order",
        );
    }

    /// **A cohort of one sample is called, not emptied.**
    ///
    /// The guard is "no sample of the run is callable", and at one sample that is one sample —
    /// so a guard written `<= 1` rather than `== 0` would drop **every** locus of a
    /// single-sample run and hand back an empty output with a count beside it. A single sample
    /// is the thinnest end of the range this caller commits to
    /// (`doc/devel/ng/spec/design_principles.md` §0) and the one no other fixture in this file
    /// exercises.
    #[test]
    fn a_cohort_of_one_sample_is_called() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_carrying("zeta", "zeta.bam", b'C');

        let called = open_over_declaring_inbreeding(
            std::slice::from_ref(&zeta),
            &reference,
            &DeclaredInbreeding::nothing_said(),
            CandidateSelectionConfig::DEFAULT,
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        assert!(
            called.loci_with_nobody_to_call.is_empty(),
            "one callable sample is somebody to call: {:?}",
            called.loci_with_nobody_to_call,
        );
        assert_eq!(called.called_loci.len(), 1, "and the locus is called");
        assert_eq!(
            called.called_loci[0].per_sample.len(),
            1,
            "over the run's one sample",
        );
    }

    /// **The same cohort at the shipped cap calls that locus normally**, which is what says the
    /// test above is about the cap rather than about the reads.
    ///
    /// At six candidate alleles nothing is cut, nobody has lost an allele they earned, and both
    /// samples are called.
    #[test]
    fn the_same_cohort_at_the_shipped_cap_calls_that_locus() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) =
            sample_carrying_two_alternatives("zeta", "zeta.bam", b'C', b'G', &[5]);
        let (_mu_dir, mu) = sample_carrying_two_alternatives("mu", "mu.bam", b'C', b'G', &[5]);

        let called = open_over_declaring_inbreeding(
            &[zeta, mu],
            &reference,
            &DeclaredInbreeding::nothing_said(),
            CandidateSelectionConfig::DEFAULT,
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        assert!(
            called.loci_with_nobody_to_call.is_empty(),
            "nothing was cut, so nobody lost an allele they earned: {:?}",
            called.loci_with_nobody_to_call,
        );
        assert_eq!(called.called_loci.len(), 1, "the locus is called");
        assert_eq!(
            called.called_loci[0].alleles().len(),
            3,
            "over the reference and both alternatives",
        );
    }

    /// **A cohort where only some samples lose their allele is called, not counted**, so the
    /// new list is the *nobody* case and not a rename of the set-aside one.
    ///
    /// `nu` alone carries the lower-ranked `G`; `zeta` and `mu` carry the `C` that survives.
    /// The cap sets `nu` aside and the other two are called, which is the ordinary behaviour
    /// the ruling above does not change.
    #[test]
    fn a_locus_where_only_some_samples_lose_their_allele_is_still_called() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_carrying("zeta", "zeta.bam", b'C');
        let (_nu_dir, nu) = sample_carrying("nu", "nu.bam", b'G');
        let (_mu_dir, mu) = sample_carrying("mu", "mu.bam", b'C');

        let called = open_over_declaring_inbreeding(
            &[zeta, nu, mu],
            &reference,
            &DeclaredInbreeding::nothing_said(),
            a_cap_of_one_alternative(),
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("calls");

        assert!(
            called.loci_with_nobody_to_call.is_empty(),
            "two of the three samples keep the allele they earned, so there is somebody to call",
        );
        assert_eq!(called.called_loci.len(), 1);
        assert_eq!(
            called.called_loci[0]
                .per_sample
                .iter()
                .map(|call| call.is_missing())
                .collect::<Vec<_>>(),
            vec![false, true, false],
            "nu alone is set aside — which is a different fact from the locus having nobody",
        );
    }
}

#[cfg(test)]
mod records_handed_over_as_the_run_finishes_them {
    //! **The run's own path: records leave one at a time and none is kept.**
    //!
    //! [`AlignedFilesVariantCaller::call_cohort`] is the oracle here — same reads, same
    //! genotyper, same settings — and what these check is the difference: that the records
    //! describe exactly the loci it calls minus the ones that establish no variant, that they
    //! arrive in genome order, and that a run whose output refuses a record stops and says
    //! where.

    use super::calling_joined_to_the_merge::the_shipped_genotyper;
    use super::the_merge_over_walkers::{
        open_over, read_of, sample_showing, sample_showing_on_reads, shipped_calling_loop,
    };
    use super::*;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference_from_its_index, header, indexed_named_bam, matching_contigs,
    };
    use crate::ng::types::{ContigId, Ploidy, Position};
    use crate::ng::vcf::header::{HeaderContig, VcfHeaderMetadata};
    use crate::ng::vcf::writer::VcfWriter;
    use crate::ng::vcf::{SampleCall, VcfRecord};
    use noodles_sam::alignment::RecordBuf;
    use std::path::PathBuf;
    use tempfile::TempDir;

    /// Collect every record a run hands over, and the counts beside them.
    fn records_of(caller: AlignedFilesVariantCaller) -> (Vec<VcfRecord>, WrittenCohort) {
        let mut records = Vec::new();
        let written = caller
            .call_cohort_handing_each_record_over(&the_shipped_genotyper(), &mut |record| {
                records.push(record.clone());
                Ok::<(), std::io::Error>(())
            })
            .expect("the fixture cohort calls and every record is taken");
        (records, written)
    }

    /// **The records are the loci `call_cohort` calls, less the ones that establish no
    /// variant** — same spans, same genotypes, same order.
    ///
    /// The oracle is the entry point every Milestone D test is written against, so a record
    /// path that called differently — a scratch it reset differently, a sample order it read
    /// differently — shows here as a genotype that does not match.
    #[test]
    fn the_records_describe_the_loci_call_cohort_calls() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let (_mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);
        let paths = [zeta, alpha, mu];

        let called = open_over(&paths, &reference)
            .call_cohort(&the_shipped_genotyper())
            .expect("calls");
        let (records, written) = records_of(open_over(&paths, &reference));

        assert_eq!(
            records.iter().map(VcfRecord::region).collect::<Vec<_>>(),
            called
                .called_loci
                .iter()
                .map(|locus| locus.region)
                .collect::<Vec<_>>(),
            "every locus of this fixture carries an alternative somebody was called on, so the \
             two lists are the same spans in the same order",
        );
        assert_eq!(written.records_written, records.len() as u64);
        assert_eq!(written.loci_called(), called.called_loci.len() as u64);
        assert_eq!(
            written.loci_called_but_not_written, 0,
            "nothing here was called all-reference",
        );

        for (record, locus) in records.iter().zip(&called.called_loci) {
            assert_eq!(
                record
                    .sample_columns()
                    .iter()
                    .map(|column| matches!(column.call, SampleCall::Called { .. }))
                    .collect::<Vec<_>>(),
                locus
                    .per_sample
                    .iter()
                    .map(|call| !call.is_missing())
                    .collect::<Vec<_>>(),
                "the samples the file writes a genotype for at {}",
                record.region(),
            );
        }
    }

    /// A sample carrying `first` on one read and `second` on another at `chr1:15`, with four
    /// reads matching the reference — thirty bases from `chr1:10`.
    ///
    /// **The shape of a locus that is built and then establishes nothing.** The merge keeps a
    /// position on the cohort's *pooled* non-reference reads, which two of them reach;
    /// candidate selection then asks each sequence separately, and one read apiece is below its
    /// floor of two, so both are dropped and the locus is called over the reference alone.
    /// Measured at more than one built locus in four on both benchmarks
    /// (`doc/devel/ng/spec/candidate_alleles.md` §6.2), so it is the ordinary case.
    fn sample_showing_two_sequences_once_each(
        sample: &str,
        file_name: &str,
        first: u8,
        second: u8,
    ) -> (TempDir, PathBuf) {
        let with = |alt: u8| {
            let mut bases = [b'A'; 30];
            bases[5] = alt;
            bases
        };
        let records: Vec<RecordBuf> = vec![
            read_of(&format!("{sample}-a"), 10, &with(first)),
            read_of(&format!("{sample}-b"), 10, &with(second)),
            read_of(&format!("{sample}-r0"), 10, &[b'A'; 30]),
            read_of(&format!("{sample}-r1"), 10, &[b'A'; 30]),
            read_of(&format!("{sample}-r2"), 10, &[b'A'; 30]),
            read_of(&format!("{sample}-r3"), 10, &[b'A'; 30]),
        ];
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
            file_name,
        )
    }

    /// **A locus called over the reference alone establishes no variant and is not written**,
    /// and the run counts it rather than losing it (spec §9).
    ///
    /// There is no gVCF and no reference block, so a record's absence is the file saying
    /// *nothing here*. A run that wrote this locus would emit `ALT .` with every sample `0/0`,
    /// which spec §5 admits only for a filtered repeat-tract record.
    #[test]
    fn a_locus_called_over_the_reference_alone_is_counted_and_not_written() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) =
            sample_showing_two_sequences_once_each("zeta", "zeta.bam", b'C', b'G');

        let called = open_over(std::slice::from_ref(&zeta), &reference)
            .call_cohort(&the_shipped_genotyper())
            .expect("calls");
        let (records, written) = records_of(open_over(std::slice::from_ref(&zeta), &reference));

        assert_eq!(
            called.called_loci.len(),
            1,
            "the merge builds the position on the cohort's two pooled non-reference reads",
        );
        assert_eq!(
            called.called_loci[0].alleles().len(),
            1,
            "and candidate selection leaves the reference alone, each sequence having one read \
             against a floor of two",
        );
        assert!(records.is_empty(), "so nothing is written");
        assert_eq!(written.records_written, 0);
        assert_eq!(
            written.loci_called_but_not_written, 1,
            "the locus is counted, not lost: called and establishing nothing",
        );
        assert_eq!(written.loci_called(), 1);
    }

    /// **A cohort of one sample writes its records**, which is the end of the range the caller
    /// commits to (`doc/devel/ng/spec/design_principles.md` §0) and the shape most likely to
    /// have a guard written `<= 1` where it meant `== 0`.
    #[test]
    fn a_cohort_of_one_sample_writes_its_records() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5, 20]);

        let (records, written) = records_of(open_over(std::slice::from_ref(&zeta), &reference));

        assert_eq!(
            written.records_written, 2,
            "the two positions zeta varies at"
        );
        assert_eq!(records.len(), 2);
        for record in &records {
            assert_eq!(
                record.sample_columns().len(),
                1,
                "one column, for the run's one sample",
            );
        }
    }

    /// **A run whose output refuses a record stops there and names the locus it was writing.**
    ///
    /// A full disk is the ordinary way to reach this, and the cause — an `io::Error` — says
    /// nothing about where the file stopped. A consumer holding the two knows how much of its
    /// output is complete.
    #[test]
    fn an_output_that_refuses_a_record_stops_the_run_naming_the_locus() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5, 20]);

        let mut taken = 0;
        let stopped = open_over(std::slice::from_ref(&zeta), &reference)
            .call_cohort_handing_each_record_over(&the_shipped_genotyper(), &mut |_record| {
                taken += 1;
                Err(std::io::Error::other("the disk is full"))
            })
            .expect_err("a refused record ends the run");

        assert_eq!(taken, 1, "the sink is not called again after it refuses");
        match stopped {
            RunError::RecordNotWritten { locus, source } => {
                assert_eq!(
                    locus,
                    GenomeRegion {
                        contig: ContigId(0),
                        start: Position(15),
                        end: Position(15),
                    },
                    "the first of the two positions, which is where it stopped",
                );
                assert!(source.to_string().contains("the disk is full"));
            }
            other => panic!("expected the record refusal, got {other:?}"),
        }
    }

    /// A sample whose reads carry a `bases`-long insertion after `chr1:14`, thirty reference
    /// bases from `chr1:10`.
    ///
    /// **An insertion rather than a deletion, and the fixture reference is why.** Every module
    /// here shares a reference of a hundred `A`s, so a deletion inside it sits in one
    /// homopolymer and left-alignment slides it off the record — measured at D1, where three
    /// reads carrying a five-base deletion produced no cohort locus at all. An insertion of
    /// `C`s introduces bases the reference does not have, so it survives left-alignment and
    /// reaches the merge.
    fn sample_carrying_an_insertion(
        sample: &str,
        file_name: &str,
        inserted: usize,
    ) -> (TempDir, PathBuf) {
        use noodles_sam::alignment::record::cigar::op::{Kind, Op};
        use noodles_sam::alignment::record_buf::Sequence;
        let carrying = |name: &str| {
            let mut record = read_of(name, 10, &[b'A'; 30]);
            let mut sequence = vec![b'A'; 5];
            sequence.extend(std::iter::repeat_n(b'C', inserted));
            sequence.extend(std::iter::repeat_n(b'A', 25));
            *record.sequence_mut() = Sequence::from(sequence.clone());
            *record.cigar_mut() = [
                Op::new(Kind::Match, 5),
                Op::new(Kind::Insertion, inserted),
                Op::new(Kind::Match, 25),
            ]
            .into_iter()
            .collect();
            *record.quality_scores_mut() =
                noodles_sam::alignment::record_buf::QualityScores::from(vec![
                    30_u8;
                    sequence.len()
                ]);
            record
        };
        let records: Vec<RecordBuf> = (0..3)
            .map(|read| carrying(&format!("{sample}-ins{read}")))
            .collect();
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
            file_name,
        )
    }

    /// **An insertion goes through the whole path and its record needs no padding base.**
    ///
    /// This is the claim three documents now rest on, as a test rather than as prose: the
    /// generic mint anchors its indels — an insertion's reference span is its anchor base alone
    /// (`ReadEvent::footprint_span`) — so the record is `REF A` against `ALT ACC`, **no allele
    /// is empty**, and `POS` does not move. A mint that instead emitted the inserted bases with
    /// an empty reference would need a padding base here, and `VcfRecord::new` would refuse the
    /// record without one.
    ///
    /// It is also the only fixture in this module whose reads are not all substitutions, which
    /// is what the correctness review found missing: discarding the fetched padding base passed
    /// every other test in the crate.
    #[test]
    fn an_insertion_is_written_as_a_record_that_needs_no_padding_base() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_carrying_an_insertion("zeta", "zeta.bam", 2);

        let (records, written) = records_of(open_over(std::slice::from_ref(&zeta), &reference));

        assert_eq!(
            written.records_written,
            1,
            "the insertion is one record; got {} record(s) at {:?}",
            written.records_written,
            records.iter().map(VcfRecord::region).collect::<Vec<_>>(),
        );
        let record = &records[0];
        assert_eq!(
            record.padding_base(),
            None,
            "the mint anchors the insertion, so the reference allele spells its anchor base and \
             no allele of the record is empty",
        );
        assert_eq!(
            record.alleles()[0].as_ref(),
            b"A",
            "REF is the anchor base alone — the insertion's reference span is 1",
        );
        assert!(
            record
                .alternatives()
                .iter()
                .any(|allele| allele.as_ref() == b"ACC"),
            "and the alternative is the anchor plus the two inserted bases, got {:?}",
            record
                .alternatives()
                .iter()
                .map(|allele| String::from_utf8_lossy(allele).into_owned())
                .collect::<Vec<_>>(),
        );
    }

    /// **A run's records go through the real writer and come back as a readable VCF.**
    ///
    /// The artefact a person sees, end to end: the header the run states, one line per record,
    /// the samples in the run's own order. What this catches that the record assertions above
    /// cannot is a record the writer refuses — an order it will not accept, a column it cannot
    /// encode — which is exactly what a run discovers only when it writes.
    #[test]
    fn a_runs_records_become_a_vcf_a_reader_can_parse() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_zeta_dir, zeta) = sample_showing("zeta", "zeta.bam", &[5]);
        let (_alpha_dir, alpha) = sample_showing("alpha", "alpha.bam", &[5, 20]);
        let paths = [zeta, alpha];
        let caller = open_over(&paths, &reference);
        let sample_names: Vec<String> = caller.sample_names().map(str::to_owned).collect();
        let metadata = VcfHeaderMetadata::try_new(
            vec![HeaderContig {
                name: "chr1".to_owned(),
                length: 100,
                md5: None,
            }],
            sample_names.clone(),
            "pop_var_caller_exp call-from-alignments".to_owned(),
            "reference.fa".to_owned(),
            String::new(),
        )
        .expect("a header this run can state");

        let out = tempfile::tempdir().expect("a temporary directory");
        let path = out.path().join("calls.vcf");
        let mut writer = VcfWriter::create(&path, metadata, Ploidy::try_new(2).expect("a diploid"))
            .expect("the output opens");
        let written = caller
            .call_cohort_handing_each_record_over(&the_shipped_genotyper(), &mut |record| {
                writer.write_record(record)
            })
            .expect("calls");
        writer.finish().expect("the file is renamed into place");

        let text = std::fs::read_to_string(&path).expect("the VCF is there");
        let heading = text
            .lines()
            .find(|line| line.starts_with("#CHROM"))
            .expect("a #CHROM line");
        assert!(
            heading.ends_with(&format!("\t{}\t{}", sample_names[0], sample_names[1])),
            "the samples in the run's own order, and got {heading}",
        );
        let records: Vec<&str> = text.lines().filter(|line| !line.starts_with('#')).collect();
        assert_eq!(records.len(), written.records_written as usize);
        assert_eq!(
            records
                .iter()
                .map(|line| line.split('\t').take(2).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
            vec![vec!["chr1", "15"], vec!["chr1", "30"]],
            "the two positions the cohort varies at, by contig name and 1-based position",
        );
        for line in &records {
            let columns: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                columns.len(),
                11,
                "nine fixed columns and two samples: {line}"
            );
            assert_eq!(columns[8], "GT:GQ:DP:AD", "the generic FORMAT string");
        }
    }

    // ---------------------------------------------------------------
    // Milestone E2 — concurrency invariance (spec §12.2; the plan's E2).
    //
    // The record path's whole claim is the same VCF at every worker count, and since E1 a
    // worker count is a rayon thread count: the cover's sweep is the one thing that threads.
    // The oracle is the serial caller, on a cohort whose loci differ in kind — a shared SNP,
    // a sample-private SNP, an insertion, a locus that is called and establishes no variant,
    // a sample with no evidence at all — over ground with a repeat tract interleaved between
    // two ordinary stretches, so the walk crosses every kind of ground the caller routes.
    //
    // Spec §8's calling-scratch trap is what this oracle exists to catch, and today it cannot
    // fire: no built arrangement reorders the loci a scratch sees, because assembly and
    // genotyping stay on the merge thread in genome order. The oracle is built now so that
    // whatever first threads the calling is caught by a test that predates it.
    //
    // Seeded by the E1 determinism review's probe suite, which ran this shape (without the
    // tract, the mixed cohort and the width sweep) at pools of 1, 2, 4 and 8 and found the
    // bytes identical.
    // ---------------------------------------------------------------

    /// The analysed ground E2 calls over: an ordinary stretch, a repeat tract, an ordinary
    /// stretch — so ordinary sites and a tract are interleaved, as the plan's oracle asks.
    ///
    /// The tract routes to an unfilled generator slot and is charged to *not built yet*
    /// (every E2 read lies inside the first stretch), which is exactly what a run over
    /// tract-bearing ground does today; what the fixture pins is that the routing and the
    /// accounting are identical at every thread count, not that a tract is called.
    fn ground_with_a_tract_interleaved() -> Segmentation {
        use crate::ng::region_typing::segment_criteria::SsrSegment;
        use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion};
        use crate::ng::repeat_catalog::StrRepeatCriteria;
        use crate::ng::types::Motif;
        use crate::regions::ContigBounds;

        let chr1 = |start: u64, end: u64| GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        };
        let bounds = [ContigBounds {
            name: "chr1",
            length: 100,
        }];
        let segments = vec![
            TypedRegion {
                region: chr1(1, 40),
                kind: RegionKind::Generic,
            },
            TypedRegion {
                region: chr1(41, 52),
                kind: RegionKind::SsrSegment(
                    SsrSegment::new("chr1".into(), 41, 52, Motif::new(b"AT").unwrap(), 1.0)
                        .expect("a twelve-base AT tract inside the contig"),
                ),
            },
            TypedRegion {
                region: chr1(53, 100),
                kind: RegionKind::Generic,
            },
        ];
        Segmentation::build(
            segments.into_iter().map(Ok),
            GenomeRegions::whole_contigs(&bounds),
            super::tests::catalog_header(),
            StrRepeatCriteria::default(),
            PathBuf::from("/genomes/test.catalog.parquet"),
        )
        .expect("a clean stream builds")
    }

    /// A caller over `paths` on the tract-interleaved ground, with the shipped settings and
    /// the merge's knobs named — E2 sweeps the building-region width so that the run is made
    /// to divide the same ground two different ways while answering the same file. See
    /// [`the_record_path_is_byte_identical_at_every_thread_count`] for what that sweep does
    /// and does not buy; it is not what makes a broken parallel reduction visible, which was
    /// measured rather than assumed.
    fn open_over_the_tract_ground_with(
        paths: &[PathBuf],
        reference: &crate::ng::read::input::reference::OpenReference,
        merge: MergeParameters,
    ) -> AlignedFilesVariantCaller {
        use crate::ng::calling::parameters_file::DeclaredInbreeding;
        use crate::ng::read::input::read_groups::build_read_groups;

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
                read_filters: ReadFilterConfig::default(),
                build_index_if_missing: false,
                locus_generator_settings: PileupGeneratorConfig::default(),
                reference_with_checksums: reference.info(),
            },
            ground_with_a_tract_interleaved(),
            parameters,
            shipped_calling_loop(),
            crate::ng::calling::allele_candidates::CandidateSelectionConfig::DEFAULT,
            merge,
        )
        .expect("five readable samples over a readable reference open")
    }

    /// The two building-region widths E2 runs everything at: the shipped default, where the
    /// fixture's first ordinary stretch is a single building region, and a seven-base width,
    /// which cuts that stretch into six — so the run is made to cover, evict and build its
    /// ground in six steps instead of one while answering the same file.
    fn the_two_widths() -> [MergeParameters; 2] {
        [
            MergeParameters::DEFAULT,
            MergeParameters {
                cohort_locus_builder_regions_len:
                    crate::ng::run::cohort_merge::CohortLocusBuilderRegionsLen(
                        std::num::NonZeroU32::new(7).expect("seven is non-zero"),
                    ),
                ..MergeParameters::DEFAULT
            },
        ]
    }

    /// E2's cohort: every kind of locus the built caller can produce, in one run.
    ///
    /// zeta and alpha carry a shared SNP at `chr1:17` (three and four reads, so the samples
    /// are numerically distinguishable); alpha alone carries `chr1:35`; iota carries a
    /// two-base insertion anchored at `chr1:14`; kappa shows two sequences one read each at
    /// `chr1:15`, which the merge builds on the pooled two and candidate selection then
    /// empties — the called-but-not-written locus; and mu has reads that all match the
    /// reference, so it is genotyped from coverage alone at every locus.
    ///
    /// **⚑ Every locus this cohort produces is one reference position wide**, and that is a
    /// limitation rather than a choice: the fixture reference is a hundred identical bases,
    /// so two substitutions never share a base to chain on, an insertion's reference span is
    /// its anchor alone, and a deletion left-aligns off the record (measured at D1). Two
    /// samples departing at adjacent positions were tried and closed as two separate loci,
    /// not one. `the_mixed_cohorts_records_describe_the_serial_callers_loci` asserts the
    /// one-position width, so the limitation is checked rather than assumed — it is what
    /// stops this module from pinning the parallel cover's chain-following, which
    /// `cohort_merge`'s own fixtures do instead.
    fn the_mixed_cohort() -> (Vec<TempDir>, Vec<PathBuf>) {
        let (zeta_dir, zeta) = sample_showing_on_reads("zeta", "zeta.bam", &[7], 3);
        let (alpha_dir, alpha) = sample_showing_on_reads("alpha", "alpha.bam", &[7, 25], 4);
        let (mu_dir, mu) = sample_showing("mu", "mu.bam", &[]);
        let (iota_dir, iota) = sample_carrying_an_insertion("iota", "iota.bam", 2);
        let (kappa_dir, kappa) =
            sample_showing_two_sequences_once_each("kappa", "kappa.bam", b'C', b'G');
        (
            vec![zeta_dir, alpha_dir, mu_dir, iota_dir, kappa_dir],
            vec![zeta, alpha, mu, iota, kappa],
        )
    }

    /// Write the mixed cohort's VCF inside a pool of `threads`; the file's bytes and the
    /// run's answer come back for comparison.
    fn mixed_cohort_vcf_in_a_pool(
        threads: usize,
        paths: &[PathBuf],
        reference: &crate::ng::read::input::reference::OpenReference,
        merge: MergeParameters,
    ) -> (Vec<u8>, WrittenCohort) {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a fixture pool");
        pool.install(|| {
            let caller = open_over_the_tract_ground_with(paths, reference, merge);
            let sample_names: Vec<String> = caller.sample_names().map(str::to_owned).collect();
            let metadata = VcfHeaderMetadata::try_new(
                vec![HeaderContig {
                    name: "chr1".to_owned(),
                    length: 100,
                    md5: None,
                }],
                sample_names,
                "pop_var_caller_exp call-from-alignments".to_owned(),
                "reference.fa".to_owned(),
                String::new(),
            )
            .expect("a header this run can state");
            let out = tempfile::tempdir().expect("a temporary directory");
            let path = out.path().join("calls.vcf");
            let mut writer =
                VcfWriter::create(&path, metadata, Ploidy::try_new(2).expect("a diploid"))
                    .expect("the output opens");
            let written = caller
                .call_cohort_handing_each_record_over(&the_shipped_genotyper(), &mut |record| {
                    writer.write_record(record)
                })
                .expect("the mixed cohort calls");
            writer.finish().expect("the file is renamed into place");
            (std::fs::read(&path).expect("the VCF is there"), written)
        })
    }

    /// One sample's walk as the thread-count sweep compares it: its name, and the three
    /// tally structs whole rather than fields picked out of them.
    type SampleWalkTalliesForComparison = (
        String,
        LocusCounts,
        Vec<(Option<ReadGroupId>, ReadFilterCounts)>,
        Option<PileupGeneratorCounts>,
    );

    /// **Each sample's walk, named and counted** — the projection the thread-count sweep
    /// compares. The sample's own name rides along with its counts, so a run that produced
    /// the right multiset of tallies but paired them with the wrong samples fails; written
    /// once and called on both sides, so the two cannot drift into comparing different things.
    fn walk_tallies_of(run: &WrittenCohort) -> Vec<SampleWalkTalliesForComparison> {
        run.walk
            .per_sample
            .iter()
            .map(|walk| {
                (
                    walk.sample_name.clone(),
                    walk.regions.clone(),
                    walk.read_filters.clone(),
                    walk.snp_indel,
                )
            })
            .collect()
    }

    /// **Each sample's admitted-read count under its own name, against absolute numbers** —
    /// the guard the thread-count sweep structurally cannot provide, because that sweep
    /// compares the record path against itself.
    ///
    /// **⚑ This closes a hole the record path had and the oracle path did not.** Reversing
    /// the sample-name list before it is paired with the walkers in
    /// [`AlignedFilesVariantCaller::call_cohort_handing_each_record_over`] passed the whole
    /// library — measured, 437 of 437 green — while the same reversal in
    /// [`AlignedFilesVariantCaller::call_cohort`] is killed by two tests. A run report would
    /// then have carried one sample's read-drop rates under another sample's name: a wrong
    /// report rather than a crash, which is the failure `CohortWalkTallies::of`'s own length
    /// check exists to prevent and cannot see, both lists being the right length.
    ///
    /// The counts are what make the pairing checkable: the five samples walk **3, 4, 3, 3 and
    /// 6 reads**, so no two of them are interchangeable and any permutation shows.
    fn walk_counts_by_name(run: &WrittenCohort) -> Vec<(&str, u64)> {
        run.walk
            .per_sample
            .iter()
            .map(|walk| {
                (
                    walk.sample_name.as_str(),
                    walk.snp_indel
                        .as_ref()
                        .map_or(0, |counts| counts.reads_admitted),
                )
            })
            .collect()
    }

    /// **The written VCF is byte-identical at every thread count, at both widths, and across
    /// the widths** — spec §12.2's oracle, taken all the way to the artefact a person diffs,
    /// on the cohort and ground described above. A pool of one takes the serial-sweep
    /// fallback; pools of 2, 4, 8 and 16 take the Jacobi sweep, each run three times because
    /// a schedule-dependent divergence shows only under some interleavings. The seven-base
    /// width runs the same cohort over six building regions instead of one, so the file must
    /// not depend on where the merge cuts its ground either. The counts beside the file —
    /// written, called-but-not-written, refused, nobody-callable, and each sample's whole walk
    /// tallies — must agree too, because a run report built from them must not depend on the
    /// thread count either. **The two refusal lists are empty on this fixture**, so those two
    /// comparisons cannot fail today and are kept for the fixture that fills them; the walk
    /// tallies are compared as whole structs rather than as chosen fields, so a count added to
    /// `LocusCounts`, `ReadFilterCounts` or `PileupGeneratorCounts` later joins the
    /// comparison instead of silently dropping out of it.
    ///
    /// **⚑ What this reaches and what it does not, measured by three mutations to
    /// [`ObservationCache::cover_in_parallel`] rather than argued.**
    ///
    /// | mutation | effect on the sweep | these two tests |
    /// |---|---|---|
    /// | drop the last sample from `par_iter_mut` | a sample is never drawn forward | **both fail** |
    /// | `max` → `min` in the `try_reduce` | a cover stops at the least reach any sample grew to | both pass |
    /// | `break` after the first iteration | the fixpoint never iterates | both pass |
    ///
    /// So the oracle has real power over **who** the sweep draws and none at all over **how
    /// far** it keeps drawing. The reason is the fixture reference, and it is not fixable
    /// here: it is a hundred identical bases, so no observation ever reaches past a building
    /// region into where another sample's begins — two substitutions share no base to chain
    /// on, an insertion's reference span is its anchor alone, and a deletion left-aligns off
    /// the record (measured at D1). Two extra samples departing at adjacent positions were
    /// tried during E2 to give a cover something to chain, and closed as two separate loci.
    /// With no chain there is nothing for a second sweep to find, so a cover that stops early
    /// loses nothing *on this ground*.
    ///
    /// The fixpoint is pinned a layer down, on fixtures minted in memory that can hold a
    /// 26-base observation: the same two mutations are killed by
    /// `the_parallel_cover_gives_the_serial_drivers_answer` and five others in
    /// [`cohort_merge`](super::super::cohort_merge). **So the layering is deliberate — the
    /// cover's fixpoint at the merge, the end-to-end tie here**, and the one-position locus
    /// width this fixture is limited to is asserted by
    /// [`the_mixed_cohorts_records_describe_the_serial_callers_loci`] so it cannot silently
    /// stop being true.
    #[test]
    fn the_record_path_is_byte_identical_at_every_thread_count() {
        let (_dirs, paths) = the_mixed_cohort();
        let (_reference_dir, reference) = fixture_reference_from_its_index();

        let (default_width_bytes, baseline) =
            mixed_cohort_vcf_in_a_pool(1, &paths, &reference, MergeParameters::DEFAULT);
        assert_eq!(
            baseline.records_written, 3,
            "the shared SNP, alpha's own SNP and iota's insertion reach the file",
        );
        assert_eq!(
            baseline.loci_called_but_not_written, 1,
            "kappa's two-singletons locus is called and establishes nothing",
        );
        // **Every sample, not the first one.** A run hands one segmentation to the whole
        // cohort, so all five walk the same three segments and handle all three; reading only
        // `per_sample[0]` would pass just as well if the tallies were paired with the wrong
        // samples, or if four of the five were left empty.
        //
        // **The second number was 1 until the tract slot was filled** (C2, 2026-09-02): the
        // fixture's tract was a region whose generator did not exist, and is now a region with
        // one. What the run cannot do with the *locus* it builds there is a different count —
        // `tract_loci_set_aside` below — and keeping the two apart is the point: one says
        // ground nobody looked at, the other says loci nobody scored.
        assert_eq!(
            baseline
                .walk
                .per_sample
                .iter()
                .map(|walk| (
                    walk.regions.regions_in,
                    walk.regions.unhandled_not_implemented
                ))
                .collect::<Vec<_>>(),
            vec![(3, 0); paths.len()],
            "each sample walks the fixture's three segments and now handles all three, on \
             every thread count",
        );
        assert_eq!(
            baseline.tracts.built(),
            0,
            "no sample of this cohort varies inside the tract, so the merge finds it too quiet \
             to build and there is no locus to set aside — \
             `a_tract_a_sample_varies_at_is_built_and_set_aside_uncalled` is the fixture that \
             does vary there",
        );
        assert_eq!(
            walk_counts_by_name(&baseline),
            [
                ("zeta", 3),
                ("alpha", 4),
                ("mu", 3),
                ("iota", 3),
                ("kappa", 6),
            ],
            "each sample's admitted reads come back under that sample's own name",
        );

        for merge in the_two_widths() {
            let width = merge.cohort_locus_builder_regions_len.get();
            let (baseline_bytes, _) = mixed_cohort_vcf_in_a_pool(1, &paths, &reference, merge);
            // At the default width this compares a second serial run against the first, which
            // is worth its cost for a different reason than the name suggests: it is the only
            // place the whole path is run twice with everything held fixed, so a run that
            // depended on the order two temporary directories happened to be created in, or on
            // any other per-run state, fails here rather than being blamed on a thread count.
            assert_eq!(
                baseline_bytes, default_width_bytes,
                "a serial run at a {width}-base building region must give the file the shipped \
                 width gave",
            );
            for threads in [2, 4, 8, 16] {
                for repetition in 0..3 {
                    let (bytes, again) =
                        mixed_cohort_vcf_in_a_pool(threads, &paths, &reference, merge);
                    assert_eq!(
                        bytes, baseline_bytes,
                        "the VCF differs between a pool of 1 and a pool of {threads} at a \
                         {width}-base building region (repetition {repetition})",
                    );
                    assert_eq!(again.records_written, baseline.records_written);
                    assert_eq!(
                        again.loci_called_but_not_written,
                        baseline.loci_called_but_not_written
                    );
                    assert_eq!(
                        again.loci_too_wide_to_assemble,
                        baseline.loci_too_wide_to_assemble
                    );
                    assert_eq!(
                        again.loci_with_nobody_to_call,
                        baseline.loci_with_nobody_to_call
                    );
                    assert_eq!(
                        again.tracts, baseline.tracts,
                        "the tract loci set aside at a pool of {threads}",
                    );
                    assert_eq!(
                        walk_tallies_of(&again),
                        walk_tallies_of(&baseline),
                        "each sample's walk tallies at a pool of {threads}",
                    );
                }
            }
        }
    }

    /// One sample whose reads vary **inside** [`ground_with_a_tract_interleaved`]'s tract, so
    /// the merge has a tract locus to build.
    ///
    /// Forty-five bases from chr1:25, so 25–69 — the whole 41–52 tract with the fifteen bases
    /// of flank the generator fetches on each side, so its reads pin a repeat length rather
    /// than coming back as partials. The changed base is at 45, inside the tract. Four reads,
    /// where the merge's keep rule asks for two.
    ///
    /// **The mixed cohort does not vary there**, which is why its own set-aside count is zero
    /// and why this fixture exists.
    fn a_sample_varying_inside_the_tract() -> (TempDir, PathBuf) {
        let mut bases = [b'A'; 45];
        bases[20] = b'C';
        let records: Vec<RecordBuf> = (0..4)
            .map(|read| read_of(&format!("varies-r{read}"), 25, &bases))
            .collect();
        indexed_named_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("varies"))],
            ),
            &records,
            "varies.bam",
        )
    }

    /// **A cohort locus at a repeat tract is built, merged, and now called** — the dispatch
    /// this milestone builds, replacing the guard that stood in its place.
    ///
    /// **This test asserted the opposite until 2026-09-02 and the reversal is the step.** The
    /// guard set every tract aside because `call_one_generic_locus` knows nothing about
    /// repeats: it would have taken the distinct tract lengths the reads showed as ordinary
    /// alleles, ranked them by read support, and emitted a record whose `REF` and `ALT` are
    /// whole tract sequences scored under a substitution model with no stutter term. What
    /// replaces the guard is not that path — it is `select_ssr` and the tract model, reached
    /// through `call_one_cohort_locus`'s branch on the observation's kind.
    ///
    /// So the two halves invert together: **nothing is set aside** — the count now holds
    /// bundles alone — and **a locus over the tract's own ground reaches the caller**. A run
    /// that dispatched but scored the tract as ordinary sequence would pass the first half and
    /// is what the record-level fixtures below are for.
    ///
    /// The sample's reads span the tract with room either side, and carry one changed base
    /// inside it, so the merge has something to build.
    #[test]
    fn a_tract_a_sample_varies_at_is_built_and_called_through_the_tract_path() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();

        let (_bam_dir, bam) = a_sample_varying_inside_the_tract();

        let called = open_over_the_tract_ground_with(
            std::slice::from_ref(&bam),
            &reference,
            MergeParameters::DEFAULT,
        )
        .call_cohort(&the_shipped_genotyper())
        .expect("the fixture cohort calls");

        assert_eq!(
            called.tracts.bundles_set_aside, 0,
            "a repeat tract is dispatched to the tract path, and only a bundle is set aside",
        );
        let over_the_tract: Vec<&LocusInference> = called
            .called_loci
            .iter()
            .filter(|locus| locus.region.start.get() >= 41 && locus.region.end.get() <= 52)
            .collect();
        assert_eq!(
            over_the_tract.len(),
            1,
            "the tract's own ground is called once, and the caller produced {:?}",
            called
                .called_loci
                .iter()
                .map(|locus| locus.region)
                .collect::<Vec<_>>(),
        );
        // **Which model called it, and not merely that something did.** The SNP/indel path
        // would also produce a record over this ground — that is what the guard this dispatch
        // replaces existed to prevent — and the only thing in the answer that tells the two
        // apart is the kind selection stamped on the candidate table.
        assert!(
            matches!(over_the_tract[0].alleles().kind(), LocusKind::Ssr(_)),
            "the tract was called through the tract path, and its candidates say {:?}",
            over_the_tract[0].alleles().kind(),
        );
    }

    /// **Each of selection's verdicts lands in its own outcome**, and a truncated tract is not
    /// counted among the cleanly called ones.
    ///
    /// The counting cannot be read off the file: two of the three leave no record — a tract
    /// refused as `notPeriodic` is called over the reference alone, so it establishes no variant
    /// and is left out. So the mapping is asserted here rather than against an output.
    #[test]
    fn each_selection_verdict_lands_in_its_own_tract_outcome() {
        let outcome_of = |verdict| {
            let mut tracts = TractOutcomes::default();
            count_this_tract(verdict, &mut tracts);
            tracts
        };
        assert_eq!(outcome_of(SelectionVerdict::Selected).called, 1);
        assert_eq!(outcome_of(SelectionVerdict::NotPeriodic).not_periodic, 1);
        assert_eq!(
            outcome_of(SelectionVerdict::Truncated { dropped: 4 }).too_many_alleles,
            1,
        );
        assert_eq!(
            outcome_of(SelectionVerdict::Truncated { dropped: 4 }).called,
            0,
            "a truncated tract is called over what the cap kept, and is not counted among the \
             tracts nothing was cut from",
        );
        // Each of the three adds exactly one locus to the partition, and to one part of it.
        for verdict in [
            SelectionVerdict::Selected,
            SelectionVerdict::NotPeriodic,
            SelectionVerdict::Truncated { dropped: 1 },
        ] {
            assert_eq!(outcome_of(verdict).built(), 1);
        }
    }

    /// **A candidate carrying no whole motif copy stops the tract rather than being scored as
    /// one repeat** — the conversion that is the refusal.
    ///
    /// Selection floors a candidate's length by the period, so a sequence shorter than one
    /// copy of the unit comes back as zero whole repeats. The stutter ladder is written in
    /// whole repeats and has no rung below its bottom one; admitting such a candidate as
    /// `NonZeroU32::new(1)` would put it on the ladder's first rung, which is a different
    /// allele from the one the reads showed.
    #[test]
    fn a_candidate_with_no_whole_repeat_is_not_convertible_for_the_tract_model() {
        assert_eq!(
            repeat_counts_the_tract_model_can_take(&[7, 6, 5])
                .expect("every candidate carries a whole repeat")
                .iter()
                .map(|count| count.get())
                .collect::<Vec<u32>>(),
            vec![7, 6, 5],
        );
        assert!(
            repeat_counts_the_tract_model_can_take(&[7, 0]).is_none(),
            "one candidate below the ladder's bottom rung stops the locus",
        );
        assert!(
            repeat_counts_the_tract_model_can_take(&[]).is_some(),
            "an empty list is not a candidate with no repeats — it is no candidates, which \
             selection cannot produce, since the reference tract is always admitted",
        );
    }

    /// **And the records are the serial caller's loci** — the oracle the plan names for E2.
    /// `call_cohort` never touches the parallel cover, so agreement here ties the whole
    /// thread-count sweep above back to the one driver whose schedule has never changed.
    ///
    /// **⚑ The record path runs inside a pool of eight on purpose.** Left to the ambient
    /// pool it would take whatever `rayon::current_num_threads()` happened to be, and on a
    /// one-CPU runner — or under `RAYON_NUM_THREADS=1` — the driver takes its serial-sweep
    /// fallback, so this would quietly become the serial cover compared against itself and
    /// still pass. Naming the pool makes the comparison the one the test claims.
    #[test]
    fn the_mixed_cohorts_records_describe_the_serial_callers_loci() {
        let (_dirs, paths) = the_mixed_cohort();
        let (_reference_dir, reference) = fixture_reference_from_its_index();

        let called = open_over_the_tract_ground_with(&paths, &reference, MergeParameters::DEFAULT)
            .call_cohort(&the_shipped_genotyper())
            .expect("the serial oracle calls");
        let mut records = Vec::new();
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(8)
            .build()
            .expect("a fixture pool");
        let written = pool.install(|| {
            assert!(
                rayon::current_num_threads() > 1,
                "the record path must take the parallel sweep here, not its serial fallback",
            );
            open_over_the_tract_ground_with(&paths, &reference, MergeParameters::DEFAULT)
                .call_cohort_handing_each_record_over(&the_shipped_genotyper(), &mut |record| {
                    records.push(record.clone());
                    Ok::<(), std::io::Error>(())
                })
                .expect("the record path calls")
        });

        assert_eq!(
            written.loci_called(),
            called.called_loci.len() as u64,
            "written plus establishing-nothing is exactly what the serial caller called",
        );
        let written_regions: Vec<GenomeRegion> =
            records.iter().map(|record| record.region()).collect();
        let called_regions: Vec<GenomeRegion> = called
            .called_loci
            .iter()
            .map(|locus| locus.region)
            .collect();
        assert!(
            written_regions
                .iter()
                .all(|region| called_regions.contains(region)),
            "every record's span is a span the serial caller called: {written_regions:?} \
             against {called_regions:?}",
        );
        assert_eq!(
            called_regions.len() - written_regions.len(),
            written.loci_called_but_not_written as usize,
            "and the difference is exactly the called-but-not-written count",
        );

        // **The fixture's limitation, asserted rather than described.** Every locus this
        // cohort produces is one reference position wide, because a reference of a hundred
        // identical bases gives substitutions no shared base to chain on and slides deletions
        // off the record. That is why this module cannot pin the parallel cover's
        // chain-following and `cohort_merge`'s in-memory fixtures can. Asserting it keeps the
        // claim honest in both directions: if a later change to the mint or to candidate
        // selection ever does produce a wider locus here, this fails and the paragraph above
        // has to be rewritten rather than quietly becoming false.
        let widest = written_regions
            .iter()
            .map(|region| region.end.0.saturating_sub(region.start.0) + 1)
            .max()
            .expect("the fixture writes records");
        assert_eq!(
            widest, 1,
            "every locus this reference can express is one position wide: {written_regions:?}",
        );
    }

    /// **A cohort of one sample gives the same file at one thread and at eight.**
    ///
    /// The single low-coverage sample is the hardest end of the range this caller commits to
    /// (`CLAUDE.md`, spec §7.2), and it is the shape where a sweep over samples has the least
    /// to do: one sample means the Jacobi reduction folds a single value, so anything that
    /// only works because a second sample happened to widen the reach shows here and nowhere
    /// else in this module. `a_cohort_of_one_sample_writes_its_records` already proves such a
    /// run reaches the file; what it does not do is vary the thread count, because it runs
    /// under whatever pool the harness has.
    #[test]
    fn a_cohort_of_one_sample_is_byte_identical_at_every_thread_count() {
        let (_dirs, paths) = the_mixed_cohort();
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let alone = &paths[..1];

        for merge in the_two_widths() {
            let width = merge.cohort_locus_builder_regions_len.get();
            let (serial_bytes, serial) = mixed_cohort_vcf_in_a_pool(1, alone, &reference, merge);
            assert_eq!(
                serial.records_written, 1,
                "zeta alone still carries its SNP at chr1:17",
            );
            for threads in [2, 8] {
                for repetition in 0..3 {
                    let (bytes, again) =
                        mixed_cohort_vcf_in_a_pool(threads, alone, &reference, merge);
                    assert_eq!(
                        bytes, serial_bytes,
                        "one sample's VCF differs between a pool of 1 and a pool of \
                         {threads} at a {width}-base building region (repetition {repetition})",
                    );
                    assert_eq!(walk_tallies_of(&again), walk_tallies_of(&serial));
                }
            }
        }
    }

    /// How many records the written VCF holds over the fixture tract's own ground.
    ///
    /// **The tract is `chr1:41-52` in the fixture's segmentation**, which
    /// `a_tract_a_sample_varies_at_is_built_and_called_through_the_tract_path` names too; a
    /// record's `POS` is its first base, one-based, so the window is asked inclusively at both
    /// ends.
    fn records_over_the_tract(vcf: &[u8]) -> usize {
        String::from_utf8_lossy(vcf)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .filter_map(|line| line.split('\t').nth(1)?.parse::<u64>().ok())
            .filter(|position| (41..=52).contains(position))
            .count()
    }

    /// **The tract path is thread-invariant too, and on a cohort where it actually fires** —
    /// C4's half of the E2 oracle.
    ///
    /// The sweep above runs the tract generator at every thread count, but no sample of that
    /// cohort varies inside its tract, so the merge finds the tract too quiet to build and
    /// nothing the tract path does reaches the file at any pool size: **a comparison of
    /// files that hold no tract record, which any implementation passes.** This runs the same
    /// sweep on the sample that does vary there, so the bytes and the walk tallies are
    /// compared at a value the tract path produced.
    ///
    /// **The anchor moved with the dispatch (2026-09-02).** It used to be the set-aside count,
    /// which is now zero at a tract; what stands in its place is a record over the tract's own
    /// ground, which is the thing the sweep would otherwise be comparing the absence of.
    ///
    /// **What could go wrong here that the sweep above cannot see.** The tract generator holds
    /// a cursor and its own reference accessors, one set per sample, and a run's walkers cross
    /// threads under the merge's parallel cover. A generator that shared a window or a cursor
    /// between workers would answer differently depending on which worker reached it first —
    /// and what would move is the tract's own observation, which is exactly what the fixture
    /// above holds constant at nothing.
    #[test]
    fn the_tract_path_is_byte_identical_at_every_thread_count() {
        let (_reference_dir, reference) = fixture_reference_from_its_index();
        let (_bam_dir, bam) = a_sample_varying_inside_the_tract();
        let alone = [bam];

        for merge in the_two_widths() {
            let width = merge.cohort_locus_builder_regions_len.get();
            let (serial_bytes, serial) = mixed_cohort_vcf_in_a_pool(1, &alone, &reference, merge);
            assert!(
                records_over_the_tract(&serial_bytes) > 0,
                "the fixture varies inside the tract, so the tract path writes a record there \
                 — without this the comparisons below are of files holding none",
            );
            assert_eq!(
                serial.tracts.bundles_set_aside, 0,
                "a tract is called rather than set aside; the count is bundles alone",
            );
            for threads in [2, 4, 8, 16] {
                for repetition in 0..3 {
                    let (bytes, again) =
                        mixed_cohort_vcf_in_a_pool(threads, &alone, &reference, merge);
                    assert_eq!(
                        bytes, serial_bytes,
                        "the VCF differs between a pool of 1 and a pool of {threads} at a \
                         {width}-base building region (repetition {repetition})",
                    );
                    assert_eq!(
                        again.tracts, serial.tracts,
                        "the tract loci set aside at a pool of {threads}",
                    );
                    assert_eq!(
                        records_over_the_tract(&bytes),
                        records_over_the_tract(&serial_bytes),
                        "the records over the tract's own ground at a pool of {threads}",
                    );
                    assert_eq!(
                        walk_tallies_of(&again),
                        walk_tallies_of(&serial),
                        "the walk tallies at a pool of {threads} — the tract generator's own \
                         read-filter counts among them",
                    );
                }
            }
        }
    }
}
