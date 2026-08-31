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
use crate::ng::read::input::read_groups::ReadGroups;
use crate::ng::read::input::reference::OpenReference;
use crate::ng::run::cohort_merge::{CohortLocusBuilderRegionsLen, MaxCohortLocusSpan, MinAltReads};

use super::RunError;
use super::segments::Segmentation;

/// The alignment files a run reads, and how it reads them.
///
/// **Grouped rather than passed one by one**, because they are answers to one question — what
/// this run's reads are — and because the four travel together into every sample's open.
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

        Ok(Self {
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
    fn bam_for(sample: &str, file_name: &str) -> (TempDir, PathBuf) {
        let (dir, path) = unindexed_bam_for(sample, file_name);
        index(&path);
        (dir, path)
    }

    /// The same file with no index beside it — what a run refuses when it was not asked to
    /// build one.
    fn unindexed_bam_for(sample: &str, file_name: &str) -> (TempDir, PathBuf) {
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
    fn bam_holding_two_samples(first: &str, second: &str, file_name: &str) -> (TempDir, PathBuf) {
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

    fn catalog_header() -> RepeatCatalogHeader {
        RepeatCatalogHeader {
            contigs: Vec::new(),
            reference_md5: [7; 16],
            built_under: StrRepeatCriteria::default(),
            scan: ScanParams::default(),
            tool_version: "test".to_string(),
            longest_tract_bp: Vec::new(),
        }
    }

    /// A segmentation over one short contig, with one ordinary stretch in it.
    fn segmentation() -> Segmentation {
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
            catalog_header(),
            StrRepeatCriteria::default(),
            PathBuf::from("/genomes/test.catalog.parquet"),
        )
        .expect("a clean stream builds")
    }

    // -----------------------------------------------------------------
    // Settings that are NOT their type's default
    //
    // A test that hands in a default and asserts the default back cannot tell "held what it
    // was given" from "replaced with the default" — the shape that let four mutations survive
    // an earlier draft of this suite. Everything below differs from what ships.
    // -----------------------------------------------------------------

    fn unusual_read_filters() -> ReadFilterConfig {
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
