//! psp mode's walk stage: one sample's observations, gathered into a psp.
//!
//! **The gatherer is the walk half of psp mode** (`doc/devel/ng/spec/run_streaming.md` §5.2;
//! `doc/devel/ng/arch/run_streaming.md` §3.3): one sample's alignment files in, that sample's
//! observations out, in genome order, ready to be written to a psp. It is the same machinery
//! direct mode walks — [`SampleLocusObservationsIterator`] over [`RunSegments`], the chain
//! [`AlignedFilesVariantCaller`](super::AlignedFilesVariantCaller) builds one of per sample —
//! wrapped the other way: direct mode's walker answers a merge one observation at a time,
//! where the gatherer is an [`Iterator`] a psp writer drains to the end.
//!
//! **One gatherer, one sample, serial within the sample** (spec §5.2). The walk stage is a
//! loop of these, one per sample, and a cohort is parallelised by running invocations — so
//! nothing here holds more than one sample's files, and the observations leave in the order
//! one deterministic walk produces them, which is what makes the written file reproducible
//! byte for byte (spec §12.1).
//!
//! **The gatherer builds the header, and it is the first production code that does.** Every
//! field a calling run will later check — the analysed regions and segmentation inputs, the
//! walk-local read-group table, the observation reach ceiling, the read filters as provenance
//! (spec §6.1) — is known at construction, before a record exists, so the header is built at
//! [`open`](SampleObservationGatherer::open) and a file whose header cannot be built is
//! refused before any walking starts.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
use crate::ng::locus_generation::{
    LocusCounts, SampleLocusObservations, SampleLocusObservationsIterator,
};
use crate::ng::psp::{
    ContigIdentity, FORMAT_VERSION, Header, Manifest, PspWriter, ReadGroupIdentity,
    ReferenceIdentity, WriteStats, WriterProvenance,
};
use crate::ng::read::ReadFilterConfig;
use crate::ng::read::input::read_groups::{ReadGroups, build_read_groups};
use crate::ng::read::input::reference::OpenReference;
use crate::ng::read::input::{IngestError, SampleReads};
use crate::ng::reference_info::{ContigInfo, ReferenceInfo};
use crate::ng::types::Bp;

use super::walker::{WalkReference, generic_path_generators};
use super::{RunError, RunSegments, Segmentation, WalkProgress};

/// What one sample's walk opens and how it filters — the per-sample counterpart of direct
/// mode's [`AlignmentInputs`](super::AlignmentInputs), minus everything that only exists at
/// cohort scope (the run-wide read-group table, the checksummed reference the assembly check
/// compares against — both stay with the command that loops the samples).
pub struct SampleWalkInputs<'a> {
    /// One sample's alignment files. **Every path must name the same sample**: the gatherer
    /// builds its read-group table over exactly these, so files naming two samples are
    /// refused ([`RunError::FilesNotFromOneSample`]).
    pub alignments: &'a [PathBuf],
    /// The reference the walk reads against. It must hold bases — a `.fai`-only reference
    /// is refused ([`RunError::ReferenceHasNoBases`]), exactly as for a calling run.
    pub reference: &'a OpenReference,
    /// Which reads the walk admits — recorded in the header's provenance, never compared
    /// (spec §6.1).
    pub read_filters: ReadFilterConfig,
    /// The five knobs the locus generator walks with. `max_record_span` is also the
    /// header's `observation_reach_ceiling_bp`: **the configured value, not a measured
    /// one** — it is a setting, known before the first record (spec §6.1;
    /// `psp_file_format.md` §3.1).
    pub locus_generator_settings: PileupGeneratorConfig,
    /// Whether a missing alignment index is built or refused.
    pub build_index_if_missing: bool,
}

/// One sample's observations in genome order, ready to become a psp (arch §3.3).
///
/// An [`Iterator`] over the sample's [`SampleLocusObservations`], plus the [`Header`] the
/// file they belong in will carry — built at [`open`](Self::open), before any record exists.
/// [`write_psp`](Self::write_psp) is the ordinary way to consume one: it drains the iterator
/// into a [`PspWriter`] record by record.
///
/// **Not `Clone`**, for the walker's reason: it owns one sample's open files and its
/// generators' accumulated state, and a second gatherer over one sample would decode the
/// same ground twice.
///
/// **The census accumulator is not here yet.** Arch §3.3 puts it at the yield point and
/// `finish()` hands the evidence over; that wiring is the plan's Milestone G, and the
/// constructor will grow the census's inputs when it lands
/// (`doc/devel/ng/impl_plan/run_driver_psp_mode.md`, G1).
pub struct SampleObservationGatherer {
    /// The header the sample's psp will carry — everything spec §6.1 asks for, fixed
    /// before the first record. The sample's name lives here (`header.sample`), not in a
    /// second field: one struct, one copy of the fact.
    header: Header,
    /// How far the walk has got — the second half of locating a failure (spec §9).
    reached: WalkProgress,
    loci: SampleLocusObservationsIterator<RunSegments>,
}

impl SampleObservationGatherer {
    /// Open one sample's files and fix the header its psp will carry.
    ///
    /// **Named `open` rather than `new` for the run's reason**: the sample's files are
    /// opened, validated and index-checked here, before a single read flows.
    ///
    /// The read-group table is built over exactly this sample's files, so its numbering
    /// starts at zero — the walk-local numbering the header records and the calling stage
    /// later merges run-wide (spec §6.1, §6.2).
    ///
    /// `provenance` is what only the caller can know — the tool, its version, the
    /// subcommand, the command line, when the run started. The gatherer overwrites what it
    /// knows better: the input basenames from the files it actually opens, and the read
    /// filters recorded as provenance parameters
    /// ([`ReadFilterConfig::provenance_parameters`]).
    ///
    /// # Defaults
    ///
    /// The header's encoding manifest — block geometry and field list — is fixed to
    /// [`Manifest::as_this_build_writes_it`]; there is no knob here to change it. The
    /// chosen values are readable from [`header`](Self::header) before anything is
    /// written, and are recorded in the file itself, which is what drives every reader.
    ///
    /// # Errors
    ///
    /// An empty file list is [`RunError::NoAlignmentFiles`]; impossible generator settings
    /// are [`RunError::LocusGeneratorSettings`]; a reference holding no bases is
    /// [`RunError::ReferenceHasNoBases`]; files that cannot be read as one sample's —
    /// unreadable read-group headers, or two samples named — are
    /// [`RunError::FilesNotFromOneSample`]; and a sample whose files will not open is
    /// [`RunError::OpeningSample`], naming it.
    pub fn open(
        inputs: SampleWalkInputs<'_>,
        segmentation: Arc<Segmentation>,
        provenance: WriterProvenance,
    ) -> Result<Self, RunError> {
        if inputs.alignments.is_empty() {
            return Err(RunError::NoAlignmentFiles);
        }
        // The same door-checks a calling run makes, for the same reason: a walk whose
        // settings are impossible or whose reference holds no bases is wrong at the door,
        // and a refusal at the first locus would arrive after the files were opened.
        inputs
            .locus_generator_settings
            .check()
            .map_err(|source| RunError::LocusGeneratorSettings { source })?;
        let walk_reference = WalkReference::of(inputs.reference)?;

        // One sample's table, built over exactly this sample's files. The classification —
        // one sample, none, or a mismatch listed per read group — is the table's own
        // (`ReadGroups::only_sample`), shared with `SampleReads::open_only_sample`; what
        // this call site keeps that one throws away is the table, which the header needs.
        let read_groups = build_read_groups(inputs.alignments).map_err(|source| {
            RunError::FilesNotFromOneSample {
                source: Box::new(IngestError::ReadGroups(source)),
            }
        })?;
        let only_sample = read_groups.only_sample().map_err(|source| match source {
            IngestError::NoFiles => RunError::NoAlignmentFiles,
            source => RunError::FilesNotFromOneSample {
                source: Box::new(source),
            },
        })?;
        let reads = SampleReads::open(
            only_sample,
            &read_groups,
            inputs.reference,
            inputs.read_filters,
            inputs.build_index_if_missing,
        )
        .map_err(|source| RunError::OpeningSample {
            sample: only_sample.sample.to_string(),
            source: Box::new(source),
        })?;

        let header = header_for(
            &inputs,
            reads.sample_name().to_string(),
            &read_groups,
            &walk_reference,
            &segmentation,
            provenance,
        );
        let generators = generic_path_generators(
            &walk_reference,
            inputs.locus_generator_settings,
            // The generators derive the tract generator's bundle radius from here — the
            // criteria the ground was actually cut with, the same rule as direct mode's
            // by construction.
            segmentation.inputs(),
        )?;
        Ok(Self {
            header,
            reached: WalkProgress::NothingYet,
            loci: SampleLocusObservationsIterator::new(
                RunSegments::of(segmentation),
                reads,
                generators,
            ),
        })
    }

    /// The header the sample's psp will carry — fixed at [`open`](Self::open).
    #[must_use]
    pub fn header(&self) -> &Header {
        &self.header
    }

    /// The individual this walk is of.
    #[must_use]
    pub fn sample_name(&self) -> &str {
        &self.header.sample
    }

    /// How far this walk has got: `NothingYet` until the first observation is yielded, and
    /// the last base of the last one after that.
    #[must_use]
    pub fn reached(&self) -> WalkProgress {
        self.reached
    }

    /// The walk's running tally — current at any point, final once the walk is spent.
    /// [`write_psp`](Self::write_psp) hands the final tally back, since it consumes the
    /// gatherer.
    #[must_use]
    pub fn counts(&self) -> &LocusCounts {
        self.loci.counts()
    }

    /// Drain this walk into a psp at `path`: header, every observation as one record, seal —
    /// the trailer is sealed **empty**, a legal trailer whose contents stay opaque until
    /// something needs them (`psp_file_format.md` §3.4). Returns the store's write totals
    /// and the walk's final tally, side by side — the two halves of a per-sample report.
    ///
    /// This is spec §5.2's "psp writer consumes the gatherer", whole. The walk is serial and
    /// the records leave in genome order, so the file is the walk — which is what lets
    /// gathering the same sample twice produce the same bytes apart from the header's
    /// timestamp (spec §12.1).
    ///
    /// **⚠ An existing file at `path` is truncated** — [`PspWriter::create`]'s contract. A
    /// caller that must not destroy a finished psp checks for one first
    /// (`generate-psps` does, plan step C3).
    ///
    /// # Errors
    ///
    /// A walk failure names the sample and how far it got ([`RunError::SourceFailed`]); a
    /// record the store refuses names its locus ([`RunError::RecordNotWritten`]); a file
    /// that cannot be created or sealed names its path ([`RunError::PspNotWritten`]). After
    /// any of these the file at `path` is not whole, and a reader will refuse it as
    /// interrupted rather than read it as complete (`psp_file_format.md` §10).
    pub fn write_psp(mut self, path: &Path) -> Result<(WriteStats, LocusCounts), RunError> {
        let mut writer = PspWriter::create(path, self.header.clone()).map_err(|source| {
            RunError::PspNotWritten {
                path: path.to_path_buf(),
                source: Box::new(source),
            }
        })?;
        for observation in &mut self {
            let observation = observation?;
            writer
                .push(&observation)
                .map_err(|source| RunError::RecordNotWritten {
                    locus: observation.region,
                    source: Box::new(source),
                })?;
        }
        let counts = self.loci.counts().clone();
        let stats = writer
            .finish(&[])
            .map_err(|source| RunError::PspNotWritten {
                path: path.to_path_buf(),
                source: Box::new(source),
            })?;
        Ok((stats, counts))
    }
}

/// The header spec §6.1 asks for, assembled from what `open` established. Split out so
/// `open`'s own body stays the control flow — which checks refuse, in what order.
fn header_for(
    inputs: &SampleWalkInputs<'_>,
    sample: String,
    read_groups: &ReadGroups,
    walk_reference: &WalkReference,
    segmentation: &Segmentation,
    mut provenance: WriterProvenance,
) -> Header {
    // The FASTA's basename is the printable identity the header wants — never the
    // directory, which would leak the producer's layout (spec §6.1). The path comes from
    // the walk reference, whose construction already refused a reference without one.
    let reference_basename = walk_reference
        .fasta_path()
        .file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .expect("WalkReference::of accepted this reference, so its FASTA path names a file");
    provenance.input_alignments = inputs
        .alignments
        .iter()
        .map(|path| {
            path.file_name()
                .expect("a path that opened as an alignment file has a final component")
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    provenance.input_reference = reference_basename.clone();
    provenance.record_parameters(inputs.read_filters.provenance_parameters());

    // Exhaustive destructures of the *source* types, so a field added to either must be
    // dispositioned here — recorded, or discarded by name — rather than silently skipped
    // by the only production code that builds psp headers.
    let ReferenceInfo {
        md5,
        contigs,
        fasta_path: _,
    } = inputs.reference.info();
    Header {
        format_version: FORMAT_VERSION,
        sample,
        reference: ReferenceIdentity {
            name: reference_basename,
            md5: *md5,
        },
        contigs: contigs
            .iter()
            .map(|contig| {
                // The `.fai` geometry stays out of the header deliberately
                // (`ReferenceIdentity`'s doc): it does not survive a round trip and would
                // leak the producer's file layout.
                let ContigInfo {
                    name,
                    length,
                    md5,
                    offset: _,
                    line_bases: _,
                    line_width: _,
                } = contig;
                ContigIdentity {
                    name: name.clone(),
                    length: *length,
                    md5: *md5,
                }
            })
            .collect(),
        read_groups: read_groups
            .iter()
            .map(|(walk_local_id, group)| ReadGroupIdentity {
                id: group.id.to_string(),
                library: group.library.value.to_string(),
                walk_local_id,
            })
            .collect(),
        observation_reach_ceiling_bp: Bp(u64::from(
            inputs.locus_generator_settings.max_record_span,
        )),
        writer: provenance,
        segmentation_inputs: segmentation.inputs().clone(),
        manifest: Manifest::as_this_build_writes_it(),
    }
}

/// The sample's name and how far it has got — not its reads. A derived `Debug` would print
/// every open file and every generator's accumulated state.
impl std::fmt::Debug for SampleObservationGatherer {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // No `..`: a new field must be added here, printed or withheld by name.
        let Self {
            header,
            reached,
            loci: _,
        } = self;
        formatter
            .debug_struct("SampleObservationGatherer")
            .field("sample", &header.sample)
            .field("reached", reached)
            .finish_non_exhaustive()
    }
}

impl Iterator for SampleObservationGatherer {
    type Item = Result<SampleLocusObservations, RunError>;

    /// The next observation this sample has, or `None` once its ground is walked.
    ///
    /// The same failure shape as direct mode's walker, because the two are one machinery:
    /// an error names the sample and where the walk had reached (spec §9), and exhaustion
    /// is final — the wrapped iterator latches on both.
    fn next(&mut self) -> Option<Self::Item> {
        match self.loci.next()? {
            Ok(observation) => {
                // `reach_position` rather than `region.end`, for the walker's reason: the
                // crate keeps the reach rule in one place.
                self.reached = WalkProgress::After(observation.reach_position());
                Some(Ok(observation))
            }
            Err(source) => Some(Err(RunError::SourceFailed {
                sample: self.header.sample.clone(),
                reached: self.reached,
                source: Box::new(source),
            })),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::psp::{ParameterValue, PspReader, ZSTD_COMPRESSION_LEVEL};
    use crate::ng::read::input::test_fixtures::{
        FixtureReadGroup, fixture_reference, header_with_read_groups, matching_contigs, named_bam,
        read_named_with_length_in_read_group,
    };
    use crate::ng::region_typing::{GenomeRegions, RegionKind, TypedRegion};
    use crate::ng::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};
    use crate::ng::tandem_repeat::ScanParams;
    use crate::ng::types::{ContigId, GenomeRegion, MapQual, Position, ReadGroupId};
    use crate::regions::ContigBounds;
    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record::MappingQuality;
    use std::collections::BTreeMap;
    use tempfile::TempDir;

    /// The two fixture contigs, as the fixture reference declares them.
    fn fixture_bounds() -> [ContigBounds<'static>; 2] {
        [
            ContigBounds {
                name: "chr1",
                length: 100,
            },
            ContigBounds {
                name: "chr2",
                length: 200,
            },
        ]
    }

    /// A segmentation over `segments`, analysing `analysed`, under the catalog header and
    /// criteria every fixture here shares — **the one copy of that scaffold**, so a change
    /// to `Segmentation::build` or to the catalog header's fields touches one site.
    fn build_segmentation(
        segments: Vec<TypedRegion>,
        analysed: GenomeRegions,
    ) -> Arc<Segmentation> {
        Arc::new(
            Segmentation::build(
                segments.into_iter().map(Ok),
                analysed,
                RepeatCatalogHeader {
                    contigs: Vec::new(),
                    reference_md5: [7; 16],
                    built_under: StrRepeatCriteria::default(),
                    scan: ScanParams::default(),
                    tool_version: "test".to_string(),
                    longest_tract_bp: Vec::new(),
                },
                StrRepeatCriteria::default(),
                PathBuf::from("/genomes/test.catalog.parquet"),
            )
            .expect("a clean stream builds"),
        )
    }

    /// A whole-contig generic segment.
    fn generic_segment(contig: u32, length: u64) -> TypedRegion {
        TypedRegion {
            region: GenomeRegion {
                contig: ContigId(contig),
                start: Position(1),
                end: Position(length),
            },
            kind: RegionKind::Generic,
        }
    }

    /// One generic segment per fixture contig — the whole of the fixture reference's ground,
    /// so every read the BAM fixtures place is on the walk.
    fn segmentation() -> Arc<Segmentation> {
        build_segmentation(
            vec![generic_segment(0, 100), generic_segment(1, 200)],
            GenomeRegions::whole_contigs(&fixture_bounds()),
        )
    }

    /// What only the caller knows, with one parameter of its own — so a test can see the
    /// gatherer *add* to the map rather than replace it — and a fixed timestamp, because a
    /// gatherer stamps nothing itself.
    fn provenance() -> WriterProvenance {
        WriterProvenance {
            tool: "ng".to_string(),
            version: "0.0.0-test".to_string(),
            subcommand: "generate-psps".to_string(),
            input_alignments: vec!["to-be-overwritten".to_string()],
            input_reference: "to-be-overwritten".to_string(),
            command_line: "ng generate-psps --test".to_string(),
            parameters: BTreeMap::from([("depth-cap".to_string(), ParameterValue::Integer(300))]),
            created: "2026-09-03T00:00:00Z".parse().expect("a datetime"),
        }
    }

    // -----------------------------------------------------------------
    // Settings that are NOT their type's default — a gatherer that dropped what it was
    // given, walking with the shipped constants while recording the configured ones (or
    // the reverse), must be visible. The fixture is built to discriminate BOTH sides:
    // `the_fixture_tells_applied_settings_from_the_defaults` proves it can.
    // -----------------------------------------------------------------

    fn unusual_read_filters() -> ReadFilterConfig {
        ReadFilterConfig {
            min_mapq: Some(MapQual(37)),
            ..ReadFilterConfig::default()
        }
    }

    fn unusual_locus_generator_settings() -> PileupGeneratorConfig {
        PileupGeneratorConfig {
            max_record_span: 4_321,
            // Two, against the fixture's three stacked reads — what makes an applied
            // default-substitution visible in the walk's output, not only in the header.
            max_snp_column_depth: 2,
            ..PileupGeneratorConfig::default()
        }
    }

    /// A read the two filter policies disagree about: MAPQ 30 sits between the default
    /// floor (20) and `unusual_read_filters`' 37.
    fn read_at_mapq_30(qname: &str, start: usize, read_group_id: &str) -> RecordBuf {
        let mut record = read_named_with_length_in_read_group(qname, 0, start, 30, read_group_id);
        *record.mapping_quality_mut() = MappingQuality::new(30);
        record
    }

    /// An indexed BAM whose one read group carries `sample` and `library`, with reads the
    /// walk yields — including three stacked at one position (the depth cap's prey) and one
    /// at MAPQ 30 (the filter floor's) so applied settings are visible in the output.
    fn bam_with_library(sample: &str, library: &str, file_name: &str) -> (TempDir, PathBuf) {
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records = vec![
            read_named_with_length_in_read_group(&format!("{stem}-r0"), 0, 5, 30, "rg1"),
            read_named_with_length_in_read_group(&format!("{stem}-r1"), 0, 12, 30, "rg1"),
            read_named_with_length_in_read_group(&format!("{stem}-s0"), 0, 50, 30, "rg1"),
            read_named_with_length_in_read_group(&format!("{stem}-s1"), 0, 50, 30, "rg1"),
            read_named_with_length_in_read_group(&format!("{stem}-s2"), 0, 50, 30, "rg1"),
            read_at_mapq_30(&format!("{stem}-low"), 55, "rg1"),
        ];
        bam_in_read_group(sample, library, "rg1", file_name, records)
    }

    /// The same shape under a different `@RG ID` — so a two-file sample's table has two rows.
    fn second_bam_with_library(sample: &str, library: &str, file_name: &str) -> (TempDir, PathBuf) {
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records = vec![read_named_with_length_in_read_group(
            &format!("{stem}-r0"),
            0,
            40,
            30,
            "rg2",
        )];
        bam_in_read_group(sample, library, "rg2", file_name, records)
    }

    /// One read on each contig, so a failure can land after progress.
    fn bam_with_reads_on_both_contigs(sample: &str, file_name: &str) -> (TempDir, PathBuf) {
        let stem = file_name.split('.').next().unwrap_or(file_name);
        let records = vec![
            read_named_with_length_in_read_group(&format!("{stem}-c1"), 0, 5, 30, "rg1"),
            read_named_with_length_in_read_group(&format!("{stem}-c2"), 1, 10, 30, "rg1"),
        ];
        bam_in_read_group(sample, "libA", "rg1", file_name, records)
    }

    fn bam_in_read_group(
        sample: &str,
        library: &str,
        read_group_id: &str,
        file_name: &str,
        records: Vec<RecordBuf>,
    ) -> (TempDir, PathBuf) {
        let sam_header = header_with_read_groups(
            Some("coordinate"),
            &matching_contigs(),
            &[FixtureReadGroup {
                id: read_group_id,
                sample: Some(sample),
                library: Some(library),
                platform: None,
            }],
        );
        let (dir, path) = named_bam(&sam_header, &records, file_name);
        index(&path);
        (dir, path)
    }

    fn index(path: &PathBuf) {
        crate::bam::index_preflight::preflight_alignment_indexes(std::slice::from_ref(path), true)
            .expect("build index");
    }

    /// A gatherer over `segmentation` at the unusual settings — **the one place those
    /// settings are spelled**, so the "every setting differs from its default" discipline
    /// cannot drift between the tests that rely on it.
    fn open_gatherer_over(
        alignments: &[PathBuf],
        reference: &OpenReference,
        segmentation: Arc<Segmentation>,
    ) -> Result<SampleObservationGatherer, RunError> {
        SampleObservationGatherer::open(
            SampleWalkInputs {
                alignments,
                reference,
                read_filters: unusual_read_filters(),
                locus_generator_settings: unusual_locus_generator_settings(),
                build_index_if_missing: false,
            },
            segmentation,
            provenance(),
        )
    }

    /// The same, over the module's shared two-contig ground.
    fn open_gatherer(
        alignments: &[PathBuf],
        reference: &OpenReference,
    ) -> Result<SampleObservationGatherer, RunError> {
        open_gatherer_over(alignments, reference, segmentation())
    }

    /// A fixture walk through the bare chain direct mode builds — the oracle a gatherer's
    /// yield is compared against, over whatever ground and settings the caller names.
    fn direct_walk(
        paths: &[PathBuf],
        reference: &OpenReference,
        segmentation: Arc<Segmentation>,
        read_filters: ReadFilterConfig,
        settings: PileupGeneratorConfig,
    ) -> Vec<SampleLocusObservations> {
        let reads = SampleReads::open_only_sample(paths, reference, read_filters, false)
            .expect("the direct open");
        let generators = generic_path_generators(
            &WalkReference::of(reference).expect("bases"),
            settings,
            segmentation.inputs(),
        )
        .expect("generators");
        SampleLocusObservationsIterator::new(RunSegments::of(segmentation), reads, generators)
            .collect::<Result<_, _>>()
            .expect("a clean walk")
    }

    /// **The header holds what spec §6.1 asks for, from the values actually given** — every
    /// asserted setting differs from its default, so a gatherer that recorded the shipped
    /// constants instead fails here.
    #[test]
    fn open_fixes_the_header_a_calling_run_will_check() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let (_b_dir, bam_b) = second_bam_with_library("NA12878", "libB", "b.bam");

        let gatherer = open_gatherer(&[bam_a, bam_b], &reference).expect("one sample's files open");
        let header = gatherer.header();

        assert_eq!(gatherer.sample_name(), "NA12878");
        assert_eq!(header.sample, "NA12878");

        // The walk-local read-group table: numbered from zero, in file order, each row's
        // library the one its file declared.
        assert_eq!(
            header.read_groups,
            vec![
                ReadGroupIdentity {
                    id: "rg1".to_string(),
                    library: "libA".to_string(),
                    walk_local_id: ReadGroupId(0),
                },
                ReadGroupIdentity {
                    id: "rg2".to_string(),
                    library: "libB".to_string(),
                    walk_local_id: ReadGroupId(1),
                },
            ],
        );

        // The reach ceiling is the CONFIGURED cap — a setting, not a measurement, and not
        // the default.
        assert_eq!(header.observation_reach_ceiling_bp, Bp(4_321));
        assert_ne!(
            unusual_locus_generator_settings().max_record_span,
            PileupGeneratorConfig::default().max_record_span,
            "the fixture must differ from the default, or this test proves nothing",
        );

        // The segmentation's own record, whole.
        assert_eq!(&header.segmentation_inputs, segmentation().inputs());

        // The reference as the run opened it: basename + whole-assembly digest, one contig
        // row per `.fai` entry with its own digest.
        let info = reference.info();
        let reference_basename = info
            .fasta_path
            .as_deref()
            .and_then(Path::file_name)
            .and_then(|name| name.to_str());
        assert!(info.md5.is_some(), "the digested fixture carries an md5");
        assert_eq!(header.reference.md5, info.md5);
        assert_eq!(Some(header.reference.name.as_str()), reference_basename);
        assert_eq!(header.contigs.len(), info.contigs.len());
        for (identity, contig) in header.contigs.iter().zip(&info.contigs) {
            assert_eq!(identity.name, contig.name);
            assert_eq!(identity.length, contig.length);
            assert_eq!(identity.md5, contig.md5);
        }

        // Provenance: the caller's fields survive, the gatherer's own are overwritten from
        // what it opened, and the read filters land beside the caller's parameters.
        assert_eq!(header.writer.subcommand, "generate-psps");
        assert_eq!(header.writer.input_alignments, vec!["a.bam", "b.bam"]);
        assert_eq!(
            Some(header.writer.input_reference.as_str()),
            reference_basename,
        );
        assert_eq!(
            header.writer.parameters.get("depth-cap"),
            Some(&ParameterValue::Integer(300)),
            "the caller's own parameter survives the read-filter recording",
        );
        assert_eq!(
            header.writer.parameters.get("read-filter-min-mapq"),
            Some(&ParameterValue::Integer(37)),
            "the applied filter is recorded, not the default",
        );

        assert_eq!(header.format_version, crate::ng::psp::FORMAT_VERSION);
        assert_eq!(header.manifest, Manifest::as_this_build_writes_it());
    }

    /// **The fixture can tell applied settings from the defaults** — the guard that makes
    /// the parity test below able to see a gatherer walking with shipped constants while
    /// recording the configured ones. Each dimension is varied alone: if either
    /// `assert_ne!` fails, the fixture went blind on that axis and the parity test proves
    /// nothing about it.
    #[test]
    fn the_fixture_tells_applied_settings_from_the_defaults() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let paths = vec![bam_a];

        let unusual = direct_walk(
            &paths,
            &reference,
            segmentation(),
            unusual_read_filters(),
            unusual_locus_generator_settings(),
        );
        let default_filters = direct_walk(
            &paths,
            &reference,
            segmentation(),
            ReadFilterConfig::default(),
            unusual_locus_generator_settings(),
        );
        let default_generator = direct_walk(
            &paths,
            &reference,
            segmentation(),
            unusual_read_filters(),
            PileupGeneratorConfig::default(),
        );

        assert_ne!(
            unusual, default_filters,
            "the MAPQ-30 read must separate the two filter policies",
        );
        assert_ne!(
            unusual, default_generator,
            "the three stacked reads must separate the two depth caps",
        );
    }

    /// **The gatherer is the direct walk** — the same fixture through the bare iterator chain
    /// direct mode builds gives the same observations, in the same order. With the guard
    /// above, this is also the applied-settings pin: a gatherer that walked with defaults
    /// would differ from this oracle.
    #[test]
    fn the_gatherer_yields_what_the_direct_walk_yields() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let paths = vec![bam_a];

        let mut gatherer = open_gatherer(&paths, &reference).expect("opens");
        let gathered: Vec<SampleLocusObservations> = gatherer
            .by_ref()
            .collect::<Result<_, _>>()
            .expect("a clean walk");

        // The spent walk's tally is final and consistent with what was yielded.
        assert_eq!(gatherer.counts().loci_emitted, gathered.len() as u64);
        assert_eq!(
            gatherer.counts().regions_in,
            2,
            "the fixture segmentation has two segments",
        );

        let walked = direct_walk(
            &paths,
            &reference,
            segmentation(),
            unusual_read_filters(),
            unusual_locus_generator_settings(),
        );
        assert!(
            !walked.is_empty(),
            "the fixture must yield observations, or this test compares two empty walks",
        );
        assert_eq!(gathered, walked);
    }

    /// `reached()` moves to `After` the last observation's reach — the success half of the
    /// progress contract a failure report depends on (spec §9).
    #[test]
    fn reached_advances_to_the_last_observations_reach() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let mut gatherer = open_gatherer(&[bam_a], &reference).expect("opens");

        assert_eq!(gatherer.reached(), WalkProgress::NothingYet);
        let first = gatherer
            .next()
            .expect("the fixture yields")
            .expect("a clean first draw");
        assert_eq!(
            gatherer.reached(),
            WalkProgress::After(first.reach_position()),
        );
    }

    /// A mid-walk failure names the sample and how far the walk had got — not "no
    /// observations yet" (spec §9). The FASTA is opened lazily per contig, so deleting it
    /// after chr1's draws makes chr2's fetch the failure.
    #[test]
    fn a_midwalk_failure_names_the_sample_and_its_progress() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_bam_dir, bam) = bam_with_reads_on_both_contigs("NA12878", "both.bam");
        let mut gatherer = open_gatherer(&[bam], &reference).expect("opens");

        let first = gatherer
            .next()
            .expect("chr1 yields")
            .expect("a clean first draw");
        let fasta = reference
            .info()
            .fasta_path
            .clone()
            .expect("the digested fixture has bases");
        std::fs::remove_file(&fasta).expect("remove the fixture FASTA");

        // chr1's window is already resident, so its remaining observations still come out
        // cleanly; the failure lands at chr2's fetch. The error's progress must be the LAST
        // clean observation's reach, whichever that turns out to be.
        let mut last_reach = first.reach_position();
        let error = loop {
            match gatherer
                .next()
                .expect("the walk must end in chr2's failure, not exhaustion")
            {
                Ok(observation) => last_reach = observation.reach_position(),
                Err(error) => break error,
            }
        };
        match error {
            RunError::SourceFailed {
                sample, reached, ..
            } => {
                assert_eq!(sample, "NA12878");
                assert_eq!(reached, WalkProgress::After(last_reach));
            }
            other => panic!("expected SourceFailed, got {other:?}"),
        }
    }

    /// **The file is the walk** (smoke at the fixture scale — B2 is the full oracle): what
    /// `write_psp` puts on disk reads back record for record as the walk streamed, under the
    /// header the gatherer fixed plus the one parameter the store adds at `create`.
    #[test]
    fn write_psp_round_trips_the_walk() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let paths = vec![bam_a];
        let psp_dir = TempDir::new().expect("a scratch dir");
        let psp_path = psp_dir.path().join("NA12878.psp");

        let (stats, counts) = open_gatherer(&paths, &reference)
            .expect("opens")
            .write_psp(&psp_path)
            .expect("the walk writes");

        let expected: Vec<SampleLocusObservations> = open_gatherer(&paths, &reference)
            .expect("opens again")
            .collect::<Result<_, _>>()
            .expect("a clean walk");
        assert!(
            !expected.is_empty(),
            "the fixture must yield observations, or this test proves an empty file",
        );
        assert_eq!(stats.records, expected.len() as u64);
        assert_eq!(
            counts.loci_emitted, stats.records,
            "the tally handed back is the walk's: one locus emitted per record written",
        );

        let mut reader = PspReader::open(&psp_path).expect("a finished psp opens");
        let mut expected_header = open_gatherer(&paths, &reference)
            .expect("opens again")
            .header()
            .clone();
        expected_header.writer.parameters.insert(
            "zstd-compression-level".to_string(),
            ParameterValue::Integer(i64::from(ZSTD_COMPRESSION_LEVEL)),
        );
        assert_eq!(reader.header(), &expected_header);

        let read_back: Vec<SampleLocusObservations> = reader
            .records()
            .expect("the walk starts")
            .map(|streamed| {
                streamed
                    .expect("a clean read")
                    .record
                    .expect("records() builds every body")
            })
            .collect();
        assert_eq!(read_back, expected);
    }

    /// **Gathering the same sample twice gives byte-identical files** (spec §12.1) — the
    /// walk is serial and deterministic and nothing in the gatherer reads a clock, so the
    /// file is a pure function of its inputs. The timestamp §12.1 exempts is the caller's
    /// to supply; with the fixture's fixed one, identity holds over the *whole* file,
    /// timestamp included — a stronger check than the exemption needs.
    #[test]
    fn gathering_the_same_sample_twice_gives_identical_bytes() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let paths = vec![bam_a];
        let psp_dir = TempDir::new().expect("a scratch dir");

        let first_psp = psp_dir.path().join("first.psp");
        let second_psp = psp_dir.path().join("second.psp");
        let (first_stats, _) = open_gatherer(&paths, &reference)
            .expect("opens")
            .write_psp(&first_psp)
            .expect("the first gather writes");
        let (second_stats, _) = open_gatherer(&paths, &reference)
            .expect("opens again")
            .write_psp(&second_psp)
            .expect("the second gather writes");

        assert!(
            first_stats.records > 0,
            "the fixture must yield records, or identical empty files prove nothing",
        );
        assert_eq!(first_stats.records, second_stats.records);
        let first_bytes = std::fs::read(&first_psp).expect("read the first file");
        let second_bytes = std::fs::read(&second_psp).expect("read the second file");
        assert_eq!(
            first_bytes, second_bytes,
            "two gathers of one sample must be the same bytes",
        );

        // **And the timestamp is the caller's**, which is what makes the comparison above
        // mean anything: a gatherer that stamped its own clock instead would pass it
        // whenever two gathers land in the same second, and produce irreproducible files
        // the moment they straddle one.
        let reader = PspReader::open(&first_psp).expect("the first file opens");
        assert_eq!(
            reader.header().writer.created,
            provenance().created,
            "the file must carry the caller's timestamp, not one the gatherer minted",
        );
    }

    /// A walk that fails mid-file must surface the failure, not seal the psp: a swallowed
    /// error here would leave a footer-complete file every reader accepts, holding none of
    /// the sample's records — the one failure shape the format cannot catch after the fact.
    #[test]
    fn write_psp_propagates_a_walk_failure() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let gatherer = open_gatherer(&[bam_a], &reference).expect("opens");
        // The FASTA is opened lazily per contig, so deleting it now fails the first draw.
        let fasta = reference
            .info()
            .fasta_path
            .clone()
            .expect("the digested fixture has bases");
        std::fs::remove_file(&fasta).expect("remove the fixture FASTA");

        let psp_dir = TempDir::new().expect("a scratch dir");
        let psp_path = psp_dir.path().join("NA12878.psp");
        let result = gatherer.write_psp(&psp_path);
        assert!(
            matches!(result, Err(RunError::SourceFailed { .. })),
            "{result:?}",
        );
        assert!(
            PspReader::open(&psp_path).is_err(),
            "the interrupted file must not read back as whole",
        );
    }

    /// A psp path that cannot be created is refused with the path named.
    #[test]
    fn write_psp_names_the_path_when_the_file_cannot_be_created() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let dir = TempDir::new().expect("scratch");
        let unwritable = dir.path().join("missing").join("NA12878.psp");

        let error = open_gatherer(&[bam_a], &reference)
            .expect("opens")
            .write_psp(&unwritable)
            .expect_err("a missing parent directory refuses the create");
        match error {
            RunError::PspNotWritten { path, .. } => assert_eq!(path, unwritable),
            other => panic!("expected PspNotWritten, got {other:?}"),
        }
    }

    /// One walk, one sample: files naming two samples are refused, and the refusal pairs
    /// **each file with the sample it claims** — the pairing, not just the names, is what
    /// lets an operator pick the stray file out of a hand-assembled list.
    #[test]
    fn files_naming_two_samples_are_refused() {
        let (_reference_dir, reference) = fixture_reference(true);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let (_b_dir, bam_b) = second_bam_with_library("NA12891", "libB", "b.bam");

        let error = open_gatherer(&[bam_a.clone(), bam_b.clone()], &reference)
            .expect_err("two samples' files must not open as one walk");

        assert!(
            matches!(error, RunError::FilesNotFromOneSample { .. }),
            "{error:?}"
        );
        let rendered = format_error_chain(&error);
        assert!(
            rendered.contains(&format!("'{}' names 'NA12878'", bam_a.display())),
            "{rendered}",
        );
        assert!(
            rendered.contains(&format!("'{}' names 'NA12891'", bam_b.display())),
            "{rendered}",
        );
    }

    /// The other route into the same refusal: files whose read-group headers cannot be read
    /// at all — here, a path that does not exist.
    #[test]
    fn open_refuses_files_whose_headers_cannot_be_read() {
        let (_reference_dir, reference) = fixture_reference(true);
        let missing = PathBuf::from("/no/such/place/ghost.bam");

        let error = open_gatherer(std::slice::from_ref(&missing), &reference)
            .expect_err("an unreadable file cannot establish a sample");
        assert!(
            matches!(error, RunError::FilesNotFromOneSample { .. }),
            "{error:?}"
        );
        let rendered = format_error_chain(&error);
        assert!(rendered.contains("ghost.bam"), "{rendered}");
    }

    /// A sample whose headers read fine but whose files will not open — a missing index the
    /// walk was told not to build — is refused naming the sample.
    #[test]
    fn open_names_the_sample_when_a_file_will_not_open() {
        let (_reference_dir, reference) = fixture_reference(true);
        let stem_reads = vec![read_named_with_length_in_read_group(
            "u-r0", 0, 5, 30, "rg1",
        )];
        let sam_header = header_with_read_groups(
            Some("coordinate"),
            &matching_contigs(),
            &[FixtureReadGroup {
                id: "rg1",
                sample: Some("NA12878"),
                library: Some("libA"),
                platform: None,
            }],
        );
        let (_dir, unindexed) = named_bam(&sam_header, &stem_reads, "unindexed.bam");

        let error = open_gatherer(&[unindexed], &reference)
            .expect_err("a missing index the walk may not build refuses the open");
        match error {
            RunError::OpeningSample { sample, .. } => assert_eq!(sample, "NA12878"),
            other => panic!("expected OpeningSample, got {other:?}"),
        }
    }

    /// **Analysed-but-empty ground round-trips** (spec §12.9): a walk over ground with no
    /// reads writes a psp with zero records, and the file still knows what ground was
    /// analysed — which is what lets a reader tell "analysed, nothing there" from "never
    /// looked", the distinction the header's analysed regions exist for.
    #[test]
    fn analysed_but_empty_ground_round_trips() {
        let (_reference_dir, reference) = fixture_reference(true);
        // Every read the fixture BAM holds is on chr1; the analysed ground is chr2 alone.
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let analysed = GenomeRegions::from_normalized_spans(
            vec![crate::regions::Region {
                chrom_id: 1,
                start: 1,
                end: 200,
            }],
            &fixture_bounds(),
        )
        .expect("one whole-contig span is normalized");
        let chr2_only_ground = build_segmentation(vec![generic_segment(1, 200)], analysed.clone());

        let psp_dir = TempDir::new().expect("a scratch dir");
        let psp_path = psp_dir.path().join("NA12878.psp");
        let (stats, counts) =
            open_gatherer_over(std::slice::from_ref(&bam_a), &reference, chr2_only_ground)
                .expect("opens")
                .write_psp(&psp_path)
                .expect("an empty walk still writes a whole file");

        assert_eq!(stats.records, 0);
        assert_eq!(counts.regions_in, 1, "chr2's one segment was dispatched");
        assert_eq!(
            counts.regions_handled, 1,
            "and it was genuinely handed to a generator rather than skipped for holding \
             no reads — the psp is byte-identical either way, so only the tally can tell \
             \"analysed, nothing there\" from \"never looked\"",
        );
        assert_eq!(counts.loci_emitted, 0);

        let mut reader = PspReader::open(&psp_path).expect("an empty psp is whole, not refused");
        assert_eq!(
            reader.header().segmentation_inputs.analysed_regions,
            analysed,
            "the file records the ground that was analysed and found empty",
        );
        assert_eq!(
            reader
                .records()
                .expect("the walk over no blocks starts")
                .count(),
            0,
        );
    }

    /// **Ground with a repeat tract in it gathers and round-trips** — the segmentation
    /// routes a tract segment to the tract generator, and whatever the walk yields over it,
    /// the file reads back equal to the walk. The guard asserts tract-kind records are
    /// actually in the stream, so this cannot silently degrade into a generic-only test.
    ///
    /// **Limitation: the fixture reference is a homopolymer**, so this round trip cannot see
    /// a store defect that is left/right-symmetric — a flank swap, a byte reversal — because
    /// the tract's two flanks are identical bytes here (measured: a mutant swapping them on
    /// write passed this test). That class is pinned where the sequence differs:
    /// `psp::record`'s `every_locus_kind_round_trips` uses asymmetric flanks, and
    /// `examples/ng_psp_gather_oracle.rs` runs on real sequence.
    #[test]
    fn tract_bearing_ground_round_trips_the_walk() {
        use crate::ng::locus_generation::LocusKind;
        use crate::ng::region_typing::segment_criteria::SsrSegment;
        use crate::ng::types::Motif;

        let (_reference_dir, reference) = fixture_reference(true);
        // Reads crossing the declared tract at 41-52: one anchored before it, one over it,
        // one after — the middle read is what gives the tract generator its evidence.
        let records = vec![
            read_named_with_length_in_read_group("t-r0", 0, 10, 30, "rg1"),
            read_named_with_length_in_read_group("t-r1", 0, 35, 30, "rg1"),
            read_named_with_length_in_read_group("t-r2", 0, 60, 30, "rg1"),
        ];
        let (_bam_dir, bam) = bam_in_read_group("NA12878", "libA", "rg1", "tract.bam", records);

        let chr1 = |start: u64, end: u64| GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        };
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
        let tract_ground = build_segmentation(
            segments,
            GenomeRegions::whole_contigs(&fixture_bounds()[..1]),
        );
        let open_over_tract_ground = || {
            open_gatherer_over(
                std::slice::from_ref(&bam),
                &reference,
                Arc::clone(&tract_ground),
            )
            .expect("opens")
        };

        let psp_dir = TempDir::new().expect("a scratch dir");
        let psp_path = psp_dir.path().join("NA12878.psp");
        let (stats, _counts) = open_over_tract_ground()
            .write_psp(&psp_path)
            .expect("the tract walk writes");

        let walked: Vec<SampleLocusObservations> = open_over_tract_ground()
            .collect::<Result<_, _>>()
            .expect("a clean walk");
        assert!(
            walked
                .iter()
                .any(|observation| matches!(observation.kind, LocusKind::Ssr(_))),
            "the tract must put tract-kind records on the walk, or this fixture is \
             generic-only and proves nothing about the tract path",
        );
        assert_eq!(stats.records, walked.len() as u64);

        // **And over tract ground too, the gatherer is the bare chain** — the equality the
        // real-data oracle's in-memory side rests on, which until now was pinned only over
        // generic segments.
        assert_eq!(
            walked,
            direct_walk(
                std::slice::from_ref(&bam),
                &reference,
                Arc::clone(&tract_ground),
                unusual_read_filters(),
                unusual_locus_generator_settings(),
            ),
        );

        let read_back: Vec<SampleLocusObservations> = PspReader::open(&psp_path)
            .expect("a finished psp opens")
            .records()
            .expect("the walk starts")
            .map(|streamed| {
                streamed
                    .expect("a clean read")
                    .record
                    .expect("records() builds every body")
            })
            .collect();
        assert_eq!(read_back, walked);
    }

    #[test]
    fn an_empty_file_list_is_refused() {
        let (_reference_dir, reference) = fixture_reference(true);
        let error = open_gatherer(&[], &reference).expect_err("no files, no walk");
        assert!(matches!(error, RunError::NoAlignmentFiles), "{error:?}");
    }

    /// A `.fai`-only reference holds no bases, and a walk cannot fetch loci from geometry.
    #[test]
    fn a_reference_without_bases_is_refused() {
        let (_reference_dir, reference) = fixture_reference(false);
        let (_a_dir, bam_a) = bam_with_library("NA12878", "libA", "a.bam");
        let error = open_gatherer(&[bam_a], &reference)
            .expect_err("a bases-less reference must refuse the walk");
        assert!(matches!(error, RunError::ReferenceHasNoBases), "{error:?}");
    }
}
