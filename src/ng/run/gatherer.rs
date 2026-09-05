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

use crate::fasta::ContigList;
use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
use crate::ng::locus_generation::ssr::DEFAULT_SSR_MAX_READS_PER_LOCUS;
use crate::ng::locus_generation::{
    LocusCounts, SampleLocusObservations, SampleLocusObservationsIterator,
};
use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use crate::ng::parameter_estimation::joint::census::{CensusWriter, DepthCap, ReadCap};
use crate::ng::parameter_estimation::joint::census_file::{PileupIdentity, write_census};
use crate::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLoci, ReferenceDigest, RegionSetDigest, SelectableRegions,
    SelectionError, SelectionTerms, select_kept_loci,
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
use crate::ng::region_typing::GenomeRegions;
use crate::ng::region_typing::RegionKind;
use crate::ng::repeat_catalog::{RepeatCatalog, StrRepeatCriteria};
use crate::ng::types::{Bp, ContigId};

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
/// **The census is fed at the yield point** (arch §3.3, spec §5.2), where a locus is handed
/// over — so the alignment files are read exactly once and produce both files. A gatherer
/// opened without a [`CensusPlan`] feeds nothing and writes nothing, which is what every test
/// that is about the psp alone wants.
pub struct SampleObservationGatherer {
    /// The header the sample's psp will carry — everything spec §6.1 asks for, fixed
    /// before the first record. The sample's name lives here (`header.sample`), not in a
    /// second field: one struct, one copy of the fact.
    header: Header,
    /// How far the walk has got — the second half of locating a failure (spec §9).
    reached: WalkProgress,
    /// **The census being accumulated, or nothing** — fed at the yield point below.
    census: Option<CensusWriter>,
    loci: SampleLocusObservationsIterator<RunSegments>,
}

/// **What a run needs to build a census beside each sample's psp**, shared by every sample of
/// the run.
///
/// **The selection is the run's and not the sample's**, which is why this is one value handed to
/// every gatherer rather than something each one computes. Which positions and which tracts are
/// kept is a function of the seed, the reference, the analysed regions and the catalog — none of
/// them per-sample — so two samples that selected separately would be selecting the same set
/// twice, and a run whose samples selected *differently* could not be fitted at all
/// (`parameter_prepass_census_sites.md` §3).
///
/// **Cheap to clone**: the loci and the contigs are behind handles, and the rest is a handful of
/// numbers and digests. `generate-psps` builds one and hands it to every sample's walk.
#[derive(Clone)]
pub struct CensusPlan {
    /// Which positions and which tracts this run keeps.
    pub loci: Arc<CensusLoci>,
    /// What the selection was made under — recorded in every census file, and what a later fit
    /// compares before pooling two samples.
    pub terms: SelectionTerms,
    /// The run's contigs, so that a tract named in the catalog becomes an id in the records.
    pub contigs: Arc<ContigList>,
    /// The most reads one repeat tract's evidence counts.
    pub read_cap: ReadCap,
    /// Where a position's reads stop being counted one by one and its allele counts are thinned
    /// proportionally.
    pub depth_cap: DepthCap,
}

/// **What a run chooses about its census**, as a person or a command supplies it.
///
/// Two of the three numbers are settled by the design and carry their figure here; the third,
/// the seed, is what makes one run's selection a different set from another's.
#[derive(Debug, Clone, Copy)]
pub struct CensusSelection {
    /// **The run's selection seed.** Every sample of a cohort must use the same one: which
    /// positions are kept is `hash(contig, position, seed) < threshold`, so two invocations that
    /// seeded differently keep disjoint sets and their samples cannot be pooled
    /// (`parameter_prepass_census_sites.md` §3). It travels in the census file, so the
    /// disagreement is named at the fit rather than averaged over.
    pub seed: u64,
    /// **How many generic positions the run aims to keep.** A budget, not a threshold: the
    /// threshold is computed from it and how much ground there is.
    pub generic_target: u64,
    /// **How many repeat tracts a stratum keeps at most.**
    pub ssr_cap: usize,
}

impl CensusSelection {
    /// **About two million positions, five thousand tracts a stratum, and a fixed seed.**
    ///
    /// The two counts are the design's own figures rather than numbers chosen here:
    /// *"the knob is a number of positions, and about two million of them is the default"*
    /// (`parameter_prepass_census_sites.md` §5.1), sized to yield about ten thousand
    /// segregating sites; and five thousand tracts a stratum is where
    /// `parameter_prepass_joint_loci.md` §6's first question closed, measured on a tomato
    /// archive.
    ///
    /// **The seed is a constant and not a clock**, which is what makes psp mode work at all: a
    /// cohort is walked by separate invocations, so a seed that differed between them would keep
    /// disjoint sets of positions and the samples could not be pooled. Two invocations agree
    /// here by construction rather than by somebody typing the same number twice.
    pub const SHIPPED: Self = Self {
        seed: 0x5EED_C0FF_EE15_0000,
        generic_target: 2_000_000,
        ssr_cap: 5_000,
    };
}

impl std::fmt::Debug for CensusPlan {
    /// **The sizes, not the sets.** A derived `Debug` would print every kept position — two
    /// million of them on a real run, in a message somebody is reading to find out which run
    /// this is.
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("CensusPlan")
            .field("generic_positions", &self.loci.generic().len())
            .field("ssr_loci", &self.loci.ssr_stratum_counts().total())
            .field("seed", &self.terms.seed)
            .finish_non_exhaustive()
    }
}

impl CensusPlan {
    /// **Choose this run's census positions and tracts, once, for every sample of it.**
    ///
    /// The selection is a pure function of the seed, the reference, the analysed regions and the
    /// catalog — none of them per-sample — which is why it is made here and shared rather than
    /// recomputed per walk (`parameter_prepass_census_sites.md` §3).
    ///
    /// **`unambiguous` is where the reference is sequence at all**, the runs of `A`, `C`, `G` and
    /// `T` collected while the FASTA was read. Generic positions are drawn from the analysed
    /// ground *intersected* with it: a position in a run of `N` has no reference base to compare
    /// a read against, and keeping one would put a permanent hole in every sample's records.
    ///
    /// # Errors
    ///
    /// [`RunError::CensusNotPlanned`] when the reference carries no digest to bind the selection
    /// to, when the analysed regions overlap, or when the catalog cannot serve the criteria this
    /// run routed under.
    pub fn of_run(
        selection: CensusSelection,
        catalog: &RepeatCatalog,
        analysed: &GenomeRegions,
        unambiguous: &SelectableRegions,
        reference: &ReferenceInfo,
        criteria: &StrRepeatCriteria,
    ) -> Result<Self, RunError> {
        let planned = |source: SelectionError| RunError::CensusNotPlanned {
            source: Box::new(source),
        };
        let analysed = SelectableRegions::new(analysed.iter().collect()).map_err(planned)?;
        let terms = SelectionTerms {
            seed: selection.seed,
            reference: ReferenceDigest::of(reference).map_err(planned)?,
            analysed_regions: RegionSetDigest::of(&analysed),
            catalog_built_under: CatalogBuildSettings::of(catalog),
            ssr_criteria: criteria.clone(),
            generic_target: selection.generic_target,
            ssr_cap: selection.ssr_cap,
        };
        let loci = select_kept_loci(&terms, catalog, &analysed, unambiguous).map_err(planned)?;
        Ok(Self {
            loci: Arc::new(loci),
            terms,
            contigs: Arc::new(reference.contig_list()),
            // **The tract read cap is the locus generator's own**, so the census counts the
            // reads the walk counted rather than a second number nobody set.
            read_cap: ReadCap(DEFAULT_SSR_MAX_READS_PER_LOCUS),
            depth_cap: DEPTH_CAP,
        })
    }
}

/// **Where a position's reads stop being counted one by one**, and its allele counts are thinned
/// to that many proportionally.
///
/// **Not the depth ladder's top, and the two are separate knobs.** The depth code records what
/// the position held, all the way to the ladder's ceiling near 1,500; this is where the *reads*
/// stop being enumerated, and the fractions the alleles showed survive the thinning.
const DEPTH_CAP: DepthCap = DepthCap::new(124);

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
        census: Option<&CensusPlan>,
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
        // **Built before the walk, and every generic stretch marked walked before a locus
        // arrives.** Without the marking a position no read reached is indistinguishable from a
        // region the run never opened, because the generic generator emits no locus where there
        // is no read — measured on tomato SRR7279482 at 25×, 1 in 21 kept positions, every one
        // of them data reported as a defect (`census.rs`'s own note on `mark_walked`).
        let census = census.map(|plan| {
            let contigs = Arc::clone(&plan.contigs);
            let contig_of = move |name: &str| {
                contigs
                    .entries
                    .iter()
                    .position(|entry| entry.name == name)
                    .map(|index| ContigId(index as u32))
            };
            let mut writer = CensusWriter::new(
                header.sample.clone(),
                &plan.loci,
                read_groups.iter().map(|(id, _)| id).collect(),
                &contig_of,
                plan.terms.clone(),
                // **The census's own ladder, by name.** `DepthBinEdges::for_census` is the one
                // ladder a census is recorded on — exact depths to 124 and ten widening rungs
                // above — so it is called here rather than carried in the plan, where it would
                // read as a knob a run may set and is not.
                DepthBinEdges::for_census(),
                plan.read_cap,
                plan.depth_cap,
            );
            for region in segmentation.segments() {
                if region.kind == RegionKind::Generic {
                    writer.mark_walked(region.region);
                }
            }
            writer
        });
        Ok(Self {
            header,
            reached: WalkProgress::NothingYet,
            census,
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
    /// **⚠ An existing file at `path` is truncated** — [`PspWriter::create`]'s contract, so a
    /// caller that must not destroy a finished psp has to keep it out of the way itself.
    /// `generate-psps` writes each sample beside its psp and renames once the file is whole,
    /// which is what stops a stopped re-walk from destroying the psp it was replacing;
    /// refusing an overwrite outright is plan step C3's.
    ///
    /// # Errors
    ///
    /// A walk failure names the sample and how far it got ([`RunError::SourceFailed`]); a
    /// record the store refuses names its locus ([`RunError::RecordNotWritten`]); a file
    /// that cannot be created or sealed names its path ([`RunError::PspNotWritten`]). After
    /// any of these the file at `path` is not whole, and a reader will refuse it as
    /// interrupted rather than read it as complete (`psp_file_format.md` §10).
    pub fn write_psp(
        mut self,
        path: &Path,
        census: Option<&Path>,
    ) -> Result<(WriteStats, LocusCounts), RunError> {
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
        self.write_census_beside(census, &stats)?;
        Ok((stats, counts))
    }

    /// **The census file, written once the psp is whole**, and named by that psp.
    ///
    /// **The identity is built from the psp's own header and its record count**, which is this
    /// value's first real construction site: `PileupIdentity::of_header` takes the count as its
    /// own argument precisely because a header written before the first record cannot carry one
    /// (spec §6.1's ruling of 2026-09-03). So it is built here, after `finish`, from the
    /// [`WriteStats`] that call returns.
    ///
    /// **A census asked for and not written fails the sample's walk** (spec §2, plan step G2).
    /// The two files are one product: a psp whose census is missing forces the sample to be
    /// walked again, which is the one thing psp mode exists to avoid — so a run that could not
    /// write it says so rather than leaving a psp that looks finished.
    ///
    /// **Asking for a census from a gatherer that was not given a plan is a defect and says
    /// so.** The pair is settled at `open`: a command that names a census path has already
    /// handed over the plan that fills it, and silently writing nothing there would leave an
    /// absent file that reads exactly like a walk that was never asked.
    fn write_census_beside(self, path: Option<&Path>, stats: &WriteStats) -> Result<(), RunError> {
        let (Some(path), Some(census)) = (path, self.census) else {
            assert!(
                path.is_none(),
                "a census was asked for at {} and this walk was opened without a plan to build \
                 one. This is a defect in ng rather than anything about the data.",
                path.map(Path::display)
                    .map(|it| it.to_string())
                    .unwrap_or_default(),
            );
            return Ok(());
        };
        // **The digest the writer handed back, not one taken from this gatherer's own header.**
        // `PspWriter::create` records the compression level into the header before encoding it,
        // so the header this object holds is not the header in the file — measured, and the
        // difference is one line of TOML that changes every byte of the digest. A census naming
        // a header no psp carries would make every freshness check say *rebuild* for ever, which
        // is exactly the failure naming the pileup exists to prevent.
        let identity = PileupIdentity {
            header: stats.header_digest,
            records: stats.records,
        };
        let evidence = census.finish();
        let mut file =
            std::fs::File::create(path).map_err(|source| RunError::CensusNotWritten {
                path: path.to_path_buf(),
                source: Box::new(source),
            })?;
        write_census(&evidence, Some(identity), &mut file).map_err(|source| {
            RunError::CensusNotWritten {
                path: path.to_path_buf(),
                source: Box::new(source),
            }
        })
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
            census,
            loci: _,
        } = self;
        formatter
            .debug_struct("SampleObservationGatherer")
            .field("sample", &header.sample)
            .field("reached", reached)
            .field("census", &census.is_some())
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
                // **The census is fed here, at the yield point, and not by the psp writer**
                // (arch §3.3, spec §5.2). Feeding it from the writer would make the census a
                // function of what was *stored* rather than of what the walk *saw*, and would
                // leave a consumer that iterates a gatherer without writing a psp — which the
                // suite does — building no census at all. This is also what makes spec §2's
                // promise true: the alignment files are read once and produce both files.
                if let Some(census) = self.census.as_mut() {
                    census.add_locus(&observation);
                }
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
    use crate::ng::run::test_fixtures::{build_segmentation, index, unusual_read_filters};
    use crate::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId};
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
    pub(super) fn provenance() -> WriterProvenance {
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
            None,
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
            .write_psp(&psp_path, None)
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
            .write_psp(&first_psp, None)
            .expect("the first gather writes");
        let (second_stats, _) = open_gatherer(&paths, &reference)
            .expect("opens again")
            .write_psp(&second_psp, None)
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
        let result = gatherer.write_psp(&psp_path, None);
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
            .write_psp(&unwritable, None)
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
                .write_psp(&psp_path, None)
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
            .write_psp(&psp_path, None)
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

#[cfg(test)]
mod census_tests {
    //! **The census beside the psp** (plan step G1): fed at the yield point, written once the
    //! psp is whole, and named by that psp.
    //!
    //! **These use the binary namespace's on-disk cohort fixture**, which is the only one in the
    //! tree that builds a real catalog file — and a real catalog is what the selection needs.
    //! The alternative was a fourth copy of the reference-and-catalog builder, which is the debt
    //! Milestone C recorded and F1 paid off.

    use super::*;
    use crate::ng::parameter_estimation::joint::census_file::open_census;
    use crate::ng::psp::PspReader;
    use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
    use crate::ng::region_typing::segment_criteria::{
        DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
    };
    use crate::pop_var_caller_exp::run_ground::{self, GroundRequest, RepeatRouting};
    use crate::pop_var_caller_exp::test_fixtures::{ACohortOnDisk, a_cohort_on_disk};

    /// The fixture cohort, its ground, and a census plan over it.
    fn a_cohort_with_a_census_plan() -> (ACohortOnDisk, Arc<Segmentation>, CensusPlan) {
        use crate::ng::parameter_estimation::joint::loci::UnambiguousRuns;
        use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
        use crate::ng::repeat_catalog::RepeatCatalog;

        let cohort = a_cohort_on_disk();
        // **The reference is read with an observer**, because the selection needs to know where
        // the genome is sequence at all: a position inside a run of `N` has no reference base to
        // compare a read against.
        let mut callable = UnambiguousRuns::default();
        let reference = read_reference_info_observing(
            ReferenceSource::Fasta {
                fasta: cohort.reference.clone(),
                fai: None,
            },
            &mut callable,
        )
        .expect("the fixture's reference reads");
        let unambiguous = callable
            .into_selectable()
            .expect("maximal runs are disjoint");

        let request = GroundRequest {
            reference: &cohort.reference,
            catalog: Some(&cohort.catalog),
            regions: None,
            routing: RepeatRouting {
                min_copies: MinCopies::default(),
                min_period: DEFAULT_MIN_PERIOD,
                max_period: DEFAULT_MAX_PERIOD,
                max_str_len: DEFAULT_MAX_STR_LEN,
                min_purity: DEFAULT_MIN_PURITY,
            },
        };
        let analysed = run_ground::analysed_regions(&request, &reference.contig_list())
            .expect("the whole reference");
        let segmentation =
            Arc::new(run_ground::segments_over(&request, &analysed, &reference).expect("it types"));
        let catalog = RepeatCatalog::open_checking_against_reference(&cohort.catalog, &reference)
            .expect("the fixture's catalog is this reference's");
        // **A target of one position per base, so the fixture keeps some.** The shipped budget
        // is two million positions and this genome is 300 bases; at the shipped number the
        // threshold keeps everything, which is what a test of the wiring wants — but saying so
        // is better than relying on it.
        let plan = CensusPlan::of_run(
            CensusSelection {
                generic_target: 300,
                ..CensusSelection::SHIPPED
            },
            &catalog,
            &analysed,
            &unambiguous,
            &reference,
            &segmentation.inputs().repeat_tract_criteria,
        )
        .expect("the fixture's ground can be selected from");
        (cohort, segmentation, plan)
    }

    /// Open a gatherer over one of the fixture's samples.
    fn gatherer_over(
        cohort: &ACohortOnDisk,
        which: usize,
        segmentation: &Arc<Segmentation>,
        plan: Option<&CensusPlan>,
    ) -> SampleObservationGatherer {
        use crate::ng::read::input::reference::OpenReference;
        use crate::ng::reference_info::{
            ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
        };

        let cache = Arc::new(ReferenceInfoCache::new());
        let (info, _) = read_reference_verifying_or_creating_fai(
            &cache,
            cohort.reference.clone(),
            ReferenceCheck::TrustIndexWithoutChecking,
        )
        .expect("the reference reads");
        let reference = OpenReference::new(info);
        SampleObservationGatherer::open(
            SampleWalkInputs {
                alignments: std::slice::from_ref(&cohort.alignments[which]),
                reference: &reference,
                read_filters: ReadFilterConfig::default(),
                locus_generator_settings: PileupGeneratorConfig::default(),
                build_index_if_missing: false,
            },
            Arc::clone(segmentation),
            super::tests::provenance(),
            plan,
        )
        .expect("the sample opens")
    }

    /// **The walk writes both files, and the census names the psp it was built from.**
    ///
    /// The identity is the psp header's digest and its record count, which is what a later fit
    /// compares before trusting a census: two censuses built from different reads are otherwise
    /// indistinguishable, and that is the whole reason the pileup is named
    /// (`census_file.rs`'s `Freshness`).
    ///
    /// **The expected value is rebuilt from the psp on disk**, which is the check that matters:
    /// the writer amends the header before writing it, so a census built from the header the
    /// gatherer *holds* names a file that does not exist. This test failed exactly that way
    /// before `WriteStats` began carrying the digest of what was written.
    #[test]
    fn the_walk_writes_a_census_beside_the_psp_and_names_that_psp_in_it() {
        let (cohort, segmentation, plan) = a_cohort_with_a_census_plan();
        let psp = cohort.directory.path().join("zeta.psp");
        let census = cohort.directory.path().join("zeta.census");

        let (stats, _) = gatherer_over(&cohort, 0, &segmentation, Some(&plan))
            .write_psp(&psp, Some(&census))
            .expect("the walk writes both files");

        assert!(census.is_file(), "the census is at {census:?}");
        let (_evidence, named) = open_census(&census).expect("this build's own census");
        let named = named.expect("the census names the psp it was built from");
        let expected = PileupIdentity::of_header(
            &PspReader::open(&psp)
                .expect("the psp opens")
                .header()
                .encode()
                .expect("the header it was written with re-encodes"),
            stats.records,
        );
        assert_eq!(
            named, expected,
            "the identity is the psp's own header and its own record count",
        );
    }

    /// **A walk opened without a plan writes no census**, which is what every test about the psp
    /// alone wants — and what a run that has not asked for one must do.
    #[test]
    fn a_walk_with_no_census_plan_writes_no_census() {
        let (cohort, segmentation, _plan) = a_cohort_with_a_census_plan();
        let psp = cohort.directory.path().join("zeta.psp");

        let _ = gatherer_over(&cohort, 0, &segmentation, None)
            .write_psp(&psp, None)
            .expect("the walk writes its psp");

        assert!(psp.is_file());
        assert!(
            !cohort.directory.path().join("zeta.census").exists(),
            "nothing asked for a census, so none was written",
        );
    }

    /// **A census that cannot be written fails the sample's walk** (spec §2, plan step G2).
    ///
    /// The two files are one product: a psp whose census is missing forces the sample to be
    /// walked again, which is what psp mode exists to avoid. A run that reported this and
    /// carried on would be storing that re-walk for somebody to find later.
    #[test]
    fn a_census_that_cannot_be_written_fails_the_walk() {
        let (cohort, segmentation, plan) = a_cohort_with_a_census_plan();
        let psp = cohort.directory.path().join("zeta.psp");
        let census = cohort
            .directory
            .path()
            .join("no-such-directory")
            .join("zeta.census");

        let refused = gatherer_over(&cohort, 0, &segmentation, Some(&plan))
            .write_psp(&psp, Some(&census))
            .expect_err("there is nowhere to write the census");

        let rendered = crate::error_render::format_error_chain(&refused);
        assert!(
            rendered.contains("census") && rendered.contains("no-such-directory"),
            "the refusal names the census and where it would have gone, and got: {rendered}",
        );
    }

    /// **The census is fed by the walk, not by the psp writer** — so what it holds is what the
    /// sample showed, and a sample that showed nothing over the ground is told apart from ground
    /// nobody walked.
    ///
    /// `zeta` carries three reads and `alpha` none, over the same ground. Both censuses cover
    /// the same positions — the selection is the run's — and they must differ in what those
    /// positions say, or the accumulator is not being fed at all.
    #[test]
    fn what_the_census_holds_is_what_the_sample_showed() {
        let (cohort, segmentation, plan) = a_cohort_with_a_census_plan();
        let mut written = Vec::new();
        for (which, sample) in ["zeta", "alpha"].iter().enumerate() {
            let psp = cohort.directory.path().join(format!("{sample}.psp"));
            let census = cohort.directory.path().join(format!("{sample}.census"));
            let _ = gatherer_over(&cohort, which, &segmentation, Some(&plan))
                .write_psp(&psp, Some(&census))
                .expect("the walk writes both files");
            written.push(std::fs::read(&census).expect("the census reads"));
        }

        assert_ne!(
            written[0], written[1],
            "one sample carried reads and the other none, so their censuses cannot be the same \
             bytes — if they are, nothing is being accumulated",
        );
    }
}
