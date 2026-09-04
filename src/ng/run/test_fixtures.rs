//! Fixtures every `ng::run` test module shares — the catalog header, the segmentation
//! scaffold, an alignment index, and the read filters that are deliberately not the default.
//!
//! **Why they live here rather than in each module's own `mod tests`.** Each of
//! [`callers`](super::callers), [`walker`](super::walker) and [`gatherer`](super::gatherer)
//! grew its own copy, because a `pub(super)` helper inside a `#[cfg(test)] mod tests`
//! resolves to *that module*, not to `run` — so a sibling file cannot reach it and starts
//! over. Three copies of two of these had already appeared before this module existed.
//!
//! **Two of the values here are load-bearing and must not be edited in one place alone**,
//! which is the whole reason for the hoist:
//!
//! - the catalog's `reference_md5` of `[7; 16]` — a run whose reference carries a checksum
//!   refuses a catalog built on a different one, so every fixture that hands over a
//!   checksummed reference has to hand over a catalog agreeing with it;
//! - [`unusual_read_filters`]' `MapQual(37)` — several tests prove a run *applied* what it
//!   was given by showing the result differs from the shipped default's, which is only a
//!   proof while this value is not the default.

use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::ng::read::ReadFilterConfig;
use crate::ng::region_typing::{GenomeRegions, TypedRegion};
use crate::ng::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};
use crate::ng::tandem_repeat::ScanParams;
use crate::ng::types::MapQual;

use super::Segmentation;

/// The catalog header every fixture segmentation is built under, claiming the digest
/// `[7; 16]` — the one a checksummed fixture reference must agree with.
pub(crate) fn catalog_header() -> RepeatCatalogHeader {
    RepeatCatalogHeader {
        contigs: Vec::new(),
        reference_md5: [7; 16],
        built_under: StrRepeatCriteria::default(),
        scan: ScanParams::default(),
        tool_version: "test".to_string(),
        longest_tract_bp: Vec::new(),
    }
}

/// The same, claiming `reference_md5` instead — for a test that hands over a reference whose
/// checksum is not the fixture's.
pub(crate) fn catalog_header_built_on(reference_md5: [u8; 16]) -> RepeatCatalogHeader {
    RepeatCatalogHeader {
        reference_md5,
        ..catalog_header()
    }
}

/// A segmentation over `segments`, analysing `analysed`, under [`catalog_header`] and the
/// default repeat criteria — **the one copy of that scaffold**, so a change to
/// [`Segmentation::build`]'s signature or to the catalog header's fields touches one site.
pub(crate) fn build_segmentation(
    segments: Vec<TypedRegion>,
    analysed: GenomeRegions,
) -> Arc<Segmentation> {
    build_segmentation_under(segments, analysed, catalog_header())
}

/// [`build_segmentation`], with the catalog header named — for the tests that are *about*
/// the catalog's identity.
pub(crate) fn build_segmentation_under(
    segments: Vec<TypedRegion>,
    analysed: GenomeRegions,
    catalog: RepeatCatalogHeader,
) -> Arc<Segmentation> {
    Arc::new(
        Segmentation::build(
            segments.into_iter().map(Ok),
            analysed,
            catalog,
            StrRepeatCriteria::default(),
            PathBuf::from("/genomes/test.catalog.parquet"),
        )
        .expect("a clean stream builds"),
    )
}

/// Build a `.bai`/`.crai` beside `path`, as a run that was asked to would.
pub(crate) fn index(path: &PathBuf) {
    crate::bam::index_preflight::preflight_alignment_indexes(std::slice::from_ref(path), true)
        .expect("build index");
}

/// Read filters that are **not** the shipped default, so a run that dropped what it was
/// given and applied the constants instead is visible.
///
/// `MapQual(37)` sits above the default floor of 20, which is what lets a fixture read at
/// MAPQ 30 separate the two policies.
pub(crate) fn unusual_read_filters() -> ReadFilterConfig {
    ReadFilterConfig {
        min_mapq: Some(MapQual(37)),
        ..ReadFilterConfig::default()
    }
}

// ---------------------------------------------------------------------
// The census plan, and a walk over one sample of the fixture cohort
// ---------------------------------------------------------------------

/// **The ground of a fixture cohort, and a census selection over it.**
///
/// Both producers of a census need this — the one that builds it while the reads are walked
/// ([`gatherer`](super::gatherer)) and the one that builds it afterwards from the stored psp
/// ([`census_from_psp`](super::census_from_psp)) — and the byte-for-byte comparison between
/// them is only a statement about the psp while both are selecting the same loci. So the
/// selection is made here, once, rather than in each test module.
///
/// **The reference is read with an observer**, because the selection has to know where the
/// genome is sequence at all: a position inside a run of `N` has no reference base to compare a
/// read against, and keeping one would leave a permanent hole in every sample's records.
///
/// **A target of one position per base**, so a 300-base fixture keeps some. The shipped budget
/// is two million positions; at that number the threshold keeps everything here, which is what
/// a test of the wiring wants — but saying so is better than relying on it.
pub(crate) fn a_census_plan_over(
    reference_fasta: &Path,
    catalog_path: &Path,
) -> (Arc<Segmentation>, super::CensusPlan) {
    use crate::ng::parameter_estimation::joint::loci::UnambiguousRuns;
    use crate::ng::reference_info::{ReferenceSource, read_reference_info_observing};
    use crate::ng::region_typing::DEFAULT_MAX_STR_LEN;
    use crate::ng::region_typing::segment_criteria::{
        DEFAULT_MAX_PERIOD, DEFAULT_MIN_PERIOD, DEFAULT_MIN_PURITY, MinCopies,
    };
    use crate::ng::repeat_catalog::RepeatCatalog;
    use crate::ng::run::{CensusPlan, CensusSelection};
    use crate::pop_var_caller_exp::run_ground::{self, GroundRequest, RepeatRouting};

    let mut callable = UnambiguousRuns::default();
    let reference = read_reference_info_observing(
        ReferenceSource::Fasta {
            fasta: reference_fasta.to_path_buf(),
            fai: None,
        },
        &mut callable,
    )
    .expect("the fixture's reference reads");
    let unambiguous = callable
        .into_selectable()
        .expect("maximal runs are disjoint");

    let request = GroundRequest {
        reference: reference_fasta,
        catalog: Some(catalog_path),
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
    let catalog = RepeatCatalog::open_checking_against_reference(catalog_path, &reference)
        .expect("the fixture's catalog is this reference's");
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
    (segmentation, plan)
}

/// A walk over one sample of a fixture cohort, with a census plan or without one.
pub(crate) fn gatherer_over(
    alignment: &PathBuf,
    reference: &Path,
    segmentation: &Arc<Segmentation>,
    plan: Option<&super::CensusPlan>,
) -> super::SampleObservationGatherer {
    use crate::ng::locus_generation::pileup::PileupGeneratorConfig;
    use crate::ng::read::ReadFilterConfig;
    use crate::ng::read::input::reference::OpenReference;
    use crate::ng::reference_info::{
        ReferenceCheck, ReferenceInfoCache, read_reference_verifying_or_creating_fai,
    };
    use crate::ng::run::{SampleObservationGatherer, SampleWalkInputs};

    let cache = Arc::new(ReferenceInfoCache::new());
    let (info, _) = read_reference_verifying_or_creating_fai(
        &cache,
        reference.to_path_buf(),
        ReferenceCheck::TrustIndexWithoutChecking,
    )
    .expect("the reference reads");
    let open_reference = OpenReference::new(info);
    SampleObservationGatherer::open(
        SampleWalkInputs {
            alignments: std::slice::from_ref(alignment),
            reference: &open_reference,
            read_filters: ReadFilterConfig::default(),
            locus_generator_settings: PileupGeneratorConfig::default(),
            build_index_if_missing: false,
        },
        Arc::clone(segmentation),
        walk_provenance(),
        plan,
    )
    .expect("the sample opens")
}

/// What a fixture walk records about the program that produced its psp.
pub(crate) fn walk_provenance() -> crate::ng::psp::WriterProvenance {
    use crate::ng::psp::{ParameterValue, WriterProvenance};
    use std::collections::BTreeMap;

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
