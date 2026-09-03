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

use std::path::PathBuf;
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
