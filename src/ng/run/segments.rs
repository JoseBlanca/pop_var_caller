//! The ground a run analyses: its segments, and the record of what they were computed from.
//!
//! A **segment** is one stretch of the reference the typed-region generator has classified —
//! a repeat tract, or a stretch of ordinary sequence between tracts. Every sample of a run
//! walks the same segments in the same order, because they are a function of the reference,
//! the repeat catalog, the criteria the run asked with and the regions it was asked to
//! analyse, and of no sample's reads (`doc/devel/ng/spec/run_streaming.md` §4.2).
//!
//! That is why they are computed once and shared: k samples over one segmentation, never k
//! segmentations that happen to agree.

use std::path::PathBuf;

use crate::ng::region_typing::{GenomeRegions, TypedRegion};
use crate::ng::repeat_catalog::{RepeatCatalogError, RepeatCatalogHeader, StrRepeatCriteria};
use crate::ng::types::GenomeRegion;

use super::RunError;

/// Re-exported from its own module so `ng::run`'s call sites keep their import path.
/// [`SegmentationInputs`] moved to a top-level home when the psp header started carrying it
/// ([`crate::ng::psp`]) — an interchange type both stages hold belongs to neither
/// (`doc/devel/ng/impl_plan/run_driver_psp_mode.md`, Checkpoint A).
pub use crate::ng::segmentation_inputs::SegmentationInputs;

/// The run's segments, in genome order, beside the values they were computed from.
///
/// **Held whole rather than streamed**, because every sample needs the same list and the
/// samples do not advance through it together — the merge draws each one forward at its own
/// pace, so a stream would mean re-reading the catalog once per sample
/// (`doc/devel/ng/arch/run_streaming.md` §1, §2).
///
/// **What that costs is unmeasured, and it is the one term of a run that grows with the
/// genome rather than with the cohort.** A `TypedRegion` is a span, a kind, and — at a repeat
/// tract — an owned motif and contig name. How many segments a genome has at the criteria a
/// run asks with has never been counted; spec §11's question 1 is that measurement. Beside 63
/// open alignment files it is small; at one sample it may not be.
///
/// **Neither `Clone` nor a deriving `Debug`.** A clone would deep-copy the whole list for a
/// second holder that the design says should not exist, and a derived `Debug` would print
/// every segment.
pub struct Segmentation {
    inputs: SegmentationInputs,
    segments: Vec<TypedRegion>,
    /// The analysed regions as the merge takes them — a slice, where
    /// [`SegmentationInputs::analysed_regions`] is the set the compatibility checks compare.
    ///
    /// **Both spellings of one fact, built together here** so they cannot come apart: the walk
    /// advances over the segments and the merge over these, and a run whose two halves
    /// disagreed about the ground would call positions no sample walked.
    analysed_regions: Vec<GenomeRegion>,
}

impl Segmentation {
    /// Consume the typed-region generator's stream once, and record what produced it.
    ///
    /// **What the recording does and does not promise.** `analysed` is both recorded and used
    /// here, so the regions this reports and the regions the merge walks are the same value.
    /// `catalog` and `repeat_tract_criteria` are recorded as given: the caller read the
    /// segment stream, and only the caller can guarantee it asked the catalog with the
    /// criteria it passes here. Handing in one set and reading the stream with another would
    /// make every later compatibility check compare the wrong thing (spec §6.2), so the two
    /// belong in one call site.
    ///
    /// `catalog_path` is carried only so a failure can name the file. Four of the catalog's
    /// own failures name no path — a digest mismatch, over-permissive criteria, differing scan
    /// weights, a differing tool version — and a person with several catalogs on disk needs to
    /// know which one spoke.
    pub fn build(
        segments: impl Iterator<Item = Result<TypedRegion, RepeatCatalogError>>,
        analysed: GenomeRegions,
        catalog: RepeatCatalogHeader,
        repeat_tract_criteria: StrRepeatCriteria,
        catalog_path: PathBuf,
    ) -> Result<Self, RunError> {
        let segments =
            segments
                .collect::<Result<Vec<_>, _>>()
                .map_err(|source| RunError::Catalog {
                    path: catalog_path,
                    source,
                })?;
        let analysed_regions = analysed.iter().collect();
        Ok(Self {
            inputs: SegmentationInputs {
                catalog,
                repeat_tract_criteria,
                analysed_regions: analysed,
            },
            segments,
            analysed_regions,
        })
    }

    /// What this segmentation was computed from — the operand of every compatibility check
    /// (spec §6.2).
    #[must_use]
    pub fn inputs(&self) -> &SegmentationInputs {
        &self.inputs
    }

    /// The segments, in genome order. A segment never crosses a contig and is never cut
    /// (spec §4.3), which is the generator's guarantee and not re-checked here.
    #[must_use]
    pub fn segments(&self) -> &[TypedRegion] {
        &self.segments
    }

    /// The regions the run was asked to analyse, as the merge takes them.
    ///
    /// The merge advances over the *analysed regions*, not over the segments: it cuts them
    /// into its own building regions and asks each sample for what falls inside
    /// ([`merge_cohort_through_cache`](crate::ng::run::cohort_merge::serial::merge_cohort_through_cache)).
    /// A sample's walk advances over the segments. Both come out of this one object, so the
    /// two cannot be built from different ground.
    #[must_use]
    pub fn analysed_regions(&self) -> &[GenomeRegion] {
        &self.analysed_regions
    }
}

/// **The sizes, not the contents.** A derived `Debug` would print every segment — millions of
/// them on a genome — in a line someone is reading to find out what ground a run covers.
impl std::fmt::Debug for Segmentation {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("Segmentation")
            .field("segments", &self.segments.len())
            .field("analysed_regions", &self.analysed_regions.len())
            .finish_non_exhaustive()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::region_typing::RegionKind;
    use crate::ng::tandem_repeat::ScanParams;
    use crate::ng::types::{ContigId, Position};
    use crate::regions::ContigBounds;

    fn contigs() -> Vec<ContigBounds<'static>> {
        vec![
            ContigBounds {
                name: "chr1",
                length: 1_000,
            },
            ContigBounds {
                name: "chr2",
                length: 500,
            },
        ]
    }

    fn header() -> RepeatCatalogHeader {
        RepeatCatalogHeader {
            contigs: Vec::new(),
            reference_md5: [7; 16],
            built_under: StrRepeatCriteria::default(),
            scan: ScanParams::default(),
            tool_version: "test".to_string(),
            longest_tract_bp: Vec::new(),
        }
    }

    /// **Criteria that are not the default**, so a test asserting they were recorded cannot
    /// pass by the recorder substituting a default.
    fn unusual_criteria() -> StrRepeatCriteria {
        let mut criteria = StrRepeatCriteria::default();
        criteria.classification.min_purity = 0.93;
        criteria
    }

    fn catalog_path() -> PathBuf {
        PathBuf::from("/genomes/chosen.catalog.parquet")
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    fn generic(start: u64, end: u64) -> TypedRegion {
        TypedRegion {
            region: region(start, end),
            kind: RegionKind::Generic,
        }
    }

    /// The stream is consumed once and kept in the order it arrived.
    #[test]
    fn build_keeps_the_generators_segments_in_order() {
        let given = vec![generic(1, 10), generic(11, 20), generic(21, 30)];
        let segmentation = Segmentation::build(
            given.clone().into_iter().map(Ok),
            GenomeRegions::whole_contigs(&contigs()),
            header(),
            unusual_criteria(),
            catalog_path(),
        )
        .expect("a clean stream builds");

        assert_eq!(segmentation.segments(), given.as_slice());
    }

    /// What a segmentation reports about itself is what it was built from — the property
    /// every compatibility check rests on.
    ///
    /// **Every input here differs from its type's default**, so recording a default instead of
    /// the argument fails rather than passing by coincidence.
    #[test]
    fn build_records_the_inputs_it_was_given() {
        let analysed = GenomeRegions::whole_contigs(&contigs());
        let segmentation = Segmentation::build(
            std::iter::empty(),
            analysed.clone(),
            header(),
            unusual_criteria(),
            catalog_path(),
        )
        .expect("an empty stream is a segmentation over no segments");

        assert_eq!(segmentation.inputs().analysed_regions, analysed);
        assert_eq!(segmentation.inputs().catalog, header());
        assert_eq!(
            segmentation.inputs().repeat_tract_criteria,
            unusual_criteria()
        );
        assert_ne!(
            unusual_criteria(),
            StrRepeatCriteria::default(),
            "the fixture must differ from the default, or this test proves nothing",
        );
        assert!(segmentation.segments().is_empty());
    }

    /// A failure part-way through the catalog's stream is the run's failure, and **the whole
    /// chain reaches the person**: which catalog file, and what the catalog said about it.
    #[test]
    fn a_failing_catalog_stream_fails_the_build_naming_the_file() {
        let stream = vec![
            Ok(generic(1, 10)),
            Err(RepeatCatalogError::ToolVersionDiffers {
                built: "0.4.1".to_string(),
                running: "0.5.0".to_string(),
            }),
        ];

        let error = Segmentation::build(
            stream.into_iter(),
            GenomeRegions::whole_contigs(&contigs()),
            header(),
            unusual_criteria(),
            catalog_path(),
        )
        .expect_err("the stream failed, so the segmentation cannot be built");

        let rendered = format_error_chain(&error);
        assert!(
            rendered.contains("/genomes/chosen.catalog.parquet"),
            "the failure names the catalog the run read: {rendered}",
        );
        assert!(
            rendered.contains("0.4.1"),
            "the catalog's own reason survives the wrapper: {rendered}",
        );
    }

    /// The merge and the walk read their ground out of the same object.
    #[test]
    fn the_analysed_regions_come_back_as_the_merge_takes_them() {
        let segmentation = Segmentation::build(
            std::iter::empty(),
            GenomeRegions::whole_contigs(&contigs()),
            header(),
            unusual_criteria(),
            catalog_path(),
        )
        .expect("builds");

        let analysed = segmentation.analysed_regions();
        assert_eq!(analysed.len(), 2, "one span per contig: {analysed:?}");
    }

    /// The debug rendering is sizes, not the segments themselves.
    #[test]
    fn the_debug_rendering_of_a_segmentation_is_sizes() {
        let segmentation = Segmentation::build(
            vec![generic(1, 10), generic(11, 20)].into_iter().map(Ok),
            GenomeRegions::whole_contigs(&contigs()),
            header(),
            unusual_criteria(),
            catalog_path(),
        )
        .expect("builds");

        let rendered = format!("{segmentation:?}");
        assert!(rendered.contains("segments: 2"), "{rendered}");
        assert!(
            !rendered.contains("Generic"),
            "the segments themselves are not printed: {rendered}",
        );
    }

    // `SegmentationInputs::first_difference` is tested where the type now lives
    // (`crate::ng::segmentation_inputs`).
}
