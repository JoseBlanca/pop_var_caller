//! The values a segmentation is a function of — the operand of every compatibility check.
//!
//! A run's segments are computed from the reference, the repeat catalog, the criteria the run
//! asked with and the regions it was asked to analyse
//! (`doc/devel/ng/spec/run_streaming.md` §4.2). This record is those inputs, kept so two
//! segmentations can be compared without recomputing either.
//!
//! **A top-level module because both sides of the psp seam hold it.** The run side records it
//! when a segmentation is built ([`crate::ng::run::Segmentation`]) and compares it when a
//! cohort of stored files is opened; the psp header carries it in every file
//! ([`crate::ng::psp`]). Living under either of those would make the other reach across a
//! pipeline-stage boundary for its interchange type.

use crate::ng::region_typing::GenomeRegions;
use crate::ng::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};

/// The values a segmentation is a function of.
///
/// Two different questions are asked of this record, and each mode asks a different one
/// (spec `run_streaming.md` §6.2): whether a stored file was written over the same ground the
/// run analyses, and whether it was written under the same catalog and the same repeat-tract
/// criteria. Both refusals name the field that differs, because "these two disagree" leaves a
/// user nothing to fix.
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentationInputs {
    /// The catalog file's own header, **reused whole rather than restated**: it already
    /// carries the whole-reference digest, the criteria the catalog was built under, the
    /// scan weights and the tool version.
    pub catalog: RepeatCatalogHeader,
    /// The criteria the *reader* asked with, which is what decides where a segment ends: which
    /// repeat periods count, how pure a tract must be, how close two tracts must be to be
    /// bundled.
    ///
    /// Not the same value as `catalog.built_under`: the catalog is built below every floor a
    /// reader might ask with, so that a reader filters rather than re-scans.
    pub repeat_tract_criteria: StrRepeatCriteria,
    /// The regions the run was asked to analyse — the field a user actually changes between
    /// runs, and the one compared across a cohort of stored files.
    pub analysed_regions: GenomeRegions,
}

impl SegmentationInputs {
    /// The name of the first field that differs, or `None` when the two agree.
    ///
    /// **A name rather than a `bool`**: a refusal that says only "these two segmentations
    /// differ" leaves the user nothing to act on (spec §6.1). The names are written to be read
    /// inside a sentence — arch §5's refusal renders them as "written under a different
    /// {field}" — so they are noun phrases in the user's vocabulary, not field identifiers.
    ///
    /// **The order is the order a person should fix them in**, and it is deliberate: the
    /// catalog carries the reference's identity, so a catalog difference makes the other two
    /// comparisons meaningless. Criteria come next because they decide where segments end;
    /// the analysed regions last, because they are the one a user changes on purpose.
    #[must_use]
    pub fn first_difference(&self, other: &Self) -> Option<&'static str> {
        if self.catalog != other.catalog {
            return Some("repeat catalog");
        }
        if self.repeat_tract_criteria != other.repeat_tract_criteria {
            return Some("set of repeat-tract criteria");
        }
        if self.analysed_regions != other.analysed_regions {
            return Some("set of analysed regions");
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::tandem_repeat::ScanParams;
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

    /// **Criteria that are not the default**, so a fixture varying them cannot agree with one
    /// that kept them.
    fn unusual_criteria() -> StrRepeatCriteria {
        let mut criteria = StrRepeatCriteria::default();
        criteria.classification.min_purity = 0.93;
        criteria
    }

    fn inputs() -> SegmentationInputs {
        SegmentationInputs {
            catalog: header(),
            repeat_tract_criteria: StrRepeatCriteria::default(),
            analysed_regions: GenomeRegions::whole_contigs(&contigs()),
        }
    }

    fn with_other_catalog(mut inputs: SegmentationInputs) -> SegmentationInputs {
        inputs.catalog.reference_md5 = [9; 16];
        inputs
    }

    fn with_other_criteria(mut inputs: SegmentationInputs) -> SegmentationInputs {
        inputs.repeat_tract_criteria = unusual_criteria();
        inputs
    }

    fn with_other_regions(mut inputs: SegmentationInputs) -> SegmentationInputs {
        inputs.analysed_regions = GenomeRegions::whole_contigs(&contigs()[..1]);
        inputs
    }

    #[test]
    fn identical_inputs_have_no_first_difference() {
        assert_eq!(inputs().first_difference(&inputs()), None);
    }

    #[test]
    fn a_different_catalog_is_named() {
        assert_eq!(
            inputs().first_difference(&with_other_catalog(inputs())),
            Some("repeat catalog"),
        );
    }

    #[test]
    fn different_repeat_tract_criteria_are_named() {
        assert_eq!(
            inputs().first_difference(&with_other_criteria(inputs())),
            Some("set of repeat-tract criteria"),
        );
    }

    #[test]
    fn different_analysed_regions_are_named() {
        assert_eq!(
            inputs().first_difference(&with_other_regions(inputs())),
            Some("set of analysed regions"),
        );
    }

    /// **When two inputs differ at once, the catalog is what the person is told about.**
    ///
    /// Only a fixture differing in two fields can see the order at all, and the order is
    /// load-bearing: the catalog carries the reference's identity, so under a different
    /// catalog the other two comparisons are about different genomes. A run told "a different
    /// set of analysed regions" would be sent to change `--regions`, which cannot help.
    #[test]
    fn the_catalog_is_named_before_the_criteria_and_the_regions() {
        let other = with_other_regions(with_other_criteria(with_other_catalog(inputs())));

        assert_eq!(inputs().first_difference(&other), Some("repeat catalog"));
    }

    /// And with the catalog agreeing, the criteria are named before the regions — they decide
    /// where a segment ends, where the regions only decide which ground is looked at.
    #[test]
    fn the_criteria_are_named_before_the_regions() {
        let other = with_other_regions(with_other_criteria(inputs()));

        assert_eq!(
            inputs().first_difference(&other),
            Some("set of repeat-tract criteria"),
        );
    }
}
