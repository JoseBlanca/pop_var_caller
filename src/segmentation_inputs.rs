//! The values a segmentation is a function of — the operand of every compatibility check.
//!
//! A run's segments are computed from the reference, the setting for repeat ground, the repeat
//! catalog and the criteria the run asked it with, and the regions it was asked to analyse
//! (`doc/devel/ng/spec/run_streaming.md` §4.2, `command_line.md` §4.4.1). This record is those
//! inputs, kept so two segmentations can be compared without recomputing either.
//!
//! **A top-level module because both sides of the psp seam hold it.** The run side records it
//! when a segmentation is built ([`crate::run::Segmentation`]) and compares it when a
//! cohort of stored files is opened; the psp header carries it in every file
//! ([`crate::psp`]). Living under either of those would make the other reach across a
//! pipeline-stage boundary for its interchange type.

use crate::region_typing::GenomeRegions;
use crate::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};

/// The values a segmentation is a function of.
///
/// Two different questions are asked of this record, and each mode asks a different one
/// (spec `run_streaming.md` §6.2): whether a stored file was written over the same ground the
/// run analyses, and whether it was written under the same setting, catalog and repeat-tract
/// criteria. Every refusal names the field that differs, because "these two disagree" leaves a
/// user nothing to fix.
#[derive(Debug, Clone, PartialEq)]
pub struct SegmentationInputs {
    /// Which path each kind of repeat ground went to, and what the repeats were found with.
    pub repeat_ground: RepeatGround,
    /// The regions the run was asked to analyse — the field a user actually changes between
    /// runs, and the one compared across a cohort of stored files.
    pub analysed_regions: GenomeRegions,
}

/// Which path a run sent each kind of repeat ground to — one of the three settings of
/// `command_line.md` §4.4.1 — with the catalog and criteria that found the repeats, under the two
/// settings that read a catalog.
///
/// | Setting | STR tracts and bundles | Satellites |
/// |---|---|---|
/// | [`StrPaths`](Self::StrPaths), the default | the STR paths | not called |
/// | [`CallSatellites`](Self::CallSatellites) | the STR paths | the generic path |
/// | [`GenericPathOnly`](Self::GenericPathOnly) | the generic path | the generic path |
///
/// A *satellite* is a tract longer than `--max-str-len` or a bundle longer than
/// `--max-str-bundle-len`. **An enum rather than a flag beside the catalog**, so a run that read
/// no catalog cannot record one.
#[derive(Debug, Clone, PartialEq)]
pub enum RepeatGround {
    /// Tracts to the STR path, bundles to the STR-bundle path, satellites not called.
    StrPaths(RepeatSorting),
    /// As [`StrPaths`](Self::StrPaths), and satellites to the generic path.
    CallSatellites(RepeatSorting),
    /// Every repeat to the generic path. No catalog was read, so none is recorded.
    GenericPathOnly,
}

/// What a run sorted repeat ground with: the catalog it read and the criteria it asked it with.
#[derive(Debug, Clone, PartialEq)]
pub struct RepeatSorting {
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
}

impl RepeatGround {
    /// The catalog and criteria the repeats were found with, or `None` under
    /// [`GenericPathOnly`](Self::GenericPathOnly), which read no catalog.
    #[must_use]
    pub fn sorting(&self) -> Option<&RepeatSorting> {
        match self {
            Self::StrPaths(sorting) | Self::CallSatellites(sorting) => Some(sorting),
            Self::GenericPathOnly => None,
        }
    }

    /// The setting's name as the command line and the psp header spell it: `str-paths`,
    /// `call-satellites` or `generic-path-only`.
    #[must_use]
    pub fn setting_name(&self) -> &'static str {
        match self {
            Self::StrPaths(_) => Self::STR_PATHS,
            Self::CallSatellites(_) => Self::CALL_SATELLITES,
            Self::GenericPathOnly => Self::GENERIC_PATH_ONLY,
        }
    }

    /// The default setting's name, [`StrPaths`](Self::StrPaths).
    pub const STR_PATHS: &'static str = "str-paths";
    /// [`CallSatellites`](Self::CallSatellites)'s name, the option `--call-satellites`.
    pub const CALL_SATELLITES: &'static str = "call-satellites";
    /// [`GenericPathOnly`](Self::GenericPathOnly)'s name, the option `--generic-path-only`.
    pub const GENERIC_PATH_ONLY: &'static str = "generic-path-only";
}

impl SegmentationInputs {
    /// The criteria the repeats were found with, or `None` when the run sorted none
    /// ([`RepeatGround::GenericPathOnly`]).
    #[must_use]
    pub fn repeat_tract_criteria(&self) -> Option<&StrRepeatCriteria> {
        self.repeat_ground
            .sorting()
            .map(|sorting| &sorting.repeat_tract_criteria)
    }

    /// The catalog the repeats were found in, or `None` when the run read none
    /// ([`RepeatGround::GenericPathOnly`]).
    #[must_use]
    pub fn catalog(&self) -> Option<&RepeatCatalogHeader> {
        self.repeat_ground.sorting().map(|sorting| &sorting.catalog)
    }

    /// What [`first_difference`](Self::first_difference) calls the regions a run was asked to
    /// analyse.
    ///
    /// **A constant because one caller acts on which field it was.** A cohort whose files
    /// disagree about the ground is refused in its own words — the samples were walked over
    /// different ground, so they cannot be called as one
    /// ([`RunError::AnalysedRegionsDiffer`](crate::run::RunError::AnalysedRegionsDiffer)) —
    /// and a cohort disagreeing about the catalog or the criteria is refused in others. Matching
    /// on the string spelled out at the call site would put the two names a file apart, where
    /// changing one of them is a silent change of behaviour rather than a compile error.
    pub const ANALYSED_REGIONS: &'static str = "set of analysed regions";

    /// The name of the first field that differs, or `None` when the two agree.
    ///
    /// **A name rather than a `bool`**: a refusal that says only "these two segmentations
    /// differ" leaves the user nothing to act on (spec §6.1). The names are written to be read
    /// inside a sentence — arch §5's refusal renders them as "written under a different
    /// {field}" — so they are noun phrases in the user's vocabulary, not field identifiers.
    ///
    /// **The order is the order a person should fix them in**, and it is deliberate. The setting
    /// comes first: it decides whether there is a catalog to compare at all, and a psp written
    /// with `--call-satellites` beside one written without it holds different loci whatever
    /// else agrees. The catalog next, because it carries the reference's identity, so a catalog
    /// difference makes the rest meaningless. Criteria next, because they decide where segments
    /// end; the analysed regions last, because they are the one a user changes on purpose.
    #[must_use]
    pub fn first_difference(&self, other: &Self) -> Option<&'static str> {
        if self.repeat_ground.setting_name() != other.repeat_ground.setting_name() {
            return Some("setting for repeat ground");
        }
        if let (Some(ours), Some(theirs)) =
            (self.repeat_ground.sorting(), other.repeat_ground.sorting())
        {
            if ours.catalog != theirs.catalog {
                return Some("repeat catalog");
            }
            if ours.repeat_tract_criteria != theirs.repeat_tract_criteria {
                return Some("set of repeat-tract criteria");
            }
        }
        if self.analysed_regions != other.analysed_regions {
            return Some(Self::ANALYSED_REGIONS);
        }
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::regions::ContigBounds;
    use crate::tandem_repeat::ScanParams;

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

    fn sorting() -> RepeatSorting {
        RepeatSorting {
            catalog: header(),
            repeat_tract_criteria: StrRepeatCriteria::default(),
        }
    }

    fn inputs() -> SegmentationInputs {
        SegmentationInputs {
            repeat_ground: RepeatGround::StrPaths(sorting()),
            analysed_regions: GenomeRegions::whole_contigs(&contigs()),
        }
    }

    /// The sorting of a fixture that has one, for the fixtures that vary it.
    fn sorting_of(inputs: &mut SegmentationInputs) -> &mut RepeatSorting {
        match &mut inputs.repeat_ground {
            RepeatGround::StrPaths(sorting) | RepeatGround::CallSatellites(sorting) => sorting,
            RepeatGround::GenericPathOnly => panic!("the fixture sorts its repeats with a catalog"),
        }
    }

    fn with_other_catalog(mut inputs: SegmentationInputs) -> SegmentationInputs {
        sorting_of(&mut inputs).catalog.reference_md5 = [9; 16];
        inputs
    }

    fn with_other_criteria(mut inputs: SegmentationInputs) -> SegmentationInputs {
        sorting_of(&mut inputs).repeat_tract_criteria = unusual_criteria();
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

    /// Two psps typed with different bundle caps are two typings, and are refused together
    /// as any two are (`str_bundle_loci.md` §5.2): a bundle one of them calls the other
    /// typed as a satellite.
    #[test]
    fn a_different_bundle_cap_is_a_different_set_of_criteria() {
        let mut capped_lower = inputs();
        sorting_of(&mut capped_lower)
            .repeat_tract_criteria
            .max_str_bundle_len_bp = crate::types::Bp(60);
        assert_eq!(
            inputs().first_difference(&capped_lower),
            Some("set of repeat-tract criteria"),
        );
    }

    /// **Through the constant, which is what a caller matches on.** A second spelling of the
    /// name here would let the two drift apart, and a cohort opener that routes the analysed
    /// regions to their own refusal would then route them to the general one instead — with
    /// every test still green.
    #[test]
    fn different_analysed_regions_are_named() {
        assert_eq!(
            inputs().first_difference(&with_other_regions(inputs())),
            Some(SegmentationInputs::ANALYSED_REGIONS),
        );
        assert_eq!(
            SegmentationInputs::ANALYSED_REGIONS,
            "set of analysed regions",
            "the words a refusal is read in",
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

    /// **A psp written with `--call-satellites` and one written without it disagree on the setting,
    /// and the setting is named before anything else** — even with the catalog and regions also
    /// differing, since under different settings the two hold different loci whatever else agrees
    /// (`command_line.md` §4.4.1: a cohort mixing settings is refused).
    #[test]
    fn a_different_setting_for_repeat_ground_is_named_first() {
        let mut other = with_other_regions(with_other_catalog(inputs()));
        let RepeatGround::StrPaths(sorting) = other.repeat_ground else {
            panic!("the fixture is the default setting");
        };
        other.repeat_ground = RepeatGround::CallSatellites(sorting);

        assert_eq!(
            inputs().first_difference(&other),
            Some("setting for repeat ground")
        );
    }

    /// Two `--generic-path-only` segmentations have no catalog to compare, so only their regions
    /// can differ — and against a setting that read a catalog, the setting is what differs.
    #[test]
    fn the_generic_path_only_setting_compares_its_regions_alone() {
        let generic_only = |inputs: SegmentationInputs| SegmentationInputs {
            repeat_ground: RepeatGround::GenericPathOnly,
            ..inputs
        };

        assert_eq!(
            generic_only(inputs()).first_difference(&generic_only(inputs())),
            None
        );
        assert_eq!(
            generic_only(inputs()).first_difference(&generic_only(with_other_regions(inputs()))),
            Some(SegmentationInputs::ANALYSED_REGIONS)
        );
        assert_eq!(
            inputs().first_difference(&generic_only(inputs())),
            Some("setting for repeat ground")
        );
    }
}
