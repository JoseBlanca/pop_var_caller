//! Which tandem repeats count as STR loci — the one policy value the catalog is built
//! under and every reader asks with, plus the comparison that decides whether a file can
//! answer a reader at all.

use crate::ng::region_typing::segment_criteria::{MinCopies, SsrSegmentCriteria};
use crate::ng::tandem_repeat::PeriodRange;
use crate::ng::types::Bp;

/// The copy floors the catalog is built at, per period 1 to 6.
///
/// **Below every floor a caller routes on** — the measured stutter onsets are
/// `[8, 6, 6, 6, 5, 4]` ([`MinCopies::default`]) — by three repeats at period 1 and one at
/// period 6, so a caller can move its routing floor anywhere inside that gap by filtering
/// the file rather than re-scanning the genome. High enough that the one- to four-copy
/// tracts nothing would ever route never reach the file (spec §4.1).
pub const CATALOG_MIN_COPIES: [u32; 6] = [5, 5, 4, 4, 4, 3];

/// The copy floor for a period beyond the table — the same fallback production uses.
pub const CATALOG_MIN_COPIES_BEYOND_TABLE: u32 = 3;

/// Sequence required on each side of a tract, to the contig's own end (bp).
///
/// What makes a tract addressable at all: with less than this beside it, no read can be
/// anchored to it and no caller can use it. Below every flank a caller works with — ng's
/// short-read default is 30 — so a reader asking for more is served by arithmetic on
/// stored coordinates, and only a reader asking for less is refused (spec §1, §4.1).
pub const CATALOG_MIN_FLANK_BP: u64 = 15;

/// The length at which a tract stops being a locus and becomes a satellite (bp).
///
/// A **read-time filter over stored spans**, not a property of the file: the builder scans
/// each contig whole, so a tract of any length is recorded whole and a reader may put this
/// boundary wherever it likes (spec §4.2). 500 against step 3's calling default of 100, so
/// the boundary stays movable in both directions.
pub const CATALOG_MAX_STR_LEN_BP: u64 = 500;

/// The widest period range the scanner can serve, and what the catalog is built at:
/// homopolymers through hexamers.
pub const CATALOG_MIN_PERIOD: u8 = 1;
/// See [`CATALOG_MIN_PERIOD`].
pub const CATALOG_MAX_PERIOD: u8 = 6;

/// Which tandem repeats count as STR loci: how short they may be, how long, and how much
/// room they need beside them.
///
/// **One value, used on both sides of the file** — the builder is given one and records it
/// in the header, every reader passes one in, and the permissiveness rule is a single
/// comparison between the two ([`serves`](Self::serves), spec §4.3).
///
/// It **wraps** [`SsrSegmentCriteria`] rather than restating its five fields, so
/// classification's own admission rules keep one home and two copies of a policy cannot
/// drift apart. The two extra fields are the ones step 3 has no concept of: it drops a
/// locus only when its flank clamps to zero, and its satellite cap doubles as a window
/// margin it does not have here.
///
/// **Only three of its values are honoured by the builder** — the period range, the copy
/// floors and the minimum flank. The rest are pure filters over stored fields, so a reader
/// may change its purity floor, score floor, satellite cap or bundle radius without
/// needing a new file (spec §4.2).
#[derive(Debug, Clone, PartialEq)]
pub struct StrRepeatCriteria {
    /// Periods, copy floors, purity floor, score floor and bundle radius — step 3's own
    /// admission rules, unchanged.
    pub classification: SsrSegmentCriteria,
    /// Sequence required on each side of a tract, to the contig's end.
    ///
    /// **Not** [`SsrSegmentCriteria::bundle_threshold`], which is a bundling distance and
    /// drops a locus only when its flank clamps to nothing. This is a floor, and the file
    /// is built at [`CATALOG_MIN_FLANK_BP`].
    pub min_flank_bp: Bp,
    /// A tract longer than this is a satellite, not a locus — a read-time filter over
    /// stored spans (see [`CATALOG_MAX_STR_LEN_BP`]).
    pub max_str_len_bp: Bp,
}

impl Default for StrRepeatCriteria {
    /// **The catalog's settings, not step 3's calling defaults**: the file is built to
    /// outlive any one calling policy, so its floors sit below every one of them (spec
    /// §4.1). The purity floor, score floor and bundle radius keep classification's own
    /// defaults, since the builder applies none of them.
    fn default() -> Self {
        // `SsrSegmentCriteria::default()` is the short-read *calling* policy; the two
        // fields named here are the ones the catalog deliberately widens, and the purity
        // floor, score floor and bundle radius it inherits are never applied by the
        // builder anyway.
        let classification = SsrSegmentCriteria {
            periods: PeriodRange::new(CATALOG_MIN_PERIOD, CATALOG_MAX_PERIOD)
                .expect("the catalog period range 1..=6 is valid"),
            min_copies: MinCopies::new(CATALOG_MIN_COPIES, CATALOG_MIN_COPIES_BEYOND_TABLE),
            ..SsrSegmentCriteria::default()
        };

        Self {
            classification,
            min_flank_bp: Bp(CATALOG_MIN_FLANK_BP),
            max_str_len_bp: Bp(CATALOG_MAX_STR_LEN_BP),
        }
    }
}

impl StrRepeatCriteria {
    /// Whether a catalog built under `self` can answer a reader asking with `wanted`.
    ///
    /// `Ok(())`, or the **first** axis on which the reader is more permissive than the
    /// file, carrying both values so the error names them. Total, and does no I/O.
    ///
    /// Only the three axes the builder honours can refuse (spec §4.1): a copy floor below
    /// the built table at any period, a flank floor below the built one, or a period range
    /// reaching outside the built one. The purity floor, score floor, satellite cap and
    /// bundle radius are filters over stored fields, so **any** value of them is
    /// serviceable from any catalog and none of them appears here (spec §4.2).
    ///
    /// This is the only place the permissiveness rule is expressed; the read methods call
    /// it rather than re-implementing it.
    pub fn serves(&self, wanted: &StrRepeatCriteria) -> Result<(), CriteriaRefusal> {
        let built_periods = self.classification.periods;
        let wanted_periods = wanted.classification.periods;
        if wanted_periods.min() < built_periods.min() || wanted_periods.max() > built_periods.max()
        {
            return Err(CriteriaRefusal::PeriodRange {
                built_min: built_periods.min(),
                built_max: built_periods.max(),
                wanted_min: wanted_periods.min(),
                wanted_max: wanted_periods.max(),
            });
        }

        // Only the periods the file actually holds are worth comparing: a floor outside
        // the built range cannot be served whatever its value, and the range check above
        // already refused that case.
        for period in built_periods.min()..=built_periods.max() {
            let built = self.classification.min_copies.for_period(period);
            let asked = wanted.classification.min_copies.for_period(period);
            if asked < built {
                return Err(CriteriaRefusal::CopyFloor {
                    period,
                    built,
                    wanted: asked,
                });
            }
        }

        if wanted.min_flank_bp.get() < self.min_flank_bp.get() {
            return Err(CriteriaRefusal::MinFlank {
                built: self.min_flank_bp.get(),
                wanted: wanted.min_flank_bp.get(),
            });
        }

        Ok(())
    }
}

/// Why a reader cannot be served from a given catalog: it is asking about tracts the file
/// does not contain.
///
/// One variant per **bounded** axis, each carrying both numbers, because an error that
/// says only "not permissive enough" leaves the reader to guess which knob to move.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum CriteriaRefusal {
    /// The reader would keep tracts with fewer copies than the builder ever wrote down.
    #[error(
        "period {period}: catalog holds tracts of {built} copies and up, reader asked for {wanted}"
    )]
    CopyFloor {
        /// The period whose floor is too low.
        period: u8,
        /// The floor the catalog was built at.
        built: u32,
        /// The floor the reader asked for.
        wanted: u32,
    },

    /// The reader would keep tracts nearer a contig's end than the builder recorded.
    #[error("catalog holds tracts with {built} bp of flank and up, reader asked for {wanted} bp")]
    MinFlank {
        /// The flank floor the catalog was built at, bp.
        built: u64,
        /// The flank floor the reader asked for, bp.
        wanted: u64,
    },

    /// The reader wants periods the scan never looked for.
    #[error(
        "catalog holds periods {built_min}..={built_max}, reader asked for {wanted_min}..={wanted_max}"
    )]
    PeriodRange {
        /// Smallest period the catalog holds.
        built_min: u8,
        /// Largest period the catalog holds.
        built_max: u8,
        /// Smallest period the reader asked for.
        wanted_min: u8,
        /// Largest period the reader asked for.
        wanted_max: u8,
    },
}

/// What a caller holding step 3's own policy is asking the catalog for.
///
/// **The bridge every consumer crosses when it stops scanning.** A scan states its policy as
/// a `TypedRegionConfig`; a reader states it as a [`StrRepeatCriteria`]. The two carry the
/// same admission rules — the conversion moves them across and supplies the one field a scan
/// has no concept of.
///
/// **That field is the flank floor, and it is the catalog's own 15 bp.** A scan applies no
/// flank floor at all: it keeps any tract with a non-zero flank, which is *less* permissive
/// than nothing but *more* permissive than the file, so asking with a scan's own rule would
/// be refused. Fifteen is what the file was built at and what ng's STR locus generator
/// fetches by default, so a caller at the default is served exactly. The consequence is the
/// one stated difference between the two paths: a tract 1 to 14 bases from a contig's end is
/// a locus to a scan and is not in the file.
impl From<&crate::ng::region_typing::TypedRegionConfig> for StrRepeatCriteria {
    fn from(config: &crate::ng::region_typing::TypedRegionConfig) -> Self {
        Self {
            classification: config.criteria.clone(),
            min_flank_bp: Bp(CATALOG_MIN_FLANK_BP),
            max_str_len_bp: config.max_str_len,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The catalog's floors sit strictly below the calling floors on every period — the
    /// property the whole file rests on, and one that a later edit to either table could
    /// break silently (spec §4.1).
    #[test]
    fn the_catalog_floors_are_below_the_calling_floors_everywhere() {
        let calling = MinCopies::default();
        let catalog = MinCopies::new(CATALOG_MIN_COPIES, CATALOG_MIN_COPIES_BEYOND_TABLE);
        for period in CATALOG_MIN_PERIOD..=CATALOG_MAX_PERIOD {
            assert!(
                catalog.for_period(period) < calling.for_period(period),
                "period {period}: catalog floor {} must leave room below the calling floor {}",
                catalog.for_period(period),
                calling.for_period(period)
            );
        }
    }

    #[test]
    fn a_reader_asking_for_exactly_what_was_built_is_served() {
        let built = StrRepeatCriteria::default();
        assert_eq!(built.serves(&StrRepeatCriteria::default()), Ok(()));
    }

    /// Build criteria differing from the catalog's defaults only in `classification`.
    fn asking(classification: SsrSegmentCriteria) -> StrRepeatCriteria {
        StrRepeatCriteria {
            classification,
            ..StrRepeatCriteria::default()
        }
    }

    /// Build criteria differing from the catalog's defaults only in the copy floors.
    fn asking_floors(by_period: [u32; 6]) -> StrRepeatCriteria {
        asking(SsrSegmentCriteria {
            min_copies: MinCopies::new(by_period, CATALOG_MIN_COPIES_BEYOND_TABLE),
            ..StrRepeatCriteria::default().classification
        })
    }

    #[test]
    fn a_stricter_reader_is_served() {
        let built = StrRepeatCriteria::default();
        let wanted = StrRepeatCriteria {
            classification: SsrSegmentCriteria {
                // the calling floors, higher than the catalog's on every period
                min_copies: MinCopies::default(),
                periods: PeriodRange::new(2, 6).expect("valid"),
                ..StrRepeatCriteria::default().classification
            },
            min_flank_bp: Bp(30),
            ..StrRepeatCriteria::default()
        };
        assert_eq!(built.serves(&wanted), Ok(()));
    }

    #[test]
    fn a_lower_copy_floor_is_refused_naming_the_period_and_both_numbers() {
        let built = StrRepeatCriteria::default();
        let mut floors = CATALOG_MIN_COPIES;
        floors[2] = 3; // period 3, one below the catalog's 4
        let wanted = asking_floors(floors);

        assert_eq!(
            built.serves(&wanted),
            Err(CriteriaRefusal::CopyFloor {
                period: 3,
                built: 4,
                wanted: 3,
            })
        );
    }

    #[test]
    fn a_lower_flank_is_refused_naming_both_numbers() {
        let built = StrRepeatCriteria::default();
        let wanted = StrRepeatCriteria {
            min_flank_bp: Bp(5),
            ..StrRepeatCriteria::default()
        };

        assert_eq!(
            built.serves(&wanted),
            Err(CriteriaRefusal::MinFlank {
                built: 15,
                wanted: 5,
            })
        );
    }

    /// Step 3's own rule today keeps any tract with a non-zero flank, which is *less than
    /// 15*, so it is refused by this same test until it adopts a flank floor (spec §5.3).
    #[test]
    fn a_reader_keeping_every_non_empty_flank_is_refused() {
        let built = StrRepeatCriteria::default();
        let wanted = StrRepeatCriteria {
            min_flank_bp: Bp(1),
            ..StrRepeatCriteria::default()
        };

        assert!(matches!(
            built.serves(&wanted),
            Err(CriteriaRefusal::MinFlank { .. })
        ));
    }

    #[test]
    fn a_wider_period_range_is_refused_at_either_end() {
        let periods = |min, max| {
            asking(SsrSegmentCriteria {
                periods: PeriodRange::new(min, max).expect("valid"),
                ..StrRepeatCriteria::default().classification
            })
        };
        let built = periods(2, 5);

        assert!(matches!(
            built.serves(&periods(1, 5)),
            Err(CriteriaRefusal::PeriodRange { .. })
        ));
        assert!(matches!(
            built.serves(&periods(2, 6)),
            Err(CriteriaRefusal::PeriodRange { .. })
        ));
    }

    /// The mirror case, and it is what §4.2 turns on: the four unbounded axes are filters
    /// over stored fields, so any value of them is serviceable from any catalog.
    #[test]
    fn the_unbounded_axes_never_refuse() {
        let built = StrRepeatCriteria::default();

        let wanted = StrRepeatCriteria {
            classification: SsrSegmentCriteria {
                min_purity: 0.0, // far more permissive than the built 0.8
                min_score: i32::MIN,
                bundle_threshold: 10_000,
                ..StrRepeatCriteria::default().classification
            },
            max_str_len_bp: Bp(1_000_000),
            ..StrRepeatCriteria::default()
        };

        assert_eq!(
            built.serves(&wanted),
            Ok(()),
            "purity, score, bundle radius and the satellite cap bound nothing"
        );
    }

    /// A copy floor outside the built period range must not be compared: the range check
    /// owns that case, and comparing a fallback against a table entry would refuse for the
    /// wrong reason.
    #[test]
    fn a_floor_beyond_the_built_range_does_not_refuse_on_its_own() {
        let narrow = PeriodRange::new(1, 4).expect("valid");
        let built = asking(SsrSegmentCriteria {
            periods: narrow,
            ..StrRepeatCriteria::default().classification
        });

        // Period 5's floor is below the catalog's, but period 5 is outside the built
        // range, so it says nothing about what the file holds.
        let mut floors = CATALOG_MIN_COPIES;
        floors[4] = 1;
        let wanted = asking(SsrSegmentCriteria {
            periods: narrow,
            min_copies: MinCopies::new(floors, CATALOG_MIN_COPIES_BEYOND_TABLE),
            ..StrRepeatCriteria::default().classification
        });

        assert_eq!(built.serves(&wanted), Ok(()));
    }
}
