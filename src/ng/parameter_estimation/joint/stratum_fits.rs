//! **The one borrow of slippage numbers that crosses the calling seam.**
//!
//! The caller scoring a repeat-tract candidate needs three numbers — how often a read reports a
//! tract length other than its allele's, which way the slips go, and how fast multi-repeat slips
//! fall off — for the read group that produced the read and the stratum the candidate sits in
//! (`doc/devel/ng/arch/read_likelihoods.md` §4.2). The pre-pass has already fitted them. What it
//! has not had is a way to ask for one, and this module is that way and nothing more.
//!
//! # What it deliberately does not do
//!
//! **It does not blend, and re-blending here would be a defect rather than a duplication.** The
//! level on a [`StratumOutcome`] is already the emitted one, by one of three routes, and its
//! [`LevelProvenance`] says which:
//!
//! - a stratum fitted on its own tracts has had its level weighed against its period's curve —
//!   [`fit_strata`](super::ssr_fit::fit_strata) draws the curves after every stratum has its own
//!   answer, then applies [`blend_level`](super::slippage_curve::blend_level) in place;
//! - a stratum too thin to be fitted takes its period's curve whole, never through `blend_level`;
//! - a run with curves switched off keeps every cell's own answer.
//!
//! Feeding any of those back through `blend_level` would weigh the curve against a number the
//! curve is already inside. `spec/str_slippage_level_curve.md` §5.1 does not name this act — what
//! it forbids in so many words is a *curve* fitted from blended values, "otherwise each round of
//! smoothing fits a curve to the previous round's curve, and the cells stop being evidence". This
//! is the same circularity one step downstream, and the measurement is in the report: re-blending
//! the five blended strata of a small real fit moves their levels by 0.6% to 4.1%, while leaving
//! the `curve_weight` in their provenance unchanged — so a consumer inspecting the provenance
//! would see nothing wrong and only the number would have moved.
//!
//! **It does not decide anything about a stratum with no answer.** Four different absences reach
//! a lookup and [`NoSlippage`] keeps them apart, because a caller that cannot tell *this run
//! never named that read group* from *that library put no read in this stratum* will report a
//! quiet tract as an unsequenced one.
//!
//! # The grain, and why the key has a read group in it
//!
//! Slippage is a property of the chemistry, so it is fitted per read group
//! (`spec/parameter_prepass_joint_fit.md` §4). A run may declare that several of its read groups
//! slip alike by naming them in one **slippage group**, and one group per read group is the
//! **specified grain** — which is not the same as what happens by default: the only builder of
//! that map in this tree pools every read group into one set unless told otherwise
//! (`examples/ng_joint_records_walk.rs`). So a lookup takes the read group the caller has and
//! this module translates, rather than making every caller carry the translation.

use std::collections::BTreeMap;

use crate::ng::types::ReadGroupId;

use super::census::Stratum;
use super::ssr_fit::{LevelProvenance, SharesProvenance, Slippage, StratumOutcome};

/// What one `(read group, stratum)` cell answers with.
///
/// **The three numbers and where they came from, together.** A level fitted from 8,000 slipped
/// reads and one read off a curve through four cells are the same `f64`, and a consumer that
/// weighs them alike is treating an interpolation as a measurement
/// (`str_slippage_level_curve.md` §8).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedSlippage {
    /// How often a read reports a tract length other than its allele's, which way the slips go,
    /// and how fast multi-repeat slips fall off — **the emitted numbers**. The level is the one
    /// the fit settled on, which the field beside it says was the cell's, the curve's, or a blend
    /// of the two; it is not always a blend, and this type never makes one.
    pub slippage: Slippage,
    /// Where the level came from: the stratum's own fit, its period's curve, or a blend of the
    /// two with the share the curve carried.
    pub level: LevelProvenance,
    /// Where the direction split and the fall-off came from. Separate from the level's, because
    /// the three numbers are smoothed on their own curves and a stratum can take its level from
    /// a curve while keeping its own shares.
    pub shares: Option<SharesProvenance>,
}

/// Why a lookup has no numbers.
///
/// **Four absences and not one**, and two of them are ordinary while two say the run is not what
/// it claims. The structural twin is [`NotIdentifiedReason`], which does the same for a
/// contamination fraction: several named reasons behind one absence, because a caller told only
/// "no number" would act on it.
///
/// [`NotIdentifiedReason`]: super::contamination::NotIdentifiedReason
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NoSlippage {
    /// **Ordinary.** No stratum of this period and repeat count is in the fit. Either the cohort
    /// holds no kept tract of that shape, or every one of them was refused. A candidate several
    /// repeats from its reference tract's length can land here on perfectly good data, so a
    /// caller has to have an answer for it.
    NoSuchStratum,
    /// **Ordinary.** The stratum is in the fit, and this read group's slippage group put no read
    /// in it. The group is named because it is what was looked under, and because two read
    /// groups pooled into one group share this answer.
    GroupPutNoReadHere { slippage_group: u32 },
    /// **The run is not what it claims.** This run's slippage fit never named this read group,
    /// so there is no slippage group to look its numbers up under — a library present at calling
    /// time that the pre-pass did not know existed.
    UnknownReadGroup,
    /// **The run is not what it claims.** The read group's slippage group is past the end of the
    /// fit's own rows, so the map this type was built with names more groups than the fit was run
    /// over. Not a quiet library: the map and the fit came from different runs.
    GroupNotInTheFit {
        slippage_group: u32,
        groups_fitted: usize,
    },
}

/// One stratum's numbers, per slippage group — the three vectors a [`StratumOutcome`] exposes,
/// held together so a lookup indexes once.
#[derive(Debug, Clone, PartialEq)]
struct StratumRow {
    slippage: Vec<Option<Slippage>>,
    level: Vec<Option<LevelProvenance>>,
    shares: Vec<Option<SharesProvenance>>,
}

/// **Every stratum's slippage numbers, indexed by the key the caller has.**
///
/// Built once per run from what [`fit_strata`](super::ssr_fit::fit_strata) returned, then read
/// unchanged for the whole of calling: slippage is a frozen parameter, and nothing about it is
/// re-estimated per locus (`arch/calling_em_loop.md` §2).
#[derive(Debug, Clone, PartialEq)]
pub struct StratumFits {
    /// Which set of slippage numbers each read group's reads are drawn under — the run's own
    /// declaration, the same map [`gather_strata`](super::ssr_fit::gather_strata) was given.
    slippage_group_of: BTreeMap<ReadGroupId, u32>,
    by_stratum: BTreeMap<Stratum, StratumRow>,
}

impl StratumFits {
    /// Gather what the fit produced.
    ///
    /// **A refused stratum contributes nothing rather than an empty row.** Its three per-group
    /// accessors return empty slices, so a row built from one would answer every read group with
    /// [`NoSlippage::GroupPutNoReadHere`] — which would say a library was silent where the truth
    /// is that the stratum has no answer for anybody. Leaving it out makes the lookup say
    /// [`NoSuchStratum`](NoSlippage::NoSuchStratum), which is the honest one of the two.
    ///
    /// **A fitted stratum and a derived one are kept alike**, exactly as
    /// [`StratumOutcome::slippage`] keeps them: the read likelihood reads three numbers and the
    /// provenance beside them, and where the fit stopped and the curve started is in the
    /// provenance rather than in which door the numbers came through.
    ///
    /// # Panics
    ///
    /// **When two outcomes name one stratum**, which would otherwise lose one of them without a
    /// word — a map insert keeps the last and says nothing, and the two levels can differ by a
    /// factor of five.
    ///
    /// **The guarantee that they do not belongs to
    /// [`gather_strata`](super::ssr_fit::gather_strata), not to `fit_strata`.** `fit_strata`
    /// returns one outcome per *evidence* it was handed, and `derive_thin_strata` rewrites those
    /// in place rather than adding any; what makes the strata distinct is that `gather_strata`
    /// keys its evidence off a map. `fit_strata` is public and three examples and a benchmark
    /// build `StratumEvidence` by hand, so a caller that assembled its own list could reach this.
    /// **Release-level, deliberately**: the cost is one comparison per stratum, of which a run has
    /// tens, and the alternative is a caller scoring every tract of one shape against another
    /// shape's polymerase.
    #[must_use]
    pub fn over(
        outcomes: &[StratumOutcome],
        slippage_group_of: BTreeMap<ReadGroupId, u32>,
    ) -> Self {
        let mut by_stratum = BTreeMap::new();
        for outcome in outcomes {
            if matches!(outcome, StratumOutcome::Refused { .. }) {
                continue;
            }
            let stratum = outcome.stratum();
            let row = StratumRow {
                slippage: outcome.slippage().to_vec(),
                level: outcome.level_provenance().to_vec(),
                shares: outcome.shares_provenance().to_vec(),
            };
            // **Checked once here rather than at every lookup.** Both of the fit's paths build
            // the three vectors from one mask in one pass, so they are the same length by
            // construction — but every field of `StratumFit` and `DerivedStratum` is public, so
            // a caller that assembled its own outcome could hand over a short one, and the
            // failure would then be an index panic inside a lookup rather than a sentence
            // naming the stratum. Build time is where a caller can act on it.
            assert!(
                row.slippage.len() == row.level.len() && row.slippage.len() == row.shares.len(),
                "the outcome for period {}, {} repeats holds {} slippage groups, {} level \
                 provenances and {} shares provenances, where the fit builds all three from one \
                 mask and they are always the same length",
                stratum.period,
                stratum.reference_repeats,
                row.slippage.len(),
                row.level.len(),
                row.shares.len(),
            );
            let displaced = by_stratum.insert(stratum, row);
            assert!(
                displaced.is_none(),
                "two of the fit's outcomes are for period {}, {} repeats, and one of them would \
                 be lost without a word — look at how the evidence handed to `fit_strata` was \
                 assembled, since `gather_strata` cannot produce a repeat",
                stratum.period,
                stratum.reference_repeats,
            );
        }
        Self {
            slippage_group_of,
            by_stratum,
        }
    }

    /// The numbers for one read group at one stratum.
    ///
    /// **Fill the stratum from the *candidate*, not from the tract.** A read's chance of slipping
    /// is a property of the tract it was copied from, and that is the candidate allele
    /// (`spec/read_likelihoods.md` §4.4): a candidate of 6 repeats and one of 12 at the same
    /// locus are drawn from different strata and slip at measurably different rates — slippage
    /// rises about 1.3-fold per repeat count over the measured range. **So the stutter parameters
    /// cannot be hoisted out of the candidate loop**, and a caller that looked one up per locus
    /// from the reference tract's own length would score every candidate there against one
    /// polymerase model.
    ///
    /// **The two numbers are taken by name rather than as a [`Stratum`], and that is the whole
    /// reason.** `Stratum`'s own field is `reference_repeats` — the right word on the fit's side
    /// of the seam, where a stratum is the bin a *reference* tract was sorted into so that
    /// tracts of one shape could be pooled. A caller handed that type has a `Stratum` for the
    /// locus already in its hand and would pass it, which is the wrong number and nothing would
    /// say so. Naming the argument `candidate_repeats` makes the mistake one somebody has to
    /// type on purpose. The bins are the same bins; only which length picks one differs.
    ///
    /// # Errors
    ///
    /// [`NoSlippage`], which names which of the four absences this is.
    pub fn at(
        &self,
        read_group: ReadGroupId,
        period: u8,
        candidate_repeats: u64,
    ) -> Result<FittedSlippage, NoSlippage> {
        let stratum = Stratum {
            period,
            reference_repeats: candidate_repeats,
        };
        let group = *self
            .slippage_group_of
            .get(&read_group)
            .ok_or(NoSlippage::UnknownReadGroup)?;
        let row = self
            .by_stratum
            .get(&stratum)
            .ok_or(NoSlippage::NoSuchStratum)?;
        let index = group as usize;
        // **A group past the end of the row is not a quiet library**, and answering as though it
        // were would hide the thing worth knowing: the map this type was built with names more
        // groups than the fit was run over, so the two were assembled from different runs. It is
        // the same class of fact as [`NoSlippage::UnknownReadGroup`] and gets its own answer.
        if index >= row.slippage.len() {
            return Err(NoSlippage::GroupNotInTheFit {
                slippage_group: group,
                groups_fitted: row.slippage.len(),
            });
        }
        // Indexed rather than `get`-ed from here on: `over` has checked that the three vectors
        // are the same length, and the bound above is that length.
        let slippage = row.slippage[index].ok_or(NoSlippage::GroupPutNoReadHere {
            slippage_group: group,
        })?;
        let level = row.level[index].expect(
            "a slippage group with numbers at a stratum has a level provenance beside them",
        );
        Ok(FittedSlippage {
            slippage,
            level,
            shares: row.shares[index],
        })
    }

    /// Which slippage group a read group's reads are drawn under, for a caller reporting what it
    /// looked up rather than looking one up.
    #[must_use]
    pub fn slippage_group_of(&self, read_group: ReadGroupId) -> Option<u32> {
        self.slippage_group_of.get(&read_group).copied()
    }

    /// How many strata carry an answer — what a run summary reports.
    ///
    /// **It cannot tell a cohort with no repeat tracts from one where every stratum was
    /// refused**, and it is not meant to: both are zero, and both mean the caller gets no
    /// slippage anywhere. A run that needs to tell them apart reads the refusals off the
    /// outcomes, which carry their own reason.
    #[must_use]
    pub fn strata(&self) -> usize {
        self.by_stratum.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::joint::share_curve::ShareSource;
    use crate::ng::parameter_estimation::joint::slippage_curve::{
        LevelSource, RiseShape, SlippageCurve,
    };
    use crate::ng::parameter_estimation::joint::ssr_fit::{
        DerivedStratum, ShareProvenance, StratumRefusal,
    };

    fn stratum(period: u8, reference_repeats: u64) -> Stratum {
        Stratum {
            period,
            reference_repeats,
        }
    }

    fn slippage(level: f64) -> Slippage {
        Slippage {
            level,
            shorter_share: 0.83,
            fall_off: 0.25,
        }
    }

    /// A curve whose value at 10 repeats is 0.11 — the shape a level provenance carries, so a
    /// test can ask the curve itself what it would have said.
    fn a_curve() -> SlippageCurve {
        SlippageCurve {
            rise_shape: RiseShape::new(1.0).expect("a rise shape of one is on the grid"),
            intercept: 0.01,
            slope: 0.01,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.077,
            cells: 5,
        }
    }

    /// A level provenance carrying a distinguishable `slipped_reads`, so that a test can tell
    /// one slippage group's provenance from another's.
    ///
    /// **The two are told apart by a number rather than by the variant**, because the variant is
    /// what a mutation would most easily preserve: reading group 0's provenance for every group
    /// leaves every `source` right and every count wrong.
    fn from_the_cell(slipped_reads: f64) -> LevelProvenance {
        LevelProvenance {
            source: LevelSource::Cell,
            curve: None,
            reach: None,
            slipped_reads: Some(slipped_reads),
        }
    }

    /// A shares provenance with real content — **the shape the fit actually emits**, where the
    /// helper below used to hand out `None`. `derive_thin_strata` and the pooled fit both set a
    /// shares provenance wherever they set a slippage number, so a fixture without one is a
    /// shape no run can produce.
    fn shares_from_a_curve(slipped_reads: f64) -> SharesProvenance {
        let from_a_curve = ShareProvenance {
            source: ShareSource::Curve,
            curve: None,
            reach: None,
        };
        SharesProvenance {
            slipped_reads: Some(slipped_reads),
            shorter_share: from_a_curve,
            fall_off: from_a_curve,
        }
    }

    /// A stratum nothing was fitted at, whose numbers came from its period's curves — the shape
    /// [`StratumOutcome::Derived`] carries, built by hand so the tests do not pay for a fit.
    ///
    /// Every group with a slippage number gets a level provenance and a shares provenance beside
    /// it, which is the invariant both of the fit's paths hold.
    fn derived(
        at: Stratum,
        slippage: Vec<Option<Slippage>>,
        level: Vec<Option<LevelProvenance>>,
    ) -> StratumOutcome {
        let shares = slippage
            .iter()
            .zip(&level)
            .map(|(numbers, level)| {
                numbers
                    .and(level.map(|level| shares_from_a_curve(level.slipped_reads.unwrap_or(1.0))))
            })
            .collect();
        StratumOutcome::Derived(Box::new(DerivedStratum {
            stratum: at,
            slippage,
            level_provenance: level,
            shares_provenance: shares,
            tracts_of_its_own: 4,
            reads_crossing: 40,
        }))
    }

    fn one_group() -> BTreeMap<ReadGroupId, u32> {
        BTreeMap::from([(ReadGroupId(0), 0)])
    }

    /// **The numbers come back keyed by the pair the caller has** — the read group that produced
    /// the read, and the candidate's own motif period and repeat count.
    #[test]
    fn a_lookup_answers_with_the_stratum_and_groups_own_numbers() {
        let fits = StratumFits::over(
            &[
                derived(
                    stratum(2, 9),
                    vec![Some(slippage(0.08)), Some(slippage(0.11))],
                    vec![Some(from_the_cell(400.0)), Some(from_the_cell(9_999.0))],
                ),
                derived(
                    stratum(2, 10),
                    vec![Some(slippage(0.10)), Some(slippage(0.13))],
                    vec![Some(from_the_cell(401.0)), Some(from_the_cell(9_998.0))],
                ),
                derived(
                    stratum(3, 9),
                    vec![Some(slippage(0.02)), Some(slippage(0.03))],
                    vec![Some(from_the_cell(402.0)), Some(from_the_cell(9_997.0))],
                ),
            ],
            BTreeMap::from([(ReadGroupId(7), 0), (ReadGroupId(9), 1)]),
        );
        let level_of = |group, period, repeats| {
            fits.at(group, period, repeats)
                .expect("the stratum has numbers")
                .slippage
                .level
        };

        // Both halves of the key are load-bearing: the same read group at two strata, and the
        // same stratum for two read groups, all four differ.
        assert_eq!(level_of(ReadGroupId(7), 2, 9), 0.08);
        assert_eq!(
            level_of(ReadGroupId(7), 2, 10),
            0.10,
            "the same library at one repeat count more",
        );
        assert_eq!(
            level_of(ReadGroupId(7), 3, 9),
            0.02,
            "the same library at the same repeat count of a longer motif",
        );
        assert_eq!(
            level_of(ReadGroupId(9), 2, 9),
            0.11,
            "the other library, which the run put in its own slippage group",
        );

        // **Each group's provenance is its own**, told apart by a count rather than by a
        // variant: reading group 0's provenance for every group would leave every `source`
        // right and every number wrong.
        let answer = fits
            .at(ReadGroupId(9), 2, 9)
            .expect("the second group has numbers");
        assert_eq!(answer.level.slipped_reads, Some(9_999.0));
        assert_eq!(
            answer
                .shares
                .expect("a group with numbers has a shares provenance")
                .slipped_reads,
            Some(9_999.0),
            "and so is its shares provenance",
        );
    }

    /// **Two candidates at one locus with different repeat counts get different numbers** — the
    /// case `spec/read_likelihoods.md` §4.4 names, and the reason the stutter parameters cannot
    /// be hoisted out of the candidate loop.
    ///
    /// A caller that keyed on the tract's *reference* length would score both against one
    /// polymerase model, and at tomato dinucleotides that is the difference between about 6 %
    /// of reads slipping and about 15 %.
    #[test]
    fn two_candidates_at_one_tract_are_two_strata() {
        let fits = StratumFits::over(
            &[
                derived(
                    stratum(2, 6),
                    vec![Some(slippage(0.06))],
                    vec![Some(from_the_cell(400.0))],
                ),
                derived(
                    stratum(2, 12),
                    vec![Some(slippage(0.15))],
                    vec![Some(from_the_cell(400.0))],
                ),
            ],
            one_group(),
        );

        let six = fits
            .at(ReadGroupId(0), 2, 6)
            .expect("six repeats is fitted");
        let twelve = fits
            .at(ReadGroupId(0), 2, 12)
            .expect("twelve repeats is fitted");
        assert_eq!(six.slippage.level, 0.06);
        assert_eq!(twelve.slippage.level, 0.15);
        assert_ne!(
            six.slippage.level, twelve.slippage.level,
            "one tract, two candidates, two slippage rates",
        );
    }

    /// **Two read groups a run declares alike share one answer**, which is what a slippage group
    /// is for: a run that knows two libraries ran on one machine may pool them, and one that
    /// pools everything is saying it cannot tell them apart.
    #[test]
    fn read_groups_pooled_into_one_slippage_group_get_one_answer() {
        let fits = StratumFits::over(
            &[derived(
                stratum(2, 9),
                vec![Some(slippage(0.08))],
                vec![Some(from_the_cell(400.0))],
            )],
            BTreeMap::from([(ReadGroupId(3), 0), (ReadGroupId(4), 0)]),
        );

        // **Unwrapped on both sides**, because two absences also compare equal: written as a
        // comparison of two `Result`s this assertion passed even under a mutation that made
        // every lookup fail.
        let one = fits.at(ReadGroupId(3), 2, 9).expect("the first library");
        let other = fits.at(ReadGroupId(4), 2, 9).expect("the second library");
        assert_eq!(one, other);
        assert_eq!(one.slippage.level, 0.08);
        assert_eq!(fits.slippage_group_of(ReadGroupId(4)), Some(0));
        assert_eq!(
            fits.slippage_group_of(ReadGroupId(5)),
            None,
            "and a read group the run never named has no group at all",
        );
    }

    /// **The four ways a lookup can come back empty are four different answers**, and a caller
    /// that could not tell them apart would read a library the fit never saw as a library that
    /// was quiet here.
    #[test]
    fn the_four_absences_are_told_apart() {
        let fits = StratumFits::over(
            &[
                derived(
                    stratum(2, 9),
                    vec![Some(slippage(0.08)), None],
                    vec![Some(from_the_cell(400.0)), None],
                ),
                StratumOutcome::Refused {
                    stratum: stratum(2, 20),
                    tracts: 3,
                    reason: StratumRefusal::BelowTheFloor {
                        tracts: 3,
                        floor: 50,
                    },
                },
            ],
            BTreeMap::from([
                (ReadGroupId(7), 0),
                (ReadGroupId(9), 1),
                (ReadGroupId(11), 4),
            ]),
        );

        assert_eq!(
            fits.at(ReadGroupId(5), 2, 9),
            Err(NoSlippage::UnknownReadGroup),
            "a read group this run's fit never named",
        );
        assert_eq!(
            fits.at(ReadGroupId(7), 2, 11),
            Err(NoSlippage::NoSuchStratum),
            "a candidate repeat count no kept tract of that period occupies",
        );
        assert_eq!(
            fits.at(ReadGroupId(9), 2, 9),
            Err(NoSlippage::GroupPutNoReadHere { slippage_group: 1 }),
            "a stratum with an answer, and a library that put no read in it",
        );
        assert_eq!(
            fits.at(ReadGroupId(11), 2, 9),
            Err(NoSlippage::GroupNotInTheFit {
                slippage_group: 4,
                groups_fitted: 2,
            }),
            "a group past the end of the fit's own rows — the map and the fit disagree, which \
             is not a quiet library",
        );
        assert_eq!(
            fits.at(ReadGroupId(7), 2, 20),
            Err(NoSlippage::NoSuchStratum),
            "a refused stratum has no answer for anybody, which is not the same claim as one \
             library being silent — so it is left out rather than carried as an empty row",
        );
        assert_eq!(
            fits.strata(),
            1,
            "the refusal is not a stratum with numbers"
        );
    }

    /// **The lookup returns the level the fit emitted, not one it recomputed from the curve.**
    ///
    /// This is the property the module exists to hold, and the fixture is built so that it can
    /// fail. The stratum's emitted level is a blend of its own fit with its period's curve, so
    /// it sits **between** the two and equals neither — `assert_ne!` pins that, because a
    /// fixture whose level happened to equal the curve's value would make every recomputation
    /// look correct and this test would prove nothing.
    ///
    /// Measured on a small real fit: putting the emitted level back through `blend_level` moves
    /// the five blended strata by 0.6 % to 4.1 %, and leaves the `curve_weight` in their
    /// provenance unchanged — so the number moves and the provenance does not say so.
    #[test]
    fn the_level_is_the_one_the_fit_emitted_and_not_the_curves_own() {
        let curve = a_curve();
        let repeats = 10;
        // A cell at 0.08 weighed against this curve's 0.11: between the two, equal to neither.
        let emitted = 0.084_711_129_860_078_91;
        assert_ne!(
            emitted,
            curve.level_at(repeats),
            "the fixture must not be a level the curve would also produce, or nothing it \
             asserts can fail",
        );
        let fits = StratumFits::over(
            &[derived(
                stratum(2, repeats),
                vec![Some(slippage(emitted))],
                vec![Some(LevelProvenance {
                    source: LevelSource::Blend {
                        curve_weight: 0.179_7,
                    },
                    curve: Some(curve),
                    reach: None,
                    slipped_reads: Some(400.0),
                })],
            )],
            one_group(),
        );

        let answer = fits
            .at(ReadGroupId(0), 2, repeats)
            .expect("the stratum has numbers");
        assert_eq!(
            answer.slippage.level, emitted,
            "the emitted level, not the curve's own value and not a second blend of the two",
        );
        assert_eq!(
            answer
                .level
                .curve
                .expect("the curve is carried")
                .level_at(repeats),
            0.11,
            "and the curve travels beside it, so a consumer can see what it would have said",
        );
        assert_eq!(
            answer.slippage.shorter_share, 0.83,
            "the two shape numbers are the ones the row carried",
        );
        assert_eq!(answer.slippage.fall_off, 0.25);
    }

    /// **Two outcomes naming one stratum is refused rather than silently halved.** A map insert
    /// keeps the last and says nothing, and the two levels can differ by a factor of five.
    #[test]
    #[should_panic(expected = "two of the fit's outcomes are for period 2, 10 repeats")]
    fn two_outcomes_for_one_stratum_are_refused() {
        let _ = StratumFits::over(
            &[
                derived(
                    stratum(2, 10),
                    vec![Some(slippage(0.05))],
                    vec![Some(from_the_cell(400.0))],
                ),
                derived(
                    stratum(2, 10),
                    vec![Some(slippage(0.99))],
                    vec![Some(from_the_cell(400.0))],
                ),
            ],
            one_group(),
        );
    }
}
