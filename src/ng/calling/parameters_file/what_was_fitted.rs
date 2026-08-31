//! **Which of the file's numbers rest on a measurement of this run's reads, and which do not** —
//! the question a geneticist asks first of a file they did not watch being produced
//! (`doc/devel/ng/spec/parameters_file.md` §7, §8).
//!
//! # Why the file has to answer it in one line
//!
//! **A file a run fitted and a file a run guessed are the same shape**, by design: spec §7's
//! "one writer, three sources — and after assembly the run cannot tell them apart" is what makes
//! a defaults run auditable at all. The cost is that the two look alike on the page. Every number
//! carries a `warrant`, so the answer is *in* the file — but only as a hundred-odd separate
//! answers a reader has to gather themselves, and a reader who does not know to gather them
//! reads a defaults run's file as a fit's.
//!
//! **Measured, on the two files this module's own fixtures produce.** A fitted run's carries 7 of 7
//! groups and a defaults run's carries 0 of 7 — and before this line existed, **both opened with the
//! same 39 lines of prose, and the first thing in either that said which run it was is the
//! `warrant` on line 105**, in the first row of `base_quality_calibration`. Everything above that
//! line was identical but for the cohort's own names.
//!
//! # What a *group* is, and why it is not a number
//!
//! **A group is one thing a fit either did or did not do**, which is the grain at which a reader
//! can act. There is no useful sense in which a run "fitted 4,211 of 5,102 numbers": the
//! substitution-rate table alone grows with `(read group × stratum × ploidy)`, so a count of
//! numbers would be a report about the cohort's size. A count of groups is a report about the
//! run.
//!
//! **Two of the file's sections are not groups, because no fit can produce them**, and putting
//! them in the denominator would mean no run could ever reach it:
//!
//! - `[sequencing_batches]` is *declared* by the operator (§3.4) — the run is told, or it assumes
//!   one batch; nothing measures it;
//! - `[stated_constants]` is the section for numbers nothing fits (§3.8), which is what its name
//!   says.
//!
//! `[fitted_from]` and `ploidy` are identity rather than numbers, and `format_version` is about
//! the format.

use super::{GroupOfNumbers, ParametersFile, SeedRung, Warrant, WhatTheRunFitted};

impl Warrant {
    /// **Whether a number under this warrant was measured from the reads of the run that wrote
    /// the file** — true for `fitted_here` and `borrowed`, false for `supplied` and `defaulted`.
    ///
    /// **`borrowed` counts as fitted and that is the whole of the subtlety.** A borrowed number
    /// is the mean of the sample's *other* read groups (spec §2), so it is a measurement of this
    /// cohort's reads taken at a neighbouring grain — the ladder's own stated grounds for ranking
    /// it above `supplied`. A file whose calibrations are all borrowed did fit them; a file whose
    /// calibrations are all supplied did not.
    pub(super) fn was_measured_from_the_runs_reads(self) -> bool {
        match self {
            Self::FittedHere | Self::Borrowed => true,
            Self::Supplied | Self::Defaulted => false,
        }
    }
}

impl ParametersFile {
    /// **Which groups of numbers this file's run fitted from its own reads, and which it did
    /// not.**
    ///
    /// **Derived from the file every time rather than recorded in it**, which is the rule step E3
    /// set for the note about missing slippage and the reason it holds here too: a recorded count
    /// is a second statement of something the numbers already say, and the two can come to
    /// disagree — most easily when a person edits a value and its warrant, which is exactly what
    /// the file invites them to do (§1.2 goal 3).
    ///
    /// **⚑ *Fitted* here means the file says a fit produced these numbers, not that this run's
    /// reads produced them.** Only the three groups [`GroupOfNumbers::states_whose_reads`] names
    /// can answer the second question, so a file demoted under spec §2.1 still counts its
    /// slippage, its spectra, its contamination and its prior's seed among the fitted. The file's
    /// opening line says so; the count does not silently mean more than it can.
    #[must_use]
    pub fn what_the_run_fitted(&self) -> WhatTheRunFitted {
        let mut fitted = Vec::new();
        let mut not_fitted = Vec::new();
        for group in GroupOfNumbers::EVERY {
            if self.was_fitted(group) {
                fitted.push(group);
            } else {
                not_fitted.push(group);
            }
        }
        WhatTheRunFitted { fitted, not_fitted }
    }

    /// Whether one group rests on a measurement of the run's own reads.
    ///
    /// **Three shapes of answer, because the file states the three groups' provenance three
    /// different ways**, and every one of them is the file's own statement rather than a guess
    /// from the value:
    ///
    /// - **a warrant on each row** — the calibration and the inbreeding coefficients. Fitted
    ///   where *any* row was measured, since a cohort in which one plant's coefficient could be
    ///   fitted and another's could not did fit the group.
    /// - **a warrant on each row** — and the repeat-tract substitution rates, whose rows carry one
    ///   too.
    /// - **the row's presence** — contamination, the slippage rows and the two length-spectrum
    ///   tables. A row exists only where the fit had something to say, which is spec §5's rule;
    ///   contamination needs the extra step of asking whether any row carries a `measurement`,
    ///   because a run where *some* library was measured writes rows for the ones that were not
    ///   (§3.4).
    /// - **a rung** — the ordinary-site prior, whose bottom rung is a stated species-range
    ///   heterozygosity and whose other three all rest on a fitted moment (§3.6).
    ///
    /// **⚑ Only the first shape can see spec §2.1's demotion**, and that is a limit of the
    /// format rather than of this function — see [`GroupOfNumbers::states_whose_reads`].
    fn was_fitted(&self, group: GroupOfNumbers) -> bool {
        match group {
            GroupOfNumbers::BaseQualityCalibration => self
                .base_quality_calibration
                .by_read_group
                .iter()
                .any(|row| {
                    row.error_probability_multiplier
                        .warrant
                        .was_measured_from_the_runs_reads()
                }),
            GroupOfNumbers::Contamination => self.contamination.as_ref().is_some_and(|table| {
                table
                    .by_read_group
                    .iter()
                    .any(|row| row.measurement.is_some())
            }),
            GroupOfNumbers::Inbreeding => self.inbreeding.by_sample.iter().any(|row| {
                row.inbreeding_coefficient
                    .warrant
                    .was_measured_from_the_runs_reads()
            }),
            GroupOfNumbers::OrdinarySitePrior => {
                self.ordinary_site_prior.rung != SeedRung::StatedHeterozygosity
            }
            GroupOfNumbers::RepeatTractSlippage => {
                !self.repeat_tracts.slippage_by_stratum_and_group.is_empty()
            }
            // **Either rung of the tract ladder that a fit populates**, and they are one group
            // because they are one answer to the reader's question: a run that fitted its
            // periods' curves and furnished every stratum from them fitted the length spectra.
            GroupOfNumbers::RepeatTractLengthSpectra => {
                !self.repeat_tracts.length_spectrum_by_stratum.is_empty()
                    || !self.repeat_tracts.length_spectrum_by_period.is_empty()
            }
            // **A warrant, not the row's presence** — a substitution rate is the one repeat-tract
            // number that carries one, and the demotion of spec §2.1 moves it
            // (`bindings::demoted_to_no_better_than_supplied`). Reading `is_empty()` here counted
            // a demoted file's rates as this run's own fit, which is the defect three reviewers
            // found on 2026-08-31.
            GroupOfNumbers::RepeatTractSubstitutionRates => self
                .repeat_tracts
                .substitution_rate_by_stratum
                .iter()
                .any(|row| row.rate.warrant.was_measured_from_the_runs_reads()),
        }
    }
}

impl GroupOfNumbers {
    /// Every group, in the order the file writes them.
    ///
    /// **⚑ Maintained by hand, and this is the one place a new group can be forgotten.** Adding a
    /// variant breaks the four exhaustive `match`es — `was_fitted`, `states_whose_reads`,
    /// `in_the_readers_words`, `key` — so it cannot be added silently; but nothing makes the
    /// compiler notice a variant missing from *this array*, and a group missing here drops out of
    /// the denominator with every test in the module still green, because they all iterate it.
    /// **Add the variant here first.**
    pub const EVERY: [Self; 7] = [
        Self::BaseQualityCalibration,
        Self::Contamination,
        Self::Inbreeding,
        Self::OrdinarySitePrior,
        Self::RepeatTractSlippage,
        Self::RepeatTractLengthSpectra,
        Self::RepeatTractSubstitutionRates,
    ];

    /// **What this group is, in the words a geneticist would use for it** — a noun phrase that
    /// reads inside a sentence listing several of them, and not a key path.
    ///
    /// The key path is [`Self::key`], and the two are both given because they answer different
    /// questions: this one tells a reader what was not measured, and that one tells them where to
    /// look.
    #[must_use]
    pub fn in_the_readers_words(self) -> &'static str {
        match self {
            Self::BaseQualityCalibration => "the base-quality calibration",
            Self::Contamination => "contamination",
            Self::Inbreeding => "the inbreeding coefficients",
            Self::OrdinarySitePrior => "the ordinary-site prior's seed",
            Self::RepeatTractSlippage => "repeat-tract slippage",
            Self::RepeatTractLengthSpectra => "repeat-tract length spectra",
            Self::RepeatTractSubstitutionRates => "repeat-tract substitution rates",
        }
    }

    /// **Whether this group's numbers say *whose* reads they came from.**
    ///
    /// **Three of the seven do and four do not, and that is a property of the format.** A
    /// `warrant` is the word that answers it, and only the base-quality calibration, the
    /// inbreeding coefficients and the repeat-tract substitution rates carry one on every number.
    /// A slippage row carries a *smoothing origin*, a length spectrum a *rung*, a contamination
    /// row *whose reads it was fitted from*, and the prior's seed a *rung* — each says how the
    /// number was arrived at, and none has a state meaning *somebody handed this over*.
    ///
    /// **So spec §2.1's demotion cannot reach four of the seven.**
    /// [`ParametersFile::demoted_to_no_better_than_supplied`] moves five warrants and nothing
    /// else, by design; a file demoted because it was fitted under another census therefore still
    /// shows its slippage, its spectra, its contamination and its prior's rung as fitted. The
    /// file's own opening line says which groups can answer and which cannot rather than papering
    /// over it, because the reader is who pays for the gap.
    ///
    /// **⚑ Closing it is a change to what the file records and is the owner's** — it is the same
    /// gap `demoted_to_no_better_than_supplied` records for `SeedRung::FittedCurve`, which after a
    /// demotion still reads *"both moments came off the run's own fitted population curve"*.
    #[must_use]
    pub fn states_whose_reads(self) -> bool {
        match self {
            Self::BaseQualityCalibration
            | Self::Inbreeding
            | Self::RepeatTractSubstitutionRates => true,
            Self::Contamination
            | Self::OrdinarySitePrior
            | Self::RepeatTractSlippage
            | Self::RepeatTractLengthSpectra => false,
        }
    }

    /// **Where in the file to look** — the key path, in the file's own spelling.
    #[must_use]
    pub fn key(self) -> &'static str {
        match self {
            Self::BaseQualityCalibration => "base_quality_calibration",
            Self::Contamination => "contamination",
            Self::Inbreeding => "inbreeding",
            Self::OrdinarySitePrior => "ordinary_site_prior",
            Self::RepeatTractSlippage => "repeat_tracts.slippage_by_stratum_and_group",
            // **Both tables, because the group is both** — a run whose spectra came off its
            // periods' curves has rows in the second and none in the first, so naming only the
            // first sends that reader to an empty table.
            Self::RepeatTractLengthSpectra => {
                "repeat_tracts.length_spectrum_by_stratum and .length_spectrum_by_period"
            }
            Self::RepeatTractSubstitutionRates => "repeat_tracts.substitution_rate_by_stratum",
        }
    }
}

impl WhatTheRunFitted {
    /// The groups that rest on a measurement of the run's own reads.
    #[must_use]
    pub fn fitted(&self) -> &[GroupOfNumbers] {
        &self.fitted
    }

    /// The groups that do not — every number in them a compiled-in constant or a value somebody
    /// handed the run.
    #[must_use]
    pub fn not_fitted(&self) -> &[GroupOfNumbers] {
        &self.not_fitted
    }

    /// How many groups there are in all — the denominator, and the same number for every file.
    #[must_use]
    pub fn groups(&self) -> usize {
        self.fitted.len() + self.not_fitted.len()
    }

    /// **Whether this run fitted nothing at all** — the defaults run of spec §8, and the state a
    /// reader most needs told rather than left to work out.
    #[must_use]
    pub fn nothing_was_fitted(&self) -> bool {
        self.fitted.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::super::tests::a_file_using_every_shape;
    use super::super::{GroupOfNumbers, ParametersFile, SeedRung, Warrant};

    /// Blank the one group `group` names, leaving every other group of the file untouched.
    ///
    /// **The point of the helper is that each arm is the *smallest* edit that unfits its group**,
    /// so a predicate reading the wrong table shows up as two groups moving rather than one.
    fn with_nothing_fitted_for(group: GroupOfNumbers) -> ParametersFile {
        let mut file = a_file_using_every_shape();
        match group {
            GroupOfNumbers::BaseQualityCalibration => {
                for row in &mut file.base_quality_calibration.by_read_group {
                    row.error_probability_multiplier.warrant = Warrant::Supplied;
                }
            }
            // **Every row kept and every `measurement` dropped**, which is the run in which
            // nothing could be measured anywhere — not the absent section, which is a different
            // state (spec §5's first row) and would test the `is_some` rather than the `any`.
            GroupOfNumbers::Contamination => {
                for row in &mut file
                    .contamination
                    .as_mut()
                    .expect("the fixture is contaminated")
                    .by_read_group
                {
                    row.measurement = None;
                }
            }
            GroupOfNumbers::Inbreeding => {
                for row in &mut file.inbreeding.by_sample {
                    row.inbreeding_coefficient.warrant = Warrant::Defaulted;
                    row.inbreeding_coefficient.observations = None;
                }
            }
            GroupOfNumbers::OrdinarySitePrior => {
                file.ordinary_site_prior.rung = SeedRung::StatedHeterozygosity;
            }
            GroupOfNumbers::RepeatTractSlippage => {
                file.repeat_tracts.slippage_by_stratum_and_group.clear();
            }
            GroupOfNumbers::RepeatTractLengthSpectra => {
                file.repeat_tracts.length_spectrum_by_stratum.clear();
                file.repeat_tracts.length_spectrum_by_period.clear();
            }
            GroupOfNumbers::RepeatTractSubstitutionRates => {
                file.repeat_tracts.substitution_rate_by_stratum.clear();
            }
        }
        file
    }

    #[test]
    fn a_file_a_fitted_run_wrote_says_every_group_was_fitted() {
        let what = a_file_using_every_shape().what_the_run_fitted();
        assert_eq!(what.groups(), 7);
        assert_eq!(
            what.fitted().len(),
            7,
            "not fitted: {:?}",
            what.not_fitted()
        );
        assert!(what.not_fitted().is_empty());
        assert!(!what.nothing_was_fitted());
    }

    /// **The test that tells the seven predicates apart.** Each group is unfitted on its own and
    /// nothing else may move — which is what a predicate reading a neighbouring table fails, and
    /// what a fixture that unfitted everything at once could not see.
    #[test]
    fn unfitting_one_group_moves_that_group_and_no_other() {
        for group in GroupOfNumbers::EVERY {
            let what = with_nothing_fitted_for(group).what_the_run_fitted();
            assert_eq!(
                what.not_fitted(),
                &[group],
                "unfitting {} moved something else too",
                group.key()
            );
            assert_eq!(what.fitted().len(), 6);
        }
    }

    /// **A contamination section that is absent and one in which nobody could be measured are
    /// the same answer to this question and different states of the file** (spec §5's first two
    /// rows), so both are pinned rather than only the one the helper above builds.
    #[test]
    fn an_absent_contamination_section_is_a_group_that_was_not_fitted() {
        let mut file = a_file_using_every_shape();
        file.contamination = None;
        assert_eq!(
            file.what_the_run_fitted().not_fitted(),
            &[GroupOfNumbers::Contamination]
        );
    }

    /// **A borrowed number was measured — from the sample's other read groups.** The group is
    /// fitted where every row is borrowed, which is the arm of `Warrant`'s predicate that no
    /// other test here reaches on its own: the fixture's calibration mixes one fitted row in.
    #[test]
    fn a_group_whose_every_row_is_borrowed_was_fitted() {
        let mut file = a_file_using_every_shape();
        for row in &mut file.base_quality_calibration.by_read_group {
            row.error_probability_multiplier.warrant = Warrant::Borrowed;
        }
        assert!(file.what_the_run_fitted().not_fitted().is_empty());
    }

    /// **The three rungs above the bottom one all count as fitted**, since each rests on a moment
    /// the pre-pass measured — including `ZeroDiversity`, which is a cohort measured and found
    /// invariant rather than a cohort nothing could be measured on.
    #[test]
    fn only_the_priors_bottom_rung_is_a_seed_nothing_was_fitted_for() {
        for rung in [
            SeedRung::FittedCurve,
            SeedRung::NeutralShape,
            SeedRung::ZeroDiversity,
        ] {
            let mut file = a_file_using_every_shape();
            file.ordinary_site_prior.rung = rung;
            assert!(
                file.what_the_run_fitted().not_fitted().is_empty(),
                "{rung:?} is a fitted seed"
            );
        }
    }

    /// **Either rung of the tract ladder counts**, so a run that fitted its periods' curves and
    /// furnished every stratum from them fitted the length spectra — the state spec §5's fourth
    /// row describes, where a stratum has no row of its own *because* its period had a curve.
    #[test]
    fn a_length_spectrum_furnished_from_a_period_curve_is_still_fitted() {
        let mut file = a_file_using_every_shape();
        file.repeat_tracts.length_spectrum_by_stratum.clear();
        assert!(file.what_the_run_fitted().not_fitted().is_empty());
    }

    /// **The three groups that can say whose reads, and the four that cannot** — pinned, because
    /// the file's own opening line is derived from this and a wrong answer there is a claim about
    /// somebody's data.
    ///
    /// The rule is not a list to be kept in step by hand: a group answers exactly where every one
    /// of its rows carries a `Warrant`. Checked against the shape below.
    #[test]
    fn exactly_the_groups_whose_rows_carry_a_warrant_say_whose_reads() {
        let answers: Vec<GroupOfNumbers> = GroupOfNumbers::EVERY
            .into_iter()
            .filter(|group| group.states_whose_reads())
            .collect();
        assert_eq!(
            answers,
            vec![
                GroupOfNumbers::BaseQualityCalibration,
                GroupOfNumbers::Inbreeding,
                GroupOfNumbers::RepeatTractSubstitutionRates,
            ]
        );

        // And the shape agrees: those three rows hold a `WarrantedValue` and the other four hold
        // an origin, a rung or a measurement.
        let file = a_file_using_every_shape();
        assert_eq!(
            file.base_quality_calibration.by_read_group[0]
                .error_probability_multiplier
                .warrant,
            Warrant::FittedHere
        );
        assert_eq!(
            file.inbreeding.by_sample[0].inbreeding_coefficient.warrant,
            Warrant::FittedHere
        );
        assert_eq!(
            file.repeat_tracts.substitution_rate_by_stratum[0]
                .rate
                .warrant,
            Warrant::Borrowed
        );
    }

    /// **A demoted file's substitution rates stop counting as fitted, and its slippage does not.**
    ///
    /// The half of spec §2.1's demotion this step could close and the half it could not, asserted
    /// together so that neither can move without the other being looked at. Until 2026-08-31 the
    /// substitution-rate predicate read `is_empty()`, so a demoted file counted another cohort's
    /// rates among the groups fitted from the reader's reads.
    ///
    /// **The four that stay are the format's limit, not this function's** — their rows carry no
    /// word for *handed over*, so `demoted_to_no_better_than_supplied` has nothing to move. The
    /// file says so in its own prose (`to_toml`), and closing it is the owner's.
    #[test]
    fn a_demoted_file_stops_counting_the_groups_that_carry_a_warrant() {
        let demoted = a_file_using_every_shape().demoted_to_no_better_than_supplied();
        let what = demoted.what_the_run_fitted();

        assert_eq!(
            what.not_fitted(),
            &[
                GroupOfNumbers::BaseQualityCalibration,
                GroupOfNumbers::Inbreeding,
                GroupOfNumbers::RepeatTractSubstitutionRates,
            ],
            "exactly the three that carry a warrant"
        );
        assert_eq!(
            what.fitted(),
            &[
                GroupOfNumbers::Contamination,
                GroupOfNumbers::OrdinarySitePrior,
                GroupOfNumbers::RepeatTractSlippage,
                GroupOfNumbers::RepeatTractLengthSpectra,
            ],
            "and exactly the four that carry none — the gap the file's prose discloses"
        );
    }

    /// The seven names are seven names, and the seven key paths seven key paths — a copy-paste
    /// that gave two groups one name would print a list naming one of them twice.
    #[test]
    fn every_group_is_named_once_in_each_vocabulary() {
        let mut words: Vec<_> = GroupOfNumbers::EVERY
            .iter()
            .map(|group| group.in_the_readers_words())
            .collect();
        words.sort_unstable();
        words.dedup();
        assert_eq!(words.len(), 7);

        let mut keys: Vec<_> = GroupOfNumbers::EVERY
            .iter()
            .map(|group| group.key())
            .collect();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), 7);
    }
}
