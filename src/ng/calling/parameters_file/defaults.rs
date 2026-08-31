//! **What a run puts in each of its parameters when nothing measured one** — the compiled-in
//! defaults, in one place, together with the two parameters that have no default and the different
//! reasons they have none.
//!
//! **The design is `doc/devel/ng/spec/parameters_file.md` §8**, and every bare section number
//! below is that document's unless another one is named.
//!
//! **Why the defaults are in the binary rather than in a shipped file.** A user choosing to run
//! without a fit should not have to find a file on disk first, so *run with defaults* is a flag and
//! never a path (§8, owner's decision of 2026-08-28). What a shipped file would have bought —
//! defaults a person can see and edit — is bought instead by §7: **every run writes out the
//! parameters it used**, so a defaults run still produces the file, and each number below arrives
//! in it with `warrant = "defaulted"` beside it.
//!
//! # The whole list, and it is not one list
//!
//! Seven things a run needs and no fit gave it. §8 sorts the first four by what their default *is*;
//! the last three are the ones §8's three cases have no slot for, and they are why this table has
//! seven rows rather than four.
//!
//! | what a run needs | what it takes with no fit | what that is |
//! |---|---|---|
//! | the base-quality multiplier, per read group | [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`](crate::ng::calling::likelihood::DEFAULT_ERROR_PROBABILITY_MULTIPLIER), one | the value at which the model does nothing |
//! | the repeat-tract outlier weight, one per run | [`DEFAULT_OUTLIER_WEIGHT`](crate::ng::calling::likelihood::ssr::DEFAULT_OUTLIER_WEIGHT), 0.01 | another caller's number, never measured here |
//! | the tract ladder's fallback concentration | [`STATED_FLAT_CONCENTRATION`](crate::ng::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION), one | a stated uninformative prior |
//! | contamination, per read group | **absence** — no `[contamination]` section | a model state, not a guess |
//! | the repeat-tract substitution rate, per (read group × stratum) | [`DEFAULT_SSR_SUBSTITUTION_RATE`](crate::ng::calling::inference::repeat_tract_parameters::DEFAULT_SSR_SUBSTITUTION_RATE), 0.001 | a default taken at the tract, not written in the file |
//! | the slippage numbers, per (stratum × slippage group) | **no row**, and the tract falls back to another caller's shipped model | a default that is owed a measurement |
//! | the inbreeding coefficient, per sample | **nothing** | the one parameter that may not be defaulted |
//!
//! The prior's seed is deliberately absent from the table: it has a fallback
//! ([`ExpectedHeterozygosity::SPECIES_FALLBACK`](crate::ng::types::ExpectedHeterozygosity::SPECIES_FALLBACK)),
//! but what records it is [`SeedRegime::FallbackDiversity`](crate::ng::calling::genotype_prior::SeedRegime),
//! a rung in the file's `[ordinary_site_prior]` section, and not a `warrant`. So it is marked, and
//! it is not marked the way the four above are — which is a difference §8 does not discuss.
//!
//! # What each row means for a run
//!
//! - **A multiplier of one declines to recalibrate; it does not abstain from a claim.** It leaves
//!   every read's error probability at what the instrument minted, which asserts the instrument was
//!   right — the assumption `read_likelihoods.md` §3.2 says the calibration exists to remove, and
//!   one that fitted multipliers refute routinely (a multiplier above one is common and says the
//!   instrument was optimistic). What is true is that it is the value at which the model does
//!   nothing, so it is the one default that cannot push a call in a particular direction. **It is
//!   also the only per-read-group default**, so a cohort can be half calibrated and half not.
//!
//!   **⚑ And it is the one row of this table whose `defaulted` warrant does not mean the run took
//!   the number in the middle column.** The other three warranted numbers are *taken* from a
//!   constant, so `defaulted` fixes the value; a multiplier is a fitted error **rate** divided by
//!   the geometric mean of that read group's minted error, and
//!   [`from_fitted_rate`](crate::ng::calling::likelihood::ReadGroupCalibration::from_fitted_rate)
//!   copies the **rate's** warrant onto the ratio. The pre-pass's error-rate ladder has a
//!   `Defaulted` bottom rung of its own —
//!   [`DEFAULT_ERROR_RATE`](crate::ng::parameter_estimation::generic::DEFAULT_ERROR_RATE) at
//!   0.001, taken by a read group with too few sites to fit, no sibling to borrow from and nothing
//!   supplied — so a run can legitimately write a `defaulted` multiplier of `0.001 / that
//!   library's mean minted error`, which is one only by coincidence. **So the file's reader cannot
//!   hold this key to its constant the way it holds the other two**, and `validate` says at length
//!   why not.
//!
//!   **And that is deliberate — owner's ruling, 2026-08-31.** A library's real error rate is never
//!   its reported sequencing quality: the quality scores describe base calling, and the reads also
//!   carry mismapping, chimeras and damage. So a read group the fit could not measure is charged a
//!   stated rate rather than taken at its word, and on any real library that pushes the reads the
//!   *conservative* way — on HG002's mean minted error of 2.9055 × 10⁻⁴ the multiplier is 3.44,
//!   5.4 Phred less confident than the instrument claimed. Spec §5's third row says such a read
//!   group gets "scale 1.0" and is the sentence to correct; `DEFAULT_ERROR_RATE` is itself a
//!   placeholder until it is fitted from GIAB. Recorded in `PROJECT_STATUS.md`.
//! - **The outlier weight is another caller's number**, inherited from the existing caller at 0.01
//!   and never measured here (§3.8). It is the share of a repeat tract's reads the model expects to
//!   have come from somewhere it cannot explain — a chimera, a paralogous tract, a mismapped read.
//!   **Too low and a stray read has nowhere to go but into a genotype**, so a tract with one
//!   aberrant read is called over-confidently; too high and every repeat-tract genotype loses
//!   evidence to a term that explains nothing. Nothing here measures which way 0.01 errs; what a
//!   run can look at is how its repeat-tract calls move when the number is edited, which is why
//!   §3.8 puts it in the file: *marking a number soft is the point of writing it down.*
//! - **The fallback concentration is a stated uninformative prior, not an inherited one.** It is
//!   one chromosome's worth of belief spread flat over a tract's candidate lengths — the same
//!   quantity and the same reading `ALPHA_REF` carries on the ordinary-site path — and at one
//!   chromosome the reads move the prior from the first read onward. It is reached only where a run
//!   fitted no stratum at all; a run that fitted any takes its own median (§3.7).
//! - **Contamination's default is absence, and absence is a real answer.** A run told nothing about
//!   contamination is *scored as* uncontaminated: the read likelihood has no fraction to mix in, so
//!   it computes its plain two-term formula rather than the three-term mixture. That is a modelling
//!   default and not a finding about the samples — the file says the first ("nobody identified any
//!   contamination") and never the second. There is no constant here at all: the file writes no
//!   `[contamination]` section, which is §5's first row and the absence a reader most easily
//!   collapses into a table of zeros.
//! - **The substitution rate's default is taken at the tract and never written down.** The pre-pass
//!   emits the rate as `FittedHere` or not at all, so a `(read group × stratum)` it never
//!   accumulated has no row in the file; the cell that reaches it takes
//!   [`DEFAULT_SSR_SUBSTITUTION_RATE`](crate::ng::calling::inference::repeat_tract_parameters::DEFAULT_SSR_SUBSTITUTION_RATE)
//!   and `Provenance::Defaulted`, and
//!   [`TractScoringFits`](crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits)
//!   counts how many cells did. **So `defaulted` is not a warrant this file's substitution-rate
//!   rows can legitimately carry**, and nothing checks that they do not.
//!
//! # The two parameters with no default, and the two reasons differ
//!
//! **The slippage numbers are owed a measurement; the inbreeding coefficient is forbidden a
//! default.** They look alike in the file — both are simply absent — and they are not alike.
//!
//! - **Slippage, per (stratum × slippage group).** §8's third bullet decides these are to be fitted
//!   from the GIAB HG002 alignments and compiled in like the rest, and §12 question 1 records that
//!   the measurement does not exist. So a run with no slippage fit writes no slippage rows, and
//!   **the gap is filled one level down rather than left open**:
//!   [`repeat_tract_parameters`](crate::ng::calling::inference::repeat_tract_parameters) gives such
//!   a cell [`StutterModel::hipstr_shipped`](crate::ng::alignment::StutterModel::hipstr_shipped)
//!   with `Provenance::Defaulted` and counts it. **Those are HipSTR's shipped constants and not a
//!   fit**: one read in twenty comes back a whole repeat short and one in twenty a whole repeat
//!   long — symmetric, where `StutterModel::hipstr_shipped`'s own documentation records that
//!   HipSTR's *fitted* values are contraction-biased. So a run without slippage numbers scores its
//!   tracts under a symmetric guess taken from another caller, on no organism in particular, and
//!   [`cells_with_no_fitted_slippage`](crate::ng::calling::inference::repeat_tract_parameters::TractScoringFits::cells_with_no_fitted_slippage)
//!   is what says how much of the run that was.
//! - **The inbreeding coefficient, per sample.** `parameter_estimation::generic::fallback`'s own
//!   header states the rule: *"The inbreeding coefficient has one rung and it is not a default …
//!   it is the parameter that differs most between an outcrosser and a selfing landrace, and a
//!   cohort's diversity divides by `1 − F`, so a wrong constant would be amplified rather than
//!   absorbed."* It is fitted, or it is supplied, or the run fails — and §3.5 requires at least one
//!   row. **So a run that fitted nothing and was told nothing has no coefficient to write**, and
//!   `to_toml`'s `origins::INBREEDING_COEFFICIENT` says as much in the file: *a run should not be
//!   able to write this line.* Where a defaults run gets its coefficients from is step E2's
//!   question and it is the owner's, not this module's.
//!
//! # Where each constant lives, and why not here
//!
//! **Beside the code that reads it, and named there once.** A constant re-declared here would be a
//! second spelling of a number the caller already reads, and the two could then disagree — which is
//! the failure this file exists to make visible rather than one to introduce. So this module
//! documents the set and pins its behaviour; the numbers stay with their readers:
//! [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`](crate::ng::calling::likelihood::DEFAULT_ERROR_PROBABILITY_MULTIPLIER)
//! beside [`ReadGroupCalibration`](crate::ng::calling::likelihood::ReadGroupCalibration), whose
//! `defaulted` constructor is the only thing that *reads* it — the projection in from a file builds
//! the struct literally, which is the path the reader's check below guards;
//! [`DEFAULT_OUTLIER_WEIGHT`](crate::ng::calling::likelihood::ssr::DEFAULT_OUTLIER_WEIGHT) beside
//! [`RepeatTractOutlierWeight`](crate::ng::calling::likelihood::ssr::RepeatTractOutlierWeight);
//! [`STATED_FLAT_CONCENTRATION`](crate::ng::parameter_estimation::joint::stratum_fits::STATED_FLAT_CONCENTRATION)
//! beside [`StratumFits`](crate::ng::parameter_estimation::joint::stratum_fits::StratumFits);
//! [`DEFAULT_SSR_SUBSTITUTION_RATE`](crate::ng::calling::inference::repeat_tract_parameters::DEFAULT_SSR_SUBSTITUTION_RATE)
//! beside the assembly that takes it.
//!
//! **The sentence a reader is shown beside each defaulted number is `to_toml`'s `origins`**, which
//! is the list this module reconciles against — one place, so a number's origin and the comment
//! above it in the file cannot drift apart.
//!
//! # What "the warrant says `defaulted`" is worth
//!
//! **A defaulted value and a fitted one are the same number.** A multiplier of one that nobody
//! fitted and a multiplier of one a fit arrived at multiply every read's error probability
//! identically, so nothing downstream of the arithmetic can tell them apart — the warrant beside
//! the value is the whole of the difference, and §5's third row is the requirement that it survive
//! a write and a read. The tests below hold each of the three constants to its warrant at the point
//! of use, and hold the file's reader to refusing a `defaulted` value that is not the number this
//! caller holds.

#[cfg(test)]
mod tests {
    use crate::ng::calling::likelihood::ssr::{DEFAULT_OUTLIER_WEIGHT, RepeatTractOutlierWeight};
    use crate::ng::calling::likelihood::{
        DEFAULT_ERROR_PROBABILITY_MULTIPLIER, MIN_BASE_ERROR, ReadGroupCalibration,
    };
    use crate::ng::calling::parameters_file::Warrant;
    use crate::ng::calling::parameters_file::tests::a_file_using_every_shape;
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::joint::stratum_fits::{
        LengthSpectrumRung, STATED_FLAT_CONCENTRATION, StratumFits,
    };
    use std::collections::BTreeMap;

    /// **The multiplier a read group with no fitted rate is charged under leaves its reads
    /// exactly as the instrument minted them**, and says it was defaulted.
    ///
    /// The second half is what the assertion on the arithmetic is for: `charged_error` is
    /// `scale · exp(q_sum / n)`, so at a scale of one it is the geometric mean of the reads' own
    /// error probabilities and nothing else. Compared on `to_bits` as the module's other float
    /// assertions are — it is the same verdict as `==` on an ordinary finite double, and the
    /// habit is what step C3 established so that a lost sign on a zero cannot pass.
    #[test]
    fn a_read_group_with_no_fitted_rate_is_charged_what_its_reads_were_minted_with() {
        let defaulted = ReadGroupCalibration::defaulted();
        assert_eq!(defaulted.scale, DEFAULT_ERROR_PROBABILITY_MULTIPLIER);
        assert_eq!(defaulted.provenance, Provenance::Defaulted);

        // Three reads spanning Phred 40 to Phred 13, as the calibration's own fixture does.
        let minted = [1e-4_f64, 10f64.powf(-2.0), 10f64.powf(-1.3)];
        let q_sum: f64 = minted.iter().map(|error| error.ln()).sum();
        let reads = u32::try_from(minted.len()).expect("three reads fit in a u32");
        let geometric_mean = (q_sum / f64::from(reads)).exp();
        assert!(
            geometric_mean > MIN_BASE_ERROR,
            "the fixture must sit above the floor, or this asserts the floor instead"
        );
        assert_eq!(
            defaulted.charged_error(q_sum, reads).to_bits(),
            geometric_mean.to_bits(),
            "a defaulted calibration changes no read's error probability"
        );
        // `ln 1 = 0` exactly, so the log-space form adds nothing either.
        assert_eq!(defaulted.log_scale(), 0.0);
    }

    /// **The outlier weight a run inherits is the existing caller's 0.01 and says it was
    /// defaulted** — the one number here that is a guess at a quantity nothing in this project has
    /// measured.
    #[test]
    fn the_outlier_weight_a_run_inherits_is_the_existing_callers_number() {
        let inherited = RepeatTractOutlierWeight::defaulted();
        assert_eq!(inherited.value(), DEFAULT_OUTLIER_WEIGHT);
        assert_eq!(inherited.provenance(), Provenance::Defaulted);
        assert_eq!(DEFAULT_OUTLIER_WEIGHT, 0.01);
    }

    /// **A run that fitted no stratum states the flat concentration, says it was defaulted, and
    /// seeds every tract from the ladder's bottom rung.**
    ///
    /// All three matter together: the number alone cannot say where it came from — a run whose
    /// own strata's median happened to be 1.0 would carry the same number — and the rung is what
    /// a tract's own record shows.
    #[test]
    fn a_run_that_fitted_no_stratum_seeds_every_tract_from_the_flat_rung() {
        let nothing_fitted = StratumFits::over(&[], BTreeMap::new());
        assert_eq!(
            nothing_fitted.stated_concentration(),
            STATED_FLAT_CONCENTRATION
        );
        assert_eq!(
            nothing_fitted.stated_concentration_warrant(),
            Provenance::Defaulted
        );
        assert_eq!(STATED_FLAT_CONCENTRATION, 1.0);

        let spectrum = nothing_fitted.length_spectrum_at(2, 11);
        assert_eq!(spectrum.rung(), LengthSpectrumRung::StatedFlat);
        assert_eq!(spectrum.concentration(), STATED_FLAT_CONCENTRATION);
        assert!(
            spectrum.fitted_weights().is_none(),
            "the bottom rung has no shape of its own to hand out"
        );
    }

    /// **The two parameters with no default are absent from the file a run with nothing fitted
    /// writes, and the file's own reader accepts that.**
    ///
    /// **Contamination's absence is already held elsewhere** — `to_run_parameters`'s
    /// `an_absent_contamination_table_is_not_a_table_of_zeros` is spec §5's first-row fixture and
    /// covers what the run then scores and reports, and `to_toml` pins that no section is
    /// written. What no test held is the pair above: that a file naming **no slippage row and no
    /// substitution-rate row at all** is a legal file rather than a broken one, which is what a
    /// defaults run will produce and what step E2 depends on.
    #[test]
    fn a_file_with_nothing_fitted_for_repeat_tracts_is_still_a_legal_file() {
        let mut file = a_file_using_every_shape();
        file.contamination = None;
        file.repeat_tracts.slippage_by_stratum_and_group.clear();
        file.repeat_tracts.slippage_group_by_read_group.clear();
        file.repeat_tracts.substitution_rate_by_stratum.clear();
        file.repeat_tracts.length_spectrum_by_stratum.clear();
        file.repeat_tracts.length_spectrum_by_period.clear();
        // The bottom rung is then the only thing left saying anything about a tract, and its
        // warrant has to be the one a run that fitted nothing carries.
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .warrant = Warrant::Defaulted;
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .value = STATED_FLAT_CONCENTRATION;
        file.repeat_tracts
            .fallback_length_spectrum_concentration
            .observations = None;

        let projected = file
            .to_run_parameters()
            .expect("a run that fitted no repeat tracts writes a file its own reader accepts");
        let fits = projected.parameters.ssr_slippage_fits();
        assert_eq!(fits.strata(), 0);
        assert_eq!(fits.stated_concentration(), STATED_FLAT_CONCENTRATION);
        assert_eq!(fits.stated_concentration_warrant(), Provenance::Defaulted);
        assert_eq!(projected.parameters.ssr_substitution_rate().count(), 0);
        assert!(projected.parameters.view().contamination_is_absent());
    }

    /// **The file's reader holds the two keys it can to their own constant**: a value marked
    /// `defaulted` that is not the number this caller holds is refused, naming the key, quoting
    /// both numbers as the file spells them, and saying what to type instead.
    ///
    /// This is the edit spec §7's third bullet invites — copy the file your run wrote and change
    /// one line — made to a number whose warrant the reader forgot to move. **Both refusals share
    /// one closing clause**, which is what the last assertion holds: a reader who has met one of
    /// them can act on the other without re-reading.
    ///
    /// **The base-quality multiplier is deliberately not here**, and the module header says why:
    /// its `defaulted` warrant is copied from the error rate the multiplier was built from, so a
    /// legitimate run writes one at a value that is not
    /// [`DEFAULT_ERROR_PROBABILITY_MULTIPLIER`](crate::ng::calling::likelihood::DEFAULT_ERROR_PROBABILITY_MULTIPLIER).
    /// The test below is the one that pins it accepted.
    ///
    /// **The concentration's edited value is a whole number on purpose.** `Display` and `Debug`
    /// agree on 3.5 and differ on 3.0, so a fixture at 3.5 leaves the `{:?}` this message is
    /// written with unpinned — E1's review measured the revert to `{}` surviving its whole
    /// mutation suite. The outlier weight has no whole number to use, since a legal weight is
    /// strictly inside zero and one, so its `{:?}` is convention rather than something a test
    /// can hold.
    #[test]
    fn a_defaulted_value_that_is_not_the_binarys_own_number_is_refused() {
        let weight_edited = {
            let mut file = a_file_using_every_shape();
            file.stated_constants.repeat_tract_outlier_weight.warrant = Warrant::Defaulted;
            file.stated_constants.repeat_tract_outlier_weight.value = DEFAULT_OUTLIER_WEIGHT * 2.0;
            file.stated_constants
                .repeat_tract_outlier_weight
                .observations = None;
            file
        };
        let concentration_edited = {
            let mut file = a_file_using_every_shape();
            let rung = &mut file.repeat_tracts.fallback_length_spectrum_concentration;
            rung.warrant = Warrant::Defaulted;
            rung.value = STATED_FLAT_CONCENTRATION + 2.0;
            rung.observations = None;
            file
        };

        for (edited, key, constant, edited_to) in [
            (
                weight_edited,
                "stated_constants.repeat_tract_outlier_weight",
                DEFAULT_OUTLIER_WEIGHT,
                DEFAULT_OUTLIER_WEIGHT * 2.0,
            ),
            (
                concentration_edited,
                "repeat_tracts.fallback_length_spectrum_concentration",
                STATED_FLAT_CONCENTRATION,
                STATED_FLAT_CONCENTRATION + 2.0,
            ),
        ] {
            let refusal = edited
                .validate()
                .expect_err("a `defaulted` value that is not the constant is refused")
                .to_string();
            assert!(
                refusal.contains(key),
                "the refusal must name the key to edit; got {refusal}"
            );
            // **Both numbers spelled as the file spells them**, which is `Debug` for a float —
            // the writer formats every value with it, so `3.0` is the string a reader can search
            // their own file for and `3` is not.
            assert!(
                refusal.contains(&format!("{constant:?}")),
                "the refusal must quote the number this caller holds, as the file spells it; \
                 got {refusal}"
            );
            assert!(
                refusal.contains(&format!("{edited_to:?}")),
                "the refusal must quote the number in the file, as the file spells it; \
                 got {refusal}"
            );
            // **One closing clause, word for word, on both.** Two sentences that each *mention*
            // `supplied` are two messages a reader has to parse separately; a reader who has
            // acted on one of these has acted on the other. The clause is the outlier weight's,
            // which was the first written and is the one that names the claim and the fix.
            assert!(
                refusal.ends_with(
                    "a number you changed is one the run was handed, so change the warrant \
                     beside it to `supplied`"
                ),
                "the two refusals must close with one clause, so a reader meets one shape of \
                 message rather than two; got {refusal}"
            );
        }
    }

    /// **Each of the three constants is accepted beside a `defaulted` warrant**, so the test
    /// above refuses the edit rather than the shape it was made in — and the multiplier, which
    /// nothing checks, is accepted at its constant too.
    #[test]
    fn each_of_the_three_constants_is_accepted_beside_a_defaulted_warrant() {
        let mut file = a_file_using_every_shape();

        let row = &mut file.base_quality_calibration.by_read_group[0];
        row.error_probability_multiplier.warrant = Warrant::Defaulted;
        row.error_probability_multiplier.value = DEFAULT_ERROR_PROBABILITY_MULTIPLIER;
        row.error_probability_multiplier.observations = None;

        file.stated_constants.repeat_tract_outlier_weight.warrant = Warrant::Defaulted;
        file.stated_constants.repeat_tract_outlier_weight.value = DEFAULT_OUTLIER_WEIGHT;
        file.stated_constants
            .repeat_tract_outlier_weight
            .observations = None;

        let rung = &mut file.repeat_tracts.fallback_length_spectrum_concentration;
        rung.warrant = Warrant::Defaulted;
        rung.value = STATED_FLAT_CONCENTRATION;
        rung.observations = None;

        file.validate()
            .expect("the three constants are what a `defaulted` warrant claims");
    }
}
