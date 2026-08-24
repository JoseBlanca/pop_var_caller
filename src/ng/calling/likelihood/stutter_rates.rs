//! The seven shares of the stutter distribution, from the three numbers the fit reports.
//!
//! The parameters fit gives three numbers per read group per stratum: how often a read
//! slips at all, which way it slips when it does, and how fast bigger slips fall away
//! ([`Slippage`]). The stutter distribution wants seven
//! ([`StutterRates`] plus the same-length share it derives).
//! [`stutter_rates_for`] is the conversion, and it is the only place in ng that performs it.
//!
//! **The emission itself is not here.** `SsrEmissionModel` and Model A are Milestone F's,
//! in `ssr_emission.rs` beside this file; this module is only the parameter adapter, which
//! is why it is named for what it produces.
//!
//! # Why the conversion lives here and not beside the distribution
//!
//! [`crate::ng::alignment::stutter`] holds the distribution, and its own contract says the
//! type **fits nothing** — it reads seven numbers frozen. The fit lives in
//! [`crate::ng::parameter_estimation`]. The two modules are siblings and neither imports the
//! other; putting the adapter in the alignment module would make the shared *distribution*
//! depend on the *fitting*, which is the direction its doc disclaims. This module already
//! depends on both, and `doc/devel/ng/spec/read_likelihoods.md` §7 puts the reading of frozen
//! parameters on this side of the boundary.
//!
//! # Two of the seven are placeholders, and they are named as such
//!
//! [`PART_REPEAT_SHARE_OF_WHOLE`] and the tying of the two one-step shares are **not fitted
//! results**. Production carries both, ng inherits both, and spec §4.2 requires them recorded
//! as placeholders rather than mistaken for estimates. §10 gives the second a home; the first
//! is waiting for an owner to bin part-repeat reads separately in the pre-pass, as HipSTR
//! does.

use crate::ng::alignment::{StutterModel, StutterRates};
use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;

/// How much stutter mass goes to **part-repeat** changes, as a fraction of the whole-repeat
/// mass — production's `OUT_FRAME_REL`
/// ([`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs)).
///
/// **A placeholder, not an estimate.** Production's own comment calls it "the Step-4 declared
/// estimator … pinned to a real per-period estimate in Step 5", and that estimate was never
/// made. The follow-up that would replace it — binning part-repeat reads separately in the
/// parameter pre-pass, as HipSTR does — is recorded in the alignment spec's §5.2 and has no
/// owner; the parameters fit sends its part-repeat reads to a bucket it calls a diagnostic
/// rather than a parameter (spec §4.2, §10).
///
/// **What relying on it costs:** every part-repeat score in ng is a fixed twentieth of the
/// whole-repeat one, per read group and stratum, whatever the motif. Any comparison that
/// turns on part-repeat reads inherits that.
pub const PART_REPEAT_SHARE_OF_WHOLE: f64 = 0.05;

/// The seven shares of the stutter distribution, from the three numbers the fit produces for
/// one read group in one stratum.
///
/// # The complement, which is the trap this function exists to contain
///
/// [`Slippage::fall_off`] is the probability of **carrying on** to a further step:
/// `read_probabilities` weights step *n* as `(1 − fall_off) · fall_off^(n − 1)`. A
/// [`StutterRates`] one-step share is the probability of **stopping at one step**. **They are
/// complements**, so the conversion is `one_step_share = 1 − fall_off`, and getting it
/// backwards inverts the size distribution — large slips become the common ones — with
/// nothing crashing (spec §4.2's first trap).
///
/// Production makes the same conversion in the same words
/// ([`hipstr.rs`](../../../../src/ssr/cohort/read_model/hipstr.rs): *"B's `decay` is the
/// geometric continuation probability … so the matching conversion is `geom = 1 − decay`"*).
/// Two tests here fail if it is dropped: `the_one_step_share_is_the_complement_of_the_fall_off`
/// and `a_fitted_row_yields_a_distribution_that_sums_to_one`. **The sums-to-one tripwire in
/// `alignment::stutter` does not**, and an earlier version of this sentence said it did —
/// that test builds its rates from literals, so it re-spells this conversion instead of
/// calling it. Measured: dropping the complement fails exactly two tests, both in this
/// module, and the tripwire stays green.
///
/// # How the level is split four ways
///
/// The level is the whole-repeat mass; the part-repeat mass is
/// [`PART_REPEAT_SHARE_OF_WHOLE`] of it **added on top**, not carved out of it — production's
/// shape, so the total slip mass is `1.05 × level`. Each is then split by the fitted
/// direction share. **The two one-step shares are tied to one value**, which HipSTR keeps
/// independent; that is the second placeholder.
///
/// # Clamping is not done here
///
/// [`StutterModel::new`] holds the one-step shares strictly inside `(0, 1)` and floors the
/// derived same-length share. This function does the arithmetic and hands over; a second
/// clamp here would be a second place for the contract to live.
#[must_use]
pub fn stutter_rates_for(slippage: &Slippage) -> StutterRates {
    let whole_repeat_mass = slippage.level;
    let part_repeat_mass = slippage.level * PART_REPEAT_SHARE_OF_WHOLE;
    let shorter_share = slippage.shorter_share;
    let longer_share = 1.0 - shorter_share;

    // The complement — see the doc above. This is the one line the trap is about.
    let one_step_share = 1.0 - slippage.fall_off;

    StutterRates {
        whole_repeat_longer_share: whole_repeat_mass * longer_share,
        whole_repeat_shorter_share: whole_repeat_mass * shorter_share,
        whole_repeat_one_step_share: one_step_share,
        part_repeat_longer_share: part_repeat_mass * longer_share,
        part_repeat_shorter_share: part_repeat_mass * shorter_share,
        // Tied to the whole-repeat share — the second placeholder (spec §10).
        part_repeat_one_step_share: one_step_share,
    }
}

/// The stutter distribution for one read group in one stratum, built from that cell's fit.
///
/// A convenience over [`stutter_rates_for`] plus [`StutterModel::new`], so a caller that
/// wants the distribution cannot accidentally build the model from unconverted numbers.
#[must_use]
pub fn stutter_model_for(slippage: &Slippage) -> StutterModel {
    StutterModel::new(stutter_rates_for(slippage))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A fit whose three numbers are **all different and none of them a half**, so a
    /// transposition or a dropped complement changes an answer.
    fn a_fit() -> Slippage {
        Slippage {
            level: 0.02,
            shorter_share: 0.83,
            fall_off: 0.35,
        }
    }

    /// **The one-step share is the complement of the fall-off** (spec §4.2's first trap).
    /// A fall-off of 0.35 — the chance of carrying on to a further step — is a one-step share
    /// of 0.65. Using 0.35 directly inverts the size distribution and nothing crashes, which
    /// is why the value here is deliberately not 0.5: at a half the two are equal and the
    /// mistake is invisible.
    #[test]
    fn the_one_step_share_is_the_complement_of_the_fall_off() {
        let rates = stutter_rates_for(&a_fit());
        assert!((rates.whole_repeat_one_step_share - 0.65).abs() < 1e-15);
        assert!((rates.part_repeat_one_step_share - 0.65).abs() < 1e-15);

        // And the distribution it builds puts more mass on one repeat than on two, in the
        // ratio the share states — the observable consequence of getting it right.
        let model = stutter_model_for(&a_fit());
        let period = std::num::NonZeroU8::new(3).expect("a valid period");
        let one = model.probability(-3, period);
        let two = model.probability(-6, period);
        assert!((one / two - 1.0 / 0.35).abs() < 1e-9, "{one} against {two}");
    }

    /// **The level is the whole-repeat mass, split by the fitted direction share**, and the
    /// part-repeat mass is a twentieth of it **added on top** rather than carved out — which
    /// is production's shape, so the four direction shares total `1.05 × level`.
    #[test]
    fn the_level_becomes_four_direction_shares_totalling_one_and_a_twentieth_of_it() {
        let fit = a_fit();
        let rates = stutter_rates_for(&fit);

        assert!((rates.whole_repeat_shorter_share - 0.02 * 0.83).abs() < 1e-15);
        assert!((rates.whole_repeat_longer_share - 0.02 * 0.17).abs() < 1e-15);
        assert!((rates.part_repeat_shorter_share - 0.02 * 0.05 * 0.83).abs() < 1e-15);
        assert!((rates.part_repeat_longer_share - 0.02 * 0.05 * 0.17).abs() < 1e-15);

        let total = rates.whole_repeat_longer_share
            + rates.whole_repeat_shorter_share
            + rates.part_repeat_longer_share
            + rates.part_repeat_shorter_share;
        assert!((total - 0.02 * 1.05).abs() < 1e-15, "total {total}");

        // Contraction-biased, because the fit said so — 0.83 shorter against 0.17 longer is
        // tomato's dinucleotide split, 2,438 reads against 501.
        assert!(rates.whole_repeat_shorter_share > rates.whole_repeat_longer_share);
        assert!(rates.part_repeat_shorter_share > rates.part_repeat_longer_share);
    }

    /// **The part-repeat mass is a fixed twentieth of the whole-repeat mass**, and that is a
    /// placeholder rather than a fitted result. Pinned so that replacing it with a real
    /// estimator is a visible change rather than a silent one.
    #[test]
    fn the_part_repeat_share_is_the_declared_placeholder() {
        assert_eq!(PART_REPEAT_SHARE_OF_WHOLE, 0.05);

        let rates = stutter_rates_for(&a_fit());
        let whole = rates.whole_repeat_longer_share + rates.whole_repeat_shorter_share;
        let part = rates.part_repeat_longer_share + rates.part_repeat_shorter_share;
        assert!((part / whole - PART_REPEAT_SHARE_OF_WHOLE).abs() < 1e-15);
    }

    /// **The two one-step shares are tied to one number**, which HipSTR keeps independent.
    /// The second declared placeholder (spec §10): pinned here so untying them is a change
    /// someone makes on purpose.
    #[test]
    fn the_two_one_step_shares_are_tied_to_one_number() {
        for fall_off in [0.05, 0.35, 0.6, 0.9] {
            let rates = stutter_rates_for(&Slippage {
                level: 0.02,
                shorter_share: 0.83,
                fall_off,
            });
            assert_eq!(
                rates.whole_repeat_one_step_share, rates.part_repeat_one_step_share,
                "at a fall-off of {fall_off}"
            );
        }
    }

    /// **A fitted row sums to one over the whole support.**
    ///
    /// The tripwire in `alignment::stutter` spells this conversion out by hand — it builds
    /// its rates from literals — so **no mistake in `stutter_rates_for` can fail it**.
    /// Measured: dropping the complement fails exactly two tests, both in this module, and
    /// leaves the tripwire green. This one goes through the conversion, which is what makes
    /// it a guard on the conversion.
    ///
    /// The fall-offs here are slow enough that the correct conversion leaves only
    /// `fall_off^10` of each branch past the cutoff, while reading the fall-off *as* the
    /// one-step share leaves `(1 − fall_off)^10` — six tenths of a branch at a fall-off of
    /// 0.05.
    #[test]
    fn a_fitted_row_yields_a_distribution_that_sums_to_one() {
        for fall_off in [0.01, 0.05, 0.1] {
            for level in [0.002, 0.02, 0.2] {
                let model = stutter_model_for(&Slippage {
                    level,
                    shorter_share: 0.83,
                    fall_off,
                });
                for period_bases in 2..=6u8 {
                    let period = std::num::NonZeroU8::new(period_bases).expect("a valid period");
                    let total: f64 = (-600i64..=600)
                        .map(|bp_diff| model.probability(bp_diff, period))
                        .sum();
                    assert!(
                        (total - 1.0).abs() < 1e-9,
                        "period {period_bases}, level {level}, fall-off {fall_off}: {total}"
                    );
                }
            }
        }
    }

    /// **The fit's own top level, which nothing covered.**
    ///
    /// A slippage curve may report a level as high as `LEVEL_CEILING` (0.999). Because the
    /// part-repeat mass is added on top rather than carved out, the four direction shares
    /// then total `1.05 × level` — **past one**, at any level above 0.9524. Every individual
    /// share is still under one, so `StutterModel::new`'s debug assertion does not fire; the
    /// same-length share floors at `GEOM_MIN` instead, and `unreachable_mass` reports
    /// **exactly no loss** for a row that has over-allocated six parts in a hundred.
    ///
    /// **This is a row the parameter fit is entitled to emit**, not a hostile hand-built one,
    /// which is why it is pinned here rather than left to the floored-model test in the
    /// stutter module. Whether reporting nothing there is the wanted answer is a decision for
    /// whoever owns the fit's ceiling; this test makes it visible either way.
    #[test]
    fn the_fits_top_level_floors_the_same_length_share_and_reports_no_loss() {
        use crate::ng::parameter_estimation::joint::slippage_curve::LEVEL_CEILING;

        let rates = stutter_rates_for(&Slippage {
            level: LEVEL_CEILING,
            shorter_share: 0.83,
            fall_off: 0.35,
        });
        let four_directions = rates.whole_repeat_longer_share
            + rates.whole_repeat_shorter_share
            + rates.part_repeat_longer_share
            + rates.part_repeat_shorter_share;
        assert!(
            (four_directions - LEVEL_CEILING * 1.05).abs() < 1e-12,
            "{four_directions}"
        );
        assert!(four_directions > 1.0, "{four_directions}");

        let model = StutterModel::new(rates);
        assert_eq!(
            model.same_length_share(),
            crate::ng::calling::likelihood::GEOM_MIN
        );
        let five = model.same_length_share() + four_directions;
        assert!((five - 1.058_95).abs() < 1e-9, "{five}");

        let period = std::num::NonZeroU8::new(3).expect("a valid period");
        let thirty_repeats = std::num::NonZeroU32::new(30).expect("a valid repeat count");
        assert_eq!(model.unreachable_mass(period, thirty_repeats), 0.0);
    }

    /// **Clamping belongs to the model, not here.** A fit at the edges — no slippage at all,
    /// or a fall-off of one — produces rates outside the model's contractual range, and it is
    /// [`StutterModel::new`] that holds them inside it. Asserted so that adding a second
    /// clamp here would be a visible change of ownership.
    #[test]
    fn the_edges_of_the_fit_are_clamped_by_the_model_rather_than_the_conversion() {
        let never_slips = stutter_rates_for(&Slippage {
            level: 0.0,
            shorter_share: 0.83,
            fall_off: 1.0,
        });
        // The conversion passes the edge through untouched...
        assert_eq!(never_slips.whole_repeat_one_step_share, 0.0);
        assert_eq!(never_slips.whole_repeat_longer_share, 0.0);

        // ...and the model holds it inside the range its contract promises.
        let model = StutterModel::new(never_slips);
        assert_eq!(
            model.whole_repeat_one_step_share(),
            crate::ng::calling::likelihood::GEOM_MIN
        );
        // No slippage at all means every read shows the allele's own length.
        assert!((model.same_length_share() - 1.0).abs() < 1e-15);
    }
}
