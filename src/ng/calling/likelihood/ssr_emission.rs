//! The STR emission — and, first, the seven shares it is built from.
//!
//! The parameters fit gives three numbers per read group per stratum: how often a read
//! slips at all, which way it slips when it does, and how fast bigger slips fall away
//! ([`Slippage`]). The stutter distribution wants seven
//! ([`StutterRates`] plus the same-length share it derives).
//! [`stutter_rates_for`] is the conversion, and it is the only place in ng that performs it.
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
/// `the_one_step_share_is_the_complement_of_the_fall_off` fails if it is dropped, and
/// `the_distribution_sums_to_one_over_its_whole_support` in the stutter module fails too,
/// because a complemented share loses six tenths of a branch's mass past the cutoff.
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
    let shorter = slippage.shorter_share;
    let longer = 1.0 - shorter;

    // The complement — see the doc above. This is the one line the trap is about.
    let one_step_share = 1.0 - slippage.fall_off;

    StutterRates {
        whole_repeat_longer_share: whole_repeat_mass * longer,
        whole_repeat_shorter_share: whole_repeat_mass * shorter,
        whole_repeat_one_step_share: one_step_share,
        part_repeat_longer_share: part_repeat_mass * longer,
        part_repeat_shorter_share: part_repeat_mass * shorter,
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
