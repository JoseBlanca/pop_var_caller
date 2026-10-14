//! **The between-locus share F, and the only place it is converted to and from a concentration.**
//!
//! A repeat-tract class's allele length variation is a Dirichlet over length offsets. Its
//! *concentration* `c` says how many chromosomes' worth of belief the class's shares carry: a
//! large `c` means every locus of the class shows nearly the same spread of allele lengths, a
//! small one that loci differ widely. The fit hands calling a concentration, and calling works in
//! it. The parameters file carries the same number on a scale from 0 to 1, the **between-locus
//! share**
//! `F = 1 / (1 + c)`: the share of the variation in allele lengths that lies between loci rather
//! than within one (`doc/devel/ng/spec/str_fit_to_calling.md` §2.1).
//!
//! **The file's F is converted here and nowhere else.** No other code writes F into the file or
//! turns the file's F back into a concentration, so the two directions cannot drift apart. The
//! fits compute `1 / (1 + c)` internally as a working quantity (`penalised_fit.rs`, which
//! optimises on the logit of F, and `composite_fit.rs`); that is not the file's value. Three
//! functions cross the file's edge:
//!
//! | function | what it converts |
//! |---|---|
//! | [`between_locus_share_of`] | a concentration into F, refusing one F cannot hold |
//! | [`concentration_of`] | F back into a concentration |
//! | [`between_locus_share_error_of`] | a concentration's standard error into F's |
//!
//! **Not every positive concentration has an F, and not every F between 0 and 1 has a
//! concentration.** Floating-point arithmetic bounds both ends:
//!
//! - a concentration at or below 2⁻⁵³ ≈ 1.1×10⁻¹⁶ gives an F that rounds to exactly 1, which
//!   stands for a concentration of zero;
//! - an F at or below about 5.6×10⁻³⁰⁹ (the reciprocal of the largest `f64`) converts back to an
//!   infinite concentration, which calling refuses (`AlleleLengthVariation::fitted` asserts a
//!   finite one).
//!
//! [`BetweenLocusShare`] holds only an F that is inside both bounds.

/// **F, checked: strictly between 0 and 1, and its concentration `(1 − F) / F` finite.**
///
/// The file's rows and its edge hold this rather than a bare `f64`, so that F and a concentration
/// cannot be swapped without the compiler noticing. Built by [`BetweenLocusShare::new`] or
/// [`between_locus_share_of`]; both refuse a value outside the two bounds the module comment
/// gives.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetweenLocusShare(f64);

impl BetweenLocusShare {
    /// The share, or `None` when it is not strictly between 0 and 1 (NaN included) or so close
    /// to 0 that its concentration `(1 − F) / F` overflows to infinity.
    pub fn new(share: f64) -> Option<Self> {
        let strictly_between_0_and_1 = share > 0.0 && share < 1.0;
        if strictly_between_0_and_1 && ((1.0 - share) / share).is_finite() {
            Some(Self(share))
        } else {
            None
        }
    }

    /// The share as a number strictly between 0 and 1.
    pub fn get(self) -> f64 {
        self.0
    }
}

/// **The between-locus share of a concentration: `F = 1 / (1 + c)`.**
///
/// `None` for a concentration that is not finite and strictly positive, for one at or below
/// 2⁻⁵³ ≈ 1.1×10⁻¹⁶, whose F rounds to exactly 1, and for one so large that its F's own
/// concentration overflows (`f64::MAX` is one: its F is the reciprocal of the largest `f64`, and
/// converting that back gives infinity).
pub fn between_locus_share_of(concentration: f64) -> Option<BetweenLocusShare> {
    if !(concentration.is_finite() && concentration > 0.0) {
        return None;
    }
    BetweenLocusShare::new(1.0 / (1.0 + concentration))
}

/// **The concentration a between-locus share stands for: `c = (1 − F) / F`.**
///
/// Finite and strictly positive for every [`BetweenLocusShare`], by its constructor's check.
///
/// **A round trip is not exact.** `1 − F` loses the digits of `F` that lie below `F`'s own
/// precision, so a concentration converted to F and back comes back with an error of about
/// 10⁻¹⁶ relative to `1 + c`, not to `c`: negligible for a concentration of 1 or more, and
/// growing as the concentration nears the 1.1×10⁻¹⁶ bound.
pub fn concentration_of(share: BetweenLocusShare) -> f64 {
    (1.0 - share.0) / share.0
}

/// **The standard error of F, from the concentration and its standard error:
/// `se(F) = se(c) / (1 + c)²`.**
///
/// The fit and calling carry the concentration's error; the file carries F's. The conversion is
/// the delta method — the error is scaled by how fast F changes with `c`, which is
/// `dF/dc = −1 / (1 + c)²` — so it holds while the error is small against `1 + c`
/// (`doc/devel/ng/spec/str_fit_to_calling.md` §2.1). It does not check its inputs: a caller
/// holds a concentration that [`between_locus_share_of`] accepted and an error that is finite
/// and not negative.
///
/// It divides by `1 + c` twice rather than once by its square, because the square overflows to
/// infinity for a concentration above about 1.3×10¹⁵⁴ and would turn a positive error into 0,
/// which the file's validation refuses. Dividing twice keeps the error positive until F's error
/// itself is too small for an `f64`.
pub fn between_locus_share_error_of(concentration: f64, concentration_error: f64) -> f64 {
    let one_plus_concentration = 1.0 + concentration;
    concentration_error / one_plus_concentration / one_plus_concentration
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 2⁻⁵³, the largest concentration whose F rounds to exactly 1: `1 + 2⁻⁵³` rounds to 1.
    const LARGEST_CONCENTRATION_WITHOUT_AN_F: f64 = 1.1102230246251565e-16;

    #[test]
    fn a_concentration_converted_to_f_and_back_comes_back_within_rounding() {
        for concentration in [1e-12, 0.01, 0.5, 1.0, 4.0, 30.0, 500.0, 1e6, 1e300] {
            let share = between_locus_share_of(concentration)
                .unwrap_or_else(|| panic!("{concentration} is a valid concentration"));
            let back = concentration_of(share);
            // The error is relative to 1 + c, because 1 − F has only the digits F leaves it.
            assert!(
                (back - concentration).abs() <= 4.0 * f64::EPSILON * (1.0 + concentration),
                "{concentration} came back as {back}"
            );
        }
    }

    #[test]
    fn an_f_converted_to_a_concentration_and_back_comes_back_within_rounding() {
        for share in [1e-300, 1e-6, 0.1, 0.5, 0.9, 0.999_999] {
            let checked = BetweenLocusShare::new(share).expect("a share inside both bounds");
            let back = between_locus_share_of(concentration_of(checked))
                .expect("the concentration of a checked share has a share");
            assert!(
                (back.get() - share).abs() <= 4.0 * f64::EPSILON * share,
                "{share} came back as {}",
                back.get()
            );
        }
    }

    #[test]
    fn known_values_convert_exactly() {
        assert_eq!(
            between_locus_share_of(1.0).map(BetweenLocusShare::get),
            Some(0.5)
        );
        assert_eq!(
            between_locus_share_of(3.0).map(BetweenLocusShare::get),
            Some(0.25)
        );
        assert_eq!(concentration_of(BetweenLocusShare::new(0.2).unwrap()), 4.0);
    }

    /// **The mean of two Fs converted back is the closed form the stated concentration uses.**
    ///
    /// Calling's stated concentration at an even count of classes is the mean of the middle two
    /// Fs converted back. It is computed in the closed form `(c1 + c2 + 2·c1·c2) / (2 + c1 + c2)`
    /// (`median_concentration`, `parameter_estimation/joint/stratum_fits.rs`), which cannot
    /// cancel to 0 for tiny concentrations. Concentrations 4 and 6 have Fs 1/5 and 1/7, whose
    /// mean 6/35 is a concentration of 29/6 — the value that function's test expects.
    #[test]
    fn the_mean_of_two_shares_converted_back_matches_the_closed_form() {
        for (c1, c2) in [(4.0, 6.0), (1.0, 6.0), (0.3, 250.0), (1e-3, 1e3)] {
            let mean_share = (between_locus_share_of(c1).unwrap().get()
                + between_locus_share_of(c2).unwrap().get())
                / 2.0;
            let round_trip = concentration_of(BetweenLocusShare::new(mean_share).unwrap());
            let closed_form = (c1 + c2 + 2.0 * c1 * c2) / (2.0 + c1 + c2);
            assert!(
                (round_trip - closed_form).abs() <= 1e-12 * (1.0 + closed_form),
                "({c1}, {c2}): round trip {round_trip}, closed form {closed_form}"
            );
        }
        let mean_share = (0.2 + 1.0 / 7.0) / 2.0;
        let round_trip = concentration_of(BetweenLocusShare::new(mean_share).unwrap());
        assert!((round_trip - 29.0 / 6.0).abs() < 1e-14, "got {round_trip}");
    }

    #[test]
    fn an_f_of_0_or_1_or_outside_them_is_refused() {
        for share in [
            0.0,
            1.0,
            -0.1,
            1.5,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ] {
            assert_eq!(BetweenLocusShare::new(share), None, "{share} is not an F");
        }
    }

    /// **An F below about 5.6×10⁻³⁰⁹ passes "between 0 and 1" but has no finite
    /// concentration.** The reciprocal of the largest `f64`, 5.562684646268003e-309, converts
    /// back to infinity, as does 5.56e-309; 5.57e-309 converts back to about 1.795e308.
    #[test]
    fn an_f_whose_concentration_overflows_is_refused() {
        let reciprocal_of_the_largest = 1.0 / f64::MAX;
        assert_eq!(reciprocal_of_the_largest, 5.562684646268003e-309);
        assert_eq!(BetweenLocusShare::new(reciprocal_of_the_largest), None);
        assert_eq!(BetweenLocusShare::new(5.56e-309), None);
        // The smallest positive `f64`, 5e-324, a subnormal.
        assert_eq!(BetweenLocusShare::new(f64::from_bits(1)), None);

        let just_above = BetweenLocusShare::new(5.57e-309).expect("5.57e-309 has a concentration");
        assert!(concentration_of(just_above).is_finite());
        // The smallest normal `f64`, 2.2e-308, is above the bound.
        assert!(BetweenLocusShare::new(f64::MIN_POSITIVE).is_some());
    }

    /// **A concentration at or below 2⁻⁵³ ≈ 1.1×10⁻¹⁶ has no F**, because `1 / (1 + c)` rounds
    /// to exactly 1. The next `f64` above 2⁻⁵³, 1.1102230246251568e-16, gives
    /// 0.9999999999999998 and is accepted.
    #[test]
    fn a_concentration_whose_f_rounds_to_1_is_refused() {
        // 2⁻⁵³ built from its bits: exponent field 1023 − 53, mantissa zero.
        assert_eq!(
            LARGEST_CONCENTRATION_WITHOUT_AN_F,
            f64::from_bits((1023 - 53) << 52)
        );
        assert_eq!(
            between_locus_share_of(LARGEST_CONCENTRATION_WITHOUT_AN_F),
            None
        );
        assert_eq!(between_locus_share_of(1e-17), None);

        let next_up = f64::from_bits(LARGEST_CONCENTRATION_WITHOUT_AN_F.to_bits() + 1);
        assert_eq!(next_up, 1.1102230246251568e-16);
        assert_eq!(
            between_locus_share_of(next_up).map(BetweenLocusShare::get),
            Some(0.9999999999999998)
        );
    }

    #[test]
    fn a_concentration_that_is_not_finite_and_positive_or_too_large_is_refused() {
        for concentration in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
            assert_eq!(
                between_locus_share_of(concentration),
                None,
                "{concentration} has no F"
            );
        }
        assert!(between_locus_share_of(1e308).is_some());
    }

    /// At `c = 1`, `dF/dc = −1/4`, so an error of 0.4 on the concentration is 0.1 on F; at
    /// `c = 3` it is a sixteenth.
    #[test]
    fn the_error_is_scaled_by_the_slope_of_f() {
        assert_eq!(between_locus_share_error_of(1.0, 0.4), 0.1);
        assert_eq!(between_locus_share_error_of(3.0, 1.6), 0.1);
        assert_eq!(between_locus_share_error_of(4.0, 0.0), 0.0);
    }

    /// **A huge concentration keeps a positive error.** At `c = 1e200`, `(1 + c)²` would be 1e400,
    /// past the largest `f64`, and an error of 1e200 would come out as 0. Divided twice it comes
    /// out as 1e-200, the true value.
    #[test]
    fn the_error_at_a_huge_concentration_does_not_underflow_to_0() {
        let error = between_locus_share_error_of(1e200, 1e200);
        assert!(error > 0.0, "got {error}");
        assert!(((error - 1e-200) / 1e-200).abs() < 1e-12, "got {error}");
    }
}
