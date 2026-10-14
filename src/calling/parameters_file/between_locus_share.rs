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
//! | [`concentration_error_of`] | F's standard error back into the concentration's |
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
//!
//! # What survives a trip through calling
//!
//! A run handed a file reads each F as a concentration, scores with it, and writes it back beside
//! its VCF as F. **Plain arithmetic does not bring every F back**: `1 / (1 + (1 − F) / F)` is
//! the neighbouring `f64` for about one F in five drawn uniformly from (0, 1), and for fewer of
//! the Fs converted from concentrations: about one in 37 for concentrations drawn log-uniform
//! over 10⁻¹² to 10⁸, and one in 24 over 10⁻³ to 10⁴, the range the new fit clamps a
//! concentration to (a million draws each, with numpy). So the two conversions *into* the file
//! pick, among the few `f64`s next to the plain answer, one that the reverse conversion maps back
//! exactly, and they pick it the same way every time (nearest the plain answer, the lower of two
//! equally near). Measured by this module's tests:
//!
//! - **A file this build writes comes back byte for byte.** Write F from a concentration, read it
//!   back, write it again: the same F, for every one of 100,000 concentrations drawn log-uniform
//!   over 10⁻¹² to 10⁸, and the same for F's standard error.
//! - **The concentration calling reads back is the fitted one only where some F maps back to it
//!   exactly.** Over the same draws that is 35,369 of 100,000, and 25,553 of 50,000 drawn over
//!   10⁻³ to 10⁴, the range the new fit clamps a concentration to. It depends on the size of the
//!   concentration (measured on 20,000 draws over each hundredfold range, with the same rule written in
//!   Python): almost never below 0.01, about 15 in 100 between 0.01 and 1, and 79 to 86 in 100
//!   above 1. A concentration below 0.01 has an F near 1, where the `f64`s are too sparse for one
//!   to land on it. Where none maps back, the concentration read back is the one the written F
//!   stands for, within about 10⁻¹⁶ relative to `1 + c` ([`concentration_of`]). F's standard
//!   error reads back exactly in 92,040 of the 100,000.
//! - **An F typed by hand comes back unchanged unless a neighbouring F stands for the same
//!   concentration.** Of 100,000 drawn uniformly from (0, 1), 97,437 came back unchanged; each of
//!   the others shares its concentration with the F next to it, and calling, which holds only the
//!   concentration, cannot tell the two apart.

use serde::{Deserialize, Serialize};

/// **F, checked: strictly between 0 and 1, and its concentration `(1 − F) / F` finite.**
///
/// The file's rows and its edge hold this rather than a bare `f64`, so that F and a concentration
/// cannot be swapped without the compiler noticing. Built by [`BetweenLocusShare::new`] or
/// [`between_locus_share_of`]; both refuse a value outside the two bounds the module comment
/// gives.
///
/// **The file spells it as a bare number, and reading one checks it.** Deserialising goes through
/// [`BetweenLocusShare::new`], so a row whose F is 0, 1 or below the lower bound is refused as the
/// file is read, at the line it is written on, and no unchecked F reaches calling.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(try_from = "f64", into = "f64")]
pub struct BetweenLocusShare(f64);

impl TryFrom<f64> for BetweenLocusShare {
    type Error = NotABetweenLocusShare;

    fn try_from(share: f64) -> Result<Self, Self::Error> {
        Self::new(share).ok_or(NotABetweenLocusShare(share))
    }
}

impl From<BetweenLocusShare> for f64 {
    fn from(share: BetweenLocusShare) -> Self {
        share.0
    }
}

/// **A number read where a between-locus share belongs that is not one.** Its message is what a
/// person editing the file reads, at the line the number is on.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NotABetweenLocusShare(pub f64);

impl std::fmt::Display for NotABetweenLocusShare {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            out,
            "{} is not a between-locus share: it lies strictly between 0 and 1 (0 would mean \
             every locus of the class carries its whole variation, 1 that each locus is fixed at \
             one length), and above about 5.6e-309, below which the concentration it stands for \
             is too large to hold",
            self.0
        )
    }
}

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

/// How many `f64`s on each side of a plain conversion are searched for one that converts back
/// exactly. In a Python model of the rule, searching two, four or eight on each side found the
/// same values; four is margin.
const NEIGHBOURS_SEARCHED: u64 = 4;

/// **The positive `f64` nearest `anchor`, at most [`NEIGHBOURS_SEARCHED`] steps away, that
/// `maps_back` accepts** — the anchor itself first, then one step down, one step up, two down and
/// so on, so the choice depends on the anchor and the test alone.
fn the_nearest_that_maps_back(anchor: f64, maps_back: impl Fn(f64) -> bool) -> Option<f64> {
    if !(anchor.is_finite() && anchor > 0.0) {
        return None;
    }
    // For a positive `f64` the bit pattern counts the values in order, so one step in the bits is
    // one step to the next value.
    let bits = anchor.to_bits();
    std::iter::once(bits)
        .chain((1..=NEIGHBOURS_SEARCHED).flat_map(|step| {
            [bits.checked_sub(step), bits.checked_add(step)]
                .into_iter()
                .flatten()
        }))
        .map(f64::from_bits)
        .find(|&candidate| candidate.is_finite() && candidate > 0.0 && maps_back(candidate))
}

/// The F nearest the plain `1 / (1 + c)` whose [`concentration_of`] is exactly `c`, if one is near.
fn the_share_that_maps_back_to(concentration: f64) -> Option<BetweenLocusShare> {
    the_nearest_that_maps_back(1.0 / (1.0 + concentration), |share| {
        BetweenLocusShare::new(share).is_some_and(|share| concentration_of(share) == concentration)
    })
    .map(BetweenLocusShare)
}

/// **The between-locus share of a concentration: `F = 1 / (1 + c)`**, chosen among its nearest
/// `f64`s so that a file survives a trip through calling (the module comment, which gives what
/// was measured).
///
/// - **Where an F near `1 / (1 + c)` converts back to exactly `c`**, it is that F — the nearest
///   such, and the lower of two equally near.
/// - **Where none does**, it is the F the same rule picks for the concentration the plain
///   `1 / (1 + c)` stands for. That F reads back to the same concentration as the plain one, so
///   nothing calling scores with changes, and writing that concentration again gives the same F.
/// - **Where that fails too**, it is the plain `1 / (1 + c)`. No measurement reached this case:
///   not one of 7 million concentrations and Fs drawn across the whole range (the review of this
///   step, which ran the rule on plain `f64`s).
///
/// `None` for a concentration that is not finite and strictly positive, for one at or below
/// 2⁻⁵³ ≈ 1.1×10⁻¹⁶, whose F rounds to exactly 1, and for one so large that its F's own
/// concentration overflows (`f64::MAX` is one: its F is the reciprocal of the largest `f64`, and
/// converting that back gives infinity).
pub fn between_locus_share_of(concentration: f64) -> Option<BetweenLocusShare> {
    if !(concentration.is_finite() && concentration > 0.0) {
        return None;
    }
    let plain = BetweenLocusShare::new(1.0 / (1.0 + concentration))?;
    Some(
        the_share_that_maps_back_to(concentration)
            .or_else(|| the_share_that_maps_back_to(concentration_of(plain)))
            .unwrap_or(plain),
    )
}

/// **The concentration a between-locus share stands for: `c = (1 − F) / F`.**
///
/// Finite and strictly positive for every [`BetweenLocusShare`], by its constructor's check.
///
/// **A concentration converted to F and back is not always the same `f64`.** `1 − F` loses the
/// digits of `F` that lie below `F`'s own precision, so where no F maps back to `c` exactly
/// ([`between_locus_share_of`]) it comes back with an error of about 10⁻¹⁶ relative to `1 + c`,
/// not to `c`: negligible for a concentration of 1 or more, and growing as the concentration
/// nears the 1.1×10⁻¹⁶ bound.
pub fn concentration_of(share: BetweenLocusShare) -> f64 {
    (1.0 - share.0) / share.0
}

/// **The standard error of F, from the concentration and its standard error:
/// `se(F) = se(c) / (1 + c)²`**, written as `se(c) · F · F` with the F
/// [`between_locus_share_of`] writes for `c`, and chosen among its nearest `f64`s as that F is.
///
/// The fit and calling carry the concentration's error; the file carries F's. The conversion is
/// the delta method — the error is scaled by how fast F changes with `c`, which is
/// `dF/dc = −1 / (1 + c)²` — so it holds while the error is small against `1 + c`
/// (`doc/devel/ng/spec/str_fit_to_calling.md` §2.1).
///
/// **Chosen so a file survives a trip through calling**, by the same rule as F: the nearest value
/// whose [`concentration_error_of`], with the F written beside it, is exactly `se(c)`; where none
/// is near, the one the rule picks for the error the plain value reads back as; where that fails
/// too, the plain value. The anchor depends only on F and the error, so reading the file and
/// writing it again gives the same error (the module comment).
///
/// It multiplies by F twice rather than dividing by `(1 + c)²`, because the square overflows to
/// infinity for a concentration above about 1.3×10¹⁵⁴ and would turn a positive error into 0,
/// which the file's validation refuses. It does not check its inputs beyond that: a caller holds
/// a concentration that [`between_locus_share_of`] accepted and an error that is finite and not
/// negative. An error of 0, or a concentration with no F, is converted plainly.
pub fn between_locus_share_error_of(concentration: f64, concentration_error: f64) -> f64 {
    let Some(share) = between_locus_share_of(concentration) else {
        let one_plus_concentration = 1.0 + concentration;
        return concentration_error / one_plus_concentration / one_plus_concentration;
    };
    let the_error_that_maps_back_to = |error_on_the_concentration: f64| {
        the_nearest_that_maps_back(
            error_on_the_concentration * share.0 * share.0,
            |candidate| concentration_error_of(share, candidate) == error_on_the_concentration,
        )
    };
    let plain = concentration_error * share.0 * share.0;
    the_error_that_maps_back_to(concentration_error)
        .or_else(|| the_error_that_maps_back_to(concentration_error_of(share, plain)))
        .unwrap_or(plain)
}

/// **The concentration's standard error, from F and F's standard error:
/// `se(c) = se(F) / F²`** — the reverse of [`between_locus_share_error_of`], by the same delta
/// method, since `dc/dF = −1 / F²`.
///
/// Calling carries the concentration's error and the file F's, so a file read back into a run
/// converts its errors here. Like its partner it does not check its inputs. It divides by F twice
/// rather than once by its square, which loses digits for an F below about 1.5×10⁻¹⁵⁴ and is 0
/// below about 2.2×10⁻¹⁶². For a very small F the result can be infinite; the file's validation
/// refuses such an error.
pub fn concentration_error_of(share: BetweenLocusShare, share_error: f64) -> f64 {
    share_error / share.0 / share.0
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
        // The plain 1 / (1 + c) is 0.9999999999999998; the F one step above it is the one that
        // converts back to exactly this concentration, so that is the one written.
        let share = between_locus_share_of(next_up).expect("just above the bound has an F");
        assert_eq!(share.get(), 0.9999999999999999);
        assert_eq!(concentration_of(share), next_up);
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

    /// At `F = 1/2` an error of 0.1 on F is 0.4 on the concentration, the reverse of the test
    /// above; at `F = 1/4` it is sixteen times F's.
    #[test]
    fn the_concentrations_error_is_fs_scaled_back() {
        let half = BetweenLocusShare::new(0.5).unwrap();
        assert_eq!(concentration_error_of(half, 0.1), 0.4);
        let quarter = BetweenLocusShare::new(0.25).unwrap();
        assert_eq!(concentration_error_of(quarter, 0.1), 1.6);
        for (concentration, error) in [(0.3, 0.05), (4.0, 1.5), (250.0, 30.0)] {
            let share = between_locus_share_of(concentration).unwrap();
            let back =
                concentration_error_of(share, between_locus_share_error_of(concentration, error));
            assert!(
                ((back - error) / error).abs() < 1e-12,
                "({concentration}, {error}) came back as {back}"
            );
        }
    }

    /// **The file spells F as a bare number, and reading one that is not an F refuses it.**
    #[test]
    fn f_is_written_as_a_number_and_read_back_checked() {
        #[derive(Debug, PartialEq, Serialize, Deserialize)]
        struct Row {
            share: BetweenLocusShare,
        }
        let row = Row {
            share: BetweenLocusShare::new(0.25).unwrap(),
        };
        let text = toml::to_string(&row).expect("serialises");
        assert_eq!(text.trim(), "share = 0.25");
        assert_eq!(toml::from_str::<Row>(&text).unwrap(), row);
        for refused in [
            "share = 0.0",
            "share = 1.0",
            "share = 1e-320",
            "share = -0.5",
        ] {
            let error = toml::from_str::<Row>(refused).expect_err(refused);
            assert!(
                error.message().contains("is not a between-locus share"),
                "{refused}: {error}"
            );
        }
    }

    /// A fixed-seed generator for the measurements below (SplitMix64), so their counts are the
    /// same on every run and can be pinned.
    struct Draws(u64);

    impl Draws {
        /// A number uniform in [0, 1).
        fn uniform(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            (z >> 11) as f64 / (1_u64 << 53) as f64
        }

        /// A number log-uniform between `low` and `high`.
        fn log_uniform(&mut self, low: f64, high: f64) -> f64 {
            // The crate's own `powf` and `log10`, which round alike on every platform, so the
            // pinned counts below are the same everywhere.
            let (low, high) = (crate::float::log10(low), crate::float::log10(high));
            crate::float::powf(10.0, low + self.uniform() * (high - low))
        }
    }

    /// **A file this build writes comes back byte for byte**, F and its error, for every one of
    /// 100,000 concentrations drawn log-uniform over 10⁻¹² to 10⁸, each with an error between a
    /// thousandth of the concentration and ten times it: write F from the concentration, read it
    /// back as calling does, write it again. And how often calling reads back the fitted number
    /// itself, pinned at what was measured: the concentration in 35,369 draws of 100,000 (in
    /// 25,553 of 50,000 drawn over 10⁻³ to 10⁴, the range the new fit clamps to), the error in
    /// 92,040.
    #[test]
    fn a_file_this_build_writes_comes_back_byte_for_byte() {
        let mut draws = Draws(20_261_014);
        let (mut concentration_exact, mut error_exact) = (0_u32, 0_u32);
        for _ in 0..100_000 {
            let concentration = draws.log_uniform(1e-12, 1e8);
            let error = concentration * draws.log_uniform(1e-3, 10.0);

            let written = between_locus_share_of(concentration).expect("inside the bounds");
            let written_error = between_locus_share_error_of(concentration, error);
            let read = concentration_of(written);
            let read_error = concentration_error_of(written, written_error);
            let again = between_locus_share_of(read).expect("a read F has a concentration");
            let again_error = between_locus_share_error_of(read, read_error);

            assert_eq!(
                again.get().to_bits(),
                written.get().to_bits(),
                "{concentration:e}: F {} came back as {}",
                written.get(),
                again.get()
            );
            assert_eq!(
                again_error.to_bits(),
                written_error.to_bits(),
                "{concentration:e} ± {error:e}: error {written_error:e} came back as \
                 {again_error:e}"
            );
            concentration_exact += u32::from(read == concentration);
            error_exact += u32::from(read_error == error);
        }

        let mut in_the_clamped_range = 0_u32;
        for _ in 0..50_000 {
            let concentration = draws.log_uniform(1e-3, 1e4);
            let written = between_locus_share_of(concentration).expect("inside the bounds");
            in_the_clamped_range += u32::from(concentration_of(written) == concentration);
        }
        assert_eq!(
            (concentration_exact, error_exact, in_the_clamped_range),
            (35_369, 92_040, 25_553),
            "how often calling reads back the fitted concentration, its error, and the \
             concentration over 1e-3..1e4"
        );
    }

    /// **An F typed by hand comes back unchanged unless the F next to it stands for the same
    /// concentration.** Of 100,000 drawn uniformly from (0, 1), 97,437 come back unchanged;
    /// each of the others comes back as a neighbour with the same concentration, which is all
    /// calling holds, so no conversion could tell them apart.
    #[test]
    fn an_f_typed_by_hand_comes_back_unless_its_neighbour_shares_its_concentration() {
        let mut draws = Draws(7);
        let mut unchanged = 0_u32;
        for _ in 0..100_000 {
            let Some(typed) = BetweenLocusShare::new(draws.uniform()) else {
                continue;
            };
            let back = between_locus_share_of(concentration_of(typed)).expect("an F's own");
            if back == typed {
                unchanged += 1;
                continue;
            }
            assert_eq!(
                concentration_of(back),
                concentration_of(typed),
                "{} came back as {}, which stands for another concentration",
                typed.get(),
                back.get()
            );
            assert!(
                back.get().to_bits().abs_diff(typed.get().to_bits()) <= NEIGHBOURS_SEARCHED,
                "{} came back as {}, which is not its neighbour",
                typed.get(),
                back.get()
            );
        }
        assert_eq!(unchanged, 97_437);
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
