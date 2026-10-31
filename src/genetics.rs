//! The population-genetics arithmetic three parts of the caller share: how likely each
//! genotype is before any read is looked at, and the two numerical primitives that
//! calculation needs.
//!
//! What lives here is what more than one module uses. The genotype priors
//! ([`calling::genotype_prior`](crate::calling::genotype_prior)) build a site's prior from
//! `alpha_from_diversity` and [`wright_genotype_log_priors`]; the site quality
//! ([`calling::quality`](crate::calling::quality)) needs [`lgamma`] for the same
//! Beta-Binomial the prior is built on; and the hidden-duplication filter
//! ([`paralog`](crate::paralog)) walks a grid of allele frequencies with
//! [`sfs_grid_point`] and scores each one with [`wright_genotype_log_priors`]. A name three
//! modules share is vocabulary, which is why it sits here beside [`types`](crate::types)
//! rather than inside whichever module happened to need it first.
//!
//! **Two words a reader needs before the rest.** *Diversity* (θ) is how often two randomly
//! chosen chromosomes from the population differ at a base — about 1 in 1,000 in humans.
//! The *inbreeding coefficient* (F) is how much likelier than chance it is that a sample's
//! two copies of a locus are identical by descent: 0 in a randomly mating population, and
//! approaching 1 in a selfing line. The first says how much variation to expect at all; the
//! second says how much of it shows up as homozygotes.
//!
//! # This is production's arithmetic, copied
//!
//! Every item below is `src/genetics.rs` unchanged, doc comments included, because these are
//! the numbers ng was graded against while both callers existed. What is not copied is
//! production's `dirichlet_multinomial_log_priors`: ng has its own in
//! [`calling::genotype_prior::dirichlet_multinomial`](crate::calling::genotype_prior),
//! and production's survives only as that port's test oracle.
//!
//! One item came later and from elsewhere: [`regularised_incomplete_beta`], moved here unchanged
//! from the site quality's artifact correction on 2026-10-17, when the joint fit's report of what
//! the data determine needed it too.

use crate::float;

/// A tiny floor clamping an analytically-positive probability away from
/// `ln(0) = −∞`. The clamp is equivalent to `f64` precision for any realistic
/// probability and never bites on the default grids.
pub const PROBABILITY_FLOOR: f64 = 1e-300;

/// Natural log of the gamma function, `ln Γ(x)`, for `x > 0`.
///
/// The Dirichlet-multinomial genotype prior needs `ln Γ` at **non-integer**
/// arguments (the concentration `α` is derived from the continuous diversity
/// `θ`), which the integer-only factorial helpers cannot supply. This is a thin
/// wrapper over `libm::lgamma` (the rust-lang port of musl's libm), kept here so
/// the rest of the caller depends on one project-owned name and the accuracy
/// tests live beside it.
///
/// Only ever called with `x > 0` in this crate (Dirichlet concentrations and
/// `α + k` are strictly positive). `libm::lgamma` is defined for other inputs,
/// but this wrapper documents the contract the callers actually rely on.
///
/// The `x > 0` precondition is checked with a debug assertion (zero-cost in
/// release): `libm::lgamma` does not panic on a non-positive argument, so a
/// wiring regression that let `α` reach `0` or go negative would otherwise
/// silently produce a corrupt log-prior. The bare `libm` call in release keeps
/// the hot path allocation- and branch-free.
pub fn lgamma(x: f64) -> f64 {
    debug_assert!(x > 0.0, "lgamma requires x > 0, got {x}");
    libm::lgamma(x)
}

/// **The regularised incomplete Beta function `I_x(a, b)`**: the probability that a value drawn
/// from `Beta(a, b)` is at most `x`. Two modules read it: the site quality's artifact correction,
/// for which it is the exact binomial cumulative distribution
/// ([`calling::quality`](crate::calling::quality)), and the joint fit's report of what the data
/// determine, for which one minus it is the share of a Beta above a stated frequency
/// (`doc/devel/ng/spec/parameter_prepass_joint_fit.md` §2.1.3).
///
/// `x` at or below 0 gives 0 and at or above 1 gives 1; both shapes must be above zero. It is
/// Numerical Recipes' `betai` (Press et al., §6.4), over [`lgamma`] and the continued fraction
/// below, each evaluated on the side of the Beta's bulk where the fraction converges fast — for `x`
/// past `(a + 1)/(a + b + 2)` through the symmetry `I_x(a, b) = 1 − I_{1−x}(b, a)`.
///
/// **Moved here from the artifact correction on 2026-10-17, its arithmetic unchanged**, so that the
/// site quality's bits do not move. Production's `lgamma` was a hand-written Lanczos approximation
/// and this one is `libm`'s, so the two are not bit-identical by construction; they agree to about
/// one part in `1e-13`, and what the agreement costs at the end of the whole correction is measured
/// by the differential the correction's port is judged against. `ln` and `exp` go through
/// [`crate::float`], so the same arguments give the same bits on every platform.
pub fn regularised_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    debug_assert!(
        a > 0.0 && b > 0.0,
        "Beta shapes must be above zero, got {a} and {b}"
    );
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let log_front =
        lgamma(a + b) - lgamma(a) - lgamma(b) + a * float::ln(x) + b * float::ln(1.0 - x);
    let front = float::exp(log_front);
    // Each branch is the one the continued fraction converges quickly on; the other is reached
    // through the symmetry `I_x(a, b) = 1 − I_{1−x}(b, a)`.
    if x < (a + 1.0) / (a + b + 2.0) {
        front * incomplete_beta_continued_fraction(a, b, x) / a
    } else {
        1.0 - front * incomplete_beta_continued_fraction(b, a, 1.0 - x) / b
    }
}

/// Lentz's continued fraction for the incomplete beta (Numerical Recipes' `betacf`).
///
/// It converges to about one part in `1e-15` within a few dozen iterations **whatever the
/// magnitude of `a` and `b`**, which is what makes the tail above cost the same at 30 reads and
/// at 30,000.
fn incomplete_beta_continued_fraction(a: f64, b: f64, x: f64) -> f64 {
    const MAX_ITERATIONS: usize = 300;
    const CONVERGED: f64 = 1e-15;
    /// Smallest magnitude the recurrence is allowed to carry, so a denominator that lands on
    /// zero is nudged off it rather than producing an infinity.
    const FLOOR: f64 = 1e-300;

    let a_plus_b = a + b;
    let a_plus_one = a + 1.0;
    let a_minus_one = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - a_plus_b * x / a_plus_one;
    if d.abs() < FLOOR {
        d = FLOOR;
    }
    d = 1.0 / d;
    let mut fraction = d;
    for iteration in 1..=MAX_ITERATIONS {
        let step = iteration as f64;
        let two_steps = 2.0 * step;

        // The even half of the recurrence.
        let term = step * (b - step) * x / ((a_minus_one + two_steps) * (a + two_steps));
        d = 1.0 + term * d;
        if d.abs() < FLOOR {
            d = FLOOR;
        }
        c = 1.0 + term / c;
        if c.abs() < FLOOR {
            c = FLOOR;
        }
        d = 1.0 / d;
        fraction *= d * c;

        // The odd half.
        let term =
            -(a + step) * (a_plus_b + step) * x / ((a + two_steps) * (a_plus_one + two_steps));
        d = 1.0 + term * d;
        if d.abs() < FLOOR {
            d = FLOOR;
        }
        c = 1.0 + term / c;
        if c.abs() < FLOOR {
            c = FLOOR;
        }
        d = 1.0 / d;
        let delta = d * c;
        fraction *= delta;
        if (delta - 1.0).abs() < CONVERGED {
            break;
        }
    }
    fraction
}

/// The `i`-th of `n` points on the inclusive linear grid `[lo, hi]`. `n <= 1`
/// yields the midpoint (a degenerate grid still returns a usable point).
pub fn linear_grid_point(i: usize, n: usize, lo: f64, hi: f64) -> f64 {
    if n <= 1 {
        return 0.5 * (lo + hi);
    }
    lo + (hi - lo) * (i as f64) / ((n - 1) as f64)
}

/// The `i`-th of `n` folded-SFS grid points on `[grid_inset, 1−grid_inset]`,
/// uniform in `p`. A single-point grid collapses to the midpoint `0.5`.
/// `grid_inset` is the amount the endpoints are pulled in from `0` and `1`
/// (typically `1/(2n)`) so no grid point sits at a degenerate frequency.
pub fn sfs_grid_point(i: usize, n: usize, grid_inset: f64) -> f64 {
    linear_grid_point(i, n, grid_inset, 1.0 - grid_inset)
}

/// The Wright inbreeding-adjusted Hardy–Weinberg genotype log-priors
/// `(hom-ref, het, hom-alt)` at ALT frequency `p` and inbreeding coefficient
/// `f`: `P(het) = 2pq(1−f)`, homozygotes `q²+f·pq` / `p²+f·pq` (`q = 1−p`).
/// Each probability is floored at [`PROBABILITY_FLOOR`] before the `ln` so a
/// zero-probability genotype yields a finite, very negative log-prior rather
/// than `−∞`.
pub fn wright_genotype_log_priors(p: f64, f: f64) -> (f64, f64, f64) {
    let q = 1.0 - p;
    let het = 2.0 * p * q * (1.0 - f);
    let hom_ref = q * q + f * p * q;
    let hom_alt = p * p + f * p * q;
    (
        float::ln(hom_ref.max(PROBABILITY_FLOOR)),
        float::ln(het.max(PROBABILITY_FLOOR)),
        float::ln(hom_alt.max(PROBABILITY_FLOOR)),
    )
}

#[cfg(test)]
/// The reference-allele Dirichlet concentration `α_ref`. Fixed at `1`, the value
/// that makes the biallelic-diploid het:hom-alt ratio `2·α_ref/(α_alt+1)`
/// approach the defensible **2:1** as `α_alt → 0`. It doubles as the
/// monomorphic-site weight (arch §9.2): with the small `α_alt = θ̂` from
/// [`alpha_from_diversity`] the Dirichlet-multinomial's hom-ref probability comes
/// out at the genetically-correct `1 − 3θ/2`, so no separate invariant mass is
/// needed at the default.
pub const ALPHA_REF: f64 = 1.0;

/// A tiny positive floor for each ALT concentration, so the
/// Dirichlet-multinomial's `lgamma(α_alt)` stays finite when the estimated
/// diversity is exactly zero (a fully invariant cohort, or `--diversity 0`). It
/// sits far below any real diversity (human `θ ≈ 1e-3`), so it never perturbs a
/// genuine estimate; at `θ = 0` it yields an effectively-certain hom-ref prior,
/// matching the biallelic grid path's `θ = 0` behaviour.
pub const MIN_ALT_CONCENTRATION: f64 = 1e-12;

#[cfg(test)]
/// The Dirichlet concentration `α = (α_ref, α_alt(1), …, α_alt(k−1))` for the SFS
/// genotype prior at estimated diversity `theta` (`θ̂`).
///
/// - `α_ref = ALPHA_REF = 1`.
/// - The total ALT concentration is `θ̂`, split evenly across the `n_alleles − 1`
///   ALT alleles: `α_alt(a) = θ̂ / (n_alleles − 1)`, floored at
///   [`MIN_ALT_CONCENTRATION`] so it stays strictly positive. Splitting keeps a
///   site's total polymorphism `θ̂` independent of how many ALT alleles it
///   carries.
///
/// Fed to ng's Dirichlet-multinomial
/// ([`calling::genotype_prior::dirichlet_multinomial`](crate::calling::genotype_prior)),
/// this yields the clean
/// population-genetics marginals for a biallelic-diploid site (`F = 0`): `P(het)
/// ≈ θ`, `P(hom-alt) ≈ θ/2`, monomorphic weight `≈ 1 − 3θ/2`, and a het:hom-alt
/// ratio that stays `≈ 2:1` at every realistic diversity (because `θ̂` — hence
/// `α_alt` — is always small). Per-sample inbreeding `F` is applied on top by the
/// engine's Wright mixture, not here. See the SFS-prior architecture doc §9.2 for
/// why this is the settled mapping (no calibration constant).
///
/// # Preconditions
///
/// `n_alleles >= 1` is a hard assertion (a zero-allele shape is impossible — every
/// site has a reference allele — and would flow a wrong-length `α` into the
/// Dirichlet-multinomial). A monomorphic shape (`n_alleles == 1`) has no ALT to
/// carry diversity and returns `[ALPHA_REF]`. `theta` finite and `>= 0` is a
/// `debug_assert` (a bad θ degrades the prior but cannot mis-shape `α`).
pub fn alpha_from_diversity(n_alleles: usize, theta: f64) -> Vec<f64> {
    assert!(n_alleles >= 1, "n_alleles must be >= 1");
    debug_assert!(
        theta.is_finite() && theta >= 0.0,
        "theta must be finite and non-negative, got {theta}"
    );

    let n_alt = n_alleles - 1;
    if n_alt == 0 {
        return vec![ALPHA_REF];
    }
    let per_alt = (theta / n_alt as f64).max(MIN_ALT_CONCENTRATION);
    let mut alpha = Vec::with_capacity(n_alleles);
    alpha.push(ALPHA_REF);
    alpha.resize(n_alleles, per_alt);
    alpha
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `I_x(a, b)` pinned from an independent implementation: SciPy 1.18.1's
    /// `scipy.special.betainc(a, b, x)`, run on 2026-10-17 (`tmp/impl_CD/betainc_values.py`).
    /// `x` = 0.01 at shapes 0.02, 1 and 50 — the joint fit's two bounds and the middle — in every
    /// pairing, plus the shapes the fits measured in the spec reach (0.158 and 1.408, a carrier
    /// Beta; 17.9, a density with every segregating position above 1%) and four points past the
    /// Beta's bulk, which go through the symmetry.
    const SCIPY_BETAINC: [(f64, f64, f64, f64); 16] = [
        (0.01, 0.02, 0.02, 0.45638517158324815),
        (0.01, 0.02, 1.0, 0.9120108393559098),
        (0.01, 0.02, 4.28, 0.9466555925383774),
        (0.01, 0.02, 50.0, 0.9886315492736919),
        (0.01, 1.0, 0.02, 0.0002009865165733764),
        (0.01, 1.0, 1.0, 0.01),
        (0.01, 1.0, 50.0, 0.39499393286246337),
        (0.01, 50.0, 0.02, 4.4159013553768126e-104),
        (0.01, 50.0, 1.0, 1.000000000000001e-100),
        (0.01, 50.0, 50.0, 3.113021933716422e-72),
        (0.01, 0.158, 1.408, 0.520803752881817),
        (0.01, 17.9, 3.0, 2.9242835052820803e-34),
        (0.5, 2.0, 3.0, 0.6875),
        (0.99, 0.02, 4.28, 0.9999999999865253),
        (0.3, 50.0, 50.0, 1.554779976144308e-05),
        (0.7, 50.0, 50.0, 0.9999844522002386),
    ];

    /// **The incomplete Beta agrees with SciPy's to 10⁻¹² relative** at every pinned point, values
    /// from 10⁻¹⁰⁴ to one.
    #[test]
    fn the_incomplete_beta_agrees_with_scipy() {
        for (x, a, b, expected) in SCIPY_BETAINC {
            let got = regularised_incomplete_beta(x, a, b);
            assert!(
                (got - expected).abs() <= 1e-12 * expected,
                "I_{x}({a}, {b}) = {got:e}, SciPy {expected:e}"
            );
        }
    }

    /// **`I_x(a, b) = 1 − I_{1−x}(b, a)`**, on both sides of the bulk where the two terms are
    /// computed by different branches, to 10⁻¹⁴.
    #[test]
    fn the_incomplete_beta_is_symmetric() {
        for (x, a, b) in [
            (0.01, 0.02, 4.28),
            (0.2, 2.0, 3.0),
            (0.5, 0.02, 0.02),
            (0.9, 50.0, 1.0),
            (0.45, 50.0, 50.0),
            (0.99, 1.0, 50.0),
        ] {
            let direct = regularised_incomplete_beta(x, a, b);
            let mirrored = 1.0 - regularised_incomplete_beta(1.0 - x, b, a);
            assert!(
                (direct - mirrored).abs() < 1e-14,
                "I_{x}({a}, {b}) = {direct:e} against 1 − I_(1−x)({b}, {a}) = {mirrored:e}"
            );
        }
    }

    /// **At and past the ends**: 0 at `x` = 0 and below, 1 at `x` = 1 and above, whatever the
    /// shapes; not a number for a `x` that is not one.
    #[test]
    fn the_incomplete_beta_is_zero_and_one_at_the_ends() {
        for (a, b) in [(0.02, 0.02), (0.02, 50.0), (50.0, 0.02), (1.0, 1.0)] {
            assert_eq!(regularised_incomplete_beta(0.0, a, b), 0.0);
            assert_eq!(regularised_incomplete_beta(-0.5, a, b), 0.0);
            assert_eq!(regularised_incomplete_beta(1.0, a, b), 1.0);
            assert_eq!(regularised_incomplete_beta(1.5, a, b), 1.0);
            assert!(regularised_incomplete_beta(f64::NAN, a, b).is_nan());
        }
    }

    /// **Just below the fit's lower shape bound it is still defined**: a central difference at a
    /// shape on its bound of 0.02 evaluates it at 0.02 less a step of 10⁻⁵ of it.
    #[test]
    fn the_incomplete_beta_is_defined_just_below_the_shape_bound() {
        let below = regularised_incomplete_beta(0.01, 0.02 * (1.0 - 1e-5), 4.28);
        let at = regularised_incomplete_beta(0.01, 0.02, 4.28);
        assert!(
            below.is_finite() && below > at && below < 1.0,
            "{below} against {at}"
        );
    }

    /// `alpha_from_diversity` sets `α_ref = 1` and splits `θ` across the ALTs:
    /// biallelic → `[1, θ]`, triallelic → `[1, θ/2, θ/2]`.
    #[test]
    fn alpha_from_diversity_sets_ref_one_and_splits_theta() {
        assert_eq!(alpha_from_diversity(2, 1e-3), vec![1.0, 1e-3]);
        let tri = alpha_from_diversity(3, 2e-3);
        assert_eq!(tri[0], 1.0);
        assert!((tri[1] - 1e-3).abs() < 1e-18 && (tri[2] - 1e-3).abs() < 1e-18);
        // Total ALT concentration is θ regardless of allele count.
        assert!((tri[1] + tri[2] - 2e-3).abs() < 1e-15);
    }

    /// A monomorphic shape (`n_alleles = 1`) has no ALT to carry diversity and
    /// returns just `[α_ref]`.
    #[test]
    fn alpha_from_diversity_monomorphic_shape_is_ref_only() {
        assert_eq!(alpha_from_diversity(1, 1e-3), vec![1.0]);
    }

    /// A legitimately tiny — but non-zero — diversity passes through unfloored:
    /// `MIN_ALT_CONCENTRATION = 1e-12` sits far below any real cohort's θ, so a
    /// low-diversity cohort at θ = 1e-8 keeps `α_alt = 1e-8`. Guards the doc
    /// claim that the floor "never perturbs a genuine estimate".
    #[test]
    fn alpha_from_diversity_tiny_real_theta_is_not_floored() {
        let alpha = alpha_from_diversity(2, 1e-8);
        assert_eq!(alpha[1], 1e-8, "a real θ=1e-8 was clamped by the floor");
        assert!(alpha[1] > MIN_ALT_CONCENTRATION);
    }

    /// A nonsensical θ > 1 (reachable via `--diversity`) is passed through as
    /// `[1, θ]` rather than rejected — documents the map has no upper clamp (the
    /// value precondition is only θ ≥ 0).
    #[test]
    fn alpha_from_diversity_theta_above_one_passes_through() {
        assert_eq!(alpha_from_diversity(2, 3.0), vec![1.0, 3.0]);
    }

    /// A zero-allele shape is impossible and trips the hard assertion rather than
    /// silently returning a wrong-length `α` (which would then mis-shape the
    /// Dirichlet-multinomial).
    #[test]
    #[should_panic(expected = "n_alleles must be >= 1")]
    fn alpha_from_diversity_panics_on_zero_alleles() {
        alpha_from_diversity(0, 1e-3);
    }

    /// A degenerate (`n <= 1`) grid returns the interval midpoint rather than
    /// dividing by `n − 1 = 0`.
    #[test]
    fn linear_grid_point_degenerate_grid_returns_midpoint() {
        assert_eq!(linear_grid_point(0, 1, 0.2, 0.8), 0.5);
        assert_eq!(linear_grid_point(0, 0, 0.2, 0.8), 0.5);
        assert_eq!(sfs_grid_point(0, 1, 0.01), 0.5);
    }

    /// `lgamma` reproduces the log-factorial identity `ln Γ(n+1) = ln n!` at
    /// integer arguments and the classic `ln Γ(1/2) = ln √π` at a non-integer
    /// one — the case the integer factorial helpers cannot cover.
    #[test]
    fn lgamma_matches_known_values() {
        // ln Γ(1) = ln 0! = 0, ln Γ(2) = ln 1! = 0.
        assert!(lgamma(1.0).abs() < 1e-12);
        assert!(lgamma(2.0).abs() < 1e-12);
        // ln Γ(n+1) = ln n! for a few n.
        for (n, fact) in [(3u32, 6.0_f64), (5, 120.0), (6, 720.0)] {
            let got = lgamma(n as f64 + 1.0);
            assert!(
                (got - float::ln(fact)).abs() < 1e-10,
                "lgamma({}) = {got}, want ln {fact}",
                n + 1
            );
        }
        // Half-integer absolute anchors — closed forms with √π, so a shared
        // systematic error (which the relative recurrence test cannot see)
        // would show up here.
        let half_ln_pi = float::ln(std::f64::consts::PI) / 2.0;
        // ln Γ(1/2) = ln √π.
        assert!((lgamma(0.5) - half_ln_pi).abs() < 1e-12, "Γ(1/2)");
        // ln Γ(3/2) = ½ln π − ln 2.
        assert!(
            (lgamma(1.5) - (half_ln_pi - float::ln(2.0))).abs() < 1e-12,
            "Γ(3/2) = {}",
            lgamma(1.5)
        );
        // ln Γ(5/2) = ln(3/4) + ½ln π.
        assert!(
            (lgamma(2.5) - (float::ln(3.0 / 4.0) + half_ln_pi)).abs() < 1e-12,
            "Γ(5/2) = {}",
            lgamma(2.5)
        );
    }

    /// The `x > 0` contract is enforced in debug builds: a non-positive argument
    /// trips the debug assertion rather than silently returning a corrupt value.
    #[test]
    #[should_panic(expected = "lgamma requires x > 0")]
    fn lgamma_panics_on_non_positive_in_debug() {
        lgamma(0.0);
    }

    /// `lgamma` satisfies the recurrence `ln Γ(x+1) = ln x + ln Γ(x)` at the
    /// small non-integer arguments the prior actually evaluates (`α ≈ θ`, small).
    #[test]
    fn lgamma_satisfies_recurrence_at_small_args() {
        for &x in &[1e-3, 0.01, 0.3, 1.7] {
            let lhs = lgamma(x + 1.0);
            let rhs = float::ln(x) + lgamma(x);
            assert!((lhs - rhs).abs() < 1e-10, "x={x}: {lhs} vs {rhs}");
        }
    }

    /// The Wright genotype priors are a proper distribution: `hom-ref + het +
    /// hom-alt = 1` for every `(p, F)` — a coefficient typo would break this.
    #[test]
    fn wright_genotype_priors_sum_to_one() {
        for &p in &[0.01, 0.2, 0.5, 0.9] {
            for &f in &[0.0, 0.3, 0.99] {
                let (a, b, c) = wright_genotype_log_priors(p, f);
                let sum = float::exp(a) + float::exp(b) + float::exp(c);
                assert!((sum - 1.0).abs() < 1e-12, "p={p} F={f} sum={sum}");
            }
        }
    }
}
