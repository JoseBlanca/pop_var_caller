//! The population-genetics arithmetic three parts of the caller share: how likely each
//! genotype is before any read is looked at, and the two numerical primitives that
//! calculation needs.
//!
//! What lives here is what more than one module uses. The genotype priors
//! ([`calling::genotype_prior`](crate::calling::genotype_prior)) build a site's prior from
//! [`alpha_from_diversity`] and [`wright_genotype_log_priors`]; the site quality
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
        hom_ref.max(PROBABILITY_FLOOR).ln(),
        het.max(PROBABILITY_FLOOR).ln(),
        hom_alt.max(PROBABILITY_FLOOR).ln(),
    )
}

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
                (got - fact.ln()).abs() < 1e-10,
                "lgamma({}) = {got}, want ln {fact}",
                n + 1
            );
        }
        // Half-integer absolute anchors — closed forms with √π, so a shared
        // systematic error (which the relative recurrence test cannot see)
        // would show up here.
        let half_ln_pi = std::f64::consts::PI.ln() / 2.0;
        // ln Γ(1/2) = ln √π.
        assert!((lgamma(0.5) - half_ln_pi).abs() < 1e-12, "Γ(1/2)");
        // ln Γ(3/2) = ½ln π − ln 2.
        assert!(
            (lgamma(1.5) - (half_ln_pi - 2.0_f64.ln())).abs() < 1e-12,
            "Γ(3/2) = {}",
            lgamma(1.5)
        );
        // ln Γ(5/2) = ln(3/4) + ½ln π.
        assert!(
            (lgamma(2.5) - ((3.0_f64 / 4.0).ln() + half_ln_pi)).abs() < 1e-12,
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
            let rhs = x.ln() + lgamma(x);
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
                let sum = a.exp() + b.exp() + c.exp();
                assert!((sum - 1.0).abs() < 1e-12, "p={p} F={f} sum={sum}");
            }
        }
    }
}
