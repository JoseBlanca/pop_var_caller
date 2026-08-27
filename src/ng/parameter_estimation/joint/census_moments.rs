//! **The population's two numbers, averaged over the census positions the fit walked** — how
//! often a chromosome drawn at random carries something other than the reference base, and how
//! often two chromosomes drawn at random differ.
//!
//! These are the two numbers the SNP/indel genotype prior is built from
//! (`doc/devel/ng/spec/calling_priors.md` §4). The caller can also read both off the *fitted
//! curve*, in closed form and with no census pass
//! (`FrequencyDensity::expected_alternative_frequency` and `::expected_heterozygosity`), and that
//! is what it does today. **What averaging the census positions buys instead is an error that goes
//! to zero as the cohort grows**: a curve converges on the best-fitting member of its family, and
//! a population outside that family stays a few per cent away however much data arrives. Measured
//! at 63 individuals on a population with two frequency peaks, the census average returns 1.000 of
//! the panel's own genotypes at every depth where the integrated curve settles 5 to 7% off
//! (`doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md` §9.1).
//!
//! Design: `doc/devel/ng/spec/ordinary_site_prior_moments.md` §3.
//!
//! **Nothing calls this yet.** Wiring it to a run is that spec's §5 and the implementation plan's
//! Milestone C; what is here is the reduction and the type it returns.

/// **How many alternative copies a position carries across the panel, and how uncertain that
/// count is** — one position's contribution, before it is averaged with the others.
///
/// `expected_copies` is `E[k]` and `copy_count_variance` is the sum of the samples' own posterior
/// variances. They are separated because the two moments need them differently: the mean frequency
/// is linear in `k` and uses only the first, while the heterozygosity is quadratic and needs both
/// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §3.1).
#[derive(Copy, Clone, PartialEq, Debug)]
struct AlternativeCopiesAtAPosition {
    expected_copies: f64,
    copy_count_variance: f64,
}

/// **The two moments, averaged over the census positions.**
///
/// Both describe the population rather than the panel: the finite-panel corrections are already
/// applied, so a run of ten individuals and a run of a thousand are estimating the same two
/// numbers.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct CensusMoments {
    /// **How often a chromosome drawn at random from the population carries something other than
    /// the reference base** — the mean over census positions of `k / 2N`, for `k` alternative
    /// copies among the panel's `2N` chromosomes.
    ///
    /// **The estimator is linear in `k`, so substituting `k`'s expected count is exact.** That is
    /// what makes this the easy one of the two: no curvature term, no correction, and unbiased at
    /// every panel size measured, from one individual to a thousand
    /// (`doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md` §3.1).
    pub mean_alternative_frequency: f64,
    /// **How often two chromosomes drawn at random from the population differ** — the mean over
    /// census positions of `2 k (2N − k) / (2N (2N − 1))`.
    ///
    /// **The `2N − 1` rather than `2N` is what makes this a property of the population rather than
    /// of the panel.** Dropping it returns the heterozygosity `1/(2N)` low, which is 50% at one
    /// individual (`doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md` §10).
    ///
    /// **⚠ Two things this number still owes, and both are named in the plan rather than hidden.**
    /// It substitutes `E[k]` into a formula that is quadratic in `k`, which returns it **high** by
    /// the variance — at one sample and three reads a position, measured, 2.538 ± 0.165 times the
    /// truth rather than 1.219 ± 0.152 (report §4.1). And it applies no inbreeding correction,
    /// which is a further 80% at one individual with an inbreeding coefficient of 0.8 (spec §4).
    /// Plan steps B2 and B3 are those two terms.
    pub heterozygosity: f64,
}

impl CensusMoments {
    /// **Average the two moments over every census position the fit kept.**
    ///
    /// `genotype_posterior` is `JointFit::genotype_posterior`: three numbers a sample a position,
    /// in position order — the posterior that the sample is heterozygous there, that both its
    /// copies are non-reference, and that it carries an extra copy of the position. **The third is
    /// not a genotype and takes no part here**: a sample carrying more copies of a position than
    /// the reference does is a mapping fact, not an allele count, and the fit scores it as its own
    /// class precisely so that it does not have to be read as a heterozygote.
    ///
    /// ## Why `k` is an expectation and never a count
    ///
    /// **At three reads a position a heterozygote often shows only one of its two alleles, and a
    /// sequencing error often looks like a third.** The first hides variation and the second
    /// invents it, and they do not cancel. So the alternative-copy count at a position is taken
    /// as an expectation under the read model, from the fit's own converged per-position
    /// posteriors, rather than counted off called genotypes
    /// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §3.1).
    ///
    /// ## Every kept position counts, unweighted
    ///
    /// **A position that looks mismapped is not weighted down**, and that is a decision rather
    /// than an omission. The fit already scores every position under both the ordinary and the
    /// mismapped error rate and weights the two by their posteriors, so a position whose reads
    /// look mismapped has already had its genotype posteriors computed under an error rate that
    /// explains them without a heterozygote. Weighting a second time removes real variation as
    /// well as artefact: measured over eighteen cells, the weighting is worse in thirteen, better
    /// in four, unchanged in one, and every difference is inside its own error bar except one
    /// loss — at one sample and three reads it moves the estimate from 0.842 of the truth to
    /// **0.674** (spec §3.2, report §6).
    ///
    /// # Panics
    ///
    /// **On a posterior array that is not `positions × samples × 3` long**, held in release. The
    /// reduction reads it by computed offset, so a length disagreement does not fail — it reads
    /// one position's numbers as another's and returns a plausible answer. It also refuses a run
    /// of no samples, which has no chromosomes for a frequency to be a share of.
    #[must_use]
    pub fn from_posteriors(genotype_posterior: &[f32], samples: usize, positions: usize) -> Self {
        assert!(
            samples > 0,
            "a census average over no samples has no chromosomes to be a share of"
        );
        assert_eq!(
            genotype_posterior.len(),
            positions * samples * 3,
            "the fit writes three numbers a sample a position, so {positions} positions over \
             {samples} samples is {} values; the array holds {}",
            positions * samples * 3,
            genotype_posterior.len()
        );
        let chromosomes = 2.0 * samples as f64;
        let mut frequency = 0.0_f64;
        let mut heterozygosity = 0.0_f64;
        for position in 0..positions {
            let copies = alternative_copies_at(genotype_posterior, samples, position);
            frequency += copies.expected_copies / chromosomes;
            heterozygosity += nei_heterozygosity(copies.expected_copies, chromosomes);
        }
        // A census of no positions leaves both sums at zero and the divide undefined. It is a run
        // whose fit kept nothing, which the caller's own fallback ladder handles by regime; what
        // this must not do is return a `NaN` that looks like a number.
        let positions = positions.max(1) as f64;
        Self {
            mean_alternative_frequency: frequency / positions,
            heterozygosity: heterozygosity / positions,
        }
    }
}

/// One position's expected alternative-copy count across the panel, and the sum of the samples'
/// own posterior variances.
///
/// ```text
/// E[k]    =  Σ over samples of   P(het)  +  2 · P(both copies non-reference)
/// Var(k) ≈=  Σ over samples of ( P(het) + 4 · P(both)  −  [P(het) + 2 · P(both)]² )
/// ```
///
/// **The square is inside the sum, per sample.** Squaring the whole sum instead is a different and
/// much larger number, and it is the first thing to check when the heterozygosity comes back
/// wrong (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §3.1).
///
/// **⚠ This variance is not the whole variance, and the shortfall has a known sign.** The exact
/// quantity is `Σᵢ Varᵢ + ΣΣ_{i≠j} Covᵢⱼ`; the samples at a position are coupled through the
/// frequency they share, so the covariance is **positive** and this sum is an under-estimate —
/// which makes the heterozygosity come back slightly high. **Its size has not been measured**; the
/// whole variance term is between 1.6 and 2.2 parts in a hundred at ten samples and three reads,
/// which bounds the residual above by that and no lower (spec §3.1, §8's first open question).
fn alternative_copies_at(
    genotype_posterior: &[f32],
    samples: usize,
    position: usize,
) -> AlternativeCopiesAtAPosition {
    let base = position * samples * 3;
    let mut expected_copies = 0.0_f64;
    let mut copy_count_variance = 0.0_f64;
    for sample in 0..samples {
        let heterozygous = f64::from(genotype_posterior[base + sample * 3]);
        let both_non_reference = f64::from(genotype_posterior[base + sample * 3 + 1]);
        let copies = heterozygous + 2.0 * both_non_reference;
        let copies_squared = heterozygous + 4.0 * both_non_reference;
        expected_copies += copies;
        copy_count_variance += copies_squared - copies * copies;
    }
    AlternativeCopiesAtAPosition {
        expected_copies,
        copy_count_variance,
    }
}

/// **How often two of the panel's chromosomes drawn at random differ at one position** — Nei's
/// gene diversity with the finite-panel correction, `2 k (2N − k) / (2N (2N − 1))`.
///
/// The `2N − 1` in the denominator is what makes the answer a property of the population rather
/// than of the panel: it is the count of *other* chromosomes the first draw can be paired with.
/// Writing `2N` there returns the heterozygosity `1/(2N)` low, which is 50% at one individual and
/// 0.8% at 63.
fn nei_heterozygosity(expected_copies: f64, chromosomes: f64) -> f64 {
    debug_assert!(
        chromosomes >= 2.0,
        "two chromosomes are needed before a pair of them can differ; got {chromosomes}"
    );
    2.0 * expected_copies * (chromosomes - expected_copies) / (chromosomes * (chromosomes - 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Posteriors that are point masses on known genotypes: `copies` alternative copies for each
    /// sample at each position, laid out the way the fit writes them.
    ///
    /// `copies_per_sample` is one entry a sample a position, in position order — so a caller can
    /// build a census where different positions carry different counts, which a single number
    /// could not.
    fn point_masses(copies_per_sample: &[&[u8]]) -> Vec<f32> {
        let mut out = Vec::new();
        for position in copies_per_sample {
            for &copies in position.iter() {
                let (het, both) = match copies {
                    0 => (0.0, 0.0),
                    1 => (1.0, 0.0),
                    2 => (0.0, 1.0),
                    other => panic!("a diploid carries 0, 1 or 2 alternative copies, not {other}"),
                };
                // The third number is the carrier posterior, which takes no part.
                out.extend_from_slice(&[het, both, 0.0]);
            }
        }
        out
    }

    /// **Both moments return the census's own exact values from point-mass posteriors**, at one
    /// individual and at a thousand — `doc/devel/ng/spec/ordinary_site_prior_moments.md` §9's
    /// first test, minus the inbreeding half, which arrives at plan step B3.
    ///
    /// **The two panel sizes are not a formality.** At a thousand individuals a `2N` written where
    /// `2N − 1` belongs is a 0.05% error and sits inside any tolerance a test would set; at one
    /// individual it is 50%. Only the small panel sees it, and only the large one would see a
    /// correction that overshoots.
    #[test]
    fn point_mass_posteriors_return_the_census_s_own_moments() {
        // One individual, four positions: reference, heterozygous, both non-reference,
        // heterozygous. Alternative copies 0, 1, 2, 1 out of 2 chromosomes.
        let one = point_masses(&[&[0], &[1], &[2], &[1]]);
        let moments = CensusMoments::from_posteriors(&one, 1, 4);
        // Mean of 0/2, 1/2, 2/2, 1/2.
        assert!(
            (moments.mean_alternative_frequency - 0.5).abs() < 1e-15,
            "got {}",
            moments.mean_alternative_frequency
        );
        // 2k(2N−k)/(2N(2N−1)) at 2N = 2 is k(2−k): 0, 1, 0, 1 — a mean of a half.
        assert!(
            (moments.heterozygosity - 0.5).abs() < 1e-15,
            "got {}",
            moments.heterozygosity
        );

        // A thousand individuals at one position: one heterozygote and 999 reference samples, so
        // one alternative copy among 2,000 chromosomes.
        let mut thousand_at_one_position = vec![0u8; 1_000];
        thousand_at_one_position[0] = 1;
        let thousand = point_masses(&[&thousand_at_one_position]);
        let moments = CensusMoments::from_posteriors(&thousand, 1_000, 1);
        assert!(
            (moments.mean_alternative_frequency - 1.0 / 2_000.0).abs() < 1e-15,
            "got {}",
            moments.mean_alternative_frequency
        );
        // 2 · 1 · 1999 / (2000 · 1999) = 1/1000. Writing 2N for 2N − 1 gives 3998/4_000_000,
        // which is 0.9995 in 1,000 — 0.05% out, and this tolerance is tighter than that.
        let truth = 1.0 / 1_000.0;
        assert!(
            (moments.heterozygosity / truth - 1.0).abs() < 1e-5,
            "got {} against {truth}",
            moments.heterozygosity
        );
    }

    /// **A census where every sample carries every copy is a heterozygosity of zero, not of one.**
    ///
    /// Two chromosomes drawn from a panel that is fixed for the alternative allele never differ.
    /// This is the mirror of an all-reference census and it is the one a formula that forgot the
    /// `(2N − k)` factor would get wrong: that formula returns the frequency instead, which is 1.
    #[test]
    fn a_panel_fixed_for_either_allele_has_no_heterozygosity() {
        let all_alternative = point_masses(&[&[2, 2, 2], &[2, 2, 2]]);
        let moments = CensusMoments::from_posteriors(&all_alternative, 3, 2);
        assert_eq!(moments.mean_alternative_frequency, 1.0);
        assert_eq!(moments.heterozygosity, 0.0);

        let all_reference = point_masses(&[&[0, 0, 0], &[0, 0, 0]]);
        let moments = CensusMoments::from_posteriors(&all_reference, 3, 2);
        assert_eq!(moments.mean_alternative_frequency, 0.0);
        assert_eq!(moments.heterozygosity, 0.0);
    }

    /// **The two moments are different questions and this fixture tells them apart.**
    ///
    /// At a panel of two individuals with one alternative copy among the four chromosomes, the
    /// frequency is 1 in 4 and the heterozygosity is 2·1·3/(4·3) = 1 in 2 — so a reduction that
    /// returned one where the other belongs is visible. **No fixture in which the two agree can
    /// say that**, and at 2N = 2 they do agree, which is why this panel is two individuals and
    /// not one.
    #[test]
    fn the_frequency_and_the_heterozygosity_are_not_the_same_number() {
        let one_copy_in_four = point_masses(&[&[1, 0]]);
        let moments = CensusMoments::from_posteriors(&one_copy_in_four, 2, 1);
        assert!((moments.mean_alternative_frequency - 0.25).abs() < 1e-15);
        assert!((moments.heterozygosity - 0.5).abs() < 1e-15);
    }

    /// **The carrier posterior — the fit's third number a sample — takes no part.**
    ///
    /// A sample carrying more copies of a position than the reference does is a mapping fact
    /// rather than an allele count, and the fit scores it as its own class so that it need not be
    /// read as a heterozygote. Setting it to one at every sample must move neither moment.
    #[test]
    fn the_carrier_posterior_moves_neither_moment() {
        let plain = point_masses(&[&[1, 0], &[2, 1]]);
        let mut with_carriers = plain.clone();
        for sample in 0..4 {
            with_carriers[sample * 3 + 2] = 1.0;
        }
        assert_eq!(
            CensusMoments::from_posteriors(&plain, 2, 2),
            CensusMoments::from_posteriors(&with_carriers, 2, 2)
        );
    }

    /// **The per-sample variance is a sum of per-sample squares, not the square of the sum.**
    ///
    /// Nothing reads it until plan step B2, so this pins the quantity where it is computed. At
    /// posteriors of a half on heterozygous and a half on both-non-reference, one sample's
    /// expected copies are `0.5 + 1.0 = 1.5` and its `E[k²]` is `0.5 + 2.0 = 2.5`, so its variance
    /// is `2.5 − 2.25 = 0.25`. **Two such samples give 0.5**, where squaring the whole sum instead
    /// would give `5.0 − 9.0`, a negative number — which is the check.
    #[test]
    fn the_copy_count_variance_is_summed_per_sample() {
        let midway = vec![0.5_f32, 0.5, 0.0, 0.5, 0.5, 0.0];
        let copies = alternative_copies_at(&midway, 2, 0);
        assert!((copies.expected_copies - 3.0).abs() < 1e-15);
        assert!(
            (copies.copy_count_variance - 0.5).abs() < 1e-15,
            "got {}",
            copies.copy_count_variance
        );
    }

    /// **A point-mass genotype has no variance**, whichever of the three it is — the sanity check
    /// that the variance is about uncertainty and not about the genotype.
    #[test]
    fn a_certain_genotype_has_no_copy_count_variance() {
        for copies in [0u8, 1, 2] {
            let certain = point_masses(&[&[copies]]);
            let at = alternative_copies_at(&certain, 1, 0);
            assert_eq!(
                at.copy_count_variance, 0.0,
                "a sample certain of carrying {copies} copies has no uncertainty about it"
            );
        }
    }

    /// **A posterior array of the wrong length is refused rather than read by offset.**
    ///
    /// The reduction indexes by `position × samples × 3`, so a disagreement between the length and
    /// the two counts does not fail — it reads one position's numbers as another's and returns a
    /// plausible answer.
    #[test]
    #[should_panic(expected = "three numbers a sample a position")]
    fn a_posterior_array_of_the_wrong_length_is_refused() {
        let two_positions = point_masses(&[&[1], &[0]]);
        let _ = CensusMoments::from_posteriors(&two_positions, 1, 3);
    }

    /// **A run of no samples is refused**, because a frequency is a share of the panel's
    /// chromosomes and there are none.
    #[test]
    #[should_panic(expected = "no chromosomes to be a share of")]
    fn a_census_over_no_samples_is_refused() {
        let _ = CensusMoments::from_posteriors(&[], 0, 0);
    }

    /// **A census of no positions returns zeros rather than `NaN`.**
    ///
    /// It is a run whose fit kept nothing, and what must not happen is a division by zero
    /// producing two numbers that look like measurements.
    #[test]
    fn a_census_of_no_positions_returns_zeros() {
        let moments = CensusMoments::from_posteriors(&[], 3, 0);
        assert_eq!(moments.mean_alternative_frequency, 0.0);
        assert_eq!(moments.heterozygosity, 0.0);
    }
}
