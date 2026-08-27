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

use crate::ng::types::InbreedingF;

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
    /// **The probability that this position segregates in this panel** — that the panel's `2N`
    /// chromosomes are not all reference and not all alternative.
    ///
    /// See `probability_that_the_panel_segregates` below for how it is formed and what it assumes.
    segregating: f64,
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
    /// **`k` is an expectation and `k (2N − k)` is quadratic, so the curvature term is carried
    /// rather than dropped**: `E[k(2N − k)] = 2N·E[k] − E[k]² − Var(k)`. Without it this comes
    /// back high by exactly `Var(k)`, which at one sample and three reads a position is most of
    /// the answer — 2.538 ± 0.165 times the truth against 1.219 ± 0.152 with it (report §4.1).
    /// See [`nei_heterozygosity`].
    ///
    /// **And the panel's inbreeding is divided back out**, by `1 − F/(2N − 1)`. A pair of
    /// chromosomes drawn at random from the panel comes from the same individual with probability
    /// `1/(2N − 1)`, and such a pair cannot differ if it is one ancestral copy counted twice — so
    /// an inbred panel shows fewer differences than the population has. See [`inbreeding_factor`],
    /// where the size at each cohort size is: 80% at one individual and 0.04% at a thousand, at
    /// tomato's fitted `F` of 0.8.
    ///
    /// **⚠ What this number still owes is one term, and it is bounded rather than unknown.** The
    /// variance above is the sum of the samples' own and ignores the positive coupling between
    /// samples at a position, so this comes back slightly **high** — see [`alternative_copies_in`].
    /// Nothing has measured the size; the whole variance term is 1.6 to 2.2 parts in a hundred at
    /// ten samples and three reads, which bounds the residual above by that and no lower (spec
    /// §3.1, §8's first open question).
    pub heterozygosity: f64,
    /// **How many census positions the fit walked** — the denominator both moments were averaged
    /// over, not the count of positions that carry any variation.
    pub positions: u64,
    /// **How many of those positions segregate in this panel**, as a **soft count**:
    /// `Σ over positions of P(the position segregates)`.
    ///
    /// **Both moments' precision rests on this rather than on [`Self::positions`]**: a position
    /// where the population carries one allele contributes zero to the heterozygosity and tells
    /// you nothing about the frequency. A two-million-position census on a panel segregating at 1
    /// in 200 carries about ten thousand segregating positions; one run over a small `--regions`
    /// BED may carry a few hundred (spec §6.2).
    ///
    /// **It is a soft count and not a count above a threshold**, and
    /// `probability_that_the_panel_segregates` in this module says why: the obvious hard version
    /// reports 100% of positions.
    ///
    /// **No floor is applied to it anywhere.** Spec §6.2 forbids picking one until it is measured
    /// and the measurement needs a real census, so the run reports the count and takes no action —
    /// which is distinguishable in the output from a floor that never fires.
    pub segregating_positions: f64,
    /// **A lower bound on how far [`Self::mean_alternative_frequency`] would move if the census
    /// were drawn again** — the standard error of the mean across positions.
    ///
    /// **A floor and not an interval.** See [`standard_error_of_the_mean`]: census positions are
    /// linked, so a spread computed as though they were independent is too narrow by a factor
    /// `parameter_prepass_census_sites.md` §5 puts between 3 and 16.
    pub frequency_standard_error_floor: f64,
    /// The same for [`Self::heterozygosity`], and a floor for the same reason. It carries the
    /// inbreeding correction, because the correction is a constant divide and scales the spread
    /// with the number it belongs to.
    pub heterozygosity_standard_error_floor: f64,
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
    /// ## The panel's inbreeding coefficient, and where it belongs
    ///
    /// `panel_inbreeding` is the probability that an individual's two copies of a position are one
    /// ancestral copy counted twice, **averaged over the panel's samples, unweighted** — a sample
    /// with more census positions covered must not count for more, because with a per-individual
    /// `Fᵢ` the estimator's expectation is `π · (1 − F̄/(2N − 1))` with `F̄` the plain mean
    /// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §4).
    ///
    /// **It corrects the heterozygosity and not the frequency**, because the frequency is linear
    /// in the copy count and inbreeding rearranges copies between individuals without changing
    /// how many there are. See [`Self::heterozygosity`] for the size at each end of the cohort
    /// range.
    ///
    /// # Panics
    ///
    /// **On a posterior array that is not `positions × samples × 3` long**, held in release. The
    /// reduction reads it by computed offset, so a length disagreement does not fail — it reads
    /// one position's numbers as another's and returns a plausible answer. It also refuses a run
    /// of no samples, which has no chromosomes for a frequency to be a share of.
    #[must_use]
    pub fn from_posteriors(
        genotype_posterior: &[f32],
        samples: usize,
        positions: usize,
        panel_inbreeding: InbreedingF,
    ) -> Self {
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
        let mut running = CensusMomentSums::over(samples);
        let mut one_position = vec![0.0_f64; samples * 3];
        for position in 0..positions {
            let base = position * samples * 3;
            for (slot, value) in one_position
                .iter_mut()
                .zip(&genotype_posterior[base..base + samples * 3])
            {
                *slot = f64::from(*value);
            }
            running.add_position(&one_position);
        }
        running.finish(panel_inbreeding)
    }
}

/// **The two moments as running sums, so a run never stores the array they are summed from.**
///
/// One of these lives on the fit's per-chunk statistics, is fed one position at a time from the
/// expectation step's own scratch buffer, and is merged chunk into chunk the way every other sum
/// there is. What it holds is a handful of `f64`s, whatever the census and whatever the cohort.
///
/// ## Why this and not the stored array
///
/// The fit can be asked to keep every position's genotype posteriors
/// (`JointFitConfig::genotype_posteriors`), and [`CensusMoments::from_posteriors`] reduces that
/// array. **It is off by default because it weighs twelve bytes a position a sample** — three
/// `f32`s — which is **1.5 GB for fifty samples over the shipped two-million-position census**,
/// held to be summed once and thrown away. The design asks for the sums instead
/// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §5), and this is them.
///
/// **The array route is kept and is not dead**: it is what the measurement harnesses reduce, and
/// it is the oracle these sums are checked against — the two must agree to the last bit on the
/// same posteriors, which is a test.
#[derive(Clone, PartialEq, Debug)]
pub struct CensusMomentSums {
    /// The panel's size in diploid individuals. Held rather than derived from the chromosome count
    /// because [`Self::add_position`] runs once per census position per pass of the fit — on the
    /// shipped census that is two million times a pass — and a stored count is one field read
    /// where a derived one is a divide and a cast.
    samples: usize,
    chromosomes: f64,
    positions: u64,
    frequency: f64,
    heterozygosity: f64,
    /// The sums of squares, which are what a spread across positions needs beside the sums.
    frequency_squares: f64,
    heterozygosity_squares: f64,
    /// `Σ P(the position segregates in this panel)` — the soft count of
    /// [`probability_that_the_panel_segregates`].
    segregating: f64,
}

impl CensusMomentSums {
    /// Start a run's sums over a panel of `samples` diploid individuals.
    ///
    /// # Panics
    ///
    /// On a panel of no samples, which has no chromosomes for a frequency to be a share of.
    #[must_use]
    pub fn over(samples: usize) -> Self {
        assert!(
            samples > 0,
            "a census average over no samples has no chromosomes to be a share of"
        );
        Self {
            samples,
            chromosomes: 2.0 * samples as f64,
            positions: 0,
            frequency: 0.0,
            heterozygosity: 0.0,
            frequency_squares: 0.0,
            heterozygosity_squares: 0.0,
            segregating: 0.0,
        }
    }

    /// Add one census position, from **three numbers a sample in sample order** — the posterior
    /// that the sample is heterozygous there, that both its copies are non-reference, and that it
    /// carries an extra copy of the position.
    ///
    /// That is the expectation step's own per-position scratch buffer, unchanged and uncopied.
    ///
    /// # Panics
    ///
    /// On a buffer that is not three long per sample. It is read by computed offset, so a
    /// disagreement would read one sample's numbers as another's rather than fail.
    pub fn add_position(&mut self, position_genotype: &[f64]) {
        assert_eq!(
            position_genotype.len(),
            self.samples * 3,
            "one position carries three numbers a sample, so {} samples is {} values; the buffer \
             holds {}",
            self.samples,
            self.samples * 3,
            position_genotype.len()
        );
        let copies = alternative_copies_in(position_genotype, self.samples);
        let frequency = copies.expected_copies / self.chromosomes;
        let heterozygosity = nei_heterozygosity(
            copies.expected_copies,
            copies.copy_count_variance,
            self.chromosomes,
        );
        self.positions += 1;
        self.frequency += frequency;
        self.heterozygosity += heterozygosity;
        self.frequency_squares += frequency * frequency;
        self.heterozygosity_squares += heterozygosity * heterozygosity;
        self.segregating += copies.segregating;
    }

    /// Fold another chunk's sums in. Every field is a sum over positions, so this is addition.
    ///
    /// # Panics
    ///
    /// On sums taken over a different panel, which cannot be added: the per-position terms are
    /// already divided by the panel's chromosome count.
    pub fn merge(&mut self, other: &Self) {
        assert_eq!(
            self.samples, other.samples,
            "two chunks of one run share a panel, so their sums are over the same sample count; \
             got {} and {}",
            self.samples, other.samples
        );
        self.positions += other.positions;
        self.frequency += other.frequency;
        self.heterozygosity += other.heterozygosity;
        self.frequency_squares += other.frequency_squares;
        self.heterozygosity_squares += other.heterozygosity_squares;
        self.segregating += other.segregating;
    }

    /// How many census positions have been added.
    #[must_use]
    pub fn positions(&self) -> u64 {
        self.positions
    }

    /// Divide through by the positions and apply the panel's inbreeding correction.
    ///
    /// **A census of no positions returns zeros rather than `NaN`.** It is a run whose fit kept
    /// nothing; the caller's own fallback ladder handles that by regime, and what this must not do
    /// is return two numbers that look like measurements.
    #[must_use]
    pub fn finish(&self, panel_inbreeding: InbreedingF) -> CensusMoments {
        let positions = self.positions.max(1) as f64;
        let correction = inbreeding_factor(panel_inbreeding, self.chromosomes);
        CensusMoments {
            mean_alternative_frequency: self.frequency / positions,
            heterozygosity: self.heterozygosity / positions / correction,
            positions: self.positions,
            segregating_positions: self.segregating,
            // **The inbreeding correction is a constant divide, so it scales the spread with the
            // number.** Applying it to one and not the other would make the two describe different
            // quantities.
            frequency_standard_error_floor: standard_error_of_the_mean(
                self.frequency,
                self.frequency_squares,
                self.positions,
            ),
            heterozygosity_standard_error_floor: standard_error_of_the_mean(
                self.heterozygosity,
                self.heterozygosity_squares,
                self.positions,
            ) / correction,
        }
    }
}

/// **How far the mean over positions would move if the census were drawn again — if the positions
/// were independent, which they are not.**
///
/// The plain standard error of the mean: the spread across positions divided by the square root of
/// how many there are, using the one-pass form `Σx² / n − (Σx / n)²` for the variance.
///
/// **What it is for is being labelled a floor.** Census positions are scattered across the genome
/// but they are linked — a run of homozygosity or a shared haplotype makes neighbours carry the
/// same information twice — so a spread computed as though they were independent counts the same
/// evidence more than once and comes back **too narrow**.
/// `parameter_prepass_census_sites.md` §5 puts that factor **between 3 and 16**. So the honest
/// reading is *the true spread is at least this*, never *the answer is within this*.
///
/// **Below two positions it is zero**, because one position has no spread to measure and the
/// `n − 1` a sample variance divides by is zero there. That is a run whose fit kept almost
/// nothing, and the segregating count beside it is what says so.
///
/// The variance is clamped at zero: the one-pass form can return a tiny negative from
/// cancellation when every position carries nearly the same value, which is exactly what a census
/// of a nearly invariant cohort is.
fn standard_error_of_the_mean(sum: f64, sum_of_squares: f64, count: u64) -> f64 {
    if count < 2 {
        return 0.0;
    }
    let n = count as f64;
    let mean = sum / n;
    let variance = (sum_of_squares / n - mean * mean).max(0.0) * n / (n - 1.0);
    (variance / n).sqrt()
}

/// **What the panel's inbreeding takes off the heterozygosity before the correction puts it
/// back** — `1 − F/(2N − 1)`.
///
/// A pair of chromosomes drawn at random from the panel is drawn from the *same individual* with
/// probability `1/(2N − 1)`, and with probability `F` such a pair is one ancestral copy counted
/// twice and cannot differ. So an inbred panel shows fewer differences than the population has,
/// by exactly this factor, and dividing by it recovers the population's own
/// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §3, §4).
///
/// **The size at the two ends of the committed cohort range is not the same order of thing**, at
/// tomato's fitted `F` of 0.8:
///
/// ```text
///   individuals      1      2      3     10     63   1000
///   the shortfall  80%    27%    16%     4%   0.6%  0.04%
/// ```
///
/// Measured across four populations at nine panel sizes: of 36 cells, 21 sit within one standard
/// error of the value this factor predicts, 33 within two and all 36 within three — which is what
/// a correct formula and 36 draws give (report §3.3).
///
/// **It never divides by zero**, because [`InbreedingF`] admits `[0, 1)` and the factor is
/// smallest at one individual, where it is `1 − F`. **At one individual and an `F` near one it is
/// near zero, and the corrected heterozygosity is correspondingly enormous — which is the factor
/// meaning what it says**: a fully inbred individual shows no heterozygotes at all, so no observed
/// rate can identify the population's diversity. `ordinary_site_prior_moments.md` §4.2 records that
/// nothing reading positions one at a time can estimate `F` from a single genome either, which is
/// why a user is given the coefficient to override.
fn inbreeding_factor(panel_inbreeding: InbreedingF, chromosomes: f64) -> f64 {
    debug_assert!(
        chromosomes >= 2.0,
        "a pair of chromosomes is drawn from the same individual with probability 1/(2N − 1), \
         which needs at least two; got {chromosomes}"
    );
    1.0 - panel_inbreeding.get() / (chromosomes - 1.0)
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
fn alternative_copies_in(
    position_genotype: &[f64],
    samples: usize,
) -> AlternativeCopiesAtAPosition {
    let mut expected_copies = 0.0_f64;
    let mut copy_count_variance = 0.0_f64;
    let mut all_reference = 1.0_f64;
    let mut all_alternative = 1.0_f64;
    for sample in 0..samples {
        let heterozygous = position_genotype[sample * 3];
        let both_non_reference = position_genotype[sample * 3 + 1];
        let copies = heterozygous + 2.0 * both_non_reference;
        let copies_squared = heterozygous + 4.0 * both_non_reference;
        expected_copies += copies;
        copy_count_variance += copies_squared - copies * copies;
        all_reference *= (1.0 - heterozygous - both_non_reference).max(0.0);
        all_alternative *= both_non_reference;
    }
    AlternativeCopiesAtAPosition {
        expected_copies,
        copy_count_variance,
        segregating: (1.0 - all_reference - all_alternative).clamp(0.0, 1.0),
    }
}

/// **The probability that a position segregates in this panel** — that its `2N` chromosomes are
/// neither all reference nor all alternative.
///
/// ```text
/// P(segregates)  =  1  −  Π over samples of P(no alternative copy)
///                      −  Π over samples of P(both copies non-reference)
/// ```
///
/// ## Why it has to be a probability and not a count
///
/// **`doc/devel/ng/spec/ordinary_site_prior_moments.md` §6.2 records an earlier draft of that spec
/// getting this wrong**, and the shape of the error is worth keeping: it defined the segregating
/// positions as *those whose expected alternative-copy count is above zero*. Posteriors from reads
/// are continuous, so at a census of two million positions essentially every one of them has an
/// expected count above zero and the run would report **100% of positions segregating** — a number
/// that looks like a measurement and is a property of arithmetic.
///
/// A soft count needs no threshold constant and degrades smoothly where the reads are thin, where
/// a hard one steps.
///
/// ## What it assumes
///
/// **That the samples' genotypes are independent given the position's posteriors, which they are
/// not.** They are coupled through the frequency they share, exactly as `Var(k)` is
/// (`ordinary_site_prior_moments.md` §3.1). So this is an approximation, and the direction is
/// knowable: positive coupling makes *all-reference* and *all-alternative* more likely than the
/// products above, so the true probability of segregating is a little **lower** and this count runs
/// a little high. **Its size has not been measured and this claims none.**
///
/// **The carrier posterior takes no part**, for the same reason it takes no part in `E[k]`: a
/// sample carrying an extra copy of the position contributes no alternative copies to `k`, so a
/// position where every sample is a carrier reads as all-reference here — consistent with the
/// count it is a count of.
///
/// **What the count is for**: a run over a small `--regions` BED may carry a few hundred
/// segregating positions where the shipped census carries about ten thousand, and both moments'
/// precision rests on that number rather than on how many positions were walked. **No floor is
/// applied to it** — spec §6.2 forbids picking one until it is measured, and the measurement needs
/// a real census. The run reports it and takes no action.
///
/// **The run's own path does not call this**, because [`alternative_copies_in`] forms the same
/// number in the loop it is already making over the samples and a second pass would double the
/// hot path's work. This is that field under the name it deserves, so the argument above has
/// somewhere to live and the tests have something to call.
#[cfg(test)]
fn probability_that_the_panel_segregates(position_genotype: &[f64], samples: usize) -> f64 {
    alternative_copies_in(position_genotype, samples).segregating
}

/// **How often two of the panel's chromosomes drawn at random differ at one position** — Nei's
/// gene diversity with the finite-panel correction, `2 k (2N − k) / (2N (2N − 1))`.
///
/// The `2N − 1` in the denominator is what makes the answer a property of the population rather
/// than of the panel: it is the count of *other* chromosomes the first draw can be paired with.
/// Writing `2N` there returns the heterozygosity `1/(2N)` low, which is 50% at one individual and
/// 0.8% at 63.
///
/// ## Why the variance is an argument, and what happens without it
///
/// **`k` is never counted — it is an expectation under the read model — and `k (2N − k)` is
/// quadratic, so substituting `E[k]` is not the same as taking the expectation:**
///
/// ```text
///   E[k (2N − k)]  =  2N · E[k]  −  E[k]²  −  Var(k)
/// ```
///
/// **Leaving `Var(k)` out returns the heterozygosity high by exactly it**, and at low depth that
/// is not a correction but most of the answer. Measured on drawn cohorts read through the fit, at
/// one sample and three reads a position: **2.538 ± 0.165 times the truth without the term and
/// 1.219 ± 0.152 with it** (`doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md`
/// §4.1).
///
/// **A cohort test cannot see this**, which is why the fixture that pins it runs at one sample.
/// `Var(k)` grows with the panel and `E[k]²` grows with its square, so the term's share of the
/// answer falls like `1/N`: at 63 samples of the same per-sample uncertainty it is under 1%, and a
/// real 63-sample fit is more certain per sample than that, which is where the report's *three
/// decimals* comes from.
///
/// **The variance passed in is the sum of the samples' own and is an under-estimate** — see
/// [`alternative_copies_in`] — so what comes back is still slightly high, in the same direction
/// and much smaller.
fn nei_heterozygosity(expected_copies: f64, copy_count_variance: f64, chromosomes: f64) -> f64 {
    debug_assert!(
        chromosomes >= 2.0,
        "two chromosomes are needed before a pair of them can differ; got {chromosomes}"
    );
    2.0 * (chromosomes * expected_copies - expected_copies * expected_copies - copy_count_variance)
        / (chromosomes * (chromosomes - 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A panel of unrelated individuals: an inbreeding coefficient of zero, where the correction
    /// is exactly one and every test above is measuring something else.
    fn outbred() -> InbreedingF {
        InbreedingF::try_new(0.0).expect("zero is a legal coefficient")
    }

    /// **The inbreeding correction at one individual, where it is `1 − F` and its absence is an
    /// 80% error** — `doc/devel/ng/spec/ordinary_site_prior_moments.md` §9's first test, its
    /// inbreeding half.
    ///
    /// At one individual there is one pair of chromosomes to draw and it is inside that
    /// individual, so `1/(2N − 1)` is exactly 1 and the factor collapses to `1 − F`. At `F = 0.8`
    /// that is 0.2: **a run without the correction reports a fifth of the population's diversity**,
    /// and the panel it measured genuinely is that much less variable than the population it came
    /// from.
    ///
    /// The fixture is a single genome heterozygous at one census position in four. Its observed
    /// heterozygosity is 0.25; the population's is 1.25.
    #[test]
    fn at_one_individual_the_correction_is_one_minus_f() {
        let genome = point_masses(&[&[1], &[0], &[0], &[0]]);
        let selfing = InbreedingF::try_new(0.8).expect("tomato's fitted range");

        let observed = CensusMoments::from_posteriors(&genome, 1, 4, outbred()).heterozygosity;
        assert!(
            (observed - 0.25).abs() < 1e-15,
            "one heterozygous position in four is an observed heterozygosity of 0.25; got \
             {observed}"
        );

        let corrected = CensusMoments::from_posteriors(&genome, 1, 4, selfing).heterozygosity;
        assert!(
            (corrected - 1.25).abs() < 1e-15,
            "at F = 0.8 and one individual the factor is 1 − F = 0.2, so 0.25 becomes 1.25; got \
             {corrected}"
        );
        assert!(
            (observed / corrected - 0.2).abs() < 1e-15,
            "without the correction a run reports a fifth of the population's diversity; got {}",
            observed / corrected
        );
    }

    /// **The same coefficient at a thousand individuals moves the answer by 4 parts in 10,000**,
    /// which is the other end of the committed cohort range and the reason the fixture above runs
    /// at one.
    ///
    /// `1 − 0.8/1999` is 0.99960, so the panel shows 4.002 parts in 10,000 less than the population
    /// and putting it back lifts the answer by 4.004 parts in 10,000 — a rounding term at a large
    /// panel and most of the answer at a single genome. **A test written here would pass with the correction
    /// deleted** at any tolerance loose enough to survive the census's own scatter.
    ///
    /// It also pins the `2N − 1`: writing `2N` gives `1 − 0.8/2000`, which differs from the truth
    /// by 2 parts in ten million — invisible at this panel, and the reason the *shape* of the
    /// denominator is checked at one individual instead, where `2N − 1 = 1` and `2N = 2` are a
    /// factor of two apart.
    #[test]
    fn at_a_thousand_individuals_the_correction_is_four_parts_in_ten_thousand() {
        let mut one_heterozygote = vec![0u8; 1_000];
        one_heterozygote[0] = 1;
        let panel = point_masses(&[&one_heterozygote]);
        let selfing = InbreedingF::try_new(0.8).expect("tomato's fitted range");

        let observed = CensusMoments::from_posteriors(&panel, 1_000, 1, outbred()).heterozygosity;
        let corrected = CensusMoments::from_posteriors(&panel, 1_000, 1, selfing).heterozygosity;
        // **The shortfall and the lift are the same fact from the two sides, and they are not the
        // same number.** The panel shows `1 − F/(2N − 1)` of the population's diversity — a
        // shortfall of `0.8/1999`, 4.002 parts in 10,000, which is spec §4's table entry. Putting
        // it back divides by that factor, so the lift is `(0.8/1999)/(1 − 0.8/1999)`: 4.004 parts
        // in 10,000. Asserting the first where the second belongs is the slip this comment exists
        // to stop.
        let shortfall = 0.8 / 1_999.0;
        let lift = corrected / observed - 1.0;
        assert!(
            (lift - shortfall / (1.0 - shortfall)).abs() < 1e-15,
            "at a thousand individuals the factor is 1 − 0.8/1999, so the lift is 4.004 parts in \
             10,000; got {lift}"
        );
        assert!(
            (shortfall - 4.002e-4).abs() < 1e-7 && (lift - 4.004e-4).abs() < 1e-7,
            "shortfall {shortfall}, lift {lift}"
        );
        // 4 parts in 10,000 against the 80% at one individual — a four-order-of-magnitude range
        // across the cohort sizes this caller commits to.
        assert!(lift < 1e-3);
    }

    /// **The correction moves the heterozygosity and never the mean frequency.**
    ///
    /// Inbreeding rearranges copies between an individual's two chromosomes; it does not change
    /// how many alternative copies the panel holds. The frequency is linear in that count, so it
    /// is untouched — and a correction applied to both would be wrong in a way no ratio test would
    /// see, because both would move together.
    #[test]
    fn the_correction_leaves_the_mean_frequency_alone() {
        let panel = point_masses(&[&[1, 0, 2], &[0, 1, 1]]);
        let selfing = InbreedingF::try_new(0.8).expect("tomato's fitted range");
        let plain = CensusMoments::from_posteriors(&panel, 3, 2, outbred());
        let corrected = CensusMoments::from_posteriors(&panel, 3, 2, selfing);
        assert_eq!(
            plain.mean_alternative_frequency,
            corrected.mean_alternative_frequency
        );
        assert!(corrected.heterozygosity > plain.heterozygosity);
    }
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
        let moments = CensusMoments::from_posteriors(&one, 1, 4, outbred());
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
        let moments = CensusMoments::from_posteriors(&thousand, 1_000, 1, outbred());
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
        let moments = CensusMoments::from_posteriors(&all_alternative, 3, 2, outbred());
        assert_eq!(moments.mean_alternative_frequency, 1.0);
        assert_eq!(moments.heterozygosity, 0.0);

        let all_reference = point_masses(&[&[0, 0, 0], &[0, 0, 0]]);
        let moments = CensusMoments::from_posteriors(&all_reference, 3, 2, outbred());
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
        let moments = CensusMoments::from_posteriors(&one_copy_in_four, 2, 1, outbred());
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
            CensusMoments::from_posteriors(&plain, 2, 2, outbred()),
            CensusMoments::from_posteriors(&with_carriers, 2, 2, outbred())
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
        let midway = vec![0.5_f64, 0.5, 0.0, 0.5, 0.5, 0.0];
        let copies = alternative_copies_in(&midway, 2);
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
            let certain = as_f64(&point_masses(&[&[copies]]));
            let at = alternative_copies_in(&certain, 1);
            assert_eq!(
                at.copy_count_variance, 0.0,
                "a sample certain of carrying {copies} copies has no uncertainty about it"
            );
        }
    }

    /// Posteriors that are **midway between genotypes** — the shape low depth produces, where the
    /// reads have not decided which of the three a sample carries.
    ///
    /// `(reference, heterozygous, both non-reference)` per sample, the same three at every sample
    /// and every position.
    fn midway(
        reference: f64,
        heterozygous: f64,
        both: f64,
        samples: usize,
        positions: usize,
    ) -> Vec<f32> {
        assert!(
            (reference + heterozygous + both - 1.0).abs() < 1e-12,
            "a sample's three genotype posteriors are a distribution"
        );
        let mut out = Vec::new();
        for _ in 0..positions * samples {
            // The third number is the carrier posterior, which takes no part.
            out.extend_from_slice(&[heterozygous as f32, both as f32, 0.0]);
        }
        out
    }

    /// The heterozygosity this reduction would return if it substituted `E[k]` and stopped there —
    /// **the defect plan step B2 exists to prevent**, written out so a test can measure the gap
    /// rather than assert it.
    fn without_the_variance_term(
        genotype_posterior: &[f32],
        samples: usize,
        positions: usize,
    ) -> f64 {
        let chromosomes = 2.0 * samples as f64;
        let wide = as_f64(genotype_posterior);
        let mut total = 0.0_f64;
        for position in 0..positions {
            let base = position * samples * 3;
            let copies = alternative_copies_in(&wide[base..base + samples * 3], samples);
            total += nei_heterozygosity(copies.expected_copies, 0.0, chromosomes);
        }
        total / positions as f64
    }

    /// The fit writes `f32`; every helper here reads `f64`. Widening is exact.
    fn as_f64(narrow: &[f32]) -> Vec<f64> {
        narrow.iter().map(|value| f64::from(*value)).collect()
    }

    /// **At one individual the heterozygosity is exactly the posterior that it is heterozygous**,
    /// and that is the oracle the variance term has to reproduce.
    ///
    /// It is what the question means: an individual's two chromosomes differ at a position exactly
    /// when it is heterozygous there, so averaging over the posterior gives `P(het)` and nothing
    /// else. **The algebra collapses to it and shares nothing with the formula** — with `h` the
    /// heterozygous posterior and `d` the both-non-reference one, `E[k] = h + 2d`,
    /// `Var(k) = h + 4d − (h + 2d)²`, and
    ///
    /// ```text
    ///   2 E[k] − E[k]² − Var(k)  =  2(h + 2d) − (h + 2d)² − h − 4d + (h + 2d)²  =  h
    /// ```
    ///
    /// so every term but `h` cancels. **Drop `Var(k)` and it does not**: what is left is
    /// `2(h + 2d) − (h + 2d)²`, which is a different number at every `d` above zero.
    #[test]
    fn at_one_individual_the_heterozygosity_is_the_heterozygous_posterior() {
        for (reference, heterozygous, both) in [
            (0.3, 0.4, 0.3),
            (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
            (0.8, 0.15, 0.05),
            (0.05, 0.05, 0.90),
            (0.5, 0.5, 0.0),
        ] {
            let posteriors = midway(reference, heterozygous, both, 1, 1);
            let moments = CensusMoments::from_posteriors(&posteriors, 1, 1, outbred());
            // **Against what the array actually holds, not against the literal.** The fit writes
            // `f32`, so 0.4 arrives as 0.40000000596; comparing to the `f64` literal would need a
            // tolerance a thousand million times looser than the arithmetic deserves, and would
            // then admit a real error of that size.
            let in_the_array = f64::from(heterozygous as f32);
            assert!(
                (moments.heterozygosity - in_the_array).abs() < 1e-15,
                "at posteriors ({reference}, {heterozygous}, {both}) the heterozygosity is {} \
                 where the sample's own heterozygous posterior is {in_the_array}",
                moments.heterozygosity
            );
        }
    }

    /// **The variance term, at one sample — where dropping it returns two and a half times the
    /// truth.** `doc/devel/ng/spec/ordinary_site_prior_moments.md` §9's second test.
    ///
    /// The posteriors are `(0.3, 0.4, 0.3)` over reference, heterozygous and both non-reference:
    /// reads that have barely decided, which is the shape three reads a position produces. **The
    /// truth is 0.400** — the heterozygous posterior, by the test above. Substituting `E[k]` and
    /// stopping gives `2 E[k] − E[k]²`, and `E[k]` is exactly 1 here, so it gives **1.000**:
    /// **2.5 times the truth**, against the 2.538 ± 0.165 the report measures through a whole fit.
    ///
    /// **A cohort fixture cannot catch this, and the second half of this test shows why rather
    /// than asserting it.** `Var(k)` grows with the panel while `E[k]²` grows with its square, so
    /// the term's share of the answer falls like `1/N`: at 63 samples of the *same* per-sample
    /// uncertainty the two differ by **0.96%**, and a real 63-sample fit is far more certain
    /// per sample than this — which is where the report's *agree to three decimals* comes from.
    #[test]
    fn dropping_the_variance_term_returns_two_and_a_half_times_the_truth_at_one_sample() {
        let one_sample = midway(0.3, 0.4, 0.3, 1, 1);
        let with_term = CensusMoments::from_posteriors(&one_sample, 1, 1, outbred()).heterozygosity;
        let without = without_the_variance_term(&one_sample, 1, 1);
        // The fit writes `f32`, so the tolerances here are about what the array holds rather than
        // about the arithmetic: 0.4 arrives as 0.40000000596 and 0.3 as 0.30000001192.
        assert!(
            (with_term - 0.4).abs() < 1e-7,
            "the truth at these posteriors is the heterozygous one, 0.4; got {with_term}"
        );
        assert!(
            (without - 1.0).abs() < 1e-7,
            "substituting E[k] = 1 gives 2·1 − 1² = 1; got {without}"
        );
        assert!(
            (without / with_term - 2.5).abs() < 1e-6,
            "dropping the term must return 2.5 times the truth here; got {}",
            without / with_term
        );

        // **And the same posteriors at 63 samples, which is the point.** A test written here
        // instead would pass with the term deleted.
        let sixty_three = midway(0.3, 0.4, 0.3, 63, 1);
        let with_term =
            CensusMoments::from_posteriors(&sixty_three, 63, 1, outbred()).heterozygosity;
        let without = without_the_variance_term(&sixty_three, 63, 1);
        // 63 samples of these posteriors: `E[k] = 63`, `Var(k) = 63 · 0.6 = 37.8`, so the term
        // moves 0.5040 to 0.4992 — an inflation of 0.0048/0.4992, or 0.96%. **Pinned rather than
        // bounded**, because a bound of "under 1%" is also satisfied by zero, and zero is exactly
        // what a deleted term gives.
        let inflation = without / with_term - 1.0;
        assert!(
            (inflation - 0.009_615).abs() < 1e-5,
            "at 63 samples of the same per-sample uncertainty the term is 0.96% of the answer, \
             which is why the fixture that pins it runs at one; got {inflation}"
        );
    }

    /// **A census of certain genotypes is unmoved by the term**, which is what makes the fixture
    /// above about uncertainty rather than about the arithmetic.
    ///
    /// Point-mass posteriors have no variance, so `E[k(2N − k)]` and `2N·E[k] − E[k]²` are the
    /// same number — and the point-mass tests above therefore say nothing about the term. Stated
    /// here so nobody reads them as covering it.
    #[test]
    fn certain_genotypes_leave_the_variance_term_at_zero() {
        let certain = point_masses(&[&[1, 0, 2], &[2, 2, 0]]);
        let with_term = CensusMoments::from_posteriors(&certain, 3, 2, outbred()).heterozygosity;
        let without = without_the_variance_term(&certain, 3, 2);
        assert_eq!(with_term, without);
    }

    /// **The soft count is not a count of positions whose expected copy count is above zero, and
    /// this is the difference** — `doc/devel/ng/spec/ordinary_site_prior_moments.md` §6.2 records
    /// an earlier draft of that spec defining it the second way.
    ///
    /// The fixture is the shape a real census has: five samples whose reads leave each of them a
    /// **1 in 100** posterior of being heterozygous and no more. Every such position has an
    /// expected alternative-copy count of 0.05, which is above zero — so a hard count calls **every
    /// one of them segregating**, and a run over two million positions reports 100% segregating,
    /// which is a property of arithmetic and not a measurement.
    ///
    /// The soft count asks what it means instead: the panel is all-reference with probability
    /// `0.99⁵ = 0.951`, so the position segregates with probability **0.049**. Over a hundred such
    /// positions the soft count reports about five and the hard one reports a hundred.
    #[test]
    fn the_segregating_count_is_soft_and_the_hard_version_reports_every_position() {
        let samples = 5;
        let barely = [0.01_f64, 0.0, 0.0].repeat(samples);
        let soft = probability_that_the_panel_segregates(&barely, samples);
        let expected_copies = alternative_copies_in(&barely, samples).expected_copies;

        assert!(
            expected_copies > 0.0,
            "the hard version's test is `expected copies above zero`, and it passes here at \
             {expected_copies}"
        );
        let all_reference = 0.99_f64.powi(5);
        assert!(
            (soft - (1.0 - all_reference)).abs() < 1e-12,
            "the panel is all-reference at 0.99^5 = {all_reference}, so it segregates at {}; got \
             {soft}",
            1.0 - all_reference
        );
        assert!(
            (soft - 0.049).abs() < 5e-4,
            "about 5 positions in 100, not 100 in 100; got {soft}"
        );

        // And over a census of a hundred such positions, that is what the run reports.
        let mut sums = CensusMomentSums::over(samples);
        for _ in 0..100 {
            sums.add_position(&barely);
        }
        let moments = sums.finish(outbred());
        assert_eq!(moments.positions, 100);
        assert!(
            (moments.segregating_positions - 4.9).abs() < 0.1,
            "about 5 of the 100 positions segregate; the hard count would say 100. Got {}",
            moments.segregating_positions
        );
    }

    /// **The two ends: a panel that is certainly all-reference or certainly all-alternative does
    /// not segregate, and one certain heterozygote makes it certain that it does.**
    ///
    /// The all-alternative end is the one a formula that only subtracted the all-reference term
    /// would get wrong, and no fixture in which some sample is reference can see it.
    #[test]
    fn a_panel_fixed_for_either_allele_does_not_segregate() {
        let all_reference = as_f64(&point_masses(&[&[0, 0, 0]]));
        assert_eq!(
            probability_that_the_panel_segregates(&all_reference, 3),
            0.0
        );

        let all_alternative = as_f64(&point_masses(&[&[2, 2, 2]]));
        assert_eq!(
            probability_that_the_panel_segregates(&all_alternative, 3),
            0.0
        );

        let one_heterozygote = as_f64(&point_masses(&[&[1, 0, 0]]));
        assert_eq!(
            probability_that_the_panel_segregates(&one_heterozygote, 3),
            1.0
        );
        // A panel of one homozygous-alternative sample beside two reference ones segregates too —
        // the alternative allele is there and so is the reference one.
        let one_homozygous_alternative = as_f64(&point_masses(&[&[2, 0, 0]]));
        assert_eq!(
            probability_that_the_panel_segregates(&one_homozygous_alternative, 3),
            1.0
        );
    }

    /// **The carrier posterior does not make a position segregate**, because it contributes no
    /// alternative copies to `k` — the same reason it takes no part in either moment.
    #[test]
    fn a_panel_of_carriers_does_not_segregate() {
        let carriers = [0.0_f64, 0.0, 1.0].repeat(4);
        assert_eq!(probability_that_the_panel_segregates(&carriers, 4), 0.0);
    }

    /// **The spread across positions is the standard error of the mean, and a census that says the
    /// same thing everywhere has none.**
    ///
    /// Two arms. A census of identical positions has a spread of exactly zero, whatever the value —
    /// which is what says the spread is measuring disagreement between positions and not the size
    /// of the numbers. And a census alternating between two values has the sample standard
    /// deviation of those two over `√n`, which is written out here rather than taken from the same
    /// one-pass form the code uses.
    #[test]
    fn the_spread_is_the_standard_error_of_the_mean_across_positions() {
        let heterozygous = as_f64(&point_masses(&[&[1, 0]]));
        let mut identical = CensusMomentSums::over(2);
        for _ in 0..8 {
            identical.add_position(&heterozygous);
        }
        let moments = identical.finish(outbred());
        assert_eq!(moments.frequency_standard_error_floor, 0.0);
        assert_eq!(moments.heterozygosity_standard_error_floor, 0.0);
        assert!(moments.mean_alternative_frequency > 0.0);

        // Four positions at a frequency of 1/4 and four at 0: a mean of 1/8, a sample standard
        // deviation of √(4·(1/8)² · 2/7) = 0.13363, and a standard error of that over √8.
        let reference = as_f64(&point_masses(&[&[0, 0]]));
        let mut alternating = CensusMomentSums::over(2);
        for _ in 0..4 {
            alternating.add_position(&heterozygous);
            alternating.add_position(&reference);
        }
        let moments = alternating.finish(outbred());
        assert!((moments.mean_alternative_frequency - 0.125).abs() < 1e-15);
        let deviation = (8.0 * 0.125_f64 * 0.125 / 7.0).sqrt();
        let truth = deviation / 8.0_f64.sqrt();
        assert!(
            (moments.frequency_standard_error_floor - truth).abs() < 1e-15,
            "got {} against {truth}",
            moments.frequency_standard_error_floor
        );
    }

    /// **A census of one position reports no spread rather than a `NaN`**, because one position has
    /// nothing to disagree with and the `n − 1` a sample variance divides by is zero there.
    #[test]
    fn one_position_has_no_spread() {
        let mut sums = CensusMomentSums::over(2);
        sums.add_position(&as_f64(&point_masses(&[&[1, 0]])));
        let moments = sums.finish(outbred());
        assert_eq!(moments.positions, 1);
        assert_eq!(moments.frequency_standard_error_floor, 0.0);
        assert_eq!(moments.heterozygosity_standard_error_floor, 0.0);
    }

    /// **The heterozygosity's spread carries the inbreeding correction and the frequency's does
    /// not** — because the correction scales the heterozygosity and leaves the frequency alone, and
    /// a spread that did not travel with its own number would describe a different quantity.
    #[test]
    fn the_inbreeding_correction_scales_the_heterozygosity_and_its_spread_together() {
        let mut sums = CensusMomentSums::over(2);
        for copies in [&[1u8, 0][..], &[0, 0], &[2, 1], &[0, 1]] {
            sums.add_position(&as_f64(&point_masses(&[copies])));
        }
        let plain = sums.finish(outbred());
        let selfing = sums.finish(InbreedingF::try_new(0.6).expect("a legal coefficient"));
        // At two individuals the factor is 1 − F/3 = 0.8.
        let lift = 1.0 / 0.8;
        assert!((selfing.heterozygosity / plain.heterozygosity - lift).abs() < 1e-12);
        assert!(
            (selfing.heterozygosity_standard_error_floor
                / plain.heterozygosity_standard_error_floor
                - lift)
                .abs()
                < 1e-12
        );
        assert_eq!(
            selfing.frequency_standard_error_floor,
            plain.frequency_standard_error_floor
        );
        // And the count of segregating positions is untouched by any of it.
        assert_eq!(selfing.segregating_positions, plain.segregating_positions);
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
        let _ = CensusMoments::from_posteriors(&two_positions, 1, 3, outbred());
    }

    /// **A run of no samples is refused**, because a frequency is a share of the panel's
    /// chromosomes and there are none.
    #[test]
    #[should_panic(expected = "no chromosomes to be a share of")]
    fn a_census_over_no_samples_is_refused() {
        let _ = CensusMoments::from_posteriors(&[], 0, 0, outbred());
    }

    /// **A census of no positions returns zeros rather than `NaN`.**
    ///
    /// It is a run whose fit kept nothing, and what must not happen is a division by zero
    /// producing two numbers that look like measurements.
    #[test]
    fn a_census_of_no_positions_returns_zeros() {
        let moments = CensusMoments::from_posteriors(&[], 3, 0, outbred());
        assert_eq!(moments.mean_alternative_frequency, 0.0);
        assert_eq!(moments.heterozygosity, 0.0);
    }
}
