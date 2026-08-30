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
//! **Where this is reached from.** The joint fit's expectation step feeds [`CensusMomentSums`] one
//! position at a time and `JointFit` carries the result; [`CensusMomentsReport::of`] turns that
//! plus the panel's inbreeding coefficient into what a run reports.
//!
//! **What is not built is anything downstream of that.** `RunParameters::assemble` — which takes
//! the seed these two numbers imply — has no caller outside its own tests, so the pre-pass to
//! calling handover does not exist as a whole. Nothing here waits on it.

use super::fit::JointFit;
use crate::ng::parameter_estimation::generic::runs::MIN_WINDOWS_TO_FIT_INBREEDING;
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
        // **A panel of no samples is refused by `CensusMomentSums::over` below, not here.** This
        // function used to carry its own guard with the *same message*, so a `#[should_panic]`
        // test naming that message was satisfied by whichever of the two fired — measured, the
        // guard could be deleted outright and the suite stayed green. One check, one message, one
        // site that a test can pin.
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
/// `f32`s — which over the shipped two-million-position census is **1.2 GB at fifty samples and
/// 1.5 GB at the tomato cohort's sixty-three**, held to be summed once and thrown away. The design
/// asks for the sums instead (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §5), and this is
/// them.
///
/// **⚑ Spec §5 pairs the 1.5 GB figure with fifty samples, and 1.5 GB is the sixty-three-sample
/// number**; `12 × 50 × 2×10⁶` is 1.2 GB and `12 × 63 × 2×10⁶` is 1.512 GB. The same pairing is on
/// `JointFitConfig::genotype_posteriors` and predates this module. Corrected here rather than
/// propagated.
///
/// **And the array's peak is about twice whichever figure applies.** The pass that collects it
/// gathers every chunk's own vector and *then* absorbs them into an exactly-reserved merged one,
/// so both are live at once — the collect is deliberate, because `reduce` may join chunks in any
/// order and the array is in position order.
///
/// **The array route is kept and is not dead**: it is what the sums are checked against, at a
/// relative `1e-6` and **not** to the last bit — the array passes through `f32` on the way out
/// while the sums stay in `f64` throughout
/// ([`fit::whole_fit_tests::the_summed_moments_and_the_stored_array_agree`](super::fit)). It is
/// also what `examples/ng_prior_moments_from_reads.rs` sets the flag for.
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

    /// How many diploid individuals the panel holds — the `N` both moments' finite-panel
    /// corrections are written in.
    #[must_use]
    pub fn samples(&self) -> usize {
        self.samples
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

/// **Where the panel's inbreeding coefficient came from**, which a run must say because the three
/// sources are not interchangeable and one of them carries a circularity
/// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §4.1, §7).
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum InbreedingSource {
    /// **The runs-of-homozygosity estimator** — a two-state model over genome windows returning
    /// the share of the analysable genome lying in stretches where an individual's two copies
    /// descend from one recent ancestor (`parameter_estimation::generic::runs`).
    ///
    /// **This is the source §4.1 prefers**, on three reasons of which the first is decisive: it
    /// reads the *distribution* of heterozygosity along a genome and needs no population
    /// expectation, so nothing about it depends on the diversity this correction is computing. It
    /// is also what the derivation asks for — realized autozygosity is what a run of homozygosity
    /// *is* — and it works at one sample.
    ///
    /// `windows` is how many the model was fitted over, which spec §7 asks a run to print because
    /// the estimator's own floor is
    /// [`MIN_WINDOWS_TO_FIT_INBREEDING`](crate::ng::parameter_estimation::generic::runs::MIN_WINDOWS_TO_FIT_INBREEDING)
    /// — 3,000 — below which what it returns is its own noise.
    ///
    /// **⚠ The below-the-floor warning this report carries cannot fire on a coefficient that came
    /// from `fit_inbreeding`, and that is worth knowing rather than discovering.** That function
    /// **refuses** below the floor — `ParameterEstimationError::InbreedingNotFittable`, naming the
    /// window count and the floor — so a run never gets a thin coefficient to report in the first
    /// place. What the warning guards is a report assembled by hand, or by some future route that
    /// reaches the coefficient without going through that refusal.
    ///
    /// **⚑ What spec §4.2 actually asks for is a warning when the count is *near* the floor, and
    /// it names no number.** This report does not invent one: `runs::resolution_at` exists and its
    /// own documentation forbids being used as a threshold — measured, three fits above it on
    /// genomes with no runs at all would have been called detections. So the count is printed and
    /// the reader judges, which is the same treatment the segregating count and the two-route gap
    /// get. **Where "near" goes is the owner's.**
    RunsOfHomozygosity { windows: u32 },
    /// **The joint fit's homozygote excess** — how much less heterozygous an individual is than
    /// the fitted frequencies predict.
    ///
    /// **⚠ It is circular here and the output must not hide that.** This quantity is measured
    /// against a population expectation the same fit produced, and the correction it feeds divides
    /// a diversity by `1 − F`. `parameter_prepass_generic.md` §6.3 states the rule in as many
    /// words: *"Do not take `F` from the ratio estimator and then compute the cohort's diversity
    /// from it… the ratio estimator needs a diversity to produce `F`, so feeding its `F` back in
    /// returns whatever was assumed."*
    ///
    /// **The joint fit is not the pure ratio estimator** — from two samples up it also sees how
    /// many samples carry the allele at each position, which is real information the ratio has
    /// not got, and it recovers 0.80 from a truth of 0.8 there (report §3.5). **But the direction
    /// of the dependence is the wrong way round**, and the runs estimator has no such dependence
    /// at all.
    ///
    /// It also absorbs population structure: a cohort that is really two subpopulations looks
    /// homozygote-excessive for reasons no individual's parents caused.
    JointFitHomozygoteExcess,
    /// **The user said so** — per sample or one value for the whole panel, including zero.
    ///
    /// A user who knows how their material was bred knows it whatever the cohort size, so this is
    /// not a single-sample feature: a fitted coefficient at three samples is worth overriding for
    /// the same reason it is at one (§4.2, owner's decision of 2026-08-27).
    User,
}

/// **What a run had to correct its heterozygosity with, before the panel's one number is taken
/// from it** — the three shapes a run can arrive in, and the decision of which source it is
/// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §4.1, §4.2).
///
/// **The panel's value is the plain mean over samples, unweighted, whichever source supplied
/// it.** With a per-individual `Fᵢ` the estimator's expectation is `π · (1 − F̄/(2N − 1))` with
/// `F̄` the unweighted mean, so a sample with more census positions covered must not count for
/// more. **⚠ Nothing has tested this on a panel of mixed coefficients** — every drawn panel behind
/// the design shares one, so a weighted rule would have passed every arm of every measurement
/// (spec §4.1, §8's sixth open question).
#[derive(Clone, PartialEq, Debug)]
pub enum PanelInbreeding<'a> {
    /// **The user gave one value for the whole panel, or one per sample.** Not a single-sample
    /// feature: a user who knows how their material was bred knows it whatever the cohort size,
    /// and a fitted coefficient at three samples is worth overriding for the same reason it is at
    /// one (§4.2, owner's decision of 2026-08-27).
    Supplied(&'a [InbreedingF]),
    /// **The per-sample route's runs-of-homozygosity coefficients**, one entry per sample that has
    /// one, each with the count of windows holding sites it was fitted over.
    ///
    /// **A sample can legitimately be missing from this list**: `fit_inbreeding_if_diploid`
    /// declines above and below two genome copies — above two, `F` needs several
    /// identity-by-descent coefficients; below two there are no heterozygotes to be short of — so
    /// a haploid sample contributes nothing rather than a zero. **Averaging a zero in for it would
    /// invent a coefficient**, which is why the mean is over the entries present.
    FittedFromRuns(&'a [(InbreedingF, u32)]),
    /// **Nothing produced a runs coefficient, so the joint fit's own homozygote excess stands
    /// in** — per sample, and circular (see [`InbreedingSource::JointFitHomozygoteExcess`]).
    HomozygoteExcess(&'a [InbreedingF]),
}

impl PanelInbreeding<'_> {
    /// The panel's one coefficient and where it came from.
    ///
    /// # Panics
    ///
    /// **On an empty list, in every arm.** A run claiming a source and supplying nothing from it
    /// is a run-assembly defect, and the alternative — quietly meaning an empty list to zero —
    /// would report a panel of unrelated individuals, which is the answer that makes the
    /// correction vanish.
    #[must_use]
    pub fn for_the_panel(&self) -> (InbreedingF, InbreedingSource) {
        match self {
            Self::Supplied(coefficients) => {
                (mean_of(coefficients, "the user"), InbreedingSource::User)
            }
            Self::FittedFromRuns(fitted) => {
                assert!(
                    !fitted.is_empty(),
                    "a run whose coefficient came from the runs estimator has at least one \
                     sample it was fitted on; an empty list is a run whose per-sample \
                     coefficients went missing between the two routes"
                );
                let coefficients: Vec<InbreedingF> =
                    fitted.iter().map(|(coefficient, _)| *coefficient).collect();
                // **The thinnest sample's window count, not the mean of them.** What the count is
                // printed for is telling a reader whether any of the panel's coefficients rests on
                // too little, and a mean hides one thin sample among sixty-two whole genomes.
                let windows = fitted
                    .iter()
                    .map(|(_, windows)| *windows)
                    .min()
                    .expect("the list is not empty");
                (
                    mean_of(&coefficients, "the runs estimator"),
                    InbreedingSource::RunsOfHomozygosity { windows },
                )
            }
            Self::HomozygoteExcess(coefficients) => (
                mean_of(coefficients, "the joint fit's homozygote excess"),
                InbreedingSource::JointFitHomozygoteExcess,
            ),
        }
    }
}

/// The unweighted mean of a panel's coefficients, refusing an empty panel.
fn mean_of(coefficients: &[InbreedingF], source: &str) -> InbreedingF {
    assert!(
        !coefficients.is_empty(),
        "a run whose inbreeding coefficient came from {source} has at least one value; an empty \
         list would mean to zero, which is the answer that makes the correction vanish"
    );
    let mean = coefficients
        .iter()
        .map(|coefficient| coefficient.get())
        .sum::<f64>()
        / coefficients.len() as f64;
    // Every entry is in `[0, 1)`, so their mean is too — this cannot fail, and saying so is
    // cheaper than a `Result` no caller could act on.
    InbreedingF::try_new(mean).expect("the mean of coefficients each below one is below one")
}
/// **Everything a run owes its reader about the two numbers its SNP/indel prior was built from** —
/// `doc/devel/ng/spec/ordinary_site_prior_moments.md` §7's list, assembled so that a run that used
/// different information cannot look like one that did not.
///
/// **It judges nothing.** Two of the numbers here have thresholds that would be useful and neither
/// has been measured: how few segregating positions is too few (§6.2), and how far apart the two
/// heterozygosity estimates may drift before a fit is suspect (§8's fourth open question). Both
/// measurements need a real census, which this checkout cannot rebuild. So the report prints and
/// the reader decides — and *no threshold applied* is distinguishable in the output from a
/// threshold that never fired.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct CensusMomentsReport {
    /// The two moments the census average produced, with their spreads and the counts behind them.
    pub measured: CensusMoments,
    /// **The same heterozygosity by the other route** — `∫ π(f) · 2f(1−f) df` off the fitted
    /// curve, with no census pass and no inbreeding coefficient
    /// (`JointFit::expected_heterozygosity`).
    ///
    /// **Two routes to one quantity, and that is why both are printed.** They are computed from
    /// the same converged fit: one averages the per-position posteriors, the other integrates the
    /// curve. Measured at 63 samples across three populations and four depths, the curve's number
    /// sits between **1.1% below and 10.7% above** the census average's — and which population is
    /// the wide one is *not* predicted by whether the curve can hold its shape (§7).
    pub curve_heterozygosity: f64,
    /// The coefficient the heterozygosity was corrected by, and where it came from.
    pub panel_inbreeding: f64,
    pub inbreeding_source: InbreedingSource,
    /// How many samples the panel holds — which the one-sample warning below needs, and which
    /// nothing else here branches on.
    pub samples: usize,
}

impl CensusMomentsReport {
    /// **Assemble the report from a converged joint fit and whatever the run had to correct its
    /// heterozygosity with** — the join `doc/devel/ng/spec/ordinary_site_prior_moments.md` §8's
    /// fifth open question is about, taken the way the owner ruled on 2026-08-27.
    ///
    /// ## Where the coefficient comes from, and why nothing new was built for it
    ///
    /// **The coefficient already exists, per sample, and this takes it.** §4.1 prefers the
    /// runs-of-homozygosity estimator over the fit's own homozygote excess, and §8's fifth
    /// question asked how a run on this route obtains one — the runs estimator walks genome
    /// windows in the per-sample histogram route, and this route walks census positions. **The
    /// answer is that both routes run**: `parameter_prepass.md`'s step table gives `F` its own
    /// row — computed after each sample's walk, needing that sample's windowed histogram and
    /// nothing else — and marks it as *not* one of the quantities the two routes produce competing
    /// estimates of. So there is one coefficient, fitted once per sample, and the only thing that
    /// was missing is this: taking those values, meaning them unweighted over the panel, and
    /// handing the mean to [`CensusMomentSums::finish`].
    ///
    /// **Fitting runs from the census positions themselves is a later piece of work and not a
    /// fallback** (owner, 2026-08-27). The arithmetic that says the signal is there is in §8's
    /// fifth question; nothing here anticipates it.
    ///
    /// ## What a run does when there is no runs coefficient
    ///
    /// It passes [`PanelInbreeding::HomozygoteExcess`] and the report says so, with the
    /// circularity warning that source carries and — at one sample — the second warning that the
    /// excess there is 0.000 whatever the truth. **That path stays reachable**: `fit_inbreeding`
    /// *refuses* below its 3,000-window floor rather than returning a thin coefficient, so a run
    /// whose samples are region-restricted has no runs coefficient at all rather than a bad one.
    #[must_use]
    pub fn of(fit: &JointFit, panel_inbreeding: PanelInbreeding<'_>) -> Self {
        let (coefficient, inbreeding_source) = panel_inbreeding.for_the_panel();
        Self {
            measured: fit.census_moments.finish(coefficient),
            curve_heterozygosity: fit.expected_heterozygosity,
            panel_inbreeding: coefficient.get(),
            inbreeding_source,
            samples: fit.census_moments.samples(),
        }
    }

    /// **The curve's heterozygosity over the census average's** — the two-route ratio §7 asks a
    /// run to print.
    ///
    /// **A threshold on it cannot be set at a tenth**: a converged, healthy fit already reaches
    /// 10.7% on one of the three populations measured, and the widest cell is the rare-allele
    /// population, whose shape the curve *can* hold exactly. Calibrating one needs a fit that
    /// genuinely failed to converge, and nothing in this work produced one (§8's fourth open
    /// question). Returns `None` where the census average is zero, which is a cohort with no
    /// variation rather than a disagreement.
    #[must_use]
    pub fn curve_over_census(&self) -> Option<f64> {
        (self.measured.heterozygosity > 0.0)
            .then(|| self.curve_heterozygosity / self.measured.heterozygosity)
    }

    /// **What share of the census positions segregate**, as the soft count over the walked count.
    /// `None` where no position was walked.
    #[must_use]
    pub fn segregating_share(&self) -> Option<f64> {
        (self.measured.positions > 0)
            .then(|| self.measured.segregating_positions / self.measured.positions as f64)
    }

    /// **The warnings this run owes, in the order a reader needs them.** Empty is the ordinary
    /// case and is not itself a claim that anything is well founded.
    #[must_use]
    pub fn warnings(&self) -> Vec<String> {
        let mut warnings = Vec::new();
        match self.inbreeding_source {
            InbreedingSource::RunsOfHomozygosity { windows } => {
                if (windows as usize) < MIN_WINDOWS_TO_FIT_INBREEDING {
                    warnings.push(format!(
                        "the inbreeding coefficient was fitted over {windows} genome windows, \
                         below the {MIN_WINDOWS_TO_FIT_INBREEDING} its estimator needs before what \
                         it returns is a measurement rather than its own noise"
                    ));
                }
            }
            InbreedingSource::JointFitHomozygoteExcess => {
                warnings.push(
                    "the inbreeding coefficient is the joint fit's own homozygote excess, which \
                     is measured against a population expectation this same fit produced — so the \
                     heterozygosity below has been divided by a number that depends on it. The \
                     runs-of-homozygosity estimator has no such dependence and is what the design \
                     prefers"
                        .to_owned(),
                );
                if self.samples == 1 {
                    // **The multiplier a user needs depends on what was already divided out**, and
                    // at one sample the factor applied was `1 − F_printed`. So the residual is
                    // `(1 − F_printed)/(1 − F_true)`, which is `1/(1 − F_true)` only where the
                    // printed coefficient is zero. The design says it *is* zero here — a single
                    // genome's totals cannot identify the excess — but this type's fields are
                    // public and `of` takes whatever it is handed, so the two cases are told
                    // apart rather than assumed. An earlier version printed the zero-coefficient
                    // advice unconditionally and would have told a user holding 0.4 to multiply
                    // by 5 where the right factor is 3.
                    if self.panel_inbreeding == 0.0 {
                        warnings.push(
                            "and at one sample that excess is 0.000 whatever the truth, because a \
                             single genome's totals cannot identify it — so the heterozygosity \
                             below is uncorrected and is a floor. A user who knows their material \
                             is selfing at F should multiply by 1/(1 − F): at F = 0.8 the \
                             population's diversity is five times what is printed"
                                .to_owned(),
                        );
                    } else {
                        warnings.push(format!(
                            "and at one sample the joint fit's excess is 0.000 whatever the truth, \
                             because a single genome's totals cannot identify it — so this {:.3} \
                             did not come from a single genome's own excess. The heterozygosity \
                             below has been divided by 1 − {:.3}; a user who knows their material \
                             is selfing at F should multiply by (1 − {:.3})/(1 − F)",
                            self.panel_inbreeding, self.panel_inbreeding, self.panel_inbreeding
                        ));
                    }
                }
            }
            InbreedingSource::User => {}
        }
        warnings
    }
}

impl std::fmt::Display for CensusMomentsReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let measured = &self.measured;
        writeln!(
            f,
            "the population's two numbers, averaged over {} census positions",
            measured.positions
        )?;
        writeln!(
            f,
            "  mean alternative-allele frequency  {:.4e}   at least ±{:.1e} (a floor, not an \
             interval)",
            measured.mean_alternative_frequency, measured.frequency_standard_error_floor
        )?;
        writeln!(
            f,
            "  heterozygosity                     {:.4e}   at least ±{:.1e} (a floor, not an \
             interval)",
            measured.heterozygosity, measured.heterozygosity_standard_error_floor
        )?;
        match self.segregating_share() {
            Some(share) => writeln!(
                f,
                "  positions that segregate           {:.1} of {} ({:.3}%), a soft count",
                measured.segregating_positions,
                measured.positions,
                100.0 * share
            )?,
            None => writeln!(f, "  positions that segregate           none walked")?,
        }
        writeln!(
            f,
            "  the same heterozygosity off the fitted curve: {:.4e}{}",
            self.curve_heterozygosity,
            match self.curve_over_census() {
                Some(ratio) => format!(" ({ratio:.3}x the census average's)"),
                None => String::new(),
            }
        )?;
        writeln!(
            f,
            "  inbreeding coefficient             {:.3}, from {}",
            self.panel_inbreeding,
            match self.inbreeding_source {
                InbreedingSource::RunsOfHomozygosity { windows } =>
                    format!("the runs-of-homozygosity estimator over {windows} genome windows"),
                InbreedingSource::JointFitHomozygoteExcess =>
                    "the joint fit's own homozygote excess".to_owned(),
                InbreedingSource::User => "the user".to_owned(),
            }
        )?;
        for warning in self.warnings() {
            writeln!(f, "  ⚠ {warning}")?;
        }
        Ok(())
    }
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

    /// **The heterozygosity formula itself, against arithmetic done by hand** — both of its two
    /// parts, at numbers chosen so that each part's absence is visible.
    ///
    /// Every other fixture reaches [`nei_heterozygosity`] through a posterior array, and the two
    /// whole-fit tests reach it through *both* routes at once, so a defect inside it moves every
    /// number in sight by the same factor and nothing fails. Measured: dropping the variance term
    /// and writing the denominator's `2N − 1` as `2N` left the whole suite green.
    ///
    /// Two samples, so four chromosomes. With `E[k] = 1.5` and `Var(k) = 0.75`:
    ///
    /// ```text
    ///   2 (4·1.5 − 1.5² − 0.75) / (4·3)  =  2 (6 − 2.25 − 0.75) / 12  =  6/12  =  0.5
    /// ```
    ///
    /// Drop the variance and it reads `0.625`; write `4·4` in the denominator and it reads
    /// `0.375`. All three are different numbers, which is the point of the fixture.
    #[test]
    fn the_heterozygosity_formula_matches_arithmetic_done_by_hand() {
        assert_eq!(nei_heterozygosity(1.5, 0.75, 4.0), 0.5);
        // The two ways it has been seen to break, named so the numbers above cannot be tuned to
        // whatever the code happens to return.
        assert_ne!(nei_heterozygosity(1.5, 0.0, 4.0), 0.5, "the variance term");
        assert_ne!(
            2.0 * (4.0 * 1.5 - 1.5 * 1.5 - 0.75) / (4.0 * 4.0),
            0.5,
            "the chromosome the pair is not drawn against twice"
        );
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

    /// **The panel's coefficient is the plain mean over samples, and a sample with more census
    /// positions covered does not count for more** — spec §4.1, whose derivation gives the
    /// estimator's expectation as `π · (1 − F̄/(2N − 1))` with `F̄` the *unweighted* mean.
    ///
    /// **⚠ Nothing has tested this on a panel of mixed coefficients but this**, and it is a
    /// fixture rather than a measurement: every drawn panel behind the design shares one
    /// coefficient across its individuals, so a weighted rule would have passed every arm of every
    /// sweep (§8's sixth open question). What this pins is that the code means them, not that
    /// meaning them is right on real data.
    #[test]
    fn the_panels_coefficient_is_the_unweighted_mean_over_its_samples() {
        let coefficients: Vec<InbreedingF> = [0.0, 0.4, 0.8]
            .into_iter()
            .map(|f| InbreedingF::try_new(f).expect("a legal coefficient"))
            .collect();
        let (mean, source) = PanelInbreeding::Supplied(&coefficients).for_the_panel();
        assert!((mean.get() - 0.4).abs() < 1e-15, "got {}", mean.get());
        assert_eq!(source, InbreedingSource::User);
    }

    /// **A haploid sample contributes no coefficient rather than a zero**, and the mean is over
    /// the samples that have one.
    ///
    /// `fit_inbreeding_if_diploid` declines above and below two genome copies — above two, `F`
    /// needs several identity-by-descent coefficients; below two there are no heterozygotes to be
    /// short of. **Averaging a zero in for such a sample would invent a coefficient**, and on a
    /// panel of two where one is haploid it would halve the correction.
    #[test]
    fn a_sample_with_no_coefficient_is_absent_rather_than_zero() {
        let one_fitted = [(InbreedingF::try_new(0.8).expect("legal"), 8_004)];
        let (mean, _) = PanelInbreeding::FittedFromRuns(&one_fitted).for_the_panel();
        assert!((mean.get() - 0.8).abs() < 1e-15, "got {}", mean.get());

        // The same panel with a zero averaged in for the sample that has none would give 0.4.
        let with_an_invented_zero = [
            (InbreedingF::try_new(0.8).expect("legal"), 8_004),
            (InbreedingF::try_new(0.0).expect("legal"), 8_004),
        ];
        let (halved, _) = PanelInbreeding::FittedFromRuns(&with_an_invented_zero).for_the_panel();
        assert!((halved.get() - 0.4).abs() < 1e-15);
    }

    /// **The window count reported for a panel is the thinnest sample's, not the mean of them** —
    /// because what the count is printed for is telling a reader whether *any* of the panel's
    /// coefficients rests on too little evidence, and a mean hides one thin sample among
    /// sixty-two whole genomes.
    ///
    /// Here: sixty-two samples at a tomato genome's 8,004 windows and one at 3,100. The mean is
    /// 7,926 and says nothing; the minimum is 3,100, which clears the estimator's 3,000 floor by
    /// 3.3%.
    #[test]
    fn the_reported_window_count_is_the_thinnest_samples() {
        let mut fitted: Vec<(InbreedingF, u32)> =
            vec![(InbreedingF::try_new(0.8).expect("legal"), 8_004); 62];
        fitted.push((InbreedingF::try_new(0.8).expect("legal"), 3_100));
        let (_, source) = PanelInbreeding::FittedFromRuns(&fitted).for_the_panel();
        assert_eq!(
            source,
            InbreedingSource::RunsOfHomozygosity { windows: 3_100 }
        );
    }

    /// **A run claiming a source and supplying nothing from it is refused, not meant to zero** —
    /// and zero is the answer that makes the correction vanish, which is why this is a panic and
    /// not a default.
    #[test]
    #[should_panic(expected = "makes the correction vanish")]
    fn a_panel_with_no_coefficients_at_all_is_refused() {
        let _ = PanelInbreeding::Supplied(&[]).for_the_panel();
    }

    /// **The same, on the runs arm, with its own message** — because a run that fitted
    /// coefficients per sample and arrived here with none lost them between the two routes, which
    /// is a different defect from a user supplying an empty list.
    #[test]
    #[should_panic(expected = "went missing between the two routes")]
    fn a_runs_source_with_no_samples_is_refused() {
        let _ = PanelInbreeding::FittedFromRuns(&[]).for_the_panel();
    }
    /// A report over a small census, so the tests below can name every number in it.
    fn a_report(inbreeding_source: InbreedingSource, samples: usize) -> CensusMomentsReport {
        let mut sums = CensusMomentSums::over(samples);
        for copies in [1u8, 0, 0, 2] {
            sums.add_position(&as_f64(&point_masses(&[&vec![copies; samples]])));
        }
        CensusMomentsReport {
            measured: sums.finish(outbred()),
            curve_heterozygosity: 7.0e-4,
            panel_inbreeding: 0.0,
            inbreeding_source,
            samples,
        }
    }

    /// **A bigger panel finds more positions segregating, from the same per-sample uncertainty** —
    /// which is the whole reason the count is a sum of probabilities rather than a tally.
    ///
    /// Every other fixture that exercises the count is at **point masses**, where each position
    /// either segregates or does not and the answer is the same at any panel size: `a_report`'s
    /// four positions give a share of exactly one in four whether the panel holds two samples or
    /// two thousand. So nothing showed the count doing the one thing it exists to do.
    ///
    /// Here each sample is heterozygous with probability 0.2 and never carries two copies, so it
    /// shows no alternative copy with probability 0.8 and the panel shows none with `0.8^N`:
    ///
    /// ```text
    ///   two samples:  1 − 0.8²  = 0.36        four samples: 1 − 0.8⁴ = 0.5904
    /// ```
    ///
    /// **This is also the soft count the plan warns about not being** — a tally of positions whose
    /// expected alternative-copy count is above zero would call all four of these segregating, at
    /// either panel size.
    #[test]
    fn a_bigger_panel_finds_more_of_the_same_positions_segregating() {
        let share_at = |samples: usize| {
            let mut sums = CensusMomentSums::over(samples);
            for _ in 0..4 {
                sums.add_position(&as_f64(&midway(0.8, 0.2, 0.0, samples, 1)));
            }
            let moments = sums.finish(outbred());
            moments.segregating_positions / moments.positions as f64
        };
        assert!((share_at(2) - 0.36).abs() < 1e-7, "got {}", share_at(2));
        assert!((share_at(4) - 0.5904).abs() < 1e-7, "got {}", share_at(4));
        // And none of these positions is certain either way, which is what a tally would miss.
        assert!(share_at(2) > 0.0 && share_at(4) < 1.0);
    }

    /// **The report prints both routes to the heterozygosity and the ratio between them, and
    /// judges neither** — `doc/devel/ng/spec/ordinary_site_prior_moments.md` §7.
    ///
    /// **A threshold on that ratio cannot be set at a tenth**: a converged, healthy fit already
    /// shows the curve's number 10.7% above the census average's on one of the three populations
    /// measured, and that population is the one whose shape the curve *can* hold exactly. So the
    /// number is printed and nothing branches on it — which this test pins by checking that a
    /// report whose two routes disagree by a factor of two still produces no warning.
    #[test]
    fn the_two_heterozygosity_routes_are_both_printed_and_neither_is_judged() {
        let mut report = a_report(InbreedingSource::User, 2);
        report.curve_heterozygosity = 2.0 * report.measured.heterozygosity;
        assert!((report.curve_over_census().expect("it segregates") - 2.0).abs() < 1e-12);
        assert!(
            report.warnings().is_empty(),
            "no threshold on this gap has been measured, so a run that shows one must not be \
             warned about it; got {:?}",
            report.warnings()
        );
        let printed = report.to_string();
        assert!(printed.contains("2.000x the census average's"), "{printed}");
        assert!(printed.contains("off the fitted curve"), "{printed}");
    }

    /// **The segregating count and both spreads reach the output, and the spreads say they are
    /// floors** — the one word that stops a reader taking one for an interval.
    #[test]
    fn the_report_carries_the_counts_and_calls_the_spreads_floors() {
        let report = a_report(InbreedingSource::User, 2);
        let printed = report.to_string();
        assert!(printed.contains("a soft count"), "{printed}");
        assert!(printed.contains("positions that segregate"), "{printed}");
        assert_eq!(
            printed.matches("a floor, not an interval").count(),
            2,
            "both spreads are floors and both must say so: {printed}"
        );
        // **One of the four positions segregates, not two.** The fixture's samples are identical,
        // so the position where both carry two alternative copies is fixed *for the alternative*
        // and does not segregate any more than the two all-reference ones do. An earlier version
        // of this assertion said two, which is the mistake the count exists to prevent read from
        // the other side.
        assert!(
            (report.segregating_share().expect("positions were walked") - 0.25).abs() < 1e-12,
            "got {:?}",
            report.segregating_share()
        );
    }

    /// **A coefficient taken from the joint fit's own homozygote excess is circular, and the
    /// output must not hide it** — §4.1, and the one row of §7's table that exists because of a
    /// rule stated elsewhere in as many words.
    ///
    /// The fit measures that excess against a population expectation it produced itself, and the
    /// correction divides a diversity by `1 − F`. So the warning is unconditional on that source,
    /// whatever the cohort size.
    #[test]
    fn the_homozygote_excess_source_warns_about_its_own_circularity() {
        let report = a_report(InbreedingSource::JointFitHomozygoteExcess, 10);
        let warnings = report.warnings();
        assert_eq!(warnings.len(), 1, "got {warnings:?}");
        assert!(warnings[0].contains("measured against a population expectation this same fit"));
        assert!(report.to_string().contains("⚠"));

        // And the other two sources are not circular, so neither warns on that account.
        assert!(a_report(InbreedingSource::User, 10).warnings().is_empty());
        assert!(
            a_report(InbreedingSource::RunsOfHomozygosity { windows: 8_004 }, 10)
                .warnings()
                .is_empty()
        );
    }

    /// **At one sample the homozygote excess is 0.000 whatever the truth, so the run says its
    /// coefficient is a floor and what the diversity would be at a stated one** — §4.2, where a
    /// single genome's census is shown to be drawn from the identical distribution under two
    /// populations whose diversities stand in the ratio `1 − F`.
    ///
    /// **The warning is for one sample and no other**, and that is a measurement rather than a
    /// taste: the fit's coefficient goes from 0.000 at one sample to 0.833 at two against a truth
    /// of 0.8, and stays within 0.03 of the truth from three samples to sixty-three (report §3.5).
    #[test]
    fn one_sample_on_the_homozygote_excess_is_told_its_diversity_is_a_floor() {
        let one = a_report(InbreedingSource::JointFitHomozygoteExcess, 1);
        let warnings = one.warnings();
        assert_eq!(warnings.len(), 2, "got {warnings:?}");
        assert!(warnings[1].contains("is a floor"), "{:?}", warnings[1]);
        assert!(
            warnings[1].contains("five times what is printed"),
            "{:?}",
            warnings[1]
        );

        // Two samples is where the fit starts recovering the coefficient, so the second warning
        // stops there.
        let two = a_report(InbreedingSource::JointFitHomozygoteExcess, 2);
        assert_eq!(two.warnings().len(), 1, "got {:?}", two.warnings());

        // **And a one-sample report carrying a non-zero coefficient is told a different thing**,
        // because the advice above is only right when nothing has been divided out yet. The
        // heterozygosity has already been divided by `1 − 0.4`, so the residual factor a user
        // needs is `(1 − 0.4)/(1 − F)` and not `1/(1 − F)` — telling them to multiply by 5 at a
        // true F of 0.8 where the right factor is 3 is the defect this branch prevents. The fields
        // are public, so this state is constructible even though the design says a single genome's
        // own excess is 0.000.
        let mismatched = CensusMomentsReport {
            panel_inbreeding: 0.4,
            ..a_report(InbreedingSource::JointFitHomozygoteExcess, 1)
        };
        let warnings = mismatched.warnings();
        assert_eq!(warnings.len(), 2, "got {warnings:?}");
        assert!(
            warnings[1].contains("(1 − 0.400)/(1 − F)"),
            "{:?}",
            warnings[1]
        );
        assert!(
            !warnings[1].contains("five times what is printed"),
            "the zero-coefficient advice must not be given here: {:?}",
            warnings[1]
        );
    }

    /// **The runs estimator warns below its own floor of 3,000 windows and not above it** — below
    /// that, `parameter_prepass_generic.md` §6.1 records that what it returns is its own noise.
    ///
    /// **⚠ This branch cannot be reached by a coefficient that came from `fit_inbreeding`**, which
    /// refuses below the floor rather than returning a thin estimate. So what is tested here is a
    /// guard on a hand-assembled report, and the fixture is hand-assembled to match. Said plainly
    /// because a test whose subject cannot arise on the shipped path reads as coverage it is not.
    #[test]
    fn the_runs_estimator_warns_below_its_own_window_floor() {
        let thin = a_report(InbreedingSource::RunsOfHomozygosity { windows: 1_200 }, 4);
        let warnings = thin.warnings();
        assert_eq!(warnings.len(), 1, "got {warnings:?}");
        assert!(
            warnings[0].contains("1200 genome windows"),
            "{:?}",
            warnings[0]
        );

        // A tomato genome is 8,004 windows, which is well above the floor.
        let whole_genome = a_report(InbreedingSource::RunsOfHomozygosity { windows: 8_004 }, 4);
        assert!(whole_genome.warnings().is_empty());
        // And exactly at the floor is not below it.
        let at_the_floor = a_report(
            InbreedingSource::RunsOfHomozygosity {
                windows: MIN_WINDOWS_TO_FIT_INBREEDING as u32,
            },
            4,
        );
        assert!(at_the_floor.warnings().is_empty());
    }

    /// **Where the coefficient came from reaches the printed output**, distinguishably — which is
    /// the whole requirement §7 restates from `calling_priors.md` §4: two runs that used different
    /// information must not look the same.
    #[test]
    fn each_inbreeding_source_prints_differently() {
        let printed: Vec<String> = [
            InbreedingSource::User,
            InbreedingSource::JointFitHomozygoteExcess,
            InbreedingSource::RunsOfHomozygosity { windows: 8_004 },
        ]
        .into_iter()
        .map(|source| a_report(source, 4).to_string())
        .collect();
        assert!(printed[0].contains("from the user"));
        assert!(printed[1].contains("from the joint fit's own homozygote excess"));
        assert!(printed[2].contains("runs-of-homozygosity estimator over 8004 genome windows"));
        assert_ne!(printed[0], printed[1]);
        assert_ne!(printed[1], printed[2]);
    }

    /// **Two chunks merged are one walk** — field for field, not only on the two fields a whole-fit
    /// test happens to read.
    ///
    /// **`merge` is live production arithmetic**: the expectation step splits the census across
    /// cores and `Statistics::absorb` folds every chunk into a fresh empty one, so a field the
    /// merge forgets comes back as **zero** on any run with more than one chunk. Measured before
    /// this test existed: dropping the two sums of squares and the segregating count from `merge`
    /// left the whole 4,910-test suite green, and a real run would then have reported both spreads
    /// and the segregating count as zero — three numbers that look like measurements.
    ///
    /// Only `positions`, `frequency` and `heterozygosity` were pinned at all, and only indirectly,
    /// through the two whole-fit tests.
    #[test]
    fn merging_two_chunks_is_the_same_as_one_walk() {
        let positions = [
            as_f64(&point_masses(&[&[1, 0, 2]])),
            as_f64(&point_masses(&[&[0, 0, 0]])),
            as_f64(&midway(0.3, 0.4, 0.3, 3, 1)),
            as_f64(&point_masses(&[&[2, 2, 1]])),
        ];

        let mut whole = CensusMomentSums::over(3);
        for position in &positions {
            whole.add_position(position);
        }

        let mut first = CensusMomentSums::over(3);
        let mut second = CensusMomentSums::over(3);
        first.add_position(&positions[0]);
        first.add_position(&positions[1]);
        second.add_position(&positions[2]);
        second.add_position(&positions[3]);
        first.merge(&second);

        // **The whole struct, so a field added later cannot escape the check** — this is the
        // assertion the field-by-field version would have let through.
        assert_eq!(first, whole);

        // And the finished numbers agree too, which is what a run reads.
        let selfing = InbreedingF::try_new(0.5).expect("a legal coefficient");
        assert_eq!(first.finish(selfing), whole.finish(selfing));
        assert!(whole.finish(selfing).segregating_positions > 0.0);
        assert!(whole.finish(selfing).frequency_standard_error_floor > 0.0);
    }

    /// **Sums taken over two different panels cannot be added**, because each position's terms are
    /// already divided by the panel's chromosome count.
    ///
    /// Nothing merged two accumulators at all until the test above, so this release-held check was
    /// unreached — a review demoted it and no test noticed.
    #[test]
    #[should_panic(expected = "their sums are over the same sample count")]
    fn merging_sums_over_different_panels_is_refused() {
        let mut two = CensusMomentSums::over(2);
        let three = CensusMomentSums::over(3);
        two.merge(&three);
    }

    /// **A position buffer of the wrong length is refused**, which is the check
    /// [`CensusMomentSums::add_position`]'s own `# Panics` note argues for and which nothing
    /// reached until now: it is read by computed offset, so a short buffer would read one sample's
    /// numbers as another's rather than fail.
    #[test]
    #[should_panic(expected = "three numbers a sample, so")]
    fn a_position_buffer_of_the_wrong_length_is_refused() {
        let mut sums = CensusMomentSums::over(3);
        sums.add_position(&[0.5, 0.0, 0.0, 0.5, 0.0, 0.0]);
    }

    /// **Posteriors that do not leave room for a homozygous-reference genotype are clamped, not
    /// carried through negative** — the guard inside the all-reference product.
    ///
    /// Every other fixture here has `P(het) + P(both) ≤ 1` by construction, so the clamp was
    /// unreachable and could be deleted with the suite green. A sample at `(0.6, 0.6)` has no room
    /// left: `1 − 0.6 − 0.6` is `−0.2`, and **two negatives multiply to a positive** — the
    /// unclamped all-reference term for two such samples is `+0.04`, a positive probability of a
    /// state neither sample can be in, subtracted from the segregating count.
    ///
    /// Both samples also carry `P(both non-reference) = 0.6`, so the all-alternative term is
    /// `0.36` whatever the clamp does. The clamp is therefore worth exactly the `0.04`:
    /// **0.64 with it, 0.60 without**.
    ///
    /// The posteriors are not reachable from a converged fit; what this guards is a buffer built by
    /// hand, and the point is that the failure is silent rather than loud.
    ///
    /// **⚠ The outer `clamp(0.0, 1.0)` on the returned probability is still unreachable** — with
    /// the inner guard in place the two products are each in `[0, 1]` and cannot sum above one, so
    /// only floating-point rounding can move it. Said rather than left for the next reviewer to
    /// find again.
    #[test]
    fn posteriors_leaving_no_room_for_a_reference_genotype_are_clamped() {
        let impossible = [0.6_f64, 0.6, 0.0, 0.6, 0.6, 0.0];
        let segregating = probability_that_the_panel_segregates(&impossible, 2);
        assert!(
            (0.0..=1.0).contains(&segregating),
            "the probability that a position segregates is a probability; got {segregating}"
        );
        // 1 − 0 − 0.36 with the clamp; 1 − 0.04 − 0.36 without it.
        assert!((segregating - 0.64).abs() < 1e-15, "got {segregating}");
    }

    /// **A census that says almost the same thing at every position has a spread of zero rather
    /// than the square root of a negative number.**
    ///
    /// The one-pass variance `Σx²/n − (Σx/n)²` cancels to a small *negative* residue when every
    /// position carries nearly the same value, which is exactly what a nearly invariant cohort is —
    /// and `sqrt` of it is `NaN`. Every other fixture here sits at values where the cancellation is
    /// exact, so the clamp was unreachable and could be deleted with the suite green.
    #[test]
    fn a_nearly_invariant_census_has_a_spread_of_zero_rather_than_a_nan() {
        let mut sums = CensusMomentSums::over(4);
        let barely = [1e-7_f64, 0.0, 0.0].repeat(4);
        for _ in 0..10_000 {
            sums.add_position(&barely);
        }
        let moments = sums.finish(outbred());
        assert!(
            moments.frequency_standard_error_floor.is_finite()
                && moments.frequency_standard_error_floor >= 0.0,
            "got {}",
            moments.frequency_standard_error_floor
        );
        assert!(
            moments.heterozygosity_standard_error_floor.is_finite()
                && moments.heterozygosity_standard_error_floor >= 0.0,
            "got {}",
            moments.heterozygosity_standard_error_floor
        );
    }

    /// **A cohort with no variation at all has no two-route ratio and a run over nothing has no
    /// segregating share** — both come back `None` rather than dividing by zero.
    ///
    /// Nothing constructed either state, so both guards could be replaced by unconditional `Some`
    /// with the suite green — and the `Display` arm that prints *none walked* was dead code.
    #[test]
    fn a_report_with_nothing_to_divide_by_returns_none_rather_than_a_nan() {
        let invariant = CensusMomentsReport {
            measured: CensusMomentSums::over(2).finish(outbred()),
            curve_heterozygosity: 7.0e-4,
            panel_inbreeding: 0.0,
            inbreeding_source: InbreedingSource::User,
            samples: 2,
        };
        assert_eq!(invariant.measured.positions, 0);
        assert_eq!(invariant.measured.heterozygosity, 0.0);
        assert_eq!(invariant.curve_over_census(), None);
        assert_eq!(invariant.segregating_share(), None);
        let printed = invariant.to_string();
        assert!(printed.contains("none walked"), "{printed}");
        assert!(
            !printed.contains("NaN") && !printed.contains("inf"),
            "{printed}"
        );
    }

    /// **Both moments' own values reach the printed output, on their own lines.**
    ///
    /// Nothing asserted that either number was printed at all: the `Display` tests checked fixed
    /// substrings and the source labels, so **exchanging the two value expressions left the suite
    /// green** — a report calling the heterozygosity a mean frequency and the frequency a
    /// heterozygosity.
    #[test]
    fn the_printed_report_carries_each_moment_on_its_own_line() {
        // A fixture whose two moments differ, so the swap is visible: one heterozygote in a panel
        // of two is a frequency of 1 in 4 and a heterozygosity of 1 in 2.
        let mut sums = CensusMomentSums::over(2);
        sums.add_position(&as_f64(&point_masses(&[&[1, 0]])));
        sums.add_position(&as_f64(&point_masses(&[&[0, 0]])));
        let report = CensusMomentsReport {
            measured: sums.finish(outbred()),
            curve_heterozygosity: 7.0e-4,
            panel_inbreeding: 0.0,
            inbreeding_source: InbreedingSource::User,
            samples: 2,
        };
        let printed = report.to_string();
        let frequency_line = printed
            .lines()
            .find(|line| line.contains("mean alternative-allele frequency"))
            .expect("the frequency has a line");
        let heterozygosity_line = printed
            .lines()
            .find(|line| line.trim_start().starts_with("heterozygosity"))
            .expect("the heterozygosity has a line");
        assert!(
            frequency_line.contains(&format!(
                "{:.4e}",
                report.measured.mean_alternative_frequency
            )),
            "the frequency line does not carry the frequency: {frequency_line}"
        );
        assert!(
            heterozygosity_line.contains(&format!("{:.4e}", report.measured.heterozygosity)),
            "the heterozygosity line does not carry the heterozygosity: {heterozygosity_line}"
        );
        // And the two really are different numbers here, so the check above can fail.
        assert!(
            report.measured.mean_alternative_frequency < report.measured.heterozygosity,
            "a fixture where the two agree cannot see them swapped"
        );
    }

    /// **Samples that disagree with each other each contribute their own uncertainty**, and their
    /// own share of the all-reference product.
    ///
    /// Every other fixture with non-point-mass posteriors gives **every sample the same three
    /// numbers**, and every fixture whose samples differ is at point masses, where the variance is
    /// identically zero. So the per-sample loop was only ever exercised over identical uncertain
    /// samples or over differing certain ones — never over differing uncertain ones, which is what
    /// a real census at low depth is.
    ///
    /// Here: one sample certain of being heterozygous, one at `(0.5, 0.25)`, one certain of
    /// carrying nothing. `E[k] = 1 + 1.0 + 0 = 2`; the variances are `0`, `0.5 + 4(0.25) − 1 = 0.5`
    /// and `0`, so `Var(k) = 0.5`.
    #[test]
    fn samples_that_disagree_each_contribute_their_own_uncertainty() {
        let mixed = [1.0_f64, 0.0, 0.0, 0.5, 0.25, 0.0, 0.0, 0.0, 0.0];
        let copies = alternative_copies_in(&mixed, 3);
        assert!(
            (copies.expected_copies - 2.0).abs() < 1e-15,
            "got {}",
            copies.expected_copies
        );
        assert!(
            (copies.copy_count_variance - 0.5).abs() < 1e-15,
            "got {}",
            copies.copy_count_variance
        );
        // The all-reference product is `0 · 0.25 · 1 = 0` and the all-alternative one
        // `0 · 0.25 · 0 = 0`, so this position certainly segregates — and it does so because of
        // the *first* sample, which no product over one shared posterior could show.
        assert_eq!(copies.segregating, 1.0);
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
