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
//! position at a time and `JointFit` carries the result; `CensusMomentsReport::of` turns that
//! plus the panel's inbreeding coefficient into what a run reports.
//!
//! **What is not built is anything downstream of that.** `RunParameters::assemble` — which takes
//! the seed these two numbers imply — has no caller outside its own tests, so the pre-pass to
//! calling handover does not exist as a whole. Nothing here waits on it.

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

/// **The two moments as running sums, so a run never stores the array they are summed from.**
///
/// One of these lives on the fit's per-chunk statistics, is fed one position at a time from the
/// expectation step's own scratch buffer, and is merged chunk into chunk the way every other sum
/// there is. What it holds is a handful of `f64`s, whatever the census and whatever the cohort.
///
/// ## Why this and not the stored array
///
/// The fit can be asked to keep every position's genotype posteriors
/// (`JointFitConfig::genotype_posteriors`), and `CensusMoments::from_posteriors` reduces that
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

    /// The fit writes `f32`; every helper here reads `f64`. Widening is exact.
    fn as_f64(narrow: &[f32]) -> Vec<f64> {
        narrow.iter().map(|value| f64::from(*value)).collect()
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
}
