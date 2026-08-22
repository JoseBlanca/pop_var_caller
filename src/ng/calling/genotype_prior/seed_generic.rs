//! The SNP/indel starting point: the run's two concentration numbers — the chromosomes'
//! worth of prior belief attached to the reference allele and to the alternatives — read
//! off the pre-pass's fitted frequency spectrum.
//!
//! At an ordinary site most alternative alleles are rare, so the one chromosome the reference
//! records is almost always the common one and the reference's number is the larger — but only
//! just, and *how much* larger is fitted rather than fixed
//! (`doc/devel/ng/spec/calling_priors.md` §4.1).
//!
//! **This file holds the prediction; the fit that uses it is plan step D2.** Given a candidate
//! pair, [`fill_expected_spectrum`] says what a panel's allele counts would look like. The fit
//! then searches for the pair whose prediction matches the spectrum the pre-pass measured.

use crate::genetics::lgamma;
use crate::ng::types::InbreedingF;

/// The largest concentration — chromosomes' worth of prior belief — this projection will predict a
/// spectrum for.
///
/// **Not a modelling limit but an arithmetic one.** The sum's terms are ratios of gamma functions
/// at arguments that grow with the concentration, and past this point they stop agreeing with each
/// other: measured at 63 individuals and `F = 0.9`, the classes total 1.0000016 at `1e9`, 1.0022 at
/// `1e12` and 1,107 at `1e15`, every entry still finite and non-negative. A run that reached those
/// values would get a spectrum that looks like one.
///
/// `1e6` chromosomes is five hundred thousand diploid individuals, two orders of magnitude past
/// the several thousand this caller commits to, and the total is within 2e-11 of one there.
pub const MAX_PROJECTION_CONCENTRATION: f64 = 1e6;

/// How far below the most likely one a branch split's probability may fall before it is dropped
/// from the sum.
///
/// **Not an approximation in the usual sense: the error has a bound rather than an estimate.** For
/// a fixed number of inbred individuals the classes are themselves a distribution summing to one,
/// so dropping that split moves no class by more than its own probability, and the whole spectrum
/// by no more than the dropped tail.
///
/// **Measured against the term-by-term sum, at tomato's fitted diversity and `F = 0.8`, this value
/// costs nothing at all**: the worst class-by-class disagreement is 4e-13 at 800 and 1,600
/// individuals and 1e-12 at 3,200, which is the same disagreement the untrimmed recurrence has —
/// floating-point accumulation, not the trim. Loosening it does start to cost, and in proportion:
/// `1e-8` gives 1e-9 and `1e-6` gives 1e-7. It buys 3 to 4-fold over not trimming at all
/// (`examples/ng_spectrum_projection_cost.rs`).
const BRANCH_TAIL_TOLERANCE: f64 = 1e-18;

/// A branch split whose probability is below this cannot move any class by more than itself, so it
/// is skipped whatever [`BRANCH_TAIL_TOLERANCE`] says — it is below what a `f64` total can carry.
const NEGLIGIBLE_BRANCH_WEIGHT: f64 = 1e-300;

/// Fill `out` with **what fraction of sites carry the alternative allele on exactly `j` of a
/// panel's chromosomes**, for `j = 0` up to every chromosome, if the population's allele
/// frequencies were drawn from the concentration pair given.
///
/// This is the *frequency spectrum* — the shape the pre-pass measures on real data, and the thing
/// step D2 matches a candidate pair against. `out` holds one entry per class and they sum to one:
/// entry 0 is the fraction of sites where no chromosome carries the alternative allele, entry 1
/// singletons, entry 2 doubletons, and so on.
///
/// ## Where the shape comes from
///
/// Two draws, in order. First the locus's alternative-allele frequency `p` is drawn from the
/// concentration pair — a Beta distribution with those two numbers. Then the panel's genotypes are
/// drawn at that `p` under the **same two-branch model the genotype prior itself uses**
/// (`doc/devel/ng/spec/calling_priors.md` §3.2): with probability `F` an individual's two copies
/// are one ancestral copy counted twice, otherwise they are two independent draws. Averaging over
/// `p` gives the numbers below.
///
/// **Nothing is simulated.** The whole calculation is a finite sum of Beta-function ratios and
/// binomial coefficients — no draws, no quadrature — so the same inputs give the same answer to
/// the last bit at any panel size. The one departure from exactness is the branch-tail trim
/// described under cost, whose error is bounded and, at the tolerance shipped, not measurable.
///
/// **That is what lets step D2's tests state a target rather than only a tolerance** (spec §12
/// tests 5–7): a neutral panel's spectrum can be written down and the fit asked to return `(1, θ)`
/// from it. Note what this does *not* buy — the same function supplies D2's objective, so the
/// tests are not an independent check of the fit's mathematics, only of its search. The
/// independent checks live here, in this file's own oracles.
///
/// ## Why the two-branch model and not `2N` independent chromosomes
///
/// **Because a selfer's panel does not produce an independent-chromosome spectrum at any
/// diversity, and on tomato that is the dominant feature rather than a correction.** In the cohort
/// VCF over 26 accessions, 10,786 sites carry the alternative allele on exactly two chromosomes
/// against 5,142 on exactly one — doubletons outnumbering singletons 2.1 to 1. No
/// independent-chromosome spectrum can do that at any `θ`, because its classes fall monotonically
/// like `θ/k`. Inbreeding is the shape of that spectrum, not a second-order effect on it (spec
/// §4.1). Treating the chromosomes as independent biases `α_ref` **down** with a fixed sign — 12 to
/// 14% at tomato's fitted `F` of 0.8 to 0.9, and 8.6% as far down as `F = 0.6` — which is the
/// number the whole of step D is about.
///
/// ## How it is computed
///
/// Split the panel by which branch each individual took. With `M` of the `N` individuals
/// identical-by-descent, the panel holds `2N − M` **distinct** chromosomes — the other `M` are
/// copies of one of them — and each distinct chromosome carries the alternative allele
/// independently at `p`. So the class probability is a sum over `M`, over how many distinct
/// chromosomes are alternative, and over how many of the duplicated ones are among them.
///
/// **Every term is a product of non-negative factors, and that is the reason for this shape.**
/// Written instead as a polynomial in `p`, the coefficients alternate in sign and grow about
/// 8.3-fold per individual — measured at `F = 0`, they reach 9.5e9 by twelve individuals — so
/// they would cancel catastrophically against a total that must come out at one. These do not
/// cancel at all.
///
/// ## What it costs
///
/// **Paid once per *fit*, not once per run**: this function is the objective step D2 searches, so a
/// multistart fit evaluates it on the order of a hundred times. Measured in release at tomato's
/// fitted diversity and `F = 0.8`, one prediction:
///
/// ```text
///   samples      400     800    1600    3200
///   this        5.8ms  29.9ms   179ms   960ms
///   term-by-term 43.8ms  340ms   2.1s   12.1s
/// ```
///
/// About `N^2.45`, against `N^2.95` for the straight term-by-term sum — so a fit at the top of the
/// committed cohort range is minutes rather than half an hour. Three things buy that, and none of
/// them changes the model:
///
/// - the log-factorials are tabulated once per call rather than rebuilt from `lgamma` per term;
/// - each term is written as **a beta-binomial weight times a hypergeometric one**, both genuine
///   probabilities, so the hypergeometric can be stepped by an exact ratio instead of exponentiated
///   — one exponential per `(split, draw count)` pair rather than one per term;
/// - branch splits far out in their own tail are dropped, at [`BRANCH_TAIL_TOLERANCE`], whose
///   error is bounded rather than estimated.
///
/// **The hypergeometric walk starts at its mode and goes out both ways, and that is not a
/// refinement.** Started at the low end, the first weight underflows to zero long before the ones
/// at the mode become small, and since the walk is multiplicative the whole row then vanishes.
/// Measured at 1,600 individuals when it was written that way: one class came back 5.7e-16 against
/// its true 6.1e-7, and the spectrum lost 3 parts in 10,000 of its mass — with every entry still
/// finite and non-negative, so nothing downstream could have told.
///
/// ## Preconditions and edges
///
/// `out` must hold exactly `2 × individuals + 1` entries, held in release. The panel is diploid:
/// the two-branch model has no state for an individual with two of four copies identical by
/// descent, which spec §3.3 defers to a spec of its own.
///
/// **An alternative concentration of exactly zero is a real answer**, not a degenerate one — a
/// fully invariant cohort — and it puts all the mass in class 0. It is taken early, because the
/// general sum would otherwise reach `lgamma(0)`.
///
/// **Both concentrations are bounded above**, at [`MAX_PROJECTION_CONCENTRATION`]. Past it the sum
/// stops being computable and says nothing about it: measured at 63 individuals and `F = 0.9`, the
/// classes total 1.0000016 at a reference concentration of `1e9`, 1.0022 at `1e12`, and **1,107 at
/// `1e15`** — with every entry still finite and non-negative, so nothing downstream could tell.
/// Step D2 searches over exactly these two axes, which is why this is a release check and not a
/// remark.
pub fn fill_expected_spectrum(
    alpha_ref: f64,
    alpha_alt: f64,
    individuals: u32,
    inbreeding: InbreedingF,
    out: &mut [f64],
) {
    let n = individuals as usize;
    assert_eq!(
        out.len(),
        2 * n + 1,
        "one entry per allele-count class: a panel of {individuals} diploid individuals has \
         {} classes, and `out` holds {}",
        2 * n + 1,
        out.len()
    );
    assert!(
        alpha_ref.is_finite() && alpha_ref > 0.0 && alpha_ref <= MAX_PROJECTION_CONCENTRATION,
        "the reference concentration must be finite, strictly positive and at most \
         {MAX_PROJECTION_CONCENTRATION:e} chromosomes, got {alpha_ref}"
    );
    assert!(
        alpha_alt.is_finite() && (0.0..=MAX_PROJECTION_CONCENTRATION).contains(&alpha_alt),
        "the alternative concentration must be finite, non-negative and at most \
         {MAX_PROJECTION_CONCENTRATION:e} chromosomes, got {alpha_alt}"
    );

    out.fill(0.0);
    if alpha_alt == 0.0 {
        // No alternative allele exists to be drawn, so every site is monomorphic. Taken here
        // because the general sum would otherwise reach `lgamma(0)`.
        out[0] = 1.0;
        return;
    }

    let f = inbreeding.get();
    let concentration_total = alpha_ref + alpha_alt;
    let log_pair_constant = lgamma(concentration_total) - lgamma(alpha_alt) - lgamma(alpha_ref);
    // `ln k!` for every count this call can reach, filled once. The sum asks for binomial
    // coefficients with heavily repeated arguments, and reading them from here rather than calling
    // `lgamma` three times apiece is bit-identical and several times faster. It is the one
    // allocation in this file, and it is per fit rather than per sample per pass — the
    // no-allocation rule of spec §8 governs the row, not the seed.
    let log_factorial: Vec<f64> = (0..=2 * n + 1).map(|k| lgamma(k as f64 + 1.0)).collect();
    let log_binomial = |top: usize, chosen: usize| {
        log_factorial[top] - log_factorial[chosen] - log_factorial[top - chosen]
    };

    // How likely each split of the panel into inbred and outbred individuals is, and the likeliest
    // of them — the tail is measured against that rather than against zero.
    let splits: Vec<Option<f64>> = (0..=n)
        .map(|identical_by_descent| log_branch_split(n, identical_by_descent, f))
        .collect();
    let tail_floor = splits
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max)
        + BRANCH_TAIL_TOLERANCE.ln();

    for (identical_by_descent, split) in splits.iter().enumerate() {
        let Some(log_branch_weight) = *split else {
            continue;
        };
        if log_branch_weight < tail_floor || log_branch_weight < NEGLIGIBLE_BRANCH_WEIGHT.ln() {
            continue;
        }
        let branch_weight = log_branch_weight.exp();
        // The panel's distinct chromosomes: one for each identical-by-descent individual, two for
        // each of the others.
        let distinct = 2 * n - identical_by_descent;
        let single_chromosomes = distinct - identical_by_descent; // 2 × the outbred individuals
        let log_draw_constant = log_pair_constant - lgamma(concentration_total + distinct as f64);

        for alternative_draws in 0..=distinct {
            // How likely it is that exactly this many of the distinct chromosomes carry the
            // alternative allele — a beta-binomial weight, and a probability in its own right.
            let draw_weight = (log_binomial(distinct, alternative_draws)
                + lgamma(alpha_alt + alternative_draws as f64)
                + lgamma(alpha_ref + (distinct - alternative_draws) as f64)
                + log_draw_constant)
                .exp();
            if draw_weight == 0.0 {
                continue;
            }

            // How many of the duplicated chromosomes are among those draws. Each one contributes a
            // second copy, so it moves the class up by one. Given the draw count this is a
            // hypergeometric weight, walked out from its mode by an exact ratio — see the cost
            // section for why the mode and not the low end.
            let lowest = alternative_draws.saturating_sub(single_chromosomes);
            let highest = identical_by_descent.min(alternative_draws);
            let mode = (((alternative_draws + 1) * (identical_by_descent + 1)) / (distinct + 2))
                .clamp(lowest, highest);
            let at_mode = (log_binomial(identical_by_descent, mode)
                + log_binomial(single_chromosomes, alternative_draws - mode)
                - log_binomial(distinct, alternative_draws))
            .exp();
            let scale = branch_weight * draw_weight;
            out[alternative_draws + mode] += scale * at_mode;

            let step_up = |doubled: usize| {
                let rise =
                    ((identical_by_descent - doubled) * (alternative_draws - doubled)) as f64;
                let fall = ((doubled + 1)
                    * ((single_chromosomes + doubled + 1) - alternative_draws))
                    as f64;
                rise / fall
            };
            let mut climbing = at_mode;
            for doubled in mode..highest {
                climbing *= step_up(doubled);
                if climbing == 0.0 {
                    break;
                }
                out[alternative_draws + doubled + 1] += scale * climbing;
            }
            let mut falling = at_mode;
            for doubled in (lowest..mode).rev() {
                falling /= step_up(doubled);
                if falling == 0.0 {
                    break;
                }
                out[alternative_draws + doubled] += scale * falling;
            }
        }
    }
}

/// `ln` of the chance that exactly `identical_by_descent` of `individuals` took the
/// identical-by-descent branch — a binomial at the panel's `F`.
///
/// `None` where that chance is exactly zero, which is every count but one at `F = 0` and at
/// `F = 1`. Returning `None` rather than `−∞` is what keeps `0 × ln 0` — a `NaN` — out of the sum
/// at both ends of the range, where the caller spends most of its life: `F = 0` is an outbred
/// panel and the pre-pass's own default.
fn log_branch_split(
    individuals: usize,
    identical_by_descent: usize,
    inbreeding: f64,
) -> Option<f64> {
    if inbreeding == 0.0 {
        return (identical_by_descent == 0).then_some(0.0);
    }
    if inbreeding == 1.0 {
        return (identical_by_descent == individuals).then_some(0.0);
    }
    let outbred = individuals - identical_by_descent;
    Some(
        lgamma(individuals as f64 + 1.0)
            - lgamma(identical_by_descent as f64 + 1.0)
            - lgamma(outbred as f64 + 1.0)
            + identical_by_descent as f64 * inbreeding.ln()
            + outbred as f64 * (1.0 - inbreeding).ln(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `ln C(n, k)` for the tests' own oracles, written from `lgamma` so it shares no table with
    /// the function under test.
    fn log_binomial(n: usize, k: usize) -> f64 {
        lgamma(n as f64 + 1.0) - lgamma(k as f64 + 1.0) - lgamma((n - k) as f64 + 1.0)
    }

    fn spectrum(alpha_ref: f64, alpha_alt: f64, individuals: u32, inbreeding: f64) -> Vec<f64> {
        let mut out = vec![f64::NAN; 2 * individuals as usize + 1];
        fill_expected_spectrum(
            alpha_ref,
            alpha_alt,
            individuals,
            InbreedingF::try_new(inbreeding).unwrap(),
            &mut out,
        );
        out
    }

    /// **Every site falls in exactly one class, so the classes sum to one** — across panel sizes,
    /// diversities and inbreeding coefficients.
    ///
    /// The property the whole construction rests on. It is not free: the sum runs over three
    /// nested indices with binomial coefficients that reach `10^36.8` at the largest panel here,
    /// and any bound that is off by one, or any branch weight that does not itself sum to one,
    /// shows up as a total away from 1.
    #[test]
    fn the_classes_of_a_panel_sum_to_one() {
        let mut worst = 0.0_f64;
        for individuals in [1_u32, 2, 3, 7, 26, 63] {
            for alpha_alt in [1e-4, 6e-4, 1e-2, 0.5, 2.0] {
                for alpha_ref in [0.5, 1.0, 3.0] {
                    for inbreeding in [0.0, 0.25, 0.6, 0.9, 1.0] {
                        let out = spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                        let total: f64 = out.iter().sum();
                        worst = worst.max((total - 1.0).abs());
                        assert!(
                            (total - 1.0).abs() < 1e-9,
                            "{individuals} individuals, α_ref {alpha_ref}, α_alt {alpha_alt}, F \
                             {inbreeding}: the classes carry {total} rather than one unit of \
                             probability"
                        );
                        assert!(
                            out.iter().all(|c| c.is_finite() && *c >= 0.0),
                            "{individuals} individuals, α_ref {alpha_ref}, α_alt {alpha_alt}, F \
                             {inbreeding}: a class came out negative or not finite: {out:?}"
                        );
                    }
                }
            }
        }
        // Recorded so the `1e-9` above is known to be a bound rather than a hope: the worst
        // departure over this whole grid is 1.7e-13, about 6,000 times inside it. Asserted, so
        // that a change which quietly costs accuracy fails here rather than being absorbed by the
        // looser budget.
        assert!(worst < 1e-12, "worst departure from one was {worst}");
    }

    /// **At no inbreeding the panel is `2N` independent chromosomes, and the spectrum is exactly
    /// the beta-binomial** — the closed form written a second way, sharing no code with the sum
    /// under test.
    ///
    /// This is the oracle for the whole construction at one end of the `F` range: the three nested
    /// sums collapse to a single term per class, and that term is a beta-binomial probability,
    /// which is written here from `lgamma` directly.
    #[test]
    fn at_no_inbreeding_the_spectrum_is_the_beta_binomial() {
        for individuals in [1_u32, 3, 12, 26] {
            for alpha_alt in [1e-3, 1e-2, 0.4] {
                for alpha_ref in [0.5, 1.0, 2.5] {
                    let out = spectrum(alpha_ref, alpha_alt, individuals, 0.0);
                    let chromosomes = 2 * individuals as usize;
                    let total = alpha_ref + alpha_alt;
                    for (class, got) in out.iter().enumerate() {
                        let want = (log_binomial(chromosomes, class)
                            + lgamma(alpha_alt + class as f64)
                            + lgamma(alpha_ref + (chromosomes - class) as f64)
                            + lgamma(total)
                            - lgamma(alpha_alt)
                            - lgamma(alpha_ref)
                            - lgamma(total + chromosomes as f64))
                        .exp();
                        assert!(
                            (got - want).abs() < 1e-12,
                            "{individuals} individuals, α_ref {alpha_ref}, α_alt {alpha_alt}, \
                             class {class}: got {got}, the beta-binomial says {want}"
                        );
                    }
                }
            }
        }
    }

    /// **At full inbreeding every individual is one chromosome counted twice**, so no odd class
    /// can hold anything and the even ones are a beta-binomial over `N` draws rather than `2N`.
    ///
    /// The oracle at the other end of the `F` range, and the check that the doubling is a doubling
    /// rather than a reweighting: a panel of 26 selfers carries the same information as 26
    /// chromosomes, not 52.
    #[test]
    fn at_full_inbreeding_only_the_even_classes_can_hold_anything() {
        for individuals in [1_u32, 4, 26] {
            for alpha_alt in [1e-3, 0.2] {
                let out = spectrum(1.0, alpha_alt, individuals, 1.0);
                for (class, weight) in out.iter().enumerate() {
                    if class % 2 == 1 {
                        assert_eq!(
                            *weight, 0.0,
                            "{individuals} individuals at F = 1: class {class} is odd and must be \
                             empty, got {weight}"
                        );
                    }
                }
                let n = individuals as usize;
                let total = 1.0 + alpha_alt;
                for pairs in 0..=n {
                    let want = (log_binomial(n, pairs)
                        + lgamma(alpha_alt + pairs as f64)
                        + lgamma(1.0 + (n - pairs) as f64)
                        + lgamma(total)
                        - lgamma(alpha_alt)
                        - lgamma(1.0)
                        - lgamma(total + n as f64))
                    .exp();
                    assert!(
                        (out[2 * pairs] - want).abs() < 1e-12,
                        "{individuals} individuals at F = 1, class {}: got {}, a beta-binomial \
                         over {n} draws says {want}",
                        2 * pairs,
                        out[2 * pairs]
                    );
                }
            }
        }
    }

    /// **Doubletons outnumbering singletons is the signature no independent-chromosome spectrum
    /// can produce**, and it is what this function exists to get right (spec §4.1).
    ///
    /// An independent-chromosome spectrum falls monotonically — about `θ/k`, so singletons always
    /// beat doubletons at every diversity. An inbred panel does the opposite, because a selfer
    /// that carries the allele at all usually carries it twice. On the tomato cohort VCF over 26
    /// accessions the real counts are 10,786 doubletons against 5,142 singletons, 2.1 to 1.
    ///
    /// Measured here at that panel size and tomato's fitted diversity: the ratio is 0.50 at
    /// `F = 0`, passes 1 between `F = 0.5` and `F = 0.6`, and reaches 4.7 at `F = 0.9`. So the
    /// real panel's 2.1 sits inside the range this model spans, which the independent-chromosome
    /// one cannot reach at any diversity — checked here across four decades of it.
    #[test]
    fn doubletons_beat_singletons_only_once_the_panel_is_inbred() {
        let individuals = 26;
        let alpha_alt = 6e-4;

        // Independent chromosomes: singletons always win, whatever the diversity.
        for alpha_alt in [1e-5, 1e-4, 1e-3, 1e-2, 0.1] {
            let out = spectrum(1.0, alpha_alt, individuals, 0.0);
            assert!(
                out[2] < out[1],
                "at F = 0 and α_alt {alpha_alt} the doubleton class reached {} against the \
                 singleton's {} — an independent-chromosome spectrum cannot do that",
                out[2],
                out[1]
            );
        }

        // Inbred: the ratio rises with F and crosses one.
        let mut previous = 0.0_f64;
        for inbreeding in [0.0, 0.25, 0.5, 0.6, 0.8, 0.9] {
            let out = spectrum(1.0, alpha_alt, individuals, inbreeding);
            let ratio = out[2] / out[1];
            assert!(
                ratio > previous,
                "the doubleton-to-singleton ratio fell from {previous} to {ratio} when F rose to \
                 {inbreeding}"
            );
            previous = ratio;
        }
        assert!(
            previous > 4.0,
            "at F = 0.9 doubletons should outnumber singletons several times over, got {previous}"
        );
        // And the panel's own measured 2.1 is reached inside the fitted range, not beyond it.
        let at_half = spectrum(1.0, alpha_alt, individuals, 0.5);
        let at_six = spectrum(1.0, alpha_alt, individuals, 0.6);
        assert!(
            at_half[2] / at_half[1] < 1.0 && at_six[2] / at_six[1] > 1.0,
            "the crossing should sit between F = 0.5 and F = 0.6, got {} and {}",
            at_half[2] / at_half[1],
            at_six[2] / at_six[1]
        );
    }

    /// **The neutral shape is `θ/k`, and this reproduces it in the limit the spec names** — small
    /// diversity, independent chromosomes, `α_ref = 1`.
    ///
    /// The spec is careful that the two are *not* an identity: the gap is the panel's own chance
    /// of being polymorphic, about `θ · H(2N)` with `H` the harmonic number. Measured here at
    /// tomato's fitted diversity of 6 in 10,000 over 52 chromosomes, that predicted gap is 3 in a
    /// thousand, and the largest departure of any class from `θ/k` is 0.272% — inside it.
    ///
    /// **This is why the fit's targets are built from this function and never from `θ/k`**
    /// (spec §12 test 5): at a diversity of 1 in 100 the same departure is 4.4%, which is larger
    /// than the effect step D2 is trying to measure.
    #[test]
    fn the_neutral_shape_appears_in_the_small_diversity_limit() {
        let individuals = 26_u32;
        let chromosomes = 2 * individuals as usize;
        let harmonic: f64 = (1..chromosomes).map(|k| 1.0 / k as f64).sum();

        for (alpha_alt, budget) in [(6e-4, 4e-3), (1e-2, 6e-2)] {
            let out = spectrum(1.0, alpha_alt, individuals, 0.0);
            let mut worst = 0.0_f64;
            for (k, weight) in out.iter().enumerate().take(chromosomes).skip(1) {
                let neutral = alpha_alt / k as f64;
                worst = worst.max((weight - neutral).abs() / neutral);
            }
            assert!(
                worst < budget,
                "at θ {alpha_alt} the worst class departed from θ/k by {worst}, above the {budget} \
                 the panel's own polymorphism rate θ·H(2N) = {} allows",
                alpha_alt * harmonic
            );
            // And the departure is of the size the spec predicts, not merely small: it tracks
            // θ·H(2N) rather than being an order of magnitude off it.
            assert!(
                worst > alpha_alt * harmonic / 10.0,
                "the departure {worst} is far below θ·H(2N) = {}, which would mean the limit is \
                 being reached for the wrong reason",
                alpha_alt * harmonic
            );
        }
    }

    /// **A cohort with no alternative allele at all is a real answer**, and it is every site in
    /// class 0 rather than a division by zero.
    #[test]
    fn a_cohort_with_no_alternative_concentration_puts_every_site_in_the_monomorphic_class() {
        for individuals in [1_u32, 5, 26] {
            for inbreeding in [0.0, 0.7, 1.0] {
                let out = spectrum(1.0, 0.0, individuals, inbreeding);
                assert_eq!(out[0], 1.0);
                assert!(out[1..].iter().all(|c| *c == 0.0), "{out:?}");
            }
        }
    }

    /// **One individual is two chromosomes and the spectrum has three classes**, the low end of
    /// the committed cohort range. At `F = 0` they are the beta-binomial's; at `F = 1` the
    /// heterozygous class is empty, because the individual's two copies are one copy.
    #[test]
    fn a_single_individual_still_has_a_spectrum() {
        let outbred = spectrum(1.0, 1e-3, 1, 0.0);
        assert_eq!(outbred.len(), 3);
        assert!((outbred.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!(
            outbred[1] > 0.0,
            "an outbred individual can be heterozygous"
        );

        let selfed = spectrum(1.0, 1e-3, 1, 1.0);
        assert_eq!(selfed[1], 0.0, "at F = 1 there is no heterozygous class");
        assert!((selfed.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    /// **A mis-sized output buffer is refused in release.** The class count is `2N + 1` and a
    /// caller reusing one buffer across panel sizes is the likeliest way to get it wrong; a short
    /// one would leave the top classes carrying the previous panel's weights.
    #[test]
    #[should_panic(expected = "one entry per allele-count class")]
    fn a_mis_sized_class_buffer_is_refused() {
        let mut out = [f64::NAN; 4];
        fill_expected_spectrum(1.0, 1e-3, 3, InbreedingF::try_new(0.0).unwrap(), &mut out);
    }

    /// **An oracle strictly inside the inbreeding range**, which neither exact test above reaches:
    /// both of those sit at `F = 0` or `F = 1`, where the triple sum collapses to one term per
    /// class and a whole family of wrong models still passes.
    ///
    /// Measured, on this code: replacing `class = draws + doubled` with `class = draws` — deleting
    /// the doubling, which is a different model — still passes the sum-to-one test, the
    /// beta-binomial test *and* the neutral-limit test, because the inner sum is a Vandermonde
    /// identity either way.
    ///
    /// What separates them is the first two moments, which the two-branch model fixes exactly:
    ///
    /// ```text
    /// E[j]  = 2N · E[p]                                   — free of F
    /// E[j²] = 4N²·E[p²] + 2N(1 + F)·(E[p] − E[p²])        — carries F linearly
    /// ```
    ///
    /// The mean is `F`-free because inbreeding rearranges copies between individuals without
    /// creating or destroying any; the second moment is where it shows, because a selfer's two
    /// copies agree. `E[p]` and `E[p²]` are the Beta's own moments, written here from the
    /// concentrations directly.
    ///
    /// **It works at one individual**, where the doubleton test has no two classes to compare —
    /// which is the low end of the committed cohort range.
    #[test]
    fn the_first_two_moments_match_the_two_branch_model_at_every_inbreeding_coefficient() {
        for individuals in [1_u32, 2, 5, 17, 40] {
            for (alpha_ref, alpha_alt) in [(1.0, 6e-4), (1.0, 1e-2), (0.5, 0.5), (3.0, 1.0)] {
                for inbreeding in [1e-6, 0.05, 0.25, 0.5, 0.6, 0.9, 0.999] {
                    let out = spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    let n = f64::from(individuals);
                    let total = alpha_ref + alpha_alt;
                    let mean_p = alpha_alt / total;
                    let mean_p_squared = alpha_alt * (alpha_alt + 1.0) / (total * (total + 1.0));

                    let got_mean: f64 = out.iter().enumerate().map(|(j, w)| j as f64 * w).sum();
                    let got_second: f64 = out
                        .iter()
                        .enumerate()
                        .map(|(j, w)| (j * j) as f64 * w)
                        .sum();

                    let want_mean = 2.0 * n * mean_p;
                    let want_second = 4.0 * n * n * mean_p_squared
                        + 2.0 * n * (1.0 + inbreeding) * (mean_p - mean_p_squared);

                    assert!(
                        (got_mean - want_mean).abs() <= 1e-9 * want_mean.max(1.0),
                        "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F {inbreeding}: \
                         mean copies {got_mean}, the model says {want_mean}"
                    );
                    assert!(
                        (got_second - want_second).abs() <= 1e-9 * want_second.max(1.0),
                        "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F {inbreeding}: \
                         second moment {got_second}, the model says {want_second}"
                    );
                }
            }
        }
    }

    /// **The whole model, written a second way and compared class by class** — every assignment of
    /// a branch and a genotype to every individual, enumerated, with each one's chance averaged
    /// over the frequency exactly.
    ///
    /// `5^N` panels, so it is only affordable at a handful of individuals; that is the point of
    /// having it beside the moment identity above, which is cheap at any size but pins two numbers
    /// rather than the whole shape. This shares the Beta-moment formula with the code under test
    /// and shares nothing else — no split by how many individuals were inbred, no hypergeometric
    /// counting, no binomial coefficients at all.
    #[test]
    fn every_panel_enumerated_one_by_one_gives_the_same_spectrum() {
        for individuals in [1_usize, 2, 3, 5] {
            for (alpha_ref, alpha_alt) in [(1.0, 6e-4), (1.0, 0.3), (2.0, 1.5)] {
                for inbreeding in [0.0, 0.3, 0.75, 1.0] {
                    let got = spectrum(alpha_ref, alpha_alt, individuals as u32, inbreeding);
                    let want = enumerated_spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    for (class, (got, want)) in got.iter().zip(&want).enumerate() {
                        assert!(
                            (got - want).abs() <= 1e-12 * want.max(1e-12),
                            "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F \
                             {inbreeding}, class {class}: got {got}, enumeration says {want}"
                        );
                    }
                }
            }
        }
    }

    /// Every `(branch, genotype)` assignment for `individuals` people, each contributing a
    /// monomial in `p` and `1 − p`, averaged over the frequency exactly.
    ///
    /// One individual has five ways to be: identical by descent and alternative (`p`, two copies)
    /// or reference (`1 − p`, none), or outbred with two, one or no alternative copies
    /// (`p²`, `2p(1 − p)`, `(1 − p)²`). The five are walked as base-5 digits.
    fn enumerated_spectrum(
        alpha_ref: f64,
        alpha_alt: f64,
        individuals: usize,
        inbreeding: f64,
    ) -> Vec<f64> {
        // (weight, powers of p, powers of 1 − p, alternative copies)
        let ways: [(f64, usize, usize, usize); 5] = [
            (inbreeding, 1, 0, 2),
            (inbreeding, 0, 1, 0),
            (1.0 - inbreeding, 1 + 1, 0, 2),
            (2.0 * (1.0 - inbreeding), 1, 1, 1),
            (1.0 - inbreeding, 0, 2, 0),
        ];
        let total = alpha_ref + alpha_alt;
        let mut out = vec![0.0; 2 * individuals + 1];
        for panel in 0..5_usize.pow(individuals as u32) {
            let mut weight = 1.0;
            let (mut alt_powers, mut ref_powers, mut class) = (0, 0, 0);
            let mut rest = panel;
            for _ in 0..individuals {
                let (w, a, b, copies) = ways[rest % 5];
                rest /= 5;
                weight *= w;
                alt_powers += a;
                ref_powers += b;
                class += copies;
            }
            if weight == 0.0 {
                continue;
            }
            // The exact average of p^a (1 − p)^b over Beta(α_alt, α_ref).
            let moment = (lgamma(alpha_alt + alt_powers as f64) - lgamma(alpha_alt)
                + lgamma(alpha_ref + ref_powers as f64)
                - lgamma(alpha_ref)
                + lgamma(total)
                - lgamma(total + (alt_powers + ref_powers) as f64))
            .exp();
            out[class] += weight * moment;
        }
        out
    }

    /// **A reference concentration at or below zero is refused in release.** Nothing downstream
    /// could catch it: measured with the check removed, `α_ref = −0.5` at four individuals returns
    /// nine finite, non-negative class weights totalling 1.0097 — a spectrum that looks like a
    /// spectrum. `lgamma`'s own guard is a `debug_assert!`, so release has no other.
    #[test]
    #[should_panic(expected = "the reference concentration must be finite")]
    fn a_reference_concentration_at_zero_is_refused() {
        let mut out = [f64::NAN; 5];
        fill_expected_spectrum(0.0, 1e-3, 2, InbreedingF::try_new(0.5).unwrap(), &mut out);
    }

    /// **A negative alternative concentration is refused in release**, for the same reason: with
    /// the check removed it returns a plausible spectrum summing to 1.0036.
    #[test]
    #[should_panic(expected = "the alternative concentration must be finite")]
    fn a_negative_alternative_concentration_is_refused() {
        let mut out = [f64::NAN; 5];
        fill_expected_spectrum(1.0, -1e-3, 2, InbreedingF::try_new(0.5).unwrap(), &mut out);
    }

    /// **A concentration past the point where the sum is computable is refused**, rather than
    /// returning a spectrum that no longer sums to one. Step D2 searches over these two axes, so
    /// the ceiling has to be a check and not a comment.
    #[test]
    #[should_panic(expected = "at most 1e6 chromosomes")]
    fn a_concentration_past_the_computable_range_is_refused() {
        let mut out = [f64::NAN; 5];
        fill_expected_spectrum(1e15, 1e-3, 2, InbreedingF::try_new(0.9).unwrap(), &mut out);
    }

    /// **The ceiling is where the arithmetic is still good, not where it has already gone.** At the
    /// largest concentration the function accepts, the classes must still sum to one.
    #[test]
    fn the_spectrum_still_sums_to_one_at_the_largest_concentration_accepted() {
        for inbreeding in [0.0, 0.5, 0.9] {
            let out = spectrum(
                MAX_PROJECTION_CONCENTRATION,
                MAX_PROJECTION_CONCENTRATION * 1e-3,
                26,
                inbreeding,
            );
            let total: f64 = out.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-9,
                "at the ceiling and F {inbreeding} the classes carry {total}"
            );
        }
    }

    /// **Skipping the rarest branch splits drops no mass**, which is what makes it a saving rather
    /// than an approximation.
    ///
    /// At an inbreeding coefficient of one in a million, almost every split with more than a
    /// handful of inbred individuals has a probability below `1e-300`, so the skip fires on nearly
    /// the whole range of `M`. Two things must survive that, and neither is about speed:
    ///
    /// - **the classes still sum to one** — a skip that dropped real mass would show here first;
    /// - **the spectrum is within order `F` of the outbred one**, because that is how far apart the
    ///   two genuinely are. Measured at three individuals, the largest class moves by 3.0e-10
    ///   between `F = 0` and `F = 1e-6`, so the bound below is loose by four orders of magnitude
    ///   and is a sanity check rather than the oracle.
    ///
    /// **The oracle for the skip is the moment identity above**, which now runs at `F = 1e-6`: it
    /// pins two numbers exactly, and dropping a branch split that carried real weight moves them.
    #[test]
    fn a_vanishing_inbreeding_coefficient_drops_no_mass() {
        for individuals in [3_u32, 26, 63] {
            let outbred = spectrum(1.0, 6e-4, individuals, 0.0);
            let nearly = spectrum(1.0, 6e-4, individuals, 1e-6);
            let total: f64 = nearly.iter().sum();
            assert!(
                (total - 1.0).abs() < 1e-12,
                "{individuals} individuals at F = 1e-6: the classes carry {total}"
            );
            for (class, (a, b)) in outbred.iter().zip(&nearly).enumerate() {
                assert!(
                    (a - b).abs() <= 1e-4 * a.max(1e-12),
                    "{individuals} individuals, class {class}: F = 0 gives {a}, F = 1e-6 gives {b}"
                );
            }
        }
    }

    /// **What one prediction costs, at the panel sizes the caller commits to.**
    ///
    /// Ignored by default because it measures wall-clock and would be flaky in a suite; run it
    /// with `cargo test --release --lib -- --ignored --nocapture seed_generic::tests::cost`. The
    /// numbers it prints are what the function's doc quotes, and re-running it is how they stay
    /// true.
    ///
    /// **Read it as cost per objective evaluation, not per run.** Step D2 searches this function,
    /// so a fit pays it on the order of a hundred times over.
    #[test]
    #[ignore = "measures wall-clock; run explicitly in release"]
    fn cost_of_one_prediction_by_panel_size() {
        for individuals in [26_u32, 63, 200, 400, 800] {
            let mut out = vec![0.0; 2 * individuals as usize + 1];
            let started = std::time::Instant::now();
            fill_expected_spectrum(
                1.0,
                6e-4,
                individuals,
                InbreedingF::try_new(0.8).unwrap(),
                &mut out,
            );
            let elapsed = started.elapsed();
            println!(
                "COST individuals={individuals} terms≈{} elapsed={:?}",
                (individuals as u64).pow(3) / 3,
                elapsed
            );
            assert!((out.iter().sum::<f64>() - 1.0).abs() < 1e-9);
        }
    }

    /// The same sum, written the slow obvious way: **one exponential per term, nothing stepped and
    /// nothing skipped**.
    ///
    /// This is what the shipped function computed before it was made fast, and it is kept because
    /// it is the version whose agreement with three outside oracles was established — enumeration,
    /// quadrature and a generating function. Every term is built independently from logarithms, so
    /// a term that underflows was genuinely negligible; nothing here can lose a whole row the way
    /// a multiplicative walk can.
    fn term_by_term_spectrum(
        alpha_ref: f64,
        alpha_alt: f64,
        individuals: u32,
        inbreeding: f64,
    ) -> Vec<f64> {
        let n = individuals as usize;
        let mut out = vec![0.0; 2 * n + 1];
        if alpha_alt == 0.0 {
            out[0] = 1.0;
            return out;
        }
        let concentration_total = alpha_ref + alpha_alt;
        let log_pair_constant = lgamma(concentration_total) - lgamma(alpha_alt) - lgamma(alpha_ref);
        for identical_by_descent in 0..=n {
            let Some(log_branch_weight) = log_branch_split(n, identical_by_descent, inbreeding)
            else {
                continue;
            };
            let distinct = 2 * n - identical_by_descent;
            let singles = distinct - identical_by_descent;
            let log_draw_constant =
                log_pair_constant - lgamma(concentration_total + distinct as f64);
            for alternative_draws in 0..=distinct {
                let log_frequency_weight = lgamma(alpha_alt + alternative_draws as f64)
                    + lgamma(alpha_ref + (distinct - alternative_draws) as f64)
                    + log_draw_constant;
                let lowest = alternative_draws.saturating_sub(singles);
                let highest = identical_by_descent.min(alternative_draws);
                for doubled in lowest..=highest {
                    out[alternative_draws + doubled] += (log_branch_weight
                        + log_binomial(identical_by_descent, doubled)
                        + log_binomial(singles, alternative_draws - doubled)
                        + log_frequency_weight)
                        .exp();
                }
            }
        }
        out
    }

    /// **The fast sum and the slow obvious one agree class by class**, over panel sizes, both
    /// concentrations and the whole inbreeding range.
    ///
    /// The shipped function steps its innermost factor by a ratio, tabulates its factorials and
    /// drops branch splits far out in their own tail. None of that changes the model, and this is
    /// where that claim is checked rather than asserted. Measured worst disagreement on this grid
    /// is 2.5e-14 relative.
    ///
    /// **Panel sizes stop at 120 because the slow sum is cubic and this runs in debug.** The sizes
    /// where the difference between the two actually bites are covered by the release-only test
    /// below and by `examples/ng_spectrum_projection_cost.rs`.
    #[test]
    fn the_fast_sum_matches_the_term_by_term_sum() {
        let mut worst = 0.0_f64;
        for individuals in [1_u32, 2, 5, 26, 63, 120] {
            for (alpha_ref, alpha_alt) in [(1.0, 6e-4), (1.0, 1e-2), (0.5, 0.5), (3.0, 2.0)] {
                for inbreeding in [0.0, 1e-6, 0.25, 0.8, 0.999, 1.0] {
                    let fast = spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    let slow = term_by_term_spectrum(alpha_ref, alpha_alt, individuals, inbreeding);
                    for (class, (a, b)) in fast.iter().zip(&slow).enumerate() {
                        let gap = if *b > 0.0 { (a - b).abs() / b } else { a.abs() };
                        worst = worst.max(gap);
                        assert!(
                            gap < 1e-11,
                            "{individuals} individuals, α ({alpha_ref}, {alpha_alt}), F \
                             {inbreeding}, class {class}: fast {a}, term-by-term {b}"
                        );
                    }
                }
            }
        }
        assert!(worst < 1e-12, "worst disagreement was {worst}");
    }

    /// **A panel large enough that the hypergeometric walk's low end underflows.**
    ///
    /// The walk is multiplicative, so a first weight of zero makes the whole row zero. Started at
    /// the low end rather than the mode, that is what happens above about a thousand individuals —
    /// measured at 1,600, one class came back 5.7e-16 against its true 6.1e-7 and the spectrum lost
    /// 3 parts in 10,000 of its mass, with every entry still finite and non-negative.
    ///
    /// Ignored by default because a panel that large costs seconds in a debug build; run it with
    /// `cargo test --release --lib -- --ignored --nocapture seed_generic::tests::the_hypergeometric`.
    #[test]
    #[ignore = "needs a panel of 1,600; run explicitly in release"]
    fn the_hypergeometric_walk_survives_a_panel_that_underflows_its_low_end() {
        let out = spectrum(1.0, 6e-4, 1600, 0.8);
        let mass: f64 = out.iter().sum();
        assert!(
            (mass - 1.0).abs() < 1e-9,
            "the classes carry {mass}; a row lost to underflow shows here first"
        );
        // Every class a panel this size can reach carries something: the walk reaching zero early
        // is what empties one.
        let emptiest = out[..=2 * 1600]
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        assert!(
            emptiest > 1e-30,
            "the emptiest class holds {emptiest}, which is a row the walk failed to fill"
        );
    }
}
