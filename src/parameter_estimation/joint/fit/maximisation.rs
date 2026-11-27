//! The maximisation step: every parameter's own maximisation over the counts one pass
//! accumulated, with the reads untouched.

use crate::types::Ploidy;

use super::beta_quadrature::{BetaQuadrature, trigamma};
use super::{
    BETA_SHAPE_BOUNDS, CLEAN_ERROR_BOUNDS, DUPLICATED_SHARE_BOUNDS, HOM_EXCESS_BOUNDS,
    JointFitConfig, NOISY_ERROR_BOUNDS, NOISY_SHARE_BOUNDS, P_FIXED_ALT_BOUNDS, P_INVARIANT_BOUNDS,
    Parameters, ReadTally, Statistics, candidate_read_probability, count_times_ln,
    genotype_frequencies, reference_read_probability,
};

/// Every parameter's own maximisation, over the counts one pass accumulated. Returns the
/// largest relative move any parameter made — what the per-pass trace ([`PassSummary`]) records;
/// since plan step B1 it no longer decides when the fit stops.
pub(super) fn maximisation(
    parameters: &mut Parameters,
    statistics: &Statistics,
    config: &JointFitConfig,
) -> f64 {
    let moved = std::cell::Cell::new(0.0_f64);
    let note = |before: f64, after: f64| {
        let scale = before.abs().max(1e-6);
        moved.set(moved.get().max((after - before).abs() / scale));
    };
    // **A share is judged against a floor rather than against itself.** A class the data does
    // not want shrinks geometrically — halving on every pass — and halving is a relative move
    // of one half however small the number has become, so the move of such a share would read as
    // large forever. One position in ten thousand is the point below which a share stops being a
    // quantity anyone reads.
    let note_share = |before: f64, after: f64| {
        let scale = before.abs().max(1e-4);
        moved.set(moved.get().max((after - before).abs() / scale));
    };

    // The share of positions that are mismapped, and the density's two masses: closed form.
    let positions = statistics.positions.max(1.0);
    let before = parameters.noisy_share;
    parameters.noisy_share =
        (statistics.noisy / positions).clamp(NOISY_SHARE_BOUNDS.0, NOISY_SHARE_BOUNDS.1);
    note_share(before, parameters.noisy_share);

    let branch_total = (statistics.invariant + statistics.fixed_alt + statistics.segregating)
        .max(f64::MIN_POSITIVE);
    let before = parameters.density.p_invariant;
    parameters.density.p_invariant =
        (statistics.invariant / branch_total).clamp(P_INVARIANT_BOUNDS.0, P_INVARIANT_BOUNDS.1);
    note(before, parameters.density.p_invariant);
    let before = parameters.density.p_fixed_alt;
    parameters.density.p_fixed_alt =
        (statistics.fixed_alt / branch_total).clamp(P_FIXED_ALT_BOUNDS.0, P_FIXED_ALT_BOUNDS.1);
    note(before, parameters.density.p_fixed_alt);

    // The duplicated class's share of positions, and the Beta its carrier frequency is drawn
    // from. **Bounded well below a half**: a class that grows without bound starts explaining
    // real variants, which is the failure mode measured when duplications are mostly private.
    if let Some(duplicated) = parameters.duplicated.as_mut() {
        let before = duplicated.share;
        duplicated.share = (statistics.duplicated / positions)
            .clamp(DUPLICATED_SHARE_BOUNDS.0, DUPLICATED_SHARE_BOUNDS.1);
        note_share(before, duplicated.share);
        if statistics.duplicated > 0.0 {
            let (a, b) = step_beta_shapes(
                statistics.carrier_shape_slopes,
                statistics.duplicated,
                duplicated.carrier_a,
                duplicated.carrier_b,
            );
            note(duplicated.carrier_a, a);
            note(duplicated.carrier_b, b);
            duplicated.carrier_a = a;
            duplicated.carrier_b = b;
        }
    }

    // The Beta's two shapes, one step along the log-likelihood's own slope in them.
    if statistics.segregating > 0.0 {
        let (a, b) = step_beta_shapes(
            statistics.density_shape_slopes,
            statistics.segregating,
            parameters.density.a,
            parameters.density.b,
        );
        note(parameters.density.a, a);
        note(parameters.density.b, b);
        parameters.density.a = a;
        parameters.density.b = b;
    }

    // Each read group's two error rates, each over nine accumulated read counts.
    for group in 0..parameters.clean.len() {
        for (class, bounds) in [(0_usize, CLEAN_ERROR_BOUNDS), (1, NOISY_ERROR_BOUNDS)] {
            let tallies = &statistics.reads[group][class];
            let current = if class == 0 {
                parameters.clean[group]
            } else {
                parameters.noisy[group]
            };
            let fitted = maximise_error_rate(tallies, config.ploidy, bounds, current);
            note(current, fitted);
            if class == 0 {
                parameters.clean[group] = fitted;
            } else {
                parameters.noisy[group] = fitted;
            }
        }
    }
    // **The two rates are not put in order** (`parameter_prepass_joint_fit.md` §3.3, decided
    // 2026-10-17). A library whose mismapped-position rate comes out below its ordinary one keeps
    // both: swapping them alone, with the share and the other libraries' rates unswapped, handed
    // the ordinary positions a rate fitted to the few mismapped ones and lowered the
    // log-likelihood by tens of units a pass. The cap of 0.5 on the mismapped share is what makes
    // that class the minority one.

    // Each sample's homozygote excess, over its expected genotype counts at each node.
    let quadrature = BetaQuadrature::new(
        parameters.density.a,
        parameters.density.b,
        config.quadrature_nodes,
    );
    // **At one sample the homozygote excess is not identified and is not moved.** Nothing
    // separates "this individual is inbred" from "the population's frequencies are what they
    // are" when there is one individual, so a fit that searched it anyway would wander without
    // converging and hand back a plausible number. `fit_jointly` marks it as not fitted.
    if fits_homozygote_excess(samples_with_reads(statistics)) {
        for (s, counts) in statistics.genotypes.iter().enumerate() {
            let fitted = maximise_hom_excess(counts, &quadrature, parameters.hom_excess[s]);
            note(parameters.hom_excess[s], fitted);
            parameters.hom_excess[s] = fitted;
        }
    }
    moved.get()
}

/// The error rate that best explains the expected read counts, by golden-section search over
/// a concave one-dimensional objective.
///
/// **This is where the cost of the route is contained.** The objective reads nine numbers, so
/// the twenty evaluations a search costs are twenty passes over nine numbers rather than
/// twenty passes over the reads.
fn maximise_error_rate(
    tallies: &[ReadTally; 3],
    ploidy: Ploidy,
    bounds: (f64, f64),
    current: f64,
) -> f64 {
    let score = |rate: f64| {
        let mut total = 0.0;
        for (j, tally) in tallies.iter().enumerate() {
            let on_candidate = candidate_read_probability(j as u8, ploidy, rate);
            let on_reference = reference_read_probability(j as u8, ploidy, rate);
            total += count_times_ln(tally.candidate, on_candidate)
                + count_times_ln(tally.neither, rate / 3.0)
                + count_times_ln(tally.reference, on_reference);
        }
        total
    };
    let observations: f64 = tallies
        .iter()
        .map(|t| t.candidate + t.neither + t.reference)
        .sum();
    if observations <= 0.0 {
        return current;
    }
    golden_section(&score, bounds.0, bounds.1)
}

/// **Whether the fit moves each sample's homozygote excess**: only with two samples or more that
/// carry reads at an ordinary position. At one sample nothing separates an inbred individual from a
/// population whose frequencies are what they are, so the excess is held where it starts — and,
/// holding it, the fit has no standard error for it either ([`standard_errors`]). **A sample with no
/// such read is not counted**: one sample with reads beside one whose census holds repeat tracts
/// only is a one-sample cohort, and counting both had the first sample's excess come back fitted,
/// with an error, where alone it is held.
pub(super) fn fits_homozygote_excess(samples_with_reads: usize) -> bool {
    samples_with_reads >= 2
}

/// How many samples have at least one read at an ordinary position — what
/// [`fits_homozygote_excess`] counts.
pub(super) fn samples_with_reads(statistics: &Statistics) -> usize {
    statistics
        .with_reads
        .iter()
        .filter(|&&reads| reads > 0)
        .count()
}

/// The homozygote excess that best explains one sample's expected genotype counts.
fn maximise_hom_excess(counts: &[[f64; 3]], quadrature: &BetaQuadrature, current: f64) -> f64 {
    let observations: f64 = counts.iter().flatten().sum();
    if observations <= 0.0 {
        return current;
    }
    let score = |excess: f64| {
        let mut total = 0.0;
        for (node, row) in counts.iter().enumerate() {
            let prior = genotype_frequencies(quadrature.nodes[node], excess);
            for (j, count) in row.iter().enumerate() {
                total += count_times_ln(*count, prior[j]);
            }
        }
        total
    };
    golden_section(&score, HOM_EXCESS_BOUNDS.0, HOM_EXCESS_BOUNDS.1)
}

/// **How narrow the golden section's bracket must be for it to stop**, as a share of one plus the
/// size of its upper end. It returns the bracket's middle, so a maximum on an end of the interval
/// comes back up to half this beside the end rather than on it — which the box solve of the Newton
/// step treats as on it (`standard_errors::BESIDE_AN_END`).
pub(super) const GOLDEN_SECTION_WIDTH: f64 = 1e-9;

/// The maximum of a unimodal function on `[low, high]`, to about six figures.
pub(super) fn golden_section(score: &dyn Fn(f64) -> f64, low: f64, high: f64) -> f64 {
    const INVERSE_PHI: f64 = 0.618_033_988_749_895;
    let (mut low, mut high) = (low, high);
    let mut c = high - INVERSE_PHI * (high - low);
    let mut d = low + INVERSE_PHI * (high - low);
    let (mut fc, mut fd) = (score(c), score(d));
    for _ in 0..60 {
        if high - low < GOLDEN_SECTION_WIDTH * (1.0 + high.abs()) {
            break;
        }
        if fc > fd {
            high = d;
            d = c;
            fd = fc;
            c = high - INVERSE_PHI * (high - low);
            fc = score(c);
        } else {
            low = c;
            c = d;
            fc = fd;
            d = low + INVERSE_PHI * (high - low);
            fd = score(d);
        }
    }
    0.5 * (low + high)
}

/// **A Beta's two shapes, one Newton step along the log-likelihood's own slope in them.**
///
/// `slope` is the log-likelihood's slope in `[a, b]` summed over the positions — the slope of
/// the likelihood the quadrature rule computes, nodes and weights moving with the shapes
/// ([`Statistics::density_shape_slopes`]). The step divides it by the curvature `positions` Beta
/// draws would give if each position's frequency were seen — `positions` times the matrix of
/// trigamma terms `[[ψ′(a) − ψ′(a+b), −ψ′(a+b)], [−ψ′(a+b), ψ′(b) − ψ′(a+b)]]`, with `positions`
/// the posterior count of positions the Beta covers. That curvature is the complete-data one, so
/// this is one step of the EM gradient algorithm (Lange 1995), **at rest only where the
/// likelihood's slope is zero** — where the update it replaces solved the digamma form, whose
/// integral on the rule misses that slope (`information` module doc, "The Beta shapes' slopes"),
/// and came to rest 0.65 to 2.3 standard errors from the likelihood's flat point
/// (`fit_precision_a4_2026-09-27.md` §2.3).
///
/// **It takes more passes than the update it replaces**: one step a pass, where that update solved
/// its two equations to convergence each pass. Measured on drawn cohorts of 20 samples at 3 reads,
/// an accelerated start converges in 18 to 30 passes where it took 12 to 18
/// (`fit_precision_a7_2026-09-28.md` §3).
///
/// **At a bound the step is projected.** A shape whose step would leave its interval is held at the
/// bound it crosses, and the other takes the step the same curvature gives with that shape held —
/// so the free shape still comes to rest only where its own slope is zero. Clamping each shape of
/// the joint step instead left the free one where its slope was not zero: a joint step's cross term
/// assumes both shapes move (plan step A7's review, where a density wanting its second shape above
/// 50 stopped its first 4.1 units from the maximum and lost log-likelihood from pass to pass). A
/// slope that is not a number moves nothing.
pub(super) fn step_beta_shapes(slope: [f64; 2], positions: f64, a: f64, b: f64) -> (f64, f64) {
    let (low, high) = BETA_SHAPE_BOUNDS;
    let within = |shape: f64| (low..=high).contains(&shape);
    if !(slope[0].is_finite() && slope[1].is_finite()) {
        return (a, b);
    }
    let ab = a + b;
    let h11 = positions * (trigamma(a) - trigamma(ab));
    let h22 = positions * (trigamma(b) - trigamma(ab));
    let h12 = -positions * trigamma(ab);
    let determinant = h11 * h22 - h12 * h12;
    if !(determinant > 0.0 && determinant.is_finite()) {
        return (a, b);
    }
    let joint = (
        a + (slope[0] * h22 - slope[1] * h12) / determinant,
        b + (slope[1] * h11 - slope[0] * h12) / determinant,
    );
    if within(joint.0) && within(joint.1) {
        return joint;
    }
    // One shape held at the bound its step crosses; the other's step with it held.
    let a_given = |held_b: f64| a + (slope[0] - h12 * (held_b - b)) / h11;
    let b_given = |held_a: f64| b + (slope[1] - h12 * (held_a - a)) / h22;
    if !within(joint.1) {
        let held_b = joint.1.clamp(low, high);
        let free_a = a_given(held_b);
        if within(free_a) {
            return (free_a, held_b);
        }
    }
    if !within(joint.0) {
        let held_a = joint.0.clamp(low, high);
        let free_b = b_given(held_a);
        if within(free_b) {
            return (held_a, free_b);
        }
    }
    (joint.0.clamp(low, high), joint.1.clamp(low, high))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The shapes' step is a Newton step with the complete-data curvature, at rest where the
    /// slope is zero, and projected at a bound.** A slope built as the curvature times a known move
    /// lands on that move exactly; a zero slope moves nothing, to the bit; a slope far beyond both
    /// bounds stops at them; a slope that is not a number moves nothing. **With the second shape at
    /// its bound and its slope pushing out**, the second stays and the first takes its own step,
    /// `slope / curvature` — so it moves not at all where its own slope is zero, which clamping each
    /// shape of the joint step does not give.
    #[test]
    fn the_shapes_step_follows_the_slope_by_the_complete_data_curvature() {
        let (a, b, weight) = (0.7, 2.5, 1_000.0);
        let (to_a, to_b) = (0.8, 2.3);
        let ab = a + b;
        let h11 = weight * (trigamma(a) - trigamma(ab));
        let h22 = weight * (trigamma(b) - trigamma(ab));
        let h12 = -weight * trigamma(ab);
        let slope = [
            h11 * (to_a - a) + h12 * (to_b - b),
            h12 * (to_a - a) + h22 * (to_b - b),
        ];
        let (stepped_a, stepped_b) = step_beta_shapes(slope, weight, a, b);
        assert!(
            (stepped_a - to_a).abs() < 1e-12 && (stepped_b - to_b).abs() < 1e-12,
            "stepped to ({stepped_a}, {stepped_b})"
        );
        assert_eq!(step_beta_shapes([0.0, 0.0], weight, a, b), (a, b));
        // A slope pointing at (100, 100), past the upper bound of both shapes, stops at it.
        let (far_a, far_b) = (100.0, 100.0);
        let towards_far = [
            h11 * (far_a - a) + h12 * (far_b - b),
            h12 * (far_a - a) + h22 * (far_b - b),
        ];
        assert_eq!(
            step_beta_shapes(towards_far, weight, a, b),
            (BETA_SHAPE_BOUNDS.1, BETA_SHAPE_BOUNDS.1)
        );
        assert_eq!(step_beta_shapes([f64::NAN, 1.0], weight, a, b), (a, b));
        // A count of positions so small that the curvature's determinant underflows to zero.
        assert_eq!(step_beta_shapes([1.0, -1.0], 1e-170, a, b), (a, b));
        // One step worked by hand rather than from the code's own trigamma: at (1, 1), one position's
        // curvature is [[1, −c], [−c, 1]] with c = ψ′(2) = π²/6 − 1, so a slope of (1, 0) moves the
        // shapes by (1, c) / (1 − c²).
        let c = std::f64::consts::PI * std::f64::consts::PI / 6.0 - 1.0;
        let (by_hand_a, by_hand_b) = step_beta_shapes([1.0, 0.0], 1.0, 1.0, 1.0);
        assert!(
            (by_hand_a - (1.0 + 1.0 / (1.0 - c * c))).abs() < 1e-6
                && (by_hand_b - (1.0 + c / (1.0 - c * c))).abs() < 1e-6,
            "from (1, 1) the step reached ({by_hand_a}, {by_hand_b})"
        );

        let (a, b) = (0.7, BETA_SHAPE_BOUNDS.1);
        let ab = a + b;
        let h11 = weight * (trigamma(a) - trigamma(ab));
        for slope_in_a in [0.0, 30.0, -12.0] {
            let (stepped_a, stepped_b) = step_beta_shapes([slope_in_a, 500.0], weight, a, b);
            assert_eq!(
                stepped_b, BETA_SHAPE_BOUNDS.1,
                "the second shape held at its bound"
            );
            assert!(
                (stepped_a - (a + slope_in_a / h11)).abs() < 1e-12,
                "slope {slope_in_a} in the first shape: stepped to {stepped_a}"
            );
        }
    }
}
