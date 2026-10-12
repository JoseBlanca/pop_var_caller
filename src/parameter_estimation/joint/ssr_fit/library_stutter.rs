//! **Each library's stutter multiplier, estimated with the repeat classes' other numbers held**
//! (`doc/devel/ng/spec/str_library_stutter.md` §3.2–3.4) — a prototype.
//!
//! Each read a library put on a tract either slipped or reported its allele's length. Given what a
//! class's fit says about the tract — its allele-length frequencies, how a slipped read is spread,
//! each sample's homozygote excess — every read is shared between those two outcomes by how likely
//! each is ([`expected_slips`]). Summed over a library's tracts in a class, that gives the expected
//! slipped reads `S` and exactly-read reads `N` of the spec's §3.2, and from them:
//!
//! - **each library's multiplier** with the class rates held ([`fit_multiplier`]);
//! - **each class's rate** with the multipliers held ([`fit_class_rate`]) — the prototype's stand-in
//!   for refitting the whole class, which takes the climb an hour on a 50-sample cohort.
//!
//! Alternating the two, each time recounting `S` and `N` at the new numbers, is expectation
//! maximisation in the class rates and the multipliers together, with every other number of a class
//! held where the pooled fit left it ([`fit_rates_and_multipliers`]).

use super::*;

/// **The numbers of one repeat class held while rates and multipliers are fitted**: the class's
/// slippage numbers for each slippage group, its allele-length frequencies and their spread.
#[derive(Debug, Clone, PartialEq)]
pub struct HeldClass {
    /// Per slippage group, the class's three slippage numbers. The level is the class rate; the
    /// direction and fall-off are held.
    pub slippage: Vec<Slippage>,
    /// How the class's chromosomes are spread over the allele lengths, `-allele_span ..= +allele_span`.
    pub allele_length_variation: Vec<f64>,
    /// How monomorphic the class's tracts are.
    pub concentration: f64,
}

/// **One class's expected slipped and exactly-read reads, for each library** of
/// [`StratumEvidence::libraries`] — zero for a library with no read in the class.
#[derive(Debug, Clone, PartialEq)]
pub struct SlipCounts {
    pub slipped: Vec<f64>,
    pub exact: Vec<f64>,
}

/// **Shares every read of one class between "slipped" and "read exactly"**, by how likely each is
/// under each sample's genotype posterior at the tract (spec §3.2).
///
/// The genotype posterior is the class's own likelihood's: the sample's reads under each genotype,
/// times the genotype's probability under the tract's allele-frequency integral and the sample's
/// homozygote excess — the same integral the climb scores, taken over the same points.
#[must_use]
pub fn expected_slips(
    evidence: &StratumEvidence,
    held: &HeldClass,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
) -> SlipCounts {
    let classes = allele_classes_of(config);
    assert_eq!(
        held.allele_length_variation.len(),
        classes,
        "period {} at {} repeats: an allele length variation of {} classes held against an allele span of {}",
        evidence.stratum.period,
        evidence.stratum.reference_repeats,
        held.allele_length_variation.len(),
        config.allele_span
    );
    let genotypes = genotype_pairs(classes);
    let distributions = ReadDistributionsByLevel::of(
        &held.slippage,
        &evidence.libraries,
        config.allele_span,
        evidence.read_span,
    );
    let quadrature = dirichlet_points(
        &held.allele_length_variation,
        held.concentration,
        config.quadrature_points,
        &genotypes,
        WhereTheThreadsGo::AcrossStrata,
    );
    // Per library, the probability a read of it reports its allele's length: one minus its level.
    let exact_share: Vec<f64> = evidence
        .libraries
        .iter()
        .map(|library| {
            1.0 - library_level(
                held.slippage[library.slippage_group as usize].level,
                library.stutter_multiplier,
            )
        })
        .collect();
    let slot_of = |class: usize| {
        let offset = class as i32 - config.allele_span;
        (offset.clamp(-evidence.read_span, evidence.read_span) + evidence.read_span) as usize
    };
    let slots: Vec<usize> = (0..classes).map(slot_of).collect();

    let per_tract: Vec<SlipCounts> = evidence
        .tracts
        .iter()
        .map(|tract| {
            let mut counts = SlipCounts {
                slipped: vec![0.0; evidence.libraries.len()],
                exact: vec![0.0; evidence.libraries.len()],
            };
            if tract.samples.is_empty() {
                return counts;
            }
            let likelihoods =
                TractLikelihoods::of(tract, &distributions, &genotypes, homozygote_excess);
            let posteriors = genotype_posteriors(&likelihoods, &quadrature, &genotypes);
            for (row, sample) in tract.samples.iter().enumerate() {
                let posterior = &posteriors[row * genotypes.len()..][..genotypes.len()];
                for (library, reads) in &sample.by_library {
                    let per_allele = distributions.of_library(*library);
                    let stay = exact_share[library.get()];
                    for (bucket, count) in reads.iter().enumerate() {
                        if *count == 0 {
                            continue;
                        }
                        // P(slipped | this bucket), averaged over the sample's genotypes.
                        let mut slipped_share = 0.0;
                        for ((first, second), weight) in genotypes.iter().zip(posterior) {
                            if *weight < 1e-12 {
                                continue;
                            }
                            let reported =
                                0.5 * (per_allele[*first][bucket] + per_allele[*second][bucket]);
                            let unslipped = 0.5
                                * stay
                                * (f64::from(u8::from(slots[*first] == bucket))
                                    + f64::from(u8::from(slots[*second] == bucket)));
                            if reported > 0.0 {
                                slipped_share += weight * (1.0 - unslipped / reported).max(0.0);
                            }
                        }
                        let count = f64::from(*count);
                        counts.slipped[library.get()] += count * slipped_share;
                        counts.exact[library.get()] += count * (1.0 - slipped_share);
                    }
                }
            }
            counts
        })
        .collect();
    let mut total = SlipCounts {
        slipped: vec![0.0; evidence.libraries.len()],
        exact: vec![0.0; evidence.libraries.len()],
    };
    for counts in per_tract {
        for (into, from) in total.slipped.iter_mut().zip(&counts.slipped) {
            *into += from;
        }
        for (into, from) in total.exact.iter_mut().zip(&counts.exact) {
            *into += from;
        }
    }
    total
}

/// **Each sample's posterior over the genotypes at one tract**, a row a sample with reads, laid out
/// as [`TractLikelihoods::scaled`] is.
///
/// At each point of the integral, a sample's genotype is weighted by its prior there times its
/// reads' likelihood; the points are weighted by how well they explain every sample at once.
fn genotype_posteriors(
    likelihoods: &TractLikelihoods,
    quadrature: &Quadrature,
    genotypes: &[(usize, usize)],
) -> Vec<f64> {
    let width = genotypes.len();
    let rows = likelihoods.homozygote_excess.len();
    let points = quadrature.frequencies.len() / quadrature.classes;
    // Per point, each row's total over genotypes, and the point's log weight.
    let mut sums = vec![0.0; points * rows];
    let mut ln_weights = vec![f64::NEG_INFINITY; points];
    for (index, point) in quadrature
        .frequencies
        .chunks_exact(quadrature.classes)
        .enumerate()
    {
        let independent = &quadrature.independent[index * width..][..width];
        let mut ln_weight = quadrature.ln_weight;
        for (row, (values, excess)) in likelihoods
            .scaled
            .chunks_exact(width)
            .zip(&likelihoods.homozygote_excess)
            .enumerate()
        {
            let at_random: f64 = independent.iter().zip(values).map(|(w, v)| w * v).sum();
            let by_descent: f64 = quadrature
                .diagonal
                .iter()
                .enumerate()
                .map(|(class, slot)| point[class] * values[*slot])
                .sum();
            let sum = excess * by_descent + (1.0 - excess) * at_random;
            sums[index * rows + row] = sum;
            ln_weight += float::ln(sum.max(1e-300));
        }
        ln_weights[index] = ln_weight;
    }
    let largest = ln_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut point_weights: Vec<f64> = ln_weights.iter().map(|w| float::exp(w - largest)).collect();
    normalise(&mut point_weights);

    let mut posteriors = vec![0.0; rows * width];
    for (index, point) in quadrature
        .frequencies
        .chunks_exact(quadrature.classes)
        .enumerate()
    {
        let point_weight = point_weights[index];
        if point_weight < 1e-15 {
            continue;
        }
        let independent = &quadrature.independent[index * width..][..width];
        for (row, (values, excess)) in likelihoods
            .scaled
            .chunks_exact(width)
            .zip(&likelihoods.homozygote_excess)
            .enumerate()
        {
            let sum = sums[index * rows + row];
            if sum <= 0.0 {
                continue;
            }
            let into = &mut posteriors[row * width..][..width];
            for (slot, value) in values.iter().enumerate() {
                into[slot] += point_weight * (1.0 - excess) * independent[slot] * value / sum;
            }
            for (class, slot) in quadrature.diagonal.iter().enumerate() {
                into[*slot] += point_weight * excess * point[class] * values[*slot] / sum;
            }
        }
    }
    posteriors
}

/// **The multiplier one library's counts say, with every class rate held** (spec §3.2–3.4), and its
/// standard error.
///
/// `by_class` holds, for each class the library read in, its expected slipped reads, exactly-read
/// reads and the class rate. The objective is `Σ [S ln(m L) + N ln(1 − m L)] + P (ln m − m)` — the
/// second term `pull` pseudo-reads that hold thin evidence near one — maximised by Newton's method
/// inside `(0, MULTIPLIED_LEVEL_CEILING / max L)`.
#[must_use]
pub fn fit_multiplier(by_class: &[(f64, f64, f64)], pull: f64) -> (f64, f64) {
    let largest = by_class
        .iter()
        .map(|(_, _, rate)| *rate)
        .fold(0.0, f64::max);
    let upper = if largest > 0.0 {
        MULTIPLIED_LEVEL_CEILING / largest
    } else {
        f64::INFINITY
    };
    let slope_and_curvature = |m: f64| {
        let mut slope = pull / m - pull;
        let mut curvature = -pull / (m * m);
        for (slipped, exact, rate) in by_class {
            slope += slipped / m - exact * rate / (1.0 - m * rate);
            curvature += -slipped / (m * m) - exact * rate * rate / float::powi(1.0 - m * rate, 2);
        }
        (slope, curvature)
    };
    let (mut low, mut high) = (1e-6_f64, upper.min(1e6) * (1.0 - 1e-9));
    let mut m = 1.0_f64.clamp(low, high);
    for _ in 0..100 {
        let (slope, curvature) = slope_and_curvature(m);
        if slope > 0.0 {
            low = m;
        } else {
            high = m;
        }
        let newton = m - slope / curvature;
        let next = if curvature < 0.0 && newton > low && newton < high {
            newton
        } else {
            0.5 * (low + high)
        };
        if (next - m).abs() <= 1e-12 * m {
            m = next;
            break;
        }
        m = next;
    }
    let (_, curvature) = slope_and_curvature(m);
    (m, (-1.0 / curvature).sqrt())
}

/// **The rate one class's counts say, with every library's multiplier held**: the `L` maximising
/// `Σ_l [S_l ln(m_l L) + N_l ln(1 − m_l L)]`, from `(S_l, N_l, m_l)` for each library that read in
/// the class. `None` where no read slipped or none read exactly.
#[must_use]
pub fn fit_class_rate(by_library: &[(f64, f64, f64)]) -> Option<f64> {
    let slipped: f64 = by_library.iter().map(|(s, _, _)| s).sum();
    let exact: f64 = by_library.iter().map(|(_, n, _)| n).sum();
    if slipped <= 0.0 || exact <= 0.0 {
        return None;
    }
    let largest = by_library.iter().map(|(_, _, m)| *m).fold(0.0, f64::max);
    let (mut low, mut high) = (1e-12_f64, (1.0 / largest) * (1.0 - 1e-9));
    let slope_and_curvature = |rate: f64| {
        let mut slope = slipped / rate;
        let mut curvature = -slipped / (rate * rate);
        for (_, n, m) in by_library {
            slope -= n * m / (1.0 - m * rate);
            curvature -= n * m * m / float::powi(1.0 - m * rate, 2);
        }
        (slope, curvature)
    };
    let mut rate = (slipped / (slipped + exact)).clamp(low, high);
    for _ in 0..100 {
        let (slope, curvature) = slope_and_curvature(rate);
        if slope > 0.0 {
            low = rate;
        } else {
            high = rate;
        }
        let newton = rate - slope / curvature;
        let next = if newton > low && newton < high {
            newton
        } else {
            0.5 * (low + high)
        };
        if (next - rate).abs() <= 1e-14 * rate {
            return Some(next);
        }
        rate = next;
    }
    Some(rate)
}

/// What [`fit_rates_and_multipliers`] returns.
#[derive(Debug, Clone, PartialEq)]
pub struct RatesAndMultipliers {
    /// Per library, its multiplier, scaled so the read-weighted geometric mean is one (spec §2.2).
    pub multipliers: Vec<f64>,
    /// Per library, the multiplier's standard error, on the same scale.
    pub standard_errors: Vec<f64>,
    /// Per library, its expected slipped reads over every class, at the answer.
    pub slipped_reads: Vec<f64>,
    /// Per library, its reads over every class — the weights of the scaling.
    pub reads: Vec<f64>,
    /// Per class, in the order the strata were handed in, the rate of slippage group 0.
    pub class_rates: Vec<f64>,
    /// How many rounds ran, and the largest move of any multiplier in its own errors in the last.
    pub rounds: usize,
    pub last_move_in_errors: f64,
}

/// **The class rates and the libraries' multipliers fitted together** (spec §3.1), every other
/// number of a class held at `held`: rounds of counting the expected slips at the current numbers,
/// then each multiplier with the rates held, the scaling, then each rate with the multipliers held,
/// until no multiplier moves by more than `SETTLED_FRACTION` of its error or `max_rounds` pass.
///
/// **One slippage group only** — the prototype's cohorts have one.
///
/// # Panics
///
/// When `strata` and `held` differ in length, or a class names more than one slippage group.
#[must_use]
pub fn fit_rates_and_multipliers(
    strata: &mut [StratumEvidence],
    held: &mut [HeldClass],
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    pull: f64,
    max_rounds: usize,
    mut report: impl FnMut(usize, &[f64], f64),
) -> RatesAndMultipliers {
    assert_eq!(strata.len(), held.len());
    assert!(held.iter().all(|class| class.slippage.len() == 1));
    let libraries = strata[0].libraries.len();
    let mut multipliers = vec![1.0; libraries];
    let mut errors = vec![f64::INFINITY; libraries];
    let mut slipped_reads = vec![0.0; libraries];
    let mut read_weights = vec![0.0; libraries];
    let mut rounds = 0;
    let mut last_move = f64::INFINITY;
    while rounds < max_rounds {
        rounds += 1;
        let list: Arc<[Library]> = strata[0]
            .libraries
            .iter()
            .zip(&multipliers)
            .map(|(library, m)| Library {
                stutter_multiplier: *m,
                ..*library
            })
            .collect();
        for evidence in strata.iter_mut() {
            evidence.libraries = Arc::clone(&list);
        }
        let counts: Vec<SlipCounts> = strata
            .par_iter()
            .zip(held.par_iter())
            .map(|(evidence, class)| expected_slips(evidence, class, homozygote_excess, config))
            .collect();

        // Each multiplier, the rates held.
        let mut fitted = vec![1.0; libraries];
        let mut fitted_errors = vec![f64::INFINITY; libraries];
        let mut reads = vec![0.0; libraries];
        for library in 0..libraries {
            let by_class: Vec<(f64, f64, f64)> = counts
                .iter()
                .zip(held.iter())
                .filter(|(count, _)| count.slipped[library] + count.exact[library] > 0.0)
                .map(|(count, class)| {
                    (
                        count.slipped[library],
                        count.exact[library],
                        class.slippage[0].level,
                    )
                })
                .collect();
            reads[library] = by_class.iter().map(|(s, n, _)| s + n).sum();
            slipped_reads[library] = by_class.iter().map(|(s, _, _)| s).sum();
            if by_class.is_empty() {
                continue;
            }
            let (m, error) = fit_multiplier(&by_class, pull);
            fitted[library] = m;
            fitted_errors[library] = error;
        }
        // The scale: a read-weighted geometric mean of one, the rates taking the factor.
        let total: f64 = reads.iter().sum();
        let ln_mean: f64 = fitted
            .iter()
            .zip(&reads)
            .map(|(m, r)| r * float::ln(*m))
            .sum::<f64>()
            / total.max(1e-300);
        let scale = float::exp(ln_mean);
        for (m, e) in fitted.iter_mut().zip(fitted_errors.iter_mut()) {
            *m /= scale;
            *e /= scale;
        }
        // Each rate, the multipliers held — counts from this round, multipliers just fitted.
        for (count, class) in counts.iter().zip(held.iter_mut()) {
            let by_library: Vec<(f64, f64, f64)> = (0..libraries)
                .filter(|l| count.slipped[*l] + count.exact[*l] > 0.0)
                .map(|l| (count.slipped[l], count.exact[l], fitted[l]))
                .collect();
            if let Some(rate) = fit_class_rate(&by_library) {
                class.slippage[0].level = rate;
            }
        }
        last_move = multipliers
            .iter()
            .zip(&fitted)
            .zip(&fitted_errors)
            .map(|((old, new), error)| (new - old).abs() / error)
            .fold(0.0, f64::max);
        multipliers = fitted;
        errors = fitted_errors;
        read_weights = reads;
        report(rounds, &multipliers, last_move);
        if last_move < SETTLED_FRACTION {
            break;
        }
    }
    RatesAndMultipliers {
        multipliers,
        standard_errors: errors,
        slipped_reads,
        reads: read_weights,
        class_rates: held.iter().map(|class| class.slippage[0].level).collect(),
        rounds,
        last_move_in_errors: last_move,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameter_estimation::joint::ssr_fit::bench_fixtures::{draw_stratum, spectrum_of};

    /// **A library's multiplier is its own slip share over the class rate when it reads one class,
    /// and the pull holds a library with no evidence at one.**
    #[test]
    fn the_multiplier_is_the_librarys_own_ratio_and_one_without_reads() {
        let (m, _) = fit_multiplier(&[(30.0, 970.0, 0.015)], 0.0);
        assert!((m * 0.015 - 0.03).abs() < 1e-9, "{m}");
        let (m, _) = fit_multiplier(&[(0.0, 0.0, 0.015)], 20.0);
        assert!((m - 1.0).abs() < 1e-9, "{m}");
    }

    /// **The class rate is the pooled share where every multiplier is one.**
    #[test]
    fn the_class_rate_is_the_pooled_share_at_multipliers_of_one() {
        let rate = fit_class_rate(&[(10.0, 90.0, 1.0), (30.0, 170.0, 1.0)]).expect("a rate");
        assert!((rate - 40.0 / 300.0).abs() < 1e-12, "{rate}");
    }

    /// **The expected slips over a drawn class recover its level**: the slipped share of every read,
    /// counted under the true numbers, is the level the reads were drawn at, within sampling error.
    #[test]
    fn the_expected_slips_at_the_truth_give_the_drawn_level() {
        let truth = Slippage {
            level: 0.08,
            shorter_share: 0.8,
            fall_off: 0.25,
        };
        let config = SsrFitConfig {
            allele_span: 1,
            ..SsrFitConfig::default()
        };
        let spectrum = spectrum_of(3);
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 300, 20, 10, 1, 11);
        let held = HeldClass {
            slippage: vec![truth],
            allele_length_variation: spectrum,
            concentration: 0.5,
        };
        let counts = expected_slips(&evidence, &held, &[0.4; 20], &config);
        let slipped: f64 = counts.slipped.iter().sum();
        let exact: f64 = counts.exact.iter().sum();
        let share = slipped / (slipped + exact);
        assert!((share - 0.08).abs() < 0.01, "slipped share {share}");
    }
}
