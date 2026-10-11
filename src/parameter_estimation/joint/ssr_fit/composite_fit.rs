//! **A repeat class's stutter and allele-length shares by expectation maximisation, without the
//! average over tract frequencies** — the owner's routes 1 and 2 together (2026-10-11), a prototype.
//!
//! The fit proper averages each tract's likelihood over how its alleles could be shared among its
//! chromosomes, and that average, over 256 points, is several log-likelihood units a tract wrong at
//! thirteen allele lengths (`fit_precision_quadrature_average_2026-10-01.md`). It also makes the fit
//! slow. This estimates the same stutter and the class's allele-length shares from a simpler
//! objective: a sample's alleles at a random tract of the class are drawn from the class's shares,
//! whatever that tract's own frequencies, so
//!
//! `Σ over (tract, sample) of ln Σ_g prior(g | shares, F) L(g | stutter)`
//!
//! is maximised by the class's numbers. It ignores that samples at one tract share that tract's
//! frequencies — which costs precision, not consistency — and so says nothing about how much tracts
//! differ (the concentration).
//!
//! Each round shares every read between "read exactly" and "slipped this far, this way", given each
//! sample's genotype probabilities, and every sample's alleles among the lengths; the slip rate, the
//! shorter share, the fall-off and the shares then follow from the expected counts.

use super::*;

/// What [`fit_class_by_counting`] returns.
#[derive(Debug, Clone, PartialEq)]
pub struct CountedClass {
    /// The class's stutter.
    pub slippage: Slippage,
    /// How the class's chromosomes are spread over the allele lengths, `-allele_span ..= +allele_span`.
    pub length_spectrum: Vec<f64>,
    /// How much the class's tracts differ from one another: the Dirichlet concentration, read off
    /// how much more often a sample's two copies match than its homozygote excess alone explains.
    pub concentration: f64,
    /// How many rounds ran, and the largest change of any number in the last.
    pub rounds: usize,
    pub last_change: f64,
    /// The objective at the answer.
    pub composite_log_likelihood: f64,
}

/// The expected counts one round gathers over a class.
#[derive(Debug, Clone, Default)]
struct Counts {
    /// Per allele class, how many independent allele draws.
    draws: Vec<f64>,
    /// Of each sample's two copies not identical by descent, how many pairs, and how many of those
    /// were one draw from the tract's frequencies.
    pairs_not_by_descent: f64,
    shared_in_tract: f64,
    exact: f64,
    slipped: f64,
    /// Slipped reads whose bucket is not an end bucket, by direction, and their steps.
    shorter: f64,
    longer: f64,
    steps: f64,
    stepped: f64,
    log_likelihood: f64,
}

impl Counts {
    fn new(classes: usize) -> Self {
        Self {
            draws: vec![0.0; classes],
            ..Self::default()
        }
    }

    fn add(&mut self, other: &Counts) {
        for (into, from) in self.draws.iter_mut().zip(&other.draws) {
            *into += from;
        }
        self.exact += other.exact;
        self.pairs_not_by_descent += other.pairs_not_by_descent;
        self.shared_in_tract += other.shared_in_tract;
        self.slipped += other.slipped;
        self.shorter += other.shorter;
        self.longer += other.longer;
        self.steps += other.steps;
        self.stepped += other.stepped;
        self.log_likelihood += other.log_likelihood;
    }
}

/// **Fits one class's stutter and shares by counting**, from `start`'s stutter and flat shares,
/// until no number moves by more than `tolerance` or `max_rounds` pass. Every library is scored at
/// the class's stutter (multipliers are not fitted here).
#[must_use]
pub fn fit_class_by_counting(
    evidence: &StratumEvidence,
    start: Slippage,
    homozygote_excess: &[f64],
    config: &SsrFitConfig,
    max_rounds: usize,
    tolerance: f64,
) -> CountedClass {
    let classes = allele_classes_of(config);
    let genotypes = genotype_pairs(classes);
    let span = config.allele_span;
    let read_span = evidence.read_span;
    let last_bucket = (2 * read_span) as usize;
    let slot =
        |class: usize| ((class as i32 - span).clamp(-read_span, read_span) + read_span) as usize;
    let mut slippage = start;
    let mut shares = vec![1.0 / classes as f64; classes];
    let mut concentration = 1.0_f64;
    let mut rounds = 0;
    let mut last_change = f64::INFINITY;
    let mut log_likelihood = f64::NEG_INFINITY;
    let libraries: Vec<Library> = evidence
        .libraries
        .iter()
        .map(|library| Library {
            slippage_group: 0,
            stutter_multiplier: 1.0,
            ..*library
        })
        .collect();
    while rounds < max_rounds {
        rounds += 1;
        let distributions = ReadDistributionsByLevel::of(&[slippage], &libraries, span, read_span);
        let per_allele = &distributions.read_distribution_by_allele[0];
        let stay = 1.0 - slippage.level;
        let censored_extra = slippage.fall_off / (1.0 - slippage.fall_off);
        // Each genotype's prior, given a sample's homozygote excess, split into the part where its two
        // copies are one draw (identical by descent) and the part where they are two.
        let one_draw: Vec<f64> = genotypes
            .iter()
            .map(|(a, b)| if a == b { shares[*a] } else { 0.0 })
            .collect();
        let two_draws: Vec<f64> = genotypes
            .iter()
            .map(|(a, b)| {
                if a == b {
                    shares[*a] * shares[*a]
                } else {
                    2.0 * shares[*a] * shares[*b]
                }
            })
            .collect();
        let per_tract: Vec<Counts> = evidence
            .tracts
            .par_iter()
            .map(|tract| {
                let mut counts = Counts::new(classes);
                for sample in &tract.samples {
                    let excess = homozygote_excess[sample.sample as usize];
                    // Pool the sample's libraries: one stutter for all of them here.
                    let mut reads = vec![0_u32; last_bucket + 1];
                    for (_, library_reads) in &sample.by_library {
                        for (into, count) in reads.iter_mut().zip(library_reads) {
                            *into += count;
                        }
                    }
                    // Each genotype's log-likelihood, then its posterior.
                    let log_likelihoods: Vec<f64> = genotypes
                        .iter()
                        .map(|(a, b)| {
                            let mut total = 0.0;
                            add_ln_reads(&mut total, per_allele, &reads, *a, *b);
                            total
                        })
                        .collect();
                    let largest = log_likelihoods
                        .iter()
                        .copied()
                        .fold(f64::NEG_INFINITY, f64::max);
                    // A sample's two copies are one draw if identical by descent, or, failing that, if they
                    // drew the same chromosome lineage from the tract's own frequencies: chance
                    // 1 / (1 + concentration) under a Dirichlet.
                    let from_the_tract = (1.0 - excess) / (1.0 + concentration);
                    let one_draw_chance = excess + from_the_tract;
                    let two_draw_chance = (1.0 - excess) * concentration / (1.0 + concentration);
                    let joint_one: Vec<f64> = log_likelihoods
                        .iter()
                        .zip(&one_draw)
                        .map(|(ll, prior)| one_draw_chance * prior * float::exp(ll - largest))
                        .collect();
                    let joint_two: Vec<f64> = log_likelihoods
                        .iter()
                        .zip(&two_draws)
                        .map(|(ll, prior)| two_draw_chance * prior * float::exp(ll - largest))
                        .collect();
                    let total: f64 = joint_one.iter().sum::<f64>() + joint_two.iter().sum::<f64>();
                    if total <= 0.0 {
                        continue;
                    }
                    counts.log_likelihood += float::ln(total) + largest;
                    for (index, (a, b)) in genotypes.iter().enumerate() {
                        let one = joint_one[index] / total;
                        let two = joint_two[index] / total;
                        let weight = one + two;
                        if weight < 1e-12 {
                            continue;
                        }
                        counts.draws[*a] += one + two;
                        counts.draws[*b] += two;
                        let shared = one * from_the_tract / one_draw_chance.max(1e-300);
                        counts.shared_in_tract += shared;
                        counts.pairs_not_by_descent += shared + two;
                        for (bucket, count) in reads.iter().enumerate() {
                            if *count == 0 {
                                continue;
                            }
                            let reported = 0.5 * (per_allele[*a][bucket] + per_allele[*b][bucket]);
                            if reported <= 0.0 {
                                continue;
                            }
                            for allele in [*a, *b] {
                                let from_allele = 0.5 * per_allele[allele][bucket] / reported;
                                let share = weight * f64::from(*count) * from_allele;
                                let unslipped = if slot(allele) == bucket {
                                    stay / per_allele[allele][bucket]
                                } else {
                                    0.0
                                };
                                let exact = share * unslipped.min(1.0);
                                let slipped = share - exact;
                                counts.exact += exact;
                                counts.slipped += slipped;
                                // Which way a slip went is read from any bucket; how far, only
                                // from one the recording did not saturate.
                                let from = slot(allele);
                                if slipped > 0.0 && bucket < from {
                                    counts.shorter += slipped;
                                } else if slipped > 0.0 && bucket > from {
                                    counts.longer += slipped;
                                }
                                // How far: a slip into an end bucket went at least that far, and a
                                // geometric fall-off has no memory, so its expected length is the
                                // distance plus f / (1 − f) at the current fall-off.
                                if slipped > 0.0 && bucket != from {
                                    let distance = bucket.abs_diff(from) as f64;
                                    let censored = bucket == 0 || bucket == last_bucket;
                                    let expected = if censored {
                                        distance + censored_extra
                                    } else {
                                        distance
                                    };
                                    counts.steps += slipped * expected;
                                    counts.stepped += slipped;
                                }
                            }
                        }
                    }
                }
                counts
            })
            .collect();
        let mut counts = Counts::new(classes);
        for tract in &per_tract {
            counts.add(tract);
        }
        log_likelihood = counts.log_likelihood;
        let draws: f64 = counts.draws.iter().sum();
        let new_shares: Vec<f64> = counts.draws.iter().map(|d| d / draws).collect();
        let new_concentration = if counts.shared_in_tract > 0.0 {
            (counts.pairs_not_by_descent / counts.shared_in_tract - 1.0).clamp(1e-3, 1e4)
        } else {
            1e4
        };
        let new_slippage = Slippage {
            level: counts.slipped / (counts.slipped + counts.exact),
            shorter_share: if counts.shorter + counts.longer > 0.0 {
                counts.shorter / (counts.shorter + counts.longer)
            } else {
                slippage.shorter_share
            },
            // A step k has probability (1 − f) f^(k − 1), whose mean is 1 / (1 − f).
            fall_off: if counts.steps > 0.0 {
                (1.0 - counts.stepped / counts.steps).clamp(0.0, 0.99)
            } else {
                slippage.fall_off
            },
        };
        last_change = new_shares
            .iter()
            .zip(&shares)
            .map(|(new, old)| (new - old).abs())
            .chain([
                (new_slippage.level - slippage.level).abs(),
                (new_slippage.shorter_share - slippage.shorter_share).abs(),
                (new_slippage.fall_off - slippage.fall_off).abs(),
                (float::ln(new_concentration) - float::ln(concentration)).abs() * 0.01,
            ])
            .fold(0.0, f64::max);
        shares = new_shares;
        slippage = new_slippage;
        concentration = new_concentration;
        if last_change < tolerance {
            break;
        }
    }
    CountedClass {
        slippage,
        length_spectrum: shares,
        concentration,
        rounds,
        last_change,
        composite_log_likelihood: log_likelihood,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parameter_estimation::joint::ssr_fit::bench_fixtures::{draw_stratum, spectrum_of};

    /// **On a drawn class the counting fit returns the stutter and the shares it was drawn at**,
    /// within sampling error, from a start far from them. At a concentration of 0.5 each tract is nearly fixed
    /// at one length, so 2,000 tracts pin a share to about ±0.011.
    #[test]
    fn counting_recovers_a_drawn_class() {
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
        let evidence = draw_stratum(truth, &spectrum, 0.5, 0.4, 2000, 20, 10, 1, 13);
        let start = Slippage {
            level: 0.2,
            shorter_share: 0.5,
            fall_off: 0.5,
        };
        let fitted = fit_class_by_counting(&evidence, start, &[0.4; 20], &config, 500, 1e-9);
        assert!((fitted.slippage.level - 0.08).abs() < 0.01, "{fitted:?}");
        assert!(
            (fitted.slippage.shorter_share - 0.8).abs() < 0.05,
            "{fitted:?}"
        );
        for (fitted, truth) in fitted.length_spectrum.iter().zip(&spectrum) {
            assert!((fitted - truth).abs() < 0.035, "{fitted} against {truth}");
        }
        assert!(
            (0.3..0.8).contains(&fitted.concentration),
            "concentration {} against 0.5",
            fitted.concentration
        );
    }
}
