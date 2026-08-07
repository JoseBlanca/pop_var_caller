//! The coupled fit: each library's error rate and the individual's genotype frequencies,
//! fitted together by alternating between the two tables that hold them.
//!
//! **Why they cannot be fitted apart.** A higher error rate explains the same alternative
//! reads as less real variation, so `ε` and the genotype frequencies trade off inside one
//! likelihood — and they are read off two different tables, the rates from the read-group
//! one and the frequencies from the whole-sample one
//! (`spec/parameter_prepass_generic.md` §5.1).
//!
//! **What one iteration is**, and it is the alternation the research harness measured
//! rather than the one an earlier draft of the architecture described:
//!
//! 1. **The frequencies**, climbed on the whole-sample table at the rates the previous
//!    iteration produced — one set per ploidy present, because a haploid region has two
//!    genotype classes and a diploid three.
//! 2. **Each read group's rate**, scanned over the error-rate ladder on that group's own
//!    table, at the frequencies step 1 just produced and **without re-climbing them**.
//!    That is `read_group_error_rate::fit_read_group_error_rates`, and it is
//!    `examples/ng_multilib_key_harness.rs`'s `fit_eps_on_read_group(space, freqs)`.
//!
//! **It stops when every read group's winning rung is the one it had last iteration.** The
//! scan returns a rung index, so "moves by less than one rung" and "does not move" are the
//! same condition and only the second is testable — a loop oscillating between two adjacent
//! rungs would satisfy a movement tolerance forever. Rung stability is also the right
//! resolution: worlds that ran the harness's full 200 iterations without meeting a movement
//! tolerance of 10⁻¹² were already at the truth to better than a thousandth of a rung
//! (research note §2.6).
//!
//! **This is a fixed point of two estimating equations, not a climb on one objective**, so
//! it is capped, the best-scoring iterate is kept rather than the last, and how it ended is
//! reported. That it lands on the truth at all is measured rather than argued: from a start
//! at three times the true rates and half the true frequencies, the fixed point is the
//! truth in all 25 of the harness's worlds — 0.000 rungs on every error rate and 0.000% on
//! both frequencies (research note §2.6).
//!
//! **How much it can matter.** At one read group the two tables are the same table
//! (`spec/parameter_prepass_generic.md` §12.6), so the alternation is plain coordinate
//! ascent on a single objective and reaches the same joint maximum the profile scan
//! returns. That is 1,550 of the 1,707 samples in the tomato archive survey; the coupling
//! bites only on the 157 multi-library ones, and on neither cohort in hand.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §5.2, spec §5.1.

use std::collections::{BTreeMap, BTreeSet};

use smallvec::SmallVec;

use crate::ng::parameter_estimation::fitting::ladder_scan::fit_by_fixed_frequency_scan;
use crate::ng::parameter_estimation::fitting::mixture_weights::{
    GenotypeLikelihoodTable, fit_mixture_weights,
};
use crate::ng::parameter_estimation::fitting::{FitTermination, NoiseModel};
use crate::ng::parameter_estimation::generic::accumulators::GenericAccumulators;
use crate::ng::parameter_estimation::generic::histogram::{Cell, DepthAltHistogram};
use crate::ng::parameter_estimation::generic::noise_model::{
    LibraryNoise, SampleLibraryNoise, SubstitutionNoiseModel,
};
use crate::ng::parameter_estimation::generic::read_group_error_rate::fit_read_group_error_rates;
use crate::ng::parameter_estimation::generic::{
    CoupledFit, DEFAULT_ERROR_RATE, MAX_COUPLED_FIT_ITERATIONS, MIN_SITES_TO_FIT, SampleRates,
};
use crate::ng::parameter_estimation::{Estimate, ParameterEstimationError, Provenance};
use crate::ng::types::{ErrorRate, GenotypeFrequency, LogProb, Ploidy, ReadGroupId};

/// Fit a sample's error rates and genotype frequencies together, from its accumulators.
///
/// The thin door: it pulls the two tables out and hands them to
/// [`fit_coupled_from_tables`], which is where the alternation is and what the tests drive.
/// A fit that could only be reached through an accumulator could only be tested through a
/// locus stream, and this plan's rule is that the fits are proven before a locus is read.
///
/// # Errors
///
/// [`ParameterEstimationError::GenotypeFrequenciesNotFittable`] when some ploidy's
/// whole-sample table holds fewer than [`MIN_SITES_TO_FIT`] sites. There is no fallback:
/// a sample has one heterozygosity, so there is no sibling to borrow it from and no
/// constant worth inventing, because it is the biology.
///
/// # Panics
///
/// If `ladder` is empty, or if the sample has no read group with reads.
pub fn fit_coupled(
    sample: &str,
    accumulators: &GenericAccumulators,
    ladder: &[ErrorRate],
) -> Result<CoupledFit, ParameterEstimationError> {
    let ploidies: BTreeSet<Ploidy> = accumulators
        .windowed_histograms()
        .keys()
        .map(|&(_, ploidy)| ploidy)
        .collect();
    let whole_sample: BTreeMap<Ploidy, DepthAltHistogram<u64>> = ploidies
        .into_iter()
        .map(|ploidy| (ploidy, accumulators.whole_sample_histogram(ploidy)))
        .collect();

    fit_coupled_from_tables(
        sample,
        accumulators.read_group_histograms(),
        &whole_sample,
        ladder,
    )
}

/// [`fit_coupled`] over tables handed in directly, which is how it is proven.
///
/// `whole_sample` is one table per ploidy — the windows of that ploidy folded together,
/// which is what [`GenericAccumulators::whole_sample_histogram`] returns.
///
/// **Every read group starts at [`DEFAULT_ERROR_RATE`]**, or rather at the rung of `ladder`
/// nearest it. Where the alternation starts does not change where it ends — that is what
/// the harness's deliberately-wrong start measures — but it does change how many iterations
/// it takes, and a start at the middle of the range chemistry actually occupies is the
/// cheapest honest guess.
///
/// # Errors
///
/// As [`fit_coupled`].
///
/// # Panics
///
/// As [`fit_coupled`].
pub(crate) fn fit_coupled_from_tables(
    sample: &str,
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    whole_sample: &BTreeMap<Ploidy, DepthAltHistogram<u64>>,
    ladder: &[ErrorRate],
) -> Result<CoupledFit, ParameterEstimationError> {
    let start = rung_nearest(ladder, DEFAULT_ERROR_RATE);
    let shares = library_shares(read_group_histograms);
    let start: BTreeMap<ReadGroupId, usize> = shares.keys().map(|&group| (group, start)).collect();

    alternate(
        sample,
        read_group_histograms,
        whole_sample,
        ladder,
        &start,
        MAX_COUPLED_FIT_ITERATIONS,
    )
    .map(|(fit, _)| fit)
}

/// The alternation, with the start and the cap named rather than taken from constants.
///
/// **Both are parameters for the same reason [`climb_with_cap`] takes its cap: so that what
/// they cost is a test rather than a recompile.** The start is what the harness's oracle
/// moves — three times the true rates — and the cap is what a test needs to see a fit that
/// ran out of iterations rather than settled.
///
/// Returns the fit and **every iterate's score in order**, which is how a test can assert
/// that the iterate kept is the best-scoring one rather than the last.
///
/// [`climb_with_cap`]: crate::ng::parameter_estimation::fitting::mixture_weights
fn alternate(
    sample: &str,
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    whole_sample: &BTreeMap<Ploidy, DepthAltHistogram<u64>>,
    ladder: &[ErrorRate],
    start: &BTreeMap<ReadGroupId, usize>,
    max_iterations: u32,
) -> Result<(CoupledFit, Vec<f64>), ParameterEstimationError> {
    assert!(!ladder.is_empty(), "a scan needs at least one rung to try");
    assert!(
        max_iterations > 0,
        "an alternation that runs no iterations has nothing to report"
    );

    let shares = library_shares(read_group_histograms);
    assert!(
        !shares.is_empty(),
        "sample {sample} has no read group with reads, so there is nothing to fit"
    );

    // The whole-sample cells, materialised once: the frequency climb walks them at every
    // iteration and so does the score, and both would otherwise re-derive the attributed
    // arm's keys twenty times over.
    let mut cells_of_ploidy: BTreeMap<Ploidy, Vec<Cell>> = BTreeMap::new();
    for (&ploidy, table) in whole_sample {
        let sites = table.total_loci();
        if sites < MIN_SITES_TO_FIT {
            return Err(ParameterEstimationError::GenotypeFrequenciesNotFittable {
                sample: sample.to_string(),
                ploidy,
                sites,
                floor: MIN_SITES_TO_FIT,
            });
        }
        cells_of_ploidy.insert(ploidy, table.cells(ploidy));
    }
    let all_cells: Vec<Cell> = cells_of_ploidy.values().flatten().cloned().collect();

    let mut rungs = start.clone();
    let mut scores: Vec<f64> = Vec::new();
    let mut best: Option<Iterate> = None;
    let mut iterations = 0;
    let mut converged = false;

    while iterations < max_iterations {
        iterations += 1;
        let noise = noise_from(&shares, &rungs, ladder);

        // Step 1 — the frequencies, from the whole-sample table at the current rates.
        let genotype_frequencies = climb_frequencies(&cells_of_ploidy, &noise);

        // Step 2 — each read group's rate, from its own table at those frequencies.
        let fitted =
            fit_read_group_error_rates(read_group_histograms, &genotype_frequencies, ladder);
        let next_rungs: BTreeMap<ReadGroupId, usize> = fitted
            .iter()
            .map(|(&group, fit)| (group, fit.rung))
            .collect();

        // The iterate is the pair (rates, frequencies) this round arrived at, and its score
        // is the whole-sample table's likelihood **at that pair** — one objective on one
        // table, which is what makes "best-scoring" a defined comparison between rounds.
        // Neither block's own score is: step 1's belongs to the previous rates and step 2's
        // to a different table.
        let noise = noise_from(&shares, &next_rungs, ladder);
        let score = whole_sample_score(&all_cells, &noise, &genotype_frequencies);
        scores.push(score.get());

        let settled = next_rungs == rungs;
        if best.as_ref().is_none_or(|kept| score.get() > kept.score) {
            best = Some(Iterate {
                rungs: next_rungs.clone(),
                genotype_frequencies,
                score: score.get(),
            });
        }
        rungs = next_rungs;
        if settled {
            converged = true;
            break;
        }
    }

    let best = best.expect("the loop runs at least once, so some iterate was scored");
    let termination = FitTermination {
        iterations,
        converged,
    };

    Ok((
        assemble(
            read_group_histograms,
            whole_sample,
            ladder,
            &best,
            termination,
        )?,
        scores,
    ))
}

/// One round's answer, and what it scored — kept so the loop can return the best rather
/// than the last.
struct Iterate {
    rungs: BTreeMap<ReadGroupId, usize>,
    genotype_frequencies: BTreeMap<Ploidy, SmallVec<[f64; 3]>>,
    /// The whole-sample table's weighted log-likelihood at this iterate's rates **and** its
    /// frequencies.
    score: f64,
}

/// Turn the winning iterate into the fit a caller reads, attaching each number's warrant.
fn assemble(
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    whole_sample: &BTreeMap<Ploidy, DepthAltHistogram<u64>>,
    ladder: &[ErrorRate],
    best: &Iterate,
    termination: FitTermination,
) -> Result<CoupledFit, ParameterEstimationError> {
    // **Reads and not sites**, because an error rate is per read
    // (`arch/parameter_prepass_generic.md` §2.4) and the two differ by the mean depth.
    let mut reads_of_group: BTreeMap<ReadGroupId, u64> = BTreeMap::new();
    for (&(group, _), table) in read_group_histograms {
        *reads_of_group.entry(group).or_default() += table.total_reads();
    }

    let error_rate = best
        .rungs
        .iter()
        .map(|(&group, &rung)| {
            (
                group,
                Estimate {
                    value: ladder[rung],
                    provenance: Provenance::FittedHere,
                    observations: reads_of_group.get(&group).copied().unwrap_or_default(),
                },
            )
        })
        .collect();

    let mut rates = BTreeMap::new();
    for (&ploidy, frequencies) in &best.genotype_frequencies {
        let by_alt_copies: SmallVec<[GenotypeFrequency; 5]> = frequencies
            .iter()
            .map(|&frequency| {
                GenotypeFrequency::try_new(frequency)
                    .expect("the climb returns a point on the simplex")
            })
            .collect();
        rates.insert(
            ploidy,
            Estimate {
                value: SampleRates::try_new(ploidy, by_alt_copies)?,
                provenance: Provenance::FittedHere,
                observations: whole_sample[&ploidy].total_loci(),
            },
        );
    }

    Ok(CoupledFit {
        error_rate,
        rates,
        termination,
    })
}

/// Each library's share of the sample's reads — the `w_g` of
/// `spec/parameter_prepass_generic.md` §1, summing to one over the sample.
///
/// **Reads and not sites, and not loci.** A library that contributed a tenth of the reads
/// contributed a tenth of the chances to misread one, whatever fraction of the sites it
/// happened to cover; the share weights a per-read probability.
///
/// A group whose table holds no reads is left out rather than given a share of zero, which
/// [`SampleLibraryNoise::new`] would refuse: a library with no reads has no rate to fit and
/// contributes no term to the mixture.
fn library_shares(
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
) -> BTreeMap<ReadGroupId, f64> {
    let mut reads_of_group: BTreeMap<ReadGroupId, u64> = BTreeMap::new();
    for (&(group, _), table) in read_group_histograms {
        *reads_of_group.entry(group).or_default() += table.total_reads();
    }
    reads_of_group.retain(|_, reads| *reads > 0);

    let total: u64 = reads_of_group.values().sum();
    reads_of_group
        .into_iter()
        .map(|(group, reads)| (group, reads as f64 / total as f64))
        .collect()
}

/// The noise parameters at one point of the alternation: every library's share paired with
/// the rate its rung names.
///
/// **The pairing is made here, in one place, from two maps keyed by read group** — never
/// from two collections indexed by position. A rule with two libraries' rates swapped
/// between their read groups is still a probability over the cell space, so none of the
/// scoring rule's identities can see it (spec §12.8); the only thing that prevents it is
/// that a share and a rate reach [`LibraryNoise`] under the same key.
///
/// # Panics
///
/// If a library has a share but no rung, or a rung but no share.
fn noise_from(
    shares: &BTreeMap<ReadGroupId, f64>,
    rungs: &BTreeMap<ReadGroupId, usize>,
    ladder: &[ErrorRate],
) -> SampleLibraryNoise {
    assert_eq!(
        shares.keys().collect::<Vec<_>>(),
        rungs.keys().collect::<Vec<_>>(),
        "the libraries with a share of the reads and the libraries with a fitted rate are \
         different sets"
    );
    SampleLibraryNoise::new(
        shares
            .iter()
            .map(|(&read_group, &share_of_reads)| LibraryNoise {
                read_group,
                share_of_reads,
                error_rate: ladder[rungs[&read_group]],
            }),
    )
}

/// Step 1 — the genotype frequencies that best explain the whole-sample table at these
/// noise parameters, one set per ploidy.
///
/// **Once per ploidy and not once for the table**, because a haploid cell has two genotypes
/// to mix and a diploid three, and a single vector would mean picking one of them and
/// dropping the rest.
fn climb_frequencies(
    cells_of_ploidy: &BTreeMap<Ploidy, Vec<Cell>>,
    noise: &SampleLibraryNoise,
) -> BTreeMap<Ploidy, SmallVec<[f64; 3]>> {
    let model = SubstitutionNoiseModel;
    let mut ln_likelihood_row_major: Vec<f64> = Vec::new();
    let mut cell_weights: Vec<f64> = Vec::new();
    let mut climbed = BTreeMap::new();

    for (&ploidy, cells) in cells_of_ploidy {
        let genotypes = model.genotypes(ploidy);
        ln_likelihood_row_major.clear();
        cell_weights.clear();
        for cell in cells {
            model.append_genotype_likelihoods(cell, noise, ploidy, &mut ln_likelihood_row_major);
            cell_weights.push(cell.sites as f64);
        }
        let table = GenotypeLikelihoodTable::from_natural_logs(&ln_likelihood_row_major, genotypes);
        climbed.insert(ploidy, fit_mixture_weights(table, &cell_weights));
    }
    climbed
}

/// The whole-sample table's weighted log-likelihood at one pair of rates and frequencies,
/// summed over the ploidies.
///
/// **A one-rung scan**, so that scoring a point and scoring a ladder go through exactly the
/// same checks — the width the model declares, the values it wrote, and the refusal of a
/// frequency set that explains none of the cells. The rail flag a one-rung ladder reports
/// is meaningless and is discarded.
fn whole_sample_score(
    cells: &[Cell],
    noise: &SampleLibraryNoise,
    genotype_frequencies: &BTreeMap<Ploidy, SmallVec<[f64; 3]>>,
) -> LogProb {
    fit_by_fixed_frequency_scan(
        &SubstitutionNoiseModel,
        cells,
        std::slice::from_ref(noise),
        genotype_frequencies,
    )
    .log_likelihood
}

/// The rung of `ladder` whose rate is closest to `rate`, in log space — the ladder is
/// geometric, so "closest" there and "closest in Phred" are the same question and
/// "closest in probability" is not.
///
/// # Panics
///
/// If `ladder` is empty.
fn rung_nearest(ladder: &[ErrorRate], rate: f64) -> usize {
    assert!(!ladder.is_empty(), "a scan needs at least one rung to try");
    let wanted = rate.ln();
    ladder
        .iter()
        .enumerate()
        .min_by(|(_, left), (_, right)| {
            let left = (left.get().ln() - wanted).abs();
            let right = (right.get().ln() - wanted).abs();
            left.total_cmp(&right)
        })
        .map(|(rung, _)| rung)
        .expect("the ladder is not empty")
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use crate::ng::parameter_estimation::fitting::ladder_scan::fit_by_profile_scan;
    use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
    use crate::ng::parameter_estimation::generic::error_rate_ladder;
    use crate::ng::parameter_estimation::generic::expected_counts::{
        alternative_read_probability, table_generated_at,
    };

    /// Phred 30 — rung 80 of the adopted ladder, which starts at Phred 10 and steps by
    /// 0.25. Also [`DEFAULT_ERROR_RATE`], so it is where the alternation starts by default.
    const RUNG_AT_PHRED_30: usize = 80;
    /// Phred 26, rung 64 — the second library, four Phred worse than the first.
    const RUNG_AT_PHRED_26: usize = 64;

    /// Tomato-like: heterozygous at 1.5 sites in a thousand, homozygous non-reference at
    /// 0.5.
    const TRUTH: [f64; 3] = [0.998, 0.0015, 0.0005];
    /// Each library reads every site this deep, so the sample's own depth is twice it.
    const PER_LIBRARY_DEPTH: u32 = 20;
    const SITES: f64 = 200_000.0;

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("a positive copy number")
    }

    /// A two-library world whose truth is known exactly: both libraries read every one of
    /// `SITES` sites [`PER_LIBRARY_DEPTH`] deep, at their own error rates, so their shares
    /// are a half each and the sample's own table sits at twice that depth.
    ///
    /// **The whole-sample table is generated at the share-weighted rate**, which is what
    /// the scoring rule reduces a pooled cell to (`noise_model`'s `share_weighted_rates`).
    /// That is not an approximation of the split — it is the exact sum of the attributed
    /// form over the splits a pooled key forgets — so the truth below really is a fixed
    /// point of both blocks rather than of one.
    struct TwoLibraryWorld {
        read_group_histograms: BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
        whole_sample: BTreeMap<Ploidy, DepthAltHistogram<u64>>,
        ladder: Vec<ErrorRate>,
    }

    impl TwoLibraryWorld {
        fn build() -> Self {
            let edges = Arc::new(DepthBinEdges::new());
            let ladder = error_rate_ladder();
            let diploid = ploidy(2);
            let rates = [
                ladder[RUNG_AT_PHRED_30].get(),
                ladder[RUNG_AT_PHRED_26].get(),
            ];

            let read_group_histograms = (1u32..=2)
                .map(|group| {
                    (
                        (ReadGroupId(group), diploid),
                        table_generated_at(
                            &edges,
                            PER_LIBRARY_DEPTH,
                            rates[group as usize - 1],
                            diploid,
                            &TRUTH,
                            SITES,
                        ),
                    )
                })
                .collect();

            // The sample's own table: one entry per site at the total depth, its
            // alternative reads drawn at the half-and-half mixture of the two rates. Built
            // by generating at each dosage's share-weighted `p_j` — which is what
            // `table_generated_at` does when handed a rate, so the rate handed in is the
            // one whose `p_j` equals the mixture at every dosage. Both libraries carry the
            // same ploidy, so that rate is the plain mean.
            let pooled_rate = 0.5 * rates[0] + 0.5 * rates[1];
            for (alt_copies, dosage) in [(0u8, 0), (1, 1), (2, 2)] {
                let mixture = 0.5
                    * (alternative_read_probability(alt_copies, diploid, rates[0])
                        + alternative_read_probability(alt_copies, diploid, rates[1]));
                let single = alternative_read_probability(alt_copies, diploid, pooled_rate);
                assert!(
                    (mixture - single).abs() < 1e-15,
                    "dosage {dosage}: the share-weighted rate is not the mean rate's"
                );
            }

            let whole_sample = BTreeMap::from([(
                diploid,
                table_generated_at(
                    &edges,
                    2 * PER_LIBRARY_DEPTH,
                    pooled_rate,
                    diploid,
                    &TRUTH,
                    SITES,
                ),
            )]);

            Self {
                read_group_histograms,
                whole_sample,
                ladder,
            }
        }

        /// The starting rungs the harness's oracle uses: **three times** each library's
        /// true rate, so that a fixed point at the truth is a result rather than a
        /// starting condition.
        fn three_times_the_truth(&self) -> BTreeMap<ReadGroupId, usize> {
            [RUNG_AT_PHRED_30, RUNG_AT_PHRED_26]
                .into_iter()
                .enumerate()
                .map(|(index, rung)| {
                    (
                        ReadGroupId(index as u32 + 1),
                        rung_nearest(&self.ladder, 3.0 * self.ladder[rung].get()),
                    )
                })
                .collect()
        }

        fn fit_from(
            &self,
            start: &BTreeMap<ReadGroupId, usize>,
            max_iterations: u32,
        ) -> (CoupledFit, Vec<f64>) {
            alternate(
                "world",
                &self.read_group_histograms,
                &self.whole_sample,
                &self.ladder,
                start,
                max_iterations,
            )
            .expect("the world holds enough sites to fit")
        }
    }

    fn rung_of(fit: &CoupledFit, ladder: &[ErrorRate], group: u32) -> usize {
        let rate = fit.error_rate[&ReadGroupId(group)].value;
        ladder
            .iter()
            .position(|&rung| rung == rate)
            .expect("the fitted rate is a rung of the ladder")
    }

    /// **E2's oracle: from a deliberately wrong start, the fixed point is the truth.**
    /// Three times each library's true rate going in, and both come back on their own rung
    /// — 80 and 64, sixteen apart — with the genotype frequencies at the truth.
    ///
    /// This is the harness's own experiment reduced to two libraries: from a start at three
    /// times the true rates and half the true frequencies, `ng_multilib_key_harness.rs`
    /// lands on the truth in all 25 worlds, to 0.000 rungs and 0.000% (research note §2.6).
    /// The start here moves only the rates, because the loop's first step derives the
    /// frequencies from them — so the frequency start is not a free parameter of this
    /// procedure, which is itself a difference from the harness worth knowing.
    #[test]
    fn from_three_times_the_true_rates_the_fixed_point_is_the_truth() {
        let world = TwoLibraryWorld::build();
        let start = world.three_times_the_truth();
        assert_ne!(
            start[&ReadGroupId(1)],
            RUNG_AT_PHRED_30,
            "the start has to be somewhere else, or there is nothing to converge from"
        );

        let (fit, _) = world.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);

        assert_eq!(rung_of(&fit, &world.ladder, 1), RUNG_AT_PHRED_30);
        assert_eq!(rung_of(&fit, &world.ladder, 2), RUNG_AT_PHRED_26);
        assert!(fit.termination.converged, "{:?}", fit.termination);

        let rates = &fit.rates[&ploidy(2)].value;
        for (dosage, (&fitted, &truth)) in rates
            .by_alt_copies()
            .iter()
            .map(|f| f.get())
            .collect::<Vec<_>>()
            .iter()
            .zip(&TRUTH)
            .enumerate()
        {
            assert!(
                (fitted - truth).abs() < 0.01 * truth,
                "dosage {dosage}: fitted {fitted}, truth {truth}"
            );
        }
    }

    /// **The rates and the frequencies really are coupled here**: a loop that fitted the
    /// rates once and stopped would keep its start. Same world, same procedure, started at
    /// a *different* wrong place — the ladder's own default rung for both libraries —
    /// and it lands on the same answer.
    ///
    /// Two starts landing together is what "fixed point" means, and it is what one start
    /// cannot say: a procedure that ignored its start entirely would also pass the test
    /// above, and a procedure that returned its start would pass neither.
    #[test]
    fn two_different_starts_reach_the_same_fixed_point() {
        let world = TwoLibraryWorld::build();
        let default_start: BTreeMap<ReadGroupId, usize> = (1u32..=2)
            .map(|group| {
                (
                    ReadGroupId(group),
                    rung_nearest(&world.ladder, DEFAULT_ERROR_RATE),
                )
            })
            .collect();

        let (from_default, _) = world.fit_from(&default_start, MAX_COUPLED_FIT_ITERATIONS);
        let (from_three_times, _) =
            world.fit_from(&world.three_times_the_truth(), MAX_COUPLED_FIT_ITERATIONS);

        // **The rates are compared exactly and the frequencies are not**, and the asymmetry
        // is the design's: a rate is a rung of the ladder, so two starts either land on the
        // same rung or they do not. The frequencies come off a climb whose pass count
        // depends on the rates it was handed at every round, so two starts agree to about
        // ten significant figures rather than to the bit — measured, 0.99800499 either way.
        // Convergence here is linear (research note §2.6), which is exactly the shape that
        // gives a great many correct digits and no exact equality.
        assert_eq!(from_default.error_rate, from_three_times.error_rate);
        for (dosage, (left, right)) in from_default.rates[&ploidy(2)]
            .value
            .by_alt_copies()
            .iter()
            .zip(from_three_times.rates[&ploidy(2)].value.by_alt_copies())
            .enumerate()
        {
            assert!(
                (left.get() - right.get()).abs() < 1e-6 * left.get(),
                "dosage {dosage}: {} from the default start, {} from three times the truth",
                left.get(),
                right.get()
            );
        }
    }

    /// **The amended second oracle: at one read group the alternation reaches the profile
    /// scan's answer.**
    ///
    /// The plan asked for "at one read group it terminates after one iteration", which is
    /// true of a profile scan and **false** of this alternation — with one library the
    /// alternation is plain coordinate ascent on a single objective, so it iterates. What
    /// it costs is iterations and not answers: each block is an exact maximisation of one
    /// objective, and at one read group the two tables are the same table
    /// (`spec/parameter_prepass_generic.md` §12.6), so both procedures converge to the same
    /// joint maximum. That is the stronger property, and it is what is asserted here.
    ///
    /// It is also the only consumer the profiling scan has on this path, which is the other
    /// reason to keep the comparison.
    #[test]
    fn at_one_read_group_the_alternation_reaches_the_profile_scans_answer() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let table = table_generated_at(
            &edges,
            PER_LIBRARY_DEPTH,
            ladder[RUNG_AT_PHRED_26].get(),
            diploid,
            &TRUTH,
            SITES,
        );
        // One library, so the two tables coincide — which is the premise, not a shortcut.
        let cells = table.cells(diploid);
        let read_group_histograms = BTreeMap::from([(
            (ReadGroupId(7), diploid),
            table_generated_at(
                &edges,
                PER_LIBRARY_DEPTH,
                ladder[RUNG_AT_PHRED_26].get(),
                diploid,
                &TRUTH,
                SITES,
            ),
        )]);
        let whole_sample = BTreeMap::from([(diploid, table)]);

        let alternated = fit_coupled_from_tables(
            "one-library",
            &read_group_histograms,
            &whole_sample,
            &ladder,
        )
        .expect("enough sites");

        let noise_ladder: Vec<SampleLibraryNoise> = ladder
            .iter()
            .map(|&rate| SampleLibraryNoise::single(ReadGroupId(7), rate))
            .collect();
        let profiled = fit_by_profile_scan(&SubstitutionNoiseModel, &cells, &noise_ladder);

        assert_eq!(
            alternated.error_rate[&ReadGroupId(7)].value,
            profiled.noise.libraries()[0].error_rate,
            "the alternation and the profile scan disagree about the rate"
        );
        for (dosage, (&climbed, fitted)) in profiled.genotype_frequencies[&diploid]
            .iter()
            .zip(alternated.rates[&diploid].value.by_alt_copies())
            .enumerate()
        {
            assert!(
                (climbed - fitted.get()).abs() < 1e-6,
                "dosage {dosage}: the scan climbed to {climbed}, the alternation to {}",
                fitted.get()
            );
        }
    }

    /// **The iterate kept is the best-scoring one, not the last** — which is the whole
    /// reason the loop keeps one at all, since a fixed point of two estimating equations
    /// has no monotonicity to lean on.
    ///
    /// Asserted by re-scoring the answer against every iterate's score: the fit's own score
    /// must be the largest of them.
    ///
    /// **What this test cannot say, stated rather than implied.** On this world the
    /// alternation settles in two iterations and the second scores higher than the first,
    /// so "keep the best" and "keep the last" agree here — the assertion below separates
    /// both from "keep an arbitrary iterate" and from a fit that reports numbers it never
    /// scored, and it does not separate them from each other. No fixture in hand makes a
    /// later iterate worse than an earlier one, which is the case the rule exists for; the
    /// real 583-cell tables of Milestone F2 are where such a trace could turn up.
    #[test]
    fn the_iterate_kept_is_the_best_scoring_one() {
        let world = TwoLibraryWorld::build();
        let (fit, scores) = world.fit_from(&world.three_times_the_truth(), 3);
        assert_eq!(
            scores.len(),
            2,
            "there has to be more than one iterate for 'best' to mean anything: {scores:?}"
        );

        let shares = library_shares(&world.read_group_histograms);
        let rungs: BTreeMap<ReadGroupId, usize> = fit
            .error_rate
            .keys()
            .map(|&group| (group, rung_of(&fit, &world.ladder, group.get())))
            .collect();
        let noise = noise_from(&shares, &rungs, &world.ladder);
        let frequencies: BTreeMap<Ploidy, SmallVec<[f64; 3]>> = fit
            .rates
            .iter()
            .map(|(&ploidy, estimate)| {
                (
                    ploidy,
                    estimate
                        .value
                        .by_alt_copies()
                        .iter()
                        .map(|f| f.get())
                        .collect(),
                )
            })
            .collect();
        let cells: Vec<Cell> = world
            .whole_sample
            .iter()
            .flat_map(|(&ploidy, table)| table.cells(ploidy))
            .collect();

        let kept = whole_sample_score(&cells, &noise, &frequencies).get();
        let best = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        assert!(
            (kept - best).abs() < 1e-6,
            "the fit scores {kept} where the best iterate scored {best}, of {scores:?}"
        );
    }

    /// A loop that ran out of iterations says so, and one that settled says so — the
    /// distinction a caller would otherwise consume as though it had settled.
    ///
    /// **Two iterations from a start nineteen rungs away, and that is the fewest a start
    /// away from the answer can take**: the first reaches the truth's rungs and the second
    /// is what observes that they did not move. Pinned as an equality rather than as
    /// "more than one", so that a loop silently running to its cap of twenty is a failure
    /// rather than a passing test with a slow fixture.
    #[test]
    fn a_capped_alternation_reports_that_it_did_not_settle() {
        let world = TwoLibraryWorld::build();
        let start = world.three_times_the_truth();

        let (capped, _) = world.fit_from(&start, 1);
        let (settled, _) = world.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);

        assert_eq!(capped.termination.iterations, 1);
        assert!(!capped.termination.converged);
        assert!(settled.termination.converged);
        assert_eq!(settled.termination.iterations, 2);
    }

    /// **A library's share is its share of the sample's *reads*, not of its sites.** Two
    /// libraries covering the same number of sites at four times the depth are not equal
    /// partners in a per-read error rate, and a share taken from site counts would say they
    /// were.
    #[test]
    fn a_librarys_share_is_read_weighted_and_not_site_weighted() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let rate = ladder[RUNG_AT_PHRED_30].get();
        let histograms = BTreeMap::from([
            (
                (ReadGroupId(1), diploid),
                table_generated_at(&edges, 40, rate, diploid, &TRUTH, 10_000.0),
            ),
            (
                (ReadGroupId(2), diploid),
                table_generated_at(&edges, 10, rate, diploid, &TRUTH, 10_000.0),
            ),
        ]);

        let shares = library_shares(&histograms);

        assert_eq!(shares.len(), 2);
        assert!(
            (shares[&ReadGroupId(1)] - 0.8).abs() < 1e-3,
            "the deep library's share is {}, where equal site counts would give 0.5",
            shares[&ReadGroupId(1)]
        );
        assert!((shares.values().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    /// A ploidy whose whole-sample table is too thin fails rather than emitting: a sample
    /// has one heterozygosity, so there is nothing to borrow it from.
    #[test]
    fn a_sample_too_thin_to_fit_its_frequencies_is_refused() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let rate = ladder[RUNG_AT_PHRED_30].get();
        let thin = 1_000.0;
        assert!(thin < MIN_SITES_TO_FIT as f64);

        let error = fit_coupled_from_tables(
            "thin",
            &BTreeMap::from([(
                (ReadGroupId(1), diploid),
                table_generated_at(&edges, PER_LIBRARY_DEPTH, rate, diploid, &TRUTH, thin),
            )]),
            &BTreeMap::from([(
                diploid,
                table_generated_at(&edges, PER_LIBRARY_DEPTH, rate, diploid, &TRUTH, thin),
            )]),
            &ladder,
        )
        .expect_err("a thousand sites is below the floor");

        assert!(
            matches!(
                error,
                ParameterEstimationError::GenotypeFrequenciesNotFittable { ref sample, .. }
                    if sample == "thin"
            ),
            "{error}"
        );
        assert!(error.to_string().contains("10000"), "{error}");
    }

    /// The default start is the ladder rung at [`DEFAULT_ERROR_RATE`], which on the adopted
    /// ladder is rung 80 exactly — Phred 30. Pinned because "nearest in log space" and
    /// "nearest in probability" are different rungs on a geometric ladder, and only the
    /// first is the one a Phred-spaced ladder means.
    #[test]
    fn the_default_start_is_the_rung_at_the_default_error_rate() {
        let ladder = error_rate_ladder();

        assert_eq!(rung_nearest(&ladder, DEFAULT_ERROR_RATE), RUNG_AT_PHRED_30);
        assert!((ladder[RUNG_AT_PHRED_30].get() - DEFAULT_ERROR_RATE).abs() < 1e-15);
        // A rate between two rungs goes to the nearer one in log space, which is the
        // higher-Phred rung when it is nearer there — the half-step lands on the ladder's
        // own geometric midpoint, so this asks for the rung just below.
        let between = (ladder[80].get() * ladder[81].get()).sqrt();
        assert!(matches!(rung_nearest(&ladder, between * 1.001), 80));
    }

    /// The pairing of a share to a rate goes through one key per library, so a set of
    /// shares and a set of rungs that name different libraries is a fault rather than a
    /// silent mismatch — the failure spec §12.8's identities cannot see, because a rule
    /// with two libraries' rates swapped is still a probability.
    #[test]
    #[should_panic(expected = "different sets")]
    fn shares_and_rungs_that_name_different_libraries_are_refused() {
        let ladder = error_rate_ladder();
        let shares = BTreeMap::from([(ReadGroupId(1), 0.5), (ReadGroupId(2), 0.5)]);
        let rungs = BTreeMap::from([(ReadGroupId(1), 80), (ReadGroupId(3), 64)]);

        let _ = noise_from(&shares, &rungs, &ladder);
    }

    /// Each library's rate reaches the library that produced it. Two libraries whose rungs
    /// are sixteen apart, and the noise set carries each rate under its own read group —
    /// the one thing `LibraryNoise` holding all three fields in a struct is for.
    #[test]
    fn each_librarys_rate_is_paired_with_its_own_read_group() {
        let ladder = error_rate_ladder();
        let shares = BTreeMap::from([(ReadGroupId(1), 0.25), (ReadGroupId(2), 0.75)]);
        let rungs = BTreeMap::from([
            (ReadGroupId(1), RUNG_AT_PHRED_30),
            (ReadGroupId(2), RUNG_AT_PHRED_26),
        ]);

        let noise = noise_from(&shares, &rungs, &ladder);

        let libraries = noise.libraries();
        assert_eq!(libraries.len(), 2);
        assert_eq!(libraries[0].read_group, ReadGroupId(1));
        assert_eq!(libraries[0].error_rate, ladder[RUNG_AT_PHRED_30]);
        assert!((libraries[0].share_of_reads - 0.25).abs() < 1e-15);
        assert_eq!(libraries[1].read_group, ReadGroupId(2));
        assert_eq!(libraries[1].error_rate, ladder[RUNG_AT_PHRED_26]);
        assert!((libraries[1].share_of_reads - 0.75).abs() < 1e-15);
    }

    /// An error rate's warrant is a count of **reads** and a genotype-frequency set's is a
    /// count of **sites** — the two differ by the mean depth, so an estimate carrying the
    /// wrong one under `observations` overstates or understates its own evidence by a
    /// factor of twenty here.
    #[test]
    fn each_estimate_carries_the_observation_count_its_quantity_is_per() {
        let world = TwoLibraryWorld::build();
        let (fit, _) = world.fit_from(&world.three_times_the_truth(), MAX_COUPLED_FIT_ITERATIONS);

        let diploid = ploidy(2);
        for group in 1u32..=2 {
            let estimate = &fit.error_rate[&ReadGroupId(group)];
            assert_eq!(
                estimate.observations,
                world.read_group_histograms[&(ReadGroupId(group), diploid)].total_reads()
            );
            assert_eq!(estimate.provenance, Provenance::FittedHere);
        }
        assert_eq!(
            fit.rates[&diploid].observations,
            world.whole_sample[&diploid].total_loci()
        );
        assert!(
            fit.error_rate[&ReadGroupId(1)].observations > 15 * fit.rates[&diploid].observations,
            "reads and sites differ by the depth, so a swap is not a small error"
        );
    }
}
