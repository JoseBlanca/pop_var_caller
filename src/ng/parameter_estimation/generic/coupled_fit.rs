//! The coupled fit: each library's error rate and the individual's genotype frequencies,
//! fitted together by alternating between the two tables that hold them.
//!
//! **Why they cannot be fitted apart.** A higher error rate explains the same alternative
//! reads as less real variation, so `ε` and the genotype frequencies trade off inside one
//! likelihood — and they are read off two different tables, the rates from the read-group
//! one and the frequencies from the whole-sample one
//! (`spec/parameter_prepass_generic.md` §5.1).
//!
//! **What one iteration is**, in the order the research harness runs it:
//!
//! 1. **The frequencies**, climbed on the whole-sample table at the rates the previous
//!    iteration produced — one set per ploidy present, because a haploid region has two
//!    genotype classes and a diploid three.
//! 2. **Each read group's rate**, scanned over the error-rate ladder on that group's own
//!    table, at the frequencies step 1 just produced and **without re-climbing them**.
//!    That is `read_group_error_rate::fit_read_group_error_rates`, and it is
//!    `examples/ng_multilib_key_harness.rs`'s `fit_eps_on_read_group(space, freqs)`.
//!
//! **Of those two, only the second is a property the tests defend.** Not re-climbing is
//! what makes this the estimator that was measured, and E1's
//! `the_frequencies_handed_in_move_the_fitted_rate` is what holds it. The *order* is
//! followed because the harness follows it, and nothing here can tell the two orders
//! apart: they share the same fixed point, and what differs is only which half of a round
//! an unconverged answer was reported from. Reversing it left every test green.
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
//! (`spec/parameter_prepass_generic.md` §1), so the alternation is plain coordinate ascent
//! on a single objective and reaches the same joint maximum the profile scan returns. That
//! is 1,550 of the 1,707 samples in the tomato archive survey; the coupling bites only on
//! the 157 multi-library ones, and on neither cohort in hand.
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
use crate::ng::parameter_estimation::generic::fallback::resolve_error_rates;
use crate::ng::parameter_estimation::generic::histogram::{Cell, DepthAltHistogram};
use crate::ng::parameter_estimation::generic::noise_model::{
    LibraryNoise, SampleLibraryNoise, SubstitutionNoiseModel,
};
use crate::ng::parameter_estimation::generic::read_group_error_rate::{
    ReadGroupErrorRateFit, fit_read_group_error_rates,
};
use crate::ng::parameter_estimation::generic::{
    CoupledFit, DEFAULT_ERROR_RATE, MAX_COUPLED_FIT_ITERATIONS, MIN_SITES_TO_FIT, SampleRates,
};
use crate::ng::parameter_estimation::{Estimate, ParameterEstimationError, Provenance};
use crate::ng::types::{ErrorRate, GenotypeFrequency, LogProb, Ploidy, ReadGroupId};

/// Fit a sample's error rates and genotype frequencies together, from its accumulators.
///
/// The thin door: it pulls the two tables out and hands them to
/// [`fit_coupled_from_tables`], which is where the alternation is and what most of the
/// tests drive, because a table built cell by cell is a fixture whose right answer can be
/// stated and a walked accumulator is not.
///
/// **This function has its own test even so**, and it needed one: the ploidy collection and
/// the fold below are the only things it does, and gutting them left all 3,147 tests green.
/// An accumulator is built in memory from hand-written loci, so reaching it costs no more
/// than reaching the tables.
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
    supplied: &BTreeMap<ReadGroupId, ErrorRate>,
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
        supplied,
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
    supplied: &BTreeMap<ReadGroupId, ErrorRate>,
) -> Result<CoupledFit, ParameterEstimationError> {
    let start = nearest_rung(ladder, DEFAULT_ERROR_RATE);
    let shares = library_shares(read_group_histograms);
    let start: BTreeMap<ReadGroupId, usize> = shares.keys().map(|&group| (group, start)).collect();

    fit_by_alternation(
        sample,
        read_group_histograms,
        whole_sample,
        ladder,
        supplied,
        &start,
        MAX_COUPLED_FIT_ITERATIONS,
    )
    .map(|(fit, _)| fit)
}

/// The alternation, with the start and the cap named rather than taken from constants.
///
/// **Both are parameters for the same reason `mixture_weights`' `climb_with_cap` takes its
/// cap: so that what they cost is a test rather than a recompile.** The start is what the
/// harness's oracle moves — three times the true rates — and the cap is what a test needs to
/// see a fit that ran out of iterations rather than settled.
///
/// Returns the fit and **every iterate in order**, so that a test can assert the fit is the
/// argmax of the trace rather than its last entry. The scores alone are not enough for
/// that: on a converged world every iterate scores the same to the bit, and the rungs and
/// frequencies are what separate them.
fn fit_by_alternation(
    sample: &str,
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    whole_sample: &BTreeMap<Ploidy, DepthAltHistogram<u64>>,
    ladder: &[ErrorRate],
    supplied: &BTreeMap<ReadGroupId, ErrorRate>,
    start: &BTreeMap<ReadGroupId, usize>,
    max_iterations: u32,
) -> Result<(CoupledFit, Vec<ScoredIterate>), ParameterEstimationError> {
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
    let mut trace: Vec<ScoredIterate> = Vec::new();
    // The last round's per-group fits, kept for their **site counts**: how much evidence
    // stood behind a group is what Milestone E4's fallback ladder gates on, and it does
    // not move between rounds.
    let mut fitted_sites: BTreeMap<ReadGroupId, ReadGroupErrorRateFit> = BTreeMap::new();
    let mut best: Option<ScoredIterate> = None;
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
        fitted_sites = fitted;

        // The iterate is the pair (rates, frequencies) this round arrived at, and its score
        // is the whole-sample table's likelihood **at that pair** — one objective on one
        // table, which is what makes "best-scoring" a defined comparison between rounds.
        // Neither block's own score is: step 1's belongs to the previous rates and step 2's
        // to a different table.
        let noise = noise_from(&shares, &next_rungs, ladder);
        let score = whole_sample_score(&all_cells, &noise, &genotype_frequencies);

        let settled = next_rungs == rungs;
        let scored = ScoredIterate {
            rungs: next_rungs.clone(),
            genotype_frequencies,
            score: score.get(),
        };
        // `>=` and not `>`: a tie keeps the **later** iterate, the same positional rule
        // the ladder scan uses for a tied rung. It is not cosmetic — on a converged world
        // every iterate scores the same `f64` to the bit, because they differ only at the
        // eleventh significant figure of the frequencies and a log-likelihood over 200,000
        // sites cannot resolve that. Keeping the first would then report the *less*
        // settled of two equally-scoring answers, and the rule exists to guard against a
        // worse last iterate, not against a more converged one.
        if best.as_ref().is_none_or(|kept| scored.score >= kept.score) {
            best = Some(scored.clone());
        }
        trace.push(scored);
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
        into_coupled_fit(
            read_group_histograms,
            whole_sample,
            ladder,
            supplied,
            &fitted_sites,
            &best,
            termination,
        )?,
        trace,
    ))
}

/// One round's answer, and what it scored — kept so the loop can return the best rather
/// than the last.
#[derive(Clone, PartialEq, Debug)]
struct ScoredIterate {
    rungs: BTreeMap<ReadGroupId, usize>,
    genotype_frequencies: BTreeMap<Ploidy, SmallVec<[f64; 3]>>,
    /// The whole-sample table's weighted log-likelihood at this iterate's rates **and** its
    /// frequencies.
    score: f64,
}

/// Turn the winning iterate into the fit a caller reads, attaching each number's warrant.
fn into_coupled_fit(
    read_group_histograms: &BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
    whole_sample: &BTreeMap<Ploidy, DepthAltHistogram<u64>>,
    ladder: &[ErrorRate],
    supplied: &BTreeMap<ReadGroupId, ErrorRate>,
    fitted_sites: &BTreeMap<ReadGroupId, ReadGroupErrorRateFit>,
    best: &ScoredIterate,
    termination: FitTermination,
) -> Result<CoupledFit, ParameterEstimationError> {
    // **Reads and not sites**, because an error rate is per read
    // (`arch/parameter_prepass_generic.md` §2.4) and the two differ by the mean depth.
    let mut reads_of_group: BTreeMap<ReadGroupId, u64> = BTreeMap::new();
    for (&(group, _), table) in read_group_histograms {
        *reads_of_group.entry(group).or_default() += table.total_reads();
    }

    // **The winning iterate's rungs, not the last round's fits**, which is why this rebuilds
    // the per-group answer rather than passing the loop's own `fitted` straight on: the
    // rate a group gets is the one the best-scoring iterate chose, and only its site count
    // comes from the fit.
    // Both maps come from `fit_read_group_error_rates` over the same histograms, so their
    // key sets are the same every round. Stated as a check rather than assumed: a group in
    // one and not the other would silently take a lower rung of the fallback ladder — a
    // borrowed or defaulted rate where it had a fitted one — with nothing to show for it.
    assert!(
        fitted_sites.keys().eq(best.rungs.keys()),
        "the winning iterate names read groups {:?} and the last round's fits name {:?}",
        best.rungs.keys().collect::<Vec<_>>(),
        fitted_sites.keys().collect::<Vec<_>>()
    );
    let at_the_winning_rungs: BTreeMap<ReadGroupId, ReadGroupErrorRateFit> = fitted_sites
        .iter()
        .filter_map(|(&group, fit)| {
            best.rungs.get(&group).map(|&rung| {
                (
                    group,
                    ReadGroupErrorRateFit {
                        error_rate: ladder[rung],
                        rung,
                        ..*fit
                    },
                )
            })
        })
        .collect();
    let error_rate = resolve_error_rates(
        &at_the_winning_rungs,
        &reads_of_group,
        supplied,
        MIN_SITES_TO_FIT,
    );

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
pub(super) fn library_shares(
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
/// scoring rule's identities can see it (spec §12, check 8); the only thing that prevents it is
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
fn nearest_rung(ladder: &[ErrorRate], rate: f64) -> usize {
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
    ///
    /// **What this world cannot see, and it is most of what the loop does.** At twenty
    /// reads a library a heterozygote shows about twenty alternative reads and a
    /// sequencing error nought or one, so the two classes never overlap: the frequency
    /// climb returns the same answer whatever error rate it is handed, and the read-group
    /// scan returns rungs 80 and 64 from **every** start on the ladder — including rung 0,
    /// a hundred times the true rate. Measured: a whole-sample table claiming thirty
    /// heterozygotes in a hundred still gives rungs 80 and 64 here. So the arrow from
    /// block 1 to block 2 carries no information on this world, and four mutations
    /// survived it, one of them the very estimator the design rejects.
    ///
    /// [`CoupledWorld`] is the fixture where the two blocks genuinely trade off, and every
    /// claim about the *coupling* is asserted there. This one is for the fixed point.
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
            // alternative reads drawn at the half-and-half mixture of the two rates.
            //
            // **`p_j(ε) = j/P + ε·(1 − 4j/3P)` is affine in `ε`**, so the mean of the two
            // libraries' per-read probabilities *is* the probability at the mean rate — for
            // any rates, shares, ploidy and dosage. The check below therefore cannot fail
            // for the reason a reader would guess; what it pins is that `pooled_rate` was
            // built with the shares this world actually has (a half each) and not some
            // other weighting.
            //
            // The consequence worth carrying: **the whole-sample table identifies only the
            // share-weighted mean rate and never the two libraries separately.** That is
            // why a transposition of the two libraries' rates has to be caught by the
            // direct test on `noise_from` rather than by any fit on this world.
            let pooled_rate = 0.5 * rates[0] + 0.5 * rates[1];
            for (alt_copies, dosage) in [(0u8, 0), (1, 1), (2, 2)] {
                let mixture = 0.5
                    * (alternative_read_probability(alt_copies, diploid, rates[0])
                        + alternative_read_probability(alt_copies, diploid, rates[1]));
                let single = alternative_read_probability(alt_copies, diploid, pooled_rate);
                assert!(
                    (mixture - single).abs() < 1e-15,
                    "dosage {dosage}: the pooled rate was not built from this world's shares"
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
                        nearest_rung(&self.ladder, 3.0 * self.ladder[rung].get()),
                    )
                })
                .collect()
        }

        fn fit_from(
            &self,
            start: &BTreeMap<ReadGroupId, usize>,
            max_iterations: u32,
        ) -> (CoupledFit, Vec<ScoredIterate>) {
            fit_by_alternation(
                "world",
                &self.read_group_histograms,
                &self.whole_sample,
                &self.ladder,
                &BTreeMap::new(),
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
                    nearest_rung(&world.ladder, DEFAULT_ERROR_RATE),
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
    /// (`spec/parameter_prepass_generic.md` §1), so both procedures converge to the same
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
            &BTreeMap::new(),
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
    /// **A comparison of scores cannot say this on this world, and the first version of
    /// this test did not know that.** Measured: the two iterates score the same `f64` to
    /// the bit — `−55460.001665379226` both — because they differ only at the eleventh
    /// significant figure of the frequencies, which a log-likelihood over 200,000 sites
    /// cannot resolve. A test asserting "the fit's score is the largest of them" therefore
    /// passed a loop rewritten to keep the *first* iterate.
    ///
    /// So the assertion is on the **iterate**, not on its score: the fit must carry the
    /// rungs and the frequencies of the trace's argmax, compared exactly. That separates
    /// keep-the-best from keep-the-first at the eleventh figure, where the scores are
    /// equal.
    ///
    /// **What it still cannot say**: with every iterate scoring alike, keep-the-best and
    /// keep-the-last are the same rule *on this world*. They are separated by
    /// `the_reported_rates_are_the_winning_iterates_and_not_the_last_rounds`, on the quiet
    /// `CoupledWorld`, whose trace really does end on a round scoring worse than an earlier
    /// one. An earlier version of this paragraph said no fixture in hand produced such a
    /// trace; one three functions down the file does.
    #[test]
    fn the_iterate_kept_is_the_argmax_of_the_trace() {
        let world = TwoLibraryWorld::build();
        let (fit, trace) = world.fit_from(&world.three_times_the_truth(), 3);
        assert_eq!(
            trace.len(),
            2,
            "there has to be more than one iterate for an argmax to mean anything"
        );

        let argmax = trace
            .iter()
            .max_by(|left, right| left.score.total_cmp(&right.score))
            .expect("a non-empty trace");
        let reported_rungs: BTreeMap<ReadGroupId, usize> = fit
            .error_rate
            .keys()
            .map(|&group| (group, rung_of(&fit, &world.ladder, group.get())))
            .collect();
        let reported_frequencies: BTreeMap<Ploidy, SmallVec<[f64; 3]>> = fit
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

        assert_eq!(reported_rungs, argmax.rungs);
        assert_eq!(
            reported_frequencies, argmax.genotype_frequencies,
            "the fit reports frequencies from an iterate that is not the argmax"
        );
        // The premise of the whole test, stated so that a future fixture whose scores do
        // separate does not quietly change what is being proven here.
        assert_eq!(
            trace[0].score, trace[1].score,
            "on this world every iterate scores alike, which is why the assertion above is \
             on the iterate and not on its score"
        );
        assert_ne!(
            trace[0].genotype_frequencies, trace[1].genotype_frequencies,
            "and they differ in the frequencies, which is what makes it discriminating"
        );
    }

    /// **An iterate's score is taken at its own rates**, not at the rates it started the
    /// round with — the distinction the module doc's "one objective on one table" argument
    /// rests on.
    ///
    /// **Asserted on a capped run, and it has to be.** At a converged iterate the round's
    /// starting rates and its finishing rates are equal by definition, so the two candidate
    /// scores coincide and no converged fixture can tell them apart. Capped at one
    /// iteration the gap is 12,552 nats.
    #[test]
    fn an_iterates_score_is_taken_at_its_own_rates() {
        let world = TwoLibraryWorld::build();
        let start = world.three_times_the_truth();
        let (capped, trace) = world.fit_from(&start, 1);
        assert_eq!(trace.len(), 1);
        assert!(!capped.termination.converged);

        let shares = library_shares(&world.read_group_histograms);
        let rungs: BTreeMap<ReadGroupId, usize> = capped
            .error_rate
            .keys()
            .map(|&group| (group, rung_of(&capped, &world.ladder, group.get())))
            .collect();
        assert_ne!(
            rungs, start,
            "the round has to have moved the rates, or there is nothing to tell apart"
        );
        let frequencies: BTreeMap<Ploidy, SmallVec<[f64; 3]>> = capped
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

        let at_its_own = whole_sample_score(
            &cells,
            &noise_from(&shares, &rungs, &world.ladder),
            &frequencies,
        )
        .get();
        let at_the_start = whole_sample_score(
            &cells,
            &noise_from(&shares, &start, &world.ladder),
            &frequencies,
        )
        .get();
        assert!(
            (at_the_start - at_its_own).abs() > 1.0,
            "the two candidate scores have to differ for this to discriminate: \
             {at_its_own} against {at_the_start}"
        );
        assert!(
            (trace[0].score - at_its_own).abs() < 1e-6,
            "the reported score {} is not the score at the iterate's own rates \
             ({at_its_own}); the start's rates would give {at_the_start}",
            trace[0].score
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
            &BTreeMap::new(),
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

        assert_eq!(nearest_rung(&ladder, DEFAULT_ERROR_RATE), RUNG_AT_PHRED_30);
        assert!((ladder[RUNG_AT_PHRED_30].get() - DEFAULT_ERROR_RATE).abs() < 1e-15);

        // Which space the nearness is measured in is
        // `nearest_rung_breaks_the_tie_in_log_space_and_not_in_probability`'s, below —
        // `DEFAULT_ERROR_RATE` sits exactly on rung 80, at distance zero in either space,
        // so this test cannot say it.
    }

    /// The pairing of a share to a rate goes through one key per library, so a set of
    /// shares and a set of rungs that name different libraries is a fault rather than a
    /// silent mismatch — the failure spec §12's eighth check cannot see, because a rule
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
        // **The line with independent force**, and the one that catches a `total_reads`
        // summing site counts: the two assertions above compare each estimate against the
        // same function that produced it, so a wrong `total_reads` moves both sides
        // together. Each library reads every site twenty deep here, so an error rate's
        // warrant is twenty times its frequencies' — pooling the two libraries would give
        // forty, which is neither estimate's number.
        assert!(
            fit.error_rate[&ReadGroupId(1)].observations > 15 * fit.rates[&diploid].observations,
            "reads and sites differ by the depth, so a swap is not a small error"
        );
    }

    /// **The fallback ladder reaches the fit a caller reads.** A second library with 500
    /// sites — a twentieth of [`MIN_SITES_TO_FIT`] — is fitted like any other inside the
    /// loop, and comes out of it marked `Borrowed`, carrying the deep library's rate and
    /// the deep library's reads as its warrant.
    ///
    /// Without this the ladder is a function nothing calls: `CoupledFit` would report a
    /// rate fitted from 500 sites as `FittedHere`, which is the provenance record saying
    /// the opposite of what happened.
    #[test]
    fn a_thin_library_comes_out_of_the_coupled_fit_marked_borrowed() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let thin = 500.0;
        assert!(thin < MIN_SITES_TO_FIT as f64);

        let read_group_histograms = BTreeMap::from([
            (
                (ReadGroupId(1), diploid),
                table_generated_at(
                    &edges,
                    PER_LIBRARY_DEPTH,
                    ladder[RUNG_AT_PHRED_30].get(),
                    diploid,
                    &TRUTH,
                    SITES,
                ),
            ),
            (
                (ReadGroupId(2), diploid),
                table_generated_at(
                    &edges,
                    PER_LIBRARY_DEPTH,
                    ladder[RUNG_AT_PHRED_26].get(),
                    diploid,
                    &TRUTH,
                    thin,
                ),
            ),
        ]);
        let whole_sample = BTreeMap::from([(
            diploid,
            table_generated_at(
                &edges,
                PER_LIBRARY_DEPTH,
                ladder[RUNG_AT_PHRED_30].get(),
                diploid,
                &TRUTH,
                SITES,
            ),
        )]);

        let fit = fit_coupled_from_tables(
            "one-thin",
            &read_group_histograms,
            &whole_sample,
            &ladder,
            &BTreeMap::new(),
        )
        .expect("the whole-sample table holds enough sites");

        assert_eq!(
            fit.error_rate[&ReadGroupId(1)].provenance,
            Provenance::FittedHere
        );
        let thin_group = &fit.error_rate[&ReadGroupId(2)];
        assert_eq!(thin_group.provenance, Provenance::Borrowed);
        assert_eq!(
            thin_group.value,
            fit.error_rate[&ReadGroupId(1)].value,
            "with one lender the borrowed rate is that lender's"
        );
        assert_eq!(
            thin_group.observations,
            read_group_histograms[&(ReadGroupId(1), diploid)].total_reads(),
            "the warrant is the lender's reads, not the 500 sites that were too few"
        );
    }

    /// **What the run supplied reaches the fit a caller reads.** With *both* libraries below
    /// [`MIN_SITES_TO_FIT`] there is no lender, so the ladder falls past `Borrowed` to its
    /// two lower rungs: the group the run supplied a rate for comes out `Supplied` carrying
    /// that rate, and the group it said nothing about comes out `Defaulted`.
    ///
    /// This is the only test that hands either public door a non-empty supplied map, so it
    /// is what stops the parameter being threaded in and then dropped — a fit that ignored
    /// it would report `Defaulted` for both and no other test would notice.
    ///
    /// The whole-sample table still holds [`SITES`], because the frequencies' floor is a
    /// different floor from the per-group one this test is about.
    #[test]
    fn a_supplied_rate_reaches_the_coupled_fit_when_no_group_can_lend() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let thin = 500.0;
        assert!(thin < MIN_SITES_TO_FIT as f64);

        let read_group_histograms: BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>> = (1u32
            ..=2)
            .map(|group| {
                (
                    (ReadGroupId(group), diploid),
                    table_generated_at(
                        &edges,
                        PER_LIBRARY_DEPTH,
                        ladder[RUNG_AT_PHRED_30].get(),
                        diploid,
                        &TRUTH,
                        thin,
                    ),
                )
            })
            .collect();
        let whole_sample = BTreeMap::from([(
            diploid,
            table_generated_at(
                &edges,
                PER_LIBRARY_DEPTH,
                ladder[RUNG_AT_PHRED_30].get(),
                diploid,
                &TRUTH,
                SITES,
            ),
        )]);
        let supplied = BTreeMap::from([(ReadGroupId(1), ladder[RUNG_AT_PHRED_26])]);

        let fit = fit_coupled_from_tables(
            "supplied",
            &read_group_histograms,
            &whole_sample,
            &ladder,
            &supplied,
        )
        .expect("the whole-sample table holds enough sites");

        let told = &fit.error_rate[&ReadGroupId(1)];
        assert_eq!(told.provenance, Provenance::Supplied);
        assert_eq!(
            told.value, ladder[RUNG_AT_PHRED_26],
            "the rate the run handed in, not the one the thin table fitted"
        );
        assert_eq!(
            told.observations, 0,
            "a value nothing in this sample stood behind carries no observations"
        );
        assert_eq!(
            fit.error_rate[&ReadGroupId(2)].provenance,
            Provenance::Defaulted,
            "the group the run said nothing about falls to the last rung"
        );
    }

    /// **A library that contributed no reads is left out rather than given a share of
    /// zero**, which [`SampleLibraryNoise::new`] refuses — a library with no reads has no
    /// rate to fit and contributes no term to the mixture.
    #[test]
    fn a_library_with_no_reads_takes_no_share() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let histograms = BTreeMap::from([
            (
                (ReadGroupId(1), diploid),
                table_generated_at(
                    &edges,
                    PER_LIBRARY_DEPTH,
                    ladder[RUNG_AT_PHRED_30].get(),
                    diploid,
                    &TRUTH,
                    SITES,
                ),
            ),
            (
                (ReadGroupId(2), diploid),
                DepthAltHistogram::new(Arc::clone(&edges)),
            ),
        ]);

        let shares = library_shares(&histograms);

        assert_eq!(shares.keys().copied().collect::<Vec<_>>(), [ReadGroupId(1)]);
        assert!((shares[&ReadGroupId(1)] - 1.0).abs() < 1e-15);
    }

    /// **The public door agrees with the one the tests drive.** `fit_coupled` pulls the
    /// ploidies out of the windowed table and folds each one's windows into a whole-sample
    /// table; nothing else in the crate calls it, so without this its wiring could be wrong
    /// in any way at all and no test would move.
    ///
    /// The accumulator is built in memory from hand-written loci, which is what makes this
    /// cost the same as everything else here — the fit does not need a locus stream to be
    /// reached, only an accumulator.
    #[test]
    fn the_public_door_agrees_with_the_tables_it_folds() {
        use crate::ng::locus_generation::{
            LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
        };
        use crate::ng::parameter_estimation::generic::accumulators::{
            ConstantPloidy, GenericAccumulators, InbreedingMode,
        };
        use crate::ng::types::{ContigId, GenomeRegion, Position};

        let edges = Arc::new(DepthBinEdges::new());
        let diploid = ploidy(2);
        let mut accumulators = GenericAccumulators::new(
            Arc::clone(&edges),
            &[ReadGroupId(1)],
            Arc::new(ConstantPloidy(diploid)),
            InbreedingMode::Fitted,
        );

        let observation = |bases: &[u8], reads: u32| SequenceObservation {
            bases: bases.into(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(1),
            num_obs: reads,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        };

        // Enough sites to clear `MIN_SITES_TO_FIT`, with one alternative read at every
        // thousandth of them so the fit has something to explain.
        for site in 0..MIN_SITES_TO_FIT + 1_000 {
            let alt = site % 1_000 == 0;
            let observations = if alt {
                vec![
                    observation(b"C", PER_LIBRARY_DEPTH - 1),
                    observation(b"A", 1),
                ]
            } else {
                vec![observation(b"C", PER_LIBRARY_DEPTH)]
            };
            accumulators.add_locus(&SampleLocusObservations {
                kind: LocusKind::Generic,
                region: GenomeRegion {
                    contig: ContigId(0),
                    start: Position(site + 1),
                    end: Position(site + 1),
                },
                reference_bases: b"C".as_slice().into(),
                observations,
                reads_without_observation: 0,
                reads_discarded_by_cap: 0,
            });
        }

        let ladder = error_rate_ladder();
        let through_the_door = fit_coupled("walked", &accumulators, &ladder, &BTreeMap::new())
            .expect("the accumulator holds enough sites");
        let from_the_tables = fit_coupled_from_tables(
            "walked",
            accumulators.read_group_histograms(),
            &BTreeMap::from([(diploid, accumulators.whole_sample_histogram(diploid))]),
            &ladder,
            &BTreeMap::new(),
        )
        .expect("the same tables");

        assert_eq!(through_the_door.error_rate, from_the_tables.error_rate);
        assert_eq!(through_the_door.rates, from_the_tables.rates);
        assert_eq!(through_the_door.termination, from_the_tables.termination);
    }

    // ------------------------------------------------------------------
    // Review additions (tmp/review_e2_a) — the tests that kill the
    // mutations `TwoLibraryWorld` cannot see.
    // ------------------------------------------------------------------

    /// Phred 15, rung 20 — noisy enough that an error and a heterozygote are the same
    /// observation at the depths below, which is the regime the alternation exists for.
    const RUNG_AT_PHRED_15: usize = 20;
    /// Phred 17, rung 28 — the second library of the coupled world.
    const RUNG_AT_PHRED_17: usize = 28;
    /// A variable individual: 7.5 sites in a hundred heterozygous, 2.5 homozygous
    /// non-reference.
    const COUPLED_TRUTH: [f64; 3] = [0.90, 0.075, 0.025];
    /// Three reads per library per site, so the sample's own table sits at six.
    const COUPLED_DEPTH: u32 = 3;
    /// What the coupled world costs the loop, measured: four rounds, against the two
    /// `TwoLibraryWorld` takes from any start at all.
    const COUPLED_WORLD_ROUNDS: u32 = 4;

    /// **A world where the two blocks are genuinely coupled**, which
    /// [`TwoLibraryWorld`] is not.
    ///
    /// At [`TwoLibraryWorld`]'s depth — twenty reads per library, forty in the sample's
    /// own table — a heterozygote shows about twenty alternative reads and a sequencing
    /// error shows nought or one, so the two classes never overlap and the frequency
    /// climb returns the same answer whatever error rate it is handed. Measured on that
    /// world: the frequencies climbed at three times the true rates and at the true rates
    /// agree to ten significant figures (0.998004990024 against 0.998004990011), and the
    /// read-group scan returns rungs 80 and 64 from **every** start on the ladder,
    /// including rung 0 — a hundred times the true rate. Nothing there can tell a coupled
    /// loop from two independent fits.
    ///
    /// Here each library reads three deep at Phred 15 and 17 and the individual is
    /// variable, so an error and a heterozygote are the same observation. Climbing the
    /// frequencies at three times the true rate moves them by 0.075 — the whole
    /// heterozygous rate — and moves the rung the read-group scan then picks by ten,
    /// which is 2.5 Phred.
    struct CoupledWorld {
        read_group_histograms: BTreeMap<(ReadGroupId, Ploidy), DepthAltHistogram<u64>>,
        whole_sample: BTreeMap<Ploidy, DepthAltHistogram<u64>>,
        ladder: Vec<ErrorRate>,
    }

    impl CoupledWorld {
        fn build_with(truth_in_the_whole_sample_table: [f64; 3]) -> Self {
            let edges = Arc::new(DepthBinEdges::new());
            let ladder = error_rate_ladder();
            let diploid = ploidy(2);
            let rates = [
                ladder[RUNG_AT_PHRED_15].get(),
                ladder[RUNG_AT_PHRED_17].get(),
            ];

            let read_group_histograms = (1u32..=2)
                .map(|group| {
                    (
                        (ReadGroupId(group), diploid),
                        table_generated_at(
                            &edges,
                            COUPLED_DEPTH,
                            rates[group as usize - 1],
                            diploid,
                            &COUPLED_TRUTH,
                            SITES,
                        ),
                    )
                })
                .collect();

            let pooled_rate = 0.5 * rates[0] + 0.5 * rates[1];
            let whole_sample = BTreeMap::from([(
                diploid,
                table_generated_at(
                    &edges,
                    2 * COUPLED_DEPTH,
                    pooled_rate,
                    diploid,
                    &truth_in_the_whole_sample_table,
                    SITES,
                ),
            )]);

            Self {
                read_group_histograms,
                whole_sample,
                ladder,
            }
        }

        fn build() -> Self {
            Self::build_with(COUPLED_TRUTH)
        }

        fn fit_from(
            &self,
            start: &BTreeMap<ReadGroupId, usize>,
            max_iterations: u32,
        ) -> (CoupledFit, Vec<ScoredIterate>) {
            fit_by_alternation(
                "coupled",
                &self.read_group_histograms,
                &self.whole_sample,
                &self.ladder,
                &BTreeMap::new(),
                start,
                max_iterations,
            )
            .expect("the world holds enough sites to fit")
        }
    }

    fn start_at(first: usize, second: usize) -> BTreeMap<ReadGroupId, usize> {
        BTreeMap::from([(ReadGroupId(1), first), (ReadGroupId(2), second)])
    }

    /// **The frequencies the whole-sample table produced are what the read-group scan is
    /// handed** — the wire that makes this a coupled fit rather than two independent ones.
    ///
    /// Two samples with **identical** read-group tables and whole-sample tables that
    /// disagree about the individual: one says 7.5 sites in a hundred are heterozygous,
    /// the other 1.5 in a hundred. If step 2 saw those frequencies, the rates it returns
    /// must differ; if step 2 re-climbed its own frequencies from the read-group table —
    /// the estimator arch §5.2 describes and the harness rejected — the two would be
    /// identical, because nothing else about the two samples differs.
    ///
    /// This is the one assertion in the file that fails when the two blocks are
    /// disconnected. `TwoLibraryWorld` cannot make it: there, a whole-sample table
    /// claiming 30% heterozygotes still returns rungs 80 and 64.
    #[test]
    fn the_whole_sample_tables_frequencies_reach_the_fitted_rates() {
        let variable = CoupledWorld::build_with(COUPLED_TRUTH);
        let quiet = CoupledWorld::build_with([0.98, 0.015, 0.005]);
        for (&(group, at), table) in &variable.read_group_histograms {
            assert_eq!(
                table.cells(at),
                quiet.read_group_histograms[&(group, at)].cells(at),
                "the two samples must differ only in their whole-sample table"
            );
        }

        let start = start_at(RUNG_AT_PHRED_15, RUNG_AT_PHRED_17);
        let (from_variable, _) = variable.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);
        let (from_quiet, _) = quiet.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);

        assert_ne!(
            rung_of(&from_variable, &variable.ladder, 1),
            rung_of(&from_quiet, &quiet.ladder, 1),
            "read group 1 got the same rate from two different whole-sample tables, so \
             the frequency block is not feeding the rate block"
        );
        assert_ne!(
            rung_of(&from_variable, &variable.ladder, 2),
            rung_of(&from_quiet, &quiet.ladder, 2),
            "read group 2 got the same rate from two different whole-sample tables"
        );
    }

    /// **The rates a caller reads are the winning iterate's, not the last round's** — and
    /// unlike `the_iterate_kept_is_the_argmax_of_the_trace`, this fixture can tell the two
    /// apart. That test says in its own doc comment that no fixture in hand produces a trace
    /// whose last iterate is not the argmax; this one does. The coupled world whose
    /// whole-sample table claims 1.5 heterozygotes in a hundred, started at the *other*
    /// world's answer, ends on a round that scores worse than an earlier one, so
    /// keep-the-best and keep-the-last name different rungs.
    ///
    /// It is the assertion `into_coupled_fit`'s rebuild exists for: the fits the loop hands
    /// on belong to the last round, and only their site counts may be used — the rate has to
    /// come from `best`.
    #[test]
    fn the_reported_rates_are_the_winning_iterates_and_not_the_last_rounds() {
        let world = CoupledWorld::build_with([0.98, 0.015, 0.005]);
        let start = start_at(RUNG_AT_PHRED_15, RUNG_AT_PHRED_17);

        let (fit, trace) = world.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);

        // `max_by` keeps the later of two equal maxima, which is the loop's own `>=` rule.
        let argmax = trace
            .iter()
            .max_by(|left, right| left.score.total_cmp(&right.score))
            .expect("a non-empty trace");
        let last = trace.last().expect("a non-empty trace");
        assert_ne!(
            argmax.rungs, last.rungs,
            "the premise of this test: on this world the last round's rungs differ from the \
             best-scoring round's, which is what lets it separate the two rules"
        );

        let reported: BTreeMap<ReadGroupId, usize> = fit
            .error_rate
            .keys()
            .map(|&group| (group, rung_of(&fit, &world.ladder, group.get())))
            .collect();

        assert_eq!(
            reported, argmax.rungs,
            "the fit reports the last round's rates rather than the best-scoring round's"
        );
    }

    /// **On the coupled world the loop lands on the truth too, and it needs more than the
    /// two rounds `TwoLibraryWorld` takes.** Pinned as an equality for the same reason the
    /// two-round count is: a loop that reaches its cap must fail rather than pass slowly.
    #[test]
    fn the_coupled_world_reaches_its_truth_and_takes_more_than_two_rounds() {
        let world = CoupledWorld::build();
        let start = start_at(
            nearest_rung(&world.ladder, 3.0 * world.ladder[RUNG_AT_PHRED_15].get()),
            nearest_rung(&world.ladder, 3.0 * world.ladder[RUNG_AT_PHRED_17].get()),
        );

        let (fit, _) = world.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);

        assert_eq!(rung_of(&fit, &world.ladder, 1), RUNG_AT_PHRED_15);
        assert_eq!(rung_of(&fit, &world.ladder, 2), RUNG_AT_PHRED_17);
        assert!(fit.termination.converged, "{:?}", fit.termination);
        assert_eq!(fit.termination.iterations, COUPLED_WORLD_ROUNDS);
    }

    /// **A fit started at its own answer settles in one round.** The stopping rule is that
    /// no read group's rung moved, so the round that observes it is the first one — a rule
    /// that counted rounds instead would report two.
    #[test]
    fn a_fit_started_at_its_answer_settles_in_one_round() {
        let world = TwoLibraryWorld::build();

        let (fit, trace) = world.fit_from(
            &start_at(RUNG_AT_PHRED_30, RUNG_AT_PHRED_26),
            MAX_COUPLED_FIT_ITERATIONS,
        );

        assert_eq!(fit.termination.iterations, 1, "{:?}", fit.termination);
        assert!(fit.termination.converged);
        assert_eq!(trace.len(), 1);
        assert_eq!(rung_of(&fit, &world.ladder, 1), RUNG_AT_PHRED_30);
        assert_eq!(rung_of(&fit, &world.ladder, 2), RUNG_AT_PHRED_26);
    }

    /// **One library already at its answer does not stop the loop for the other.** The
    /// rule is *every* read group's rung, and a rule that watched one of them — the first,
    /// the last, any single one — would stop a round early here and report a converged fit
    /// whose second library had just moved nineteen rungs.
    #[test]
    fn one_library_already_home_does_not_settle_the_loop() {
        let world = TwoLibraryWorld::build();
        let start = start_at(
            RUNG_AT_PHRED_30,
            nearest_rung(&world.ladder, 3.0 * world.ladder[RUNG_AT_PHRED_26].get()),
        );
        assert_eq!(start[&ReadGroupId(1)], RUNG_AT_PHRED_30);
        assert_ne!(start[&ReadGroupId(2)], RUNG_AT_PHRED_26);

        let (fit, _) = world.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);

        assert_eq!(
            fit.termination.iterations, 2,
            "read group 1 started on its answer and read group 2 did not, so the first \
             round cannot be the one that observes nothing moved: {:?}",
            fit.termination
        );
        assert_eq!(rung_of(&fit, &world.ladder, 2), RUNG_AT_PHRED_26);
    }

    /// **[`fit_coupled_from_tables`] really does start every group at
    /// [`DEFAULT_ERROR_RATE`]'s rung**, which no test reached before: every start on the
    /// ladder reaches the same answer on these worlds, so only the round count can see
    /// where the loop began. A one-library sample whose true rate *is* the default rung
    /// settles in one round; from anywhere else it takes two.
    #[test]
    fn fit_coupled_from_tables_starts_at_the_default_error_rates_rung() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let diploid = ploidy(2);
        let at_the_default = |rung: usize| {
            let generate = || {
                table_generated_at(
                    &edges,
                    PER_LIBRARY_DEPTH,
                    ladder[rung].get(),
                    diploid,
                    &TRUTH,
                    SITES,
                )
            };
            fit_coupled_from_tables(
                "one-library",
                &BTreeMap::from([((ReadGroupId(1), diploid), generate())]),
                &BTreeMap::from([(diploid, generate())]),
                &ladder,
                &BTreeMap::new(),
            )
            .expect("enough sites")
        };

        // The default rung is the answer, so the first round observes that nothing moved.
        assert_eq!(at_the_default(RUNG_AT_PHRED_30).termination.iterations, 1);
        // Four Phred away, so it takes a round to get there and a round to see it stayed.
        assert_eq!(at_the_default(RUNG_AT_PHRED_26).termination.iterations, 2);
    }

    /// **A sample with two ploidies fits a frequency set for each**, which no test reached
    /// before: every fixture in this file is diploid throughout, so
    /// [`climb_frequencies`]'s per-ploidy loop ran exactly once and neither its genotype
    /// count nor its buffer reuse was ever exercised twice.
    ///
    /// A haploid region has two genotype classes and a diploid three, so a climb that used
    /// one ploidy's count for both, or carried one ploidy's likelihood rows into the
    /// next's, would be caught here and nowhere else.
    #[test]
    fn a_sample_with_a_haploid_and_a_diploid_region_fits_both() {
        let edges = Arc::new(DepthBinEdges::new());
        let ladder = error_rate_ladder();
        let haploid = ploidy(1);
        let diploid = ploidy(2);
        let rate = ladder[RUNG_AT_PHRED_30].get();
        let haploid_truth = [0.997, 0.003];

        let read_group_histograms = BTreeMap::from([
            (
                (ReadGroupId(1), haploid),
                table_generated_at(
                    &edges,
                    PER_LIBRARY_DEPTH,
                    rate,
                    haploid,
                    &haploid_truth,
                    SITES,
                ),
            ),
            (
                (ReadGroupId(1), diploid),
                table_generated_at(&edges, PER_LIBRARY_DEPTH, rate, diploid, &TRUTH, SITES),
            ),
        ]);
        let whole_sample = BTreeMap::from([
            (
                haploid,
                table_generated_at(
                    &edges,
                    PER_LIBRARY_DEPTH,
                    rate,
                    haploid,
                    &haploid_truth,
                    SITES,
                ),
            ),
            (
                diploid,
                table_generated_at(&edges, PER_LIBRARY_DEPTH, rate, diploid, &TRUTH, SITES),
            ),
        ]);

        let fit = fit_coupled_from_tables(
            "two-ploidies",
            &read_group_histograms,
            &whole_sample,
            &ladder,
            &BTreeMap::new(),
        )
        .expect("enough sites at both ploidies");

        // One rate for the read group, spanning both ploidies.
        assert_eq!(fit.error_rate.len(), 1);
        assert_eq!(
            fit.error_rate[&ReadGroupId(1)].value,
            ladder[RUNG_AT_PHRED_30]
        );

        // Two frequency sets, each as wide as its own ploidy's genotype set.
        assert_eq!(fit.rates.len(), 2);
        let fitted_haploid = &fit.rates[&haploid].value;
        let fitted_diploid = &fit.rates[&diploid].value;
        assert_eq!(fitted_haploid.by_alt_copies().len(), 2);
        assert_eq!(fitted_diploid.by_alt_copies().len(), 3);
        for (dosage, (fitted, truth)) in fitted_haploid
            .by_alt_copies()
            .iter()
            .zip(&haploid_truth)
            .enumerate()
        {
            assert!(
                (fitted.get() - truth).abs() < 0.01 * truth,
                "haploid dosage {dosage}: fitted {}, truth {truth}",
                fitted.get()
            );
        }
        for (dosage, (fitted, truth)) in fitted_diploid
            .by_alt_copies()
            .iter()
            .zip(&TRUTH)
            .enumerate()
        {
            assert!(
                (fitted.get() - truth).abs() < 0.01 * truth,
                "diploid dosage {dosage}: fitted {}, truth {truth}",
                fitted.get()
            );
        }
    }

    /// **Every iterate is scored at its own rates**, not at the rates it started the round
    /// with — which is what makes the scores a comparison between iterates rather than a
    /// mixture of two rounds.
    ///
    /// Checked by re-deriving each score from the outside: replay the alternation, and at
    /// each round score the whole-sample table at the rungs that round *produced*. Nothing
    /// in the file pinned this, and a loop scoring at the previous rungs is a plausible
    /// off-by-one — it converges to the same place and reports different scores.
    #[test]
    fn every_iterates_score_is_taken_at_the_rates_that_round_produced() {
        let world = TwoLibraryWorld::build();
        let start = world.three_times_the_truth();
        let (_, trace) = world.fit_from(&start, MAX_COUPLED_FIT_ITERATIONS);

        let shares = library_shares(&world.read_group_histograms);
        let cells_of_ploidy: BTreeMap<Ploidy, Vec<Cell>> = world
            .whole_sample
            .iter()
            .map(|(&ploidy, table)| (ploidy, table.cells(ploidy)))
            .collect();
        let all_cells: Vec<Cell> = cells_of_ploidy.values().flatten().cloned().collect();

        let mut rungs = start;
        for (round, reported) in trace.iter().map(|iterate| iterate.score).enumerate() {
            let genotype_frequencies = climb_frequencies(
                &cells_of_ploidy,
                &noise_from(&shares, &rungs, &world.ladder),
            );
            let fitted = fit_read_group_error_rates(
                &world.read_group_histograms,
                &genotype_frequencies,
                &world.ladder,
            );
            let next: BTreeMap<ReadGroupId, usize> = fitted
                .iter()
                .map(|(&group, fit)| (group, fit.rung))
                .collect();

            let at_this_rounds_rates = whole_sample_score(
                &all_cells,
                &noise_from(&shares, &next, &world.ladder),
                &genotype_frequencies,
            )
            .get();
            let at_last_rounds_rates = whole_sample_score(
                &all_cells,
                &noise_from(&shares, &rungs, &world.ladder),
                &genotype_frequencies,
            )
            .get();

            assert!(
                (reported - at_this_rounds_rates).abs() < 1e-6,
                "round {round}: the loop reported {reported}, and the whole-sample table \
                 scores {at_this_rounds_rates} at that round's own rates"
            );
            if round == 0 {
                assert!(
                    (at_this_rounds_rates - at_last_rounds_rates).abs() > 1.0,
                    "round 0 has to separate the two rate sets or this asserts nothing: \
                     {at_this_rounds_rates} against {at_last_rounds_rates}"
                );
            }
            rungs = next;
        }
    }

    /// **`nearest_rung` breaks a tie in log space and not in probability**, and the case
    /// that separates them is narrow enough to be worth writing down: on a ladder stepping
    /// by a quarter of a Phred, two adjacent rungs' geometric and arithmetic midpoints sit
    /// 4 parts in 10,000 apart — a six-hundredth of a rung. A rate inside that band is
    /// nearest rung 80 in log space and rung 81 in probability.
    ///
    /// The existing test's example — the geometric midpoint times 1.001 — is *above* the
    /// band, where both metrics answer 80, so it does not separate them.
    #[test]
    fn nearest_rung_breaks_the_tie_in_log_space_and_not_in_probability() {
        let ladder = error_rate_ladder();
        let geometric = (ladder[80].get() * ladder[81].get()).sqrt();
        let arithmetic = 0.5 * (ladder[80].get() + ladder[81].get());
        assert!(geometric < arithmetic);

        let inside_the_band = geometric * 1.0002;
        assert!(inside_the_band < arithmetic, "the band is 4e-4 wide");

        let in_probability = ladder
            .iter()
            .enumerate()
            .min_by(|(_, left), (_, right)| {
                (left.get() - inside_the_band)
                    .abs()
                    .total_cmp(&(right.get() - inside_the_band).abs())
            })
            .map(|(rung, _)| rung)
            .expect("the ladder is not empty");

        assert_eq!(nearest_rung(&ladder, inside_the_band), 80);
        assert_eq!(in_probability, 81, "the two metrics have to disagree here");
    }
}
