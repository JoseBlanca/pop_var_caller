//! **What does it cost to read the prior's two numbers off reads instead of off genotypes — and
//! does the population curve the fit assumes pull them toward itself?**
//!
//! # The question, and why it is the one that can kill the proposal
//!
//! `doc/devel/ng/research/ordinary_site_prior_moments.md` proposes reading the genotype prior's two
//! numbers — the population's **mean alternative-allele frequency** and its **heterozygosity** —
//! straight off the census positions:
//!
//! ```text
//! heterozygosity   pi  =  mean over positions of   2 k (2N - k) / (2N (2N - 1))
//! mean frequency   f   =  mean over positions of   k / 2N
//! ```
//!
//! **Nobody can count `k`.** At three reads a position a heterozygote often shows only one of its
//! two alleles, and a sequencing error often looks like a third; the first hides variation and the
//! second invents it, and they do not cancel. So `k` has to be an *expected* count under the read
//! model, taken from the joint fit's own per-position posteriors.
//!
//! Those posteriors are computed under the fit's description of the population: **one Beta over
//! what segregates, plus a spike of positions carrying only the reference base and a spike carrying
//! only a non-reference one**. A moment computed from them therefore inherits a pull toward that
//! description. At a hundred reads a position the reads decide and the pull vanishes; at three they
//! do not.
//!
//! **A population drawn from that same family cannot show the pull**, because there is nothing to
//! be pulled toward — which is why half the populations here are outside it.
//!
//! # The three ways of computing each moment, and what separates them
//!
//! Every cell computes both moments four times over **the same drawn cohort**:
//!
//! 1. **From the genotypes the cohort was drawn with.** No reads in it. This is the ceiling: it is
//!    what the estimator would give if sequencing were perfect.
//! 2. **From the posteriors, the obvious way.** Replace `k` by its posterior mean and evaluate the
//!    two formulas. The second formula is not linear in `k`, so this is not the posterior mean of
//!    the formula — see 3.
//! 3. **From the posteriors, with the term the second formula's curvature needs.**
//!    `E[k(2N - k)] = 2N·E[k] - E[k]² - Var(k)`, and dropping `Var(k)` makes the heterozygosity
//!    come back **high**. The variance used here is the sum of the samples' own, which is exact
//!    only if the samples' posteriors are independent given the reads; they are not, so a residual
//!    remains and its size is what column 3 minus column 1 shows.
//! 4. **The fit's own population curve, integrated in closed form.** Its mean frequency is
//!    `p_fixed_alt + p_segregating · a/(a + b)` and its heterozygosity
//!    `p_segregating · 2ab/((a + b)(a + b + 1))` — two integrals, no census pass and no panel.
//!    **This is what the caller does since 2026-08-27**, and the column that matters is how many
//!    calls it moves against the census average.
//!
//! **A fifth arm used to be here and its machinery is deleted**: the curve evaluated into the
//! panel's `2N + 1` allele-count classes and searched for a two-parameter pair, which is the
//! detour this work removed. Its figures are in report §9.
//!
//! **The gap between 1 and 3 is what the read model costs.** The gap between 3 and 1 *on the
//! populations outside the fit's family, at low depth, where it does not appear on the populations
//! inside it* is the pull this program exists to find.
//!
//! # What each arm is scored against
//!
//! **Every arm is scored against the moments of the positions that cohort was itself drawn at** —
//! not against the population the positions came from. Scoring against the population would add
//! the census's own scatter, which is 3% to 16% of the mean at these sizes, to every error column
//! and call it a property of the estimator. The census's scatter is measured separately, in
//! `ng_prior_moment_estimators.rs`.
//!
//! # What is held fixed
//!
//! Within one population and one depth the positions, the frequencies and every sample's genotypes
//! and reads are drawn **once**, at the largest panel; each panel-size arm refits the first `N`
//! samples of that same draw. Nothing moves across the panel-size arms for a reason unrelated to
//! panel size.
//!
//! # Mismapped positions
//!
//! A position where two stretches of genome the reference holds once both pile their reads up
//! reads part non-reference **in every sample**, which is what a heterozygous position looks like.
//! The fit already produces a per-position posterior that a position is that. A share of positions
//! is planted as such in the second half of this program, and the heterozygosity is computed with
//! and without weighting each position by one minus that posterior.
//!
//! Run: `./scripts/dev.sh cargo run --release --example ng_prior_moments_from_reads`

use std::collections::BTreeMap;
use std::env;
use std::time::Instant;

use pop_var_caller::ng::calling::genotype_prior::{GenotypePriorModel, MarginalizedDirichletPrior};
use pop_var_caller::ng::calling::genotype_prior::{PriorRow, SpectrumSeed};
use pop_var_caller::ng::calling::genotype_prior::{
    VariantClass, fill_locus_concentration, seed_from_population_moments,
};
use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, NamedReadGroup, ObservedAllele, PackedDepthCodes, ReadCap, RecordingTerms,
    SampleCensusEvidence, Section, SectionKey, SelectionTermsDigest,
};
use pop_var_caller::ng::parameter_estimation::joint::census_moments::CensusMoments;
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFitConfig, StartingPoint, fit_jointly,
};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::{
    AlleleId, ExpectedAlternativeFrequency, ExpectedHeterozygosity, InbreedingF, LogProb,
    ReadGroupId,
};

/// How often a read misreads a base at an ordinary position.
const CLEAN_ERROR_RATE: f64 = 0.002;

/// **The census keeps at most this many reads at a position, and a run at a hundred reads a
/// position meets it.** The store's allele counts are one byte and its depth ladder is exact to
/// 124 reads and coarse above (`generic::depth_bins`), so a run subsamples each position down to
/// the cap before recording it (`DepthCap`). Modelling that here is not a convenience: without it
/// the drawn depths at 100 reads a position spill into the ladder's first widening bin, which
/// stands for 35 depths, and the fit refuses a code that wide.
///
/// Drawing exactly this many reads is the same distribution as drawing more and thinning down to
/// it, because thinning a binomial sample gives a binomial with the same probability.
const CENSUS_DEPTH_CAP: u32 = 124;

/// Panel sizes the read sweep runs at. **Shorter than the estimator sweep's**, because every arm
/// here is a full joint fit rather than an average over an array.
const PANELS: [usize; 3] = [1, 10, 63];

/// Reads a sample at a position — the committed range's two ends and two points between
/// (`CLAUDE.md`, *what this caller has to work on*).
const DEPTHS: [f64; 4] = [3.0, 8.0, 20.0, 100.0];

fn main() {
    let mut args = env::args().skip(1);
    let positions: usize = args.next().map_or(20_000, |a| a.parse().expect("a count"));
    let inbreeding: f64 = args
        .next()
        .map_or(0.15, |a| a.parse().expect("a coefficient"));
    let mismapped_share: f64 = args.next().map_or(0.01, |a| a.parse().expect("a share"));
    let replicates: usize = args.next().map_or(3, |a| a.parse().expect("a count"));
    // Which populations to run, as indices into `populations()`, and whether to run the mismapped
    // half. **Both exist so a re-run can add one population without repeating twenty minutes of
    // arithmetic whose answer is already recorded**; the seeds are fixed, so a repeated cell
    // returns the identical numbers.
    let chosen: Vec<usize> = args.next().map_or_else(
        || (0..populations().len()).collect(),
        |a| {
            a.split(',')
                .map(|part| part.parse().expect("a population index"))
                .collect()
        },
    );
    // `mismapped` (the default) runs both halves, `nomismapped` only the first, `mismappedonly`
    // only the second — so a re-run that only needs the planted half does not repeat twenty
    // minutes of arithmetic whose answer is already recorded.
    let half = args.next().unwrap_or_else(|| "mismapped".to_string());
    let run_main = half != "mismappedonly";
    let run_mismapped = half != "nomismapped";
    let mismapped_error_rate: f64 = args
        .next()
        .map_or(DEFAULT_MISMAPPED_ERROR_RATE, |a| a.parse().expect("a rate"));

    println!("# Reading the prior's two numbers off reads instead of off genotypes");
    println!();
    println!(
        "{positions} census positions and {replicates} independently drawn cohorts a cell, \
         homozygote excess {inbreeding} in every drawn individual, panel sizes {PANELS:?}, depths \
         {DEPTHS:?} reads a position."
    );
    println!(
        "Error rate {CLEAN_ERROR_RATE} at an ordinary position. Every arm is scored against the \
         moments of the positions that cohort was drawn at."
    );

    for (index, population) in populations().into_iter().enumerate() {
        if !run_main || !chosen.contains(&index) {
            continue;
        }
        println!();
        println!("## {}", population.name);
        println!();
        println!(
            "{} — heterozygosity {:.6} ({:.2} per kilobase), mean frequency {:.6}.",
            population.family_note(),
            population.heterozygosity(),
            1_000.0 * population.heterozygosity(),
            population.mean_frequency()
        );
        for &depth in &DEPTHS {
            run_depth(
                &population,
                depth,
                inbreeding,
                positions,
                replicates,
                0.0,
                mismapped_error_rate,
            );
        }
    }

    if !run_mismapped {
        return;
    }
    println!();
    println!("# Mismapped positions");
    println!();
    println!(
        "The same two populations with {:.1} in 100 positions planted as mismapped: at those, \
         every sample's reads disagree with the reference at {mismapped_error_rate} instead of \
         {CLEAN_ERROR_RATE}, whatever its genotype. The heterozygosity is then computed twice — \
         once over all positions, once weighting each position by one minus the fit's own \
         posterior that it is mismapped.",
        mismapped_share * 100.0
    );
    for (index, population) in populations().into_iter().enumerate() {
        if !run_mismapped || !chosen.contains(&index) {
            continue;
        }
        println!();
        println!("## {} — with mismapped positions", population.name);
        for &depth in &[3.0_f64, 20.0] {
            run_depth(
                &population,
                depth,
                inbreeding,
                positions,
                replicates,
                mismapped_share,
                mismapped_error_rate,
            );
        }
    }
}

/// How often a read at a mismapped position disagrees with the reference, unless the command line
/// says otherwise.
///
/// **Two stretches of genome piling up at one place differ from each other at a few percent of
/// their bases, so this is the mild end of the phenomenon**, and it is the rate
/// `examples/ng_joint_sample_count_sweep.rs` plants. At three reads a position it leaves a
/// disagreeing read at about one planted position in six, which is not far from what an ordinary
/// position shows. **The rate is a command-line argument because the question this half asks —
/// whether weighting positions away by the fit's own mismapped posterior recovers the
/// heterozygosity — only has teeth at a rate high enough that the unweighted estimate fails.**
const DEFAULT_MISMAPPED_ERROR_RATE: f64 = 0.06;

// ---------------------------------------------------------------------------------------------
// The populations
// ---------------------------------------------------------------------------------------------

/// One Beta over the segregating positions, at a share of the segregating whole.
#[derive(Clone, Copy)]
struct Component {
    weight: f64,
    a: f64,
    b: f64,
}

/// A population, as this program draws one.
///
/// **Two of them, and the second is the load-bearing one**: a mixture of two Betas is a shape the
/// joint fit's single Beta cannot hold, so a posterior computed under the fit is being pulled
/// toward something the population is not. On the first population there is nothing to pull
/// toward and any such pull is invisible.
struct Population {
    name: &'static str,
    invariant: f64,
    fixed_alt: f64,
    segregating: Vec<Component>,
    inside_the_fitted_family: bool,
}

impl Population {
    fn share_segregating(&self) -> f64 {
        1.0 - self.invariant - self.fixed_alt
    }

    fn family_note(&self) -> &'static str {
        if self.inside_the_fitted_family {
            "A shape the joint fit can hold exactly"
        } else {
            "A shape the joint fit cannot hold: two Betas, where it has one"
        }
    }

    fn mean_frequency(&self) -> f64 {
        self.fixed_alt
            + self.share_segregating()
                * self
                    .segregating
                    .iter()
                    .map(|c| c.weight * c.a / (c.a + c.b))
                    .sum::<f64>()
    }

    fn heterozygosity(&self) -> f64 {
        self.share_segregating()
            * self
                .segregating
                .iter()
                .map(|c| c.weight * 2.0 * c.a * c.b / ((c.a + c.b) * (c.a + c.b + 1.0)))
                .sum::<f64>()
    }

    fn draw_frequency(&self, rng: &mut Rng) -> f64 {
        match rng.pick(&[self.invariant, self.fixed_alt, self.share_segregating()]) {
            0 => 0.0,
            1 => 1.0,
            _ => {
                let mut cut = rng.uniform();
                for component in &self.segregating {
                    cut -= component.weight;
                    if cut <= 0.0 {
                        return rng.beta(component.a, component.b);
                    }
                }
                let last = self.segregating.last().expect("a population segregates");
                rng.beta(last.a, last.b)
            }
        }
    }
}

/// The two populations.
///
/// **Both segregate at 2 positions in 100, which is more than tomato's 1 in 200 and is a
/// deliberate choice about what this program can afford.** What sets how well the fit resolves the
/// population is the number of *segregating* positions it sees, and a real run gets about ten
/// thousand of them — two million census positions at tomato's rate
/// (`parameter_prepass_census_sites.md` §5.1). Every arm here is a full joint fit, so two million
/// positions is out of reach; at 2 in 100 a twenty-thousand-position census carries four hundred
/// segregating positions, which is a twenty-fifth of a real run's rather than a two-hundred-and-
/// fiftieth. **The diversity level is not what this program measures** — what it measures is the
/// gap between reading the moments off genotypes and reading them off reads, and both arms see the
/// same positions. The level's own effect is in `ng_prior_moment_estimators.rs`, at tomato's rate.
fn populations() -> Vec<Population> {
    vec![
        Population {
            name: "nearly all alternative alleles rare — a shape the fit can hold",
            invariant: 0.9780,
            fixed_alt: 0.0020,
            segregating: vec![Component {
                weight: 1.0,
                a: 0.20,
                b: 1.00,
            }],
            inside_the_fitted_family: true,
        },
        Population {
            name: "two peaks, off centre — a shape the fit cannot hold",
            invariant: 0.9790,
            fixed_alt: 0.0010,
            segregating: vec![
                Component {
                    weight: 0.70,
                    a: 30.0,
                    b: 70.0,
                },
                Component {
                    weight: 0.30,
                    a: 90.0,
                    b: 10.0,
                },
            ],
            inside_the_fitted_family: false,
        },
        // **The control for the population above, and without it that comparison proves nothing.**
        // The two-peaked population differs from the rare-allele one in two ways at once: it is
        // outside the fit's family, and its alternative alleles sit at frequencies of 0.3 and 0.9
        // where a few reads settle a genotype easily, rather than piled up near zero where they do
        // not. A gap between those two could be either. **This population keeps the frequencies
        // and drops the second peak**: one `Beta(30, 70)`, squarely inside the family, at a
        // segregating share chosen so its heterozygosity lands within 5% of the two-peaked one's.
        // A pull toward the fit's own family must show up as a gap between *these* two.
        Population {
            name: "one peak at the same place — the control, inside the family",
            invariant: 0.9816,
            fixed_alt: 0.0010,
            segregating: vec![Component {
                weight: 1.0,
                a: 30.0,
                b: 70.0,
            }],
            inside_the_fitted_family: true,
        },
    ]
}

// ---------------------------------------------------------------------------------------------
// One depth: draw a cohort, refit at each panel size, repeat
// ---------------------------------------------------------------------------------------------

/// **Every cell is several independently drawn cohorts, and every ratio is printed with the
/// spread across them.** A single drawn cohort at these sizes carries about four hundred
/// segregating positions, and the ratios wander by ten to fifteen parts in a hundred from one
/// draw to the next — enough to make a depth trend appear that is not there. The `+-` figure
/// beside each ratio is how precisely this run pins it: one standard deviation over the
/// replicates, divided by the square root of their number.
fn run_depth(
    population: &Population,
    depth: f64,
    inbreeding: f64,
    positions: usize,
    replicates: usize,
    mismapped_share: f64,
    mismapped_error_rate: f64,
) {
    println!();
    println!("### {depth} reads a position");
    println!();
    println!(
        "**The columns marked `/gt` are the number the plan asks for**: what reading the moment \
         off reads gives, divided by what the same positions' own genotypes give. A 1.000 there \
         means sequencing cost nothing. `genotypes /drawn` is the genotype arm against the \
         positions it was drawn at — the panel's own sampling, which is \
         `ng_prior_moment_estimators.rs`'s subject and not this one's."
    );
    println!();
    println!(
        "The genotype columns call one sample's genotype at {GENOTYPE_COMPARISON_LOCI} freshly \
         drawn loci under each seed, counted over the loci that segregate. **`trebled` is the \
         control**: the same comparison with the direct seed's alternative concentration \
         multiplied by three. If that column is also zero the comparison cannot detect movement \
         and the others mean nothing."
    );
    println!();
    println!(
        "| individuals | freq: genotypes /drawn | freq: posteriors /gt | het: genotypes /drawn | \
         het: posteriors plain /gt | het: posteriors + variance /gt | het: the curve's own /gt | \
         freq: the curve's own /gt | **calls moved, curve vs census** | calls moved, trebled | \
         alpha_alt curve | alpha_alt direct | seconds |"
    );
    println!("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|");

    let largest = *PANELS.last().expect("a panel");
    let mut cells: Vec<Cell> = PANELS.iter().map(|_| Cell::default()).collect();
    let mut planted_cells: Vec<MismappedCell> =
        PANELS.iter().map(|_| MismappedCell::default()).collect();
    let started = Instant::now();

    for replicate in 0..replicates {
        // **The positions and every sample's reads are drawn once a replicate**, and each
        // panel-size arm refits the first `N` samples of that same draw, so nothing moves across
        // the arms for a reason unrelated to panel size.
        let cohort = draw(
            population,
            largest,
            positions,
            depth,
            inbreeding,
            mismapped_share,
            mismapped_error_rate,
            mix(
                0x9E37_79B9_7F4A_7C15,
                replicate as u64,
                (depth as u64) * 7919 + population.name.len() as u64,
            ),
        );
        if replicate == 0 {
            println!(
                "\nFirst replicate: mean frequency {:.6}, heterozygosity {:.6}, {} of the \
                 {positions} positions segregating.\n",
                cohort.drawn_mean_frequency,
                cohort.drawn_heterozygosity,
                cohort.segregating_positions
            );
        }

        for (index, &individuals) in PANELS.iter().enumerate() {
            let subset: Vec<SampleCensusEvidence> = cohort.samples[..individuals].to_vec();
            let mut evidence =
                CohortCensusEvidence::new(subset).expect("a drawn cohort records one way");
            let config = JointFitConfig {
                quadrature_nodes: 12,
                starting_points: StartingPoint::spanning_the_class_separation(),
                genotype_posteriors: true,
                ..JointFitConfig::default()
            };
            let fit = fit_jointly(&mut evidence, &config).expect("a drawn cohort pools");

            let chromosomes = 2.0 * individuals as f64;
            let panel_inbreeding =
                InbreedingF::try_new(inbreeding.min(0.99)).expect("a coefficient");

            // ---- from the genotypes the cohort was drawn with ---------------------------------
            let mut oracle_frequency = 0.0_f64;
            let mut oracle_heterozygosity = 0.0_f64;
            let mut oracle_over_unplanted_total = 0.0_f64;
            let mut unplanted = 0_u64;
            for position in 0..positions {
                let copies: u32 = cohort.genotypes[position][..individuals]
                    .iter()
                    .map(|g| u32::from(*g))
                    .sum();
                oracle_frequency += f64::from(copies) / chromosomes;
                let here = nei_heterozygosity(f64::from(copies), 0.0, chromosomes);
                oracle_heterozygosity += here;
                if !cohort.planted_mismapped[position] {
                    oracle_over_unplanted_total += here;
                    unplanted += 1;
                }
            }
            oracle_frequency /= positions as f64;
            oracle_heterozygosity /= positions as f64;
            let oracle_over_unplanted = oracle_over_unplanted_total / unplanted.max(1) as f64;

            let posterior =
                MomentsFromPosteriors::of(&fit.genotype_posterior, individuals, positions);
            // **The shipped estimator, over the same posteriors, checked against this program's
            // own reduction.** The implementation plan's verification table asks that this harness
            // reproduce the report's figures *against the wired-in estimators*, and until this
            // check existed it did not: the reduction below is a copy, so a change to the library's
            // arithmetic would have moved nothing here and every number this program prints would
            // have gone on describing code that no longer ran.
            //
            // The copy is kept rather than replaced, because two of the three columns it produces
            // are quantities the library deliberately does not compute — the heterozygosity
            // *without* the curvature term (report §4.1's subject) and the one weighted by each
            // position's mismapped posterior (§6's). What is checked is the third, which is the
            // one the library ships.
            check_against_the_shipped_estimator(
                &fit.genotype_posterior,
                individuals,
                positions,
                &posterior,
            );
            // **Both numbers integrated off the fitted curve in closed form**, with no projection,
            // no search and no census average. This is the two-line repair, and the column it
            // feeds asks the question the recommendation turns on — whether choosing between it
            // and the census average moves a genotype at all.
            let curve = FromTheCurve::of(&fit.density.value, fit.expected_heterozygosity);
            let direct_seed =
                seed_from_moments(posterior.frequency, posterior.heterozygosity_with_variance);
            let curve_against_census = genotype_calls(
                curve.seed,
                direct_seed,
                panel_inbreeding,
                population,
                depth,
                inbreeding,
                GENOTYPE_COMPARISON_LOCI,
                replicate,
                individuals,
            );

            let cell = &mut cells[index];
            cell.oracle_frequency_over_drawn
                .add(oracle_frequency / cohort.drawn_mean_frequency);
            cell.posterior_frequency_over_oracle
                .add(posterior.frequency / oracle_frequency);
            cell.oracle_heterozygosity_over_drawn
                .add(oracle_heterozygosity / cohort.drawn_heterozygosity);
            cell.posterior_plain_over_oracle
                .add(posterior.heterozygosity_plain / oracle_heterozygosity);
            cell.posterior_with_variance_over_oracle
                .add(posterior.heterozygosity_with_variance / oracle_heterozygosity);
            cell.curve_heterozygosity_over_oracle
                .add(curve.heterozygosity / oracle_heterozygosity);
            cell.curve_frequency_over_oracle
                .add(curve.frequency / oracle_frequency);
            cell.calls_moved_curve_against_census
                .add(curve_against_census.moved_in_a_hundred());
            cell.alpha_alt_curve.add(curve.seed.alpha_alt_total());
            cell.calls_moved_when_trebled
                .add(curve_against_census.moved_when_trebled_in_a_hundred());
            cell.alpha_alt_direct.add(direct_seed.alpha_alt_total());

            if mismapped_share > 0.0 {
                let weighted = MomentsFromPosteriors::weighted_by_being_ordinary(
                    &fit.genotype_posterior,
                    &fit.noisy_posterior,
                    individuals,
                    positions,
                );
                // **Both scored against the same panel's own genotypes, and that is a correction
                // a review caught.** The numerator is Nei's panel estimator, whose expectation is
                // the population value times `1 - F/(2N-1)`; the drawn positions' `2 f (1 - f)`
                // carries no such factor. Dividing one by the other left 15% of that factor
                // standing at one sample with a homozygote excess of 0.15 — larger than anything
                // mismapping does here, and it was being read as the cost of mismapping. Against
                // the oracle it cancels exactly, as it does in every other table in this program.
                let planted_cell = &mut planted_cells[index];
                planted_cell
                    .unweighted_over_oracle
                    .add(posterior.heterozygosity_with_variance / oracle_heterozygosity);
                planted_cell
                    .weighted_over_oracle
                    .add(weighted.heterozygosity_with_variance / oracle_heterozygosity);
                // **A control on the plant itself.** A planted position here carries real
                // variation drawn from the same population and only noisier reads, so the
                // genotypes over every position and over the unplanted ones should agree. If they
                // do not, the plant moved the truth rather than the reading of it.
                planted_cell
                    .oracle_over_all_against_oracle_over_unplanted
                    .add(oracle_heterozygosity / oracle_over_unplanted);
            }
        }
    }

    let seconds = started.elapsed().as_secs_f64() / (replicates * PANELS.len()) as f64;
    for (index, &individuals) in PANELS.iter().enumerate() {
        let cell = &cells[index];
        println!(
            "| {individuals} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.2e} | \
             {:.2e} | {seconds:.0} |",
            cell.oracle_frequency_over_drawn.as_ratio(),
            cell.posterior_frequency_over_oracle.as_ratio(),
            cell.oracle_heterozygosity_over_drawn.as_ratio(),
            cell.posterior_plain_over_oracle.as_ratio(),
            cell.posterior_with_variance_over_oracle.as_ratio(),
            cell.curve_heterozygosity_over_oracle.as_ratio(),
            cell.curve_frequency_over_oracle.as_ratio(),
            cell.calls_moved_curve_against_census.as_percentage(),
            cell.calls_moved_when_trebled.as_percentage(),
            cell.alpha_alt_curve.mean,
            cell.alpha_alt_direct.mean,
        );
    }

    if mismapped_share > 0.0 {
        println!();
        println!(
            "**Both columns divided by the same panel's own genotypes**, so these rows compare \
             directly with the no-plant sweep and the inbreeding factor cancels. The last column \
             is a control on the plant: the genotypes over every position against the genotypes \
             over the unplanted ones, which is 1.000 if planting changed the reads and not the \
             truth."
        );
        println!();
        println!(
            "| individuals | over all positions /gt | weighted by one minus the mismapped \
             posterior /gt | control: genotypes, all positions over unplanted |"
        );
        println!("|---:|---:|---:|---:|");
        for (index, &individuals) in PANELS.iter().enumerate() {
            println!(
                "| {individuals} | {} | {} | {} |",
                planted_cells[index].unweighted_over_oracle.as_ratio(),
                planted_cells[index].weighted_over_oracle.as_ratio(),
                planted_cells[index]
                    .oracle_over_all_against_oracle_over_unplanted
                    .as_ratio(),
            );
        }
    }
}

/// One panel size's numbers, accumulated over the replicates.
///
/// **The planted-positions half of this program keeps its three numbers in [`MismappedCell`]
/// rather than borrowing three of these fields.** An earlier version borrowed
/// `posterior_with_variance_over_oracle` and `curve_heterozygosity_over_oracle` to hold quantities
/// that were neither over the oracle nor anything to do with the caller's own path, which told
/// a reader auditing the scoring — the one audit this program most needs — the wrong thing.
#[derive(Default)]
struct Cell {
    oracle_frequency_over_drawn: Welford,
    posterior_frequency_over_oracle: Welford,
    oracle_heterozygosity_over_drawn: Welford,
    posterior_plain_over_oracle: Welford,
    posterior_with_variance_over_oracle: Welford,
    /// The heterozygosity integrated off the fitted curve, over the same panel's genotypes.
    curve_heterozygosity_over_oracle: Welford,
    /// The mean frequency integrated off the fitted curve in closed form, over the same panel's
    /// genotypes — the route that keeps the curve and drops only the projection and the search.
    curve_frequency_over_oracle: Welford,
    /// **The recommendation's own question**: how many calls come out differently under the two
    /// routes that both remove the search — the curve integrated in closed form, and the census
    /// average.
    calls_moved_curve_against_census: Welford,
    alpha_alt_curve: Welford,
    calls_moved_when_trebled: Welford,
    alpha_alt_direct: Welford,
}

/// The planted-positions half's three numbers at one panel size, over the replicates.
#[derive(Default)]
struct MismappedCell {
    /// The heterozygosity from the posteriors over every position, over the same panel's genotypes.
    unweighted_over_oracle: Welford,
    /// …and with each position weighted by one minus the fit's posterior that it is mismapped.
    weighted_over_oracle: Welford,
    /// The control: the genotypes over every position against the genotypes over the unplanted
    /// ones. Near 1.000 means the plant changed the reads and not the truth.
    oracle_over_all_against_oracle_over_unplanted: Welford,
}

/// Running mean and spread, so a sweep does not hold every replicate's number.
#[derive(Clone, Copy, Default)]
struct Welford {
    count: f64,
    mean: f64,
    sum_of_squared_deviations: f64,
}

impl Welford {
    fn add(&mut self, value: f64) {
        self.count += 1.0;
        let delta = value - self.mean;
        self.mean += delta / self.count;
        self.sum_of_squared_deviations += delta * (value - self.mean);
    }

    /// How precisely this run pins the mean: one standard deviation over the replicates, divided
    /// by the square root of their number. **A departure from 1.000 smaller than this is a
    /// departure the run cannot see.**
    fn uncertainty_of_the_mean(&self) -> f64 {
        if self.count < 2.0 {
            return f64::NAN;
        }
        (self.sum_of_squared_deviations / (self.count - 1.0) / self.count).sqrt()
    }

    fn as_ratio(&self) -> String {
        format!("{:.3}±{:.3}", self.mean, self.uncertainty_of_the_mean())
    }

    fn as_percentage(&self) -> String {
        format!("{:.2}±{:.2}%", self.mean, self.uncertainty_of_the_mean())
    }
}

/// Fold three numbers into one seed. Neighbouring seeds in this generator's state produce
/// correlated first draws, so replicates are separated by mixing rather than by adding.
fn mix(base: u64, first: u64, second: u64) -> u64 {
    let mut state = base ^ first.wrapping_mul(0xFF51_AFD7_ED55_8CCD);
    state ^= state >> 33;
    state = state.wrapping_mul(0xC4CE_B9FE_1A85_EC53) ^ second.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    state ^= state >> 29;
    state | 1
}

/// How many single-sample loci the genotype comparison draws.
const GENOTYPE_COMPARISON_LOCI: usize = 200_000;

// ---------------------------------------------------------------------------------------------
// The moments, three ways
// ---------------------------------------------------------------------------------------------

/// **How often two chromosomes drawn at random from the panel differ**, at one position, from the
/// expected alternative-copy count and its variance.
///
/// Of the `2N (2N - 1)` ordered pairs of distinct chromosomes, `2 k (2N - k)` have one carrying the
/// alternative allele and the other not. When `k` is known this is `2 k (2N - k)`. When only its
/// posterior mean is known it is **not**: the expression is quadratic, so
/// `E[k (2N - k)] = 2N·E[k] - E[k]² - Var(k)`, and passing `0.0` for the variance is the plain
/// substitution — which comes back **high**, by exactly the variance the reads left behind.
fn nei_heterozygosity(expected_copies: f64, variance: f64, chromosomes: f64) -> f64 {
    debug_assert!(
        chromosomes >= 2.0,
        "two chromosomes are needed before a pair of them can differ"
    );
    2.0 * (chromosomes * expected_copies - expected_copies * expected_copies - variance)
        / (chromosomes * (chromosomes - 1.0))
}

/// The two moments, computed from the fit's per-position per-sample genotype posteriors.
struct MomentsFromPosteriors {
    frequency: f64,
    /// `k` replaced by its posterior mean in both formulas.
    heterozygosity_plain: f64,
    /// …and with the variance term the second formula's curvature needs.
    heterozygosity_with_variance: f64,
}

impl MomentsFromPosteriors {
    /// The posteriors are three numbers a sample a position, in position order: heterozygous,
    /// both copies non-reference, and carrying an extra copy of the position
    /// (`JointFit::genotype_posterior`). A carrier is not a genotype and takes no part here.
    fn of(genotype_posterior: &[f32], individuals: usize, positions: usize) -> Self {
        Self::over(genotype_posterior, individuals, positions, |_| 1.0)
    }

    /// The same, with each position weighted by one minus the fit's posterior that it is
    /// mismapped — so a position where every sample reads part non-reference counts for almost
    /// nothing.
    fn weighted_by_being_ordinary(
        genotype_posterior: &[f32],
        noisy_posterior: &[f32],
        individuals: usize,
        positions: usize,
    ) -> Self {
        Self::over(genotype_posterior, individuals, positions, |position| {
            1.0 - f64::from(noisy_posterior[position])
        })
    }

    fn over(
        genotype_posterior: &[f32],
        individuals: usize,
        positions: usize,
        weight_of: impl Fn(usize) -> f64,
    ) -> Self {
        assert_eq!(
            genotype_posterior.len(),
            positions * individuals * 3,
            "the fit returns three numbers a sample a position: {positions} positions × \
             {individuals} samples × 3 is {}, and the fit returned {}",
            positions * individuals * 3,
            genotype_posterior.len()
        );
        let chromosomes = 2.0 * individuals as f64;
        let mut frequency = 0.0_f64;
        let mut plain = 0.0_f64;
        let mut with_variance = 0.0_f64;
        let mut total_weight = 0.0_f64;
        for position in 0..positions {
            let weight = weight_of(position);
            let base = position * individuals * 3;
            let mut expected_copies = 0.0_f64;
            let mut variance = 0.0_f64;
            for sample in 0..individuals {
                let heterozygous = f64::from(genotype_posterior[base + sample * 3]);
                let homozygous_alt = f64::from(genotype_posterior[base + sample * 3 + 1]);
                let copies = heterozygous + 2.0 * homozygous_alt;
                let copies_squared = heterozygous + 4.0 * homozygous_alt;
                expected_copies += copies;
                variance += copies_squared - copies * copies;
            }
            frequency += weight * expected_copies / chromosomes;
            plain += weight * nei_heterozygosity(expected_copies, 0.0, chromosomes);
            with_variance += weight * nei_heterozygosity(expected_copies, variance, chromosomes);
            total_weight += weight;
        }
        let total_weight = total_weight.max(f64::MIN_POSITIVE);
        Self {
            frequency: frequency / total_weight,
            heterozygosity_plain: plain / total_weight,
            heterozygosity_with_variance: with_variance / total_weight,
        }
    }
}

/// **What the caller reads off the fit: the two integrals of its own population curve**, and the
/// seed they imply.
///
/// **Until 2026-08-27 this held a third number** — the mean frequency of the pair the caller got by
/// evaluating the curve into the panel's `2N + 1` allele-count classes and searching for a
/// two-parameter match. That search is deleted
/// (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §5) and its figures are in report §9.
struct FromTheCurve {
    /// **The mean allele frequency integrated straight off the fitted curve**, with no projection
    /// and no search: `p_fixed_alt + p_segregating · a/(a + b)`.
    ///
    /// **The prior takes two numbers and the fitted curve carries four, and the argument that
    /// started this work was that compressing four into two must bias the prior.** It need not:
    /// the two numbers the prior wants *are* two integrals of the curve, and both have a closed
    /// form — this one, and `FrequencyDensity::expected_heterozygosity` for the other.
    frequency: f64,
    heterozygosity: f64,
    seed: SpectrumSeed,
}

impl FromTheCurve {
    fn of(density: &FrequencyDensity, expected_heterozygosity: f64) -> Self {
        let diversity = ExpectedHeterozygosity::try_new(expected_heterozygosity.clamp(0.0, 0.5))
            .expect("a clamped heterozygosity is a probability");
        let frequency = density.expected_alternative_frequency();
        let seed = seed_from_population_moments(
            ExpectedAlternativeFrequency::try_new(frequency).ok(),
            Some(diversity),
        );
        Self {
            frequency,
            heterozygosity: expected_heterozygosity,
            seed,
        }
    }
}

/// **Refuse to print a number this program computed itself if the shipped estimator disagrees with
/// it.**
///
/// `MomentsFromPosteriors` below is this program's own reduction, and it exists because two of the
/// three quantities it produces are ones the library deliberately does not: the heterozygosity with
/// the curvature term *dropped*, which is what §4.1 measures the cost of, and the one weighted by
/// each position's mismapped posterior, which is §6's subject. The library computes only the third.
///
/// **Without this check the copy would drift silently.** Every figure this program prints would go
/// on describing arithmetic that no longer ran anywhere, which is exactly what the implementation
/// plan's *"a change that moves these numbers and leaves the tests green is a change whose effect
/// nobody has looked at"* is about.
///
/// The tolerance is relative and `1e-9`: both routes read the same `f32` array into `f64` and do
/// the same arithmetic in the same order, so what separates them is summation order alone.
///
/// **Run at `F = 0`, because the copy applies no inbreeding correction** — the correction is a
/// constant divide the harness's own arms do not want, so comparing at any other coefficient would
/// be comparing two different quantities.
fn check_against_the_shipped_estimator(
    genotype_posterior: &[f32],
    individuals: usize,
    positions: usize,
    ours: &MomentsFromPosteriors,
) {
    let outbred = InbreedingF::try_new(0.0).expect("a legal coefficient");
    let shipped =
        CensusMoments::from_posteriors(genotype_posterior, individuals, positions, outbred);
    let apart = |ours: f64, shipped: f64| {
        if ours == 0.0 && shipped == 0.0 {
            0.0
        } else {
            (shipped / ours - 1.0).abs()
        }
    };
    let frequency = apart(ours.frequency, shipped.mean_alternative_frequency);
    let heterozygosity = apart(ours.heterozygosity_with_variance, shipped.heterozygosity);
    assert!(
        frequency < 1e-9 && heterozygosity < 1e-9,
        "this program's own reduction and the shipped `CensusMoments` disagree — frequency by \
         {frequency:e}, heterozygosity by {heterozygosity:e}. Every number below is then about \
         code that no longer runs; fix the copy or the library before reading any of it."
    );
}

/// Turn a mean frequency and a heterozygosity into the concentration pair the genotype prior
/// takes.
///
/// **The same arithmetic the shipped seam uses once it has a shape** (`ordinary_site_seed.md` §3):
/// a pair of expected frequency `f` and total `A` implies a heterozygosity of
/// `2 f (1 - f) A / (A + 1)`, so the total that reproduces a measured `pi` is
/// `A = pi / (2 f (1 - f) - pi)`. A measurement no pair can reach — `pi` at or above `2f(1-f)` —
/// falls back to the neutral pair, which is what the shipped seam does too.
fn seed_from_moments(frequency: f64, heterozygosity: f64) -> SpectrumSeed {
    use pop_var_caller::ng::calling::genotype_prior::SeedRegime;
    let frequency = frequency.clamp(1e-12, 1.0 - 1e-12);
    let ceiling = 2.0 * frequency * (1.0 - frequency);
    if heterozygosity <= 0.0 || heterozygosity >= ceiling {
        return SpectrumSeed::new(
            1.0,
            heterozygosity.clamp(1e-12, 0.5),
            SeedRegime::NeutralShape,
        );
    }
    let total = heterozygosity / (ceiling - heterozygosity);
    SpectrumSeed::new(
        total * (1.0 - frequency),
        total * frequency,
        SeedRegime::NeutralShape,
    )
}

// ---------------------------------------------------------------------------------------------
// Does the difference move a genotype?
// ---------------------------------------------------------------------------------------------

/// **The share of single-sample genotype calls that come out differently under the two seeds.**
///
/// One locus is one drawn position with two alleles: the reference and one alternative. The
/// sample's genotype is drawn from the population at that position, its reads are drawn at the
/// depth, and the call is whichever of the three genotypes has the largest posterior — the
/// caller's own prior (`MarginalizedDirichletPrior`) times the reads' likelihood.
///
/// **This is a genotype call and not the caller**: it has no candidate-allele selection, no
/// quality gate and no cohort step. What it isolates is the one thing the two seeds differ in,
/// which is what the question asks.
#[allow(
    clippy::too_many_arguments,
    reason = "two seeds, the population they are compared on, and the cell they belong to"
)]
fn genotype_calls(
    todays: SpectrumSeed,
    direct: SpectrumSeed,
    inbreeding: InbreedingF,
    population: &Population,
    depth: f64,
    drawn_inbreeding: f64,
    loci: usize,
    replicate: usize,
    individuals: usize,
) -> GenotypeComparison {
    let todays_prior = genotype_log_priors(todays, inbreeding);
    let direct_prior = genotype_log_priors(direct, inbreeding);
    // **The control.** The same comparison against a seed whose alternative concentration is three
    // times the direct one — far outside anything the two routes disagree by. If this column is
    // zero too, the comparison cannot see a prior at all and the columns beside it say nothing.
    let trebled_prior = genotype_log_priors(
        SpectrumSeed::new(
            direct.alpha_ref(),
            direct.alpha_alt_total() * 3.0,
            direct.regime(),
        ),
        inbreeding,
    );
    // **Keyed on the replicate and the panel size as well as the depth, and that is a defect a
    // review caught.** Seeded by depth alone, every replicate and every panel size at one depth
    // called the identical 200,000 loci with the identical reads, so the spread printed beside
    // `calls moved` was a spread over seed pairs on one fixed locus set rather than over
    // independent draws.
    let mut rng = Rng(mix(
        0x2545_F491_4F6C_DD1D,
        (depth as u64).wrapping_mul(0x9E37_79B9) ^ replicate as u64,
        individuals as u64,
    ));
    let mut result = GenotypeComparison {
        segregating_loci: 0,
        moved_at_segregating_loci: 0,
        moved_when_trebled: 0,
        wrong_under_todays_path: 0,
        wrong_under_the_direct_moments: 0,
    };
    for _ in 0..loci {
        let frequency = population.draw_frequency(&mut rng);
        let genotype = draw_genotype(&mut rng, frequency, drawn_inbreeding);
        let reads = rng.poisson(depth).min(CENSUS_DEPTH_CAP);
        let mut alternative_reads = 0_u32;
        let carried = f64::from(genotype) / 2.0;
        let on_alternative =
            carried * (1.0 - CLEAN_ERROR_RATE) + (1.0 - carried) * CLEAN_ERROR_RATE / 3.0;
        for _ in 0..reads {
            if rng.uniform() < on_alternative {
                alternative_reads += 1;
            }
        }
        let likelihood = read_log_likelihoods(alternative_reads, reads);
        let under_todays = call(&todays_prior, &likelihood);
        let under_direct = call(&direct_prior, &likelihood);
        let under_trebled = call(&trebled_prior, &likelihood);
        // **Counted over the positions that segregate, and not over every position.** A
        // population that is invariant at 98 positions in 100 makes any share over all of them a
        // statement about how rare variants are, which nobody needs from this program.
        if frequency <= 0.0 || frequency >= 1.0 {
            continue;
        }
        result.segregating_loci += 1;
        if under_todays != under_direct {
            result.moved_at_segregating_loci += 1;
        }
        if under_todays != under_trebled {
            result.moved_when_trebled += 1;
        }
        if under_todays != genotype as usize {
            result.wrong_under_todays_path += 1;
        }
        if under_direct != genotype as usize {
            result.wrong_under_the_direct_moments += 1;
        }
    }
    result
}

/// What the two seeds do to a genotype, counted over the loci that segregate.
struct GenotypeComparison {
    segregating_loci: u64,
    moved_at_segregating_loci: u64,
    /// The control: how many move when the direct seed's alternative concentration is trebled.
    moved_when_trebled: u64,
    wrong_under_todays_path: u64,
    wrong_under_the_direct_moments: u64,
}

impl GenotypeComparison {
    fn moved_in_a_hundred(&self) -> f64 {
        100.0 * self.moved_at_segregating_loci as f64 / self.segregating_loci.max(1) as f64
    }

    fn moved_when_trebled_in_a_hundred(&self) -> f64 {
        100.0 * self.moved_when_trebled as f64 / self.segregating_loci.max(1) as f64
    }

    /// How often the call disagrees with the genotype the locus was drawn with. **Kept because a
    /// run where both seeds are equally wrong and a run where both are right look identical in
    /// the `moved` column**, and only this tells them apart.
    #[allow(
        dead_code,
        reason = "read by the assertion below and by a reader of the source"
    )]
    fn wrong_in_a_hundred(&self, wrong: u64) -> f64 {
        100.0 * wrong as f64 / self.segregating_loci.max(1) as f64
    }
}

/// The caller's own genotype prior at a biallelic diploid locus, as three log probabilities in the
/// order homozygous reference, heterozygous, homozygous alternative.
fn genotype_log_priors(seed: SpectrumSeed, inbreeding: InbreedingF) -> [f64; 3] {
    let mut concentration = [0.0_f64; 2];
    let concentration =
        fill_locus_concentration(seed, VariantClass::Substitution, 2, &mut concentration);
    // Diploid, two alleles: two copies of the reference, one of each, two of the alternative.
    let genotype_allele_counts: [u32; 6] = [2, 0, 1, 1, 0, 2];
    let log_multinomial_coefficients = [0.0, std::f64::consts::LN_2, 0.0];
    let homozygous_allele_for = [Some(AlleleId(0)), None, Some(AlleleId(1))];
    let mut scratch = [0.0_f64; 2];
    let mut out = [LogProb(0.0); 3];
    let mut row = PriorRow::new(
        concentration,
        &genotype_allele_counts,
        &log_multinomial_coefficients,
        &homozygous_allele_for,
        &mut scratch,
        &mut out,
    );
    MarginalizedDirichletPrior.fill_genotype_log_priors(&mut row, inbreeding);
    [out[0].get(), out[1].get(), out[2].get()]
}

/// `ln P(reads | genotype)` for each of the three genotypes, dropping the binomial coefficient the
/// three share.
fn read_log_likelihoods(alternative_reads: u32, reads: u32) -> [f64; 3] {
    let mut out = [0.0_f64; 3];
    for (copies, slot) in out.iter_mut().enumerate() {
        let carried = copies as f64 / 2.0;
        let on_alternative = (carried * (1.0 - CLEAN_ERROR_RATE)
            + (1.0 - carried) * CLEAN_ERROR_RATE / 3.0)
            .clamp(1e-12, 1.0 - 1e-12);
        *slot = f64::from(alternative_reads) * on_alternative.ln()
            + f64::from(reads - alternative_reads) * (1.0 - on_alternative).ln();
    }
    out
}

/// Which genotype has the largest posterior.
fn call(prior: &[f64; 3], likelihood: &[f64; 3]) -> usize {
    let mut best = 0;
    let mut best_score = f64::NEG_INFINITY;
    for index in 0..3 {
        let score = prior[index] + likelihood[index];
        if score > best_score {
            best_score = score;
            best = index;
        }
    }
    best
}

// ---------------------------------------------------------------------------------------------
// Drawing a cohort into the records the fit reads
// ---------------------------------------------------------------------------------------------

struct Drawn {
    samples: Vec<SampleCensusEvidence>,
    /// Each sample's alternative-copy count at each position — the truth the oracle arm reads.
    genotypes: Vec<Vec<u8>>,
    drawn_mean_frequency: f64,
    drawn_heterozygosity: f64,

    /// Which positions were planted mismapped. **Kept so the genotype-derived oracle can be
    /// restricted to the same positions the estimate is asked about** — scoring a panel estimator
    /// against a population quantity leaves the inbreeding factor uncancelled, which at one sample
    /// is 15% and swamps what this half measures.
    planted_mismapped: Vec<bool>,
    /// How many of the drawn positions carry the alternative allele at a frequency strictly
    /// between zero and one — the count everything about the population rests on.
    segregating_positions: u64,
}

#[allow(
    clippy::too_many_arguments,
    reason = "the drawn cohort's own parameters"
)]
fn draw(
    population: &Population,
    samples: usize,
    positions: usize,
    depth: f64,
    inbreeding: f64,
    mismapped_share: f64,
    mismapped_error_rate: f64,
    seed: u64,
) -> Drawn {
    let mut rng = Rng(seed);
    let edges = DepthBinEdges::for_census();
    let mut codes: Vec<PackedDepthCodes> = (0..samples)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
    let mut genotypes: Vec<Vec<u8>> = Vec::with_capacity(positions);
    let mut frequency_total = 0.0_f64;
    let mut heterozygosity_total = 0.0_f64;
    let mut segregating_positions = 0_u64;
    let mut planted_mismapped: Vec<bool> = Vec::with_capacity(positions);

    for index in 0..positions {
        let mismapped = rng.uniform() < mismapped_share;
        planted_mismapped.push(mismapped);
        let rate = if mismapped {
            mismapped_error_rate
        } else {
            CLEAN_ERROR_RATE
        };
        let frequency = population.draw_frequency(&mut rng);
        frequency_total += frequency;
        heterozygosity_total += 2.0 * frequency * (1.0 - frequency);
        if frequency > 0.0 && frequency < 1.0 {
            segregating_positions += 1;
        }
        // **Code 0 is the reference base by construction**, so the three candidates are codes 1
        // to 3 and the fit's own sum over which one segregates has to find the right one.
        let allele = 1 + ((rng.uniform() * 3.0) as usize).min(2);
        let mut at_position = vec![0_u8; samples];
        for sample in 0..samples {
            let reads = rng.poisson(depth).min(CENSUS_DEPTH_CAP);
            let genotype = draw_genotype(&mut rng, frequency, inbreeding);
            at_position[sample] = genotype as u8;
            let carried = f64::from(genotype) / 2.0;
            let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
            let on_reference = (1.0 - carried) * (1.0 - rate) + carried * rate / 3.0;
            let mut counts = [0_u32; 5];
            for _ in 0..reads {
                let u = rng.uniform();
                let code = if u < on_candidate {
                    allele
                } else if u < on_candidate + on_reference {
                    0
                } else if u < on_candidate + on_reference + rate / 3.0 {
                    (allele % 3) + 1
                } else {
                    ((allele + 1) % 3) + 1
                };
                counts[code] += 1;
            }
            codes[sample].set(index, DepthCode::Binned(edges.bin_for(reads)));
            for (code, count) in counts.iter().enumerate() {
                if code == 0 || *count == 0 {
                    continue;
                }
                sparse[sample].push(AlleleObservation {
                    index: index as u32,
                    allele: match code {
                        1 => ObservedAllele::C,
                        2 => ObservedAllele::G,
                        _ => ObservedAllele::T,
                    },
                    reads: u8::try_from((*count).min(255))
                        .expect("a clamped count fits the census's one-byte field"),
                });
            }
        }
        genotypes.push(at_position);
    }

    let terms = RecordingTerms {
        selection: SelectionTermsDigest::of(&SelectionTerms {
            seed,
            reference: ReferenceDigest([7; 16]),
            analysed_regions: RegionSetDigest([9; 16]),
            catalog_built_under: CatalogBuildSettings {
                criteria: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "0.1.0".to_string(),
            },
            ssr_criteria: StrRepeatCriteria::default(),
            generic_target: positions as u64,
            ssr_cap: 1_000,
        }),
        kept_loci: CensusLociDigester::new().finish(),
        ssr_stratum_counts: Default::default(),
        read_cap: ReadCap(1_000),
        depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
        depth_cap: DepthCap::new(CENSUS_DEPTH_CAP),
    };
    let records = (0..samples)
        .map(|s| {
            SampleCensusEvidence::resident(
                format!("s{s:03}"),
                terms.clone(),
                NamedReadGroup::drawn_for(&format!("s{s:03}"), [ReadGroupId(s as u32)]),
                BTreeMap::new(),
                BTreeMap::from([(
                    SectionKey::Generic(ReadGroupId(s as u32)),
                    Section::Generic(GenericEvidence::from_parts(
                        std::mem::replace(&mut codes[s], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[s]),
                    )),
                )]),
            )
        })
        .collect();

    Drawn {
        samples: records,
        genotypes,
        drawn_mean_frequency: frequency_total / positions as f64,
        drawn_heterozygosity: heterozygosity_total / positions as f64,
        segregating_positions,
        planted_mismapped,
    }
}

/// Draw one diploid individual's count of alternative copies at a position of frequency `f`: with
/// probability `F` the two copies are one ancestral copy counted twice, otherwise two independent
/// draws.
fn draw_genotype(rng: &mut Rng, frequency: f64, inbreeding: f64) -> u32 {
    if rng.uniform() < inbreeding {
        if rng.uniform() < frequency { 2 } else { 0 }
    } else {
        u32::from(rng.uniform() < frequency) + u32::from(rng.uniform() < frequency)
    }
}

/// The stream every drawn number comes from — the same xorshift the sibling harnesses use.
struct Rng(u64);

impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn pick(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut cut = self.uniform() * total;
        for (index, weight) in weights.iter().enumerate() {
            cut -= weight;
            if cut <= 0.0 {
                return index;
            }
        }
        weights.len() - 1
    }

    fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        (x / (x + y)).clamp(1e-9, 1.0 - 1e-9)
    }

    fn gamma(&mut self, shape: f64) -> f64 {
        if shape < 1.0 {
            let u = self.uniform().max(1e-12);
            return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let z = self.normal();
            let v = (1.0 + c * z).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.uniform().max(1e-12);
            if u.ln() < 0.5 * z * z + d - d * v + d * v.ln() {
                return d * v;
            }
        }
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-12);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    fn poisson(&mut self, mean: f64) -> u32 {
        let limit = (-mean).exp();
        let mut product = self.uniform();
        let mut count = 0;
        while product > limit && count < 400 {
            count += 1;
            product *= self.uniform();
        }
        count
    }
}
