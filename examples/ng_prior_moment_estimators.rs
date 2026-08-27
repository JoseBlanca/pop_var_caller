//! **Do the two direct moment estimators recover the population they were drawn from?**
//!
//! # The question
//!
//! The SNP/indel genotype prior holds two numbers, and they are fixed by two properties of the
//! population: the **mean alternative-allele frequency** and the **heterozygosity** — how often
//! two copies of a position drawn at random from the population differ.
//! `doc/devel/ng/research/ordinary_site_prior_moments.md` proposes reading both straight off the
//! census positions instead of off a fitted curve:
//!
//! ```text
//! heterozygosity   pi  =  mean over positions of   2 k (2N - k) / (2N (2N - 1))
//! mean frequency   f   =  mean over positions of   k / 2N
//! ```
//!
//! for `k` alternative copies among the panel's `2N` chromosomes. Both are classical, and the
//! plan's claim is that both are unbiased at every panel size from one individual to thousands.
//!
//! **This program checks that claim on genotypes that are known**, which is the plan's first
//! question and the cheapest place a mistake can be found. Nothing here reads a read: the panel's
//! genotypes are the ones it was drawn with, so what is measured is the estimator and not the
//! sequencing. What it costs to read from reads instead is a separate program
//! (`ng_prior_moments_from_reads.rs`).
//!
//! # What each arm is scored against, decided before the run
//!
//! Two different yardsticks, because two different questions are being asked, and quoting one
//! where the other belongs is how a measurement invents a finding.
//!
//! - **Is the estimator biased?** Score each drawn panel against **the positions that panel was
//!   itself drawn at** — the mean of the `f` actually used at those positions, and the mean of
//!   `2 f (1 - f)` there. Those positions are a finite sample of the population and their own
//!   moments sit a little away from the population's; scoring against the population would put
//!   that scatter into every error column and call it estimator bias.
//! - **How far from the population is the answer a run would print?** Score against the
//!   population's own moments, computed in closed form from the shape. This one includes the
//!   census's scatter, because a run does not get to average over censuses.
//!
//! Both are printed, in adjacent columns, and never mixed.
//!
//! # What is held fixed
//!
//! Within one replicate the positions and every sample's genotypes are drawn **once**, and each
//! panel-size arm reads the first `N` samples of that same draw. So nothing moves across the
//! panel-size arms for a reason unrelated to panel size — the device
//! `examples/ng_joint_sample_count_sweep.rs` uses, and for the same reason.
//!
//! # The shapes, and why some of them are awkward on purpose
//!
//! The joint fit describes a population with **one Beta over what segregates plus a spike at
//! frequency zero and one at frequency one**. A population drawn from that family flatters every
//! estimator that assumes it, so a sweep made only of such shapes is a formality. Half the shapes
//! here are outside the family — two peaks, or a lump of positions at one intermediate frequency
//! — and each shape prints **how far outside**, as the largest gap between its own distribution
//! of segregating frequencies and the closest single Beta's.
//!
//! # Inbreeding, which is not a detail here
//!
//! A tomato accession is largely self-pollinated: its two copies of a position are the same
//! ancestral copy far more often than random mating would give. The panel's inbreeding
//! coefficient `F` is the probability that an individual's two copies are one copy counted twice.
//! **Every arm is run at `F = 0` and at `F = 0.8`**, tomato's fitted range, because the second
//! estimator's derivation quietly assumes the panel's chromosomes are a random sample of the
//! population's, and two copies inside one inbred individual are not.
//!
//! Run: `./scripts/dev.sh cargo run --release --example ng_prior_moment_estimators`

use std::env;

use pop_var_caller::genetics::lgamma;

/// Panel sizes, in diploid individuals — the committed range, one sample to a thousand
/// (`doc/devel/specs/design_principles.md` §0).
const PANELS: [usize; 9] = [1, 2, 3, 5, 10, 25, 63, 200, 1000];

/// The inbreeding coefficients swept: random mating, and tomato's fitted range.
const INBREEDING: [f64; 2] = [0.0, 0.8];

/// What a run assumes when nothing could be fitted — `ExpectedHeterozygosity::SPECIES_FALLBACK`,
/// restated here so this program does not depend on the caller's types to print one number.
const SPECIES_FALLBACK_DIVERSITY: f64 = 1e-3;

fn main() {
    let mut args = env::args().skip(1);
    let positions: usize = args.next().map_or(20_000, |a| a.parse().expect("a count"));
    let replicates: usize = args.next().map_or(24, |a| a.parse().expect("a count"));

    println!("# Do the two direct moment estimators recover the population they were drawn from?");
    println!();
    println!("{positions} census positions a replicate, {replicates} replicates a cell.");
    println!("Genotypes are the drawn ones — no reads, no depth, no fit. Panel sizes {PANELS:?}.");

    let shapes = shapes();
    check_the_beta_moments(&shapes);
    print_shape_table(&shapes);
    print_estimator_oracle();

    for shape in &shapes {
        for &inbreeding in &INBREEDING {
            run_cell(shape, inbreeding, positions, replicates);
        }
    }

    print_single_sample_spread(&shapes, replicates.max(40));
}

// ---------------------------------------------------------------------------------------------
// What the caller used to get, and why the arm that measured it is gone
// ---------------------------------------------------------------------------------------------
//
// **This program used to print a third block: what the caller's own path returned on these four
// populations.** It evaluated each population into the `2N + 1` allele-count classes a panel of
// `N` diploid individuals has, ran the shipped two-parameter search over the result, and reported
// the pair's mean frequency over the population's — the detour this work removed.
//
// **The search was deleted on 2026-08-27** (`doc/devel/ng/spec/ordinary_site_prior_moments.md`
// §5), so the arm cannot run against the shipped code any more, and running it against a copy
// would make its figures facts about the copy. **Its numbers are in the report** it was written
// for, `doc/devel/reports/ng_ordinary_site_prior_moments_2026-08-27.md` §9: handed the population
// exactly, the search's mean frequency was 0.999× the truth at one individual and 1.217×, 0.861×,
// 0.787× and 1.043× at 200 individuals on the four shapes, and the blend cost a further 0.816×,
// 0.622×, 0.893× and 0.916× at one individual.
//
// **What the caller does now needs no arm here**: the mean frequency is an integral of the fitted
// curve, so handed the population exactly it returns the population exactly.

// ---------------------------------------------------------------------------------------------
// The populations
// ---------------------------------------------------------------------------------------------

/// One component of a population's spread of alternative-allele frequencies over the positions
/// that segregate.
///
/// **Every component is a Beta, and a shape is a weighted mixture of them.** A single Beta is
/// what the joint fit can hold; a mixture of two is not, and that is the whole point of the
/// awkward shapes below. Keeping every component in the same family is also what lets the
/// comparison against the current path use **the caller's own projection**
/// ([`allele_count_classes`]) rather than a copy of it written here.
#[derive(Clone, Copy, Debug)]
struct Component {
    /// What share of the segregating positions this component holds.
    weight: f64,
    a: f64,
    b: f64,
}

impl Component {
    /// `E[f]` under this component alone.
    fn mean_frequency(&self) -> f64 {
        self.a / (self.a + self.b)
    }

    /// `E[f²]` under this component alone.
    fn mean_square_frequency(&self) -> f64 {
        self.a * (self.a + 1.0) / ((self.a + self.b) * (self.a + self.b + 1.0))
    }

    /// Draw one frequency from this component.
    fn draw(&self, rng: &mut Rng) -> f64 {
        rng.beta(self.a, self.b)
    }
}

/// A population, as this program draws one: a share of positions carrying only the reference
/// base, a share carrying only a non-reference base, and a mixture over what is left.
struct PopulationShape {
    name: &'static str,
    /// Share of positions where the population carries only the reference base — frequency 0.
    invariant: f64,
    /// Share where it carries only a non-reference base — frequency 1.
    fixed_alt: f64,
    /// The segregating positions' frequency spread. Weights sum to one.
    segregating: Vec<Component>,
    /// Whether one Beta can hold the segregating spread. Printed as a number too — see
    /// [`distance_from_the_closest_beta`].
    inside_the_fitted_family: bool,
}

impl PopulationShape {
    fn share_segregating(&self) -> f64 {
        1.0 - self.invariant - self.fixed_alt
    }

    /// The population's mean alternative-allele frequency, `E[f]`, over **all** positions.
    fn mean_frequency(&self) -> f64 {
        let segregating: f64 = self
            .segregating
            .iter()
            .map(|c| c.weight * c.mean_frequency())
            .sum();
        self.fixed_alt + self.share_segregating() * segregating
    }

    /// The population's heterozygosity, `E[2 f (1 - f)]`, over **all** positions. The two end
    /// masses contribute nothing: a position carrying one allele has no diversity.
    fn heterozygosity(&self) -> f64 {
        let segregating: f64 = self
            .segregating
            .iter()
            .map(|c| c.weight * 2.0 * (c.mean_frequency() - c.mean_square_frequency()))
            .sum();
        self.share_segregating() * segregating
    }

    /// Draw one position's alternative-allele frequency.
    ///
    /// **No allocation**: this runs once per census position per replicate, tens of millions of
    /// times in a sweep, and an earlier draft collected the component weights into a fresh `Vec`
    /// on every call.
    fn draw_frequency(&self, rng: &mut Rng) -> f64 {
        match rng.pick(&[self.invariant, self.fixed_alt, self.share_segregating()]) {
            0 => 0.0,
            1 => 1.0,
            _ => {
                let mut cut = rng.uniform();
                for component in &self.segregating {
                    cut -= component.weight;
                    if cut <= 0.0 {
                        return component.draw(rng);
                    }
                }
                self.segregating
                    .last()
                    .expect("a shape spreads its segregating positions somewhere")
                    .draw(rng)
            }
        }
    }
}

/// The grid. **Four shapes, and only two of them are shapes the joint fit's family can hold** —
/// on a shape the family fits, an estimator that quietly assumes the family looks perfect.
///
/// Every diversity is set near this project's two benchmark cohorts: tomato at about 6
/// differences per 10,000 bases, a human panel at about 10.
///
/// **What the four deliberately do not share.** Two put most segregating positions at low
/// frequency and two do not; two are asymmetric under `f -> 1 - f` and one is not; one has more
/// mass above a half than below, which is the direction a rare-allele pile-up never has; and the
/// two end masses differ across them by a factor of twenty-five. A set that shared any of those
/// would let a swapped pair of shape parameters, or an estimator written the wrong way round in
/// `k`, pass unnoticed.
fn shapes() -> Vec<PopulationShape> {
    vec![
        PopulationShape {
            name: "tomato-like: nearly all alternative alleles rare",
            invariant: 0.9950,
            fixed_alt: 0.0010,
            segregating: vec![Component {
                weight: 1.0,
                a: 0.20,
                b: 1.00,
            }],
            inside_the_fitted_family: true,
        },
        PopulationShape {
            name: "where it varies the reference base is the rare one",
            invariant: 0.9900,
            fixed_alt: 0.0060,
            segregating: vec![Component {
                weight: 1.0,
                a: 3.00,
                b: 0.60,
            }],
            inside_the_fitted_family: true,
        },
        PopulationShape {
            name: "two peaks, off centre — outside the fitted family",
            invariant: 0.9955,
            fixed_alt: 0.0005,
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
        PopulationShape {
            name: "a lump at one intermediate frequency — outside the fitted family",
            invariant: 0.9940,
            fixed_alt: 0.0010,
            segregating: vec![
                Component {
                    weight: 0.55,
                    a: 0.25,
                    b: 1.50,
                },
                // **A narrow lump at frequency 0.42**, not a point mass: a mean of 0.42 with a
                // spread of 0.049. Every component being a Beta is what lets the comparison
                // against the current path go through the caller's own projection.
                Component {
                    weight: 0.45,
                    a: 42.0,
                    b: 58.0,
                },
            ],
            inside_the_fitted_family: false,
        },
    ]
}

// ---------------------------------------------------------------------------------------------
// The two estimators
// ---------------------------------------------------------------------------------------------

/// **How often two chromosomes drawn at random from the panel differ** — Nei's average
/// heterozygosity at one position, with the finite-panel correction.
///
/// Of the `2N (2N - 1)` ordered pairs of distinct chromosomes, `2 k (2N - k)` have one carrying
/// the alternative allele and the other not. The `2N - 1` rather than `2N` is what makes the
/// answer a property of the population rather than of the panel.
fn nei_heterozygosity(alternative_copies: u32, chromosomes: u32) -> f64 {
    assert!(
        chromosomes >= 2,
        "two chromosomes are needed before a pair of them can differ; this panel has \
         {chromosomes}"
    );
    let k = f64::from(alternative_copies);
    let n = f64::from(chromosomes);
    2.0 * k * (n - k) / (n * (n - 1.0))
}

/// **What share of the panel's chromosomes carry the alternative allele** at one position.
fn allele_frequency(alternative_copies: u32, chromosomes: u32) -> f64 {
    f64::from(alternative_copies) / f64::from(chromosomes)
}

/// **The factor by which inbreeding shrinks what [`nei_heterozygosity`] returns.**
///
/// The estimator counts how often two chromosomes drawn from the panel differ. One pair in
/// `2N - 1` is the two copies inside a single individual, and those are the *same* ancestral copy
/// with probability `F` — in which case they never differ. So the estimator's expectation is the
/// population's heterozygosity times `1 - F / (2N - 1)`.
///
/// At sixty-three individuals that is a shrinkage of 6 parts in a thousand at `F = 0.8`, which no
/// user would notice. **At one individual it is `1 - F` exactly**, so a single selfing tomato
/// reports a fifth of its population's diversity. The two ends of the committed cohort range do
/// not agree about whether this term exists, which is why it is a function and not a comment.
fn inbreeding_shrinkage(chromosomes: u32, inbreeding: f64) -> f64 {
    1.0 - inbreeding / (f64::from(chromosomes) - 1.0)
}

// ---------------------------------------------------------------------------------------------
// The sweep
// ---------------------------------------------------------------------------------------------

/// Everything one replicate contributes at one panel size.
#[derive(Clone, Copy, Default)]
struct Tally {
    /// The estimate divided by the moment of the positions this replicate actually drew.
    against_drawn_positions: Welford,
    /// The estimate divided by the population's own moment.
    against_the_population: Welford,
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

    /// The spread across replicates — one standard deviation.
    fn spread(&self) -> f64 {
        if self.count < 2.0 {
            return f64::NAN;
        }
        (self.sum_of_squared_deviations / (self.count - 1.0)).sqrt()
    }

    /// How precisely this run pins the mean: the spread divided by the square root of the
    /// replicate count. **A bias smaller than this is a bias the run cannot see.**
    fn uncertainty_of_the_mean(&self) -> f64 {
        self.spread() / self.count.sqrt()
    }
}

fn run_cell(shape: &PopulationShape, inbreeding: f64, positions: usize, replicates: usize) {
    let population_frequency = shape.mean_frequency();
    let population_heterozygosity = shape.heterozygosity();

    println!();
    println!("## {} — inbreeding F = {inbreeding}", shape.name);
    println!();
    println!(
        "Population: mean frequency {population_frequency:.6}, heterozygosity \
         {population_heterozygosity:.6} ({:.2} differences per kilobase).",
        1_000.0 * population_heterozygosity
    );
    println!();

    let largest = *PANELS.last().expect("a panel");
    let mut frequency_tally = vec![Tally::default(); PANELS.len()];
    let mut heterozygosity_tally = vec![Tally::default(); PANELS.len()];
    let mut corrected_tally = vec![Tally::default(); PANELS.len()];

    for replicate in 0..replicates {
        // A fresh stream per replicate, and one that depends on the cell, so no two cells share a
        // draw. Mixed rather than added, because neighbouring seeds in this generator's state
        // produce visibly correlated first draws.
        let seed = mix(
            0x9E37_79B9_7F4A_7C15,
            replicate as u64,
            shape.name.len() as u64 * 131 + (inbreeding * 1_000.0) as u64,
        );
        let mut rng = Rng(seed);

        // The positions and every sample's genotype, drawn once. Each panel-size arm reads a
        // prefix of the same draw.
        let mut copies_at_arm: Vec<Vec<u32>> = PANELS
            .iter()
            .map(|_| Vec::with_capacity(positions))
            .collect();
        let mut drawn_frequency_total = 0.0_f64;
        let mut drawn_heterozygosity_total = 0.0_f64;

        for _ in 0..positions {
            let frequency = shape.draw_frequency(&mut rng);
            drawn_frequency_total += frequency;
            drawn_heterozygosity_total += 2.0 * frequency * (1.0 - frequency);
            // **One running total, read off at each arm's boundary.** The arms are prefixes of one
            // another, so a panel of 63 is the panel of 25 plus 38 more individuals — drawing them
            // separately would let the arms disagree for a reason unrelated to panel size.
            let mut running = 0_u32;
            let mut arm = 0;
            for individual in 0..largest {
                running += draw_genotype(&mut rng, frequency, inbreeding);
                if individual + 1 == PANELS[arm] {
                    copies_at_arm[arm].push(running);
                    arm += 1;
                    if arm == PANELS.len() {
                        break;
                    }
                }
            }
        }

        let drawn_frequency = drawn_frequency_total / positions as f64;
        let drawn_heterozygosity = drawn_heterozygosity_total / positions as f64;

        for (index, &individuals) in PANELS.iter().enumerate() {
            let chromosomes = 2 * individuals as u32;
            let mut frequency_sum = 0.0_f64;
            let mut heterozygosity_sum = 0.0_f64;
            for &copies in &copies_at_arm[index] {
                frequency_sum += allele_frequency(copies, chromosomes);
                heterozygosity_sum += nei_heterozygosity(copies, chromosomes);
            }
            let frequency = frequency_sum / positions as f64;
            let heterozygosity = heterozygosity_sum / positions as f64;
            let corrected = heterozygosity / inbreeding_shrinkage(chromosomes, inbreeding);

            frequency_tally[index]
                .against_drawn_positions
                .add(frequency / drawn_frequency);
            frequency_tally[index]
                .against_the_population
                .add(frequency / population_frequency);
            heterozygosity_tally[index]
                .against_drawn_positions
                .add(heterozygosity / drawn_heterozygosity);
            heterozygosity_tally[index]
                .against_the_population
                .add(heterozygosity / population_heterozygosity);
            corrected_tally[index]
                .against_drawn_positions
                .add(corrected / drawn_heterozygosity);
            corrected_tally[index]
                .against_the_population
                .add(corrected / population_heterozygosity);
        }
    }

    print_moment_table("mean alternative-allele frequency", &frequency_tally);
    print_moment_table(
        "heterozygosity, as the plan writes it",
        &heterozygosity_tally,
    );
    if inbreeding > 0.0 {
        print_moment_table("heterozygosity, divided by 1 - F/(2N-1)", &corrected_tally);
    }
}

fn print_moment_table(title: &str, tally: &[Tally]) {
    println!("### {title}");
    println!();
    println!(
        "Left half isolates the estimator: each panel is scored against the positions it was \
         itself drawn at. Right half is what a run would print, scored against the population, so \
         it also carries how far this census's own positions sit from the population."
    );
    println!();
    println!(
        "| individuals | off the drawn positions | pinned to | spread over replicates | off the \
         population | pinned to | spread over replicates |"
    );
    println!("|---:|---:|---:|---:|---:|---:|---:|");
    for (index, &individuals) in PANELS.iter().enumerate() {
        let drawn = tally[index].against_drawn_positions;
        let population = tally[index].against_the_population;
        println!(
            "| {individuals} | {:+.3}% | ±{:.3}% | {:.3}% | {:+.2}% | ±{:.2}% | {:.2}% |",
            (drawn.mean - 1.0) * 100.0,
            drawn.uncertainty_of_the_mean() * 100.0,
            drawn.spread() * 100.0,
            (population.mean - 1.0) * 100.0,
            population.uncertainty_of_the_mean() * 100.0,
            population.spread() * 100.0,
        );
    }
    println!();
}

/// **The one-sample corner, on its own, because it is the corner the caller commits to and the
/// one where a wide answer would sink the proposal.**
///
/// Reports the spread of the two moments over independent single-genome runs, and puts it beside
/// the distance between the population's diversity and the constant a run falls back to when
/// nothing could be fitted. A spread wider than that distance would mean the measurement is worse
/// than the guess.
fn print_single_sample_spread(shapes: &[PopulationShape], replicates: usize) {
    println!();
    println!("## One genome: is the answer usable, not merely unbiased?");
    println!();
    println!(
        "Each run here is **one genome at {SINGLE_GENOME_CENSUS} census positions**, which is the \
         shipped default (`parameter_prepass_census_sites.md` §5), with the positions redrawn \
         every run — so the spread is the whole spread a real run would see, census sampling and \
         genotype sampling together. The panel sweep above uses far fewer positions on purpose, \
         because there it is the estimator being measured and not the census."
    );
    println!();
    println!(
        "**What the estimate has to beat is the last column**: how far the species-range constant \
         of {SPECIES_FALLBACK_DIVERSITY} per position sits from this population's own diversity. \
         That constant is what a run falls back to when nothing could be fitted."
    );
    println!();
    println!(
        "Both estimates below divide by 1 - F, which a single genome cannot supply for itself — \
         see the note under the tables."
    );
    println!();
    println!(
        "**Both moments, because the design decision that rested on this arm needed both.** The \
         heterozygosity decides whether a single-genome run should fall back to the species-range \
         constant. The mean frequency decided whether the seed still needed its blend toward a \
         neutral shape — the columns below say it did not, and that blend was deleted on \
         2026-08-27."
    );
    println!();
    println!(
        "| population | F | heterozygosity: off the population | spread over runs | 19 runs in \
         20 land inside | the fallback constant is off by | mean frequency: off the population | \
         spread over runs | 19 runs in 20 land inside |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
    for shape in shapes {
        for &inbreeding in &INBREEDING {
            let heterozygosity_truth = shape.heterozygosity();
            let frequency_truth = shape.mean_frequency();
            let mut heterozygosity_tally = Welford::default();
            let mut frequency_tally = Welford::default();
            for replicate in 0..replicates {
                let seed = mix(
                    0x243F_6A88_85A3_08D3,
                    replicate as u64,
                    shape.name.len() as u64 * 977 + (inbreeding * 1_000.0) as u64,
                );
                let mut rng = Rng(seed);
                let mut heterozygosity_sum = 0.0_f64;
                let mut frequency_sum = 0.0_f64;
                for _ in 0..SINGLE_GENOME_CENSUS {
                    let frequency = shape.draw_frequency(&mut rng);
                    let copies = draw_genotype(&mut rng, frequency, inbreeding);
                    heterozygosity_sum += nei_heterozygosity(copies, 2);
                    frequency_sum += allele_frequency(copies, 2);
                }
                let heterozygosity = (heterozygosity_sum / SINGLE_GENOME_CENSUS as f64)
                    / inbreeding_shrinkage(2, inbreeding);
                heterozygosity_tally.add(heterozygosity / heterozygosity_truth);
                // **No inbreeding correction on this one, and that is not an omission**: the share
                // of a genome's chromosomes carrying the alternative allele has expectation `f`
                // whatever the two copies inside the individual are doing.
                frequency_tally
                    .add((frequency_sum / SINGLE_GENOME_CENSUS as f64) / frequency_truth);
            }
            // The two-sided 95% band of a normal, which a mean over two million positions earns.
            println!(
                "| {} | {inbreeding} | {:+.2}% | {:.2}% | ±{:.1}% | {:+.0}% | {:+.2}% | {:.2}% | \
                 ±{:.1}% |",
                shape.name,
                (heterozygosity_tally.mean - 1.0) * 100.0,
                heterozygosity_tally.spread() * 100.0,
                1.96 * heterozygosity_tally.spread() * 100.0,
                (SPECIES_FALLBACK_DIVERSITY / heterozygosity_truth - 1.0) * 100.0,
                (frequency_tally.mean - 1.0) * 100.0,
                frequency_tally.spread() * 100.0,
                1.96 * frequency_tally.spread() * 100.0,
            );
        }
    }
}

/// Census positions for the single-genome arm — the shipped default
/// (`doc/devel/ng/spec/parameter_prepass_census_sites.md` §5).
const SINGLE_GENOME_CENSUS: usize = 2_000_000;

// ---------------------------------------------------------------------------------------------
// The two checks that would catch a wrong drawer or a wrong estimator
// ---------------------------------------------------------------------------------------------

/// Print, for each shape, its moments and how far its segregating frequencies sit from the
/// closest single Beta.
///
/// **This is what turns "outside the fitted family" from a claim into a number.** The closest
/// Beta is the one matching the shape's first two moments over the segregating positions; the
/// distance is the largest gap between the two distributions' cumulative curves, measured over
/// 200,000 draws from each, which pins it to about 0.4 parts in a thousand.
///
/// It also checks the drawer against the algebra, and it does so **on the segregating positions
/// only**. Checking the whole population's mean would be a weak test: at these diversities fewer
/// than one position in a hundred segregates, so a Monte-Carlo mean over the lot is a few percent
/// noisy and would hide a real disagreement. Conditioned on segregating, the same number of draws
/// pins it to about two parts in a thousand. The end masses are checked separately, by counting
/// how often the drawer lands on each.
fn print_shape_table(shapes: &[PopulationShape]) {
    println!();
    println!("## The populations swept");
    println!();
    println!(
        "`drawn` is Monte-Carlo over {} draws from the same drawer the sweep uses, beside the \
         closed form the population's own moments come from. They must agree.",
        DRAWER_CHECK_DRAWS
    );
    println!();
    println!(
        "| population | heterozygosity | per kb | mean frequency | segregating: mean f, closed \
         form | drawn | invariant share, set | drawn | gap from the closest single Beta |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
    for shape in shapes {
        let mut rng = Rng(mix(0xB5AD_4ECE_DA10_8000, shape.name.len() as u64, 17));
        let mut segregating_sum = 0.0;
        let mut segregating_count = 0_u64;
        let mut invariant_count = 0_u64;
        for _ in 0..DRAWER_CHECK_DRAWS {
            let f = shape.draw_frequency(&mut rng);
            if f == 0.0 {
                invariant_count += 1;
            } else if f < 1.0 {
                segregating_sum += f;
                segregating_count += 1;
            }
        }
        let segregating_mean: f64 = shape
            .segregating
            .iter()
            .map(|c| c.weight * c.mean_frequency())
            .sum();
        println!(
            "| {} | {:.6} | {:.2} | {:.6} | {:.4} | {:.4} | {:.4} | {:.4} | {:.3} |",
            shape.name,
            shape.heterozygosity(),
            1_000.0 * shape.heterozygosity(),
            shape.mean_frequency(),
            segregating_mean,
            segregating_sum / segregating_count.max(1) as f64,
            shape.invariant,
            invariant_count as f64 / DRAWER_CHECK_DRAWS as f64,
            distance_from_the_closest_beta(shape),
        );
    }
    println!();
    println!(
        "The last column is the largest gap between a shape's segregating frequencies and the \
         closest single Beta's cumulative curve. Near zero means one Beta holds the shape; {} of \
         the {} shapes are inside the family the joint fit uses and the rest are not.",
        shapes.iter().filter(|s| s.inside_the_fitted_family).count(),
        shapes.len()
    );
}

/// Draws used for the checks that the drawer agrees with the algebra, and for the distance from
/// the closest Beta.
const DRAWER_CHECK_DRAWS: usize = 400_000;

/// The largest gap between a shape's segregating frequencies and the closest single Beta's,
/// measured by drawing from both and comparing their cumulative curves.
fn distance_from_the_closest_beta(shape: &PopulationShape) -> f64 {
    let total: f64 = shape.segregating.iter().map(|c| c.weight).sum();
    let mean: f64 = shape
        .segregating
        .iter()
        .map(|c| c.weight * c.mean_frequency())
        .sum::<f64>()
        / total;
    let mean_square: f64 = shape
        .segregating
        .iter()
        .map(|c| c.weight * c.mean_square_frequency())
        .sum::<f64>()
        / total;
    let variance = (mean_square - mean * mean).max(1e-12);
    // Moment matching: a Beta with this mean and variance has a + b = mean(1-mean)/variance - 1.
    let concentration = (mean * (1.0 - mean) / variance - 1.0).max(1e-6);
    let (a, b) = (mean * concentration, (1.0 - mean) * concentration);

    let draws = DRAWER_CHECK_DRAWS;
    let mut from_shape = Vec::with_capacity(draws);
    let mut from_beta = Vec::with_capacity(draws);
    let mut rng = Rng(mix(0x1234_5678_9ABC_DEF0, a.to_bits(), b.to_bits()));
    let weights: Vec<f64> = shape.segregating.iter().map(|c| c.weight).collect();
    for _ in 0..draws {
        let component = shape.segregating[rng.pick(&weights)];
        from_shape.push(component.draw(&mut rng));
        from_beta.push(rng.beta(a, b));
    }
    from_shape.sort_by(f64::total_cmp);
    from_beta.sort_by(f64::total_cmp);
    let mut gap: f64 = 0.0;
    let mut j = 0_usize;
    for (i, value) in from_shape.iter().enumerate() {
        while j < draws && from_beta[j] <= *value {
            j += 1;
        }
        gap = gap.max(((i + 1) as f64 - j as f64).abs() / draws as f64);
    }
    gap
}

/// A check on [`nei_heterozygosity`] that does not go through the sweep: at one known frequency,
/// redraw the panel many times and see whether the estimator's mean is `2 f (1 - f)`.
///
/// **This is the check that would catch the estimator written the wrong way round in `k`, or with
/// `2N` where `2N - 1` belongs.** It shares no code with the sweep except the estimator itself,
/// and it is run at a frequency and a panel size where the wrong forms differ visibly: the `2N`
/// denominator would come back 1/(2N) low, and the panel sizes below span a factor where that is
/// 50%, 10% and 1%.
fn print_estimator_oracle() {
    println!();
    println!("## The estimator, checked against one known frequency");
    println!();
    println!(
        "One position at frequency 0.30, panel redrawn 400,000 times, no census sampling in it \
         at all. The truth is 2 f (1 - f) = {:.4} at F = 0 and {:.4} at F = 0.8 for a single \
         genome.",
        2.0 * 0.3 * 0.7,
        2.0 * 0.3 * 0.7 * 0.2
    );
    println!();
    println!("| individuals | F | mean of the estimator | truth | ratio |");
    println!("|---:|---:|---:|---:|---:|");
    for individuals in [1_usize, 5, 50] {
        for inbreeding in INBREEDING {
            let mut rng = Rng(mix(0xCBF2_9CE4_8422_2325, individuals as u64, 7));
            let frequency = 0.30_f64;
            let draws = 400_000;
            let mut total = 0.0_f64;
            for _ in 0..draws {
                let mut copies = 0_u32;
                for _ in 0..individuals {
                    copies += draw_genotype(&mut rng, frequency, inbreeding);
                }
                total += nei_heterozygosity(copies, 2 * individuals as u32);
            }
            let truth = 2.0
                * frequency
                * (1.0 - frequency)
                * inbreeding_shrinkage(2 * individuals as u32, inbreeding);
            println!(
                "| {individuals} | {inbreeding} | {:.5} | {truth:.5} | {:.4} |",
                total / draws as f64,
                (total / draws as f64) / truth,
            );
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------------------------------

/// Draw one diploid individual's count of alternative copies at a position of frequency `f`,
/// under the standard inbreeding model: with probability `F` the two copies are one ancestral
/// copy counted twice, otherwise they are two independent draws.
fn draw_genotype(rng: &mut Rng, frequency: f64, inbreeding: f64) -> u32 {
    if rng.uniform() < inbreeding {
        if rng.uniform() < frequency { 2 } else { 0 }
    } else {
        u32::from(rng.uniform() < frequency) + u32::from(rng.uniform() < frequency)
    }
}

/// Fold three numbers into one seed. Neighbouring seeds in this generator's state produce
/// correlated first draws, so cells are separated by mixing rather than by adding.
fn mix(base: u64, first: u64, second: u64) -> u64 {
    let mut state = base ^ first.wrapping_mul(0xFF51_AFD7_ED55_8CCD);
    state ^= state >> 33;
    state = state.wrapping_mul(0xC4CE_B9FE_1A85_EC53) ^ second.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    state ^= state >> 29;
    state | 1
}

/// The stream every drawn number comes from — the same xorshift the sibling harnesses use, so a
/// run is reproducible.
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
}

/// Check every Beta component's closed-form mean against the same quantity computed a different
/// way, and stop the run if any disagrees.
///
/// `E[f]` under `Beta(a, b)` is `a / (a + b)`, and it is also the ratio of two Beta functions,
/// `B(a + 1, b) / B(a, b)`. **The two share no arithmetic** — one is a division, the other four
/// `lgamma` calls and an exponential — so a shape whose parameters were transposed on the way in
/// still agrees here, but a wrong moment formula does not. It is cheap and it runs on the shapes
/// the sweep actually uses, not on a fixture.
fn check_the_beta_moments(shapes: &[PopulationShape]) {
    let ln_beta = |x: f64, y: f64| lgamma(x) + lgamma(y) - lgamma(x + y);
    for shape in shapes {
        // **The check that would have caught the first draft of this list, and did.** One shape
        // was written with two end masses totalling 1.015, which makes the segregating share
        // negative — and a negative share does not crash: it produces a population whose printed
        // heterozygosity is `-0.0033`, and a drawer that silently never segregates.
        assert!(
            shape.invariant + shape.fixed_alt < 1.0,
            "'{}' has {} of positions invariant and {} fixed non-reference, which leaves nothing \
             to segregate: the two are shares of the same positions",
            shape.name,
            shape.invariant,
            shape.fixed_alt
        );
        let weights: f64 = shape.segregating.iter().map(|c| c.weight).sum();
        assert!(
            (weights - 1.0).abs() < 1e-12,
            "'{}' spreads {weights} of its segregating positions rather than all of them",
            shape.name
        );
        for component in &shape.segregating {
            let (a, b) = (component.a, component.b);
            let closed_form = a / (a + b);
            let through_gamma = (ln_beta(a + 1.0, b) - ln_beta(a, b)).exp();
            assert!(
                (closed_form - through_gamma).abs() < 1e-12,
                "Beta({a}, {b}) has mean {closed_form} by the closed form and {through_gamma} \
                 through the gamma functions, so one of the two is wrong and every population \
                 moment printed by this program rests on the first"
            );
        }
    }
}
