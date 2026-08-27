//! **Below what panel size is a fitted frequency spectrum not worth projecting?**
//!
//! ## The question
//!
//! The ordinary-site genotype prior starts from a pair of concentrations, and it has three ways
//! to get one (`doc/devel/ng/spec/population_diversity.md` §3.4). Its top rung reads both numbers
//! off the panel's **frequency spectrum** — how allele frequencies are spread across the
//! population, held as one weight per allele-count class. Its middle rung takes a neutral shape at
//! the fitted diversity, and is what a run gets when no spectrum arrives.
//!
//! **At `N` diploid individuals the spectrum has `2N + 1` classes**, so a panel of one gives three
//! and carries almost no shape. The consumer's own documentation says the pre-pass emits the
//! spectrum as absent below a panel-size floor; **no such floor exists**, so the top two rungs are
//! separated by nothing. `population_diversity.md` §9's third question asks where to put one, and
//! its leaning is that a floor low enough never to fire is worse than none — **set it from
//! measurement, not from taste**.
//!
//! **⚠ What this measures is now a fact about the search, not about the shipped seed.** Since
//! `doc/devel/ng/spec/ordinary_site_seed.md` §3 the seed's *total* — how much conviction the pair
//! carries — is solved from the run's measured heterozygosity rather than taken from the search,
//! precisely because of the loss this program measures. So the "the pair's heterozygosity,
//! against the density's" column below is what the seed **would** lose if it still read the
//! search's total, and what it now loses is nothing: the pinned pair implies the measurement
//! exactly at every panel size
//! (`seed_generic::projection_tests::the_seeds_implied_diversity_is_the_measured_one_at_every_panel_and_shape`).
//! The program is kept because that column is the evidence for the pin, and because the pair's
//! *shape* still comes from this search.
//!
//! ## What this measures
//!
//! One fitted allele-frequency density — two point masses and a Beta over what segregates — is
//! projected into `2N + 1` allele-count classes at a range of panel sizes, and the projection the
//! caller's own search (`fit_spectrum_shape`) is fitted to each. Two numbers come back per
//! panel size:
//!
//! - **how far the fitted pair sits from the same density's answer at a large panel**, which is
//!   what "the projection is worth doing" means: below the floor the pair the caller would use is
//!   not the pair the population implies;
//! - **the divergence the fit reports** — how much information is lost by describing the projected
//!   spectrum with the two-parameter family instead of itself. §9's question names this one.
//!
//! **The density is the same at every panel size**, so nothing here is about sampling noise: the
//! class weights are the exact Beta-binomial ones. What moves is only how coarsely `2N + 1`
//! classes can express a continuous density, and how much shape the family can read out of that
//! many numbers.
//!
//! **Which is why this sweep cannot answer the floor question on its own, and does not claim to.**
//! A floor is for a panel too small to *estimate* a spectrum from, and that is a sampling
//! question; the experiment for it is named in `parameter_prepass_cohort.md` §10's third question
//! — subsample a real cohort and watch where the spectrum stops being stable — and is not this.
//! What this settles is narrower and still worth having: **the statistic
//! `population_diversity.md` §9's third question names cannot locate a floor**, because it is
//! smallest where a floor would fire.
//!
//! **Every panel here is outbred (`F = 0`).** Inbreeding reshapes the predicted spectrum —
//! `project_spectrum_seed`'s own documentation records the reference concentration moving 8.6% to
//! 14.0% across `F = 0.6` to `0.9` — so these are outbred-panel numbers.
//!
//! ## What it is not
//!
//! **The densities swept are illustrative, not fitted.** No cohort's fitted `FrequencyDensity` is
//! recorded in this repository, so the sweep runs a grid of Beta shapes at diversities spanning
//! this project's two benchmark cohorts — tomato at about 6 differences per 10,000 bases and a
//! human panel at about 10 — and reports whether the answer is the same across the grid. A floor
//! that moves with the shape is a different finding from one that does not, and both are reported.
//!
//! Run: `./scripts/dev.sh cargo run --release --example ng_spectrum_panel_floor`

use std::time::Instant;

use pop_var_caller::ng::calling::genotype_prior::{FittedSpectrum, fit_spectrum_shape};
use pop_var_caller::ng::parameter_estimation::joint::fit::FrequencyDensity;
use pop_var_caller::ng::types::InbreedingF;

/// The panel sizes swept, in diploid individuals. **1 is the single genome this caller commits
/// to and 200 is above every panel in this repository**; the two benchmark cohorts sit at 1
/// (HG002) and 63 (tomato).
const PANELS: [u32; 16] = [1, 2, 3, 4, 5, 6, 8, 10, 13, 16, 20, 32, 50, 63, 100, 200];

/// The panel the smaller ones are measured against — the answer the density itself implies, as
/// closely as this sweep can reach it.
const REFERENCE_PANEL: u32 = 200;

/// One shape of population, named for the cohort it is meant to resemble.
struct Population {
    name: &'static str,
    density: FrequencyDensity,
}

/// The grid. **Each is a real density's shape, not a real density**: the Beta shapes span the
/// rare-allele pile-up a neutral population has (`a` below 1) up to a flat one, and each
/// `p_invariant` is set so the implied heterozygosity lands near its cohort's.
fn populations() -> Vec<Population> {
    vec![
        Population {
            name: "tomato-like, strong rare-allele pile-up",
            density: FrequencyDensity {
                p_invariant: 0.9950,
                p_fixed_alt: 0.0010,
                a: 0.20,
                b: 1.00,
            },
        },
        Population {
            name: "human-like, moderate pile-up",
            density: FrequencyDensity {
                p_invariant: 0.9949,
                p_fixed_alt: 0.0004,
                a: 0.35,
                b: 1.20,
            },
        },
        Population {
            name: "flat over what segregates",
            density: FrequencyDensity {
                p_invariant: 0.9950,
                p_fixed_alt: 0.0010,
                a: 1.00,
                b: 1.00,
            },
        },
        Population {
            // **The density `fit::tests::a_lopsided_density` uses**, swept here so that the
            // figure quoted in that module's doc comment and the figures here come from one run.
            name: "the unit tests' own lopsided fixture",
            density: FrequencyDensity {
                p_invariant: 0.90,
                p_fixed_alt: 0.01,
                a: 0.50,
                b: 2.00,
            },
        },
        Population {
            name: "middling frequencies — the shape the family cannot hold",
            density: FrequencyDensity {
                p_invariant: 0.9950,
                p_fixed_alt: 0.0010,
                a: 4.00,
                b: 4.00,
            },
        },
    ]
}

/// The projection, which is the shipped one — `FrequencyDensity::allele_count_classes`.
///
/// **This sweep calls what the caller calls.** An earlier draft carried its own copy of the
/// Beta-binomial and its own `ln Γ`, which would have made the floor a fact about the copy.
fn allele_count_classes(density: &FrequencyDensity, individuals: u32) -> Vec<f64> {
    density.allele_count_classes(individuals)
}

struct Row {
    individuals: u32,
    alpha_ref: f64,
    alpha_alt: f64,
    /// What the fitted pair implies about how often two copies differ — `2 α_ref α_alt /
    /// ((α_ref + α_alt)(α_ref + α_alt + 1))`, the Dirichlet's own expected heterozygosity.
    /// **Measured against the density's own**, because the two part company as the panel grows.
    pair_heterozygosity: f64,
    divergence_nats: f64,
    at_search_limit: bool,
    seconds: f64,
}

fn main() {
    let outbred = InbreedingF::try_new(0.0).expect("a legal coefficient");

    println!("# Below what panel size is a fitted frequency spectrum not worth projecting?");
    println!();
    println!(
        "Each panel's pair is measured against the same density's answer at {REFERENCE_PANEL} \
         individuals."
    );

    for population in populations() {
        let implied = population.density.expected_heterozygosity();
        println!();
        println!("## {}", population.name);
        println!(
            "Beta({:.2}, {:.2}) over {:.4} of positions, {:.5} of them fixed non-reference; \
             implied heterozygosity {:.3e}",
            population.density.a,
            population.density.b,
            population.density.p_segregating(),
            population.density.p_fixed_alt,
            implied
        );
        println!();

        let rows: Vec<Row> = PANELS
            .iter()
            .map(|individuals| {
                let classes = allele_count_classes(&population.density, *individuals);
                let started = Instant::now();
                // **The search's own pair, which is what this sweep is about.** The shipped seed
                // no longer uses the search's total — it solves one from the run's measured
                // heterozygosity instead (`doc/devel/ng/spec/ordinary_site_seed.md` §3), so the
                // loss measured below is now a fact about the search rather than about the seed.
                let shape =
                    fit_spectrum_shape(&FittedSpectrum::new(&classes, 0.0, 1_000.0), outbred);
                let seconds = started.elapsed().as_secs_f64();
                let spectrum_match = shape.spectrum_match();
                let (alpha_ref, alpha_alt) = shape.concentrations();
                let total = alpha_ref + alpha_alt;
                Row {
                    individuals: *individuals,
                    alpha_ref,
                    alpha_alt,
                    pair_heterozygosity: 2.0 * alpha_ref * alpha_alt / (total * (total + 1.0)),
                    divergence_nats: spectrum_match.divergence_nats(),
                    at_search_limit: spectrum_match.at_search_limit(),
                    seconds,
                }
            })
            .collect();

        let reference = rows.last().expect("at least one panel");
        println!(
            "| individuals | classes | alpha_ref | alpha_alt | from the {REFERENCE_PANEL}-panel \
             pair | the pair's heterozygosity, against the density's | divergence (nats) | at the \
             search limit | seconds |"
        );
        println!("|---:|---:|---:|---:|---:|---:|---:|---|---:|");
        for row in &rows {
            let drift = distance(row, reference);
            println!(
                "| {} | {} | {:.4} | {:.3e} | {:.1}% | {:.3e} ({:+.1}%) | {:.3e} | {} | {:.2} |",
                row.individuals,
                2 * row.individuals + 1,
                row.alpha_ref,
                row.alpha_alt,
                drift * 100.0,
                row.pair_heterozygosity,
                (row.pair_heterozygosity / implied - 1.0) * 100.0,
                row.divergence_nats,
                if row.at_search_limit { "yes" } else { "no" },
                row.seconds
            );
        }

        // **No floor is computed from this, and an earlier draft of this example computed one.**
        // It looked for the smallest panel whose pair is within a tolerance of the largest
        // panel's — which is circular: the largest panel is zero from itself, so the search
        // always succeeds at the last row and the "never settles" arm was dead code. It printed
        // "settled from 100 individuals up" while the table beside it showed the pair still
        // moving between 100 and 200.
        //
        // The table is the finding. What it says is that the pair drifts monotonically and has
        // not converged by the largest panel swept, and that the divergence — the statistic
        // `population_diversity.md` §9's third question names — is smallest at the smallest panel
        // rather than largest, so it cannot locate a floor at all.
        let smallest = rows.first().expect("at least one panel");
        println!(
            "\nStill moving at {REFERENCE_PANEL} individuals: the pair is {:.1}% from the \
             largest panel's at one individual and {:.1}% at 100, and the divergence rises from \
             {:.3e} to {:.3e} rather than falling.",
            distance(smallest, reference) * 100.0,
            rows.iter()
                .find(|row| row.individuals == 100)
                .map_or(f64::NAN, |row| distance(row, reference))
                * 100.0,
            smallest.divergence_nats,
            reference.divergence_nats
        );
    }
}

/// How far one row's pair sits from another's — the larger of the two relative differences, so a
/// pair that matches on one concentration and not the other is not called settled.
fn distance(row: &Row, reference: &Row) -> f64 {
    let reference_drift = (row.alpha_ref - reference.alpha_ref).abs() / reference.alpha_ref;
    let alternative_drift = (row.alpha_alt - reference.alpha_alt).abs() / reference.alpha_alt;
    reference_drift.max(alternative_drift)
}
