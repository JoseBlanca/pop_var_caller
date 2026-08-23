//! **How wrong is the caller's starting belief when the panel's inbreeding is wrong?**
//!
//! The genotype prior starts every locus from two numbers — one attached to the reference base,
//! one shared among the alternatives — and both are now read off what the panel actually showed
//! (`doc/devel/ng/spec/calling_priors.md` §4.1). The panel's inbreeding coefficient enters that
//! reading twice, and this measures both:
//!
//! 1. **Directly, in the prediction.** Working out what allele counts a candidate pair of numbers
//!    would produce needs a model of how the panel's chromosomes came about, and in a selfer an
//!    individual's two copies are usually one ancestral copy inherited twice. Give that model the
//!    wrong coefficient and the pair it recovers is wrong. Measured here by building the exact
//!    allele-count distribution of a panel at a known coefficient and then fitting it at a wrong
//!    one.
//!
//! 2. **Indirectly, through the diversity.** Where the panel is too small for the pre-pass to
//!    emit an allele-count distribution at all — one sample, or a cohort below its floor — the
//!    pair falls back to `(1, θ)`, and `θ` is each sample's observed heterozygosity divided by
//!    `(1 − F)` (spec §4). **In a selfer that divides by a small number.** At an inbreeding
//!    coefficient of 0.85 only about 15 alternative copies in 100 sit in heterozygotes, so the
//!    diversity is measured through the thinnest available channel and multiplied by 6.7 to
//!    recover it. That factor is exact arithmetic rather than a simulation, and it is printed
//!    beside the fitted results.
//!
//! **Read the one-sample rows as a bound, not as a prediction.** At one sample no site can vary
//! *across* the panel, so the distribution the fit is handed is the pre-pass's own neutral prior,
//! built from the same coefficient the fit then uses. A coefficient that is wrong is therefore
//! wrong at both ends and the two errors cancel exactly in the reference number — what survives
//! is the diversity, by the factor in the second table. The rows below hold the coefficient wrong
//! at one end only, which is what a **cohort** meets: it has one coefficient standing for samples
//! that each have their own, and no document yet says how that one is arrived at.
//!
//! **And what an error costs, in the units that move a call.** The reference number sets the
//! prior odds on a heterozygote against a homozygous-variant call before any read is looked at —
//! 2:1 when it is 1. Every row below carries those odds, computed by running the shipped prior,
//! so an error in the concentration can be read as an error in the call.
//!
//! Run with:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example ng_inbreeding_sensitivity
//! ```

use pop_var_caller::ng::calling::GenotypeTable;
use pop_var_caller::ng::calling::genotype_prior::seed_generic::{
    FittedSpectrum, fill_expected_spectrum, project_spectrum_seed,
};
use pop_var_caller::ng::calling::genotype_prior::{
    Concentration, GenotypePriorModel, MarginalizedDirichletPrior, PriorRow,
};
use pop_var_caller::ng::types::{ExpectedHeterozygosity, InbreedingF, LogProb, Ploidy};

/// Tomato's fitted diversity: 6 differences per 10,000 bases (spec §4.1).
const THETA: f64 = 6e-4;

/// The exact allele-count distribution a panel of this many diploid individuals would show if its
/// population frequencies came from the concentration pair given — the same closed form the fit
/// searches, so the truth here is known rather than simulated.
fn exact_spectrum(alpha_ref: f64, alpha_alt: f64, individuals: u32, inbreeding: f64) -> Vec<f64> {
    let mut out = vec![0.0; 2 * individuals as usize + 1];
    fill_expected_spectrum(
        alpha_ref,
        alpha_alt,
        individuals,
        InbreedingF::try_new(inbreeding).unwrap(),
        &mut out,
    );
    out
}

/// The prior odds on a heterozygote against a homozygous-variant call, at a biallelic diploid
/// site, before any read — read off the shipped prior rather than off a formula.
fn heterozygote_odds(alpha_ref: f64, alpha_alt: f64) -> f64 {
    let table = GenotypeTable::build(Ploidy::try_new(2).unwrap(), 2);
    let view = table.view();
    let concentration = [alpha_ref, alpha_alt];
    let mut scratch = vec![0.0; 2];
    let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
    let mut row = PriorRow::new(
        Concentration::new(&concentration),
        view.genotype_allele_counts(),
        view.log_multinomial_coeffs(),
        view.homozygous_alleles(),
        &mut scratch,
        &mut out,
    );
    // At F = 0, so what the row shows is the concentration's own doing and not the mixture's.
    MarginalizedDirichletPrior
        .fill_genotype_log_priors(&mut row, InbreedingF::try_new(0.0).unwrap());
    (out[1].get() - out[2].get()).exp()
}

fn main() {
    println!(
        "How far the run's two starting numbers move when the panel's inbreeding coefficient is\n\
         wrong. The panel really is at `F true`; the fit is told `F used`. Truth is (1, {THETA:e}).\n"
    );
    println!(
        "{:>6}  {:>6}  {:>6}  {:>10}  {:>10}  {:>9}  {:>8}",
        "people", "F true", "F used", "alpha_ref", "alpha_alt", "alt / th", "het:hom"
    );

    for individuals in [1u32, 26, 63] {
        for f_true in [0.6f64, 0.85] {
            let weights = exact_spectrum(1.0, THETA, individuals, f_true);
            for delta in [-0.10f64, -0.05, 0.0, 0.05, 0.10] {
                let f_used = (f_true + delta).clamp(0.0, 0.999);
                let seed = project_spectrum_seed(
                    Some(FittedSpectrum::new(&weights, 10.0, 3_000.0)),
                    Some(ExpectedHeterozygosity::try_new(THETA).unwrap()),
                    InbreedingF::try_new(f_used).unwrap(),
                );
                println!(
                    "{individuals:>6}  {f_true:>6.2}  {f_used:>6.2}  {:>10.4}  {:>10.3e}  {:>9.4}  {:>7.3}:1",
                    seed.alpha_ref(),
                    seed.alpha_alt_total(),
                    seed.alpha_alt_total() / THETA,
                    heterozygote_odds(seed.alpha_ref(), seed.alpha_alt_total()),
                );
            }
            println!();
        }
    }

    println!(
        "Where the pre-pass emits no allele-count distribution — one sample, or a cohort below its\n\
         floor — the pair is (1, th) instead, and th is each sample's observed heterozygosity\n\
         divided by (1 - F). An error in F then moves th by dF / (1 - F), exactly:\n"
    );
    println!(
        "{:>6}  {:>22}  {:>28}",
        "F", "copies in het (of 100)", "th error from dF = 0.05"
    );
    for f in [0.0f64, 0.6, 0.85, 0.9] {
        // Under Hardy-Weinberg with inbreeding and a rare allele, heterozygotes carry a fraction
        // (1 - F) of the alternative copies and homozygotes the rest.
        let in_heterozygotes = 100.0 * (1.0 - f);
        let theta_error = 0.05 / (1.0 - f);
        println!(
            "{f:>6.2}  {in_heterozygotes:>22.0}  {:>27.0}%",
            100.0 * theta_error
        );
    }
}
