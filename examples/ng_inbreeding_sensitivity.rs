//! **How wrong is the caller's diversity when the panel's inbreeding coefficient is wrong?**
//!
//! Where the pre-pass supplies a diversity rather than a fitted population curve — the per-sample
//! histogram route — the run's `θ` is each sample's *observed* heterozygosity divided by `(1 − F)`
//! (`doc/devel/ng/spec/calling_priors.md` §4). **In a selfer that divides by a small number.** At
//! an inbreeding coefficient of 0.85 only about 15 alternative copies in 100 sit in heterozygotes,
//! so the diversity is measured through the thinnest available channel and multiplied by 6.7 to
//! recover it; an error `dF` in the coefficient moves `θ` by about `dF / (1 − F)`.
//!
//! **That is the first-order size and the table below prints it; the exact one is
//! `dF / (1 − F − dF)`, and the two part company where they matter most.** At `F = 0.9` with a
//! coefficient 0.05 too high, the column says 50% and the true inflation is **100%** — because
//! `1 − F` has halved rather than shrunk by a twentieth. The first-order figure is what a reader
//! can carry between rows; the exact one is what a run at the wrong `F` actually gets.
//!
//! That factor is arithmetic rather than a simulation, and printing it is all this program now
//! does.
//!
//! **⛔ It used to measure a second, larger route, and that route is gone.** Until 2026-08-27 the
//! seed's shape was fitted to the panel's own allele-count classes under a two-branch model that
//! takes an inbreeding coefficient, so a wrong `F` reached the run's two starting numbers directly.
//! This program built a panel's exact allele-count distribution at a known coefficient, refitted it
//! at a wrong one, and reported how far the pair moved. **The seed's numbers are now two integrals
//! of the fitted population curve** (`doc/devel/ng/spec/ordinary_site_prior_moments.md` §2), which
//! carry no panel and no coefficient, so there is nothing left for that arm to move.
//!
//! **What that arm found, kept because it is the reason the route was safe to remove:** with the
//! total already pinned to the measured heterozygosity, an `F` wrong by ±0.10 moved the shipped
//! pair **not at all at one individual**, to five digits, while the search's own reference
//! concentration moved by a factor of three; at 63 individuals it moved the shipped number by at
//! most 1.4% against 4% for the search's. The full tables are in
//! `doc/devel/ng/reports/inbreeding_sensitivity_of_the_seed_2026-08-23.md`.
//!
//! Run with:
//!
//! ```text
//! ./scripts/dev.sh cargo run --release --example ng_inbreeding_sensitivity
//! ```

fn main() {
    println!(
        "Where the pre-pass supplies a diversity rather than a fitted population curve, the pair\n\
         is (1, th), and th is each sample's observed heterozygosity divided by (1 - F). An error\n\
         in F then moves th by about dF / (1 - F). The exact size is dF / (1 - F - dF), which\n\
         at F = 0.9 and dF = 0.05 is 100% against the 50% printed here:\n"
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
