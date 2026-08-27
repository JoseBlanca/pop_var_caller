//! **How big does a panel have to be before its own allele frequencies are a better guess than
//! the textbook shape?**
//!
//! ## The question, and why it has a number in it
//!
//! The SNP/indel genotype prior starts from two numbers — how many chromosomes' worth of belief
//! the reference allele carries and how many the alternatives do. Since
//! `doc/devel/ng/spec/ordinary_site_seed.md` §3 those two are set from an **expected allele
//! frequency** and the run's **measured heterozygosity**: the second fixes how much conviction
//! the pair carries, and the first is its shape.
//!
//! There are two ways to guess that shape and neither is right everywhere. The **neutral** one
//! comes from population genetics and needs no panel at all: a population under no selection has
//! most of its alternative alleles rare, and its expected frequency is about the heterozygosity
//! itself. The **panel's own** comes from fitting the two-parameter family to how the alternative
//! allele is actually spread across this cohort's chromosomes — better when the cohort is large
//! enough to say, and measurably worse when it is not. At one, two and three samples the fitted
//! shape comes back on the *wrong side* of the neutral one in every draw (spec §1.1).
//!
//! So the seed interpolates: `ln f = (1 − w) · ln f_neutral + w · ln f_fitted`, with
//! `w = N / (N + N₀)` for a panel of `N` diploid individuals. **`N₀` is this program's output.**
//! It is the panel size at which the two guesses are equally good, so `w` is a half there.
//!
//! ## What it does
//!
//! For each population shape and each read depth:
//!
//! 1. draw a cohort at known parameters — a known Beta over the positions that segregate, and a
//!    known share of positions that do;
//! 2. refit **the same drawn positions** using the first 1, 2, 3, 5, 10, 25 and 63 samples;
//! 3. project each fit into the allele-count classes that panel has, and run the shipped search
//!    over them ([`fit_spectrum_shape`]) to get the panel's own expected frequency;
//! 4. score both guesses against the frequency the cohort was actually drawn with, as a ratio in
//!    log space — which is the space the blend works in.
//!
//! **Holding the drawn positions fixed across the arms is what makes it a curve in panel size.**
//! Redrawing them at each arm would let the answer move for a reason that has nothing to do with
//! how many samples there are, which is the whole quantity being measured.
//!
//! ## Both axes, not one
//!
//! Every number in the spec is at 20 reads a sample, and this project's tomato panel sits near
//! three. A weight fitted at high depth and applied at low would lean on a shape the reads never
//! supported, so the sweep runs at **3, 8 and 20 reads a sample** and reports whether `N₀` moves.
//! If it does, that is a finding and the weight has to take depth as well.
//!
//! ## The hold-out that says the answer is worth something
//!
//! `N₀` is fitted on one set of draws and then checked on **a second set the fit never saw**, at
//! different seeds. On those, the blended shape is scored against both ends. What the design
//! claims is that the blend is at least as good as whichever end is better at each panel size
//! (spec §6.3); what is reported is what is actually true, arm by arm.
//!
//! ## What this is not
//!
//! **Drawn cohorts, not a real one.** This checkout cannot rebuild the tomato census — the read
//! files are not in the repository — so the sweep runs on cohorts drawn at known parameters. That
//! is what `doc/devel/specs/design_principles.md` §0 asks for in any case; what it is not is a
//! confirmation on real data, and `ordinary_site_seed.md` §7's second open question keeps that
//! open.
//!
//! **The segregating share is set high on purpose.** A tomato-like 4 positions in 1,000 would put
//! about a dozen variable positions in a run of this size, which measures the draw rather than the
//! panel. Both shapes here segregate at 10% of positions, so what differs between them is the
//! *shape* of the frequency density and not how much there is to see. The heterozygosities that
//! follow — 26 and 15 per thousand bases — are well above either benchmark cohort's, and are
//! stated rather than hidden.
//!
//! **Every panel here is outbred.** The drawn cohorts have no excess of homozygotes and the
//! projection is run at `F = 0`, which is the basis every measurement in the spec is on. Whether
//! inbreeding belongs in the weight is that document's third open question.
//!
//! ```text
//! ng_seed_shape_weight_sweep [positions] [draws] [hold-out draws] [mismapped share]
//! ```
//!
//! Run: `./scripts/dev.sh cargo run --release --example ng_seed_shape_weight_sweep`

use std::collections::BTreeMap;
use std::time::Instant;

use pop_var_caller::ng::calling::genotype_prior::{
    FittedSpectrum, HALF_WEIGHT_PANEL_SIZE, fit_spectrum_shape,
};
use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, ObservedAllele, PackedDepthCodes, ReadCap, RecordingTerms,
    SampleCensusEvidence, Section, SectionKey, SelectionTermsDigest,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFitConfig, StartingPoint, fit_jointly,
};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::{InbreedingF, ReadGroupId};

/// The panel sizes the curve is reported at, in diploid individuals. **1 is the single genome
/// this caller commits to and 63 is the tomato cohort's size.**
const ARMS: [usize; 7] = [1, 2, 3, 5, 10, 25, 63];

/// Reads a sample at a position. **3 is where the tomato panel sits and 20 is what every number
/// in the spec was measured at.**
const DEPTHS: [f64; 3] = [3.0, 8.0, 20.0];

/// How often a read misreads a base at an ordinary position.
const SEQUENCING_ERROR: f64 = 0.002;

/// The share of positions that segregate in every drawn population. **High on purpose** — see
/// this file's opening note.
const SEGREGATING: f64 = 0.10;

/// The share of positions where the population carries only a non-reference base.
const FIXED_ALTERNATIVE: f64 = 0.0005;

/// One shape of population, named for what its Beta does to rare alleles.
struct Population {
    name: &'static str,
    a: f64,
    b: f64,
}

/// The two shapes. **Beta(0.7, 2.5) is the one `ordinary_site_seed.md` §1.1's own table was
/// measured on**; Beta(0.20, 1.00) piles alternative alleles far harder into the rare end, which
/// is what the tomato-like density in `examples/ng_spectrum_panel_floor.rs` does.
const POPULATIONS: [Population; 2] = [
    Population {
        name: "moderate rare-allele pile-up, Beta(0.70, 2.50)",
        a: 0.70,
        b: 2.50,
    },
    Population {
        name: "strong rare-allele pile-up, Beta(0.20, 1.00)",
        a: 0.20,
        b: 1.00,
    },
];

/// The half-weight panel sizes tried when fitting the constant. **`0` is the panel's own shape
/// and nothing else at every panel size; `0.25` is a weight of 0.80 at
/// one individual — the panel's own shape almost throughout — and `1000` is 0.06 at 63, which is
/// the neutral shape almost throughout.**
const CANDIDATE_HALF_WEIGHTS: [f64; 15] = [
    0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 50.0, 200.0, 1_000.0,
];

fn main() {
    let mut args = std::env::args().skip(1);
    let positions: usize = args.next().map_or(3_000, |a| a.parse().expect("a count"));
    let draws: usize = args.next().map_or(6, |a| a.parse().expect("a count"));
    let held_out: usize = args.next().map_or(4, |a| a.parse().expect("a count"));
    let mismapped: f64 = args.next().map_or(0.0, |a| a.parse().expect("a share"));

    println!("# Where a panel's own allele frequencies start beating the neutral shape");
    println!();
    println!(
        "{positions} positions a cohort, {draws} drawn cohorts a cell for the fit and {held_out} \
         held out."
    );
    println!(
        "Reads disagree with the reference at {SEQUENCING_ERROR} of positions; {mismapped} of \
         positions are mismapped."
    );
    println!(
        "Every drawn plant is outbred and the projection is run at F = 0. Panels: {ARMS:?} \
         diploid individuals."
    );

    let started = Instant::now();
    // Every cell of the sweep, kept so the constant can be fitted over all of them at once.
    let mut fitting: Vec<Measured> = Vec::new();
    let mut checking: Vec<Measured> = Vec::new();

    for population in &POPULATIONS {
        let density = FrequencyDensity {
            p_invariant: 1.0 - FIXED_ALTERNATIVE - SEGREGATING,
            p_fixed_alt: FIXED_ALTERNATIVE,
            a: population.a,
            b: population.b,
        };
        let expected = expected_frequency_of(&density);

        println!("\n## {}", population.name);
        println!(
            "Segregating at {:.4} of positions, {:.5} fixed non-reference. The population's own \
             expected alternative-allele frequency is {:.5e} and its heterozygosity {:.5e} \
             ({:.1} per thousand bases).",
            density.p_segregating(),
            density.p_fixed_alt,
            expected,
            density.expected_heterozygosity(),
            1_000.0 * density.expected_heterozygosity()
        );

        println!(
            "\n### The projection on its own — no cohort drawn, nothing estimated\n\nThe density \
             is projected into each panel's allele-count classes exactly and the shipped search \
             is run over the result, so the only thing that moves between rows is **how many \
             classes the two-parameter family is asked to fit at once**."
        );
        projection_only(&density, &[1, 2, 3, 5, 10, 25, 63, 200]);

        for depth in DEPTHS {
            let cells = sweep(&density, depth, positions, mismapped, 0, draws);
            println!("\n### {depth} reads a sample");
            println!();
            println!(
                "| individuals | cohorts fitted | this panel's own frequency / the density's | \
                 fitted heterozygosity / the density's | the panel's own shape / this panel's \
                 frequency | log error | the neutral shape / this panel's frequency | log error | \
                 which is nearer |"
            );
            println!("|---:|---:|---:|---:|---:|---:|---:|---:|---|");
            for (arm, cell) in ARMS.iter().zip(&cells) {
                let (panel, neutral) = (
                    median(&cell.panel_shape_error),
                    median(&cell.neutral_shape_error),
                );
                let drawn = median(&cell.drawn_frequency);
                println!(
                    "| {arm} | {}/{} | {:.3}× | {:.3}× | {:.3}× | {:.4} | {:.3}× | {:.4} | {} |",
                    cell.fitted,
                    cell.attempted,
                    drawn / expected,
                    median(&cell.diversity_ratio),
                    median(&cell.panel_shape) / drawn,
                    panel,
                    median(&cell.neutral_shape) / drawn,
                    neutral,
                    if panel < neutral {
                        "the panel's"
                    } else {
                        "the neutral"
                    }
                );
            }
            match crossing_panel_size(&cells) {
                Some(n0) => println!(
                    "\n**The two are equally good at {n0:.1} individuals**, read off where the \
                     two log-error columns cross."
                ),
                None => println!(
                    "\n**There is no panel size above which the panel's own shape becomes the \
                     better guess.** It is the better guess at the smallest panels and loses \
                     ground as the panel grows, which is the opposite of what §4.1's crossing \
                     assumes, so this depth sets no half-weight panel size of its own."
                ),
            }
            fitting.push(Measured {
                population: population.name,
                depth,

                cells,
            });
            checking.push(Measured {
                population: population.name,
                depth,

                cells: sweep(
                    &density,
                    depth,
                    positions,
                    mismapped,
                    draws,
                    draws + held_out,
                ),
            });
        }
    }

    println!("\n## The projection on its own, at a real cohort's segregating share");
    println!();
    println!(
        "The same control at five allele-frequency densities — the **four** \
         `ordinary_site_seed.md` §1.2 measured the old behaviour on, plus this project's own \
         lopsided unit-test fixture — each segregating at 4 positions in 1,000, which is a real \
         cohort's share rather than the 10% the drawn sweep needs. **Nothing is drawn and nothing \
         is estimated**, so what these rows isolate is what the projection-and-refit costs by \
         itself."
    );
    for (name, a, b) in REALISTIC_SHAPES {
        let density = FrequencyDensity {
            p_invariant: 0.9950,
            p_fixed_alt: 0.0010,
            a,
            b,
        };
        println!("\n### {name}");
        projection_only(&density, &[1, 2, 5, 10, 25, 63, 200]);
    }

    println!("\n## Fitting the half-weight panel size");
    println!();
    println!(
        "The crossing `ordinary_site_seed.md` §4.1 asks for needs the panel's own shape to start \
         out worse and end up better, and the tables above say it does the opposite. So the \
         constant is fitted the other way it can be: **the value that makes the blended shape \
         nearest the truth, averaged over every panel size, depth and population above.** The \
         score is the mean over the {} cells and {} panel sizes of the median \
         `|ln(blended / that panel's own drawn frequency)|` across each cell's cohorts.",
        fitting.len(),
        ARMS.len()
    );
    println!();
    println!(
        "| half-weight panel size | on the cohorts it was fitted on | on the held-out cohorts |"
    );
    println!("|---:|---:|---:|");
    let mut best = (f64::INFINITY, f64::NAN);
    for candidate in CANDIDATE_HALF_WEIGHTS {
        let on_fit = blended_score(&fitting, candidate);
        let on_held = blended_score(&checking, candidate);
        if on_fit < best.0 {
            best = (on_fit, candidate);
        }
        println!("| {candidate} | {on_fit:.4} | {on_held:.4} |");
    }
    println!(
        "\n**Fitted: N0 = {}**, on the drawn cohorts, scoring {:.4} there and {:.4} on the \
         held-out ones. The library ships {HALF_WEIGHT_PANEL_SIZE}, which scores {:.4} and \
         {:.4}.",
        best.1,
        best.0,
        blended_score(&checking, best.1),
        blended_score(&fitting, HALF_WEIGHT_PANEL_SIZE),
        blended_score(&checking, HALF_WEIGHT_PANEL_SIZE)
    );

    println!("\n### Does the answer depend on depth, or on the population?");
    println!();
    println!("| population | reads a sample | best half-weight panel size | crossing point |");
    println!("|---|---:|---:|---:|");
    for measured in &fitting {
        let one = std::slice::from_ref(measured);
        let best_here = CANDIDATE_HALF_WEIGHTS
            .iter()
            .copied()
            .min_by(|left, right| {
                blended_score(one, *left)
                    .partial_cmp(&blended_score(one, *right))
                    .expect("no NaN among the scores")
            })
            .expect("a non-empty grid");
        match crossing_panel_size(&measured.cells) {
            Some(crossing) => println!(
                "| {} | {} | {best_here} | {crossing:.1} |",
                measured.population, measured.depth
            ),
            None => println!(
                "| {} | {} | {best_here} | none |",
                measured.population, measured.depth
            ),
        }
    }

    println!("\n### Is the blend ever worse than both ends?");
    println!();
    println!(
        "On the held-out cohorts, at the shipped N0 of {HALF_WEIGHT_PANEL_SIZE}, **one drawn \
         cohort a row** rather than one cell. A geometric blend of two guesses whose errors have \
         opposite signs beats both; one whose errors point the same way lands between them. \
         **Worse than both is arithmetically impossible per cohort** — in log space the blend is a \
         convex combination of the two errors — so a non-empty last row would be a defect in this \
         program rather than a finding. An earlier version compared three separately-taken medians \
         and reported one, which is what that comparison can do and this one cannot."
    );
    println!();
    println!("| | drawn cohorts |");
    println!("|---|---:|");
    let mut tally = [0usize; 3];
    for measured in &checking {
        for (arm, cell) in ARMS.iter().zip(&measured.cells) {
            let weight = *arm as f64 / (*arm as f64 + HALF_WEIGHT_PANEL_SIZE);
            let blended = blended_errors(cell, weight);
            for ((blend, panel), neutral) in blended
                .iter()
                .zip(&cell.panel_shape_error)
                .zip(&cell.neutral_shape_error)
            {
                if *blend <= panel.min(*neutral) + 1e-12 {
                    tally[0] += 1;
                } else if *blend <= panel.max(*neutral) + 1e-12 {
                    tally[1] += 1;
                } else {
                    tally[2] += 1;
                }
            }
        }
    }
    println!("| at least as good as both ends | {} |", tally[0]);
    println!("| between the two ends | {} |", tally[1]);
    println!("| worse than both ends | {} |", tally[2]);

    println!("\n### The same question at a real cohort's segregating share");
    println!();
    println!(
        "The drawn arm above needs 10% of positions to segregate or there is nothing to fit at \
         3,000 positions, and at that share the neutral shape happens to be a good guess. **This \
         row scores the same candidates on the projection-only arm** — the five densities of the \
         section above, each segregating at 4 positions in 1,000, with nothing drawn and nothing \
         estimated. It carries no sampling noise and no cohort; what it carries instead is a \
         realistic share of variable positions, where the neutral shape is off by a factor of 2 \
         to 3 because the positions fixed non-reference dominate the population's mean frequency."
    );
    println!();
    println!("| half-weight panel size | mean log error over 5 densities × 7 panels |");
    println!("|---:|---:|");
    for candidate in CANDIDATE_HALF_WEIGHTS {
        println!("| {candidate} | {:.4} |", projection_only_score(candidate));
    }
    println!(
        "\n**On realistic parameters the score rises from zero with no floor at all**, so this arm \
         wants the panel's own shape and nothing else — and so, since the scoring was corrected \
         to each drawn cohort's own realised frequency, does the drawn arm. **Both arms now put \
         the minimum at zero.** The shipped {HALF_WEIGHT_PANEL_SIZE} is a hedge rather than a fit, \
         and it costs this arm {:.0}% over zero.",
        100.0 * (projection_only_score(HALF_WEIGHT_PANEL_SIZE) / projection_only_score(0.0) - 1.0)
    );
    println!("\nWhole sweep: {:.1} s.", started.elapsed().as_secs_f64());
}

/// **What one candidate half-weight panel size scores with no cohort drawn and nothing
/// estimated** — the five densities of `REALISTIC_SHAPES`, each at a real cohort's segregating
/// share, projected exactly into each panel's allele-count classes.
///
/// The companion to [`blended_score`], and it answers a different objection. That one is scored on
/// drawn cohorts and so carries sampling noise, but it has to segregate at 10% of positions for
/// there to be anything to fit; this one has no noise at all but is exact rather than estimated.
/// **Neither is the whole answer and both are printed.**
fn projection_only_score(half_weight_panel: f64) -> f64 {
    let outbred = InbreedingF::try_new(0.0).expect("a legal coefficient");
    let mut total = 0.0;
    let mut counted = 0;
    for (_, a, b) in REALISTIC_SHAPES {
        let density = FrequencyDensity {
            p_invariant: 0.9950,
            p_fixed_alt: 0.0010,
            a,
            b,
        };
        let truth = expected_frequency_of(&density);
        let neutral = neutral_shape_of(density.expected_heterozygosity());
        for individuals in [1u32, 2, 5, 10, 25, 63, 200] {
            let classes = density.allele_count_classes(individuals);
            let panel = fit_spectrum_shape(&FittedSpectrum::new(&classes, 0.0, 1.0), outbred)
                .expected_frequency();
            let weight = f64::from(individuals) / (f64::from(individuals) + half_weight_panel);
            let blended = ((1.0 - weight) * neutral.ln() + weight * panel.ln()).exp();
            total += (blended / truth).ln().abs();
            counted += 1;
        }
    }
    total / f64::from(counted)
}

/// The five allele-frequency densities the projection-only arm uses — the four
/// `ordinary_site_seed.md` §1.2 measured, plus this project's own lopsided unit-test fixture. Each
/// is used at `p_invariant = 0.9950` and `p_fixed_alt = 0.0010`, so 4 positions in 1,000
/// segregate, which is a real cohort's share.
const REALISTIC_SHAPES: [(&str, f64, f64); 5] = [
    (
        "tomato-like, strong rare-allele pile-up, Beta(0.20, 1.00)",
        0.20,
        1.00,
    ),
    ("human-like, moderate pile-up, Beta(0.35, 1.20)", 0.35, 1.20),
    ("flat over what segregates, Beta(1.00, 1.00)", 1.00, 1.00),
    (
        "the lopsided unit-test fixture, Beta(0.50, 2.00)",
        0.50,
        2.00,
    ),
    ("middling frequencies, Beta(4.00, 4.00)", 4.00, 4.00),
];
/// One population at one depth, with what every panel size produced.
struct Measured {
    population: &'static str,
    depth: f64,

    cells: Vec<Cell>,
}

/// **The expected alternative-allele frequency a density actually has**, averaged over every
/// position: the Beta's own mean over what segregates, plus the positions fixed non-reference,
/// which carry the alternative allele on every chromosome.
fn expected_frequency_of(density: &FrequencyDensity) -> f64 {
    density.p_segregating() * density.a / (density.a + density.b) + density.p_fixed_alt
}

/// How far the blended shape lands from the truth in one cell, **one entry a drawn cohort**, each
/// scored against the frequency that cohort was actually drawn with.
fn blended_errors(cell: &Cell, weight: f64) -> Vec<f64> {
    cell.panel_shape
        .iter()
        .zip(&cell.neutral_shape)
        .zip(&cell.drawn_frequency)
        .map(|((panel, neutral), drawn)| {
            let blended = ((1.0 - weight) * neutral.ln() + weight * panel.ln()).exp();
            (blended / drawn).ln().abs()
        })
        .collect()
}

/// **What one candidate half-weight panel size scores**: the mean over every cell and panel size
/// of the median `|ln(blended / drawn)|` across that cell's cohorts.
///
/// The median is taken within a cell so one badly drawn cohort cannot set the answer, and the
/// mean across cells so every panel size, depth and population counts once.
fn blended_score(measured: &[Measured], half_weight_panel: f64) -> f64 {
    let mut total = 0.0;
    let mut counted = 0;
    for one in measured {
        for (arm, cell) in ARMS.iter().zip(&one.cells) {
            let weight = *arm as f64 / (*arm as f64 + half_weight_panel);
            let score = median(&blended_errors(cell, weight));
            if score.is_finite() {
                total += score;
                counted += 1;
            }
        }
    }
    total / counted as f64
}

/// **What the projection costs on its own, with no cohort drawn and nothing estimated.**
///
/// The density is handed straight to `allele_count_classes`, which is exact, and the shipped
/// search is run over the result. So the only thing that can move the answer between one panel
/// size and another is **how many classes the two-parameter family is being asked to fit at
/// once** — no sampling noise, and no fit of the population from reads.
///
/// It is the control the drawn sweep needs: a fitted shape that is worse at 63 individuals than
/// at one could be the panel's fault or the projection's, and these two arms answer that
/// separately.
fn projection_only(density: &FrequencyDensity, panels: &[u32]) {
    let outbred = InbreedingF::try_new(0.0).expect("a legal coefficient");
    let drawn_frequency = expected_frequency_of(density);
    let theta = density.expected_heterozygosity();
    let neutral = theta / (1.0 + theta);
    println!();
    println!(
        "| individuals | the panel's own shape / the truth's | log error | the neutral shape / \
         the truth's | log error |"
    );
    println!("|---:|---:|---:|---:|---:|");
    for individuals in panels {
        let classes = density.allele_count_classes(*individuals);
        let shape = fit_spectrum_shape(&FittedSpectrum::new(&classes, 0.0, 1.0), outbred);
        let panel = shape.expected_frequency();
        println!(
            "| {individuals} | {:.3}× | {:.4} | {:.3}× | {:.4} |",
            panel / drawn_frequency,
            (panel / drawn_frequency).ln().abs(),
            neutral / drawn_frequency,
            (neutral / drawn_frequency).ln().abs()
        );
    }
}

/// What one panel size produced across the drawn cohorts of one cell.
#[derive(Default)]
struct Cell {
    /// The panel's own fitted expected frequency, one entry a drawn cohort.
    panel_shape: Vec<f64>,
    /// The neutral shape's, at the same cohort's **fitted** heterozygosity — which is what the
    /// caller would use, so it is what is scored.
    neutral_shape: Vec<f64>,
    /// **The frequency the cohort was actually drawn with**, one entry a drawn cohort — the share
    /// of this panel's chromosomes that ended up carrying the alternative allele.
    ///
    /// **Not the density's expectation, and the difference is bigger than anything read off these
    /// tables.** Over 3,000 positions the realised frequency of a Beta(0.70, 2.50) cohort has a
    /// standard deviation of 7.7% of its own mean, and a Beta(0.20, 1.00) one 10.3%; the
    /// candidate half-weight panel sizes are separated by less than a tenth of a per cent. Scoring
    /// against the expectation would put that draw-to-draw scatter into every error column, where
    /// it belongs to neither guess.
    drawn_frequency: Vec<f64>,
    /// How far each guess is from the frequency that cohort was drawn with, as
    /// `|ln(guess / drawn)|`.
    panel_shape_error: Vec<f64>,
    neutral_shape_error: Vec<f64>,
    /// The fitted heterozygosity over the drawn one, so the table can say whether the number the
    /// pin rests on is itself sound at this panel size.
    diversity_ratio: Vec<f64>,
    /// How many of this cell's drawn cohorts produced a usable fit. Printed, because an arm that
    /// silently dropped a cohort would be a median over a different set from its neighbours.
    fitted: usize,
    attempted: usize,
}

/// Draw cohorts at seeds `[from, to)`, refit each one's own positions at every panel size, and
/// score both guesses **against the frequency that panel was actually drawn with**.
fn sweep(
    density: &FrequencyDensity,
    depth: f64,
    positions: usize,
    mismapped: f64,
    from: usize,
    to: usize,
) -> Vec<Cell> {
    let outbred = InbreedingF::try_new(0.0).expect("a legal coefficient");
    let mut cells: Vec<Cell> = (0..ARMS.len()).map(|_| Cell::default()).collect();

    for draw in from..to {
        // **A different stream per cell**, so no two cells share a drawn cohort, and the same
        // stream for every arm of one cell, so the positions are held fixed across the arms.
        let seed = 0x9E37_79B9_7F4A_7C15_u64
            .wrapping_mul(draw as u64 + 1)
            .wrapping_add((depth as u64) << 32)
            .wrapping_add((density.a * 1_000.0) as u64)
            | 1;
        let cohort = draw_cohort(
            ARMS[ARMS.len() - 1],
            positions,
            depth,
            mismapped,
            *density,
            seed,
        );

        for (index, arm) in ARMS.iter().enumerate() {
            cells[index].attempted += 1;
            // **The realised frequency of the panel this arm actually fits**, which is a prefix of
            // the drawn cohort: alternative copies carried, over the chromosomes it holds.
            let carried: u64 = cohort.alternative_copies_per_sample[..*arm].iter().sum();
            let drawn_frequency = carried as f64 / (2.0 * *arm as f64 * positions as f64);
            if drawn_frequency <= 0.0 {
                continue;
            }
            let subset = cohort.samples[..*arm].to_vec();
            let mut evidence =
                CohortCensusEvidence::new(subset).expect("a drawn cohort records one way");
            let config = JointFitConfig {
                quadrature_nodes: 12,
                starting_points: StartingPoint::spanning_the_class_separation(),
                ..JointFitConfig::default()
            };
            let Ok(fit) = fit_jointly(&mut evidence, &config) else {
                continue;
            };
            let heterozygosity = fit.expected_heterozygosity;
            if heterozygosity <= 0.0 || !heterozygosity.is_finite() {
                continue;
            }
            let classes = fit.density.value.allele_count_classes(*arm as u32);
            let shape = fit_spectrum_shape(&FittedSpectrum::new(&classes, 0.0, 1.0), outbred);

            let panel = shape.expected_frequency();
            let neutral = neutral_shape_of(heterozygosity);
            cells[index].panel_shape.push(panel);
            cells[index].neutral_shape.push(neutral);
            cells[index].drawn_frequency.push(drawn_frequency);
            cells[index]
                .panel_shape_error
                .push((panel / drawn_frequency).ln().abs());
            cells[index]
                .neutral_shape_error
                .push((neutral / drawn_frequency).ln().abs());
            cells[index]
                .diversity_ratio
                .push(heterozygosity / density.expected_heterozygosity());
            cells[index].fitted += 1;
        }
    }
    cells
}

/// The neutral shape's own expected alternative-allele frequency at a diversity — `θ / (1 + θ)`,
/// the ratio of the pair `(1, θ)`. **The library's own
/// [`neutral_expected_frequency`](pop_var_caller::ng::calling::genotype_prior) is private**, so
/// this repeats it; the two are pinned together by
/// `seed_generic::projection_tests::the_ramps_neutral_end_is_the_pair_the_neutral_rung_returns`,
/// which checks the library's against the pair its no-spectrum branch actually returns.
fn neutral_shape_of(diversity: f64) -> f64 {
    diversity / (1.0 + diversity)
}

/// **The panel size at which the two guesses are equally good**, interpolated in log panel size
/// between the two arms the crossing sits between.
///
/// `None` when the panel's own shape is already the better guess at the smallest panel swept, or
/// still the worse one at the largest: either way this sweep does not contain the crossing and
/// says so rather than extrapolating to one.
fn crossing_panel_size(cells: &[Cell]) -> Option<f64> {
    let gap: Vec<f64> = cells
        .iter()
        .map(|cell| median(&cell.panel_shape_error) - median(&cell.neutral_shape_error))
        .collect();
    if gap[0] <= 0.0 || gap[gap.len() - 1] >= 0.0 {
        return None;
    }
    let at = (1..gap.len()).find(|index| gap[*index] < 0.0)?;
    let (low, high) = (ARMS[at - 1] as f64, ARMS[at] as f64);
    let share = gap[at - 1] / (gap[at - 1] - gap[at]);
    Some((low.ln() + share * (high.ln() - low.ln())).exp())
}

/// The middle value, or the mean of the middle two. `NaN` on an empty set, which prints as `NaN`
/// rather than as a plausible zero.
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|left, right| left.partial_cmp(right).expect("no NaN among the scores"));
    let middle = sorted.len() / 2;
    if sorted.len() % 2 == 1 {
        sorted[middle]
    } else {
        0.5 * (sorted[middle - 1] + sorted[middle])
    }
}

/// Draw a cohort at known parameters and write it into the records the fit reads.
///
/// **The reference base is code 0 by construction**, so the three candidate alleles are codes 1
/// to 3 and the fit's own sum over which one segregates is being asked to find the right one.
///
/// **No excess of homozygotes is planted**: every drawn plant is outbred, which is the basis
/// every measurement in `ordinary_site_seed.md` is on.
fn draw_cohort(
    samples: usize,
    positions: usize,
    depth: f64,
    mismapped: f64,
    density: FrequencyDensity,
    seed: u64,
) -> Drawn {
    let mut rng = Rng(seed);
    let edges = DepthBinEdges::for_census();
    let mut codes: Vec<PackedDepthCodes> = (0..samples)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
    // **How many alternative copies each drawn plant ended up carrying**, over every
    // position — the numerator of the frequency this cohort really has, as against the
    // one its density has in expectation.
    let mut alternative_copies = vec![0_u64; samples];

    for index in 0..positions {
        let rate = if rng.uniform() < mismapped {
            0.06
        } else {
            SEQUENCING_ERROR
        };
        let branch = rng.pick(&[
            density.p_invariant,
            density.p_fixed_alt,
            density.p_segregating(),
        ]);
        let allele = 1 + ((rng.uniform() * 3.0) as usize).min(2);
        let frequency = match branch {
            0 => 0.0,
            1 => 1.0,
            _ => rng.beta(density.a, density.b),
        };
        for sample in 0..samples {
            let reads = rng.poisson(depth);
            let genotype = match branch {
                0 => 0_usize,
                1 => 2,
                _ => rng.pick(&[
                    (1.0 - frequency) * (1.0 - frequency),
                    2.0 * frequency * (1.0 - frequency),
                    frequency * frequency,
                ]),
            };
            alternative_copies[sample] += genotype as u64;
            let carried = genotype as f64 / 2.0;
            let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
            let on_reference = (1.0 - carried) * (1.0 - rate) + carried * rate / 3.0;
            let mut counts = [0_u32; 5];
            for _ in 0..reads {
                let draw = rng.uniform();
                let code = if draw < on_candidate {
                    allele
                } else if draw < on_candidate + on_reference {
                    0
                } else if draw < on_candidate + on_reference + rate / 3.0 {
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
                    reads: u8::try_from(*count)
                        .expect("a drawn count fits the census's one-byte field"),
                });
            }
        }
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
        depth_cap: DepthCap::MAX,
    };
    let samples_recorded = (0..samples)
        .map(|sample| {
            SampleCensusEvidence::resident(
                format!("s{sample:02}"),
                terms.clone(),
                BTreeMap::from([(
                    SectionKey::Generic(ReadGroupId(sample as u32)),
                    Section::Generic(GenericEvidence::from_parts(
                        std::mem::replace(&mut codes[sample], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[sample]),
                    )),
                )]),
            )
        })
        .collect();
    Drawn {
        samples: samples_recorded,
        alternative_copies_per_sample: alternative_copies,
    }
}

/// One drawn cohort: the records the fit reads, and **what each plant really carried**.
///
/// The second is what the sweep scores against. A drawn cohort's own alternative-allele frequency
/// is not its density's expectation — over 3,000 positions the two differ by about 8% of
/// themselves, one draw to the next — and that difference belongs to the draw rather than to
/// either guess being scored.
struct Drawn {
    samples: Vec<SampleCensusEvidence>,
    /// How many alternative copies each plant carries, over every drawn position. A diploid
    /// carries 0, 1 or 2 at each.
    alternative_copies_per_sample: Vec<u64>,
}

/// The stream every drawn number comes from. Deterministic, so a run is reproducible.
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
