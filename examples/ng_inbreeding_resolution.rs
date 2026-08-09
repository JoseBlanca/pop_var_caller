//! **What is the runs model's noise floor a function of, and where do its two states stop
//! being tellable apart?** — the harness for the two numbers Milestone E3 shipped with a
//! constant it could not justify.
//!
//! `resolution_at` answers *what `F` comes back as on a genome with no runs at all*, and a
//! consumer reads a fitted `F` against it to decide whether anything was detected. Research
//! note §3.6 measured it at four **window counts** and nothing else, so the shipped function
//! is a function of the window count alone. That is not what the floor depends on. A
//! two-state chain classifies a window on how many heterozygotes it holds, and that is
//! `sites × heterozygosity`, not the number of windows — §3.6's genomes carried 100,000
//! sites a window at one heterozygote per kilobase, which is **100 heterozygotes a window**,
//! and a run over a region subset or a shallow sample carries far fewer.
//!
//! The same sweep answers a second question. `MAX_IDENTIFIED_STATE_RATIO` refuses a fit
//! whose two states came out too close to tell apart — at coincident states the likelihood
//! is *exactly* flat in `F` (research note §3.1), so what a fit returns is an accident. E3
//! set it at 0.9 from ten fitted ratios on one fixture shape. Whether 0.9 is right, or even
//! whether one constant can serve every shape, needs the ratio measured **on both sides**:
//! on genomes with no runs, where the states genuinely coincide and every fit should be
//! refused; and on genomes with real runs, where they do not and none should be.
//!
//! **This harness simulates rather than computing a bias exactly**, as
//! `ng_inbreeding_harness.rs` does and for the same reason: the runs model's objective is a
//! chain over windows, not a sum over independent cells, so there is no infinite-genome
//! table to weight. Every number below is a mean over seeds with its worst case beside it.
//!
//! ```text
//! cargo run --release --example ng_inbreeding_resolution
//! ```

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::sync::Arc;

use pop_var_caller::ng::parameter_estimation::generic::accumulators::WindowKey;
use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::generic::histogram::{
    DepthAltHistogram, DepthAndAltReads,
};
use pop_var_caller::ng::parameter_estimation::generic::noise_model::SampleLibraryNoise;
use pop_var_caller::ng::parameter_estimation::generic::runs::{
    MAX_IDENTIFIED_STATE_RATIO, RunsModelStarts, fit_inbreeding, resolution_at,
};
use pop_var_caller::ng::parameter_estimation::generic::{SampleRates, WindowIndex};
use pop_var_caller::ng::types::{
    Bp, ContigId, ErrorRate, GenotypeFrequency, InbreedingF, Ploidy, ReadGroupId,
};

const GENOTYPES: usize = 3;
const ERROR_RATE: f64 = 0.001;
/// The homozygous-non-reference rate, the same in both states — a landrace far from the
/// reference, as `ng_inbreeding_harness.rs`'s scenarios have it.
const HOM_ALT: f64 = 0.02;
/// One window in thirty ends a run, so a run is 3 Mb at 100 kb windows.
const LEAVE_RUN: f64 = 1.0 / 30.0;
/// The inside state's heterozygote rate as a fraction of the outside one — the real
/// separation a drawn genome has, before any fit looks at it.
const TRUE_SEPARATION: f64 = 0.04;
/// How many seeds each cell of the grid is averaged over.
const SEEDS: u64 = 5;

/// SplitMix64 — `ng_inbreeding_harness.rs`'s generator, so a genome drawn here and one
/// drawn there are drawn the same way.
struct SplitMix64(u64);

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn binomial(&mut self, n: u64, p: f64) -> u64 {
        if p <= 0.0 {
            return 0;
        }
        if p >= 1.0 {
            return n;
        }
        (0..n).filter(|_| self.unit() < p).count() as u64
    }
}

/// `p_j(ε)` — the chance one read shows something other than the reference base.
fn alternative_read_probability(alt_copies: u8, error_rate: f64) -> f64 {
    let carried = f64::from(alt_copies) / 2.0;
    carried * (1.0 - error_rate / 3.0) + (1.0 - carried) * error_rate
}

/// One point of the grid: how big the genome is and how much evidence each window carries.
#[derive(Copy, Clone)]
struct Shape {
    windows: usize,
    sites_per_window: u64,
    depth: u32,
    outside_het: f64,
}

impl Shape {
    /// **The quantity the floor should actually depend on**: how many heterozygous sites a
    /// window holds. A two-state chain decides a window on that count, and its sampling
    /// noise is `√h` on a mean of `h`, so what a window can resolve goes as `1/√h`.
    fn heterozygotes_per_window(self) -> f64 {
        self.sites_per_window as f64 * self.outside_het
    }
}

/// A drawn genome, and the autozygous fraction it actually realised.
struct Drawn {
    windowed: BTreeMap<(WindowKey, Ploidy), DepthAltHistogram<u32>>,
    realised_f: f64,
}

fn draw(shape: Shape, nominal_f: f64, seed: u64) -> Drawn {
    let edges = Arc::new(DepthBinEdges::new());
    let diploid = Ploidy::try_new(2).expect("a positive copy number");
    let mut rng = SplitMix64::new(seed);
    let enter = if nominal_f <= 0.0 {
        0.0
    } else {
        LEAVE_RUN * nominal_f / (1.0 - nominal_f)
    };

    let frequencies = |inside: bool| -> [f64; GENOTYPES] {
        let het = if inside {
            shape.outside_het * TRUE_SEPARATION
        } else {
            shape.outside_het
        };
        [1.0 - het - HOM_ALT, het, HOM_ALT]
    };
    let cell_probabilities = |inside: bool| -> Vec<f64> {
        let frequencies = frequencies(inside);
        (0..=shape.depth)
            .map(|alt_reads| {
                let mut total = 0.0;
                for (dosage, &frequency) in frequencies.iter().enumerate() {
                    let p = alternative_read_probability(dosage as u8, ERROR_RATE);
                    let mut term = (1.0 - p).powi(shape.depth as i32);
                    for step in 1..=alt_reads {
                        term *= f64::from(shape.depth - step + 1) / f64::from(step) * p / (1.0 - p);
                    }
                    total += frequency * term;
                }
                total
            })
            .collect()
    };
    let inside_cells = cell_probabilities(true);
    let outside_cells = cell_probabilities(false);

    // Twelve contigs, so the chain restarts eleven times, as a real genome makes it.
    let contigs = 12u32;
    let per_contig = (shape.windows as u32).div_ceil(contigs);
    let mut windowed = BTreeMap::new();
    let mut inside_positions = 0u64;
    let mut total_positions = 0u64;

    for contig in 0..contigs {
        let mut inside = enter > 0.0 && rng.unit() < enter / (enter + LEAVE_RUN);
        for window in 0..per_contig {
            let cells = if inside {
                &inside_cells
            } else {
                &outside_cells
            };
            let mut table = DepthAltHistogram::new(Arc::clone(&edges));

            let mut left = shape.sites_per_window;
            let mut remaining: f64 = cells.iter().sum();
            for (alt_reads, &probability) in cells.iter().enumerate() {
                if left == 0 {
                    break;
                }
                let in_cell = if remaining > 0.0 {
                    rng.binomial(left, (probability / remaining).clamp(0.0, 1.0))
                } else {
                    0
                };
                remaining -= probability;
                left -= in_cell;
                for _ in 0..in_cell {
                    table.add_site(DepthAndAltReads::new(shape.depth, alt_reads as u32), Bp(1));
                }
            }

            total_positions += table.total_covered_positions();
            if inside {
                inside_positions += table.total_covered_positions();
            }
            windowed.insert(
                (
                    WindowKey::InWindow {
                        contig: ContigId(contig),
                        window: WindowIndex(window),
                    },
                    diploid,
                ),
                table,
            );

            inside = if inside {
                rng.unit() >= LEAVE_RUN
            } else {
                enter > 0.0 && rng.unit() < enter
            };
        }
    }

    Drawn {
        windowed,
        realised_f: inside_positions as f64 / total_positions.max(1) as f64,
    }
}

/// What one fit returned, reduced to the two numbers this harness is about.
struct Outcome {
    realised_f: f64,
    /// `None` when the fit refused — which is itself a result, not a gap.
    fitted_f: Option<f64>,
    /// The fitted inside heterozygote rate as a fraction of the outside one, whether or not
    /// the fit was accepted. **Reported even for a refused fit**, because the question is
    /// where the two populations of ratios sit relative to each other.
    state_ratio: Option<f64>,
    /// **How far apart the nine starting points landed**, largest `F` minus smallest.
    ///
    /// Spec §6.5 names this as the one thing separating a genuine answer from a search that
    /// found nothing — "a run where every start returned the same `F` at the same score has
    /// not measured zero autozygosity; it has failed to find anything" — and E3 reports it
    /// in `starts_tried` without ever reading it. If the failures disagree across starts
    /// and the real fits agree, it is the criterion the ratio check could not be.
    across_starts: Option<f64>,
}

fn fit(shape: Shape, nominal_f: f64, seed: u64) -> Outcome {
    let diploid = Ploidy::try_new(2).expect("a positive copy number");
    let drawn = draw(shape, nominal_f, seed);
    let noise = SampleLibraryNoise::single(
        ReadGroupId(1),
        ErrorRate::try_new(ERROR_RATE).expect("a probability"),
    );
    let rates = SampleRates::try_new(
        diploid,
        [
            1.0 - shape.outside_het - HOM_ALT,
            shape.outside_het,
            HOM_ALT,
        ]
        .into_iter()
        .map(|f| GenotypeFrequency::try_new(f).expect("a frequency"))
        .collect(),
    )
    .expect("a well-formed diploid set");

    match fit_inbreeding(
        "sweep",
        &drawn.windowed,
        &noise,
        diploid,
        &rates,
        &RunsModelStarts::default(),
    ) {
        Ok((estimate, model)) => {
            let spread = model
                .starts_tried
                .iter()
                .map(|start| start.inbreeding)
                .fold(f64::NEG_INFINITY, f64::max)
                - model
                    .starts_tried
                    .iter()
                    .map(|start| start.inbreeding)
                    .fold(f64::INFINITY, f64::min);
            Outcome {
                realised_f: drawn.realised_f,
                fitted_f: Some(estimate.value.get()),
                state_ratio: Some(model.inside_het_floor / model.outside_het),
                across_starts: Some(spread),
            }
        }
        Err(_) => Outcome {
            realised_f: drawn.realised_f,
            fitted_f: None,
            state_ratio: None,
            across_starts: None,
        },
    }
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

fn worst(values: &[f64]) -> f64 {
    values.iter().copied().fold(0.0f64, f64::max)
}

/// The grid: three genome sizes crossed with four levels of evidence per window.
///
/// The evidence levels are chosen to bracket what the design cares about. **100
/// heterozygotes a window is research note §3.6's own genome** — 100,000 sites at one
/// heterozygote per kilobase, which is also tomato at 100 kb windows. Twenty is Milestone
/// E3's unit-test fixture. Five and two are what a shallow sample or a region-restricted
/// run leaves, and they are where the shipped `resolution_at` has never been checked.
fn grid() -> Vec<Shape> {
    let mut shapes = Vec::new();
    for windows in [3_000usize, 6_000, 12_000] {
        for &(sites_per_window, outside_het) in &[
            (200u64, 0.01f64), // 2 heterozygotes a window
            (200, 0.025),      // 5
            (400, 0.05),       // 20
            (2_000, 0.05),     // 100 — §3.6's own evidence level
        ] {
            shapes.push(Shape {
                windows,
                sites_per_window,
                depth: 10,
                outside_het,
            });
        }
    }
    shapes
}

fn main() {
    let mut report = String::new();
    let _ = writeln!(
        report,
        "# The runs model's noise floor, and where its two states stop being tellable apart\n"
    );
    let _ = writeln!(
        report,
        "Generated by `cargo run --release --example ng_inbreeding_resolution`. Every row is \
         {SEEDS} seeds.\n"
    );

    // ------------------------------------------------------------------
    let _ = writeln!(
        report,
        "\n## 1. The floor on a genome with **no runs at all**\n"
    );
    let _ = writeln!(
        report,
        "`F` must come back at nothing. What it actually comes back at is the estimator's \
         resolution, and the question is what that is a function of. `shipped` is what \
         `resolution_at(windows)` reports today — a function of the window count alone.\n"
    );
    let _ = writeln!(
        report,
        "{:>8} {:>8} {:>7} {:>8} {:>8} {:>9} {:>13} {:>11}",
        "windows",
        "het/win",
        "refused",
        "mean F",
        "worst F",
        "shipped",
        "worst/shipped",
        "starts spread"
    );

    let mut floor_rows: Vec<(Shape, usize, f64, f64)> = Vec::new();
    for shape in grid() {
        let outcomes: Vec<Outcome> = (0..SEEDS)
            .map(|seed| fit(shape, 0.0, 0xF100_0000 + seed * 7919 + shape.windows as u64))
            .collect();
        let answered: Vec<f64> = outcomes.iter().filter_map(|o| o.fitted_f).collect();
        let refused = outcomes.len() - answered.len();
        let shipped = resolution_at(shape.windows);
        let (m, w) = (mean(&answered), worst(&answered));
        let spreads: Vec<f64> = outcomes.iter().filter_map(|o| o.across_starts).collect();
        let _ = writeln!(
            report,
            "{:>8} {:>8.0} {:>7} {:>8.4} {:>8.4} {:>9.4} {:>13.2} {:>11.4}",
            shape.windows,
            shape.heterozygotes_per_window(),
            refused,
            m,
            w,
            shipped,
            w / shipped,
            worst(&spreads)
        );
        floor_rows.push((shape, refused, m, w));
    }

    // ------------------------------------------------------------------
    let _ = writeln!(
        report,
        "\n## 2. The two states' fitted separation, on both sides\n"
    );
    let _ = writeln!(
        report,
        "The ratio is the fitted inside heterozygote rate over the fitted outside one, \
         reported whether or not the fit was accepted. A genome with no runs has no second \
         state to find, so its ratio should sit near one; a genome 30% covered by runs \
         has one, drawn at {TRUE_SEPARATION}. **The threshold has to separate the two \
         columns at every shape, or one constant cannot serve them all.** Shipped: \
         {MAX_IDENTIFIED_STATE_RATIO}.\n"
    );
    let _ = writeln!(
        report,
        "{:>8} {:>8} {:>21} {:>21} {:>8}",
        "windows", "het/win", "no runs: mean (worst)", "with runs: mean (worst)", "gap"
    );

    for shape in grid() {
        let ratios = |nominal: f64, salt: u64| -> Vec<f64> {
            (0..SEEDS)
                .filter_map(|seed| {
                    fit(shape, nominal, salt + seed * 7919 + shape.windows as u64).state_ratio
                })
                .collect()
        };
        // A refused fit still reports its ratio, so both columns are over all seeds.
        let none = ratios(0.0, 0xF100_0000);
        let some = ratios(0.30, 0xF200_0000);
        let (none_mean, none_worst) = (mean(&none), worst(&none));
        let (some_mean, some_worst) = (mean(&some), worst(&some));
        let _ = writeln!(
            report,
            "{:>8} {:>8.0} {:>12.3} ({:.3}) {:>12.3} ({:.3}) {:>8.3}",
            shape.windows,
            shape.heterozygotes_per_window(),
            none_mean,
            none_worst,
            some_mean,
            some_worst,
            none_mean - some_worst
        );
    }

    // ------------------------------------------------------------------
    let _ = writeln!(
        report,
        "\n## 3. Does a genome that **has** runs still get its answer?\n"
    );
    let _ = writeln!(
        report,
        "The floor is only half the question: a threshold tight enough to refuse every \
         unidentified fit must still answer every identified one. Drawn 30% covered by runs.\n"
    );
    let _ = writeln!(
        report,
        "{:>8} {:>8} {:>7} {:>9} {:>9} {:>9} {:>13}",
        "windows", "het/win", "refused", "realised", "mean F", "worst err", "starts spread"
    );
    for shape in grid() {
        let outcomes: Vec<Outcome> = (0..SEEDS)
            .map(|seed| {
                fit(
                    shape,
                    0.30,
                    0xF200_0000 + seed * 7919 + shape.windows as u64,
                )
            })
            .collect();
        let refused = outcomes.iter().filter(|o| o.fitted_f.is_none()).count();
        let realised = mean(&outcomes.iter().map(|o| o.realised_f).collect::<Vec<_>>());
        let answered: Vec<f64> = outcomes.iter().filter_map(|o| o.fitted_f).collect();
        let errors: Vec<f64> = outcomes
            .iter()
            .filter_map(|o| o.fitted_f.map(|f| (f - o.realised_f).abs()))
            .collect();
        let spreads: Vec<f64> = outcomes.iter().filter_map(|o| o.across_starts).collect();
        let _ = writeln!(
            report,
            "{:>8} {:>8.0} {:>7} {:>9.4} {:>9.4} {:>9.4} {:>13.4}",
            shape.windows,
            shape.heterozygotes_per_window(),
            refused,
            realised,
            mean(&answered),
            worst(&errors),
            worst(&spreads)
        );
    }

    // A sanity line so a reader can tell a run that measured nothing from one that did.
    let _ = writeln!(
        report,
        "\n{} grid points × {SEEDS} seeds × 3 questions; `InbreedingF` accepts \
         [0, 1] and every fitted value above was inside it.",
        grid().len()
    );
    let _ = InbreedingF::try_new(0.5).expect("a fraction");

    print!("{report}");
}
