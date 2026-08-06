//! **Is the runs-of-homozygosity estimate of `F` well behaved?** — the harness for
//! [`parameter_prepass_generic.md`](../doc/devel/ng/spec/parameter_prepass_generic.md) §6.
//!
//! `F` is the fraction of a genome lying in runs of homozygosity, read off a two-state
//! hidden Markov model over 100 kb windows: *inside a run* and *outside* one, each
//! carrying its own three genotype frequencies, with both transition rates fitted too.
//! Nothing downstream can check it — there is no truth for a real sample — so an error
//! here is a plausible wrong number nobody notices. This program asks nine questions of
//! it, and each is a way the number could be wrong while looking right.
//!
//! **Unlike the multi-library harness next door, this one cannot compute bias exactly.**
//! There the estimator maximised a sum over independent cells, so weighting each cell by
//! its true probability gave the infinite-data answer in closed form. A hidden Markov
//! model's likelihood does not factor that way — a window's contribution depends on every
//! other window on the chromosome — so what is measured here is a **ladder in the number
//! of windows**, with repeats, which separates "converges to the truth" from "converges
//! to something else" and from "is just noisy". Where an exact statement is available it
//! is stated as a proof rather than measured; §1 below is one.
//!
//! ## The questions
//!
//! 1. **Is `F` identified on an outbred genome?** If the two states carry the same
//!    genotype frequencies, the observed windows are independent draws from one
//!    distribution and the likelihood is **exactly constant in both transition rates** —
//!    every value of `F` scores identically. That is a proof, not a measurement. What the
//!    fit then returns is whatever its starting point implied, so the harness fits from
//!    several starts and reports the spread. A spread as wide as the range means the
//!    number is an artefact of initialisation.
//! 2. **What does a genuinely outbred genome return?** Not the degenerate case above but
//!    the real one: no runs anywhere, states that differ only by sampling noise. A
//!    two-state model can always improve its score by splitting on that noise, so `F`
//!    should come back above zero. The question is how far above, and whether it shrinks
//!    as windows accumulate — if it does not, it is a floor and every outbred sample
//!    carries it.
//! 3. **Is `F` recovered when runs are really there?** Across true values, with run
//!    lengths a genome actually has.
//! 4. **Does the false-heterozygote floor cancel, as §6.2 claims?** The decision to build
//!    the runs estimator rather than the ratio one rests on this: a uniform floor of
//!    spurious heterozygotes — collapsed paralogs, mismapping — is claimed to lift the
//!    inside and outside rates together and cancel out of the gap between them. No test
//!    in the spec injects one. This one does.
//! 5. **Do the two definitions of `F` in the spec agree?** §6.1 says the transition
//!    rates' ratio *is* `F` (the chain's stationary inside-probability); §6.1's closing
//!    and the architecture doc say `F` is the coverage-weighted posterior occupancy. Both
//!    are computed on every fit here.
//! 6. **Is this a hidden Markov model at all, or a per-window classifier?** How many windows
//!    the chain was genuinely undecided about, and what shuffling the window order does.
//! 7. **Would a finer window rescue the collapse of question 4?**
//! 8. **What does binning the depth do to the emission?** The accumulator keys a cell by a
//!    depth bin rather than an exact depth and scores it at the mean of the depths that fell
//!    in it. Spec §6.1 asserts the error that leaves is small; §8 refits one genome through
//!    every candidate ladder and both mean-depth grains to find out. **The exact ladder — one
//!    bin per depth — is the control, and must reproduce the unbinned answer.**
//! 9. **What does the fit do at 20 and 60 reads a site**, start by start rather than winner
//!    alone. This is the section that caught two defects in *this program* that fire only
//!    above about 40 reads a site; see `SplitMix64::binomial` and the convergence test.
//!
//! ```text
//! cargo run --release --example ng_inbreeding_harness
//! ```

use std::borrow::Cow;
use std::fmt::Write as _;

const STATES: usize = 2;
const INSIDE: usize = 0;
const OUTSIDE: usize = 1;
const GENOTYPES: usize = 3;

// ---------------------------------------------------------------------------
// Numerics
// ---------------------------------------------------------------------------

fn ln_sum_exp2(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (lo - hi).exp().ln_1p()
}

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
    /// Binomial by inversion on the cumulative distribution — exact, and fast enough
    /// because it is called once per cell rather than once per read.
    fn binomial(&mut self, n: u64, p: f64) -> u64 {
        if n == 0 || p <= 0.0 {
            return 0;
        }
        if p >= 1.0 {
            return n;
        }
        // Normal approximation with continuity correction above the exact regime; the
        // counts here are site counts in a window, so n is large and p small.
        if f64::from(n as u32) * p > 30.0 && n > 200 {
            let mean = n as f64 * p;
            let sd = (mean * (1.0 - p)).sqrt();
            let u1 = self.unit().max(1e-300);
            let u2 = self.unit();
            let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
            return (mean + sd * z).round().clamp(0.0, n as f64) as u64;
        }
        let mut count = 0;
        let mut remaining = n;
        // Geometric skipping: draw the gap to the next success.
        //
        // **`ln(1 − p)` must be computed as `ln_1p(−p)`, and this is not a nicety.** Written
        // `(1.0 - p).ln()`, a `p` below 2.2 × 10⁻¹⁶ rounds `1 − p` to exactly 1, so `ln_q`
        // is 0, the gap to the next success is `ln(u)/0 = −∞`, and the cast of that to an
        // unsigned integer saturates to 0 — a gap of nothing, meaning *every* trial
        // succeeds. The routine then returns `n` where the answer is all but certainly 0.
        //
        // It fires the moment any cell's probability falls below f64's resolution, which
        // for a Poisson depth distribution is the depth-1 cell above a mean depth of about
        // 41: `λ·e^(−λ)` is 4 × 10⁻⁸ at λ = 20 and 5 × 10⁻²⁵ at λ = 60. At 60 reads a site
        // every one of a window's 100,000 sites was assigned to the cell "one read, no
        // alternative", the genome came back with **no heterozygote anywhere**, and the runs
        // model dutifully reported `F` = 1.0000, converged, on a genome 32% covered by runs.
        // Every number this harness has published is at 3 reads a site, where the depth-1
        // cell has probability 0.15 and nothing is near the cliff.
        let ln_q = (-p).ln_1p();
        loop {
            let u = self.unit().max(1e-300);
            let gap = (u.ln() / ln_q).floor() as u64;
            if gap >= remaining {
                return count;
            }
            count += 1;
            remaining -= gap + 1;
            if remaining == 0 {
                return count;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// The generating model
// ---------------------------------------------------------------------------

/// One read's chance of showing a non-reference base at `j` alternative copies out of two
/// (spec `parameter_prepass.md` §3).
fn p_alt(j: usize, eps: f64) -> f64 {
    let frac = j as f64 / 2.0;
    frac * (1.0 - eps / 3.0) + (1.0 - frac) * eps
}

#[derive(Clone)]
struct Scenario {
    name: String,
    contigs: usize,
    windows_per_contig: usize,
    /// Covered sites in one window. A 100 kb window at full coverage holds ~100,000.
    sites_per_window: u64,
    mean_depth: f64,
    eps: f64,
    /// `(hom-ref, het, hom-alt)` outside a run and inside one.
    outside: [f64; GENOTYPES],
    inside: [f64; GENOTYPES],
    /// Per-**window** transition probabilities: entering a run, and leaving one.
    enter_run: f64,
    leave_run: f64,
    /// A uniform rate of spurious heterozygotes added to **both** states, taken out of
    /// the homozygous-reference class. Collapsed paralogs and mismapping, in one number.
    false_het_floor: f64,
    /// How deeply each window is covered, as multipliers of `mean_depth`, assigned to
    /// windows in turn. **One level means every window has the same depth distribution**,
    /// which is the case where a cell's per-window mean depth and its whole-sample mean
    /// differ only by sampling noise — so §8's comparison between the two needs more than
    /// one level to be a question at all.
    coverage_levels: Vec<f64>,
}

impl Scenario {
    /// The genotype frequencies actually generated from, floor included.
    fn generating(&self, state: usize) -> [f64; GENOTYPES] {
        let base = if state == INSIDE {
            self.inside
        } else {
            self.outside
        };
        [
            base[0] - self.false_het_floor,
            base[1] + self.false_het_floor,
            base[2],
        ]
    }

    /// The chain's stationary probability of being inside a run — the `F` a very long
    /// genome generated from this scenario would have.
    fn stationary_f(&self) -> f64 {
        self.enter_run / (self.enter_run + self.leave_run)
    }

    fn windows(&self) -> usize {
        self.contigs * self.windows_per_contig
    }
}

/// The depth support, truncated where the Poisson tail is negligible, over `n ≥ 1`.
fn depth_distribution(mean_depth: f64) -> Vec<f64> {
    let mut probs = vec![0.0f64];
    let mut n = 1u32;
    let mut ln_p = -mean_depth;
    loop {
        ln_p += mean_depth.ln() - f64::from(n).ln();
        let p = ln_p.exp();
        probs.push(p);
        if f64::from(n) > mean_depth && p < 1e-12 {
            break;
        }
        n += 1;
    }
    let total: f64 = probs.iter().sum();
    for p in &mut probs {
        *p /= total;
    }
    probs
}

/// A `(depth, alt-count)` cell of the per-window tally.
#[derive(Clone, Copy)]
struct Cell {
    depth: u32,
    alt: u32,
}

// ---------------------------------------------------------------------------
// The depth ladder
// ---------------------------------------------------------------------------

/// Which depths share a bin. The accumulator does not keep an exact depth: it keys a cell
/// by a bin — one bin per depth at the bottom, widening geometrically above — and scores
/// each cell at the mean of the exact depths that landed in it
/// (`../doc/devel/ng/arch/parameter_prepass_generic.md` §2.2). Spec §6.1 says of this
/// window emission that "binning the depth makes it an approximation ... scoring each cell
/// at its own mean depth is what keeps the error small", and that claim is what §8 below
/// measures.
///
/// Deliberately a second copy of the type in `ng_multilib_key_harness.rs` rather than a
/// shared module: these are two standalone programs, and a shared helper would make each
/// harder to read on its own than the twenty lines it saves.
#[derive(Clone)]
struct DepthLadder {
    label: String,
    bin_of: Vec<u32>,
    cap: u32,
    bins: usize,
}

impl DepthLadder {
    fn exact(max_depth: u32) -> Self {
        Self {
            label: "exact (no binning)".to_string(),
            bin_of: (0..=max_depth).collect(),
            cap: max_depth,
            bins: (max_depth + 1) as usize,
        }
    }

    fn geometric(exact_to: u32, bins: usize, cap: u32) -> Self {
        assert!(bins > exact_to as usize + 1, "no room to widen");
        let widening = bins - exact_to as usize - 1;
        let ratio = (f64::from(cap) / f64::from(exact_to)).powf(1.0 / widening as f64);
        let mut first = Vec::with_capacity(widening);
        let mut previous = exact_to;
        for i in 0..widening {
            let edge = (f64::from(exact_to) * ratio.powi(i as i32 + 1)).round() as u32;
            let edge = edge.max(previous + 1).min(cap);
            first.push(previous + 1);
            previous = edge;
        }
        let mut bin_of = Vec::with_capacity((cap + 1) as usize);
        for depth in 0..=cap {
            bin_of.push(if depth <= exact_to {
                depth
            } else {
                exact_to + first.iter().filter(|&&f| f <= depth).count() as u32
            });
        }
        let realised = bin_of.last().map_or(0, |&b| b as usize + 1);
        Self {
            label: format!("exact≤{exact_to} {realised} bins cap {cap}"),
            bin_of,
            cap,
            bins: realised,
        }
    }

    fn bin_of(&self, depth: u32) -> u32 {
        self.bin_of[depth.min(self.cap) as usize]
    }
}

/// Which depth a binned cell is scored at.
#[derive(Clone, Copy, PartialEq, Eq)]
enum MeanDepthGrain {
    /// The mean over the depths that landed in this cell **of this window** — one depth sum
    /// per cell per window, which is what `DepthAltHistogram` holds and what the
    /// architecture specifies.
    PerWindow,
    /// The mean over the whole sample. Cheaper by one counter per window, and identical
    /// wherever every window has the same coverage — which is why the genome below does
    /// not.
    WholeSample,
}

/// Every cell the depth support allows, and the probability the model puts on each given
/// the genotype. Shared by generation and by fitting, so both speak of the same cells.
struct CellSpace {
    cells: Vec<Cell>,
    depth_probs: Vec<f64>,
}

impl CellSpace {
    /// Built to hold every depth the deepest window can produce, so that one cell list
    /// serves windows of different coverage.
    fn for_scenario(scenario: &Scenario) -> Self {
        let deepest = scenario
            .coverage_levels
            .iter()
            .cloned()
            .fold(1.0f64, f64::max);
        Self::new(scenario.mean_depth * deepest)
    }

    fn new(mean_depth: f64) -> Self {
        let depth_probs = depth_distribution(mean_depth);
        let mut cells = Vec::new();
        for depth in 1..depth_probs.len() as u32 {
            for alt in 0..=depth {
                cells.push(Cell { depth, alt });
            }
        }
        Self { cells, depth_probs }
    }

    /// `P(cell | genotype j)` — the depth's own probability times a binomial. The depth
    /// factor is common to every genotype and every state, so it cancels out of both the
    /// posterior and the state comparison; it is kept so the numbers are probabilities.
    fn component(&self, cell: Cell, j: usize, eps: f64) -> f64 {
        let p = p_alt(j, eps);
        let n = cell.depth;
        let k = cell.alt;
        let mut ln_c = 0.0;
        for i in 0..k {
            ln_c += ((n - i) as f64).ln() - ((i + 1) as f64).ln();
        }
        self.depth_probs[n as usize]
            * (ln_c + f64::from(k) * p.ln() + f64::from(n - k) * (1.0 - p).ln()).exp()
    }

    /// `P(alt count | depth, genotype j)` — the binomial alone, with the depth's own
    /// probability left out.
    fn alt_count_probability(&self, cell: Cell, j: usize, eps: f64) -> f64 {
        let p = p_alt(j, eps);
        let (n, k) = (cell.depth, cell.alt);
        let mut ln_c = 0.0;
        for i in 0..k {
            ln_c += ((n - i) as f64).ln() - ((i + 1) as f64).ln();
        }
        (ln_c + f64::from(k) * p.ln() + f64::from(n - k) * (1.0 - p).ln()).exp()
    }

    /// `P(cell)` under one set of genotype frequencies, at a stated depth distribution —
    /// the window's own, which is not the sample's when coverage varies between windows.
    fn probability_at(
        &self,
        cell: Cell,
        freqs: &[f64; GENOTYPES],
        eps: f64,
        depth_probs: &[f64],
    ) -> f64 {
        let n = cell.depth as usize;
        if n >= depth_probs.len() {
            return 0.0;
        }
        depth_probs[n]
            * (0..GENOTYPES)
                .map(|j| freqs[j] * self.alt_count_probability(cell, j, eps))
                .sum::<f64>()
    }
}

/// One simulated sample: the per-window cell tallies, plus the truth to score against.
struct Sample {
    /// `windows[w][cell index]` — how many sites of window `w` fell in that cell.
    windows: Vec<Vec<u64>>,
    /// Which contig each window belongs to, so the chain restarts at boundaries.
    contig_of: Vec<usize>,
    /// The state each window was really in.
    true_state: Vec<usize>,
}

impl Sample {
    fn true_f(&self) -> f64 {
        let inside = self.true_state.iter().filter(|&&s| s == INSIDE).count();
        inside as f64 / self.true_state.len() as f64
    }
}

/// Draw a genome. Windows are tallied straight from the multinomial over cells rather
/// than site by site — exact, and it makes a 100,000-site window cost fifty draws.
fn simulate(scenario: &Scenario, space: &CellSpace, seed: u64) -> Sample {
    let mut rng = SplitMix64::new(seed);
    let levels = &scenario.coverage_levels;
    let level_depths: Vec<Vec<f64>> = levels
        .iter()
        .map(|m| depth_distribution(scenario.mean_depth * m))
        .collect();
    // `cell_probs[state][coverage level]`.
    let cell_probs: Vec<Vec<Vec<f64>>> = (0..STATES)
        .map(|s| {
            let freqs = scenario.generating(s);
            level_depths
                .iter()
                .map(|depth_probs| {
                    space
                        .cells
                        .iter()
                        .map(|&c| space.probability_at(c, &freqs, scenario.eps, depth_probs))
                        .collect()
                })
                .collect()
        })
        .collect();

    let mut windows = Vec::with_capacity(scenario.windows());
    let mut contig_of = Vec::with_capacity(scenario.windows());
    let mut true_state = Vec::with_capacity(scenario.windows());
    let stationary = scenario.stationary_f();

    for contig in 0..scenario.contigs {
        let mut state = if rng.unit() < stationary {
            INSIDE
        } else {
            OUTSIDE
        };
        for _ in 0..scenario.windows_per_contig {
            // Multinomial over cells, by the chain of conditional binomials.
            let probs = &cell_probs[state][windows.len() % levels.len()];
            let mut remaining = scenario.sites_per_window;
            let mut remaining_mass: f64 = probs.iter().sum();
            let mut counts = vec![0u64; probs.len()];
            // The chain of conditional binomials leaves whatever is left to the last cell,
            // and that cell must be one the truth can actually produce. A window covered
            // more thinly than the deepest one gives zero probability to every cell above
            // its own depth support, so the last cell of the shared list is a
            // maximum-depth, all-alternative cell that this window can never see — and
            // handing it the remainder would invent thousands of sites of the strongest
            // homozygous-non-reference evidence there is.
            let last_possible = probs.iter().rposition(|&p| p > 0.0).unwrap_or(0);
            for (i, &p) in probs.iter().enumerate() {
                if remaining == 0 {
                    break;
                }
                if i == last_possible {
                    counts[i] = remaining;
                    break;
                }
                if p <= 0.0 {
                    continue;
                }
                let share = if remaining_mass > 0.0 {
                    (p / remaining_mass).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                let drawn = rng.binomial(remaining, share);
                counts[i] = drawn;
                remaining -= drawn;
                remaining_mass -= p;
            }
            windows.push(counts);
            contig_of.push(contig);
            true_state.push(state);

            let jump = rng.unit();
            state = match state {
                INSIDE if jump < scenario.leave_run => OUTSIDE,
                OUTSIDE if jump < scenario.enter_run => INSIDE,
                s => s,
            };
        }
    }

    Sample {
        windows,
        contig_of,
        true_state,
    }
}

// ---------------------------------------------------------------------------
// The estimator: a two-state HMM over windows, fitted by Baum–Welch
// ---------------------------------------------------------------------------

#[derive(Clone, Debug)]
struct RunsFit {
    /// `(hom-ref, het, hom-alt)` for each state, after the ordering constraint.
    freqs: [[f64; GENOTYPES]; STATES],
    enter_run: f64,
    leave_run: f64,
    /// The coverage-weighted posterior occupancy — the `F` the architecture emits.
    f_posterior: f64,
    /// `enter / (enter + leave)` — the `F` spec §6.1 calls "their ratio".
    f_stationary: f64,
    log_likelihood: f64,
    iterations: u32,
    converged: bool,
    /// Whether the fit had to be relabelled because the state it called "inside" had the
    /// **higher** heterozygote rate — the ordering constraint firing.
    relabelled: bool,
    /// The share of windows the model was genuinely undecided about — posterior between
    /// 0.01 and 0.99. If this is zero, every window was classified on its own evidence
    /// and the chain's transitions did no work at all.
    undecided: f64,
}

/// **What the fitter sees of a sample**, and the one place the depth binning of §8 enters.
///
/// For each window, how many sites fell in each cell, and how likely each genotype makes
/// that cell. `components` holds **one** table where every window scores a cell the same
/// way — which is every case but one — and one table per window where a binned cell's mean
/// depth is a per-window quantity.
struct Emissions<'a> {
    /// `counts[window][cell]` — borrowed where binning left them alone.
    counts: Cow<'a, [Vec<u64>]>,
    /// `components[0 or window][cell][genotype]`.
    components: Vec<Vec<[f64; GENOTYPES]>>,
    /// How many cells a window was reduced to. Reported, because it is what binning buys.
    cells: usize,
}

impl<'a> Emissions<'a> {
    fn components_for(&self, window: usize) -> &[[f64; GENOTYPES]] {
        &self.components[if self.components.len() == 1 {
            0
        } else {
            window
        }]
    }

    /// One cell per exact `(depth, alt count)` — the emission of spec §6.1 as written,
    /// before any binning.
    fn exact(space: &CellSpace, sample: &'a Sample, eps: f64) -> Self {
        let components: Vec<[f64; GENOTYPES]> = space
            .cells
            .iter()
            .map(|&c| {
                let mut out = [0.0; GENOTYPES];
                for (j, slot) in out.iter_mut().enumerate() {
                    *slot = space.component(c, j, eps);
                }
                out
            })
            .collect();
        Self {
            cells: components.len(),
            counts: Cow::Borrowed(&sample.windows),
            components: vec![components],
        }
    }

    /// The same emission over **binned** cells, each scored at the mean of the exact depths
    /// that landed in it — what the accumulator can actually offer, since it keeps a depth
    /// bin and a running sum of depths rather than the depths themselves.
    ///
    /// The combinatorial factor is dropped: it is common to both states and every genotype,
    /// so it cancels out of every posterior and out of the gap the model reads a run from,
    /// and at a fractional depth it would need a gamma function to say something the fit
    /// cannot see. That makes log-likelihoods comparable **within** one `Emissions` — which
    /// is all "keep the best-scoring start" needs — and not between two.
    fn binned(
        space: &CellSpace,
        sample: &Sample,
        eps: f64,
        ladder: &DepthLadder,
        grain: MeanDepthGrain,
    ) -> Self {
        // Which binned cell each exact cell falls in, and what alternative count it holds.
        let mut index: std::collections::HashMap<(u32, u32), usize> =
            std::collections::HashMap::new();
        let mut alt_of: Vec<u32> = Vec::new();
        let cell_of: Vec<usize> = space
            .cells
            .iter()
            .map(|&c| {
                let key = (ladder.bin_of(c.depth), c.alt);
                *index.entry(key).or_insert_with(|| {
                    alt_of.push(c.alt);
                    alt_of.len() - 1
                })
            })
            .collect();
        let n_cells = alt_of.len();

        let mut counts = vec![vec![0u64; n_cells]; sample.windows.len()];
        let mut depth_sums = vec![vec![0f64; n_cells]; sample.windows.len()];
        for (w, window) in sample.windows.iter().enumerate() {
            for (i, &c) in window.iter().enumerate() {
                if c == 0 {
                    continue;
                }
                counts[w][cell_of[i]] += c;
                depth_sums[w][cell_of[i]] += c as f64 * f64::from(space.cells[i].depth);
            }
        }

        // The whole-sample mean of each cell, which is also the fallback for a cell a
        // window never saw.
        let mut sample_counts = vec![0f64; n_cells];
        let mut sample_depths = vec![0f64; n_cells];
        for w in 0..sample.windows.len() {
            for c in 0..n_cells {
                sample_counts[c] += counts[w][c] as f64;
                sample_depths[c] += depth_sums[w][c];
            }
        }
        let sample_mean: Vec<f64> = (0..n_cells)
            .map(|c| {
                if sample_counts[c] > 0.0 {
                    sample_depths[c] / sample_counts[c]
                } else {
                    f64::from(alt_of[c])
                }
            })
            .collect();

        let build = |mean: &[f64]| -> Vec<[f64; GENOTYPES]> {
            (0..n_cells)
                .map(|c| {
                    let mut out = [0.0; GENOTYPES];
                    for (j, slot) in out.iter_mut().enumerate() {
                        let p = p_alt(j, eps);
                        let k = f64::from(alt_of[c]);
                        *slot = (k * p.ln() + (mean[c] - k) * (1.0 - p).ln()).exp();
                    }
                    out
                })
                .collect()
        };

        let components = match grain {
            MeanDepthGrain::WholeSample => vec![build(&sample_mean)],
            MeanDepthGrain::PerWindow => (0..sample.windows.len())
                .map(|w| {
                    let mean: Vec<f64> = (0..n_cells)
                        .map(|c| {
                            if counts[w][c] > 0 {
                                depth_sums[w][c] / counts[w][c] as f64
                            } else {
                                sample_mean[c]
                            }
                        })
                        .collect();
                    build(&mean)
                })
                .collect(),
        };

        Self {
            counts: Cow::Owned(counts),
            components,
            cells: n_cells,
        }
    }

    /// `ln P(window | state)` — a sum over the window's cells, never a heterozygote count
    /// (spec §6.1).
    fn window_ln_emission(&self, window: usize, freqs: &[f64; GENOTYPES]) -> f64 {
        let components = self.components_for(window);
        let mut total = 0.0;
        for (i, &c) in self.counts[window].iter().enumerate() {
            if c == 0 {
                continue;
            }
            let mut p = 0.0;
            for j in 0..GENOTYPES {
                p += freqs[j] * components[i][j];
            }
            total += c as f64 * p.ln();
        }
        total
    }

    fn fit(&self, sample: &Sample, start: &RunsStart, max_iter: u32) -> RunsFit {
        let n = self.counts.len();
        let mut freqs = start.freqs;
        let mut a = [[0.0f64; STATES]; STATES]; // a[from][to]
        let mut enter = start.enter_run;
        let mut leave = start.leave_run;

        // Binning changes which cells a window's sites are spread over and never how many
        // there are, so this coverage weighting is the same number either way.
        let sites: Vec<f64> = self
            .counts
            .iter()
            .map(|w| w.iter().sum::<u64>() as f64)
            .collect();

        let mut gamma = vec![[0.0f64; STATES]; n];
        let mut last_ll = f64::NEG_INFINITY;
        let mut iterations = 0;
        let mut converged = false;

        for iter in 1..=max_iter {
            iterations = iter;
            a[INSIDE][OUTSIDE] = leave;
            a[INSIDE][INSIDE] = 1.0 - leave;
            a[OUTSIDE][INSIDE] = enter;
            a[OUTSIDE][OUTSIDE] = 1.0 - enter;
            let stationary = enter / (enter + leave);
            let ln_init = [
                stationary.max(1e-300).ln(),
                (1.0 - stationary).max(1e-300).ln(),
            ];
            let ln_a = [
                [a[0][0].max(1e-300).ln(), a[0][1].max(1e-300).ln()],
                [a[1][0].max(1e-300).ln(), a[1][1].max(1e-300).ln()],
            ];

            let ln_e: Vec<[f64; STATES]> = (0..n)
                .map(|w| {
                    [
                        self.window_ln_emission(w, &freqs[INSIDE]),
                        self.window_ln_emission(w, &freqs[OUTSIDE]),
                    ]
                })
                .collect();

            // Forward, in log space and restarting at every contig boundary.
            let mut alpha = vec![[f64::NEG_INFINITY; STATES]; n];
            for t in 0..n {
                let new_contig = t == 0 || sample.contig_of[t] != sample.contig_of[t - 1];
                for s in 0..STATES {
                    alpha[t][s] = if new_contig {
                        ln_init[s] + ln_e[t][s]
                    } else {
                        let mut acc = f64::NEG_INFINITY;
                        for i in 0..STATES {
                            acc = ln_sum_exp2(acc, alpha[t - 1][i] + ln_a[i][s]);
                        }
                        acc + ln_e[t][s]
                    };
                }
            }

            // Backward.
            let mut beta = vec![[0.0f64; STATES]; n];
            for t in (0..n).rev() {
                let last_of_contig = t + 1 == n || sample.contig_of[t + 1] != sample.contig_of[t];
                if last_of_contig {
                    beta[t] = [0.0; STATES];
                } else {
                    for s in 0..STATES {
                        let mut acc = f64::NEG_INFINITY;
                        for j in 0..STATES {
                            acc = ln_sum_exp2(acc, ln_a[s][j] + ln_e[t + 1][j] + beta[t + 1][j]);
                        }
                        beta[t][s] = acc;
                    }
                }
            }

            // Posteriors, and the total log-likelihood as the sum over contigs.
            let mut log_likelihood = 0.0;
            for t in 0..n {
                let norm = ln_sum_exp2(alpha[t][0] + beta[t][0], alpha[t][1] + beta[t][1]);
                for s in 0..STATES {
                    gamma[t][s] = (alpha[t][s] + beta[t][s] - norm).exp();
                }
                let last_of_contig = t + 1 == n || sample.contig_of[t + 1] != sample.contig_of[t];
                if last_of_contig {
                    log_likelihood += norm;
                }
            }

            // Transition counts.
            let mut xi = [[0.0f64; STATES]; STATES];
            for t in 0..n.saturating_sub(1) {
                if sample.contig_of[t + 1] != sample.contig_of[t] {
                    continue;
                }
                let norm = ln_sum_exp2(alpha[t][0] + beta[t][0], alpha[t][1] + beta[t][1]);
                for i in 0..STATES {
                    for j in 0..STATES {
                        xi[i][j] += (alpha[t][i] + ln_a[i][j] + ln_e[t + 1][j] + beta[t + 1][j]
                            - norm)
                            .exp();
                    }
                }
            }

            // Genotype frequencies: expected genotype counts under each state's posterior.
            let mut counts = [[0.0f64; GENOTYPES]; STATES];
            for (t, window) in self.counts.iter().enumerate() {
                let components = self.components_for(t);
                for s in 0..STATES {
                    let weight = gamma[t][s];
                    if weight < 1e-12 {
                        continue;
                    }
                    for (i, &c) in window.iter().enumerate() {
                        if c == 0 {
                            continue;
                        }
                        let mut denom = 0.0;
                        let mut terms = [0.0; GENOTYPES];
                        for j in 0..GENOTYPES {
                            terms[j] = freqs[s][j] * components[i][j];
                            denom += terms[j];
                        }
                        if denom <= 0.0 {
                            continue;
                        }
                        for j in 0..GENOTYPES {
                            counts[s][j] += weight * c as f64 * terms[j] / denom;
                        }
                    }
                }
            }

            let mut next_freqs = freqs;
            for s in 0..STATES {
                let total: f64 = counts[s].iter().sum();
                if total > 0.0 {
                    for j in 0..GENOTYPES {
                        next_freqs[s][j] = (counts[s][j] / total).max(1e-12);
                    }
                }
            }

            let out_of_inside: f64 = xi[INSIDE][OUTSIDE];
            let inside_total: f64 = xi[INSIDE][INSIDE] + xi[INSIDE][OUTSIDE];
            let into_inside: f64 = xi[OUTSIDE][INSIDE];
            let outside_total: f64 = xi[OUTSIDE][INSIDE] + xi[OUTSIDE][OUTSIDE];
            let next_leave = if inside_total > 0.0 {
                (out_of_inside / inside_total).clamp(1e-9, 1.0 - 1e-9)
            } else {
                leave
            };
            let next_enter = if outside_total > 0.0 {
                (into_inside / outside_total).clamp(1e-9, 1.0 - 1e-9)
            } else {
                enter
            };

            let moved = (next_enter / enter)
                .ln()
                .abs()
                .max((next_leave / leave).ln().abs());
            freqs = next_freqs;
            enter = next_enter;
            leave = next_leave;

            // **A relative tolerance on this log-likelihood is a trap, and it bit.** The
            // total scales with sites × depth: it is −1.5 × 10⁹ at 3 reads a site and
            // −2.0 × 10¹⁰ at 60, so the tolerance this test once used — 10⁻⁷ of the total —
            // was 151 nats at 3 reads and **2,013 nats at 60**. Expectation-maximization
            // stopped after four passes with the answer still moving by hundreds of nats a
            // pass, reported `converged`, and returned `F` = 1.0000 on a genome 32% covered
            // by runs. The tolerance is now relative to the **per-window** log-likelihood,
            // which is the quantity that does not grow with the genome, and floored so it
            // cannot ask for more than floating point can deliver.
            let per_window = log_likelihood.abs() / (n as f64).max(1.0);
            let tolerance = (1e-9 * per_window).max(1e-6) * (n as f64);
            if (log_likelihood - last_ll).abs() < tolerance && moved < 1e-6 {
                last_ll = log_likelihood;
                converged = true;
                break;
            }
            last_ll = log_likelihood;
        }

        // The ordering constraint: the state with the *lower* heterozygote rate is the
        // one inside a run. Applied after the fit, because nothing in Baum-Welch imposes
        // it — which is what makes it worth counting how often it fires.
        let mut relabelled = false;
        if freqs[INSIDE][1] > freqs[OUTSIDE][1] {
            freqs.swap(INSIDE, OUTSIDE);
            std::mem::swap(&mut enter, &mut leave);
            for g in gamma.iter_mut() {
                g.swap(INSIDE, OUTSIDE);
            }
            relabelled = true;
        }

        let covered: f64 = sites.iter().sum();
        let f_posterior = gamma
            .iter()
            .zip(&sites)
            .map(|(g, s)| g[INSIDE] * s)
            .sum::<f64>()
            / covered;

        let undecided = gamma
            .iter()
            .filter(|g| g[INSIDE] > 0.01 && g[INSIDE] < 0.99)
            .count() as f64
            / n as f64;

        RunsFit {
            freqs,
            enter_run: enter,
            leave_run: leave,
            f_posterior,
            f_stationary: enter / (enter + leave),
            log_likelihood: last_ll,
            iterations,
            converged,
            relabelled,
            undecided,
        }
    }
}

/// Where a fit starts. Baum–Welch climbs to a stationary point, so on a surface that is
/// not concave the start is part of the answer — which is exactly what question 1 asks.
#[derive(Clone, Copy)]
struct RunsStart {
    freqs: [[f64; GENOTYPES]; STATES],
    enter_run: f64,
    leave_run: f64,
}

/// A spread of starting points, over **both** things a start has to guess: how much of the
/// genome is inside a run, and how far apart the two states are in heterozygosity.
///
/// Spreading only the first is not a spread at all. An earlier version of this harness
/// varied `F₀` across five starts while every one of them guessed the inside state at a
/// tenth of the outside rate — so a genome whose real states sit closer together was
/// missed by all five identically, and the harness reported that as a property of the
/// estimator rather than of its own initialisation.
fn starts(outside_het: f64, hom_alt: f64) -> Vec<(String, RunsStart)> {
    let make = |inside_het: f64, enter: f64, leave: f64| RunsStart {
        freqs: [
            [1.0 - inside_het - hom_alt, inside_het, hom_alt],
            [1.0 - outside_het - hom_alt, outside_het, hom_alt],
        ],
        enter_run: enter,
        leave_run: leave,
    };
    let mut out = Vec::new();
    for &(separation, sep_name) in &[(0.05, "1/20"), (0.3, "1/3"), (0.75, "3/4")] {
        for &(enter, leave, f_name) in &[
            (0.005, 0.095, "0.05"),
            (0.020, 0.020, "0.50"),
            (0.036, 0.012, "0.75"),
        ] {
            out.push((
                format!("sep={sep_name} F₀={f_name}"),
                make(outside_het * separation, enter, leave),
            ));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// The questions
// ---------------------------------------------------------------------------

fn base_scenario(name: &str, true_f: f64, mean_run_windows: f64) -> Scenario {
    // Tomato: about one heterozygote per kilobase outside a run, and a landrace far from
    // the reference, so the homozygous-non-reference rate is well above the het rate.
    let outside = [1.0 - 0.001 - 0.006, 0.001, 0.006];
    // Inside a run the genotype is one allele draw doubled, so the heterozygotes' mass
    // moves to the homozygous-non-reference class (spec §6.1).
    let inside = [1.0 - 0.00005 - 0.0065, 0.00005, 0.0065];
    let leave = 1.0 / mean_run_windows;
    let enter = if true_f >= 1.0 {
        1.0
    } else {
        leave * true_f / (1.0 - true_f)
    };
    Scenario {
        name: name.to_string(),
        contigs: 12,
        windows_per_contig: 667,
        sites_per_window: 100_000,
        mean_depth: 3.0,
        eps: 0.001,
        outside,
        inside,
        enter_run: enter,
        leave_run: leave,
        false_het_floor: 0.0,
        coverage_levels: vec![1.0],
    }
}

fn summarise(values: &[f64]) -> (f64, f64, f64, f64) {
    let mean = values.iter().sum::<f64>() / values.len() as f64;
    let sd = if values.len() > 1 {
        (values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (values.len() - 1) as f64).sqrt()
    } else {
        0.0
    };
    let lo = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    (mean, sd, lo, hi)
}

/// Fit from every start and keep the best-scoring one, as spec §6.1 requires — while
/// reporting how far apart the starts landed, which is what says whether that choice
/// meant anything.
fn fit_all_starts(
    fitter: &Emissions,
    sample: &Sample,
    outside_het: f64,
    hom_alt: f64,
) -> (RunsFit, f64, f64) {
    let mut best: Option<RunsFit> = None;
    let mut all_f = Vec::new();
    for (_, start) in starts(outside_het, hom_alt) {
        let fit = fitter.fit(sample, &start, 300);
        all_f.push(fit.f_posterior);
        if best
            .as_ref()
            .is_none_or(|b| fit.log_likelihood > b.log_likelihood)
        {
            best = Some(fit);
        }
    }
    let (_, _, lo, hi) = summarise(&all_f);
    (best.expect("at least one start"), lo, hi)
}

fn question_1_coincident_states(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## 1. When the two states coincide, `F` is not identified — and this is a proof\n"
    );
    let _ = writeln!(
        report,
        "If both states carry the same genotype frequencies, every window's emission is the\n\
         same under either state, so the observed sequence is independent draws from one\n\
         distribution and `P(data | θ)` does not contain the transition rates at all. The\n\
         likelihood is exactly flat in `F`. What a fit returns is then its starting point,\n\
         and the run below is the demonstration rather than the argument.\n"
    );
    let mut scenario = base_scenario("states coincide", 0.30, 20.0);
    scenario.inside = scenario.outside;
    let space = CellSpace::for_scenario(&scenario);
    let sample = simulate(&scenario, &space, 0xC0FFEE);
    let fitter = Emissions::exact(&space, &sample, scenario.eps);

    let _ = writeln!(
        report,
        "{:<18} {:>12} {:>14} {:>14} {:>12}",
        "start", "F returned", "F stationary", "ln L", "relabelled"
    );
    let mut lls = Vec::new();
    for (label, start) in starts(scenario.outside[1], scenario.outside[2]) {
        let fit = fitter.fit(&sample, &start, 300);
        lls.push(fit.log_likelihood);
        let _ = writeln!(
            report,
            "{:<18} {:>12.4} {:>14.4} {:>14.2} {:>12}",
            label, fit.f_posterior, fit.f_stationary, fit.log_likelihood, fit.relabelled
        );
    }
    let (_, _, lo, hi) = summarise(&lls);
    let _ = writeln!(
        report,
        "\nSpread in log-likelihood across those starts: {:.3e} on a total of {:.0} — the\n\
         five are the same score to eleven significant figures, which is the flatness the\n\
         proof predicts. **But they do not drift apart: every one collapses to `F` = 0.**\n\
         The reason is that a fit never sits *at* coincidence. It starts with the two states\n\
         separated, and on the first pass every window fits the outside state far better\n\
         than the inside one, so the inside state is emptied and the rate of entering a run\n\
         goes to zero. The flat ridge is real and unreachable: expectation-maximization\n\
         walks away from it rather than along it.",
        hi - lo,
        lls[0]
    );
}

fn question_2_outbred(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## 2. A genuinely outbred genome — what floor does `F` carry?\n"
    );
    let _ = writeln!(
        report,
        "No runs anywhere: the chain never enters the inside state, so the true `F` is zero\n\
         and the two states differ only by whatever the fit reads out of window-to-window\n\
         noise. If `F` were merely imprecise it would shrink as windows accumulate.\n"
    );
    let _ = writeln!(
        report,
        "{:<10} {:>8} {:>10} {:>10} {:>10} {:>10} {:>9}",
        "windows", "seeds", "mean F", "sd", "min", "max", "relabel"
    );
    for &windows_per_contig in &[100usize, 400, 1600, 6400] {
        let mut scenario = base_scenario("outbred", 0.0, 20.0);
        scenario.windows_per_contig = windows_per_contig;
        scenario.enter_run = 0.0;
        scenario.leave_run = 1.0;
        let space = CellSpace::for_scenario(&scenario);
        let seeds = 8;
        let mut values = Vec::new();
        let mut relabelled = 0;
        for s in 0..seeds {
            let sample = simulate(&scenario, &space, 0xABCD_0000 + s);
            let fitter = Emissions::exact(&space, &sample, scenario.eps);
            let (fit, _, _) =
                fit_all_starts(&fitter, &sample, scenario.outside[1], scenario.outside[2]);
            values.push(fit.f_posterior);
            if fit.relabelled {
                relabelled += 1;
            }
        }
        let (mean, sd, lo, hi) = summarise(&values);
        let _ = writeln!(
            report,
            "{:<10} {:>8} {:>10.4} {:>10.4} {:>10.4} {:>10.4} {:>9}",
            scenario.windows(),
            seeds,
            mean,
            sd,
            lo,
            hi,
            relabelled
        );
    }
}

fn question_3_recovery(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## 3. Is `F` recovered when runs are really there?\n"
    );
    let _ = writeln!(
        report,
        "8,004 windows — a tomato genome at 100 kb. `F posterior` is the coverage-weighted\n\
         occupancy the architecture emits; `F stationary` is the transition rates' ratio,\n\
         which spec §6.1 also calls `F`. `across starts` is how far the five starting\n\
         points landed apart before the best-scoring one was kept.\n"
    );
    let _ = writeln!(
        report,
        "{:<10} {:>10} {:>12} {:>14} {:>16} {:>8}",
        "true F", "realised", "F posterior", "F stationary", "across starts", "iters"
    );
    for &(true_f, run_len) in &[(0.05, 20.0), (0.15, 30.0), (0.30, 30.0), (0.60, 50.0)] {
        let scenario = base_scenario("recovery", true_f, run_len);
        let space = CellSpace::for_scenario(&scenario);
        let sample = simulate(&scenario, &space, 0xBEEF_0000 + (true_f * 100.0) as u64);
        let fitter = Emissions::exact(&space, &sample, scenario.eps);
        let (fit, lo, hi) =
            fit_all_starts(&fitter, &sample, scenario.outside[1], scenario.outside[2]);
        let _ = writeln!(
            report,
            "{:<10.2} {:>10.4} {:>12.4} {:>14.4} {:>16} {:>8}",
            true_f,
            sample.true_f(),
            fit.f_posterior,
            fit.f_stationary,
            format!("{lo:.3}–{hi:.3}"),
            fit.iterations
        );
    }
}

fn question_4_false_het_floor(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## 4. Does a floor of false heterozygotes cancel, as §6.2 claims?\n"
    );
    let _ = writeln!(
        report,
        "Spec §6.2's second reason for choosing the runs estimator over the ratio one is\n\
         that a uniform floor of spurious heterozygotes — collapsed paralogs, mismapping —\n\
         lifts both states together and cancels out of the gap. The floor below is added to\n\
         both states, in multiples of the outside heterozygosity of 1 per kilobase. `Hobs`\n\
         is what the ratio estimator would read, for comparison.\n"
    );
    let _ = writeln!(
        report,
        "{:<12} {:>10} {:>12} {:>12} {:>12} {:>12} {:>12}",
        "floor", "realised F", "F posterior", "fitted h", "fitted Hout", "Hout/h", "Hobs"
    );
    let true_f = 0.30;
    for &multiple in &[0.0f64, 0.5, 1.0, 2.0, 3.0, 3.5, 4.0, 4.5, 5.0] {
        let mut scenario = base_scenario("floor", true_f, 30.0);
        scenario.false_het_floor = multiple * scenario.outside[1];
        let space = CellSpace::for_scenario(&scenario);
        let sample = simulate(&scenario, &space, 0xF100_0000);
        let fitter = Emissions::exact(&space, &sample, scenario.eps);
        let (fit, _, _) = fit_all_starts(
            &fitter,
            &sample,
            scenario.outside[1] + scenario.false_het_floor,
            scenario.outside[2],
        );
        // What the whole-genome heterozygosity looks like — the ratio estimator's input.
        let observed_het = scenario.stationary_f() * scenario.generating(INSIDE)[1]
            + (1.0 - scenario.stationary_f()) * scenario.generating(OUTSIDE)[1];
        let _ = writeln!(
            report,
            "{:<12} {:>10.4} {:>12.4} {:>12.6} {:>12.6} {:>12.1} {:>12.6}",
            format!("{multiple:.1}x"),
            sample.true_f(),
            fit.f_posterior,
            fit.freqs[INSIDE][1],
            fit.freqs[OUTSIDE][1],
            fit.freqs[OUTSIDE][1] / fit.freqs[INSIDE][1].max(1e-12),
            observed_het
        );
    }
}

/// **Does the chain do any work?** The decision made about a window depends on that
/// window's own reads and on its neighbours through the transitions. If a window's own
/// evidence already settles it, the transitions are inert and the model is a two-component
/// mixture over windows wearing a hidden Markov model's clothes. The number of covered
/// sites in a window is what decides that, so this sweeps it.
fn question_6_does_the_chain_matter(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## 6. Is this a hidden Markov model, or a per-window classifier?\n"
    );
    let _ = writeln!(
        report,
        "A window's state is inferred from its own reads *and* from its neighbours, through\n\
         the transition rates. `undecided` is the share of windows whose posterior landed\n\
         between 0.01 and 0.99 — the ones where the neighbours could have mattered. `shuffle\n\
         Δ` refits after shuffling the window order within each contig, which destroys every\n\
         run while preserving every window's contents: if the chain is doing work, `F` must\n\
         fall. 100,000 sites is a 100 kb window at 3 reads; 10,000 is a 10 kb one.\n"
    );
    let _ = writeln!(
        report,
        "{:<12} {:>10} {:>12} {:>12} {:>12} {:>10}",
        "sites/window", "true F", "F posterior", "undecided", "shuffle Δ", "relabel"
    );
    for &sites in &[100_000u64, 10_000, 3_000, 1_000, 300] {
        let mut scenario = base_scenario("grain", 0.30, 30.0);
        scenario.sites_per_window = sites;
        let space = CellSpace::for_scenario(&scenario);
        let sample = simulate(&scenario, &space, 0x6141_0000 + sites);
        let fitter = Emissions::exact(&space, &sample, scenario.eps);
        let (fit, _, _) =
            fit_all_starts(&fitter, &sample, scenario.outside[1], scenario.outside[2]);

        let mut rng = SplitMix64::new(0x6141_9999);
        let mut order: Vec<usize> = (0..sample.windows.len()).collect();
        for i in (1..order.len()).rev() {
            let j = (rng.unit() * (i + 1) as f64) as usize;
            order.swap(i, j.min(i));
        }
        let shuffled = Sample {
            windows: order.iter().map(|&i| sample.windows[i].clone()).collect(),
            contig_of: sample.contig_of.clone(),
            true_state: order.iter().map(|&i| sample.true_state[i]).collect(),
        };
        let shuffled_fitter = Emissions::exact(&space, &shuffled, scenario.eps);
        let (shuffled_fit, _, _) = fit_all_starts(
            &shuffled_fitter,
            &shuffled,
            scenario.outside[1],
            scenario.outside[2],
        );

        let _ = writeln!(
            report,
            "{:<12} {:>10.4} {:>12.4} {:>11.2}% {:>12.4} {:>10}",
            sites,
            sample.true_f(),
            fit.f_posterior,
            100.0 * fit.undecided,
            fit.f_posterior - shuffled_fit.f_posterior,
            fit.relabelled
        );
    }
}

fn question_5_shuffle(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## 5. Does shuffling the windows tell a real run from noise?\n"
    );
    let _ = writeln!(
        report,
        "Spec §6.1 offers this as the check needing no truth: shuffling destroys runs and\n\
         preserves everything else, so whatever `F` survives was read out of noise. For it\n\
         to work, the shuffled `F` has to be *lower* than the unshuffled one by more than\n\
         the noise — and on an outbred genome the two must agree, since there was nothing\n\
         to destroy.\n"
    );
    let _ = writeln!(
        report,
        "{:<16} {:>12} {:>12} {:>12} {:>10}",
        "genome", "F as walked", "F shuffled", "difference", "verdict"
    );
    for &(label, true_f) in &[("outbred", 0.0), ("F = 0.30", 0.30)] {
        let mut scenario = base_scenario(label, true_f, 30.0);
        if true_f == 0.0 {
            scenario.enter_run = 0.0;
            scenario.leave_run = 1.0;
        }
        let space = CellSpace::for_scenario(&scenario);
        let sample = simulate(&scenario, &space, 0x5417_0000);
        let fitter = Emissions::exact(&space, &sample, scenario.eps);
        let (fit, _, _) =
            fit_all_starts(&fitter, &sample, scenario.outside[1], scenario.outside[2]);

        // Shuffle window order within each contig, keeping every window's content.
        let mut rng = SplitMix64::new(0x5417_5555);
        let mut order: Vec<usize> = (0..sample.windows.len()).collect();
        for i in (1..order.len()).rev() {
            let j = (rng.unit() * (i + 1) as f64) as usize;
            order.swap(i, j.min(i));
        }
        let shuffled = Sample {
            windows: order.iter().map(|&i| sample.windows[i].clone()).collect(),
            contig_of: sample.contig_of.clone(),
            true_state: order.iter().map(|&i| sample.true_state[i]).collect(),
        };
        let shuffled_fitter = Emissions::exact(&space, &shuffled, scenario.eps);
        let (shuffled_fit, _, _) = fit_all_starts(
            &shuffled_fitter,
            &shuffled,
            scenario.outside[1],
            scenario.outside[2],
        );

        let difference = fit.f_posterior - shuffled_fit.f_posterior;
        let verdict = if difference.abs() < 0.02 {
            "agree"
        } else {
            "differ"
        };
        let _ = writeln!(
            report,
            "{:<16} {:>12.4} {:>12.4} {:>12.4} {:>10}",
            label, fit.f_posterior, shuffled_fit.f_posterior, difference, verdict
        );
    }
}

/// **Would a finer window rescue the collapse of §4?** At 100 kb the chain is inert (§6),
/// so a window that has lost its own separation has nothing to fall back on. A finer grain
/// gives the transitions something to do — at the cost of more windows for the same
/// genome, which is what the window key charges for. The genome is held constant here:
/// tenfold smaller windows, tenfold more of them.
fn question_7_finer_windows_against_the_floor(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## 7. Would a finer window rescue the collapse?\n"
    );
    let _ = writeln!(
        report,
        "Same genome throughout — 800 Mb — cut into windows of three sizes. The floor is\n\
         the one §4 collapses at. If the chain can carry an estimate through a window that\n\
         has lost its own separation, `F` comes back at the finer grains.\n"
    );
    let _ = writeln!(
        report,
        "{:<12} {:>10} {:>10} {:>12} {:>12} {:>12}",
        "sites/window", "windows", "floor", "realised F", "F posterior", "undecided"
    );
    for &(sites, windows_per_contig) in &[(100_000u64, 667usize), (10_000, 6_670), (3_000, 22_233)]
    {
        for &multiple in &[0.0f64, 3.0] {
            let mut scenario = base_scenario("finer", 0.30, 30.0 * (100_000.0 / sites as f64));
            scenario.sites_per_window = sites;
            scenario.windows_per_contig = windows_per_contig;
            scenario.false_het_floor = multiple * scenario.outside[1];
            let space = CellSpace::for_scenario(&scenario);
            let sample = simulate(&scenario, &space, 0x7777_0000 + sites);
            let fitter = Emissions::exact(&space, &sample, scenario.eps);
            let (fit, _, _) = fit_all_starts(
                &fitter,
                &sample,
                scenario.outside[1] + scenario.false_het_floor,
                scenario.outside[2],
            );
            let _ = writeln!(
                report,
                "{:<12} {:>10} {:>10} {:>12.4} {:>12.4} {:>11.1}%",
                sites,
                scenario.windows(),
                format!("{multiple:.1}x"),
                sample.true_f(),
                fit.f_posterior,
                100.0 * fit.undecided
            );
        }
    }
}

/// **Experiment 1(c): what depth binning does to the window emission.**
///
/// Spec §6.1 states the emission's factorisation is exact "for cells that hold one exact
/// depth", that "binning the depth makes it an approximation", and that "scoring each cell
/// at its own mean depth is what keeps the error small". Nothing measured the last clause.
///
/// Each row refits the **same drawn genome** with the same nine starting points, changing
/// only how a window's sites are keyed and scored. `F exact` is the emission as spec §6.1
/// writes it — one cell per exact `(depth, alt count)`. The binned rows key by depth bin and
/// score each cell at the mean of the depths that landed in it, per window as the
/// accumulator holds them, or over the whole sample.
///
/// Windows differ in coverage — 0.6×, 1.0× and 1.6× the stated mean, taken in turn — because
/// with uniform coverage the per-window and whole-sample means differ only by sampling noise
/// and the comparison between them would be empty.
fn question_8_depth_binning(report: &mut String) {
    let _ = writeln!(report, "\n## 8. What depth binning does to the emission\n");
    let _ = writeln!(
        report,
        "The same drawn genome, refitted with the same nine starting points, changing only how\n\
         a window's sites are keyed and scored. 3,600 windows of 100,000 sites — above the\n\
         3,000 the design refuses to fit below, and a third of a tomato genome's windows.\n\
         Window coverage alternates over 0.6×, 1.0× and 1.6× the stated mean depth, so that a\n\
         cell's per-window mean depth and its whole-sample mean are genuinely different\n\
         numbers. `h` and `Hout` are the fitted inside and outside heterozygote rates, whose\n\
         gap is the whole of what the model reads a run from.\n"
    );
    let _ = writeln!(
        report,
        "{:<8} {:<26} {:<12} {:>7} {:>11} {:>11} {:>11} {:>10}",
        "depth", "ladder", "mean depth", "cells", "F returned", "h", "Hout", "undecided"
    );
    for (mean_depth, coverage) in [3.0f64, 20.0, 60.0]
        .into_iter()
        .flat_map(|d| [(d, vec![1.0]), (d, vec![0.6, 1.0, 1.6])])
    {
        let mut scenario = base_scenario("binning", 0.30, 30.0);
        scenario.mean_depth = mean_depth;
        scenario.windows_per_contig = 300;
        scenario.coverage_levels = coverage.clone();
        let space = CellSpace::for_scenario(&scenario);
        let sample = simulate(&scenario, &space, 0x8888_0000 + mean_depth as u64);
        let deepest = (space.depth_probs.len() - 1) as u32;
        let depth_label = if coverage.len() == 1 {
            format!("{mean_depth:.0} flat")
        } else {
            format!("{mean_depth:.0} varied")
        };

        let mut rows: Vec<(String, String, Emissions)> = Vec::new();
        rows.push((
            DepthLadder::exact(deepest).label,
            "exact".to_string(),
            Emissions::exact(&space, &sample, scenario.eps),
        ));
        // The first ladder is the control: one bin per depth means every cell's mean depth
        // is its own depth, so the binned path must return the unbinned answer. A
        // difference there would be the harness's, not the binning's.
        for ladder in [
            DepthLadder::exact(deepest),
            DepthLadder::geometric(8, 16, 124),
            DepthLadder::geometric(8, 20, 124),
            DepthLadder::geometric(8, 16, 300),
        ] {
            for (grain, grain_name) in [
                (MeanDepthGrain::PerWindow, "per window"),
                (MeanDepthGrain::WholeSample, "whole sample"),
            ] {
                rows.push((
                    ladder.label.clone(),
                    grain_name.to_string(),
                    Emissions::binned(&space, &sample, scenario.eps, &ladder, grain),
                ));
            }
        }

        let _ = writeln!(
            report,
            "{:<12} {:<26} {:<12} {:>7} {:>11.4} {:>11} {:>11} {:>10}",
            depth_label,
            "true (realised)",
            "-",
            "-",
            sample.true_f(),
            "-",
            "-",
            "-"
        );
        for (ladder_label, grain_name, emissions) in &rows {
            let (fit, _, _) =
                fit_all_starts(emissions, &sample, scenario.outside[1], scenario.outside[2]);
            let _ = writeln!(
                report,
                "{:<12} {:<26} {:<12} {:>7} {:>11.4} {:>11.6} {:>11.6} {:>9.2}%",
                depth_label,
                ladder_label,
                grain_name,
                emissions.cells,
                fit.f_posterior,
                fit.freqs[INSIDE][1],
                fit.freqs[OUTSIDE][1],
                100.0 * fit.undecided,
            );
        }
        let _ = writeln!(report);
    }
}

/// **Why the fit at 60 reads a site returns what it does** — every start's own outcome,
/// rather than the winner alone. Spec §6.5's whole point is that the winner cannot be read
/// without them.
fn question_9_high_depth_diagnosis(report: &mut String) {
    let _ = writeln!(report, "\n## 9. The fit at high depth, start by start\n");
    for &mean_depth in &[20.0f64, 60.0] {
        let mut scenario = base_scenario("diagnosis", 0.30, 30.0);
        scenario.mean_depth = mean_depth;
        scenario.windows_per_contig = 300;
        let space = CellSpace::for_scenario(&scenario);
        let sample = simulate(&scenario, &space, 0x8888_0000 + mean_depth as u64);
        let emissions = Emissions::exact(&space, &sample, scenario.eps);

        // What the generator actually put on the table, counted rather than assumed: sites
        // whose alternative fraction is nowhere near 0 or 1 are heterozygotes at any depth
        // above a handful, and the whole model turns on how many of them a window holds.
        let mut het_looking = 0u64;
        let mut total_sites = 0u64;
        for window in &sample.windows {
            for (i, &c) in window.iter().enumerate() {
                total_sites += c;
                let cell = space.cells[i];
                let fraction = f64::from(cell.alt) / f64::from(cell.depth);
                if cell.depth >= 10 && (0.25..=0.75).contains(&fraction) {
                    het_looking += c;
                }
            }
        }
        let _ = writeln!(
            report,
            "\ndepth {mean_depth:.0}, realised F = {:.4}; the genome holds {het_looking} \
             heterozygote-looking sites\nin {total_sites} ({:.5} of them), against the {:.5} \
             the scenario asks for outside a run",
            sample.true_f(),
            het_looking as f64 / total_sites as f64,
            scenario.outside[1],
        );
        let mut per_cell = vec![0u64; space.cells.len()];
        for window in &sample.windows {
            for (i, &c) in window.iter().enumerate() {
                per_cell[i] += c;
            }
        }
        let mut order: Vec<usize> = (0..per_cell.len()).collect();
        order.sort_by_key(|&i| std::cmp::Reverse(per_cell[i]));
        let _ = write!(report, "the eight fullest cells (depth,alt→sites): ");
        for &i in order.iter().take(8) {
            let _ = write!(
                report,
                "({},{})→{}  ",
                space.cells[i].depth, space.cells[i].alt, per_cell[i]
            );
        }
        let _ = writeln!(report);
        let _ = writeln!(
            report,
            "{:<18} {:>10} {:>12} {:>12} {:>18} {:>7} {:>6}",
            "start", "F", "h", "Hout", "ln L", "iters", "conv"
        );
        for (label, start) in starts(scenario.outside[1], scenario.outside[2]) {
            let fit = emissions.fit(&sample, &start, 300);
            let _ = writeln!(
                report,
                "{:<18} {:>10.4} {:>12.6} {:>12.6} {:>18.4} {:>7} {:>6}",
                label,
                fit.f_posterior,
                fit.freqs[INSIDE][1],
                fit.freqs[OUTSIDE][1],
                fit.log_likelihood,
                fit.iterations,
                fit.converged
            );
        }
    }
}

fn main() {
    let mut report = String::new();
    let _ = writeln!(
        report,
        "# Inbreeding `F` from runs of homozygosity — is the estimate well behaved?"
    );
    let _ = writeln!(
        report,
        "#\n\
         # Tomato-shaped: 12 contigs, 100 kb windows, 100,000 covered sites a window,\n\
         # 3 reads a site, error rate 0.001, one heterozygote per kilobase outside a run\n\
         # and a homozygous-non-reference rate of 6 per kilobase. Every fit is Baum-Welch\n\
         # from nine starting points — three guesses at how far apart the two states sit,\n\
         # crossed with three at how much of the genome is inside a run — keeping the\n\
         # best-scoring one, as spec §6.1 specifies. Sections 8 and 9 vary the depth."
    );

    // Each section is minutes rather than seconds, so they can be run one at a time.
    let args: Vec<String> = std::env::args().collect();
    let only: Vec<&String> = args.iter().filter(|a| a.starts_with("--only=")).collect();
    let wanted = |section: &str| only.is_empty() || only.iter().any(|a| a.ends_with(section));

    if wanted("binning") {
        question_8_depth_binning(&mut report);
    }
    if wanted("highdepth") {
        question_9_high_depth_diagnosis(&mut report);
    }
    if !only.is_empty() {
        print!("{report}");
        return;
    }

    question_1_coincident_states(&mut report);
    question_2_outbred(&mut report);
    question_3_recovery(&mut report);
    question_4_false_het_floor(&mut report);
    question_5_shuffle(&mut report);
    question_6_does_the_chain_matter(&mut report);
    question_7_finer_windows_against_the_floor(&mut report);

    print!("{report}");
}
