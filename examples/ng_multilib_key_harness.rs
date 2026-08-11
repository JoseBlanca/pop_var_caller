//! **Does a multi-library cell key give an unbiased estimate?** — the harness that
//! decides, rather than the sweep that compares one draw against another.
//!
//! The question a review keeps re-opening is which scoring rule to use when a cell key
//! has thrown away *which library produced which read*. Three properties are wanted of
//! whatever is chosen: the estimate should be unbiased, it should approach the truth as
//! the data grow, and it should be as precise as the key allows. This program measures
//! all three, and separates them, which the previous comparison could not.
//!
//! **The central measurement needs no simulation.** An estimator maximises a sum over
//! cells. Replace the observed cell counts with the cell's *exact probability under the
//! truth* and you get the objective the estimator is climbing with an infinite genome:
//!
//! ```text
//! Q(θ)  =  Σ over cells   P_true(cell) · ln L_rule(cell ; θ)
//! θ*    =  argmax Q(θ)              the value the estimate converges to
//! bias  =  θ* − θ_true              a fixed number, not a sampling error
//! ```
//!
//! `θ*` is the pseudo-true value a misspecified maximum-likelihood fit converges to
//! (White 1982) — the parameter whose model sits closest in Kullback–Leibler divergence
//! to the truth. Computing it takes no random draws, no seeds and no repeats, so
//! `bias = 0` and `bias ≠ 0` are decided exactly rather than to within Monte Carlo
//! noise. **`bias = 0` is exactly the statement that the estimator is consistent.**
//!
//! Monte Carlo is still run, for two things the analytic calculation cannot give: a
//! check that the implementation actually converges to `θ*`, and the *variance*, which
//! is what "as precise as possible" means.
//!
//! ## The candidates
//!
//! A candidate is a **coarsening** (what the key keeps) crossed with a **scoring rule**
//! (how a coarsened observation is scored).
//!
//! Coarsenings, parameterised by `K` = how many alternative reads keep their library
//! attribution, and `D` = the total depth at or below which the whole per-library
//! breakdown is kept:
//!
//! - `K=0, D=0` — today's pooled key: total depth and total alternative count;
//! - `K=4, D=0` — attribution of the alternative reads only;
//! - `K=4, D=4` — attribution, plus the whole breakdown at shallow sites;
//! - `exact` — nothing coarsened. The oracle, and the precision ceiling.
//!
//! Scoring rules for a coarsened cell:
//!
//! - `plug-avg` — each library's depth is taken as its average share, `n̂_g = w_g·n`,
//!   with `max(n̂_g − k_g, 0)` where that goes negative. **This is what the design
//!   proposes today.**
//! - `plug-res` — the alternative reads are placed first and only the remainder is
//!   split by share, `n̂_g = k_g + w_g·(n − k)`. Never negative, so no clamp.
//! - `marginal` — the likelihood of the coarsened observation: the sum of the exact
//!   probabilities of every exact outcome the key maps to. There is a closed form (see
//!   `ln_component_*` below), so this costs the same as the plug-ins.
//!
//! ## The checks that need no simulation at all
//!
//! Run before any fitting, because they are identities:
//!
//! 1. **Proper likelihood** — the rule sums to one over the cell space at any `θ`. A
//!    rule that does not is not the likelihood of anything and no consistency result
//!    covers it.
//! 2. **Non-negative exponents** — `n̂_g ≥ k_g` at every reachable cell, so no factor
//!    exceeds one. Reported as the probability mass on which it fails.
//! 3. **Equal error rates** — with every `ε_g` equal, every rule must reproduce the
//!    exact likelihood to floating point, since there is then nothing to attribute.
//! 4. **One library** — every rule collapses to the pooled key.
//!
//! ## Two further questions the same machinery answers
//!
//! Both are misspecifications rather than coarsenings — the fit is handed a model that is
//! not the one the data came from — and both are priced the same exact way.
//!
//! - **What depth binning costs.** The accumulator keys a cell by a depth *bin* and scores
//!   it at the mean of the exact depths that fell in it, so the closed form above is handed
//!   a fractional `n`. Sweeps the number of bins, how far up the exact-per-depth region
//!   runs, and the cap. **The exact ladder — one bin per depth — is the control, and must
//!   return the unbinned answer.**
//! - **What assuming a heterozygote is a half costs.** Generate at a true allele balance
//!   `b`, fit with the model that assumes `½ + ε/3`, and read the bias in `ε` and in the two
//!   genotype frequencies. **`b` = 0.50 is the control, and must return exactly zero.**
//!
//! ```text
//! cargo run --release --example ng_multilib_key_harness
//! cargo run --release --example ng_multilib_key_harness -- --monte-carlo
//! cargo run --release --example ng_multilib_key_harness -- --only=binning
//! cargo run --release --example ng_multilib_key_harness -- --only=balance
//! ```

// `j` here is a **genotype** — how many of the individual's copies are non-reference —
// and the arrays it steps through are indexed by it. Rewriting those loops as
// `iter_mut().enumerate()` would hide the quantity behind an iterator position, on the
// one file whose arithmetic is the oracle every later change to the scoring rule is
// checked against.
#![allow(clippy::needless_range_loop)]

use std::collections::HashMap;
use std::fmt::Write as _;

const PLOIDY: u32 = 2;
const GENOTYPES: usize = (PLOIDY + 1) as usize;

// ---------------------------------------------------------------------------
// Small numerical helpers
// ---------------------------------------------------------------------------

/// Lanczos approximation, g = 7. Accurate to about 15 digits over the range used
/// here (integers up to a few thousand), which is all the factorials need.
fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_1,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x)
    } else {
        let y = x - 1.0;
        let mut a = C[0];
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (y + i as f64);
        }
        let t = y + 7.5;
        0.5 * (2.0 * std::f64::consts::PI).ln() + (y + 0.5) * t.ln() - t + a.ln()
    }
}

fn ln_factorial(n: u32) -> f64 {
    ln_gamma(f64::from(n) + 1.0)
}

/// `ln x!` for a **fractional** `x`, which is what a cell scored at a mean depth needs.
/// Returns zero where `x` is below `−1`, since the factor is then undefined — and it is a
/// factor common to every genotype, so it cancels out of the mixture and cannot move the
/// fit. The caller counts those cells rather than letting them vanish silently.
fn ln_factorial_real(x: f64) -> f64 {
    if x <= -1.0 + 1e-9 {
        0.0
    } else {
        ln_gamma(x + 1.0)
    }
}

fn ln_binomial_coefficient(n: u32, k: u32) -> f64 {
    ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
}

fn ln_binomial_coefficient_real(n: f64, k: u32) -> f64 {
    ln_factorial_real(n) - ln_factorial(k) - ln_factorial_real(n - f64::from(k))
}

fn ln_sum_exp(values: &[f64]) -> f64 {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    max + values.iter().map(|v| (v - max).exp()).sum::<f64>().ln()
}

/// `x·ln y` with the convention `0·ln 0 = 0`, so a zero count against a zero
/// probability contributes nothing instead of a NaN.
fn x_ln_y(x: f64, y: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x * y.ln() }
}

// ---------------------------------------------------------------------------
// The model
// ---------------------------------------------------------------------------

/// One read's chance of showing a non-reference base, at `j` alternative copies out of
/// `PLOIDY` and error rate `eps` (spec `parameter_prepass.md` §3).
///
/// `balance` is the chance that a read at a **heterozygote** came off the alternative copy
/// before any misread. The model assumes `MODEL_HET_BALANCE`; a world may generate from a
/// different one, which is the misspecification spec `parameter_prepass_generic.md` §8
/// proposes spending a parameter on. At `balance = ½` this is the spec's `p_j` exactly, and
/// the two homozygous rows do not contain `balance` at all.
fn p_alt(j: u32, eps: f64, balance: f64) -> f64 {
    let frac = if j == 0 {
        0.0
    } else if j == PLOIDY {
        1.0
    } else {
        balance
    };
    frac * (1.0 - eps / 3.0) + (1.0 - frac) * eps
}

/// What every fit in this program assumes a heterozygote's reads look like: half of them
/// off each copy. Spec `parameter_prepass.md` §3's `½`.
const MODEL_HET_BALANCE: f64 = 0.5;

// ---------------------------------------------------------------------------
// The depth ladder
// ---------------------------------------------------------------------------

/// Which depths share a bin. The real accumulator does not keep an exact depth: it keys a
/// cell by a bin, one bin per depth at the bottom and widening geometrically above, and
/// scores each cell at the **mean of the exact depths that landed in it**
/// (`../doc/devel/ng/arch/parameter_prepass_generic.md` §2.2).
///
/// The exact ladder — one bin per depth — is the control this program checks first: under
/// it every cell's mean depth is its own depth, so binned scoring must reproduce the
/// unbinned answer to floating point, and any bias it reports is the harness's.
#[derive(Clone)]
struct DepthLadder {
    label: String,
    /// `bin_of[d]` for `d` in `0..=cap`.
    bin_of: Vec<u32>,
    /// The deepest depth the ladder has a bin for. A site above it would be subsampled
    /// down to it by the accumulator; no world here reaches it, and the check in
    /// `binning_checks` reports the truth mass that does.
    cap: u32,
    bins: usize,
}

impl DepthLadder {
    /// One bin per depth to `max_depth`. The null case.
    fn exact(max_depth: u32) -> Self {
        Self {
            label: "exact (no binning)".to_string(),
            bin_of: (0..=max_depth).collect(),
            cap: max_depth,
            bins: (max_depth + 1) as usize,
        }
    }

    /// Exact integers to `exact_to`, then `bins − exact_to − 1` geometrically widening bins
    /// up to `cap`.
    fn geometric(exact_to: u32, bins: usize, cap: u32) -> Self {
        assert!(bins > exact_to as usize + 1, "no room to widen");
        let widening = bins - exact_to as usize - 1;
        let ratio = (f64::from(cap) / f64::from(exact_to)).powf(1.0 / widening as f64);
        // The first depth of each widening bin, forced strictly increasing.
        let mut first: Vec<u32> = Vec::with_capacity(widening);
        let mut previous = exact_to;
        for i in 0..widening {
            let edge = (f64::from(exact_to) * ratio.powi(i as i32 + 1)).round() as u32;
            let edge = edge.max(previous + 1).min(cap);
            first.push(previous + 1);
            previous = edge;
        }
        let mut bin_of = Vec::with_capacity((cap + 1) as usize);
        for depth in 0..=cap {
            let bin = if depth <= exact_to {
                depth
            } else {
                let above = first.iter().filter(|&&f| f <= depth).count() as u32;
                exact_to + above
            };
            bin_of.push(bin);
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

    /// The widest bin's span, as a readable summary of how coarse the ladder gets.
    fn widest_bin(&self) -> u32 {
        let mut widths = vec![0u32; self.bins];
        for (depth, &bin) in self.bin_of.iter().enumerate() {
            widths[bin as usize] = widths[bin as usize].max(depth as u32);
        }
        let mut lows = vec![u32::MAX; self.bins];
        for (depth, &bin) in self.bin_of.iter().enumerate() {
            lows[bin as usize] = lows[bin as usize].min(depth as u32);
        }
        (0..self.bins)
            .map(|b| widths[b] - lows[b] + 1)
            .max()
            .unwrap_or(1)
    }
}

/// Which depth a binned cell is scored at.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum DepthScoring {
    /// The mean of the exact depths that landed in **this cell** — what
    /// `DepthAltHistogram::mean_depth_in_cell` computes, and what the arch doc argues is a
    /// correctness requirement rather than a refinement.
    PerCell,
    /// The mean over the whole **bin**, alternative counts pooled. The variant the arch doc
    /// says diverges; measured here rather than asserted.
    PerBin,
}

/// Everything that defines one simulated world: the libraries, their error rates and
/// read shares, the depth distribution, and the genotype frequencies.
#[derive(Clone)]
struct World {
    name: String,
    /// One error rate per library.
    eps: Vec<f64>,
    /// Each library's share of the reads. Sums to 1.
    weights: Vec<f64>,
    mean_depth: f64,
    pi_het: f64,
    pi_hom_alt: f64,
    /// The chance a read at a heterozygote came off the alternative copy, **in the world**.
    /// Every fit assumes `MODEL_HET_BALANCE`; a world that generates from something else is
    /// the misspecification of spec §8.
    het_balance: f64,
    /// **A second class of site**: with probability `share`, a site's reads disagree with the
    /// reference at `rate` instead of at their library's own `ε`
    /// (`research/noise_model_overdispersion_2026-08-10.md`). `None` is a world with one
    /// class, which is what every world here was before the noise-model milestone.
    ///
    /// **One rate, shared by every library at a noisy site, while `ε` stays per library** —
    /// which is the production model's own asymmetry and the reason this harness can say
    /// something no fixture in `src/` can: whether the estimate stays unbiased once the cell
    /// key has thrown away which library produced each alternative read.
    site_noise: Option<(f64, f64)>,
}

impl World {
    fn libraries(&self) -> usize {
        self.eps.len()
    }

    fn genotype_freqs(&self) -> [f64; GENOTYPES] {
        [
            1.0 - self.pi_het - self.pi_hom_alt,
            self.pi_het,
            self.pi_hom_alt,
        ]
    }

    /// The deepest site the world produces, past the Poisson truncation.
    fn max_depth(&self) -> u32 {
        (self.depth_distribution().len() - 1) as u32
    }

    /// The depth support and its probabilities, truncated where the Poisson tail falls
    /// below `1e-13` and renormalised over `n ≥ 1` — a site with no reads is not a site.
    fn depth_distribution(&self) -> Vec<f64> {
        let lambda = self.mean_depth;
        let mut n_max = lambda.ceil() as u32 + 1;
        loop {
            let ln_p = -lambda + f64::from(n_max) * lambda.ln() - ln_factorial(n_max);
            if ln_p.exp() < 1e-13 && f64::from(n_max) > lambda {
                break;
            }
            n_max += 1;
        }
        let mut probs = vec![0.0; (n_max + 1) as usize];
        for (n, slot) in probs.iter_mut().enumerate() {
            let n = n as u32;
            *slot = (-lambda + f64::from(n) * lambda.ln() - ln_factorial(n)).exp();
        }
        probs[0] = 0.0;
        let total: f64 = probs.iter().sum();
        for p in &mut probs {
            *p /= total;
        }
        probs
    }
}

/// The parameters being fitted: one error rate per library plus the genotype
/// frequencies. Everything else — the read shares, the depth distribution — is known.
#[derive(Clone, Debug)]
struct Params {
    eps: Vec<f64>,
    freqs: [f64; GENOTYPES],
}

// ---------------------------------------------------------------------------
// Cells: what a key keeps about one site
// ---------------------------------------------------------------------------

/// One cell of a keyed table. `Whole` and `Attributed` and `Pooled` are the three arms
/// of the proposed key; `Whole` with no depth bound is the exact per-library oracle.
///
/// **The two coarsened arms are keyed by a depth *bin*, not by a depth**, which is what the
/// accumulator does. Under `DepthLadder::exact` a bin is one depth and the key is the
/// unbinned one, so the same code covers both and the exact ladder is a control rather than
/// a separate path.
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
enum Cell {
    /// The whole per-library breakdown, exact depths included. Never binned: it is the
    /// oracle, and the point of it is that nothing has been thrown away.
    Whole(Vec<(u32, u32)>),
    /// Depth bin, and which library each alternative read came from.
    Attributed { bin: u32, alt: Vec<u32> },
    /// Depth bin and total alternative count.
    Pooled { bin: u32, alt: u32 },
}

/// How much a key keeps: alternative reads up to `max_attributed` keep their library,
/// and sites of total depth at most `exact_below_depth` keep everything.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Coarsening {
    max_attributed: u32,
    exact_below_depth: u32,
    /// The oracle: nothing is coarsened at all.
    exact: bool,
}

impl Coarsening {
    fn pooled() -> Self {
        Self {
            max_attributed: 0,
            exact_below_depth: 0,
            exact: false,
        }
    }
    fn attributed(max_attributed: u32, exact_below_depth: u32) -> Self {
        Self {
            max_attributed,
            exact_below_depth,
            exact: false,
        }
    }
    fn oracle() -> Self {
        Self {
            max_attributed: 0,
            exact_below_depth: 0,
            exact: true,
        }
    }

    fn key_of(&self, per_lib: &[(u32, u32)], ladder: &DepthLadder) -> Cell {
        let depth: u32 = per_lib.iter().map(|&(n, _)| n).sum();
        let alt: u32 = per_lib.iter().map(|&(_, k)| k).sum();
        if self.exact || depth <= self.exact_below_depth {
            Cell::Whole(per_lib.to_vec())
        } else if alt <= self.max_attributed {
            Cell::Attributed {
                bin: ladder.bin_of(depth),
                alt: per_lib.iter().map(|&(_, k)| k).collect(),
            }
        } else {
            Cell::Pooled {
                bin: ladder.bin_of(depth),
                alt,
            }
        }
    }
}

/// How a coarsened cell is scored.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Rule {
    /// `n̂_g = w_g·n`, clamped at zero. What the design proposes today.
    PlugAverageShare,
    /// `n̂_g = k_g + w_g·(n − k)`. Never negative.
    PlugResidualShare,
    /// The sum of the exact probabilities of every outcome the key maps to.
    Marginal,
}

impl Rule {
    fn label(self) -> &'static str {
        match self {
            Rule::PlugAverageShare => "plug-avg",
            Rule::PlugResidualShare => "plug-res",
            Rule::Marginal => "marginal",
        }
    }
}

// ---------------------------------------------------------------------------
// Component likelihoods — the per-genotype probability of one cell
// ---------------------------------------------------------------------------

/// `ln P(cell | genotype j)` for the whole breakdown. Exact by construction: the
/// per-library depths are observed, so this is the multinomial depth split times one
/// binomial per library.
fn ln_component_whole(
    per_lib: &[(u32, u32)],
    weights: &[f64],
    eps: &[f64],
    balance: f64,
    j: u32,
) -> f64 {
    let depth: u32 = per_lib.iter().map(|&(n, _)| n).sum();
    let mut total = ln_factorial(depth);
    for (g, &(n_g, k_g)) in per_lib.iter().enumerate() {
        let p = p_alt(j, eps[g], balance);
        total += -ln_factorial(n_g) + x_ln_y(f64::from(n_g), weights[g]);
        total += ln_binomial_coefficient(n_g, k_g)
            + x_ln_y(f64::from(k_g), p)
            + x_ln_y(f64::from(n_g - k_g), 1.0 - p);
    }
    total
}

/// `ln P(cell | genotype j)` for the attributed key, in closed form.
///
/// Each read independently picks a library (share `w_g`) and then shows the alternative
/// allele (`p_j(ε_g)`) or not. The cell records how many alternative reads came from
/// each library and how many reads showed the reference in total, so the outcome is a
/// multinomial over `G + 1` categories: one "alternative from library g" per library,
/// and one pooled "showed the reference".
///
/// **`depth` is a real number, because a binned cell is scored at the mean of the exact
/// depths that fell in it.** Everything except the two factorial prefactors is already
/// continuous in the depth; the prefactors are common to every genotype, so they cancel out
/// of the mixture and cannot move the fit whatever is done with them.
fn ln_component_attributed(
    depth: f64,
    alt_by_lib: &[u32],
    weights: &[f64],
    eps: &[f64],
    balance: f64,
    j: u32,
) -> f64 {
    let alt: u32 = alt_by_lib.iter().sum();
    let reference_reads = depth - f64::from(alt);
    let mut total = ln_factorial_real(depth) - ln_factorial_real(reference_reads);
    let mut p_reference = 0.0;
    for (g, &k_g) in alt_by_lib.iter().enumerate() {
        let p = p_alt(j, eps[g], balance);
        total += -ln_factorial(k_g) + x_ln_y(f64::from(k_g), weights[g] * p);
        p_reference += weights[g] * (1.0 - p);
    }
    total + reference_reads * p_reference.ln()
}

/// `ln P(cell | genotype j)` for the pooled key, in closed form.
///
/// The same argument collapses further: with the library forgotten, each read shows the
/// alternative allele with the share-weighted rate `Σ_g w_g·p_j(ε_g)`, so the pooled
/// count is one binomial at that rate. **No per-library depth split appears at all** —
/// which is worth noticing, because the plug-in rules invent one here.
fn ln_component_pooled(
    depth: f64,
    alt: u32,
    weights: &[f64],
    eps: &[f64],
    balance: f64,
    j: u32,
) -> f64 {
    let mut p_bar = 0.0;
    for (g, &e) in eps.iter().enumerate() {
        p_bar += weights[g] * p_alt(j, e, balance);
    }
    ln_binomial_coefficient_real(depth, alt)
        + x_ln_y(f64::from(alt), p_bar)
        + (depth - f64::from(alt)) * (1.0 - p_bar).ln()
}

/// A plug-in score: invent a per-library depth, then score as if it had been observed.
///
/// `residual` selects which invention. Everything except the reference-read term is
/// identical to `ln_component_attributed`, deliberately — the alternative reads are
/// observed per library under both, so the only thing a plug-in changes is how the
/// `n − k` reads that showed the reference are charged:
///
/// ```text
/// marginal  (n−k) · ln( Σ_g w_g·(1−p_g) )        the log of a mean
/// plug-res  (n−k) · Σ_g w_g·ln(1−p_g)            the mean of a log   — Jensen, so lower
/// plug-avg  Σ_g max(w_g·n − k_g, 0) · ln(1−p_g)  same, plus a clamp
/// ```
///
/// **Takes no heterozygote balance, unlike the three functions above it.** Those are used to
/// build the truth as well as to score it, so they have to be able to speak of a world where
/// a heterozygote is not a half. A plug-in is a scoring rule and never a generator, so it
/// assumes `MODEL_HET_BALANCE` and there is nowhere for a second value to come from.
fn ln_component_plug_in(
    depth: f64,
    alt_by_lib: &[u32],
    weights: &[f64],
    eps: &[f64],
    j: u32,
    residual: bool,
    clamped: &mut bool,
) -> f64 {
    let alt: u32 = alt_by_lib.iter().sum();
    let mut total = ln_factorial_real(depth) - ln_factorial_real(depth - f64::from(alt));
    for &k_g in alt_by_lib {
        total -= ln_factorial(k_g);
    }
    for (g, &k_g) in alt_by_lib.iter().enumerate() {
        let p = p_alt(j, eps[g], MODEL_HET_BALANCE);
        let n_hat = if residual {
            f64::from(k_g) + weights[g] * (depth - f64::from(alt))
        } else {
            weights[g] * depth
        };
        let reference_reads = n_hat - f64::from(k_g);
        if reference_reads < -1e-12 {
            *clamped = true;
        }
        total += x_ln_y(f64::from(k_g), weights[g] * p) + reference_reads.max(0.0) * (1.0 - p).ln();
    }
    total
}

/// The `GENOTYPES` component log-likelihoods of one cell under one rule, scored at
/// `depth` — the cell's own depth where the ladder is exact, and the mean of the exact
/// depths that landed in it where it is not.
///
/// **`site_noise` is the second class of site, and it wraps the whole rule rather than
/// entering it.** A site is clean with probability `1 − share` and noisy with probability
/// `share`, and the reads at a noisy site disagree with the reference at `rate` whichever
/// library produced them — so the cell's likelihood is the same rule evaluated twice, at the
/// libraries' own rates and at the noisy rate, and averaged by the share. **That is why the
/// multi-library closed form needs no rewriting**: the site's class is a property of the
/// *site* and the library split is a property of the *reads*, so the sum over which library
/// produced each alternative read happens inside each branch.
fn ln_components(
    cell: &Cell,
    depth: f64,
    rule: Rule,
    weights: &[f64],
    eps: &[f64],
    site_noise: Option<(f64, f64)>,
    clamped: &mut bool,
) -> [f64; GENOTYPES] {
    let clean = ln_components_at(cell, depth, rule, weights, eps, clamped);
    let Some((share, rate)) = site_noise else {
        return clean;
    };
    // Every library at a noisy site reads at the same rate, which is the one place the two
    // classes differ in shape: `eps` is per library and this is not.
    let noisy_eps = vec![rate; eps.len()];
    let noisy = ln_components_at(cell, depth, rule, weights, &noisy_eps, clamped);
    let mut out = [0.0; GENOTYPES];
    for (j, slot) in out.iter_mut().enumerate() {
        *slot = ln_sum_exp(&[(1.0 - share).ln() + clean[j], share.ln() + noisy[j]]);
    }
    out
}

/// The rule itself, at one class of site.
fn ln_components_at(
    cell: &Cell,
    depth: f64,
    rule: Rule,
    weights: &[f64],
    eps: &[f64],
    clamped: &mut bool,
) -> [f64; GENOTYPES] {
    let balance = MODEL_HET_BALANCE;
    let mut out = [0.0; GENOTYPES];
    for (j, slot) in out.iter_mut().enumerate() {
        let j = j as u32;
        *slot = match cell {
            // The whole breakdown is exact under every rule: nothing was thrown away.
            Cell::Whole(per_lib) => ln_component_whole(per_lib, weights, eps, balance, j),
            Cell::Attributed { alt, .. } => match rule {
                Rule::Marginal => ln_component_attributed(depth, alt, weights, eps, balance, j),
                Rule::PlugAverageShare => {
                    ln_component_plug_in(depth, alt, weights, eps, j, false, clamped)
                }
                Rule::PlugResidualShare => {
                    ln_component_plug_in(depth, alt, weights, eps, j, true, clamped)
                }
            },
            Cell::Pooled { alt, .. } => match rule {
                Rule::Marginal => ln_component_pooled(depth, *alt, weights, eps, balance, j),
                // With the library forgotten, the plug-ins split the alternative count
                // by share too — the previous comparison's own baseline.
                Rule::PlugAverageShare | Rule::PlugResidualShare => {
                    let mut total = ln_binomial_coefficient_real(depth, *alt);
                    for (g, &w) in weights.iter().enumerate() {
                        let p = p_alt(j, eps[g], balance);
                        let k_hat = w * f64::from(*alt);
                        let n_hat = w * depth;
                        total += x_ln_y(k_hat, p) + (n_hat - k_hat).max(0.0) * (1.0 - p).ln();
                    }
                    total
                }
            },
        };
    }
    out
}

// ---------------------------------------------------------------------------
// The cell space, with each cell's exact probability under the truth
// ---------------------------------------------------------------------------

/// Every cell a coarsening can produce, with the exact probability the truth puts on it.
struct CellSpace {
    cells: Vec<Cell>,
    /// `P_true(cell)`. Sums to 1 up to the depth truncation.
    mass: Vec<f64>,
    /// `Σ n · P_true(n, cell)` — the numerator of the cell's mean depth.
    depth_sum: Vec<f64>,
    index: HashMap<Cell, usize>,
}

/// Enumerate every per-library breakdown of `depth` across `libraries`, calling `visit`
/// with each `(n_g, k_g)` vector.
fn for_each_breakdown(depth: u32, libraries: usize, visit: &mut impl FnMut(&[(u32, u32)])) {
    let mut buffer = vec![(0u32, 0u32); libraries];
    fn recurse(
        g: usize,
        remaining_depth: u32,
        libraries: usize,
        buffer: &mut Vec<(u32, u32)>,
        visit: &mut impl FnMut(&[(u32, u32)]),
    ) {
        if g + 1 == libraries {
            for k in 0..=remaining_depth {
                buffer[g] = (remaining_depth, k);
                visit(buffer);
            }
            return;
        }
        for n_g in 0..=remaining_depth {
            for k in 0..=n_g {
                buffer[g] = (n_g, k);
                recurse(g + 1, remaining_depth - n_g, libraries, buffer, visit);
            }
        }
    }
    recurse(0, depth, libraries, &mut buffer, visit);
}

/// Enumerate every way `total` alternative reads split across `libraries`.
fn for_each_alt_split(total: u32, libraries: usize, visit: &mut impl FnMut(&[u32])) {
    let mut buffer = vec![0u32; libraries];
    fn recurse(
        g: usize,
        remaining: u32,
        libraries: usize,
        buffer: &mut Vec<u32>,
        visit: &mut impl FnMut(&[u32]),
    ) {
        if g + 1 == libraries {
            buffer[g] = remaining;
            visit(buffer);
            return;
        }
        for k in 0..=remaining {
            buffer[g] = k;
            recurse(g + 1, remaining - k, libraries, buffer, visit);
        }
    }
    recurse(0, total, libraries, &mut buffer, visit);
}

impl CellSpace {
    /// Build the space under `ladder`. **The truth is always exact**: each cell's
    /// probability and its mean depth are sums over the *exact* depths that map into it,
    /// scored with the world's own error rates and heterozygote balance. Only the fit is
    /// binned and only the fit assumes `MODEL_HET_BALANCE`, which is what makes the bias
    /// this program reports the bias of the estimator rather than of the generator.
    fn build(world: &World, coarsening: Coarsening, ladder: &DepthLadder) -> Self {
        let depth_probs = world.depth_distribution();
        let freqs = world.genotype_freqs();
        let libraries = world.libraries();
        let balance = world.het_balance;
        let n_max = (depth_probs.len() - 1) as u32;

        let mut cells = Vec::new();
        let mut mass = Vec::new();
        let mut depth_sum = Vec::new();
        let mut index: HashMap<Cell, usize> = HashMap::new();

        let mut push = |cell: Cell,
                        p: f64,
                        depth: u32,
                        cells: &mut Vec<Cell>,
                        mass: &mut Vec<f64>,
                        depth_sum: &mut Vec<f64>| {
            if p <= 0.0 {
                return;
            }
            match index.get(&cell) {
                Some(&i) => {
                    mass[i] += p;
                    depth_sum[i] += p * f64::from(depth);
                }
                None => {
                    index.insert(cell.clone(), cells.len());
                    cells.push(cell);
                    mass.push(p);
                    depth_sum.push(p * f64::from(depth));
                }
            }
        };

        // **The world's own second class of site, applied to the truth rather than to a
        // rule.** A site is clean with probability `1 − share` and noisy with probability
        // `share`, and at a noisy site every library reads at the same `rate` — so a cell's
        // probability under one genotype is the same closed form evaluated at two rate
        // vectors and averaged. `None` leaves every world in this file exactly as it was.
        let noisy_eps: Option<(f64, Vec<f64>)> = world
            .site_noise
            .map(|(share, rate)| (share, vec![rate; libraries]));
        let over_classes = |at: &dyn Fn(&[f64]) -> f64| -> f64 {
            let clean = at(&world.eps);
            match &noisy_eps {
                None => clean,
                Some((share, eps)) => {
                    ln_sum_exp(&[(1.0 - share).ln() + clean, share.ln() + at(eps)])
                }
            }
        };

        // Probability of a cell = P(depth) × Σ_j π_j · P(cell | depth, j).
        let mixture = |components: &[f64; GENOTYPES]| -> f64 {
            let mut terms = [0.0; GENOTYPES];
            for j in 0..GENOTYPES {
                terms[j] = freqs[j].ln() + components[j];
            }
            ln_sum_exp(&terms).exp()
        };

        for depth in 1..=n_max {
            let p_depth = depth_probs[depth as usize];
            if p_depth <= 0.0 {
                continue;
            }
            let whole_arm = coarsening.exact || depth <= coarsening.exact_below_depth;
            if whole_arm {
                for_each_breakdown(depth, libraries, &mut |per_lib| {
                    let mut components = [0.0; GENOTYPES];
                    for j in 0..GENOTYPES {
                        components[j] = over_classes(&|eps| {
                            ln_component_whole(per_lib, &world.weights, eps, balance, j as u32)
                        });
                    }
                    push(
                        Cell::Whole(per_lib.to_vec()),
                        p_depth * mixture(&components),
                        depth,
                        &mut cells,
                        &mut mass,
                        &mut depth_sum,
                    );
                });
                continue;
            }
            let bin = ladder.bin_of(depth);
            // Attributed arm: every split of every alt count up to the bound.
            for alt_total in 0..=coarsening.max_attributed.min(depth) {
                for_each_alt_split(alt_total, libraries, &mut |split| {
                    let mut components = [0.0; GENOTYPES];
                    for j in 0..GENOTYPES {
                        components[j] = over_classes(&|eps| {
                            ln_component_attributed(
                                f64::from(depth),
                                split,
                                &world.weights,
                                eps,
                                balance,
                                j as u32,
                            )
                        });
                    }
                    push(
                        Cell::Attributed {
                            bin,
                            alt: split.to_vec(),
                        },
                        p_depth * mixture(&components),
                        depth,
                        &mut cells,
                        &mut mass,
                        &mut depth_sum,
                    );
                });
            }
            // Pooled arm: everything above the attribution bound.
            for alt in (coarsening.max_attributed + 1)..=depth {
                let mut components = [0.0; GENOTYPES];
                for j in 0..GENOTYPES {
                    components[j] = over_classes(&|eps| {
                        ln_component_pooled(
                            f64::from(depth),
                            alt,
                            &world.weights,
                            eps,
                            balance,
                            j as u32,
                        )
                    });
                }
                push(
                    Cell::Pooled { bin, alt },
                    p_depth * mixture(&components),
                    depth,
                    &mut cells,
                    &mut mass,
                    &mut depth_sum,
                );
            }
        }

        Self {
            cells,
            mass,
            depth_sum,
            index,
        }
    }

    #[expect(
        dead_code,
        reason = "kept for ad-hoc runs that check the cell space sums to one"
    )]
    fn total_mass(&self) -> f64 {
        self.mass.iter().sum()
    }

    /// The depth each cell is scored at.
    ///
    /// `PerCell` is `Σ n·P(n, cell) / Σ P(n, cell)` — the infinite-genome limit of what
    /// `mean_depth_in_cell` computes from a real sample's depth sums. `PerBin` pools the
    /// alternative counts first, which is the variant the architecture doc says diverges.
    fn score_depths(&self, scoring: DepthScoring) -> Vec<f64> {
        let per_cell: Vec<f64> = self
            .mass
            .iter()
            .zip(&self.depth_sum)
            .map(|(&m, &s)| if m > 0.0 { s / m } else { 0.0 })
            .collect();
        match scoring {
            DepthScoring::PerCell => per_cell,
            DepthScoring::PerBin => {
                let mut bin_mass: HashMap<u32, (f64, f64)> = HashMap::new();
                for (i, cell) in self.cells.iter().enumerate() {
                    if let Some(bin) = cell_bin(cell) {
                        let entry = bin_mass.entry(bin).or_insert((0.0, 0.0));
                        entry.0 += self.mass[i];
                        entry.1 += self.depth_sum[i];
                    }
                }
                self.cells
                    .iter()
                    .enumerate()
                    .map(|(i, cell)| match cell_bin(cell) {
                        Some(bin) => {
                            let (m, s) = bin_mass[&bin];
                            if m > 0.0 { s / m } else { per_cell[i] }
                        }
                        None => per_cell[i],
                    })
                    .collect()
            }
        }
    }
}

/// The depth bin a cell belongs to, or `None` for the unbinned whole-breakdown arm.
fn cell_bin(cell: &Cell) -> Option<u32> {
    match cell {
        Cell::Whole(_) => None,
        Cell::Attributed { bin, .. } | Cell::Pooled { bin, .. } => Some(*bin),
    }
}

/// How many alternative reads a cell holds — the count its score depth must not fall below.
fn cell_alt(cell: &Cell) -> u32 {
    match cell {
        Cell::Whole(per_lib) => per_lib.iter().map(|&(_, k)| k).sum(),
        Cell::Attributed { alt, .. } => alt.iter().sum(),
        Cell::Pooled { alt, .. } => *alt,
    }
}

// ---------------------------------------------------------------------------
// The read-group table, and the coupled fit
// ---------------------------------------------------------------------------

/// **One library's own table**, which is the other object step 4 accumulates: a site
/// covered by two libraries enters *twice*, once per library, each entry holding only
/// that library's depth and alternative count.
///
/// The entry's marginal distribution is exact and simple. A library's reads at a site are
/// a thinning of the site's total depth, so at Poisson total depth `λ` its own depth is
/// Poisson at `λ·w_g`, truncated to at least one read because a library that covered
/// nothing enters nothing. The genotype is still drawn once for the site, so it appears
/// here through the same mixture:
///
/// ```text
/// P(m, k)  =  PoissonTruncated(λ·w_g ; m) · Σ_j π_j · Binom(k ; m, p_j(ε_g))
/// ```
///
/// **This entry's marginal is correctly specified** — nothing about splitting a site
/// between two entries makes it wrong. What is lost is that the two entries of one site
/// share a genotype, and scoring them as independent throws that dependence away. In the
/// standard vocabulary the product over entries is a **composite likelihood**: each factor
/// is a true marginal, so the truth still maximises it in expectation, and what the split
/// costs is precision rather than correctness. This section measures whether that holds.
struct ReadGroupSpace {
    library: usize,
    /// `(depth, alt reads)` for this library alone.
    cells: Vec<(u32, u32)>,
    /// The exact probability of each entry under the truth, summing to one.
    mass: Vec<f64>,
}

fn build_read_group_spaces(world: &World) -> Vec<ReadGroupSpace> {
    let freqs = world.genotype_freqs();
    (0..world.libraries())
        .map(|g| {
            let lambda = world.mean_depth * world.weights[g];
            let thinned = World {
                mean_depth: lambda,
                ..world.clone()
            };
            let depth_probs = thinned.depth_distribution();
            let mut cells = Vec::new();
            let mut mass = Vec::new();
            for m in 1..depth_probs.len() as u32 {
                for k in 0..=m {
                    // The same two classes of site as the whole-sample table: a noisy site
                    // is noisy for every library that read it, so this group's own table
                    // carries the sample's share and the sample's noisy rate.
                    let (share, noisy_rate) = world.site_noise.unwrap_or((0.0, world.eps[g]));
                    let mut p = 0.0;
                    for (j, &pi) in freqs.iter().enumerate() {
                        let at = |rate: f64| {
                            let q = p_alt(j as u32, rate, world.het_balance);
                            (ln_binomial_coefficient(m, k)
                                + x_ln_y(f64::from(k), q)
                                + x_ln_y(f64::from(m - k), 1.0 - q))
                            .exp()
                        };
                        p += pi * ((1.0 - share) * at(world.eps[g]) + share * at(noisy_rate));
                    }
                    let p = depth_probs[m as usize] * p;
                    if p > 0.0 {
                        cells.push((m, k));
                        mass.push(p);
                    }
                }
            }
            ReadGroupSpace {
                library: g,
                cells,
                mass,
            }
        })
        .collect()
}

/// `Σ over entries  P_true(entry) · ln L(entry ; ε_g, π, site noise)` — the objective the
/// error-rate fit climbs on one library's table.
///
/// `site_noise` is the **fitted** pair, and it enters exactly as it does in production's
/// per-read-group scan: every candidate rate is scored with the second class beside it.
/// Leaving it out is what made production's clean rate the one-class rate on every sample it
/// had ever seen (`reports/implementations/ng_noise_model_extension_n5_fix_2026-08-10.md`).
fn score_read_group(
    space: &ReadGroupSpace,
    eps_g: f64,
    freqs: &[f64; GENOTYPES],
    site_noise: Option<(f64, f64)>,
) -> f64 {
    let (share, noisy_rate) = site_noise.unwrap_or((0.0, eps_g));
    let mut total = 0.0;
    for (&(m, k), &w) in space.cells.iter().zip(&space.mass) {
        let mut terms = [0.0; GENOTYPES];
        for (j, slot) in terms.iter_mut().enumerate() {
            let at = |rate: f64| {
                let q = p_alt(j as u32, rate, MODEL_HET_BALANCE);
                ln_binomial_coefficient(m, k)
                    + x_ln_y(f64::from(k), q)
                    + x_ln_y(f64::from(m - k), 1.0 - q)
            };
            *slot = freqs[j].max(1e-300).ln()
                + ln_sum_exp(&[
                    (1.0 - share).max(1e-300).ln() + at(eps_g),
                    share.max(1e-300).ln() + at(noisy_rate),
                ]);
        }
        total += w * ln_sum_exp(&terms);
    }
    total
}

/// Maximise one library's error rate on its own table, with the genotype frequencies held
/// where the caller put them. Golden section in `ln ε`.
fn fit_eps_on_read_group(
    space: &ReadGroupSpace,
    freqs: &[f64; GENOTYPES],
    site_noise: Option<(f64, f64)>,
) -> f64 {
    let inverse_phi = 0.5 * (5f64.sqrt() - 1.0);
    let (mut lo, mut hi) = ((1e-7f64).ln(), (0.3f64).ln());
    let mut c = hi - inverse_phi * (hi - lo);
    let mut d = lo + inverse_phi * (hi - lo);
    let mut fc = score_read_group(space, c.exp(), freqs, site_noise);
    let mut fd = score_read_group(space, d.exp(), freqs, site_noise);
    for _ in 0..80 {
        if fc > fd {
            hi = d;
            d = c;
            fd = fc;
            c = hi - inverse_phi * (hi - lo);
            fc = score_read_group(space, c.exp(), freqs, site_noise);
        } else {
            lo = c;
            c = d;
            fc = fd;
            d = lo + inverse_phi * (hi - lo);
            fd = score_read_group(space, d.exp(), freqs, site_noise);
        }
        if (hi - lo) < 1e-9 {
            break;
        }
    }
    (0.5 * (lo + hi)).exp()
}

/// **Does alternating between the two tables land on the truth?**
///
/// The design reads the error rates off the read-group table and the genotype frequencies
/// off the whole-sample one, and resolves the coupling between them by alternating. That is
/// not coordinate ascent on any single objective — it is a fixed point of two estimating
/// equations, one per table — so whether it is even consistent has to be checked rather
/// than assumed. Both objectives are exact expectations under a known truth here, so what
/// the loop converges to is what an infinite genome would give.
fn question_coupled_fit(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## The coupled fit: alternating between two tables\n"
    );
    let _ = writeln!(
        report,
        "The error rates are read off the read-group table, the genotype frequencies off the\n\
         whole-sample one, and the two are coupled. Alternating between them is a fixed point\n\
         of two estimating equations rather than a climb on one objective, so consistency is\n\
         not automatic. Both tables are weighted by their exact probabilities under a known\n\
         truth, so a departure below is bias with no sampling noise in it. Error rates are in\n\
         rungs of the ladder; the frequencies in relative error.\n"
    );
    let _ = writeln!(
        report,
        "{:<34} {:>9} {:>9} {:>10} {:>11} {:>7}",
        "world", "ε₁ rungs", "ε₂ rungs", "π_het", "π_hom_alt", "iters"
    );

    for world in worlds() {
        if world.mean_depth > 20.0 {
            continue; // three depths are enough; the cell spaces above 20 are slow.
        }
        let truth = world.genotype_freqs();
        let ladder = DepthLadder::exact(world.max_depth());
        let whole = CellSpace::build(&world, Coarsening::attributed(4, 0), &ladder);
        let depths = whole.score_depths(DepthScoring::PerCell);
        let weights = whole.mass.clone();
        let groups = build_read_group_spaces(&world);

        // Start away from the truth on both blocks, so a fixed point at the truth is a
        // result rather than a starting condition.
        let mut eps: Vec<f64> = world.eps.iter().map(|e| e * 3.0).collect();
        let mut freqs = [
            1.0 - 0.5 * (truth[1] + truth[2]),
            0.5 * truth[1],
            0.5 * truth[2],
        ];
        let mut iterations = 0;
        for iter in 1..=200 {
            iterations = iter;
            // Block 1: the frequencies, from the whole-sample table at the current rates.
            let next_freqs = score_at(&whole, &depths, Rule::Marginal, &world, &eps, &weights).1;
            // Block 2: each library's rate, from its own table at the current frequencies.
            let next_eps: Vec<f64> = groups
                .iter()
                .map(|g| fit_eps_on_read_group(g, &next_freqs, None))
                .collect();
            let moved = next_eps
                .iter()
                .zip(&eps)
                .map(|(a, b)| (a / b).ln().abs())
                .fold(0.0f64, f64::max)
                .max(
                    next_freqs
                        .iter()
                        .zip(&freqs)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0f64, f64::max),
                );
            eps = next_eps;
            freqs = next_freqs;
            if moved < 1e-12 {
                break;
            }
        }

        let last = world.libraries() - 1;
        let _ = writeln!(
            report,
            "{:<34} {:>9.3} {:>9.3} {:>9.3}% {:>10.3}% {:>7}",
            world.name,
            rungs(eps[0] / world.eps[0]),
            rungs(eps[last] / world.eps[last]),
            100.0 * (freqs[1] / truth[1] - 1.0),
            100.0 * (freqs[2] / truth[2] - 1.0),
            iterations
        );
    }

    let _ = writeln!(
        report,
        "\nAnd the same worlds fitted **entirely on the read-group table** — both the rates and\n\
         the frequencies — which says whether splitting a site between two entries costs\n\
         anything beyond precision.\n"
    );
    let _ = writeln!(
        report,
        "{:<34} {:>9} {:>9} {:>10} {:>11} {:>7}",
        "world", "ε₁ rungs", "ε₂ rungs", "π_het", "π_hom_alt", "iters"
    );
    for world in worlds() {
        if world.mean_depth > 20.0 {
            continue;
        }
        let truth = world.genotype_freqs();
        let groups = build_read_group_spaces(&world);
        let mut eps: Vec<f64> = world.eps.iter().map(|e| e * 3.0).collect();
        let mut freqs = [
            1.0 - 0.5 * (truth[1] + truth[2]),
            0.5 * truth[1],
            0.5 * truth[2],
        ];
        let mut iterations = 0;
        for iter in 1..=200 {
            iterations = iter;
            // The frequencies, climbed on the pooled read-group evidence: every library's
            // entries are cells of one mixture with the same weights.
            let mut components: Vec<[f64; GENOTYPES]> = Vec::new();
            let mut cell_weights: Vec<f64> = Vec::new();
            for g in &groups {
                for (&(m, k), &w) in g.cells.iter().zip(&g.mass) {
                    let mut comp = [0.0; GENOTYPES];
                    for (j, slot) in comp.iter_mut().enumerate() {
                        let q = p_alt(j as u32, eps[g.library], MODEL_HET_BALANCE);
                        *slot = ln_binomial_coefficient(m, k)
                            + x_ln_y(f64::from(k), q)
                            + x_ln_y(f64::from(m - k), 1.0 - q);
                    }
                    components.push(comp);
                    cell_weights.push(w);
                }
            }
            let next_freqs = climb_frequencies(&components, &cell_weights);
            let next_eps: Vec<f64> = groups
                .iter()
                .map(|g| fit_eps_on_read_group(g, &next_freqs, None))
                .collect();
            let moved = next_eps
                .iter()
                .zip(&eps)
                .map(|(a, b)| (a / b).ln().abs())
                .fold(0.0f64, f64::max);
            eps = next_eps;
            freqs = next_freqs;
            if moved < 1e-12 {
                break;
            }
        }
        let last = world.libraries() - 1;
        let _ = writeln!(
            report,
            "{:<34} {:>9.3} {:>9.3} {:>9.3}% {:>10.3}% {:>7}",
            world.name,
            rungs(eps[0] / world.eps[0]),
            rungs(eps[last] / world.eps[last]),
            100.0 * (freqs[1] / truth[1] - 1.0),
            100.0 * (freqs[2] / truth[2] - 1.0),
            iterations
        );
    }
}

// ---------------------------------------------------------------------------
// Fitting
// ---------------------------------------------------------------------------

/// The genotype frequencies that maximise the weighted log-likelihood at fixed error
/// rates. The surface is concave in the frequencies (spec `parameter_prepass.md` §3.1),
/// so expectation-maximization reaches the maximum from any interior start.
fn climb_frequencies(components: &[[f64; GENOTYPES]], weights: &[f64]) -> [f64; GENOTYPES] {
    let total: f64 = weights.iter().sum();
    let mut freqs = [1.0 / GENOTYPES as f64; GENOTYPES];
    for _ in 0..400 {
        let mut next = [0.0; GENOTYPES];
        for (cell, &w) in components.iter().zip(weights) {
            if w == 0.0 {
                continue;
            }
            let mut terms = [0.0; GENOTYPES];
            for j in 0..GENOTYPES {
                terms[j] = freqs[j].ln() + cell[j];
            }
            let norm = ln_sum_exp(&terms);
            for j in 0..GENOTYPES {
                next[j] += w * (terms[j] - norm).exp();
            }
        }
        let mut moved: f64 = 0.0;
        for j in 0..GENOTYPES {
            next[j] /= total;
            moved = moved.max((next[j] - freqs[j]).abs());
        }
        freqs = next;
        if moved < 1e-13 {
            break;
        }
    }
    freqs
}

fn score_at(
    space: &CellSpace,
    depths: &[f64],
    rule: Rule,
    world: &World,
    eps: &[f64],
    weights: &[f64],
) -> (f64, [f64; GENOTYPES]) {
    let mut clamped = false;
    let components: Vec<[f64; GENOTYPES]> = space
        .cells
        .iter()
        .zip(depths)
        .map(|(c, &d)| ln_components(c, d, rule, &world.weights, eps, None, &mut clamped))
        .collect();
    let freqs = climb_frequencies(&components, weights);
    let mut score = 0.0;
    for (cell, &w) in components.iter().zip(weights) {
        if w == 0.0 {
            continue;
        }
        let mut terms = [0.0; GENOTYPES];
        for j in 0..GENOTYPES {
            terms[j] = freqs[j].ln() + cell[j];
        }
        score += w * ln_sum_exp(&terms);
    }
    (score, freqs)
}

/// Maximise over the error rates by coordinate-wise golden-section search in `ln ε`,
/// with the frequencies climbed to their optimum at every trial.
///
/// The search is continuous rather than on the design's quarter-Phred ladder, because
/// the bias being measured is a property of the objective and not of the search
/// resolution — quantising it to the ladder would hide anything below 6%.
fn fit(
    space: &CellSpace,
    depths: &[f64],
    rule: Rule,
    world: &World,
    weights: &[f64],
    start: &[f64],
) -> Params {
    let mut eps: Vec<f64> = start.to_vec();
    let inverse_phi = 0.5 * (5f64.sqrt() - 1.0);

    for _sweep in 0..6 {
        let mut moved: f64 = 0.0;
        for g in 0..eps.len() {
            let (mut lo, mut hi) = ((1e-7f64).ln(), (0.3f64).ln());
            let mut c = hi - inverse_phi * (hi - lo);
            let mut d = lo + inverse_phi * (hi - lo);
            let evaluate = |x: f64, eps: &mut Vec<f64>| {
                let keep = eps[g];
                eps[g] = x.exp();
                let s = score_at(space, depths, rule, world, eps, weights).0;
                eps[g] = keep;
                s
            };
            let mut fc = evaluate(c, &mut eps);
            let mut fd = evaluate(d, &mut eps);
            for _ in 0..60 {
                if fc > fd {
                    hi = d;
                    d = c;
                    fd = fc;
                    c = hi - inverse_phi * (hi - lo);
                    fc = evaluate(c, &mut eps);
                } else {
                    lo = c;
                    c = d;
                    fc = fd;
                    d = lo + inverse_phi * (hi - lo);
                    fd = evaluate(d, &mut eps);
                }
                if (hi - lo) < 1e-8 {
                    break;
                }
            }
            let best = (0.5 * (lo + hi)).exp();
            moved = moved.max((best / eps[g]).ln().abs());
            eps[g] = best;
        }
        if moved < 1e-8 {
            break;
        }
    }
    let (_, freqs) = score_at(space, depths, rule, world, &eps, weights);
    Params { eps, freqs }
}

// ---------------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------------

/// One rung of the design's error-rate ladder is a quarter-Phred, so a ratio of
/// `10^0.025` — about 5.9% in probability. Reporting the bias in rungs keeps it on the
/// scale the spec already uses to decide what a caller can feel.
fn rungs(ratio: f64) -> f64 {
    ratio.ln() / (10f64.powf(0.025)).ln()
}

struct Candidate {
    label: &'static str,
    coarsening: Coarsening,
    rule: Rule,
}

fn candidates() -> Vec<Candidate> {
    vec![
        Candidate {
            label: "exact per-library (oracle)",
            coarsening: Coarsening::oracle(),
            rule: Rule::Marginal,
        },
        Candidate {
            label: "pooled K=0 D=0  plug-avg (today)",
            coarsening: Coarsening::pooled(),
            rule: Rule::PlugAverageShare,
        },
        Candidate {
            label: "pooled K=0 D=0  marginal",
            coarsening: Coarsening::pooled(),
            rule: Rule::Marginal,
        },
        Candidate {
            label: "attrib K=4 D=0  plug-avg (proposed)",
            coarsening: Coarsening::attributed(4, 0),
            rule: Rule::PlugAverageShare,
        },
        Candidate {
            label: "attrib K=4 D=4  plug-avg (proposed)",
            coarsening: Coarsening::attributed(4, 4),
            rule: Rule::PlugAverageShare,
        },
        Candidate {
            label: "attrib K=4 D=0  plug-res",
            coarsening: Coarsening::attributed(4, 0),
            rule: Rule::PlugResidualShare,
        },
        Candidate {
            label: "attrib K=4 D=0  marginal",
            coarsening: Coarsening::attributed(4, 0),
            rule: Rule::Marginal,
        },
        Candidate {
            label: "attrib K=2 D=0  marginal",
            coarsening: Coarsening::attributed(2, 0),
            rule: Rule::Marginal,
        },
    ]
}

// ---------------------------------------------------------------------------
// The second class of site — the question `src/` cannot ask
// ---------------------------------------------------------------------------

/// One point of the profile: its score, the rates and frequencies that settled around it,
/// and the pair itself.
type ProfilePoint = (f64, Vec<f64>, [f64; GENOTYPES], (f64, f64));

/// The whole-sample objective at one point: `Σ_c P_true(c) · ln Σ_j π_j · L(c | j, θ)`.
fn score_whole_sample(
    space: &CellSpace,
    depths: &[f64],
    rule: Rule,
    world: &World,
    eps: &[f64],
    freqs: &[f64; GENOTYPES],
    site_noise: Option<(f64, f64)>,
) -> f64 {
    let mut clamped = false;
    let mut total = 0.0;
    for ((cell, &depth), &w) in space.cells.iter().zip(depths).zip(&space.mass) {
        if w == 0.0 {
            continue;
        }
        let comp = ln_components(
            cell,
            depth,
            rule,
            &world.weights,
            eps,
            site_noise,
            &mut clamped,
        );
        let mut terms = [0.0; GENOTYPES];
        for j in 0..GENOTYPES {
            terms[j] = freqs[j].max(1e-300).ln() + comp[j];
        }
        total += w * ln_sum_exp(&terms);
    }
    total
}

/// The share of noisy sites that best explains the table, with everything else held — the
/// same expectation-maximisation production climbs, and concave in the share for a fixed
/// noisy rate because both branches are then constants.
fn climb_noisy_share(
    space: &CellSpace,
    depths: &[f64],
    rule: Rule,
    world: &World,
    eps: &[f64],
    freqs: &[f64; GENOTYPES],
    noisy_rate: f64,
) -> f64 {
    let mut clamped = false;
    let noisy_eps = vec![noisy_rate; eps.len()];
    // The two branches, marginalised over the genotypes once: neither moves as the share does.
    let mut branches = Vec::with_capacity(space.cells.len());
    for (cell, &depth) in space.cells.iter().zip(depths) {
        let at = |rates: &[f64], clamped: &mut bool| {
            let comp = ln_components_at(cell, depth, rule, &world.weights, rates, clamped);
            let mut terms = [0.0; GENOTYPES];
            for j in 0..GENOTYPES {
                terms[j] = freqs[j].max(1e-300).ln() + comp[j];
            }
            ln_sum_exp(&terms)
        };
        branches.push((at(eps, &mut clamped), at(&noisy_eps, &mut clamped)));
    }

    let total: f64 = space.mass.iter().sum();
    let mut share: f64 = 0.5;
    for _ in 0..400 {
        let mut responsibility = 0.0;
        for ((clean, noisy), &w) in branches.iter().zip(&space.mass) {
            let (from_clean, from_noisy) = (
                (1.0 - share).max(1e-300).ln() + clean,
                share.max(1e-300).ln() + noisy,
            );
            let larger = from_clean.max(from_noisy);
            let (a, b) = ((from_clean - larger).exp(), (from_noisy - larger).exp());
            responsibility += w * b / (a + b);
        }
        let next = responsibility / total;
        let settled = (next - share).abs() < 1e-14;
        share = next;
        if settled {
            break;
        }
    }
    share
}

/// **The coupled fit with a second class of site, fitted the way production fits it.** The
/// noisy rate is *profiled* — held at each point of a grid while the rates, the frequencies
/// and the share are refitted around it — and the best-scoring point wins.
///
/// **Why the profile and not a third block in the alternation**, and it is production's own
/// measurement rather than a preference: fitting the pair once at whatever the previous round
/// settled on leaves the answer at a point that is optimal in each block and not jointly. On
/// tomato SRR7279481 that cost 209 nats
/// (`reports/implementations/ng_noise_model_extension_n5_fix_2026-08-10.md`).
///
/// The grid is continuous rather than the design's quarter-Phred ladder, for this file's
/// standing reason: the bias being measured is a property of the objective, and quantising
/// the search would hide anything below a rung.
fn fit_coupled_with_site_noise(
    whole: &CellSpace,
    depths: &[f64],
    rule: Rule,
    world: &World,
    groups: &[ReadGroupSpace],
    start_eps: &[f64],
    start_freqs: [f64; GENOTYPES],
) -> (Vec<f64>, [f64; GENOTYPES], Option<(f64, f64)>) {
    // **Started where the caller says and then warm-started**, which is what makes this
    // affordable: the profile refits everything at every point of its grid, and from a cold
    // start at three times the truth that is a full alternation each time. Every point after
    // the first begins at the previous point's answer, and the rate step inside is an
    // exhaustive search rather than a local one, so a warm start cannot pin an answer — only
    // reach it sooner. Measured: the whole question goes from over 40 minutes to under two.
    let settle = |site_noise: Option<(f64, f64)>,
                  from_eps: &[f64],
                  from_freqs: [f64; GENOTYPES]|
     -> (Vec<f64>, [f64; GENOTYPES]) {
        let mut eps = from_eps.to_vec();
        let mut freqs = from_freqs;
        for _ in 0..60 {
            let next_freqs = {
                let mut clamped = false;
                let components: Vec<[f64; GENOTYPES]> = whole
                    .cells
                    .iter()
                    .zip(depths)
                    .map(|(c, &d)| {
                        ln_components(c, d, rule, &world.weights, &eps, site_noise, &mut clamped)
                    })
                    .collect();
                climb_frequencies(&components, &whole.mass)
            };
            let next_eps: Vec<f64> = groups
                .iter()
                .map(|g| fit_eps_on_read_group(g, &next_freqs, site_noise))
                .collect();
            let moved = next_eps
                .iter()
                .zip(&eps)
                .map(|(a, b)| (a / b).ln().abs())
                .fold(0.0f64, f64::max)
                .max(
                    next_freqs
                        .iter()
                        .zip(&freqs)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0f64, f64::max),
                );
            eps = next_eps;
            freqs = next_freqs;
            if moved < 1e-12 {
                break;
            }
        }
        (eps, freqs)
    };

    // The control: one class of site, which is what a world without a second class must come
    // back to and the score every candidate has to beat.
    let (one_eps, one_freqs) = settle(None, start_eps, start_freqs);
    let one_score = score_whole_sample(whole, depths, rule, world, &one_eps, &one_freqs, None);

    // One point of the profile: the noisy rate is held and everything else refitted around
    // it. Returns `None` where the share collapses, which is this rate buying no class at all.
    let mut warm_eps = one_eps.clone();
    let mut warm_freqs = one_freqs;
    let at_rate = |rate: f64, warm_eps: &mut Vec<f64>, warm_freqs: &mut [f64; GENOTYPES]| {
        let (mut eps, mut freqs) = (warm_eps.clone(), *warm_freqs);
        let mut share = 0.0;
        for _ in 0..12 {
            let next = climb_noisy_share(whole, depths, rule, world, &eps, &freqs, rate);
            let settled = (next - share).abs() < 1e-13;
            share = next;
            let (e, f) = settle(Some((share, rate)), &eps, freqs);
            eps = e;
            freqs = f;
            if settled {
                break;
            }
        }
        *warm_eps = eps.clone();
        *warm_freqs = freqs;
        if share <= 0.0 {
            return None;
        }
        let score = score_whole_sample(
            whole,
            depths,
            rule,
            world,
            &eps,
            &freqs,
            Some((share, rate)),
        );
        Some((score, eps, freqs, (share, rate)))
    };

    // **A grid, then a continuous refinement between the winner's neighbours.** The grid alone
    // is not enough and the reason is arithmetic rather than judgement: 25 points spanning
    // 10⁻³ to 0.3 step by a factor of 1.27, which is **4.3 rungs of the ladder**, so the noisy
    // rate could only ever come back within about two rungs of the truth — and it did, at 1.0
    // to 2.3 rungs, which reads as bias and is nothing but the spacing. Golden section in
    // `ln rate` between the neighbours of the best point removes the grid from the answer, for
    // the same reason every other search in this file is continuous: what is being measured is
    // a property of the objective, and quantising the search hides anything finer than a step.
    const STEPS: usize = 24;
    let rate_at = |step: f64| (1e-3f64.ln() + (0.3f64.ln() - 1e-3f64.ln()) * step / 24.0).exp();
    let mut best: Option<ProfilePoint> = None;
    let mut best_step = 0usize;
    for step in 0..=STEPS {
        let Some(point) = at_rate(rate_at(step as f64), &mut warm_eps, &mut warm_freqs) else {
            continue;
        };
        if best.as_ref().is_none_or(|(kept, ..)| point.0 > *kept) {
            best = Some(point);
            best_step = step;
        }
    }

    if best.is_some() {
        let inverse_phi = 0.5 * (5f64.sqrt() - 1.0);
        let (mut lo, mut hi) = (
            rate_at(best_step.saturating_sub(1) as f64).ln(),
            rate_at(best_step.saturating_add(1).min(STEPS) as f64).ln(),
        );
        let mut c = hi - inverse_phi * (hi - lo);
        let mut d = lo + inverse_phi * (hi - lo);
        let mut point_c = at_rate(c.exp(), &mut warm_eps, &mut warm_freqs);
        let mut point_d = at_rate(d.exp(), &mut warm_eps, &mut warm_freqs);
        for _ in 0..40 {
            let (score_c, score_d) = (
                point_c.as_ref().map_or(f64::NEG_INFINITY, |p| p.0),
                point_d.as_ref().map_or(f64::NEG_INFINITY, |p| p.0),
            );
            if score_c > score_d {
                hi = d;
                d = c;
                point_d = point_c;
                c = hi - inverse_phi * (hi - lo);
                point_c = at_rate(c.exp(), &mut warm_eps, &mut warm_freqs);
            } else {
                lo = c;
                c = d;
                point_c = point_d;
                d = lo + inverse_phi * (hi - lo);
                point_d = at_rate(d.exp(), &mut warm_eps, &mut warm_freqs);
            }
            if (hi - lo) < 1e-9 {
                break;
            }
        }
        for point in [point_c, point_d].into_iter().flatten() {
            if best.as_ref().is_none_or(|(kept, ..)| point.0 > *kept) {
                best = Some(point);
            }
        }
    }

    match best {
        // **Production's three-nat floor does not belong here, and putting it here made every
        // two-class sample decline its own second class.** That floor is a likelihood-ratio
        // test on a real table, whose weights are counts of sites — χ²(2) at p ≈ 0.05 is
        // 5.99, a gain of 3 nats over half a million sites. The weights in this file are
        // *probabilities summing to one*, so a nat here is a nat **per site**, and three of
        // them is a threshold five orders of magnitude above any gain a real second class
        // buys. What this file measures is the argmax at infinite data, where a floor has no
        // meaning at all: the only thing worth rejecting is arithmetic noise.
        Some((score, eps, freqs, pair)) if score - one_score > 1e-12 => (eps, freqs, Some(pair)),
        _ => (one_eps, one_freqs, None),
    }
}

/// The worlds that have a second class of site, at the shares and rates HG002 and the tomato
/// samples actually returned — crossed with one and two libraries, because **the second class
/// is per sample while `ε` is per library, and no fixture in `src/` can ask what that costs a
/// key that has thrown the library away.**
fn two_class_worlds() -> Vec<World> {
    let mut out = Vec::new();
    for &(share, noisy, label) in &[
        (0.0088, 5.31e-2, "hg002-30x"),
        (0.0127, 4.22e-2, "hg002-300x"),
        (0.0142, 6.31e-2, "tomato"),
    ] {
        for &(libraries, split_name) in &[(1usize, "one library"), (2, "two libraries")] {
            let (eps, weights) = if libraries == 1 {
                (vec![1.9e-3], vec![1.0])
            } else {
                (vec![1.9e-3, 7.6e-3], vec![0.5, 0.5])
            };
            out.push(World {
                name: format!("{label} {split_name}"),
                eps,
                weights,
                mean_depth: 10.0,
                pi_het: 1e-3,
                pi_hom_alt: 6e-4,
                het_balance: MODEL_HET_BALANCE,
                site_noise: Some((share, noisy)),
            });
        }
    }
    out
}

/// **Is the second class of site recovered when the cell key has thrown away which library
/// produced each alternative read?**
///
/// Everything else about this milestone was decided on single-library evidence: both cohorts
/// carry one read group, and the research note's two synthetic worlds have one rate. The
/// question that leaves open is this one, and it is the question this harness exists for —
/// the cells are weighted by their exact probability under a known truth, so what the fit
/// converges to is what an infinite genome would give and a departure is bias with no
/// sampling noise in it.
///
/// **A control that must return nothing**: the same worlds with the second class removed.
fn question_the_second_class_of_site(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## The second class of site, under a coarsened key\n"
    );
    let _ = writeln!(
        report,
        "A site is clean with probability `1 − w` and noisy with probability `w`, and at a\n\
         noisy site every library reads at `ε_noisy` — one rate for the sample, where `ε` is\n\
         one rate per library. The fit is handed none of the three. Bias is exact: every cell\n\
         carries its probability under the truth, so a non-zero entry is bias and not noise.\n\
         Error rates are in rungs of the ladder; `w` and the frequencies in relative error.\n"
    );
    let _ = writeln!(
        report,
        "{:<28} {:<22} {:>9} {:>9} {:>9} {:>9} {:>10}",
        "world", "key", "ε₁ rungs", "ε_noisy", "w", "π_het", "π_hom_alt"
    );

    for world in two_class_worlds() {
        let truth = world.genotype_freqs();
        let (true_share, true_noisy) = world.site_noise.expect("a two-class world");
        let ladder = DepthLadder::exact(world.max_depth());
        let groups = build_read_group_spaces(&world);
        let start_eps: Vec<f64> = world.eps.iter().map(|e| e * 3.0).collect();
        let start_freqs = [
            1.0 - 0.5 * (truth[1] + truth[2]),
            0.5 * truth[1],
            0.5 * truth[2],
        ];

        for (label, coarsening) in [
            ("pooled K=0 D=0", Coarsening::pooled()),
            ("attributed K=4 D=0", Coarsening::attributed(4, 0)),
        ] {
            let space = CellSpace::build(&world, coarsening, &ladder);
            let depths = space.score_depths(DepthScoring::PerCell);
            let (eps, freqs, pair) = fit_coupled_with_site_noise(
                &space,
                &depths,
                Rule::Marginal,
                &world,
                &groups,
                &start_eps,
                start_freqs,
            );
            match pair {
                Some((share, noisy)) => {
                    let _ = writeln!(
                        report,
                        "{:<28} {:<22} {:>9.3} {:>9.3} {:>8.2}% {:>8.2}% {:>9.2}%",
                        world.name,
                        label,
                        rungs(eps[0] / world.eps[0]),
                        rungs(noisy / true_noisy),
                        100.0 * (share / true_share - 1.0),
                        100.0 * (freqs[1] / truth[1] - 1.0),
                        100.0 * (freqs[2] / truth[2] - 1.0),
                    );
                }
                None => {
                    let _ = writeln!(
                        report,
                        "{:<28} {:<22} {:>9.3} {:>9} {:>9} {:>9.2} {:>10.2}",
                        world.name,
                        label,
                        rungs(eps[0] / world.eps[0]),
                        "DECLINED",
                        "-",
                        100.0 * (freqs[1] / truth[1] - 1.0),
                        100.0 * (freqs[2] / truth[2] - 1.0),
                    );
                }
            }
        }
    }

    let _ = writeln!(
        report,
        "\nAnd the control — the same worlds with the second class taken out of the truth.\n\
         Every row must decline it: the two-class model contains the one-class model, so on\n\
         a world generated by one rate the maximum is the truth and the extra pair buys\n\
         nothing above the three-nat floor.\n"
    );
    let _ = writeln!(
        report,
        "{:<28} {:<22} {:>9} {:>12} {:>9} {:>10}",
        "world", "key", "ε₁ rungs", "second class", "π_het", "π_hom_alt"
    );
    for world in two_class_worlds() {
        let control = World {
            site_noise: None,
            ..world.clone()
        };
        let truth = control.genotype_freqs();
        let ladder = DepthLadder::exact(control.max_depth());
        let groups = build_read_group_spaces(&control);
        let start_eps: Vec<f64> = control.eps.iter().map(|e| e * 3.0).collect();
        let start_freqs = [
            1.0 - 0.5 * (truth[1] + truth[2]),
            0.5 * truth[1],
            0.5 * truth[2],
        ];
        let space = CellSpace::build(&control, Coarsening::pooled(), &ladder);
        let depths = space.score_depths(DepthScoring::PerCell);
        let (eps, freqs, pair) = fit_coupled_with_site_noise(
            &space,
            &depths,
            Rule::Marginal,
            &control,
            &groups,
            &start_eps,
            start_freqs,
        );
        let _ = writeln!(
            report,
            "{:<28} {:<22} {:>9.3} {:>12} {:>8.3}% {:>9.3}%",
            control.name,
            "pooled K=0 D=0",
            rungs(eps[0] / control.eps[0]),
            match pair {
                None => "declined".to_string(),
                Some((share, rate)) => format!("{:.3}% at {rate:.1e}", 100.0 * share),
            },
            100.0 * (freqs[1] / truth[1] - 1.0),
            100.0 * (freqs[2] / truth[2] - 1.0),
        );
    }
}

fn worlds() -> Vec<World> {
    let mut out = Vec::new();
    for &(ratio, ratio_name) in &[(1.0, "1"), (4.0, "4"), (10.0, "10")] {
        // Depths 6 and 10 are the band the design's two mechanisms both miss: too deep
        // for the whole-breakdown arm (`D = 4`) and too shallow for the average-share
        // depth to exceed four alternative reads in the minor library.
        for &depth in &[3.0f64, 6.0, 10.0, 20.0, 60.0] {
            for &(split_name, w0) in &[("even", 0.5), ("skew90", 0.9)] {
                out.push(World {
                    name: format!("ratio={ratio_name} depth={depth:.0} split={split_name}"),
                    eps: vec![1e-3, 1e-3 * ratio],
                    weights: vec![w0, 1.0 - w0],
                    mean_depth: depth,
                    pi_het: 1e-2,
                    pi_hom_alt: 6e-3,
                    het_balance: MODEL_HET_BALANCE,
                    site_noise: None,
                });
            }
        }
    }
    out.push(World {
        name: "four libraries depth=20".to_string(),
        eps: vec![5e-4, 1e-3, 2e-3, 4e-3],
        weights: vec![0.25, 0.25, 0.25, 0.25],
        mean_depth: 20.0,
        pi_het: 1e-2,
        pi_hom_alt: 6e-3,
        het_balance: MODEL_HET_BALANCE,
        site_noise: None,
    });
    out
}

// ---------------------------------------------------------------------------
// The checks that need no simulation
// ---------------------------------------------------------------------------

fn algebraic_checks(world: &World) -> String {
    let mut report = String::new();
    let ladder = DepthLadder::exact(world.max_depth());
    let coarsenings: [(&str, Coarsening); 3] = [
        ("pooled K=0 D=0", Coarsening::pooled()),
        ("attrib K=4 D=0", Coarsening::attributed(4, 0)),
        ("attrib K=4 D=4", Coarsening::attributed(4, 4)),
    ];

    for (name, coarsening) in coarsenings {
        let space = CellSpace::build(world, coarsening, &ladder);
        let depths = space.score_depths(DepthScoring::PerCell);

        // 1. Proper likelihood: does the rule sum to one over the cell space?
        for rule in [
            Rule::Marginal,
            Rule::PlugResidualShare,
            Rule::PlugAverageShare,
        ] {
            let mut clamped = false;
            let mut total = 0.0;
            let depth_probs = world.depth_distribution();
            for (cell, &depth) in space.cells.iter().zip(&depths) {
                let comp = ln_components(
                    cell,
                    depth,
                    rule,
                    &world.weights,
                    &world.eps,
                    None,
                    &mut clamped,
                );
                let mut terms = [0.0; GENOTYPES];
                let freqs = world.genotype_freqs();
                for j in 0..GENOTYPES {
                    terms[j] = freqs[j].ln() + comp[j];
                }
                total += depth_probs[depth.round() as usize] * ln_sum_exp(&terms).exp();
            }
            let _ = writeln!(
                report,
                "  {name:<15} {:<9} sums to {total:.9}   {}",
                rule.label(),
                if (total - 1.0).abs() < 1e-9 {
                    "PASS"
                } else {
                    "FAIL"
                }
            );
        }

        // 2. Non-negative exponents: how much truth mass sits on a clamped cell?
        let mut clamped_mass = 0.0;
        for ((cell, &m), &depth) in space.cells.iter().zip(&space.mass).zip(&depths) {
            if let Cell::Attributed { alt, .. } = cell {
                for (g, &k_g) in alt.iter().enumerate() {
                    if world.weights[g] * depth < f64::from(k_g) - 1e-12 {
                        clamped_mass += m;
                        break;
                    }
                }
            }
        }
        let _ = writeln!(
            report,
            "  {name:<15} plug-avg  clamp fires on {:.4}% of sites   {}",
            100.0 * clamped_mass,
            if clamped_mass < 1e-12 { "PASS" } else { "FAIL" }
        );

        // 3. Equal error rates: every rule must reproduce the exact likelihood.
        //
        // Compared after subtracting each cell's own first component, so that a
        // per-cell factor common to every genotype — which cancels out of the mixture
        // and cannot move the fit — is not counted as a difference. What survives is
        // exactly what the estimator can see.
        let mut null = world.clone();
        let flat = null.eps[0];
        null.eps = vec![flat; null.libraries()];
        let null_space = CellSpace::build(&null, coarsening, &ladder);
        let null_depths = null_space.score_depths(DepthScoring::PerCell);
        for rule in [
            Rule::Marginal,
            Rule::PlugResidualShare,
            Rule::PlugAverageShare,
        ] {
            let mut clamped = false;
            let mut worst: f64 = 0.0;
            for (cell, &depth) in null_space.cells.iter().zip(&null_depths) {
                let got = ln_components(
                    cell,
                    depth,
                    rule,
                    &null.weights,
                    &null.eps,
                    None,
                    &mut clamped,
                );
                let want = ln_components(
                    cell,
                    depth,
                    Rule::Marginal,
                    &null.weights,
                    &null.eps,
                    None,
                    &mut clamped,
                );
                for j in 1..GENOTYPES {
                    let d = (got[j] - got[0]) - (want[j] - want[0]);
                    if d.is_finite() {
                        worst = worst.max(d.abs());
                    }
                }
            }
            let _ = writeln!(
                report,
                "  {name:<15} {:<9} equal-ε max |Δ ln L| = {worst:.3e}   {}",
                rule.label(),
                if worst < 1e-9 { "PASS" } else { "FAIL" }
            );
        }
    }

    // 4. One library: every rule collapses to the pooled key.
    let mut single = world.clone();
    single.eps = vec![world.eps[0]];
    single.weights = vec![1.0];
    let space = CellSpace::build(&single, Coarsening::attributed(4, 4), &ladder);
    let single_depths = space.score_depths(DepthScoring::PerCell);
    let mut worst: f64 = 0.0;
    let mut clamped = false;
    for (cell, &depth) in space.cells.iter().zip(&single_depths) {
        let alt = cell_alt(cell);
        let mut want = [0.0; GENOTYPES];
        for j in 0..GENOTYPES {
            want[j] = ln_component_pooled(
                depth,
                alt,
                &single.weights,
                &single.eps,
                MODEL_HET_BALANCE,
                j as u32,
            );
        }
        for rule in [
            Rule::Marginal,
            Rule::PlugResidualShare,
            Rule::PlugAverageShare,
        ] {
            let got = ln_components(
                cell,
                depth,
                rule,
                &single.weights,
                &single.eps,
                None,
                &mut clamped,
            );
            for j in 1..GENOTYPES {
                let d = (got[j] - got[0]) - (want[j] - want[0]);
                worst = worst.max(d.abs());
            }
        }
    }
    let _ = writeln!(
        report,
        "  one library     all rules collapse to pooled, max |Δ ln L| = {worst:.3e}   {}",
        if worst < 1e-9 { "PASS" } else { "FAIL" }
    );

    report
}

// ---------------------------------------------------------------------------
// Experiment 1 — what depth binning costs
// ---------------------------------------------------------------------------

/// The candidate ladders. The first is the control: one bin per depth, so every cell's
/// mean depth is its own depth and the binned machinery must reproduce the unbinned
/// answer exactly. A bias it reported would be the harness's, not the binning's.
fn ladders(max_depth: u32) -> Vec<DepthLadder> {
    vec![
        DepthLadder::exact(max_depth),
        DepthLadder::geometric(4, 16, 124),
        DepthLadder::geometric(8, 16, 124),
        DepthLadder::geometric(16, 20, 124),
        DepthLadder::geometric(8, 20, 124),
        DepthLadder::geometric(12, 20, 124),
        DepthLadder::geometric(8, 24, 124),
        DepthLadder::geometric(8, 16, 300),
        DepthLadder::geometric(8, 20, 300),
    ]
}

/// The worlds the binning question is asked in. Depth 3 is tomato, 10 and 30 bracket an
/// ordinary whole-genome run, and 60 is a deep one. **No world reaches any ladder's cap**:
/// at mean depth 60 the deepest site the Poisson truncation keeps is 125, against caps of
/// 124 and 300. So the cap enters these numbers only through how far it stretches the
/// geometric region — what a site *deeper* than the cap costs is the subsampling rule, and
/// this harness does not implement one.
fn binning_worlds() -> Vec<World> {
    let mut out = Vec::new();
    for &(ratio, ratio_name) in &[(1.0, "1"), (4.0, "4")] {
        for &depth in &[3.0f64, 10.0, 20.0, 30.0, 60.0] {
            for &(split_name, w0) in &[("even", 0.5), ("skew90", 0.9)] {
                out.push(World {
                    name: format!("ratio={ratio_name} depth={depth:.0} {split_name}"),
                    eps: vec![1e-3, 1e-3 * ratio],
                    weights: vec![w0, 1.0 - w0],
                    mean_depth: depth,
                    pi_het: 1e-2,
                    pi_hom_alt: 6e-3,
                    het_balance: MODEL_HET_BALANCE,
                    site_noise: None,
                });
            }
        }
    }
    out
}

/// The checks that reject a broken binning before any fit is run.
fn binning_checks(report: &mut String, world: &World) {
    let _ = writeln!(report, "\n## Depth binning — the checks that need no fit\n");
    let _ = writeln!(
        report,
        "On {}. The two `moves` columns are how far a site's own contribution to the objective\n\
         moves when its exact depth is replaced by its cell's mean — the truth-mass-weighted\n\
         mean of |Δ| in nats per site, split by whether the site showed any alternative read,\n\
         and **exactly zero on the exact ladder**, which is the control. They are taken on the\n\
         mixture over genotypes rather than on its components, because a component that shifts\n\
         by hundreds of nats while sitting far below its neighbours shifts nothing the fit can\n\
         see. `n̄ < k` is the truth mass on cells scored at a depth below their own alternative\n\
         count: those cells charge a **negative number of reference reads**, and their term\n\
         `(ε/3)^(n−k)` grows as ε falls. `Σ L` is what the rule sums to over the cell space.\n",
        world.name
    );
    let _ = writeln!(
        report,
        "{:<26} {:>7} {:>9} {:>10} {:>10} {:>13} {:>13} {:>10} {:>10}",
        "ladder",
        "bins",
        "widest",
        "moves k=0",
        "moves k≥1",
        "n̄<k per-cell",
        "n̄<k per-bin",
        "Σ L cell",
        "cells"
    );

    let exact_ladder = DepthLadder::exact(world.max_depth());
    let exact_space = CellSpace::build(world, Coarsening::attributed(4, 0), &exact_ladder);
    let exact_depths = exact_space.score_depths(DepthScoring::PerCell);
    let mut clamped = false;
    let exact_components: HashMap<Cell, [f64; GENOTYPES]> = exact_space
        .cells
        .iter()
        .zip(&exact_depths)
        .map(|(c, &d)| {
            (
                c.clone(),
                ln_components(
                    c,
                    d,
                    Rule::Marginal,
                    &world.weights,
                    &world.eps,
                    None,
                    &mut clamped,
                ),
            )
        })
        .collect();

    for ladder in ladders(world.max_depth()) {
        let space = CellSpace::build(world, Coarsening::attributed(4, 0), &ladder);
        let per_cell = space.score_depths(DepthScoring::PerCell);
        let per_bin = space.score_depths(DepthScoring::PerBin);

        // 1. How far does a site's score move? Compare each *exact* cell's score against
        //    the score its binned cell now receives, after subtracting the
        //    genotype-independent part — what survives is exactly what the estimator can
        //    see — and weight by how often such a site occurs.
        let mut moved = [0.0f64; 2];
        let mut moved_mass = [0.0f64; 2];
        for depth in 1..=world.max_depth() {
            let bin = ladder.bin_of(depth);
            for alt_total in 0..=4u32.min(depth) {
                for_each_alt_split(alt_total, world.libraries(), &mut |split| {
                    let binned_cell = Cell::Attributed {
                        bin,
                        alt: split.to_vec(),
                    };
                    let exact_cell = Cell::Attributed {
                        bin: depth,
                        alt: split.to_vec(),
                    };
                    let (Some(&i), Some(&j_exact), Some(want)) = (
                        space.index.get(&binned_cell),
                        exact_space.index.get(&exact_cell),
                        exact_components.get(&exact_cell),
                    ) else {
                        return;
                    };
                    let got = ln_components(
                        &binned_cell,
                        per_cell[i],
                        Rule::Marginal,
                        &world.weights,
                        &world.eps,
                        None,
                        &mut clamped,
                    );
                    // What the fit actually sums is the *mixture* over genotypes, and each
                    // cell's own prefactor cancels out of it, so both are taken relative to
                    // the homozygous-reference component. A component that moves by
                    // hundreds of nats while sitting 10⁻⁴⁰⁰ below its neighbours moves
                    // nothing the estimator can see.
                    let freqs = world.genotype_freqs();
                    let mixture = |c: &[f64; GENOTYPES]| -> f64 {
                        let terms: Vec<f64> = (0..GENOTYPES)
                            .map(|j| freqs[j].ln() + c[j] - c[0])
                            .collect();
                        ln_sum_exp(&terms)
                    };
                    let d = mixture(&got) - mixture(want);
                    if d.is_finite() {
                        let mass = exact_space.mass[j_exact];
                        let slot = if alt_total == 0 { 0 } else { 1 };
                        moved[slot] += mass * d.abs();
                        moved_mass[slot] += mass;
                    }
                });
            }
        }

        // 2. Cells scored below their own alternative count, under each depth rule.
        let mass_below = |depths: &[f64]| -> f64 {
            space
                .cells
                .iter()
                .zip(&space.mass)
                .zip(depths)
                .filter(|((cell, _), d)| **d + 1e-12 < f64::from(cell_alt(cell)))
                .map(|((_, &m), _)| m)
                .sum::<f64>()
        };

        // 3. Sum over the cell space, at the truth.
        let freqs = world.genotype_freqs();
        let mut total = 0.0;
        for (i, cell) in space.cells.iter().enumerate() {
            let comp = ln_components(
                cell,
                per_cell[i],
                Rule::Marginal,
                &world.weights,
                &world.eps,
                None,
                &mut clamped,
            );
            let mut terms = [0.0; GENOTYPES];
            for j in 0..GENOTYPES {
                terms[j] = freqs[j].ln() + comp[j];
            }
            // Each binned cell carries the mass of every exact depth that maps into it, so
            // the depth factor is that bin's total probability rather than one depth's.
            let depth_probs = world.depth_distribution();
            let bin = cell_bin(cell).unwrap_or(0);
            let bin_mass: f64 = (1..depth_probs.len())
                .filter(|&d| ladder.bin_of(d as u32) == bin)
                .map(|d| depth_probs[d])
                .sum();
            total += bin_mass * ln_sum_exp(&terms).exp();
        }

        let mean_move = |i: usize| {
            if moved_mass[i] > 0.0 {
                moved[i] / moved_mass[i]
            } else {
                0.0
            }
        };
        let _ = writeln!(
            report,
            "{:<26} {:>7} {:>9} {:>10.2e} {:>10.2e} {:>12.4}% {:>12.4}% {:>10.4} {:>10}",
            ladder.label,
            ladder.bins,
            ladder.widest_bin(),
            mean_move(0),
            mean_move(1),
            100.0 * mass_below(&per_cell),
            100.0 * mass_below(&per_bin),
            total,
            space.cells.len()
        );
    }
    let _ = writeln!(
        report,
        "\n`Σ L cell` is not expected to be exactly one away from the exact ladder: a binned\n\
         cell's score is one binomial at a fractional depth, and summing those over a bin's\n\
         alternative counts is not the bin's probability. Its distance from one is the size of\n\
         the approximation, not a pass/fail."
    );
}

/// **Experiment 1 (a) and (b): the asymptotic bias of the binned score.**
fn question_depth_binning(report: &mut String) {
    let _ = writeln!(report, "\n## Depth binning — asymptotic bias\n");
    let _ = writeln!(
        report,
        "Every cell is weighted by its exact probability under the truth and scored at the\n\
         exact mean of the depths that landed in it, `Σ n·P(n, cell) / Σ P(n, cell)` — the\n\
         infinite-genome limit of what `mean_depth_in_cell` computes from a sample's depth\n\
         sums. So a departure below is bias with no sampling noise in it. Error rates in rungs\n\
         of the quarter-Phred ladder (~5.9% each); the two frequencies in relative error.\n\
         `spread` is how far the answer for the second library's rate moves when the search\n\
         starts at 3× and at ⅓ of the truth rather than at it.\n"
    );

    for world in binning_worlds() {
        let truth = world.genotype_freqs();
        let true_mean_eps: f64 = world
            .eps
            .iter()
            .zip(&world.weights)
            .map(|(e, w)| e * w)
            .sum();
        let last = world.libraries() - 1;
        let _ = writeln!(
            report,
            "\n### {}   ε = {:?}, shares = {:?}",
            world.name, world.eps, world.weights
        );
        let _ = writeln!(
            report,
            "{:<26} {:>8} {:>9} {:>9} {:>8} {:>9} {:>11} {:>7}",
            "ladder", "ε̄", "ε first", "ε last", "spread", "π_het", "π_hom_alt", "cells"
        );
        for ladder in ladders(world.max_depth()) {
            let space = CellSpace::build(&world, Coarsening::attributed(4, 0), &ladder);
            let depths = space.score_depths(DepthScoring::PerCell);
            let weights = space.mass.clone();
            let mut fits = Vec::new();
            for scale in [1.0, 3.0, 1.0 / 3.0] {
                let start: Vec<f64> = world.eps.iter().map(|e| e * scale).collect();
                fits.push(fit(
                    &space,
                    &depths,
                    Rule::Marginal,
                    &world,
                    &weights,
                    &start,
                ));
            }
            let fitted = &fits[0];
            let spread = fits
                .iter()
                .map(|f| rungs(f.eps[last] / world.eps[last]))
                .fold(f64::NEG_INFINITY, f64::max)
                - fits
                    .iter()
                    .map(|f| rungs(f.eps[last] / world.eps[last]))
                    .fold(f64::INFINITY, f64::min);
            let fitted_mean_eps: f64 = fitted
                .eps
                .iter()
                .zip(&world.weights)
                .map(|(e, w)| e * w)
                .sum();
            let _ = writeln!(
                report,
                "{:<26} {:>8.3} {:>9.3} {:>9.3} {:>8.3} {:>8.3}% {:>10.3}% {:>7}",
                ladder.label,
                rungs(fitted_mean_eps / true_mean_eps),
                rungs(fitted.eps[0] / world.eps[0]),
                rungs(fitted.eps[last] / world.eps[last]),
                spread,
                100.0 * (fitted.freqs[1] / truth[1] - 1.0),
                100.0 * (fitted.freqs[2] / truth[2] - 1.0),
                space.cells.len(),
            );
        }
    }
}

/// **The per-bin mean depth**, which the architecture doc says diverges. Measured rather
/// than asserted, on the worlds where a bin is widest.
fn question_per_bin_mean_depth(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## Scoring a cell at its bin's mean depth instead of its own\n"
    );
    let _ = writeln!(
        report,
        "The architecture doc (§2.2) argues that a **bin** mean is not merely coarser but\n\
         unbounded: a cell whose alternative count exceeds its bin's mean depth is charged a\n\
         negative number of reference reads, and its homozygous-non-reference term\n\
         `(ε/3)^(n−k)` then grows without limit as ε falls, railing the fit to the ladder's\n\
         floor. Each row below fits the same world twice, changing only which mean the cell is\n\
         scored at.\n"
    );
    let _ = writeln!(
        report,
        "{:<26} {:<22} {:>9} {:>9} {:>10} {:>11}",
        "world", "ladder", "n̄<k mass", "ε̄ rungs", "π_het", "π_hom_alt"
    );
    for world in binning_worlds()
        .into_iter()
        .filter(|w| w.mean_depth >= 20.0)
    {
        let truth = world.genotype_freqs();
        let true_mean_eps: f64 = world
            .eps
            .iter()
            .zip(&world.weights)
            .map(|(e, w)| e * w)
            .sum();
        for ladder in [
            DepthLadder::geometric(8, 16, 124),
            DepthLadder::geometric(8, 16, 300),
        ] {
            let space = CellSpace::build(&world, Coarsening::attributed(4, 0), &ladder);
            let depths = space.score_depths(DepthScoring::PerBin);
            let weights = space.mass.clone();
            let below: f64 = space
                .cells
                .iter()
                .zip(&space.mass)
                .zip(&depths)
                .filter(|((cell, _), d)| **d + 1e-12 < f64::from(cell_alt(cell)))
                .map(|((_, &m), _)| m)
                .sum();
            let fitted = fit(
                &space,
                &depths,
                Rule::Marginal,
                &world,
                &weights,
                &world.eps,
            );
            let fitted_mean_eps: f64 = fitted
                .eps
                .iter()
                .zip(&world.weights)
                .map(|(e, w)| e * w)
                .sum();
            let _ = writeln!(
                report,
                "{:<26} {:<22} {:>8.4}% {:>9.3} {:>9.2}% {:>10.2}%",
                world.name,
                ladder.label,
                100.0 * below,
                rungs(fitted_mean_eps / true_mean_eps),
                100.0 * (fitted.freqs[1] / truth[1] - 1.0),
                100.0 * (fitted.freqs[2] / truth[2] - 1.0),
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Experiment 2 — what assuming ½ costs
// ---------------------------------------------------------------------------

/// **What a fit that assumes `½ + ε/3` returns when the truth is not a half.**
///
/// Spec `parameter_prepass_generic.md` §8 leans toward replacing the `½` with a fitted
/// per-read-group constant, on the grounds that reads carrying the alternative allele map
/// slightly less often than reference-carrying ones, so a true heterozygote sits nearer
/// 0.47–0.49. §11.3 makes the decision turn on whether a fitted value departs from a half by
/// more than its standard error on real data — which needs a pipeline that does not exist.
/// The prior question does not: **if the misspecification costs nothing, the parameter is not
/// worth building whatever real data would say about it.**
///
/// Generate at balance `b`, fit with `b = ½`, and read the exact bias. Depth binning is
/// deliberately not mixed in: every world here runs on the exact ladder.
fn question_het_allele_balance(report: &mut String) {
    let _ = writeln!(report, "\n## Assuming a heterozygote is a half\n");
    let _ = writeln!(
        report,
        "One library, so nothing about the multi-library key is in these numbers. `b` is the\n\
         chance a read at a heterozygote came off the alternative copy before any misread; the\n\
         fit assumes ½ at every row. **`b = 0.50` is the control and must read exactly zero.**\n\
         `ε` in rungs of the quarter-Phred ladder, the two frequencies in relative error,\n\
         `spread` the movement across starts at 3× and ⅓ of the true rate.\n"
    );
    for &(pi_het, pi_hom_alt, rates_name) in &[
        (
            1e-2,
            6e-3,
            "10 het/kb, 6 hom-alt/kb — the other worlds' rates",
        ),
        (
            1e-3,
            6e-3,
            "1 het/kb, 6 hom-alt/kb — tomato's measured rate",
        ),
    ] {
        let _ = writeln!(report, "\n**{rates_name}**\n");
        let _ = writeln!(
            report,
            "{:<8} {:>7} {:>10} {:>10} {:>10} {:>11} {:>8}",
            "depth", "b", "ε rungs", "ε rel", "π_het", "π_hom_alt", "spread"
        );
        for &depth in &[3.0f64, 6.0, 10.0, 20.0, 60.0] {
            for &b in &[0.50f64, 0.49, 0.48, 0.47, 0.46, 0.45, 0.44] {
                let world = World {
                    name: format!("balance depth={depth:.0} b={b}"),
                    eps: vec![1e-3],
                    weights: vec![1.0],
                    mean_depth: depth,
                    pi_het,
                    pi_hom_alt,
                    het_balance: b,
                    site_noise: None,
                };
                let truth = world.genotype_freqs();
                let ladder = DepthLadder::exact(world.max_depth());
                let space = CellSpace::build(&world, Coarsening::pooled(), &ladder);
                let depths = space.score_depths(DepthScoring::PerCell);
                let weights = space.mass.clone();
                let mut fits = Vec::new();
                for scale in [1.0, 3.0, 1.0 / 3.0] {
                    let start: Vec<f64> = world.eps.iter().map(|e| e * scale).collect();
                    fits.push(fit(
                        &space,
                        &depths,
                        Rule::Marginal,
                        &world,
                        &weights,
                        &start,
                    ));
                }
                let fitted = &fits[0];
                let spread = fits
                    .iter()
                    .map(|f| rungs(f.eps[0] / world.eps[0]))
                    .fold(f64::NEG_INFINITY, f64::max)
                    - fits
                        .iter()
                        .map(|f| rungs(f.eps[0] / world.eps[0]))
                        .fold(f64::INFINITY, f64::min);
                let _ = writeln!(
                    report,
                    "{:<8.0} {:>7.2} {:>10.3} {:>9.2}% {:>9.2}% {:>10.2}% {:>8.3}",
                    depth,
                    b,
                    rungs(fitted.eps[0] / world.eps[0]),
                    100.0 * (fitted.eps[0] / world.eps[0] - 1.0),
                    100.0 * (fitted.freqs[1] / truth[1] - 1.0),
                    100.0 * (fitted.freqs[2] / truth[2] - 1.0),
                    spread,
                );
            }
            let _ = writeln!(report);
        }
    }
}

// ---------------------------------------------------------------------------
// Monte Carlo: the sample-size ladder
// ---------------------------------------------------------------------------

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
    fn binomial(&mut self, n: u32, p: f64) -> u32 {
        (0..n).filter(|_| self.unit() < p).count() as u32
    }
    fn poisson(&mut self, mean: f64) -> u32 {
        let limit = (-mean).exp();
        let mut k = 0u32;
        let mut prod = self.unit();
        while prod > limit {
            k += 1;
            prod *= self.unit();
            if k > 100_000 {
                break;
            }
        }
        k
    }
}

/// Draw `sites` sites and tally them under `coarsening`, returning the cells, their counts,
/// and the sum of the exact depths that landed in each — which is what the accumulator's
/// `depth_sums` holds and what a cell's mean depth is computed from.
fn simulate(
    world: &World,
    coarsening: Coarsening,
    ladder: &DepthLadder,
    sites: u64,
    seed: u64,
) -> (Vec<Cell>, Vec<f64>, Vec<f64>) {
    let mut rng = SplitMix64::new(seed);
    let freqs = world.genotype_freqs();
    let libraries = world.libraries();
    let mut table: HashMap<Cell, (f64, f64)> = HashMap::new();

    for _ in 0..sites {
        let u = rng.unit();
        let j = if u < freqs[0] {
            0
        } else if u < freqs[0] + freqs[1] {
            1
        } else {
            2
        };
        let depth = rng.poisson(world.mean_depth);
        if depth == 0 {
            continue;
        }
        let mut remaining = depth;
        let mut remaining_weight = 1.0;
        let mut per_lib = Vec::with_capacity(libraries);
        for g in 0..libraries {
            let n_g = if g + 1 == libraries {
                remaining
            } else {
                let p = (world.weights[g] / remaining_weight).clamp(0.0, 1.0);
                let drawn = rng.binomial(remaining, p);
                remaining -= drawn;
                remaining_weight -= world.weights[g];
                drawn
            };
            let k_g = rng.binomial(n_g, p_alt(j, world.eps[g], world.het_balance));
            per_lib.push((n_g, k_g));
        }
        let entry = table
            .entry(coarsening.key_of(&per_lib, ladder))
            .or_insert((0.0, 0.0));
        entry.0 += 1.0;
        entry.1 += f64::from(depth);
    }

    let mut cells = Vec::with_capacity(table.len());
    let mut counts = Vec::with_capacity(table.len());
    let mut depth_sums = Vec::with_capacity(table.len());
    let mut entries: Vec<_> = table.into_iter().collect();
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    for (cell, (count, depth_sum)) in entries {
        cells.push(cell);
        counts.push(count);
        depth_sums.push(depth_sum);
    }
    (cells, counts, depth_sums)
}

/// The sample-size ladder: does the estimate approach the truth as the data grow, and
/// how precise is it?
///
/// Heterozygosity is the headline because every key identifies it — the individual
/// error rates are not identified by the pooled key at all (see the `spread` column
/// above), so an error bar on those would be an error bar on where a search stopped.
/// `ε̄` is the share-weighted mean error rate, which every key does identify.
///
/// `sites× ` is the extra genome a key needs to match the exact per-library oracle's
/// precision: the ratio of variances. 1.00× is free.
fn monte_carlo_ladder(world: &World) {
    println!("\n## Sample-size ladder — {}\n", world.name);
    println!(
        "{:<38} {:>10} {:>11} {:>10} {:>11} {:>10} {:>8}",
        "candidate", "sites", "π_het bias", "π_het sd", "ε̄ bias", "ε̄ sd", "sites×"
    );

    let ladder: [u64; 4] = [10_000, 100_000, 1_000_000, 10_000_000];
    let repeats = 10;
    let true_het = world.genotype_freqs()[1];
    let true_mean_eps: f64 = world
        .eps
        .iter()
        .zip(&world.weights)
        .map(|(e, w)| e * w)
        .sum();
    let mut oracle_sd: HashMap<u64, f64> = HashMap::new();

    for candidate in candidates() {
        for &sites in &ladder {
            let mut het_errors = Vec::with_capacity(repeats);
            let mut eps_errors = Vec::with_capacity(repeats);
            for r in 0..repeats {
                let seed = 0x5EED_0000_u64
                    .wrapping_add(r as u64)
                    .wrapping_mul(0x9E37_79B9)
                    ^ sites;
                let ladder = DepthLadder::exact(world.max_depth());
                let (cells, counts, depth_sums) =
                    simulate(world, candidate.coarsening, &ladder, sites, seed);
                let total: f64 = counts.iter().sum();
                let space = CellSpace {
                    index: HashMap::new(),
                    cells,
                    mass: counts.iter().map(|c| c / total).collect(),
                    depth_sum: depth_sums.iter().map(|s| s / total).collect(),
                };
                let depths = space.score_depths(DepthScoring::PerCell);
                let weights = space.mass.clone();
                let fitted = fit(&space, &depths, candidate.rule, world, &weights, &world.eps);
                het_errors.push(fitted.freqs[1] / true_het - 1.0);
                let mean_eps: f64 = fitted
                    .eps
                    .iter()
                    .zip(&world.weights)
                    .map(|(e, w)| e * w)
                    .sum();
                eps_errors.push(rungs(mean_eps / true_mean_eps));
            }
            let summarise = |v: &[f64]| {
                let mean = v.iter().sum::<f64>() / v.len() as f64;
                let sd = (v.iter().map(|e| (e - mean).powi(2)).sum::<f64>() / (v.len() - 1) as f64)
                    .sqrt();
                (mean, sd)
            };
            let (het_bias, het_sd) = summarise(&het_errors);
            let (eps_bias, eps_sd) = summarise(&eps_errors);
            if candidate.coarsening.exact {
                oracle_sd.insert(sites, het_sd);
            }
            let cost = oracle_sd
                .get(&sites)
                .map(|o| format!("{:.2}x", (het_sd / o).powi(2)))
                .unwrap_or_else(|| "-".to_string());
            println!(
                "{:<38} {:>10} {:>10.2}% {:>9.2}% {:>11.3} {:>10.3} {:>8}",
                candidate.label,
                sites,
                100.0 * het_bias,
                100.0 * het_sd,
                eps_bias,
                eps_sd,
                cost
            );
        }
        println!();
    }
}

// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let run_monte_carlo = args.iter().any(|a| a == "--monte-carlo");
    // Each section is minutes rather than seconds, so they can be run one at a time.
    let wanted = |section: &str| {
        let selected: Vec<&String> = args.iter().filter(|a| a.starts_with("--only=")).collect();
        selected.is_empty() || selected.iter().any(|a| a.ends_with(section))
    };

    if wanted("binning") {
        let mut report = String::new();
        for probe_name in ["ratio=4 depth=20 skew90", "ratio=4 depth=60 skew90"] {
            let probe = binning_worlds()
                .into_iter()
                .find(|w| w.name == probe_name)
                .expect("the probe world is in the list");
            binning_checks(&mut report, &probe);
        }
        question_depth_binning(&mut report);
        question_per_bin_mean_depth(&mut report);
        print!("{report}");
    }
    if wanted("balance") {
        let mut report = String::new();
        question_het_allele_balance(&mut report);
        print!("{report}");
    }
    // Its own section, and printed the moment it is done: it profiles a rate at 25 points
    // and refits everything at each, so it is the slowest question here by an order of
    // magnitude. Buffered behind another section, a timeout would take both.
    if wanted("noise") {
        let mut report = String::new();
        question_the_second_class_of_site(&mut report);
        print!("{report}");
    }
    // Print-only, and reachable only by naming it: `--only=oracle`. It computes nothing
    // this program measures and changes no output any other section produces. It exists
    // because `ln_component_attributed` below is the oracle the caller's own multi-library
    // scoring rule is checked against (`impl_plan/parameter_prepass_generic.md` D2), and a
    // unit test cannot call into an example binary — so the numbers are dumped here and
    // pasted there, with this section as their provenance.
    if args.iter().any(|a| a == "--only=oracle") {
        dump_attributed_oracle();
        return;
    }

    if !wanted("key") {
        return;
    }

    println!("# Multi-library cell key — is the estimate unbiased, consistent, precise?");
    println!("#");
    println!("# Bias is computed exactly: each cell is weighted by its probability under");
    println!("# the truth, so the fit below is what an infinite genome would return. Units");
    println!("# are rungs of the design's error-rate ladder (a quarter-Phred, ~5.9% each)");
    println!("# for the error rates, and relative error for the two genotype frequencies.");
    println!("# A bias of exactly 0 means the estimator is consistent for that parameter.");

    let all = worlds();

    println!("\n## Checks that need no simulation\n");
    let probe = all
        .iter()
        .find(|w| w.name == "ratio=4 depth=20 split=skew90")
        .expect("the probe world is in the list");
    println!("(on {})\n", probe.name);
    print!("{}", algebraic_checks(probe));

    println!("\n## Asymptotic bias — what the estimate converges to\n");
    println!("Columns. `ε̄` is the share-weighted mean error rate `Σ_g w_g·ε_g`; `ε first`");
    println!("and `ε last` are the individual libraries' rates. `spread` is how far the");
    println!("answer for `ε last` moves when the search is started at 3× and at 1/3 of the");
    println!("truth instead of at it — a large spread means that rate is not identified by");
    println!("this key at all, and the zero beside it would be an artefact of starting on");
    println!("the answer. `clamp` is the share of sites where `n̂_g < k_g`.\n");

    for world in &all {
        let true_freqs = world.genotype_freqs();
        let true_mean_eps: f64 = world
            .eps
            .iter()
            .zip(&world.weights)
            .map(|(e, w)| e * w)
            .sum();
        println!(
            "### {}   ε = {:?}, shares = {:?}",
            world.name, world.eps, world.weights
        );
        println!(
            "{:<38} {:>8} {:>8} {:>8} {:>8} {:>9} {:>10} {:>8} {:>7}",
            "candidate", "ε̄", "ε first", "ε last", "spread", "π_het", "π_hom_alt", "clamp", "cells"
        );
        for candidate in candidates() {
            // The oracle enumerates every per-library breakdown, so its cell space grows
            // as the fourth power of the depth cap. It is run only where that is cheap —
            // its job here is to confirm the harness returns zero bias when the key
            // throws nothing away, and one depth settles that.
            if candidate.coarsening.exact && world.mean_depth > 5.0 {
                println!(
                    "{:<38} {:>8} {:>8} {:>8} {:>8} {:>9} {:>10} {:>8} {:>7}",
                    candidate.label, "-", "-", "-", "-", "-", "-", "-", "skipped"
                );
                continue;
            }
            let ladder = DepthLadder::exact(world.max_depth());
            let space = CellSpace::build(world, candidate.coarsening, &ladder);
            let depths = space.score_depths(DepthScoring::PerCell);
            let weights = space.mass.clone();
            let last = world.libraries() - 1;

            // Three starts, so a flat ridge in the likelihood cannot masquerade as an
            // unbiased answer.
            let starts: [f64; 3] = [1.0, 3.0, 1.0 / 3.0];
            let mut fits = Vec::with_capacity(starts.len());
            for scale in starts {
                let start: Vec<f64> = world.eps.iter().map(|e| e * scale).collect();
                fits.push(fit(
                    &space,
                    &depths,
                    candidate.rule,
                    world,
                    &weights,
                    &start,
                ));
            }
            let fitted = &fits[0];
            let spread = fits
                .iter()
                .map(|f| rungs(f.eps[last] / world.eps[last]))
                .fold(f64::NEG_INFINITY, f64::max)
                - fits
                    .iter()
                    .map(|f| rungs(f.eps[last] / world.eps[last]))
                    .fold(f64::INFINITY, f64::min);

            let fitted_mean_eps: f64 = fitted
                .eps
                .iter()
                .zip(&world.weights)
                .map(|(e, w)| e * w)
                .sum();

            // How much truth mass sits on a cell the plug-in has to clamp.
            let mut clamped_mass = 0.0;
            if candidate.rule == Rule::PlugAverageShare {
                for ((cell, &m), &depth) in space.cells.iter().zip(&space.mass).zip(&depths) {
                    if let Cell::Attributed { alt, .. } = cell
                        && alt
                            .iter()
                            .enumerate()
                            .any(|(g, &k)| world.weights[g] * depth < f64::from(k) - 1e-12)
                    {
                        clamped_mass += m;
                    }
                }
            }

            println!(
                "{:<38} {:>8.3} {:>8.3} {:>8.3} {:>8.3} {:>8.2}% {:>9.2}% {:>7.3}% {:>7}",
                candidate.label,
                rungs(fitted_mean_eps / true_mean_eps),
                rungs(fitted.eps[0] / world.eps[0]),
                rungs(fitted.eps[last] / world.eps[last]),
                spread,
                100.0 * (fitted.freqs[1] / true_freqs[1] - 1.0),
                100.0 * (fitted.freqs[2] / true_freqs[2] - 1.0),
                100.0 * clamped_mass,
                space.cells.len(),
            );
        }
        println!();
    }

    let mut coupled = String::new();
    question_coupled_fit(&mut coupled);
    print!("{coupled}");

    if run_monte_carlo {
        // Tomato's regime, and the one where the analytic bias is largest — also the
        // only depth where the exact per-library oracle is cheap enough to stand as the
        // precision ceiling.
        for world in all
            .iter()
            .filter(|w| w.name == "ratio=4 depth=3 split=skew90")
        {
            monte_carlo_ladder(world);
        }
    } else {
        println!("(re-run with --monte-carlo for the sample-size ladder and the variances)");
    }
}

/// Dump `ln_component_attributed` on a handful of cells, as Rust literals, so the
/// caller's own multi-library scoring rule can be checked against this program's.
///
/// **Reachable only by `--only=oracle`, and it measures nothing.** The world is
/// `ratio=4 depth=6 split=skew90` — two libraries at `ε` = 0.001 and 0.004 with a 90/10
/// split of the reads, which is the regime the attributed key exists for.
fn dump_attributed_oracle() {
    let weights = [0.9, 0.1];
    let eps = [1e-3, 4e-3];
    let splits: [[u32; 2]; 5] = [[1, 0], [0, 1], [2, 1], [1, 3], [0, 4]];

    println!(
        "// Generated by: cargo run --release --example ng_multilib_key_harness -- --only=oracle"
    );
    println!("// World: two libraries, eps = {eps:?}, shares = {weights:?}, ploidy 2.");
    println!("// Each row: depth, alternative reads from each library, then ln L at 0, 1");
    println!("// and 2 alternative copies.");
    for depth in [3.0_f64, 6.0, 6.5, 20.0] {
        for split in splits {
            let alt: u32 = split.iter().sum();
            if f64::from(alt) > depth {
                continue;
            }
            let ln: Vec<String> = (0..=PLOIDY)
                .map(|j| {
                    format!(
                        "{:e}",
                        ln_component_attributed(
                            depth,
                            &split,
                            &weights,
                            &eps,
                            MODEL_HET_BALANCE,
                            j
                        )
                    )
                })
                .collect();
            println!(
                "({depth:.1}, [{}, {}], [{}]),",
                split[0],
                split[1],
                ln.join(", ")
            );
        }
    }
}
