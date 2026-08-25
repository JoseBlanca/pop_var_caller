//! **Does the joint fit recover what it was given, and which way of describing the
//! allele frequency does it best?**
//!
//! The per-site route weights each sample's genotype by *that locus's* allele frequency
//! and integrates the frequency away
//! (`spec/parameter_prepass_joint_fit.md` §2.1). How that frequency is described is the
//! one open choice in the estimator, and this program measures the three candidates
//! against a truth it drew itself.
//!
//! ## The central measurement needs no simulation
//!
//! The method is `ng_multilib_key_harness.rs`'s, one level up. An estimator maximises a
//! sum over **the patterns a locus can show** — every sample's depth and alternative-read
//! count together. Replace the observed pattern counts with each pattern's *exact
//! probability under the truth* and you get the objective the estimator climbs with an
//! infinite genome:
//!
//! ```text
//! Q(θ)  =  Σ over patterns  P_true(pattern) · ln L_model(pattern ; θ)
//! θ*    =  argmax Q(θ)          the value the estimate converges to
//! bias  =  θ* − θ_true          a fixed number, not a sampling error
//! ```
//!
//! `θ*` is the pseudo-true value a misspecified fit converges to (White 1982), so
//! `bias = 0` and `bias ≠ 0` are decided exactly. **`bias = 0` is the statement that the
//! estimator is consistent.** The pattern space grows as `(depths × counts)^samples`, so
//! the exact arm runs at two or three samples; the drawn arm runs at fifty and uses the
//! identical objective with empirical weights in place of exact ones.
//!
//! ## The three candidates
//!
//! - **`density-4`** — the adopted one (spec §2.1.2): a mass on *the population carries
//!   only the reference base*, a mass on *only a non-reference base*, and a Beta over what
//!   segregates. Four numbers, and the grid it is evaluated on is quadrature rather than
//!   freedom.
//! - **`grid-free`** — the same frequency form with **one free weight per grid cell**,
//!   `2N + 1` of them. This is what an earlier draft of the spec described. Recovering a
//!   frequency density from genotypes means undoing a binomial blur, and an unregularised
//!   deconvolution onto a grid this fine is expected to collapse onto spikes; whether that
//!   costs anything a caller would feel is what this arm measures.
//! - **`count-sfs`** — the *count* form, one free weight per allele count in the panel,
//!   with the convolution that keeps the samples' genotypes adding up to that count
//!   (ANGSD's `realSFS`; Nielsen et al. 2012). **It cannot carry a per-sample inbreeding
//!   coefficient** — the cancellation that makes it work needs each genotype to be a pair
//!   of independent draws (spec §2.1.1) — so it is fitted with inbreeding held at zero and
//!   is only comparable where the truth has none. It is here because it is the
//!   best-conditioned of the three when it applies, and that is worth knowing.
//!
//! ## Running it
//!
//! ```text
//! ng_joint_fit_harness [exact|drawn|both|sweep|structured] [samples] [depth] [loci]
//! ```
//!
//! `sweep` refits the adopted candidate at 1, 2, 5, 10, 25 and 50 samples and reports
//! where each parameter stops being estimable (spec §12.5); `structured` fits a truth
//! whose frequency density the fitted shape cannot reach (spec §12.9).
//!
//! Diploid throughout: the Hardy–Weinberg-with-inbreeding genotype prior is diploid-only
//! and the spec defers its restatement above `P = 2` (§10).

use rayon::prelude::*;
use std::time::Instant;

// ---------------------------------------------------------------------
// Small numerics
// ---------------------------------------------------------------------

fn ln_gamma(x: f64) -> f64 {
    // Lanczos, g = 7, n = 9 — the coefficients every implementation uses.
    //
    // Written to the digit count the published set is published with, so a reader
    // checking them against a reference sees the same characters. 1.98's clippy asks
    // for the last digit of the first one to go; it parses to the same `f64` either
    // way, and matching the reference is worth more here than matching the lint.
    #[allow(clippy::excessive_precision)]
    const C: [f64; 9] = [
        0.999_999_999_999_809_93,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection, so the fit may probe shape parameters near zero without care.
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().abs().ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + 7.5;
        for (i, c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

fn ln_binomial_coefficient(n: u32, k: u32) -> f64 {
    ln_gamma(f64::from(n) + 1.0) - ln_gamma(f64::from(k) + 1.0) - ln_gamma(f64::from(n - k) + 1.0)
}

/// `ln` of the binomial probability of `k` successes in `n` trials at `p`.
fn ln_binomial(k: u32, n: u32, p: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mut out = ln_binomial_coefficient(n, k);
    out += if k > 0 { f64::from(k) * p.ln() } else { 0.0 };
    out += if n > k {
        f64::from(n - k) * (1.0 - p).ln()
    } else {
        0.0
    };
    out
}

fn ln_sum_exp(values: &[f64]) -> f64 {
    let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    max + values.iter().map(|v| (v - max).exp()).sum::<f64>().ln()
}

/// The regularised incomplete beta `I_x(a, b)` — the Beta's own cumulative, which is what
/// turns a grid of frequency cells into exact cell masses.
///
/// Continued fraction, the standard Lentz evaluation. Needed because the Beta's shape is
/// fitted, so no fixed quadrature rule's nodes can be precomputed for it.
fn regularised_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    // The fraction converges quickly only on one side of the mean; reflect otherwise.
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularised_incomplete_beta(1.0 - x, b, a);
    }
    let tiny = 1e-300;
    let mut f = 1.0_f64;
    let mut c = 1.0_f64;
    let mut d = 0.0_f64;
    for i in 0..=300 {
        let m = i / 2;
        let numerator = if i == 0 {
            1.0
        } else if i % 2 == 0 {
            let m = f64::from(m);
            (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        } else {
            let m = f64::from(m);
            -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        };
        d = 1.0 + numerator * d;
        if d.abs() < tiny {
            d = tiny;
        }
        d = 1.0 / d;
        c = 1.0 + numerator / c;
        if c.abs() < tiny {
            c = tiny;
        }
        let step = c * d;
        f *= step;
        if (1.0 - step).abs() < 1e-12 {
            break;
        }
    }
    front * (f - 1.0) / a
}

/// A deterministic generator, so a run is reproducible without a seed file.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        // splitmix64.
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1_u64 << 53) as f64
    }
    fn poisson(&mut self, mean: f64) -> u32 {
        // Knuth; the means here are single digits.
        let limit = (-mean).exp();
        let mut k = 0;
        let mut p = 1.0;
        loop {
            p *= self.uniform();
            if p <= limit {
                return k;
            }
            k += 1;
            if k > 1_000 {
                return k;
            }
        }
    }
    fn binomial(&mut self, n: u32, p: f64) -> u32 {
        (0..n).filter(|_| self.uniform() < p).count() as u32
    }
    fn weighted_index(&mut self, weights: &[f64]) -> usize {
        let mut u = self.uniform() * weights.iter().sum::<f64>();
        for (i, w) in weights.iter().enumerate() {
            u -= w;
            if u <= 0.0 {
                return i;
            }
        }
        weights.len() - 1
    }
    /// A Beta draw by the ratio of two Gammas, each by Marsaglia–Tsang with the
    /// boost that makes shapes below one legal.
    fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        if x + y == 0.0 { 0.5 } else { x / (x + y) }
    }
    fn gamma(&mut self, shape: f64) -> f64 {
        if shape < 1.0 {
            let u = self.uniform().max(1e-300);
            return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let x = self.normal();
            let v = (1.0 + c * x).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.uniform().max(1e-300);
            if u.ln() < 0.5 * x * x + d - d * v + d * v.ln() {
                return d * v;
            }
        }
    }
    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------------------------------------------------------------
// The model's vocabulary
// ---------------------------------------------------------------------

const PLOIDY: u32 = 2;

/// The chance one read shows something other than the reference base, for an individual
/// carrying `j` alternative copies of `PLOIDY` — `spec/parameter_prepass.md` §3, including
/// the `ε/3` that an earlier version of that table got wrong.
fn p_alt(j: u32, eps: f64) -> f64 {
    let share = f64::from(j) / f64::from(PLOIDY);
    share * (1.0 - eps / 3.0) + (1.0 - share) * eps
}

/// The diploid genotype frequencies at population frequency `f` with inbreeding `inbreed`
/// — Hardy–Weinberg with the inbreeding correction, which is the factorisation the count
/// form cannot use (spec §2.1.1).
fn genotype_frequencies(f: f64, inbreed: f64) -> [f64; 3] {
    let hetero = 2.0 * f * (1.0 - f) * (1.0 - inbreed);
    let shared = inbreed * f * (1.0 - f);
    [(1.0 - f) * (1.0 - f) + shared, hetero, f * f + shared]
}

/// The frequency density: two masses and a Beta (spec §2.1.2).
#[derive(Clone, Copy, Debug)]
struct Density {
    p_invariant: f64,
    p_fixed_alt: f64,
    a: f64,
    b: f64,
    /// A **second bump** over the segregating frequencies: `(its share, a₂, b₂)`.
    ///
    /// Three more fitted numbers, not two — the share is one of them. Present so the cost
    /// of adding shape can be measured rather than argued: what it buys is a truth with
    /// two bumps in it, and what it risks is that two similar bumps are interchangeable,
    /// so the fit lands wherever it started.
    second: Option<(f64, f64, f64)>,
}

impl Density {
    fn p_segregating(&self) -> f64 {
        1.0 - self.p_invariant - self.p_fixed_alt
    }

    /// The population's expected heterozygosity, `∫ π(f) · 2f(1−f) df` (spec §5.3). For a
    /// Beta that integral is closed form: `2·a·b / ((a+b)(a+b+1))`.
    fn expected_heterozygosity(&self) -> f64 {
        let of = |a: f64, b: f64| 2.0 * a * b / ((a + b) * (a + b + 1.0));
        let shape = match self.second {
            None => of(self.a, self.b),
            Some((share, a2, b2)) => (1.0 - share) * of(self.a, self.b) + share * of(a2, b2),
        };
        self.p_segregating() * shape
    }
}

/// The two grids a run holds, and they are different kinds of thing.
///
/// - `quadrature` is **accuracy**: where the Beta's integral is evaluated. Its cell count
///   is an accuracy knob and adds no freedom to the fit (spec §2.1.2).
/// - `free` is **freedom**: the frequencies `grid-free` puts one weight on each of, at
///   `c / 2N` for `c = 0..2N` — a faithful reconstruction of the earlier draft's *"one
///   weight per possible allele count in the panel"*, endpoints included so an invariant
///   position has a component to live in.
struct Grids {
    quadrature: FrequencyGrid,
    free: Vec<f64>,
}

impl Grids {
    fn new(samples: usize, quadrature_cells: usize) -> Self {
        let chromosomes = 2 * samples;
        Self {
            quadrature: FrequencyGrid::new(quadrature_cells),
            free: (0..=chromosomes)
                .map(|c| c as f64 / chromosomes as f64)
                .collect(),
        }
    }
}

/// The frequency grid the segregating mass is evaluated on.
///
/// **Cell edges are fixed; only the masses move.** Arcsine spacing clusters cells at both
/// ends, which is where a Beta with shapes below one puts its mass and where a rare allele
/// lives. Each cell's mass and its conditional mean are exact, from the Beta's own
/// cumulative — so this is quadrature, not a discretisation of the parameter.
struct FrequencyGrid {
    edges: Vec<f64>,
}

impl FrequencyGrid {
    fn new(cells: usize) -> Self {
        let edges = (0..=cells)
            .map(|e| 0.5 * (1.0 - (std::f64::consts::PI * e as f64 / cells as f64).cos()))
            .collect();
        Self { edges }
    }

    fn cells(&self) -> usize {
        self.edges.len() - 1
    }

    /// Each cell's `(mass, representative frequency)` under `Beta(a, b)`.
    /// The cells under a density with one bump or two — the mixture is taken on the
    /// **cell masses**, which is exact, rather than on the representatives.
    fn under_density(&self, density: &Density) -> Vec<(f64, f64)> {
        let first = self.under_beta(density.a, density.b);
        let Some((share, a2, b2)) = density.second else {
            return first;
        };
        let second = self.under_beta(a2, b2);
        first
            .iter()
            .zip(&second)
            .map(|(&(m1, f1), &(m2, f2))| {
                let mass = (1.0 - share) * m1 + share * m2;
                // The cell's representative is the two components' conditional means
                // averaged by the mass each actually put in this cell.
                let representative = if mass > 1e-300 {
                    ((1.0 - share) * m1 * f1 + share * m2 * f2) / mass
                } else {
                    0.5 * (f1 + f2)
                };
                (mass, representative)
            })
            .collect()
    }

    fn under_beta(&self, a: f64, b: f64) -> Vec<(f64, f64)> {
        let mean = a / (a + b);
        let mut lower = 0.0;
        let mut lower_shifted = 0.0;
        let mut out = Vec::with_capacity(self.cells());
        for e in 1..self.edges.len() {
            let upper = regularised_incomplete_beta(self.edges[e], a, b);
            let upper_shifted = regularised_incomplete_beta(self.edges[e], a + 1.0, b);
            let mass = (upper - lower).max(0.0);
            // The cell's conditional mean, which is the representative a first-order
            // quadrature wants — never the midpoint, since the density inside a cell at
            // the ends of the range is anything but flat.
            let representative = if mass > 1e-300 {
                (mean * (upper_shifted - lower_shifted) / mass)
                    .clamp(self.edges[e - 1], self.edges[e])
            } else {
                0.5 * (self.edges[e - 1] + self.edges[e])
            };
            out.push((mass, representative));
            lower = upper;
            lower_shifted = upper_shifted;
        }
        let total: f64 = out.iter().map(|(m, _)| m).sum();
        if total > 0.0 {
            for cell in &mut out {
                cell.0 /= total;
            }
        }
        out
    }
}

// ---------------------------------------------------------------------
// The truth, and the data
// ---------------------------------------------------------------------

#[derive(Clone, Debug)]
struct Truth {
    samples: usize,
    depth_mean: f64,
    eps_clean: f64,
    eps_noisy: f64,
    w_noisy: f64,
    density: Density,
    inbreeding: Vec<f64>,
    /// When set, the segregating frequencies are drawn from **two** Beta components in
    /// equal parts rather than one — `[a1, b1, a2, b2]`.
    ///
    /// **This is a truth the fitted shape cannot reach**, which is the point: a panel that
    /// is landraces from several regions has a frequency density with more than one bump,
    /// and a single Beta has one. A fit graded only on data drawn from its own family
    /// reports that its arithmetic is right and nothing about its model (spec §12.1).
    /// Drawn-arm only — the exact arm's truth has to be the same object the model
    /// evaluates.
    bimodal: Option<[f64; 4]>,
}

impl Truth {
    /// The population's expected heterozygosity under whichever density this truth holds.
    fn expected_heterozygosity(&self) -> f64 {
        match self.bimodal {
            None => self.density.expected_heterozygosity(),
            Some([a1, b1, a2, b2]) => {
                let of = |a: f64, b: f64| 2.0 * a * b / ((a + b) * (a + b + 1.0));
                self.density.p_segregating() * 0.5 * (of(a1, b1) + of(a2, b2))
            }
        }
    }

    fn tomato_like(samples: usize, depth_mean: f64, inbreeding: f64) -> Self {
        Self {
            samples,
            depth_mean,
            // The rates the histogram route fits on HG002, which is the only place they
            // have been measured against a truth set
            // (`research/noise_model_overdispersion_2026-08-10.md`).
            eps_clean: 1.895e-3,
            eps_noisy: 5.29e-2,
            w_noisy: 0.0088,
            density: Density {
                p_invariant: 0.99,
                p_fixed_alt: 0.002,
                // `a < 1` is the rare-allele pile-up a neutral population has.
                a: 0.3,
                b: 1.2,
                second: None,
            },
            inbreeding: vec![inbreeding; samples],
            bimodal: None,
        }
    }

    /// The same truth with a two-bump frequency density — one component around 0.2 and
    /// one around 0.8, in equal parts, which is the shape two diverged subpopulations
    /// leave behind.
    fn structured(samples: usize, depth_mean: f64, inbreeding: f64) -> Self {
        Self {
            bimodal: Some([6.0, 24.0, 24.0, 6.0]),
            ..Self::tomato_like(samples, depth_mean, inbreeding)
        }
    }
}

/// One locus's evidence: `(depth, alternative reads)` for every sample, in sample order.
type Pattern = Vec<(u32, u32)>;

/// Patterns with the weight each carries — an exact probability in the exact arm, a count
/// in the drawn one. **One objective reads both**, which is what makes the drawn arm a
/// check on the exact one rather than a separate program.
struct Weighted {
    /// Kept only where a candidate has to score every sample one at a time — the count
    /// form does, and so does a fit with an inbreeding coefficient per sample. **Dropped
    /// above a few dozen samples**, where a pair per sample per locus is gigabytes and the
    /// fast path never reads it.
    patterns: Option<Vec<Pattern>>,
    shapes: Vec<PatternShape>,
    weights: Vec<f64>,
    /// The deepest any sample went, which sizes the by-depth tables.
    deepest: usize,
    /// How many of the drawn loci segregate **in the population** — a property of the
    /// truth, and the same at every panel size.
    segregating: usize,
    /// How many segregate **in the panel**: at least one of the drawn chromosomes carries
    /// the allele and at least one does not.
    ///
    /// **These two counts come apart as the panel grows, and it is the second that
    /// contamination needs** — a marker is only usable where the cohort's own alleles
    /// differ there. A panel of two thousand samples sees alleles down to a frequency of
    /// about one in four thousand; fifty samples sees nothing below one in a hundred. So a
    /// locus budget set by contamination should fall as the cohort grows, where one set by
    /// the noisy-locus share should not.
    panel_segregating: usize,
}

/// Every sample's per-genotype read likelihood at one locus, `ln` scale, for each site
/// class. Precomputed per pattern because every candidate reads the same thing.
struct PatternLikelihoods {
    /// `[class][sample][genotype]`
    ln_read: Vec<Vec<[f64; 3]>>,
}

impl PatternLikelihoods {
    fn new(pattern: &Pattern, classes: &[f64]) -> Self {
        let ln_read = classes
            .iter()
            .map(|&eps| {
                pattern
                    .iter()
                    .map(|&(depth, alt)| {
                        [
                            ln_binomial(alt, depth, p_alt(0, eps)),
                            ln_binomial(alt, depth, p_alt(1, eps)),
                            ln_binomial(alt, depth, p_alt(2, eps)),
                        ]
                    })
                    .collect()
            })
            .collect();
        Self { ln_read }
    }
}

/// A locus's samples split into the ones that said nothing and the ones that did.
///
/// **This is what makes a large panel cost what a small one does.** A sample with no
/// alternative read contributes a factor that depends only on its **depth** and on its
/// inbreeding coefficient — not on which sample it is — so where every sample shares one
/// coefficient, the product over a thousand quiet samples is a table of a dozen numbers
/// raised to counts. Only the samples that showed an alternative read have to be scored
/// one at a time, and at three reads a site and an error rate of 0.002 that is about six
/// samples in a thousand plus whoever really carries the variant.
///
/// Computed once per locus and independent of every parameter, so it survives the search.
struct PatternShape {
    /// How many silent samples sat at each depth. Index is the depth.
    quiet_at_depth: Vec<u32>,
    /// The samples that showed at least one alternative read, as `(depth, alternative
    /// reads)`. **Which sample each came from is not kept**: where one inbreeding
    /// coefficient is shared across the panel, a loud observation's identity changes
    /// nothing about the likelihood.
    loud: Vec<(u32, u32)>,
}

impl PatternShape {
    fn new(pattern: &Pattern) -> Self {
        let deepest = pattern.iter().map(|&(d, _)| d).max().unwrap_or(0) as usize;
        let mut quiet_at_depth = vec![0_u32; deepest + 1];
        let mut loud = Vec::new();
        for &(depth, alt) in pattern {
            if alt == 0 {
                quiet_at_depth[depth as usize] += 1;
            } else {
                loud.push((depth, alt));
            }
        }
        Self {
            quiet_at_depth,
            loud,
        }
    }

    /// The loud samples' per-genotype read likelihoods, `[class][loud index][genotype]` —
    /// the only likelihoods the fast path evaluates.
    fn loud_likelihoods(&self, classes: &[f64]) -> PatternLikelihoods {
        PatternLikelihoods {
            ln_read: classes
                .iter()
                .map(|&eps| {
                    self.loud
                        .iter()
                        .map(|&(depth, alt)| {
                            [
                                ln_binomial(alt, depth, p_alt(0, eps)),
                                ln_binomial(alt, depth, p_alt(1, eps)),
                                ln_binomial(alt, depth, p_alt(2, eps)),
                            ]
                        })
                        .collect()
                })
                .collect(),
        }
    }
}

/// Everything a silent sample can contribute, tabulated by depth so that no locus has to
/// evaluate it per sample.
///
/// Built once per objective evaluation — two site classes by a couple of dozen frequency
/// cells by a dozen depths is a few hundred entries, against one evaluation per sample per
/// locus without it. **Only legal where every sample shares one inbreeding coefficient**,
/// which is what `shared` in the fit's signature means; with a coefficient per sample the
/// quiet samples stop being interchangeable and the slow path is the only one.
struct QuietTables {
    /// `[class][cell][depth]` — a silent sample's contribution to the segregating term.
    segregating: Vec<Vec<Vec<f64>>>,
    /// `[class][depth]` — its contribution where the population carries only the reference
    /// base, and only a non-reference one.
    all_reference: Vec<Vec<f64>>,
    all_alternative: Vec<Vec<f64>>,
}

impl QuietTables {
    fn new(classes: &[f64], cells: &[(f64, f64)], inbreeding: f64, deepest: usize) -> Self {
        let depths = || 0..=deepest as u32;
        Self {
            segregating: classes
                .iter()
                .map(|&eps| {
                    cells
                        .iter()
                        .map(|&(_, f)| {
                            let priors = genotype_frequencies(f, inbreeding);
                            depths()
                                .map(|depth| {
                                    ln_sum_exp(&[
                                        priors[0].max(1e-300).ln()
                                            + ln_binomial(0, depth, p_alt(0, eps)),
                                        priors[1].max(1e-300).ln()
                                            + ln_binomial(0, depth, p_alt(1, eps)),
                                        priors[2].max(1e-300).ln()
                                            + ln_binomial(0, depth, p_alt(2, eps)),
                                    ])
                                })
                                .collect()
                        })
                        .collect()
                })
                .collect(),
            all_reference: classes
                .iter()
                .map(|&eps| depths().map(|d| ln_binomial(0, d, p_alt(0, eps))).collect())
                .collect(),
            all_alternative: classes
                .iter()
                .map(|&eps| depths().map(|d| ln_binomial(0, d, p_alt(2, eps))).collect())
                .collect(),
        }
    }
}

/// [`ln_mixture`], with the silent samples read off [`QuietTables`] instead of one at a
/// time. **Returns the same number** — which its own test pins, because a fast path that
/// quietly disagrees with the slow one would look like a modelling result.
fn ln_mixture_shared(
    shape: &PatternShape,
    likelihoods: &PatternLikelihoods,
    class_weights: &[f64],
    density: &Density,
    cells: &[(f64, f64)],
    ln_priors: &[[f64; 3]],
    quiet: &QuietTables,
) -> f64 {
    let mut per_class = Vec::with_capacity(class_weights.len());
    for (class, &class_weight) in class_weights.iter().enumerate() {
        if class_weight <= 0.0 {
            per_class.push(f64::NEG_INFINITY);
            continue;
        }
        let ln_read = &likelihoods.ln_read[class];
        let quiet_sum = |table: &[f64]| -> f64 {
            shape
                .quiet_at_depth
                .iter()
                .enumerate()
                .filter(|&(_, &count)| count > 0)
                .map(|(depth, &count)| f64::from(count) * table[depth])
                .sum()
        };
        let mut terms = Vec::with_capacity(cells.len() + 2);
        if density.p_invariant > 0.0 {
            let mut total = density.p_invariant.ln() + quiet_sum(&quiet.all_reference[class]);
            total += ln_read.iter().map(|g| g[0]).sum::<f64>();
            terms.push(total);
        }
        if density.p_fixed_alt > 0.0 {
            let mut total = density.p_fixed_alt.ln() + quiet_sum(&quiet.all_alternative[class]);
            total += ln_read.iter().map(|g| g[2]).sum::<f64>();
            terms.push(total);
        }
        let segregating = density.p_segregating();
        if segregating > 0.0 {
            for (cell, &(mass, _)) in cells.iter().enumerate() {
                if mass <= 0.0 {
                    continue;
                }
                let mut total =
                    segregating.ln() + mass.ln() + quiet_sum(&quiet.segregating[class][cell]);
                let ln_prior = &ln_priors[cell];
                for ln_g in ln_read.iter() {
                    total += ln_sum_exp(&[
                        ln_prior[0] + ln_g[0],
                        ln_prior[1] + ln_g[1],
                        ln_prior[2] + ln_g[2],
                    ]);
                }
                terms.push(total);
            }
        }
        per_class.push(class_weight.ln() + ln_sum_exp(&terms));
    }
    ln_sum_exp(&per_class)
}

// ---------------------------------------------------------------------
// Building the data
// ---------------------------------------------------------------------

/// Every pattern a locus can show at depths up to `depth_cap`, with its exact probability
/// under the truth. **This is the whole of the exact arm's input.**
fn exact_patterns(truth: &Truth, depth_cap: u32, grids: &Grids) -> Weighted {
    let per_sample: Vec<(u32, u32)> = (0..=depth_cap)
        .flat_map(|d| (0..=d).map(move |k| (d, k)))
        .collect();
    let mut patterns = Vec::new();
    let mut current = Vec::with_capacity(truth.samples);
    enumerate(&per_sample, truth.samples, &mut current, &mut patterns);

    // Depth is drawn independently of everything else, so it multiplies in.
    //
    // **Truncated at the cap and renormalised, deliberately.** Leaving the Poisson's tail
    // outside the enumerated patterns would make the exact arm's data a *truncated* sample
    // of the truth — and truncation removes the deepest sites, which are the most
    // informative about an error rate, so what came back would be a bias of the harness
    // rather than of the estimator. Renormalising makes the truth here exactly "a Poisson
    // conditioned on landing at or below the cap", which is a truth the model can hit.
    let ln_poisson: Vec<f64> = (0..=depth_cap)
        .map(|d| {
            f64::from(d) * truth.depth_mean.ln() - truth.depth_mean - ln_gamma(f64::from(d) + 1.0)
        })
        .collect();
    let retained: f64 = ln_poisson.iter().map(|l| l.exp()).sum();
    let ln_depth: Vec<f64> = ln_poisson.iter().map(|l| l - retained.ln()).collect();

    let classes = [truth.eps_clean, truth.eps_noisy];
    let class_weights = [1.0 - truth.w_noisy, truth.w_noisy];
    let cells = grids.quadrature.under_density(&truth.density);
    let ln_priors = genotype_priors(&cells, &truth.inbreeding);
    let weights = patterns
        .iter()
        .map(|pattern| {
            let ln_depth_term: f64 = pattern.iter().map(|&(d, _)| ln_depth[d as usize]).sum();
            let likelihoods = PatternLikelihoods::new(pattern, &classes);
            let ln = ln_mixture(
                &likelihoods,
                &class_weights,
                &truth.density,
                &cells,
                &ln_priors,
            );
            (ln_depth_term + ln).exp()
        })
        .collect::<Vec<f64>>();

    // **Prune, and report what was dropped.** The pattern space is dominated by
    // combinations no genome would ever show — every sample at its deepest and every read
    // disagreeing — and each carries a weight near zero while costing the same to score.
    // The retained mass is printed, so a cut that mattered would be visible rather than
    // silent.
    let total: f64 = weights.iter().sum();
    let floor = total * 1e-13;
    let (patterns, weights): (Vec<Pattern>, Vec<f64>) = patterns
        .into_iter()
        .zip(weights)
        .filter(|(_, w)| *w > floor)
        .unzip();
    let shapes: Vec<PatternShape> = patterns.iter().map(PatternShape::new).collect();
    let deepest = depth_cap as usize;
    Weighted {
        patterns: Some(patterns),
        shapes,
        weights,
        deepest,
        // Every pattern is enumerated here rather than drawn, so "how many segregate" is
        // not a property of this arm at all.
        segregating: 0,
        panel_segregating: 0,
    }
}

fn enumerate(
    per_sample: &[(u32, u32)],
    remaining: usize,
    current: &mut Pattern,
    out: &mut Vec<Pattern>,
) {
    if remaining == 0 {
        out.push(current.clone());
        return;
    }
    for &cell in per_sample {
        current.push(cell);
        enumerate(per_sample, remaining - 1, current, out);
        current.pop();
    }
}

/// `ln` of each sample's genotype prior at each grid cell — `[cell][sample][genotype]`.
///
/// **Hoisted out of the pattern loop deliberately.** It depends on the cell's frequency and
/// the sample's inbreeding and on nothing about the reads, so computing it inside the loop
/// over loci would repeat two million identical evaluations per objective call. This is
/// the whole difference between a fit that runs in seconds and one that does not.
fn genotype_priors(cells: &[(f64, f64)], inbreeding: &[f64]) -> Vec<Vec<[f64; 3]>> {
    cells
        .iter()
        .map(|&(_, f)| {
            inbreeding
                .iter()
                .map(|&inbreed| {
                    let p = genotype_frequencies(f, inbreed);
                    [
                        p[0].max(1e-300).ln(),
                        p[1].max(1e-300).ln(),
                        p[2].max(1e-300).ln(),
                    ]
                })
                .collect()
        })
        .collect()
}

/// The locus likelihood every frequency-form candidate evaluates: sum over site class,
/// then the three frequency terms — invariant, fixed-alternative, and the segregating
/// integral over the grid (spec §3.1).
fn ln_mixture(
    likelihoods: &PatternLikelihoods,
    class_weights: &[f64],
    density: &Density,
    cells: &[(f64, f64)],
    ln_priors: &[Vec<[f64; 3]>],
) -> f64 {
    let mut per_class = Vec::with_capacity(class_weights.len());
    for (class, &class_weight) in class_weights.iter().enumerate() {
        if class_weight <= 0.0 {
            per_class.push(f64::NEG_INFINITY);
            continue;
        }
        let ln_read = &likelihoods.ln_read[class];
        let mut terms = Vec::with_capacity(cells.len() + 2);
        if density.p_invariant > 0.0 {
            terms.push(density.p_invariant.ln() + ln_read.iter().map(|g| g[0]).sum::<f64>());
        }
        if density.p_fixed_alt > 0.0 {
            terms.push(density.p_fixed_alt.ln() + ln_read.iter().map(|g| g[2]).sum::<f64>());
        }
        let segregating = density.p_segregating();
        if segregating > 0.0 {
            for (cell, &(mass, _)) in cells.iter().enumerate() {
                if mass <= 0.0 {
                    continue;
                }
                let mut total = segregating.ln() + mass.ln();
                for (sample, ln_g) in ln_read.iter().enumerate() {
                    let ln_prior = &ln_priors[cell][sample];
                    total += ln_sum_exp(&[
                        ln_prior[0] + ln_g[0],
                        ln_prior[1] + ln_g[1],
                        ln_prior[2] + ln_g[2],
                    ]);
                }
                terms.push(total);
            }
        }
        per_class.push(class_weight.ln() + ln_sum_exp(&terms));
    }
    ln_sum_exp(&per_class)
}

/// Draw `loci` loci from the truth — the arm that runs at fifty samples, where the exact
/// pattern space does not fit in a universe.
fn drawn_patterns(truth: &Truth, loci: usize, seed: u64) -> Weighted {
    // Above a few dozen samples the full patterns are gigabytes and nothing reads them:
    // the count form and the per-sample-inbreeding path are the only callers, and both
    // are small-panel arms.
    let keep_patterns = truth.samples <= 64;
    let mut rng = Rng(seed);
    let mut patterns = Vec::with_capacity(if keep_patterns { loci } else { 0 });
    let mut shapes = Vec::with_capacity(loci);
    let mut deepest = 0_usize;
    let mut segregating = 0_usize;
    let mut panel_segregating = 0_usize;
    for _ in 0..loci {
        let eps = if rng.uniform() < truth.w_noisy {
            truth.eps_noisy
        } else {
            truth.eps_clean
        };
        let u = rng.uniform();
        let f = if u < truth.density.p_invariant {
            0.0
        } else if u < truth.density.p_invariant + truth.density.p_fixed_alt {
            1.0
        } else {
            match truth.bimodal {
                None => rng.beta(truth.density.a, truth.density.b),
                Some([a1, b1, a2, b2]) => {
                    if rng.uniform() < 0.5 {
                        rng.beta(a1, b1)
                    } else {
                        rng.beta(a2, b2)
                    }
                }
            }
        };
        if f > 0.0 && f < 1.0 {
            segregating += 1;
        }
        let mut alternative_copies = 0_u32;
        let pattern: Pattern = (0..truth.samples)
            .map(|sample| {
                let priors = genotype_frequencies(f, truth.inbreeding[sample]);
                let genotype = rng.weighted_index(&priors) as u32;
                alternative_copies += genotype;
                let depth = rng.poisson(truth.depth_mean);
                (depth, rng.binomial(depth, p_alt(genotype, eps)))
            })
            .collect();
        if alternative_copies > 0 && alternative_copies < PLOIDY * truth.samples as u32 {
            panel_segregating += 1;
        }
        deepest = deepest.max(pattern.iter().map(|&(d, _)| d as usize).max().unwrap_or(0));
        shapes.push(PatternShape::new(&pattern));
        if keep_patterns {
            patterns.push(pattern);
        }
    }
    let weights = vec![1.0 / loci as f64; loci];
    Weighted {
        patterns: keep_patterns.then_some(patterns),
        shapes,
        weights,
        deepest,
        segregating,
        panel_segregating,
    }
}

// ---------------------------------------------------------------------
// The candidates
// ---------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
enum Candidate {
    /// Two bumps over the segregating frequencies instead of one — three more fitted
    /// numbers, and the arm that measures what buying shape costs.
    Density7,
    /// Two masses and a Beta — the adopted form (spec §2.1.2).
    Density4,
    /// One free weight per grid cell, `2N + 1` cells.
    GridFree,
    /// One free weight per panel allele count, with the convolution — inbreeding held at
    /// zero, since the form cannot carry it (spec §2.1.1).
    CountSpectrum,
}

impl Candidate {
    fn name(self) -> &'static str {
        match self {
            Candidate::Density4 => "density-4",
            Candidate::Density7 => "density-7",
            Candidate::GridFree => "grid-free",
            Candidate::CountSpectrum => "count-sfs",
        }
    }
}

/// What a fit holds and moves.
#[derive(Clone)]
struct Fitted {
    eps_clean: f64,
    eps_noisy: f64,
    w_noisy: f64,
    density: Density,
    inbreeding: Vec<f64>,
    /// `grid-free`'s own masses over the grid cells; `count-sfs`'s over allele counts.
    weights: Vec<f64>,
}

impl Fitted {
    /// The truth's own values, so the objective can be evaluated there. **`Q(θ_true)` is
    /// the reference every fitted score is read against**: a fitted score *below* it is an
    /// optimiser that stopped early, which is a different failure from a biased estimator
    /// and looks identical without this number.
    fn at_truth(truth: &Truth, grids: &Grids) -> Self {
        Self {
            eps_clean: truth.eps_clean,
            eps_noisy: truth.eps_noisy,
            w_noisy: truth.w_noisy,
            density: truth.density,
            inbreeding: truth.inbreeding.clone(),
            weights: vec![1.0 / grids.free.len() as f64; grids.free.len()],
        }
    }

    /// A starting point. **`variant` spans the separation between the clean and the noisy
    /// class**, which is what the per-sample route's two-state model needed after a start
    /// that put them close together emptied one and reported convergence
    /// (`parameter_prepass_generic.md` §6.5).
    fn start(samples: usize, grids: &Grids, candidate: Candidate, variant: usize) -> Self {
        let cells = grids.free.len();
        let (eps_clean, eps_noisy, w_noisy) = match variant {
            0 => (6e-3, 2e-2, 0.05),
            1 => (5e-4, 2e-1, 0.002),
            _ => (2e-2, 6e-2, 0.20),
        };
        // Deliberately wrong starts, three times the truth's error rate and a flat
        // frequency shape — the shape `parameter_prepass_generic.md` §5.1 uses, so a fit
        // that only works from the answer is caught.
        let weights = match candidate {
            Candidate::CountSpectrum => {
                let n = 2 * samples + 1;
                let mut w = vec![1.0 / (n as f64 * 100.0); n];
                w[0] = 1.0 - (n - 1) as f64 / (n as f64 * 100.0);
                w
            }
            _ => {
                // Most of the mass on "invariant", which cell zero is, and the rest
                // spread — a flat vector would start the fit believing half the genome
                // segregates at a half.
                let mut w = vec![0.01 / cells as f64; cells];
                w[0] = 0.99;
                w
            }
        };
        Self {
            eps_clean,
            eps_noisy,
            w_noisy,
            density: Density {
                p_invariant: 0.95,
                p_fixed_alt: 0.01,
                a: 1.0,
                b: 1.0,
                // A second bump only where one is being fitted, and **placed differently
                // by each start** — two bumps are interchangeable, so a fit that lands
                // in the same place from every start is the evidence that matters.
                second: match (candidate, variant) {
                    (Candidate::Density7, 0) => Some((0.5, 0.5, 0.5)),
                    (Candidate::Density7, 1) => Some((0.2, 3.0, 8.0)),
                    (Candidate::Density7, _) => Some((0.8, 8.0, 3.0)),
                    _ => None,
                },
            },
            inbreeding: vec![0.1; samples],
            weights,
        }
    }
}

/// The objective: `Σ weight · ln L(pattern; θ)`, which is `Q(θ)` in the exact arm and the
/// mean log-likelihood in the drawn one.
/// The whole frequency description a candidate is currently holding, as the pair
/// `ln_mixture` reads: the masses-and-representatives, and the genotype priors they imply.
///
/// **Built once per objective evaluation, never per locus** — the representatives move
/// only when the Beta's shape does.
struct FrequencyView {
    density: Density,
    cells: Vec<(f64, f64)>,
    ln_priors: Vec<Vec<[f64; 3]>>,
}

fn frequency_view(fitted: &Fitted, candidate: Candidate, grids: &Grids) -> FrequencyView {
    match candidate {
        Candidate::Density4 | Candidate::Density7 => {
            let cells = grids.quadrature.under_density(&fitted.density);
            let ln_priors = genotype_priors(&cells, &fitted.inbreeding);
            FrequencyView {
                density: fitted.density,
                cells,
                ln_priors,
            }
        }
        // `grid-free` puts its whole description on the grid, endpoints included, so it
        // has no separate invariant mass: component zero *is* that mass.
        _ => {
            let total: f64 = fitted.weights.iter().sum();
            let cells: Vec<(f64, f64)> = fitted
                .weights
                .iter()
                .zip(&grids.free)
                .map(|(w, &f)| (w / total, f))
                .collect();
            let ln_priors = genotype_priors(&cells, &fitted.inbreeding);
            FrequencyView {
                density: Density {
                    p_invariant: 0.0,
                    p_fixed_alt: 0.0,
                    a: 1.0,
                    b: 1.0,
                    second: None,
                },
                cells,
                ln_priors,
            }
        }
    }
}

fn objective(
    data: &Weighted,
    cache: &[PatternLikelihoods],
    fitted: &Fitted,
    candidate: Candidate,
    grids: &Grids,
    // `shared` says whether `cache` holds only the loud samples. It must agree with how
    // the cache was built, so it is passed rather than inferred.
    shared: bool,
) -> f64 {
    let class_weights = [1.0 - fitted.w_noisy, fitted.w_noisy];
    if candidate == Candidate::CountSpectrum {
        let terms: Vec<f64> = data
            .weights
            .par_iter()
            .enumerate()
            .map(|(i, weight)| {
                weight * ln_pattern_under_count_spectrum(&cache[i], &class_weights, &fitted.weights)
            })
            .collect();
        return terms.iter().sum();
    }
    let view = frequency_view(fitted, candidate, grids);
    // Collected and then summed in order, never reduced in parallel: a tree reduction over
    // floats depends on how the work was split, so the objective would move with the
    // thread count and so would every fitted value (spec §7, determinism).
    let terms: Vec<f64> = if shared {
        // **The path that makes a large panel affordable.** One inbreeding coefficient for
        // everyone means the silent samples are interchangeable, so what a locus costs is
        // a table lookup per depth plus one evaluation per sample that actually showed an
        // alternative read — a dozen or so, whatever the panel size.
        let quiet = QuietTables::new(
            &[fitted.eps_clean, fitted.eps_noisy],
            &view.cells,
            fitted.inbreeding[0],
            data.deepest,
        );
        let ln_priors: Vec<[f64; 3]> = view.ln_priors.iter().map(|per| per[0]).collect();
        data.weights
            .par_iter()
            .enumerate()
            .map(|(i, weight)| {
                weight
                    * ln_mixture_shared(
                        &data.shapes[i],
                        &cache[i],
                        &class_weights,
                        &view.density,
                        &view.cells,
                        &ln_priors,
                        &quiet,
                    )
            })
            .collect()
    } else {
        data.weights
            .par_iter()
            .enumerate()
            .map(|(i, weight)| {
                weight
                    * ln_mixture(
                        &cache[i],
                        &class_weights,
                        &view.density,
                        &view.cells,
                        &view.ln_priors,
                    )
            })
            .collect()
    };
    terms.iter().sum()
}

/// The count form: `Σ_c sfs(c) · P(pattern | c)`, with `P(pattern | c)` from the
/// convolution that keeps the genotypes adding up to `c` (spec §2.1.1).
fn ln_pattern_under_count_spectrum(
    likelihoods: &PatternLikelihoods,
    class_weights: &[f64],
    sfs: &[f64],
) -> f64 {
    let chromosomes = (sfs.len() - 1) as u32;
    let mut per_class = Vec::with_capacity(class_weights.len());
    let total_sfs: f64 = sfs.iter().sum();
    for (class, &class_weight) in class_weights.iter().enumerate() {
        if class_weight <= 0.0 {
            per_class.push(f64::NEG_INFINITY);
            continue;
        }
        let ln_read = &likelihoods.ln_read[class];
        // `h[k]` is the `ln` of the partial sum over the samples seen so far, of the
        // read likelihoods times the number of ways their genotypes place `k` alleles.
        let mut h = vec![f64::NEG_INFINITY; sfs.len()];
        h[0] = 0.0;
        for ln_g in ln_read {
            let mut next = vec![f64::NEG_INFINITY; sfs.len()];
            for (k, &value) in h.iter().enumerate() {
                if !value.is_finite() {
                    continue;
                }
                for genotype in 0..=PLOIDY {
                    let target = k + genotype as usize;
                    if target >= next.len() {
                        continue;
                    }
                    let ways = ln_binomial_coefficient(PLOIDY, genotype);
                    let term = value + ways + ln_g[genotype as usize];
                    next[target] = ln_sum_exp(&[next[target], term]);
                }
            }
            h = next;
        }
        let terms: Vec<f64> = sfs
            .iter()
            .enumerate()
            .filter(|&(_, &w)| w > 0.0)
            .map(|(c, &w)| {
                (w / total_sfs).ln() + h[c] - ln_binomial_coefficient(chromosomes, c as u32)
            })
            .collect();
        per_class.push(class_weight.ln() + ln_sum_exp(&terms));
    }
    ln_sum_exp(&per_class)
}

// ---------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------

/// Golden-section maximisation of a one-dimensional slice, in a transformed coordinate so
/// that a rate never leaves `(0, 1)` and a shape never leaves `(0, ∞)`.
fn climb_scalar(mut score: impl FnMut(f64) -> f64, start: f64, span: f64) -> f64 {
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let (mut lo, mut hi) = (start - span, start + span);
    let mut c = hi - phi * (hi - lo);
    let mut d = lo + phi * (hi - lo);
    let (mut fc, mut fd) = (score(c), score(d));
    for _ in 0..25 {
        if fc < fd {
            lo = c;
            c = d;
            fc = fd;
            d = lo + phi * (hi - lo);
            fd = score(d);
        } else {
            hi = d;
            d = c;
            fd = fc;
            c = hi - phi * (hi - lo);
            fc = score(c);
        }
        if (hi - lo).abs() < 1e-7 {
            break;
        }
    }
    0.5 * (lo + hi)
}

fn logit(p: f64) -> f64 {
    let p = p.clamp(1e-9, 1.0 - 1e-9);
    (p / (1.0 - p)).ln()
}

fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// One expectation-maximization step for a free weight vector: each pattern's posterior
/// over the components, averaged with the pattern's weight. Monotone in the objective,
/// which a coordinate climb over a hundred correlated weights would not be.
fn expectation_step(
    data: &Weighted,
    cache: &[PatternLikelihoods],
    fitted: &Fitted,
    candidate: Candidate,
    grids: &Grids,
) -> Vec<f64> {
    let class_weights = [1.0 - fitted.w_noisy, fitted.w_noisy];
    let view = frequency_view(fitted, candidate, grids);
    let mut accumulated = vec![0.0; fitted.weights.len()];
    let total_weight: f64 = fitted.weights.iter().sum();
    for (i, &pattern_weight) in data.weights.iter().enumerate() {
        let ln_component: Vec<f64> = match candidate {
            Candidate::GridFree => (0..view.cells.len())
                .map(|cell| {
                    let mut terms = Vec::with_capacity(class_weights.len());
                    for (class, &cw) in class_weights.iter().enumerate() {
                        let ln_read = &cache[i].ln_read[class];
                        let mut total = cw.ln();
                        for (sample, ln_g) in ln_read.iter().enumerate() {
                            let ln_prior = &view.ln_priors[cell][sample];
                            total += ln_sum_exp(&[
                                ln_prior[0] + ln_g[0],
                                ln_prior[1] + ln_g[1],
                                ln_prior[2] + ln_g[2],
                            ]);
                        }
                        terms.push(total);
                    }
                    (fitted.weights[cell] / total_weight).max(1e-300).ln() + ln_sum_exp(&terms)
                })
                .collect(),
            _ => ln_count_components(&cache[i], &class_weights, &fitted.weights),
        };
        let normaliser = ln_sum_exp(&ln_component);
        for (component, ln) in ln_component.iter().enumerate() {
            accumulated[component] += pattern_weight * (ln - normaliser).exp();
        }
    }
    let total: f64 = accumulated.iter().sum();
    accumulated.iter().map(|a| a / total).collect()
}

/// `ln` of each allele count's contribution to one pattern — the count form's components,
/// which the expectation step needs separately rather than summed.
fn ln_count_components(
    likelihoods: &PatternLikelihoods,
    class_weights: &[f64],
    sfs: &[f64],
) -> Vec<f64> {
    let chromosomes = (sfs.len() - 1) as u32;
    let total_sfs: f64 = sfs.iter().sum();
    let mut out = vec![f64::NEG_INFINITY; sfs.len()];
    for (class, &class_weight) in class_weights.iter().enumerate() {
        if class_weight <= 0.0 {
            continue;
        }
        let ln_read = &likelihoods.ln_read[class];
        let mut h = vec![f64::NEG_INFINITY; sfs.len()];
        h[0] = 0.0;
        for ln_g in ln_read {
            let mut next = vec![f64::NEG_INFINITY; sfs.len()];
            for (k, &value) in h.iter().enumerate() {
                if !value.is_finite() {
                    continue;
                }
                for genotype in 0..=PLOIDY {
                    let target = k + genotype as usize;
                    if target >= next.len() {
                        continue;
                    }
                    let term =
                        value + ln_binomial_coefficient(PLOIDY, genotype) + ln_g[genotype as usize];
                    next[target] = ln_sum_exp(&[next[target], term]);
                }
            }
            h = next;
        }
        for (c, slot) in out.iter_mut().enumerate() {
            if sfs[c] <= 0.0 {
                continue;
            }
            let term = class_weight.ln() + (sfs[c] / total_sfs).ln() + h[c]
                - ln_binomial_coefficient(chromosomes, c as u32);
            *slot = ln_sum_exp(&[*slot, term]);
        }
    }
    out
}

// ---------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------

/// `Hexp` for a candidate that holds free weights rather than a Beta: the mean of
/// `2f(1−f)` over the fitted weights.
fn expected_heterozygosity_of(fitted: &Fitted, candidate: Candidate, grids: &Grids) -> f64 {
    match candidate {
        Candidate::Density4 | Candidate::Density7 => fitted.density.expected_heterozygosity(),
        Candidate::GridFree => {
            let total: f64 = fitted.weights.iter().sum();
            fitted
                .weights
                .iter()
                .zip(&grids.free)
                .map(|(w, f)| w / total * 2.0 * f * (1.0 - f))
                .sum()
        }
        Candidate::CountSpectrum => {
            // From the panel allele-count spectrum, with the finite-sample correction
            // that turns a count spectrum into a population quantity (spec §5.3).
            let chromosomes = (fitted.weights.len() - 1) as f64;
            let total: f64 = fitted.weights.iter().sum();
            fitted
                .weights
                .iter()
                .enumerate()
                .map(|(c, w)| {
                    let c = c as f64;
                    w / total * 2.0 * c * (chromosomes - c) / (chromosomes * (chromosomes - 1.0))
                })
                .sum()
        }
    }
}

// Eight arguments against clippy's seven, and they are eight columns of one printed
// line: bundling them into a struct would name the same eight fields one indirection
// away from the `println!` that formats them.
#[allow(clippy::too_many_arguments)]
fn report(
    label: &str,
    truth: &Truth,
    fitted: &Fitted,
    candidate: Candidate,
    grids: &Grids,
    score: f64,
    passes: usize,
    elapsed: f64,
) {
    let relative = |got: f64, want: f64| {
        if want.abs() < 1e-12 {
            f64::NAN
        } else {
            100.0 * (got - want) / want
        }
    };
    let hexp_truth = truth.density.expected_heterozygosity();
    let hexp = expected_heterozygosity_of(fitted, candidate, grids);
    let mean_inbreeding = fitted.inbreeding.iter().sum::<f64>() / fitted.inbreeding.len() as f64;
    let truth_inbreeding = truth.inbreeding.iter().sum::<f64>() / truth.inbreeding.len() as f64;
    println!(
        "{label:>10} {:>9} | eps_clean {:>10.3e} {:>+7.1}% | eps_noisy {:>9.3e} {:>+7.1}% | \
         w {:>8.4} {:>+7.1}% | Hexp {:>9.3e} {:>+7.1}% | F {:>6.3} {:>+7.1}% | lnL {:>12.6} | \
         {passes} passes, {elapsed:.1}s",
        candidate.name(),
        fitted.eps_clean,
        relative(fitted.eps_clean, truth.eps_clean),
        fitted.eps_noisy,
        relative(fitted.eps_noisy, truth.eps_noisy),
        fitted.w_noisy,
        relative(fitted.w_noisy, truth.w_noisy),
        hexp,
        relative(hexp, hexp_truth),
        mean_inbreeding,
        relative(mean_inbreeding, truth_inbreeding),
        score,
    );
}

/// How spiky a fitted weight vector is: the share of its mass on the segregating cells
/// that sits in its single largest one. **A deconvolution's failure looks like this** —
/// the maximiser of a nonparametric mixing distribution is discrete, so a value near 1
/// means the fit has replaced a density with a handful of atoms.
fn concentration(weights: &[f64]) -> f64 {
    let segregating: Vec<f64> = weights[1..weights.len().saturating_sub(1)].to_vec();
    let total: f64 = segregating.iter().sum();
    if total <= 0.0 {
        return f64::NAN;
    }
    segregating.iter().copied().fold(0.0, f64::max) / total
}

/// **Where does each parameter stop being estimable as the panel shrinks?** — the number
/// a user needs and nothing currently states (spec §12.5).
///
/// The adopted candidate only, and **one inbreeding coefficient shared by the panel**,
/// since the truth gives every sample the same one. That is what makes fifty samples
/// affordable: a per-sample block is `N` climbs a pass where everything else is seven.
fn panel_sweep(depth: f64, loci: usize) {
    println!(
        "\n=== panel sweep: density-4, {depth} reads a site, {loci} loci, one shared \
         inbreeding coefficient ===\n\
         samples |  eps_clean |  eps_noisy |          w |       Hexp |          F | seconds"
    );
    for &inbreeding in &[0.0, 0.6] {
        println!("--- truth inbreeding {inbreeding} ---");
        for &samples in &[1_usize, 2, 5, 10, 25, 50] {
            let truth = Truth::tomato_like(samples, depth, inbreeding);
            let grids = Grids::new(samples, 24);
            let data = drawn_patterns(&truth, loci, 20_260_812);
            let started = Instant::now();
            // **Three starting points, best taken.** An earlier version of this sweep ran
            // one, and the start turned out to be worth seven percentage points on the
            // clean error rate at ten samples — so a one-start sweep measures the search
            // and reports it as the estimator.
            let (fitted, _, _) = (0..3)
                .map(|variant| {
                    fit_rebuilding(
                        &data,
                        Candidate::Density4,
                        &grids,
                        samples,
                        10,
                        variant,
                        true,
                    )
                })
                .max_by(|a, b| a.1.partial_cmp(&b.1).expect("finite scores"))
                .expect("three starts");
            let relative = |got: f64, want: f64| 100.0 * (got - want) / want;
            let hexp = fitted.density.expected_heterozygosity();
            println!(
                "{samples:>7} | {:>+9.1}% | {:>+9.1}% | {:>+9.1}% | {:>+9.1}% | {:>10.3} | {:>7.1}",
                relative(fitted.eps_clean, truth.eps_clean),
                relative(fitted.eps_noisy, truth.eps_noisy),
                relative(fitted.w_noisy, truth.w_noisy),
                relative(hexp, truth.expected_heterozygosity()),
                fitted.inbreeding[0],
                started.elapsed().as_secs_f64(),
            );
        }
    }
    println!(
        "\nThe inbreeding column is the fitted value itself, not an error: the truth is \
         0.000 in the first block and 0.600 in the second."
    );
}

/// **Is a single Beta enough shape for a panel that is landraces from several regions?**
/// (spec §11.6, §12.9).
///
/// The truth's frequency density has two bumps and the fitted one has room for a single
/// hump, so this is the arm that grades the *model* rather than the arithmetic. What it
/// reports is which parameters feel it: the noise rates read the density only as a
/// background and should barely move, while `Hexp` is read straight off it.
fn structured_panel(samples: usize, depth: f64, loci: usize) {
    println!(
        "\n=== a frequency density the fitted shape cannot reach: two components at 0.2 \
         and 0.8, {samples} samples, {depth} reads a site, {loci} loci ==="
    );
    println!(
        "Each row is one starting point. **The spread down a block is the number that \
         decides whether a second bump is worth having**: a shape that lands in the same \
         place from every start is describing the data, and one that does not is \
         describing where it began."
    );
    let grids = Grids::new(samples, 24);
    for (truth_label, truth) in [
        ("one bump", Truth::tomato_like(samples, depth, 0.0)),
        ("two bumps", Truth::structured(samples, depth, 0.0)),
    ] {
        let data = drawn_patterns(&truth, loci, 20_260_812);
        let hexp_truth = truth.expected_heterozygosity();
        println!(
            "\n--- truth: {truth_label}, Hexp {hexp_truth:.3e} ---\n\
             {:>10} {:>5} | {:>10} | {:>10} | {:>10} |       lnL",
            "fitted", "start", "eps_clean", "Hexp", "F"
        );
        for candidate in [Candidate::Density4, Candidate::Density7] {
            let mut recovered = Vec::new();
            for variant in 0..3 {
                let (fitted, score, _) =
                    fit_rebuilding(&data, candidate, &grids, samples, 10, variant, true);
                let relative = |got: f64, want: f64| 100.0 * (got - want) / want;
                let hexp = fitted.density.expected_heterozygosity();
                recovered.push(hexp);
                println!(
                    "{:>10} {variant:>5} | {:>+9.1}% | {:>+9.1}% | {:>10.3} | {score:>12.6}",
                    candidate.name(),
                    relative(fitted.eps_clean, truth.eps_clean),
                    relative(hexp, hexp_truth),
                    fitted.inbreeding[0],
                );
            }
            let lowest = recovered.iter().copied().fold(f64::INFINITY, f64::min);
            let highest = recovered.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            println!(
                "{:>10} {:>5} |            | Hexp spans {:>+6.1}% of its own value across \
                 the three starts",
                "",
                "",
                100.0 * (highest - lowest) / hexp_truth
            );
        }
    }
}

/// **The fast path and the slow path must return the same number**, and this is where
/// that is checked rather than asserted in a comment.
///
/// A fast path that quietly disagreed would not crash: it would return a slightly
/// different likelihood, the fit would land somewhere slightly different, and the result
/// would be reported as a property of the estimator. So every run that uses the fast path
/// prints this first.
fn check_fast_path_agrees(grids: &Grids) {
    let truth = Truth::tomato_like(6, 3.0, 0.4);
    let data = drawn_patterns(&truth, 4_000, 20_260_812);
    let fitted = Fitted::at_truth(&truth, grids);
    let classes = [fitted.eps_clean, fitted.eps_noisy];
    let full: Vec<PatternLikelihoods> = data
        .patterns
        .as_ref()
        .expect("six samples keeps its patterns")
        .iter()
        .map(|p| PatternLikelihoods::new(p, &classes))
        .collect();
    let loud: Vec<PatternLikelihoods> = data
        .shapes
        .iter()
        .map(|s| s.loud_likelihoods(&classes))
        .collect();
    let slow = objective(&data, &full, &fitted, Candidate::Density4, grids, false);
    let fast = objective(&data, &loud, &fitted, Candidate::Density4, grids, true);
    println!(
        "fast path against slow path: {slow:.12} vs {fast:.12} — {}",
        if (slow - fast).abs() < 1e-9 {
            "agree"
        } else {
            "DISAGREE, and nothing below can be trusted"
        }
    );
}

/// **How do the errors respond to the number of sampled loci?**
///
/// The knob's own experiment ([`parameter_prepass_joint_loci.md`] §4.3): refit at each
/// budget on the same truth and report every parameter's error separately, because they do
/// not degrade together — the noise rates go as the square root of the budget while the
/// frequency density is carried by the few loci that actually segregate, and one
/// "precision" column would call a budget adequate that had already lost the density.
///
/// **The segregating count is reported beside each budget**, since it is what the density
/// and the diversity actually rest on, and it is an outcome of the panel rather than a
/// target.
fn locus_sweep(samples: usize, depth: f64, largest: usize) {
    let grids = Grids::new(samples, 24);
    check_fast_path_agrees(&grids);
    println!(
        "\n=== locus sweep: density-4, {samples} samples, {depth} reads a site, one shared \
         inbreeding coefficient, three starts ===\n\
            loci | in population | in panel |  eps_clean |  eps_noisy |          w |       Hexp |      F | seconds"
    );
    for &inbreeding in &[0.0, 0.6] {
        println!("--- truth inbreeding {inbreeding} ---");
        let truth = Truth::tomato_like(samples, depth, inbreeding);
        let mut budget = largest;
        let mut budgets = Vec::new();
        while budget >= 2_000 {
            budgets.push(budget);
            budget /= 4;
        }
        budgets.reverse();
        for &loci in &budgets {
            // A fresh draw per budget rather than a prefix of one draw, so the rows are
            // independent measurements and a lucky draw cannot make a whole column look
            // good.
            let data = drawn_patterns(&truth, loci, 20_260_812 + loci as u64);
            let started = Instant::now();
            let (fitted, _, _) = (0..3)
                .map(|variant| {
                    fit_rebuilding(
                        &data,
                        Candidate::Density4,
                        &grids,
                        samples,
                        10,
                        variant,
                        true,
                    )
                })
                .max_by(|a, b| a.1.partial_cmp(&b.1).expect("finite scores"))
                .expect("three starts");
            let relative = |got: f64, want: f64| 100.0 * (got - want) / want;
            println!(
                "{loci:>8} | {:>13} | {:>8} | {:>+9.1}% | {:>+9.1}% | {:>+9.1}% | {:>+9.1}% | {:>6.3} | {:>7.1}",
                data.segregating,
                data.panel_segregating,
                relative(fitted.eps_clean, truth.eps_clean),
                relative(fitted.eps_noisy, truth.eps_noisy),
                relative(fitted.w_noisy, truth.w_noisy),
                relative(
                    fitted.density.expected_heterozygosity(),
                    truth.expected_heterozygosity()
                ),
                fitted.inbreeding[0],
                started.elapsed().as_secs_f64(),
            );
        }
    }
    println!(
        "\nThe inbreeding column is the fitted value, not an error: the truth is 0.000 in \
         the first block and 0.600 in the second."
    );
}

fn main() {
    let mut args = std::env::args().skip(1);
    let arm = args.next().unwrap_or_else(|| "both".to_string());
    let samples: usize = args
        .next()
        .map_or(2, |a| a.parse().expect("a sample count"));
    let depth: f64 = args
        .next()
        .map_or(3.0, |a| a.parse().expect("a mean depth"));
    let loci: usize = args
        .next()
        .map_or(200_000, |a| a.parse().expect("a locus count"));

    if arm == "loci" {
        locus_sweep(if samples == 0 { 50 } else { samples }, depth, loci);
        return;
    }
    if arm == "sweep" {
        panel_sweep(depth, loci);
        return;
    }
    if arm == "structured" {
        structured_panel(samples, depth, loci);
        return;
    }

    // Twenty-four quadrature cells for the Beta, which is accuracy; the free grid is
    // `2N + 1` frequencies, which is freedom (see `Grids`).
    let grids = Grids::new(samples, 24);
    let candidates = [
        Candidate::Density4,
        Candidate::Density7,
        Candidate::GridFree,
        Candidate::CountSpectrum,
    ];

    for &inbreeding in &[0.0, 0.6] {
        let truth = Truth::tomato_like(samples, depth, inbreeding);
        println!(
            "\n=== truth: {samples} samples, {depth} reads a site, inbreeding {inbreeding} ===\n\
             eps_clean {:.3e}  eps_noisy {:.3e}  w {:.4}  Hexp {:.3e}  \
             (invariant {:.3}, fixed-alt {:.3}, Beta({}, {}))",
            truth.eps_clean,
            truth.eps_noisy,
            truth.w_noisy,
            truth.density.expected_heterozygosity(),
            truth.density.p_invariant,
            truth.density.p_fixed_alt,
            truth.density.a,
            truth.density.b,
        );

        let mut arms: Vec<(&str, Weighted)> = Vec::new();
        if arm != "drawn" {
            // The pattern space is `(depths × counts)^samples`, so the cap is what keeps
            // the exact arm finite. At three reads a site a cap of 6 leaves under one
            // locus in a thousand outside it.
            let depth_cap = 6;
            let exact = exact_patterns(&truth, depth_cap, &grids);
            let mass: f64 = exact.weights.iter().sum();
            println!(
                "exact arm: {} patterns at depths up to {depth_cap}, holding {:.4} of the \
                 probability",
                exact.shapes.len(),
                mass
            );
            arms.push(("exact", exact));
        }
        if arm != "exact" {
            arms.push(("drawn", drawn_patterns(&truth, loci, 20_260_812)));
        }

        for (label, data) in &arms {
            // The objective at the truth itself. **A fitted score below this is an
            // optimiser that stopped early, not a biased estimator** — the two are
            // indistinguishable in a table of parameter errors, and this is the one
            // number that tells them apart.
            let at_truth = Fitted::at_truth(&truth, &grids);
            let cache: Vec<PatternLikelihoods> = data
                .patterns
                .as_ref()
                .expect("the small-panel arms keep their patterns")
                .iter()
                .map(|p| PatternLikelihoods::new(p, &[truth.eps_clean, truth.eps_noisy]))
                .collect();
            let truth_score =
                objective(data, &cache, &at_truth, Candidate::Density4, &grids, false);
            println!("{label:>10} {:>9} | lnL {truth_score:>12.6}", "at truth");

            for &candidate in &candidates {
                if candidate == Candidate::CountSpectrum && inbreeding > 0.0 {
                    println!(
                        "{label:>10} {:>9} | not comparable: the count form cannot carry a \
                         per-sample inbreeding coefficient (spec §2.1.1)",
                        candidate.name()
                    );
                    continue;
                }
                let started = Instant::now();
                // **Three starting points, spanning the separation between the clean and
                // the noisy class** — the shape the per-sample route's two-state model
                // needed, and the best of the three is the estimate.
                // Three starts on the exact arm, where a pass is cheap and a trapped climb
                // is what is being looked for; one on the drawn arm, where the same
                // question has already been answered and the cost is fifteen times higher.
                let starts = if *label == "exact" { 3 } else { 1 };
                let (fitted, score, passes) = (0..starts)
                    .map(|variant| {
                        fit_rebuilding(data, candidate, &grids, samples, 10, variant, false)
                    })
                    .max_by(|a, b| a.1.partial_cmp(&b.1).expect("finite scores"))
                    .expect("at least one start");
                report(
                    label,
                    &truth,
                    &fitted,
                    candidate,
                    &grids,
                    score,
                    passes,
                    started.elapsed().as_secs_f64(),
                );
                if candidate != Candidate::Density4 {
                    println!(
                        "{:>21} largest single component holds {:.1}% of the segregating mass",
                        "",
                        100.0 * concentration(&fitted.weights)
                    );
                }
            }
        }
    }
}

/// Alternate: the noise numbers, then the frequency description, then each sample's
/// inbreeding (spec §3.3). Coordinate-wise, because the objective is smooth and every
/// coordinate has a natural transform that keeps it legal.
///
/// **The read likelihoods depend on the error rates being tried**, so they are rebuilt
/// whenever those move rather than cached across the fit — which is what makes the exact
/// arm's pattern cap matter.
fn fit_rebuilding(
    data: &Weighted,
    candidate: Candidate,
    grids: &Grids,
    samples: usize,
    passes: usize,
    variant: usize,
    shared_inbreeding: bool,
) -> (Fitted, f64, usize) {
    // The cache must match the path the objective will take: only the loud samples where
    // the panel shares one inbreeding coefficient, every sample where it does not.
    let build = |eps_clean: f64, eps_noisy: f64| -> Vec<PatternLikelihoods> {
        let classes = [eps_clean, eps_noisy];
        if shared_inbreeding {
            data.shapes
                .iter()
                .map(|s| s.loud_likelihoods(&classes))
                .collect()
        } else {
            data.patterns
                .as_ref()
                .expect("a per-sample inbreeding fit needs the whole pattern")
                .iter()
                .map(|p| PatternLikelihoods::new(p, &classes))
                .collect()
        }
    };
    let mut fitted = Fitted::start(samples, grids, candidate, variant);
    let mut previous = f64::NEG_INFINITY;
    let mut used = 0;
    for pass in 0..passes {
        used = pass + 1;
        for which in 0..3 {
            let base = fitted.clone();
            let start = match which {
                0 => logit(base.eps_clean),
                1 => logit(base.eps_noisy),
                _ => logit(base.w_noisy),
            };
            let best = climb_scalar(
                |x| {
                    let mut trial = base.clone();
                    match which {
                        0 => trial.eps_clean = expit(x),
                        1 => trial.eps_noisy = expit(x),
                        _ => trial.w_noisy = expit(x),
                    }
                    let cache = build(trial.eps_clean, trial.eps_noisy);
                    objective(data, &cache, &trial, candidate, grids, shared_inbreeding)
                },
                start,
                3.0,
            );
            match which {
                0 => fitted.eps_clean = expit(best),
                1 => fitted.eps_noisy = expit(best),
                _ => fitted.w_noisy = expit(best),
            }
        }
        let cache = build(fitted.eps_clean, fitted.eps_noisy);
        match candidate {
            Candidate::Density4 | Candidate::Density7 => {
                let numbers = if candidate == Candidate::Density7 {
                    7
                } else {
                    4
                };
                for which in 0..numbers {
                    let base = fitted.clone();
                    let second = base.density.second.unwrap_or((0.5, 1.0, 1.0));
                    let start = match which {
                        0 => logit(base.density.p_invariant),
                        1 => logit(base.density.p_fixed_alt),
                        2 => base.density.a.ln(),
                        3 => base.density.b.ln(),
                        4 => logit(second.0),
                        5 => second.1.ln(),
                        _ => second.2.ln(),
                    };
                    let best = climb_scalar(
                        |x| {
                            let mut trial = base.clone();
                            let mut second = second;
                            match which {
                                0 => trial.density.p_invariant = expit(x),
                                1 => trial.density.p_fixed_alt = expit(x),
                                2 => trial.density.a = x.exp(),
                                3 => trial.density.b = x.exp(),
                                4 => second.0 = expit(x),
                                5 => second.1 = x.exp(),
                                _ => second.2 = x.exp(),
                            }
                            if which >= 4 {
                                trial.density.second = Some(second);
                            }
                            if trial.density.p_segregating() <= 1e-9 {
                                return f64::NEG_INFINITY;
                            }
                            objective(data, &cache, &trial, candidate, grids, shared_inbreeding)
                        },
                        start,
                        3.0,
                    );
                    let mut second = second;
                    match which {
                        0 => fitted.density.p_invariant = expit(best),
                        1 => fitted.density.p_fixed_alt = expit(best),
                        2 => fitted.density.a = best.exp(),
                        3 => fitted.density.b = best.exp(),
                        4 => second.0 = expit(best),
                        5 => second.1 = best.exp(),
                        _ => second.2 = best.exp(),
                    }
                    if which >= 4 {
                        fitted.density.second = Some(second);
                    }
                }
            }
            _ => {
                for _ in 0..20 {
                    fitted.weights = expectation_step(data, &cache, &fitted, candidate, grids);
                }
            }
        }
        if candidate != Candidate::CountSpectrum {
            if shared_inbreeding {
                // **One number for the whole panel.** Only legitimate where the truth
                // gave every sample the same coefficient, which the panel sweep does —
                // and it is what makes that sweep affordable, since a per-sample block is
                // `N` climbs a pass and everything else is seven.
                let base = fitted.clone();
                let best = climb_scalar(
                    |x| {
                        let mut trial = base.clone();
                        trial.inbreeding.iter_mut().for_each(|f| *f = expit(x));
                        objective(data, &cache, &trial, candidate, grids, shared_inbreeding)
                    },
                    logit(base.inbreeding[0].clamp(1e-6, 1.0 - 1e-6)),
                    3.0,
                );
                fitted.inbreeding.iter_mut().for_each(|f| *f = expit(best));
            } else {
                for sample in 0..samples {
                    let base = fitted.clone();
                    let best = climb_scalar(
                        |x| {
                            let mut trial = base.clone();
                            trial.inbreeding[sample] = expit(x);
                            objective(data, &cache, &trial, candidate, grids, shared_inbreeding)
                        },
                        logit(base.inbreeding[sample].clamp(1e-6, 1.0 - 1e-6)),
                        3.0,
                    );
                    fitted.inbreeding[sample] = expit(best);
                }
            }
        }
        let score = objective(data, &cache, &fitted, candidate, grids, shared_inbreeding);
        if (score - previous).abs() < 1e-10 {
            previous = score;
            break;
        }
        previous = score;
    }
    (fitted, previous, used)
}
