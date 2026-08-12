//! **Does the repeat-tract half of the joint fit recover what it was given, and is the
//! "how monomorphic are loci" number identified at all?**
//!
//! `spec/parameter_prepass_joint_fit.md` §4.1 and §4.2 were written in one sitting by
//! analogy with the ordinary-position path, with no program behind them. They say two
//! things:
//!
//! - **§4.1** — a locus's length frequencies are a latent vector drawn from a Dirichlet
//!   whose mean is the stratum's own length spectrum and whose **concentration `κ` says how
//!   monomorphic loci are**. A large `κ` is meant to return the per-stratum model exactly,
//!   so that "per locus or per stratum" becomes a comparison of one fitted number.
//! - **§4.2** — the integral over that vector is done by summing over the configurations a
//!   locus can be in: fixed for one length, or segregating two.
//!
//! **Those two are not the same model, and this program is where that stops being an
//! argument.** A Dirichlet is a density on the interior of the simplex, and every
//! configuration §4.2 sums over is a face of it — a set the Dirichlet gives probability
//! zero. Where the two disagree is exactly where §4.1's large-`κ` limit lives: a locus
//! carrying all of the stratum's lengths at the stratum's own frequencies is in the
//! interior, and §4.2 cannot represent it.
//!
//! So three descriptions of a locus are fitted here against one truth, and none of them is
//! assumed:
//!
//! - **`per-stratum`** — every locus's length frequencies *are* the stratum's spectrum. This
//!   is what the per-sample route does today, and it is §4.1's `κ → ∞` limit.
//! - **`dirichlet`** — §4.1 taken literally, integrated over the whole simplex by nested
//!   quantile quadrature. Fits `κ`.
//! - **`faces`** — §4.2's support: a mass on each *fixed for one length* corner and a mass
//!   on each *segregating two lengths* edge, with the frequency along an edge on a
//!   quadrature. **Stated as its own model rather than as an approximation of the
//!   Dirichlet**, because a face has no Dirichlet density to inherit. Fits a monomorphic
//!   share and a shape number in place of `κ`, which is the same four-number shape §2.1.2
//!   fits on the ordinary-position path.
//!
//! ## What a read does
//!
//! A read reports the tract length of one of the sample's two alleles, or slips. The three
//! slippage numbers are `parameter_prepass_ssr.md` §3's: **how often** a read slips, **which
//! way** — reads showing a shorter tract outnumber longer ones by 4.9 times at tomato
//! dinucleotides — and **how far**, geometrically. A step that would land outside the
//! recorded range saturates into the end class, as the record does.
//!
//! ## What is measured, and what is not
//!
//! The truth is drawn, so every error below carries sampling scatter as well as bias; the
//! `seeds` mode refits the same truth from several draws so the two can be told apart. The
//! inbreeding coefficient is **supplied**, as the spec has it arriving from the
//! ordinary-position path, so nothing here says whether it could be fitted jointly.
//!
//! ```text
//! ng_joint_str_harness [fit|faces-truth|kappa|classes|seeds] [samples] [depth] [loci] [classes] [kappa]
//! ```

use std::time::Instant;

// ---------------------------------------------------------------------
// Small numerics
// ---------------------------------------------------------------------

/// A reproducible stream, the same one the other harnesses use.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1_u64 << 53) as f64
    }

    /// A Gamma(shape, 1) draw — Marsaglia and Tsang, with the shape-below-one boost.
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
            if u.ln() < 0.5 * x * x + d - d * v + d * (v.ln()) {
                return d * v;
            }
        }
    }

    fn normal(&mut self) -> f64 {
        // Box–Muller; one of the pair is enough and the other costs nothing to drop.
        let u1 = self.uniform().max(1e-300);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// A draw from `Dirichlet(alpha)`.
    fn dirichlet(&mut self, alpha: &[f64]) -> Vec<f64> {
        let draws: Vec<f64> = alpha.iter().map(|a| self.gamma(*a).max(1e-300)).collect();
        let total: f64 = draws.iter().sum();
        draws.into_iter().map(|d| d / total).collect()
    }

    /// Which class, given the weights (which must sum to one).
    fn categorical(&mut self, weights: &[f64]) -> usize {
        let mut u = self.uniform();
        for (index, weight) in weights.iter().enumerate() {
            u -= weight;
            if u <= 0.0 {
                return index;
            }
        }
        weights.len() - 1
    }
}

fn ln_sum_exp(values: &[f64]) -> f64 {
    let largest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !largest.is_finite() {
        return largest;
    }
    largest + values.iter().map(|v| (v - largest).exp()).sum::<f64>().ln()
}

fn ln_gamma(x: f64) -> f64 {
    // Lanczos, g = 7, n = 9 — the same coefficients the other harnesses carry.
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
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
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let mut series = COEFFICIENTS[0];
    for (index, coefficient) in COEFFICIENTS.iter().enumerate().skip(1) {
        series += coefficient / (x + index as f64);
    }
    let t = x + 7.5;
    0.5 * std::f64::consts::TAU.ln() + (x + 0.5) * t.ln() - t + series.ln()
}

/// `I_x(a, b)`, by its continued fraction — enough for a Beta quantile by bisection.
fn regularised_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front = (a * x.ln() + b * (1.0 - x).ln() + ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b))
        .exp()
        / a;
    if x > (a + 1.0) / (a + b + 2.0) {
        return 1.0 - regularised_incomplete_beta(1.0 - x, b, a);
    }
    let (mut f, mut c, mut d) = (1.0_f64, 1.0_f64, 0.0_f64);
    for index in 0..=200 {
        let m = index / 2;
        let numerator = if index == 0 {
            1.0
        } else if index % 2 == 0 {
            let m = m as f64;
            (m * (b - m) * x) / ((a + 2.0 * m - 1.0) * (a + 2.0 * m))
        } else {
            let m = m as f64;
            -((a + m) * (a + b + m) * x) / ((a + 2.0 * m) * (a + 2.0 * m + 1.0))
        };
        d = 1.0 + numerator * d;
        if d.abs() < 1e-30 {
            d = 1e-30;
        }
        d = 1.0 / d;
        c = 1.0 + numerator / c;
        if c.abs() < 1e-30 {
            c = 1e-30;
        }
        let step = c * d;
        f *= step;
        if (1.0 - step).abs() < 1e-12 {
            break;
        }
    }
    front * (f - 1.0)
}

/// The `p`-th quantile of `Beta(a, b)`, by bisection on the CDF.
fn beta_quantile(p: f64, a: f64, b: f64) -> f64 {
    let (mut low, mut high) = (0.0_f64, 1.0_f64);
    for _ in 0..60 {
        let middle = 0.5 * (low + high);
        if regularised_incomplete_beta(middle, a, b) < p {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * (low + high)
}

/// `n` equally probable nodes of `Beta(a, b)`, each standing for `1/n` of its mass.
///
/// **Equal probability rather than equal spacing**, because a Dirichlet with a small
/// concentration piles its mass against the ends: an evenly spaced grid puts almost every
/// node where there is nothing and misses the spike, and the spike is the whole phenomenon
/// `κ` describes.
fn beta_nodes(a: f64, b: f64, n: usize) -> Vec<f64> {
    (0..n)
        .map(|index| beta_quantile((index as f64 + 0.5) / n as f64, a, b))
        .collect()
}

// ---------------------------------------------------------------------
// What a read does
// ---------------------------------------------------------------------

/// The three numbers `parameter_prepass_ssr.md` §3 fits, per read group and stratum.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Slippage {
    /// How often a read reports a length other than its allele's.
    level: f64,
    /// Of the reads that slip, the share showing a **shorter** tract. Tomato dinucleotides
    /// sit near 0.83 — 2,438 shorter against 501 longer.
    down_share: f64,
    /// How fast two-repeat slips fall off against one-repeat slips.
    fall_off: f64,
}

impl Slippage {
    /// `P(read reports class c | allele is class a)`, over `classes` length classes.
    ///
    /// **A step past the end saturates into the end class**, as the record's end buckets do
    /// — dropping it would make the level fitted from an edge allele too small, and wrapping
    /// it would put mass at the far end of the tract.
    fn read_probabilities(&self, allele: usize, classes: usize) -> Vec<f64> {
        let mut out = vec![0.0; classes];
        out[allele] = 1.0 - self.level;
        for step in 1..=classes {
            // The geometric weight of a step of exactly `step` repeats, and — at the last
            // one — of every step at least that long, which is what the end class absorbs.
            let weight = if step == classes {
                self.fall_off.powi(step as i32 - 1)
            } else {
                (1.0 - self.fall_off) * self.fall_off.powi(step as i32 - 1)
            };
            let shorter = allele.saturating_sub(step);
            out[shorter] += self.level * self.down_share * weight;
            let longer = (allele + step).min(classes - 1);
            out[longer] += self.level * (1.0 - self.down_share) * weight;
        }
        out
    }
}

/// One sample at one locus: how many reads reported each length class.
type ReadCounts = Vec<u32>;

/// `ln P(this sample's reads | genotype)`, the multinomial coefficient dropped — it is the
/// same for every genotype and every candidate, so it cancels out of both the fit and the
/// comparison.
fn ln_reads_given_genotype(counts: &[u32], per_allele: &[Vec<f64>], first: usize, second: usize) -> f64 {
    let mut total = 0.0;
    for (class, count) in counts.iter().enumerate() {
        if *count == 0 {
            continue;
        }
        let probability = 0.5 * (per_allele[first][class] + per_allele[second][class]);
        total += f64::from(*count) * probability.max(1e-300).ln();
    }
    total
}

// ---------------------------------------------------------------------
// The truth, and the data it draws
// ---------------------------------------------------------------------

/// Where a locus's length frequencies come from, in the truth.
///
/// **Both are fitted by both candidates**, because a comparison in which one candidate is
/// asked to recover its own family and the other is not measures the family and not the
/// candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TruthModel {
    /// §4.1's: a Dirichlet over the whole simplex.
    Dirichlet,
    /// §4.2's: fixed for one length, or segregating two.
    Faces,
}

#[derive(Debug, Clone)]
struct Truth {
    model: TruthModel,
    classes: usize,
    /// The stratum's own length spectrum — the Dirichlet's mean.
    spectrum: Vec<f64>,
    /// How monomorphic loci are. Small is monomorphic.
    kappa: f64,
    slippage: Slippage,
    inbreeding: f64,
    samples: usize,
    depth: u32,
}

impl Truth {
    /// A stratum shaped like a real one: most loci at the reference length, one repeat
    /// either side carrying most of the rest.
    fn default_for(classes: usize, samples: usize, depth: u32, kappa: f64) -> Self {
        let middle = classes / 2;
        let mut spectrum = vec![0.0; classes];
        for (class, weight) in spectrum.iter_mut().enumerate() {
            let distance = (class as i32 - middle as i32).abs();
            *weight = 0.55_f64.powi(distance);
        }
        let total: f64 = spectrum.iter().sum();
        for weight in &mut spectrum {
            *weight /= total;
        }
        Self {
            model: TruthModel::Dirichlet,
            classes,
            spectrum,
            kappa,
            slippage: Slippage {
                level: 0.08,
                down_share: 0.83,
                fall_off: 0.25,
            },
            inbreeding: 0.4,
            samples,
            depth,
        }
    }

    /// One locus's length frequencies, under whichever model this truth carries.
    ///
    /// The faces truth reuses `κ` rather than taking a parameter of its own: its monomorphic
    /// share is `1 / (1 + κ)`, so it runs from nearly every locus fixed at small `κ` to
    /// nearly every locus segregating at large `κ`, in the same direction the Dirichlet's
    /// concentration does. **It is not the same distribution** — that is the whole point of
    /// running it — only the same knob pointing the same way.
    fn draw_frequencies(&self, rng: &mut Rng) -> Vec<f64> {
        match self.model {
            TruthModel::Dirichlet => {
                let alpha: Vec<f64> = self.spectrum.iter().map(|q| self.kappa * q).collect();
                rng.dirichlet(&alpha)
            }
            TruthModel::Faces => {
                let mut frequencies = vec![0.0; self.classes];
                let monomorphic = 1.0 / (1.0 + self.kappa);
                if rng.uniform() < monomorphic {
                    frequencies[rng.categorical(&self.spectrum)] = 1.0;
                    return frequencies;
                }
                let first = rng.categorical(&self.spectrum);
                let mut second = rng.categorical(&self.spectrum);
                while second == first {
                    second = rng.categorical(&self.spectrum);
                }
                // A symmetric Beta(1/2, 1/2), drawn from two Gammas.
                let a = rng.gamma(0.5);
                let b = rng.gamma(0.5);
                let share = a / (a + b).max(1e-300);
                frequencies[first] = share;
                frequencies[second] = 1.0 - share;
                frequencies
            }
        }
    }

    /// One locus's reads for every sample, and the frequencies they were drawn under.
    fn draw_locus(&self, rng: &mut Rng) -> (Vec<ReadCounts>, Vec<f64>) {
        let frequencies = self.draw_frequencies(rng);
        let per_allele: Vec<Vec<f64>> = (0..self.classes)
            .map(|allele| self.slippage.read_probabilities(allele, self.classes))
            .collect();

        let mut panel = Vec::with_capacity(self.samples);
        for _ in 0..self.samples {
            let first = rng.categorical(&frequencies);
            // Inbreeding: with probability F the two copies are identical by descent.
            let second = if rng.uniform() < self.inbreeding {
                first
            } else {
                rng.categorical(&frequencies)
            };
            let mut counts = vec![0_u32; self.classes];
            for _ in 0..self.depth {
                let allele = if rng.uniform() < 0.5 { first } else { second };
                counts[rng.categorical(&per_allele[allele])] += 1;
            }
            panel.push(counts);
        }
        (panel, frequencies)
    }

    fn draw(&self, loci: usize, seed: u64) -> Vec<Vec<ReadCounts>> {
        let mut rng = Rng(seed);
        (0..loci).map(|_| self.draw_locus(&mut rng).0).collect()
    }
}

// ---------------------------------------------------------------------
// The three descriptions of a locus
// ---------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Candidate {
    /// Every locus carries the stratum's own spectrum — `κ → ∞`, and what the per-sample
    /// route does today.
    PerStratum,
    /// §4.1 taken literally: a Dirichlet over the whole simplex.
    Dirichlet,
    /// §4.2's support: fixed for one length, or segregating two.
    Faces,
    /// **Not in the spec.** The same Dirichlet as [`Dirichlet`](Self::Dirichlet), integrated
    /// over a fixed low-discrepancy point set instead of a full tensor grid — so the cost is
    /// the number of points rather than `nodes^(classes − 1)`, and it does not grow with the
    /// number of length classes at all.
    DirichletPoints,
}

impl Candidate {
    fn name(self) -> &'static str {
        match self {
            Self::PerStratum => "per-stratum",
            Self::Dirichlet => "dirichlet",
            Self::Faces => "faces",
            Self::DirichletPoints => "dirichlet-pts",
        }
    }
}

/// What a fit holds. `kappa` is the Dirichlet's; `monomorphic_share` and `edge_shape` are
/// the faces model's replacement for it, and each candidate reads only what it uses.
#[derive(Debug, Clone)]
struct Parameters {
    spectrum: Vec<f64>,
    kappa: f64,
    monomorphic_share: f64,
    edge_shape: f64,
    slippage: Slippage,
}

impl Parameters {
    fn start(classes: usize, level: f64, kappa: f64, monomorphic_share: f64) -> Self {
        Self {
            spectrum: vec![1.0 / classes as f64; classes],
            kappa,
            monomorphic_share,
            edge_shape: 0.5,
            slippage: Slippage {
                level,
                down_share: 0.5,
                fall_off: 0.3,
            },
        }
    }
}

/// One point of the simplex the locus's frequencies are integrated over, with the weight it
/// stands for.
struct Support {
    frequencies: Vec<Vec<f64>>,
    ln_weights: Vec<f64>,
}

/// How many quantile nodes each stick-breaking dimension gets.
const DIRICHLET_NODES: usize = 24;
/// How many quadrature nodes the frequency along an edge gets.
const EDGE_NODES: usize = 12;

/// The support the candidate integrates a locus's frequencies over.
///
/// **This is where the three candidates actually differ**, and it is the whole comparison:
/// one point, a grid over the interior, or the corners and edges.
fn support(candidate: Candidate, parameters: &Parameters, classes: usize) -> Support {
    match candidate {
        Candidate::PerStratum => Support {
            frequencies: vec![parameters.spectrum.clone()],
            ln_weights: vec![0.0],
        },
        Candidate::Dirichlet => dirichlet_support(parameters, classes),
        Candidate::DirichletPoints => dirichlet_point_support(parameters, classes),
        Candidate::Faces => faces_support(parameters, classes),
    }
}

/// Stick-breaking quadrature: `π₁ = x₁`, `π₂ = (1 − x₁)x₂`, and so on, each stick a Beta
/// whose nodes carry equal probability.
///
/// The node count is `DIRICHLET_NODES^(classes − 1)`, which is what makes this candidate the
/// expensive one — 576 points at three classes and 331,776 at five. **That cost is §4.2's
/// entire motivation**, and it is reported rather than asserted.
fn dirichlet_support(parameters: &Parameters, classes: usize) -> Support {
    let alpha: Vec<f64> = parameters
        .spectrum
        .iter()
        .map(|q| (parameters.kappa * q).max(1e-3))
        .collect();
    let mut points: Vec<(Vec<f64>, f64)> = vec![(Vec::new(), 1.0)];
    for stick in 0..classes - 1 {
        let a = alpha[stick];
        let b: f64 = alpha[stick + 1..].iter().sum();
        let nodes = beta_nodes(a, b, DIRICHLET_NODES);
        let mut next = Vec::with_capacity(points.len() * DIRICHLET_NODES);
        for (prefix, remaining) in &points {
            for node in &nodes {
                let mut extended = prefix.clone();
                extended.push(remaining * node);
                next.push((extended, remaining * (1.0 - node)));
            }
        }
        points = next;
    }
    let count = points.len() as f64;
    let mut frequencies = Vec::with_capacity(points.len());
    for (mut prefix, remaining) in points {
        prefix.push(remaining);
        frequencies.push(prefix);
    }
    Support {
        ln_weights: vec![-count.ln(); frequencies.len()],
        frequencies,
    }
}

/// How many low-discrepancy points [`Candidate::DirichletPoints`] integrates over.
const QMC_POINTS: usize = 256;

/// The `index`-th term of the van der Corput sequence in `base` — the coordinate of a Halton
/// point.
fn van_der_corput(mut index: usize, base: usize) -> f64 {
    let (mut fraction, mut result) = (1.0 / base as f64, 0.0);
    while index > 0 {
        result += (index % base) as f64 * fraction;
        index /= base;
        fraction /= base as f64;
    }
    result
}

/// The same Dirichlet, over a **fixed** Halton point set pushed through the stick-breaking
/// quantiles.
///
/// **The cost stops depending on the number of length classes.** A tensor grid needs
/// `nodes^(classes − 1)` points — 576 at three classes and 331,776 at five — where this needs
/// [`QMC_POINTS`] whatever the class count is.
///
/// **The points are fixed and the quantile map is continuous**, so the objective is a smooth
/// function of the concentration rather than a jittery one: the same uniforms are reused at
/// every value the climb tries, which is what stops the search chasing sampling noise.
fn dirichlet_point_support(parameters: &Parameters, classes: usize) -> Support {
    const BASES: [usize; 8] = [2, 3, 5, 7, 11, 13, 17, 19];
    let alpha: Vec<f64> = parameters
        .spectrum
        .iter()
        .map(|q| (parameters.kappa * q).max(1e-3))
        .collect();
    let mut frequencies = Vec::with_capacity(QMC_POINTS);
    for point in 0..QMC_POINTS {
        let mut remaining = 1.0;
        let mut piece = Vec::with_capacity(classes);
        for stick in 0..classes - 1 {
            let a = alpha[stick];
            let b: f64 = alpha[stick + 1..].iter().sum();
            let u = van_der_corput(point + 1, BASES[stick.min(BASES.len() - 1)]);
            let share = beta_quantile(u, a, b);
            piece.push(remaining * share);
            remaining *= 1.0 - share;
        }
        piece.push(remaining);
        frequencies.push(piece);
    }
    Support {
        ln_weights: vec![-(QMC_POINTS as f64).ln(); frequencies.len()],
        frequencies,
    }
}

/// The corners and the edges: a locus fixed for one length, or segregating two.
///
/// **Its own model, not a discretised Dirichlet.** A face carries no Dirichlet density, so
/// the weights are stated: a corner's is the monomorphic share times the spectrum, an edge's
/// is the rest times the two lengths' spectrum weights, and along an edge the minor
/// frequency is a symmetric Beta on a quantile quadrature. This is the same shape §2.1.2
/// fits on the ordinary-position path, with `classes` corners instead of two.
fn faces_support(parameters: &Parameters, classes: usize) -> Support {
    let mut frequencies = Vec::new();
    let mut ln_weights = Vec::new();

    for class in 0..classes {
        let mut point = vec![0.0; classes];
        point[class] = 1.0;
        frequencies.push(point);
        ln_weights
            .push((parameters.monomorphic_share * parameters.spectrum[class]).max(1e-300).ln());
    }

    let pair_total: f64 = (0..classes)
        .flat_map(|first| (first + 1..classes).map(move |second| (first, second)))
        .map(|(first, second)| parameters.spectrum[first] * parameters.spectrum[second])
        .sum();
    let nodes = beta_nodes(parameters.edge_shape, parameters.edge_shape, EDGE_NODES);
    for first in 0..classes {
        for second in first + 1..classes {
            let pair = parameters.spectrum[first] * parameters.spectrum[second]
                / pair_total.max(1e-300);
            for node in &nodes {
                let mut point = vec![0.0; classes];
                point[first] = *node;
                point[second] = 1.0 - node;
                frequencies.push(point);
                ln_weights.push(
                    ((1.0 - parameters.monomorphic_share) * pair / EDGE_NODES as f64)
                        .max(1e-300)
                        .ln(),
                );
            }
        }
    }
    Support {
        frequencies,
        ln_weights,
    }
}

/// One locus's read likelihoods, in a form the support can be swept over cheaply.
///
/// **Rescaled out of log space, once.** The inner loop of every candidate is *for each point
/// of the support, for each sample, sum the genotypes* — so doing it in logs costs one `exp`
/// per genotype per sample **per point**, which at the Dirichlet arm's hundreds of points is
/// where the whole program's time goes. Subtracting each sample's largest log-likelihood
/// once makes the sum a plain dot product, and the offsets are added back at the end.
struct LocusLikelihoods {
    /// `scaled[sample][genotype]`, in the genotype order [`genotype_pairs`] gives.
    scaled: Vec<Vec<f64>>,
    /// Σ over samples of the log-likelihood each row was divided by.
    ln_offset: f64,
}

fn genotype_pairs(classes: usize) -> Vec<(usize, usize)> {
    (0..classes)
        .flat_map(|first| (first..classes).map(move |second| (first, second)))
        .collect()
}

impl LocusLikelihoods {
    fn of(panel: &[ReadCounts], per_allele: &[Vec<f64>], genotypes: &[(usize, usize)]) -> Self {
        let mut scaled = Vec::with_capacity(panel.len());
        let mut ln_offset = 0.0;
        for counts in panel {
            let logs: Vec<f64> = genotypes
                .iter()
                .map(|(first, second)| ln_reads_given_genotype(counts, per_allele, *first, *second))
                .collect();
            let largest = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            ln_offset += largest;
            scaled.push(logs.into_iter().map(|log| (log - largest).exp()).collect());
        }
        Self { scaled, ln_offset }
    }
}

/// `ln P(this locus's panel | parameters)`.
fn ln_locus(
    likelihoods: &LocusLikelihoods,
    support: &Support,
    genotypes: &[(usize, usize)],
    inbreeding: f64,
) -> f64 {
    let mut terms = Vec::with_capacity(support.frequencies.len());
    let mut priors = vec![0.0_f64; genotypes.len()];
    for (point, ln_weight) in support.frequencies.iter().zip(&support.ln_weights) {
        for (slot, (first, second)) in genotypes.iter().enumerate() {
            priors[slot] = if first == second {
                inbreeding * point[*first] + (1.0 - inbreeding) * point[*first] * point[*first]
            } else {
                (1.0 - inbreeding) * 2.0 * point[*first] * point[*second]
            };
        }
        let mut total = *ln_weight;
        for row in &likelihoods.scaled {
            let mut sum = 0.0;
            for (prior, scaled) in priors.iter().zip(row) {
                sum += prior * scaled;
            }
            if sum <= 0.0 {
                total = f64::NEG_INFINITY;
                break;
            }
            total += sum.ln();
        }
        terms.push(total);
    }
    ln_sum_exp(&terms) + likelihoods.ln_offset
}

/// The data, with the read likelihoods already computed for one slippage value.
///
/// **Held between objective evaluations**, because the climb moves the spectrum, the
/// concentration and the monomorphic share far more often than it moves slippage, and none
/// of those changes a read likelihood.
struct Prepared {
    slippage: Slippage,
    loci: Vec<LocusLikelihoods>,
}

impl Prepared {
    fn for_slippage(data: &[Vec<ReadCounts>], slippage: Slippage, classes: usize) -> Self {
        let genotypes = genotype_pairs(classes);
        let per_allele: Vec<Vec<f64>> = (0..classes)
            .map(|allele| slippage.read_probabilities(allele, classes))
            .collect();
        Self {
            slippage,
            loci: data
                .iter()
                .map(|panel| LocusLikelihoods::of(panel, &per_allele, &genotypes))
                .collect(),
        }
    }

    fn refresh(&mut self, data: &[Vec<ReadCounts>], slippage: Slippage, classes: usize) {
        if self.slippage != slippage {
            *self = Self::for_slippage(data, slippage, classes);
        }
    }
}

fn objective(
    prepared: &Prepared,
    candidate: Candidate,
    parameters: &Parameters,
    inbreeding: f64,
    classes: usize,
) -> f64 {
    let support = support(candidate, parameters, classes);
    let genotypes = genotype_pairs(classes);
    prepared
        .loci
        .iter()
        .map(|locus| ln_locus(locus, &support, &genotypes, inbreeding))
        .sum::<f64>()
        / prepared.loci.len() as f64
}

// ---------------------------------------------------------------------
// Climbing
// ---------------------------------------------------------------------

fn logit(p: f64) -> f64 {
    let p = p.clamp(1e-9, 1.0 - 1e-9);
    (p / (1.0 - p)).ln()
}

fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// The data plus the read likelihoods for whichever slippage was last asked about.
///
/// **A scorer rather than a free function** so the cache has somewhere to live: the climb
/// asks about hundreds of spectra for every slippage value it tries, and recomputing a read
/// likelihood for each of them is the difference between seconds and hours.
struct Scorer<'a> {
    data: &'a [Vec<ReadCounts>],
    prepared: Prepared,
    classes: usize,
    inbreeding: f64,
    evaluations: u64,
}

impl<'a> Scorer<'a> {
    fn new(data: &'a [Vec<ReadCounts>], slippage: Slippage, classes: usize, inbreeding: f64) -> Self {
        Self {
            data,
            prepared: Prepared::for_slippage(data, slippage, classes),
            classes,
            inbreeding,
            evaluations: 0,
        }
    }

    fn score(&mut self, candidate: Candidate, parameters: &Parameters) -> f64 {
        self.prepared
            .refresh(self.data, parameters.slippage, self.classes);
        self.evaluations += 1;
        objective(
            &self.prepared,
            candidate,
            parameters,
            self.inbreeding,
            self.classes,
        )
    }
}

/// Golden-section on one coordinate, over a bracket of `span` either side of `start`.
fn climb_scalar(mut score: impl FnMut(f64) -> f64, start: f64, span: f64) -> f64 {
    const PHI: f64 = 0.618_033_988_749_895;
    let (mut low, mut high) = (start - span, start + span);
    let (mut c, mut d) = (high - PHI * (high - low), low + PHI * (high - low));
    let (mut fc, mut fd) = (score(c), score(d));
    for _ in 0..16 {
        if fc > fd {
            high = d;
            d = c;
            fd = fc;
            c = high - PHI * (high - low);
            fc = score(c);
        } else {
            low = c;
            c = d;
            fc = fd;
            d = low + PHI * (high - low);
            fd = score(d);
        }
    }
    if fc > fd { c } else { d }
}

/// Coordinate ascent over everything the candidate fits.
///
/// **Several starting points, best score taken** — the ordinary-position path measured that a
/// start putting its two noise classes close together collapses them and reports convergence,
/// costing 46% of the clean error rate, and nothing says this path is better behaved.
fn fit(
    data: &[Vec<ReadCounts>],
    candidate: Candidate,
    inbreeding: f64,
    classes: usize,
    starts: &[Parameters],
) -> (Parameters, f64) {
    let mut best: Option<(Parameters, f64)> = None;
    let mut scorer = Scorer::new(data, starts[0].slippage, classes, inbreeding);
    for start in starts {
        let mut parameters = start.clone();
        let mut score = scorer.score(candidate, &parameters);
        for _ in 0..5 {
            let before = score;

            // The three slippage numbers, on their own scales.
            for which in 0..3 {
                let current = match which {
                    0 => parameters.slippage.level,
                    1 => parameters.slippage.down_share,
                    _ => parameters.slippage.fall_off,
                };
                let moved = climb_scalar(
                    |x| {
                        let mut trial = parameters.clone();
                        let value = expit(x);
                        match which {
                            0 => trial.slippage.level = value,
                            1 => trial.slippage.down_share = value,
                            _ => trial.slippage.fall_off = value,
                        }
                        scorer.score(candidate, &trial)
                    },
                    logit(current),
                    3.0,
                );
                let value = expit(moved);
                match which {
                    0 => parameters.slippage.level = value,
                    1 => parameters.slippage.down_share = value,
                    _ => parameters.slippage.fall_off = value,
                }
            }

            // The spectrum, one class at a time on a log scale, renormalised each time.
            for class in 0..classes {
                let current = parameters.spectrum[class].max(1e-9).ln();
                let moved = climb_scalar(
                    |x| {
                        let mut trial = parameters.clone();
                        trial.spectrum[class] = x.exp();
                        let total: f64 = trial.spectrum.iter().sum();
                        for weight in &mut trial.spectrum {
                            *weight /= total;
                        }
                        scorer.score(candidate, &trial)
                    },
                    current,
                    2.0,
                );
                parameters.spectrum[class] = moved.exp();
                let total: f64 = parameters.spectrum.iter().sum();
                for weight in &mut parameters.spectrum {
                    *weight /= total;
                }
            }

            match candidate {
                Candidate::Dirichlet | Candidate::DirichletPoints => {
                    let moved = climb_scalar(
                        |x| {
                            let mut trial = parameters.clone();
                            trial.kappa = x.exp();
                            scorer.score(candidate, &trial)
                        },
                        parameters.kappa.ln(),
                        2.5,
                    );
                    parameters.kappa = moved.exp();
                }
                Candidate::Faces => {
                    let moved = climb_scalar(
                        |x| {
                            let mut trial = parameters.clone();
                            trial.monomorphic_share = expit(x);
                            scorer.score(candidate, &trial)
                        },
                        logit(parameters.monomorphic_share),
                        3.0,
                    );
                    parameters.monomorphic_share = expit(moved);
                    let moved = climb_scalar(
                        |x| {
                            let mut trial = parameters.clone();
                            trial.edge_shape = x.exp();
                            scorer.score(candidate, &trial)
                        },
                        parameters.edge_shape.ln(),
                        2.0,
                    );
                    parameters.edge_shape = moved.exp();
                }
                Candidate::PerStratum => {}
            }

            score = scorer.score(candidate, &parameters);
            if score - before < 1e-6 {
                break;
            }
        }
        if best.as_ref().is_none_or(|(_, b)| score > *b) {
            best = Some((parameters, score));
        }
    }
    best.expect("at least one start")
}

fn starting_points(classes: usize) -> Vec<Parameters> {
    // Spread over the slippage level and over how monomorphic loci are assumed to be, since
    // those are the two the objective is least likely to be unimodal in.
    vec![
        Parameters::start(classes, 0.02, 0.3, 0.9),
        Parameters::start(classes, 0.10, 3.0, 0.5),
        Parameters::start(classes, 0.30, 30.0, 0.1),
    ]
}

// ---------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------

fn relative_error(fitted: f64, truth: f64) -> String {
    if truth == 0.0 {
        return format!("{fitted:>10.4}");
    }
    format!("{:>+9.1}%", 100.0 * (fitted - truth) / truth)
}

fn report(truth: &Truth, candidate: Candidate, fitted: &Parameters, score: f64, seconds: f64) {
    println!(
        "  {:<12}  level {} ({:.4})   shorter-share {} ({:.3})   fall-off {} ({:.3})",
        candidate.name(),
        relative_error(fitted.slippage.level, truth.slippage.level),
        fitted.slippage.level,
        relative_error(fitted.slippage.down_share, truth.slippage.down_share),
        fitted.slippage.down_share,
        relative_error(fitted.slippage.fall_off, truth.slippage.fall_off),
        fitted.slippage.fall_off,
    );
    let spectrum: Vec<String> = fitted.spectrum.iter().map(|q| format!("{q:.3}")).collect();
    let shape = match candidate {
        Candidate::Dirichlet | Candidate::DirichletPoints => format!(
            "kappa {} ({:.3}, truth {:.3})",
            relative_error(fitted.kappa, truth.kappa),
            fitted.kappa,
            truth.kappa
        ),
        Candidate::Faces => format!(
            "monomorphic share {:.3}, edge shape {:.3}",
            fitted.monomorphic_share, fitted.edge_shape
        ),
        Candidate::PerStratum => "no shape number".to_string(),
    };
    println!(
        "                spectrum [{}]   {shape}",
        spectrum.join(", ")
    );
    println!("                objective {score:.6} nats a locus, {seconds:.1} s");
}

/// How often a locus drawn under this truth carries three or more lengths — the population
/// §4.2's support cannot represent.
fn interior_share(truth: &Truth, loci: usize, seed: u64, panel_chromosomes: usize) -> (f64, f64) {
    let mut rng = Rng(seed);
    let (mut monomorphic, mut over_two) = (0_u64, 0_u64);
    for _ in 0..loci {
        let frequencies = truth.draw_frequencies(&mut rng);
        let mut seen = std::collections::BTreeSet::new();
        for _ in 0..panel_chromosomes {
            seen.insert(rng.categorical(&frequencies));
        }
        if seen.len() == 1 {
            monomorphic += 1;
        }
        if seen.len() > 2 {
            over_two += 1;
        }
    }
    (
        monomorphic as f64 / loci as f64,
        over_two as f64 / loci as f64,
    )
}

fn run_fit(truth: &Truth, loci: usize, seed: u64, candidates: &[Candidate]) {
    println!(
        "\ntruth        {:?} frequencies, {} classes, {} samples, {} reads a locus, kappa {:.3}, \
         F {:.2}",
        truth.model, truth.classes, truth.samples, truth.depth, truth.kappa, truth.inbreeding
    );
    let spectrum: Vec<String> = truth.spectrum.iter().map(|q| format!("{q:.3}")).collect();
    println!("             spectrum [{}]", spectrum.join(", "));
    println!(
        "             slippage level {:.4}, shorter-share {:.3}, fall-off {:.3}",
        truth.slippage.level, truth.slippage.down_share, truth.slippage.fall_off
    );

    let (monomorphic, over_two) =
        interior_share(truth, 20_000, seed ^ 0x5EED, 2 * truth.samples);
    println!(
        "             of loci, {:.1}% carry one length in the panel and {:.1}% carry three or more",
        100.0 * monomorphic,
        100.0 * over_two
    );
    println!(
        "             — the second is the population the corners-and-edges support cannot hold"
    );

    let started = Instant::now();
    let data = truth.draw(loci, seed);
    println!(
        "\ndrawn        {loci} loci in {:.1} s\n",
        started.elapsed().as_secs_f64()
    );

    let starts = starting_points(truth.classes);
    for candidate in candidates {
        let started = Instant::now();
        let (fitted, score) = fit(&data, *candidate, truth.inbreeding, truth.classes, &starts);
        report(
            truth,
            *candidate,
            &fitted,
            score,
            started.elapsed().as_secs_f64(),
        );
        let points = support(*candidate, &fitted, truth.classes).frequencies.len();
        println!("                {points} points of support a locus");
    }

    // The objective at the truth, which separates a biased estimator from a climb that
    // stopped early — the diagnostic the ordinary-position harness found load-bearing.
    let at_truth = Parameters {
        spectrum: truth.spectrum.clone(),
        kappa: truth.kappa,
        monomorphic_share: monomorphic,
        edge_shape: 0.5,
        slippage: truth.slippage,
    };
    let mut scorer = Scorer::new(&data, truth.slippage, truth.classes, truth.inbreeding);
    for candidate in candidates {
        println!(
            "  {:<12}  objective at the true values {:.6}",
            candidate.name(),
            scorer.score(*candidate, &at_truth)
        );
    }
}

fn kappa_sweep(classes: usize, samples: usize, depth: u32, loci: usize) {
    println!(
        "\nHow monomorphic loci are, swept — {classes} classes, {samples} samples, {depth} reads \
         a locus, {loci} loci"
    );
    println!(
        "\n kappa   one length   3+ lengths | slippage level, error against truth 0.0800"
    );
    println!("                                  | per-stratum   dirichlet-pts       faces");
    for kappa in [0.05_f64, 0.2, 1.0, 5.0, 50.0] {
        let truth = Truth::default_for(classes, samples, depth, kappa);
        let (monomorphic, over_two) = interior_share(&truth, 20_000, 7, 2 * samples);
        let data = truth.draw(loci, 11);
        let starts = starting_points(classes);
        let mut cells = Vec::new();
        for candidate in [
            Candidate::PerStratum,
            Candidate::DirichletPoints,
            Candidate::Faces,
        ] {
            let (fitted, _) = fit(&data, candidate, truth.inbreeding, classes, &starts);
            cells.push(format!(
                "{:>7} ({:.4})",
                relative_error(fitted.slippage.level, truth.slippage.level),
                fitted.slippage.level
            ));
        }
        println!(
            "{kappa:>6.2}  {:>9.1}%  {:>10.1}% | {}",
            100.0 * monomorphic,
            100.0 * over_two,
            cells.join("  ")
        );
    }
}

fn seed_scatter(truth: &Truth, loci: usize, candidate: Candidate) {
    println!(
        "\nThe same truth from five draws — {} under {}, {loci} loci",
        truth.samples,
        candidate.name()
    );
    println!("  seed      level    shorter-share   fall-off");
    for seed in 1..=5_u64 {
        let data = truth.draw(loci, seed * 977);
        let (fitted, _) = fit(&data, candidate, truth.inbreeding, truth.classes, &starting_points(truth.classes));
        println!(
            "  {seed:>4}   {:>8.4}      {:>8.3}   {:>8.3}",
            fitted.slippage.level, fitted.slippage.down_share, fitted.slippage.fall_off
        );
    }
    println!(
        "  truth  {:>8.4}      {:>8.3}   {:>8.3}",
        truth.slippage.level, truth.slippage.down_share, truth.slippage.fall_off
    );
}

fn main() {
    let mut args = std::env::args().skip(1);
    let mode = args.next().unwrap_or_else(|| "fit".to_string());
    let samples: usize = args.next().map_or(20, |a| a.parse().expect("a sample count"));
    let depth: u32 = args.next().map_or(6, |a| a.parse().expect("a read depth"));
    let loci: usize = args.next().map_or(4_000, |a| a.parse().expect("a locus count"));
    let classes: usize = args.next().map_or(3, |a| a.parse().expect("a class count"));
    let kappa: f64 = args.next().map_or(0.5, |a| a.parse().expect("a concentration"));

    match mode.as_str() {
        "kappa" => kappa_sweep(classes, samples, depth, loci),
        "classes" => {
            for classes in [3_usize, 5] {
                let truth = Truth::default_for(classes, samples, depth, kappa);
                run_fit(
                    &truth,
                    loci,
                    11,
                    &[
                        Candidate::PerStratum,
                        Candidate::DirichletPoints,
                        Candidate::Faces,
                    ],
                );
            }
        }
        "seeds" => {
            let truth = Truth::default_for(classes, samples, depth, kappa);
            seed_scatter(&truth, loci, Candidate::Faces);
            seed_scatter(&truth, loci, Candidate::PerStratum);
        }
        "faces-truth" => {
            let mut truth = Truth::default_for(classes, samples, depth, kappa);
            truth.model = TruthModel::Faces;
            run_fit(
                &truth,
                loci,
                11,
                &[
                    Candidate::PerStratum,
                    Candidate::DirichletPoints,
                    Candidate::Faces,
                ],
            );
        }
        _ => {
            let truth = Truth::default_for(classes, samples, depth, kappa);
            run_fit(
                &truth,
                loci,
                11,
                &[
                    Candidate::PerStratum,
                    Candidate::Dirichlet,
                    Candidate::DirichletPoints,
                    Candidate::Faces,
                ],
            );
        }
    }
}
