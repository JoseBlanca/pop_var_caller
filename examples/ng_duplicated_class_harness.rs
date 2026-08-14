//! **What does a duplicated position cost the joint fit's three headline numbers, and can the
//! cohort recognise one without a per-sample coverage summary?**
//!
//! A plant that carries two copies of a stretch the reference holds once has both copies' reads
//! landing at the same place. Wherever the two copies differ, about half the reads disagree with the
//! reference — which is what a heterozygote looks like. The joint parameters fit
//! (`spec/parameter_prepass_joint_fit.md` §2.2) answers that with a **third class of site**, whose
//! discriminator is the sample's local relative coverage: a window carrying twice the sample's
//! normal depth is carrying two copies. That discriminator is the whole reason for a per-sample
//! coverage-by-window summary, which in a two-phase run costs a full extra pass over every sample's
//! pileup (`spec/parameter_prepass_joint_records.md` §4).
//!
//! This program fits **one drawn cohort three ways** and reports what each way gets:
//!
//! - **`coverage`** — the third class, recognised from each sample's local relative coverage;
//! - **`pattern`** — the same third class, recognised only from the pattern across the cohort:
//!   a real variant at frequency a half leaves about a quarter of the samples homozygous for the
//!   non-reference allele, and a duplication leaves nobody there;
//! - **`no-class`** — no third class at all, which is what the fit does today.
//!
//! Everything else is identical: same drawn panel, same starting points, same search. The three
//! differ only in what evidence the third class's per-sample carrier state is allowed to read, and
//! `no-class` holds that class's weight at zero.
//!
//! ## The two questions it answers
//!
//! **Does ignoring the class cancel out of the inbreeding coefficient?** A duplicated position
//! raises observed heterozygosity, because the fit's only home for it is *heterozygous*. But it also
//! enters the fitted allele-frequency density as a mid-frequency variant, and expected heterozygosity
//! is read off that density (spec §5.3). The inbreeding coefficient is `1 − Hobs/Hexp`, so the two
//! inflations might cancel there. Each of the three numbers is reported separately against the drawn
//! truth, which is the only way to see whether they do.
//!
//! **Can the coverage summary be dropped?** That is `pattern` against `coverage`, plus the two counts
//! only a drawn truth can supply: how many drawn duplicated loci `pattern` finds, and how many drawn
//! variants it wrongly calls duplicated.
//!
//! ## How the panel is drawn
//!
//! Three kinds of locus, from `reports/duplicated_locus_probe_2026-08-12.md`:
//!
//! - **ordinary loci** — a population frequency from a monomorphic mass, a fixed-alternative mass and
//!   a Beta, then a genotype per sample under Hardy–Weinberg with inbreeding (the truth
//!   `ng_joint_fit_harness.rs` already uses);
//! - **duplicated loci** — a carrier frequency per locus, then a carrier indicator per sample. A
//!   carrier gets **twice the depth** and about **half its reads disagreeing**; a non-carrier is
//!   homozygous reference at normal depth;
//! - and every sample at every locus carries a **relative coverage** reading, drawn around its own
//!   copy number with the scatter a window of that many aligned bases really has.
//!
//! Run it:
//!
//! ```text
//! ng_duplicated_class_harness [samples] [depth] [loci] [inbreeding] [window_aligned_bases] \
//!                             [collapse_share] [carrier_a] [carrier_b]
//! ```
//!
//! A depth of `0` runs both three reads a position and twenty-five. The last three set how the
//! carriers are spread — the share of duplicated positions every sample carries, then the two shapes
//! of the Beta the rest draw their carrier frequency from. **That spread is what the answer turns
//! on**, so it is a knob rather than a constant.
//!
//! Diploid throughout, as the genotype prior is (spec §10).

use rayon::prelude::*;
use std::time::Instant;

// ---------------------------------------------------------------------
// Small numerics — the same ones `ng_joint_fit_harness.rs` uses.
// ---------------------------------------------------------------------

fn ln_gamma(x: f64) -> f64 {
    // Lanczos, g = 7, n = 9 — the coefficients every implementation uses.
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

fn ln_add(a: f64, b: f64) -> f64 {
    if a > b {
        a + (b - a).exp().ln_1p()
    } else if b.is_finite() {
        b + (a - b).exp().ln_1p()
    } else {
        a
    }
}

/// The regularised incomplete beta `I_x(a, b)` — the Beta's own cumulative, which turns a grid of
/// frequency cells into exact cell masses. Continued fraction, the standard Lentz evaluation.
fn regularised_incomplete_beta(x: f64, a: f64, b: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
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

/// A deterministic generator, so a run reproduces without a seed file. splitmix64.
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
    fn poisson(&mut self, mean: f64) -> u32 {
        // Knuth for the small means, a normal approximation above 30 where Knuth's product
        // underflows and the approximation is already good to a fraction of a read.
        if mean > 30.0 {
            let draw = mean + mean.sqrt() * self.normal();
            return draw.round().max(0.0) as u32;
        }
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
    fn from_weights(&mut self, weights: &[f64]) -> usize {
        let mut u = self.uniform() * weights.iter().sum::<f64>();
        for (i, w) in weights.iter().enumerate() {
            u -= w;
            if u <= 0.0 {
                return i;
            }
        }
        weights.len() - 1
    }
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

/// The chance one read shows something other than the reference base, for an individual carrying
/// `j` alternative copies of `PLOIDY` (`spec/parameter_prepass.md` §3).
///
/// **A duplication carrier reads at `p_alt(1, eps)` too** — half its reads come from the copy the
/// reference does not hold. That is exactly why the read counts alone cannot separate the two.
fn p_alt(j: u32, eps: f64) -> f64 {
    let share = f64::from(j) / f64::from(PLOIDY);
    share * (1.0 - eps / 3.0) + (1.0 - share) * eps
}

/// Diploid genotype frequencies at population frequency `f` with inbreeding `inbreed` —
/// Hardy–Weinberg with the inbreeding correction.
fn genotype_frequencies(f: f64, inbreed: f64) -> [f64; 3] {
    let hetero = 2.0 * f * (1.0 - f) * (1.0 - inbreed);
    let shared = inbreed * f * (1.0 - f);
    [(1.0 - f) * (1.0 - f) + shared, hetero, f * f + shared]
}

/// `∫ Beta(a, b) · 2f(1−f) df`, in closed form.
fn beta_heterozygosity(a: f64, b: f64) -> f64 {
    2.0 * a * b / ((a + b) * (a + b + 1.0))
}

/// A grid of frequency cells with fixed edges, arcsine-spaced so cells cluster at both ends —
/// where a Beta with shapes below one puts its mass, and where a rare allele lives. Each cell's
/// mass and conditional mean come from the Beta's own cumulative, so this is quadrature rather
/// than a discretisation of the parameter.
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

    /// Each cell's `(mass, representative frequency)` under `Beta(a, b)`.
    fn under_beta(&self, a: f64, b: f64) -> Vec<(f64, f64)> {
        let mean = a / (a + b);
        let mut lower = 0.0;
        let mut lower_shifted = 0.0;
        let mut out = Vec::with_capacity(self.edges.len() - 1);
        for e in 1..self.edges.len() {
            let upper = regularised_incomplete_beta(self.edges[e], a, b);
            let upper_shifted = regularised_incomplete_beta(self.edges[e], a + 1.0, b);
            let mass = (upper - lower).max(0.0);
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
// How well a window's mean depth tells one copy from two
// ---------------------------------------------------------------------

/// **What a window's relative coverage reading looks like when the sample really carries `copy`
/// copies there.** Relative coverage is the window's mean depth over the sample's own
/// depth-against-GC curve, rescaled so the median window sits at 1.0 — so a two-copy window reads
/// about 2.0, blurred by two kinds of scatter.
///
/// Both are calibrated against
/// [`reports/duplicated_locus_probe_2026-08-12.md`](../doc/devel/ng/reports/duplicated_locus_probe_2026-08-12.md):
///
/// - **read-count scatter**, which shrinks as the window collects more aligned bases. At tomato
///   SRR7279482's 25× over 500 bp — 12,600 aligned bases — 86% of windows land inside 0.6 to 1.4 and
///   0.86% reach 1.6 or above, which puts the standard deviation at about 0.25;
/// - **a window-to-window floor** that no amount of depth removes, from mappability and from the
///   GC curve being fitted rather than known. At 3.6× over 500 bp — 1,800 aligned bases — 9.7% of
///   positions land in the two-copy band, which puts the standard deviation there at about 0.46.
///
/// Those two points fix `sd² = copy²·floor² + copy·scatter/aligned_bases` at a floor of 0.194 and a
/// scatter of 313 aligned bases. The scatter is twice what independent 150 bp reads would give,
/// which is the usual amount of clustering.
#[derive(Clone, Copy)]
struct CoverageReading {
    aligned_bases_per_window: f64,
}

impl CoverageReading {
    const FLOOR: f64 = 0.194;
    const SCATTER: f64 = 313.0;

    fn standard_deviation(&self, copy: f64) -> f64 {
        (copy * copy * Self::FLOOR * Self::FLOOR
            + copy * Self::SCATTER / self.aligned_bases_per_window)
            .sqrt()
    }

    fn draw(&self, rng: &mut Rng, copy: f64) -> f64 {
        (copy + self.standard_deviation(copy) * rng.normal()).max(0.02)
    }

    /// `ln P(reading | two copies) − ln P(reading | one copy)` — the only way the coverage summary
    /// enters the likelihood, and the one line that separates the `coverage` fit from the
    /// `pattern` fit.
    fn ln_two_copies_over_one(&self, reading: f64) -> f64 {
        let ln_normal = |mean: f64| {
            let sd = self.standard_deviation(mean);
            -0.5 * ((reading - mean) / sd).powi(2) - sd.ln()
        };
        ln_normal(2.0) - ln_normal(1.0)
    }
}

// ---------------------------------------------------------------------
// The truth
// ---------------------------------------------------------------------

/// **How many duplicated carrier positions a sample has, per genuinely heterozygous position.**
///
/// On tomato SRR7279482 at 25×, positions that are both in a window near two copies and reading
/// between 35% and 65% alternative number 150 to 590 per two million across eight samples, against
/// about 668 genuinely near-half positions per two million in the same sample — so the artefact is
/// about a third of the real signal. **The ratio is what the bias depends on**, not either count on
/// its own, so it is what the drawn panel holds fixed.
const DUPLICATED_PER_HETEROZYGOUS: f64 = 1.0 / 3.0;

/// **The share of duplicated loci that every sample carries** — the reference's own collapse rather
/// than copy number segregating in the panel. Of 84 windows read near two copies by at least one of
/// eight tomato samples, 11 are read that way by seven or eight of them; allowing for the windows a
/// segregating duplication leaves with no carrier at all among eight, that is 9.0% of duplicated
/// loci (probe report §5).
const REFERENCE_COLLAPSE_SHARE: f64 = 0.090;

/// **The carrier frequency of the rest**, as `Beta(a, b)`. Fitted to the same eight-sample counts:
/// of the 73 windows that are not the reference's own collapse, 40 are read near two copies by
/// exactly one sample and 33 by two to six, which `Beta(1.19, 9.55)` reproduces to within a
/// thousandth. Its mean carrier frequency is 0.111.
const CARRIER_BETA: (f64, f64) = (1.1895, 9.5477);

#[derive(Clone)]
struct Truth {
    samples: usize,
    depth_mean: f64,
    eps_clean: f64,
    eps_noisy: f64,
    w_noisy: f64,
    p_invariant: f64,
    p_fixed_alt: f64,
    beta_a: f64,
    beta_b: f64,
    inbreeding: f64,
    /// The share of loci that are duplicated in at least one sample's genome — derived, so that the
    /// drawn panel holds [`DUPLICATED_PER_HETEROZYGOUS`].
    duplicated_share: f64,
    /// How the carriers are spread: a share carried by every sample, and a `Beta` over the carrier
    /// frequency of the rest. **This is the input the cohort-pattern answer is most sensitive to**,
    /// so it is a run parameter rather than a constant.
    collapse_share: f64,
    carrier_a: f64,
    carrier_b: f64,
    coverage: CoverageReading,
}

impl Truth {
    /// The population's expected heterozygosity over the **ordinary** loci — the quantity `Hexp`
    /// estimates. A duplicated locus is not population polymorphism and contributes nothing to it.
    fn expected_heterozygosity(&self) -> f64 {
        self.p_segregating() * beta_heterozygosity(self.beta_a, self.beta_b)
    }

    fn p_segregating(&self) -> f64 {
        1.0 - self.p_invariant - self.p_fixed_alt
    }

    /// The mean carrier frequency of a duplicated locus, over both components.
    fn mean_carrier_frequency(&self) -> f64 {
        self.collapse_share
            + (1.0 - self.collapse_share) * self.carrier_a / (self.carrier_a + self.carrier_b)
    }

    fn tomato_like(
        samples: usize,
        depth_mean: f64,
        inbreeding: f64,
        aligned_bases_per_window: f64,
        carriers: (f64, f64, f64),
    ) -> Self {
        let mut truth = Self {
            samples,
            depth_mean,
            // The rates the per-sample histogram route fits on HG002, the only place they have been
            // measured against a truth set (`research/noise_model_overdispersion_2026-08-10.md`).
            eps_clean: 1.895e-3,
            eps_noisy: 5.29e-2,
            w_noisy: 0.0088,
            p_invariant: 0.99,
            p_fixed_alt: 0.002,
            // `a < 1` is the rare-allele pile-up a neutral population has.
            beta_a: 0.3,
            beta_b: 1.2,
            inbreeding,
            duplicated_share: 0.0,
            collapse_share: carriers.0,
            carrier_a: carriers.1,
            carrier_b: carriers.2,
            coverage: CoverageReading {
                aligned_bases_per_window,
            },
        };
        // A sample is heterozygous at `Hexp · (1 − F)` of the loci; give it a third as many
        // duplicated carrier positions, which is what tomato has.
        let heterozygous_rate = truth.expected_heterozygosity() * (1.0 - inbreeding);
        truth.duplicated_share =
            heterozygous_rate * DUPLICATED_PER_HETEROZYGOUS / truth.mean_carrier_frequency();
        truth
    }
}

// ---------------------------------------------------------------------
// The drawn panel
// ---------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq)]
enum LocusKind {
    Monomorphic,
    FixedAlternative,
    Segregating,
    Duplicated,
}

/// One drawn cohort: every sample's reads and coverage reading at every locus, plus what each locus
/// really is.
struct Panel {
    samples: usize,
    loci: usize,
    /// `[locus * samples + sample]`.
    depth: Vec<u16>,
    alt: Vec<u16>,
    /// `ln P(reading | two copies) − ln P(reading | one copy)`, per sample per locus.
    coverage_ln_ratio: Vec<f32>,
    /// Per locus: how many silent samples sat at each depth, and the ones that spoke.
    quiet_at_depth: Vec<Vec<u32>>,
    loud: Vec<Vec<(u16, u16)>>,
    kind: Vec<LocusKind>,
    /// Per locus: how many samples carry the duplication (0 elsewhere).
    carriers: Vec<u32>,
    /// Per locus: how many samples are truly heterozygous (0 at a duplicated locus).
    heterozygotes: Vec<u32>,
    deepest: usize,
    /// The drawn panel's own expected heterozygosity, `mean over loci of 2f(1−f)`.
    drawn_expected_heterozygosity: f64,
}

impl Panel {
    fn draw(truth: &Truth, loci: usize, seed: u64) -> Self {
        let samples = truth.samples;
        let mut rng = Rng(seed);
        let mut depth = vec![0_u16; loci * samples];
        let mut alt = vec![0_u16; loci * samples];
        let mut coverage_ln_ratio = vec![0.0_f32; loci * samples];
        let mut quiet_at_depth = Vec::with_capacity(loci);
        let mut loud = Vec::with_capacity(loci);
        let mut kind = Vec::with_capacity(loci);
        let mut carriers = vec![0_u32; loci];
        let mut heterozygotes = vec![0_u32; loci];
        let mut deepest = 0_usize;
        let mut frequency_sum = 0.0;

        for locus in 0..loci {
            let eps = if rng.uniform() < truth.w_noisy {
                truth.eps_noisy
            } else {
                truth.eps_clean
            };
            let duplicated = rng.uniform() < truth.duplicated_share;
            let mut this_kind = LocusKind::Segregating;
            let mut frequency = 0.0;
            let mut carrier_frequency = 0.0;
            if duplicated {
                this_kind = LocusKind::Duplicated;
                carrier_frequency = if rng.uniform() < truth.collapse_share {
                    1.0
                } else {
                    rng.beta(truth.carrier_a, truth.carrier_b)
                };
            } else {
                let u = rng.uniform();
                if u < truth.p_invariant {
                    this_kind = LocusKind::Monomorphic;
                    frequency = 0.0;
                } else if u < truth.p_invariant + truth.p_fixed_alt {
                    this_kind = LocusKind::FixedAlternative;
                    frequency = 1.0;
                } else {
                    frequency = rng.beta(truth.beta_a, truth.beta_b);
                }
                frequency_sum += 2.0 * frequency * (1.0 - frequency);
            }

            let base = locus * samples;
            let mut quiet = vec![0_u32; 1];
            let mut loud_here = Vec::new();
            for sample in 0..samples {
                let (copies, read_alt_probability, heterozygous, carrier) = if duplicated {
                    if rng.uniform() < carrier_frequency {
                        // Two copies' reads at one position: twice the depth, and about half the
                        // reads disagreeing wherever the copies differ.
                        (2.0, p_alt(1, eps), false, true)
                    } else {
                        (1.0, p_alt(0, eps), false, false)
                    }
                } else {
                    let priors = genotype_frequencies(frequency, truth.inbreeding);
                    let genotype = rng.from_weights(&priors) as u32;
                    (1.0, p_alt(genotype, eps), genotype == 1, false)
                };
                if heterozygous {
                    heterozygotes[locus] += 1;
                }
                if carrier {
                    carriers[locus] += 1;
                }
                let d = rng.poisson(truth.depth_mean * copies);
                let k = rng.binomial(d, read_alt_probability);
                let reading = truth.coverage.draw(&mut rng, copies);
                depth[base + sample] = d as u16;
                alt[base + sample] = k as u16;
                coverage_ln_ratio[base + sample] =
                    truth.coverage.ln_two_copies_over_one(reading) as f32;
                if quiet.len() <= d as usize {
                    quiet.resize(d as usize + 1, 0);
                }
                if k == 0 {
                    quiet[d as usize] += 1;
                } else {
                    loud_here.push((d as u16, k as u16));
                }
                deepest = deepest.max(d as usize);
            }
            quiet_at_depth.push(quiet);
            loud.push(loud_here);
            kind.push(this_kind);
        }

        Self {
            samples,
            loci,
            depth,
            alt,
            coverage_ln_ratio,
            quiet_at_depth,
            loud,
            kind,
            carriers,
            heterozygotes,
            deepest,
            drawn_expected_heterozygosity: frequency_sum / loci as f64,
        }
    }

    /// The drawn panel's own observed heterozygosity: the share of (locus, sample) pairs at which
    /// the sample really is heterozygous. **A duplication carrier is not counted** — it is not
    /// heterozygous, and that is the whole point.
    fn drawn_observed_heterozygosity(&self) -> f64 {
        let total: u64 = self.heterozygotes.iter().map(|&h| u64::from(h)).sum();
        total as f64 / (self.loci as f64 * self.samples as f64)
    }

    /// The share of (locus, sample) pairs at which the sample carries a duplication — the quantity
    /// the tomato measurement counts at 150 to 590 per two million.
    fn drawn_carrier_rate(&self) -> f64 {
        let total: u64 = self.carriers.iter().map(|&c| u64::from(c)).sum();
        total as f64 / (self.loci as f64 * self.samples as f64)
    }
}

// ---------------------------------------------------------------------
// The fitted model
// ---------------------------------------------------------------------

/// Which evidence the third class's per-sample carrier state is allowed to read.
#[derive(Clone, Copy, PartialEq)]
enum Arm {
    /// The class is there and each sample's local relative coverage is evidence for it.
    Coverage,
    /// The class is there and the only evidence is the pattern of alternative reads across the
    /// cohort.
    Pattern,
    /// No third class at all — the fit as it stands today.
    NoClass,
}

impl Arm {
    fn name(self) -> &'static str {
        match self {
            Arm::Coverage => "coverage",
            Arm::Pattern => "pattern",
            Arm::NoClass => "no-class",
        }
    }
    fn has_third_class(self) -> bool {
        self != Arm::NoClass
    }
    fn reads_coverage(self) -> bool {
        self == Arm::Coverage
    }
}

#[derive(Clone)]
struct Fitted {
    eps_clean: f64,
    eps_noisy: f64,
    w_noisy: f64,
    p_invariant: f64,
    p_fixed_alt: f64,
    beta_a: f64,
    beta_b: f64,
    inbreeding: f64,
    /// The third class's weight, and the Beta its carrier frequency is drawn from.
    duplicated_share: f64,
    carrier_a: f64,
    carrier_b: f64,
}

impl Fitted {
    fn p_segregating(&self) -> f64 {
        (1.0 - self.p_invariant - self.p_fixed_alt).max(0.0)
    }

    fn expected_heterozygosity(&self) -> f64 {
        self.p_segregating() * beta_heterozygosity(self.beta_a, self.beta_b)
    }

    /// A deliberately wrong start: three times the truth's error rate and a flat frequency shape.
    /// **The three variants span the separation between the clean and the noisy class**, which
    /// `spec/parameter_prepass_joint_fit.md` §11 question 3 measured to be worth more than any other
    /// choice in this fit — a start that puts the two classes near each other collapses them into
    /// one and reports convergence.
    fn start(variant: usize, arm: Arm) -> Self {
        let (eps_clean, eps_noisy, w_noisy) = match variant {
            0 => (5e-4, 2e-1, 0.002),
            1 => (6e-3, 2e-2, 0.05),
            _ => (2e-2, 6e-2, 0.20),
        };
        Self {
            eps_clean,
            eps_noisy,
            w_noisy,
            p_invariant: 0.95,
            p_fixed_alt: 0.01,
            beta_a: 1.0,
            beta_b: 1.0,
            inbreeding: 0.1,
            duplicated_share: if arm.has_third_class() { 1e-3 } else { 0.0 },
            carrier_a: 1.0,
            carrier_b: 4.0,
        }
    }
}

/// `ln` of the binomial read likelihood for every `(genotype, error class, depth, alternative
/// reads)` the panel can show.
///
/// **Rebuilt whenever an error rate moves and looked up everywhere else.** Without it the third
/// class's per-sample product would call `ln_gamma` fifty times a locus per objective evaluation,
/// which is the whole cost of the fit.
struct BinomialTable {
    deepest: usize,
    /// `[(genotype * 2 + class) * (deepest+1) * (deepest+1) + depth * (deepest+1) + alt]`
    values: Vec<f64>,
}

impl BinomialTable {
    fn new(deepest: usize, eps: [f64; 2]) -> Self {
        let stride = deepest + 1;
        let mut values = vec![f64::NEG_INFINITY; 3 * 2 * stride * stride];
        for genotype in 0..3_u32 {
            for (class, &e) in eps.iter().enumerate() {
                let p = p_alt(genotype, e);
                let block = (genotype as usize * 2 + class) * stride * stride;
                for d in 0..=deepest {
                    for k in 0..=d {
                        values[block + d * stride + k] = ln_binomial(k as u32, d as u32, p);
                    }
                }
            }
        }
        Self { deepest, values }
    }

    #[inline]
    fn get(&self, genotype: usize, class: usize, depth: u16, alt: u16) -> f64 {
        let stride = self.deepest + 1;
        self.values
            [(genotype * 2 + class) * stride * stride + depth as usize * stride + alt as usize]
    }
}

/// Everything a silent sample can contribute to the ordinary-variant class, tabulated by depth so
/// that no locus has to evaluate it per sample. Only legal because every sample shares one
/// inbreeding coefficient here.
struct QuietTables {
    /// `[class][cell][depth]`
    segregating: Vec<Vec<Vec<f64>>>,
    /// `[class][depth]`
    all_reference: Vec<Vec<f64>>,
    all_alternative: Vec<Vec<f64>>,
    /// `[class][carrier cell][depth]` — the third class's silent samples, which are interchangeable
    /// only where no coverage reading distinguishes them.
    duplicated: Option<Vec<Vec<Vec<f64>>>>,
}

impl QuietTables {
    fn new(
        binomial: &BinomialTable,
        cells: &[(f64, f64)],
        carrier_cells: &[(f64, f64)],
        inbreeding: f64,
        arm: Arm,
    ) -> Self {
        let deepest = binomial.deepest as u16;
        let depths = || 0..=deepest;
        let segregating = (0..2)
            .map(|class| {
                cells
                    .iter()
                    .map(|&(_, f)| {
                        let priors = genotype_frequencies(f, inbreeding);
                        depths()
                            .map(|depth| {
                                ln_sum_exp(&[
                                    priors[0].max(1e-300).ln() + binomial.get(0, class, depth, 0),
                                    priors[1].max(1e-300).ln() + binomial.get(1, class, depth, 0),
                                    priors[2].max(1e-300).ln() + binomial.get(2, class, depth, 0),
                                ])
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect();
        let all_reference = (0..2)
            .map(|class| depths().map(|d| binomial.get(0, class, d, 0)).collect())
            .collect();
        let all_alternative = (0..2)
            .map(|class| depths().map(|d| binomial.get(2, class, d, 0)).collect())
            .collect();
        let duplicated = (arm.has_third_class() && !arm.reads_coverage()).then(|| {
            (0..2)
                .map(|class| {
                    carrier_cells
                        .iter()
                        .map(|&(_, q)| {
                            depths()
                                .map(|depth| {
                                    ln_add(
                                        (1.0 - q).max(1e-300).ln()
                                            + binomial.get(0, class, depth, 0),
                                        q.max(1e-300).ln() + binomial.get(1, class, depth, 0),
                                    )
                                })
                                .collect()
                        })
                        .collect()
                })
                .collect()
        });
        Self {
            segregating,
            all_reference,
            all_alternative,
            duplicated,
        }
    }
}

/// The whole frequency description a fit is currently holding, built once per objective evaluation.
struct View {
    cells: Vec<(f64, f64)>,
    ln_genotype_prior: Vec<[f64; 3]>,
    carrier_cells: Vec<(f64, f64)>,
}

fn build_view(fitted: &Fitted, variant_grid: &FrequencyGrid, carrier_grid: &FrequencyGrid) -> View {
    let cells = variant_grid.under_beta(fitted.beta_a, fitted.beta_b);
    let ln_genotype_prior = cells
        .iter()
        .map(|&(_, f)| {
            let p = genotype_frequencies(f, fitted.inbreeding);
            [
                p[0].max(1e-300).ln(),
                p[1].max(1e-300).ln(),
                p[2].max(1e-300).ln(),
            ]
        })
        .collect();
    let carrier_cells = carrier_grid.under_beta(fitted.carrier_a, fitted.carrier_b);
    View {
        cells,
        ln_genotype_prior,
        carrier_cells,
    }
}

/// Buffers one worker thread reuses across loci.
///
/// **Every one of them would otherwise be an allocation per locus per objective evaluation** — two
/// hundred thousand of them a pass — and the fit makes thousands of passes.
#[derive(Default)]
struct Scratch {
    /// The mixture components of one class at one locus, before they are summed.
    terms: Vec<f64>,
    /// Each sample's two branches under the third class, exponentiated once and reused at every
    /// carrier-frequency cell: `(homozygous reference, duplication carrier)`.
    branches: Vec<(f64, f64)>,
}

/// One locus's log-likelihood.
///
/// **The coverage readings factor almost all the way out.** Under every class except the third one,
/// every sample is single-copy, so its coverage reading contributes `P(reading | one copy)` — the
/// same factor at every class and independent of every parameter. Dividing it out leaves the
/// coverage entering in exactly one place: the ratio `P(reading | two copies) / P(reading | one
/// copy)` multiplying the carrier branch of the third class. So the `coverage` and `pattern` fits
/// differ by that ratio and by nothing else.
#[allow(clippy::too_many_arguments)]
fn ln_locus(
    panel: &Panel,
    locus: usize,
    binomial: &BinomialTable,
    fitted: &Fitted,
    view: &View,
    quiet: &QuietTables,
    arm: Arm,
    scratch: &mut Scratch,
) -> f64 {
    let class_weights = [1.0 - fitted.w_noisy, fitted.w_noisy];
    let mut per_class = [f64::NEG_INFINITY; 2];
    for (class, &class_weight) in class_weights.iter().enumerate() {
        if class_weight <= 0.0 {
            continue;
        }
        let ordinary =
            ln_ordinary_class(panel, locus, binomial, fitted, view, quiet, class, scratch);
        let total = if arm.has_third_class() && fitted.duplicated_share > 0.0 {
            let duplicated =
                ln_duplicated_class(panel, locus, binomial, view, quiet, class, scratch);
            ln_add(
                (1.0 - fitted.duplicated_share).max(1e-300).ln() + ordinary,
                fitted.duplicated_share.ln() + duplicated,
            )
        } else {
            ordinary
        };
        per_class[class] = class_weight.ln() + total;
    }
    ln_sum_exp(&per_class)
}

/// The two masses and the Beta — the ordinary-variant class, with the silent samples read off the
/// by-depth tables rather than one at a time.
#[allow(clippy::too_many_arguments)]
fn ln_ordinary_class(
    panel: &Panel,
    locus: usize,
    binomial: &BinomialTable,
    fitted: &Fitted,
    view: &View,
    quiet: &QuietTables,
    class: usize,
    scratch: &mut Scratch,
) -> f64 {
    let loud = &panel.loud[locus];
    let quiet_at_depth = &panel.quiet_at_depth[locus];
    let quiet_sum = |table: &[f64]| -> f64 {
        quiet_at_depth
            .iter()
            .enumerate()
            .filter(|&(_, &count)| count > 0)
            .map(|(depth, &count)| f64::from(count) * table[depth])
            .sum()
    };
    let terms = &mut scratch.terms;
    terms.clear();
    if fitted.p_invariant > 0.0 {
        let mut total = fitted.p_invariant.ln() + quiet_sum(&quiet.all_reference[class]);
        for &(d, k) in loud {
            total += binomial.get(0, class, d, k);
        }
        terms.push(total);
    }
    if fitted.p_fixed_alt > 0.0 {
        let mut total = fitted.p_fixed_alt.ln() + quiet_sum(&quiet.all_alternative[class]);
        for &(d, k) in loud {
            total += binomial.get(2, class, d, k);
        }
        terms.push(total);
    }
    let segregating = fitted.p_segregating();
    if segregating > 0.0 {
        for (cell, &(mass, _)) in view.cells.iter().enumerate() {
            if mass <= 0.0 {
                continue;
            }
            let mut total =
                segregating.ln() + mass.ln() + quiet_sum(&quiet.segregating[class][cell]);
            let prior = &view.ln_genotype_prior[cell];
            for &(d, k) in loud {
                total += ln_sum_exp(&[
                    prior[0] + binomial.get(0, class, d, k),
                    prior[1] + binomial.get(1, class, d, k),
                    prior[2] + binomial.get(2, class, d, k),
                ]);
            }
            terms.push(total);
        }
    }
    ln_sum_exp(terms)
}

/// The third class: a carrier frequency from its own Beta, then every sample either a carrier —
/// half its reads disagreeing — or homozygous reference.
///
/// **What makes it a different shape from an ordinary variant at the same frequency is what it has
/// no room for**: a sample homozygous for the non-reference allele. An ordinary variant at a
/// frequency of a half leaves about a quarter of the panel there; a duplication leaves none.
fn ln_duplicated_class(
    panel: &Panel,
    locus: usize,
    binomial: &BinomialTable,
    view: &View,
    quiet: &QuietTables,
    class: usize,
    scratch: &mut Scratch,
) -> f64 {
    // **Where a coverage reading is in play, each sample's two branches are exponentiated once and
    // reused at every carrier-frequency cell.** Written the obvious way the cell loop would take two
    // exponentials and a logarithm per sample per cell, and that loop is the whole cost of this fit.
    let mut common_scale = 0.0;
    if quiet.duplicated.is_none() {
        let base = locus * panel.samples;
        scratch.branches.clear();
        for sample in 0..panel.samples {
            let d = panel.depth[base + sample];
            let k = panel.alt[base + sample];
            let ln_reference = binomial.get(0, class, d, k);
            let ln_carrier =
                binomial.get(1, class, d, k) + f64::from(panel.coverage_ln_ratio[base + sample]);
            let scale = ln_reference.max(ln_carrier);
            common_scale += scale;
            scratch
                .branches
                .push(((ln_reference - scale).exp(), (ln_carrier - scale).exp()));
        }
    }
    let mut terms = std::mem::take(&mut scratch.terms);
    terms.clear();
    for (cell, &(mass, q)) in view.carrier_cells.iter().enumerate() {
        if mass <= 0.0 {
            continue;
        }
        let not_carrier = (1.0 - q).max(1e-300);
        let carrier = q.max(1e-300);
        let mut total = mass.ln();
        match &quiet.duplicated {
            // No coverage reading: silent samples at one depth are interchangeable, so only the
            // samples that showed an alternative read are scored one at a time.
            Some(table) => {
                for (depth, &count) in panel.quiet_at_depth[locus].iter().enumerate() {
                    if count > 0 {
                        total += f64::from(count) * table[class][cell][depth];
                    }
                }
                for &(d, k) in &panel.loud[locus] {
                    total += ln_add(
                        not_carrier.ln() + binomial.get(0, class, d, k),
                        carrier.ln() + binomial.get(1, class, d, k),
                    );
                }
            }
            // With one: every sample carries its own evidence, so every sample is scored.
            None => {
                total += common_scale;
                for &(reference, carrier_branch) in &scratch.branches {
                    total += (not_carrier * reference + carrier * carrier_branch).ln();
                }
            }
        }
        terms.push(total);
    }
    let out = ln_sum_exp(&terms);
    scratch.terms = terms;
    out
}

/// The objective: the mean log-likelihood over the drawn loci.
///
/// Collected and then summed in order, never reduced in parallel — a tree reduction over floats
/// depends on how the work was split, so the objective would move with the thread count and so would
/// every fitted value.
fn objective(
    panel: &Panel,
    binomial: &BinomialTable,
    fitted: &Fitted,
    variant_grid: &FrequencyGrid,
    carrier_grid: &FrequencyGrid,
    arm: Arm,
) -> f64 {
    let view = build_view(fitted, variant_grid, carrier_grid);
    let quiet = QuietTables::new(
        binomial,
        &view.cells,
        &view.carrier_cells,
        fitted.inbreeding,
        arm,
    );
    let terms: Vec<f64> = (0..panel.loci)
        .into_par_iter()
        .map_init(Scratch::default, |scratch, locus| {
            ln_locus(panel, locus, binomial, fitted, &view, &quiet, arm, scratch)
        })
        .collect();
    terms.iter().sum::<f64>() / panel.loci as f64
}

// ---------------------------------------------------------------------
// The search
// ---------------------------------------------------------------------

fn logit(p: f64) -> f64 {
    let p = p.clamp(1e-9, 1.0 - 1e-9);
    (p / (1.0 - p)).ln()
}

fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Golden-section maximisation of a one-dimensional slice, in a transformed coordinate so that a
/// rate never leaves `(0, 1)` and a shape never leaves `(0, ∞)`.
fn climb_scalar(mut score: impl FnMut(f64) -> f64, start: f64, span: f64) -> f64 {
    let phi = (5.0_f64.sqrt() - 1.0) / 2.0;
    let (mut lo, mut hi) = (start - span, start + span);
    let mut c = hi - phi * (hi - lo);
    let mut d = lo + phi * (hi - lo);
    let (mut fc, mut fd) = (score(c), score(d));
    for _ in 0..16 {
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
        if (hi - lo).abs() < 1e-5 {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// Which fitted number a coordinate climb is moving.
#[derive(Clone, Copy)]
enum Knob {
    EpsClean,
    EpsNoisy,
    NoisyShare,
    Invariant,
    FixedAlternative,
    BetaA,
    BetaB,
    DuplicatedShare,
    CarrierA,
    CarrierB,
    Inbreeding,
}

impl Knob {
    fn read(self, fitted: &Fitted) -> f64 {
        match self {
            Knob::EpsClean => logit(fitted.eps_clean),
            Knob::EpsNoisy => logit(fitted.eps_noisy),
            Knob::NoisyShare => logit(fitted.w_noisy),
            Knob::Invariant => logit(fitted.p_invariant),
            Knob::FixedAlternative => logit(fitted.p_fixed_alt),
            Knob::BetaA => fitted.beta_a.ln(),
            Knob::BetaB => fitted.beta_b.ln(),
            Knob::DuplicatedShare => logit(fitted.duplicated_share),
            Knob::CarrierA => fitted.carrier_a.ln(),
            Knob::CarrierB => fitted.carrier_b.ln(),
            Knob::Inbreeding => logit(fitted.inbreeding.clamp(1e-6, 1.0 - 1e-6)),
        }
    }
    fn write(self, fitted: &mut Fitted, x: f64) {
        match self {
            Knob::EpsClean => fitted.eps_clean = expit(x),
            Knob::EpsNoisy => fitted.eps_noisy = expit(x),
            Knob::NoisyShare => fitted.w_noisy = expit(x),
            Knob::Invariant => fitted.p_invariant = expit(x),
            Knob::FixedAlternative => fitted.p_fixed_alt = expit(x),
            Knob::BetaA => fitted.beta_a = x.exp(),
            Knob::BetaB => fitted.beta_b = x.exp(),
            Knob::DuplicatedShare => fitted.duplicated_share = expit(x),
            Knob::CarrierA => fitted.carrier_a = x.exp(),
            Knob::CarrierB => fitted.carrier_b = x.exp(),
            Knob::Inbreeding => fitted.inbreeding = expit(x),
        }
    }
    /// Whether moving this one invalidates the read-likelihood table.
    fn moves_error_rates(self) -> bool {
        matches!(self, Knob::EpsClean | Knob::EpsNoisy)
    }
}

/// Alternate over the coordinates until the objective stops moving, as `spec §3.3` does.
fn fit(
    panel: &Panel,
    arm: Arm,
    variant_grid: &FrequencyGrid,
    carrier_grid: &FrequencyGrid,
    variant: usize,
    passes: usize,
) -> (Fitted, f64) {
    let mut fitted = Fitted::start(variant, arm);
    let mut knobs = vec![
        Knob::EpsClean,
        Knob::EpsNoisy,
        Knob::NoisyShare,
        Knob::Invariant,
        Knob::FixedAlternative,
        Knob::BetaA,
        Knob::BetaB,
    ];
    if arm.has_third_class() {
        knobs.extend([Knob::DuplicatedShare, Knob::CarrierA, Knob::CarrierB]);
    }
    knobs.push(Knob::Inbreeding);

    let mut previous = f64::NEG_INFINITY;
    let mut table = BinomialTable::new(panel.deepest, [fitted.eps_clean, fitted.eps_noisy]);
    for _ in 0..passes {
        for &knob in &knobs {
            let base = fitted.clone();
            let best = climb_scalar(
                |x| {
                    let mut trial = base.clone();
                    knob.write(&mut trial, x);
                    if trial.p_segregating() <= 1e-9 {
                        return f64::NEG_INFINITY;
                    }
                    if knob.moves_error_rates() {
                        let table =
                            BinomialTable::new(panel.deepest, [trial.eps_clean, trial.eps_noisy]);
                        objective(panel, &table, &trial, variant_grid, carrier_grid, arm)
                    } else {
                        objective(panel, &table, &trial, variant_grid, carrier_grid, arm)
                    }
                },
                knob.read(&fitted),
                3.0,
            );
            knob.write(&mut fitted, best);
            if knob.moves_error_rates() {
                table = BinomialTable::new(panel.deepest, [fitted.eps_clean, fitted.eps_noisy]);
            }
        }
        let score = objective(panel, &table, &fitted, variant_grid, carrier_grid, arm);
        if (score - previous).abs() < 1e-11 {
            previous = score;
            break;
        }
        previous = score;
    }
    (fitted, previous)
}

// ---------------------------------------------------------------------
// What the fit says afterwards
// ---------------------------------------------------------------------

/// What a converged fit reports: the three headline numbers, and how it sorted the loci.
struct Verdict {
    observed_heterozygosity: f64,
    expected_heterozygosity: f64,
    inbreeding: f64,
    /// Drawn duplicated loci with at least one carrier, and how many of them the fit gives a
    /// posterior above a half of being duplicated.
    duplicated_loci: usize,
    duplicated_found: usize,
    /// The same, split by how many samples carry the duplication.
    found_by_carriers: [(usize, usize); 4],
    /// **The count that says how much of the damage is left**: how many (position, sample) pairs
    /// where a sample really carries a duplication sit at loci the fit does *not* call duplicated,
    /// out of all such pairs. A locus missed because one sample of fifty carries it leaves one pair
    /// behind; a locus missed because thirty do leaves thirty.
    missed_carrier_positions: u64,
    total_carrier_positions: u64,
    /// Drawn segregating loci, and how many the fit wrongly calls duplicated.
    segregating_loci: usize,
    segregating_called_duplicated: usize,
    /// Drawn monomorphic or fixed-alternative loci wrongly called duplicated.
    invariant_called_duplicated: usize,
    fitted: Fitted,
    score: f64,
}

/// The bands the found-count is split into: carried by one sample, by a handful, by a large
/// minority, by more than half the panel.
fn carrier_band(carriers: u32, samples: usize) -> usize {
    let share = f64::from(carriers) / samples as f64;
    if carriers <= 1 {
        0
    } else if share < 0.1 {
        1
    } else if share < 0.5 {
        2
    } else {
        3
    }
}

/// How to label each band, in carrier counts. A band can be empty on a small panel — under a tenth
/// of ten samples is one sample, which the first band already holds — and then it says so rather
/// than printing a backwards range.
fn carrier_band_labels(samples: usize) -> [String; 4] {
    let first = |share: f64| {
        (0..=samples)
            .find(|&c| c as f64 / samples as f64 >= share)
            .unwrap_or(samples)
    };
    let (tenth, half) = (first(0.1), first(0.5));
    let range = |low: usize, high: usize| {
        if low > high {
            "(none possible)".to_string()
        } else if low == high {
            format!("{low}")
        } else {
            format!("{low} to {high}")
        }
    };
    [
        "1".to_string(),
        range(2, tenth.saturating_sub(1)),
        range(tenth.max(2), half.saturating_sub(1)),
        format!("{} or more", half.max(2)),
    ]
}

/// **Observed heterozygosity is a mean of genotype posteriors, not a count** (spec §3.2), so it has
/// to be computed from the converged fit rather than read off anything.
///
/// The third class contributes nothing to it: under that class a sample is either a duplication
/// carrier — which is not a genotype — or homozygous reference. That is exactly how recognising the
/// class removes the inflation.
fn verdict(
    panel: &Panel,
    fitted: &Fitted,
    score: f64,
    variant_grid: &FrequencyGrid,
    carrier_grid: &FrequencyGrid,
    arm: Arm,
) -> Verdict {
    let view = build_view(fitted, variant_grid, carrier_grid);
    let table = BinomialTable::new(panel.deepest, [fitted.eps_clean, fitted.eps_noisy]);
    let quiet = QuietTables::new(
        &table,
        &view.cells,
        &view.carrier_cells,
        fitted.inbreeding,
        arm,
    );
    let class_weights = [1.0 - fitted.w_noisy, fitted.w_noisy];

    // Per locus: the summed heterozygous posterior over the samples, and the posterior that the
    // locus is duplicated.
    let per_locus: Vec<(f64, f64)> = (0..panel.loci)
        .into_par_iter()
        .map_init(Scratch::default, |scratch: &mut Scratch, locus| {
            let base = locus * panel.samples;
            let mut ln_terms = Vec::new();
            let mut ln_duplicated_terms = Vec::new();
            // `heterozygous[sample]` accumulates, in linear space, the unnormalised weight of the
            // sample being heterozygous.
            let mut heterozygous = vec![f64::NEG_INFINITY; panel.samples];
            for (class, &class_weight) in class_weights.iter().enumerate() {
                if class_weight <= 0.0 {
                    continue;
                }
                let ln_class = class_weight.ln()
                    + if arm.has_third_class() && fitted.duplicated_share > 0.0 {
                        (1.0 - fitted.duplicated_share).max(1e-300).ln()
                    } else {
                        0.0
                    };
                // The segregating branch, cell by cell, with a leave-one-out so each sample's own
                // heterozygous weight can be lifted out of the product.
                let segregating = fitted.p_segregating();
                if segregating > 0.0 {
                    for (cell, &(mass, _)) in view.cells.iter().enumerate() {
                        if mass <= 0.0 {
                            continue;
                        }
                        let prior = &view.ln_genotype_prior[cell];
                        let mut inner = Vec::with_capacity(panel.samples);
                        let mut total = 0.0;
                        for sample in 0..panel.samples {
                            let d = panel.depth[base + sample];
                            let k = panel.alt[base + sample];
                            let value = ln_sum_exp(&[
                                prior[0] + table.get(0, class, d, k),
                                prior[1] + table.get(1, class, d, k),
                                prior[2] + table.get(2, class, d, k),
                            ]);
                            inner.push(value);
                            total += value;
                        }
                        let branch = ln_class + segregating.ln() + mass.ln() + total;
                        ln_terms.push(branch);
                        for sample in 0..panel.samples {
                            let d = panel.depth[base + sample];
                            let k = panel.alt[base + sample];
                            let het = branch - inner[sample] + prior[1] + table.get(1, class, d, k);
                            heterozygous[sample] = ln_add(heterozygous[sample], het);
                        }
                    }
                }
                // The two masses carry no heterozygotes at all.
                if fitted.p_invariant > 0.0 {
                    let mut total = ln_class + fitted.p_invariant.ln();
                    for (depth, &count) in panel.quiet_at_depth[locus].iter().enumerate() {
                        if count > 0 {
                            total += f64::from(count) * quiet.all_reference[class][depth];
                        }
                    }
                    for &(d, k) in &panel.loud[locus] {
                        total += table.get(0, class, d, k);
                    }
                    ln_terms.push(total);
                }
                if fitted.p_fixed_alt > 0.0 {
                    let mut total = ln_class + fitted.p_fixed_alt.ln();
                    for (depth, &count) in panel.quiet_at_depth[locus].iter().enumerate() {
                        if count > 0 {
                            total += f64::from(count) * quiet.all_alternative[class][depth];
                        }
                    }
                    for &(d, k) in &panel.loud[locus] {
                        total += table.get(2, class, d, k);
                    }
                    ln_terms.push(total);
                }
                if arm.has_third_class() && fitted.duplicated_share > 0.0 {
                    let branch = class_weight.ln()
                        + fitted.duplicated_share.ln()
                        + ln_duplicated_class(panel, locus, &table, &view, &quiet, class, scratch);
                    ln_terms.push(branch);
                    ln_duplicated_terms.push(branch);
                }
            }
            let normaliser = ln_sum_exp(&ln_terms);
            let summed: f64 = heterozygous.iter().map(|&h| (h - normaliser).exp()).sum();
            let duplicated = if ln_duplicated_terms.is_empty() {
                0.0
            } else {
                (ln_sum_exp(&ln_duplicated_terms) - normaliser).exp()
            };
            (summed, duplicated)
        })
        .collect();

    let observed_heterozygosity =
        per_locus.iter().map(|&(h, _)| h).sum::<f64>() / (panel.loci as f64 * panel.samples as f64);
    let expected_heterozygosity = fitted.expected_heterozygosity();

    let mut duplicated_loci = 0;
    let mut duplicated_found = 0;
    let mut found_by_carriers = [(0_usize, 0_usize); 4];
    let mut segregating_loci = 0;
    let mut segregating_called_duplicated = 0;
    let mut invariant_called_duplicated = 0;
    let mut missed_carrier_positions = 0_u64;
    let mut total_carrier_positions = 0_u64;
    for (locus, &(_, duplicated_posterior)) in per_locus.iter().enumerate() {
        let called = duplicated_posterior > 0.5;
        match panel.kind[locus] {
            LocusKind::Duplicated if panel.carriers[locus] > 0 => {
                duplicated_loci += 1;
                total_carrier_positions += u64::from(panel.carriers[locus]);
                let band = carrier_band(panel.carriers[locus], panel.samples);
                found_by_carriers[band].0 += 1;
                if called {
                    duplicated_found += 1;
                    found_by_carriers[band].1 += 1;
                } else {
                    missed_carrier_positions += u64::from(panel.carriers[locus]);
                }
            }
            LocusKind::Segregating => {
                segregating_loci += 1;
                if called {
                    segregating_called_duplicated += 1;
                }
            }
            _ => {
                if called {
                    invariant_called_duplicated += 1;
                }
            }
        }
    }

    Verdict {
        observed_heterozygosity,
        expected_heterozygosity,
        inbreeding: 1.0 - observed_heterozygosity / expected_heterozygosity,
        duplicated_loci,
        duplicated_found,
        found_by_carriers,
        missed_carrier_positions,
        total_carrier_positions,
        segregating_loci,
        segregating_called_duplicated,
        invariant_called_duplicated,
        fitted: fitted.clone(),
        score,
    }
}

// ---------------------------------------------------------------------
// Reporting
// ---------------------------------------------------------------------

fn relative(got: f64, want: f64) -> f64 {
    if want.abs() < 1e-15 {
        f64::NAN
    } else {
        100.0 * (got - want) / want
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let samples: usize = args
        .next()
        .map_or(50, |a| a.parse().expect("a sample count"));
    let depth: f64 = args
        .next()
        .map_or(0.0, |a| a.parse().expect("a mean depth"));
    let loci: usize = args
        .next()
        .map_or(200_000, |a| a.parse().expect("a locus count"));
    let inbreeding: f64 = args
        .next()
        .map_or(0.6, |a| a.parse().expect("an inbreeding coefficient"));
    let window: f64 = args
        .next()
        .map_or(12_000.0, |a| a.parse().expect("aligned bases per window"));
    // How the carriers are spread: the share of duplicated loci every sample carries, then the two
    // shapes of the Beta the rest draw their carrier frequency from. The defaults are what the eight
    // tomato samples' window counts imply.
    let collapse_share: f64 = args
        .next()
        .map_or(REFERENCE_COLLAPSE_SHARE, |a| a.parse().expect("a share"));
    let carrier_a: f64 = args
        .next()
        .map_or(CARRIER_BETA.0, |a| a.parse().expect("a Beta shape"));
    let carrier_b: f64 = args
        .next()
        .map_or(CARRIER_BETA.1, |a| a.parse().expect("a Beta shape"));

    // Twenty-four cells for the allele-frequency Beta and ten for the carrier-frequency one — both
    // accuracy, not freedom. Eight is enough for the carrier frequency because nothing is read off its
    // shape; the arcsine grid's top cell runs from 0.96 to 1.0, which is where the reference's own
    // collapse sits.
    let variant_grid = FrequencyGrid::new(24);
    let carrier_grid = FrequencyGrid::new(8);

    let depths: Vec<f64> = if depth > 0.0 {
        vec![depth]
    } else {
        vec![3.0, 25.0]
    };

    for depth_mean in depths {
        let truth = Truth::tomato_like(
            samples,
            depth_mean,
            inbreeding,
            window,
            (collapse_share, carrier_a, carrier_b),
        );
        let started = Instant::now();
        let panel = Panel::draw(&truth, loci, 20_260_813);
        let hobs_truth = panel.drawn_observed_heterozygosity();
        let hexp_truth = panel.drawn_expected_heterozygosity;
        let inbreeding_truth = 1.0 - hobs_truth / hexp_truth;
        let carrier_rate = panel.drawn_carrier_rate();

        let duplicated_loci: usize = (0..loci)
            .filter(|&l| panel.kind[l] == LocusKind::Duplicated && panel.carriers[l] > 0)
            .count();
        let labels = carrier_band_labels(samples);
        let mut carrier_bands = [0_usize; 4];
        for locus in 0..loci {
            if panel.kind[locus] == LocusKind::Duplicated && panel.carriers[locus] > 0 {
                carrier_bands[carrier_band(panel.carriers[locus], samples)] += 1;
            }
        }

        println!(
            "\n==================================================================================\n\
             {samples} samples, {depth_mean} reads a site, {loci} loci, {window:.0} aligned bases \
             a coverage window\n\
             ==================================================================================\n\
             drawn truth   Hobs {hobs_truth:.4e}   Hexp {hexp_truth:.4e}   \
             inbreeding coefficient {inbreeding_truth:.4}\n\
             duplicated carrier positions {carrier_rate:.4e} — {:.0} per two million, against \
             {:.0} heterozygous per two million ({:.2} of them)\n\
             {duplicated_loci} duplicated loci carried by someone, by how many samples carry each: \
             {} carried by {}, {} by {}, {} by {}, {} by {}\n\
             carrier spread: {collapse_share} of duplicated loci carried by every sample, the rest \
             Beta({carrier_a}, {carrier_b})\n\
             drawn in {:.1}s",
            carrier_rate * 2e6,
            hobs_truth * 2e6,
            carrier_rate / hobs_truth,
            carrier_bands[0],
            labels[0],
            carrier_bands[1],
            labels[1],
            carrier_bands[2],
            labels[2],
            carrier_bands[3],
            labels[3],
            started.elapsed().as_secs_f64(),
        );

        println!(
            "\n{:>9} | {:>10} {:>8} | {:>10} {:>8} | {:>8} {:>8} | {:>10} | {:>7}",
            "fit", "Hobs", "vs truth", "Hexp", "vs truth", "F", "vs truth", "lnL", "seconds"
        );

        let mut verdicts = Vec::new();
        for arm in [Arm::Coverage, Arm::Pattern, Arm::NoClass] {
            let started = Instant::now();
            // Three starting points, best score taken — spec §11 question 3 measured the start to
            // be worth 46% of the clean error rate, more than any other choice here.
            let (fitted, score) = (0..3)
                .map(|variant| fit(&panel, arm, &variant_grid, &carrier_grid, variant, 10))
                .max_by(|a, b| a.1.partial_cmp(&b.1).expect("finite scores"))
                .expect("three starts");
            let verdict = verdict(&panel, &fitted, score, &variant_grid, &carrier_grid, arm);
            println!(
                "{:>9} | {:>10.4e} {:>+7.1}% | {:>10.4e} {:>+7.1}% | {:>8.4} {:>+7.1}% | \
                 {:>10.5} | {:>7.1}",
                arm.name(),
                verdict.observed_heterozygosity,
                relative(verdict.observed_heterozygosity, hobs_truth),
                verdict.expected_heterozygosity,
                relative(verdict.expected_heterozygosity, hexp_truth),
                verdict.inbreeding,
                relative(verdict.inbreeding, inbreeding_truth),
                verdict.score,
                started.elapsed().as_secs_f64(),
            );
            verdicts.push((arm, verdict));
        }

        println!(
            "\n{:>9} | {:>10} | {:>10} | {:>10} | {:>8} | {:>8} | {:>8}",
            "fit", "eps_clean", "eps_noisy", "w_noisy", "dup share", "carrier a", "carrier b"
        );
        for (arm, v) in &verdicts {
            println!(
                "{:>9} | {:>+9.1}% | {:>+9.1}% | {:>+9.1}% | {:>8.2e} | {:>8.3} | {:>8.3}",
                arm.name(),
                relative(v.fitted.eps_clean, truth.eps_clean),
                relative(v.fitted.eps_noisy, truth.eps_noisy),
                relative(v.fitted.w_noisy, truth.w_noisy),
                v.fitted.duplicated_share,
                v.fitted.carrier_a,
                v.fitted.carrier_b,
            );
        }
        println!(
            "the third class's true weight is {:.2e}",
            truth.duplicated_share
        );

        println!(
            "\nwhich loci each fit calls duplicated (posterior above a half)\n\
             {:>9} | {:>16} | {:>18} | {:>18}",
            "fit", "duplicated found", "variants called dup", "invariants called dup"
        );
        for (arm, v) in &verdicts {
            if !arm.has_third_class() {
                continue;
            }
            println!(
                "{:>9} | {:>7} of {:<6} | {:>8} of {:<7} | {:>21}",
                arm.name(),
                v.duplicated_found,
                v.duplicated_loci,
                v.segregating_called_duplicated,
                v.segregating_loci,
                v.invariant_called_duplicated,
            );
        }

        println!(
            "\nfound, split by how many samples carry the duplication\n\
             {:>9} | {:>14} | {:>14} | {:>14} | {:>14}",
            "fit", labels[0], labels[1], labels[2], labels[3],
        );
        for (arm, v) in &verdicts {
            if !arm.has_third_class() {
                continue;
            }
            println!(
                "{:>9} | {:>6} of {:<5} | {:>6} of {:<5} | {:>6} of {:<5} | {:>6} of {:<5}",
                arm.name(),
                v.found_by_carriers[0].1,
                v.found_by_carriers[0].0,
                v.found_by_carriers[1].1,
                v.found_by_carriers[1].0,
                v.found_by_carriers[2].1,
                v.found_by_carriers[2].0,
                v.found_by_carriers[3].1,
                v.found_by_carriers[3].0,
            );
        }

        println!(
            "\nhow much of the artefact each fit leaves behind\n\
             {:>9} | {:>34}",
            "fit", "carrier positions at missed loci"
        );
        for (arm, v) in &verdicts {
            if !arm.has_third_class() {
                continue;
            }
            println!(
                "{:>9} | {:>10} of {:<8} ({:>5.1}%)",
                arm.name(),
                v.missed_carrier_positions,
                v.total_carrier_positions,
                100.0 * v.missed_carrier_positions as f64 / v.total_carrier_positions as f64,
            );
        }
    }
}
