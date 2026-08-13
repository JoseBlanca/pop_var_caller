//! **Does an uncontaminated panel of landraces come back uncontaminated?**
//!
//! `spec/parameter_prepass_joint_fit.md` §3.4 fits one contamination fraction per sample:
//! the share of a sample's reads that came from another individual. It is identified by
//! three signatures, and the load-bearing one is that **a contaminant's alleles are alleles
//! the population carries** — where a sequencing error's wrong base is one of three at
//! random. So the estimate leans on the fit's own answer to *which loci vary and at what
//! frequency*.
//!
//! §3.4.2 records the hazard, from `verifyBamID2`: if the frequency used is the **pooled**
//! one, a sample from a diverged subpopulation carries alleles the pooled spectrum calls
//! rare — and rare alleles turning up in a sample is exactly the contamination signature.
//! **Structure would be read as contamination.** The tomato panel is landraces from several
//! regions, so it has precisely that structure.
//!
//! **Nothing in the spec has a number behind it, and the number that matters is not that a
//! contaminated sample returns its own fraction. It is that a structured, clean panel
//! returns zero.** This program measures both.
//!
//! ## The two arms
//!
//! - **`pooled`** — every sample's alleles are scored against one frequency per locus,
//!   estimated from the whole panel. This is what a route that fits one spectrum does.
//! - **`by-group`** — each sample's alleles are scored against the frequency of **its own
//!   subpopulation**, estimated from that subpopulation's members. This is an *oracle*: it
//!   is handed the group membership that `verifyBamID2` has to infer from principal
//!   components. **It is therefore a ceiling on the fix rather than the fix** — if even this
//!   cannot return zero, no principal-component model will.
//!
//! ## How a structured panel is drawn
//!
//! Balding–Nichols: each locus has an ancestral allele frequency, and each subpopulation's
//! own frequency is drawn around it with a spread set by `F_st` — 0 is one panmictic
//! population, 0.15 is about as diverged as landraces of one crop get.
//!
//! ```text
//! ng_joint_contamination_harness [null|spike|budget] [samples] [depth] [loci] [groups] [fst]
//! ```

use std::time::Instant;

// ---------------------------------------------------------------------
// Small numerics
// ---------------------------------------------------------------------

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

    fn gamma(&mut self, shape: f64) -> f64 {
        if shape < 1.0 {
            let u = self.uniform().max(1e-300);
            return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let u1 = self.uniform().max(1e-300);
            let u2 = self.uniform();
            let x = (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos();
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

    fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        x / (x + y).max(1e-300)
    }

    fn binomial(&mut self, n: u32, p: f64) -> u32 {
        (0..n).filter(|_| self.uniform() < p).count() as u32
    }
}

fn ln_binomial(k: u32, n: u32, p: f64) -> f64 {
    let p = p.clamp(1e-12, 1.0 - 1e-12);
    f64::from(k) * p.ln() + f64::from(n - k) * (1.0 - p).ln()
}

fn ln_sum_exp(values: &[f64]) -> f64 {
    let largest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !largest.is_finite() {
        return largest;
    }
    largest + values.iter().map(|v| (v - largest).exp()).sum::<f64>().ln()
}

// ---------------------------------------------------------------------
// The panel
// ---------------------------------------------------------------------

/// What one sample showed at one locus.
#[derive(Copy, Clone)]
struct Observation {
    depth: u32,
    alternative: u32,
}

struct Panel {
    /// `data[locus][sample]`.
    data: Vec<Vec<Observation>>,
    /// Which subpopulation each sample belongs to.
    group_of: Vec<usize>,
    groups: usize,
    /// The contamination fraction each sample was drawn with.
    truth: Vec<f64>,
    /// Loci at which the panel's chromosomes are not all the same allele — **the markers
    /// contamination can actually use**, since a contaminant carries the same allele as its
    /// host wherever the population is fixed.
    segregating: usize,
    /// `true_frequencies[locus][group]` — what the genotypes were actually drawn at.
    /// **Not available to any real fit**; it is here so that a frequency that is right but
    /// noisily estimated can be told from one that is wrong.
    true_frequencies: Vec<Vec<f64>>,
}

#[derive(Clone)]
struct Truth {
    samples: usize,
    depth: u32,
    loci: usize,
    groups: usize,
    /// How diverged the subpopulations are. 0 is one panmictic population.
    fst: f64,
    error: f64,
    inbreeding: f64,
    /// The contamination fraction of each sample; all zero in the test that matters.
    contamination: Vec<f64>,
    /// How many samples each subpopulation holds. `None` splits them evenly.
    ///
    /// **The unequal case is the one to watch**: the axes of a decomposition go to the
    /// largest groups first, so a small divergent set can fail to get one and fall back to
    /// something near the panel average.
    group_sizes: Option<Vec<usize>>,
}

impl Truth {
    fn clean(samples: usize, depth: u32, loci: usize, groups: usize, fst: f64) -> Self {
        Self {
            samples,
            depth,
            loci,
            groups,
            fst,
            error: 0.002,
            inbreeding: 0.5,
            contamination: vec![0.0; samples],
            group_sizes: None,
        }
    }

    /// Which subpopulation each sample belongs to.
    fn group_of(&self) -> Vec<usize> {
        match &self.group_sizes {
            None => (0..self.samples)
                .map(|i| i * self.groups / self.samples)
                .collect(),
            Some(sizes) => {
                let mut out = Vec::with_capacity(self.samples);
                for (group, count) in sizes.iter().enumerate() {
                    out.extend(std::iter::repeat_n(group, *count));
                }
                assert_eq!(
                    out.len(),
                    self.samples,
                    "the group sizes must cover the panel"
                );
                out
            }
        }
    }

    /// The probability a read shows the alternative allele, given a genotype's alternative
    /// copy count.
    fn read_probability(&self, copies: u32) -> f64 {
        match copies {
            0 => self.error,
            1 => 0.5,
            _ => 1.0 - self.error,
        }
    }

    fn draw(&self, seed: u64) -> Panel {
        let mut rng = Rng(seed);
        let group_of = self.group_of();
        let mut data = Vec::with_capacity(self.loci);
        let mut true_frequencies = Vec::with_capacity(self.loci);
        let mut segregating = 0_usize;

        for _ in 0..self.loci {
            // Balding–Nichols: an ancestral frequency, then one frequency per subpopulation
            // drawn around it with a spread set by F_st.
            let ancestral = rng.beta(0.3, 1.2);
            let group_frequencies: Vec<f64> = (0..self.groups)
                .map(|_| {
                    if self.fst <= 0.0 {
                        ancestral
                    } else {
                        let scale = (1.0 - self.fst) / self.fst;
                        rng.beta(
                            (ancestral * scale).max(1e-3),
                            ((1.0 - ancestral) * scale).max(1e-3),
                        )
                    }
                })
                .collect();

            let genotype = |rng: &mut Rng, frequency: f64| -> u32 {
                // Hardy–Weinberg with inbreeding.
                if rng.uniform() < self.inbreeding {
                    if rng.uniform() < frequency { 2 } else { 0 }
                } else {
                    u32::from(rng.uniform() < frequency) + u32::from(rng.uniform() < frequency)
                }
            };

            let mut row = Vec::with_capacity(self.samples);
            let mut copies_seen = 0_u32;
            let mut chromosomes = 0_u32;
            for sample in 0..self.samples {
                let frequency = group_frequencies[group_of[sample]];
                let own = genotype(&mut rng, frequency);
                let alpha = self.contamination[sample];
                // The contaminant is another individual of the same subpopulation — a second
                // plant in the tube or a neighbouring library on the same run, which is the
                // case `verifyBamID2`'s equal-ancestry model covers.
                let other = if alpha > 0.0 {
                    genotype(&mut rng, frequency)
                } else {
                    own
                };
                let mut alternative = 0;
                for _ in 0..self.depth {
                    let from = if rng.uniform() < alpha { other } else { own };
                    alternative += rng.binomial(1, self.read_probability(from));
                }
                copies_seen += own;
                chromosomes += 2;
                row.push(Observation {
                    depth: self.depth,
                    alternative,
                });
            }
            if copies_seen > 0 && copies_seen < chromosomes {
                segregating += 1;
            }
            data.push(row);
            true_frequencies.push(group_frequencies);
        }

        Panel {
            data,
            group_of,
            groups: self.groups,
            truth: self.contamination.clone(),
            segregating,
            true_frequencies,
        }
    }
}

// ---------------------------------------------------------------------
// Fitting one contamination fraction per sample
// ---------------------------------------------------------------------

#[derive(Copy, Clone, PartialEq, Eq)]
enum Frequencies {
    /// One frequency per locus, from the whole panel.
    Pooled,
    /// One frequency per locus **per subpopulation**, from that subpopulation's members.
    /// It is handed the membership `verifyBamID2` infers from principal components — but it
    /// still has to *estimate* the frequency, from a twelfth of the panel.
    ByGroup,
    /// The frequencies the genotypes were actually drawn at. **No fit can have these**; the
    /// arm exists to separate a frequency that is wrong from one that is right and noisy.
    TrueByGroup,
    /// One frequency per **sample** per locus, fitted as a straight line in the samples'
    /// ancestry coordinates. **The thing being measured**: no groups, no assignment, and
    /// every locus's line fitted from the whole panel.
    PcRegression {
        components: usize,
        /// Shrink each locus's slopes towards zero — towards the pooled frequency —
        /// by how much of that locus's dosage spread the line actually explains.
        ///
        /// **A line fitted from fifty samples is itself an estimate**, and an error in the
        /// frequency inflates `α`: contamination is the parameter that absorbs reads which
        /// do not fit the genotype the prior expected, so a frequency that is noisily wrong
        /// manufactures some. Shrinking trades that noise for a little of the pooled
        /// frequency's bias, and §5 measures which is worth more.
        shrink: bool,
    },
}

/// The allele frequency each sample's genotype is scored against, per locus.
///
/// **Estimated from the reads, not supplied**, because that is what the fit has: the moment
/// estimator inverts the error rate out of the alternative-read share, which is what a
/// spectrum fitted over the same loci converges to.
fn estimate_frequencies(panel: &Panel, which: Frequencies, error: f64) -> Vec<Vec<f64>> {
    if let Frequencies::PcRegression { components, shrink } = which {
        return pc_regression_frequencies(panel, error, components, shrink);
    }
    let samples = panel.data[0].len();
    let mut out = Vec::with_capacity(panel.data.len());
    for row in &panel.data {
        let mut per_group = vec![(0_u64, 0_u64); panel.groups.max(1)];
        let (mut all_alt, mut all_depth) = (0_u64, 0_u64);
        for (sample, observation) in row.iter().enumerate() {
            all_alt += u64::from(observation.alternative);
            all_depth += u64::from(observation.depth);
            let group = panel.group_of[sample];
            per_group[group].0 += u64::from(observation.alternative);
            per_group[group].1 += u64::from(observation.depth);
        }
        let invert = |alt: u64, depth: u64| -> f64 {
            if depth == 0 {
                return 0.0;
            }
            let share = alt as f64 / depth as f64;
            ((share - error) / (1.0 - 2.0 * error)).clamp(1e-4, 1.0 - 1e-4)
        };
        let per_group: Vec<f64> = match which {
            Frequencies::Pooled => vec![invert(all_alt, all_depth); panel.groups.max(1)],
            Frequencies::ByGroup => per_group
                .iter()
                .map(|(alt, depth)| invert(*alt, *depth))
                .collect(),
            Frequencies::TrueByGroup => panel.true_frequencies[out.len()]
                .iter()
                .map(|f| f.clamp(1e-4, 1.0 - 1e-4))
                .collect(),
            Frequencies::PcRegression { .. } => unreachable!("handled above"),
        };
        out.push(
            (0..samples)
                .map(|sample| per_group[panel.group_of[sample]])
                .collect(),
        );
    }
    out
}

/// One expected allele frequency per sample per locus, fitted as a straight line in the
/// samples' own ancestry coordinates.
///
/// **This is the thing being measured**, and it has three steps.
///
/// 1. **A dosage per sample per locus** — the posterior mean number of alternative copies
///    under a pooled-frequency prior. At three reads a site a raw read fraction is far too
///    noisy to decompose; the prior is what makes it usable, and it is the same bootstrap
///    `PCAngsd` starts from.
/// 2. **Coordinates**, from the eigenvectors of the samples' own similarity matrix, on
///    Patterson-normalised dosages. A few numbers per sample and nothing else is kept: no
///    thresholds on the axes, no groups.
/// 3. **A straight line per locus**, `dosage ≈ b₀ + Σ bₖ·coordinateₖ`, fitted across **all**
///    samples. Each sample's own frequency is that line's height at its own coordinates,
///    halved. So a sample alone at one end of an axis borrows the *slope* — measured from
///    the whole panel — rather than any neighbour's allele counts.
fn pc_regression_frequencies(
    panel: &Panel,
    error: f64,
    components: usize,
    shrink: bool,
) -> Vec<Vec<f64>> {
    let loci = panel.data.len();
    let samples = panel.data[0].len();
    let (dosage, pooled, coordinates) = dosages_and_coordinates(panel, error, components);
    let _ = (loci, samples);
    pc_lines(&dosage, &pooled, &coordinates, components, shrink)
}

/// Steps 1 and 2 of [`pc_regression_frequencies`], kept separate because the **leverage**
/// each sample has on the fitted lines depends only on the coordinates — and it is one
/// number per sample for the whole run, not one per locus.
fn dosages_and_coordinates(
    panel: &Panel,
    error: f64,
    components: usize,
) -> (Vec<Vec<f64>>, Vec<f64>, Vec<Vec<f64>>) {
    let loci = panel.data.len();
    let samples = panel.data[0].len();

    // ---- 1. dosages ----------------------------------------------------------------
    let mut dosage = vec![vec![0.0_f64; samples]; loci];
    let mut pooled = vec![0.0_f64; loci];
    for (locus, row) in panel.data.iter().enumerate() {
        let (mut alt, mut depth) = (0_u64, 0_u64);
        for observation in row {
            alt += u64::from(observation.alternative);
            depth += u64::from(observation.depth);
        }
        let share = if depth == 0 {
            0.0
        } else {
            alt as f64 / depth as f64
        };
        let frequency = ((share - error) / (1.0 - 2.0 * error)).clamp(1e-3, 1.0 - 1e-3);
        pooled[locus] = frequency;
        // Hardy-Weinberg prior at the pooled frequency; inbreeding is left out here
        // deliberately, since the dosages only have to order the samples along an axis.
        let priors = [
            (1.0 - frequency) * (1.0 - frequency),
            2.0 * frequency * (1.0 - frequency),
            frequency * frequency,
        ];
        for (sample, observation) in row.iter().enumerate() {
            let mut weights = [0.0_f64; 3];
            for (copies, weight) in weights.iter_mut().enumerate() {
                let p = match copies {
                    0 => error,
                    1 => 0.5,
                    _ => 1.0 - error,
                };
                *weight = priors[copies]
                    * ln_binomial(observation.alternative, observation.depth, p).exp();
            }
            let total: f64 = weights.iter().sum();
            dosage[locus][sample] = if total > 0.0 {
                (weights[1] + 2.0 * weights[2]) / total
            } else {
                2.0 * frequency
            };
        }
    }

    // ---- 2. coordinates ------------------------------------------------------------
    // The samples' similarity matrix, on dosages centred and scaled the way a population
    // structure analysis scales them — dividing by sqrt(p(1-p)) so a rare allele's
    // contribution is not swamped by a common one.
    let mut gram = vec![vec![0.0_f64; samples]; samples];
    for locus in 0..loci {
        let mean = 2.0 * pooled[locus];
        let scale = (pooled[locus] * (1.0 - pooled[locus])).sqrt().max(1e-6);
        let centred: Vec<f64> = (0..samples)
            .map(|s| (dosage[locus][s] - mean) / scale)
            .collect();
        for i in 0..samples {
            for j in i..samples {
                gram[i][j] += centred[i] * centred[j];
            }
        }
    }
    for i in 0..samples {
        for j in 0..i {
            gram[i][j] = gram[j][i];
        }
    }
    let coordinates = leading_eigenvectors(&gram, components);
    (dosage, pooled, coordinates)
}

/// Step 3: a line per locus.
fn pc_lines(
    dosage: &[Vec<f64>],
    pooled: &[f64],
    coordinates: &[Vec<f64>],
    components: usize,
    shrink: bool,
) -> Vec<Vec<f64>> {
    let loci = dosage.len();
    let samples = coordinates.len();
    let width = components + 1;
    let mut out = Vec::with_capacity(loci);
    for locus in 0..loci {
        // Normal equations for `dosage ~ 1 + coordinates`, over all samples.
        let mut xtx = vec![vec![0.0_f64; width]; width];
        let mut xty = vec![0.0_f64; width];
        for sample in 0..samples {
            let mut design = vec![1.0_f64; width];
            for k in 0..components {
                design[k + 1] = coordinates[sample][k];
            }
            for a in 0..width {
                for b in 0..width {
                    xtx[a][b] += design[a] * design[b];
                }
                xty[a] += design[a] * dosage[locus][sample];
            }
        }
        // A ridge term, small beside the 50 samples on the diagonal, so a locus whose
        // dosages are constant cannot make the system singular.
        for (a, row) in xtx.iter_mut().enumerate() {
            row[a] += 1e-6;
        }
        let mut beta = solve(xtx, xty).unwrap_or_else(|| {
            let mut fallback = vec![0.0; width];
            fallback[0] = 2.0 * pooled[locus];
            fallback
        });
        if shrink && components > 0 {
            // The positive-part James-Stein factor: how much of the dosage spread the line
            // explains, against how much noise alone would explain. A locus whose slopes
            // are indistinguishable from noise keeps only its intercept, which is the
            // pooled frequency — so shrinking never does worse than not modelling
            // structure at that locus.
            let mean: f64 = dosage[locus].iter().sum::<f64>() / samples as f64;
            let mut explained = 0.0;
            let mut residual = 0.0;
            for sample in 0..samples {
                let mut fitted = beta[0];
                for k in 0..components {
                    fitted += beta[k + 1] * coordinates[sample][k];
                }
                explained += (fitted - mean) * (fitted - mean);
                residual += (dosage[locus][sample] - fitted) * (dosage[locus][sample] - fitted);
            }
            let noise = residual / (samples - width) as f64;
            let factor = if explained > 0.0 {
                (1.0 - components as f64 * noise / explained).max(0.0)
            } else {
                0.0
            };
            for slope in beta.iter_mut().skip(1) {
                *slope *= factor;
            }
            // The intercept moves back to the panel mean by whatever the slopes gave up.
            let mut centre = 0.0;
            for sample in 0..samples {
                for k in 0..components {
                    centre += beta[k + 1] * coordinates[sample][k];
                }
            }
            beta[0] = mean - centre / samples as f64;
        }
        out.push(
            (0..samples)
                .map(|sample| {
                    let mut fitted = beta[0];
                    for k in 0..components {
                        fitted += beta[k + 1] * coordinates[sample][k];
                    }
                    (fitted / 2.0).clamp(1e-4, 1.0 - 1e-4)
                })
                .collect(),
        );
    }
    out
}

/// **How much of its own fitted allele frequency each sample supplies.**
///
/// The line at every locus is fitted against the same coordinates, so this is one number per
/// sample for the whole run rather than one per locus, and it can be computed before a single
/// locus is touched. It runs from `(components + 1) / samples` — a sample pulling its fair
/// share — up towards 1, where the line at that sample's position is determined by that
/// sample alone and its "expected" frequency is really its own noisy reading.
///
/// **This is the number that says whose contamination estimate to trust.** An accession sitting
/// alone at the end of an axis has a fitted frequency that is mostly its own echo, and §5's
/// mechanism then applies: a noisy frequency manufactures contamination.
fn coordinate_leverage(coordinates: &[Vec<f64>], components: usize) -> Vec<f64> {
    let samples = coordinates.len();
    let width = components + 1;
    let design = |sample: usize| -> Vec<f64> {
        let mut row = vec![1.0_f64; width];
        for k in 0..components {
            row[k + 1] = coordinates[sample][k];
        }
        row
    };
    let mut xtx = vec![vec![0.0_f64; width]; width];
    for sample in 0..samples {
        let row = design(sample);
        for a in 0..width {
            for b in 0..width {
                xtx[a][b] += row[a] * row[b];
            }
        }
    }
    for (a, row) in xtx.iter_mut().enumerate() {
        row[a] += 1e-9;
    }
    (0..samples)
        .map(|sample| {
            let row = design(sample);
            solve(xtx.clone(), row.clone())
                .map(|z| row.iter().zip(&z).map(|(x, z)| x * z).sum())
                .unwrap_or(1.0)
        })
        .collect()
}

/// The `wanted` eigenvectors of largest eigenvalue of a small symmetric matrix, by cyclic
/// Jacobi rotation. The matrix is one row per sample, so 50 by 50 — this is nothing.
fn leading_eigenvectors(matrix: &[Vec<f64>], wanted: usize) -> Vec<Vec<f64>> {
    let n = matrix.len();
    let mut a: Vec<Vec<f64>> = matrix.to_vec();
    let mut v: Vec<Vec<f64>> = (0..n)
        .map(|i| (0..n).map(|j| if i == j { 1.0 } else { 0.0 }).collect())
        .collect();
    for _ in 0..60 {
        let mut off = 0.0;
        for i in 0..n {
            for j in i + 1..n {
                off += a[i][j] * a[i][j];
            }
        }
        if off < 1e-12 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                if a[p][q].abs() < 1e-14 {
                    continue;
                }
                let theta = 0.5 * (a[q][q] - a[p][p]) / a[p][q];
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|x, y| a[*y][*y].partial_cmp(&a[*x][*x]).expect("no NaN"));
    (0..n)
        .map(|sample| {
            order
                .iter()
                .take(wanted)
                .map(|axis| v[sample][*axis])
                .collect()
        })
        .collect()
}

/// Gaussian elimination with partial pivoting, for the `(components + 1)`-wide normal
/// equations. `None` where the system is singular.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    for column in 0..n {
        let pivot = (column..n).max_by(|x, y| {
            a[*x][column]
                .abs()
                .partial_cmp(&a[*y][column].abs())
                .expect("no NaN")
        })?;
        if a[pivot][column].abs() < 1e-12 {
            return None;
        }
        a.swap(column, pivot);
        b.swap(column, pivot);
        for row in column + 1..n {
            let factor = a[row][column] / a[column][column];
            for k in column..n {
                a[row][k] -= factor * a[column][k];
            }
            b[row] -= factor * b[column];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let mut total = b[row];
        for k in row + 1..n {
            total -= a[row][k] * x[k];
        }
        x[row] = total / a[row][row];
    }
    Some(x)
}

fn genotype_priors(frequency: f64, inbreeding: f64) -> [f64; 3] {
    let (p, q) = (frequency, 1.0 - frequency);
    [
        q * q + inbreeding * p * q,
        2.0 * p * q * (1.0 - inbreeding),
        p * p + inbreeding * p * q,
    ]
}

/// `ln P(this sample's reads at this locus | α)`, both genotypes summed over.
///
/// **The sequence-only two-genotype mixture** of `verifyBamID` (Jun et al. 2012): a read
/// comes from the intended sample's genotype with probability `1 − α` and from the
/// contaminant's with probability `α`, and neither genotype is known.
fn ln_locus(
    observation: Observation,
    frequency: f64,
    inbreeding: f64,
    error: f64,
    alpha: f64,
) -> f64 {
    let priors = genotype_priors(frequency, inbreeding);
    let read_probability = |copies: usize| match copies {
        0 => error,
        1 => 0.5,
        _ => 1.0 - error,
    };
    let mut terms = Vec::with_capacity(9);
    for own in 0..3 {
        if priors[own] <= 0.0 {
            continue;
        }
        for other in 0..3 {
            if priors[other] <= 0.0 {
                continue;
            }
            let p = (1.0 - alpha) * read_probability(own) + alpha * read_probability(other);
            terms.push(
                priors[own].ln()
                    + priors[other].ln()
                    + ln_binomial(observation.alternative, observation.depth, p),
            );
        }
    }
    ln_sum_exp(&terms)
}

/// Fit one sample's contamination fraction.
///
/// **The search is restricted to `α ≤ ½`**: the sequence-only likelihood is symmetric, so it
/// cannot tell a 20% contaminated sample from an 80% one, and a sample swap is invisible to
/// it by construction (§3.4.1).
fn fit_alpha(
    panel: &Panel,
    frequencies: &[Vec<f64>],
    sample: usize,
    inbreeding: f64,
    error: f64,
) -> f64 {
    let score = |alpha: f64| -> f64 {
        panel
            .data
            .iter()
            .zip(frequencies)
            .map(|(row, per_sample)| {
                ln_locus(row[sample], per_sample[sample], inbreeding, error, alpha)
            })
            .sum::<f64>()
    };
    // Golden section on [0, 0.5]. A grid first, because the likelihood is flat near zero and
    // a bracket that opens in the middle can walk away from a true zero.
    const PHI: f64 = 0.618_033_988_749_895;
    let (mut low, mut high) = (0.0_f64, 0.5_f64);
    let (mut c, mut d) = (high - PHI * (high - low), low + PHI * (high - low));
    let (mut fc, mut fd) = (score(c), score(d));
    for _ in 0..40 {
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
    let best = if fc > fd { c } else { d };
    // The bracket cannot reach 0 exactly; if the boundary scores higher, that is the answer.
    if score(0.0) >= score(best) { 0.0 } else { best }
}

struct Fitted {
    alphas: Vec<f64>,
}

impl Fitted {
    fn of(panel: &Panel, which: Frequencies, inbreeding: f64, error: f64) -> Self {
        let frequencies = estimate_frequencies(panel, which, error);
        Self {
            alphas: (0..panel.data[0].len())
                .map(|sample| fit_alpha(panel, &frequencies, sample, inbreeding, error))
                .collect(),
        }
    }

    fn mean(&self) -> f64 {
        self.alphas.iter().sum::<f64>() / self.alphas.len() as f64
    }

    fn max(&self) -> f64 {
        self.alphas.iter().copied().fold(0.0, f64::max)
    }

    /// How many samples a 1% flagging threshold would call contaminated.
    fn flagged(&self, threshold: f64) -> usize {
        self.alphas.iter().filter(|a| **a >= threshold).count()
    }
}

// ---------------------------------------------------------------------
// The runs
// ---------------------------------------------------------------------

/// The arms every run reports, in one place so the header and the cells cannot drift apart.
fn arms() -> Vec<Frequencies> {
    vec![
        Frequencies::Pooled,
        Frequencies::PcRegression {
            components: 4,
            shrink: false,
        },
        Frequencies::PcRegression {
            components: 4,
            shrink: true,
        },
        Frequencies::PcRegression {
            components: 8,
            shrink: true,
        },
        Frequencies::TrueByGroup,
    ]
}

fn arm_name(which: Frequencies) -> String {
    match which {
        Frequencies::Pooled => "pooled".to_string(),
        Frequencies::ByGroup => "by-group".to_string(),
        Frequencies::TrueByGroup => "true group".to_string(),
        Frequencies::PcRegression { components, shrink } => {
            format!("{components} axes{}", if shrink { ", shrunk" } else { "" })
        }
    }
}

fn arm_header() -> String {
    arms()
        .into_iter()
        .map(|which| format!("{:>21}", format!("{}: fitted / others", arm_name(which))))
        .collect::<Vec<_>>()
        .join(" |")
}

fn arm_header_wide() -> String {
    arms()
        .into_iter()
        .map(|which| format!("{:>24}", format!("{}: mean max flagged", arm_name(which))))
        .collect::<Vec<_>>()
        .join(" |")
}

/// **The measurement that matters**: a clean panel, swept over how diverged its
/// subpopulations are.
fn null_sweep(samples: usize, depth: u32, loci: usize, groups: usize) {
    println!(
        "\nA clean panel — every sample's true contamination is zero.\n\
         {samples} samples in {groups} subpopulations, {depth} reads a site, {loci} loci."
    );
    println!("\n  F_st |{}", arm_header_wide());
    for fst in [0.0_f64, 0.02, 0.05, 0.10, 0.20] {
        let truth = Truth::clean(samples, depth, loci, groups, fst);
        let panel = truth.draw(4242);
        let mut cells = Vec::new();
        for which in arms() {
            let fitted = Fitted::of(&panel, which, truth.inbreeding, truth.error);
            cells.push(format!(
                "{:>8.4} {:>7.4} {:>7}",
                fitted.mean(),
                fitted.max(),
                format!("{}/{}", fitted.flagged(0.01), samples)
            ));
        }
        println!("  {fst:>4.2} |{}", cells.join(" |"));
    }
    println!(
        "\n  A sample flagged here is a clean sample called contaminated. Under a 1% threshold\n  \
         that is a false positive, and under §3.4.5 it is a run the user is told to repeat."
    );
}

/// One genuinely contaminated sample in an otherwise clean panel: does its own fraction come
/// back, and do the others stay at zero?
fn spike(samples: usize, depth: u32, loci: usize, groups: usize, fst: f64) {
    println!(
        "\nOne contaminated sample among {samples}, {groups} subpopulations at F_st {fst:.2}, \
         {depth} reads a site, {loci} loci."
    );
    println!("\n  true alpha |{}", arm_header());
    for alpha in [0.01_f64, 0.03, 0.10] {
        let mut truth = Truth::clean(samples, depth, loci, groups, fst);
        truth.contamination[0] = alpha;
        let panel = truth.draw(4242);
        assert_eq!(
            panel.truth[0], alpha,
            "the spiked sample is the one that was spiked"
        );
        let mut cells = Vec::new();
        for which in arms() {
            let fitted = Fitted::of(&panel, which, truth.inbreeding, truth.error);
            cells.push(format!(
                "{:>10.4} {:>10.4}",
                fitted.alphas[0],
                fitted.alphas[1..].iter().copied().fold(0.0_f64, f64::max)
            ));
        }
        println!("  {alpha:>10.3} | {}", cells.join(" |"));
    }
}

/// How many loci the estimate needs — the sweep that prices §3.4.4's budget.
fn budget_sweep(samples: usize, depth: u32, groups: usize, fst: f64) {
    println!(
        "\nHow many loci a contamination estimate needs — {samples} samples, {groups} \
         subpopulations at F_st {fst:.2}, {depth} reads a site."
    );
    println!(
        "\n  The frequency each sample is scored against is its own subpopulation's, correct \
         rather than estimated —\n  the ceiling of §3.4.2's fix, so this prices the budget the \
         estimator needs and not the frequency's own error.\n"
    );
    println!(
        "     loci  segregating | clean panel: mean   max | one sample at 3%: fitted   others' max"
    );
    for loci in [5_000_usize, 20_000, 80_000, 320_000] {
        let started = Instant::now();
        let clean = Truth::clean(samples, depth, loci, groups, fst);
        let panel = clean.draw(4242);
        let panel_segregating = panel.segregating;
        let null = Fitted::of(
            &panel,
            Frequencies::TrueByGroup,
            clean.inbreeding,
            clean.error,
        );

        let mut spiked = clean.clone();
        spiked.contamination[0] = 0.03;
        let panel = spiked.draw(909);
        let fitted = Fitted::of(
            &panel,
            Frequencies::TrueByGroup,
            spiked.inbreeding,
            spiked.error,
        );
        println!(
            "  {loci:>7} {:>12} |          {:>8.4} {:>6.4} |             {:>8.4} {:>13.4}   ({:.0} s)",
            panel_segregating,
            null.mean(),
            null.max(),
            fitted.alphas[0],
            fitted.alphas[1..].iter().copied().fold(0.0_f64, f64::max),
            started.elapsed().as_secs_f64(),
        );
    }
}

/// **The case the owner raised**: subpopulations of wildly different size, and the
/// contaminated sample sitting in the smallest of them.
///
/// A decomposition gives its leading axes to the largest groups, so the two samples off on
/// their own are the ones whose ancestry a truncated set of axes can miss — and by §4's
/// measurement, a sample whose frequency falls back towards the panel average is a sample
/// whose contamination is underestimated.
fn unbalanced(depth: u32, loci: usize, fst: f64) {
    let sizes = vec![40_usize, 5, 3, 2];
    let samples: usize = sizes.iter().sum();
    println!(
        "\nUnbalanced subpopulations {sizes:?} — {samples} samples at F_st {fst:.2}, {depth} reads \
         a site, {loci} loci.\n  The contaminated sample is the last one, which sits in the group \
         of {}.",
        sizes[sizes.len() - 1]
    );
    println!("\n  true alpha |{}", arm_header());
    for alpha in [0.03_f64, 0.10] {
        let mut truth = Truth::clean(samples, depth, loci, sizes.len(), fst);
        truth.group_sizes = Some(sizes.clone());
        *truth.contamination.last_mut().expect("a sample") = alpha;
        let panel = truth.draw(4242);
        let spiked = samples - 1;
        let mut cells = Vec::new();
        for which in arms() {
            let fitted = Fitted::of(&panel, which, truth.inbreeding, truth.error);
            let others = fitted.alphas[..spiked]
                .iter()
                .copied()
                .fold(0.0_f64, f64::max);
            cells.push(format!("{:>10.4} {:>10.4}", fitted.alphas[spiked], others));
        }
        println!("  {alpha:>10.3} | {}", cells.join(" |"));
    }

    println!("\n  And the same panel with nobody contaminated, by group:");
    let mut truth = Truth::clean(samples, depth, loci, sizes.len(), fst);
    truth.group_sizes = Some(sizes.clone());
    let panel = truth.draw(4242);
    let group_of = truth.group_of();
    for which in arms() {
        let fitted = Fitted::of(&panel, which, truth.inbreeding, truth.error);
        let mut by_group = vec![0.0_f64; sizes.len()];
        for (sample, alpha) in fitted.alphas.iter().enumerate() {
            by_group[group_of[sample]] = by_group[group_of[sample]].max(*alpha);
        }
        let cells: Vec<String> = by_group.iter().map(|a| format!("{a:>8.4}")).collect();
        println!(
            "    {:>14}  worst in each group: {}",
            arm_name(which),
            cells.join(" ")
        );
    }

    println!(
        "\n  How much of its own fitted frequency each group supplies (its leverage).\n           A fair share would be {:.3}; 1.000 means the line at that sample is that sample.",
        5.0 / samples as f64
    );
    let (_, _, coordinates) = dosages_and_coordinates(&panel, truth.error, 4);
    let leverage = coordinate_leverage(&coordinates, 4);
    let mut worst = vec![0.0_f64; sizes.len()];
    for (sample, value) in leverage.iter().enumerate() {
        worst[group_of[sample]] = worst[group_of[sample]].max(*value);
    }
    let cells: Vec<String> = worst.iter().map(|h| format!("{h:>8.3}")).collect();
    println!(
        "    {:>14}  worst in each group: {}",
        "4 axes",
        cells.join(" ")
    );
}

fn main() {
    let mut args = std::env::args().skip(1);
    let mode = args.next().unwrap_or_else(|| "null".to_string());
    let samples: usize = args
        .next()
        .map_or(50, |a| a.parse().expect("a sample count"));
    let depth: u32 = args.next().map_or(3, |a| a.parse().expect("a read depth"));
    let loci: usize = args
        .next()
        .map_or(20_000, |a| a.parse().expect("a locus count"));
    let groups: usize = args.next().map_or(4, |a| a.parse().expect("a group count"));
    let fst: f64 = args.next().map_or(0.1, |a| a.parse().expect("an F_st"));

    match mode.as_str() {
        "spike" => spike(samples, depth, loci, groups, fst),
        "unbalanced" => unbalanced(depth, loci, fst),
        "budget" => budget_sweep(samples, depth, groups, fst),
        _ => null_sweep(samples, depth, loci, groups),
    }
}
