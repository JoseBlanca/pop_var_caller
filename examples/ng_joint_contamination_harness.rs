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
        let group_of: Vec<usize> = (0..self.samples)
            .map(|i| i * self.groups / self.samples)
            .collect();
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
}

/// The allele frequency each sample's genotype is scored against, per locus.
///
/// **Estimated from the reads, not supplied**, because that is what the fit has: the moment
/// estimator inverts the error rate out of the alternative-read share, which is what a
/// spectrum fitted over the same loci converges to.
fn estimate_frequencies(panel: &Panel, which: Frequencies, error: f64) -> Vec<Vec<f64>> {
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
        out.push(match which {
            Frequencies::Pooled => vec![invert(all_alt, all_depth); panel.groups.max(1)],
            Frequencies::ByGroup => per_group
                .iter()
                .map(|(alt, depth)| invert(*alt, *depth))
                .collect(),
            Frequencies::TrueByGroup => panel.true_frequencies[out.len()]
                .iter()
                .map(|f| f.clamp(1e-4, 1.0 - 1e-4))
                .collect(),
        });
    }
    out
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
    let group = panel.group_of[sample];
    let score = |alpha: f64| -> f64 {
        panel
            .data
            .iter()
            .zip(frequencies)
            .map(|(row, per_group)| {
                ln_locus(row[sample], per_group[group], inbreeding, error, alpha)
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

/// **The measurement that matters**: a clean panel, swept over how diverged its
/// subpopulations are.
fn null_sweep(samples: usize, depth: u32, loci: usize, groups: usize) {
    println!(
        "\nA clean panel — every sample's true contamination is zero.\n\
         {samples} samples in {groups} subpopulations, {depth} reads a site, {loci} loci."
    );
    println!(
        "\n  F_st |      pooled: mean     max  flagged |    by-group: mean     max  flagged \
         | true group: mean     max  flagged"
    );
    for fst in [0.0_f64, 0.02, 0.05, 0.10, 0.20] {
        let truth = Truth::clean(samples, depth, loci, groups, fst);
        let panel = truth.draw(4242);
        let mut cells = Vec::new();
        for which in [
            Frequencies::Pooled,
            Frequencies::ByGroup,
            Frequencies::TrueByGroup,
        ] {
            let fitted = Fitted::of(&panel, which, truth.inbreeding, truth.error);
            cells.push(format!(
                "{:>10.4} {:>7.4} {:>8}",
                fitted.mean(),
                fitted.max(),
                format!("{}/{}", fitted.flagged(0.01), samples)
            ));
        }
        println!("  {fst:>4.2} | {}", cells.join(" | "));
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
    println!(
        "\n  true alpha |    pooled: fitted  others' max |  by-group: fitted  others' max \
         | true group: fitted  others' max"
    );
    for alpha in [0.01_f64, 0.03, 0.10] {
        let mut truth = Truth::clean(samples, depth, loci, groups, fst);
        truth.contamination[0] = alpha;
        let panel = truth.draw(4242);
        assert_eq!(
            panel.truth[0], alpha,
            "the spiked sample is the one that was spiked"
        );
        let mut cells = Vec::new();
        for which in [
            Frequencies::Pooled,
            Frequencies::ByGroup,
            Frequencies::TrueByGroup,
        ] {
            let fitted = Fitted::of(&panel, which, truth.inbreeding, truth.error);
            cells.push(format!(
                "{:>16.4} {:>12.4}",
                fitted.alphas[0],
                fitted.alphas[1..].iter().copied().fold(0.0_f64, f64::max)
            ));
        }
        println!("  {alpha:>10.3} | {}", cells.join(" | "));
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
        "budget" => budget_sweep(samples, depth, groups, fst),
        _ => null_sweep(samples, depth, loci, groups),
    }
}
