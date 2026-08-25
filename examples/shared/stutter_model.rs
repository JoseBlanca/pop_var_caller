//! **The STR stutter model and the search over it — shared, so one control guards both users.**
//!
//! `doc/devel/ng/spec/parameter_prepass_ssr.md` §3 and §4 settle a noise model for repeat tracts:
//! how often a read shows a length other than its allele's, which way it moves, and how far. This
//! module is that model and the maximisation over it, and **nothing here knows where its numbers
//! came from** — a caller hands it, per locus, how the reads split across whole-repeat offsets, and
//! gets back the three parameters that best explain the whole collection.
//!
//! ## Where the check that guards this code lives, and why it is not in the harness
//!
//! [`ng_str_stutter_harness.rs`](../ng_str_stutter_harness.rs) measures this estimator's bias
//! exactly: it enumerates **every possible locus shape**, weights each by its probability under a
//! known truth, and maximises — so what it reports has no sampling noise in it and "unbiased" is
//! decided rather than estimated. Its control reads exactly zero.
//!
//! **That control cannot be made to guard a real-data program, and the reason is structural.**
//! Enumerating the shape space is only affordable at the depths and bucket counts the harness uses;
//! a real locus at 30 reads across nine buckets has more shapes than there are loci in the genome.
//! So a program fitting real data has to build its per-genotype likelihoods **directly from the
//! shapes it observed**, which is different code from the harness's enumeration. Sharing the
//! optimiser would not close that gap, because the gap is in the component builder.
//!
//! **So the control travels with the program that needs it.**
//! [`ng_str_stutter_rate.rs`](../ng_str_stutter_rate.rs) carries its own exact-bias check, run
//! through its own component builder and this module's kernel and optimiser, and refuses to report
//! a real measurement until it passes. That is the lesson of a survey that ran for two days and
//! could not answer the question it was built for: the instrument was never asked to reproduce an
//! answer somebody already knew.
//!
//! What this module holds is the part both programs genuinely share — the slip kernel, the concave
//! climb over genotype frequencies, the multi-start search, and the algebraic gates. The harness
//! predates it and still carries its own copy of the kernel; reconciling the two is worth doing and
//! is not required for the checking story above.

#![allow(dead_code)]

// ---------------------------------------------------------------------------
// Arithmetic helpers
// ---------------------------------------------------------------------------

/// Lanczos approximation, g = 7 — accurate to about 15 digits over the small integers these
/// factorials use.
pub fn ln_gamma(x: f64) -> f64 {
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

pub fn ln_factorial(n: u32) -> f64 {
    ln_gamma(f64::from(n) + 1.0)
}

pub fn ln_sum_exp(values: &[f64]) -> f64 {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    max + values.iter().map(|v| (v - max).exp()).sum::<f64>().ln()
}

/// `x·ln y` with the convention `0·ln 0 = 0`, so an empty bucket against an impossible outcome
/// contributes nothing instead of a NaN.
pub fn x_ln_y(x: f64, y: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x * y.ln() }
}

/// The logit and its inverse — the scale the two share-shaped parameters are searched on, so a
/// search cannot walk out of `(0, 1)`.
pub fn logit(p: f64) -> f64 {
    (p / (1.0 - p)).ln()
}

pub fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

// ---------------------------------------------------------------------------
// The noise model
// ---------------------------------------------------------------------------

/// The furthest a read is allowed to slip. The kernel is renormalised over `1..=this`, so it is a
/// proper distribution; at the fall-offs measured on real data (7 to 12 reads in 100 take a second
/// step) the mass beyond eight steps is below 1e-8 of the slipped reads.
pub const MAX_SLIP_STEP: i32 = 8;

/// **How a read's repeat count moves away from the allele it came off.** Three numbers, each
/// answering a different question:
///
/// - `level` — how often a read slips at all, `P(|d| ≥ 1)`;
/// - `up_share` — of the reads that slipped, the share that *gained* repeats. Measured at 0.17 on
///   tomato dinucleotides (2,438 losses against 501 gains), so strongly asymmetric;
/// - `falloff` — of the reads that slipped, the chance of a second step given a first. Measured at
///   0.065 to 0.12, and taken to be the same in both directions (spec §3).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Slip {
    pub level: f64,
    pub up_share: f64,
    pub falloff: f64,
}

impl Slip {
    /// `P(a read shows exactly `d` whole repeats more than its allele)`.
    pub fn p(&self, d: i32) -> f64 {
        if d == 0 {
            return 1.0 - self.level;
        }
        let steps = d.abs();
        if steps > MAX_SLIP_STEP {
            return 0.0;
        }
        let direction = if d > 0 {
            self.up_share
        } else {
            1.0 - self.up_share
        };
        // Geometric over 1..=MAX_SLIP_STEP, renormalised so the truncation loses no mass.
        let tail = 1.0 - self.falloff.powi(MAX_SLIP_STEP);
        self.level * direction * (1.0 - self.falloff) * self.falloff.powi(steps - 1) / tail
    }
}

/// How the end buckets of a saturating offset range are scored.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EdgeScoring {
    /// The bucket's probability is the sum over every offset it absorbs. **The adopted rule**:
    /// measured exactly unbiased at every range tried (spec §4.1).
    Marginal,
    /// The bucket is scored as though every read in it sat exactly on the edge. Improper — the
    /// bucket probabilities do not sum to one — and kept only so the gate can reject it.
    PlugAtEdge,
    /// The same plug-in, then rescaled so the buckets sum to one. Proper, so only bias is at stake,
    /// and it costs +33% of the slippage level where slippage is heaviest.
    PlugAtEdgeRenormalised,
}

/// How offsets are recorded: a saturating range around the reference tract length.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Keying {
    /// Offsets are recorded over `−half_range ..= +half_range`, the ends saturating.
    pub half_range: i32,
    pub edges: EdgeScoring,
}

impl Keying {
    pub fn buckets(&self) -> usize {
        (2 * self.half_range + 1) as usize
    }

    /// Which bucket an offset falls in. The two end buckets are saturating.
    pub fn bucket_of(&self, offset: i32) -> usize {
        (offset.clamp(-self.half_range, self.half_range) + self.half_range) as usize
    }
}

/// **The distribution of one read's bucket, given the allele it came off.**
///
/// The read slips by `d` with the kernel above and lands at offset `allele + d`, which the key
/// clamps into a bucket. The end buckets therefore absorb every offset beyond the range, and how
/// that absorbed mass is scored is [`EdgeScoring`].
pub fn read_bucket_probs(slip: &Slip, allele_offset: i32, keying: &Keying) -> Vec<f64> {
    let mut out = vec![0.0; keying.buckets()];
    match keying.edges {
        EdgeScoring::Marginal => {
            for step in -MAX_SLIP_STEP..=MAX_SLIP_STEP {
                out[keying.bucket_of(allele_offset + step)] += slip.p(step);
            }
        }
        EdgeScoring::PlugAtEdge | EdgeScoring::PlugAtEdgeRenormalised => {
            for (bucket, slot) in out.iter_mut().enumerate() {
                let offset = bucket as i32 - keying.half_range;
                *slot = slip.p(offset - allele_offset);
            }
            if keying.edges == EdgeScoring::PlugAtEdgeRenormalised {
                let total: f64 = out.iter().sum();
                if total > 0.0 {
                    for p in &mut out {
                        *p /= total;
                    }
                }
            }
        }
    }
    out
}

/// The bucket distribution for a whole genotype: each read picks one of the two copies with equal
/// chance, then slips.
pub fn genotype_bucket_probs(
    slip: &Slip,
    allele_offsets: &[i32],
    genotype: (usize, usize),
    keying: &Keying,
) -> Vec<f64> {
    let first = read_bucket_probs(slip, allele_offsets[genotype.0], keying);
    let second = read_bucket_probs(slip, allele_offsets[genotype.1], keying);
    first
        .iter()
        .zip(&second)
        .map(|(a, b)| 0.5 * a + 0.5 * b)
        .collect()
}

/// Unordered pairs of allele indices — one genotype each.
pub fn genotype_pairs(alleles: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for i in 0..alleles {
        for j in i..alleles {
            out.push((i, j));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Fitting
// ---------------------------------------------------------------------------

/// **One stratum's per-genotype likelihoods, in log space and one flat block.**
///
/// `ln[entry * genotypes + g]` is the log-likelihood of that entry under that genotype.
///
/// **Logs rather than probabilities, and flat rather than nested, because this is the hot loop.**
/// The climb below evaluates `ln(freq) + ln(component)` once per entry per genotype per iteration;
/// holding the component already logged removes a `ln` from the innermost loop, and a single
/// allocation removes one per entry per candidate — of which a fit tries a few hundred.
pub struct LnComponents {
    pub ln: Vec<f64>,
    pub genotypes: usize,
}

impl LnComponents {
    pub fn entries(&self) -> usize {
        self.ln.len().checked_div(self.genotypes).unwrap_or(0)
    }
    pub fn row(&self, entry: usize) -> &[f64] {
        &self.ln[entry * self.genotypes..(entry + 1) * self.genotypes]
    }
}

/// **The genotype frequencies that best explain the table, at fixed slippage parameters.**
///
/// With the noise parameters held still the surface in the frequencies is concave, so this climb
/// cannot stop on a false summit — the one part of this fit with a proof behind it
/// (`parameter_prepass.md` §3.1). Expectation-maximization.
///
/// `weights[entry]` is how much of the data sits at that entry — a locus count on real data, an
/// exact probability where a truth is being recovered.
///
/// **`start` is a warm start and it is a speed device, not a modelling one.** The surface is
/// concave, so every start reaches the same summit; beginning from the previous candidate's answer
/// simply arrives in a few passes instead of hundreds, because neighbouring candidates have nearly
/// the same frequencies. It changes the answer only if the climb stops before converging, which is
/// what the iteration cap and the self-check between them are for.
pub fn climb_genotype_frequencies(
    component: &LnComponents,
    weights: &[f64],
    iterations: usize,
    start: Option<&[f64]>,
) -> Vec<f64> {
    let genotypes = component.genotypes;
    let total: f64 = weights.iter().sum();
    let mut freqs = match start {
        Some(warm) if warm.len() == genotypes => warm.to_vec(),
        _ => vec![1.0 / genotypes as f64; genotypes],
    };
    if total == 0.0 || genotypes == 0 {
        return freqs;
    }
    let mut ln_freqs = vec![0.0; genotypes];
    let mut terms = vec![0.0; genotypes];
    let mut next = vec![0.0; genotypes];
    for _ in 0..iterations {
        for g in 0..genotypes {
            ln_freqs[g] = freqs[g].max(1e-300).ln();
        }
        next.iter_mut().for_each(|v| *v = 0.0);
        for (entry, &weight) in weights.iter().take(component.entries()).enumerate() {
            if weight == 0.0 {
                continue;
            }
            let row = component.row(entry);
            for g in 0..genotypes {
                terms[g] = ln_freqs[g] + row[g];
            }
            let norm = ln_sum_exp(&terms);
            for g in 0..genotypes {
                next[g] += weight * (terms[g] - norm).exp();
            }
        }
        let mut moved: f64 = 0.0;
        for g in 0..genotypes {
            next[g] /= total;
            moved = moved.max((next[g] - freqs[g]).abs());
        }
        freqs.copy_from_slice(&next);
        if moved < 1e-12 {
            break;
        }
    }
    freqs
}

/// The marginal log-likelihood of the whole table at these frequencies: over entries, the weight
/// times the log of the genotype-averaged likelihood.
pub fn score_at(component: &LnComponents, weights: &[f64], freqs: &[f64]) -> f64 {
    let genotypes = component.genotypes;
    let mut ln_freqs = vec![0.0; genotypes];
    for g in 0..genotypes {
        ln_freqs[g] = freqs[g].max(1e-300).ln();
    }
    let mut terms = vec![0.0; genotypes];
    let mut score = 0.0;
    for (entry, &weight) in weights.iter().take(component.entries()).enumerate() {
        if weight == 0.0 {
            continue;
        }
        let row = component.row(entry);
        for g in 0..genotypes {
            terms[g] = ln_freqs[g] + row[g];
        }
        score += weight * ln_sum_exp(&terms);
    }
    score
}

/// How finely each axis is resolved, on its own transformed scale, and how many coordinate sweeps
/// the search takes.
///
/// **These are resolution knobs, not accuracy ones, and they are set by what a consumer can feel.**
/// A tolerance of 0.01 on the rate's log scale is 1% of the rate; the shared spec argues a few
/// percent between neighbouring values is already finer than a genotype likelihood can register,
/// since shifting a prior by that much moves a call by far less than one more read would. An
/// earlier version resolved to 1e-7 — seven digits on a quantity known to two — and spent about
/// four times the work to do it.
///
/// The self-check runs at exactly these values, so if they are ever too coarse the exact control
/// stops reading zero and says so.
#[derive(Clone, Copy, Debug)]
pub struct SearchPrecision {
    /// How finely an axis is resolved, on its own transformed scale. On the rate's log scale this
    /// is a relative resolution: 0.01 is 1% of the rate.
    pub tolerance: f64,
    pub max_axis_steps: usize,
    pub max_sweeps: usize,
}

impl SearchPrecision {
    /// **For a table with many entries — resolution matched to what a consumer can feel.** 1% of
    /// the rate, because the shared spec argues a few percent between neighbouring values is
    /// already finer than a genotype likelihood registers: shifting a prior by that much moves a
    /// call far less than one more read would.
    pub fn fast() -> Self {
        Self {
            tolerance: 0.01,
            max_axis_steps: 16,
            max_sweeps: 5,
        }
    }

    /// **For the exact control, where the point is sharpness rather than speed.** The control's
    /// job is to say the estimator is unbiased, and it can only say so to the resolution the
    /// search itself has — at `fast()` a genuine 1% bias and the search's own step are the same
    /// size, so the control would be reporting its own grid. The enumerated table it runs on is a
    /// couple of hundred entries, so this costs nothing worth counting.
    pub fn fine() -> Self {
        Self {
            tolerance: 1e-5,
            max_axis_steps: 40,
            max_sweeps: 8,
        }
    }
}

/// Maximise any score over the three slippage parameters, by coordinate-wise golden-section search
/// on the scale each parameter lives on: a rate on a log scale, two shares on a logit one, so a
/// search cannot walk out of the range a parameter has.
///
/// The search is continuous rather than on a ladder, because on the harness the bias being measured
/// is a property of the objective and not of the search resolution.
pub fn maximise_slip(
    score: impl Fn(&Slip) -> f64,
    start: Slip,
    precision: SearchPrecision,
) -> (Slip, f64) {
    let axes: [(f64, f64); 3] = [
        ((1e-5f64).ln(), (0.6f64).ln()),
        (logit(0.005), logit(0.995)),
        (logit(0.002), logit(0.95)),
    ];
    let inverse_phi = 0.5 * (5f64.sqrt() - 1.0);
    let mut current = start;
    let put = |slip: &mut Slip, axis: usize, x: f64| match axis {
        0 => slip.level = x.exp(),
        1 => slip.up_share = expit(x),
        _ => slip.falloff = expit(x),
    };

    for _sweep in 0..precision.max_sweeps {
        let mut moved: f64 = 0.0;
        for (axis, &(low, high)) in axes.iter().enumerate() {
            let (mut lo, mut hi) = (low, high);
            let mut c = hi - inverse_phi * (hi - lo);
            let mut d = lo + inverse_phi * (hi - lo);
            let evaluate = |x: f64, base: &Slip| {
                let mut trial = *base;
                put(&mut trial, axis, x);
                score(&trial)
            };
            let mut fc = evaluate(c, &current);
            let mut fd = evaluate(d, &current);
            for _ in 0..precision.max_axis_steps {
                if fc > fd {
                    hi = d;
                    d = c;
                    fd = fc;
                    c = hi - inverse_phi * (hi - lo);
                    fc = evaluate(c, &current);
                } else {
                    lo = c;
                    c = d;
                    fc = fd;
                    d = lo + inverse_phi * (hi - lo);
                    fd = evaluate(d, &current);
                }
                if (hi - lo) < precision.tolerance {
                    break;
                }
            }
            let before = match axis {
                0 => current.level.ln(),
                1 => logit(current.up_share),
                _ => logit(current.falloff),
            };
            let best = 0.5 * (lo + hi);
            moved = moved.max((best - before).abs());
            put(&mut current, axis, best);
        }
        if moved < precision.tolerance {
            break;
        }
    }
    let final_score = score(&current);
    (current, final_score)
}

/// **The starting points, spread over every axis the fit can stick on** (spec §4.2).
///
/// Starts that disagree only about the headline parameter are not a spread: on ng's inbreeding
/// estimator five starts sharing one guess at a nuisance axis returned a confident zero on a genome
/// 29% covered by runs. So each start here disagrees about the level, the direction and the
/// fall-off at once.
///
/// **The level is a multiplier on a starting estimate, not an absolute value**, because a stratum's
/// rate spans twenty-two-fold across repeat counts and a fixed ladder of absolute rates would begin
/// every stratum in the wrong place.
pub fn starting_points(level_estimate: f64) -> Vec<Slip> {
    let level = level_estimate.clamp(1e-4, 0.3);
    vec![
        Slip {
            level: (level * 3.0).min(0.5),
            up_share: 0.20,
            falloff: 0.03,
        },
        Slip {
            level: (level / 3.0).max(1e-5),
            up_share: 0.80,
            falloff: 0.40,
        },
        Slip {
            level,
            up_share: 0.50,
            falloff: 0.15,
        },
        Slip {
            level: (level * 0.3).max(1e-5),
            up_share: 0.35,
            falloff: 0.08,
        },
    ]
}

/// What a multi-start search returned: the winner, and enough beside it to tell a fitted number
/// from a stopped search.
#[derive(Clone, Debug)]
pub struct MultistartFit {
    pub slip: Slip,
    pub genotype_freqs: Vec<f64>,
    pub score: f64,
    /// The ratio between the highest and lowest level any start reached. **The column that says
    /// whether the fit found something**, and the one that is not sufficient on its own — a
    /// deterministic search returns the same point from every start wherever the objective is
    /// flat, so it is paired with a score comparison wherever a reference point exists.
    pub level_spread: f64,
    pub starts: usize,
}

/// Fit from every start and keep the best-scoring, with the spread reported beside it.
///
/// `score_and_freqs` is the caller's objective: given candidate slippage parameters it returns the
/// marginal log-likelihood and the genotype frequencies it climbed to. Keeping it a closure is what
/// lets the harness score an enumerated cell space and the real program score observed counts,
/// through one optimiser.
pub fn fit_from_starts(
    score_and_freqs: impl Fn(&Slip) -> (f64, Vec<f64>),
    starts: &[Slip],
    precision: SearchPrecision,
) -> MultistartFit {
    let mut reached: Vec<(Slip, f64)> = starts
        .iter()
        .map(|start| maximise_slip(|trial| score_and_freqs(trial).0, *start, precision))
        .collect();
    let levels: Vec<f64> = reached.iter().map(|(s, _)| s.level).collect();
    let high = levels.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let low = levels.iter().cloned().fold(f64::INFINITY, f64::min);
    reached.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
    let (slip, score) = reached[0];
    let (_, genotype_freqs) = score_and_freqs(&slip);
    MultistartFit {
        slip,
        genotype_freqs,
        score,
        level_spread: if low > 0.0 { high / low } else { f64::INFINITY },
        starts: starts.len(),
    }
}

// ---------------------------------------------------------------------------
// The gates — checks that need no fit, and reject a broken rule in one line
// ---------------------------------------------------------------------------

/// **Does the scoring rule sum to one over the buckets, for every genotype?** A rule that does not
/// is not the likelihood of anything, and no consistency result covers it. The un-rescaled plug-in
/// sums to 0.9488 at a half-range of 1, which is how it is rejected without fitting anything.
///
/// Returns the worst deviation from one, over the genotypes tried.
pub fn gate_sums_to_one(slip: &Slip, alleles: &[i32], keying: &Keying) -> f64 {
    let mut worst: f64 = 0.0;
    for genotype in genotype_pairs(alleles.len()) {
        let probs = genotype_bucket_probs(slip, alleles, genotype, keying);
        let total: f64 = probs.iter().sum();
        worst = worst.max((total - 1.0).abs());
    }
    worst
}

/// **With the slippage level at zero, does every locus's reads land on its own alleles?** A kernel
/// that puts mass anywhere else at a level of zero is describing movement that cannot happen.
///
/// Returns the largest probability found off the genotype's own allele buckets.
pub fn gate_silent_kernel(alleles: &[i32], keying: &Keying) -> f64 {
    let silent = Slip {
        level: 0.0,
        up_share: 0.5,
        falloff: 0.1,
    };
    let mut worst: f64 = 0.0;
    for genotype in genotype_pairs(alleles.len()) {
        let probs = genotype_bucket_probs(&silent, alleles, genotype, keying);
        let own = [
            keying.bucket_of(alleles[genotype.0]),
            keying.bucket_of(alleles[genotype.1]),
        ];
        for (bucket, &p) in probs.iter().enumerate() {
            if !own.contains(&bucket) {
                worst = worst.max(p);
            }
        }
    }
    worst
}

/// **Is every bucket probability a probability?** Negative or greater-than-one mass is the
/// arithmetic equivalent of charging a bucket a negative number of reads.
///
/// Returns the worst violation, zero when there is none.
pub fn gate_probabilities_in_range(slip: &Slip, alleles: &[i32], keying: &Keying) -> f64 {
    let mut worst: f64 = 0.0;
    for genotype in genotype_pairs(alleles.len()) {
        for p in genotype_bucket_probs(slip, alleles, genotype, keying) {
            if p < 0.0 {
                worst = worst.max(-p);
            } else if p > 1.0 {
                worst = worst.max(p - 1.0);
            }
        }
    }
    worst
}
