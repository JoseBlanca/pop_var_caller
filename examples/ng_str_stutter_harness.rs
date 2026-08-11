//! **Is the STR path's stutter accumulator able to give an unbiased answer?** — the
//! harness that decides, rather than the argument that asserts.
//!
//! ng step 4's STR path estimates, per (read group × stratum), how often a read shows a
//! different number of repeats than the allele it came off, which way it moves, how far,
//! and a per-base substitution rate. **Nothing downstream can check any of them**, so a
//! wrong answer here is a plausible number nobody notices. This program measures the
//! errors that survive an infinite genome.
//!
//! ## The method — bias without simulation
//!
//! The estimator maximises a sum over cells. Replace each cell's observed count with the
//! cell's *exact probability under a known truth* and the sum becomes the objective the
//! estimator climbs with an infinite genome:
//!
//! ```text
//! Q(θ)  =  Σ over cells   P_true(cell) · ln L_model(cell ; θ)
//! θ*    =  argmax Q(θ)              what the estimate converges to
//! bias  =  θ* − θ_true              a fixed number, with no sampling error in it
//! ```
//!
//! `θ*` is the value a misspecified maximum-likelihood fit converges to — the parameter
//! whose model sits closest in Kullback–Leibler divergence to the truth (White 1982). No
//! draws, no seeds, no repeats: `bias = 0` and `bias ≠ 0` are decided rather than
//! estimated, and `bias = 0` is exactly the statement that the estimator is consistent.
//! Same shape as [`ng_multilib_key_harness`](ng_multilib_key_harness.rs), a different
//! accumulator.
//!
//! ## What is varied
//!
//! One **truth** — a slippage kernel, an allele-length spectrum, a depth distribution —
//! and several **keys**, each being one answer to a design question the spec leaves open:
//!
//! - **What the cell is.** One tally over every read in the stratum, which is what
//!   `parameter_prepass_ssr.md` §4.1 specifies; or one entry per locus, holding that
//!   locus's reads together. The genotype belongs to a locus, so this decides whether it
//!   can be summed over at all.
//! - **Where the offsets are measured from** — the reference tract length, or each
//!   locus's own modal observed length (§4.1's `OPEN`).
//! - **How a saturated end bucket is scored** — as the sum over everything it absorbs, or
//!   by plugging in its edge offset.
//!
//! ```text
//! cargo run --release --example ng_str_stutter_harness
//! cargo run --release --example ng_str_stutter_harness -- --only=origin
//! ```
//!
//! Sections: `gates`, `identify`, `origin`, `saturation`, `shape`, `strata`, `composition`.

use std::collections::HashMap;
use std::fmt::Write as _;

// ---------------------------------------------------------------------------
// Numerical helpers
// ---------------------------------------------------------------------------

/// Lanczos approximation, g = 7 — accurate to about 15 digits over the small integers
/// these factorials use.
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

fn ln_sum_exp(values: &[f64]) -> f64 {
    let max = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    max + values.iter().map(|v| (v - max).exp()).sum::<f64>().ln()
}

/// `x·ln y` with the convention `0·ln 0 = 0`, so an empty bucket against an impossible
/// outcome contributes nothing instead of a NaN.
fn x_ln_y(x: f64, y: f64) -> f64 {
    if x == 0.0 { 0.0 } else { x * y.ln() }
}

/// The logit and its inverse — the scale the two share-shaped parameters are searched on,
/// so a search cannot walk out of `(0, 1)`.
fn logit(p: f64) -> f64 {
    (p / (1.0 - p)).ln()
}
fn expit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

// ---------------------------------------------------------------------------
// The noise model
// ---------------------------------------------------------------------------

/// The furthest a read is allowed to slip. The kernel is renormalised over `1..=this`, so
/// it is a proper distribution; at the fall-offs measured on real data (7 to 12 reads in
/// 100 take a second step) the mass beyond eight steps is below 1e-8 of the slipped reads.
const MAX_SLIP_STEP: i32 = 8;

/// **How a read's repeat count moves away from the allele it came off.** Three numbers,
/// and each answers a different question:
///
/// - `level` — how often a read slips at all, `P(|d| ≥ 1)`;
/// - `up_share` — of the reads that slipped, the share that *gained* repeats. Measured at
///   0.17 on tomato dinucleotides (2,438 losses against 501 gains), so strongly asymmetric;
/// - `falloff` — of the reads that slipped, the chance of a second step given a first.
///   Measured at 0.065 to 0.12, and taken to be the same in both directions
///   (`parameter_prepass_ssr.md` §3).
#[derive(Clone, Copy, Debug, PartialEq)]
struct Slip {
    level: f64,
    up_share: f64,
    falloff: f64,
}

impl Slip {
    /// `P(a read shows exactly `d` whole repeats more than its allele)`.
    fn p(&self, d: i32) -> f64 {
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

/// Where a locus's offsets are measured from.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Origin {
    /// The reference tract length — a property of the locus, identical in every sample.
    Reference,
    /// The locus's own modal observed length, in this sample's reads. Data-dependent.
    LocusMode(TieRule),
}

/// What the modal origin does when two lengths tie for most common — which at three reads
/// a site is not a corner case but the ordinary case.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum TieRule {
    /// The shortest of the tied lengths.
    Shortest,
    /// The tied length closest to the reference.
    NearestReference,
}

/// How the end buckets of a saturating offset range are scored.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum EdgeScoring {
    /// The bucket's probability is the sum over every offset it absorbs.
    Marginal,
    /// The bucket is scored as though every read in it sat exactly on the edge. Improper:
    /// the resulting bucket probabilities do not sum to one.
    PlugAtEdge,
    /// The same plug-in, then rescaled so the buckets sum to one. A proper likelihood, so
    /// only bias is at stake.
    PlugAtEdgeRenormalised,
}

/// What the accumulator keeps about a locus.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Keying {
    /// Offsets are recorded over `−half_range ..= +half_range`, ends saturating.
    half_range: i32,
    origin: Origin,
    edges: EdgeScoring,
}

impl Keying {
    fn buckets(&self) -> usize {
        (2 * self.half_range + 1) as usize
    }
    /// Which bucket an offset falls in.
    fn bucket_of(&self, offset: i32) -> usize {
        (offset.clamp(-self.half_range, self.half_range) + self.half_range) as usize
    }
}

/// Everything that defines one truth: the slippage, the allele-length spectrum the loci
/// are drawn from, how correlated a locus's two alleles are, and the depth.
#[derive(Clone, Debug)]
struct World {
    /// What the world is called in the report.
    #[allow(dead_code)]
    name: String,
    slip: Slip,
    /// The allele lengths present, as whole-repeat offsets from the reference tract length.
    allele_offsets: Vec<i32>,
    /// How common each of them is. Sums to one.
    allele_probs: Vec<f64>,
    /// How much more often the two alleles of a locus are the same one than chance —
    /// the inbreeding coefficient, which ng's generic path fits per sample. Zero is
    /// Hardy-Weinberg.
    inbreeding: f64,
    mean_depth: f64,
}

impl World {
    /// The genotypes: unordered pairs of allele indices.
    fn genotypes(&self) -> Vec<(usize, usize)> {
        let n = self.allele_offsets.len();
        let mut out = Vec::new();
        for i in 0..n {
            for j in i..n {
                out.push((i, j));
            }
        }
        out
    }

    /// How common each genotype is: `F·p_a + (1−F)·p_a²` for a homozygote, `(1−F)·2p_a·p_b`
    /// otherwise.
    fn genotype_probs(&self) -> Vec<f64> {
        let f = self.inbreeding;
        self.genotypes()
            .into_iter()
            .map(|(i, j)| {
                let (pi, pj) = (self.allele_probs[i], self.allele_probs[j]);
                if i == j {
                    f * pi + (1.0 - f) * pi * pi
                } else {
                    (1.0 - f) * 2.0 * pi * pj
                }
            })
            .collect()
    }

    /// How often a locus's two alleles differ — the quantity the fit has to keep apart
    /// from slippage.
    fn heterozygosity(&self) -> f64 {
        self.genotypes()
            .into_iter()
            .zip(self.genotype_probs())
            .filter(|((i, j), _)| i != j)
            .map(|(_, p)| p)
            .sum()
    }

    /// Depth per locus: Poisson at `mean_depth`, conditioned on at least one read and
    /// truncated at `max_depth`. Returns `probs[n]` for `n` in `0..=max_depth`.
    fn depth_probs(&self, max_depth: u32) -> Vec<f64> {
        let lambda = self.mean_depth;
        let mut probs = vec![0.0; (max_depth + 1) as usize];
        for (n, slot) in probs.iter_mut().enumerate().skip(1) {
            let n = n as u32;
            *slot = (-lambda + f64::from(n) * lambda.ln() - ln_factorial(n)).exp();
        }
        let total: f64 = probs.iter().sum();
        for p in &mut probs {
            *p /= total;
        }
        probs
    }

    /// The mean depth the truncated distribution actually realises. **The truncated
    /// Poisson is the world**, not an approximation of an untruncated one, so this is what
    /// a report should label a row with — the nominal `mean_depth` is only how the world
    /// was specified.
    fn realised_mean_depth(&self, max_depth: u32) -> f64 {
        self.depth_probs(max_depth)
            .iter()
            .enumerate()
            .map(|(n, p)| n as f64 * p)
            .sum()
    }
}

/// **The distribution of one read's bucket, given the allele it came off.**
///
/// The read slips by `d` with the kernel above and lands at offset `allele + d`, which the
/// key clamps into a bucket. The end buckets therefore absorb every offset beyond the
/// range, and how that absorbed mass is scored is `EdgeScoring`.
fn read_bucket_probs(slip: &Slip, allele_offset: i32, keying: &Keying) -> Vec<f64> {
    let mut out = vec![0.0; keying.buckets()];
    match keying.edges {
        EdgeScoring::Marginal => {
            for step in -MAX_SLIP_STEP..=MAX_SLIP_STEP {
                out[keying.bucket_of(allele_offset + step)] += slip.p(step);
            }
        }
        EdgeScoring::PlugAtEdge | EdgeScoring::PlugAtEdgeRenormalised => {
            // Every bucket is scored at one offset — its own — so the mass the end buckets
            // absorb is simply not counted.
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

/// The bucket distribution for a whole genotype: each read picks one of the two copies
/// with equal chance, then slips.
fn genotype_bucket_probs(
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

// ---------------------------------------------------------------------------
// Cells
// ---------------------------------------------------------------------------

/// **One entry of the accumulator, holding a whole locus.** The reads of one locus share a
/// genotype, so keeping them together is what lets the genotype be summed over; a key that
/// pools reads across loci has thrown that away (see `question_identification`).
#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
struct LocusCell {
    depth: u32,
    counts: Vec<u16>,
}

/// Every way `total` reads split across `buckets` buckets.
fn for_each_split(total: u32, buckets: usize, visit: &mut impl FnMut(&[u16])) {
    fn recurse(
        index: usize,
        left: u32,
        buckets: usize,
        buffer: &mut Vec<u16>,
        visit: &mut impl FnMut(&[u16]),
    ) {
        if index + 1 == buckets {
            buffer[index] = left as u16;
            visit(buffer);
            return;
        }
        for taken in 0..=left {
            buffer[index] = taken as u16;
            recurse(index + 1, left - taken, buckets, buffer, visit);
        }
    }
    let mut buffer = vec![0u16; buckets];
    recurse(0, total, buckets, &mut buffer, visit);
}

/// Re-key a locus's bucket counts onto its own modal bucket.
fn recentre_on_mode(counts: &[u16], keying: &Keying, tie: TieRule) -> Vec<u16> {
    let mut best = 0usize;
    for bucket in 1..counts.len() {
        let better = counts[bucket] > counts[best];
        let tied = counts[bucket] == counts[best];
        let prefer_on_tie = match tie {
            TieRule::Shortest => false, // the lowest index already wins
            TieRule::NearestReference => {
                let here = (bucket as i32 - keying.half_range).abs();
                let there = (best as i32 - keying.half_range).abs();
                here < there
            }
        };
        if better || (tied && prefer_on_tie) {
            best = bucket;
        }
    }
    let mut out = vec![0u16; counts.len()];
    for (bucket, &count) in counts.iter().enumerate() {
        let shifted = bucket as i32 - best as i32;
        out[keying.bucket_of(shifted)] += count;
    }
    out
}

/// **The accumulator's cells, with each genotype's probability of producing each one.**
///
/// `component[cell][genotype]` is `P(this locus's depth) · P(these bucket counts | depth,
/// genotype)`. It depends on the slippage parameters and on the key, and **not** on how
/// common each genotype is — which is what lets the genotype frequencies be climbed
/// separately at fixed noise parameters, the one part of the surface with a proof behind
/// it.
struct CellTable {
    cells: Vec<LocusCell>,
    index: HashMap<LocusCell, usize>,
    component: Vec<Vec<f64>>,
    /// What the components sum to over the whole cell space, per genotype — one for a
    /// proper likelihood.
    mass_by_genotype: Vec<f64>,
}

fn build_cells(
    slip: &Slip,
    allele_offsets: &[i32],
    genotypes: &[(usize, usize)],
    depth_probs: &[f64],
    keying: &Keying,
) -> CellTable {
    let buckets = keying.buckets();
    let key_origin = keying.origin;
    let mut cells: Vec<LocusCell> = Vec::new();
    let mut index: HashMap<LocusCell, usize> = HashMap::new();
    let mut component: Vec<Vec<f64>> = Vec::new();

    let bucket_probs: Vec<Vec<f64>> = genotypes
        .iter()
        .map(|&g| genotype_bucket_probs(slip, allele_offsets, g, keying))
        .collect();

    for (depth, &p_depth) in depth_probs.iter().enumerate() {
        if p_depth <= 0.0 {
            continue;
        }
        let depth = depth as u32;
        for_each_split(depth, buckets, &mut |counts| {
            let keyed = match key_origin {
                Origin::Reference => counts.to_vec(),
                Origin::LocusMode(tie) => recentre_on_mode(counts, keying, tie),
            };
            let cell = LocusCell {
                depth,
                counts: keyed,
            };
            let slot = *index.entry(cell.clone()).or_insert_with(|| {
                cells.push(cell.clone());
                component.push(vec![0.0; genotypes.len()]);
                cells.len() - 1
            });
            let mut ln_coefficient = ln_factorial(depth);
            for &count in counts {
                ln_coefficient -= ln_factorial(u32::from(count));
            }
            for (g, probs) in bucket_probs.iter().enumerate() {
                let mut ln_p = ln_coefficient;
                for (bucket, &count) in counts.iter().enumerate() {
                    ln_p += x_ln_y(f64::from(count), probs[bucket]);
                }
                component[slot][g] += p_depth * ln_p.exp();
            }
        });
    }

    let mut mass_by_genotype = vec![0.0; genotypes.len()];
    for row in &component {
        for (g, &value) in row.iter().enumerate() {
            mass_by_genotype[g] += value;
        }
    }
    CellTable {
        cells,
        index,
        component,
        mass_by_genotype,
    }
}

/// The truth's weight on each cell: how often a locus of that shape occurs.
fn truth_mass(table: &CellTable, genotype_probs: &[f64]) -> Vec<f64> {
    table
        .component
        .iter()
        .map(|row| {
            row.iter()
                .zip(genotype_probs)
                .map(|(c, p)| c * p)
                .sum::<f64>()
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Fitting
// ---------------------------------------------------------------------------

/// **The genotype frequencies that best explain the table, at fixed slippage parameters.**
///
/// With the noise parameters held still the surface in the frequencies is concave, so this
/// climb cannot stop on a false summit — the one part of this fit with a proof behind it
/// (`parameter_prepass.md` §3.1). Expectation-maximization, from a uniform start.
fn climb_genotype_frequencies(
    component: &[Vec<f64>],
    weights: &[f64],
    genotypes: usize,
) -> Vec<f64> {
    climb_genotype_frequencies_for(component, weights, genotypes, CLIMB_ITERATIONS.get())
}

thread_local! {
    /// How many expectation-maximization passes the inner climb takes. **A knob because it
    /// turned out to matter**: the climb converges at a rate set by how much the genotypes
    /// overlap, so a keying that makes a locus's shape less informative about its genotype
    /// makes it slow — and a climb that stops short leaves the outer search a tilted
    /// objective, which is indistinguishable from bias unless it is checked.
    /// `question_convergence` is that check.
    static CLIMB_ITERATIONS: std::cell::Cell<usize> = const { std::cell::Cell::new(200) };
}

fn climb_genotype_frequencies_for(
    component: &[Vec<f64>],
    weights: &[f64],
    genotypes: usize,
    iterations: usize,
) -> Vec<f64> {
    let total: f64 = weights.iter().sum();
    let mut freqs = vec![1.0 / genotypes as f64; genotypes];
    for _ in 0..iterations {
        let mut next = vec![0.0; genotypes];
        for (row, &weight) in component.iter().zip(weights) {
            if weight == 0.0 {
                continue;
            }
            let mut terms = vec![0.0; genotypes];
            for g in 0..genotypes {
                terms[g] = freqs[g].max(1e-300).ln() + row[g].max(1e-300).ln();
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
        freqs = next;
        if moved < 1e-12 {
            break;
        }
    }
    freqs
}

/// Everything held fixed while a candidate slippage is scored.
struct FitInputs<'a> {
    world: &'a World,
    /// The cells the truth put weight on, and how much.
    truth_cells: &'a [LocusCell],
    truth_weights: &'a [f64],
    /// How the model believes the accumulator was keyed. Deliberately separate from how it
    /// **was** keyed: a model that scores a mode-centred table as though it were
    /// reference-centred is the naive implementation, and the gap is what this measures.
    model_keying: Keying,
    depth_probs: &'a [f64],
    /// **The allele lengths the fit is allowed to place mass on**, which need not be the
    /// ones the world used. Wider than the recorded offset range is the ordinary case, and
    /// *narrower than the world's alleles* is the case `question_allele_support` measures:
    /// a locus whose allele the fit cannot represent has its reads explained the only way
    /// left, which is slippage. `None` means "the same alleles the world has", the
    /// well-specified default every other question uses.
    model_alleles: Option<&'a [i32]>,
}

impl FitInputs<'_> {
    fn model_alleles(&self) -> &[i32] {
        self.model_alleles.unwrap_or(&self.world.allele_offsets)
    }
}

/// Unordered pairs of allele indices — one genotype each.
fn genotype_pairs(alleles: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for i in 0..alleles {
        for j in i..alleles {
            out.push((i, j));
        }
    }
    out
}

/// The objective at one candidate slippage: climb the genotype frequencies, then score.
fn score_slip(inputs: &FitInputs, slip: &Slip) -> (f64, Vec<f64>) {
    let alleles = inputs.model_alleles();
    let genotypes = genotype_pairs(alleles.len());
    let model = build_cells(
        slip,
        alleles,
        &genotypes,
        inputs.depth_probs,
        &inputs.model_keying,
    );
    // Line the model's components up with the cells the truth actually produced.
    let mut component = Vec::with_capacity(inputs.truth_cells.len());
    for cell in inputs.truth_cells {
        match model.index.get(cell) {
            Some(&slot) => component.push(model.component[slot].clone()),
            None => component.push(vec![1e-300; genotypes.len()]),
        }
    }
    let freqs = climb_genotype_frequencies(&component, inputs.truth_weights, genotypes.len());
    let mut score = 0.0;
    for (row, &weight) in component.iter().zip(inputs.truth_weights) {
        if weight == 0.0 {
            continue;
        }
        let mut terms = vec![0.0; genotypes.len()];
        for g in 0..genotypes.len() {
            terms[g] = freqs[g].max(1e-300).ln() + row[g].max(1e-300).ln();
        }
        score += weight * ln_sum_exp(&terms);
    }
    (score, freqs)
}

/// One fit's answer.
#[derive(Clone, Debug)]
struct Fitted {
    slip: Slip,
    genotype_freqs: Vec<f64>,
    score: f64,
}

impl Fitted {
    /// How often the fitted frequencies say a locus carries two different alleles. Read off
    /// the **fit's** genotype list, which is the model's allele support and not necessarily
    /// the world's.
    fn heterozygosity(&self, alleles: usize) -> f64 {
        genotype_pairs(alleles)
            .into_iter()
            .zip(&self.genotype_freqs)
            .filter(|((i, j), _)| i != j)
            .map(|(_, &p)| p)
            .sum()
    }
}

/// Maximise any score over the three slippage parameters, by coordinate-wise golden-section
/// search on the scale each parameter lives on: a rate on a log scale, two shares on a
/// logit one, so a search cannot walk out of the range a parameter has.
///
/// The search is continuous rather than on a ladder, because the bias being measured is a
/// property of the objective and not of the search resolution.
fn maximise_slip(score: impl Fn(&Slip) -> f64, start: Slip) -> (Slip, f64) {
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

    for _sweep in 0..8 {
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
            for _ in 0..34 {
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
                if (hi - lo) < 1e-7 {
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
        if moved < 1e-7 {
            break;
        }
    }
    let final_score = score(&current);
    (current, final_score)
}

/// The per-locus fit: the same search, with the genotype frequencies climbed to their
/// optimum at every trial.
fn fit_slip(inputs: &FitInputs, start: Slip) -> Fitted {
    let (slip, score) = maximise_slip(|trial| score_slip(inputs, trial).0, start);
    let (_, freqs) = score_slip(inputs, &slip);
    Fitted {
        slip,
        genotype_freqs: freqs,
        score,
    }
}

/// **The starting points, spread over every axis the fit can stick on.**
///
/// Starts that disagree only about the headline parameter are not a spread: on ng's
/// inbreeding estimator five starts sharing one guess at a nuisance axis returned a
/// confident zero on a genome 29% covered by runs. So each start here disagrees about the
/// level, the direction and the fall-off at once, and none of them begins at the truth.
fn starting_points(truth: &Slip) -> Vec<Slip> {
    vec![
        Slip {
            level: truth.level * 3.0,
            up_share: 0.20,
            falloff: 0.03,
        },
        Slip {
            level: truth.level / 3.0,
            up_share: 0.80,
            falloff: 0.40,
        },
        Slip {
            level: 0.05,
            up_share: 0.50,
            falloff: 0.15,
        },
        Slip {
            level: 0.005,
            up_share: 0.35,
            falloff: 0.08,
        },
    ]
}

/// How far apart the starting points left the answer — the column that says whether a fit
/// found something or merely stopped.
struct Spread {
    /// The ratio between the highest and lowest slippage level returned.
    level: f64,
    /// The range of fall-offs returned, on the fall-off's own scale. **Kept and not
    /// reported, deliberately**: it read exactly zero on a fit whose fall-off was 0.02 out,
    /// because a deterministic search returns the same point from every start wherever the
    /// objective is nearly flat. `question_convergence`'s score-at-the-truth is the
    /// diagnostic that works; this one is here so nobody adds it back as if it did.
    #[allow(dead_code)]
    falloff: f64,
}

/// Fit from every start and keep the best-scoring, which is what the design does. The
/// spread across starts is reported beside it, because a zero from a search that never
/// looked elsewhere is not a measurement.
fn fit_from_spread(inputs: &FitInputs) -> (Fitted, Spread) {
    let starts = starting_points(&inputs.world.slip);
    let mut fits: Vec<Fitted> = starts
        .into_iter()
        .map(|start| fit_slip(inputs, start))
        .collect();
    let range = |values: Vec<f64>| {
        (
            values.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
            values.iter().cloned().fold(f64::INFINITY, f64::min),
        )
    };
    let (level_high, level_low) = range(fits.iter().map(|f| f.slip.level).collect());
    let (falloff_high, falloff_low) = range(fits.iter().map(|f| f.slip.falloff).collect());
    fits.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
    (
        fits.remove(0),
        Spread {
            level: level_high / level_low,
            falloff: falloff_high - falloff_low,
        },
    )
}

// ---------------------------------------------------------------------------
// The worlds
// ---------------------------------------------------------------------------

/// A stratum shaped like tomato's dinucleotides at six or more repeats: 2 reads in 100
/// slip, a read is about five times as likely to lose a repeat as gain one, and 9 in 100
/// of the slipped reads take a second step (`parameter_prepass_ssr.md` §3, §5).
fn measured_slip() -> Slip {
    Slip {
        level: 0.020,
        up_share: 0.17,
        falloff: 0.087,
    }
}

/// A locus population: allele lengths within one repeat of the reference, with the
/// heterozygosity set by how spread the spectrum is.
fn world(name: &str, slip: Slip, spread: f64, inbreeding: f64, mean_depth: f64) -> World {
    World {
        name: name.to_string(),
        slip,
        allele_offsets: vec![-1, 0, 1],
        allele_probs: vec![spread / 2.0, 1.0 - spread, spread / 2.0],
        inbreeding,
        mean_depth,
    }
}

// ---------------------------------------------------------------------------
// The checks that need no fit
// ---------------------------------------------------------------------------

fn question_gates(report: &mut String) {
    let _ = writeln!(report, "\n## Checks that need no fit\n");
    let _ = writeln!(
        report,
        "Three identities, one line each, run before anything is fitted. A rule that fails\n\
         any of them cannot be unbiased, so each rejects a design without a measurement.\n"
    );

    let truth = measured_slip();
    let probe = world("probe", truth, 0.30, 0.0, 3.0);
    let depth_probs = probe.depth_probs(10);
    let genotypes = probe.genotypes();

    let _ = writeln!(
        report,
        "**1. Does the scoring rule sum to one over the cell space?** Per genotype, since a\n\
         cell's probability is conditional on it. Reported as the worst genotype's total.\n"
    );
    let _ = writeln!(
        report,
        "{:<12} {:<26} {:>14} {:>8}",
        "range", "edge scoring", "worst Σ L", "verdict"
    );
    for half_range in [1, 2, 3] {
        for edges in [
            EdgeScoring::Marginal,
            EdgeScoring::PlugAtEdge,
            EdgeScoring::PlugAtEdgeRenormalised,
        ] {
            let keying = Keying {
                half_range,
                origin: Origin::Reference,
                edges,
            };
            let table = build_cells(
                &truth,
                &probe.allele_offsets,
                &genotypes,
                &depth_probs,
                &keying,
            );
            let worst = table
                .mass_by_genotype
                .iter()
                .map(|m| (m - 1.0).abs())
                .fold(0.0f64, f64::max);
            let reported = table
                .mass_by_genotype
                .iter()
                .cloned()
                .fold(1.0, |acc: f64, m| {
                    if (m - 1.0).abs() > (acc - 1.0).abs() {
                        m
                    } else {
                        acc
                    }
                });
            let _ = writeln!(
                report,
                "±{:<11} {:<26} {:>14.9} {:>8}",
                half_range,
                format!("{edges:?}"),
                reported,
                if worst < 1e-9 { "PASS" } else { "FAIL" }
            );
        }
    }

    let _ = writeln!(
        report,
        "\n**2. Is any bucket charged a negative number of reads?** The origin bucket holds\n\
         the depth less every other bucket, so a key whose buckets can exceed the depth\n\
         would charge one. Reported as the worst cell's origin count.\n"
    );
    let keying = Keying {
        half_range: 2,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };
    let table = build_cells(
        &truth,
        &probe.allele_offsets,
        &genotypes,
        &depth_probs,
        &keying,
    );
    let bad = table
        .cells
        .iter()
        .filter(|cell| u32::from(cell.counts.iter().sum::<u16>()) != cell.depth)
        .count();
    let _ = writeln!(
        report,
        "  every cell's buckets sum to its depth: {} cells, {bad} violations   {}",
        table.cells.len(),
        if bad == 0 { "PASS" } else { "FAIL" }
    );

    let _ = writeln!(
        report,
        "\n**3. Does it reduce exactly in the null case?** With the slippage level at zero a\n\
         locus's reads must all land on its own alleles, so a homozygous genotype puts all\n\
         its mass on the single cell that holds every read at that allele's offset.\n"
    );
    let silent = Slip {
        level: 0.0,
        ..truth
    };
    let table = build_cells(
        &silent,
        &probe.allele_offsets,
        &genotypes,
        &depth_probs,
        &keying,
    );
    let mut worst_leak: f64 = 0.0;
    for (slot, cell) in table.cells.iter().enumerate() {
        for (g, &(i, j)) in genotypes.iter().enumerate() {
            if i != j {
                continue;
            }
            let expected_bucket = keying.bucket_of(probe.allele_offsets[i]);
            let all_there = cell.counts[expected_bucket] as u32 == cell.depth;
            if !all_there {
                worst_leak = worst_leak.max(table.component[slot][g]);
            }
        }
    }
    let _ = writeln!(
        report,
        "  mass a silent kernel puts anywhere but the allele's own bucket: {worst_leak:.3e}   {}",
        if worst_leak < 1e-12 { "PASS" } else { "FAIL" }
    );

    let _ = writeln!(
        report,
        "\n**4. The control the whole method rests on.** Key the accumulator by locus with the\n\
         reference origin, generate and fit under the same key, and the bias must be exactly\n\
         zero. A number other than zero here is the harness's, not the estimator's.\n"
    );
    let keying = Keying {
        half_range: 2,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };
    let table = build_cells(
        &truth,
        &probe.allele_offsets,
        &genotypes,
        &depth_probs,
        &keying,
    );
    let weights = truth_mass(&table, &probe.genotype_probs());
    let inputs = FitInputs {
        world: &probe,
        truth_cells: &table.cells,
        truth_weights: &weights,
        model_keying: keying,
        depth_probs: &depth_probs,
        model_alleles: None,
    };
    let (best, spread) = fit_from_spread(&inputs);
    let _ = writeln!(
        report,
        "  level {:+.3}%   up-share {:+.4}   fall-off {:+.4}   heterozygosity {:+.3}%   \
         spread across starts {:.3}×",
        100.0 * (best.slip.level / truth.level - 1.0),
        best.slip.up_share - truth.up_share,
        best.slip.falloff - truth.falloff,
        100.0 * (best.heterozygosity(probe.allele_offsets.len()) / probe.heterozygosity() - 1.0),
        spread.level
    );
}

// ---------------------------------------------------------------------------
// Question 1 — can a per-read tally identify anything?
// ---------------------------------------------------------------------------

/// **The accumulator `parameter_prepass_ssr.md` §4.1 specifies pools reads across loci.**
/// One tally per (read group, period, repeat count), counting how many *reads* showed each
/// whole-repeat offset. This asks what such a tally can still say.
///
/// A read that pools with reads from other loci carries no genotype: it drew one allele
/// from the stratum's spectrum and then slipped. So the observed offset distribution is the
/// **convolution** of the allele spectrum with the slippage kernel, and separating the two
/// is a deconvolution with both halves unknown.
fn question_identification(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## Pooling reads across loci: what a per-read tally still identifies\n"
    );
    let truth = measured_slip();
    let _ = writeln!(
        report,
        "The tally holds how many **reads** landed at each whole-repeat offset, pooled over\n\
         every locus in the stratum. A read carries no genotype — it drew one allele and\n\
         slipped — so what the tally holds is the allele spectrum **convolved** with the\n\
         slippage kernel, and recovering the kernel from that means undoing a convolution\n\
         with both halves unknown.\n\n\
         `alleles` is the share of loci carrying an allele other than the reference length.\n\
         `spectrum` says whether the fit has to find the allele spectrum too, or was handed\n\
         it. `spread` is the ratio between the highest and lowest answer across four starting\n\
         points that disagree about all three parameters at once — the column that says\n\
         whether the fit found an answer or merely stopped. Truth: {:.1}% of reads slip.\n",
        100.0 * truth.level
    );
    let _ = writeln!(
        report,
        "{:<12} {:>9} {:<11} {:>12} {:>11} {:>11} {:>9}",
        "allele range", "alleles", "spectrum", "level", "up-share", "fall-off", "spread"
    );

    let keying = Keying {
        half_range: 3,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };

    for (allele_offsets, spread_of_alleles) in [
        (vec![0i32], 0.0f64),
        (vec![-1, 0, 1], 0.05),
        (vec![-1, 0, 1], 0.30),
        (vec![-2, -1, 0, 1, 2], 0.30),
        (vec![-3, -2, -1, 0, 1, 2, 3], 0.30),
    ] {
        for supplied in [true, false] {
            let alleles = allele_offsets.clone();
            let allele_probs: Vec<f64> = if alleles.len() == 1 {
                vec![1.0]
            } else {
                alleles
                    .iter()
                    .map(|&a| {
                        if a == 0 {
                            1.0 - spread_of_alleles
                        } else {
                            spread_of_alleles / (alleles.len() - 1) as f64
                        }
                    })
                    .collect()
            };
            // The per-read cell space: a read's bucket, with the allele spectrum as the
            // mixing weights. There is no genotype anywhere in it.
            let bucket_component = |slip: &Slip| -> Vec<Vec<f64>> {
                (0..keying.buckets())
                    .map(|bucket| {
                        alleles
                            .iter()
                            .map(|&a| read_bucket_probs(slip, a, &keying)[bucket])
                            .collect()
                    })
                    .collect()
            };
            let truth_component = bucket_component(&truth);
            let weights: Vec<f64> = truth_component
                .iter()
                .map(|row| {
                    row.iter()
                        .zip(&allele_probs)
                        .map(|(c, p)| c * p)
                        .sum::<f64>()
                })
                .collect();

            let score = |slip: &Slip| -> f64 {
                let component = bucket_component(slip);
                let freqs = if supplied {
                    allele_probs.clone()
                } else {
                    climb_genotype_frequencies(&component, &weights, alleles.len())
                };
                let mut total = 0.0;
                for (row, &weight) in component.iter().zip(&weights) {
                    let mut terms = vec![0.0; alleles.len()];
                    for a in 0..alleles.len() {
                        terms[a] = freqs[a].max(1e-300).ln() + row[a].max(1e-300).ln();
                    }
                    total += weight * ln_sum_exp(&terms);
                }
                total
            };

            let mut answers: Vec<(Slip, f64)> = starting_points(&truth)
                .into_iter()
                .map(|start| maximise_slip(score, start))
                .collect();
            let lowest = answers
                .iter()
                .map(|a| a.0.level)
                .fold(f64::INFINITY, f64::min);
            let highest = answers
                .iter()
                .map(|a| a.0.level)
                .fold(f64::NEG_INFINITY, f64::max);
            answers.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            let best = answers[0].0;
            let _ = writeln!(
                report,
                "{:<12} {:>8.0}% {:<11} {:>+11.1}% {:>+11.4} {:>+11.4} {:>8.1}×",
                format!("{}..{}", alleles[0], alleles[alleles.len() - 1]),
                100.0 * spread_of_alleles,
                if supplied { "supplied" } else { "fitted" },
                100.0 * (best.level / truth.level - 1.0),
                best.up_share - truth.up_share,
                best.falloff - truth.falloff,
                highest / lowest
            );
        }
    }

    let _ = writeln!(
        report,
        "\nThe first two rows are the control: with every locus at the reference length there\n\
         is nothing for slippage to be confused with, and the tally must return the truth\n\
         whether the spectrum is handed over or not."
    );
}

// ---------------------------------------------------------------------------
// Question 2 — the offset origin
// ---------------------------------------------------------------------------

/// **What the offsets are measured from** — the reference tract length, which is a property
/// of the locus, or the locus's own modal observed length, which is a property of this
/// sample's reads.
///
/// The spec's claim is that centring on the mode is safe because the fit marginalises over
/// the genotype, so a heterozygous locus's second allele is explained by the genotype term
/// rather than charged to slippage. Two ways the claim can fail, and both are measured
/// here: the model may score a mode-centred table as if it were reference-centred (the
/// naive implementation), or it may model the centring correctly and still lose.
fn question_origin(report: &mut String) {
    let _ = writeln!(report, "\n## Where the offsets are measured from\n");
    let _ = writeln!(
        report,
        "Every row generates from the same truth — {:.1}% of reads slip, {:.0} in 100 of\n\
         those gain rather than lose, {:.0} in 100 take a second step — and differs only in\n\
         how the accumulator was keyed and how the fit scored it. `hets` is the share of\n\
         loci carrying two different alleles. Bias is what an infinite genome returns.\n",
        100.0 * measured_slip().level,
        100.0 * measured_slip().up_share,
        100.0 * measured_slip().falloff
    );

    let truth = measured_slip();
    let reference = Keying {
        half_range: 2,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };
    let modal = Keying {
        half_range: 2,
        origin: Origin::LocusMode(TieRule::Shortest),
        edges: EdgeScoring::Marginal,
    };
    let modal_near = Keying {
        half_range: 2,
        origin: Origin::LocusMode(TieRule::NearestReference),
        edges: EdgeScoring::Marginal,
    };

    for &depth in &[3.0f64, 6.0, 12.0] {
        for &spread in &[0.05f64, 0.30] {
            let w = world("origin", truth, spread, 0.0, depth);
            let depth_probs = w.depth_probs(12);
            let genotypes = w.genotypes();
            let realised = w.realised_mean_depth(12);
            let _ = writeln!(
                report,
                "\n### {realised:.1} reads a locus, {:.0} loci in 100 heterozygous",
                100.0 * w.heterozygosity(),
            );
            let _ = writeln!(
                report,
                "{:<44} {:>11} {:>11} {:>11} {:>10} {:>9} {:>12}",
                "accumulator keyed / model scores it as",
                "level",
                "up-share",
                "fall-off",
                "hets",
                "level ×",
                "beats truth"
            );
            let arms: [(&str, Keying, Keying); 4] = [
                ("reference origin / reference", reference, reference),
                ("modal origin / reference (naive)", modal, reference),
                ("modal origin / modal (marginal)", modal, modal),
                (
                    "modal origin, ties to reference / modal",
                    modal_near,
                    modal_near,
                ),
            ];
            for (label, truth_keying, model_keying) in arms {
                let table = build_cells(
                    &truth,
                    &w.allele_offsets,
                    &genotypes,
                    &depth_probs,
                    &truth_keying,
                );
                let weights = truth_mass(&table, &w.genotype_probs());
                let inputs = FitInputs {
                    world: &w,
                    truth_cells: &table.cells,
                    truth_weights: &weights,
                    model_keying,
                    depth_probs: &depth_probs,
                    model_alleles: None,
                };
                let (best, spread) = fit_from_spread(&inputs);
                // **How much better than the truth the fitted answer scores**, in nats per
                // locus. A correctly specified model cannot beat the truth, so a value at
                // zero beside a non-zero bias means the objective is *flat* in that
                // direction — the answer is not identified rather than displaced. The
                // spread across starts cannot say so on its own: golden-section search is
                // deterministic, so on an exactly flat direction every start returns the
                // same point and the spread reads zero.
                let beats_truth = best.score - score_slip(&inputs, &truth).0;
                let _ = writeln!(
                    report,
                    "{:<44} {:>+10.1}% {:>+11.4} {:>+11.4} {:>+9.1}% {:>8.2}× {:>12.2e}",
                    label,
                    100.0 * (best.slip.level / truth.level - 1.0),
                    best.slip.up_share - truth.up_share,
                    best.slip.falloff - truth.falloff,
                    100.0
                        * (best.heterozygosity(w.allele_offsets.len()) / w.heterozygosity() - 1.0),
                    spread.level,
                    beats_truth
                );
            }
        }
    }
}

/// **Is the mode-centred fit's residual a flat direction, or a climb that stopped short?**
///
/// The mode-centred table scored by summing over its own centring returns the slippage
/// level and the direction split right, and the fall-off 0.02 to 0.04 high. A correctly
/// specified model cannot be biased, so one of two things is happening: the objective is
/// flat along the fall-off, or the inner climb over the genotype frequencies has not
/// converged and the outer search is walking a tilted surface.
///
/// **The spread across starting points cannot tell them apart** — golden-section search is
/// deterministic, so on a flat direction every start returns the same point and the spread
/// reads zero either way. What tells them apart is the score at the truth, and how it moves
/// as the climb is given more iterations.
fn question_convergence(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## The mode-centred residual: flat direction, or a climb that stopped short?\n"
    );
    let _ = writeln!(
        report,
        "One world, one arm — the accumulator keyed on each locus's modal length and the fit\n\
         scoring it by summing over the centring, which is the arm that returns the level right\n\
         and the fall-off high. `beats truth` is how much better the fitted answer scores than\n\
         the truth, in nats per locus: **a correctly specified model cannot exceed zero**, so a\n\
         positive value is the climb's, and it should fall as the climb is given more passes.\n"
    );
    let _ = writeln!(
        report,
        "{:<12} {:>14} {:>16} {:>16} {:>14}",
        "origin", "climb passes", "score at truth", "score at fitted", "fitted − truth"
    );

    let truth = measured_slip();
    // The fitted fall-off that came back at 200 passes on this world, and the point whose
    // score is being compared against the truth's.
    let fitted = Slip {
        falloff: truth.falloff + 0.0231,
        ..truth
    };
    let w = world("converge", truth, 0.30, 0.0, 3.0);
    let depth_probs = w.depth_probs(12);
    let genotypes = w.genotypes();

    for (label, keying) in [
        (
            "modal",
            Keying {
                half_range: 2,
                origin: Origin::LocusMode(TieRule::Shortest),
                edges: EdgeScoring::Marginal,
            },
        ),
        (
            "reference",
            Keying {
                half_range: 2,
                origin: Origin::Reference,
                edges: EdgeScoring::Marginal,
            },
        ),
    ] {
        let table = build_cells(&truth, &w.allele_offsets, &genotypes, &depth_probs, &keying);
        let weights = truth_mass(&table, &w.genotype_probs());
        let inputs = FitInputs {
            world: &w,
            truth_cells: &table.cells,
            truth_weights: &weights,
            model_keying: keying,
            depth_probs: &depth_probs,
            model_alleles: None,
        };
        for passes in [200usize, 2_000, 20_000, 200_000] {
            CLIMB_ITERATIONS.with(|c| c.set(passes));
            let at_truth = score_slip(&inputs, &truth).0;
            let at_fitted = score_slip(&inputs, &fitted).0;
            let _ = writeln!(
                report,
                "{label:<12} {passes:>14} {at_truth:>16.10} {at_fitted:>16.10} {:>14.2e}",
                at_fitted - at_truth
            );
        }
    }
    CLIMB_ITERATIONS.with(|c| c.set(200));
    let _ = writeln!(
        report,
        "\nTwo scores per row, not a fit: the objective at the truth and at the point the\n\
         200-pass search returned, with only the number of expectation-maximization passes\n\
         changing. **A correctly specified model puts its maximum at the truth**, so the last\n\
         column must be at or below zero once the climb has converged. The reference-origin\n\
         rows are the control — there the climb converges at 200 passes and the column is\n\
         already negative."
    );
}

/// **What a locus whose allele the fit cannot represent costs.**
///
/// The fit places mass on a bounded list of allele lengths. Real tract lengths have a long tail:
/// on HG002 at 300×, 88.9 loci in 100 sit exactly at the reference length but about one in 200
/// sits further than six repeats away (research note §6.8). A locus outside the list has its reads
/// explained the only way left — as slippage — so the question is what that costs the three
/// parameters, and it is a different question from the saturating end buckets, which measured
/// alleles outside the *recorded range* while still inside the fitted support.
///
/// The world's alleles reach ±4 here and the fit's support is narrowed a step at a time. The row
/// where the support covers the world is the control and must read exactly zero.
fn question_allele_support(report: &mut String) {
    let _ = writeln!(report, "\n## Alleles the fit is not allowed to place\n");
    let truth = measured_slip();
    let _ = writeln!(
        report,
        "Loci carry alleles out to ±4 repeats from the reference. `support` is how far the fit may\n\
         place one; `outside` is the share of **loci** whose allele it therefore cannot represent.\n\
         Truth: {:.1}% of reads slip, direction split {:.2}, fall-off {:.3}. The last row is the\n\
         control — a support covering every allele the world has must read exactly zero.\n",
        100.0 * truth.level,
        truth.up_share,
        truth.falloff
    );
    let _ = writeln!(
        report,
        "{:<10} {:>10} {:>12} {:>12} {:>12} {:>10}",
        "support", "outside", "level", "up-share", "fall-off", "spread"
    );

    // A spectrum with a tail: most loci at the reference, and a geometric fall away from it, which
    // is the shape the real distribution has (research note §6.8).
    let offsets: Vec<i32> = (-4..=4).collect();
    let probs: Vec<f64> = {
        let mut raw: Vec<f64> = offsets
            .iter()
            .map(|&a| {
                if a == 0 {
                    1.0
                } else {
                    0.09 * 0.45f64.powi(a.abs() - 1)
                }
            })
            .collect();
        let total: f64 = raw.iter().sum();
        for p in &mut raw {
            *p /= total;
        }
        raw
    };
    let world = World {
        name: "allele support".to_string(),
        slip: truth,
        allele_offsets: offsets.clone(),
        allele_probs: probs.clone(),
        inbreeding: 0.0,
        mean_depth: 6.0,
    };
    let depth_probs = world.depth_probs(8);
    let genotypes = world.genotypes();
    // **Recorded at ±2 while the alleles reach ±4**, deliberately: nine allele lengths against nine
    // offset buckets makes the cell space forty-five genotypes over ninety thousand cells, which is
    // hours a fit. §6.4 already measured that a range narrower than the alleles costs nothing under
    // the marginal rule, so folding them into the end buckets isolates the question this asks
    // rather than confounding it — and the control row is what confirms that.
    let keying = Keying {
        half_range: 2,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };
    let table = build_cells(
        &truth,
        &world.allele_offsets,
        &genotypes,
        &depth_probs,
        &keying,
    );
    let weights = truth_mass(&table, &world.genotype_probs());

    for limit in [0i32, 1, 2, 3, 4] {
        let model_alleles: Vec<i32> = (-limit..=limit).collect();
        // How much of the locus population carries an allele the fit cannot place. A locus is
        // outside if **either** copy is.
        let inside: f64 = world
            .genotypes()
            .into_iter()
            .zip(world.genotype_probs())
            .filter(|((i, j), _)| offsets[*i].abs() <= limit && offsets[*j].abs() <= limit)
            .map(|(_, p)| p)
            .sum();
        let inputs = FitInputs {
            world: &world,
            truth_cells: &table.cells,
            truth_weights: &weights,
            model_keying: keying,
            depth_probs: &depth_probs,
            model_alleles: Some(&model_alleles),
        };
        let (best, spread) = fit_from_spread(&inputs);
        let _ = writeln!(
            report,
            "±{:<9} {:>9.2}% {:>+11.1}% {:>+12.4} {:>+12.4} {:>9.2}×",
            limit,
            100.0 * (1.0 - inside),
            100.0 * (best.slip.level / truth.level - 1.0),
            best.slip.up_share - truth.up_share,
            best.slip.falloff - truth.falloff,
            spread.level
        );
    }
}

/// **Does the scoring rule stay unbiased above twelve reads a locus?**
///
/// Every other question here stops at a cap of twelve, because a locus's cell space is every way
/// its reads split across the buckets and that grows as the depth to the power of the bucket count.
/// That left `MAX_LOCUS_READS` justified by where the evidence stopped rather than by anything
/// about the data — an unsatisfying reason for a constant that decides how many reads are thrown
/// away at a deep locus.
///
/// **Narrowing the recorded range is what buys the depth**, and it is free to do: three buckets
/// instead of nine takes the cell space at 60 reads from eight million to 39,710, and the
/// saturating-bucket measurement already showed a range of ±1 scored by its marginal is exactly
/// unbiased. So the rule can be checked at real depths, on a key that is narrower but not weaker.
fn question_depth(report: &mut String) {
    let _ = writeln!(report, "\n## Whether the scoring rule survives depth\n");
    let truth = measured_slip();
    let _ = writeln!(
        report,
        "Loci keyed by locus, offsets recorded over **±1** with the ends scored by their marginal —\n\
         narrow on purpose, because that is what keeps the cell space affordable to 60 reads a\n\
         locus. `beats truth` is how much better the fitted answer scores than the truth, which a\n\
         correctly specified model cannot make positive; a value that grows with depth is the\n\
         inner climb slowing down, not the estimator failing.\n"
    );
    let _ = writeln!(
        report,
        "{:>12} {:>8} {:>11} {:>12} {:>12} {:>9} {:>12}",
        "reads/locus", "cells", "level", "up-share", "fall-off", "spread", "beats truth"
    );

    let keying = Keying {
        half_range: 1,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };
    for &(mean_depth, cap) in &[
        (3.0f64, 12u32),
        (6.0, 20),
        (12.0, 30),
        (20.0, 40),
        (30.0, 55),
        (45.0, 75),
    ] {
        let w = world("depth", truth, 0.10, 0.0, mean_depth);
        let depth_probs = w.depth_probs(cap);
        let genotypes = w.genotypes();
        let table = build_cells(&truth, &w.allele_offsets, &genotypes, &depth_probs, &keying);
        let weights = truth_mass(&table, &w.genotype_probs());
        let inputs = FitInputs {
            world: &w,
            truth_cells: &table.cells,
            truth_weights: &weights,
            model_keying: keying,
            depth_probs: &depth_probs,
            model_alleles: None,
        };
        let (best, spread) = fit_from_spread(&inputs);
        let beats_truth = best.score - score_slip(&inputs, &truth).0;
        let _ = writeln!(
            report,
            "{:>12.1} {:>8} {:>+10.2}% {:>+12.4} {:>+12.4} {:>8.2}× {:>12.2e}",
            w.realised_mean_depth(cap),
            table.cells.len(),
            100.0 * (best.slip.level / truth.level - 1.0),
            best.slip.up_share - truth.up_share,
            best.slip.falloff - truth.falloff,
            spread.level,
            beats_truth
        );
    }
}

// ---------------------------------------------------------------------------
// Question 3 — the saturating end buckets
// ---------------------------------------------------------------------------

/// **How a saturated end bucket is scored.** The accumulator records offsets over a small
/// range with the ends absorbing everything beyond; treating "at least this far" as
/// "exactly this far" is a plug-in, and the marginal is the sum over everything the bucket
/// takes in.
fn question_saturation(report: &mut String) {
    let _ = writeln!(report, "\n## Scoring a saturated end bucket\n");
    let _ = writeln!(
        report,
        "The end buckets hold every read that moved at least that far. `Marginal` scores a\n\
         bucket by the sum over everything it absorbs; `PlugAtEdge` scores it as though every\n\
         read sat exactly on the edge, and `PlugAtEdgeRenormalised` rescales that to a proper\n\
         distribution. The narrower the range, the more the end buckets absorb.\n"
    );
    let truth = measured_slip();
    // A fall-off three times the measured one, so that a meaningful share of the slipped
    // reads reaches the second and third step and the end buckets have something to absorb.
    let heavy = Slip {
        falloff: 0.30,
        ..truth
    };
    // The third block is the case the reference origin makes unavoidable: with offsets
    // measured from the reference tract length, the range has to cover the **allele**
    // spectrum and not only the slippage, so the end buckets absorb whole alleles.
    for (label, slip, allele_offsets) in [
        ("measured fall-off", truth, vec![-1i32, 0, 1]),
        ("fall-off 0.30", heavy, vec![-1, 0, 1]),
        (
            "alleles wider than the range",
            truth,
            vec![-3, -2, -1, 0, 1, 2, 3],
        ),
    ] {
        let _ = writeln!(
            report,
            "\n### {label} — {:.1}% of reads slip, {:.0} in 100 of those take a second step, \
             alleles {}..{}\n",
            100.0 * slip.level,
            100.0 * slip.falloff,
            allele_offsets[0],
            allele_offsets[allele_offsets.len() - 1]
        );
        let _ = writeln!(
            report,
            "{:<10} {:<26} {:>11} {:>11} {:>11} {:>10}",
            "range", "edge scoring", "level", "up-share", "fall-off", "hets"
        );
        // The wide-allele world enumerates 28 genotypes against seven allele lengths, so it
        // is run shallower and at the two ranges the question is about — a range wider than
        // the alleles is not the case in doubt.
        let wide = allele_offsets.len() > 3;
        let ranges: &[i32] = if wide { &[1, 2] } else { &[1, 2, 3] };
        let mean_depth = if wide { 3.0 } else { 6.0 };
        let depth_cap = if wide { 8 } else { 10 };
        for &half_range in ranges {
            for edges in [
                EdgeScoring::Marginal,
                EdgeScoring::PlugAtEdge,
                EdgeScoring::PlugAtEdgeRenormalised,
            ] {
                let mut w = world("saturation", slip, 0.10, 0.0, mean_depth);
                if wide {
                    let spread = 0.30;
                    w.allele_offsets = allele_offsets.clone();
                    w.allele_probs = allele_offsets
                        .iter()
                        .map(|&a| {
                            if a == 0 {
                                1.0 - spread
                            } else {
                                spread / (allele_offsets.len() - 1) as f64
                            }
                        })
                        .collect();
                }
                let depth_probs = w.depth_probs(depth_cap);
                let genotypes = w.genotypes();
                let truth_keying = Keying {
                    half_range,
                    origin: Origin::Reference,
                    edges: EdgeScoring::Marginal,
                };
                let model_keying = Keying {
                    half_range,
                    origin: Origin::Reference,
                    edges,
                };
                let table = build_cells(
                    &slip,
                    &w.allele_offsets,
                    &genotypes,
                    &depth_probs,
                    &truth_keying,
                );
                let weights = truth_mass(&table, &w.genotype_probs());
                let inputs = FitInputs {
                    world: &w,
                    truth_cells: &table.cells,
                    truth_weights: &weights,
                    model_keying,
                    depth_probs: &depth_probs,
                    model_alleles: None,
                };
                let (best, _) = fit_from_spread(&inputs);
                let _ = writeln!(
                    report,
                    "±{:<9} {:<26} {:>+10.2}% {:>+11.4} {:>+11.4} {:>+9.2}%",
                    half_range,
                    format!("{edges:?}"),
                    100.0 * (best.slip.level / slip.level - 1.0),
                    best.slip.up_share - slip.up_share,
                    best.slip.falloff - slip.falloff,
                    100.0
                        * (best.heterozygosity(w.allele_offsets.len()) / w.heterozygosity() - 1.0),
                );
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Question 4 — the shape of the surface the search walks
// ---------------------------------------------------------------------------

/// **Does the score curve have one hump?** The design steps through the noise parameters
/// end to end precisely because nobody has shown it does. At one parameter that costs 161
/// scores; at three it costs 4.2 million, per (read group × stratum), so whether the answer
/// is one hump decides whether the scan can be replaced.
fn question_shape(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## The shape of the surface, and what a scan would cost\n"
    );
    let _ = writeln!(
        report,
        "The curve below is a **profile**: at each slippage level the other two parameters and\n\
         the genotype frequencies are maximised out, leaving a curve in the level alone. A\n\
         second hump would mean a scan coarse enough to afford could step over the answer.\n"
    );
    let truth = measured_slip();
    let keying = Keying {
        half_range: 2,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };
    for &(depth, spread) in &[(3.0f64, 0.30f64), (12.0, 0.05)] {
        let w = world("shape", truth, spread, 0.0, depth);
        let depth_probs = w.depth_probs(12);
        let genotypes = w.genotypes();
        let table = build_cells(&truth, &w.allele_offsets, &genotypes, &depth_probs, &keying);
        let weights = truth_mass(&table, &w.genotype_probs());
        let inputs = FitInputs {
            world: &w,
            truth_cells: &table.cells,
            truth_weights: &weights,
            model_keying: keying,
            depth_probs: &depth_probs,
            model_alleles: None,
        };

        // A ladder in the level, with the other two axes maximised at each rung.
        let rungs = 41;
        let (lo, hi) = ((1e-4f64).ln(), (0.3f64).ln());
        let mut curve = Vec::with_capacity(rungs);
        for step in 0..rungs {
            let level = (lo + (hi - lo) * step as f64 / (rungs - 1) as f64).exp();
            // Maximise the two shares at this level, from two starts.
            let mut best = f64::NEG_INFINITY;
            for start in [(0.20, 0.03), (0.80, 0.40)] {
                let mut current = Slip {
                    level,
                    up_share: start.0,
                    falloff: start.1,
                };
                let inverse_phi = 0.5 * (5f64.sqrt() - 1.0);
                for _sweep in 0..3 {
                    for axis in 1..3 {
                        let (mut a, mut b) = if axis == 1 {
                            (logit(0.005), logit(0.995))
                        } else {
                            (logit(0.002), logit(0.95))
                        };
                        let mut c = b - inverse_phi * (b - a);
                        let mut d = a + inverse_phi * (b - a);
                        let evaluate = |x: f64, base: &Slip| {
                            let mut trial = *base;
                            if axis == 1 {
                                trial.up_share = expit(x);
                            } else {
                                trial.falloff = expit(x);
                            }
                            score_slip(&inputs, &trial).0
                        };
                        let mut fc = evaluate(c, &current);
                        let mut fd = evaluate(d, &current);
                        for _ in 0..26 {
                            if fc > fd {
                                b = d;
                                d = c;
                                fd = fc;
                                c = b - inverse_phi * (b - a);
                                fc = evaluate(c, &current);
                            } else {
                                a = c;
                                c = d;
                                fc = fd;
                                d = a + inverse_phi * (b - a);
                                fd = evaluate(d, &current);
                            }
                            if (b - a) < 1e-6 {
                                break;
                            }
                        }
                        let chosen = 0.5 * (a + b);
                        if axis == 1 {
                            current.up_share = expit(chosen);
                        } else {
                            current.falloff = expit(chosen);
                        }
                    }
                }
                best = best.max(score_slip(&inputs, &current).0);
            }
            curve.push((level, best));
        }
        let peak = curve
            .iter()
            .cloned()
            .fold((0.0f64, f64::NEG_INFINITY), |acc, x| {
                if x.1 > acc.1 { x } else { acc }
            });
        let humps = curve
            .windows(3)
            .filter(|w| w[1].1 > w[0].1 + 1e-12 && w[1].1 > w[2].1 + 1e-12)
            .count();
        let _ = writeln!(
            report,
            "\n### mean depth {depth:.0}, hets {:.0} in 100",
            100.0 * w.heterozygosity()
        );
        let _ = writeln!(
            report,
            "  {rungs} rungs from 0.0001 to 0.3; interior local maxima: **{humps}**; \
             peak at level {:.4} against a truth of {:.4}",
            peak.0, truth.level
        );
        // A readable slice of the curve, relative to the peak.
        let _ = writeln!(report, "\n{:<12} {:>16}", "level", "score − peak");
        for (level, score) in curve.iter().step_by(4) {
            let _ = writeln!(report, "{level:<12.5} {:>16.6}", score - peak.1);
        }
    }
}

// ---------------------------------------------------------------------------
// Question 5 — merging two strata
// ---------------------------------------------------------------------------

/// **What pooling two strata costs.** Thin strata borrow their neighbours' value, and a
/// fitted sequence that fails to rise with repeat count is merged and refitted. Both change
/// the estimate; neither has ever had a bias measured against it.
fn question_strata(report: &mut String) {
    let _ = writeln!(report, "\n## Borrowing and merging across strata\n");
    let _ = writeln!(
        report,
        "Two strata whose true slippage differs by `ratio`, pooled and fitted as one. The\n\
         merged answer is what a merge-and-refit produces when it fires, and the two columns\n\
         say what each stratum then carries. `share` is the first stratum's share of the loci.\n"
    );
    let _ = writeln!(
        report,
        "{:<10} {:>8} {:>12} {:>12} {:>14} {:>14}",
        "ratio", "share", "truth low", "truth high", "merged fit", "worst error"
    );
    let base = measured_slip();
    let keying = Keying {
        half_range: 2,
        origin: Origin::Reference,
        edges: EdgeScoring::Marginal,
    };
    for &ratio in &[1.0f64, 1.5, 2.0, 4.0] {
        for &share in &[0.5f64, 0.8] {
            let low = Slip {
                level: base.level,
                ..base
            };
            let high = Slip {
                level: base.level * ratio,
                ..base
            };
            let w = world("strata", base, 0.10, 0.0, 3.0);
            let depth_probs = w.depth_probs(12);
            let genotypes = w.genotypes();
            let table_low = build_cells(&low, &w.allele_offsets, &genotypes, &depth_probs, &keying);
            let table_high =
                build_cells(&high, &w.allele_offsets, &genotypes, &depth_probs, &keying);
            let genotype_probs = w.genotype_probs();
            // Pool the two strata into one cell table, weighted by their share of the loci.
            let mut pooled: HashMap<LocusCell, f64> = HashMap::new();
            for (table, weight) in [(&table_low, share), (&table_high, 1.0 - share)] {
                for (slot, cell) in table.cells.iter().enumerate() {
                    let mass: f64 = table.component[slot]
                        .iter()
                        .zip(&genotype_probs)
                        .map(|(c, p)| c * p)
                        .sum();
                    *pooled.entry(cell.clone()).or_insert(0.0) += weight * mass;
                }
            }
            let mut cells: Vec<LocusCell> = pooled.keys().cloned().collect();
            cells.sort();
            let weights: Vec<f64> = cells.iter().map(|c| pooled[c]).collect();
            let inputs = FitInputs {
                world: &w,
                truth_cells: &cells,
                truth_weights: &weights,
                model_keying: keying,
                depth_probs: &depth_probs,
                model_alleles: None,
            };
            let (best, _) = fit_from_spread(&inputs);
            let worst = ((best.slip.level / low.level - 1.0).abs())
                .max((best.slip.level / high.level - 1.0).abs());
            let _ = writeln!(
                report,
                "{ratio:<10.1} {:>7.0}% {:>11.3}% {:>11.3}% {:>13.3}% {:>13.0}%",
                100.0 * share,
                100.0 * low.level,
                100.0 * high.level,
                100.0 * best.slip.level,
                100.0 * worst
            );
        }
    }
    let _ = writeln!(
        report,
        "\nThe first two rows are the control: with the two strata identical a merge must cost\n\
         exactly nothing."
    );
}

// ---------------------------------------------------------------------------
// Question 6 — the composition channel
// ---------------------------------------------------------------------------

/// **`ε` from two scalars per stratum.** The accumulator keeps how many bases were compared
/// and how many mismatched, and nothing else. This asks what that coarse a key costs.
fn question_composition(report: &mut String) {
    let _ = writeln!(
        report,
        "\n## The substitution rate from two running counts\n"
    );
    let _ = writeln!(
        report,
        "Each read is compared against the tract at the length **that read shows**, so a\n\
         mismatch is a substitution and not a slip. Two consequences, both arithmetic rather\n\
         than measurement, and both checked below.\n"
    );

    // 1. Separability: the mismatch count is conditionally independent of the length
    //    outcome, so the joint likelihood factorises and ε has a closed form.
    let bases = 40u32;
    let eps_true = 0.003f64;
    let mut best_eps = 0.0;
    let mut best_score = f64::NEG_INFINITY;
    for step in 0..=4000 {
        let eps = 1e-5 + step as f64 * (0.05 - 1e-5) / 4000.0;
        // The exact expectation of the pooled log-likelihood under the truth.
        let mut score = 0.0;
        for k in 0..=bases {
            let ln_choose = ln_factorial(bases) - ln_factorial(k) - ln_factorial(bases - k);
            let p_true = (ln_choose
                + x_ln_y(f64::from(k), eps_true)
                + x_ln_y(f64::from(bases - k), 1.0 - eps_true))
            .exp();
            score += p_true * (x_ln_y(f64::from(k), eps) + x_ln_y(f64::from(bases - k), 1.0 - eps));
        }
        if score > best_score {
            best_score = score;
            best_eps = eps;
        }
    }
    let _ = writeln!(
        report,
        "**Pooling every read's bases into two counters is exact.** A read's mismatches are\n\
         binomial at `ε` whatever its length, so the pooled counts are a sufficient statistic\n\
         and the maximum-likelihood value is mismatches over bases compared. Truth\n\
         {eps_true:.4}, recovered {best_eps:.4} — a division, not a search."
    );

    // 2. What pooling reads of two different true rates returns.
    let _ = writeln!(
        report,
        "\n**What it returns when the reads are not all alike.** Reads inside a repeat tract\n\
         carry worse base quality than reads outside one, so a stratum may hold two\n\
         populations. The pooled counter returns the **base-weighted mean** of their rates,\n\
         which is the right answer for a model with one rate — the same result the generic\n\
         path's shared error rate has (research note §2.2).\n"
    );
    let _ = writeln!(
        report,
        "{:<28} {:>12} {:>12} {:>14}",
        "two populations", "rate A", "rate B", "pooled returns"
    );
    for &(share_a, rate_a, rate_b) in &[
        (0.5f64, 0.001f64, 0.001f64),
        (0.5, 0.001, 0.004),
        (0.9, 0.001, 0.010),
    ] {
        let pooled = share_a * rate_a + (1.0 - share_a) * rate_b;
        let _ = writeln!(
            report,
            "{:<28} {:>12.4} {:>12.4} {:>14.5}",
            format!("{:.0}% at A", 100.0 * share_a),
            rate_a,
            rate_b,
            pooled
        );
    }
    let _ = writeln!(
        report,
        "\nThe first row is the control: with the two populations identical the pooled rate\n\
         must be that rate exactly."
    );

    // 3. The scan arithmetic that follows.
    let _ = writeln!(
        report,
        "\n**What this removes from the search.** `ε` is fitted from its own two counters in\n\
         closed form, so it is not an axis of the scan. What is left to search is the\n\
         slippage: how often a read slips, which way, and how far — three axes, not the\n\
         error rate plus two."
    );
}

// ---------------------------------------------------------------------------

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let wanted = |section: &str| {
        let selected: Vec<&String> = args.iter().filter(|a| a.starts_with("--only=")).collect();
        selected.is_empty() || selected.iter().any(|a| a.ends_with(section))
    };

    println!("# The STR stutter accumulator — is it able to give an unbiased answer?");
    println!("#");
    println!("# Every cell is weighted by its exact probability under a known truth, so each");
    println!("# number below is what an infinite genome returns: bias with no sampling error");
    println!("# in it, and a zero that means the estimator is consistent rather than lucky.");

    let mut report = String::new();
    if wanted("gates") {
        question_gates(&mut report);
    }
    if wanted("identify") {
        question_identification(&mut report);
    }
    if wanted("origin") {
        question_origin(&mut report);
    }
    if wanted("converge") {
        question_convergence(&mut report);
    }
    if wanted("support") {
        question_allele_support(&mut report);
    }
    if wanted("depth") {
        question_depth(&mut report);
    }
    if wanted("saturation") {
        question_saturation(&mut report);
    }
    if wanted("shape") {
        question_shape(&mut report);
    }
    if wanted("strata") {
        question_strata(&mut report);
    }
    if wanted("composition") {
        question_composition(&mut report);
    }
    print!("{report}");
}
