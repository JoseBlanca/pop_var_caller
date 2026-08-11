//! **A stratum's loci, and the per-genotype likelihood of each — built from observed shapes.**
//!
//! This is the half of the STR stutter fit that a real-data program cannot borrow from the
//! exact-bias harness. The harness enumerates every way a locus's reads could split across the
//! offset buckets and weights each by its probability; that space is `C(depth + 8, 8)` shapes at one
//! depth — 220 at three reads, 294,000 at twelve, and beyond counting at the thirty a deep library
//! reaches. A program reading real alignments instead **tallies the shapes it saw** and evaluates
//! the model only there.
//!
//! ## What an entry is, and why it is a locus
//!
//! One entry is one locus's shape: how many of that locus's reads showed each whole-repeat offset
//! from the **reference** tract length, plus how many showed a length that is not a whole number of
//! copies. Loci with identical shapes are counted together.
//!
//! Keeping a locus's reads together is what makes the fit possible at all. A read carries no
//! genotype — it drew one of the locus's two alleles and then slipped — so a tally that pools reads
//! across loci holds the allele spectrum convolved with the slippage kernel, and the fitted rate
//! then moves **333-fold depending only on where the search starts** (spec §4.1). A real allele
//! shifts all of a locus's reads; slippage shifts some of them. That is the whole separation.
//!
//! ## The likelihood of one entry
//!
//! Given a genotype — an unordered pair of allele offsets — each read picks one of the two copies
//! and then slips, so the per-bucket probability is the average of the two copies' slip kernels and
//! the shape's probability is one multinomial over the buckets.
//!
//! **Two factors are dropped and the reason is that they cancel exactly.** The multinomial
//! coefficient and the probability of the observed depth depend on the entry alone, not on the
//! genotype and not on the slippage parameters. Genotype-independent factors cancel in the climb
//! over the frequencies, which normalises per entry; and factors that do not depend on the
//! candidate parameters shift every candidate's score by the same amount, so they cannot move the
//! maximum. **What this costs is that a score from here is a log-likelihood up to an additive
//! constant** — differences between candidates are exact, absolute values are not comparable with
//! the harness's.

#![allow(dead_code)]

use std::collections::BTreeMap;

#[path = "stutter_model.rs"]
mod stutter_model;

pub use stutter_model::*;

/// One locus's reads, laid out across the offset buckets.
///
/// `counts.iter().sum() + not_whole_repeat == depth`, always — the invariant the "no bucket is
/// charged a negative number of reads" gate checks. Ordered and hashable so it can key a table and
/// so iteration order is fixed, which is the whole of the determinism requirement.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct LocusShape {
    /// Reads at each bucket, in bucket order. The middle index is the reference length.
    pub counts: Vec<u16>,
    /// Reads differing from the reference by something that is not a whole number of copies.
    ///
    /// **A guard, not a parameter.** It is modelled as an independent per-read outcome, so the
    /// likelihood splits exactly into *how many reads were non-whole-repeat* times *how the rest
    /// fell across the offsets* — nothing about the slippage parameters is estimated from it and
    /// nothing about it disturbs them (spec §4.1). It is carried so the fitter can report it.
    pub not_whole_repeat: u16,
}

impl LocusShape {
    /// Reads that landed in an offset bucket — the depth the length channel sees.
    pub fn scored_depth(&self) -> u32 {
        self.counts.iter().map(|&c| u32::from(c)).sum()
    }
}

/// One stratum's evidence: every distinct locus shape and how many loci had it.
///
/// Sparse rather than dense, and that is forced rather than chosen — the possible shape space is
/// vast and only a small, data-dependent corner of it is ever occupied. A `BTreeMap` rather than a
/// `HashMap` because every fit is a floating-point sum over entries and floating-point addition is
/// not associative, so a fixed iteration order is what makes two runs agree.
#[derive(Default, Clone, Debug)]
pub struct StratumTable {
    loci_by_shape: BTreeMap<LocusShape, u64>,
    loci: u64,
    reads: u64,
    /// Reads that landed off the reference bucket, whole-repeat or not — the moment estimate the
    /// search starts from.
    reads_off_reference: u64,
    reads_not_whole_repeat: u64,
}

impl StratumTable {
    pub fn add_locus(&mut self, shape: LocusShape) {
        let depth = u64::from(shape.scored_depth()) + u64::from(shape.not_whole_repeat);
        if depth == 0 {
            return;
        }
        let middle = shape.counts.len() / 2;
        let at_reference = u64::from(shape.counts[middle]);
        self.loci += 1;
        self.reads += depth;
        self.reads_off_reference += depth - at_reference;
        self.reads_not_whole_repeat += u64::from(shape.not_whole_repeat);
        *self.loci_by_shape.entry(shape).or_insert(0) += 1;
    }

    pub fn loci(&self) -> u64 {
        self.loci
    }
    pub fn reads(&self) -> u64 {
        self.reads
    }
    pub fn entries(&self) -> usize {
        self.loci_by_shape.len()
    }

    /// Of the reads differing from the reference length, the share that did so by something other
    /// than a whole number of copies. §5's diagnostic: above one in ten, this noise model does not
    /// describe the stratum and a rate fitted from it is mostly mis-modelled ordinary indel.
    pub fn not_whole_repeat_share(&self) -> f64 {
        if self.reads_off_reference == 0 {
            0.0
        } else {
            self.reads_not_whole_repeat as f64 / self.reads_off_reference as f64
        }
    }

    /// **Where the search starts from: the share of reads sitting off the reference bucket.**
    ///
    /// An over-estimate of the slippage level, because it counts reads that correctly measured a
    /// non-reference allele as well as slipped ones — which is why the starting points spread below
    /// it as well as above (spec §4.2).
    pub fn level_estimate(&self) -> f64 {
        if self.reads == 0 {
            0.01
        } else {
            (self.reads_off_reference as f64 / self.reads as f64).clamp(1e-4, 0.3)
        }
    }

    /// Every distinct shape and how many loci had it, materialised once. A `Vec` rather than an
    /// iterator because the search re-walks these once per candidate.
    pub fn shapes(&self) -> Vec<(&LocusShape, u64)> {
        self.loci_by_shape.iter().map(|(s, &n)| (s, n)).collect()
    }
}

/// The allele offsets the fit may place mass on, for a stratum of this reference repeat count.
///
/// `±limit` around the reference length, **clipped at the low end because an allele cannot be
/// shorter than nothing**: a stratum at 4 repeats reaches only −4. So the genotype count is a
/// per-stratum quantity rather than a constant.
///
/// The limit matters and the recorded bucket range does not: measured, a locus whose allele the fit
/// cannot represent has its reads explained the only way left, as slippage, and the cost is nothing
/// and then everything — 2.5% of loci outside costs 0.1% of the rate, 19% costs +499% with the
/// direction asymmetry destroyed (spec §8.1).
pub fn allele_support(reference_repeats: u32, limit: i32) -> Vec<i32> {
    let lowest = -(limit.min(reference_repeats as i32));
    (lowest..=limit).collect()
}

/// **The per-genotype likelihood of every observed shape, at these slippage parameters.**
///
/// `component[entry][genotype]`, and the entries are in `table.shapes()` order. The multinomial
/// coefficient and the depth's own probability are omitted for the reason the module doc gives:
/// both are genotype-independent and candidate-independent, so they cancel exactly.
pub fn component_matrix(
    shapes: &[(&LocusShape, u64)],
    slip: &Slip,
    alleles: &[i32],
    keying: &Keying,
) -> LnComponents {
    let genotypes = genotype_pairs(alleles.len());
    // The bucket distribution per genotype, logged once and reused across every entry — a shape's
    // likelihood is a sum of `count · ln p` over the buckets, so the logs are what the loop wants.
    let ln_probs: Vec<Vec<f64>> = genotypes
        .iter()
        .map(|&g| {
            genotype_bucket_probs(slip, alleles, g, keying)
                .into_iter()
                .map(|p| p.max(1e-300).ln())
                .collect()
        })
        .collect();
    let mut ln = Vec::with_capacity(shapes.len() * genotypes.len());
    for (shape, _) in shapes {
        for probs in &ln_probs {
            let mut total = 0.0;
            for (&count, &ln_p) in shape.counts.iter().zip(probs) {
                if count > 0 {
                    total += f64::from(count) * ln_p;
                }
            }
            ln.push(total);
        }
    }
    LnComponents {
        ln,
        genotypes: genotypes.len(),
    }
}

/// What one stratum's fit returned.
#[derive(Clone, Debug)]
pub struct StratumFit {
    pub slip: Slip,
    /// The ratio between the highest and lowest level any starting point reached.
    pub level_spread: f64,
    pub score: f64,
    /// How often the fitted genotype frequencies say a locus carries two different alleles. A
    /// by-product of the fit rather than something this path emits, and useful as a sanity read.
    pub heterozygosity: f64,
    pub loci: u64,
    pub reads: u64,
    pub not_whole_repeat_share: f64,
}

/// **Fit the three slippage parameters for one stratum.**
///
/// The genotype is summed over rather than chosen, so the fit alternates: hold the slippage
/// parameters, climb to the allele-pair frequencies that best explain the stratum's loci — the half
/// with a concavity proof behind it, so it cannot stop on a false summit — then move the slippage
/// parameters and repeat, from several starting points.
pub fn fit_stratum(
    table: &StratumTable,
    reference_repeats: u32,
    allele_limit: i32,
    keying: &Keying,
    climb_iterations: usize,
    precision: SearchPrecision,
) -> StratumFit {
    let alleles = allele_support(reference_repeats, allele_limit);
    let shapes = table.shapes();
    let weights: Vec<f64> = shapes.iter().map(|(_, n)| *n as f64).collect();

    // **The climb starts cold every time, and that is not an oversight.** Seeding it from the
    // previous candidate's answer is an obvious speed-up and it is wrong here: the outer search
    // compares two candidates by their scores, so the objective has to be a function of the
    // candidate alone. Warm-started, it also depends on which candidates were tried before, and
    // golden-section then compares two numbers computed under different conditions — a search that
    // wanders rather than converges. Measured: with warm starts the exact control's fitted point
    // scored *below* the truth and moved 1.3% when the climb's cap was raised tenfold, both of
    // which say the maximum was never reached.
    let objective = |slip: &Slip| {
        let component = component_matrix(&shapes, slip, &alleles, keying);
        let freqs = climb_genotype_frequencies(&component, &weights, climb_iterations, None);
        let score = score_at(&component, &weights, &freqs);
        (score, freqs)
    };

    let fit = fit_from_starts(
        &objective,
        &starting_points(table.level_estimate()),
        precision,
    );
    let heterozygosity = genotype_pairs(alleles.len())
        .into_iter()
        .zip(&fit.genotype_freqs)
        .filter(|((i, j), _)| i != j)
        .map(|(_, &p)| p)
        .sum();

    StratumFit {
        slip: fit.slip,
        level_spread: fit.level_spread,
        score: fit.score,
        heterozygosity,
        loci: table.loci(),
        reads: table.reads(),
        not_whole_repeat_share: table.not_whole_repeat_share(),
    }
}

/// Score a stratum's table at a *given* slippage, climbing only the frequencies.
///
/// **The check that a spread across starting points cannot make.** A correctly specified model
/// cannot be beaten at the truth, so where a truth exists this is evaluated there and at the fitted
/// point: a fitted score above the truth's is a defect in the code, not a finding about the
/// estimator. A deterministic search returns the same point from every start wherever the objective
/// is flat, which is how a fall-off 0.02 out once passed four agreeing starts (spec §4.2).
pub fn score_stratum_at(
    table: &StratumTable,
    slip: &Slip,
    reference_repeats: u32,
    allele_limit: i32,
    keying: &Keying,
    climb_iterations: usize,
) -> f64 {
    let alleles = allele_support(reference_repeats, allele_limit);
    let shapes = table.shapes();
    let weights: Vec<f64> = shapes.iter().map(|(_, n)| *n as f64).collect();
    let component = component_matrix(&shapes, slip, &alleles, keying);
    let freqs = climb_genotype_frequencies(&component, &weights, climb_iterations, None);
    score_at(&component, &weights, &freqs)
}
