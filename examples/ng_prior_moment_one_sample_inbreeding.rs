//! **Two populations that differ five-fold in diversity and that one genome cannot tell apart.**
//!
//! # Why this exists
//!
//! `examples/ng_prior_moment_estimators.rs` finds that the heterozygosity estimator the research
//! plan proposes comes back **one fifth of the truth at one self-pollinated individual**, because
//! the estimator counts how often two chromosomes drawn from the panel differ and at one
//! individual the only pair available is the two copies inside it — which are the same ancestral
//! copy with probability `F`, the inbreeding coefficient. The fix is to divide by `1 − F`.
//!
//! **That raises the obvious question: can `F` be estimated from the one genome?** It cannot, and
//! this program shows it in the strongest available form — not by fitting and finding the answer
//! unstable, but by constructing **two populations whose single-genome data have exactly the same
//! distribution**, position for position and read for read, while their diversities differ by a
//! factor of five.
//!
//! **If the data distributions are identical then no estimator can separate them**, so the
//! conclusion covers the direct estimate, the fitted curve the caller uses today, and anything
//! anyone writes later. What this program measures is only that the construction is right and that
//! the shipped fit behaves as the argument says it must.
//!
//! # The construction
//!
//! A population is a share of positions carrying only the reference base, a share carrying only a
//! non-reference base, and a `Beta(a, b)` over what segregates. One genome from it shows exactly
//! two things across the census: how often it is heterozygous, and how often both its copies are
//! non-reference. Writing `s` for the segregating share, `q` for the fixed-non-reference share and
//! `m1`, `m2` for the Beta's first two moments:
//!
//! ```text
//! heterozygous          (1 - F) · s · 2 (m1 - m2)
//! both copies non-ref   q + s · m2 + F · s · (m1 - m2)
//! ```
//!
//! Hold the Beta fixed and both of those equal between an outbred population and a selfing one,
//! and the two equations solve:
//!
//! ```text
//! s_outbred = (1 - F) · s_selfing
//! q_outbred = q_selfing + F · s_selfing · m1
//! ```
//!
//! The population diversity is `s · 2 (m1 - m2)`, which is in the ratio `1 - F`. At `F = 0.8`
//! that is a factor of five.
//!
//! Run: `./scripts/dev.sh cargo run --release --example ng_prior_moment_one_sample_inbreeding`

use std::collections::BTreeMap;
use std::env;

use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, NamedReadGroup, ObservedAllele, PackedDepthCodes, ReadCap, RecordingTerms,
    SampleCensusEvidence, Section, SectionKey, SelectionTermsDigest,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    JointFitConfig, StartingPoint, fit_jointly,
};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::ReadGroupId;

/// The selfing population's inbreeding coefficient — tomato's fitted range.
const SELFING_COEFFICIENT: f64 = 0.8;

/// The Beta over the segregating positions. **The same in both populations**, so the only thing
/// that differs is how many positions segregate and how many are fixed non-reference.
const BETA_A: f64 = 0.20;
const BETA_B: f64 = 1.00;

/// The selfing population's shares.
const SELFING_SEGREGATING: f64 = 0.05;
const SELFING_FIXED_ALT: f64 = 0.001;

/// How often a read misreads a base.
const CLEAN_ERROR_RATE: f64 = 0.002;

/// The census's per-position read ceiling — see `ng_prior_moments_from_reads.rs`.
const CENSUS_DEPTH_CAP: u32 = 124;

fn main() {
    let mut args = env::args().skip(1);
    let positions: usize = args.next().map_or(200_000, |a| a.parse().expect("a count"));

    let m1 = BETA_A / (BETA_A + BETA_B);
    let m2 = BETA_A * (BETA_A + 1.0) / ((BETA_A + BETA_B) * (BETA_A + BETA_B + 1.0));

    let selfing = Population {
        name: "selfing",
        inbreeding: SELFING_COEFFICIENT,
        segregating: SELFING_SEGREGATING,
        fixed_alt: SELFING_FIXED_ALT,
    };
    let outbred = Population {
        name: "outbred",
        inbreeding: 0.0,
        segregating: (1.0 - SELFING_COEFFICIENT) * SELFING_SEGREGATING,
        fixed_alt: SELFING_FIXED_ALT + SELFING_COEFFICIENT * SELFING_SEGREGATING * m1,
    };

    println!("# Two populations one genome cannot tell apart");
    println!();
    println!(
        "Both spread their segregating positions as Beta({BETA_A}, {BETA_B}); they differ only in \
         how many positions segregate, how many are fixed non-reference, and how inbred the \
         individuals are."
    );
    println!();
    println!(
        "| population | F | segregating | fixed non-reference | population diversity | per kb |"
    );
    println!("|---|---:|---:|---:|---:|---:|");
    for population in [&selfing, &outbred] {
        println!(
            "| {} | {} | {:.5} | {:.5} | {:.6} | {:.2} |",
            population.name,
            population.inbreeding,
            population.segregating,
            population.fixed_alt,
            population.diversity(m1, m2),
            1_000.0 * population.diversity(m1, m2),
        );
    }
    println!();
    println!(
        "**The selfing population is {:.2} times as diverse as the outbred one.**",
        selfing.diversity(m1, m2) / outbred.diversity(m1, m2)
    );

    println!();
    println!("## What one genome from each shows, in closed form");
    println!();
    println!("| population | heterozygous | both copies non-reference |");
    println!("|---|---:|---:|");
    for population in [&selfing, &outbred] {
        println!(
            "| {} | {:.8} | {:.8} |",
            population.name,
            population.heterozygous_rate(m1, m2),
            population.homozygous_alt_rate(m1, m2),
        );
    }
    // **Held rather than printed and trusted.** If the two rows above ever stop agreeing, every
    // sentence this program exists to support is void, and a reader comparing eight decimal places
    // by eye is not a check.
    let het_gap = (selfing.heterozygous_rate(m1, m2) - outbred.heterozygous_rate(m1, m2)).abs();
    let alt_gap = (selfing.homozygous_alt_rate(m1, m2) - outbred.homozygous_alt_rate(m1, m2)).abs();
    assert!(
        het_gap < 1e-15 && alt_gap < 1e-15,
        "the two populations were built so that one genome shows the same two rates in both, and \
         they differ by {het_gap} and {alt_gap}: the construction in this file's header is wrong \
         and nothing below it means anything"
    );
    println!();
    println!(
        "The two rows agree to within {:.1e} and {:.1e}, which is floating-point noise. **A single \
         genome's census is drawn from the same distribution under both populations**, so no \
         estimator can prefer one over the other.",
        het_gap.max(1e-18),
        alt_gap.max(1e-18)
    );

    println!();
    println!("## What the shipped fit returns from one genome of each");
    println!();
    println!(
        "One genome, {positions} census positions, at three depths. The fit is the run's own \
         (`fit_jointly`), asked for nothing special."
    );
    println!();
    println!(
        "| depth | population | truth | fitted diversity | fitted / truth | fitted homozygote \
         excess | truth |"
    );
    println!("|---:|---|---:|---:|---:|---:|---:|");
    for depth in [3.0_f64, 20.0, 100.0] {
        for population in [&selfing, &outbred] {
            let truth = population.diversity(m1, m2);
            let fitted = fit_one_genome(population, positions, depth, m1);
            println!(
                "| {depth} | {} | {truth:.6} | {:.6} | {:.3}× | {:.3} | {} |",
                population.name,
                fitted.diversity,
                fitted.diversity / truth,
                fitted.homozygote_excess,
                population.inbreeding,
            );
        }
    }
    println!();
    println!(
        "**Read the two rows of each depth against each other, not against the truth.** They are \
         two draws from one distribution, so whatever the fit returns it returns for both — and \
         the two truths differ by a factor of five, so at most one of the two rows can be right."
    );

    print_cohort_sweep(&selfing, m1, m2);
}

/// **From how many samples up does the fit recover the inbreeding coefficient?**
///
/// The section above shows one genome cannot, and shows *why*: one genome's data is the same under
/// two populations five-fold apart in diversity. **That argument stops applying the moment there
/// are two genomes**, because the samples at a position share that position's frequency — how many
/// of them carry the allele says what the frequency is, and the excess of homozygotes over what
/// that frequency predicts is the coefficient. What the argument does not say is **how many samples
/// it takes before the fit actually finds it**, and the whole correction rests on that number.
///
/// So: one selfing population at `F = 0.8`, drawn at a range of cohort sizes, and what the fit
/// returns for the coefficient and for the diversity at each.
fn print_cohort_sweep(selfing: &Population, m1: f64, m2: f64) {
    let truth = selfing.diversity(m1, m2);
    println!();
    println!("## From how many samples up does the fit find the coefficient?");
    println!();
    println!(
        "The selfing population above — inbreeding coefficient **{SELFING_COEFFICIENT}**, diversity \
         {truth:.6} ({:.2} per kilobase) — drawn at each cohort size over {COHORT_SWEEP_POSITIONS} \
         census positions. **The coefficient column is what a run would use to correct its \
         diversity, and the diversity column is what the fit reports today**, which reads it off \
         the fitted curve and applies no correction at all.",
        1_000.0 * truth
    );
    println!();
    println!(
        "| depth | individuals | fitted coefficient (truth {SELFING_COEFFICIENT}) | smallest and \
         largest over the samples | fitted diversity | over the truth |"
    );
    println!("|---:|---:|---:|---|---:|---:|");
    for depth in COHORT_SWEEP_DEPTHS {
        for individuals in COHORT_SWEEP_PANELS {
            let fitted = fit_cohort(selfing, individuals, COHORT_SWEEP_POSITIONS, depth, m1);
            println!(
                "| {depth} | {individuals} | {:.3} | {:.3} to {:.3} | {:.6} | {:.3}× |",
                fitted.homozygote_excess,
                fitted.smallest_coefficient,
                fitted.largest_coefficient,
                fitted.diversity,
                fitted.diversity / truth,
            );
        }
    }
}

/// Cohort sizes for the sweep above. **One is in it on purpose**: it is the row the section before
/// proves cannot work, and a sweep that started at two would leave the reader to take that on
/// trust.
const COHORT_SWEEP_PANELS: [usize; 7] = [1, 2, 3, 5, 10, 25, 63];

/// Depths for the sweep: the shallow end this caller commits to, and one where the reads settle
/// each genotype outright.
const COHORT_SWEEP_DEPTHS: [f64; 2] = [3.0, 20.0];

/// Census positions for the sweep. **Fewer than the one-genome section's two hundred thousand**,
/// because a 63-sample fit over that many is minutes rather than seconds; at 5 positions in 100
/// segregating this still leaves 2,500 segregating positions a cohort.
const COHORT_SWEEP_POSITIONS: usize = 50_000;

struct Population {
    name: &'static str,
    inbreeding: f64,
    segregating: f64,
    fixed_alt: f64,
}

impl Population {
    /// `E[2 f (1 - f)]` over all positions — what the genotype prior wants, and what no single
    /// genome can supply.
    fn diversity(&self, m1: f64, m2: f64) -> f64 {
        self.segregating * 2.0 * (m1 - m2)
    }

    /// How often one genome from this population is heterozygous at a census position.
    fn heterozygous_rate(&self, m1: f64, m2: f64) -> f64 {
        (1.0 - self.inbreeding) * self.segregating * 2.0 * (m1 - m2)
    }

    /// …and how often both its copies are non-reference.
    fn homozygous_alt_rate(&self, m1: f64, m2: f64) -> f64 {
        self.fixed_alt + self.segregating * m2 + self.inbreeding * self.segregating * (m1 - m2)
    }

    fn invariant(&self) -> f64 {
        1.0 - self.segregating - self.fixed_alt
    }
}

struct Fitted {
    diversity: f64,
    /// The mean over the cohort's samples — the panel coefficient the correction takes.
    homozygote_excess: f64,
    smallest_coefficient: f64,
    largest_coefficient: f64,
}

fn fit_one_genome(population: &Population, positions: usize, depth: f64, m1: f64) -> Fitted {
    fit_cohort(population, 1, positions, depth, m1)
}

/// Draw `individuals` genomes from a population, record them as a census, and fit.
fn fit_cohort(
    population: &Population,
    individuals: usize,
    positions: usize,
    depth: f64,
    m1: f64,
) -> Fitted {
    // **Mixed from the name's bytes and not its length, and that is a defect a review caught.**
    // `"selfing"` and `"outbred"` are both seven bytes, so keying on the length gave the two
    // populations the identical generator state at every depth. They are meant to be independent
    // draws — the whole reading of the table below is that two independent draws from the same
    // distribution come back the same — so a shared stream would have made that reading circular.
    let mut name_key = 0xCBF2_9CE4_8422_2325_u64;
    for byte in population.name.bytes() {
        name_key = (name_key ^ u64::from(byte)).wrapping_mul(0x1000_0000_01B3);
    }
    // **The cohort size is mixed in too**, so the panel-size sweep's arms are independent draws
    // rather than nested prefixes of one. Nesting is right where the question is *what changes with
    // panel size on fixed positions* (`ng_prior_moment_estimators.rs`); here the question is what
    // the fit returns at each size, and each size should get its own cohort.
    let seed = 0xA076_1D64_78BD_642F
        ^ (depth as u64).wrapping_mul(0x9E37_79B9)
        ^ name_key
        ^ (individuals as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let mut rng = Rng(seed | 1);
    let edges = DepthBinEdges::for_census();
    let mut codes: Vec<PackedDepthCodes> = (0..individuals)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); individuals];

    for index in 0..positions {
        let frequency = match rng.pick(&[
            population.invariant(),
            population.fixed_alt,
            population.segregating,
        ]) {
            0 => 0.0,
            1 => 1.0,
            _ => rng.beta(BETA_A, BETA_B),
        };
        let allele = 1 + ((rng.uniform() * 3.0) as usize).min(2);
        // **Every sample at a position shares that position's frequency**, which is the whole
        // reason a cohort can separate the inbreeding coefficient from the diversity and one
        // genome cannot: how many of the samples carry the allele says what the frequency is,
        // and the excess of homozygotes over what that frequency predicts is the coefficient.
        for sample in 0..individuals {
            let reads = rng.poisson(depth).min(CENSUS_DEPTH_CAP);
            let genotype = if rng.uniform() < population.inbreeding {
                if rng.uniform() < frequency { 2 } else { 0 }
            } else {
                u32::from(rng.uniform() < frequency) + u32::from(rng.uniform() < frequency)
            };
            let carried = f64::from(genotype) / 2.0;
            let on_candidate =
                carried * (1.0 - CLEAN_ERROR_RATE) + (1.0 - carried) * CLEAN_ERROR_RATE / 3.0;
            let on_reference =
                (1.0 - carried) * (1.0 - CLEAN_ERROR_RATE) + carried * CLEAN_ERROR_RATE / 3.0;
            let mut counts = [0_u32; 5];
            for _ in 0..reads {
                let u = rng.uniform();
                let code = if u < on_candidate {
                    allele
                } else if u < on_candidate + on_reference {
                    0
                } else if u < on_candidate + on_reference + CLEAN_ERROR_RATE / 3.0 {
                    (allele % 3) + 1
                } else {
                    ((allele + 1) % 3) + 1
                };
                counts[code] += 1;
            }
            codes[sample].set(index, DepthCode::Binned(edges.bin_for(reads)));
            for (code, count) in counts.iter().enumerate() {
                if code == 0 || *count == 0 {
                    continue;
                }
                sparse[sample].push(AlleleObservation {
                    index: index as u32,
                    allele: match code {
                        1 => ObservedAllele::C,
                        2 => ObservedAllele::G,
                        _ => ObservedAllele::T,
                    },
                    reads: u8::try_from(*count).expect("a capped count fits one byte"),
                });
            }
        }
    }
    let _ = m1;

    let terms = RecordingTerms {
        selection: SelectionTermsDigest::of(&SelectionTerms {
            seed,
            reference: ReferenceDigest([7; 16]),
            analysed_regions: RegionSetDigest([9; 16]),
            catalog_built_under: CatalogBuildSettings {
                criteria: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "0.1.0".to_string(),
            },
            ssr_criteria: StrRepeatCriteria::default(),
            generic_target: positions as u64,
            ssr_cap: 1_000,
        }),
        kept_loci: CensusLociDigester::new().finish(),
        ssr_stratum_counts: Default::default(),
        read_cap: ReadCap(1_000),
        depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
        depth_cap: DepthCap::new(CENSUS_DEPTH_CAP),
    };
    let records: Vec<SampleCensusEvidence> = (0..individuals)
        .map(|sample| {
            SampleCensusEvidence::resident(
                format!("s{sample:03}"),
                terms.clone(),
                NamedReadGroup::drawn_for(&format!("s{sample:03}"), [ReadGroupId(sample as u32)]),
                BTreeMap::from([(
                    SectionKey::Generic(ReadGroupId(sample as u32)),
                    Section::Generic(GenericEvidence::from_parts(
                        std::mem::replace(&mut codes[sample], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[sample]),
                    )),
                )]),
            )
        })
        .collect();
    let mut evidence = CohortCensusEvidence::new(records).expect("a drawn cohort records one way");
    let config = JointFitConfig {
        quadrature_nodes: 12,
        starting_points: StartingPoint::spanning_the_class_separation(),
        ..JointFitConfig::default()
    };
    let fit = fit_jointly(&mut evidence, &config).expect("one genome pools");
    let coefficients: Vec<f64> = fit
        .hom_excess
        .values()
        .map(|estimate| estimate.value.get())
        .collect();
    Fitted {
        diversity: fit.expected_heterozygosity,
        homozygote_excess: coefficients.iter().sum::<f64>() / coefficients.len().max(1) as f64,
        smallest_coefficient: coefficients.iter().copied().fold(f64::INFINITY, f64::min),
        largest_coefficient: coefficients
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max),
    }
}

/// The stream every drawn number comes from — the same xorshift the sibling harnesses use.
struct Rng(u64);

impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }

    fn pick(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut cut = self.uniform() * total;
        for (index, weight) in weights.iter().enumerate() {
            cut -= weight;
            if cut <= 0.0 {
                return index;
            }
        }
        weights.len() - 1
    }

    fn beta(&mut self, a: f64, b: f64) -> f64 {
        let x = self.gamma(a);
        let y = self.gamma(b);
        (x / (x + y)).clamp(1e-9, 1.0 - 1e-9)
    }

    fn gamma(&mut self, shape: f64) -> f64 {
        if shape < 1.0 {
            let u = self.uniform().max(1e-12);
            return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
        }
        let d = shape - 1.0 / 3.0;
        let c = 1.0 / (9.0 * d).sqrt();
        loop {
            let z = self.normal();
            let v = (1.0 + c * z).powi(3);
            if v <= 0.0 {
                continue;
            }
            let u = self.uniform().max(1e-12);
            if u.ln() < 0.5 * z * z + d - d * v + d * v.ln() {
                return d * v;
            }
        }
    }

    fn normal(&mut self) -> f64 {
        let u1 = self.uniform().max(1e-12);
        let u2 = self.uniform();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    fn poisson(&mut self, mean: f64) -> u32 {
        let limit = (-mean).exp();
        let mut product = self.uniform();
        let mut count = 0;
        while product > limit && count < 400 {
            count += 1;
            product *= self.uniform();
        }
        count
    }
}
