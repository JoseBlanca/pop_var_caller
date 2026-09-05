//! Where does having many samples at one position start paying?
//!
//! The whole argument for ng's second route to the calling parameters is that **several
//! samples at one position separate a mismapped position from a heterozygous one**. A
//! mismapped position — two stretches of genome that the reference holds once, so both
//! stretches' reads pile up at one place — reads part non-reference in *every* sample. A real
//! heterozygote reads part non-reference only in the samples that carry the allele. One sample
//! at that position cannot tell the two apart and the estimator falls back to a blind share:
//! *some* fraction of positions are mismapped, and every position pays a slice of it.
//!
//! **How many samples it takes before that stops being a slogan is the number a user needs,
//! and nothing states it.** `parameter_prepass_joint_fit.md` §8's third measurement asks for
//! it as a curve in sample count rather than a single ratio, and this program is that curve.
//!
//! # What it does
//!
//! One cohort is drawn at known parameters, with a share of positions planted as mismapped —
//! at those positions every sample's reads disagree with the reference far more often than the
//! sequencer's error rate allows. The same drawn positions are then refitted using the first
//! 2, 3, 5, 10, 25 and 50 samples, and what is reported is **the fitted heterozygosity against
//! the heterozygosity the samples were actually drawn with**.
//!
//! **Holding the positions fixed across the arms is what makes it a curve in sample count.**
//! Redrawing them would let the answer move for a reason that has nothing to do with the
//! mechanism being measured.
//!
//! # The control that says the curve means something
//!
//! Every arm is run twice: once on the cohort as drawn, and once on **the same cohort with no
//! mismapped positions in it at all**, everything else identical. If the excess falls with
//! sample count in the first and there is no excess in the second, the curve is measuring the
//! mechanism. If the second shows an excess too, the curve is measuring something else — a
//! shallow-cohort bias in the estimator — and would have been reported as the mechanism.
//!
//! ```text
//! ng_joint_sample_count_sweep [positions] [depth] [noisy-share] [het-per-kb] [noisy-rate]
//! ```

use std::collections::BTreeMap;
use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, NamedReadGroup, ObservedAllele, PackedDepthCodes, ReadCap, RecordingTerms,
    SampleCensusEvidence, Section, SectionKey, SelectionTermsDigest,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFitConfig, StartingPoint, fit_jointly,
};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::ReadGroupId;

/// The sample counts the curve is reported at.
const ARMS: [usize; 6] = [2, 3, 5, 10, 25, 50];

/// How often a read misreads a base at an ordinary position.
const CLEAN: f64 = 0.002;

/// How inbred each drawn plant is — a departure from random mating, which the fit also has to
/// find and which trades against the frequency density.
const HOM_EXCESS: f64 = 0.15;

fn main() {
    let mut args = std::env::args().skip(1);
    let positions: usize = args.next().map_or(60_000, |a| a.parse().expect("a count"));
    let depth: f64 = args.next().map_or(8.0, |a| a.parse().expect("a depth"));
    let noisy_share: f64 = args.next().map_or(0.01, |a| a.parse().expect("a share"));
    // **Heterozygosity is the axis the mechanism lives on**, so it is an input rather than a
    // consequence of the density's masses: at tomato's floor of 0.149 per kilobase at least 97
    // positions in 100 carrying an alternative read are artefact, and at 12 per kilobase
    // almost none are (`parameter_prepass_joint_fit.md` §3.2).
    let het_per_kb: f64 = args.next().map_or(12.0, |a| a.parse().expect("a rate"));
    let noisy: f64 = args.next().map_or(0.06, |a| a.parse().expect("a rate"));

    // A Beta(0.7, 2.5) population is heterozygous at 2ab/((a+b)(a+b+1)) of its segregating
    // positions, and inbreeding takes another slice, so the segregating share is what is left
    // when the wanted heterozygosity is divided by both.
    let (a, b) = (0.7_f64, 2.5_f64);
    let per_segregating = 2.0 * a * b / ((a + b) * (a + b + 1.0)) * (1.0 - HOM_EXCESS);
    let p_fixed_alt = 0.0005;
    let density = FrequencyDensity {
        p_invariant: 1.0 - p_fixed_alt - (het_per_kb / 1_000.0) / per_segregating,
        p_fixed_alt,
        a,
        b,
    };

    println!("positions        {positions}");
    println!("depth            {depth} reads a position");
    println!(
        "mismapped        {:.3} of positions, where a read disagrees at {noisy} instead of \
         {CLEAN}",
        noisy_share
    );
    println!(
        "segregating      {:.5} of positions, for a drawn heterozygosity near {het_per_kb} per \
         kilobase",
        1.0 - density.p_invariant - density.p_fixed_alt
    );
    println!("inbreeding       {HOM_EXCESS} homozygote excess in every drawn plant");

    for (label, planted) in [
        ("with mismapped positions", noisy_share),
        ("none planted", 0.0),
    ] {
        println!("\n=== {label} ===");
        let cohort = draw(
            ARMS[ARMS.len() - 1],
            positions,
            depth,
            planted,
            noisy,
            density,
            0x9E37_79B9_7F4A_7C15,
        );
        println!(
            "  drawn heterozygosity {:.4} per position ({:.3} per kilobase)",
            cohort.heterozygosity,
            1_000.0 * cohort.heterozygosity
        );
        println!(
            "\n  {:>8}{:>14}{:>14}{:>10}{:>14}{:>12}",
            "samples", "het drawn", "het fitted", "ratio", "error rate", "mismapped"
        );
        for &arm in &ARMS {
            let at = Instant::now();
            // The first `arm` samples as their own cohort. **Cloned**, because the drawn
            // positions are held fixed across the arms and each arm fits its own prefix.
            let subset = cohort.samples[..arm].to_vec();
            let mut fitted_cohort =
                CohortCensusEvidence::new(subset.clone()).expect("a drawn cohort records one way");
            let config = JointFitConfig {
                quadrature_nodes: 12,
                starting_points: StartingPoint::spanning_the_class_separation(),
                ..JointFitConfig::default()
            };
            let fit = fit_jointly(&mut fitted_cohort, &config).expect("a drawn cohort pools");
            let drawn: f64 = cohort.heterozygous_per_sample[..arm].iter().sum::<f64>() / arm as f64;
            let fitted: f64 = subset
                .iter()
                .map(|sample| fit.rates[&sample.sample].value.heterozygous)
                .sum::<f64>()
                / arm as f64;
            println!(
                "  {:>8}{:>14.5}{:>14.5}{:>9.2}×{:>14.5}{:>12.4}   {}{:.1} s",
                arm,
                drawn,
                fitted,
                fitted / drawn,
                // The mean over the libraries: one a sample, since a library belongs to one plant.
                fit.noise.values().map(|n| n.value.clean).sum::<f64>() / fit.noise.len() as f64,
                fit.noisy_share,
                if fit.converged { "" } else { "RAN OUT, " },
                at.elapsed().as_secs_f64()
            );
        }
    }
}

struct Drawn {
    samples: Vec<SampleCensusEvidence>,
    /// The share of positions at which each sample was actually drawn heterozygous.
    heterozygous_per_sample: Vec<f64>,
    heterozygosity: f64,
}

/// Draw a cohort at known parameters and write it into records the fit will read.
///
/// **The reference base is code 0 by construction**, so the three candidate alleles are codes
/// 1 to 3 and the fit's own sum over which one segregates is being asked to find the right one.
#[allow(
    clippy::too_many_arguments,
    reason = "the drawn cohort's own parameters"
)]
fn draw(
    samples: usize,
    positions: usize,
    depth: f64,
    noisy_share: f64,
    noisy: f64,
    density: FrequencyDensity,
    seed: u64,
) -> Drawn {
    let mut rng = Rng(seed);
    let edges = DepthBinEdges::for_census();
    let mut codes: Vec<PackedDepthCodes> = (0..samples)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
    let mut heterozygous = vec![0_u64; samples];

    for index in 0..positions {
        let rate = if rng.uniform() < noisy_share {
            noisy
        } else {
            CLEAN
        };
        let branch = rng.pick(&[
            density.p_invariant,
            density.p_fixed_alt,
            1.0 - density.p_invariant - density.p_fixed_alt,
        ]);
        let allele = 1 + ((rng.uniform() * 3.0) as usize).min(2);
        let f = match branch {
            0 => 0.0,
            1 => 1.0,
            _ => rng.beta(density.a, density.b),
        };
        for sample in 0..samples {
            let reads = rng.poisson(depth);
            let genotype = match branch {
                0 => 0_usize,
                1 => 2,
                _ => {
                    let het = 2.0 * f * (1.0 - f) * (1.0 - HOM_EXCESS);
                    let shift = HOM_EXCESS * f * (1.0 - f);
                    rng.pick(&[(1.0 - f) * (1.0 - f) + shift, het, f * f + shift])
                }
            };
            if genotype == 1 {
                heterozygous[sample] += 1;
            }
            let carried = genotype as f64 / 2.0;
            let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
            let on_reference = (1.0 - carried) * (1.0 - rate) + carried * rate / 3.0;
            let mut counts = [0_u32; 5];
            for _ in 0..reads {
                let u = rng.uniform();
                let code = if u < on_candidate {
                    allele
                } else if u < on_candidate + on_reference {
                    0
                } else if u < on_candidate + on_reference + rate / 3.0 {
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
                    reads: u8::try_from(*count)
                        .expect("a drawn count fits the census's one-byte field"),
                });
            }
        }
    }

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
        depth_cap: DepthCap::MAX,
    };
    let records = (0..samples)
        .map(|s| {
            SampleCensusEvidence::resident(
                format!("s{s:02}"),
                terms.clone(),
                NamedReadGroup::drawn_for(&format!("s{s:02}"), [ReadGroupId(s as u32)]),
                BTreeMap::from([(
                    SectionKey::Generic(ReadGroupId(s as u32)),
                    Section::Generic(GenericEvidence::from_parts(
                        std::mem::replace(&mut codes[s], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[s]),
                    )),
                )]),
            )
        })
        .collect();
    let per_sample: Vec<f64> = heterozygous
        .iter()
        .map(|n| *n as f64 / positions as f64)
        .collect();
    let mean = per_sample.iter().sum::<f64>() / per_sample.len() as f64;
    Drawn {
        samples: records,
        heterozygous_per_sample: per_sample,
        heterozygosity: mean,
    }
}

/// The stream every drawn number comes from. Deterministic, so a run is reproducible.
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
