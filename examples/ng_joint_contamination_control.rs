//! **Does excluding the mismapped positions leave a real 3% contamination still findable?**
//!
//! On the 63 tomato accessions the contamination estimate came back with a median of 6.5% and a
//! highest of 12.5%, which cannot be true of sixty-three plants from a public archive: it is a
//! floor, not a measurement. The named cause is that **a position where two stretches of genome
//! the reference holds once both pile reads up produces the contamination signature in every
//! sample at once** — a small share of reads carrying an allele the sample should not have —
//! and nothing excluded such positions from the ones contamination was measured over.
//!
//! The fit already computes, for every position, the probability that it is one of those. This
//! program is the control that says whether using it works, and it needs both halves:
//!
//! - **the floor** — a drawn panel *with* mismapped positions planted and nobody contaminated,
//!   which must come back near zero once they are excluded and does not today;
//! - **the signal** — the same panel with one sample genuinely contaminated at 3%, which must
//!   still be found. An estimator that has been made blind returns zero on both, and that is
//!   the trap this project has already fallen into once.
//!
//! Everything runs through the library: the cohort is drawn as the records a walk would have
//! written, [`fit_jointly`] fits it and produces the per-position probabilities, and
//! [`fit_contamination`] is then called once per arm on the same fit. So what is compared is
//! the shipped code under different settings, not a re-implementation of it.
//!
//! # The arms
//!
//! Two settings crossed. **Which positions are markers**: all of them, as before, against
//! dropping those the fit says are more likely mismapped than not. And **where the sample is
//! taken to stand on the panel's axes of variation** while its fraction is searched for: as the
//! decomposition read it off the sample's own reads, or divided by `1 − α` — which undoes
//! exactly the drag a fraction `α` of stray reads causes — or each axis searched freely beside
//! `α`.
//!
//! ```text
//! ng_joint_contamination_control [positions] [samples] [depth] [groups] [fst] [noisy-share] [noisy-rate] [alpha]
//! ```

use std::collections::BTreeMap;
use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, DepthCap, DepthCode, DepthLadderDigest, GenericEvidence, ObservedAllele,
    PackedDepthCodes, ReadCap, RecordingTerms, SampleCensusEvidence,
};
use pop_var_caller::ng::parameter_estimation::joint::contamination::{
    ContaminationConfig, ContaminationEstimate, OwnCoordinates, fit_contamination,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFit, JointFitConfig, fit_jointly,
};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::ReadGroupId;

/// How often a read misreads a base at an ordinary position.
const CLEAN: f64 = 0.003;

/// How much less heterozygous than random mating each drawn plant is — a selfing crop.
const HOM_EXCESS: f64 = 0.6;

fn main() {
    let mut args = std::env::args().skip(1);
    let positions: usize = args.next().map_or(400_000, |a| a.parse().expect("a count"));
    let samples: usize = args.next().map_or(40, |a| a.parse().expect("a count"));
    let depth: f64 = args.next().map_or(3.0, |a| a.parse().expect("a depth"));
    let groups: usize = args.next().map_or(4, |a| a.parse().expect("a count"));
    let fst: f64 = args.next().map_or(0.20, |a| a.parse().expect("an F_st"));
    let noisy_share: f64 = args.next().map_or(0.033, |a| a.parse().expect("a share"));
    let noisy: f64 = args.next().map_or(0.024, |a| a.parse().expect("a rate"));
    let alpha: f64 = args.next().map_or(0.03, |a| a.parse().expect("a fraction"));
    let components: usize = args.next().map_or(4, |a| a.parse().expect("a count"));

    println!("positions        {positions}");
    println!("samples          {samples} in {groups} subpopulations at F_st {fst}");
    println!("depth            {depth} reads a position");
    println!(
        "mismapped        {noisy_share} of positions, where a read disagrees at {noisy} instead \
         of {CLEAN}"
    );
    println!("contaminated     sample 0, at {alpha} of its reads, by a plant drawn from the panel");
    println!(
        "\nBoth halves are needed. A panel with nobody contaminated says whether the floor is\n\
         gone; the same panel with one sample at {alpha} says whether the estimator can still\n\
         see anything at all."
    );

    let mut scenarios = vec![
        (noisy_share, alpha),
        (noisy_share, 0.0),
        (0.0, alpha),
        (0.0, 0.0),
    ];
    scenarios.dedup();
    for (planted, spike) in scenarios {
        let label = format!(
            "{}, {}",
            if planted > 0.0 {
                format!("{planted} of positions mismapped")
            } else {
                "no mismapped positions at all".to_string()
            },
            if spike > 0.0 {
                format!("sample 0 contaminated at {spike}")
            } else {
                "nobody contaminated".to_string()
            }
        );
        println!("\n=== {label} ===");
        let at = Instant::now();
        let drawn = draw(
            samples,
            positions,
            depth,
            planted,
            noisy,
            groups,
            fst,
            spike,
            0x9E37_79B9_7F4A_7C15,
        );
        let config = JointFitConfig {
            quadrature_nodes: 12,
            // **The third class is off here**, and this harness is not the place to judge it:
            // it is measured in `ng_joint_duplicated_in_fit.rs`. Left on, every fit below pays
            // for a second frequency integral and the numbers move for a reason that has
            // nothing to do with contamination.
            duplicated_positions: false,
            ..JointFitConfig::default()
        };
        let fit = fit_jointly(&drawn.samples, &config).expect("a drawn cohort pools");
        println!(
            "  the fit: {} of positions booked mismapped against {planted} planted, error rate \
             {:.5} against {CLEAN}, {} passes, {:.0} s",
            format_args!("{:.4}", fit.noisy_share),
            fit.noise[&ReadGroupId(0)].value.clean,
            fit.passes,
            at.elapsed().as_secs_f64()
        );
        let condemned = fit.noisy_posterior.iter().filter(|p| **p > 0.5).count();
        println!(
            "  {condemned} of {} positions are more likely mismapped than not",
            fit.noisy_posterior.len()
        );

        println!(
            "\n  {:<50}{:>10}{:>10}{:>10}{:>10}",
            "", "markers", "sample 0", "median", "worst other"
        );
        for keep_mismapped in [true, false] {
            for (own, integrate) in [
                (OwnCoordinates::AsRead, false),
                (OwnCoordinates::AsRead, true),
                (OwnCoordinates::UndoneByAlpha, true),
                (OwnCoordinates::MaximisedFreely, true),
            ] {
                let settings = ContaminationConfig {
                    components,
                    max_noisy_posterior: if keep_mismapped { 1.0 } else { 0.5 },
                    weight_by_posterior: !keep_mismapped,
                    own_coordinates: own,
                    integrate_over_depth_bin: integrate,
                };
                let (markers, alphas) = run(&drawn.samples, &fit, &settings);
                let mut others: Vec<f64> = alphas[1..].to_vec();
                others.sort_by(f64::total_cmp);
                println!(
                    "  {:<50}{markers:>10}{:>10.4}{:>10.4}{:>10.4}",
                    format!(
                        "{}, {}{}",
                        if keep_mismapped {
                            "every position"
                        } else {
                            "mismapped dropped"
                        },
                        match own {
                            OwnCoordinates::AsRead => "coordinates as read",
                            OwnCoordinates::UndoneByAlpha => "coordinates undone by α",
                            OwnCoordinates::MaximisedFreely => "coordinates free",
                        },
                        if integrate { ", depth summed over" } else { "" }
                    ),
                    alphas[0],
                    others[others.len() / 2],
                    others.last().copied().unwrap_or(0.0),
                );
            }
        }
    }
}

/// One arm: how many markers survived, and every sample's fraction in order.
fn run(
    samples: &[SampleCensusEvidence],
    fit: &JointFit,
    settings: &ContaminationConfig,
) -> (u64, Vec<f64>) {
    let error: Vec<f64> = samples
        .iter()
        .map(|_| fit.noise[&ReadGroupId(0)].value.clean)
        .collect();
    let excess: Vec<f64> = samples
        .iter()
        .map(|sample| fit.hom_excess[&sample.sample].value.get())
        .collect();
    let estimates = fit_contamination(
        samples,
        &DepthBinEdges::new(),
        &error,
        &excess,
        &fit.noisy_posterior,
        settings,
    );
    let markers = estimates
        .iter()
        .find_map(|estimate| match estimate {
            ContaminationEstimate::Estimated { markers, .. } => Some(*markers),
            ContaminationEstimate::NotIdentified { .. } => None,
        })
        .unwrap_or(0);
    let alphas = estimates
        .iter()
        .map(|estimate| match estimate {
            ContaminationEstimate::Estimated { alpha, .. } => *alpha,
            ContaminationEstimate::NotIdentified { .. } => f64::NAN,
        })
        .collect();
    (markers, alphas)
}

struct Drawn {
    samples: Vec<SampleCensusEvidence>,
}

/// Draw a structured cohort and write it into the records the fit reads.
///
/// **The contaminant is drawn from the whole panel**, not from the contaminated sample's own
/// subpopulation, because a neighbouring library on a plate is whoever else was on the plate.
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
    groups: usize,
    fst: f64,
    alpha: f64,
    seed: u64,
) -> Drawn {
    // A Beta(0.7, 2.5) population, and the same three branches the fit's own density has: most
    // positions carry one allele, a few are fixed on a non-reference one, the rest segregate.
    let (a, b) = (0.7_f64, 2.5_f64);
    let density = FrequencyDensity {
        p_invariant: 0.97,
        p_fixed_alt: 0.0005,
        a,
        b,
    };
    let mut rng = Rng(seed);
    let edges = DepthBinEdges::new();
    let mut codes: Vec<PackedDepthCodes> = (0..samples)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
    let group_of = |s: usize| s * groups / samples;

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
        // Balding–Nichols: one ancestral frequency, then a frequency per subpopulation drawn
        // around it with a spread set by F_st.
        let ancestral = match branch {
            0 => 0.0,
            1 => 1.0,
            _ => rng.beta(a, b),
        };
        let per_group: Vec<f64> = (0..groups)
            .map(|_| {
                if branch != 2 || fst <= 0.0 {
                    ancestral
                } else {
                    let scale = (1.0 - fst) / fst;
                    rng.beta(
                        (ancestral * scale).max(1e-3),
                        ((1.0 - ancestral) * scale).max(1e-3),
                    )
                }
            })
            .collect();

        let genotype = |rng: &mut Rng, f: f64| -> usize {
            let het = 2.0 * f * (1.0 - f) * (1.0 - HOM_EXCESS);
            let shift = HOM_EXCESS * f * (1.0 - f);
            rng.pick(&[(1.0 - f) * (1.0 - f) + shift, het, f * f + shift])
        };
        let genotypes: Vec<usize> = (0..samples)
            .map(|s| match branch {
                0 => 0,
                1 => 2,
                _ => genotype(&mut rng, per_group[group_of(s)]),
            })
            .collect();

        for sample in 0..samples {
            let reads = rng.poisson(depth);
            let mut counts = [0_u32; 5];
            for _ in 0..reads {
                let from = if sample == 0 && alpha > 0.0 && rng.uniform() < alpha {
                    genotypes[(rng.uniform() * samples as f64) as usize % samples]
                } else {
                    genotypes[sample]
                };
                let carried = from as f64 / 2.0;
                let on_candidate = carried * (1.0 - rate) + (1.0 - carried) * rate / 3.0;
                let on_reference = (1.0 - carried) * (1.0 - rate) + carried * rate / 3.0;
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
                    reads: *count,
                });
            }
        }
    }

    let terms = RecordingTerms {
        selection: SelectionTerms {
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
        },
        kept_loci: CensusLociDigester::new().finish(),
        ssr_stratum_counts: Default::default(),
        read_cap: ReadCap(1_000),
        depth_ladder: DepthLadderDigest::of(&DepthBinEdges::new()),
        depth_cap: DepthCap(124),
        coverage_window: None,
    };
    Drawn {
        samples: (0..samples)
            .map(|s| SampleCensusEvidence {
                sample: format!("s{s:02}"),
                generic: [(
                    ReadGroupId(0),
                    GenericEvidence::from_parts(
                        std::mem::replace(&mut codes[s], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[s]),
                    ),
                )]
                .into_iter()
                .collect(),
                ssr: BTreeMap::new(),
                coverage: None,
                terms: terms.clone(),
            })
            .collect(),
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
