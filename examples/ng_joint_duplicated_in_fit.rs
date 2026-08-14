//! **Does the third class, now inside the estimator, do what the drawn measurement said it would?**
//!
//! A position where a plant carries two copies of a stretch the reference holds once reads about
//! half non-reference in that plant and not in the others. The two-class model has nowhere to put
//! it but *heterozygous*, and
//! [`duplicated_class_identification_2026-08-13.md`](../doc/devel/ng/reports/duplicated_class_identification_2026-08-13.md)
//! measured what that costs on a purpose-built program: on a fifty-sample selfing panel at three
//! reads a position, observed heterozygosity **50.6% above the truth**, expected heterozygosity
//! only 10.6% above it, and the fitted homozygote excess **0.4471 where the truth is 0.5942**.
//!
//! That program is not the estimator. This one puts the same cohort through
//! `parameter_estimation::joint::fit`, which now carries the class, and asks three things:
//!
//! 1. **Does the shipped estimator reproduce the cost of ignoring the class?** If it does not,
//!    either the drawing here or the implementation disagrees with that measurement, and either
//!    way the number cannot be carried over.
//! 2. **Does the class recover it, identified from the cohort alone?** The evidence is the
//!    absence of samples homozygous for the non-reference allele, which needs a panel rather
//!    than depth.
//! 3. **What does the class cost on a cohort that has none?** This is the control, and it is the
//!    one that decides whether the class can be on by default: a class that invents duplications
//!    where there are none buys the first two answers at the price of every other cohort.
//!
//! Each panel size is fitted three ways: with the class off, with it on and identified from the
//! cohort pattern, and with it on and handed a coverage reading per sample per position.
//!
//! **The coverage arm is an upper bound rather than a method.** It is handed readings drawn from
//! exactly the scatter model it scores them with, and each position's reading is independent of
//! its neighbours', where a real window shares one reading across every position inside it.
//!
//! ```text
//! ng_joint_duplicated_in_fit [positions] [depth] [inbreeding] [aligned-bases-a-window]
//! ```

use std::collections::BTreeMap;
use std::sync::Arc;
use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, ObservedAllele, PackedDepthCodes, ReadCap, RecordingTerms,
    SampleCensusEvidence, Section, SectionKey,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{
    FrequencyDensity, JointFitConfig, fit_jointly,
};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::ReadGroupId;

/// The panel sizes the answer is reported at. **The cohort pattern needs samples, not reads**, so
/// this is the axis it lives on.
const PANELS: [usize; 3] = [10, 25, 50];

/// How often a read misreads a base at an ordinary position, and at a mismapped one, and how many
/// positions are mismapped.
const CLEAN: f64 = 0.0019;
const NOISY: f64 = 0.053;
const NOISY_SHARE: f64 = 0.01;

/// **How many duplicated carrier positions a sample has, per genuinely heterozygous position** —
/// 150 to 590 per two million against about 668 real ones, on eight tomato samples. The bias this
/// causes is a ratio, so the ratio is what the drawn panel holds fixed.
const DUPLICATED_PER_HETEROZYGOUS: f64 = 1.0 / 3.0;

/// The share of duplicated positions the whole panel carries — the reference's own collapse.
const COLLAPSE_SHARE: f64 = 0.090;

/// The carrier frequency of the rest, fitted to the same eight-sample counts.
const CARRIER: (f64, f64) = (1.1895, 9.5477);

fn main() {
    let mut args = std::env::args().skip(1);
    let positions: usize = args.next().map_or(100_000, |a| a.parse().expect("a count"));
    let depth: f64 = args.next().map_or(3.0, |a| a.parse().expect("a depth"));
    let inbreeding: f64 = args
        .next()
        .map_or(0.6, |a| a.parse().expect("a coefficient"));
    let aligned: f64 = args
        .next()
        .map_or(12_000.0, |a| a.parse().expect("a count"));

    let (a, b) = (0.3_f64, 1.2_f64);
    let density = FrequencyDensity {
        p_invariant: 0.990,
        p_fixed_alt: 0.002,
        a,
        b,
    };

    println!("positions        {positions}");
    println!("depth            {depth} reads a position");
    println!("inbreeding       {inbreeding} homozygote excess in every drawn plant");
    println!(
        "duplicated       one carrier position for every {:.1} heterozygous ones, {COLLAPSE_SHARE} \
         of them carried by the whole panel",
        1.0 / DUPLICATED_PER_HETEROZYGOUS
    );
    println!("coverage         a window collecting {aligned} aligned bases");

    for samples in PANELS {
        for planted in [true, false] {
            let drawn = draw(
                samples,
                positions,
                depth,
                inbreeding,
                density,
                planted,
                aligned,
                0x9E37_79B9_7F4A_7C15,
            );
            println!(
                "\n=== {samples} samples, {} ===",
                if planted {
                    format!(
                        "{} duplicated positions somebody carries, {} carrier pairs",
                        drawn.duplicated_positions, drawn.carrier_pairs
                    )
                } else {
                    "no duplicated positions at all".to_string()
                }
            );
            println!(
                "  drawn: heterozygosity {:.5}, expected heterozygosity {:.5}, homozygote excess \
                 {inbreeding}",
                drawn.heterozygosity, drawn.expected_heterozygosity
            );
            println!(
                "\n  {:<26}{:>12}{:>9}{:>12}{:>9}{:>10}{:>12}",
                "", "het", "", "expected", "", "excess", "class weight"
            );
            for arm in ["the class off", "cohort pattern", "pattern and coverage"] {
                let at = Instant::now();
                let config = JointFitConfig {
                    quadrature_nodes: 12,
                    max_passes: 120,
                    duplicated_positions: arm != "the class off",
                    coverage_odds: if arm == "pattern and coverage" {
                        drawn.coverage_odds.clone()
                    } else {
                        Vec::new()
                    },
                    ..JointFitConfig::default()
                };
                let mut cohort = CohortCensusEvidence::new(drawn.samples.clone())
                    .expect("a drawn cohort records one way");
                let fit = fit_jointly(&mut cohort, &config).expect("a drawn cohort pools");
                let het: f64 = drawn
                    .samples
                    .iter()
                    .map(|sample| fit.rates[&sample.sample].value.heterozygous)
                    .sum::<f64>()
                    / samples as f64;
                let excess: f64 = drawn
                    .samples
                    .iter()
                    .map(|sample| fit.hom_excess[&sample.sample].value.get())
                    .sum::<f64>()
                    / samples as f64;
                println!(
                    "  {arm:<26}{het:>12.5}{:>8.1}%{:>12.5}{:>8.1}%{excess:>10.4}{:>12.5}   {}{:.0} s",
                    100.0 * (het / drawn.heterozygosity - 1.0),
                    fit.expected_heterozygosity,
                    100.0 * (fit.expected_heterozygosity / drawn.expected_heterozygosity - 1.0),
                    fit.duplicated.as_ref().map_or(0.0, |d| d.value.share),
                    if fit.converged { "" } else { "RAN OUT, " },
                    at.elapsed().as_secs_f64()
                );
            }
        }
    }
}

struct Drawn {
    samples: Vec<SampleCensusEvidence>,
    /// Per sample, `ln P(coverage | two copies) − ln P(coverage | one)` at every position.
    coverage_odds: Vec<Arc<[f32]>>,
    /// The share of positions at which a sample was drawn genuinely heterozygous, averaged.
    heterozygosity: f64,
    /// What random mating in the drawn population would give.
    expected_heterozygosity: f64,
    duplicated_positions: usize,
    /// How many (position, sample) pairs really carry a duplication.
    carrier_pairs: usize,
}

/// The scatter a window's relative coverage reading really has, from the tomato measurement:
/// `√(copies² × 0.194² + copies × 313 / aligned bases)`.
fn coverage_sd(copies: f64, aligned: f64) -> f64 {
    (copies * copies * 0.194 * 0.194 + copies * 313.0 / aligned).sqrt()
}

#[allow(
    clippy::too_many_arguments,
    reason = "the drawn cohort's own parameters"
)]
fn draw(
    samples: usize,
    positions: usize,
    depth: f64,
    inbreeding: f64,
    density: FrequencyDensity,
    planted: bool,
    aligned: f64,
    seed: u64,
) -> Drawn {
    let mut rng = Rng(seed);
    let edges = DepthBinEdges::new();
    let mut codes: Vec<PackedDepthCodes> = (0..samples)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
    let mut odds: Vec<Vec<f32>> = vec![vec![0.0; positions]; samples];
    let mut heterozygous = vec![0_u64; samples];
    let mut expected = 0.0_f64;
    let (mut duplicated_positions, mut carrier_pairs) = (0_usize, 0_usize);

    // The share of positions drawn from the duplicated class, so that a sample ends up with one
    // carrier position for every three heterozygous ones.
    let per_segregating = 2.0 * density.a * density.b
        / ((density.a + density.b) * (density.a + density.b + 1.0))
        * (1.0 - inbreeding);
    let heterozygous_rate = density.p_segregating() * per_segregating;
    let mean_carrier = CARRIER.0 / (CARRIER.0 + CARRIER.1);
    let duplicated_share = if planted {
        heterozygous_rate * DUPLICATED_PER_HETEROZYGOUS
            / (COLLAPSE_SHARE + (1.0 - COLLAPSE_SHARE) * mean_carrier)
    } else {
        0.0
    };

    for index in 0..positions {
        let rate = if rng.uniform() < NOISY_SHARE {
            NOISY
        } else {
            CLEAN
        };
        let allele = 1 + ((rng.uniform() * 3.0) as usize).min(2);
        let is_duplicated = rng.uniform() < duplicated_share;
        // A duplication is carried by everybody, or by a share drawn from the same Beta the eight
        // tomato samples' counts were fitted to.
        let carrier_frequency = if rng.uniform() < COLLAPSE_SHARE {
            1.0
        } else {
            rng.beta(CARRIER.0, CARRIER.1)
        };
        let branch = rng.pick(&[
            density.p_invariant,
            density.p_fixed_alt,
            density.p_segregating(),
        ]);
        let f = match branch {
            0 => 0.0,
            1 => 1.0,
            _ => rng.beta(density.a, density.b),
        };
        if !is_duplicated && branch == 2 {
            expected += 2.0 * f * (1.0 - f);
        }
        let mut carried_here = 0;
        for sample in 0..samples {
            // **A carrier collects twice the reads**, which is the signal the coverage summary
            // reads and which neither discriminator is allowed to see position by position.
            let copies = if is_duplicated && rng.uniform() < carrier_frequency {
                2.0
            } else {
                1.0
            };
            if copies > 1.0 {
                carried_here += 1;
            }
            let reads = rng.poisson(depth * copies);
            let genotype = if is_duplicated {
                // A carrier reads about half non-reference; a non-carrier is homozygous
                // reference. **There is no third state**, and that is the class's whole
                // signature.
                if copies > 1.0 { 1_usize } else { 0 }
            } else {
                match branch {
                    0 => 0,
                    1 => 2,
                    _ => {
                        let het = 2.0 * f * (1.0 - f) * (1.0 - inbreeding);
                        let shift = inbreeding * f * (1.0 - f);
                        rng.pick(&[(1.0 - f) * (1.0 - f) + shift, het, f * f + shift])
                    }
                }
            };
            if genotype == 1 && !is_duplicated {
                heterozygous[sample] += 1;
            }
            let reading = (copies + coverage_sd(copies, aligned) * rng.normal()).max(0.02);
            let ln_normal = |mean: f64| {
                let sd = coverage_sd(mean, aligned);
                -0.5 * ((reading - mean) / sd).powi(2) - sd.ln()
            };
            odds[sample][index] = (ln_normal(2.0) - ln_normal(1.0)) as f32;

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
        if carried_here > 0 {
            duplicated_positions += 1;
            carrier_pairs += carried_here;
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
        depth_cap: DepthCap::new(124),
    };
    Drawn {
        samples: (0..samples)
            .map(|s| {
                SampleCensusEvidence::resident(
                    format!("s{s:02}"),
                    terms.clone(),
                    BTreeMap::from([(
                        SectionKey::Generic(ReadGroupId(0)),
                        Section::Generic(GenericEvidence::from_parts(
                            std::mem::replace(&mut codes[s], PackedDepthCodes::never_walked(0)),
                            std::mem::take(&mut sparse[s]),
                        )),
                    )]),
                )
            })
            .collect(),
        coverage_odds: odds.into_iter().map(Arc::from).collect(),
        heterozygosity: heterozygous.iter().sum::<u64>() as f64 / (samples * positions) as f64,
        expected_heterozygosity: expected / positions as f64,
        duplicated_positions,
        carrier_pairs,
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
