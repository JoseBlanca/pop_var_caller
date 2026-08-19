//! **Does the duplicated class still eat the heterozygosity when the truth is drawn from
//! outside the model the fit assumes?**
//!
//! The fit describes a population's allele frequencies with four numbers: the share of
//! positions carrying only the reference, the share carrying only a non-reference base, and
//! the two shapes of a Beta for everything in between. Every drawn cohort measured so far drew
//! its frequencies from exactly that Beta, so those runs grade the arithmetic and cannot say
//! whether the description fits a real population — and a class free to claim positions at any
//! frequency is exactly what would absorb the mismatch if it did not.
//!
//! This program draws the frequencies from **drift instead**, which is where a real
//! population's come from and which is not a Beta:
//!
//! * A segregating position is one whose alternative allele is carried by `k` of the panel's
//!   `2n` chromosomes, with `k` drawn from `1/k` — the frequency spectrum a neutral,
//!   constant-sized population leaves in a sample, so most variants are carried by one or two
//!   chromosomes.
//! * The `k` copies are then **dealt out among the panel's chromosomes**, so the samples are
//!   not independent draws at a shared frequency: one sample carrying the allele makes another
//!   slightly less likely to. The fit assumes the opposite.
//! * A duplicated stretch drifts the same way, rather than being handed a carrier frequency
//!   from a Beta fitted to eight tomato accessions.
//!
//! Three cohorts are fitted, each with the class on and off:
//!
//! 1. **Drift, with duplications planted** at the rate the duplication probe measured on real
//!    alignments — 300 carrier positions in two million per accession.
//! 2. **Drift, with no duplications at all.** This is the control that decides H3: if the class
//!    still takes the heterozygosity here, it is absorbing the difference between drift and a
//!    Beta rather than finding duplications.
//! 3. **A Beta, no duplications**, drawn from the fit's own family at the same heterozygosity,
//!    which is what says how much of arm 2 belongs to the spectrum's shape.
//!
//! ```text
//! ng_joint_duplicated_drifting [positions] [panel sizes, comma separated]
//! ```

use std::collections::BTreeMap;
use std::time::Instant;

use pop_var_caller::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
use pop_var_caller::ng::parameter_estimation::joint::census::{
    AlleleObservation, CohortCensusEvidence, DepthCap, DepthCode, DepthLadderDigest,
    GenericEvidence, ObservedAllele, PackedDepthCodes, ReadCap, RecordingTerms,
    SampleCensusEvidence, Section, SectionKey, SelectionTermsDigest,
};
use pop_var_caller::ng::parameter_estimation::joint::fit::{JointFitConfig, fit_jointly};
use pop_var_caller::ng::parameter_estimation::joint::loci::{
    CatalogBuildSettings, CensusLociDigester, ReferenceDigest, RegionSetDigest, SelectionTerms,
};
use pop_var_caller::ng::repeat_catalog::StrRepeatCriteria;
use pop_var_caller::ng::tandem_repeat::ScanParams;
use pop_var_caller::ng::types::ReadGroupId;

/// How often a read misreads a base at an ordinary position, and at a mismapped one, and how
/// many positions mismap. The tomato panel's own fitted values.
const CLEAN: f64 = 0.0034;
const NOISY: f64 = 0.0254;
const NOISY_SHARE: f64 = 0.03;

/// How much less heterozygous every drawn plant is than random mating in the panel predicts —
/// the tomato panel's median with the class off.
const HOMOZYGOTE_EXCESS: f64 = 0.79;

/// The share of positions that segregate at all. Tuned so that a plant comes out near one
/// heterozygous position per kilobase, which is what a rough SNP caller measured on tomato.
const SEGREGATING_SHARE: f64 = 0.030;

/// **How many positions an accession carries an extra copy of, per two million.** The
/// duplication probe measured 150 to 590 on eight tomato accessions
/// (`duplicated_locus_probe_2026-08-12.md` §6); this is the middle of that range.
const CARRIER_POSITIONS_PER_TWO_MILLION: f64 = 300.0;

/// The Beta the third arm draws its frequencies from — **the tomato panel's own fitted shapes
/// with the duplicated class off**, so the comparison is against the shape the fit believes the
/// real cohort has rather than against an arbitrary one.
const BETA_SHAPES: (f64, f64) = (0.555, 6.151);

/// How many more of its positions the Beta arm segregates at. The two spectra put different
/// mass at intermediate frequencies, so the same share of segregating positions would leave the
/// Beta arm at 0.69 heterozygous positions per kilobase against drift's 1.74 — and a comparison
/// between two arms at different heterozygosities is not the one this program is for.
const BETA_SEGREGATING_SCALE: f64 = 2.5;

/// The shallowest and deepest accession in the tomato panel, in reads a position. Every drawn
/// panel spreads its samples geometrically between them, because a class that behaves at one
/// depth and not another would otherwise be graded at whichever depth this program picked.
const DEPTH_RANGE: (f64, f64) = (2.4, 30.6);

fn main() {
    let mut args = std::env::args().skip(1);
    let positions: usize = args.next().map_or(200_000, |a| a.parse().expect("a count"));
    let panels: Vec<usize> = args.next().map_or_else(
        || vec![25, 63],
        |a| {
            a.split(',')
                .map(|part| part.parse().expect("a panel size"))
                .collect()
        },
    );

    println!("positions        {positions}");
    println!(
        "depth            {} to {} reads a position, spread geometrically across the panel",
        DEPTH_RANGE.0, DEPTH_RANGE.1
    );
    println!(
        "the population   {:.1}% of positions segregate, the alternative allele carried by k of \
         the panel's chromosomes with k drawn from 1/k",
        100.0 * SEGREGATING_SHARE
    );
    println!("every plant      {HOMOZYGOTE_EXCESS} less heterozygous than random mating predicts");
    println!(
        "duplications     {CARRIER_POSITIONS_PER_TWO_MILLION:.0} carrier positions per two \
         million per accession, drifting the same way"
    );

    for samples in &panels {
        for arm in [Arm::DriftWithDuplications, Arm::Drift, Arm::Beta] {
            let drawn = draw(*samples, positions, arm, 0x9E37_79B9_7F4A_7C15);
            println!("\n=== {samples} samples, {} ===", arm.describe());
            println!(
                "  drawn: heterozygosity {:.5} ({:.3} per kilobase), homozygote excess \
                 {HOMOZYGOTE_EXCESS}, {} carrier positions ({} of them somebody carries)",
                drawn.heterozygosity,
                1_000.0 * drawn.heterozygosity,
                drawn.carrier_pairs,
                drawn.duplicated_positions
            );
            println!(
                "\n  {:<20}{:>12}{:>12}{:>14}{:>16}{:>10}",
                "", "het/kb", "against", "less het by", "class weight", "time"
            );
            for with_class in [false, true] {
                let at = Instant::now();
                let config = JointFitConfig {
                    quadrature_nodes: 12,
                    max_passes: 200,
                    duplicated_positions: with_class,
                    ..JointFitConfig::default()
                };
                let mut cohort = CohortCensusEvidence::new(drawn.samples.clone())
                    .expect("a drawn cohort records one way");
                let fit = fit_jointly(&mut cohort, &config).expect("a drawn cohort pools");
                let mut heterozygosity: Vec<f64> = drawn
                    .samples
                    .iter()
                    .map(|sample| fit.rates[&sample.sample].value.heterozygous)
                    .collect();
                let mut excess: Vec<f64> = drawn
                    .samples
                    .iter()
                    .map(|sample| fit.hom_excess[&sample.sample].value.get())
                    .collect();
                heterozygosity.sort_by(f64::total_cmp);
                excess.sort_by(f64::total_cmp);
                let median = heterozygosity[heterozygosity.len() / 2];
                println!(
                    "  {:<20}{:>12.3}{:>11.1}%{:>14.3}{:>16.5}{:>8.0} s   {}",
                    if with_class {
                        "the class on"
                    } else {
                        "the class off"
                    },
                    1_000.0 * median,
                    100.0 * (median / drawn.heterozygosity - 1.0),
                    excess[excess.len() / 2],
                    fit.duplicated.as_ref().map_or(0.0, |d| d.value.share),
                    at.elapsed().as_secs_f64(),
                    if fit.converged { "" } else { "RAN OUT" },
                );
            }
        }
    }
}

#[derive(Copy, Clone, PartialEq)]
enum Arm {
    DriftWithDuplications,
    Drift,
    Beta,
}

impl Arm {
    fn describe(self) -> &'static str {
        match self {
            Arm::DriftWithDuplications => {
                "frequencies from drift, duplications planted at the probe's rate"
            }
            Arm::Drift => "frequencies from drift, no duplications at all",
            Arm::Beta => "frequencies from the fit's own Beta, no duplications at all",
        }
    }
}

struct Drawn {
    samples: Vec<SampleCensusEvidence>,
    /// The share of positions at which a sample was drawn genuinely heterozygous, averaged over
    /// the panel — **the truth the fitted heterozygosity is compared against**.
    heterozygosity: f64,
    duplicated_positions: usize,
    /// (position, sample) pairs that really carry an extra copy.
    carrier_pairs: usize,
}

fn draw(samples: usize, positions: usize, arm: Arm, seed: u64) -> Drawn {
    let mut rng = Rng(seed);
    let edges = DepthBinEdges::new();
    let mut codes: Vec<PackedDepthCodes> = (0..samples)
        .map(|_| PackedDepthCodes::never_walked(positions))
        .collect();
    let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
    let mut heterozygous = vec![0_u64; samples];
    let (mut duplicated_positions, mut carrier_pairs) = (0_usize, 0_usize);

    // Each accession's own depth, spread geometrically over the tomato panel's range.
    let depths: Vec<f64> = (0..samples)
        .map(|s| {
            let t = if samples > 1 {
                s as f64 / (samples - 1) as f64
            } else {
                0.0
            };
            DEPTH_RANGE.0 * (DEPTH_RANGE.1 / DEPTH_RANGE.0).powf(t)
        })
        .collect();

    // **The spectrum, and it is where this program differs from every other drawn cohort.**
    // `k` alternative chromosomes out of `2n`, with `k` drawn from `1/k` — so a position
    // carried by one chromosome is `2n − 2` times as likely as one carried by half the panel.
    let chromosomes = 2 * samples;
    let spectrum: Vec<f64> = (1..chromosomes).map(|k| 1.0 / k as f64).collect();

    // How often a position is duplicated, so that an accession ends up carrying the number the
    // duplication probe measured. A sample carries where either of its chromosomes does.
    let mean_carrier: f64 = {
        let total: f64 = spectrum.iter().sum();
        spectrum
            .iter()
            .enumerate()
            .map(|(index, weight)| {
                let f = (index + 1) as f64 / chromosomes as f64;
                weight / total * (1.0 - (1.0 - f) * (1.0 - f))
            })
            .sum()
    };
    let duplicated_share = if arm == Arm::DriftWithDuplications {
        CARRIER_POSITIONS_PER_TWO_MILLION / 2e6 / mean_carrier
    } else {
        0.0
    };

    let mut chromosome_carries = vec![false; chromosomes];
    for index in 0..positions {
        let rate = if rng.uniform() < NOISY_SHARE {
            NOISY
        } else {
            CLEAN
        };
        let allele = 1 + ((rng.uniform() * 3.0) as usize).min(2);
        let duplicated = rng.uniform() < duplicated_share;
        let segregating_share = SEGREGATING_SHARE
            * if arm == Arm::Beta {
                BETA_SEGREGATING_SCALE
            } else {
                1.0
            };
        let segregating = !duplicated && rng.uniform() < segregating_share;

        // Which chromosomes carry the alternative allele — or the extra copy. Both are dealt
        // out among the panel rather than drawn independently per sample, which is what makes
        // a fixed allele count and what the fit's own model does not have.
        chromosome_carries.fill(false);
        if duplicated {
            let count = 1 + rng.pick(&spectrum);
            rng.deal(&mut chromosome_carries, count);
        } else if segregating {
            match arm {
                // The fit's own family: one frequency for the population, and every chromosome
                // an independent toss at it.
                Arm::Beta => {
                    let f = rng.beta(BETA_SHAPES.0, BETA_SHAPES.1);
                    for slot in chromosome_carries.iter_mut() {
                        *slot = rng.uniform() < f;
                    }
                }
                _ => {
                    let count = 1 + rng.pick(&spectrum);
                    rng.deal(&mut chromosome_carries, count);
                }
            }
        }
        let mut carried_here = 0;
        for sample in 0..samples {
            // **Inbreeding is drawn as autozygosity**: with probability `HOMOZYGOTE_EXCESS` a
            // plant's two chromosomes are copies of one, which is what makes it less
            // heterozygous than the panel's frequencies predict.
            let first = chromosome_carries[sample * 2];
            let second = if rng.uniform() < HOMOZYGOTE_EXCESS {
                first
            } else {
                chromosome_carries[sample * 2 + 1]
            };
            let copies = if duplicated && (first || second) {
                carried_here += 1;
                2.0
            } else {
                1.0
            };
            let genotype = if duplicated {
                // A carrier reads about half non-reference and a non-carrier is homozygous
                // reference. **There is no third state**, and that absence is the class's whole
                // evidence across a cohort.
                usize::from(copies > 1.0)
            } else {
                usize::from(first) + usize::from(second)
            };
            if genotype == 1 && !duplicated {
                heterozygous[sample] += 1;
            }
            let reads = rng.poisson(depths[sample] * copies);

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
            codes[sample].set(index, DepthCode::Binned(edges.bin_for(reads.min(124))));
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
        heterozygosity: heterozygous.iter().sum::<u64>() as f64 / (samples * positions) as f64,
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

    /// Set exactly `count` of the slots, chosen uniformly among them — the deal that makes an
    /// allele count fixed rather than a per-sample coin toss.
    fn deal(&mut self, slots: &mut [bool], count: usize) {
        let mut remaining = count.min(slots.len());
        let mut left = slots.len();
        for slot in slots.iter_mut() {
            if remaining == 0 {
                break;
            }
            if self.uniform() < remaining as f64 / left as f64 {
                *slot = true;
                remaining -= 1;
            }
            left -= 1;
        }
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
