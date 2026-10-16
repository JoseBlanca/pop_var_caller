//! **A whole motif length's repeat classes, drawn from a known model at the level of read counts**
//! — the count-level simulator of the fitting research plan
//! (`doc/devel/ng/research/str_fit_research_plan.md` §3.1).
//!
//! The census hands the repeat-tract fit, per class (one motif length and one reference repeat
//! count), every sample's reads at every locus, counted by how many repeat units each read's tract
//! differs from the reference. This module draws exactly those counts from a model whose every
//! number is known, without simulating reads or aligning anything, so a 2,000-sample cohort is drawn
//! in seconds and a fit's answer can be held against the truth.
//!
//! # The model, as the draw runs
//!
//! For each class, every number comes from a **curve across the reference repeat count**
//! ([`ParameterCurve`]), so a whole motif length is described by a handful of curves:
//!
//! 1. **The class's allele-length shares**: how the class's chromosomes are spread over allele
//!    lengths, from `-reach` to `+reach` repeat units around the reference
//!    ([`CohortModel::true_allele_reach`], 14 by default). Either a six-number shape
//!    ([`AlleleLengthShape`]) or explicit share vectors.
//! 2. **Each locus's own allele frequencies**, drawn from a Dirichlet distribution centred on those
//!    shares with concentration `c`. The model is described by the between-locus share
//!    `F = 1 / (1 + c)`, the Balding–Nichols F_ST with loci in place of populations: at 0 every
//!    locus carries the class's shares exactly, at 1 every locus is fixed at one length.
//! 3. **Each sample's genotype at the locus**: with probability equal to the sample's homozygote
//!    excess its two copies are one draw from the locus's frequencies, otherwise two independent
//!    draws.
//! 4. **Each sample's reads at the locus**: a Poisson number with mean [`CohortModel::depth`], each
//!    read from either copy with even odds, its length reported through the stutter of the sample's
//!    library ([`Slippage::read_probabilities`] at [`library_level`]). A sample with no read at a
//!    locus is left out of that locus, as the census leaves it out.
//!
//! Reads are recorded to `±read_span` (the census's [`RECORDED_OFFSET_RANGE`], 8) with the end
//! buckets saturating, while alleles reach ±14 — so the fit's outermost allele classes (±9, which
//! mean "this far or further") are exercised by alleles past them.
//!
//! # What comes back
//!
//! [`DrawnPeriod`]: one [`StratumEvidence`] per class, exactly as the census would build it, the
//! per-sample homozygote excess, the truth per class ([`ClassTruth`]), which
//! [`write_truth_tsv`] writes out, and the genotype every sample drew at every locus.
//!
//! **The draw is deterministic.** Each class is drawn from its own random stream, seeded from the
//! model's seed and the class's repeat count, so a class's draw does not change when other classes
//! are added or removed, and the classes are drawn in parallel.
//!
//! # The default truth
//!
//! [`CohortModel::tomato_like`] gives numbers that look like the tomato cohorts; each curve's
//! comment says where its numbers come from.
//!
//! # Population structure
//!
//! One population only. The plan defers structure (§6, decision 3); [`PopulationStructure`] is the
//! place it will go.

use std::io::Write;

use serde::Deserialize;

use super::*;

/// **How far either side of the reference the simulator places alleles by default: 14 repeat
/// units**, five past the fit's outermost allele class (±9), so that class's "this far or further"
/// meaning is tested by alleles that really are further.
pub const DEFAULT_TRUE_ALLELE_REACH: i32 = 14;

// ---------------------------------------------------------------------
// The curves
// ---------------------------------------------------------------------

/// **One number of the model as a function of the reference repeat count.**
///
/// Written in a model file as `{ constant = 0.7 }` or `{ knots = [[8, 0.005], [16, 0.27]] }`.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum ParameterCurve {
    /// The same value at every repeat count.
    Constant(f64),
    /// `(repeat count, value)` points, joined by straight lines and held flat past the first and
    /// last point. The points must be in increasing repeat count.
    Knots(Vec<(f64, f64)>),
}

impl ParameterCurve {
    /// The curve's value at `repeats`.
    ///
    /// # Panics
    ///
    /// When a curve of knots has none.
    #[must_use]
    pub fn at(&self, repeats: u64) -> f64 {
        match self {
            Self::Constant(value) => *value,
            Self::Knots(knots) => {
                assert!(
                    !knots.is_empty(),
                    "a curve of knots needs at least one knot"
                );
                let x = repeats as f64;
                let (first, last) = (knots[0], knots[knots.len() - 1]);
                if x <= first.0 {
                    return first.1;
                }
                if x >= last.0 {
                    return last.1;
                }
                let right = knots
                    .iter()
                    .position(|(at, _)| *at >= x)
                    .unwrap_or(knots.len() - 1);
                let (x0, y0) = knots[right - 1];
                let (x1, y1) = knots[right];
                y0 + (y1 - y0) * (x - x0) / (x1 - x0)
            }
        }
    }

    /// The same curve with every value multiplied by `factor`.
    #[must_use]
    pub fn scaled(&self, factor: f64) -> Self {
        match self {
            Self::Constant(value) => Self::Constant(value * factor),
            Self::Knots(knots) => Self::Knots(
                knots
                    .iter()
                    .map(|(at, value)| (*at, value * factor))
                    .collect(),
            ),
        }
    }
}

/// The three stutter numbers of a library at multiplier one, each a curve.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StutterTruth {
    /// How often a read reports a length other than its allele's.
    pub level: ParameterCurve,
    /// Of the reads that slip, the share that report a shorter tract.
    pub shorter_share: ParameterCurve,
    /// How fast two-unit slips fall off against one-unit slips.
    pub fall_off: ParameterCurve,
}

/// **A six-number allele-length shape**, the one the kimura shape comparison tested
/// (`scripts/research/str_allele_length_shapes.py`, `ref_step_power6`), with each number a curve.
///
/// In words: a share at the reference length; the rest split between shorter and longer alleles;
/// on each side a share one repeat away; and the remainder of that side spread over two or more
/// repeats away, falling as a power of the distance.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AlleleLengthShape {
    /// The share of chromosomes at the reference length.
    pub reference_share: ParameterCurve,
    /// Of the chromosomes not at the reference length, the share that are shorter.
    pub shorter_split: ParameterCurve,
    /// Of the shorter chromosomes, the share exactly one repeat shorter.
    pub shorter_one_step_share: ParameterCurve,
    /// How fast the shorter tail falls beyond one repeat: the weight at `d` repeats is `d^-power`.
    pub shorter_power: ParameterCurve,
    /// Of the longer chromosomes, the share exactly one repeat longer.
    pub longer_one_step_share: ParameterCurve,
    /// How fast the longer tail falls beyond one repeat.
    pub longer_power: ParameterCurve,
}

impl AlleleLengthShape {
    /// The shares at `repeats`, over `-reach ..= +reach`.
    #[must_use]
    pub fn shares_at(&self, repeats: u64, reach: i32) -> Vec<f64> {
        let reference = self.reference_share.at(repeats).clamp(0.0, 1.0);
        let shorter = self.shorter_split.at(repeats).clamp(0.0, 1.0);
        let side = |one_step: f64, power: f64| -> Vec<f64> {
            // Distances 1 ..= reach: the first gets `one_step`, the rest share the remainder as
            // d^-power.
            let one_step = one_step.clamp(0.0, 1.0);
            let tail: Vec<f64> = (2..=reach)
                .map(|distance| float::powf(f64::from(distance), -power))
                .collect();
            let tail_total: f64 = tail.iter().sum();
            let mut out = vec![if reach >= 2 { one_step } else { 1.0 }];
            out.extend(
                tail.iter()
                    .map(|weight| (1.0 - one_step) * weight / tail_total),
            );
            out
        };
        let below = side(
            self.shorter_one_step_share.at(repeats),
            self.shorter_power.at(repeats),
        );
        let above = side(
            self.longer_one_step_share.at(repeats),
            self.longer_power.at(repeats),
        );
        let mut shares = vec![0.0; (2 * reach + 1) as usize];
        shares[reach as usize] = reference;
        for distance in 1..=reach {
            let index = (distance - 1) as usize;
            shares[(reach - distance) as usize] = (1.0 - reference) * shorter * below[index];
            shares[(reach + distance) as usize] =
                (1.0 - reference) * (1.0 - shorter) * above[index];
        }
        normalise(&mut shares);
        shares
    }
}

/// One class's explicit allele-length shares.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExplicitShares {
    pub repeats: u64,
    /// Over `-reach ..= +reach`; normalised when drawn.
    pub shares: Vec<f64>,
}

/// How the classes' allele-length shares are given.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum AlleleLengthVariationTruth {
    /// A six-number shape whose numbers are curves across repeat count.
    Shape(AlleleLengthShape),
    /// A share vector per class; a class not listed takes the nearest listed one.
    Explicit(Vec<ExplicitShares>),
}

impl AlleleLengthVariationTruth {
    fn shares_at(&self, repeats: u64, reach: i32) -> Vec<f64> {
        match self {
            Self::Shape(shape) => shape.shares_at(repeats, reach),
            Self::Explicit(listed) => {
                let nearest = listed
                    .iter()
                    .min_by_key(|class| class.repeats.abs_diff(repeats))
                    .expect("explicit shares list at least one class");
                assert_eq!(
                    nearest.shares.len(),
                    (2 * reach + 1) as usize,
                    "explicit shares run from -reach to +reach, {} values",
                    2 * reach + 1
                );
                let mut shares = nearest.shares.clone();
                normalise(&mut shares);
                shares
            }
        }
    }
}

// ---------------------------------------------------------------------
// The cohort
// ---------------------------------------------------------------------

/// How many loci each class holds.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum LociPerClass {
    /// The same number in every class.
    Each(usize),
    /// As many as tomato2 has in that class ([`TOMATO2_LOCI`]), optionally capped.
    Tomato2 { at_most: Option<usize> },
    /// `(repeat count, loci)` pairs; a class not listed has none.
    Listed(Vec<(u64, usize)>),
}

/// Each sample's homozygote excess: the chance its two copies at a locus are one draw.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(rename_all = "snake_case", deny_unknown_fields)]
pub enum HomozygoteExcessTruth {
    /// Every sample the same.
    Each(f64),
    /// Each sample drawn uniformly between the two.
    Uniform { low: f64, high: f64 },
}

/// **How the samples are related.** One population today.
///
/// The research plan defers population structure until the one-population model is settled (§6,
/// decision 3). TODO(plan §6.3): a `TwoSubpopulations` variant, each locus's frequencies drawn
/// around the class's shares and then each subpopulation's around the locus's, with a second
/// between-subpopulation share.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum PopulationStructure {
    #[default]
    OnePopulation,
}

/// Which default stutter to draw from: a PCR library, or a PCR-free one.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StutterRegime {
    Pcr,
    /// Slippage levels ten times lower than [`Pcr`](Self::Pcr)'s, everything else the same.
    PcrFree,
}

/// **Everything the simulator draws a motif length from.** Built with [`CohortModel::tomato_like`]
/// and adjusted field by field, or read from a TOML file whose keys are these fields.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CohortModel {
    /// The motif length, in bases.
    pub period: u8,
    /// The reference repeat counts to draw. `None` takes every class [`LociPerClass`] gives loci
    /// to (for [`LociPerClass::Each`], every tomato2 class of this period).
    #[serde(default)]
    pub repeat_counts: Option<Vec<u64>>,
    pub loci: LociPerClass,
    pub stutter: StutterTruth,
    /// The between-locus share `F = 1 / (1 + c)` of each class, `c` the Dirichlet concentration.
    pub between_locus_share: ParameterCurve,
    pub allele_length_variation: AlleleLengthVariationTruth,
    /// How far either side of the reference alleles are placed ([`DEFAULT_TRUE_ALLELE_REACH`]).
    #[serde(default = "default_reach")]
    pub true_allele_reach: i32,
    /// How many samples.
    pub samples: usize,
    /// The mean number of reads spanning a locus, per sample per locus (Poisson). Spanning reads,
    /// not coverage: at a long tract fewer reads span it than cover a position.
    pub depth: f64,
    pub homozygote_excess: HomozygoteExcessTruth,
    /// One stutter multiplier per library (`[1.0]` by default); sample `s` is in library
    /// `s mod libraries`.
    #[serde(default = "default_multipliers")]
    pub library_multipliers: Vec<f64>,
    /// Whether the evidence's library list carries the true multipliers (as if a previous fit had
    /// found them) or one for every library (as if none had been fitted). The truth's class levels
    /// are always at multiplier one.
    #[serde(default)]
    pub multipliers_told: bool,
    #[serde(default)]
    pub structure: PopulationStructure,
    pub seed: u64,
}

fn default_reach() -> i32 {
    DEFAULT_TRUE_ALLELE_REACH
}

fn default_multipliers() -> Vec<f64> {
    vec![1.0]
}

/// **Loci per class in tomato2** (59 accessions, about three reads a position): period, reference
/// repeat count, and the loci with at least one spanning read. Copied from the counting fit's run
/// at ±9 on 2026-10-11 (`pop_var_caller-stutter-multiplier/tmp/window/tomato2_counted_span9_fixed.tsv`,
/// column `tracts`), which lists the classes with at least 8 such loci.
pub const TOMATO2_LOCI: &[(u8, u64, usize)] = &[
    (1, 8, 4947),
    (1, 9, 4045),
    (1, 10, 2014),
    (1, 11, 821),
    (1, 12, 428),
    (1, 13, 267),
    (1, 14, 154),
    (1, 15, 88),
    (1, 16, 60),
    (1, 17, 32),
    (1, 18, 22),
    (1, 19, 30),
    (1, 20, 23),
    (1, 21, 16),
    (2, 6, 426),
    (2, 7, 181),
    (2, 8, 142),
    (2, 9, 101),
    (2, 10, 71),
    (2, 11, 60),
    (2, 12, 44),
    (2, 13, 33),
    (2, 14, 36),
    (2, 15, 20),
    (2, 16, 34),
    (2, 17, 25),
    (2, 18, 25),
    (2, 19, 19),
    (2, 20, 21),
    (2, 21, 13),
    (2, 22, 17),
    (2, 23, 16),
    (2, 24, 11),
    (2, 25, 21),
    (2, 27, 12),
    (3, 6, 120),
    (3, 7, 64),
    (3, 8, 31),
    (3, 9, 16),
    (3, 10, 15),
    (3, 11, 14),
    (3, 12, 10),
    (3, 13, 9),
    (4, 6, 15),
    (6, 4, 20),
];

impl CohortModel {
    /// **A truth that looks like the tomato cohorts**, for period 1 or 2, at 60 samples and three
    /// reads, loci as tomato2 has them. Where each number comes from:
    ///
    /// - **Slippage level**: tomato2's counting fit at ±9 (`TOMATO2_LOCI`'s source), class by class,
    ///   lightly smoothed: one-base runs from 0.005 at 8 repeats to 0.27 at 16, two-base runs from
    ///   0.0075 at 6 to 0.17 at 14 and 0.21 at 18. Held flat past the last knot.
    ///   [`StutterRegime::PcrFree`] divides every level by ten.
    /// - **Shorter share**: the same run gives 0.53–0.62 for one-base runs and 0.74–0.90 for
    ///   two-base runs; the curves run 0.55 → 0.60 and 0.85 → 0.78.
    /// - **Fall-off**: the same run gives 0.05–0.21 for one-base runs up to 16 repeats and
    ///   0.02–0.22 for two-base runs up to 13; the curves rise from 0.10 to 0.30 at 18 repeats.
    /// - **Between-locus share F**: 0.7 throughout. tomato2's one-base runs give concentrations of
    ///   0.45–0.81 (F 0.55–0.69); its two-base runs give 1.0–2.5 (F 0.29–0.50), lower than this.
    /// - **Allele-length shape**: the kimura refits of the six-number shape at ±6
    ///   (`pop_var_caller/tmp/str_shapes_kimura_2026-10-11.md`, question 2): reference share 0.93
    ///   at 8 repeats falling to 0.50 at 16 (one-base); the one-repeat shares falling from
    ///   0.96/0.83 to 0.33/0.60; tail powers held near 1.5–3, where the refits above 18 repeats
    ///   gave 0, an artefact of the ±6 window. The two-base curves start at kimura's only two-base
    ///   class (6 repeats: reference share 0.96, shorter split 0.26) and fall to a reference share
    ///   of 0.5 by 16 repeats.
    /// - **Homozygote excess**: uniform between 0.85 and 0.99 a sample; tomato accessions are
    ///   largely inbred lines.
    ///
    /// # Panics
    ///
    /// For a period other than 1 or 2: no default truth is written for them; give a model file.
    #[must_use]
    pub fn tomato_like(period: u8, regime: StutterRegime) -> Self {
        let knots = |points: &[(f64, f64)]| ParameterCurve::Knots(points.to_vec());
        let (stutter, shape) = match period {
            1 => (
                StutterTruth {
                    level: knots(&[
                        (8.0, 0.005),
                        (9.0, 0.016),
                        (10.0, 0.040),
                        (11.0, 0.078),
                        (12.0, 0.123),
                        (13.0, 0.160),
                        (14.0, 0.190),
                        (15.0, 0.230),
                        (16.0, 0.270),
                    ]),
                    shorter_share: knots(&[(8.0, 0.55), (16.0, 0.60)]),
                    fall_off: knots(&[(8.0, 0.10), (18.0, 0.30)]),
                },
                AlleleLengthShape {
                    reference_share: knots(&[
                        (8.0, 0.93),
                        (10.0, 0.77),
                        (12.0, 0.60),
                        (14.0, 0.53),
                        (16.0, 0.50),
                    ]),
                    shorter_split: knots(&[(8.0, 0.35), (9.0, 0.50), (16.0, 0.60)]),
                    shorter_one_step_share: knots(&[
                        (8.0, 0.96),
                        (12.0, 0.60),
                        (16.0, 0.42),
                        (23.0, 0.33),
                    ]),
                    shorter_power: knots(&[(8.0, 2.0), (11.0, 3.0), (16.0, 1.5)]),
                    longer_one_step_share: knots(&[(8.0, 0.83), (10.0, 0.62), (16.0, 0.60)]),
                    longer_power: ParameterCurve::Constant(2.0),
                },
            ),
            2 => (
                StutterTruth {
                    level: knots(&[
                        (6.0, 0.0075),
                        (7.0, 0.023),
                        (8.0, 0.044),
                        (9.0, 0.076),
                        (10.0, 0.110),
                        (11.0, 0.130),
                        (12.0, 0.145),
                        (13.0, 0.157),
                        (14.0, 0.170),
                        (18.0, 0.210),
                    ]),
                    shorter_share: knots(&[(6.0, 0.85), (14.0, 0.78)]),
                    fall_off: knots(&[(6.0, 0.10), (18.0, 0.30)]),
                },
                AlleleLengthShape {
                    reference_share: knots(&[
                        (6.0, 0.96),
                        (8.0, 0.85),
                        (10.0, 0.72),
                        (12.0, 0.60),
                        (14.0, 0.53),
                        (16.0, 0.50),
                    ]),
                    shorter_split: knots(&[(6.0, 0.26), (10.0, 0.50)]),
                    shorter_one_step_share: knots(&[(6.0, 0.90), (14.0, 0.50)]),
                    shorter_power: knots(&[(6.0, 2.6), (14.0, 1.8)]),
                    longer_one_step_share: knots(&[(6.0, 0.93), (14.0, 0.55)]),
                    longer_power: knots(&[(6.0, 2.7), (14.0, 1.8)]),
                },
            ),
            other => panic!(
                "no default truth is written for period {other}; give the simulator a model file"
            ),
        };
        let stutter = match regime {
            StutterRegime::Pcr => stutter,
            StutterRegime::PcrFree => StutterTruth {
                level: stutter.level.scaled(0.1),
                ..stutter
            },
        };
        Self {
            period,
            repeat_counts: None,
            loci: LociPerClass::Tomato2 { at_most: None },
            stutter,
            between_locus_share: ParameterCurve::Constant(0.7),
            allele_length_variation: AlleleLengthVariationTruth::Shape(shape),
            true_allele_reach: DEFAULT_TRUE_ALLELE_REACH,
            samples: 60,
            depth: 3.0,
            homozygote_excess: HomozygoteExcessTruth::Uniform {
                low: 0.85,
                high: 0.99,
            },
            library_multipliers: default_multipliers(),
            multipliers_told: false,
            structure: PopulationStructure::OnePopulation,
            seed: 1,
        }
    }

    /// The classes to draw and how many loci each holds, in increasing repeat count; classes with
    /// no loci are left out.
    #[must_use]
    pub fn classes(&self) -> Vec<(u64, usize)> {
        let tomato2: Vec<(u64, usize)> = TOMATO2_LOCI
            .iter()
            .filter(|(period, _, _)| *period == self.period)
            .map(|(_, repeats, loci)| (*repeats, *loci))
            .collect();
        let mut classes: Vec<(u64, usize)> = match &self.loci {
            LociPerClass::Each(loci) => {
                let repeats: Vec<u64> = self
                    .repeat_counts
                    .clone()
                    .unwrap_or_else(|| tomato2.iter().map(|(repeats, _)| *repeats).collect());
                repeats
                    .into_iter()
                    .map(|repeats| (repeats, *loci))
                    .collect()
            }
            LociPerClass::Tomato2 { at_most } => tomato2
                .iter()
                .filter(|(repeats, _)| {
                    self.repeat_counts
                        .as_ref()
                        .is_none_or(|wanted| wanted.contains(repeats))
                })
                .map(|(repeats, loci)| (*repeats, at_most.map_or(*loci, |cap| (*loci).min(cap))))
                .collect(),
            LociPerClass::Listed(listed) => listed
                .iter()
                .filter(|(repeats, _)| {
                    self.repeat_counts
                        .as_ref()
                        .is_none_or(|wanted| wanted.contains(repeats))
                })
                .copied()
                .collect(),
        };
        classes.retain(|(_, loci)| *loci > 0);
        classes.sort_unstable();
        classes.dedup_by_key(|(repeats, _)| *repeats);
        classes
    }

    /// The truth of the class at `repeats` holding `loci` loci.
    #[must_use]
    pub fn truth_at(&self, repeats: u64, loci: usize) -> ClassTruth {
        let between_locus_share = self.between_locus_share.at(repeats).clamp(1e-9, 1.0 - 1e-9);
        ClassTruth {
            stratum: Stratum {
                period: self.period,
                reference_repeats: repeats,
            },
            loci,
            slippage: Slippage {
                level: self.stutter.level.at(repeats),
                shorter_share: self.stutter.shorter_share.at(repeats),
                fall_off: self.stutter.fall_off.at(repeats),
            },
            concentration: 1.0 / between_locus_share - 1.0,
            true_allele_reach: self.true_allele_reach,
            true_shares: self
                .allele_length_variation
                .shares_at(repeats, self.true_allele_reach),
        }
    }
}

/// **The truth of one drawn class**, at multiplier one.
#[derive(Debug, Clone, PartialEq)]
pub struct ClassTruth {
    pub stratum: Stratum,
    /// Loci drawn; those where no sample drew a read are in the evidence with no samples.
    pub loci: usize,
    pub slippage: Slippage,
    /// The Dirichlet concentration `c`.
    pub concentration: f64,
    /// How far either side of the reference the shares reach.
    pub true_allele_reach: i32,
    /// The class's allele-length shares, `-true_allele_reach ..= +true_allele_reach`.
    pub true_shares: Vec<f64>,
}

impl ClassTruth {
    /// `F = 1 / (1 + c)`.
    #[must_use]
    pub fn between_locus_share(&self) -> f64 {
        1.0 / (1.0 + self.concentration)
    }

    /// **The shares as a fit at `allele_span` sees them**: `-allele_span ..= +allele_span`, the two
    /// outermost classes holding every allele that far or further.
    ///
    /// The concentration needs no such folding: merging classes of a Dirichlet leaves its
    /// concentration unchanged.
    #[must_use]
    pub fn shares_in_window(&self, allele_span: i32) -> Vec<f64> {
        let reach = self.true_allele_reach;
        let mut out = vec![0.0; (2 * allele_span + 1) as usize];
        for (index, share) in self.true_shares.iter().enumerate() {
            let offset = (index as i32 - reach).clamp(-allele_span, allele_span);
            out[(offset + allele_span) as usize] += share;
        }
        out
    }
}

/// **What a draw returns**: the evidence the census would hand the fit, and the truth behind it.
#[derive(Debug, Clone)]
pub struct DrawnPeriod {
    /// One class each, in increasing repeat count.
    pub strata: Vec<StratumEvidence>,
    /// The truth of each class, in the same order.
    pub truth: Vec<ClassTruth>,
    /// Each sample's homozygote excess, as drawn; the fit takes it as given.
    pub homozygote_excess: Vec<f64>,
    /// **The genotype every sample drew**: per class, per locus, per sample listed at the locus
    /// (in the order of the locus's samples in `strata`), its two copies as repeat offsets from the
    /// reference, `-true_allele_reach ..= +true_allele_reach`, in the order drawn. What the fit
    /// never sees; the bench's proxy for calling scores calls against it. An offset is one byte,
    /// so [`draw_period`] panics if it draws an allele more than 127 repeats out (the default
    /// reach is 14).
    pub genotypes: Vec<Vec<Vec<(i8, i8)>>>,
}

// ---------------------------------------------------------------------
// Drawing
// ---------------------------------------------------------------------

/// How many reads a sample draws at a locus.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) enum Depth {
    /// Exactly this many, every sample at every locus, and a sample with none is still listed —
    /// [`draw_stratum`]'s rule.
    Fixed(u32),
    /// A Poisson number with this mean, and a sample with none is left out — the census's rule.
    Poisson(f64),
}

/// What drawing one locus needs besides its frequencies.
pub(super) struct LocusDrawing<'a> {
    /// Per library, per allele class, the probability of each read bucket.
    pub per_library: &'a [Vec<Vec<f64>>],
    /// How many libraries the samples cycle through.
    pub libraries: usize,
    pub homozygote_excess: &'a [f64],
    pub depth: Depth,
    pub buckets: usize,
}

/// **One locus's reads, every sample**, from the locus's allele `frequencies`.
///
/// The random numbers are taken in the order [`draw_stratum`] has always taken them — for each
/// sample, its first copy, the coin for one draw or two, its second copy, then for each read the
/// copy and the bucket — so that drawing at a fixed depth reproduces it exactly. A Poisson depth
/// takes its number first.
pub(super) fn draw_locus(
    rng: &mut Rng,
    frequencies: &[f64],
    drawing: &LocusDrawing<'_>,
) -> TractReads {
    draw_locus_and_genotypes(rng, frequencies, drawing).0
}

/// [`draw_locus`], and the genotype each sample listed at the locus drew, in the order of the
/// locus's samples: its two copies as indices into `frequencies`. The random numbers are the same.
pub(super) fn draw_locus_and_genotypes(
    rng: &mut Rng,
    frequencies: &[f64],
    drawing: &LocusDrawing<'_>,
) -> (TractReads, Vec<(usize, usize)>) {
    let mut reads = TractReads::default();
    let mut genotypes = Vec::new();
    for sample in 0..drawing.homozygote_excess.len() {
        let (depth, listed_without_reads) = match drawing.depth {
            Depth::Fixed(depth) => (depth, true),
            Depth::Poisson(mean) => (rng.poisson(mean), false),
        };
        if depth == 0 && !listed_without_reads {
            continue;
        }
        let first = rng.categorical(frequencies);
        let second = if rng.uniform() < drawing.homozygote_excess[sample] {
            first
        } else {
            rng.categorical(frequencies)
        };
        let library = sample % drawing.libraries;
        let per_allele = &drawing.per_library[library];
        let mut counts = vec![0_u32; drawing.buckets];
        for _ in 0..depth {
            let allele = if rng.uniform() < 0.5 { first } else { second };
            counts[rng.categorical(&per_allele[allele])] += 1;
        }
        reads.samples.push(SampleTractReads {
            sample: sample as u32,
            by_library: vec![(LibraryIndex(library as u32), counts)],
        });
        genotypes.push((first, second));
    }
    (reads, genotypes)
}

/// A locus's allele frequencies: one Dirichlet draw with the given weights (`concentration ×
/// shares`).
pub(super) fn draw_locus_frequencies(rng: &mut Rng, weights: &[f64]) -> Vec<f64> {
    rng.dirichlet(weights)
}

/// A seed for one part of a draw, mixed from the model's seed and a label, so that parts of a draw
/// do not share a stream.
fn seed_for(seed: u64, label: u64) -> u64 {
    let mut mixer = Rng(seed ^ label.wrapping_mul(0xD1B5_4A32_D192_ED03));
    mixer.next_u64()
}

/// **Draw every class of `model`**, reads recorded to `±read_span` with the ends saturating.
///
/// # Panics
///
/// When the model has no samples, a negative depth, no libraries, or a non-positive multiplier.
#[must_use]
pub fn draw_period(model: &CohortModel, read_span: i32) -> DrawnPeriod {
    assert!(model.samples > 0, "a cohort has at least one sample");
    assert!(
        model.depth >= 0.0,
        "depth is a mean number of reads, not negative"
    );
    assert!(
        !model.library_multipliers.is_empty()
            && model
                .library_multipliers
                .iter()
                .all(|m| *m > 0.0 && m.is_finite()),
        "every library has a finite positive stutter multiplier"
    );
    let PopulationStructure::OnePopulation = model.structure;

    let mut excess_rng = Rng(seed_for(model.seed, u64::MAX));
    let homozygote_excess: Vec<f64> = (0..model.samples)
        .map(|_| match model.homozygote_excess {
            HomozygoteExcessTruth::Each(value) => value,
            HomozygoteExcessTruth::Uniform { low, high } => {
                low + (high - low) * excess_rng.uniform()
            }
        })
        .collect();
    let libraries: Arc<[Library]> = model
        .library_multipliers
        .iter()
        .enumerate()
        .map(|(index, multiplier)| Library {
            stutter_multiplier: if model.multipliers_told {
                *multiplier
            } else {
                1.0
            },
            ..Library::in_group(ReadGroupId(index as u32), 0)
        })
        .collect();

    let classes = model.classes();
    type DrawnClass = (StratumEvidence, ClassTruth, Vec<Vec<(i8, i8)>>);
    let drawn: Vec<DrawnClass> = classes
        .par_iter()
        .map(|(repeats, loci)| {
            let truth = model.truth_at(*repeats, *loci);
            let reach = model.true_allele_reach;
            let per_library: Vec<Vec<Vec<f64>>> = model
                .library_multipliers
                .iter()
                .map(|multiplier| {
                    let slippage = Slippage {
                        level: library_level(truth.slippage.level, *multiplier),
                        ..truth.slippage
                    };
                    (-reach..=reach)
                        .map(|allele| slippage.read_probabilities(allele, read_span))
                        .collect()
                })
                .collect();
            let drawing = LocusDrawing {
                per_library: &per_library,
                libraries: model.library_multipliers.len(),
                homozygote_excess: &homozygote_excess,
                depth: Depth::Poisson(model.depth),
                buckets: (2 * read_span + 1) as usize,
            };
            let weights: Vec<f64> = truth
                .true_shares
                .iter()
                .map(|share| truth.concentration * share)
                .collect();
            let mut rng = Rng(seed_for(
                model.seed,
                (u64::from(model.period) << 32) | repeats,
            ));
            let offset = |index: usize| {
                i8::try_from(index as i32 - reach).expect("an allele within 127 repeats")
            };
            let (tracts, genotypes): (Vec<TractReads>, Vec<Vec<(i8, i8)>>) = (0..*loci)
                .map(|_| {
                    let frequencies = draw_locus_frequencies(&mut rng, &weights);
                    let (reads, genotypes) =
                        draw_locus_and_genotypes(&mut rng, &frequencies, &drawing);
                    let genotypes = genotypes
                        .into_iter()
                        .map(|(first, second)| (offset(first), offset(second)))
                        .collect();
                    (reads, genotypes)
                })
                .unzip();
            let evidence = StratumEvidence {
                stratum: truth.stratum,
                tracts,
                read_span,
                groups: 1,
                libraries: Arc::clone(&libraries),
                tracts_over_guard_threshold: 0,
                reads_reaching_not_crossing: 0,
                guard_reads: 0,
                bases_compared: 0,
                mismatching_bases: 0,
            };
            (evidence, truth, genotypes)
        })
        .collect();
    let mut strata = Vec::with_capacity(drawn.len());
    let mut truth = Vec::with_capacity(drawn.len());
    let mut genotypes = Vec::with_capacity(drawn.len());
    for (evidence, class_truth, class_genotypes) in drawn {
        strata.push(evidence);
        truth.push(class_truth);
        genotypes.push(class_genotypes);
    }
    DrawnPeriod {
        strata,
        truth,
        homozygote_excess,
        genotypes,
    }
}

/// **Write the truth, one row a class**: every number, the between-locus share, the shares as a
/// fit at `allele_span` sees them (`share_<offset>`), and the shares as drawn (`true_share_<offset>`).
///
/// # Errors
///
/// When the writer fails.
pub fn write_truth_tsv(
    truth: &[ClassTruth],
    allele_span: i32,
    mut out: impl Write,
) -> std::io::Result<()> {
    let Some(first) = truth.first() else {
        return Ok(());
    };
    let reach = first.true_allele_reach;
    let mut header: Vec<String> = [
        "period",
        "repeats",
        "loci",
        "level",
        "shorter_share",
        "fall_off",
        "concentration",
        "between_locus_share",
    ]
    .map(str::to_owned)
    .to_vec();
    header.extend((-allele_span..=allele_span).map(|offset| format!("share_{offset}")));
    header.extend((-reach..=reach).map(|offset| format!("true_share_{offset}")));
    writeln!(out, "{}", header.join("\t"))?;
    for class in truth {
        let mut row = vec![
            class.stratum.period.to_string(),
            class.stratum.reference_repeats.to_string(),
            class.loci.to_string(),
            format!("{:.6}", class.slippage.level),
            format!("{:.6}", class.slippage.shorter_share),
            format!("{:.6}", class.slippage.fall_off),
            format!("{:.6}", class.concentration),
            format!("{:.6}", class.between_locus_share()),
        ];
        row.extend(
            class
                .shares_in_window(allele_span)
                .iter()
                .map(|share| format!("{share:.8}")),
        );
        row.extend(class.true_shares.iter().map(|share| format!("{share:.8}")));
        writeln!(out, "{}", row.join("\t"))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_curve_joins_its_knots_and_holds_flat_past_them() {
        let curve = ParameterCurve::Knots(vec![(8.0, 0.1), (10.0, 0.3)]);
        assert_eq!(curve.at(6), 0.1);
        assert!((curve.at(9) - 0.2).abs() < 1e-12);
        assert_eq!(curve.at(12), 0.3);
    }

    #[test]
    fn the_shape_gives_its_reference_share_and_split() {
        let model = CohortModel::tomato_like(1, StutterRegime::Pcr);
        let AlleleLengthVariationTruth::Shape(shape) = &model.allele_length_variation else {
            unreachable!()
        };
        let shares = shape.shares_at(12, 14);
        assert!((shares.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((shares[14] - 0.60).abs() < 1e-12);
        let shorter: f64 = shares[..14].iter().sum();
        // At 12 repeats the split is 0.50 + 0.10 × 3/7.
        let split = 0.5 + 0.1 * 3.0 / 7.0;
        assert!((shorter - 0.4 * split).abs() < 1e-12);
        // At 16 repeats, where the tails fall slowest, alleles past ±9 carry about 2.5% of the
        // chromosomes, so the window's edge classes are exercised.
        let long = shape.shares_at(16, 14);
        let beyond: f64 = long[..5].iter().chain(&long[24..]).sum();
        assert!((0.015..0.04).contains(&beyond), "{beyond}");
    }

    #[test]
    fn folding_into_the_window_puts_the_tails_in_the_edge_classes() {
        let truth = CohortModel::tomato_like(1, StutterRegime::Pcr).truth_at(16, 10);
        let folded = truth.shares_in_window(9);
        assert_eq!(folded.len(), 19);
        assert!((folded.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let beyond: f64 = truth.true_shares[..6].iter().sum();
        assert!((folded[0] - beyond).abs() < 1e-15);
    }

    /// **The Poisson draws have the mean and variance they should**, at the two ends of the depth
    /// range.
    #[test]
    fn poisson_draws_have_their_mean_and_variance() {
        for mean in [0.5, 3.0, 30.0, 300.0] {
            let mut rng = Rng(7);
            let draws: Vec<f64> = (0..200_000).map(|_| f64::from(rng.poisson(mean))).collect();
            let average = draws.iter().sum::<f64>() / draws.len() as f64;
            let variance = draws
                .iter()
                .map(|d| float::powi(d - average, 2))
                .sum::<f64>()
                / (draws.len() - 1) as f64;
            let standard_error = (mean / draws.len() as f64).sqrt();
            assert!(
                (average - mean).abs() < 5.0 * standard_error,
                "mean {mean}: drew {average}"
            );
            assert!(
                (variance / mean - 1.0).abs() < 0.03,
                "mean {mean}: variance {variance}"
            );
        }
    }

    /// **Each locus's frequencies have the spread a Dirichlet of the class's concentration
    /// gives**: across loci, the frequency of a length has mean `share` and variance
    /// `share (1 − share) / (1 + c)`.
    #[test]
    fn locus_frequencies_have_the_dirichlet_mean_and_variance() {
        let truth = CohortModel::tomato_like(1, StutterRegime::Pcr).truth_at(12, 1);
        let weights: Vec<f64> = truth
            .true_shares
            .iter()
            .map(|share| truth.concentration * share)
            .collect();
        let mut rng = Rng(3);
        let draws = 100_000;
        let frequencies: Vec<Vec<f64>> = (0..draws)
            .map(|_| draw_locus_frequencies(&mut rng, &weights))
            .collect();
        for class in [13_usize, 14, 15, 11] {
            let share = truth.true_shares[class];
            let values: Vec<f64> = frequencies.iter().map(|f| f[class]).collect();
            let mean = values.iter().sum::<f64>() / draws as f64;
            let variance =
                values.iter().map(|v| float::powi(v - mean, 2)).sum::<f64>() / (draws - 1) as f64;
            let expected = share * (1.0 - share) / (1.0 + truth.concentration);
            assert!(
                (mean - share).abs() < 5.0 * (expected / draws as f64).sqrt(),
                "class {class}: mean {mean} against {share}"
            );
            assert!(
                (variance / expected - 1.0).abs() < 0.05,
                "class {class}: variance {variance} against {expected}"
            );
        }
    }

    /// **The reads reproduce their truth**: over many loci, the share of reads in each bucket is
    /// the class's shares mixed through the stutter, `Σ_a share_a · P(bucket | a)` — each read's
    /// allele is, averaged over loci and genotypes, a draw from the class's shares. Checked with a
    /// second library at three times the stutter, whose reads must follow its own level.
    #[test]
    fn the_reads_reproduce_the_shares_through_the_stutter() {
        let mut model = CohortModel::tomato_like(2, StutterRegime::Pcr);
        model.repeat_counts = Some(vec![12]);
        model.loci = LociPerClass::Each(4_000);
        model.samples = 20;
        model.depth = 5.0;
        model.library_multipliers = vec![1.0, 3.0];
        let drawn = draw_period(&model, RECORDED_OFFSET_RANGE);
        let truth = &drawn.truth[0];
        let read_span = RECORDED_OFFSET_RANGE;
        let reach = truth.true_allele_reach;
        let buckets = (2 * read_span + 1) as usize;
        for (library, multiplier) in [(0_u32, 1.0), (1, 3.0)] {
            // Per locus, this library's reads in each bucket: reads at one locus share its
            // alleles, so the error is taken over loci, which are independent.
            let per_locus: Vec<Vec<f64>> = drawn.strata[0]
                .tracts
                .iter()
                .map(|tract| {
                    let mut counts_here = vec![0.0; buckets];
                    for sample in &tract.samples {
                        for (index, counts) in &sample.by_library {
                            if index.0 == library {
                                for (into, count) in counts_here.iter_mut().zip(counts) {
                                    *into += f64::from(*count);
                                }
                            }
                        }
                    }
                    counts_here
                })
                .collect();
            let observed: Vec<f64> = (0..buckets)
                .map(|bucket| per_locus.iter().map(|locus| locus[bucket]).sum())
                .collect();
            let reads: f64 = observed.iter().sum();
            let slippage = Slippage {
                level: library_level(truth.slippage.level, multiplier),
                ..truth.slippage
            };
            let mut expected = vec![0.0; observed.len()];
            for (index, share) in truth.true_shares.iter().enumerate() {
                let from = slippage.read_probabilities(index as i32 - reach, read_span);
                for (into, probability) in expected.iter_mut().zip(from) {
                    *into += share * probability;
                }
            }
            for (bucket, (seen, probability)) in observed.iter().zip(&expected).enumerate() {
                let share = seen / reads;
                // The ratio estimator's error over independent loci.
                let spread: f64 = per_locus
                    .iter()
                    .map(|locus| {
                        let at_locus: f64 = locus.iter().sum();
                        float::powi(locus[bucket] - share * at_locus, 2)
                    })
                    .sum();
                let standard_error = spread.sqrt() / reads;
                assert!(
                    (share - probability).abs() < 5.0 * standard_error.max(1e-5),
                    "library {library}, bucket {bucket}: {share} against {probability}, \
                     error {standard_error}"
                );
            }
        }
    }

    #[test]
    fn a_draw_is_deterministic_and_a_class_does_not_depend_on_its_neighbours() {
        let mut model = CohortModel::tomato_like(1, StutterRegime::Pcr);
        model.loci = LociPerClass::Each(30);
        model.samples = 5;
        model.repeat_counts = Some(vec![9, 10, 11]);
        let all = draw_period(&model, RECORDED_OFFSET_RANGE);
        model.repeat_counts = Some(vec![10]);
        let one = draw_period(&model, RECORDED_OFFSET_RANGE);
        assert_eq!(all.strata[1], one.strata[0]);
        assert_eq!(all.homozygote_excess, one.homozygote_excess);
    }

    /// **[`draw_stratum`] draws what it drew before the simulator was built around it**: a
    /// checksum of one draw, taken before the change.
    #[test]
    fn draw_stratum_draws_what_it_always_drew() {
        let drawn = draw_stratum(
            Slippage {
                level: 0.05,
                shorter_share: 0.6,
                fall_off: 0.3,
            },
            &spectrum_of(5),
            0.8,
            0.3,
            40,
            7,
            4,
            2,
            11,
        );
        assert_eq!(checksum(&drawn), DRAW_STRATUM_CHECKSUM);
    }

    const DRAW_STRATUM_CHECKSUM: u64 = 4_336_400_559_775_743_465;

    fn checksum(evidence: &StratumEvidence) -> u64 {
        let mut hash: u64 = 0xCBF2_9CE4_8422_2325;
        let mut mix = |value: u64| {
            hash ^= value;
            hash = hash.wrapping_mul(0x0100_0000_01B3);
        };
        for tract in &evidence.tracts {
            mix(tract.samples.len() as u64);
            for sample in &tract.samples {
                mix(u64::from(sample.sample));
                for (library, counts) in &sample.by_library {
                    mix(u64::from(library.0));
                    for count in counts {
                        mix(u64::from(*count));
                    }
                }
            }
        }
        hash
    }
}
