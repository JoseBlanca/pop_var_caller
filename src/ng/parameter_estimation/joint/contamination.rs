//! How much of a sample's DNA came from a different plant.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_fit.md` §3.4. Measured before it was
//! built: `doc/devel/ng/reports/joint_contamination_2026-08-12.md`, whose numbers this module
//! is a port of rather than a reinvention.
//!
//! # What contamination is, and what makes it hard
//!
//! A second seedling in the tube, or a neighbouring library on the same sequencing run, puts a
//! small share of another individual's reads into this sample's. It is **invisible wherever the
//! two plants carry the same allele**, and where they differ it shows as a *small* share of
//! reads carrying the other allele. That smallness is what tells it from a heterozygote, whose
//! two alleles are balanced.
//!
//! To call a share of reads "unexpected" the estimator has to know what to expect, which is
//! **the allele frequency of the population the contaminant came from**. Use one frequency for
//! the whole panel and a diverged accession's own alleles look unexpected everywhere — and the
//! measured consequence is not a false alarm but the opposite:
//!
//! > On a panel of four subpopulations at `F_st` 0.20, a pooled frequency returns **0.005 for a
//! > sample truly contaminated at 3%**, and exactly zero for one at 1%. Both pass any threshold
//! > as clean. **Structure does not invent contamination; it hides it.**
//!
//! # What is done instead, and the one thing not to do
//!
//! **Every sample gets its own allele frequency at every position, and the panel is never
//! partitioned to get it.** Splitting the panel into groups and estimating within each is worse
//! than ignoring structure altogether — it adds about 0.015 to every sample's estimate and puts
//! 41 to 47 of 50 clean samples over a 1% threshold — because twelve samples are too few to
//! estimate a frequency from, and a *noisy* frequency manufactures contamination for the mirror
//! image of the reason a *smoothly wrong* one hides it.
//!
//! So each sample gets a few coordinates in the cohort's own axes of variation, and at each
//! position the allele frequency is fitted as **a straight line in those coordinates, using
//! every sample**. A sample alone at one end of an axis gets its own frequency, and what it
//! borrows from the panel is the *slope* — how fast frequency changes along the axis — never a
//! neighbour's allele counts. There is no threshold anywhere, so an admixed accession simply
//! gets an intermediate frequency.
//!
//! # The refusal, and why it is free
//!
//! **A sample sitting alone at the end of an axis has a fitted frequency that is mostly its own
//! echo**, and by the mechanism above that manufactures contamination: on a panel of 40, 5, 3
//! and 2 accessions with nobody contaminated, the group of two returns a spurious 0.031.
//!
//! How much of its own fitted frequency a sample supplies depends **only on the coordinates**,
//! so it is one number per sample for the whole run and can be computed before a single position
//! is fitted. It tracked the damage exactly on that panel — 0.027, 0.307, 0.429 and 0.857 across
//! the four groups. A sample supplying more than about half of its own frequency is told *not
//! identified* rather than given a number.
//!
//! # What this finds, and what it does not — measured 2026-08-13
//!
//! On forty samples in four subpopulations at `F_st` 0.20 and three reads a position, with one
//! sample contaminated at 3%:
//!
//! | | 12,000 varying positions | 60,000 |
//! |---|---:|---:|
//! | the contaminated sample | 0.0166 | 0.0163 |
//! | the worst of the thirty-nine clean ones | 0.0032 | 0.0004 |
//!
//! **It finds the sample and it understates the fraction.** The separation is what a threshold
//! needs and it improves with positions — forty times the noise floor at 60,000 — but the value
//! itself sits at about half the truth and **does not move with more positions**, so it is a
//! bias rather than noise.
//!
//! **The cause is named in the specification and is not implemented here.** The contaminated
//! sample's own coordinates are estimated from its own reads, which the contamination has
//! already pulled towards the panel average; its fitted frequency therefore sits closer to the
//! contaminant's than it should, and the difference the estimator lives on shrinks.
//! `verifyBamID2` maximises over `α` **and the intended sample's coordinates together** for
//! exactly this reason. Until that is done, read `α` as *this sample stands out from the panel*
//! and not as *this sample is 1.7% contaminated*.

use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;

use super::records::{DepthCode, SampleRecords};

/// The default number of axes of variation each sample's frequency is a line in.
///
/// **Four by inheritance from `verifyBamID2`, and the specification says so is not an
/// argument** — how many a plant panel needs is set by the panel. It is a parameter for that
/// reason.
pub const DEFAULT_COMPONENTS: usize = 4;

/// How much of its own fitted frequency a sample may supply before its estimate is refused.
///
/// A sample pulling its fair share supplies `(components + 1) / samples` — 0.10 on a panel of
/// fifty with four axes. At 0.857 the fitted line at that sample's position is that sample's
/// own reading, and the estimate it produces is a reading of its own noise.
pub const MAX_LEVERAGE: f64 = 0.5;

/// A position the cohort varies at, and every sample's reads there.
///
/// **Only positions that vary carry information about who differs from whom**, and at a
/// position where the population carries one allele the two-genotype mixture has nothing to
/// separate.
struct Marker {
    /// Reads on the cohort's most common non-reference allele, per sample.
    alternative: Vec<u32>,
    depth: Vec<u32>,
    /// The panel-wide frequency, with the error rate inverted out of the read share.
    pooled: f64,
    /// The posterior mean number of alternative copies, per sample, under a prior at `pooled`.
    ///
    /// **A raw read fraction is far too noisy to decompose at three reads a position**; the
    /// prior is what makes it usable, and it is the same bootstrap `PCAngsd` starts from.
    dosage: Vec<f64>,
}

/// A position enters only if this many samples put a read on it.
const MIN_SAMPLES_WITH_DATA: usize = 8;

/// …and only if the panel-wide frequency is inside this band.
const MIN_FREQUENCY: f64 = 0.02;

/// The fraction of a sample's reads that came from another individual, or the reason there is
/// no number.
///
/// **An enum and not an `Option<f64>`**: *not identified* and *zero* are different answers, and
/// a caller told "no contamination" would act on it.
#[derive(Clone, PartialEq, Debug)]
pub enum ContaminationEstimate {
    Estimated {
        /// **Restricted to a half.** The likelihood reads sequence only, so it is symmetric: it
        /// cannot tell a sample 20% contaminated from one 80% contaminated, and a swap of two
        /// samples is invisible to it by construction. A number above a half would not be a
        /// stronger claim but a mirror image.
        alpha: f64,
        /// How many positions stood behind it — what says how far to trust the value.
        markers: u64,
        /// How much of its own fitted allele frequency this sample supplied, against the
        /// `(components + 1) / samples` a sample pulling its fair share would.
        leverage: f64,
    },
    NotIdentified {
        reason: NotIdentifiedReason,
    },
}

#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum NotIdentifiedReason {
    /// Fewer than two samples: there is no panel, so there is no frequency to be surprised by.
    NoPanel,
    /// Too few positions where the cohort varies.
    TooFewMarkers,
    /// This sample supplies most of its own fitted frequency, so an estimate from it would be a
    /// reading of its own noise (`MAX_LEVERAGE`).
    OwnFrequencyIsItsOwnEcho,
}

impl std::fmt::Display for NotIdentifiedReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let text = match self {
            Self::NoPanel => "there is no panel to compare against",
            Self::TooFewMarkers => "too few positions where the cohort varies",
            Self::OwnFrequencyIsItsOwnEcho => {
                "this sample supplies most of its own fitted allele frequency"
            }
        };
        f.write_str(text)
    }
}

/// Every sample's contamination fraction, in the order `samples` are given.
///
/// `error_rate` is per sample — the rate at which a read misreads a base, which the fit has
/// already produced — and `hom_excess` how much less heterozygous each sample is than random
/// mating in the panel predicts.
pub fn fit_contamination(
    samples: &[SampleRecords],
    edges: &DepthBinEdges,
    error_rate: &[f64],
    hom_excess: &[f64],
    components: usize,
) -> Vec<ContaminationEstimate> {
    let count = samples.len();
    if count < 2 {
        return vec![
            ContaminationEstimate::NotIdentified {
                reason: NotIdentifiedReason::NoPanel,
            };
            count
        ];
    }
    let mean_error = error_rate.iter().sum::<f64>() / count as f64;
    let markers = markers(samples, edges, mean_error);
    if markers.len() < 100 {
        return vec![
            ContaminationEstimate::NotIdentified {
                reason: NotIdentifiedReason::TooFewMarkers,
            };
            count
        ];
    }

    let components = components.min(count.saturating_sub(2)).max(1);
    let coordinates = ancestry_coordinates(&markers, count, components);
    let leverage = coordinate_leverage(&coordinates, components);
    let frequencies = fitted_frequencies(&markers, &coordinates, components);

    (0..count)
        .map(|sample| {
            if leverage[sample] > MAX_LEVERAGE {
                return ContaminationEstimate::NotIdentified {
                    reason: NotIdentifiedReason::OwnFrequencyIsItsOwnEcho,
                };
            }
            ContaminationEstimate::Estimated {
                alpha: fit_alpha(
                    &markers,
                    &frequencies,
                    sample,
                    hom_excess[sample],
                    error_rate[sample],
                ),
                markers: markers.len() as u64,
                leverage: leverage[sample],
            }
        })
        .collect()
}

/// The positions the cohort varies at, with each sample's reads and its dosage there.
fn markers(samples: &[SampleRecords], edges: &DepthBinEdges, error: f64) -> Vec<Marker> {
    let count = samples.len();
    let positions = samples
        .first()
        .and_then(|s| s.generic.values().next())
        .map_or(0, |g| g.depth().len());

    // Which non-reference allele the cohort carries at each position: the one the most reads
    // across the whole panel fell on.
    let mut totals = vec![[0_u32; 5]; positions];
    for sample in samples {
        for group in sample.generic.values() {
            for observation in group.non_reference() {
                totals[observation.index as usize][usize::from(observation.allele.code())] +=
                    observation.reads;
            }
        }
    }
    let major: Vec<u8> = totals
        .iter()
        .map(|counts| {
            counts
                .iter()
                .take(4)
                .enumerate()
                .max_by_key(|(_, n)| **n)
                .map(|(allele, _)| allele as u8)
                .unwrap_or(0)
        })
        .collect();

    let mut alternative = vec![vec![0_u32; positions]; count];
    let mut depth = vec![vec![0_u32; positions]; count];
    for (s, sample) in samples.iter().enumerate() {
        for group in sample.generic.values() {
            for (index, code) in group.depth().iter().enumerate() {
                if let DepthCode::Binned(bin) = code {
                    depth[s][index] = depth[s][index]
                        .saturating_add(edges.representative_depth(bin).round() as u32);
                }
            }
            for observation in group.non_reference() {
                let index = observation.index as usize;
                if observation.allele.code() == major[index] {
                    alternative[s][index] = alternative[s][index].saturating_add(observation.reads);
                }
            }
        }
    }

    let mut out = Vec::new();
    for index in 0..positions {
        let covered = (0..count).filter(|&s| depth[s][index] > 0).count();
        if covered < MIN_SAMPLES_WITH_DATA {
            continue;
        }
        let (alt, total): (u64, u64) = (0..count).fold((0, 0), |(a, d), s| {
            (
                a + u64::from(alternative[s][index]),
                d + u64::from(depth[s][index]),
            )
        });
        if total == 0 {
            continue;
        }
        let share = alt as f64 / total as f64;
        // The read share is the frequency blurred by the error rate; inverting it is what a
        // spectrum fitted over the same positions converges to.
        let pooled = ((share - error) / (1.0 - 2.0 * error)).clamp(1e-3, 1.0 - 1e-3);
        if !(MIN_FREQUENCY..=1.0 - MIN_FREQUENCY).contains(&pooled) {
            continue;
        }
        let priors = genotype_priors(pooled, 0.0);
        let dosage: Vec<f64> = (0..count)
            .map(|s| {
                let mut weights = [0.0_f64; 3];
                for (copies, weight) in weights.iter_mut().enumerate() {
                    *weight = priors[copies]
                        * binomial(
                            alternative[s][index],
                            depth[s][index],
                            read_share(copies, error),
                        );
                }
                let total: f64 = weights.iter().sum();
                if total > 0.0 {
                    (weights[1] + 2.0 * weights[2]) / total
                } else {
                    2.0 * pooled
                }
            })
            .collect();
        out.push(Marker {
            alternative: (0..count).map(|s| alternative[s][index]).collect(),
            depth: (0..count).map(|s| depth[s][index]).collect(),
            pooled,
            dosage,
        });
    }
    out
}

/// Each sample's coordinates in the cohort's own axes of variation.
///
/// The samples' similarity matrix on dosages centred and scaled the way a population-structure
/// analysis scales them — divided by `√(p(1−p))`, so a rare allele's contribution is not
/// swamped by a common one — and then its leading eigenvectors. **A few numbers per sample and
/// nothing else is kept: no thresholds on the axes, and no groups.**
fn ancestry_coordinates(markers: &[Marker], samples: usize, components: usize) -> Vec<Vec<f64>> {
    let mut gram = vec![0.0_f64; samples * samples];
    let mut centred = vec![0.0_f64; samples];
    for marker in markers {
        let mean = 2.0 * marker.pooled;
        let scale = (marker.pooled * (1.0 - marker.pooled)).sqrt().max(1e-6);
        for (slot, dosage) in centred.iter_mut().zip(&marker.dosage) {
            *slot = (dosage - mean) / scale;
        }
        for i in 0..samples {
            for j in i..samples {
                gram[i * samples + j] += centred[i] * centred[j];
            }
        }
    }
    for i in 0..samples {
        for j in 0..i {
            gram[i * samples + j] = gram[j * samples + i];
        }
    }
    leading_eigenvectors(&gram, samples, components)
}

/// **How much of its own fitted allele frequency each sample supplies.**
///
/// The line at every position is fitted against the same coordinates, so this is one number per
/// sample for the whole run and it can be computed before a single position is touched. It runs
/// from `(components + 1) / samples` — a sample pulling its fair share — up towards one, where
/// the line at that sample's position is determined by that sample alone.
fn coordinate_leverage(coordinates: &[Vec<f64>], components: usize) -> Vec<f64> {
    let samples = coordinates.len();
    let width = components + 1;
    let design = |sample: usize| {
        let mut row = vec![1.0_f64; width];
        row[1..].copy_from_slice(&coordinates[sample][..components]);
        row
    };
    let mut xtx = vec![0.0_f64; width * width];
    for sample in 0..samples {
        let row = design(sample);
        for a in 0..width {
            for b in 0..width {
                xtx[a * width + b] += row[a] * row[b];
            }
        }
    }
    for a in 0..width {
        xtx[a * width + a] += 1e-9;
    }
    (0..samples)
        .map(|sample| {
            let row = design(sample);
            solve(&xtx, &row, width).map_or(1.0, |z| row.iter().zip(&z).map(|(x, z)| x * z).sum())
        })
        .collect()
}

/// One allele frequency per sample per position, from a straight line in the coordinates.
///
/// **The slopes are shrunk.** A position whose slopes are indistinguishable from noise keeps
/// only its intercept — the panel-wide frequency — so modelling structure is never worse than
/// not modelling it. Unshrunk, the same fit returned 0.0443 where the truth was 0.030, and on an
/// unbalanced panel it walked to the search boundary.
fn fitted_frequencies(
    markers: &[Marker],
    coordinates: &[Vec<f64>],
    components: usize,
) -> Vec<Vec<f64>> {
    let samples = coordinates.len();
    let width = components + 1;
    let design: Vec<Vec<f64>> = (0..samples)
        .map(|sample| {
            let mut row = vec![1.0_f64; width];
            row[1..].copy_from_slice(&coordinates[sample][..components]);
            row
        })
        .collect();
    let mut xtx = vec![0.0_f64; width * width];
    for row in &design {
        for a in 0..width {
            for b in 0..width {
                xtx[a * width + b] += row[a] * row[b];
            }
        }
    }
    // A ridge term, small beside the samples on the diagonal, so a position whose dosages are
    // constant cannot make the system singular.
    for a in 0..width {
        xtx[a * width + a] += 1e-6;
    }

    markers
        .iter()
        .map(|marker| {
            let mut xty = vec![0.0_f64; width];
            for (sample, row) in design.iter().enumerate() {
                for a in 0..width {
                    xty[a] += row[a] * marker.dosage[sample];
                }
            }
            let mut beta = solve(&xtx, &xty, width).unwrap_or_else(|| {
                let mut fallback = vec![0.0; width];
                fallback[0] = 2.0 * marker.pooled;
                fallback
            });
            let mean: f64 = marker.dosage.iter().sum::<f64>() / samples as f64;
            if components > 0 && samples > width {
                // The positive-part James–Stein factor: how much of the dosage spread the line
                // explains, against how much noise alone would explain.
                let (mut explained, mut residual) = (0.0, 0.0);
                for (sample, row) in design.iter().enumerate() {
                    let fitted: f64 = row.iter().zip(&beta).map(|(x, b)| x * b).sum();
                    explained += (fitted - mean) * (fitted - mean);
                    residual += (marker.dosage[sample] - fitted) * (marker.dosage[sample] - fitted);
                }
                let noise = residual / (samples - width) as f64;
                let factor = if explained > 0.0 {
                    (1.0 - components as f64 * noise / explained).max(0.0)
                } else {
                    0.0
                };
                for slope in beta.iter_mut().skip(1) {
                    *slope *= factor;
                }
                // The intercept takes back whatever the slopes gave up.
                let centre: f64 = design
                    .iter()
                    .map(|row| {
                        row[1..]
                            .iter()
                            .zip(&beta[1..])
                            .map(|(x, b)| x * b)
                            .sum::<f64>()
                    })
                    .sum();
                beta[0] = mean - centre / samples as f64;
            }
            design
                .iter()
                .map(|row| {
                    let fitted: f64 = row.iter().zip(&beta).map(|(x, b)| x * b).sum();
                    (fitted / 2.0).clamp(1e-4, 1.0 - 1e-4)
                })
                .collect()
        })
        .collect()
}

/// One sample's contamination fraction.
///
/// **The two-genotype mixture, reading sequence only**: a read comes from this sample's own
/// genotype with probability `1 − α` and from the contaminant's with probability `α`, and
/// neither genotype is known, so both are summed over.
fn fit_alpha(
    markers: &[Marker],
    frequencies: &[Vec<f64>],
    sample: usize,
    hom_excess: f64,
    error: f64,
) -> f64 {
    let score = |alpha: f64| -> f64 {
        markers
            .iter()
            .zip(frequencies)
            .map(|(marker, per_sample)| {
                ln_marker(
                    marker.alternative[sample],
                    marker.depth[sample],
                    per_sample[sample],
                    marker.pooled,
                    hom_excess,
                    error,
                    alpha,
                )
            })
            .sum()
    };
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
    // The bracket never reaches zero, and zero is the answer that matters most: a clean sample
    // must be able to come back clean rather than at the smallest value the search can express.
    if score(0.0) >= score(best) { 0.0 } else { best }
}

/// `ln P(this sample's reads here | α)`, both genotypes summed over.
///
/// **The two genotypes are drawn against two different frequencies, and that is the design
/// rather than an oversight.** The sample's own genotype is drawn at *its* frequency — the line
/// fitted at its own coordinates — because that is a statement about its ancestry. The
/// contaminant's is drawn at the frequency of **whoever else was sequenced beside it**, which by
/// default is the whole panel, because a neighbouring library on a plate is not chosen for
/// ancestry (spec §3.4.3). Scoring both against the sample's own frequency inflates `α` on a
/// diverged panel: a contaminant from a different group carries alleles the sample's own
/// frequency calls rare, and rare alleles turning up is the contamination signature.
fn ln_marker(
    alternative: u32,
    depth: u32,
    own_frequency: f64,
    batch_frequency: f64,
    hom_excess: f64,
    error: f64,
    alpha: f64,
) -> f64 {
    if depth == 0 {
        return 0.0;
    }
    let priors = genotype_priors(own_frequency, hom_excess);
    // The contaminant is somebody else, so it carries no inbreeding of this sample's.
    let contaminant = genotype_priors(batch_frequency, 0.0);
    let mut terms = [f64::NEG_INFINITY; 9];
    for own in 0..3 {
        for other in 0..3 {
            let weight = priors[own] * contaminant[other];
            if weight <= 0.0 {
                continue;
            }
            let p = (1.0 - alpha) * read_share(own, error) + alpha * read_share(other, error);
            terms[own * 3 + other] =
                weight.ln() + binomial(alternative, depth, p).max(f64::MIN_POSITIVE).ln();
        }
    }
    ln_sum_exp(&terms)
}

/// The chance one read shows the alternative allele from a genotype carrying `copies` of it.
fn read_share(copies: usize, error: f64) -> f64 {
    match copies {
        0 => error,
        1 => 0.5,
        _ => 1.0 - error,
    }
}

fn genotype_priors(frequency: f64, hom_excess: f64) -> [f64; 3] {
    let (p, q) = (frequency, 1.0 - frequency);
    [
        q * q + hom_excess * p * q,
        2.0 * p * q * (1.0 - hom_excess),
        p * p + hom_excess * p * q,
    ]
}

/// `P(k of n reads showed the allele)`, without the binomial coefficient — it is the same for
/// every genotype pair and every `α`, so it cancels out of both the sum and the search.
fn binomial(k: u32, n: u32, p: f64) -> f64 {
    let p = p.clamp(1e-12, 1.0 - 1e-12);
    p.powi(k as i32) * (1.0 - p).powi((n - k.min(n)) as i32)
}

fn ln_sum_exp(values: &[f64]) -> f64 {
    let largest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if largest == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    largest + values.iter().map(|v| (v - largest).exp()).sum::<f64>().ln()
}

/// The `wanted` eigenvectors of largest eigenvalue of a small symmetric matrix, by cyclic
/// Jacobi rotation. One row per sample, so this is fifty by fifty at most.
fn leading_eigenvectors(matrix: &[f64], n: usize, wanted: usize) -> Vec<Vec<f64>> {
    let mut a = matrix.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..n)
            .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
            .map(|(i, j)| a[i * n + j] * a[i * n + j])
            .sum();
        if off < 1e-18 {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                if a[p * n + q].abs() < 1e-15 {
                    continue;
                }
                let theta = (a[q * n + q] - a[p * n + p]) / (2.0 * a[p * n + q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..n {
                    let (akp, akq) = (a[k * n + p], a[k * n + q]);
                    a[k * n + p] = c * akp - s * akq;
                    a[k * n + q] = s * akp + c * akq;
                }
                for k in 0..n {
                    let (apk, aqk) = (a[p * n + k], a[q * n + k]);
                    a[p * n + k] = c * apk - s * aqk;
                    a[q * n + k] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let (vkp, vkq) = (v[k * n + p], v[k * n + q]);
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        a[j * n + j]
            .partial_cmp(&a[i * n + i])
            .expect("no NaN in a similarity matrix")
    });
    (0..n)
        .map(|sample| {
            order[..wanted.min(n)]
                .iter()
                .map(|&axis| v[sample * n + axis])
                .collect()
        })
        .collect()
}

/// Gaussian elimination with partial pivoting on a `width × width` system. `None` where it is
/// singular.
fn solve(a: &[f64], b: &[f64], width: usize) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut y = b.to_vec();
    for column in 0..width {
        let pivot = (column..width)
            .max_by(|&i, &j| {
                m[i * width + column]
                    .abs()
                    .partial_cmp(&m[j * width + column].abs())
                    .expect("no NaN in normal equations")
            })
            .expect("at least one row");
        if m[pivot * width + column].abs() < 1e-12 {
            return None;
        }
        for k in 0..width {
            m.swap(column * width + k, pivot * width + k);
        }
        y.swap(column, pivot);
        for row in (column + 1)..width {
            let factor = m[row * width + column] / m[column * width + column];
            for k in column..width {
                m[row * width + k] -= factor * m[column * width + k];
            }
            y[row] -= factor * y[column];
        }
    }
    let mut out = vec![0.0; width];
    for row in (0..width).rev() {
        let mut value = y[row];
        for k in (row + 1)..width {
            value -= m[row * width + k] * out[k];
        }
        out[row] = value / m[row * width + row];
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **`(components + 1) / samples` is the panel's *mean* leverage, not every sample's.**
    /// Twenty samples spread evenly along one axis average exactly 0.10 and range from 0.05 at
    /// the middle to 0.19 at the ends — a sample at the end of an axis always supplies more of
    /// its own line, and that is a property of a straight line rather than a defect. The
    /// refusal is set well above the ends of an evenly filled axis for that reason.
    #[test]
    fn an_evenly_spread_panel_averages_its_fair_share_and_none_of_it_is_refused() {
        let coordinates: Vec<Vec<f64>> = (0..20).map(|s| vec![(s as f64 - 9.5) / 10.0]).collect();
        let leverage = coordinate_leverage(&coordinates, 1);
        let mean: f64 = leverage.iter().sum::<f64>() / leverage.len() as f64;
        let fair = 2.0 / 20.0;
        assert!((mean - fair).abs() < 1e-6, "{mean} against {fair}");
        let largest = leverage.iter().copied().fold(0.0_f64, f64::max);
        assert!(
            largest < MAX_LEVERAGE,
            "the worst placed sample is at {largest}"
        );
        assert!(
            largest > mean,
            "an end of the axis supplies more than the average"
        );
    }

    /// **The refusal's own oracle.** A sample sitting alone at the end of an axis supplies most
    /// of its own fitted frequency, and that is what the refusal is keyed on.
    #[test]
    fn a_sample_alone_at_the_end_of_an_axis_supplies_most_of_its_own_frequency() {
        let mut coordinates: Vec<Vec<f64>> = (0..19).map(|_| vec![0.0]).collect();
        coordinates.push(vec![1.0]);
        let leverage = coordinate_leverage(&coordinates, 1);
        assert!(
            leverage[19] > MAX_LEVERAGE,
            "the lone sample supplies {} of its own line",
            leverage[19]
        );
        assert!(
            leverage[0] < 0.2,
            "the nineteen together supply {} each",
            leverage[0]
        );
    }

    #[test]
    fn the_mixture_is_flat_where_the_two_genotypes_agree() {
        // Both plants homozygous reference: no share of reads from the other one shows.
        let clean = ln_marker(0, 20, 0.001, 0.001, 0.0, 0.002, 0.0);
        let contaminated = ln_marker(0, 20, 0.001, 0.001, 0.0, 0.002, 0.2);
        assert!(
            (clean - contaminated).abs() < 0.05,
            "contamination is invisible where the two plants carry the same allele"
        );
    }

    #[test]
    fn a_small_share_of_stray_reads_scores_better_contaminated_than_clean() {
        // Two reads in fifty on the other allele, at a position where the allele is common:
        // too few for a heterozygote, too many for the error rate.
        let clean = ln_marker(2, 50, 0.3, 0.3, 0.0, 0.002, 0.0);
        let contaminated = ln_marker(2, 50, 0.3, 0.3, 0.0, 0.002, 0.04);
        assert!(
            contaminated > clean,
            "clean {clean}, contaminated {contaminated}"
        );
    }

    // ---- the whole estimator, on a panel drawn with known contamination ----------

    use crate::ng::parameter_estimation::joint::loci::{
        CatalogBuildSettings, KeptLociDigester, ReferenceDigest, RegionSetDigest, SelectionIdentity,
    };
    use crate::ng::parameter_estimation::joint::records::{
        AlleleObservation, DepthCap, DepthLadderDigest, GenericRecords, ObservedAllele,
        PackedDepthCodes, ReadCap, RecordIdentity,
    };
    use crate::ng::repeat_catalog::StrRepeatCriteria;
    use crate::ng::tandem_repeat::ScanParams;
    use crate::ng::types::ReadGroupId;
    use std::collections::BTreeMap;

    struct Rng(u64);

    impl Rng {
        fn uniform(&mut self) -> f64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 >> 11) as f64 / (1_u64 << 53) as f64
        }

        fn gamma(&mut self, shape: f64) -> f64 {
            if shape < 1.0 {
                let u = self.uniform().max(1e-12);
                return self.gamma(shape + 1.0) * u.powf(1.0 / shape);
            }
            let d = shape - 1.0 / 3.0;
            let c = 1.0 / (9.0 * d).sqrt();
            loop {
                let u1 = self.uniform().max(1e-12);
                let u2 = self.uniform();
                let z = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
                let v = (1.0 + c * z).powi(3);
                if v <= 0.0 {
                    continue;
                }
                if self.uniform().max(1e-12).ln() < 0.5 * z * z + d - d * v + d * v.ln() {
                    return d * v;
                }
            }
        }

        fn beta(&mut self, a: f64, b: f64) -> f64 {
            let x = self.gamma(a);
            let y = self.gamma(b);
            (x / (x + y)).clamp(1e-4, 1.0 - 1e-4)
        }

        fn binomial(&mut self, n: u32, p: f64) -> u32 {
            (0..n).filter(|_| self.uniform() < p).count() as u32
        }

        fn poisson(&mut self, mean: f64) -> u32 {
            let limit = (-mean).exp();
            let mut product = self.uniform();
            let mut count = 0;
            while product > limit && count < 200 {
                count += 1;
                product *= self.uniform();
            }
            count
        }
    }

    /// A panel of `groups` subpopulations diverged by `fst`, one sample of which is
    /// contaminated at `alpha` by a plant drawn from the whole panel.
    ///
    /// **The contaminant is drawn from the panel and not from the contaminated sample's own
    /// subpopulation**, because that is the harder case and the one a sequencing batch
    /// produces: a neighbouring library on the same run is whoever else was on the plate.
    fn structured_panel(
        samples: usize,
        markers: usize,
        depth: f64,
        groups: usize,
        fst: f64,
        spiked: Option<(usize, f64)>,
        seed: u64,
    ) -> Vec<SampleRecords> {
        const ERROR: f64 = 0.002;
        let mut rng = Rng(seed);
        let edges = DepthBinEdges::new();
        let mut codes: Vec<PackedDepthCodes> = (0..samples)
            .map(|_| PackedDepthCodes::never_walked(markers))
            .collect();
        let mut sparse: Vec<Vec<AlleleObservation>> = vec![Vec::new(); samples];
        let group_of = |s: usize| s * groups / samples;

        for index in 0..markers {
            let ancestral = 0.05 + 0.90 * rng.uniform();
            let scale = (1.0 - fst) / fst.max(1e-9);
            let per_group: Vec<f64> = (0..groups)
                .map(|_| {
                    if fst <= 0.0 {
                        ancestral
                    } else {
                        rng.beta(ancestral * scale, (1.0 - ancestral) * scale)
                    }
                })
                .collect();
            let genotypes: Vec<u32> = (0..samples)
                .map(|s| rng.binomial(2, per_group[group_of(s)]))
                .collect();
            for s in 0..samples {
                let reads = rng.poisson(depth);
                let mut alternative = 0;
                for _ in 0..reads {
                    let from = match spiked {
                        Some((which, alpha)) if which == s && rng.uniform() < alpha => {
                            genotypes[(rng.uniform() * samples as f64) as usize % samples]
                        }
                        _ => genotypes[s],
                    };
                    if rng.uniform() < read_share(from as usize, ERROR) {
                        alternative += 1;
                    }
                }
                codes[s].set(index, DepthCode::Binned(edges.bin_for(reads)));
                if alternative > 0 {
                    sparse[s].push(AlleleObservation {
                        index: index as u32,
                        allele: ObservedAllele::C,
                        reads: alternative,
                    });
                }
            }
        }

        let identity = RecordIdentity {
            selection: SelectionIdentity {
                seed,
                reference: ReferenceDigest([7; 16]),
                analysed_regions: RegionSetDigest([9; 16]),
                catalog_built_under: CatalogBuildSettings {
                    criteria: StrRepeatCriteria::default(),
                    scan: ScanParams::default(),
                    tool_version: "0.1.0".to_string(),
                },
                ssr_criteria: StrRepeatCriteria::default(),
                generic_target: markers as u64,
                ssr_cap: 1_000,
            },
            kept_loci: KeptLociDigester::new().finish(),
            ssr_stratum_counts: Default::default(),
            read_cap: ReadCap(1_000),
            depth_ladder: DepthLadderDigest::of(&DepthBinEdges::new()),
            depth_cap: DepthCap(124),
            coverage_window: None,
        };
        (0..samples)
            .map(|s| SampleRecords {
                sample: format!("s{s:02}"),
                generic: [(
                    ReadGroupId(0),
                    GenericRecords::from_parts(
                        std::mem::replace(&mut codes[s], PackedDepthCodes::never_walked(0)),
                        std::mem::take(&mut sparse[s]),
                    ),
                )]
                .into_iter()
                .collect(),
                ssr: BTreeMap::new(),
                coverage: None,
                identity: identity.clone(),
            })
            .collect()
    }

    fn alphas(estimates: &[ContaminationEstimate]) -> Vec<f64> {
        estimates
            .iter()
            .map(|estimate| match estimate {
                ContaminationEstimate::Estimated { alpha, .. } => *alpha,
                ContaminationEstimate::NotIdentified { .. } => f64::NAN,
            })
            .collect()
    }

    /// **The test that says the whole thing works, and the one a pooled frequency fails.**
    /// A panel of four diverged subpopulations with one sample contaminated at 3%: that sample
    /// has to be found, and it has to stand clear of the cleanest-looking of the others.
    ///
    /// A clean structured panel is *passed* by an estimator using the pooled frequency, which
    /// returns zero for everybody — so the test that catches a broken frequency is a spiked
    /// panel whose contaminated sample must be found, not a clean one whose samples must all
    /// read zero (`joint_contamination_2026-08-12.md` §3).
    #[test]
    fn a_contaminated_sample_in_a_diverged_panel_is_found_and_stands_clear() {
        let samples = 40;
        let panel = structured_panel(samples, 12_000, 3.0, 4, 0.20, Some((0, 0.03)), 0x51ED2709);
        let estimates = fit_contamination(
            &panel,
            &DepthBinEdges::new(),
            &vec![0.002; samples],
            &vec![0.0; samples],
            4,
        );
        let alpha = alphas(&estimates);
        let spiked = alpha[0];
        let worst_clean = alpha[1..].iter().copied().fold(0.0_f64, f64::max);
        eprintln!("spiked at 0.030 came back {spiked:.4}; worst of the 39 clean {worst_clean:.4}");
        // **Detection, not magnitude** — see the module doc's note on attenuation. Measured
        // here: 0.0166 for a truth of 0.030, against a worst clean sample of 0.0032, and the
        // estimate does not move with more markers (0.0163 at 60,000) while the floor does
        // (0.0004). So what this asserts is that the contaminated sample is found and stands
        // clear, which is what a threshold has to do.
        assert!(
            spiked > 0.012,
            "the 3% sample came back at {spiked}, which no sensible threshold would catch"
        );
        assert!(
            spiked > 4.0 * worst_clean.max(1e-4),
            "the 3% sample at {spiked} does not stand clear of the worst clean one at \
             {worst_clean}"
        );
    }

    /// The other half of the pair: **nobody contaminated must come back as somebody
    /// contaminated.** A false positive means telling a lab to repeat a sequencing run.
    #[test]
    fn a_clean_diverged_panel_flags_nobody_badly() {
        let samples = 40;
        let panel = structured_panel(samples, 12_000, 3.0, 4, 0.20, None, 0xA5A5_1234);
        let estimates = fit_contamination(
            &panel,
            &DepthBinEdges::new(),
            &vec![0.002; samples],
            &vec![0.0; samples],
            4,
        );
        let alpha = alphas(&estimates);
        let worst = alpha.iter().copied().fold(0.0_f64, f64::max);
        let mean = alpha.iter().sum::<f64>() / alpha.len() as f64;
        eprintln!("clean panel: worst {worst:.4}, mean {mean:.4}");
        assert!(
            worst < 0.03,
            "a clean panel's worst sample came back at {worst}"
        );
    }

    #[test]
    fn the_solver_returns_none_on_a_singular_system() {
        let singular = [1.0, 2.0, 2.0, 4.0];
        assert!(solve(&singular, &[1.0, 2.0], 2).is_none());
        let ordinary = [2.0, 0.0, 0.0, 4.0];
        let answer = solve(&ordinary, &[2.0, 8.0], 2).expect("not singular");
        assert!((answer[0] - 1.0).abs() < 1e-12 && (answer[1] - 2.0).abs() < 1e-12);
    }
}
