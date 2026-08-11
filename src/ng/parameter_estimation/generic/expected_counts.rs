//! Test support: cell tables holding **the counts an infinite genome would produce**.
//!
//! The device the two research harnesses use (research note §1): replace each cell's
//! observed count with its probability under a known truth, and the answer a fit returns
//! has no sampling noise in it — so a departure from the generating parameters is bias,
//! decided rather than estimated.
//!
//! **Shared by every fit's tests rather than restated in each**, because two copies of a
//! fixture generator are two chances for one of them to drift into agreeing with the code
//! it is supposed to check.
//!
//! Not compiled outside tests.

use std::sync::Arc;

use super::depth_bins::DepthBinEdges;
use super::histogram::{Cell, DepthAltHistogram, DepthAndAltReads, SiteKey};
use crate::genetics::lgamma;
use crate::ng::types::{Bp, Ploidy};

/// `p_j(ε)` — the chance one read shows something other than the reference base, at `j`
/// alternative copies of `ploidy`.
///
/// **The same expression `noise_model` scores with, restated here rather than reached
/// for**, so that the fixture and the model are two statements of the rule and a fit that
/// recovers its own parameters is a claim about the **fit** rather than about the
/// expression. Whether the expression itself is right is Milestone D2's four identities.
///
/// **What that division of labour rests on, stated because it is easy to over-read.** On a
/// table of expected counts the score is `N · Σ_c p_c(θ₀) · ln p_c(θ)`, which Gibbs'
/// inequality maximises at `θ = θ₀` for *any* rule whose cell probabilities sum to one over
/// the cell space. So recovery of the generating parameters cannot catch an expression that
/// is wrong self-consistently — and it is not vacuous either: it fails the moment a fit
/// mislabels, misgathers, misreports or misranks. The sum-to-one it leans on holds only
/// because every site in a table sits at one exact depth, which is why
/// [`table_generated_at`] takes a single `depth`.
pub(super) fn alternative_read_probability(alt_copies: u8, ploidy: Ploidy, error_rate: f64) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    carried * (1.0 - error_rate / 3.0) + (1.0 - carried) * error_rate
}

/// A table of `sites` sites, every one at `depth`: cell `k` gets
/// `sites · Σ_j π_j · Binomial(k; depth, p_j(ε))`, rounded.
///
/// One depth, so the cell's mean depth is that depth exactly, the binning rule contributes
/// nothing to the answer, and the cell probabilities sum to one — which is what
/// [`alternative_read_probability`]'s note says the recovery claim needs.
///
/// Each cell is rounded on its own, so the table holds a site or two either side of
/// `sites`. Nothing should assert the total against `sites`; compare against
/// [`DepthAltHistogram::total_loci`] of the table itself.
///
/// # Panics
///
/// If `depth` is above the ladder's cap, or if `genotype_frequencies` is not one entry per
/// dosage of `ploidy`.
pub(super) fn table_generated_at(
    edges: &Arc<DepthBinEdges>,
    depth: u32,
    error_rate: f64,
    ploidy: Ploidy,
    genotype_frequencies: &[f64],
    sites: f64,
) -> DepthAltHistogram<u64> {
    assert_eq!(
        genotype_frequencies.len(),
        usize::from(ploidy.get()) + 1,
        "a ploidy-{ploidy} truth needs one frequency per dosage"
    );

    let mut histogram = DepthAltHistogram::new(Arc::clone(edges));
    for alt_reads in 0..=depth {
        let mut probability = 0.0;
        for (alt_copies, &frequency) in genotype_frequencies.iter().enumerate() {
            let p = alternative_read_probability(alt_copies as u8, ploidy, error_rate);
            probability += frequency * binomial_probability(alt_reads, depth, p);
        }
        let in_cell = (sites * probability).round() as u64;
        for _ in 0..in_cell {
            histogram.add_site(DepthAndAltReads::new(depth, alt_reads), Bp(1));
        }
    }
    histogram
}

/// `Binomial(k; n, p)`, computed through logs.
///
/// **An earlier version built the term up from `k = 0` by repeated multiplication**, starting
/// at `(1 − p)^n`, on the stated grounds that it "keeps every intermediate inside `f64` at the
/// depths used here". It did — until F2 raised the depth. At `n = 124` and a homozygous
/// non-reference genotype, `p` is `1 − ε/3` and the starting term is `0.00033^124`, which is
/// about `10⁻⁴³¹` and **underflows to exactly zero**; every later multiplication keeps it
/// there. The whole genotype then contributes nothing to any cell, so the table silently
/// omits its homozygous non-reference sites and the fit correctly reports a frequency of
/// 0.0000 for a class the fixture never generated.
///
/// That is the failure mode this file exists to avoid, in the file every fit's tests share: a
/// fixture that quietly asks a different question than the one written down. Logs have no
/// such cliff — `ln` of the same term is −992, an ordinary number — and `exp` at the end
/// returns a genuinely tiny probability as a tiny probability rather than as nothing.
fn binomial_probability(successes: u32, trials: u32, p: f64) -> f64 {
    if p <= 0.0 {
        return f64::from(u8::from(successes == 0));
    }
    if p >= 1.0 {
        return f64::from(u8::from(successes == trials));
    }
    let ln_choose = lgamma(f64::from(trials) + 1.0)
        - lgamma(f64::from(successes) + 1.0)
        - lgamma(f64::from(trials - successes) + 1.0);
    (ln_choose + f64::from(successes) * p.ln() + f64::from(trials - successes) * (1.0 - p).ln())
        .exp()
}

/// The depth distribution of a real 30× human walk: `(depth, loci)` for the 551,843 generic
/// loci ng produced over HG002's 100 GIAB confident regions.
///
/// **Measured, not invented, and kept here because it cannot be re-derived from anything in
/// this repository.** It came from walking
/// `benchmarks/giab/per_sample/bam/30x/HG002.30x.seed42.bam` — an alignment no worktree
/// carries — under `ReadFilterConfig::default()`, and counting the depth
/// `count_whole_site` gave each locus. 69 distinct depths, mean 30.12.
///
/// **What it is for.** A fixture at one depth cannot see anything the depth distribution
/// decides, and the second class of site is exactly such a thing: the measured excess
/// heterozygosity it exists to remove falls monotonically with depth, from 27.7 sites per
/// 10,000 loci at depth 11–15 to 1.0 above depth 40
/// (`research/noise_model_overdispersion_2026-08-10.md`). A 30× sample is not 30 reads
/// everywhere; it is a long shoulder reaching down to 11, and the shoulder is where the
/// two explanations of an alternative read stop being distinguishable.
pub(super) const REAL_DEPTH_DISTRIBUTION: &[(u32, u64)] = &[
    (1, 72),
    (2, 192),
    (3, 379),
    (4, 340),
    (5, 380),
    (6, 343),
    (7, 372),
    (8, 430),
    (9, 553),
    (10, 694),
    (11, 812),
    (12, 1217),
    (13, 1425),
    (14, 1765),
    (15, 2773),
    (16, 3809),
    (17, 5554),
    (18, 7953),
    (19, 10349),
    (20, 12739),
    (21, 14989),
    (22, 16918),
    (23, 19382),
    (24, 21944),
    (25, 23339),
    (26, 25612),
    (27, 27464),
    (28, 29740),
    (29, 30112),
    (30, 30988),
    (31, 30321),
    (32, 28695),
    (33, 26217),
    (34, 24070),
    (35, 22863),
    (36, 20201),
    (37, 18248),
    (38, 15892),
    (39, 13822),
    (40, 11422),
    (41, 9444),
    (42, 7974),
    (43, 6628),
    (44, 5375),
    (45, 4295),
    (46, 3479),
    (47, 2727),
    (48, 2168),
    (49, 1572),
    (50, 1032),
    (51, 812),
    (52, 527),
    (53, 353),
    (54, 327),
    (55, 236),
    (56, 139),
    (57, 91),
    (58, 64),
    (59, 51),
    (60, 43),
    (61, 39),
    (62, 34),
    (63, 9),
    (64, 11),
    (65, 9),
    (66, 4),
    (67, 6),
    (68, 3),
    (69, 1),
];

/// The cells an **infinite** genome of this depth distribution would produce: one cell per
/// `(depth, alternative reads)` pair, holding that pair's exact expected number of sites.
///
/// **One cell per depth rather than per depth *bin*, which is deliberate and is not what the
/// accumulator does.** Sharing a bin makes several depths one cell scored at their mean, and
/// what that costs is the binning bias the ladder was chosen to bound — 0.054 rungs and 0.3%
/// (research note 2026-08-06 §4.3). Keeping the depths apart takes that term out of the
/// answer, so a departure here is the *fit's* bias and nothing else. It is the setting the
/// two controls of `research/noise_model_overdispersion_2026-08-10.md` were measured in.
///
/// `site_noise` gives the sites two classes; `None` generates a world with one error rate,
/// which is the control a richer model must leave alone.
///
/// # Panics
///
/// If `genotype_frequencies` is not one entry per dosage of `ploidy`.
pub(super) fn cells_over_a_real_depth_distribution(
    edges: &Arc<DepthBinEdges>,
    error_rate: f64,
    site_noise: Option<(f64, f64)>,
    ploidy: Ploidy,
    genotype_frequencies: &[f64],
) -> Vec<Cell> {
    assert_eq!(
        genotype_frequencies.len(),
        usize::from(ploidy.get()) + 1,
        "a ploidy-{ploidy} truth needs one frequency per dosage"
    );
    let (noisy_fraction, noisy_rate) = site_noise.unwrap_or((0.0, error_rate));

    let mut cells = Vec::new();
    for &(depth, loci) in REAL_DEPTH_DISTRIBUTION {
        for alt_reads in 0..=depth {
            let mut probability = 0.0;
            for (alt_copies, &frequency) in genotype_frequencies.iter().enumerate() {
                let clean = alternative_read_probability(alt_copies as u8, ploidy, error_rate);
                let noisy = alternative_read_probability(alt_copies as u8, ploidy, noisy_rate);
                probability += frequency
                    * ((1.0 - noisy_fraction) * binomial_probability(alt_reads, depth, clean)
                        + noisy_fraction * binomial_probability(alt_reads, depth, noisy));
            }
            let sites = probability * loci as f64;
            // Below this a cell holds less than a millionth of a site and contributes
            // nothing any fit can resolve, while the cell count trebles.
            if sites > 1e-6 {
                cells.push(Cell {
                    key: SiteKey::pooled(edges.bin_for(depth), alt_reads),
                    ploidy,
                    // The weight is a count of sites and the type is integral, so the
                    // fractional part is lost. Scaling by a million first keeps five
                    // significant figures of the rarest cell that matters, and every fit
                    // here is invariant to a common factor on the weights.
                    sites: (sites * 1e6).round() as u64,
                    mean_depth: f64::from(depth),
                });
            }
        }
    }
    cells
}

/// The same world as [`cells_over_a_real_depth_distribution`], as a **table** rather than as
/// cells — which is what the whole coupled fit takes, where `fit_site_noise` alone takes
/// cells.
///
/// The difference matters for exactly one question: `fit_site_noise` is handed the clean rate,
/// so a test built on cells asks whether the *second* class is recovered given the first,
/// while the coupled fit is handed nothing and has to find both. On real alignments the clean
/// rate came back on the same rung the one-class fit chose, twice, which is a claim only a
/// table-shaped world can check.
///
/// **A count of sites and not a scaled weight**, so cells below half a site round away
/// entirely — a table is what an accumulator built and holds whole loci. That is the one way
/// this world is coarser than its cells-shaped sibling, and it is why the site total is the
/// depth distribution's own 551,843 rather than a millionfold weight.
///
/// # Panics
///
/// If `genotype_frequencies` is not one entry per dosage of `ploidy`.
pub(super) fn table_over_a_real_depth_distribution(
    edges: &Arc<DepthBinEdges>,
    error_rate: f64,
    site_noise: Option<(f64, f64)>,
    ploidy: Ploidy,
    genotype_frequencies: &[f64],
) -> DepthAltHistogram<u64> {
    assert_eq!(
        genotype_frequencies.len(),
        usize::from(ploidy.get()) + 1,
        "a ploidy-{ploidy} truth needs one frequency per dosage"
    );
    let (noisy_fraction, noisy_rate) = site_noise.unwrap_or((0.0, error_rate));

    let mut histogram = DepthAltHistogram::new(Arc::clone(edges));
    for &(depth, loci) in REAL_DEPTH_DISTRIBUTION {
        for alt_reads in 0..=depth {
            let mut probability = 0.0;
            for (alt_copies, &frequency) in genotype_frequencies.iter().enumerate() {
                let clean = alternative_read_probability(alt_copies as u8, ploidy, error_rate);
                let noisy = alternative_read_probability(alt_copies as u8, ploidy, noisy_rate);
                probability += frequency
                    * ((1.0 - noisy_fraction) * binomial_probability(alt_reads, depth, clean)
                        + noisy_fraction * binomial_probability(alt_reads, depth, noisy));
            }
            let in_cell = (probability * loci as f64).round() as u64;
            for _ in 0..in_cell {
                histogram.add_site(DepthAndAltReads::new(depth, alt_reads), Bp(1));
            }
        }
    }
    histogram
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::generic::histogram::Attribution;

    fn edges() -> Arc<DepthBinEdges> {
        Arc::new(DepthBinEdges::new())
    }

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a positive copy number")
    }

    /// **The fixture's own mass check, and it exists because the first version of this
    /// experiment did not have one.** Written in Python while the model was being chosen,
    /// it enumerated alternative counts only to 25 — which silently discarded almost all of
    /// the homozygous-non-reference probability, since such a site shows alternative reads
    /// at nearly its whole depth. Every model then came back with that class 73% low, and
    /// the numbers looked like a finding about the models rather than a bug in the
    /// generator.
    ///
    /// Each depth's cells must carry that depth's whole complement of loci.
    ///
    /// **The tolerance is 1 in 10,000 and the worst real loss is 5 in a million**, at depth
    /// 69, where the distribution holds a single locus and the builder's prune drops cells
    /// worth a fraction of a micro-site.
    ///
    /// **What it must catch is smaller than it sounds, and the first version of this note
    /// got it wrong.** Truncating at 25 alternative reads destroys 73% of the
    /// *homozygous-non-reference class* — that is the number the Python run reported — but
    /// only **0.1% of the total probability mass**, because that class is a thousandth of
    /// the sites. Measured against this fixture: 9.9 × 10⁻⁴ missing at depth 26, ten times
    /// the tolerance rather than the ten thousand times "73%" would imply. The check that
    /// sees the truncation loudly is the next one, which asks *where* the mass went.
    #[test]
    fn every_depths_cells_hold_that_depths_whole_complement_of_loci() {
        for site_noise in [None, Some((0.0088, 5.29e-2))] {
            let cells = cells_over_a_real_depth_distribution(
                &edges(),
                1.0e-3,
                site_noise,
                diploid(),
                &[0.9885, 0.0105, 0.0010],
            );
            for &(depth, loci) in REAL_DEPTH_DISTRIBUTION {
                let held: u64 = cells
                    .iter()
                    .filter(|cell| (cell.mean_depth - f64::from(depth)).abs() < 1e-9)
                    .map(|cell| cell.sites)
                    .sum();
                let want = loci * 1_000_000;
                let lost = want.abs_diff(held) as f64 / want as f64;
                assert!(
                    lost < 1e-4,
                    "depth {depth}: the cells hold {held} against {want} — {:.3e} of the \
                     probability is missing, so a fit over this world is answering a \
                     question about a truncated cell space",
                    lost
                );
            }
        }
    }

    /// The homozygous-non-reference class is present and lands where it belongs — near the
    /// top of each depth's range, not near the bottom. This is the class the Python
    /// truncation destroyed, and a mass check alone would not have said *where* the mass
    /// went.
    #[test]
    fn the_homozygous_non_reference_class_lands_at_nearly_full_alternative_counts() {
        let cells = cells_over_a_real_depth_distribution(
            &edges(),
            1.0e-3,
            None,
            diploid(),
            &[0.9885, 0.0105, 0.0010],
        );
        let deep: Vec<&Cell> = cells
            .iter()
            .filter(|cell| (cell.mean_depth - 30.0).abs() < 1e-9)
            .collect();
        let at_full: u64 = deep
            .iter()
            .filter(|cell| cell.key.alt_reads() >= 28)
            .map(|cell| cell.sites)
            .sum();
        let expected = (0.0010 * 30_988.0 * 1e6) as u64;
        assert!(
            at_full > expected / 2,
            "at depth 30 the cells with 28 or more alternative reads hold {at_full}, and a \
             frequency of 0.0010 over 30,988 loci should put about {expected} there"
        );
    }

    /// A world given no second class of site is generated by the one-rate expression, so
    /// that the control a richer model must leave alone really is the world today's model
    /// describes.
    #[test]
    fn no_site_noise_generates_the_same_world_as_a_zero_share() {
        let plain = cells_over_a_real_depth_distribution(
            &edges(),
            1e-3,
            None,
            diploid(),
            &[0.98, 0.015, 0.005],
        );
        let zero_share = cells_over_a_real_depth_distribution(
            &edges(),
            1e-3,
            Some((0.0, 9.9e-2)),
            diploid(),
            &[0.98, 0.015, 0.005],
        );
        assert_eq!(plain.len(), zero_share.len());
        for (a, b) in plain.iter().zip(&zero_share) {
            assert_eq!(a.sites, b.sites);
            assert!((a.mean_depth - b.mean_depth).abs() < 1e-12);
        }
    }

    /// Every cell is pooled, because a world with one library has nothing to attribute —
    /// stated by a test so that a later multi-library oracle cannot quietly reuse this one.
    #[test]
    fn the_generated_cells_are_all_pooled() {
        let cells = cells_over_a_real_depth_distribution(
            &edges(),
            1e-3,
            None,
            diploid(),
            &[0.98, 0.015, 0.005],
        );
        assert!(!cells.is_empty());
        assert!(
            cells
                .iter()
                .all(|cell| matches!(cell.key.attribution(), Attribution::Pooled)),
            "a single-library world attributes nothing"
        );
    }
}
