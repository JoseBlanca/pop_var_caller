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
use super::histogram::{DepthAltHistogram, DepthAndAltReads};
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
            // The binomial term, built up from `k = 0` rather than through a factorial,
            // which at the depths used here keeps every intermediate inside `f64`.
            let mut term = (1.0 - p).powi(depth as i32);
            for step in 1..=alt_reads {
                term *= f64::from(depth - step + 1) / f64::from(step) * p / (1.0 - p);
            }
            probability += frequency * term;
        }
        let in_cell = (sites * probability).round() as u64;
        for _ in 0..in_cell {
            histogram.add_site(DepthAndAltReads::new(depth, alt_reads), Bp(1));
        }
    }
    histogram
}
