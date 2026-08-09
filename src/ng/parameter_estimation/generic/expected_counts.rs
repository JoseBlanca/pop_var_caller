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
