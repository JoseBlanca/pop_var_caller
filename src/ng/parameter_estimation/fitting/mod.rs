//! The mathematics, with the domain taken out.
//!
//! Given a table of cells — each with a count of how many sites looked like it — and
//! a model saying how likely each genotype makes a cell, this finds the noise
//! parameters and genotype frequencies that best explain the table. It knows nothing
//! about markers, loci or windows; the ladder of candidate noise parameters it steps
//! through is handed to it by the path that owns one.
//!
//! **The one genuine swappable seam in step 4.** The SNP/indel path and the STR path
//! run the same procedure over two different models of what can go wrong with a read:
//! a base miscalled here, a repeat unit gained or lost there
//! (`spec/parameter_prepass.md` §3.2). This path is the first consumer; the STR path
//! is the second, and whether the seam was cut in the right place is a question its
//! plan answers, not this one.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §4. The climb over the
//! genotype frequencies is [`mixture_weights`]; the scan over the noise parameters
//! joins it in Milestone D. The two result types below are Milestone A.

pub mod mixture_weights;

use smallvec::SmallVec;

use crate::ng::types::{LogProb, Ploidy};

/// What a path assumes can go wrong with a read.
///
/// **The one seam between the two paths of step 4**: the same procedure over two models.
/// A base miscalled here, a repeat unit gained or lost there
/// (`spec/parameter_prepass.md` §3.2). The SNP/indel path implements it in
/// `parameter_estimation::generic::noise_model`; the STR path is the second implementor.
///
/// **Static dispatch, deliberately** — the scan is written `<M: NoiseModel>` and not
/// `&dyn NoiseModel`, so the compiler emits one specialised copy per model with
/// `M::Cell` substituted. Sharing the procedure across the two paths then costs nothing
/// at run time in a loop that runs about 75,000 times per fit
/// (`arch/parameter_prepass_generic.md` §4.2).
pub trait NoiseModel {
    /// The cell type this model's histogram is keyed on — a depth and an alternative
    /// count on the SNP/indel path, a table of repeat-length offsets on the STR path.
    type Cell;
    /// The noise parameters being scanned — error rates on the SNP/indel path, three
    /// slippage parameters on the STR path.
    type NoiseParams;

    /// How likely each genotype makes this cell, at these noise parameters, as natural
    /// logarithms.
    ///
    /// **Appends `ploidy + 1` entries and clears nothing**, one per number of
    /// alternative copies, ascending from zero. The name says so because the contract is
    /// what the profile scan is built on: the scan clears one flat buffer per rung and
    /// calls this once per cell, so what comes out is the row-major table
    /// [`mixture_weights::GenotypeLikelihoodTable`] borrows, with no per-cell row and no
    /// copy. A model that cleared instead would silently leave the scan holding one
    /// cell's row.
    ///
    /// `−∞` is a legal entry and says this genotype cannot have produced this cell.
    fn append_genotype_likelihoods(
        &self,
        cell: &Self::Cell,
        noise: &Self::NoiseParams,
        ploidy: Ploidy,
        out: &mut Vec<f64>,
    );
}

/// What one scan over a ladder of noise parameters returned.
///
/// Generic over the noise parameters, because the two paths scan different things: one
/// error rate here, three stutter parameters on the STR path.
#[derive(Clone, PartialEq, Debug)]
pub struct ScanResult<P> {
    /// The winning rung.
    pub noise: P,
    /// The genotype frequencies climbed to at that rung. On the error-rate scan these
    /// are a means rather than an output — the scan is run for `noise` and they are
    /// discarded — while the sample's own rates come from a scan run for these.
    pub frequencies: SmallVec<[f64; 3]>,
    /// What makes "the best-scoring iterate" a defined comparison in an alternating
    /// fit. A [`LogProb`] rather than a bare `f64` because comparing is the only thing
    /// it is for, and `LogProb` carries `ln(0)` as `-∞` — the score of a rung nothing
    /// could have produced — where a linear probability would reach `0` and be
    /// indistinguishable from a value that merely got too small to represent.
    pub log_likelihood: LogProb,
    /// **Whether the answer sat on the ladder's edge, and it is not decoration.** A
    /// read group whose true rate lies outside the ladder — a bad run, heavy
    /// contamination, or any of the arithmetic failures a scan can suffer — has its
    /// answer silently clamped to an endpoint and emitted as though it were fitted,
    /// with a large observation count behind it. This one bit is what stands between a
    /// railed fit and a plausible-looking number.
    pub argmax_at_ladder_end: bool,
}

/// How an alternating fit ended.
///
/// **Emitted rather than discarded**, because a fit that ran out of iterations is still
/// a number a caller would otherwise consume as though it had settled. This is the
/// **outer** alternation between the two tables, which has no convergence proof at all;
/// the inner climb over the genotype frequencies is provably concave, so it cannot get
/// stuck — but it is capped too, because concavity says nothing about how *fast* it
/// arrives and one of `mixture_weights`' own fixtures takes 1,234 passes. The inner cap
/// is not reported to a consumer and this one is: the outer loop can end either way at
/// any iteration count and keeps its best-scoring iterate, so which happened is not
/// derivable from the count.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct FitTermination {
    pub iterations: u32,
    pub converged: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The rail flag is the field with teeth, so it is stated rather than left to a
    /// reader to notice: a scan that railed reports the same shape as one that did not.
    #[test]
    fn a_scan_result_reports_whether_its_answer_sat_on_the_ladders_edge() {
        let railed = ScanResult {
            noise: 0.1_f64,
            frequencies: SmallVec::from_slice(&[0.98, 0.015, 0.005]),
            log_likelihood: LogProb(-1.2e9),
            argmax_at_ladder_end: true,
        };

        assert!(railed.argmax_at_ladder_end);
        assert_eq!(railed.frequencies.len(), 3);
    }

    /// A single-library sample settles after one pass, because the two tables the
    /// alternation reads are then the same table. That is 1,550 of the 1,707 samples in
    /// the archive survey, which is why the coupled fit is low-risk.
    #[test]
    fn a_fit_that_settled_and_one_that_ran_out_are_distinguishable() {
        let settled = FitTermination {
            iterations: 1,
            converged: true,
        };
        let ran_out = FitTermination {
            iterations: 20,
            converged: false,
        };

        assert_ne!(settled, ran_out);
        assert!(settled.converged && !ran_out.converged);
    }
}
