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
//! genotype frequencies is [`mixture_weights`]; the scan over the noise parameters is
//! [`ladder_scan`], which joined it in Milestone D.
//!
//! **What this file holds is the seam and nothing else**: the two traits a path
//! implements to be fitted — [`NoiseModel`] and [`WeightedCell`] — and
//! [`FitTermination`], which belongs to the alternating fit of Milestone E rather than
//! to either half. Each fit's result type sits with the function that builds it:
//! [`ladder_scan::ScanResult`] and [`ladder_scan::FixedFrequencyScanResult`] are the two
//! scans'.

pub mod ladder_scan;
pub mod mixture_weights;

use crate::ng::types::Ploidy;

/// What the scan needs to know about a cell **besides how to score it**: which genotype
/// set it belongs to, and how much of the table it accounts for.
///
/// **A trait rather than the `(Cell, Ploidy, u64)` the architecture sketches**
/// (`arch/parameter_prepass_generic.md` §4.2), for the reason that replaced the same
/// tuple in `generic/histogram.rs`: nothing in a tuple says which of the three members
/// is which, and two of them are integers. Both paths' cell types carry these two facts
/// already, so a tuple would also mean restating them beside a cell that already knows
/// them — one more place for the two to disagree.
pub trait WeightedCell {
    /// The genotype set this cell is scored against. It travels with the cell because
    /// one noise parameter is fitted across every ploidy its reads covered, so a single
    /// scan sees cells of more than one.
    fn ploidy(&self) -> Ploidy;
    /// How much this cell counts for — how many sites landed in it.
    fn sites(&self) -> u64;
}

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
/// at run time in a loop that runs tens of thousands of times per fit — the generic
/// path's ladder is 161 rungs and its cell table holds up to 583 cells, so at most
/// 93,863 evaluations (`arch/parameter_prepass_generic.md` §4.2, which rounds it to
/// "~75,000" against a table that is rarely full).
pub trait NoiseModel {
    /// The cell type this model's histogram is keyed on — a depth and an alternative
    /// count on the SNP/indel path, a table of repeat-length offsets on the STR path.
    ///
    /// **It must also be a [`WeightedCell`]**, which completes the seam in one place: a
    /// path plugs into `fitting/` by supplying a model *and* a cell that knows its ploidy
    /// and its site count. Stated here rather than as a bound on
    /// [`ladder_scan::fit_by_profile_scan`], where a model whose cell knew neither would
    /// compile and fail only at the one call site that scans it.
    type Cell: WeightedCell;
    /// The noise parameters being scanned — error rates on the SNP/indel path, three
    /// slippage parameters on the STR path.
    type NoiseParams;

    /// How many entries [`Self::append_genotype_likelihoods`] appends per cell, and so
    /// how wide the row-major table the scan builds is.
    ///
    /// **Asked of the model rather than derived from the ploidy, because a genotype is
    /// not the same object on the two paths.** On the SNP/indel path a genotype is a
    /// dosage and there are `ploidy + 1` of them. On the STR path it is an unordered
    /// tuple of allele lengths, so a diploid stratum whose alleles span nine lengths has
    /// **45**, not three (`spec/parameter_prepass_ssr.md` §4.2). A scan that computed
    /// the width as `ploidy + 1` would hand
    /// [`mixture_weights::GenotypeLikelihoodTable`] a width the STR path's rows do not
    /// have.
    fn genotypes(&self, ploidy: Ploidy) -> usize;

    /// How likely each genotype makes this cell, at these noise parameters, as natural
    /// logarithms.
    ///
    /// **Appends [`Self::genotypes`] entries and clears nothing**, in the model's own
    /// genotype order — ascending number of alternative copies on the SNP/indel path.
    /// The name says so because the contract is
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
