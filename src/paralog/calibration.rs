//! **What pass two settles for pass three: how common hidden duplications are in this run,
//! and where to cut.**
//!
//! The per-locus score ([`super::score_locus_for_paralogy`]) says which of the two stories
//! fits a locus better; it does not say how *likely* that locus is to be a hidden
//! duplication. That needs the rate at which they occur in this run, which is estimated from
//! the run's own scores by the EM in [`super::prior`]. With it, each locus's likelihood ratio
//! becomes a probability and a false-discovery q-value, and the operator's target FDR
//! resolves to a cut.
//!
//! **Copied from production's `src/var_calling/paralog_filter/calibrate.rs`, from its first item
//! down to but not including the cohort inbreeding coefficient** — the one place in this module
//! where ng takes a span of a larger file rather than a whole one, because production keeps this
//! type in the same file as its spill-streaming driver, which ng does not want. Everything below
//! the `use` is that span, line for line and byte for byte but for **seven** lines a textual guard
//! (`copy_fidelity.rs`, deleted at promotion step C3) declared and checked: one path into
//! production, one rustdoc link to the driver ng leaves behind, and five items whose `pub(crate)`
//! is widened to `pub` because ng re-exports them from a public module and a `pub use` cannot
//! re-export a crate-private item. The file was guarded there like every other copy. Since
//! 2026-09-14 its `exp` and `ln` also go through [`crate::float`], so the posterior is the same on
//! macOS and Linux.
//!
//! **ng's `F` is per sample where production's is one cohort number** — the parameters file
//! carries a fitted coefficient for each sample (`doc/devel/ng/spec/hidden_paralog_filter.md`
//! §3.2), and the copied scorer already takes it as a per-sample slice. So production's
//! `cohort_inbreeding` helper, which sits immediately after this span in its file and
//! fills a slice with one coefficient repeated for every sample, is deliberately **not**
//! copied: it is the piece ng replaces rather than ports.

use super::{EmConfig, ParalogFdrCurve, ParalogPrior};
use crate::float;

/// Default fallback `π` used when the EM does not converge — a rare-paralog
/// prior deliberately low so a pathological dataset yields a conservative
/// (few-drops) calibration rather than trusting a junk estimate. Chosen to
/// match [`crate::paralog::prior::DEFAULT_EM_START`] (`0.03`); they are
/// independent constants, so this comment is the only thing coupling them.
/// The operator is warned when this is used (S6).
pub const DEFAULT_FALLBACK_PARALOG_PRIOR: f64 = 0.03;

/// Tuning knobs for the calibration. [`Default`] is the production configuration.
///
/// There is deliberately **no minimum-sample gate**: an under-powered locus
/// self-gates in the likelihood ratio (when the coverage can't tell 1× from 2×,
/// H1 and H2 fit the data equally → LR≈0 → kept), so a count threshold would be
/// redundant. See `doc/devel/architecture/hidden_paralog_single_sample_scoring.md`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CalibrationConfig {
    /// EM configuration for the prior estimate.
    pub em: EmConfig,
    /// Fallback `π` when the EM does not converge.
    pub fallback_prior: f64,
}

impl Default for CalibrationConfig {
    fn default() -> Self {
        Self {
            em: EmConfig::default(),
            fallback_prior: DEFAULT_FALLBACK_PARALOG_PRIOR,
        }
    }
}

/// The global calibration derived from the spill: the prior, the FDR curve, and
/// the resolved LR cut for the target FDR.
#[derive(Debug, Clone)]
pub struct ParalogCalibration {
    /// The empirical-Bayes prior `π` (from the EM, or the fallback). Its
    /// `converged` flag is `false` when the fallback was used.
    pub prior: ParalogPrior,
    /// The tail-FDR curve `q_of_lr`, built over the same histogram bins the LRs
    /// were folded into.
    pub curve: ParalogFdrCurve,
    /// The least-stringent LR at/above which a locus's tail FDR is `<=
    /// target_fdr`; `None` if the target is unachievable (nothing is flagged).
    /// Recorded for VCF-header provenance (S5).
    pub lr_threshold: Option<f64>,
    /// The operator's target false-discovery rate — the operating knob.
    pub target_fdr: f64,
}

impl ParalogCalibration {
    /// Whether a locus with likelihood ratio `lr` is flagged (dropped): its
    /// tail-FDR q-value is at or below the target. Uses the curve directly, so
    /// an unachievable target flags nothing (q never reaches it). Equivalent to
    /// `lr >= lr_threshold` by the curve's monotonicity, but stated on the FDR
    /// so the guarantee is explicit.
    ///
    /// A non-finite `lr` is never flagged: the histogram itself refuses to fold
    /// a non-finite LR (it is not a valid Bayes factor), so it contributed to
    /// neither `π` nor the curve — classifying on it here would be treating the
    /// same value by a different rule than the calibration did. Screened
    /// explicitly because `q_of_lr` would otherwise saturate a `NaN` into bin 0.
    pub fn flags(&self, lr: f64) -> bool {
        lr.is_finite() && self.curve.q_of_lr(lr) <= self.target_fdr
    }

    /// The empirical-Bayes posterior `P(paralog | data) = σ(LR + logit π)` for a
    /// locus's likelihood ratio, using the calibration's estimated prior `π`.
    /// This is the same posterior the prototype reports; the caller emits it as
    /// the `PARALOG_POST` INFO field.
    ///
    /// `None` for a non-finite `lr` (an unscored locus never entered the
    /// calibration), or a degenerate prior (`π ∉ (0, 1)`) where the log-odds
    /// offset is undefined — the field is then omitted rather than emitting a
    /// saturated `0`/`1`.
    pub fn posterior(&self, lr: f64) -> Option<f64> {
        let pi = self.prior.prior_probability;
        if !lr.is_finite() || pi <= 0.0 || pi >= 1.0 {
            return None;
        }
        let log_odds = lr + float::ln(pi / (1.0 - pi)); // LR + logit(π)
        Some(1.0 / (1.0 + float::exp(-log_odds))) // σ(·)
    }
}
