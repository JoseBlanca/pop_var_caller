//! **The neighbour-evidence pass** — re-calling each sample's genotypes from its own likelihoods and
//! the samples it shares stretches of chromosome with (`doc/devel/ng/spec/neighbour_evidence_pass.md`).
//!
//! - [`model`] — one sample's walk along a window: the copying model with an autozygous state (§3).
//! - [`records`] — one record as the pass holds it: each sample's likelihoods and allele shares,
//!   and the background template's shares (§3.2, §3.3).
//! - [`templates`] — which other samples a sample copies from in a window (§3.3).
//! - [`fit`] — each sample's walk parameters, fitted on the first windows and then frozen (§3.5).
//! - [`windows`] — where each window's core and margins end, under the memory budget (§2.2, §6).
//! - [`refine`] — re-calling every sample in one window, and the rule for samples without reads
//!   (§4.2).
//! - [`line`] — reading a finished VCF line into the pass and writing it back re-called (§4.1).
//! - [`stage`] — the pass as a stage in front of the writer, and whether a run uses it (§2, §5).

pub mod fit;
pub mod line;
pub mod model;
pub mod records;
pub mod refine;
pub mod stage;
pub mod templates;
pub mod windows;

#[cfg(test)]
mod test_cohort;
