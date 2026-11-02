//! **The neighbour-evidence pass** — re-calling each sample's genotypes from its own likelihoods and
//! the samples it shares stretches of chromosome with (`doc/devel/ng/spec/neighbour_evidence_pass.md`).
//!
//! - [`model`] — one sample's walk along a window: the copying model with an autozygous state (§3);
//!   with phasing, the path drawn from it and the haplotypes re-estimated along it (phasing spec §2).
//! - [`random`] — the seeded uniforms phasing draws from (phasing spec §5).
//! - [`records`] — one record as the pass holds it: each sample's likelihoods and allele shares,
//!   and the background template's shares (§3.2, §3.3).
//! - [`templates`] — which other samples a sample copies from in a window (§3.3), and which other
//!   samples' haplotypes in a phasing round (phasing spec §2.3).
//! - [`fit`] — each sample's walk parameters, fitted on the first windows and then frozen (§3.5).
//! - [`windows`] — where each window's core and margins end, under the memory budget (§2.2, §6).
//! - [`phase`] — the phase the rounds estimate, written in blocks of confident phase (phase output
//!   spec).
//! - [`refine`] — re-calling every sample in one window, in phasing's rounds, and the rule for
//!   samples without reads and the fill check (§4.2, §4.3; phasing spec §4).
//! - [`line`](mod@line) — reading a finished VCF line into the pass and writing it back re-called (§4.1).
//! - [`stage`] — the pass as a stage in front of the writer, and whether a run uses it (§2, §5).

pub mod fit;
pub mod line;
pub mod model;
pub mod phase;
pub(crate) mod random;
pub mod records;
pub mod refine;
pub mod stage;
pub mod templates;
pub mod windows;

#[cfg(test)]
pub(crate) mod test_cohort;
