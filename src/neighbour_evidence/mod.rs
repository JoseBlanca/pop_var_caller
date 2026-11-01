//! **The neighbour-evidence pass** — re-calling each sample's genotypes from its own likelihoods and
//! the samples it shares stretches of chromosome with (`doc/devel/ng/spec/neighbour_evidence_pass.md`).
//!
//! - [`model`] — one sample's walk along a window: the copying model with an autozygous state (§3).

pub mod model;
