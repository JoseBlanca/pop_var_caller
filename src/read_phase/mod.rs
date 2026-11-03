//! **Phase from reads** — which allele each of a sample's reads showed at each written record, kept
//! from the cohort merge to whichever step writes phase
//! (`doc/devel/ng/spec/neighbour_evidence_read_phase.md`).
//!
//! A read, or a read pair, comes from one molecule, so two heterozygous sites it covers sit on one
//! chromosome copy as the read shows them. This module holds what that evidence is made of, from
//! the moment the merge still knows each read by its chain id.
//!
//! - [`capture`] — one sample's read alleles at one cohort locus, taken in the merge before the
//!   chain ids are dropped, and settled: one per read, at most 32 a sample (§2.2, §2.4).

pub mod capture;
