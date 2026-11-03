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
//! - [`record_reads`] — every sample's read alleles at one written record, numbered as the
//!   record numbers its alleles (§2.3).
//! - [`keeper`] — on the ordered thread, the read alleles of reads seen at two records or more,
//!   numbered as per-sample fragments; the rest dropped (§2.4).
//! - [`link`] — what the fragments say about one link between two heterozygous records: its
//!   log-odds, and the thresholds both phase writers judge it by (§2.5).
//! - [`stage`] — the writer's stage that phases from reads alone where refinement does not run
//!   (§2.6).

pub mod capture;
pub mod keeper;
pub mod link;
pub mod record_reads;
pub mod stage;
