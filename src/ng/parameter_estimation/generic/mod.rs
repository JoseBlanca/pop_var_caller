//! The SNP/indel path: two tallies of what a sample's sites looked like, and the
//! four numbers fitted from them.
//!
//! Two accumulators, differing only in how a site is keyed. The **read-group** one
//! enters a site once per read group that covered it, because an error rate describes
//! the chemistry and two libraries of one sample can genuinely differ. The
//! **windowed** one enters that same site once, at its total depth, because
//! heterozygosity describes the individual — one genome has one heterozygosity
//! however many libraries were used to read it. Neither is derivable from the other
//! once a sample has two read groups
//! (`arch/parameter_prepass_generic.md` §3).
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_generic.md` and its architecture
//! companion. Built across Milestones B, C, E and F.

pub mod depth_and_alt_reads;
pub mod histogram;
pub mod runs;
