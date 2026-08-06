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
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §4. Implemented in
//! Milestone D.

pub mod mixture_weights;
