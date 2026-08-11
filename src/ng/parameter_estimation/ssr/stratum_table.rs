//! One stratum's evidence: every distinct locus shape and how many loci had it.
//!
//! **An entry is a locus, not a read**, and that is the design's central choice rather than
//! a storage detail. A read carries no genotype — it drew one of its locus's alleles and
//! then slipped — so a tally that pools reads across loci holds the allele spectrum
//! convolved with the slippage kernel, and recovering the kernel from that means undoing a
//! convolution with both halves unknown. Measured, the fitted slippage level then moves
//! **333-fold depending only on where the search starts**; keyed by locus the same fit is
//! exactly unbiased (`spec/parameter_prepass_ssr.md` §4.1).
//!
//! What separates two loci in this table is only their shape — how their reads fell across
//! the offset buckets — so loci that looked alike are counted together and which loci they
//! were is never asked again.
//!
//! Empty until Milestone B; the buckets it will count are in [`super`].
