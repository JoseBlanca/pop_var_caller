//! One STR locus, reduced to the entry a stratum's table is keyed on.
//!
//! This is the only place that decides what a read's repeat count *is*, which the locus
//! type cannot answer: a locus carries its reference bases, its motif and the sequence each
//! read showed, and turning those into "this read sits two copies short of the reference"
//! is a division that has to be done somewhere and done once.
//!
//! Its own file because the shaping of data and the mathematics on it never live together
//! (`arch/parameter_prepass_ssr.md` §Module home) — [`super::slippage`] is the mathematics.
//!
//! Empty until Milestone C; the vocabulary it will produce is in [`super`].
