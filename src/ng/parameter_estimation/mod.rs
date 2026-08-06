//! ng step 4 — the parameters the caller runs on, measured from the sample's own
//! loci before anything is called.
//!
//! Four numbers per sample come out of the SNP/indel path: a per-read-group error
//! rate, the sample's heterozygosity, its homozygous-non-reference rate, and its
//! inbreeding coefficient. They are measured from **every** covered position,
//! including the overwhelming majority that show no alternative allele at all —
//! which is what separates this step from production's estimator. Production writes
//! the pure-reference columns; it is production's *heterozygosity accumulator* that
//! never looks at them (`spec/parameter_prepass.md` §2.1), so the loss is in the
//! estimator rather than in the data — and what is lost is the strongest evidence
//! there is about the error rate.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_generic.md` (the design and its
//! rationale), `doc/devel/ng/spec/parameter_prepass.md` (the shared framing), and
//! `doc/devel/ng/arch/parameter_prepass_generic.md` (types and interfaces).
//!
//! Two sub-units, split so that the shaping of data and the mathematics on it never
//! live in one file:
//!
//! - [`fitting`] — the mathematics. Knows nothing about markers, loci or windows: it
//!   is given a table of numbers and returns the values that best explain them. A
//!   folder rather than a file because it is the one genuine swappable seam — one
//!   trait, an implementation on this path and a second on the STR path.
//! - [`generic`] — the SNP/indel path: the two accumulators, the cell table, the
//!   vocabulary they are keyed on, and what each of the four numbers is fitted from.
//!
//! A third sub-unit for the STR path joins them later, which is why this file holds
//! no vocabulary of its own: an error-rate ladder in per-base probabilities and a
//! window size for runs of homozygosity are the SNP/indel path's, and they live in
//! [`generic`] where the STR path will not inherit them.

pub mod fitting;
pub mod generic;
