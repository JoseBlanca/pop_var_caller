//! The joint fit — every parameter estimated once, over every sample at the same loci.
//!
//! ng step 4's **second** route to the parameters
//! ([`generic`](super::generic) is the first). Where that route folds each sample's loci
//! into histograms and fits from those, this one keeps raw evidence at a bounded set of
//! loci — **the same loci in every sample** — and fits everything once when every sample
//! is in. What that buys is the locus's own allele frequency in the cohort, a quantity a
//! histogram has forgotten by construction.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_fit.md` (what the route is),
//! `parameter_prepass_joint_loci.md` (which loci), `parameter_prepass_joint_records.md`
//! (what is recorded at each). Types: the three `arch/parameter_prepass_joint_*.md`.
//!
//! **Under construction.** [`loci`] settles which loci every sample keeps, [`records`] what
//! each writes down there, [`coverage`] fills the window summary they travel with, and
//! [`fit`] is the estimator — **the ordinary-position half of it**. The repeat-tract half,
//! contamination and the duplicated-stretch class are not in it yet, and its own docs say so.

pub mod coverage;
pub mod fit;
pub mod loci;
pub mod records;
