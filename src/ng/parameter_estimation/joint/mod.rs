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
//! **Under construction.** [`loci`] settles which loci every sample keeps, [`census`] what
//! each writes down there, [`fit`] is
//! the estimator over ordinary positions — including the class for positions a sample carries
//! more copies of than the reference does — [`ssr_fit`] the estimator over repeat tracts, and
//! [`contamination`] the share of a sample's reads that came from another individual, which no
//! per-sample route can produce at all. Each module's own docs say what it does and does not do.
//!
//! [`slippage_curve`] holds how the repeat-tract slippage *level* rises with repeat count, which
//! is what lets a thin cell read a line drawn through all its neighbours instead of copying one.
//! [`share_curve`] does the same for the other two slippage numbers — which way a read slipped
//! and how far — whose shapes differ from the level's and from each other.
//!
//! [`sequencing_batches`] holds who was sequenced beside whom, as the run was told — the
//! grouping the contaminating population is drawn from, and the one thing in this module that is
//! declared rather than fitted.
//!
//! **The two halves of the estimator run in that order and depend on each other in one direction
//! only**: [`ssr_fit`] takes each sample's homozygote excess from [`fit`] and gives nothing back,
//! so a run may drop the ordinary-position records before it reads a single tract.

pub mod census;
pub mod census_file;
pub mod contamination;
pub mod fit;
pub mod loci;
pub mod sequencing_batches;
pub mod share_curve;
pub mod slippage_curve;
pub mod ssr_fit;
pub mod stratum_fits;
