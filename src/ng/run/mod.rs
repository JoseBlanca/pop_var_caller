//! ng's calling run — the machinery the two variant callers drive.
//!
//! A calling run reads every sample's locus observations in coordinate order and
//! emits called variants in coordinate order. `doc/devel/ng/spec/run_streaming.md`
//! owns that outer shape (the caller objects, the sources, the look-ahead and the
//! VCF writing); this module is where its parts land as they are built.
//!
//! Landed so far: [`cohort_merge`]'s three run parameters. The stage itself — which
//! turns k samples' observations into one stream of cohort observations, in parallel —
//! lands beside them as the plan's milestones complete.

pub mod cohort_merge;
