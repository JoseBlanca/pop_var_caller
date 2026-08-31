//! ng's calling run — the machinery the two variant callers drive.
//!
//! A calling run reads every sample's locus observations in coordinate order and emits called
//! variants in coordinate order. `doc/devel/ng/spec/run_streaming.md` owns that outer shape
//! (the caller objects, the sources, the VCF writing) and
//! `doc/devel/ng/arch/run_streaming.md` the types; this module is where its parts land as
//! they are built.
//!
//! **Landed so far:** [`cohort_merge`], which turns k samples' observations into one stream of
//! cohort observations; [`segments`]'s [`Segmentation`], the ground every sample of a run
//! walks; and [`callers`]'s [`AlignedFilesVariantCaller`], constructed but not yet iterating.

pub mod callers;
pub mod cohort_merge;
pub mod segments;

pub use callers::{AlignedFilesVariantCaller, AlignmentInputs, MergeParameters};
pub use segments::{Segmentation, SegmentationInputs};

use std::path::PathBuf;

use crate::ng::read::input::IngestError;
use crate::ng::repeat_catalog::RepeatCatalogError;

/// What can go wrong driving a run.
///
/// **Every variant names what a person can act on**, and each one names the *sample* or the
/// *file* the trouble is in — a run over thousands of samples that says only "it failed"
/// leaves nobody anywhere to look (spec §9).
///
/// **The reason comes from the cause, not from the top line**, so a command reporting one of
/// these must render the whole chain with
/// [`format_error_chain`](crate::error_render::format_error_chain), never `Display` alone. A
/// bare `Display` says which sample would not open; the chain says its index is missing and
/// where it was looked for.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum RunError {
    /// The run's segments could not be read out of the repeat catalog.
    ///
    /// **The path is carried here because the catalog's own error often does not have it**:
    /// a digest mismatch, over-permissive criteria, differing scan weights and a differing
    /// tool version all describe the file without naming it, and a person with several
    /// catalogs on disk cannot act on that.
    #[error("the run's segments could not be read from the repeat catalog {}", path.display())]
    Catalog {
        /// The catalog the run read.
        path: PathBuf,
        /// What reading it hit.
        #[source]
        source: RepeatCatalogError,
    },

    /// One sample's alignment files could not be opened.
    ///
    /// **The sample's name is a field of its own because the wrapped error does not carry
    /// it**: an open failure knows which file it was, not which individual the file holds.
    /// The cause is boxed to keep this type small, since it travels in every `Result` a run
    /// returns.
    #[error("sample {sample}: its alignment files could not be opened")]
    OpeningSample {
        /// The sample as its read groups name it.
        sample: String,
        /// What opening it hit — the file, and why.
        source: Box<IngestError>,
    },
}
