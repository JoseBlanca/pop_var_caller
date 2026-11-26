//! Errors raised while finding or building an alignment file's index.

use std::path::PathBuf;

use thiserror::Error;

/// Errors raised by [`crate::bam::index_preflight::preflight_alignment_indexes`].
///
/// The variants here are deliberately CLI-vocabulary-free (no `--flag`
/// names). The CLI layer wraps them into user-facing error variants
/// that name the responsible flag — keeps `crate::bam` independent of
/// any single caller's argument surface.
#[derive(Error, Debug)]
#[non_exhaustive]
pub enum AlignmentIndexError {
    /// An input alignment file has no index next to it.
    #[error(
        "input alignment file '{path}' has no index \
         (looked for '{expected_index_path}')"
    )]
    MissingAlignmentIndex {
        path: PathBuf,
        expected_index_path: PathBuf,
    },

    /// Index construction or write failed (commonly a read-only
    /// directory next to the source file).
    #[error("failed to build alignment index for '{path}': {source}")]
    BuildFailed {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// Loading a previously-built alignment index from disk
    /// failed (corrupted index, permission denied, …). Surfaced
    /// for both CRAM (`.crai`) and BAM (`.csi`/`.bai`) inputs.
    #[error("failed to load alignment index '{index_path}' for '{path}': {source}")]
    LoadFailed {
        path: PathBuf,
        index_path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    /// Input's extension is neither `.cram` nor `.bam`. Pre-flight
    /// does not sniff file magic — an explicit extension is required.
    #[error(
        "input alignment file '{path}' has an unsupported extension \
         (expected .cram or .bam)"
    )]
    UnsupportedExtension { path: PathBuf },
}
