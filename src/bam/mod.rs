//! Mapped-read file input — BAM and CRAM: the shared read vocabulary and the index check.
//!
//! The reader itself is [`crate::read::input`]. What stays here is what it and the rest of the
//! crate share:
//!
//! - [`alignment_input`] — the vocabulary a decoded read is spelled in: [`CigarOp`]
//!   (used crate-wide), the SAM `FLAG_*` bits, the read filters' `DEFAULT_*` thresholds,
//!   [`MappedRead`], and the per-read helpers the filters and the decoder share.
//! - [`errors`] — [`errors::AlignmentIndexError`], raised by the index check.
//! - [`index_preflight`] — `.crai` / `.bai` / `.csi` detection and opt-in building, for the
//!   commands that need per-contig random access.
//!
//! A decoded read becomes a [`MappedRead`] and then ng's
//! [`PreparedRead`](crate::read::PreparedRead) in read preparation.
//!
//! [`CigarOp`]: alignment_input::CigarOp
//! [`MappedRead`]: alignment_input::MappedRead

pub mod alignment_input;
// Test fixtures: synthetic FASTA and CRAM files written to a tempdir. Used by the read input's
// tests, `fasta`'s and the commands'.
#[cfg(test)]
pub(crate) mod cram_files;
pub mod errors;
pub mod index_preflight;
