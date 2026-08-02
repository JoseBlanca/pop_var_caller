//! **The alignment cursor — a reader that stays where it is.** So far only its errors:
//! [`CursorError`]. The cursor itself lands in Milestone B, and the record readers beneath it
//! in the next step; the rest of this module doc is the problem being built for.
//!
//! Reading a sorted alignment file today opens a fresh query per region: the index is
//! consulted, the reader seeks, and the block it lands in is decompressed and decoded — and
//! then the next region, 390 bases along, does all of it again. Measured on chromosome 21
//! of a tandem-repeat-targeted HG002 run, **82 % of the seeks land in the block the reader
//! already holds**, and the same 35,228 records — that is the probe's **whole-contig** mode;
//! the typed-region walk counts 34,633, and spec §11.5 requires the mode to be named because
//! both figures are in circulation — are decoded 1,067,729 times. Consecutive
//! queries overlap by about 93 %, because the caller widens every region by 5,000 bases
//! against regions averaging 390.
//!
//! A cursor is the reader kept between regions instead: positioned in one chromosome of one
//! file, holding the reads it has already decoded *and filtered*, and handing them back when
//! the next region can use them.
//!
//! Design: `doc/devel/ng/spec/alignment_cursor.md` (what and why) and
//! `doc/devel/ng/arch/alignment_cursor.md` (types and interfaces). Build order in
//! `doc/devel/ng/impl_plan/alignment_cursor.md`.
//!
//! # What is here so far
//!
//! Only the errors. The record readers beneath the cursor arrive in the next step of this
//! milestone; the cursor itself and the sample-level merge above it in later ones. Milestone A
//! is types and the instrument, and changes no behaviour.

use std::path::Path;
use std::sync::Arc;

use crate::ng::types::ContigId;

/// What can go wrong once a cursor exists.
///
/// Two conditions, and they are not peers. `Io` is the file failing under a read that was
/// asked for correctly; [`WrongChromosome`](Self::WrongChromosome) is a **caller bug** — a
/// guard, not a step in normal control flow. Correct code compares against the cursor's own
/// chromosome first and never sees it.
///
/// There is deliberately **no ordering variant**. Within its chromosome a cursor answers any
/// region in any order and the answer is always right; only the *speed* depends on how close
/// the region is to the last one (spec §4). Requiring regions to move forward was weighed
/// and rejected: a backward jump costs a seek and a block, which is what every request costs
/// today, so the restriction would protect against almost nothing while putting error
/// handling at every call site. Requiring a new cursor per chromosome was kept, because
/// there the cost prevented is chromosome-sized — on CRAM, re-reading hundreds of megabytes
/// of reference bases — and because nothing in a cursor survives a chromosome change anyway.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CursorError {
    /// A region on a chromosome this cursor does not cover. Make a cursor for that
    /// chromosome; this one is unharmed and still good for its own.
    ///
    /// **The contigs are reported as the numbers they are, and the file is named so the
    /// numbers mean something.** A cursor holds no name table, so it cannot say `chr21` — but
    /// an index is only interpretable against a particular table, and the path says which:
    /// [`AlignmentFile::contigs`](super::open_bam::AlignmentFile::contigs) on that file turns
    /// both numbers into names. Without it, a run holding one cursor per chromosome per
    /// generator per worker — 32 of them at one file per sample, 320 at ten — reports two bare
    /// integers and no way to tell which of the 320 produced them.
    #[error(
        "cursor on '{}' covers contig {} but the region is on contig {}",
        path.display(),
        cursor_contig.get(),
        requested_contig.get()
    )]
    WrongChromosome {
        path: Arc<Path>,
        /// The chromosome this cursor was made for.
        cursor_contig: ContigId,
        /// The chromosome the region the caller asked for lies on.
        requested_contig: ContigId,
    },

    /// Reading the next record from the file failed.
    #[error("reading alignment file '{}' failed", path.display())]
    ReadRecord {
        path: Arc<Path>,
        #[source]
        source: std::io::Error,
    },
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The message names both contigs, names them apart, **and names the file** — two bare
    /// integers are only meaningful against a particular contig table, and a run holds up to
    /// 320 cursors over many files.
    #[test]
    fn the_wrong_chromosome_message_names_the_file_and_both_contigs() {
        let error = CursorError::WrongChromosome {
            path: Arc::from(Path::new("/data/sample.bam")),
            cursor_contig: ContigId(20),
            requested_contig: ContigId(7),
        };

        assert_eq!(
            error.to_string(),
            "cursor on '/data/sample.bam' covers contig 20 but the region is on contig 7",
        );
    }

    /// The path is rendered, not debug-printed: `Path` has no `Display`, so the naive
    /// `{path}` does not compile and the naive `{path:?}` prints quotes and escapes into an
    /// operator-facing message.
    #[test]
    fn the_read_failure_message_renders_the_path_and_keeps_the_cause() {
        let error = CursorError::ReadRecord {
            path: Arc::from(Path::new("/data/sample.bam")),
            source: std::io::Error::new(std::io::ErrorKind::UnexpectedEof, "truncated block"),
        };

        assert_eq!(
            error.to_string(),
            "reading alignment file '/data/sample.bam' failed",
        );
        // The cause survives as a `#[source]`, so the renderer that walks the chain reaches
        // it — without it, "reading … failed" would be the whole story.
        let source = std::error::Error::source(&error).expect("the io error is the source");
        assert_eq!(source.to_string(), "truncated block");
    }
}
