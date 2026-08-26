//! The **psp store**: everything one sample's reads showed at every reference position a
//! run analysed, written once by the locus generator and read back by the cohort merge.
//!
//! A caller opens one file per sample and holds them all open for the whole run, so what
//! one open file costs is multiplied by the cohort size. That is the requirement the
//! format is shaped around: **an open file costs no more than 500 kB of resident memory,
//! and that figure does not grow with the block size, the read depth or the length of the
//! genome** (`doc/devel/ng/spec/psp_file_format.md` §1.1).
//!
//! The way it gets there is to untie two numbers production's `.psp` ties together. There,
//! a block is at once *how far back the compressor may look for a repeated pattern* and
//! *how much a reader must inflate before it can hand out its first record*, so shrinking a
//! block to save memory costs compression and multiplies the index. Here the block stays
//! large — for its ratio and for a small index — while a separate declared number caps the
//! compressor's look-back, and the reader decodes incrementally instead of inflating a
//! whole block.
//!
//! ```text
//! +---------------------------+
//! | header                    |  plain text: magic, length, TOML body, sentinel
//! +---------------------------+
//! | psp block 0               |  records, one compressed stream, each record
//! | psp block 1               |    opening with the head that lets a reader skip it
//! | ...                       |
//! +---------------------------+
//! | block index               |  one entry per psp block
//! +---------------------------+
//! | trailer                   |  the writer's closing payload; may be empty
//! +---------------------------+
//! | footer                    |  fixed size: offsets, counts, checksum, magic last
//! +---------------------------+
//! ```
//!
//! **Three words this module uses differently from production's `src/psp/`**, and a coder
//! moving between the two will get them wrong: what production calls its *metadata
//! section* is the **trailer** here, what production calls its *trailer* is the **footer**
//! here, and *block* here always means a **psp block** — a span of reference and the
//! records in it — never one of zstd's own internal subdivisions.
//!
//! **This module is not a pipeline step.** It is infrastructure beside
//! [`crate::ng::ref_seq`], used by the locus generator to write and by the cohort merge to
//! read, so there is no trait and no competing implementation.
//!
//! Design authority: `doc/devel/ng/spec/psp_file_format.md` (the container),
//! `doc/devel/ng/spec/psp_record_encoding.md` (what a record holds),
//! `doc/devel/ng/spec/psp_chain_id_encoding.md` (the chain ids), and
//! `doc/devel/ng/arch/psp_file_format.md` (the code shape).

use std::path::PathBuf;

use crate::ng::types::GenomeRegion;

pub mod footer;
pub mod header;
pub mod index;
pub mod record;

pub use footer::{FOOTER_BYTES, Footer};
pub use header::{
    ContigEntry, FieldEncoding, FieldName, FieldSpec, Header, Manifest, ParameterValue,
    ReferenceIdentity, WriterProvenance,
};
pub use index::BlockIndexEntry;
pub use record::RecordHead;

// ---------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------

/// Everything that can go wrong reading a psp.
///
/// **Every variant is an input problem, not a bug.** A corrupt or truncated file is data a
/// run was handed, so none of these is a panic and none of them may reach a caller as a
/// half-built record (spec §6.7).
///
/// The variants are separate because **the instruction to whoever sees them differs**:
/// rebuild the file, upgrade the reader, raise a limit, the data is damaged. Collapsing
/// them into one error loses exactly that.
///
/// **Same name as production's `psp::errors::PspReadError`, a different type in a
/// different module.** ng's surface is its own; the two must not be confused at a `use`
/// site.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum PspReadError {
    /// No valid footer, so the writer never finished: the run was interrupted and the file
    /// holds however many blocks reached disk. **Refusing it is the point** — read short,
    /// it is indistinguishable from a sample covering less of the genome.
    #[error("{} has no valid footer — the writer did not finish", path.display())]
    Incomplete { path: PathBuf },

    /// Written by a newer format than this reader understands. Upgrade the reader; the
    /// file is fine.
    ///
    /// **Reached by parsing the version and nothing else**, which is why the header stays
    /// plain text: a reader must be able to learn the version of a file it cannot
    /// otherwise read (spec §3.1).
    #[error(
        "{} is format {}.{}; this reader understands up to {}.{}",
        path.display(), found.0, found.1, supported.0, supported.1
    )]
    UnsupportedVersion {
        path: PathBuf,
        found: (u16, u16),
        supported: (u16, u16),
    },

    /// The file declares a compressor look-back window larger than this reader budgeted
    /// for, so its decoder would have to hold more than it is allowed to.
    ///
    /// **A variant of its own because the fix is a knob, not a rebuild**, and because
    /// zstd's own error here names a code rather than a number anyone can act on.
    #[error(
        "{} needs a {needed_bytes}-byte look-back window; this reader allows {allowed_bytes}",
        path.display()
    )]
    WindowTooLarge {
        path: PathBuf,
        needed_bytes: u64,
        allowed_bytes: u64,
    },

    /// The header is not a header this reader can make sense of: the magic is wrong, the
    /// declared body length disagrees with the sentinel, the TOML does not parse, or a
    /// field in it breaks a rule the format requires — an empty contig list, an MD5 that is
    /// not 32 hex characters, a field encoding with no width.
    ///
    /// **One class rather than the dozen production distinguishes**, because they share an
    /// instruction: the file is damaged, rebuild it. `reason` carries which rule broke.
    ///
    /// **Not the same thing as [`UnsupportedVersion`](Self::UnsupportedVersion)**, and the
    /// order the two are checked in is the point: a file written by a newer format must say
    /// so, not come back as unparseable TOML (spec §3.1).
    #[error("{}: {reason}", path.display())]
    MalformedHeader { path: PathBuf, reason: String },

    /// A block failed to decompress, or a record ran past the end of its block. The file
    /// is damaged.
    #[error("{}: block {block} is corrupt", path.display())]
    CorruptBlock {
        path: PathBuf,
        block: u64,
        #[source]
        source: std::io::Error,
    },

    #[error(transparent)]
    Io(#[from] std::io::Error),
}

/// Everything that can go wrong writing a psp.
///
/// **Same name as production's `psp::errors::PspWriteError`, a different type in a
/// different module** — see [`PspReadError`].
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum PspWriteError {
    /// A record was pushed that starts before the one before it. **Refused rather than
    /// accepted**: coordinate order is what the block index and every seek rest on, and a
    /// file that breaks it seeks wrongly rather than failing (spec §6.3).
    #[error("record at {offered} starts before the previous record at {previous}")]
    OutOfOrder {
        previous: GenomeRegion,
        offered: GenomeRegion,
    },

    /// A header field the writer was handed cannot be written: an empty contig list, a
    /// non-finite parameter value, a version this writer does not produce.
    #[error("header field {field} is not writable: {reason}")]
    InvalidHeaderField { field: String, reason: String },

    /// An append was asked to extend a file whose manifest this writer cannot honour.
    /// **Appending does not rewrite the header**, so the added records must use the
    /// encodings the file already declares (spec §6.4).
    #[error("{} declares a manifest this writer cannot honour: {reason}", path.display())]
    UnsupportedManifest { path: PathBuf, reason: String },

    /// Reopening a finished file failed. **Append and trailer replacement read before they
    /// write** — the footer for the offsets, the header for the manifest — so both surface
    /// the read-side classes (spec §6.7) rather than restating them here.
    #[error(transparent)]
    Reopen(#[from] PspReadError),

    #[error(transparent)]
    Io(#[from] std::io::Error),
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position};

    /// The messages are the contract spec §6.7 tabulates, so they are pinned rather than
    /// left to whatever `thiserror` last rendered. Each has to say what happened *and*
    /// carry the number whoever sees it must act on.
    #[test]
    fn read_errors_name_the_file_and_the_limit_that_was_exceeded() {
        let incomplete = PspReadError::Incomplete {
            path: PathBuf::from("SRR7279481.psp"),
        };
        assert_eq!(
            incomplete.to_string(),
            "SRR7279481.psp has no valid footer — the writer did not finish"
        );

        let version = PspReadError::UnsupportedVersion {
            path: PathBuf::from("SRR7279481.psp"),
            found: (2, 0),
            supported: (1, 0),
        };
        assert_eq!(
            version.to_string(),
            "SRR7279481.psp is format 2.0; this reader understands up to 1.0"
        );

        let window = PspReadError::WindowTooLarge {
            path: PathBuf::from("SRR7279481.psp"),
            needed_bytes: 524_288,
            allowed_bytes: 32_768,
        };
        assert_eq!(
            window.to_string(),
            "SRR7279481.psp needs a 524288-byte look-back window; this reader allows 32768"
        );
    }

    /// An out-of-order push names both regions, because the useful question is which pair
    /// of records disagreed and not merely that some pair did.
    #[test]
    fn an_out_of_order_push_names_both_regions() {
        let at = |start: u64, end: u64| GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        };
        let refused = PspWriteError::OutOfOrder {
            previous: at(1_000, 1_000),
            offered: at(900, 900),
        };
        assert_eq!(
            refused.to_string(),
            "record at contig 0:900-900 starts before the previous record at contig 0:1000-1000"
        );
    }

    /// Reopening for an append or a trailer rewrite reads the footer first, so a file with
    /// no footer must surface as the read-side class through `?` without a second variant
    /// restating it.
    #[test]
    fn a_reopen_carries_the_read_side_error_through() {
        let interrupted = PspReadError::Incomplete {
            path: PathBuf::from("half-written.psp"),
        };
        let as_write_error: PspWriteError = interrupted.into();
        assert_eq!(
            as_write_error.to_string(),
            "half-written.psp has no valid footer — the writer did not finish"
        );
    }
}
