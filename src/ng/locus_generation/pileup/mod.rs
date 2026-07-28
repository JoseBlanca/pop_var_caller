//! ng's generic locus generator — the pileup walk, copied from production.
//!
//! A **folder**, unlike its `ssr.rs` sibling, because it holds the whole copied
//! walker plus (later) the generator that wraps it. Production is not edited:
//! ng copies rather than reaches in, so there is no visibility lift and no field
//! on a frozen type (`doc/devel/ng/spec/locus_generation_pileup.md` §3,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*).
//!
//! # Seven files, verbatim
//!
//! [`genome_walk`], [`open_record`], [`cigar_cursor`], [`decompose`],
//! [`active_read_set`], [`chain_id_allocator`] and [`errors`] are transcribed
//! from `src/pileup/walker/` **unchanged**, and are proven to compute exactly
//! what production computes before a line of them is edited (spec §3, the
//! stage-1 differential). The rule that paid three times on this branch is
//! *transcribe first, change second*: a copy that is provably production is the
//! baseline every later change is measured against, and without it the
//! generator's deliberate divergences could not be told from transcription
//! slips.
//!
//! So the copies still emit production's [`PileupRecord`](crate::pileup_record::PileupRecord),
//! not ng's `SampleLocusObservations`. That changes in the next plan.
//!
//! **One file is renamed on the way in: `driver.rs` → [`genome_walk`]** — it is
//! the only one of the seven named for a *role* rather than for what it owns,
//! and "driver" answers *driver of what?* with nothing. `genome_walk` names the
//! one job that file has and the others do not: advancing a position cursor
//! along genome coordinates over an active read set. The **type** keeps
//! production's `PileupWalker`, so the differential reads as a straight
//! comparison. *(A walk covers one region, not the genome; `genome_walk` names
//! the axis it advances along, not the extent.)*
//!
//! # What this module re-exports, and why
//!
//! The copied files reach their shared vocabulary through `super::` — that is
//! how they were written against `pileup/walker/mod.rs`, and leaving those paths
//! alone is what keeps the transcription verbatim. This module therefore stands
//! in for production's `walker/mod.rs`, drawing each name from wherever it now
//! lives:
//!
//! - [`PreparedRead`], [`MateRole`], [`ReadLengthError`] — **ng's**, copied and
//!   extended with `read_group` (spec §6). They are the reason the whole walker
//!   is copied: every one of the seven names `PreparedRead` in its signatures.
//! - [`CigarOp`], [`WalkerConfig`] and the `DEFAULT_*` constants — **production's,
//!   reused as-is**. ng does not modify them, so it does not copy them; the
//!   constants are reached by name rather than by literal so there is one source
//!   of truth until ng deliberately diverges.

mod active_read_set;
mod chain_id_allocator;
mod cigar_cursor;
mod decompose;
mod errors;
mod genome_walk;
mod open_record;

// The vocabulary the copied files resolve through `super::`. Production's own
// `walker/mod.rs` declares exactly these names around exactly these modules.
pub use crate::ng::read::prepared_read::{MateRole, PreparedRead, ReadLengthError};
pub use crate::pileup::walker::{
    CigarOp, DEFAULT_MATE_LOOKUP_WINDOW, DEFAULT_MAX_INDEL_COLUMN_DEPTH, DEFAULT_MAX_RECORD_SPAN,
    DEFAULT_MAX_SNP_COLUMN_DEPTH, WalkerConfig,
};

pub use chain_id_allocator::DEFAULT_MAX_ACTIVE_READS;
pub use errors::WalkerError;
pub use genome_walk::{PileupWalker, RunSummary, run};

use std::io;

use crate::fasta::{ChromRefFetchError, MultiChromRefFetcher};
use crate::ng::ref_seq::{RefSeq, RefSeqError};
use crate::ng::types::ContigId;

// ---------------------------------------------------------------------
// The reference shim
// ---------------------------------------------------------------------

/// Serve ng's [`RefSeq`] through the [`MultiChromRefFetcher`] the copied walker
/// asks for.
///
/// **Semantically empty, and that is a checked claim rather than a hopeful one.**
/// The walker's contract is "uppercase ASCII over `{A,C,G,T,N}`, canonicalised by
/// the fetcher implementation" ([`MultiChromRefFetcher::fetch`]), and that is
/// exactly what `RefSeq::fetch_into` promises — it folds the stored bytes to
/// canonical, so a soft-masked or ambiguity-coded reference arrives at the walk
/// already normalised. Read preparation's case divergence (`read_preparation.md`
/// §6) does **not** recur here: this shim moves bytes and decides nothing
/// (`locus_generation_pileup.md` §3). `the_shim_canonicalises_what_it_serves`
/// pins it against a reference built to need every fold.
///
/// It is a newtype rather than a blanket `impl<R: RefSeq> MultiChromRefFetcher
/// for R` because the two traits are both ours to implement and a blanket impl
/// would claim every present and future `RefSeq` for a role most of them will
/// never play.
///
/// # What is lost, and it is only error text
///
/// [`ChromRefFetchError`] names contigs by **name** where [`RefSeqError`] names
/// them by [`ContigId`] — and `RefSeq` alone cannot resolve one to the other
/// (the contig table is a separate capability, deliberately). So the id is
/// rendered in the name's place. Nothing computes on it: the walker wraps every
/// fetch failure in `WalkerError::Fasta { chrom_id, .. }`, which carries the id
/// already. The `u64 → u32` narrowings in the out-of-bounds arm are the same
/// kind of thing — message fields on a window that has already been rejected.
#[derive(Debug, Clone, Copy)]
pub struct RefSeqFetcher<R: RefSeq>(pub R);

impl<R: RefSeq> MultiChromRefFetcher for RefSeqFetcher<R> {
    fn fetch(
        &self,
        chrom_id: u32,
        start_1based: u32,
        length: u32,
    ) -> Result<Vec<u8>, ChromRefFetchError> {
        self.0
            .fetch(
                ContigId(chrom_id),
                u64::from(start_1based),
                u64::from(length),
            )
            .map_err(|error| to_fetch_error(chrom_id, error))
    }
}

/// Translate a [`RefSeqError`] into the walker's [`ChromRefFetchError`], variant
/// for variant.
///
/// An exhaustive `match`, so a variant added to either enum stops this compiling
/// rather than falling into a catch-all that reports the wrong failure mode. The
/// one arm that is not a rename is `UnknownContig`, which production's trait
/// documents as belonging to `Io` ("`Io` on any underlying FASTA I/O failure or
/// unknown `chrom_id`") — so it is reported there, with a synthesised
/// `NotFound` source that says which id was unknown.
fn to_fetch_error(chrom_id: u32, error: RefSeqError) -> ChromRefFetchError {
    // Contig names are what `ChromRefFetchError` wants and what `RefSeq` cannot
    // give; the id stands in for it. Message text only — see the type's note.
    let chrom_name = || format!("chrom_id {chrom_id}");
    // The window was rejected, so these are diagnostic strings, not arithmetic.
    // Saturating rather than truncating: a wrapped number would misreport the
    // very window that was refused.
    let narrow = |value: u64| u32::try_from(value).unwrap_or(u32::MAX);

    match error {
        RefSeqError::OutOfBounds {
            contig: _,
            contig_length,
            start,
            end,
        } => ChromRefFetchError::OutOfBounds {
            chrom_name: chrom_name(),
            chrom_length: narrow(contig_length),
            start: narrow(start),
            end: narrow(end),
        },
        RefSeqError::InvalidStart => ChromRefFetchError::InvalidStart,
        RefSeqError::UnknownContig(contig) => ChromRefFetchError::Io {
            chrom_name: chrom_name(),
            source: io::Error::new(
                io::ErrorKind::NotFound,
                format!("no contig with id {} in the reference", contig.0),
            ),
        },
        RefSeqError::Io { contig: _, source } => ChromRefFetchError::Io {
            chrom_name: chrom_name(),
            source,
        },
    }
}

// `shim_tests`, not `tests`: the copied `walker/tests.rs` lands as this module's
// `tests` child (A4), and mirroring production's module names is what makes the
// two suites comparable name for name.
#[cfg(test)]
mod shim_tests {
    use super::*;
    use crate::ng::ref_seq::InMemoryRefSeq;

    /// Two contigs, the second one written the way a real assembly writes
    /// repeats and ambiguity: lower case, and a code that is neither `ACGT` nor
    /// `N`. Both are what `RefSeq` folds and what the walker must never see.
    fn reference() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![b"ACGTACGTAC".to_vec(), b"acgtRYacgt".to_vec()])
    }

    #[test]
    fn the_shim_serves_the_window_it_was_asked_for() {
        let fetcher = RefSeqFetcher(reference());
        assert_eq!(fetcher.fetch(0, 1, 4).expect("in range"), b"ACGT");
        // 1-based and inclusive of `start`: base 5 is the fifth, not the sixth.
        assert_eq!(fetcher.fetch(0, 5, 3).expect("in range"), b"ACG");
        assert_eq!(fetcher.fetch(0, 10, 1).expect("in range"), b"C");
    }

    /// **The claim that this shim is semantically empty, checked.** The walker's
    /// contract is canonical uppercase `{A,C,G,T,N}`; a reference that is neither
    /// uppercase nor `ACGTN` comes out as both. Without this, "the two contracts
    /// already agree" would be a doc comment and nothing else.
    #[test]
    fn the_shim_canonicalises_what_it_serves() {
        let fetcher = RefSeqFetcher(reference());
        let bases = fetcher.fetch(1, 1, 10).expect("in range");
        assert_eq!(bases, b"ACGTNNACGT", "lower case up, non-ACGT folded to N");
        assert!(
            bases
                .iter()
                .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')),
            "every byte is inside the walker's alphabet"
        );
    }

    /// A window past the contig's end is refused, and reported as the walker's
    /// own out-of-bounds — not folded into an I/O failure, which is what the
    /// walker's error routing distinguishes.
    #[test]
    fn a_window_past_the_contig_end_is_out_of_bounds() {
        let fetcher = RefSeqFetcher(reference());
        let error = fetcher.fetch(0, 9, 5).expect_err("only 10 bases exist");
        assert!(
            matches!(
                error,
                ChromRefFetchError::OutOfBounds {
                    chrom_length: 10,
                    start: 9,
                    end: 14,
                    ..
                }
            ),
            "unexpected: {error:?}"
        );
    }

    /// The 1-based contract, carried across unchanged: `start_1based == 0` is
    /// its own variant on both sides, never an out-of-bounds.
    #[test]
    fn a_zero_start_keeps_its_own_variant() {
        let fetcher = RefSeqFetcher(reference());
        let error = fetcher.fetch(0, 0, 4).expect_err("1-based coordinates");
        assert!(
            matches!(error, ChromRefFetchError::InvalidStart),
            "unexpected: {error:?}"
        );
    }

    /// An unknown contig is an `Io` failure, which is where production's trait
    /// documentation puts it — and the message says which id was unknown, since
    /// the name it would rather print does not exist.
    #[test]
    fn an_unknown_contig_is_reported_as_io_naming_the_id() {
        let fetcher = RefSeqFetcher(reference());
        let error = fetcher.fetch(7, 1, 4).expect_err("only two contigs exist");
        let ChromRefFetchError::Io { source, .. } = &error else {
            panic!("unexpected: {error:?}");
        };
        assert_eq!(source.kind(), io::ErrorKind::NotFound);
        assert!(
            source.to_string().contains('7'),
            "the message names the id: {source}"
        );
    }
}
