//! What is left of the reference-FASTA fetchers: the error a reference read raises
//! ([`ChromRefFetchError`]) and the rule that canonicalises a FASTA byte ([`canonicalise`]).
//!
//! The fetchers themselves (`StreamingChromRefFetcher`, `ManualEvictChromRefFetcher`,
//! `RepositoryRefFetcher`) and the `ChromRefFetcher` trait served the deleted production caller
//! and are gone. The readers this caller uses are [`crate::raw_chrom_reader`] and the
//! [`crate::ref_seq`] implementations; both raise this error, and `ref_seq` canonicalises with
//! this function so its bytes match what the old fetchers returned.

use std::io;

/// Why a reference read failed.
// `#[non_exhaustive]` so a new variant breaks no caller: every in-tree match already carries the
// wildcard arm a new variant would need.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ChromRefFetchError {
    /// The requested `[start, start+length)` window exceeds the
    /// bound contig's total length.
    #[error("fetch [{start}, {end}) past contig {chrom_name} length {chrom_length}")]
    OutOfBounds {
        chrom_name: String,
        chrom_length: u32,
        start: u32,
        end: u32,
    },
    /// `start_1based` was 0 — 1-based coordinate contract violated.
    #[error("start_1based must be >= 1")]
    InvalidStart,
    /// A read before the start of a streaming fetcher's sliding buffer. Only the deleted
    /// `StreamingChromRefFetcher` raised it; no reader in this tree does.
    #[error(
        "fetch at base {requested_start} lies before the streamer's current buffer (origin = {buffer_origin}). \
         StreamingChromRefFetcher requires monotonic non-decreasing access"
    )]
    OutOfPattern {
        requested_start: u32,
        buffer_origin: u32,
    },
    /// Underlying I/O failure: file open, seek, read, or `.fai`
    /// parse.
    #[error("I/O failure on contig {chrom_name}: {source}")]
    Io {
        chrom_name: String,
        #[source]
        source: io::Error,
    },
}

/// Canonicalise a single FASTA byte to SAM-spec uppercase `ACGT`/`N`.
///
/// Uppercasing alone would leave IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`,
/// `D`, `H`, `V`), the RNA `U`, gap characters (`-`, `.`) and any other ASCII byte unchanged,
/// and a consumer comparing bases against `ACGT` is not built for them. So this function folds
/// any non-`ACGT` byte to `N` *after* uppercasing: `a`/`c`/`g`/`t` (soft-masked) map to
/// `A`/`C`/`G`/`T`, and everything else (`R`, `n`, `N`, `U`, `-`, `.`, stray whitespace bytes)
/// maps to `N`.
#[inline]
pub(crate) fn canonicalise(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => b.to_ascii_uppercase(),
        _ => b'N',
    }
}
