//! Reference-FASTA fetchers. Three types live here:
//!
//! - `StreamingChromRefFetcher` — per-contig, sliding-buffer reader.
//!   The workhorse for monotonic-forward access (cohort var-calling,
//!   PerGroupMerger). Implements the single-contig typed-error
//!   `ChromRefFetcher` trait. Bound to one contig at construction;
//!   serves `fetch(start, length)` and `iter_bases()` through a 1 MB
//!   sliding buffer.
//! - `ManualEvictChromRefFetcher` — per-contig, caller-managed
//!   buffer lifecycle. Designed for genuinely random-access patterns
//!   within a known region (Stage 1 BAQ per-worker instance).
//! - `RepositoryRefFetcher` — multi-contig fetcher used by the
//!   Stage 1 pileup walker. Reads windows out of the shared noodles
//!   `fasta::Repository` the reader already keeps resident (for CRAM
//!   decode and the per-read F1/F3 fetch), so the walker does not open a
//!   second reader over the same FASTA. Implements the multi-contig
//!   typed-error `MultiChromRefFetcher` trait; `Send + Sync` (holds
//!   only an `Arc` clone of the repository).
//!
//! M12 (2026-05-23 code review): the legacy `RefSeqFetcher` trait
//! (`io::Error`-returning, multi-chrom) was retired. All consumers
//! migrated to the typed-error
//! `MultiChromRefFetcher` trait.

use std::io;

// ---------------------------------------------------------------------
// StreamingChromRefFetcher — per-worker, sliding-buffer reader.
// ---------------------------------------------------------------------

// ---------------------------------------------------------------------
// New unified ChromRefFetcher API (in-flight migration).
// ---------------------------------------------------------------------
//
// The trait below is the single-contig typed-error fetcher API
// used by per-chrom consumers (DUST, PerGroupMerger). The
// multi-contig sibling is [`MultiChromRefFetcher`]
// (used by the Stage 1 walker). Both use [`ChromRefFetchError`].
//
// Design choices:
//
// - **No `chrom_id` parameter.** The fetcher is bound to one contig
//   at construction; multi-contig consumers build one fetcher per
//   contig rather than passing chrom ids per call.
// - **Typed error.** `ChromRefFetchError` distinguishes "I/O
//   failure", "out of contig bounds", and "out of access pattern"
//   so consumers can route them differently.
// - **`OutOfPattern` is the load-bearing contract.** A streaming
//   impl errors when a `fetch` start lies before the current sliding
//   buffer's origin. The contract is "monotonic non-decreasing
//   `start` across `fetch` calls within a phase"; `iter_bases()`
//   defines its own phase (resets the buffer state on construction
//   and on Drop). If the error fires in real code, that consumer
//   needs a random-access impl rather than a streaming one.
// - **`length()`** is exposed on the trait so consumers can use
//   the contig's declared length without threading it through
//   constructor calls.

/// Why a `ChromRefFetcher` call failed.
// M2 (2026-05-23 code review): #[non_exhaustive] so adding a future
// variant (e.g. an InvalidIndex variant for B1's validation if it
// graduates from io::Error to a typed variant) doesn't break callers.
// The enum lives behind a `pub` surface but the crate is unpublished;
// in-tree matchers must accept the marker today so they pay the
// add-a-variant cost upfront.
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
    /// Streaming impl only. The requested `start_1based` lies
    /// before the current sliding buffer's origin; the consumer's
    /// access pattern is not monotonic non-decreasing. Production
    /// callers (DUST + PerGroupMerger) satisfy the contract by
    /// construction; this error firing in real code means a
    /// non-monotonic consumer that needs a random-access fetcher
    /// impl with the same trait.
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

// ---------------------------------------------------------------------
// ManualEvictChromRefFetcher — caller-managed buffer lifecycle.
// ---------------------------------------------------------------------
//
// Designed for consumers whose access pattern is genuinely random
// access *within a known region* — Stage 1 BAQ is the motivating
// example. Each rayon worker owns its own instance via `map_init`
// and processes reads that may visit positions in any order within
// the current chunk's coord range. The worker calls `evict_before`
// at deterministic points to keep memory bounded.
//
// Contract:
//
// - `fetch(start, length)` returns a borrowed `&[u8]` slice. If the
//   range isn't covered, the buffer extends in whichever direction
//   is needed — forward by appending fresh bytes, backward by
//   prepending. Never shrinks on its own.
// - `evict_before(pos)` is the *only* shrinking operation. The
//   caller decides when bytes are no longer needed. CRAM is
//   coordinate-sorted, so calling `evict_before(read.pos)` after
//   each read keeps the buffer at "current-position window plus
//   whatever the next read extends forward."
//
// Trade-off vs [`StreamingChromRefFetcher`]:
//
// - Streaming forces monotonic forward access (small look-back
//   within the buffer is fine; large backward jumps are
//   `OutOfPattern`). Good for cohort var-calling.
// - ManualEvict tolerates any access pattern within the resident
//   range and trusts the caller to keep memory bounded. Good for
//   rayon-parallel random-access-within-chunk consumers. Worst case
//   (no `evict_before` calls + adversarial input) is one whole
//   contig resident — acceptable because the consumer controls it.
//
// Not `Sync`. Each thread holds its own instance; concurrent
// fetches from multiple threads on a single instance are forbidden.
// `Send` because every field is `Send`.

/// Read `n_bases` uppercased, newline-stripped bases from `reader`
/// (currently positioned at the file offset for those bases) and
/// append them to `dst`. Returns an `UnexpectedEof` if the source
/// runs out before `n_bases` bases are gathered. Shared by both
/// streaming and manual-evict fetchers.
/// Canonicalise a single FASTA byte to SAM-spec uppercase `ACGT`/`N`.
///
/// M26 (2026-05-23 code review): the `ChromRefFetcher::fetch` trait
/// docs promise output is "uppercased bases, SAM-spec canonical
/// (`A`/`C`/`G`/`T`/`N`) regardless of how the FASTA is masked on
/// disk". The previous implementation used only
/// `b.to_ascii_uppercase()`, which leaves IUPAC ambiguity codes
/// (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`), the RNA `U`,
/// gap characters (`-`, `.`), and any other ASCII byte unchanged.
/// Downstream DUST + BAQ + PerGroupMerger then saw raw IUPAC bytes
/// that they were not built to handle.
///
/// This function folds any non-`ACGT` byte to `N` AFTER
/// uppercasing — so `a`/`c`/`g`/`t` (soft-masked) map to `A`/`C`/
/// `G`/`T`, and everything else (`R`, `n`, `N`, `U`, `-`, `.`,
/// stray whitespace bytes that survived the newline filter) maps
/// to `N`.
#[inline]
pub(crate) fn canonicalise(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' | b'C' | b'G' | b'T' => b.to_ascii_uppercase(),
        _ => b'N',
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {

    // ---------------------------------------------------------------
    // StreamingChromRefFetcher tests
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // ChromRefFetcher (new API) tests
    // ---------------------------------------------------------------
    //
    // These tests exercise the new-API methods: `for_contig`
    // construction, `length()`, `fetch(start, length)` (without
    // chrom_id), `iter_bases()` (with buffer reset semantics), and
    // the `OutOfPattern` contract on backward jumps beyond the
    // buffer.

    // ---------------------------------------------------------------
    // ManualEvictChromRefFetcher tests
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // B1 (2026-05-23 code review): .fai validation tests
    // ---------------------------------------------------------------
    //
    // The .fai file format permits attacker-influenced values
    // (noodles_fasta::fai::Record parses without validating
    // `line_bases > 0` or `line_width >= line_bases`). The
    // `--reference` CLI argument is a real attacker surface.
    // Pre-B1 a malformed .fai with `line_bases = 0` would parse,
    // and the first `fetch` call would panic inside
    // ContigFai::base_to_file_offset on integer division by zero.

    // ---------------------------------------------------------------
    // M26 (2026-05-23 code review): IUPAC + soft-mask canonicalisation
    // ---------------------------------------------------------------

    // ---------------------------------------------------------------
    // B3 + M22 + M24 + M25 (2026-05-23 code review): contract /
    // edge-case tests that were missing.
    // ---------------------------------------------------------------

    // -----------------------------------------------------------------
    // RepositoryRefFetcher
    // -----------------------------------------------------------------
}
