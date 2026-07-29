//! ng's generic locus generator — the pileup walk, copied from production.
//!
//! A **folder**, unlike its `ssr.rs` sibling, because it holds the whole copied
//! walker plus (later) the generator that wraps it. Production is not edited:
//! ng copies rather than reaches in, so there is no visibility lift and no field
//! on a frozen type (`doc/devel/ng/spec/locus_generation_pileup.md` §3,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*).
//!
//! # Eight files, verbatim
//!
//! [`genome_walk`], [`open_record`], [`cigar_cursor`], [`decompose`],
//! [`active_read_set`], [`chain_id_allocator`] and [`errors`] — plus `tests.rs`,
//! production's own end-to-end suite, which is Milestone A's gate — are
//! transcribed from `src/pileup/walker/` **unchanged**. `copy_fidelity.rs` checks
//! that textually, over all eight, so "verbatim" is a property with a test rather
//! than a claim in this comment. The rule that paid three times on this branch is
//! *transcribe first, change second*: a copy that is provably production is the
//! baseline every later change is measured against, and without it the
//! generator's deliberate divergences could not be told from transcription
//! slips.
//!
//! What the copy is proven to *compute* is Milestone B's differential; what it is
//! proven to *be* is `copy_fidelity.rs`. The two are different claims and only the
//! second one has landed.
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
//! - [`CigarOp`], [`WalkerConfig`] and two of the `DEFAULT_*` constants —
//!   **production's, reused as-is**. ng does not modify them, so it does not copy
//!   them; they are reached by name rather than by literal so there is one source
//!   of truth until ng deliberately diverges.
//!   **One constant is ng's, forced by the verbatim rule:**
//!   [`DEFAULT_MAX_ACTIVE_READS`] is declared *inside* `chain_id_allocator.rs`, so
//!   the copy brought its own — two definitions of `4096` now exist, and
//!   `chain_id_allocator.rs`'s `Self::with_caps(DEFAULT_MAX_ACTIVE_READS,
//!   super::DEFAULT_MATE_LOOKUP_WINDOW)` reaches ng's for the first and production's
//!   for the second. `the_copied_active_reads_cap_is_still_productions` pins the two
//!   equal; when they are deliberately allowed to differ, that test is what says so.
//!
//! The vocabulary is bound `pub(crate)`, not `pub`: it is an internal aid for the
//! copies, not an ng-flavoured public alias for production's types. A consumer that
//! wants `CigarOp` should say `crate::pileup::walker::CigarOp` and see whose it is —
//! which matters from plan 3 on, when ng's walker starts to diverge and two live
//! paths to one name would stop being harmless.

use std::io;

use crate::fasta::{ChromRefFetchError, MultiChromRefFetcher};
use crate::ng::ref_seq::{RefSeq, RefSeqError};
use crate::ng::types::ContigId;

mod active_read_set;
mod chain_id_allocator;
mod cigar_cursor;
mod decompose;
mod errors;
mod genome_walk;
mod open_record;

/// Production's own end-to-end suite, copied verbatim and run against the copy —
/// Milestone A's gate (spec §12). `pub(crate)` because production's declaration is,
/// so the two module trees stay comparable; its `MockFasta` / `snp_read` /
/// `paired_snp_reads` fixtures are reachable from B1's parity harness either way,
/// that being a descendant of this module.
#[cfg(test)]
pub(crate) mod tests;

/// **ng's, not a copy** — the textual check that the eight copies are still
/// production's, from outside the files it checks (spec §3).
#[cfg(test)]
mod copy_fidelity;

/// **ng's, not a copy** — the stage-1 differential: the two walkers compute the
/// same thing, over one read stream (spec §3, §13.1). `copy_fidelity` says the
/// copy *is* production's; this says it *does* what production's does.
#[cfg(test)]
mod parity;

// The vocabulary the copied files resolve through `super::`. Production's own
// `walker/mod.rs` declares the same names, around the same modules bar one: its
// `pub(crate) mod indel_norm` is not copied, because none of the seven reaches it —
// its consumers are production's `read_processor.rs` and ng's own `src/ng/alignment/`,
// neither of which is part of the walk.
//
// `DEFAULT_MAX_SNP_COLUMN_DEPTH` / `DEFAULT_MAX_INDEL_COLUMN_DEPTH` are deliberately
// absent: no copied file reaches them through `super::` (production's `WalkerConfig`
// holds the column caps and comes with its own defaults), and an unused re-export is
// invisible to the compiler, so carrying them would be shape without substance.
// `PileupGeneratorConfig` names them in plan 3, from production directly.
pub(crate) use crate::ng::read::prepared_read::{MateRole, PreparedRead, ReadLengthError};
pub(crate) use crate::pileup::walker::{CigarOp, DEFAULT_MAX_RECORD_SPAN, WalkerConfig};
// `#[cfg(test)]` because that is where the one copy reaching it lives:
// `ChainIdAllocator::new` is itself `#[cfg(test)]` (production code calls `with_caps`
// from `run`). Left ungated it is an unused import in a non-test build — which the
// previous `pub use` hid, since a `pub` re-export is never reported unused.
#[cfg(test)]
pub(crate) use crate::pileup::walker::DEFAULT_MATE_LOOKUP_WINDOW;

// This module's own surface, as against the vocabulary above.
pub use chain_id_allocator::DEFAULT_MAX_ACTIVE_READS;
pub use errors::WalkerError;
pub use genome_walk::{PileupWalker, RunSummary, run};

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
///
/// # What is *not* adapted
///
/// [`MultiChromRefFetcher::iter_bases`] is left on its inherited default, which
/// materialises a whole contig via `fetch(chrom_id, 1, length)`. Nothing in the walk
/// calls it — `open_record.rs`'s two `fetch` sites are the only consumers — so it is
/// unadapted and untested rather than considered and kept. It is a trap for the next
/// caller: reached through a [`WindowedRefSeq`](crate::ng::ref_seq::WindowedRefSeq) it
/// would produce exactly the whole-contig residency that reference exists to avoid,
/// with no error. Override it against `RefSeq::fetch_into` before wiring one up.
///
/// `pub(crate)`, matching arch §1.3's private newtype: B1's parity harness is a sibling
/// inside `pileup/` and needs nothing wider, and a `pub` tuple field would let a caller
/// reach past the shim to the reference underneath — the opposite of "moves bytes and
/// decides nothing". No `Copy`: no `RefSeq` implementation in the tree is `Copy` (they
/// own `Vec`s, `RefCell`s or an `Arc`), so the derive would be inert while advertising
/// bitwise-copy semantics for an accessor of unbounded size. The `R: RefSeq` bound lives
/// on the impl below, not here, so it is stated where it does work.
///
/// # Temporary — plan 3's A0 deletes this (owner, 2026-07-29)
///
/// This type exists for one reason: the copies were transcribed **verbatim**, so their
/// signatures are production's and ask for `MultiChromRefFetcher`. That is a consequence of
/// the copy, not a design choice, and it stops being true the moment the two walkers
/// diverge. Plan 3's **first** step — while the stage-1 differential can still prove a
/// refactor free — has `open_record.rs` take a `RefSeq` directly, at which point this type,
/// [`to_chrom_ref_fetch_error`] and both of the lossy spots documented above go away, and
/// `fetch_into` replaces `fetch` (arch §4's allocation note).
///
/// The name is also one this codebase **deliberately retired** for a different concept
/// (`fasta/fetcher.rs:20-23`, a 2026-05-23 review); deleting it settles that too, where a
/// rename would only have moved it.
///
/// **Dead outside tests until plan 3**, which is what the `expect` records: this
/// milestone builds the shim and proves it, and `PileupGenerator` is the first
/// non-test consumer. The attribute is a compiler backstop that forces its own
/// removal the moment that lands — an `allow` would simply go stale.
#[derive(Debug, Clone)]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "plan 3's PileupGenerator is the first
    non-test consumer; this expect must be deleted when it lands"
    )
)]
pub(crate) struct RefSeqFetcher<R>(pub(crate) R);

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
            .map_err(|error| to_chrom_ref_fetch_error(chrom_id, error))
    }
}

/// Translate a [`RefSeqError`] into the walker's [`ChromRefFetchError`], variant
/// for variant.
///
/// An exhaustive `match` **on [`RefSeqError`]**, so a variant added *there* stops this
/// compiling rather than falling into a catch-all that reports the wrong failure mode.
/// *The constructed side is not checked the same way, and cannot be:* an into-mapping
/// owes no variant coverage. `ChromRefFetchError` already carries a variant this never
/// mints — `OutOfPattern`, which belongs to the streaming fetcher — and a new one would
/// go unproduced in silence.
///
/// The one arm that is not a rename is `UnknownContig`, which production's trait
/// documents as belonging to `Io` ("`Io` on any underlying FASTA I/O failure or
/// unknown `chrom_id`") — so it is reported there, with a synthesised
/// `NotFound` source that says which id was unknown.
///
/// Dead outside tests for as long as [`RefSeqFetcher`] is — same backstop, same reason.
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "reached only through RefSeqFetcher,
    which plan 3's PileupGenerator is the first non-test consumer of"
    )
)]
fn to_chrom_ref_fetch_error(chrom_id: u32, error: RefSeqError) -> ChromRefFetchError {
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

// `ref_seq_fetcher_tests`, not `tests`: the copied `walker/tests.rs` lands as this
// module's `tests` child (A4), and mirroring production's module names is what makes
// the two suites comparable name for name. Named for its subject rather than for the
// pattern, following `mod baq_tests`.
#[cfg(test)]
mod ref_seq_fetcher_tests {
    use super::*;
    use crate::ng::ref_seq::InMemoryRefSeq;

    /// Two contigs, the second one written the way a real assembly writes
    /// repeats and ambiguity: lower case, and a code that is neither `ACGT` nor
    /// `N`. Both are what `RefSeq` folds and what the walker must never see.
    fn reference() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![b"ACGTACGTAC".to_vec(), b"acgtRYacgt".to_vec()])
    }

    /// A reference whose every read fails — the one [`RefSeqError`] arm no fixture
    /// built from bytes in memory can produce, and the one a user actually meets
    /// (a truncated or absent FASTA).
    struct BrokenRefSeq;
    impl RefSeq for BrokenRefSeq {
        fn fetch_into(
            &self,
            contig: ContigId,
            _start_1based: u64,
            _length: u64,
            _dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            Err(RefSeqError::Io {
                contig,
                source: io::Error::new(io::ErrorKind::UnexpectedEof, "truncated fasta"),
            })
        }
    }

    #[test]
    fn the_shim_serves_the_window_it_was_asked_for() {
        let fetcher = RefSeqFetcher(reference());
        assert_eq!(fetcher.fetch(0, 1, 4).expect("in range"), b"ACGT");
        // 1-based and inclusive of `start`: base 5 is the fifth, not the sixth.
        assert_eq!(fetcher.fetch(0, 5, 3).expect("in range"), b"ACG");
        assert_eq!(fetcher.fetch(0, 10, 1).expect("in range"), b"C");
    }

    /// A zero-length window is a legal question with an empty answer — including at
    /// the position one past the contig's last base, which is where the walker's
    /// **exclusive** footprint ends sit. `open_record.rs`'s `widen` passes exactly
    /// such an end as a 1-based start, so a shim that refused it would fail the walk
    /// at every record that grows to the contig's edge.
    #[test]
    fn a_zero_length_window_is_served_rather_than_refused() {
        let fetcher = RefSeqFetcher(reference());
        assert_eq!(fetcher.fetch(0, 5, 0).expect("in range"), b"");
        assert_eq!(
            fetcher.fetch(0, 11, 0).expect("one past the last base"),
            b""
        );
        // But one *base* past the end is still refused.
        assert!(fetcher.fetch(0, 11, 1).is_err());
    }

    /// **The claim that this shim is semantically empty, checked.** The walker's
    /// contract is canonical uppercase `{A,C,G,T,N}`; a reference that is neither
    /// uppercase nor `ACGTN` comes out as both. Without this, "the two contracts
    /// already agree" would be a doc comment and nothing else.
    ///
    /// It overlaps `ref_seq.rs`'s own canonicalisation test by design and is not a
    /// duplicate of it: what this pins is that **the shim reads the canonical
    /// surface**, so it fails against the one mutation that matters — repointing it
    /// at [`RawRefSeq`](crate::ng::ref_seq::RawRefSeq), which would hand the walker
    /// raw bytes.
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
        // `chrom_name: _` rather than `..`: a bare `..` would also absorb a field
        // added to the variant, and the synthesised name is asserted in its own test.
        assert!(
            matches!(
                error,
                ChromRefFetchError::OutOfBounds {
                    chrom_length: 10,
                    start: 9,
                    end: 14,
                    chrom_name: _,
                }
            ),
            "unexpected: {error:?}"
        );
    }

    /// **The narrowing saturates rather than wraps.** `start + length` can exceed
    /// `u32` even though both arrive as `u32`, and a truncating cast would report an
    /// `end` *below* the `start` it belongs to — a message that lies about the very
    /// window it refused. The doc comment on `narrow` states this; without a test it
    /// would be a claim about an unreachable-looking path that is in fact reachable
    /// through the shim's own signature.
    #[test]
    fn an_out_of_bounds_window_wider_than_u32_saturates_rather_than_wrapping() {
        let fetcher = RefSeqFetcher(reference());
        let error = fetcher
            .fetch(0, u32::MAX, 10)
            .expect_err("far past a 10-base contig");
        let ChromRefFetchError::OutOfBounds {
            chrom_name,
            chrom_length,
            start,
            end,
        } = &error
        else {
            panic!("unexpected: {error:?}");
        };
        assert_eq!(*start, u32::MAX);
        assert_eq!(*end, u32::MAX, "saturated, not wrapped");
        assert!(end >= start, "the reported window must not read backwards");
        assert_eq!(*chrom_length, 10, "the contig's real length");
        assert!(
            chrom_name.contains('0'),
            "the id stands in for the name: {chrom_name}"
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

    /// A broken reference read stays an I/O failure and **keeps its source**. The
    /// walker wraps it in `WalkerError::Fasta`, which renders the source and nothing
    /// else about the cause — so a dropped source is a run that fails without saying
    /// why. This is the arm no in-memory fixture can reach and the one a user will.
    #[test]
    fn a_broken_reference_read_is_reported_as_io_with_its_source_intact() {
        let error = RefSeqFetcher(BrokenRefSeq)
            .fetch(3, 1, 4)
            .expect_err("every fetch on this reference fails");
        let ChromRefFetchError::Io { chrom_name, source } = &error else {
            panic!("unexpected: {error:?}");
        };
        assert_eq!(source.kind(), io::ErrorKind::UnexpectedEof);
        assert_eq!(source.to_string(), "truncated fasta");
        assert!(
            chrom_name.contains('3'),
            "the id stands in for the name: {chrom_name}"
        );
    }

    /// **The shim, composed with the walker it exists for.**
    ///
    /// Every other test here calls `RefSeqFetcher::fetch` directly, and the 113
    /// inherited tests run the walk over `MockFasta` — so without this, the one piece
    /// of code this milestone *authored* is never on a walk. Same reads, same bases,
    /// two fetchers: any difference is the shim's.
    ///
    /// The deletion is not decoration. It puts
    /// [`OpenPileupRecordTable::widen`](super::open_record) on the path, which fetches
    /// `(chrom_id, old_end, extra_len)` — a record's **exclusive** footprint end passed
    /// as a **1-based start**. That is the one call shape where a coordinate-convention
    /// disagreement between ng's reference and the walker would produce wrong REF bases
    /// rather than an error.
    ///
    /// `PileupRecord` has no `PartialEq`, so its `Debug` rendering is the oracle —
    /// exact, and it names the field that moved when it fails.
    #[test]
    fn a_walk_over_the_shim_matches_the_same_walk_over_the_mock_reference() {
        use super::tests::{MockFasta, snp_read};
        use crate::ng::types::ReadGroupId;
        use std::sync::Arc;

        const BASES: &str = "ACGTACGTAC";
        let reads = || {
            vec![
                snp_read("ref", 1, b"ACGTACGTAC", &[30; 10]),
                snp_read("snp", 1, b"ACGTTCGTAC", &[30; 10]),
                // M4 D3 M3: ref 1–4 matched, 5–7 deleted, 8–10 matched. The deletion
                // widens the record at position 4 past the reads that opened it.
                PreparedRead {
                    chrom_id: 0,
                    alignment_start: 1,
                    alignment_end: 10,
                    cigar: vec![CigarOp::Match(4), CigarOp::Deletion(3), CigarOp::Match(3)],
                    seq: b"ACGTTAC".to_vec(),
                    bq_baq: vec![30; 7],
                    mq_log_err: -3.0,
                    mapq: 60,
                    is_reverse_strand: false,
                    qname: Arc::from("del"),
                    mate_role: MateRole::Solo,
                    adaptor_boundary: None,
                    read_group: ReadGroupId(0),
                },
            ]
        };

        let records: Vec<_> = run(reads(), MockFasta::new(BASES), &WalkerConfig::default())
            .map(|record| record.expect("the mock walk succeeds"))
            .collect();

        // The fixture has to *reach* the path this test exists for, or it compares two
        // walks that both did nothing interesting. Asserted on the records themselves
        // rather than on their rendering: a REF allele longer than one base is a record
        // `widen` grew, and growing it is the fetch shape (`old_end` as a 1-based start)
        // that the direct-`fetch` tests above never exercise.
        assert!(!records.is_empty(), "the fixture must emit records");
        let widened = records
            .iter()
            .find(|record| record.alleles[0].seq.len() > 1)
            .expect("the deletion must widen a record, or `widen` is not on the walk");
        assert_eq!(
            widened.alleles.len(),
            3,
            "the widened record should carry REF, the SNP read's haplotype and the deletion"
        );

        let over_mock: Vec<String> = records.iter().map(|record| format!("{record:?}")).collect();
        let over_shim: Vec<String> = run(
            reads(),
            RefSeqFetcher(InMemoryRefSeq::from_contigs(vec![
                BASES.as_bytes().to_vec(),
            ])),
            &WalkerConfig::default(),
        )
        .map(|record| format!("{:?}", record.expect("the shim walk succeeds")))
        .collect();

        assert_eq!(over_shim, over_mock, "the shim must serve the same walk");
    }

    /// The one `DEFAULT_*` the copy had to bring with it, pinned against production's.
    ///
    /// `DEFAULT_MAX_ACTIVE_READS` is declared *inside* `chain_id_allocator.rs`, so the
    /// verbatim copy forked it — where the other constants are reached from production
    /// by name and cannot drift. Production's value is the baseline Milestone B's
    /// differential runs against, so a retune on either side must fail here rather than
    /// silently in the walk.
    #[test]
    fn the_copied_active_reads_cap_is_still_productions() {
        assert_eq!(
            DEFAULT_MAX_ACTIVE_READS,
            crate::pileup::walker::DEFAULT_MAX_ACTIVE_READS,
            "ng's copied chain_id_allocator has drifted from production's cap"
        );
    }
}
