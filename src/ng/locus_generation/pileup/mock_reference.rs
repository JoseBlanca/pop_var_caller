//! **ng's, not a copy** — production's `MockFasta` served through ng's
//! [`RefSeq`], so the copied test suite can drive ng's walker (A0).
//!
//! # Why this exists, and why it is not the shim coming back
//!
//! A0 deletes `RefSeqFetcher`, the newtype that presented **ng's reference** as
//! **production's `MultiChromRefFetcher`** so a verbatim copy of production's
//! walker would compile. That direction is gone: ng's walk is bound by
//! [`RefSeq`] and imports neither of production's fetcher types.
//!
//! This is the *other* direction, and it is a fixture concern rather than an
//! architectural one. `MockFasta` — production's in-memory test double — lives
//! inside [`tests.rs`](super::tests), which is one of the files
//! [`copy_fidelity`](super::copy_fidelity) still guards as verbatim, and its
//! `chromosomes` field is private to that module. So ng cannot make it a
//! [`RefSeq`] *in place* without editing a file the plan keeps frozen, and cannot
//! read its bytes from outside either. What it can do is call the one surface
//! `MockFasta` exposes and re-present the answer, which is what this file is.
//!
//! It dies with `tests.rs`, when the copied suite is adapted to ng's own locus
//! type (arch, *Module home*).
//!
//! # The one semantic decision, stated rather than left implicit
//!
//! `MockFasta::fetch` is a **raw passthrough** — it serves whatever bytes the
//! fixture string held. [`RefSeq`], by contract, serves canonical uppercase
//! `{A,C,G,T,N}`. So this impl canonicalises, using the very function every
//! file-backed `RefSeq` uses ([`canonicalise`]), rather than honouring the
//! double's laxity and handing ng's walk bytes its trait says cannot arrive.
//!
//! On every inherited fixture the two coincide: the copied suite's references are
//! uppercase `ACGT` throughout, so canonicalisation is the identity and the walk
//! sees exactly what production's walk sees. That coincidence is *pinned*, not
//! assumed — see `parity.rs`'s
//! `both_sides_of_the_differential_are_served_the_same_bytes`, which fails if a
//! fixture ever grows a lower-case or ambiguity-coded base. Without that pin such
//! a fixture would surface as "the two walkers disagree" and be chased into the
//! walk.

use crate::fasta::MultiChromRefFetcher;
use crate::fasta::fetcher::{ChromRefFetchError, canonicalise};
use crate::ng::ref_seq::{RefSeq, RefSeqError};
use crate::ng::types::ContigId;

use super::tests::MockFasta;

impl RefSeq for MockFasta {
    /// Delegate to the double's own `fetch` and canonicalise on the way out.
    ///
    /// The **window arithmetic is the double's**, deliberately: production's walk
    /// and ng's then agree on which fetches succeed, not merely on the bytes the
    /// successful ones return. Reimplementing the bounds check here would make
    /// the two sides of the differential disagree at the edges for a reason that
    /// belongs to neither walker.
    ///
    /// `u32` narrowing rather than a `try_from`: the walker's own coordinates are
    /// `u32` (that is the width the copies keep, arch §3), so a request that
    /// cannot be expressed as one cannot have come from the walk. The saturation
    /// turns such a request into an out-of-bounds against a fixture contig of at
    /// most a few hundred bases, which is what it is.
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        let narrow = |value: u64| u32::try_from(value).unwrap_or(u32::MAX);
        let raw =
            MultiChromRefFetcher::fetch(self, contig.get(), narrow(start_1based), narrow(length))
                .map_err(|error| to_ref_seq_error(contig, start_1based, length, error))?;
        // Only on success, so `dst` is left untouched on error — `RefSeq`'s
        // contract, and the property `fetch_into_leaves_dst_unchanged_on_error`
        // pins for the real implementations.
        dst.clear();
        dst.extend(raw.into_iter().map(canonicalise));
        Ok(())
    }
}

/// The same, through a shared reference.
///
/// `MultiChromRefFetcher` carries a blanket `impl … for &T` and [`RefSeq`] does not, so
/// the copied suite's three `run(reads, &fa, …)` call sites would otherwise not compile —
/// and `tests.rs` is still one of the verbatim copies, so they cannot be changed.
///
/// **The general form belongs in `ref_seq.rs`, beside the `Arc<T>` forward it mirrors**,
/// and is deliberately not written here: that would widen a shared ng trait's surface
/// from inside a step whose subject is the pileup's own reference plumbing. Raised at
/// Checkpoint A instead.
impl RefSeq for &MockFasta {
    fn fetch_into(
        &self,
        contig: ContigId,
        start_1based: u64,
        length: u64,
        dst: &mut Vec<u8>,
    ) -> Result<(), RefSeqError> {
        RefSeq::fetch_into(*self, contig, start_1based, length, dst)
    }
}

/// Translate the double's [`ChromRefFetchError`] into ng's [`RefSeqError`],
/// variant for variant.
///
/// Exhaustive on the **matched** side, so a variant added to production's enum
/// stops this compiling rather than falling into a catch-all that reports the
/// wrong failure mode. `OutOfPattern` belongs to the streaming fetcher and
/// `MockFasta` never mints it; it is mapped rather than ignored so the match can
/// stay exhaustive without a wildcard.
///
/// `chrom_name` is discarded: it is a string the double synthesises from the id
/// this function already has, and `RefSeqError` names contigs by [`ContigId`].
/// The requested window is taken from the **caller's** `u64` arguments rather
/// than from the narrowed ones the double echoes back, so a request too wide for
/// `u32` is reported as the window that was actually asked for.
fn to_ref_seq_error(
    contig: ContigId,
    start_1based: u64,
    length: u64,
    error: ChromRefFetchError,
) -> RefSeqError {
    match error {
        ChromRefFetchError::OutOfBounds { chrom_length, .. } => RefSeqError::OutOfBounds {
            contig,
            contig_length: u64::from(chrom_length),
            start: start_1based,
            end: start_1based.saturating_add(length),
        },
        ChromRefFetchError::InvalidStart => RefSeqError::InvalidStart,
        // The double reports an unknown chromosome as `Io { NotFound }`; ng has a
        // variant that says exactly that, so the round trip is lossless rather
        // than folded into I/O. Any other `Io` is a genuine read failure.
        ChromRefFetchError::Io { source, .. } if source.kind() == std::io::ErrorKind::NotFound => {
            RefSeqError::UnknownContig(contig)
        }
        ChromRefFetchError::Io { source, .. } => RefSeqError::Io { contig, source },
        ChromRefFetchError::OutOfPattern { .. } => RefSeqError::Io {
            contig,
            source: std::io::Error::other("MockFasta cannot report OutOfPattern"),
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::ref_seq::InMemoryRefSeq;

    /// Two contigs, the second one written the way a real assembly writes repeats
    /// and ambiguity: lower case, and a code that is neither `ACGT` nor `N`.
    fn reference() -> MockFasta {
        MockFasta::with_chromosomes(&["ACGTACGTAC", "acgtRYacgt"])
    }

    #[test]
    fn it_serves_the_window_it_was_asked_for() {
        let reference = reference();
        assert_eq!(
            RefSeq::fetch(&reference, ContigId(0), 1, 4).expect("in range"),
            b"ACGT"
        );
        // 1-based and inclusive of `start`: base 5 is the fifth, not the sixth.
        assert_eq!(
            RefSeq::fetch(&reference, ContigId(0), 5, 3).expect("in range"),
            b"ACG"
        );
        assert_eq!(
            RefSeq::fetch(&reference, ContigId(0), 10, 1).expect("in range"),
            b"C"
        );
    }

    /// A zero-length window is a legal question with an empty answer — including
    /// at the position one past the contig's last base, which is where the
    /// walker's **exclusive** footprint ends sit. `open_record.rs`'s `widen`
    /// passes exactly such an end as a 1-based start, so a reference that refused
    /// it would fail the walk at every record that grows to the contig's edge.
    #[test]
    fn a_zero_length_window_is_served_rather_than_refused() {
        let reference = reference();
        assert_eq!(
            RefSeq::fetch(&reference, ContigId(0), 5, 0).expect("in range"),
            b""
        );
        assert_eq!(
            RefSeq::fetch(&reference, ContigId(0), 11, 0).expect("one past the last base"),
            b""
        );
        // But one *base* past the end is still refused.
        assert!(RefSeq::fetch(&reference, ContigId(0), 11, 1).is_err());
    }

    /// **The canonicalisation, checked.** The double passes its bytes through;
    /// `RefSeq` promises uppercase `{A,C,G,T,N}`. Without this the promise would
    /// be a doc comment and nothing else — and the mutation it fails against is
    /// dropping the `.map(canonicalise)`.
    #[test]
    fn it_canonicalises_what_the_double_passes_through() {
        let reference = reference();
        let bases = RefSeq::fetch(&reference, ContigId(1), 1, 10).expect("in range");
        assert_eq!(bases, b"ACGTNNACGT", "lower case up, non-ACGT folded to N");
        // The double itself does not: this is what says the fold is ours.
        assert_eq!(
            MultiChromRefFetcher::fetch(&reference, 1, 1, 10).expect("in range"),
            b"acgtRYacgt",
        );
    }

    /// A window past the contig's end is refused as ng's own out-of-bounds — not
    /// folded into an I/O failure, which is what the walker's error routing
    /// distinguishes. The reported length is the contig's real one.
    #[test]
    fn a_window_past_the_contig_end_is_out_of_bounds() {
        let reference = reference();
        let error = RefSeq::fetch(&reference, ContigId(0), 9, 5).expect_err("only 10 bases exist");
        assert!(
            matches!(
                error,
                RefSeqError::OutOfBounds {
                    contig: ContigId(0),
                    contig_length: 10,
                    start: 9,
                    end: 14,
                }
            ),
            "unexpected: {error:?}"
        );
    }

    /// The 1-based contract, carried across unchanged: `start_1based == 0` is its
    /// own variant on both sides, never an out-of-bounds.
    #[test]
    fn a_zero_start_keeps_its_own_variant() {
        let reference = reference();
        let error = RefSeq::fetch(&reference, ContigId(0), 0, 4).expect_err("1-based coordinates");
        assert!(
            matches!(error, RefSeqError::InvalidStart),
            "unexpected: {error:?}"
        );
    }

    /// An unknown contig reaches ng's own variant rather than being folded into
    /// I/O — the double reports it as `Io { NotFound }`, and the translation
    /// recovers the meaning.
    #[test]
    fn an_unknown_contig_keeps_its_own_variant() {
        let reference = reference();
        let error =
            RefSeq::fetch(&reference, ContigId(7), 1, 4).expect_err("only two contigs exist");
        assert!(
            matches!(error, RefSeqError::UnknownContig(ContigId(7))),
            "unexpected: {error:?}"
        );
    }

    /// **A start too large for the walker's own `u32` is reported at its real
    /// value, not at the saturated one.** The narrowing saturates on the way in,
    /// so the double refuses the window; the error is then built from the
    /// caller's `u64` arguments, and a report built from the narrowed ones would
    /// name a start of `u32::MAX` for a window nobody asked about.
    ///
    /// The length is `0` because the double adds `start + length` in `u32` and a
    /// debug build panics on the overflow before this impl's translation is
    /// reached — `tests.rs` is still verbatim, so that arithmetic is not ours to
    /// change. A zero-length window at a start past the contig is refused all the
    /// same, which is the path under test.
    #[test]
    fn a_start_wider_than_u32_is_reported_at_its_real_value() {
        let reference = reference();
        let error = RefSeq::fetch(&reference, ContigId(0), u64::from(u32::MAX) + 5, 0)
            .expect_err("far past a 10-base contig");
        let RefSeqError::OutOfBounds {
            contig_length,
            start,
            end,
            ..
        } = error
        else {
            panic!("unexpected error variant");
        };
        assert_eq!(start, u64::from(u32::MAX) + 5, "not narrowed in the report");
        assert_eq!(
            end,
            u64::from(u32::MAX) + 5,
            "start + length, both un-narrowed"
        );
        assert!(end >= start, "the reported window must not read backwards");
        assert_eq!(contig_length, 10, "the contig's real length");
    }

    /// **The double, driven through the walk it exists to feed.**
    ///
    /// Every other test here calls `fetch` directly, so without this the walker
    /// never sees this impl at a coordinate it computed itself. Same reads, same
    /// bases, two `RefSeq` implementations — the double and ng's real in-memory
    /// one — so any difference is this file's.
    ///
    /// The deletion is not decoration. It puts
    /// [`OpenPileupRecordTable::widen`](super::super::open_record) on the path,
    /// which fetches `(chrom_id, old_end, extra_len)` — a record's **exclusive**
    /// footprint end passed as a **1-based start**. That is the one call shape
    /// where a coordinate-convention disagreement produces wrong REF bases rather
    /// than an error.
    ///
    /// Compared through the `Debug` rendering rather than by value: it names the field
    /// that moved when it fails. (`PileupRecord` *does* have a hand-written `PartialEq`
    /// — [pileup_record.rs:208](crate::pileup_record) — which an earlier version of this
    /// comment denied; `parity.rs` relies on it. The `Debug` comparison is a
    /// diagnostics choice here, not a necessity.)
    #[test]
    fn a_walk_over_the_double_matches_the_same_walk_over_a_real_reference() {
        use super::super::tests::snp_read;
        use super::super::{CigarOp, MateRole, PreparedRead, WalkerConfig, run};
        use crate::ng::types::ReadGroupId;
        use std::sync::Arc;

        const BASES: &str = "ACGTACGTAC";
        let reads = || {
            vec![
                snp_read("ref", 1, b"ACGTACGTAC", &[30; 10]),
                snp_read("snp", 1, b"ACGTTCGTAC", &[30; 10]),
                // M4 D3 M3: ref 1–4 matched, 5–7 deleted, 8–10 matched. The
                // deletion widens the record at position 4 past the reads that
                // opened it.
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
            .map(|record| record.expect("the walk over the double succeeds"))
            .collect();

        // The fixture has to *reach* the path this test exists for, or it compares
        // two walks that both did nothing interesting. Asserted on the records
        // themselves rather than on their rendering: a REF allele longer than one
        // base is a record `widen` grew, and growing it is the fetch shape
        // (`old_end` as a 1-based start) the direct-`fetch` tests never exercise.
        assert!(!records.is_empty(), "the fixture must emit records");
        // The reach check reads ng's own type (D1): the locus's reference bytes are a field,
        // and a row per read that folded. It used to detour through the `PileupRecord`
        // back-projection to say "REF is `alleles[0]`" — a statement about production's shape
        // that this walk has not produced since B2.
        let widened = records
            .iter()
            .find(|locus| locus.reference_bases.len() > 1)
            .expect("the deletion must widen a record, or `widen` is not on the walk");
        assert_eq!(
            widened.observations.len(),
            3,
            "the widened locus should carry a row each for the reference-matching read, the \
             SNP read's haplotype and the deletion"
        );

        let over_double: Vec<String> = records.iter().map(|record| format!("{record:?}")).collect();
        let over_real: Vec<String> = run(
            reads(),
            InMemoryRefSeq::from_contigs(vec![BASES.as_bytes().to_vec()]),
            &WalkerConfig::default(),
        )
        .map(|record| {
            format!(
                "{:?}",
                record.expect("the walk over the real reference succeeds")
            )
        })
        .collect();

        assert_eq!(
            over_real, over_double,
            "the double must serve the same walk as a real reference"
        );
    }
}
