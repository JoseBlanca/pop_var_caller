//! Step 2's v1 implementation — **pass-through + left-alignment**, producing production's
//! `PreparedRead`.
//!
//! A read whose CIGAR carries no indel is built straight from the record: nothing to shift, no
//! reference fetched, no buffer touched. A read that does carry one has its indels rewritten into
//! their leftmost equivalent spelling, so that equivalent variants get an identical one and the
//! reads supporting them pool into a single candidate instead of scattering across several weak
//! ones. Bases and qualities are copied through untouched; only the CIGAR changes.
//!
//! **The shifting is not implemented here.** It is production's `left_align_indels` (itself a port
//! of GATK's `AlignmentUtils.leftAlignIndels`), reached through the
//! [`AlignmentNormalizer`](crate::ng::alignment::AlignmentNormalizer) trait in
//! [`alignment`](crate::ng::alignment) — this module supplies the reference window, the round-trip
//! into and out of an `Alignment`, and the policy of when to fetch at all.
//!
//! Design: `doc/devel/ng/spec/read_preparation.md`, `doc/devel/ng/arch/read_preparation.md`.

use super::{ReadPrepError, ReadPreparer};
use crate::bam::alignment_input::{MappedRead, cigar_ref_span};
use crate::ng::alignment::{Alignment, AlignmentNormalizer, DefaultAlignmentNormalizer};
use crate::ng::ref_seq::RefSeq;
use crate::ng::types::ContigId;
use crate::pileup::per_sample::baq_engine::prepare_passthrough;
use crate::pileup::walker::{CigarOp, PreparedRead};

/// Whether this line-up carries anything left-alignment could move.
///
/// **It must stay exactly production's `is_indel`** — insertion or deletion, nothing else
/// ([indel_norm.rs](../../../pileup/walker/indel_norm.rs)). This predicate is what decides
/// whether a read fetches a reference window at all, so a narrower one here would silently skip
/// normalization for a read production *would* shift: a mis-placed indel, with no error and
/// nothing in the output to show for it.
///
/// Production reaches the same conclusion one layer down — `left_align_indels` early-returns on
/// a CIGAR with no indel before it touches the reference. ng asks first instead, because unlike
/// production it has nothing else to fetch the window *for*: `F1` moved to step 1, so
/// left-alignment is the reference's only consumer here (spec §5).
fn cigar_has_indel(cigar: &[CigarOp]) -> bool {
    cigar
        .iter()
        .any(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_)))
}

/// Reused buffers for [`LeftAlignPreparer`] — allocated once per worker, never per read.
///
/// It holds one thing, and holding it is the whole reason the scratch exists:
/// [`ReadPreparer::prepare_read`](super::ReadPreparer::prepare_read) takes `&self`, and
/// [`RefSeq::fetch_into`] writes the window into a buffer the *caller* owns — so without a
/// scratch every indel-bearing read would allocate one. **Reads with no indel never touch it**:
/// they need no reference at all.
///
/// Buffers only. Nothing here may change a result, which is what makes `Default` a safe way to
/// build it: whatever `Default` picks is invisible at every call site. When the gated re-align
/// mode lands, its alignment matrices join this struct.
#[derive(Debug, Default)]
pub struct LeftAlignScratch {
    /// The reference window for the read being prepared, refilled per indel-bearing read.
    reference_window: Vec<u8>,
}

/// Step 2's v1 preparer: **pass through, or canonicalise the indels** a read already carries.
///
/// # What it is generic over, and why both are visible
///
/// `R` is the reference it fetches windows from. It is held as a field rather than passed per
/// call because a preparer that needs no reference should not have to carry an unused one — the
/// pass-through-only arm holds nothing, and "needs no reference" stays a compile-time fact.
///
/// `N` is the normalizer that does the shifting. The alignment module ships three so they can be
/// *compared*, and names [`DefaultAlignmentNormalizer`] to record which one won — so burying a
/// concrete choice in here would hide it at the use site and make that comparison unrepeatable
/// through this step. It costs nothing: the default normalizer is a unit struct.
///
/// # Per worker, never shared
///
/// One instance belongs to one worker thread. The reference accessors carry state — a resident
/// window and a reader sliding forward through the contig — and sharing that state across threads
/// would mean synchronising a cache on a per-read path, to save an amount of memory that does not
/// matter (arch §4). So the driver builds one of these inside each worker.
///
/// Deliberately **not** `Clone` or `Copy`: no real reference accessor is either — a windowed one
/// holds an open file and a sliding buffer — and a derive that can never fire would still suggest
/// the wrong ownership model, that a preparer is something to hand around rather than something a
/// worker owns.
#[derive(Debug)]
pub struct LeftAlignPreparer<R: RefSeq, N: AlignmentNormalizer = DefaultAlignmentNormalizer> {
    /// The reference this preparer fetches its windows from — **the uppercased, canonical view**
    /// ([`RefSeq::fetch_into`]), not production's raw bytes. An alignment has no business
    /// treating a base differently because a masking tool lowercased it (spec §6).
    reference: R,
    /// The normalizer the canonicalize path calls. Stateless; one value serves every read.
    normalizer: N,
}

impl<R: RefSeq> LeftAlignPreparer<R, DefaultAlignmentNormalizer> {
    /// Build over `reference` with the recommended normalizer — algorithm 1a, the
    /// GATK/production left-aligner.
    ///
    /// This is the constructor a caller wants unless it is running a bake-off; the alternative
    /// normalizers exist as comparators, not as choices to make by taste.
    pub fn with_default_normalizer(reference: R) -> Self {
        Self::new(reference, DefaultAlignmentNormalizer::default())
    }
}

impl<R: RefSeq, N: AlignmentNormalizer> LeftAlignPreparer<R, N> {
    /// Build over `reference`, naming the normalizer explicitly — the bake-off constructor.
    pub fn new(reference: R, normalizer: N) -> Self {
        Self {
            reference,
            normalizer,
        }
    }

    /// The reference this preparer reads from. **Test/diagnostic only** — a driver holds the
    /// reference before it builds the preparer, so nothing in production needs it back, and
    /// eviction needs `&mut` which this cannot give.
    #[cfg(test)]
    pub fn reference(&self) -> &R {
        &self.reference
    }

    /// Rewrite this read's indels into their leftmost equivalent spelling, in place.
    ///
    /// Only the CIGAR changes: bases, qualities and the placement start are all left where they
    /// are. The default normalizer wraps production's `left_align_indels`, which production calls
    /// with end-deletion stripping *off*, so a deletion that slides to the read's start stays a
    /// leading deletion rather than moving `alignment_start` (spec §5).
    fn canonicalize(
        &self,
        read: &mut MappedRead,
        scratch: &mut LeftAlignScratch,
    ) -> Result<(), ReadPrepError> {
        let reference_span = u64::from(cigar_ref_span(&read.cigar));
        // A read that consumes no reference — an all-insertion CIGAR — has nothing to align
        // against, and asking for a zero-length window would be asking a question with no
        // answer. Production guards the same case by skipping F3 entirely when its span is zero
        // (`fetch_raw_slice` returns `None`), so skipping here keeps the two in step.
        if reference_span == 0 {
            return Ok(());
        }

        // `ref_id` indexes the `u32` contig table, so a value that did not fit would be a corrupt
        // record rather than a legal one — step 1 makes the same conversion with the same reason.
        let contig = ContigId(u32::try_from(read.ref_id).expect("ref_id fits u32"));
        // The **uppercased** view, not production's raw bytes: an alignment must not treat a base
        // differently because a masking tool lowercased it. Production compares raw reference
        // bytes against read bases that were uppercased at decode, so its left-alignment does
        // nothing at all on a soft-masked reference (spec §6).
        self.reference
            .fetch_into(
                contig,
                read.pos,
                reference_span,
                &mut scratch.reference_window,
            )
            .map_err(ReadPrepError::Reference)?;
        // The window must cover the read's whole footprint, and `fetch_into` guarantees it: on
        // success it writes exactly `length` bases, on failure it returns `Err`. Worth asserting
        // because the failure mode of a short window is silent — the normalizer fails *safe*,
        // handing back the CIGAR untouched, so normalization would simply not happen and nothing
        // would say so.
        debug_assert_eq!(
            scratch.reference_window.len() as u64,
            reference_span,
            "a successful fetch must fill the whole window"
        );

        // The normalizer works on an `Alignment`, not on a bare CIGAR, so the CIGAR moves out of
        // the read and back — no clone. `reference_offset` is 0 because the window *starts* at
        // the read's first aligned base; the offset exists for callers holding a larger stretch,
        // and a non-zero value here would place the read against the wrong bases.
        let mut alignment = Alignment {
            reference_offset: 0,
            cigar: std::mem::take(&mut read.cigar),
        };
        self.normalizer
            .normalize(&mut alignment, &read.seq, &scratch.reference_window);
        read.cigar = alignment.cigar;
        Ok(())
    }
}

impl<R: RefSeq, N: AlignmentNormalizer> ReadPreparer for LeftAlignPreparer<R, N> {
    type Scratch = LeftAlignScratch;

    /// Prepare one read: canonicalise its indels if it has any, otherwise hand back what the
    /// record already says.
    ///
    /// **A read with no indel costs one CIGAR scan.** No reference is fetched, the scratch is
    /// never touched, and nothing can fail — which is most reads. That is not an optimisation
    /// bolted on afterwards; it is what makes the reference's *only* consumer here the
    /// left-alignment of the minority of reads that carry a gap (spec §5).
    fn prepare_read(
        &self,
        mut read: MappedRead,
        scratch: &mut LeftAlignScratch,
    ) -> Result<Option<PreparedRead>, ReadPrepError> {
        if cigar_has_indel(&read.cigar) {
            self.canonicalize(&mut read, scratch)?;
        }
        Ok(Some(into_prepared(read)))
    }
}

/// Build the `PreparedRead` from a record whose CIGAR is final.
///
/// The whole field wiring — `alignment_end` from the CIGAR, `mate_role` and strand from the flag,
/// `mq_log_err` from MAPQ, `qname`, the adaptor boundary, and the raw qualities copied into
/// `bq_baq` uncapped — is production's `prepare_passthrough`, its `--no-baq` path, called rather
/// than re-derived. That is what makes the parity fixture exact on every field this step does not
/// itself compute (spec §9, §11).
fn into_prepared(read: MappedRead) -> PreparedRead {
    // `ref_id` indexes the `u32` contig table, so a value that did not fit would be a corrupt
    // record rather than a legal one — the same conversion, with the same reasoning, that step 1
    // and production's own passthrough arm make.
    let chrom_id = u32::try_from(read.ref_id).expect("ref_id fits u32");
    prepare_passthrough(read, chrom_id)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::alignment::left_align_repeated::RepeatedLeftAligner;
    use crate::ng::ref_seq::{InMemoryRefSeq, RefSeqError};

    fn reference() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![b"ACGTACGTAC".to_vec()])
    }

    /// A reference that **fails the test if anything fetches from it**.
    ///
    /// This is how "a read with no indel never touches the reference" is pinned as behaviour
    /// rather than left as a claim in a doc comment: the assertion is that the fetch is not
    /// reached at all, which no amount of inspecting the returned read could show.
    struct NeverFetched;
    impl RefSeq for NeverFetched {
        fn fetch_into(
            &self,
            _contig: ContigId,
            _start_1based: u64,
            _length: u64,
            _dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            panic!("a read with no indel must not fetch a reference window");
        }
    }

    /// A reference whose every fetch fails — for proving a broken fetch ends the run.
    struct AlwaysFailsToFetch;
    impl RefSeq for AlwaysFailsToFetch {
        fn fetch_into(
            &self,
            contig: ContigId,
            _start_1based: u64,
            _length: u64,
            _dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            Err(RefSeqError::UnknownContig(contig))
        }
    }

    /// A read at 1-based `pos` with `cigar` over `seq`.
    fn read_with(pos: u64, cigar: Vec<CigarOp>, seq: &[u8]) -> MappedRead {
        MappedRead {
            qname: b"read1".to_vec(),
            flag: 0,
            ref_id: 0,
            pos,
            mapq: 60,
            cigar,
            seq: seq.to_vec(),
            qual: vec![30; seq.len()],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            source_file_index: 0,
        }
    }

    /// The scratch starts empty and allocates nothing: a preparer that never meets an
    /// indel-bearing read never grows it.
    #[test]
    fn a_fresh_scratch_holds_no_window() {
        let scratch = LeftAlignScratch::default();
        assert!(scratch.reference_window.is_empty());
        assert_eq!(
            scratch.reference_window.capacity(),
            0,
            "Default must not pre-allocate — it decides nothing"
        );
    }

    /// The default constructor reaches the recommended normalizer, and the preparer keeps the
    /// reference it was handed.
    #[test]
    fn the_default_constructor_holds_its_reference() {
        let preparer = LeftAlignPreparer::with_default_normalizer(reference());
        let mut window = Vec::new();
        preparer
            .reference()
            .fetch_into(ContigId(0), 1, 4, &mut window)
            .expect("the fixture contig is 10 bases long");
        assert_eq!(window, b"ACGT");
    }

    /// A second normalizer can be named at construction — the bake-off surface the type
    /// parameter exists for. If this stops compiling, the comparison has been foreclosed.
    #[test]
    fn an_alternative_normalizer_can_be_named() {
        let preparer = LeftAlignPreparer::new(reference(), RepeatedLeftAligner::new());
        let mut window = Vec::new();
        preparer
            .reference()
            .fetch_into(ContigId(0), 5, 4, &mut window)
            .expect("the fixture contig is 10 bases long");
        assert_eq!(window, b"ACGT");
    }

    #[test]
    fn the_indel_predicate_is_exactly_production_s() {
        assert!(cigar_has_indel(&[CigarOp::Match(4), CigarOp::Insertion(2)]));
        assert!(cigar_has_indel(&[CigarOp::Match(4), CigarOp::Deletion(2)]));
        assert!(!cigar_has_indel(&[CigarOp::Match(10)]));
        // Everything else a CIGAR can carry is *not* an indel — a narrower predicate here would
        // skip normalization for a read production would shift.
        assert!(!cigar_has_indel(&[
            CigarOp::SoftClip(3),
            CigarOp::Match(4),
            CigarOp::Skip(5),
            CigarOp::HardClip(2),
            CigarOp::Padding(1),
            CigarOp::SeqMatch(4),
            CigarOp::SeqMismatch(1),
        ]));
    }

    /// **The conditional fetch, pinned.** A read with no indel is built straight from the record:
    /// the reference is never consulted, which is why this passes a reference that panics.
    #[test]
    fn a_read_with_no_indel_never_fetches_a_reference_window() {
        let preparer = LeftAlignPreparer::with_default_normalizer(NeverFetched);
        let mut scratch = LeftAlignScratch::default();
        let prepared = preparer
            .prepare_read(read_with(3, vec![CigarOp::Match(4)], b"ACGT"), &mut scratch)
            .expect("no fetch, so nothing can fail")
            .expect("v1 never declines a read");

        assert_eq!(prepared.cigar, vec![CigarOp::Match(4)]);
        assert_eq!(prepared.seq, b"ACGT");
        assert_eq!(prepared.alignment_start, 3);
        assert!(
            scratch.reference_window.is_empty(),
            "the scratch is not touched either"
        );
    }

    /// Qualities pass through **uncapped** — production's `--no-baq` behaviour, and all of v1's
    /// quality handling, since BAQ is deferred.
    #[test]
    fn base_qualities_pass_through_uncapped() {
        let preparer = LeftAlignPreparer::with_default_normalizer(NeverFetched);
        let mut scratch = LeftAlignScratch::default();
        let mut read = read_with(1, vec![CigarOp::Match(4)], b"ACGT");
        read.qual = vec![2, 40, 11, 37];

        let prepared = preparer
            .prepare_read(read, &mut scratch)
            .expect("no fetch")
            .expect("prepared");
        assert_eq!(prepared.bq_baq, vec![2, 40, 11, 37]);
    }

    /// **The transform itself.** A deletion sitting inside a homopolymer run can be written at
    /// several equivalent positions; left-alignment picks the leftmost, so that two reads showing
    /// the same event spell it the same way and pool into one candidate.
    ///
    /// Reference `TAAAAG`, read `TAAAG` — one `A` deleted from the run. The mapper spelled the
    /// deletion at the run's right end (`M4 D1 M1`); leftmost is its left end (`M1 D1 M4`).
    #[test]
    fn an_indel_in_a_homopolymer_is_moved_to_its_leftmost_spelling() {
        let reference = InMemoryRefSeq::from_contigs(vec![b"TAAAAG".to_vec()]);
        let preparer = LeftAlignPreparer::with_default_normalizer(reference);
        let mut scratch = LeftAlignScratch::default();

        let prepared = preparer
            .prepare_read(
                read_with(
                    1,
                    vec![CigarOp::Match(4), CigarOp::Deletion(1), CigarOp::Match(1)],
                    b"TAAAG",
                ),
                &mut scratch,
            )
            .expect("the window is in range")
            .expect("prepared");

        assert_eq!(
            prepared.cigar,
            vec![CigarOp::Match(1), CigarOp::Deletion(1), CigarOp::Match(4)],
            "the deletion should sit at the leftmost equivalent position"
        );
        // The placement start does not move: production strips no end deletions, so a shifted
        // indel never changes where the read begins.
        assert_eq!(prepared.alignment_start, 1);
        assert_eq!(prepared.seq, b"TAAAG", "bases are untouched");
    }

    /// A failed fetch ends the run — it is never folded into `Ok(None)`, which would look like a
    /// read that simply produced nothing.
    #[test]
    fn a_failed_fetch_is_an_error_not_a_decline() {
        let preparer = LeftAlignPreparer::with_default_normalizer(AlwaysFailsToFetch);
        let mut scratch = LeftAlignScratch::default();
        let outcome = preparer.prepare_read(
            read_with(
                1,
                vec![CigarOp::Match(2), CigarOp::Deletion(1), CigarOp::Match(2)],
                b"ACGT",
            ),
            &mut scratch,
        );
        assert!(matches!(
            outcome,
            Err(ReadPrepError::Reference(RefSeqError::UnknownContig(_)))
        ));
    }

    /// A read that consumes no reference has nothing to align against. Production skips `F3`
    /// for it; so does this, rather than asking for a zero-length window.
    #[test]
    fn an_all_insertion_read_is_passed_through_without_fetching() {
        let preparer = LeftAlignPreparer::with_default_normalizer(NeverFetched);
        let mut scratch = LeftAlignScratch::default();
        let prepared = preparer
            .prepare_read(
                read_with(1, vec![CigarOp::Insertion(4)], b"ACGT"),
                &mut scratch,
            )
            .expect("no fetch is attempted")
            .expect("prepared");
        assert_eq!(prepared.cigar, vec![CigarOp::Insertion(4)]);
    }
}
