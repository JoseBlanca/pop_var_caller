//! Step 2's v1 implementation — **pass-through + left-alignment**, producing ng's
//! [`PreparedRead`].
//!
//! The read it hands back is production's `prepare_passthrough` output plus the read group —
//! ng owns the type precisely because production's has nowhere to put that field
//! (`doc/devel/ng/spec/locus_generation_pileup.md` §6).
//!
//! A read whose CIGAR carries no indel is built straight from the record: nothing to shift, no
//! reference fetched, no buffer touched. A read that does carry one has its indels rewritten into
//! their leftmost equivalent spelling, so that equivalent variants get an identical one and the
//! reads supporting them pool into a single candidate instead of scattering across several weak
//! ones. Bases and qualities are copied through untouched; only the CIGAR changes.
//!
//! **The shifting is not implemented here.** It is
//! [`left_align_indels`], a port of GATK's
//! `AlignmentUtils.leftAlignIndels`, called directly — this module supplies the reference window
//! and the policy of when to fetch at all.
//!
//! Design: `doc/devel/ng/spec/read_preparation.md`, `doc/devel/ng/arch/read_preparation.md`.

use super::{ReadPrepError, ReadPreparer};
use crate::alignment::indel_norm::left_align_indels;
use crate::bam::alignment_input::CigarOp;
use crate::bam::alignment_input::cigar_ref_span;
use crate::read::aligned_read::AlignedRead;
use crate::read::prepared_read::PreparedRead;
use crate::read::prepared_read::prepare_passthrough;
use crate::ref_seq::{EvictableRefSeq, RefSeq};
use crate::types::ContigId;

/// Whether this line-up carries anything left-alignment could move.
///
/// **It must stay exactly `is_indel` in `alignment/indel_norm.rs`** — insertion or deletion,
/// nothing else. This predicate is what decides whether a read fetches a reference window at
/// all, so a narrower one here would silently skip normalization for a read `left_align_indels`
/// *would* shift: a mis-placed indel, with no error and nothing in the output to show for it.
///
/// The aligner reaches the same conclusion one layer down — `left_align_indels` early-returns on
/// a CIGAR with no indel before it touches the reference. ng asks first instead, because unlike
/// production it has nothing else to fetch the window *for*: the mismatch-fraction filter moved
/// to step 1, so
/// left-alignment is the reference's only consumer here (spec §5).
fn cigar_has_indel(cigar: &[CigarOp]) -> bool {
    cigar
        .iter()
        .any(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_)))
}

/// Reused buffers for [`LeftAlignPreparer`] — allocated once per worker, never per read.
///
/// It holds one thing, and holding it is the whole reason the scratch exists:
/// [`ReadPreparer::prepare_read`] takes `&self`, and
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
/// `R` is the reference it fetches windows from. It is held as a field rather than passed per
/// call because a preparer that needs no reference should not have to carry an unused one — the
/// pass-through-only arm holds nothing, and "needs no reference" stays a compile-time fact.
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
pub struct LeftAlignPreparer<R: RefSeq> {
    /// The reference this preparer fetches its windows from — **the uppercased, canonical view**
    /// ([`RefSeq::fetch_into`]), not production's raw bytes. An alignment has no business
    /// treating a base differently because a masking tool lowercased it (spec §6).
    reference: R,
}

impl<R: RefSeq> LeftAlignPreparer<R> {
    /// Build over `reference`.
    pub fn new(reference: R) -> Self {
        Self { reference }
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
    /// are. [`left_align_indels`] runs with end-deletion stripping *off*, so a deletion that
    /// slides to the read's start stays a leading deletion rather than moving `alignment_start`
    /// (spec §5).
    fn canonicalize(
        &self,
        read: &mut AlignedRead,
        scratch: &mut LeftAlignScratch,
    ) -> Result<(), ReadPrepError> {
        let reference_span = u64::from(cigar_ref_span(&read.cigar));
        // A read that consumes no reference — an all-insertion CIGAR — has nothing to align
        // against, and asking for a zero-length window would be asking a question with no
        // answer. Production skips left-alignment entirely when its span is zero
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
        // because the failure mode of a short window is silent — `left_align_indels` fails
        // *safe*, handing back the CIGAR untouched, so normalization would simply not happen and
        // nothing would say so.
        debug_assert_eq!(
            scratch.reference_window.len() as u64,
            reference_span,
            "a successful fetch must fill the whole window"
        );

        // The window *starts* at the read's first aligned base, which is the reference
        // `left_align_indels` expects (its `ref_seq[0]`).
        left_align_indels(&mut read.cigar, &read.seq, &scratch.reference_window);
        Ok(())
    }
}

impl<R: RefSeq + EvictableRefSeq> ReadPreparer for LeftAlignPreparer<R> {
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
        mut read: AlignedRead,
        scratch: &mut LeftAlignScratch,
    ) -> Result<Option<PreparedRead>, ReadPrepError> {
        if cigar_has_indel(&read.cigar) {
            self.canonicalize(&mut read, scratch)?;
        }
        // The rest of the build — `alignment_end`, mate role, strand, `mq_log_err`, `qname`, the
        // raw qualities as `bq_baq` and the read group — is production's `--no-baq` arm, called
        // rather than re-derived, which is what makes the parity fixture exact on every field this
        // step does not itself compute (spec §9, §11).
        Ok(Some(prepare_passthrough(read)))
    }

    /// Release the reference bases already left behind — see the trait method.
    ///
    /// This preparer fetches only for reads that carry an indel, so its window is sparse and
    /// grows more slowly than the walk's. It still grows: two indels a kilobase apart leave
    /// the kilobase between them resident.
    fn evict_reference_before(&self, pos: u64) {
        self.reference.evict_before(pos);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::PLACEHOLDER_READ_GROUP;
    use crate::ref_seq::{InMemoryRefSeq, RefSeqError};
    use crate::types::ReadGroupId;

    fn reference() -> InMemoryRefSeq {
        InMemoryRefSeq::from_contigs(vec![b"ACGTACGTAC".to_vec()])
    }

    /// A reference that **fails the test if anything fetches from it**.
    ///
    /// This is how "a read with no indel never touches the reference" is pinned as behaviour
    /// rather than left as a claim in a doc comment: the assertion is that the fetch is not
    /// reached at all, which no amount of inspecting the returned read could show.
    struct NeverFetched;
    impl EvictableRefSeq for NeverFetched {
        /// Holds no window, so there is nothing to release.
        fn evict_before(&self, _pos: u64) {}
    }

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
    impl EvictableRefSeq for AlwaysFailsToFetch {
        /// Holds no window, so there is nothing to release.
        fn evict_before(&self, _pos: u64) {}
    }

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
    fn read_with(pos: u64, cigar: Vec<CigarOp>, seq: &[u8]) -> AlignedRead {
        AlignedRead {
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
            read_group: PLACEHOLDER_READ_GROUP,
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

    /// The preparer keeps the reference it was handed.
    #[test]
    fn the_constructor_holds_its_reference() {
        let preparer = LeftAlignPreparer::new(reference());
        let mut window = Vec::new();
        preparer
            .reference()
            .fetch_into(ContigId(0), 1, 4, &mut window)
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
        let preparer = LeftAlignPreparer::new(NeverFetched);
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
        let preparer = LeftAlignPreparer::new(NeverFetched);
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
        let preparer = LeftAlignPreparer::new(reference);
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
        let preparer = LeftAlignPreparer::new(AlwaysFailsToFetch);
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

    /// **The read group survives preparation — on both paths.**
    ///
    /// This is the whole point of ng owning the read type: production's `prepare_passthrough`
    /// drops the group, so a preparer that simply returned its output would leave the pileup with
    /// nothing to key an observation's group cell on. Both arms are driven because they build the
    /// read through different code — the indel arm rewrites the CIGAR first — and **two different
    /// non-zero ids** are used: a non-zero id rules out a defaulted field (`ReadGroupId(0)` reads
    /// the same either way), and two different ones additionally rule out a preparer that latched
    /// the first group it saw.
    #[test]
    fn the_read_group_rides_through_both_paths() {
        let preparer =
            LeftAlignPreparer::new(InMemoryRefSeq::from_contigs(vec![b"TAAAAG".to_vec()]));
        let mut scratch = LeftAlignScratch::default();

        let mut no_indel = read_with(1, vec![CigarOp::Match(5)], b"TAAAA");
        no_indel.read_group = ReadGroupId(7);
        let prepared = preparer
            .prepare_read(no_indel, &mut scratch)
            .expect("in range")
            .expect("prepared");
        assert_eq!(prepared.read_group, ReadGroupId(7), "pass-through path");

        let mut with_indel = read_with(
            1,
            vec![CigarOp::Match(4), CigarOp::Deletion(1), CigarOp::Match(1)],
            b"TAAAG",
        );
        with_indel.read_group = ReadGroupId(9);
        let prepared = preparer
            .prepare_read(with_indel, &mut scratch)
            .expect("in range")
            .expect("prepared");
        assert_eq!(
            prepared.cigar,
            vec![CigarOp::Match(1), CigarOp::Deletion(1), CigarOp::Match(4)],
            "the fixture must actually take the canonicalize arm"
        );
        assert_eq!(prepared.read_group, ReadGroupId(9), "canonicalize path");
    }

    /// Each read keeps **its own** group across one preparer and one scratch.
    ///
    /// The scratch is where a latching bug would be introduced — it is the one piece of state that
    /// survives between reads — and a fixture using one id twice cannot see it. Distinct non-zero
    /// ids, compared as a sequence, so a preparer that stamped the first group it saw on every
    /// subsequent read fails here.
    #[test]
    fn two_reads_from_different_groups_keep_their_own() {
        let preparer = LeftAlignPreparer::new(reference());
        let mut scratch = LeftAlignScratch::default();
        let groups: Vec<ReadGroupId> = [ReadGroupId(3), ReadGroupId(11)]
            .into_iter()
            .map(|group| {
                let mut read = read_with(1, vec![CigarOp::Match(4)], b"ACGT");
                read.read_group = group;
                preparer
                    .prepare_read(read, &mut scratch)
                    .expect("no fetch")
                    .expect("prepared")
                    .read_group
            })
            .collect();
        assert_eq!(groups, vec![ReadGroupId(3), ReadGroupId(11)]);
    }

    /// A read that consumes no reference has nothing to align against. Production skips
    /// left-alignment for it; so does this, rather than asking for a zero-length window.
    #[test]
    fn an_all_insertion_read_is_passed_through_without_fetching() {
        let preparer = LeftAlignPreparer::new(NeverFetched);
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
