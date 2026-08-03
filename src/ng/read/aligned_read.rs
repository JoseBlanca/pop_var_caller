//! One read from an alignment file, in its two states, and the conversion
//! between them.
//!
//! - [`RawAlignedRead`] is the read as it comes off the file, undecoded: a flag
//!   and a mapping quality readable without unpacking anything else.
//!   [`NoodlesRawAlignedRead`] is the one implementation, a view of a reused
//!   noodles `RecordBuf`.
//! - [`AlignedRead`] is the decoded one: the fields the caller works with, plus
//!   the read group it came from.
//!
//! They live together because they are one thing in two states, and
//! [`decode_record`] — the conversion — was already here. Deciding *which* reads
//! get converted belongs to `read/filtering.rs` and stays there: this module
//! holds no keep-or-drop rule.
//!
//! Design: `doc/devel/ng/spec/read_filtering_stages.md` §6,
//! `doc/devel/ng/arch/read_filtering_stages.md` §1, §3.1.
//!
//! **ng-owned, modelled on production's `MappedRead`**
//! ([`crate::bam::alignment_input::MappedRead`]). ng carried production's type
//! until read groups arrived; the identifier had to ride on the read, and adding
//! it to production's struct would have made production carry an ng concept for
//! ng's benefit alone. The difference is exactly one field: production's
//! `source_file_index` (which file, positionally, within one sample) is replaced
//! by `read_group` (which `@RG`, run-wide) — a strictly better answer to the
//! same question, since a read group knows its own file.
//!
//! The decode is a copy of production's assembly, **not** of its logic:
//! `compute_adaptor_boundary` and `cigar_to_ops` stay production's and are
//! called from here. Field-for-field parity against production's decode of the
//! same record is a test in this module.
//!
//! Design: `doc/devel/ng/spec/read_groups.md` §8,
//! `doc/devel/ng/arch/read_groups.md` §1.4.

use std::io;

use noodles_sam::alignment::RecordBuf;

use crate::bam::alignment_input::{cigar_to_ops, compute_adaptor_boundary};
use crate::ng::types::{MapQual, ReadGroupId};
use crate::pileup::walker::CigarOp;

/// One alignment record as it comes off the file, undecoded — the seam that lets
/// the flag/MAPQ cascade (#1–#6) run *before* the record is decoded.
/// [`Self::flag`] and [`Self::mapq`] are cheap field reads on the still-packed
/// record; [`Self::decode`] is the expensive phase (base/quality decode +
/// adaptor-boundary annotation → [`AlignedRead`]), run only on the reads that
/// clear the pre-decode gate. It borrows `&self`, not `self`, so the underlying
/// buffer stays reusable for the next read.
///
/// **An unmapped read is one of these.** Filter #5 rejects it, so it exists here
/// before it is rejected, and SAM calls every line an alignment record — unmapped
/// ones included. [`AlignedRead`] has no such case: [`decode_record`] refuses a
/// record with no reference sequence id or no alignment start.
pub trait RawAlignedRead {
    /// The SAM flag bitfield — the same `u16` [`AlignedRead::flag`] carries; read
    /// by filters #1, #3–#6.
    fn flag(&self) -> u16;
    /// The SAM mapping quality; filter #2. "Unavailable" (SAM `0xFF`) resolves to
    /// [`MapQual`]`(0)`, matching production, so any non-zero minimum drops it.
    fn mapq(&self) -> MapQual;
    /// Decode the record into an owned [`AlignedRead`] (filters #7–#9 read the
    /// result). Fallible because the decode is: a record reaching here has passed
    /// the pre-decode cascade (so it is mapped), but a corrupt record — the unmapped
    /// flag clear yet no position — surfaces here as an `Err` that, like the #8
    /// fetch and the read off the file, is **fatal to the run** rather than a
    /// per-read drop or a panic (spec §7). (The spec's illustrative signature
    /// showed this infallible; it is fallible in practice because the reused
    /// production decoder is.)
    fn decode(&self) -> io::Result<AlignedRead>;
    /// Which read group this record belongs to, once its source has resolved it.
    ///
    /// Read by the tally, not by any filter — a drop has to be charged to the
    /// read group it came from, and the pre-decode filters run before any
    /// `AlignedRead` exists to carry it. `None` means the source has not stamped
    /// the buffer, which [`Self::decode`] refuses; the tally charges it to no
    /// read group rather than to an arbitrary one.
    fn read_group(&self) -> Option<ReadGroupId>;
}

/// One decoded, filtered read.
///
/// Owned and free of lifetimes: a read outlives the reused buffer it was
/// decoded from, and the merge holds several at once.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignedRead {
    pub qname: Vec<u8>,
    pub flag: u16,
    /// Index into the canonical contig table — equal to the file's `ref_id`,
    /// which the open gate proved.
    pub ref_id: usize,
    /// 1-based leftmost mapped position.
    pub pos: u64,
    pub mapq: u8,
    pub cigar: Vec<CigarOp>,
    /// Uppercase ACGTN.
    pub seq: Vec<u8>,
    /// Phred 0-93, raw base quality before BAQ capping.
    pub qual: Vec<u8>,
    pub mate_ref_id: Option<usize>,
    pub mate_pos: Option<u64>,
    /// Reference position of the first base on this read that lies inside the
    /// mate-pair adaptor (1-based), or `None` where the boundary cannot be
    /// computed. Semantics and the computation are production's.
    pub adaptor_boundary: Option<u32>,
    /// Which read group this read came from — the library it was prepared in,
    /// and through the run's table, the file it was read from.
    pub read_group: ReadGroupId,
}

/// Decode one record into an [`AlignedRead`], attributing it to `read_group`.
///
/// Fallible for the reason production's is: a record that cleared the
/// pre-decode gate can still be corrupt — unmapped flag clear yet no position,
/// or no reference sequence id — and that is fatal to the run rather than a
/// per-read drop.
pub(crate) fn decode_record(
    record: &RecordBuf,
    read_group: ReadGroupId,
) -> io::Result<AlignedRead> {
    let qname = record
        .name()
        .map(|name| AsRef::<[u8]>::as_ref(name).to_vec())
        .unwrap_or_default();
    let flag = record.flags().bits();

    // Both failures name the read. Nothing above this point can: the error
    // travels up as `ReadFilterError::Decode` and then `AlignmentFileError::
    // Filter`, neither of which knows which record was being decoded, so a
    // message built here without the name leaves a run-killing error with
    // nothing to look up. The name is already in hand, one line above.
    let named = || String::from_utf8_lossy(&qname).into_owned();
    let ref_id = record.reference_sequence_id().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("read '{}' has no reference sequence id", named()),
        )
    })?;
    let pos = record
        .alignment_start()
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "read '{}' on reference sequence {ref_id} has no alignment start",
                    named()
                ),
            )
        })?
        .get() as u64;
    let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
    let cigar = cigar_to_ops(record.cigar());
    let seq: Vec<u8> = record
        .sequence()
        .as_ref()
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect();
    let qual = record.quality_scores().as_ref().to_vec();
    let mate_ref_id = record.mate_reference_sequence_id();
    let mate_pos = record.mate_alignment_start().map(|pos| pos.get() as u64);
    let adaptor_boundary = compute_adaptor_boundary(
        flag,
        ref_id,
        pos,
        seq.len() as u32,
        mate_ref_id,
        mate_pos,
        record.template_length(),
    );

    Ok(AlignedRead {
        qname,
        flag,
        ref_id,
        pos,
        mapq,
        cigar,
        seq,
        qual,
        mate_ref_id,
        mate_pos,
        adaptor_boundary,
        read_group,
    })
}

/// An ng-owned adapter viewing a noodles [`RecordBuf`] as a [`RawAlignedRead`].
/// Whatever reads the file refills the inner `record` in place;
/// [`RawAlignedRead::flag`]/[`RawAlignedRead::mapq`] are cheap field reads on the
/// undecoded record, and [`RawAlignedRead::decode`] runs ng's own
/// [`decode_record`] (which still calls production's `compute_adaptor_boundary`
/// and `cigar_to_ops`). `read_group` is the read group the decoded read is
/// attributed to; whatever fills the buffer sets it before each read.
///
/// ng-owned by design (spec §7, decision a): the dependency points ng → existing
/// code, so production never learns about ng.
///
/// Both fields are `pub(crate)` internal state, not a caller-set surface:
/// `read_group` is cleared and (re)stamped on every read by whatever fills the
/// buffer, and `record` is refilled in place — a caller-written value would just
/// be overwritten. No `Clone`: the whole point is to reuse one buffer, and
/// cloning would deep-copy the `RecordBuf` it exists to avoid copying.
#[derive(Debug)]
pub struct NoodlesRawAlignedRead {
    /// The reused undecoded record buffer.
    pub(crate) record: RecordBuf,
    /// The read group this record belongs to, stamped onto the decoded
    /// [`AlignedRead`].
    ///
    /// **Cleared when the buffer is handed out and set again after each record
    /// is filled** — in that order, and it matters. The buffer is reused, so a
    /// source that failed to stamp would otherwise leave the *previous*
    /// record's group in place and attribute this read to it, silently. Clearing
    /// first turns that into a refusal at [`RawAlignedRead::decode`], which is
    /// what the `Option` is for; without it the guard only ever caught a missed
    /// stamp on the very first record.
    pub(crate) read_group: Option<ReadGroupId>,
}

/// Hand-written rather than derived, because a fresh buffer has **no** read
/// group and must not pretend to one.
///
/// `ReadGroupId(0)` would be the obvious filler and is the trap: zero is not a
/// sentinel, it is the first identifier the run's table mints, so a source that
/// failed to stamp would attribute every read to a real library — silently, and
/// on exactly the grouping the error model is fitted per. `None` makes that a
/// loud failure at the one place that reads the field.
impl Default for NoodlesRawAlignedRead {
    fn default() -> Self {
        Self {
            record: RecordBuf::default(),
            read_group: None,
        }
    }
}

impl RawAlignedRead for NoodlesRawAlignedRead {
    fn flag(&self) -> u16 {
        self.record.flags().bits()
    }

    fn mapq(&self) -> MapQual {
        // `map(u8::from).unwrap_or(0)`: SAM `0xFF` (unavailable) decodes to `None`
        // → treated as MAPQ 0, matching production `classify_pre_decode`.
        MapQual(self.record.mapping_quality().map(u8::from).unwrap_or(0))
    }

    fn read_group(&self) -> Option<ReadGroupId> {
        self.read_group
    }

    fn decode(&self) -> io::Result<AlignedRead> {
        // A buffer reaching decode without a read group means whatever filled it
        // skipped its stamp. Refusing is the point: the alternative is a read
        // silently attributed to whichever library the filler happened to name.
        let read_group = self.read_group.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "a record reached decode with no read group: its source did not stamp one",
            )
        })?;
        decode_record(&self.record, read_group)
    }
}

impl AlignedRead {
    /// The same read as production's [`MappedRead`](crate::bam::alignment_input::MappedRead).
    ///
    /// **The read group does not survive this conversion**, because production's
    /// type has nowhere to put it. It exists for one boundary: the generic
    /// path's read preparation hands the read to production's
    /// `prepare_passthrough`, which returns a production `PreparedRead`. So on
    /// that path the identifier stops here — the deferred item in spec §12, and
    /// the reason an ng-owned prepared read is the next thing this design wants.
    /// The STR path never calls this: it works on `AlignedRead` directly.
    ///
    /// `source_file_index` is `0` and means nothing. Production sets it from a
    /// sample's file list, ng no longer has one on the read, and the only
    /// consumer of this conversion ignores it.
    pub(crate) fn into_mapped_read(self) -> crate::bam::alignment_input::MappedRead {
        // Destructured, not field-by-field off `self`: a field added to
        // `AlignedRead` later then fails to compile here instead of being
        // silently dropped on the way to production's type.
        let Self {
            qname,
            flag,
            ref_id,
            pos,
            mapq,
            cigar,
            seq,
            qual,
            mate_ref_id,
            mate_pos,
            adaptor_boundary,
            read_group: _,
        } = self;

        crate::bam::alignment_input::MappedRead {
            qname,
            flag,
            ref_id,
            pos,
            mapq,
            cigar,
            seq,
            qual,
            mate_ref_id,
            mate_pos,
            adaptor_boundary,
            source_file_index: 0,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::{
        FLAG_DUPLICATE, FLAG_MATE_REVERSE_STRAND, FLAG_PAIRED, FLAG_REVERSE_STRAND,
        record_buf_to_mapped_read,
    };
    use crate::ng::read::input::test_fixtures::{read_named, read_named_with_length};

    use noodles_core::Position as RecordPosition;
    use noodles_sam::alignment::record::cigar::Op;
    use noodles_sam::alignment::record::cigar::op::Kind;
    use noodles_sam::alignment::record::{Flags, MappingQuality};
    use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

    /// A read whose mate makes `compute_adaptor_boundary` return `Some`.
    ///
    /// Every precondition has to hold at once — paired, both mates mapped, a
    /// non-zero template length, mates on opposite strands and the same
    /// reference, and a molecule *shorter* than the read so part of it is
    /// adaptor. Miss one and the function returns `None`, which is what made
    /// the first version of the parity test compare `None` with `None` and
    /// prove nothing.
    fn readthrough_pair() -> RecordBuf {
        RecordBuf::builder()
            .set_name(b"pair-1")
            .set_flags(Flags::from(FLAG_PAIRED | FLAG_MATE_REVERSE_STRAND))
            .set_reference_sequence_id(0)
            .set_alignment_start(RecordPosition::try_from(20).unwrap())
            .set_mate_reference_sequence_id(0)
            .set_mate_alignment_start(RecordPosition::try_from(24).unwrap())
            .set_template_length(6)
            .set_mapping_quality(MappingQuality::new(60).unwrap())
            .set_cigar([Op::new(Kind::Match, 10)].into_iter().collect())
            .set_sequence(Sequence::from(b"ACGTACGTAC".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 10]))
            .build()
    }

    /// The reverse-strand half of the same pair — the *other* branch of
    /// `compute_adaptor_boundary`, which computes the boundary from the mate's
    /// start rather than from this read's.
    fn readthrough_pair_reverse() -> RecordBuf {
        let mut record = readthrough_pair();
        *record.flags_mut() = Flags::from(FLAG_PAIRED | FLAG_REVERSE_STRAND);
        *record.template_length_mut() = -6;
        record
    }

    /// Lower-case bases, a clipped and gapped CIGAR, no mapping quality and no
    /// name — the arms the uniform fixtures never reach: `to_ascii_uppercase`,
    /// every `CigarOp` variant but `Match`, `unwrap_or(0)` and
    /// `unwrap_or_default()`.
    fn awkward_record() -> RecordBuf {
        RecordBuf::builder()
            .set_flags(Flags::empty())
            .set_reference_sequence_id(1)
            .set_alignment_start(RecordPosition::try_from(31).unwrap())
            .set_cigar(
                [
                    Op::new(Kind::SoftClip, 2),
                    Op::new(Kind::Match, 3),
                    Op::new(Kind::Insertion, 1),
                    Op::new(Kind::Deletion, 2),
                    Op::new(Kind::Match, 4),
                ]
                .into_iter()
                .collect(),
            )
            .set_sequence(Sequence::from(b"acgtnACGTn".to_vec()))
            .set_quality_scores(QualityScores::from(vec![2u8; 10]))
            .build()
    }

    /// **The parity oracle for the copy.** `decode_record` is production's
    /// assembly rewritten; a field silently dropped or mis-mapped would produce
    /// wrong reads rather than a crash, and no other test would notice. Every
    /// field is compared against production's decode of the same record, one by
    /// one — `assert_eq!` on the whole struct is impossible (the types differ)
    /// and would hide which field moved even if it were.
    ///
    /// **The fixtures have to make each field differ from its default**, or the
    /// comparison proves nothing: an earlier version drove three unpaired
    /// single-op reads, so `mate_ref_id`, `mate_pos` and `adaptor_boundary`
    /// were `None` on both sides every time and deleting the
    /// `compute_adaptor_boundary` call outright would have kept it green.
    #[test]
    fn every_field_matches_productions_decode_of_the_same_record() {
        let records = [
            ("plain", read_named("r1", 0, 1)),
            ("second contig", read_named("r2", 1, 40)),
            ("long", read_named_with_length("r3", 0, 7, 40)),
            ("readthrough, forward", readthrough_pair()),
            ("readthrough, reverse", readthrough_pair_reverse()),
            (
                "clipped, gapped, lower-case, unmapped-quality",
                awkward_record(),
            ),
        ];

        for (case, record) in records {
            let ours = decode_record(&record, ReadGroupId(3))
                .unwrap_or_else(|error| panic!("{case} decodes: {error}"));
            let production = record_buf_to_mapped_read(&record, 0)
                .unwrap_or_else(|error| panic!("{case} decodes in production: {error}"));

            assert_eq!(ours.qname, production.qname, "{case}: qname");
            assert_eq!(ours.flag, production.flag, "{case}: flag");
            assert_eq!(ours.ref_id, production.ref_id, "{case}: ref_id");
            assert_eq!(ours.pos, production.pos, "{case}: pos");
            assert_eq!(ours.mapq, production.mapq, "{case}: mapq");
            assert_eq!(ours.cigar, production.cigar, "{case}: cigar");
            assert_eq!(ours.seq, production.seq, "{case}: seq");
            assert_eq!(ours.qual, production.qual, "{case}: qual");
            assert_eq!(
                ours.mate_ref_id, production.mate_ref_id,
                "{case}: mate_ref_id"
            );
            assert_eq!(ours.mate_pos, production.mate_pos, "{case}: mate_pos");
            assert_eq!(
                ours.adaptor_boundary, production.adaptor_boundary,
                "{case}: adaptor_boundary"
            );
            assert_eq!(
                ours.read_group,
                ReadGroupId(3),
                "{case}: the one field production has no counterpart for"
            );
        }
    }

    /// The fixtures above are only worth anything if they reach the branches
    /// they were built for. This asserts they do — without it, a later edit
    /// could quietly return them to the all-`None` shape that made the oracle
    /// vacuous, and every parity assertion would still pass.
    #[test]
    fn the_parity_fixtures_reach_the_branches_they_exist_for() {
        let forward = decode_record(&readthrough_pair(), ReadGroupId(0)).expect("decodes");
        assert_eq!(
            forward.adaptor_boundary,
            Some(26),
            "a forward readthrough puts the boundary at start + |tlen|"
        );
        assert_eq!(forward.mate_ref_id, Some(0));
        assert_eq!(forward.mate_pos, Some(24));

        let reverse = decode_record(&readthrough_pair_reverse(), ReadGroupId(0)).expect("decodes");
        assert_eq!(
            reverse.adaptor_boundary,
            Some(23),
            "a reverse readthrough puts it just before the mate's start"
        );

        let awkward = decode_record(&awkward_record(), ReadGroupId(0)).expect("decodes");
        assert_eq!(awkward.seq, b"ACGTNACGTN".to_vec(), "bases are upper-cased");
        assert_eq!(awkward.mapq, 0, "an absent mapping quality reads as 0");
        assert!(awkward.qname.is_empty(), "an absent name reads as empty");
        assert!(
            awkward.cigar.len() > 1,
            "the CIGAR reaches ops other than Match"
        );
    }

    /// A corrupt record is fatal, not a silent default — the same contract
    /// production's decode has — and the message names the read, because
    /// nothing above this point can.
    #[test]
    fn a_record_with_no_position_fails_naming_the_read() {
        let mut record = read_named("r1", 0, 1);
        *record.alignment_start_mut() = None;

        let error = decode_record(&record, ReadGroupId(0)).expect_err("a corrupt record is fatal");

        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        let message = error.to_string();
        assert!(
            message.contains("r1"),
            "the message names the read: {message}"
        );
        assert!(
            message.contains("alignment start"),
            "and says what was missing: {message}"
        );
    }

    /// A record with the flag and mapping quality a reader would have filled the
    /// buffer with, and a body that decodes cleanly. Not `read_named`: the tests
    /// below are about the adapter reading *those two fields specifically*, so
    /// both have to differ from that fixture's defaults.
    ///
    /// **Everything else it sets is asserted somewhere**, which is the point: a
    /// fixture field nothing reads looks like it is being proved and is not.
    fn record_with_mapq_and_flags(mapq: u8, flags: Flags) -> RecordBuf {
        RecordBuf::builder()
            .set_name(b"r1")
            .set_reference_sequence_id(0)
            .set_flags(flags)
            .set_mapping_quality(MappingQuality::new(mapq).expect("mapq in range"))
            .set_alignment_start(RecordPosition::try_from(5).unwrap())
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
            .set_sequence(Sequence::from(vec![b'A'; 4]))
            .set_quality_scores(QualityScores::from(vec![30u8; 4]))
            .build()
    }

    #[test]
    fn noodles_raw_aligned_read_reads_flag_mapq_and_decodes() {
        let record = record_with_mapq_and_flags(42, Flags::from(FLAG_DUPLICATE));
        let raw = NoodlesRawAlignedRead {
            record,
            read_group: Some(ReadGroupId(7)),
        };
        // Cheap pre-decode reads.
        assert_eq!(raw.flag(), FLAG_DUPLICATE);
        assert_eq!(raw.mapq(), MapQual(42));
        // Decode produces the AlignedRead the cascade's phase two consumes.
        let read = raw.decode().unwrap();
        assert_eq!(read.ref_id, 0);
        assert_eq!(read.pos, 5);
        assert_eq!(read.mapq, 42);
        assert_eq!(read.seq, b"AAAA");
        assert_eq!(read.read_group, ReadGroupId(7));
        // The two fields the fixture set that nothing used to read. Without
        // these the helper's name and quality scores are inert, so the adapter's
        // pass-through of them is asserted by nothing.
        assert_eq!(read.qname, b"r1");
        assert_eq!(read.qual, vec![30u8; 4]);
    }

    #[test]
    fn noodles_raw_aligned_read_maps_unavailable_mapq_to_zero() {
        // A record with no mapping quality (SAM 0xFF) → MapQual(0).
        let record = RecordBuf::builder()
            .set_reference_sequence_id(0)
            .set_flags(Flags::from(FLAG_PAIRED))
            .set_alignment_start(RecordPosition::try_from(1usize).unwrap())
            .set_cigar([Op::new(Kind::Match, 4)].into_iter().collect())
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30u8; 4]))
            .build();
        let raw = NoodlesRawAlignedRead {
            record,
            read_group: Some(ReadGroupId(0)),
        };
        assert_eq!(raw.mapq(), MapQual(0));
    }

    /// A corrupt record is fatal through the adapter too, and the error that
    /// surfaces is the **decoder's** — the buffer is stamped, so the read-group
    /// guard cannot fire first.
    ///
    /// **This test used to drive `NoodlesRawAlignedRead::default()`**, whose
    /// `read_group` is `None`, so it returned at the read-group guard and never
    /// reached `decode_record` at all — while its name, its comment and the
    /// impl report all said it exercised the missing-position path. It asserted
    /// only `err.kind()`, which both paths produce, so it passed either way.
    /// Stamping the buffer is what makes the name true; asserting the message is
    /// what keeps it true.
    #[test]
    fn noodles_raw_aligned_read_decode_errors_on_a_stamped_record_with_no_position() {
        let mut record = record_with_mapq_and_flags(42, Flags::from(FLAG_DUPLICATE));
        *record.alignment_start_mut() = None;
        let raw = NoodlesRawAlignedRead {
            record,
            read_group: Some(ReadGroupId(7)),
        };

        let error = raw.decode().expect_err("a corrupt record is fatal");

        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(
            error.to_string().contains("alignment start"),
            "the decoder's own error surfaces, not the read-group guard's: {error}"
        );
    }

    /// **The refusal the `Option<ReadGroupId>` exists for.** A buffer that
    /// reaches `decode` without a stamp is refused, because `ReadGroupId(0)` is
    /// the first id the run's table mints — a real library — so a reader that
    /// forgot to stamp would otherwise attribute every one of its reads to it,
    /// silently, on exactly the grouping the error model is fitted per.
    ///
    /// **The fixture has to decode cleanly or this proves nothing** — which is
    /// how the guard went unpinned: the test that looked like it covered this
    /// used a buffer that was invalid in several respects at once, so it could
    /// not tell the read-group refusal from the decoder's own complaint.
    /// Replacing the whole guard with `unwrap_or(ReadGroupId(0))` left all 2,837
    /// tests green.
    #[test]
    fn noodles_raw_aligned_read_decode_refuses_a_buffer_with_no_read_group() {
        let record = record_with_mapq_and_flags(42, Flags::from(FLAG_DUPLICATE));
        assert!(
            decode_record(&record, ReadGroupId(0)).is_ok(),
            "the fixture must be decodable, or the missing stamp is not what fails"
        );
        let raw = NoodlesRawAlignedRead {
            record,
            read_group: None,
        };

        let error = raw.decode().expect_err("an unstamped buffer is refused");

        assert_eq!(error.kind(), io::ErrorKind::InvalidData);
        assert!(
            error.to_string().contains("read group"),
            "and the message names the missed stamp: {error}"
        );
    }

    /// A fresh buffer reports **no** read group — the whole argument of the
    /// hand-written [`Default`], asserted rather than left to its doc comment.
    ///
    /// Read through the trait accessor, not the field: the accessor is what the
    /// tally calls, and it was equally unpinned. Both mutations this kills —
    /// `Default` filling `Some(ReadGroupId(0))`, and the accessor rewritten so
    /// it can never answer `None` — passed the whole 2,837-test suite before it
    /// existed. The neighbouring tests assert on the field, which is why.
    #[test]
    fn noodles_raw_aligned_read_default_reports_no_read_group() {
        assert_eq!(
            RawAlignedRead::read_group(&NoodlesRawAlignedRead::default()),
            None,
            "ReadGroupId(0) is a real id, so a fresh buffer must not pretend to one"
        );
    }

    /// The conversion to production's type is exhaustive by construction (it
    /// destructures), but nothing checks that it maps each field to its
    /// counterpart rather than to a neighbour of the same type.
    #[test]
    fn the_conversion_to_productions_read_moves_every_field_to_its_counterpart() {
        let ours = decode_record(&readthrough_pair(), ReadGroupId(7)).expect("decodes");
        let expected = ours.clone();
        let theirs = ours.into_mapped_read();

        assert_eq!(theirs.qname, expected.qname);
        assert_eq!(theirs.flag, expected.flag);
        assert_eq!(theirs.ref_id, expected.ref_id);
        assert_eq!(theirs.pos, expected.pos);
        assert_eq!(theirs.mapq, expected.mapq);
        assert_eq!(theirs.cigar, expected.cigar);
        assert_eq!(theirs.seq, expected.seq);
        assert_eq!(theirs.qual, expected.qual);
        assert_eq!(theirs.mate_ref_id, expected.mate_ref_id);
        assert_eq!(theirs.mate_pos, expected.mate_pos);
        assert_eq!(theirs.adaptor_boundary, expected.adaptor_boundary);
        assert_eq!(
            theirs.source_file_index, 0,
            "the read group has nowhere to go, and this field means nothing here"
        );
    }
}
