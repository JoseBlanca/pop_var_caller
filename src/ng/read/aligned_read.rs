//! [`AlignedRead`] — ng's decoded read.
//!
//! One record from an alignment file, decoded into the fields the caller works
//! with, plus the read group it came from.
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

use noodles_sam as sam;

use crate::bam::alignment_input::{cigar_to_ops, compute_adaptor_boundary};
use crate::ng::types::ReadGroupId;
use crate::pileup::walker::CigarOp;

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
    record: &sam::alignment::RecordBuf,
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
    use crate::bam::alignment_input::record_buf_to_mapped_read;
    use crate::ng::read::input::test_fixtures::{read_named, read_named_with_length};

    use noodles_core::Position as RecordPosition;
    use noodles_sam::alignment::RecordBuf;
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
        const FLAG_PAIRED: u16 = 0x1;
        const FLAG_MATE_REVERSE: u16 = 0x20;

        RecordBuf::builder()
            .set_name(b"pair-1")
            .set_flags(Flags::from(FLAG_PAIRED | FLAG_MATE_REVERSE))
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
        const FLAG_PAIRED: u16 = 0x1;
        const FLAG_REVERSE: u16 = 0x10;

        let mut record = readthrough_pair();
        *record.flags_mut() = Flags::from(FLAG_PAIRED | FLAG_REVERSE);
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
