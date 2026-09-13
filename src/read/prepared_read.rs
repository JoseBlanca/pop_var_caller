//! [`PreparedRead`] — ng's prepared read, the pileup walk's input.
//!
//! **ng-owned, copied from production's `pileup/walker/mod.rs`** and extended
//! with one field: `read_group`. The difference is the whole reason the copy
//! exists — a read's library membership is a property of the read, exactly like
//! its MAPQ or its strand, and production's type predates read groups and has
//! nowhere to put it. Adding the field there would mean editing a frozen module
//! for ng's benefit alone; copying reaches the same place and edits production
//! zero times (`doc/devel/ng/spec/locus_generation_pileup.md` §6).
//!
//! Copying the *read* is also what settles copying the *walker*: every module
//! under `pileup/walker/` names `PreparedRead` in its signatures, so an ng-owned
//! read type reaches all of them whatever else is decided (spec §3).
//!
//! The fields, their invariants and [`PreparedRead::length`] are production's,
//! transcribed rather than re-derived. **So is the construction path, which ng now
//! carries**: [`prepare_passthrough`] is production's `--no-baq` arm copied here on
//! 2026-09-12, deriving the same four values from the same flags and CIGAR and returning
//! ng's read with its group attached. Until then ng minted production's read and converted
//! it; that conversion survived as a test-only bridge for the parity oracles until promotion
//! step C26 deleted it.
//! Either way ng re-derives none of the per-field wiring, and the read-preparation parity
//! fixture is what says so.
//!
//! Design: `doc/devel/ng/spec/locus_generation_pileup.md` §6,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*.

use std::sync::Arc;

use crate::bam::alignment_input::{
    FLAG_FIRST_OF_PAIR, FLAG_PAIRED, FLAG_REVERSE_STRAND, MappedRead,
};
use crate::types::ReadGroupId;
// The decoder's, which is where a CIGAR run is built and read back. It sat in
// `pileup::walker` until 2026-09-12, and the old path was the misleading one: this is
// crate-wide read vocabulary, not the walk's.
use crate::bam::alignment_input::CigarOp;

// ---------------------------------------------------------------------
// Minting one
// ---------------------------------------------------------------------

/// Turn a decoded read into the walk's input: derive what the walk needs, carry the rest
/// across, and attach the read group.
///
/// **"Passthrough" is about the base qualities.** Production's function of this name is its
/// `--no-baq` arm — the one that hands the read's own qualities to the walk instead of
/// replacing them with BAQ-adjusted ones. ng defers BAQ entirely
/// (`doc/devel/ng/spec/read_preparation.md` §1), so this is the only arm there is, and the
/// name is kept because the parity fixture compares against the arm it names.
///
/// Four things are derived and the rest is moved unchanged: where the alignment ends, the
/// read's name as a shared string, its mate role, and its mapping quality on a log scale.
///
/// # ng's copy of production's `baq_engine::prepare_passthrough`
///
/// Copied on 2026-09-12 with its four helpers, so that ng mints its own read rather than
/// minting production's and converting. The arithmetic is production's: the same reference
/// span, the same flag tests, the same `ln` of the same Phred. What it does differently is
/// take `read_group` and return ng's [`PreparedRead`] — which is the one field production's
/// type has nowhere to put, and the reason ng owns this type at all.
pub(crate) fn prepare_passthrough(
    read: MappedRead,
    chrom_id: u32,
    read_group: ReadGroupId,
) -> PreparedRead {
    let bq_baq = read.qual.clone();
    let alignment_start = read.pos as u32;
    let alignment_end = alignment_end_of(alignment_start, &read.cigar);
    let qname = qname_to_arc(&read.qname);
    let mate_role = mate_role_of(read.flag);
    let is_reverse_strand = read.flag & FLAG_REVERSE_STRAND != 0;
    let mq_log_err = phred_to_ln_perr(read.mapq);
    let mapq = read.mapq;
    PreparedRead {
        chrom_id,
        alignment_start,
        alignment_end,
        cigar: read.cigar,
        seq: read.seq,
        bq_baq,
        mq_log_err,
        mapq,
        is_reverse_strand,
        qname,
        mate_role,
        adaptor_boundary: read.adaptor_boundary,
        read_group,
    }
}

/// The last reference base the alignment covers, 1-based and inclusive.
///
/// Sums the CIGAR steps that consume reference. **Every variant is named rather than caught
/// by a wildcard**, so a step added to [`CigarOp`] has to be classified deliberately instead
/// of silently contributing nothing and shortening every read that carries it.
fn alignment_end_of(start_1: u32, cigar: &[CigarOp]) -> u32 {
    let mut ref_span: u32 = 0;
    for op in cigar {
        match *op {
            CigarOp::Match(l)
            | CigarOp::SeqMatch(l)
            | CigarOp::SeqMismatch(l)
            | CigarOp::Deletion(l)
            | CigarOp::Skip(l) => ref_span = ref_span.saturating_add(l),
            CigarOp::Insertion(_)
            | CigarOp::SoftClip(_)
            | CigarOp::HardClip(_)
            | CigarOp::Padding(_) => {}
        }
    }
    start_1.saturating_add(ref_span).saturating_sub(1)
}

/// The read's name, shared rather than copied — the walk keeps one in its pending-mate table
/// beside the read itself.
///
/// SAM says a name is printable ASCII, so the ordinary path is one allocation. A name that is
/// not valid UTF-8 goes through `from_utf8_lossy`, whose owned string is stolen rather than
/// copied again.
fn qname_to_arc(qname: &[u8]) -> Arc<str> {
    match std::str::from_utf8(qname) {
        Ok(s) => Arc::<str>::from(s),
        Err(_) => Arc::<str>::from(String::from_utf8_lossy(qname)),
    }
}

/// Which half of a pair this read is, from its SAM flags.
///
/// Unpaired if `0x1` is clear. Otherwise first-of-pair if `0x40` is set and second-of-pair if
/// it is not — which treats `0x80` as redundant, so a malformed record with both bits set or
/// both clear is read as second-of-pair rather than refused here. Rejecting those is the
/// input stage's job.
fn mate_role_of(flag: u16) -> MateRole {
    if flag & FLAG_PAIRED == 0 {
        MateRole::Solo
    } else if flag & FLAG_FIRST_OF_PAIR != 0 {
        MateRole::FirstOfPair
    } else {
        MateRole::SecondOfPair
    }
}

/// A Phred mapping quality as the natural log of the probability the placement is wrong:
/// `ln(10^(-Q/10))`.
///
/// `Q = 0` means *no information about the placement*, so the probability is 1 and its log is
/// 0 — which the general formula gives anyway, and which is written out because the caller
/// reads it as a special case.
fn phred_to_ln_perr(q: u8) -> f64 {
    if q == 0 {
        return 0.0;
    }
    -(q as f64) * std::f64::consts::LN_10 / 10.0
}

// ---------------------------------------------------------------------
// MateRole
// ---------------------------------------------------------------------

/// SAM-flag-derived mate role for a [`PreparedRead`]. Collapses the
/// SAM `0x1` (paired) and `0x40` (first segment) bits into a single
/// field so combinations the walker never wants to see — notably
/// "solo + first-of-pair" — cannot be constructed.
///
/// "First" / "Second" here is the SAM-flag sense (which segment of
/// the template the read came from), not the walker arrival-order
/// sense — under coordinate-sorted input the second-of-pair can be
/// the first to reach the walker.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MateRole {
    /// SAM flag `0x1` unset. Treated as solo by the walker and never
    /// registered in `pending_mates`.
    Solo,
    /// SAM flag `0x1` set, flag `0x40` set — first segment in the
    /// template.
    FirstOfPair,
    /// SAM flag `0x1` set, flag `0x40` unset — last segment in the
    /// template.
    SecondOfPair,
}

impl MateRole {
    /// True for both `FirstOfPair` and `SecondOfPair` — i.e. SAM
    /// flag `0x1` is set.
    ///
    /// An exhaustive `match` where production writes `!matches!`: **this enum is ng's
    /// to extend and production's is frozen**, so a fourth variant here must be
    /// classified deliberately rather than defaulted into "paired". Behaviour on the
    /// three current variants is production's exactly, which
    /// `the_mate_role_predicates_are_productions` pins.
    pub fn is_paired(self) -> bool {
        match self {
            MateRole::Solo => false,
            MateRole::FirstOfPair | MateRole::SecondOfPair => true,
        }
    }

    /// True only for `FirstOfPair`. Used by mate-overlap tie-breaks
    /// on equal-BQ positions; on `Solo` and `SecondOfPair` it is
    /// false.
    ///
    /// Exhaustive for the reason [`is_paired`](Self::is_paired) gives — and this is the
    /// predicate with teeth: it is the deterministic tie-break on equal-BQ mate-overlap
    /// positions, so a variant defaulted to `false` here is a silently wrong tie-break.
    pub fn is_first_of_pair(self) -> bool {
        match self {
            MateRole::FirstOfPair => true,
            MateRole::Solo | MateRole::SecondOfPair => false,
        }
    }
}

// ---------------------------------------------------------------------
// PreparedRead
// ---------------------------------------------------------------------

/// A read after the upstream filter and BAQ stages: every field the
/// walker needs is owned and pre-decoded. The walker treats `bq_baq`
/// as opaque per-base BQ — it does not re-apply BAQ — so for tests
/// or until the BAQ stage lands, raw BAM BQ flows through unchanged.
///
/// All bytes are uppercase ASCII over `{A,C,G,T,N}` (`seq`) or Phred
/// 0–93 (`bq_baq`). `qname` is shared as `Arc<str>` so cheap clones
/// can sit in the `pending_mates` map alongside the `ActiveRead`
/// without the bytes being duplicated.
///
/// # Construction: every field, every time
///
/// Deliberately **not** `#[non_exhaustive]`, and deliberately **without a `Default`
/// impl** — production's own stated reason, which ng inherits along with the type.
/// `PreparedRead` is the input contract: callers populate it with concrete bytes per the
/// field-level docs, and adding a new field *should* force every caller (test, bench,
/// constructor) to update its literal explicitly. `#[non_exhaustive]` would push callers
/// toward `..Default::default()`, which is exactly the silent-absorb hazard the
/// refactor-safety rule guards against — and a `Default` impl would hand them that
/// escape hatch directly. The 28 `read_group` lines this port added to the copied
/// fixtures are that decision working.
///
/// *(Production states the same reason in a `//` comment, which rustdoc drops. It is a
/// doc comment here because a caller meeting a compile error on a new field should be
/// able to read that the error is the point.)*
#[derive(Debug, Clone)]
pub struct PreparedRead {
    /// Index into the merged `ContigList`.
    pub chrom_id: u32,
    /// 1-based reference position of the first reference base the
    /// read covers (CIGAR `M`/`=`/`X` aligned).
    pub alignment_start: u32,
    /// 1-based reference position of the last reference base the
    /// read covers, inclusive. Cached at decode time so the walker
    /// does not re-walk the CIGAR to compute it.
    pub alignment_end: u32,
    pub cigar: Vec<CigarOp>,
    /// Read bases. Length matches the read-consuming CIGAR ops
    /// (M/=/X/I/S).
    pub seq: Vec<u8>,
    /// BAQ-capped per-base quality, Phred. Same length as `seq`.
    /// Outside test contexts this is `min(BQ, BAQ)`; inside tests
    /// (or before BAQ lands) raw BAM BQ flows through unchanged.
    pub bq_baq: Vec<u8>,
    /// `ln(P_err)` derived from MAPQ once at filter time. Cached
    /// because every event-folding step uses it.
    pub mq_log_err: f64,
    /// Raw BAM mapping quality, preserved alongside `mq_log_err` so
    /// per-allele `mapq_sum` / `mapq_sum_sq` can be accumulated in
    /// the pileup fold. Reads that pass `--min-mapq` carry their
    /// original Phred MAPQ here (BWA-MEM caps at 60 in practice).
    pub mapq: u8,
    pub is_reverse_strand: bool,
    pub qname: Arc<str>,
    /// SAM-flag-derived mate role. Collapses the SAM `0x1` (paired)
    /// and `0x40` (first segment) bits into one field so the
    /// previously-unrepresentable "solo + first-of-pair" combination
    /// is gone at the type level. The walker consults
    /// [`MateRole::is_paired`] to decide whether to register the
    /// qname for second-mate lookup and [`MateRole::is_first_of_pair`]
    /// as the deterministic tie-breaker on equal-BQ mate-overlap
    /// positions.
    pub mate_role: MateRole,
    /// 1-based reference position of the first read base that lies
    /// inside the mate-pair adaptor. `None` when the boundary cannot
    /// be reliably computed (single-end, mate unmapped, geometry
    /// inconsistent, or molecule longer than the read).
    ///
    /// On the **forward** strand any read base at `ref_pos >=
    /// adaptor_boundary` was sequenced *through* the molecule's 3′
    /// end into the far adaptor. On the **reverse** strand any read
    /// base at `ref_pos <= adaptor_boundary` was sequenced through
    /// the 5′ end into the near adaptor. The cursor's Match-emit
    /// sites apply this test direction-aware. See finding `G1` in
    /// `doc/devel/reports/reviews/pileup_gatk_comparison_2026-05-08.md`.
    ///
    /// # Default
    /// `None` *disables the G1 adaptor filter* for this read: no
    /// read base is treated as past the boundary. This is the
    /// safe choice when the upstream cannot trust its mate
    /// geometry — a false-positive filter would drop real
    /// evidence. Callers that need strict filtering must compute
    /// and set a boundary themselves. Mi17 in
    /// `doc/devel/reports/reviews/pileup_2026-05-11.md`.
    pub adaptor_boundary: Option<u32>,
    /// Which read group this read came from — **the field production's type
    /// does not have, and the reason ng owns this one.**
    ///
    /// One `@RG`, i.e. one lane: the finest grain available, carried so the
    /// consumer picks library, experiment or read group rather than inheriting
    /// this step's guess. It rides straight through from
    /// [`AlignedRead::read_group`](super::AlignedRead), untouched by
    /// preparation — nothing here reads it (spec §6).
    pub read_group: ReadGroupId,
}

// ---------------------------------------------------------------------
// ReadLengthError
// ---------------------------------------------------------------------

/// Failure modes for [`PreparedRead::length`]. Carries the raw
/// lengths only — locus context (qname, chrom_id, pos) is added by
/// whichever caller maps this into its own error type (the walker
/// wraps it in `WalkerError::MalformedRead`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReadLengthError {
    /// `seq.len()` and `bq_baq.len()` disagree.
    SeqBqMismatch { seq_len: usize, bq_baq_len: usize },
    /// The CIGAR's read-consuming ops sum to `cigar_consumed` but
    /// the `seq` field has length `seq_len`.
    CigarSeqMismatch { cigar_consumed: u64, seq_len: usize },
}

impl PreparedRead {
    /// The read's length (in `u64`, matching the cigar-sum width),
    /// verified against the three internal representations that
    /// must agree on it:
    ///
    /// 1. `seq.len() == bq_baq.len()`.
    /// 2. The CIGAR's read-consuming ops (M/=/X/I/S) sum to
    ///    `seq.len()`.
    ///
    /// On success the three lengths agree and the common value is
    /// returned. On failure a typed [`ReadLengthError`] reports
    /// which invariant broke.
    ///
    /// M21 in `doc/devel/reports/reviews/pileup_2026-05-09.md`. The cigar cursor
    /// and `decompose` both index `seq[..]` / `bq_baq[..]` using
    /// offsets derived from the CIGAR; a mismatch would otherwise
    /// panic with `slice index out of bounds` and kill the run on
    /// the offending read. The accessor surfaces a typed error
    /// instead so callers can attach locus context and continue.
    pub fn length(&self) -> Result<u64, ReadLengthError> {
        if self.seq.len() != self.bq_baq.len() {
            return Err(ReadLengthError::SeqBqMismatch {
                seq_len: self.seq.len(),
                bq_baq_len: self.bq_baq.len(),
            });
        }
        let cigar_consumed: u64 = self
            .cigar
            .iter()
            .map(|op| match *op {
                CigarOp::Match(n)
                | CigarOp::SeqMatch(n)
                | CigarOp::SeqMismatch(n)
                | CigarOp::Insertion(n)
                | CigarOp::SoftClip(n) => n as u64,
                CigarOp::Deletion(_)
                | CigarOp::Skip(_)
                | CigarOp::HardClip(_)
                | CigarOp::Padding(_) => 0,
            })
            .sum();
        if cigar_consumed != self.seq.len() as u64 {
            return Err(ReadLengthError::CigarSeqMismatch {
                cigar_consumed,
                seq_len: self.seq.len(),
            });
        }
        Ok(cigar_consumed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read::PLACEHOLDER_READ_GROUP;

    /// A read whose three length representations agree: 4 matched bases, 4
    /// qualities, a CIGAR consuming 4.
    fn consistent_read() -> PreparedRead {
        PreparedRead {
            chrom_id: 0,
            alignment_start: 10,
            alignment_end: 13,
            cigar: vec![CigarOp::Match(4)],
            seq: b"ACGT".to_vec(),
            bq_baq: vec![30, 30, 30, 30],
            mq_log_err: -6.0,
            mapq: 60,
            is_reverse_strand: false,
            qname: Arc::from("read1"),
            mate_role: MateRole::Solo,
            adaptor_boundary: None,
            read_group: PLACEHOLDER_READ_GROUP,
        }
    }

    #[test]
    fn a_consistent_read_reports_its_length() {
        assert_eq!(consistent_read().length(), Ok(4));
    }

    /// The degenerate read: no ops, no bases, no qualities. All three
    /// representations agree on zero, so it is a length and not an error — the
    /// `0`/empty boundary the other cases skip past.
    #[test]
    fn a_read_with_no_bases_reports_a_length_of_zero() {
        let mut read = consistent_read();
        read.cigar = Vec::new();
        read.seq = Vec::new();
        read.bq_baq = Vec::new();
        assert_eq!(read.length(), Ok(0));
    }

    /// Only the read-consuming ops count: a deletion and a skip add reference
    /// span, not read bases, and a hard clip and padding add neither. A copy
    /// that folded any of them in would report a length no `seq` could match.
    #[test]
    fn only_the_read_consuming_ops_count_toward_the_length() {
        let mut read = consistent_read();
        read.cigar = vec![
            CigarOp::HardClip(3),
            CigarOp::SoftClip(1),
            CigarOp::Match(1),
            CigarOp::Insertion(1),
            CigarOp::Deletion(2),
            CigarOp::Skip(5),
            CigarOp::Padding(1),
            CigarOp::SeqMatch(1),
            CigarOp::SeqMismatch(1),
        ];
        // S1 + M1 + I1 + =1 + X1 = 5 read bases.
        read.seq = b"ACGTA".to_vec();
        read.bq_baq = vec![30; 5];
        assert_eq!(read.length(), Ok(5));
    }

    #[test]
    fn a_seq_and_quality_of_different_lengths_is_an_error() {
        let mut read = consistent_read();
        read.bq_baq = vec![30, 30, 30];
        assert_eq!(
            read.length(),
            Err(ReadLengthError::SeqBqMismatch {
                seq_len: 4,
                bq_baq_len: 3,
            })
        );
    }

    /// The check the cursor's slice indexing depends on: a CIGAR that consumes
    /// a different number of read bases than `seq` holds is caught once, here,
    /// rather than panicking deep inside the walk.
    #[test]
    fn a_cigar_disagreeing_with_seq_is_an_error() {
        let mut read = consistent_read();
        read.cigar = vec![CigarOp::Match(6)];
        assert_eq!(
            read.length(),
            Err(ReadLengthError::CigarSeqMismatch {
                cigar_consumed: 6,
                seq_len: 4,
            })
        );
    }

    /// The `seq`/`bq_baq` check runs first: a read failing both invariants
    /// reports the length mismatch, not the CIGAR one. Pinned because the
    /// walker's error message differs between the two and a reordered copy
    /// would change what a malformed read reports.
    #[test]
    fn the_length_check_precedes_the_cigar_check() {
        let mut read = consistent_read();
        read.bq_baq = vec![30, 30];
        read.cigar = vec![CigarOp::Match(9)];
        assert!(matches!(
            read.length(),
            Err(ReadLengthError::SeqBqMismatch { .. })
        ));
    }

    #[test]
    fn the_mate_role_predicates_are_productions() {
        assert!(!MateRole::Solo.is_paired());
        assert!(MateRole::FirstOfPair.is_paired());
        assert!(MateRole::SecondOfPair.is_paired());

        assert!(!MateRole::Solo.is_first_of_pair());
        assert!(MateRole::FirstOfPair.is_first_of_pair());
        assert!(!MateRole::SecondOfPair.is_first_of_pair());
    }

    /// **The transcription, checked against its original.** `length()` was copied out of
    /// production's `walker/mod.rs` by hand — the one thing in this milestone that was,
    /// everything else being a byte copy — and the tests above pin ng against ng. This
    /// pins ng against production's recorded answers, over every op class and the
    /// CIGAR-length failure — the quality-length failure is the next test's — so a slip
    /// in the op classification or in the order of the two checks fails here rather than
    /// surviving as a silent divergence.
    ///
    /// Production's answers were recorded at commit `d9e7b076` (promotion step C26) and are
    /// the third column. Its `ReadLengthError` was a distinct type with the same variants and
    /// fields as ng's, so its errors are recorded, and compared, through `Debug`.
    #[test]
    fn length_agrees_with_productions_on_every_op_mix() {
        // (cigar, seq_len, production's answer) — the last row is a failure mode.
        let cases: Vec<(Vec<CigarOp>, usize, Result<u64, &str>)> = vec![
            (vec![], 0, Ok(0)),
            (vec![CigarOp::Match(4)], 4, Ok(4)),
            (vec![CigarOp::Deletion(4)], 0, Ok(0)),
            (vec![CigarOp::Skip(4), CigarOp::Match(2)], 2, Ok(2)),
            (
                vec![
                    CigarOp::HardClip(3),
                    CigarOp::SoftClip(2),
                    CigarOp::Match(1),
                ],
                3,
                Ok(3),
            ),
            (vec![CigarOp::Padding(2), CigarOp::Insertion(3)], 3, Ok(3)),
            (
                vec![CigarOp::SeqMatch(2), CigarOp::SeqMismatch(2)],
                4,
                Ok(4),
            ),
            (
                vec![CigarOp::Match(4)],
                5,
                Err("CigarSeqMismatch { cigar_consumed: 4, seq_len: 5 }"),
            ),
        ];
        for (cigar, seq_len, production) in cases {
            let ours = PreparedRead {
                cigar: cigar.clone(),
                seq: vec![b'A'; seq_len],
                bq_baq: vec![30; seq_len],
                ..consistent_read()
            };
            assert_eq!(
                ours.length().map_err(|error| format!("{error:?}")),
                production.map_err(str::to_string),
                "cigar {cigar:?} against {seq_len} bases",
            );
        }
    }

    /// The `seq`/`bq_baq` check runs before the CIGAR check on **both** sides — the one
    /// ordering the table above cannot see, since it never breaks two invariants at once.
    /// Production's answer on this read, recorded at commit `d9e7b076` (promotion step C26),
    /// was the quality-length error.
    #[test]
    fn the_check_order_agrees_with_productions_when_both_invariants_break() {
        let ours = PreparedRead {
            cigar: vec![CigarOp::Match(9)],
            seq: b"ACGT".to_vec(),
            bq_baq: vec![30, 30],
            ..consistent_read()
        };
        assert_eq!(
            ours.length().map_err(|error| format!("{error:?}")),
            Err("SeqBqMismatch { seq_len: 4, bq_baq_len: 2 }".to_string()),
        );
    }
}
