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
//! transcribed rather than re-derived. The construction path is production's
//! too — [`PreparedRead::from_production`] takes the output of
//! `prepare_passthrough` and re-attaches the group — so ng computes none of the
//! per-field wiring itself and the read-preparation parity fixture stays exact.
//!
//! Design: `doc/devel/ng/spec/locus_generation_pileup.md` §6,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*.

use std::sync::Arc;

use crate::ng::types::ReadGroupId;
// Production's, reused rather than re-minted — the same call, for the same reason, as
// `ng::alignment`'s (see its "# Reusing production's `CigarOp`"). The import path
// misleads: `CigarOp` is crate-wide CIGAR vocabulary that happens to live in the walker.
use crate::pileup::walker::CigarOp;
// Aliased, so "ours" and "production's" read at a glance instead of being carried by a
// four-segment path at every one of the dozen sites below.
use crate::pileup::walker::MateRole as ProductionMateRole;
use crate::pileup::walker::PreparedRead as ProductionPreparedRead;

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

    /// The same role, on ng's copy of the enum.
    ///
    /// Written as an exhaustive `match` rather than a cast, so a variant added to
    /// **production's** enum stops this compiling instead of being silently mapped to a
    /// neighbour. *That is the only direction the compiler checks here:* a variant added
    /// to **ng's** would simply be unreachable from any production read, in silence —
    /// [`to_production`](Self::to_production) is the `#[cfg(test)]` witness that closes
    /// it.
    ///
    /// An inherent `pub(crate)` function rather than a `From` impl, matching its
    /// neighbour [`PreparedRead::from_production`]: trait impls carry no visibility, so a
    /// `From` would publish an ng→production coupling that nothing outside the crate can
    /// use.
    pub(crate) fn from_production(role: ProductionMateRole) -> Self {
        match role {
            ProductionMateRole::Solo => MateRole::Solo,
            ProductionMateRole::FirstOfPair => MateRole::FirstOfPair,
            ProductionMateRole::SecondOfPair => MateRole::SecondOfPair,
        }
    }

    /// The reverse direction, **`#[cfg(test)]` only**: nothing needs it at run time.
    /// See also [`PreparedRead::into_production`], which the walker parity harness needs
    /// for the same reason — to hand one prepared stream to both walkers.
    ///
    /// It exists so that an exhaustive `match` over *ng's* enum also has to be updated
    /// when ng's enum grows — which is what makes
    /// [`from_production`](Self::from_production)'s claim true in both directions instead
    /// of half-true.
    #[cfg(test)]
    pub(crate) fn to_production(self) -> ProductionMateRole {
        match self {
            MateRole::Solo => ProductionMateRole::Solo,
            MateRole::FirstOfPair => ProductionMateRole::FirstOfPair,
            MateRole::SecondOfPair => ProductionMateRole::SecondOfPair,
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

    /// Production's prepared read plus the read group it had nowhere to put.
    ///
    /// **The per-field wiring stays production's.** Read preparation builds its
    /// output by calling `prepare_passthrough` — production's `--no-baq` arm —
    /// rather than re-deriving `alignment_end`, `mate_role`, `mq_log_err` and
    /// the rest; this is the one line where that output becomes ng's type. Every
    /// field is moved across unchanged, and `read_group` is the only thing added.
    ///
    /// Written as a **destructure** so a field added to production's type stops
    /// this compiling instead of being silently dropped on the way into ng's —
    /// the same reason `AlignedRead::into_mapped_read` destructures on the way
    /// out.
    pub(crate) fn from_production(read: ProductionPreparedRead, read_group: ReadGroupId) -> Self {
        let ProductionPreparedRead {
            chrom_id,
            alignment_start,
            alignment_end,
            cigar,
            seq,
            bq_baq,
            mq_log_err,
            mapq,
            is_reverse_strand,
            qname,
            mate_role,
            adaptor_boundary,
        } = read;

        Self {
            chrom_id,
            alignment_start,
            alignment_end,
            cigar,
            seq,
            bq_baq,
            mq_log_err,
            mapq,
            is_reverse_strand,
            qname,
            mate_role: MateRole::from_production(mate_role),
            adaptor_boundary,
            read_group,
        }
    }

    /// ng's prepared read as production's, **dropping the read group** — the one field
    /// production has nowhere to put.
    ///
    /// `#[cfg(test)]`, and it exists for exactly one caller: the walker parity harness,
    /// which must hand **one** prepared stream to both walkers. Preparing twice would
    /// compare two different inputs and call the result parity. Reads are prepared once,
    /// by ng's preparer, and converted down here.
    ///
    /// Destructured for the reason `from_production` is: a field added to ng's type must
    /// stop this compiling, so that whoever adds it decides whether it has a production
    /// counterpart rather than discovering later that the parity harness quietly stopped
    /// carrying it.
    #[cfg(test)]
    pub(crate) fn into_production(self) -> ProductionPreparedRead {
        let Self {
            chrom_id,
            alignment_start,
            alignment_end,
            cigar,
            seq,
            bq_baq,
            mq_log_err,
            mapq,
            is_reverse_strand,
            qname,
            mate_role,
            adaptor_boundary,
            // Production's type has no counterpart; this is the whole reason ng owns its
            // own read. The walk does not read it, which is what makes dropping it safe
            // *here* and only here.
            read_group: _,
        } = self;

        ProductionPreparedRead {
            chrom_id,
            alignment_start,
            alignment_end,
            cigar,
            seq,
            bq_baq,
            mq_log_err,
            mapq,
            is_reverse_strand,
            qname,
            mate_role: mate_role.to_production(),
            adaptor_boundary,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::read::PLACEHOLDER_READ_GROUP;

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

    /// Each production role maps to its own counterpart — not to a neighbour.
    /// A conversion that collapsed two roles would silently disable the
    /// mate-overlap tie-break rather than fail.
    ///
    /// Both directions, because the two catch different mistakes: `from_production`'s
    /// `match` is what the compiler checks when *production* grows a variant, and
    /// `to_production` is what it checks when *ng* does. Asserting the round trip on top
    /// is what rules out a pair of conversions that are each exhaustive and disagree.
    #[test]
    fn every_production_mate_role_maps_to_its_counterpart() {
        for (theirs, ours) in [
            (ProductionMateRole::Solo, MateRole::Solo),
            (ProductionMateRole::FirstOfPair, MateRole::FirstOfPair),
            (ProductionMateRole::SecondOfPair, MateRole::SecondOfPair),
        ] {
            assert_eq!(MateRole::from_production(theirs), ours);
            assert_eq!(ours.to_production(), theirs);
        }
    }

    /// **The transcription, checked against its original.** `length()` was copied out of
    /// production's `walker/mod.rs` by hand — the one thing in this milestone that was,
    /// everything else being a byte copy — and the tests above pin ng against ng. This
    /// pins ng against production, over every op class and both failure modes, so a slip
    /// in the op classification or in the order of the two checks fails here rather than
    /// surviving as a silent divergence.
    ///
    /// The two `ReadLengthError` types are distinct but structurally identical, so the
    /// comparison goes through `Debug`.
    #[test]
    fn length_agrees_with_productions_on_every_op_mix() {
        // (cigar, seq_len) — the last two rows are the two failure modes.
        let cases: Vec<(Vec<CigarOp>, usize)> = vec![
            (vec![], 0),
            (vec![CigarOp::Match(4)], 4),
            (vec![CigarOp::Deletion(4)], 0),
            (vec![CigarOp::Skip(4), CigarOp::Match(2)], 2),
            (
                vec![
                    CigarOp::HardClip(3),
                    CigarOp::SoftClip(2),
                    CigarOp::Match(1),
                ],
                3,
            ),
            (vec![CigarOp::Padding(2), CigarOp::Insertion(3)], 3),
            (vec![CigarOp::SeqMatch(2), CigarOp::SeqMismatch(2)], 4),
            (vec![CigarOp::Match(4)], 5),
        ];
        for (cigar, seq_len) in cases {
            let production = ProductionPreparedRead {
                chrom_id: 0,
                alignment_start: 1,
                alignment_end: 1,
                cigar: cigar.clone(),
                seq: vec![b'A'; seq_len],
                bq_baq: vec![30; seq_len],
                mq_log_err: -6.0,
                mapq: 60,
                is_reverse_strand: false,
                qname: Arc::from("read1"),
                mate_role: ProductionMateRole::Solo,
                adaptor_boundary: None,
            };
            let ours = PreparedRead::from_production(production.clone(), ReadGroupId(1));
            assert_eq!(
                ours.length().map_err(|error| format!("{error:?}")),
                production.length().map_err(|error| format!("{error:?}")),
                "cigar {cigar:?} against {seq_len} bases",
            );
        }
    }

    /// The `seq`/`bq_baq` check runs before the CIGAR check on **both** sides — the one
    /// ordering the table above cannot see, since it never breaks two invariants at once.
    #[test]
    fn the_check_order_agrees_with_productions_when_both_invariants_break() {
        let production = ProductionPreparedRead {
            chrom_id: 0,
            alignment_start: 1,
            alignment_end: 1,
            cigar: vec![CigarOp::Match(9)],
            seq: b"ACGT".to_vec(),
            bq_baq: vec![30, 30],
            mq_log_err: -6.0,
            mapq: 60,
            is_reverse_strand: false,
            qname: Arc::from("read1"),
            mate_role: ProductionMateRole::Solo,
            adaptor_boundary: None,
        };
        let ours = PreparedRead::from_production(production.clone(), ReadGroupId(1));
        assert_eq!(
            ours.length().map_err(|error| format!("{error:?}")),
            production.length().map_err(|error| format!("{error:?}")),
        );
    }

    /// **The reverse conversion moves every field too.** The mirror of the test below, and
    /// the one that needs saying out loud: `into_production`'s only caller is the
    /// `#[ignore]`d walker parity harness, so without this it is compiled and never
    /// executed — zero runtime coverage on the conversion the entire real-data half of the
    /// oracle rests on. A lossy one would surface as "the two walkers disagree", which is
    /// the wrong diagnosis and would be chased into the walker.
    #[test]
    fn into_production_moves_every_field_to_its_counterpart() {
        let ours = PreparedRead {
            chrom_id: 2,
            alignment_start: 101,
            alignment_end: 140,
            cigar: vec![CigarOp::SoftClip(1), CigarOp::Match(3)],
            seq: b"TACG".to_vec(),
            bq_baq: vec![11, 22, 33, 44],
            mq_log_err: -13.5,
            mapq: 37,
            is_reverse_strand: true,
            qname: Arc::from("frag/1"),
            mate_role: MateRole::SecondOfPair,
            adaptor_boundary: Some(137),
            read_group: ReadGroupId(4),
        };

        let theirs = ours.clone().into_production();

        assert_eq!(theirs.chrom_id, 2);
        assert_eq!(theirs.alignment_start, 101);
        assert_eq!(theirs.alignment_end, 140);
        assert_eq!(theirs.cigar, ours.cigar);
        assert_eq!(theirs.seq, b"TACG".to_vec());
        assert_eq!(theirs.bq_baq, vec![11, 22, 33, 44]);
        assert_eq!(theirs.mq_log_err, -13.5);
        assert_eq!(theirs.mapq, 37);
        assert!(theirs.is_reverse_strand);
        assert_eq!(&*theirs.qname, "frag/1");
        assert_eq!(theirs.mate_role, ProductionMateRole::SecondOfPair);
        assert_eq!(theirs.adaptor_boundary, Some(137));
    }

    /// **The claim the real-data differential rests on: the round trip is lossless for
    /// everything the walk reads.**
    ///
    /// One prepared stream is handed to both walkers by converting ng's *down* to
    /// production's, so a field that does not survive the trip means the two walkers saw
    /// different inputs — and the divergence would be blamed on the walk. The two
    /// conversions are individually plausible and could still disagree; this is what rules
    /// that out, as `every_production_mate_role_maps_to_its_counterpart` does for the role.
    ///
    /// `read_group` is the one field *meant* to be lost, so the round trip re-attaches it
    /// and everything else must match.
    #[test]
    fn the_two_conversions_round_trip_every_field_the_walk_reads() {
        for (role, boundary) in [
            (ProductionMateRole::Solo, None),
            (ProductionMateRole::FirstOfPair, Some(137)),
            (ProductionMateRole::SecondOfPair, Some(1)),
        ] {
            let original = ProductionPreparedRead {
                chrom_id: 2,
                alignment_start: 101,
                alignment_end: 140,
                cigar: vec![CigarOp::SoftClip(1), CigarOp::Match(3)],
                seq: b"TACG".to_vec(),
                bq_baq: vec![11, 22, 33, 44],
                mq_log_err: -13.5,
                mapq: 37,
                is_reverse_strand: true,
                qname: Arc::from("frag/1"),
                mate_role: role,
                adaptor_boundary: boundary,
            };
            let back =
                PreparedRead::from_production(original.clone(), ReadGroupId(4)).into_production();
            // Production's type has no `PartialEq`, so `Debug` is the oracle — total over
            // these fields, and it names the one that moved.
            assert_eq!(
                format!("{back:?}"),
                format!("{original:?}"),
                "the round trip lost or altered a field on {role:?}",
            );
        }
    }

    /// **The conversion moves every field to its counterpart, not to a
    /// neighbour of the same type.** The destructure makes it exhaustive; only
    /// this makes it correct. The fixture gives each field a distinct value for
    /// that reason — `chrom_id` and `alignment_start` are both `u32`, and
    /// swapping them would compile.
    #[test]
    fn the_conversion_from_productions_read_moves_every_field() {
        let production = ProductionPreparedRead {
            chrom_id: 2,
            alignment_start: 101,
            alignment_end: 140,
            cigar: vec![CigarOp::SoftClip(1), CigarOp::Match(3)],
            seq: b"TACG".to_vec(),
            bq_baq: vec![11, 22, 33, 44],
            mq_log_err: -13.5,
            mapq: 37,
            is_reverse_strand: true,
            qname: Arc::from("frag/1"),
            mate_role: ProductionMateRole::SecondOfPair,
            adaptor_boundary: Some(137),
        };

        let ours = PreparedRead::from_production(production.clone(), ReadGroupId(4));

        assert_eq!(ours.chrom_id, production.chrom_id);
        assert_eq!(ours.alignment_start, production.alignment_start);
        assert_eq!(ours.alignment_end, production.alignment_end);
        assert_eq!(ours.cigar, production.cigar);
        assert_eq!(ours.seq, production.seq);
        assert_eq!(ours.bq_baq, production.bq_baq);
        assert_eq!(ours.mq_log_err, production.mq_log_err);
        assert_eq!(ours.mapq, production.mapq);
        assert_eq!(ours.is_reverse_strand, production.is_reverse_strand);
        assert_eq!(ours.qname, production.qname);
        assert_eq!(ours.mate_role, MateRole::SecondOfPair);
        assert_eq!(ours.adaptor_boundary, production.adaptor_boundary);
        assert_eq!(
            ours.read_group,
            ReadGroupId(4),
            "the one field production has no counterpart for"
        );
    }
}
