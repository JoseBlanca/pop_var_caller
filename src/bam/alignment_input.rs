//! The vocabulary every decoded read is spelled in: [`CigarOp`], the SAM
//! `FLAG_*` bits, the read filters' `DEFAULT_*` thresholds, [`MappedRead`],
//! and the per-read helpers the filters and the decoder share
//! ([`cigar_to_ops`], [`cigar_ref_span`], [`cigar_is_bad`],
//! [`read_exceeds_mismatch_fraction`], [`compute_adaptor_boundary`]).
//!
//! The reader that uses them is `crate::read::input`; the filter cascade is
//! `crate::read::filtering`.

#[cfg(test)]
use crate::fasta::ContigEntry;
use noodles_sam as sam;

/// One step of a read's CIGAR: what the aligner did there, and for how many bases.
///
/// A read's placement against the reference is spelled out as a run of these — so many bases
/// lined up, so many inserted, so many deleted — and the length each carries counts
/// *reference* bases or *read* bases depending on which step it is. `Match`, `SeqMatch` and
/// `SeqMismatch` consume both; `Insertion` and `SoftClip` consume read only; `Deletion` and
/// `Skip` consume reference only; `HardClip` and `Padding` consume neither.
///
/// The set is SAM's, so htslib's and noodles' too. **`Match` is SAM's `M`, which means
/// *aligned*, not *equal*** — a mismatched base is still `M`, and most aligners emit nothing
/// else. `SeqMatch` and `SeqMismatch` are SAM's `=` and `X`, which do distinguish the two.
///
/// **It lives here because this module is where a CIGAR becomes one**: [`cigar_to_ops`] builds
/// the run from a record, [`MappedRead::cigar`] holds it, and [`cigar_ref_span`] and
/// [`cigar_is_bad`] read it back. Consumers downstream — the pileup walk, the left-aligner,
/// ng's read filters and alignment algorithms — all name the type this decoder produces.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// SAM `M`: aligned to the reference, matching or not. Consumes both.
    Match(u32),
    /// SAM `I`: bases the read has and the reference does not. Consumes read only.
    Insertion(u32),
    /// SAM `D`: bases the reference has and the read does not. Consumes reference only.
    Deletion(u32),
    /// SAM `N`: reference the read skips over — an intron, in spliced alignment. Consumes
    /// reference only, and differs from [`Deletion`](Self::Deletion) in what it means rather
    /// than in what it consumes.
    Skip(u32),
    /// SAM `S`: read bases left out of the alignment but still carried in the record.
    /// Consumes read only.
    SoftClip(u32),
    /// SAM `H`: read bases the aligner cut from the record altogether. Consumes neither.
    HardClip(u32),
    /// SAM `P`: padding against a multiple alignment. Consumes neither.
    Padding(u32),
    /// SAM `=`: aligned *and* equal to the reference. Consumes both.
    SeqMatch(u32),
    /// SAM `X`: aligned and *not* equal to the reference. Consumes both.
    SeqMismatch(u32),
}

// ---------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------

/// Reads with MAPQ strictly below this are dropped. Matches bcftools'
/// default
///
/// Reads with MAPQ unavailable (SAM 0xFF / noodles `mapping_quality()`
/// returning `None`) are treated as MAPQ 0 and therefore rejected
/// under any non-zero minimum. Matches the bcftools/freebayes
/// convention.
pub const DEFAULT_MIN_MAPQ: u8 = 20;

/// Default for `ReadFilterConfig::max_read_mismatch_fraction` (`crate::read::filtering`).
/// 10% — comfortably outside the ~0.5-1% mismatch rate of normal
/// Illumina at MAPQ 20+, inside the regime where adapter-runthrough,
/// contamination, and chimeric tails live. See finding `F1` in
/// `ia/reviews/pileup_freebayes_comparison_2026-05-08.md`.
pub const DEFAULT_MAX_READ_MISMATCH_FRACTION: f32 = 0.10;

/// Default for `ReadFilterConfig::mismatch_bq_floor`. Matches
/// freebayes' `BQL2`. Mismatches whose raw base quality is below this
/// floor do not count toward the mismatch fraction — keeps the filter
/// from firing on genuinely low-quality bases (which the downstream
/// likelihood already de-emphasises). `0` disables the floor entirely.
pub const DEFAULT_MISMATCH_BQ_FLOOR: u8 = 10;

/// Decoded SEQ length below this is dropped. Reads shorter than this
/// rarely contribute reliable alignments at the project's coverage
/// targets.
pub const DEFAULT_MIN_READ_LENGTH: u32 = 30;

// ---------------------------------------------------------------------
// MappedRead
// ---------------------------------------------------------------------

/// A read decoded out of an alignment file (CRAM or BAM). Every
/// byte (qname, sequence, qualities, CIGAR) is owned by the
/// `MappedRead` itself rather than borrowed from the noodles
/// record it was decoded from. That lets downstream stages hold
/// onto reads without being tied to the source reader's lifetime.
///
/// Every field is public — downstream stages (BAQ, the pileup walker)
/// need to read them directly.
///
/// A borrowing design would save one allocation per read at decode
/// time, but each downstream stage would then have to copy the bytes
/// for itself before using them, which costs more in total than
/// allocating once upfront and handing the read along.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MappedRead {
    pub qname: Vec<u8>,
    pub flag: u16,
    /// Index into the canonical `ContigList`.
    pub ref_id: usize,
    /// 1-based leftmost mapped position.
    pub pos: u64,
    pub mapq: u8,
    pub cigar: Vec<CigarOp>,
    /// Uppercase ACGTN.
    pub seq: Vec<u8>,
    /// Phred 0-93, raw BQ before BAQ capping.
    pub qual: Vec<u8>,
    pub mate_ref_id: Option<usize>,
    pub mate_pos: Option<u64>,
    /// Reference position of the first base on this read that lies
    /// inside the mate-pair adaptor (1-based, inclusive on the read's
    /// 3′ side; for reverse-strand reads, inclusive on the 5′ side).
    /// `None` when the boundary cannot be reliably computed (single-end,
    /// mate unmapped, mates on different contigs, TLEN=0, mates on the
    /// same strand, geometry inconsistent, or the molecule is at least
    /// as long as the read so no readthrough is possible). See finding
    /// `G1` in `ia/reviews/pileup_gatk_comparison_2026-05-08.md`.
    pub adaptor_boundary: Option<u32>,
    pub source_file_index: usize,
}

// SAM/BAM flag bit constants — used both inside this module and by
// tests building synthetic records.
pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE_STRAND: u16 = 0x10;
pub const FLAG_MATE_REVERSE_STRAND: u16 = 0x20;
pub const FLAG_FIRST_OF_PAIR: u16 = 0x40;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_QC_FAIL: u16 = 0x200;
pub const FLAG_DUPLICATE: u16 = 0x400;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

// ---------------------------------------------------------------------
// Per-record conversion
// ---------------------------------------------------------------------

/// Finding `G1` (adaptor-region per-base filter) — pure-logic
/// boundary computation. Returns the 1-based reference position of
/// the first base that lies *inside* the mate-pair adaptor, or
/// `None` when the read is single-end, the mate is unreliable, the
/// fragment geometry is inconsistent, or the molecule is at least
/// as long as the read (so no read base could have been sequenced
/// past the molecule end).
///
/// Mirrors GATK's
/// [`ReadUtils.getAdaptorBoundary`](../../gatk/src/main/java/org/broadinstitute/hellbender/utils/read/ReadUtils.java#L527)
/// with one principled deviation: instead of GATK's hardcoded
/// `|tlen| > DEFAULT_ADAPTOR_SIZE = 100` cap (calibrated for 100bp
/// reads), we gate on `|tlen| < seq_len`. The molecule has to be
/// shorter than the *read's own* sequenced length for adaptor
/// readthrough to be possible — that condition is read-aware and
/// stays correct for ancient-DNA short fragments and for modern
/// 150bp+ reads alike.
///
/// On the forward strand the boundary is `read.start + |tlen|`
/// (first base past the molecule's 3′ end). On the reverse strand
/// it is `mate.start - 1` (last base before the molecule's 5′ end);
/// any read base at or before this position is inside the adaptor.
/// The walker's emit sites apply the direction-aware test using
/// `is_reverse_strand`.
pub(crate) fn compute_adaptor_boundary(
    flag: u16,
    ref_id: usize,
    pos: u64,
    seq_len: u32,
    mate_ref_id: Option<usize>,
    mate_pos: Option<u64>,
    template_length: i32,
) -> Option<u32> {
    // Required preconditions, in order of cheapness:
    if flag & FLAG_PAIRED == 0 {
        return None;
    }
    if flag & FLAG_UNMAPPED != 0 || flag & FLAG_MATE_UNMAPPED != 0 {
        return None;
    }
    if template_length == 0 {
        return None;
    }
    let is_reverse = flag & FLAG_REVERSE_STRAND != 0;
    let mate_is_reverse = flag & FLAG_MATE_REVERSE_STRAND != 0;
    // Same-strand pair → not a normal FR/RF pair, geometry is
    // not the simple readthrough case.
    if is_reverse == mate_is_reverse {
        return None;
    }
    let mate_ref_id = mate_ref_id?;
    if mate_ref_id != ref_id {
        return None;
    }
    let mate_pos = mate_pos?;

    let abs_tlen = template_length.unsigned_abs();
    // Per-read gate: only check when the read could physically have
    // run off the end of the molecule. If the molecule is at least
    // as long as the read's sequence, no base sequenced from this
    // read is in adaptor.
    if abs_tlen >= seq_len {
        return None;
    }

    if is_reverse {
        // Reverse-strand read: molecule's 5′ end is at the mate's
        // 1-based start; the first base inside adaptor (on the
        // read's 5′ side, which is the *higher* read offset since
        // the read is reversed) is the position immediately before
        // the mate's start. Skip when ref_pos <= boundary.
        if mate_pos == 0 {
            // Defensive: 1-based positions are non-zero in valid
            // SAM, but a corrupt record could place mate at 0.
            return None;
        }
        Some(mate_pos as u32 - 1)
    } else {
        // Forward-strand read: molecule's 3′ end is at
        // `read.start + |tlen|`. Skip when ref_pos >= boundary.
        Some(pos as u32 + abs_tlen)
    }
}

/// Reference-span sum of `M`-, `=`-, `X`-, `D`-, and `N`-op lengths
/// — i.e. the count of reference positions the read covers. Helper
/// for the F1 mismatch-fraction filter, which needs the slice
/// `[pos, pos + ref_span)` of the reference to compare against.
pub(crate) fn cigar_ref_span(cigar: &[CigarOp]) -> u32 {
    cigar
        .iter()
        .map(|op| match *op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Deletion(n)
            | CigarOp::Skip(n) => n,
            _ => 0,
        })
        .sum()
}

/// Pure-logic implementation of finding `G2` (`GoodCigar`-style
/// read-level rejection). Returns `true` when the CIGAR is
/// "ill-formed" by either of two rules — both sentinel-grade:
///
/// 1. **Adjacent `I`/`D` pair, in either order.** No biological
///    event produces an immediate insertion-then-deletion (or
///    deletion-then-insertion) in a single read; when the aligner
///    emits one, the alignment is genuinely confused and the read's
///    other M-events near that region are also unreliable.
/// 2. **Starts or ends with a deletion** (after stripping leading
///    soft- or hard-clips). A deletion at a CIGAR boundary has no
///    flanking evidence on the missing side; the aligner ran out of
///    read bases at the wrong moment.
///
/// The boundary-deletion rule must be applied to the **original**
/// CIGAR — i.e. before F3 left-alignment runs. F3 deliberately
/// shifts indels through homopolymers/tandem repeats and can in
/// principle land an indel at the read's edge as part of normal
/// canonicalisation; rejecting that case would punish F3 for doing
/// its job. Pre-F3 boundary deletions are the case where the
/// aligner itself produced the suspect alignment.
///
/// Rule 1 (consecutive `I`/`D`) targets the **aligner's** original
/// CIGAR — an adjacent I/D pair as emitted is the suspect signal. It
/// must run before left-alignment: the GATK port left-alignment can
/// *itself* merge two colliding indels into a canonical adjacent
/// `D`/`I` pair, which is a legitimate (parsimonious) representation,
/// not the aligner artefact this rule rejects. Running rule 1 first
/// checks only what the aligner produced.
pub(crate) fn cigar_is_bad(cigar: &[CigarOp]) -> bool {
    // Rule 1 — adjacent I/D in either order.
    for window in cigar.windows(2) {
        let (a, b) = (window[0], window[1]);
        let a_is_indel = matches!(a, CigarOp::Insertion(_) | CigarOp::Deletion(_));
        let b_is_indel = matches!(b, CigarOp::Insertion(_) | CigarOp::Deletion(_));
        // Adjacent same-kind indels (II, DD) shouldn't occur in
        // canonical CIGARs either, but they are not the GATK
        // GoodCigar pattern — only mixed I/D pairs trigger this rule.
        let a_is_ins = matches!(a, CigarOp::Insertion(_));
        let b_is_ins = matches!(b, CigarOp::Insertion(_));
        if a_is_indel && b_is_indel && a_is_ins != b_is_ins {
            return true;
        }
    }

    // Rule 2 — first or last op (after stripping clips) is a
    // deletion. `iter().find` past the leading clips, and `iter()
    // .rev().find` past the trailing clips, give us the first/last
    // *non-clip* op. A read consisting entirely of clips has no
    // such op and is not a boundary-deletion case.
    let first_non_clip = cigar
        .iter()
        .find(|op| !matches!(op, CigarOp::SoftClip(_) | CigarOp::HardClip(_)));
    if matches!(first_non_clip, Some(CigarOp::Deletion(_))) {
        return true;
    }
    let last_non_clip = cigar
        .iter()
        .rev()
        .find(|op| !matches!(op, CigarOp::SoftClip(_) | CigarOp::HardClip(_)));
    if matches!(last_non_clip, Some(CigarOp::Deletion(_))) {
        return true;
    }

    false
}

/// Pure-logic implementation of finding `F1` (per-read
/// mismatch-fraction filter). Returns `true` when the read should
/// be dropped — i.e. its `M`-op mismatch fraction strictly exceeds
/// `threshold`.
///
/// - **Numerator:** count of `M`-op (or `=`/`X`-op) positions where
///   the read base differs from the reference base **and** the raw
///   base quality clears `bq_floor`. Positions where either base is
///   `N` (or non-ATGC) are skipped from both numerator and
///   denominator — they carry no usable signal.
/// - **Denominator:** count of `M`-op positions counted in the
///   numerator's denominator (i.e. ATGC-on-both-sides).
/// - If the denominator is `0` (e.g. a read consisting entirely of
///   soft-clips, or all `M` positions hit `N` bases), the function
///   returns `false` — there is no signal on which to base a drop
///   decision, and the conservative choice is to keep the read and
///   let downstream filters speak.
///
/// `seq`, `qual`, and `cigar` are the read's owned fields after
/// decoding; `ref_seq` is the reference slice covering
/// `[read.pos, read.pos + cigar_ref_span(cigar))`. Indexing into
/// these slices uses the standard CIGAR semantics — see
/// `src/locus_generation/pileup/decompose.rs` for the
/// reference walk pattern this function mirrors.
pub(crate) fn read_exceeds_mismatch_fraction(
    cigar: &[CigarOp],
    seq: &[u8],
    qual: &[u8],
    ref_seq: &[u8],
    bq_floor: u8,
    threshold: f32,
) -> bool {
    let mut read_pos: usize = 0;
    let mut ref_pos: usize = 0;
    let mut mismatches: u32 = 0;
    let mut comparable_bases: u32 = 0;

    for op in cigar {
        match *op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                let n = n as usize;

                // **A position past the end of the read or of the reference contributes to
                // neither count, so the loop stops there rather than stepping through it.**
                // Both used to be read as `N`, and `N` is in neither the read's ATGC set nor
                // the reference's, so such a position could never reach either counter.
                // Cutting the range at the shorter of the two therefore leaves both counts
                // exactly as they were — and lets the body below index three slices taken
                // once per CIGAR run instead of bounds-checking each of them per base.
                // **The starts are clamped as well as the lengths, and that is not
                // belt-and-braces.** A malformed record whose CIGAR claims more read bases
                // than it carries leaves `read_pos` past the end of `seq`, and Rust rejects
                // `&seq[read_pos..read_pos]` there even though the range is empty — so
                // clamping only the length turns a record the old loop merely counted
                // nothing for into a panic that stops the run. Where `usable` is non-zero the
                // clamps are identities, so nothing else sees them.
                let read_from = read_pos.min(seq.len());
                let ref_from = ref_pos.min(ref_seq.len());
                let score_from = read_pos.min(qual.len());
                let usable = n.min(seq.len() - read_from).min(ref_seq.len() - ref_from);
                // The quality scores are the one slice a short tail of does **not** make the
                // position uncountable: an absent score was read as 0, and a position scored
                // 0 is still comparable. So the zipped loop is cut at the qualities too and
                // the remainder finished after it at quality 0 — a tail that is empty for
                // every record whose bases and qualities are the same length.
                let scored = usable.min(qual.len() - score_from);

                let bases = &seq[read_from..read_from + scored];
                let reference = &ref_seq[ref_from..ref_from + scored];
                let scores = &qual[score_from..score_from + scored];

                for ((&base, &reference_base), &score) in bases.iter().zip(reference).zip(scores) {
                    // **`& 0xdf` is exactly `to_ascii_uppercase` on the bytes that matter,
                    // and it needs no branch.** Clearing bit five maps `a`→`A` … `t`→`T`,
                    // and the only bytes it maps *into* `ACGT` are those eight — so testing
                    // the folded byte accepts precisely the eight-case set this arm used to
                    // test the reference against, and the two `matches!` become symmetric.
                    let folded = reference_base & 0xdf;
                    let read_is_acgt = matches!(base, b'A' | b'C' | b'G' | b'T');
                    let reference_is_acgt = matches!(folded, b'A' | b'C' | b'G' | b'T');

                    // Counted rather than branched on, so the body is straight-line.
                    let comparable = read_is_acgt & reference_is_acgt;
                    comparable_bases += u32::from(comparable);
                    mismatches += u32::from(comparable & (base != folded) & (score >= bq_floor));
                }

                // The positions past the end of the quality scores. An absent score was read
                // as 0, and 0 clears the floor exactly when the floor is 0 — which is the
                // whole of `score >= bq_floor` here, hoisted out of the loop.
                let absent_score_clears_the_floor = bq_floor == 0;
                for k in scored..usable {
                    let base = seq[read_from + k];
                    let folded = ref_seq[ref_from + k] & 0xdf;
                    let comparable = matches!(base, b'A' | b'C' | b'G' | b'T')
                        & matches!(folded, b'A' | b'C' | b'G' | b'T');
                    comparable_bases += u32::from(comparable);
                    mismatches +=
                        u32::from(comparable & (base != folded) & absent_score_clears_the_floor);
                }

                read_pos += n;
                ref_pos += n;
            }
            CigarOp::Insertion(n) => {
                read_pos += n as usize;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                ref_pos += n as usize;
            }
            CigarOp::SoftClip(n) => {
                read_pos += n as usize;
            }
            CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }

    if comparable_bases == 0 {
        return false;
    }
    let fraction = (mismatches as f32) / (comparable_bases as f32);
    fraction > threshold
}

/// `pub(crate)` for ng's own decode, which reuses this rather than copying the
/// CIGAR translation (`src/read/aligned_read.rs`).
pub(crate) fn cigar_to_ops(cigar: &sam::alignment::record_buf::Cigar) -> Vec<CigarOp> {
    use sam::alignment::record::cigar::op::Kind;
    cigar
        .as_ref()
        .iter()
        .map(|op| {
            let len = op.len() as u32;
            match op.kind() {
                Kind::Match => CigarOp::Match(len),
                Kind::Insertion => CigarOp::Insertion(len),
                Kind::Deletion => CigarOp::Deletion(len),
                Kind::Skip => CigarOp::Skip(len),
                Kind::SoftClip => CigarOp::SoftClip(len),
                Kind::HardClip => CigarOp::HardClip(len),
                Kind::Pad => CigarOp::Padding(len),
                Kind::SequenceMatch => CigarOp::SeqMatch(len),
                Kind::SequenceMismatch => CigarOp::SeqMismatch(len),
            }
        })
        .collect()
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // --- Pure-type tests ---------------------------------------------

    #[test]
    fn p1_contig_list_md5_wildcard_equality() {
        let with_md5 = ContigEntry {
            name: "chr1".into(),
            length: 100,
            md5: Some([1u8; 16]),
        };
        let same_md5 = with_md5.clone();
        let no_md5 = ContigEntry {
            name: "chr1".into(),
            length: 100,
            md5: None,
        };
        let different_md5 = ContigEntry {
            name: "chr1".into(),
            length: 100,
            md5: Some([2u8; 16]),
        };
        let different_name = ContigEntry {
            name: "chr2".into(),
            length: 100,
            md5: None,
        };
        let different_length = ContigEntry {
            name: "chr1".into(),
            length: 200,
            md5: None,
        };

        assert_eq!(with_md5, same_md5);
        assert_eq!(with_md5, no_md5);
        assert_eq!(no_md5, with_md5);
        assert_ne!(with_md5, different_md5);
        assert_ne!(with_md5, different_name);
        assert_ne!(with_md5, different_length);
    }

    /// **The rewritten loop against the one it replaced, on 20,000 random reads.**
    ///
    /// This filter is shared with the frozen production caller
    /// (`src/pileup/per_sample/read_processor.rs`), so the loop that stopped bounds-checking
    /// three slices per base had to be *exactly* the loop it replaced rather than close to
    /// it. Two of its steps are arguments rather than rewrites — that a position past the end
    /// of the read or of the reference can reach neither counter, and that `& 0xdf` accepts
    /// the same eight reference bytes the two-case list did — and an argument is what a
    /// randomised comparison is for.
    ///
    /// The generator aims at exactly the cases the arguments turn on: sequences, qualities
    /// and reference deliberately allowed to run short of what the CIGAR claims; every base
    /// drawn from `ACGTNacgtn` so both cases and the excluded byte appear; qualities
    /// straddling the floor, and a floor of 0 so that the short-quality tail's `0 >= bq_floor`
    /// is exercised in both directions.
    #[test]
    fn the_rewritten_mismatch_loop_agrees_with_the_one_it_replaced() {
        /// The loop exactly as it stood before the rewrite, kept here as the oracle and
        /// nowhere else.
        fn as_it_was(
            cigar: &[CigarOp],
            seq: &[u8],
            qual: &[u8],
            ref_seq: &[u8],
            bq_floor: u8,
            threshold: f32,
        ) -> bool {
            let mut read_pos: usize = 0;
            let mut ref_pos: usize = 0;
            let mut mismatches: u32 = 0;
            let mut comparable_bases: u32 = 0;
            for op in cigar {
                match *op {
                    CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                        let n = n as usize;
                        for k in 0..n {
                            let r = seq.get(read_pos + k).copied().unwrap_or(b'N');
                            let g = ref_seq.get(ref_pos + k).copied().unwrap_or(b'N');
                            let q = qual.get(read_pos + k).copied().unwrap_or(0);
                            let r_atgc = matches!(r, b'A' | b'C' | b'G' | b'T');
                            let g_atgc =
                                matches!(g, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't');
                            if r_atgc && g_atgc {
                                comparable_bases += 1;
                                if r != g.to_ascii_uppercase() && q >= bq_floor {
                                    mismatches += 1;
                                }
                            }
                        }
                        read_pos += n;
                        ref_pos += n;
                    }
                    CigarOp::Insertion(n) | CigarOp::SoftClip(n) => read_pos += n as usize,
                    CigarOp::Deletion(n) | CigarOp::Skip(n) => ref_pos += n as usize,
                    CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
                }
            }
            if comparable_bases == 0 {
                return false;
            }
            (mismatches as f32) / (comparable_bases as f32) > threshold
        }

        // A fixed generator rather than a crate: the sequence is the same on every machine
        // and every run, so a disagreement can be reproduced from the seed alone.
        // The seed is the date this fixture was written, 2026-09-08, so that it reads as a
        // seed rather than as a number somebody tuned. Grouped in fours because that is what
        // `clippy::unusual_byte_groupings` accepts, which still leaves the year and the
        // month-day legible either side of the underscore.
        let mut state: u64 = 0x2026_0908_u64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        const ALPHABET: &[u8] = b"ACGTNacgtn";

        let mut compared_something = 0usize;
        for _ in 0..20_000 {
            let ops = 1 + (next() % 4) as usize;
            let mut cigar = Vec::with_capacity(ops);
            let mut read_span = 0usize;
            let mut ref_span = 0usize;
            for _ in 0..ops {
                let n = 1 + (next() % 12) as u32;
                match next() % 7 {
                    0 => {
                        cigar.push(CigarOp::Insertion(n));
                        read_span += n as usize;
                    }
                    1 => {
                        cigar.push(CigarOp::Deletion(n));
                        ref_span += n as usize;
                    }
                    2 => {
                        cigar.push(CigarOp::SoftClip(n));
                        read_span += n as usize;
                    }
                    3 => cigar.push(CigarOp::HardClip(n)),
                    op => {
                        cigar.push(match op {
                            4 => CigarOp::SeqMatch(n),
                            5 => CigarOp::SeqMismatch(n),
                            _ => CigarOp::Match(n),
                        });
                        read_span += n as usize;
                        ref_span += n as usize;
                    }
                }
            }

            // Deliberately allowed to fall short of the span the CIGAR claims — a
            // well-formed record never does, and the truncation argument is about exactly
            // the records that do.
            let short_by = |limit: usize, draw: u64| limit.saturating_sub((draw % 4) as usize);
            let seq: Vec<u8> = (0..short_by(read_span, next()))
                .map(|_| ALPHABET[(next() % ALPHABET.len() as u64) as usize])
                .collect();
            let qual: Vec<u8> = (0..short_by(read_span, next()))
                .map(|_| (next() % 60) as u8)
                .collect();
            let reference: Vec<u8> = (0..short_by(ref_span, next()))
                .map(|_| ALPHABET[(next() % ALPHABET.len() as u64) as usize])
                .collect();
            // 0 among the floors so the short-quality tail's `0 >= bq_floor` is true
            // sometimes and false sometimes.
            let bq_floor = [0u8, 1, 20, 30][(next() % 4) as usize];
            let threshold = [0.0f32, 0.05, 0.10, 0.5][(next() % 4) as usize];

            if !seq.is_empty() && !reference.is_empty() {
                compared_something += 1;
            }
            assert_eq!(
                read_exceeds_mismatch_fraction(
                    &cigar, &seq, &qual, &reference, bq_floor, threshold
                ),
                as_it_was(&cigar, &seq, &qual, &reference, bq_floor, threshold),
                "the two loops disagreed on cigar {cigar:?}, seq {seq:?}, qual {qual:?}, \
                 reference {reference:?}, floor {bq_floor}, threshold {threshold}",
            );
        }
        // A guard on the generator, not on the loops: a CIGAR drawn as all clips, insertions
        // and deletions gives a read or a reference of nothing, and a run of those would let
        // this test pass while comparing two functions that both returned `false`. It lands at
        // about 13,900 of 20,000.
        assert!(
            compared_something > 10_000,
            "the generator produced too many empty cases to be comparing anything: only \
             {compared_something} of 20,000 had both a read and a reference",
        );
    }

    // --- F1 mismatch-fraction helper: pure-logic tests ---------------
    //
    // These exercise `read_exceeds_mismatch_fraction` directly, with
    // synthetic inputs (no Repository, no CRAM). Integration of the
    // filter into the alignment_input pipeline is covered by the
    // `default_config_*` tests below.

    fn m_only(len: u32) -> Vec<CigarOp> {
        vec![CigarOp::Match(len)]
    }

    #[test]
    fn f1_zero_mismatches_passes() {
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"ACGTACGTAC";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_high_mismatch_drops() {
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"AGGGAGGTAG"; // 4 mismatches in 10 → 40 %
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_threshold_boundary_is_strict_greater_than() {
        // Exactly at threshold: should NOT drop (filter is `>`, not `>=`).
        // 1 mismatch in 10 = 0.10 fraction; threshold = 0.10 → keep.
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"AGGTACGTAC"; // 1 mismatch at pos 1
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
        // Just above threshold: 2/10 = 0.20.
        let ref_seq_2mm = b"AGGAACGTAC"; // 2 mismatches
        assert!(read_exceeds_mismatch_fraction(
            &cigar,
            seq,
            &qual,
            ref_seq_2mm,
            10,
            0.10
        ));
    }

    #[test]
    fn f1_bq_floor_excludes_low_quality_mismatches() {
        // 5 mismatches in 10 bases (= 50 %), but all with BQ=5 (< floor=10).
        // Filter should ignore them all → keep.
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![5u8; 10];
        let ref_seq = b"AGGGAGGGAC"; // 5 mismatches
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
        // Same input with BQ=20 → mismatches count → drop.
        let qual_high = vec![20u8; 10];
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual_high, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_bq_floor_zero_counts_every_mismatch() {
        // floor=0 → every mismatch counts regardless of BQ.
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![0u8; 10];
        let ref_seq = b"AGGGAGGGAC"; // 5 mismatches
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 0, 0.10
        ));
    }

    #[test]
    fn f1_n_in_ref_skipped_from_both_numerator_and_denominator() {
        // Read = ACGTACGTAC, Ref = NNNNNNNNNN. All M positions hit
        // N-in-ref → comparable_bases = 0 → keep (no signal).
        let cigar = m_only(10);
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b"NNNNNNNNNN";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_n_in_read_skipped_from_both_numerator_and_denominator() {
        // 5 read bases are N (skipped); the other 5 all match
        // → 0 mismatches in 5 comparable → keep.
        let cigar = m_only(10);
        let seq = b"NNNNNACGTA";
        let qual = vec![30u8; 10];
        let ref_seq = b"AAAAAACGTA";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_soft_clips_excluded_from_denominator() {
        // 5S5M: only the 5 M-op bases count. Two of them mismatch
        // → 2/5 = 40 % → drop. The freebayes formula (with seqlen
        // denominator = 10) would give 2/10 = 20 % — still drop here
        // but this test pins the M-op-only behaviour we documented.
        let cigar = vec![CigarOp::SoftClip(5), CigarOp::Match(5)];
        let seq = b"AAAAAACGTA"; // first 5 = soft clip
        let qual = vec![30u8; 10];
        let ref_seq = b"GGGTA"; // ref covers only the 5 M bases
        assert!(read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_indels_advance_cursors_but_dont_count() {
        // 4M 2I 4M with all matches in M-ops and an insertion of
        // garbage in the middle. M-op bases = 8, mismatches = 0.
        // Insertion does not contribute to either count → keep.
        let cigar = vec![CigarOp::Match(4), CigarOp::Insertion(2), CigarOp::Match(4)];
        let seq = b"ACGTNNACGT"; // last 4 align after the 2-base ins
        let qual = vec![30u8; 10];
        let ref_seq = b"ACGTACGT";
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_no_m_op_returns_false() {
        // All-soft-clip read (shouldn't happen post-MAPQ-filter but
        // could). comparable_bases = 0 → keep (conservative).
        let cigar = vec![CigarOp::SoftClip(10)];
        let seq = b"ACGTACGTAC";
        let qual = vec![30u8; 10];
        let ref_seq = b""; // ref_span = 0
        assert!(!read_exceeds_mismatch_fraction(
            &cigar, seq, &qual, ref_seq, 10, 0.10
        ));
    }

    #[test]
    fn f1_default_config_has_filter_enabled_at_ten_percent() {
        // Pin the defaults: threshold 0.10, floor 10. That the read filter's
        // default config uses them is pinned in `crate::read::filtering`.
        assert_eq!(DEFAULT_MAX_READ_MISMATCH_FRACTION, 0.10);
        assert_eq!(DEFAULT_MISMATCH_BQ_FLOOR, 10);
    }

    // --- G1 adaptor-boundary helper: pure-logic tests ----------------
    //
    // These exercise `compute_adaptor_boundary` directly, with synthetic
    // SAM-flag/TLEN/mate-pos inputs. Its use in decoding a record is
    // covered by `crate::read::aligned_read`'s tests.

    /// Build the standard "properly paired, FR orientation" flag for a
    /// forward-strand read whose mate is reverse-strand. Tests pass
    /// this directly so the SAM bit fiddling is centralised.
    const FLAG_FR_FWD: u16 = FLAG_PAIRED | FLAG_MATE_REVERSE_STRAND;
    /// Same as `FLAG_FR_FWD` but for the reverse-strand mate.
    const FLAG_FR_REV: u16 = FLAG_PAIRED | FLAG_REVERSE_STRAND;

    #[test]
    fn g1_single_end_returns_none() {
        // No FLAG_PAIRED → boundary is undefined.
        assert_eq!(
            compute_adaptor_boundary(0, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_mate_unmapped_returns_none() {
        let flag = FLAG_FR_FWD | FLAG_MATE_UNMAPPED;
        assert_eq!(
            compute_adaptor_boundary(flag, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_self_unmapped_returns_none() {
        let flag = FLAG_FR_FWD | FLAG_UNMAPPED;
        assert_eq!(
            compute_adaptor_boundary(flag, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_zero_tlen_returns_none() {
        // TLEN=0 means the aligner could not determine the insert size,
        // typically because the mates are on different references or
        // the placement is ambiguous. Don't filter.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(200), 0),
            None
        );
    }

    #[test]
    fn g1_same_strand_pair_returns_none() {
        // Same-strand FF or RR: improperly oriented, geometry is not
        // the simple readthrough case.
        let same_strand_flag = FLAG_PAIRED; // both forward (no FLAG_MATE_REVERSE_STRAND)
        assert_eq!(
            compute_adaptor_boundary(same_strand_flag, 0, 100, 150, Some(0), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_mate_on_different_contig_returns_none() {
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(1), Some(200), 130),
            None
        );
    }

    #[test]
    fn g1_molecule_at_least_as_long_as_read_returns_none() {
        // 150bp read in a 200bp insert: fragment ≥ read, no readthrough.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(250), 200),
            None
        );
        // Edge case: |tlen| == seq_len, also no readthrough.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(250), 150),
            None
        );
    }

    #[test]
    fn g1_forward_strand_short_insert_boundary_is_start_plus_tlen() {
        // 150bp read, insert 50bp: read overruns by 100bp into the
        // 3′ adaptor. Boundary at read.start + |tlen| = 100 + 50 = 150.
        // Any base at ref_pos >= 150 is in adaptor.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(149), 50),
            Some(150)
        );
    }

    #[test]
    fn g1_reverse_strand_short_insert_boundary_is_mate_start_minus_one() {
        // The reverse-strand read in an FR pair: mate's start is the
        // forward-strand mate's leftmost ref pos (the molecule's 5′
        // end on the forward strand). Boundary = mate.start - 1.
        // ancient-DNA shape: 50bp molecule, mate starts at 100, our
        // reverse read end placed beyond.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_REV, 0, 200, 150, Some(0), Some(100), -50),
            Some(99)
        );
    }

    #[test]
    fn g1_negative_tlen_takes_absolute_value() {
        // Forward-strand mate that happens to be the second-in-pair
        // (TLEN convention is reverse-mate TLEN, forward-mate -TLEN).
        // |tlen| should be used so a negative sign does not break
        // the calculation.
        assert_eq!(
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 100, 150, Some(0), Some(149), -50),
            Some(150)
        );
    }

    #[test]
    fn g1_ancient_dna_short_fragment_yields_boundary_inside_read_span() {
        // 70bp ancient-DNA molecule sequenced with 100bp reads:
        // mate starts at 1000, |tlen| = 70.
        // Forward mate (start 1000): boundary at 1000 + 70 = 1070.
        // Read covers ref [1000, 1100), so positions 1070..1100
        // are filtered (30 of 100 bases).
        let boundary_fwd =
            compute_adaptor_boundary(FLAG_FR_FWD, 0, 1000, 100, Some(0), Some(1030), 70);
        assert_eq!(boundary_fwd, Some(1070));
        // Reverse mate (the same molecule's other end): boundary
        // at mate.start - 1 = 1000 - 1 = 999. The reverse mate
        // sits at ref [1030, 1130), but its leftward overrun is
        // anywhere ≤ 999 — outside its own placement, which means
        // its leftmost soft-clipped/`M` bases past 999 are flagged.
        let boundary_rev =
            compute_adaptor_boundary(FLAG_FR_REV, 0, 1030, 100, Some(0), Some(1000), -70);
        assert_eq!(boundary_rev, Some(999));
    }

    #[test]
    fn g1_default_record_has_no_boundary() {
        // A record without any of the paired flags / TLEN must yield
        // None — the no-op case the cursor relies on.
        assert_eq!(
            compute_adaptor_boundary(0, 0, 100, 150, None, None, 0),
            None
        );
    }

    // --- G2 cigar_is_bad helper: pure-logic tests --------------------
    //
    // These exercise `cigar_is_bad` directly. Its use in the read
    // filter cascade is covered by `crate::read::filtering`'s tests.

    #[test]
    fn g2_well_formed_cigar_is_not_bad() {
        // Standard shapes that should pass: pure-match, single
        // interior indel, soft-clipped ends, multi-indel CIGARs
        // separated by M ops.
        assert!(!cigar_is_bad(&[CigarOp::Match(100)]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Insertion(2),
            CigarOp::Match(50),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Deletion(3),
            CigarOp::Match(50),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::SoftClip(5),
            CigarOp::Match(90),
            CigarOp::SoftClip(5),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Match(40),
            CigarOp::Deletion(1),
            CigarOp::Match(20),
        ]));
    }

    #[test]
    fn g2_consecutive_id_is_bad() {
        // I followed by D.
        assert!(cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Deletion(1),
            CigarOp::Match(20),
        ]));
        // D followed by I (the symmetric case).
        assert!(cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Deletion(1),
            CigarOp::Insertion(2),
            CigarOp::Match(20),
        ]));
    }

    #[test]
    fn g2_consecutive_same_kind_indels_are_not_bad() {
        // II or DD is non-canonical but not the GoodCigar pattern —
        // those should be merged or normalised by the aligner; we
        // don't reject on them here because the GATK rule is only
        // about *mixed* I/D pairs.
        assert!(!cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Insertion(2),
            CigarOp::Insertion(1),
            CigarOp::Match(20),
        ]));
        assert!(!cigar_is_bad(&[
            CigarOp::Match(20),
            CigarOp::Deletion(2),
            CigarOp::Deletion(1),
            CigarOp::Match(20),
        ]));
    }

    #[test]
    fn g2_first_op_deletion_is_bad() {
        assert!(cigar_is_bad(&[CigarOp::Deletion(2), CigarOp::Match(50),]));
    }

    #[test]
    fn g2_last_op_deletion_is_bad() {
        assert!(cigar_is_bad(&[CigarOp::Match(50), CigarOp::Deletion(2),]));
    }

    #[test]
    fn g2_first_op_deletion_after_clip_is_bad() {
        // GATK's rule explicitly says "with or without preceding
        // clips": a deletion after only soft/hard clips is still a
        // boundary deletion.
        assert!(cigar_is_bad(&[
            CigarOp::SoftClip(5),
            CigarOp::Deletion(2),
            CigarOp::Match(50),
        ]));
        assert!(cigar_is_bad(&[
            CigarOp::HardClip(5),
            CigarOp::SoftClip(3),
            CigarOp::Deletion(2),
            CigarOp::Match(50),
        ]));
    }

    #[test]
    fn g2_last_op_deletion_after_trailing_clip_is_bad() {
        assert!(cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Deletion(2),
            CigarOp::SoftClip(5),
        ]));
        assert!(cigar_is_bad(&[
            CigarOp::Match(50),
            CigarOp::Deletion(2),
            CigarOp::SoftClip(3),
            CigarOp::HardClip(2),
        ]));
    }

    #[test]
    fn g2_first_op_insertion_is_not_bad() {
        // Boundary insertions are not the GoodCigar rule. The walker
        // already drops first/last-op insertions as events (no
        // flanking evidence), but the read itself stays in the active
        // set and contributes Match events from the surrounding ops.
        assert!(!cigar_is_bad(&[CigarOp::Insertion(2), CigarOp::Match(50),]));
        assert!(!cigar_is_bad(&[CigarOp::Match(50), CigarOp::Insertion(2),]));
        assert!(!cigar_is_bad(&[
            CigarOp::SoftClip(3),
            CigarOp::Insertion(2),
            CigarOp::Match(50),
        ]));
    }

    #[test]
    fn g2_only_clips_is_not_bad() {
        // Pathological but not the GoodCigar pattern: a CIGAR of just
        // clips has no first/last *non-clip* op, so the boundary-
        // deletion check finds no candidate and returns false.
        // Reads like this are dropped earlier by the min-read-length
        // / unmapped-flag checks, so this is just a defensive case.
        assert!(!cigar_is_bad(&[CigarOp::SoftClip(50)]));
        assert!(!cigar_is_bad(&[
            CigarOp::HardClip(10),
            CigarOp::SoftClip(40),
        ]));
    }

    #[test]
    fn g2_empty_cigar_is_not_bad() {
        // An empty CIGAR has no ops at all — `windows(2)` yields
        // nothing, `find` returns None, neither rule fires. The
        // upstream record decode would normally reject empty CIGARs
        // before this point; this test pins the no-panic behaviour.
        assert!(!cigar_is_bad(&[]));
    }
}
