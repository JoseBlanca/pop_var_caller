//! Classify one read against an STR locus: delimit its repeat with the generator's
//! [`RepeatDelimiter`], recover a long allele by widening the read region, gate the tract's base
//! quality, and map the result to an observation or to the reason there is none.
//!
//! The per-read pipeline follows the deleted caller's `classify_read`, over the new delimiters.
//! Its result, [`Classified`], is what [`tally`](super::tally) folds.

use super::read_region::{
    Footprint, MIN_REGION_Q1, extract_region, flank_truncated, passes_quality_gate, read_footprint,
    tract_span_from_input, widen_region,
};
use super::{MeasuredGround, RepeatDelimiter, SsrLocus};
use crate::alignment::{ReadBases, RepeatContext, RepeatGeometry, RepeatSpan, StutterModel};
use crate::bam::alignment_input::CigarOp;
use crate::locus_generation::{FlankDeletion, HeldFlank, LocusLen, ReadWitness};
use crate::read::aligned_read::AlignedRead;
use crate::types::Bp;
use std::ops::Range;

/// Why a read yielded no usable observation — the tally increments the matching
/// `SsrGeneratorCounts` reason (and `reads_without_observation`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum NoObservationReason {
    /// Neither flank anchored (the read lies wholly inside the repeat), or a malformed read.
    NoBorderAnchored,
    /// The delimited tract failed the base-quality gate.
    LowQuality,
    /// The allele stayed flank-truncated even after widening.
    WindowTruncated,
    /// **The read is not in this locus.** It anchored a flank and crossed no tract
    /// position at all, so the "lower bound" it would supply is *at least zero* — no
    /// evidence about the repeat's length, on a read that never entered the repeat.
    ///
    /// It arises because the fetch queries the tract **plus its margin**, so a read
    /// overlapping only the flank is delimited too; the aligner then reports
    /// [`RepeatSpan::FromLeft`] or
    /// `FromRight` with an empty span, the rejected side having fallen back to the
    /// read's own edge. Measured on tomato chr01 of `SRR7279503`: **6,704 such reads
    /// against 7,085 genuine partials**, at a median window overlap of 16 bases against
    /// a 30-base flank.
    ///
    /// Those bases belong to the SNP/indel path, which analyses them; the STR path
    /// discards them and counts them (owner, 2026-07-30). Before that decision they
    /// became observations with **empty bases** and a witness covering zero locus
    /// positions — 3,180 rows of the STR dump — contributing to `num_obs` while adding
    /// nothing to depth anywhere along the locus.
    OutsideTract,
}

/// What one read contributes to a locus: an observed tract (complete or partial), or a
/// reason it observed nothing.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum Classified {
    Observed {
        bases: Box<[u8]>,
        read_witness: ReadWitness,
        /// The flank bases the read's deletion took beside the tract, when the read is
        /// spelled from its own alignment (`str_edge_deletion_alleles.md` §5.1); `NONE`
        /// otherwise.
        flank_deletion: FlankDeletion,
        /// Which flank a partial read anchored, where its witness cannot say it: set at an
        /// STR segment or bundle, on a witness covering the whole locus (`partial`); `None`
        /// otherwise.
        held_flank: Option<HeldFlank>,
        /// The read's base-quality error mass over the tract, `Σ ln(P_err)`, in log-error
        /// space (freebayes' `q_sum` convention). The **BQ** support moment the spec names
        /// (spec §3, "strand/BQ/MAPQ moments") — computed here because the tract base
        /// qualities are already sliced, a free by-product; the tally folds it into
        /// [`SequenceObservation::q_sum`](crate::locus_generation::SequenceObservation). MAPQ
        /// is carried separately, off the [`AlignedRead`], in the tally. **Soft** — filled,
        /// unconsumed today.
        q_sum: f64,
    },
    NoObservation(NoObservationReason),
}

/// Classify one read against `locus` using `aligner` (the chosen [`RepeatDelimiter`]). `stutter`
/// feeds the context — algorithm 3 (flat-gap) ignores it, algorithm 4 (unit-slip) prices its
/// slips from it; `align_scratch` (the aligner's own scratch) and `qual_buffer` are reused
/// across reads.
///
/// The pipeline mirrors production: extract the read's slice over the locus window, delimit
/// it, and — if a *complete* tract looks truncated by the window (a mapper-collapsed long
/// allele) — widen the slice and re-delimit, giving up as `WindowTruncated` if it stays
/// truncated. A complete tract is quality-gated; **partials are kept as lower bounds**
/// without the gate (spec §3, the new behaviour production discards).
///
/// **A complete observation's bases come from the read's input alignment where it can
/// serve** ([`tract_span_from_input`]), and from the delimiter's own span only as the
/// fallback for reads the input alignment does not carry across the tract (clips, long
/// alleles). The delimiter still decides completeness, anchoring and the widen-and-retry;
/// what it must not decide is the spelling, because its stutter-priced re-alignment
/// re-spells junction variation toward a pure motif run (tract-accuracy program, L4).
pub(super) fn classify_read<A: RepeatDelimiter>(
    read: &AlignedRead,
    locus: &SsrLocus,
    aligner: &A,
    stutter: &StutterModel,
    align_scratch: &mut A::Scratch,
    qual_buffer: &mut Vec<u8>,
) -> Classified {
    // A malformed record (qual length ≠ seq length) cannot be sliced safely.
    if read.qual.len() != read.seq.len() {
        return Classified::NoObservation(NoObservationReason::NoBorderAnchored);
    }
    let read_len = read.seq.len();
    let reference = &*locus.tract_with_margin_bases;
    let left = locus.left_flank_len();
    let right = locus.right_flank_len();
    let geometry = RepeatGeometry {
        left_flank_len: Bp(left as u64),
        right_flank_len: Bp(right as u64),
        motif: locus.repeat.delimiter_motif,
    };
    let context = RepeatContext {
        geometry: &geometry,
        stutter,
    };

    let fp = read_footprint(&read.cigar, read.pos);
    let window_start = (locus.margin_start.get() - 1) as u32; // 0-based
    let window_len = reference.len() as u32;

    // **A read whose own alignment never reaches the tract measures nothing about it**
    // (owner, 2026-10-10). The fetch takes the tract plus its margin, so reads lying wholly
    // in a flank arrive here, and the delimiter cannot always tell them from a read that
    // crossed: where a flank's last bases repeat the other flank's first — `CCACCA` on both
    // sides of a 6-base repeat on the stuttering benchmark — a read stopping exactly at the
    // tract matched its last six bases to the far flank and was measured as a complete
    // observation of an **empty** tract. At 50 samples and 15× that was 122 reads at 34
    // segments, and an empty candidate stopped the locus. The read's CIGAR settles it before
    // any alignment: ending at or before the tract's first base, or starting at or after its
    // last, with no unaligned bases on the side facing it — a soft clip or an insertion
    // there could be tract the mapper did not place, and goes to the delimiter as before.
    let tract_start = window_start + left as u32;
    let tract_end = tract_start + locus.repeat.tract_len() as u32;
    let tract_on_reference = tract_start..tract_end;
    let (unaligned_before, unaligned_after) = unaligned_bases_at_the_ends(&read.cigar);
    let ends_before_the_tract = fp.ref_end <= tract_start && unaligned_after == 0;
    let starts_after_the_tract = fp.ref_start >= tract_end && unaligned_before == 0;
    if ends_before_the_tract || starts_after_the_tract {
        return Classified::NoObservation(NoObservationReason::OutsideTract);
    }

    let region = extract_region(&read.cigar, fp, read_len, window_start, window_len);

    // **A read that copies the reference across the whole window has nothing to align.**
    // Its repeat count is the reference's, and the tract sits where the reference's tract
    // sits — so the delimiter would spend a full dynamic-programming pass rediscovering a
    // span this function already knows. On tomato at three reads a position that is 4 reads
    // in every 9 reaching the aligner (57,507 of 132,069 over 24 spans of chromosome 1).
    //
    // The base-quality gate below still runs: skipping the aligner skips *where the tract
    // is*, not *whether the read is good enough to be counted*, which is a separate question
    // and is what turns the span into an observation.
    if let Some(tract) = tract_needing_no_alignment(read, &region, locus, reference) {
        let tract = to_usize(&tract);
        return complete_or_low_quality(
            read,
            region.start + tract.start..region.start + tract.end,
            FlankDeletion::NONE,
            qual_buffer,
        );
    }

    // **The read's own input alignment, mapped over the tract's reference positions —
    // but only where that alignment actually says something.** A CIGAR with an
    // insertion or deletion is the mapper's genome-wide claim about where the read's
    // length variation sits, and it is information the delimiter's window-local,
    // stutter-priced re-alignment destroys at junctions (see `tract_span_from_input`).
    // An all-Match CIGAR claims nothing the bases don't: it is exactly how a mapper
    // spells a collapsed long allele — the length change hidden as flank mismatches —
    // and rescuing that read is what the delimiter and its widen-and-retry are *for*,
    // so an indel-free read is left entirely to them.
    let cigar: &[CigarOp] = &read.cigar;
    let input_span = if left > 0
        && right > 0
        && cigar
            .iter()
            .any(|op| matches!(op, CigarOp::Insertion(_) | CigarOp::Deletion(_)))
    {
        tract_span_from_input(
            cigar,
            &read.seq,
            fp,
            window_start + left as u32,
            locus.repeat.tract_len() as u32,
            reference[left - 1],
            reference[left + locus.repeat.tract_len() as usize - 1],
        )
    } else {
        None
    };

    // **A read whose own alignment carries a deletion across the tract's edge, and is
    // aligned on both sides beyond the tract, crossed the tract** (owner, 2026-10-10;
    // spec §3). The delimiter calls a read complete only where both flanks match intact,
    // and a deletion that takes the first bases of a flank breaks that match, so every read
    // carrying it came back a partial and the deleted allele was never a candidate: at 50
    // samples and 15× on the stuttering benchmark, 33 of the 62 such alleles the cohort
    // carries, 563 carrier genotypes, none right. The read's own alignment already says it
    // crossed, and is what spells it (above). Measured: +1,122 right genotypes at 50 samples
    // over four runs for 16 more false calls in 321,116, +76 at one sample for none, and
    // +138 on HG002's tandem-repeat benchmark at 30× and 50×, where one tract a depth went
    // from right to wrong.
    //
    // **And the flank bases its deletion took are recorded beside the spelling** (owner,
    // 2026-10-10; `str_edge_deletion_alleles.md` §5.1), so the whole deletion is one allele
    // of the tract's record rather than a tract piece here and a flank piece on the generic
    // path. A deletion reaching past the flank this locus holds is not recorded, since the
    // record could not be widened over it.
    //
    // **In its smallest form** (`str_edge_deletion_alleles.md` §4): the mapper chooses where a
    // deletion sits in a run of equal bases, so the same sequence can come as different
    // splits between spelling and flank deletion, and a split that can give flank bases back
    // does. A deletion that fits wholly inside the tract — a stutter contraction, whichever
    // way the mapper or a left-aligner placed it — comes out as no flank deletion at all.
    let edge_deleted = flank_bases_deleted_at_the_edges(cigar, fp.ref_start, &tract_on_reference);
    let flank_deletion = flank_deletion_of(edge_deleted)
        .filter(|deleted| deleted.left as usize <= left && deleted.right as usize <= right)
        .unwrap_or(FlankDeletion::NONE);
    let (input_span, flank_deletion) = match input_span {
        Some(span) if !flank_deletion.is_none() => {
            let tract_len = locus.repeat.tract_len() as usize;
            let (span, flank_deletion) = through_the_flank_indels(
                cigar,
                fp,
                span,
                flank_deletion,
                &tract_on_reference,
                guarded_flank_widths(left, right),
            );
            let (span, deleted) = in_smallest_form(
                &read.seq,
                span,
                flank_deletion,
                (
                    &reference[left - flank_deletion.left as usize..left],
                    &reference[left + tract_len..left + tract_len + flank_deletion.right as usize],
                ),
                locus.repeat.delimiter_motif.period(),
            );
            (Some(span), deleted)
        }
        other => (other, FlankDeletion::NONE),
    };
    if let Some(span) = &input_span {
        let (aligned_before, aligned_after) =
            aligned_bases_beyond_the_tract(cigar, fp.ref_start, &tract_on_reference);
        if edge_deleted != (0, 0)
            && aligned_before >= MIN_ALIGNED_FLANK_BASES
            && aligned_after >= MIN_ALIGNED_FLANK_BASES
        {
            return complete_or_low_quality(read, span.clone(), flank_deletion, qual_buffer);
        }
    }

    let span = delimit(aligner, read, &region, reference, context, align_scratch);
    if locus.repeat.ground == MeasuredGround::BundleHull
        && let Some((spelled, deleted)) = flank_deletion_the_delimiter_implies(
            read,
            &span,
            &region,
            fp,
            window_start..window_start + window_len,
            tract_start..tract_end,
            reference,
            (left, right),
        )
    {
        return complete_or_low_quality(read, spelled, deleted, qual_buffer);
    }
    match span {
        RepeatSpan::Between(tract)
            if flank_truncated(&region, &to_usize(&tract), read_len, left, right) =>
        {
            // A complete tract whose flanks look eaten by the window: widen and retry.
            // **The delimiter's own spelling, always** — a read in this arm is the
            // mapper-collapsed long allele whose input alignment is exactly the account
            // not to trust.
            let wide = widen_region(region, read_len, left, right);
            let wspan = delimit(aligner, read, &wide, reference, context, align_scratch);
            match wspan {
                RepeatSpan::Between(wtract)
                    if !flank_truncated(&wide, &to_usize(&wtract), read_len, left, right) =>
                {
                    let wtract = to_usize(&wtract);
                    complete_or_low_quality(
                        read,
                        wide.start + wtract.start..wide.start + wtract.end,
                        FlankDeletion::NONE,
                        qual_buffer,
                    )
                }
                _ => Classified::NoObservation(NoObservationReason::WindowTruncated),
            }
        }
        RepeatSpan::Between(tract) => {
            let tract = to_usize(&tract);
            // **The input spelling wins whenever it is available — including when it
            // reads as the reference.** The alternative — deferring reference-spelling
            // inputs to the delimiter — was built and measured (arm `l4_v2`,
            // tract-accuracy program L4): it repaired the handful of loci where a
            // mapper places an ambiguous indel outside the tract span
            // (`chr1:194,220,455`) and re-broke twice as many where the mapper's
            // reference account of the tract is exactly right — at
            // `chr1:104,262,346` every read carries a flank-`C` deletion the truth
            // confirms (`GC→G`, hom), the input span correctly reads eleven reference
            // `T`s, and the delimiter re-absorbs the flank deletion as a spurious
            // one-unit contraction, double-counting the deleted base.
            match input_span {
                Some(spelled) => {
                    complete_or_low_quality(read, spelled, flank_deletion, qual_buffer)
                }
                None => complete_or_low_quality(
                    read,
                    region.start + tract.start..region.start + tract.end,
                    FlankDeletion::NONE,
                    qual_buffer,
                ),
            }
        }
        // **An empty one-sided span is a read outside the locus, not a lower bound of
        // zero** — and `partial` is where that is decided, in one place.
        //
        // This arm used to hold a `if tract.is_empty()` guard of its own, added
        // when `ReadWitness::from_left`/`from_right` were infallible and nothing else
        // could reject the case. They now return an `Option`, so `partial` has to handle
        // "no position covered" anyway and answers exactly this — which made the guard
        // **unfalsifiable**: the Milestone C review deleted it outright and the whole
        // suite stayed green, because no input can tell the two paths apart. Two
        // spellings of one decision, one of them untestable, is worse than one.
        //
        // **Unless the read's own alignment puts an indel in the flank bases the delimiter
        // may not open a gap in** (2026-11-01). The delimiter's junction guard closes the
        // true path for such a read, and where the read's tract also differs from the
        // reference the cheapest path left ends the tract early and lays the read's flank
        // across the reference's with mismatches, which the anchor rule then rejects: at
        // SL4.0ch05:51017179 all ten of one sample's reads came back partial and tied two
        // equal-length candidates, though eight spelled one of them base for base. The
        // read's own alignment opened the gap the guard forbids, so it is what spells it.
        RepeatSpan::FromLeft(_) | RepeatSpan::FromRight(_)
            if input_span.is_some()
                && indel_in_the_guarded_flank(
                    cigar,
                    fp.ref_start,
                    &tract_on_reference,
                    guarded_flank_widths(left, right),
                )
                && {
                    let (before, after) =
                        aligned_bases_beyond_the_tract(cigar, fp.ref_start, &tract_on_reference);
                    before >= MIN_ALIGNED_FLANK_BASES && after >= MIN_ALIGNED_FLANK_BASES
                } =>
        {
            let spelled = input_span.expect("the guard checked it is Some");
            complete_or_low_quality(read, spelled, flank_deletion, qual_buffer)
        }
        RepeatSpan::FromLeft(tract) => partial(read, &region, &tract, HeldFlank::Left, locus),
        RepeatSpan::FromRight(tract) => partial(read, &region, &tract, HeldFlank::Right, locus),
        RepeatSpan::Unanchored => Classified::NoObservation(NoObservationReason::NoBorderAnchored),
    }
}

/// **The flank deletion an STR bundle read's delimiter spelling implies**, as a complete
/// observation in its smallest form — or `None`, and the read is measured as at a segment
/// (found in review, 2026-10-11; `str_bundle_loci.md` §6.1).
///
/// A deletion across a hull's outer edge is a sequence the mapper can write without any
/// deletion crossing the edge: as substitutions (`CTCCGA|AAAC…` read as `CTCAAC|AAAC…`), or
/// as a deletion inside the flank and an insertion at the hull's first base. Then the
/// own-alignment rule above does not fire, and the flat-gap delimiter — which prices nothing
/// inside the hull and a stiff gap in the flanks — either deletes the flank bases and spells
/// the read's extra bases into the hull, or fails to anchor that flank at all. The first
/// wrote the event twice, once as hull bases in the bundle record and once as the flank
/// deletion the generic path calls beside it; the second made every carrier read a partial,
/// which calling leaves out at a bundle. On the stuttering benchmark (50 samples, PCR, 30×)
/// that was 24 of the 51 wrong bundle cells. (A segment's own delimiter prices an
/// insertion inside the repeat that is not whole units, so it keeps such a flank in place
/// and leaves the change to the generic path, once.)
///
/// So at a bundle, **where the delimiter itself says a flank lost bases** — it anchored a
/// complete read with fewer read bases on that side than the flank holds, or it could not
/// anchor that side of a read whose own alignment deletes bases in that flank — the read's
/// bases over the window are read as the flank's longest matching stretch from its far end,
/// the deleted bases between it and the hull, and the spelling. Both flanks must keep at
/// least [`MIN_ALIGNED_FLANK_BASES`] matching bases, the read's own alignment must cover the
/// whole window and reach that far beyond the hull on each side, and the side the delimiter
/// anchored in full keeps its junction. A flank the delimiter left intact never yields a
/// deletion here, so a substituted flank base still reads as it does at a segment
/// (`str_edge_deletion_alleles.md` §4, the rule tried and not adopted).
#[allow(clippy::too_many_arguments)]
fn flank_deletion_the_delimiter_implies(
    read: &AlignedRead,
    span: &RepeatSpan,
    region: &Range<usize>,
    fp: Footprint,
    window: Range<u32>,
    tract: Range<u32>,
    reference: &[u8],
    (left, right): (usize, usize),
) -> Option<(Range<usize>, FlankDeletion)> {
    let region_len = region.end - region.start;
    let cigar: &[CigarOp] = &read.cigar;
    let (deleted_before, deleted_after) =
        reference_bases_deleted_in_the_flanks(cigar, fp.ref_start, &window, &tract);
    // Each side: (whether the delimiter left it short, the junction it kept, region-relative).
    let (left_short, left_junction, right_short, right_junction) = match span {
        RepeatSpan::Between(tract_span) => {
            let tract_span = to_usize(tract_span);
            (
                tract_span.start < left,
                tract_span.start,
                region_len.saturating_sub(tract_span.end) < right,
                tract_span.end,
            )
        }
        RepeatSpan::FromRight(tract_span) => {
            (deleted_before > 0, 0, false, to_usize(tract_span).end)
        }
        RepeatSpan::FromLeft(tract_span) => (
            false,
            to_usize(tract_span).start,
            deleted_after > 0,
            region_len,
        ),
        RepeatSpan::Unanchored => return None,
    };
    if !(left_short || right_short) {
        return None;
    }
    if fp.ref_start > window.start || fp.ref_end < window.end {
        return None;
    }
    let (aligned_before, aligned_after) =
        aligned_bases_beyond_the_tract(cigar, fp.ref_start, &tract);
    if aligned_before < MIN_ALIGNED_FLANK_BASES || aligned_after < MIN_ALIGNED_FLANK_BASES {
        return None;
    }
    let bases = &read.seq[region.clone()];
    let left_flank = &reference[..left];
    let right_flank = &reference[reference.len() - right..];
    let kept_left = bases
        .iter()
        .zip(left_flank)
        .take_while(|(read_base, flank_base)| read_base == flank_base)
        .count();
    let kept_right = bases
        .iter()
        .rev()
        .zip(right_flank.iter().rev())
        .take_while(|(read_base, flank_base)| read_base == flank_base)
        .count();
    let floor = MIN_ALIGNED_FLANK_BASES as usize;
    if (left_short && kept_left < floor) || (right_short && kept_right < floor) {
        return None;
    }
    let (start, deleted_left) = if left_short {
        (kept_left, left - kept_left)
    } else {
        (left_junction, 0)
    };
    let (end, deleted_right) = if right_short {
        (region_len - kept_right, right - kept_right)
    } else {
        (right_junction, 0)
    };
    if start > end {
        return None;
    }
    let deleted = FlankDeletion {
        left: u8::try_from(deleted_left).ok()?,
        right: u8::try_from(deleted_right).ok()?,
    };
    if deleted.is_none() {
        return None;
    }
    let tract_len = reference.len() - left - right;
    Some(in_smallest_form(
        &read.seq,
        region.start + start..region.start + end,
        deleted,
        (
            &left_flank[left - deleted_left..],
            &reference[left + tract_len..left + tract_len + deleted_right],
        ),
        0,
    ))
}

/// **How many reference bases the read's deletions take inside each flank of the window**,
/// as `(left flank, right flank)` — the flanks being `window` less `tract`. Coordinates as in
/// [`flank_bases_deleted_at_the_edges`].
fn reference_bases_deleted_in_the_flanks(
    cigar: &[CigarOp],
    ref_start: u32,
    window: &Range<u32>,
    tract: &Range<u32>,
) -> (u32, u32) {
    let overlap = |from: u32, to: u32, lo: u32, hi: u32| to.min(hi).saturating_sub(from.max(lo));
    let (mut before, mut after, mut at) = (0, 0, ref_start);
    for op in cigar {
        match op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Skip(n) => at += n,
            CigarOp::Deletion(n) => {
                before += overlap(at, at + n, window.start, tract.start);
                after += overlap(at, at + n, tract.end, window.end);
                at += n;
            }
            _ => {}
        }
    }
    (before, after)
}

/// **The fewest aligned bases a read must carry beyond the tract on each side** for its own
/// alignment, deletion and all, to count as having crossed the tract. Five keeps a read
/// whose alignment merely touches a flank on the delimiter's path. It was not swept: it is
/// the value the measurement in [`classify_read`]'s comment was taken at.
const MIN_ALIGNED_FLANK_BASES: u32 = 5;

/// **How many flank bases the read's deletions take at the tract's edges**, as `(before the
/// tract, after it)` — `(0, 0)` when no deletion takes any. A deletion counts on the left when
/// it crosses the tract's first base or ends at it, taking the flank bases before it; on the
/// right when it crosses the tract's last base or starts just after it, taking the flank bases
/// after; one covering the whole tract counts on both sides. Where two deletions count on one
/// side, the later one's count stands. `tract` is the tract's 0-based reference positions,
/// the end exclusive; `ref_start` is 0-based too.
///
/// Classification uses the one count two ways: whether it is `(0, 0)` decides whether a
/// read's own alignment, deletion and all, counts as having crossed the tract, and
/// [`flank_deletion_of`] turns it into the flank deletion recorded beside the spelling.
///
/// **A deletion that touches the edge counts** (`str_edge_deletion_alleles.md` §5.1, amended
/// 2026-10-18). Where a run of one base is followed by a different base, a read that misreads
/// that base as the run's carries a run one base too long; minimap2 writes such a read as the
/// run at full length and the neighbouring base deleted, a deletion that ends at the edge
/// without crossing it. Left out, the read was spelled over the tract's positions alone and
/// counted as the reference length: at HG002's homozygous one-base-deletion homopolymers at
/// 30x, 565 of the 926 reference-length reads were such reads. Recorded, the read keeps its
/// flank deletion, and calling compares it with each allele over the wider extent (§5.6), where
/// it matches the deletion allele but for the misread base.
fn flank_bases_deleted_at_the_edges(
    cigar: &[CigarOp],
    ref_start: u32,
    tract: &Range<u32>,
) -> (u32, u32) {
    let (tract_start, tract_end) = (tract.start, tract.end);
    let (mut left, mut right, mut at) = (0_u32, 0_u32, ref_start);
    for op in cigar {
        match op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Skip(n) => at += n,
            CigarOp::Deletion(n) => {
                let (from, to) = (at, at + n);
                if from < tract_start && to >= tract_start {
                    left = tract_start - from;
                }
                if from <= tract_end && to > tract_end {
                    right = to - tract_end;
                }
                at = to;
            }
            _ => {}
        }
    }
    (left, right)
}

/// **How many flank bases beside each edge the delimiter will not open a gap in**, as `(left,
/// right)`: the default junction guard's flank width, cut to half the flank as the delimiter
/// cuts it ([`JunctionGuard`]'s `flank_width`).
///
/// [`JunctionGuard`]: crate::alignment::ssr_unit_robust::JunctionGuard
fn guarded_flank_widths(left_flank_len: usize, right_flank_len: usize) -> (u32, u32) {
    let width = crate::alignment::ssr_unit_robust::JunctionGuard::default().flank_columns;
    (
        width.min(left_flank_len / 2) as u32,
        width.min(right_flank_len / 2) as u32,
    )
}

/// **A read's flank deletion run on through every indel its alignment has in the guarded flank
/// bases**, with the spelling taking the read's bases up to where its alignment rejoins the
/// reference (2026-11-01).
///
/// A mapper can write one event beside a tract as several: a deletion at the tract's edge, a
/// few matches and an insertion. The deletion alone then made the read's allele, and the
/// rest went to the generic path, which wrote the flank the read actually carries; merged,
/// the two disagreed on the flank bases and the STR record won, so the sample got an allele
/// nobody carried (simulated, `ng_tract_simulator` with flank deletions: one repeat
/// inserted and `GAACT` deleted one base from the tract came back as the tract unchanged and
/// two flank bases deleted). Taken as one event, the allele is the read's own sequence there.
///
/// On each side, the flank deletion grows to the farthest reference base an indel in the
/// guarded bases reaches, and the spelling to the read's bases up to it. A flank that grows
/// past the window's own is left as it was, as is everything when no indel lies there.
fn through_the_flank_indels(
    cigar: &[CigarOp],
    fp: Footprint,
    span: Range<usize>,
    deleted: FlankDeletion,
    tract: &Range<u32>,
    widths: (u32, u32),
) -> (Range<usize>, FlankDeletion) {
    let (tract_start, tract_end) = (tract.start, tract.end);
    // The farthest flank bases an indel in the guard reaches: `first` on the left (the
    // first base it takes or sits before), `last` on the right (one past the last).
    let mut first = tract_start - u32::from(deleted.left);
    let mut last = tract_end + u32::from(deleted.right);
    let left_guard = tract_start.saturating_sub(widths.0);
    let right_guard = tract_end + widths.1;
    let mut at = fp.ref_start;
    for op in cigar {
        match op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Skip(n) => at += n,
            CigarOp::Deletion(n) => {
                let (from, to) = (at, at + n);
                if from < tract_start && to > left_guard {
                    first = first.min(from);
                }
                if to > tract_end && from < right_guard {
                    last = last.max(to);
                }
                at = to;
            }
            CigarOp::Insertion(_) => {
                if at < tract_start && at > left_guard {
                    first = first.min(at);
                }
                if at > tract_end && at <= right_guard {
                    last = last.max(at);
                }
            }
            _ => {}
        }
    }
    let grown = FlankDeletion {
        left: u8::try_from(tract_start - first).unwrap_or(deleted.left),
        right: u8::try_from(last - tract_end).unwrap_or(deleted.right),
    };
    if grown == deleted
        || grown.left as u32 > widths.0.max(u32::from(deleted.left))
        || grown.right as u32 > widths.1.max(u32::from(deleted.right))
    {
        return (span, deleted);
    }
    // The read position before reference base `first` (and any insertion before it), and
    // the one after the insertions before base `last`.
    let (mut start, mut end) = (None, None);
    let (mut ref_cur, mut read_cur) = (fp.ref_start, fp.leading_clip as usize);
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                if start.is_none() && (ref_cur..ref_cur + n).contains(&first) {
                    start = Some(read_cur + (first - ref_cur) as usize);
                }
                if end.is_none() && (ref_cur..ref_cur + n).contains(&last) {
                    end = Some(read_cur + (last - ref_cur) as usize);
                }
                ref_cur += n;
                read_cur += *n as usize;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                if start.is_none() && (ref_cur..ref_cur + n).contains(&first) {
                    start = Some(read_cur);
                }
                if end.is_none() && (ref_cur..ref_cur + n).contains(&last) {
                    end = Some(read_cur);
                }
                ref_cur += n;
            }
            CigarOp::Insertion(n) => {
                if start.is_none() && ref_cur == first {
                    start = Some(read_cur);
                }
                read_cur += *n as usize;
            }
            _ => {}
        }
    }
    match (start, end) {
        (Some(from), Some(to)) if from <= span.start && to >= span.end => (from..to, grown),
        _ => (span, deleted),
    }
}

/// **Does an insertion or deletion in the CIGAR sit in the guarded flank bases** — the
/// `widths.0` reference bases before the tract or the `widths.1` after it? A deletion counts
/// when it takes any of those bases, an insertion when it sits between two of them or at the
/// edge. Coordinates as in [`flank_bases_deleted_at_the_edges`].
fn indel_in_the_guarded_flank(
    cigar: &[CigarOp],
    ref_start: u32,
    tract: &Range<u32>,
    widths: (u32, u32),
) -> bool {
    let (tract_start, tract_end) = (tract.start, tract.end);
    let left_guard = tract_start.saturating_sub(widths.0)..tract_start;
    let right_guard = tract_end..tract_end + widths.1;
    let mut at = ref_start;
    for op in cigar {
        match op {
            CigarOp::Match(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n)
            | CigarOp::Skip(n) => at += n,
            CigarOp::Deletion(n) => {
                let (from, to) = (at, at + n);
                let takes = |guard: &Range<u32>| from < guard.end && to > guard.start;
                if takes(&left_guard) || takes(&right_guard) {
                    return true;
                }
                at = to;
            }
            // An insertion sits before reference base `at`.
            CigarOp::Insertion(_)
                if (at > left_guard.start && at <= left_guard.end)
                    || (at >= right_guard.start && at < right_guard.end) =>
            {
                return true;
            }
            _ => {}
        }
    }
    false
}

/// **The flank deletion [`flank_bases_deleted_at_the_edges`] counted, as a
/// [`FlankDeletion`]**, or `None` when it counted none on either side, or more than 255 on
/// one (which no flank this path holds is long enough for).
fn flank_deletion_of((left, right): (u32, u32)) -> Option<FlankDeletion> {
    let deleted = FlankDeletion {
        left: u8::try_from(left).ok()?,
        right: u8::try_from(right).ok()?,
    };
    (!deleted.is_none()).then_some(deleted)
}

/// **A read's spelling and flank deletion in their smallest form** — the split that gives
/// back every flank base it can (`str_edge_deletion_alleles.md` §4).
///
/// With `left` flank bases missing, the read's sequence in place of the reference from `left`
/// bases before the tract to its end is the spelling. If the spelling begins with the
/// first of those flank bases, the same sequence is the spelling without that base in place
/// of one flank base fewer — the deletion placed one base further right. So while the spelling
/// begins with the first deleted flank base, that base moves back into the flank. The right
/// side is the mirror image: while the spelling ends with the last deleted flank base. The
/// result does not depend on where in a run of equal bases the mapper put the deletion, and a
/// deletion that fits wholly inside the tract — stutter — comes out as none.
///
/// `span` is the spelling in read coordinates; `left_flank` and `right_flank` are the
/// reference bases the deletion took, the tract's side last.
///
/// **No flank base is given back that would leave the spelling shorter than `shortest`**
/// (2026-11-01): a segment's period, since an allele without a whole copy of the motif is one
/// the read model cannot score, and candidate selection drops it; zero at a bundle, whose
/// alleles can hold no copy. At chr20:2,472,054 (HG002) a 9-base deletion took seven flank
/// bases, `AAAAAAG`, and two of an `A` ×8 tract: given back without the floor, every carrier
/// read became an empty spelling with the `G` deleted, the allele was dropped, and only the
/// generic path's shortened deletion was written. With it the reads spell `A` with `AG`
/// deleted — the same sequence.
fn in_smallest_form(
    seq: &[u8],
    mut span: Range<usize>,
    mut deleted: FlankDeletion,
    (left_flank, right_flank): (&[u8], &[u8]),
    shortest: usize,
) -> (Range<usize>, FlankDeletion) {
    while deleted.left > 0 && span.len() > shortest {
        let first_deleted = left_flank[left_flank.len() - deleted.left as usize];
        if seq[span.start] != first_deleted {
            break;
        }
        span.start += 1;
        deleted.left -= 1;
    }
    while deleted.right > 0 && span.len() > shortest {
        let last_deleted = right_flank[deleted.right as usize - 1];
        if seq[span.end - 1] != last_deleted {
            break;
        }
        span.end -= 1;
        deleted.right -= 1;
    }
    (span, deleted)
}

/// **How many aligned bases a read has before the tract and after it**, as `(before,
/// after)` — match operations only, so deleted, skipped, clipped and inserted bases do not
/// count. Coordinates as in [`flank_bases_deleted_at_the_edges`].
fn aligned_bases_beyond_the_tract(
    cigar: &[CigarOp],
    ref_start: u32,
    tract: &Range<u32>,
) -> (u32, u32) {
    let (tract_start, tract_end) = (tract.start, tract.end);
    let (mut before, mut after, mut at) = (0, 0, ref_start);
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                let end = at + n;
                before += end.min(tract_start).saturating_sub(at);
                after += end.saturating_sub(at.max(tract_end));
                at = end;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => at += n,
            _ => {}
        }
    }
    (before, after)
}

/// **How many read bases sit unaligned at each end of a read** — a soft clip, or an
/// insertion with no aligned base beyond it — as `(at the start, at the end)`. Hard clips are
/// not in the read and are stepped over. These are the bases a long allele's extra tract can
/// hide in, which is why a read with any of them facing a tract is never ruled out of it by
/// its footprint alone.
fn unaligned_bases_at_the_ends(cigar: &[CigarOp]) -> (u32, u32) {
    fn run(ops: &mut dyn Iterator<Item = &CigarOp>) -> u32 {
        let mut bases = 0;
        for op in ops {
            match op {
                CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
                CigarOp::SoftClip(n) | CigarOp::Insertion(n) => bases += n,
                _ => break,
            }
        }
        bases
    }
    (run(&mut cigar.iter()), run(&mut cigar.iter().rev()))
}

/// The tract's span inside `region`, for a read the aligner has nothing to decide about.
///
/// **The condition is that the read reproduces the reference across the whole window.** Then
/// the read carries the reference's repeat count — a different count would have to appear in
/// the CIGAR as an insertion or a deletion — and the tract sits at the reference's own
/// offsets, so the answer is arithmetic rather than a search.
///
/// Four things are required, and each rules out a way the CIGAR could be lying:
///
/// - **No insertion or deletion anywhere in the read.** One outside the window still shifts
///   every reference-to-read offset after it, and `region` is derived from those offsets.
/// - **No clip at either end.** A soft clip is the one place a long allele's extra units can
///   hide from the aligned span — the deleted caller's gate admitted a clipped read
///   unconditionally for exactly this reason — so a clipped read goes
///   to the aligner however clean the rest of its CIGAR looks.
/// - **The aligned span brackets the whole window**, or the read is a partial observation and
///   belongs on the aligner's path, which is what decides how far into the tract it reached.
/// - **The bases match the reference over the whole window**, not merely their count. The
///   delimiter scores emissions by base quality, so *equal length* alone would not let this
///   function predict its answer; *equal bases* does, because then the best path is the
///   diagonal at any quality — every alternative either mismatches a base or pays a slip's
///   open cost, and both are strictly worse than a run of exact matches.
///
/// The last is why a `SeqMismatch` in the CIGAR is not enough to disqualify a read on its own
/// and is not tested for: what matters is the bases, and they are compared directly.
///
/// **A fifth condition is about the locus rather than the read, and it was found by
/// measurement rather than by argument.** Where the reference's own repeat run continues past
/// the boundary region typing drew, the delimiter answers a *lower bound* and not a length —
/// on a tomato run it returned `FromLeft(15..42)` for a 12-base `A` tract whose right flank is
/// more `A`, because the run reaches the window edge and might go further. That is the right
/// answer and this function must not overrule it, so a locus whose run is extendable in the
/// reference is left to the aligner entirely. It is checked against the reference and so is
/// the same verdict for every read at that locus.
fn tract_needing_no_alignment(
    read: &AlignedRead,
    region: &Range<usize>,
    locus: &SsrLocus,
    reference: &[u8],
) -> Option<Range<u64>> {
    if reference_run_escapes_the_tract(locus, reference) {
        return None;
    }
    let straight = read.cigar.iter().all(|op| {
        matches!(
            op,
            CigarOp::Match(_) | CigarOp::SeqMatch(_) | CigarOp::SeqMismatch(_)
        )
    });
    if !straight {
        return None;
    }
    // With no indel and no clip, `extract_region` returns exactly the window's width when the
    // aligned span brackets it, and something shorter when it does not — so this one length
    // test covers the bracketing requirement.
    if region.len() != reference.len() {
        return None;
    }
    if &read.seq[region.clone()] != reference {
        return None;
    }
    let left = locus.left_flank_len() as u64;
    Some(left..left + locus.repeat.tract_len())
}

/// Does the reference's own repeat run carry on past the tract region typing drew?
///
/// A run of period `p` continues one base to the left when that base repeats the one `p`
/// further along, and one base to the right when it repeats the one `p` back. Either way the
/// delimiter sees a run that does not end inside the window, so it answers a lower bound
/// rather than a length, and [`tract_needing_no_alignment`] must stand aside.
///
/// A flank shorter than the period cannot be tested and is treated as escaping, which costs
/// an alignment and never a wrong answer.
///
/// **Each edge is asked about the run that reaches it** (`str_bundle_loci.md` §6.1): a
/// segment's one motif at both, a bundle's hull the run at that edge — and nothing where
/// fixed sequence reaches the edge, since no run is there to carry on.
fn reference_run_escapes_the_tract(locus: &SsrLocus, reference: &[u8]) -> bool {
    let left = locus.left_flank_len();
    let tract_end = left + locus.repeat.tract_len() as usize;
    let periods = locus.repeat.edge_run_periods;
    let escapes_left = periods.left.is_some_and(|period| {
        period == 0
            || left < period
            || reference.len() < left + period
            || reference[left - 1] == reference[left - 1 + period]
    });
    let escapes_right = periods.right.is_some_and(|period| {
        period == 0
            || reference.len() < tract_end + period
            || tract_end < period
            || reference[tract_end] == reference[tract_end - period]
    });
    escapes_left || escapes_right
}

/// Align the read's `region` slice against the reference frame.
fn delimit<A: RepeatDelimiter>(
    aligner: &A,
    read: &AlignedRead,
    region: &Range<usize>,
    reference: &[u8],
    context: RepeatContext<'_>,
    scratch: &mut A::Scratch,
) -> RepeatSpan {
    let bases = ReadBases::try_new(&read.seq[region.clone()], &read.qual[region.clone()])
        .expect("seq and qual are the same slice, hence equal length");
    aligner.align(bases, reference, context, scratch)
}

fn to_usize(range: &Range<u64>) -> Range<usize> {
    range.start as usize..range.end as usize
}

/// A complete tract, if it clears the base-quality gate; else `LowQuality`. `tract` is in
/// **absolute read coordinates** — the callers translate whichever frame they hold
/// (region-relative for the delimiter's spans, absolute already for the input-alignment
/// span), so this function cannot be handed a range in the wrong frame silently.
fn complete_or_low_quality(
    read: &AlignedRead,
    tract: Range<usize>,
    flank_deletion: FlankDeletion,
    qual_buffer: &mut Vec<u8>,
) -> Classified {
    let tract_qual = &read.qual[tract.clone()];
    if passes_quality_gate(tract_qual, MIN_REGION_Q1, qual_buffer) {
        Classified::Observed {
            bases: read.seq[tract].into(),
            read_witness: ReadWitness::Complete,
            flank_deletion,
            held_flank: None,
            q_sum: ln_p_err_sum(tract_qual),
        }
    } else {
        Classified::NoObservation(NoObservationReason::LowQuality)
    }
}

/// `Σ ln(P_err)` over base qualities — the per-read base-quality error mass on the tract, in
/// log-error space (freebayes' `q_sum` convention). Negative and monotone: higher quality →
/// less error mass. The BQ support moment (spec §3); there is no BAQ on this path and MAPQ is
/// carried separately, so this is the base-quality term alone.
fn ln_p_err_sum(quals: &[u8]) -> f64 {
    const LN10_OVER_10: f64 = std::f64::consts::LN_10 / 10.0;
    quals.iter().map(|&q| -(q as f64) * LN10_OVER_10).sum()
}

/// A partial (lower-bound) observation: the tract bases the read showed, as the run of locus
/// positions it witnessed, anchored at whichever border held — `border`, the half of
/// `RepeatSpan`'s answer that decides where the witnessed run sits inside the locus.
fn partial(
    read: &AlignedRead,
    region: &Range<usize>,
    tract: &Range<u64>,
    border: HeldFlank,
    locus: &SsrLocus,
) -> Classified {
    // **The span has to be ordered, and only a `debug_assert` upstream says so.**
    // `classify`'s own ordering check (`ssr_best_path_flat_gap.rs`) is debug-only on a
    // struct with public fields — kept public for the delimiter parity harness — and
    // this repo has recorded twice that debug-only guards compile out of the build it
    // actually runs. An inverted span underflows `tract.end - tract.start` below: a
    // panic in debug, and in release a wrap that clamps to `u16::MAX` and then dies on
    // the slice index, naming neither the cause nor the read. The `tract.is_empty()`
    // guard removed from `classify_read` used to swallow this shape too, since
    // `Range::is_empty` is `start >= end` (Milestone C review).
    debug_assert!(
        tract.start <= tract.end,
        "inverted tract span {tract:?} reached `partial`",
    );
    let region_seq = &read.seq[region.clone()];
    let region_qual = &read.qual[region.clone()];
    let tract = to_usize(tract);
    // `reach` is the observed tract length in **read** coordinates, which diverge from locus
    // positions under stutter — an expanded allele reaches further in read bases than the
    // reference tract has positions.
    //
    // # The convention that turns a read length into reference positions (owner, 2026-07-31)
    //
    // A witness is a set of **reference** positions, so `reach` has to be placed on the
    // reference before it can be stored, and inside a repeat that placement needs a rule:
    // if a read shows 12 repeat bases where the reference tract has 10 positions, *which*
    // two are the extra copies is not determined by the sequence.
    //
    // **The rule is: lay the read's repeat down from the border it anchored.** A read that
    // held the left flank starts its repeat at the tract's left border; one that held the
    // right flank ends its repeat at the right border. Bases beyond the tract's far border
    // are then extra copies — an insertion — rather than positions of the tract. It is the
    // same rule indel left-alignment uses, with the anchored side choosing the direction.
    //
    // Two consequences, and the clamp in `ReadWitness::from_left` / `from_right` is exactly
    // this rule implemented:
    //
    // 1. a reach **shorter** than the tract covers that many positions from the anchored
    //    border, which is what the constructors build directly;
    // 2. a reach **at or past** the tract length covers the tract end to end — every
    //    reference position — with the surplus falling outside as inserted copies. That is
    //    what the clamp produces, and it is correct rather than a saturation artefact.
    //
    // **Covering every position is still not a measurement.** The read anchored one border
    // and then ran out; the allele can continue past what it showed, so the evidence stays a
    // lower bound. That is why this path builds a *partial* witness in both cases and only a
    // read holding **both** flanks is `Complete` — and why the dumps spell case 2
    // `partial:both` rather than folding it into `partial:left` (spec §8; 2,530 of 6,216
    // partial observations on chr01 of tomato SRR7279503).
    let reach = (tract.end - tract.start).min(u16::MAX as usize) as u16;
    let locus_len = LocusLen::from_positions(locus.repeat.tract_len());
    // **The read is outside the locus, and this is where that is decided.** The
    // constructors answer `None` when the clamped run covers no position, which is
    // exactly the read that anchored a flank and crossed no tract position: the
    // delimiter's rejected side falls back to the read's own edge, so an alignment that
    // never entered the repeat comes back as `tract_start == tract_end`. The bound it
    // would supply is *at least zero* — no evidence — and it arrives in bulk, because
    // the fetch queries the tract plus its margin so every read clipping the flank is
    // delimited. Those bases are the SNP/indel path's to analyse (see
    // `NoObservationReason::OutsideTract` for the measurement).
    //
    // A second `tract.is_empty()` guard stood in `classify_read` until the
    // Milestone C review showed no input could distinguish the two — the `Option` had
    // made this branch answer first. One decision, one place.
    let Some(read_witness) = (match border {
        HeldFlank::Left => ReadWitness::from_left(reach, locus_len),
        HeldFlank::Right => ReadWitness::from_right(reach, locus_len),
    }) else {
        return Classified::NoObservation(NoObservationReason::OutsideTract);
    };
    // **The side the clamp loses is written down** (`str_bundle_loci.md` §7.6, the review's
    // finding M3; segments since the merge of main's partial-read phase fix, 2026-10-12). A
    // witness that covers the whole locus no longer says which flank the read anchored, and
    // both a segment and a bundle compare a left-held read with the start of the repeat and a
    // right-held one with its end (`read_likelihoods.md` §5.2), so two such reads of the same
    // bases from the two sides stay two observations.
    let held_flank =
        (read_witness.is_flush_left() && read_witness.is_flush_right(locus_len)).then_some(border);
    // **At a bundle, a read ending inside the hull is kept as a partial too**: it counts toward
    // the hull's depth, which the hidden-duplication filter's coverage windows read, and calling
    // scores it as a lower bound from the flank it held (`str_bundle_loci.md` §6.1, §7.6).
    Classified::Observed {
        flank_deletion: crate::locus_generation::FlankDeletion::NONE,
        bases: region_seq[tract.clone()].into(),
        read_witness,
        held_flank,
        // Partials are kept without the quality gate (spec §3), but still carry the BQ moment.
        q_sum: ln_p_err_sum(&region_qual[tract]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::PerQualityEmission;
    use crate::alignment::ssr_best_path_flat_gap::{SsrFlatGapAligner, ViterbiScratch};
    use crate::bam::alignment_input::CigarOp;
    use crate::locus_generation::ssr::{EdgeRunPeriods, MeasuredRepeat};
    use crate::region_typing::segment_criteria::{Motif, SsrSegment};
    use crate::types::{Position, ReadGroupId};

    // Reference frame: 6-base flanks around a CACACA tract → "GGGGGGCACACATTTTTT".
    const FRAME: &[u8] = b"GGGGGGCACACATTTTTT";

    /// The flank deletion classification records for a read with this CIGAR.
    fn flank_bases_deleted(
        cigar: &[CigarOp],
        ref_start: u32,
        tract_start: u32,
        tract_end: u32,
    ) -> Option<FlankDeletion> {
        flank_deletion_of(flank_bases_deleted_at_the_edges(
            cigar,
            ref_start,
            &(tract_start..tract_end),
        ))
    }

    /// Whether classification reads this CIGAR as a deletion taking flank bases at an edge.
    fn deletion_takes_flank_bases_at_a_tract_edge(
        cigar: &[CigarOp],
        ref_start: u32,
        tract_start: u32,
        tract_end: u32,
    ) -> bool {
        flank_bases_deleted_at_the_edges(cigar, ref_start, &(tract_start..tract_end)) != (0, 0)
    }

    fn locus() -> SsrLocus {
        SsrLocus {
            repeat: SsrSegment::new("chr1".into(), 7, 12, Motif::new(b"CA").unwrap(), 1.0)
                .unwrap()
                .into(),
            tract_with_margin_bases: FRAME.into(),
            margin_start: Position(1),
        }
    }

    fn read(seq: &[u8], qual_value: u8) -> AlignedRead {
        AlignedRead {
            qname: b"r".to_vec(),
            flag: 0,
            ref_id: 0,
            pos: 1,
            mapq: 60,
            cigar: vec![CigarOp::Match(seq.len() as u32)],
            seq: seq.to_vec(),
            qual: vec![qual_value; seq.len()],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            read_group: ReadGroupId(0),
        }
    }

    /// Classify `read` against an STR bundle's hull `hull` between `left` and `right`
    /// flanks, as the bundle generator measures it: the flat-gap delimiter with the
    /// segments' firm-anchor test.
    fn classify_at_a_bundle(
        read: &AlignedRead,
        (left, hull, right): (&[u8], &[u8], &[u8]),
        edges: EdgeRunPeriods,
    ) -> Classified {
        use crate::alignment::ssr_unit_robust::AnchorRule;
        let frame = [left, hull, right].concat();
        let start = left.len() as u64 + 1;
        let bundle = SsrLocus {
            repeat: MeasuredRepeat::bundle_hull(
                start,
                start + hull.len() as u64 - 1,
                Motif::new(b"A").unwrap(),
                edges,
            ),
            tract_with_margin_bases: frame.into(),
            margin_start: Position(1),
        };
        classify_read(
            read,
            &bundle,
            &SsrFlatGapAligner::with_anchor_rule(PerQualityEmission::new(), AnchorRule::default()),
            &StutterModel::hipstr_shipped(),
            &mut ViterbiScratch::new(),
            &mut Vec::new(),
        )
    }

    /// **A deletion across a bundle's edge that the mapper wrote as substitutions is one
    /// allele with a flank deletion** (found in review, 2026-10-11). The stuttering
    /// benchmark's `bun00119`: the left flank's last three bases `CGA` are gone and the
    /// `(AAAC)6` run carries one copy more, so the read is the reference's length and the
    /// mapper writes it with no indel. The flat-gap delimiter deleted two flank bases and
    /// spelled `AC` into the hull, while the generic path called the flank change beside
    /// it — the event written twice. Now the read is the hull's spelling `AAC (AAAC)6 …` with
    /// three flank bases deleted, the allele the truth holds.
    #[test]
    fn a_bundle_read_whose_delimiter_deletes_flank_bases_records_them_as_a_flank_deletion() {
        let left = b"CGGACGCCGCTCCGA";
        let hull = b"AAACAAACAAACAAACAAACAAACGTCAATTAATGGTCTCTCTCTCTC";
        let right = b"AACGGTGCTATCCGC";
        let spelled = [&b"AAC"[..], hull].concat();
        let carrier = read(&[&left[..12], &spelled, right].concat(), 40);
        let edges = EdgeRunPeriods {
            left: Some(4),
            right: Some(2),
        };
        match classify_at_a_bundle(&carrier, (left, hull, right), edges) {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                flank_deletion,
                ..
            } => {
                assert_eq!(&*bases, spelled.as_slice());
                assert_eq!(flank_deletion, FlankDeletion { left: 3, right: 0 });
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
        // The reference read is measured as before: no flank deletion.
        let reference = read(&[&left[..], hull, right].concat(), 40);
        match classify_at_a_bundle(&reference, (left, hull, right), edges) {
            Classified::Observed {
                bases,
                flank_deletion,
                ..
            } => {
                assert_eq!(&*bases, hull.as_slice());
                assert!(flank_deletion.is_none());
            }
            other => panic!("expected an observation, got {other:?}"),
        }
    }

    /// **A deletion across a bundle's edge that the mapper wrote inside the flank, with an
    /// insertion at the hull's first base, is one allele with a flank deletion** (found in
    /// review, 2026-10-11). The benchmark's `bun00018`: `TCCGA` is gone from the left flank
    /// and the `(AAAT)6` run gained a copy, written `10M 4D 2M 6I 61M`. No deletion crosses
    /// the hull's edge, so the own-alignment rule does not fire, and the delimiter cannot
    /// anchor a flank whose last five bases are missing; the read was a partial, which
    /// calling leaves out at a bundle.
    #[test]
    fn a_bundle_read_with_a_deletion_inside_the_flank_records_it_as_a_flank_deletion() {
        let left = b"GGACATTGTGTCCGA";
        let hull = b"AAATAAATAAATAAATAAATAAATTCCGGCGAGCCTATATATATATA";
        let right = b"ACGGTAGTCCAGCGG";
        let spelled = [&b"AAT"[..], &b"AAAT".repeat(7), &hull[24..]].concat();
        let mut carrier = read(&[&left[..10], &spelled, right].concat(), 40);
        carrier.cigar = vec![
            CigarOp::Match(10),
            CigarOp::Deletion(4),
            CigarOp::Match(2),
            CigarOp::Insertion(6),
            CigarOp::Match(61),
        ];
        let edges = EdgeRunPeriods {
            left: Some(4),
            right: Some(2),
        };
        match classify_at_a_bundle(&carrier, (left, hull, right), edges) {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                flank_deletion,
                ..
            } => {
                assert_eq!(&*bases, spelled.as_slice());
                assert_eq!(flank_deletion, FlankDeletion { left: 5, right: 0 });
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    fn classify(read: &AlignedRead) -> Classified {
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let mut scratch = ViterbiScratch::new();
        let mut qual_buffer = Vec::new();
        classify_read(
            read,
            &locus(),
            &aligner,
            &stutter,
            &mut scratch,
            &mut qual_buffer,
        )
    }

    /// A read spanning the whole frame measures the tract exactly — a complete observation
    /// of the reference allele.
    #[test]
    fn a_spanning_read_yields_a_complete_observation() {
        match classify(&read(FRAME, 40)) {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                flank_deletion: FlankDeletion::NONE,
                held_flank: None,
                q_sum,
            } => {
                assert_eq!(&*bases, b"CACACA");
                // Six Q40 bases: Σ ln(P_err) = 6 · 40 · (−ln10/10), all negative.
                assert!(q_sum < 0.0, "the BQ error mass is negative, got {q_sum}");
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    // A homopolymer frame for the input-alignment spelling tests: 6-base flanks around
    // a 6-base poly-A. The left flank ends in `T` and the right flank starts with `C`,
    // so neither continues the run and every junction rule below tests what it claims.
    const FRAME_A: &[u8] = b"GGGTTTAAAAAACCCGGG";

    fn locus_a() -> SsrLocus {
        SsrLocus {
            repeat: SsrSegment::new("chr1".into(), 7, 12, Motif::new(b"A").unwrap(), 1.0)
                .unwrap()
                .into(),
            tract_with_margin_bases: FRAME_A.into(),
            margin_start: Position(1),
        }
    }

    fn classify_at_a(read: &AlignedRead) -> Classified {
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let mut scratch = ViterbiScratch::new();
        let mut qual_buffer = Vec::new();
        classify_read(
            read,
            &locus_a(),
            &aligner,
            &stutter,
            &mut scratch,
            &mut qual_buffer,
        )
    }

    fn observed_bases(classified: Classified) -> Vec<u8> {
        match classified {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                ..
            } => bases.into_vec(),
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    /// **The `chr3:33,877,690` shape, at unit-test size: a real flank deletion beside
    /// the junction plus a substituted first tract base.** The mapper spells it as a
    /// one-base deletion in the flank and a `C` against the tract's first `A`; the
    /// honest tract account is `CAAAAA` — same length as the reference's, substituted
    /// edge. The delimiter's stutter-priced re-alignment absorbs the `C` into the flank
    /// and reports a bare 5-`A` run instead (the L4 defect, sized at 232 of 463
    /// never-offered tracts on HG002 at 30×); the input-alignment spelling keeps it.
    #[test]
    fn a_flank_deletion_beside_the_junction_keeps_the_tract_edge() {
        let mut read = read(b"GGGTTCAAAAACCCGGG", 40);
        read.cigar = vec![CigarOp::Match(5), CigarOp::Deletion(1), CigarOp::Match(12)];
        assert_eq!(observed_bases(classify_at_a(&read)), b"CAAAAA");
    }

    /// A right-anchored whole-unit expansion (`M I M`, the simulator's and some mappers'
    /// spelling) measures at the expanded length: the inserted `A`s left-align into the
    /// tract, so they are the tract's.
    #[test]
    fn a_right_anchored_expansion_is_the_tracts() {
        let mut read = read(b"GGGTTTAAAAAAAACCCGGG", 40);
        read.cigar = vec![CigarOp::Match(12), CigarOp::Insertion(2), CigarOp::Match(6)];
        assert_eq!(observed_bases(classify_at_a(&read)), b"AAAAAAAA");
    }

    /// A right-junction insertion that cannot left-align into the tract (its last base
    /// is not the tract's) is flank content and stays out of the observation. The tract
    /// carries a substitution so the input account is non-reference and the junction
    /// rule — not the reference-deferral to the delimiter — is what decides.
    #[test]
    fn a_right_junction_foreign_insertion_is_the_flanks() {
        let mut read = read(b"GGGTTTAACAAAGGCCCGGG", 40);
        read.cigar = vec![CigarOp::Match(12), CigarOp::Insertion(2), CigarOp::Match(6)];
        assert_eq!(observed_bases(classify_at_a(&read)), b"AACAAA");
    }

    /// A left-junction insertion that cannot left-align into the flank belongs to the
    /// tract — the boundary convention: an insertion at a junction goes with what
    /// follows, unless left-alignment carries it across.
    #[test]
    fn a_left_junction_foreign_insertion_is_the_tracts() {
        let mut read = read(b"GGGTTTAGAAAAAACCCGGG", 40);
        read.cigar = vec![CigarOp::Match(6), CigarOp::Insertion(2), CigarOp::Match(12)];
        assert_eq!(observed_bases(classify_at_a(&read)), b"AGAAAAAA");
    }

    /// A left-junction insertion whose bases left-align into the flank (`TT` beside a
    /// flank ending in `T`) is flank content, and the tract reads without it. The tract
    /// carries a substitution so the input account is non-reference and the junction
    /// rule — not the reference-deferral to the delimiter — is what decides.
    #[test]
    fn a_left_junction_flank_run_insertion_is_the_flanks() {
        let mut read = read(b"GGGTTTTTAACAAACCCGGG", 40);
        read.cigar = vec![CigarOp::Match(6), CigarOp::Insertion(2), CigarOp::Match(12)];
        assert_eq!(observed_bases(classify_at_a(&read)), b"AACAAA");
    }

    /// A read whose tract quality is below the gate is `LowQuality`, not an observation.
    #[test]
    fn a_low_quality_tract_is_rejected() {
        // Whole read at quality 5 (< MIN_REGION_Q1 = 15).
        let classified = classify(&read(FRAME, 5));
        assert_eq!(
            classified,
            Classified::NoObservation(NoObservationReason::LowQuality)
        );
    }

    /// A read that anchors the left flank but ends inside the tract is a partial (left)
    /// observation — a lower bound production would have discarded.
    #[test]
    fn a_read_running_off_the_right_is_a_left_partial() {
        // left flank + 4 tract bases, then the read ends (no right flank).
        let classified = classify(&read(b"GGGGGGCACA", 40));
        match classified {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Partial { positions },
                ..
            } => {
                assert_eq!(&*bases, b"CACA");
                // Flush with the left border: the run starts at locus position 0 and
                // covers four. Half-open, so it ends at 4.
                assert_eq!(positions.runs().collect::<Vec<_>>(), vec![(0, 4)]);
            }
            other => panic!("expected a left partial, got {other:?}"),
        }
    }

    /// **A segment records the held flank of a partial over the whole tract**, from either
    /// side: the witness covers every position and no longer says it, and a segment compares
    /// a left-held read with the tract's start and a right-held one with its end
    /// (`read_likelihoods.md` §5.2; `str_bundle_loci.md` §7.6).
    #[test]
    fn a_segment_partial_over_the_whole_tract_records_its_flank() {
        let mut right_held = read(b"CACACACACATTTTTT", 40);
        right_held.pos = 3;
        for (read, flank) in [
            (read(b"GGGGGGCACACACACA", 40), HeldFlank::Left),
            (right_held, HeldFlank::Right),
        ] {
            match classify(&read) {
                Classified::Observed {
                    read_witness: ReadWitness::Partial { positions },
                    held_flank,
                    ..
                } => {
                    assert_eq!(positions.runs().collect::<Vec<_>>(), vec![(0, 6)]);
                    assert_eq!(held_flank, Some(flank));
                }
                other => panic!("expected a partial over the whole tract, got {other:?}"),
            }
        }
    }

    /// A read that anchors the right flank but begins inside the tract is a partial
    /// (right) — the mirror of the left case, which a left/right swap would fail.
    ///
    /// **The offset is the assertion that matters here.** Since the reshape the side is not
    /// a variant but a derivation from where the run sits, so `offset_in_locus == 2` on a
    /// 6-base tract is what says "flush with the right border" — and it is the one number a
    /// mint that forgot the locus length could not produce.
    #[test]
    fn a_read_running_off_the_left_is_a_right_partial() {
        // 4 tract bases + full right flank, mapped at the tract's 3rd base (0-based 8).
        let mut read = read(b"CACATTTTTT", 40);
        read.pos = 9;
        match classify(&read) {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Partial { positions },
                ..
            } => {
                assert_eq!(&*bases, b"CACA");
                // The tract is "CACACA" — 6 positions — so a 4-position run flush with the
                // right border starts at 2 and ends at the locus's own end.
                assert_eq!(positions.runs().collect::<Vec<_>>(), vec![(6 - 4, 6)]);
            }
            other => panic!("expected a right partial, got {other:?}"),
        }
    }

    /// **A read that stops at the tract's edge is outside the locus, not a partial of
    /// length zero.**
    ///
    /// It covers the left flank and none of `CACACA`. The delimiter anchors the left
    /// flank and, having no right anchor, falls the far end back to the read's own edge
    /// — which is the same offset, so the span is empty. Before this was classified it
    /// became an observation with **empty bases** and a witness covering zero positions,
    /// and it was not a corner: 6,704 reads against 7,085 genuine partials on tomato
    /// chr01 of `SRR7279503`, at a median window overlap of 16 bases against a 30-base
    /// flank. Those bases belong to the SNP/indel path (owner, 2026-07-30).
    ///
    /// Both borders are asserted because the two arrive by mirror-image routes — the
    /// left case through `tract_start == region_len`, the right through
    /// `tract_end == 0`. **They now share one decision point** (`partial`'s handling of
    /// the constructors' `None`), so neither can be fixed without the other; asserting
    /// both is what keeps that true if the routes ever separate again. The mutation that
    /// discriminates is on the constructor: make `from_left` saturate a zero-length run
    /// to one position and both cases here fail, along with
    /// `a_constructor_asked_for_no_positions_answers_none`.
    #[test]
    fn a_read_covering_only_a_flank_is_outside_the_tract() {
        // The left flank alone: the read ends where the tract begins.
        assert_eq!(
            classify(&read(b"GGGGGG", 40)),
            Classified::NoObservation(NoObservationReason::OutsideTract),
            "a read that stops at the tract's first base witnessed no tract position",
        );

        // The right flank alone, mapped at the base just past the tract.
        let mut right_only = read(b"TTTTTT", 40);
        right_only.pos = 13;
        assert_eq!(
            classify(&right_only),
            Classified::NoObservation(NoObservationReason::OutsideTract),
            "and the mirror case, which arrives through the other end of the span",
        );
    }

    /// **A read that stops where the tract begins is outside it, even where the flanks echo
    /// each other** (owner, 2026-10-10). The stuttering benchmark's segment: the left flank
    /// ends `CCACCA` and the right flank begins `CCACCA`, so the production delimiter
    /// matched a flank-only read's last six bases to the right flank and measured an empty
    /// tract — a complete observation that stopped the locus. Without the check this
    /// 30-base read fares no better: it comes back as a *partial* observation of the tract
    /// spelled from 24 flank bases. The read's own CIGAR says it never reached the tract,
    /// and that is now decided before any alignment.
    #[test]
    fn a_read_stopping_at_the_tract_is_outside_it_when_the_flanks_echo() {
        use crate::alignment::ssr_unit_robust::{SsrUnitRobustAligner, UnitRobustScratch};
        let left_flank = b"TAAGCTGTACTCTTCCAGTCACCGCCACCA";
        let tract = b"AACAAT".repeat(11);
        let right_flank = b"CCACCAACGTGACCGTGACCGTACTAGCTC";
        let frame = [left_flank.as_slice(), &tract, right_flank].concat();
        let echo_locus = SsrLocus {
            repeat: SsrSegment::new(
                "chr1".into(),
                31,
                30 + tract.len() as u64,
                Motif::new(b"AACAAT").unwrap(),
                1.0,
            )
            .unwrap()
            .into(),
            tract_with_margin_bases: frame.into(),
            margin_start: Position(1),
        };
        let flank_only = read(left_flank, 40);
        let classified = classify_read(
            &flank_only,
            &echo_locus,
            &SsrUnitRobustAligner::new(PerQualityEmission::new()),
            &StutterModel::hipstr_shipped(),
            &mut UnitRobustScratch::default(),
            &mut Vec::new(),
        );
        assert_eq!(
            classified,
            Classified::NoObservation(NoObservationReason::OutsideTract)
        );
    }

    /// **A read carrying a deletion across the tract's edge crossed the tract** (owner,
    /// 2026-10-10). The stuttering benchmark's shape: a 13-base deletion taking the last nine
    /// bases of an `AAGAG` repeat and the first four of the right flank. The delimiter could
    /// not match the damaged flank and returned a partial; the read's own alignment says it
    /// crossed, and spells the tract as the 21 bases the deletion left.
    #[test]
    fn a_read_with_a_deletion_across_the_tract_edge_is_complete() {
        use crate::alignment::ssr_unit_robust::{SsrUnitRobustAligner, UnitRobustScratch};
        let left_flank = b"GTCAGCATTGCAGGTCATCGGACTTAGCCA";
        let tract = b"AAGAG".repeat(6);
        let right_flank = b"CCACGATTGCAGTCAGCTTGACCATAGCAG";
        let frame = [left_flank.as_slice(), &tract, right_flank].concat();
        let locus = SsrLocus {
            repeat: SsrSegment::new("chr1".into(), 31, 60, Motif::new(b"AAGAG").unwrap(), 1.0)
                .unwrap()
                .into(),
            tract_with_margin_bases: frame.into(),
            margin_start: Position(1),
        };
        let seq = [left_flank.as_slice(), &tract[..21], &right_flank[4..]].concat();
        let mut carrier = read(&seq, 40);
        carrier.cigar = vec![
            CigarOp::Match(51),
            CigarOp::Deletion(13),
            CigarOp::Match(26),
        ];
        let classified = classify_read(
            &carrier,
            &locus,
            &SsrUnitRobustAligner::new(PerQualityEmission::new()),
            &StutterModel::hipstr_shipped(),
            &mut UnitRobustScratch::default(),
            &mut Vec::new(),
        );
        match classified {
            Classified::Observed {
                bases,
                read_witness,
                ..
            } => {
                assert_eq!(&*bases, &tract[..21]);
                assert_eq!(read_witness, ReadWitness::Complete);
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
        // The four flank bases the deletion took are recorded beside the spelling (owner,
        // 2026-10-10): the deletion cannot slide left, since the tract's `A` before it is not
        // the flank's `C` at its end.
        match classify_read(
            &carrier,
            &locus,
            &SsrUnitRobustAligner::new(PerQualityEmission::new()),
            &StutterModel::hipstr_shipped(),
            &mut UnitRobustScratch::default(),
            &mut Vec::new(),
        ) {
            Classified::Observed { flank_deletion, .. } => {
                assert_eq!(flank_deletion, FlankDeletion { left: 0, right: 4 })
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    /// **Reads of one sequence get one key, whatever split the mapper chose** (owner,
    /// 2026-10-10; `str_edge_deletion_alleles.md` §4). Two `A`s deleted from `…GA|AAGAAG…`
    /// fit wholly inside the repeat: a mapper that put the gap across the edge, taking the
    /// flank's last `A`, still yields the tract without its first two bases and no flank
    /// deletion — the same observation as a read with the gap inside. Where the flank's last
    /// base is a `C`, the same CIGAR takes a base the repeat cannot give back, and the read
    /// records one flank base missing on the left.
    #[test]
    fn a_deletion_across_the_edge_is_recorded_in_its_smallest_form() {
        use crate::alignment::ssr_unit_robust::{SsrUnitRobustAligner, UnitRobustScratch};
        let tract = b"AAG".repeat(10);
        let right_flank = b"CCACGATTGCAGTCAGCTTGACCATAGCAG";
        for (left_flank, expected_bases, expected_deleted) in [
            (
                b"GTCAGCATTGCAGGTCATCGGACTTAGCGA",
                &tract[2..],
                FlankDeletion::NONE,
            ),
            (
                b"GTCAGCATTGCAGGTCATCGGACTTAGCGC",
                &tract[1..],
                FlankDeletion { left: 1, right: 0 },
            ),
        ] {
            let frame = [left_flank.as_slice(), &tract, right_flank].concat();
            let locus = SsrLocus {
                repeat: SsrSegment::new("chr1".into(), 31, 60, Motif::new(b"AAG").unwrap(), 1.0)
                    .unwrap()
                    .into(),
                tract_with_margin_bases: frame.into(),
                margin_start: Position(1),
            };
            let seq = [&left_flank[..29], &tract[1..], right_flank].concat();
            let mut carrier = read(&seq, 40);
            carrier.cigar = vec![CigarOp::Match(29), CigarOp::Deletion(2), CigarOp::Match(59)];
            match classify_read(
                &carrier,
                &locus,
                &SsrUnitRobustAligner::new(PerQualityEmission::new()),
                &StutterModel::hipstr_shipped(),
                &mut UnitRobustScratch::default(),
                &mut Vec::new(),
            ) {
                Classified::Observed {
                    bases,
                    read_witness,
                    flank_deletion,
                    ..
                } => {
                    assert_eq!(&*bases, expected_bases);
                    assert_eq!(read_witness, ReadWitness::Complete);
                    assert_eq!(flank_deletion, expected_deleted);
                }
                other => panic!("expected a complete observation, got {other:?}"),
            }
        }
    }

    /// A 15-base-flank locus around `tract`, its reference flanks those of
    /// SL4.0ch05:51017179, the read of `seq` aligned at the window's first base by `cigar`,
    /// classified by the shipped delimiter.
    fn classify_at_ch05(tract: &[u8], seq: &[u8], cigar: Vec<CigarOp>) -> Classified {
        use crate::alignment::ssr_unit_robust::{SsrUnitRobustAligner, UnitRobustScratch};
        let frame = [&b"AAATGAGGATTAAGG"[..], tract, b"TATGAATAGAAAATA"].concat();
        let locus = SsrLocus {
            repeat: SsrSegment::new(
                "chr1".into(),
                16,
                15 + tract.len() as u64,
                Motif::new(&tract[..1]).unwrap(),
                1.0,
            )
            .unwrap()
            .into(),
            tract_with_margin_bases: frame.into(),
            margin_start: Position(1),
        };
        let mut carrier = read(seq, 40);
        carrier.cigar = cigar;
        classify_read(
            &carrier,
            &locus,
            &SsrUnitRobustAligner::new(PerQualityEmission::new()),
            &StutterModel::hipstr_shipped(),
            &mut UnitRobustScratch::default(),
            &mut Vec::new(),
        )
    }

    /// **A read whose own alignment opens an indel in the flank bases the delimiter's guard
    /// closes is spelled by that alignment** (2026-11-01). SL4.0ch05:51017179: the reads
    /// carry `ATAAAGAA` over an `A` ×8 tract and delete `ATAGA` five bases into the right
    /// flank. The guard forbids that gap, the delimiter's best path left ends the tract early
    /// with the flank laid on mismatches, the anchor rule drops the right side, and the
    /// read came back a partial that tied every 8-base candidate.
    #[test]
    fn a_flank_indel_the_guard_forbids_is_spelled_by_the_reads_alignment() {
        use CigarOp::{Deletion, Match};
        let seq = [&b"AAATGAGGATTAAGG"[..], b"ATAAAGAA", b"TATGA", b"AAATA"].concat();
        match classify_at_ch05(b"AAAAAAAA", &seq, vec![Match(28), Deletion(5), Match(5)]) {
            Classified::Observed {
                bases,
                read_witness,
                flank_deletion,
                ..
            } => {
                assert_eq!(&*bases, b"ATAAAGAA");
                assert_eq!(read_witness, ReadWitness::Complete);
                // The deletion does not touch the tract: the generic path calls it.
                assert_eq!(flank_deletion, FlankDeletion::NONE);
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    /// **A deletion at the tract's edge runs on through the indels the read's alignment has in
    /// the guarded flank bases** (2026-11-01). A sample carrying one more `AAGC` and the
    /// flank's `GAACT` deleted one base past the tract: minimap2 wrote its reads as the
    /// reference's eight repeats, the flank's `CG` deleted, two matches, the `G` of the
    /// extra repeat inserted, and a mismatch. Taken alone, the edge deletion made the allele
    /// the reference tract with two flank bases gone, which no read carried.
    #[test]
    fn an_edge_deletion_runs_on_through_the_flank_indels_beside_it() {
        use crate::alignment::ssr_unit_robust::{SsrUnitRobustAligner, UnitRobustScratch};
        use CigarOp::{Deletion, Insertion, Match};
        let left_flank = b"GTCAGCATTGCAGGTCATCGGACTTAGCCA";
        let tract = b"AAGC".repeat(8);
        let right_flank = b"CGAACTAGTGGAGCAGTTGACCATAGCAGT";
        let frame = [left_flank.as_slice(), &tract, right_flank].concat();
        let locus = SsrLocus {
            repeat: SsrSegment::new("chr1".into(), 31, 62, Motif::new(b"AAGC").unwrap(), 1.0)
                .unwrap()
                .into(),
            tract_with_margin_bases: frame.into(),
            margin_start: Position(1),
        };
        // The truth: nine repeats, the flank's `C`, then the flank after `GAACT`.
        let seq = [
            left_flank.as_slice(),
            &b"AAGC".repeat(9),
            b"C",
            &right_flank[6..],
        ]
        .concat();
        let mut carrier = read(&seq, 40);
        carrier.cigar = vec![
            Match(62),
            Deletion(2),
            Match(2),
            Insertion(1),
            Match(seq.len() as u32 - 65),
        ];
        match classify_read(
            &carrier,
            &locus,
            &SsrUnitRobustAligner::new(PerQualityEmission::new()),
            &StutterModel::hipstr_shipped(),
            &mut UnitRobustScratch::default(),
            &mut Vec::new(),
        ) {
            Classified::Observed {
                bases,
                read_witness,
                flank_deletion,
                ..
            } => {
                // The read's bases up to where its alignment rejoins the reference after the
                // insertion, in place of the tract and the four flank bases `CGAA`; the
                // mismatch after it is the generic path's.
                assert_eq!(&*bases, [tract.as_slice(), b"AAG"].concat().as_slice());
                assert_eq!(read_witness, ReadWitness::Complete);
                assert_eq!(flank_deletion, FlankDeletion { left: 0, right: 4 });
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    /// **A deletion taking a run of flank bases equal to the tract's keeps one copy in the
    /// spelling** (2026-11-01; HG002 chr20:2,472,054). `…TC AAAAAAG | A ×8`, and the reads
    /// delete `AAAAAAGAA`: seven flank bases and two of the tract. Smallest form gave every
    /// spelled `A` back to the flank, the spelling came out empty, and candidate selection
    /// dropped an allele without a copy of the motif.
    #[test]
    fn smallest_form_keeps_a_whole_copy_in_a_segments_spelling() {
        use crate::alignment::ssr_unit_robust::{SsrUnitRobustAligner, UnitRobustScratch};
        use CigarOp::{Deletion, Match};
        let left_flank = [&b"GTCAGCATTGCAGGTCATCGGTC"[..], b"AAAAAAG"].concat();
        let tract = b"A".repeat(8);
        let right_flank = b"CCACGATTGCAGTCAGCTTGACCATAGCAG";
        let frame = [left_flank.as_slice(), &tract, right_flank].concat();
        let locus = SsrLocus {
            repeat: SsrSegment::new("chr1".into(), 31, 38, Motif::new(b"A").unwrap(), 1.0)
                .unwrap()
                .into(),
            tract_with_margin_bases: frame.into(),
            margin_start: Position(1),
        };
        let seq = [&left_flank[..23], &tract[2..], right_flank].concat();
        let mut carrier = read(&seq, 40);
        carrier.cigar = vec![Match(23), Deletion(9), Match(36)];
        match classify_read(
            &carrier,
            &locus,
            &SsrUnitRobustAligner::new(PerQualityEmission::new()),
            &StutterModel::hipstr_shipped(),
            &mut UnitRobustScratch::default(),
            &mut Vec::new(),
        ) {
            Classified::Observed {
                bases,
                read_witness,
                flank_deletion,
                ..
            } => {
                assert_eq!(&*bases, b"A");
                assert_eq!(read_witness, ReadWitness::Complete);
                assert_eq!(flank_deletion, FlankDeletion { left: 2, right: 0 });
            }
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    /// Which indels sit in the guarded flank bases, and how far an edge deletion grows
    /// through them: the tract at 30–60, eight guarded bases each side.
    #[test]
    fn indels_in_the_guarded_flank_are_found_and_an_edge_deletion_grows_through_them() {
        use CigarOp::{Deletion, Insertion, Match};
        let widths = (8, 8);
        for (cigar, expected) in [
            (vec![Match(65), Deletion(5), Match(20)], true),
            (vec![Match(68), Deletion(5), Match(20)], false),
            (vec![Match(18), Deletion(3), Match(70)], false),
            (vec![Match(24), Deletion(3), Match(70)], true),
            (vec![Match(62), Insertion(2), Match(20)], true),
            (vec![Match(68), Insertion(2), Match(20)], false),
            (vec![Match(40), Deletion(3), Match(40)], false),
        ] {
            assert_eq!(
                indel_in_the_guarded_flank(&cigar, 0, &(30..60), widths),
                expected,
                "{cigar:?}"
            );
        }
        let fp = Footprint {
            ref_start: 0,
            ref_end: 100,
            leading_clip: 0,
            trailing_clip: 0,
        };
        // An edge deletion of two flank bases, two matches, an insertion: the flank deletion
        // grows to the insertion, four bases, and the spelling to the read base after it.
        let cigar = vec![Match(60), Deletion(2), Match(2), Insertion(1), Match(36)];
        assert_eq!(
            through_the_flank_indels(
                &cigar,
                fp,
                30..60,
                FlankDeletion { left: 0, right: 2 },
                &(30..60),
                widths,
            ),
            (30..63, FlankDeletion { left: 0, right: 4 })
        );
        // Nothing else in the guard: unchanged.
        let cigar = vec![Match(60), Deletion(2), Match(38)];
        assert_eq!(
            through_the_flank_indels(
                &cigar,
                fp,
                30..60,
                FlankDeletion { left: 0, right: 2 },
                &(30..60),
                widths,
            ),
            (30..60, FlankDeletion { left: 0, right: 2 })
        );
    }

    /// The smallest form on both sides: a spelling gives a flank base back while it starts
    /// (or ends) with the deleted flank base nearest it, and stops at the first that differs.
    /// And the flank bases a deletion takes are counted on each side of the tract.
    #[test]
    fn the_smallest_form_gives_back_flank_bases_and_deletions_are_counted() {
        use CigarOp::{Deletion, Match};
        // Spelling "AAT" with two left flank bases "AA" missing: both given back.
        assert_eq!(
            in_smallest_form(
                b"AAT",
                0..3,
                FlankDeletion { left: 2, right: 0 },
                (b"GAA", b""),
                0,
            ),
            (2..3, FlankDeletion::NONE)
        );
        // "CAT" with flank "GA" missing on the left: the first deleted base is `G`, kept.
        assert_eq!(
            in_smallest_form(
                b"CAT",
                0..3,
                FlankDeletion { left: 2, right: 0 },
                (b"GA", b""),
                0,
            ),
            (0..3, FlankDeletion { left: 2, right: 0 })
        );
        // Right side: "AGT" missing "TC" after the tract: the last deleted base is `C`, kept;
        // "AGC" missing "TC": `C` given back, then `T` against `G`, kept.
        assert_eq!(
            in_smallest_form(
                b"AGT",
                0..3,
                FlankDeletion { left: 0, right: 2 },
                (b"", b"TC"),
                0,
            ),
            (0..3, FlankDeletion { left: 0, right: 2 })
        );
        assert_eq!(
            in_smallest_form(
                b"AGC",
                0..3,
                FlankDeletion { left: 0, right: 2 },
                (b"", b"TC"),
                0,
            ),
            (0..2, FlankDeletion { left: 0, right: 1 })
        );
        // chr20:2,472,054: six `A`s spelled with `AAAAAAG` missing on the left. With no floor
        // every `A` is given back and the spelling is empty; with an `A` tract's floor of one
        // base, one `A` stays and `AG` is the deletion — the same sequence.
        assert_eq!(
            in_smallest_form(
                b"AAAAAA",
                0..6,
                FlankDeletion { left: 7, right: 0 },
                (b"AAAAAAG", b""),
                0,
            ),
            (6..6, FlankDeletion { left: 1, right: 0 })
        );
        assert_eq!(
            in_smallest_form(
                b"AAAAAA",
                0..6,
                FlankDeletion { left: 7, right: 0 },
                (b"AAAAAAG", b""),
                1,
            ),
            (5..6, FlankDeletion { left: 2, right: 0 })
        );
        // Tract at [30, 60).
        assert_eq!(
            flank_bases_deleted(&[Match(27), Deletion(6), Match(40)], 0, 30, 60),
            Some(FlankDeletion { left: 3, right: 0 })
        );
        assert_eq!(
            flank_bases_deleted(&[Match(51), Deletion(13), Match(26)], 0, 30, 60),
            Some(FlankDeletion { left: 0, right: 4 })
        );
        assert_eq!(
            flank_bases_deleted(&[Match(28), Deletion(34), Match(26)], 0, 30, 60),
            Some(FlankDeletion { left: 2, right: 2 }),
            "a deletion over the whole tract takes flank bases on both sides"
        );
        assert_eq!(
            flank_bases_deleted(&[Match(40), Deletion(5), Match(40)], 0, 30, 60),
            None
        );
        assert_eq!(
            flank_bases_deleted(&[Match(29), Deletion(1), Match(40)], 0, 30, 60),
            Some(FlankDeletion { left: 1, right: 0 }),
            "a deletion ending at the tract's first base takes the flank base before it"
        );
        assert_eq!(
            flank_bases_deleted(&[Match(60), Deletion(2), Match(20)], 0, 30, 60),
            Some(FlankDeletion { left: 0, right: 2 }),
            "a deletion starting just after the tract's last base takes the flank bases after it"
        );
    }

    /// The two CIGAR measurements behind that rule: a deletion crosses an edge only when it
    /// takes bases from both sides of it, and aligned bases are counted outside the tract
    /// and outside deletions.
    #[test]
    fn a_deletion_crossing_an_edge_and_the_aligned_bases_beyond_are_read_off_the_cigar() {
        use CigarOp::{Deletion, Match, SoftClip};
        // Tract at [30, 60).
        let across_the_end = [Match(51), Deletion(13), Match(26)];
        assert!(deletion_takes_flank_bases_at_a_tract_edge(
            &across_the_end,
            0,
            30,
            60
        ));
        assert_eq!(
            aligned_bases_beyond_the_tract(&across_the_end, 0, &(30..60)),
            (30, 26)
        );
        let across_the_start = [Match(27), Deletion(6), Match(40)];
        assert!(deletion_takes_flank_bases_at_a_tract_edge(
            &across_the_start,
            0,
            30,
            60
        ));
        let inside = [Match(40), Deletion(5), Match(40)];
        assert!(!deletion_takes_flank_bases_at_a_tract_edge(
            &inside, 0, 30, 60
        ));
        let starting_at_the_first_base = [Match(30), Deletion(5), Match(40)];
        assert!(
            !deletion_takes_flank_bases_at_a_tract_edge(&starting_at_the_first_base, 0, 30, 60),
            "a deletion that starts at the tract's first base takes nothing from the flank"
        );
        let ending_at_the_last_base = [Match(55), Deletion(5), Match(40)];
        assert!(
            !deletion_takes_flank_bases_at_a_tract_edge(&ending_at_the_last_base, 0, 30, 60),
            "a deletion that ends at the tract's last base takes nothing from the flank"
        );
        let touching_the_start = [Match(29), Deletion(1), Match(40)];
        assert!(
            deletion_takes_flank_bases_at_a_tract_edge(&touching_the_start, 0, 30, 60),
            "a flank deletion ending at the tract's first base touches the edge"
        );
        let touching_the_end = [Match(60), Deletion(1), Match(40)];
        assert!(
            deletion_takes_flank_bases_at_a_tract_edge(&touching_the_end, 0, 30, 60),
            "a flank deletion starting just after the tract's last base touches the edge"
        );
        assert_eq!(
            aligned_bases_beyond_the_tract(&[SoftClip(10), Match(80)], 10, &(30..60)),
            (20, 30)
        );
    }

    /// Unaligned bases at a read's ends: soft clips and edge insertions count, hard clips are
    /// stepped over, and an aligned base ends the run.
    #[test]
    fn unaligned_bases_at_the_ends_counts_clips_and_edge_insertions() {
        use CigarOp::{HardClip, Insertion, Match, SoftClip};
        assert_eq!(unaligned_bases_at_the_ends(&[Match(50)]), (0, 0));
        assert_eq!(
            unaligned_bases_at_the_ends(&[HardClip(3), SoftClip(4), Match(50), Insertion(6)]),
            (4, 6)
        );
        assert_eq!(
            unaligned_bases_at_the_ends(&[Match(10), Insertion(6), Match(10)]),
            (0, 0),
            "an insertion with aligned bases beyond it is not at an end"
        );
    }

    /// A read wholly inside the tract anchors neither flank — no per-read fact, so no
    /// observation.
    #[test]
    fn a_read_wholly_inside_the_tract_is_unanchored() {
        let mut read = read(b"CACACA", 40); // exactly the tract, no flanks
        read.pos = 7;
        assert_eq!(
            classify(&read),
            Classified::NoObservation(NoObservationReason::NoBorderAnchored)
        );
    }

    /// A long-allele read the mapper laid down all-Match pushes its far flank out of the
    /// ref-sized window; the first delimit sees a truncated flank, and widening the read
    /// slice recovers the full long allele (spec §2, the long-allele window recovery).
    #[test]
    fn a_window_truncated_long_allele_is_recovered_by_widening() {
        // 5 CA units (10 bp tract) instead of the reference's 3, mapped 22M.
        match classify(&read(b"GGGGGGCACACACACATTTTTT", 40)) {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                ..
            } => assert_eq!(&*bases, b"CACACACACA"),
            other => panic!("expected the recovered long allele, got {other:?}"),
        }
    }

    // --- An STR bundle's hull, measured as one repeat (`str_bundle_loci.md` §6.1) ---

    // 8-base flanks around a 19-base hull: (CA)4, a fixed TT, (GAT)3. Neither flank
    // continues the run at its edge.
    const BUNDLE_LEFT: &[u8] = b"GGATCCTG";
    const BUNDLE_HULL: &[u8] = b"CACACACATTGATGATGAT";
    const BUNDLE_RIGHT: &[u8] = b"CCTAGGAC";

    fn bundle_frame() -> Vec<u8> {
        [BUNDLE_LEFT, BUNDLE_HULL, BUNDLE_RIGHT].concat()
    }

    fn bundle_locus(edges: super::super::EdgeRunPeriods) -> SsrLocus {
        SsrLocus {
            repeat: super::super::MeasuredRepeat::bundle_hull(
                9,
                27,
                Motif::new(b"CA").unwrap(),
                edges,
            ),
            tract_with_margin_bases: bundle_frame().into(),
            margin_start: Position(1),
        }
    }

    fn runs_at_both_edges() -> super::super::EdgeRunPeriods {
        super::super::EdgeRunPeriods {
            left: Some(2),
            right: Some(3),
        }
    }

    fn classify_in_the_bundle(read: &AlignedRead) -> Classified {
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let mut scratch = ViterbiScratch::new();
        let mut qual_buffer = Vec::new();
        classify_read(
            read,
            &bundle_locus(runs_at_both_edges()),
            &aligner,
            &stutter,
            &mut scratch,
            &mut qual_buffer,
        )
    }

    /// A read crossing the hull is a complete observation spelled over the whole hull.
    #[test]
    fn a_read_crossing_the_hull_is_spelled_over_all_of_it() {
        match classify_in_the_bundle(&read(&bundle_frame(), 40)) {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                flank_deletion: FlankDeletion::NONE,
                ..
            } => assert_eq!(&*bases, BUNDLE_HULL),
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    /// A read one `CA` short in the first run, the mapper's deletion inside the hull, is
    /// spelled from its own alignment over the hull, both runs and the fixed bases between.
    #[test]
    fn a_read_with_a_contraction_in_one_run_is_spelled_over_the_hull() {
        let seq = [BUNDLE_LEFT, b"CACACA", b"TTGATGATGAT", BUNDLE_RIGHT].concat();
        let mut contracted = read(&seq, 40);
        contracted.cigar = vec![
            CigarOp::Match(8 + 6),
            CigarOp::Deletion(2),
            CigarOp::Match(11 + 8),
        ];
        match classify_in_the_bundle(&contracted) {
            Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                ..
            } => assert_eq!(&*bases, b"CACACATTGATGATGAT"),
            other => panic!("expected a complete observation, got {other:?}"),
        }
    }

    /// **A read that ends inside the hull is kept as a partial, as at a segment** — so it
    /// counts toward the hull's depth and is scored by calling (`str_bundle_loci.md`
    /// §6.1, §7.6; changed during implementation, 2026-10-11, where it was first dropped).
    #[test]
    fn a_read_ending_inside_the_hull_is_kept_as_a_partial() {
        let seq = [BUNDLE_LEFT, b"CACACACATTGA"].concat();
        match classify_in_the_bundle(&read(&seq, 40)) {
            Classified::Observed {
                read_witness: ReadWitness::Partial { .. },
                ..
            } => {}
            other => panic!("a bundle keeps the read as a partial, got {other:?}"),
        }
    }

    /// **A partial whose witness covers the whole hull carries the flank it held; one whose
    /// witness says its side does not** (`str_bundle_loci.md` §7.6; the review's finding
    /// M3). A read from the left flank showing the hull and five bases more, with no right
    /// flank to anchor, is clamped to every position of the hull, so only the recorded flank
    /// says it was read from the left.
    #[test]
    fn a_partial_over_the_whole_hull_records_the_flank_it_held() {
        let past_the_hull = [BUNDLE_LEFT, BUNDLE_HULL, b"GATGA"].concat();
        match classify_in_the_bundle(&read(&past_the_hull, 40)) {
            Classified::Observed {
                read_witness: ReadWitness::Partial { positions },
                held_flank,
                ..
            } => {
                assert_eq!(positions.first_run(), (0, BUNDLE_HULL.len() as u16));
                assert_eq!(held_flank, Some(HeldFlank::Left));
            }
            other => panic!("expected a partial over the whole hull, got {other:?}"),
        }
        let short = [BUNDLE_LEFT, b"CACACACATTGA"].concat();
        match classify_in_the_bundle(&read(&short, 40)) {
            Classified::Observed {
                read_witness: ReadWitness::Partial { .. },
                held_flank,
                ..
            } => assert_eq!(held_flank, None, "the witness says the side"),
            other => panic!("expected a partial, got {other:?}"),
        }
    }

    /// **The mirror image: a read from the right flank that ran past the hull's first base
    /// records the right flank** — five more `CA` bases where the left flank would be, so it
    /// anchors only the right flank, and its witness is clamped to every position of the hull
    /// as the left-held read's is. Recording `Left` for it, or nothing, would score it from the
    /// wrong end or as the half-and-half average (this review, 2026-10-12).
    #[test]
    fn a_partial_over_the_whole_hull_from_the_right_records_the_right_flank() {
        let past_the_hull = [&b"ACACA"[..], BUNDLE_HULL, BUNDLE_RIGHT].concat();
        let mut from_the_right = read(&past_the_hull, 40);
        from_the_right.pos = 1 + (BUNDLE_LEFT.len() - 5) as u64;
        match classify_in_the_bundle(&from_the_right) {
            Classified::Observed {
                read_witness: ReadWitness::Partial { positions },
                held_flank,
                ..
            } => {
                assert_eq!(positions.first_run(), (0, BUNDLE_HULL.len() as u16));
                assert_eq!(held_flank, Some(HeldFlank::Right));
            }
            other => panic!("expected a partial over the whole hull, got {other:?}"),
        }
    }

    /// A read lying wholly in a flank is outside the hull, as it is outside a segment —
    /// not counted as a read that ended inside.
    #[test]
    fn a_read_in_the_flank_alone_is_outside_the_hull() {
        assert_eq!(
            classify_in_the_bundle(&read(BUNDLE_LEFT, 40)),
            Classified::NoObservation(NoObservationReason::OutsideTract)
        );
    }

    /// **The edge-escape check asks about the run at each edge, and about nothing where
    /// fixed sequence reaches it** (`str_bundle_loci.md` §6.1). The left flank ending in
    /// `CA` continues a `CA` run at the hull's first base; with the run there the reference
    /// run escapes, and with fixed sequence there it does not.
    #[test]
    fn the_edge_escape_check_follows_the_run_at_each_edge() {
        let frame = [b"GGATCCCA".as_slice(), BUNDLE_HULL, BUNDLE_RIGHT].concat();
        let at = |edges| {
            let mut locus = bundle_locus(edges);
            locus.tract_with_margin_bases = frame.clone().into();
            reference_run_escapes_the_tract(&locus, &frame)
        };
        assert!(
            at(runs_at_both_edges()),
            "the CA run carries on into the flank"
        );
        assert!(
            !at(super::super::EdgeRunPeriods {
                left: None,
                right: Some(3),
            }),
            "fixed sequence at the left edge: nothing to carry on"
        );
        assert!(!at(super::super::EdgeRunPeriods {
            left: None,
            right: None,
        }));
        // The untouched frame escapes at neither edge.
        let locus = bundle_locus(runs_at_both_edges());
        assert!(!reference_run_escapes_the_tract(&locus, &bundle_frame()));
    }

    /// The mirror at the right edge alone: a right flank starting `GAT` continues the
    /// `(GAT)3` run, and escapes only when a period-3 run is said to reach that edge.
    #[test]
    fn the_edge_escape_check_asks_the_right_edge_on_its_own() {
        let frame = [BUNDLE_LEFT, BUNDLE_HULL, b"GATAGGAC"].concat();
        let at = |edges| {
            let mut locus = bundle_locus(edges);
            locus.tract_with_margin_bases = frame.clone().into();
            reference_run_escapes_the_tract(&locus, &frame)
        };
        assert!(
            at(super::super::EdgeRunPeriods {
                left: None,
                right: Some(3),
            }),
            "the GAT run carries on into the right flank"
        );
        assert!(
            !at(super::super::EdgeRunPeriods {
                left: Some(2),
                right: None,
            }),
            "fixed sequence at the right edge, and the left run does not carry on"
        );
    }
}
