//! The read-region mapping: from a read's CIGAR, the slice of the read that covers a locus
//! window, and the base-quality gate a delimited tract must pass.
//!
//! Adapted from the STR pileup of the caller deleted on 2026-09-13, to take window coordinates
//! directly and a plain scratch buffer. Used by [`classify`](super::classify).

use crate::bam::alignment_input::CigarOp;
use std::ops::Range;

/// The lower-quartile base-quality floor a delimited tract must clear (production's
/// `MIN_REGION_Q1`).
pub(super) const MIN_REGION_Q1: u8 = 15;

/// A read's reference footprint from its CIGAR — the aligned span (0-based, half-open)
/// and the soft-clip on each end (the optimistic long-allele reach).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Footprint {
    pub(super) ref_start: u32,
    pub(super) ref_end: u32,
    pub(super) leading_clip: u32,
    pub(super) trailing_clip: u32,
}

/// First soft-clip length at an end (skipping a hard-clip), or `0` if the end is aligned.
fn end_soft_clip(op: &CigarOp) -> Option<u32> {
    match op {
        CigarOp::HardClip(_) => None,
        CigarOp::SoftClip(n) => Some(*n),
        _ => Some(0),
    }
}

/// The read's reference footprint from its CIGAR and 1-based mapping position.
pub(super) fn read_footprint(cigar: &[CigarOp], pos: u64) -> Footprint {
    let ref_start = pos.saturating_sub(1) as u32;
    let ref_span: u32 = cigar
        .iter()
        .map(|op| match op {
            CigarOp::Match(n)
            | CigarOp::Deletion(n)
            | CigarOp::Skip(n)
            | CigarOp::SeqMatch(n)
            | CigarOp::SeqMismatch(n) => *n,
            _ => 0,
        })
        .sum();
    Footprint {
        ref_start,
        ref_end: ref_start + ref_span,
        leading_clip: cigar.iter().find_map(end_soft_clip).unwrap_or(0),
        trailing_clip: cigar.iter().rev().find_map(end_soft_clip).unwrap_or(0),
    }
}

/// Read coordinate of a 0-based reference position `target` inside the aligned span. A
/// dual-cursor CIGAR walk; a `target` inside a deletion maps to the read position the
/// deletion sits at.
fn ref_to_read(cigar: &[CigarOp], ref_start: u32, leading_clip: u32, target: u32) -> usize {
    let mut ref_cur = ref_start;
    let mut read_cur = leading_clip as usize;
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                if target < ref_cur + n {
                    return read_cur + (target - ref_cur) as usize;
                }
                ref_cur += n;
                read_cur += *n as usize;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                if target < ref_cur + n {
                    return read_cur;
                }
                ref_cur += n;
            }
            CigarOp::Insertion(n) => read_cur += *n as usize,
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }
    read_cur
}

/// The read-coordinate span covering the locus window `[window_start, window_start +
/// window_len)` (0-based reference), where the whole soft-clip is included on a side the
/// window opens past the aligned span — that is where a long allele's extra tract lives,
/// and the aligner realigns within the grabbed bases (spec §2).
pub(super) fn extract_region(
    cigar: &[CigarOp],
    fp: Footprint,
    read_len: usize,
    window_start: u32,
    window_len: u32,
) -> Range<usize> {
    let w_end = window_start + window_len;
    let r_start = if window_start < fp.ref_start {
        0 // window opens left of the alignment → take the full leading clip
    } else {
        ref_to_read(cigar, fp.ref_start, fp.leading_clip, window_start)
    };
    let r_end = if w_end > fp.ref_end {
        read_len // window closes right of the alignment → take the full trailing clip
    } else {
        ref_to_read(cigar, fp.ref_start, fp.leading_clip, w_end)
    };
    // Clamp into `[0, read_len]` — defense-in-depth against a length-inconsistent CIGAR.
    let r_start = r_start.min(read_len);
    r_start..r_end.min(read_len).max(r_start)
}

/// The read-coordinate span of the tract per the read's **own input alignment** — what
/// the mapper said the read's account of the tract's reference positions is.
///
/// # Why a complete observation is spelled from here and not from the delimiter
///
/// The delimiter re-aligns the window under a stutter-priced model, and that model makes
/// a whole-unit slip far cheaper (≈ −2.9 nats at the shipped rates) than an out-of-frame
/// gap (≈ −10.4), so wherever the true variant sits at a junction — a real flank indel
/// beside the tract, or a unit of a neighbouring repeat left-aligned into it — the
/// winning re-alignment re-spells the read toward a pure motif run and the non-motif
/// evidence is destroyed. Measured on the HG002 tandem-repeat benchmark at 30×: 232 of
/// the 463 tracts whose true sequence was never offered as a candidate miss a spelling
/// that is not a pure motif run (tract-accuracy program, L4;
/// `doc/devel/ng/research/tract_accuracy_program_report.md`). The input aligner has no
/// cheap-slip incentive, so its account keeps interruptions and junction variation; the
/// delimiter keeps the jobs it is uniquely right for — anchoring, complete-versus-partial,
/// the widen-and-retry for clipped long alleles — and rules *whether* the read measures
/// the tract; this function only rules *what the measured tract says*.
///
/// # The junction convention
///
/// An inserted run at a junction belongs to the tract **iff left-alignment cannot move it
/// across that junction** — the same left-align-first convention the truth comparison
/// uses, applied symmetrically:
///
/// - at the **left** junction (anchored between the flank's last base and the tract's
///   first), the insertion is the tract's unless its last base equals the flank base
///   beside it — then left-alignment carries it into the flank and it is flank content;
/// - at the **right** junction, the insertion is the flank's unless its last base equals
///   the tract's last base — then left-alignment carries it into the tract. This is how a
///   simulator or mapper that *right*-anchors an expansion (`M I M`) still measures at the
///   expanded length.
///
/// Insertions strictly inside the tract are the tract's wherever the mapper put them; a
/// deletion needs no rule, because it contributes no bases either way.
///
/// `None` when the read's aligned footprint does not bracket the tract plus one base on
/// each side — a clipped or partial read — and the caller falls back to the delimiter's
/// own spelling, which is the only account such a read has. **The caller also asks only
/// for reads whose CIGAR carries an indel**: an all-Match CIGAR is how a mapper spells a
/// collapsed long allele (the length change hidden as flank mismatches), and rescuing
/// that shape is the delimiter's one irreplaceable job — `classify_read` states the rule
/// beside the call.
pub(super) fn tract_span_from_input(
    cigar: &[CigarOp],
    seq: &[u8],
    fp: Footprint,
    tract_start: u32,
    tract_len: u32,
    flank_base_before: u8,
    tract_last_base: u8,
) -> Option<Range<usize>> {
    let tract_end = tract_start + tract_len; // exclusive
    if fp.ref_start + 1 > tract_start || fp.ref_end < tract_end + 1 {
        return None;
    }
    let mut ref_cur = fp.ref_start;
    let mut read_cur = fp.leading_clip as usize;
    let mut start: Option<usize> = None;
    let mut end: Option<usize> = None;
    for op in cigar {
        match op {
            CigarOp::Match(n) | CigarOp::SeqMatch(n) | CigarOp::SeqMismatch(n) => {
                if start.is_none() && (ref_cur..ref_cur + n).contains(&tract_start) {
                    start = Some(read_cur + (tract_start - ref_cur) as usize);
                }
                if end.is_none() && tract_end > ref_cur && tract_end <= ref_cur + n {
                    end = Some(read_cur + (tract_end - ref_cur) as usize);
                }
                ref_cur += n;
                read_cur += *n as usize;
            }
            CigarOp::Deletion(n) | CigarOp::Skip(n) => {
                if start.is_none() && (ref_cur..ref_cur + n).contains(&tract_start) {
                    start = Some(read_cur);
                }
                if end.is_none() && tract_end > ref_cur && tract_end <= ref_cur + n {
                    end = Some(read_cur);
                }
                ref_cur += n;
            }
            CigarOp::Insertion(n) => {
                let count = *n as usize;
                let last_inserted = seq.get(read_cur + count - 1).copied()?;
                if start.is_none() && ref_cur == tract_start && last_inserted != flank_base_before {
                    start = Some(read_cur);
                }
                if ref_cur == tract_end
                    && end == Some(read_cur)
                    && crate::locus_generation::right_junction_insertion_extends_the_repeat(
                        last_inserted,
                        tract_last_base,
                    )
                {
                    end = Some(read_cur + count);
                }
                read_cur += count;
            }
            CigarOp::SoftClip(_) | CigarOp::HardClip(_) | CigarOp::Padding(_) => {}
        }
    }
    match (start, end) {
        (Some(from), Some(to)) if from <= to && to <= seq.len() => Some(from..to),
        _ => None,
    }
}

/// Whether a delimited tract looks truncated by the extraction window rather than by the
/// read genuinely ending — the long-allele recovery trigger. A side is window-bounded when
/// [`extract_region`] did not reach the read edge there; the tract is suspicious when a
/// window-bounded side carries fewer flank bytes than the locus declares. `tract` is
/// region-relative.
pub(super) fn flank_truncated(
    region: &Range<usize>,
    tract: &Range<usize>,
    read_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
) -> bool {
    let region_len = region.end - region.start;
    let left_flank_bytes = tract.start;
    let right_flank_bytes = region_len - tract.end;
    let left_window_bounded = region.start > 0;
    let right_window_bounded = region.end < read_len;
    (left_window_bounded && left_flank_bytes < left_flank_len)
        || (right_window_bounded && right_flank_bytes < right_flank_len)
}

/// Widen a read-coordinate region by one full reference flank each side, clamped to the
/// read — in read, not reference, coordinates, sidestepping the unreliable CIGAR of a
/// mis-aligned long-allele read.
pub(super) fn widen_region(
    region: Range<usize>,
    read_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
) -> Range<usize> {
    let start = region.start.saturating_sub(left_flank_len);
    let end = (region.end + right_flank_len).min(read_len);
    start..end
}

/// Whether the tract's base qualities clear the gate: the nearest-rank lower quartile (the
/// element at sorted index `⌊(n-1)/4⌋`) is at least `threshold`. `buffer` is reused
/// scratch. An empty tract passes vacuously.
pub(super) fn passes_quality_gate(quals: &[u8], threshold: u8, buffer: &mut Vec<u8>) -> bool {
    if quals.is_empty() {
        return true;
    }
    buffer.clear();
    buffer.extend_from_slice(quals);
    let k = (buffer.len() - 1) / 4;
    let (_, q1, _) = buffer.select_nth_unstable(k);
    *q1 >= threshold
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A window fully inside an all-Match read maps to the same read offsets as reference
    /// offsets (shifted by the mapping position).
    #[test]
    fn extract_region_maps_a_window_inside_an_all_match_read() {
        // read at pos 11 (0-based 10), 30M; window [15, 25) (0-based).
        let cigar = vec![CigarOp::Match(30)];
        let fp = read_footprint(&cigar, 11);
        let region = extract_region(&cigar, fp, 30, 15, 10);
        // ref 15 → read 5 (15 - 10), ref 25 → read 15.
        assert_eq!(region, 5..15);
    }

    /// An internal deletion shifts the window's right edge left in read coordinates.
    #[test]
    fn extract_region_maps_across_an_internal_deletion() {
        // pos 11 (0-based 10): 10M 4D 16M — ref span 30, read len 26. Window [15, 28).
        let cigar = vec![CigarOp::Match(10), CigarOp::Deletion(4), CigarOp::Match(16)];
        let fp = read_footprint(&cigar, 11);
        assert_eq!(fp.ref_end, 40);
        let region = extract_region(&cigar, fp, 26, 15, 13);
        // ref 15 → read 5; ref 28 is past the deletion (ref 20..24 deleted), read = 10 + (28-24) = 14.
        assert_eq!(region, 5..14);
    }

    /// A window opening left of the alignment takes the full leading soft-clip.
    #[test]
    fn extract_region_takes_the_leading_softclip_when_the_window_opens_left() {
        // pos 21 (0-based 20): 5S 20M. Window [18, 30) opens left of ref_start 20.
        let cigar = vec![CigarOp::SoftClip(5), CigarOp::Match(20)];
        let fp = read_footprint(&cigar, 21);
        assert_eq!(fp.leading_clip, 5);
        let region = extract_region(&cigar, fp, 25, 18, 12);
        assert_eq!(region.start, 0, "the full leading clip is grabbed");
    }

    #[test]
    fn flank_truncated_flags_a_window_bounded_short_flank_only() {
        // region [3, 40) within a 50-base read (both sides window-bounded), tract [5, 32)
        // → left flank 5, right flank region_len(37) - 32 = 5.
        assert!(
            flank_truncated(&(3..40), &(5..32), 50, 6, 6),
            "a flank of 5 below the declared 6 is truncated"
        );
        assert!(
            !flank_truncated(&(3..40), &(5..32), 50, 5, 5),
            "flanks exactly the declared length are not truncated"
        );
        // A region reaching both read edges is never flagged (genuine allele ≥ read length).
        assert!(!flank_truncated(&(0..50), &(0..45), 50, 6, 6));
    }

    #[test]
    fn widen_region_extends_a_flank_each_side_clamped_to_the_read() {
        assert_eq!(widen_region(20..40, 100, 10, 10), 10..50);
        // Clamp at both ends.
        assert_eq!(widen_region(5..95, 100, 10, 10), 0..100);
    }

    #[test]
    fn the_quality_gate_keys_on_the_nearest_rank_lower_quartile() {
        let mut buffer = Vec::new();
        assert!(
            passes_quality_gate(&[], MIN_REGION_Q1, &mut buffer),
            "empty passes"
        );
        // len 1 → k = 0 (the lowest); 15 passes, 14 fails.
        assert!(passes_quality_gate(&[15], MIN_REGION_Q1, &mut buffer));
        assert!(!passes_quality_gate(&[14], MIN_REGION_Q1, &mut buffer));
        // len 8 → k = ⌊7/4⌋ = 1 (the 2nd-lowest): one base below the floor still passes.
        assert!(passes_quality_gate(
            &[10, 20, 20, 20, 20, 20, 20, 20],
            MIN_REGION_Q1,
            &mut buffer
        ));
        assert!(!passes_quality_gate(
            &[10, 10, 20, 20, 20, 20, 20, 20],
            MIN_REGION_Q1,
            &mut buffer
        ));
    }
}
