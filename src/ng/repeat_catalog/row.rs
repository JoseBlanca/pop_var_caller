//! Turning one detected interval into one catalog row — the edge where the detector's
//! 0-based half-open offsets become ng's 1-based inclusive coordinates, and where the two
//! values that need the bases (the whole-motif cut and the purity over it) are computed
//! once so no reader ever has to.
//!
//! **This is the file's one silent-failure site.** Every other stage fails loudly: a wrong
//! criteria comparison refuses, a broken file will not open. An off-by-one here produces a
//! row that reads perfectly and describes the wrong bases, and every consumer inherits it.

use crate::ng::region_typing::segment_criteria::{
    SsrSegmentCriteria, minimal_trim, recompute_purity, upper,
};
use crate::ng::repeat_catalog::criteria::StrRepeatCriteria;
use crate::ng::repeat_catalog::{FoundRepeat, TractSpan};
use crate::ng::tandem_repeat::RepeatInterval;
use crate::ng::types::{Bp, ContigId, Motif, Position};

/// Why a detected interval did not become a catalog row.
///
/// Not an error type: these are the builder's own admission decisions, and each one is a
/// row the file deliberately does not hold. The builder tallies them so a run can say what
/// it dropped rather than leaving the count of what it kept unexplained.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowRejection {
    /// Fewer motif copies in the **detected** span than the catalog's floor for that
    /// period.
    ///
    /// **Measured on the detected span, as `prefilter` measures it** — deliberately not on
    /// the trimmed one. Bundling runs on detected coordinates before any trim, so a row
    /// dropped here must be a row a reader's own pre-screen would have dropped too;
    /// measuring the trimmed span instead would remove tracts that are still a neighbour's
    /// bundle member, and the neighbour would then be classified differently from a live
    /// scan (spec §3.1).
    CopyFloor,
    /// Closer to a contig's end than the catalog's flank floor, on one side or both. Such
    /// a tract cannot be anchored by a read, so no caller can use it (spec §1).
    TooCloseToContigEnd,
    /// The tract is shorter than one whole motif, so it has no motif to name. Charged
    /// separately from the copy floor because it is a malformed interval rather than a
    /// short one; the detector does not emit these.
    NoMotif,
    /// The period is outside 1..=`MAX_MOTIF_LEN`, so `Motif` cannot hold the unit. Also
    /// not reachable from a scan whose period range the criteria bound.
    PeriodOutOfRange,
}

/// Build the catalog row for one detected interval, or say why there is none.
///
/// `interval` is the detector's own output: **0-based half-open offsets into `bases`**,
/// which are `contig`'s bases from position 1 (the builder scans a contig whole, so the
/// slice starts at the contig's first base — spec §2.3). `contig_len` is that contig's own
/// length, which is what the flank floor is measured against.
///
/// What it does, in order, and the order is [`crate::ng::region_typing::segment_criteria::classify`]'s:
///
/// 1. the motif — the tract's first `period` bases, uppercased, verbatim and phase-faithful;
/// 2. the copy floor on the **detected** span (see [`RowRejection::CopyFloor`]);
/// 3. the whole-motif cut, `None` when the tract has no clean one — **the row survives**,
///    because such a tract is still a bundle member for its neighbour (spec §3.1);
/// 4. the purity over the cut tract, `None` exactly when the cut is;
/// 5. the flank floor against the contig's ends, measured on the cut tract when there is
///    one and on the detected span when there is not.
///
/// It applies **no purity floor, no satellite cap and no bundling**: those are the reader's
/// (spec §1.3).
pub fn row_for_interval(
    contig: ContigId,
    contig_len: Bp,
    bases: &[u8],
    interval: &RepeatInterval,
    criteria: &StrRepeatCriteria,
) -> Result<FoundRepeat, RowRejection> {
    let period = usize::from(interval.period);
    if interval.end <= interval.start || interval.end as usize > bases.len() {
        return Err(RowRejection::NoMotif);
    }

    let tract = upper(&bases[interval.start as usize..interval.end as usize]);
    let motif_bytes = tract.get(..period).ok_or(RowRejection::NoMotif)?.to_vec();
    let motif = Motif::new(&motif_bytes).map_err(|_| RowRejection::PeriodOutOfRange)?;

    // **1-based inclusive, once.** `[start, end)` at offset `start` into a contig-length
    // slice is `[start + 1, end]` in ng's coordinates — the only conversion in the module,
    // and the reason `TractSpan` exists as a named pair rather than two loose integers.
    let detected = TractSpan {
        start: Position(interval.start + 1),
        end: Position(interval.end),
    };

    let copy_floor = criteria
        .classification
        .min_copies
        .for_period(interval.period);
    if detected.repeat_count(interval.period) < u64::from(copy_floor) {
        return Err(RowRejection::CopyFloor);
    }

    // The cut is offsets into the tract, so it rebases by the tract's own start.
    let cut = minimal_trim(&tract, &motif_bytes);
    let trimmed = cut.map(|(start_offset, end_offset)| TractSpan {
        start: Position(interval.start + start_offset as u64 + 1),
        end: Position(interval.start + end_offset as u64),
    });
    let purity = cut.map(|(start_offset, end_offset)| {
        recompute_purity(&tract[start_offset..end_offset], &motif_bytes)
    });

    // The flank is measured against the **locus** — the cut tract where there is one,
    // since that is what a read would be anchored around.
    let anchored = trimmed.unwrap_or(detected);
    let flank = criteria.min_flank_bp.get();
    let room_left = anchored.start.get() - 1;
    let room_right = contig_len.get().saturating_sub(anchored.end.get());
    if room_left < flank || room_right < flank {
        return Err(RowRejection::TooCloseToContigEnd);
    }

    Ok(FoundRepeat {
        contig,
        detected,
        trimmed,
        period: interval.period,
        score: interval.score,
        motif,
        purity,
    })
}

/// The catalog's own copy floor, applied to a **detected** span, exactly as
/// [`crate::ng::region_typing::segment_criteria::prefilter`] applies a policy's.
///
/// Exposed so the builder can pre-screen intervals before it pays for the trim, and so the
/// property that keeps bundling identical has a name a test can pin.
#[inline]
pub fn clears_detected_copy_floor(
    interval: &RepeatInterval,
    criteria: &SsrSegmentCriteria,
) -> bool {
    let copies = (interval.end - interval.start) / u64::from(interval.period);
    copies >= u64::from(criteria.min_copies.for_period(interval.period))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A contig with `left` bases of filler, then `copies` copies of `motif`, then `right`
    /// bases of filler. Filler is `G`, which tiles nothing used here.
    fn contig_with_tract(left: usize, motif: &str, copies: usize, right: usize) -> Vec<u8> {
        let mut bases = vec![b'G'; left];
        for _ in 0..copies {
            bases.extend_from_slice(motif.as_bytes());
        }
        bases.extend(std::iter::repeat_n(b'G', right));
        bases
    }

    fn interval(start: u64, end: u64, period: u8) -> RepeatInterval {
        RepeatInterval {
            start,
            end,
            period,
            score: 42,
        }
    }

    fn criteria() -> StrRepeatCriteria {
        StrRepeatCriteria::default()
    }

    /// The conversion, on a tract whose coordinates are checkable by hand: 40 bases of
    /// filler, then `(CAG)6` at offsets 40..58, then 40 more. In ng's 1-based inclusive
    /// coordinates that tract is 41..=58.
    #[test]
    fn the_detected_span_is_the_detector_offsets_made_one_based_inclusive() {
        let bases = contig_with_tract(40, "CAG", 6, 40);
        let row = row_for_interval(
            ContigId(0),
            Bp(bases.len() as u64),
            &bases,
            &interval(40, 58, 3),
            &criteria(),
        )
        .expect("a 6-copy trimer 40 bp from either end is a row");

        assert_eq!(row.detected.start, Position(41));
        assert_eq!(row.detected.end, Position(58));
        assert_eq!(row.detected.len_bp(), 18);
        assert_eq!(row.detected.repeat_count(3), 6);
        assert_eq!(row.motif.as_bytes(), b"CAG");
        assert_eq!(row.period, 3);
        assert_eq!(row.score, 42);
    }

    /// A perfect tract needs no trimming, so the cut is the tract and the purity is 1.
    #[test]
    fn a_perfect_tract_trims_to_itself_at_purity_one() {
        let bases = contig_with_tract(40, "CAG", 6, 40);
        let row = row_for_interval(
            ContigId(0),
            Bp(bases.len() as u64),
            &bases,
            &interval(40, 58, 3),
            &criteria(),
        )
        .expect("a row");

        assert_eq!(row.trimmed, Some(row.detected));
        assert_eq!(row.purity, Some(1.0));
        assert_eq!(row.stratum(), Some((3, 6)));
    }

    /// The detector's span often overshoots by a partial copy at each end. The cut moves
    /// **inside** it, and the stratum follows the cut, not the detection.
    #[test]
    fn a_ragged_tract_cuts_back_to_whole_copies() {
        // `AG` + (CAG)6 + `CA`: the detector reports the whole 22-base run, the cut keeps
        // the 18 bases that are whole copies.
        let mut bases = vec![b'G'; 40];
        bases.extend_from_slice(b"AG");
        for _ in 0..6 {
            bases.extend_from_slice(b"CAG");
        }
        bases.extend_from_slice(b"CA");
        bases.extend(std::iter::repeat_n(b'T', 40));

        let row = row_for_interval(
            ContigId(0),
            Bp(bases.len() as u64),
            &bases,
            &interval(40, 62, 3),
            &criteria(),
        )
        .expect("a row");

        assert_eq!(
            row.detected,
            TractSpan {
                start: Position(41),
                end: Position(62)
            }
        );
        let trimmed = row.trimmed.expect("a ragged tract still has a clean cut");
        assert!(
            trimmed.start >= row.detected.start && trimmed.end <= row.detected.end,
            "the cut lies inside the detection: {trimmed:?} in {:?}",
            row.detected
        );
        assert_eq!(
            trimmed.len_bp() % 3,
            0,
            "the cut tract is whole copies: {trimmed:?}"
        );
        assert_eq!(row.stratum(), Some((3, trimmed.len_bp() / 3)));
    }

    /// A tract with no clean cut keeps its row: it is still a bundle member for its
    /// neighbour, and dropping it would change how that neighbour is classified (spec
    /// §3.1). Its stratum, its trim and its purity are all absent.
    #[test]
    fn a_tract_with_no_clean_cut_is_still_a_row() {
        // Six interrupted copies: every third copy is broken, so no two adjacent whole
        // copies exist at either end.
        let mut bases = vec![b'G'; 40];
        bases.extend_from_slice(b"CAGCATCAGCATCAGCAT");
        bases.extend(std::iter::repeat_n(b'T', 40));

        let row = row_for_interval(
            ContigId(0),
            Bp(bases.len() as u64),
            &bases,
            &interval(40, 58, 3),
            &criteria(),
        )
        .expect("an untrimmable tract is still recorded");

        assert_eq!(row.trimmed, None);
        assert_eq!(row.purity, None);
        assert_eq!(row.stratum(), None);
        assert_eq!(row.detected.repeat_count(3), 6, "the detection is intact");
    }

    /// The floor is measured on the **detected** span. A tract that clears it there but
    /// whose cut falls below it is kept — the rule that keeps bundling identical to a live
    /// scan (spec §3.1).
    #[test]
    fn a_tract_below_the_floor_after_trimming_is_kept() {
        // `(CAG)3 CAT`: 12 detected bases at period 3 is 4 copies, exactly the catalog's
        // floor, so the row is admitted. The cut ends after the third whole copy, leaving
        // 9 bases — **3 copies, below that floor** — and the row still stands, which is the
        // rule: the floor is measured where `prefilter` measures it, so that a tract this
        // one might bundle with sees the same neighbour a live scan sees.
        let mut bases = vec![b'G'; 40];
        bases.extend_from_slice(b"CAGCAGCAGCAT");
        bases.extend(std::iter::repeat_n(b'G', 40));

        let row = row_for_interval(
            ContigId(0),
            Bp(bases.len() as u64),
            &bases,
            &interval(40, 52, 3),
            &criteria(),
        )
        .expect("4 detected copies clear the period-3 floor of 4");

        let trimmed = row.trimmed.expect("a clean cut exists");
        assert_eq!(row.detected.repeat_count(3), 4, "detected: 4 copies");
        assert_eq!(trimmed.repeat_count(3), 3, "trimmed: 3 copies");
        assert!(
            u64::from(criteria().classification.min_copies.for_period(3)) > trimmed.repeat_count(3),
            "the fixture must actually fall below the floor after trimming, or it proves nothing"
        );
    }

    /// A tract below the floor **on the detected span** is not a row at all.
    #[test]
    fn a_tract_below_the_detected_copy_floor_is_rejected() {
        // Period 3's catalog floor is 4 copies; offer 3.
        let bases = contig_with_tract(40, "CAG", 3, 40);
        assert_eq!(
            row_for_interval(
                ContigId(0),
                Bp(bases.len() as u64),
                &bases,
                &interval(40, 49, 3),
                &criteria()
            ),
            Err(RowRejection::CopyFloor)
        );
    }

    /// The flank floor is 15 bp on each side, so 15 is kept and 14 is not — the pair that
    /// pins the boundary rather than merely exercising it (spec §1, §4.1).
    #[test]
    fn fifteen_bases_of_flank_is_kept_and_fourteen_is_not() {
        for (left, expected_kept) in [(15usize, true), (14, false)] {
            let bases = contig_with_tract(left, "CAG", 6, 40);
            let start = left as u64;
            let outcome = row_for_interval(
                ContigId(0),
                Bp(bases.len() as u64),
                &bases,
                &interval(start, start + 18, 3),
                &criteria(),
            );
            assert_eq!(
                outcome.is_ok(),
                expected_kept,
                "{left} bp of left flank: expected kept={expected_kept}, got {outcome:?}"
            );
            if !expected_kept {
                assert_eq!(outcome, Err(RowRejection::TooCloseToContigEnd));
            }
        }
    }

    /// The same boundary at the contig's right end, which is the side measured against the
    /// contig's length rather than against position 1.
    #[test]
    fn the_right_flank_is_measured_against_the_contigs_own_end() {
        for (right, expected_kept) in [(15usize, true), (14, false)] {
            let bases = contig_with_tract(40, "CAG", 6, right);
            let outcome = row_for_interval(
                ContigId(0),
                Bp(bases.len() as u64),
                &bases,
                &interval(40, 58, 3),
                &criteria(),
            );
            assert_eq!(
                outcome.is_ok(),
                expected_kept,
                "{right} bp of right flank: expected kept={expected_kept}, got {outcome:?}"
            );
        }
    }

    /// Soft-masked reference sequence is real (tomato SL4.0 carries 227,170 lowercase
    /// bases), and a motif read verbatim off it would compare unequal to everything.
    #[test]
    fn a_soft_masked_tract_yields_an_upper_case_motif() {
        let bases = contig_with_tract(40, "cag", 6, 40);
        let row = row_for_interval(
            ContigId(0),
            Bp(bases.len() as u64),
            &bases,
            &interval(40, 58, 3),
            &criteria(),
        )
        .expect("a row");

        assert_eq!(row.motif.as_bytes(), b"CAG");
        assert_eq!(row.purity, Some(1.0), "purity is measured on the same case");
    }

    /// The pre-screen and the row builder must agree about the copy floor, or the builder
    /// would pay for trims on intervals it then drops — and, worse, the two could disagree
    /// about which rows exist.
    #[test]
    fn the_pre_screen_agrees_with_the_row_builder_on_the_floor() {
        let criteria = criteria();
        let bases = contig_with_tract(40, "CAG", 3, 40);
        let short = interval(40, 49, 3);
        assert!(!clears_detected_copy_floor(
            &short,
            &criteria.classification
        ));
        assert_eq!(
            row_for_interval(
                ContigId(0),
                Bp(bases.len() as u64),
                &bases,
                &short,
                &criteria
            ),
            Err(RowRejection::CopyFloor)
        );

        let bases = contig_with_tract(40, "CAG", 6, 40);
        let long = interval(40, 58, 3);
        assert!(clears_detected_copy_floor(&long, &criteria.classification));
        assert!(
            row_for_interval(
                ContigId(0),
                Bp(bases.len() as u64),
                &bases,
                &long,
                &criteria
            )
            .is_ok()
        );
    }
}
