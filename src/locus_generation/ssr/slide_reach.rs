//! **How far an indel inside a repeat can slide past the repeat's edges**, and what that does to
//! the reads measured across it.
//!
//! A tract is the stretch region typing drew, one motif from edge to edge. An indel inside it can
//! sometimes be written further along than the tract reaches: at GIAB HG002 chr22:18094763 the
//! tract is nine `T`s, the reference after it is `CTTTTCTTTTCTTTTCTTC`, and inserting `TTTTC`
//! after the fifth `T` gives the same haplotype as inserting it anywhere up to 18 bases past the
//! tract's last base. The array that carries it is too short for the catalog to list, so nothing
//! joins it to the tract.
//!
//! Two things then go wrong, both measured on simulated reads from that site's homozygote at 200×
//! (every read carrying the insertion), where the tract came out `0/1` at GQ 99 with 47 reads
//! counted as the reference:
//!
//! - **a read that stops inside the array is the reference's bases**, whichever haplotype it came
//!   from, and its flank matches, so it was measured as a complete reference observation (18 of
//!   the 47, called on their own);
//! - **a mapper that wrote the insertion inside the array**, rather than at its left-most place
//!   in the tract, gave a read whose own alignment spells the tract as the reference, and that
//!   spelling is what a complete observation keeps (17 of the 47, called on their own).
//!
//! This module answers the first from the reference alone, so every sample of a cohort gets the
//! same answer at a locus ([`SlideReach`], [`demoted_if_short_of_the_slide`]): all 18 are gone.
//!
//! **The second is not answered.** Moving such a read's indel to its left-most placement before
//! the tract is spelled was built and measured, and rescued none of the 17: each carries a
//! sequencing error inside the inserted copy, which pins the mapper's placement where no exact move
//! can reach. On the GIAB per-sample benchmark at 30× and 300× and on 63 tomato accessions it
//! changed no byte of the calls, so it was removed (`doc/devel/reports/
//! ng_ref_reads_in_repeats_2026-10-17.md` §3.1).

use super::SsrLocus;
use super::classify::{Classified, NoObservationReason};
use super::read_region::read_footprint;
use crate::locus_generation::{HeldFlank, LocusLen, ReadWitness};
use crate::read::aligned_read::AlignedRead;

/// The longest indel unit followed. Six is the longest motif region typing calls a short tandem
/// repeat; a longer unit sliding out of a tract is rarer than anything this measures.
pub(super) const MAX_INDEL_UNIT: usize = 6;

/// The furthest a slide is followed past an edge, in bases. Only keeps the reference fetch
/// bounded; at GRCh38 chr22 no tract's slide reaches 100 bases past its edge.
pub(super) const MAX_SLIDE: usize = 100;

/// How many bases past each edge of a tract an indel starting inside it can still be placed.
/// Zero on a side means no indel in the tract can be written past that edge.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct SlideReach {
    pub(super) left: u64,
    pub(super) right: u64,
}

impl SlideReach {
    /// Measured on `window`, which holds the tract at `first..=last` (0-based) and the bases
    /// around it — [`MAX_SLIDE`] + [`MAX_INDEL_UNIT`] either side, or fewer at a contig end.
    pub(super) fn measure(window: &[u8], first: usize, last: usize) -> Self {
        let right = past_the_right_edge(window, first, last);
        let mirrored: Vec<u8> = window.iter().rev().copied().collect();
        let n = window.len();
        let left = past_the_right_edge(&mirrored, n - 1 - last, n - 1 - first);
        Self {
            left: left as u64,
            right: right as u64,
        }
    }

    pub(super) fn is_none(self) -> bool {
        self.left == 0 && self.right == 0
    }
}

/// Bases past `last` that an indel with its first changed base inside `first..=last` can reach.
///
/// An indel of `p` bases placed after base `a` covers `a+1 ..= a+p`, and can be moved one base
/// right for every following base that equals the base `p` before it. So past the tract it
/// reaches the end of a run of such repeats that starts at most `p` bases after the tract —
/// `a+1` being inside it. **A run is counted only when it holds for at least `p` bases**, one more
/// whole copy: a single base that happens to equal the one `p` back would otherwise give almost
/// every tract a reach, and the indel that would use it — a duplication of the flank's next bases
/// placed across the junction — is the flank's by the junction convention, not the tract's.
fn past_the_right_edge(window: &[u8], first: usize, last: usize) -> usize {
    let mut reach = 0;
    for unit in 1..=MAX_INDEL_UNIT {
        let earliest = (last + 1).max(first + unit + 1);
        for start in earliest..=last + unit {
            if start >= window.len() {
                break;
            }
            // Only the start of a run: an earlier start already followed it.
            if start > last + 1 && window[start - 1] == window[start - 1 - unit] {
                continue;
            }
            let mut j = start;
            while j < window.len()
                && j <= last + MAX_SLIDE
                && window[j] != b'N'
                && window[j] == window[j - unit]
            {
                j += 1;
            }
            if j - start >= unit {
                reach = reach.max(j - 1 - last);
            }
        }
    }
    reach
}

/// **A complete reference observation from a read that stopped inside the slide becomes a partial
/// one**, held from the flank it did reach; every other outcome is returned as it came.
///
/// A read has to show at least one base past where the indel could still sit before it can tell a
/// carrier from the reference. Short of that, and spelling the tract as the reference, it is
/// exactly what a carrier's read would show. **Only reference spellings are demoted**: a read that
/// stopped in the same place but spelled the tract some other way — one copy short, say — saw the
/// tract end where a carrier's would not, and is evidence against both.
///
/// What the demotion costs is the read's evidence *against longer alleles of the tract itself*: as
/// a partial it says "at least this many bases", so it no longer counts against an expansion of
/// the tract. It still counts against contractions.
pub(super) fn demoted_if_short_of_the_slide(
    outcome: Classified,
    read: &AlignedRead,
    locus: &SsrLocus,
    reach: SlideReach,
) -> Classified {
    if reach.is_none() {
        return outcome;
    }
    let Classified::Observed {
        bases,
        read_witness: ReadWitness::Complete,
        flank_deletion,
        held_flank: None,
        q_sum,
    } = outcome
    else {
        return outcome;
    };
    let left = locus.left_flank_len();
    let tract_len = locus.repeat.tract_len() as usize;
    let reference_tract = &locus.tract_with_margin_bases[left..left + tract_len];
    if !flank_deletion.is_none() || &*bases != reference_tract {
        return Classified::Observed {
            bases,
            read_witness: ReadWitness::Complete,
            flank_deletion,
            held_flank: None,
            q_sum,
        };
    }
    let footprint = read_footprint(&read.cigar, read.pos);
    let first_aligned = u64::from(footprint.ref_start) + 1;
    let last_aligned = u64::from(footprint.ref_end);
    let short_right = reach.right > 0 && last_aligned < locus.repeat.end() + reach.right + 1;
    let short_left = reach.left > 0 && first_aligned + reach.left + 1 > locus.repeat.start();
    let locus_len = LocusLen::from_positions(locus.repeat.tract_len());
    let covered = u16::try_from(tract_len).unwrap_or(u16::MAX);
    let (read_witness, held) = match (short_left, short_right) {
        (false, false) => {
            return Classified::Observed {
                bases,
                read_witness: ReadWitness::Complete,
                flank_deletion,
                held_flank: None,
                q_sum,
            };
        }
        (true, true) => return Classified::NoObservation(NoObservationReason::WindowTruncated),
        (false, true) => (ReadWitness::from_left(covered, locus_len), HeldFlank::Left),
        (true, false) => (
            ReadWitness::from_right(covered, locus_len),
            HeldFlank::Right,
        ),
    };
    let Some(read_witness) = read_witness else {
        return Classified::NoObservation(NoObservationReason::WindowTruncated);
    };
    Classified::Observed {
        bases,
        read_witness,
        flank_deletion,
        held_flank: Some(held),
        q_sum,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::CigarOp;
    use crate::locus_generation::FlankDeletion;
    use crate::region_typing::segment_criteria::{Motif, SsrSegment};
    use crate::types::{Position, ReadGroupId};

    /// GRCh38 chr22:18094730-18094820. The tract is the nine `T`s at 0-based 29..=37; inserting
    /// `TTTTC` inside it slides through the `(CTTTT)` array to 18094785, 18 bases past the tract.
    const HG002_CHR22_18094763: &[u8] =
        b"CCATAGCTAGTTTCTTATTTTCCTTTGGCTTTTTTTTTCTTTTCTTTTCTTTTCTTCTTTTTTTTTTTGAGATGGAATCTTGCTCTTGTCA";

    #[test]
    fn a_tract_followed_by_an_array_too_short_to_list_reaches_to_the_arrays_end() {
        let reach = SlideReach::measure(HG002_CHR22_18094763, 29, 37);
        assert_eq!(reach.right, 18);
        assert_eq!(reach.left, 0);
    }

    /// The STR locus at 18094763: the tract at 1-based 30..=38 of [`HG002_CHR22_18094763`], with
    /// 15-base flanks — the slide reaches 18 past the tract, to 56.
    fn locus_at_18094763() -> SsrLocus {
        SsrLocus {
            repeat: SsrSegment::new("chr1".into(), 30, 38, Motif::new(b"T").unwrap(), 1.0)
                .expect("a valid segment")
                .into(),
            tract_with_margin_bases: HG002_CHR22_18094763[14..53].into(),
            margin_start: Position(15),
        }
    }

    /// A read aligned from 1 to `last`, matching the reference.
    fn read_to(last: u64) -> AlignedRead {
        let seq = HG002_CHR22_18094763[..last as usize].to_vec();
        AlignedRead {
            qname: b"r".to_vec(),
            flag: 0,
            ref_id: 0,
            pos: 1,
            mapq: 60,
            cigar: vec![CigarOp::Match(last as u32)],
            qual: vec![30; seq.len()],
            seq,
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            read_group: ReadGroupId(0),
        }
    }

    fn complete(bases: &[u8]) -> Classified {
        Classified::Observed {
            bases: bases.into(),
            read_witness: ReadWitness::Complete,
            flank_deletion: FlankDeletion::NONE,
            held_flank: None,
            q_sum: -1.0,
        }
    }

    /// **A reference read that stops inside the slide becomes a partial held from the left**;
    /// one that shows a base past the slide stays complete.
    #[test]
    fn a_reference_read_stopping_inside_the_slide_is_a_partial() {
        let locus = locus_at_18094763();
        let reach = SlideReach::measure(HG002_CHR22_18094763, 29, 37);
        let tract = b"TTTTTTTTT";

        let stopped = demoted_if_short_of_the_slide(complete(tract), &read_to(50), &locus, reach);
        let Classified::Observed {
            read_witness,
            held_flank,
            ..
        } = stopped
        else {
            panic!("still an observation: {stopped:?}");
        };
        assert!(!matches!(read_witness, ReadWitness::Complete));
        assert_eq!(held_flank, Some(HeldFlank::Left));

        // One base past the slide's last, 56, is enough.
        assert_eq!(
            demoted_if_short_of_the_slide(complete(tract), &read_to(57), &locus, reach),
            complete(tract),
        );
        assert!(!matches!(
            demoted_if_short_of_the_slide(complete(tract), &read_to(56), &locus, reach),
            Classified::Observed {
                read_witness: ReadWitness::Complete,
                ..
            }
        ));
    }

    /// **A read that spells the tract otherwise is not demoted**, however short: it saw the tract
    /// end where a carrier's read would not.
    #[test]
    fn a_non_reference_read_stopping_inside_the_slide_stays_complete() {
        let locus = locus_at_18094763();
        let reach = SlideReach::measure(HG002_CHR22_18094763, 29, 37);
        let one_short = b"TTTTTTTT";
        assert_eq!(
            demoted_if_short_of_the_slide(complete(one_short), &read_to(50), &locus, reach),
            complete(one_short),
        );
    }

    #[test]
    fn a_tract_in_plain_sequence_reaches_nowhere() {
        //                 tract: AAAAAAAA at 6..=13
        let window = b"CGTCAGAAAAAAAACGTGCATCG";
        assert!(SlideReach::measure(window, 6, 13).is_none());
    }

    #[test]
    fn one_coincidental_base_after_the_tract_is_not_a_reach() {
        // `TTTTTTTT` then `GAGC`: the second `G` repeats the first, two bases back, once — and
        // past where a two-base indel from inside the tract could start.
        let window = b"CAGCTTTTTTTTGAGCACGTA";
        assert_eq!(SlideReach::measure(window, 4, 11).right, 0);
    }
}
