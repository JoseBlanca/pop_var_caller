//! **Partial reads that single out one sequence count toward candidate selection** at STR
//! segments and STR bundles (`doc/devel/ng/spec/candidate_alleles_ssr.md` §8.1, adopted
//! 2026-10-13).
//!
//! **The rule in plain words.** A read that ended inside the repeat cannot spell an allele, but
//! it can rule alleles out: its bases, read from the flank it held, either fit an allele (the
//! allele, then the reference flank past it, starts with them) or do not. **When the read fits
//! exactly one of the sequences in the locus's table, it counts as one more read of that
//! sequence**, for its own sample, in both halves of the bar — the sequence's reads and the
//! sample's total. A partial read that fits two sequences or more says nothing about which, and
//! one that fits none (a sequencing error, or an allele no crossing read spelled) says nothing
//! either; both count nowhere.
//!
//! **A read need not reach the flank to fit one sequence alone.** One that ran past the end of
//! every shorter sequence into the flank fits only the sequences of the length it crossed; but
//! one that ended inside the repeat, longer than every sequence of the table but one, fits that
//! one only, and counts for it. The table is all the rule knows: a longer allele no crossing read
//! spelled would fit the read as well, and is not there to be fitted (spec §8.1's first limit).
//!
//! **Why this lets alleles in without letting stutter in.** A partial read singles out one
//! sequence about as often for a stutter product as for the allele it came from, and because the
//! deciding partials enter the denominator as well as the numerator, a stutter sequence's share
//! of a sample's reads is unchanged at depth, where the share decides. What changes is the floor
//! of two reads at low depth: one crossing read and one deciding partial now reach it.
//!
//! **Precisely.** For a partial read of sample `s` with bases `w` that held flank `h`
//! ([`HeldFlank::of_witness`]), and an allele `a` of the merge's table with flank deletion
//! `(dl, dr)`:
//!
//! ```text
//! fits(w, a)  ⟺  h = Left  and  a ++ right_flank[dr..]  starts with w
//!              or  h = Right and  left_flank[..len − dl] ++ a  ends with w
//! the read counts for a  ⟺  a is the only allele of the table that w fits
//! ```
//!
//! and each such read is added to `s`'s support row for `(a, the read's group)` before
//! selection runs. A read whose side neither its witness nor the classifier recorded counts for
//! nothing.
//!
//! **Selection's bar only.** The copy with the counted reads feeds nomination, admission and the
//! cap's "this sample earned it" test; the leftover's read count and error mass, which reach the
//! record's `DP`, are taken from the merge's own rows ([`leftover_of`](super::leftover_of)), and
//! the read likelihood is handed the merge's own rows and partials, so a partial read is never
//! scored as complete and never counted twice.

use crate::locus_generation::{FlankDeletion, HeldFlank, LocusLen, ReadWitness};
use crate::run::cohort_merge::build::{
    AlleleSupport, CohortObservation, PartialObservation, SampleSupport, SupportedAllele,
};

/// **The locus with every deciding partial read added to its sample's support rows**, or `None`
/// where no partial read of any sample singles out one allele — the caller then selects over the
/// observation as it is, without a copy.
///
/// The copy's samples carry **no partial reads**: the ones that counted are in the support rows
/// now, and selection reads nothing else of them, so carrying them would cost one allocation per
/// partial read and leave each counted read in the copy twice.
///
/// `flanks` are the reference flanks either side of the segment or the hull — the bases a read
/// that ran past the locus's far end shows, and that the bundle read model compares it with too.
pub(super) fn count_deciding_partial_reads(
    observation: &CohortObservation,
    flanks: (&[u8], &[u8]),
) -> Option<CohortObservation> {
    // **No counted reads where some allele is a deletion across the locus's edge.** A partial
    // read records no flank deletion, so a carrier's partial read is spelled as if its flank
    // were whole, and it can fit exactly a crossing read that the delimiter spelled the same
    // wrong way. Measured on the stuttering benchmark (50 samples, PCR, 30×): one `AGC` segment
    // whose edge deletion was carried by 18 samples; 85 partial reads of carriers singled out a
    // misspelling shown by one crossing read, admitted it, and turned all 18 right calls wrong.
    if observation
        .allele_flank_deletions
        .iter()
        .any(|deletion| !deletion.is_none())
    {
        return None;
    }
    // The axis the witnesses are on, derived as the read likelihood derives it.
    let locus_len = LocusLen::of_region(observation.region);
    let deciding_reads: Vec<Vec<(usize, &PartialObservation)>> = observation
        .per_sample
        .iter()
        .map(|sample| {
            sample
                .partials
                .iter()
                .filter_map(|partial| {
                    the_one_allele_a_partial_read_fits(observation, partial, flanks, locus_len)
                        .map(|allele| (allele, partial))
                })
                .collect()
        })
        .collect();
    if deciding_reads.iter().all(Vec::is_empty) {
        return None;
    }
    let per_sample = observation
        .per_sample
        .iter()
        .zip(&deciding_reads)
        .map(|(sample, deciding)| {
            let SampleSupport {
                sample,
                supported,
                partials: _,
                reads_without_observation,
                reads_removed_as_evidence,
                reads_composed_across_records,
            } = sample;
            let mut supported = supported.clone();
            add_deciding_reads(&mut supported, deciding);
            SampleSupport {
                sample: *sample,
                supported,
                partials: Vec::new(),
                reads_without_observation: *reads_without_observation,
                reads_removed_as_evidence: *reads_removed_as_evidence,
                reads_composed_across_records: *reads_composed_across_records,
            }
        })
        .collect();
    // Field by field: the merge's type is not `Clone`, and this is the one place that needs a
    // copy of it.
    Some(CohortObservation {
        region: observation.region,
        alleles: observation.alleles.clone(),
        per_sample,
        kind: observation.kind.clone(),
        window_coverage: observation.window_coverage.clone(),
        allele_flank_deletions: observation.allele_flank_deletions.clone(),
    })
}

/// **Add one sample's deciding partial reads to its support rows**, each to the row of the
/// allele it singled out and its own read group, a row made where the sample had none — and
/// leave the rows in the merge's ascending `(allele, read group)` order, which every walk over
/// them asserts.
fn add_deciding_reads(
    supported: &mut Vec<SupportedAllele>,
    deciding: &[(usize, &PartialObservation)],
) {
    if deciding.is_empty() {
        return;
    }
    for &(allele, partial) in deciding {
        let read_group = partial.read_group;
        match supported
            .iter_mut()
            .find(|row| row.allele == allele && row.read_group == read_group)
        {
            Some(row) => {
                row.support.num_reads = row.support.num_reads.saturating_add(partial.num_reads);
                row.support.q_sum += partial.q_sum;
            }
            // Only the read count and the error mass: selection reads nothing else of a row,
            // and a partial read carries nothing else (`PartialObservation::q_sum`).
            None => supported.push(SupportedAllele {
                allele,
                read_group,
                support: AlleleSupport {
                    num_reads: partial.num_reads,
                    q_sum: partial.q_sum,
                    ..AlleleSupport::default()
                },
            }),
        }
    }
    supported.sort_by_key(|row| (row.allele, row.read_group));
}

/// **The one allele of the table a partial read's bases fit, if exactly one does.**
fn the_one_allele_a_partial_read_fits(
    observation: &CohortObservation,
    partial: &PartialObservation,
    flanks: (&[u8], &[u8]),
    locus_len: LocusLen,
) -> Option<usize> {
    // Reads the calling run re-read after the merge carry their own comparison rule. None reach
    // an STR locus today — the re-read is the generic path's — and if one did it would count
    // for nothing rather than be compared by the wrong rule.
    if partial.sequence_may_run_on_past.is_some() {
        return None;
    }
    let witness = ReadWitness::Partial {
        positions: partial.witnessed_in_locus.clone(),
    };
    let held_flank = HeldFlank::of_witness(&witness, locus_len, partial.held_flank)?;
    let mut fitting = (0..observation.alleles.len()).filter(|&allele| {
        partial_read_fits(
            &partial.bases,
            held_flank,
            &observation.alleles[allele],
            observation.flank_deletion_of(allele),
            flanks,
        )
    });
    let only_fitting_allele = fitting.next()?;
    fitting.next().is_none().then_some(only_fitting_allele)
}

/// **Whether a partial read's bases `read`, from flank `held_flank`, are what a carrier of
/// `allele` would show** — the allele followed by the far reference flank (less the bases the
/// allele's own edge deletion takes on that side), starting with them from the left, or ending
/// with them from the right. Exact: a read carrying a sequencing error fits nothing and so counts
/// for nothing. A read showing more flank bases than the stored flank holds fits nothing either;
/// it would have anchored on them and been a crossing read.
fn partial_read_fits(
    read: &[u8],
    held_flank: HeldFlank,
    allele: &[u8],
    flank_deletion: FlankDeletion,
    (left_flank, right_flank): (&[u8], &[u8]),
) -> bool {
    let bases_in_allele = read.len().min(allele.len());
    match held_flank {
        HeldFlank::Left => {
            let far_flank =
                &right_flank[usize::from(flank_deletion.right).min(right_flank.len())..];
            read[..bases_in_allele] == allele[..bases_in_allele]
                && far_flank.starts_with(&read[bases_in_allele..])
        }
        HeldFlank::Right => {
            let far_flank = &left_flank
                [..left_flank.len() - usize::from(flank_deletion.left).min(left_flank.len())];
            let bases_in_flank = read.len() - bases_in_allele;
            read[bases_in_flank..] == allele[allele.len() - bases_in_allele..]
                && far_flank.ends_with(&read[..bases_in_flank])
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::*;
    use super::*;
    use crate::calling::allele_candidates::SelectionScratch;
    use crate::calling::allele_candidates::ssr::{
        SsrSelectionConfig, select_ssr, select_ssr_bundle,
    };
    use crate::locus_generation::{LocusKind, SsrBundleDetail, SsrDetail, WitnessedLocusPositions};
    use crate::types::{Motif, Ploidy, Position, ReadGroupId};

    const LEFT_FLANK: &[u8] = b"GGTCAGTCAG";
    const RIGHT_FLANK: &[u8] = b"TTCAGTCAGG";

    fn ca(copies: usize) -> Vec<u8> {
        b"CA".repeat(copies)
    }

    /// A locus whose region is as long as its reference allele, which is what the merge builds
    /// and what the witnesses are laid against.
    fn locus_spanning_its_reference(
        alleles: &[&[u8]],
        per_sample: Vec<crate::run::cohort_merge::build::SampleSupport>,
    ) -> CohortObservation {
        let mut observation = locus_of(alleles, per_sample);
        observation.region.end =
            Position(observation.region.start.get() + alleles[0].len() as u64 - 1);
        observation
    }

    fn segment(
        alleles: &[&[u8]],
        per_sample: Vec<crate::run::cohort_merge::build::SampleSupport>,
    ) -> CohortObservation {
        let mut observation = locus_spanning_its_reference(alleles, per_sample);
        observation.kind = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"CA").unwrap(),
            left_flank: Box::from(LEFT_FLANK),
            right_flank: Box::from(RIGHT_FLANK),
        });
        observation
    }

    /// `num_reads` reads that held `held` and showed `bases`, on a locus of `locus_len` positions,
    /// with an error mass of `-1` per row. A read reaching the far border has a whole-locus
    /// witness and its side recorded, as the mint records it; a shorter one has its side in its
    /// witness alone.
    fn partial_read(
        bases: &[u8],
        held: HeldFlank,
        locus_len: u16,
        num_reads: u32,
    ) -> PartialObservation {
        let reach = (bases.len() as u16).min(locus_len);
        let witnessed = match held {
            HeldFlank::Left => WitnessedLocusPositions::one_run_from_offset_and_length(0, reach),
            HeldFlank::Right => {
                WitnessedLocusPositions::one_run_from_offset_and_length(locus_len - reach, reach)
            }
        }
        .expect("a witness of at least one position");
        PartialObservation {
            witnessed_in_locus: witnessed,
            read_group: ReadGroupId(0),
            bases: Box::from(bases),
            num_reads,
            q_sum: -1.0,
            sequence_may_run_on_past: None,
            held_flank: (usize::from(reach) == usize::from(locus_len)).then_some(held),
        }
    }

    fn config() -> SsrSelectionConfig {
        SsrSelectionConfig::at_ploidy(Ploidy::try_new(2).unwrap())
    }

    fn admitted(selection: &crate::calling::allele_candidates::LocusSelection) -> Vec<Vec<u8>> {
        selection.alleles().iter().map(<[u8]>::to_vec).collect()
    }

    /// A heterozygote at 10 and 12 `CA` copies, the 12 seen by one crossing read: without its
    /// partial read the floor of two refuses it. A left-held read that ran past the 12's end into
    /// the flank fits the 12 alone, so it is the second read, and the 12 is a candidate.
    #[test]
    fn a_partial_read_reaching_the_flank_lifts_a_once_seen_allele_over_the_floor() {
        let reference = ca(10);
        let longer = ca(12);
        let reaching = [&longer[..], &RIGHT_FLANK[..3]].concat();
        let sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
        let without = segment(&[&reference, &longer], vec![sample.clone()]);
        let mut sample = sample;
        sample.partials = vec![partial_read(&reaching, HeldFlank::Left, 20, 1)];
        let observation = segment(&[&reference, &longer], vec![sample]);

        let refused = select_ssr(&without, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&refused.selection), vec![reference.clone()]);
        let lifted = select_ssr(&observation, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&lifted.selection), vec![reference, longer]);
    }

    /// **The same for a contraction**, whose partial reads are shorter than the locus, so their
    /// side is in the witness alone: a left-held read of 8 copies and one flank base, and a
    /// right-held one of the flank's last base and 8 copies, each lift the 8.
    #[test]
    fn a_partial_read_reaching_the_flank_lifts_a_once_seen_contraction() {
        let reference = ca(10);
        let shorter = ca(8);
        for (bases, held) in [
            ([&shorter[..], &RIGHT_FLANK[..1]].concat(), HeldFlank::Left),
            (
                [&LEFT_FLANK[LEFT_FLANK.len() - 1..], &shorter[..]].concat(),
                HeldFlank::Right,
            ),
        ] {
            let read = partial_read(&bases, held, 20, 1);
            assert_eq!(read.held_flank, None, "the witness alone must say the side");
            let mut sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
            sample.partials = vec![read];
            let observation = segment(&[&reference, &shorter], vec![sample]);
            let selected = select_ssr(&observation, &config(), &mut SelectionScratch::new());
            assert_eq!(
                admitted(&selected.selection),
                vec![reference.clone(), shorter.clone()],
                "held {held:?}"
            );
        }
    }

    /// At a locus where some allele is a deletion across either edge, no partial read counts.
    #[test]
    fn no_partial_read_counts_where_an_allele_deletes_across_the_edge() {
        let reference = ca(10);
        let longer = ca(12);
        let shorter = ca(9);
        let reaching = [&longer[..], &RIGHT_FLANK[..3]].concat();
        for deletion in [
            FlankDeletion { left: 2, right: 0 },
            FlankDeletion { left: 0, right: 2 },
        ] {
            let mut sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
            sample.partials = vec![partial_read(&reaching, HeldFlank::Left, 20, 1)];
            let mut observation = segment(&[&reference, &longer, &shorter], vec![sample]);
            observation.allele_flank_deletions =
                vec![FlankDeletion::NONE, FlankDeletion::NONE, deletion];
            assert!(
                count_deciding_partial_reads(&observation, (LEFT_FLANK, RIGHT_FLANK)).is_none(),
                "{deletion:?}"
            );
        }
    }

    /// The same from the right flank: the read shows the flank's last bases, then the allele's
    /// end.
    #[test]
    fn a_right_held_partial_read_counts_from_the_end() {
        let reference = ca(10);
        let longer = ca(12);
        let reaching = [&LEFT_FLANK[LEFT_FLANK.len() - 2..], &longer[..]].concat();
        let mut sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
        sample.partials = vec![partial_read(&reaching, HeldFlank::Right, 20, 1)];
        let observation = segment(&[&reference, &longer], vec![sample]);
        let selected = select_ssr(&observation, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&selected.selection), vec![reference, longer]);
    }

    /// A partial read that ended inside the repeat short of every allele's end fits every allele,
    /// so it decides nothing; one carrying a wrong base fits none. Neither counts.
    #[test]
    fn a_partial_read_fitting_two_alleles_or_none_counts_nowhere() {
        let reference = ca(10);
        let longer = ca(12);
        let inside = ca(9);
        let mut misread = [&longer[..], &RIGHT_FLANK[..3]].concat();
        misread[4] = b'G';
        for bases in [inside, misread] {
            let mut sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
            sample.partials = vec![partial_read(&bases, HeldFlank::Left, 20, 5)];
            let observation = segment(&[&reference, &longer], vec![sample]);
            assert!(
                count_deciding_partial_reads(&observation, (LEFT_FLANK, RIGHT_FLANK)).is_none()
            );
            let selected = select_ssr(&observation, &config(), &mut SelectionScratch::new());
            assert_eq!(admitted(&selected.selection), vec![reference.clone()]);
        }
    }

    /// **A read that ended inside the repeat, past the end of every allele but the longest,
    /// counts for the longest** — it fits that one alone, and the rule reads only the table
    /// (spec §8.1). Pinned so that requiring a flank base is a decision someone takes, not a
    /// drift.
    #[test]
    fn a_read_ending_inside_the_repeat_past_every_shorter_allele_counts_for_the_longest() {
        let reference = ca(10);
        let longer = ca(12);
        let inside = ca(11);
        let mut sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
        sample.partials = vec![partial_read(&inside, HeldFlank::Left, 20, 1)];
        let observation = segment(&[&reference, &longer], vec![sample]);
        let counted =
            count_deciding_partial_reads(&observation, (LEFT_FLANK, RIGHT_FLANK)).unwrap();
        assert_eq!(counted.per_sample[0].pooled_support_for(1).num_reads, 2);
    }

    /// **Deciding partials count in the denominator too**, so a stutter sequence's share is not
    /// raised by them. 30 crossing reads at 10 copies and 2 at 11 (a tenth of 32 is 4: refused),
    /// plus 30 deciding partials of the 10 and 2 of the 11: 4 of 64 against a bar of 7, still
    /// refused. The partial rows' error masses join the crossing rows' (`-1` a row here).
    #[test]
    fn deciding_partials_leave_a_stutter_sequences_share_where_it_was() {
        let reference = ca(10);
        let stutter = ca(11);
        let of_reference = [&reference[..], &RIGHT_FLANK[..3]].concat();
        let of_stutter = [&stutter[..], &RIGHT_FLANK[..3]].concat();
        let mut sample = sample_showing(0, vec![row(0, 30, -30.0), row(1, 2, -2.0)]);
        let mut stutter_partial = partial_read(&of_stutter, HeldFlank::Left, 20, 2);
        stutter_partial.q_sum = -0.5;
        sample.partials = vec![
            partial_read(&of_reference, HeldFlank::Left, 20, 30),
            stutter_partial,
        ];
        let observation = segment(&[&reference, &stutter], vec![sample]);
        let counted =
            count_deciding_partial_reads(&observation, (LEFT_FLANK, RIGHT_FLANK)).unwrap();
        let counted_sample = &counted.per_sample[0];
        assert_eq!(counted_sample.pooled_support_for(0).num_reads, 60);
        assert_eq!(counted_sample.pooled_support_for(1).num_reads, 4);
        assert_eq!(counted_sample.pooled_support_for(0).q_sum, -31.0);
        assert_eq!(counted_sample.pooled_support_for(1).q_sum, -2.5);
        assert!(
            counted_sample.partials.is_empty(),
            "the counted reads are in the rows, not carried a second time"
        );
        let selected = select_ssr(&observation, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&selected.selection), vec![reference]);
    }

    /// **The leftover is the merge's, not the copy's.** The stutter sequence above is refused,
    /// and its 2 crossing reads are what the record's `DP` adds as "no written allele explains
    /// these" — the record adds the 32 partial reads separately, so the 2 that counted for the
    /// stutter must not be in the leftover as well (they were, and `DP` came out 66 of 64 reads).
    #[test]
    fn the_leftover_holds_the_merges_reads_and_not_the_counted_partials() {
        let reference = ca(10);
        let stutter = ca(11);
        let of_reference = [&reference[..], &RIGHT_FLANK[..3]].concat();
        let of_stutter = [&stutter[..], &RIGHT_FLANK[..3]].concat();
        let mut sample = sample_showing(0, vec![row(0, 30, -30.0), row(1, 2, -2.0)]);
        sample.partials = vec![
            partial_read(&of_reference, HeldFlank::Left, 20, 30),
            partial_read(&of_stutter, HeldFlank::Left, 20, 2),
        ];
        let observation = segment(&[&reference, &stutter], vec![sample]);
        let selected = select_ssr(&observation, &config(), &mut SelectionScratch::new());
        let leftover = selected.selection.unmatched()[0];
        assert_eq!(leftover.num_reads, 2);
        assert_eq!(leftover.q_sum, -2.0);
        assert_eq!(leftover.earned_reads_cut_by_the_cap, 0);
    }

    /// Counted reads land on the sample whose reads they are: two samples, a deciding read in the
    /// second only.
    #[test]
    fn deciding_reads_count_only_for_the_sample_whose_reads_they_are() {
        let reference = ca(10);
        let longer = ca(12);
        let reaching = [&longer[..], &RIGHT_FLANK[..1]].concat();
        let first = sample_showing(0, vec![row(0, 5, -5.0)]);
        let mut second = sample_showing(1, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
        second.partials = vec![partial_read(&reaching, HeldFlank::Left, 20, 1)];
        let observation = segment(&[&reference, &longer], vec![first, second]);
        let counted =
            count_deciding_partial_reads(&observation, (LEFT_FLANK, RIGHT_FLANK)).unwrap();
        assert_eq!(counted.per_sample[0].pooled_support_for(1).num_reads, 0);
        assert_eq!(counted.per_sample[0].pooled_support_for(0).num_reads, 5);
        assert_eq!(counted.per_sample[1].pooled_support_for(1).num_reads, 2);
    }

    /// **A sample with no crossing read at all can nominate on two deciding partial reads**, if
    /// the sequence they single out is in the table — spelled by a crossing read of another
    /// sample. Without them the second sample shows nothing and nominates nothing.
    #[test]
    fn a_sample_without_crossing_reads_nominates_on_two_deciding_partial_reads() {
        let reference = ca(10);
        let longer = ca(12);
        let reaching = [&longer[..], &RIGHT_FLANK[..2]].concat();
        let first = sample_showing(0, vec![row(0, 6, -6.0), row(1, 1, -1.0)]);
        let mut second = sample_showing(1, Vec::new());
        second.partials = vec![partial_read(&reaching, HeldFlank::Left, 20, 2)];
        let without = segment(
            &[&reference, &longer],
            vec![first.clone(), sample_showing(1, Vec::new())],
        );
        let observation = segment(&[&reference, &longer], vec![first, second]);
        let refused = select_ssr(&without, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&refused.selection), vec![reference.clone()]);
        let selected = select_ssr(&observation, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&selected.selection), vec![reference, longer]);
    }

    /// At a bundle the read is compared with the hull's allele and the far flank, as at a
    /// segment: a read reaching past a longer `TTG` run into the flank lifts that allele.
    #[test]
    fn a_bundle_partial_read_lifts_a_once_seen_allele() {
        use crate::region_typing::reference_runs::ReferenceRun;
        let hull = |ca_copies: usize, ttg: usize| {
            [&b"GG"[..], &ca(ca_copies), b"GT", &b"TTG".repeat(ttg), b"C"].concat()
        };
        let reference = hull(8, 5);
        let longer = hull(8, 7);
        let detail = SsrBundleDetail {
            left_flank: Box::from(LEFT_FLANK),
            right_flank: Box::from(RIGHT_FLANK),
            runs: vec![
                ReferenceRun::new(2, 16, Motif::new(b"CA").unwrap()).unwrap(),
                ReferenceRun::new(20, 15, Motif::new(b"TTG").unwrap()).unwrap(),
            ]
            .into_boxed_slice(),
        };
        let reaching = [&longer[..], &RIGHT_FLANK[..2]].concat();
        let sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0)]);
        let mut without =
            locus_spanning_its_reference(&[&reference, &longer], vec![sample.clone()]);
        without.kind = LocusKind::SsrBundle(detail.clone());
        let mut sample = sample;
        sample.partials = vec![partial_read(
            &reaching,
            HeldFlank::Left,
            reference.len() as u16,
            1,
        )];
        let mut observation = locus_spanning_its_reference(&[&reference, &longer], vec![sample]);
        observation.kind = LocusKind::SsrBundle(detail);
        let refused = select_ssr_bundle(&without, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&refused.selection), vec![reference.clone()]);
        let lifted = select_ssr_bundle(&observation, &config(), &mut SelectionScratch::new());
        assert_eq!(admitted(&lifted.selection), vec![reference, longer]);
    }

    /// The merge's row order holds after a counted read creates a row: a read for allele 1 from a
    /// second read group lands between allele 1's first-group row and allele 2's.
    #[test]
    fn counted_rows_keep_the_merges_order() {
        let reference = ca(10);
        let longer = ca(12);
        let longest = ca(14);
        let reaching = [&longer[..], &RIGHT_FLANK[..3]].concat();
        let mut sample = sample_showing(0, vec![row(0, 3, -3.0), row(1, 1, -1.0), row(2, 1, -1.0)]);
        let mut deciding = partial_read(&reaching, HeldFlank::Left, 20, 1);
        deciding.read_group = ReadGroupId(1);
        sample.partials = vec![deciding];
        let observation = segment(&[&reference, &longer, &longest], vec![sample]);
        // 12 copies then flank fits the 12 alone: the 14 holds `CA` where the flank's `TTC` is.
        let counted =
            count_deciding_partial_reads(&observation, (LEFT_FLANK, RIGHT_FLANK)).unwrap();
        let keys: Vec<(usize, u32)> = counted.per_sample[0]
            .supported
            .iter()
            .map(|row| (row.allele, row.read_group.0))
            .collect();
        assert_eq!(keys, vec![(0, 0), (1, 0), (1, 1), (2, 0)]);
    }

    /// **The far flank loses the bases the allele's own edge deletion took**, on the side the
    /// read reads into: two bases off the right flank for a left-held read, two off the end of
    /// the left flank for a right-held one — and not on the other side.
    #[test]
    fn partial_read_fits_takes_the_alleles_own_flank_deletion_from_the_far_side() {
        let flanks = (LEFT_FLANK, RIGHT_FLANK);
        let allele = ca(10);
        let past_a_right_deletion = [&allele[..], &RIGHT_FLANK[2..4]].concat();
        let right_two = FlankDeletion { left: 0, right: 2 };
        let left_two = FlankDeletion { left: 2, right: 0 };
        assert!(partial_read_fits(
            &past_a_right_deletion,
            HeldFlank::Left,
            &allele,
            right_two,
            flanks
        ));
        assert!(!partial_read_fits(
            &past_a_right_deletion,
            HeldFlank::Left,
            &allele,
            FlankDeletion::NONE,
            flanks
        ));
        assert!(!partial_read_fits(
            &past_a_right_deletion,
            HeldFlank::Left,
            &allele,
            left_two,
            flanks
        ));
        let end = LEFT_FLANK.len();
        let past_a_left_deletion = [&LEFT_FLANK[end - 4..end - 2], &allele[..]].concat();
        assert!(partial_read_fits(
            &past_a_left_deletion,
            HeldFlank::Right,
            &allele,
            left_two,
            flanks
        ));
        assert!(!partial_read_fits(
            &past_a_left_deletion,
            HeldFlank::Right,
            &allele,
            FlankDeletion::NONE,
            flanks
        ));
        assert!(!partial_read_fits(
            &past_a_left_deletion,
            HeldFlank::Right,
            &allele,
            right_two,
            flanks
        ));
    }

    /// The edges of the comparison: an empty read fits anything; a read exactly the allele fits
    /// it; one showing the whole stored flank fits, one base more does not (it would have
    /// anchored), and a flank clamped short by a contig end is no panic.
    #[test]
    fn partial_read_fits_at_its_boundaries() {
        let flanks = (LEFT_FLANK, RIGHT_FLANK);
        let allele = ca(3);
        let none = FlankDeletion::NONE;
        assert!(partial_read_fits(
            b"",
            HeldFlank::Left,
            &allele,
            none,
            flanks
        ));
        assert!(partial_read_fits(
            &allele,
            HeldFlank::Right,
            &allele,
            none,
            flanks
        ));
        let whole_flank = [&allele[..], RIGHT_FLANK].concat();
        assert!(partial_read_fits(
            &whole_flank,
            HeldFlank::Left,
            &allele,
            none,
            flanks
        ));
        let past_the_flank = [&whole_flank[..], b"A"].concat();
        assert!(!partial_read_fits(
            &past_the_flank,
            HeldFlank::Left,
            &allele,
            none,
            flanks
        ));
        let clamped = (&LEFT_FLANK[..2], &RIGHT_FLANK[..2]);
        let three_flank_bases = [&allele[..], &RIGHT_FLANK[..3]].concat();
        assert!(!partial_read_fits(
            &three_flank_bases,
            HeldFlank::Left,
            &allele,
            none,
            clamped
        ));
        assert!(!partial_read_fits(
            &[&LEFT_FLANK[..3], &allele[..]].concat(),
            HeldFlank::Right,
            &allele,
            FlankDeletion { left: 9, right: 0 },
            clamped
        ));
    }
}
