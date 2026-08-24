//! **The ordinary SNP/indel path: narrowing one locus's allele table to the sequences worth
//! calling over** (`doc/devel/ng/spec/candidate_alleles.md` §3, §4; arch §3.1).
//!
//! One function, [`select_generic`]. The repeat-tract path is a sibling that takes different
//! evidence and returns different extras, and which of the two runs is decided by the locus's
//! kind rather than by a swappable recipe — see this module's parent for why that is two
//! functions and not two implementations of one trait.

use super::{
    AlleleRemap, CandidateSelectionConfig, LocusSelection, SelectionScratch, SelectionVerdict,
    UnmatchedSupport, summarise_alleles,
};
use crate::ng::calling::CandidateAlleles;
use crate::ng::locus_generation::LocusKind;
use crate::ng::run::cohort_merge::build::CohortObservation;
use crate::ng::types::AlleleId;

/// **Narrow one locus's allele table to the sequences worth calling over** (spec §3, §4).
///
/// **Two walks of the locus's rows, and neither allocates beyond the answer.** The first is
/// [`summarise_alleles`], which asks each sample separately whether its own reads lent each
/// sequence enough — `max(2 reads, ceil(share × that sample's compared reads))` — and one
/// sample reaching it admits the sequence for the whole cohort (spec §3.2). Then the survivors
/// are admitted, in the merge table's own order, and where each of the merge's alleles ended
/// up is recorded as they go.
///
/// **The order the steps have to run in is fixed, and arch §3.1's sentence about it cannot be
/// implemented as written.** That sentence says one pass "admits the survivors in table order,
/// applies the cap, and fills the leftover"; the cap has to run *before* admission, because
/// admission needs to know which alleles survived it, and the leftover has to run *after*,
/// because it is per sample where admission is per allele and it reads the finished remapping.
/// So the shape is: fold the rows, cap the survivor list (step C2), admit, then walk the rows
/// again for the leftover (step C3).
///
/// **This step is the middle two of those four and returns placeholders for the others.** The
/// verdict is always [`SelectionVerdict::Selected`], and the leftover is one zeroed entry per
/// covering sample — the right length in the right order, which is what C3 fills.
///
/// **The reference is admitted first and is asked nothing.** It is exempt from the rule and
/// from the cap (spec §6.1), so it is seeded structurally rather than by reading whether it
/// passed — which it may well have, since the fold asks it like any other allele.
///
/// **A locus can select down to the reference alone, and that is a first-class outcome, not
/// an error.** The merge builds a locus when some sample's non-reference reads *pooled* reach
/// its rule, and two reads split one and one across two alternatives clear that while clearing
/// neither allele's own bar. Measured at more than one built locus in four on both benchmarks —
/// 27.4% on the 63-accession tomato panel, 27.3% on the GIAB trio at 30× and 28.0% at 300×
/// (spec §6.2). What the run does with such a locus is emission's business.
///
/// **It relies on the merge's allele table holding each sequence once**, which
/// [`CohortObservation::alleles`] states and `AlleleTable` enforces by interning on the bytes
/// themselves. Two copies of one sequence would be admitted as two candidates here, and the
/// read likelihood would then split one allele's evidence across both — a genotype that looks
/// ordinary. Held by a `debug_assert!` rather than in release, because the check is a scan of
/// the table where the invariant belongs to the producer.
///
/// # Panics
///
/// On a locus whose allele table is empty. The merge always interns the reference at index 0
/// ([`CohortObservation::alleles`]), so this is a caller bug rather than an input. Without the
/// assertion the failure is `index out of bounds` two statements later, where the reference's
/// bases are read — measured, not assumed: a reviewer deleted the assertion and got exactly
/// that, which points at the wrong line and says nothing about the locus.
///
/// **And, until the cap of step C2 lands, on a locus where more than 65,535 alternatives clear
/// the rule**: this step admits every one of them, and [`CandidateAlleles::admit`] refuses an
/// admission no [`AlleleId`] could name. Nothing upstream bounds a locus's allele table at that
/// width — see [`SelectionVerdict::Truncated`]'s own documentation, which records a review
/// building a locus of 70,001 alternatives. C2 makes it unreachable and this paragraph goes
/// with it.
pub fn select_generic(
    observation: &CohortObservation,
    config: &CandidateSelectionConfig,
    scratch: &mut SelectionScratch,
) -> LocusSelection {
    assert!(
        !observation.alleles.is_empty(),
        "a cohort locus always holds at least its reference allele, and the one at {} holds none",
        observation.region
    );
    debug_assert!(
        observation
            .alleles
            .iter()
            .collect::<std::collections::HashSet<_>>()
            .len()
            == observation.alleles.len(),
        "the merge interns each sequence once, and two copies of one would be admitted as two \
         candidates at {} — one allele's evidence split across both",
        observation.region
    );
    summarise_alleles(observation, config.min_allele_support, scratch);

    let SelectionScratch {
        per_allele,
        ranked_table_indices,
    } = scratch;
    // Every alternative some sample's reads earned, in the merge table's own order — which is
    // also the order they are admitted in, so nothing has to sort it. Step C2's cap is what
    // first reorders this buffer, and it puts it back before admission.
    ranked_table_indices.extend(
        (1..observation.alleles.len())
            .filter(|&index| per_allele[index].cleared_the_bar())
            .map(|index| {
                u32::try_from(index).expect("a merge table narrower than four billion alleles")
            }),
    );

    let mut alleles = CandidateAlleles::new(observation.alleles[0].clone(), LocusKind::Generic);
    let mut remap = AlleleRemap::with_all_dropped(observation.alleles.len());
    remap.admit(0, AlleleId::REFERENCE);
    for &table_index in ranked_table_indices.iter() {
        let table_index = table_index as usize;
        let candidate = alleles.admit(observation.alleles[table_index].clone());
        remap.admit(table_index, candidate);
    }

    let covering_samples = observation.per_sample.len();
    LocusSelection::new(
        alleles,
        SelectionVerdict::Selected,
        vec![UnmatchedSupport::default(); covering_samples],
        remap,
        covering_samples,
    )
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::*;
    use super::*;
    use crate::ng::run::cohort_merge::MinAltReads;
    use crate::ng::run::cohort_merge::build::SampleSupport;
    use crate::ng::types::ReadGroupId;

    /// Narrow `observation` under a support rule of `floor` reads or `share`, and the default
    /// cap.
    fn selection_of(
        observation: &CohortObservation,
        min_allele_support: MinAltReads,
    ) -> LocusSelection {
        let config = CandidateSelectionConfig {
            min_allele_support,
            ..CandidateSelectionConfig::DEFAULT
        };
        select_generic(observation, &config, &mut SelectionScratch::new())
    }

    /// The surviving sequences, in candidate-id order.
    fn surviving_bases(selection: &LocusSelection) -> Vec<Vec<u8>> {
        selection.alleles().iter().map(<[u8]>::to_vec).collect()
    }

    /// **The plan's oracle for this step: the round trip, with a hole in the middle.**
    ///
    /// Five sequences, and the middle one is the only alternative no sample earned — one
    /// sample showed it a single read where the rule asks two. So the merge's indices 1, 3 and
    /// 4 survive and 2 does not, and the survivors take the dense candidate ids 1, 2 and 3.
    ///
    /// **The hole is what the test is for.** An off-by-one in the remapping hands the calling
    /// loop a real but *wrong* allele's evidence, which is a wrong genotype rather than a
    /// panic, and a table with no gap in it cannot tell a correct remapping from one that
    /// simply counts. Both directions are asserted: every surviving merge index answers with
    /// its dense id, and the dropped one answers `None`.
    #[test]
    fn every_surviving_merge_index_answers_with_its_dense_id_and_the_dropped_one_with_none() {
        let observation = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 4, -4.0),
                    row(1, 3, -3.0),
                    row(2, 1, -1.0),
                    row(3, 2, -2.0),
                    row(4, 2, -2.0),
                ],
            )],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));

        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"C".to_vec(), b"T".to_vec(), b"AA".to_vec()],
            "the reference and the three alternatives some sample's reads earned"
        );
        let remap = selection.remap();
        assert_eq!(remap.candidate_for(0), Some(AlleleId::REFERENCE));
        assert_eq!(remap.candidate_for(1), Some(AlleleId(1)));
        assert_eq!(
            remap.candidate_for(2),
            None,
            "the one allele the rule dropped, and the hole the ids close over"
        );
        assert_eq!(remap.candidate_for(3), Some(AlleleId(2)));
        assert_eq!(remap.candidate_for(4), Some(AlleleId(3)));
        assert_eq!(remap.table_len(), 5);
        assert_eq!(remap.num_admitted(), 4);
    }

    /// **The hand-off arch §3.2 fixes, reproduced exactly.** The calling loop builds each
    /// sample's evidence by keeping the rows whose allele the remapping still knows and
    /// re-keying them onto the candidate ids, so this walks that recipe and checks the rows it
    /// yields — which is the only thing that shows the remapping is usable rather than merely
    /// self-consistent.
    ///
    /// The sample here shows four of the five sequences, one of them from two read groups, and
    /// one of them dropped. **The two read-group rows must both survive and stay apart**: they
    /// are pooled to ask the support rule and never pooled for a likelihood, because two lanes
    /// have different error rates (`doc/devel/ng/spec/read_likelihoods.md` §2.3).
    #[test]
    fn the_evidence_hand_off_keeps_every_row_of_a_surviving_allele_and_re_keys_it() {
        let observation = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 4, -4.0),
                    row_from_group(1, ReadGroupId(0), 2, -2.0),
                    row_from_group(1, ReadGroupId(1), 1, -1.0),
                    row(2, 1, -1.0),
                    row(3, 2, -2.0),
                ],
            )],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        let remap = selection.remap();

        let re_keyed: Vec<(AlleleId, ReadGroupId, u32)> = observation.per_sample[0]
            .supported
            .iter()
            .filter_map(|supported| {
                remap
                    .candidate_for(supported.allele)
                    .map(|id| (id, supported.read_group, supported.support.num_reads))
            })
            .collect();
        assert_eq!(
            re_keyed,
            vec![
                (AlleleId::REFERENCE, ReadGroupId(0), 4),
                (AlleleId(1), ReadGroupId(0), 2),
                (AlleleId(1), ReadGroupId(1), 1),
                (AlleleId(2), ReadGroupId(0), 2),
            ],
            "the dropped allele's row is gone, both of allele 1's lanes are kept and stay \
             separate, and allele 3 of the merge is candidate 2"
        );
    }

    /// **A locus can select down to the reference alone, and it is `Selected`, not an error.**
    ///
    /// The merge built this locus because the sample's non-reference reads *pooled* reach its
    /// rule — two of them — and they are split one and one across two alternatives, so neither
    /// alternative clears a bar of two on its own. Measured at more than one built locus in
    /// four on both benchmarks (spec §6.2), so it is the ordinary case and not a corner.
    #[test]
    fn a_locus_whose_alternatives_all_failed_the_rule_keeps_the_reference_and_is_selected() {
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            vec![sample_showing(
                0,
                vec![row(0, 8, -8.0), row(1, 1, -1.0), row(2, 1, -1.0)],
            )],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));

        assert_eq!(surviving_bases(&selection), vec![b"A".to_vec()]);
        assert_eq!(selection.verdict(), SelectionVerdict::Selected);
        assert_eq!(
            selection.alternative_allele_count(),
            0,
            "the number the genotype prior divides its concentration by"
        );
        assert_eq!(selection.remap().candidate_for(1), None);
        assert_eq!(selection.remap().candidate_for(2), None);
    }

    /// **The survivors are admitted in the merge table's order, not in rank order** (arch
    /// §3.1), which is the order that reaches the VCF's `ALT` column. The ranking decides
    /// *which* alleles survive; it never decides what order they appear in.
    ///
    /// Here the second alternative is much the better evidenced — 6 of one sample's 10 reads
    /// against 2 — so a pass that admitted in rank order would swap them.
    #[test]
    fn the_survivors_are_admitted_in_the_merge_tables_order_and_not_by_rank() {
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            vec![sample_showing(
                0,
                vec![row(0, 2, -2.0), row(1, 2, -2.0), row(2, 6, -6.0)],
            )],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"C".to_vec(), b"G".to_vec()]
        );
        assert_eq!(selection.remap().candidate_for(1), Some(AlleleId(1)));
        assert_eq!(selection.remap().candidate_for(2), Some(AlleleId(2)));
    }

    /// **Adding a sample that shows only reference reads changes nothing** — spec §3.2's
    /// principle at the level of the whole answer rather than of the fold, and the plan's
    /// standing property.
    ///
    /// **It covers a cohort term in the *denominator* and not in the numerator**, and an
    /// earlier version of this comment claimed both. The added sample lends the alternative no
    /// reads, so it moves any rule counting *samples* — a majority-of-covering-samples rule
    /// fails here — and moves nothing that pools the alternative's reads across the cohort.
    /// `an_alternative_no_single_sample_earned_is_refused_however_many_samples_showed_it` is
    /// the one that covers the numerator, and the two are needed together.
    #[test]
    fn a_sample_showing_only_reference_reads_changes_neither_the_list_nor_the_remapping() {
        let carrier = || sample_showing(0, vec![row(0, 1, -1.0), row(1, 2, -2.0)]);
        let alone = selection_of(
            &locus_of(&[b"A", b"C"], vec![carrier()]),
            support_rule_of(2, 0.5),
        );
        let with_a_bystander = selection_of(
            &locus_of(
                &[b"A", b"C"],
                vec![carrier(), sample_showing(1, vec![row(0, 60, -60.0)])],
            ),
            support_rule_of(2, 0.5),
        );
        assert_eq!(alone.alleles(), with_a_bystander.alleles());
        assert_eq!(alone.remap(), with_a_bystander.remap());
    }

    /// The surviving list is always a subset of the merge's table and always holds the
    /// reference first — the plan's other standing property, asserted here on a locus where
    /// some alternatives survive and some do not.
    #[test]
    fn the_survivors_are_a_subset_of_the_merge_table_with_the_reference_first() {
        let observation = locus_of(
            &[b"A", b"C", b"G", b"T"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 5, -5.0),
                    row(1, 3, -3.0),
                    row(2, 1, -1.0),
                    row(3, 4, -4.0),
                ],
            )],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        let table: Vec<Vec<u8>> = observation.alleles.iter().map(|a| a.to_vec()).collect();
        let survivors = surviving_bases(&selection);
        assert_eq!(survivors[0], table[0], "the reference is candidate 0");
        assert!(survivors.iter().all(|bases| table.contains(bases)));
        assert!(survivors.len() <= table.len());
        assert_eq!(
            selection.alleles().kind(),
            &LocusKind::Generic,
            "this path builds generic loci and nothing else"
        );
    }

    /// The leftover runs parallel to the locus's *covering* samples, which is not the run's
    /// sample order: a run of 63 accessions can produce a locus four of them covered, and
    /// indexing by the run's order would put four leftovers at four scattered positions with
    /// 59 zeroed rows between them — indistinguishable from a covering sample that dropped
    /// nothing.
    ///
    /// Step C3 fills these; C1 only has to get the length and the order right.
    #[test]
    fn the_leftover_has_one_entry_per_covering_sample() {
        let covering: Vec<SampleSupport> = [3_usize, 17, 40]
            .into_iter()
            .map(|sample| sample_showing(sample, vec![row(0, 2, -2.0), row(1, 2, -2.0)]))
            .collect();
        let selection = selection_of(&locus_of(&[b"A", b"C"], covering), support_rule_of(2, 0.0));
        assert_eq!(selection.unmatched().len(), 3);
        assert!(
            selection
                .unmatched()
                .iter()
                .all(|leftover| *leftover == UnmatchedSupport::default()),
            "C1 drops nothing, so every leftover is zeroed — and a non-zero \
             `earned_reads_cut_by_the_cap` on this line would emit every covering sample at \
             every locus as a missing genotype, which is the slip C3 is in a position to make"
        );
    }

    /// **A covering sample whose reads all stopped inside the locus still owns a leftover.**
    /// It has partials and no support rows at all, and the merge does build it —
    /// `per_sample` holds the samples that *covered* the span, not the ones that spanned it.
    ///
    /// The three-sample fixture above cannot show this, because all three of its samples have
    /// rows: a length that counted only the samples with rows passes it and fails here. The
    /// cost of getting it wrong is not a short vector but a **shifted** one — the leftover is
    /// parallel to `per_sample` by position, so every later sample's leftover slides onto its
    /// neighbour, and once C3 fills these that is a missing genotype for a sample that lost
    /// nothing and an invented one for the sample that did.
    #[test]
    fn the_leftover_counts_a_covering_sample_whose_reads_all_stopped_inside_the_locus() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 4, -4.0), row(1, 3, -3.0)]),
                sample_with_only_partials(1, 5),
                sample_showing(2, vec![row(0, 4, -4.0)]),
            ],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(selection.unmatched().len(), 3);
    }

    /// A locus whose table is the reference alone — no alternative was ever interned — comes
    /// back as the reference alone, with a remapping of one entry.
    #[test]
    fn a_reference_only_table_selects_to_the_reference() {
        let observation = locus_of(&[b"AC"], vec![sample_showing(0, vec![row(0, 3, -3.0)])]);
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(surviving_bases(&selection), vec![b"AC".to_vec()]);
        assert_eq!(selection.remap().table_len(), 1);
        assert_eq!(selection.remap().num_admitted(), 1);
    }

    /// A locus no sample covered still selects — to the reference alone, with an empty
    /// leftover. The merge does not build such a locus, and the assertion this replaces would
    /// have been the wrong shape: it is legal input to a pure function of one locus.
    #[test]
    fn a_locus_no_sample_covers_selects_to_the_reference_with_an_empty_leftover() {
        let selection = selection_of(
            &locus_of(&[b"A", b"C"], Vec::new()),
            support_rule_of(2, 0.0),
        );
        assert_eq!(surviving_bases(&selection), vec![b"A".to_vec()]);
        assert!(selection.unmatched().is_empty());
    }

    /// An empty allele table is a caller bug, and it is refused here rather than three lines
    /// later by [`CandidateAlleles::new`] complaining about an empty reference allele — a
    /// message that would send the reader looking in the wrong place.
    #[test]
    #[should_panic(expected = "always holds at least its reference allele")]
    fn a_locus_with_no_alleles_at_all_is_refused() {
        selection_of(&locus_of(&[], Vec::new()), support_rule_of(2, 0.0));
    }

    /// The scratch is reused across loci, so a second narrowing must carry nothing of the
    /// first — in particular the buffer of surviving indices, which is appended to rather than
    /// assigned.
    #[test]
    fn narrowing_a_second_locus_with_the_same_scratch_carries_nothing_from_the_first() {
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.0),
            ..CandidateSelectionConfig::DEFAULT
        };
        let mut scratch = SelectionScratch::new();
        let wide = locus_of(
            &[b"A", b"C", b"G"],
            vec![sample_showing(
                0,
                vec![row(0, 2, -2.0), row(1, 3, -3.0), row(2, 3, -3.0)],
            )],
        );
        let first = select_generic(&wide, &config, &mut scratch);
        assert_eq!(first.alternative_allele_count(), 2);

        let narrow = locus_of(
            &[b"A", b"T"],
            vec![sample_showing(0, vec![row(0, 2, -2.0), row(1, 3, -3.0)])],
        );
        let second = select_generic(&narrow, &config, &mut scratch);
        assert_eq!(
            surviving_bases(&second),
            vec![b"A".to_vec(), b"T".to_vec()],
            "the first locus's alleles must not reappear in the second's table"
        );
        assert_eq!(second.remap().table_len(), 2);
    }

    /// **Two samples lending one read each pool to two, and the rule refuses both** — the
    /// input that separates spec §3.2's per-sample rule from any rule that pools the
    /// alternative's reads across the cohort.
    ///
    /// **Without it the whole admission rule can be replaced by `cohort_reads >= 2` and every
    /// other test here stays green**, because no other fixture has two samples each lending an
    /// alternative *less* than the floor. A cohort term in this rule makes a sample's candidate
    /// list depend on who else is in the run — the property this module exists to keep — and it
    /// admits error alleles in proportion to cohort size: at 63 accessions or at a thousand,
    /// one error read per sample pools past any fixed floor.
    #[test]
    fn an_alternative_no_single_sample_earned_is_refused_however_many_samples_showed_it() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 1, -1.0), row(1, 1, -1.0)]),
                sample_showing(1, vec![row(0, 1, -1.0), row(1, 1, -1.0)]),
            ],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec()],
            "the two reads pool to the floor of two, and neither sample reached it alone"
        );
        assert_eq!(selection.remap().candidate_for(1), None);
    }

    /// **The share the caller configured reaches the fold whole**, shown in the only regime
    /// where it can be: above 41 compared reads, where `ceil(0.05 × 300) = 15` asks for more
    /// than a floor of 2.
    ///
    /// **Every other fixture in this file is a handful of reads, where the floor decides
    /// whatever the share is** — so without this one the share could be dropped, or the whole
    /// configured rule replaced by the shipped default, and the suite would not notice. The
    /// regime it covers is not a corner: the GIAB trio runs at 30× and 300×, so it is where
    /// the high-depth benchmark spends all of its time, and a rule that degraded to its floor
    /// there would admit sequencing error as a candidate — 10 reads in 300 is about the error
    /// rate — with a longer `ALT` list and no crash.
    #[test]
    fn the_configured_share_refuses_an_allele_the_floor_alone_would_admit() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(
                0,
                vec![row(0, 290, -290.0), row(1, 10, -10.0)],
            )],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.05));
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec()],
            "10 of 300 compared reads clears the floor of 2 and misses the share's 15"
        );
        assert_eq!(selection.remap().candidate_for(1), None);
    }

    /// The scratch reused the other way round — a **narrow** locus and then a much wider one.
    /// The sibling test above only shrinks, so a per-allele buffer that grew without being
    /// refilled would pass it; this one admits twelve alleles after two.
    #[test]
    fn narrowing_a_wider_locus_after_a_narrow_one_admits_every_survivor() {
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.0),
            ..CandidateSelectionConfig::DEFAULT
        };
        let mut scratch = SelectionScratch::new();
        let narrow = locus_of(
            &[b"A", b"C"],
            vec![sample_showing(0, vec![row(0, 4, -4.0), row(1, 3, -3.0)])],
        );
        assert_eq!(
            select_generic(&narrow, &config, &mut scratch)
                .alleles()
                .len(),
            2
        );

        let bases: Vec<Vec<u8>> = (0..12_u8).map(|i| vec![b'A' + i]).collect();
        let table: Vec<&[u8]> = bases.iter().map(Vec::as_slice).collect();
        let rows: Vec<_> = (0..12).map(|allele| row(allele, 3, -3.0)).collect();
        let wide = locus_of(&table, vec![sample_showing(0, rows)]);
        let second = select_generic(&wide, &config, &mut scratch);
        assert_eq!(second.alleles().len(), 12);
        assert_eq!(second.remap().candidate_for(11), Some(AlleleId(11)));
    }

    /// A merge table holding one sequence twice is a producer bug the read likelihood would
    /// turn into a genotype rather than a crash: both copies become candidates and one
    /// allele's evidence is split across them. Held in debug builds only, because the check
    /// is a scan of the table and the invariant belongs to the merge, which interns on the
    /// bytes themselves.
    #[test]
    #[should_panic(expected = "interns each sequence once")]
    #[cfg(debug_assertions)]
    fn a_merge_table_holding_one_sequence_twice_is_refused() {
        selection_of(
            &locus_of(
                &[b"A", b"C", b"C"],
                vec![sample_showing(
                    0,
                    vec![row(0, 3, -3.0), row(1, 3, -3.0), row(2, 3, -3.0)],
                )],
            ),
            support_rule_of(2, 0.0),
        );
    }
}
