//! **The ordinary SNP/indel path: narrowing one locus's allele table to the sequences worth
//! calling over** (`doc/devel/ng/spec/candidate_alleles.md` §3, §4; arch §3.1).
//!
//! One function, [`select_generic`]. The repeat-tract path is a sibling that takes different
//! evidence and returns different extras, and which of the two runs is decided by the locus's
//! kind rather than by a swappable recipe — see this module's parent for why that is two
//! functions and not two implementations of one trait.

use super::{
    AlleleRemap, CandidateSelectionConfig, LocusSelection, RankedAlternative, SelectionScratch,
    SelectionVerdict, compare_best_first, leftover_of, summarise_alleles,
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
/// again for the leftover — which is what [`leftover_of`] does, per sample.
///
/// **Above the cap the list is cut and the locus is still called; it is never refused**
/// (spec §4.1). What is cut is the worst-evidenced by [`compare_best_first`], and the count is
/// reported as [`SelectionVerdict::Truncated`] — not how many alleles were dropped in total,
/// since an alternative that failed the admission rule was never a candidate for the cap.
/// **Refusing the whole locus was the alternative and it loses more**: it is what HipSTR does
/// above 1,000 haplotypes and what production's repeat-tract path does above 24 candidates
/// (`src/ssr/cohort/candidate_set.rs`),
/// and at 63 accessions it costs 62 samples a locus they were called at perfectly well because
/// one accession carried something rare.
///
/// **The reference is admitted first and is asked nothing.** It is exempt from the rule and
/// from the cap (spec §6.1), so it is seeded structurally rather than by reading whether it
/// passed — which it may well have, since the fold asks it like any other allele.
///
/// **A locus can select down to the reference alone, and that is a first-class outcome, not
/// an error.** The merge builds a locus when some sample's non-reference reads *pooled* reach
/// its rule, and two reads split one and one across two alternatives clear that while clearing
/// neither allele's own bar. Measured at more than one built locus in four on both benchmarks —
/// 27.4% on the 63-accession tomato panel, 27.3% on HG002 at 30× and 28.0% at 300×
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
/// **And on a merge table of more than four billion alleles**, where the two `u32` conversions
/// here refuse rather than wrap. That width is not reachable through any input this project can
/// build; it is stated because the conversions are visible in the code and a reader should know
/// they are the only remaining ones.
///
/// **What can no longer panic is a wide locus.** The cap holds the candidate table at
/// [`MaxCandidateAlleles`](super::MaxCandidateAlleles), six by default, so
/// [`CandidateAlleles::admit`]'s refusal at 65,536 alleles is unreachable however many
/// alternatives the merge interned — and the count of cut alternatives is a `u32`, which is why
/// [`SelectionVerdict::Truncated`] carries one.
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

    // `ladder` is named and ignored rather than covered by a `..`, so that a buffer added to the
    // scratch later is a compile error here — the same rule `SelectionScratch::reset_for` states.
    // The ordinary path builds no ladder: a locus's alleles are an unordered set here.
    let SelectionScratch {
        per_allele,
        ranked_table_indices,
        ladder: _,
        sample_reads_per_rung: _,
        promoted_rungs: _,
        rung_is_promoted: _,
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

    // The cap, and the two orders the buffer holds in turn. First the ranking, so that what is
    // cut is the worst-evidenced; then the merge table's own index order, because that is the
    // order the survivors are admitted in (arch §3.1) and so the order that reaches the VCF's
    // `ALT` column. The ranking decides *which* alleles survive and never what order they
    // appear in.
    let allowed_alternatives = usize::from(config.max_candidate_alleles.alternatives());
    let verdict = if ranked_table_indices.len() <= allowed_alternatives {
        SelectionVerdict::Selected
    } else {
        let alternative_of = |table_index: u32| RankedAlternative {
            summary: per_allele[table_index as usize],
            bases: &observation.alleles[table_index as usize],
        };
        ranked_table_indices.sort_unstable_by(|&left, &right| {
            compare_best_first(alternative_of(left), alternative_of(right))
        });
        let dropped = ranked_table_indices.len() - allowed_alternatives;
        ranked_table_indices.truncate(allowed_alternatives);
        ranked_table_indices.sort_unstable();
        SelectionVerdict::Truncated {
            dropped: u32::try_from(dropped).expect("a merge table narrower than four billion"),
        }
    };

    let mut alleles = CandidateAlleles::new(observation.alleles[0].clone(), LocusKind::Generic);
    let mut remap = AlleleRemap::with_all_dropped(observation.alleles.len());
    remap.admit(0, AlleleId::REFERENCE);
    for &table_index in ranked_table_indices.iter() {
        let table_index = table_index as usize;
        let candidate = alleles.admit(observation.alleles[table_index].clone());
        remap.admit(table_index, candidate);
    }

    // The leftover, last: it is per sample where everything above is per allele, and it reads
    // the finished remapping to know what was dropped.
    let covering_samples = observation.per_sample.len();
    let leftovers = observation
        .per_sample
        .iter()
        .map(|sample| {
            leftover_of(
                sample,
                observation.region,
                &remap,
                config.min_allele_support,
            )
        })
        .collect();
    LocusSelection::new(alleles, verdict, leftovers, remap, covering_samples)
}

#[cfg(test)]
mod tests {
    use super::super::MaxCandidateAlleles;
    use super::super::UnmatchedSupport;
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
    /// regime it covers is not a corner: the human benchmark runs at 30× and 300×, so it is where
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
        // A cap wide enough not to bind, so what this test measures is the scratch and not the
        // cap — under the default cap of six the twelve-allele locus is truncated. Thirteen and
        // not twelve, so the cap misses binding by one rather than by nothing.
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.0),
            max_candidate_alleles: cap_of(13),
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

    // ---- the cap (step C2) --------------------------------------------------------------

    /// Narrow `observation` under a support rule of 2 reads and a cap of `cap` alleles
    /// **counting the reference**.
    fn selection_capped_at(
        observation: &CohortObservation,
        max_candidate_alleles: MaxCandidateAlleles,
    ) -> LocusSelection {
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.0),
            max_candidate_alleles,
        };
        select_generic(observation, &config, &mut SelectionScratch::new())
    }

    /// A cap of `alleles` **counting the reference**, said in the type so no call site has to
    /// respell it in words.
    fn cap_of(alleles: u16) -> MaxCandidateAlleles {
        MaxCandidateAlleles::new(alleles).expect("a cap of at least two")
    }

    /// One sample showing every allele at `num_reads[i]` reads, the reference first.
    fn one_sample_showing(num_reads: &[u32]) -> Vec<SampleSupport> {
        let rows = num_reads
            .iter()
            .enumerate()
            .map(|(allele, &num_reads)| row(allele, num_reads, -f64::from(num_reads)))
            .collect();
        vec![sample_showing(0, rows)]
    }

    /// **Below the cap nothing is cut and the verdict is `Selected`** — including at exactly
    /// the cap, which is the boundary a `<` in place of a `<=` would move.
    #[test]
    fn a_locus_at_or_below_the_cap_keeps_everything_that_cleared_the_rule() {
        let three = locus_of(&[b"A", b"C", b"G"], one_sample_showing(&[4, 3, 2]));
        let selection = selection_capped_at(&three, cap_of(3));
        assert_eq!(selection.verdict(), SelectionVerdict::Selected);
        assert_eq!(selection.alleles().len(), 3, "exactly the cap, nothing cut");

        let two = locus_of(&[b"A", b"C"], one_sample_showing(&[4, 3]));
        let selection = selection_capped_at(&two, cap_of(3));
        assert_eq!(selection.verdict(), SelectionVerdict::Selected);
        assert_eq!(selection.alleles().len(), 2);
    }

    /// **Above the cap the list is cut to the best and the locus is still called.** Four
    /// alternatives, a cap of three alleles counting the reference, so two alternatives fit and
    /// two are cut.
    ///
    /// The table holds five alternatives at 8, 3, 6, 4 and 1 reads of one sample's 24, and the
    /// last of them never cleared a rule of two reads. So **four** are candidates for the cap,
    /// two fit under it, and the two kept are the 8 and the 6 — the ranking's first key is the
    /// largest share of one sample's compared reads, and one sample is all there is here.
    ///
    /// **The cut count is two, not three**: the 1-read alternative was dropped by the admission
    /// rule and was never a candidate for the cap, and spec §4.1 makes that distinction because
    /// it is what step C3's second count keys on.
    #[test]
    fn above_the_cap_the_worst_evidenced_alternatives_are_cut_and_the_locus_is_still_called() {
        let observation = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA", b"AC"],
            one_sample_showing(&[2, 8, 3, 6, 4, 1]),
        );
        let selection = selection_capped_at(&observation, cap_of(3));
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 2 },
            "four alternatives cleared the rule, two fit under the cap, and the fifth never \
             cleared it so was never a candidate for it"
        );
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"C".to_vec(), b"T".to_vec()],
            "the reference, and the two best-evidenced alternatives"
        );
    }

    /// **The survivors of a binding cap are admitted in the merge table's order, not in rank
    /// order** — the whole `ALT` column depends on it, and no test before this one could see it.
    ///
    /// A reviewer deleted the sort that puts the kept prefix back into index order and every
    /// test in this file still passed while the column came out permuted. Here the best-ranked
    /// alternative is the *last* of the table's four, so admitting in rank order gives
    /// `AA, C` where the table's order gives `C, AA`.
    #[test]
    fn the_survivors_of_a_binding_cap_are_admitted_in_the_merge_tables_order() {
        let observation = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA"],
            one_sample_showing(&[2, 4, 2, 3, 9]),
        );
        let selection = selection_capped_at(&observation, cap_of(3));
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 2 }
        );
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"C".to_vec(), b"AA".to_vec()],
            "the best-ranked alternative is the table's last, and it stays last"
        );
        assert_eq!(selection.remap().candidate_for(1), Some(AlleleId(1)));
        assert_eq!(selection.remap().candidate_for(4), Some(AlleleId(2)));
        assert_eq!(
            selection.remap().candidate_for(3),
            None,
            "the 3-read alternative is what the cap cut"
        );
    }

    /// **The reference is exempt from the cap as well as from the rule** (spec §6.1): at the
    /// smallest cap a locus can carry — two alleles, one alternative — the reference is still
    /// there and one alternative fits beside it.
    #[test]
    fn the_reference_survives_the_tightest_cap_there_is() {
        // The better-evidenced alternative is the table's *second*, so keeping the merge
        // table's leading prefix would give `C` and fail here.
        let observation = locus_of(&[b"A", b"C", b"G"], one_sample_showing(&[2, 3, 5]));
        let selection = selection_capped_at(&observation, cap_of(MaxCandidateAlleles::SMALLEST));
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 }
        );
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"G".to_vec()],
            "the reference and the better-evidenced of the two alternatives, which is the \
             table's second"
        );
        assert_eq!(selection.alternative_allele_count(), 1);
    }

    /// **The default cap is six alleles counting the reference**, so five alternatives fit and it
    /// first bites at six — the width production's own constant is set at, and which spec §4.2
    /// measures binding at 23 of 53,935 tomato loci and none of HG002's.
    #[test]
    fn the_default_cap_first_bites_at_six_alternatives() {
        let six = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA", b"AC"],
            one_sample_showing(&[2, 7, 6, 5, 4, 3]),
        );
        let selection = select_generic(
            &six,
            &CandidateSelectionConfig::DEFAULT,
            &mut SelectionScratch::new(),
        );
        assert_eq!(selection.verdict(), SelectionVerdict::Selected);
        assert_eq!(selection.alleles().len(), 6);

        // **Every alternative here has to clear the bar comfortably, or this stops being a
        // test of the cap.** An earlier version gave the weakest alternative 2 reads of a
        // sample's 29, and when the shipped share moved from 5 in 100 to 10 the bar — not the
        // cap — removed it: six alternatives became five, five fit, and the verdict came back
        // `Selected` at a fixture whose whole purpose is a cap that bites. Every alternative
        // below takes at least 5 of the sample's 47 reads, where the rule asks for
        // `max(2, ceil(0.10 × 47)) = 5`, so all six reach the cap and the cap is what cuts one.
        let seven = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA", b"AC", b"AG"],
            one_sample_showing(&[2, 5, 8, 7, 6, 9, 10]),
        );
        let selection = select_generic(
            &seven,
            &CandidateSelectionConfig::DEFAULT,
            &mut SelectionScratch::new(),
        );
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 }
        );
        assert_eq!(
            selection.alleles().len(),
            6,
            "six counting the reference, which is five alternatives"
        );
        assert_eq!(
            selection.remap().candidate_for(1),
            None,
            "the 5-read alternative is what the cap cut, and it is the table's first rather \
             than its last — so keeping the table's leading prefix fails here"
        );
    }

    /// **A cap that binds does not change which alleles cleared the rule**, only how many of
    /// them are kept — so the same locus under a wide cap and a tight one gives the tight one's
    /// survivors as a subset of the wide one's, and the same relative order.
    #[test]
    fn the_cap_keeps_a_subset_of_what_a_wider_cap_would_keep() {
        let observation = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA"],
            one_sample_showing(&[2, 9, 3, 7, 5]),
        );
        let wide = surviving_bases(&selection_capped_at(&observation, cap_of(5)));
        let tight = surviving_bases(&selection_capped_at(&observation, cap_of(3)));
        assert_eq!(
            wide,
            vec![
                b"A".to_vec(),
                b"C".to_vec(),
                b"G".to_vec(),
                b"T".to_vec(),
                b"AA".to_vec()
            ]
        );
        assert_eq!(tight, vec![b"A".to_vec(), b"C".to_vec(), b"T".to_vec()]);
        assert!(
            tight.iter().all(|bases| wide.contains(bases)),
            "the cap removes, it never admits"
        );
    }

    /// **The cap is what makes a wide locus safe, at the cohort size it exists for.** Four
    /// hundred samples each carrying a different private allele, every one of them earning it
    /// on its own reads: the merge's table is 401 sequences, the default cap leaves six, and
    /// the verdict counts the 395 it cut.
    ///
    /// **What that costs is not the alleles, it is the samples** — each of the 395 earned the
    /// allele the cap took away, so **395 of the 400 come back with a missing genotype**, which
    /// this test asserts rather than leaves to be discovered. Spec §4.1 measures the cap binding
    /// at 23 of 53,935 tomato loci and reassures that "only the samples that earned a cut allele
    /// are affected"; at several hundred samples that sentence means almost everybody. Raised
    /// with the owner at Checkpoint C.
    #[test]
    fn four_hundred_private_alleles_are_cut_to_the_cap_and_the_locus_is_still_called() {
        // A sequence nobody earns, shown once by every sample, so that all 400 have a non-empty
        // pool while only 395 lose an allele the cap took. Without it, "has a pool" and "the cap
        // cost it something" pick out the same 395 samples and the assertion below cannot tell
        // the two rules apart at the cohort size where the claim is weakest.
        let bases: Vec<Vec<u8>> = std::iter::once(b"A".to_vec())
            .chain((0..400).map(|i| format!("{i:03}").into_bytes()))
            .chain(std::iter::once(b"NOISE".to_vec()))
            .collect();
        let table: Vec<&[u8]> = bases.iter().map(Vec::as_slice).collect();
        let per_sample = (0..400)
            .map(|sample| {
                sample_showing(
                    sample,
                    vec![row(0, 4, -4.0), row(sample + 1, 6, -6.0), row(401, 1, -0.5)],
                )
            })
            .collect();
        let observation = locus_of(&table, per_sample);
        let selection = select_generic(
            &observation,
            &CandidateSelectionConfig::DEFAULT,
            &mut SelectionScratch::new(),
        );
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 395 },
            "400 alternatives cleared the rule and five fit under the default cap; the noise \
             sequence never cleared it and was never a candidate for the cap"
        );
        assert_eq!(selection.alleles().len(), 6);
        assert_eq!(
            surviving_bases(&selection),
            vec![
                b"A".to_vec(),
                b"000".to_vec(),
                b"001".to_vec(),
                b"002".to_vec(),
                b"003".to_vec(),
                b"004".to_vec(),
            ],
            "every alternative ties on all three numeric keys — a share of six in eleven, one \
             clearing sample, six cohort reads — so the bases are the only thing separating 400 \
             of them, and the five kept are the five the bases rank first"
        );
        assert_eq!(selection.unmatched().len(), 400);
        assert_eq!(
            selection
                .unmatched()
                .iter()
                .filter(|leftover| leftover.num_reads > 0)
                .count(),
            400,
            "every sample has a non-empty pool, because every sample showed the noise sequence"
        );
        let missing = selection
            .unmatched()
            .iter()
            .filter(|leftover| leftover.genotype_must_be_missing())
            .count();
        assert_eq!(
            missing, 395,
            "but only the 395 whose own allele the cap cut are emitted as missing — the five \
             whose alleles survived keep their genotypes although their pools are non-empty, \
             which is what separates the cap rule from the pool rule at this cohort size"
        );
    }

    /// **The first ranking key is a share of one sample's reads, not the cohort's total** (spec
    /// §4.1), on the only kind of input that can tell the two apart.
    ///
    /// One sample is heterozygous for `C`, showing it at 15 of its 30 reads — a share of a half.
    /// Ten other samples each show `G` at 2 reads in 30, which is a low-level artefact: a share
    /// of one in fifteen each, but **20 reads across the cohort against the heterozygote's 15**.
    /// At a cap of one alternative, this ranking keeps `C` and production's cohort-read-total
    /// ranking keeps `G`.
    ///
    /// **Every other cap fixture in this file has one sample**, where the share and the cohort
    /// total are the same number over a constant and no fixture can separate the two rankings —
    /// so without this test the whole first key could be replaced by production's and nothing
    /// would fail. That is the key spec §4.1 spends four paragraphs defending, and getting it
    /// wrong cuts the real allele, which step C3 then turns into a missing genotype for the one
    /// sample that carries it.
    #[test]
    fn the_cap_ranks_by_the_within_sample_share_and_not_the_cohort_total() {
        let heterozygote = sample_showing(0, vec![row(0, 15, -15.0), row(1, 15, -15.0)]);
        let artefact_carriers =
            (1..11).map(|sample| sample_showing(sample, vec![row(0, 28, -28.0), row(2, 2, -2.0)]));
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            std::iter::once(heterozygote)
                .chain(artefact_carriers)
                .collect(),
        );
        let selection = selection_capped_at(&observation, cap_of(2));
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 }
        );
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"C".to_vec()],
            "half of one sample's reads beats 20 cohort reads spread one in fifteen"
        );
    }

    /// **At about three reads a position the shares tie and the count of samples decides** — the
    /// tomato panel's whole regime, and spec §4.1's stated reason for the tie-break order.
    ///
    /// Both alternatives take exactly two thirds of their carriers' reads, so the first key
    /// ties. Three samples cleared the rule for `C` and one for `G` — and `G` carries 20 cohort
    /// reads to `C`'s 6, so a ranking that skipped the second key and fell through to the third
    /// would keep the wrong one.
    #[test]
    fn the_cap_falls_through_to_the_sample_count_before_the_cohort_total() {
        let shallow =
            (0..3).map(|sample| sample_showing(sample, vec![row(0, 1, -1.0), row(1, 2, -2.0)]));
        let deep = sample_showing(3, vec![row(0, 10, -10.0), row(2, 20, -20.0)]);
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            shallow.chain(std::iter::once(deep)).collect(),
        );
        let selection = selection_capped_at(&observation, cap_of(2));
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 }
        );
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"C".to_vec()],
            "both shares are two thirds; three samples cleared the rule for C and one for G, \
             and G's larger cohort total is the third key, which never gets asked"
        );
    }

    /// **Where all three numeric keys tie, the bases decide, ascending** — and that is the
    /// regime the cap actually meets at scale, since a cohort of private alleles at one depth
    /// ties on every number the fold records.
    ///
    /// It is also **the only shape in which a summary paired with the wrong allele's bases is
    /// visible**, because the bases are read nowhere else. `RankedAlternative` exists to make
    /// that mis-pairing impossible at a call site, and this test is what checks the call site
    /// the cap actually wrote.
    #[test]
    fn the_cap_at_a_tie_on_every_number_keeps_the_alleles_the_bases_rank_first() {
        let observation = locus_of(
            &[b"A", b"CG", b"CA", b"CT"],
            one_sample_showing(&[4, 3, 3, 3]),
        );
        let selection = selection_capped_at(&observation, cap_of(3));
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 }
        );
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"CG".to_vec(), b"CA".to_vec()],
            "each alternative takes 3 of the sample's 13 reads for one sample, so all three \
             numbers tie and the bases break it ascending: CA and CG survive, CT is cut"
        );
    }

    /// **A scratch reused across a locus the cap truncated.** The cap sorts and truncates the
    /// buffer of surviving indices **in place**, so after a binding cap it holds a short,
    /// reordered list — and the next locus inherits the buffer.
    ///
    /// What makes that safe is `SelectionScratch::reset_for`'s `clear`, which lives in the
    /// parent module and had no test tying it to the in-place truncation the cap introduced. A
    /// stale index surviving into the second locus is either an out-of-range panic or, worse, an
    /// in-range index naming a *different* allele of the new locus — a wrong `ALT` entry with no
    /// other symptom.
    #[test]
    fn a_scratch_reused_after_a_truncated_locus_carries_no_index_into_the_next() {
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.0),
            max_candidate_alleles: cap_of(3),
        };
        let mut scratch = SelectionScratch::new();
        let truncated = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA"],
            one_sample_showing(&[2, 9, 3, 7, 5]),
        );
        let first = select_generic(&truncated, &config, &mut scratch);
        assert_eq!(first.verdict(), SelectionVerdict::Truncated { dropped: 2 });

        let narrow = locus_of(&[b"A", b"CC"], one_sample_showing(&[3, 4]));
        let second = select_generic(&narrow, &config, &mut scratch);
        assert_eq!(second.verdict(), SelectionVerdict::Selected);
        assert_eq!(
            surviving_bases(&second),
            vec![b"A".to_vec(), b"CC".to_vec()],
            "the second locus's own alleles, and nothing left over from the first"
        );
        assert_eq!(second.remap().table_len(), 2);
    }

    // ---- the leftover (step C3) ----------------------------------------------------------

    /// **The plan's first oracle: the pool is the sum of the dropped rows' own `q_sum`, to the
    /// bit** — never re-derived from a count and a rate, which is why the assertion is an
    /// equality against the rows written out in the order the walk visits them.
    ///
    /// One sample, five sequences. A rule of three reads keeps `C` and drops `G`, `T` and `AA`;
    /// `AA` is shown from two read groups, so the pool has to add both of its rows. The
    /// reference is never dropped and its mass is never in the pool.
    ///
    /// **What this equality does and does not pin.** It pins the *set* of rows summed — that the
    /// reference is out, that both of `AA`'s lanes are in, that nothing was scaled — and it
    /// pins the value to the bit. It does **not** pin the addition order, and an earlier version
    /// of this comment claimed it did: the four masses here are exact binary fractions, so every
    /// ordering gives exactly −5.25. Pinning the order would need masses that do not sum
    /// associatively, and there is nothing this step does differently at those.
    #[test]
    fn the_pool_is_the_dropped_rows_own_quality_mass_and_not_a_rate() {
        let observation = locus_of(
            &[b"A", b"C", b"G", b"T", b"AA"],
            vec![sample_showing(
                0,
                vec![
                    row(0, 6, -6.0),
                    row(1, 4, -4.0),
                    row(2, 1, -3.0),
                    row(3, 1, -1.5),
                    row_from_group(4, ReadGroupId(0), 1, -0.5),
                    row_from_group(4, ReadGroupId(1), 1, -0.25),
                ],
            )],
        );
        let selection = selection_of(&observation, support_rule_of(3, 0.0));
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"C".to_vec()],
            "only C reaches three reads — AA's two lanes pool to two, which is one short"
        );
        let leftover = selection.unmatched()[0];
        assert_eq!(
            leftover.num_reads, 4,
            "one read each on G and T, and two on AA across its two read groups"
        );
        assert_eq!(leftover.q_sum, -3.0 + -1.5 + (-0.5 + -0.25));
        assert_eq!(
            leftover.earned_reads_cut_by_the_cap, 0,
            "the rule dropped these, not the cap, so the sample is still callable"
        );
        assert!(!leftover.genotype_must_be_missing());
    }

    /// **A sample with nothing dropped gets a zero pool**, which is the default value and not a
    /// computed one: the walk skips every row, so nothing is added and nothing is rounded.
    /// Every allele it showed survived, so the walk skips every row and the leftover is the
    /// default — which matters because a non-zero third field would emit the sample as missing.
    #[test]
    fn a_sample_whose_alleles_all_survived_has_an_untouched_leftover() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 4, -4.0), row(1, 4, -4.0)]),
                sample_showing(1, vec![row(0, 9, -9.0), row(1, 1, -1.0)]),
            ],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(selection.unmatched()[0], UnmatchedSupport::default());
        assert_eq!(
            selection.unmatched()[1].num_reads,
            0,
            "the second sample's one read is on an allele the first sample's reads admitted, \
             so nothing of its is in the pool either"
        );
    }

    /// **The plan's second oracle, and the one that separates the count from the pool.**
    ///
    /// A locus above the cap. Two samples earned `C` on their own reads — 6 of their 10 each —
    /// and a third had a single error read on it. `G`, which one deeply-sampled accession shows
    /// at 9 of its 10 reads, outranks `C` on the largest within-sample share, so at a cap of one
    /// alternative the cap cuts `C`.
    ///
    /// **The two carriers must come back missing and the third must not.** All three have a
    /// non-zero pool — the third's is its one error read — so a rule keyed on the pool would
    /// no-call all three, and at nearly every locus in a real run it would no-call everybody:
    /// the admission rule drops sequencing error at 13,166 of 15,474 alternatives on the GIAB
    /// HG002 at 300× (spec §3.3), and every sample carries a few error reads almost everywhere.
    /// **The condition is the cap and not the pool** (spec §4.1, §5).
    #[test]
    fn only_the_samples_that_earned_the_cut_allele_are_emitted_as_missing() {
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            vec![
                sample_showing(0, vec![row(0, 4, -4.0), row(1, 6, -6.0)]),
                sample_showing(1, vec![row(0, 4, -4.0), row(1, 6, -6.5)]),
                sample_showing(2, vec![row(0, 9, -9.0), row(1, 1, -1.25)]),
                sample_showing(3, vec![row(0, 1, -1.0), row(2, 9, -9.0)]),
            ],
        );
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.0),
            max_candidate_alleles: cap_of(2),
        };
        let selection = select_generic(&observation, &config, &mut SelectionScratch::new());
        assert_eq!(
            selection.verdict(),
            SelectionVerdict::Truncated { dropped: 1 }
        );
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"G".to_vec()],
            "G takes nine tenths of one sample's reads and C only six tenths, so the cap keeps G"
        );

        let leftovers = selection.unmatched();
        assert_eq!(
            leftovers[0],
            UnmatchedSupport {
                num_reads: 6,
                q_sum: -6.0,
                earned_reads_cut_by_the_cap: 6,
            }
        );
        assert_eq!(
            leftovers[1],
            UnmatchedSupport {
                num_reads: 6,
                q_sum: -6.5,
                earned_reads_cut_by_the_cap: 6,
            }
        );
        assert_eq!(
            leftovers[2],
            UnmatchedSupport {
                num_reads: 1,
                q_sum: -1.25,
                earned_reads_cut_by_the_cap: 0,
            }
        );
        assert_eq!(leftovers[3], UnmatchedSupport::default());

        let missing: Vec<bool> = leftovers
            .iter()
            .map(UnmatchedSupport::genotype_must_be_missing)
            .collect();
        assert_eq!(
            missing,
            vec![true, true, false, false],
            "the two that earned the cut allele are missing; the one with an error read on it \
             is genotyped although its pool is non-zero, and so is the carrier of the survivor"
        );
    }

    /// **Partial reads are not in the pool** (spec §5.1). A read that stopped inside the locus
    /// does not say what the sample carries, it says the sample carries *at least* this, and it
    /// is scored on its own axis — it is not a read matching no candidate, it is a read matching
    /// a *set* of candidates.
    ///
    /// **Nor are the merge's other two counts**: reads that covered the locus and showed nothing,
    /// and reads it removed as evidence. Neither carries a quality sum and neither was ever in
    /// the allele table, so neither can join a pool of error mass.
    #[test]
    fn the_pool_holds_no_partial_no_silent_and_no_removed_read() {
        let mut sample = sample_showing(0, vec![row(0, 6, -6.0), row(1, 1, -1.5)]);
        sample.partials = vec![partial_of(7)];
        sample.reads_without_observation = 5;
        sample.reads_removed_as_evidence = 3;
        let selection = selection_of(
            &locus_of(&[b"A", b"C"], vec![sample]),
            support_rule_of(2, 0.0),
        );
        assert_eq!(
            selection.unmatched()[0],
            UnmatchedSupport {
                num_reads: 1,
                q_sum: -1.5,
                earned_reads_cut_by_the_cap: 0,
            },
            "the one dropped read and its mass, and nothing from the other three axes"
        );
    }

    /// **The reference's reads are never in the pool**, because the reference is never dropped
    /// (spec §5.1, §6.1) — asserted at a locus where the reference carries far more mass than
    /// anything else, so including it would be unmissable.
    #[test]
    fn the_references_own_reads_are_never_in_the_pool() {
        let selection = selection_of(
            &locus_of(
                &[b"A", b"C"],
                vec![sample_showing(0, vec![row(0, 90, -90.0), row(1, 1, -1.0)])],
            ),
            support_rule_of(2, 0.0),
        );
        assert_eq!(selection.unmatched()[0].num_reads, 1);
        assert_eq!(selection.unmatched()[0].q_sum, -1.0);
    }

    /// **A sample's two read groups are pooled before the rule is asked of the leftover, exactly
    /// as they are in the fold.** Here the sample shows the cut allele 3 reads from one lane and
    /// 2 from another: pooled that is 5 and it clears a rule of five, so the sample is missing;
    /// asked per row it would be 3 and 2 and the sample would be genotyped against a table that
    /// does not hold what it carries.
    ///
    /// The two walks must agree because they ask the same rule of the same reads —
    /// `one_run_per_allele` is shared for that reason, and this is the fixture that would notice
    /// if they came apart.
    #[test]
    fn the_leftover_pools_a_samples_read_groups_before_asking_the_rule() {
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            vec![
                sample_showing(
                    0,
                    vec![
                        row(0, 5, -5.0),
                        row_from_group(1, ReadGroupId(0), 3, -3.0),
                        row_from_group(1, ReadGroupId(1), 2, -2.0),
                    ],
                ),
                sample_showing(1, vec![row(0, 1, -1.0), row(2, 9, -9.0)]),
            ],
        );
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.5),
            max_candidate_alleles: cap_of(2),
        };
        let selection = select_generic(&observation, &config, &mut SelectionScratch::new());
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"G".to_vec()],
            "G takes nine tenths of its sample's reads and C a half, so the cap keeps G"
        );
        assert_eq!(
            selection.unmatched()[0],
            UnmatchedSupport {
                num_reads: 5,
                q_sum: -5.0,
                earned_reads_cut_by_the_cap: 5,
            },
            "five reads across two lanes clear a rule of five; asked per lane, three and two \
             would clear nothing and the sample would be genotyped against a table missing \
             what it carries"
        );
        assert!(selection.unmatched()[0].genotype_must_be_missing());
    }
    /// **The rule the leftover re-asks is asked against the sample's *compared reads*, not
    /// against the allele's own** — the same denominator the fold used, which is the whole point
    /// of asking the same rule twice.
    ///
    /// A cut allele, and two samples that show it very differently. The first has 300 compared
    /// reads and 10 on it: a share of 5 in 100 asks for 15, so that sample never earned it and
    /// **must still be genotyped**. The second has 20 compared reads and 15 on it, earns it, and
    /// must come back missing.
    ///
    /// **Asked against the allele's own reads the first sample would be missing too** — 10 reads
    /// against `max(2, ceil(0.05 × 10)) = 2` clears easily — and every other fixture here is
    /// built at a share of zero or a depth where the two denominators agree, so this is the only
    /// one that separates them. At HG002's 30× and 300× that wrong denominator would
    /// no-call every sample with a handful of error reads on a cut allele.
    #[test]
    fn the_leftover_asks_the_rule_against_the_samples_compared_reads() {
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            vec![
                sample_showing(0, vec![row(0, 290, -290.0), row(1, 10, -10.0)]),
                sample_showing(1, vec![row(0, 5, -5.0), row(1, 15, -15.0)]),
                sample_showing(2, vec![row(0, 1, -1.0), row(2, 19, -19.0)]),
            ],
        );
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.05),
            max_candidate_alleles: cap_of(2),
        };
        let selection = select_generic(&observation, &config, &mut SelectionScratch::new());
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"G".to_vec()],
            "G takes 19 of its sample's 20 reads and C 15 of 20, so the cap keeps G"
        );
        assert_eq!(
            selection.unmatched()[0],
            UnmatchedSupport {
                num_reads: 10,
                q_sum: -10.0,
                earned_reads_cut_by_the_cap: 0,
            },
            "ten reads in three hundred miss the share's fifteen, so this sample never earned \
             the cut allele and keeps its genotype"
        );
        assert_eq!(
            selection.unmatched()[1].earned_reads_cut_by_the_cap,
            15,
            "fifteen reads in twenty is an allele this sample really carries"
        );
        let missing: Vec<bool> = selection
            .unmatched()
            .iter()
            .map(UnmatchedSupport::genotype_must_be_missing)
            .collect();
        assert_eq!(missing, vec![false, true, false]);
    }

    /// **A quality mass that is not a number is refused** — the third of the three caller bugs
    /// spec §8 names as assertions this step holds in release, and the one that only becomes
    /// reachable here, because nothing before C3 read `q_sum` at all.
    ///
    /// It is not a crash waiting to happen, which is why it has to be caught rather than left:
    /// the mass flows into the pool, the pool into every genotype's data likelihood, and a
    /// non-finite likelihood prefers no genotype over any other. The locus would come out called
    /// with nothing chosen and nothing failed.
    #[test]
    #[should_panic(expected = "which is not a number the arithmetic can carry")]
    fn a_dropped_allele_whose_quality_mass_is_not_a_number_is_refused() {
        selection_of(
            &locus_of(
                &[b"A", b"C", b"G"],
                vec![sample_showing(
                    0,
                    vec![row(0, 6, -6.0), row(1, 4, -4.0), row(2, 1, f64::NAN)],
                )],
            ),
            support_rule_of(2, 0.0),
        );
    }

    /// **Adding a sample that shows only reference reads changes no other sample's leftover** —
    /// the second half of the plan's standing scale-freedom property, which only became
    /// checkable at C3. The first half, that it changes neither the surviving list nor the
    /// remapping, is asserted at
    /// [`a_sample_showing_only_reference_reads_changes_neither_the_list_nor_the_remapping`].
    ///
    /// It is the property that fails first if a cohort term ever creeps into the leftover — if
    /// the pool were shared, or the missing-genotype condition asked anything about the cohort
    /// rather than about this sample.
    #[test]
    fn a_sample_showing_only_reference_reads_changes_no_other_samples_leftover() {
        let carriers = || {
            vec![
                sample_showing(0, vec![row(0, 6, -6.0), row(1, 4, -4.0), row(2, 1, -1.5)]),
                sample_showing(1, vec![row(0, 8, -8.0), row(2, 1, -0.5)]),
            ]
        };
        let alone = selection_of(
            &locus_of(&[b"A", b"C", b"G"], carriers()),
            support_rule_of(2, 0.0),
        );
        let mut with_a_bystander_samples = carriers();
        with_a_bystander_samples.push(sample_showing(2, vec![row(0, 60, -60.0)]));
        let with_a_bystander = selection_of(
            &locus_of(&[b"A", b"C", b"G"], with_a_bystander_samples),
            support_rule_of(2, 0.0),
        );
        assert_eq!(
            alone.unmatched(),
            &with_a_bystander.unmatched()[..2],
            "the two carriers' leftovers are what they were before the third sample arrived"
        );
        assert_eq!(
            with_a_bystander.unmatched()[2],
            UnmatchedSupport::default(),
            "and the bystander showed only the reference, which is never dropped"
        );
    }

    /// **The pool and the earned count are two different numbers, and this is the only fixture
    /// where they differ.** The sample carries one error read on `T`, which the admission rule
    /// drops, and six on `C`, which it earns and the cap cuts: its pool is seven reads and the
    /// reads the cap took from it are six.
    ///
    /// Every other leftover fixture gives the affected sample **exactly one** dropped allele,
    /// where the two totals coincide — so without this one the earned count could be the running
    /// pool total and nothing would fail. And `T` sits *before* `C` in the merge table on
    /// purpose: after `C` the running total and the right answer coincide again.
    #[test]
    fn the_earned_count_is_this_alleles_reads_and_not_the_pool_so_far() {
        let observation = locus_of(
            &[b"A", b"T", b"C", b"G"],
            vec![
                sample_showing(0, vec![row(0, 4, -4.0), row(1, 1, -1.0), row(2, 6, -6.0)]),
                sample_showing(1, vec![row(0, 1, -1.0), row(3, 9, -9.0)]),
            ],
        );
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.0),
            max_candidate_alleles: cap_of(2),
        };
        let selection = select_generic(&observation, &config, &mut SelectionScratch::new());
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"G".to_vec()],
            "G takes nine tenths of its sample's reads and C six elevenths, so the cap keeps G"
        );
        assert_eq!(
            selection.unmatched()[0],
            UnmatchedSupport {
                num_reads: 7,
                q_sum: -7.0,
                earned_reads_cut_by_the_cap: 6,
            },
            "seven reads in the pool, six of them on the allele this sample earned"
        );
    }

    /// **The leftover is parallel to `per_sample` by position, and a partial-only sample holds
    /// its place.** The samples either side drop the same one read of `C` but carry different
    /// error mass, so a leftover that slid by one position shows up in the values.
    ///
    /// **The length assertion in `LocusSelection::new` cannot see that**, because a fill that
    /// skips a sample and pads the tail is still the right length — and the failure it hides is
    /// the wrong sample emitted as missing while the sample that lost an allele is genotyped
    /// against a table that does not hold it.
    #[test]
    fn the_leftover_stays_aligned_when_a_partial_only_sample_sits_between_two_others() {
        let observation = locus_of(
            &[b"A", b"C"],
            vec![
                sample_showing(0, vec![row(0, 5, -5.0), row(1, 1, -1.0)]),
                sample_with_only_partials(1, 4),
                sample_showing(2, vec![row(0, 5, -5.0), row(1, 1, -7.0)]),
            ],
        );
        let selection = selection_of(&observation, support_rule_of(2, 0.0));
        assert_eq!(
            selection.unmatched(),
            &[
                UnmatchedSupport {
                    num_reads: 1,
                    q_sum: -1.0,
                    earned_reads_cut_by_the_cap: 0,
                },
                UnmatchedSupport::default(),
                UnmatchedSupport {
                    num_reads: 1,
                    q_sum: -7.0,
                    earned_reads_cut_by_the_cap: 0,
                },
            ],
            "the middle sample has only partials and drops nothing, and the third sample's mass \
             must not slide onto it"
        );
    }

    /// **The denominator is every compared read, the cut allele's own included.** The first
    /// sample holds 6 reads of `C` in 16 compared and the rule asks 8, so it never earned `C`
    /// and keeps its genotype. Against the 10 reads it has *left* after `C` is cut, the same
    /// share asks only 5 and the sample would be no-called.
    ///
    /// **That wrong denominator is systematically smaller than the right one** — it is the right
    /// one minus whatever selection dropped — so it asks less of every sample and no-calls
    /// samples the rule never meant to touch. Every other leftover fixture is built at a share of
    /// zero or at a depth where the two denominators give the same answer.
    #[test]
    fn the_leftover_asks_the_rule_against_every_compared_read_including_the_cut_alleles_own() {
        let observation = locus_of(
            &[b"A", b"C", b"G"],
            vec![
                sample_showing(0, vec![row(0, 10, -10.0), row(1, 6, -6.0)]),
                sample_showing(1, vec![row(0, 2, -2.0), row(1, 8, -8.0)]),
                sample_showing(2, vec![row(0, 1, -1.0), row(2, 9, -9.0)]),
            ],
        );
        let config = CandidateSelectionConfig {
            min_allele_support: support_rule_of(2, 0.5),
            max_candidate_alleles: cap_of(2),
        };
        let selection = select_generic(&observation, &config, &mut SelectionScratch::new());
        assert_eq!(
            surviving_bases(&selection),
            vec![b"A".to_vec(), b"G".to_vec()],
            "G takes nine tenths of its sample's reads and C four fifths, so the cap keeps G"
        );
        assert_eq!(
            selection.unmatched()[0],
            UnmatchedSupport {
                num_reads: 6,
                q_sum: -6.0,
                earned_reads_cut_by_the_cap: 0,
            },
            "six reads in sixteen compared misses the share's eight"
        );
        assert!(
            selection.unmatched()[1].genotype_must_be_missing(),
            "eight reads in ten earns it, and the cap took it away"
        );
    }
}
