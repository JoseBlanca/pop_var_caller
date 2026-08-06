//! One locus reduced to one cell key: how many reads covered the site, and how many
//! of those showed something other than the reference.
//!
//! **The only place that decides what counts as an alternative read**, which is why
//! it is its own file rather than a method on the locus type — a locus cannot know a
//! model's answer to that question
//! (`arch/parameter_prepass_generic.md` §2.3).
//!
//! The answer here is byte equality against the locus's own reference bases, and it
//! reduces to production's `allele_index == 0` rule exactly where the two can be
//! compared. Three things follow from it, and each is a decision rather than an
//! implementation detail — see [`count_whole_site`].

use smallvec::SmallVec;

use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::ReadGroupId;

use super::histogram::DepthAndAltReads;

/// Count one generic locus's reads, taking the site **whole** — the spec's word for the
/// un-split site (`spec/parameter_prepass_generic.md` §4). This is the pair the windowed
/// histogram enters, because a site's genotype is one thing however many read groups
/// covered it.
///
/// **What counts as an alternative read.** A complete witness whose bases differ from
/// the locus's reference bases, byte for byte. That is the per-read form of production's
/// `allele_index == 0` rule and it needs no special cases: `bases` is in *read*
/// coordinates, so an insertion makes it longer than the reference slice and a deletion
/// shorter, and either way it cannot compare equal — which is the right answer for both.
/// The reference bases arrive canonical `{A, C, G, T, N}` from
/// [`RefSeq::fetch_into`](crate::ng::ref_seq::RefSeq::fetch_into), so a soft-masked
/// lowercase reference cannot make every read at a repeat look alternative; the verbatim
/// reader that preserves soft-masking is a different capability, used by the typed-region
/// catalog and not on this path.
///
/// **Complete witnesses only.** A read that spanned part of the locus witnessed neither
/// the reference allele nor an alternative one at the positions it missed, and scoring it
/// as either is the mis-scoring `complete_observations` exists to guard. Loci one base
/// wide — the overwhelming majority on this path — lose nothing, since a read either
/// covers the base or does not.
///
/// **`reads_without_observation` does not enter the depth.** Those reads covered the
/// locus and witnessed nothing, so counting them would assert they showed the reference,
/// and they did not show anything. The depth here is reads whose evidence is present,
/// which is what every likelihood in this design conditions on. They are counted
/// separately by the accumulator instead: a locus where many reads say nothing is one
/// whose depth is not what it appears, and that is a fact about the mapping rather than
/// about the genotype. It also has to be this way for
/// `spec/parameter_prepass_generic.md` §12.6 to hold at all — the field is a bare scalar
/// with no read-group attribution, so including it would make the two histograms differ
/// by construction on a single-library sample.
///
/// **`reads_discarded_by_cap` does not skip the locus.** An earlier draft skipped any
/// locus whose reads an upstream cap had already subsampled. That is depth-dependent
/// selection of exactly the kind this step exists to remove: ng's generic pileup caps
/// indel-bearing columns at 250 reads against 8,000 for the rest, so at 300× every
/// indel-bearing column would be dropped and at 30× none would, and the coverage-invariance
/// anchor of §9 would be measuring the skip rule rather than the estimator. The locus is
/// entered at the depth observed and the accumulator counts it
/// (`AccumulationCounts::loci_with_upstream_subsample`).
///
/// # Panics
///
/// If the locus's supports sum past `u32`, which no pileup cap allows.
#[must_use]
pub fn count_whole_site(locus: &SampleLocusObservations) -> DepthAndAltReads {
    let mut depth = 0u32;
    let mut alt_reads = 0u32;

    for observation in locus.complete_observations() {
        depth = add_support(depth, observation.num_obs, "depth");
        if *observation.bases != *locus.reference_bases {
            alt_reads = add_support(alt_reads, observation.num_obs, "alternative reads");
        }
    }

    DepthAndAltReads::new(depth, alt_reads)
}

/// The same count, split between the read groups that covered the site — one entry per
/// group, ascending by [`ReadGroupId`].
///
/// This is what the read-group histogram enters: a site with 20 reads from one library
/// and 10 from another becomes two entries, one at depth 20 and one at depth 10, because
/// an error rate describes the chemistry and two libraries of one sample can genuinely
/// differ. The split costs that fit nothing, since what it counts is reads.
///
/// **It is also where the windowed histogram's library attribution comes from**, whose
/// depth is the *sum* of these — exact rather than approximate, since these are raw
/// counts at one position (`spec/parameter_prepass_generic.md` §4).
///
/// `out` is scratch: cleared, then filled. It is an argument rather than a return value
/// because this runs once per covered position over hundreds of millions of them, and a
/// fresh `Vec` per locus would be an allocation on the hottest path in the pre-pass.
///
/// What counts as an alternative read, which witnesses are read, and what is left out of
/// the depth are all [`count_whole_site`]'s doc — the two functions differ only in grain.
///
/// # Panics
///
/// If any group's supports sum past `u32`, which no pileup cap allows.
pub fn count_by_read_group(
    locus: &SampleLocusObservations,
    out: &mut Vec<(ReadGroupId, DepthAndAltReads)>,
) {
    // Depth and alternative count are accumulated as a bare pair and only made a
    // `DepthAndAltReads` at the end, because that type's constructor checks
    // `alt_reads <= depth` and a half-accumulated group has not earned the check yet.
    // Inline for two groups, which is every sample that has more than one library at
    // all: 1,550 of the 1,707 in the tomato archive survey have exactly one.
    let mut running: SmallVec<[(ReadGroupId, u32, u32); 2]> = SmallVec::new();

    for observation in locus.complete_observations() {
        let group = observation.read_group;
        let at = match running.iter().position(|&(seen, _, _)| seen == group) {
            Some(at) => at,
            None => {
                running.push((group, 0, 0));
                running.len() - 1
            }
        };
        let entry = &mut running[at];
        entry.1 = add_support(entry.1, observation.num_obs, "depth");
        if *observation.bases != *locus.reference_bases {
            entry.2 = add_support(entry.2, observation.num_obs, "alternative reads");
        }
    }

    running.sort_unstable_by_key(|&(group, _, _)| group);
    out.clear();
    out.extend(
        running
            .into_iter()
            .map(|(group, depth, alt_reads)| (group, DepthAndAltReads::new(depth, alt_reads))),
    );
}

/// Add one observation's support, refusing to wrap.
///
/// The release profile leaves `overflow-checks` off, so a wrapped depth would come back
/// as a small number and be scored as a shallow site. It cannot happen through the
/// generic pileup, which caps a column at 8,000 reads — this is the module's rule that a
/// counter says so rather than wrapping quietly.
fn add_support(running: u32, support: u32, what: &str) -> u32 {
    running.checked_add(support).unwrap_or_else(|| {
        panic!("a locus's {what} passed u32 at {running}, with {support} more to add")
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, LocusLen, ReadWitness, SequenceObservation};
    use crate::ng::types::{ContigId, GenomeRegion, Position};

    fn group(id: u32) -> ReadGroupId {
        ReadGroupId(id)
    }

    /// One observation, with only the fields this step reads set to anything meaningful.
    fn observation(
        bases: &[u8],
        witness: ReadWitness,
        read_group: ReadGroupId,
        num_obs: u32,
    ) -> SequenceObservation {
        SequenceObservation {
            bases: bases.into(),
            read_witness: witness,
            read_group,
            num_obs,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn locus(
        reference_bases: &[u8],
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        let span = reference_bases.len().max(1) as u64;
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(1),
                end: Position(span),
            },
            reference_bases: reference_bases.into(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    fn by_group(locus: &SampleLocusObservations) -> Vec<(ReadGroupId, DepthAndAltReads)> {
        let mut out = Vec::new();
        count_by_read_group(locus, &mut out);
        out
    }

    /// **The one-base site, which is the overwhelming majority of this path.** Eleven
    /// reads showed the reference base and three showed something else, so the site is
    /// entered at depth 14 with 3 alternative reads.
    #[test]
    fn a_one_base_site_counts_the_reads_that_disagreed_with_the_reference() {
        let site = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(0), 11),
                observation(b"C", ReadWitness::Complete, group(0), 3),
            ],
        );

        let counted = count_whole_site(&site);
        assert_eq!(counted.depth(), 14);
        assert_eq!(counted.alt_reads(), 3);
    }

    /// **An insertion and a deletion are alternative reads without either being a
    /// special case.** `bases` is in read coordinates, so an inserted base makes the
    /// observation longer than the reference slice and a deletion shorter; neither can
    /// compare equal, which is the right answer for both. A rule that compared only
    /// same-length sequences, or that compared position by position, would need to name
    /// them — and would then have to decide what a length mismatch means.
    #[test]
    fn an_insertion_and_a_deletion_are_both_alternative_reads() {
        let widened = locus(
            b"ACGT",
            vec![
                observation(b"ACGT", ReadWitness::Complete, group(0), 20),
                observation(b"ACCGT", ReadWitness::Complete, group(0), 4), // an insertion
                observation(b"AGT", ReadWitness::Complete, group(0), 6),   // a deletion
            ],
        );

        let counted = count_whole_site(&widened);
        assert_eq!(counted.depth(), 30);
        assert_eq!(counted.alt_reads(), 10);
    }

    /// **A read that witnessed part of the locus witnessed neither allele at the
    /// positions it missed**, so it enters neither the depth nor the alternative count.
    /// Scoring it as the reference would assert it agreed where it saw nothing; scoring
    /// it as alternative would assert the opposite. Both are the mis-scoring
    /// `complete_observations` exists to guard.
    #[test]
    fn a_partial_witness_enters_neither_the_depth_nor_the_alternative_count() {
        // Two of the four locus positions, witnessed from the left — a read that ran out
        // before the locus's right border.
        let saw_two_of_four =
            ReadWitness::from_left(2, LocusLen::from_positions(4)).expect("a non-empty run");
        let widened = locus(
            b"ACGT",
            vec![
                observation(b"ACGT", ReadWitness::Complete, group(0), 9),
                observation(b"AC", saw_two_of_four.clone(), group(0), 5),
                observation(b"AG", saw_two_of_four, group(0), 7),
            ],
        );

        let counted = count_whole_site(&widened);
        assert_eq!(counted.depth(), 9, "only the complete witness");
        assert_eq!(counted.alt_reads(), 0);
    }

    /// A locus every read covered without witnessing anything is entered at depth zero,
    /// not skipped and not counted as reference. `reads_without_observation` is a bare
    /// scalar with no read-group attribution, so it could not be split between the two
    /// histograms even if it were wanted — including it would break the cell-for-cell
    /// equality of `spec/parameter_prepass_generic.md` §12.6 by construction.
    #[test]
    fn reads_that_witnessed_nothing_do_not_enter_the_depth() {
        let mut silent = locus(b"A", Vec::new());
        silent.reads_without_observation = 30;

        let counted = count_whole_site(&silent);
        assert_eq!(counted.depth(), 0);
        assert_eq!(counted.alt_reads(), 0);
        assert!(by_group(&silent).is_empty(), "no group covered it");
    }

    /// A locus whose reads an upstream cap already subsampled is entered **at the depth
    /// observed**, not skipped. Skipping it would be depth-dependent selection: ng's
    /// pileup caps indel-bearing columns at 250 reads and everything else at 8,000, so
    /// at 300× every indel-bearing column would go and at 30× none would.
    #[test]
    fn a_locus_an_upstream_cap_subsampled_is_entered_at_the_depth_observed() {
        let mut capped = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(0), 240),
                observation(b"T", ReadWitness::Complete, group(0), 10),
            ],
        );
        capped.reads_discarded_by_cap = 1_400;

        let counted = count_whole_site(&capped);
        assert_eq!(
            counted.depth(),
            250,
            "the reads that are here, not the 1,650"
        );
        assert_eq!(counted.alt_reads(), 10);
    }

    /// **The read-group split, and the entries come back in read-group order** whatever
    /// order the observations arrived in — the histogram keys on this, and an unsorted
    /// listing would split one site's cell in two.
    #[test]
    fn the_read_group_split_comes_back_in_read_group_order() {
        let site = locus(
            b"A",
            vec![
                observation(b"T", ReadWitness::Complete, group(7), 2),
                observation(b"A", ReadWitness::Complete, group(3), 18),
                observation(b"A", ReadWitness::Complete, group(7), 8),
                observation(b"T", ReadWitness::Complete, group(3), 2),
            ],
        );

        let split = by_group(&site);
        assert_eq!(split.len(), 2);
        assert_eq!(split[0].0, group(3));
        assert_eq!((split[0].1.depth(), split[0].1.alt_reads()), (20, 2));
        assert_eq!(split[1].0, group(7));
        assert_eq!((split[1].1.depth(), split[1].1.alt_reads()), (10, 2));
    }

    /// **The two grains agree, which is the property the two histograms rest on.** The
    /// windowed table enters a site once at its total depth and the read-group table
    /// once per group; if the split did not sum to the whole, the two tables would
    /// describe different data and §12.6's cell-for-cell equality on a single-library
    /// sample could not hold.
    ///
    /// **Every rule the two functions share has to appear in at least one of these
    /// loci, or this cannot see them drifting apart.** Making `count_whole_site` alone
    /// read the partial witnesses left it green while the two grains disagreed by
    /// twelve reads, because no fixture had one — so the last locus carries a partial
    /// witness, a widened span and two libraries at once.
    #[test]
    fn the_read_group_split_sums_to_the_whole_site() {
        let saw_two_of_four =
            ReadWitness::from_left(2, LocusLen::from_positions(4)).expect("a non-empty run");
        for site in [
            locus(
                b"A",
                vec![
                    observation(b"A", ReadWitness::Complete, group(0), 11),
                    observation(b"C", ReadWitness::Complete, group(1), 3),
                    observation(b"A", ReadWitness::Complete, group(1), 7),
                    observation(b"G", ReadWitness::Complete, group(4), 1),
                ],
            ),
            locus(
                b"A",
                vec![observation(b"A", ReadWitness::Complete, group(2), 5)],
            ),
            locus(b"ACGT", Vec::new()),
            locus(
                b"ACGT",
                vec![
                    observation(b"ACGT", ReadWitness::Complete, group(0), 9),
                    observation(b"ACCGT", ReadWitness::Complete, group(1), 4),
                    observation(b"AC", saw_two_of_four.clone(), group(0), 5),
                    observation(b"AG", saw_two_of_four.clone(), group(1), 7),
                ],
            ),
        ] {
            let whole = count_whole_site(&site);
            let split = by_group(&site);

            let depth: u32 = split.iter().map(|(_, counted)| counted.depth()).sum();
            let alt_reads: u32 = split.iter().map(|(_, counted)| counted.alt_reads()).sum();
            assert_eq!(depth, whole.depth());
            assert_eq!(alt_reads, whole.alt_reads());
        }
    }

    /// A single-library site produces exactly one entry, which is the case that matters
    /// most: 1,550 of the 1,707 samples in the tomato archive survey carry one library,
    /// and for them the read-group table and the windowed one hold the same numbers.
    #[test]
    fn a_single_library_site_produces_one_entry_equal_to_the_whole() {
        let site = locus(
            b"A",
            vec![
                observation(b"A", ReadWitness::Complete, group(6), 11),
                observation(b"C", ReadWitness::Complete, group(6), 3),
            ],
        );

        let split = by_group(&site);
        assert_eq!(split.len(), 1);
        assert_eq!(split[0].0, group(6));
        assert_eq!(split[0].1, count_whole_site(&site));
    }

    /// The scratch vector is cleared before it is filled, so a caller reusing one across
    /// loci does not accumulate the previous locus's groups.
    #[test]
    fn the_scratch_vector_is_cleared_between_loci() {
        let mut out = Vec::new();
        let first = locus(
            b"A",
            vec![observation(b"C", ReadWitness::Complete, group(9), 4)],
        );
        count_by_read_group(&first, &mut out);
        assert_eq!(out.len(), 1);

        let second = locus(
            b"A",
            vec![observation(b"A", ReadWitness::Complete, group(1), 6)],
        );
        count_by_read_group(&second, &mut out);
        assert_eq!(out.len(), 1, "the previous locus's group is gone");
        assert_eq!(out[0].0, group(1));
        assert_eq!(out[0].1.alt_reads(), 0);
    }

    /// **The reference bases arrive canonical, and this is what would break if they did
    /// not.** ng has two reference readers: `RefSeq::fetch_into` returns canonical
    /// `{A, C, G, T, N}`, and the raw reader preserves soft-masked lowercase for the
    /// typed-region catalog's byte oracle. A locus built from the second would call
    /// every read at a soft-masked position alternative — about half the human genome,
    /// and the repeat-rich half. This pins the comparison as byte equality so the
    /// dependency is visible rather than implied.
    #[test]
    fn byte_equality_is_the_rule_and_a_lowercase_reference_would_break_it() {
        let canonical = locus(
            b"A",
            vec![observation(b"A", ReadWitness::Complete, group(0), 12)],
        );
        assert_eq!(count_whole_site(&canonical).alt_reads(), 0);

        let soft_masked = locus(
            b"a",
            vec![observation(b"A", ReadWitness::Complete, group(0), 12)],
        );
        assert_eq!(
            count_whole_site(&soft_masked).alt_reads(),
            12,
            "every read would read as alternative — which is why the generic locus \
             generator fetches through RefSeq and not the raw reader"
        );
    }
}
