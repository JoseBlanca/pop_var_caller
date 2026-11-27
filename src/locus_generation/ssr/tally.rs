//! The per-locus tally: fold each kept read's [`Classified`] into the locus's observed
//! sequences, one per distinct [`ObservationIdentity`], summing the support moments and counting
//! the reasons reads yielded no observation.
//!
//! It follows the deleted caller's `tally`, extended with partial observations and the
//! strand, base-quality and mapping-quality moments the shared locus type carries (spec §3, §6).

use super::SsrGeneratorCounts;
use super::classify::{Classified, NoObservationReason};
use crate::bam::alignment_input::FLAG_REVERSE_STRAND;
use crate::locus_generation::{
    FlankDeletion, HeldFlank, ReadWitness, SequenceObservation, read_placed_left,
};
use crate::read::aligned_read::AlignedRead;
use crate::read::prepared_read::alignment_end_of;
use crate::types::ReadGroupId;
use std::collections::HashMap;

/// What makes two reads one observation: their bases, witness, read group, flank deletion
/// and the flank they held where the witness does not say it.
type ObservationIdentity = (
    Box<[u8]>,
    ReadWitness,
    ReadGroupId,
    FlankDeletion,
    Option<HeldFlank>,
);

/// The per-locus tally the generator folds onto the `SampleLocusObservations`: the deduped
/// observed sequences and how many reads reached the aligner but yielded nothing.
/// (`reads_fetched` / `reads_discarded_by_cap` are the caller's to set — they come from the
/// cap, not the outcomes.)
pub(super) struct SsrTally {
    pub(super) observations: Vec<SequenceObservation>,
    pub(super) reads_without_observation: u32,
}

/// Accumulated support for one distinct [`ObservationIdentity`] — the moments summed as reads
/// fold in, materialised into a [`SequenceObservation`] at the end.
#[derive(Default)]
struct Support {
    num_obs: u32,
    num_fwd: u32,
    q_sum: f64,
    mapq_sum: u32,
    mapq_sum_sq: u64,
    placed_left: u32,
}

/// Fold each kept read's classification into the locus tally.
///
/// Observations dedup by **[`ObservationIdentity`]** — bases, witness, read group, flank
/// deletion and held flank. So a `Complete` and a partial of the same bases are different
/// evidence and stay separate observations, two identical partials from one read group merge,
/// and the **same allele seen from two read groups is two observations** (spec §3, §6). Each
/// bucket accumulates the strand (`num_fwd` off the reverse-strand flag), BQ (`q_sum`, off the
/// read's tract error mass), MAPQ and `placed_left` moments; `chain_ids` stays empty because
/// the STR path does not phase. `observations` is sorted by all five fields of that identity,
/// so — like the deleted caller's `tally` — the bases, the counts
/// and the integer moments are independent of the order reads were folded. `q_sum` is the one
/// exception: it sums in `f64`, which is commutative but not associative, so a bucket's
/// `q_sum` can differ in its low bits under a different fold order. That is immaterial while
/// `q_sum` is soft and unconsumed, and the parity oracle checks only bytes and counts (spec
/// §6).
///
/// `locus_start` and `locus_end` are the tract's first and last base, 1-based and inclusive,
/// needed only for `placed_left` — the count of supporting reads with at least as much of
/// their span left of the tract as right of it ([`read_placed_left`]). The whole tract is
/// the locus here, unlike the generic path's anchor, because it is fixed before any read
/// is measured.
///
/// `counts` carries the **run-level** totals (complete/partial observations and the four
/// no-observation reasons); the returned `reads_without_observation` is this locus's own
/// total.
pub(super) fn tally<'a>(
    reads_and_outcomes: impl IntoIterator<Item = (&'a AlignedRead, Classified)>,
    locus_start: u64,
    locus_end: u64,
    counts: &mut SsrGeneratorCounts,
) -> SsrTally {
    let mut buckets: HashMap<ObservationIdentity, Support> = HashMap::new();
    let mut reads_without_observation = 0u32;
    for (read, outcome) in reads_and_outcomes {
        match outcome {
            Classified::Observed {
                bases,
                read_witness,
                flank_deletion,
                held_flank,
                q_sum,
            } => {
                match read_witness {
                    ReadWitness::Complete => counts.observations_complete += 1,
                    ReadWitness::Partial { .. } => counts.observations_partial += 1,
                }
                let support = buckets
                    .entry((
                        bases,
                        read_witness,
                        read.read_group,
                        flank_deletion,
                        held_flank,
                    ))
                    .or_default();
                support.num_obs += 1;
                if read.flag & FLAG_REVERSE_STRAND == 0 {
                    support.num_fwd += 1;
                }
                support.q_sum += q_sum;
                let mapq = u32::from(read.mapq);
                support.mapq_sum += mapq;
                support.mapq_sum_sq += u64::from(mapq) * u64::from(mapq);
                // `pos` fits in 32 bits on any contig the generic path accepts, which
                // converts it the same way (`prepare_passthrough`).
                let alignment_end = alignment_end_of(read.pos as u32, &read.cigar);
                support.placed_left += u32::from(read_placed_left(
                    read.pos,
                    u64::from(alignment_end),
                    locus_start,
                    locus_end,
                ));
            }
            Classified::NoObservation(reason) => {
                reads_without_observation += 1;
                match reason {
                    NoObservationReason::NoBorderAnchored => counts.no_border_anchored += 1,
                    NoObservationReason::LowQuality => counts.low_quality += 1,
                    NoObservationReason::WindowTruncated => counts.window_truncated += 1,
                    NoObservationReason::OutsideTract => counts.outside_tract += 1,
                }
            }
        }
    }

    let mut observations: Vec<SequenceObservation> = buckets
        .into_iter()
        .map(
            |((bases, read_witness, read_group, flank_deletion, held_flank), support)| {
                SequenceObservation {
                    bases,
                    read_witness,
                    read_group,
                    num_obs: support.num_obs,
                    num_fwd: support.num_fwd,
                    q_sum: crate::types::SummedLogError::from_nats(support.q_sum),
                    mapq_sum: support.mapq_sum,
                    mapq_sum_sq: support.mapq_sum_sq,
                    placed_left: support.placed_left,
                    // The STR path has no artifact penalties to feed (`low_base_quality.md`).
                    num_low_quality: 0,
                    // The STR path does not phase — there is no chain id to fold (spec §3).
                    chain_ids: Vec::new(),
                    flank_deletion,
                    held_flank,
                }
            },
        )
        .collect();
    // The read group joins the sort key because it joined the bucket key: without it two
    // observations differing only by group would tie, and `HashMap` iteration order is seeded per
    // process — the output would be non-deterministic run to run on any multi-group sample.
    observations.sort_unstable_by(|a, b| {
        a.bases
            .cmp(&b.bases)
            .then_with(|| a.read_witness.sort_key().cmp(&b.read_witness.sort_key()))
            .then_with(|| a.read_group.cmp(&b.read_group))
            .then_with(|| a.flank_deletion.cmp(&b.flank_deletion))
            .then_with(|| a.held_flank.cmp(&b.held_flank))
    });

    SsrTally {
        observations,
        reads_without_observation,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam::alignment_input::CigarOp;
    use crate::locus_generation::LocusLen;
    use crate::types::ReadGroupId;

    /// An `AlignedRead` with a given strand flag and MAPQ — the only fields the tally reads
    /// off the read; the sequence and qualities live in the `Classified` handed alongside.
    fn read(flag: u16, mapq: u8) -> AlignedRead {
        AlignedRead {
            qname: b"r".to_vec(),
            flag,
            ref_id: 0,
            pos: 1,
            mapq,
            cigar: vec![CigarOp::Match(6)],
            seq: b"CACACA".to_vec(),
            qual: vec![40; 6],
            mate_ref_id: None,
            mate_pos: None,
            adaptor_boundary: None,
            read_group: ReadGroupId(0),
        }
    }

    /// The same, in a named read group and starting at `pos` — the two fields the
    /// `(…, read_group)` bucket key and `placed_left` read.
    fn read_in_group(group: u32, pos: u64) -> AlignedRead {
        AlignedRead {
            pos,
            read_group: ReadGroupId(group),
            ..read(0, 60)
        }
    }

    fn observed(bases: &[u8], read_witness: ReadWitness, q_sum: f64) -> Classified {
        Classified::Observed {
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            held_flank: None,
            bases: bases.into(),
            read_witness,
            q_sum,
        }
    }

    /// **The read group is part of the identity**: one allele seen from two read groups is
    /// two observations, and their `num_obs` sum to what a single-group tally would have reported
    /// (spec §6). This is the check that the split is *computed* and not defaulted — with
    /// `read_group` left at a constant both reads would merge into one observation of two.
    #[test]
    fn one_allele_from_two_read_groups_is_two_observations_that_sum_back() {
        let rg0 = read_in_group(0, 10);
        let rg1 = read_in_group(1, 10);
        let outcomes = vec![
            (&rg0, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            (&rg0, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            (&rg1, observed(b"CACACA", ReadWitness::Complete, -1.0)),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let observations = tally(outcomes, 10, 10, &mut counts).observations;

        assert_eq!(
            observations.len(),
            2,
            "one observation per (allele, read group)"
        );
        assert!(
            observations
                .iter()
                .all(|observation| &*observation.bases == b"CACACA"),
            "both observations are the same allele"
        );
        assert_eq!(
            observations
                .iter()
                .map(|observation| observation.read_group)
                .collect::<Vec<_>>(),
            vec![ReadGroupId(0), ReadGroupId(1)],
            "and they are ordered by group, which is what keeps the output deterministic"
        );
        assert_eq!(
            observations
                .iter()
                .map(|observation| observation.num_obs)
                .sum::<u32>(),
            3,
            "collapsing the group axis recovers the single-group total exactly"
        );
    }

    /// **The group tie-break is asserted on the sort key, not on a lucky fold order.**
    ///
    /// The obvious form of this test — fold two groups and check the observation order — rides on
    /// `HashMap` iteration, which is seeded per process: with the `then_with` deleted it
    /// passes roughly six runs in ten, so the regression it exists to catch would reach a
    /// green CI most of the time. Enough observations that a wrong order cannot be a coin flip is
    /// the fix: six groups over two alleles, whose sorted order is fully determined.
    #[test]
    fn observations_sort_by_group_within_an_allele_deterministically() {
        let reads: Vec<AlignedRead> = (0..6).map(|g| read_in_group(g, 10)).collect();
        let mut outcomes = Vec::new();
        for (i, r) in reads.iter().enumerate() {
            // Alternate the alleles so group order and fold order disagree.
            let bases: &[u8] = if i % 2 == 0 { b"AA" } else { b"CC" };
            outcomes.push((r, observed(bases, ReadWitness::Complete, -1.0)));
        }
        let mut counts = SsrGeneratorCounts::default();
        let observations = tally(outcomes, 10, 10, &mut counts).observations;

        let order: Vec<(&[u8], u32)> = observations
            .iter()
            .map(|observation| (observation.bases.as_ref(), observation.read_group.get()))
            .collect();
        assert_eq!(
            order,
            vec![
                (b"AA".as_ref(), 0),
                (b"AA".as_ref(), 2),
                (b"AA".as_ref(), 4),
                (b"CC".as_ref(), 1),
                (b"CC".as_ref(), 3),
                (b"CC".as_ref(), 5),
            ],
            "bases first, then read group — ascending, and never fold order"
        );
    }

    /// **An expanded allele merges the two sides into one observation, and that is a real
    /// behaviour change the reshape brought.**
    ///
    /// `reach` is measured in *read* bases; on an allele longer than the reference tract it
    /// exceeds the locus length, so `from_left` and `from_right` both clamp to the whole
    /// locus and produce the **same run**. Two reads anchored at opposite borders with the
    /// same bases then share a bucket key and merge — where `PartialLeft(n)` and
    /// `PartialRight(n)` kept them as two observations.
    ///
    /// It is arguably the right answer (identical constraints are one observation) but it is not
    /// the pre-reshape answer, it is invisible on any fixture whose reads are exact
    /// reference slices, and the plan's stated equivalence
    /// `PartialRight(n) ⇔ Partial { len - n, n }` silently stops holding at `n = len`.
    /// Pinned here so it is a decision on the record rather than a surprise in a dump.
    ///
    /// *Since 2026-10-12 the classifier writes the side down on such a read
    /// ([`SequenceObservation::held_flank`]), at segments as at bundles, so two real reads
    /// from opposite borders no longer merge
    /// (`partials_held_from_two_flanks_stay_two_observations`). What this still pins is the
    /// witness: once saturated, it alone cannot tell the sides apart.*
    #[test]
    fn an_expanded_allele_merges_the_two_sides_into_one_observation() {
        let locus_len = LocusLen::from_positions(6);
        let left =
            ReadWitness::from_left(9, locus_len).expect("a run covering at least one position");
        let right =
            ReadWitness::from_right(9, locus_len).expect("a run covering at least one position");
        assert_eq!(
            left, right,
            "the two sides denote the same run once saturated"
        );

        let r = read(0, 60);
        let outcomes = vec![
            (&r, observed(b"CACACACACA", left, -1.0)),
            (&r, observed(b"CACACACACA", right, -1.0)),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let observations = tally(outcomes, 1, 1, &mut counts).observations;

        assert_eq!(
            observations.len(),
            1,
            "one observation, not two — the sides are indistinguishable once the run saturates"
        );
        assert_eq!(observations[0].num_obs, 2, "and both reads support it");
        assert_eq!(
            counts.observations_partial, 2,
            "the run-level per-read tally is unaffected: two reads, two partials"
        );
    }

    /// `placed_left` counts supporting reads with at least as many aligned bases left of
    /// the tract as right of it — freebayes' `placedLeft` — and not reads that merely
    /// started left of it, which is nearly every read that covers a tract.
    ///
    /// Every read here starts left of the tract `[20, 25]`, so the old rule counted all
    /// four. Two of them lie mostly to its right. The two in the middle are the boundary:
    /// ten bases each side is a tie and goes left; nine against eleven goes right.
    #[test]
    fn placed_left_counts_reads_lying_mostly_left_of_the_tract() {
        fn spanning(pos: u64, len: u32) -> AlignedRead {
            AlignedRead {
                cigar: vec![CigarOp::Match(len)],
                ..read_in_group(0, pos)
            }
        }
        let mostly_right = spanning(18, 30); // 2 bases left, 22 right
        let mostly_left = spanning(1, 30); // 19 left, 5 right
        let tie = spanning(10, 26); // 10 left, 10 right
        let one_short = spanning(11, 26); // 9 left, 11 right
        let outcomes = vec![
            (
                &mostly_right,
                observed(b"CACACA", ReadWitness::Complete, -1.0),
            ),
            (
                &mostly_left,
                observed(b"CACACA", ReadWitness::Complete, -1.0),
            ),
            (&tie, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            (&one_short, observed(b"CACACA", ReadWitness::Complete, -1.0)),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let observations = tally(outcomes, 20, 25, &mut counts).observations;

        assert_eq!(
            observations.len(),
            1,
            "one allele, one group — one observation"
        );
        assert_eq!(observations[0].num_obs, 4);
        assert_eq!(
            observations[0].placed_left, 2,
            "the read with 19 bases left and the tie are placed left; the two with more \
             bases right are not, though all four started left of the tract"
        );
    }

    /// A `Complete` and a partial run of the **same** bases are different evidence, so they
    /// stay as two separate observations (spec §3) — the property the [`ObservationIdentity`]
    /// dedup key rests on.
    #[test]
    fn a_complete_and_a_partial_of_the_same_bases_stay_separate() {
        let fwd = read(0, 60);
        let outcomes = vec![
            (&fwd, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            (
                &fwd,
                observed(
                    b"CACACA",
                    ReadWitness::from_left(6, LocusLen::from_positions(6))
                        .expect("a run covering at least one position"),
                    -1.0,
                ),
            ),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let result = tally(outcomes, 1, 1, &mut counts);
        assert_eq!(result.observations.len(), 2);
        assert_eq!(counts.observations_complete, 1);
        assert_eq!(counts.observations_partial, 1);
    }

    /// Two identical partials (same bases, same witness) are the identical constraint, so
    /// they merge into one observation with `num_obs == 2` (spec §3).
    #[test]
    fn two_identical_partials_merge_into_one_count() {
        let fwd = read(0, 60);
        let outcomes = vec![
            (
                &fwd,
                observed(
                    b"CACA",
                    ReadWitness::from_left(4, LocusLen::from_positions(6))
                        .expect("a run covering at least one position"),
                    -1.0,
                ),
            ),
            (
                &fwd,
                observed(
                    b"CACA",
                    ReadWitness::from_left(4, LocusLen::from_positions(6))
                        .expect("a run covering at least one position"),
                    -1.0,
                ),
            ),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let result = tally(outcomes, 1, 1, &mut counts);
        assert_eq!(result.observations.len(), 1);
        assert_eq!(result.observations[0].num_obs, 2);
    }

    /// **The held flank is part of the identity**: one sequence over a whole bundle hull read
    /// from the left and from the right is two pieces of evidence, kept apart and carried
    /// onto the observations (`str_bundle_loci.md` §7.6).
    #[test]
    fn partials_held_from_two_flanks_stay_two_observations() {
        let fwd = read(0, 60);
        let whole_hull = || {
            ReadWitness::from_left(8, LocusLen::from_positions(6))
                .expect("a run covering at least one position")
        };
        let held = |flank| Classified::Observed {
            bases: b"CACACACA"[..].into(),
            read_witness: whole_hull(),
            flank_deletion: FlankDeletion::NONE,
            held_flank: Some(flank),
            q_sum: -1.0,
        };
        let outcomes = vec![
            (&fwd, held(HeldFlank::Right)),
            (&fwd, held(HeldFlank::Left)),
            (&fwd, held(HeldFlank::Left)),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let result = tally(outcomes, 1, 1, &mut counts);
        let flanks: Vec<_> = result
            .observations
            .iter()
            .map(|observation| (observation.held_flank, observation.num_obs))
            .collect();
        assert_eq!(
            flanks,
            vec![(Some(HeldFlank::Left), 2), (Some(HeldFlank::Right), 1)]
        );
    }

    /// The strand, BQ and MAPQ moments accumulate across the reads folded into a bucket:
    /// `num_fwd` counts the forward-strand reads, `q_sum` sums the per-read BQ mass,
    /// `mapq_sum` / `mapq_sum_sq` sum MAPQ and MAPQ². `chain_ids` stays empty.
    #[test]
    fn the_support_moments_accumulate() {
        let fwd = read(0, 60);
        let rev = read(FLAG_REVERSE_STRAND, 30);
        let outcomes = vec![
            (&fwd, observed(b"CACACA", ReadWitness::Complete, -2.0)),
            (&rev, observed(b"CACACA", ReadWitness::Complete, -3.0)),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let result = tally(outcomes, 1, 1, &mut counts);
        let obs = &result.observations[0];
        assert_eq!(obs.num_obs, 2);
        assert_eq!(obs.num_fwd, 1, "one forward, one reverse");
        assert_eq!(obs.q_sum.nats(), -5.0, "the BQ masses sum");
        assert_eq!(obs.mapq_sum, 90, "60 + 30");
        assert_eq!(obs.mapq_sum_sq, 60 * 60 + 30 * 30);
        assert!(obs.chain_ids.is_empty(), "the STR path does not phase");
    }

    /// Each no-observation reason lands in its own run-level counter, and every such read is
    /// counted in the locus's `reads_without_observation` total.
    #[test]
    fn no_observation_reasons_are_counted_by_reason_and_in_total() {
        let r = read(0, 60);
        let outcomes = vec![
            (
                &r,
                Classified::NoObservation(NoObservationReason::NoBorderAnchored),
            ),
            (
                &r,
                Classified::NoObservation(NoObservationReason::LowQuality),
            ),
            (
                &r,
                Classified::NoObservation(NoObservationReason::WindowTruncated),
            ),
            (
                &r,
                Classified::NoObservation(NoObservationReason::LowQuality),
            ),
            // The outside-the-tract reason, the largest of the four on real data — its tally arm went
            // untested until the Milestone C review said so.
            (
                &r,
                Classified::NoObservation(NoObservationReason::OutsideTract),
            ),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let result = tally(outcomes, 1, 1, &mut counts);
        assert!(result.observations.is_empty());
        assert_eq!(result.reads_without_observation, 5);
        assert_eq!(counts.no_border_anchored, 1);
        assert_eq!(counts.low_quality, 2);
        assert_eq!(counts.window_truncated, 1);
        assert_eq!(counts.outside_tract, 1);
    }

    /// `observations` is sorted by bytes, then by witness — so the record is identical
    /// regardless of the order reads folded in (production's order-independence, extended to
    /// the witness tie-break).
    #[test]
    fn observations_are_sorted_by_bases_then_witness() {
        let r = read(0, 60);
        let outcomes = vec![
            (&r, observed(b"GG", ReadWitness::Complete, -1.0)),
            (
                &r,
                observed(
                    b"AA",
                    ReadWitness::from_right(2, LocusLen::from_positions(6))
                        .expect("a run covering at least one position"),
                    -1.0,
                ),
            ),
            (&r, observed(b"AA", ReadWitness::Complete, -1.0)),
            (
                &r,
                observed(
                    b"AA",
                    ReadWitness::from_left(2, LocusLen::from_positions(6))
                        .expect("a run covering at least one position"),
                    -1.0,
                ),
            ),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let result = tally(outcomes, 1, 1, &mut counts);
        let order: Vec<(&[u8], ReadWitness)> = result
            .observations
            .iter()
            .map(|o| (o.bases.as_ref(), o.read_witness.clone()))
            .collect();
        assert_eq!(
            order,
            vec![
                (b"AA".as_ref(), ReadWitness::Complete),
                (
                    b"AA".as_ref(),
                    ReadWitness::from_left(2, LocusLen::from_positions(6))
                        .expect("a run covering at least one position")
                ),
                (
                    b"AA".as_ref(),
                    ReadWitness::from_right(2, LocusLen::from_positions(6))
                        .expect("a run covering at least one position")
                ),
                (b"GG".as_ref(), ReadWitness::Complete),
            ]
        );
    }

    /// **Which of the tie-break's two components dominates** — the claim
    /// [`ReadWitness::sort_key`]'s doc makes, and which the test above cannot check.
    ///
    /// `observations_are_sorted_by_bases_then_witness` builds both of its partials with the
    /// **same** `positions_covered` (`from_left(2, len 6)` = `{0,2}`,
    /// `from_right(2, len 6)` = `{4,2}`), so exchanging the two components of the sort key
    /// maps them to `(1,2,0)` and `(1,2,4)` and the asserted order survives by accident.
    /// The review made exactly that mutation and the whole suite stayed green.
    ///
    /// Here the long run is the left-flush one and the short run the right-flush one, so
    /// "offset outranks length" and "shortest first" disagree — which is the only shape
    /// that can fail. It matters because on the STR path the two are different genetic
    /// constraints, a prefix and a suffix, and their emission order is what a cohort merge
    /// reads.
    #[test]
    fn tally_orders_two_partials_of_one_sequence_by_offset_before_length() {
        let r = read(0, 60);
        let len = LocusLen::from_positions(6);
        let outcomes = vec![
            (
                &r,
                observed(
                    b"AA",
                    ReadWitness::from_right(2, len).expect("a run covering at least one position"),
                    -1.0,
                ),
            ),
            (
                &r,
                observed(
                    b"AA",
                    ReadWitness::from_left(4, len).expect("a run covering at least one position"),
                    -1.0,
                ),
            ),
        ];
        let mut counts = SsrGeneratorCounts::default();
        let result = tally(outcomes, 1, 1, &mut counts);
        let order: Vec<ReadWitness> = result
            .observations
            .iter()
            .map(|o| o.read_witness.clone())
            .collect();
        assert_eq!(
            order,
            vec![
                ReadWitness::from_left(4, len).expect("a run covering at least one position"),
                ReadWitness::from_right(2, len).expect("a run covering at least one position"),
            ],
            "a left-flush run must precede a right-flush one whatever their lengths",
        );
    }

    /// The integer moments of a **multi-read bucket** are order-independent: the same reads
    /// (differing strand and MAPQ) folded into one bucket in two orders
    /// give the same `num_obs` / `num_fwd` / `mapq_sum` / `mapq_sum_sq`. This is the case the
    /// singleton buckets of `tally_is_order_independent` never exercise, and it isolates
    /// `q_sum` as the sole order-sensitive field (its f64 sum is not asserted here).
    #[test]
    fn a_multi_read_bucket_accumulates_integer_moments_order_independently() {
        let fwd = read(0, 60);
        let rev = read(FLAG_REVERSE_STRAND, 30);
        let moments = |outcomes: Vec<(&AlignedRead, Classified)>| {
            let mut counts = SsrGeneratorCounts::default();
            let obs = tally(outcomes, 1, 1, &mut counts).observations;
            assert_eq!(obs.len(), 1, "all reads share one bucket");
            let o = &obs[0];
            (o.num_obs, o.num_fwd, o.mapq_sum, o.mapq_sum_sq)
        };
        let forward_first = moments(vec![
            (&fwd, observed(b"CACACA", ReadWitness::Complete, -2.0)),
            (&rev, observed(b"CACACA", ReadWitness::Complete, -3.0)),
        ]);
        let reverse_first = moments(vec![
            (&rev, observed(b"CACACA", ReadWitness::Complete, -3.0)),
            (&fwd, observed(b"CACACA", ReadWitness::Complete, -2.0)),
        ]);
        assert_eq!(forward_first, reverse_first);
        assert_eq!(forward_first, (2, 1, 90, 60 * 60 + 30 * 30));
    }

    /// The tally is order-independent: the same multiset of outcomes folded in two different
    /// orders yields the same observed sequences and the same counts. (Its buckets are
    /// singletons, so no `q_sum` f64 reordering occurs — the multi-read integer case is
    /// `a_multi_read_bucket_accumulates_integer_moments_order_independently`.)
    #[test]
    fn tally_is_order_independent() {
        let r = read(0, 60);
        let run = |outcomes: Vec<(&AlignedRead, Classified)>| {
            let mut counts = SsrGeneratorCounts::default();
            let result = tally(outcomes, 1, 1, &mut counts);
            (result.observations, counts)
        };
        let a = run(vec![
            (&r, observed(b"CACACA", ReadWitness::Complete, -1.0)),
            (&r, observed(b"CACA", ReadWitness::Complete, -1.0)),
            (
                &r,
                Classified::NoObservation(NoObservationReason::LowQuality),
            ),
        ]);
        let b = run(vec![
            (
                &r,
                Classified::NoObservation(NoObservationReason::LowQuality),
            ),
            (&r, observed(b"CACA", ReadWitness::Complete, -1.0)),
            (&r, observed(b"CACACA", ReadWitness::Complete, -1.0)),
        ]);
        assert_eq!(a, b);
    }
}
