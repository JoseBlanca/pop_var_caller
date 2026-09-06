//! Which reference positions one record reports depth at, and how much (spec §3.1).
//!
//! The accumulator next door takes one covered position at a time. This is what turns a
//! sample's records into that stream: given one record as the merge's cache drew it, it calls
//! back once per position the record speaks for, with the depth there.
//!
//! **Depth here is observation depth** — the number of reads whose observation covers the
//! position. It is below read depth by the reads that covered the ground and produced no
//! observation and by any a depth cap discarded, consistently on both sides of the filter's
//! division, which is what spec §6's first trap is about: the yardstick and the measurement
//! have to be the same quantity.
//!
//! **Three rules:**
//!
//! - **A record spanning one base reports `reads_compared_with_reference`** — the reads whose
//!   whole sequence over the locus was compared against the reference. Every record's summary
//!   carries it: psp mode reads it from the stored record's head, direct mode derives it from
//!   the record the walker just minted. At a one-base locus no read's witness can stop inside
//!   the locus, so every observation is whole and that count is the sum of `num_obs` — which is
//!   why psp mode need not decode the evidence here at all. **That equality is a claim about
//!   real data, not a theorem** (spec §6 trap 7); a counter-example — a read whose evidence
//!   stops inside a one-base locus — moves this rule to "build every body". Plan step B2
//!   checked **8,784,182 one-base records and found no disagreement and no such read**, on two
//!   tomato stores: six accessions over 200 kb of one chromosome at a mean 14.4 reads compared
//!   with the reference, and one accession over 8 Mb of all twelve at 10.3. **This rule decides
//!   992 positions in every 1,000 there, and the build rule the other 8** — a share that moves
//!   with indel density and with how much of the ground is repeat tract, so it is a fact about
//!   those files rather than about the caller.
//! - **A generic record spanning more than one base reports at its first base only.** The
//!   generic walk emits one record per covered position, and a record widened by a deletion is
//!   anchored at its first base while the interior positions have records of their own — so
//!   spreading its depth over its span would count that interior twice (spec §6 trap 2).
//! - **A repeat-tract record reports at every position of its span.** Region typing partitions
//!   the reference and a tract is one segment, so a tract's record is the only record over its
//!   ground and nothing else will report those positions.
//!
//! The last two need the evidence, which direct mode already has in hand and psp mode builds
//! on demand from the bytes it kept. Which of the three applies is decided by the record's
//! **span**, never by which shape the draw arrived in, so both modes compute the same number
//! from the same rule (spec §1.1 goal 5).

use core::ops::Range;

use crate::ng::locus_generation::{LocusKind, SampleLocusObservations};
use crate::ng::run::cohort_merge::observation_cache::{Drawn, LocusSummary};
use crate::ng::types::{GenomePosition, Position};

/// Call `report` once for every reference position `drawn` speaks for, with this sample's
/// depth there, in increasing position order and all on the record's own contig.
///
/// `summary` must be `drawn`'s own — **checked in release, not only in debug**, because a
/// summary belonging to another record would silently attribute one record's ground to
/// another's and the accumulator downstream would fill the wrong window. The parameter exists
/// so that a built record's sequences are not walked again here: the cache derives the summary
/// as it draws and holds it, and deriving it once more is a walk of every observation the
/// record carries. A kept record needs no such economy — its summary travels inside the draw —
/// so that one is used and the parameter is only checked against it.
///
/// `build` is how the evidence is obtained for a record that arrived without it, and it is
/// **called at most once, and only for a record spanning more than one base**: a one-base
/// record is answered from its summary, so psp mode never decodes the bytes it kept for one —
/// on the stores B2 measured, all but about one record in a thousand (the module doc above).
///
/// A position with a record but no reads is reported at depth `0`. That is a covered position
/// — the sample has a record there — and dropping it would leave the window's denominator
/// counting only the positions that happened to have depth.
///
/// **A record whose region ends before it starts reports nothing at all.** Such a record
/// describes no ground: `GenomeRegion::len` saturates to zero for it and the depth vector comes
/// back empty. Nothing mints one today; the behaviour is pinned by a test rather than left to
/// be discovered.
///
/// # Errors
///
/// Whatever `build` refuses, unchanged: this adds nothing to it, because the source that
/// failed is the one that can name itself. Nothing is reported for a record whose evidence
/// could not be built, so a caller's accumulator is never left holding half of one.
pub fn for_each_reported_depth<E>(
    drawn: &Drawn,
    summary: LocusSummary,
    build: impl FnOnce(Range<usize>) -> Result<SampleLocusObservations, E>,
    mut report: impl FnMut(GenomePosition, u32),
) -> Result<(), E> {
    // **A release check and not a `debug_assert!`** — the release profile is the one this repo
    // runs (`Cargo.toml`'s `[profile.release]` leaves `debug-assertions` off, and
    // `[profile.soak]` exists to turn them back on), and the same reasoning put a release
    // assert on the cache's own ordering check one file away. A mismatch here is not a caught
    // error but a silent one: the positions would come from this summary and the depths from
    // another record.
    let summary = match drawn {
        Drawn::Built(record) => {
            assert_eq!(
                record.region, summary.region,
                "the summary handed in is not this record's",
            );
            summary
        }
        // The summary a kept draw carries is this record's by construction, so it is the one
        // used; the parameter is checked against it whole, which covers the head count as well
        // as the ground — and the head count is the number a one-base record reports.
        Drawn::Kept {
            summary: its_own, ..
        } => {
            assert_eq!(
                *its_own, summary,
                "the summary handed in is not this record's",
            );
            *its_own
        }
    };

    let first_base = summary.start_position();
    if first_base.position == summary.reach() {
        report(first_base, summary.reads_compared_with_reference);
        return Ok(());
    }

    match drawn {
        Drawn::Built(record) => report_body_depths(record, first_base, &mut report),
        // `summary: _` rather than `..`: the summary this arm has already used is right there
        // in the variant, and naming it is what makes ignoring it here a decision. A third
        // field on `Drawn::Kept` then fails to compile at the one site that walks every record.
        Drawn::Kept { body, summary: _ } => {
            let record = build(body.clone())?;
            // A fact about the file, in psp mode: a stored head that claims other ground than
            // the body behind it. At plan step C2, where the source doing the building is a psp
            // reader, this becomes a `RunError` naming the sample, beside the other
            // producer-guarantee checks the merge already makes.
            assert_eq!(
                record.region, summary.region,
                "the body built for this record covers other ground than its head claimed",
            );
            report_body_depths(&record, first_base, &mut report);
        }
    }
    Ok(())
}

/// The positions a built record reports, and the depth at each, decided by its kind.
///
/// `first_base` is the record's first position, taken from its summary — **one derivation, not
/// two**, so the position a depth is reported at cannot come from a different reading of the
/// region than the decision to report it at all. The asserts above are what let this trust it:
/// by the time this is called, the record's region and the summary's are equal.
///
/// **The depth is always `num_obs_along_locus`'s**, even where only its first entry is kept, so
/// that a generic record's anchor and a tract's positions are the same quantity derived by the
/// same code. Taking the anchor's depth straight from the observations would save an allocation
/// and a write at every position of the span for every observation, and would leave two
/// derivations of "depth" that agree until a witness kind neither author thought about — trap
/// 1's failure one level down.
fn report_body_depths(
    record: &SampleLocusObservations,
    first_base: GenomePosition,
    report: &mut impl FnMut(GenomePosition, u32),
) {
    let depth_along_the_locus = record.num_obs_along_locus();
    match record.kind {
        // Its interior positions have records of their own.
        LocusKind::Generic => {
            if let Some(&depth_at_the_anchor) = depth_along_the_locus.first() {
                report(first_base, depth_at_the_anchor);
            }
        }
        // Nothing else reports this ground. `SsrBundle` regions emit no loci today
        // (`locus_generation.md` §11), so no bundle record reaches this arm; it is grouped
        // with the tract because a bundle is a repeat segment too — one record over ground
        // that is partitioned to it — and the alternative, reporting a bundle at its anchor
        // alone, would silently under-cover every window overlapping it if such records ever
        // start arriving.
        LocusKind::Ssr(_) | LocusKind::SsrBundle => {
            for (offset, &depth) in depth_along_the_locus.iter().enumerate() {
                // The add cannot overflow: `offset` is below the region's length, so the
                // largest position produced is the region's own last base.
                report(
                    GenomePosition {
                        contig: first_base.contig,
                        position: Position(first_base.position.get() + offset as u64),
                    },
                    depth,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{
        LocusLen, ReadWitness, SequenceObservation, SsrDetail, WitnessedLocusPositions,
    };
    use crate::ng::types::{ContigId, GenomeRegion, Motif, ReadGroupId, SummedLogError};

    /// The contig every fixture here sits on; which one is never the point.
    const CONTIG: ContigId = ContigId(3);

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: CONTIG,
            start: Position(start),
            end: Position(end),
        }
    }

    fn at(position: u64) -> GenomePosition {
        GenomePosition {
            contig: CONTIG,
            position: Position(position),
        }
    }

    /// An observation of `num_obs` reads at a given witness. Every field the depth rule does
    /// not read is fixed, so a fixture's numbers are only the ones under test.
    fn obs(read_witness: ReadWitness, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(&b"A"[..]),
            read_witness,
            read_group: ReadGroupId(0),
            num_obs,
            num_fwd: 0,
            q_sum: SummedLogError::from_nats(0.0),
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn generic_record(
        region: GenomeRegion,
        observed: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region,
            reference_bases: Box::from(&b""[..]),
            observations: observed,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    fn tract_record(
        region: GenomeRegion,
        observed: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            kind: LocusKind::Ssr(SsrDetail {
                motif: Motif::new(b"AT").expect("a two-base motif"),
                left_flank: Box::from(&b"CCC"[..]),
                right_flank: Box::from(&b"GGG"[..]),
            }),
            ..generic_record(region, observed)
        }
    }

    /// The summary a record's head carries, with the head count said out loud rather than
    /// derived — which is the point at a one-base record, where the rule reads the summary and
    /// never the evidence.
    fn summary_of(region: GenomeRegion, reads_compared_with_reference: u32) -> LocusSummary {
        LocusSummary {
            region,
            non_reference_reads: 0,
            reads_compared_with_reference,
        }
    }

    /// A `build` that must not be called, and says which record it was asked for if it is.
    fn build_refused(body: Range<usize>) -> Result<SampleLocusObservations, &'static str> {
        panic!(
            "the evidence was built for a record that should have been answered from its summary: {body:?}"
        );
    }

    /// Collect what the rule reports for a drawn record whose evidence is already in hand.
    fn collect_reported(drawn: &Drawn, summary: LocusSummary) -> Vec<(GenomePosition, u32)> {
        collect_reported_building_with(drawn, summary, |body| -> Result<_, &'static str> {
            panic!("no body was kept for this record, yet one was asked for: {body:?}")
        })
    }

    fn collect_reported_building_with<E: core::fmt::Debug>(
        drawn: &Drawn,
        summary: LocusSummary,
        build: impl FnOnce(Range<usize>) -> Result<SampleLocusObservations, E>,
    ) -> Vec<(GenomePosition, u32)> {
        let mut reported = Vec::new();
        for_each_reported_depth(drawn, summary, build, |position, depth| {
            reported.push((position, depth));
        })
        .expect("this fixture's build does not fail");
        reported
    }

    /// The one-base rule: the summary's count is reported, and the evidence is not consulted
    /// even when it is in hand and says something else. The record's observations sum to 9;
    /// the summary says 7; 7 is what a one-base record reports.
    #[test]
    fn a_record_spanning_one_base_reports_the_head_count_and_not_its_evidence() {
        let record = generic_record(
            region(42, 42),
            vec![obs(ReadWitness::Complete, 7), obs(ReadWitness::Complete, 2)],
        );
        assert_eq!(record.num_obs_along_locus(), vec![9]);

        let reported = collect_reported(&Drawn::Built(record), summary_of(region(42, 42), 7));
        assert_eq!(reported, vec![(at(42), 7)]);
    }

    /// The same rule on the shape that has no evidence in hand: a kept one-base record must
    /// not reach for its body, which is what leaves psp mode's stored bytes undecoded wherever
    /// a record covers a single base.
    #[test]
    fn a_kept_record_spanning_one_base_builds_no_body() {
        let summary = summary_of(region(42, 42), 7);
        let drawn = Drawn::Kept {
            summary,
            body: 100..180,
        };
        let mut reported = Vec::new();
        for_each_reported_depth(&drawn, summary, build_refused, |position, depth| {
            reported.push((position, depth));
        })
        .expect("the summary answers this record");
        assert_eq!(reported, vec![(at(42), 7)]);
    }

    /// A generic record widened by a deletion reports at its first base alone: the interior
    /// positions have records of their own, and reporting them here would count them twice.
    #[test]
    fn a_generic_record_spanning_more_than_one_base_reports_at_its_anchor_only() {
        let record = generic_record(
            region(10, 12),
            vec![
                obs(ReadWitness::Complete, 3),
                obs(
                    ReadWitness::from_left(1, LocusLen::from_positions(3))
                        .expect("a run covering at least one position"),
                    5,
                ),
            ],
        );
        assert_eq!(record.num_obs_along_locus(), vec![8, 3, 3]);

        // The head count is deliberately none of 8, 3 or their sum, so that a rule reading the
        // summary where it should read the body fails this test rather than passing by
        // coincidence.
        let reported = collect_reported(&Drawn::Built(record), summary_of(region(10, 12), 99));
        assert_eq!(reported, vec![(at(10), 8)]);
    }

    /// The anchor's depth is the depth **at** the anchor, not the record's largest. A read
    /// whose evidence starts inside a widened record leaves the first position shallower than
    /// the interior, and the fixture above cannot tell the two rules apart because its anchor
    /// is also its deepest position.
    #[test]
    fn a_generic_records_anchor_depth_is_its_first_position_and_not_its_deepest() {
        let record = generic_record(
            region(10, 12),
            vec![obs(
                ReadWitness::from_right(2, LocusLen::from_positions(3))
                    .expect("a run covering at least one position"),
                5,
            )],
        );
        assert_eq!(record.num_obs_along_locus(), vec![0, 5, 5]);

        let reported = collect_reported(&Drawn::Built(record), summary_of(region(10, 12), 99));
        assert_eq!(reported, vec![(at(10), 0)]);
    }

    /// A repeat tract is the only record over its ground, so every position of its span is
    /// reported — and at the depth the witnesses give there, not at one number spread flat.
    #[test]
    fn a_repeat_tract_reports_every_position_of_its_span() {
        let record = tract_record(
            region(10, 12),
            vec![
                obs(ReadWitness::Complete, 3),
                obs(
                    ReadWitness::from_right(2, LocusLen::from_positions(3))
                        .expect("a run covering at least one position"),
                    4,
                ),
            ],
        );
        assert_eq!(record.num_obs_along_locus(), vec![3, 7, 7]);

        let reported = collect_reported(&Drawn::Built(record), summary_of(region(10, 12), 99));
        assert_eq!(reported, vec![(at(10), 3), (at(11), 7), (at(12), 7)]);
    }

    /// The evidence a kept record left behind is built through the source, and the bytes it
    /// is asked for are the ones the draw named.
    #[test]
    fn a_kept_record_spanning_more_than_one_base_is_built_through_the_source() {
        let summary = summary_of(region(10, 12), 99);
        let drawn = Drawn::Kept {
            summary,
            body: 4_096..4_224,
        };
        let mut asked_for = None;
        let reported = collect_reported_building_with(&drawn, summary, |body| {
            asked_for = Some(body);
            Ok::<_, &'static str>(tract_record(
                region(10, 12),
                vec![obs(ReadWitness::Complete, 6)],
            ))
        });
        assert_eq!(asked_for, Some(4_096..4_224));
        assert_eq!(reported, vec![(at(10), 6), (at(11), 6), (at(12), 6)]);
    }

    /// A **generic** record that arrived kept reports at its anchor alone, exactly as the built
    /// shape does. Without this, psp mode could spread a deletion-widened record's depth over
    /// its span — double-counting every deletion interior (spec §6 trap 2) at exactly the loci
    /// direct mode gets right — and the rest of the suite would stay green.
    #[test]
    fn a_kept_generic_record_spanning_more_than_one_base_reports_at_its_anchor_only() {
        let summary = summary_of(region(10, 12), 99);
        let record = || {
            generic_record(
                region(10, 12),
                vec![
                    obs(ReadWitness::Complete, 3),
                    obs(
                        ReadWitness::from_left(1, LocusLen::from_positions(3))
                            .expect("a run covering at least one position"),
                        5,
                    ),
                ],
            )
        };
        let kept = collect_reported_building_with(
            &Drawn::Kept {
                summary,
                body: 0..1,
            },
            summary,
            |_| Ok::<_, &'static str>(record()),
        );
        assert_eq!(kept, vec![(at(10), 8)]);
        assert_eq!(kept, collect_reported(&Drawn::Built(record()), summary));
    }

    /// **Mode equivalence in miniature.** The same record reported through the shape a walker
    /// over alignment files hands over and through the shape a reader over stored files hands
    /// over gives the same positions and the same depths — which is what stops direct mode and
    /// psp mode training two different yardsticks (spec §1.1 goal 5).
    #[test]
    fn both_draw_shapes_report_the_same_depths_for_the_same_record() {
        let summary = summary_of(region(10, 12), 99);
        let record = || {
            tract_record(
                region(10, 12),
                vec![
                    obs(ReadWitness::Complete, 3),
                    obs(
                        ReadWitness::from_left(2, LocusLen::from_positions(3))
                            .expect("a run covering at least one position"),
                        4,
                    ),
                ],
            )
        };

        let built = collect_reported(&Drawn::Built(record()), summary);
        let kept = collect_reported_building_with(
            &Drawn::Kept {
                summary,
                body: 0..1,
            },
            summary,
            |_| Ok::<_, &'static str>(record()),
        );
        assert_eq!(built, kept);
        assert_eq!(built, vec![(at(10), 7), (at(11), 7), (at(12), 3)]);
    }

    /// A body that cannot be built stops the walk and hands its own error back untouched.
    /// Absorbed instead, a psp decode failure would take the record's positions out of the
    /// window's denominator and out of the sample's histogram while the run reported success.
    #[test]
    fn a_build_that_fails_stops_the_walk_and_returns_its_error() {
        let summary = summary_of(region(10, 12), 99);
        let drawn = Drawn::Kept {
            summary,
            body: 0..1,
        };
        let mut reported = Vec::new();
        let outcome = for_each_reported_depth(
            &drawn,
            summary,
            |_| Err::<SampleLocusObservations, _>("the source could not read this record"),
            |position, depth| reported.push((position, depth)),
        );
        assert_eq!(outcome, Err("the source could not read this record"));
        assert!(
            reported.is_empty(),
            "a record whose evidence failed to build reported {reported:?}, so a caller's \
             accumulator would hold half of it",
        );
    }

    /// A record with no reads behind it is still a covered position, reported at depth zero.
    /// Dropping it would leave the window's denominator counting only positions that happened
    /// to have depth, and the mean would be over a set nothing else knows about.
    #[test]
    fn a_record_with_no_reads_is_reported_at_depth_zero() {
        let one_base = collect_reported(
            &Drawn::Built(generic_record(region(7, 7), Vec::new())),
            summary_of(region(7, 7), 0),
        );
        assert_eq!(one_base, vec![(at(7), 0)]);

        let tract = collect_reported(
            &Drawn::Built(tract_record(region(7, 9), Vec::new())),
            summary_of(region(7, 9), 0),
        );
        assert_eq!(tract, vec![(at(7), 0), (at(8), 0), (at(9), 0)]);
    }

    /// A witness that stops inside a tract raises only the positions it saw. This is the
    /// difference between observation depth along a locus and one pooled count spread flat,
    /// and it is why a tract's depth is taken from the body rather than from the summary.
    #[test]
    fn a_partial_witness_raises_only_the_positions_it_witnessed() {
        let record = tract_record(
            region(100, 105),
            vec![obs(
                ReadWitness::from_right(2, LocusLen::from_positions(6))
                    .expect("a run covering at least one position"),
                9,
            )],
        );
        let reported = collect_reported(&Drawn::Built(record), summary_of(region(100, 105), 0));
        assert_eq!(
            reported,
            vec![
                (at(100), 0),
                (at(101), 0),
                (at(102), 0),
                (at(103), 0),
                (at(104), 9),
                (at(105), 9),
            ],
        );
    }

    /// A repeat cluster covers its ground the way a tract does. No such record is emitted
    /// today, and this pins which arm it would take if one ever were.
    #[test]
    fn a_repeat_bundle_reports_every_position_of_its_span() {
        let record = SampleLocusObservations {
            kind: LocusKind::SsrBundle,
            ..generic_record(region(20, 22), vec![obs(ReadWitness::Complete, 4)])
        };
        let reported = collect_reported(&Drawn::Built(record), summary_of(region(20, 22), 0));
        assert_eq!(reported, vec![(at(20), 4), (at(21), 4), (at(22), 4)]);
    }

    /// The rule is decided by the record's span and never by which shape it arrived in, so a
    /// tract that spans one base is answered from its summary like any other one-base record.
    /// Nothing emits such a tract; the test exists because the alternative — dispatching on
    /// kind first — would make the two modes' answers depend on whether a body was in hand.
    #[test]
    fn the_span_decides_the_rule_and_not_the_records_kind() {
        let record = tract_record(region(50, 50), vec![obs(ReadWitness::Complete, 12)]);
        assert_eq!(record.num_obs_along_locus(), vec![12]);

        let reported = collect_reported(&Drawn::Built(record), summary_of(region(50, 50), 5));
        assert_eq!(reported, vec![(at(50), 5)]);
    }

    /// A summary belonging to another record is refused loudly, **in the release profile this
    /// repo runs and not only in a debug build**. Silently accepted, it would report one
    /// record's depths at another record's coordinates.
    #[test]
    #[should_panic(expected = "the summary handed in is not this record's")]
    fn a_summary_that_is_not_the_records_own_is_refused() {
        let record = tract_record(region(9_000, 9_002), vec![obs(ReadWitness::Complete, 4)]);
        let _ = collect_reported(&Drawn::Built(record), summary_of(region(10, 12), 99));
    }

    /// The same refusal on a kept draw, where the summary the draw carries is the one that is
    /// used and the parameter is the one checked against it.
    #[test]
    #[should_panic(expected = "the summary handed in is not this record's")]
    fn a_kept_draw_whose_own_summary_disagrees_with_the_parameter_is_refused() {
        let drawn = Drawn::Kept {
            summary: summary_of(region(10, 12), 99),
            body: 0..1,
        };
        let _ = collect_reported_building_with(&drawn, summary_of(region(50, 52), 99), |_| {
            Ok::<_, &'static str>(tract_record(
                region(10, 12),
                vec![obs(ReadWitness::Complete, 4)],
            ))
        });
    }

    /// A stored head that claims other ground than the body behind it is refused rather than
    /// reported at the head's coordinates — a fact about the file, in psp mode.
    #[test]
    #[should_panic(expected = "covers other ground than its head claimed")]
    fn a_built_body_that_covers_other_ground_than_its_head_claimed_is_refused() {
        let summary = summary_of(region(10, 12), 99);
        let drawn = Drawn::Kept {
            summary,
            body: 0..1,
        };
        let _ = collect_reported_building_with(&drawn, summary, |_| {
            Ok::<_, &'static str>(tract_record(
                region(9_000, 9_002),
                vec![obs(ReadWitness::Complete, 4)],
            ))
        });
    }

    /// A region whose end precedes its start describes no ground, and reports nothing rather
    /// than a position that is not there. Nothing mints such a record today — this pins the
    /// behaviour rather than leaving it to be discovered, because the two derivations involved
    /// disagree about it: the span test reads the region as the ground it names, and
    /// `GenomeRegion::len` saturates to zero.
    #[test]
    fn a_region_that_ends_before_it_starts_reports_no_position() {
        let record = tract_record(region(12, 10), vec![obs(ReadWitness::Complete, 4)]);
        assert_eq!(record.num_obs_along_locus(), Vec::<u32>::new());

        let reported = collect_reported(&Drawn::Built(record), summary_of(region(12, 10), 7));
        assert_eq!(reported, Vec::new());
    }

    proptest::proptest! {
        /// The post-conditions the doc states, over the whole small domain rather than the
        /// dozen points the fixtures above reach: the two draw shapes agree, positions
        /// strictly increase, every one lies inside the record's region, and how many there
        /// are is decided by the span and the kind. The first is the one the plan cares about
        /// most — it is what stops a later edit making psp mode and direct mode differ at a
        /// span or a witness layout nobody wrote a fixture for.
        #[test]
        fn every_reported_position_is_inside_the_record_and_the_two_draw_shapes_agree(
            span in 1u64..=8,
            witness_from in 0u64..8,
            witness_covers in 1u64..=8,
            is_tract in proptest::bool::ANY,
            reads in 0u32..5,
        ) {
            let region = region(1_000, 1_000 + span - 1);
            // A run of `witness_covers` positions from `witness_from`, both clamped into the
            // locus, so the sweep covers whole and partial witnesses alike.
            let from = witness_from.min(span - 1);
            let covers = witness_covers.min(span - from);
            let witness = ReadWitness::Partial {
                positions: WitnessedLocusPositions::one_run_from_offset_and_length(
                    from as u16,
                    covers as u16,
                )
                .expect("a run covering at least one position"),
            };
            let observed = vec![obs(witness, reads)];
            let record = || {
                if is_tract {
                    tract_record(region, observed.clone())
                } else {
                    generic_record(region, observed.clone())
                }
            };
            let summary = summary_of(region, reads);

            let built = collect_reported(&Drawn::Built(record()), summary);
            let kept = collect_reported_building_with(
                &Drawn::Kept { summary, body: 0..1 },
                summary,
                |_| Ok::<_, &'static str>(record()),
            );
            proptest::prop_assert_eq!(&built, &kept);

            let positions_expected = if span == 1 || !is_tract { 1 } else { span as usize };
            proptest::prop_assert_eq!(built.len(), positions_expected);
            for pair in built.windows(2) {
                proptest::prop_assert!(pair[0].0 < pair[1].0);
            }
            for (position, _) in &built {
                proptest::prop_assert_eq!(position.contig, region.contig);
                proptest::prop_assert!(position.position >= region.start);
                proptest::prop_assert!(position.position <= region.end);
            }
        }
    }
}
