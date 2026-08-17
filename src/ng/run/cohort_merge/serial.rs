//! The whole analysed stretch, merged by one builder on one thread — **the oracle**.
//!
//! This is the simplest thing that produces the right answer: one builder, no cache, no
//! organiser, no threads, holding every sample's observations and walking the analysed
//! regions in order. Everything the later milestones add is about speed and memory, so
//! anything that changes this output is a defect rather than a trade
//! (`doc/devel/ng/spec/cohort_merge.md` §15; the plan's C2).
//!
//! **It is deliberately thin.** With the whole stretch in hand, merging is
//! [`build_region`](super::build::build_region) over each analysed region and nothing more:
//! the loci a builder closes depend on the observations it can see, and here it can see all
//! of them. What the parallel arrangement adds is a *narrower view* per builder — a window
//! drawn from the cache (milestone D) and regions handed out in parallel (E) — and the loci
//! that differ between the two are exactly what the organiser's overlap resolution exists to
//! settle (spec §6.1). That is why this stands as the oracle: it has no window and no
//! resolution to get wrong.

use super::build::{RegionOutcome, build_region};
use super::{MaxCohortLocusSpan, MinAltObs};
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::types::GenomeRegion;

/// Merge the cohort over `analysed`, one region at a time, in the order given.
///
/// `observations_per_sample` is one slice per sample in the run's sample order, each in
/// coordinate order, holding **everything over the analysed stretch** — this driver hands the
/// same whole slices to every region rather than windowing them, which is the point of it.
///
/// The analysed regions are the run's own (a contig, a `--regions` interval); they are walked
/// in the order given and each locus is owned by the one its first position falls in, so a
/// locus that reaches from one into the next is built once, whole, by the first
/// ([`build_region`]). Where the analysed regions are disjoint and ascending — which is what
/// a run supplies — the observations come out in genome order and the failed spans with them.
///
/// **No overlap resolution happens here, and none is needed.** Two builders resolve overlaps
/// because each saw a different window; here every region is built from the same complete
/// view, so the loci closed over one region are the loci closed over any other, and the
/// ownership rule alone keeps each exactly once (spec §6.1, §9).
///
/// **The analysed regions must be disjoint and ascending, and that is checked rather than
/// assumed.** Two regions covering the same ground would each own the loci opening in the
/// overlap and the run would carry them twice — not out of order, but *twice*, which no
/// consumer can tell from a cohort that really varied there. A run's own regions are
/// normalised before they reach here, but this signature takes a bare slice, and a
/// user-supplied BED is exactly the shape that arrives overlapping. Release-level, at one
/// comparison per region against a walk that closes every locus in the genome.
///
/// **The cost of dividing the stretch finely is real, and this driver does not pay it for
/// you.** Every call closes the loci from the beginning of the observations it is given and
/// discards those before its own ground ([`build_region`]), so the same 20,000 observations
/// cost 5.4 ms in one region and 184 ms in a thousand — 34 times the work for the same
/// answer, measured in a release build by the C2 review. **Hand it the run's own analysed
/// regions**, which are few and long; short building regions are what the observation cache
/// exists to make affordable (milestone D), and that is where the parallel arrangement will
/// get them.
pub fn merge_cohort_serially(
    analysed: &[GenomeRegion],
    observations_per_sample: &[&[SampleLocusObservations]],
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_obs: MinAltObs,
) -> RegionOutcome {
    let mut merged = RegionOutcome::default();

    for pair in analysed.windows(2) {
        let (earlier, later) = (pair[0], pair[1]);
        assert!(
            (earlier.contig, earlier.end) < (later.contig, later.start),
            "the analysed regions {earlier} and {later} are not disjoint and ascending, so \
             the loci opening in the ground they share would be built — and carried — twice",
        );
    }

    for region in analysed {
        let outcome = build_region(
            *region,
            observations_per_sample,
            max_cohort_locus_span,
            min_alt_obs,
        );
        merged
            .cohort_observations
            .extend(outcome.cohort_observations);
        merged.failed_locus_spans.extend(outcome.failed_locus_spans);
    }

    merged
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::types::{ContigId, Position, ReadGroupId};

    fn region_on(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        region_on(0, start, end)
    }

    /// One sample's record over `region`, showing `observed_bases` to three reads.
    fn member(
        region: GenomeRegion,
        reference_bases: &[u8],
        observed_bases: &[u8],
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region,
            reference_bases: Box::from(reference_bases),
            observations: vec![SequenceObservation {
                bases: Box::from(observed_bases),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 3,
                num_fwd: 3,
                q_sum: -6.0,
                mapq_sum: 180,
                mapq_sum_sq: 10_800,
                placed_left: 0,
                chain_ids: vec![1, 2, 3],
            }],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// The observations come out in genome order across several analysed regions, and each
    /// locus is built exactly once — by the region its first position falls in.
    #[test]
    fn the_analysed_regions_are_walked_in_order_and_each_locus_built_once() {
        let sample = [
            member(region(12, 12), b"G", b"T"),
            member(region(48, 52), b"ACGTA", b"A"),
            member(region(77, 77), b"G", b"T"),
        ];
        let analysed = [region(1, 50), region(51, 100)];

        let merged = merge_cohort_serially(
            &analysed,
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(12, 12), region(48, 52), region(77, 77)],
            "the deletion opening at 48 belongs to the first region and reaches into the \
             second, which does not build it again",
        );
    }

    /// A run over several contigs comes out contig by contig, and a locus is never claimed by
    /// another contig's region however its positions compare.
    #[test]
    fn several_contigs_come_out_in_the_order_they_were_analysed() {
        let first = [member(region_on(0, 20, 20), b"G", b"T")];
        let second = [member(region_on(1, 20, 20), b"G", b"C")];
        let analysed = [region_on(0, 1, 100), region_on(1, 1, 100)];

        let merged = merge_cohort_serially(
            &analysed,
            &[&first, &second],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region_on(0, 20, 20), region_on(1, 20, 20)],
        );
        assert_eq!(
            merged.cohort_observations[0].per_sample.len(),
            1,
            "and each locus holds only the sample that covered it",
        );
    }

    /// The failed spans travel with the observations, region by region, and a locus that was
    /// merely too quiet appears in neither — the distinction the run summary rests on.
    ///
    /// **One failed locus in each of the two regions**, deliberately: with only the first
    /// region's, a driver that stopped after one region or walked them backwards would pass
    /// this test unchanged.
    #[test]
    fn the_failed_spans_are_gathered_across_regions_and_the_quiet_are_not() {
        let wide = [
            member(region(20, 40), &[b'A'; 21], b"A"),
            member(region(60, 85), &[b'A'; 26], b"A"),
        ];
        let quiet = [SampleLocusObservations {
            observations: vec![SequenceObservation {
                num_obs: 1,
                ..member(region(95, 95), b"A", b"C").observations.remove(0)
            }],
            ..member(region(95, 95), b"A", b"C")
        }];
        let analysed = [region(1, 50), region(51, 100)];

        let merged = merge_cohort_serially(
            &analysed,
            &[&wide, &quiet],
            MaxCohortLocusSpan(std::num::NonZeroU32::new(20).expect("20 is non-zero")),
            MinAltObs::DEFAULT,
        );

        assert!(merged.cohort_observations.is_empty());
        assert_eq!(
            merged.failed_locus_spans,
            vec![region(20, 40), region(60, 85)],
            "21 and 26 bases against a bound of 20 are refused and counted, one from each \
             region and in genome order; the single read at 95 reaches neither vector",
        );
    }

    // ---------------------------------------------------------------
    // The milestone's own claim: cohort observations from observations the generator
    // actually minted, over reads on disk — one builder, one thread (Checkpoint C).
    // ---------------------------------------------------------------

    /// Reads prepared as they arrive: no left-alignment, no BAQ. What this fixture is about
    /// is the merge, and the generator's own tests already pin what preparation does.
    ///
    /// **It matters that this is stated rather than defaulted.** Unification by byte equality
    /// is sound only because indels are left-aligned upstream (spec §4.2), and these reads
    /// carry their deletion at one placement in every read, so the fixture does not depend on
    /// what this preparer skips.
    struct ReadsAsTheyArrive;

    impl crate::ng::read::ReadPreparer for ReadsAsTheyArrive {
        type Scratch = ();

        fn prepare_read(
            &self,
            read: crate::ng::read::AlignedRead,
            _scratch: &mut Self::Scratch,
        ) -> Result<Option<crate::ng::read::PreparedRead>, crate::ng::read::ReadPrepError> {
            let read_group = read.read_group;
            let chrom_id = u32::try_from(read.ref_id).expect("a fixture contig id fits u32");
            Ok(Some(crate::ng::read::PreparedRead::from_production(
                crate::pileup::per_sample::baq_engine::prepare_passthrough(
                    read.into_mapped_read(),
                    chrom_id,
                ),
                read_group,
            )))
        }
    }

    /// Mint one sample's observations over `region` from `records` on disk, through the real
    /// generic locus generator.
    ///
    /// The two temporary directories come back with the loci because they hold the BAM and
    /// the reference: dropped early, the sample would be reading files that no longer exist.
    fn minted_over(
        records: &[noodles_sam::alignment::RecordBuf],
        region: GenomeRegion,
    ) -> (
        tempfile::TempDir,
        tempfile::TempDir,
        Vec<SampleLocusObservations>,
    ) {
        use crate::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
        use crate::ng::read::filtering::ReadFilterConfig;
        use crate::ng::read::input::SampleReads;
        use crate::ng::read::input::test_fixtures::{
            fixture_reference, fixture_reference_bases, header, indexed_bam, matching_contigs,
        };
        use std::sync::Arc;

        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, bam) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            records,
        );
        let reads =
            SampleReads::open_only_sample(&[bam], &reference, ReadFilterConfig::default(), false)
                .expect("the fixture sample opens");

        let mut generator = PileupGenerator::new(
            Arc::new(fixture_reference_bases()),
            fixture_reference_bases,
            ReadsAsTheyArrive,
            PileupGeneratorConfig::default(),
        )
        .expect("the generator's configuration is the default one");

        generator.begin_segment(region);
        let mut loci = Vec::new();
        while let Some(locus) = generator
            .next_locus(&reads)
            .expect("the walk over a fixture BAM succeeds")
        {
            loci.push(locus);
        }
        (reference_dir, bam_dir, loci)
    }

    /// A 30-base read at `start` on chr2, all reference except one base.
    fn read_with_a_substitution(
        name: &str,
        start: usize,
        offset_into_read: usize,
    ) -> noodles_sam::alignment::RecordBuf {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        use noodles_sam::alignment::record_buf::Sequence;

        let mut record = read_named_with_length(name, 1, start, 30);
        let mut bases = vec![b'A'; 30];
        bases[offset_into_read] = b'C';
        *record.sequence_mut() = Sequence::from(bases);
        record
    }

    /// A read that matches for `before` bases, deletes five, and matches for `after` more.
    fn read_with_a_five_base_deletion(
        name: &str,
        start: usize,
        before: usize,
        after: usize,
    ) -> noodles_sam::alignment::RecordBuf {
        use crate::ng::read::input::test_fixtures::read_named_with_length;
        use noodles_sam::alignment::record::cigar::Op;
        use noodles_sam::alignment::record::cigar::op::Kind;
        use noodles_sam::alignment::record_buf::{QualityScores, Sequence};

        let mut record = read_named_with_length(name, 1, start, before + after);
        *record.cigar_mut() = [
            Op::new(Kind::Match, before),
            Op::new(Kind::Deletion, 5),
            Op::new(Kind::Match, after),
        ]
        .into_iter()
        .collect();
        *record.sequence_mut() = Sequence::from(vec![b'A'; before + after]);
        *record.quality_scores_mut() = QualityScores::from(vec![30u8; before + after]);
        record
    }

    /// **A cohort observation built from what the generator actually minted** — the
    /// milestone's own claim, and the first time this module meets records it did not
    /// fabricate.
    ///
    /// Two samples over chr2 of the fixture reference, which is all `A`. Both samples' reads
    /// start at position 95:
    ///
    /// - **sample 0** carries a substitution at position **112** — the eighteenth base of a
    ///   thirty-base read — and nothing else;
    /// - **sample 1** carries a five-base deletion over **110–114**, which the mint records
    ///   as one record spanning 109–114: the anchor base and the five it deleted.
    ///
    /// The deletion covers the substitution's position, so the two samples' records chain
    /// into **one** cohort locus — the case the whole of milestone B is built around, here
    /// on records minted from reads on disk rather than written by hand.
    ///
    /// **One of sample 0's reads carries a bad base at 110**, quality 5 against 30
    /// everywhere else. Without it every one of that sample's six records inside the locus
    /// would carry the same quality sum, and "the weakest of the six" would equal the
    /// strongest — an assertion any rule satisfies. With it, taking the best sighting
    /// instead of the worst gives a different answer here as well as in `build.rs`.
    #[test]
    fn a_cohort_observation_is_built_from_minted_observations() {
        let analysed = region_on(1, 80, 140);

        let substituted: Vec<_> = (0..3)
            .map(|read| {
                let mut record = read_with_a_substitution(&format!("sub{read}"), 95, 17);
                if read == 0 {
                    // Position 110 is the sixteenth base of a read starting at 95, and this
                    // is what makes the six records' qualities differ — see the doc above.
                    let mut qualities = vec![30u8; 30];
                    qualities[15] = 5;
                    *record.quality_scores_mut() =
                        noodles_sam::alignment::record_buf::QualityScores::from(qualities);
                }
                record
            })
            .collect();
        let deleted: Vec<_> = (0..3)
            .map(|read| read_with_a_five_base_deletion(&format!("del{read}"), 95, 15, 20))
            .collect();

        let (_substituted_reference, _substituted_bam, first) = minted_over(&substituted, analysed);
        let (_deleted_reference, _deleted_bam, second) = minted_over(&deleted, analysed);

        assert!(
            !first.is_empty() && !second.is_empty(),
            "the generator minted nothing: {} / {}",
            first.len(),
            second.len(),
        );

        let merged = merge_cohort_serially(
            &[analysed],
            &[&first, &second],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert_eq!(merged.cohort_observations.len(), 1, "one locus, not two");
        let observed = &merged.cohort_observations[0];

        assert_eq!(
            observed.region,
            region_on(1, 109, 114),
            "the deletion's own record, which the substitution chained into",
        );
        assert_eq!(
            observed
                .alleles
                .iter()
                .map(|allele| String::from_utf8_lossy(allele).into_owned())
                .collect::<Vec<_>>(),
            vec!["AAAAAA", "AAACAA", "A"],
            "the reference over six bases, the substitution widened onto them, and the \
             deletion — which the reference is all `A`, so the substituted base is the only \
             letter in it",
        );

        // **The substitution sample is the multi-record case, on minted records.** The
        // generic mint writes a record at every covered position, so inside this six-base
        // locus that sample has six of them, and each of its three reads is named at all six
        // — so its allele is composed across them and its quality sums are divided.
        let substituting = &observed.per_sample[0];
        assert_eq!(substituting.sample, 0);
        assert_eq!(substituting.reads_composed_across_records, 3);
        assert_eq!(substituting.reads_removed_as_evidence, 0);
        assert_eq!(
            substituting.support_for(1).num_reads,
            3,
            "three reads showed the substitution across the whole locus",
        );

        let weakest_of_its_records = first
            .iter()
            .filter(|record| {
                record.region.start >= observed.region.start
                    && record.region.end <= observed.region.end
            })
            .map(|record| {
                let sequence = &record.observations[0];
                sequence.q_sum / f64::from(sequence.num_obs)
            })
            .fold(f64::NEG_INFINITY, f64::max);
        assert_eq!(
            substituting.support_for(1).q_sum,
            weakest_of_its_records * 3.0,
            "each read takes the weakest of the six positions it was seen at",
        );

        // **The deletion sample is the one-record case**, so its numbers are the mint's own,
        // undivided — asserted against the record itself rather than against a constant.
        let deleting = &observed.per_sample[1];
        let minted = &second
            .iter()
            .find(|record| record.region == observed.region)
            .expect("the deletion's record is the locus")
            .observations[0];
        assert_eq!(deleting.sample, 1);
        assert_eq!(deleting.reads_composed_across_records, 0);
        assert_eq!(
            deleting.support_for(2),
            crate::ng::run::cohort_merge::build::AlleleSupport {
                num_reads: minted.num_obs,
                num_fwd: minted.num_fwd,
                q_sum: minted.q_sum,
                mapq_sum: minted.mapq_sum,
                mapq_sum_sq: minted.mapq_sum_sq,
                placed_left: minted.placed_left,
            },
            "one record, so every sum is the one the generator wrote",
        );

        assert_eq!(
            substituting.support_for(0),
            crate::ng::run::cohort_merge::build::AlleleSupport::default(),
            "neither sample's reads showed the reference over the whole locus",
        );
        assert_eq!(
            deleting.support_for(0),
            crate::ng::run::cohort_merge::build::AlleleSupport::default(),
        );
    }

    /// **The same observations give the same output however the analysed ground is divided**
    /// — the property spec §15 calls this component's regression anchor, and the reason this
    /// driver is the oracle at all.
    ///
    /// Sixty loci over six hundred bases, built as one region, as six, as sixty and as a
    /// hundred and twenty. A rule that lost a locus on a boundary, claimed one twice, or
    /// reordered them would show here at one width and not another.
    #[test]
    fn the_same_loci_come_out_however_the_analysed_ground_is_divided() {
        let sample: Vec<_> = (0..60)
            .map(|locus| {
                let at = 10 * locus + 1;
                member(region(at, at), b"G", b"T")
            })
            .collect();

        let built_in = |regions: u64| {
            let width = 600 / regions;
            let analysed: Vec<_> = (0..regions)
                .map(|piece| region(piece * width + 1, (piece + 1) * width))
                .collect();
            merge_cohort_serially(
                &analysed,
                &[&sample],
                MaxCohortLocusSpan::DEFAULT,
                MinAltObs::DEFAULT,
            )
            .cohort_observations
            .iter()
            .map(|observed| observed.region)
            .collect::<Vec<_>>()
        };

        let whole = built_in(1);
        assert_eq!(whole.len(), 60, "every locus is built once in one region");
        for regions in [6, 60, 120] {
            assert_eq!(
                built_in(regions),
                whole,
                "dividing the ground into {regions} regions changed the answer",
            );
        }
    }

    /// **Analysed regions that overlap are refused**, rather than building the loci in the
    /// ground they share twice. The duplicate is the dangerous outcome: nothing downstream
    /// can tell one locus carried twice from a cohort that really varied at two places.
    #[test]
    #[should_panic(expected = "not disjoint and ascending")]
    fn analysed_regions_that_overlap_are_refused() {
        let sample = [member(region(45, 45), b"G", b"T")];

        let _ = merge_cohort_serially(
            &[region(1, 60), region(40, 100)],
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );
    }

    /// And regions in descending order are refused by the same check: they would come out in
    /// the order they were given, which is not genome order, and every consumer of this
    /// stream reads it as genome order.
    #[test]
    #[should_panic(expected = "not disjoint and ascending")]
    fn analysed_regions_out_of_order_are_refused() {
        let sample = [member(region(45, 45), b"G", b"T")];

        let _ = merge_cohort_serially(
            &[region(61, 100), region(1, 60)],
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );
    }

    /// **The keep threshold a run sets is the one the builders use.** Every other test here
    /// runs at the default, so a driver that dropped the argument would ship green.
    #[test]
    fn the_keep_threshold_reaches_the_builders() {
        let sample = [member(region(12, 12), b"G", b"T")];
        let analysed = [region(1, 50)];

        let built = merge_cohort_serially(
            &analysed,
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );
        let too_quiet = merge_cohort_serially(
            &analysed,
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs(std::num::NonZeroU32::new(4).expect("4 is non-zero")),
        );

        assert_eq!(built.cohort_observations.len(), 1, "three reads reach two");
        assert!(too_quiet.cohort_observations.is_empty(), "and not four");
    }

    /// Nothing analysed is not an error: an empty run yields an empty outcome rather than a
    /// failure, and the caller's own emptiness check is where a zero-sample run is refused
    /// (spec §7.2).
    #[test]
    fn analysing_nothing_yields_nothing() {
        let sample = [member(region(12, 12), b"G", b"T")];

        let merged = merge_cohort_serially(
            &[],
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltObs::DEFAULT,
        );

        assert!(merged.cohort_observations.is_empty() && merged.failed_locus_spans.is_empty());
    }
}
