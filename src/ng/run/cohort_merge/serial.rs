//! Two serial drivers of the same merge: the oracle, and the same merge read through the cache.
//!
//! **Vocabulary, because the whole file turns on it.** An *analysed region* is the run's own
//! interval — a contig, a `--regions` interval — and they are few and long. A *building region*
//! is the short stretch one builder owns, `cohort_locus_builder_regions_len` bases, 20 by
//! default (spec §6.1); building regions divide an analysed region.
//!
//! [`merge_cohort_serially`] is **the oracle**: the simplest thing that produces the right
//! answer — one builder, no cache, no organiser, no threads, holding every sample's
//! observations and walking the analysed regions in order. Everything the later milestones add
//! is about speed and memory, so anything that changes this output is a defect rather than a
//! trade (`doc/devel/ng/spec/cohort_merge.md` §15; the plan's C2).
//!
//! [`merge_cohort_through_cache`] does the same job through **one forward reader per sample**,
//! and its whole claim is that its output is the oracle's byte for byte (the plan's D2). It is
//! still one builder on one thread; what changes is the *view* — a window over each builder's
//! own ground instead of the whole stretch, which is what makes short building regions
//! affordable. The builders run at the same time as each other in [`parallel`](super::parallel),
//! over this same cache; the loci that differ once builders see different windows are what the
//! organiser's overlap resolution settles (spec §6.1).

use super::build::{
    CohortObservation, RegionOutcome, build_region, build_region_handing_over_windowed,
};
use crate::ng::locus_generation::SampleLocusObservations;
use super::observation_cache::{ObservationCache, ObservationSource, building_regions_of};
use super::timing;
use super::{
    CohortLocusBuilderRegionsLen, MaxCohortLocusSpan, MinAltReads,
    refuse_malformed_analysed_regions,
};
use crate::ng::types::{GenomePosition, GenomeRegion};

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
pub fn merge_cohort_serially<'a>(
    analysed: &[GenomeRegion],
    observations_per_sample: &[&'a [SampleLocusObservations]],
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
) -> RegionOutcome {
    let mut merged = RegionOutcome::default();
    refuse_malformed_analysed_regions(analysed);

    for region in analysed {
        let outcome = build_region(
            *region,
            observations_per_sample,
            max_cohort_locus_span,
            min_alt_reads,
        );
        merged
            .cohort_observations
            .extend(outcome.cohort_observations);
        merged.failed_locus_spans.extend(outcome.failed_locus_spans);
    }

    merged
}

/// Merge the cohort over `analysed`, reading each sample through the observation cache and
/// building in regions of `cohort_locus_builder_regions_len` bases.
///
/// **The output is [`merge_cohort_serially`]'s, and that is the whole claim of this step**
/// (the plan's D2). The oracle holds every sample's observations at once and hands the whole
/// stretch to every builder; this holds one forward reader per sample and hands each builder a
/// window over its own ground. Nothing else differs, and nothing else may: the cache exists to
/// make the merge affordable, not to change what it answers
/// (`doc/devel/ng/spec/cohort_merge.md` §15).
///
/// **What it buys is that the analysed ground can now be divided finely.** Handing a builder
/// the whole stretch costs it every locus that opened before its own ground, closed in full and
/// then discarded — about 3.3 µs per prefix base at 63 samples (the C1 review), which is why
/// the oracle must be given the run's own long analysed regions and this need not. Here the
/// building regions are as short as the caller asks for, which is what the parallel arrangement
/// (milestone E) will hand out. **What it costs in exchange is one cover per building region**,
/// `sweeps × (samples + held)` work each — 616 µs at 1,000 samples and 2.87 ms at 3,000 on
/// 20-base regions, measured by the D1 review; [`ObservationCache::cover`] carries the table.
///
/// **The cache comes in from the caller rather than being made here**, because it outlives one
/// merge: it holds one reader per sample *for the whole run* (spec §6.4), and at milestone E it
/// belongs to the organiser. It also makes what this driver does to memory visible — a test can
/// ask the cache afterwards how much it still holds
/// ([`held_observations_len`](ObservationCache::held_observations_len)), which is the only way to tell a
/// driver that evicts from one that does not: their outputs are identical.
///
/// Its sources are one per sample **in the run's sample order**, each yielding that sample's
/// observations in coordinate order — the shape `run_streaming.md` arch §2's `ObservationSource`
/// will have. The analysed regions must be disjoint and ascending, checked as in the oracle.
///
/// **The building regions tile each analysed region exactly, and the last one is clamped to
/// it.** A locus belongs to the builder whose region holds its first position, so a building
/// region running past the analysed ground would claim loci the oracle never builds — the
/// division has to end where the run's own interval ends.
///
/// **Eviction happens before each cover, at that region's first base.** What ends before it can
/// reach no locus a later builder can own (see [`ObservationCache::evict_before`]), and holding
/// it would make this driver's memory the whole stretch again. Evicting *after* the cover
/// instead would give the same output and the same final window — what it costs is the sweep,
/// which re-reads the whole held window: measured on a record every ten bases at 20-base
/// regions, the pre-cover eviction takes the window from three records to one at every region.
///
/// **A failure ends the merge and leaves the cache where it stopped.** The observations built so
/// far are dropped — the caller gets the source's error and nothing else — and the cache keeps
/// the readers' position and the window the last cover drew. So **the same cache cannot be used
/// to try the same ground again**: the ground behind the failure has already been evicted, and a
/// second merge over it comes back short and says `Ok`. A run that means to retry builds a new
/// cache over new sources; a run that does not, abandons the stretch. Making that unrepresentable
/// belongs with the organiser (milestone E), which owns the cache.
///
/// **Not every source-side failure arrives as `Err`.** A source whose observations go backwards
/// trips the cache's coordinate-order check and panics ([`ObservationCache::cover`]), because the
/// cache has no error of its own to mint — `E` is the source's type and this driver only passes
/// it through. That is right while observations come from this crate's generator, where going
/// backwards is a bug; it stops being right when they are decoded from a psp file, and
/// `organise.rs` records that it owes the change to `RunError`.
pub fn merge_cohort_through_cache<S, E>(
    analysed: &[GenomeRegion],
    cache: &mut ObservationCache<S>,
    cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen,
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
) -> Result<RegionOutcome, E>
where
    S: ObservationSource<Error = E>,
{
    let mut merged = RegionOutcome::default();
    // Two disjoint borrows of one outcome, so that collecting the loci is this driver's
    // streaming form with `Vec::push` for a sink rather than a second copy of it.
    let RegionOutcome {
        cohort_observations,
        failed_locus_spans,
    } = &mut merged;
    merge_cohort_handing_each_locus_over(
        analysed,
        cache,
        cohort_locus_builder_regions_len,
        max_cohort_locus_span,
        min_alt_reads,
        &mut |built| cohort_observations.push(built),
        failed_locus_spans,
    )?;
    Ok(merged)
}

/// Merge the cohort over `analysed` as [`merge_cohort_through_cache`] does, **handing each
/// surviving locus to `keep` where it is built** instead of collecting them.
///
/// The two are one driver: `merge_cohort_through_cache` is this with `Vec::push` for a sink,
/// so nothing about the division of the ground, the eviction or the covering differs between
/// them, and the merge's oracles check both at once.
///
/// **The parallel-cover form below is what a calling run drives since E1; this serial form is
/// its oracle** — one driver has to keep the schedule the fixtures were reasoned about under.
/// What handing-over buys either way is that the buffer holds what the sink made of each locus
/// rather than the locus itself: a called locus is one genotype per sample, where a cohort
/// observation is every covering sample's reads folded onto the locus's alleles, and
/// collecting the observations for a whole run is what spec §5.1's bound forbids.
///
/// **The failure behaviour is `merge_cohort_through_cache`'s, and what the sink has already
/// been handed is kept**: a caller that owns the sink's buffer still holds every locus built
/// before the failing cover. The collecting form drops them, because its buffer is the value
/// it does not return.
pub fn merge_cohort_handing_each_locus_over<S, E>(
    analysed: &[GenomeRegion],
    cache: &mut ObservationCache<S>,
    cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen,
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
    keep: &mut impl FnMut(CohortObservation),
    refused: &mut Vec<GenomeRegion>,
) -> Result<(), E>
where
    S: ObservationSource<Error = E>,
{
    merge_handing_each_locus_over_with(
        analysed,
        cache,
        cohort_locus_builder_regions_len,
        max_cohort_locus_span,
        min_alt_reads,
        keep,
        refused,
        &mut ObservationCache::cover,
    )
}

/// [`merge_cohort_handing_each_locus_over`], with each cover's samples drawn forward
/// **concurrently** instead of one after another.
///
/// **This is where a calling run's parallelism went, and the measurement that put it here is
/// the plan's own.** On 63 tomato accessions, drawing the readers forward — every sample's
/// walk, which is the reads being decoded — is **88% of `call_cohort` and runs on one
/// thread**, while assembling the loci and genotyping them together are 11%
/// (`doc/devel/reports/implementations/ng_run_driver_e1_2026-09-01.md`, which extends D3's
/// 3–24-sample table to the whole cohort). A pool that genotypes several loci at once therefore cannot buy a
/// run more than a few percent; sweeping the cohort's samples concurrently inside each cover
/// is the arrangement that reaches the 88%.
///
/// **Nothing about the answer changes, and that is the cover's own guarantee, not this
/// driver's.** [`ObservationCache::cover_in_parallel`] reaches the same fixpoint as the
/// serial sweep by a different schedule and leaves the same held window (its documentation
/// carries the argument; the parallel merge's whole oracle battery rests on it). Everything
/// downstream of the cover — eviction, the builders, the sink — runs on the calling thread
/// exactly as in the serial form, in the same order, so the loci handed to `keep` are
/// byte-identical at every thread count. `the_parallel_cover_gives_the_serial_drivers_answer`
/// pins it here; the run's own concurrency-invariance oracle (spec §12.2, the plan's E2) is
/// what will pin it end to end, and is not built yet.
///
/// **On a pool of one thread it takes the serial sweep**, the same way the parallel merge
/// asks `rayon::current_num_threads` before handing eviction to workers: the Jacobi schedule
/// buys nothing with nobody to share the sweep, and can cost one extra sweep per chain link.
///
/// **One failure shape is less determined than the serial form's**: when two samples' sources
/// fail during one sweep, which error comes back depends on the schedule — the serial sweep
/// always reports the first in the run's sample order. The parallel merge has the same
/// property. A run stops either way, naming a sample that really failed.
pub fn merge_cohort_handing_each_locus_over_covering_samples_in_parallel<S, E>(
    analysed: &[GenomeRegion],
    cache: &mut ObservationCache<S>,
    cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen,
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
    keep: &mut impl FnMut(CohortObservation),
    refused: &mut Vec<GenomeRegion>,
) -> Result<(), E>
where
    S: ObservationSource<Error = E> + Send,
    E: Send,
{
    if rayon::current_num_threads() > 1 {
        merge_handing_each_locus_over_with(
            analysed,
            cache,
            cohort_locus_builder_regions_len,
            max_cohort_locus_span,
            min_alt_reads,
            keep,
            refused,
            &mut ObservationCache::cover_in_parallel,
        )
    } else {
        merge_cohort_handing_each_locus_over(
            analysed,
            cache,
            cohort_locus_builder_regions_len,
            max_cohort_locus_span,
            min_alt_reads,
            keep,
            refused,
        )
    }
}

/// The one body behind both handing-over drivers: everything but how a cover sweeps.
///
/// `cover` is [`ObservationCache::cover`] or [`ObservationCache::cover_in_parallel`] and
/// nothing else — the division of the ground, the eviction, the building and the sink are
/// written once here, so the two public forms cannot drift on anything but the sweep
/// schedule, which is the one thing the cover's own fixpoint argument covers.
#[expect(
    clippy::too_many_arguments,
    reason = "the seven are the two public drivers' shared signature plus the cover; \
              grouping them would make a struct whose only purpose is this private call"
)]
fn merge_handing_each_locus_over_with<S, E>(
    analysed: &[GenomeRegion],
    cache: &mut ObservationCache<S>,
    cohort_locus_builder_regions_len: CohortLocusBuilderRegionsLen,
    max_cohort_locus_span: MaxCohortLocusSpan,
    min_alt_reads: MinAltReads,
    keep: &mut impl FnMut(CohortObservation),
    refused: &mut Vec<GenomeRegion>,
    cover: &mut impl FnMut(&mut ObservationCache<S>, GenomeRegion) -> Result<(), E>,
) -> Result<(), E>
where
    S: ObservationSource<Error = E>,
{
    // Zero-sized and doing nothing without `--features merge-timing` (`super::timing`), and
    // named the same parts as the parallel driver's so the two breakdowns can be compared
    // line for line. A region is this driver's round.
    let whole_merge = timing::Stopwatch::start();
    refuse_malformed_analysed_regions(analysed);
    // The sink is timed apart from the assembling it happens inside, because in a calling run
    // it *is* the genotyping and the two together are one undifferentiated builder time.
    let mut timed_keep = |built| {
        let after_assembly = timing::Stopwatch::start();
        keep(built);
        after_assembly.add_to(&timing::AFTER_ASSEMBLY_NANOS);
    };

    for analysed_region in analysed {
        for building_region in
            building_regions_of(*analysed_region, cohort_locus_builder_regions_len)
        {
            let evicting = timing::Stopwatch::start();
            cache.evict_before(GenomePosition {
                contig: building_region.contig,
                position: building_region.start,
            });
            evicting.add_to(&timing::EVICT_NANOS);
            let covering = timing::Stopwatch::start();
            cover(cache, building_region)?;
            covering.add_to(&timing::COVER_NANOS);
            timing::ROUNDS.add(1);
            timing::REGIONS.add(1);
            let builder = timing::Stopwatch::start();
            cache.with_observations(building_region, |window| {
                build_region_handing_over_windowed(
                    building_region,
                    window,
                    max_cohort_locus_span,
                    min_alt_reads,
                    &mut timed_keep,
                    refused,
                );
            });
            // **The sink's own work is inside this**, where the collecting form's push was
            // too: a run that calls each locus here charges the genotyping to the builder,
            // which is what the milestone-E measurement has to be able to see.
            let busy = builder.elapsed_nanos();
            timing::BUILDER_BUSY_NANOS.add(busy);
            timing::SLOWEST_BUILDER_NANOS.add(busy);
            timing::ROUND_WALL_NANOS.add(busy);
        }
    }

    whole_merge.add_to(&timing::MERGE_WALL_NANOS);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::super::build::build_region_windowed;
    use crate::ng::locus_generation::SampleLocusObservations;
    use super::super::fixtures::{
        SourceFailed, in_flight, member, refuse_any_difference, region, region_on, render,
        source_of, three_samples_over_six_hundred_bases, width,
    };
    use super::super::organise::{Organiser, RegionIndex};
    use super::super::parallel::merge_cohort_in_parallel;
    use super::super::{MinAltObs, MinAltReadShare};
    use super::*;
    use crate::ng::locus_generation::{LocusKind, SequenceObservation};
    use crate::ng::types::{ContigId, Position};

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
            MinAltReads::DEFAULT,
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
            MinAltReads::DEFAULT,
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
            MinAltReads::DEFAULT,
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
            MinAltReads::DEFAULT,
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
            substituting.pooled_support_for(1).num_reads,
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
                sequence.q_sum.nats() / f64::from(sequence.num_obs)
            })
            .fold(f64::NEG_INFINITY, f64::max);
        assert_eq!(
            substituting.pooled_support_for(1).q_sum,
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
            deleting.pooled_support_for(2),
            crate::ng::run::cohort_merge::build::AlleleSupport {
                num_reads: minted.num_obs,
                num_fwd: minted.num_fwd,
                q_sum: minted.q_sum.nats(),
                mapq_sum: minted.mapq_sum,
                mapq_sum_sq: minted.mapq_sum_sq,
                placed_left: minted.placed_left,
            },
            "one record, so every sum is the one the generator wrote",
        );

        assert_eq!(
            substituting.pooled_support_for(0),
            crate::ng::run::cohort_merge::build::AlleleSupport::default(),
            "neither sample's reads showed the reference over the whole locus",
        );
        assert_eq!(
            deleting.pooled_support_for(0),
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
                MinAltReads::DEFAULT,
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
            MinAltReads::DEFAULT,
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
            MinAltReads::DEFAULT,
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
            MinAltReads::DEFAULT,
        );
        let too_quiet = merge_cohort_serially(
            &analysed,
            &[&sample],
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads {
                floor: MinAltObs(std::num::NonZeroU32::new(4).expect("4 is non-zero")),
                share: MinAltReadShare::DEFAULT,
            },
        );

        assert_eq!(built.cohort_observations.len(), 1, "three reads reach two");
        assert!(too_quiet.cohort_observations.is_empty(), "and not four");
    }

    /// **The two drivers agree on layouts no fixture enumerates.** The five fixtures above are
    /// hand-written; this walks 200 random ones — 1 to 6 samples, 1 or 2 contigs, records 1 to
    /// 60 bases with one in ten wide, several analysed regions with gaps between them, building
    /// widths 1 to 30, span bounds 5 to 64 and keep thresholds 1 to 4 — and compares the whole
    /// outcome of each.
    ///
    /// Written and run by the D2 review at 600 layouts with **no disagreement**, and it is not
    /// decoration: it killed two of that review's mutations. Its standing value is the regime
    /// the fixtures do not reach — a record that opens inside an analysed region and ends past
    /// its last base, which is where the cache must draw past the analysed ground and the oracle
    /// simply already has it. The counters keep the test from passing vacuously if the generator
    /// is ever narrowed.
    #[test]
    fn the_two_drivers_agree_on_random_layouts() {
        /// A seeded linear congruential generator — no dependency, and the seed is in the
        /// failure message.
        struct Seeded(u64);
        impl Seeded {
            fn next(&mut self, bound: u64) -> u64 {
                self.0 = self
                    .0
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                (self.0 >> 33) % bound
            }
        }

        let ground_end = 400u64;
        let (mut with_observations, mut with_failed, mut with_a_straddling_record) = (0, 0, 0);

        for seed in 0..200u64 {
            let mut draw = Seeded(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0x5EED);
            let samples = 1 + draw.next(6) as usize;
            let contigs = 1 + draw.next(2) as u32;
            let bound = MaxCohortLocusSpan(
                std::num::NonZeroU32::new(5 + u32::try_from(draw.next(60)).expect("small"))
                    .expect("at least five"),
            );
            let keep = MinAltReads {
                floor: MinAltObs(
                    std::num::NonZeroU32::new(1 + u32::try_from(draw.next(4)).expect("small"))
                        .expect("at least one"),
                ),
                share: MinAltReadShare::DEFAULT,
            };

            // Every reference base is `A`: two members of one locus that disagree on the
            // reference are refused by `build_region`, in both drivers alike.
            let mut layouts: Vec<Vec<SampleLocusObservations>> = Vec::new();
            for _ in 0..samples {
                let mut records = Vec::new();
                for contig in 0..contigs {
                    let mut at_base = 1 + draw.next(10);
                    while at_base <= ground_end {
                        let bases = match draw.next(10) {
                            0 => 1 + draw.next(60),
                            _ => 1 + draw.next(4),
                        };
                        let end = at_base + bases - 1;
                        let width = usize::try_from(end - at_base + 1).expect("small");
                        records.push(SampleLocusObservations {
                            reference_bases: vec![b'A'; width].into_boxed_slice(),
                            ..member(region_on(contig, at_base, end), b"A", b"T")
                        });
                        at_base = end + 1 + draw.next(6);
                    }
                }
                layouts.push(records);
            }

            let analysed_width = 20 + draw.next(80);
            let mut analysed = Vec::new();
            for contig in 0..contigs {
                let mut at_base = 1u64;
                while at_base <= ground_end {
                    let end = at_base + analysed_width - 1;
                    if draw.next(4) != 0 {
                        analysed.push(region_on(contig, at_base, end));
                    }
                    at_base = end + 1;
                }
            }
            let building_width = width(1 + u32::try_from(draw.next(30)).expect("a small width"));

            if layouts.iter().flatten().any(|record| {
                analysed.iter().any(|ground| {
                    ground.contig == record.region.contig
                        && record.region.start >= ground.start
                        && record.region.start <= ground.end
                        && record.reach() > ground.end
                })
            }) {
                with_a_straddling_record += 1;
            }

            let per_sample: Vec<&[SampleLocusObservations]> =
                layouts.iter().map(Vec::as_slice).collect();
            let merged = the_outcome_both_drivers_agree_on(
                &analysed,
                &per_sample,
                building_width,
                bound,
                keep,
            );
            if !merged.cohort_observations.is_empty() {
                with_observations += 1;
            }
            if !merged.failed_locus_spans.is_empty() {
                with_failed += 1;
            }
        }

        assert!(
            with_observations > 100 && with_failed > 50 && with_a_straddling_record > 100,
            "the generator stopped producing the shapes this test is for: {with_observations} \
             layouts with observations, {with_failed} with failed spans, \
             {with_a_straddling_record} with a record straddling an analysed edge",
        );
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
            MinAltReads::DEFAULT,
        );

        assert!(merged.cohort_observations.is_empty() && merged.failed_locus_spans.is_empty());
    }

    // ---------------------------------------------------------------
    // D2 — the same merge, read through the observation cache. Everything below asserts
    // one thing: the cache changes nothing (the plan's Checkpoint D).
    // ---------------------------------------------------------------

    /// Run both drivers over the same observations, refuse any difference, and hand back the
    /// outcome they agree on.
    ///
    /// **The comparison is on the `Debug` rendering**, which is what "byte-identical" means
    /// here: `CohortObservation` has no `PartialEq`, and a comparison written field by field
    /// would silently stop covering a field added later. Two distinct `f64` sums render as
    /// distinct strings, so a quality divided differently shows.
    ///
    /// **Rendered locus by locus rather than whole**, so that a disagreement names the first
    /// locus that differs: the widest fixture here renders as 26 kB, and two of those inside an
    /// `assert_eq!` message is not something a reader can diff by eye.
    ///
    /// **It asserts rather than returning the material for an assertion.** A call site that
    /// forgot to compare would look right and check nothing, which is the one failure mode this
    /// step's tests exist to prevent. What it returns is the outcome itself, so a test can also
    /// say what its fixture reached — an agreement between two empty outcomes would otherwise
    /// pass for a proof.
    fn the_outcome_both_drivers_agree_on(
        analysed: &[GenomeRegion],
        per_sample: &[&[SampleLocusObservations]],
        building_region_width: CohortLocusBuilderRegionsLen,
        max_cohort_locus_span: MaxCohortLocusSpan,
        min_alt_reads: MinAltReads,
    ) -> RegionOutcome {
        let oracle =
            merge_cohort_serially(analysed, per_sample, max_cohort_locus_span, min_alt_reads);
        let mut cache =
            ObservationCache::over(per_sample.iter().map(|sample| source_of(sample)).collect());
        let through_cache = merge_cohort_through_cache(
            analysed,
            &mut cache,
            building_region_width,
            max_cohort_locus_span,
            min_alt_reads,
        )
        .expect("the fixture sources hold");

        let (from_oracle, from_cache) = (render(&oracle), render(&through_cache));
        let bases = building_region_width.get();
        // **Before the comparison, not after.** Once the two drivers are known equal, the
        // cached output *is* the oracle's, whose loci one walk per analysed region makes
        // disjoint by construction — so a guard placed below could only fire on an input the
        // oracle itself had overlapped, which is to say never. Measured by E2's review: with
        // the cached driver deliberately broken, eight tests failed and none of the eight was
        // this guard; moved here, the same break trips it by name.
        refuse_overlapping_ground(&through_cache, bases);
        refuse_displaced_loci(
            analysed,
            per_sample,
            building_region_width,
            max_cohort_locus_span,
            min_alt_reads,
        );
        refuse_a_parallel_difference(
            analysed,
            per_sample,
            building_region_width,
            max_cohort_locus_span,
            min_alt_reads,
            &oracle,
        );
        if let Some(first) = from_oracle
            .iter()
            .zip(&from_cache)
            .position(|(oracle_entry, cache_entry)| oracle_entry != cache_entry)
        {
            panic!(
                "building regions of {bases} bases changed the merge, first at entry {first}:\
                 \n  oracle:        {}\n  through cache: {}",
                from_oracle[first], from_cache[first],
            );
        }
        assert_eq!(
            from_cache.len(),
            from_oracle.len(),
            "building regions of {bases} bases changed how many loci the merge produced",
        );

        oracle
    }

    /// **The parallel driver's answer is the oracle's, on every fixture this helper sees** —
    /// milestone E's claim (spec §15), made here rather than only on the parallel file's own
    /// fixtures so that it covers what those cannot: the two hundred random layouts, and the
    /// locus built from observations the generic generator actually minted.
    ///
    /// **One, four and sixteen regions in flight, not the full sweep**, because the helper's
    /// busiest caller runs it two hundred times. One is the round that is not a round, which
    /// makes every building-region boundary a round boundary too. Four puts several builders in
    /// one round at every width this file uses **except 600**, where the analysed stretch is a
    /// single building region and every count gives the same one-builder round. Sixteen is here
    /// because E4's review measured what the smaller two miss: a defect confined to a round's
    /// fifth region or later fails five tests and **none** of them is in this file.
    ///
    /// The exhaustive width × count sweep is
    /// `super::super::parallel::tests::the_parallel_merge_is_the_oracles_at_every_width_and_count`;
    /// what this adds is reach — the two hundred random layouts, and the locus built from
    /// observations the generic generator actually minted, neither of which the parallel file
    /// has.
    fn refuse_a_parallel_difference(
        analysed: &[GenomeRegion],
        per_sample: &[&[SampleLocusObservations]],
        building_region_width: CohortLocusBuilderRegionsLen,
        max_cohort_locus_span: MaxCohortLocusSpan,
        min_alt_reads: MinAltReads,
        oracle: &RegionOutcome,
    ) {
        for regions in [1, 4, 16] {
            let mut cache =
                ObservationCache::over(per_sample.iter().map(|sample| source_of(sample)).collect());
            let in_parallel = merge_cohort_in_parallel(
                analysed,
                &mut cache,
                building_region_width,
                in_flight(regions),
                max_cohort_locus_span,
                min_alt_reads,
            )
            .expect("the fixture sources hold");
            refuse_any_difference(
                &format!(
                    "{regions} regions in flight on {}-base regions",
                    building_region_width.get()
                ),
                oracle,
                &in_parallel,
            );
        }
    }

    /// **A real merge never displaces a locus** — the same claim as
    /// [`refuse_overlapping_ground`], made where it counts: through the builders and into a
    /// real [`Organiser`], rather than on the merged output afterwards.
    ///
    /// This is `merge_cohort_through_cache`'s own loop with the organiser wired in — evict at
    /// the building region's first base, cover it, build it, submit it — which is as close to
    /// the parallel arrangement as one thread gets. Without it nothing in the suite ever hands
    /// the organiser an outcome a builder produced, and the counter that is meant to be the
    /// alarm would first be read on real data.
    ///
    /// **It re-runs the driver's loop rather than calling the driver**, because
    /// `merge_cohort_through_cache` returns one merged outcome and the organiser needs them
    /// region by region. So a defect in the *driver's* own loop is caught by the byte-identity
    /// comparison above and not here; what this catches is a builder producing overlapping loci
    /// from a window drawn as the design says to draw it.
    ///
    /// **The eviction point is the discipline it pins**, and E2's review is why that matters:
    /// the safety-net argument reads as though it rests on the cache, and it rests on where
    /// whoever draws the cache forward chooses to evict. Measured here: move this loop's
    /// eviction from the building region's first base to its last — one line — and the suite
    /// reports `at building regions of 29 bases a builder produced a locus on ground an earlier
    /// locus already owned`. The reviewer measured the same change on a harness of its own at
    /// 1,170 displacements over 4,000 random layouts. E3 hands eviction to the organiser with
    /// several builders in flight, where the safe point is the earliest live region's first
    /// base rather than the latest.
    fn refuse_displaced_loci(
        analysed: &[GenomeRegion],
        per_sample: &[&[SampleLocusObservations]],
        building_region_width: CohortLocusBuilderRegionsLen,
        max_cohort_locus_span: MaxCohortLocusSpan,
        min_alt_reads: MinAltReads,
    ) {
        let mut cache =
            ObservationCache::over(per_sample.iter().map(|sample| source_of(sample)).collect());
        let mut organiser = Organiser::new();
        let mut next_index = 0u64;

        for analysed_region in analysed {
            for building_region in building_regions_of(*analysed_region, building_region_width) {
                cache.evict_before(GenomePosition {
                    contig: building_region.contig,
                    position: building_region.start,
                });
                cache
                    .cover(building_region)
                    .expect("the fixture sources hold");
                let outcome = cache.with_observations(building_region, |window| {
                    build_region_windowed(
                        building_region,
                        window,
                        max_cohort_locus_span,
                        min_alt_reads,
                    )
                });
                organiser.submit(RegionIndex(next_index), outcome);
                next_index += 1;
                assert!(organiser.drain_ready().count() < usize::MAX);
            }
        }

        assert_eq!(
            organiser.displaced_locus_count(),
            0,
            "at building regions of {} bases a builder produced a locus on ground an earlier \
             locus already owned",
            building_region_width.get(),
        );
    }

    /// **No two loci a merge produces may overlap** — the claim the organiser's displacement
    /// rule is the safety net for (`super::organise::Organiser`, and spec §6.1).
    ///
    /// Asserted inside [`the_outcome_both_drivers_agree_on`] rather than in a test of its own,
    /// so that it holds over every fixture routed through that helper rather than over one case
    /// written to demonstrate it — **six of this file's twenty-eight tests**, among them the two
    /// hundred random layouts and the 305–330 deletion that a building region's boundary falls
    /// inside at four of its five widths. **A test that calls a driver directly is not checked**,
    /// and `a_locus_reaching_across_a_building_region_boundary_is_built_once_and_whole` is one
    /// such. The organiser's own tests say what the rule *does*; this says it is not needed on
    /// the fixtures it sees.
    ///
    /// It is checked on the **cached** driver's output, because that is the one that divides
    /// the ground into short regions: the oracle hands a whole analysed region to one builder,
    /// where the loci are disjoint by construction and the claim is vacuous.
    ///
    /// A failed span counts as ground exactly as an emitted locus does (spec §3.2), so the two
    /// are sorted together and checked as one sequence. **Disjointness is the whole of what it
    /// asserts**: the sort is how the two vectors are merged, and it also throws away the order
    /// the driver produced them in — which the byte-identity comparison above is what checks.
    fn refuse_overlapping_ground(outcome: &RegionOutcome, bases: u32) {
        let mut ground: Vec<GenomeRegion> = outcome
            .cohort_observations
            .iter()
            .map(|observed| observed.region)
            .chain(outcome.failed_locus_spans.iter().copied())
            .collect();
        ground.sort_by_key(|span| (span.contig, span.start));

        for pair in ground.windows(2) {
            let [earlier, later] = pair else { continue };
            let (earlier, later) = (*earlier, *later);
            assert!(
                earlier.contig != later.contig || earlier.end < later.start,
                "at building regions of {bases} bases the merge produced {earlier} and \
                 {later}, which share ground — the second was built by a builder that never \
                 saw what opened in the first",
            );
        }
    }

    /// The disjointness guard above asserts a property every fixture in this file has, so
    /// nothing here can tell a working guard from a vacuous one. These two feed it output no
    /// merge produces, which is the only way to know it would fire.
    #[test]
    #[should_panic(expected = "which share ground")]
    fn the_disjointness_guard_refuses_two_overlapping_loci() {
        refuse_overlapping_ground(
            &RegionOutcome {
                cohort_observations: Vec::new(),
                failed_locus_spans: vec![region(10, 40), region(25, 60)],
            },
            20,
        );
    }

    /// And a failed span counts as ground, so an emitted locus opening inside one is refused
    /// exactly as two emitted loci would be — which is why the two vectors are sorted together
    /// rather than checked apart.
    #[test]
    #[should_panic(expected = "which share ground")]
    fn the_disjointness_guard_weighs_a_failed_span_as_ground() {
        refuse_overlapping_ground(
            &RegionOutcome {
                cohort_observations: vec![locus_over(region(40, 44))],
                failed_locus_spans: vec![region(10, 200)],
            },
            20,
        );
    }

    /// **The boundary is one shared base.** `GenomeRegion` is inclusive at both ends, so
    /// 10-40 and 40-60 share base 40 — the case an `earlier.end <= later.start` slip would
    /// call disjoint.
    #[test]
    #[should_panic(expected = "which share ground")]
    fn the_disjointness_guard_refuses_two_loci_sharing_one_base() {
        refuse_overlapping_ground(
            &RegionOutcome {
                cohort_observations: Vec::new(),
                failed_locus_spans: vec![region(10, 40), region(40, 60)],
            },
            20,
        );
    }

    /// And the other side of that boundary: adjacent is not overlapping, so 10-39 beside
    /// 40-60 passes.
    #[test]
    fn the_disjointness_guard_allows_two_adjacent_loci() {
        refuse_overlapping_ground(
            &RegionOutcome {
                cohort_observations: Vec::new(),
                failed_locus_spans: vec![region(10, 39), region(40, 60)],
            },
            20,
        );
    }

    /// **The sort is by contig first, and that carries weight.** Ordered on position alone
    /// these three interleave, every adjacent pair straddles a contig, and the overlap on
    /// contig 1 is never compared.
    #[test]
    #[should_panic(expected = "which share ground")]
    fn the_disjointness_guard_finds_an_overlap_interleaved_across_contigs() {
        refuse_overlapping_ground(
            &RegionOutcome {
                cohort_observations: Vec::new(),
                failed_locus_spans: vec![
                    region_on(1, 10, 80),
                    region_on(0, 20, 25),
                    region_on(1, 30, 40),
                ],
            },
            20,
        );
    }

    /// A cohort observation over `span`, for the guard's own tests — what it reads is the
    /// ground, and nothing else.
    fn locus_over(span: GenomeRegion) -> crate::ng::run::cohort_merge::build::CohortObservation {
        crate::ng::run::cohort_merge::build::CohortObservation {
            region: span,
            alleles: Vec::new(),
            per_sample: Vec::new(),
            kind: LocusKind::Generic,
        }
    }

    /// **The cache changes nothing, at every building-region width** — the milestone's claim.
    ///
    /// Three samples over six hundred bases: single-base loci every ten bases, a deletion at
    /// 305–330 that carries a locus across whatever boundary a width puts near it, and a
    /// sample whose one record sits inside that deletion, so the two chain into one locus
    /// however finely the ground is divided.
    ///
    /// The widths run from one base — where a building region is narrower than most loci — to
    /// six hundred, the whole analysed stretch as one region, which is what the oracle does.
    #[test]
    fn the_cache_changes_nothing_at_every_building_region_width() {
        let layouts = three_samples_over_six_hundred_bases();
        let (dotted, deleting, inside) = (&layouts[0], &layouts[1], &layouts[2]);
        let analysed = [region(1, 600)];
        let per_sample: [&[SampleLocusObservations]; 3] = [dotted, deleting, inside];

        for bases in [1, 3, 20, 47, 600] {
            let merged = the_outcome_both_drivers_agree_on(
                &analysed,
                &per_sample,
                width(bases),
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
            );
            assert_eq!(
                merged.cohort_observations.len(),
                59,
                "the fixture's own shape: 60 dotted loci, two of which the deletion \
                 swallowed into one locus with the sample at 310",
            );
            // **The fixture's one partial reaches the comparison, and that is a premise rather
            // than a detail.** Every driver comparison in this module is on the whole rendering
            // of the outcome, so a field that is empty in every entry is a field two drivers
            // agree on by both building nothing. This is the assertion that keeps the fixture
            // from drifting back to one where `partials` is that field.
            assert_eq!(
                merged
                    .cohort_observations
                    .iter()
                    .flat_map(|observed| observed.per_sample.iter())
                    .filter(|sample| !sample.partials.is_empty())
                    .count(),
                1,
                "the deleting sample's record carries the fixture's only partial",
            );
        }
    }

    /// **A locus that opens in one building region and reaches into the next is built once,
    /// whole, by the region its first position falls in** — the case a window makes possible
    /// to get wrong, since the builder that owns it has to read past its own end.
    ///
    /// The deletion opens at 305 and reaches 330; at twenty-base regions the boundary at 320
    /// falls inside it, and the sample at 310 is what makes the locus wider than the deletion
    /// alone would need. The assertion is on the span, not merely on the count: a locus cut at
    /// the boundary would still be one locus.
    #[test]
    fn a_locus_reaching_across_a_building_region_boundary_is_built_once_and_whole() {
        let deleting = [member(region(305, 330), &[b'A'; 26], b"A")];
        let inside = [member(region(310, 310), b"A", b"T")];

        let mut cache = ObservationCache::over(vec![source_of(&deleting), source_of(&inside)]);
        let merged = merge_cohort_through_cache(
            &[region(1, 600)],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        )
        .expect("the fixture sources hold");

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(305, 330)],
            "one locus over its whole span, owned by the region holding base 305",
        );
        assert_eq!(
            merged.cohort_observations[0].per_sample.len(),
            2,
            "and both samples' evidence is in it, though they sit either side of base 320",
        );

        // **The shape the plan's E2 names**, taken all the way into a real organiser: a wide
        // deletion beginning before a building region and reaching into it is exactly where a
        // later builder would work from a partial picture. Nothing is displaced, because the
        // builder that owns base 305 follows the locus to 330 and every later builder skips
        // the chain it opened.
        refuse_displaced_loci(
            &[region(1, 600)],
            &[&deleting, &inside],
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );
    }

    /// **The building regions stop where the analysed ground stops.** A locus opening at 55 is
    /// outside an analysed region ending at 50, so neither driver builds it — but a division
    /// that let the last building region run to its full width would run to 60 and claim it.
    #[test]
    fn a_locus_past_the_analysed_ground_is_not_claimed_by_the_last_building_region() {
        let sample = [
            member(region(45, 45), b"G", b"T"),
            member(region(55, 55), b"G", b"T"),
        ];
        let analysed = [region(1, 50)];

        let merged = the_outcome_both_drivers_agree_on(
            &analysed,
            &[&sample],
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region(45, 45)],
            "45 is inside the analysed ground and 55 is not",
        );
    }

    /// The failed spans and the quiet ground come through the cache exactly as they come
    /// through the oracle, across several analysed regions and on both sides of a boundary.
    #[test]
    fn the_failed_spans_come_through_the_cache_unchanged() {
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
        let bound = MaxCohortLocusSpan(std::num::NonZeroU32::new(20).expect("20 is non-zero"));

        let merged = the_outcome_both_drivers_agree_on(
            &analysed,
            &[&wide, &quiet],
            width(7),
            bound,
            MinAltReads::DEFAULT,
        );

        assert_eq!(
            merged.failed_locus_spans,
            vec![region(20, 40), region(60, 85)],
            "the fixture is the one that refuses two loci, so the comparison has something \
             to compare — and the single read at 95 reaches neither vector",
        );
        assert!(merged.cohort_observations.is_empty());
    }

    /// **Several contigs, read through one cache.** Each sample's reader crosses the contig
    /// boundary forward, and the eviction at the next contig's first base drops the previous
    /// contig's window — the case where a position comparison blind to the contig would keep
    /// everything.
    #[test]
    fn several_contigs_come_through_the_cache_unchanged() {
        let first = [
            member(region_on(0, 20, 20), b"G", b"T"),
            member(region_on(1, 900, 900), b"G", b"C"),
        ];
        let second = [
            member(region_on(0, 21, 21), b"G", b"C"),
            member(region_on(1, 40, 40), b"G", b"T"),
        ];
        let analysed = [region_on(0, 1, 100), region_on(1, 1, 1000)];

        let merged = the_outcome_both_drivers_agree_on(
            &analysed,
            &[&first, &second],
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![
                region_on(0, 20, 20),
                region_on(0, 21, 21),
                region_on(1, 40, 40),
                region_on(1, 900, 900),
            ],
            "the fixture reaches the second contig, and 900 is past every position on the \
             first",
        );
    }

    /// **The cache changes nothing on observations the generator actually minted** — the same
    /// two samples on disk as Checkpoint C's own fixture, at a building-region width narrower
    /// than the locus they chain into.
    #[test]
    fn the_cache_changes_nothing_on_minted_observations() {
        let analysed = region_on(1, 80, 140);

        let substituted: Vec<_> = (0..3)
            .map(|read| {
                let mut record = read_with_a_substitution(&format!("sub{read}"), 95, 17);
                if read == 0 {
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

        let merged = the_outcome_both_drivers_agree_on(
            &[analysed],
            &[&first, &second],
            width(4),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );

        assert_eq!(
            merged
                .cohort_observations
                .iter()
                .map(|observed| observed.region)
                .collect::<Vec<_>>(),
            vec![region_on(1, 109, 114)],
            "the fixture is Checkpoint C's own locus, so there is something to compare",
        );
        assert_eq!(
            merged.cohort_observations[0]
                .alleles
                .iter()
                .map(|allele| String::from_utf8_lossy(allele).into_owned())
                .collect::<Vec<_>>(),
            vec!["AAAAAA", "AAACAA", "A"],
            "and it is the six-base locus the two samples chained into, alleles and all",
        );
    }

    /// A source's failure ends the merge and comes back unchanged, rather than yielding a
    /// short answer that looks like a cohort with nothing to say.
    #[test]
    fn a_failing_source_ends_the_merge_through_the_cache() {
        let sample = [member(region(12, 12), b"G", b"T")];
        let failing = vec![
            Ok(member(region(12, 12), b"G", b"T")),
            Err(SourceFailed("the block would not decode")),
        ]
        .into_iter();

        let mut cache = ObservationCache::over(vec![source_of(&sample), failing]);
        let outcome = merge_cohort_through_cache(
            &[region(1, 600)],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );

        // `RegionOutcome` has no `PartialEq` — deliberately, since every comparison of one is
        // made through its `Debug` rendering — so the failure is matched rather than equated.
        assert_eq!(
            outcome.err(),
            Some(SourceFailed("the block would not decode")),
        );
    }

    /// Analysed regions that overlap are refused here for the same reason the oracle refuses
    /// them: the loci in the ground they share would be carried twice.
    #[test]
    #[should_panic(expected = "not disjoint and ascending")]
    fn analysed_regions_that_overlap_are_refused_through_the_cache() {
        let sample = [member(region(45, 45), b"G", b"T")];

        let mut cache = ObservationCache::over(vec![source_of(&sample)]);
        let _ = merge_cohort_through_cache(
            &[region(1, 60), region(40, 100)],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );
    }

    /// **The driver evicts as it goes, and the merge's own output cannot say so** — a driver
    /// that never evicted would produce exactly these observations while holding the whole
    /// stretch. So the claim is made against the cache: after merging six hundred bases at
    /// twenty-base regions, with a record every ten bases, the window holds **2 of the 60**.
    ///
    /// **Two is what the last building region holds**, not a bound a forward reader pays: the
    /// records at 581 and 591 both lie inside 581–600, and the source is spent after 591, so
    /// nothing was drawn past the analysed ground. On this same fixture at five-base regions the
    /// window ends holding none — the number moves with the width and the record spacing. The
    /// load-bearing comparison is against the sixty the stretch holds when nothing is evicted.
    ///
    /// **The order — evict, then cover — is not what this pins.** Evicting after the cover ends
    /// with the same window and the same output; what it costs is the sweep, which re-reads the
    /// held window (`ObservationCache::cover`), and on this fixture the pre-cover eviction takes
    /// the window from three records to one at every region.
    #[test]
    fn the_driver_evicts_as_it_goes_and_the_window_stays_short() {
        let dotted: Vec<_> = (0..60)
            .map(|locus| {
                let at = 10 * locus + 1;
                member(region(at, at), b"A", b"T")
            })
            .collect();
        let mut cache = ObservationCache::over(vec![source_of(&dotted)]);

        let merged = merge_cohort_through_cache(
            &[region(1, 600)],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        )
        .expect("the fixture source holds");

        assert_eq!(
            merged.cohort_observations.len(),
            60,
            "every locus was built"
        );
        assert_eq!(
            cache.held_observations_len(),
            2,
            "the records at 581 and 591, not the sixty the stretch holds",
        );
    }

    /// **The width the caller asks for is the width the driver builds in** — which the output
    /// cannot show, since both answers are identical. What differs is the cache: at twenty-base
    /// regions the window ends at two records, and handing the whole six hundred bases to one
    /// builder leaves all sixty held.
    ///
    /// Found by the D2 review: mutating the driver's *call site* to build each analysed region
    /// as one is invisible to the test that checks the division itself, because that test calls
    /// the divider directly and nothing said the driver used it.
    #[test]
    fn the_width_the_caller_asks_for_is_the_width_the_driver_builds_in() {
        let dotted: Vec<_> = (0..60)
            .map(|locus| {
                let at = 10 * locus + 1;
                member(region(at, at), b"A", b"T")
            })
            .collect();

        let held_at = |bases: u32| {
            let mut cache = ObservationCache::over(vec![source_of(&dotted)]);
            merge_cohort_through_cache(
                &[region(1, 600)],
                &mut cache,
                width(bases),
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
            )
            .expect("the fixture source holds");
            cache.held_observations_len()
        };

        assert_eq!(
            held_at(600),
            60,
            "one region for the stretch holds all of it"
        );
        assert_eq!(
            held_at(20),
            2,
            "and twenty-base regions hold only the last region's records",
        );
    }

    /// **The window is short at the moment a merge fails, not merely at the end of one** — the
    /// only place from outside where the driver's *drawing pace* is visible. A driver that
    /// covered the whole analysed region instead of each building region would evict down to the
    /// same window by the end, and every other test here would still pass.
    #[test]
    fn the_window_stays_short_up_to_a_failure() {
        let mut records: Vec<Result<SampleLocusObservations, SourceFailed>> = (0..30)
            .map(|locus| {
                let at = 10 * locus + 1;
                Ok(member(region(at, at), b"A", b"T"))
            })
            .collect();
        records.push(Err(SourceFailed("the block would not decode")));
        let mut cache = ObservationCache::over(vec![records.into_iter()]);

        let outcome = merge_cohort_through_cache(
            &[region(1, 600)],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );

        assert!(outcome.is_err());
        assert_eq!(
            cache.held_observations_len(),
            2,
            "the window at the moment of failure, not the thirty records behind it",
        );
    }

    /// **An analysed region whose ends are the wrong way round is refused**, rather than read
    /// one way by the division and another by the builder: the division orders the two ends and
    /// `build_region` does not, so `50-1` builds nothing through the oracle and the whole of
    /// 1–50 through the cache. Found by the D2 review, and it is the byte-identity this module
    /// claims, broken by an input nothing was checking.
    #[test]
    #[should_panic(expected = "the wrong way round")]
    fn an_analysed_region_with_inverted_ends_is_refused() {
        let sample = [member(region(45, 45), b"G", b"T")];
        let inverted = GenomeRegion {
            contig: ContigId(0),
            start: Position(50),
            end: Position(1),
        };
        let mut cache = ObservationCache::over(vec![source_of(&sample)]);

        let _ = merge_cohort_through_cache(
            &[inverted],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );
    }

    /// **Regions sharing exactly one base are refused**, which is the comparison the guard turns
    /// on: at `<=` in place of `<` this pair is accepted and the locus at 50 is built by both
    /// regions and carried twice. The other refusal fixtures overlap by 21 bases, where a
    /// loosened comparison still refuses.
    #[test]
    #[should_panic(expected = "not disjoint and ascending")]
    fn analysed_regions_sharing_one_base_are_refused() {
        let sample = [member(region(50, 50), b"G", b"T")];
        let mut cache = ObservationCache::over(vec![source_of(&sample)]);

        let _ = merge_cohort_through_cache(
            &[region(1, 50), region(50, 100)],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );
    }

    /// Nothing analysed is not an error here either, and nothing is drawn — so a driver that
    /// covered or evicted before the loop would show.
    #[test]
    fn merging_no_analysed_regions_through_the_cache_yields_nothing() {
        let sample = [member(region(12, 12), b"G", b"T")];
        let mut cache = ObservationCache::over(vec![source_of(&sample)]);

        let merged = merge_cohort_through_cache(
            &[],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        )
        .expect("nothing can fail");

        assert!(merged.cohort_observations.is_empty() && merged.failed_locus_spans.is_empty());
        assert_eq!(cache.held_observations_len(), 0, "and nothing was drawn");
    }

    /// A cohort of no samples merges to nothing rather than failing — the bottom of the range
    /// this caller commits to (spec §7.2), and the shape `ObservationCache::over(Vec::new())`
    /// produces.
    #[test]
    fn merging_no_samples_through_the_cache_yields_nothing() {
        let mut cache: ObservationCache<
            std::vec::IntoIter<Result<SampleLocusObservations, SourceFailed>>,
        > = ObservationCache::over(Vec::new());

        let merged = merge_cohort_through_cache(
            &[region(1, 100)],
            &mut cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        )
        .expect("nothing can fail");

        assert!(merged.cohort_observations.is_empty() && merged.failed_locus_spans.is_empty());
    }

    /// The keep threshold reaches the builders through the cache too — every other test here
    /// runs at the default, so a driver that dropped the argument would ship green.
    #[test]
    fn the_keep_threshold_reaches_the_builders_through_the_cache() {
        let sample = [member(region(12, 12), b"G", b"T")];

        let mut built_cache = ObservationCache::over(vec![source_of(&sample)]);
        let built = merge_cohort_through_cache(
            &[region(1, 50)],
            &mut built_cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        )
        .expect("the fixture source holds");
        let mut quiet_cache = ObservationCache::over(vec![source_of(&sample)]);
        let too_quiet = merge_cohort_through_cache(
            &[region(1, 50)],
            &mut quiet_cache,
            width(20),
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads {
                floor: MinAltObs(std::num::NonZeroU32::new(4).expect("4 is non-zero")),
                share: MinAltReadShare::DEFAULT,
            },
        )
        .expect("the fixture source holds");

        assert_eq!(built.cohort_observations.len(), 1, "three reads reach two");
        assert!(too_quiet.cohort_observations.is_empty(), "and not four");
    }

    /// Drive the parallel-cover form over `layouts` inside a pool of `threads`, collecting
    /// what it hands over — the shape every parallel-cover test here compares.
    fn merge_with_parallel_cover_in_a_pool(
        threads: usize,
        analysed: &[GenomeRegion],
        layouts: &[Vec<SampleLocusObservations>],
        building_region_width: CohortLocusBuilderRegionsLen,
    ) -> RegionOutcome {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .expect("a fixture pool");
        pool.install(|| {
            let mut cache = ObservationCache::over(
                layouts
                    .iter()
                    .map(|sample| source_of(sample))
                    .collect::<Vec<_>>(),
            );
            let mut merged = RegionOutcome::default();
            let RegionOutcome {
                cohort_observations,
                failed_locus_spans,
            } = &mut merged;
            merge_cohort_handing_each_locus_over_covering_samples_in_parallel(
                analysed,
                &mut cache,
                building_region_width,
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
                &mut |built| cohort_observations.push(built),
                failed_locus_spans,
            )
            .expect("the fixture sources hold");
            merged
        })
    }

    /// **The parallel cover gives the serial drivers' answer, at every width and every pool
    /// size** — Milestone E1's claim at the driver, on the fixture the module's byte-identity
    /// rests on, which carries a locus chaining two samples across region boundaries and a
    /// span the width bound refuses at the narrow widths.
    ///
    /// The pool of one exercises the fallback branch (a one-thread pool takes the serial
    /// sweep); the pools of two and eight exercise the Jacobi sweep with genuinely concurrent
    /// samples. **Whether the sweep truly occupies several threads is untested here**, exactly
    /// as the parallel merge's own module records for its builders — the claim under test is
    /// that the schedule cannot change the answer.
    #[test]
    fn the_parallel_cover_gives_the_serial_drivers_answer() {
        let layouts = {
            let mut layouts = three_samples_over_six_hundred_bases();
            // A span the 50-base default bound refuses, so the refused list is compared too.
            layouts.push(vec![member(region(420, 510), &[b'A'; 91], b"A")]);
            layouts
        };
        let analysed = [region(1, 600)];
        let per_sample: Vec<&[SampleLocusObservations]> =
            layouts.iter().map(Vec::as_slice).collect();
        let oracle = merge_cohort_serially(
            &analysed,
            &per_sample,
            MaxCohortLocusSpan::DEFAULT,
            MinAltReads::DEFAULT,
        );
        assert_eq!(
            oracle.failed_locus_spans,
            vec![region(420, 510)],
            "the fixture must carry a locus the width bound refuses",
        );

        for threads in [1, 2, 8] {
            for bases in [1, 3, 20, 47, 600] {
                refuse_any_difference(
                    &format!("the parallel cover on {bases}-base regions in a pool of {threads}"),
                    &oracle,
                    &merge_with_parallel_cover_in_a_pool(
                        threads,
                        &analysed,
                        &layouts,
                        width(bases),
                    ),
                );
            }
        }
    }

    /// A reader that fails ends the merge under the parallel cover, and its own error comes
    /// back untouched — the same contract as both serial drivers and the parallel merge.
    #[test]
    fn a_failing_source_ends_the_merge_under_the_parallel_cover() {
        let failing = vec![
            Ok(member(region(5, 5), b"G", b"T")),
            Err(SourceFailed("the walk failed")),
        ];
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("a fixture pool");

        let merged = pool.install(|| {
            let mut cache = ObservationCache::over(vec![failing.into_iter()]);
            let mut refused = Vec::new();
            merge_cohort_handing_each_locus_over_covering_samples_in_parallel(
                &[region(1, 600)],
                &mut cache,
                width(20),
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
                &mut |_built| {},
                &mut refused,
            )
        });

        assert_eq!(merged, Err(SourceFailed("the walk failed")));
    }
}
