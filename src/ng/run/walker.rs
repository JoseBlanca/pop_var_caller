//! One sample's alignment files behind the merge's source interface.
//!
//! **A source answers one question: what did this sample see next?** The merge asks each of its
//! k samples in turn and never asks twice, so what it needs from a sample is one observation at
//! a time, in coordinate order, forward only, for the whole run
//! (`doc/devel/ng/arch/run_streaming.md` §2). The two calling modes differ in exactly this one
//! place — direct mode reads the observations out of alignment files, psp mode decodes them from
//! a file — and nothing above the trait can tell which it is holding.
//!
//! This file is direct mode's half: the walker. **It is not new machinery.**
//! [`SampleLocusObservationsIterator`] already drives the typed-region generators over one
//! sample's reads and yields loci in genome order; what the walker adds is the two things the
//! merge's trait asks for that a plain iterator cannot give.
//!
//! - **A failure that names the sample and where the walk had reached.** The merge adds nothing
//!   to a source's error and passes it through, so the error has to locate itself: in a run over
//!   a thousand samples, *reading a region failed* names neither the individual to look at nor
//!   the place to look ([`RunError::SourceFailed`], spec §9).
//! - **The offer of a spent record back for reuse.** [`ObservationSource::next_observation`]
//!   hands the walker a record the merge will not read again. **This walker drops it**, which
//!   the trait explicitly permits — the offer is not an obligation. A later step of the plan
//!   (`doc/devel/ng/impl_plan/run_driver_direct_mode.md`, step G1) refills it instead, and that
//!   is where the measurement motivating it lives. **What this step buys that step** is that
//!   the trait is implemented by this type rather than picked up from the blanket
//!   implementation: the blanket one covers *every* iterator of one sample's observations, so
//!   it can never be specialised for this walker alone, and a step that wants to refill the
//!   spare needs an implementation of its own to change.
//!
//! **One walker per sample for the whole run** — not one per worker and not one per segment
//! (spec §3.4). The merge is its only consumer and the merge only moves forward, so each
//! stretch of a file is decoded once and the backward jump the cursors are capable of never
//! happens. It is also why the walker holds its generators for the run: a locus generator
//! carries state across segments (spec §8) — the generic one keeps a read cursor per
//! chromosome, which `end_walk` does not clear
//! ([`locus_generation/mod.rs`](../../../../src/ng/locus_generation/mod.rs), the `Drop`
//! commentary) — and a fresh generator at each of a run's segments would re-open and re-seek
//! per segment where the cursor answers from where it already is
//! ([`read/input/mod.rs`](../../../../src/ng/read/input/mod.rs), `SampleReads::cursor`).
//!
//! **The walker is deliberately not an [`Iterator`].** Every iterator of one sample's
//! observations is already a source, through a blanket implementation that drops the spare
//! ([`observation_cache`](super::cohort_merge::observation_cache)) — so a type that was both
//! would implement the trait twice, and Rust refuses the overlap. Nothing is lost: the walk is
//! the same machinery as [`SampleLocusObservationsIterator`], so the differential that checks a
//! source against the walk builds that iterator directly over the same fixture.

use std::convert::Infallible;

use crate::ng::locus_generation::{
    GeneratorSet, LocusCounts, LocusGenerationError, SampleLocusObservations,
    SampleLocusObservationsIterator,
};
use crate::ng::read::input::SampleReads;
use crate::ng::region_typing::TypedRegion;

use super::cohort_merge::observation_cache::ObservationSource;
use super::{RunError, Segmentation, WalkProgress};

/// One sample's observations, read from its open alignment files.
///
/// Constructed once per sample and advanced by the merge alone. See the module's own
/// documentation for why it is one per sample and why it is not an [`Iterator`].
///
/// **Generic over the region stream, exactly as the iterator it wraps is**, so a test can hand
/// it a `Vec` of segments where a run hands it [`RunSegments`] over the whole segmentation.
///
/// **Not `Clone`**: it owns one sample's open files and its generators' accumulated state, and
/// a second walker over one sample would decode the same ground twice while each told the merge
/// a different story about how far it had got.
///
/// **Neither `Send` nor `Sync`, and it cannot be made so from here.** Arch §2's last contract
/// clause says a source needs those two only for the parallel merge
/// ([`merge_cohort_in_parallel`](super::cohort_merge::parallel), which bounds
/// `S: Sync + Send`); the single-threaded merge this is built against does not. What blocks it
/// is one layer down: [`GeneratorSlot::Generator`](crate::ng::locus_generation::GeneratorSlot)
/// holds a `Box<dyn LocusGenerator<S>>` with no auto-trait bound, a deliberate omission its own
/// documentation records. So a walker and the merge that draws from it stay on one thread, and
/// putting a walker under the parallel merge means widening that trait object first — stated
/// here so it is read rather than met as a compiler error.
pub struct AlignmentFilesWalker<T> {
    /// The individual this walk is of. **Held rather than read back from the reads**, because
    /// the iterator owns the [`SampleReads`] and does not hand it out, and a failure has to name
    /// the sample after the walk is under way.
    sample: String,
    /// How far the walk has got — the second half of locating a failure (spec §9).
    reached: WalkProgress,
    loci: SampleLocusObservationsIterator<T>,
}

impl<T> AlignmentFilesWalker<T> {
    /// A walker over `regions`, reading `reads` through `generators`.
    ///
    /// **The sample's name is taken from the reads rather than passed in**, so the name a
    /// failure carries and the files it failed on cannot come apart.
    pub fn new(regions: T, reads: SampleReads, generators: GeneratorSet) -> Self {
        Self {
            sample: reads.sample_name().to_string(),
            reached: WalkProgress::NothingYet,
            loci: SampleLocusObservationsIterator::new(regions, reads, generators),
        }
    }

    /// The individual this walk is of.
    #[must_use]
    pub fn sample_name(&self) -> &str {
        &self.sample
    }

    /// How far this walk has got. `NothingYet` until the first observation is yielded, and the last
    /// base of the last one after that — **not** the position the walk is decoding, which is ahead
    /// of it and which nothing outside the generators knows.
    #[must_use]
    pub fn reached(&self) -> WalkProgress {
        self.reached
    }

    /// The walk's running tally — current at any point, final once the walk is spent. What the
    /// run report says about regions handled and regions refused comes from here.
    #[must_use]
    pub fn counts(&self) -> &LocusCounts {
        self.loci.counts()
    }

    /// The generator set, for the per-generator counts [`counts`](Self::counts) does not carry.
    ///
    /// **Forwarded because the walker would otherwise close the only door to them.** Once a run
    /// owns its walkers, nothing else holds these generators, and
    /// [`GeneratorSet::generic_counts`](crate::ng::locus_generation::GeneratorSet::generic_counts)
    /// and its siblings — the nine counts that explain a covered region emitting nothing — would
    /// have no reader outside a test. Handed out by `&`, so a caller can read every slot and
    /// change none.
    ///
    /// **The read-filter tallies are still out of reach, and that is a gap this cannot close.**
    /// They belong to a cursor from the moment it is made (spec §8), and cursors live inside the
    /// generators; neither this type nor [`SampleReads`] hands one out. Whatever builds the run
    /// report needs a route to them that does not exist yet.
    #[must_use]
    pub fn generators(&self) -> &GeneratorSet {
        self.loci.generators()
    }
}

impl<'a> AlignmentFilesWalker<RunSegments<'a>> {
    /// A walker over **the whole run's ground**: every segment of the segmentation, in genome
    /// order, once.
    ///
    /// This is the shape a run builds. Every sample of the run gets the same segments from the
    /// same object, which is what makes "k samples over one segmentation" true rather than "k
    /// segmentations that happen to agree" (spec §4.2).
    pub fn over(
        segmentation: &'a Segmentation,
        reads: SampleReads,
        generators: GeneratorSet,
    ) -> Self {
        Self::new(RunSegments::of(segmentation), reads, generators)
    }
}

/// **The sample's name and how far it has got, not its reads.** A derived `Debug` would print
/// every open file and every generator's accumulated state.
impl<T> std::fmt::Debug for AlignmentFilesWalker<T> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("AlignmentFilesWalker")
            .field("sample", &self.sample)
            .field("reached", &self.reached)
            .finish_non_exhaustive()
    }
}

impl<T, E> ObservationSource for AlignmentFilesWalker<T>
where
    T: Iterator<Item = Result<TypedRegion, E>>,
    E: Into<LocusGenerationError>,
{
    type Error = RunError;

    /// The next observation this sample has, or `None` once its ground is walked.
    ///
    /// **The spare is dropped**, which the trait permits; a later step of the plan refills it
    /// instead (see the module documentation). It is released explicitly at the top rather than
    /// left to fall out of scope, so its buffers are back with the allocator before the walk's
    /// next draw asks for any.
    ///
    /// **⚑ What the tests here pin is that the spare does not come back out as an observation,
    /// not that it is released.** A walker that stashed every offered record and never freed one
    /// passes all of them — measured, that mutation survived a 21-mutation pass that killed the
    /// other twenty. Nothing can pin the release while there is no pool to count, so the test
    /// belongs to the step that adds one: at 63 samples it is unbounded growth of exactly the
    /// records the reuse hook exists to stop allocating, so the step that starts keeping them
    /// must bound how many it keeps and assert the bound.
    ///
    /// **Exhaustion is final**, which is the wrapped iterator's guarantee rather than a fresh
    /// one: it latches a `done` flag and is a [`FusedIterator`](std::iter::FusedIterator). A
    /// source that yielded `Some` after a `None` would be drawn in behind the merge's window and
    /// so silently out of coordinate order.
    ///
    /// **⚑ A failed walk is spent, and the trait's contract says a failure leaves a source
    /// live.** The two are not the same, and this implementation takes the first: the wrapped
    /// iterator latches `done` on an error, so a consumer that swallowed the error and asked
    /// again would get `None`. The trait's own note says what live-after-failure is for — it is
    /// "what lets a cover be made again"
    /// ([`observation_cache`](super::cohort_merge::observation_cache)) — and **nothing in the
    /// merge does that**: `ObservationCache::draw_next` propagates the error without marking the
    /// source spent, and both drivers abandon the cache rather than retry. So the deviation is
    /// unreachable today. What it would cost if something did retry is worth stating, because it
    /// is silent: the cache would read the `None` as exhaustion, mark the sample spent, and go on
    /// building cohort loci without it — wrong genotypes, not an error. Anything that adds a
    /// retry must fix this first.
    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, RunError>> {
        drop(spare);
        match self.loci.next()? {
            Ok(observation) => {
                // **`reach_position` rather than `region.end`**, because the crate keeps that
                // rule in one place: `reach` is `end.max(start)`, and `GenomeRegion` has public
                // fields and no constructor, so an inverted region read straight off `end` would
                // put an observation's reach before its own first base. The merge's own cache
                // keys on the same call.
                self.reached = WalkProgress::After(observation.reach_position());
                Some(Ok(observation))
            }
            Err(source) => Some(Err(RunError::SourceFailed {
                sample: self.sample.clone(),
                reached: self.reached,
                source: Box::new(source),
            })),
        }
    }
}

/// The run's segments as a walker reads them: a cursor over [`Segmentation`]'s list, in genome
/// order, handing out one segment at a time.
///
/// **A borrow rather than a copy.** The list grows with the genome, not with the cohort:
/// 100,171 segments over the 80 BED regions of `benchmarks/tomato1`, which that benchmark's 63
/// accessions all share. A copy per sample would multiply one list by the cohort size for no
/// gain, since every walker reads the same segments in the same order.
///
/// **Its item is `Result<_, Infallible>` because reading it cannot fail.** A run's segments were
/// read out of the repeat catalog once, at [`Segmentation::build`], and every catalog failure
/// was reported there. The wrapped iterator takes a fallible stream because the catalog's own
/// reader is one; this says, in the type, that this particular stream has nothing left to fail
/// at.
pub struct RunSegments<'a> {
    remaining: std::slice::Iter<'a, TypedRegion>,
}

impl<'a> RunSegments<'a> {
    /// Every segment of `segmentation`, in genome order.
    #[must_use]
    pub fn of(segmentation: &'a Segmentation) -> Self {
        Self {
            remaining: segmentation.segments().iter(),
        }
    }
}

impl Iterator for RunSegments<'_> {
    type Item = Result<TypedRegion, Infallible>;

    /// **Cloned rather than borrowed**, because a generator is *given* the region it begins and
    /// holds it while it drains.
    ///
    /// What a clone costs: the span and the kind copy inline — a
    /// [`Motif`](crate::ng::types::Motif) is a six-byte buffer and a length, not a heap
    /// string — and two of the four kinds allocate. A repeat tract copies its contig name
    /// (`SsrSegment`'s `chrom`, a `Box<str>`), and a bundle copies its tract list
    /// (`RegionKind::SsrBundle`'s `Box<[RepeatInterval]>`, two or more members). So at most two
    /// allocations a segment, paid once per segment rather than once per locus.
    fn next(&mut self) -> Option<Self::Item> {
        self.remaining.next().cloned().map(Ok)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::locus_generation::{GeneratorSlot, LocusGenerator, LocusKind, UnhandledReason};
    use crate::ng::read::filtering::ReadFilterConfig;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference, header, indexed_bam, matching_contigs, read_named_with_length,
    };
    use crate::ng::region_typing::{GenomeRegions, RegionKind};
    use crate::ng::repeat_catalog::{RepeatCatalogHeader, StrRepeatCriteria};
    use crate::ng::tandem_repeat::ScanParams;
    use crate::ng::types::{ContigId, GenomePosition, GenomeRegion, Position};
    use crate::regions::ContigBounds;
    use std::path::PathBuf;

    // -----------------------------------------------------------------
    // Fixtures
    //
    // **Everything here differs from everything else it could be confused with.** Two samples
    // with different names, segments on two different contigs, loci whose start and end differ,
    // and a locus count that differs per segment — because a walker over one segment of one
    // sample cannot tell "advances to the next segment" from "hands back everything it has",
    // and one whose loci are all at the same position cannot tell "the end of the last
    // observation" from "the start of the first".
    // -----------------------------------------------------------------

    /// A one-read indexed BAM naming `sample`, opened as its own cohort.
    ///
    /// The temp dirs come back because they must outlive the reads.
    fn reads_named(sample: &str) -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        let (reference_dir, reference) = fixture_reference(false);
        let records = vec![read_named_with_length("r0", 0, 1, 30)];
        let (bam_dir, bam_path) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some(sample))],
            ),
            &records,
        );
        let reads = SampleReads::open_only_sample(
            &[bam_path],
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    /// The segments as a walker's region stream. **`Infallible`**, because a list already in
    /// hand has nothing left to fail at — the same statement `RunSegments` makes about a run's
    /// own segments.
    fn stream(segments: Vec<TypedRegion>) -> std::vec::IntoIter<Result<TypedRegion, Infallible>> {
        segments.into_iter().map(Ok).collect::<Vec<_>>().into_iter()
    }

    fn generic_segment(contig: u32, start: u64, end: u64) -> TypedRegion {
        TypedRegion {
            region: GenomeRegion {
                contig: ContigId(contig),
                start: Position(start),
                end: Position(end),
            },
            kind: RegionKind::Generic,
        }
    }

    /// What one segment does when the scripted generator reaches it.
    #[derive(Clone, Copy)]
    struct SegmentScript {
        /// Loci to emit before the segment ends or fails.
        loci: u32,
        /// Fail instead of ending cleanly, once those loci are out.
        then_fails: bool,
    }

    /// A generic-slot generator that follows a script, one entry per segment begun.
    ///
    /// **Its loci are placed from the segment's own start plus how many it has emitted**, and
    /// each spans three bases (the regions are 1-based and inclusive), so every locus of a walk has a
    /// distinct region and a start that
    /// differs from its end. A test can therefore say exactly which observation it is looking
    /// at, and cannot confuse an end with a start.
    struct ScriptedGenerator {
        remaining_segments: std::vec::IntoIter<SegmentScript>,
        current: Option<SegmentScript>,
        segment: GenomeRegion,
        emitted: u32,
    }

    impl ScriptedGenerator {
        fn following(script: Vec<SegmentScript>) -> Self {
            Self {
                remaining_segments: script.into_iter(),
                current: None,
                segment: GenomeRegion {
                    contig: ContigId(0),
                    start: Position(1),
                    end: Position(1),
                },
                emitted: 0,
            }
        }
    }

    impl LocusGenerator<()> for ScriptedGenerator {
        fn begin_segment(&mut self, region: GenomeRegion) {
            self.segment = region;
            self.current = self.remaining_segments.next();
            self.emitted = 0;
        }

        fn next_locus(
            &mut self,
            _segment: &(),
            _reads: &SampleReads,
        ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
            let Some(script) = self.current else {
                return Ok(None);
            };
            if self.emitted < script.loci {
                let start = self.segment.start.get() + u64::from(self.emitted);
                self.emitted += 1;
                return Ok(Some(SampleLocusObservations {
                    region: GenomeRegion {
                        contig: self.segment.contig,
                        start: Position(start),
                        end: Position(start + 2),
                    },
                    reference_bases: Box::from(&b"AAA"[..]),
                    observations: Vec::new(),
                    reads_without_observation: 0,
                    reads_discarded_by_cap: 0,
                    kind: LocusKind::Generic,
                }));
            }
            if script.then_fails {
                // **A stand-in failure, chosen because it needs nothing but a region.** The
                // walker treats every `LocusGenerationError` alike — it wraps it, names the
                // sample and says how far the walk had got — so which one a fixture raises
                // decides only what the rendered chain says underneath.
                return Err(LocusGenerationError::ForeignSample {
                    region: self.segment,
                });
            }
            Ok(None)
        }
    }

    fn generators_following(script: Vec<SegmentScript>) -> GeneratorSet {
        GeneratorSet::new(
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Generator(Box::new(ScriptedGenerator::following(script))),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        )
    }

    fn emits(loci: u32) -> SegmentScript {
        SegmentScript {
            loci,
            then_fails: false,
        }
    }

    fn emits_then_fails(loci: u32) -> SegmentScript {
        SegmentScript {
            loci,
            then_fails: true,
        }
    }

    /// No fixture here scripts more than a handful of loci, so a walk that reaches this many
    /// draws is not walking — it is repeating.
    ///
    /// **Every draining loop in this module is bounded, and the bound is not decoration.** A
    /// walker that handed back the record it was offered instead of dropping it would yield for
    /// ever, and an unbounded loop turns that into a test binary killed by the allocator: signal
    /// 9, no assertion, no test named. Measured — the one mutation that produced it took about
    /// four minutes against about one for every other, and reported nothing. Bounded, it names
    /// the test and says what went wrong.
    const A_WALK_THAT_IS_REPEATING: usize = 100;

    /// Every observation a walker yields, and the error if it stops on one.
    fn drain(
        walker: &mut impl ObservationSource<Error = RunError>,
    ) -> (Vec<GenomeRegion>, Option<RunError>) {
        let mut regions = Vec::new();
        for _ in 0..A_WALK_THAT_IS_REPEATING {
            match walker.next_observation(None) {
                None => return (regions, None),
                Some(Ok(observation)) => regions.push(observation.region),
                Some(Err(error)) => return (regions, Some(error)),
            }
        }
        panic!(
            "this walk yielded {A_WALK_THAT_IS_REPEATING} observations without ending; \
             no fixture here scripts that many, so the walk is repeating itself"
        );
    }

    // -----------------------------------------------------------------
    // The walk itself
    // -----------------------------------------------------------------

    /// **A walker crosses its segments and yields each one's loci, in order.**
    ///
    /// Three segments emitting 2, 0 and 3 loci, the middle one empty and the third on another
    /// contig. A walker that began only its first segment would yield two; one that ignored
    /// segment boundaries could not place the loci on the right contigs.
    #[test]
    fn a_walk_yields_the_loci_of_every_segment_in_genome_order() {
        let (_reference_dir, _bam_dir, reads) = reads_named("NA12878");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![
                generic_segment(0, 10, 20),
                generic_segment(0, 30, 40),
                generic_segment(1, 50, 60),
            ]),
            reads,
            generators_following(vec![emits(2), emits(0), emits(3)]),
        );

        let (regions, failure) = drain(&mut walker);
        assert!(failure.is_none(), "a clean script does not fail");
        assert_eq!(
            regions,
            vec![
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(10),
                    end: Position(12)
                },
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(11),
                    end: Position(13)
                },
                GenomeRegion {
                    contig: ContigId(1),
                    start: Position(50),
                    end: Position(52)
                },
                GenomeRegion {
                    contig: ContigId(1),
                    start: Position(51),
                    end: Position(53)
                },
                GenomeRegion {
                    contig: ContigId(1),
                    start: Position(52),
                    end: Position(54)
                },
            ],
        );
    }

    /// **The walk's tally is the generators' own**, so a run report reads one number rather
    /// than counting again. Three segments in, three handled, five loci out — which the walk
    /// above emits and nothing else here does.
    #[test]
    fn the_walk_reports_the_counts_its_generators_kept() {
        let (_reference_dir, _bam_dir, reads) = reads_named("NA12878");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![
                generic_segment(0, 10, 20),
                generic_segment(0, 30, 40),
                generic_segment(1, 50, 60),
            ]),
            reads,
            generators_following(vec![emits(2), emits(0), emits(3)]),
        );
        drain(&mut walker);

        let counts = walker.counts();
        assert_eq!(counts.regions_in, 3);
        assert_eq!(counts.regions_handled, 3);
        assert_eq!(counts.loci_emitted, 5);
    }

    // -----------------------------------------------------------------
    // How far it had got
    // -----------------------------------------------------------------

    /// **Nothing has been reached before the first draw.**
    #[test]
    fn a_fresh_walker_has_reached_nothing() {
        let (_reference_dir, _bam_dir, reads) = reads_named("NA12878");
        let walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            reads,
            generators_following(vec![emits(1)]),
        );
        assert_eq!(walker.reached(), WalkProgress::NothingYet);
    }

    /// **Where it reached is the *end* of the *last* observation, on that observation's own
    /// contig.**
    ///
    /// The walk finishes on contig 1 at a locus spanning 51–53, and every other candidate the
    /// code could have taken is a different number: the start of that locus (51), the first
    /// observation of the walk (contig 0, 10–12), and the last segment's own end (60).
    #[test]
    fn where_a_walk_reached_is_the_end_of_its_last_observation() {
        let (_reference_dir, _bam_dir, reads) = reads_named("NA12878");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20), generic_segment(1, 50, 60)]),
            reads,
            generators_following(vec![emits(1), emits(2)]),
        );
        drain(&mut walker);

        assert_eq!(
            walker.reached(),
            WalkProgress::After(GenomePosition {
                contig: ContigId(1),
                position: Position(53),
            }),
        );
    }

    /// **It advances with the walk rather than only at the end**, so a failure part-way names
    /// where the walk was and not where it started.
    #[test]
    fn where_a_walk_reached_advances_observation_by_observation() {
        let (_reference_dir, _bam_dir, reads) = reads_named("NA12878");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            reads,
            generators_following(vec![emits(3)]),
        );

        // Bounded for the reason `drain` states: an endless walk must name a test, not the
        // allocator.
        let mut reached_after_each = Vec::new();
        for _ in 0..A_WALK_THAT_IS_REPEATING {
            let Some(next) = walker.next_observation(None) else {
                break;
            };
            next.expect("a clean script does not fail");
            reached_after_each.push(walker.reached());
        }
        assert!(
            reached_after_each.len() < A_WALK_THAT_IS_REPEATING,
            "this walk is repeating itself rather than ending"
        );

        assert_eq!(
            reached_after_each,
            vec![
                WalkProgress::After(GenomePosition {
                    contig: ContigId(0),
                    position: Position(12)
                }),
                WalkProgress::After(GenomePosition {
                    contig: ContigId(0),
                    position: Position(13)
                }),
                WalkProgress::After(GenomePosition {
                    contig: ContigId(0),
                    position: Position(14)
                }),
            ],
        );
    }

    // -----------------------------------------------------------------
    // Failure
    // -----------------------------------------------------------------

    /// **A failure names this walker's sample and how far this walk had got.**
    ///
    /// Two samples are walked to make the name a claim rather than a coincidence: with one, a
    /// hard-coded name or a name read from the wrong walker would pass. The second sample's
    /// walk also fails at a different place, so the position is its own too.
    #[test]
    fn a_failure_names_the_sample_that_failed_and_where_it_had_reached() {
        let (_reference_dir, _bam_dir, first) = reads_named("zeta");
        let (_other_reference_dir, _other_bam_dir, second) = reads_named("alpha");

        let mut walking_zeta = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            first,
            generators_following(vec![emits_then_fails(2)]),
        );
        let mut walking_alpha = AlignmentFilesWalker::new(
            stream(vec![generic_segment(1, 70, 80)]),
            second,
            generators_following(vec![emits_then_fails(1)]),
        );

        let (_, zetas_failure) = drain(&mut walking_zeta);
        let (_, alphas_failure) = drain(&mut walking_alpha);

        let zeta = zetas_failure.expect("the script fails after two loci");
        assert_eq!(
            zeta.to_string(),
            "sample zeta: reading its observations failed; \
             its last complete observation ended at contig 0 position 13",
        );
        let alpha = alphas_failure.expect("the script fails after one locus");
        assert_eq!(
            alpha.to_string(),
            "sample alpha: reading its observations failed; \
             its last complete observation ended at contig 1 position 72",
        );

        // **The cause says what went wrong**, which the top line deliberately does not: it says
        // which sample and where, and the chain says the rest.
        let chain = format_error_chain(&zeta);
        assert!(
            chain.contains("opened for another sample"),
            "the cause reaches the rendered chain: {chain}",
        );
    }

    /// **A failure on the first draw says so rather than naming a position it never reached.**
    ///
    /// An operator sent to contig 0:1 by a walk that had read nothing would be looking at an
    /// innocent locus.
    #[test]
    fn a_failure_before_any_observation_names_no_position() {
        let (_reference_dir, _bam_dir, reads) = reads_named("zeta");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            reads,
            generators_following(vec![emits_then_fails(0)]),
        );

        let (regions, failure) = drain(&mut walker);
        assert!(regions.is_empty(), "nothing was yielded");
        assert_eq!(
            failure.expect("the script fails at once").to_string(),
            "sample zeta: reading its observations failed; it had produced no observations yet",
        );
    }

    /// **A walk that has failed is spent**, so a merge that ignored the error and asked again
    /// gets `None` rather than a second walk of ground it has already been given.
    #[test]
    fn a_walk_that_failed_yields_nothing_more() {
        let (_reference_dir, _bam_dir, reads) = reads_named("zeta");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            reads,
            generators_following(vec![emits_then_fails(1)]),
        );

        let (_, failure) = drain(&mut walker);
        assert!(failure.is_some(), "the script fails");
        assert!(
            walker.next_observation(None).is_none(),
            "asking a failed walk again yields nothing",
        );
    }

    /// **Exhaustion is final.** A source that yielded `Some` after a `None` would be drawn in
    /// behind the merge's window and so silently out of coordinate order.
    #[test]
    fn a_walk_that_answered_none_keeps_answering_none() {
        let (_reference_dir, _bam_dir, reads) = reads_named("zeta");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            reads,
            generators_following(vec![emits(1)]),
        );

        drain(&mut walker);
        for _ in 0..3 {
            assert!(walker.next_observation(None).is_none());
        }
    }

    // -----------------------------------------------------------------
    // The spare
    // -----------------------------------------------------------------

    /// **A record offered back does not come out again as the next observation.**
    ///
    /// This walker drops the spare (G1 refills it instead), and what must stay true either way
    /// is that what comes back is the *next* observation. A spare handed back unchanged would
    /// put the merge at a position the sample never reported — here, an unmistakable one on a
    /// contig neither segment is on.
    #[test]
    fn a_spare_offered_back_never_comes_out_as_an_observation() {
        let (_reference_dir, _bam_dir, reads) = reads_named("zeta");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            reads,
            generators_following(vec![emits(2)]),
        );

        let spare = || SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(1),
                start: Position(199),
                end: Position(200),
            },
            reference_bases: Box::from(&b"CG"[..]),
            observations: Vec::new(),
            reads_without_observation: 7,
            reads_discarded_by_cap: 7,
            kind: LocusKind::Generic,
        };

        // **Bounded, and this is the loop that made the bound necessary.** A walker that handed
        // the offered record straight back would satisfy the `while let` for ever; measured, that
        // mutation killed the test binary with signal 9 after about four minutes and named
        // nothing.
        let mut regions = Vec::new();
        for _ in 0..A_WALK_THAT_IS_REPEATING {
            let Some(next) = walker.next_observation(Some(spare())) else {
                break;
            };
            regions.push(next.expect("a clean script does not fail").region);
        }
        assert!(
            regions.len() < A_WALK_THAT_IS_REPEATING,
            "this walk kept yielding: the spare it was offered is coming back out"
        );

        assert_eq!(
            regions,
            vec![
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(10),
                    end: Position(12)
                },
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(11),
                    end: Position(13)
                },
            ],
            "the walk's own observations, not the record it was handed",
        );
    }

    // -----------------------------------------------------------------
    // The run's segments
    // -----------------------------------------------------------------

    fn segmentation_over(segments: Vec<TypedRegion>) -> Segmentation {
        let bounds = [
            ContigBounds {
                name: "chr1",
                length: 100,
            },
            ContigBounds {
                name: "chr2",
                length: 200,
            },
        ];
        Segmentation::build(
            segments.into_iter().map(Ok),
            GenomeRegions::whole_contigs(&bounds),
            RepeatCatalogHeader {
                contigs: Vec::new(),
                reference_md5: [7; 16],
                built_under: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "test".to_string(),
                longest_tract_bp: Vec::new(),
            },
            StrRepeatCriteria::default(),
            PathBuf::from("/genomes/test.catalog.parquet"),
        )
        .expect("a clean stream builds")
    }

    /// **Every segment of the run's ground reaches the walk, once, in genome order.**
    ///
    /// A stream that dropped the first, stopped at a contig change, or handed one out twice
    /// would leave a stretch of the genome unwalked or double-counted, and nothing downstream
    /// would say so.
    #[test]
    fn every_run_segment_is_handed_out_exactly_once_in_genome_order() {
        let segments = vec![
            generic_segment(0, 1, 10),
            generic_segment(0, 40, 50),
            generic_segment(1, 5, 9),
        ];
        let segmentation = segmentation_over(segments.clone());

        let handed_out: Vec<TypedRegion> = RunSegments::of(&segmentation)
            .map(|segment| segment.expect("reading a built list cannot fail"))
            .collect();

        assert_eq!(handed_out, segments);
    }

    /// **A walker built over a segmentation walks that segmentation's segments**, which is what
    /// makes "k samples over one segmentation" true rather than "k lists that happen to agree".
    #[test]
    fn a_walker_over_a_segmentation_walks_its_segments() {
        let (_reference_dir, _bam_dir, reads) = reads_named("NA12878");
        let segmentation =
            segmentation_over(vec![generic_segment(0, 10, 20), generic_segment(1, 50, 60)]);
        let mut walker = AlignmentFilesWalker::over(
            &segmentation,
            reads,
            generators_following(vec![emits(1), emits(1)]),
        );

        let (regions, failure) = drain(&mut walker);
        assert!(failure.is_none());
        assert_eq!(
            regions,
            vec![
                GenomeRegion {
                    contig: ContigId(0),
                    start: Position(10),
                    end: Position(12)
                },
                GenomeRegion {
                    contig: ContigId(1),
                    start: Position(50),
                    end: Position(52)
                },
            ],
        );
        assert_eq!(walker.counts().regions_in, 2);
    }

    /// **The sample a walker names is the one whose reads it was given.** Two walkers, two
    /// samples, two names — a name hard-coded or read from anywhere else fails here.
    #[test]
    fn a_walker_names_the_sample_whose_reads_it_holds() {
        let (_reference_dir, _bam_dir, zeta) = reads_named("zeta");
        let (_other_reference_dir, _other_bam_dir, alpha) = reads_named("alpha");
        let segmentation = segmentation_over(vec![generic_segment(0, 10, 20)]);

        let walking_zeta =
            AlignmentFilesWalker::over(&segmentation, zeta, generators_following(vec![emits(1)]));
        let walking_alpha =
            AlignmentFilesWalker::over(&segmentation, alpha, generators_following(vec![emits(1)]));

        assert_eq!(walking_zeta.sample_name(), "zeta");
        assert_eq!(walking_alpha.sample_name(), "alpha");
    }

    /// **A walker's `Debug` says which sample and how far, and does not print its reads.**
    #[test]
    fn a_walkers_debug_names_the_sample_and_not_its_files() {
        let (_reference_dir, _bam_dir, reads) = reads_named("zeta");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![generic_segment(0, 10, 20)]),
            reads,
            generators_following(vec![emits(1)]),
        );
        drain(&mut walker);

        let rendered = format!("{walker:?}");
        assert!(rendered.contains("zeta"), "{rendered}");
        assert!(rendered.contains("After"), "{rendered}");
        assert!(!rendered.contains(".bam"), "{rendered}");
    }
}
