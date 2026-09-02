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
use std::sync::Arc;

use crate::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
use crate::ng::locus_generation::{
    GeneratorSet, GeneratorSlot, LocusCounts, LocusGenerationError, SampleLocusObservations,
    SampleLocusObservationsIterator, UnhandledReason,
};
use crate::ng::read::input::SampleReads;
use crate::ng::read::input::reference::OpenReference;
use crate::ng::read::left_align::LeftAlignPreparer;
use crate::ng::ref_seq::WindowedRefSeq;
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
/// **`Send`, since 2026-09-01, and that is what lets a calling run cover its samples in
/// parallel** — the merge's parallel cover sweeps the cohort's walkers from a worker pool, so
/// a walker must be able to cross a thread (one thread at a time; nothing here is `Sync`, and
/// nothing shares a walker). Three things below this type were widened to get there — two at
/// sites whose own documentation had reserved the change, one (`WindowedRefSeq`) whose doc had
/// only recorded the per-worker ownership that makes the widening safe:
/// [`GeneratorSlot::Generator`](crate::ng::locus_generation::GeneratorSlot)'s trait object
/// gained `+ Send`; the pileup generator's read-preparation cell traded `Rc<RefCell<_>>` for
/// `Arc<Mutex<_>>`; and `WindowedRefSeq`'s resident window traded `RefCell` for `Mutex`, so
/// the `Arc` handles the generator keeps are `Send`. Every lock is uncontended — ownership
/// stays per worker — and `a_run_walker_can_cross_a_thread` is the compile-time proof.
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

impl AlignmentFilesWalker<RunSegments> {
    /// A walker over **the whole run's ground**: every segment of the segmentation, in genome
    /// order, once.
    ///
    /// This is the shape a run builds. Every sample of the run reads the same segments from the
    /// same object — one `Arc::clone` a sample, not one copy of the list — which is what makes
    /// "k samples over one segmentation" true rather than "k segmentations that happen to
    /// agree" (spec §4.2).
    ///
    /// **The type it produces carries no lifetime**, and that is the point rather than a detail:
    /// a run has to hold its walkers and the segmentation they read in one object, and a walker
    /// that borrowed the segmentation could not be stored beside it.
    pub fn over(
        segmentation: Arc<Segmentation>,
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
pub struct RunSegments {
    /// **Shared ownership, not a borrow, and that is what makes a run possible.**
    ///
    /// A run holds one walker per sample for the whole run, and it holds the segmentation those
    /// walkers read. With a borrow the two cannot live in one object: the walkers would borrow a
    /// field of the struct that holds them, which is self-referential and which safe Rust cannot
    /// express. Cloning is not the escape — [`Segmentation`] is deliberately not `Clone`, for the
    /// reason this type exists — and neither is minting a walker per draw, which breaks the
    /// one-source-per-sample rule outright. So the handle is shared, and the genome-sized list is
    /// still stored once however many samples read it.
    ///
    /// **`Arc` rather than `Rc`**, and since 2026-09-01 the choice is load-bearing rather than
    /// insurance: the walker is `Send` and crosses threads under the merge's parallel cover, so
    /// an `Rc` here would take that back. (When this was written the walker was `!Send` through
    /// the generator set's then-unbounded trait object, and the `Arc` was to avoid adding a
    /// second blocker; the blocker was lifted and the insurance paid off.)
    segmentation: Arc<Segmentation>,
    /// How far through the list this stream is. An index rather than a slice iterator because a
    /// borrowing iterator is what the shared handle exists to avoid.
    next: usize,
}

impl RunSegments {
    /// Every segment of `segmentation`, in genome order.
    #[must_use]
    pub fn of(segmentation: Arc<Segmentation>) -> Self {
        Self {
            segmentation,
            next: 0,
        }
    }
}

impl Iterator for RunSegments {
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
        let segment = self.segmentation.segments().get(self.next)?;
        self.next += 1;
        Some(Ok(segment.clone()))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::error_render::format_error_chain;
    use crate::ng::locus_generation::{GeneratorSlot, LocusGenerator, LocusKind, UnhandledReason};
    use crate::ng::read::filtering::ReadFilterConfig;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference_from_its_index, header, indexed_bam, matching_contigs,
        read_named_with_length,
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
        let (reference_dir, reference) = fixture_reference_from_its_index();
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

    /// A segment no generator will take. **Refused permanently rather than as unbuilt**, and the
    /// dispatcher decides that from the kind alone, whatever the slots hold — so it is the one
    /// region kind a scripted fixture can use to fill the out-of-scope counter.
    fn satellite_segment(contig: u32, start: u64, end: u64) -> TypedRegion {
        TypedRegion {
            kind: RegionKind::Satellite,
            ..generic_segment(contig, start, end)
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
                    reference_bases: b"AAA".to_vec(),
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
    /// than counting again. Four segments in, three handled, five loci out.
    ///
    /// **A satellite is in the fixture so the two kinds of nothing land in different counters.**
    /// A region a generator looked at and found nothing in, and a region no generator would take,
    /// are different facts, and with every region handled both refusal counters are zero — a
    /// tally that booked one to the other would read as correct.
    #[test]
    fn the_walk_reports_the_counts_its_generators_kept() {
        let (_reference_dir, _bam_dir, reads) = reads_named("NA12878");
        let mut walker = AlignmentFilesWalker::new(
            stream(vec![
                generic_segment(0, 10, 20),
                generic_segment(0, 30, 40),
                satellite_segment(0, 42, 44),
                generic_segment(1, 50, 60),
            ]),
            reads,
            generators_following(vec![emits(2), emits(0), emits(3)]),
        );
        drain(&mut walker);

        let counts = walker.counts();
        assert_eq!(counts.regions_in, 4);
        assert_eq!(counts.regions_handled, 3, "the three generic segments");
        assert_eq!(counts.loci_emitted, 5);
        assert_eq!(
            counts.unhandled_out_of_scope, 1,
            "the satellite, which is refused permanently"
        );
        assert_eq!(
            counts.unhandled_not_implemented, 0,
            "nothing here is refused as merely unbuilt: the generic slot is filled and no \
             segment asks for another"
        );
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
            reference_bases: b"CG".to_vec(),
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

    fn segmentation_over(segments: Vec<TypedRegion>) -> Arc<Segmentation> {
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
        .map(Arc::new)
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

        let handed_out: Vec<TypedRegion> = RunSegments::of(Arc::clone(&segmentation))
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
            Arc::clone(&segmentation),
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

    /// **A run's walker can cross a thread**, which is what the merge's parallel cover does
    /// with it — each sweep hands every sample's walker to a pool worker, one thread at a
    /// time. This fails at the compiler if anything below the walker loses `Send` again:
    /// the generator slot's trait object, the pileup generator's preparation cell, or the
    /// reference accessor's window (see the type's own note for the three). It holds for
    /// every filling of the generator slots at once, because the slot's box carries the
    /// bound — a generator that is not `Send` is refused where it is boxed, not here.
    #[test]
    fn a_run_walker_can_cross_a_thread() {
        fn assert_send<T: Send>() {}
        assert_send::<super::AlignmentFilesWalker<super::RunSegments>>();
    }

    /// **The sample a walker names is the one whose reads it was given.** Two walkers, two
    /// samples, two names — a name hard-coded or read from anywhere else fails here.
    #[test]
    fn a_walker_names_the_sample_whose_reads_it_holds() {
        let (_reference_dir, _bam_dir, zeta) = reads_named("zeta");
        let (_other_reference_dir, _other_bam_dir, alpha) = reads_named("alpha");
        let segmentation = segmentation_over(vec![generic_segment(0, 10, 20)]);

        let walking_zeta = AlignmentFilesWalker::over(
            Arc::clone(&segmentation),
            zeta,
            generators_following(vec![emits(1)]),
        );
        let walking_alpha = AlignmentFilesWalker::over(
            Arc::clone(&segmentation),
            alpha,
            generators_following(vec![emits(1)]),
        );

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

    /// **A run can hold its walkers and the segmentation they read, in one object** — which is
    /// the whole reason the segments are shared rather than borrowed.
    ///
    /// This is the shape the next step needs and the shape a borrow cannot take: a struct whose
    /// walkers borrowed its own segmentation would be self-referential, and safe Rust would
    /// refuse it. **The test is that this compiles**, so a change back to a borrow fails here at
    /// the compiler rather than three steps later at the wiring.
    ///
    /// **And the list is stored once, not once per sample.** Three walkers over one segmentation
    /// leave four holders of the same allocation — the three of them and the run — where three
    /// copies would leave four allocations of a genome-sized list.
    #[test]
    fn a_run_can_hold_its_walkers_beside_the_segmentation_they_read() {
        struct ARunHoldingBoth {
            segmentation: Arc<Segmentation>,
            walkers: Vec<AlignmentFilesWalker<RunSegments>>,
        }

        let samples = ["zeta", "alpha", "mu"];
        let opened: Vec<_> = samples.iter().map(|name| reads_named(name)).collect();
        let segmentation = segmentation_over(vec![generic_segment(0, 10, 20)]);

        let mut run = ARunHoldingBoth {
            walkers: opened
                .into_iter()
                .map(|(_reference_dir, _bam_dir, reads)| {
                    AlignmentFilesWalker::over(
                        Arc::clone(&segmentation),
                        reads,
                        generators_following(vec![emits(1)]),
                    )
                })
                .collect(),
            segmentation,
        };

        assert_eq!(
            Arc::strong_count(&run.segmentation),
            4,
            "the run and its three walkers, all on one list"
        );
        assert_eq!(
            run.walkers
                .iter()
                .map(AlignmentFilesWalker::sample_name)
                .collect::<Vec<_>>(),
            samples,
        );
        for walker in &mut run.walkers {
            let (regions, failure) = drain(walker);
            assert!(failure.is_none());
            assert_eq!(regions.len(), 1, "each walker walks the shared segment");
        }
    }
}

/// **A source yields exactly what the walk yields** — the observations-equal-the-walk oracle of
/// `doc/devel/ng/spec/run_streaming.md` §12, built as step B2 of
/// `doc/devel/ng/impl_plan/run_driver_direct_mode.md`.
///
/// The tests above prove the *adapter*: the ordering, how far the walk reached, the failure, the
/// segment stream. They cannot prove the *walk*, because every one of them drives a scripted
/// generator that never touches the reads it was handed. This module closes that: it runs the
/// real generic locus generator over a real indexed BAM, twice — once through
/// [`SampleLocusObservationsIterator`] directly, which is the machinery that existed before this
/// step, and once through [`AlignmentFilesWalker`] behind the merge's trait — and compares the
/// two, observation for observation.
///
/// **The oracle is the iterator, not another walker**, which is what makes this a differential
/// rather than a self-check. If the two ever disagree, the walker has changed what a sample
/// reports, and no test above would say so.
#[cfg(test)]
mod walking_the_real_generator {
    use super::*;
    use crate::ng::locus_generation::pileup::{PileupGenerator, PileupGeneratorConfig};
    use crate::ng::locus_generation::{
        GeneratorCounts, GeneratorSlot, LocusCounts, UnhandledReason,
    };
    use crate::ng::read::filtering::ReadFilterConfig;
    use crate::ng::read::input::test_fixtures::{
        fixture_reference_bases, fixture_reference_from_its_index, header, indexed_bam,
        matching_contigs, read_named_with_length,
    };
    use crate::ng::read::left_align::LeftAlignPreparer;
    use crate::ng::region_typing::RegionKind;
    use crate::ng::types::{ContigId, GenomeRegion, Position};
    use std::sync::Arc;

    use noodles_sam::alignment::RecordBuf;
    use noodles_sam::alignment::record_buf::Sequence;

    /// A read of `bases` at `start`, so the walk has something other than the reference to
    /// report.
    ///
    /// **The fixture reference is all `A`** ([`fixture_reference_bases`]), so a read that
    /// carried its default sequence would agree with it everywhere. Giving a read `C`s puts
    /// distinguishable evidence in the records the two walks are compared on: against an all-`A`
    /// reference every one of the 62 loci would otherwise hold a single observation whose bases
    /// equal the reference, and a defect that mangled `bases` or dropped an observation would be
    /// invisible.
    fn read_of(qname: &str, contig: usize, start: usize, bases: &[u8]) -> RecordBuf {
        let mut record = read_named_with_length(qname, contig, start, bases.len());
        *record.sequence_mut() = Sequence::from(bases.to_vec());
        record
    }

    /// The reads both walks see: on `chr1` one that matches the reference everywhere and one
    /// carrying a `C` at two positions, so a locus has both a matching and a non-matching
    /// witness; on `chr2` one more, so a walk that stopped at the first contig change is visible.
    ///
    /// **Thirty bases each, because the shipped read filter drops anything shorter**
    /// (`DEFAULT_MIN_READ_LENGTH`, 30). A fixture of ten-base reads reaches the generator as no
    /// reads at all and every walk comes back empty. That is what the first draft did, and it is
    /// why the locus count below is asserted rather than described.
    fn the_fixture_reads() -> Vec<RecordBuf> {
        vec![
            read_of("chr1-ref", 0, 10, b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"),
            read_of("chr1-alt", 0, 15, b"AACAAAAAAAAAAAAAAAAACAAAAAAAAA"),
            read_of("chr2-alt", 1, 30, b"AAAAACAAAAAAAAAAAAAAAAAAAAAAAA"),
        ]
    }

    /// The ground both walks cover: two generic stretches on `chr1` separated by a satellite that
    /// no generator handles, a generic stretch on `chr2` — and one more on `chr2` that no read
    /// reaches.
    ///
    /// A source that quietly widened a segment, or that handled a region the dispatcher refuses,
    /// would emit loci the iterator does not, and with one uninterrupted segment neither could be
    /// seen. **The last segment is analysed and empty**, which is a different state from
    /// unanalysed and the only one of the five that a filled generator looks at and finds nothing
    /// in.
    fn the_fixture_segments() -> Vec<TypedRegion> {
        vec![
            segment(RegionKind::Generic, 0, 5, 25),
            segment(RegionKind::Satellite, 0, 26, 28),
            segment(RegionKind::Generic, 0, 29, 50),
            segment(RegionKind::Generic, 1, 25, 65),
            segment(RegionKind::Generic, 1, 100, 120),
        ]
    }

    fn segment(kind: RegionKind, contig: u32, start: u64, end: u64) -> TypedRegion {
        TypedRegion {
            region: GenomeRegion {
                contig: ContigId(contig),
                start: Position(start),
                end: Position(end),
            },
            kind,
        }
    }

    /// One sample opened over the fixture reads, and a generator set with the **real** generic
    /// generator in it.
    ///
    /// Built fresh on each call because neither piece can be shared: `SampleReads` is not
    /// `Clone` — deliberately — and a generator carries the state of the walk it has done. So
    /// each arm writes and opens its own copy of the same three records.
    fn a_sample_and_its_generators() -> (
        tempfile::TempDir,
        tempfile::TempDir,
        SampleReads,
        GeneratorSet,
    ) {
        let (reference_dir, reference) = fixture_reference_from_its_index();
        let (bam_dir, bam) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            &the_fixture_reads(),
        );
        let reads =
            SampleReads::open_only_sample(&[bam], &reference, ReadFilterConfig::default(), false)
                .expect("the fixture sample opens");

        let preparer = LeftAlignPreparer::with_default_normalizer(fixture_reference_bases());
        // An ordinary `Arc`: `InMemoryRefSeq` is `Send + Sync`, as `WindowedRefSeq` also is
        // since 2026-09-01. (An `arc_with_non_send_sync` waiver contrast used to live here;
        // the lint stopped firing anywhere when the file-backed accessor became `Sync`.)
        let shared = Arc::new(fixture_reference_bases());
        let generator = PileupGenerator::new(
            shared,
            fixture_reference_bases,
            preparer,
            PileupGeneratorConfig::default(),
        )
        .expect("the generic generator builds against the fixture reference");
        let generators = GeneratorSet::new(
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Generator(Box::new(generator)),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        );
        (reference_dir, bam_dir, reads, generators)
    }

    /// The walk as it existed before this step: the iterator, driven directly.
    fn walked_directly() -> (Vec<SampleLocusObservations>, LocusCounts) {
        let (_reference_dir, _bam_dir, reads, generators) = a_sample_and_its_generators();
        let mut iterator = SampleLocusObservationsIterator::new(
            the_fixture_segments().into_iter().map(Ok::<_, Infallible>),
            reads,
            generators,
        );
        let observations: Vec<_> = (&mut iterator)
            .collect::<Result<_, _>>()
            .expect("the fixture walk succeeds");
        let counts = iterator.counts().clone();
        (observations, counts)
    }

    /// The same ground through the merge's interface. `spare` is what the merge offers back on
    /// every draw, so both the no-reuse and the reuse-offered cases go through one function.
    fn drawn_through_the_source(
        spare: impl Fn() -> Option<SampleLocusObservations>,
    ) -> (Vec<SampleLocusObservations>, LocusCounts) {
        let (_reference_dir, _bam_dir, reads, generators) = a_sample_and_its_generators();
        let mut walker = AlignmentFilesWalker::new(
            the_fixture_segments().into_iter().map(Ok::<_, Infallible>),
            reads,
            generators,
        );
        let mut observations = Vec::new();
        for _ in 0..A_WALK_THAT_IS_REPEATING {
            let Some(next) = walker.next_observation(spare()) else {
                break;
            };
            observations.push(next.expect("the fixture walk succeeds"));
        }
        assert!(
            observations.len() < A_WALK_THAT_IS_REPEATING,
            "this walk is repeating itself rather than ending"
        );
        let counts = walker.counts().clone();
        (observations, counts)
    }

    /// This module's own bound, and a different number from the sibling module's because the
    /// fixtures are: a scripted generator emits a handful of loci where this one emits 62. What
    /// the bound is for is the same — a walk that reaches this many draws is not walking, it is
    /// repeating, and an unbounded loop reports that as a killed test binary rather than as a
    /// named failure.
    const A_WALK_THAT_IS_REPEATING: usize = 1_000;

    /// **The observations a source yields are exactly the observations the walk yields, in the
    /// same order.**
    ///
    /// Compared whole — region, reference bases, every sequence observation and its support, the
    /// reads that witnessed nothing, the reads a cap discarded, and the locus kind — because
    /// `SampleLocusObservations` is `PartialEq` and anything less would let a dropped field
    /// through.
    #[test]
    fn a_source_yields_exactly_what_the_walk_yields() {
        let (walked, _) = walked_directly();
        let (drawn, _) = drawn_through_the_source(|| None);

        // **What the fixture actually covers, asserted rather than described.** Its first draft
        // used ten-base reads, which the shipped filter drops, so every walk came back empty and
        // every comparison below passed on two empty vectors. An asserted count is what says so.
        //
        // 62 is the covered positions of the generic segments, and the regions are 1-based and
        // inclusive: `chr1` 5–25 sees the two reads that start at 10 and 15, so positions 10–25,
        // sixteen; `chr1` 29–50 sees them out to where they end at 39 and 44, positions 29–44,
        // sixteen again; `chr2` 25–65 sees one read over 30–59, thirty. The satellite and the
        // read-free `chr2` segment contribute none.
        assert_eq!(
            walked.len(),
            62,
            "the covered positions of the generic segments"
        );
        // `non_reference_reads` rather than an open-coded comparison: it is the one place the
        // crate writes that test, and it counts only the observations that witnessed the whole
        // locus, which a partial run's shorter bases would otherwise fail against.
        let carrying_a_non_reference_base = walked
            .iter()
            .filter(|locus| locus.non_reference_reads() > 0)
            .count();
        assert_eq!(
            carrying_a_non_reference_base, 3,
            "the three `C`s the fixture reads carry against an all-`A` reference — without them \
             all 62 loci hold one reference-matching observation and the comparison below cannot \
             see a mangled or dropped one"
        );

        assert_eq!(drawn, walked);
    }

    /// **The two walks report the same tallies, not just the same observations.** The
    /// observations could match while the tallies diverged — a region accounted to the wrong
    /// counter emits nothing either way — and the run report reads the tallies.
    #[test]
    fn a_source_and_the_walk_account_for_the_same_regions() {
        let (_, walked) = walked_directly();
        let (_, drawn) = drawn_through_the_source(|| None);

        assert_eq!(drawn, walked);
        assert_eq!(walked.regions_in, 5, "the fixture's five segments");
        assert_eq!(
            walked.unhandled_out_of_scope, 1,
            "the satellite — and out of scope because `begin_region` routes that kind there \
             whatever the slots hold, not because a slot was left unfilled"
        );
        assert_eq!(
            walked.regions_handled, 4,
            "the four generic segments, the read-free one included: a region a filled generator \
             looked at and found nothing in is handled, not unhandled"
        );
    }

    /// One segment walked alone against the same span inside the whole walk — the
    /// segment-independence oracle of `doc/devel/ng/spec/run_streaming.md` §12, compared on
    /// everything but the chain-id numbering (see [`with_chain_ids_renumbered`]).
    ///
    /// **A segment is reached in one of two states, and they are not the same test.** On a
    /// contig the walk has just entered, nothing is carried: the generic generator mints a
    /// cursor, a reference window and a per-chromosome walker afresh at every contig change. On
    /// a contig it is already inside, the cursor has advanced past the earlier segments and the
    /// reference window has released the bases behind it. Both are checked below.
    ///
    /// If either differed, what a sample reports would depend on what ground was walked before
    /// it, and the failure is a stretch of genome missing from the output rather than a crash.
    /// **This is not the thirds-chopping failure spec §4.3 records** — that one cut a segment 74
    /// bases inside a 91-base deletion and lost the bases past the cut. Nothing here cuts a
    /// segment.
    fn walked_alone_matches_the_whole_walk(alone: TypedRegion, expected_loci: usize) {
        // **The whole-walk arm is the iterator's, not another walker's**, so this stays a
        // differential rather than becoming a walker compared with itself. Measured: with both
        // arms drawn through the source, three defects that mangled every yielded observation —
        // cleared chain ids, zeroed support counts, blanked reference bases — passed this test
        // while failing the one above it, because the same defect was applied to both sides.
        let (whole_walk, _) = walked_directly();
        let inside_the_whole_walk: Vec<_> = whole_walk
            .into_iter()
            .filter(|locus| {
                locus.region.contig == alone.region.contig
                    && locus.region.start >= alone.region.start
                    && locus.region.end <= alone.region.end
            })
            .collect();

        let (_reference_dir, _bam_dir, reads, generators) = a_sample_and_its_generators();
        let mut walker = AlignmentFilesWalker::new(
            vec![alone].into_iter().map(Ok::<_, Infallible>),
            reads,
            generators,
        );
        let mut walked_alone = Vec::new();
        while let Some(next) = walker.next_observation(None) {
            walked_alone.push(next.expect("the fixture walk succeeds"));
            assert!(
                walked_alone.len() < A_WALK_THAT_IS_REPEATING,
                "this walk is repeating itself rather than ending"
            );
        }

        assert_eq!(walked_alone.len(), expected_loci);
        assert_eq!(
            with_chain_ids_renumbered(walked_alone),
            with_chain_ids_renumbered(inside_the_whole_walk),
        );
    }

    /// **A segment on a contig the walk has just entered**: `chr2` 25–65, which is reached
    /// fourth inside the whole walk and first when walked alone. Since the generator enters a
    /// contig fresh either way, what this pins is that arriving fourth changes nothing.
    #[test]
    fn a_first_segment_on_a_contig_emits_the_same_loci_alone_and_in_the_whole_walk() {
        walked_alone_matches_the_whole_walk(segment(RegionKind::Generic, 1, 25, 65), 30);
    }

    /// **A segment behind earlier ones on the same contig**: `chr1` 29–50. Inside the whole walk
    /// the cursor has already crossed 5–25 and a satellite, the reference window has released
    /// the bases behind it, and both reads are being met for the second time. **This is where
    /// carried state could change what is emitted**, and the first case cannot see it.
    #[test]
    fn a_later_segment_on_a_walked_contig_emits_the_same_loci_alone_and_in_the_whole_walk() {
        walked_alone_matches_the_whole_walk(segment(RegionKind::Generic, 0, 29, 50), 16);
    }

    /// The same observations with every chain id replaced by the order in which it first
    /// appears — so two walks are compared on **which reads were grouped together**, not on
    /// what those reads were called.
    ///
    /// **Chain ids are walk-relative and the type says so**: "an id names a read within one
    /// walk" ([`SequenceObservation::chain_ids`](crate::ng::locus_generation::SequenceObservation)).
    /// The `chr2` read is id 0 when its segment is walked alone and id 4 when two `chr1`
    /// segments and a satellite came first, because the allocator counts up across a whole walk.
    /// Comparing the raw numbers would therefore fail on a property nobody claims, while
    /// comparing the grouping still catches what matters — a read split into two, two reads
    /// merged into one, or a locus that lost its witnesses.
    ///
    /// **This is the one field of spec §12's fourth oracle that is not literally equal.** The
    /// oracle says a segment walked alone emits *exactly* what it emits inside a whole walk;
    /// measured, everything is equal but this, and this cannot be, by the design of the
    /// allocator. Recorded, not worked around.
    fn with_chain_ids_renumbered(
        mut loci: Vec<SampleLocusObservations>,
    ) -> Vec<SampleLocusObservations> {
        let mut first_seen: Vec<u64> = Vec::new();
        for locus in &mut loci {
            for observed in &mut locus.observations {
                for id in &mut observed.chain_ids {
                    let rank = first_seen
                        .iter()
                        .position(|seen| seen == id)
                        .unwrap_or_else(|| {
                            first_seen.push(*id);
                            first_seen.len() - 1
                        });
                    *id = rank as u64;
                }
            }
        }
        loci
    }

    /// **The spare the merge offers back does not reach the output.**
    ///
    /// Today the walker drops it, so this restates the sibling module's scripted check over the
    /// real generator. It becomes load-bearing at the step that refills the record, where one
    /// stale field carried from a reused record would show up here — against 62 real
    /// observations — and nowhere else.
    #[test]
    fn offering_a_spare_does_not_change_what_a_source_yields() {
        let (walked, _) = walked_directly();
        let (drawn, _) = drawn_through_the_source(|| {
            Some(SampleLocusObservations {
                region: GenomeRegion {
                    contig: ContigId(1),
                    start: Position(199),
                    end: Position(200),
                },
                reference_bases: b"CG".to_vec(),
                observations: Vec::new(),
                reads_without_observation: 7,
                reads_discarded_by_cap: 7,
                kind: crate::ng::locus_generation::LocusKind::Generic,
            })
        });

        // **Its own guard, not one inherited from a sibling test.** Comparing two empty vectors
        // passes, so a fixture that stopped producing anything would leave this green while
        // proving nothing.
        assert_eq!(walked.len(), 62, "the fixture still produces its loci");
        assert_eq!(drawn, walked);
    }

    /// **The generators are reachable through the walker**, so the per-generator counters
    /// [`AlignmentFilesWalker::counts`] does not carry still have a reader once a run owns its
    /// walkers.
    ///
    /// Asserted on a **running, non-zero** count read back through the walker: an accessor wired
    /// to a fresh set, or to a set nothing drove, passes an `is_some` and fails this.
    ///
    /// **Five admissions from three reads, and the difference is the point of the number.**
    /// `reads_admitted` counts admissions, not distinct reads: a read is admitted once per
    /// segment it is met in, and the two `chr1` reads span both `chr1` segments. Asserting three
    /// would be asserting a count this field does not keep.
    #[test]
    fn the_generators_counts_are_reachable_through_the_walker() {
        let (_reference_dir, _bam_dir, reads, generators) = a_sample_and_its_generators();
        let mut walker = AlignmentFilesWalker::new(
            the_fixture_segments().into_iter().map(Ok::<_, Infallible>),
            reads,
            generators,
        );
        while let Some(next) = walker.next_observation(None) {
            next.expect("the fixture walk succeeds");
        }

        let Some(GeneratorCounts::Pileup(counts)) = walker.generators().generic_counts() else {
            panic!("the generic slot is filled, so its counts must be reachable");
        };
        assert_eq!(
            counts.reads_admitted, 5,
            "the `chr2` read once, and each `chr1` read once per `chr1` segment"
        );
        assert!(
            walker.generators().ssr_counts().is_none(),
            "an unfilled slot counts nothing"
        );
    }
}

/// The reference the walk fetches its bases from, opened once for the whole run.
///
/// **Two accessors are held per sample, a third is minted per file per chromosome, and none of
/// them may be shared.** The locus generator keeps one for the walk's own REF fetches and the
/// read preparer keeps a second — neither rebuilt per segment, because a fresh accessor at every
/// boundary throws away the sliding buffer and re-pays the `.fai` parse. The third is a
/// *factory*: each of a sample's files gets its own accessor every time a cursor is made, which
/// is once per file per chromosome. `WindowedRefSeq` holds an open per-contig reader, and the
/// input layer takes a factory because sharing one accessor across cursors would collapse them
/// onto one file position and one window (spec §8's trap). The type stopped *forbidding* the
/// sharing on 2026-09-01, when it became `Sync` so a walker could cross threads — the reason
/// not to share is the window's ownership, not the auto-trait, and it is unchanged.
///
/// **What is shared instead is the index and the contig table.** An accessor that parses the
/// `.fai` for itself costs about 189 µs on a GRCh38-shaped reference of 2,580 records; sharing
/// the index brings that to 52 µs, of which 34 is cloning the contig table; sharing the table
/// too leaves about 18 µs, which is the `open(2)` that giving each reader its own cursor
/// actually means (measured, [`WindowedRefSeq::with_shared_index`]). **And the parse is paid at
/// every contig open, not once per accessor** ([`WindowedRefSeq::new`]) — so on a 63-sample
/// cohort over a dozen contigs, sharing is the difference between one parse and several hundred.
pub(crate) struct WalkReference {
    fasta: std::path::PathBuf,
    contigs: Arc<crate::fasta::ContigList>,
    index: Arc<noodles_fasta::fai::Index>,
}

impl WalkReference {
    /// Open the run's reference for walking, parsing its index once.
    ///
    /// **Refuses a reference that has no bases.** A reference read from a `.fai` alone describes
    /// a genome's geometry and holds no sequence, so the walk has nothing to fetch REF bases
    /// from — and the refusal belongs here rather than at the first locus, where it would arrive
    /// after every file was opened and a genome's worth of setup was done.
    pub(crate) fn of(reference: &OpenReference) -> Result<Self, RunError> {
        let fasta = reference
            .info()
            .fasta_path
            .clone()
            .ok_or(RunError::ReferenceHasNoBases)?;
        let index = WindowedRefSeq::read_index(&fasta).map_err(|source| {
            RunError::ReferenceIndexUnreadable {
                reference: fasta.clone(),
                source,
            }
        })?;
        Ok(Self {
            fasta,
            contigs: Arc::new(reference.info().contig_list()),
            index,
        })
    }

    /// One accessor of its own, over the shared index and contig table.
    #[must_use]
    pub(crate) fn accessor(&self) -> WindowedRefSeq {
        WindowedRefSeq::with_shared_index(
            self.fasta.clone(),
            Arc::clone(&self.contigs),
            Arc::clone(&self.index),
        )
    }
}

/// **The sizes and the path, not the bases.** A derived `Debug` would print a contig table.
impl std::fmt::Debug for WalkReference {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("WalkReference")
            .field("fasta", &self.fasta)
            .field("contigs", &self.contigs.entries.len())
            .finish_non_exhaustive()
    }
}

/// The generator set a run builds today: **the generic path filled, and both repeat-tract slots
/// refused as unbuilt**.
///
/// A segment routed to an unfilled slot emits no locus and is counted against
/// `unhandled_not_implemented`, which is how a run says *this ground was analysed and this
/// caller cannot yet speak for it* — as opposed to the satellite's permanent refusal. **So a run
/// over ground with repeat tracts in it is not wrong, it is short**, and the tally says by how
/// much. Candidate selection at a tract is specified and unbuilt
/// (`doc/devel/ng/impl_plan/candidate_alleles_ssr.md`); this is where the second slot is filled
/// when it exists.
///
/// One set per sample: a locus generator carries state across segments and cannot be shared
/// (spec §8).
pub(crate) fn generic_path_generators(
    reference: &WalkReference,
    config: PileupGeneratorConfig,
) -> Result<GeneratorSet, RunError> {
    let make_reference = {
        let reference = WalkReference {
            fasta: reference.fasta.clone(),
            contigs: Arc::clone(&reference.contigs),
            index: Arc::clone(&reference.index),
        };
        move || reference.accessor()
    };
    // **The `Arc` is the generator's constructor asking for one.** Each of these three
    // accessors is its own — the input layer takes a factory precisely so nothing shares a
    // window — and nothing shares this handle either; the `Arc` is the constructor's shape.
    // *(A `clippy::arc_with_non_send_sync` waiver stood here while `WindowedRefSeq` was
    // `!Sync`; it became `Sync` on 2026-09-01 so a walker can cross threads under the
    // merge's parallel cover, and the lint no longer fires.)*
    let shared = Arc::new(reference.accessor());
    let generator = PileupGenerator::new(
        shared,
        make_reference,
        LeftAlignPreparer::with_default_normalizer(reference.accessor()),
        config,
    )
    .map_err(|source| RunError::LocusGeneratorSettings { source })?;
    Ok(GeneratorSet::new(
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        GeneratorSlot::Generator(Box::new(generator)),
        GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
    ))
}
