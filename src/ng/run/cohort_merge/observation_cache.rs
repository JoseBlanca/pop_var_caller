//! The observation cache — one forward reader per sample, and the window the builders read.
//!
//! Upstream produces a sample's observations in one forward pass; the builders want them in
//! several places at once (`doc/devel/ng/spec/cohort_merge.md` §6.4). Nothing can serve that
//! by seeking, and giving every builder its own reader would mean as many readers per sample
//! as there are builders. So there is **one reader per sample for the whole run**, advancing
//! forward only, and this is what sits between it and the builders: it draws the readers
//! forward until a region's ground is covered, hands the covering observations out for the
//! length of a call, and drops what nothing can reach any more.
//!
//! **A *region* here is a building region** — the stretch of genome assigned to one builder,
//! `cohort_locus_builder_regions_len` bases (spec §6.1). The run's own intervals are called
//! *analysed regions*, as in [`merge_cohort_serially`](super::serial::merge_cohort_serially),
//! and are much longer.
//!
//! **Why the merge needs this rather than merely liking it.** A builder closes loci from the
//! beginning of whatever observations it is handed and discards those that opened before its
//! own ground ([`build_region`](super::build::build_region)), so handing every builder the
//! whole analysed stretch costs each of them the whole prefix — about **3.3 µs per prefix
//! base at 63 samples**, measured in a release build by the C1 review
//! (`doc/devel/reports/reviews/ng_cohort_merge_c1_2026-08-17.md`). The same effect end to
//! end, measured on **one sample** by the C2 review
//! (`doc/devel/reports/reviews/ng_cohort_merge_c2_2026-08-17.md`): 20,000 observations cost
//! **5.4 ms merged as one analysed region and 184 ms as a thousand**. Short building regions
//! are only affordable when each builder is handed a window over its own ground, and that
//! window is what this file produces.
//!
//! **Builders read it and never write it** (spec §6.4, goal 1):
//! [`cover`](ObservationCache::cover) and [`evict_before`](ObservationCache::evict_before)
//! take `&mut self` and belong to the organiser, while
//! [`with_observations`](ObservationCache::with_observations) takes `&self` and is a builder's
//! only way in.
//!
//! **The organiser is [`super::organise`], not here** — the two shared a file until the owner
//! split them at Checkpoint E. The argument for keeping them together had been that the
//! organiser would become the cache's only writer, at which point
//! [`cover`](ObservationCache::cover) and [`evict_before`](ObservationCache::evict_before)
//! could turn private to that file. E1's review found the premise false: the cached serial
//! driver calls both from a sibling module and nothing in the plan removes it. So the
//! reachable narrowing is `pub(super)` — taken here — which a file of its own gets equally,
//! and against it stood one file of two thousand-odd lines covering two types that share no
//! field and no function.

use super::CohortLocusBuilderRegionsLen;
use crate::ng::locus_generation::SampleLocusObservations;
use crate::ng::ref_seq::{ContigTable, EvictableRefSeq, RefSeq, RefSeqError};
use crate::ng::types::{ContigId, GenomePosition, GenomeRegion, Position};
use crate::ng::window_coverage::depth::{EvidenceForOneRecord, for_each_reported_depth};
use crate::ng::window_coverage::{
    self, SampleHistogram, WindowCoverage, WindowCoverageAccumulator, WindowCoverageConfig,
};

/// What the merge asks of a reference: bases, and the release of what it has walked past.
///
/// **Both halves, because a merge walks a whole contig forward.** A windowed reader extends its
/// buffer while each request lands near the last — which is exactly the merge's pattern, one
/// building region after the next — and shrinks it only when it is told to. A reader that is
/// never told ends a contig holding every base the merge passed: about 250 MB on human
/// chromosome 1, against a walk that otherwise peaks near 25 MB
/// ([`ref_seq`](crate::ng::ref_seq)'s `a_forward_walk_holds_the_whole_span_unless_it_releases`).
/// Asking for [`RefSeq`] alone would make releasing unavailable rather than merely unused.
///
/// **And the contig table, because the look-ahead has to stop at a contig's end.** A cover draws
/// half a window past the region it was asked for ([`cover`](ObservationCache::cover)), and the
/// ground it then reads runs that far too — so on the last region of a contig the read would run
/// off the end and the fetch would refuse it. The lengths are the only thing that says where to
/// stop, and every reference this merge is given already carries them.
pub trait MergeReference: RefSeq + EvictableRefSeq + ContigTable {}

impl<T: RefSeq + EvictableRefSeq + ContigTable + ?Sized> MergeReference for T {}

/// The reference bases over ground a cover holds could not be read.
///
/// **A failure of the cache and not of a source**, which is why it is its own type rather than
/// one of a source's errors: the cache is generic over what its readers refuse, and this is
/// something the cache itself hit. A caller converts it into whatever it fails with — a run
/// converts it into [`RunError::WindowCoverageGroundUnreadable`](crate::ng::run::RunError), and
/// the bound that says so sits on [`cover`](ObservationCache::cover) rather than on the type,
/// so a caller that never covers never has to name it.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
#[error("the reference bases over {region} could not be read")]
pub struct ReferenceUnreadable {
    /// The ground the fetch was over: every position of every record the cover holds on that
    /// contig, widened to the region it was asked for.
    pub region: GenomeRegion,
    /// What the reference fetch hit.
    #[source]
    pub source: RefSeqError,
}

/// What the merge knows about one sample's observation **before** it decides to build
/// anything.
///
/// **These are the whole of what closing a cohort locus reads *per observation*.** Where the
/// locus opens and how far it reaches come from [`region`](Self::region); whether any sample
/// varied enough for the locus to be worth calling comes from the two counts, which are the
/// numerator and the denominator of the keep rule
/// ([`MinAltReads::reached_by`](super::MinAltReads::reached_by)). Nothing else about an
/// observation is consulted until the locus has survived both verdicts and is being assembled.
///
/// **The locus's kind is deliberately not here, because it is not a fact about an
/// observation.** The width bound applies to generic loci and not to repeat tracts, whose span
/// the reference fixes — but that verdict is passed once on the closed locus, reading the kind
/// the opening observation carried, not once per member. The only per-member read of a kind is
/// a release assertion that a locus never mixes the two, which is structurally guaranteed by
/// segments being the reference's own partition. Putting a kind in this type would have cost
/// every observation a field to serve a check that cannot fire, and would have obliged a
/// summary drawn from a psp record's head to answer something that record's head does not
/// carry (the owner's ruling of 2026-09-04 keeps the tag in the body).
///
/// **Why that set has a name.** A sample whose observations are stored in a psp can answer all
/// four without decoding the evidence they describe, because a stored record's head carries
/// them (`spec/psp_file_format.md` §4.3) — so a run over stored files can decide which loci are
/// worth building before it builds any of them, and at about one position in a hundred varying,
/// most are never built at all (`spec/cohort_merge_psp_path.md` §2). While the closing walk
/// reached into each record for each fact separately, no source could supply the facts without
/// supplying the record; naming them is what lets one exist.
///
/// **Direct mode derives them from the record it is already holding**, at the cost it already
/// paid — [`of`](Self::of) is the same pair of reads and the same single walk over the
/// sequences that the closing walk made inline before this type existed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocusSummary {
    /// The reference the observation covers, first base to last.
    pub region: GenomeRegion,
    /// The reads that showed something other than the reference — the keep rule's numerator.
    pub non_reference_reads: u32,
    /// The reads whose whole sequence over the locus was compared against the reference — the
    /// keep rule's denominator, and neither read depth nor the reads that merely covered the
    /// ground.
    pub reads_compared_with_reference: u32,
}

impl LocusSummary {
    /// The summary of a record already in hand — direct mode's only path, and the psp path's
    /// oracle.
    ///
    /// **Both counts come from one walk of the record's sequences**, which is why they are
    /// taken together here rather than read one at a time: they filter on the same witness and
    /// read the same support, so asking separately walks the observations twice. What that
    /// saves grows with the reads at a position rather than with the cohort.
    #[must_use]
    pub fn of(observation: &SampleLocusObservations) -> Self {
        let (non_reference_reads, reads_compared_with_reference) =
            observation.non_reference_and_compared_reads();
        Self {
            region: observation.region,
            non_reference_reads,
            reads_compared_with_reference,
        }
    }

    /// Where the observation begins, as a whole-genome position.
    #[must_use]
    pub fn start_position(self) -> GenomePosition {
        GenomePosition {
            contig: self.region.contig,
            position: self.region.start.min(self.region.end),
        }
    }

    /// The last reference position the observation covers, as a whole-genome position.
    #[must_use]
    pub fn reach_position(self) -> GenomePosition {
        GenomePosition {
            contig: self.region.contig,
            position: self.reach(),
        }
    }

    /// The last reference position the observation covers.
    ///
    /// Saturating for the reason [`SampleLocusObservations::reach`] is: a region ending at the
    /// coordinate ceiling must not wrap a comparison the closing walk makes on every record.
    #[must_use]
    pub fn reach(self) -> Position {
        self.region.end.max(self.region.start)
    }
}

/// One sample's observations in coordinate order, **and the place a record goes when the
/// merge has finished with it**.
///
/// **The second half is why this is a trait and not an `Iterator` bound.** The merge frees
/// far more than it allocates — 6.4 million blocks a round against 216 thousand on the
/// tomato panel, counted — because the records it walks were allocated by the stage upstream
/// and are released as it passes them. Measured by making the merge leak instead of free,
/// that is **63% of the eight-thread merge** and 23% of the single-threaded one. No
/// scheduling change reaches it: the work is real, it is just work nobody needs done.
///
/// So a source is offered its own spent records back. A source that mints them — the
/// generator, or a psp reader decoding into buffers — can fill the record it is handed
/// instead of allocating a new one, and then neither side allocates per position after the
/// first window. A source that cannot simply ignores the offer, which is what the blanket
/// implementation below does for every plain iterator, so nothing that exists today has to
/// change.
///
/// **The spare is an offer and not an obligation**, and that is the whole of the contract: a
/// source may fill it, drop it, or keep it for later, and the record it returns need have
/// nothing to do with it. Making it an obligation would mean a source could not decide per
/// record whether reuse is possible — which a decoder must, since a record whose buffers are
/// the wrong size is cheaper to allocate than to reshape.
/// What a source hands over when the cache draws from it.
///
/// **Two shapes, because two sources answer differently and the merge must not care.** A walker
/// over alignment files has just minted the record and gives it whole; a reader over a stored
/// sample gives the summary its head carried and leaves the evidence in its own arena, to be
/// built only if a locus survives (`spec/cohort_merge_psp_path.md` §3.1). At about one position
/// in a hundred surviving, the second shape is what makes the other ninety-nine free.
pub enum Drawn {
    /// The record itself, already built.
    Built(SampleLocusObservations),
    /// Its summary, with the evidence kept by the source under `body`.
    Kept {
        /// Everything the merge decides on before assembling.
        summary: LocusSummary,
        /// Where the source is holding the evidence — meaningful only to that source.
        body: core::ops::Range<usize>,
    },
}

impl Drawn {
    /// What the merge decides on, whichever shape arrived.
    #[must_use]
    pub fn summary(&self) -> LocusSummary {
        match self {
            Self::Built(record) => LocusSummary::of(record),
            Self::Kept { summary, .. } => *summary,
        }
    }
}

pub trait ObservationSource {
    /// What a failed read is. The cache adds nothing to it and passes it through, so it must
    /// name the sample it came from (arch §5).
    type Error;

    /// The next observation, or `None` once this sample is spent.
    ///
    /// `spare` is a record the merge will not read again, offered for reuse. It is `None`
    /// when the cache has none to hand back.
    ///
    /// **Once this answers `None` it is never called again**, which the cache guarantees with
    /// a flag of its own — the guard [`Fuse`](std::iter::Fuse) used to give, and the reason it
    /// matters is unchanged: a source that yielded `Some` after a `None` would be drawn in
    /// behind the window's own right edge and so silently out of coordinate order. A
    /// *failure* is `Some(Err(_))` and leaves the source live, which is what lets a cover be
    /// made again.
    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, Self::Error>>;

    /// The next observation in whichever shape this source gives, **and the one the cache
    /// actually calls**.
    ///
    /// The default builds every record, which is what a walker over alignment files does
    /// anyway and what every fixture in this module expects. A source that can answer more
    /// cheaply — a reader over a stored sample, which has the summary in the record's head and
    /// need not decode the body behind it — overrides this and keeps its evidence, at the cost
    /// of implementing [`build`](Self::build).
    fn next_drawn(&mut self, spare: Option<SampleLocusObservations>) -> Option<Result<Drawn, Self::Error>> {
        Some(self.next_observation(spare)?.map(Drawn::Built))
    }

    /// Build the evidence this source kept at `body`.
    ///
    /// **Called from a shared reference, on whichever thread got there, in any order** — which
    /// a psp reader can serve because a stored body is decoded from its own bytes alone
    /// (`spec/cohort_merge_psp_path.md` §3.2). A source whose [`next_drawn`](Self::next_drawn)
    /// only ever returns [`Drawn::Built`] is never asked, and the default says so.
    ///
    /// # Errors
    ///
    /// Whatever building the evidence refuses.
    fn build(&self, body: core::ops::Range<usize>) -> Result<SampleLocusObservations, Self::Error> {
        unreachable!(
            "a source that keeps no evidence was asked to build bytes {body:?}: only a source \
             whose `next_drawn` returns `Drawn::Kept` is ever asked, and this one does not"
        )
    }
}

/// Every iterator of one sample's observations is a source that does not reuse.
///
/// This is what keeps the trait from being a migration: the fixtures, the probes and the
/// direct path all hand the cache a plain iterator and go on working, paying exactly the
/// allocator traffic they paid before. Reuse is something a source opts into by implementing
/// the trait itself.
impl<I, E> ObservationSource for I
where
    I: Iterator<Item = Result<SampleLocusObservations, E>>,
{
    type Error = E;

    fn next_observation(
        &mut self,
        spare: Option<SampleLocusObservations>,
    ) -> Option<Result<SampleLocusObservations, E>> {
        drop(spare);
        self.next()
    }
}

/// Every sample's observations over the ground currently assigned to builders (spec §6.4).
///
/// **One forward reader per sample, and a window of what it has drawn.** The cache never
/// seeks and never goes back: the only ways the window moves are [`cover`](Self::cover),
/// which draws it forward, and [`evict_before`](Self::evict_before), which drops from its
/// left edge.
///
/// **This is the module's dominant memory** (spec §8), and it is what a building region's
/// width is paid for in: the ground it spans is `builders × cohort_locus_builder_regions_len`
/// — 3,200 bases at 16 builders on the default 200-base regions — plus the tail of
/// observations reaching past it.
///
/// **A source is any iterator of one sample's observations in coordinate order, and its error
/// type is its own.** The architecture writes `cover(…) -> Result<(), RunError>` and names
/// `run_streaming.md` arch §2's `ObservationSource` as what the cache is handed; neither
/// exists yet — both belong to the run's own document, which is out of this plan's scope — so
/// the cache is generic over the source and passes its failures through untouched. When
/// `ObservationSource` lands, its `observations_in` iterator is exactly this shape and
/// `RunError` is exactly this `E`.
///
/// **`E` must identify the sample it came from.** This cache knows which sample's reader
/// failed and passes the error on without adding to it, so a source whose error does not name
/// its own sample produces a failure an operator cannot act on at a cohort of thousands —
/// arch §5's two variants both carry `sample` for exactly that reason.
pub struct ObservationCache<S> {
    /// One per sample, in the run's sample order — the order every consumer indexes by.
    samples: Vec<SampleWindow<S>>,
    /// **One reference accessor for the whole run** (`spec/window_coverage.md` §3.2). It slides
    /// forward with the merge and releases what the merge has walked past — the release is
    /// [`fetch_the_ground_this_cover_reached`](Self::fetch_the_ground_this_cover_reached)'s last
    /// statement, and [`MergeReference`] is the bound that makes it available at all. It is
    /// never shared with another thread; `Send + Sync` is what keeps the cache itself `Send`
    /// where its sources are, and what lets the parallel cover borrow the cache shared.
    ///
    /// A trait object rather than a type parameter: the accessor is read **once per cover**
    /// and not once per position, so the indirection costs nothing measurable, and every
    /// signature that names this cache would otherwise gain a second parameter for it — while
    /// the fixtures, which want a hand-built reference rather than a FASTA on disk, would gain
    /// it too.
    reference: Box<dyn MergeReference + Send + Sync>,
    /// The reference bases over the ground the last [`cover`](Self::cover) reached, fetched
    /// once and read by every sample at its own offset.
    ///
    /// **A buffer the cover refills rather than a per-position fetch**: the window-coverage
    /// measurement needs a base at every position of every sample's records, and asking the
    /// reference for each one would be one file read per position per sample.
    reference_bases: Vec<u8>,
    /// Which position [`reference_bases`](Self::reference_bases)`[0]` is, or `None` before the
    /// first cover. Read through [`reference_base_at`](Self::reference_base_at), which is the
    /// only thing that should do the arithmetic.
    reference_bases_from: Option<GenomePosition>,
    /// **Whether any source has kept its evidence rather than building it**, latched on the
    /// first such draw and never cleared.
    ///
    /// A window's records cannot be inferred from a sample's holding none: a sample that has
    /// drawn nothing yet, or whose window was just emptied, holds no records either way. Read
    /// that as "records present" and a builder indexes an empty slice for a locus whose
    /// members came from summaries — which is a panic, and was one.
    keeps_evidence: bool,
    /// How far a **successful** [`cover`](Self::cover) has drawn, genome-wide.
    ///
    /// A cover that failed does not move it, which is what lets
    /// [`with_observations`](Self::with_observations) refuse ground no reader reached rather
    /// than hand out a window that is short — and short is a locus closed over the wrong
    /// ground, which is a wrong answer rather than a failure.
    covered_to: Option<GenomePosition>,
}

/// One sample's reader and the observations drawn from it that have not been evicted.
struct SampleWindow<S> {
    /// The forward reader. Never seeks, never rewinds.
    source: S,
    /// Whether the source has already answered `None`.
    ///
    /// **This is the guard [`Fuse`](std::iter::Fuse) used to be**, kept because the reason
    /// for it is unchanged and moved here because a source is now a trait rather than an
    /// iterator: one that yielded `Some` after a `None` would be drawn in behind the window's
    /// own right edge, and so silently out of coordinate order. Set only on `None`, so a
    /// *failure* — `Some(Err(_))` — leaves the source live and a cover can be made again.
    spent: bool,
    /// Records this sample's window has finished with, offered back to its source.
    ///
    /// **Capped at what the sample currently holds**, which is the round's own ground and is
    /// the term spec §8 already prices — so recycling cannot make the cache's memory a
    /// different shape from the one that section bounds. What does not fit is freed as
    /// before.
    spare: Vec<SampleLocusObservations>,
    /// What has been drawn and not yet evicted, in coordinate order.
    ///
    /// **A `Vec` and not a `VecDeque`**, because a builder is handed a contiguous slice of it
    /// ([`ObservationCache::with_observations`]) and a deque's two halves are not one.
    /// Eviction pays a move of what survives, which is the window — short by construction.
    held_observations: Vec<SampleLocusObservations>,
    /// One summary per held observation, at the same index, drawn and evicted with it.
    ///
    /// **Kept beside the records rather than derived when wanted, because deriving one walks
    /// the record's sequences** and the closing walk wants every summary in its window on
    /// every pass. Held here, that walk happens once per record drawn instead of once per
    /// record per cover.
    ///
    /// **And because a summary is what a run over stored files will have *instead* of a
    /// record.** A psp answers all of it from a record's head without decoding the evidence
    /// behind it, so this is the array that path fills while `held_observations` stays empty
    /// until a locus survives (`spec/cohort_merge_psp_path.md` §3.1). Direct mode fills both,
    /// which is why the two are parallel arrays rather than one array of pairs: only one of
    /// them is always present.
    held_summaries: Vec<LocusSummary>,
    /// Where each held observation's evidence sits in its source, when the source kept it
    /// rather than building it. Empty for a source that builds.
    held_bodies: Vec<core::ops::Range<usize>>,
    /// Whether this source has ever kept a body — latched, so an emptied window still says so.
    keeps_evidence: bool,
    /// Where the last observation drawn from `source` began — the ordering check's memory.
    last_drawn: Option<GenomePosition>,
    /// This sample's window-coverage measurement, at whatever stage the pass has reached it.
    ///
    /// **One per sample and never shared**, because a window is over *that sample's* covered
    /// positions: the GC is averaged over the positions that contributed the depth, so two
    /// samples covering different ground have different windows at the same centre.
    window_coverage: WindowCoverageInProgress,
}

/// One sample's window-coverage measurement while the pass is still running: the accumulator
/// being fed, the windows it has finalised, and how far it has read.
///
/// **The three are one thing and move together**, which is why they are a type rather than
/// three fields on [`SampleWindow`]: the cursor is only meaningful against the accumulator it
/// guards, eviction touches only the finalised windows, and the step that finishes the pass
/// (plan step C5) moves the whole of this out by value.
struct WindowCoverageInProgress {
    /// Every covered position this sample has a record at, fed in coordinate order, giving
    /// back the mean depth and GC of the window centred on each
    /// (`spec/window_coverage.md` §3.3).
    accumulator: WindowCoverageAccumulator,
    /// The windows the accumulator has finalised and nothing has evicted, oldest first — the
    /// values a builder reads at a locus.
    ///
    /// **A `Vec` drained from the front, for the reason `held_summaries` is one**: a consumer
    /// is handed a contiguous slice of it, and a deque's two halves are not one.
    finalised: Vec<(GenomePosition, WindowCoverage)>,
    /// Where the record the accumulator last saw began, or `None` before the first.
    ///
    /// **The cursor that makes a record observed exactly once.** A record is held across as
    /// many covers as it reaches into, and every one of those covers reads ground that contains
    /// it; without this, its positions would be counted once per cover and the sample's depth
    /// would rise with how many covers happened to overlap it. Records reach the accumulator in
    /// ascending order and never repeat a start — a sample's records are disjoint and ascending
    /// — so a start is enough to say what has been seen.
    observed_through: Option<GenomePosition>,
}

impl WindowCoverageInProgress {
    /// A measurement that has seen nothing, configured from the run's constants.
    ///
    /// **The configuration is spelled here and nowhere else**, which is what
    /// `WindowCoverageConfig`'s own doc asks for: it has no `Default`, so the run's values live
    /// at the one site that builds an accumulator.
    fn new() -> Self {
        Self {
            accumulator: WindowCoverageAccumulator::new(WindowCoverageConfig {
                window_bp: window_coverage::WINDOW_BP,
                gc_bins: window_coverage::GC_BINS,
                depth_bins: window_coverage::DEPTH_BINS,
                depth_scale_windows: window_coverage::DEPTH_SCALE_WINDOWS,
                depth_range_in_medians: window_coverage::DEPTH_RANGE_IN_MEDIANS,
                min_window_positions: window_coverage::MIN_WINDOW_POSITIONS,
            }),
            finalised: Vec::new(),
            observed_through: None,
        }
    }

    /// Where in `held` the records this has not seen begin.
    ///
    /// A sample's records are disjoint and ascending, so what has not been seen is a suffix and
    /// where it starts is a partition rather than a scan of the whole window.
    fn first_not_yet_observed(&self, held: &[LocusSummary]) -> usize {
        self.observed_through.map_or(0, |seen| {
            held.partition_point(|summary| summary.start_position() <= seen)
        })
    }

    /// The summaries in `held` this has not seen.
    fn summaries_not_yet_observed<'h>(&self, held: &'h [LocusSummary]) -> &'h [LocusSummary] {
        &held[self.first_not_yet_observed(held)..]
    }

    /// Feed one record's reported positions to the accumulator, against `bases`, and keep
    /// whatever windows that finalises.
    ///
    /// `summary` must be `evidence`'s own, and `bases` must hold every position the record
    /// reports — the ground a cover reads is computed from the very records it then observes
    /// ([`ObservationCache::ground_this_cover_holds_on`]), so a position with no base is a
    /// defect in that computation rather than an absent value, and it says so.
    ///
    /// # Errors
    ///
    /// Whatever building a kept record's evidence refuses, unchanged.
    fn observe<E>(
        &mut self,
        summary: LocusSummary,
        evidence: EvidenceForOneRecord<'_>,
        build: impl FnOnce(core::ops::Range<usize>) -> Result<SampleLocusObservations, E>,
        bases: &[u8],
        bases_from: GenomePosition,
    ) -> Result<(), E> {
        let Self {
            accumulator,
            finalised,
            observed_through,
        } = self;
        for_each_reported_depth(summary, evidence, build, |at, depth| {
            let base = base_in(bases, bases_from, at).unwrap_or_else(|| {
                panic!(
                    "the cover read no reference base at {at:?}, where this sample holds a \
                     record — the ground a cover reads is computed from the records it holds, \
                     so this is that computation and not the data"
                )
            });
            accumulator.observe(at.contig, at.position, base, depth);
        })?;
        *observed_through = Some(summary.start_position());
        while let Some(window) = accumulator.pop_ready() {
            finalised.push(window);
        }
        Ok(())
    }

    /// Drop the finalised windows centred before `position`.
    ///
    /// **By the windows' own coordinates and not by the records'**: a window is centred on a
    /// position, and the record that position came from may already be gone. Dropped ones can
    /// never be asked for again — a builder reads a window at the locus it is building, and the
    /// organiser releases ground only once nothing can reach back into it.
    ///
    /// **A prefix drain, for the reason the summaries' is one**: the windows are ascending, so
    /// the survivors are a suffix, and the cost is proportional to what is dropped rather than
    /// to the window.
    fn evict_before(&mut self, position: GenomePosition) {
        let first_survivor = self
            .finalised
            .partition_point(|(centre, _)| *centre < position);
        self.finalised.drain(..first_survivor);
    }
}

impl<S> SampleWindow<S> {
    /// Drop everything this sample holds that ends before `position`, offering the records back
    /// to its own spare list and pushing what will not fit onto `dead`.
    ///
    /// **One body, called by both evictors** — the serial one frees `dead` itself and the
    /// parallel one hands it to the caller's graveyard. They differ in what becomes of the
    /// overflow and in nothing else, and a statement added to one and not the other would
    /// leak a sample's window without failing to compile.
    fn evict_before(&mut self, position: GenomePosition, dead: &mut Vec<SampleLocusObservations>) {
        let first_survivor = first_reaching_summary(&self.held_summaries, position);
        // **Destructured so that the drain and the spare list are two borrows, not one.** Both
        // live on the same `SampleWindow`, and a method call inside the drain would borrow the
        // whole of it a second time.
        // **Every field named, none absorbed by a `..`.** A field added to `SampleWindow` that
        // an eviction has to drain would otherwise compile here and leak.
        let Self {
            held_observations,
            held_summaries,
            held_bodies,
            spare,
            window_coverage,
            source: _,
            spent: _,
            keeps_evidence: _,
            last_drawn: _,
        } = self;
        // The finalised windows go with the records, by their own coordinates.
        window_coverage.evict_before(position);
        held_summaries.drain(..first_survivor);
        held_bodies.drain(..first_survivor.min(held_bodies.len()));
        let room = held_observations.len();
        for record in held_observations.drain(..first_survivor.min(room)) {
            if spare.len() < room {
                spare.push(record);
            } else {
                dead.push(record);
            }
        }
    }
}

/// The cohort's window over one stretch of ground: what every sample recorded there, in the
/// run's sample order, in the two forms the merge reads it in.
///
/// **Two views of one sequence, not two sequences.** Index *i* of a sample's summaries
/// describes index *i* of its observations, so a range of members decided from the summaries
/// addresses the records without a second search.
///
/// **Why the merge is handed both at once.** Closing a cohort locus reads only
/// [`summaries`](Self::summaries) — where each observation sits, and the keep rule's two counts
/// — and assembly reads only [`observations`](Self::observations). Keeping them apart is what
/// lets a run over stored files fill the first from records' heads and leave the second empty
/// until a locus has survived both verdicts (`spec/cohort_merge_psp_path.md` §3.1); direct mode
/// fills both, and pays for the second what it always paid.
#[derive(Debug, Clone, Copy)]
pub struct WindowedCohort<'a> {
    /// One slice per sample: the evidence, in coordinate order — **or `None` when the caller
    /// has summaries and the evidence is still compressed**.
    ///
    /// Absent, nothing before assembly is affected: the closing walk reads summaries and the
    /// two verdicts are passed on those. What changes is who materialises a surviving locus's
    /// members, which becomes the source that kept them
    /// (`spec/cohort_merge_psp_path.md` §3.1).
    pub observations: Option<&'a [&'a [SampleLocusObservations]]>,
    /// One slice per sample at the same indices, or **`None` when the caller has records and
    /// no summaries beside them**.
    ///
    /// Absent, a summary is derived from the record when the walk asks for it — which is
    /// exactly what the closing walk did inline before summaries were held at all, at the
    /// same cost. Present, the walk reads it and never touches the record. The distinction
    /// exists because a run over stored files will have the summaries and *not* the records
    /// (`spec/cohort_merge_psp_path.md` §3.1), and every fixture in this module has the
    /// records and no summaries; both must be expressible or one of them needs a second
    /// merge.
    pub summaries: Option<&'a [&'a [LocusSummary]]>,
    /// One slice per sample at the same index again: the windows that sample's coverage
    /// accumulator has finalised and nothing has evicted, oldest first.
    ///
    /// **Not parallel to the other two** — a window is centred on a position, not on a record,
    /// and a sample has one per covered position while it has one record per locus. It is read
    /// by position, through [`window_coverage_at`](Self::window_coverage_at).
    ///
    /// `None` for a window built from records alone, which is what every fixture that predates
    /// the measurement hands over.
    pub finalised_windows: Option<&'a [&'a [(GenomePosition, WindowCoverage)]]>,
}

impl<'a> WindowedCohort<'a> {
    /// The summary of sample `sample`'s observation at `index`, read or derived.
    #[must_use]
    pub fn summary_at(&self, sample: usize, index: usize) -> Option<LocusSummary> {
        match self.summaries {
            Some(summaries) => summaries[sample].get(index).copied(),
            None => self
                .observations
                .expect("a window with neither summaries nor records describes nothing")[sample]
                .get(index)
                .map(LocusSummary::of),
        }
    }

    /// How many samples the window covers.
    #[must_use]
    pub fn samples(&self) -> usize {
        match (self.summaries, self.observations) {
            (Some(summaries), _) => summaries.len(),
            (None, Some(observations)) => observations.len(),
            (None, None) => 0,
        }
    }

    /// Where sample `sample`'s first held observation begins, if it holds one.
    #[must_use]
    pub fn first_start(&self, sample: usize) -> Option<GenomePosition> {
        self.summary_at(sample, 0).map(LocusSummary::start_position)
    }

    /// Sample `sample`'s window centred on `at`, or `None` where it has none there.
    ///
    /// **`None` is an ordinary answer and not a gap in the measurement.** A sample has a window
    /// at a position only if it has a record there — coverage begins a base later, a repeat
    /// tract emits no loci, no read reached — and the filter that reads this skips a sample it
    /// has no window for, which is what its scorer already does for an absent sample.
    ///
    /// A binary search, because the windows are in ascending coordinate order: a sample holds
    /// one per covered position, so a scan would be linear in the window's ground.
    #[must_use]
    pub fn window_coverage_at(&self, sample: usize, at: GenomePosition) -> Option<WindowCoverage> {
        let per_sample = self.finalised_windows?;
        let windows = per_sample.get(sample)?;
        let found = windows
            .binary_search_by_key(&at, |(centre, _)| *centre)
            .ok()?;
        windows.get(found).map(|(_, window)| *window)
    }
}

/// A borrowed window is a window — so a caller handed one by the cache can pass it on.
impl<'a> From<&WindowedCohort<'a>> for WindowedCohort<'a> {
    fn from(window: &WindowedCohort<'a>) -> Self {
        *window
    }
}

impl WindowedCohort<'static> {
    /// **A window that answers nothing**, and the one definition of it.
    ///
    /// It covers no sample, holds no record and has no window at any position — which is what a
    /// caller assembling a cohort locus from records it already holds has to hand over
    /// ([`CohortObservation::over`](super::build::CohortObservation::over)), and what every
    /// fixture written before the coverage measurement existed means. A locus assembled through
    /// it gets an absent pair for every covering sample.
    ///
    /// **Named rather than `Default`**, because it is not a generally usable view: `samples()`
    /// answers zero for it and `summary_at` panics on it. It is a legal answer to one question —
    /// what was this sample's window here — and the name is what says so.
    #[must_use]
    pub fn with_nothing_measured() -> Self {
        Self {
            observations: None,
            summaries: None,
            finalised_windows: None,
        }
    }
}

/// Records alone are a window: the summaries come from them as they are asked for.
impl<'a> From<&'a [&'a [SampleLocusObservations]]> for WindowedCohort<'a> {
    fn from(observations: &'a [&'a [SampleLocusObservations]]) -> Self {
        Self {
            observations: Some(observations),
            summaries: None,
            finalised_windows: None,
        }
    }
}

impl<S> ObservationCache<S> {
    /// A cache over one source per sample, in the run's sample order, and one reference
    /// accessor for the whole merge.
    ///
    /// Zero samples is not an error here: refusing a zero-sample *run* happens where the run
    /// is configured (spec §7.2), and a cache over an empty cohort covers nothing and hands
    /// nothing out.
    ///
    /// **The accessor is the cache's own and is never shared.** Both callers mint one beside
    /// the padding accessor the run already holds (`walk_reference.accessor()`), because an
    /// accessor shared with another reader would have two callers sliding one window in two
    /// directions — the cost `WindowedRefSeq`'s own documentation records at 14.6 GB of peak
    /// resident memory when it was rebuilt per region instead.
    pub fn over(sources: Vec<S>, reference: Box<dyn MergeReference + Send + Sync>) -> Self {
        Self {
            reference,
            reference_bases: Vec::new(),
            reference_bases_from: None,
            samples: sources
                .into_iter()
                .map(|source| SampleWindow {
                    source,
                    spent: false,
                    spare: Vec::new(),
                    held_observations: Vec::new(),
                    held_summaries: Vec::new(),
                    held_bodies: Vec::new(),
                    keeps_evidence: false,
                    last_drawn: None,
                    window_coverage: WindowCoverageInProgress::new(),
                })
                .collect(),
            covered_to: None,
            keeps_evidence: false,
        }
    }

    /// **The sources back and each sample's finished coverage histogram**, both in the run's
    /// sample order — the same order the sources were handed over in, so entry `i` of either list
    /// is the sample the cohort's `i`th entry describes.
    ///
    /// **Why a cache hands its readers back at all.** A source is not only a reader: a run's
    /// walker carries what its walk saw — the regions it handled, the regions it could not,
    /// and its generators' per-slot counts, which are what explain a covered region emitting
    /// nothing. The cache owns the sources for the whole merge, so without this the merge's
    /// return is where all of that is dropped, and a run report has nothing to say about a
    /// sample beyond its genotypes (arch §3.4).
    ///
    /// **It says nothing about how far the readers got**, and a caller must not infer it: a
    /// merge that failed leaves its sources wherever they stopped, and one that succeeded
    /// leaves them spent. Which of the two happened is the merge's return value, not this.
    ///
    /// **This is the only place a sample's histogram is finished**, and finishing is not a
    /// formality: [`WindowCoverageAccumulator::finish`] closes the centres no cover could — the
    /// last half-window of the sample's own stream, which has no later position to complete it —
    /// and fits the depth axis for a sample whose pass ended before the scale sample was full.
    /// A run that dropped the accumulators instead would produce histograms short of their tails
    /// and, at a sample with fewer than `depth_scale_windows` windows, no histogram at all.
    /// (`spec/window_coverage.md` §3.5.)
    ///
    /// **The tail windows `finish` hands back are dropped here**, deliberately: a window is read
    /// at the locus a builder is building, and by the time this runs every builder has run. They
    /// are the same centres the whole-store recomputation finalises and a run does not
    /// (`cohort_merge::recorded_windows`), and the histogram carries them either way.
    ///
    /// **There is deliberately no sources-only form.** One existed until the histograms did, and
    /// it would now be a shorter name that finishes every accumulator and throws the result away
    /// — and the histograms are the one thing here nothing downstream can recompute without a
    /// second pass over the reads. A caller with no use for them says `.0`, which says so.
    pub fn into_sources_and_histograms(self) -> (Vec<S>, Vec<SampleHistogram>) {
        self.samples
            .into_iter()
            .map(|window| {
                // **Both levels destructured, so a field either type gains has to be answered
                // for here**: this is where the per-sample state stops existing, and one dropped
                // silently is one nothing downstream can ask for. Reaching the accumulator
                // through `.accumulator` would have left a field added beside it invisible.
                let SampleWindow {
                    source,
                    window_coverage,
                    spent: _,
                    spare: _,
                    held_observations: _,
                    held_summaries: _,
                    held_bodies: _,
                    keeps_evidence: _,
                    last_drawn: _,
                } = window;
                let WindowCoverageInProgress {
                    accumulator,
                    // The windows this sample had finalised and nothing had evicted. They are
                    // read at the locus a builder is building, and every builder has run.
                    finalised: _,
                    observed_through: _,
                } = window_coverage;
                let (_tail_windows, histogram) = accumulator.finish();
                (source, histogram)
            })
            .unzip()
    }

    /// **How many samples this cache merges** — the cohort's size, and the length of every
    /// per-sample list a consumer indexes by.
    ///
    /// It exists so that a caller driving the merge can read the cohort's size off the thing it
    /// is merging rather than carry a second copy beside it: the two would be the same number
    /// until the day they were not, and the failure is a genotyper told one cohort size while a
    /// different number of samples is drawn.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    /// [`over`](Self::over) against the fixtures' hand-built reference — four contigs of
    /// `ACGT` repeated, which every fixture in this module and its neighbours sits inside.
    ///
    /// **It exists so that a test says what it is about**: passing a reference at every one of
    /// thirty construction sites would bury the two tests that are about the reference among
    /// twenty-eight that are not, and the alternative — a cache that may hold no reference —
    /// is the one shape that could measure no coverage without saying so.
    #[cfg(test)]
    pub(super) fn over_fixture(sources: Vec<S>) -> Self {
        Self::over(sources, super::fixtures::a_reference())
    }

    /// The reference base at `at`, or `None` where the last cover's ground does not hold it.
    ///
    /// **Canonical `{A,C,G,T,N}`**, as every `RefSeq` fetch is, so a soft-masked FASTA and an
    /// uppercase one give the same answer.
    ///
    /// **`None` means the last cover did not read that position, and there are two ways to meet
    /// it**: a position on a contig the merge has left, and a position no cover has reached yet.
    /// What it does *not* mean is a record without a base — the ground a cover reads covers
    /// every record every sample holds on that contig
    /// ([`fetch_the_ground_this_cover_reached`](Self::fetch_the_ground_this_cover_reached)), so
    /// a consumer walking the held records finds a base at each of their positions.
    ///
    /// It is `pub(super)` and not `pub` because [`cover`](Self::cover) is: outside the merge
    /// nothing can fill the buffer this reads, so a wider accessor could only ever answer
    /// `None`.
    #[must_use]
    // **`not(test)`, because the tests below are still its only callers.** C2 was expected to
    // make it live and did not: the measurement reads the buffer through [`base_in`], which is
    // this accessor's own arithmetic lifted out so the two cannot drift. What the accessor adds
    // is the `&self` form the tests pin the buffer's edges with; it turns live at plan step C4,
    // where a builder reads a base at a locus.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the cover reads the buffer through `base_in`; this `&self` form is what \
                      the module's tests pin the buffer's edges with, and plan step C4's \
                      builder is its first non-test caller"
        )
    )]
    pub(super) fn reference_base_at(&self, at: GenomePosition) -> Option<u8> {
        base_in(&self.reference_bases, self.reference_bases_from?, at)
    }
}

impl<S> ObservationCache<S> {
    /// Every sample's held observations from the first one that reaches into `span`, for the
    /// length of the call — a builder's only way in, and read-only by construction.
    ///
    /// The slices are handed out in the run's sample order, which is what
    /// [`build_region`](super::build::build_region) indexes by.
    ///
    /// **`span`'s two ends do different jobs, and neither is the obvious one.** Its start
    /// selects the window's left edge. Its end does **not** bound the right — trimming there
    /// would cut a locus that opens inside `span` and reaches past it, precisely the deletion
    /// the ownership rule exists to keep whole (spec §6.1) — but neither is it ignored: it is
    /// what this checks against the ground [`cover`](Self::cover) actually reached, so a
    /// builder asking for ground no reader has drawn is refused rather than handed a window
    /// one observation short.
    ///
    /// **Observations that end before `span` are skipped rather than evicted**, so a builder
    /// pays no prefix for ground the organiser has not released yet. Skipping them changes no
    /// locus the caller can own: a locus whose first position is at or after `span.start` has
    /// every member at or after it too, and a locus that chains back past `span.start` does so
    /// through an observation that reaches into `span` — which is kept, so the chain is seen
    /// and the locus still opens before `span`, and is still skipped by the ownership rule.
    pub fn with_observations<R>(
        &self,
        span: GenomeRegion,
        f: impl FnOnce(&WindowedCohort<'_>) -> R,
    ) -> R {
        // `min`/`max` for the reason `cover` and `SampleLocusObservations::reach` use them:
        // `GenomeRegion` has public fields and no constructor ordering them, and these two
        // methods must read an inverted region the same way, or one draws ground the other
        // then trims away.
        let left_edge = GenomePosition {
            contig: span.contig,
            position: span.start.min(span.end),
        };
        let last_base = GenomePosition {
            contig: span.contig,
            position: span.start.max(span.end),
        };

        // **A release check, for the same reason the coordinate-order one is**: a window over
        // ground no cover reached is short, and short is a locus closed over the wrong ground
        // rather than a failure. What it says is that the ground was *drawn* — what is still
        // held is `evict_before`'s caller's business.
        assert!(
            self.covered_to
                .is_some_and(|covered_to| last_base <= covered_to),
            "the window was asked for {span} but cover has only reached {:?} — a cover that \
             failed, or was never made, leaves a window short of what a locus there reaches",
            self.covered_to,
        );

        // One slice reference per sample, one allocation per call. The closer copies them
        // again into its own array (`LocusCloser::over`); doing it twice costs nothing beside
        // the walk, and the alternative — a scratch buffer on the cache — would need interior
        // mutability on a `&self` method that several builders will hold at once.
        let windowing = super::timing::Stopwatch::start();
        // **Both views start at the same index**, because the two arrays are one sequence held
        // twice (`SampleWindow::held_summaries`): the summary at index *i* is the summary of
        // the observation at index *i*, so a member range decided on one addresses the other.
        let mut observations_per_sample: Vec<&[SampleLocusObservations]> =
            Vec::with_capacity(self.samples.len());
        let mut summaries_per_sample: Vec<&[LocusSummary]> = Vec::with_capacity(self.samples.len());
        // **The whole of what each sample holds, not the part from `left_edge`.** A window is
        // centred on a position and a builder asks for it by position, so trimming a prefix
        // would only make the search shorter while giving the caller a slice whose start is not
        // where it thinks; eviction is what bounds this, as it bounds the other two.
        let mut finalised_windows_per_sample: Vec<&[(GenomePosition, WindowCoverage)]> =
            Vec::with_capacity(self.samples.len());
        for sample in &self.samples {
            let from = first_reaching_summary(&sample.held_summaries, left_edge);
            summaries_per_sample.push(&sample.held_summaries[from..]);
            finalised_windows_per_sample.push(&sample.window_coverage.finalised);
            if !self.keeps_evidence {
                observations_per_sample.push(&sample.held_observations[from..]);
            }
        }
        windowing.add_to(&super::timing::WINDOW_NANOS);
        let records = (!self.keeps_evidence).then_some(&observations_per_sample[..]);
        f(&WindowedCohort {
            observations: records,
            summaries: Some(&summaries_per_sample),
            finalised_windows: Some(&finalised_windows_per_sample),
        })
    }


    /// How many observations are held, summed across samples — the size of the window this
    /// cache is the memory of (spec §8).
    ///
    /// **It exists so that "eviction keeps up" can be asserted rather than assumed.** Nothing
    /// in a merge's output shows whether [`evict_before`](Self::evict_before) was ever called:
    /// a driver that never evicted would produce exactly the right answer and hold the whole
    /// stretch while doing it. This is what a test — and later the run's memory report — reads
    /// to tell those two apart.
    pub fn held_observations_len(&self) -> usize {
        self.samples
            .iter()
            .map(|sample| sample.held_observations.len())
            .sum()
    }

    /// Drop everything that ends before `position`. Called once the organiser has released
    /// every locus that could have started there.
    ///
    /// **An observation that began before `position` and reaches into it stays**, which is
    /// what makes eviction safe rather than merely cheap: that observation is what chains a
    /// locus across the evicted point, and dropping it would let a locus that really opens
    /// before `position` appear to open at it — and so be claimed and built by a builder that
    /// never saw its first bases.
    ///
    /// **It is a prefix drain, and on legal input that is a distinction without a
    /// difference.** A sample's records are disjoint and ascending — `build_region` asserts
    /// it — so reach is monotone across one sample's window and the first survivor is the last
    /// non-survivor's successor: nothing behind it could have been kept, nothing after it
    /// could have been dropped. The prefix form is chosen for its cost, which is proportional
    /// to what it drops rather than to the window, and not for a difference in what it keeps.
    pub(super) fn evict_before(&mut self, position: GenomePosition) {
        // The records that will not fit in a sample's spare list. The serial evictor has
        // nowhere to send them, so it frees each sample's before moving to the next — which is
        // where they were freed before this walk and the parallel one shared a body.
        let mut dead = Vec::new();
        for sample in &mut self.samples {
            sample.evict_before(position, &mut dead);
            dead.clear();
        }
    }

    /// [`evict_before`](Self::evict_before) with the samples evicted at the same time as each
    /// other, and what it drops moved out rather than freed here.
    ///
    /// **Eviction is per sample and shares nothing between them**, so this is that walk on a pool:
    /// each sample finds its own first survivor, drains its own prefix into its own spare list,
    /// and pushes what will not fit into a list of its own. Those come back as one buffer, so the
    /// caller still has a single thing to free beside the next round's builders.
    ///
    /// **What the graveyard is for.** What eviction costs is not deciding what to drop — that is
    /// one search per sample — but returning every dropped record's buffers to the allocator.
    /// Moving them out first makes that free a job of its own that the caller can put beside its
    /// builders rather than in front of them; the caller owns the buffer and decides when it dies,
    /// and one that never empties it holds every record the run ever evicted.
    ///
    /// **Why it is worth a method rather than left serial**: it runs between rounds, when every
    /// worker but the organiser's is idle. Measured on 63 tomato accessions over 100 kb it was
    /// 13.6% of the merge and the largest part no thread but one ever ran; putting it here halves
    /// it — 225 ms to 105 over eight merges — and takes 4.4% off the merge. About a third of the
    /// saving comes back as a slower cover, because a record is now freed by whichever worker
    /// drained its sample and the recycled ones cross threads more often on the next draw.
    pub(super) fn evict_before_in_parallel(
        &mut self,
        position: GenomePosition,
        graveyard: &mut Vec<SampleLocusObservations>,
    ) where
        S: Send,
    {
        use rayon::prelude::*;

        let dropped: Vec<Vec<SampleLocusObservations>> = self
            .samples
            .par_iter_mut()
            .map(|sample| {
                let mut dead = Vec::new();
                sample.evict_before(position, &mut dead);
                dead
            })
            .collect();
        for mut sample in dropped {
            graveyard.append(&mut sample);
        }
    }
}

impl<S, E> ObservationCache<S>
where
    S: ObservationSource<Error = E>,
{

    /// Build sample `sample`'s evidence at `index`, which its source kept.
    ///
    /// **Takes `&self`, so several builders may call it at once** — sound because a stored
    /// body is decoded from its own bytes alone, sharing nothing mutable
    /// (`spec/cohort_merge_psp_path.md` §3.2).
    ///
    /// # Errors
    ///
    /// Whatever the sample's source refuses — a damaged block, or a body that will not parse.
    pub(super) fn build_at(
        &self,
        sample: usize,
        index: usize,
    ) -> Result<SampleLocusObservations, E> {
        let window = &self.samples[sample];
        window.source.build(window.held_bodies[index].clone())
    }
    /// Draw every sample forward until `region` is covered, and far enough past it to hold two
    /// things: what a locus starting inside it can reach (spec §6.4), and half a window more,
    /// which is what the coverage measurement's last centres need (`window_coverage.md` §3.3).
    ///
    /// **The first of those is the whole difficulty, and it is a fixpoint across samples.** A
    /// locus opening inside `region` closes only when the next observation begins beyond its
    /// reach, and each observation the reach pulls in may push it further — so the chain's
    /// reach starts half a window past `region`'s last base
    /// ([`where_the_chain_reaches_before_any_sweep`](Self::where_the_chain_reaches_before_any_sweep))
    /// and grows with every observation that begins at or
    /// before it, which is the same `<=` the closer chains on (spec §4.1). One sample's
    /// deletion is what makes another sample's later observation part of the locus, so the
    /// samples are swept repeatedly until a whole sweep moves the reach no further. A single
    /// sweep stops short whenever the widening sample is swept after the sample it widens
    /// onto, and three samples can need three sweeps
    /// (`a_chain_that_needs_a_third_sweep_is_drawn_whole`): the sweep count is a property of
    /// the data, not of the code.
    ///
    /// **How far past the look-ahead this can go is not bounded here.** The reach grows only
    /// through observations that chain into it, so what limits it is the widest observation the
    /// generator can mint — the reach ceiling, `max_record_span` (spec §1.3) — times the
    /// length of the chain. On ground where observations overlap wall to wall (spec §7.1) one
    /// cover can draw a whole segment. **Nothing here reads or checks that ceiling**, and that
    /// is unchanged now the psp header carries it (2026-09-04, `run_streaming.md` plan step
    /// E4): the ceiling is a sizing fact and this cache has no capacity to size — it grows on
    /// demand. What landed is that the number is now *available* to a run over stored files,
    /// which takes the maximum over its cohort's files at open
    /// ([`PspVariantCaller::observation_reach_ceiling`](crate::ng::run::PspVariantCaller)),
    /// where before it existed only as the configured `max_record_span` of a walker that a
    /// calling run over psps does not have. Spec §13's deferral is discharged; no refusal
    /// accompanies it, as §13 says none is needed.
    ///
    /// **What a cover costs.** Each sweep visits every sample and re-reads its held window, so
    /// a cover is `sweeps × (samples + held)`. Two sweeps is the ordinary count; the worst
    /// case is `samples` sweeps, when a chain of overlapping observations runs through the
    /// cohort in decreasing sample order. Measured by the D1 review on a synthetic cohort in a
    /// release build, over 20-base regions: one extra sweep costs 3.1 µs at 3,000 samples
    /// against 2.87 ms for the whole cover, and that worst case costs 28 ms for one 11-base
    /// region at 3,000 samples. The `held` term stays short **only while the organiser evicts
    /// at the pace it releases ground**, which this module cannot enforce — `evict_before` is
    /// the organiser's call (milestone E). At 1,000 samples the same walk costs 616 µs a cover
    /// with 4 observations held per sample and 1,028 µs with 200. **The look-ahead adds to the
    /// `held` term rather than to the sweep count**: half a window more of held summaries per
    /// sample, which `window_coverage.md` §3.3 prices at about 250 at three reads a position and
    /// some tens of kilobytes across a thousand samples.
    ///
    /// **The window overshoots by at most one observation per sample**, and that is what a
    /// forward reader costs: the only way to know whether the next observation begins beyond
    /// the chain's reach is to draw it, and once drawn it is held rather than thrown away.
    ///
    /// Failures from a source end the cover and are passed through unchanged. The window keeps
    /// whatever was drawn before the failure and the cover can be made again — **which asks of
    /// a source that it may be polled after yielding `Err`**, something `Iterator` does not
    /// grant on its own. A source that cannot honour it must yield `None` after failing, which
    /// this reads as a spent sample rather than as a retry.
    pub(super) fn cover(&mut self, region: GenomeRegion) -> Result<(), E>
    where
        E: From<ReferenceUnreadable>,
    {
        let mut chain_reach = self.where_the_chain_reaches_before_any_sweep(region);

        // The fixpoint: sweep until a whole sweep moves nothing.
        while self.sweep(&mut chain_reach)? {}

        self.close_the_cover(region, chain_reach)
    }

    /// How far the chain reaches before a single sweep has run: **half a window past the
    /// region's last base**, stopped at the contig's end. The sweeps then widen it.
    ///
    /// **The half window is the accumulator's look-ahead, and getting it wrong is silent** (spec
    /// `window_coverage.md` §3.3, §6 trap 3). A window centred at `p` is finalised once that
    /// sample's stream **reaches** `p + WINDOW_BP / 2` — the accumulator closes a centre at
    /// `centre + half <= position`, not past it, which is what makes half a window exactly
    /// enough rather than one base short. A cover that stopped at its region's last base would
    /// leave the region's last centres unfinalised: the builder handed that region would find no
    /// window there and every sample would read as absent at exactly the loci nearest each region
    /// boundary. No test fails, no VCF byte moves, and the only thing that can see it is the
    /// whole-store recomputation (`cohort_merge::recorded_windows`, and the probe that reads what
    /// it writes).
    ///
    /// **Drawing that far is not the same as finalising**, and the difference is what makes this
    /// enough rather than merely necessary. What closes a centre is a *position* at or beyond
    /// `p + WINDOW_BP / 2`, not a reach — but the draw stops at the first record starting past
    /// the chain ([`draw_to`](SampleWindow::draw_to)) and **holds** it, and every held record no
    /// accumulator has seen is observed by this same cover
    /// ([`read_the_ground_and_measure_coverage_over_it`](Self::read_the_ground_and_measure_coverage_over_it)).
    /// So a sample with any record at all past the chain has had it observed here. What is left
    /// over is a sample with no such record — the last half-window of its last records on its last
    /// contig — whose centres a contig change closes, or, at the very end of a sample's stream,
    /// nothing in the merge does: `WindowCoverageAccumulator::finish` has no caller until plan
    /// step C5.
    ///
    /// **The contig's end is a hard stop**, because the ground a cover reads runs to the chain
    /// ([`ground_this_cover_holds_on`](Self::ground_this_cover_holds_on)) and a fetch past a
    /// contig's last base is refused rather than truncated
    /// ([`RefSeqError::OutOfBounds`](crate::ng::ref_seq::RefSeqError)). It costs nothing where it
    /// does not bind: nothing sits past a contig's end for a chain to reach.
    ///
    /// `max` on the region's own ends for the reason `SampleLocusObservations::reach` uses it:
    /// `GenomeRegion` has public fields and no constructor enforcing `start <= end`, and an
    /// inverted region must not put the chain's reach before the ground it is meant to cover.
    fn where_the_chain_reaches_before_any_sweep(&self, region: GenomeRegion) -> GenomePosition {
        let last_base = region.end.max(region.start).get();
        // **A contig the table does not name gets no clamp**, and needs none: the fetch over its
        // ground refuses it by name (`RefSeqError::UnknownContig`) whatever the reach was, which
        // is what `a_failed_fetch_names_the_ground_the_cover_held` covers.
        let contig_length = self
            .reference
            .contigs()
            .entries
            .get(region.contig.get() as usize)
            .map_or(u64::MAX, |contig| contig.length);
        // **The clamp may not pull the reach *back*.** A region already at or past its contig's
        // last base — which a caller may hand over, since nothing here checks a region against the
        // table — would otherwise be given a chain behind the ground its own builder was handed,
        // and `with_observations` would then refuse that builder's window.
        let no_further_than = contig_length.max(last_base);
        let with_the_look_ahead =
            last_base.saturating_add(u64::from(window_coverage::WINDOW_BP / 2));
        GenomePosition {
            contig: region.contig,
            position: Position(with_the_look_ahead.min(no_further_than)),
        }
    }

    /// [`cover`](Self::cover), with each sweep's samples swept concurrently.
    ///
    /// **Same fixpoint, a different schedule.** The serial sweep is Gauss-Seidel — sample `j`
    /// sees the reach sample `i < j` widened inside the same sweep — and this is Jacobi: every
    /// sample is drawn against the reach the last sweep ended on, and the sweep's answer is the
    /// widest of theirs. Drawing is monotone in the reach it is given and the reach only ever
    /// grows, so both schedules climb to the same least fixpoint above the region's last base,
    /// and the held window is the same because every sample's last draw is against that
    /// fixpoint. What differs is the sweep count: a chain that runs through the cohort in
    /// decreasing sample order costs one sweep per link either way, but a chain the serial form
    /// follows within one sweep costs this one a sweep per link.
    pub(super) fn cover_in_parallel(&mut self, region: GenomeRegion) -> Result<(), E>
    where
        S: Send,
        E: Send + From<ReferenceUnreadable>,
    {
        use rayon::prelude::*;

        let mut chain_reach = self.where_the_chain_reaches_before_any_sweep(region);
        loop {
            let snapshot = chain_reach;
            let widest = self
                .samples
                .par_iter_mut()
                // The return type is spelled out because `E` now carries a `From` bound, and
                // without it inference reads the `?` below as converting into that bound's
                // type rather than into the source's own error.
                .map(|sample| -> Result<GenomePosition, E> {
                    let drawing = super::timing::Stopwatch::start();
                    let mut reach = snapshot;
                    sample.draw_to(&mut reach)?;
                    drawing.add_to(&super::timing::COVER_BUSY_NANOS);
                    Ok(reach)
                })
                .try_reduce(|| snapshot, |left, right| Ok(left.max(right)))?;
            super::timing::COVER_SWEEPS.add(1);
            if widest == chain_reach {
                break;
            }
            chain_reach = widest;
        }

        self.close_the_cover(region, chain_reach)
    }

    /// What both fixpoints do once the drawing has stopped: read the ground, then record how
    /// far the cover got.
    ///
    /// **One copy, because the two covers must not drift.** They differ only in how they run
    /// the fixpoint, and the standing oracle for this module is that they produce the same
    /// answer; a step that added a statement to one tail and not the other would break that
    /// without failing to compile.
    ///
    /// **The fetch happens before either fact is recorded**, so a cover that could not read its
    /// ground is not counted as covered: `covered_to` is what
    /// [`with_observations`](Self::with_observations) refuses ground against, and advancing it
    /// over ground whose bases were never read would hand a builder a window this cache cannot
    /// describe.
    ///
    /// # Errors
    ///
    /// Whatever the reference fetch refuses, converted by the caller.
    fn close_the_cover(
        &mut self,
        region: GenomeRegion,
        chain_reach: GenomePosition,
    ) -> Result<(), E>
    where
        E: From<ReferenceUnreadable>,
    {
        self.read_the_ground_and_measure_coverage_over_it(region, chain_reach)?;
        self.keeps_evidence |= self.samples.iter().any(|sample| sample.keeps_evidence);
        self.covered_to = Some(
            self.covered_to
                .map_or(chain_reach, |reached| reached.max(chain_reach)),
        );
        Ok(())
    }

    /// Read the reference over the ground this cover holds and feed every record it added to
    /// its sample's accumulator — **one contig at a time**, in ascending order, ending on the
    /// region's own.
    ///
    /// **A cover crosses contigs, and a buffer is one contig's.** Drawing stops at the first
    /// record past the reach, and a reach on a later contig is past every position of an
    /// earlier one — so the cover that first reaches contig *n* also draws whatever the sample
    /// still had on contig *n − 1*. Those records are real covered positions of that sample and
    /// their bases are on the contig they sit on, so this reads each contig the cover added
    /// records on before observing the records there.
    ///
    /// **The buffer is left holding the region's own ground, and that takes a second fetch on the
    /// covers that crossed.** The contigs are observed *ascending*, because that is the order the
    /// accumulators demand; the region's own is last only while no contig sorts after it, and one
    /// regularly does — a sample is drawn one record past the reach, that record can be the first
    /// of the next contig, and the next contig sorts after. So a cover that ended on another
    /// contig reads the region's ground again before returning, observing nothing, because every
    /// record on it was observed on the first visit. Without it the buffer a builder reads through
    /// [`reference_base_at`](Self::reference_base_at) is the *next* contig's, and that accessor
    /// answers `None` at every position of the region the builder owns.
    ///
    /// **Each record exactly once**, which the per-sample cursor is for: a record is held across
    /// every cover it reaches into, and each of those covers reads ground containing it.
    ///
    /// **A record spanning more than one base has its evidence built here in psp mode**, which
    /// is the one place this measurement costs a run over stored files anything. Plan step B2
    /// measured how often that is: 1 record in 871 on the six-accession slice and 1 in 1,221 on the
    /// wider store.
    ///
    /// # Errors
    ///
    /// A failed reference fetch, named by the ground it was over — and whatever building a kept
    /// record's evidence refuses, unchanged.
    fn read_the_ground_and_measure_coverage_over_it(
        &mut self,
        region: GenomeRegion,
        chain_reach: GenomePosition,
    ) -> Result<(), E>
    where
        E: From<ReferenceUnreadable>,
    {
        let contigs = self.contigs_this_cover_must_read(region);
        // Which contig the buffer will be left holding, so the refetch below can say whether it
        // is the region's.
        let last_contig_read = contigs.last().copied();
        for contig in contigs {
            let ground = self.ground_this_cover_holds_on(contig, region, chain_reach);
            self.fetch_the_ground(ground)?;
            // **Destructured so the buffer and the samples are two borrows.** The bases belong
            // to the cache and the accumulators to the samples, and a method call would borrow
            // the whole of `self` for both.
            let Self {
                samples,
                reference_bases,
                reference_bases_from,
                reference: _,
                keeps_evidence: _,
                covered_to: _,
            } = self;
            // **The buffer's own first position, not the ground's**, so that which contig the
            // bases are on and where they start cannot be two answers: `fetch_the_ground` is
            // what set it, and it is `None` only when the fetch failed and `?` already left.
            let bases_from = reference_bases_from
                .expect("a fetch that returned `Ok` recorded where its bases start");
            for sample in &mut *samples {
                sample.measure_coverage_on(reference_bases, bases_from)?;
            }
        }
        // **The region's ground back in the buffer**, where observing ascending did not leave it
        // there. Nothing is observed against it: every record on this contig was seen on its own
        // turn above, and each sample's cursor is past them.
        if last_contig_read != Some(region.contig) {
            let ground = self.ground_this_cover_holds_on(region.contig, region, chain_reach);
            self.fetch_the_ground(ground)?;
        }
        Ok(())
    }

    /// The contigs this cover must read the reference over, ascending: every contig a record no
    /// accumulator has seen sits on, plus the region's.
    ///
    /// **Three at most in practice** — the one the cover is leaving, the region's own, and the one
    /// the overshoot draw reached ahead of it — so this is a short sorted list rather than a set.
    fn contigs_this_cover_must_read(&self, region: GenomeRegion) -> Vec<ContigId> {
        let mut contigs = vec![region.contig];
        for sample in &self.samples {
            let unobserved = sample
                .window_coverage
                .summaries_not_yet_observed(&sample.held_summaries);
            for summary in unobserved {
                let contig = summary.region.contig;
                if !contigs.contains(&contig) {
                    contigs.push(contig);
                }
            }
        }
        contigs.sort_unstable();
        contigs
    }

    /// The reference bases over one stretch of one contig, into the buffer every sample reads
    /// by offset.
    ///
    /// **One fetch a contig a cover, whatever the cohort's size.** The bases are the same for
    /// every sample; only which of them a sample's records fall on differs.
    ///
    /// # Errors
    ///
    /// A failed fetch, named by the ground it was over. It ends the cover rather than being
    /// absorbed: one failure costs every sample its coverage over that stretch, and a run that
    /// carried on would report a window whose denominator quietly excluded it. **A ground whose
    /// first base is 0 fails here too** — there is no position 0 on a contig, and the fetch
    /// refuses rather than translating a caller's bug into a base.
    fn fetch_the_ground(&mut self, ground: GenomeRegion) -> Result<(), E>
    where
        E: From<ReferenceUnreadable>,
    {
        // **The old ground stops being answerable the moment a new one is asked for.** A failed
        // fetch leaves the buffer as it was (`RefSeq::fetch_into`'s contract), so without this
        // a caller that retried a cover would read the previous cover's bases as this one's.
        self.reference_bases_from = None;
        self.reference
            .fetch_into(
                ground.contig,
                ground.start.get(),
                ground.len(),
                &mut self.reference_bases,
            )
            .map_err(|source| {
                E::from(ReferenceUnreadable {
                    region: ground,
                    source,
                })
            })?;
        self.reference_bases_from = Some(GenomePosition {
            contig: ground.contig,
            position: ground.start,
        });
        // **Release what the merge has walked past**, exactly as the run's padding accessor
        // does at the record (`run::records`). A windowed reader extends its buffer forward
        // and only ever shrinks here, so a merge that never released would end a contig
        // holding every base it had covered — about 250 MB on human chromosome 1, against a
        // walk that otherwise peaks near 25 MB (`ref_seq.rs`'s
        // `a_forward_walk_holds_the_whole_span_unless_it_releases`).
        self.reference.evict_before(ground.start.get());
        Ok(())
    }

    /// The stretch of `contig` this cover's records lie on: every record any sample still holds
    /// there, widened to the region when `contig` is the region's own.
    ///
    /// **The ground is what the cache holds, not what the region asked for**, and the two are
    /// different at both ends. Past the region: an observation chaining beyond it widens the
    /// reach, and every sample is drawn one observation *past* that reach — the only way to
    /// know a record is beyond the reach is to draw it, and once drawn it is held rather than
    /// thrown away. Before the region: a record drawn by an earlier cover and still held is
    /// handed to the builders again, and this buffer is refilled rather than appended to, so
    /// the bases the cover that drew it fetched are gone.
    fn ground_this_cover_holds_on(
        &self,
        contig: ContigId,
        region: GenomeRegion,
        chain_reach: GenomePosition,
    ) -> GenomeRegion {
        // `min`/`max` for the reason `SampleLocusObservations::reach` uses them: `GenomeRegion`
        // has public fields and no constructor enforcing `start <= end`.
        let on_the_regions_own_contig = contig == region.contig;
        let mut first_base = None;
        let mut last_base = None;
        if on_the_regions_own_contig {
            first_base = Some(region.start.min(region.end));
            last_base = Some(chain_reach.position.max(region.start.min(region.end)));
        }
        for sample in &self.samples {
            // **The first base is an end, not a scan.** A sample's held summaries ascend in
            // `(contig, position)` (`draw_next`'s ordering check), so the ones on `contig` are a
            // contiguous run and the earliest of them is the first. Scanning instead would cost
            // this `samples × held` per contig on every cover.
            let held = summaries_on(&sample.held_summaries, contig);
            if let Some(first) = held.first() {
                let start = first.start_position().position;
                first_base = Some(first_base.map_or(start, |first: Position| first.min(start)));
            }
            // **The last base is an end only on the region's own contig**, and the difference is
            // what a source out of order can turn into a wrong answer. There, the chain already
            // dominates every held record: each one starting at or before it widened the chain to
            // its own reach, and at most one record begins past it, which is the last. On a contig
            // the cover is *leaving* there is no chain reach to lean on, and `held.last()` is the
            // furthest only while reaches are monotone — which nothing here checks, since
            // `draw_next` compares starts. A record reaching past its successor would then have
            // ground the fetch never read, and the measurement would panic on a missing base.
            let furthest = if on_the_regions_own_contig {
                held.last().map(|summary| summary.reach())
            } else {
                held.iter().map(|summary| summary.reach()).max()
            };
            if let Some(reach) = furthest {
                last_base = Some(last_base.map_or(reach, |last: Position| last.max(reach)));
            }
        }
        // A contig reaches this function only because the region names it or a record sits on
        // it, so both ends are set; the fallback is the region's own ground and cannot fire.
        GenomeRegion {
            contig,
            start: first_base.unwrap_or(region.start.min(region.end)),
            end: last_base.unwrap_or(region.end.max(region.start)),
        }
    }

    /// One sweep of every sample against `chain_reach`. Answers whether any of them moved it.
    fn sweep(&mut self, chain_reach: &mut GenomePosition) -> Result<bool, E> {
        let mut reach_grew = false;
        super::timing::COVER_SWEEPS.add(1);
        for sample in &mut self.samples {
            // Not `reach_grew |= sample.draw_to(…)?`: that reads as though the call could be
            // skipped, and a later `||` in its place would skip it.
            let drawing = super::timing::Stopwatch::start();
            let grew = sample.draw_to(chain_reach)?;
            drawing.add_to(&super::timing::COVER_BUSY_NANOS);
            if grew {
                reach_grew = true;
            }
        }
        Ok(reach_grew)
    }
}

impl<S, E> SampleWindow<S>
where
    S: ObservationSource<Error = E>,
{
    /// Widen `chain_reach` with this sample's observations while they begin at or before it,
    /// drawing from the source as the held ones run out. Answers whether the reach moved.
    ///
    /// **It starts from the window's first observation every time, rather than remembering
    /// where the last sweep stopped.** Re-reading one that is already inside the reach cannot
    /// move it — the reach only ever grows — so the scan is idempotent, and keeping no mark
    /// means there is no mark for eviction to correct or for a later region to inherit: a
    /// cover is a function of the window, the source and the region it was asked for. A mark
    /// carried across covers is not merely stale but can point past the end of a window that
    /// lost two entries at once, which `a_survivor_of_an_eviction_still_widens_the_next_reach`
    /// and `two_evicted_at_once_leave_the_window_sound` exist to catch.
    fn draw_to(&mut self, chain_reach: &mut GenomePosition) -> Result<bool, E> {
        let mut reach_grew = false;
        let mut considered = 0;
        loop {
            if considered == self.held_summaries.len() {
                // Nothing left in this sample: it cannot widen the reach again.
                let Some(observation) = self.draw_next()? else {
                    break;
                };
                match observation {
                    Drawn::Built(record) => {
                        self.held_summaries.push(LocusSummary::of(&record));
                        self.held_observations.push(record);
                    }
                    Drawn::Kept { summary, body } => {
                        self.keeps_evidence = true;
                        self.held_summaries.push(summary);
                        self.held_bodies.push(body);
                    }
                }
            }
            let observation = self.held_summaries[considered];
            if observation.start_position() > *chain_reach {
                // Held, but beyond the reach of any locus this cover can see. It stays in the
                // window, and the next cover reconsiders it against a later reach.
                break;
            }
            considered += 1;
            let reach = observation.reach_position();
            if reach > *chain_reach {
                *chain_reach = reach;
                reach_grew = true;
            }
        }
        Ok(reach_grew)
    }

    /// Feed this sample's accumulator the records it holds on `bases_from`'s contig that the
    /// accumulator has not seen, against the bases just read for that contig.
    ///
    /// `bases` is that contig's ground and `bases_from` where its first byte sits — **one
    /// argument for both, so which contig the bases are on cannot be two answers.** Every
    /// record this walks lies inside them, because the ground is computed from these very
    /// records ([`ObservationCache::ground_this_cover_holds_on`]).
    ///
    /// **Records on a later contig are left for that contig's turn**, and the cursor stops with
    /// them: the caller walks the contigs in ascending order, so a record skipped here is
    /// reached by a later call of the same cover.
    ///
    /// # Errors
    ///
    /// Whatever building a kept record's evidence refuses, unchanged.
    fn measure_coverage_on(&mut self, bases: &[u8], bases_from: GenomePosition) -> Result<(), E> {
        // Destructured so that the source can be asked to build while the measurement beside
        // it is borrowed mutably: they are separate fields, and only a destructuring says so.
        let Self {
            source,
            held_summaries,
            held_observations,
            held_bodies,
            window_coverage,
            spent: _,
            spare: _,
            keeps_evidence: _,
            last_drawn: _,
        } = self;
        let first_new = window_coverage.first_not_yet_observed(held_summaries);
        for index in first_new..held_summaries.len() {
            let summary = held_summaries[index];
            if summary.region.contig != bases_from.contig {
                // Ascending, so this and everything after it belongs to a later contig's turn.
                break;
            }
            let evidence = match held_observations.get(index) {
                Some(record) => EvidenceForOneRecord::InHand(record),
                // The psp path holds bytes where direct mode holds records, at the same index.
                None => EvidenceForOneRecord::Kept(held_bodies[index].clone()),
            };
            window_coverage.observe(
                summary,
                evidence,
                |body| source.build(body),
                bases,
                bases_from,
            )?;
        }
        Ok(())
    }

    /// The next observation from the source, or `None` once it is spent.
    fn draw_next(&mut self) -> Result<Option<Drawn>, E> {
        if self.spent {
            return Ok(None);
        }
        // **The offer, and it is the whole of the recycling.** A source that mints records
        // fills this one instead of allocating; a source that cannot drops it and allocates,
        // which is what every plain iterator does.
        let spare = self.spare.pop();
        let Some(next) = self.source.next_drawn(spare).transpose()? else {
            self.spent = true;
            return Ok(None);
        };
        // **Every record the merge will later free passes here**, which is why the count that
        // prices Milestone G is taken here rather than at the mint: the generator's own
        // records include ones the region clamp discards, and those the merge never owns.
        super::timing::RECORDS_DRAWN.add(1);
        if let Drawn::Built(record) = &next {
            super::timing::OBSERVATIONS_DRAWN.add(record.observations.len() as u64);
        }
        let summary = next.summary();

        // **A release check, not a `debug_assert!`** — the release profile is the one this
        // repo runs, and a source that goes backwards is silent otherwise: the draw loop stops
        // at the first observation beyond the reach, so an out-of-order one would end a cover
        // early and hand a builder a locus cut short, which is a wrong answer rather than a
        // failure. The closer's own ordering check (`LocusCloser::over`) compares only within
        // one locus and would not see it. One comparison per observation drawn.
        //
        // **When observations are decoded from a psp file this must become a `RunError`**
        // beside `ObservationExceedsReachCeiling` (arch §5), like the producer-guarantee
        // checks in `build.rs` — a source out of coordinate order is then a fact about the
        // file rather than a bug in this crate, and this is the first such check the psp path
        // reaches.
        // **Strictly after, not at or after, and the measurement is what tightened it.** The
        // per-sample cursor that makes a record observed exactly once is a *start*, compared with
        // `<=` ([`WindowCoverageInProgress::first_not_yet_observed`]) — so two records sharing a
        // start are one record to it. Inside a single cover both are still observed, because the
        // observation walks by index; split across two covers the second falls out of the
        // not-yet-observed suffix and **never reaches the accumulator at all**, leaving that
        // sample short of one record's positions with nothing to say so. Refusing the shape here
        // is what keeps the cursor's representation adequate to what the window can hold.
        //
        // **It is still not the disjointness the window is really read by.** Two records that
        // start apart and overlap pass this, and `build_region`'s own producer-guarantee
        // assertion is what catches them — deliberately, since
        // `parallel.rs`'s `a_builder_that_panics_leaves_the_cache_in_the_callers_hands_and_advanced`
        // passes exactly such a pair through this cache to reach it.
        let start = summary.start_position();
        if let Some(previous) = self.last_drawn {
            assert!(
                start > previous,
                "this sample's source is not in coordinate order: {} does not begin after \
                 {previous:?}",
                summary.region,
            );
        }
        self.last_drawn = Some(start);

        Ok(Some(next))
    }
}

/// One analysed region divided into the regions single builders own — `building_region_width`
/// bases each, adjacent, in genome order, the last one clamped to the analysed region's own
/// last base.
///
/// **It lives here, beside the cache, because it is the organiser's geometry**: milestone E
/// hands these regions out to builders, and a second derivation of the clamp there is exactly
/// the defect this function's tests exist to catch. Its first caller is the serial driver that
/// reads through the cache (`super::serial::merge_cohort_through_cache`).
///
/// **It is a function of its own so that the division can be checked**, which no merge's output
/// can do: dividing the ground changes nothing a caller sees — that is the whole claim of the
/// cached driver — so a driver that ignored the width and built each analysed region as one
/// would give the right answer and lose everything the division buys.
///
/// **The clamp is load-bearing.** A locus belongs to the builder whose region holds its first
/// position, so a last building region running past the analysed ground would claim loci the
/// run never analysed.
///
/// An analysed region whose ends are the wrong way round is read as the ground it names, the
/// same defence `SampleLocusObservations::reach` makes — belt-and-braces, since the drivers
/// refuse such a region before they get here
/// (`super::serial::refuse_malformed_analysed_regions`).
pub fn building_regions_of(
    analysed_region: GenomeRegion,
    building_region_width: CohortLocusBuilderRegionsLen,
) -> impl Iterator<Item = GenomeRegion> {
    let first_base = analysed_region.start.min(analysed_region.end);
    let last_base = analysed_region.end.max(analysed_region.start);
    let bases_per_region = u64::from(building_region_width.get());

    let building_region_from = move |start: Position| GenomeRegion {
        contig: analysed_region.contig,
        start,
        // Saturating, then clamped: the arithmetic must not wrap at the top of the coordinate
        // space, and the region must not reach past the ground the run analysed. The `- 1`
        // cannot underflow — the width wraps a `NonZeroU32`.
        end: Position(
            start
                .0
                .saturating_add(bases_per_region - 1)
                .min(last_base.0),
        ),
    };

    std::iter::successors(Some(building_region_from(first_base)), move |previous| {
        // `checked_add`, because a region ending on the last base of the coordinate space has
        // no successor and the addition would wrap to zero — and with the saturating arithmetic
        // above, a wrap would divide the genome for ever rather than panic.
        let next_first_base = previous.end.0.checked_add(1)?;
        (next_first_base <= last_base.0).then(|| building_region_from(Position(next_first_base)))
    })
}

/// The first held observation that reaches `position` or beyond — the window's left edge for
/// a call at `position`, and the window's length when every one of them ends before it.
/// [`first_reaching_index`] over summaries, which every window has and only some have records
/// for.
/// The base at `at` in a buffer whose first byte is `bases_from`, or `None` where the buffer
/// does not hold that position.
///
/// **One derivation of the offset, read by both the cache's accessor and the measurement.**
/// The two would agree until the day one of them was changed, and the failure is a window
/// filled from bases one position over.
fn base_in(bases: &[u8], bases_from: GenomePosition, at: GenomePosition) -> Option<u8> {
    if bases_from.contig != at.contig || at.position < bases_from.position {
        return None;
    }
    // The guard above is what makes this subtraction total.
    let offset = at.position.get() - bases_from.position.get();
    bases.get(usize::try_from(offset).ok()?).copied()
}

/// The stretch of `held` that sits on `contig` — a contiguous run, because a sample's summaries
/// ascend in `(contig, position)` (`SampleWindow::draw_next`'s ordering check).
fn summaries_on(held: &[LocusSummary], contig: ContigId) -> &[LocusSummary] {
    let from = held.partition_point(|summary| summary.region.contig < contig);
    let to = held.partition_point(|summary| summary.region.contig <= contig);
    &held[from..to]
}

fn first_reaching_summary(held: &[LocusSummary], position: GenomePosition) -> usize {
    held.partition_point(|summary| summary.reach_position() < position)
}


#[cfg(test)]
mod tests {
    use super::super::build::build_region_windowed;
    use super::super::fixtures::{
        ReferenceCountingItsReleases, ReferenceFailureRecorded, SourceFailed, position_on, region,
        region_on,
    };
    use super::*;
    use crate::ng::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::ng::types::{ContigId, Position, ReadGroupId};
    use std::cell::Cell;
    use std::rc::Rc;
    use std::sync::Arc;

    /// One sample's record over `region`. What it showed is irrelevant to the cache — which
    /// reads nothing but where an observation begins and how far it reaches — so every fixture
    /// observation carries the same single sighting. It is a well-formed one, because the
    /// differential test at the end of this file hands these to a real builder.
    fn observation_over(region: GenomeRegion) -> SampleLocusObservations {
        let width = usize::try_from(region.end.0.saturating_sub(region.start.0) + 1)
            .expect("a fixture region fits in memory");
        SampleLocusObservations {
            region,
            reference_bases: vec![b'A'; width].into_boxed_slice(),
            observations: vec![SequenceObservation {
                bases: Box::from(&b"C"[..]),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 3,
                num_fwd: 3,
                q_sum: crate::ng::types::SummedLogError::from_nats(-6.0),
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

    /// One sample's reader, counting what was drawn out of it. The count is what says whether
    /// the cache stopped drawing where it claims to.
    struct CountingSource {
        remaining: std::vec::IntoIter<Result<SampleLocusObservations, SourceFailed>>,
        drawn: Rc<Cell<usize>>,
    }

    impl Iterator for CountingSource {
        type Item = Result<SampleLocusObservations, SourceFailed>;

        fn next(&mut self) -> Option<Self::Item> {
            let next = self.remaining.next();
            if next.is_some() {
                self.drawn.set(self.drawn.get() + 1);
            }
            next
        }
    }

    fn width(bases: u32) -> CohortLocusBuilderRegionsLen {
        CohortLocusBuilderRegionsLen(
            std::num::NonZeroU32::new(bases).expect("a fixture width is non-zero"),
        )
    }

    /// One sample's reader over records already in hand, with the failure type these tests
    /// name spelled once rather than at every fixture.
    fn a_source(
        records: Vec<Result<SampleLocusObservations, SourceFailed>>,
    ) -> std::vec::IntoIter<Result<SampleLocusObservations, SourceFailed>> {
        records.into_iter()
    }

    /// The base the fixture reference holds at a 1-based position: it repeats `ACGT`, so
    /// position 1 is `A`, position 5 is `A` again, and a test can say what it expects without
    /// counting.
    fn fixture_base_at(position: u64) -> u8 {
        b"ACGT"[usize::try_from((position - 1) % 4).expect("a remainder of four fits")]
    }

    /// A cover reads the reference over the ground it drew, and every sample reads it from the
    /// same buffer at its own offset. This is the first half: the bases are the reference's,
    /// at the positions they belong to.
    #[test]
    fn a_cover_holds_the_reference_bases_over_the_ground_it_drew() {
        let source = a_source(vec![Ok(observation_over(region(40, 40)))]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source and reference hold");

        for position in 40..=50 {
            assert_eq!(
                cache.reference_base_at(position_on(0, position)),
                Some(fixture_base_at(position)),
                "the base at {position} is not the reference's",
            );
        }
    }

    /// **The ground is the cover's reach and not the region's end**, which is the whole reason
    /// the fetch happens after the fixpoint rather than before it: an observation chaining past
    /// the region widens what was drawn, and every position of it needs a base.
    ///
    /// The observation reaches well past the look-ahead's own half window
    /// ([`where_the_chain_reaches_before_any_sweep`](ObservationCache::where_the_chain_reaches_before_any_sweep)), so that what the
    /// ground follows here is the chain and not the constant.
    #[test]
    fn the_ground_read_reaches_as_far_as_the_cover_drew_and_no_further() {
        // One observation opens inside the region and reaches far past its end — past the
        // look-ahead too, so widening it is what decides the ground.
        let source = a_source(vec![Ok(observation_over(region(48, 400)))]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source and reference hold");

        assert_eq!(
            cache.reference_base_at(position_on(0, 400)),
            Some(fixture_base_at(400)),
            "the last base the cover drew to has no reference base",
        );
        assert_eq!(
            cache.reference_base_at(position_on(0, 401)),
            None,
            "a base past the cover's reach was fetched, so the ground is not the one drawn",
        );
        assert_eq!(
            cache.reference_base_at(position_on(0, 39)),
            None,
            "a base before the region's own start was fetched",
        );
    }

    /// **A cover reads half a window past its region even where nothing chains there**, which is
    /// the look-ahead (`spec/window_coverage.md` §3.3): a window centred on the region's last
    /// base is finalised only once that sample's stream has passed it by half a window, and a
    /// cover stopping at the region's end would leave the builder nothing there.
    ///
    /// **Its failure is silent** — no panic, no moved VCF byte — so this is a test rather than an
    /// oracle's business, and the number comes from the constant rather than being retyped.
    #[test]
    fn a_cover_reads_half_a_window_past_its_region_with_nothing_chaining_there() {
        let half_a_window = u64::from(crate::ng::window_coverage::WINDOW_BP / 2);
        let source = a_source(vec![Ok(observation_over(region(40, 45)))]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source and reference hold");

        assert_eq!(
            cache.reference_base_at(position_on(0, 50 + half_a_window)),
            Some(fixture_base_at(50 + half_a_window)),
            "the cover stopped short of half a window past its region, so that region's last \
             centres will be finalised after their builder has run",
        );
        assert_eq!(
            cache.reference_base_at(position_on(0, 50 + half_a_window + 1)),
            None,
            "the look-ahead is half a window and no more",
        );
    }

    /// **A cover that reached onto the next contig still leaves the region's ground in the
    /// buffer.** The contigs are read ascending, because the accumulators demand it, and the
    /// record a sample is drawn one *past* the reach can be the first of the next contig — which
    /// sorts after the region's. A cover that stopped there would leave a builder an accessor
    /// answering `None` at every position of its own region, and with the look-ahead clamped at a
    /// contig's end this is the ordinary case at every contig boundary rather than a corner.
    #[test]
    fn a_record_held_on_the_next_contig_does_not_take_the_buffer_with_it() {
        let source = a_source(vec![
            Ok(observation_over(region_on(0, 40, 45))),
            Ok(observation_over(region_on(1, 10, 10))),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region_on(0, 40, 50))
            .expect("the fixture source and reference hold");

        assert_eq!(
            cache.reference_base_at(position_on(0, 45)),
            Some(fixture_base_at(45)),
            "the cover ended holding another contig's bases, so a builder over its own region \
             finds none",
        );
        assert_eq!(
            cache.reference_base_at(position_on(1, 10)),
            None,
            "the buffer is one contig's, and the one it holds is the region's",
        );
    }

    /// **A region already past its contig's end keeps its own ground**, rather than having the
    /// clamp pull the chain back behind it. `GenomeRegion` has public fields and nothing here
    /// checks a region against the contig table, so this is representable; a chain behind the
    /// region would leave `with_observations` refusing the very builder the cover was made for.
    #[test]
    fn a_region_past_its_contigs_end_is_not_clamped_back_behind_itself() {
        let source = a_source(vec![Ok(observation_over(region(1_900, 1_900)))]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        // The fixture reference's contigs are 2,000 bases; this region ends past that.
        cache
            .cover(region(1_900, 2_400))
            .expect_err("ground past the contig's last base cannot be read");

        // The chain was not pulled back: had it been, the ground would have ended at 2,000 and
        // the fetch would have succeeded rather than refusing.
        assert_eq!(
            cache.reference_base_at(position_on(0, 1_950)),
            None,
            "a cover that could not read its ground left bases readable",
        );
    }

    /// **The look-ahead stops at the contig's end**, because the ground a cover reads runs to the
    /// chain and a fetch past a contig's last base is refused rather than truncated. The fixture
    /// reference's contigs are 2,000 bases, so a region ending inside half a window of that is
    /// the case: without the clamp the fetch asks for base 2,150 and the cover fails.
    #[test]
    fn the_look_ahead_stops_at_the_contigs_last_base() {
        let source = a_source(vec![Ok(observation_over(region(1_890, 1_895)))]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(1_890, 1_900))
            .expect("a region within half a window of the contig's end is covered, not refused");

        assert_eq!(
            cache.reference_base_at(position_on(0, 2_000)),
            Some(fixture_base_at(2_000)),
            "the ground stops short of the contig's last base",
        );
        assert_eq!(
            cache.reference_base_at(position_on(0, 2_001)),
            None,
            "the ground ran past the contig's last base",
        );
    }

    /// A cover refills the buffer rather than adding to it, so the ground of the cover before
    /// it is gone — including the whole of another contig. Reading a stale base would give a
    /// sample the wrong contig's GC at every position of the region.
    ///
    /// **The two covers overlap in coordinate**, which is what makes this about the contig
    /// rather than about the position: positions restart at 1 on every contig, so a fixture
    /// whose second cover sat at higher coordinates would come back absent under either rule
    /// and could not tell the two guards apart.
    #[test]
    fn a_position_on_the_contig_left_behind_is_not_read_from_the_new_ones_buffer() {
        let first = a_source(vec![Ok(observation_over(region_on(0, 1, 1)))]);
        let mut cache = ObservationCache::over_fixture(vec![first]);
        cache
            .cover(region_on(0, 1, 900))
            .expect("the fixture source and reference hold");
        assert_eq!(
            cache.reference_base_at(position_on(0, 500)),
            Some(fixture_base_at(500)),
        );

        cache
            .cover(region_on(1, 1, 900))
            .expect("the fixture source and reference hold");
        assert_eq!(
            cache.reference_base_at(position_on(1, 500)),
            Some(fixture_base_at(500)),
            "the second cover did not read its own contig",
        );
        assert_eq!(
            cache.reference_base_at(position_on(0, 500)),
            None,
            "a position on the contig the merge has left was answered out of the new contig's \
             buffer, which is a base read off the wrong chromosome",
        );
    }

    /// **Before the first cover there is no ground**, and asking for a base is `None` rather
    /// than an empty buffer's byte or a panic.
    #[test]
    fn a_cache_that_has_not_covered_yet_has_no_reference_base_anywhere() {
        let source = a_source(vec![Ok(observation_over(region(40, 40)))]);
        let cache = ObservationCache::over_fixture(vec![source]);
        assert_eq!(cache.reference_base_at(position_on(0, 40)), None);
    }

    /// One sample's records over `1..=bases`, one per position, each seen by `reads` reads —
    /// the shape the generic walk emits, and long enough that windows in the middle of it
    /// finalise: a centre closes only once the stream has passed 250 bases beyond it.
    fn one_record_a_position(
        bases: u64,
        reads: u32,
    ) -> std::vec::IntoIter<Result<SampleLocusObservations, SourceFailed>> {
        (1..=bases)
            .map(|at| {
                let mut record = observation_over(region(at, at));
                record.observations[0].num_obs = reads;
                Ok(record)
            })
            .collect::<Vec<_>>()
            .into_iter()
    }

    /// **The measurement, end to end through the cache.** Each of 600 positions carries three
    /// reads, and the reference repeats `ACGT`, so the window centred on 300 — the covered
    /// positions of 50 through 550 — is a mean depth of exactly 3 over 501 positions, 251 of
    /// which are `G` or `C`.
    #[test]
    fn a_builder_reads_the_window_a_sample_has_at_a_position() {
        let mut cache = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        cache
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");

        let window = cache.with_observations(region(1, 600), |window| {
            window.window_coverage_at(0, position_on(0, 300))
        });
        let window = window.expect("the sample has a window at 300");
        assert_eq!(
            window.mean_depth, 3.0,
            "the window's depth is not the three reads every position carries",
        );
        assert!(
            (f64::from(window.gc_fraction) - 251.0 / 501.0).abs() < 1e-6,
            "the window's GC is {} where the reference gives 251 of 501",
            window.gc_fraction,
        );
    }

    /// **A record is observed once, however many covers hold it.** A record is held across
    /// every cover it reaches into, and each of those covers reads ground containing it; a
    /// second observation would fold the same centre twice, so the check is that every centre a
    /// sample holds is distinct and ascending — which is also what `window_coverage_at`'s search assumes.
    #[test]
    fn a_record_held_across_two_covers_is_observed_by_one_of_them() {
        let mut over_two = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        over_two
            .cover(region(1, 300))
            .expect("the fixture source and reference hold");
        over_two
            .cover(region(301, 600))
            .expect("the fixture source and reference hold");

        let centres: Vec<GenomePosition> = over_two.samples[0]
            .window_coverage
            .finalised
            .iter()
            .map(|(centre, _)| *centre)
            .collect();
        let mut distinct = centres.clone();
        distinct.dedup();
        assert_eq!(
            centres, distinct,
            "a centre was finalised twice, so a record was observed by both covers holding it",
        );
        assert!(
            centres.windows(2).all(|pair| pair[0] < pair[1]),
            "the finalised centres are not ascending, which `window_coverage_at`'s search assumes",
        );

        // And the same ground covered in one go gives the same windows, which is what says the
        // division into covers is not visible in the measurement.
        let mut in_one = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        in_one
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");
        assert_eq!(
            over_two.samples[0].window_coverage.finalised,
            in_one.samples[0].window_coverage.finalised,
            "covering the same ground in two goes gave different windows from covering it in one",
        );
    }

    /// **The two covers finalise the same windows**, which nothing else in this module checks.
    ///
    /// The standing oracle for the merge is that its serial and parallel paths do not drift, and
    /// it is a comparison of *output* — but no output carries a window yet, so a parallel cover
    /// that observed a record twice, or missed one, would move nothing a VCF or a driver's
    /// agreement test can see. This compares the windows directly.
    #[test]
    fn the_two_covers_finalise_the_same_windows() {
        let three_covers = [region(1, 200), region(201, 400), region(401, 600)];
        let mut serially = ObservationCache::over_fixture(vec![
            one_record_a_position(600, 3),
            one_record_a_position(600, 5),
        ]);
        let mut in_parallel = ObservationCache::over_fixture(vec![
            one_record_a_position(600, 3),
            one_record_a_position(600, 5),
        ]);
        for ground in three_covers {
            serially
                .cover(ground)
                .expect("the fixture source and reference hold");
            in_parallel
                .cover_in_parallel(ground)
                .expect("the fixture source and reference hold");
        }

        for sample in 0..2 {
            assert_eq!(
                serially.samples[sample].window_coverage.finalised,
                in_parallel.samples[sample].window_coverage.finalised,
                "sample {sample}'s windows differ between the two covers",
            );
        }
        assert!(
            !serially.samples[0].window_coverage.finalised.is_empty(),
            "no window was finalised, so this compared two empty lists",
        );
    }

    /// **A sample's histogram comes out of the cache finished, and finishing is what closes the
    /// tail.** The centres in the last half-window of a sample's stream have no later position to
    /// complete them, so no cover can finalise them; only `finish` does, and this is the one place
    /// it is called. A cache that handed its sources back without it would give every sample a
    /// histogram short of its last windows.
    #[test]
    fn the_histograms_come_out_finished_with_the_tail_no_cover_could_close() {
        let mut cache = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        cache
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");
        // What the covers themselves finalised: everything up to half a window before the end.
        let finalised_by_the_covers = cache.samples[0].window_coverage.finalised.len();

        let (_sources, histograms) = cache.into_sources_and_histograms();

        let [SampleHistogram::Fitted(histogram)] = &histograms[..] else {
            panic!(
                "one sample with 600 covered positions has a fitted histogram, not {histograms:?}"
            );
        };
        assert_eq!(
            histogram.windows_folded, 600,
            "every covered position finalises one window, and the floor silenced none here",
        );
        assert!(
            u64::try_from(finalised_by_the_covers).expect("a small count")
                < histogram.windows_folded,
            "the covers finalised {finalised_by_the_covers} of the 600, so `finish` is what closed \
             the rest — and a histogram equal to what the covers had would mean it was never called",
        );
        assert_eq!(
            histogram
                .counts
                .iter()
                .map(|count| u64::from(*count))
                .sum::<u64>(),
            histogram.windows_folded,
            "the cells hold every window the histogram says it folded",
        );
    }

    /// **A sample the pass reached nothing of says which silence it was**, rather than coming back
    /// with an empty histogram that a consumer would read as a measurement of zero.
    #[test]
    fn a_sample_with_no_covered_position_comes_back_with_the_reason_it_has_no_histogram() {
        let (source, _drawn) = source_over(&[]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");

        let (_sources, histograms) = cache.into_sources_and_histograms();

        assert_eq!(histograms, vec![SampleHistogram::NoWindowFinalised]);
    }

    /// **The two evictors drop the same windows**, for the reason above: eviction touches the
    /// finalised windows, and nothing downstream can yet show that one path dropped more than the
    /// other.
    #[test]
    fn the_two_evictors_drop_the_same_windows() {
        let mut serially = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        let mut in_parallel = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        serially
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");
        in_parallel
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");
        let before = serially.samples[0].window_coverage.finalised.len();

        serially.evict_before(position_on(0, 200));
        let mut graveyard = Vec::new();
        in_parallel.evict_before_in_parallel(position_on(0, 200), &mut graveyard);

        assert_eq!(
            serially.samples[0].window_coverage.finalised,
            in_parallel.samples[0].window_coverage.finalised,
            "the two evictors left different windows behind",
        );
        assert!(
            serially.samples[0].window_coverage.finalised.len() < before,
            "nothing was evicted, so this compared two untouched lists",
        );
    }

    /// **A cover that crosses a contig reads both.** Drawing stops at the first record past the
    /// reach, and a reach on a later contig is past every position of an earlier one — so the
    /// cover that first reaches contig 1 also draws whatever the sample still had on contig 0.
    /// Those are real covered positions of that sample, and their bases are on the contig they
    /// sit on.
    #[test]
    fn a_cover_that_crosses_a_contig_observes_the_records_left_on_the_one_it_leaves() {
        let mut records: Vec<Result<SampleLocusObservations, SourceFailed>> = (1..=600)
            .map(|at| Ok(observation_over(region_on(0, at, at))))
            .collect();
        records.extend((1..=600).map(|at| Ok(observation_over(region_on(1, at, at)))));
        let mut cache = ObservationCache::over_fixture(vec![records.into_iter()]);

        // A first cover that stops well short of contig 0's records, so plenty are left for the
        // cover that moves to contig 1 to draw.
        cache
            .cover(region_on(0, 1, 100))
            .expect("the fixture source and reference hold");
        cache
            .cover(region_on(1, 1, 600))
            .expect("the fixture source and reference hold");

        let held = &cache.samples[0].window_coverage.finalised;
        assert!(
            held.iter().any(|(centre, _)| centre.contig == ContigId(0)),
            "the records left on the contig the cover left were never observed",
        );
        assert!(
            held.iter().any(|(centre, _)| centre.contig == ContigId(1)),
            "the records on the contig the cover moved to were never observed",
        );
    }

    /// Eviction drops the windows behind it and keeps the rest — by the windows' own
    /// coordinates, because a window is centred on a position and the record that position came
    /// from may already be gone.
    #[test]
    fn eviction_drops_the_windows_behind_it_and_keeps_the_rest() {
        let mut cache = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        cache
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");
        let before = cache.samples[0].window_coverage.finalised.len();
        assert!(before > 0, "the cover finalised no window to evict");

        cache.evict_before(position_on(0, 200));

        let after = &cache.samples[0].window_coverage.finalised;
        assert!(
            after.len() < before,
            "eviction dropped no window, so the deque grows with the contig",
        );
        assert!(
            after
                .iter()
                .all(|(centre, _)| centre.position >= Position(200)),
            "a window behind the evicted position survived",
        );
    }

    /// A sample with no record at a position has no window there, and that is an ordinary
    /// answer: the filter that reads this skips a sample it has no window for.
    #[test]
    fn a_sample_with_no_record_at_a_position_has_no_window_there() {
        let mut cache = ObservationCache::over_fixture(vec![one_record_a_position(600, 3)]);
        cache
            .cover(region(1, 600))
            .expect("the fixture source and reference hold");

        let window = cache.with_observations(region(1, 600), |window| {
            window.window_coverage_at(0, position_on(0, 599))
        });
        assert_eq!(
            window, None,
            "a centre the stream has not passed by half a window is not finalised, so there is \
             no window to read there yet",
        );
    }

    /// **Every record the cache holds has a base**, which is the property the whole step
    /// exists to provide and the one the next step reads. Two shapes make it non-trivial, and
    /// both occur on ordinary input: each sample is drawn **one record past** the reach — the
    /// only way to know a record is beyond it is to draw it — and a record drawn by an earlier
    /// cover is still held and handed out again, with its head before this region's start.
    #[test]
    fn every_position_of_every_held_record_has_a_reference_base() {
        // 30..44 is drawn by an earlier cover and still reaches into this one; 60..60 is the
        // record this cover draws past its own reach.
        let source = a_source(vec![
            Ok(observation_over(region(30, 44))),
            Ok(observation_over(region(45, 45))),
            Ok(observation_over(region(60, 60))),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(20, 35))
            .expect("the fixture source and reference hold");
        cache
            .cover(region(40, 50))
            .expect("the fixture source and reference hold");

        for position in [30, 40, 44, 45, 50, 60] {
            assert_eq!(
                cache.reference_base_at(position_on(0, position)),
                Some(fixture_base_at(position)),
                "the cache holds a record covering {position} and no base for it",
            );
        }
    }

    /// The parallel cover reads its ground too. **It is the production path on any multi-core
    /// machine** — the driver picks it whenever the pool has more than one thread — so a fetch
    /// that ran only in the serial cover would leave the whole measurement absent from a real
    /// run and present in the oracle.
    #[test]
    fn the_parallel_cover_holds_the_reference_bases_over_the_ground_it_drew() {
        let source = a_source(vec![Ok(observation_over(region(48, 60)))]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover_in_parallel(region(40, 50))
            .expect("the fixture source and reference hold");

        assert_eq!(
            cache.reference_base_at(position_on(0, 45)),
            Some(fixture_base_at(45)),
        );
        assert_eq!(
            cache.reference_base_at(position_on(0, 60)),
            Some(fixture_base_at(60)),
            "the parallel cover did not read as far as it drew",
        );
    }

    /// **The merge lets go of what it has walked past.** A windowed reader extends its buffer
    /// while each request lands near the last, which is exactly the merge's pattern, and shrinks
    /// it only when told — so a merge that never told it would end a contig holding every base
    /// it had covered. An in-memory reference has nothing to release, so what this pins is that
    /// the release is asked for, and from where.
    #[test]
    fn a_cover_releases_the_reference_bases_the_merge_has_walked_past() {
        let reference = ReferenceCountingItsReleases::new();
        let source = a_source(vec![
            Ok(observation_over(region(40, 40))),
            Ok(observation_over(region(140, 140))),
        ]);
        let mut cache = ObservationCache::over(vec![source], Box::new(Arc::clone(&reference)));
        cache
            .cover(region(40, 50))
            .expect("the fixture source and reference hold");
        cache.evict_before(position_on(0, 130));
        cache
            .cover(region(130, 150))
            .expect("the fixture source and reference hold");

        let released = reference
            .released_before
            .lock()
            .expect("no test holds this lock");
        assert_eq!(
            *released,
            vec![40, 130],
            "each cover releases the bases before the ground it went on to read",
        );
    }

    /// A failed fetch names the ground the cover held, not the region it was asked for — the
    /// ground is what tells an operator which stretch of FASTA went unreadable, and it is the
    /// field the run's own error prints.
    #[test]
    fn a_failed_fetch_names_the_ground_the_cover_held() {
        let records: Vec<Result<SampleLocusObservations, ReferenceFailureRecorded>> =
            vec![Ok(observation_over(region_on(9, 48, 60)))];
        let mut cache = ObservationCache::over_fixture(vec![records.into_iter()]);
        let refused = cache
            .cover(region_on(9, 40, 50))
            .expect_err("contig 9 is past the fixture reference's four");
        // **Half a window past the region, not the observation's own last base**: the
        // look-ahead is what the cover drew to, so it is what the failure has to name.
        let half_a_window = u64::from(crate::ng::window_coverage::WINDOW_BP / 2);
        assert_eq!(
            refused,
            ReferenceFailureRecorded(Some(region_on(9, 40, 50 + half_a_window))),
            "the failure names other ground than the cover held",
        );
    }

    /// **A cover that could not read its ground leaves none readable**, rather than answering
    /// out of the cover before it. The type says a cover can be made again, so a caller that
    /// catches a failure and retries must not be handed the previous cover's bases as this
    /// one's.
    #[test]
    fn a_failed_fetch_leaves_no_ground_readable_and_a_retry_restores_it() {
        let source = a_source(vec![
            Ok(observation_over(region_on(0, 40, 40))),
            Ok(observation_over(region_on(0, 60, 60))),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region_on(0, 40, 50))
            .expect("the fixture source and reference hold");
        assert!(cache.reference_base_at(position_on(0, 45)).is_some());

        cache
            .cover(region_on(9, 100, 110))
            .expect_err("contig 9 is past the fixture reference's four");
        assert_eq!(
            cache.reference_base_at(position_on(0, 45)),
            None,
            "after a failed fetch the previous cover's bases are still readable",
        );

        cache
            .cover(region_on(0, 51, 70))
            .expect("the fixture source and reference hold");
        assert_eq!(
            cache.reference_base_at(position_on(0, 60)),
            Some(fixture_base_at(60)),
            "a cover made again after a failure did not restore the ground",
        );
    }

    /// A reference that cannot serve the ground ends the cover, naming the ground rather than
    /// the position inside it that failed. **It is not absorbed**: one failed fetch costs every
    /// sample its coverage over that stretch, and a run carrying on would report windows whose
    /// denominators quietly excluded it.
    #[test]
    fn a_reference_that_cannot_serve_the_cover_ground_ends_the_cover_naming_it() {
        // Contig 9 is past the fixture reference's four, so the fetch cannot be served.
        let source = a_source(Vec::new());
        let mut cache = ObservationCache::over_fixture(vec![source]);
        let refused = cache
            .cover(region_on(9, 40, 50))
            .expect_err("a reference with no contig 9 cannot serve this cover");
        assert_eq!(
            refused,
            SourceFailed("the fixture reference could not serve the cover's ground"),
            "the cover did not fail through the reference conversion",
        );
    }

    /// A source over the given regions, with the counter the test reads afterwards.
    fn source_over(regions: &[GenomeRegion]) -> (CountingSource, Rc<Cell<usize>>) {
        let drawn = Rc::new(Cell::new(0));
        let items: Vec<_> = regions
            .iter()
            .map(|region| Ok(observation_over(*region)))
            .collect();
        (
            CountingSource {
                remaining: items.into_iter(),
                drawn: Rc::clone(&drawn),
            },
            drawn,
        )
    }

    /// The regions each sample's window hands out over `span`, in sample order.
    fn handed_out<S>(cache: &ObservationCache<S>, span: GenomeRegion) -> Vec<Vec<GenomeRegion>> {
        cache.with_observations(span, |per_sample| {
            per_sample
                .observations
                .expect("the fixture's sources build every record")
                .iter()
                .map(|observations| {
                    observations
                        .iter()
                        .map(|observation| observation.region)
                        .collect()
                })
                .collect()
        })
    }

    /// The plain case: covering a region draws the observations over it, and they come back.
    #[test]
    fn covering_a_region_draws_the_observations_over_it() {
        let (source, _drawn) = source_over(&[region(45, 45), region(48, 48)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(45, 45), region(48, 48)]],
        );
    }

    /// **A window covering a region also holds an observation that started before it and
    /// reaches in** — the plan's own D1 case. The deletion opening at 30 is what a locus over
    /// this ground chains back to, so a cache that handed out only what begins inside the
    /// region would let a builder claim a locus whose first bases it never saw.
    #[test]
    fn an_observation_that_began_earlier_and_reaches_in_is_handed_out() {
        let (source, _drawn) = source_over(&[region(30, 44), region(45, 45)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(30, 44), region(45, 45)]],
            "30–44 ends inside the region, so it is part of what a locus there can chain to",
        );
    }

    /// And one that ends before the span is not handed out, though the forward reader had to
    /// pass over it. That is the prefix a builder does not pay for.
    ///
    /// **The two assertions sit either side of the one base that decides it.** 20–39 ends one
    /// base before a span opening at 40 and is trimmed away; asked again over a span opening
    /// on 39 — its own last base — it comes back, because a locus opening there chains into
    /// it. An edge drawn one base further right would lose the locus's true first position.
    #[test]
    fn an_observation_that_ends_before_the_span_is_drawn_but_not_handed_out() {
        let (source, drawn) = source_over(&[region(20, 39), region(41, 45), region(46, 46)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(41, 45), region(46, 46)]],
            "20–39 ends one base before the span and can reach no locus inside it",
        );
        assert_eq!(
            handed_out(&cache, region(39, 50)),
            vec![vec![region(20, 39), region(41, 45), region(46, 46)]],
            "and it comes back the moment the span opens on its own last base",
        );
        assert_eq!(
            drawn.get(),
            3,
            "it was still read: a forward reader passes over it to reach 41",
        );
    }

    /// **The chain's reach follows a widening in a later sample.** Sample 1's deletion, from
    /// 48 to 70, is what makes sample 0's observations part of a locus opening inside the
    /// region — and only a second sweep, after that deletion has moved the reach, draws them.
    ///
    /// **The second of sample 0's records is what makes the sweep visible.** One observation
    /// beyond the reach is held whatever the cover does, because the draw that discovers it
    /// keeps it; what a single sweep loses is the observation *after* it. Measured by the D1
    /// review: with only the first, this test passes under a single-sweep cover.
    #[test]
    fn the_chain_reach_follows_a_widening_in_a_later_sample() {
        let (far, far_drawn) = source_over(&[region(60, 60), region(65, 65)]);
        let (widening, _widening_drawn) = source_over(&[region(48, 70)]);
        let mut cache = ObservationCache::over_fixture(vec![far, widening]);

        cache
            .cover(region(40, 55))
            .expect("the fixture source holds");

        assert_eq!(far_drawn.get(), 2, "a single sweep would draw only 60");
        assert_eq!(
            handed_out(&cache, region(40, 55)),
            vec![vec![region(60, 60), region(65, 65)], vec![region(48, 70)]],
            "the deletion reaching to 70 is what draws 65 into the window; 60 was already \
             kept by the draw that discovered it",
        );
    }

    /// **The fixpoint needs a third sweep here, and a cover that stopped after two would hand
    /// out a window one observation short.** Sample 2's deletion carries the reach to 70; only
    /// then is sample 1's 60–120 inside it, which carries it to 120; only then is sample 0's
    /// 100–300 inside it, which carries it to 300 — and 260 is in the window only because of
    /// that third widening.
    ///
    /// Found by the D1 review: capping the sweeps at two passed every other test in this file,
    /// and disagreed with a whole-stretch builder on 410 of 600 random layouts.
    #[test]
    fn a_chain_that_needs_a_third_sweep_is_drawn_whole() {
        let (far, _far_drawn) = source_over(&[region(100, 300)]);
        let (middle, _middle_drawn) =
            source_over(&[region(60, 120), region(250, 250), region(260, 260)]);
        let (opening, _opening_drawn) = source_over(&[region(50, 70)]);
        let mut cache = ObservationCache::over_fixture(vec![far, middle, opening]);

        cache
            .cover(region(40, 55))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 55)),
            vec![
                vec![region(100, 300)],
                vec![region(60, 120), region(250, 250), region(260, 260)],
                vec![region(50, 70)],
            ],
            "260 is in the window only because the third sweep carried the reach to 300",
        );
    }

    /// **An observation beginning on the reach's own base is inside it, and what it reaches is
    /// drawn too.** This is the same `<=` the closer chains on (spec §4.1), and it is the
    /// boundary that decides whether the window is complete: the third sample's second
    /// observation is in the window only because sample 0's record — which begins exactly
    /// where sample 1's deletion ends — carried the reach out to 200.
    ///
    /// **The observation at the boundary alone would not show it.** The one draw that
    /// discovers an observation is beyond the reach leaves it in the window either way, so a
    /// cache that stopped at the boundary would still hand out sample 0's record; what it
    /// would not hold is what that record reaches, which is why the late sample carries two.
    #[test]
    fn an_observation_beginning_on_the_reach_is_drawn_and_widens_it() {
        let (widened, widened_drawn) = source_over(&[region(70, 200), region(500, 500)]);
        let (widening, _widening_drawn) = source_over(&[region(48, 70)]);
        let (late, late_drawn) = source_over(&[region(150, 150), region(160, 160)]);
        let mut cache = ObservationCache::over_fixture(vec![widened, widening, late]);

        cache
            .cover(region(40, 55))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 55)),
            vec![
                vec![region(70, 200), region(500, 500)],
                vec![region(48, 70)],
                vec![region(150, 150), region(160, 160)],
            ],
            "70 begins on the reach the deletion left at 70, so its own reach to 200 is what \
             puts 160 in the window",
        );
        assert_eq!(widened_drawn.get(), 2, "500 is the one draw past the chain");
        assert_eq!(
            late_drawn.get(),
            2,
            "both of the late sample's observations were drawn, not merely the first",
        );
    }

    /// **The window overshoots by exactly one observation per sample, and stops.** The
    /// observation at 800 is beyond the reach and had to be drawn to find that out; the one at
    /// 900 was never read at all.
    ///
    /// 800 is past the 750 the region's own look-ahead reaches
    /// ([`where_the_chain_reaches_before_any_sweep`](ObservationCache::where_the_chain_reaches_before_any_sweep)),
    /// so the draw that ends the chain is the second record and not the third.
    #[test]
    fn drawing_stops_one_observation_past_the_reach() {
        let (source, drawn) = source_over(&[region(450, 450), region(800, 800), region(900, 900)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        cache
            .cover(region(400, 500))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            2,
            "450 is inside the region and 800 is the one draw that says the chain has ended",
        );
        assert_eq!(
            handed_out(&cache, region(400, 500)),
            vec![vec![region(450, 450), region(800, 800)]],
            "the overshoot is held rather than thrown away, and a later cover reconsiders it",
        );
    }

    /// A second cover carries on from where the first stopped, and re-reads nothing.
    #[test]
    fn a_second_cover_carries_on_without_re_reading() {
        let (source, drawn) = source_over(&[region(45, 45), region(60, 60), region(80, 80)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache
            .cover(region(51, 70))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            3,
            "the second cover drew only 80 — 60 was already held from the first",
        );
        assert_eq!(
            handed_out(&cache, region(51, 70)),
            vec![vec![region(60, 60), region(80, 80)]],
            "and 45 is behind the span, so the second builder is not handed it",
        );
    }

    /// **Eviction drops nothing a live region can still reach** — the plan's second D1 case.
    /// The observation at 20–21 cannot touch the ground from 40 on and goes.
    ///
    /// **The second observation ends exactly on the evicted base**, which is the case the rule
    /// turns on: 25–40 still covers position 40, so a locus opening there chains back to 25,
    /// and a cache that dropped it would let that locus appear to open at 40 — claimed and
    /// built by a builder that never saw its first fifteen bases.
    ///
    /// The window is asked over ground reaching back to base 1, so a survivor that eviction
    /// had merely hidden would show here.
    #[test]
    fn eviction_keeps_what_reaches_past_the_evicted_point() {
        let (source, _drawn) = source_over(&[
            region(20, 21),
            region(25, 40),
            region(41, 44),
            region(45, 45),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        cache.evict_before(position_on(0, 40));

        assert_eq!(
            handed_out(&cache, region(1, 50)),
            vec![vec![region(25, 40), region(41, 44), region(45, 45)]],
            "20–21 is gone; 25–40 ends on the evicted base and stays",
        );
    }

    /// **Every sample's window is evicted, not merely the first.** Both samples here have a
    /// record behind the evicted point and one ahead of it, so a drop that reached only sample
    /// 0 would leave sample 1's 22–23 in the window — and the cache's whole memory bound is
    /// that this loop runs to the end (spec §8). Nothing in the output would show it.
    #[test]
    fn eviction_drops_from_every_sample_not_only_the_first() {
        let (first, _first_drawn) = source_over(&[region(20, 21), region(45, 45)]);
        let (second, _second_drawn) = source_over(&[region(22, 23), region(46, 46)]);
        let mut cache = ObservationCache::over_fixture(vec![first, second]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        cache.evict_before(position_on(0, 40));

        assert_eq!(
            handed_out(&cache, region(1, 50)),
            vec![vec![region(45, 45)], vec![region(46, 46)]],
            "both samples lost the record that ends before 40",
        );
    }

    /// **A contig boundary is past every position on the contig before it.** The held record
    /// at 900 on contig 0 is a much higher *number* than the position evicted at, 10 on contig
    /// 1, and it still goes: nothing on contig 1 can chain back to it.
    #[test]
    fn eviction_at_a_later_contig_drops_the_previous_contigs_window() {
        let (source, _drawn) = source_over(&[region_on(0, 900, 900), region_on(1, 45, 45)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region_on(1, 40, 50))
            .expect("the fixture source holds");

        cache.evict_before(position_on(1, 10));

        assert_eq!(
            handed_out(&cache, region_on(0, 1, 999)),
            vec![vec![region_on(1, 45, 45)]],
            "contig 0's record is behind the evicted point though 900 is ahead of 10",
        );
    }

    /// Eviction leaves the window usable: the next cover draws forward from where it stood,
    /// and what survived is still handed out.
    #[test]
    fn a_cover_after_an_eviction_draws_forward_and_keeps_the_survivor() {
        let (source, drawn) = source_over(&[region(20, 21), region(30, 44), region(60, 60)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 40));

        cache
            .cover(region(51, 70))
            .expect("the fixture source holds");

        assert_eq!(drawn.get(), 3, "everything the fixture holds, and no more");
        assert_eq!(
            handed_out(&cache, region(1, 70)),
            vec![vec![region(30, 44), region(60, 60)]],
            "the survivor of the eviction is still there, and the new ground was drawn",
        );
    }

    /// **The window's own observations are re-read at every cover, so an eviction leaves no
    /// mark to correct.** After the eviction, 300–800 is still what a locus in 510–540 chains
    /// through, so it must widen the reach again.
    ///
    /// **The positions are wide apart because the look-ahead is 250 bases**
    /// ([`where_the_chain_reaches_before_any_sweep`](ObservationCache::where_the_chain_reaches_before_any_sweep)): a chain that
    /// reaches less than half a window past its region widens nothing a cover had not already
    /// asked for, and the widening this test is about would be invisible. 300–800 reaches 800,
    /// which is past the second cover's own 790.
    #[test]
    fn a_survivor_of_an_eviction_still_widens_the_next_reach() {
        let (source, drawn) = source_over(&[
            region(200, 210),
            region(300, 800),
            region(900, 900),
            region(1_500, 1_500),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(400, 500))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 400));

        cache
            .cover(region(510, 540))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            3,
            "300–800 reaches to 800, and 900 is the draw that ends the chain there",
        );
        assert_eq!(
            handed_out(&cache, region(1, 800)),
            vec![vec![region(300, 800), region(900, 900)]],
        );
    }

    /// Two observations evicted at once, with more still held — the case where an index
    /// remembered across covers would be not merely stale but out of range.
    #[test]
    fn two_evicted_at_once_leave_the_window_sound() {
        // Wide apart for the reason `a_survivor_of_an_eviction_still_widens_the_next_reach` is:
        // a cover draws half a window past its region whatever chains there, so a fixture
        // packed inside 250 bases is drawn whole by the first cover and has nothing to evict.
        let (source, _drawn) = source_over(&[
            region(100, 110),
            region(200, 210),
            region(300, 800),
            region(1_000, 1_500),
            region(1_600, 1_600),
            region(1_800, 1_800),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(400, 500))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 400));

        cache
            .cover(region(510, 540))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(1, 800)),
            vec![vec![region(300, 800), region(1_000, 1_500)]],
            "both of the evicted pair are gone and the window behind them is intact",
        );
    }

    /// **A region whose ends are the wrong way round covers the ground it names**, rather than
    /// a reach before its own first base. `GenomeRegion` has public fields and no constructor
    /// ordering them, and a reach at the lower end would end the cover before the ground was
    /// drawn — a window short of what a locus there reaches, which is a wrong answer and not a
    /// failure. The same defence `SampleLocusObservations::reach` makes, and the second
    /// assertion is what says `with_observations` reads an inverted span the same way.
    #[test]
    fn a_region_given_end_first_still_covers_its_ground() {
        let (source, drawn) = source_over(&[region(45, 45), region(60, 60)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        let inverted = GenomeRegion {
            contig: ContigId(0),
            start: Position(50),
            end: Position(40),
        };

        cache.cover(inverted).expect("the fixture source holds");

        assert_eq!(drawn.get(), 2, "45 was drawn, and 60 ended the chain");
        assert_eq!(
            handed_out(&cache, inverted),
            vec![vec![region(45, 45), region(60, 60)]],
            "asked back over the same inverted region, the window is the ground it names",
        );
        assert_eq!(
            cache.reference_base_at(position_on(0, 40)),
            Some(fixture_base_at(40)),
            "the ground read does not start at the inverted region's lower bound",
        );
    }

    /// A contig boundary is beyond every reach on the contig before it, whatever the positions
    /// say: covering a region on contig 1 draws past contig 0's observations and stops at
    /// contig 2's.
    #[test]
    fn a_cover_stops_at_the_next_contig() {
        let (source, drawn) = source_over(&[
            region_on(0, 90, 90),
            region_on(1, 45, 45),
            region_on(2, 10, 10),
            region_on(2, 20, 20),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        cache
            .cover(region_on(1, 40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            drawn.get(),
            3,
            "contig 0's observation, contig 1's, and the one draw that ends the chain",
        );
        assert_eq!(
            handed_out(&cache, region_on(1, 40, 50)),
            vec![vec![region_on(1, 45, 45), region_on(2, 10, 10)]],
            "contig 0's observation is behind the span even though position 90 is ahead of 40",
        );
    }

    /// A source's failure ends the cover and comes back unchanged.
    #[test]
    fn a_failing_source_ends_the_cover() {
        let drawn = Rc::new(Cell::new(0));
        let source = CountingSource {
            remaining: vec![
                Ok(observation_over(region(45, 45))),
                Err(SourceFailed("the block would not decode")),
            ]
            .into_iter(),
            drawn: Rc::clone(&drawn),
        };
        let mut cache = ObservationCache::over_fixture(vec![source]);

        let outcome = cache.cover(region(40, 50));

        assert_eq!(outcome, Err(SourceFailed("the block would not decode")));
        assert_eq!(drawn.get(), 2, "it failed on its second observation");
    }

    /// **A window over ground no cover reached is refused**, which is the state a failed cover
    /// leaves behind. Handing one out instead would give a builder a window short of what a
    /// locus there reaches — a locus closed over the wrong ground, and no failure anywhere.
    #[test]
    #[should_panic(expected = "cover has only reached")]
    fn a_window_over_ground_no_cover_reached_is_refused() {
        // 1,000 onwards is past the look-ahead's own half window as well as past the chain, so
        // the ground below is ground no cover reached however far the region alone would carry.
        let (source, _drawn) = source_over(&[region(450, 450), region(1_200, 1_200)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(400, 500))
            .expect("the fixture source holds");

        let _ = handed_out(&cache, region(1_000, 1_300));
    }

    /// And a window is refused for reaching **past** covered ground, not merely for starting
    /// past it: the span below opens at 450, which was drawn, and ends at 1,300, which was not.
    /// That is the shape a builder handed too wide a region takes, and the loci it would lose
    /// are the ones in the ground the reader never reached.
    #[test]
    #[should_panic(expected = "cover has only reached")]
    fn a_window_reaching_past_the_covered_ground_is_refused() {
        let (source, _drawn) = source_over(&[region(450, 450), region(1_200, 1_200)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);
        cache
            .cover(region(400, 500))
            .expect("the fixture source holds");

        let _ = handed_out(&cache, region(450, 1_300));
    }

    /// **A cover that failed can be made again**, and the sample whose source failed is drawn
    /// on from where it stopped — which is why nothing latches on `Err`. The reach is rebuilt
    /// by re-reading the held window, so the retry still reaches 70 without re-reading the
    /// source, and what was drawn before the failure is still there.
    #[test]
    fn a_cover_can_be_made_again_after_a_failure() {
        let (steady, _steady_drawn) = source_over(&[region(45, 45)]);
        let recovering = CountingSource {
            remaining: vec![
                Ok(observation_over(region(48, 70))),
                Err(SourceFailed("the block would not decode")),
                Ok(observation_over(region(60, 60))),
            ]
            .into_iter(),
            drawn: Rc::new(Cell::new(0)),
        };
        let mut cache = ObservationCache::over_fixture(vec![steady, recovering]);

        let first = cache.cover(region(40, 50));
        assert_eq!(first, Err(SourceFailed("the block would not decode")));

        cache.cover(region(40, 50)).expect("the retry holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![vec![region(45, 45)], vec![region(48, 70), region(60, 60)]],
            "the failed sample was drawn on, and the re-read window still reached 70",
        );
    }

    /// A spent source is not drawn from again. An iterator is free to yield `Some` after a
    /// `None`, and a cache that kept asking would draw observations behind the window's own
    /// reach, and so silently out of order.
    #[test]
    fn a_spent_source_is_not_drawn_from_again() {
        struct Resurrecting {
            polls: Rc<Cell<usize>>,
        }

        impl Iterator for Resurrecting {
            type Item = Result<SampleLocusObservations, SourceFailed>;

            fn next(&mut self) -> Option<Self::Item> {
                self.polls.set(self.polls.get() + 1);
                match self.polls.get() {
                    1 => Some(Ok(observation_over(region(45, 45)))),
                    2 => None,
                    _ => Some(Ok(observation_over(region(10, 10)))),
                }
            }
        }

        let polls = Rc::new(Cell::new(0));
        let mut cache = ObservationCache::over_fixture(vec![Resurrecting {
            polls: Rc::clone(&polls),
        }]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache
            .cover(region(51, 60))
            .expect("the fixture source holds");

        assert_eq!(
            polls.get(),
            2,
            "the second cover did not poll a spent source"
        );
        assert_eq!(
            handed_out(&cache, region(1, 60)),
            vec![vec![region(45, 45)]]
        );
    }

    /// A source whose observations go backwards is refused rather than quietly ending the
    /// cover where the walk cannot see it.
    #[test]
    #[should_panic(expected = "not in coordinate order")]
    fn a_source_that_goes_backwards_is_refused() {
        let (source, _drawn) = source_over(&[region(45, 45), region(30, 30)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        let _ = cache.cover(region(40, 50));
    }

    /// **A source whose records overlap is refused by the measurement's own order check, and it
    /// is refused for the right reason.**
    ///
    /// Two things had to be right for that. The ground a cover reads on a contig it is *leaving*
    /// has no chain reach to lean on, so it is the furthest reach of what that contig holds —
    /// taking the last held record's would be right only while reaches are monotone, which
    /// nothing checks, since `draw_next` compares starts. Here a repeat tract reaches to 1,500
    /// and the record after it stops at 400, so a ground ending at 400 leaves the tract's own
    /// positions with no base, and the cover panics **blaming its own ground computation** for
    /// what is a fact about the source. With the ground right, the stream reaches the
    /// accumulator and its order check says what is really wrong — and says it in release, which
    /// is the profile this repository ships.
    ///
    /// **Nothing mints such a pair today** — the generic walk's records are disjoint — but a run
    /// over stored files reads them from a file, and `draw_next`'s own comment already says an
    /// out-of-order source is a fact about the file rather than a bug in this crate.
    #[test]
    #[should_panic(expected = "observed out of order")]
    fn a_source_whose_records_overlap_is_refused_naming_the_order_and_not_the_ground() {
        let mut tract = observation_over(region_on(0, 300, 1_500));
        tract.kind = LocusKind::Ssr(crate::ng::locus_generation::SsrDetail {
            motif: crate::ng::types::Motif::new(b"AT").expect("a two-base motif"),
            left_flank: Box::from(&b"CCC"[..]),
            right_flank: Box::from(&b"GGG"[..]),
        });
        // A cover whose region is on contig 1 draws every contig-0 record the sample still has,
        // all of them unobserved, and then reads contig 0 as a contig it is leaving.
        let source = a_source(vec![
            Ok(tract),
            Ok(observation_over(region_on(0, 400, 400))),
            Ok(observation_over(region_on(1, 10, 10))),
        ]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        let _ = cache.cover(region_on(1, 10, 20));
    }

    /// **Two records at one start are refused**, because the measurement's cursor cannot tell
    /// them apart: it is a *start*, so the second falls out of the not-yet-observed suffix the
    /// moment a cover boundary separates the two, and never reaches the accumulator at all —
    /// leaving that sample short of one record's positions with nothing to say so.
    #[test]
    #[should_panic(expected = "does not begin after")]
    fn a_source_with_two_records_at_one_start_is_refused() {
        let (source, _drawn) = source_over(&[region(45, 45), region(45, 47)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        let _ = cache.cover(region(40, 50));
    }

    /// And one that goes back to an **earlier contig** is refused too, though its positions
    /// rise. `GenomePosition`'s order is contig first, so this is the same check — but a
    /// source merging per-contig files is where it actually fires, and a comparison made only
    /// within a contig would let it through.
    #[test]
    #[should_panic(expected = "does not begin after")]
    fn a_source_that_goes_back_a_contig_is_refused() {
        let (source, _drawn) = source_over(&[region_on(1, 45, 45), region_on(0, 90, 90)]);
        let mut cache = ObservationCache::over_fixture(vec![source]);

        let _ = cache.cover(region_on(1, 40, 50));
    }

    /// A cohort of no samples covers nothing and hands out nothing, rather than failing:
    /// refusing a zero-sample run happens where the run is configured (spec §7.2).
    #[test]
    fn a_cache_over_no_samples_covers_nothing() {
        let mut cache: ObservationCache<CountingSource> =
            ObservationCache::over_fixture(Vec::new());

        cache.cover(region(40, 50)).expect("nothing can fail");

        assert!(handed_out(&cache, region(40, 50)).is_empty());
    }

    /// A sample that has run out does not stop the others: the cohort's ground is still
    /// covered, and the spent sample simply hands out nothing.
    #[test]
    fn a_sample_that_has_run_out_does_not_stop_the_cover() {
        let (spent, _spent_drawn) = source_over(&[region(10, 10)]);
        let (covering, _covering_drawn) = source_over(&[region(45, 45)]);
        let mut cache = ObservationCache::over_fixture(vec![spent, covering]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");

        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![Vec::new(), vec![region(45, 45)]],
        );
    }

    /// A sample with no records at all is not a special case either — an ordinary thing at 63
    /// samples, where one accession may have nothing on a contig.
    #[test]
    fn covering_a_region_with_an_empty_source_hands_out_nothing() {
        let (empty, empty_drawn) = source_over(&[]);
        let (covering, _covering_drawn) = source_over(&[region(45, 45)]);
        let mut cache = ObservationCache::over_fixture(vec![empty, covering]);

        cache
            .cover(region(40, 50))
            .expect("the fixture source holds");
        cache.evict_before(position_on(0, 40));

        assert_eq!(
            empty_drawn.get(),
            0,
            "an empty source yields nothing to count",
        );
        assert_eq!(
            handed_out(&cache, region(40, 50)),
            vec![Vec::new(), vec![region(45, 45)]],
        );
    }

    /// **The held count is the whole cohort's, not the first sample's**, and it is zero before
    /// anything is drawn. A count over one sample would under-report the cache by the cohort
    /// size, which is what this accessor exists to measure (spec §8).
    #[test]
    fn the_held_count_sums_every_samples_window() {
        let (first, _first_drawn) = source_over(&[region(45, 45)]);
        let (second, _second_drawn) = source_over(&[region(46, 46)]);
        let mut cache = ObservationCache::over_fixture(vec![first, second]);

        assert_eq!(
            cache.held_observations_len(),
            0,
            "nothing has been drawn yet"
        );
        cache
            .cover(region(40, 50))
            .expect("the fixture sources hold");

        assert_eq!(
            cache.held_observations_len(),
            2,
            "one from each sample, not one in all",
        );
    }

    /// And a cache over no samples holds nothing rather than failing.
    #[test]
    fn the_held_count_is_zero_over_no_samples() {
        let cache: ObservationCache<CountingSource> = ObservationCache::over_fixture(Vec::new());

        assert_eq!(cache.held_observations_len(), 0);
    }

    /// **The building regions tile the analysed ground exactly** — adjacent, in genome order,
    /// none of them empty, and the last one clamped to the analysed region's own last base.
    ///
    /// **This is what a merge's output cannot check.** Dividing the ground changes nothing a
    /// caller sees, so a driver that ignored the width and built each analysed region as one
    /// would pass every byte-identity test in `serial.rs` while losing the whole point of the
    /// cache. Fifty bases divided five ways, including two widths that do not divide it.
    #[test]
    fn the_building_regions_tile_the_analysed_ground_exactly() {
        for bases in [1, 7, 20, 50, 600] {
            let divided: Vec<_> = building_regions_of(region(1, 50), width(bases)).collect();

            assert_eq!(divided[0].start, Position(1), "at {bases} bases");
            assert_eq!(
                divided.last().expect("at least one region").end,
                Position(50),
                "the last region is clamped to the analysed ground, at {bases} bases",
            );
            for pair in divided.windows(2) {
                assert_eq!(
                    pair[1].start.0,
                    pair[0].end.0 + 1,
                    "the regions are adjacent, at {bases} bases",
                );
            }
            for one in &divided {
                assert!(one.start <= one.end, "no empty region, at {bases} bases");
                // The width of an inclusive region, which is what the `+ 1` is for.
                let width = one.end.0 - one.start.0 + 1;
                assert!(
                    width <= u64::from(bases),
                    "no region wider than the width asked for, at {bases} bases",
                );
            }
            assert_eq!(
                divided.len(),
                usize::try_from(50u64.div_ceil(u64::from(bases))).expect("a small count"),
                "as many regions as the width divides fifty bases into",
            );
        }
    }

    /// A building region ending on the last base of the coordinate space has no successor, and
    /// the division stops rather than wrapping to zero and dividing the genome for ever.
    ///
    /// **Bounded, because the failure this guards against is unbounded.** With the wrap in
    /// place the iterator never ends, and an unbounded `collect` would take the whole test
    /// binary down with it rather than print a difference.
    #[test]
    fn the_division_stops_at_the_coordinate_ceiling() {
        let at_the_ceiling = GenomeRegion {
            contig: ContigId(0),
            start: Position(u64::MAX - 1),
            end: Position(u64::MAX),
        };

        let divided: Vec<_> = building_regions_of(at_the_ceiling, width(20))
            .take(3)
            .collect();

        assert_eq!(
            divided,
            vec![at_the_ceiling],
            "one region, and no second lap"
        );
    }

    /// **An analysed region whose ends are the wrong way round is read as the ground it
    /// names** — belt-and-braces behind the drivers' own refusal
    /// (`super::serial::refuse_malformed_analysed_regions`), and the reason the division cannot
    /// yield one empty inverted region instead.
    #[test]
    fn the_division_reads_an_inverted_region_as_the_ground_it_names() {
        let inverted = GenomeRegion {
            contig: ContigId(0),
            start: Position(50),
            end: Position(1),
        };

        assert_eq!(
            building_regions_of(inverted, width(20)).collect::<Vec<_>>(),
            vec![region(1, 20), region(21, 40), region(41, 50)],
        );
    }

    /// **The window is sufficient**, which is the property this whole file exists for and the
    /// one no single fixture can state: a builder fed from the cache closes exactly the loci it
    /// would close given every observation in the genome
    /// ([`build_region`](super::super::build::build_region)'s input contract).
    ///
    /// Randomised over 200 layouts from a seeded generator, so a disagreement is reproducible
    /// from the seed it prints. Records are disjoint and ascending within a sample, which is
    /// what the generator mints and what `build_region` asserts; one observation in ten is a
    /// deletion up to 150 bases wide, so chains across samples are common; the building regions
    /// have gaps, so the cache is asked to jump forward over ground nobody builds.
    ///
    /// Written by the D1 review, which measured it: it agrees on 600 of 600 layouts as the code
    /// stands, and disagrees on 410 of 600 under a cover capped at two sweeps.
    ///
    /// **The evict/cover/hand-out loop is written out here rather than calling the serial
    /// driver**, deliberately: this test leaves *gaps* between its builder regions, and the
    /// driver tiles the analysed ground without gaps. Unifying the two would quietly drop the
    /// gaps, which are what make the cache jump forward over ground nobody builds.
    #[test]
    fn a_builder_fed_from_the_cache_closes_the_loci_a_whole_stretch_would() {
        use super::super::build::build_region;
        use super::super::{MaxCohortLocusSpan, MinAltReads};

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

        let max_span = MaxCohortLocusSpan::DEFAULT;
        let keep = MinAltReads::DEFAULT;
        let ground_end = 400u64;
        let mut disagreements = Vec::new();

        for seed in 0..200u64 {
            let mut draw = Seeded(seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xDEAD_BEEF);
            let samples = 1 + draw.next(8) as usize;
            let contigs = 1 + draw.next(2) as u32;

            let mut layouts: Vec<Vec<GenomeRegion>> = Vec::new();
            for _ in 0..samples {
                let mut regions = Vec::new();
                for contig in 0..contigs {
                    let mut at_base = 1 + draw.next(10);
                    while at_base <= ground_end {
                        let width = match draw.next(10) {
                            0 => 1 + draw.next(150),
                            1 | 2 => 1 + draw.next(20),
                            _ => 1 + draw.next(4),
                        };
                        let end = at_base + width - 1;
                        regions.push(region_on(contig, at_base, end));
                        at_base = end + 1 + draw.next(6);
                    }
                }
                layouts.push(regions);
            }

            let width = 1 + draw.next(12);
            let mut builder_regions = Vec::new();
            for contig in 0..contigs {
                let mut at_base = 1u64;
                while at_base <= ground_end {
                    if draw.next(4) != 0 {
                        builder_regions.push(region_on(contig, at_base, at_base + width - 1));
                    }
                    at_base += width;
                }
            }

            // The oracle: every builder is handed the whole stretch.
            let whole: Vec<Vec<SampleLocusObservations>> = layouts
                .iter()
                .map(|regions| regions.iter().map(|at| observation_over(*at)).collect())
                .collect();
            let whole_slices: Vec<&[SampleLocusObservations]> =
                whole.iter().map(Vec::as_slice).collect();
            // **Rendered through `fixtures::render`, not `Debug` on the whole outcome**, which
            // is what leaves the window coverage out: only the cached side takes that
            // measurement, so comparing it here would fail on a difference that is not a
            // disagreement. What checks that field is
            // `serial::tests::the_cached_driver_measures_windows_where_the_in_memory_one_has_none`.
            let oracle: Vec<String> = builder_regions
                .iter()
                .map(|at| {
                    super::super::fixtures::render(&build_region(
                        *at,
                        &whole_slices,
                        max_span,
                        keep,
                    ))
                    .join(" | ")
                })
                .collect();

            // Under test: every builder is handed the cache's window.
            let sources: Vec<CountingSource> = layouts
                .iter()
                .map(|regions| source_over(regions).0)
                .collect();
            let mut cache = ObservationCache::over_fixture(sources);
            let mut through_cache = Vec::new();
            for at in &builder_regions {
                cache.evict_before(GenomePosition {
                    contig: at.contig,
                    position: at.start,
                });
                cache.cover(*at).expect("the fixture sources hold");
                through_cache.push(cache.with_observations(*at, |windows| {
                    super::super::fixtures::render(
                        &build_region_windowed::<std::convert::Infallible>(
                            *at,
                            windows,
                            max_span,
                            keep,
                            &|_, _| unreachable!("the fixture's sources build every record"),
                        )
                        .expect("an infallible build"),
                    )
                    .join(" | ")
                }));
            }

            if oracle != through_cache {
                let first = oracle
                    .iter()
                    .zip(&through_cache)
                    .position(|(whole, windowed)| whole != windowed)
                    .expect("they differ somewhere");
                disagreements.push(format!(
                    "seed {seed}, region {}\n  whole stretch: {}\n  through cache: {}",
                    builder_regions[first], oracle[first], through_cache[first],
                ));
            }
        }

        assert!(
            disagreements.is_empty(),
            "{} of 200 layouts disagreed with the whole-stretch builder; the first:\n{}",
            disagreements.len(),
            disagreements.first().expect("not empty"),
        );
    }
}
