//! The organiser — it takes the builders' outcomes and releases their loci in genome order.
//!
//! Builders own short **building regions** and finish out of order; the consumer must see
//! genome order (`doc/devel/ng/spec/cohort_merge.md` §6.3). This file holds what puts them back
//! in it: outcomes keyed by the index the run handed the region out under, released along an
//! unbroken run of indexes, with the loci that lost their ground to an earlier locus dropped on
//! the way out (§6.1).
//!
//! **The observation cache is [`super::observation_cache`], not here.** The two lived in one
//! file until E3's review priced the argument for it: the case was that the organiser would
//! become the cache's only writer and `cover` and `evict_before` could then turn private, but
//! `super::serial`'s cached driver calls both from a sibling module and nothing in the plan
//! removes it. The reachable narrowing is `pub(super)`, which a file of its own gets equally —
//! and does.
//!
//! **What the organiser does not yet hold is that cache**, which arch §4 gives it. Drawing the
//! readers forward is the driver's today (`super::parallel::merge_cohort_in_parallel`), because
//! the cache is generic over its source's error type and the run's own `ObservationSource` and
//! `RunError` do not exist yet.

use std::collections::{BTreeMap, VecDeque};

use super::build::{CohortObservation, RegionOutcome};
use crate::ng::types::{GenomePosition, GenomeRegion};

/// Where one building region falls in the order the run hands them out — 0 for the run's
/// first, 1 for its second, counting on across every analysed region of the run.
///
/// **It is the run's numbering, not a coordinate**, and that is the point: the organiser
/// releases on a gapless run of indexes, so it can tell "region 7 has not arrived yet" from
/// "region 7 found nothing", which two genome positions cannot say. The regions themselves
/// are [`building_regions_of`]'s, taken in the order that iterator yields them.
#[derive(Copy, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct RegionIndex(pub u64);

impl std::fmt::Display for RegionIndex {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{}", self.0)
    }
}

/// A run that ended owing output — the release-level guard production keeps as
/// `WriterError::MissingChunks` (`var_calling/vcf_writer.rs:152-158`), and what arch §5 will
/// fold into `RunError` when that type exists.
///
/// **The whole of E1 is that a gap must be an error and not a truncation** (spec §6.3): the
/// output would simply stop early, and the failed-locus total would be short by everything the
/// lost regions refused — the one number that says the width bound is charging more than
/// expected (spec §3.3). Both are silent, and a partly-lost run answers exactly like a complete
/// one.
///
/// **An enum, because the two ways to end short can happen together and can happen apart, and
/// no combination of them means "nothing was lost".** Arch §5 sketches one struct with one
/// count, which production can afford because it emits from inside its own submit and has no
/// second step to forget; here taking the released loci is the caller's own call, so a run can
/// also end with loci released and never taken.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum RunEndedShort {
    /// Regions the run handed out whose loci never reached the caller: the first index that
    /// stalled the drain, and how many regions from it onwards were never released. **The
    /// count covers both halves of a stall** — the region that never delivered, and every
    /// region behind it that delivered and could not be let out.
    #[error(
        "region {first_stalled} never delivered its result, so {regions} region(s) were \
         handed out and never released"
    )]
    RegionsNeverReleased {
        first_stalled: RegionIndex,
        regions: u64,
    },

    /// Loci released in order and never taken from the organiser. Nothing stalled; the caller
    /// stopped draining.
    #[error("{loci} released locus/loci were never taken from the organiser")]
    LociNeverDrained { loci: u64 },

    /// Both at once, which one count could not express.
    #[error(
        "region {first_stalled} never delivered its result, so {regions} region(s) were \
         handed out and never released, and {loci} released locus/loci were never taken"
    )]
    RegionsNeverReleasedAndLociNeverDrained {
        first_stalled: RegionIndex,
        regions: u64,
        loci: u64,
    },
}

/// What a finished run has to say about the loci it did not emit.
///
/// **Both counts leave with the organiser**, because `finish` consumes it: a caller that read
/// them afterwards could not, and one that had to read them first would lose them by taking the
/// obvious order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MergeTally {
    /// Loci refused by the width bound — spec §3.3's number, the only signal that the bound is
    /// charging more than expected.
    pub failed_loci: u64,
    /// Loci dropped because an earlier locus already owned their first base (spec §6.1).
    /// **Expected to be zero**; [`Organiser`] says why.
    pub displaced_loci: u64,
}

/// Holds the builders' outcomes and releases their loci in genome order (spec §6.1, §6.3).
///
/// **Builders finish out of order and the consumer must see genome order**, so the outcomes
/// are keyed by region index and drained while the head equals the next expected index — the
/// reorder map production's VCF writer drains on `next_expected`, carried whole
/// (`var_calling/vcf_writer.rs:168-176`).
///
/// **Every region delivers exactly one outcome, including the empty ones** (spec §6.3). The
/// drain advances only along an unbroken run of indexes, so a region that found nothing must
/// still submit — an empty [`RegionOutcome`] — or every region behind it is held for ever,
/// and what those regions found is lost with nothing to say so. That is what
/// [`RunEndedShort`] refuses at the end of the run.
///
/// **Waiting for the predecessor is not merely about order.** A region's loci can only be
/// confirmed once the region before it has been resolved, because that is what says whether a
/// locus owned earlier already covers this region's ground (spec §6.3). So the release point
/// is also the resolution point: each region's loci and failed spans are taken in genome
/// order, and one whose first base falls on ground an earlier locus already owns is dropped —
/// the earlier owner stands, whether it was emitted or failed
/// ([`resolve_and_release`](Self::resolve_and_release)).
///
/// **That rule is a safety net, not a live one, and saying which matters.** Under
/// [`build_region`](super::build::build_region)'s input contract two loci owned by different
/// regions cannot overlap at all, so nothing in a healthy run reaches it. Three terms carry
/// the argument and are worth naming first: a locus's **members** are the per-sample
/// observations it was closed over; a builder's **window** is the observations it is handed
/// for its own region; and its members are **chained** when each begins at or before the
/// furthest reach of those already open, which is what keeps one locus open across them
/// (spec §4.1, [`super::close::LocusCloser`]). The argument in full, because the code is
/// otherwise a branch nobody can explain:
///
/// - a builder is handed every observation overlapping its own ground, including those that
///   opened earlier. **This is a discipline of whoever draws the cache forward, not a property
///   of the cache**: [`ObservationCache::evict_before`] drops what ends before whatever
///   position it is handed, and the driver chooses that position to be the building region's
///   own first base (`super::serial::merge_cohort_through_cache`), so everything reaching into
///   the region survives. Hand it the region's *last* base instead and the argument fails at
///   once — measured, and the regression test is
///   `super::serial::tests::refuse_displaced_loci`. With several builders in flight the choice
///   is `super::parallel::merge_cohort_in_parallel`'s, and it evicts at the first base of the
///   round's first region — the earliest that is live — which keeps the argument as it stands;
/// - so if a locus L owned by an earlier region reaches into this one, every member of L that
///   reaches this region's first base is in this builder's window. Take any member of L
///   starting inside this region: it reaches at least its own start, so it is in the window,
///   and so is the member it overlaps, and so on backwards — the sub-chain is unbroken and it
///   begins before this region;
/// - so this builder closes that ground as one locus starting before its own first base, and
///   skips it as an earlier region's — the `locus.region.start < builder_region.start` arm of
///   [`build_region`](super::build::build_region). Every locus it does own begins after that
///   chain ends, which is after L's own reach.
///
/// The rule is kept because it is spec §6.1's, because the contract it rests on belongs to
/// whoever feeds the builders rather than to this file, and because it costs one comparison
/// per locus. **What it does not have is a fixture that reaches it through a real merge.** No
/// driver builds an organiser yet, so every test below hands it fabricated outcomes, and
/// [`displaced_locus_count`](Self::displaced_locus_count) is how a run would say the argument
/// had failed. The same claim is made from the other end in `super::serial`, whose
/// `tests::refuse_overlapping_ground` asserts that the loci and failed spans the **cached**
/// driver produces never share ground, and whose `tests::refuse_displaced_loci` drives the
/// builders' own outcomes into a real organiser and asserts that none is displaced. Both run
/// over six of that file's twenty-eight tests — among them the two hundred random layouts —
/// and over the 305–330 deletion that a building region's boundary falls inside, which is the
/// shape the plan named. What the first asserts is disjointness alone: it sorts the spans
/// together before comparing them, so the order they came out in is not checked there.
///
/// **It does not yet hold the observation cache**, which arch §4 gives it. The cache is drawn
/// forward and evicted by whoever hands regions out, and that is the parallel arrangement's
/// shape to settle (the plan's E3); until then the two live side by side in this file rather
/// than one inside the other.
#[derive(Debug)]
pub struct Organiser {
    /// The next region index that may be released. Everything below it has been released
    /// already; nothing at or above it has.
    next_expected_region: RegionIndex,
    /// Outcomes that arrived before their turn, keyed by region index.
    held_outcomes: BTreeMap<RegionIndex, RegionOutcome>,
    /// Loci already released, in genome order, waiting for the caller to take them. They are
    /// held rather than handed straight on because releasing and taking are two different
    /// moments: a region becomes releasable when its predecessor arrives, which is not when
    /// the caller next asks.
    released_loci: VecDeque<CohortObservation>,
    /// Summed over the regions **released so far**, not over those submitted — the number
    /// spec §3.3 requires to reach the run summary. Counting at release rather than at
    /// submission is what lets a failed locus an earlier one displaced be dropped without the
    /// total having counted it already.
    failed_locus_count: u64,
    /// The last base owned by a locus already resolved, or `None` before the first one — the
    /// **frontier**, which is the word the comments and tests below use for it. A locus whose
    /// first position is at or before it belongs to ground an earlier locus already owns, and
    /// is dropped (spec §6.1).
    owned_through: Option<GenomePosition>,
    /// How often that rule fired. **Nothing in a healthy run should reach it** — see
    /// [`Organiser`]'s note on why the rule is a safety net — so a run that reports a non-zero
    /// figure here has had the exclusion argument broken somewhere upstream, and that is worth
    /// being able to see rather than infer.
    displaced_locus_count: u64,
}

impl Default for Organiser {
    /// [`Organiser::new`], and only that: a second construction path would be a second place
    /// to forget a field, which is the thing `new` is written out longhand to prevent.
    fn default() -> Self {
        Self::new()
    }
}

impl Organiser {
    /// An organiser expecting region index 0 and holding nothing.
    ///
    /// Written out field by field rather than deferring to `Default`, so that a field E2 or E3
    /// adds has to be given a starting value here instead of taking whatever its own `Default`
    /// happens to be.
    #[must_use]
    pub fn new() -> Self {
        Self {
            next_expected_region: RegionIndex(0),
            held_outcomes: BTreeMap::new(),
            released_loci: VecDeque::new(),
            failed_locus_count: 0,
            owned_through: None,
            displaced_locus_count: 0,
        }
    }

    /// Take one region's outcome, and release everything the arrival makes releasable.
    ///
    /// **Releasing here rather than in [`drain_ready`](Self::drain_ready) is production's
    /// shape** (`var_calling/vcf_writer.rs:246`): a region becomes releasable the moment the
    /// index in front of it arrives, which has nothing to do with when the caller next asks
    /// for loci. Keeping the two apart is also what makes
    /// [`failed_locus_count`](Self::failed_locus_count) meaningful before anything has been
    /// drained.
    ///
    /// **Panics** when `index` was submitted already or has been released. Both are bugs in
    /// whoever hands the regions out rather than facts about the data, and both are caught
    /// mid-flight, where the release order is already wrong and nothing coherent can follow —
    /// which is what separates them from [`RunEndedShort`], a bug too, but one caught at the
    /// end of a run where what was lost can still be named and reported.
    pub fn submit(&mut self, index: RegionIndex, outcome: RegionOutcome) {
        assert!(
            index >= self.next_expected_region,
            "the outcome of region {index} arrived after it was released (next expected {}), \
             so its loci would be released a second time, out of order",
            self.next_expected_region,
        );
        // **Checked before the insert, not inside it.** Writing the map inside the assertion's
        // condition would put the module's only insertion in an expression a later edit to
        // `debug_assert!` would stop evaluating — and this crate's release profile leaves debug
        // assertions off, so every region's outcome would be dropped in the shipped binary and
        // in no test. Checking first also keeps the first outcome, rather than replacing it
        // and then panicking.
        assert!(
            !self.held_outcomes.contains_key(&index),
            "region {index} delivered a second outcome, and a region owns its loci exactly once",
        );
        self.held_outcomes.insert(index, outcome);
        self.release_regions_in_turn();
    }

    /// Everything released and not yet taken, in genome order.
    ///
    /// **What is left when the iterator is dropped early stays here**, to be handed out by the
    /// next call: this takes one locus at a time from the front rather than emptying the buffer
    /// into the iterator, so a caller that stops halfway loses nothing. Loci still here at the
    /// end of a run are what [`finish`](Self::finish) refuses.
    #[must_use = "the loci stay in the organiser until the iterator is consumed"]
    pub fn drain_ready(&mut self) -> impl Iterator<Item = CohortObservation> + '_ {
        let released_loci = &mut self.released_loci;
        std::iter::from_fn(move || released_loci.pop_front())
    }

    /// The failed loci of every region released so far, summed (spec §3.3).
    #[must_use]
    pub fn failed_locus_count(&self) -> u64 {
        self.failed_locus_count
    }

    /// How many loci — emitted and failed alike — were dropped because an earlier locus
    /// already owned the ground they started on (spec §6.1).
    ///
    /// **Expected to be zero, and that is the point.** [`Organiser`]'s note gives the argument
    /// for why no builder working under `build_region`'s input contract can produce an
    /// overlapping pair; this counter is what would show the argument failing on real data
    /// instead of leaving it to be inferred from a locus quietly missing.
    #[must_use]
    pub fn displaced_locus_count(&self) -> u64 {
        self.displaced_locus_count
    }

    /// Nothing outstanding: every region the run handed out has been released, and every
    /// released locus taken (arch §4).
    ///
    /// **`regions_handed_out` is how many building regions the run dealt out**, their indexes
    /// being `0..regions_handed_out`. Without it the organiser cannot see a gap at the *tail*
    /// of a run — indexes that never submitted, with no later index behind them to hold — and
    /// a run that lost its last regions would look exactly like one that finished.
    #[must_use]
    pub fn is_finished(&self, regions_handed_out: u64) -> bool {
        self.next_expected_region.0 >= regions_handed_out
            && self.held_outcomes.is_empty()
            && self.released_loci.is_empty()
    }

    /// End the run: what it did not emit, or a refusal naming what would have been lost.
    ///
    /// Consuming, because there is nothing to ask an organiser afterwards — production's
    /// `VcfWriter::finish` has the same shape and the same reason
    /// (`var_calling/vcf_writer.rs:256`). `regions_handed_out` is
    /// [`is_finished`](Self::is_finished)'s, and this returns `Ok` on exactly the runs that
    /// method calls finished.
    ///
    /// **Panics** when an outcome was submitted for an index the run says it never handed out,
    /// which is the same class of hand-out bug [`submit`](Self::submit) refuses.
    pub fn finish(self, regions_handed_out: u64) -> Result<MergeTally, RunEndedShort> {
        // Destructured rather than field-accessed, so that anything the organiser gains at a
        // later step has to be answered for here — drained by the end of the run, or
        // deliberately not — instead of being left behind in silence.
        let Self {
            next_expected_region,
            held_outcomes,
            released_loci,
            failed_locus_count,
            owned_through: _,
            displaced_locus_count,
        } = self;

        if let Some((last_held, _)) = held_outcomes.last_key_value() {
            assert!(
                last_held.0 < regions_handed_out,
                "region {last_held} delivered an outcome though the run handed out only \
                 {regions_handed_out} region(s)",
            );
        }

        // Every index from the cursor to the last one handed out: the region that stalled the
        // drain, and every region behind it, whether it delivered or not. `saturating_sub`
        // because a caller that under-reports its own hand-out is the panic above, not a wrap.
        let regions = regions_handed_out.saturating_sub(next_expected_region.0);
        let loci = released_loci.len() as u64;

        match (regions, loci) {
            (0, 0) => Ok(MergeTally {
                failed_loci: failed_locus_count,
                displaced_loci: displaced_locus_count,
            }),
            (0, loci) => Err(RunEndedShort::LociNeverDrained { loci }),
            (regions, 0) => Err(RunEndedShort::RegionsNeverReleased {
                first_stalled: next_expected_region,
                regions,
            }),
            (regions, loci) => Err(RunEndedShort::RegionsNeverReleasedAndLociNeverDrained {
                first_stalled: next_expected_region,
                regions,
                loci,
            }),
        }
    }

    /// Release the unbroken run of regions now at the head, oldest first.
    fn release_regions_in_turn(&mut self) {
        while let Some(outcome) = self.held_outcomes.remove(&self.next_expected_region) {
            // Destructured for the reason `build.rs` gives where it consumes a comparable
            // type: anything `RegionOutcome` gains has to be answered for here — carried or
            // dropped deliberately — rather than vanishing.
            let RegionOutcome {
                cohort_observations,
                failed_locus_spans,
            } = outcome;
            self.resolve_and_release(cohort_observations, failed_locus_spans);
            // PANIC-FREE: the cursor rises once per released region, and a run cannot hand out
            // more than `u64::MAX` building regions.
            self.next_expected_region = RegionIndex(
                self.next_expected_region
                    .0
                    .checked_add(1)
                    .expect("a run cannot hand out more than u64::MAX building regions"),
            );
        }
    }

    /// One region's loci and failed spans, taken in genome order, each claiming its ground or
    /// losing it to a locus that got there first (spec §6.1).
    ///
    /// **The two vectors are merged rather than handled apart**, because the rule is about
    /// ground and a failed locus owns ground exactly as an emitted one does: it wins against
    /// what overlaps it, and the only things that differ are at the end of the line — nothing
    /// is emitted for it, and the run counts it (spec §3.2). A separate pass over each would
    /// let a failed locus be resolved against a frontier that had already run past it.
    ///
    /// **What the two vectors have to be, and who owes it.** Interleaved by first base they
    /// must form one ascending sequence of spans that do not overlap. That is what one walk
    /// over one region produces — it closes disjoint loci in genome order and sends each to
    /// whichever vector its verdict picks ([`RegionOutcome`]) — but *neither vector says it
    /// alone*, and nothing here checks it: it is the submitter's to keep, and a submitter that
    /// breaks it sees the rule fire *inside* one region. Given it, only the join between two
    /// regions can reach the rule at all, and [`Organiser`]'s note says why even that does not,
    /// today.
    fn resolve_and_release(&mut self, loci: Vec<CohortObservation>, failed: Vec<GenomeRegion>) {
        // Driven off `next_if` rather than a peek-then-take pair, so that no branch has to
        // assert that the value it just looked at is still there.
        let mut loci = loci.into_iter().peekable();
        for span in failed {
            let span_first = first_base_of(span);
            while let Some(locus) = loci.next_if(|locus| first_base_of(locus.region) <= span_first)
            {
                self.claim_and_release(locus);
            }
            if self.claim(span) {
                // Saturating where the region cursor is checked, and the difference is
                // deliberate: a saturated total is a truer answer than a wrap for a count,
                // while a cursor that stopped advancing would release one region for ever.
                self.failed_locus_count = self.failed_locus_count.saturating_add(1);
            }
        }
        for locus in loci {
            self.claim_and_release(locus);
        }
    }

    /// Release `locus` if it claims its ground, and drop it if an earlier locus owns it.
    fn claim_and_release(&mut self, locus: CohortObservation) {
        if self.claim(locus.region) {
            self.released_loci.push_back(locus);
        }
    }

    /// Claim `locus_span`'s ground for the locus that covers it, or refuse it to a locus that
    /// got there first — spec §6.1's one rule, with no special case for a failed locus. The
    /// span is a **locus's** own ground, not a building region's, which is the other thing this
    /// file calls a region.
    ///
    /// **The test is the first base, not the overlap.** A locus belongs to the builder whose
    /// region holds its first position, so the earlier *start* is the earlier owner; asking
    /// instead whether the two spans intersect would give the same answer here and a different
    /// one on any pair the walk could not have produced.
    ///
    /// **"Earlier" means claimed first, not lower-numbered.** Under `build_region`'s input
    /// contract the two are the same, since a builder starts a locus only inside its own
    /// region and the regions are resolved in order. They come apart only where that contract
    /// has already failed — a later region delivering a locus that starts before the standing
    /// owner's — and there the standing owner keeps its ground.
    #[must_use = "a locus that did not claim its ground is dropped, and the caller decides how"]
    fn claim(&mut self, locus_span: GenomeRegion) -> bool {
        if self
            .owned_through
            .is_some_and(|owned| first_base_of(locus_span) <= owned)
        {
            self.displaced_locus_count = self.displaced_locus_count.saturating_add(1);
            return false;
        }
        // **Assigned, not maxed, and the frontier cannot retreat**: a span reaches here only if
        // its first base is past the frontier, and `last_base_of` never returns anything before
        // `first_base_of`, so its last base is past the frontier too. A `max` here would be a
        // branch no input can take, which reads as a hazard someone guarded against and leaves
        // the next reader looking for the case.
        self.owned_through = Some(last_base_of(locus_span));
        true
    }
}

/// A locus span's first base, genome-wide — what decides which of two loci is the earlier
/// owner.
///
/// **The ends are ordered rather than trusted**, as [`ObservationCache::with_observations`] and
/// [`ObservationCache::cover`] order theirs and for the same reason: `GenomeRegion` has public
/// fields and no constructor putting them in order, and an inverted span would otherwise walk
/// the frontier backwards — past its own first base, releasing the locus behind it.
fn first_base_of(locus_span: GenomeRegion) -> GenomePosition {
    GenomePosition {
        contig: locus_span.contig,
        position: locus_span.start.min(locus_span.end),
    }
}

/// A locus span's last base, genome-wide — where the frontier stands once the span is claimed.
/// Ordered for [`first_base_of`]'s reason, and never before it.
fn last_base_of(locus_span: GenomeRegion) -> GenomePosition {
    GenomePosition {
        contig: locus_span.contig,
        position: locus_span.end.max(locus_span.start),
    }
}

#[cfg(test)]
mod tests {
    use super::super::fixtures::{region, region_on};
    use super::*;

    /// A cohort locus over `region`. The organiser reads nothing but where a locus sits, so
    /// the allele table and the per-sample support are left empty: what these tests check is
    /// which loci come out and in what order, and a fabricated one carries that in its region.
    fn locus_at(region: GenomeRegion) -> CohortObservation {
        CohortObservation {
            region,
            alleles: Vec::new(),
            per_sample: Vec::new(),
        }
    }

    /// One region's outcome: loci over `locus_regions`, and failures over `failed_regions`.
    ///
    /// **Both lists must be ascending, and the two must not overlap each other**, because that
    /// is what one walk over one region produces: every locus it closes is disjoint from the
    /// last, and each one is either built or failed. A fixture that breaks it makes the overlap
    /// rule fire *inside* a single region, which no builder can do — five of the fixtures
    /// written before the rule existed did exactly that, and writing the rule is what found
    /// them.
    fn outcome_of(
        locus_regions: &[GenomeRegion],
        failed_regions: &[GenomeRegion],
    ) -> RegionOutcome {
        RegionOutcome {
            cohort_observations: locus_regions.iter().copied().map(locus_at).collect(),
            failed_locus_spans: failed_regions.to_vec(),
        }
    }

    /// The regions of everything the organiser will now part with, in the order it parts
    /// with them.
    fn drained_regions(organiser: &mut Organiser) -> Vec<GenomeRegion> {
        organiser
            .drain_ready()
            .map(|observation| observation.region)
            .collect()
    }

    /// The plain case: regions arriving in their own order release as they arrive.
    #[test]
    fn a_region_arriving_in_its_turn_releases_at_once() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert_eq!(drained_regions(&mut organiser), vec![region(1, 3)]);

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        assert_eq!(drained_regions(&mut organiser), vec![region(21, 24)]);
    }

    /// The reorder buffer's whole job: a region cannot pass the one in front of it.
    #[test]
    fn a_region_that_arrives_early_waits_for_the_one_before_it() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        assert!(
            drained_regions(&mut organiser).is_empty(),
            "region 1 was released before region 0 had arrived, so its loci were never \
             offered the chance to be displaced by an earlier owner's",
        );

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(1, 3), region(21, 24)],
        );
    }

    /// A whole convoy held behind one missing index comes out at once, in index order, when
    /// that index lands — and the order is the regions' own, not the order they arrived in.
    #[test]
    fn a_convoy_held_behind_a_gap_releases_in_region_order_when_the_gap_closes() {
        let mut organiser = Organiser::new();

        for (index, locus) in [
            (3, region(61, 63)),
            (1, region(21, 24)),
            (2, region(41, 41)),
        ] {
            organiser.submit(RegionIndex(index), outcome_of(&[locus], &[]));
        }
        assert!(drained_regions(&mut organiser).is_empty());

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(1, 3), region(21, 24), region(41, 41), region(61, 63)],
        );
    }

    /// **Exactly one outcome per region, empty ones included** (spec §6.3). The empty region
    /// carries no loci and no failures, and delivering it is the only thing that lets the
    /// region behind it out.
    #[test]
    fn an_empty_region_still_lets_the_region_behind_it_out() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        organiser.submit(RegionIndex(0), outcome_of(&[], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(21, 24)]);
    }

    /// The other half of the rule above: a region that never submits holds back every region
    /// after it, for ever. This is what a builder dropping its result costs.
    ///
    /// **The three held regions carry one, two and three loci**, so the refusal's count of
    /// *regions* (3) differs from the loci they hold (6): a count that added up the held loci
    /// instead would read 6 and be caught here. The reviewer found that mutation alive against
    /// an earlier fixture of three one-locus regions, where the two rules give the same number.
    #[test]
    fn a_region_that_never_submits_holds_back_every_region_after_it() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        organiser.submit(
            RegionIndex(2),
            outcome_of(&[region(31, 34), region(36, 36)], &[]),
        );
        organiser.submit(
            RegionIndex(3),
            outcome_of(&[region(41, 44), region(46, 46), region(48, 48)], &[]),
        );

        assert!(drained_regions(&mut organiser).is_empty());
        assert!(!organiser.is_finished(4));
        assert_eq!(
            organiser.finish(4),
            Err(RunEndedShort::RegionsNeverReleased {
                first_stalled: RegionIndex(0),
                regions: 4,
            }),
        );
    }

    /// **A gap at the tail of the run is a gap too**, and it is the one the organiser cannot
    /// see from its own state: nothing is held, because the regions that went missing have no
    /// later index behind them. Five regions dealt out, two heard from.
    #[test]
    fn finish_refuses_a_run_whose_last_regions_never_submitted() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3)], &[region(10, 15)]),
        );
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        drained_regions(&mut organiser);

        assert!(
            !organiser.is_finished(5),
            "three of the five regions handed out never delivered an outcome",
        );
        assert_eq!(
            organiser.finish(5),
            Err(RunEndedShort::RegionsNeverReleased {
                first_stalled: RegionIndex(2),
                regions: 3,
            }),
        );
    }

    /// Within one region the builder's order is the genome's, and the organiser keeps it.
    #[test]
    fn the_loci_of_one_region_come_out_in_the_order_the_builder_gave_them() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7), region(11, 18)], &[]),
        );

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(1, 3), region(7, 7), region(11, 18)],
        );
    }

    /// The count spec §3.3 needs is the sum over the regions, and it reaches it whatever
    /// order the regions arrive in.
    #[test]
    fn the_failed_locus_count_sums_every_released_region() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[], &[region(45, 60)]));
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[], &[region(1, 20), region(25, 40)]),
        );
        organiser.submit(RegionIndex(2), outcome_of(&[], &[]));

        assert_eq!(organiser.failed_locus_count(), 3);
    }

    /// **The count is of what has been released, not of what has been submitted**, so that
    /// the next step can drop a failed locus an earlier one displaced without the total
    /// having counted it already. Region 1 is in the organiser's hands and uncounted until
    /// region 0 lets it out.
    #[test]
    fn the_failed_locus_count_ignores_a_region_still_held_behind_a_gap() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(1),
            outcome_of(&[], &[region(40, 60), region(65, 90)]),
        );
        assert_eq!(organiser.failed_locus_count(), 0);

        organiser.submit(RegionIndex(0), outcome_of(&[], &[region(1, 30)]));
        assert_eq!(organiser.failed_locus_count(), 3);
    }

    /// **Releasing happens when the region arrives, not when the caller asks.** Nothing has
    /// been drained here, and the failed loci of both regions are already counted — which is
    /// what says the release ran inside `submit`.
    #[test]
    fn a_region_is_released_on_arrival_rather_than_at_the_next_drain() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3)], &[region(10, 60)]),
        );
        organiser.submit(RegionIndex(1), outcome_of(&[], &[region(80, 120)]));

        assert_eq!(organiser.failed_locus_count(), 2);
    }

    /// A drain the caller stops halfway leaves the rest where it was: the loci are taken one
    /// at a time from the front, not emptied into the iterator and dropped with it.
    #[test]
    fn a_drain_stopped_halfway_leaves_the_rest_for_the_next_call() {
        let mut organiser = Organiser::new();
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7), region(11, 18)], &[]),
        );

        let first = organiser
            .drain_ready()
            .next()
            .expect("the region's first locus");
        assert_eq!(first.region, region(1, 3));

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(7, 7), region(11, 18)],
        );
    }

    /// Nothing outstanding, which is what a run asserts at its end — and `finish` returns `Ok`
    /// on exactly those runs.
    #[test]
    fn an_organiser_is_finished_once_every_region_is_released_and_drained() {
        let mut organiser = Organiser::new();
        assert!(Organiser::new().is_finished(0));

        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        assert!(
            !organiser.is_finished(1),
            "a locus released and not yet taken is a locus the run still owes its output",
        );

        drained_regions(&mut organiser);
        assert!(organiser.is_finished(1));
        assert_eq!(
            organiser.finish(1),
            Ok(MergeTally {
                failed_loci: 0,
                displaced_loci: 0,
            }),
        );
    }

    /// Loci released and never taken truncate the output exactly as a missing region does,
    /// so `finish` refuses them too — and says so without inventing a gap.
    #[test]
    fn finish_refuses_a_run_with_loci_released_and_never_drained() {
        let mut organiser = Organiser::new();
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7)], &[]),
        );

        assert_eq!(
            organiser.finish(1),
            Err(RunEndedShort::LociNeverDrained { loci: 2 }),
        );
    }

    /// **Both ways of ending short at once**, each named by its own count: region 1 stalled
    /// the drain and region 2 is held behind it — two regions handed out and never released,
    /// though they hold three loci between them — while region 0's two loci were released and
    /// never taken. A refusal that could carry one count at a time would drop half of this.
    #[test]
    fn finish_names_both_counts_when_a_stall_and_undrained_loci_coincide() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 7)], &[]),
        );
        organiser.submit(
            RegionIndex(2),
            outcome_of(&[region(41, 41), region(45, 45), region(51, 51)], &[]),
        );

        assert_eq!(
            organiser.finish(3),
            Err(RunEndedShort::RegionsNeverReleasedAndLociNeverDrained {
                first_stalled: RegionIndex(1),
                regions: 2,
                loci: 2,
            }),
        );
    }

    /// The refusal reads back the way a person will read it, each count against its own noun.
    /// The two counts are deliberately different, so swapping them in the message cannot pass.
    #[test]
    fn the_refusal_names_each_count_against_its_own_noun() {
        let both = RunEndedShort::RegionsNeverReleasedAndLociNeverDrained {
            first_stalled: RegionIndex(4),
            regions: 3,
            loci: 7,
        };
        assert_eq!(
            both.to_string(),
            "region 4 never delivered its result, so 3 region(s) were handed out and never \
             released, and 7 released locus/loci were never taken",
        );

        assert_eq!(
            RunEndedShort::LociNeverDrained { loci: 7 }.to_string(),
            "7 released locus/loci were never taken from the organiser",
            "a run that lost nothing to a stall must not have one described to it",
        );
    }

    /// The end of a healthy run: the failed-locus total comes back, and nothing is owed.
    #[test]
    fn finish_returns_the_failed_locus_total_of_a_run_that_ended_clean() {
        let mut organiser = Organiser::new();
        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3)], &[region(10, 60)]),
        );
        organiser.submit(RegionIndex(1), outcome_of(&[], &[]));
        drained_regions(&mut organiser);

        assert_eq!(
            organiser.finish(2),
            Ok(MergeTally {
                failed_loci: 1,
                displaced_loci: 0,
            }),
        );
    }

    /// A run over no regions at all is finished and owes nothing.
    #[test]
    fn an_organiser_that_was_given_nothing_finishes_clean() {
        assert_eq!(
            Organiser::new().finish(0),
            Ok(MergeTally {
                failed_loci: 0,
                displaced_loci: 0,
            }),
        );
    }

    /// One outcome per region, whether or not its loci have gone out yet.
    #[test]
    #[should_panic(expected = "delivered a second outcome")]
    fn a_region_submitted_twice_is_refused() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));
    }

    /// The duplicate is refused **and the first outcome is the one kept** — the check runs
    /// before the map is written, so a second delivery cannot displace what the region
    /// already owns on its way to the panic.
    #[test]
    fn a_refused_second_outcome_does_not_displace_the_first() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(1), outcome_of(&[region(21, 24)], &[]));

        let refused = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            organiser.submit(RegionIndex(1), outcome_of(&[region(91, 94)], &[]));
        }));
        assert!(
            refused.is_err(),
            "a second outcome for region 1 must be refused"
        );

        organiser.submit(RegionIndex(0), outcome_of(&[], &[]));
        assert_eq!(drained_regions(&mut organiser), vec![region(21, 24)]);
    }

    /// An index that has already been released cannot be submitted again either — the drain
    /// has passed it, so its loci would come out behind loci that are already downstream.
    #[test]
    #[should_panic(expected = "arrived after it was released")]
    fn a_region_submitted_after_it_was_released_is_refused() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
        organiser.submit(RegionIndex(0), outcome_of(&[region(1, 3)], &[]));
    }

    // ---------------------------------------------------------------
    // Overlap resolution. **Every fixture below is fabricated**, because no builder working
    // under `build_region`'s input contract can produce an overlapping pair — the argument is
    // on [`Organiser`]. These state what the rule does if it ever is reached; the fixtures
    // `super::super::serial::tests::refuse_overlapping_ground` sees state, from the other end,
    // that it is not needed.
    // ---------------------------------------------------------------

    /// The rule itself: of two loci that overlap, the one that starts earlier stands.
    #[test]
    fn a_locus_starting_on_ground_an_earlier_locus_owns_is_dropped() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 40)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(25, 60)], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(10, 40)]);
        assert_eq!(organiser.displaced_locus_count(), 1);
        assert_eq!(
            organiser.finish(2),
            Ok(MergeTally {
                failed_loci: 0,
                displaced_loci: 1,
            }),
            "the tally is what a run reads, and it must carry the displacement too",
        );
    }

    /// **The test is the first base, not the overlap.** A locus that opens on the very last
    /// base an earlier one owns is dropped; one that opens on the next base is kept, even
    /// though the two are adjacent.
    #[test]
    fn the_boundary_is_the_earlier_locus_last_base() {
        let mut on_the_last_base = Organiser::new();
        on_the_last_base.submit(RegionIndex(0), outcome_of(&[region(10, 40)], &[]));
        on_the_last_base.submit(RegionIndex(1), outcome_of(&[region(40, 55)], &[]));
        assert_eq!(
            drained_regions(&mut on_the_last_base),
            vec![region(10, 40)],
            "a locus opening on base 40 starts inside ground owned through base 40",
        );

        let mut on_the_next_base = Organiser::new();
        on_the_next_base.submit(RegionIndex(0), outcome_of(&[region(10, 40)], &[]));
        on_the_next_base.submit(RegionIndex(1), outcome_of(&[region(41, 55)], &[]));
        assert_eq!(
            drained_regions(&mut on_the_next_base),
            vec![region(10, 40), region(41, 55)],
        );
        assert_eq!(on_the_next_base.displaced_locus_count(), 0);
    }

    /// **A failed locus displaces exactly as an emitted one does** (spec §3.2, §6.1). This is
    /// the case the rule exists for: a wide deletion refused by the width bound owns ground
    /// that a later builder, which never saw what opened there, would build from a partial
    /// picture.
    #[test]
    fn a_failed_locus_displaces_what_starts_inside_its_span() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[], &[region(10, 200)]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(120, 130)], &[]));
        organiser.submit(RegionIndex(2), outcome_of(&[region(210, 214)], &[]));

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(210, 214)],
            "the locus built inside the refused span survived with nothing to displace it",
        );
        assert_eq!(organiser.failed_locus_count(), 1);
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// And the same the other way round: a failed locus starting inside an emitted one's
    /// ground is dropped, and — because the count is taken at release — is never counted.
    #[test]
    fn a_failed_locus_displaced_by_an_earlier_one_is_not_counted() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 90)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[], &[region(50, 300)]));

        assert_eq!(drained_regions(&mut organiser), vec![region(10, 90)]);
        assert_eq!(
            organiser.failed_locus_count(),
            0,
            "a refusal the run never owned must not reach the run summary",
        );
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// **The two lists are resolved as one sequence, not one after the other.** Region 0's
    /// failure at 10-200 and its locus at 220-224 interleave with region 1's locus at 100-110;
    /// resolving all the loci first and then all the failures would let the failure be judged
    /// against a frontier that had already run to 224.
    #[test]
    fn loci_and_failed_spans_are_resolved_in_one_genome_order() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(220, 224)], &[region(10, 200)]),
        );
        organiser.submit(RegionIndex(1), outcome_of(&[region(100, 110)], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(220, 224)]);
        assert_eq!(organiser.failed_locus_count(), 1);
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// One wide locus can displace several later ones, and the frontier does not retreat to
    /// the last one kept.
    #[test]
    fn one_wide_locus_displaces_every_locus_that_opens_inside_it() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 300)], &[]));
        for (index, locus) in [
            (1, region(40, 44)),
            (2, region(120, 124)),
            (3, region(280, 284)),
        ] {
            organiser.submit(RegionIndex(index), outcome_of(&[locus], &[]));
        }
        organiser.submit(RegionIndex(4), outcome_of(&[region(301, 305)], &[]));

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(10, 300), region(301, 305)],
        );
        assert_eq!(organiser.displaced_locus_count(), 3);
    }

    /// **The frontier is the latest ground owned, not the first.** A locus kept after another
    /// owns what it covers just as the first did, so the one that opens inside *it* is
    /// displaced — which a frontier that stopped at the first locus would let through.
    #[test]
    fn a_locus_kept_after_another_owns_its_own_ground_too() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 20)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(30, 100)], &[]));
        organiser.submit(RegionIndex(2), outcome_of(&[region(50, 60)], &[]));

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(10, 20), region(30, 100)],
        );
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// **A displaced locus owns nothing** — the frontier stays on the standing owner's last
    /// base, not the losing locus's. Region 1's locus reaches to 100 and loses to region 0's,
    /// so the ground from 41 to 100 is nobody's and region 2's locus there stands. Every other
    /// fixture here hides the difference, because in each the displaced locus ends inside the
    /// standing owner's span.
    #[test]
    fn a_displaced_locus_does_not_own_the_ground_it_lost() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 40)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(25, 100)], &[]));
        organiser.submit(RegionIndex(2), outcome_of(&[region(50, 60)], &[]));

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(10, 40), region(50, 60)],
        );
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// **"Earlier" means claimed first, not lower-numbered.** A later region delivering a locus
    /// that starts before the standing owner's cannot happen under `build_region`'s input
    /// contract; where it does, the standing owner keeps its ground. Nothing else in the suite
    /// separates the two readings of spec §6.1's "earlier".
    #[test]
    fn a_later_region_delivering_an_earlier_start_still_loses() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(30, 40)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(10, 35)], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(30, 40)]);
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// **A displaced locus is not a failed one.** The tally keeps them apart: the failed total
    /// is what the width bound refused, and displacement — which here drops one locus and one
    /// refusal — adds to neither it nor the ground the run reports as refused.
    #[test]
    fn the_tally_keeps_the_failed_and_the_displaced_apart() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 300)], &[]));
        organiser.submit(
            RegionIndex(1),
            outcome_of(&[region(50, 60)], &[region(70, 80)]),
        );
        drained_regions(&mut organiser);

        assert_eq!(
            organiser.finish(2),
            Ok(MergeTally {
                failed_loci: 0,
                displaced_loci: 2,
            }),
        );
    }

    /// A span whose ends are the wrong way round is read as the ground it names, which is how
    /// [`ObservationCache::cover`] and [`ObservationCache::with_observations`] read one too.
    /// Taken raw, its last base would put the frontier *behind* its own first base, and the
    /// locus that follows would be released whether or not it overlapped.
    #[test]
    fn a_span_with_its_ends_the_wrong_way_round_still_owns_its_ground() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(50, 40)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region(45, 60)], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(50, 40)]);
        assert_eq!(
            organiser.displaced_locus_count(),
            1,
            "45 falls inside the ground 40-50 names, whichever way round it was written",
        );

        // And the other end: an inverted span *arriving* on owned ground must lose it. Read
        // raw, this one's first base would be 50 — past the frontier at 45 — and it would
        // claim ground 40-45 that region 0 already owns.
        let mut arriving_inverted = Organiser::new();
        arriving_inverted.submit(RegionIndex(0), outcome_of(&[region(10, 45)], &[]));
        arriving_inverted.submit(RegionIndex(1), outcome_of(&[region(50, 40)], &[]));

        assert_eq!(
            drained_regions(&mut arriving_inverted),
            vec![region(10, 45)],
        );
        assert_eq!(arriving_inverted.displaced_locus_count(), 1);
    }

    /// Ownership does not cross a contig: a locus at the start of the next contig is kept
    /// however far the last locus of the previous one reached.
    #[test]
    fn a_locus_on_the_next_contig_is_never_displaced() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region_on(0, 10, 4_000)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[region_on(1, 5, 9)], &[]));

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region_on(0, 10, 4_000), region_on(1, 5, 9)],
        );
        assert_eq!(organiser.displaced_locus_count(), 0);
    }

    /// **The frontier is the organiser's, not one region's**, so it carries across the join
    /// between two regions and across a region that had nothing in it.
    #[test]
    fn the_frontier_carries_across_an_empty_region() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 300)], &[]));
        organiser.submit(RegionIndex(1), outcome_of(&[], &[]));
        organiser.submit(RegionIndex(2), outcome_of(&[region(250, 260)], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(10, 300)]);
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// The resolution runs in **region order, not arrival order** — which is the reason the
    /// release waits for the predecessor at all (spec §6.3). Region 1 arrives first and would
    /// win if the frontier were built as outcomes landed; it loses because region 0 owns the
    /// ground when the two are resolved.
    #[test]
    fn resolution_follows_region_order_rather_than_arrival_order() {
        let mut organiser = Organiser::new();

        organiser.submit(RegionIndex(1), outcome_of(&[region(25, 60)], &[]));
        organiser.submit(RegionIndex(0), outcome_of(&[region(10, 40)], &[]));

        assert_eq!(drained_regions(&mut organiser), vec![region(10, 40)]);
        assert_eq!(organiser.displaced_locus_count(), 1);
    }

    /// Nothing is displaced in a run where nothing overlaps, which is every run the builders
    /// can actually produce.
    #[test]
    fn nothing_is_displaced_when_no_two_loci_overlap() {
        let mut organiser = Organiser::new();

        organiser.submit(
            RegionIndex(0),
            outcome_of(&[region(1, 3), region(7, 9)], &[region(20, 40)]),
        );
        organiser.submit(
            RegionIndex(1),
            outcome_of(&[region(41, 44)], &[region(50, 70)]),
        );

        assert_eq!(
            drained_regions(&mut organiser),
            vec![region(1, 3), region(7, 9), region(41, 44)],
        );
        assert_eq!(organiser.failed_locus_count(), 2);
        assert_eq!(organiser.displaced_locus_count(), 0);
    }

    /// An outcome for a region the run says it never handed out is the same class of hand-out
    /// bug as a duplicate, and is refused where it is first visible.
    #[test]
    #[should_panic(expected = "though the run handed out only")]
    fn an_outcome_for_a_region_never_handed_out_is_refused() {
        let mut organiser = Organiser::new();
        organiser.submit(RegionIndex(7), outcome_of(&[region(1, 3)], &[]));
        let _ = organiser.finish(3);
    }
}
