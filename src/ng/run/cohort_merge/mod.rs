//! The cohort merge — k samples' locus observations into cohort observations, in
//! parallel.
//!
//! Every sample's observations arrive in coordinate order. This module groups the
//! positions the cohort varied at into **cohort loci**, judges each one, and
//! assembles the survivors into the unit the caller consumes. It is the stage that
//! ran on a single thread in production and was the wall floor there
//! (`cohort_integration.rs`); here the genome is dealt out in short regions and the
//! builders working them share nothing but their results.
//!
//! Design: `doc/devel/ng/spec/cohort_merge.md` (what and why),
//! `doc/devel/ng/arch/cohort_merge.md` (types and contracts).
//!
//! **What has landed:** this file's four parameters; [`close`]'s walk and its two
//! verdicts; [`build`]'s assembly of a survivor — every member projected onto the locus
//! span, unified into one allele table, with each covering sample's support against it;
//! [`organise`]'s observation cache — the window one builder is handed — the division of the
//! analysed ground into the regions single builders own, the organiser's **ordered release**,
//! which takes the builders' outcomes by region index and lets their loci out along an
//! unbroken run of indexes, and its resolution of the overlaps between neighbouring builders;
//! [`serial`]'s **two** single-threaded drivers, the oracle every later milestone must
//! reproduce and the same merge read through the cache, byte for byte; and [`parallel`]'s
//! round-by-round arrangement of the builders, whose output is those two byte for byte at any
//! number of regions in flight and any region width.
//!
//! **Still to come:** the caller objects that will own these parameters and this cache, and
//! the run summary that reports what the merge refused (spec §13).
//!
//! **`pub`, though the architecture calls this crate-private machinery.** The two
//! caller objects that will own it do not exist yet, so `pub(crate)` items here would
//! have no consumer and the crate's `-D warnings` gate would reject them as dead code;
//! ng's probes also live in `examples/`, which are separate crate targets and see only
//! `pub` items. The intent is the architecture's — narrow this when the caller objects
//! land.

pub mod build;
pub mod close;
pub mod organise;
pub mod parallel;
pub mod serial;

/// The fixtures the module's test suites share — the coordinates every test writes and the
/// failure every fake source yields.
///
/// **One home, because the copies had started to multiply.** `region` and `region_on` are
/// written out in all four of this module's files, and D2's review found `SourceFailed`
/// becoming the fifth such copy. [`organise`] and [`serial`] read them from here; [`build`] and
/// [`close`] still carry their own and are the next two to fold in.
#[cfg(test)]
pub(super) mod fixtures {
    use crate::ng::locus_generation::{
        LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
    };
    use crate::ng::types::{ContigId, GenomePosition, GenomeRegion, Position, ReadGroupId};

    /// A region on the named contig, both ends inclusive.
    pub(super) fn region_on(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    /// A region on contig 0 — the one most fixtures need.
    pub(super) fn region(start: u64, end: u64) -> GenomeRegion {
        region_on(0, start, end)
    }

    /// One base, genome-wide.
    pub(super) fn position_on(contig: u32, position: u64) -> GenomePosition {
        GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        }
    }

    /// One sample's record over `region`, showing `observed_bases` to three reads.
    ///
    /// **Here rather than in one test module, because three drivers want it** — the two serial
    /// ones, and the parallel one that must reproduce them. A fixture that differs between the
    /// files comparing their outputs is a difference nobody would look for.
    pub(super) fn member(
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

    /// One sample's forward reader over the observations it minted.
    pub(super) fn source_of(
        observations: &[SampleLocusObservations],
    ) -> std::vec::IntoIter<Result<SampleLocusObservations, SourceFailed>> {
        observations
            .iter()
            .cloned()
            .map(Ok)
            .collect::<Vec<_>>()
            .into_iter()
    }

    /// A building-region width in reference bases.
    pub(super) fn width(bases: u32) -> super::CohortLocusBuilderRegionsLen {
        super::CohortLocusBuilderRegionsLen(
            std::num::NonZeroU32::new(bases).expect("a fixture width is non-zero"),
        )
    }

    /// How many building regions a fixture works at once.
    pub(super) fn in_flight(regions: usize) -> super::CohortLocusBuilderRegionsInFlight {
        super::CohortLocusBuilderRegionsInFlight(
            std::num::NonZeroUsize::new(regions).expect("a fixture region count is non-zero"),
        )
    }

    /// **The fixture the module's byte-identity claims rest on**, and the reason it is here:
    /// three drivers compare their outputs on it, and a copy that drifted would leave each
    /// file passing while they merged different ground.
    ///
    /// Sixty single-base loci every ten bases, a deletion at 305–330 that carries a locus
    /// across whatever boundary a region width puts near it, and a third sample whose one
    /// record sits inside that deletion — so two samples chain into one locus however finely
    /// the ground is divided. The whole fixture's reference is `A`, because a locus gathers
    /// its reference from its members and two members that disagree over shared ground are
    /// refused: the deletion's ground covers two of the dotted sample's positions.
    pub(super) fn three_samples_over_six_hundred_bases() -> Vec<Vec<SampleLocusObservations>> {
        let dotted: Vec<_> = (0..60)
            .map(|locus| {
                let at = 10 * locus + 1;
                member(region(at, at), b"A", b"T")
            })
            .collect();
        vec![
            dotted,
            vec![member(region(305, 330), &[b'A'; 26], b"A")],
            vec![member(region(310, 310), b"A", b"C")],
        ]
    }

    /// A whole outcome rendered entry by entry — what "the same answer" means across this
    /// module's drivers.
    ///
    /// **The comparison is on the `Debug` rendering**: `CohortObservation` has no `PartialEq`,
    /// a comparison written field by field would silently stop covering a field added later,
    /// and two distinct `f64` sums render as distinct strings, so a quality divided differently
    /// shows.
    pub(super) fn render(outcome: &super::build::RegionOutcome) -> Vec<String> {
        outcome
            .cohort_observations
            .iter()
            .map(|observed| format!("{observed:?}"))
            .chain(
                outcome
                    .failed_locus_spans
                    .iter()
                    .map(|span| format!("failed {span}")),
            )
            .collect()
    }

    /// Refuse any difference between two rendered outcomes, naming the **first** entry that
    /// differs and `what_changed` as the setting that produced it.
    ///
    /// **Entry by entry rather than whole**: the widest fixture in this module renders as 26 kB,
    /// and two of those inside one `assert_eq!` message is not something a reader can diff by
    /// eye.
    pub(super) fn refuse_any_difference(
        what_changed: &str,
        expected: &super::build::RegionOutcome,
        actual: &super::build::RegionOutcome,
    ) {
        let (expected, actual) = (render(expected), render(actual));
        if let Some(first) = expected
            .iter()
            .zip(&actual)
            .position(|(one, other)| one != other)
        {
            panic!(
                "{what_changed} changed the merge, first at entry {first}:\
                 \n  expected: {}\n  actual:   {}",
                expected[first], actual[first],
            );
        }
        assert_eq!(
            actual.len(),
            expected.len(),
            "{what_changed} changed how many entries the merge produced",
        );
    }

    /// What a source failure looks like in a test. The run's own error type does not exist
    /// yet, so the cache is generic over the source's and both drivers pass it through
    /// (`doc/devel/ng/arch/run_streaming.md` §2, §5).
    #[derive(Debug, PartialEq, Eq)]
    pub(super) struct SourceFailed(pub &'static str);
}

use std::num::{NonZeroU32, NonZeroUsize};

use crate::ng::types::GenomeRegion;

/// The widest cohort locus the caller undertakes to build, in reference bases.
///
/// A locus wider than this is **failed**: not assembled, not emitted, and counted in
/// the run summary, while its ground still displaces the loci that overlap it
/// (spec §3.2, the owner's 2026-08-17 ruling). It is a policy bound and not a fact
/// about the data, which is why it is the operator's to set: a run over long reads
/// is expected to raise it, since the widest event worth merging into one locus
/// grows with the reads.
///
/// **A command-line parameter of a calling run**, default
/// [`DEFAULT_MAX_COHORT_LOCUS_SPAN`]. It is never recorded in a psp file — those hold
/// what the generator minted — so re-calling under a new value needs no second walk
/// over the alignments (spec §3.1).
///
/// **Owed:** the effective value belongs in the run's output beside the failed-locus
/// count, because it decides which ground was refused and two runs over the same
/// records under different values are otherwise indistinguishable (arch §1; spec §3.1,
/// §3.3). Where the run summary lives is the emission step's (spec §13), so nothing
/// here writes it yet.
///
/// It governs **generic** loci. An STR locus's span is its reference tract, which the
/// segmentation defines and which may be wider (spec §3.1).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct MaxCohortLocusSpan(pub NonZeroU32);

impl MaxCohortLocusSpan {
    /// The default bound, [`DEFAULT_MAX_COHORT_LOCUS_SPAN`] reference bases — typed, so
    /// a call site that wants to name what it is passing can.
    pub const DEFAULT: Self = Self(non_zero_default(DEFAULT_MAX_COHORT_LOCUS_SPAN));

    /// The bound in reference bases.
    #[inline]
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}

impl Default for MaxCohortLocusSpan {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// 50 bases — the owner's number, unmeasured and soft (spec §14 question 3).
///
/// Cheap to revisit: re-calling under a different bound needs no re-walk (spec §3.1),
/// so the measurement that would settle it — how much real signal sits just above 50
/// — can be made against records that already exist.
pub const DEFAULT_MAX_COHORT_LOCUS_SPAN: u32 = 50;

/// How many non-reference **reads** a cohort locus needs, summed across it, to be built
/// at all (spec §4.3).
///
/// The `Obs` is production's word for a read's allele observation, and **not** this
/// module's `observation`, which is one sample's whole record over a stretch of genome
/// (spec §1.3). One observation can carry many of these reads.
///
/// Below the threshold the locus is **dropped**: nothing is assembled, nothing is
/// emitted, and — unlike a locus that failed [`MaxCohortLocusSpan`] — nothing is
/// counted. A failure is ground the caller refused; this is ground it judged empty, and
/// conflating the two would stop the failed count meaning anything.
///
/// **A command-line parameter**, default [`DEFAULT_MIN_ALT_OBS`], and the reason it
/// exists is measured rather than aesthetic: in production, under production's rule
/// (below), it removes a large number of very low-quality variants and improves
/// performance substantially. Its cost is that a variant seen once, in one sample, is
/// unrecoverable — nothing downstream is emitted for that locus at all.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct MinAltObs(pub NonZeroU32);

impl MinAltObs {
    /// The default threshold, [`DEFAULT_MIN_ALT_OBS`] non-reference reads.
    pub const DEFAULT: Self = Self(non_zero_default(DEFAULT_MIN_ALT_OBS));

    /// The threshold in non-reference reads.
    #[inline]
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}

impl Default for MinAltObs {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// 2 reads — production's number, over a rule that is not production's.
///
/// **The rule differs, so the number does not mean quite the same thing.** Production
/// sums, across a group's positions, the *maximum over samples* of that sample's
/// non-reference observations (`max_nonref_obs`, `cohort_integration.rs:64-78`, summed
/// in `derive_is_kept`); ng sums every sample's non-reference reads, which spec §4.3
/// chose deliberately and spec §15 pins with a test inverting production's. A maximum
/// is never larger than a sum, so at the same threshold ng keeps everything production
/// keeps and more: identically at one sample, and by a widening margin as the cohort
/// grows — at 63 samples, one non-reference read in each of two samples at one position
/// reaches ng's 2 and never reaches production's. The performance claim above was
/// measured under production's rule, so it is inherited evidence rather than evidence
/// about ng.
///
/// The value is copied and the name is not. Production spells it
/// [`DEFAULT_MIN_ALT_OBS_PER_SAMPLE`], the default of `--min-alt-obs-per-sample`, and
/// feeds that one number both to the cohort keep rule ng's descends from
/// (`derive_is_kept`) and to a per-sample pre-EM filter (`passes_min_alt_obs`,
/// `variant_caller.rs`). ng's threshold is only the first, so reaching the constant by
/// name would let a retune of production's per-sample filter move ng's cohort keep.
///
/// [`DEFAULT_MIN_ALT_OBS_PER_SAMPLE`]: crate::var_calling::DEFAULT_MIN_ALT_OBS_PER_SAMPLE
pub const DEFAULT_MIN_ALT_OBS: u32 = 2;

/// How many reference bases one builder's region covers (spec §6.1).
///
/// **A command-line parameter**, default [`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`].
/// It is deliberately *not* derived from [`MaxCohortLocusSpan`]: what a region's width
/// really costs is the observation cache, which has to cover every region in play at
/// once, so `builders × this` is the ground held resident (spec §6.4, §8). How wide a
/// locus may be has nothing to do with it.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CohortLocusBuilderRegionsLen(pub NonZeroU32);

impl CohortLocusBuilderRegionsLen {
    /// The default width, [`DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN`] reference bases.
    pub const DEFAULT: Self = Self(non_zero_default(DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN));

    /// The region width in reference bases.
    #[inline]
    pub const fn get(self) -> u32 {
        self.0.get()
    }
}

impl Default for CohortLocusBuilderRegionsLen {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// 20 bases — the owner's starting value, unmeasured (spec §14 question 1).
///
/// The sweep that would settle it trades two things against each other: wider regions
/// mean fewer joins between builders and so less overlapping work discarded, while
/// narrower ones shrink the ground the observation cache must cover, which is this
/// module's main memory.
pub const DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN: u32 = 20;

/// How many of those regions the merge works at once (spec §6.2).
///
/// **A count of regions in flight, not of threads.** The threads come from rayon's pool; what
/// this number sets is the ground the observation cache must hold, which is
/// `cohort_locus_builder_regions_in_flight × cohort_locus_builder_regions_len` bases plus the
/// tail of the observations reaching past it (spec §6.4, §8) — 320 bases at 16 regions in
/// flight and this module's default width. One region in flight gives the serial driver's
/// memory and its answer.
///
/// **A command-line parameter, and the only one of this module's four with no constant
/// default**: what it should be depends on the machine's cores and on how much memory the
/// cohort's width leaves, neither of which is knowable where the other three defaults are
/// written. A run given no value takes [`one_per_worker_thread`](Self::one_per_worker_thread).
///
/// **Owed:** the resolved value belongs in the run's output beside the other three, because
/// two runs differing only in it differ in memory and in nothing a reader of the output can
/// see.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CohortLocusBuilderRegionsInFlight(pub NonZeroUsize);

impl CohortLocusBuilderRegionsInFlight {
    /// One region in flight per thread in rayon's pool — what a run takes when the operator
    /// names no value.
    ///
    /// Never zero, and not by an assertion: rayon's pool has at least one thread, and the
    /// fallback says what to do if that ever stops being true rather than panicking part-way
    /// through a run.
    #[must_use]
    pub fn one_per_worker_thread() -> Self {
        Self(NonZeroUsize::new(rayon::current_num_threads()).unwrap_or(NonZeroUsize::MIN))
    }

    /// How many regions are worked at once.
    #[inline]
    pub const fn get(self) -> usize {
        self.0.get()
    }
}

/// Refuse the analysed regions no driver in this module can merge — see
/// [`merge_cohort_serially`](serial::merge_cohort_serially), whose doc says what a duplicate
/// would cost.
///
/// **Two rules, and the second was found by D2's review.** A region must not overlap or precede
/// the one before it, or the loci in the ground they share are built — and carried — twice. And
/// a region's own ends must be the right way round: the division orders them
/// ([`building_regions_of`](organise::building_regions_of)) and
/// [`build_region`](build::build_region) does not, so on `50-1` the cached driver builds the
/// whole of 1–50 and the oracle builds nothing. That is the byte-identity this module claims,
/// broken by an input nothing was checking.
///
/// **Here rather than in one driver's file, because it is the contract all three share.** Every
/// driver opens by calling it, so a reader of any of them finds the rule where the module's
/// other shared terms are, not in the file named for a different arrangement.
pub(super) fn refuse_malformed_analysed_regions(analysed: &[GenomeRegion]) {
    for region in analysed {
        assert!(
            region.start <= region.end,
            "the analysed region {region} has its ends the wrong way round, and the drivers \
             do not read such a region the same way",
        );
    }
    for pair in analysed.windows(2) {
        let (earlier, later) = (pair[0], pair[1]);
        assert!(
            (earlier.contig, earlier.end) < (later.contig, later.start),
            "the analysed regions {earlier} and {later} are not disjoint and ascending, so \
             the loci opening in the ground they share would be built — and carried — twice",
        );
    }
}

/// A default's value as the [`NonZeroU32`] the newtypes hold.
///
/// **A zero default is a build error, not a panic**, and it is the *call* that makes it
/// one: every caller is a `pub const DEFAULT: Self` item, which is evaluated when the
/// crate is compiled, so the `None` arm below never survives into a running binary.
/// Anything added later gets the same guarantee from the same line — there is no
/// separate assertion to remember to write.
///
/// Writing the defaults as `u32` rather than as `NonZeroU32` keeps them readable at
/// their declaration, which is where an operator reading the source looks for them.
const fn non_zero_default(default_value: u32) -> NonZeroU32 {
    match NonZeroU32::new(default_value) {
        Some(value) => value,
        None => panic!("a cohort-merge default must be non-zero"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Both spellings of each default — the constant a command line will advertise and
    /// the typed value a run will use — pinned to the documented number, so the two
    /// cannot drift into printing one and running the other.
    #[test]
    fn the_defaults_are_the_documented_values() {
        assert_eq!(DEFAULT_MAX_COHORT_LOCUS_SPAN, 50);
        assert_eq!(DEFAULT_MIN_ALT_OBS, 2);
        assert_eq!(DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, 20);

        assert_eq!(MaxCohortLocusSpan::default().get(), 50);
        assert_eq!(MinAltObs::default().get(), 2);
        assert_eq!(CohortLocusBuilderRegionsLen::default().get(), 20);
    }

    /// The fourth parameter's default is a rule rather than a number, so the rule is what is
    /// pinned: one region in flight per thread in rayon's pool, whatever that machine's pool
    /// turns out to be. A default that ignored the pool — one, or a constant 16 — would pass
    /// every other test in this module, since the merge's answer is the same at every count.
    #[test]
    fn the_regions_in_flight_default_is_one_per_worker_thread() {
        assert_eq!(
            CohortLocusBuilderRegionsInFlight::one_per_worker_thread().get(),
            rayon::current_num_threads().max(1),
        );
    }

    /// An operator-set value is the only case that tells reading the field apart from
    /// returning the default — and it is the case these parameters exist for, since all
    /// three are set from a command line and spec §3.1 expects a long-read run to raise
    /// the bound. Without this, an accessor that ignored its argument and answered with
    /// its own default would run every cohort at 50/2/20 whatever the operator asked
    /// for, and no test would notice.
    #[test]
    fn get_returns_the_wrapped_value_not_the_default() {
        assert_eq!(MaxCohortLocusSpan(NonZeroU32::new(200).unwrap()).get(), 200);
        assert_eq!(MinAltObs(NonZeroU32::new(7).unwrap()).get(), 7);
        assert_eq!(
            CohortLocusBuilderRegionsLen(NonZeroU32::new(100).unwrap()).get(),
            100
        );

        // Both ends of the type the newtypes wrap.
        assert_eq!(MinAltObs(NonZeroU32::MIN).get(), 1);
        assert_eq!(MaxCohortLocusSpan(NonZeroU32::MAX).get(), u32::MAX);
    }
}
