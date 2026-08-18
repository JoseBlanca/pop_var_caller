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
//! **What has landed:** this file's three parameters; [`close`]'s walk and its two
//! verdicts; [`build`]'s assembly of a survivor — every member projected onto the locus
//! span, unified into one allele table, with each covering sample's support against it;
//! [`organise`]'s observation cache — the window one builder is handed — the division of the
//! analysed ground into the regions single builders own, and the organiser's **ordered
//! release**, which takes the builders' outcomes by region index and lets their loci out
//! along an unbroken run of indexes; and [`serial`]'s **two** single-threaded drivers: the
//! oracle every later milestone must reproduce, and the same merge read through the cache,
//! byte for byte. **Still to come:** the organiser's resolution of the overlaps between
//! neighbouring builders, and the parallel arrangement around it
//! (`doc/devel/ng/impl_plan/cohort_merge.md`, E2 and E3).
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
    use crate::ng::types::{ContigId, GenomePosition, GenomeRegion, Position};

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

    /// What a source failure looks like in a test. The run's own error type does not exist
    /// yet, so the cache is generic over the source's and both drivers pass it through
    /// (`doc/devel/ng/arch/run_streaming.md` §2, §5).
    #[derive(Debug, PartialEq, Eq)]
    pub(super) struct SourceFailed(pub &'static str);
}

use std::num::NonZeroU32;

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
