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
//! **What has landed:** this file's five parameters; [`close`]'s walk and its two
//! verdicts; [`build`]'s assembly of a survivor — every member projected onto the locus
//! span, unified into one allele table, with each covering sample's support against it;
//! [`organise`]'s observation cache — the window one builder is handed — the division of the
//! analysed ground into the regions single builders own, the organiser's **ordered release**,
//! which takes the builders' outcomes by region index and lets their loci out along an
//! unbroken run of indexes, and its resolution of the overlaps between neighbouring builders;
//! [`serial`]'s **two** single-threaded drivers — the oracle every later milestone must
//! reproduce and the same merge read through the cache, byte for byte — plus, since the run
//! driver's E1, a third whose only departure is that each cover sweeps the samples
//! concurrently; and [`parallel`]'s
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
pub mod observation_cache;
pub mod organise;
pub mod parallel;
pub mod serial;
pub mod timing;

/// The fixtures the module's test suites share — the coordinates every test writes and the
/// failure every fake source yields.
///
/// **Two of these are read from outside the merge**, and deliberately: `render` and
/// `refuse_any_difference` are what "the same answer" means for a `RegionOutcome`, and the run
/// driver's own differential asks the same question of the same type. A second spelling of it
/// there is exactly what `render`'s destructuring exists to prevent — a field the outcome gains
/// would have to be answered for twice, and would silently drop out of one.
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
            reference_bases: reference_bases.to_vec(),
            observations: vec![SequenceObservation {
                bases: observed_bases.to_vec(),
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
        // **The deletion's record carries a partial sequence as well as its deletion**, and it
        // is the only one in the fixture. Without it no driver comparison in this module ever
        // looks at `SampleSupport::partials`: every other record here is minted `Complete`, so
        // the field renders empty in every entry and a driver that built it differently would
        // agree with the others by both building nothing. It is the field where such a
        // divergence would now hide, because a sample whose observations are all partial
        // contributes to no other one. Adding it changes no locus's existence — the keep rule
        // counts complete observations only — which is why it can go on an existing record.
        let mut deletion = member(region(305, 330), &[b'A'; 26], b"A");
        let mut ran_out = deletion.observations[0].clone();
        ran_out.bases = b"AAAAAT".to_vec();
        ran_out.num_obs = 2;
        ran_out.chain_ids = vec![4, 5];
        ran_out.q_sum = crate::ng::types::SummedLogError::from_nats(-4.5);
        ran_out.read_witness = ReadWitness::from_left(6, deletion.locus_len())
            .expect("six positions inside a twenty-six-base record");
        deletion.observations.push(ran_out);

        vec![
            dotted,
            vec![deletion],
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
    pub(in crate::ng::run) fn render(outcome: &super::build::RegionOutcome) -> Vec<String> {
        // Destructured, not field-accessed: this function is what "the same answer" means for
        // every comparison in the module, so a field `RegionOutcome` gains has to be answered
        // for here or it silently drops out of all of them.
        let super::build::RegionOutcome {
            cohort_observations,
            failed_locus_spans,
        } = outcome;
        cohort_observations
            .iter()
            .map(|observed| format!("{observed:?}"))
            .chain(
                failed_locus_spans
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
    pub(in crate::ng::run) fn refuse_any_difference(
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

/// The fewest non-reference **reads** one sample may show at a cohort locus and still
/// have it built — the floor half of [`MinAltReads`] (spec §4.3).
///
/// The `Obs` is production's word for a read's allele observation, and **not** this
/// module's `observation`, which is one sample's whole record over a stretch of genome
/// (spec §1.3). One observation can carry many of these reads.
///
/// Below the threshold in *every* sample the locus is **dropped**: nothing is assembled,
/// nothing is emitted, and — unlike a locus that failed [`MaxCohortLocusSpan`] — nothing
/// is counted. A failure is ground the caller refused; this is ground it judged empty,
/// and conflating the two would stop the failed count meaning anything.
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

/// 2 reads — production's number, over a rule that is now production's shape.
///
/// **ng asks the same question production asks, and asks it of a share as well.**
/// Production takes, across a group's positions, the *maximum over samples* of that
/// sample's non-reference observations (`max_nonref_obs`,
/// `cohort_integration.rs:64-78`, summed in `derive_is_kept`) — a per-sample test. ng
/// summed every sample's non-reference reads instead until 2026-08-19, and the owner
/// retired that rule when it was measured: because the bar was a fixed count and the sum
/// grows with the cohort, it stopped filtering as samples were added. Over one 100 kb
/// tomato interval it discarded 997 loci in 1,000 at one accession and only 439 in 1,000
/// at 63, where the merge then assembled a cohort observation at 546 positions per
/// kilobase.
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

/// What share of one sample's compared reads must be non-reference — the share half of
/// [`MinAltReads`] (spec §4.3).
///
/// A fraction of 1, not a percentage: 0.02 is two reads in a hundred.
///
/// **The denominator is the sample's own reads, and that is the whole design.** A share
/// taken of the *cohort's* reads would grow its bar as samples are added while a rare
/// allele's evidence stays where it is — one carrier's handful of reads — so a large
/// study could not call the rare variation it was run to find (spec §4.3, and the
/// derivation there).
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct MinAltReadShare(f64);

impl MinAltReadShare {
    /// The default share, [`DEFAULT_MIN_ALT_READ_SHARE`].
    pub const DEFAULT: Self = Self(DEFAULT_MIN_ALT_READ_SHARE);

    /// Whether `share` is a fraction of one — at least nothing, at most everything, and
    /// a real number.
    ///
    /// **One spelling for both constructors**, so that the fallible one and the `const`
    /// one cannot come to disagree about what a legal share is. That is not a
    /// hypothetical: a first draft of the pair wrote the test out twice, and a review's
    /// mutation of one copy left the other's tests green.
    ///
    /// The comparison is written out rather than asked of `(0.0..=1.0).contains(&share)`
    /// because `RangeInclusive::contains` is not a `const fn`. It refuses everything the
    /// range plus an `is_finite()` call would refuse: **a negative share and `-∞` fail
    /// the lower comparison, a share above one and `+∞` fail the upper, and `NaN` fails
    /// both.** `-0.0` is a share of nothing and passes, as it does through the range.
    const fn is_a_fraction_of_one(share: f64) -> bool {
        share >= 0.0 && share <= 1.0
    }

    /// The share, or `None` if it is not a fraction of 1 — negative, above one, or not a
    /// number. Rejecting rather than clamping: a mistyped share is a run that quietly
    /// keeps everything or nothing, and neither is visible in the output.
    pub fn new(share: f64) -> Option<Self> {
        Self::is_a_fraction_of_one(share).then_some(Self(share))
    }

    /// The same share, for a `const` that has to name one — and it **panics** where
    /// [`new`](Self::new) returns `None`.
    ///
    /// **The severity is the point, not a lapse.** It is written for `const`
    /// declarations of a rule the whole run is judged against, where a value outside
    /// `0..=1` fails the compile rather than the run and there is nobody to hand an
    /// `Option` to. Written for
    /// [`DEFAULT_MIN_ALLELE_SUPPORT`](crate::ng::calling::allele_candidates::DEFAULT_MIN_ALLELE_SUPPORT),
    /// which needs a share of 5 in 100 and cannot reach this type's private field.
    ///
    /// **Called at runtime it aborts the process** — `[profile.release]` sets
    /// `panic = "abort"` — and nothing in the signature stops that, because a `const fn`
    /// is an ordinary function too. **A share that is not a literal in this source tree —
    /// one an operator typed, one read from a file, one computed — goes through
    /// [`new`](Self::new)**, which refuses rather than aborts.
    pub const fn new_or_panic(share: f64) -> Self {
        assert!(
            Self::is_a_fraction_of_one(share),
            "a minimum non-reference read share is a fraction of one"
        );
        Self(share)
    }

    /// The share, as a fraction of 1.
    #[inline]
    pub const fn get(self) -> f64 {
        self.0
    }
}

impl Default for MinAltReadShare {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// 2 reads in a hundred — the owner's number, chosen against two measurements a
/// hundredfold apart in depth and not yet against a truth set.
///
/// **It is the half of the rule that only high depth can exercise.** At the tomato
/// panel's 11 reads a sample at a locus, two in a hundred asks for 0.22 reads and
/// [`MinAltObs`]'s floor of 2 decides everything — the share changed 12,034 kept loci to
/// 12,033. At the GIAB trio's 313 reads a sample it asks for 7, and there it does all
/// the work: 53,796 kept loci become 2,721 over 136 kb. A true heterozygote at that
/// depth shows about 156 non-reference reads, so the bar sits a factor of 20 below what
/// a real carrier produces.
pub const DEFAULT_MIN_ALT_READ_SHARE: f64 = 0.02;

/// **How much non-reference evidence one sample must show before the cohort builds a
/// locus** — the merge's keep rule, and the two numbers it is made of (spec §4.3).
///
/// A locus is built when **some single sample** shows at least
/// [`required_of`](Self::required_of) its compared reads as non-reference. Not the sum
/// across the cohort: a real variant concentrates its evidence in the samples that carry
/// it, so asking one sample is what tracks a carrier, and it is the only form of the
/// question whose answer does not drift as samples are added.
///
/// The two numbers answer the two ends of the depth axis, which is why neither alone
/// would do:
///
/// - the **floor** ([`MinAltObs`], 2 reads) is what decides at low coverage, where a
///   share of three reads rounds to nothing;
/// - the **share** ([`MinAltReadShare`], 2 in 100) is what decides at high coverage,
///   where two reads out of 300 is the error rate rather than an allele.
///
/// At one sample the rule is exactly what the cohort-sum rule it replaced was, because a
/// sum over one sample is that sample's own count.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Default)]
pub struct MinAltReads {
    /// The fewest non-reference reads that can ever suffice, whatever the depth.
    pub floor: MinAltObs,
    /// The share of a sample's compared reads asked for once depth makes it the larger.
    pub share: MinAltReadShare,
}

impl MinAltReads {
    /// The default rule: [`MinAltObs::DEFAULT`] reads or [`MinAltReadShare::DEFAULT`] of
    /// the sample's own, whichever is more.
    pub const DEFAULT: Self = Self {
        floor: MinAltObs::DEFAULT,
        share: MinAltReadShare::DEFAULT,
    };

    /// How many non-reference reads a sample that compared `reads_compared_with_reference`
    /// against the reference must show — the floor, or the share of those reads rounded
    /// **up**, whichever is larger.
    ///
    /// Rounded up because rounding down would let the share ask for less than it says: at
    /// 200 reads and a share of 2 in 100, 4 reads is the bar and 3 is not.
    ///
    /// The cast cannot go wrong at either end. The product is finite and never negative,
    /// so a share of 1 at the largest possible read count lands exactly on `u32::MAX`;
    /// anything past that saturates rather than wrapping, which is Rust's rule for
    /// float-to-integer casts since 1.45.
    #[inline]
    pub fn required_of(self, reads_compared_with_reference: u32) -> u32 {
        let asked = (self.share.get() * f64::from(reads_compared_with_reference)).ceil();
        self.floor.get().max(asked as u32)
    }

    /// Whether one sample that showed `non_reference_reads` out of
    /// `reads_compared_with_reference` reaches the rule — the question the merge asks of
    /// every contributing sample at every locus.
    ///
    /// **The numerator is the caller's, and the callers count different reads.** The
    /// merge passes a sample's non-reference reads pooled over the whole allele table,
    /// which is what the argument is named for. Candidate selection passes *one
    /// allele's* reads, so that the question "is this sequence worth calling over"
    /// is the merge's own question asked one level down
    /// (`doc/devel/ng/spec/candidate_alleles.md` §3). A discovery round passes a
    /// narrower count still — the reads the converged model is explaining as slippage
    /// (§3.4). The denominator is the same in all three: that sample's compared reads
    /// at the locus. **There is deliberately no wrapper per caller**: a second spelling
    /// of one rule is how two rules become different rules.
    ///
    /// **The floor is tested first because it is the cheap half and it cannot be
    /// skipped.** [`required_of`](Self::required_of) is never below the floor, so a
    /// sample under the floor fails whatever its depth — and answering it by integer
    /// comparison keeps the multiply off the path taken at almost every position, where
    /// the sample showed no non-reference read at all.
    #[inline]
    pub fn reached_by(self, non_reference_reads: u32, reads_compared_with_reference: u32) -> bool {
        non_reference_reads >= self.floor.get()
            && non_reference_reads >= self.required_of(reads_compared_with_reference)
    }
}

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

/// 500 bases — **raised from 200 on 2026-08-28**, on a sweep taken at 63 accessions rather than
/// the 16 the earlier number came from
/// ([`cohort_merge_parallel_cost_2026-08-28.md`](../../../../doc/devel/ng/research/cohort_merge_parallel_cost_2026-08-28.md) §5.1).
///
/// **What a region's width really decides is whether threads help at all.** Measured on the
/// tomato benchmark's 63 accessions over 100 kb of SL4.0 — real reads through the generic locus
/// generator, one record per covered position per sample — the merge on **one thread barely
/// notices the width**: 656 ms at 20 bases against 615 at 100, 616 at 200, 624 at 500 and 636 at
/// 1,000. On eight threads it notices a great deal, because what a narrow region adds is
/// per-round cost: at 20 bases, launching each round's builders is 12.0% of the merge and waiting
/// for its slowest is 5.0%; at 1,000 bases those are 0.5% and 3.9%.
///
/// **Eight threads at 63 accessions:** 579 ms at 20 bases, 406 at 100, 260 at 200, 220 at 500,
/// 219 at 1,000. **The ordering is what this number rests on, not those milliseconds** — five
/// sweeps under different host loads span a factor of two in absolute terms, and every one of
/// them puts 200 slower than 500 and 500 level with or slightly slower than 1,000.
///
/// **What it costs is the ground the observation cache holds**, and that was measured too: peak
/// resident 3.148 GB at 200 bases, 3.244 at 500, 3.343 at 1,000, 3.590 at 2,000. **500 buys the
/// speed for 3% more memory**; 1,000 is not reliably better and 2,000 is worse. With 16 regions
/// in flight the cache spans 8,000 bases rather than 3,200.
///
/// **The 16-sample sweep this replaces put the optimum at 100–200 bases**, and it was not wrong
/// about its own corner: fewer samples make each region's fixed cost smaller relative to the
/// building, which moves the optimum narrower. Denser ground moves it wider for the same reason.
/// So this number is a fact about tens of samples at about one record per covered base, and a
/// cohort far from that should sweep it again — it is a run parameter precisely so that costs no
/// code.
///
/// **The fabricated fixture that first argued for the old number argued too strongly**, and the
/// note it left is kept as a caution. Over ground with a record every hundred bases it made 200
/// look 1.7 times better than 20 and wider still better again; real observations arrive about
/// **one per covered base**, a hundred times denser.
pub const DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN: u32 = 500;

/// How many of those regions the merge works at once (spec §6.2).
///
/// **A count of regions in flight, not of threads.** The threads come from rayon's pool; what
/// this number sets is the ground the observation cache must hold, which is
/// `cohort_locus_builder_regions_in_flight × cohort_locus_builder_regions_len` bases plus the
/// tail of the observations reaching past it (spec §6.4, §8) — 3,200 bases at 16 regions in
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
        assert_eq!(DEFAULT_MIN_ALT_READ_SHARE, 0.02);
        assert_eq!(DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, 500);

        assert_eq!(MaxCohortLocusSpan::default().get(), 50);
        assert_eq!(MinAltObs::default().get(), 2);
        assert_eq!(MinAltReadShare::default().get(), 0.02);
        assert_eq!(CohortLocusBuilderRegionsLen::default().get(), 500);

        assert_eq!(MinAltReads::default(), MinAltReads::DEFAULT);
        assert_eq!(MinAltReads::DEFAULT.floor.get(), 2);
        assert_eq!(MinAltReads::DEFAULT.share.get(), 0.02);
    }

    /// **The two halves of the keep rule, each at the depth where it decides.** The
    /// crossing point is where the share first exceeds the floor, and it is the number
    /// this rule's behaviour hinges on: below it every sample is asked for the floor, and
    /// above it the bar climbs with the sample's own depth.
    #[test]
    fn the_floor_decides_at_low_depth_and_the_share_at_high() {
        let rule = MinAltReads::DEFAULT;

        // At and below 100 reads, 2% asks for no more than the floor of 2.
        assert_eq!(
            rule.required_of(0),
            2,
            "no reads at all still asks the floor"
        );
        assert_eq!(rule.required_of(3), 2, "2% of 3 rounds up to 1");
        assert_eq!(rule.required_of(100), 2, "2% of 100 is exactly the floor");

        // Above it the share takes over, and it rounds up rather than down: 2% of 150 is
        // 3 exactly, and 2% of 151 is 3.02, which asks for 4.
        assert_eq!(
            rule.required_of(101),
            3,
            "the first depth the share decides"
        );
        assert_eq!(rule.required_of(150), 3);
        assert_eq!(rule.required_of(151), 4, "rounded up, not down");
        assert_eq!(rule.required_of(313), 7, "the GIAB trio's depth a sample");
    }

    /// `reached_by` is `required_of` with the comparison folded in — pinned separately
    /// because it carries a short-circuit on the floor, and a short-circuit that answered
    /// differently from the arithmetic it is skipping would be invisible in the merge's
    /// output.
    #[test]
    fn reaching_the_rule_agrees_with_what_it_asks_for() {
        let rule = MinAltReads::DEFAULT;

        for compared in [0u32, 1, 3, 50, 100, 101, 150, 313, 1_000] {
            for alt in 0..=12u32 {
                assert_eq!(
                    rule.reached_by(alt, compared),
                    alt >= rule.required_of(compared),
                    "{alt} non-reference of {compared} compared"
                );
            }
        }
    }

    /// A share is a fraction of one, and anything else is refused rather than clamped —
    /// a run that silently kept everything or nothing would look like a run that found
    /// everything or nothing.
    #[test]
    fn a_share_outside_a_fraction_of_one_is_refused() {
        assert_eq!(
            MinAltReadShare::new(0.0).map(MinAltReadShare::get),
            Some(0.0)
        );
        assert_eq!(
            MinAltReadShare::new(1.0).map(MinAltReadShare::get),
            Some(1.0)
        );
        assert_eq!(
            MinAltReadShare::new(0.02).map(MinAltReadShare::get),
            Some(0.02)
        );

        assert!(MinAltReadShare::new(-0.01).is_none(), "below zero");
        assert!(MinAltReadShare::new(1.01).is_none(), "above one");
        assert!(MinAltReadShare::new(f64::NAN).is_none(), "not a number");
        assert!(MinAltReadShare::new(f64::INFINITY).is_none(), "not finite");
        assert!(
            MinAltReadShare::new(f64::NEG_INFINITY).is_none(),
            "not finite, the other way"
        );
    }

    /// **The `const` constructor refuses exactly what [`MinAltReadShare::new`] refuses**,
    /// which is the whole claim `is_a_fraction_of_one` exists to make true — one test,
    /// both constructors, over the same ten values.
    ///
    /// `new_or_panic` panics where `new` returns `None`, so the two are compared through
    /// [`std::panic::catch_unwind`] rather than by reading them side by side. The values
    /// that matter here are the ones no earlier fixture reached: **a negative share**,
    /// which does not crash anything downstream but silently deletes the share half of
    /// the rule — `required_of` casts a negative product to `u32`, which saturates to
    /// zero, so `max(floor, 0)` is the floor at every depth — and **both ends of the
    /// closed range**, which must be accepted rather than refused.
    #[test]
    fn the_const_share_refuses_exactly_what_the_fallible_one_refuses() {
        for share in [
            -1.0,
            -0.05,
            -0.0,
            0.0,
            0.05,
            1.0,
            1.5,
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::NAN,
        ] {
            let accepted_by_new = MinAltReadShare::new(share).is_some();
            let accepted_by_new_or_panic =
                std::panic::catch_unwind(|| MinAltReadShare::new_or_panic(share)).is_ok();
            assert_eq!(
                accepted_by_new, accepted_by_new_or_panic,
                "the two constructors disagree about a share of {share}"
            );
        }
    }

    /// The `const` constructor's message, so the two ends of its range are pinned as
    /// *refusals* and not merely as "something went wrong".
    ///
    /// **A negative share is the one that has to be loud.** Above one and `NaN` are
    /// visibly nonsense; a negative share reads as a plausible typo and, unrefused,
    /// produces a run that keeps every allele the floor admits at every depth with
    /// nothing in the output saying so.
    #[test]
    #[should_panic(expected = "a fraction of one")]
    fn a_const_share_below_zero_fails_rather_than_clamping() {
        let _ = MinAltReadShare::new_or_panic(-0.05);
    }

    #[test]
    #[should_panic(expected = "a fraction of one")]
    fn a_const_share_above_one_fails_rather_than_clamping() {
        let _ = MinAltReadShare::new_or_panic(1.5);
    }

    #[test]
    #[should_panic(expected = "a fraction of one")]
    fn a_const_share_that_is_not_a_number_fails() {
        let _ = MinAltReadShare::new_or_panic(f64::NAN);
    }

    /// Both ends of the closed range are legal shares — no read at all, and every read.
    #[test]
    fn a_const_share_may_sit_on_either_end_of_the_range() {
        assert_eq!(MinAltReadShare::new_or_panic(0.0).get(), 0.0);
        assert_eq!(MinAltReadShare::new_or_panic(1.0).get(), 1.0);
    }

    /// **A share of one at the largest read count lands on `u32::MAX` and does not
    /// wrap.** The cast in `required_of` is the only arithmetic in this module that could
    /// answer a number smaller than the truth, and a keep rule that asked for *less* at
    /// the extreme would build every locus in a run that had asked for none.
    #[test]
    fn the_share_at_its_extreme_saturates_rather_than_wrapping() {
        let everything = MinAltReads {
            floor: MinAltObs::DEFAULT,
            share: MinAltReadShare::new(1.0).expect("a fraction of one"),
        };

        assert_eq!(everything.required_of(u32::MAX), u32::MAX);
        assert!(!everything.reached_by(u32::MAX - 1, u32::MAX));
        assert!(everything.reached_by(u32::MAX, u32::MAX));
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
    /// its own default would run every cohort at 50/2/200 whatever the operator asked
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
