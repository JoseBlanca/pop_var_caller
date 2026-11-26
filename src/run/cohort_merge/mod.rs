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
//! concurrently; and `parallel`'s
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
pub mod recorded_windows;
pub mod rounds;
pub mod serial;
pub mod stream;

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
    use crate::locus_generation::{
        LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
    };
    use crate::types::{ContigId, GenomePosition, GenomeRegion, Position, ReadGroupId};

    use super::observation_cache;

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
            chain_id_base: 0,
            region,
            reference_bases: Box::from(reference_bases),
            observations: vec![SequenceObservation {
                flank_deletion: crate::locus_generation::FlankDeletion::NONE,
                held_flank: None,
                bases: Box::from(observed_bases),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 3,
                num_fwd: 3,
                q_sum: crate::types::SummedLogError::from_nats(-6.0),
                mapq_sum: 180,
                mapq_sum_sq: 10_800,
                placed_left: 0,
                num_low_quality: 0,
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

    /// **One sample's reader that hands over summaries and keeps its records**, the way a
    /// reader over a stored `.psp` does — and the shape every other fixture in this module
    /// lacks.
    ///
    /// # Why this exists, and it is not for coverage
    ///
    /// [`source_of`] returns each record whole, so the cache's `keeps_evidence` is false, its
    /// windows carry the records themselves and
    /// [`ObservationSource::build`](observation_cache::ObservationSource::build) is never
    /// called. **Every byte-for-byte comparison this module makes between its drivers was made
    /// over sources of that shape**, and a defect reachable only through `build` therefore
    /// passed all of them: the member indices a closed locus carries are numbered inside the
    /// window the cache hands out, while `build_at` indexes a sample's whole held list, and the
    /// two agree only where the window starts at the head. A driver evicting at each region's
    /// own first base gets that for free; one evicting once for a round of regions does not,
    /// and every region after a round's first would have built its members from the region
    /// before it. It was found by running, not by the suite.
    ///
    /// So this source keeps its records and hands out `Drawn::Kept`, with the body range
    /// naming the record's own index — `i..i+1`, which is meaningful to nobody but itself,
    /// which is exactly what the trait says a body range is.
    pub(super) struct KeepingSource {
        records: Vec<SampleLocusObservations>,
        handed_out: usize,
    }

    impl KeepingSource {
        /// A reader over `observations`, keeping every one of them.
        pub(super) fn over(observations: &[SampleLocusObservations]) -> Self {
            Self {
                records: observations.to_vec(),
                handed_out: 0,
            }
        }
    }

    impl observation_cache::ObservationSource for KeepingSource {
        type Error = SourceFailed;

        /// Never called: this source overrides `next_drawn`, and the cache calls that.
        fn next_observation(
            &mut self,
            _spare: Option<SampleLocusObservations>,
        ) -> Option<Result<SampleLocusObservations, Self::Error>> {
            unreachable!("a keeping source is drawn through `next_drawn`")
        }

        fn next_drawn(
            &mut self,
            _spare: Option<SampleLocusObservations>,
        ) -> Option<Result<observation_cache::Drawn, Self::Error>> {
            let at = self.handed_out;
            let record = self.records.get(at)?;
            self.handed_out += 1;
            Some(Ok(observation_cache::Drawn::Kept {
                summary: observation_cache::LocusSummary::of(record),
                body: at..at + 1,
            }))
        }

        /// **Nothing is released**, which the trait permits and which this fixture wants: a
        /// source that dropped what the merge had passed could not answer a `build` the driver
        /// got wrong, and answering it wrongly is the failure this fixture exists to catch.
        fn build(
            &self,
            body: core::ops::Range<usize>,
        ) -> Result<SampleLocusObservations, Self::Error> {
            Ok(self.records[body.start].clone())
        }

        /// **Seals only what a real source could still build**: the records drawn so far, from
        /// `from_body` on. A driver that asked a seal for a body before its round, or one not
        /// drawn yet, gets a panic rather than a record.
        fn seal(
            &mut self,
            from_body: usize,
        ) -> Option<std::sync::Arc<dyn observation_cache::BuildsKeptEvidence<SourceFailed>>>
        {
            let records = self
                .records
                .iter()
                .enumerate()
                .map(|(at, record)| {
                    (at >= from_body && at < self.handed_out).then(|| record.clone())
                })
                .collect();
            Some(std::sync::Arc::new(SealedKeeping { records }))
        }
    }

    /// What a [`KeepingSource`] sealed: its records drawn by then from the seal's first body on.
    struct SealedKeeping {
        records: Vec<Option<SampleLocusObservations>>,
    }

    impl observation_cache::BuildsKeptEvidence<SourceFailed> for SealedKeeping {
        fn build(
            &self,
            body: core::ops::Range<usize>,
        ) -> Result<SampleLocusObservations, SourceFailed> {
            Ok(self.records[body.start]
                .clone()
                .unwrap_or_else(|| panic!("body {body:?} was not sealed with this round")))
        }
    }

    /// **The width bound these fixtures were drawn against: 50 bases**, named rather than taken
    /// from [`super::MaxCohortLocusSpan::DEFAULT`].
    ///
    /// The drivers' comparisons rest on a fixture locus the bound refuses — a 91-base span, so
    /// the refused list is compared as well as the built loci. When the default rose to 100
    /// (2026-10-04) that span stopped being refused and three tests found their fixture no
    /// longer carried what they compare. The bound a fixture exercises is a fact about the
    /// fixture, not about the shipped default.
    pub(super) const FIXTURE_BOUND: super::MaxCohortLocusSpan =
        super::MaxCohortLocusSpan(super::non_zero_default(50));

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
        ran_out.bases = Box::from(&b"AAAAAT"[..]);
        ran_out.num_obs = 2;
        ran_out.chain_ids = vec![4, 5];
        ran_out.q_sum = crate::types::SummedLogError::from_nats(-4.5);
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
    /// **Every field is rendered by `Debug` and the fields are named one at a time**, which is
    /// what lets one of them be left out. `CohortObservation` has no `PartialEq`, and two distinct
    /// `f64` sums render as distinct strings, so a quality divided differently shows. The
    /// destructure is what keeps the exclusion from spreading: a field this type gains is a
    /// compile error here, and whoever adds it has to say whether the drivers must agree on it.
    ///
    /// **The window coverage is the field left out, and leaving it out is a decision.** Only the
    /// driver reading through the observation cache takes that measurement; the in-memory driver
    /// it is compared against holds every record at once and takes none, so its every pair is
    /// absent. Rendering it would make these comparisons fail on a difference that is not a
    /// disagreement. **The tempting repair is worse**: asserting only that both sides are absent
    /// passes whatever the cached driver did, including measuring nothing at all. What checks the
    /// field is
    /// `serial::tests::the_cached_driver_measures_windows_where_the_in_memory_one_has_none`,
    /// which says out loud that one side has numbers and the other does not.
    pub(in crate::run) fn render(outcome: &super::build::RegionOutcome) -> Vec<String> {
        // Destructured, not field-accessed: this function is what "the same answer" means for
        // every comparison in the module, so a field `RegionOutcome` gains has to be answered
        // for here or it silently drops out of all of them.
        let super::build::RegionOutcome {
            cohort_observations,
            failed_locus_spans,
        } = outcome;
        cohort_observations
            .iter()
            .map(|observed| {
                // Destructured for the same reason, one level down: a field `CohortObservation`
                // gains is a compile error here rather than a comparison that quietly stops
                // covering it.
                let crate::cohort_observation::CohortObservation {
                    region,
                    alleles,
                    per_sample,
                    kind,
                    window_coverage: _,
                    allele_flank_deletions,
                } = observed;
                format!("{region:?} {alleles:?} {allele_flank_deletions:?} {per_sample:?} {kind:?}")
            })
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
    pub(in crate::run) fn refuse_any_difference(
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

    /// A fixture's source failure can also carry the cache's own — the reference fetch it makes
    /// once per cover — because [`cover`](super::observation_cache::ObservationCache::cover)
    /// asks every caller for the conversion.
    ///
    /// **A fixture reference is built in memory, so the only way to reach this is to ask it for
    /// ground it does not hold** — a contig past its four, or a position past its 2,000 bases —
    /// which the cache's own reference tests do deliberately and no other fixture does by
    /// accident. What the ground *was* is thrown away here, so a test that cares which ground
    /// failed uses [`ReferenceFailureRecorded`] instead.
    impl From<super::observation_cache::ReferenceUnreadable> for SourceFailed {
        fn from(_: super::observation_cache::ReferenceUnreadable) -> Self {
            Self("the fixture reference could not serve the cover's ground")
        }
    }

    /// A source failure that keeps the ground a failed fetch named, for the one test that is
    /// about the naming.
    ///
    /// **[`SourceFailed`] cannot carry it** — it holds a `&'static str` — so a test asserting on
    /// the region would be asserting on a constant, and any change to the ground a failure
    /// names would pass unnoticed.
    #[derive(Debug, PartialEq, Eq)]
    pub(super) struct ReferenceFailureRecorded(pub Option<GenomeRegion>);

    impl From<super::observation_cache::ReferenceUnreadable> for ReferenceFailureRecorded {
        fn from(failure: super::observation_cache::ReferenceUnreadable) -> Self {
            Self(Some(failure.region))
        }
    }

    /// A reference the cache's fixtures cover against: four contigs of 2,000 bases, the widest
    /// any fixture here reaches being 600.
    ///
    /// **The bases repeat `ACGT`**, so a window over any stretch of it has a GC fraction of a
    /// half — a value a test can predict without counting, and one no arithmetic slip lands on
    /// by accident the way `0` or `1` would.
    pub(super) fn a_reference() -> Box<dyn super::observation_cache::MergeReference + Send + Sync> {
        let bases: Vec<u8> = b"ACGT".iter().copied().cycle().take(2_000).collect();
        Box::new(crate::ref_seq::InMemoryRefSeq::from_contigs(vec![bases; 4]))
    }

    /// A reference that answers like [`a_reference`] and counts what it was told to release —
    /// the only way a test can see that the merge lets go of what it has walked past, since an
    /// in-memory reference holds nothing to release.
    pub(super) struct ReferenceCountingItsReleases {
        bases: crate::ref_seq::InMemoryRefSeq,
        pub(super) released_before: std::sync::Mutex<Vec<u64>>,
    }

    impl crate::ref_seq::RefSeq for ReferenceCountingItsReleases {
        fn fetch_into(
            &self,
            contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), crate::ref_seq::RefSeqError> {
            self.bases.fetch_into(contig, start_1based, length, dst)
        }
    }

    /// The lengths the cover's look-ahead stops at, forwarded from the bases behind this
    /// counter so that the two cannot describe different contigs.
    impl crate::ref_seq::ContigTable for ReferenceCountingItsReleases {
        fn contigs(&self) -> &crate::fasta::ContigList {
            self.bases.contigs()
        }
    }

    impl crate::ref_seq::EvictableRefSeq for ReferenceCountingItsReleases {
        fn evict_before(&self, pos: u64) {
            self.released_before
                .lock()
                .expect("a fixture reference is used by one test at a time")
                .push(pos);
        }

        fn resident_bases(&self) -> usize {
            0
        }
    }

    impl ReferenceCountingItsReleases {
        pub(super) fn new() -> std::sync::Arc<Self> {
            let bases: Vec<u8> = b"ACGT".iter().copied().cycle().take(2_000).collect();
            std::sync::Arc::new(Self {
                bases: crate::ref_seq::InMemoryRefSeq::from_contigs(vec![bases; 4]),
                released_before: std::sync::Mutex::new(Vec::new()),
            })
        }
    }
}

use std::num::{NonZeroU32, NonZeroUsize};

use crate::cohort_observation::non_zero_default;
use crate::types::GenomeRegion;

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

/// 100 bases since 2026-10-04; it was 50, the owner's unmeasured number (spec §14 question 3).
///
/// **Raised on two measurements, one at each end of the depth range** (owner, 2026-10-04):
///
/// - **63 tomato accessions at about 3×**, tomato1's 2 Mb. At 50, 278 loci were refused; at 100,
///   267 of them were built and 162 records added, with the same run time (37 s) and peak memory
///   (about 290 MB). Against GATK run on the same accessions over those loci, 1,506 of the
///   variant genotypes they add are ones GATK also calls (same allele, either zygosity) and 394
///   are not: 21 in 100 of the new calls disputed, against 31 in 100 at ordinary loci 20 to 50
///   bases wide.
/// - **GIAB HG002, HG003 and HG004 at 300×**: the one locus refused at 50 was a true 31-base
///   heterozygous deletion, one copy of a two-copy tandem duplication, which can sit anywhere
///   along 62 bases. At 100 it is called, and nothing else changes.
///
/// A deletion of one copy of a two-copy duplication needs a locus twice its length, so 50 ruled
/// out every such deletion longer than 25 bases. Re-calling under a different bound needs no
/// re-walk (spec §3.1). Details: `doc/devel/ng/research/giab_unexplained_fn_2026-10-03.md`.
pub const DEFAULT_MAX_COHORT_LOCUS_SPAN: u32 = 100;

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
/// written. A run given no value takes `one_per_worker_thread`.
///
/// **Owed:** the resolved value belongs in the run's output beside the other three, because
/// two runs differing only in it differ in memory and in nothing a reader of the output can
/// see.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct CohortLocusBuilderRegionsInFlight(pub NonZeroUsize);

impl CohortLocusBuilderRegionsInFlight {
    /// One region at a time — the streaming driver's arrangement, and every default's.
    pub const ONE: Self = Self(NonZeroUsize::MIN);

    #[cfg(test)]
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cohort_observation::{
        DEFAULT_MIN_ALT_OBS, DEFAULT_MIN_ALT_READ_SHARE, MinAltObs, MinAltReadShare, MinAltReads,
    };

    /// Both spellings of each default — the constant a command line will advertise and
    /// the typed value a run will use — pinned to the documented number, so the two
    /// cannot drift into printing one and running the other.
    #[test]
    fn the_defaults_are_the_documented_values() {
        assert_eq!(DEFAULT_MAX_COHORT_LOCUS_SPAN, 100);
        assert_eq!(DEFAULT_MIN_ALT_OBS, 2);
        assert_eq!(DEFAULT_MIN_ALT_READ_SHARE, 0.02);
        assert_eq!(DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN, 500);

        assert_eq!(MaxCohortLocusSpan::default().get(), 100);
        assert_eq!(MinAltObs::default().get(), 2);
        assert_eq!(MinAltReadShare::default().get(), 0.02);
        assert_eq!(CohortLocusBuilderRegionsLen::default().get(), 500);

        assert_eq!(MinAltReads::default(), MinAltReads::DEFAULT);
        assert_eq!(MinAltReads::DEFAULT.floor.get(), 2);
        assert_eq!(MinAltReads::DEFAULT.share.get(), 0.02);
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
