//! ng's generic locus generator — the pileup walk over a non-STR region.
//!
//! A **folder**, unlike its `ssr` sibling's single entry file, because it holds a whole walker
//! plus the generator that wraps it (`doc/devel/ng/spec/locus_generation_pileup.md`,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*).
//!
//! # Where the walker came from
//!
//! [`genome_walk`], [`open_record`], [`cigar_cursor`], [`decompose`], [`active_read_set`],
//! [`chain_id_allocator`] and [`errors`] — plus `tests.rs`, that walker's own end-to-end
//! suite — were transcribed unchanged from the pileup walker of the caller deleted on
//! 2026-09-13 (`src/pileup/walker/` at commit `d9e7b076`). The rule that paid on this branch
//! was *transcribe first, change second*: a copy that provably matched the old walker was the
//! baseline every later change was measured against, so the generator's deliberate
//! divergences could be told from transcription slips. Every one of those files has changed
//! since, and each says so in its own header. Only one was renamed on the way in: `driver.rs`
//! became [`genome_walk`], which names the axis the walk advances along; its walker type keeps
//! the old name, `PileupWalker`.
//!
//! [`generator`], [`fast_column`], [`witnessed_ref`] and [`read_sampling`] are ng's own.
//! `parity.rs` is what is left of the comparison against the old walker: determinism across
//! processes, the generated cases' coverage floor, and the old walker's error stream, frozen.
//!
//! The walk emits ng's own [`SampleLocusObservations`](crate::locus_generation::SampleLocusObservations).
//! The inherited tests assert on that type through two accessors that **panic** where the old
//! walker's positional idiom has no ng answer (`tests::Locus::reference_observation`,
//! `first_alt_observation`), so a test landing on one of those cases has to be rewritten
//! rather than quietly reinterpreted.
//!
//! # What this module re-exports, and why
//!
//! The transcribed files reach their shared vocabulary through `super::`, as they did in the
//! old walker. This module therefore declares or re-exports each name they reach that way:
//!
//! - [`PreparedRead`], [`MateRole`], [`ReadLengthError`] — **ng's**, extended with
//!   `read_group` (spec §6). Every transcribed file names `PreparedRead` in its signatures.
//! - [`CigarOp`] — **the decoder's** ([`crate::bam::alignment_input`]): the input stage both
//!   builds a CIGAR run and reads it back, so it owns the type. Re-exported here only so the
//!   transcribed files' `super::CigarOp` resolves.
//! - [`WalkerConfig`] and the `DEFAULT_*` limits — **ng's own**, declared below. All but one
//!   keep the old walker's values. **[`DEFAULT_MAX_ACTIVE_READS`] is the one that diverged**,
//!   and it is declared in `chain_id_allocator.rs` because the transcription put it there: ng
//!   holds 32,768 reads open against the old walker's 4,096, and
//!   `the_copied_active_reads_cap_is_still_productions` records the pair rather than pinning
//!   them equal.
//!
//! The vocabulary is bound `pub(crate)`, not `pub`: it is an internal aid for the
//! transcribed files, not a second public path to a name declared elsewhere. A consumer
//! outside this module that wants `CigarOp` should say
//! [`crate::bam::alignment_input::CigarOp`] and name the decoder that produced it.
//!
//! # The reference: ng's own, with no adaptor between
//!
//! The walk fetches through [`RefSeq`](crate::ref_seq::RefSeq) directly. It used to go
//! through a newtype that presented ng's reference as the old walker's
//! `MultiChromRefFetcher`; that newtype, its error translation and both of that translation's
//! lossy spots (a contig *name* rendered as an id; a `u64 → u32` narrowing) are **deleted**,
//! and no file in this module imports `MultiChromRefFetcher` or `ChromRefFetchError` outside
//! `#[cfg(test)]` code.

mod active_read_set;
mod chain_id_allocator;
mod cigar_cursor;
mod decompose;
mod errors;
mod fast_column;
mod generator;
mod genome_walk;
mod open_record;
/// **ng's, not a copy** — the deterministic per-read number both depth caps select on,
/// so which reads survive a cap is a fact about the reads and not about the container
/// that happens to hold them.
mod read_sampling;
mod witnessed_ref;

/// The deleted caller's end-to-end walker suite, transcribed and run against ng's walk
/// (spec §12). Its `MockFasta` / `snp_read` / `paired_snp_reads` fixtures are what
/// `parity.rs` builds on.
#[cfg(test)]
pub(crate) mod tests;

/// **Two walks emitted the same evidence, up to what their chain ids are called** — the
/// comparison the region-tiling tests make since the owner's ruling of 2026-08-17 put an id
/// on every observation.
///
/// **A chain id is allocated in admission order, so it is a label and not a name.** Two
/// walks over the same ground admit the same reads in the same order *relative to each
/// other* but not from the same starting number: splitting one region into two adjacent
/// ones re-admits the reads that straddle the join, so everything after it shifts. Before
/// the ruling the shift was invisible, because the reads that carry the loci over most of
/// the genome agree with the reference and carried no id at all; now it shows in every
/// observation, and byte-equality between two walks is the wrong test.
///
/// So this asserts the property that is actually meant: the evidence is identical, and
/// **one consistent renaming** carries the first walk's ids onto the second's across
/// everything it is given. A walk that merged two reads into one identity, split one into
/// two, or renumbered a read partway through fails it.
///
/// **Give it one region at a time, because a renaming survives exactly one region.** A read
/// that straddles a join is admitted in *both* regions of a split walk and allocated an id in
/// each: on this module's tiling fixture the read that the whole-region walk calls 1
/// throughout is 1 before position 50 and **3** after it. One read, two identities, and no
/// way for the walk to know otherwise — it meets the read twice. Handing a whole split walk
/// to one call would therefore fail; handing it region by region is what asserts that the
/// renumbering happens **only** at the join, which is the property that matters. A map
/// rebuilt per locus instead would accept a walk that renumbered its reads mid-region —
/// measured: adding 100 to every id from position 60 on passed under that weaker form.
///
/// **That is safe for everything downstream, and only because of where walks are cut.** The
/// cohort merge compares ids within one sample's stream inside one cohort locus; a locus
/// never crosses a segment boundary, and a segment is never cut
/// (`doc/devel/ng/spec/run_streaming.md` §4.3, the owner's rule of 2026-08-09), so the two
/// identities of a straddling read always fall in different loci. A run that ever cut inside
/// a segment would break the merge's read-linking silently, and this is the test that would
/// stop saying so.
#[cfg(test)]
fn assert_same_evidence_up_to_chain_renaming(
    left: &[crate::locus_generation::SampleLocusObservations],
    right: &[crate::locus_generation::SampleLocusObservations],
    what: &str,
) {
    use crate::locus_generation::SampleLocusObservations;
    use crate::types::ChainId;
    use std::collections::HashMap;

    let without_ids = |locus: &SampleLocusObservations| {
        let mut stripped = locus.clone();
        // The region's first id goes with the ids: two walked regions are two numberings.
        stripped.chain_id_base = 0;
        for observation in &mut stripped.observations {
            observation.chain_ids.clear();
        }
        stripped
    };
    let left_without_ids: Vec<_> = left.iter().map(without_ids).collect();
    let right_without_ids: Vec<_> = right.iter().map(without_ids).collect();
    assert_eq!(
        left_without_ids, right_without_ids,
        "{what}: the evidence differs, chain ids aside",
    );

    // One map for the whole call, which is why the caller passes one region at a time.
    let mut onto: HashMap<ChainId, ChainId> = HashMap::new();
    let mut back: HashMap<ChainId, ChainId> = HashMap::new();
    for (left_locus, right_locus) in left.iter().zip(right) {
        for (left_observation, right_observation) in left_locus
            .observations
            .iter()
            .zip(&right_locus.observations)
        {
            assert_eq!(
                left_observation.chain_ids.len(),
                right_observation.chain_ids.len(),
                "{what}: {} observations name a different number of reads",
                left_locus.region,
            );
            for (&left_id, &right_id) in left_observation
                .chain_ids
                .iter()
                .zip(&right_observation.chain_ids)
            {
                assert_eq!(
                    *onto.entry(left_id).or_insert(right_id),
                    right_id,
                    "{what}: at {} chain id {left_id} stands for two different reads",
                    left_locus.region,
                );
                assert_eq!(
                    *back.entry(right_id).or_insert(left_id),
                    left_id,
                    "{what}: at {} two reads collapsed onto chain id {right_id}",
                    left_locus.region,
                );
            }
        }
    }
}

/// The comparison above is a test's only guard against a walk that renumbers its reads, so
/// it is itself tested: **it has to reject a renaming that changes partway through.**
///
/// Written because the first version rebuilt its map per locus, which accepts exactly that —
/// measured on the tiling fixture, where adding 100 to every id from position 60 on passed.
#[cfg(test)]
mod chain_renaming_tests {
    use crate::locus_generation::{
        LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
    };
    use crate::types::ChainId;
    use crate::types::{ContigId, GenomeRegion, Position, ReadGroupId};

    /// One position covered by one read, named `chain`.
    fn locus_named(position: u64, chain: ChainId) -> SampleLocusObservations {
        SampleLocusObservations {
            chain_id_base: 0,
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(position),
                end: Position(position),
            },
            reference_bases: Box::from(&b"A"[..]),
            observations: vec![SequenceObservation {
                flank_deletion: crate::locus_generation::FlankDeletion::NONE,
                held_flank: None,
                bases: Box::from(&b"A"[..]),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 1,
                num_fwd: 1,
                q_sum: crate::types::SummedLogError::from_nats(0.0),
                mapq_sum: 60,
                mapq_sum_sq: 3600,
                placed_left: 0,
                num_low_quality: 0,
                chain_ids: vec![chain],
            }],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// One read across two positions, called 1 by one walk and 7 by the other, is the same
    /// read consistently — which is what the two walks of the tiling tests actually differ by.
    #[test]
    fn a_renaming_that_holds_throughout_is_accepted() {
        let one_walk = [locus_named(10, 1), locus_named(11, 1)];
        let another = [locus_named(10, 7), locus_named(11, 7)];

        super::assert_same_evidence_up_to_chain_renaming(&one_walk, &another, "consistent");
    }

    /// **The same read renumbered at the second position is refused.** That is one read
    /// becoming two identities inside a stretch the caller said was one region, which is the
    /// state the cohort merge cannot work in — and the per-locus map this helper started with
    /// accepted it.
    #[test]
    #[should_panic(expected = "stands for two different reads")]
    fn a_renaming_that_changes_partway_through_is_refused() {
        let one_walk = [locus_named(10, 1), locus_named(11, 1)];
        let renumbered = [locus_named(10, 7), locus_named(11, 100)];

        super::assert_same_evidence_up_to_chain_renaming(&one_walk, &renumbered, "renumbered");
    }

    /// And two reads collapsing onto one identity is refused from the other side.
    #[test]
    #[should_panic(expected = "two reads collapsed onto chain id")]
    fn two_reads_given_one_identity_are_refused() {
        let one_walk = [locus_named(10, 1), locus_named(11, 2)];
        let collapsed = [locus_named(10, 7), locus_named(11, 7)];

        super::assert_same_evidence_up_to_chain_renaming(&one_walk, &collapsed, "collapsed");
    }
}

/// **ng's, not a copy** — [`RefSeq`](crate::ref_seq::RefSeq) over `tests.rs`'s
/// `MockFasta`, so production's copied suite can drive ng's walker (A0).
#[cfg(test)]
mod mock_reference;

/// **ng's, not a copy** — what is left of the stage-1 differential: determinism across
/// processes, the generated cases' coverage floor, and production's error stream
/// frozen (spec §3, §13.1).
#[cfg(test)]
mod parity;

// The vocabulary the transcribed files resolve through `super::` — the names the old
// walker's `walker/mod.rs` declared, bar its `indel_norm` module, which none of the seven
// reached.
pub(crate) use crate::bam::alignment_input::CigarOp;
pub(crate) use crate::read::prepared_read::{MateRole, PreparedRead, ReadLengthError};

// ---------------------------------------------------------------------
// What the walk is allowed to hold, and how far it will look
// ---------------------------------------------------------------------

/// How wide a single record may get before the walk refuses it, in reference bases.
///
/// Bounds what one open record costs. Reads wider than this are dropped by the read filters
/// upstream, so in an ordinary run the walk's own check never fires; it is here because a
/// record that grew past it would otherwise grow without limit.
pub(crate) const DEFAULT_MAX_RECORD_SPAN: u32 = 5_000;

/// How far past a first mate's start the walk keeps looking for its partner, in reference
/// bases, before giving up and treating the read as unpaired.
///
/// Sized for ordinary Illumina paired-end inserts of 200 to 500 bases with room to spare. A
/// protocol whose mates sit further apart than this would need it raised, and would otherwise
/// lose the pairing quietly rather than loudly.
pub(crate) const DEFAULT_MATE_LOOKUP_WINDOW: u32 = 10_000;

/// How many reads **of one read group** the walk will fold into one position — the default of
/// `--max-reads-per-position`.
///
/// **1,000 since 2026-09-29; it was samtools' `MPLP_MAX_DEPTH`, 8,000, per sample.** At 8,000
/// the cap never fired at the collapsed-repeat pile-up that stopped the 2,169-sample tomato run
/// (a median of 5,669 reads a position in each 3× sample), so every one of those reads reached
/// the stored files and the calling that read them. 1,000 is the owner's choice: far above any
/// position's honest depth at the coverages the caller is for, low enough to bound a pile-up.
pub(crate) const DEFAULT_MAX_READS_PER_POSITION: u32 = 1_000;

// **There is no second, tighter cap at a position where a read has an insertion or deletion.**
// There was one until 2026-10-03 — 250 reads, samtools' `MPLP_MAX_INDEL_DEPTH` — and it biased
// every indel the walk recorded at depth. An indel's record spans more than one position: the
// deletion's reads are folded at the position the deletion is anchored on, where the tighter
// cap applied, but the reference reads are also folded at the positions the deletion covers,
// where it did not. So above 250 reads the deletion lost reads and the reference kept nearly
// all of its own. At 300× a true heterozygous deletion in 47% of fragments was stored at 31%,
// and the hidden-duplication filter then removed it as a collapsed repeat
// (`doc/devel/ng/research/giab_unexplained_fn_2026-10-03.md` §3). The tighter cap saved no
// measurable work there either: the psps came out 0.2% smaller without it and the calling took
// the same time.

/// The four limits the walk reads as it runs.
///
/// **The truncation is per position and per read group, not per allele**, and that is the
/// property worth stating: when one read group has more reads at a position than the cap
/// allows, the walk keeps `cap` of them and drops the rest. Clipping each allele separately
/// instead would bias the allele frequency — a position where 99 reads show one base and 1
/// shows another would come back as roughly 71 against 29 at a cap of 250.
///
/// **Which reads are kept is a function of the reads**: the smallest sampling keys, a hash of
/// the query name ([`read_sampling`]), so the same reads survive at neighbouring positions.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub(crate) struct WalkerConfig {
    /// Reads of one read group folded at one position, whatever the reads show there. See
    /// [`DEFAULT_MAX_READS_PER_POSITION`].
    pub max_reads_per_position: u32,
    /// How wide one record may get. See [`DEFAULT_MAX_RECORD_SPAN`].
    pub max_record_span: u32,
    /// How far to look for a mate. See [`DEFAULT_MATE_LOOKUP_WINDOW`].
    pub mate_lookup_window: u32,
    /// How many reads may be open at once; exceeding it stops the walk rather than growing
    /// memory without limit. See [`DEFAULT_MAX_ACTIVE_READS`].
    pub max_active_reads: u32,
}

impl Default for WalkerConfig {
    /// **`max_active_reads` is ng's 32,768 and not production's 4,096**, which is the one
    /// limit of the five ng deliberately moved (`chain_id_allocator`'s own header says why).
    /// Nothing a run executes reads this `Default` — [`PileupGeneratorConfig`] carries the
    /// limits a run walks under and builds its own — so it is the fixtures that see it.
    fn default() -> Self {
        Self {
            max_reads_per_position: DEFAULT_MAX_READS_PER_POSITION,
            max_record_span: DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: DEFAULT_MATE_LOOKUP_WINDOW,
            max_active_reads: DEFAULT_MAX_ACTIVE_READS,
        }
    }
}

// This module's own surface, as against the vocabulary above.
//
// **Kept to what a consumer outside the module needs.** The walk itself — `run`,
// `PileupWalker`, `RunSummary` — stays private, and `run` mattered most:
// `to_walker_config` is `pub(super)` precisely so every caller passes through
// `PileupGeneratorConfig::check()` and its `MAX_RECORD_SPAN_CEILING`. With `run` public, an
// external caller could build a `WalkerConfig` from its `Default` and walk straight past
// that check; measured from an external example when `run` was still public:
// `max_record_span = 983_025`, **15× the ceiling `check()` enforces**.
//
// The surface was audited by demoting each name in turn and rebuilding with
// `cargo check --all-targets --all-features`, **without** `-D warnings`: demoting
// `PileupGenerator` makes the lib target treat this module as dead, and that cascade stops
// cargo ever reaching the code that consumes it.
pub use errors::WalkerError;
pub use generator::{
    MAX_RECORD_SPAN_CEILING, PileupGenerator, PileupGeneratorConfig, PileupGeneratorConfigError,
    PileupGeneratorCounts,
};
/// **How wrong one read is, as the walk mints it** — the worse of its base quality and its
/// mapping quality, in log space.
///
/// Re-exported because spec §3.2 makes it a requirement rather than a convenience: the
/// quantity the calibration accumulator averages and the quantity the read likelihood charges
/// must be computed by *the same function*, or the scale calibrates against a different
/// definition of "how wrong is this read" than the one it is applied to. The likelihood's
/// tests mint their fixtures through this rather than re-deriving `10^(−q/10)`, which differs
/// from the crate's own table in the last place or two and would be a second definition of
/// the quantity all the same.
///
/// `#[cfg(test)]` because that is where the consumer is, and for the reason
/// [`DEFAULT_MATE_LOOKUP_WINDOW`] above carries: left ungated it is an unused import in a
/// non-test build. **The shipping path never needs it** — a row charges the `q_sum` the walk
/// already summed, and minting is the walk's own business.
#[cfg(test)]
pub(crate) use open_record::minted_ln_read_error;

// **Not `pub`: no consumer outside this module, and it is production's name.** Kept in
// scope rather than deleted because `generator.rs` reaches it as
// `super::DEFAULT_MAX_ACTIVE_READS` and this module's own doc links it — an intra-doc
// link resolves against private imports, so demoting costs nothing while deleting would
// make `cargo doc` gain a thirteenth error.
use chain_id_allocator::DEFAULT_MAX_ACTIVE_READS;

// **Demoted, not deleted, and `#[cfg(test)]` because every consumer is a test one.**
// `run` and `RunSummary` are reached as `super::run` / `super::RunSummary` by four
// `#[cfg(test)]` descendants (`tests.rs`, `open_record.rs`'s tests, `mock_reference.rs`,
// `parity.rs`) — a private `use` still serves them, because privacy in Rust reaches
// descendants. Non-test code does not go through here at all: `generator.rs` names
// `super::genome_walk::{PileupWalker, RunSummary}` directly, which is why `PileupWalker`
// is gone from this list rather than demoted into it.
//
// Ungated it would be an unused import in a non-test build — the same trap the
// `DEFAULT_MATE_LOOKUP_WINDOW` note above records, and the reason it went unnoticed for
// both: **a `pub` re-export is never reported unused**, so `pub` was hiding the fact that
// nothing outside wanted these.
#[cfg(test)]
use genome_walk::{RunSummary, run};

// `walker_vocabulary_tests`, not `tests`: the copied `walker/tests.rs` lands as this
// module's `tests` child (A4), and mirroring production's module names is what makes
// the two suites comparable name for name. Named for its subject rather than for the
// pattern, following `mod baq_tests`.
#[cfg(test)]
mod walker_vocabulary_tests {
    use super::*;

    /// The one `DEFAULT_*` the copy had to bring with it — **and the first one ng has
    /// deliberately moved.**
    ///
    /// `DEFAULT_MAX_ACTIVE_READS` is declared *inside* `chain_id_allocator.rs`, so the
    /// verbatim copy forked it, where the other constants are reached from production by
    /// name and cannot drift. This test used to pin the two equal, and its own doc said
    /// what to do when they were allowed to differ: *"when they are deliberately allowed
    /// to differ, that test is what says so."* This is that.
    ///
    /// ng holds **eight times** as many reads as production. Production's 4,096 was
    /// refusing 19,725 reads of 113,629,764 on one ~130× tomato chromosome, and a read
    /// refused at the door contributes at no position — so positions ended up with less
    /// coverage than the input had for them (owner, 2026-08-05). Both numbers were asserted
    /// until promotion step C18; production's is now only this comment, since there is no
    /// production constant left to retune.
    #[test]
    fn the_copied_active_reads_cap_is_still_productions() {
        assert_eq!(
            DEFAULT_MAX_ACTIVE_READS, 32_768,
            "ng's ceiling on reads held open at once"
        );
    }
}
