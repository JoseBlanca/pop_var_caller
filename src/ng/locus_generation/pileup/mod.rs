//! ng's generic locus generator — the pileup walk, copied from production.
//!
//! A **folder**, unlike its `ssr.rs` sibling, because it holds the whole copied
//! walker plus (later) the generator that wraps it. Production is not edited:
//! ng copies rather than reaches in, so there is no visibility lift and no field
//! on a frozen type (`doc/devel/ng/spec/locus_generation_pileup.md` §3,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*).
//!
//! # Eight files were transcribed; three are still verbatim
//!
//! [`genome_walk`], [`open_record`], [`cigar_cursor`], [`decompose`],
//! [`active_read_set`], [`chain_id_allocator`] and [`errors`] — plus `tests.rs`,
//! production's own end-to-end suite — were transcribed from
//! `src/pileup/walker/` **unchanged**, and `copy_fidelity.rs` checks that
//! textually rather than leaving it a claim in this comment. The rule that paid
//! three times on this branch is *transcribe first, change second*: a copy that is
//! provably production is the baseline every later change is measured against, and
//! without it the generator's deliberate divergences could not be told from
//! transcription slips.
//!
//! **Plan 3 spends that baseline, file by file.** Each step releases the file it
//! changes from `copy_fidelity.rs`'s checked set and says so in that file's own
//! header, so the guard keeps protecting what is still a copy instead of being
//! switched off wholesale at the first change. A0 released [`genome_walk`],
//! [`open_record`] and [`errors`] — the reference adaptor's removal — and B2 released
//! `tests.rs`, whose assertions had to move to ng's own locus type (spec §12), and D2
//! released [`active_read_set`] for the per-read "ever contributed" flag. **Three are
//! still guarded**, and `copy_fidelity.rs`'s release table is the list that stays true.
//!
//! What the copy was proven to *compute* is the stage-1 differential
//! (`parity.rs`); what it is proven to *be* is `copy_fidelity.rs`. The two are
//! different claims. Both are named as files rather than linked: they are
//! `#[cfg(test)]` modules, so an intra-doc link to them breaks `cargo doc`.
//!
//! **The walk emits ng's own [`SampleLocusObservations`](crate::ng::locus_generation::SampleLocusObservations)
//! from B2 on**, not production's `PileupRecord` — and since Milestone D **nothing in this
//! module sees it through anything else.** B2 adapted the 44 inherited tests through
//! `to_pileup_record`, one reviewed back-projection rather than 67 hand-edited assertions
//! (spec §12); that projection merged the observations ng splits and dropped the three fields
//! production has no counterpart for, and at Milestone B its losses hid three live surfaces
//! from the review. D1 removed the differential's need for it by projecting production
//! *forward* instead, and the same commit that answers Checkpoint D removed the suite's:
//! the inherited tests now assert on ng's type through two accessors that **panic** where
//! production's positional idiom has no ng answer (`tests::Locus::reference_observation`,
//! `first_alt_observation`), so a test landing on one of those cases has to be rewritten rather
//! than quietly reinterpreted. `to_pileup_record` is deleted.
//!
//! **One file is renamed on the way in: `driver.rs` → [`genome_walk`]** — it is
//! the only one of the seven named for a *role* rather than for what it owns,
//! and "driver" answers *driver of what?* with nothing. `genome_walk` names the
//! one job that file has and the others do not: advancing a position cursor
//! along genome coordinates over an active read set. The **type** keeps
//! production's `PileupWalker`, so the differential reads as a straight
//! comparison. *(`genome_walk` names the axis it advances along, not the extent:
//! since D1 a walker is pointed at one region after another and lives for a
//! chromosome, and before that it covered a single region — the name was right
//! for both.)*
//!
//! # What this module re-exports, and why
//!
//! The copied files reach their shared vocabulary through `super::` — that is
//! how they were written against `pileup/walker/mod.rs`, and leaving those paths
//! alone is what keeps the transcription verbatim. This module therefore stands
//! in for production's `walker/mod.rs`, drawing each name from wherever it now
//! lives:
//!
//! - [`PreparedRead`], [`MateRole`], [`ReadLengthError`] — **ng's**, copied and
//!   extended with `read_group` (spec §6). They are the reason the whole walker
//!   is copied: every one of the seven names `PreparedRead` in its signatures.
//! - [`CigarOp`], [`WalkerConfig`] and two of the `DEFAULT_*` constants —
//!   **production's, reused as-is**. ng does not modify them, so it does not copy
//!   them; they are reached by name rather than by literal so there is one source
//!   of truth until ng deliberately diverges.
//!   **One constant is ng's, forced by the verbatim rule:**
//!   [`DEFAULT_MAX_ACTIVE_READS`] is declared *inside* `chain_id_allocator.rs`, so
//!   the copy brought its own — two definitions of `4096` now exist, and
//!   `chain_id_allocator.rs`'s `Self::with_caps(DEFAULT_MAX_ACTIVE_READS,
//!   super::DEFAULT_MATE_LOOKUP_WINDOW)` reaches ng's for the first and production's
//!   for the second. `the_copied_active_reads_cap_is_still_productions` pins the two
//!   equal; when they are deliberately allowed to differ, that test is what says so.
//!
//! The vocabulary is bound `pub(crate)`, not `pub`: it is an internal aid for the
//! copies, not an ng-flavoured public alias for production's types. A consumer that
//! wants `CigarOp` should say `crate::pileup::walker::CigarOp` and see whose it is —
//! which matters from plan 3 on, when ng's walker starts to diverge and two live
//! paths to one name would stop being harmless.
//!
//! # The reference: ng's own, with no adaptor between
//!
//! The walk fetches through [`RefSeq`](crate::ng::ref_seq::RefSeq) directly
//! (A0). It used to go through a `RefSeqFetcher` newtype presenting ng's
//! reference as production's `MultiChromRefFetcher` — a consequence of the
//! transcription being verbatim, never a design choice, and one that stopped
//! being true the moment the two walkers began to diverge. The newtype, its
//! error translation and both of that translation's lossy spots (a contig *name*
//! rendered as an id; a `u64 → u32` narrowing) are **deleted**, and **no file in
//! this module** imports `MultiChromRefFetcher` or `ChromRefFetchError` outside
//! `#[cfg(test)]` code. Scoped to this module deliberately: an earlier wording
//! claimed it of ng's non-test code as a whole, which is false —
//! [`raw_chrom_reader`](crate::ng::raw_chrom_reader) names `MultiChromRefFetcher`
//! in public signatures, and always did.

mod active_read_set;
mod chain_id_allocator;
mod cigar_cursor;
mod decompose;
mod errors;
mod fast_column;
mod generator;
mod genome_walk;
// **Measurement scaffolding, not part of the walk** — per-read-group sums of the minted
// per-read error in both the shapes an average can be taken in, so the gap between the
// geometric mean the scale uses and the arithmetic one the spec first asked for can be
// measured on real reads. Off unless `PVC_MINTED_ERROR_CENSUS=1`; read by
// `ng_minted_error_means`. The module's own header is the documentation; a doc comment
// here as well would be resolved in *this* module's scope, where the names it links are
// not in scope, and `cargo doc` counts each of those as an error.
pub mod minted_error_census;
mod open_record;
/// **ng's, not a copy** — the deterministic per-read number both depth caps select on,
/// so which reads survive a cap is a fact about the reads and not about the container
/// that happens to hold them.
mod read_sampling;
mod witnessed_ref;

/// Production's own end-to-end suite, copied verbatim and run against the copy —
/// Milestone A's gate (spec §12). `pub(crate)` because production's declaration is,
/// so the two module trees stay comparable; its `MockFasta` / `snp_read` /
/// `paired_snp_reads` fixtures are reachable from B1's parity harness either way,
/// that being a descendant of this module.
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
    left: &[crate::ng::locus_generation::SampleLocusObservations],
    right: &[crate::ng::locus_generation::SampleLocusObservations],
    what: &str,
) {
    use crate::ng::locus_generation::SampleLocusObservations;
    use crate::pileup_record::ChainId;
    use std::collections::HashMap;

    let without_ids = |locus: &SampleLocusObservations| {
        let mut stripped = locus.clone();
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
    use crate::ng::locus_generation::{
        LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
    };
    use crate::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId};
    use crate::pileup_record::ChainId;

    /// One position covered by one read, named `chain`.
    fn locus_named(position: u64, chain: ChainId) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(position),
                end: Position(position),
            },
            reference_bases: Box::from(&b"A"[..]),
            observations: vec![SequenceObservation {
                bases: Box::from(&b"A"[..]),
                read_witness: ReadWitness::Complete,
                read_group: ReadGroupId(0),
                num_obs: 1,
                num_fwd: 1,
                q_sum: crate::ng::types::SummedLogError::from_nats(0.0),
                mapq_sum: 60,
                mapq_sum_sq: 3600,
                placed_left: 0,
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

/// **ng's, not a copy** — the textual check that the still-untouched copies are
/// production's, from outside the files it checks (spec §3).
#[cfg(test)]
mod copy_fidelity;

/// **ng's, not a copy** — [`RefSeq`](crate::ng::ref_seq::RefSeq) over `tests.rs`'s
/// `MockFasta`, so production's copied suite can drive ng's walker (A0).
#[cfg(test)]
mod mock_reference;

/// **ng's, not a copy** — the stage-1 differential: the two walkers compute the
/// same thing, over one read stream (spec §3, §13.1). `copy_fidelity` says the
/// copy *is* production's; this says it *does* what production's does.
#[cfg(test)]
mod parity;

// The vocabulary the copied files resolve through `super::`. Production's own
// `walker/mod.rs` declares the same names, around the same modules bar one: its
// `pub(crate) mod indel_norm` is not copied, because none of the seven reaches it —
// its consumers are production's `read_processor.rs` and ng's own `src/ng/alignment/`,
// neither of which is part of the walk.
//
// `DEFAULT_MAX_SNP_COLUMN_DEPTH` / `DEFAULT_MAX_INDEL_COLUMN_DEPTH` are deliberately
// absent: no copied file reaches them through `super::` (production's `WalkerConfig`
// holds the column caps and comes with its own defaults), and an unused re-export is
// invisible to the compiler, so carrying them would be shape without substance.
// `PileupGeneratorConfig` names them in plan 3, from production directly.
pub(crate) use crate::ng::read::prepared_read::{MateRole, PreparedRead, ReadLengthError};
pub(crate) use crate::pileup::walker::{CigarOp, DEFAULT_MAX_RECORD_SPAN, WalkerConfig};
// `#[cfg(test)]` because that is where the one copy reaching it lives:
// `ChainIdAllocator::new` is itself `#[cfg(test)]` (production code calls `with_caps`
// from `run`). Left ungated it is an unused import in a non-test build — which the
// previous `pub use` hid, since a `pub` re-export is never reported unused.
#[cfg(test)]
pub(crate) use crate::pileup::walker::DEFAULT_MATE_LOOKUP_WINDOW;

// This module's own surface, as against the vocabulary above.
//
// **Audited by demoting each name in turn and rebuilding**, not by grepping: production
// declares identically-named `WalkerError`, `PileupWalker`, `RunSummary`, `run` and
// `DEFAULT_MAX_ACTIVE_READS`, so a bare-name search answers about the wrong crate. Four
// of the ten were load-bearing (`E0603` on demotion), two have no build consumer but a
// reason to stay, and four were surplus — the reverse of Checkpoint C's guess that nine
// of ten had no consumer. The answering build is `cargo check --all-targets
// --all-features` **without** `-D warnings`: demoting `PileupGenerator` makes the lib
// target treat this module as dead, and that cascade stops cargo ever reaching the
// example that actually consumes it.
//
// **Why the four went, and why `run` mattered most.** They put production's walker
// vocabulary into ng's *public* API under production's own names — the thing the
// vocabulary note above forbids, for the reason it gives: from plan 3 on, two live paths
// to one name stop being harmless. And `run` had teeth. `to_walker_config` is
// `pub(super)` precisely so every caller passes through
// `PileupGeneratorConfig::check()` and its `MAX_RECORD_SPAN_CEILING`; with `run` public,
// an external caller builds a `WalkerConfig` from production's `pub` `Default` and walks
// straight past that. Measured from an external example before the change:
// `max_record_span = 983_025`, **15× the ceiling `check()` enforces**.
//
// `generator.rs` imports `PileupWalker` and `RunSummary` from `genome_walk` directly, and
// `parity.rs` now names `genome_walk::run` the same way, so nothing needed a replacement
// binding.
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

/// **Measurement scaffolding, not part of the walk.** Process-global tallies of how many
/// columns the walk saw and how many of them were the *ordinary* column — one covered base,
/// every contributor showing a plain reference-anchored match, nothing to reconcile. Read by
/// `ng_generic_walk_probe`; incremented from
/// [`genome_walk`](crate::ng::locus_generation::pileup) only when `PVC_COLUMN_CENSUS=1`.
pub mod column_census {
    use std::sync::atomic::{AtomicU64, Ordering};

    macro_rules! counters {
        ($($name:ident),* $(,)?) => {
            $(pub static $name: AtomicU64 = AtomicU64::new(0);)*
            /// Every counter's current value, in declaration order, with its name.
            pub fn snapshot() -> Vec<(&'static str, u64)> {
                vec![$((stringify!($name), $name.load(Ordering::Relaxed))),*]
            }
        };
    }

    /// **Always counted, unlike everything below** — how many columns the ordinary-column
    /// path answered. A fast lane that never fires and one that works look identical in a
    /// timing, so this one is not behind the env gate. One relaxed increment per ordinary
    /// column, against the ~90 reads that column would otherwise fold.
    pub static FAST_COLUMNS: AtomicU64 = AtomicU64::new(0);

    counters! {
        COLUMNS,
        COLUMNS_ORDINARY,
        CONTRIBUTORS,
        CONTRIBUTORS_ORDINARY,
        REJECT_RECORD_ALREADY_OPEN,
        REJECT_INDEL_EVENT,
        REJECT_READ_HAS_DELETION,
        REJECT_MATE_OVERLAP,
        REJECT_DEPTH_CAP,
        REJECT_MULTI_READ_GROUP,
        REJECT_READ_HAS_INDEL,
        COLUMNS_SIMPLE,
        CONTRIBUTORS_SIMPLE,
        COLUMNS_SIMPLE_WITH_MATE,
        CONTRIBUTORS_SIMPLE_WITH_MATE,
    }

    /// Whether the census is armed for this process. Read once; the walk pays a relaxed
    /// atomic load per column when it is off.
    pub fn enabled() -> bool {
        static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
        *ON.get_or_init(|| std::env::var_os("PVC_COLUMN_CENSUS").is_some_and(|v| v == "1"))
    }

    pub(super) fn add(counter: &AtomicU64, n: u64) {
        counter.fetch_add(n, Ordering::Relaxed);
    }
}

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
    /// coverage than the input had for them (owner, 2026-08-05). Both numbers are
    /// asserted, so a retune on either side lands here as a decision rather than as
    /// drift.
    #[test]
    fn the_copied_active_reads_cap_is_still_productions() {
        assert_eq!(
            DEFAULT_MAX_ACTIVE_READS, 32_768,
            "ng's ceiling on reads held open at once"
        );
        assert_eq!(
            crate::pileup::walker::DEFAULT_MAX_ACTIVE_READS,
            4096,
            "production's, which ng no longer follows"
        );
    }
}
