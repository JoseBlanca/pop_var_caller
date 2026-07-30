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
//! production's positional idiom has no ng answer (`tests::Locus::reference_row`,
//! `first_alt_row`), so a test landing on one of those cases has to be rewritten rather
//! than quietly reinterpreted. `to_pileup_record` is deleted.
//!
//! **One file is renamed on the way in: `driver.rs` → [`genome_walk`]** — it is
//! the only one of the seven named for a *role* rather than for what it owns,
//! and "driver" answers *driver of what?* with nothing. `genome_walk` names the
//! one job that file has and the others do not: advancing a position cursor
//! along genome coordinates over an active read set. The **type** keeps
//! production's `PileupWalker`, so the differential reads as a straight
//! comparison. *(A walk covers one region, not the genome; `genome_walk` names
//! the axis it advances along, not the extent.)*
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
mod generator;
mod genome_walk;
mod open_record;

/// Production's own end-to-end suite, copied verbatim and run against the copy —
/// Milestone A's gate (spec §12). `pub(crate)` because production's declaration is,
/// so the two module trees stay comparable; its `MockFasta` / `snp_read` /
/// `paired_snp_reads` fixtures are reachable from B1's parity harness either way,
/// that being a descendant of this module.
#[cfg(test)]
pub(crate) mod tests;

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

    /// The one `DEFAULT_*` the copy had to bring with it, pinned against production's.
    ///
    /// `DEFAULT_MAX_ACTIVE_READS` is declared *inside* `chain_id_allocator.rs`, so the
    /// verbatim copy forked it — where the other constants are reached from production
    /// by name and cannot drift. Production's value is the baseline the stage-1
    /// differential runs against, so a retune on either side must fail here rather than
    /// silently in the walk.
    #[test]
    fn the_copied_active_reads_cap_is_still_productions() {
        assert_eq!(
            DEFAULT_MAX_ACTIVE_READS,
            crate::pileup::walker::DEFAULT_MAX_ACTIVE_READS,
            "ng's copied chain_id_allocator has drifted from production's cap"
        );
    }
}
