//! ng's generic locus generator — the pileup walk, copied from production.
//!
//! A **folder**, unlike its `ssr.rs` sibling, because it holds the whole copied
//! walker plus (later) the generator that wraps it. Production is not edited:
//! ng copies rather than reaches in, so there is no visibility lift and no field
//! on a frozen type (`doc/devel/ng/spec/locus_generation_pileup.md` §3,
//! `doc/devel/ng/arch/locus_generation_pileup.md` *Module home*).
//!
//! # Eight files were transcribed; five are still verbatim
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
//! [`open_record`] and [`errors`] — the reference adaptor's removal; the other
//! five are still guarded.
//!
//! What the copy was proven to *compute* is the stage-1 differential
//! (`parity.rs`); what it is proven to *be* is `copy_fidelity.rs`. The two are
//! different claims. Both are named as files rather than linked: they are
//! `#[cfg(test)]` modules, so an intra-doc link to them breaks `cargo doc`.
//!
//! The copies still emit production's [`PileupRecord`](crate::pileup_record::PileupRecord),
//! not ng's `SampleLocusObservations`. That changes in Milestone B.
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
pub use chain_id_allocator::DEFAULT_MAX_ACTIVE_READS;
pub use errors::WalkerError;
pub use genome_walk::{PileupWalker, RunSummary, run};

/// Lay a finished locus back out as production's [`PileupRecord`] — **test-only, and the
/// regression floor depends on it being read as what it is.**
///
/// B2 changed what the walk emits. The 44 inherited tests that came with the copy are the
/// walk's regression floor, and spec §12 sanctions adapting them mechanically to the new
/// type. Routing them through one projection *is* that adaptation, and it is the safer form
/// of it: 67 hand-translated assertions are 67 chances to re-express a test slightly weaker
/// than it was, which is the failure mode this branch has hit in six consecutive
/// milestones. One function can be reviewed once.
///
/// # What it cannot reproduce, and why each is a deliberate change rather than a gap
///
/// - **`placed_start` comes back as `0`.** ng does not carry it at all from B2 on (spec §6);
///   no model consumes it, and it is a pure function of the read's start against the anchor,
///   so a later consumer re-derives it without touching the fold. The three inherited
///   assertions on it change, which spec §12 already lists.
/// - **Allele order is ng's**, which is sorted by `(bases, read_coverage, read_group)`,
///   where production's is bucket-creation order. Spec §12 lists order-asserting tests as
///   must-change. REF is placed first regardless, because `alleles[0] == REF` is an
///   invariant of the type being projected onto, not an ordering preference.
/// - **Chain ids follow ng's per-read rule**, so a read that agreed with the reference over
///   everything it witnessed carries none even where production would have given it one.
///
/// Rows that split by coverage or read group are **merged back together** here. That is the
/// point: the inherited suite tests the *walk* — which reads folded where, with what
/// evidence — not the shape of the emitted type, which has its own tests (`observation_rows`)
/// and the differential.
#[cfg(test)]
pub(crate) fn as_pileup_record(
    locus: &crate::ng::locus_generation::SampleLocusObservations,
) -> crate::pileup_record::PileupRecord {
    use crate::pileup_record::{AlleleObservation, AlleleSupportStats};

    let mut alleles: Vec<AlleleObservation> = Vec::new();
    // REF first and always present: production creates `alleles[0]` at record open
    // regardless of support, and ng emits no row for a reference nobody witnessed.
    alleles.push(AlleleObservation {
        seq: locus.reference_bases.to_vec(),
        support: AlleleSupportStats::default(),
        chain_ids: Vec::new(),
    });
    for observation in &locus.observed_sequences {
        let index = match alleles
            .iter()
            .position(|allele| allele.seq.as_slice() == observation.bases.as_ref())
        {
            Some(index) => index,
            None => {
                alleles.push(AlleleObservation {
                    seq: observation.bases.to_vec(),
                    support: AlleleSupportStats::default(),
                    chain_ids: Vec::new(),
                });
                alleles.len() - 1
            }
        };
        let support = &mut alleles[index].support;
        support.num_obs += observation.num_obs;
        support.q_sum += observation.q_sum;
        support.fwd += observation.num_fwd;
        support.placed_left += observation.placed_left;
        support.mapq_sum += observation.mapq_sum;
        support.mapq_sum_sq += observation.mapq_sum_sq;
        alleles[index]
            .chain_ids
            .extend_from_slice(&observation.chain_ids);
    }
    for allele in &mut alleles {
        allele.chain_ids.sort_unstable();
        allele.chain_ids.dedup();
    }
    crate::pileup_record::PileupRecord {
        chrom_id: locus.region.contig.0,
        pos: locus.region.start.get() as u32,
        alleles,
        // The walker cannot compute the centred window (it needs look-ahead); the
        // pileup→psp seam fills these from the sliding-window accumulator before writing.
        windowed_gc: f32::NAN,
        windowed_coverage: f32::NAN,
    }
}

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
