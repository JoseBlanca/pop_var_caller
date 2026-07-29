//! The generic locus generator itself — the knobs it takes, the run-level counts
//! it keeps, and the state it holds across segments.
//!
//! **ng's own file, not a copy.** Everything beside it in this module was
//! transcribed from `src/pileup/walker/`; this is the part production has no
//! counterpart for, because production drives its walker from a per-chromosome
//! pipeline stage rather than from a region-at-a-time generator
//! (`doc/devel/ng/arch/locus_generation_pileup.md` §1.1, §2.1).
//!
//! C1 lands the types and the one invariant that cannot be inherited: the
//! ceiling on [`PileupGeneratorConfig::max_record_span`]. The walk over a region
//! ([`begin_segment`](super::super::LocusGenerator::begin_segment) and the halo)
//! is C2's, and the allocator's lifetime across segments is C3's.

use crate::ng::read::ReadPreparer;
use crate::ng::ref_seq::RefSeq;

use super::WalkerConfig;
use super::chain_id_allocator::ChainIdAllocator;

/// The widest `max_record_span` this generator accepts: 65,535 reference
/// positions, the widest footprint a [`ReadCoverage`] run can describe.
///
/// [`ReadCoverage::Observed`] carries `offset_in_locus` and `positions_covered`
/// as `u16`s, minted through [`LocusLen::from_positions`], which **saturates**
/// rather than failing. A record footprint wider than this therefore makes a
/// partial witness report a *truncated* `positions_covered` — a wrong number, no
/// error, at exactly the long-deletion loci this generator exists to get right.
///
/// [`ReadCoverage`]: crate::ng::locus_generation::ReadCoverage
/// [`ReadCoverage::Observed`]: crate::ng::locus_generation::ReadCoverage::Observed
/// [`LocusLen::from_positions`]: crate::ng::locus_generation::LocusLen::from_positions
pub const MAX_RECORD_SPAN_CEILING: u32 = u16::MAX as u32;

/// **Production's default `max_record_span` is inside ng's ceiling**, which is
/// why [`PileupGeneratorConfig::default`] can inherit the knob by name without
/// shipping a configuration [`PileupGenerator::new`] would reject.
///
/// A compile-time check rather than a test, because it compares two constants:
/// should production ever raise its default past 65,535, this breaks the build
/// and the divergence becomes a decision rather than a runtime surprise.
const _: () = assert!(
    crate::pileup::walker::DEFAULT_MAX_RECORD_SPAN <= MAX_RECORD_SPAN_CEILING,
    "production's default max_record_span no longer fits a ReadCoverage run: ng must either \
     widen the run or stop inheriting the default",
);

/// This generator's knobs — owned, taken at construction, and **production's
/// values, inherited and never measured by ng** (spec §7).
///
/// Raw `u32` rather than `Bp`, because the copied walker speaks production's
/// integer widths and a port must not change types under itself (spec §3). The
/// five are exactly production's [`WalkerConfig`] fields, reached from
/// production's own `pub const`s **by name** in [`Default`], so there is one
/// source of truth until ng deliberately diverges and the divergence shows up as
/// a diff.
///
/// # One knob is not simply production's
///
/// [`max_record_span`](Self::max_record_span) is capped at
/// [`MAX_RECORD_SPAN_CEILING`] and [`check`](Self::check) rejects more.
/// Production's `--max-record-span` is an unbounded `u32`; inheriting the knob by
/// name would inherit a silent truncation ng's coverage runs cannot survive. The
/// cap is not a constraint in practice — production's own default of 5,000 is
/// already unreachable with Illumina reads, and the ceiling is thirteen times
/// that again (owner, 2026-07-29).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PileupGeneratorConfig {
    /// Reads folded at a position with no indel anchored there. Production: 8,000.
    pub max_snp_column_depth: u32,
    /// Reads folded at a position where any read has an indel. Production: 250.
    pub max_indel_column_depth: u32,
    /// Widest record footprint before the walk fails. Production: 5,000; ng
    /// additionally rejects anything above [`MAX_RECORD_SPAN_CEILING`].
    pub max_record_span: u32,
    /// How far a first mate stays available for pairing. Production: 10,000.
    pub mate_lookup_window: u32,
    /// Active-read ceiling. Production: 4,096.
    pub max_active_reads: u32,
}

impl Default for PileupGeneratorConfig {
    /// Production's five constants, **by name** — never by literal, so a retune
    /// on production's side reaches ng as a diff rather than as drift.
    fn default() -> Self {
        Self {
            max_snp_column_depth: crate::pileup::walker::DEFAULT_MAX_SNP_COLUMN_DEPTH,
            max_indel_column_depth: crate::pileup::walker::DEFAULT_MAX_INDEL_COLUMN_DEPTH,
            max_record_span: crate::pileup::walker::DEFAULT_MAX_RECORD_SPAN,
            mate_lookup_window: crate::pileup::walker::DEFAULT_MATE_LOOKUP_WINDOW,
            // ng's copy of the constant, which `walker_vocabulary_tests` pins equal
            // to production's — the one `DEFAULT_*` the verbatim copy forked, because
            // it is declared inside `chain_id_allocator.rs`.
            max_active_reads: super::DEFAULT_MAX_ACTIVE_READS,
        }
    }
}

impl PileupGeneratorConfig {
    /// Reject a configuration a [`ReadCoverage`] run could not describe.
    ///
    /// Called by [`PileupGenerator::new`], so a bad knob never reaches a locus.
    /// `coverage_of` carries a `debug_assert` stating the same envelope; this is
    /// what makes it provable rather than hopeful, in release builds too.
    ///
    /// [`ReadCoverage`]: crate::ng::locus_generation::ReadCoverage
    pub fn check(&self) -> Result<(), PileupGeneratorConfigError> {
        if self.max_record_span > MAX_RECORD_SPAN_CEILING {
            return Err(PileupGeneratorConfigError::RecordSpanExceedsCoverageRun {
                max_record_span: self.max_record_span,
                ceiling: MAX_RECORD_SPAN_CEILING,
            });
        }
        Ok(())
    }

    /// The same five knobs in the shape the copied walker reads them.
    ///
    /// Written as an exhaustive struct literal rather than with `..default()`:
    /// a knob production adds to [`WalkerConfig`] is then a compile error here,
    /// which is a decision to make rather than a default to inherit silently.
    pub fn to_walker_config(&self) -> WalkerConfig {
        WalkerConfig {
            max_snp_column_depth: self.max_snp_column_depth,
            max_indel_column_depth: self.max_indel_column_depth,
            max_record_span: self.max_record_span,
            mate_lookup_window: self.mate_lookup_window,
            max_active_reads: self.max_active_reads,
        }
    }
}

/// A configuration this generator cannot walk under. `#[non_exhaustive]`; raised
/// at construction, never mid-walk.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum PileupGeneratorConfigError {
    /// `max_record_span` is wider than a coverage run can describe, so a partial
    /// witness of a maximally-wide record would report a saturated
    /// `positions_covered` — silently, since [`LocusLen::from_positions`] clamps
    /// rather than failing.
    ///
    /// [`LocusLen::from_positions`]: crate::ng::locus_generation::LocusLen::from_positions
    #[error(
        "max_record_span ({max_record_span}) exceeds the widest footprint a coverage run can \
         describe ({ceiling}): a partially-witnessed record that wide would report a truncated \
         positions_covered rather than an error"
    )]
    RecordSpanExceedsCoverageRun { max_record_span: u32, ceiling: u32 },
}

/// Run-level counts for this generator, kept alongside the shared
/// [`LocusCounts`](super::super::LocusCounts) the dispatcher owns.
///
/// The first seven mirror production's [`RunSummary`](super::RunSummary) field for
/// field (spec §7); the last two are ng's. **`records_emitted` is deliberately not
/// mirrored**: what a caller wants is loci *kept*, which is the walk's emissions
/// minus [`records_outside_region`](Self::records_outside_region), and the kept
/// count is already `LocusCounts::loci_emitted`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PileupGeneratorCounts {
    /// Reads the walk admitted, across every segment.
    pub reads_admitted: u64,
    /// Times an open record's footprint grew under a longer deletion.
    pub record_widen_events: u64,
    /// Positions where two overlapping mates were reconciled.
    pub mate_overlap_positions: u64,
    /// Chain ids allocated — fragments the walk gave an identity.
    pub chain_allocations: u64,
    /// The largest number of concurrently-active reads seen. A **max**, not a
    /// sum: regions are walked one at a time, so the run's peak is the largest
    /// single region's peak.
    pub active_reads_high_water: u32,
    /// First mates evicted from `pending_mates` unpaired, having waited longer
    /// than `mate_lookup_window`.
    pub mate_lookup_evictions: u64,
    /// Columns where the contributor list was truncated by the applicable
    /// per-column depth cap.
    pub column_depth_truncations: u64,
    /// Reads silent over a whole record footprint — every base `N` or
    /// adaptor-masked, so never a contributor at any position and invisible to
    /// the per-locus `reads_without_observation` tally (spec §6).
    ///
    /// **Nothing increments this yet.** It needs a per-active-read "ever
    /// contributed" flag in the walk, which no step through C4 adds; the dump
    /// tool's read-accounting assertion (spec §13, plan D2) is what forces it.
    /// Until then the field reads zero, which is not the same claim as "no read
    /// was silent".
    pub reads_silent_over_footprint: u64,
    /// Records the region clamp dropped, because their anchor fell outside the
    /// region being walked. Observably zero-sum across neighbouring regions,
    /// which is how the gap-free tiling argument stays checkable (spec §7).
    pub records_outside_region: u64,
}

/// ng's generic locus generator: a streaming pileup walk over one `Generic`
/// region, emitting one
/// [`SampleLocusObservations`](super::super::SampleLocusObservations) per covered
/// position.
///
/// # What it holds, and why it holds it
///
/// A generator owns its accessors as fields (`locus_generation.md` §2). **Two
/// reference accessors, not one:** `preparer` carries its own (read
/// preparation's rule), and `reference` serves the walk's REF fetches. Neither is
/// rebuilt per segment — a fresh accessor per region throws away the sliding
/// buffer at every boundary and re-pays a `.fai` parse plus two `open(2)`s, which
/// is the ~564k-opens trap the STR side already paid for (spec §8).
///
/// The chain-id allocator likewise lives here rather than inside the walk: ng
/// walks one region where production walks a chromosome, and a fresh allocator
/// per segment would give two fragments of different regions the same id (spec
/// §8). What that costs — `reset()` between segments, and counters that must be
/// folded as deltas because `reset()` preserves them — is C3's.
pub struct PileupGenerator<R: RefSeq, P: ReadPreparer> {
    /// The reference the walk fetches REF bases from. Built once, for the run.
    #[expect(
        dead_code,
        reason = "C2's walk is the reader; C1 is the state that walk will run on"
    )]
    reference: R,
    /// Canonicalises each read the query returns before the walk sees it.
    #[expect(
        dead_code,
        reason = "C2's walk is the reader; C1 is the state that walk will run on"
    )]
    preparer: P,
    /// The preparer's reusable buffers — allocated once for the generator, never
    /// per read and never per segment.
    #[expect(
        dead_code,
        reason = "C2's walk is the reader; C1 is the state that walk will run on"
    )]
    prep_scratch: P::Scratch,
    /// Lives across segments so `next_id` never repeats.
    #[expect(
        dead_code,
        reason = "C3 resets it between segments and folds its counters as deltas"
    )]
    chain_ids: ChainIdAllocator,
    config: PileupGeneratorConfig,
    counts: PileupGeneratorCounts,
}

impl<R: RefSeq, P: ReadPreparer> PileupGenerator<R, P> {
    /// Build a generator over `reference` (the walk's REF fetches) and `preparer`
    /// (per-read canonicalisation), with `config` checked before anything is
    /// held.
    ///
    /// Fails only on a configuration a coverage run could not describe — see
    /// [`PileupGeneratorConfig::check`].
    pub fn new(
        reference: R,
        preparer: P,
        config: PileupGeneratorConfig,
    ) -> Result<Self, PileupGeneratorConfigError> {
        config.check()?;
        Ok(Self {
            reference,
            preparer,
            prep_scratch: P::Scratch::default(),
            chain_ids: ChainIdAllocator::with_caps(
                config.max_active_reads,
                config.mate_lookup_window,
            ),
            config,
            counts: PileupGeneratorCounts::default(),
        })
    }

    /// The knobs this generator was built with.
    pub fn config(&self) -> &PileupGeneratorConfig {
        &self.config
    }

    /// The run-level counts accumulated across every segment walked so far.
    pub fn counts(&self) -> &PileupGeneratorCounts {
        &self.counts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::ng::locus_generation::LocusLen;
    use crate::ng::read::{AlignedRead, PreparedRead, ReadPrepError};
    use crate::ng::ref_seq::InMemoryRefSeq;

    /// A preparer that prepares nothing — the generator's constructor never calls
    /// it, and C1 has no walk to feed.
    struct PreparesNothing;

    impl ReadPreparer for PreparesNothing {
        type Scratch = ();
        fn prepare_read(
            &self,
            _read: AlignedRead,
            _scratch: &mut Self::Scratch,
        ) -> Result<Option<PreparedRead>, ReadPrepError> {
            Ok(None)
        }
    }

    fn a_generator(
        config: PileupGeneratorConfig,
    ) -> Result<PileupGenerator<InMemoryRefSeq, PreparesNothing>, PileupGeneratorConfigError> {
        PileupGenerator::new(
            InMemoryRefSeq::from_contigs(vec![b"ACGTACGTAC".to_vec()]),
            PreparesNothing,
            config,
        )
    }

    /// **The defaults are production's, read from production's own constants.**
    ///
    /// Asserted against the `pub const`s rather than against literals: a literal
    /// here would let production retune a knob and ng silently keep the old
    /// value, which is the drift the "by name, not by literal" rule exists to
    /// prevent (arch §1.1).
    #[test]
    fn the_default_knobs_are_productions_five_constants() {
        let config = PileupGeneratorConfig::default();
        assert_eq!(
            config.max_snp_column_depth,
            crate::pileup::walker::DEFAULT_MAX_SNP_COLUMN_DEPTH
        );
        assert_eq!(
            config.max_indel_column_depth,
            crate::pileup::walker::DEFAULT_MAX_INDEL_COLUMN_DEPTH
        );
        assert_eq!(
            config.max_record_span,
            crate::pileup::walker::DEFAULT_MAX_RECORD_SPAN
        );
        assert_eq!(
            config.mate_lookup_window,
            crate::pileup::walker::DEFAULT_MATE_LOOKUP_WINDOW
        );
        assert_eq!(
            config.max_active_reads,
            crate::pileup::walker::DEFAULT_MAX_ACTIVE_READS
        );
    }

    /// The whole config, not just the one knob C1 constrains: production's
    /// defaults must all be walkable, or the `Default` impl ships a generator
    /// that cannot be constructed.
    #[test]
    fn the_default_config_is_accepted() {
        assert!(PileupGeneratorConfig::default().check().is_ok());
        assert!(a_generator(PileupGeneratorConfig::default()).is_ok());
    }

    /// **The ceiling is exactly where the narrowing starts lying.** This is the
    /// reason the knob is capped, asserted rather than described: at the ceiling
    /// `LocusLen` still reports the span it was given; one position past it, it
    /// reports a smaller number and no error.
    #[test]
    fn a_span_one_past_the_ceiling_is_where_a_coverage_run_starts_lying() {
        let at_ceiling = u64::from(MAX_RECORD_SPAN_CEILING);
        assert_eq!(
            u64::from(LocusLen::from_positions(at_ceiling).get()),
            at_ceiling,
            "a run at the ceiling must still describe itself exactly"
        );
        assert!(
            u64::from(LocusLen::from_positions(at_ceiling + 1).get()) < at_ceiling + 1,
            "one position past the ceiling the run saturates — silently, which is what \
             `check` exists to make unreachable"
        );
    }

    /// The boundary is inclusive: a span **at** the ceiling is describable, so it
    /// is accepted.
    #[test]
    fn a_max_record_span_at_the_ceiling_is_accepted() {
        let config = PileupGeneratorConfig {
            max_record_span: MAX_RECORD_SPAN_CEILING,
            ..PileupGeneratorConfig::default()
        };
        assert!(config.check().is_ok());
        assert!(a_generator(config).is_ok());
    }

    /// One position past it is rejected, and the error names both numbers —
    /// production's own `--max-record-span` is an unbounded `u32`, so this is the
    /// knob a user is most likely to raise past what ng can describe.
    #[test]
    fn a_max_record_span_past_the_ceiling_is_rejected() {
        let config = PileupGeneratorConfig {
            max_record_span: MAX_RECORD_SPAN_CEILING + 1,
            ..PileupGeneratorConfig::default()
        };
        let error = config
            .check()
            .expect_err("a span wider than a coverage run must be rejected");
        let PileupGeneratorConfigError::RecordSpanExceedsCoverageRun {
            max_record_span,
            ceiling,
        } = error;
        assert_eq!(max_record_span, MAX_RECORD_SPAN_CEILING + 1);
        assert_eq!(ceiling, MAX_RECORD_SPAN_CEILING);
    }

    /// **The rejection reaches the constructor**, which is the only door the
    /// generator has: a `check` nobody calls would leave the invariant a comment.
    #[test]
    fn the_constructor_rejects_a_span_no_coverage_run_could_describe() {
        let error = a_generator(PileupGeneratorConfig {
            max_record_span: MAX_RECORD_SPAN_CEILING + 1,
            ..PileupGeneratorConfig::default()
        })
        .err()
        .expect("the constructor must reject what `check` rejects");
        assert!(
            error.to_string().contains("65535"),
            "the message must name the ceiling the user has to come back under, got: {error}"
        );
    }

    /// **Every knob reaches the walker, and reaches the right field.** Five
    /// deliberately distinct, non-default values: with the defaults, or with two
    /// knobs sharing a value, a transposed pair of fields in `to_walker_config`
    /// would pass.
    #[test]
    fn every_knob_reaches_the_walker_config_it_names() {
        let config = PileupGeneratorConfig {
            max_snp_column_depth: 11,
            max_indel_column_depth: 22,
            max_record_span: 33,
            mate_lookup_window: 44,
            max_active_reads: 55,
        };
        let walker_config = config.to_walker_config();
        assert_eq!(walker_config.max_snp_column_depth, 11);
        assert_eq!(walker_config.max_indel_column_depth, 22);
        assert_eq!(walker_config.max_record_span, 33);
        assert_eq!(walker_config.mate_lookup_window, 44);
        assert_eq!(walker_config.max_active_reads, 55);
    }

    /// A fresh generator has counted nothing — the baseline every later delta is
    /// folded onto (C3).
    #[test]
    fn a_fresh_generator_has_counted_nothing() {
        let generator = a_generator(PileupGeneratorConfig::default()).expect("default config");
        assert_eq!(*generator.counts(), PileupGeneratorCounts::default());
        assert_eq!(*generator.config(), PileupGeneratorConfig::default());
    }
}
