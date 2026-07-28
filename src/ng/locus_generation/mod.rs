//! ng locus generation — typed regions → a sample's loci (the shared shape).
//!
//! The middle arrow of the caller spine: region typing says *what the reference is*
//! at every position, read ingestion says *what reads a sample has*; this step joins
//! them into **loci**, a locus being a stretch of genome with one sample's evidence
//! attached. The evidence is the same kind of thing for every kind of locus — the
//! distinct sequences the reads showed, each with its support — so one type
//! ([`SampleLocusObservations`]) serves a candidate SNP, a microsatellite, and a
//! repeat cluster alike; [`LocusKind`] tags which.
//!
//! This module owns the **shared shape**: the locus type, the `LocusGenerator`
//! contract, the dispatcher, and the count-only `NoLoci` generator (landing across
//! this plan's milestones). Each real generator plugs in from its own module —
//! `ssr.rs` (STR), `pileup/` (generic). See `doc/devel/ng/spec/locus_generation.md`
//! (design) and `doc/devel/ng/arch/locus_generation.md` (types & interfaces).

pub mod ssr;

use crate::ng::read::input::{IngestError, SampleReads};
use crate::ng::ref_seq::RefSeqError;
use crate::ng::region_typing::segment_criteria::{Motif, SsrSegment};
use crate::ng::region_typing::{RegionKind, TypedRegion, TypedRegionError};
use crate::ng::types::{GenomeRegion, ReadGroupId};
use crate::pileup_record::ChainId;

/// One sample's locus: the stretch of genome it covers, and what that sample's reads
/// showed there.
///
/// **Owned, no lifetimes** — a cohort stage merges these across samples and a future
/// artifact writes them, so it must outlive every buffer it was built from
/// (spec §3). The evidence is uniform across kinds; [`kind`](Self::kind) is what names
/// the locus, not which fields are populated.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleLocusObservations {
    /// The stretch this locus covers — one base for a candidate SNP, several for an
    /// indel, the tract for a microsatellite.
    pub region: GenomeRegion,
    /// The reference (REF) bases over `region` — what a wider-span projection needs
    /// when samples merge.
    pub reference_bases: Box<[u8]>,
    /// The distinct sequences the reads showed, each with its support. **Observations,
    /// not alleles** — they become alleles when something calls them.
    pub observed_sequences: Vec<ObservedSequence>,
    /// Reads that covered this locus but produced no observation at all. A scalar with
    /// no positions: *no coverage* and *coverage that said nothing* are different
    /// states, and only one means "look at the mapping" (spec §3).
    pub reads_without_observation: u32,
    /// Reads a depth cap discarded. **Non-zero means the support counts are a
    /// subsample, not the depth** (spec §3).
    pub reads_discarded_by_cap: u32,
    /// What kind of locus this is, and the extras that kind carries — read off the
    /// type, never inferred from which fields are populated.
    pub kind: LocusKind,
}

impl SampleLocusObservations {
    /// Read depth at each position of `region`, in order — **derived, not stored**.
    ///
    /// A [`Complete`](ReadCoverage::Complete) observation counts its `num_obs` at every
    /// position; an [`Observed`](ReadCoverage::Observed) run counts it over the stretch
    /// it witnessed. The returned vector has exactly `region.len()` entries.
    ///
    /// This is *observation* depth and only exact per locus: it omits reads that
    /// covered the tract but anchored no border (they are in
    /// [`reads_without_observation`](Self::reads_without_observation), a scalar with no
    /// positions), and overlapping loci on the generic path double-count if summed. The
    /// paralog filter owns both caveats (spec §3, §11).
    pub fn num_obs_along_locus(&self) -> Vec<u32> {
        let len = self.region.len() as usize;
        let mut depth = vec![0u32; len];
        for obs in &self.observed_sequences {
            // A partial's covered extent cannot exceed the locus length — that is a
            // producer invariant, enforced where `ReadCoverage` is minted (the STR
            // generator). Here we *clamp* rather than `debug_assert`: this is a
            // consumer-side derivation run over whole cohorts, and a debug-only guard
            // compiles out of the release build this repo actually runs (a trap it has
            // recorded hitting twice). Clamping keeps the derivation panic-free on any
            // input; a bad extent is caught at the producer, not here.
            let (from, to) = match obs.read_coverage {
                ReadCoverage::Complete => (0, len),
                ReadCoverage::Observed {
                    offset_in_locus,
                    positions_covered,
                } => {
                    let from = (offset_in_locus as usize).min(len);
                    (
                        from,
                        from.saturating_add(positions_covered as usize).min(len),
                    )
                }
            };
            for slot in &mut depth[from..to] {
                *slot = slot.saturating_add(obs.num_obs);
            }
        }
        depth
    }

    /// This locus's length, for the [`ReadCoverage`] predicates that need it.
    ///
    /// **The one source a consumer should use.** The mint derives the same quantity from
    /// the segment it is building the locus from, before the locus exists; every reader
    /// afterwards should ask the locus rather than re-deriving from `region`, so the two
    /// cannot drift apart unnoticed. That they agree is pinned by
    /// `ssr::tests::the_mints_locus_length_is_the_emitted_regions`.
    pub fn locus_len(&self) -> LocusLen {
        LocusLen::of_region(self.region)
    }

    /// The observations a likelihood may score directly — the
    /// [`Complete`](ReadCoverage::Complete) ones.
    ///
    /// A partial is a lower bound that mis-scores as a *short* allele until a censored
    /// likelihood models it (step 7), so reaching the partials is a deliberate act:
    /// this iterator is the guard (spec §3).
    pub fn complete_observations(&self) -> impl Iterator<Item = &ObservedSequence> + '_ {
        self.observed_sequences
            .iter()
            .filter(|obs| obs.read_coverage == ReadCoverage::Complete)
    }
}

/// One distinct sequence the reads showed at a locus, with its support.
///
/// The fields between `num_obs` and `chain_ids` are the per-read moments the SNP
/// filters read (strand bias, base-quality error, the MAPQ multi-mapper test); an STR
/// model reduces to `num_obs` alone. Modelled on production's per-allele shape
/// (`AlleleObservation` + `AlleleSupportStats`), minus `placed_start`, which no model
/// consumes (spec §6).
///
/// **This is a table of cells, not a table of sequences.** The identity is
/// `(bases, read_coverage, read_group)` — three axes, not one — so a consumer that
/// wants per-allele totals must aggregate over coverage *and* group, and one that
/// treats each entry as an allele will count the same allele several times. The
/// aggregation is exact: every support field is additive, and the merged cells share
/// their `bases` and `read_coverage` by construction (spec §6).
#[derive(Debug, Clone, PartialEq)]
pub struct ObservedSequence {
    /// The observed bases — allele content, in **read** coordinates.
    pub bases: Box<[u8]>,
    /// How much of the locus a read of this sequence spanned. **Part of the
    /// identity**: a [`Complete`](ReadCoverage::Complete) and an
    /// [`Observed`](ReadCoverage::Observed) run of the same `bases` are different
    /// evidence and stay separate entries (spec §3).
    pub read_coverage: ReadCoverage,
    /// Which read group — one `@RG`, i.e. one lane — these reads came from. **Part of
    /// the identity**, so an allele supported from several groups is several rows.
    ///
    /// Carried because a per-chemistry model needs the allele × group cross **with its
    /// quality moments**: a per-group count beside one merged observation gives the
    /// first and loses the second. The near-term consumer is the STR path, whose
    /// stutter level and per-base `ε` are already fit per sample group — off groups it
    /// currently has to *infer* ("data-driven soft clusters") because the evidence did
    /// not carry the real one. The read group is it (spec §6).
    ///
    /// At the **finest** grain available, deliberately: library and experiment stay a
    /// downstream fold, so picking a grain remains the modeller's decision and not this
    /// step's guess. Free where a sample has one read group, which is most of them.
    pub read_group: ReadGroupId,
    /// How many reads showed this sequence. The whole support on the STR path, and the
    /// one field every model on both paths reduces to.
    pub num_obs: u32,
    /// Reads on the forward strand — strand bias.
    pub num_fwd: u32,
    /// Σ per-read log-error over the supporting reads — the freebayes per-read error
    /// term (production's `q_sum`).
    pub q_sum: f64,
    /// Σ MAPQ over the supporting reads. With `mapq_sum_sq` and `num_obs` it recovers
    /// the mean and variance the MAPQ Welch's-t multi-mapper filter reads.
    pub mapq_sum: u32,
    /// Σ MAPQ² over the supporting reads (see `mapq_sum`).
    pub mapq_sum_sq: u64,
    /// How many supporting reads started **strictly left** of the locus's anchor —
    /// freebayes' `placedLeft`, and the read-position-bias term production subtracts
    /// from QUAL (`vcf/qual_refine.rs`).
    ///
    /// Carried because dropping it would forfeit the ability to reproduce production's
    /// QUAL, which outranks tidiness. Its sibling `placed_start` is **not** carried: no
    /// model consumes it, and it is a pure function of the read's start against the
    /// anchor, so a later consumer can re-derive it without changing the fold (spec §6).
    pub placed_left: u32,
    /// Phase-chain ids of the reads folded here — what lets a later step chain
    /// observations at neighbouring loci into a haplotype.
    pub chain_ids: Vec<ChainId>,
}

/// How much of a locus a single read spanned — **one read's span, not depth**.
///
/// `Complete` means the read reached **both** borders of the locus; anything else is
/// the one **run** of locus positions the read actually witnessed, in **locus**
/// coordinates (spec §3). A partial run is a *censored* observation: the sequence is
/// at least this long, but not how long.
///
/// **One `Observed` run replaces the earlier `PartialLeft`/`PartialRight` pair
/// (owner, 2026-07-28).** Two side-tagged variants cannot describe what a read
/// witnesses once the *events*, not the alignment span, define it: a read can be blind
/// in the **middle** of a footprint (an interior `N`, a ref-skip) or blind at either
/// end, and a widened record can be wider than the read on both sides. Prefix-versus-
/// suffix survives as a derivation — [`is_flush_left`](Self::is_flush_left) /
/// [`is_flush_right`](Self::is_flush_right) — so the STR path's "a prefix and a suffix
/// are different constraints" is preserved, not lost. `Complete` is kept rather than
/// folded in: it is the overwhelmingly common case and it keeps
/// [`complete_observations`](SampleLocusObservations::complete_observations) a cheap
/// equality instead of arithmetic against the footprint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReadCoverage {
    /// The read reached both borders of the locus.
    Complete,
    /// The stretch the read did witness, in **locus positions** — the axis `bases` is
    /// not on (that is allele content, in read coordinates).
    Observed {
        /// Locus positions between the locus's left border and the first one
        /// witnessed. `0` = flush with the left border, i.e. a prefix constraint.
        offset_in_locus: u16,
        /// How many locus positions were witnessed, running from `offset_in_locus`.
        positions_covered: u16,
    },
}

/// A locus's length in reference positions — the axis a [`ReadCoverage`] run lives on.
///
/// **A newtype because the alternative is silently wrong.** Every constructor and
/// predicate on `ReadCoverage` takes a covered extent *and* a locus length, both counts
/// of locus positions and both formerly `u16` — so `from_left(10, 4)` and
/// `from_left(4, 10)` each compiled, and the clamping the constructors do for their own
/// good would then hide the transposition rather than surface it. Two `u16`s in a row
/// with no way to tell them apart is exactly the shape that produces a wrong depth and
/// no panic.
///
/// It also gives the saturating cast **one** home. The count arrives as a `u64` (a
/// region length, a tract length) and has to be narrowed; doing that at each call site
/// spread the same `.min(u16::MAX as u64) as u16` across the mint and every dump tool.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocusLen(u16);

impl LocusLen {
    /// From a count of reference positions, saturating at `u16::MAX`.
    ///
    /// Saturation rather than an error: a locus longer than 65,535 positions cannot be
    /// described by a run either, and the clamp keeps the derivation total. A tract that
    /// long is a satellite, which the caller does not handle.
    pub fn from_positions(positions: u64) -> Self {
        Self(positions.min(u16::MAX as u64) as u16)
    }

    /// The length of `region` — the canonical source once a locus exists, and what
    /// [`SampleLocusObservations::locus_len`] returns.
    pub fn of_region(region: GenomeRegion) -> Self {
        Self::from_positions(region.len())
    }

    pub fn get(self) -> u16 {
        self.0
    }
}

impl ReadCoverage {
    /// A run flush with the locus's **left** border, `positions_covered` long — the
    /// old `PartialLeft(n)`.
    ///
    /// `locus_len` clamps the run into the locus. It is needed because a producer's
    /// reach can be measured in *read* bases, which diverge from locus positions under
    /// stutter; a run must never claim positions the locus does not have.
    pub fn from_left(positions_covered: u16, locus_len: LocusLen) -> Self {
        Self::Observed {
            offset_in_locus: 0,
            positions_covered: positions_covered.min(locus_len.get()),
        }
    }

    /// A run flush with the locus's **right** border, `positions_covered` long — the
    /// old `PartialRight(n)`.
    ///
    /// Clamped first, then the offset derived, so the subtraction cannot underflow:
    /// computing `locus_len - positions_covered` on an over-long reach would wrap to a
    /// huge offset, or — with a saturating subtraction — silently relabel a
    /// right-anchored read as a left-anchored one.
    ///
    /// **One consequence of the encoding, and it reaches further than labelling.** A run
    /// that covers the whole locus is flush with *both* borders, so once
    /// `positions_covered >= locus_len` this returns exactly what
    /// [`from_left`](Self::from_left) would. Three things follow, and the third is a
    /// behaviour change:
    ///
    /// 1. the run is reported as flush left by [`is_flush_left`](Self::is_flush_left);
    /// 2. every label derived from flushness calls it a *left* partial;
    /// 3. **it shares a bucket key with a left-flush run of the same bases, so the two
    ///    merge into one row** — where `PartialLeft(n)` and `PartialRight(n)` kept them
    ///    apart.
    ///
    /// That is reachable on the STR path, where the reach is measured in *read* bases:
    /// an allele longer than the reference tract gives a reach past the locus length.
    /// It is arguably the right answer — identical constraints are one cell, and a read
    /// that witnessed every position is constrained from neither side — but it is not
    /// the pre-reshape answer, and the plan's stated equivalence
    /// `PartialRight(n) ⇔ Observed { len - n, n }` stops holding at `n = len`. Pinned by
    /// `ssr::tally::tests::an_expanded_allele_merges_the_two_sides_into_one_row`.
    pub fn from_right(positions_covered: u16, locus_len: LocusLen) -> Self {
        let covered = positions_covered.min(locus_len.get());
        Self::Observed {
            offset_in_locus: locus_len.get() - covered,
            positions_covered: covered,
        }
    }

    /// Whether the run starts at the locus's left border — a **prefix** constraint on
    /// the allele. Always true of `Complete`.
    pub fn is_flush_left(&self) -> bool {
        match self {
            Self::Complete => true,
            Self::Observed {
                offset_in_locus, ..
            } => *offset_in_locus == 0,
        }
    }

    /// Whether the run ends at the locus's right border — a **suffix** constraint.
    /// Always true of `Complete`.
    pub fn is_flush_right(&self, locus_len: LocusLen) -> bool {
        match self {
            Self::Complete => true,
            Self::Observed {
                offset_in_locus,
                positions_covered,
            } => offset_in_locus.saturating_add(*positions_covered) >= locus_len.get(),
        }
    }
}

/// The kind of locus, plus whatever that kind adds to the shared evidence fields.
///
/// `#[non_exhaustive]` because a kind's payload fills in as the generator that mints
/// it is written. Shared vocabulary (`ng_step_interfaces.md` §2); authoritative in
/// spec §3.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum LocusKind {
    /// A SNP/indel candidate site. Its evidence is the observed alleles; no extras.
    Generic,
    /// A microsatellite tract — carries the motif and the flanks the read model needs.
    Ssr(SsrDetail),
    /// A repeat cluster with no clean flanks, coarser than a single tract. Its payload
    /// is the bundle generator's to decide (deferred, spec §11).
    SsrBundle,
}

/// What an [`LocusKind::Ssr`] locus carries — grouped so a repeat's motif and flanks
/// are present or absent together, never half.
///
/// The STR generator mints these, splitting the flanks out of the reference bases it
/// fetches (spec §3; `locus_generation_ssr.md` §1).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrDetail {
    /// The repeat unit.
    pub motif: Motif,
    /// The reference flank left of the tract — the read model's alignment anchor.
    /// Clamped at the contig end, so it may be shorter than the right flank.
    pub left_flank: Box<[u8]>,
    /// The reference flank right of the tract (see `left_flank`).
    pub right_flank: Box<[u8]>,
}

/// Generates a sample's loci from one segment of kind `S`, streaming **one locus at a
/// time**.
///
/// `S` is the segment payload the generator consumes — `SsrSegment` for the STR generator.
/// It is a parameter on the *contract*, not an associated type inside each implementation,
/// so two generators for the same kind stay interchangeable behind `Box<dyn
/// LocusGenerator<S>>` (spec §4). A generator holds its own accessors (reference, aligner,
/// scratch) as fields, so the only per-call context is the segment and the sample's reads.
pub trait LocusGenerator<S> {
    /// Start a new segment: reset progress. Does no gathering and cannot fail.
    fn begin_segment(&mut self, region: GenomeRegion);

    /// The next locus of the segment begun, or `None` once it has no more.
    ///
    /// Called repeatedly with the same `segment` until it returns `None`; returning `None`
    /// immediately is a normal outcome, not a failure. The `segment` must be the one whose
    /// region was passed to the preceding [`begin_segment`](Self::begin_segment) — the two
    /// calls are paired, and nothing in the types enforces it. `&mut self` because a
    /// generator owns reusable scratch (alignment matrices, sampling buffers) that must not
    /// be reallocated per segment.
    fn next_locus(
        &mut self,
        segment: &S,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError>;
}

/// A generator that produces no loci and reports why.
///
/// One implementation covers every kind, because it ignores the segment entirely — the
/// count-only fallback every region kind with no real generator resolves to, so that "we
/// produce nothing here" is a configuration with a reason attached rather than a silent
/// gap (spec §5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NoLoci {
    /// Why this kind produces no loci — the fact the dispatcher accounts by.
    pub reason: UnhandledReason,
}

/// **Why** a kind produces no loci — a boundary we chose vs. a gap not yet filled.
///
/// Not cosmetic: the two answer different questions ("what will this caller never cover?"
/// vs "how much does it not cover *yet*?") and must not be added together (spec §5, §7).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnhandledReason {
    /// Deliberately outside the caller's scope — e.g. satellite arrays. Permanent.
    OutOfScope,
    /// No generator written yet. Temporary by construction.
    NotImplemented,
}

impl<S> LocusGenerator<S> for NoLoci {
    fn begin_segment(&mut self, _region: GenomeRegion) {}

    fn next_locus(
        &mut self,
        _segment: &S,
        _reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        Ok(None)
    }
}

/// A fatal, run-level failure of locus generation.
///
/// `#[non_exhaustive]`; every variant wraps an upstream error — this step mints none of its
/// own. A read that yields no observation is a tallied per-read outcome, never an error; an
/// error means the run is broken (spec §6). A reference fetch can surface two ways — through
/// the upstream walk ([`TypedRegion`](Self::TypedRegion)) or a generator's own fetch
/// ([`Reference`](Self::Reference)) — and they stay distinct because they fail in different
/// places.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum LocusGenerationError {
    /// The upstream typed-region walk failed.
    #[error("typed-region walk failed during locus generation")]
    TypedRegion(#[from] TypedRegionError),
    /// A read query failed mid-stream, or the alignment input was malformed (the open
    /// already succeeded upstream; this fires while a generator pulls reads).
    #[error("read access failed during locus generation")]
    Reads(#[from] IngestError),
    /// A reference fetch failed — a broken reference, or a region past a contig end.
    #[error("reference fetch failed during locus generation")]
    Reference(#[from] RefSeqError),
}

/// The running tally — "no silent caps": every region and every base is accounted for, so
/// "how much genome does this caller not cover, and how much of that is temporary?" is
/// answerable from the counts alone. The base counters are the other half of why `SsrBundle`
/// and `Satellite` exist as types rather than holes (spec §7).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LocusCounts {
    /// Typed regions dispatched — the total, which partitions **exactly** into
    /// `regions_handled` plus the two unhandled counters (spec §13.2).
    pub regions_in: u64,
    /// Regions routed to a filled generator, whatever number of loci it then emitted
    /// (including zero). With the two unhandled counters this sums to `regions_in`.
    pub regions_handled: u64,
    /// Loci emitted, across every generator. **Not** a region count — one handled region
    /// yields zero, one, or many.
    pub loci_emitted: u64,
    /// Regions that produced no loci because no generator is filled for their kind.
    /// **Temporary** by construction.
    pub unhandled_not_implemented: u64,
    /// The bases those `unhandled_not_implemented` regions cover.
    pub unhandled_not_implemented_bp: u64,
    /// Regions deliberately outside scope (satellites). **Permanent.**
    pub unhandled_out_of_scope: u64,
    /// The bases those `unhandled_out_of_scope` regions cover.
    pub unhandled_out_of_scope_bp: u64,
}

impl LocusCounts {
    /// Charge one unhandled region, of `bp` bases, to the counter its reason names — the one
    /// place the two kinds of nothing are kept apart (spec §5, §7).
    fn record_unhandled(&mut self, reason: UnhandledReason, bp: u64) {
        match reason {
            UnhandledReason::NotImplemented => {
                self.unhandled_not_implemented += 1;
                self.unhandled_not_implemented_bp += bp;
            }
            UnhandledReason::OutOfScope => {
                self.unhandled_out_of_scope += 1;
                self.unhandled_out_of_scope_bp += bp;
            }
        }
    }
}

/// One region kind's generator, or the reason it has none.
///
/// A **trait object** so a generator can be swapped per run without the dispatcher being
/// generic over each kind's concrete type — the lab's `Box<dyn _>` choice
/// (`ng_step_interfaces.md` §4). `Unfilled` carries the reason the dispatcher accounts by:
/// the `NoLoci` configuration kept as data, so plugging in a real generator is a one-line
/// change at the set (spec §5).
///
/// The trait object carries **no `Send` bound**: v1 is single-threaded (arch §9). If a
/// `GeneratorSet` is ever moved onto a producer thread rather than built per thread, this
/// becomes `dyn LocusGenerator<S> + Send` — a deliberate omission now, not an oversight.
pub enum GeneratorSlot<S> {
    /// A generator supplied for this kind.
    Generator(Box<dyn LocusGenerator<S>>),
    /// No generator; account every region of this kind to the reason.
    Unfilled(UnhandledReason),
}

impl<S> GeneratorSlot<S> {
    /// Begin a region on this slot: reset a real generator, or account the region as
    /// unhandled. Returns whether a generator is filled, so the dispatcher knows whether
    /// there are loci to pull.
    fn begin(&mut self, region: GenomeRegion, bp: u64, counts: &mut LocusCounts) -> bool {
        match self {
            GeneratorSlot::Generator(generator) => {
                generator.begin_segment(region);
                true
            }
            GeneratorSlot::Unfilled(reason) => {
                counts.record_unhandled(*reason, bp);
                false
            }
        }
    }

    /// The next locus from a filled slot. An unfilled slot yields `None` — though the
    /// dispatcher never asks one, since [`begin`](Self::begin) reported it not filled.
    fn next(
        &mut self,
        segment: &S,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        match self {
            GeneratorSlot::Generator(generator) => generator.next_locus(segment, reads),
            GeneratorSlot::Unfilled(_) => Ok(None),
        }
    }
}

/// The set of generators the dispatcher routes to — one slot per region kind — plus the
/// running tally and the one-region-at-a-time cursor.
///
/// `Satellite` has no slot: it is out of scope for the whole caller and always accounted
/// `OutOfScope` (spec §5). The other three kinds each hold a [`GeneratorSlot`]. This is the
/// `GeneratorSet` the arch left as an impl-time confirmation, pinned here; the payload types
/// for `Generic` and `SsrBundle` are `()` for now and refine when those generators land.
pub struct GeneratorSet {
    ssr: GeneratorSlot<SsrSegment>,
    generic: GeneratorSlot<()>,
    ssr_bundle: GeneratorSlot<()>,
    counts: LocusCounts,
    /// The region whose generator is mid-stream, if any. `None` between regions.
    current: Option<TypedRegion>,
}

impl GeneratorSet {
    /// A set with a generator (or a reason) chosen for each kind.
    pub fn new(
        ssr: GeneratorSlot<SsrSegment>,
        generic: GeneratorSlot<()>,
        ssr_bundle: GeneratorSlot<()>,
    ) -> Self {
        Self {
            ssr,
            generic,
            ssr_bundle,
            counts: LocusCounts::default(),
            current: None,
        }
    }

    /// A set with no real generator — every kind falls back to its `NoLoci` reason, which is
    /// what this shape ships (spec §2): `SsrSegment` / `Generic` / `SsrBundle` are
    /// `NotImplemented` until a generator is supplied; `Satellite` is always `OutOfScope`.
    pub fn all_unimplemented() -> Self {
        Self::new(
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        )
    }

    /// The running tally — readable at any point, final once the stream is exhausted.
    pub fn counts(&self) -> &LocusCounts {
        &self.counts
    }

    /// Begin a region: count it, and ready its generator if one is filled. Every region is
    /// counted in `regions_in`; a handled kind also in `regions_handled`, an unfilled kind in
    /// its unhandled counter. Infallible — resetting a generator cannot fail (spec §4).
    ///
    /// Call only after the previous region is drained (`next_locus` returned `None`); calling
    /// over an undrained region silently abandons its remaining loci. The iterator upholds
    /// this, which is why it is a documented contract rather than a runtime guard.
    pub fn begin_region(&mut self, region: TypedRegion) {
        self.counts.regions_in += 1;
        let bp = region.region.len();
        // Copied out (GenomeRegion is Copy) only for readability, so `region` can still move
        // into `current` below.
        let geometry = region.region;
        let filled = match &region.kind {
            RegionKind::Satellite => {
                self.counts
                    .record_unhandled(UnhandledReason::OutOfScope, bp);
                false
            }
            RegionKind::SsrSegment(_) => self.ssr.begin(geometry, bp, &mut self.counts),
            RegionKind::Generic => self.generic.begin(geometry, bp, &mut self.counts),
            RegionKind::SsrBundle { .. } => self.ssr_bundle.begin(geometry, bp, &mut self.counts),
        };
        if filled {
            self.counts.regions_handled += 1;
        }
        self.current = filled.then_some(region);
    }

    /// The next locus of the region begun by [`begin_region`](Self::begin_region), or `None`
    /// once it — or an unfilled/absent region — has no more. After a `None` the caller pulls
    /// the next region and calls `begin_region` again. Holds **one region at a time**: no
    /// buffer of loci (spec §6).
    pub fn next_locus(
        &mut self,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        let Some(region) = self.current.take() else {
            return Ok(None);
        };
        let produced = match &region.kind {
            RegionKind::SsrSegment(segment) => self.ssr.next(segment, reads),
            RegionKind::Generic => self.generic.next(&(), reads),
            RegionKind::SsrBundle { .. } => self.ssr_bundle.next(&(), reads),
            // Satellite is never made current — begin_region reports it unfilled.
            RegionKind::Satellite => Ok(None),
        }?;
        match produced {
            Some(locus) => {
                self.counts.loci_emitted += 1;
                self.current = Some(region); // more may follow; keep driving it.
                Ok(Some(locus))
            }
            None => Ok(None),
        }
    }
}

/// Lazily turns a typed-region stream into a sample's loci — the public surface of locus
/// generation.
///
/// Holds **no buffer of loci**: it drives the current region one locus at a time via the
/// [`GeneratorSet`], and only when that region is exhausted pulls the next from the stream,
/// so exactly one locus is resident regardless of how many a region yields (spec §6, §9).
/// Loci come out in the stream's order, which is coordinate order (spec §2).
///
/// Generic over the region stream `T` so a `Vec` can stand in for the real
/// `TypedRegionIterator` in tests. The generator set is a concrete [`GeneratorSet`] — the
/// per-run swap lives in its trait-object slots, not in a type parameter.
pub struct SampleLocusObservationsIterator<T> {
    regions: T,
    // FIXME(pileup-generator): the field order below loses a step-1 tally as
    // soon as a generator holds a read stream across `next_locus` calls —
    // which is exactly what the generic (pileup) locus generator is being
    // built to do. **Delete this comment once it is fixed.**
    //
    // Since 2026-07-28 a region stream owns its files by `Arc` rather than
    // borrowing them (`read/input/`, the generic generator's Milestone A), so
    // it may outlive the `SampleReads` that made it. Rust drops fields in
    // declaration order, so `reads` dies here *before* `generators` — and any
    // stream a generator is still holding then folds its drop tally into an
    // `AlignmentFile` that only that stream owns, and which is freed in the
    // same breath. The reads already emitted are unaffected and the pooled
    // reader still goes back; what is lost is the ability to *read* that
    // query's tally through `SampleReads::counts`, so the failure is a silent
    // under-report of drop rates, not a crash.
    //
    // Two fixes, both cheap, and the choice belongs with whoever writes the
    // generic generator: declare `generators` before `reads` so the streams
    // drop first, or have the generator drop its held stream at the end of
    // each segment. Whichever is taken, it wants a test that keeps a second
    // `Arc` and asserts the tally landed — see
    // `read::input::open_bam::tests::a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally`
    // for the shape, since a tally folded into an unreachable file is
    // indistinguishable from one never folded at all.
    reads: SampleReads,
    generators: GeneratorSet,
    /// Latched on clean exhaustion or a fatal error — the fused contract.
    done: bool,
}

impl<T> SampleLocusObservationsIterator<T> {
    /// `regions` is the typed-region stream, `reads` the sample's reads, `generators` the set
    /// the dispatcher routes to (spec §6). (No `LocusConfig` yet — it lands when it has a
    /// field; an empty one would be a dormant lever.)
    pub fn new(regions: T, reads: SampleReads, generators: GeneratorSet) -> Self {
        Self {
            regions,
            reads,
            generators,
            done: false,
        }
    }

    /// The running tally — current at any point, final once the stream is exhausted.
    pub fn counts(&self) -> &LocusCounts {
        self.generators.counts()
    }
}

impl<T> Iterator for SampleLocusObservationsIterator<T>
where
    T: Iterator<Item = Result<TypedRegion, TypedRegionError>>,
{
    type Item = Result<SampleLocusObservations, LocusGenerationError>;

    /// Pull loci from the current region; when it is exhausted, take the next region and
    /// begin it. A fatal error — from a generator or the upstream walk — is yielded once as
    /// `Some(Err(_))` and then the iterator is done, so `?` makes it un-ignorable rather than
    /// a silent end of stream (spec §6).
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        loop {
            match self.generators.next_locus(&self.reads) {
                Ok(Some(locus)) => return Some(Ok(locus)),
                Err(error) => {
                    self.done = true;
                    return Some(Err(error));
                }
                Ok(None) => match self.regions.next() {
                    None => {
                        self.done = true;
                        return None;
                    }
                    Some(Err(error)) => {
                        self.done = true;
                        return Some(Err(error.into()));
                    }
                    Some(Ok(region)) => self.generators.begin_region(region),
                },
            }
        }
    }
}

impl<T> std::iter::FusedIterator for SampleLocusObservationsIterator<T> where
    T: Iterator<Item = Result<TypedRegion, TypedRegionError>>
{
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::{ContigId, Position, ReadGroupId};

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// An observation of `bases` with `num_obs` reads at a given coverage — the moment
    /// fields are irrelevant to the depth derivation, so they are fixed.
    fn obs(bases: &[u8], read_coverage: ReadCoverage, num_obs: u32) -> ObservedSequence {
        ObservedSequence {
            bases: Box::from(bases),
            read_coverage,
            read_group: ReadGroupId(0),
            num_obs,
            num_fwd: 0,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn locus(region: GenomeRegion, observed: Vec<ObservedSequence>) -> SampleLocusObservations {
        SampleLocusObservations {
            region,
            reference_bases: Box::from(&b""[..]),
            observed_sequences: observed,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// A minimal real `SampleReads` over the read-ingestion test fixtures — one indexed BAM
    /// naming one sample, opened against the fixture reference. Constructing `SampleReads`
    /// needs alignment files, so this is the cheapest honest handle; the `NoLoci` path never
    /// reads it, but the signature requires one. Returns the temp dirs so they outlive the
    /// handle.
    fn sample_reads_over_fixture() -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        use crate::ng::read::filtering::ReadFilterConfig;
        use crate::ng::read::input::test_fixtures::{
            fixture_reference, header, indexed_bam, matching_contigs, read_named_with_length,
        };

        let (reference_dir, reference) = fixture_reference(false);
        let records = vec![read_named_with_length("r0", 0, 1, 30)];
        let (bam_dir, bam_path) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            &records,
        );
        let reads = SampleReads::open_only_sample(
            &[bam_path],
            &reference,
            ReadFilterConfig::default(),
            false,
        )
        .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    /// `NoLoci` is a `LocusGenerator` for any segment type, emits no locus, and carries its
    /// reason for the dispatcher to account by (spec §5).
    #[test]
    fn no_loci_emits_nothing_and_carries_its_reason() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let mut generator = NoLoci {
            reason: UnhandledReason::OutOfScope,
        };
        // Driven over `()` as the segment — NoLoci ignores it, as it does every kind. The
        // segment type must be named because NoLoci implements the trait for *every* `S`.
        LocusGenerator::<()>::begin_segment(&mut generator, region(1, 5));
        let out = LocusGenerator::<()>::next_locus(&mut generator, &(), &reads).unwrap();
        assert!(out.is_none(), "NoLoci produces no locus");
        assert_eq!(generator.reason, UnhandledReason::OutOfScope);
    }

    fn typed(kind: RegionKind, start: u64, end: u64) -> TypedRegion {
        TypedRegion {
            region: region(start, end),
            kind,
        }
    }

    fn an_ssr_segment(start: u64, end: u64) -> RegionKind {
        RegionKind::SsrSegment(
            SsrSegment::new("chr1".into(), start, end, Motif::new(b"AT").unwrap(), 1.0).unwrap(),
        )
    }

    fn a_bundle() -> RegionKind {
        use crate::ng::tandem_repeat::RepeatInterval;
        RegionKind::SsrBundle {
            tracts: vec![RepeatInterval {
                start: 99,
                end: 160,
                period: 2,
                score: 10,
            }]
            .into_boxed_slice(),
        }
    }

    /// Drive one region through the dispatcher to exhaustion, collecting its loci.
    fn drain_region(
        set: &mut GeneratorSet,
        region: TypedRegion,
        reads: &SampleReads,
    ) -> Vec<SampleLocusObservations> {
        set.begin_region(region);
        let mut out = Vec::new();
        while let Some(locus) = set.next_locus(reads).unwrap() {
            out.push(locus);
        }
        out
    }

    /// Every kind is accounted, and the two kinds of nothing land in **different** counters
    /// — the check that §5's distinction survives contact with code (spec §13.1, §13.3).
    #[test]
    fn all_unimplemented_accounts_every_kind_and_keeps_the_two_nothings_apart() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let mut set = GeneratorSet::all_unimplemented();

        // Distinct spans so the base counters are individually checkable.
        for region in [
            typed(RegionKind::Generic, 1, 10),      // 10 bp → NotImplemented
            typed(an_ssr_segment(20, 25), 20, 25),  // 6 bp → NotImplemented
            typed(a_bundle(), 100, 160),            // 61 bp → NotImplemented
            typed(RegionKind::Satellite, 200, 400), // 201 bp → OutOfScope
        ] {
            assert!(
                drain_region(&mut set, region, &reads).is_empty(),
                "an unimplemented set emits no loci"
            );
        }

        let counts = set.counts();
        assert_eq!(counts.regions_in, 4);
        assert_eq!(counts.regions_handled, 0);
        assert_eq!(counts.loci_emitted, 0);
        assert_eq!(counts.unhandled_not_implemented, 3);
        assert_eq!(counts.unhandled_not_implemented_bp, 10 + 6 + 61);
        assert_eq!(counts.unhandled_out_of_scope, 1);
        assert_eq!(counts.unhandled_out_of_scope_bp, 201);
        // Nothing is unaccounted for: regions_in partitions exactly (spec §13.2).
        assert_eq!(
            counts.regions_in,
            counts.regions_handled
                + counts.unhandled_not_implemented
                + counts.unhandled_out_of_scope,
        );
    }

    /// A generator emitting a fixed number of loci per segment — a stand-in for a real one,
    /// so the filled-slot path (loci counted, region *handled* not *unhandled*) is exercised
    /// even though this shape ships only `NoLoci`. Generic over the segment type so it fits
    /// any kind's slot, which is what lets one routing test distinguish the three.
    struct FixedCountGenerator {
        per_segment: u32,
        remaining: u32,
    }

    impl<S> LocusGenerator<S> for FixedCountGenerator {
        fn begin_segment(&mut self, _region: GenomeRegion) {
            self.remaining = self.per_segment;
        }

        fn next_locus(
            &mut self,
            _segment: &S,
            _reads: &SampleReads,
        ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
            if self.remaining == 0 {
                return Ok(None);
            }
            self.remaining -= 1;
            Ok(Some(locus(region(1, 1), Vec::new())))
        }
    }

    fn generator(per_segment: u32) -> GeneratorSlot<()> {
        GeneratorSlot::Generator(Box::new(FixedCountGenerator {
            per_segment,
            remaining: 0,
        }))
    }

    /// Each kind reaches **its own** slot (spec §13.1). Distinguishable generators (2 / 3 / 5
    /// loci per segment) make a mis-route show up as the wrong count, which indistinguishable
    /// `NoLoci` slots could not.
    #[test]
    fn each_kind_routes_to_its_own_slot() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let mut set = GeneratorSet::new(
            GeneratorSlot::Generator(Box::new(FixedCountGenerator {
                per_segment: 2,
                remaining: 0,
            })),
            generator(3),
            generator(5),
        );

        assert_eq!(
            drain_region(&mut set, typed(an_ssr_segment(20, 25), 20, 25), &reads).len(),
            2,
            "SsrSegment → the ssr slot"
        );
        assert_eq!(
            drain_region(&mut set, typed(RegionKind::Generic, 1, 10), &reads).len(),
            3,
            "Generic → the generic slot"
        );
        assert_eq!(
            drain_region(&mut set, typed(a_bundle(), 100, 160), &reads).len(),
            5,
            "SsrBundle → the bundle slot"
        );
        assert_eq!(
            drain_region(&mut set, typed(RegionKind::Satellite, 200, 400), &reads).len(),
            0,
            "Satellite has no slot"
        );

        let counts = set.counts();
        assert_eq!(counts.regions_in, 4);
        assert_eq!(counts.regions_handled, 3);
        assert_eq!(counts.loci_emitted, 2 + 3 + 5);
        assert_eq!(counts.unhandled_out_of_scope, 1);
        assert_eq!(counts.unhandled_not_implemented, 0);
        assert_eq!(
            counts.regions_in,
            counts.regions_handled
                + counts.unhandled_not_implemented
                + counts.unhandled_out_of_scope,
        );
    }

    /// A filled slot's region is *handled*: its loci are counted and it never touches the
    /// unhandled counters — the other side of the dispatch from the NoLoci case.
    #[test]
    fn a_filled_slot_counts_its_loci_and_is_not_unhandled() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let mut set = GeneratorSet::new(
            GeneratorSlot::Generator(Box::new(FixedCountGenerator {
                per_segment: 2,
                remaining: 0,
            })),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        );

        let loci = drain_region(&mut set, typed(an_ssr_segment(20, 25), 20, 25), &reads);
        assert_eq!(loci.len(), 2, "the generator's two loci per segment");

        let counts = set.counts();
        assert_eq!(counts.regions_in, 1);
        assert_eq!(counts.regions_handled, 1);
        assert_eq!(counts.loci_emitted, 2);
        assert_eq!(
            counts.unhandled_not_implemented, 0,
            "a handled region is not unhandled"
        );
        assert_eq!(counts.unhandled_out_of_scope, 0);
    }

    /// A generator that emits `per_segment` loci per segment, each carrying the segment's own
    /// coordinates — so the output echoes region order (for the order check) and a region
    /// yielding *several* loci is exercised through the iterator (spec §13.4).
    struct EchoGenerator {
        per_segment: u32,
        remaining: u32,
        region: Option<GenomeRegion>,
    }

    impl EchoGenerator {
        fn new(per_segment: u32) -> Self {
            Self {
                per_segment,
                remaining: 0,
                region: None,
            }
        }
    }

    impl<S> LocusGenerator<S> for EchoGenerator {
        fn begin_segment(&mut self, region: GenomeRegion) {
            self.remaining = self.per_segment;
            self.region = Some(region);
        }

        fn next_locus(
            &mut self,
            _segment: &S,
            _reads: &SampleReads,
        ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
            if self.remaining == 0 {
                return Ok(None);
            }
            self.remaining -= 1;
            Ok(Some(locus(
                self.region.expect("begin_segment ran first"),
                Vec::new(),
            )))
        }
    }

    fn echo_slot(per_segment: u32) -> GeneratorSlot<()> {
        GeneratorSlot::Generator(Box::new(EchoGenerator::new(per_segment)))
    }

    /// A generator that emits one locus, then fails — so a fatal generator error *after* a
    /// locus has already been yielded is exercised (distinct from the upstream-stream error).
    struct FailAfterOneGenerator {
        emitted: bool,
    }

    impl<S> LocusGenerator<S> for FailAfterOneGenerator {
        fn begin_segment(&mut self, _region: GenomeRegion) {
            self.emitted = false;
        }

        fn next_locus(
            &mut self,
            _segment: &S,
            _reads: &SampleReads,
        ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
            if !self.emitted {
                self.emitted = true;
                return Ok(Some(locus(region(1, 1), Vec::new())));
            }
            Err(LocusGenerationError::Reads(IngestError::NoFiles))
        }
    }

    /// The iterator drains a multi-kind stream and accounts every region, yielding nothing
    /// when all kinds are unimplemented — the shape run on its own (spec §2, §13.2, §13.3).
    #[test]
    fn the_iterator_drains_and_accounts_a_multi_kind_stream() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let regions = vec![
            Ok(typed(RegionKind::Generic, 1, 10)),
            Ok(typed(an_ssr_segment(20, 25), 20, 25)),
            Ok(typed(RegionKind::Satellite, 200, 400)),
        ];
        let mut iterator = SampleLocusObservationsIterator::new(
            regions.into_iter(),
            reads,
            GeneratorSet::all_unimplemented(),
        );

        assert!(
            iterator.next().is_none(),
            "an unimplemented run emits no loci"
        );

        let counts = iterator.counts();
        assert_eq!(counts.regions_in, 3);
        assert_eq!(counts.unhandled_not_implemented, 2);
        assert_eq!(counts.unhandled_out_of_scope, 1);
        assert_eq!(
            counts.regions_in,
            counts.regions_handled
                + counts.unhandled_not_implemented
                + counts.unhandled_out_of_scope,
        );
    }

    /// Emitted loci are in coordinate order across a multi-region, multi-kind stream — the
    /// output-order contract (spec §2, §13.4).
    #[test]
    fn loci_come_out_in_coordinate_order_across_kinds() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let regions = vec![
            Ok(typed(RegionKind::Generic, 1, 10)),
            Ok(typed(an_ssr_segment(20, 25), 20, 25)),
            Ok(typed(a_bundle(), 100, 160)),
        ];
        let set = GeneratorSet::new(
            GeneratorSlot::Generator(Box::new(EchoGenerator::new(1))),
            echo_slot(1),
            echo_slot(1),
        );
        let iterator = SampleLocusObservationsIterator::new(regions.into_iter(), reads, set);

        let starts: Vec<u64> = iterator
            .map(|item| item.unwrap().region.start.get())
            .collect();
        assert_eq!(
            starts,
            vec![1, 20, 100],
            "one locus per region, in the stream's coordinate order"
        );
    }

    /// A region yielding **several** loci streams them all, in order, before the next region
    /// — the iterator's "keep driving the same region across successive polls" branch, which
    /// spec §13.4 names explicitly.
    #[test]
    fn a_region_yielding_several_loci_streams_them_all_before_advancing() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let regions = vec![
            Ok(typed(RegionKind::Generic, 1, 10)), // generic slot → 3 loci at start=1
            Ok(typed(an_ssr_segment(20, 25), 20, 25)), // ssr slot → 2 loci at start=20
        ];
        let set = GeneratorSet::new(
            GeneratorSlot::Generator(Box::new(EchoGenerator::new(2))),
            echo_slot(3),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        );
        let iterator = SampleLocusObservationsIterator::new(regions.into_iter(), reads, set);

        let starts: Vec<u64> = iterator
            .map(|item| item.unwrap().region.start.get())
            .collect();
        // The first region's 3 loci in full, then the second's 2 — none dropped, none early.
        assert_eq!(starts, vec![1, 1, 1, 20, 20]);
    }

    /// A fatal generator error — after a locus has already been yielded — is surfaced once,
    /// then the iterator fuses. Distinct from the upstream-stream error path (spec §6).
    #[test]
    fn a_generator_error_mid_region_is_fatal_and_fuses() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let regions = vec![
            Ok(typed(an_ssr_segment(20, 25), 20, 25)),
            Ok(typed(RegionKind::Generic, 1, 10)),
        ];
        let set = GeneratorSet::new(
            GeneratorSlot::Generator(Box::new(FailAfterOneGenerator { emitted: false })),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
            GeneratorSlot::Unfilled(UnhandledReason::NotImplemented),
        );
        let mut iterator = SampleLocusObservationsIterator::new(regions.into_iter(), reads, set);

        assert!(
            matches!(iterator.next(), Some(Ok(_))),
            "the one locus before the failure"
        );
        match iterator.next() {
            Some(Err(LocusGenerationError::Reads(_))) => {}
            other => panic!("expected a fatal generator error, got {other:?}"),
        }
        assert!(
            iterator.next().is_none(),
            "fused after the generator error — the second region is never reached"
        );
    }

    /// A fatal upstream error is yielded once, wrapped, then the iterator is done — a failure
    /// never looks like clean end-of-stream, and the iterator is fused (spec §6).
    #[test]
    fn a_stream_error_is_fatal_and_the_iterator_fuses() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let regions = vec![
            Ok(typed(RegionKind::Generic, 1, 10)),
            Err(TypedRegionError::MarginNarrowerThanFlank {
                max_str_len: 1,
                bundle_threshold: 2,
            }),
            Ok(typed(RegionKind::Generic, 20, 30)),
        ];
        let mut iterator = SampleLocusObservationsIterator::new(
            regions.into_iter(),
            reads,
            GeneratorSet::all_unimplemented(),
        );

        match iterator.next() {
            Some(Err(LocusGenerationError::TypedRegion(_))) => {}
            other => panic!("expected a fatal wrapped TypedRegion error, got {other:?}"),
        }
        assert!(
            iterator.next().is_none(),
            "fused: nothing after the fatal error"
        );
        assert!(iterator.next().is_none(), "still fused on a repeated poll");
    }

    /// The types compose into a locus of each kind — a smoke test that the shared
    /// shape holds together before the contract and dispatcher land on it.
    #[test]
    fn a_locus_of_each_kind_can_be_built() {
        let generic = SampleLocusObservations {
            region: region(100, 100),
            reference_bases: Box::from(&b"A"[..]),
            observed_sequences: vec![ObservedSequence {
                bases: Box::from(&b"T"[..]),
                read_coverage: ReadCoverage::Complete,
                read_group: ReadGroupId(0),
                num_obs: 9,
                num_fwd: 5,
                q_sum: -12.0,
                mapq_sum: 540,
                mapq_sum_sq: 32_400,
                placed_left: 3,
                chain_ids: vec![1, 2],
            }],
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        };
        assert_eq!(generic.region.len(), 1);
        assert_eq!(generic.observed_sequences[0].num_obs, 9);

        let ssr = SampleLocusObservations {
            region: region(10_442, 10_461),
            reference_bases: Box::from(&b"ATATATATATATATATATAT"[..]),
            observed_sequences: Vec::new(),
            reads_without_observation: 3,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(SsrDetail {
                motif: Motif::new(b"AT").unwrap(),
                left_flank: Box::from(&b"CCCGGG"[..]),
                right_flank: Box::from(&b"TTTAAA"[..]),
            }),
        };
        // Zero coverage is a real observation: the locus exists with an empty table.
        assert!(ssr.observed_sequences.is_empty());
        assert!(matches!(ssr.kind, LocusKind::Ssr(_)));

        let bundle = SampleLocusObservations {
            region: region(200, 260),
            reference_bases: Box::from(&b"N"[..]),
            observed_sequences: Vec::new(),
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::SsrBundle,
        };
        assert_eq!(bundle.kind, LocusKind::SsrBundle);
    }

    /// A complete and a partial observation of the *same* bases are distinct evidence,
    /// so they must not compare equal — the property the `(bases, read_coverage)`
    /// dedup key rests on (spec §3).
    #[test]
    fn same_bases_differ_by_read_coverage() {
        let complete = ObservedSequence {
            bases: Box::from(&b"ATATAT"[..]),
            read_coverage: ReadCoverage::Complete,
            read_group: ReadGroupId(0),
            num_obs: 1,
            num_fwd: 1,
            q_sum: 0.0,
            mapq_sum: 60,
            mapq_sum_sq: 3_600,
            placed_left: 0,
            chain_ids: Vec::new(),
        };
        let partial = ObservedSequence {
            read_coverage: ReadCoverage::from_left(6, LocusLen::from_positions(6)),
            ..complete.clone()
        };
        assert_ne!(complete, partial);
        assert_ne!(complete.read_coverage, partial.read_coverage);
    }

    /// Depth derives correctly from read-coverage (spec §13.5): the vector has
    /// `region.len()` entries; a `Complete` raises every position, a left-flush run of
    /// `n` only the leftmost `n`, a right-flush run of `n` only the rightmost `n`. A
    /// 10-position locus with one complete (×3), one left-partial reaching 4 (×2), one
    /// right-partial reaching 3 (×5).
    #[test]
    fn depth_derives_from_read_coverage() {
        let l = locus(
            region(1, 10),
            vec![
                obs(b"AAAAAAAAAA", ReadCoverage::Complete, 3),
                obs(
                    b"AAAA",
                    ReadCoverage::from_left(4, LocusLen::from_positions(10)),
                    2,
                ),
                obs(
                    b"AAA",
                    ReadCoverage::from_right(3, LocusLen::from_positions(10)),
                    5,
                ),
            ],
        );
        // positions:            1  2  3  4  5  6  7  8  9 10
        //   complete ×3:        3  3  3  3  3  3  3  3  3  3
        //   left(4)  ×2:       +2 +2 +2 +2  .  .  .  .  .  .
        //   right(3) ×5:        .  .  .  .  .  .  . +5 +5 +5
        assert_eq!(l.num_obs_along_locus(), vec![5, 5, 5, 5, 3, 3, 3, 8, 8, 8],);
    }

    /// A single-base locus (a candidate SNP) has one depth position, raised by every
    /// complete observation over it.
    #[test]
    fn single_base_locus_has_one_depth_position() {
        let l = locus(
            region(42, 42),
            vec![
                obs(b"A", ReadCoverage::Complete, 7),
                obs(b"T", ReadCoverage::Complete, 2),
            ],
        );
        assert_eq!(l.num_obs_along_locus(), vec![9]);
    }

    /// No observations → depth is zero at every position, still `region.len()` long
    /// (the zero-coverage locus is a real one, not an absent one).
    #[test]
    fn no_observations_is_all_zero_full_length() {
        assert_eq!(
            locus(region(1, 4), Vec::new()).num_obs_along_locus(),
            vec![0; 4]
        );
    }

    /// A run claiming to reach further than the locus is long is clamped, not an
    /// out-of-bounds index — the consumer-side guard, which also survives an *unclamped*
    /// producer: an unclamped `locus_len - positions_covered` wraps to a huge offset, and
    /// this test still yields a bounded window rather than panicking.
    ///
    /// **It runs one case, not two.** Before the reshape `PartialLeft(9)` and
    /// `PartialRight(9)` on a 3-position locus were distinct values and exercised two arms;
    /// they now denote the same run, which is what
    /// `from_left_and_from_right_agree_once_the_reach_covers_the_whole_locus` states
    /// directly. The doc used to claim "both ends" and no longer does.
    #[test]
    fn a_run_reaching_beyond_the_locus_is_clamped() {
        let clamped = locus(
            region(1, 3),
            vec![obs(
                b"AAA",
                ReadCoverage::from_left(9, LocusLen::from_positions(3)),
                4,
            )],
        );
        assert_eq!(clamped.num_obs_along_locus(), vec![4, 4, 4]);
    }

    /// **The constructors place the run against their own border**, and the offset is
    /// derived from the *clamped* reach, never the raw one.
    ///
    /// The right case is the one that matters: `from_right(4, 10)` must start at 6. Deriving
    /// the offset before clamping would wrap on an over-long reach, and a saturating
    /// subtraction would quietly relabel a right-anchored read as left-anchored — the
    /// silent-depth failure this step was flagged for.
    #[test]
    fn from_right_places_the_run_against_the_right_border() {
        assert_eq!(
            ReadCoverage::from_left(4, LocusLen::from_positions(10)),
            ReadCoverage::Observed {
                offset_in_locus: 0,
                positions_covered: 4
            }
        );
        assert_eq!(
            ReadCoverage::from_right(4, LocusLen::from_positions(10)),
            ReadCoverage::Observed {
                offset_in_locus: 6,
                positions_covered: 4
            }
        );
    }

    /// **Once the reach covers the whole locus the two constructors agree** — a read that
    /// witnessed every position is constrained from neither border, so there is one run and
    /// not two.
    ///
    /// Stated as a test because it is a real behaviour change, not a curiosity: on the STR
    /// path an expanded allele can give a reach longer than the reference tract, and a
    /// left-anchored and a right-anchored read of the *same bases* then land in the **same
    /// tally cell** and merge into one row — where `PartialLeft(n)` and `PartialRight(n)`
    /// kept them apart. See `ssr::tally::tests::an_expanded_allele_merges_the_two_sides`.
    #[test]
    fn from_left_and_from_right_agree_once_the_reach_covers_the_whole_locus() {
        assert_eq!(
            ReadCoverage::from_left(9, LocusLen::from_positions(3)),
            ReadCoverage::Observed {
                offset_in_locus: 0,
                positions_covered: 3
            }
        );
        assert_eq!(
            ReadCoverage::from_right(9, LocusLen::from_positions(3)),
            ReadCoverage::from_left(9, LocusLen::from_positions(3))
        );
    }

    /// The flushness predicates — the **entire** surviving representation of
    /// prefix-versus-suffix, since the reshape dropped the side-tagged variants.
    ///
    /// Both `Complete` arms are asserted here because no call site reaches them: every
    /// `coverage_label` matches `Complete` first, so inverting either arm to `false` would
    /// otherwise change nothing anywhere in the tree.
    #[test]
    fn flushness_is_derived_from_where_the_run_sits() {
        assert!(ReadCoverage::Complete.is_flush_left());
        assert!(ReadCoverage::Complete.is_flush_right(LocusLen::from_positions(10)));

        let left = ReadCoverage::from_left(4, LocusLen::from_positions(10));
        assert!(left.is_flush_left(), "a prefix constraint");
        assert!(!left.is_flush_right(LocusLen::from_positions(10)));

        let right = ReadCoverage::from_right(4, LocusLen::from_positions(10));
        assert!(!right.is_flush_left());
        assert!(
            right.is_flush_right(LocusLen::from_positions(10)),
            "a suffix constraint"
        );

        // An interior run — flush with neither border. The STR path cannot mint one, but the
        // predicates are shared and the generic path will.
        let interior = ReadCoverage::Observed {
            offset_in_locus: 3,
            positions_covered: 4,
        };
        assert!(!interior.is_flush_left());
        assert!(!interior.is_flush_right(LocusLen::from_positions(10)));
    }

    /// Depth over an **interior** run — flush with neither border. This is the case the
    /// reshape exists to represent (a read blind in the middle of a footprint: an interior
    /// `N`, a ref-skip), and the only one where a wrong window neither panics nor clamps, so
    /// it is also the only one where an off-by-one in the arm is purely silent.
    #[test]
    fn depth_over_an_interior_run_raises_only_the_witnessed_stretch() {
        let l = locus(
            region(1, 10),
            vec![obs(
                b"AAAA",
                ReadCoverage::Observed {
                    offset_in_locus: 3,
                    positions_covered: 4,
                },
                7,
            )],
        );
        // positions: 1  2  3  4  5  6  7  8  9 10
        //            .  .  .  7  7  7  7  .  .  .
        assert_eq!(l.num_obs_along_locus(), vec![0, 0, 0, 7, 7, 7, 7, 0, 0, 0]);
    }

    /// `complete_observations()` yields only the complete entries — the guard that a
    /// partial (a lower bound) is never scored as an exact allele.
    #[test]
    fn complete_observations_excludes_partials() {
        let l = locus(
            region(1, 6),
            vec![
                obs(b"ATATAT", ReadCoverage::Complete, 4),
                obs(
                    b"ATATAT",
                    ReadCoverage::from_left(6, LocusLen::from_positions(6)),
                    2,
                ),
                obs(b"ATGTAT", ReadCoverage::Complete, 3),
                obs(
                    b"ATAT",
                    ReadCoverage::from_right(4, LocusLen::from_positions(6)),
                    1,
                ),
            ],
        );
        // Both completes, and only the completes — a partial is never scored as exact.
        let complete: Vec<&[u8]> = l
            .complete_observations()
            .map(|o| o.bases.as_ref())
            .collect();
        assert_eq!(complete, vec![&b"ATATAT"[..], &b"ATGTAT"[..]]);
        assert_eq!(l.complete_observations().count(), 2);
    }
}
