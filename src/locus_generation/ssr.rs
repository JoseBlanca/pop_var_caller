//! The STR locus generator — one microsatellite tract → one locus.
//!
//! A [`LocusGenerator`]: it consumes an `SsrSegment`, fetches the reads over the tract, aligns
//! each to read off its repeat, and tallies the answers into a [`SampleLocusObservations`]. It
//! was ported from the STR pileup of the caller deleted on 2026-09-13, adapted at two seams (the
//! split of coordinates from bases, and a widened admission gate). See
//! `doc/devel/ng/spec/locus_generation_ssr.md` (design) and
//! `doc/devel/ng/arch/locus_generation_ssr.md` (types & interfaces).
//!
//! This file holds the config, the counts, the working input ([`SsrLocus`] and its margin
//! fetch) and the public surface, [`SsrGenerator`]. The four steps it runs per locus live in
//! submodules:
//!
//! - [`read_cap`] — fetch the reads over the tract and its margin, and cap their number;
//! - [`read_region`] — map each read's CIGAR to the slice of the read over the locus window;
//! - [`classify`] — delimit the read's repeat with the chosen [`RepeatDelimiter`] (by default
//!   the unit-slip aligner; the flat-gap aligner for STR bundles and the parity oracle) and
//!   turn the result into an observation or a reason for none;
//! - [`tally`] — fold the kept reads' observations into the locus's observed sequences.
//!
//! [`slide_reach`] measures how far an indel inside a repeat can slide past its edges, and
//! [`bundle`] holds the STR bundle generator, which runs the same steps over a bundle's hull.

use super::{
    CursorLedger, GeneratorCounts, LocusGenerationError, LocusGenerator, LocusKind,
    SampleLocusObservations, SequenceObservation, SsrDetail,
};
use crate::alignment::ssr_unit_robust::SsrUnitRobustAligner;
use crate::alignment::{
    BestPathAligner, PerQualityEmission, RepeatContext, RepeatSpan, StutterModel,
};
use crate::read::filtering::ReadFilterCounts;
use crate::read::input::SampleReads;
#[cfg(test)]
use crate::read::input::cursor::CursorCounts;
use crate::read::input::sample_cursor::SampleCursor;
use crate::ref_seq::{ContigTable, EvictableRefSeq, RawRefSeq, RefSeq, RefSeqError};
use crate::region_typing::segment_criteria::SsrSegment;
use crate::types::{Bp, ContigId, GenomeRegion, Motif, Position, ReadGroupId};

pub mod bundle;

/// The repeat ground one STR locus measures, as reading a read across it needs it: where it
/// lies, the motif the delimiter is handed, the period of the reference run reaching each of
/// its edges, and whether a read that ends inside it is kept.
///
/// **An STR segment and an STR bundle's hull are both one of these** (`str_bundle_loci.md`
/// §6.1): the bundle generator is the segment's generator run over the hull, and these four
/// facts are everything in which the two differ while a read is measured.
#[derive(Debug, Clone, PartialEq)]
struct MeasuredRepeat {
    /// The repeat's first base, 1-based.
    start: u64,
    /// The repeat's last base, 1-based inclusive.
    end: u64,
    /// The motif [`RepeatGeometry`](crate::alignment::RepeatGeometry) hands the delimiter. A
    /// segment's own; for a bundle, its first member's, which the flat-gap delimiter it is
    /// measured with does not read.
    delimiter_motif: Motif,
    /// The period of the reference run reaching each edge, `None` where fixed sequence does.
    edge_run_periods: EdgeRunPeriods,
    /// Whether this is an STR segment or an STR bundle's hull — which decides whether a flank
    /// deletion the delimiter's own spelling implies is recorded
    /// (`classify::flank_deletion_the_delimiter_implies`).
    ground: MeasuredGround,
}

/// Which kind of repeat ground a [`MeasuredRepeat`] is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MeasuredGround {
    /// An STR segment, measured with its own delimiter.
    Segment,
    /// An STR bundle's hull, measured with the flat-gap delimiter.
    BundleHull,
}

/// The period of the reference run that reaches each edge of a measured repeat — what decides
/// whether the reference's own repeat carries on past that edge
/// (`reference_run_escapes_the_tract`). `None` where fixed sequence reaches the edge, which no
/// run continues.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct EdgeRunPeriods {
    pub left: Option<usize>,
    pub right: Option<usize>,
}

impl MeasuredRepeat {
    /// An STR bundle's hull, `start..=end` (1-based), measured with `delimiter_motif` in the
    /// geometry and `runs` deciding its edges. A read ending inside it is kept as a partial
    /// observation, as at a segment, so that it counts toward the hull's depth; calling scores it as
    /// a lower bound from the flank it held (`str_bundle_loci.md` §6.1, §7.6; changed during
    /// implementation 2026-10-11, and scored since 2026-10-12).
    pub fn bundle_hull(
        start: u64,
        end: u64,
        delimiter_motif: Motif,
        edge_run_periods: EdgeRunPeriods,
    ) -> Self {
        Self {
            start,
            end,
            delimiter_motif,
            edge_run_periods,
            ground: MeasuredGround::BundleHull,
        }
    }

    /// The repeat's first base, 1-based.
    pub fn start(&self) -> u64 {
        self.start
    }

    /// The repeat's last base, 1-based inclusive.
    pub fn end(&self) -> u64 {
        self.end
    }

    /// The repeat's length in reference bases.
    pub fn tract_len(&self) -> u64 {
        self.end - self.start + 1
    }
}

impl From<SsrSegment> for MeasuredRepeat {
    /// A segment is one run of one motif reaching both its edges, and a read ending inside it
    /// is a partial observation.
    fn from(segment: SsrSegment) -> Self {
        Self::from(&segment)
    }
}

impl From<&SsrSegment> for MeasuredRepeat {
    fn from(segment: &SsrSegment) -> Self {
        let period = Some(segment.period());
        Self {
            start: segment.start(),
            end: segment.end(),
            delimiter_motif: segment.motif(),
            edge_run_periods: EdgeRunPeriods {
                left: period,
                right: period,
            },
            ground: MeasuredGround::Segment,
        }
    }
}

/// An STR locus ready to align against: the repeat plus the reference bases the aligner
/// aligns the reads to.
///
/// The ng counterpart of production's `Locus`, **split** so the coordinates come from region
/// typing (`SsrSegment`, or an STR bundle's hull) and the bases are fetched here (spec §2).
/// Fetching them makes the port an *adaptation, not a lift* — the most likely place for the
/// port to go subtly wrong, which is why the flank lengths are always **measured** from the
/// clamped span, never assumed to be `flank_bp`.
#[derive(Debug, Clone, PartialEq)]
struct SsrLocus {
    /// The repeat measured: its coordinates and what measuring a read across it needs.
    pub repeat: MeasuredRepeat,
    /// The tract plus its query margin, canonical `{A,C,G,T,N}` bases, **clamped at contig
    /// ends** — so this may be shorter than `2 * flank_bp + tract`, and each flank must be
    /// measured, never assumed (spec §2).
    pub tract_with_margin_bases: Box<[u8]>,
    /// 1-based position of `tract_with_margin_bases[0]`.
    pub margin_start: Position,
}

impl SsrLocus {
    /// Fetch the tract ± `flank_bp` for `repeat` (an `SsrSegment`, or a bundle's hull) into an
    /// [`SsrLocus`], clamped at the contig ends, using the reused `buffer` as scratch (spec §2).
    ///
    /// `contig` is the tract's contig (from the region handed to `begin_segment`); its length
    /// is read from `reference`'s contig table to clamp the right margin, and `fetch_into`
    /// validates the final window. The bases are canonical — the aligner compares against
    /// them, and production upper-cases the same way.
    pub fn fetch<R: RefSeq + ContigTable>(
        reference: &R,
        contig: ContigId,
        repeat: impl Into<MeasuredRepeat>,
        flank_bp: Bp,
        buffer: &mut Vec<u8>,
    ) -> Result<Self, RefSeqError> {
        let repeat = repeat.into();
        let flank = flank_bp.get();
        let tract_start = repeat.start(); // 1-based inclusive
        let tract_end = repeat.end();
        // The left margin cannot run below base 1; the right cannot run past the contig.
        let margin_start = tract_start.saturating_sub(flank).max(1);
        let margin_end = match reference.contigs().entries.get(contig.get() as usize) {
            // Clamp the right margin to the contig — but only when the tract actually fits
            // inside it. A tract reaching past the contig end is a broken segment: leaving
            // the window unclamped makes `fetch_into` reject it as out-of-bounds, rather than
            // fetching a window that is missing the tract (which would underflow the derived
            // right flank). An unknown contig (`None`) is likewise left for `fetch_into`.
            Some(entry) if tract_end <= entry.length => (tract_end + flank).min(entry.length),
            _ => tract_end + flank,
        };
        let length = (margin_end + 1).saturating_sub(margin_start);

        reference.fetch_into(contig, margin_start, length, buffer)?;
        Ok(Self {
            repeat,
            tract_with_margin_bases: buffer.as_slice().into(),
            margin_start: Position(margin_start),
        })
    }

    /// The left flank's length in bases — **measured** from the clamp (`tract_start −
    /// margin_start`), so a tract near the contig start reports a short flank, not `flank_bp`.
    pub fn left_flank_len(&self) -> usize {
        (self.repeat.start() - self.margin_start.get()) as usize
    }

    /// The right flank's length in bases — measured as what remains after the left flank and
    /// the tract, so a tract near the contig end reports a short flank.
    pub fn right_flank_len(&self) -> usize {
        self.tract_with_margin_bases.len()
            - self.left_flank_len()
            - self.repeat.tract_len() as usize
    }
}

/// The STR generator's knobs — owned and taken at construction
/// ([shared config discipline](super); spec §4). Each generator owns its own knobs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrGeneratorConfig {
    /// The flank fetched either side of the tract — the aligner's anchor and the read
    /// query margin. Must be `<= bundle_threshold`, the radius region typing guarantees
    /// repeat-free (checked by [`Self::check_flank_within`]); equal by default (spec §4).
    pub flank_bp: Bp,
    /// Reads kept per locus, reservoir-sampled. `None` = no cap (spec §4).
    pub max_reads_per_locus: Option<u32>,
}

/// ng's own per-locus read cap — **not** production's `MAX_READS_PER_LOCUS`, so the two can
/// diverge. Starts at 1000 (matching production) but is never-measured and soft: to be set
/// by experiment (spec §4).
pub const DEFAULT_SSR_MAX_READS_PER_LOCUS: u32 = 1000;

impl Default for SsrGeneratorConfig {
    /// The flank equals the bundle threshold's default (spec §4, "equal by default"), and
    /// the cap starts at [`DEFAULT_SSR_MAX_READS_PER_LOCUS`].
    fn default() -> Self {
        Self {
            flank_bp: Bp(crate::region_typing::segment_criteria::DEFAULT_BUNDLE_THRESHOLD),
            max_reads_per_locus: Some(DEFAULT_SSR_MAX_READS_PER_LOCUS),
        }
    }
}

impl SsrGeneratorConfig {
    /// Check the cross-config invariant `flank_bp <= bundle_threshold`.
    ///
    /// A wider flank than the radius region typing guarantees repeat-free would let the read
    /// query hit a neighbouring repeat, leaving the aligner's anchor no longer clean (spec
    /// §4). It is a relation between two configs, so no newtype can hold it — the generator's
    /// constructor calls this. `bundle_threshold` is
    /// [`SsrSegmentCriteria::bundle_threshold`](crate::region_typing::segment_criteria::SsrSegmentCriteria::bundle_threshold)
    /// ([`crate::region_typing::segment_criteria`]).
    pub fn check_flank_within(&self, bundle_threshold: Bp) -> Result<(), SsrGeneratorConfigError> {
        if self.flank_bp.get() > bundle_threshold.get() {
            return Err(SsrGeneratorConfigError::FlankExceedsBundleThreshold {
                flank_bp: self.flank_bp.get(),
                bundle_threshold: bundle_threshold.get(),
            });
        }
        Ok(())
    }
}

/// A malformed STR generator configuration. `#[non_exhaustive]`; raised at construction so a
/// bad knob combination never reaches a locus.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum SsrGeneratorConfigError {
    /// The flank is wider than region typing's clean-flank guarantee — the read query could
    /// reach a neighbouring repeat and spoil the aligner's anchor (spec §4).
    #[error(
        "flank_bp ({flank_bp}) exceeds bundle_threshold ({bundle_threshold}): the read query \
         would reach past the repeat-free radius region typing guarantees"
    )]
    FlankExceedsBundleThreshold {
        flank_bp: u64,
        bundle_threshold: u64,
    },
}

/// Run-level STR counts, reported alongside the shared
/// [`LocusCounts`](super::LocusCounts). The locus records *that* reads yielded nothing
/// ([`reads_without_observation`](super::SampleLocusObservations::reads_without_observation));
/// **why** is this generator's to report, because the reasons are specific to how it reads a
/// tract and mean nothing to a pileup (spec §4).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct SsrGeneratorCounts {
    /// Reads fetched over the tract-plus-margin query span. Each is counted once more, in
    /// exactly one of [`outside_tract_before_cap`](Self::outside_tract_before_cap),
    /// [`reads_discarded_by_cap`](Self::reads_discarded_by_cap) or the per-read outcomes that
    /// follow them.
    pub reads_fetched: u64,
    /// Reads whose own alignment never reaches the tract, set aside **before** the read cap so
    /// that they do not take a place in it (`classify::lies_outside_the_tract`). Before this
    /// counter existed these reads were offered to the cap and, when kept, counted under
    /// [`outside_tract`](Self::outside_tract).
    pub outside_tract_before_cap: u64,
    /// Reads the reservoir cap dropped, **out of the reads that reach the tract**: the reads
    /// counted in [`outside_tract_before_cap`](Self::outside_tract_before_cap) are never
    /// offered to it.
    pub reads_discarded_by_cap: u64,
    /// Observations pinning the tract length (both borders seen).
    pub observations_complete: u64,
    /// Observations proving a lower bound (one border, ran off the other end).
    pub observations_partial: u64,
    /// Reads that reached the aligner and anchored no border (`read_preparation_ssr.md` §4).
    pub no_border_anchored: u64,
    /// Reads dropped for low base quality.
    pub low_quality: u64,
    /// Reads whose allele stayed flank-truncated even after widening.
    pub window_truncated: u64,
    /// Reads that reached the aligner, anchored a flank and crossed **no** tract position —
    /// they overlap the locus *window* and lie outside the repeat the locus is. Logged rather
    /// than merely dropped, because it is a large population and the number is the evidence
    /// that the fetch window is wider than the locus: **6,704 reads against 7,085 real
    /// partials** on tomato chr01 of `SRR7279503` when this counter was added. That figure
    /// included the reads whose CIGAR alone shows them outside, which are now set aside before
    /// the cap and counted in [`outside_tract_before_cap`](Self::outside_tract_before_cap)
    /// instead; this counter keeps the ones only the delimiter can find.
    pub outside_tract: u64,
}

impl SsrGeneratorCounts {}

mod read_cap;

use read_cap::{fetch_capped_reads, seed_for_segment};

mod read_region;

mod slide_reach;

/// How far before a repeat's start, and after its end, the slide measurement reads the
/// reference: the furthest a slide is followed, plus one indel unit and the base beside it.
const SLIDE_FETCH_MARGIN_BP: u64 =
    (slide_reach::MAX_SLIDE + slide_reach::MAX_INDEL_UNIT + 1) as u64;

mod classify;

mod tally;

// The generator: one `SsrSegment` → one `SampleLocusObservations`. The four steps (build the
// locus, fetch and cap the reads, classify each, tally) run inside the first `next_locus`; the
// second returns `None`. One locus per segment, including zero coverage (spec §2, arch §1, §2).

/// A tract delimiter the generator can drive — any alignment-module aligner that measures a read's
/// repeat: [algorithm 3](crate::alignment::ssr_best_path_flat_gap::SsrFlatGapAligner) (the
/// flat-gap port of production's `delimit_read`, kept as the byte-parity oracle and baseline),
/// algorithm 4 (the unit-slip
/// model), or [algorithm 4u](crate::alignment::ssr_unit_robust::SsrUnitRobustAligner) (unit-slip
/// hardened by the delimiter bake-off — **the recommended default**). Any [`BestPathAligner`] whose
/// `Output` is a [`RepeatSpan`] over a [`RepeatContext`] qualifies.
///
/// It is a **trait alias** (a supertrait bundle with a blanket impl), the one bound the generator
/// and [`classify::classify_read`] repeat — so the aligner is a **type parameter**, chosen per
/// generator and monomorphised. The choice is thus **static dispatch**: `align` is a direct call in
/// the per-read loop, never a `dyn` indirection (a bake-off runs millions of reads). The only
/// dynamic dispatch is the one-per-run `Box<dyn LocusGenerator>` the dispatcher already uses to hold
/// a generator — outside the hot loop.
pub trait RepeatDelimiter:
    BestPathAligner<Output = RepeatSpan> + for<'a> BestPathAligner<Context<'a> = RepeatContext<'a>>
{
}

impl<A> RepeatDelimiter for A where
    A: BestPathAligner<Output = RepeatSpan>
        + for<'a> BestPathAligner<Context<'a> = RepeatContext<'a>>
{
}

/// The STR locus generator: turns one microsatellite tract into one locus.
///
/// Holds its own accessors and reusable scratch, the "a generator holds its own accessors"
/// convention ([`LocusGenerator`]; spec §2). The delimiter `A` is a **type
/// parameter** ([`RepeatDelimiter`]) so the algorithm is chosen per generator and dispatched
/// statically — [`with_default_aligner`](Self::with_default_aligner) builds the recommended
/// algorithm 4u (the unit-robust aligner — algorithm 4 hardened by the delimiter bake-off), while
/// [`new`](Self::new) takes any aligner, which is how a bake-off swaps in algorithm 3 (the flat-gap
/// port) and how the parity oracle pins it to unit-slip / flat-gap (spec §6).
///
/// **The reference seam (the Arc gap).** Two reference handles, because they serve two needs the
/// current [`RawRefSeq`] design keeps apart:
/// - `reference` — the margin fetch ([`SsrLocus::fetch`], [`RefSeq`] + [`ContigTable`]), a
///   persistent accessor the generator holds.
/// - `make_reference` — a **factory** the cursor needs ([`SampleReads::cursor`]): each file's
///   cursor owns its own reader, and a `RawRefSeq` is a stateful reader that is neither `Clone`
///   nor shareable by `&`, so the cursor takes `FnMut() -> R` rather than a borrow.
///
/// These are the same reference logically; a future `Arc`-shared reference would collapse them to
/// one field with a cheap per-file view. Holding both is the working stopgap the spec flags (spec
/// §8) — **and the per-chromosome cursor took most of its cost away**: the factory is called once per file per
/// *chromosome* now, where it used to run at every locus, so a file-backed `R` whose factory
/// reloads is no longer the per-locus tax it was.
pub struct SsrGenerator<R: RawRefSeq + EvictableRefSeq + Send + 'static, A: RepeatDelimiter> {
    /// The reference the margin fetch reads the flanks from.
    reference: R,
    /// Builds a fresh read reference for the cursor's mismatch-fraction filter — one per file,
    /// per chromosome (see the type doc).
    ///
    /// **Boxed, which is what keeps it off this type's parameter list**, the same trade the
    /// generic generator took (arch §3.6). It is called at chromosome boundaries and
    /// nowhere else, so the indirection is not on any path that matters.
    ///
    /// **`+ Send` since 2026-09-01, because [`GeneratorSlot`](super::GeneratorSlot) now
    /// requires it of everything boxed into a slot** — a calling run's walkers cross threads under the merge's parallel
    /// cover, so a generator that could fill a slot must be able to go with them. The caller
    /// that used to rule the bound out — `ng_ssr_cohort_stutter` handing every file a clone of
    /// one `Arc<WindowedRefSeq>` — stopped being a counter-example in the same step, when
    /// `WindowedRefSeq` became `Sync` and that `Arc` became `Send`.
    make_reference: Box<dyn FnMut() -> R + Send>,
    /// The sample's cursor and the chromosome it covers — one per chromosome, minted by
    /// the first locus on each and rebuilt at every boundary, because nothing in a cursor
    /// survives a chromosome change (`spec/alignment_cursor.md` §4).
    ///
    /// `None` before the first locus. It is *not* cleared per segment: the whole point is that
    /// it outlives them.
    cursor: Option<(ContigId, SampleCursor<R>)>,
    /// **The sample this generator reads, and what its retired chromosomes' cursors did** —
    /// one sample per generator, refused in [`cursor_for`](Self::cursor_for), and the tallies
    /// `cursor_counts` (tests only) and [`read_filter_counts`](Self::read_filter_counts)
    /// add the live cursor's to.
    ledger: CursorLedger,
    /// The tract delimiter — the chosen [`RepeatDelimiter`] (algorithm 3 or 4).
    aligner: A,
    /// Reused alignment matrices (the aligner's own scratch type), so it does not reallocate per
    /// read.
    align_scratch: A::Scratch,
    /// The stutter model feeding the aligner's context. Algorithm 3 (flat-gap) ignores it;
    /// algorithm 4 (unit-slip) prices its slips from it. The context carries it either way (spec
    /// §2, arch §5).
    stutter: StutterModel,
    /// Reused scratch for the margin fetch.
    margin_buffer: Vec<u8>,
    /// Reused scratch for the reference around a tract that [`slide_reach`] measures.
    slide_buffer: Vec<u8>,
    /// Reused scratch for the quality gate.
    qual_buffer: Vec<u8>,
    config: SsrGeneratorConfig,
    counts: SsrGeneratorCounts,
    /// The region of the segment begun, from [`begin_segment`](LocusGenerator::begin_segment) —
    /// its `ContigId` is what the margin fetch and the read query key on (`SsrSegment` carries a
    /// contig *name*, not an id).
    current_region: Option<GenomeRegion>,
    /// Latched once the segment's single locus has been produced, so the second `next_locus`
    /// returns `None` (spec §2).
    produced: bool,
}

impl<R> SsrGenerator<R, SsrUnitRobustAligner<PerQualityEmission>>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq + Send + 'static,
{
    /// A generator with the **recommended** delimiter — algorithm 4u, the *unit-robust* aligner over
    /// production's [`PerQualityEmission`] table. It is algorithm 4 (unit-slip) hardened by the
    /// delimiter bake-off with a narrow junction guard (a sequencing error near a flank/tract
    /// boundary can no longer slide the boundary) and an evidence-based anchor test (a complete is
    /// only reported when the flank was actually matched, else a lower-bound partial), the demotion
    /// capped so a partial is never lost. On HG002 it removes ~4,500 fabricated exact-lengths per
    /// chr20–22 that unit-slip reported and preserves every real partial (`ng_ssr_gain_loss`).
    ///
    /// The default for real use; a bake-off, or the production byte-parity oracle, passes a specific
    /// aligner to [`new`](Self::new) — parity is checked against unit-slip / the flat-gap port, not
    /// against this, because unit-robust deliberately departs from production's measurement.
    pub fn with_default_aligner(
        reference: R,
        make_reference: impl FnMut() -> R + Send + 'static,
        config: SsrGeneratorConfig,
        bundle_threshold: Bp,
    ) -> Result<Self, SsrGeneratorConfigError> {
        Self::new(
            reference,
            make_reference,
            SsrUnitRobustAligner::new(PerQualityEmission::new()),
            config,
            bundle_threshold,
        )
    }
}

impl<R, A> SsrGenerator<R, A>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq + Send + 'static,
    A: RepeatDelimiter,
{
    /// Build a generator over `reference` (the margin fetch), `make_reference` (the read-query
    /// factory) and `aligner` (the chosen delimiter — algorithm 3 or 4), with `config` checked
    /// against `bundle_threshold`, the region-typing radius the flank must stay within (spec §4).
    /// Fails if `flank_bp > bundle_threshold`. For the recommended default, use
    /// [`with_default_aligner`](Self::with_default_aligner).
    pub fn new(
        reference: R,
        make_reference: impl FnMut() -> R + Send + 'static,
        aligner: A,
        config: SsrGeneratorConfig,
        bundle_threshold: Bp,
    ) -> Result<Self, SsrGeneratorConfigError> {
        config.check_flank_within(bundle_threshold)?;
        Ok(Self {
            reference,
            make_reference: Box::new(make_reference),
            cursor: None,
            ledger: CursorLedger::default(),
            aligner,
            align_scratch: A::Scratch::default(),
            stutter: StutterModel::hipstr_shipped(),
            margin_buffer: Vec::new(),
            slide_buffer: Vec::new(),
            qual_buffer: Vec::new(),
            config,
            counts: SsrGeneratorCounts::default(),
            current_region: None,
            produced: false,
        })
    }

    /// The running STR counts — accumulated across every segment, readable at any point.
    pub fn counts(&self) -> &SsrGeneratorCounts {
        &self.counts
    }

    #[cfg(test)]
    /// **What the cursors did**, summed over every chromosome walked so far — see
    /// [`PileupGenerator::cursor_counts`](super::pileup::PileupGenerator::cursor_counts) for the
    /// argument, which is the same one.
    ///
    /// In short: a generator that mints a cursor per locus and one that keeps a cursor per
    /// chromosome emit **identical loci**, so nothing in the output can tell them apart, and a
    /// review proved the feature could be switched off with the whole suite green. What
    /// moves is only what the reader avoided (`spec/alignment_cursor.md` §11.5).
    pub fn cursor_counts(&self) -> CursorCounts {
        self.ledger
            .cursor_counts(self.cursor.as_ref().map(|(_, cursor)| cursor))
    }

    /// **Why this sample's reads were dropped, per read group, over the whole walk** — every
    /// chromosome, not merely the one the cursor is standing on.
    ///
    /// `run_streaming.md` §8 requires a run to sum these when it finishes, and names what
    /// skipping it costs: *"drop rates under-report … silently, since every number stays
    /// plausible"* — every figure the run prints is a number a real run could have produced,
    /// so nothing about the output says the tally is short. Retired chromosomes are added up as
    /// they go and the live one is asked directly, so the answer is current at any moment.
    ///
    /// **A drop rate is a read group's property and the entries are never summed together**: a
    /// run that went badly shows up as one library with an anomalous rate, and adding them
    /// erases exactly that signal.
    ///
    /// **The `None` key is reads that declared no read group**, which is the key the cursor
    /// itself uses; a run whose files all declare one never sees it.
    pub fn read_filter_counts(&self) -> Vec<(Option<ReadGroupId>, ReadFilterCounts)> {
        self.ledger
            .read_filter_counts(self.cursor.as_ref().map(|(_, cursor)| cursor))
    }

    /// This sample's cursor for `contig`, minting one if the last locus was on another
    /// chromosome.
    ///
    /// The cursor is what replaces a read query per locus. It stays open and positioned, so two
    /// STR loci close enough to share reads decode them once
    /// (`spec/alignment_cursor.md` §1). A chromosome change is the one thing it cannot absorb:
    /// nothing in a cursor survives it — the kept reads are useless and, on CRAM, so are the
    /// reference bases — so the boundary mints a new one rather than repositioning the old (§4).
    ///
    /// **The STR generator's own regions are far apart**, so this reuses far less than the
    /// generic generator's tiling walk does — and the ratio the Checkpoint C audit measured
    /// collapses with stride (23× dense, 3.8× at a 13× stride, 0.5× walked backwards). So the
    /// gain here is **measured, not assumed**: `ng_ssr_loci_dump` on chromosome 21 of HG002 30×
    /// goes from **17.3 s to 12.5 s of user CPU, −28 %**, with the dump byte-identical. Wall time
    /// does not move, because that tool's wall is the background whole-genome reference MD5.
    ///
    /// The two generators hold **separate** cursors on purpose: their regions interleave, and
    /// sharing one would tie their lifetimes together (spec §12).
    ///
    /// # One sample per generator, and it is checked here
    ///
    /// A cursor is opened for **one sample's files** and kept. Handing this generator a second
    /// sample would answer that sample out of the first one's files — every individual in a
    /// cohort reporting the first one's reads, with no error and output of exactly the right
    /// shape. So the sample is remembered alongside the chromosome and a foreign one is refused
    /// ([`ForeignSample`](LocusGenerationError::ForeignSample)).
    ///
    /// **Refused rather than re-opened**, and the arithmetic is why. A cohort tool walks
    /// region-major — every sample is asked about one repeat before anything moves to the next
    /// repeat — so re-opening on a change of sample would open one reader per sample *per
    /// locus*, which is worse than the per-locus query this replaced. A generator per sample
    /// opens one reader per sample per chromosome, which is what the design assumes anyway
    /// (`arch/alignment_cursor.md` §2.4 counts open files as `files × generators × workers`).
    fn cursor_for(
        &mut self,
        contig: ContigId,
        reads: &SampleReads,
        region: GenomeRegion,
    ) -> Result<&mut SampleCursor<R>, LocusGenerationError> {
        let sample = self.ledger.check_sample(reads, region)?;
        // Let both readers drop what this walk has gone past. A reference reader is a
        // sliding window over the FASTA that only ever grows unless it is told to shrink, and
        // nobody was telling it: a forward walk of a chromosome ended up holding one byte for
        // every base it had passed. The cursor's margin is one flank, which is the furthest back
        // it looks from a repeat's own start. Releasing is a hint: releasing too much costs a
        // re-read and never an answer.
        //
        // **The generator's own reader keeps more: [`SLIDE_FETCH_MARGIN_BP`], 107 bases.**
        // `slide_reach` reads that far back from the repeat's start right after this call, so
        // releasing from one flank (15 bases by default) made every STR locus re-read about 92
        // bases it had just let go of — one extra file read a locus a sample.
        let keep_from = region.start.get().saturating_sub(self.config.flank_bp.0);
        let reference_keep_from = region
            .start
            .get()
            .saturating_sub(self.config.flank_bp.0.max(SLIDE_FETCH_MARGIN_BP));
        self.reference.evict_before(reference_keep_from);
        if let Some((_, cursor)) = &self.cursor {
            cursor.evict_reference_before(keep_from);
        }

        if self.cursor.as_ref().is_none_or(|(on, _)| *on != contig) {
            // Taken before the new one is built, so a run never holds two chromosomes'
            // cursors — and their kept reads — at once. The tallies are taken with it: they
            // die with the cursor, and they are what says whether it did anything.
            if let Some((_, retiring)) = self.cursor.take() {
                self.ledger.retire(&retiring);
            }
            let cursor = reads
                .cursor(contig, &mut self.make_reference)
                .map_err(|source| LocusGenerationError::OpenReadQuery { region, source })?;
            self.cursor = Some((contig, cursor));
            // Adopted only once a cursor exists, so a failed open leaves the generator
            // unclaimed and a retry with the same sample still works.
            self.ledger.adopt(sample);
        }
        Ok(&mut self.cursor.as_mut().expect("the cursor was just ensured").1)
    }
}

impl<R, A> LocusGenerator<SsrSegment> for SsrGenerator<R, A>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq + Send + 'static,
    A: RepeatDelimiter,
{
    fn begin_segment(&mut self, region: GenomeRegion) {
        self.current_region = Some(region);
        self.produced = false;
    }

    /// The STR counts, reachable through a boxed generator — the same need the generic
    /// generator's nine had. This tool-driven path reads them off the concrete type
    /// ([`SsrGenerator::counts`]); a dispatcher-driven one cannot, and now does not have
    /// to.
    fn counts(&self) -> Option<GeneratorCounts<'_>> {
        Some(GeneratorCounts::Ssr(SsrGenerator::counts(self)))
    }

    /// **The only way the per-read-group drop tallies are reachable once the generator is
    /// boxed** — `GeneratorSlot` erases the type, so without this a run with the tract slot
    /// filled would report only the SNP/indel generator's drops and call them the sample's.
    fn read_filter_counts(&self) -> Vec<(Option<ReadGroupId>, ReadFilterCounts)> {
        SsrGenerator::read_filter_counts(self)
    }

    fn next_locus(
        &mut self,
        segment: &SsrSegment,
        reads: &SampleReads,
    ) -> Result<Option<SampleLocusObservations>, LocusGenerationError> {
        let Some(region) = self.claim_the_one_locus() else {
            return Ok(None);
        };
        let contig = region.contig;

        // The margin fetch and the read query key on `region.contig` (an id), the reservoir seed
        // on `segment.chrom()` (a name); the `LocusGenerator` contract leaves the region/segment
        // pairing unenforced, so a mis-pairing would silently fetch the wrong contig. Guard it in
        // debug/test builds against the reference's own name for that id.
        debug_assert_eq!(
            self.contig_name(contig),
            Some(segment.chrom()),
            "the region's contig and the segment must name the same contig"
        );

        // (1) Build the locus — fetch the tract ± flank, clamped at contig ends.
        let locus = self.fetch_locus(region, segment)?;

        // (2)–(4) Fetch, cap, classify and tally the reads.
        let measured = self.measure_reads(&locus, region, seed_for_segment(segment), reads)?;

        // Assemble the output: the tract coordinates, the tract bases only (flanks go in the
        // kind), and the flanks split out of the fetched margin (spec §3).
        let (left_flank, reference_bases, right_flank) = locus.split_flanks();
        Ok(Some(SampleLocusObservations {
            chain_id_base: 0,
            region: GenomeRegion {
                contig,
                start: Position(segment.start()),
                end: Position(segment.end()),
            },
            reference_bases: reference_bases.into(),
            observations: measured.observations,
            reads_without_observation: measured.reads_without_observation,
            reads_discarded_by_cap: measured.reads_discarded_by_cap,
            kind: LocusKind::Ssr(SsrDetail {
                motif: segment.motif(),
                left_flank: left_flank.into(),
                right_flank: right_flank.into(),
            }),
        }))
    }
}

/// What one locus's reads came to: the tally, and how many reads the depth cap took away.
struct MeasuredReads {
    observations: Vec<SequenceObservation>,
    reads_without_observation: u32,
    reads_discarded_by_cap: u32,
}

impl SsrLocus {
    /// The fetched bases as `(left flank, repeat, right flank)`.
    fn split_flanks(&self) -> (&[u8], &[u8], &[u8]) {
        let left = self.left_flank_len();
        let tract_len = self.repeat.tract_len() as usize;
        let bases = &*self.tract_with_margin_bases;
        (
            &bases[..left],
            &bases[left..left + tract_len],
            &bases[left + tract_len..],
        )
    }
}

/// The steps of one locus that an STR segment and an STR bundle share — the bundle generator
/// ([`bundle`]) is this generator run over a hull (`str_bundle_loci.md` §6.1), and
/// these are what it runs.
impl<R, A> SsrGenerator<R, A>
where
    R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq + Send + 'static,
    A: RepeatDelimiter,
{
    /// The region begun, the first time it is asked for and `None` after — a region yields one
    /// locus (spec §2).
    ///
    /// # Panics
    ///
    /// When no region was begun: the `LocusGenerator` contract calls `begin_segment` first.
    fn claim_the_one_locus(&mut self) -> Option<GenomeRegion> {
        if self.produced {
            return None;
        }
        self.produced = true;
        Some(
            self.current_region
                .expect("begin_segment is called before next_locus (the LocusGenerator contract)"),
        )
    }

    /// The reference's name for `contig`, which the reservoir seed is taken over.
    fn contig_name(&self, contig: ContigId) -> Option<&str> {
        self.reference
            .contigs()
            .entries
            .get(contig.get() as usize)
            .map(|entry| entry.name.as_str())
    }

    /// Step (1): the repeat ± flank, clamped at the contig ends.
    fn fetch_locus(
        &mut self,
        region: GenomeRegion,
        repeat: impl Into<MeasuredRepeat>,
    ) -> Result<SsrLocus, LocusGenerationError> {
        SsrLocus::fetch(
            &self.reference,
            region.contig,
            repeat,
            self.config.flank_bp,
            &mut self.margin_buffer,
        )
        .map_err(|source| LocusGenerationError::Reference { region, source })
    }

    /// How far an indel in `locus`'s tract can slide past its edges, measured on the reference
    /// ([`slide_reach::SlideReach`]).
    fn slide_reach(
        &mut self,
        locus: &SsrLocus,
        contig: ContigId,
        region: GenomeRegion,
    ) -> Result<slide_reach::SlideReach, LocusGenerationError> {
        let margin = SLIDE_FETCH_MARGIN_BP;
        let start = locus.repeat.start().saturating_sub(margin).max(1);
        let contig_len = self
            .reference
            .contigs()
            .entries
            .get(contig.get() as usize)
            .map_or(locus.repeat.end(), |entry| entry.length);
        let end = (locus.repeat.end() + margin).min(contig_len);
        self.reference
            .fetch_into(contig, start, end + 1 - start, &mut self.slide_buffer)
            .map_err(|source| LocusGenerationError::Reference { region, source })?;
        self.slide_buffer.make_ascii_uppercase();
        let first = (locus.repeat.start() - start) as usize;
        let last = (locus.repeat.end() - start) as usize;
        Ok(slide_reach::SlideReach::measure(
            &self.slide_buffer,
            first,
            last,
        ))
    }

    /// Steps (2) to (4): fetch the reads over the repeat-plus-margin query span, admitting on
    /// relevance (overlap, which `SampleReads` applies) — not spanning — and depth-cap them;
    /// classify each against `locus` (delimit, recover, gate); tally the outcomes into the
    /// locus's observed sequences and the run-level counts.
    fn measure_reads(
        &mut self,
        locus: &SsrLocus,
        region: GenomeRegion,
        seed: u64,
        reads: &SampleReads,
    ) -> Result<MeasuredReads, LocusGenerationError> {
        let contig = region.contig;
        let margin_start = locus.margin_start.get();
        // `margin_len >= 1` on a successful `fetch` (the tract is non-empty), so the `- 1` for
        // the inclusive end cannot underflow.
        let margin_len = locus.tract_with_margin_bases.len() as u64;
        let query_span = GenomeRegion {
            contig,
            start: Position(margin_start),
            end: Position(margin_start + margin_len - 1),
        };
        // Read out before the cursor is borrowed: `cursor_for` takes `&mut self`, and a config
        // field read in the same call would hold a shared borrow across it.
        let max_reads_per_locus = self.config.max_reads_per_locus;
        // **Reads that never reach the tract are set aside before the cap**, from their CIGAR
        // alone, so that the cap's places go to reads that can measure the tract. They are
        // still reads that covered the locus and produced no observation, so they are added to
        // `reads_without_observation` below, as they were when the cap kept them.
        let capped = fetch_capped_reads(
            self.cursor_for(contig, reads, region)?,
            query_span,
            seed,
            max_reads_per_locus,
            |read| classify::lies_outside_the_tract(read, locus),
        )?;
        // Per-locus discards cannot approach 2^32 (the cap and any real depth are far below), so
        // the `as u32` on the locus fields below is exact; the run-level counters keep the u64.
        let reads_discarded_by_cap = capped.fetched - capped.set_aside - capped.kept.len() as u64;
        self.counts.reads_fetched += capped.fetched;
        self.counts.outside_tract_before_cap += capped.set_aside;
        self.counts.reads_discarded_by_cap += reads_discarded_by_cap;

        let reach = self.slide_reach(locus, contig, region)?;
        let mut outcomes = Vec::with_capacity(capped.kept.len());
        for read in &capped.kept {
            let outcome = classify::classify_read(
                read,
                locus,
                &self.aligner,
                &self.stutter,
                &mut self.align_scratch,
                &mut self.qual_buffer,
            );
            outcomes.push(slide_reach::demoted_if_short_of_the_slide(
                outcome, read, locus, reach,
            ));
        }

        let tallied = tally::tally(
            capped.kept.iter().zip(outcomes),
            locus.repeat.start(),
            locus.repeat.end(),
            &mut self.counts,
        );
        Ok(MeasuredReads {
            observations: tallied.observations,
            reads_without_observation: tallied.reads_without_observation + capped.set_aside as u32,
            reads_discarded_by_cap: reads_discarded_by_cap as u32,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::locus_generation::LocusLen;
    use crate::ref_seq::InMemoryRefSeq;
    use crate::region_typing::segment_criteria::Motif;
    use crate::types::ReadGroupId;

    /// A 100-base contig `chr1` with a known repeating pattern, so a fetched span can be
    /// checked byte-for-byte.
    fn reference_100() -> (InMemoryRefSeq, Vec<u8>) {
        let bases: Vec<u8> = b"ACGT".iter().cycle().take(100).copied().collect();
        (
            InMemoryRefSeq::from_named_contigs(vec![("chr1".to_string(), bases.clone())]),
            bases,
        )
    }

    fn tract(start: u64, end: u64) -> SsrSegment {
        SsrSegment::new("chr1".into(), start, end, Motif::new(b"AC").unwrap(), 1.0).unwrap()
    }

    /// A tract with room either side: both flanks are exactly `flank_bp`, the margin span is
    /// `tract + 2·flank_bp`, and the fetched bytes are exactly that slice of the reference.
    #[test]
    fn a_mid_contig_tract_fetches_equal_flanks_and_the_exact_bytes() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // tract [40, 49] (10 bases), flank 10 → margin [30, 59].
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(40, 49), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(30));
        assert_eq!(locus.left_flank_len(), 10);
        assert_eq!(locus.right_flank_len(), 10);
        assert_eq!(locus.tract_with_margin_bases.len(), 30); // 10 + 10 + 10
        // margin [30, 59] 1-based == bases[29..59] 0-based.
        assert_eq!(&*locus.tract_with_margin_bases, &bases[29..59]);
    }

    /// A tract near the contig **start**: the left margin clamps to base 1, so the left flank
    /// is short while the right is full — the flanks are unequal, measured not assumed.
    #[test]
    fn a_tract_near_the_contig_start_has_a_short_left_flank() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // tract [5, 14], flank 10 → margin_start clamps to 1, margin_end 24.
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(5, 14), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(1));
        assert_eq!(locus.left_flank_len(), 4, "clamped: 5 - 1");
        assert_eq!(locus.right_flank_len(), 10, "full");
        assert_ne!(locus.left_flank_len(), locus.right_flank_len());
        assert_eq!(&*locus.tract_with_margin_bases, &bases[0..24]);
    }

    /// A tract near the contig **end**: the right margin clamps to the contig length, so the
    /// right flank is short — again unequal.
    #[test]
    fn a_tract_near_the_contig_end_has_a_short_right_flank() {
        let (reference, bases) = reference_100();
        let mut buffer = Vec::new();
        // contig length 100; tract [92, 96], flank 10 → margin [82, 100].
        let locus =
            SsrLocus::fetch(&reference, ContigId(0), tract(92, 96), Bp(10), &mut buffer).unwrap();

        assert_eq!(locus.margin_start, Position(82));
        assert_eq!(locus.left_flank_len(), 10, "full");
        assert_eq!(locus.right_flank_len(), 4, "clamped: 100 - 96");
        assert_eq!(locus.tract_with_margin_bases.len(), 19); // 10 + 5 + 4
        assert_eq!(&*locus.tract_with_margin_bases, &bases[81..100]);
    }

    /// A tract reaching past the contig end is a broken segment: `fetch` leaves the window
    /// unclamped so `fetch_into` rejects it, rather than returning a locus whose derived
    /// right flank would underflow.
    #[test]
    fn a_tract_past_the_contig_end_is_rejected_not_silently_truncated() {
        let (reference, _bases) = reference_100();
        let mut buffer = Vec::new();
        // contig length 100; tract [98, 110] runs past the end.
        assert!(
            SsrLocus::fetch(&reference, ContigId(0), tract(98, 110), Bp(10), &mut buffer).is_err()
        );
    }

    #[test]
    fn default_flank_equals_the_bundle_threshold_default_and_caps_at_1000() {
        let config = SsrGeneratorConfig::default();
        assert_eq!(
            config.flank_bp,
            Bp(crate::region_typing::segment_criteria::DEFAULT_BUNDLE_THRESHOLD),
            "flank equals the bundle threshold by default (spec §4)"
        );
        assert_eq!(
            config.max_reads_per_locus,
            Some(DEFAULT_SSR_MAX_READS_PER_LOCUS)
        );
    }

    #[test]
    fn the_flank_check_accepts_within_and_at_the_threshold() {
        // Equal is allowed (the default case), and strictly narrower.
        assert!(
            SsrGeneratorConfig {
                flank_bp: Bp(30),
                max_reads_per_locus: None,
            }
            .check_flank_within(Bp(30))
            .is_ok()
        );
        assert!(
            SsrGeneratorConfig {
                flank_bp: Bp(20),
                max_reads_per_locus: None,
            }
            .check_flank_within(Bp(30))
            .is_ok()
        );
    }

    #[test]
    fn the_flank_check_rejects_a_flank_wider_than_the_bundle_threshold() {
        let error = SsrGeneratorConfig {
            flank_bp: Bp(50),
            max_reads_per_locus: None,
        }
        .check_flank_within(Bp(30))
        .expect_err("50 > 30 must be refused");
        assert!(matches!(
            error,
            SsrGeneratorConfigError::FlankExceedsBundleThreshold {
                flank_bp: 50,
                bundle_threshold: 30,
            }
        ));
    }

    // --- the per-locus read fetch and cap -----------------------------------

    fn span(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// A `RawRefSeq` over the fixture contigs (all `A`s), for the cursor's
    /// mismatch-fraction filter — matches the fixture reference the reads are opened against.
    /// Named (`chr1` / `chr2`) so a generator's contig-name invariant holds against the
    /// `SsrSegment`'s `chrom()`; the fetch itself keys on `ContigId`, so the names are only for
    /// that check.
    fn fixture_ref_bases() -> InMemoryRefSeq {
        use crate::read::input::test_fixtures::matching_contigs;
        InMemoryRefSeq::from_named_contigs(
            matching_contigs()
                .iter()
                .map(|(name, len, _)| ((*name).to_string(), vec![b'A'; *len]))
                .collect(),
        )
    }

    /// A cursor over the fixture sample's first contig — what
    /// [`fetch_capped_reads`] is pointed at, in place of the sample plus a per-locus query.
    fn a_cursor(reads: &SampleReads) -> SampleCursor<InMemoryRefSeq> {
        reads
            .cursor(ContigId(0), fixture_ref_bases)
            .expect("the fixture sample opens a cursor")
    }

    /// Open a `SampleReads` over one indexed BAM of `records`, against the fixture reference.
    fn sample_reads_with(
        records: &[noodles_sam::alignment::RecordBuf],
    ) -> (tempfile::TempDir, tempfile::TempDir, SampleReads) {
        use crate::read::filtering::ReadFilterConfig;
        use crate::read::input::test_fixtures::{
            fixture_reference, header, indexed_bam, matching_contigs,
        };
        let (reference_dir, reference) = fixture_reference(false);
        let (bam_dir, bam) = indexed_bam(
            &header(
                Some("coordinate"),
                &matching_contigs(),
                &[("rg1", Some("NA12878"))],
            ),
            records,
        );
        let reads =
            SampleReads::open_only_sample(&[bam], &reference, ReadFilterConfig::default(), false)
                .expect("the fixture sample opens");
        (reference_dir, bam_dir, reads)
    }

    /// The cap keeps `max_reads_per_locus` of the fetched reads and counts them all; `None`
    /// keeps everything (spec §4).
    #[test]
    fn fetch_caps_reads_and_counts_the_fetched() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..5)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * 10, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);

        let mut cursor = a_cursor(&reads);
        let capped = fetch_capped_reads(&mut cursor, span(1, 100), 42, Some(3), |_| false)
            .expect("fetch succeeds");
        assert_eq!(capped.fetched, 5);
        assert_eq!(capped.kept.len(), 3, "capped to 3 of the 5 fetched");

        let mut cursor = a_cursor(&reads);
        let uncapped = fetch_capped_reads(&mut cursor, span(1, 100), 42, None, |_| false)
            .expect("fetch succeeds");
        assert_eq!(uncapped.fetched, 5);
        assert_eq!(uncapped.kept.len(), 5, "no cap keeps every read");
    }

    /// The kept subset is a deterministic function of the seed (the thread-invariance the seed
    /// buys, spec §4).
    #[test]
    fn the_capped_kept_set_is_deterministic_for_a_seed() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..20)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 1 + i * 4, 20))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let kept_positions = |seed| {
            let mut positions: Vec<u64> =
                fetch_capped_reads(&mut a_cursor(&reads), span(1, 100), seed, Some(5), |_| {
                    false
                })
                .expect("fetch succeeds")
                .kept
                .iter()
                .map(|read| read.pos)
                .collect();
            positions.sort_unstable();
            positions
        };
        assert_eq!(
            kept_positions(7),
            kept_positions(7),
            "same seed → same kept set"
        );
    }

    // --- the generator --------------------------------------------------------

    /// A generator over the all-`A` fixture reference (`chr1` = 100 bp, `chr2` = 200 bp), with
    /// `flank_bp = 10` (within a bundle threshold of 10) and no cap. The reference and the
    /// read-query factory are the same all-`A` bases, so fixture reads (also all-`A`) clear the
    /// mismatch filter — the tract content is not real STR sequence, which is fine for the
    /// structural checks here (real tract parity is the parity oracle's).
    fn ssr_generator() -> SsrGenerator<InMemoryRefSeq, SsrUnitRobustAligner<PerQualityEmission>> {
        SsrGenerator::new(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrUnitRobustAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: None,
            },
            Bp(10),
        )
        .expect("flank within the bundle threshold")
    }

    /// **The cursor is kept across loci rather than minted for each — and the test that can
    /// tell.**
    ///
    /// A generator that opens a query per locus and one that keeps a cursor per chromosome emit
    /// **identical loci**; that is the correctness requirement, and it is why a review could
    /// switch the whole feature off with 1,557 tests green. The only observable is what the
    /// reader avoided (`spec/alignment_cursor.md` §11.5), which is what `cursor_counts` exists
    /// for.
    ///
    /// Asserted on `regions_reusing`, not on `reads_decoded`: a fresh cursor has no last region
    /// to compare against, so its first move always jumps — three loci through one cursor jump
    /// once and reuse twice, and three cursors jump three times.
    #[test]
    fn the_cursor_is_kept_across_loci_rather_than_minted_for_each() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("r0", 0, 20, 30),
            read_named_with_length("r1", 0, 40, 30),
            read_named_with_length("r2", 0, 55, 30),
        ]);
        let mut generator = ssr_generator();

        for (start, end) in [(30u64, 39u64), (45, 54), (60, 69)] {
            generator.begin_segment(span(start, end));
            generator
                .next_locus(&tract(start, end), &reads)
                .expect("no fetch error")
                .expect("one locus per segment");
        }

        let counts = generator.cursor_counts();
        assert_eq!(
            (counts.regions_jumping, counts.regions_reusing),
            (1, 2),
            "three ascending loci through one cursor jump once — at the first, which has \
             nothing to reuse — and reuse twice: {counts:?}",
        );
        assert_eq!(
            (counts.reads_decoded, counts.reads_replayed),
            (3, 4),
            "the three reads are decoded once each and handed back four more times — locus 2 \
             replays r0 and r1, locus 3 replays r1 and r2. That ratio *is* the feature: \
             {counts:?}",
        );
        assert_eq!(
            counts.reads_evicted, 1,
            "and r0 ends before the third locus begins, so it is dropped — the only witness \
             that 'the kept set is bounded' is an observation rather than an argument: \
             {counts:?}",
        );
    }

    /// **Every retired chromosome's tallies survive, not just the last one's.**
    ///
    /// With a single retirement, accumulating and *overwriting* are numerically identical — so
    /// `retired += counts` mutated to `retired = counts` passed, and so did dropping four of the
    /// five fields. Both matter: `regions_reusing` and `regions_jumping` are the only observable
    /// that says whether the cursor is being kept at all, and a fold that loses them would leave
    /// the feature detectable on the first chromosome and undetectable on every other.
    ///
    /// chr1 → chr2 → chr1 is two retirements, which is what tells the two apart.
    #[test]
    fn returning_to_a_chromosome_keeps_every_retired_cursors_tallies() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("chr1-a", 0, 20, 30),
            read_named_with_length("chr2-a", 1, 20, 30),
            read_named_with_length("chr2-b", 1, 25, 30),
        ]);
        let mut generator = ssr_generator();

        let mut decoded = Vec::new();
        for (contig, name) in [(0u32, "chr1"), (1, "chr2"), (0, "chr1")] {
            let segment =
                SsrSegment::new(name.into(), 30, 39, Motif::new(b"AC").unwrap(), 1.0).unwrap();
            generator.begin_segment(GenomeRegion {
                contig: ContigId(contig),
                start: Position(30),
                end: Position(39),
            });
            generator
                .next_locus(&segment, &reads)
                .expect("no fetch error")
                .expect("one locus per segment");
            decoded.push(generator.cursor_counts().reads_decoded);
        }

        assert_eq!(
            decoded,
            vec![1, 3, 4],
            "chr1 decodes one read, chr2 two more, and returning to chr1 decodes its one \
             again — a fold that overwrote instead of accumulating would report [1, 2, 1]",
        );
    }

    /// **The generator says why its reader dropped reads, per read group, over the whole walk**
    /// — and the drops from a chromosome it has left are still in the answer.
    ///
    /// **The failure this exists for is silent.** The tally defaulted to empty, so a run with
    /// the tract slot filled reported only the SNP/indel generator's drops and called them the
    /// sample's — every number plausible, nothing to notice, and a library that had gone wrong
    /// invisible. Reading only the *live* cursor is the same failure one chromosome at a time:
    /// a whole-genome walk would report the last chromosome's drops as the run's.
    ///
    /// Two duplicate reads and two kept ones, one pair on each chromosome, so both halves of
    /// the sum are exercised and neither could be produced by the other.
    #[test]
    fn the_generator_reports_its_readers_drops_per_read_group_across_chromosomes() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let duplicate = |name: &str, contig: usize| {
            let mut read = read_named_with_length(name, contig, 20, 30);
            *read.flags_mut() = noodles_sam::alignment::record::Flags::DUPLICATE;
            read
        };
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("kept-on-chr1", 0, 20, 30),
            duplicate("dup-on-chr1", 0),
            read_named_with_length("kept-on-chr2", 1, 20, 30),
            duplicate("dup-on-chr2", 1),
        ]);
        let mut generator = ssr_generator();
        let summed = |generator: &SsrGenerator<_, _>| {
            generator
                .read_filter_counts()
                .into_iter()
                .fold((0_u64, 0_u64), |(kept, duplicate), (_, counts)| {
                    (kept + counts.kept, duplicate + counts.duplicate)
                })
        };

        fn walk<R, A>(
            generator: &mut SsrGenerator<R, A>,
            reads: &SampleReads,
            contig: u32,
            name: &str,
        ) where
            R: RefSeq + ContigTable + RawRefSeq + EvictableRefSeq + Send + 'static,
            A: RepeatDelimiter,
        {
            let segment =
                SsrSegment::new(name.into(), 30, 39, Motif::new(b"AC").unwrap(), 1.0).unwrap();
            generator.begin_segment(GenomeRegion {
                contig: ContigId(contig),
                start: Position(30),
                end: Position(39),
            });
            generator
                .next_locus(&segment, reads)
                .expect("no fetch error")
                .expect("one locus per segment");
        }

        walk(&mut generator, &reads, 0, "chr1");
        let after_first = summed(&generator);
        walk(&mut generator, &reads, 1, "chr2");
        let after_second = summed(&generator);

        assert_eq!(
            after_first,
            (1, 1),
            "the first chromosome's one kept read and one duplicate, off its live cursor",
        );
        assert_eq!(
            after_second,
            (2, 2),
            "and after the boundary both chromosomes' counts, summed — not the second's alone, \
             which is what a generator that lost the retiring cursor would report",
        );

        // **And it is reachable through the trait**, which is the only way a run can ask: a
        // generator in a `GeneratorSlot` has its type erased, and the trait's default answer is
        // an empty list. Without the override every number above would be right and no run
        // could read one.
        let boxed: &dyn LocusGenerator<SsrSegment> = &generator;
        assert_eq!(
            boxed.read_filter_counts(),
            SsrGenerator::read_filter_counts(&generator),
        );
    }

    /// **The per-read-group tallies and the aggregate one count the same reads.**
    ///
    /// Two harvests of one cursor taken at the same instant, read by different callers — the
    /// aggregate by whoever asks whether the cursor did anything, the per-group by the run
    /// report. A walk whose two answers differ has lost one of them between the cursor and the
    /// report, and nothing else would say so.
    #[test]
    fn the_tract_generators_per_group_tallies_and_its_aggregate_count_the_same_reads() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("on-chr1", 0, 20, 30),
            read_named_with_length("on-chr2", 1, 20, 30),
        ]);
        let mut generator = ssr_generator();

        for (contig, name) in [(0u32, "chr1"), (1, "chr2")] {
            let segment =
                SsrSegment::new(name.into(), 30, 39, Motif::new(b"AC").unwrap(), 1.0).unwrap();
            generator.begin_segment(GenomeRegion {
                contig: ContigId(contig),
                start: Position(30),
                end: Position(39),
            });
            generator
                .next_locus(&segment, &reads)
                .expect("no fetch error")
                .expect("one locus per segment");
        }

        let kept: u64 = generator
            .read_filter_counts()
            .iter()
            .map(|(_, counts)| counts.kept)
            .sum();
        assert_eq!(
            kept,
            generator.cursor_counts().reads_decoded,
            "the reads the filters kept, per group, are the reads the cursor decoded",
        );
        assert!(
            kept > 0,
            "a walk that read nothing could not tell the two apart",
        );
    }

    /// **A locus on a new chromosome mints a cursor, and the old one's tallies survive it.**
    ///
    /// A cursor covers one chromosome and refuses a region on any other (§4), so the generator
    /// has to notice the boundary. When a query was opened per locus nothing needed to see
    /// this: each query carried its contig.
    #[test]
    fn a_locus_on_a_new_chromosome_mints_a_cursor_and_keeps_the_old_tallies() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&[
            read_named_with_length("on-chr1", 0, 20, 30),
            read_named_with_length("on-chr2", 1, 20, 30),
        ]);
        let mut generator = ssr_generator();

        generator.begin_segment(span(30, 39));
        let first = generator
            .next_locus(&tract(30, 39), &reads)
            .expect("no fetch error")
            .expect("one locus per segment");
        let after_first = generator.cursor_counts();

        // The same coordinates on chr2 — `SsrSegment` carries a contig *name*, the region a
        // `ContigId`, and the generator's own debug assertion pins that they agree.
        let on_chr2 = SsrSegment::new("chr2".into(), 30, 39, Motif::new(b"AC").unwrap(), 1.0)
            .expect("a valid segment");
        generator.begin_segment(GenomeRegion {
            contig: ContigId(1),
            start: Position(30),
            end: Position(39),
        });
        let second = generator
            .next_locus(&on_chr2, &reads)
            .expect("chr2 must be reachable — a cursor stuck on chr1 would refuse it")
            .expect("one locus per segment");
        let after_second = generator.cursor_counts();

        // Not asserted on `region.contig`: a locus's region is built from `begin_segment`'s,
        // never from the cursor, so it re-states the test's own input and cannot fail under any
        // mutation here (review). What is load-bearing is that chr2's *read* arrived.
        // **Arrived, not observed.** The fixture's contigs are all `A`, so an indel in this tract
        // could slide to either contig end, and no read can show past it: each is set aside as
        // a read without an observation (`slide_reach`). Arriving is what this test is about.
        let arrived = |locus: &SampleLocusObservations| {
            locus.observations.len() as u32 + locus.reads_without_observation
        };
        assert!(
            arrived(&first) > 0,
            "chr1's read must reach its locus, or the comparison below is between two nothings",
        );
        assert!(
            arrived(&second) > 0,
            "chr2's read must reach the locus, not merely fail to error",
        );
        assert!(
            after_first.reads_decoded > 0,
            "the first chromosome must decode something, or this test cannot fail: \
             {after_first:?}",
        );
        assert!(
            after_second.reads_decoded > after_first.reads_decoded,
            "the retiring cursor's tallies must be taken before it is dropped, not replaced: \
             {after_first:?} then {after_second:?}",
        );
    }

    /// **A second sample through one generator is refused, not answered.**
    ///
    /// This is the shape of a real defect, caught by review before it shipped. A generator
    /// opens a reader for one sample's files and keeps it for a whole chromosome, but the
    /// `LocusGenerator` trait hands it a `&SampleReads` afresh on every call. Keyed on the
    /// chromosome alone, it answered the second sample out of the **first** sample's files —
    /// so a cohort tool asking 51 plants about one repeat would have got one plant's reads 51
    /// times, under 51 names, with no error and rows of exactly the right shape.
    ///
    /// The fixture is built so the wrong answer is unmistakable: sample A has a read on the
    /// tract, and sample B's nearest read is 45 bases away — outside the queried span
    /// `[20,49]` — so *any* observation from B is A's.
    #[test]
    fn a_second_sample_through_one_generator_is_refused() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let (_ra, _ba, sample_a) = sample_reads_with(&[read_named_with_length("a0", 0, 25, 30)]);
        let (_rb, _bb, sample_b) = sample_reads_with(&[read_named_with_length("b0", 0, 80, 30)]);
        let mut generator = ssr_generator();
        let seg = tract(30, 39);

        generator.begin_segment(span(30, 39));
        let a = generator
            .next_locus(&seg, &sample_a)
            .expect("no fetch error")
            .expect("one locus per segment");
        assert_eq!(
            a.observations.len(),
            1,
            "sample A's read must reach the tract, or this test cannot fail",
        );

        generator.begin_segment(span(30, 39));
        let refused = generator.next_locus(&seg, &sample_b);

        match refused {
            Err(LocusGenerationError::ForeignSample { region }) => {
                assert_eq!(
                    region,
                    span(30, 39),
                    "the refusal names the region asked for"
                );
            }
            Ok(Some(locus)) => panic!(
                "sample B was answered instead of refused, and with sample A's read — B has \
                 nothing within 45 bases of this tract: {:?}",
                locus.observations,
            ),
            other => panic!("expected a refusal naming the foreign sample, got {other:?}"),
        }
    }

    /// **A generator per sample is the shape that works** — the other half of the test above,
    /// because "refuses everything" would pass that one on its own.
    #[test]
    fn a_generator_per_sample_gives_each_sample_its_own_reads() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let (_ra, _ba, sample_a) = sample_reads_with(&[read_named_with_length("a0", 0, 25, 30)]);
        let (_rb, _bb, sample_b) = sample_reads_with(&[read_named_with_length("b0", 0, 80, 30)]);
        let seg = tract(30, 39);

        let observations_for = |reads: &SampleReads| {
            let mut generator = ssr_generator();
            generator.begin_segment(span(30, 39));
            generator
                .next_locus(&seg, reads)
                .expect("no fetch error")
                .expect("one locus per segment")
                .observations
                .len()
        };

        assert_eq!(
            observations_for(&sample_a),
            1,
            "sample A's read covers the tract",
        );
        assert_eq!(
            observations_for(&sample_b),
            0,
            "sample B's read is 45 bases away and must not appear — which is the whole point",
        );
    }

    /// A reference that records what it was told to release, delegating everything else.
    ///
    /// The fixture reference is in memory and holds no window, so "how much is resident" reads
    /// zero however badly the release is wired. Asking *what was released, and when* is
    /// falsifiable where asking *how much is held* is not.
    struct ReleaseSpy {
        inner: InMemoryRefSeq,
        released: std::sync::Arc<std::sync::Mutex<Vec<u64>>>,
    }

    impl EvictableRefSeq for ReleaseSpy {
        fn evict_before(&self, pos: u64) {
            self.released.lock().expect("no panic holds this").push(pos);
        }
    }

    impl RawRefSeq for ReleaseSpy {
        fn fetch_raw_into(
            &self,
            contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            self.inner.fetch_raw_into(contig, start_1based, length, dst)
        }
    }

    impl RefSeq for ReleaseSpy {
        fn fetch_into(
            &self,
            contig: ContigId,
            start_1based: u64,
            length: u64,
            dst: &mut Vec<u8>,
        ) -> Result<(), RefSeqError> {
            self.inner.fetch_into(contig, start_1based, length, dst)
        }
    }

    impl ContigTable for ReleaseSpy {
        fn contigs(&self) -> &crate::fasta::ContigList {
            self.inner.contigs()
        }
    }

    /// **Every repeat tells the reference readers what they may release, and the mark moves
    /// forward.**
    ///
    /// A reference reader is a sliding window over the FASTA. It grows as it is asked for more
    /// and shrinks only when told to — so a walk that never tells it holds one byte for every
    /// base it has passed, about 250 MB on human chromosome 1. `ref_seq`'s own tests measure
    /// that directly; this pins the two things *this* generator is responsible for: that it
    /// releases at all, and that the mark advances with the walk.
    ///
    /// Reviewers found this exact call site deletable with the whole suite green, which is why
    /// it has a test of its own rather than sharing the generic generator's.
    #[test]
    fn every_repeat_releases_the_reference_behind_it() {
        use crate::read::input::test_fixtures::read_named_with_length;
        // On the 200-base fixture contig, `chr2`: the generator's reader keeps
        // `SLIDE_FETCH_MARGIN_BP` (107) bases behind each repeat, so the repeats must start past
        // 107 for its mark to move. Kept clear of the contig's end: the last read runs 150..179.
        let records: Vec<_> = (0..4)
            .map(|index| read_named_with_length(&format!("r{index}"), 1, 105 + index * 15, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);

        // Two recorders: the margin fetch's reader, and the one the read filter holds inside
        // the cursor.
        let released = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let released_by_cursor = std::sync::Arc::new(std::sync::Mutex::new(Vec::new()));
        let spy = ReleaseSpy {
            inner: fixture_ref_bases(),
            released: std::sync::Arc::clone(&released),
        };
        let make_reference = {
            let released_by_cursor = std::sync::Arc::clone(&released_by_cursor);
            move || ReleaseSpy {
                inner: fixture_ref_bases(),
                released: std::sync::Arc::clone(&released_by_cursor),
            }
        };
        let mut generator = SsrGenerator::new(
            spy,
            make_reference,
            SsrUnitRobustAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: None,
            },
            Bp(10),
        )
        .expect("flank within the bundle threshold");

        for index in 0..4u64 {
            let start = 120 + index * 10;
            generator.begin_segment(GenomeRegion {
                contig: ContigId(1),
                ..span(start, start + 9)
            });
            let on_chr2 = SsrSegment::new(
                "chr2".into(),
                start,
                start + 9,
                Motif::new(b"AC").unwrap(),
                1.0,
            )
            .unwrap();
            generator
                .next_locus(&on_chr2, &reads)
                .expect("no fetch error");
        }

        let marks = released.lock().expect("no panic holds this").clone();
        assert_eq!(
            marks.len(),
            4,
            "one release per repeat, so a walk of any length holds one repeat's worth: \
             {marks:?}",
        );
        assert!(
            marks.windows(2).all(|pair| pair[1] > pair[0]),
            "the mark must move forward with the walk, or it releases the same nothing every \
             time: {marks:?}",
        );
        // The first repeat starts at 120; this reader keeps the slide margin, 107 bases,
        // which is further back than the flank (10).
        assert_eq!(
            marks[0],
            120 - SLIDE_FETCH_MARGIN_BP,
            "the mark is the repeat's start less the slide margin, the furthest back this \
             reader looks from a repeat's own start: {marks:?}",
        );

        // The cursor's reader has none to release at the first repeat, which opens it, and it
        // keeps one flank: the furthest back it looks.
        let by_cursor = released_by_cursor
            .lock()
            .expect("no panic holds this")
            .clone();
        assert_eq!(
            by_cursor,
            vec![120, 130, 140],
            "the cursor's reader is released too, from the repeat after the one that opens it, \
             at that repeat's start less one flank: {by_cursor:?}",
        );
    }

    /// `new` refuses a flank wider than the bundle threshold — the cross-config check at
    /// construction (spec §4).
    #[test]
    fn new_rejects_a_flank_wider_than_the_bundle_threshold() {
        // `SsrGenerator` is not `Debug` (it holds a closure and the aligner), so match on the
        // `Result` directly rather than `expect_err`.
        let result = SsrGenerator::with_default_aligner(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrGeneratorConfig {
                flank_bp: Bp(50),
                max_reads_per_locus: None,
            },
            Bp(10),
        );
        assert!(matches!(
            result,
            Err(SsrGeneratorConfigError::FlankExceedsBundleThreshold { .. })
        ));
    }

    /// A tract no read covers still yields **one** locus — present, with an empty observation
    /// table and zeroed drop counts — and the second `next_locus` returns `None`. "We looked and
    /// saw nothing" ≠ "we never looked" (spec §2). The locus carries the tract coordinates, the
    /// tract bases only, and the flanks split out of the fetched margin.
    #[test]
    fn a_zero_coverage_tract_yields_one_empty_but_present_locus() {
        // One read on chr1 at 65..94 — clear of the tract's query span (30..59), and long
        // enough (30 bp) to clear the min-read-length filter, so it is a genuine miss rather
        // than a filtered read.
        use crate::read::input::test_fixtures::read_named_with_length;
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r0", 0, 65, 30)]);
        let mut generator = ssr_generator();
        let segment = tract(40, 49); // chr1, motif AC, 10 bp

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus per segment, even at zero coverage");

        assert_eq!(locus.region, span(40, 49));
        assert_eq!(locus.region.len(), 10);
        assert_eq!(
            &*locus.reference_bases,
            &b"AAAAAAAAAA"[..],
            "the tract only"
        );
        assert!(locus.observations.is_empty(), "no read covered it");
        assert_eq!(locus.reads_without_observation, 0);
        assert_eq!(locus.reads_discarded_by_cap, 0);
        match locus.kind {
            LocusKind::Ssr(SsrDetail {
                motif,
                left_flank,
                right_flank,
            }) => {
                assert_eq!(motif, Motif::new(b"AC").unwrap());
                assert_eq!(&*left_flank, &b"AAAAAAAAAA"[..], "10 bp left flank");
                assert_eq!(&*right_flank, &b"AAAAAAAAAA"[..], "10 bp right flank");
            }
            other => panic!("expected an Ssr locus, got {other:?}"),
        }

        assert!(
            generator
                .next_locus(&segment, &reads)
                .expect("no error")
                .is_none(),
            "the second poll ends the segment"
        );
        assert_eq!(
            generator.counts().reads_fetched,
            0,
            "no read reached the tract"
        );
    }

    /// A covered tract wires fetch → classify → tally: every fetched read is accounted for,
    /// either as an observation's support or in `reads_without_observation` (no cap, so nothing
    /// is discarded). The tract content here is not real STR sequence, so *which* bucket a read
    /// lands in is E's concern; this pins the conservation the four steps must uphold.
    #[test]
    fn a_covered_tract_accounts_every_fetched_read() {
        use crate::read::input::test_fixtures::read_named_with_length;
        // Four 30 bp reads (clearing the min-read-length filter) overlapping the query span
        // (30..59) on chr1.
        let records: Vec<_> = (0..4)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 30 + i * 3, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = ssr_generator();
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus");

        let fetched = generator.counts().reads_fetched;
        assert_eq!(fetched, 4, "all four overlapping reads reached the tract");
        let supported: u32 = locus.observations.iter().map(|obs| obs.num_obs).sum();
        assert_eq!(
            supported as u64 + locus.reads_without_observation as u64,
            fetched,
            "every fetched read is either an observation or a no-observation"
        );
        assert_eq!(locus.reads_discarded_by_cap, 0, "no cap");
    }

    /// **`placed_left` is counted against the whole tract, and this is the only test that says
    /// which coordinates those are.**
    ///
    /// The fold itself (`tally`) is handed the tract's ends as bare `u64`s, and its unit test
    /// supplies them itself — so the unit test pins the comparison and nothing about *which*
    /// coordinates the generator passes. A 0-based/1-based slip, the margin's start, or one end
    /// of the tract in place of both would all be silent: `placed_left` feeds a QUAL penalty, so
    /// a wrong coordinate is a wrong QUAL with no panic.
    ///
    /// The tract is `[40, 49]` and `flank_bp` is 10, so the margin starts at 30. The reads are
    /// 30 bases long and start at 28, 29, 31 and 32: the first two lie at least half left of the
    /// tract and the other two do not. Each wrong coordinate moves that boundary past some of
    /// them: against the tract's first base only a read starting at 25 or earlier counts,
    /// against the margin's first base 15 or earlier, against the tract's last base 34 or
    /// earlier. (Each read stops inside the margin, which over this all-A reference is what
    /// gives it an observation at all.)
    #[test]
    fn placed_left_is_counted_against_the_whole_tract() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = [28usize, 29, 31, 32]
            .iter()
            .enumerate()
            .map(|(i, start)| read_named_with_length(&format!("r{i}"), 0, *start, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = ssr_generator();
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus");

        let supporting: u32 = locus.observations.iter().map(|obs| obs.num_obs).sum();
        let placed_left: u32 = locus.observations.iter().map(|obs| obs.placed_left).sum();

        assert_eq!(
            supporting, 4,
            "every read must produce an observation, or this test cannot discriminate"
        );
        assert_eq!(
            placed_left, 2,
            "the reads starting at 28 and 29 lie at least half left of the tract [40, 49]; \
             counting against its first base alone, or the margin's, would give 0, against its \
             last base alone 4, and the old started-left rule 4"
        );
    }

    /// **The length the mint clamps against is the length the emitted locus reports.**
    ///
    /// They come from different expressions — the mint uses `segment.tract_len()`, because it runs
    /// per read before the locus exists, while every consumer asks the locus
    /// (`SampleLocusObservations::locus_len`, over `region`). They agree by construction, since the
    /// region *is* `segment.start()..=segment.end()`, but nothing said so and a divergence would be
    /// silent: a run clamped against one length and tested for flushness against another reports
    /// the wrong side.
    #[test]
    fn the_mints_locus_length_is_the_emitted_regions() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let records = vec![read_named_with_length("r0", 0, 30, 30)];
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = ssr_generator();
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no fetch error")
            .expect("one locus");

        assert_eq!(
            locus.locus_len(),
            LocusLen::from_positions(segment.tract_len()),
            "the consumer's length and the mint's are the same quantity"
        );
    }

    /// The cap wires through `next_locus`: with `max_reads_per_locus = 2` over four overlapping
    /// reads, the locus reports two discarded, and the run-level counter agrees (the `as u32`
    /// path and `counts.reads_discarded_by_cap`).
    #[test]
    fn a_capped_locus_reports_the_discarded_reads() {
        use crate::read::input::test_fixtures::read_named_with_length;
        let records: Vec<_> = (0..4)
            .map(|i| read_named_with_length(&format!("r{i}"), 0, 30 + i * 3, 30))
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = SsrGenerator::new(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrUnitRobustAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: Some(2),
            },
            Bp(10),
        )
        .expect("flank within the bundle threshold");
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no error")
            .expect("one locus");

        assert_eq!(generator.counts().reads_fetched, 4);
        assert_eq!(locus.reads_discarded_by_cap, 2, "4 fetched, 2 kept");
        assert_eq!(generator.counts().reads_discarded_by_cap, 2);
    }

    /// **Reads that never reach the tract do not take places in the read cap.** 1,500 reads
    /// are fetched at one tract, 300 of them lying wholly in the left flank, under a cap of
    /// 1,000: the cap keeps 1,000 reads that cross the tract. When the flank-only reads were
    /// offered to the cap too, about a fifth of the 1,000 places went to them (they are a fifth
    /// of the fetch) and were then discarded as outside the tract.
    ///
    /// The flank-only reads still count as reads that covered the locus and showed nothing,
    /// so the locus's reads add up to the 1,500 fetched as before.
    #[test]
    fn reads_outside_the_tract_are_set_aside_before_the_cap() {
        use crate::read::input::test_fixtures::read_named_with_length;
        // Tract [40, 49], flank 10: the fetch spans [30, 59]. A 30-base read at 30 crosses the
        // tract; one at 1 ends on base 30, in the left flank. In position order, as a sorted
        // file holds them.
        let records: Vec<_> = (0..1_500)
            .map(|i| {
                let start = if i < 300 { 1 } else { 30 };
                read_named_with_length(&format!("r{i}"), 0, start, 30)
            })
            .collect();
        let (_reference_dir, _bam_dir, reads) = sample_reads_with(&records);
        let mut generator = SsrGenerator::new(
            fixture_ref_bases(),
            fixture_ref_bases as fn() -> InMemoryRefSeq,
            SsrUnitRobustAligner::new(PerQualityEmission::new()),
            SsrGeneratorConfig {
                flank_bp: Bp(10),
                max_reads_per_locus: Some(1_000),
            },
            Bp(10),
        )
        .expect("flank within the bundle threshold");
        let segment = tract(40, 49);

        generator.begin_segment(span(40, 49));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no error")
            .expect("one locus");

        let counts = generator.counts().clone();
        let outside = counts.outside_tract + counts.outside_tract_before_cap;
        let observed: u64 = locus
            .observations
            .iter()
            .map(|o| u64::from(o.num_obs))
            .sum();
        let kept_crossing = observed + u64::from(locus.reads_without_observation) - outside;
        assert_eq!(counts.reads_fetched, 1_500);
        assert_eq!(
            kept_crossing, 1_000,
            "the cap's 1,000 places all go to reads crossing the tract: {counts:?}"
        );
        assert_eq!(
            outside, 300,
            "every flank-only read is found outside: {counts:?}"
        );
        assert_eq!(
            counts.outside_tract_before_cap, 300,
            "set aside before the cap"
        );
        assert_eq!(
            locus.reads_discarded_by_cap, 200,
            "1,200 crossing, 1,000 kept"
        );
        assert_eq!(counts.reads_discarded_by_cap, 200);
        assert_eq!(
            observed
                + u64::from(locus.reads_without_observation)
                + u64::from(locus.reads_discarded_by_cap),
            1_500,
            "the locus still accounts for every read fetched"
        );
    }

    /// A tract near the contig **end** splits a short right flank end-to-end: the fetched margin
    /// is clamped, so `SsrDetail.right_flank` is shorter than the left, and `reference_bases` is
    /// still the tract only. This is the contig-end-clamp path through the generator's
    /// flank-split (proven not to misalign; pinned here as a regression anchor).
    #[test]
    fn a_near_end_tract_splits_a_short_right_flank() {
        use crate::read::input::test_fixtures::read_named_with_length;
        // A read clear of the tract's query span (82..100) — zero coverage, so the locus is
        // present-but-empty and only the flank split is under test.
        let (_reference_dir, _bam_dir, reads) =
            sample_reads_with(&[read_named_with_length("r0", 0, 1, 30)]);
        let mut generator = ssr_generator(); // flank_bp = 10
        // chr1 is 100 bp; tract [92, 96] → margin [82, 100], right flank clamped to 100 - 96 = 4.
        let segment = tract(92, 96);

        generator.begin_segment(span(92, 96));
        let locus = generator
            .next_locus(&segment, &reads)
            .expect("no error")
            .expect("one locus");

        assert_eq!(locus.region, span(92, 96));
        assert_eq!(
            &*locus.reference_bases,
            &b"AAAAA"[..],
            "the 5 bp tract only"
        );
        match locus.kind {
            LocusKind::Ssr(SsrDetail {
                left_flank,
                right_flank,
                ..
            }) => {
                assert_eq!(left_flank.len(), 10, "full left flank");
                assert_eq!(right_flank.len(), 4, "clamped right flank (100 - 96)");
            }
            other => panic!("expected an Ssr locus, got {other:?}"),
        }
    }

    // --- the parity oracle against the deleted caller ---------------------------

    /// **The north-star parity oracle** (spec §6, §9.3): every COMPLETE observation ng produces
    /// must match production's `SsrLocusObs.observed` byte for byte, in bases and count, **with the
    /// cap disabled** on a fixture shallower than any cap.
    ///
    /// **Production's side is frozen since that caller was deleted** — see the recorded block in
    /// the body. The paragraphs below describe the comparison as it was built.
    ///
    /// Both sides ran over the **same** reads and the **same** reference frame. ng runs its real
    /// per-read pipeline (`classify::classify_read`, the `SsrFlatGapAligner` over production's
    /// `PerQualityEmission` table) and its real `tally`; production ran its real delimiter
    /// (`delimit_read` over `HmmModel`) and its real `tally`. The two aligners are *different
    /// algorithms* — this is exactly the claim under test: that on the complete class of clean reads
    /// they delimit the identical tract. The generator's fetch/cap are not exercised here because,
    /// with the cap disabled, they only select reads and never change an observation's bases; that
    /// wiring is covered by the generator tests above. The production side reproduced
    /// production's *clean-read* classify path (extract → delimit → the Region tract), which for
    /// these Q40, fully-spanned reads is exactly what production's `classify_read` does — this was
    /// asserted, while production still ran, via `flank_truncated` and the quality floor.
    ///
    /// Partial observations are ng's new behaviour with **no oracle**: this checks only that they
    /// *exist* (ng's classify keeps as a partial the read production's delimiter drops as
    /// border-off-end), not their bytes (spec §6).
    ///
    /// **Scope.** The fixture is all-Match reads of the reference allele and a one-unit *contraction*
    /// (`CACA`), so it pins parity on the complete class over the delimiters' contraction and
    /// exact-match behaviour. An *expansion* allele or a soft-clipped read cannot be added without
    /// leaving this regime: both drive the tract past the reference-sized window, tripping the
    /// long-allele **widening** recovery — a branch this production-side reproduction deliberately
    /// does not model (it asserted `!flank_truncated`). Widening-path and soft-clip parity are not
    /// covered here.
    #[test]
    fn ng_complete_observations_match_frozen_production_byte_for_byte() {
        use crate::alignment::ssr_best_path_flat_gap::{
            SsrFlatGapAligner, ViterbiScratch as NgViterbiScratch,
        };
        use crate::alignment::{PerQualityEmission, StutterModel};
        use crate::bam::alignment_input::CigarOp;
        use crate::locus_generation::ReadWitness;
        use crate::read::aligned_read::AlignedRead;
        // A 6 bp G-flank + CACACA + 6 bp T-flank — production's own delimiter fixture frame.
        const FRAME: &[u8] = b"GGGGGGCACACATTTTTT";

        /// A reference-frame read: `seq` mapped all-Match at 1-based `pos`, Q40.
        fn mapped(seq: &[u8], pos: u64) -> AlignedRead {
            AlignedRead {
                qname: b"r".to_vec(),
                flag: 0,
                ref_id: 0,
                pos,
                mapq: 60,
                cigar: vec![CigarOp::Match(seq.len() as u32)],
                seq: seq.to_vec(),
                qual: vec![40u8; seq.len()],
                mate_ref_id: None,
                mate_pos: None,
                adaptor_boundary: None,
                read_group: ReadGroupId(0),
            }
        }

        // Three complete reads of the reference allele (CACACA) and two of a clean shorter allele
        // (CACA, one CA unit deleted) — two distinct observed sequences, so the parity is not a
        // single-allele triviality. Plus one partially-covering read (left flank + two CA units, no
        // right flank) — the read ng keeps as a partial and production drops as border-off-end.
        let reads = [
            mapped(FRAME, 1),
            mapped(FRAME, 1),
            mapped(FRAME, 1),
            mapped(b"GGGGGGCACATTTTTT", 1),
            mapped(b"GGGGGGCACATTTTTT", 1),
            mapped(b"GGGGGGCACA", 1),
        ];

        // --- ng: real classify + real tally ---------------------------------
        let ng_locus = SsrLocus {
            repeat: SsrSegment::new("chr1".into(), 7, 12, Motif::new(b"CA").unwrap(), 1.0)
                .unwrap()
                .into(),
            tract_with_margin_bases: FRAME.into(),
            margin_start: Position(1),
        };
        let aligner = SsrFlatGapAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let mut ng_scratch = NgViterbiScratch::new();
        let mut qual_buffer = Vec::new();
        let mut counts = SsrGeneratorCounts::default();
        let ng_outcomes: Vec<_> = reads
            .iter()
            .map(|read| {
                classify::classify_read(
                    read,
                    &ng_locus,
                    &aligner,
                    &stutter,
                    &mut ng_scratch,
                    &mut qual_buffer,
                )
            })
            .collect();
        let ng = tally::tally(
            reads.iter().zip(ng_outcomes),
            ng_locus.repeat.start(),
            ng_locus.repeat.end(),
            &mut counts,
        );
        let ng_complete: Vec<(Vec<u8>, u32)> = ng
            .observations
            .iter()
            .filter(|obs| obs.read_witness == ReadWitness::Complete)
            .map(|obs| (obs.bases.to_vec(), obs.num_obs))
            .collect();

        // --- production: recorded ---------------------------------------------
        // What production's `delimit_read` over `HmmModel`, driven through its clean-read classify
        // path, and its `tally` gave on these six reads at commit `d9e7b076`, before its deletion:
        // the two complete alleles with their counts, and the partially-covering read dropped as
        // border-off-end. This is the value the test asserted ng equal to while it still called
        // production live, and production's code did not change after `d9e7b076`. That test also
        // asserted both simplifications its glue made — no read was window-truncated, and every
        // tract base cleared production's quality floor.
        let production_observed: Vec<(Vec<u8>, u32)> =
            vec![(b"CACA".to_vec(), 2), (b"CACACA".to_vec(), 3)];

        // --- the parity assertion -------------------------------------------
        assert_eq!(
            ng_complete,
            vec![(b"CACA".to_vec(), 2), (b"CACACA".to_vec(), 3)],
            "ng's two complete alleles, sorted by bytes with their counts"
        );
        assert_eq!(
            ng_complete, production_observed,
            "ng's complete observations match production's SsrLocusObs.observed byte for byte"
        );

        // ng's classify keeps a partial where production's delimiter drops the same read.
        assert_eq!(
            counts.observations_partial, 1,
            "ng's classify keeps the partially-covering read as a partial"
        );
        // Production's delimiter dropped that same read as border-off-end (recorded above), so
        // there is no oracle for the partial's bytes.
    }
}
