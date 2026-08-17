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

pub mod pileup;
pub mod ssr;
mod witness;

/// The witness vocabulary, re-exported so no consumer's import path names the file it
/// lives in (arch *Module home*).
pub use witness::{LocusLen, ReadWitness, WitnessedLocusPositions};

use crate::ng::read::input::{IngestError, SampleReads};
use crate::ng::ref_seq::RefSeqError;
use crate::ng::region_typing::segment_criteria::{Motif, SsrSegment};
use crate::ng::region_typing::{RegionKind, TypedRegion};
use crate::ng::types::{GenomePosition, GenomeRegion, Position, ReadGroupId};
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
    pub observations: Vec<SequenceObservation>,
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
    /// A [`Complete`](ReadWitness::Complete) observation counts its `num_obs` at every
    /// position; a [`Partial`](ReadWitness::Partial) run counts it over the stretch
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
        for obs in &self.observations {
            // **This clamp is the guard, not a second one.** An earlier comment here
            // called the bound "a producer invariant, enforced where `ReadWitness` is
            // minted", which overstated it twice over: `Partial`'s fields are public, so
            // a run need not have come from `from_left`/`from_right` at all; and even one
            // that did was clamped against *some* `LocusLen`, which nothing ties to the
            // locus it ends up on. `ReadWitness` cannot know its own locus, so the
            // invariant is not expressible on the type — it can only be checked here,
            // against the region actually in hand.
            //
            // Clamping rather than `debug_assert`: this is a consumer-side derivation run
            // over whole cohorts, and a debug-only guard compiles out of the release build
            // this repo actually runs (a trap it has recorded hitting twice). Clamping
            // keeps the derivation total on any input.
            //
            // **Per run since C2, and that is the point of the set.** Clamping the run the
            // witness *encloses* would count depth straight through a hole — the read
            // credited with positions it never saw, at exactly the spliced loci the set
            // exists to describe. So each run is clamped and added on its own; a witness of
            // one run, which is every witness on DNA-seq, walks the same slice as before.
            match &obs.read_witness {
                ReadWitness::Complete => {
                    for slot in &mut depth[0..len] {
                        *slot = slot.saturating_add(obs.num_obs);
                    }
                }
                ReadWitness::Partial { positions } => {
                    for (start, end) in positions.runs() {
                        let from = (start as usize).min(len);
                        let to = (end as usize).min(len).max(from);
                        for slot in &mut depth[from..to] {
                            *slot = slot.saturating_add(obs.num_obs);
                        }
                    }
                }
            }
        }
        depth
    }

    /// This locus's length, for the [`ReadWitness`] predicates that need it.
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
    /// [`Complete`](ReadWitness::Complete) ones.
    ///
    /// A partial is a lower bound that mis-scores as a *short* allele until a censored
    /// likelihood models it (step 7), so reaching the partials is a deliberate act:
    /// this iterator is the guard (spec §3).
    pub fn complete_observations(&self) -> impl Iterator<Item = &SequenceObservation> + '_ {
        self.observations
            .iter()
            .filter(|obs| obs.read_witness == ReadWitness::Complete)
    }

    /// The last reference base this observation covers.
    ///
    /// `region.end` on every well-formed observation, and named rather than read off the
    /// field because the cohort merge groups by it: production's grouping arithmetic is
    /// `pos + max(span, 1) − 1` (`reach`, `var_calling/cohort_integration.rs`), and
    /// whoever compares ng's rule with production's should find one place where the two
    /// agree rather than an open-coded expression at each use
    /// (`doc/devel/ng/arch/cohort_merge.md` §2).
    ///
    /// **It agrees with production's answer on every region below the top of the
    /// coordinate space, and it cannot overflow.** Production's own expression saturates
    /// on both operations and cannot overflow either; what fails at the ceiling is
    /// *reaching* it from a [`GenomeRegion`], because that needs the span and
    /// [`GenomeRegion::len`] computes `end + 1` before subtracting — a debug panic at
    /// `end == u64::MAX`, and a length of 0 in the release profile, which has overflow
    /// checks off. The larger of the two ends is the same number on every input — `end`
    /// whenever the region is well formed, `start` when it is not — and reaches it in one
    /// comparison, with no span in between.
    ///
    /// **At the ceiling the two part company by one, and this form is the right one.**
    /// Production saturates the addition before subtracting, so a one-base region at
    /// `u64::MAX` reaches `u64::MAX − 1` under its expression and `u64::MAX` under this
    /// one — and the last base of that region is `u64::MAX`. Pinned by
    /// `tests::a_locus_at_the_coordinate_ceiling_reaches_its_own_end`, so nobody restores
    /// "agreement" by reintroducing the arithmetic that panics.
    ///
    /// The inverted case is worth naming because [`GenomeRegion`] has public fields and
    /// no constructor enforcing `start <= end`: reading `region.end` there would put an
    /// observation's reach *before* its own first base, and a walk keyed on "does the
    /// next position fall within the reach" would close every locus immediately.
    pub fn reach(&self) -> Position {
        self.region.end.max(self.region.start)
    }

    /// Where this observation begins, genome-wide — the key the cohort merge orders on.
    ///
    /// A [`Position`] alone does not identify a base, and every consumer of this one
    /// compares it across samples and across contigs: the merge's k-way walk keys on it
    /// (`LocusCloser`) and so does the observation cache, which must know that a later
    /// contig lies beyond every position of the one before it.
    pub fn start_position(&self) -> GenomePosition {
        GenomePosition {
            contig: self.region.contig,
            position: self.region.start,
        }
    }

    /// The last base this observation covers, genome-wide — [`reach`](Self::reach) with
    /// its contig, and the sibling of [`start_position`](Self::start_position).
    pub fn reach_position(&self) -> GenomePosition {
        GenomePosition {
            contig: self.region.contig,
            position: self.reach(),
        }
    }

    /// How many reads here showed something other than the reference — the number the
    /// cohort merge's keep rule sums across a locus
    /// (`doc/devel/ng/spec/cohort_merge.md` §4.3).
    ///
    /// **Counted over the [`Complete`](ReadWitness::Complete) observations only, and
    /// that is forced by what a partial's bases are.** A partial's `bases` cover only the
    /// stretch its read witnessed, so comparing them against this locus's whole
    /// `reference_bases` would report a partial that saw less than the whole locus as
    /// non-reference — including a read that agreed with the reference over every base it
    /// actually saw. The census writer makes the same comparison over the same subset
    /// (`parameter_estimation/joint/census.rs`, `add_generic`), which is what
    /// [`SequenceObservation::matches_reference`] exists to keep in one place.
    ///
    /// The cost is that a variant witnessed only by partial reads does not reach the
    /// keep threshold. Scoring a partial needs a censored likelihood that does not exist
    /// yet (`spec/locus_generation.md` §3), so this is where the line sits today.
    ///
    /// **The sum saturates, and it costs nothing here**: the only consumer compares the
    /// total against `min_alt_obs`, default 2, so a total pinned at `u32::MAX` gives the
    /// same verdict as the true one. The generic pre-pass panics on the same sum instead
    /// (`parameter_estimation/generic/depth_and_alt_reads.rs`) because its number is a
    /// histogram key, where a wrapped count would be scored as a shallow site.
    ///
    /// The pre-pass counts this same quantity per site under the name `alt_reads`.
    /// *Non-reference read* is the spec's word (`spec/cohort_merge.md` §1.3); *alt*
    /// survives in `min_alt_obs`, which is production's parameter name.
    pub fn non_reference_reads(&self) -> u32 {
        self.complete_observations()
            .filter(|obs| !obs.matches_reference(&self.reference_bases))
            .fold(0u32, |total, obs| total.saturating_add(obs.num_obs))
    }
}

/// One distinct `(bases, witness, read group)` the reads showed at a locus, with the
/// pooled support of every read that showed it.
///
/// The fields between `num_obs` and `chain_ids` are the per-read moments the SNP
/// filters read (strand bias, base-quality error, the MAPQ multi-mapper test); an STR
/// model reduces to `num_obs` alone. Modelled on production's per-allele shape
/// (`AlleleObservation` + `AlleleSupportStats`), minus `placed_start`, which no model
/// consumes (spec §6).
///
/// **One entry is not one allele.** The identity has three axes, not one —
/// `(bases, read_witness, read_group)` — so a consumer that wants per-allele totals
/// must aggregate over witness *and* group, and one that treats each entry as an
/// allele will count the same allele several times. The aggregation is exact: every
/// support field is additive, and the merged entries share their `bases` and
/// `read_witness` by construction (spec §6).
#[derive(Debug, Clone, PartialEq)]
pub struct SequenceObservation {
    /// The observed bases — allele content, in **read** coordinates.
    pub bases: Box<[u8]>,
    /// How much of the locus a read of this sequence spanned. **Part of the
    /// identity**: a [`Complete`](ReadWitness::Complete) and an
    /// [`Partial`](ReadWitness::Partial) run of the same `bases` are different
    /// evidence and stay separate entries (spec §3).
    pub read_witness: ReadWitness,
    /// Which read group — one `@RG`, i.e. one lane — these reads came from. **Part of
    /// the identity**, so an allele supported from several groups is several observations.
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
    /// **The ids of every read folded here** — one per read, or one per read *pair* with
    /// its mates collapsed onto a single id, and deduplicated, so this is a count of reads
    /// only where no pair overlapped itself.
    ///
    /// It answers two questions, and the second is why the reference-matching reads are
    /// named too (the owner's ruling of 2026-08-17). **Which haplotype a read came from**,
    /// which is what lets a later step chain observations at neighbouring loci into one. And
    /// **whether a read was here at all**: when a cohort locus spans several of a sample's
    /// records, a read's allele over that locus is what it showed at each of them, so a read
    /// that covered a position and agreed with the reference has to be told apart from one
    /// that never reached it. Until the ruling those were the same absence, and the merge
    /// could only invent a reference stretch the read never saw or throw the read away
    /// (`doc/devel/ng/impl_plan/cohort_merge.md` B0).
    ///
    /// **An id names a read within one walk**, and a read that straddles the boundary
    /// between two walked regions is met twice and named twice. Nothing downstream links a
    /// read across such a boundary, because a segment is never cut and no locus crosses one
    /// (`doc/devel/ng/spec/run_streaming.md` §4.3).
    ///
    /// Empty on the STR path, which does not phase and does not need to: an STR locus is one
    /// record, so [`ReadWitness`] already says whether a read spanned it.
    pub chain_ids: Vec<ChainId>,
}

impl SequenceObservation {
    /// Whether these reads showed the reference's own bases over `reference_bases`.
    ///
    /// **The one place the comparison is written** — a byte comparison, which is all it
    /// has ever been, but not open-coded at each use, because two spellings of one test
    /// are two things that can disagree (`doc/devel/ng/arch/cohort_merge.md` §2). Its
    /// callers each keep their own sum: the census writer needs the answer per read
    /// group, the generic pre-pass needs a per-site count beside a depth, and the cohort
    /// merge needs a flat total per locus. The predicate is the shared thing; the sums
    /// are not.
    ///
    /// **It is equality, not containment.** A deletion's bases are frequently a prefix of
    /// the locus's reference — `AC` where the reference reads `ACGT` — and that is a
    /// different sequence, not a shorter spelling of the same one. Pinned by
    /// `tests::matches_reference_compares_the_bases_it_is_given`, which is there because
    /// a `starts_with` in place of the `==` passed every other test in the crate.
    ///
    /// **Raw bytes, so both sides must already be canonical.** ng's reference fetch
    /// uppercases ACGT and folds everything else to N (`ref_seq.rs`), so a soft-masked
    /// reference base cannot come back lowercase and read as a variant — but that is a
    /// dependency this predicate carries rather than enforces, and it would be silent if
    /// it broke.
    ///
    /// **What still belongs to a model is which observations to ask about.** The generic
    /// pre-pass decides what counts as an alternative read for its fit — the subset, the
    /// depth cap, the read-group grain (`parameter_estimation/generic/depth_and_alt_reads.rs`,
    /// `arch/parameter_prepass_generic.md` §2.3) — and this predicate does not touch any
    /// of that. It answers one question about one observation.
    ///
    /// **`reference_bases` must cover the same stretch these `bases` do.** For a
    /// [`Complete`](ReadWitness::Complete) observation that is the whole locus's
    /// [`reference_bases`](SampleLocusObservations::reference_bases); for a
    /// [`Partial`](ReadWitness::Partial) one it is not, since a partial's bases stop
    /// where its read's witness stopped, and handing it the whole locus's reference
    /// would report a read that matched everything it saw as non-reference. Callers that
    /// cannot supply the matching stretch should stay on
    /// [`complete_observations`](SampleLocusObservations::complete_observations).
    pub fn matches_reference(&self, reference_bases: &[u8]) -> bool {
        *self.bases == *reference_bases
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

    /// What this generator has counted so far, beside the shared [`LocusCounts`] — or
    /// `None` if it counts nothing of its own.
    ///
    /// Readable at any point; final once the run is drained. **This exists because a
    /// boxed generator's own counters were otherwise unreachable**: `GeneratorSlot`
    /// erases the type, so the generic generator's nine counts had no reader that was
    /// not a test, and a walk that emitted nothing for a covered region counted the
    /// truncations that explained it into a struct nobody could see (Milestone C
    /// review).
    ///
    /// Defaulted to `None` so a generator with nothing to report — [`NoLoci`], a test
    /// fake — says so by saying nothing.
    fn counts(&self) -> Option<GeneratorCounts<'_>> {
        None
    }
}

/// What a generator counted, tagged by which generator counted it.
///
/// **The same shape [`LocusKind`] uses**, and for the same reason: a common surface
/// with a per-kind payload, where the payload's type is the kind's own. A trait
/// method cannot return an associated type through `dyn`, and a downcast would move
/// a compile-time question to run time — so the kinds are enumerated here, exactly as
/// this module already enumerates them for loci ([`LocusKind::Ssr`] carrying
/// [`SsrDetail`]) and for slots ([`GeneratorSet`] has one named field per kind).
///
/// Borrowed, not owned: a caller reads a running tally rather than taking a snapshot
/// it then has to keep fresh.
#[derive(Debug)]
#[non_exhaustive]
pub enum GeneratorCounts<'a> {
    /// The generic (SNP/indel) pileup generator's run-level counts.
    Pileup(&'a pileup::PileupGeneratorCounts),
    /// The STR generator's run-level counts.
    Ssr(&'a ssr::SsrGeneratorCounts),
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
///
/// # Every generator failure names the region it happened over (owner, 2026-07-29)
///
/// A region is the unit of this whole step: it is what `begin_segment` takes, what the
/// counters are keyed to, and what the clamp and the halo reason about. An error that did not
/// name one could be *believed of the wrong region*, and that is not hypothetical — the
/// generic generator shipped exactly that defect, a failed read charged to the next region
/// after it had emitted all of its own loci, and the missing region is why nobody noticed
/// ([Milestone C review](../../../doc/devel/reports/reviews/ng_locus_generation_pileup_generator_c_2026-07-29.md)).
///
/// So the four generator-raised variants carry a [`GenomeRegion`] and **none of them has a
/// `#[from]` conversion**. That is the enforcement, not an oversight: with a blanket `From`, a
/// bare `?` compiles and silently produces an error with no region, which is the state this
/// change exists to make unreachable. Attaching at the `?` site costs one `map_err` and is the
/// only moment the region is known for certain.
///
/// [`RepeatCatalog`](Self::RepeatCatalog) is the exception and carries none, because the
/// region *stream* is what failed — there is no region to name, that being the point.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum LocusGenerationError {
    /// The reference's repeat catalog — where a region stream comes from — could not be
    /// read: the file is missing, describes another reference, or cannot answer the policy
    /// asked for. The first case names the command that writes one, which is the only one a
    /// person can act on.
    #[error("reading the reference's repeat catalog during locus generation")]
    RepeatCatalog(#[from] crate::ng::repeat_catalog::RepeatCatalogError),
    /// The read query for a region could not be **opened** — a missing or unreadable
    /// index, a contig the file does not have, a region the planner rejects.
    ///
    /// Split from [`Reads`](Self::Reads) because the two are different operational
    /// problems: "the index query for this region could not be opened" is a broken
    /// input or a bad request, "the record stream broke 40 kb in" is a truncated or
    /// corrupt file. One variant rendered them identically, and its own doc used to
    /// deny that the first could happen at all.
    ///
    /// **Since D1 this also covers pointing a long-lived cursor at a region** — the same
    /// job, done once per region against a reader that stays open instead of once against
    /// a reader that is opened for it. One condition reaching here is *not* an open
    /// failure in the sense above: `CursorError::AfterFailure` says the cursor's file broke
    /// during an **earlier** region, which by the split's own logic belongs to
    /// [`Reads`](Self::Reads). It is unreachable through the generic generator, whose
    /// `failed` latch fires first and reports the original failure; if a caller ever meets
    /// it, the message it wants is the one already reported, not this one.
    #[error("the read query for {region} could not be opened")]
    OpenReadQuery {
        region: GenomeRegion,
        #[source]
        source: IngestError,
    },
    /// A generator was handed a **different sample** from the one it was opened for.
    ///
    /// **Each sample needs its own generator.** A generator opens a reader for one sample and
    /// keeps it positioned for a whole chromosome — that is what makes it fast — but it is
    /// lent a `&SampleReads` afresh on every call. Sharing one generator between samples would
    /// therefore answer every sample out of the first one's files: no error, no empty rows, and
    /// a cohort in which every individual looks identical to the first. This is that mistake,
    /// caught.
    ///
    /// It is a caller bug, not a data problem. Correct code builds a generator per sample and
    /// never sees it.
    #[error(
        "this generator was opened for another sample, so it cannot answer for {region}: \
         give each sample its own generator"
    )]
    ForeignSample { region: GenomeRegion },
    /// A read query failed **mid-stream**, or the alignment input was malformed: the
    /// open already succeeded and reads were flowing.
    #[error("read access over {region} failed during locus generation")]
    Reads {
        region: GenomeRegion,
        #[source]
        source: IngestError,
    },
    /// A reference fetch failed — a broken reference, or a region past a contig end.
    #[error("reference fetch over {region} failed during locus generation")]
    Reference {
        region: GenomeRegion,
        #[source]
        source: RefSeqError,
    },
    /// The pileup walk failed: a malformed read, reads out of coordinate order, a
    /// record wider than the span cap, or an exhausted chain-id space. Fatal and
    /// terminal for the walk that raised it
    /// (`locus_generation_pileup.md` §7).
    ///
    /// None of the variants above covers it — they name the *inputs* failing, and this
    /// names the walk over inputs it already accepted.
    #[error("the pileup walk over {region} failed")]
    Walker {
        region: GenomeRegion,
        #[source]
        source: pileup::WalkerError,
    },
}

impl LocusGenerationError {
    /// The region the failure is attributed to, or `None` for a failure of the region
    /// *stream* itself.
    ///
    /// The region a generator attaches is **the one it was working over** — the segment
    /// it was given, not the wider span it queried, since the segment is the unit a
    /// caller can act on. The one exception is a helper that only ever knows a span
    /// (the STR read fetch), which says so where it attaches.
    ///
    /// Here so a consumer — a log line, a per-region tally — can ask without matching
    /// every variant, which is what would rot the moment a variant is added.
    pub fn region(&self) -> Option<GenomeRegion> {
        match self {
            // A catalog failure is a failure of the region *stream*, so it has no region to
            // name — that being the point of keeping it apart from the rest.
            LocusGenerationError::RepeatCatalog(_) => None,
            LocusGenerationError::OpenReadQuery { region, .. }
            | LocusGenerationError::ForeignSample { region }
            | LocusGenerationError::Reads { region, .. }
            | LocusGenerationError::Reference { region, .. }
            | LocusGenerationError::Walker { region, .. } => Some(*region),
        }
    }
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
/// The trait object carries **no `Send` bound**: v1 is single-threaded (`locus_generation.md` §9). If a
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

    /// What the generator in this slot has counted, or `None` — for an unfilled slot,
    /// or a generator that counts nothing of its own.
    fn counts(&self) -> Option<GeneratorCounts<'_>> {
        match self {
            GeneratorSlot::Generator(generator) => generator.counts(),
            GeneratorSlot::Unfilled(_) => None,
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

    /// What the **STR** generator has counted, if one is filled and counts anything.
    pub fn ssr_counts(&self) -> Option<GeneratorCounts<'_>> {
        self.ssr.counts()
    }

    /// What the **generic** generator has counted, if one is filled and counts anything.
    ///
    /// One accessor per slot rather than one keyed by [`RegionKind`]: the kinds are
    /// already three named fields here, and a key would have to carry a payload
    /// (`RegionKind::SsrSegment` holds a segment) that has nothing to do with reading a
    /// tally.
    pub fn generic_counts(&self) -> Option<GeneratorCounts<'_>> {
        self.generic.counts()
    }

    /// What the **`SsrBundle`** generator has counted, if one is filled and counts
    /// anything. `None` today: the slot has no generator (spec §10).
    pub fn ssr_bundle_counts(&self) -> Option<GeneratorCounts<'_>> {
        self.ssr_bundle.counts()
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
/// Generic over the region stream `T` so a `Vec` can stand in for the catalog's own
/// region reader in tests. The generator set is a concrete [`GeneratorSet`] — the
/// per-run swap lives in its trait-object slots, not in a type parameter.
pub struct SampleLocusObservationsIterator<T> {
    regions: T,
    // **`generators` must be dropped before `reads`** (C4, plan 3 — the choice the
    // pileup generator's FIXME left to this step). Declaring it first is no longer
    // what enforces that: the `Drop` impl below releases it explicitly, so a future
    // edit that reorders these fields cannot silently undo the property. The order
    // is kept anyway, as the cheaper of the two mechanisms.
    //
    // Since 2026-07-28 a region stream owns its files by `Arc` rather than
    // borrowing them, so it may outlive the `SampleReads` that made it — and a
    // generator holds one such reader across `next_locus` calls. Rust drops
    // fields in declaration order: with `reads` first, a reader a generator was
    // still holding would fold its drop tally into an `AlignmentFile` that only
    // that reader owns and that is freed in the same breath. The reads already
    // emitted are unaffected and the pooled reader still goes back; what is lost
    // is the ability to *read* that tally through `SampleReads::counts` — a
    // silent under-report of drop rates, not a crash.
    //
    // **How long a generator holds one widened at D1, and the release below is
    // now the ordinary case rather than the exceptional one.** The STR generator
    // still holds a per-region stream and drops it when the region drains, so a
    // run driven to exhaustion leaves it holding nothing. The generic generator
    // holds a *cursor per chromosome*: `end_walk` clears only the region walk, so
    // the cursor survives every region and is released when the generator is.
    // Draining the run no longer empties it.
    //
    // **No test can fail if this breaks, and that is a property of the types
    // rather than an omission** — see the `Drop` impl for what it would take.
    generators: GeneratorSet,
    reads: SampleReads,
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

    /// The generator set, for the per-generator counts the shared tally does not carry
    /// ([`GeneratorSet::generic_counts`] and its siblings).
    ///
    /// The whole set rather than three more forwarding methods: it is handed out by
    /// `&`, so a caller can read every slot's tally and change none of them.
    pub fn generators(&self) -> &GeneratorSet {
        &self.generators
    }
}

impl<T, E> Iterator for SampleLocusObservationsIterator<T>
where
    T: Iterator<Item = Result<TypedRegion, E>>,
    E: Into<LocusGenerationError>,
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

impl<T, E> std::iter::FusedIterator for SampleLocusObservationsIterator<T>
where
    T: Iterator<Item = Result<TypedRegion, E>>,
    E: Into<LocusGenerationError>,
{
}

impl<T> Drop for SampleLocusObservationsIterator<T> {
    /// **Release the generators while `reads` is still alive.** A generator holds
    /// a region stream that owns its files by `Arc`, and that stream folds its
    /// drop tally into an `AlignmentFile` on the way out; if `reads` — the only
    /// other holder — has gone first, the tally lands in an object nobody can
    /// read and is freed with it. Silent: no crash, no wrong locus, just a drop
    /// rate under-reported by whatever the abandoned region had counted.
    ///
    /// Field order already gives this (Rust drops fields in declaration order,
    /// and `generators` is declared first). This is here because that made the
    /// property depend on a line's *position* in a struct, guarded by a comment —
    /// and the failure is invisible, so a reorder would not announce itself. With
    /// both, the order is an optimisation and this is the guarantee.
    ///
    /// **Nothing can test it, and the reason is worth stating rather than
    /// leaving as an absence** (owner-facing decision, 2026-07-30). Observing the
    /// tally after the drop needs a second handle onto the same files;
    /// `SampleReads` is deliberately not `Clone` and does not expose them, so the
    /// test costs a widening of that type — the shape is
    /// `read::input::open_bam::tests::a_stream_outliving_every_other_handle_still_banks_its_reader_and_tally`,
    /// one layer down, where the handle exists. Deleting this impl fails no test;
    /// it is a defence, not a checked invariant, and it is cheap because the
    /// replacement set allocates nothing.
    fn drop(&mut self) {
        drop(std::mem::replace(
            &mut self.generators,
            GeneratorSet::all_unimplemented(),
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::repeat_catalog::RepeatCatalogError;
    use crate::ng::types::{ContigId, Position, ReadGroupId};

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// A one-run partial witness, in the offset-and-length spelling these depth fixtures
    /// were written in. C2 made `Partial`'s payload a set; a fixture that says "four
    /// positions from offset three" still reads as the constraint it is.
    fn partial_run(offset_in_locus: u16, positions_covered: u16) -> ReadWitness {
        ReadWitness::Partial {
            positions: WitnessedLocusPositions::one_run_from_offset_and_length(
                offset_in_locus,
                positions_covered,
            )
            .expect("a run covering at least one position"),
        }
    }

    /// An observation of `bases` with `num_obs` reads at a given witness — the moment
    /// fields are irrelevant to the depth derivation, so they are fixed.
    fn obs(bases: &[u8], read_witness: ReadWitness, num_obs: u32) -> SequenceObservation {
        SequenceObservation {
            bases: Box::from(bases),
            read_witness,
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

    fn locus(region: GenomeRegion, observed: Vec<SequenceObservation>) -> SampleLocusObservations {
        SampleLocusObservations {
            region,
            reference_bases: Box::from(&b""[..]),
            observations: observed,
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
            Err(LocusGenerationError::Reads {
                region: region(1, 1),
                source: IngestError::NoFiles,
            })
        }
    }

    /// The iterator drains a multi-kind stream and accounts every region, yielding nothing
    /// when all kinds are unimplemented — the shape run on its own (spec §2, §13.2, §13.3).
    #[test]
    fn the_iterator_drains_and_accounts_a_multi_kind_stream() {
        let (_reference_dir, _bam_dir, reads) = sample_reads_over_fixture();
        let regions = vec![
            Ok::<_, RepeatCatalogError>(typed(RegionKind::Generic, 1, 10)),
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
            Ok::<_, RepeatCatalogError>(typed(RegionKind::Generic, 1, 10)),
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
            Ok::<_, RepeatCatalogError>(typed(RegionKind::Generic, 1, 10)), // generic slot → 3 loci at start=1
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
            Ok::<_, RepeatCatalogError>(typed(RegionKind::Generic, 1, 10)),
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
            Some(Err(LocusGenerationError::Reads { .. })) => {}
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
            Ok::<_, RepeatCatalogError>(typed(RegionKind::Generic, 1, 10)),
            Err(RepeatCatalogError::ContigTableMismatch {
                detail: "an injected stream failure".to_string(),
            }),
            Ok(typed(RegionKind::Generic, 20, 30)),
        ];
        let mut iterator = SampleLocusObservationsIterator::new(
            regions.into_iter(),
            reads,
            GeneratorSet::all_unimplemented(),
        );

        match iterator.next() {
            Some(Err(LocusGenerationError::RepeatCatalog(_))) => {}
            other => panic!("expected a fatal wrapped catalog error, got {other:?}"),
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
            observations: vec![SequenceObservation {
                bases: Box::from(&b"T"[..]),
                read_witness: ReadWitness::Complete,
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
        assert_eq!(generic.observations[0].num_obs, 9);

        let ssr = SampleLocusObservations {
            region: region(10_442, 10_461),
            reference_bases: Box::from(&b"ATATATATATATATATATAT"[..]),
            observations: Vec::new(),
            reads_without_observation: 3,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(SsrDetail {
                motif: Motif::new(b"AT").unwrap(),
                left_flank: Box::from(&b"CCCGGG"[..]),
                right_flank: Box::from(&b"TTTAAA"[..]),
            }),
        };
        // Zero coverage is a real observation: the locus exists with an empty table.
        assert!(ssr.observations.is_empty());
        assert!(matches!(ssr.kind, LocusKind::Ssr(_)));

        let bundle = SampleLocusObservations {
            region: region(200, 260),
            reference_bases: Box::from(&b"N"[..]),
            observations: Vec::new(),
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::SsrBundle,
        };
        assert_eq!(bundle.kind, LocusKind::SsrBundle);
    }

    /// A complete and a partial observation of the *same* bases are distinct evidence,
    /// so they must not compare equal — the property the `(bases, read_witness)`
    /// dedup key rests on (spec §3).
    #[test]
    fn same_bases_differ_by_read_witness() {
        let complete = SequenceObservation {
            bases: Box::from(&b"ATATAT"[..]),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: 1,
            num_fwd: 1,
            q_sum: 0.0,
            mapq_sum: 60,
            mapq_sum_sq: 3_600,
            placed_left: 0,
            chain_ids: Vec::new(),
        };
        let partial = SequenceObservation {
            read_witness: ReadWitness::from_left(6, LocusLen::from_positions(6))
                .expect("a run covering at least one position"),
            ..complete.clone()
        };
        assert_ne!(complete, partial);
        assert_ne!(complete.read_witness, partial.read_witness);
    }

    /// Depth derives correctly from the read witness (spec §13.5): the vector has
    /// `region.len()` entries; a `Complete` raises every position, a left-flush run of
    /// `n` only the leftmost `n`, a right-flush run of `n` only the rightmost `n`. A
    /// 10-position locus with one complete (×3), one left-partial reaching 4 (×2), one
    /// right-partial reaching 3 (×5).
    #[test]
    fn depth_derives_from_read_witness() {
        let l = locus(
            region(1, 10),
            vec![
                obs(b"AAAAAAAAAA", ReadWitness::Complete, 3),
                obs(
                    b"AAAA",
                    ReadWitness::from_left(4, LocusLen::from_positions(10))
                        .expect("a run covering at least one position"),
                    2,
                ),
                obs(
                    b"AAA",
                    ReadWitness::from_right(3, LocusLen::from_positions(10))
                        .expect("a run covering at least one position"),
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
                obs(b"A", ReadWitness::Complete, 7),
                obs(b"T", ReadWitness::Complete, 2),
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
                ReadWitness::from_left(9, LocusLen::from_positions(3))
                    .expect("a run covering at least one position"),
                4,
            )],
        );
        assert_eq!(clamped.num_obs_along_locus(), vec![4, 4, 4]);
    }

    /// **The clamp in `num_obs_along_locus` is the guard, and this is the input that needs
    /// it** — a witness whose runs reach past the locus *without* having come through a
    /// constructor.
    ///
    /// The test above cannot reach it: `from_left(9, LocusLen(3))` is clamped by the
    /// constructor, so the run arriving here is already `0..3` and deleting both `.min(len)`
    /// calls leaves everything green (Milestone C review). That is exactly the gap the
    /// clamp's own comment describes — `Partial`'s field is public and
    /// `WitnessedLocusPositions` cannot know which locus it ends up attached to, so a run
    /// need not have been clamped against *this* one.
    ///
    /// Unclamped, the first run indexes `depth[0..20]` on a 3-slot vector: a panic, in a
    /// release build, on a derivation run over whole cohorts. The second run is past the
    /// locus entirely and must contribute nothing rather than wrap.
    #[test]
    fn depth_clamps_a_witness_that_reaches_past_the_locus_it_is_attached_to() {
        let over_long = ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 20), (40, 44)])
                .expect("two runs, both overrunning a 3-position locus"),
        };
        let l = locus(region(1, 3), vec![obs(b"AAA", over_long, 4)]);
        assert_eq!(
            l.num_obs_along_locus(),
            vec![4, 4, 4],
            "the first run is cut at the locus's end and the second falls outside it",
        );
    }

    /// Depth over an **interior** run — flush with neither border. This is the case the
    /// reshape exists to represent (a read blind in the middle of a footprint: an interior
    /// `N`, a ref-skip), and the only one where a wrong window neither panics nor clamps, so
    /// it is also the only one where an off-by-one in the arm is purely silent.
    #[test]
    fn depth_over_an_interior_run_raises_only_the_witnessed_stretch() {
        let l = locus(region(1, 10), vec![obs(b"AAAA", partial_run(3, 4), 7)]);
        // positions: 1  2  3  4  5  6  7  8  9 10
        //            .  .  .  7  7  7  7  .  .  .
        assert_eq!(l.num_obs_along_locus(), vec![0, 0, 0, 7, 7, 7, 7, 0, 0, 0]);
    }

    /// **Depth over a witness with a hole raises the runs and not the gap** — the derivation
    /// C2's set exists for, and the one number that says whether the hole survived the trip
    /// from the fold to a depth profile.
    ///
    /// A spliced read witnessing positions 1–3 and 8–10 of a ten-position locus must leave
    /// 4–7 at zero. Summing over the extent that *encloses* its runs raises all ten instead,
    /// which is the read being credited with an intron it never sequenced — no panic, no
    /// clamp, just a depth that is 60 % too wide. That mutation is what this test fails
    /// under; nothing else in the suite notices it, because every witness the fold mints
    /// before C3 has one run and the two readings agree on those.
    #[test]
    fn depth_over_a_witness_with_a_hole_leaves_the_hole_at_zero() {
        let spliced = ReadWitness::Partial {
            positions: WitnessedLocusPositions::from_half_open_runs([(0, 3), (7, 10)])
                .expect("two runs"),
        };
        let l = locus(region(1, 10), vec![obs(b"AAACCC", spliced, 5)]);
        // positions: 1  2  3  4  5  6  7  8  9 10
        //            5  5  5  .  .  .  .  5  5  5
        assert_eq!(l.num_obs_along_locus(), vec![5, 5, 5, 0, 0, 0, 0, 5, 5, 5]);
    }

    /// `complete_observations()` yields only the complete entries — the guard that a
    /// partial (a lower bound) is never scored as an exact allele.
    #[test]
    fn complete_observations_excludes_partials() {
        let l = locus(
            region(1, 6),
            vec![
                obs(b"ATATAT", ReadWitness::Complete, 4),
                obs(
                    b"ATATAT",
                    ReadWitness::from_left(6, LocusLen::from_positions(6))
                        .expect("a run covering at least one position"),
                    2,
                ),
                obs(b"ATGTAT", ReadWitness::Complete, 3),
                obs(
                    b"ATAT",
                    ReadWitness::from_right(4, LocusLen::from_positions(6))
                        .expect("a run covering at least one position"),
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

    // ---------------------------------------------------------------
    // The two derivations the cohort merge walks on
    // (`doc/devel/ng/arch/cohort_merge.md` §2).
    // ---------------------------------------------------------------

    /// A locus with reference bases, which the depth fixtures above do not need.
    fn locus_over_reference(
        region: GenomeRegion,
        reference_bases: &[u8],
        observed: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            reference_bases: Box::from(reference_bases),
            ..locus(region, observed)
        }
    }

    /// Production's grouping arithmetic, copied from `reach` in
    /// `var_calling/cohort_integration.rs` — what [`SampleLocusObservations::reach`] has
    /// to agree with, since the cohort merge's chaining rule is production's. It
    /// saturates on both operations, exactly as production's does.
    fn production_reach(pos: u64, span: u64) -> u64 {
        pos.saturating_add(span.max(1)).saturating_sub(1)
    }

    /// The reach is the locus's last base, and it is production's answer — checked
    /// against production's own arithmetic rather than against a copy of the result, on
    /// the two shapes the merge sees: a SNP, which covers one base, and a deletion,
    /// which covers several.
    #[test]
    fn the_reach_agrees_with_production_arithmetic() {
        let snp = locus_over_reference(region(10, 10), b"A", vec![]);
        assert_eq!(snp.reach(), Position(10));
        assert_eq!(
            snp.reach().get(),
            production_reach(snp.region.start.get(), snp.region.len())
        );

        let deletion = locus_over_reference(region(10, 14), b"ACGTA", vec![]);
        assert_eq!(deletion.reach(), Position(14));
        assert_eq!(
            deletion.reach().get(),
            production_reach(deletion.region.start.get(), deletion.region.len())
        );
    }

    proptest::proptest! {
        /// The agreement is a property of the whole well-formed domain, not of the two
        /// points the fixture above names — both of which start at 10. A `reach` that
        /// happened to be right at 10–10 and 10–14 and wrong at `start == 0`, or over a
        /// long span, passes those and fails this.
        ///
        /// Bounded below the ceiling deliberately: that is the one documented input where
        /// the two forms differ, and where production's cannot be evaluated from a
        /// `GenomeRegion` at all.
        #[test]
        fn the_reach_agrees_with_production_over_every_well_formed_region(
            start in 0u64..1_000_000,
            span in 1u64..1_000,
        ) {
            let well_formed = region(start, start + span - 1);
            let l = locus_over_reference(well_formed, b"", vec![]);
            proptest::prop_assert_eq!(
                l.reach().get(),
                production_reach(start, well_formed.len())
            );
        }
    }

    /// **A locus at the top of the coordinate space answers instead of panicking**, which
    /// is what keeps the merge's arithmetic on the saturating side of spec §11's trap.
    ///
    /// **This is the one input where the two ancestors part company**, and it parts them
    /// in two ways. `GenomeRegion::len` computes `end + 1` before subtracting, so
    /// `region(u64::MAX, u64::MAX).len()` panics with "attempt to add with overflow" in a
    /// debug build — a case that method's doc does not cover either way, since it
    /// promises saturation for an *inverted* region and says nothing about the ceiling.
    /// And production's expression, handed the span it would have, saturates the addition
    /// before subtracting and lands one short of the true last base.
    ///
    /// The assertion below is what stops someone restoring "agreement" by reintroducing
    /// arithmetic that both panics and answers `u64::MAX − 1` for a base at `u64::MAX`.
    #[test]
    fn a_locus_at_the_coordinate_ceiling_reaches_its_own_end() {
        let at_the_ceiling = locus_over_reference(region(u64::MAX, u64::MAX), b"A", vec![]);
        assert_eq!(at_the_ceiling.reach(), Position(u64::MAX));

        // Production's form, given the one-base span it would have here, is short by one.
        assert_eq!(production_reach(u64::MAX, 1), u64::MAX - 1);
    }

    /// **An inverted region does not put the reach behind the start.** `GenomeRegion`'s
    /// fields are public and nothing enforces `start <= end`, and a walk keyed on "does
    /// the next position fall within the reach" would close every locus at once if one
    /// ever answered with a base before its own first. Reading `region.end` gives 4 here;
    /// production's expression gives the start, and so does this.
    #[test]
    fn an_inverted_region_reaches_its_own_start() {
        let inverted = locus_over_reference(region(10, 4), b"", vec![]);
        assert_eq!(inverted.reach(), Position(10));
        assert!(inverted.reach() >= inverted.region.start);

        // `len()` saturates to 0 here and production's `span.max(1)` puts it back to 1,
        // so the two forms land on the start together — asserted, not claimed in prose.
        assert_eq!(
            inverted.reach().get(),
            production_reach(inverted.region.start.get(), inverted.region.len())
        );
    }

    /// The predicate is a byte comparison over the stretch it is handed — the same test
    /// the census writer used to spell inline.
    ///
    /// **Bases of a different length from the reference are what earn this test its
    /// keep**, in both directions. Two containment tests — the shape a reader reaching
    /// for "handle a partial" writes — pass everything else in the crate: replacing the
    /// `==` with `reference_bases.starts_with(&self.bases)`, or with
    /// `self.bases.starts_with(reference_bases)`, left every test in this module and
    /// every test in the census green before the three assertions below were added,
    /// because no fixture anywhere compared bases of a different length against the
    /// reference. A deletion's bases are shorter and an insertion's longer,
    /// which is the whole indel half of what this caller is for, and
    /// `non_reference_reads` is what the merge's keep threshold sums — so either mutant
    /// would have stopped an indel counting, and dropped indel-only loci below the
    /// threshold with nothing objecting.
    #[test]
    fn matches_reference_compares_the_bases_it_is_given() {
        let reference = obs(b"ACGT", ReadWitness::Complete, 3);
        let snp = obs(b"ACCT", ReadWitness::Complete, 2);
        let deletion = obs(b"AT", ReadWitness::Complete, 1);

        assert!(reference.matches_reference(b"ACGT"));
        assert!(!snp.matches_reference(b"ACGT"));
        assert!(!deletion.matches_reference(b"ACGT"));

        // Equality, not containment, and both directions of it. `ACGT` starting with
        // `AC` does not make `AC` the reference; `ACGTT` starting with `ACGT` does not
        // make the insertion the reference either.
        let trailing_deletion = obs(b"AC", ReadWitness::Complete, 1);
        let insertion = obs(b"ACGTT", ReadWitness::Complete, 1);
        let left_insertion = obs(b"TACGT", ReadWitness::Complete, 1);
        assert!(!trailing_deletion.matches_reference(b"ACGT"));
        assert!(!insertion.matches_reference(b"ACGT"));
        assert!(!left_insertion.matches_reference(b"ACGT"));

        // A different stretch is a different question, and the predicate answers the one
        // it was asked: the same bases match a reference that is those bases.
        assert!(deletion.matches_reference(b"AT"));
    }

    /// The keep rule's input: reads that differ from the reference, summed, and reads
    /// that agree contributing nothing however deep they are.
    ///
    /// The 40 reference reads are what makes this discriminating — an implementation
    /// that summed every observation, or that inverted the predicate, would answer 45 or
    /// 40 here rather than 5.
    #[test]
    fn non_reference_reads_sums_only_the_reads_that_differ() {
        let l = locus_over_reference(
            region(10, 13),
            b"ACGT",
            vec![
                obs(b"ACGT", ReadWitness::Complete, 40),
                obs(b"ACCT", ReadWitness::Complete, 3),
                obs(b"AGGT", ReadWitness::Complete, 2),
            ],
        );
        assert_eq!(l.non_reference_reads(), 5);

        let quiet = locus_over_reference(
            region(10, 13),
            b"ACGT",
            vec![obs(b"ACGT", ReadWitness::Complete, 40)],
        );
        assert_eq!(quiet.non_reference_reads(), 0);
    }

    /// **A partial is not counted, and the fixture says why.** This partial's read agreed
    /// with the reference over every base it saw — it simply stopped after two — so its
    /// bases are `AC` against a locus reference of `ACGT`. Counting it would report 6
    /// non-reference reads where the reads showed 2, on a locus where nothing but the
    /// witness is unusual, and at a threshold of 2 that turns ground the cohort agreed on
    /// into a built locus.
    #[test]
    fn a_partial_that_agreed_with_the_reference_is_not_counted_against_it() {
        let l = locus_over_reference(
            region(10, 13),
            b"ACGT",
            vec![
                obs(b"ACGT", ReadWitness::Complete, 10),
                obs(b"ACCT", ReadWitness::Complete, 2),
                obs(
                    b"AC",
                    ReadWitness::from_left(2, LocusLen::from_positions(4))
                        .expect("a run covering at least one position"),
                    4,
                ),
            ],
        );
        assert_eq!(l.non_reference_reads(), 2);
    }

    /// A locus nobody covered has no non-reference reads, and that is an answer rather
    /// than an error — the fold has to be total, since the merge asks this of every
    /// locus it walks past, most of which are quiet.
    #[test]
    fn a_locus_with_no_observations_has_no_non_reference_reads() {
        let uncovered = locus_over_reference(region(10, 13), b"ACGT", vec![]);
        assert_eq!(uncovered.non_reference_reads(), 0);
    }

    /// **The cost of the complete-only rule, pinned.** This partial's read *disagreed*
    /// with the reference over the two bases it saw, and its 7 reads still contribute
    /// nothing — so a variant witnessed only by partial reads never reaches
    /// `min_alt_obs`, and nothing downstream is emitted over that locus.
    ///
    /// The neighbouring partial test cannot see this: a partial that *agreed* answers 2
    /// whether partials are excluded or compared properly against their own stretch. This
    /// one separates today's rule from a partial-aware one, which is what makes it the
    /// test that will fail, loudly, when step 7's censored likelihood changes the
    /// decision.
    #[test]
    fn a_variant_seen_only_by_partial_reads_is_not_counted() {
        let l = locus_over_reference(
            region(10, 13),
            b"ACGT",
            vec![obs(
                b"AT",
                ReadWitness::from_left(2, LocusLen::from_positions(4))
                    .expect("a run covering at least one position"),
                7,
            )],
        );
        assert_eq!(l.non_reference_reads(), 0);
    }

    /// **Saturating, not wrapping** — spec §11's trap, at the boundary that separates the
    /// two. The threshold only asks whether the total reached `min_alt_obs`, so a capped
    /// total answers that question correctly where a wrapped one answers it backwards,
    /// and `num_obs` will arrive from a file on the psp path rather than only from a walk
    /// that caps its columns.
    #[test]
    fn non_reference_reads_saturates_rather_than_wrapping() {
        let l = locus_over_reference(
            region(10, 13),
            b"ACGT",
            vec![
                obs(b"ACCT", ReadWitness::Complete, u32::MAX),
                obs(b"AGGT", ReadWitness::Complete, u32::MAX),
            ],
        );
        assert_eq!(l.non_reference_reads(), u32::MAX);
    }
}
